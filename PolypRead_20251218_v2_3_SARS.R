#!/usr/bin/env Rscript

# ============================================================
# GEO scRNA-seq Complete Processing Pipeline
# VERSION: 2.3 (2025-12-18) - FINAL STABLE
#
# Fixed in v2.3 (vs v2.2):
# P0: infer_gsm_col NA handling (prevents crash on NA-containing columns)
# Robustness: gsm_id whitespace cleaning
# Robustness: Duplicate gsm_id detection with fail-fast
#
# Previous fixes (v2.2):
# - Metadata key logic (gsm_id vs sample separation)
# - 10x prefix detection for standard naming
# - Gene name preservation in deduplication
#
# Previous fixes (v2.1):
# - File connection management
# - Sparse matrix deduplication
# - H5 ENSG→Symbol path
# - EmptyDrops-aware filtering
# ============================================================

# ===== Configuration Section =====
DATA_ROOT <- "/home/h2048/data/source/1221"
METADATA_PATH <- "/home/h2048/data/source/reference/metadata_full_atlas_20251210.csv"
OUTPUT_DIR <- "/home/h2048/data/R/1221/Polyp"
QC_PLOTS_DIR <- "/home/h2048/data/R/1221/Polyp/qc_plots"
# ===== Feature name cleanup (NEW) =====
STRIP_GRCH38_PREFIX <- TRUE
REMOVE_SARS_GENES <- TRUE

GRCH38_PREFIX_PATTERN <- "^DEPRECATED_"
SARS_PREFIX_PATTERN <- "^SARS-"

# EmptyDrops parameters
RUN_EMPTYDROPS <- TRUE
EMPTYDROPS_LOWER <- 100
EMPTYDROPS_FDR_THRESHOLD <- 0.01

# Gene ID standardization parameters
STANDARDIZE_GENE_IDS <- TRUE
SYMBOL_DEDUP_METHOD <- "sum"
KEEP_UNMAPPED_GENES <- TRUE
MIN_MAPPING_RATE <- 0.5

# QC plot settings
MAX_PLOT_POINTS <- 5000

QC_THRESHOLDS <- list(
  nFeature_min = 200,
  nFeature_max = 6000,
  nCount_min = 500,
  nCount_max = 50000,
  percent_mt_max = 20,
  percent_rb_min = 5
)

# Safety knobs
MAX_WIDE_MATRIX_CELLS <- 30000
FAIL_ON_READMM_TRUNCATED <- FALSE

# ===== Load Libraries =====
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(ggplot2)
  library(patchwork)
  library(DropletUtils)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(jsonlite)
})

cat("============================================================\n")
cat("GEO/HRA scRNA-seq Processing Pipeline v2.3 (FINAL STABLE)\n")
cat("============================================================\n")
cat("Data root:  ", DATA_ROOT, "\n")
cat("Metadata:   ", METADATA_PATH, "\n")
cat("Output:     ", OUTPUT_DIR, "\n")
cat("QC plots:   ", QC_PLOTS_DIR, "\n")
cat("EmptyDrops: ", ifelse(RUN_EMPTYDROPS, "ENABLED", "DISABLED"), "\n")
if (RUN_EMPTYDROPS) {
  cat("  Lower bound:", EMPTYDROPS_LOWER, "UMIs\n")
  cat("  FDR threshold:", EMPTYDROPS_FDR_THRESHOLD, "\n")
}
cat(
  "Gene standardization:",
  ifelse(STANDARDIZE_GENE_IDS, "ENABLED", "DISABLED"),
  "\n"
)
if (STANDARDIZE_GENE_IDS) {
  cat("  Dedup method:", SYMBOL_DEDUP_METHOD, "\n")
  cat("  Keep unmapped:", KEEP_UNMAPPED_GENES, "\n")
}
cat("============================================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

BAD_SAMPLES_TXT <- file.path(OUTPUT_DIR, "bad_samples.txt")
writeLines(
  c(
    "dataset_id\tgsm_id\tformat\treason\tfile",
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  con = BAD_SAMPLES_TXT
)

append_bad_sample <- function(
  dataset_id,
  gsm_id,
  format,
  reason,
  file = NA_character_
) {
  line <- paste(dataset_id, gsm_id, format, reason, file, sep = "\t")
  cat(line, "\n", file = BAD_SAMPLES_TXT, append = TRUE)
}

# ===== Helper Functions =====

with_text_con <- function(path, FUN) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(con), add = TRUE)
  FUN(con)
}

readMM_strict <- function(path) {
  w <- character(0)
  mat <- withCallingHandlers(
    {
      with_text_con(path, function(con) Matrix::readMM(con))
    },
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  if (
    FAIL_ON_READMM_TRUNCATED &&
      any(grepl("expected .* entries but found only", w))
  ) {
    stop(
      "READMM_TRUNCATED: ",
      paste(w[grepl("expected .* entries but found only", w)], collapse = " | ")
    )
  }
  mat
}

readLines_safe <- function(path) {
  with_text_con(path, readLines)
}

readTable_safe <- function(path, ...) {
  with_text_con(path, function(con) read.table(con, ...))
}

drop_na_columns_if_any <- function(counts, gsm_id) {
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "CsparseMatrix")
  }
  if (anyNA(counts@x)) {
    cs <- Matrix::colSums(counts)
    bad <- which(is.na(cs))
    if (length(bad) > 0) {
      cat(
        "    [",
        gsm_id,
        "] WARNING: dropping NA columns: n=",
        length(bad),
        "\n",
        sep = ""
      )
      counts <- counts[, -bad, drop = FALSE]
    }
  }
  if (anyNA(counts@x)) {
    stop("COUNTS_CONTAINS_NA: counts@x still has NA after cleanup.")
  }
  counts
}

# ===== Gene ID Standardization Functions =====

is_ensembl_id <- function(gene_ids) {
  if (length(gene_ids) == 0) {
    return(FALSE)
  }

  ensg_pattern <- "^ENS[A-Z]*G\\d{11}(\\.\\d+)?$"
  match_rate <- sum(grepl(ensg_pattern, gene_ids)) / length(gene_ids)

  match_rate > 0.7
}

convert_ensg_to_symbol <- function(
  ensg_ids,
  provided_symbols = NULL,
  gsm_id = "unknown"
) {
  cat("    Converting ENSG IDs to gene symbols...\n")

  if (
    !is.null(provided_symbols) && length(provided_symbols) == length(ensg_ids)
  ) {
    cat("      Strategy: Using provided symbols from features file\n")

    valid_symbols <- provided_symbols != "" &
      !is.na(provided_symbols) &
      provided_symbols != "NA"

    result <- data.frame(
      original_id = ensg_ids,
      symbol = ifelse(valid_symbols, provided_symbols, ensg_ids),
      source = ifelse(valid_symbols, "features_tsv", "ensg_kept"),
      stringsAsFactors = FALSE
    )

    n_mapped <- sum(valid_symbols)
    mapping_rate <- n_mapped / length(ensg_ids)

    cat(
      "      Mapped:",
      n_mapped,
      "/",
      length(ensg_ids),
      "genes (",
      round(mapping_rate * 100, 1),
      "%)\n"
    )

    if (mapping_rate < MIN_MAPPING_RATE) {
      cat(
        "      ⚠️ Warning: Low mapping rate (<",
        MIN_MAPPING_RATE * 100,
        "%)\n"
      )
    }

    return(result)
  }

  cat("      Strategy: Using org.Hs.eg.db for conversion\n")

  ensg_clean <- sub("\\.\\d+$", "", ensg_ids)

  tryCatch(
    {
      symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = ensg_clean,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )

      result <- data.frame(
        original_id = ensg_ids,
        symbol = ifelse(is.na(symbols), ensg_ids, as.character(symbols)),
        source = ifelse(is.na(symbols), "ensg_unmapped", "org.Hs.eg.db"),
        stringsAsFactors = FALSE
      )

      n_mapped <- sum(!is.na(symbols))
      mapping_rate <- n_mapped / length(ensg_ids)

      cat(
        "      Mapped:",
        n_mapped,
        "/",
        length(ensg_ids),
        "genes (",
        round(mapping_rate * 100, 1),
        "%)\n"
      )

      if (mapping_rate < MIN_MAPPING_RATE) {
        cat(
          "      ⚠️ Warning: Low mapping rate (<",
          MIN_MAPPING_RATE * 100,
          "%)\n"
        )
      }

      return(result)
    },
    error = function(e) {
      cat("      ⚠️ org.Hs.eg.db conversion failed:", conditionMessage(e), "\n")
      cat("      Keeping original ENSG IDs\n")

      return(data.frame(
        original_id = ensg_ids,
        symbol = ensg_ids,
        source = "conversion_failed",
        stringsAsFactors = FALSE
      ))
    }
  )
}

dedup_sum_sparse <- function(counts, symbols) {
  stopifnot(length(symbols) == nrow(counts))

  f <- factor(symbols, levels = unique(symbols))

  i <- seq_len(length(f))
  j <- as.integer(f)
  P <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(length(f), nlevels(f))
  )

  agg <- Matrix::t(P) %*% counts

  rownames(agg) <- levels(f)

  agg
}

deduplicate_symbols <- function(counts, gene_mapping, method = "sum") {
  symbols <- gene_mapping$symbol

  dup_symbols <- symbols[duplicated(symbols)]
  if (length(dup_symbols) == 0) {
    cat("      No duplicate symbols found\n")
    rownames(counts) <- symbols
    return(list(
      counts = counts,
      mapping = gene_mapping,
      n_dup_symbols = 0
    ))
  }

  n_unique_dups <- length(unique(dup_symbols))
  cat("      Deduplicating", n_unique_dups, "symbols with duplicates...\n")

  if (method == "sum") {
    cat("      Method: Summing counts (sparse matrix operation)\n")

    counts_dedup <- dedup_sum_sparse(counts, symbols)

    mapping_dedup <- gene_mapping[!duplicated(gene_mapping$symbol), ]
    mapping_dedup$dedup_method <- "sum"
    mapping_dedup$n_duplicates <- as.integer(table(symbols)[
      mapping_dedup$symbol
    ])

    cat("      Result:", nrow(counts), "->", nrow(counts_dedup), "genes\n")
  } else if (method == "keep_first") {
    cat("      Method: Keeping first occurrence\n")

    keep_idx <- !duplicated(symbols)
    counts_dedup <- counts[keep_idx, , drop = FALSE]
    rownames(counts_dedup) <- symbols[keep_idx]

    mapping_dedup <- gene_mapping[keep_idx, ]
    mapping_dedup$dedup_method <- "keep_first"
    mapping_dedup$n_duplicates <- as.integer(table(symbols)[
      mapping_dedup$symbol
    ])

    cat("      Result:", nrow(counts), "->", nrow(counts_dedup), "genes\n")
  }

  list(
    counts = counts_dedup,
    mapping = mapping_dedup,
    n_dup_symbols = n_unique_dups
  )
}
# ===== Feature prefix cleanup helpers (NEW) =====

strip_grch38_prefix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- as.character(x)
  sub(GRCH38_PREFIX_PATTERN, "", x)
}

clean_feature_prefixes_pre <- function(
  counts,
  provided_symbols = NULL,
  gsm_id = "unknown"
) {
  stopifnot(!is.null(counts))
  rn_raw <- rownames(counts)
  if (is.null(rn_raw)) {
    return(list(
      counts = counts,
      provided_symbols = provided_symbols,
      n_grch38_stripped = 0L
    ))
  }

  n_stripped <- 0L
  if (STRIP_GRCH38_PREFIX) {
    rn_new <- strip_grch38_prefix(rn_raw)
    n_stripped <- sum(rn_new != rn_raw, na.rm = TRUE)
    if (n_stripped > 0) {
      cat(
        "    [",
        gsm_id,
        "] Stripping GRCh38- prefix: n=",
        n_stripped,
        "\n",
        sep = ""
      )
      rownames(counts) <- rn_new
    }
    if (
      !is.null(provided_symbols) && length(provided_symbols) == length(rn_raw)
    ) {
      provided_symbols <- strip_grch38_prefix(provided_symbols)
    }
  }

  list(
    counts = counts,
    provided_symbols = provided_symbols,
    n_grch38_stripped = as.integer(n_stripped)
  )
}

remove_sars_genes_post <- function(
  counts,
  gene_mapping = NULL,
  gsm_id = "unknown"
) {
  rn <- rownames(counts)
  if (is.null(rn)) {
    return(list(
      counts = counts,
      gene_mapping = gene_mapping,
      n_sars_removed = 0L
    ))
  }

  sars_idx <- grepl(SARS_PREFIX_PATTERN, rn)
  n_sars <- sum(sars_idx)

  if (n_sars > 0) {
    cat(
      "    [",
      gsm_id,
      "] Removing SARS-* genes after EmptyDrops: n=",
      n_sars,
      "\n",
      sep = ""
    )
    counts <- counts[!sars_idx, , drop = FALSE]

    if (!is.null(gene_mapping) && nrow(gene_mapping) > 0) {
      # gene_mapping$symbol 是最终用于 Seurat 的 rownames（你 dedup 后会用 symbol）
      if ("symbol" %in% names(gene_mapping)) {
        gene_mapping <- gene_mapping[
          !grepl(SARS_PREFIX_PATTERN, gene_mapping$symbol),
          ,
          drop = FALSE
        ]
      } else if ("original_id" %in% names(gene_mapping)) {
        gene_mapping <- gene_mapping[
          !grepl(SARS_PREFIX_PATTERN, gene_mapping$original_id),
          ,
          drop = FALSE
        ]
      }
    }
  }

  list(
    counts = counts,
    gene_mapping = gene_mapping,
    n_sars_removed = as.integer(n_sars)
  )
}

standardize_gene_ids <- function(
  counts,
  gsm_id,
  provided_symbols = NULL,
  keep_unmapped_for_now = TRUE,
  dedup_method = SYMBOL_DEDUP_METHOD
) {
  if (!STANDARDIZE_GENE_IDS) {
    cat("    Gene ID standardization: DISABLED\n")
    return(list(
      counts = counts,
      mapping = data.frame(
        original_id = rownames(counts),
        symbol = rownames(counts),
        source = "standardization_disabled",
        stringsAsFactors = FALSE
      ),
      standardized = FALSE,
      n_dup_symbols = 0
    ))
  }

  cat("    Gene ID standardization: ENABLED\n")

  original_genes <- rownames(counts)

  if (!is_ensembl_id(original_genes)) {
    cat("      Input format: Gene symbols (not ENSG)\n")

    dup_genes <- original_genes[duplicated(original_genes)]
    if (length(dup_genes) > 0) {
      n_unique_dups <- length(unique(dup_genes))
      cat("      Found", n_unique_dups, "duplicate gene symbols\n")

      gene_mapping <- data.frame(
        original_id = original_genes,
        symbol = original_genes,
        source = "symbol_input",
        stringsAsFactors = FALSE
      )

      result <- deduplicate_symbols(counts, gene_mapping, dedup_method)
      result$standardized <- TRUE
      return(result)
    }

    cat("      No duplicates found, keeping original symbols\n")
    return(list(
      counts = counts,
      mapping = data.frame(
        original_id = original_genes,
        symbol = original_genes,
        source = "symbol_input",
        stringsAsFactors = FALSE
      ),
      standardized = FALSE,
      n_dup_symbols = 0
    ))
  }

  cat("      Input format: ENSG IDs detected\n")

  gene_mapping <- convert_ensg_to_symbol(
    original_genes,
    provided_symbols = provided_symbols,
    gsm_id = gsm_id
  )

  if (!keep_unmapped_for_now) {
    unmapped_idx <- gene_mapping$source %in% c("ensg_unmapped", "ensg_kept")
    n_unmapped <- sum(unmapped_idx)

    if (n_unmapped > 0) {
      cat("      Removing", n_unmapped, "unmapped genes\n")
      counts <- counts[!unmapped_idx, , drop = FALSE]
      gene_mapping <- gene_mapping[!unmapped_idx, ]
    }
  } else if (!KEEP_UNMAPPED_GENES) {
    cat("      Unmapped genes kept for EmptyDrops (will filter later)\n")
  }

  result <- deduplicate_symbols(counts, gene_mapping, dedup_method)
  result$standardized <- TRUE

  return(result)
}

filter_unmapped_genes <- function(counts, gene_mapping) {
  if (KEEP_UNMAPPED_GENES) {
    return(list(counts = counts, mapping = gene_mapping, n_removed = 0))
  }

  unmapped_idx <- gene_mapping$source %in% c("ensg_unmapped", "ensg_kept")
  n_unmapped <- sum(unmapped_idx)

  if (n_unmapped > 0) {
    cat("    Post-EmptyDrops: Removing", n_unmapped, "unmapped genes\n")
    counts <- counts[!unmapped_idx, , drop = FALSE]
    gene_mapping <- gene_mapping[!unmapped_idx, ]
  }

  list(counts = counts, mapping = gene_mapping, n_removed = n_unmapped)
}

# ===== Data Loading Functions =====

extract_10x_prefix <- function(gsm_dir) {
  files <- list.files(gsm_dir, full.names = FALSE)

  barcode_file <- grep(
    "barcodes\\.tsv(\\.gz)?$",
    files,
    value = TRUE,
    ignore.case = TRUE
  )[1]

  if (is.na(barcode_file)) {
    return(NULL)
  }

  if (grepl("^barcodes\\.tsv(\\.gz)?$", barcode_file, ignore.case = TRUE)) {
    return(NULL)
  }

  prefix <- sub("barcodes\\.tsv(\\.gz)?$", "", barcode_file, ignore.case = TRUE)

  prefix <- sub("[._-]$", "", prefix)

  if (nchar(prefix) == 0) NULL else prefix
}
# ===== Regex Escape Helper (FIX for {} / special chars) =====
escape_regex <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) {
    return(NULL)
  }

  # Escape regex metacharacters using PCRE (perl=TRUE)
  # This is robust against {}, [], \, ., ?, +, *, ^, $, (, ), | ...
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

# ===== Fixed 10x manual reader =====
read_10x_manual <- function(gsm_dir, sample_prefix = NULL) {
  files <- list.files(gsm_dir, full.names = FALSE)

  pick1 <- function(pattern) {
    out <- grep(pattern, files, value = TRUE, ignore.case = TRUE)
    if (length(out) == 0) {
      return(NA_character_)
    }
    out[1]
  }

  if (is.null(sample_prefix)) {
    sample_prefix <- extract_10x_prefix(gsm_dir)
  }

  prefix_esc <- escape_regex(sample_prefix)

  barcode_pattern <- if (is.null(prefix_esc)) {
    "barcodes\\.tsv(\\.gz)?$"
  } else {
    paste0(prefix_esc, ".*barcodes\\.tsv(\\.gz)?$")
  }

  feature_pattern <- if (is.null(prefix_esc)) {
    "(features|genes)\\.tsv(\\.gz)?$"
  } else {
    paste0(prefix_esc, ".*(features|genes)\\.tsv(\\.gz)?$")
  }

  matrix_pattern <- if (is.null(prefix_esc)) {
    "matrix\\.mtx(\\.gz)?$"
  } else {
    paste0(prefix_esc, ".*matrix\\.mtx(\\.gz)?$")
  }

  barcode_file <- pick1(barcode_pattern)
  feature_file <- pick1(feature_pattern)
  matrix_file <- pick1(matrix_pattern)

  if (is.na(barcode_file) || is.na(feature_file) || is.na(matrix_file)) {
    return(NULL)
  }

  barcode_path <- file.path(gsm_dir, barcode_file)
  feature_path <- file.path(gsm_dir, feature_file)
  matrix_path <- file.path(gsm_dir, matrix_file)

  cat("    10x(manual) files:\n")
  cat("      barcodes:", barcode_file, "\n")
  cat("      features:", feature_file, "\n")
  cat("      matrix  :", matrix_file, "\n")

  barcodes <- readLines_safe(barcode_path)

  features <- readTable_safe(
    feature_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  mat <- readMM_strict(matrix_path)
  mat <- as(mat, "CsparseMatrix")

  if (ncol(mat) != length(barcodes)) {
    stop(
      "10X_DIM_MISMATCH: ncol(mtx)=",
      ncol(mat),
      " != length(barcodes)=",
      length(barcodes)
    )
  }
  if (nrow(mat) != nrow(features)) {
    stop(
      "10X_DIM_MISMATCH: nrow(mtx)=",
      nrow(mat),
      " != nrow(features)=",
      nrow(features)
    )
  }

  colnames(mat) <- barcodes
  rownames(mat) <- features[, 1]

  gene_symbols <- if (ncol(features) >= 2) {
    features[, 2]
  } else {
    NULL
  }

  list(
    counts = mat,
    gene_ids = features[, 1],
    gene_symbols = gene_symbols
  )
}

detect_data_format <- function(sample_dir, gsm_id) {
  files <- list.files(sample_dir, full.names = FALSE)

  h5_files <- grep("\\.(h5|hdf5)$", files, value = TRUE, ignore.case = TRUE)
  if (length(h5_files) > 0) {
    return(list(
      format = "h5",
      path = file.path(sample_dir, h5_files[1])
    ))
  }

  has_barcodes <- any(grepl(
    "barcodes\\.tsv(\\.gz)?$",
    files,
    ignore.case = TRUE
  ))
  has_features <- any(grepl(
    "(features|genes)\\.tsv(\\.gz)?$",
    files,
    ignore.case = TRUE
  ))
  has_matrix <- any(grepl("matrix\\.mtx(\\.gz)?$", files, ignore.case = TRUE))

  if (has_barcodes && has_features && has_matrix) {
    return(list(
      format = "10x",
      path = sample_dir
    ))
  }

  matrix_files <- grep(
    "\\.(txt|tsv|csv|gz)$",
    files,
    value = TRUE,
    ignore.case = TRUE
  )
  matrix_files <- matrix_files[
    !grepl(
      "features|barcodes|genes|meta|anno",
      matrix_files,
      ignore.case = TRUE
    )
  ]

  if (length(matrix_files) > 0) {
    prio_idx <- grep(
      "expr|expression|count|matrix",
      matrix_files,
      ignore.case = TRUE
    )
    chosen <- if (length(prio_idx) > 0) {
      matrix_files[prio_idx[1]]
    } else {
      matrix_files[1]
    }
    return(list(
      format = "matrix",
      path = file.path(sample_dir, chosen)
    ))
  }

  list(format = "unknown", path = NULL)
}

read_h5_with_ensg <- function(h5_path) {
  counts_list <- Read10X_h5(h5_path, use.names = FALSE)

  if (is.list(counts_list) && !inherits(counts_list, "Matrix")) {
    if ("Gene Expression" %in% names(counts_list)) {
      counts <- counts_list[["Gene Expression"]]
    } else {
      counts <- counts_list[[1]]
    }
  } else {
    counts <- counts_list
  }

  tryCatch(
    {
      h5 <- H5File$new(h5_path, mode = "r")

      gene_names <- NULL

      if (h5$exists("matrix/features/name")) {
        gene_names <- h5[["matrix/features/name"]][]
      } else if (h5$exists("matrix/features/gene_name")) {
        gene_names <- h5[["matrix/features/gene_name"]][]
      } else if (h5$exists("gene_names")) {
        gene_names <- h5[["gene_names"]][]
      }

      h5$close_all()

      if (!is.null(gene_names) && length(gene_names) == nrow(counts)) {
        return(list(counts = counts, gene_symbols = gene_names))
      }
    },
    error = function(e) {}
  )

  list(counts = counts, gene_symbols = NULL)
}

load_counts_data <- function(format_info, gsm_id) {
  format_type <- format_info$format
  path <- format_info$path

  tryCatch(
    {
      counts <- NULL
      provided_symbols <- NULL

      if (format_type == "h5") {
        cat("    Loading H5 format...\n")

        h5_result <- read_h5_with_ensg(path)
        counts <- h5_result$counts
        provided_symbols <- h5_result$gene_symbols

        if (!is.null(provided_symbols)) {
          cat("    Extracted gene symbols from H5 file\n")
        } else {
          cat("    Will use org.Hs.eg.db for ENSG conversion\n")
        }
      } else if (format_type == "10x") {
        cat("    Loading 10x format...\n")

        tryCatch(
          {
            counts_list <- Read10X(path, gene.column = 1)

            if (is.list(counts_list) && !inherits(counts_list, "Matrix")) {
              if ("Gene Expression" %in% names(counts_list)) {
                counts <- counts_list[["Gene Expression"]]
              } else {
                counts <- counts_list[[1]]
              }
            } else {
              counts <- counts_list
            }

            feature_file <- list.files(
              path,
              pattern = "(features|genes)\\.tsv(\\.gz)?$",
              full.names = TRUE,
              ignore.case = TRUE
            )[1]

            if (!is.na(feature_file) && file.exists(feature_file)) {
              features <- readTable_safe(
                feature_file,
                sep = "\t",
                header = FALSE,
                stringsAsFactors = FALSE,
                quote = "",
                comment.char = ""
              )
              if (ncol(features) >= 2 && nrow(features) == nrow(counts)) {
                provided_symbols <- features[, 2]
                cat("    Gene symbols from features file\n")
              }
            }
          },
          error = function(e) {
            cat("    Seurat Read10X failed, trying manual reader...\n")
            result <- read_10x_manual(path, sample_prefix = NULL)

            if (is.null(result)) {
              stop("10X_MANUAL_READ_FAILED: ", conditionMessage(e))
            }

            counts <<- result$counts
            provided_symbols <<- result$gene_symbols
          }
        )
      } else if (format_type == "matrix") {
        cat("    Loading matrix format...\n")

        con <- with_text_con(path, function(con) {
          line1 <- readLines(con, n = 1)
          line2 <- readLines(con, n = 1)
          list(line1 = line1, line2 = line2)
        })

        line1 <- con$line1
        line2 <- con$line2

        delim <- if (grepl("\t", line1)) "\t" else ","
        tok1 <- strsplit(line1, split = delim, fixed = TRUE)[[1]]
        tok2 <- strsplit(line2, split = delim, fixed = TRUE)[[1]]

        missing_gene_placeholder <- (length(tok2) == length(tok1) + 1)

        if (length(tok1) > MAX_WIDE_MATRIX_CELLS) {
          stop(
            "WIDE_MATRIX_TOO_MANY_CELLS: n_cells=",
            length(tok1),
            " > MAX_WIDE_MATRIX_CELLS=",
            MAX_WIDE_MATRIX_CELLS
          )
        }

        if (missing_gene_placeholder) {
          barcodes <- tok1
          mat <- fread(
            path,
            sep = delim,
            header = FALSE,
            skip = 1,
            data.table = FALSE
          )
          if (ncol(mat) != (length(barcodes) + 1)) {
            stop("WIDE_MATRIX_COL_MISMATCH")
          }
          colnames(mat) <- c("gene", barcodes)
        } else {
          tokens <- tok1
          use_header <- FALSE
          if (length(tokens) > 1) {
            suppressWarnings({
              num_flags <- !is.na(as.numeric(tokens[-1]))
            })
            use_header <- !(all(num_flags))
          }
          mat <- fread(
            path,
            sep = delim,
            header = use_header,
            data.table = FALSE
          )

          if (use_header && ncol(mat) > 1) {
            cn <- colnames(mat)
            suppressWarnings({
              cn_numeric <- !is.na(as.numeric(cn[-1]))
            })
            if (all(cn_numeric)) {
              mat <- fread(
                path,
                sep = delim,
                header = FALSE,
                data.table = FALSE
              )
            }
          }
        }

        genes <- mat[, 1]
        mat2 <- mat[, -1, drop = FALSE]
        if (ncol(mat2) == 0) {
          stop("MATRIX_NO_COUNT_COLUMNS")
        }

        mat_numeric <- as.matrix(mat2)
        rownames(mat_numeric) <- genes

        if (is.null(colnames(mat2)) || any(colnames(mat2) == "")) {
          colnames(mat_numeric) <- paste0(
            gsm_id,
            "_Cell",
            seq_len(ncol(mat_numeric))
          )
        } else {
          colnames(mat_numeric) <- colnames(mat2)
        }

        counts <- as(mat_numeric, "CsparseMatrix")
        provided_symbols <- NULL
      } else {
        stop("UNKNOWN_FORMAT")
      }

      if (!inherits(counts, "dgCMatrix")) {
        counts <- as(counts, "CsparseMatrix")
      }
      counts <- drop_na_columns_if_any(counts, gsm_id)
      # >>> NEW: strip GRCh38- prefix BEFORE standardization
      pre_clean <- clean_feature_prefixes_pre(counts, provided_symbols, gsm_id)
      counts <- pre_clean$counts
      provided_symbols <- pre_clean$provided_symbols
      n_grch38_stripped <- pre_clean$n_grch38_stripped
      # <<< NEW
      cat("    Loaded:", nrow(counts), "genes x", ncol(counts), "cells\n")

      std_result <- standardize_gene_ids(
        counts,
        gsm_id,
        provided_symbols = provided_symbols,
        keep_unmapped_for_now = TRUE
      )

      counts_std <- std_result$counts
      gene_mapping <- std_result$mapping
      was_standardized <- std_result$standardized
      n_dup_symbols <- std_result$n_dup_symbols

      if (was_standardized) {
        cat(
          "    Post-standardization:",
          nrow(counts_std),
          "genes x",
          ncol(counts_std),
          "cells\n"
        )
      }

      list(
        counts = counts_std,
        gene_mapping = gene_mapping,
        standardized = was_standardized,
        n_dup_symbols = n_dup_symbols,
        n_grch38_stripped = n_grch38_stripped, # <<< NEW
        reason = NULL
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      cat("    ❌ Error loading", gsm_id, ":", msg, "\n")
      list(
        counts = NULL,
        gene_mapping = NULL,
        standardized = FALSE,
        n_dup_symbols = 0,
        n_grch38_stripped = 0, # <<< NEW
        reason = msg
      )
    }
  )
}

run_emptydrops_filter <- function(counts, gsm_id) {
  cells_before <- ncol(counts)

  if (cells_before < 100) {
    cat("    Skipping EmptyDrops (too few cells:", cells_before, ")\n")
    return(list(
      counts = counts,
      cells_before = cells_before,
      cells_after = cells_before,
      cells_removed = 0,
      method = "skipped_fewcells",
      reason = NULL
    ))
  }

  tryCatch(
    {
      counts <- drop_na_columns_if_any(counts, gsm_id)

      cat(
        "    Running EmptyDrops (lower=",
        EMPTYDROPS_LOWER,
        ", FDR<=",
        EMPTYDROPS_FDR_THRESHOLD,
        ")...\n",
        sep = ""
      )
      set.seed(42)
      e.out <- emptyDrops(counts, lower = EMPTYDROPS_LOWER)

      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "C",
          location = "emptydrops:after_call",
          message = "EmptyDrops completed",
          data = list(
            gsm_id = gsm_id,
            e_out_class = class(e.out),
            e_out_names = names(e.out),
            has_fdr = "FDR" %in% names(e.out),
            fdr_na_count = if ("FDR" %in% names(e.out)) {
              sum(is.na(e.out$FDR))
            } else {
              NA
            }
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      is.cell <- e.out$FDR <= EMPTYDROPS_FDR_THRESHOLD
      is.cell[is.na(is.cell)] <- FALSE

      high.count <- Matrix::colSums(counts) > EMPTYDROPS_LOWER
      is.cell <- is.cell | (is.na(e.out$FDR) & high.count)

      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "C",
          location = "emptydrops:after_filter",
          message = "Cell filtering completed",
          data = list(
            gsm_id = gsm_id,
            n_cells_before = length(is.cell),
            n_cells_after = sum(is.cell),
            n_na_fdr = if ("FDR" %in% names(e.out)) {
              sum(is.na(e.out$FDR))
            } else {
              NA
            }
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      counts_filtered <- counts[, is.cell, drop = FALSE]
      cells_after <- ncol(counts_filtered)
      cells_removed <- cells_before - cells_after

      cat("    EmptyDrops Results:\n")
      cat("      Cells before:", format(cells_before, big.mark = ","), "\n")
      cat("      Cells after :", format(cells_after, big.mark = ","), "\n")
      cat(
        "      Removed     :",
        format(cells_removed, big.mark = ","),
        " (",
        if (cells_before > 0) {
          round(cells_removed / cells_before * 100, 1)
        } else {
          0
        },
        "%)\n",
        sep = ""
      )

      list(
        counts = counts_filtered,
        cells_before = cells_before,
        cells_after = cells_after,
        cells_removed = cells_removed,
        method = "emptydrops",
        reason = NULL
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      cat("    ⚠️  EmptyDrops failed:", msg, "\n")
      cat("    Proceeding without EmptyDrops filtering\n")
      list(
        counts = counts,
        cells_before = cells_before,
        cells_after = cells_before,
        cells_removed = 0,
        method = "failed",
        reason = msg
      )
    }
  )
}

extract_orig_ident <- function(dataset_id) {
  dataset_id
}

calculate_qc_metrics <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^(MT-|mt-)"
  )
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^(RP[SL]|rp[sl])"
  )

  if (median(seurat_obj$percent.mt, na.rm = TRUE) < 0.1) {
    cat("    ⚠️ Warning: Very low MT%. Check gene ID conversion.\n")
  }

  seurat_obj
}

pick_group_col <- function(md, gsm_id) {
  for (k in c("sample", "sample_id", "orig.ident", "dataset")) {
    if (k %in% colnames(md) && length(unique(md[[k]])) > 1) {
      return(k)
    }
  }
  return("gsm_id_display")
}

generate_qc_plot <- function(
  seurat_obj,
  gsm_id,
  output_path,
  max_points = MAX_PLOT_POINTS
) {
  qc_data <- as.data.table(seurat_obj@meta.data)

  grp_col <- pick_group_col(seurat_obj@meta.data, gsm_id)
  if (grp_col == "gsm_id_display") {
    qc_data[, group := gsm_id]
  } else {
    qc_data[, group := get(grp_col)]
  }

  if (nrow(qc_data) > max_points) {
    sample_idx <- sample.int(nrow(qc_data), max_points)
    qc_data_sample <- qc_data[sample_idx]
  } else {
    qc_data_sample <- qc_data
  }

  p1 <- ggplot(qc_data, aes(x = group, y = nFeature_RNA)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nFeature_min,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nFeature_max,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    labs(title = "Features per Cell", x = "", y = "nFeature_RNA") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p2 <- ggplot(qc_data, aes(x = group, y = nCount_RNA)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nCount_min,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nCount_max,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    scale_y_log10() +
    labs(title = "Counts per Cell", x = "", y = "nCount_RNA (log10)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p3 <- ggplot(qc_data, aes(x = group, y = percent.mt)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$percent_mt_max,
      linetype = "dashed",
      color = "darkred",
      linewidth = 0.5
    ) +
    labs(title = "Mitochondrial %", x = "", y = "% MT genes") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p4 <- ggplot(qc_data, aes(x = group, y = percent.rb)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$percent_rb_min,
      linetype = "dashed",
      color = "darkviolet",
      linewidth = 0.5
    ) +
    labs(title = "Ribosomal %", x = "", y = "% RB genes") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p5 <- ggplot(qc_data_sample, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.8) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Feature vs Count",
      x = "nCount_RNA (log10)",
      y = "nFeature_RNA (log10)"
    ) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p6 <- ggplot(qc_data_sample, aes(x = nCount_RNA, y = percent.mt)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_x_log10() +
    labs(title = "MT% vs Count", x = "nCount_RNA (log10)", y = "% MT genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  layout <- "ABC\nDEF"
  combined <- p1 +
    p2 +
    p3 +
    p4 +
    p5 +
    p6 +
    plot_layout(design = layout) +
    plot_annotation(
      title = paste0("Quality Control Metrics: ", gsm_id),
      subtitle = paste0(
        "Cells: ",
        format(ncol(seurat_obj), big.mark = ","),
        " | Genes: ",
        format(nrow(seurat_obj), big.mark = ","),
        " | Median nFeature: ",
        format(
          round(median(qc_data$nFeature_RNA, na.rm = TRUE)),
          big.mark = ","
        ),
        " | Median nCount: ",
        format(round(median(qc_data$nCount_RNA, na.rm = TRUE)), big.mark = ","),
        if (nrow(qc_data) > max_points) {
          paste0(" | Sampled ", format(max_points, big.mark = ","), " cells")
        } else {
          ""
        }
      ),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    )

  ggsave(
    filename = output_path,
    plot = combined,
    width = 12,
    height = 8,
    dpi = 300
  )
  cat("    QC plot saved:", basename(output_path), "\n")
}

# ===== Load Metadata =====
cat("Loading metadata...\n")
if (!file.exists(METADATA_PATH)) {
  stop("Metadata file not found: ", METADATA_PATH)
}
metadata <- fread(METADATA_PATH, header = TRUE)

cat("Standardizing metadata columns...\n")
cols <- names(metadata)

if ("filtered_or_Falset" %in% cols && !"filtered_or_raw" %in% cols) {
  setnames(metadata, "filtered_or_Falset", "filtered_or_raw")
}

ref_cols <- grep("^reference", cols, value = TRUE, ignore.case = TRUE)
if (length(ref_cols) == 1 && !"reference_genome" %in% cols) {
  setnames(metadata, ref_cols, "reference_genome")
}

# P0 FIX: Robust GSM column detection with NA handling
infer_gsm_col <- function(dt) {
  # #region agent log
  log_entry <- jsonlite::toJSON(
    list(
      sessionId = "debug-session",
      runId = "run1",
      hypothesisId = "B",
      location = "infer_gsm_col:entry",
      message = "infer_gsm_col called",
      data = list(n_cols = length(names(dt)), n_rows = nrow(dt)),
      timestamp = as.numeric(Sys.time()) * 1000
    ),
    auto_unbox = TRUE
  )
  cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion

  for (k in names(dt)) {
    v <- dt[[k]]
    if (is.character(v) || is.factor(v)) {
      v_char <- as.character(v)
      # Replace NA with empty string to prevent grepl from returning NA
      v_char[is.na(v_char)] <- ""
      # Calculate match rate with na.rm=TRUE for safety
      match_rate <- mean(grepl("^(GSM|HRR)\\d+$", v_char), na.rm = TRUE)

      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "B",
          location = "infer_gsm_col:loop",
          message = "Column match rate calculated",
          data = list(
            column = k,
            match_rate = match_rate,
            n_na = sum(is.na(v_char)),
            n_total = length(v_char)
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      # Explicit NA check before comparison
      if (!is.na(match_rate) && match_rate > 0.8) {
        # #region agent log
        log_entry <- jsonlite::toJSON(
          list(
            sessionId = "debug-session",
            runId = "run1",
            hypothesisId = "B",
            location = "infer_gsm_col:return",
            message = "GSM column found",
            data = list(column = k, match_rate = match_rate),
            timestamp = as.numeric(Sys.time()) * 1000
          ),
          auto_unbox = TRUE
        )
        cat(
          log_entry,
          "\n",
          file = "/home/h2048/.cursor/debug.log",
          append = TRUE
        )
        # #endregion
        return(k)
      }
    }
  }

  # #region agent log
  log_entry <- jsonlite::toJSON(
    list(
      sessionId = "debug-session",
      runId = "run1",
      hypothesisId = "B",
      location = "infer_gsm_col:return_null",
      message = "No GSM column found",
      data = list(),
      timestamp = as.numeric(Sys.time()) * 1000
    ),
    auto_unbox = TRUE
  )
  cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion

  NULL
}

gsm_col <- infer_gsm_col(metadata)
if (!is.null(gsm_col)) {
  cat("  Detected GSM column:", gsm_col, "\n")
  metadata[, gsm_id := get(gsm_col)]

  # ROBUSTNESS FIX: Clean whitespace from gsm_id
  metadata[, gsm_id := gsub("[[:space:]]+", "", gsm_id)]

  setkey(metadata, gsm_id)

  # ROBUSTNESS FIX: Check for duplicate gsm_id (fail-fast)
  if (anyDuplicated(metadata$gsm_id)) {
    dup_ids <- metadata$gsm_id[duplicated(metadata$gsm_id)]
    stop(
      "Duplicated gsm_id detected in metadata: ",
      paste(head(dup_ids, 5), collapse = ", "),
      "\nPlease de-duplicate before running."
    )
  }
} else {
  warning("Cannot infer GSM/HRR column in metadata; metadata matching may fail")
  if ("sample" %in% names(metadata)) {
    metadata[, gsm_id := sample]
    metadata[, gsm_id := gsub("[[:space:]]+", "", gsm_id)]
    setkey(metadata, gsm_id)
    cat("  Fallback: Using 'sample' column as gsm_id\n")
  }
}

cat("  Total metadata records:", nrow(metadata), "\n")
if ("gsm_id" %in% names(metadata)) {
  cat("  Unique gsm_ids:", uniqueN(metadata$gsm_id), "\n")
}
cat("\n")

# ===== Dynamic GSE/HRA Detection =====
cat("Scanning for GSE/HRA directories...\n")
all_dirs <- list.dirs(DATA_ROOT, recursive = FALSE, full.names = FALSE)
DATASET_IDS <- all_dirs[grepl("^(GSE|HRA)\\d+$", all_dirs)]
if (length(DATASET_IDS) == 0) {
  stop("No GSE/HRA directories found in: ", DATA_ROOT)
}

cat("  Found", length(DATASET_IDS), "datasets:\n")
cat("  ", paste(DATASET_IDS, collapse = ", "), "\n\n")

# ===== Main Processing Loop =====
processing_summary_list <- list()
summary_counter <- 1

total_processed <- 0
total_failed <- 0
total_cells_raw <- 0
total_cells_filtered <- 0

for (dataset_id in DATASET_IDS) {
  dataset_dir <- file.path(DATA_ROOT, dataset_id)
  if (!dir.exists(dataset_dir)) {
    cat("⚠️  Skipping", dataset_id, "(directory not found)\n\n")
    next
  }

  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Processing Dataset:", dataset_id, "\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

  sample_dirs <- list.dirs(dataset_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("^(GSM|HRR)\\d+$", basename(sample_dirs))]

  if (length(sample_dirs) == 0) {
    cat("⚠️  No sample directories found in", dataset_id, "\n\n")
    next
  }

  cat("Found", length(sample_dirs), "samples\n\n")

  for (sample_dir in sample_dirs) {
    gsm_id <- basename(sample_dir)
    cat("📦", gsm_id, "\n")

    format_info <- detect_data_format(sample_dir, gsm_id)
    cat("  Format:", format_info$format, "\n")

    if (format_info$format == "unknown") {
      cat("  ⚠️  Skipping (unknown format)\n\n")
      append_bad_sample(
        dataset_id,
        gsm_id,
        "unknown",
        "unknown_format",
        NA_character_
      )
      processing_summary_list[[summary_counter]] <- data.table(
        Dataset = dataset_id,
        GSM_ID = gsm_id,
        Sample = NA_character_,
        Format = "unknown",
        Cells_Raw = NA_integer_,
        Cells_Filtered = NA_integer_,
        Cells_Removed = NA_integer_,
        EmptyDrops_Method = "skipped",
        Genes_Original = NA_integer_,
        Genes_Final = NA_integer_,
        Genes_Unmapped_Removed = NA_integer_,
        N_Dup_Symbols = NA_integer_,
        Dedup_Method = NA_character_,
        Gene_Standardized = FALSE,
        Median_nFeature = NA_real_,
        Median_nCount = NA_real_,
        Median_MT = NA_real_,
        Median_RB = NA_real_,
        Status = "skipped",
        Reason = "unknown_format"
      )
      summary_counter <- summary_counter + 1
      total_failed <- total_failed + 1
      next
    }

    ld <- load_counts_data(format_info, gsm_id)
    counts <- ld$counts
    gene_mapping <- ld$gene_mapping
    was_standardized <- ld$standardized
    n_dup_symbols <- ld$n_dup_symbols
    n_grch38_stripped <- ld$n_grch38_stripped

    if (is.null(counts)) {
      append_bad_sample(
        dataset_id,
        gsm_id,
        format_info$format,
        ld$reason,
        format_info$path
      )
      processing_summary_list[[summary_counter]] <- data.table(
        Dataset = dataset_id,
        GSM_ID = gsm_id,
        Sample = NA_character_,
        Format = format_info$format,
        Cells_Raw = NA_integer_,
        Cells_Filtered = NA_integer_,
        Cells_Removed = NA_integer_,
        EmptyDrops_Method = "failed_load",
        Genes_Original = NA_integer_,
        Genes_Final = NA_integer_,
        Genes_Unmapped_Removed = NA_integer_,
        N_Dup_Symbols = NA_integer_,
        Dedup_Method = NA_character_,
        Gene_Standardized = FALSE,
        Median_nFeature = NA_real_,
        Median_nCount = NA_real_,
        Median_MT = NA_real_,
        Median_RB = NA_real_,
        Status = "failed",
        Reason = ld$reason
      )
      summary_counter <- summary_counter + 1
      total_failed <- total_failed + 1
      cat("\n")
      next
    }

    genes_original <- if (!is.null(gene_mapping)) {
      nrow(gene_mapping)
    } else {
      nrow(counts)
    }

    if (RUN_EMPTYDROPS) {
      ed_result <- run_emptydrops_filter(counts, gsm_id)
      counts <- ed_result$counts
      cells_raw <- ed_result$cells_before
      cells_filtered <- ed_result$cells_after
      cells_removed <- ed_result$cells_removed
      ed_method <- ed_result$method
      ed_reason <- ed_result$reason
      if (!is.null(ed_reason)) {
        append_bad_sample(
          dataset_id,
          gsm_id,
          format_info$format,
          paste0("emptydrops_", ed_method, ":", ed_reason),
          format_info$path
        )
      }
    } else {
      cells_raw <- ncol(counts)
      cells_filtered <- ncol(counts)
      cells_removed <- 0
      ed_method <- "disabled"
      ed_reason <- NULL
    }

    genes_unmapped_removed <- 0
    if (!is.null(gene_mapping)) {
      filter_result <- filter_unmapped_genes(counts, gene_mapping)
      counts <- filter_result$counts
      gene_mapping <- filter_result$mapping
      genes_unmapped_removed <- filter_result$n_removed
    }
    # >>> NEW: remove SARS-* genes AFTER EmptyDrops (recommended)
    genes_sars_removed <- 0L
    if (REMOVE_SARS_GENES) {
      post_clean <- remove_sars_genes_post(counts, gene_mapping, gsm_id)
      counts <- post_clean$counts
      gene_mapping <- post_clean$gene_mapping
      genes_sars_removed <- post_clean$n_sars_removed
    }
    # <<< NEW

    genes_final <- nrow(counts)

    if (ncol(counts) == 0 || nrow(counts) == 0) {
      cat("  ⚠️  No valid cells or genes. Skipping.\n\n")
      append_bad_sample(
        dataset_id,
        gsm_id,
        format_info$format,
        "failed_nocells",
        format_info$path
      )

      bio_sample <- if (
        "gsm_id" %in% names(metadata) && gsm_id %in% metadata$gsm_id
      ) {
        meta_row_tmp <- metadata[.(gsm_id)]
        if ("sample" %in% names(meta_row_tmp)) {
          meta_row_tmp$sample
        } else {
          NA_character_
        }
      } else {
        NA_character_
      }

      processing_summary_list[[summary_counter]] <- data.table(
        Dataset = dataset_id,
        GSM_ID = gsm_id,
        Sample = bio_sample,
        Format = format_info$format,
        Cells_Raw = cells_raw,
        Cells_Filtered = 0,
        Cells_Removed = cells_raw,
        EmptyDrops_Method = if (RUN_EMPTYDROPS) {
          "failed_nocells"
        } else {
          "disabled_nocells"
        },
        Genes_Original = genes_original,
        Genes_Final = 0,
        Genes_Unmapped_Removed = genes_unmapped_removed,
        N_Dup_Symbols = n_dup_symbols,
        Dedup_Method = SYMBOL_DEDUP_METHOD,
        Gene_Standardized = was_standardized,
        Median_nFeature = NA_real_,
        Median_nCount = NA_real_,
        Median_MT = NA_real_,
        Median_RB = NA_real_,
        Status = "failed_nocells",
        Reason = "no_cells_or_genes_after_filter"
      )
      summary_counter <- summary_counter + 1
      total_failed <- total_failed + 1
      rm(counts)
      gc(FALSE)
      next
    }

    orig_ident <- extract_orig_ident(dataset_id)
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = orig_ident,
      min.cells = 0,
      min.features = 0
    )

    cat(
      "  Created Seurat object:",
      format(ncol(seurat_obj), big.mark = ","),
      "cells x",
      format(nrow(seurat_obj), big.mark = ","),
      "genes\n"
    )

    seurat_obj$gsm_id <- gsm_id
    seurat_obj$dataset <- dataset_id
    seurat_obj$orig.ident <- orig_ident
    seurat_obj$emptydrops_filtered <- RUN_EMPTYDROPS &&
      ed_method == "emptydrops"

    if (!is.null(gene_mapping)) {
      mapping_file <- file.path(
        OUTPUT_DIR,
        paste0(gsm_id, "_gene_mapping.csv.gz")
      )
      fwrite(gene_mapping, mapping_file)

      seurat_obj@misc$gene_standardization <- list(
        standardized = was_standardized,
        method = if (was_standardized) {
          unique(gene_mapping$source)[1]
        } else {
          "none"
        },
        dedup_method = SYMBOL_DEDUP_METHOD,
        keep_unmapped = KEEP_UNMAPPED_GENES,
        genes_original = genes_original,
        genes_final = genes_final,
        genes_unmapped_removed = genes_unmapped_removed,
        n_dup_symbols = n_dup_symbols,
        mapping_file = basename(mapping_file),
        timestamp = Sys.time()
      )
    }

    if ("gsm_id" %in% names(metadata) && gsm_id %in% metadata$gsm_id) {
      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "A",
          location = "metadata_access:before",
          message = "About to access metadata",
          data = list(gsm_id = gsm_id, metadata_cols = names(metadata)),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      meta_row <- metadata[.(gsm_id)]

      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "A",
          location = "metadata_access:after_lookup",
          message = "Metadata row retrieved",
          data = list(
            gsm_id = gsm_id,
            meta_row_cols = names(meta_row),
            n_rows = nrow(meta_row)
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      required_cols <- c(
        "sample",
        "GEO",
        "tissue",
        "tissue_level_2",
        "disease_level_1",
        "disease_level_2",
        "condition",
        "tissue_sampling_method",
        "filtered_or_raw",
        "target_cell",
        "frozen_or_fresh",
        "sex",
        "age",
        "platform",
        "reference_genome",
        "assay",
        "tissue_dissociation_protocol"
      )
      missing_cols <- setdiff(required_cols, names(meta_row))

      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "A",
          location = "metadata_access:col_check",
          message = "Checking required columns",
          data = list(
            gsm_id = gsm_id,
            missing_cols = missing_cols,
            required_cols = required_cols
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      if (length(missing_cols) > 0) {
        # #region agent log
        log_entry <- jsonlite::toJSON(
          list(
            sessionId = "debug-session",
            runId = "run1",
            hypothesisId = "A",
            location = "metadata_access:missing_cols",
            message = "Missing columns detected",
            data = list(gsm_id = gsm_id, missing_cols = missing_cols),
            timestamp = as.numeric(Sys.time()) * 1000
          ),
          auto_unbox = TRUE
        )
        cat(
          log_entry,
          "\n",
          file = "/home/h2048/.cursor/debug.log",
          append = TRUE
        )
        # #endregion
      }

      if ("sample" %in% names(meta_row)) {
        seurat_obj$sample <- meta_row$sample
      }
      if ("GEO" %in% names(meta_row)) {
        seurat_obj$GEO <- meta_row$GEO
      }
      if ("tissue" %in% names(meta_row)) {
        seurat_obj$tissue <- meta_row$tissue
      }
      if ("tissue_level_2" %in% names(meta_row)) {
        seurat_obj$tissue_level_2 <- meta_row$tissue_level_2
      }
      if ("disease_level_1" %in% names(meta_row)) {
        seurat_obj$disease_level_1 <- meta_row$disease_level_1
      }
      if ("disease_level_2" %in% names(meta_row)) {
        seurat_obj$disease_level_2 <- meta_row$disease_level_2
      }
      if ("condition" %in% names(meta_row)) {
        seurat_obj$condition <- meta_row$condition
      }
      if ("tissue_sampling_method" %in% names(meta_row)) {
        seurat_obj$tissue_sampling_method <- meta_row$tissue_sampling_method
      }
      if ("filtered_or_raw" %in% names(meta_row)) {
        seurat_obj$filtered_or_raw <- meta_row$filtered_or_raw
      }
      if ("target_cell" %in% names(meta_row)) {
        seurat_obj$target_cell <- meta_row$target_cell
      }
      if ("frozen_or_fresh" %in% names(meta_row)) {
        seurat_obj$frozen_or_fresh <- meta_row$frozen_or_fresh
      }
      if ("sex" %in% names(meta_row)) {
        seurat_obj$sex <- meta_row$sex
      }
      if ("age" %in% names(meta_row)) {
        seurat_obj$age <- meta_row$age
      }
      if ("platform" %in% names(meta_row)) {
        seurat_obj$platform <- meta_row$platform
      }
      if ("reference_genome" %in% names(meta_row)) {
        seurat_obj$reference_genome <- meta_row$reference_genome
      }
      if ("assay" %in% names(meta_row)) {
        seurat_obj$assay <- meta_row$assay
      }
      if ("tissue_dissociation_protocol" %in% names(meta_row)) {
        seurat_obj$tissue_dissociation_protocol <- meta_row$tissue_dissociation_protocol
      }

      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "A",
          location = "metadata_access:after_assign",
          message = "Metadata assignment completed",
          data = list(
            gsm_id = gsm_id,
            assigned_cols = intersect(required_cols, names(meta_row))
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      cat("  ✓ Metadata added\n")
    } else {
      # #region agent log
      log_entry <- jsonlite::toJSON(
        list(
          sessionId = "debug-session",
          runId = "run1",
          hypothesisId = "E",
          location = "metadata_access:not_found",
          message = "GSM ID not found in metadata",
          data = list(
            gsm_id = gsm_id,
            has_gsm_id_col = "gsm_id" %in% names(metadata),
            gsm_id_in_metadata = if ("gsm_id" %in% names(metadata)) {
              gsm_id %in% metadata$gsm_id
            } else {
              FALSE
            }
          ),
          timestamp = as.numeric(Sys.time()) * 1000
        ),
        auto_unbox = TRUE
      )
      cat(
        log_entry,
        "\n",
        file = "/home/h2048/.cursor/debug.log",
        append = TRUE
      )
      # #endregion

      cat("  ⚠️  No metadata found for", gsm_id, "\n")
    }

    cat("  Calculating QC metrics...\n")
    seurat_obj <- calculate_qc_metrics(seurat_obj)

    qc_stats <- seurat_obj@meta.data
    median_nfeature <- median(qc_stats$nFeature_RNA, na.rm = TRUE)
    median_ncount <- median(qc_stats$nCount_RNA, na.rm = TRUE)
    median_mt <- median(qc_stats$percent.mt, na.rm = TRUE)
    median_rb <- median(qc_stats$percent.rb, na.rm = TRUE)

    cat("  QC Summary:\n")
    cat(
      "    Median nFeature:",
      format(round(median_nfeature), big.mark = ","),
      "\n"
    )
    cat(
      "    Median nCount  :",
      format(round(median_ncount), big.mark = ","),
      "\n"
    )
    cat("    Median MT%     :", round(median_mt, 2), "\n")
    cat("    Median RB%     :", round(median_rb, 2), "\n")

    qc_plot_path <- file.path(QC_PLOTS_DIR, paste0(gsm_id, "_qc.png"))
    generate_qc_plot(seurat_obj, gsm_id, qc_plot_path)

    output_file <- file.path(OUTPUT_DIR, paste0(gsm_id, "_seurat.rds"))
    saveRDS(seurat_obj, output_file, compress = "gzip")
    cat("  ✓ Seurat object saved:", basename(output_file), "\n")

    bio_sample <- if (
      "gsm_id" %in% names(metadata) && gsm_id %in% metadata$gsm_id
    ) {
      meta_row_tmp <- metadata[.(gsm_id)]
      if ("sample" %in% names(meta_row_tmp)) {
        meta_row_tmp$sample
      } else {
        NA_character_
      }
    } else {
      NA_character_
    }

    processing_summary_list[[summary_counter]] <- data.table(
      Dataset = dataset_id,
      GSM_ID = gsm_id,
      Sample = bio_sample,
      Format = format_info$format,
      Cells_Raw = cells_raw,
      Cells_Filtered = cells_filtered,
      Cells_Removed = cells_removed,
      EmptyDrops_Method = ed_method,
      Genes_Original = genes_original,
      Genes_Final = genes_final,
      Genes_Unmapped_Removed = genes_unmapped_removed,
      N_Dup_Symbols = n_dup_symbols,
      Dedup_Method = SYMBOL_DEDUP_METHOD,
      Gene_Standardized = was_standardized,
      Median_nFeature = median_nfeature,
      Median_nCount = median_ncount,
      Median_MT = median_mt,
      Median_RB = median_rb,
      Status = "success",
      Reason = ifelse(is.null(ed_reason), "", ed_reason)
    )
    summary_counter <- summary_counter + 1

    total_processed <- total_processed + 1
    total_cells_raw <- total_cells_raw + cells_raw
    total_cells_filtered <- total_cells_filtered + cells_filtered

    rm(counts, seurat_obj, qc_stats)
    gc(FALSE)

    cat("\n")
  }
}

processing_summary <- rbindlist(processing_summary_list, fill = TRUE)

# ===== Final Summary =====
cat("============================================================\n")
cat("Processing Complete\n")
cat("============================================================\n\n")

cat("Summary Statistics:\n")
cat("  Total samples processed:", total_processed, "\n")
cat("  Total samples failed   :", total_failed, "\n")
cat("  Total cells (raw)      :", format(total_cells_raw, big.mark = ","), "\n")
cat(
  "  Total cells (filtered) :",
  format(total_cells_filtered, big.mark = ","),
  "\n"
)

if (RUN_EMPTYDROPS) {
  cells_removed_total <- total_cells_raw - total_cells_filtered
  if (total_cells_raw > 0) {
    cat(
      "  Cells removed by EmptyDrops:",
      format(cells_removed_total, big.mark = ","),
      " (",
      round(cells_removed_total / total_cells_raw * 100, 1),
      "%)\n",
      sep = ""
    )
  }
}

if (STANDARDIZE_GENE_IDS) {
  n_standardized <- sum(
    processing_summary$Gene_Standardized[
      processing_summary$Status == "success"
    ],
    na.rm = TRUE
  )
  cat(
    "  Samples with gene standardization:",
    n_standardized,
    "/",
    total_processed,
    "\n"
  )

  if (!KEEP_UNMAPPED_GENES) {
    total_unmapped_removed <- sum(
      processing_summary$Genes_Unmapped_Removed[
        processing_summary$Status == "success"
      ],
      na.rm = TRUE
    )
    cat(
      "  Total unmapped genes removed:",
      format(total_unmapped_removed, big.mark = ","),
      "\n"
    )
  }

  total_dup_symbols <- sum(
    processing_summary$N_Dup_Symbols[processing_summary$Status == "success"],
    na.rm = TRUE
  )
  if (total_dup_symbols > 0) {
    cat(
      "  Total duplicate symbols deduplicated:",
      format(total_dup_symbols, big.mark = ","),
      "\n"
    )
  }
}

cat("\n")

cat("Format Distribution:\n")
print(table(processing_summary$Format))
cat("\n")

cat("Status Distribution:\n")
print(table(processing_summary$Status))
cat("\n")

if (RUN_EMPTYDROPS) {
  cat("EmptyDrops Method Distribution:\n")
  print(table(processing_summary$EmptyDrops_Method))
  cat("\n")
}

summary_file <- file.path(OUTPUT_DIR, "processing_summary.csv")
fwrite(processing_summary, summary_file)
cat("Processing summary saved:", summary_file, "\n")
cat("Bad samples list saved  :", BAD_SAMPLES_TXT, "\n\n")

cat("Per-Dataset Statistics:\n")
dataset_summary <- processing_summary[
  Status == "success",
  .(
    Samples = .N,
    Total_Cells_Raw = as.integer(sum(Cells_Raw, na.rm = TRUE)),
    Total_Cells_Filtered = as.integer(sum(Cells_Filtered, na.rm = TRUE)),
    Pct_Removed = as.numeric(
      ifelse(
        sum(Cells_Raw, na.rm = TRUE) > 0,
        round(
          (sum(Cells_Raw, na.rm = TRUE) - sum(Cells_Filtered, na.rm = TRUE)) /
            sum(Cells_Raw, na.rm = TRUE) *
            100,
          1
        ),
        0.0 # ← 关键：使用 0.0 确保类型一致
      )
    ),
    Mean_Cells = as.numeric(round(mean(Cells_Filtered, na.rm = TRUE))),
    Median_Genes = round(median(Genes_Final, na.rm = TRUE)),
    N_Standardized = sum(Gene_Standardized),
    Total_Dup_Symbols = sum(N_Dup_Symbols),
    Genes_Unmapped_Removed = sum(Genes_Unmapped_Removed),
    Median_nFeature = round(median(Median_nFeature, na.rm = TRUE)),
    Median_nCount = round(median(Median_nCount, na.rm = TRUE)),
    Median_MT = round(median(Median_MT, na.rm = TRUE), 2),
    Median_RB = round(median(Median_RB, na.rm = TRUE), 2)
  ),
  by = Dataset
]
print(dataset_summary)

dataset_summary_file <- file.path(OUTPUT_DIR, "dataset_summary.csv")
fwrite(dataset_summary, dataset_summary_file)
cat("\nDataset summary saved:", dataset_summary_file, "\n\n")

cat("============================================================\n")
cat("✨ Pipeline Complete!\n")
cat("============================================================\n")
cat("Output locations:\n")
cat("  Seurat objects:", OUTPUT_DIR, "\n")
cat("  QC plots      :", QC_PLOTS_DIR, "\n")
cat("  Summary files :", OUTPUT_DIR, "\n")
cat("  Gene mappings :", OUTPUT_DIR, "/*_gene_mapping.csv.gz\n")
cat("============================================================\n")

rds_count <- length(list.files(OUTPUT_DIR, pattern = "_seurat\\.rds$"))
qc_count <- length(list.files(QC_PLOTS_DIR, pattern = "_qc\\.png$"))
mapping_count <- length(list.files(
  OUTPUT_DIR,
  pattern = "_gene_mapping\\.csv\\.gz$"
))

cat("\nGenerated files:\n")
cat("  RDS files    :", rds_count, "\n")
cat("  QC plots     :", qc_count, "\n")
cat("  Gene mappings:", mapping_count, "\n")

if (rds_count != qc_count || rds_count != total_processed) {
  cat("\n⚠️  Warning: File count mismatch detected\n")
  cat("  Expected:", total_processed, "files per type\n")
  cat("  Check processing_summary.csv for details\n")
}

cat("\n============================================================\n")
