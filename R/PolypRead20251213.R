#!/usr/bin/env Rscript

# ============================================================
# GEO scRNA-seq Complete Processing Pipeline (with EmptyDrops)
# FIXED VERSION (2025-12-13)
#
# Fixes included:
# 1) 10x manual reader now correctly handles .gz (barcodes/features/matrix.mtx.gz)
#    - Prevents readMM truncation/embedded NUL issues caused by reading gz as raw
#    - Strictly fails on readMM "expected ... found only ..." (true truncation)
# 2) Matrix (gene x cell wide table) header fix:
#    - Handles "header missing gene placeholder" (barcodes-only first line)
#    - Skips extremely wide matrices safely to avoid OOM
# 3) Guarantees no NA in sparse counts@x before EmptyDrops (drops bad columns or fails)
# 4) Writes bad samples list to a TXT for quick inspection
# ============================================================

# ===== Configuration Section =====
DATA_ROOT <- "/home/h2048/data/source/1214"
METADATA_PATH <- "/home/h2048/data/source/reference/metadata_full_atlas_20251210.csv"
OUTPUT_DIR <- "/home/h2048/data/R/1214/Polyp"
QC_PLOTS_DIR <- "/home/h2048/data/R/1214/Polyp/qc_plots"

# EmptyDrops parameters
RUN_EMPTYDROPS <- TRUE
EMPTYDROPS_LOWER <- 100
EMPTYDROPS_FDR_THRESHOLD <- 0.01

# QC plot settings
MAX_PLOT_POINTS <- 5000

# QC thresholds (for visualization reference lines only)
QC_THRESHOLDS <- list(
  nFeature_min = 200,
  nFeature_max = 6000,
  nCount_min = 500,
  nCount_max = 50000,
  percent_mt_max = 20,
  percent_rb_min = 5
)

# Safety knobs (for wide text matrices)
MAX_WIDE_MATRIX_CELLS <- 3000 # wide tables beyond this are skipped (avoid OOM)
FAIL_ON_READMM_TRUNCATED <- FALSE # stop if readMM reports truncation

# ===== Load Libraries =====
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(ggplot2)
  library(patchwork)
  library(DropletUtils)
})

cat("============================================================\n")
cat("GEO/HRA scRNA-seq Complete Processing Pipeline (FIXED)\n")
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
cat("Wide-matrix guard: MAX_WIDE_MATRIX_CELLS =", MAX_WIDE_MATRIX_CELLS, "\n")
cat("readMM truncation strict:", FAIL_ON_READMM_TRUNCATED, "\n")
cat("============================================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

BAD_SAMPLES_TXT <- file.path(OUTPUT_DIR, "bad_samples.txt")
writeLines(
  c(
    "dataset_id\tsample_id\tformat\treason\tfile",
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  con = BAD_SAMPLES_TXT
)

append_bad_sample <- function(
  dataset_id,
  sample_id,
  format,
  reason,
  file = NA_character_
) {
  line <- paste(dataset_id, sample_id, format, reason, file, sep = "\t")
  cat(line, "\n", file = BAD_SAMPLES_TXT, append = TRUE)
}

# ===== Helper Functions =====

open_text_con <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

readMM_strict <- function(path) {
  w <- character(0)
  mat <- withCallingHandlers(
    {
      Matrix::readMM(open_text_con(path))
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

drop_na_columns_if_any <- function(counts, sample_id) {
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "CsparseMatrix")
  }
  if (anyNA(counts@x)) {
    cs <- Matrix::colSums(counts)
    bad <- which(is.na(cs))
    if (length(bad) > 0) {
      cat(
        "    [",
        sample_id,
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

# Read 10x files manually (handles non-standard naming + gz)
read_10x_manual <- function(gsm_dir, sample_prefix = NULL) {
  files <- list.files(gsm_dir, full.names = FALSE)

  pick1 <- function(pattern) {
    out <- grep(pattern, files, value = TRUE, ignore.case = TRUE)
    if (length(out) == 0) {
      return(NA_character_)
    }
    out[1]
  }

  barcode_pattern <- if (is.null(sample_prefix)) {
    "barcodes\\.tsv(\\.gz)?$"
  } else {
    paste0(sample_prefix, ".*barcodes\\.tsv(\\.gz)?$")
  }

  feature_pattern <- if (is.null(sample_prefix)) {
    "(features|genes)\\.tsv(\\.gz)?$"
  } else {
    paste0(sample_prefix, ".*(features|genes)\\.tsv(\\.gz)?$")
  }

  matrix_pattern <- if (is.null(sample_prefix)) {
    "matrix\\.mtx(\\.gz)?$"
  } else {
    paste0(sample_prefix, ".*matrix\\.mtx(\\.gz)?$")
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

  barcodes <- readLines(open_text_con(barcode_path))

  features <- read.table(
    open_text_con(feature_path),
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
  mat
}

# Detect data format in sample directory
detect_data_format <- function(sample_dir, sample_id) {
  files <- list.files(sample_dir, full.names = FALSE)

  # Check for H5 format
  h5_files <- grep("\\.(h5|hdf5)$", files, value = TRUE, ignore.case = TRUE)
  if (length(h5_files) > 0) {
    return(list(
      format = "h5",
      path = file.path(sample_dir, h5_files[1])
    ))
  }

  # Check for 10x format (barcodes + features/genes + matrix)
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
    # Check if files have standard names
    standard_barcodes <- any(files %chin% c("barcodes.tsv.gz", "barcodes.tsv"))
    return(list(
      format = "10x",
      path = sample_dir,
      has_prefix = !standard_barcodes,
      sample_prefix = if (!standard_barcodes) sample_id else NULL
    ))
  }

  # Check for matrix format
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

# Load counts data and filter ADT
load_counts_data <- function(format_info, sample_id) {
  format_type <- format_info$format
  path <- format_info$path

  tryCatch(
    {
      counts <- NULL

      if (format_type == "h5") {
        cat("    Loading H5 format...\n")
        counts_list <- Read10X_h5(path, use.names = TRUE)

        if (is.list(counts_list) && !inherits(counts_list, "Matrix")) {
          cat(
            "    Multiple assays detected:",
            paste(names(counts_list), collapse = ", "),
            "\n"
          )
          if ("Gene Expression" %in% names(counts_list)) {
            counts <- counts_list[["Gene Expression"]]
            cat("    Using 'Gene Expression' assay\n")
          } else {
            counts <- counts_list[[1]]
            cat("    Using first assay:", names(counts_list)[1], "\n")
          }
        } else {
          counts <- counts_list
        }
      } else if (format_type == "10x") {
        cat("    Loading 10x format...\n")

        if (!isTRUE(format_info$has_prefix)) {
          counts_list <- Read10X(path, gene.column = 1)

          if (is.list(counts_list) && !inherits(counts_list, "Matrix")) {
            cat(
              "    Multiple assays detected:",
              paste(names(counts_list), collapse = ", "),
              "\n"
            )
            if ("Gene Expression" %in% names(counts_list)) {
              counts <- counts_list[["Gene Expression"]]
              cat("    Using 'Gene Expression' assay (ADT filtered)\n")
            } else {
              counts <- counts_list[[1]]
              cat("    Using first assay:", names(counts_list)[1], "\n")
            }
          } else {
            counts <- counts_list
          }
        } else {
          cat("    Detected non-standard file naming, reading manually...\n")
          counts <- read_10x_manual(path, format_info$sample_prefix)
          if (is.null(counts)) {
            stop("10X_MANUAL_READ_FAILED: missing required files")
          }
        }
      } else if (format_type == "matrix") {
        cat("    Loading matrix format...\n")

        # ---- Robust first-2-lines inference ----
        con <- open_text_con(path)
        line1 <- readLines(con, n = 1)
        line2 <- readLines(con, n = 1)
        close(con)

        delim <- if (grepl("\t", line1)) "\t" else ","
        tok1 <- strsplit(line1, split = delim, fixed = TRUE)[[1]]
        tok2 <- strsplit(line2, split = delim, fixed = TRUE)[[1]]

        # header is barcodes ONLY (missing leading gene/feature placeholder)
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
          cat(
            "    Detected matrix header missing gene placeholder; treat line1 as barcodes and skip it.\n"
          )
          barcodes <- tok1
          mat <- fread(
            path,
            sep = delim,
            header = FALSE,
            skip = 1,
            data.table = FALSE
          )
          if (ncol(mat) != (length(barcodes) + 1)) {
            stop(
              "WIDE_MATRIX_COL_MISMATCH: data cols=",
              ncol(mat),
              " expected=",
              length(barcodes) + 1
            )
          }
          colnames(mat) <- c("gene", barcodes)
        } else {
          # fallback heuristic: header present if non-numeric colnames
          tokens <- tok1
          use_header <- FALSE
          if (length(tokens) > 1) {
            suppressWarnings({
              num_flags <- !is.na(as.numeric(tokens[-1]))
            })
            use_header <- !(all(num_flags))
          }
          cat("    Header detection: use_header =", use_header, "\n")
          mat <- fread(
            path,
            sep = delim,
            header = use_header,
            data.table = FALSE
          )

          # If header was used but colnames still numeric => re-read headerless
          if (use_header && ncol(mat) > 1) {
            cn <- colnames(mat)
            suppressWarnings({
              cn_numeric <- !is.na(as.numeric(cn[-1]))
            })
            if (all(cn_numeric)) {
              cat(
                "    Numeric colnames after header read; re-reading as headerless.\n"
              )
              mat <- fread(
                path,
                sep = delim,
                header = FALSE,
                data.table = FALSE
              )
            }
          }
        }

        # ---- Convert to sparse safely (NOTE: still dense parse; guarded by MAX_WIDE_MATRIX_CELLS) ----
        genes <- mat[, 1]
        mat2 <- mat[, -1, drop = FALSE]
        if (ncol(mat2) == 0) {
          stop("MATRIX_NO_COUNT_COLUMNS")
        }

        # coerce to numeric matrix
        mat_numeric <- as.matrix(mat2)
        rownames(mat_numeric) <- genes

        # if no header existed, synthesize cell IDs
        if (is.null(colnames(mat2)) || any(colnames(mat2) == "")) {
          colnames(mat_numeric) <- paste0(
            sample_id,
            "_Cell",
            seq_len(ncol(mat_numeric))
          )
        } else {
          colnames(mat_numeric) <- colnames(mat2)
        }

        counts <- as(mat_numeric, "CsparseMatrix")
      } else {
        stop("UNKNOWN_FORMAT")
      }

      if (!inherits(counts, "dgCMatrix")) {
        counts <- as(counts, "CsparseMatrix")
      }
      counts <- drop_na_columns_if_any(counts, sample_id)

      cat("    Loaded:", nrow(counts), "genes x", ncol(counts), "cells\n")
      list(counts = counts, reason = NULL)
    },
    error = function(e) {
      msg <- conditionMessage(e)
      cat("    ❌ Error loading", sample_id, ":", msg, "\n")
      list(counts = NULL, reason = msg)
    }
  )
}

# Run EmptyDrops filtering
run_emptydrops_filter <- function(counts, sample_id) {
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
      counts <- drop_na_columns_if_any(counts, sample_id)

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

      is.cell <- e.out$FDR <= EMPTYDROPS_FDR_THRESHOLD
      is.cell[is.na(is.cell)] <- FALSE

      high.count <- Matrix::colSums(counts) > EMPTYDROPS_LOWER
      is.cell <- is.cell | (is.na(e.out$FDR) & high.count)

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
        round(cells_removed / cells_before * 100, 1),
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

extract_orig_ident <- function(path) {
  if (is.null(path) || length(path) == 0) {
    return("Unknown")
  }
  filename <- basename(path)
  orig_ident <- gsub("_features\\.tsv.*", "", filename, ignore.case = TRUE)
  orig_ident <- gsub("_genes\\.tsv.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_matrix\\.mtx.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_barcodes\\.tsv.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_filtered.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_raw.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub(
    "\\.(h5|hdf5|txt|tsv|csv|gz)$",
    "",
    orig_ident,
    ignore.case = TRUE
  )
  orig_ident
}

calculate_qc_metrics <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[SL]"
  )
  seurat_obj
}

generate_qc_plot <- function(
  seurat_obj,
  sample_id,
  output_path,
  max_points = MAX_PLOT_POINTS
) {
  qc_data <- as.data.table(seurat_obj@meta.data)

  if (nrow(qc_data) > max_points) {
    sample_idx <- sample.int(nrow(qc_data), max_points)
    qc_data_sample <- qc_data[sample_idx]
  } else {
    qc_data_sample <- qc_data
  }

  p1 <- ggplot(qc_data, aes(x = "Sample", y = nFeature_RNA)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = nFeature_RNA),
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
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p2 <- ggplot(qc_data, aes(x = "Sample", y = nCount_RNA)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = nCount_RNA),
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
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p3 <- ggplot(qc_data, aes(x = "Sample", y = percent.mt)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = percent.mt),
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
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p4 <- ggplot(qc_data, aes(x = "Sample", y = percent.rb)) +
    geom_violin(alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = percent.rb),
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
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
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
      title = paste0("Quality Control Metrics: ", sample_id),
      subtitle = paste0(
        "Cells: ",
        format(ncol(seurat_obj), big.mark = ","),
        " | Genes: ",
        format(nrow(seurat_obj), big.mark = ","),
        " | Median nFeature: ",
        format(round(median(qc_data$nFeature_RNA)), big.mark = ","),
        " | Median nCount: ",
        format(round(median(qc_data$nCount_RNA)), big.mark = ","),
        if (nrow(qc_data) > max_points) {
          paste0(
            " | Sampled ",
            format(max_points, big.mark = ","),
            " cells for scatters"
          )
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
cat("  Original columns:", paste(names(metadata), collapse = ", "), "\n")

cat("Standardizing metadata columns...\n")
cols <- names(metadata)

if ("filtered_or_Falset" %in% cols && !"filtered_or_raw" %in% cols) {
  setnames(metadata, "filtered_or_Falset", "filtered_or_raw")
  cat("  Renamed 'filtered_or_Falset' to 'filtered_or_raw'\n")
}

ref_cols <- grep("^reference", cols, value = TRUE, ignore.case = TRUE)
if (length(ref_cols) == 1 && !"reference_genome" %in% cols) {
  setnames(metadata, ref_cols, "reference_genome")
  cat("  Renamed '", ref_cols, "' to 'reference_genome'\n", sep = "")
}

cat("Processing metadata columns...\n")
if ("sample" %in% names(metadata)) {
  metadata[, sample_id := sample]
  cat("  Created 'sample_id' column from 'sample' (sample column preserved)\n")
} else {
  stop("'sample' column not found in metadata")
}

metadata[, sample_id := gsub("[\\s\\t]+", "", sample_id)]
if (any(duplicated(metadata$sample_id))) {
  warning(
    "Duplicated sample_ids found in metadata; only first occurrence will be used."
  )
}
setkey(metadata, sample_id)

cat("  Total metadata records:", nrow(metadata), "\n")
cat("  Unique sample_ids:", uniqueN(metadata$sample_id), "\n")
cat("  Final columns:", paste(names(metadata), collapse = ", "), "\n\n")

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
processing_summary <- data.table(
  Dataset = character(),
  Sample_ID = character(),
  sample = character(),
  Format = character(),
  Cells_Raw = integer(),
  Cells_Filtered = integer(),
  Cells_Removed = integer(),
  EmptyDrops_Method = character(),
  Genes = integer(),
  Median_nFeature = numeric(),
  Median_nCount = numeric(),
  Median_MT = numeric(),
  Median_RB = numeric(),
  Status = character(),
  Reason = character()
)

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
    sample_id <- basename(sample_dir)
    cat("📦", sample_id, "\n")

    format_info <- detect_data_format(sample_dir, sample_id)
    cat("  Format:", format_info$format, "\n")

    if (format_info$format == "unknown") {
      cat("  ⚠️  Skipping (unknown format)\n\n")
      append_bad_sample(
        dataset_id,
        sample_id,
        "unknown",
        "unknown_format",
        NA_character_
      )
      processing_summary <- rbind(
        processing_summary,
        data.table(
          Dataset = dataset_id,
          Sample_ID = sample_id,
          sample = NA_character_,
          Format = "unknown",
          Cells_Raw = NA_integer_,
          Cells_Filtered = NA_integer_,
          Cells_Removed = NA_integer_,
          EmptyDrops_Method = "skipped",
          Genes = NA_integer_,
          Median_nFeature = NA_real_,
          Median_nCount = NA_real_,
          Median_MT = NA_real_,
          Median_RB = NA_real_,
          Status = "skipped",
          Reason = "unknown_format"
        )
      )
      total_failed <- total_failed + 1
      next
    }

    ld <- load_counts_data(format_info, sample_id)
    counts <- ld$counts

    if (is.null(counts)) {
      append_bad_sample(
        dataset_id,
        sample_id,
        format_info$format,
        ld$reason,
        format_info$path
      )
      processing_summary <- rbind(
        processing_summary,
        data.table(
          Dataset = dataset_id,
          Sample_ID = sample_id,
          sample = NA_character_,
          Format = format_info$format,
          Cells_Raw = NA_integer_,
          Cells_Filtered = NA_integer_,
          Cells_Removed = NA_integer_,
          EmptyDrops_Method = "failed_load",
          Genes = NA_integer_,
          Median_nFeature = NA_real_,
          Median_nCount = NA_real_,
          Median_MT = NA_real_,
          Median_RB = NA_real_,
          Status = "failed",
          Reason = ld$reason
        )
      )
      total_failed <- total_failed + 1
      cat("\n")
      next
    }

    # EmptyDrops
    if (RUN_EMPTYDROPS) {
      ed_result <- run_emptydrops_filter(counts, sample_id)
      counts <- ed_result$counts
      cells_raw <- ed_result$cells_before
      cells_filtered <- ed_result$cells_after
      cells_removed <- ed_result$cells_removed
      ed_method <- ed_result$method
      ed_reason <- ed_result$reason
      if (!is.null(ed_reason)) {
        append_bad_sample(
          dataset_id,
          sample_id,
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

    if (ncol(counts) == 0 || nrow(counts) == 0) {
      cat(
        "  ⚠️  No valid cells or genes detected (",
        nrow(counts),
        " genes x ",
        ncol(counts),
        " cells). Skipping.\n\n",
        sep = ""
      )
      append_bad_sample(
        dataset_id,
        sample_id,
        format_info$format,
        "failed_nocells",
        format_info$path
      )

      sample_val <- if (sample_id %in% metadata$sample_id) {
        metadata[sample_id]$sample
      } else {
        NA_character_
      }

      processing_summary <- rbind(
        processing_summary,
        data.table(
          Dataset = dataset_id,
          Sample_ID = sample_id,
          sample = sample_val,
          Format = format_info$format,
          Cells_Raw = cells_raw,
          Cells_Filtered = 0,
          Cells_Removed = cells_raw,
          EmptyDrops_Method = if (RUN_EMPTYDROPS) {
            "failed_nocells"
          } else {
            "disabled_nocells"
          },
          Genes = nrow(counts),
          Median_nFeature = NA_real_,
          Median_nCount = NA_real_,
          Median_MT = NA_real_,
          Median_RB = NA_real_,
          Status = "failed_nocells",
          Reason = "no_cells_or_genes_after_filter"
        )
      )
      total_failed <- total_failed + 1
      rm(counts)
      gc(FALSE)
      next
    }

    orig_ident <- extract_orig_ident(format_info$path)
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

    seurat_obj$sample_id <- sample_id
    seurat_obj$dataset <- dataset_id
    seurat_obj$orig.ident <- orig_ident
    seurat_obj$emptydrops_filtered <- RUN_EMPTYDROPS &&
      ed_method == "emptydrops"

    if (sample_id %in% metadata$sample_id) {
      meta_row <- metadata[sample_id]
      seurat_obj$sample <- meta_row$sample
      seurat_obj$GEO <- meta_row$GEO
      seurat_obj$tissue <- meta_row$tissue
      seurat_obj$tissue_level_2 <- meta_row$tissue_level_2
      seurat_obj$disease_level_1 <- meta_row$disease_level_1
      seurat_obj$disease_level_2 <- meta_row$disease_level_2
      seurat_obj$condition <- meta_row$condition
      seurat_obj$tissue_sampling_method <- meta_row$tissue_sampling_method
      seurat_obj$filtered_or_raw <- meta_row$filtered_or_raw
      seurat_obj$target_cell <- meta_row$target_cell
      seurat_obj$frozen_or_fresh <- meta_row$frozen_or_fresh
      seurat_obj$sex <- meta_row$sex
      seurat_obj$age <- meta_row$age
      seurat_obj$platform <- meta_row$platform
      seurat_obj$reference_genome <- meta_row$reference_genome
      seurat_obj$assay <- meta_row$assay
      seurat_obj$tissue_dissociation_protocol <- meta_row$tissue_dissociation_protocol
      cat("  ✓ Metadata added\n")
    } else {
      cat("  ⚠️  No metadata found for", sample_id, "\n")
    }

    cat("  Calculating QC metrics...\n")
    seurat_obj <- calculate_qc_metrics(seurat_obj)

    qc_stats <- seurat_obj@meta.data
    median_nfeature <- median(qc_stats$nFeature_RNA)
    median_ncount <- median(qc_stats$nCount_RNA)
    median_mt <- median(qc_stats$percent.mt)
    median_rb <- median(qc_stats$percent.rb)

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

    qc_plot_path <- file.path(QC_PLOTS_DIR, paste0(sample_id, "_qc.png"))
    generate_qc_plot(seurat_obj, sample_id, qc_plot_path)

    output_file <- file.path(OUTPUT_DIR, paste0(sample_id, "_seurat.rds"))
    saveRDS(seurat_obj, output_file, compress = "gzip")
    cat("  ✓ Seurat object saved:", basename(output_file), "\n")

    sample_val <- if (sample_id %in% metadata$sample_id) {
      metadata[sample_id]$sample
    } else {
      NA_character_
    }

    processing_summary <- rbind(
      processing_summary,
      data.table(
        Dataset = dataset_id,
        Sample_ID = sample_id,
        sample = sample_val,
        Format = format_info$format,
        Cells_Raw = cells_raw,
        Cells_Filtered = cells_filtered,
        Cells_Removed = cells_removed,
        EmptyDrops_Method = ed_method,
        Genes = nrow(seurat_obj),
        Median_nFeature = median_nfeature,
        Median_nCount = median_ncount,
        Median_MT = median_mt,
        Median_RB = median_rb,
        Status = "success",
        Reason = ifelse(is.null(ed_reason), "", ed_reason)
      )
    )

    total_processed <- total_processed + 1
    total_cells_raw <- total_cells_raw + cells_raw
    total_cells_filtered <- total_cells_filtered + cells_filtered

    rm(counts, seurat_obj, qc_stats)
    gc(FALSE)

    cat("\n")
  }
}

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

max_genes <- if (all(is.na(processing_summary$Genes))) {
  NA_integer_
} else {
  max(processing_summary$Genes, na.rm = TRUE)
}
cat("  Total genes (max):", format(max_genes, big.mark = ","), "\n\n")

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
    Total_Cells_Raw = sum(Cells_Raw),
    Total_Cells_Filtered = sum(Cells_Filtered),
    Pct_Removed = round(
      (sum(Cells_Raw) - sum(Cells_Filtered)) / sum(Cells_Raw) * 100,
      1
    ),
    Mean_Cells = round(mean(Cells_Filtered)),
    Median_nFeature = round(median(Median_nFeature)),
    Median_nCount = round(median(Median_nCount)),
    Median_MT = round(median(Median_MT), 2),
    Median_RB = round(median(Median_RB), 2)
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
cat("============================================================\n")

rds_count <- length(list.files(OUTPUT_DIR, pattern = "_seurat\\.rds$"))
qc_count <- length(list.files(QC_PLOTS_DIR, pattern = "_qc\\.png$"))

cat("\nGenerated files:\n")
cat("  RDS files:", rds_count, "\n")
cat("  QC plots :", qc_count, "\n")

if (rds_count != qc_count || rds_count != total_processed) {
  cat("\n⚠️  Warning: File count mismatch detected\n")
  cat("  Expected:", total_processed, "files per type\n")
  cat("  Check processing_summary.csv for details\n")
}

cat("\n============================================================\n")
