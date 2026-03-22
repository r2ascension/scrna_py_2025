#!/usr/bin/env Rscript

# ============================================================
# gsm_one_by_one_debug.R
# Purpose:
#   Debug upstream input -> counts -> EmptyDrops -> Seurat QC
#   for ONE GSM/HRR sample at a time.
#
# Usage examples:
#   Rscript gsm_one_by_one_debug.R --data_root /home/h2048/data/source/1210 --sample_id GSM8499590 --outdir /home/h2048/data/R/debug/GSM8499590
#   Rscript gsm_one_by_one_debug.R --data_root /home/h2048/data/source/1210 --dataset_id GSE276503 --sample_id GSM8499590 --run_emptydrops TRUE
#   Rscript gsm_one_by_one_debug.R --rds /home/h2048/data/R/1210/Polyp/GSM8499590_seurat.rds --outdir /home/h2048/data/R/debug/GSM8499590_rds
#
# Notes:
# - If NA exists in counts@x, Seurat nCount/nFeature may become NA.
# - For dense matrix (csv/tsv) inputs, this script can optionally cap columns.
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(DropletUtils)
})

# ---------------------------
# Simple CLI parser
# ---------------------------
parse_kv_args <- function(args) {
  kv <- list()
  if (length(args) == 0) {
    return(kv)
  }
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) {
      stop("Bad arg: ", a)
    }
    key <- sub("^--", "", a)
    if (i == length(args)) {
      kv[[key]] <- TRUE
      break
    }
    val <- args[i + 1]
    if (startsWith(val, "--")) {
      kv[[key]] <- TRUE
      i <- i + 1
    } else {
      kv[[key]] <- val
      i <- i + 2
    }
  }
  kv
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  if (is.logical(x)) {
    return(x)
  }
  x <- tolower(as.character(x))
  x %in% c("true", "t", "1", "yes", "y")
}

# ---------------------------
# Utilities
# ---------------------------
cat_line <- function(...) cat(..., "\n", sep = "")

ensure_dir <- function(p) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  normalizePath(p, mustWork = FALSE)
}

write_text <- function(path, lines) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con)
}

# ---------------------------
# Locate sample dir
# ---------------------------
find_sample_dir <- function(data_root, sample_id, dataset_id = NULL) {
  if (!is.null(dataset_id)) {
    p <- file.path(data_root, dataset_id, sample_id)
    if (dir.exists(p)) {
      return(p)
    }
    return(NULL)
  }
  # scan one level: DATA_ROOT/*/SAMPLE_ID
  dsets <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)
  cand <- file.path(dsets, sample_id)
  cand <- cand[dir.exists(cand)]
  if (length(cand) == 0) {
    return(NULL)
  }
  cand[1]
}

# ---------------------------
# Fast NA inspection for dgCMatrix @x
# ---------------------------
cells_with_na_in_x <- function(m) {
  stopifnot(inherits(m, "dgCMatrix"))
  if (!anyNA(m@x)) {
    return(character())
  }
  na_pos <- which(is.na(m@x))
  # map x positions -> column index via p (column pointers)
  # findInterval returns col-1 (0-based); add 1
  col_idx <- findInterval(na_pos - 1, m@p) # 0..ncol-1
  unique(colnames(m)[col_idx])
}

genes_with_na_in_cells <- function(m, cells) {
  stopifnot(inherits(m, "dgCMatrix"))
  idx <- match(cells, colnames(m))
  idx <- idx[!is.na(idx)]
  res <- list()
  if (length(idx) == 0) {
    return(res)
  }

  for (j in idx) {
    s <- m@p[j] + 1
    e <- m@p[j + 1]
    if (s <= e) {
      xj <- m@x[s:e]
      ij <- m@i[s:e] + 1
      w <- which(is.na(xj))
      if (length(w)) {
        res[[colnames(m)[j]]] <- rownames(m)[ij[w]]
      }
    }
  }
  res
}

inspect_counts <- function(m, tag = "counts", max_show = 10) {
  stopifnot(inherits(m, "Matrix"))
  if (!inherits(m, "dgCMatrix")) {
    m <- as(m, "dgCMatrix")
  }

  cat_line(
    "    [",
    tag,
    "] class=",
    paste(class(m), collapse = "/"),
    " dim=",
    nrow(m),
    "x",
    ncol(m),
    " nnz=",
    length(m@x),
    " anyNA(x)=",
    anyNA(m@x)
  )

  if (length(m@x) > 0) {
    xr <- m@x
    cat_line(
      "    [",
      tag,
      "] x stats: min=",
      suppressWarnings(min(xr, na.rm = TRUE)),
      " median=",
      suppressWarnings(median(xr, na.rm = TRUE)),
      " max=",
      suppressWarnings(max(xr, na.rm = TRUE))
    )
    # integer-likeness
    frac <- xr[!is.na(xr)] - round(xr[!is.na(xr)])
    frac <- frac[is.finite(frac)]
    if (length(frac) > 0) {
      cat_line(
        "    [",
        tag,
        "] integer_like(~0 frac) rate=",
        round(mean(abs(frac) < 1e-8) * 100, 2),
        "%"
      )
    }
  }

  bad_cells <- cells_with_na_in_x(m)
  if (length(bad_cells) > 0) {
    cat_line(
      "    [",
      tag,
      "] cells_with_NA=",
      length(bad_cells),
      " examples=",
      paste(head(bad_cells, max_show), collapse = ", ")
    )
  }

  # sanity on library sizes
  lib <- Matrix::colSums(m, na.rm = FALSE)
  if (anyNA(lib)) {
    cat_line(
      "    [",
      tag,
      "] WARNING: colSums has NA (n=",
      sum(is.na(lib)),
      ")"
    )
  } else {
    cat_line(
      "    [",
      tag,
      "] libsize: min=",
      min(lib),
      " median=",
      median(lib),
      " max=",
      max(lib)
    )
  }
}

# ---------------------------
# Detect format
# ---------------------------
detect_data_format <- function(sample_dir, sample_id) {
  files <- list.files(sample_dir, full.names = FALSE)

  # H5
  h5 <- grep("\\.(h5|hdf5)$", files, value = TRUE, ignore.case = TRUE)
  if (length(h5) > 0) {
    return(list(format = "h5", path = file.path(sample_dir, h5[1])))
  }

  # 10x: barcodes + (features|genes) + matrix.mtx (all may be .gz)
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
  has_mtx <- any(grepl("matrix\\.mtx(\\.gz)?$", files, ignore.case = TRUE))
  if (has_barcodes && has_features && has_mtx) {
    standard_barcodes <- any(files %chin% c("barcodes.tsv.gz", "barcodes.tsv"))
    return(list(
      format = "10x",
      path = sample_dir,
      has_prefix = !standard_barcodes,
      sample_prefix = if (!standard_barcodes) sample_id else NULL
    ))
  }

  # matrix wide file
  mat <- grep("\\.(csv|tsv|txt|gz)$", files, value = TRUE, ignore.case = TRUE)
  mat <- mat[
    !grepl("features|barcodes|genes|meta|anno", mat, ignore.case = TRUE)
  ]
  if (length(mat) > 0) {
    prio <- grep("expr|expression|count|matrix", mat, ignore.case = TRUE)
    chosen <- if (length(prio) > 0) mat[prio[1]] else mat[1]
    return(list(format = "matrix", path = file.path(sample_dir, chosen)))
  }

  list(format = "unknown", path = NULL)
}

# ---------------------------
# Read 10x manually (prefix naming)
# ---------------------------
pick1 <- function(files, pattern) {
  x <- grep(pattern, files, value = TRUE, ignore.case = TRUE)
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[1]
}

open_text_con <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
}

read_lines_safe <- function(path) {
  con <- open_text_con(path)
  on.exit(close(con), add = TRUE)
  x <- readLines(con)
  x <- x[nzchar(x)]
  x
}

read_table_tsv_safe <- function(path) {
  # features/genes can be 2-3 cols; we only need col1
  con <- open_text_con(path)
  on.exit(close(con), add = TRUE)
  read.table(
    con,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

read_mtx_safe <- function(path) {
  con <- open_text_con(path)
  on.exit(close(con), add = TRUE)
  Matrix::readMM(con)
}

read_10x_manual <- function(sample_dir, sample_id, sample_prefix) {
  files <- list.files(sample_dir, full.names = FALSE)

  barcode_file <- pick1(
    files,
    paste0("^", sample_prefix, ".*barcodes\\.tsv(\\.gz)?$")
  )
  feature_file <- pick1(
    files,
    paste0("^", sample_prefix, ".*(features|genes)\\.tsv(\\.gz)?$")
  )
  matrix_file <- pick1(
    files,
    paste0("^", sample_prefix, ".*matrix\\.mtx(\\.gz)?$")
  )

  if (is.na(barcode_file) || is.na(feature_file) || is.na(matrix_file)) {
    return(NULL)
  }

  barcode_path <- file.path(sample_dir, barcode_file)
  feature_path <- file.path(sample_dir, feature_file)
  matrix_path <- file.path(sample_dir, matrix_file)

  cat_line("    10x(manual) files:")
  cat_line("      barcodes: ", barcode_file)
  cat_line("      features: ", feature_file)
  cat_line("      matrix  : ", matrix_file)

  barcodes <- read_lines_safe(barcode_path)
  feats <- read_table_tsv_safe(feature_path)
  mtx <- read_mtx_safe(matrix_path)

  if (!inherits(mtx, "dgCMatrix")) {
    mtx <- as(mtx, "dgCMatrix")
  }

  if (ncol(mtx) != length(barcodes)) {
    stop(
      "Barcode length mismatch: ncol(mtx)=",
      ncol(mtx),
      " vs length(barcodes)=",
      length(barcodes)
    )
  }
  if (nrow(mtx) != nrow(feats)) {
    stop(
      "Feature length mismatch: nrow(mtx)=",
      nrow(mtx),
      " vs nrow(features)=",
      nrow(feats)
    )
  }

  colnames(mtx) <- barcodes
  rownames(mtx) <- feats[[1]]

  mtx
}

# ---------------------------
# Load counts (h5 / 10x / matrix)
# ---------------------------
load_counts_data <- function(format_info, sample_id, matrix_max_cols = 5000) {
  fmt <- format_info$format
  path <- format_info$path

  tryCatch(
    {
      if (fmt == "h5") {
        cat_line("    Loading H5: ", basename(path))
        counts_list <- Read10X_h5(path, use.names = TRUE)
        counts <- if (
          is.list(counts_list) && !inherits(counts_list, "Matrix")
        ) {
          if ("Gene Expression" %in% names(counts_list)) {
            counts_list[["Gene Expression"]]
          } else {
            counts_list[[1]]
          }
        } else {
          counts_list
        }
      } else if (fmt == "10x") {
        cat_line("    Loading 10x directory: ", path)
        if (!isTRUE(format_info$has_prefix)) {
          counts_list <- Read10X(path, gene.column = 1)
          counts <- if (
            is.list(counts_list) && !inherits(counts_list, "Matrix")
          ) {
            if ("Gene Expression" %in% names(counts_list)) {
              counts_list[["Gene Expression"]]
            } else {
              counts_list[[1]]
            }
          } else {
            counts_list
          }
        } else {
          counts <- read_10x_manual(path, sample_id, format_info$sample_prefix)
          if (is.null(counts)) {
            stop(
              "Manual 10x read failed (prefix=",
              format_info$sample_prefix,
              ")"
            )
          }
        }
      } else if (fmt == "matrix") {
        cat_line("    Loading MATRIX file: ", basename(path))

        # Peek first line for delim + column count
        con <- open_text_con(path)
        on.exit(close(con), add = TRUE)
        first_line <- readLines(con, n = 1)
        close(con)

        delim <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
        tokens <- strsplit(first_line, split = delim, fixed = TRUE)[[1]]
        ncols_guess <- length(tokens)

        # header heuristic
        use_header <- FALSE
        if (ncols_guess > 1) {
          suppressWarnings({
            num_flags <- !is.na(as.numeric(tokens[-1]))
          })
          use_header <- !(all(num_flags))
        }
        cat_line(
          "    Detected delim=",
          ifelse(delim == "\t", "TAB", "COMMA"),
          " header=",
          use_header,
          " cols_first_line=",
          ncols_guess
        )

        # optional: cap columns for very wide matrices (debug safety)
        select_idx <- NULL
        if (ncols_guess > (matrix_max_cols + 1)) {
          select_idx <- seq_len(matrix_max_cols + 1) # gene + first N cells
          cat_line(
            "    NOTE: matrix seems very wide; reading only first ",
            matrix_max_cols,
            " cell columns for debug (plus gene column)."
          )
        }

        mat <- fread(
          path,
          sep = delim,
          header = use_header,
          data.table = FALSE,
          select = select_idx,
          showProgress = FALSE
        )

        # If header mis-detected -> re-read headerless
        if (use_header && ncol(mat) > 1) {
          cn <- colnames(mat)
          suppressWarnings({
            cn_numeric <- !is.na(as.numeric(cn[-1]))
          })
          if (all(cn_numeric)) {
            cat_line(
              "    Header likely mis-detected (numeric colnames). Re-reading without header."
            )
            mat <- fread(
              path,
              sep = delim,
              header = FALSE,
              data.table = FALSE,
              select = select_idx,
              showProgress = FALSE
            )
            use_header <- FALSE
          }
        }

        genes <- as.character(mat[[1]])
        mat <- mat[, -1, drop = FALSE]
        if (ncol(mat) == 0) {
          stop("No count columns detected after removing gene column.")
        }

        # safer numeric coercion
        mat_numeric <- data.matrix(mat) # may introduce NA if non-numeric strings exist
        rownames(mat_numeric) <- genes

        if (!use_header) {
          colnames(mat_numeric) <- paste0(
            sample_id,
            "_Cell",
            seq_len(ncol(mat_numeric))
          )
        }

        na_n <- sum(is.na(mat_numeric))
        if (na_n > 0) {
          cat_line(
            "    WARNING: NA introduced in dense numeric matrix: nNA=",
            na_n
          )
          # report first few offending columns
          na_by_col <- colSums(is.na(mat_numeric))
          bad <- which(na_by_col > 0)
          cat_line(
            "    WARNING: cols with NA=",
            length(bad),
            " examples=",
            paste(head(colnames(mat_numeric)[bad], 10), collapse = ", ")
          )
        }

        counts <- Matrix(mat_numeric, sparse = TRUE)
      } else {
        stop("Unknown format: ", fmt)
      }

      if (!inherits(counts, "dgCMatrix")) {
        counts <- as(counts, "dgCMatrix")
      }
      counts
    },
    error = function(e) {
      cat_line("    ❌ load_counts_data error: ", conditionMessage(e))
      NULL
    }
  )
}

# ---------------------------
# EmptyDrops wrapper (NA-safe)
# ---------------------------
run_emptydrops_filter <- function(counts, lower = 100, fdr = 0.01) {
  stopifnot(inherits(counts, "dgCMatrix"))
  n0 <- ncol(counts)

  if (n0 < 100) {
    cat_line("    EmptyDrops skipped (too few cells: ", n0, ")")
    return(list(
      counts = counts,
      method = "skipped_fewcells",
      before = n0,
      after = n0
    ))
  }

  cat_line("    Running EmptyDrops (lower=", lower, ", FDR<=", fdr, ")")
  set.seed(42)
  e.out <- emptyDrops(counts, lower = lower)

  is.cell <- e.out$FDR <= fdr
  is.cell[is.na(is.cell)] <- FALSE

  high.count <- Matrix::colSums(counts) > lower
  high.count[is.na(high.count)] <- FALSE

  is.cell <- is.cell | (is.na(e.out$FDR) & high.count)
  is.cell[is.na(is.cell)] <- FALSE # critical: prevent NA propagation

  kept <- counts[, is.cell, drop = FALSE]
  n1 <- ncol(kept)

  cat_line(
    "    EmptyDrops: before=",
    n0,
    " after=",
    n1,
    " removed=",
    (n0 - n1),
    " (",
    round((n0 - n1) / n0 * 100, 2),
    "%)"
  )
  list(counts = kept, method = "emptydrops", before = n0, after = n1)
}

# ---------------------------
# Seurat + QC inspect
# ---------------------------
inspect_seurat <- function(obj, tag = "seurat") {
  cat_line(
    "    [",
    tag,
    "] class=",
    paste(class(obj), collapse = "/"),
    " features=",
    nrow(obj),
    " cells=",
    ncol(obj)
  )
  md <- obj@meta.data
  cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
  cols <- cols[cols %in% colnames(md)]

  if (length(cols) > 0) {
    na_cells <- rownames(md)[apply(md[, cols, drop = FALSE], 1, function(v) {
      any(is.na(v))
    })]
    cat_line(
      "    [",
      tag,
      "] QC cols=",
      paste(cols, collapse = ", "),
      " cells_with_any_QC_NA=",
      length(na_cells)
    )
    if (length(na_cells) > 0) {
      cat_line(
        "    [",
        tag,
        "] example NA cells: ",
        paste(head(na_cells, 10), collapse = ", ")
      )
      print(md[head(na_cells, 5), cols, drop = FALSE])
    }
  } else {
    cat_line("    [", tag, "] QC cols not found.")
  }
}

# ---------------------------
# MAIN
# ---------------------------
args <- parse_kv_args(commandArgs(trailingOnly = TRUE))

DATA_ROOT <- args$data_root %||% "/home/h2048/data/source/1210"
SAMPLE_ID <- args$sample_id
DATASET_ID <- args$dataset_id
RDS_PATH <- args$rds
OUTDIR <- ensure_dir(args$outdir %||% file.path(getwd(), "gsm_debug_out"))
RUN_EMPTYDROPS <- as_bool(args$run_emptydrops, default = TRUE)
EMPTYDROPS_LOWER <- as.integer(args$emptydrops_lower %||% "100")
EMPTYDROPS_FDR <- as.numeric(args$emptydrops_fdr %||% "0.01")
MATRIX_MAX_COLS <- as.integer(args$matrix_max_cols %||% "5000")
SAVE_SEURAT <- as_bool(args$save_seurat, default = TRUE)

cat_line("============================================================")
cat_line("GSM/HRR One-by-One Debug")
cat_line("  OUTDIR: ", OUTDIR)
cat_line("  R version: ", R.version.string)
cat_line("  Seurat: ", as.character(packageVersion("Seurat")))
cat_line("  Matrix: ", as.character(packageVersion("Matrix")))
cat_line("  DropletUtils: ", as.character(packageVersion("DropletUtils")))
cat_line("============================================================")

# mode 1: inspect existing RDS
if (!is.null(RDS_PATH)) {
  if (!file.exists(RDS_PATH)) {
    stop("RDS not found: ", RDS_PATH)
  }
  cat_line("Mode: RDS inspect")
  cat_line("RDS: ", RDS_PATH)

  obj <- readRDS(RDS_PATH)
  if (!inherits(obj, "Seurat")) {
    stop("Not a Seurat object in RDS.")
  }

  # counts layer
  da <- DefaultAssay(obj)
  cat_line("Default assay: ", da)
  counts <- GetAssayData(obj, assay = da, layer = "counts")
  inspect_counts(counts, tag = "RDS:counts")

  bad_cells <- cells_with_na_in_x(as(counts, "dgCMatrix"))
  if (length(bad_cells) > 0) {
    na_genes <- genes_with_na_in_cells(as(counts, "dgCMatrix"), bad_cells)
    fwrite(
      data.table(cell = bad_cells),
      file.path(OUTDIR, "rds_cells_with_NA.csv")
    )
    # flatten gene list
    dt <- rbindlist(
      lapply(names(na_genes), function(cn) {
        data.table(cell = cn, gene = na_genes[[cn]])
      }),
      use.names = TRUE,
      fill = TRUE
    )
    fwrite(dt, file.path(OUTDIR, "rds_na_genes_by_cell.csv"))
    cat_line("Wrote: rds_cells_with_NA.csv, rds_na_genes_by_cell.csv")
  }

  inspect_seurat(obj, tag = "RDS:seurat")
  quit(save = "no", status = 0)
}

# mode 2: load from raw directory
if (is.null(SAMPLE_ID)) {
  stop("Provide --sample_id GSMxxxxxxx (or --rds /path/to.rds)")
}

sample_dir <- find_sample_dir(DATA_ROOT, SAMPLE_ID, DATASET_ID)
if (is.null(sample_dir)) {
  stop(
    "Sample directory not found under DATA_ROOT. sample_id=",
    SAMPLE_ID,
    if (!is.null(DATASET_ID)) paste0(" dataset_id=", DATASET_ID) else ""
  )
}

cat_line("Mode: RAW load")
cat_line("DATA_ROOT: ", DATA_ROOT)
if (!is.null(DATASET_ID)) {
  cat_line("DATASET_ID: ", DATASET_ID)
}
cat_line("SAMPLE_ID: ", SAMPLE_ID)
cat_line("SAMPLE_DIR: ", sample_dir)

format_info <- detect_data_format(sample_dir, SAMPLE_ID)
cat_line("Detected format: ", format_info$format)
if (format_info$format == "unknown") {
  stop("Unknown format in: ", sample_dir)
}

# load counts
counts <- load_counts_data(
  format_info,
  SAMPLE_ID,
  matrix_max_cols = MATRIX_MAX_COLS
)
if (is.null(counts)) {
  stop("Failed to load counts for ", SAMPLE_ID)
}

inspect_counts(counts, tag = paste0(SAMPLE_ID, ":after_load"))

# dump NA details if present
bad_cells0 <- cells_with_na_in_x(counts)
if (length(bad_cells0) > 0) {
  na_genes0 <- genes_with_na_in_cells(counts, bad_cells0)
  fwrite(
    data.table(cell = bad_cells0),
    file.path(OUTDIR, "after_load_cells_with_NA.csv")
  )
  dt0 <- rbindlist(
    lapply(names(na_genes0), function(cn) {
      data.table(cell = cn, gene = na_genes0[[cn]])
    }),
    use.names = TRUE,
    fill = TRUE
  )
  fwrite(dt0, file.path(OUTDIR, "after_load_na_genes_by_cell.csv"))
  cat_line(
    "Wrote: after_load_cells_with_NA.csv, after_load_na_genes_by_cell.csv"
  )
}

# emptydrops
if (RUN_EMPTYDROPS) {
  ed <- run_emptydrops_filter(
    counts,
    lower = EMPTYDROPS_LOWER,
    fdr = EMPTYDROPS_FDR
  )
  counts_ed <- ed$counts
  inspect_counts(counts_ed, tag = paste0(SAMPLE_ID, ":after_emptydrops"))

  bad_cells1 <- cells_with_na_in_x(counts_ed)
  if (length(bad_cells1) > 0) {
    na_genes1 <- genes_with_na_in_cells(counts_ed, bad_cells1)
    fwrite(
      data.table(cell = bad_cells1),
      file.path(OUTDIR, "after_emptydrops_cells_with_NA.csv")
    )
    dt1 <- rbindlist(
      lapply(names(na_genes1), function(cn) {
        data.table(cell = cn, gene = na_genes1[[cn]])
      }),
      use.names = TRUE,
      fill = TRUE
    )
    fwrite(dt1, file.path(OUTDIR, "after_emptydrops_na_genes_by_cell.csv"))
    cat_line(
      "Wrote: after_emptydrops_cells_with_NA.csv, after_emptydrops_na_genes_by_cell.csv"
    )
  }
} else {
  counts_ed <- counts
  cat_line("EmptyDrops: DISABLED")
}

# create seurat
cat_line("Creating Seurat object...")
obj <- CreateSeuratObject(
  counts = counts_ed,
  project = SAMPLE_ID,
  min.cells = 0,
  min.features = 0
)

# QC
cat_line("Computing QC metrics...")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

inspect_seurat(obj, tag = paste0(SAMPLE_ID, ":seurat"))
md <- obj@meta.data

# dump QC NA cells (if any)
qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
qc_cols <- qc_cols[qc_cols %in% colnames(md)]
if (length(qc_cols) > 0) {
  qc_na_cells <- rownames(md)[apply(
    md[, qc_cols, drop = FALSE],
    1,
    function(v) any(is.na(v))
  )]
  if (length(qc_na_cells) > 0) {
    fwrite(
      data.table(cell = qc_na_cells),
      file.path(OUTDIR, "seurat_qc_cells_with_NA.csv")
    )
    cat_line("Wrote: seurat_qc_cells_with_NA.csv")
  }
}

# save outputs
write_text(
  file.path(OUTDIR, "format_info.txt"),
  c(
    paste0("sample_id=", SAMPLE_ID),
    paste0("sample_dir=", sample_dir),
    paste0("format=", format_info$format),
    paste0("path=", format_info$path)
  )
)

if (SAVE_SEURAT) {
  out_rds <- file.path(OUTDIR, paste0(SAMPLE_ID, "_debug_seurat.rds"))
  saveRDS(obj, out_rds, compress = "gzip")
  cat_line("Saved Seurat RDS: ", out_rds)
}

cat_line("DONE. OUTDIR=", OUTDIR)
