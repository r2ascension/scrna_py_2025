#!/usr/bin/env Rscript
# =============================================================================
# RDS Folder → QC → DoubletFinder → DecontX → Smart Merge Pipeline (Seurat v5)
#
# VERSION: v4_MEMORY_OPTIMIZED (2024-12-16) ⭐ Memory-Efficient Incremental Merge
# NEW FEATURES vs v3:
#   - ⭐ Lazy loading: Only record paths, load files during merge
#   - ⭐ Incremental merge: Read-merge-free cycle to minimize memory
#   - ⭐ Cell IDs added once during save, not during merge
#   - ⭐ RNA_restored priority: Use RNA_restored as main assay if available
#   - Memory cleanup after each sample processing
#
# FEATURES from v3_STORE_RESTORED_SKIP_EXISTING:
#   - ⭐ Detects and skips already processed cleaned RDS files
#   - SKIP_EXISTING_CLEANED parameter to enable/disable skip behavior
#   - FORCE_REPROCESS parameter to force reprocessing all samples
#   - Improved logging for skipped samples
#
# FEATURES from v3_STORE_RESTORED:
#   - ⭐ Stores restored counts in new assay "RNA_restored" (optional)
#   - SAVE_RESTORED_COUNTS parameter to control storage behavior
#   - Restored counts available for inspection and downstream analysis
#
# FEATURES from v2_DENORM:
#   - De-lognormalization functions to restore raw counts from log-normalized data
#   - Auto-detection and restoration of raw counts for DecontX compatibility
#   - restore_raw_counts_log(): Core de-lognormalization function
#   - restore_counts_from_seurat(): Seurat-specific wrapper
#
# Input: Folder with multiple .rds files (each = 1 sample, EmptyDrops already done)
# Output:
#   - cleaned_samples/*.rds (per-sample cleaned objects with cell IDs)
#   - qc_stats/per_sample_qc_stats.csv
#   - merged/merged_seurat.rds
#   - merged/merged_seurat_final.rds (singlets + decontX filtered)
#
# Author: r2end
# Date: 2024-12-16
# Seurat: v5 compatible
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DoubletFinder)
  library(celda) # DecontX
  library(data.table)
  library(ggplot2)
  library(jsonlite) # For debug logging
})

options(stringsAsFactors = FALSE)
set.seed(42)

# =============================================================================
# CLI Arguments
# Usage: Rscript script.R /path/to/input_rds_dir /path/to/output_dir
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
INPUT_DIR <- if (length(args) >= 1) {
  args[1]
} else {
  stop("Please provide INPUT_DIR")
}
OUTPUT_DIR <- if (length(args) >= 2) {
  args[2]
} else {
  paste0(INPUT_DIR, "_QC_output")
}

# =============================================================================
# PARAMETERS (adjust as needed)
# =============================================================================
ASSAY_PREFERRED <- "RNA"

# QC thresholds (based on project standards)
NFEATURE_MIN <- 200
NFEATURE_MAX <- 6000
MT_PATTERN <- "^MT-"
MT_MAX <- 20
RB_PATTERN <- "^RP[SL]"
RB_MAX <- 40

# Regression outlier detection
DO_REGRESSION_OUTLIER <- TRUE
OUTLIER_PRED_LEVEL <- 0.999

# DoubletFinder
DO_DOUBLET_FINDER <- TRUE
DOUBLET_RATE <- 0.06
DOUBLET_PCS <- 30

# DecontX
DO_DECONTX <- TRUE
DECONTX_MAX_ITER <- 500
FILTER_BY_DECONTX <- TRUE
DECONTX_CONTAM_MAX <- 0.25

# De-lognormalization (v2 features)
TRY_DENORMALIZE <- TRUE  # Attempt to restore raw counts if data is log-normalized
DENORM_SCALE_FACTOR <- 10000  # Scale factor used in LogNormalize

# ⭐ NEW v3: Storage control for restored counts
SAVE_RESTORED_COUNTS <- TRUE  # Store restored counts in new assay "RNA_restored"
RESTORED_ASSAY_NAME <- "RNA_restored"  # Name for restored counts assay

# Gene filtering (final merge)
MIN_CELLS_PER_GENE <- 3

# Output
PLOT_QC_PDF <- TRUE
SAVE_INTERMEDIATE <- TRUE

# ⭐ NEW: Skip already processed samples
SKIP_EXISTING_CLEANED <- TRUE  # Skip processing if cleaned RDS already exists
FORCE_REPROCESS <- FALSE  # Set to TRUE to force reprocessing all samples

# =============================================================================
# Helper Functions
# =============================================================================

# =============================================================================
# De-lognormalization Functions (from v2)
# =============================================================================

#' Restore raw counts from log-normalized data
#'
#' This function reverses the log(counts + 1) normalization to recover raw counts.
#' Formula: counts = exp(normalized_data) - 1, then multiply by size factors
#'
#' @param normalized_data Matrix of log-normalized values
#' @param size_factors Vector of size factors (normalization factors) for each cell
#'        Length must match ncol(normalized_data)
#' @return Matrix of restored raw counts
#'
restore_raw_counts_log <- function(normalized_data, size_factors) {
  # Validate input
  if (length(size_factors) != ncol(normalized_data)) {
    stop("Size factors length does not match the number of samples in the normalized data.")
  }

  # Reverse log transformation using expm1 for numerical stability
  raw_counts <- expm1(normalized_data)

  # Multiply by size factors to restore original scale
  raw_counts <- sweep(raw_counts, 2, size_factors, "*")

  return(raw_counts)
}

#' Restore raw counts from Seurat object with LogNormalize
#'
#' Extracts log-normalized data from a Seurat object and restores it to raw counts
#' using the stored nCount (total UMI per cell) as size factors
#'
#' @param seurat_obj Seurat object with log-normalized data
#' @param assay Assay name (default: "RNA")
#' @param scale_factor Scale factor used in NormalizeData (default: 10000)
#' @return Matrix of restored raw counts (dgCMatrix)
#'
restore_counts_from_seurat <- function(seurat_obj, assay = "RNA", scale_factor = 10000) {
  # Extract normalized data (log-transformed)
  normalized_data <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "data"),
    error = function(e1) {
      tryCatch(
        GetAssayData(seurat_obj, assay = assay, slot = "data"),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(normalized_data)) {
    stop("Cannot extract normalized data from Seurat object")
  }

  # Get QC column name for nCount
  qc_cols <- get_qc_cols(assay)
  nCount_col <- qc_cols$nCount

  # Extract nCount (total UMI per cell) as size factors
  if (!nCount_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("QC column '%s' not found in metadata", nCount_col))
  }

  nCount <- seurat_obj@meta.data[[nCount_col]]

  # Calculate size factors
  size_factors <- nCount / scale_factor

  # Restore raw counts
  raw_counts <- restore_raw_counts_log(normalized_data, size_factors)

  # Round to integers and maintain sparse format
  if (inherits(raw_counts, "dgCMatrix")) {
    raw_counts@x <- round(raw_counts@x)
  } else {
    raw_counts <- round(raw_counts)
  }

  # Convert to dgCMatrix if needed
  if (!inherits(raw_counts, "dgCMatrix")) {
    raw_counts <- as(raw_counts, "dgCMatrix")
  }

  # Restore dimnames
  dimnames(raw_counts) <- dimnames(normalized_data)

  return(raw_counts)
}

# =============================================================================
# Directory Setup
# =============================================================================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(OUTPUT_DIR, "cleaned_samples"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(OUTPUT_DIR, "qc_stats"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(OUTPUT_DIR, "qc_plots"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(OUTPUT_DIR, "merged"),
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# Seurat object structure snapshot (Seurat v5/v4 compatible)
# =============================================================================
snapshot_seurat_structure <- function(obj, sample_name, outdir, assay = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }

  f_txt <- file.path(outdir, paste0(sample_name, "_structure.txt"))
  f_meta <- file.path(outdir, paste0(sample_name, "_meta_schema.csv"))
  f_asst <- file.path(outdir, paste0(sample_name, "_assay_schema.csv"))

  # --- meta schema ---
  meta <- obj@meta.data
  meta_schema <- data.frame(
    column = colnames(meta),
    class = vapply(
      meta,
      function(x) paste(class(x), collapse = "|"),
      character(1)
    ),
    n_na = vapply(meta, function(x) sum(is.na(x)), integer(1)),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(meta_schema, f_meta)

  # --- assay schema ---
  assays <- get_assay_names_safe(obj)
  assay_rows <- lapply(assays, function(a) {
    ao <- obj[[a]]
    layers <- tryCatch(Layers(ao), error = function(e) character(0))
    data.frame(
      assay = a,
      assay_class = class(ao)[1],
      layers = if (length(layers)) paste(layers, collapse = ",") else "",
      stringsAsFactors = FALSE
    )
  })
  assay_schema <- data.table::rbindlist(assay_rows, fill = TRUE)
  data.table::fwrite(assay_schema, f_asst)

  # --- counts quick sanity ---
  get_counts_quick <- function(o, a) {
    tryCatch(
      LayerData(o, assay = a, layer = "counts"),
      error = function(e1) {
        tryCatch(
          GetAssayData(o, assay = a, slot = "counts"),
          error = function(e2) NULL
        )
      }
    )
  }
  cnt <- get_counts_quick(obj, assay)
  cnt_info <- if (is.null(cnt)) {
    "counts: <NOT FOUND>"
  } else {
    nnz <- if (inherits(cnt, "dgCMatrix")) length(cnt@x) else NA_integer_
    has_na <- if (inherits(cnt, "dgCMatrix")) anyNA(cnt@x) else anyNA(cnt)
    has_rownames <- !is.null(rownames(cnt))
    has_colnames <- !is.null(colnames(cnt))
    is_integer <- is_integer_like_counts(cnt)

    sprintf(
      "counts: %d genes x %d cells | class=%s | nnz=%s | anyNA=%s | rownames=%s | colnames=%s | integer_like=%s",
      nrow(cnt),
      ncol(cnt),
      class(cnt)[1],
      ifelse(is.na(nnz), "NA", as.character(nnz)),
      as.character(has_na),
      as.character(has_rownames),
      as.character(has_colnames),
      as.character(is_integer)
    )
  }

  # --- Idents schema ---
  idt <- Idents(obj)
  idt_info <- sprintf(
    "Idents: class=%s | n_levels=%d | head=%s",
    class(idt)[1],
    length(levels(idt)),
    paste(head(as.character(idt), 5), collapse = ",")
  )

  # --- write text snapshot ---
  txt <- c(
    sprintf("sample: %s", sample_name),
    sprintf("Seurat class: %s", class(obj)[1]),
    sprintf("cells: %d | features: %d", ncol(obj), nrow(obj)),
    sprintf("DefaultAssay: %s", DefaultAssay(obj)),
    sprintf("Assays: %s", paste(assays, collapse = ",")),
    sprintf("assay_selected_for_snapshot: %s", assay),
    cnt_info,
    idt_info,
    sprintf("meta cols: %d", ncol(obj@meta.data))
  )

  writeLines(txt, f_txt)
  invisible(TRUE)
}

# =============================================================================
# Dynamic QC Column Names (for non-RNA assays)
# =============================================================================
get_qc_cols <- function(assay) {
  list(
    nFeature = paste0("nFeature_", assay),
    nCount = paste0("nCount_", assay)
  )
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(...)
  cat("\n")
}

# =============================================================================
# Safe assay name getter (Seurat v5 compatible)
# =============================================================================
get_assay_names_safe <- function(obj) {
  nm <- tryCatch(names(obj@assays), error = function(e) NULL)
  if (is.null(nm) || length(nm) == 0) {
    ao <- tryCatch(SeuratObject::Assays(obj), error = function(e) NULL)
    nm <- tryCatch(names(ao), error = function(e) character(0))
    if (length(nm) == 0) {
      nm <- tryCatch(as.character(ao), error = function(e) character(0))
    }
  }
  nm <- as.character(nm)
  nm <- nm[!is.na(nm) & nzchar(nm)]
  nm
}

# Pick preferred assay
pick_assay <- function(obj, preferred = "RNA") {
  nm <- get_assay_names_safe(obj)
  if (length(nm) == 0) {
    stop("No assays detected in Seurat object (names(obj@assays) is empty).")
  }

  preferred <- as.character(preferred)[1]
  if (!is.na(preferred) && nzchar(preferred) && preferred %in% nm) {
    return(preferred)
  }

  active <- tryCatch(obj@active.assay, error = function(e) NA_character_)
  active <- as.character(active)[1]
  if (!is.na(active) && nzchar(active) && active %in% nm) {
    return(active)
  }

  nm[1]
}

# =============================================================================
# Fix missing dimnames in counts matrix
# =============================================================================
fix_dimnames_counts <- function(counts_mat, obj, assay = "RNA") {
  if (is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
    feats <- tryCatch(
      SeuratObject::Features(obj, assay = assay),
      error = function(e) rownames(obj)
    )
    cells <- colnames(obj)
    dimnames(counts_mat) <- list(feats, cells)
  }
  counts_mat
}

# =============================================================================
# Check if counts matrix contains integer-like values
# =============================================================================
is_integer_like_counts <- function(m, tol = 1e-6, max_check = 5e5) {
  if (!inherits(m, "dgCMatrix")) {
    return(FALSE)
  }
  x <- m@x
  if (length(x) == 0) {
    return(TRUE)
  }
  if (length(x) > max_check) {
    x <- sample(x, max_check)
  }
  frac_nonint <- mean(abs(x - round(x)) > tol)
  frac_nonint < 0.01
}

# =============================================================================
# Seurat v5 Layers Handler
# =============================================================================
ensure_joined_layers <- function(seurat_obj, assay = "RNA", verbose = TRUE) {
  assay_obj <- seurat_obj[[assay]]

  if (!inherits(assay_obj, "Assay5")) {
    return(seurat_obj)
  }

  ly <- tryCatch(Layers(assay_obj), error = function(e) character(0))
  if (length(ly) <= 1) {
    return(seurat_obj)
  }

  if (verbose) {
    log_msg(sprintf(
      "  [Layers] Detected %d layers in %s → JoinLayers()",
      length(ly),
      assay
    ))
  }
  seurat_obj[[assay]] <- JoinLayers(seurat_obj[[assay]])
  return(seurat_obj)
}

# =============================================================================
# Assay5 → Assay Converter (for DoubletFinder compatibility)
# =============================================================================
rebuild_assay_from_scratch <- function(seurat_obj, assay = "RNA") {
  cat("  [Assay] Checking type...\n")

  if (!inherits(seurat_obj[[assay]], "Assay5")) {
    cat("    Already Assay (not Assay5), no rebuild needed\n")
    return(seurat_obj)
  }

  seurat_obj <- ensure_joined_layers(seurat_obj, assay = assay, verbose = TRUE)

  cat("    Extracting counts...\n")
  counts_data <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) {
      tryCatch(
        GetAssayData(seurat_obj, assay = assay, slot = "counts"),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(counts_data)) {
    stop("Unable to extract counts after JoinLayers()")
  }

  counts_data <- fix_dimnames_counts(counts_data, seurat_obj, assay = assay)

  cat(sprintf(
    "    Counts: %d genes × %d cells\n",
    nrow(counts_data),
    ncol(counts_data)
  ))

  cat("    Extracting/calculating data layer...\n")
  data_data <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )

  new_assay <- CreateAssayObject(counts = counts_data)

  if (!is.null(data_data) && nrow(data_data) > 0) {
    new_assay <- SetAssayData(new_assay, slot = "data", new.data = data_data)
    cat("    ✓ Data layer preserved\n")
  } else {
    cat("    Data layer missing → Recalculating NormalizeData()...\n")

    temp_obs <- seurat_obj@meta.data
    cat("    Converting factor columns to character...\n")
    factor_cols <- c()
    for (col in colnames(temp_obs)) {
      if (is.factor(temp_obs[[col]])) {
        temp_obs[[col]] <- as.character(temp_obs[[col]])
        factor_cols <- c(factor_cols, col)
      }
    }
    if (length(factor_cols) > 0) {
      cat(sprintf(
        "      Converted %d columns: %s\n",
        length(factor_cols),
        paste(head(factor_cols, 5), collapse = ", ")
      ))
    }

    cat("    Creating temporary Seurat object...\n")
    tmp <- CreateSeuratObject(counts = counts_data)

    cat("    Transferring metadata...\n")
    counts_cells <- colnames(counts_data)
    missing_meta <- setdiff(counts_cells, rownames(temp_obs))
    if (length(missing_meta) > 0) {
      cat("    WARNING: Metadata rownames don't fully match cell names\n")
      cat("    Attempting to align by available cells...\n")
    }

    available_cells <- counts_cells[counts_cells %in% rownames(temp_obs)]
    if (length(available_cells) == 0) {
      stop(
        "FATAL: No overlap between metadata rownames and cell names. Cannot rebuild assay safely."
      )
    }

    if (length(available_cells) < length(counts_cells)) {
      cat(sprintf(
        "    Aligned %d/%d cells (dropping %d lacking metadata)\n",
        length(available_cells),
        length(counts_cells),
        length(counts_cells) - length(available_cells)
      ))
    }

    counts_data <- counts_data[, available_cells, drop = FALSE]
    temp_obs_aligned <- temp_obs[available_cells, , drop = FALSE]
    tmp@meta.data <- temp_obs_aligned

    cat("    Running NormalizeData...\n")
    tmp <- NormalizeData(tmp, verbose = FALSE)

    cmd_key <- paste0("NormalizeData.", assay)
    if (!is.null(tmp@commands) && !is.null(tmp@commands[[cmd_key]])) {
      seurat_obj@commands[[cmd_key]] <- tmp@commands[[cmd_key]]
      cat(sprintf("    ✓ Copied NormalizeData command: %s\n", cmd_key))
    }

    new_assay <- SetAssayData(
      new_assay,
      slot = "data",
      new.data = GetAssayData(tmp, slot = "data")
    )
    rm(tmp, temp_obs, temp_obs_aligned)
    gc(verbose = FALSE)
    cat("    ✓ Normalization completed\n")
  }

  seurat_obj[[assay]] <- new_assay

  cat(sprintf("  ✓ Assay rebuilt as: %s\n", class(seurat_obj[[assay]])[1]))
  return(seurat_obj)
}

# =============================================================================
# Sample Column Fallback
# =============================================================================
ensure_sample_column <- function(seurat_obj, sample_column = "sample") {
  md <- seurat_obj@meta.data

  if (sample_column %in% colnames(md) && !all(is.na(md[[sample_column]]))) {
    Idents(seurat_obj) <- sample_column
    return(seurat_obj)
  }

  if ("orig.ident" %in% colnames(md)) {
    seurat_obj[[sample_column]] <- as.character(md$orig.ident)
    log_msg(sprintf(
      "  [Meta] Created '%s' column from orig.ident",
      sample_column
    ))
  } else {
    seurat_obj[[sample_column]] <- "all"
    log_msg(sprintf(
      "  [Meta] Created '%s' column with default 'all'",
      sample_column
    ))
  }

  Idents(seurat_obj) <- sample_column

  return(seurat_obj)
}

# =============================================================================
# Counts Matrix Getter (robust across v5/v4)
# =============================================================================
get_counts_matrix <- function(seurat_obj, assay = "RNA") {
  m <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e1) {
      tryCatch(
        GetAssayData(seurat_obj, assay = assay, slot = "counts"),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(m)) {
    stop("Cannot extract counts matrix.")
  }
  if (!inherits(m, "dgCMatrix")) {
    m <- as(m, "dgCMatrix")
  }

  m <- fix_dimnames_counts(m, seurat_obj, assay = assay)

  return(m)
}

# Drop NA-corrupted cells
drop_na_corrupted_cells <- function(counts) {
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  if (!anyNA(counts@x)) {
    return(list(dropped = 0L, bad_cells = character(0)))
  }

  cs <- suppressWarnings(Matrix::colSums(counts))
  bad <- which(is.na(cs))
  bad_cells <- colnames(counts)[bad]

  return(list(dropped = length(bad), bad_cells = bad_cells))
}

# =============================================================================
# QC Metrics Calculation
# =============================================================================
calc_qc_metrics <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  qc_cols <- get_qc_cols(assay)
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MT_PATTERN)
  }
  if (!"percent.rb" %in% colnames(obj@meta.data)) {
    obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = RB_PATTERN)
  }

  missing_qc <- setdiff(unlist(qc_cols, use.names = FALSE), colnames(obj@meta.data))
  if (length(missing_qc) > 0) {
    counts_mat <- get_counts_matrix(obj, assay = assay)
    if (!qc_cols$nCount %in% colnames(obj@meta.data)) {
      obj@meta.data[[qc_cols$nCount]] <- Matrix::colSums(counts_mat)
    }
    if (!qc_cols$nFeature %in% colnames(obj@meta.data)) {
      obj@meta.data[[qc_cols$nFeature]] <- Matrix::colSums(counts_mat > 0)
    }
  }
  return(obj)
}

# =============================================================================
# Basic Threshold Filtering
# =============================================================================
filter_basic_thresholds <- function(obj, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  qc_cols <- get_qc_cols(assay)

  before <- ncol(obj)

  if (!qc_cols$nFeature %in% colnames(obj@meta.data)) {
    stop(sprintf("QC column '%s' not found in metadata", qc_cols$nFeature))
  }

  obj <- subset(
    obj,
    subset = obj@meta.data[[qc_cols$nFeature]] >= NFEATURE_MIN &
      obj@meta.data[[qc_cols$nFeature]] <= NFEATURE_MAX &
      percent.mt < MT_MAX &
      percent.rb < RB_MAX
  )
  after <- ncol(obj)
  list(obj = obj, before = before, after = after)
}

# =============================================================================
# Regression Outlier Detection
# =============================================================================
regression_outliers <- function(obj, pred_level = 0.999, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  qc_cols <- get_qc_cols(assay)

  if (
    !qc_cols$nCount %in% colnames(obj@meta.data) ||
      !qc_cols$nFeature %in% colnames(obj@meta.data)
  ) {
    stop(sprintf(
      "QC columns '%s' or '%s' not found",
      qc_cols$nCount,
      qc_cols$nFeature
    ))
  }

  nCount <- obj@meta.data[[qc_cols$nCount]]
  nFeature <- obj@meta.data[[qc_cols$nFeature]]

  log_counts <- log10(nCount + 1)
  log_features <- log10(nFeature + 1)

  fit <- lm(log_features ~ log_counts)
  pred <- predict(
    fit,
    data.frame(log_counts = log_counts),
    interval = "prediction",
    level = pred_level
  )

  outliers <- log_features < pred[, "lwr"] | log_features > pred[, "upr"]
  keep <- colnames(obj)[!outliers]
  list(keep_cells = keep, outlier_rate = mean(outliers))
}

# =============================================================================
# QC Plot PDF
# =============================================================================
plot_qc_pdf <- function(obj, sample_name, outdir, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  qc_cols <- get_qc_cols(assay)

  pdf(
    file.path(outdir, paste0(sample_name, "_qc.pdf")),
    width = 14,
    height = 10
  )

  feats <- c(qc_cols$nFeature, qc_cols$nCount, "percent.mt", "percent.rb")
  feats <- feats[feats %in% colnames(obj@meta.data)]

  print(VlnPlot(
    obj,
    features = feats,
    ncol = min(4, length(feats)),
    pt.size = 0.1
  ))

  if (
    qc_cols$nCount %in%
      colnames(obj@meta.data) &&
      qc_cols$nFeature %in% colnames(obj@meta.data)
  ) {
    print(
      FeatureScatter(
        obj,
        feature1 = qc_cols$nCount,
        feature2 = qc_cols$nFeature,
        pt.size = 0.6
      ) +
        ggtitle(sprintf("UMIs vs Genes (%s)", assay))
    )
  }
  if (
    qc_cols$nCount %in%
      colnames(obj@meta.data) &&
      "percent.mt" %in% colnames(obj@meta.data)
  ) {
    print(
      FeatureScatter(
        obj,
        feature1 = qc_cols$nCount,
        feature2 = "percent.mt",
        pt.size = 0.6
      ) +
        ggtitle(sprintf("UMIs vs MT%% (%s)", assay))
    )
  }

  dev.off()
}

# =============================================================================
# Ensure NormalizeData command exists
# =============================================================================
ensure_normalize_command <- function(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
) {
  DefaultAssay(obj) <- assay
  cmd_key <- paste0("NormalizeData.", assay)

  if (is.null(obj@commands) || is.null(obj@commands[[cmd_key]])) {
    obj <- NormalizeData(
      obj,
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = FALSE
    )
  }
  obj
}

# =============================================================================
# DoubletFinder (per sample)
# =============================================================================
run_doubletfinder_keep_singlets <- function(
  obj,
  assay = "RNA",
  doublet_rate = 0.06,
  n_pcs = 30
) {
  DefaultAssay(obj) <- assay

  n_cells <- ncol(obj)
  if (n_cells < 50) {
    obj$doublet_class <- "Skipped_too_few_cells"
    return(list(obj = obj, kept = n_cells, status = "Skipped"))
  }

  data_mat <- tryCatch(GetAssayData(obj, slot = "data"), error = function(e) {
    NULL
  })
  if (is.null(data_mat) || nrow(data_mat) == 0) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  obj <- ensure_normalize_command(obj, assay = assay)

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

  n_pcs_use <- min(n_pcs, max(10, n_cells - 1))
  obj <- RunPCA(obj, npcs = n_pcs_use, verbose = FALSE)

  nExp <- round(doublet_rate * n_cells)

  sweep.res <- paramSweep(obj, PCs = 1:n_pcs_use, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(
    bcmvn$BCmetric
  )])))
  if (is.na(pK) || pK < 0.01 || pK > 0.3) {
    pK <- 0.09
  }

  rm(sweep.res, sweep.stats, bcmvn)
  gc(verbose = FALSE)

  obj <- doubletFinder(
    obj,
    PCs = 1:n_pcs_use,
    pN = 0.25,
    pK = pK,
    nExp = nExp,
    sct = FALSE
  )

  df_class_col <- grep(
    "DF.classifications",
    colnames(obj@meta.data),
    value = TRUE
  )
  if (length(df_class_col) == 0) {
    obj$doublet_class <- "DF_failed"
    return(list(obj = obj, kept = n_cells, status = "Failed"))
  }

  df_class_col <- df_class_col[length(df_class_col)]
  cls <- obj@meta.data[[df_class_col]]
  obj$doublet_class <- ifelse(cls == "Singlet", "Singlet", "Doublet")

  keep_cells <- colnames(obj)[obj$doublet_class == "Singlet"]
  obj2 <- subset(obj, cells = keep_cells)

  return(list(obj = obj2, kept = length(keep_cells), status = "Success"))
}

# =============================================================================
# ⭐ NEW v3: DecontX with De-lognormalization AND Storage Support
# =============================================================================
run_decontx <- function(obj, assay = "RNA", max_iter = 500,
                        save_restored = FALSE, restored_assay_name = "RNA_restored") {
  DefaultAssay(obj) <- assay
  cnt <- get_counts_matrix(obj, assay = assay)

  # Track whether we restored counts
  restored_counts_matrix <- NULL

  # Check if counts are integer-like
  if (!is_integer_like_counts(cnt)) {
    log_msg(
      "  [DecontX] ⚠️  Counts are NOT integer-like (normalized/scaled data detected)"
    )

    # Attempt de-lognormalization if enabled
    if (TRY_DENORMALIZE) {
      log_msg(
        "  [DecontX] Attempting to restore raw counts via de-lognormalization..."
      )

      restored_counts <- tryCatch(
        {
          restore_counts_from_seurat(
            obj,
            assay = assay,
            scale_factor = DENORM_SCALE_FACTOR
          )
        },
        error = function(e) {
          log_msg(sprintf(
            "  [DecontX] De-lognormalization failed: %s",
            e$message
          ))
          NULL
        }
      )

      if (!is.null(restored_counts)) {
        # Validate restored counts
        if (is_integer_like_counts(restored_counts)) {
          log_msg("  [DecontX] ✓ Successfully restored integer-like counts")
          cnt <- restored_counts
          restored_counts_matrix <- restored_counts  # Save for storage
          obj$data_type <- "restored_counts"
        } else {
          log_msg(
            "  [DecontX] ✗ Restored counts are still not integer-like → SKIPPING"
          )
          obj$decontX_contamination <- NA_real_
          obj$data_type <- "denorm_failed"
          return(obj)
        }
      } else {
        log_msg("  [DecontX] De-lognormalization returned NULL → SKIPPING")
        obj$decontX_contamination <- NA_real_
        obj$data_type <- "denorm_error"
        return(obj)
      }
    } else {
      log_msg(
        "  [DecontX] DecontX requires raw UMI counts → SKIPPING (TRY_DENORMALIZE=FALSE)"
      )
      obj$decontX_contamination <- NA_real_
      obj$data_type <- "non_raw_counts"
      return(obj)
    }
  } else {
    obj$data_type <- "raw_counts"
  }

  # Drop NA-corrupted cells
  res_na <- drop_na_corrupted_cells(cnt)
  if (res_na$dropped > 0) {
    log_msg(sprintf(
      "  [DecontX] Dropped %d NA-corrupted cells before DecontX",
      res_na$dropped
    ))
    obj <- subset(obj, cells = setdiff(colnames(obj), res_na$bad_cells))
    cnt <- get_counts_matrix(obj, assay = assay)

    # Also subset restored_counts_matrix if it exists
    if (!is.null(restored_counts_matrix)) {
      keep_cells <- setdiff(colnames(restored_counts_matrix), res_na$bad_cells)
      restored_counts_matrix <- restored_counts_matrix[, keep_cells, drop = FALSE]
    }
  }

  # Run DecontX
  set.seed(12345)
  dx <- decontX(x = cnt, maxIter = max_iter)

  # Store contamination scores
  obj$decontX_contamination <- dx$contamination

  # ⭐ NEW v3: Store restored counts in new assay if requested
  if (save_restored && !is.null(restored_counts_matrix)) {
    log_msg(sprintf(
      "  [DecontX] Storing restored counts in assay '%s'",
      restored_assay_name
    ))

    # Create new assay with restored counts
    new_assay <- tryCatch(
      {
        CreateAssayObject(counts = restored_counts_matrix)
      },
      error = function(e) {
        log_msg(sprintf(
          "  [DecontX] WARNING: Failed to create assay '%s': %s",
          restored_assay_name,
          e$message
        ))
        NULL
      }
    )

    if (!is.null(new_assay)) {
      obj[[restored_assay_name]] <- new_assay
      log_msg(sprintf("  [DecontX] ✓ Restored counts stored in '%s' assay", restored_assay_name))
    }
  }

  return(obj)
}

# =============================================================================
# ⭐ NEW v4: Cell ID Helpers
# =============================================================================
ensure_cell_prefix <- function(obj, sample_name, stage = "processing") {
  prefix <- paste0(sample_name, "_")
  cells <- colnames(obj)
  has_prefix <- startsWith(cells, prefix)

  if (all(has_prefix)) {
    log_msg(sprintf(
      "  [CellID] Prefix already present for sample '%s' (%s stage)",
      sample_name,
      stage
    ))
    return(obj)
  }

  if (any(has_prefix)) {
    stop(sprintf(
      "Inconsistent cell prefixes detected for sample '%s' during %s stage.",
      sample_name,
      stage
    ))
  }

  log_msg(sprintf(
    "  [CellID] Adding prefix '%s' to %d cells (%s stage)",
    prefix,
    length(cells),
    stage
  ))

  new_cells <- paste0(prefix, cells)
  obj <- RenameCells(obj, new.names = new_cells)

  return(obj)
}

add_cell_ids <- function(obj, sample_name) {
  log_msg(sprintf("  [CellID] Ensuring cell ID prefix for sample '%s'", sample_name))
  obj <- ensure_cell_prefix(obj, sample_name, stage = "processing")
  return(obj)
}

# =============================================================================
# ⭐ NEW v4: Promote RNA_restored to RNA assay
# =============================================================================
promote_rna_restored <- function(obj, restored_assay_name = "RNA_restored") {
  assays <- get_assay_names_safe(obj)

  if (restored_assay_name %in% assays) {
    log_msg(sprintf("  [Assay] Found '%s' assay - promoting to 'RNA'", restored_assay_name))

    # Remove old RNA assay
    if ("RNA" %in% assays) {
      obj[["RNA"]] <- NULL
      log_msg("  [Assay] Removed old 'RNA' assay")
    }

    # Get the restored assay
    restored_assay <- obj[[restored_assay_name]]

    # Create new RNA assay from restored counts
    obj[["RNA"]] <- restored_assay

    # Remove the restored assay (now it's in RNA)
    obj[[restored_assay_name]] <- NULL

    # Set RNA as default
    DefaultAssay(obj) <- "RNA"

    log_msg("  [Assay] ✓ RNA_restored promoted to RNA and set as default")
  }

  return(obj)
}

# =============================================================================
# Per-Sample Processing Function
# =============================================================================
process_one_rds <- function(rds_path, sample_name = NULL) {
  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(rds_path))
  }
  log_msg(sprintf("=== Processing Sample: %s ===", sample_name))

  log_msg("  [Step 1/11] Reading RDS file...")
  obj <- readRDS(rds_path)
  log_msg(sprintf("    Cells: %d, Features: %d", ncol(obj), nrow(obj)))

  # Pick assay
  log_msg("  [Step 2/11] Selecting assay...")
  assay_use <- pick_assay(obj, ASSAY_PREFERRED)
  DefaultAssay(obj) <- assay_use
  log_msg(sprintf(
    "    Using assay: %s (class: %s)",
    assay_use,
    class(obj[[assay_use]])[1]
  ))

  # Snapshot structure
  snapshot_seurat_structure(
    obj,
    sample_name = sample_name,
    outdir = file.path(OUTPUT_DIR, "qc_stats", "structure_snapshots"),
    assay = assay_use
  )

  # Ensure sample metadata
  log_msg("  [Step 3/11] Setting up sample metadata...")
  log_msg(sprintf("    Setting sample = '%s'", sample_name))
  obj$sample <- as.character(sample_name)

  log_msg("    Resetting Idents to sample (avoid factor issues)...")
  Idents(obj) <- "sample"

  log_msg("    Running ensure_sample_column...")
  obj <- ensure_sample_column(obj, sample_column = "sample")

  # Handle Seurat v5 layers
  log_msg("  [Step 4/11] Handling Seurat v5 layers...")
  obj <- ensure_joined_layers(obj, assay = assay_use, verbose = TRUE)

  log_msg("  [Step 5/11] Rebuilding assay (v5 compatibility)...")
  obj <- rebuild_assay_from_scratch(obj, assay = assay_use)

  # Calculate QC metrics
  log_msg("  [Step 6/11] Calculating QC metrics...")
  obj <- calc_qc_metrics(obj, assay = assay_use)

  qc_rows <- list()

  # Stage 1: Basic Thresholds
  b1 <- ncol(obj)
  fb <- filter_basic_thresholds(obj, assay = assay_use)
  obj <- fb$obj
  a1 <- fb$after
  qc_rows[[length(qc_rows) + 1]] <- data.frame(
    sample = sample_name,
    step = "basic_thresholds",
    cells_before = b1,
    cells_after = a1,
    stringsAsFactors = FALSE
  )
  log_msg(sprintf("  [QC] Basic thresholds: %d → %d", b1, a1))

  if (ncol(obj) == 0) {
    log_msg("  [QC] No cells left after basic thresholds. Skipping sample.")
    return(list(obj = NULL, qc = rbindlist(qc_rows, fill = TRUE)))
  }

  # Stage 2: Regression Outliers
  if (DO_REGRESSION_OUTLIER && ncol(obj) >= 50) {
    b2 <- ncol(obj)
    reg <- regression_outliers(
      obj,
      pred_level = OUTLIER_PRED_LEVEL,
      assay = assay_use
    )
    obj <- subset(obj, cells = reg$keep_cells)
    a2 <- ncol(obj)
    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      sample = sample_name,
      step = "regression_outliers",
      cells_before = b2,
      cells_after = a2,
      stringsAsFactors = FALSE
    )
    log_msg(sprintf(
      "  [QC] Regression outliers: %d → %d (outlier_rate=%.3f)",
      b2,
      a2,
      reg$outlier_rate
    ))
  }

  # Optional: QC Plots
  if (PLOT_QC_PDF && ncol(obj) > 0) {
    plot_qc_pdf(
      obj,
      sample_name,
      file.path(OUTPUT_DIR, "qc_plots"),
      assay = assay_use
    )
  }

  # Stage 3: DoubletFinder
  if (DO_DOUBLET_FINDER && ncol(obj) > 0) {
    b3 <- ncol(obj)
    df_res <- tryCatch(
      run_doubletfinder_keep_singlets(
        obj,
        assay = assay_use,
        doublet_rate = DOUBLET_RATE,
        n_pcs = DOUBLET_PCS
      ),
      error = function(e) {
        log_msg(sprintf("  [DoubletFinder] ERROR: %s", e$message))
        obj$doublet_class <- "DF_error"
        list(obj = obj, kept = ncol(obj), status = "Error")
      }
    )
    obj <- df_res$obj
    a3 <- ncol(obj)
    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      sample = sample_name,
      step = paste0("doubletfinder_", df_res$status),
      cells_before = b3,
      cells_after = a3,
      stringsAsFactors = FALSE
    )
    log_msg(sprintf("  [DoubletFinder] %s: %d → %d", df_res$status, b3, a3))
  }

  # Stage 4: DecontX
  if (DO_DECONTX && ncol(obj) > 0) {
    b4 <- ncol(obj)
    obj <- tryCatch(
      run_decontx(
        obj,
        assay = assay_use,
        max_iter = DECONTX_MAX_ITER,
        save_restored = SAVE_RESTORED_COUNTS,
        restored_assay_name = RESTORED_ASSAY_NAME
      ),
      error = function(e) {
        log_msg(sprintf("  [DecontX] ERROR: %s", e$message))
        obj$decontX_contamination <- NA_real_
        obj
      }
    )

    # Filter by contamination
    if (
      FILTER_BY_DECONTX && "decontX_contamination" %in% colnames(obj@meta.data)
    ) {
      cont <- obj$decontX_contamination
      if (all(is.na(cont))) {
        log_msg(
          "  [DecontX] Contamination is all NA → skip contamination filtering"
        )
      } else {
        keep <- which(!is.na(cont) & cont <= DECONTX_CONTAM_MAX)
        obj <- obj[, keep, drop = FALSE]
      }
    }

    a4 <- ncol(obj)
    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      sample = sample_name,
      step = "decontx_and_filter",
      cells_before = b4,
      cells_after = a4,
      stringsAsFactors = FALSE
    )
    log_msg(sprintf("  [DecontX] Filter: %d → %d", b4, a4))
  }

  # ⭐ NEW v4: Promote RNA_restored to RNA if it exists
  log_msg("  [Step 7/11] Checking for RNA_restored assay...")
  obj <- promote_rna_restored(obj, restored_assay_name = RESTORED_ASSAY_NAME)

  # ⭐ NEW v4: Add cell IDs BEFORE saving
  log_msg("  [Step 8/11] Adding cell ID prefixes...")
  obj <- add_cell_ids(obj, sample_name = sample_name)

  # Final Diet (keep only RNA assay)
  log_msg("  [Step 9/11] Diet Seurat (keeping only RNA assay)...")
  obj <- DietSeurat(
    obj,
    assays = "RNA",
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL
  )

  log_msg(sprintf("  [Step 10/11] Final object: %d cells, %d genes", ncol(obj), nrow(obj)))

  return(list(obj = obj, qc = rbindlist(qc_rows, fill = TRUE)))
}

# =============================================================================
# ⭐ NEW v4: Load and prepare object for merge
# =============================================================================
prepare_merge_object <- function(obj, sample_name, context = "merge") {
  obj <- promote_rna_restored(obj, restored_assay_name = RESTORED_ASSAY_NAME)
  obj <- ensure_cell_prefix(obj, sample_name, stage = context)
  DefaultAssay(obj) <- "RNA"
  obj <- ensure_sample_column(obj, sample_column = "sample")
  return(obj)
}

load_for_merge <- function(rds_path, sample_name = NULL) {
  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(rds_path))
  }

  log_msg(sprintf("  [Merge] Loading sample: %s", sample_name))

  obj <- readRDS(rds_path)

  log_msg(sprintf("    Loaded: %d cells, %d genes", ncol(obj), nrow(obj)))

  obj <- prepare_merge_object(obj, sample_name, context = "merge(load)")

  return(obj)
}

# =============================================================================
# ⭐ NEW v4: Incremental merge function
# =============================================================================
incremental_merge <- function(merge_inputs) {
  n_inputs <- length(merge_inputs)
  log_msg(sprintf("=== Starting Incremental Merge of %d Samples ===", n_inputs))

  if (n_inputs == 0) {
    stop("No samples to merge")
  }

  load_entry <- function(entry, idx, total) {

    if (is.null(entry$sample)) {
      stop(sprintf("Merge entry at index %d is missing 'sample' field", idx))
    }

    if (!is.null(entry$path)) {
      log_msg(sprintf("[%d/%d] Loading sample from disk: %s", idx, total, entry$sample))
      obj <- load_for_merge(entry$path, entry$sample)
    } else if (!is.null(entry$object)) {
      log_msg(sprintf("[%d/%d] Using in-memory sample: %s", idx, total, entry$sample))
      obj <- entry$object
      if (is.null(obj)) {
        stop(sprintf("Merge entry for sample '%s' has NULL object", entry$sample))
      }
      log_msg(sprintf("    In-memory object: %d cells, %d genes", ncol(obj), nrow(obj)))
      obj <- prepare_merge_object(obj, entry$sample, context = "merge(memory)")
    } else {
      stop(sprintf("Merge entry for sample '%s' has neither path nor object", entry$sample))
    }

    obj
  }

  if (n_inputs == 1) {
    log_msg("Only 1 sample, loading directly without merge")
    merged_obj <- load_entry(merge_inputs[[1]], 1, n_inputs)
    gc(verbose = FALSE)
    return(merged_obj)
  }

  merged_obj <- load_entry(merge_inputs[[1]], 1, n_inputs)
  gc(verbose = FALSE)

  for (i in 2:n_inputs) {
    next_obj <- load_entry(merge_inputs[[i]], i, n_inputs)

    log_msg(sprintf("    Merging: %d cells + %d cells", ncol(merged_obj), ncol(next_obj)))
    merged_obj <- merge(x = merged_obj, y = next_obj)

    rm(next_obj)
    gc(verbose = FALSE)

    log_msg(sprintf("    Cumulative total: %d cells", ncol(merged_obj)))
  }

  log_msg(sprintf("✓ Merge complete: %d total cells", ncol(merged_obj)))

  return(merged_obj)
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================
log_msg("╔══════════════════════════════════════════════════════════════╗")
log_msg("║  RDS Folder QC → DoubletFinder → DecontX → Incremental Merge ║")
log_msg("║  VERSION: v4_MEMORY_OPTIMIZED (Lazy Load + Incremental Merge)║")
log_msg("╚══════════════════════════════════════════════════════════════╝")

rds_files <- list.files(INPUT_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
  stop("No .rds files found in INPUT_DIR: ", INPUT_DIR)
}

log_msg(sprintf("Found %d RDS files in: %s", length(rds_files), INPUT_DIR))
log_msg(sprintf("Output directory: %s", OUTPUT_DIR))
log_msg(sprintf("De-lognormalization: %s", ifelse(TRY_DENORMALIZE, "ENABLED", "DISABLED")))
log_msg(sprintf("Save restored counts: %s (assay: %s)",
                ifelse(SAVE_RESTORED_COUNTS, "YES", "NO"),
                RESTORED_ASSAY_NAME))
log_msg(sprintf("RNA_restored priority: YES (will be promoted to RNA)"))
log_msg(sprintf("Skip existing cleaned files: %s", ifelse(SKIP_EXISTING_CLEANED && !FORCE_REPROCESS, "YES", "NO")))
log_msg(sprintf("Force reprocess: %s", ifelse(FORCE_REPROCESS, "YES", "NO")))
log_msg(sprintf("Memory optimization: ENABLED (lazy load + incremental merge)"))
log_msg("")

# ⭐ NEW v4: Track merge inputs (paths or in-memory objects)
merge_inputs <- list()
sample_names_valid <- character()
qc_all <- list()

record_merge_input <- function(sample_name, path = NULL, object = NULL) {
  # #region agent log
  log_entry <- list(
    sessionId = "debug-session",
    runId = "pre-fix",
    hypothesisId = "B",
    location = "rds_folder_qc_smartmerge_20251216_v4_memory_optimized.R:1419",
    message = "record_merge_input entry",
    data = list(
      sample = sample_name,
      has_path = !is.null(path),
      has_object = !is.null(object),
      path = if (!is.null(path)) path else NULL,
      object_size = if (!is.null(object)) sprintf("%d cells x %d genes", ncol(object), nrow(object)) else NULL
    ),
    timestamp = as.numeric(Sys.time()) * 1000
  )
  write(jsonlite::toJSON(log_entry, auto_unbox = TRUE), file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion

  if (is.null(sample_name) || nchar(sample_name) == 0) {
    stop("sample_name cannot be NULL or empty")
  }
  if (is.null(path) && is.null(object)) {
    stop(sprintf("Merge input for sample '%s' must have a source", sample_name))
  }
  if (!is.null(path) && !is.null(object)) {
    stop(sprintf("Sample '%s' cannot supply both path and object", sample_name))
  }
  if (!is.null(path) && !file.exists(path)) {
    stop(sprintf("Path does not exist for sample '%s': %s", sample_name, path))
  }
  if (!is.null(object) && !inherits(object, "Seurat")) {
    stop(sprintf("Object for sample '%s' is not a Seurat object", sample_name))
  }

  merge_inputs <<- append(
    merge_inputs,
    list(list(sample = sample_name, path = path, object = object))
  )
  sample_names_valid <<- c(sample_names_valid, sample_name)

  # #region agent log
  log_entry <- list(
    sessionId = "debug-session",
    runId = "pre-fix",
    hypothesisId = "B",
    location = "rds_folder_qc_smartmerge_20251216_v4_memory_optimized.R:1450",
    message = "record_merge_input exit",
    data = list(
      sample = sample_name,
      total_inputs = length(merge_inputs),
      total_names = length(sample_names_valid)
    ),
    timestamp = as.numeric(Sys.time()) * 1000
  )
  write(jsonlite::toJSON(log_entry, auto_unbox = TRUE), file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion
}

for (i in seq_along(rds_files)) {
  sample_name <- tools::file_path_sans_ext(basename(rds_files[i]))

  log_msg(sprintf(
    "(%d/%d) Processing: %s",
    i,
    length(rds_files),
    sample_name
  ))

  # ⭐ Check if cleaned RDS already exists
  out_rds <- file.path(
    OUTPUT_DIR,
    "cleaned_samples",
    paste0(sample_name, "_cleaned.rds")
  )

  if (SKIP_EXISTING_CLEANED && !FORCE_REPROCESS && file.exists(out_rds)) {
    log_msg(sprintf("  ⏭️  SKIPPED: Cleaned file already exists: %s", basename(out_rds)))

    # ⭐ NEW v4: Just record the path, don't load the object
    record_merge_input(sample_name, path = out_rds)

    # Create a QC record for skipped sample
    qc_skipped <- data.frame(
      sample = sample_name,
      step = "LOADED_FROM_EXISTING",
      cells_before = NA,
      cells_after = NA,
      stringsAsFactors = FALSE
    )
    qc_all[[i]] <- qc_skipped

    log_msg(sprintf("  ✓ Path recorded for later merge"))
    log_msg("")
    next  # Skip to next sample
  }

  # Process the sample (either new or reprocess)
  res <- tryCatch(
    {
      process_one_rds(rds_files[i], sample_name = sample_name)
    },
    error = function(e) {
      log_msg(sprintf("  ✗ FATAL ERROR processing this sample: %s", e$message))
      log_msg("  Skipping this sample and continuing with others...")

      qc_error <- data.frame(
        sample = sample_name,
        step = "FATAL_ERROR",
        cells_before = NA,
        cells_after = 0,
        stringsAsFactors = FALSE
      )
      list(obj = NULL, qc = qc_error)
    }
  )

  qc_all[[i]] <- res$qc

  if (!is.null(res$obj)) {
    cleaned_obj <- res$obj

    if (SAVE_INTERMEDIATE) {
      log_msg("  [Step 11/11] Saving cleaned object...")
      saveRDS(cleaned_obj, out_rds)
      log_msg(sprintf("  ✓ Saved: %s", basename(out_rds)))

      record_merge_input(sample_name, path = out_rds)

      rm(cleaned_obj)
      gc(verbose = FALSE)
      log_msg("  ✓ Memory cleared")
    } else {
      log_msg("  [Step 11/11] Keeping cleaned object in memory (SAVE_INTERMEDIATE=FALSE)")
      record_merge_input(sample_name, object = cleaned_obj)
    }

    rm(res)
    gc(verbose = FALSE)
  } else {
    log_msg("  ⚠ Sample dropped (no cells left or fatal error)")
  }

  log_msg("")
}

# Save QC stats
qc_df <- rbindlist(qc_all, fill = TRUE)
fwrite(qc_df, file.path(OUTPUT_DIR, "qc_stats", "per_sample_qc_stats.csv"))
log_msg("✓ Saved QC stats: qc_stats/per_sample_qc_stats.csv")
log_msg("")

# #region agent log
log_entry <- list(
  sessionId = "debug-session",
  runId = "pre-fix",
  hypothesisId = "C",
  location = "rds_folder_qc_smartmerge_20251216_v4_memory_optimized.R:1526",
  message = "Before merge check",
  data = list(
    n_merge_inputs = length(merge_inputs),
    n_sample_names = length(sample_names_valid),
    merge_inputs_summary = if (length(merge_inputs) > 0) {
      vapply(merge_inputs, function(x) {
        if (!is.null(x$path)) "path" else if (!is.null(x$object)) "object" else "unknown"
      }, character(1))
    } else character(0)
  ),
  timestamp = as.numeric(Sys.time()) * 1000
)
write(jsonlite::toJSON(log_entry, auto_unbox = TRUE), file = "/home/h2048/.cursor/debug.log", append = TRUE)
# #endregion

if (length(merge_inputs) == 0) {
  stop("No samples left after QC. Nothing to merge.")
}

if (length(merge_inputs) != length(sample_names_valid)) {
  stop(sprintf("Mismatch: %d merge inputs vs %d sample names", length(merge_inputs), length(sample_names_valid)))
}

# =============================================================================
# ⭐ NEW v4: Incremental Smart Merge
# =============================================================================
log_msg(sprintf("=== Starting Smart Merge of %d Samples ===", length(merge_inputs)))
log_msg("Sample list:")
for (i in seq_along(sample_names_valid)) {
  log_msg(sprintf("  [%d] %s", i, sample_names_valid[i]))
}
log_msg("")

merged_obj <- incremental_merge(merge_inputs)

log_msg(sprintf("Merged object: %d cells, %d genes", ncol(merged_obj), nrow(merged_obj)))

# Ensure sample column
merged_obj <- ensure_sample_column(merged_obj, sample_column = "sample")

# =============================================================================
# Save merged object before gene filtering
# =============================================================================
log_msg("Saving merged object before gene filtering (merged/merged_seurat.rds)...")
saveRDS(merged_obj, file.path(OUTPUT_DIR, "merged", "merged_seurat.rds"))
log_msg("✓ Saved: merged/merged_seurat.rds")

# =============================================================================
# Gene Filtering (min.cells)
# =============================================================================
if (MIN_CELLS_PER_GENE > 0) {
  log_msg(sprintf("Filtering genes (min.cells = %d)...", MIN_CELLS_PER_GENE))

  DefaultAssay(merged_obj) <- "RNA"

  merged_obj <- ensure_joined_layers(
    merged_obj,
    assay = "RNA",
    verbose = TRUE
  )

  cnt <- get_counts_matrix(merged_obj, assay = "RNA")

  gene_ncells <- Matrix::rowSums(cnt > 0)
  keep_genes <- names(gene_ncells[gene_ncells >= MIN_CELLS_PER_GENE])

  before_genes <- nrow(merged_obj)
  merged_obj <- subset(merged_obj, features = keep_genes)
  after_genes <- nrow(merged_obj)

  log_msg(sprintf(
    "  Genes: %d → %d (removed %d rare genes)",
    before_genes,
    after_genes,
    before_genes - after_genes
  ))
} else {
  log_msg("Skipping gene filtering (MIN_CELLS_PER_GENE <= 0)")
}

# =============================================================================
# Save final merged object
# =============================================================================
log_msg("Saving merged object after final filtering (merged/merged_seurat_final.rds)...")
saveRDS(merged_obj, file.path(OUTPUT_DIR, "merged", "merged_seurat_final.rds"))
log_msg("✓ Saved: merged/merged_seurat_final.rds")

log_msg("")
log_msg("╔══════════════════════════════════════════════════════════╗")
log_msg("║                    PIPELINE COMPLETE                     ║")
log_msg("╚══════════════════════════════════════════════════════════╝")
log_msg(sprintf("Total cells in final object: %d", ncol(merged_obj)))
log_msg(sprintf("Total genes in final object: %d", nrow(merged_obj)))
log_msg(sprintf("Number of samples: %d", length(unique(merged_obj$sample))))
log_msg(sprintf("Default assay: %s", DefaultAssay(merged_obj)))

# Report on assays
assays_present <- get_assay_names_safe(merged_obj)
log_msg(sprintf("Assays in merged object: %s", paste(assays_present, collapse = ", ")))
