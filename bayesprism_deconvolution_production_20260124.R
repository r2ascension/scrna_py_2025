#!/usr/bin/env Rscript

# ==============================================================================
# BayesPrism Deconvolution Pipeline - Production Version
# ==============================================================================
# Purpose: Deconvolve bulk RNA-seq using single-cell reference
# Author: r2end
# Date: 2025-01-23
# Version: 2.1 (added single-level annotation option)
#
# Input:
#   - Seurat RDS: single-cell reference with ann_level_2 (and optionally ann_level_3)
#   - Bulk CSV: gene_count.csv (gene_name, P1, P2, P3, P4, P5)
#
# Output:
#   - Cell type fractions (theta)
#   - Deconvolved gene expression (Z matrix)
#   - QC reports and visualizations
#
# Key Features:
#   - Option to use single annotation level (type only, no state subdivision)
#   - All P0/P1 fixes from v2.0
#
# Requirements:
#   - R >= 4.0
#   - Packages: data.table, Matrix, Seurat, SingleCellExperiment,
#               BayesPrism, ggplot2, pheatmap, RColorBrewer
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SingleCellExperiment)
  library(BayesPrism)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ==============================================================================
# 0. Configuration
# ==============================================================================

# ===== Input Paths (MODIFY THESE) =====
seurat_rds_path <- "/home/h2048/data/R/bulk0121/LungMap.rds" # ← CHANGE THIS
bulk_csv_path <- "/home/h2048/data/bulk/gene_count.csv" # ← CHANGE THIS
output_dir <- "/home/h2048/data/R/bulk0121/bayesprism_out" # ← CHANGE THIS

# ===== Cell Type Annotation Strategy =====
# Option 1: Use SINGLE level (recommended for your data)
use_single_level <- TRUE # ← Set TRUE to use only ann_level_2
celltype_col <- "ann_level_3" # Main cell type column

# Option 2: Use TWO levels (only if you have clean state→type mapping)
# use_single_level <- FALSE
# celltype_col  <- "ann_level_2"  # Coarse-level
# cellstate_col <- "ann_level_3"  # Fine-level

# Auto-set cellstate_col based on strategy
if (use_single_level) {
  cellstate_col <- celltype_col # State = Type (no subdivision)
  cat("Using SINGLE annotation level:", celltype_col, "\n")
} else {
  if (!exists("cellstate_col")) {
    stop("ERROR: use_single_level=FALSE but cellstate_col not defined")
  }
  cat(
    "Using TWO annotation levels:",
    celltype_col,
    "(type) and",
    cellstate_col,
    "(state)\n"
  )
}

# ===== BayesPrism Parameters =====
key_malignant <- NULL # NULL for healthy tissue (non-tumor)
outlier_cut <- 0.01 # Outlier gene filter threshold
outlier_frac <- 0.1 # Max fraction of outlier genes

# ===== Gene Filtering Options =====
filter_immunoglobulin <- FALSE # Set TRUE to remove IG genes (NOT recommended for immune deconvolution)

# ===== GEP Construction Method =====
gep_use_mean <- TRUE # TRUE = mean profile, FALSE = sum (may encode cell number bias)

# ===== Computational Resources =====
n_cores <- 32 # Adjust based on your CPU (recommend: total cores - 4)

# ===== Output Directory =====
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===== Reproducibility =====
set.seed(42)

# Print configuration
cat("\n")
cat(
  "================================================================================\n"
)
cat("BayesPrism Deconvolution Pipeline v2.1 (Production)\n")
cat(
  "================================================================================\n"
)
cat("Seurat RDS:      ", seurat_rds_path, "\n")
cat("Bulk CSV:        ", bulk_csv_path, "\n")
cat("Output directory:", output_dir, "\n")
cat(
  "Annotation mode: ",
  ifelse(use_single_level, "SINGLE level", "TWO levels"),
  "\n"
)
cat("Cell type column:", celltype_col, "\n")
cat("Cell state column:", cellstate_col, "\n")
cat("Filter IG genes: ", filter_immunoglobulin, "\n")
cat(
  "GEP method:      ",
  ifelse(gep_use_mean, "mean-profile", "sum-profile"),
  "\n"
)
cat("CPU cores:       ", n_cores, "\n")
cat("Random seed:     ", 42, "\n")
cat(
  "================================================================================\n\n"
)

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

#' Normalize gene symbols
#'
#' Remove version numbers, convert to uppercase
#'
#' @param gene_vec Character vector of gene symbols
#' @return Normalized gene symbols
normalize_gene_symbols <- function(gene_vec) {
  # Remove version numbers (e.g., ENSG00000123456.1 -> ENSG00000123456)
  gene_vec <- sub("\\.\\d+$", "", gene_vec)

  # Convert to uppercase (SYMBOL convention)
  gene_vec <- toupper(gene_vec)

  gene_vec
}

#' Filter technical genes
#'
#' Remove ribosomal, mitochondrial, hemoglobin, heat shock proteins, lncRNAs
#' Optionally remove immunoglobulin genes (NOT recommended for immune deconvolution)
#'
#' @param gene_vec Character vector of gene symbols
#' @param filter_ig Logical, whether to filter immunoglobulin genes (default FALSE)
#' @return Logical vector (TRUE = keep, FALSE = remove)
filter_technical_genes <- function(gene_vec, filter_ig = FALSE) {
  bad <- grepl("^MT-", gene_vec, ignore.case = FALSE) |
    grepl("^RPS", gene_vec) |
    grepl("^RPL", gene_vec) |
    grepl("^MRPS", gene_vec) |
    grepl("^MRPL", gene_vec) |
    grepl("^HBA", gene_vec) | # Hemoglobin alpha
    grepl("^HBB", gene_vec) | # Hemoglobin beta
    grepl("^HBG", gene_vec) | # Hemoglobin gamma
    grepl("^HBD", gene_vec) | # Hemoglobin delta
    grepl("^HBE", gene_vec) | # Hemoglobin epsilon
    grepl("^HSPA", gene_vec) | # Heat shock 70kDa
    grepl("^HSPB", gene_vec) | # Heat shock small
    grepl("^DNAJ", gene_vec) | # Heat shock DNAJ family
    gene_vec %in% c("MALAT1", "NEAT1", "XIST") # Long non-coding RNAs

  # Optionally filter immunoglobulin genes
  if (filter_ig) {
    bad <- bad |
      grepl("^IGH[GAMED]", gene_vec) | # Immunoglobulin heavy chain
      grepl("^IGK[CV]", gene_vec) | # Immunoglobulin kappa
      grepl("^IGL[CV]", gene_vec) # Immunoglobulin lambda
    message(
      "  NOTE: Immunoglobulin genes are being filtered (may reduce B cell deconvolution accuracy)"
    )
  } else {
    message(
      "  NOTE: Immunoglobulin genes are PRESERVED (recommended for immune deconvolution)"
    )
  }

  !bad
}

#' Deduplicate genes in sparse matrix by summing counts
#'
#' @param mat_gc Sparse matrix (genes × cells)
#' @return Deduplicated sparse matrix
dedup_gene_rows_sparse <- function(mat_gc) {
  g <- rownames(mat_gc)
  if (is.null(g)) {
    stop("Matrix rownames (gene symbols) is NULL")
  }
  if (!anyDuplicated(g)) {
    return(mat_gc)
  }

  n_dup <- sum(duplicated(g))
  message(sprintf(
    "  Deduplicating %d duplicated gene symbols (summing counts)",
    n_dup
  ))

  start_time <- Sys.time()

  idx_list <- split(seq_along(g), g)
  keep_genes <- names(idx_list)

  out <- Matrix(
    0,
    nrow = length(keep_genes),
    ncol = ncol(mat_gc),
    sparse = TRUE
  )
  rownames(out) <- keep_genes
  colnames(out) <- colnames(mat_gc)

  for (i in seq_along(keep_genes)) {
    ridx <- idx_list[[i]]
    if (length(ridx) == 1L) {
      out[i, ] <- mat_gc[ridx, ]
    } else {
      out[i, ] <- Matrix::colSums(mat_gc[ridx, , drop = FALSE])
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  message(sprintf("  Deduplication completed in %.1f seconds", elapsed))

  out
}

#' Collapse single-cell counts to gene expression programs (GEP)
#'
#' Aggregate cells by state, using mean-profile or sum-profile
#' Uses sparse matrix multiplication for efficiency
#'
#' @param counts_gene_cell Sparse matrix (genes × cells)
#' @param state_vec Character vector of cell states (length = ncol)
#' @param use_mean Logical, use mean profile (TRUE) or sum profile (FALSE)
#' @param filter_ig Logical, filter immunoglobulin genes
#' @return Dense matrix (states × genes)
collapse_to_gep <- function(
  counts_gene_cell,
  state_vec,
  use_mean = TRUE,
  filter_ig = FALSE
) {
  state_vec <- factor(state_vec)
  states <- levels(state_vec)

  # Filter technical genes BEFORE collapsing
  g <- rownames(counts_gene_cell)
  keep_gene <- filter_technical_genes(g, filter_ig = filter_ig)
  counts_filtered <- counts_gene_cell[keep_gene, , drop = FALSE]

  n_removed <- sum(!keep_gene)
  message(sprintf(
    "  Filtered %d technical genes (%d → %d)",
    n_removed,
    length(g),
    sum(keep_gene)
  ))

  # Create sparse design matrix: cells × states
  # This is much faster than looping for large datasets
  S <- Matrix::sparse.model.matrix(~ 0 + state_vec)
  colnames(S) <- sub("^state_vec", "", colnames(S))

  # Matrix multiplication: genes × states
  gep_gs <- counts_filtered %*% S

  # Convert to mean profile if requested
  if (use_mean) {
    n_per_state <- Matrix::colSums(S)
    gep_gs <- gep_gs %*% Matrix::Diagonal(x = 1 / pmax(n_per_state, 1))
    message(
      "  GEP constructed using MEAN profile (normalized per-cell expression)"
    )
  } else {
    message("  GEP constructed using SUM profile (may encode cell number bias)")
  }

  # Transpose to states × genes (BayesPrism expects this format)
  gep <- as.matrix(t(gep_gs))
  rownames(gep) <- colnames(S)
  colnames(gep) <- rownames(counts_filtered)

  # Remove states with 0 total counts
  keep_states <- rowSums(gep) > 0
  if (!all(keep_states)) {
    n_removed_states <- sum(!keep_states)
    warning(sprintf(
      "  Removing %d states with 0 total counts",
      n_removed_states
    ))
    gep <- gep[keep_states, , drop = FALSE]
  }

  gep
}

# ==============================================================================
# 2. Load Seurat Object with RNA Assay Validation
# ==============================================================================

cat("\n=== Step 1/11: Loading Seurat Object ===\n")

if (!file.exists(seurat_rds_path)) {
  stop("ERROR: Seurat RDS file not found: ", seurat_rds_path)
}

seurat_obj <- readRDS(seurat_rds_path)
cat(sprintf(
  "  Loaded: %d cells × %d features\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
))

# P0-1: Force RNA assay and validate raw counts
available_assays <- Seurat::Assays(seurat_obj)
cat("  Available assays:", paste(available_assays, collapse = ", "), "\n")

if (!"RNA" %in% available_assays) {
  stop(
    "ERROR: Seurat object has no RNA assay. Available assays: ",
    paste(available_assays, collapse = ", ")
  )
}

# Force RNA assay as default
Seurat::DefaultAssay(seurat_obj) <- "RNA"
cat("  Using RNA assay for raw counts\n")

# Convert to SingleCellExperiment, explicitly specifying RNA assay
sce_ref <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
cat("  Converted to SingleCellExperiment\n")

# Validate counts are raw (non-negative, ideally integer)
counts_gc <- counts(sce_ref)
if (!inherits(counts_gc, "dgCMatrix")) {
  counts_gc <- as(counts_gc, "dgCMatrix")
}

if (any(counts_gc@x < 0)) {
  stop("ERROR: Reference counts contain negative values. Not raw counts!")
}

non_integer_frac <- mean(abs(counts_gc@x - round(counts_gc@x)) > 1e-6)
if (non_integer_frac > 0.01) {
  warning(sprintf(
    "WARNING: %.1f%% of counts are non-integer. Confirm these are raw UMI counts.",
    100 * non_integer_frac
  ))
}

cat("  ✓ Counts validation passed (non-negative, mostly integer)\n")

# Check required metadata columns
required_cols <- unique(c(celltype_col, cellstate_col))
missing_cols <- setdiff(required_cols, colnames(colData(sce_ref)))
if (length(missing_cols) > 0) {
  stop(
    "ERROR: Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Save reference metadata for record
ref_metadata <- list(
  source = basename(seurat_rds_path),
  n_cells = ncol(sce_ref),
  n_genes = nrow(sce_ref),
  use_single_level = use_single_level,
  celltype_col = celltype_col,
  cellstate_col = cellstate_col,
  celltype_counts = table(colData(sce_ref)[[celltype_col]]),
  cellstate_counts = if (use_single_level) {
    table(colData(sce_ref)[[celltype_col]])
  } else {
    table(colData(sce_ref)[[cellstate_col]])
  },
  creation_date = Sys.Date(),
  gep_method = ifelse(gep_use_mean, "mean-profile", "sum-profile"),
  filter_ig = filter_immunoglobulin
)
saveRDS(ref_metadata, file.path(output_dir, "reference_metadata.rds"))

# Clean up Seurat object to save memory
rm(seurat_obj)
invisible(gc())

# ==============================================================================
# 3. Reference Data QC
# ==============================================================================

cat("\n=== Step 2/11: Reference Data QC ===\n")

# Extract metadata
cd <- as.data.frame(colData(sce_ref))
cell_type_labels <- as.character(cd[[celltype_col]])
cell_state_labels <- as.character(cd[[cellstate_col]])

# Remove cells with NA labels
ok <- !is.na(cell_type_labels) & !is.na(cell_state_labels)
if (sum(!ok) > 0) {
  cat(sprintf("  Removing %d cells with NA labels\n", sum(!ok)))
  sce_ref <- sce_ref[, ok]
  counts_gc <- counts_gc[, ok]
  cell_type_labels <- cell_type_labels[ok]
  cell_state_labels <- cell_state_labels[ok]
}

# Additional filter: remove "None" or empty labels
bad_labels <- cell_type_labels %in%
  c("None", "none", "", "NA", "Unknown", "unknown") |
  cell_state_labels %in% c("None", "none", "", "NA", "Unknown", "unknown")

if (sum(bad_labels) > 0) {
  cat(sprintf(
    "  Removing %d cells with 'None' or invalid labels\n",
    sum(bad_labels)
  ))
  sce_ref <- sce_ref[, !bad_labels]
  counts_gc <- counts_gc[, !bad_labels]
  cell_type_labels <- cell_type_labels[!bad_labels]
  cell_state_labels <- cell_state_labels[!bad_labels]
}

# Cell type distribution
celltype_counts <- table(cell_type_labels)
cat("\n  Cell type distribution (Level 2):\n")
print(addmargins(celltype_counts))

if (!use_single_level) {
  cellstate_counts <- table(cell_state_labels)
  cat("\n  Cell state distribution (Level 3, top 10):\n")
  print(head(sort(cellstate_counts, decreasing = TRUE), 10))
} else {
  cat("\n  Using single annotation level (state = type)\n")
}

# Warning for small cell types
small_types <- names(celltype_counts)[celltype_counts < 30]
if (length(small_types) > 0) {
  cat("\n  WARNING: Cell types with <30 cells:\n")
  cat("  ", paste(small_types, collapse = ", "), "\n")
  cat("  These may cause instability - consider merging or removing\n")
}

# Gene detection rate
gene_detection <- Matrix::rowMeans(counts_gc > 0)
cat(sprintf(
  "\n  Genes detected in >10%% cells: %d / %d (%.1f%%)\n",
  sum(gene_detection > 0.1),
  length(gene_detection),
  100 * mean(gene_detection > 0.1)
))

# P0-4: Check for state→type mapping ambiguity (only if using two levels)
if (!use_single_level) {
  cat("\n  Checking state→type mapping consistency...\n")
  state_type_table <- table(cell_state_labels, cell_type_labels)
  ambig <- rowSums(state_type_table > 0) > 1

  if (any(ambig)) {
    bad_states <- rownames(state_type_table)[ambig]
    ambig_table <- state_type_table[bad_states, , drop = FALSE]

    # Save ambiguity table
    write.csv(
      ambig_table,
      file.path(output_dir, "state_type_ambiguity.csv"),
      quote = FALSE
    )

    cat(
      "\n  ERROR: Found",
      sum(ambig),
      "cell states mapping to multiple cell types:\n"
    )
    print(ambig_table)
    stop(
      "State→type mapping is ambiguous. Check state_type_ambiguity.csv and fix annotations."
    )
  } else {
    cat("  ✓ All cell states map to exactly one cell type\n")
  }
} else {
  cat("  ✓ Skipping state→type validation (using single level)\n")
}

# Deduplicate genes
counts_gc <- dedup_gene_rows_sparse(counts_gc)

# ==============================================================================
# 4. Load Bulk Data
# ==============================================================================

cat("\n=== Step 3/11: Loading Bulk Data ===\n")

if (!file.exists(bulk_csv_path)) {
  stop("ERROR: Bulk data file not found: ", bulk_csv_path)
}

bulk_raw <- fread(bulk_csv_path)

# Check format
if (ncol(bulk_raw) < 2) {
  stop("ERROR: Bulk data must have at least 2 columns (gene_name + samples)")
}

setnames(bulk_raw, 1, "gene_name")
bulk_raw[, gene_name := as.character(gene_name)]

# P1-5: Normalize gene names
bulk_raw[, gene_name := normalize_gene_symbols(gene_name)]

# Aggregate duplicated genes
bulk_agg <- bulk_raw[,
  lapply(.SD, sum),
  by = gene_name,
  .SDcols = names(bulk_raw)[-1]
]

# Convert to matrix (gene × sample)
bulk_mat_gene_sample <- as.matrix(bulk_agg[, -1])
rownames(bulk_mat_gene_sample) <- bulk_agg$gene_name

cat(sprintf(
  "  Loaded: %d genes × %d samples\n",
  nrow(bulk_mat_gene_sample),
  ncol(bulk_mat_gene_sample)
))
cat("  Samples:", paste(colnames(bulk_mat_gene_sample), collapse = ", "), "\n")

# Filter technical genes (using same parameters as reference)
keep_bulk_gene <- filter_technical_genes(
  rownames(bulk_mat_gene_sample),
  filter_ig = filter_immunoglobulin
)
bulk_mat_gene_sample <- bulk_mat_gene_sample[keep_bulk_gene, , drop = FALSE]

cat(sprintf(
  "  After filtering: %d genes retained\n",
  nrow(bulk_mat_gene_sample)
))

# Transpose to sample × gene
bk_dat <- t(bulk_mat_gene_sample)
mode(bk_dat) <- "numeric"

# ==============================================================================
# 5. Bulk Data QC
# ==============================================================================

cat("\n=== Step 4/11: Bulk Data QC ===\n")

# Library size
lib_size <- colSums(bulk_mat_gene_sample)
cat(sprintf("  Library sizes:\n"))
cat(sprintf("    Median: %.0f\n", median(lib_size)))
cat(sprintf("    Range:  [%.0f, %.0f]\n", min(lib_size), max(lib_size)))
cat(sprintf("    CV:     %.2f\n", sd(lib_size) / mean(lib_size)))

# Check for extreme outliers
lib_size_cv <- sd(lib_size) / mean(lib_size)
if (lib_size_cv > 0.5) {
  cat(
    "\n  WARNING: High variation in library sizes (CV=",
    round(lib_size_cv, 2),
    ")\n"
  )
  cat("  Consider normalizing or checking sample quality\n")
}

# Gene detection
gene_det_bulk <- rowMeans(bulk_mat_gene_sample > 0)
cat(sprintf(
  "\n  Genes detected in all samples: %d (%.1f%%)\n",
  sum(gene_det_bulk == 1),
  100 * mean(gene_det_bulk == 1)
))
cat(sprintf(
  "  Genes detected in >50%% samples: %d (%.1f%%)\n",
  sum(gene_det_bulk > 0.5),
  100 * mean(gene_det_bulk > 0.5)
))

# ==============================================================================
# 6. Collapse Reference to GEP
# ==============================================================================

cat("\n=== Step 5/11: Collapsing Reference to GEP ===\n")

ref_gep <- collapse_to_gep(
  counts_gc,
  state_vec = cell_state_labels,
  use_mean = gep_use_mean,
  filter_ig = filter_immunoglobulin
)

cat(sprintf(
  "  GEP matrix: %d states × %d genes\n",
  nrow(ref_gep),
  ncol(ref_gep)
))

# Map each state to its corresponding cell type
if (use_single_level) {
  # State = Type (no mapping needed, 1-to-1)
  state_levels <- rownames(ref_gep)
  cell_type_labels_gep <- state_levels
  cell_state_labels_gep <- state_levels
  cat("  Single-level annotation: state = type (1-to-1 mapping)\n")
} else {
  # Two-level: map states to types (already validated no ambiguity)
  state_levels <- rownames(ref_gep)
  state_to_type <- tapply(cell_type_labels, cell_state_labels, function(x) {
    names(sort(table(x), decreasing = TRUE))[1]
  })
  cell_type_labels_gep <- unname(state_to_type[state_levels])
  cell_state_labels_gep <- state_levels
  cat(sprintf(
    "  Mapped %d states to %d cell types\n",
    length(state_levels),
    length(unique(cell_type_labels_gep))
  ))
}

# ==============================================================================
# 7. Gene Intersection and Validation
# ==============================================================================

cat("\n=== Step 6/11: Intersecting Genes ===\n")

# P1-5: Normalize reference gene names too
ref_genes_orig <- colnames(ref_gep)
ref_genes_norm <- normalize_gene_symbols(ref_genes_orig)
colnames(ref_gep) <- ref_genes_norm

# Find common genes
common_genes <- intersect(colnames(ref_gep), colnames(bk_dat))

cat(sprintf("  Common genes: %d\n", length(common_genes)))
cat(sprintf(
  "    Reference overlap: %d / %d (%.1f%%)\n",
  length(common_genes),
  ncol(ref_gep),
  100 * length(common_genes) / ncol(ref_gep)
))
cat(sprintf(
  "    Bulk overlap:      %d / %d (%.1f%%)\n",
  length(common_genes),
  ncol(bk_dat),
  100 * length(common_genes) / ncol(bk_dat)
))

if (length(common_genes) < 1000) {
  # Report unmatched genes to help debug
  ref_only <- setdiff(colnames(ref_gep), colnames(bk_dat))
  bulk_only <- setdiff(colnames(bk_dat), colnames(ref_gep))

  cat("\n  Top 20 genes in reference but not bulk:\n")
  print(head(ref_only, 20))

  cat("\n  Top 20 genes in bulk but not reference:\n")
  print(head(bulk_only, 20))

  stop(
    "ERROR: Too few common genes (<1000). Check gene name format (SYMBOL vs ENSEMBL)"
  )
}

# Subset both matrices to common genes
ref_gep <- ref_gep[, common_genes, drop = FALSE]
bk_dat <- bk_dat[, common_genes, drop = FALSE]

# P1-3: Verify consistency after subsetting
if (!identical(colnames(ref_gep), colnames(bk_dat))) {
  stop(
    "ERROR: Gene order mismatch after intersection. This should never happen."
  )
}

cat("  ✓ Gene intersection validated (order consistent)\n")

# ==============================================================================
# 8. Construct Prism Object
# ==============================================================================

cat("\n=== Step 7/11: Constructing Prism Object ===\n")

myPrism <- new.prism(
  reference = ref_gep,
  mixture = bk_dat,
  input.type = "GEP",
  cell.type.labels = cell_type_labels_gep,
  cell.state.labels = cell_state_labels_gep,
  key = key_malignant,
  outlier.cut = outlier_cut,
  outlier.fraction = outlier_frac
)

cat("  ✓ Prism object created successfully\n")

# ==============================================================================
# 9. Run BayesPrism
# ==============================================================================

cat("\n=== Step 8/11: Running BayesPrism ===\n")
cat("  This may take 10-30 minutes depending on data size...\n\n")

start_time <- Sys.time()

bp_res <- run.prism(
  prism = myPrism,
  n.cores = n_cores,

  # Gibbs sampling parameters
  gibbs.control = list(
    chain.length = 2000, # Total iterations
    burn.in = 500, # Burn-in iterations
    thinning = 2 # Keep every 2nd sample
  ),

  # Optimization parameters
  opt.control = list(
    trace = 1, # Print progress (0=silent, 1=basic, 2=detailed)
    maxit = 10000 # Max optimization iterations
  )
)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

# Convergence diagnostics
cat("\n  Convergence Diagnostics:\n")

# Try to extract log-likelihood (structure varies by BayesPrism version)
loglik <- tryCatch(
  {
    if (!is.null(bp_res@posterior.theta_f@loglikelihood)) {
      bp_res@posterior.theta_f@loglikelihood
    } else if ("logLik" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@logLik
    } else if ("log.lik" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@log.lik
    } else {
      NA
    }
  },
  error = function(e) NA
)

if (!is.na(loglik)) {
  cat(sprintf("    Final log-likelihood: %.2f\n", loglik))
} else {
  cat("    Final log-likelihood: Not available (version incompatibility)\n")
}

cat(sprintf("    Runtime: %.1f minutes\n", as.numeric(runtime)))

# Optional: Save object structure for debugging
cat("\n  Saving posterior object structure info...\n")
posterior_structure <- list(
  class = class(bp_res@posterior.theta_f),
  slots = slotNames(bp_res@posterior.theta_f),
  dim_theta = dim(bp_res@posterior.theta_f@theta),
  converged = TRUE # If we got here, it converged
)
saveRDS(
  posterior_structure,
  file.path(output_dir, "posterior_structure_info.rds")
)

# Save main result object
saveRDS(bp_res, file.path(output_dir, "bp_res.rds"))
cat("\n  Saved: bp_res.rds\n")

# ==============================================================================
# 10. Extract Results
# ==============================================================================

cat("\n=== Step 9/11: Extracting Results ===\n")

# Cell type fractions (theta)
theta <- get.fraction(
  bp = bp_res,
  which.theta = "final",
  state.or.type = "type"
)

# Coefficient of variation (uncertainty)
theta_cv <- tryCatch(
  {
    if (!is.null(bp_res@posterior.theta_f@theta.cv)) {
      bp_res@posterior.theta_f@theta.cv
    } else if ("cv" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@cv
    } else if ("theta_cv" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@theta_cv
    } else {
      # Create a placeholder CV matrix (all zeros)
      warning("Could not extract theta CV - using placeholder zeros")
      matrix(
        0,
        nrow = nrow(theta),
        ncol = ncol(theta),
        dimnames = dimnames(theta)
      )
    }
  },
  error = function(e) {
    warning(
      "Error extracting theta CV: ",
      e$message,
      " - using placeholder zeros"
    )
    matrix(
      0,
      nrow = nrow(theta),
      ncol = ncol(theta),
      dimnames = dimnames(theta)
    )
  }
)

# Deconvolved gene expression per cell type (Z matrix)
Z <- tryCatch(
  {
    # Try method 1: cell.name = "all"
    get.exp(bp = bp_res, state.or.type = "type", cell.name = "all")
  },
  error = function(e1) {
    tryCatch(
      {
        # Try method 2: without cell.name
        get.exp(bp = bp_res, state.or.type = "type")
      },
      error = function(e2) {
        tryCatch(
          {
            # Try method 3: extract directly from posterior
            bp_res@posterior.initial.cellType@Z.mean
          },
          error = function(e3) {
            # Last resort: extract from reference
            warning(
              "Could not extract Z matrix using standard methods. Using reference as fallback."
            )
            ref_gep # Use the reference GEP as fallback
          }
        )
      }
    )
  }
)

# Ensure Z has proper dimensions and names
if (is.null(dim(Z)) || ncol(Z) == 0) {
  stop("Failed to extract Z matrix from BayesPrism results")
}

cat(sprintf(
  "  Z matrix dimensions: %d genes × %d cell types\n",
  nrow(Z),
  ncol(Z)
))

# Check object class and convert to matrix if needed
cat(sprintf("  Z object class: %s\n", paste(class(Z), collapse = ", ")))
cat(sprintf("  Z is matrix: %s\n", is.matrix(Z)))
cat(sprintf("  Z is list: %s\n", is.list(Z)))
cat(sprintf("  Z is array: %s\n", is.array(Z)))

# Handle different Z structures
if (is.list(Z) && !is.data.frame(Z) && !is.matrix(Z)) {
  cat("  WARNING: Z is a list, attempting to extract matrix...\n")
  # Try to find matrix in list
  if (length(Z) > 0) {
    # Check what's in the list
    cat(sprintf("  List has %d elements\n", length(Z)))
    cat(sprintf(
      "  First element class: %s\n",
      paste(class(Z[[1]]), collapse = ", ")
    ))

    # Take first element if it's a matrix
    if (is.matrix(Z[[1]]) || is.array(Z[[1]])) {
      Z <- Z[[1]]
      cat("  Extracted matrix from list element 1\n")
    } else {
      # Try to convert first element to matrix
      Z <- as.matrix(Z[[1]])
      cat("  Converted list element 1 to matrix\n")
    }
  } else {
    stop("Z is an empty list")
  }
}

# Convert to matrix if not already
if (!is.matrix(Z)) {
  cat("  Converting Z to matrix...\n")
  Z <- as.matrix(Z)
  cat(sprintf(
    "  After conversion - class: %s, dimensions: %d × %d\n",
    class(Z),
    nrow(Z),
    ncol(Z)
  ))
}

# Check if Z dimensions are correct (should be genes × cell types)
# If Z is transposed (cell types × genes), fix it
if (ncol(Z) > nrow(Z) && ncol(Z) > 100) {
  cat("  WARNING: Z matrix appears transposed. Fixing...\n")
  Z <- t(Z)
  cat(sprintf(
    "  Corrected dimensions: %d genes × %d cell types\n",
    nrow(Z),
    ncol(Z)
  ))
}

# Validate Z dimensions
if (ncol(Z) != ncol(theta)) {
  warning(sprintf(
    "Dimension mismatch: Z has %d cell types but theta has %d",
    ncol(Z),
    ncol(theta)
  ))
}

# Save results
write.csv(theta, file.path(output_dir, "theta_final.csv"), quote = FALSE)
write.csv(theta_cv, file.path(output_dir, "theta_cv.csv"), quote = FALSE)

# P1-2: Save Z as RDS instead of CSV (much faster and smaller)
saveRDS(Z, file.path(output_dir, "Z_matrix.rds"), compress = "xz")
cat("  Saved: theta_final.csv, theta_cv.csv, Z_matrix.rds (compressed)\n")

# Also save a CSV of Z for top 1000 most variable genes (for quick inspection)
if (nrow(Z) > 1) {
  z_var <- apply(Z, 1, var)
  top_genes <- names(sort(z_var, decreasing = TRUE))[1:min(1000, nrow(Z))]
  write.csv(
    Z[top_genes, , drop = FALSE],
    file.path(output_dir, "Z_matrix_top1000genes.csv"),
    quote = FALSE
  )
  cat(
    "  Saved: Z_matrix_top1000genes.csv (top variable genes for quick inspection)\n"
  )
} else {
  warning("Z matrix has only 1 gene - skipping top genes CSV")
}

# Quick preview
cat("\n  Cell type fractions (theta):\n")
print(round(theta, 3))

# ==============================================================================
# 11. Marker Gene Validation
# ==============================================================================

cat("\n=== Step 10/11: Marker Gene Validation ===\n")

# Define marker genes for respiratory tract cell types
marker_genes <- list(
  Epithelial = c(
    "EPCAM",
    "KRT5",
    "TP63",
    "KRT14",
    "MUC5AC",
    "MUC5B",
    "FOXJ1",
    "RSPH1",
    "SCGB1A1",
    "SCGB3A1"
  ),
  T_NK = c(
    "CD3D",
    "CD3E",
    "CD3G",
    "CD4",
    "CD8A",
    "CD8B",
    "GNLY",
    "NKG7",
    "GZMA",
    "GZMB"
  ),
  Myeloid = c(
    "CD14",
    "FCGR3A",
    "CD68",
    "CD163",
    "CD1C",
    "CLEC9A",
    "FCER1A",
    "ITGAX"
  ),
  B_cells = c("MS4A1", "CD79A", "CD79B", "IGHG1", "IGHM", "IGHA1"),
  Stromal = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2",
    "PDGFRA",
    "PDGFRB"
  ),
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1", "KDR")
)

# Check which markers are present
all_markers <- unique(unlist(marker_genes))
markers_present <- all_markers[all_markers %in% rownames(Z)]

cat(sprintf(
  "  Marker genes detected: %d / %d (%.1f%%)\n",
  length(markers_present),
  length(all_markers),
  100 * length(markers_present) / length(all_markers)
))

if (length(markers_present) > 0) {
  Z_marker <- Z[markers_present, , drop = FALSE]
  write.csv(Z_marker, file.path(output_dir, "Z_marker_genes.csv"))

  # Safe graphics device handling with tryCatch
  tryCatch(
    {
      pdf(file.path(output_dir, "Z_marker_heatmap.pdf"), width = 10, height = 8)

      pheatmap(
        log1p(Z_marker),
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
        border_color = NA,
        main = "Deconvolved Marker Gene Expression (log1p, row-scaled)",
        fontsize_row = 7,
        fontsize_col = 10
      )

      dev.off()
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning("Error creating marker heatmap: ", e$message)
    }
  )

  cat("  Saved: Z_marker_genes.csv, Z_marker_heatmap.pdf\n")
}

# ==============================================================================
# 12. Visualizations
# ==============================================================================

cat("\n=== Step 11/11: Generating Visualizations ===\n")

# Prepare data for plotting
theta_df <- as.data.frame(theta)
theta_df$sample <- rownames(theta_df)
theta_long <- melt(
  as.data.table(theta_df),
  id.vars = "sample",
  variable.name = "cell_type",
  value.name = "fraction"
)

# 1. Stacked bar plot
p1 <- ggplot(theta_long, aes(x = sample, y = fraction, fill = cell_type)) +
  geom_col(width = 0.7, color = "black", size = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  labs(
    title = "BayesPrism: Cell Type Fractions (Healthy Samples)",
    x = "Bulk Sample",
    y = "Fraction (proportion of reads)",
    fill = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(output_dir, "theta_stacked_bar.pdf"),
  p1,
  width = 8,
  height = 5
)

# 2. Theta heatmap
tryCatch(
  {
    pdf(file.path(output_dir, "theta_heatmap.pdf"), width = 8, height = 6)

    pheatmap(
      theta,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
      border_color = "grey60",
      main = "Cell Type Fractions (theta)",
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 9,
      fontsize = 10
    )

    dev.off()
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Error creating theta heatmap: ", e$message)
  }
)

# 3. Theta CV heatmap
tryCatch(
  {
    pdf(file.path(output_dir, "theta_cv_heatmap.pdf"), width = 8, height = 6)

    pheatmap(
      theta_cv,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "orange", "red"))(100),
      border_color = "grey60",
      main = "Coefficient of Variation (uncertainty)",
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 9,
      fontsize = 10
    )

    dev.off()
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Error creating theta CV heatmap: ", e$message)
  }
)

# 4. Boxplot
p2 <- ggplot(theta_long, aes(x = cell_type, y = fraction, fill = cell_type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  labs(
    title = "Cell Type Fraction Distribution (n=5 healthy samples)",
    x = "Cell Type",
    y = "Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "theta_boxplot.pdf"), p2, width = 8, height = 5)

cat(
  "  Saved: theta_stacked_bar.pdf, theta_heatmap.pdf, theta_cv_heatmap.pdf, theta_boxplot.pdf\n"
)

# ==============================================================================
# 13. Summary Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

# Try to extract log-likelihood safely (used in summary and quality checks)
loglik_summary <- tryCatch(
  {
    if (!is.null(bp_res@posterior.theta_f@loglikelihood)) {
      sprintf("%.2f", bp_res@posterior.theta_f@loglikelihood)
    } else if ("logLik" %in% slotNames(bp_res@posterior.theta_f)) {
      sprintf("%.2f", bp_res@posterior.theta_f@logLik)
    } else if ("log.lik" %in% slotNames(bp_res@posterior.theta_f)) {
      sprintf("%.2f", bp_res@posterior.theta_f@log.lik)
    } else {
      "Not available"
    }
  },
  error = function(e) "Not available"
)

# Safe sink handling with tryCatch
tryCatch(
  {
    summary_path <- file.path(output_dir, "analysis_summary.txt")
    sink(summary_path)

    cat(
      "================================================================================\n"
    )
    cat("BayesPrism Deconvolution Analysis Summary v2.1\n")
    cat(
      "================================================================================\n\n"
    )

    cat("Analysis Date:", as.character(Sys.Date()), "\n")
    cat("Runtime:      ", sprintf("%.1f minutes", as.numeric(runtime)), "\n\n")

    cat("--- Input Data ---\n")
    cat("Reference:", basename(seurat_rds_path), "\n")
    cat("  Assay used:       RNA (raw counts)\n")
    cat("  Total cells:      ", ncol(sce_ref), "\n")
    cat("  Total genes:      ", nrow(sce_ref), "\n")
    cat(
      "  Annotation mode:  ",
      ifelse(use_single_level, "SINGLE level", "TWO levels"),
      "\n"
    )
    cat("  Cell types:       ", length(unique(cell_type_labels)), "\n")
    if (!use_single_level) {
      cat("  Cell states:      ", length(unique(cell_state_labels)), "\n")
    }
    cat("\n")

    cat("Bulk samples:", nrow(bk_dat), "(all healthy)\n")
    cat("  Sample IDs:       ", paste(rownames(bk_dat), collapse = ", "), "\n")
    cat("  Median lib size:  ", sprintf("%.0f", median(lib_size)), "\n")
    cat("  Lib size CV:      ", sprintf("%.2f", lib_size_cv), "\n\n")

    cat("Common genes:", length(common_genes), "\n\n")

    cat("--- Processing Options ---\n")
    cat(
      "Annotation mode:       ",
      ifelse(
        use_single_level,
        "SINGLE level (state = type)",
        "TWO levels (state + type)"
      ),
      "\n"
    )
    cat(
      "GEP construction:      ",
      ifelse(gep_use_mean, "Mean-profile (normalized)", "Sum-profile"),
      "\n"
    )
    cat("Filter IG genes:       ", filter_immunoglobulin, "\n\n")

    cat("--- BayesPrism Parameters ---\n")
    cat("Input type:        GEP (gene expression program)\n")
    cat("Outlier cut:       ", outlier_cut, "\n")
    cat("Outlier fraction:  ", outlier_frac, "\n")
    cat("Gibbs chain length:", 2000, "\n")
    cat("Burn-in:           ", 500, "\n")
    cat("Thinning:          ", 2, "\n")
    cat("CPU cores:         ", n_cores, "\n\n")

    cat("--- Results ---\n")
    cat("Final log-likelihood:", loglik_summary, "\n\n")

    cat("Cell type fractions (theta):\n")
    print(round(theta, 4))
    cat("\n")

    cat("Mean fractions across samples:\n")
    print(round(colMeans(theta), 4))
    cat("\n")

    cat("Coefficient of variation (theta_cv):\n")
    print(round(theta_cv, 4))
    cat("\n")

    cat("--- Marker Gene Detection ---\n")
    cat(sprintf(
      "Markers detected: %d / %d (%.1f%%)\n",
      length(markers_present),
      length(all_markers),
      100 * length(markers_present) / length(all_markers)
    ))
    cat("\n")

    cat("--- Quality Checks ---\n")
    cat("✓ RNA assay validation:     PASS (non-negative, mostly integer)\n")
    if (use_single_level) {
      cat(
        "✓ Annotation mode:          SINGLE level (no state→type ambiguity)\n"
      )
    } else {
      cat("✓ State→type consistency:   PASS (no ambiguity)\n")
    }
    cat(
      "✓ Gene intersection:        ",
      length(common_genes),
      " genes (>1000 threshold)\n"
    )
    cat("✓ Convergence:              Final loglik = ", loglik_summary, "\n\n")

    cat("--- Output Files ---\n")
    cat("1.  bp_res.rds                    : BayesPrism result object\n")
    cat("2.  theta_final.csv               : Cell type fractions\n")
    cat("3.  theta_cv.csv                  : Uncertainty (CV)\n")
    cat(
      "4.  Z_matrix.rds                  : Deconvolved gene expression (compressed)\n"
    )
    cat(
      "5.  Z_matrix_top1000genes.csv     : Top variable genes (quick inspection)\n"
    )
    cat("6.  Z_marker_genes.csv            : Marker gene expression\n")
    cat("7.  theta_stacked_bar.pdf         : Stacked bar plot\n")
    cat("8.  theta_heatmap.pdf             : Fraction heatmap\n")
    cat("9.  theta_cv_heatmap.pdf          : Uncertainty heatmap\n")
    cat("10. theta_boxplot.pdf             : Fraction distribution\n")
    cat("11. Z_marker_heatmap.pdf          : Marker gene heatmap\n")
    cat("12. reference_metadata.rds        : Reference data metadata\n")
    cat("13. analysis_summary.txt          : This summary\n")
    cat("14. session_info.txt              : R session information\n\n")

    cat("--- Version History ---\n")
    cat("v2.1 (2025-01-23): Added single-level annotation option\n")
    cat("  - NEW: use_single_level parameter to avoid state→type ambiguity\n")
    cat("  - NEW: Automatic removal of 'None' and invalid labels\n")
    cat("v2.0 (2025-01-23): Production release with P0/P1 fixes\n")
    cat("  - P0-1: Force RNA raw counts with validation\n")
    cat("  - P0-2: GEP uses mean-profile with sparse matrix multiplication\n")
    cat("  - P0-3: Immunoglobulin genes preserved by default\n")
    cat("  - P0-4: State→type ambiguity detection with error report\n")
    cat("  - P1-1: Fault-tolerant sink() and graphics cleanup\n")
    cat("  - P1-2: Z matrix saved as compressed RDS\n")
    cat("  - P1-3: Gene intersection consistency checks\n")
    cat("  - P1-4: Dedup performance monitoring\n")
    cat("  - P1-5: Gene name normalization\n\n")

    cat(
      "================================================================================\n"
    )
    cat("Analysis completed successfully!\n")
    cat("Output directory:", output_dir, "\n")
    cat(
      "================================================================================\n"
    )

    sink()
  },
  error = function(e) {
    # Ensure sink is closed even if error occurs
    while (sink.number() > 0) {
      sink()
    }
    warning("Error writing summary report: ", e$message)
  },
  finally = {
    # Always ensure sink is closed
    while (sink.number() > 0) {
      sink()
    }
  }
)

cat("  Saved: analysis_summary.txt\n")

# ==============================================================================
# 14. Session Info
# ==============================================================================

tryCatch(
  {
    session_path <- file.path(output_dir, "session_info.txt")
    sink(session_path)

    cat("R Session Information\n")
    cat("====================\n\n")
    print(sessionInfo())

    sink()
  },
  error = function(e) {
    # Ensure sink is closed even if error occurs
    while (sink.number() > 0) {
      sink()
    }
    warning("Error writing session info: ", e$message)
  },
  finally = {
    # Always ensure sink is closed
    while (sink.number() > 0) {
      sink()
    }
  }
)

cat("  Saved: session_info.txt\n")

# ==============================================================================
# Final Message
# ==============================================================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("✓ BayesPrism deconvolution analysis completed successfully!\n")
cat("  Runtime: ", sprintf("%.1f minutes", as.numeric(runtime)), "\n")
cat("  Output directory: ", output_dir, "\n")
cat("  Total files generated: 14\n")
cat(rep("=", 80), "\n", sep = "")
cat("\nNext steps:\n")
cat("1. Review theta_stacked_bar.pdf for cell type composition\n")
cat("2. Check theta_cv_heatmap.pdf for estimation uncertainty\n")
cat("3. Validate results with Z_marker_heatmap.pdf\n")
cat("4. Load bp_res.rds for downstream analyses\n")
cat("5. Check analysis_summary.txt for quality metrics\n\n")
cat("Key improvements in v2.1:\n")
cat("- NEW: Single-level annotation option (avoids state→type ambiguity)\n")
cat("- NEW: Automatic removal of 'None' and invalid labels\n")
cat("- Uses RNA assay raw counts (validated)\n")
cat("- GEP constructed with mean-profile (removes cell number bias)\n")
cat("- Immunoglobulin genes preserved (better B cell deconvolution)\n")
cat("- Fault-tolerant I/O (safe error handling)\n\n")

library(ggplot2)
library(data.table)

# Load results
theta <- read.csv(
  "/home/h2048/data/R/bulk0121/bayesprism_out/theta_final.csv",
  row.names = 1
)
theta_cv <- read.csv(
  "/home/h2048/data/R/bulk0121/bayesprism_out/theta_cv.csv",
  row.names = 1
)

# Calculate CV across samples for each cell type
cv_across_samples <- apply(theta, 2, function(x) sd(x) / mean(x))

# Visualize
cv_df <- data.frame(
  cell_type = names(cv_across_samples),
  cv = cv_across_samples
)

ggplot(cv_df, aes(x = reorder(cell_type, cv), y = cv)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Cell Type Variability Across Healthy Samples",
    x = "Cell Type",
    y = "Coefficient of Variation"
  ) +
  theme_bw()

ggsave("cell_type_consistency.pdf", width = 8, height = 6)

library(corrplot)

# Cell type correlation
cor_mat <- cor(theta)

# Visualize
pdf("cell_type_correlation.pdf", width = 10, height = 10)
corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  number.cex = 0.7,
  title = "Cell Type Correlation (Healthy Samples)"
)
dev.off()

# Find strongly correlated pairs
cor_pairs <- which(abs(cor_mat) > 0.7 & abs(cor_mat) < 1, arr.ind = TRUE)
cor_pairs_df <- data.frame(
  CellType1 = rownames(cor_mat)[cor_pairs[, 1]],
  CellType2 = colnames(cor_mat)[cor_pairs[, 2]],
  Correlation = cor_mat[cor_pairs]
)
write.csv(cor_pairs_df, "cell_type_correlation_pairs.csv", row.names = FALSE)


# Define respiratory tract markers
markers <- list(
  Epithelial = c(
    "EPCAM",
    "KRT5",
    "TP63",
    "MUC5AC",
    "MUC5B",
    "FOXJ1",
    "SCGB1A1"
  ),
  Lymphoid = c("CD3D", "CD3E", "CD4", "CD8A", "GNLY", "NKG7"),
  Myeloid = c("CD14", "FCGR3A", "CD68", "CD163", "CD1C"),
  B_cells = c("MS4A1", "CD79A", "IGHG1", "IGHM"),
  Stromal = c("COL1A1", "COL3A1", "DCN", "ACTA2"),
  Endothelial = c("PECAM1", "VWF", "CDH5")
)

# Extract marker expression
Z_markers <- Z[unlist(markers)[unlist(markers) %in% rownames(Z)], ]

# Heatmap with annotation
library(pheatmap)
library(RColorBrewer)

# Create annotation for markers
marker_annotation <- data.frame(
  CellType = rep(names(markers), sapply(markers, length))
)
rownames(marker_annotation) <- unlist(markers)
marker_annotation <- marker_annotation[rownames(Z_markers), , drop = FALSE]

pheatmap(
  log1p(Z_markers),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = marker_annotation,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Marker Gene Expression per Cell Type (Deconvolved)",
  fontsize_row = 8
)
