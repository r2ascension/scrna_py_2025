#!/usr/bin/env Rscript

# ==============================================================================
# BayesPrism Deconvolution Pipeline
# ==============================================================================
# Purpose: Deconvolve bulk RNA-seq using single-cell reference
# Author: r2end
# Date: 2025-01-26
# Version: 3.0
#
# Input:
#   - Seurat RDS: single-cell reference with ann_level_2 / ann_level_3
#   - Bulk CSV: gene_count.csv (gene_name, P1, P2, ...)
#
# Output:
#   - Cell type fractions (theta) + uncertainty (CV)
#   - Deconvolved gene expression (Z matrix)
#   - QC reports and visualizations
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
  library(corrplot)
})

# ==============================================================================
# 0. Configuration
# ==============================================================================

# ===== Paths (MODIFY THESE) =====
seurat_rds_path <- "/home/h2048/data/R/bulk0121/LungMap.rds"
bulk_csv_path <- "/home/h2048/data/bulk/gene_count.csv"
output_dir <- "/home/h2048/data/R/bulk0303/bayesprism_out"

# ===== Cell Type Annotation Strategy =====
# TRUE  = use single level (celltype_col only, no state subdivision, recommended)
# FALSE = use two levels (celltype_col as type, cellstate_col as state)
use_single_level <- TRUE
celltype_col <- "ann_level_3_clean"
# cellstate_col  <- "ann_level_2"   # only needed when use_single_level = FALSE

if (use_single_level) {
  cellstate_col <- celltype_col
} else {
  if (!exists("cellstate_col")) {
    stop("cellstate_col not defined for two-level mode")
  }
}

# ===== BayesPrism Parameters =====
key_malignant <- NULL # NULL for healthy / non-tumor tissue
outlier_cut <- 0.01
outlier_frac <- 0.1
filter_immunoglobulin <- FALSE # Keep IG genes for immune deconvolution
gep_use_mean <- TRUE # TRUE = mean-profile (removes cell number bias)
n_cores <- 32

# ===== Reproducibility =====
set.seed(42)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat(
  "================================================================================\n"
)
cat("BayesPrism Deconvolution Pipeline v3.0\n")
cat(
  "================================================================================\n"
)
cat("Seurat RDS:      ", seurat_rds_path, "\n")
cat("Bulk CSV:        ", bulk_csv_path, "\n")
cat("Output dir:      ", output_dir, "\n")
cat(
  "Annotation mode: ",
  ifelse(use_single_level, "SINGLE level", "TWO levels"),
  "\n"
)
cat("Cell type col:   ", celltype_col, "\n")
cat("Filter IG genes: ", filter_immunoglobulin, "\n")
cat(
  "GEP method:      ",
  ifelse(gep_use_mean, "mean-profile", "sum-profile"),
  "\n"
)
cat("CPU cores:       ", n_cores, "\n")
cat(
  "================================================================================\n\n"
)

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

# Normalize gene symbols: strip version suffix, convert to uppercase
normalize_gene_symbols <- function(gene_vec) {
  toupper(sub("\\.\\d+$", "", gene_vec))
}

# Return logical keep-vector: remove ribosomal, MT, HB, HSP, lncRNA
# Optionally remove IG genes (not recommended for immune deconvolution)
filter_technical_genes <- function(gene_vec, filter_ig = FALSE) {
  bad <- grepl("^MT-", gene_vec) |
    grepl("^RPS", gene_vec) |
    grepl("^RPL", gene_vec) |
    grepl("^MRPS", gene_vec) |
    grepl("^MRPL", gene_vec) |
    grepl("^HBA", gene_vec) |
    grepl("^HBB", gene_vec) |
    grepl("^HBG", gene_vec) |
    grepl("^HBD", gene_vec) |
    grepl("^HBE", gene_vec) |
    grepl("^HSPA", gene_vec) |
    grepl("^HSPB", gene_vec) |
    grepl("^DNAJ", gene_vec) |
    gene_vec %in% c("MALAT1", "NEAT1", "XIST")

  if (filter_ig) {
    bad <- bad |
      grepl("^IGH[GAMED]", gene_vec) |
      grepl("^IGK[CV]", gene_vec) |
      grepl("^IGL[CV]", gene_vec)
    message("  NOTE: IG genes filtered (may reduce B cell accuracy)")
  } else {
    message("  NOTE: IG genes preserved (recommended for immune deconvolution)")
  }
  !bad
}

# Deduplicate gene rows in sparse matrix by summing counts
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

  t0 <- Sys.time()
  idx_list <- split(seq_along(g), g)
  out <- Matrix(0, nrow = length(idx_list), ncol = ncol(mat_gc), sparse = TRUE)
  rownames(out) <- names(idx_list)
  colnames(out) <- colnames(mat_gc)

  for (i in seq_along(idx_list)) {
    ridx <- idx_list[[i]]
    out[i, ] <- if (length(ridx) == 1L) {
      mat_gc[ridx, ]
    } else {
      Matrix::colSums(mat_gc[ridx, , drop = FALSE])
    }
  }

  message(sprintf(
    "  Deduplication done in %.1f s",
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
  out
}

# Collapse single-cell counts to GEP: genes x cells -> states x genes
collapse_to_gep <- function(
  counts_gc,
  state_vec,
  use_mean = TRUE,
  filter_ig = FALSE
) {
  state_vec <- factor(state_vec)
  g <- rownames(counts_gc)
  keep_gene <- filter_technical_genes(g, filter_ig = filter_ig)
  counts_filtered <- counts_gc[keep_gene, , drop = FALSE]
  message(sprintf(
    "  Filtered %d technical genes (%d -> %d)",
    sum(!keep_gene),
    length(g),
    sum(keep_gene)
  ))

  S <- Matrix::sparse.model.matrix(~ 0 + state_vec)
  colnames(S) <- sub("^state_vec", "", colnames(S))
  gep_gs <- counts_filtered %*% S

  if (use_mean) {
    n_per_state <- Matrix::colSums(S)
    gep_gs <- gep_gs %*% Matrix::Diagonal(x = 1 / pmax(n_per_state, 1))
    message("  GEP: mean-profile")
  } else {
    message("  GEP: sum-profile")
  }

  gep <- as.matrix(t(gep_gs))
  rownames(gep) <- colnames(S)
  colnames(gep) <- rownames(counts_filtered)

  keep_states <- rowSums(gep) > 0
  if (!all(keep_states)) {
    warning(sprintf("  Removing %d zero-count states", sum(!keep_states)))
    gep <- gep[keep_states, , drop = FALSE]
  }
  gep
}

# Safe extraction of Z matrix from BayesPrism result object
extract_z_matrix <- function(bp_res) {
  Z_arr <- get.exp(bp = bp_res, state.or.type = "type")
  # Z_arr: samples x genes x cell_types -> average across samples -> genes x cell_types
  Z <- apply(Z_arr, c(2, 3), mean)
  cat(sprintf(
    "  Z array [%s] -> averaged to [%d genes x %d cell types]\n",
    paste(dim(Z_arr), collapse = " x "),
    nrow(Z),
    ncol(Z)
  ))
  Z
}

# ==============================================================================
# 2. Load Seurat Reference
# ==============================================================================

cat("=== Step 1/9: Loading Seurat Reference ===\n")

if (!file.exists(seurat_rds_path)) {
  stop("Seurat RDS not found: ", seurat_rds_path)
}

seurat_obj <- readRDS(seurat_rds_path)
cat(sprintf(
  "  Loaded: %d cells x %d features\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
))

# ===== 添加新的细胞类型标注列 =====
finest_col <- "ann_finest_level" # <- 修改为最细标注列名
celltype_col <- "ann_level_3_clean" # <- 新列名，后续流程使用此列
cellstate_col <- "ann_level_3_clean"
cd <- seurat_obj@meta.data
cd[["ann_level_3_clean"]] <- as.character(cd[["ann_level_3"]])

# None / Rare 替换为最细标注
replace_idx <- cd[["ann_level_3_clean"]] %in% c("None", "Rare")
cd[replace_idx, "ann_level_3_clean"] <- as.character(cd[
  replace_idx,
  "ann_finest_level"
])

# Lymphatic EC 开头统一归为 LEC
cd[["ann_level_3_clean"]] <- ifelse(
  grepl("^Lymphatic EC", cd[["ann_level_3_clean"]]),
  "LEC",
  cd[["ann_level_3_clean"]]
)

# 合并平滑肌三个亚型为 SMC
cd[["ann_level_3_clean"]] <- ifelse(
  cd[["ann_level_3_clean"]] %in%
    c("SM activated stress response", "Smooth muscle", "Smooth muscle FAM83D+"),
  "SMC",
  cd[["ann_level_3_clean"]]
)

seurat_obj@meta.data <- cd
celltype_col <- "ann_level_3_clean"
cellstate_col <- "ann_level_3_clean"

table(seurat_obj$ann_level_3_clean)

available_assays <- Seurat::Assays(seurat_obj)
cat("  Available assays:", paste(available_assays, collapse = ", "), "\n")
if (!"RNA" %in% available_assays) {
  stop(
    "No RNA assay found. Available: ",
    paste(available_assays, collapse = ", ")
  )
}

Seurat::DefaultAssay(seurat_obj) <- "RNA"
sce_ref <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
counts_gc <- counts(sce_ref)
if (!inherits(counts_gc, "dgCMatrix")) {
  counts_gc <- as(counts_gc, "dgCMatrix")
}

if (any(counts_gc@x < 0)) {
  stop("Reference counts contain negative values")
}
non_int_frac <- mean(abs(counts_gc@x - round(counts_gc@x)) > 1e-6)
if (non_int_frac > 0.01) {
  warning(sprintf(
    "%.1f%% non-integer counts - confirm raw UMI",
    100 * non_int_frac
  ))
}
cat("  Counts validation: PASS\n")

missing_cols <- setdiff(
  unique(c(celltype_col, cellstate_col)),
  colnames(colData(sce_ref))
)
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

ref_metadata <- list(
  source = basename(seurat_rds_path),
  n_cells = ncol(sce_ref),
  n_genes = nrow(sce_ref),
  use_single_level = use_single_level,
  celltype_col = celltype_col,
  cellstate_col = cellstate_col,
  celltype_counts = table(colData(sce_ref)[[celltype_col]]),
  creation_date = Sys.Date(),
  gep_method = ifelse(gep_use_mean, "mean-profile", "sum-profile"),
  filter_ig = filter_immunoglobulin
)
saveRDS(ref_metadata, file.path(output_dir, "reference_metadata.rds"))

rm(seurat_obj)
invisible(gc())

# ==============================================================================
# 3. Reference QC
# ==============================================================================

cat("\n=== Step 2/9: Reference QC ===\n")

cd <- as.data.frame(colData(sce_ref))
cell_type_labels <- as.character(cd[[celltype_col]])
cell_state_labels <- as.character(cd[[cellstate_col]])

# Remove NA labels
ok <- !is.na(cell_type_labels) & !is.na(cell_state_labels)
if (sum(!ok) > 0) {
  cat(sprintf("  Removing %d NA-label cells\n", sum(!ok)))
  sce_ref <- sce_ref[, ok]
  counts_gc <- counts_gc[, ok]
  cell_type_labels <- cell_type_labels[ok]
  cell_state_labels <- cell_state_labels[ok]
}

# Remove invalid labels
bad_labels <- cell_type_labels %in%
  c("None", "none", "", "NA", "Unknown", "unknown") |
  cell_state_labels %in% c("None", "none", "", "NA", "Unknown", "unknown")
if (sum(bad_labels) > 0) {
  cat(sprintf("  Removing %d invalid-label cells\n", sum(bad_labels)))
  sce_ref <- sce_ref[, !bad_labels]
  counts_gc <- counts_gc[, !bad_labels]
  cell_type_labels <- cell_type_labels[!bad_labels]
  cell_state_labels <- cell_state_labels[!bad_labels]
}

celltype_counts <- table(cell_type_labels)
cat("\n  Cell type distribution:\n")
print(addmargins(celltype_counts))

small_types <- names(celltype_counts)[celltype_counts < 30]
if (length(small_types) > 0) {
  cat("  WARNING: <30 cells in:", paste(small_types, collapse = ", "), "\n")
}

# Two-level mode: validate state->type mapping has no ambiguity
if (!use_single_level) {
  state_type_table <- table(cell_state_labels, cell_type_labels)
  ambig <- rowSums(state_type_table > 0) > 1
  if (any(ambig)) {
    write.csv(
      state_type_table[ambig, , drop = FALSE],
      file.path(output_dir, "state_type_ambiguity.csv"),
      quote = FALSE
    )
    stop("State->type ambiguity detected. See state_type_ambiguity.csv")
  }
  cat("  State->type mapping: PASS\n")
} else {
  cat("  Annotation mode: SINGLE level (state = type)\n")
}

gene_detection <- Matrix::rowMeans(counts_gc > 0)
cat(sprintf(
  "  Genes detected in >10%% cells: %d / %d (%.1f%%)\n",
  sum(gene_detection > 0.1),
  length(gene_detection),
  100 * mean(gene_detection > 0.1)
))
# 移除与咽鼓管不相关的细胞类型
remove_types <- c("AT1", "AT2", "Hematopoietic stem cells")
keep_cells <- !cell_type_labels %in% remove_types
sce_ref <- sce_ref[, keep_cells]
counts_gc <- counts_gc[, keep_cells]
cell_type_labels <- cell_type_labels[keep_cells]
cell_state_labels <- cell_state_labels[keep_cells]
cat(sprintf(
  "  Removed %s: remaining %d cells\n",
  paste(remove_types, collapse = ", "),
  ncol(sce_ref)
))
counts_gc <- dedup_gene_rows_sparse(counts_gc)

# ==============================================================================
# 4. Load Bulk Data
# ==============================================================================

cat("\n=== Step 3/9: Loading Bulk Data ===\n")

if (!file.exists(bulk_csv_path)) {
  stop("Bulk CSV not found: ", bulk_csv_path)
}

bulk_raw <- fread(bulk_csv_path)
setnames(bulk_raw, 1, "gene_name")
bulk_raw[, gene_name := normalize_gene_symbols(as.character(gene_name))]

bulk_agg <- bulk_raw[,
  lapply(.SD, sum),
  by = gene_name,
  .SDcols = names(bulk_raw)[-1]
]
bulk_mat <- as.matrix(bulk_agg[, -1])
rownames(bulk_mat) <- bulk_agg$gene_name

keep_bulk <- filter_technical_genes(
  rownames(bulk_mat),
  filter_ig = filter_immunoglobulin
)
bulk_mat <- bulk_mat[keep_bulk, , drop = FALSE]

bk_dat <- t(bulk_mat)
mode(bk_dat) <- "numeric"

lib_size <- colSums(bulk_mat)
lib_size_cv <- sd(lib_size) / mean(lib_size)
cat(sprintf(
  "  Loaded: %d genes x %d samples\n",
  nrow(bulk_mat),
  ncol(bulk_mat)
))
cat("  Samples:", paste(colnames(bulk_mat), collapse = ", "), "\n")
cat(sprintf(
  "  Library size: median=%.0f  CV=%.2f\n",
  median(lib_size),
  lib_size_cv
))
if (lib_size_cv > 0.5) {
  cat("  WARNING: High library size variation (CV > 0.5)\n")
}

# ==============================================================================
# 5. Collapse Reference to GEP
# ==============================================================================

cat("\n=== Step 4/9: Collapsing Reference to GEP ===\n")

ref_gep <- collapse_to_gep(
  counts_gc,
  cell_state_labels,
  gep_use_mean,
  filter_immunoglobulin
)
cat(sprintf("  GEP: %d states x %d genes\n", nrow(ref_gep), ncol(ref_gep)))

if (use_single_level) {
  cell_type_labels_gep <- rownames(ref_gep)
  cell_state_labels_gep <- rownames(ref_gep)
} else {
  state_to_type <- tapply(cell_type_labels, cell_state_labels, function(x) {
    names(sort(table(x), decreasing = TRUE))[1]
  })
  cell_type_labels_gep <- unname(state_to_type[rownames(ref_gep)])
  cell_state_labels_gep <- rownames(ref_gep)
  cat(sprintf(
    "  Mapped %d states to %d cell types\n",
    length(rownames(ref_gep)),
    length(unique(cell_type_labels_gep))
  ))
}

# ==============================================================================
# 6. Gene Intersection
# ==============================================================================

cat("\n=== Step 5/9: Intersecting Genes ===\n")

colnames(ref_gep) <- normalize_gene_symbols(colnames(ref_gep))
common_genes <- intersect(colnames(ref_gep), colnames(bk_dat))
cat(sprintf(
  "  Common genes: %d  (ref %.1f%%  bulk %.1f%%)\n",
  length(common_genes),
  100 * length(common_genes) / ncol(ref_gep),
  100 * length(common_genes) / ncol(bk_dat)
))

if (length(common_genes) < 1000) {
  cat(
    "  Top 20 in reference only:",
    paste(
      head(setdiff(colnames(ref_gep), colnames(bk_dat)), 20),
      collapse = ", "
    ),
    "\n"
  )
  cat(
    "  Top 20 in bulk only:     ",
    paste(
      head(setdiff(colnames(bk_dat), colnames(ref_gep)), 20),
      collapse = ", "
    ),
    "\n"
  )
  stop(
    "Too few common genes (<1000). Check gene name format (SYMBOL vs ENSEMBL)"
  )
}

ref_gep <- ref_gep[, common_genes, drop = FALSE]
bk_dat <- bk_dat[, common_genes, drop = FALSE]
if (!identical(colnames(ref_gep), colnames(bk_dat))) {
  stop("Gene order mismatch after intersection")
}
cat("  Gene intersection: PASS\n")

# ==============================================================================
# 7. Construct Prism & Run BayesPrism
# ==============================================================================

cat("\n=== Step 6/9: Constructing Prism Object ===\n")

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
cat("  Prism object: OK\n")

cat("\n=== Step 7/9: Running BayesPrism (may take 10-30 min) ===\n")

t0 <- Sys.time()
bp_res <- run.prism(
  prism = myPrism,
  n.cores = n_cores,
  gibbs.control = list(chain.length = 2000, burn.in = 500, thinning = 2),
  opt.control = list(trace = 1, maxit = 10000)
)
runtime <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("  Runtime: %.1f minutes\n", as.numeric(runtime)))

saveRDS(bp_res, file.path(output_dir, "bp_res.rds"))
cat("  Saved: bp_res.rds\n")

# ==============================================================================
# 8. Extract Results
# ==============================================================================

cat("\n=== Step 8/9: Extracting Results ===\n")

# Cell type fractions
theta <- get.fraction(
  bp = bp_res,
  which.theta = "final",
  state.or.type = "type"
)
cat(sprintf("  theta: %d samples x %d cell types\n", nrow(theta), ncol(theta)))

# Coefficient of variation
theta_cv <- tryCatch(
  {
    if (!is.null(bp_res@posterior.theta_f@theta.cv)) {
      bp_res@posterior.theta_f@theta.cv
    } else if ("cv" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@cv
    } else if ("theta_cv" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@theta_cv
    } else {
      warning("theta CV not found, using zero placeholder")
      matrix(0, nrow(theta), ncol(theta), dimnames = dimnames(theta))
    }
  },
  error = function(e) {
    warning("theta CV extraction error: ", e$message)
    matrix(0, nrow(theta), ncol(theta), dimnames = dimnames(theta))
  }
)

# Log-likelihood
loglik_str <- tryCatch(
  {
    ll <- if (!is.null(bp_res@posterior.theta_f@loglikelihood)) {
      bp_res@posterior.theta_f@loglikelihood
    } else if ("logLik" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@logLik
    } else if ("log.lik" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@log.lik
    } else {
      NA
    }
    if (!is.na(ll)) sprintf("%.2f", ll) else "Not available"
  },
  error = function(e) "Not available"
)
cat("  Final log-likelihood:", loglik_str, "\n")

# Z matrix
Z <- extract_z_matrix(bp_res)
cat(sprintf("  Z matrix: %d genes x %d cell types\n", nrow(Z), ncol(Z)))

# Save core outputs
write.csv(theta, file.path(output_dir, "theta_final.csv"), quote = FALSE)
write.csv(theta_cv, file.path(output_dir, "theta_cv.csv"), quote = FALSE)
saveRDS(Z, file.path(output_dir, "Z_matrix.rds"), compress = "xz")
cat("  Saved: theta_final.csv, theta_cv.csv, Z_matrix.rds\n")

if (nrow(Z) > 1) {
  z_var <- apply(Z, 1, var)
  top_genes <- names(sort(z_var, decreasing = TRUE))[1:min(1000, nrow(Z))]
  write.csv(
    Z[top_genes, , drop = FALSE],
    file.path(output_dir, "Z_matrix_top1000genes.csv"),
    quote = FALSE
  )
  cat("  Saved: Z_matrix_top1000genes.csv\n")
}

cat("\n  Cell type fractions (theta):\n")
print(round(theta, 3))

# ==============================================================================
# 9. Marker Gene Validation
# ==============================================================================

cat("\n=== Step 9/9: Marker Gene Validation & Visualizations ===\n")

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

all_markers <- unique(unlist(marker_genes))
markers_present <- all_markers[all_markers %in% rownames(Z)]
cat(sprintf(
  "  Markers detected: %d / %d (%.1f%%)\n",
  length(markers_present),
  length(all_markers),
  100 * length(markers_present) / length(all_markers)
))

if (length(markers_present) > 0) {
  Z_marker <- Z[markers_present, , drop = FALSE]
  marker_annotation <- data.frame(
    CellType = rep(names(marker_genes), sapply(marker_genes, length))
  )
  rownames(marker_annotation) <- unlist(marker_genes)
  marker_annotation <- marker_annotation[rownames(Z_marker), , drop = FALSE]
  write.csv(Z_marker, file.path(output_dir, "Z_marker_genes.csv"))

  tryCatch(
    {
      pdf(file.path(output_dir, "Z_marker_heatmap.pdf"), width = 10, height = 8)
      pheatmap(
        log1p(Z_marker),
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = marker_annotation,
        color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
        border_color = NA,
        main = "Deconvolved Marker Gene Expression (log1p, row-scaled)",
        fontsize_row = 7,
        fontsize_col = 10
      )
      dev.off()
      cat("  Saved: Z_marker_genes.csv, Z_marker_heatmap.pdf\n")
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning("Marker heatmap error: ", e$message)
    }
  )
}

# ==============================================================================
# 10. Visualizations
# ==============================================================================

theta_dt <- as.data.table(theta, keep.rownames = "sample")
theta_long <- melt(
  theta_dt,
  id.vars = "sample",
  variable.name = "cell_type",
  value.name = "fraction"
)
n_ct <- length(unique(theta_long$cell_type))
ct_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_ct)
# Stacked bar plot
tryCatch(
  {
    p_stack <- ggplot(
      theta_long,
      aes(x = sample, y = fraction, fill = cell_type)
    ) +
      geom_col(width = 0.7, color = "black", linewidth = 0.3) +
      scale_fill_manual(values = ct_colors) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_bw(base_size = 12) +
      labs(
        title = "BayesPrism: Cell Type Fractions",
        x = "Bulk Sample",
        y = "Fraction",
        fill = "Cell Type"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    ggsave(
      file.path(output_dir, "theta_stacked_bar.pdf"),
      p_stack,
      width = 8,
      height = 5
    )
    cat("  Saved: theta_stacked_bar.pdf\n")
  },
  error = function(e) warning("Stacked bar error: ", e$message)
)

# Theta heatmap
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
    cat("  Saved: theta_heatmap.pdf\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Theta heatmap error: ", e$message)
  }
)

# Theta CV heatmap
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
    cat("  Saved: theta_cv_heatmap.pdf\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("CV heatmap error: ", e$message)
  }
)

# Fraction distribution boxplot
tryCatch(
  {
    p_box <- ggplot(
      theta_long,
      aes(x = cell_type, y = fraction, fill = cell_type)
    ) +
      geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(values = ct_colors) +
      theme_bw(base_size = 12) +
      labs(
        title = sprintf(
          "Cell Type Fraction Distribution (n=%d samples)",
          nrow(theta)
        ),
        x = "Cell Type",
        y = "Fraction"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    ggsave(
      file.path(output_dir, "theta_boxplot.pdf"),
      p_box,
      width = 8,
      height = 5
    )
    cat("  Saved: theta_boxplot.pdf\n")
  },
  error = function(e) warning("Boxplot error: ", e$message)
)

# Per-cell-type CV across samples
tryCatch(
  {
    cv_across <- apply(theta, 2, function(x) sd(x) / mean(x))
    cv_df <- data.frame(cell_type = names(cv_across), cv = cv_across)
    p_cv <- ggplot(cv_df, aes(x = reorder(cell_type, cv), y = cv)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(
        title = "Cell Type Variability Across Samples",
        x = "Cell Type",
        y = "Coefficient of Variation"
      )
    ggsave(
      file.path(output_dir, "cell_type_cv.pdf"),
      p_cv,
      width = 8,
      height = 6
    )
    cat("  Saved: cell_type_cv.pdf\n")
  },
  error = function(e) warning("CV barplot error: ", e$message)
)

# Cell type correlation heatmap
tryCatch(
  {
    cor_mat <- cor(theta)
    pdf(
      file.path(output_dir, "cell_type_correlation.pdf"),
      width = 10,
      height = 10
    )
    corrplot(
      cor_mat,
      method = "color",
      type = "upper",
      tl.col = "black",
      tl.srt = 45,
      addCoef.col = "black",
      number.cex = 0.7,
      title = "Cell Type Correlation",
      mar = c(0, 0, 2, 0)
    )
    dev.off()

    cor_pairs <- which(abs(cor_mat) > 0.7 & abs(cor_mat) < 1, arr.ind = TRUE)
    if (nrow(cor_pairs) > 0) {
      cor_pairs_df <- data.frame(
        CellType1 = rownames(cor_mat)[cor_pairs[, 1]],
        CellType2 = colnames(cor_mat)[cor_pairs[, 2]],
        Correlation = cor_mat[cor_pairs]
      )
      write.csv(
        cor_pairs_df,
        file.path(output_dir, "cell_type_correlation_pairs.csv"),
        row.names = FALSE
      )
    }
    cat("  Saved: cell_type_correlation.pdf, cell_type_correlation_pairs.csv\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Correlation plot error: ", e$message)
  }
)

# ==============================================================================
# 11. Summary Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

tryCatch(
  {
    sink(file.path(output_dir, "analysis_summary.txt"))

    cat(
      "================================================================================\n"
    )
    cat("BayesPrism Deconvolution Analysis Summary v3.0\n")
    cat(
      "================================================================================\n\n"
    )
    cat("Date:    ", as.character(Sys.Date()), "\n")
    cat(sprintf("Runtime: %.1f minutes\n\n", as.numeric(runtime)))

    cat("--- Reference ---\n")
    cat("File:             ", basename(seurat_rds_path), "\n")
    cat("Cells:            ", ncol(sce_ref), "\n")
    cat("Genes:            ", nrow(sce_ref), "\n")
    cat(
      "Annotation mode:  ",
      ifelse(use_single_level, "SINGLE level", "TWO levels"),
      "\n"
    )
    cat("Cell types:       ", length(unique(cell_type_labels)), "\n\n")

    cat("--- Bulk Data ---\n")
    cat("File:             ", basename(bulk_csv_path), "\n")
    cat("Samples:          ", nrow(bk_dat), "\n")
    cat("Sample IDs:       ", paste(rownames(bk_dat), collapse = ", "), "\n")
    cat(sprintf(
      "Lib size median:  %.0f  CV: %.2f\n\n",
      median(lib_size),
      lib_size_cv
    ))

    cat("--- Processing ---\n")
    cat("Common genes:     ", length(common_genes), "\n")
    cat(
      "GEP method:       ",
      ifelse(gep_use_mean, "mean-profile", "sum-profile"),
      "\n"
    )
    cat("Filter IG genes:  ", filter_immunoglobulin, "\n")
    cat("Outlier cut:      ", outlier_cut, "\n")
    cat("Gibbs chain:      2000 (burn-in 500, thinning 2)\n")
    cat("CPU cores:        ", n_cores, "\n\n")

    cat("--- Results ---\n")
    cat("Log-likelihood:   ", loglik_str, "\n\n")

    cat("Cell type fractions (theta):\n")
    print(round(theta, 4))
    cat("\nMean fractions:\n")
    print(round(colMeans(theta), 4))
    cat("\nCoefficient of variation (theta_cv):\n")
    print(round(theta_cv, 4))

    cat(sprintf(
      "\nMarker genes detected: %d / %d (%.1f%%)\n",
      length(markers_present),
      length(all_markers),
      100 * length(markers_present) / length(all_markers)
    ))

    cat("\n--- Output Files ---\n")
    files <- c(
      "bp_res.rds                       : BayesPrism result object",
      "theta_final.csv                  : Cell type fractions",
      "theta_cv.csv                     : Uncertainty (CV)",
      "Z_matrix.rds                     : Deconvolved expression (compressed)",
      "Z_matrix_top1000genes.csv        : Top 1000 variable genes",
      "Z_marker_genes.csv               : Marker gene expression",
      "theta_stacked_bar.pdf            : Stacked bar plot",
      "theta_heatmap.pdf                : Fraction heatmap",
      "theta_cv_heatmap.pdf             : Uncertainty heatmap",
      "theta_boxplot.pdf                : Fraction distribution",
      "Z_marker_heatmap.pdf             : Marker gene heatmap",
      "cell_type_cv.pdf                 : Per-type variability",
      "cell_type_correlation.pdf        : Cell type correlation heatmap",
      "cell_type_correlation_pairs.csv  : Highly correlated pairs (|r|>0.7)",
      "reference_metadata.rds           : Reference metadata",
      "analysis_summary.txt             : This summary",
      "session_info.txt                 : R session info"
    )
    for (i in seq_along(files)) {
      cat(sprintf("%2d. %s\n", i, files[i]))
    }

    cat(
      "\n================================================================================\n"
    )
    cat("Analysis completed successfully!\n")
    cat("Output directory:", output_dir, "\n")
    cat(
      "================================================================================\n"
    )
    sink()
  },
  error = function(e) {
    while (sink.number() > 0) {
      sink()
    }
    warning("Summary report error: ", e$message)
  },
  finally = {
    while (sink.number() > 0) {
      sink()
    }
  }
)
cat("  Saved: analysis_summary.txt\n")

# Session info
tryCatch(
  {
    sink(file.path(output_dir, "session_info.txt"))
    cat("R Session Information\n====================\n\n")
    print(sessionInfo())
    sink()
  },
  error = function(e) {
    while (sink.number() > 0) {
      sink()
    }
  },
  finally = {
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
cat("BayesPrism pipeline completed successfully!\n")
cat(sprintf("Runtime: %.1f minutes\n", as.numeric(runtime)))
cat("Output:  ", output_dir, "\n")
cat(rep("=", 80), "\n", sep = "")
cat("\nKey outputs:\n")
cat("  theta_stacked_bar.pdf          <- cell type composition\n")
cat("  theta_cv_heatmap.pdf           <- estimation uncertainty\n")
cat("  Z_marker_heatmap.pdf           <- marker gene validation\n")
cat("  cell_type_cv.pdf               <- per-type variability\n")
cat("  cell_type_correlation.pdf      <- inter-type correlation\n\n")


Z_test2 <- get.exp(bp = bp_res, state.or.type = "type")
cat("class:", class(Z_test2), "\n")
cat("dim:", dim(Z_test2), "\n")
cat("is.list:", is.list(Z_test2), "\n")
if (is.list(Z_test2)) {
  cat("length:", length(Z_test2), "\n")
  cat("names:", paste(head(names(Z_test2), 5), collapse = ", "), "\n")
  cat("class[[1]]:", class(Z_test2[[1]]), "\n")
  cat("dim[[1]]:", dim(Z_test2[[1]]), "\n")
}
