# ===== Ciliated Cell Analysis with BBKNN Integration (Native R) =====
# Author: r2end
# Date: 2024-12-30
# Purpose:
# - Use native bbknnR package for batch correction (no Python dependency)
# - Leiden clustering on BBKNN graph
# - Wilcoxon marker detection with MT/RB genes filtered at marker stage only
# - Fully automated workflow

# ===== 1. Load Required Libraries =====
library(Seurat)
library(dplyr)
library(reticulate)
library(data.table)
library(ggplot2)
library(patchwork)
library(SCNT)
library(Matrix)
library(bbknnR) # Native R BBKNN implementation
use_condaenv("bbknn_env", required = TRUE)
py_config()

# ===== 2. Helper Functions =====

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

# ===== 3. Load Data and Setup Output Directory =====

H5AD_FILE <- '/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2/Ciliated/batch_aware/cnmf_analysis_k40/Ciliated_with_cnmf_k40.h5ad'


file_ext <- tools::file_ext(H5AD_FILE)

if (file_ext == "h5ad") {
  cat("Loading h5ad file...\n")
  seurat_obj <- GetSeurat(h5ad_path = H5AD_FILE, debug = TRUE)
} else if (file_ext == "rds") {
  cat("Loading RDS file...\n")
  seurat_obj <- readRDS(H5AD_FILE)
} else {
  stop(sprintf(
    "Unsupported file format: %s. Only .h5ad and .rds are supported.",
    file_ext
  ))
}

output_dir <- '/home/h2048/data/R/1230/ciliated'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

cat(sprintf("Initial cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Initial genes: %d\n", nrow(seurat_obj)))

# ===== 4. Gene Filtering (Except MT/RB - Keep for QC) =====
# Note: Keep MT/RB genes for QC metrics calculation
# They will be filtered only at the marker detection stage

all_genes <- rownames(seurat_obj)
genes_to_remove <- c()

# 1. Ribosomal pseudogenes only (not all ribosomal genes)
pseudo_genes <- grep(
  '^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$',
  all_genes,
  value = TRUE
)
genes_to_remove <- c(genes_to_remove, pseudo_genes)
cat(sprintf('Pseudogenes: %d\n', length(pseudo_genes)))

# 2. ENSG unannotated genes
ensg_genes <- grep('^ENSG[0-9]+', all_genes, value = TRUE)
genes_to_remove <- c(genes_to_remove, ensg_genes)
cat(sprintf('ENSG genes: %d\n', length(ensg_genes)))

# 3. Unannotated transcripts
unannotated_genes <- grep(
  '^(AC|AL|AP|BX|Z)[0-9]+\\.|^RP[0-9]+-|^CTD-|^CTB-|^CTC-|^LINC[0-9]+|-AS[0-9]+$|-OT[0-9]+$|^LOC[0-9]+',
  all_genes,
  value = TRUE
)
genes_to_remove <- c(genes_to_remove, unannotated_genes)
cat(sprintf('Unannotated transcripts: %d\n', length(unannotated_genes)))

# Remove duplicates and subset
genes_to_remove <- unique(genes_to_remove)
cat(sprintf(
  '\nTotal genes to remove (keeping MT/RB for QC): %d\n',
  length(genes_to_remove)
))

genes_to_keep <- setdiff(all_genes, genes_to_remove)
cat(sprintf('Genes to keep: %d\n', length(genes_to_keep)))

seurat_obj <- subset(seurat_obj, features = genes_to_keep)

# ===== 5. Additional Gene Filtering (min 3 cells) =====

cnt <- get_counts_matrix(seurat_obj)
gene_ncells <- Matrix::rowSums(cnt > 0)
keep_genes <- names(gene_ncells[gene_ncells >= 3])

before_genes <- nrow(seurat_obj)
seurat_obj <- subset(seurat_obj, features = keep_genes)
after_genes <- nrow(seurat_obj)

cat(sprintf("\nGenes before filtering (min 3 cells): %d\n", before_genes))
cat(sprintf("Genes after filtering: %d\n", after_genes))
cat(sprintf("Genes removed: %d\n", before_genes - after_genes))

# ===== 6. Calculate QC Metrics =====
# Must be done BEFORE any further gene filtering

if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
}
if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^(RPL|RPS)"
  )
}

cat(sprintf("\nQC metrics calculated:\n"))
cat(sprintf(
  "  - Mean percent.mt: %.2f%%\n",
  mean(seurat_obj$percent.mt, na.rm = TRUE)
))
cat(sprintf(
  "  - Mean percent.rb: %.2f%%\n",
  mean(seurat_obj$percent.rb, na.rm = TRUE)
))

# ===== 7. Standard Preprocessing =====

cat("\n===== Starting Standard Preprocessing =====\n")

# Normalization
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
cat("Normalization completed\n")

# Find variable features
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
cat(sprintf(
  "Variable features identified: %d\n",
  length(VariableFeatures(seurat_obj))
))

# Scale data
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
cat("Data scaling completed\n")

# PCA
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
cat("PCA completed (50 PCs)\n")

# ===== 8. BBKNN Batch Correction (Native R: bbknnR) =====

cat("\n===== Starting BBKNN Batch Correction (bbknnR) =====\n")

batch_key <- "study"
if (!batch_key %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Batch key '%s' not found in metadata", batch_key))
}

# Check batch sizes
batch_counts <- table(seurat_obj@meta.data[[batch_key]])
cat("\nCells per batch:\n")
print(batch_counts)

small_batches <- names(batch_counts[batch_counts < 3])
if (length(small_batches) > 0) {
  cat(sprintf(
    "\nWarning: %d batches have <3 cells (may cause issues)\n",
    length(small_batches)
  ))
  cat("Small batches:", paste(small_batches, collapse = ", "), "\n")
  cat("Consider filtering these batches if errors occur.\n")
}

# Run BBKNN with automatic UMAP computation
cat("\nRunning BBKNN integration...\n")
seurat_obj <- RunBBKNN(
  object = seurat_obj,
  batch_key = batch_key,
  reduction = "pca",
  n_pcs = 50L,

  # Graph naming
  graph_name = "bbknn",

  # UMAP settings
  run_UMAP = TRUE,
  UMAP_name = "umap_bbknn",
  UMAP_key = "UMAPBBKNN_",
  min_dist = 0.3,
  spread = 1.0,

  # Skip TSNE for speed
  run_TSNE = FALSE,

  # Random seed for reproducibility
  seed = 42,
  verbose = TRUE,

  # BBKNN parameters (passed to internal builder)
  neighbors_within_batch = 2,
  method = "annoy",
  metric = "euclidean",
  trim = NULL
)

cat("BBKNN integration completed\n")

# ===== 9. Clustering on BBKNN Graph =====

cat("\n===== Performing Clustering =====\n")

# bbknnR stores graph as "<assay>_<graph_name>", e.g. "RNA_bbknn"
graph_to_use <- paste0(DefaultAssay(seurat_obj), "_bbknn")

# Verify graph exists
if (!graph_to_use %in% names(seurat_obj@graphs)) {
  cat("Available graphs:\n")
  print(names(seurat_obj@graphs))
  stop(sprintf("BBKNN graph not found: %s", graph_to_use))
}

cat(sprintf("Using graph: %s\n", graph_to_use))

# Find clusters using Leiden algorithm
seurat_obj <- FindClusters(
  seurat_obj,
  graph.name = graph_to_use,
  algorithm = 4, # Leiden algorithm
  resolution = 1.0,
  verbose = FALSE
)

# Auto-detect cluster column name
cluster_col <- grep("^RNA_bbknn", colnames(seurat_obj@meta.data), value = TRUE)
if (length(cluster_col) == 0) {
  # Fallback: look for any new cluster column
  cluster_col <- grep("snn_res", colnames(seurat_obj@meta.data), value = TRUE)
}
cluster_col <- cluster_col[length(cluster_col)] # Use most recent

cat(sprintf("Cluster column: %s\n", cluster_col))
cat(sprintf(
  "Number of clusters: %d\n",
  length(unique(seurat_obj@meta.data[[cluster_col]]))
))

# Set identity
Idents(seurat_obj) <- cluster_col

cat("Clustering completed\n")

# ===== 10. Visualization =====

cat("\n===== Generating Visualizations =====\n")

pdf("01_bbknn_integration_qc.pdf", width = 14, height = 10)

# UMAP by clusters
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = cluster_col,
  label = TRUE,
  raster = TRUE
) +
  ggtitle(sprintf("BBKNN Integration - Clusters (%s)", cluster_col))

# UMAP by sample (batch)
p2 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("BBKNN Integration - Batch (Sample)")

# UMAP by dataset (if exists)
if ("dataset" %in% colnames(seurat_obj@meta.data)) {
  p3 <- DimPlot(
    seurat_obj,
    reduction = "umap_bbknn",
    group.by = "dataset",
    raster = TRUE
  ) +
    ggtitle("BBKNN Integration - Dataset")
} else {
  p3 <- NULL
}

# UMAP by tissue (if exists)
if ("tissue" %in% colnames(seurat_obj@meta.data)) {
  p4 <- DimPlot(
    seurat_obj,
    reduction = "umap_bbknn",
    group.by = "tissue",
    raster = TRUE
  ) +
    ggtitle("BBKNN Integration - Tissue")
} else {
  p4 <- NULL
}

# QC metrics
p5 <- FeaturePlot(
  seurat_obj,
  reduction = "umap_bbknn",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 2,
  raster = TRUE
)

print(p1)
print(p2)
if (!is.null(p3)) {
  print(p3)
}
if (!is.null(p4)) {
  print(p4)
}
print(p5)

dev.off()
cat("Visualization saved to 01_bbknn_integration_qc.pdf\n")

# ===== 11. QC by Cluster =====

pdf("02_QC_by_cluster.pdf", width = 14, height = 10)

# nFeature by cluster
p1 <- VlnPlot(
  seurat_obj,
  features = "nFeature_RNA",
  group.by = cluster_col,
  pt.size = 0,
  raster = TRUE
) +
  NoLegend() +
  ggtitle("nFeature_RNA by Cluster")

# nCount by cluster
p2 <- VlnPlot(
  seurat_obj,
  features = "nCount_RNA",
  group.by = cluster_col,
  pt.size = 0,
  raster = TRUE
) +
  NoLegend() +
  ggtitle("nCount_RNA by Cluster")

# percent.mt by cluster
p3 <- VlnPlot(
  seurat_obj,
  features = "percent.mt",
  group.by = cluster_col,
  pt.size = 0,
  raster = TRUE
) +
  NoLegend() +
  ggtitle("Mitochondrial % by Cluster")

# percent.rb by cluster
p4 <- VlnPlot(
  seurat_obj,
  features = "percent.rb",
  group.by = cluster_col,
  pt.size = 0,
  raster = TRUE
) +
  NoLegend() +
  ggtitle("Ribosomal % by Cluster")

print(p1)
print(p2)
print(p3)
print(p4)

dev.off()
cat("QC plots saved to 02_QC_by_cluster.pdf\n")

# ===== 12. FindAllMarkers with MT/RB Filtering =====

cat("\n===== Finding Cluster Markers (Wilcoxon Test) =====\n")

# Create marker-specific object by filtering MT/RB genes
current_genes <- rownames(seurat_obj)

# Identify MT and RB genes
mt_genes <- grep("^MT-", current_genes, value = TRUE)
rb_genes <- grep("^(RPL|RPS|MRPL|MRPS)", current_genes, value = TRUE)
tech_genes <- unique(c(mt_genes, rb_genes))

cat(sprintf("MT genes to exclude: %d\n", length(mt_genes)))
cat(sprintf("RB genes to exclude: %d\n", length(rb_genes)))
cat(sprintf("Total technical genes to exclude: %d\n", length(tech_genes)))

# Create subset for marker detection
marker_features <- setdiff(current_genes, tech_genes)
cat(sprintf("Genes for marker detection: %d\n", length(marker_features)))

seurat_marker <- subset(seurat_obj, features = marker_features)
Idents(seurat_marker) <- cluster_col

# Find all markers using Wilcoxon test
cat("Running FindAllMarkers (Wilcoxon)...\n")
markers_wilcox <- FindAllMarkers(
  seurat_marker,
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

cat(sprintf("Total markers found: %d\n", nrow(markers_wilcox)))

# Verify no MT/RB genes in results
if (any(markers_wilcox$gene %in% tech_genes)) {
  cat("Warning: Some MT/RB genes still in markers (should not happen)\n")
}

# Save all markers
write.csv(
  markers_wilcox,
  "03_all_markers_wilcox.csv",
  row.names = FALSE
)
cat("All markers saved to 03_all_markers_wilcox.csv\n")

# Get top markers per cluster
top_n <- 20
top_markers <- markers_wilcox %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  top_markers,
  "04_top20_markers_per_cluster_wilcox.csv",
  row.names = FALSE
)
cat(sprintf("Top %d markers per cluster saved\n", top_n))

# ===== 13. Marker Heatmap =====

cat("\n===== Generating Marker Heatmap =====\n")

# Get unique top markers (preserve order)
top_10_markers <- markers_wilcox %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

marker_genes <- top_10_markers$gene
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_marker)]
marker_genes <- marker_genes[!duplicated(marker_genes)]

cat(sprintf("Unique markers for heatmap: %d\n", length(marker_genes)))

# Scale only marker genes in marker object
seurat_marker <- ScaleData(
  seurat_marker,
  features = marker_genes,
  verbose = FALSE
)

# Generate heatmap
pdf("05_marker_heatmap_wilcox.pdf", width = 14, height = 16)

p <- DoHeatmap(
  seurat_marker,
  features = marker_genes,
  group.by = cluster_col,
  raster = TRUE
) +
  NoLegend() +
  ggtitle("Top 10 Markers per Cluster (Wilcoxon Test, MT/RB-free)")

print(p)
dev.off()
cat("Marker heatmap saved to 05_marker_heatmap_wilcox.pdf\n")

# ===== 14. Save Final Object =====

cat("\n===== Saving Final Seurat Object =====\n")

output_rds <- "ciliated_bbknn_final.rds"
saveRDS(seurat_obj, output_rds)
cat(sprintf("Final object saved to %s\n", output_rds))

# ===== 15. Summary Statistics =====

cat("\n===== Analysis Summary =====\n")
cat(sprintf("Final cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Final genes: %d\n", nrow(seurat_obj)))
cat(sprintf(
  "Number of clusters: %d\n",
  length(unique(seurat_obj@meta.data[[cluster_col]]))
))
cat(sprintf("Total markers detected: %d\n", nrow(markers_wilcox)))
cat(sprintf(
  "Batches integrated: %d\n",
  length(unique(seurat_obj@meta.data[[batch_key]]))
))

# Cluster sizes
cluster_sizes <- table(seurat_obj@meta.data[[cluster_col]])
cat("\nCluster sizes:\n")
print(sort(cluster_sizes, decreasing = TRUE))

# ===== 16. Session Info =====

cat("\n===== Session Information =====\n")
sessionInfo()

cat("\n===== Analysis Completed Successfully =====\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nGenerated files:\n")
cat("  - 01_bbknn_integration_qc.pdf\n")
cat("  - 02_QC_by_cluster.pdf\n")
cat("  - 03_all_markers_wilcox.csv\n")
cat("  - 04_top20_markers_per_cluster_wilcox.csv\n")
cat("  - 05_marker_heatmap_wilcox.pdf\n")
cat("  - ciliated_bbknn_final.rds\n")
