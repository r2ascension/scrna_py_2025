# ===== Ciliated Cell Analysis with BBKNN Integration =====
# Author: r2end
# Date: 2024-12-30
# Purpose:
# - Replace Harmony with BBKNN for batch correction
# - Use Wilcoxon test for marker detection
# - Filter MT/RB genes before marker analysis

# ===== 1. Load Required Libraries =====
library(Seurat)
library(reticulate)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(Matrix)

# Configure Python environment
use_condaenv("bbknn_env", required = TRUE)
py_config()

# Import Python modules
sc <- import("scanpy")
bbknn <- import("bbknn")
np <- import("numpy")

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

seurat_obj <- readRDS(
  '/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2/Ciliated/batch_aware/cnmf_analysis_k40/Ciliated_with_cnmf_k40.h5ad'
)

output_dir <- '/home/h2048/data/R/1230/ciliated'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

cat(sprintf("Initial cells: %d\n", ncol(seurat_obj)))
cat(sprintf("Initial genes: %d\n", nrow(seurat_obj)))

# ===== 4. Gene Filtering =====
# Remove mitochondrial, ribosomal, pseudogenes, ENSG, and unannotated genes

all_genes <- rownames(seurat_obj)
genes_to_remove <- c()

# Mitochondrial genes (MT-)
mt_genes <- grep('^MT-', all_genes, value = TRUE)
genes_to_remove <- c(genes_to_remove, mt_genes)
cat(sprintf('Mitochondrial genes: %d\n', length(mt_genes)))

# Ribosomal genes (RPS, RPL, MRPS, MRPL)
ribo_genes <- grep('^RPS|^RPL|^MRPS|^MRPL', all_genes, value = TRUE)
genes_to_remove <- c(genes_to_remove, ribo_genes)
cat(sprintf('Ribosomal genes: %d\n', length(ribo_genes)))

# Ribosomal pseudogenes
pseudo_genes <- grep(
  '^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$',
  all_genes,
  value = TRUE
)
genes_to_remove <- c(genes_to_remove, pseudo_genes)
cat(sprintf('Pseudogenes: %d\n', length(pseudo_genes)))

# ENSG unannotated genes
ensg_genes <- grep('^ENSG[0-9]+', all_genes, value = TRUE)
genes_to_remove <- c(genes_to_remove, ensg_genes)
cat(sprintf('ENSG genes: %d\n', length(ensg_genes)))

# Unannotated transcripts
unannotated_genes <- grep(
  '^(AC|AL|AP|BX|Z)[0-9]+\\.|^RP[0-9]+-|^CTD-|^CTB-|^CTC-|^LINC[0-9]+|-AS[0-9]+$|-OT[0-9]+$|^LOC[0-9]+',
  all_genes,
  value = TRUE
)
genes_to_remove <- c(genes_to_remove, unannotated_genes)
cat(sprintf('Unannotated transcripts: %d\n', length(unannotated_genes)))

# Remove duplicates and subset
genes_to_remove <- unique(genes_to_remove)
cat(sprintf('\nTotal genes to remove: %d\n', length(genes_to_remove)))

genes_to_keep <- setdiff(all_genes, genes_to_remove)
cat(sprintf('Genes to keep: %d\n', length(genes_to_keep)))

seurat_obj <- subset(seurat_obj, features = genes_to_keep)

# ===== 5. Add Manual Annotations =====

anno <- fread("/home/h2048/data/R/1217/Annotation.csv")
anno[, Cluster := as.character(Cluster)]

cluster_col <- "RNA_snn_res.1"
cur_cluster <- as.character(seurat_obj[[cluster_col, drop = TRUE]])

# Map annotations to cells
for (col in setdiff(colnames(anno), "Cluster")) {
  seurat_obj[[col]] <- anno[[col]][match(cur_cluster, anno$Cluster)]
}

seurat_obj <- SetIdent(seurat_obj, value = "Manual_Annotation")
cat("\nManual annotation distribution:\n")
print(table(seurat_obj$Manual_Annotation, useNA = "ifany"))

# Visualize cell types
pdf("01_initial_celltype_umap.pdf", width = 12, height = 8)
p <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "Manual_Annotation",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("Initial Cell Type Annotations")
print(p)
dev.off()

# ===== 6. Additional Gene Filtering (min 3 cells) =====

cnt <- get_counts_matrix(seurat_obj)
gene_ncells <- Matrix::rowSums(cnt > 0)
keep_genes <- names(gene_ncells[gene_ncells >= 3])

before_genes <- nrow(seurat_obj)
seurat_obj <- subset(seurat_obj, features = keep_genes)
after_genes <- nrow(seurat_obj)

cat(sprintf("\nGenes before filtering (min 3 cells): %d\n", before_genes))
cat(sprintf("Genes after filtering: %d\n", after_genes))
cat(sprintf("Genes removed: %d\n", before_genes - after_genes))

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

# ===== 8. BBKNN Batch Correction =====

cat("\n===== Starting BBKNN Batch Correction =====\n")

# Convert Seurat to AnnData
cat("Converting Seurat to AnnData...\n")

# Extract necessary components
counts_mat <- get_counts_matrix(seurat_obj)
normalized_mat <- GetAssayData(seurat_obj, slot = "data")
scaled_mat <- GetAssayData(seurat_obj, slot = "scale.data")
pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")
meta_data <- seurat_obj@meta.data

# Create AnnData object
adata <- sc$AnnData(
  X = t(normalized_mat),
  obs = meta_data,
  var = data.frame(
    gene_ids = rownames(normalized_mat),
    row.names = rownames(normalized_mat)
  )
)

# Add layers
adata$layers$update(list(
  counts = t(counts_mat),
  scaled = t(scaled_mat[, colnames(seurat_obj)])
))

# Add PCA
adata$obsm$update(list(
  X_pca = pca_embeddings
))

cat(sprintf(
  "AnnData created: %d cells x %d genes\n",
  as.integer(adata$n_obs),
  as.integer(adata$n_vars)
))

# Check batch key
batch_key <- "sample"
if (!batch_key %in% colnames(adata$obs)) {
  stop(sprintf("Batch key '%s' not found in metadata", batch_key))
}

# Count cells per batch
batch_counts <- table(adata$obs[[batch_key]])
cat("\nCells per batch:\n")
print(batch_counts)

# Filter small batches (< 3 cells)
small_batches <- names(batch_counts[batch_counts < 3])
if (length(small_batches) > 0) {
  cat(sprintf(
    "\nWarning: Removing %d small batches (< 3 cells)\n",
    length(small_batches)
  ))
  keep_mask <- !adata$obs[[batch_key]] %in% small_batches
  adata <- adata[keep_mask, ]
  cat(sprintf("Cells after filtering: %d\n", as.integer(adata$n_obs)))
}

# Run BBKNN
cat("\nRunning BBKNN integration...\n")
bbknn$bbknn(
  adata,
  batch_key = batch_key,
  neighbors_within_batch = 5L,
  n_pcs = 50L,
  trim = NULL,
  copy = FALSE
)
cat("BBKNN integration completed\n")

# Compute UMAP
cat("Computing UMAP...\n")
sc$tl$umap(adata, min_dist = 0.3, spread = 1.0)
cat("UMAP completed\n")

# ===== 9. Transfer Results Back to Seurat =====

cat("\n===== Transferring Results to Seurat =====\n")

# Filter Seurat object to match AnnData cells
valid_cells <- py_to_r(adata$obs_names$to_list())
seurat_obj <- subset(seurat_obj, cells = valid_cells)

# Add BBKNN UMAP
umap_coords <- py_to_r(adata$obsm$get("X_umap"))
rownames(umap_coords) <- valid_cells
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

seurat_obj[["umap_bbknn"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAPBBKNN_",
  assay = "RNA"
)
cat("UMAP coordinates transferred\n")

# Extract BBKNN connectivities and distances
connectivities <- py_to_r(adata$obsp$get("connectivities"))
distances <- py_to_r(adata$obsp$get("distances"))

# Convert to Graph objects
bbknn_graph <- as.Graph(connectivities)
bbknn_dist <- as.Graph(distances)

# Store in Seurat object
seurat_obj[["bbknn"]] <- bbknn_graph
seurat_obj[["bbknn_dist"]] <- bbknn_dist
cat("BBKNN graphs transferred\n")

# ===== 10. Clustering on BBKNN Graph =====

cat("\n===== Performing Clustering =====\n")

# Find clusters using BBKNN graph
seurat_obj <- FindClusters(
  seurat_obj,
  graph.name = "bbknn",
  algorithm = 4, # Leiden algorithm
  resolution = 1.0,
  verbose = FALSE
)
cat("Clustering completed\n")

# ===== 11. QC Metrics =====

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

# ===== 12. Visualization =====

cat("\n===== Generating Visualizations =====\n")

pdf("02_bbknn_integration_qc.pdf", width = 14, height = 10)

# UMAP by clusters
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = "bbknn_snn_res.1",
  label = TRUE,
  raster = TRUE
) +
  ggtitle("BBKNN Integration - Clusters (res=1.0)")

# UMAP by sample (batch)
p2 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("BBKNN Integration - Batch (Sample)")

# UMAP by dataset
p3 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = "dataset",
  raster = TRUE
) +
  ggtitle("BBKNN Integration - Dataset")

# UMAP by tissue
p4 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = "tissue",
  raster = TRUE
) +
  ggtitle("BBKNN Integration - Tissue")

# UMAP by manual annotation
p5 <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  group.by = "Manual_Annotation",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("BBKNN Integration - Manual Annotation")

# QC metrics
p6 <- FeaturePlot(
  seurat_obj,
  reduction = "umap_bbknn",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 2,
  raster = TRUE
)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

dev.off()
cat("Visualization saved to 02_bbknn_integration_qc.pdf\n")

# ===== 13. QC by Cluster =====

pdf("03_QC_by_cluster.pdf", width = 14, height = 10)

cluster_col <- "bbknn_snn_res.1"

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
cat("QC plots saved to 03_QC_by_cluster.pdf\n")

# ===== 14. FindAllMarkers with Wilcoxon Test =====

cat("\n===== Finding Cluster Markers (Wilcoxon Test) =====\n")

# Set identity to clusters
Idents(seurat_obj) <- "bbknn_snn_res.1"

# Verify no MT or RB genes remain
current_genes <- rownames(seurat_obj)
remaining_mt <- grep("^MT-", current_genes, value = TRUE)
remaining_rb <- grep("^(RPL|RPS|MRPL|MRPS)", current_genes, value = TRUE)

cat(sprintf("Remaining MT genes: %d\n", length(remaining_mt)))
cat(sprintf("Remaining RB genes: %d\n", length(remaining_rb)))

if (length(remaining_mt) > 0 | length(remaining_rb) > 0) {
  cat("Warning: Some MT/RB genes still present, filtering now...\n")
  genes_to_remove_final <- c(remaining_mt, remaining_rb)
  genes_to_keep_final <- setdiff(current_genes, genes_to_remove_final)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep_final)
  cat(sprintf("Final gene count: %d\n", nrow(seurat_obj)))
}

# Find all markers using Wilcoxon test
cat("Running FindAllMarkers (Wilcoxon)...\n")
markers_wilcox <- FindAllMarkers(
  seurat_obj,
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

cat(sprintf("Total markers found: %d\n", nrow(markers_wilcox)))

# Save all markers
write.csv(
  markers_wilcox,
  "04_all_markers_wilcox.csv",
  row.names = FALSE
)
cat("All markers saved to 04_all_markers_wilcox.csv\n")

# Get top markers per cluster
top_n <- 20
top_markers <- markers_wilcox %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  top_markers,
  "05_top20_markers_per_cluster_wilcox.csv",
  row.names = FALSE
)
cat(sprintf("Top %d markers per cluster saved\n", top_n))

# ===== 15. Marker Heatmap =====

cat("\n===== Generating Marker Heatmap =====\n")

# Get unique top markers (preserve order)
top_10_markers <- markers_wilcox %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

marker_genes <- top_10_markers$gene
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]
marker_genes <- marker_genes[!duplicated(marker_genes)]

cat(sprintf("Unique markers for heatmap: %d\n", length(marker_genes)))

# Scale only marker genes
seurat_obj <- ScaleData(
  seurat_obj,
  features = marker_genes,
  verbose = FALSE
)

# Generate heatmap
pdf("06_marker_heatmap_wilcox.pdf", width = 14, height = 16)

p <- DoHeatmap(
  seurat_obj,
  features = marker_genes,
  group.by = "bbknn_snn_res.1",
  raster = TRUE
) +
  NoLegend() +
  ggtitle("Top 10 Markers per Cluster (Wilcoxon Test)")

print(p)
dev.off()
cat("Marker heatmap saved to 06_marker_heatmap_wilcox.pdf\n")

# ===== 16. Save Final Object =====

cat("\n===== Saving Final Seurat Object =====\n")

output_rds <- "ciliated_bbknn_final.rds"
saveRDS(seurat_obj, output_rds)
cat(sprintf("Final object saved to %s\n", output_rds))

# ===== 17. Session Info =====

cat("\n===== Session Information =====\n")
sessionInfo()

cat("\n===== Analysis Completed Successfully =====\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nGenerated files:\n")
cat("  - 01_initial_celltype_umap.pdf\n")
cat("  - 02_bbknn_integration_qc.pdf\n")
cat("  - 03_QC_by_cluster.pdf\n")
cat("  - 04_all_markers_wilcox.csv\n")
cat("  - 05_top20_markers_per_cluster_wilcox.csv\n")
cat("  - 06_marker_heatmap_wilcox.pdf\n")
cat("  - ciliated_bbknn_final.rds\n")
