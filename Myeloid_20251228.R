# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
output_dir <- '/home/h2048/data/R/1228/Myeloid'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
library(reticulate)
library(harmony)
library(ggplot2)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
h5ad_file1 <- "/home/h2048/data/py/1217/cnmf_batch_production_v1_1_1/Myeloid/batch_aware/cnmf_analysis_k40_1/Myeloid_with_cnmf_k40.h5ad"
cat("Reading first h5ad file...\n")
seurat_obj <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)
seurat_obj <- NormalizeData(seurat_obj) #归一化
# ===== scANVI-based Analysis and QC =====
# Author: r2end
# Date: 2024-12-28
# Purpose: Visualize cells in scANVI space and compute markers by scANVI predictions

# ===== 1. Check Available Reductions =====
cat("Available reductions in Seurat object:\n")
print(names(seurat_obj@reductions))

# ===== 2. Set scANVI as Active Reduction =====
# Based on the log, we have: pca_, scanvi_, scvi_, umap_, umap_scanvi, umapscvi_
# We'll use umap_scanvi for visualization (UMAP computed on scANVI latent space)

if ("umap_scanvi" %in% names(seurat_obj@reductions)) {
  scanvi_reduction <- "umap_scanvi"
  cat("Using UMAP from scANVI space (umap_scanvi)\n")
} else if ("scanvi_" %in% names(seurat_obj@reductions)) {
  scanvi_reduction <- "scanvi_"
  cat("Using scANVI latent space directly (scanvi_)\n")
  # Compute UMAP on scANVI latent space if not present
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "scanvi_",
    dims = 1:30,
    reduction.name = "umap_scanvi",
    reduction.key = "umap_scanvi"
  )
  scanvi_reduction <- "umap_scanvi"
} else {
  stop("No scANVI reduction found in the object!")
}

# ===== 3. Visualize in scANVI Space =====
pdf("scanvi_visualization_qc.pdf", width = 16, height = 12)

# 3.1 scANVI predictions
p1 <- DimPlot(
  seurat_obj,
  reduction = scanvi_reduction,
  group.by = "scanvi_predictions",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("scANVI UMAP - Cell Type Predictions") +
  theme(legend.position = "right")

# 3.2 Original BBKNN clustering (for comparison)
p2 <- DimPlot(
  seurat_obj,
  reduction = scanvi_reduction,
  group.by = "leiden_bbknn_res1.0",
  label = TRUE,
  label.size = 3,
  raster = TRUE
) +
  ggtitle("scANVI UMAP - BBKNN Leiden Clusters")

# 3.3 Dataset/Batch
p3 <- DimPlot(
  seurat_obj,
  reduction = scanvi_reduction,
  group.by = "dataset",
  raster = TRUE
) +
  ggtitle("scANVI UMAP - Dataset")

# 3.4 Tissue
p4 <- DimPlot(
  seurat_obj,
  reduction = scanvi_reduction,
  group.by = "tissue",
  raster = TRUE
) +
  ggtitle("scANVI UMAP - Tissue")

print((p1 | p2) / (p3 | p4))

# 3.5 Sample distribution
p5 <- DimPlot(
  seurat_obj,
  reduction = scanvi_reduction,
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("scANVI UMAP - Sample") +
  theme(legend.position = "none")

print(p5)

# 3.6 Dominant GEP (cNMF program)
if ("dominant_gep" %in% colnames(seurat_obj@meta.data)) {
  p6 <- DimPlot(
    seurat_obj,
    reduction = scanvi_reduction,
    group.by = "dominant_gep",
    raster = TRUE
  ) +
    ggtitle("scANVI UMAP - Dominant cNMF Program")
  print(p6)
}

dev.off()

# ===== 4. Set scANVI Predictions as Identity =====
Idents(seurat_obj) <- "scanvi_predictions"

# Check cell type distribution
cat("\nCell type distribution (scANVI predictions):\n")
print(table(seurat_obj$scanvi_predictions))

# ===== 5. Quality Control Metrics by Cell Type =====
pdf("scanvi_celltype_qc_metrics.pdf", width = 14, height = 10)

# 5.1 nCount_RNA distribution
p_ncount <- VlnPlot(
  seurat_obj,
  features = "nCount_RNA",
  pt.size = 0,
  log = TRUE
) +
  ggtitle("nCount_RNA by scANVI Cell Type") +
  RotatedAxis()

# 5.2 nFeature_RNA distribution
p_nfeature <- VlnPlot(
  seurat_obj,
  features = "nFeature_RNA",
  pt.size = 0,
  log = TRUE
) +
  ggtitle("nFeature_RNA by scANVI Cell Type") +
  RotatedAxis()

# 5.3 Percent mitochondrial genes
if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
  p_pctmt <- VlnPlot(
    seurat_obj,
    features = "percent.mt",
    pt.size = 0
  ) +
    ggtitle("Percent Mitochondrial by scANVI Cell Type") +
    RotatedAxis()
  print(p_pctmt)
}

print(p_ncount)
print(p_nfeature)

# 5.4 Cell type proportions by sample
prop_table <- prop.table(
  table(seurat_obj$scanvi_predictions, seurat_obj$sample),
  margin = 2
)
prop_df <- as.data.frame(prop_table)
colnames(prop_df) <- c("CellType", "Sample", "Proportion")

p_prop <- ggplot(prop_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Type Proportions by Sample") +
  ylab("Proportion") +
  xlab("Sample")

print(p_prop)

# 5.5 Cell type proportions by tissue
if ("tissue" %in% colnames(seurat_obj@meta.data)) {
  prop_table_tissue <- prop.table(
    table(seurat_obj$scanvi_predictions, seurat_obj$tissue),
    margin = 2
  )
  prop_df_tissue <- as.data.frame(prop_table_tissue)
  colnames(prop_df_tissue) <- c("CellType", "Tissue", "Proportion")

  p_prop_tissue <- ggplot(
    prop_df_tissue,
    aes(x = Tissue, y = Proportion, fill = CellType)
  ) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Cell Type Proportions by Tissue") +
    ylab("Proportion") +
    xlab("Tissue")

  print(p_prop_tissue)
}

dev.off()

# ===== 6. Find Markers by scANVI Predictions (Wilcoxon) =====
cat("\nFinding markers for each scANVI predicted cell type...\n")
cat("Using Wilcoxon rank sum test with default parameters\n")

# Important: Make sure data is normalized
if (!"data" %in% names(seurat_obj@assays$RNA@layers)) {
  cat("Normalizing data...\n")
  seurat_obj <- NormalizeData(seurat_obj)
}

# Find all markers (one-vs-rest comparison)
all_markers <- FindAllMarkers(
  seurat_obj,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)

# Save markers
write.csv(all_markers, "scanvi_celltype_markers_wilcox.csv", row.names = FALSE)
cat("Markers saved to: scanvi_celltype_markers_wilcox.csv\n")

# Get top markers per cell type
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)

write.csv(top_markers, "scanvi_celltype_top10_markers.csv", row.names = FALSE)
cat(
  "Top 10 markers per cell type saved to: scanvi_celltype_top10_markers.csv\n"
)
# ===== 7. Enhanced Marker Visualization =====
pdf("scanvi_celltype_marker_visualization.pdf", width = 16, height = 14)

# 7.1 Prepare top markers
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)

top5_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5)

top3_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 3)

# 7.2 Heatmap of top 5 markers
if (nrow(top5_markers) > 0) {
  genes_to_plot <- unique(top5_markers$gene)
  
  cat(sprintf("Scaling %d genes for heatmap...\n", length(genes_to_plot)))
  
  # Scale only needed genes
  seurat_obj <- ScaleData(
    seurat_obj,
    features = genes_to_plot,
    verbose = FALSE
  )
  
  # Create heatmap
  p_heatmap <- DoHeatmap(
    seurat_obj,
    features = genes_to_plot,
    group.by = "scanvi_predictions",
    raster = TRUE,
    size = 3,
    draw.lines = TRUE
  ) +
    ggtitle("Top 5 Markers per scANVI Cell Type") +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "bottom"
    )
  
  print(p_heatmap)
}

# 7.3 DotPlot - Myeloid lineage markers
myeloid_markers <- list(
  "Pan_Myeloid" = c("PTPRC", "CD68", "CSF1R", "LYZ", "TYROBP"),
  "Monocyte_Classical" = c("S100A8", "S100A9", "S100A12", "FCN1", "VCAN", "IL1B"),
  "Monocyte_NonClassical" = c("FCGR3A", "MS4A7", "LST1", "IFITM3", "LGALS3"),
  "Macrophage" = c("APOE", "C1QA", "C1QB", "C1QC", "MERTK", "MAF", "FOLR2"),
  "DC_cDC1" = c("CLEC9A", "XCR1", "BATF3", "IRF8"),
  "DC_cDC2" = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  "DC_pDC" = c("GZMB", "IRF7", "IL3RA", "CLEC4C", "SERPINF1"),
  "DC_Migratory" = c("CCR7", "FSCN1", "CCL19", "CCL22", "LAMP3"),
  "Neutrophil" = c("FCGR3B", "CXCR2", "CXCR1", "CSF3R", "MMP25"),
  "Mast" = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC")
)

all_myeloid_markers <- unique(unlist(myeloid_markers))
present_markers <- intersect(all_myeloid_markers, rownames(seurat_obj))

if (length(present_markers) > 0) {
  p_dotplot <- DotPlot(
    seurat_obj,
    features = present_markers,
    assay = "RNA",
    dot.scale = 6,
    cluster.idents = FALSE
  ) +
    RotatedAxis() +
    ggtitle("Myeloid Lineage Marker Panel") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10)
    )
  
  print(p_dotplot)
}

# 7.4 Violin plots for key markers
key_violin_features <- c(
  "S100A8",    # Classical mono
  "FCGR3A",    # Non-classical mono
  "APOE",      # Macrophage
  "C1QA",      # Macrophage
  "CD1C",      # cDC2
  "IL3RA",     # pDC
  "CCR7",      # Migratory DC
  "TPSAB1"     # Mast
)

present_violin_features <- intersect(key_violin_features, rownames(seurat_obj))

if (length(present_violin_features) > 0) {
  for (feat in present_violin_features) {
    p_vln <- VlnPlot(
      seurat_obj,
      features = feat,
      group.by = "scanvi_predictions",
      pt.size = 0
    ) +
      ggtitle(paste("Expression of", feat, "by Cell Type")) +
      RotatedAxis() +
      theme(legend.position = "none")
    
    print(p_vln)
  }
}

# 7.5 Feature plots on scANVI UMAP
key_features <- c(
  "S100A8", "S100A9",      # Classical monocytes
  "FCGR3A", "MS4A7",       # Non-classical monocytes  
  "APOE", "C1QA",          # Macrophages
  "CD1C", "FCER1A",        # cDC2
  "CLEC9A", "XCR1",        # cDC1
  "GZMB", "IL3RA",         # pDC
  "CCR7", "FSCN1",         # Migratory DC
  "TPSAB1", "CPA3"         # Mast cells
)

present_key_features <- intersect(key_features, rownames(seurat_obj))

if (length(present_key_features) > 0) {
  # Plot in batches of 4
  for (i in seq(1, length(present_key_features), by = 4)) {
    feat_batch <- present_key_features[i:min(i + 3, length(present_key_features))]
    
    p_feature <- FeaturePlot(
      seurat_obj,
      features = feat_batch,
      reduction = scanvi_reduction,
      ncol = 2,
      raster = TRUE,
      order = TRUE,
      pt.size = 0.1
    ) &
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12)
      )
    
    print(p_feature)
  }
}

# 7.6 Top 3 markers overlay on UMAP (split by cell type)
if (nrow(top3_markers) > 0) {
  # For each cell type, show top 3 markers
  cell_types <- unique(top3_markers$cluster)
  
  for (ct in cell_types[1:min(3, length(cell_types))]) {  # Show first 3 cell types
    ct_markers <- top3_markers %>%
      filter(cluster == ct) %>%
      pull(gene)
    
    if (length(ct_markers) > 0) {
      p_ct <- FeaturePlot(
        seurat_obj,
        features = ct_markers,
        reduction = scanvi_reduction,
        ncol = 3,
        raster = TRUE,
        order = TRUE
      ) &
        ggtitle(paste("Top Markers -", ct)) &
        theme(plot.title = element_text(face = "bold"))
      
      print(p_ct)
    }
  }
}

dev.off()

cat("\nMarker visualization complete!\n")
# ===== 8. Statistical Summary =====
cat("\n===== Analysis Summary =====\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Total genes:", nrow(seurat_obj), "\n")
cat(
  "Number of scANVI predicted cell types:",
  length(unique(seurat_obj$scanvi_predictions)),
  "\n"
)
cat("Cell types:\n")
print(sort(table(seurat_obj$scanvi_predictions), decreasing = TRUE))

cat("\n===== Marker Gene Summary =====\n")
cat("Total significant markers:", nrow(all_markers), "\n")
marker_counts <- all_markers %>%
  group_by(cluster) %>%
  summarise(n_markers = n()) %>%
  arrange(desc(n_markers))
print(marker_counts)

# ===== 9. Save Updated Object =====
saveRDS(seurat_obj, "Myeloid_scanvi_analyzed_20251228.rds")
cat("\nUpdated Seurat object saved to: Myeloid_scanvi_analyzed_20251228.rds\n")

# ===== 10. Generate Comparison Table =====
# Compare BBKNN clusters with scANVI predictions
comparison_table <- table(
  BBKNN = seurat_obj$leiden_bbknn_res1.0,
  scANVI = seurat_obj$scanvi_predictions
)

write.csv(comparison_table, "bbknn_vs_scanvi_comparison.csv")
cat("BBKNN vs scANVI comparison saved to: bbknn_vs_scanvi_comparison.csv\n")

cat("\n===== Analysis Complete =====\n")
cat("Generated files:\n")
cat("  1. scanvi_visualization_qc.pdf\n")
cat("  2. scanvi_celltype_qc_metrics.pdf\n")
cat("  3. scanvi_celltype_markers_wilcox.csv\n")
cat("  4. scanvi_celltype_top10_markers.csv\n")
cat("  5. scanvi_celltype_marker_visualization.pdf\n")
cat("  6. bbknn_vs_scanvi_comparison.csv\n")
cat("  7. Myeloid_scanvi_analyzed_20251228.rds\n")
