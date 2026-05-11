# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
setwd("/home/h2048/data/R/1218")
library(reticulate)
library(harmony)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
h5ad_file1 <- "/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/annotation_results/adata_with_manual_annotations_raw.h5ad"
cat("Reading first h5ad file...\n")
seurat_obj <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)
seurat_obj <- NormalizeData(seurat_obj) #归一化
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
) #寻找变异基因
seurat_obj <- ScaleData(seurat_obj) #标准化


# 0) 选择聚类/分组作为身份（按需改成你的meta列名，如 "cell_type"）
Idents(seurat_obj) <- "leiden_bbknn"

# 1) FindAllMarkers
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, "all_markers.csv", row.names = FALSE)

# 2) 每个cluster取Top N marker
top_n <- 10
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE)

marker_genes <- unique(top_markers$gene)
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]

# 3) 仅对这些marker做Scale + 热图
seurat_obj <- ScaleData(seurat_obj, features = marker_genes, verbose = FALSE)

p <- DoHeatmap(
  seurat_obj,
  features = marker_genes,
  group.by = "leiden_bbknn",
  raster = TRUE
) +
  NoLegend()

# ggsave("marker_heatmap.pdf", p, width = 10, height = 12)
pdf("epi_marker.pdf", width = 12, height = 8, onefile = TRUE)
print(p)
dev.off()

# 确保有常见 QC 指标（如已存在会跳过）
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

cluster_col <- "leiden_bbknn"

pdf("QC_by_cluster.pdf", width = 12, height = 8, onefile = TRUE)

# 1) Violin：按 cluster 看 nCount / nFeature / mt / rb
p1 <- VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
  group.by = cluster_col,
  pt.size = 0.05,
  ncol = 2
)
print(p1)

# 2) Scatter：nCount vs nFeature（按 cluster 上色）
print(FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = cluster_col
))

# 3) Scatter：percent.mt vs nCount / nFeature（帮助定位低质/高线粒体群）
print(FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = cluster_col
))
print(FeatureScatter(
  seurat_obj,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt",
  group.by = cluster_col
))

dev.off()

p0 <- VlnPlot(
  seurat_obj,
  features = c('PTPRC', 'LST1', 'FCGR3B'),
  group.by = cluster_col,
  pt.size = 0.05,
  ncol = 2
)
print(p0)

pdf(
  "QC_by_cluster_oneplot_per_page.pdf",
  width = 12,
  height = 6,
  onefile = TRUE
)

qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
for (f in qc_feats) {
  p <- VlnPlot(
    seurat_obj,
    features = f,
    group.by = cluster_col,
    pt.size = 0
  )
  print(p) # 每次 print 自动新开一页
}

dev.off()

anno <- fread("/home/h2048/data/R/1217/Annotation.csv") # 你的CSV：含 Cluster, Size, Percent, Suggestion, Dominant_GEP, Manual_Annotation
anno[, Cluster := as.character(Cluster)]

cluster_col <- "leiden_bbknn" # 如果你的cluster列不是这个名，改成对应列名
cur_cluster <- as.character(seurat_obj[[cluster_col, drop = TRUE]])

# 把 CSV 里除 Cluster 外的所有列，按 cluster 映射到每个细胞
for (col in setdiff(colnames(anno), "Cluster")) {
  seurat_obj[[col]] <- anno[[col]][match(cur_cluster, anno$Cluster)]
}

# 可选：把 Manual_Annotation 设为当前 ident
seurat_obj <- SetIdent(seurat_obj, value = "Manual_Annotation")

# quick check
table(seurat_obj$Manual_Annotation, useNA = "ifany")
pdf("celltype.pdf", width = 12, height = 8, onefile = TRUE)

p <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "Manual_Annotation",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
)
print(p)

dev.off()
saveRDS(seurat_obj, 'epithelial_bbknn_raw_20251217.rds')
Idents(seurat_obj) <- "Manual_Annotation"

# （可选但推荐）明确用哪个assay做marker
DefaultAssay(seurat_obj) <- "RNA"

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, "all_markers.csv", row.names = FALSE)

top_n <- 10
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

# 保留顺序去重（而不是 unique() 打乱顺序）
marker_genes <- top_markers$gene
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]
marker_genes <- marker_genes[!duplicated(marker_genes)]

seurat_obj <- ScaleData(seurat_obj, features = marker_genes, verbose = FALSE)

p <- DoHeatmap(
  seurat_obj,
  features = marker_genes,
  group.by = "Manual_Annotation", # 也可以删掉，默认按 Idents
  raster = TRUE
) +
  NoLegend()

pdf("epi_subtype_marker.pdf", width = 12, height = 8, onefile = TRUE)
print(p)
dev.off()

seurat_obj <- readRDS(
  '/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds'
)
seurat_sub <- subset(seurat_obj, subset = Manual_Annotation %in% 'Basal')
seurat_obj <- NormalizeData(seurat_obj) #归一化
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
) #寻找变异基因
seurat_obj <- ScaleData(seurat_obj) #标准化


suppressMessages(library(ROGUE))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
OUTPUT_DIR <- "rogue_results"
GROUPING_COL <- "leiden_bbknn" # Change to desired metadata column
SAMPLE_COL <- "batch" # Sample/batch column for per-sample ROGUE
MIN_CELLS_PER_GROUP <- 10 # Filter small groups
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===== Get Cell Type Groups =====
meta_data <- seurat_obj@meta.data
cell_types <- unique(meta_data[[GROUPING_COL]])
cell_types <- cell_types[!is.na(cell_types)]

cat(sprintf("Found %d cell types in '%s'\n", length(cell_types), GROUPING_COL))

# ===== Calculate ROGUE for Each Cell Type =====
rogue_results <- list()

for (ct in cell_types) {
  cat(sprintf("\nProcessing: %s... ", ct))

  tryCatch(
    {
      # Subset cells for this cell type
      cells_keep <- meta_data[[GROUPING_COL]] == ct
      n_cells <- sum(cells_keep)

      cat(sprintf("(%d cells) ", n_cells))

      if (n_cells < MIN_CELLS_PER_GROUP) {
        cat("SKIPPED (too few cells)\n")
        next
      }

      # Extract expression matrix for this cell type only
      expr_subset <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
      expr_subset <- expr_subset[, cells_keep]
      expr_subset <- as.matrix(expr_subset)

      # Filter low-abundance genes and cells
      expr_subset <- matr.filter(expr_subset, min.cells = 10, min.genes = 200)

      # Check matrix validity
      if (ncol(expr_subset) < MIN_CELLS_PER_GROUP || nrow(expr_subset) < 100) {
        cat("SKIPPED (insufficient genes/cells after filtering)\n")
        rm(expr_subset)
        gc(verbose = FALSE)
        next
      }

      # Calculate entropy
      ent_res <- SE_fun(expr_subset)

      # Check for invalid values
      if (any(is.na(ent_res$entropy)) || any(is.infinite(ent_res$entropy))) {
        cat("SKIPPED (invalid entropy values)\n")
        rm(expr_subset, ent_res)
        gc(verbose = FALSE)
        next
      }

      # Calculate ROGUE
      rogue_val <- CalculateRogue(ent_res, platform = "UMI")

      cat(sprintf("ROGUE = %.4f\n", rogue_val))

      # Store results
      rogue_results[[ct]] <- data.frame(
        cell_type = ct,
        n_cells = n_cells,
        n_genes_used = nrow(expr_subset),
        rogue_value = rogue_val
      )

      # Clean up
      rm(expr_subset, ent_res)
      gc(verbose = FALSE)
    },
    error = function(e) {
      cat(sprintf("ERROR: %s\n", e$message))
    }
  )
}

# ===== Compile Results =====
rogue_df <- bind_rows(rogue_results)

if (nrow(rogue_df) == 0) {
  stop("No ROGUE values calculated. Check data and parameters.")
}

# Add purity classification
rogue_df <- rogue_df %>%
  mutate(
    purity_class = case_when(
      rogue_value >= 0.9 ~ "High (≥0.9)",
      rogue_value >= 0.7 ~ "Moderate (0.7-0.9)",
      TRUE ~ "Low (<0.7)"
    )
  ) %>%
  arrange(desc(rogue_value))

# Save results
write.csv(
  rogue_df,
  file.path(OUTPUT_DIR, "rogue_values_by_celltype.csv"),
  row.names = FALSE
)

# ===== Print Summary =====
cat("\n=== ROGUE Summary ===\n")
print(rogue_df, row.names = FALSE)

# ===== Visualizations =====

# 1. Bar plot with values
pdf(file.path(OUTPUT_DIR, "rogue_barplot.pdf"), width = 10, height = 6)
p1 <- ggplot(
  rogue_df,
  aes(x = reorder(cell_type, rogue_value), y = rogue_value, fill = purity_class)
) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.3f", rogue_value)),
    hjust = -0.1,
    size = 3.5
  ) +
  geom_hline(
    yintercept = c(0.7, 0.9),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "High (≥0.9)" = "forestgreen",
      "Moderate (0.7-0.9)" = "orange",
      "Low (<0.7)" = "firebrick"
    )
  ) +
  labs(
    title = "ROGUE Cluster Purity by Cell Type",
    subtitle = sprintf("Based on '%s' annotation", GROUPING_COL),
    x = "Cell Type",
    y = "ROGUE Value (Higher = More Pure)",
    fill = "Purity Class"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10), legend.position = "bottom") +
  ylim(0, 1.05)
print(p1)
dev.off()

# 2. Lollipop plot with cell counts
pdf(file.path(OUTPUT_DIR, "rogue_lollipop.pdf"), width = 10, height = 6)
p2 <- ggplot(
  rogue_df,
  aes(x = reorder(cell_type, rogue_value), y = rogue_value)
) +
  geom_segment(
    aes(xend = reorder(cell_type, rogue_value), yend = 0),
    color = "grey50"
  ) +
  geom_point(aes(size = n_cells, color = purity_class), alpha = 0.8) +
  geom_hline(
    yintercept = c(0.7, 0.9),
    linetype = "dashed",
    color = "red",
    alpha = 0.3
  ) +
  coord_flip() +
  scale_color_manual(
    values = c(
      "High (≥0.9)" = "forestgreen",
      "Moderate (0.7-0.9)" = "orange",
      "Low (<0.7)" = "firebrick"
    )
  ) +
  scale_size_continuous(range = c(3, 10), labels = scales::comma) +
  labs(
    title = "ROGUE Values with Cell Counts",
    subtitle = sprintf("Calculated from '%s' annotation", GROUPING_COL),
    x = "Cell Type",
    y = "ROGUE Value",
    color = "Purity Class",
    size = "Number of Cells"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10), legend.position = "right")
print(p2)
dev.off()

# 3. Summary table plot
pdf(file.path(OUTPUT_DIR, "rogue_table.pdf"), width = 12, height = 8)
rogue_table <- rogue_df %>%
  mutate(
    rogue_value = sprintf("%.4f", rogue_value),
    n_cells = format(n_cells, big.mark = ","),
    n_genes_used = format(n_genes_used, big.mark = ",")
  )

gridExtra::grid.table(rogue_table, rows = NULL)
dev.off()

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat(sprintf(
  "Successfully calculated ROGUE for %d/%d cell types\n",
  nrow(rogue_df),
  length(cell_types)
))

seurat_obj <- readRDS(
  '/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds'
)
seurat_sub <- subset(seurat_obj, subset = Manual_Annotation %in% 'Basal')
seurat_sub <- NormalizeData(seurat_sub) #归一化
seurat_sub <- FindVariableFeatures(
  seurat_sub,
  selection.method = "vst",
  nfeatures = 2000
) #寻找变异基因
seurat_sub <- ScaleData(seurat_sub) #标准化

# 使用高变基因进行主成分分析，降低数据维度
seurat_sub <- RunPCA(seurat_sub, npcs = 50)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
seurat_sub <- RunHarmony(
  object = seurat_sub,
  group.by.vars = c("sample"),
  theta = c(1), # Higher theta for more diverse clustering
  lambda = c(7), # Higher lambda to reduce overcorrection
  sigma = 0.1, # Lower sigma for tighter clusters
  nclust = 30, # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20,
  early_stop = TRUE,
  dims = 1:50 # More iterations for better convergence
)
seurat_sub <- RunUMAP(
  seurat_sub,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  n.trees = 500,
  min.dist = 0.4,
  # ,
  # learning.rate = 0.2, # 相对保守的学习率
  # n.epochs = 1400, # 增加迭代次数补偿较小的学习率
  # spread = 1.2,
  # repulsion.strength = 1.1,

  metric = "correlation"
)
# Find neighbors
seurat_sub <- FindNeighbors(
  seurat_sub,
  reduction = "harmony",
  dims = 1:30,
  k.param = 45
)

seurat_sub <- FindClusters(
  seurat_sub,
  algorithm = 4,
  group.singletons = FALSE,
  resolution = 3, # 多个分辨率
  verbose = TRUE
)
