# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
setwd("/home/h2048/data/R/1218")
library(reticulate)
library(cowplot)
library(harmony)
library(ggplot2)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
h5ad_file1 <- "/home/h2048/data/py/1217/bcell/bcell_with_bbknn.h5ad"
cat("Reading first h5ad file...\n")
seurat_obj <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)
table(seurat_obj$dataset)
# seurat_obj <- subset(seurat_obj,subset= dataset %in% c('Kerstin_B_Meyer_2021'))
seurat_obj <- NormalizeData(seurat_obj) #归一化
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
) #寻找变异基因
seurat_obj <- ScaleData(seurat_obj) #标准化
# Idents(seurat_obj) <- "leiden_bbknn_res1.5"
Idents(seurat_obj) <- "scanvi_predictions"

p <- DimPlot(
  seurat_obj,
  reduction = "umap_bbknn",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  NoLegend()

print(p)
# 1) FindAllMarkers
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, "b_all_markers.csv", row.names = FALSE)

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
pdf("b_marker.pdf", width = 12, height = 8, onefile = TRUE)
print(p)
dev.off()

age_associated_marker <- unique(c(
  'LST1',
  'TYROBP',
  'FCER1G',
  'S100A8',
  'S100A9',
  'LGALS3',
  'CTSS',
  'LYZ',
  'FCGR3A',
  'LILRB1',
  'MZB1',
  'JCHAIN',
  'XBP1',
  'SDC1',
  'PRDM1',
  'LST1',
  'TYROBP',
  'LYZ',
  'MS4A1',
  'CD79A',
  'CD74',
  'FCRL4',
  'FCRL5',
  'ITGAX'
))
p <- FeaturePlot(
  seurat_obj,
  features = age_associated_marker,
  raster = TRUE
)
print(p)

p1 <- DotPlot(
  seurat_obj,
  features = age_associated_marker,
  group.by = "scanvi_predictions"
)
print(p1)

p1data <- p1$data
print(p1data)

p2 <- (FeatureScatter(
  seurat_obj,
  feature1 = "MS4A1",
  feature2 = "LST1",
  group.by = 'scanvi_predictions'
))

FeatureScatter(seurat_obj, feature1 = "MS4A1", feature2 = "LST1") +
  ggplot2::theme_classic()

FeatureScatter(seurat_obj, feature1 = "CD79A", feature2 = "LST1") +
  ggplot2::theme_classic()

FeatureScatter(seurat_obj, feature1 = "CD79A", feature2 = "MS4A1") +
  ggplot2::theme_classic()

FeatureScatter(seurat_obj, feature1 = "CD79A", feature2 = "TYROBP") +
  ggplot2::theme_classic()

print(p2)

Idents(seurat_obj) <- "scanvi_predictions"
abc <- subset(seurat_obj, idents = "Age-associated B cells")

# 以表达>0 作为“阳性”阈值（你也可以用 >1 或 >2 做敏感性分析）
b_pos <- FetchData(abc, c("MS4A1", "CD79A"))
my_pos <- FetchData(abc, c("LST1", "TYROBP", "LYZ"))

prop.table(
  table(
    B_core = (b_pos$MS4A1 > 0 & b_pos$CD79A > 0),
    Myeloid = (my_pos$LST1 > 0 | my_pos$TYROBP > 0 | my_pos$LYZ > 0)
  ),
  margin = 1
)


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
getwd()
marker_genes <- c(
  'MS4A1',
  'FCER2',
  'CR2',
  'TCL1A',
  'IGHD',
  'CD27',
  'TNFRSF13B',
  'AIM2',
  'MZB1',
  'XBP1',
  'JCHAIN',
  'PRDM1',
  'SDC1',
  'HLA-DRA',
  'CD74',
  'CD83',
  'CD86'
)
Idents(seurat_obj) <- "leiden_bbknn_res1.0"
seurat_obj <- ScaleData(seurat_obj, features = marker_genes, verbose = FALSE)

p <- DoHeatmap(
  seurat_obj,
  features = marker_genes,
  group.by = "leiden_bbknn_res1.0",
  raster = TRUE
) +
  NoLegend()
p1 <- VlnPlot(seurat_obj, features = marker_genes, pt.size = 0, ncol = 3)

p2 <- DotPlot(
  seurat_obj,
  features = marker_genes,
  group.by = "leiden_bbknn_res1.0",
  split.by = NULL,
  cols = c("lightgrey", "red"),
  dot.scale = 8
)
# ggsave("marker_heatmap.pdf", p, width = 10, height = 12)
source_data <- p2$data
write.csv(source_data, "VlnPlot.csv")
pdf("epi_marker.pdf", width = 48, height = 32, onefile = TRUE)
print(p)
print(p1)
dev.off()
getwd()
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
suppressMessages(library(ROGUE))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
setwd('/home/h2048/data/R/1218/B')
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


DefaultAssay(seurat_obj) <- "RNA"

meta_cols <- colnames(seurat_obj@meta.data)
GROUP_COL <- "scanvi_predictions"

# 选择一个可用的降维（优先 umap）
red_avail <- Reductions(seurat_obj)
REDUCTION <- if ("umap" %in% red_avail) "umap" else red_avail[1]

# =========================
# 1) Myeloid signature genes（可按需增删）
# =========================
myeloid_genes <- c(
  # “骨架/经典髓系”
  "LST1",
  "TYROBP",
  "FCER1G",
  "LYZ",
  "S100A8",
  "S100A9",
  "FCGR3A",
  "LILRB1",
  # “单核/巨噬常见模块”
  "AIF1",
  "SPI1",
  "CSF1R",
  "MS4A7",
  "FCN1",
  "VCAN",
  "CSTA",
  "MNDA",
  "C1QA",
  "C1QB",
  "C1QC",
  # “抗原呈递/溶酶体”（非髓系特异，但与髓系常共变）
  "CTSS",
  "CTSD",
  "LGALS3",
  "IFITM3"
)
myeloid_genes <- intersect(myeloid_genes, rownames(seurat_obj))
if (length(myeloid_genes) < 8) {
  stop("myeloid_genes 与对象基因集交集过少，请检查基因命名。")
}

# =========================
# 2) 计算 Myeloid module score
# =========================
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(myeloid_genes),
  name = "MyeloidScore"
)
SCORE_COL <- "MyeloidScore1"

# =========================
# 3) 定义“myeloid-high”并做全局统计
#    （用全体细胞95分位做阈值；你也可改成 0.90/0.99）
# =========================
thr <- as.numeric(quantile(
  seurat_obj[[SCORE_COL]][, 1],
  probs = 0.95,
  na.rm = TRUE
))
seurat_obj$myeloid_high <- seurat_obj[[SCORE_COL]][, 1] >= thr

tab <- as.data.frame.matrix(table(
  group = seurat_obj[[GROUP_COL]][, 1],
  myeloid_high = seurat_obj$myeloid_high
))
tab$group <- rownames(tab)
rownames(tab) <- NULL
if (!all(c("FALSE", "TRUE") %in% colnames(tab))) {
  # 防止极端情况只有一个水平
  if (!("TRUE" %in% colnames(tab))) tab$
  TRUE <- 0
  if (!("FALSE" %in% colnames(tab))) tab$
  FALSE <- 0
}
colnames(tab)[colnames(tab) == "FALSE"] <- "myeloid_low"
colnames(tab)[colnames(tab) == "TRUE"] <- "myeloid_high"
tab$total <- tab$myeloid_low + tab$myeloid_high
tab$frac_myeloid_high <- ifelse(
  tab$total > 0,
  tab$myeloid_high / tab$total,
  NA_real_
)
tab <- tab[order(tab$frac_myeloid_high, decreasing = TRUE), ]

write.csv(tab, file = "myeloid_high_fraction_by_group.csv", row.names = FALSE)

# =========================
# 4) 可视化：多页PDF，每页一个图；不画单细胞散点
# =========================
pdf("myeloid_overview.pdf", width = 11, height = 8.5)

# 4.1 全局 MyeloidScore 在 UMAP 上的分布
p1 <- FeaturePlot(
  seurat_obj,
  features = SCORE_COL,
  reduction = REDUCTION,
  raster = TRUE,
  order = TRUE
)
print(p1)

# 4.2 MyeloidScore 按 cell type/cluster 的分布（不画点）
p2 <- VlnPlot(
  seurat_obj,
  features = SCORE_COL,
  group.by = GROUP_COL,
  pt.size = 0
) +
  RotatedAxis()
print(p2)

# 4.3 关键髓系骨架基因：每基因一页 FeaturePlot（更直观排查污染/混群）
key_genes <- intersect(
  c("LST1", "TYROBP", "LYZ", "FCER1G", "S100A8", "S100A9", "FCGR3A", "LILRB1"),
  rownames(seurat_obj)
)
for (g in key_genes) {
  pg <- FeaturePlot(
    seurat_obj,
    features = g,
    reduction = REDUCTION,
    raster = TRUE,
    order = TRUE
  ) +
    ggtitle(sprintf("%s on %s", g, REDUCTION))
  print(pg)
}

# 4.4 DotPlot：在所有 group 上看一组髓系基因的表达/阳性比例（整体一页）
dot_genes <- intersect(
  c(
    "LST1",
    "TYROBP",
    "LYZ",
    "FCER1G",
    "S100A8",
    "S100A9",
    "FCGR3A",
    "MS4A7",
    "FCN1",
    "C1QC",
    "CTSS",
    "LGALS3"
  ),
  rownames(seurat_obj)
)
p4 <- DotPlot(
  seurat_obj,
  features = dot_genes,
  group.by = GROUP_COL
) +
  RotatedAxis() +
  ggtitle(sprintf("Myeloid gene panel across %s", GROUP_COL))
print(p4)

print(p4$data)

dev.off()

# =========================
# 5) 控制台快速查看 top groups
# =========================
print(head(tab, 20))


DefaultAssay(seurat_obj) <- "RNA"
GROUP_COL <- "scanvi_predictions" # 你按实际替换

my_core <- intersect(
  c(
    "LST1",
    "TYROBP",
    "FCER1G",
    "FCGR3A",
    "S100A8",
    "S100A9",
    "MS4A7",
    "FCN1",
    "C1QC",
    "LYZ"
  ),
  rownames(seurat_obj)
)
apc_state <- intersect(c("CTSS", "LGALS3"), rownames(seurat_obj))

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(my_core),
  name = "MyeloidCoreScore"
)
summarise(seurat_obj$MyeloidCoreScore)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(apc_state),
  name = "APCStateScore"
)

pdf("myeloid_core_vs_state.pdf", width = 11, height = 8.5)
print(
  VlnPlot(
    seurat_obj,
    features = c("MyeloidCoreScore1", "APCStateScore1"),
    group.by = GROUP_COL,
    pt.size = 0
  ) +
    RotatedAxis()
)

pdc_markers <- intersect(
  c("IL3RA", "GZMB", "TCF4", "CLEC4C", "LILRA4", "IRF7"),
  rownames(seurat_obj)
)
if (length(pdc_markers) > 0) {
  print(
    DotPlot(seurat_obj, features = pdc_markers, group.by = GROUP_COL) +
      RotatedAxis() +
      ggtitle("pDC canonical markers")
  )
}
dev.off()

p0 <- DotPlot(seurat_obj, features = pdc_markers, group.by = GROUP_COL) +
  RotatedAxis() +
  ggtitle("pDC canonical markers")
p0$data


suppressPackageStartupMessages({
  library(Seurat)
})

DefaultAssay(seurat_obj) <- "RNA"
GROUP_COL <- "scanvi_predictions" # 按你的实际列名

# --- gene sets ---
mono_core <- intersect(
  c(
    "S100A8",
    "S100A9",
    "LST1",
    "FCN1",
    "VCAN",
    "MS4A7",
    "LILRB1",
    "FCGR3A",
    "TYROBP",
    "FCER1G",
    "CSTA",
    "MNDA",
    "AIF1",
    "CSF1R",
    "C1QA",
    "C1QB",
    "C1QC",
    "CTSS",
    "LGALS3"
  ),
  rownames(seurat_obj)
)

pdc_core <- intersect(
  c("IL3RA", "GZMB", "TCF4", "IRF7", "CLEC4C", "LILRA4"),
  rownames(seurat_obj)
)
b_core <- intersect(
  c("MS4A1", "CD79A", "CD74", "CD37", "CD22", "BANK1"),
  rownames(seurat_obj)
)

seurat_obj <- AddModuleScore(seurat_obj, list(mono_core), name = "MonoCore")
seurat_obj <- AddModuleScore(seurat_obj, list(pdc_core), name = "pDCCore")
seurat_obj <- AddModuleScore(seurat_obj, list(b_core), name = "BCore")

# thresholds (global p95)
thr_mono <- as.numeric(quantile(
  seurat_obj$MonoCore1,
  probs = 0.95,
  na.rm = TRUE
))
thr_pdc <- as.numeric(quantile(seurat_obj$pDCCore1, probs = 0.95, na.rm = TRUE))
thr_b <- as.numeric(quantile(seurat_obj$BCore1, probs = 0.95, na.rm = TRUE))

seurat_obj$mono_high <- seurat_obj$MonoCore1 >= thr_mono
seurat_obj$pdc_high <- seurat_obj$pDCCore1 >= thr_pdc
seurat_obj$b_high <- seurat_obj$BCore1 >= thr_b

# summary by group (no dplyr required)
grp <- seurat_obj[[GROUP_COL]][, 1]
df <- data.frame(
  group = grp,
  MonoCore = seurat_obj$MonoCore1,
  pDCCore = seurat_obj$pDCCore1,
  BCore = seurat_obj$BCore1,
  mono_high = seurat_obj$mono_high,
  pdc_high = seurat_obj$pdc_high,
  b_high = seurat_obj$b_high
)

summ <- do.call(
  rbind,
  lapply(split(df, df$group), function(x) {
    data.frame(
      group = x$group[1],
      n = nrow(x),
      MonoCore_mean = mean(x$MonoCore, na.rm = TRUE),
      pDCCore_mean = mean(x$pDCCore, na.rm = TRUE),
      BCore_mean = mean(x$BCore, na.rm = TRUE),
      frac_mono_high = mean(x$mono_high, na.rm = TRUE),
      frac_pdc_high = mean(x$pdc_high, na.rm = TRUE),
      frac_b_high = mean(x$b_high, na.rm = TRUE)
    )
  })
)

summ <- summ[order(summ$frac_mono_high, decreasing = TRUE), ]
print(summ)

write.csv(summ, "score_summary_by_group.csv", row.names = FALSE)
cat("Saved: score_summary_by_group.csv\n")

DefaultAssay(seurat_obj) <- "RNA"
GROUP_COL <- "scanvi_predictions" # 按实际替换

mono_strict <- intersect(
  c("S100A8", "S100A9", "FCN1", "VCAN", "LYZ"),
  rownames(seurat_obj)
)
mac_strict <- intersect(
  c("C1QA", "C1QB", "C1QC", "APOE", "LGMN", "CTSD"),
  rownames(seurat_obj)
)
pdc_strict <- intersect(
  c("IL3RA", "GZMB", "CLEC4C", "LILRA4"),
  rownames(seurat_obj)
)
b_core <- intersect(
  c("MS4A1", "CD79A", "CD74", "CD37", "CD22", "BANK1"),
  rownames(seurat_obj)
)

seurat_obj <- AddModuleScore(seurat_obj, list(mono_strict), name = "MonoStrict")
seurat_obj <- AddModuleScore(seurat_obj, list(mac_strict), name = "MacStrict")
seurat_obj <- AddModuleScore(seurat_obj, list(pdc_strict), name = "pDCStrict")
seurat_obj <- AddModuleScore(seurat_obj, list(b_core), name = "BCore")

thr_mono <- as.numeric(quantile(
  seurat_obj$MonoStrict1,
  probs = 0.95,
  na.rm = TRUE
))
thr_mac <- as.numeric(quantile(
  seurat_obj$MacStrict1,
  probs = 0.95,
  na.rm = TRUE
))
thr_pdc <- as.numeric(quantile(
  seurat_obj$pDCStrict1,
  probs = 0.95,
  na.rm = TRUE
))
thr_b <- as.numeric(quantile(seurat_obj$BCore1, probs = 0.95, na.rm = TRUE))

seurat_obj$mono_strict_high <- seurat_obj$MonoStrict1 >= thr_mono
seurat_obj$mac_strict_high <- seurat_obj$MacStrict1 >= thr_mac
seurat_obj$pdc_strict_high <- seurat_obj$pDCStrict1 >= thr_pdc
seurat_obj$b_high <- seurat_obj$BCore1 >= thr_b

grp <- seurat_obj[[GROUP_COL]][, 1]
df <- data.frame(
  group = grp,
  mono = seurat_obj$MonoStrict1,
  mac = seurat_obj$MacStrict1,
  pdc = seurat_obj$pDCStrict1,
  b = seurat_obj$BCore1,
  mono_high = seurat_obj$mono_strict_high,
  mac_high = seurat_obj$mac_strict_high,
  pdc_high = seurat_obj$pdc_strict_high,
  b_high = seurat_obj$b_high
)

summ2 <- do.call(
  rbind,
  lapply(split(df, df$group), function(x) {
    data.frame(
      group = x$group[1],
      n = nrow(x),
      MonoStrict_mean = mean(x$mono, na.rm = TRUE),
      MacStrict_mean = mean(x$mac, na.rm = TRUE),
      pDCStrict_mean = mean(x$pdc, na.rm = TRUE),
      BCore_mean = mean(x$b, na.rm = TRUE),
      frac_mono_high = mean(x$mono_high, na.rm = TRUE),
      frac_mac_high = mean(x$mac_high, na.rm = TRUE),
      frac_pdc_high = mean(x$pdc_high, na.rm = TRUE),
      frac_b_high = mean(x$b_high, na.rm = TRUE)
    )
  })
)

summ2 <- summ2[order(summ2$frac_mono_high, decreasing = TRUE), ]
print(summ2)
write.csv(summ2, "score_summary_strict_by_group.csv", row.names = FALSE)

DefaultAssay(seurat_obj) <- "RNA"
GROUP_COL <- "scanvi_predictions"

mono_nolyz <- intersect(
  c("S100A8", "S100A9", "FCN1", "VCAN"),
  rownames(seurat_obj)
)
if (length(mono_nolyz) < 2) {
  stop("mono_nolyz genes too few; check gene symbols.")
}

seurat_obj <- AddModuleScore(seurat_obj, list(mono_nolyz), name = "MonoNoLYZ")

thr_nolyz <- as.numeric(quantile(
  seurat_obj$MonoNoLYZ1,
  probs = 0.95,
  na.rm = TRUE
))
seurat_obj$mono_nolyz_high <- seurat_obj$MonoNoLYZ1 >= thr_nolyz

# 只看 Plasma cells 内部比例
pl <- subset(seurat_obj, subset = scanvi_predictions %in% c("Plasma cells"))
cat("Plasma MonoStrict_high (your previous):", mean(pl$mono_strict_high), "\n")
cat("Plasma MonoNoLYZ_high:", mean(pl$mono_nolyz_high), "\n")

# 顺手看一下 Plasma 内关键骨架的 pct（不用图）
genes_check <- intersect(
  c(
    "LYZ",
    "FCN1",
    "VCAN",
    "S100A8",
    "S100A9",
    "LST1",
    "TYROBP",
    "FCER1G",
    "MS4A7"
  ),
  rownames(pl)
)
dp <- DotPlot(pl, features = genes_check, group.by = "mono_strict_high")
print(dp$data[order(dp$data$features.plot, dp$data$id), ])

DefaultAssay(seurat_obj) <- "RNA"

# 只取 Plasma
pl <- subset(seurat_obj, subset = scanvi_predictions %in% "Plasma cells")

# Plasma core（尽量用强、相对特异的）
plasma_core <- intersect(
  c("JCHAIN", "MZB1", "XBP1", "SDC1", "PRDM1"),
  rownames(pl)
)

# Myeloid skeleton（明确排除 LYZ）
my_skel <- intersect(
  c(
    "LST1",
    "TYROBP",
    "FCER1G",
    "MS4A7",
    "FCGR3A",
    "FCN1",
    "VCAN",
    "S100A8",
    "S100A9"
  ),
  rownames(pl)
)

# Fetch
m_plasma <- FetchData(pl, plasma_core)
m_myel <- FetchData(pl, my_skel)

# 你可以把阈值从 >0 改成 >1（更严格）
plasma_pos <- rowSums(m_plasma > 0) >= 2 # 至少2个浆细胞核心阳性
myel_pos <- rowSums(m_myel > 0) >= 2 # 至少2个髓系骨架阳性（不含LYZ）

cat("Plasma core positive fraction:", mean(plasma_pos), "\n")
cat("Myeloid skeleton positive fraction:", mean(myel_pos), "\n")
cat(
  "Putative plasma-myeloid doublet fraction:",
  mean(plasma_pos & myel_pos),
  "\n"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# =========================
# Config
# =========================
SAMPLE_COL <- "dataset" # 你的“样本”列名
CELLTYPE_COL <- "scanvi_predictions" # 或 Manual_Annotation
DATASET_COL <- NULL # 如需同时保留 study/dataset 层级：填入列名；否则设 NULL

MIN_CELLS_PER_GROUP <- 50
OUTPUT_DIR <- OUTPUT_DIR # 你已有
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# meta_data 建议直接来自 seurat_obj@meta.data，且 rownames 是 cell names
meta_data <- seurat_obj@meta.data

stopifnot(SAMPLE_COL %in% colnames(meta_data))
stopifnot(CELLTYPE_COL %in% colnames(meta_data))
if (!is.null(DATASET_COL)) {
  stopifnot(DATASET_COL %in% colnames(meta_data))
}

# =========================
# Get counts ONCE (Seurat v5/v4 compatible)
# =========================
counts <- tryCatch(
  GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
)

# =========================
# Build combinations: sample × celltype (optionally dataset)
# =========================
group_cols <- c(SAMPLE_COL, CELLTYPE_COL)
if (!is.null(DATASET_COL)) {
  group_cols <- c(DATASET_COL, group_cols)
}

combo_counts <- meta_data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(across(all_of(group_cols))) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= MIN_CELLS_PER_GROUP)

cat(sprintf("Total valid combinations: %d\n", nrow(combo_counts)))
cat(sprintf(
  "Samples: %d, Cell types: %d\n",
  n_distinct(combo_counts[[SAMPLE_COL]]),
  n_distinct(combo_counts[[CELLTYPE_COL]])
))

# =========================
# Helper: safe ROGUE
# =========================
calc_rogue_safe <- function(expr_mat, min_cells = MIN_CELLS_PER_GROUP) {
  # expr_mat: genes × cells, raw counts
  # Return: list(rogue=..., n_genes=...)
  if (is.null(expr_mat) || ncol(expr_mat) < min_cells) {
    return(list(rogue = NA_real_, n_genes = NA_integer_))
  }

  # matr.filter 可能改变维度；并可能报错
  expr_f <- tryCatch(
    matr.filter(as.matrix(expr_mat), min.cells = 10, min.genes = 200),
    error = function(e) NULL
  )
  if (is.null(expr_f) || ncol(expr_f) < min_cells || nrow(expr_f) < 100) {
    return(list(
      rogue = NA_real_,
      n_genes = if (!is.null(expr_f)) nrow(expr_f) else NA_integer_
    ))
  }

  ent_res <- tryCatch(SE_fun(expr_f), error = function(e) NULL)
  if (
    is.null(ent_res) ||
      any(is.na(ent_res$entropy)) ||
      any(is.infinite(ent_res$entropy))
  ) {
    return(list(rogue = NA_real_, n_genes = nrow(expr_f)))
  }

  rogue_val <- tryCatch(
    CalculateRogue(ent_res, platform = "UMI"),
    error = function(e) NA_real_
  )
  list(rogue = rogue_val, n_genes = nrow(expr_f))
}

# =========================
# Main loop: per sample × celltype
# =========================
cat("\n=== Calculating per-sample ROGUE values ===\n")
rogue_results_list <- vector("list", nrow(combo_counts))
pb <- txtProgressBar(max = nrow(combo_counts), style = 3)

for (i in seq_len(nrow(combo_counts))) {
  n_cells <- combo_counts$n_cells[i]

  # 取出本组合的 key
  key <- combo_counts[i, group_cols, drop = FALSE]

  # 用 cell names 做 subset（避免 meta 顺序问题）
  cells_keep <- meta_data %>%
    tibble::rownames_to_column("cell") %>%
    semi_join(
      key %>% mutate(.tmp_join = 1),
      by = setNames(group_cols, group_cols)
    ) %>%
    pull(cell)

  # 进一步保证这些 cell 在 counts 里
  cells_keep <- intersect(cells_keep, colnames(counts))

  rogue_out <- tryCatch(
    {
      expr_subset <- counts[, cells_keep, drop = FALSE]
      calc_rogue_safe(expr_subset, min_cells = MIN_CELLS_PER_GROUP)
    },
    error = function(e) list(rogue = NA_real_, n_genes = NA_integer_)
  )

  # 组装输出
  row_out <- as.data.frame(key, stringsAsFactors = FALSE)
  row_out$n_cells <- length(cells_keep)
  row_out$n_genes_used <- rogue_out$n_genes
  row_out$rogue_value <- rogue_out$rogue
  rogue_results_list[[i]] <- row_out

  setTxtProgressBar(pb, i)
}
close(pb)

rogue_df <- bind_rows(rogue_results_list)

# =========================
# Save outputs
# =========================
write.csv(
  rogue_df,
  file.path(OUTPUT_DIR, "rogue_values_per_sample_long.csv"),
  row.names = FALSE
)

# sample × celltype wide matrix（如你需要 heatmap）
rogue_wide <- rogue_df %>%
  select(all_of(c(group_cols, "rogue_value"))) %>%
  tidyr::pivot_wider(
    names_from = all_of(CELLTYPE_COL),
    values_from = rogue_value
  )

write.csv(
  rogue_wide,
  file.path(OUTPUT_DIR, "rogue_values_per_sample_wide.csv"),
  row.names = FALSE
)

# Summary by cell type (across samples)
celltype_summary <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  group_by(across(all_of(CELLTYPE_COL))) %>%
  summarise(
    n_samples = n(),
    total_cells = sum(n_cells),
    mean_rogue = mean(rogue_value),
    sd_rogue = sd(rogue_value),
    median_rogue = median(rogue_value),
    min_rogue = min(rogue_value),
    max_rogue = max(rogue_value),
    .groups = "drop"
  ) %>%
  arrange(desc(median_rogue))

write.csv(
  celltype_summary,
  file.path(OUTPUT_DIR, "rogue_summary_by_celltype_per_sample.csv"),
  row.names = FALSE
)

# Summary by sample
sample_summary <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  group_by(across(all_of(SAMPLE_COL))) %>%
  summarise(
    n_celltypes = n(),
    total_cells = sum(n_cells),
    mean_rogue = mean(rogue_value),
    median_rogue = median(rogue_value),
    .groups = "drop"
  ) %>%
  arrange(desc(median_rogue))

write.csv(
  sample_summary,
  file.path(OUTPUT_DIR, "rogue_summary_by_sample.csv"),
  row.names = FALSE
)

cat("\n=== Done ===\n")
cat(sprintf("Saved to: %s\n", OUTPUT_DIR))
cat(sprintf(
  "Successfully calculated: %d/%d combinations\n",
  sum(!is.na(rogue_df$rogue_value)),
  nrow(rogue_df)
))


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(ROGUE)
})
getwd()
# =========================
# Config
# =========================
CELLTYPE_COL <- "scanvi_predictions" # 你的列名
ASSAY_USE <- "RNA"

MIN_CELLS_PER_GROUP <- 200
MAX_CELLS_PER_GROUP <- 10000 # 太大容易 as.matrix 爆内存；不想下采样就设 Inf
SEED <- 1
OUTPUT_DIR
OUTPUT_DIR <- "/home/h2048/data/R/1218"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

meta_data <- seurat_obj@meta.data
stopifnot(CELLTYPE_COL %in% colnames(meta_data))

set.seed(SEED)

# =========================
# Get counts ONCE (Seurat v5/v4 compatible)
# =========================
counts <- tryCatch(
  GetAssayData(seurat_obj, assay = ASSAY_USE, layer = "counts"),
  error = function(e) {
    GetAssayData(seurat_obj, assay = ASSAY_USE, slot = "counts")
  }
)

# =========================
# Build celltype groups
# =========================
grp_df <- meta_data %>%
  rownames_to_column("cell") %>%
  filter(!is.na(.data[[CELLTYPE_COL]]), .data[[CELLTYPE_COL]] != "") %>%
  group_by(.data[[CELLTYPE_COL]]) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells)) %>%
  filter(n_cells >= MIN_CELLS_PER_GROUP)

cat(sprintf(
  "Valid cell types (n >= %d): %d\n",
  MIN_CELLS_PER_GROUP,
  nrow(grp_df)
))

# split cells by celltype（只对通过阈值的 celltype）
cells_by_type <- meta_data %>%
  rownames_to_column("cell") %>%
  filter(.data[[CELLTYPE_COL]] %in% grp_df[[CELLTYPE_COL]]) %>%
  group_by(.data[[CELLTYPE_COL]]) %>%
  summarise(cells = list(cell), .groups = "drop")

# =========================
# Helper: safe ROGUE
# =========================
calc_rogue_safe <- function(expr_mat, min_cells = MIN_CELLS_PER_GROUP) {
  if (is.null(expr_mat) || ncol(expr_mat) < min_cells) {
    return(list(
      rogue = NA_real_,
      n_genes = NA_integer_,
      n_cells = ifelse(is.null(expr_mat), 0L, ncol(expr_mat))
    ))
  }

  # (可选) 下采样，避免 dense 爆内存
  if (is.finite(MAX_CELLS_PER_GROUP) && ncol(expr_mat) > MAX_CELLS_PER_GROUP) {
    keep <- sample(colnames(expr_mat), MAX_CELLS_PER_GROUP)
    expr_mat <- expr_mat[, keep, drop = FALSE]
  }

  # matr.filter / SE_fun 在很多实现里要求 base matrix；这里转 dense（子集后再转）
  expr_dense <- tryCatch(as.matrix(expr_mat), error = function(e) NULL)
  if (is.null(expr_dense)) {
    return(list(
      rogue = NA_real_,
      n_genes = NA_integer_,
      n_cells = ncol(expr_mat)
    ))
  }

  expr_f <- tryCatch(
    matr.filter(expr_dense, min.cells = 10, min.genes = 200),
    error = function(e) NULL
  )
  if (is.null(expr_f) || ncol(expr_f) < min_cells || nrow(expr_f) < 100) {
    return(list(
      rogue = NA_real_,
      n_genes = if (!is.null(expr_f)) nrow(expr_f) else NA_integer_,
      n_cells = ncol(expr_dense)
    ))
  }

  ent_res <- tryCatch(SE_fun(expr_f), error = function(e) NULL)
  if (
    is.null(ent_res) ||
      any(is.na(ent_res$entropy)) ||
      any(is.infinite(ent_res$entropy))
  ) {
    return(list(
      rogue = NA_real_,
      n_genes = nrow(expr_f),
      n_cells = ncol(expr_f)
    ))
  }

  rogue_val <- tryCatch(
    CalculateRogue(ent_res, platform = "UMI"),
    error = function(e) NA_real_
  )
  list(rogue = rogue_val, n_genes = nrow(expr_f), n_cells = ncol(expr_f))
}

# =========================
# Main loop
# =========================
cat("\n=== Calculating ROGUE by scanvi_predictions ===\n")
res_list <- vector("list", nrow(cells_by_type))
pb <- txtProgressBar(max = nrow(cells_by_type), style = 3)

for (i in seq_len(nrow(cells_by_type))) {
  ct <- cells_by_type[[CELLTYPE_COL]][i]
  cells_keep <- intersect(unlist(cells_by_type$cells[i]), colnames(counts))

  expr_subset <- counts[, cells_keep, drop = FALSE]
  out <- calc_rogue_safe(expr_subset)

  res_list[[i]] <- data.frame(
    cell_type = ct,
    n_cells_input = length(cells_keep),
    n_cells_used = out$n_cells, # 若发生下采样会变化
    n_genes_used = out$n_genes,
    rogue_value = out$rogue,
    stringsAsFactors = FALSE
  )

  setTxtProgressBar(pb, i)
}
close(pb)

rogue_df <- bind_rows(res_list) %>%
  arrange(desc(rogue_value))

# =========================
# Save + quick report
# =========================
write.csv(
  rogue_df,
  file.path(OUTPUT_DIR, "rogue_by_scanvi_predictions.csv"),
  row.names = FALSE
)

n_failed <- sum(is.na(rogue_df$rogue_value))
cat(sprintf("\nFailed: %d/%d\n", n_failed, nrow(rogue_df)))
cat(sprintf(
  "Saved: %s\n",
  file.path(OUTPUT_DIR, "rogue_by_scanvi_predictions.csv")
))

# =========================
# Visualization: barplot (one value per cell type)
# =========================
pdf(
  file.path(OUTPUT_DIR, "rogue_barplot_by_celltype.pdf"),
  width = 12,
  height = 6
)
p <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  ggplot(aes(x = reorder(cell_type, rogue_value), y = rogue_value)) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed") +
  labs(
    title = "ROGUE by scANVI predicted cell type",
    x = "Cell type",
    y = "ROGUE",
    subtitle = sprintf(
      "MIN_CELLS=%d; MAX_CELLS=%s",
      MIN_CELLS_PER_GROUP,
      as.character(MAX_CELLS_PER_GROUP)
    )
  ) +
  theme_bw()
print(p)
dev.off()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ROGUE) # 关键：确保 matr.filter/SE_fun/CalculateRogue 可用
})

CELLTYPE_COL <- "scanvi_predictions"
ASSAY_USE <- "RNA"

MIN_CELLS_PER_GROUP <- 200
MAX_CELLS_PER_GROUP <- 10000
SEED <- 1

OUTPUT_DIR <- "rogue_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED)

# ---- get counts once ----
counts <- tryCatch(
  GetAssayData(seurat_obj, assay = ASSAY_USE, layer = "counts"),
  error = function(e) {
    GetAssayData(seurat_obj, assay = ASSAY_USE, slot = "counts")
  }
)

meta <- seurat_obj@meta.data
stopifnot(CELLTYPE_COL %in% colnames(meta))

# ---- align meta to counts by cell names (critical) ----
common_cells <- intersect(colnames(counts), rownames(meta))
if (length(common_cells) < 1000) {
  warning(sprintf(
    "Low overlap between meta and counts: %d cells. Check cell names.",
    length(common_cells)
  ))
}
counts <- counts[, common_cells, drop = FALSE]
meta <- meta[common_cells, , drop = FALSE]

# ---- valid cell types ----
grp_df <- meta %>%
  rownames_to_column("cell") %>%
  filter(!is.na(.data[[CELLTYPE_COL]]), .data[[CELLTYPE_COL]] != "") %>%
  count(.data[[CELLTYPE_COL]], name = "n_cells") %>%
  filter(n_cells >= MIN_CELLS_PER_GROUP) %>%
  arrange(desc(n_cells))

cat(sprintf(
  "Valid cell types (n >= %d): %d\n",
  MIN_CELLS_PER_GROUP,
  nrow(grp_df)
))
if (nrow(grp_df) == 0) {
  stop("No valid cell types after MIN_CELLS_PER_GROUP filter.")
}

cells_by_type <- meta %>%
  rownames_to_column("cell") %>%
  filter(.data[[CELLTYPE_COL]] %in% grp_df[[CELLTYPE_COL]]) %>%
  group_by(.data[[CELLTYPE_COL]]) %>%
  summarise(cells = list(cell), .groups = "drop")

# ---- safe calc with explicit error capture ----
calc_rogue_safe <- function(expr_sparse) {
  # return list(rogue, n_genes, n_cells, stage, msg)
  out <- list(
    rogue = NA_real_,
    n_genes = NA_integer_,
    n_cells = ncol(expr_sparse),
    stage = "init",
    msg = ""
  )

  if (ncol(expr_sparse) < MIN_CELLS_PER_GROUP) {
    out$stage <- "ncell_check"
    out$msg <- "ncol < MIN_CELLS_PER_GROUP"
    return(out)
  }

  # downsample
  if (
    is.finite(MAX_CELLS_PER_GROUP) && ncol(expr_sparse) > MAX_CELLS_PER_GROUP
  ) {
    keep <- sample(colnames(expr_sparse), MAX_CELLS_PER_GROUP)
    expr_sparse <- expr_sparse[, keep, drop = FALSE]
  }
  out$n_cells <- ncol(expr_sparse)

  # convert to dense after downsample
  expr_dense <- tryCatch(as.matrix(expr_sparse), error = function(e) NULL)
  if (is.null(expr_dense)) {
    out$stage <- "as.matrix"
    out$msg <- "as.matrix failed (memory?)"
    return(out)
  }

  # dynamic min.genes: 防止 HVG 子集把细胞全过滤掉
  # 至少 50，最多 200，且不超过基因数的 20%
  dyn_min_genes <- max(50, min(200, floor(nrow(expr_dense) * 0.2)))

  expr_f <- tryCatch(
    ROGUE::matr.filter(expr_dense, min.cells = 10, min.genes = dyn_min_genes),
    error = function(e) {
      out$stage <<- "matr.filter"
      out$msg <<- conditionMessage(e)
      NULL
    }
  )
  if (
    is.null(expr_f) || ncol(expr_f) < MIN_CELLS_PER_GROUP || nrow(expr_f) < 100
  ) {
    if (out$stage == "init") {
      out$stage <- "post_filter"
    }
    if (out$msg == "") {
      out$msg <- sprintf(
        "after filter: nrow=%s, ncol=%s",
        nrow(expr_f),
        ncol(expr_f)
      )
    }
    out$n_genes <- if (!is.null(expr_f)) nrow(expr_f) else NA_integer_
    return(out)
  }

  out$n_genes <- nrow(expr_f)

  ent_res <- tryCatch(
    ROGUE::SE_fun(expr_f),
    error = function(e) {
      out$stage <<- "SE_fun"
      out$msg <<- conditionMessage(e)
      NULL
    }
  )
  if (
    is.null(ent_res) ||
      any(is.na(ent_res$entropy)) ||
      any(is.infinite(ent_res$entropy))
  ) {
    if (out$stage == "init") {
      out$stage <- "entropy_invalid"
    }
    if (out$msg == "") {
      out$msg <- "entropy NA/Inf"
    }
    return(out)
  }

  rv <- tryCatch(
    ROGUE::CalculateRogue(ent_res, platform = "UMI"),
    error = function(e) {
      out$stage <<- "CalculateRogue"
      out$msg <<- conditionMessage(e)
      NA_real_
    }
  )

  out$rogue <- rv
  if (is.na(out$rogue) && out$stage == "init") {
    out$stage <- "rogue_na"
    out$msg <- "CalculateRogue returned NA"
  } else if (!is.na(out$rogue)) {
    out$stage <- "ok"
  }
  out
}

# ---- main loop ----
cat("\n=== Calculating ROGUE by scanvi_predictions ===\n")
res_list <- vector("list", nrow(cells_by_type))
pb <- txtProgressBar(max = nrow(cells_by_type), style = 3)

for (i in seq_len(nrow(cells_by_type))) {
  ct <- cells_by_type[[CELLTYPE_COL]][i]
  cells_keep <- intersect(unlist(cells_by_type$cells[i]), colnames(counts))

  expr_subset <- counts[, cells_keep, drop = FALSE]
  out <- calc_rogue_safe(expr_subset)

  res_list[[i]] <- data.frame(
    cell_type = ct,
    n_cells_input = length(cells_keep),
    n_cells_used = out$n_cells,
    n_genes_used = out$n_genes,
    rogue_value = out$rogue,
    fail_stage = out$stage,
    fail_msg = out$msg,
    stringsAsFactors = FALSE
  )

  setTxtProgressBar(pb, i)
}
close(pb)

rogue_df <- bind_rows(res_list)
write.csv(
  rogue_df,
  file.path(OUTPUT_DIR, "rogue_by_scanvi_predictions_debug.csv"),
  row.names = FALSE
)

cat("\nSaved debug table:\n")
cat(file.path(OUTPUT_DIR, "rogue_by_scanvi_predictions_debug.csv"), "\n")
print(rogue_df)

set.seed(1)

CELLTYPE_COL <- "scanvi_predictions"
ASSAY_USE <- "RNA"
N_SUB <- 200 # 等细胞数（<= 最小群体的 n_cells）
N_REP <- 30 # 重复次数

counts <- tryCatch(
  GetAssayData(seurat_obj, assay = ASSAY_USE, layer = "counts"),
  error = function(e) {
    GetAssayData(seurat_obj, assay = ASSAY_USE, slot = "counts")
  }
)

meta <- seurat_obj@meta.data
common_cells <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_cells, drop = FALSE]
meta <- meta[common_cells, , drop = FALSE]

celltypes <- names(which(table(meta[[CELLTYPE_COL]]) >= N_SUB))

boot_one <- function(ct) {
  cells <- rownames(meta)[meta[[CELLTYPE_COL]] == ct]
  vals <- numeric(N_REP)
  genes_used <- integer(N_REP)
  for (r in seq_len(N_REP)) {
    keep <- sample(cells, N_SUB)
    expr <- counts[, keep, drop = FALSE]
    expr <- as.matrix(expr)
    expr_f <- matr.filter(expr, min.cells = 10, min.genes = 50) # 这里 min.genes 放宽以减少“抽样导致掉光”
    ent <- SE_fun(expr_f)
    vals[r] <- CalculateRogue(ent, platform = "UMI")
    genes_used[r] <- nrow(expr_f)
  }
  data.frame(
    cell_type = ct,
    rogue_mean = mean(vals, na.rm = TRUE),
    rogue_sd = sd(vals, na.rm = TRUE),
    genes_used_mean = mean(genes_used, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

boot_df <- dplyr::bind_rows(lapply(celltypes, boot_one)) %>%
  dplyr::arrange(desc(rogue_mean))

write.csv(
  boot_df,
  file.path(OUTPUT_DIR, "rogue_bootstrap_equalN.csv"),
  row.names = FALSE
)
print(boot_df)

pl <- subset(seurat_obj, subset = scanvi_predictions == "Plasma cells")
stopifnot("mono_strict_high" %in% colnames(pl@meta.data))

counts_pl <- tryCatch(
  GetAssayData(pl, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(pl, assay = "RNA", slot = "counts")
)

calc_rogue <- function(mat) {
  x <- as.matrix(mat)
  x <- matr.filter(x, min.cells = 10, min.genes = 50)
  ent <- SE_fun(x)
  CalculateRogue(ent, platform = "UMI")
}

r_T <- calc_rogue(counts_pl[, colnames(pl)[pl$mono_strict_high], drop = FALSE])
r_F <- calc_rogue(counts_pl[, colnames(pl)[!pl$mono_strict_high], drop = FALSE])

cat("Plasma mono_strict_high TRUE :", r_T, "\n")
cat("Plasma mono_strict_high FALSE:", r_F, "\n")

# Seurat 自带的 cell cycle genes
cc <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes))

pgc <- subset(
  seurat_obj,
  subset = scanvi_predictions == "Proliferative germinal center B cells"
)
counts_pgc <- tryCatch(
  GetAssayData(pgc, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(pgc, assay = "RNA", slot = "counts")
)

genes_keep <- setdiff(rownames(counts_pgc), intersect(cc, rownames(counts_pgc)))

x <- as.matrix(counts_pgc[genes_keep, , drop = FALSE])
x <- matr.filter(x, min.cells = 10, min.genes = 50)
ent <- SE_fun(x)
cat(
  "Prolif GC (no cell-cycle genes) ROGUE:",
  CalculateRogue(ent, platform = "UMI"),
  "\n"
)

pl <- subset(seurat_obj, subset = scanvi_predictions == "Plasma cells")
Idents(pl) <- "mono_strict_high"

# 核心：看 top markers（TRUE vs FALSE）
m_pl <- FindMarkers(
  pl,
  ident.1 = "TRUE",
  ident.2 = "FALSE",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)
head(m_pl[order(m_pl$p_val_adj), ], 30)

pl <- subset(seurat_obj, subset = scanvi_predictions == "Plasma cells")
counts_pl <- tryCatch(
  GetAssayData(pl, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(pl, assay = "RNA", slot = "counts")
)

ig <- grep("^(IGH|IGK|IGL)", rownames(counts_pl), value = TRUE)
genes_keep <- setdiff(rownames(counts_pl), ig)

x <- as.matrix(counts_pl[genes_keep, , drop = FALSE])
x <- matr.filter(x, min.cells = 10, min.genes = 50)
ent <- SE_fun(x)
cat("Plasma (no IG genes) ROGUE:", CalculateRogue(ent, platform = "UMI"), "\n")

pgc <- subset(
  seurat_obj,
  subset = scanvi_predictions == "Proliferative germinal center B cells"
)

gc_core <- list(c("BCL6", "AICDA", "MEF2B", "RGS13", "S1PR2"))
plasma_core <- list(c("XBP1", "MZB1", "JCHAIN", "PRDM1", "SDC1"))
ifn_core <- list(c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1"))

pgc <- AddModuleScore(pgc, features = gc_core, name = "GCcore")
pgc <- AddModuleScore(pgc, features = plasma_core, name = "PLcore")
pgc <- AddModuleScore(pgc, features = ifn_core, name = "IFNcore")

# 看分布（是否双峰/长尾）
VlnPlot(pgc, features = c("GCcore1", "PLcore1", "IFNcore1"), pt.size = 0)

# 看 GC vs Plasma 是否互斥/混合
FeatureScatter(pgc, feature1 = "GCcore1", feature2 = "PLcore1")


pl <- subset(seurat_obj, subset = scanvi_predictions == "Plasma cells")
stopifnot("mono_strict_high" %in% colnames(pl@meta.data))

plasma_core <- list(c("JCHAIN", "MZB1", "XBP1", "SDC1", "PRDM1"))
epi_secretory <- list(c(
  "STATH",
  "BPIFA1",
  "BPIFB1",
  "SLPI",
  "LTF",
  "PIP",
  "ZG16B",
  "EPCAM",
  "KRT19"
))

pl <- AddModuleScore(pl, features = plasma_core, name = "PLcore")
pl <- AddModuleScore(pl, features = epi_secretory, name = "EPIsec")

# 1) 看两组的平均分（不用图）
print(aggregate(
  cbind(PLcore1, EPIsec1) ~ mono_strict_high,
  data = pl@meta.data,
  mean
))

# 2) 给一个“上皮污染”阈值：比如 EPIsec1 > 0（可按你数据调整）
pl$epi_like_flag <- pl$EPIsec1 > 0

cat("mono_strict_high TRUE fraction:", mean(pl$mono_strict_high), "\n")
cat("epi_like_flag fraction:", mean(pl$epi_like_flag), "\n")
cat(
  "overlap (TRUE & epi_like):",
  mean(pl$mono_strict_high & pl$epi_like_flag),
  "\n"
)


pl2 <- subset(pl, subset = mono_strict_high == FALSE)

DefaultAssay(pl2) <- "RNA"
pl2 <- NormalizeData(pl2)
pl2 <- FindVariableFeatures(pl2, nfeatures = 3000)
pl2 <- ScaleData(pl2, features = VariableFeatures(pl2))
pl2 <- RunPCA(pl2, npcs = 30)
pl2 <- FindNeighbors(pl2, dims = 1:30)
pl2 <- FindClusters(pl2, resolution = 0.4)
pl2 <- RunUMAP(pl2, dims = 1:30)

# 每个子群 markers
m <- FindAllMarkers(
  pl2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  m,
  file.path(OUTPUT_DIR, "plasma_FALSE_subcluster_markers.csv"),
  row.names = FALSE
)

# 也可以直接看几个关键轴的分布（不用画图也行）
genes <- intersect(
  c(
    "JCHAIN",
    "MZB1",
    "XBP1",
    "SDC1",
    "PRDM1",
    "HLA-DRA",
    "ISG15",
    "IFI6",
    "XIST",
    "MKI67"
  ),
  rownames(pl2)
)
avg <- AverageExpression(
  pl2,
  features = genes,
  group.by = "seurat_clusters"
)$RNA
print(round(avg, 3))

pgc <- subset(
  seurat_obj,
  subset = scanvi_predictions == "Proliferative germinal center B cells"
)

ifn <- list(c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1"))
mhc <- list(c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1"))
stress <- list(c("HSPA1A", "HSPA1B", "DNAJB1", "FOS", "JUN"))

pgc <- AddModuleScore(pgc, ifn, name = "IFN")
pgc <- AddModuleScore(pgc, mhc, name = "MHC2")
pgc <- AddModuleScore(pgc, stress, name = "STRESS")

print(summary(pgc$IFN1))
print(summary(pgc$MHC21))
print(summary(pgc$STRESS1))

# 内部重聚类
pgc <- NormalizeData(pgc)
pgc <- FindVariableFeatures(pgc, nfeatures = 3000)
pgc <- ScaleData(pgc, features = VariableFeatures(pgc))
pgc <- RunPCA(pgc, npcs = 30)
pgc <- FindNeighbors(pgc, dims = 1:30)
pgc <- FindClusters(pgc, resolution = 0.4)

# 输出每个子群细胞数，先确认是否存在“混合的多个小群”
print(sort(table(pgc$seurat_clusters), decreasing = TRUE))

pl <- subset(seurat_obj, subset = scanvi_predictions == "Plasma cells")

epi_struct <- list(c(
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "TACSTD2",
  "CLDN1",
  "CLDN4"
))
epi_secr <- list(c("STATH", "BPIFA1", "BPIFB1", "SLPI", "LTF", "PIP", "ZG16B"))

pl <- AddModuleScore(pl, epi_struct, name = "EPIstruct")
pl <- AddModuleScore(pl, epi_secr, name = "EPIsecr")
pl <- AddModuleScore(
  pl,
  list(c("JCHAIN", "MZB1", "XBP1", "PRDM1", "SDC1")),
  name = "PLcore"
)

# 用更严格的阈值：例如 top 5% 作为“强阳性”
q95_struct <- quantile(pl$EPIstruct1, 0.95, na.rm = TRUE)
q95_secr <- quantile(pl$EPIsecr1, 0.95, na.rm = TRUE)

pl$epi_struct_hi <- pl$EPIstruct1 >= q95_struct
pl$epi_secr_hi <- pl$EPIsecr1 >= q95_secr

# 关键：区分 doublet-like vs ambient-like
pl$epi_doublet_like <- pl$epi_struct_hi # 结构基因高：更像双细胞/误标
pl$epi_ambient_like <- pl$epi_secr_hi & !pl$epi_struct_hi # 仅分泌高：更像环境RNA

print(prop.table(table(pl$epi_doublet_like)))
print(prop.table(table(pl$epi_ambient_like)))
print(prop.table(table(pl$epi_doublet_like, pl$mono_strict_high)))

# 如果你要“清理 plasma”用于下游：
pl_clean <- subset(pl, subset = !epi_doublet_like)

# pgc <- subset(seurat_obj, subset = scanvi_predictions == "Proliferative germinal center B cells")
DefaultAssay(pgc) <- "RNA"

gc_core <- list(c("BCL6", "AICDA", "MEF2B", "RGS13", "S1PR2"))
plasma_core <- list(c("XBP1", "MZB1", "JCHAIN", "PRDM1", "SDC1"))
mhc2_core <- list(c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1"))
ifn_core <- list(c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1"))
cc_core <- list(c("MKI67", "TOP2A", "HMGB2", "TUBB", "TYMS"))

pgc <- AddModuleScore(pgc, gc_core, name = "GCcore")
pgc <- AddModuleScore(pgc, plasma_core, name = "PLcore")
pgc <- AddModuleScore(pgc, mhc2_core, name = "MHC2")
pgc <- AddModuleScore(pgc, ifn_core, name = "IFN")
pgc <- AddModuleScore(pgc, cc_core, name = "CC")

# 你当前的聚类列名如果是 seurat_clusters 就用它，否则替换
tab <- pgc@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    n = dplyr::n(),
    GCcore = mean(GCcore1),
    PLcore = mean(PLcore1),
    MHC2 = mean(MHC21),
    IFN = mean(IFN1),
    CC = mean(CC1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n))

print(tab)

pgc <- subset(
  seurat_obj,
  subset = scanvi_predictions == "Proliferative germinal center B cells"
)
DefaultAssay(pgc) <- "RNA"
Idents(pgc) <- "seurat_clusters" # 你现在这张表对应的就是这个

genes <- c(
  # B-lineage core
  "MS4A1",
  "CD79A",
  "CD79B",
  "CD19",
  "PAX5",
  "CD37",
  "CD22",
  # GC program
  "BCL6",
  "AICDA",
  "RGS13",
  "S1PR2",
  "MEF2B",
  # Plasma program
  "XBP1",
  "MZB1",
  "JCHAIN",
  "PRDM1",
  "SDC1",
  # MHC-II / APC
  "CD74",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DPA1",
  "HLA-DPB1",
  # Cycling
  "MKI67",
  "TOP2A",
  "TYMS",
  # pDC
  "IL3RA",
  "TCF4",
  "GZMB",
  "CLEC4C",
  "LILRA4",
  # Myeloid
  "LST1",
  "TYROBP",
  "FCER1G",
  "LYZ",
  "S100A8",
  "S100A9",
  # Epithelial (doublet)
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19"
)
genes <- intersect(genes, rownames(pgc))

avg <- AverageExpression(
  pgc,
  features = genes,
  group.by = "seurat_clusters"
)$RNA
print(round(avg, 3))

dp <- DotPlot(pgc, features = genes) + RotatedAxis()
print(dp)

# pgc <- subset(seurat_obj, subset = scanvi_predictions == "Proliferative germinal center B cells")
Idents(pgc) <- "seurat_clusters"

# 1) 先给 subcluster 贴建议标签
pgc$pgc_refined <- as.character(pgc$seurat_clusters)
pgc$pgc_refined[pgc$seurat_clusters %in% c("0", "1")] <- "Cycling_GC_like_B"
pgc$pgc_refined[
  pgc$seurat_clusters %in% c("2")
] <- "MHC2hi_GC_like_B_noncycling"
pgc$pgc_refined[pgc$seurat_clusters %in% c("3")] <- "GC_like_low_MHC2_mixed"
pgc$pgc_refined[
  pgc$seurat_clusters %in% c("4")
] <- "Plasma_like_epithelial_doublet"
pgc$pgc_refined <- factor(pgc$pgc_refined)

table(pgc$pgc_refined)

# 2) 生成 clean 版本：剔除 g4（必要时也可连 g3 一起剔除做敏感性分析）
pgc_clean_main <- subset(
  pgc,
  subset = pgc_refined != "Plasma_like_epithelial_doublet"
)
pgc_clean_sens <- subset(
  pgc,
  subset = pgc_refined %in%
    c("Cycling_GC_like_B", "MHC2hi_GC_like_B_noncycling")
)

# 3)（推荐）比较三套的 ROGUE：原始 / main-clean / sens-clean
calc_rogue_obj <- function(obj) {
  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
  x <- as.matrix(counts)
  x <- matr.filter(x, min.cells = 10, min.genes = 50)
  ent <- SE_fun(x)
  CalculateRogue(ent, platform = "UMI")
}

cat("ProlifGC ROGUE (raw)      :", calc_rogue_obj(pgc), "\n")
cat("ProlifGC ROGUE (main)     :", calc_rogue_obj(pgc_clean_main), "\n")
cat("ProlifGC ROGUE (sens)     :", calc_rogue_obj(pgc_clean_sens), "\n")

qc_tab <- pgc@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    n = dplyr::n(),
    nCount_RNA = median(nCount_RNA),
    nFeature_RNA = median(nFeature_RNA),
    percent.mt = median(percent.mt),
    .groups = "drop"
  )
print(qc_tab)

# pgc：你已经做好的 Prolif GC 子集，并含 pgc_refined
# pgc <- subset(seurat_obj, subset = scanvi_predictions == "Proliferative germinal center B cells")

# 确保 pgc_refined 存在
stopifnot("pgc_refined" %in% colnames(pgc@meta.data))

# ---- 1) 回填到全对象：scanvi_refined ----
seurat_obj$scanvi_refined <- as.character(seurat_obj$scanvi_predictions)
seurat_obj$scanvi_refined[colnames(pgc)] <- as.character(pgc$pgc_refined)
seurat_obj$scanvi_refined <- factor(seurat_obj$scanvi_refined)

# ---- 2) 打 flag，方便后续统一过滤 ----
seurat_obj$flag_doublet_pgc_epi <- FALSE
seurat_obj$flag_doublet_pgc_epi[colnames(pgc)[
  pgc$pgc_refined == "Plasma_like_epithelial_doublet"
]] <- TRUE

seurat_obj$flag_lowqc_pgc <- FALSE
seurat_obj$flag_lowqc_pgc[colnames(pgc)[
  pgc$pgc_refined == "GC_like_low_MHC2_mixed"
]] <- TRUE

# ---- 3) 三套集合 ----
cells_pgc_cycling <- colnames(pgc)[pgc$pgc_refined == "Cycling_GC_like_B"]

cells_pgc_extended <- colnames(pgc)[
  pgc$pgc_refined %in%
    c("Cycling_GC_like_B", "MHC2hi_GC_like_B_noncycling")
]

cells_pgc_all_clean <- colnames(pgc)[
  pgc$pgc_refined != "Plasma_like_epithelial_doublet"
]

# 例：生成 clean 子对象
pgc_cycling <- subset(seurat_obj, cells = cells_pgc_cycling)
pgc_extended <- subset(seurat_obj, cells = cells_pgc_extended)
pgc_allclean <- subset(seurat_obj, cells = cells_pgc_all_clean)

table(seurat_obj$scanvi_refined)


# 1) 看样本来源
tab_src <- table(
  seurat_obj$sample[seurat_obj$flag_doublet_pgc_epi],
  useNA = "ifany"
)
print(sort(tab_src, decreasing = TRUE))

# 2) 看它们在 scanvi_predictions 里的原始标签分布（以防有其它误归类）
print(table(seurat_obj$scanvi_predictions[seurat_obj$flag_doublet_pgc_epi]))

# 3) 如果你有 UMAP：
# FeaturePlot(seurat_obj, features = "flag_doublet_pgc_epi")
# 或 DimPlot(seurat_obj, group.by = "scanvi_refined")

library(dplyr)

meta <- seurat_obj@meta.data

# 你的 B 细胞定义按需调整（这里包含你目前所有 B 相关标签）
b_levels <- c(
  "Age-associated B cells",
  "Germinal center B cells",
  "Memory B cells",
  "Naive B cells",
  "Plasma cells",
  "Cycling_GC_like_B",
  "MHC2hi_GC_like_B_noncycling",
  "GC_like_low_MHC2_mixed"
)

meta$is_B <- meta$scanvi_refined %in%
  b_levels |
  meta$scanvi_predictions %in% c("Proliferative germinal center B cells")
meta$is_doublet_known <- meta$scanvi_refined == "Plasma_like_epithelial_doublet"

rate_by_sample <- meta %>%
  filter(is_B) %>%
  group_by(sample) %>%
  summarise(
    n_B = n(),
    n_doublet = sum(is_doublet_known),
    frac_doublet = n_doublet / n_B,
    .groups = "drop"
  ) %>%
  arrange(desc(frac_doublet), desc(n_doublet))

print(rate_by_sample, n = 50)

DefaultAssay(seurat_obj) <- "RNA"

gs <- list(
  PlasmaCore = c("JCHAIN", "MZB1", "XBP1", "PRDM1", "SDC1"),
  EPIstruct = c("EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2", "CLDN4"),
  Bcore = c("MS4A1", "CD79A", "CD79B", "CD19", "PAX5", "CD37", "CD22")
)

for (nm in names(gs)) {
  genes <- intersect(gs[[nm]], rownames(seurat_obj))
  seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), name = nm)
}

pos <- function(x) pmax(0, x)
z_by_group <- function(x, g) ave(x, g, FUN = function(v) as.numeric(scale(v)))

meta <- seurat_obj@meta.data
group_col <- "sample" # 或 dataset；建议用 sample

meta$z_plasma <- z_by_group(meta$PlasmaCore1, meta[[group_col]])
meta$z_epi <- z_by_group(meta$EPIstruct1, meta[[group_col]])
meta$z_bcore <- z_by_group(meta$Bcore1, meta[[group_col]])

meta$log10_nCount <- log10(meta$nCount_RNA + 1)
meta$log10_nFeature <- log10(meta$nFeature_RNA + 1)
meta$z_nCount <- z_by_group(meta$log10_nCount, meta[[group_col]])
meta$z_nFeature <- z_by_group(meta$log10_nFeature, meta[[group_col]])

# 关键：plasma×epi 的 co-high
meta$score_plasma_epi <- pos(meta$z_plasma) * pos(meta$z_epi)

# 可选：B×epi（避免把某些非 plasma 的 B×epi doublet 漏掉）
meta$score_b_epi <- pos(meta$z_bcore) * pos(meta$z_epi)

meta$cross_score <- pmax(meta$score_plasma_epi, meta$score_b_epi)

w_qc <- 0.3
meta$qc_score <- pos(meta$z_nCount) + pos(meta$z_nFeature)

meta$doublet_score <- meta$cross_score + w_qc * meta$qc_score

seurat_obj@meta.data <- meta

meta <- seurat_obj@meta.data
meta$thr_sample <- ave(
  meta$doublet_score,
  meta[[group_col]],
  FUN = function(v) quantile(v, 0.99, na.rm = TRUE)
)

# 最终阈值：取更严格者
meta$doublet_flag_adaptive <- meta$doublet_score >= pmax(thr, meta$thr_sample)

seurat_obj@meta.data <- meta

meta <- seurat_obj@meta.data
top_cells <- rownames(meta)[order(meta$doublet_score, decreasing = TRUE)][1:200]
sub <- subset(seurat_obj, cells = top_cells)

pp <- DotPlot(
  sub,
  features = intersect(
    c(
      "EPCAM",
      "KRT19",
      "KRT8",
      "KRT18",
      "JCHAIN",
      "MZB1",
      "XBP1",
      "PRDM1",
      "MS4A1",
      "CD79A"
    ),
    rownames(sub)
  )
) +
  RotatedAxis()

pp$data

rate_by_sample2 <- rate_by_sample %>%
  mutate(flag_smallN = n_B < 50) %>%
  arrange(desc(frac_doublet))

a <- 1
b <- 99
rate_by_sample_shrunk <- rate_by_sample %>%
  mutate(frac_shrunk = (n_doublet + a) / (n_B + a + b)) %>%
  arrange(desc(frac_shrunk))

library(Matrix)

DefaultAssay(seurat_obj) <- "RNA"
counts <- tryCatch(
  GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
)

plasma_genes <- intersect(
  c("JCHAIN", "MZB1", "XBP1", "PRDM1", "SDC1"),
  rownames(counts)
)
epi_genes <- intersect(
  c("EPCAM", "TACSTD2", "CLDN4", "KRT8", "KRT18"),
  rownames(counts)
)

# counts>=2 的“共检出强度”
pl_n2 <- Matrix::colSums(counts[plasma_genes, , drop = FALSE] >= 2)
epi_n2 <- Matrix::colSums(counts[epi_genes, , drop = FALSE] >= 2)

meta <- seurat_obj@meta.data
meta$pl_n2 <- pl_n2[colnames(seurat_obj)]
meta$epi_n2 <- epi_n2[colnames(seurat_obj)]
meta$cross_count_score <- meta$pl_n2 * meta$epi_n2

# QC（样本内z）
group_col <- "sample"
z_by_group <- function(x, g) ave(x, g, FUN = function(v) as.numeric(scale(v)))
pos <- function(x) pmax(0, x)

meta$log10_nCount <- log10(meta$nCount_RNA + 1)
meta$log10_nFeature <- log10(meta$nFeature_RNA + 1)
meta$z_nCount <- z_by_group(meta$log10_nCount, meta[[group_col]])
meta$z_nFeature <- z_by_group(meta$log10_nFeature, meta[[group_col]])
meta$qc_score <- pos(meta$z_nCount) + pos(meta$z_nFeature)

w_qc <- 0.3
meta$doublet_score <- meta$cross_count_score + w_qc * meta$qc_score

seurat_obj@meta.data <- meta

pos_idx <- which(seurat_obj$scanvi_refined == "Plasma_like_epithelial_doublet")

thr <- quantile(seurat_obj$doublet_score[pos_idx], 0.10, na.rm = TRUE) # 召回约90%
seurat_obj$doublet_flag <- seurat_obj$doublet_score >= thr

cat("Global flagged rate:", mean(seurat_obj$doublet_flag, na.rm = TRUE), "\n")
cat(
  "Recall on known positives:",
  mean(seurat_obj$doublet_flag[pos_idx], na.rm = TRUE),
  "\n"
)

meta <- seurat_obj@meta.data
k <- 0.02

meta$thr_sample <- ave(
  meta$doublet_score,
  meta[[group_col]],
  FUN = function(v) quantile(v, 1 - k, na.rm = TRUE)
)

# 最终判定：同时满足全局阈值 AND 不超过样本上限（更保守）
meta$doublet_flag2 <- meta$doublet_score >= pmax(thr, meta$thr_sample)

seurat_obj@meta.data <- meta
summary(seurat_obj$cross_count_score[pos_idx])


library(Matrix)
DefaultAssay(seurat_obj) <- "RNA"

counts <- tryCatch(
  GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
)

plasma_genes <- intersect(
  c("JCHAIN", "MZB1", "XBP1", "PRDM1", "SDC1"),
  rownames(counts)
)
epi_anchor <- intersect(
  c("EPCAM", "TACSTD2", "CLDN4", "CLDN7"),
  rownames(counts)
) # 有哪些用哪些
epi_pair <- intersect(c("KRT8", "KRT18"), rownames(counts))
krt19_gene <- intersect("KRT19", rownames(counts))

meta <- seurat_obj@meta.data

# Plasma 强证据：counts>=2 的 plasma gene 个数
meta$pl_n2 <- Matrix::colSums(counts[plasma_genes, , drop = FALSE] >= 2)

# EPI 证据：counts>=1
epi_anchor_hit <- if (length(epi_anchor) > 0) {
  Matrix::colSums(counts[epi_anchor, , drop = FALSE] >= 1)
} else {
  rep(0, ncol(counts))
}
krt8_hit <- if ("KRT8" %in% rownames(counts)) {
  as.numeric(counts["KRT8", ] >= 1)
} else {
  rep(0, ncol(counts))
}
krt18_hit <- if ("KRT18" %in% rownames(counts)) {
  as.numeric(counts["KRT18", ] >= 1)
} else {
  rep(0, ncol(counts))
}
pair_hit <- krt8_hit * krt18_hit

# KRT19 仅作“弱证据”，且要求更高计数以压 ambient（比如>=5；你可调）
krt19_strong <- if (length(krt19_gene) > 0) {
  as.numeric(counts["KRT19", ] >= 5)
} else {
  rep(0, ncol(counts))
}

meta$epi_anchor_hit_n1 <- epi_anchor_hit
meta$epi_pair_hit <- pair_hit
meta$krt19_strong <- krt19_strong

# 合成 EPI “是否支持”：anchor>0 或 pair_hit==1 或 krt19_strong==1
meta$epi_support <- (meta$epi_anchor_hit_n1 > 0) |
  (meta$epi_pair_hit == 1) |
  (meta$krt19_strong == 1)

# 跨谱系证据（主信号）：pl_n2 * epi_support（epi_support 是 0/1）
meta$cross_support_score <- meta$pl_n2 * as.numeric(meta$epi_support)

# QC 仅加分
group_col <- "sample"
z_by_group <- function(x, g) ave(x, g, FUN = function(v) as.numeric(scale(v)))
pos <- function(x) pmax(0, x)

meta$log10_nCount <- log10(meta$nCount_RNA + 1)
meta$log10_nFeature <- log10(meta$nFeature_RNA + 1)
meta$z_nCount <- z_by_group(meta$log10_nCount, meta[[group_col]])
meta$z_nFeature <- z_by_group(meta$log10_nFeature, meta[[group_col]])
meta$qc_score <- pos(meta$z_nCount) + pos(meta$z_nFeature)

# 最终分数：先靠 cross_support，再少量加 QC
w_qc <- 0.1
meta$doublet_score2 <- meta$cross_support_score + w_qc * meta$qc_score

seurat_obj@meta.data <- meta
table(meta$pl_n2[pos_idx])

table(meta$epi_support[pos_idx])

library(Matrix)
meta <- seurat_obj@meta.data

counts <- tryCatch(
  GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
)

pos_idx <- which(meta$scanvi_refined == "Plasma_like_epithelial_doublet")

plasma_genes <- intersect(
  c("JCHAIN", "MZB1", "XBP1", "PRDM1", "SDC1"),
  rownames(counts)
)
epi_genes <- intersect(
  c(
    "EPCAM",
    "TACSTD2",
    "CLDN4",
    "CLDN7",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "KRT17"
  ),
  rownames(counts)
)
anchor_genes <- intersect(
  c("EPCAM", "TACSTD2", "CLDN4", "CLDN7"),
  rownames(counts)
)

pl_n1 <- Matrix::colSums(counts[plasma_genes, , drop = FALSE] >= 1)
pl_n2 <- Matrix::colSums(counts[plasma_genes, , drop = FALSE] >= 2)

epi_n1 <- Matrix::colSums(counts[epi_genes, , drop = FALSE] >= 1)
anchor_n1 <- if (length(anchor_genes) > 0) {
  Matrix::colSums(counts[anchor_genes, , drop = FALSE] >= 1)
} else {
  rep(0, ncol(counts))
}
krt19 <- if ("KRT19" %in% rownames(counts)) {
  as.numeric(counts["KRT19", ])
} else {
  rep(0, ncol(counts))
}

# 网格：你可以按需要缩放
grid <- expand.grid(
  pl_thr = c(1, 2), # 用 pl_n2 >= pl_thr
  epi_thr = c(2, 3, 4), # epi_n1 >= epi_thr（关键：至少2个epi基因共检出）
  krt19_thr = c(0, 2, 5), # krt19>=krt19_thr 作为“弱加成”而非必需
  need_anchor = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

eval_one <- function(pl_thr, epi_thr, krt19_thr, need_anchor) {
  # EPI 支持：epi_n1 达标 且（anchor可选要求） 且（KRT19可选加成）
  epi_ok <- (epi_n1 >= epi_thr) &
    (if (need_anchor) (anchor_n1 > 0) else TRUE) &
    (if (krt19_thr > 0) (krt19 >= krt19_thr | anchor_n1 > 0) else TRUE)

  cand <- (pl_n2 >= pl_thr) & epi_ok

  data.frame(
    pl_thr = pl_thr,
    epi_thr = epi_thr,
    krt19_thr = krt19_thr,
    need_anchor = need_anchor,
    global_rate = mean(cand),
    recall_pos = mean(cand[pos_idx])
  )
}

res <- do.call(
  rbind,
  Map(eval_one, grid$pl_thr, grid$epi_thr, grid$krt19_thr, grid$need_anchor)
)
res <- res[order(-res$recall_pos, res$global_rate), ]
print(head(res, 20))
