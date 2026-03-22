#!/usr/bin/env Rscript
# =============================================================================
# RNA_snn_res.2 Cluster QC + Lineage Validation (Epithelial focus) v1.0
# Output:
#  - cluster_qc_summary.csv
#  - cluster_signature_summary.csv
#  - cluster_problem_report.txt
#  - qc_violin.pdf
#  - marker_dotplot_RNA_snn_res.2.pdf
#  - signature_umap.pdf
# =============================================================================

# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
# output_dir <- "/home/h2048/data/R/0103/ciliated"
# dir.create(output_dir,recursive=TRUE)
# setwd()
library(reticulate)
library(data.table)
library(harmony)
library(ggplot2)
library(patchwork)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct

# ===== Remove Ribosomal, Mitochondrial, ENSG, and Unannotated Genes =====

# seurat_obj <- GetSeurat(
#   '/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2/Ciliated/batch_aware/cnmf_analysis_k40/Ciliated_with_cnmf_k40.h5ad'
# # )
seurat_obj <- readRDS(
  '/home/h2048/data/R/0103/ciliated/ciliated_filtered_20260103.rds'
)
output_dir <- '/home/h2048/data/R/0103/ciliated'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
library(Seurat)
library(ggplot2)

# ===== 定义关键marker基因（按功能分组）=====
markers_ciliated <- list(
  # 核心纤毛标志
  "Core Ciliated" = c("FOXJ1", "RSPH1", "DNAH5", "DNAI1", "PIFO"),

  # 纤毛发生/增殖
  "Ciliogenesis" = c("CCNO", "DEUP1", "CDC20B", "FOXN4", "MKI67", "TOP2A"),

  # 成熟纤毛亚型
  "Mature Ciliated" = c("DNAH11", "CFAP46", "CFAP54", "SPATA6L", "RSPH10B2"),

  # 干扰素响应
  "IFN Response" = c("IFITM1", "ISG15", "IFI44L", "MX1", "MX2", "OAS2"),

  # 分泌/杯状特征
  "Secretory-like" = c("MUC5AC", "SCGB1A1", "SCGB3A1", "BPIFB1", "VMO1"),

  # HLA/抗原递呈
  "HLA-high" = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "DUOX2"),

  # 炎症激活
  "Inflammatory" = c("FOS", "JUN", "IL8", "CCL20", "CXCL10", "SAA1"),

  # 代谢相关
  "Metabolic" = c("HARS", "KARS", "NARS", "ATP5F1C", "GAPDH")
)

# ===== DotPlot（经典样式）=====
DotPlot(
  seurat_obj, # 你的纤毛细胞Seurat对象
  features = unlist(markers_ciliated),
  group.by = "seurat_clusters",
  cols = c("lightgrey", "red"),
  dot.scale = 6
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(title = "Ciliated Cell Subclusters - Key Markers")

ggsave("ciliated_dotplot_all_markers.pdf", width = 18, height = 8)


# ===== DotPlot（分面展示 - 更清晰）=====
# 将marker转为数据框
marker_df <- stack(markers_ciliated)
colnames(marker_df) <- c("gene", "category")

# 按类别分别绘制
for (cat in names(markers_ciliated)) {
  genes <- markers_ciliated[[cat]]

  p <- DotPlot(
    seurat_obj,
    features = genes,
    cols = c("lightblue", "darkred"),
    dot.scale = 8
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = paste0(cat, " Markers"))

  print(p)
  ggsave(
    paste0("ciliated_dotplot_", gsub(" ", "_", cat), ".pdf"),
    width = 8,
    height = 6
  )
}


library(Seurat)
library(ggplot2)
library(dplyr)

# ===== 1. 快速定义关键亚群的top marker =====
markers_by_function <- list(
  # 纤毛发生/增殖 (Cluster 10)
  "Ciliogenesis_C10" = c(
    "CCNO",
    "CDC20B",
    "FOXN4",
    "MCIDAS",
    "EPIST",
    "NEK2",
    "POC1A",
    "HYLS1"
  ),

  # 炎症应答 (Cluster 11, 13)
  "Inflammatory_C11" = c("CYR61", "IL8", "JUNB", "FOS", "ATF3", "NR4A1"),
  "IFN_Response_C13" = c(
    "CCL20",
    "OASL",
    "IFIT2",
    "IFIT3",
    "ISG15",
    "IFITM1",
    "MX2",
    "OAS3"
  ),

  # HLA高表达/补体 (Cluster 14)
  "HLA_Complement_C14" = c(
    "C4A",
    "C4B",
    "C3",
    "SERPING1",
    "HLA-DRA",
    "HLA-DRB1",
    "HLA-DQA1",
    "CD74"
  ),

  # 分泌/杯状混合 (Cluster 1)
  "Goblet_Mixed_C1" = c(
    "MUC5AC",
    "KRT5",
    "SERPINB3",
    "SCGB1A1",
    "FABP5",
    "EPAS1"
  ),

  # 核心纤毛标志（跨cluster验证）
  "Core_Ciliated" = c("FOXJ1", "RSPH1", "DNAH5", "DNAI1", "PIFO"),

  # 代谢型 (Cluster 7)
  "Metabolic_C7" = c("ALPL", "BCAT1", "SERPINB2", "TM4SF1", "NMU"),

  # 发育相关 (Cluster 4, 8, 15)
  "Developmental" = c("ESR1", "ZBTB16", "LRP1B", "LDB2", "RBMS3")
)

# ===== 2. DotPlot - 按亚群分组 =====
# 合并所有marker
all_markers <- unique(unlist(markers_by_function))

# 基础DotPlot
DotPlot(
  seurat_obj,
  features = all_markers,
  cols = c("lightgrey", "red"),
  dot.scale = 5,
  cluster.idents = TRUE # 自动排序cluster
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 10)
  ) +
  labs(title = "Ciliated Subclusters - Functional Markers")

ggsave("ciliated_functional_dotplot.pdf", width = 16, height = 10)


# ===== 3. 分类别DotPlot（更清晰）=====
for (category in names(markers_by_function)) {
  genes <- markers_by_function[[category]]

  # 过滤掉可能不存在的基因
  genes_present <- genes[genes %in% rownames(seurat_obj)]

  if (length(genes_present) == 0) {
    next
  }

  p <- DotPlot(
    seurat_obj,
    features = genes_present,
    cols = c("lightblue", "darkred"),
    dot.scale = 8
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 11,
        face = "bold"
      ),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(title = gsub("_", " ", category))

  print(p)
  ggsave(paste0("dotplot_", category, ".pdf"), width = 8, height = 6)
}


# ===== 4. FeaturePlot - 关键marker空间分布 =====
# 4.1 纤毛发生 vs 成熟
FeaturePlot(
  seurat_obj,
  features = c(
    "CCNO",
    "CDC20B",
    "FOXN4",
    "MKI67", # 纤毛发生
    "FOXJ1",
    "RSPH1",
    "DNAH5",
    "DNAI1"
  ), # 成熟纤毛
  ncol = 4,
  pt.size = 0.5,
  order = TRUE,
  cols = c("grey90", "navy")
) &
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

ggsave("featureplot_ciliogenesis_vs_mature.pdf", width = 16, height = 8)


# 4.2 炎症相关
FeaturePlot(
  seurat_obj,
  features = c(
    "IL8",
    "CCL20",
    "CXCL10", # 趋化因子
    "IFITM1",
    "ISG15",
    "MX2"
  ), # 干扰素
  ncol = 3,
  pt.size = 0.4,
  cols = c("grey90", "red")
)

ggsave("featureplot_inflammatory.pdf", width = 12, height = 8)


# 4.3 HLA/补体
FeaturePlot(
  seurat_obj,
  features = c("HLA-DRA", "HLA-DRB1", "CD74", "C3", "C4A", "SERPING1"),
  ncol = 3,
  cols = c("grey90", "purple")
)

ggsave("featureplot_HLA_complement.pdf", width = 12, height = 8)


# ===== 5. 从你的marker表格提取Top5可视化 =====
# 定义每个cluster的Top marker（基于你的表格）
cluster_top_markers <- list(
  "1" = c("MUC5AC", "KRT5", "EPAS1", "SERPINB3", "FABP5"),
  "7" = c("ALPL", "BCAT1", "SERPINB2", "TM4SF1", "NMU"),
  "10" = c("CCNO", "CDC20B", "FOXN4", "NEK2", "EPIST"),
  "11" = c("IL8", "JUNB", "FOS", "ATF3", "NR4A1"),
  "13" = c("CCL20", "OASL", "IFIT2", "IFIT3", "ISG15"),
  "14" = c("C4A", "C4B", "C3", "HLA-DRA", "HLA-DRB1")
)

# 批量生成FeaturePlot
for (cluster_id in names(cluster_top_markers)) {
  genes <- cluster_top_markers[[cluster_id]]
  genes_present <- genes[genes %in% rownames(seurat_obj)]

  if (length(genes_present) < 4) {
    next
  }

  p <- FeaturePlot(
    seurat_obj,
    features = genes_present[1:min(6, length(genes_present))],
    ncol = 3,
    pt.size = 0.3,
    order = TRUE
  ) &
    theme(plot.title = element_text(size = 11, face = "bold"))

  print(p)
  ggsave(
    paste0("featureplot_cluster", cluster_id, "_top_markers.pdf"),
    width = 12,
    height = 8
  )
}


# ===== 6. VlnPlot - 查看表达分布 =====
# 关键功能marker的violin plot
key_features <- c(
  "CCNO", # 纤毛发生
  "FOXJ1", # 成熟纤毛
  "MKI67", # 增殖
  "IL8", # 炎症
  "IFITM1", # 干扰素
  "HLA-DRA", # 抗原递呈
  "MUC5AC", # 杯状
  "C3" # 补体
)

VlnPlot(
  seurat_obj,
  features = key_features,
  ncol = 4,
  pt.size = 0,
  cols = rainbow(length(unique(Idents(seurat_obj))))
) &
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 11, face = "bold")
  )

ggsave("vlnplot_key_features.pdf", width = 16, height = 8)


# ===== 7. DotPlot - 聚焦重点cluster对比 =====
# 对比几个重要cluster
Idents(seurat_obj) <- "seurat_clusters"
key_clusters <- c("1", "7", "10", "11", "13", "14")

if (all(key_clusters %in% levels(Idents(seurat_obj)))) {
  subset_obj <- subset(seurat_obj, idents = key_clusters)

  comparison_markers <- c(
    # 纤毛
    "FOXJ1",
    "CCNO",
    "CDC20B",
    # 炎症
    "IL8",
    "CCL20",
    "IFITM1",
    "ISG15",
    # HLA
    "HLA-DRA",
    "CD74",
    "C3",
    # 分泌
    "MUC5AC",
    "SCGB1A1",
    # 代谢
    "ALPL",
    "BCAT1"
  )

  DotPlot(
    subset_obj,
    features = comparison_markers,
    cols = c("lightgrey", "red"),
    dot.scale = 8
  ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        size = 11,
        face = "bold"
      ),
      axis.text.y = element_text(size = 11)
    ) +
    labs(title = "Key Ciliated Subclusters Comparison")

  ggsave("dotplot_key_clusters_comparison.pdf", width = 10, height = 7)
}


# ===== 8. 热图展示Top marker =====
# 如果有FindAllMarkers的结果
if (exists("ciliated_markers")) {
  top10 <- ciliated_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

  DoHeatmap(
    seurat_obj,
    features = top10$gene,
    group.by = "seurat_clusters",
    size = 3,
    angle = 90
  ) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    theme(axis.text.y = element_text(size = 6))

  ggsave("heatmap_top10_markers.pdf", width = 12, height = 16)
}


# ===== 9. 组合图 - UMAP + 关键marker =====
# UMAP + 几个关键marker的组合
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  theme(legend.position = "none") +
  labs(title = "Ciliated Cell Clusters")

p2 <- FeaturePlot(seurat_obj, features = "CCNO", pt.size = 0.5) +
  labs(title = "CCNO (Ciliogenesis)")

p3 <- FeaturePlot(seurat_obj, features = "IL8", pt.size = 0.5) +
  labs(title = "IL8 (Inflammation)")

p4 <- FeaturePlot(seurat_obj, features = "HLA-DRA", pt.size = 0.5) +
  labs(title = "HLA-DRA (Antigen Presentation)")

library(patchwork)
(p1 | p2) / (p3 | p4)

ggsave("combined_umap_features.pdf", width = 14, height = 12)
