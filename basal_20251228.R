# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
output_dir <- '/home/h2048/data/R/1228/basal'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
library(reticulate)
library(harmony)
library(ggplot2)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
h5ad_file1 <- "/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2/Basal/batch_aware/cnmf_analysis_k40/Basal_with_cnmf_k40.h5ad"
cat("Reading first h5ad file...\n")
seurat_obj <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)
# seurat_obj <- readRDS('basal_filtered_20251228.rds')
GetH5ad(seurat_obj, 'basal_filtered_20260104.h5ad')
seurat_obj <- NormalizeData(seurat_obj) #归一化
library(Seurat)
library(dplyr)

# ==== 0) 基础设置 ====
# seurat_obj <- readRDS("...")  # 你自己已加载的话可注释
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "RNA_snn_res.1.5" # 如果你的列名不同，在这里改

# # percent.mt（若没有就补上）
# if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
#   seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
#     seurat_obj,
#     pattern = "^MT-"
#   )
# }

# # ==== 1) 需要验证的基因集合 ====
# genes <- list(
#   basal_core = c("KRT5", "TP63", "KRT15", "KRT17", "ITGA6", "NGFR", "COL17A1"),
#   epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
#   endothelial = c("PECAM1", "VWF", "KDR", "EMCN"),
#   pericyte_smc = c("RGS5", "CSPG4", "ACTA2", "TAGLN", "MCAM"),
#   fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL6A1"),
#   immune = c("PTPRC", "LST1", "LYZ", "MS4A1", "FCER1A"),
#   rbc = c("HBA1", "HBA2", "HBG1")
# )

# # 只保留对象里存在的基因
# genes_present <- lapply(genes, \(x) intersect(x, rownames(seurat_obj)))
# genes_flat <- unique(unlist(genes_present))

# # ==== 2) DotPlot：一眼判断“是不是basal/是不是漏进来” ====
# pdf("validate_lineage_markers_dotplot.pdf", width = 14, height = 6)
# print(DotPlot(seurat_obj, features = genes_flat) + RotatedAxis())
# dev.off()

# # ==== 3) QC：双细胞/低质量（看 cluster 层面的偏移） ====
# meta <- seurat_obj@meta.data %>%
#   mutate(cluster = as.character(Idents(seurat_obj)))

# qc_summary <- meta %>%
#   group_by(cluster) %>%
#   summarise(
#     n_cells = n(),
#     nCount_median = median(nCount_RNA),
#     nFeature_median = median(nFeature_RNA),
#     percent_mt_median = median(percent.mt),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(nFeature_median))

# write.csv(
#   qc_summary,
#   "validate_qc_summary_by_cluster.csv",
#   row.names = FALSE,
#   quote = FALSE
# )

# pdf("validate_qc_violin.pdf", width = 12, height = 4)
# print(VlnPlot(
#   seurat_obj,
#   features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
#   pt.size = 0,
#   ncol = 3
# ))
# dev.off()

# # ==== 4) cluster 平均表达（用于表格快速判定） ====
# avg <- AverageExpression(
#   seurat_obj,
#   features = genes_flat,
#   assays = "RNA",
#   slot = "data"
# )$RNA
# avg_df <- as.data.frame(t(avg)) # rows=clusters, cols=genes
# write.csv(avg_df, "validate_marker_avgexpr_by_cluster.csv", quote = FALSE)

# # ==== 5) 可选：module score（更稳的“成套程序”证据） ====
# for (nm in names(genes_present)) {
#   if (length(genes_present[[nm]]) >= 3) {
#     seurat_obj <- AddModuleScore(
#       seurat_obj,
#       features = list(genes_present[[nm]]),
#       name = paste0("MS_", nm, "_"),
#       assay = "RNA",
#       search = FALSE
#     )
#   }
# }

# ms_cols <- grep("^MS_", colnames(seurat_obj@meta.data), value = TRUE)
# if (length(ms_cols) > 0) {
#   pdf("validate_modulescore_violin.pdf", width = 16, height = 5)
#   print(VlnPlot(seurat_obj, features = ms_cols, pt.size = 0, ncol = 4))
#   dev.off()
# }
# gene2 <- c('MYH11','CNN1','DES','LTF','BPIFA1','SCGB3A1')

# pdf("validate_lineage_markers_dotplot_2.pdf", width = 14, height = 6)
# p0 <- DotPlot(seurat_obj, features = gene2) + RotatedAxis()
# print(p0)
# dev.off()

# write.csv(p0$data,'gene2.csv')

# pdf("qc_umap_harmony.pdf", width = 12, height = 9)
# print(
#   DimPlot(seurat_obj, group.by = 'RNA_snn_res.1.5', label = TRUE, raster = TRUE) +
#     ggtitle(paste0("UMAP (Harmony) - ", 'RNA_snn_res.1.5'))
# )
# print(
#   DimPlot(seurat_obj, group.by = 'dataset', raster = TRUE) +
#     ggtitle(paste0("UMAP (Harmony) - ", 'dataset'))
# )
# print(
#   DimPlot(seurat_obj, group.by = 'tissue', raster = TRUE) +
#     ggtitle(paste0("UMAP (Harmony) - ", 'tissue'))
# )
# dev.off()

# saveRDS(seurat_obj, "basal_validated_with_scores.rds")

# # ---- A) 共表达比例：Basal(KRT5/TP63) 与 Serous(BPIFA1/LTF) 是否在同一细胞里同时高 ----
# df <- FetchData(seurat_obj, vars = c("seurat_clusters","nFeature_RNA","nCount_RNA",
#                                     "KRT5","TP63","BPIFA1","LTF","EPCAM"))

# # 用“>0”做最保守的共表达；你也可以改成更严格阈值（比如 >1）
# df <- df %>%
#   mutate(
#     basal_pos = (KRT5 > 0) & (TP63 > 0),
#     serous_pos = (BPIFA1 > 0) | (LTF > 0),
#     basal_serous_co = basal_pos & serous_pos
#   )

# co_tab <- df %>%
#   group_by(seurat_clusters) %>%
#   summarise(
#     n = n(),
#     co_rate = mean(basal_serous_co),
#     basal_rate = mean(basal_pos),
#     serous_rate = mean(serous_pos),
#     nFeature_med = median(nFeature_RNA),
#     nCount_med = median(nCount_RNA),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(co_rate))

# write.csv(co_tab, "check_basal_serous_coexpression_by_cluster.csv", row.names = FALSE)

# # 聚焦看 cluster 24
# subset(co_tab, seurat_clusters %in% c("24","18","22","9","19"))

# # ---- B) 可视化：cluster 24 的 KRT5 vs BPIFA1 是否呈“全体双阳性云团” ----
# pdf("check_cluster24_KRT5_vs_BPIFA1.pdf", width = 6, height = 5)
# FeatureScatter(subset(seurat_obj, idents = "24"), feature1 = "KRT5", feature2 = "BPIFA1")
# dev.off()

# pdf("check_allclusters_KRT5_vs_BPIFA1.pdf", width = 6, height = 5)
# FeatureScatter(seurat_obj, feature1 = "KRT5", feature2 = "BPIFA1")
# dev.off()

# # ---- C) 复杂度对比：cluster 24 是否明显高于其他 basal 主群（doublet 特征） ----
# pdf("check_nFeature_by_cluster.pdf", width = 12, height = 4)
# VlnPlot(seurat_obj, features = c("nFeature_RNA","nCount_RNA"), pt.size = 0, ncol = 2)
# dev.off()

# df <- FetchData(seurat_obj, vars = c("seurat_clusters","nFeature_RNA","nCount_RNA",
#                                     "KRT5","TP63","BPIFA1","LTF","EPCAM"))

# # 用 log-normalized(data slot) 的经验阈值：>1 约等于“不是零星泄漏”
# df2 <- df %>%
#   mutate(
#     basal_pos  = (KRT5 > 1) & (TP63 > 1),
#     serous_pos = (BPIFA1 > 1) | (LTF > 1),
#     co_pos     = basal_pos & serous_pos
#   )

# tab2 <- df2 %>%
#   group_by(seurat_clusters) %>%
#   summarise(
#     n = n(),
#     co_rate_strict = mean(co_pos),
#     basal_rate_strict = mean(basal_pos),
#     serous_rate_strict = mean(serous_pos),
#     nFeature_med = median(nFeature_RNA),
#     nCount_med   = median(nCount_RNA),
#     .groups = "drop"
#   ) %>% arrange(desc(co_rate_strict))

# write.csv(tab2, "coexpression_strict_by_cluster.csv", row.names = FALSE)
# tab2 %>% filter(seurat_clusters %in% c("24","18","9","19","22"))

# pdf("cluster24_KRT5_vs_BPIFA1.pdf", 6, 5)
# FeatureScatter(subset(seurat_obj, idents = "24"), feature1 = "KRT5", feature2 = "BPIFA1")
# dev.off()

# pdf("cluster24_qc_density.pdf", 10, 4)
# print(VlnPlot(subset(seurat_obj, idents="24"), features=c("nFeature_RNA","nCount_RNA"), pt.size=0, ncol=2))
# dev.off()

seurat_obj <- subset(
  seurat_obj,
  subset = RNA_snn_res.1.5 %in% c('20', '23'),
  invert = TRUE
)
seurat_obj$seurat_clusters <- NULL

# 1) 要删除的 metadata 列：leiden_* + (可选) RNA_snn_res.* / SCT_snn_res.*
pat_drop <- c(
  "^leiden_Epithelial", # leiden_Epithelial_res0.2/0.4/... + leiden_Epithelial
  "^leiden_bbknn", # leiden_bbknn_res1.2/...
  "^leiden_harmony", # leiden_harmony_res1.2/... + leiden_harmony
  "^RNA_snn_res\\.", # Seurat FindClusters 生成的 RNA_snn_res.X（如果存在）
  "^SCT_snn_res\\." # 如你用过 SCT（如果存在）
)

md <- seurat_obj@meta.data
cols_drop <- unique(unlist(lapply(pat_drop, \(p) {
  grep(p, colnames(md), value = TRUE)
})))

cat("Will drop metadata columns (n=", length(cols_drop), "):\n", sep = "")
print(cols_drop)

# 2) 防止当前 Idents 依赖被删列：如果你之前 Idents(seurat_obj) <- "leiden_*"，建议先切回一个保留列
# （按你的对象实际情况改，比如 "seurat_clusters" / "Manual_Annotation" / "celltype"）
if ("seurat_clusters" %in% colnames(md)) {
  Idents(seurat_obj) <- "seurat_clusters"
}

# 3) 删除 metadata 列
if (length(cols_drop) > 0) {
  seurat_obj@meta.data <- md[, setdiff(colnames(md), cols_drop), drop = FALSE]
}

cat("Current reductions:\n")
print(Reductions(seurat_obj))

# 1) Safety check
if (!"cnmf_usages" %in% Reductions(seurat_obj)) {
  stop("Reduction 'cnmf_usages' not found in this Seurat object.")
}

# 2) Drop all other reductions
keep_red <- "cnmf_usages"
drop_red <- setdiff(Reductions(seurat_obj), keep_red)

seurat_obj@reductions <- seurat_obj@reductions[keep_red]

cat("Dropped reductions:\n")
print(drop_red)
cat("Remaining reductions:\n")
print(Reductions(seurat_obj))


seurat_obj <- NormalizeData(seurat_obj, scale.factor = 1e4) #归一化
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
) #寻找变异基因
seurat_obj <- ScaleData(seurat_obj) #标准化

# 使用高变基因进行主成分分析，降低数据维度
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
seurat_obj <- RunHarmony(
  object = seurat_obj,
  group.by.vars = c("sample"),
  theta = c(8), # Higher theta for more diverse clustering
  lambda = c(2), # Higher lambda to reduce overcorrection
  sigma = 0.1, # Lower sigma for tighter clusters
  nclust = 30, # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20,
  early_stop = TRUE,
  dims = 1:30 # More iterations for better convergence
)
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 150,
  n.trees = 500,
  min.dist = 0.6,
  # ,
  # learning.rate = 0.2, # 相对保守的学习率
  # n.epochs = 1400, # 增加迭代次数补偿较小的学习率
  # spread = 1.2,
  # repulsion.strength = 1.1,

  metric = "correlation"
)
# Find neighbors
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  k.param = 150
)

seurat_obj <- FindClusters(
  seurat_obj,
  algorithm = 4,
  group.singletons = TRUE,
  resolution = 1, # 多个分辨率
  verbose = TRUE
)

pdf("qc_umap_harmony.pdf", width = 12, height = 9)
print(
  DimPlot(seurat_obj, group.by = 'RNA_snn_res.1', label = TRUE, raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'RNA_snn_res.1'))
)
print(
  DimPlot(seurat_obj, group.by = 'dataset', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'dataset'))
)
print(
  DimPlot(seurat_obj, group.by = 'tissue', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'tissue'))
)
dev.off()

# 0) 选择聚类/分组作为身份（按需改成你的meta列名，如 "cell_type"）
Idents(seurat_obj) <- "RNA_snn_res.1"

saveRDS(seurat_obj, "basal_filtered_20251228.rds")
