# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
output_dir <- '/home/h2048/data/R/0105/basal'
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
h5ad_file1 <- "/home/h2048/data/R/1228/basal/ciliated_bbknn_integrated.h5ad"
cat("Reading first h5ad file...\n")
seurat_obj <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)
# seurat_obj <- readRDS('basal_filtered_20251228.rds')
# GetH5ad(seurat_obj,'basal_filtered_20260104.h5ad')
seurat_obj <- NormalizeData(seurat_obj) #归一化


# ==== 0) 基础设置 ====
# seurat_obj <- readRDS("...")  # 你自己已加载的话可注释
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "leiden_res1.0" # 如果你的列名不同，在这里改

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

# seurat_obj <- subset(
#   seurat_obj,
#   subset = RNA_snn_res.1.5 %in% c('20', '23'),
#   invert = TRUE
# )
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
  DimPlot(seurat_obj, group.by = 'leiden_res1.0', label = TRUE, raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'leiden_res1.0'))
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

# Find all markers using Wilcoxon test
library(Seurat)
library(future)

# Enable parallel processing
plan("multisession", workers = 8) # Adjust based on your CPU cores
options(future.globals.maxSize = 8000 * 1024^2) # 8GB max object size

# Find all markers
all_markers <- FindAllMarkers(
  object = seurat_obj,
  test.use = "wilcox", # Wilcoxon rank sum test
  only.pos = TRUE, # Only positive markers
  min.pct = 0.25, # Min % cells expressing in either group
  logfc.threshold = 0.25, # Min log fold change
  verbose = TRUE
)

# Filter significant markers
sig_markers <- subset(all_markers, p_val_adj < 0.05)

cat(sprintf(
  "Found %d significant markers across %d clusters\n",
  nrow(sig_markers),
  length(unique(sig_markers$cluster))
))

# Get top markers per cluster
top_markers <- sig_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Save results
write.csv(all_markers, "all_markers_wilcox.csv", row.names = FALSE)
write.csv(top_markers, "top100_markers_per_cluster.csv", row.names = FALSE)

cat("Marker detection completed\n")

marker_genes <- unique(top_markers$gene)
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]

# 3) 仅对这些marker做Scale + 热图
seurat_obj <- ScaleData(seurat_obj, features = marker_genes, verbose = FALSE)
# c1 <- c('3','11','16')
# c2 <- c('17')
# c3 <- c('10')
# c4 <- c('4')
# c5 <- c('1','7','8')
# c6 <- c('5','6','12')
# c7 <- c('2')
# c8 <- c('15')
# c9 <- c('13')
# c10 <- c('14','19')

p <- DoHeatmap(
  seurat_obj,
  features = marker_genes,
  group.by = "leiden_res1.0",
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

cluster_col <- "leiden_res1.0"

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
saveRDS(seurat_obj, "basal_filtered_20251228.rds")


features <- unique(c(
  # Basal core
  "TP63",
  "KRT5",
  "KRT14",
  "KRT15",
  "KRT17",
  "COL17A1",
  "ITGA6",
  "LAMB3",
  "LAMC2",
  "BCAM",
  # Squamous/metaplasia
  "KRT4",
  "KRT13",
  "SPRR3",
  "CLCA2",
  "LY6D",
  "LGALS7",
  # Cycling
  "MKI67",
  "TOP2A",
  "TK1",
  "TYMS",
  # Inflammatory/IEG
  "IL8",
  "CXCL2",
  "CXCL1",
  "AREG",
  "ATF3",
  "TNFAIP3",
  # ECM/Repair
  "FN1",
  "POSTN",
  "WNT4",
  "IGFBP7",
  "CTGF",
  # Secretory/Club-like
  "SCGB1A1",
  "KRT7",
  "WFDC2",
  "CXCL17",
  "MUC1",
  # Ciliated check
  "FOXJ1",
  "TPPP3",
  'KRT8',
  'KRT18',
  'KRT19',
  'EPCAM',
  'TACSTD2',
  'SCGB3A1',
  'SCGB3A2',
  'PIGR',
  'BPIFA1',
  'SLPI'
))

p <- DotPlot(seurat_obj, features = features, group.by = "leiden_res1.0") +
  RotatedAxis()

write.csv(p$data, 'dotplot.csv')

# =========================
# Split basal / secretory / squamous and export to h5ad
# =========================

# 1) User config
# 1) User config
CLUSTER_COL <- "leiden_res1.0"
OUTDIR <- "h5ad_split_20260104"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 2) Cluster mapping
basal_ids <- as.character(c(1, 2, 3, 5, 7, 8, 9, 10))
secretory_ids <- as.character(c(0, 4, 6))
squamous_ids <- as.character(c(11))

stopifnot(CLUSTER_COL %in% colnames(seurat_obj@meta.data))

# 3) Write lineage label into metadata (cluster -> lineage)
clu <- as.character(seurat_obj[[CLUSTER_COL, drop = TRUE]])

seurat_obj$epi_lineage_3way <- "other"
seurat_obj$epi_lineage_3way[clu %in% basal_ids] <- "basal"
seurat_obj$epi_lineage_3way[clu %in% secretory_ids] <- "secretory"
seurat_obj$epi_lineage_3way[clu %in% squamous_ids] <- "squamous"

# 4) Subset by metadata
obj_basal <- subset(seurat_obj, subset = epi_lineage_3way == "basal")
obj_secretory <- subset(seurat_obj, subset = epi_lineage_3way == "secretory")
obj_squamous <- subset(seurat_obj, subset = epi_lineage_3way == "squamous")

# 5) Export
GetH5ad(obj_basal, file.path(OUTDIR, "basal_final_20260104.h5ad"))
GetH5ad(obj_secretory, file.path(OUTDIR, "secretory_final_20260104.h5ad"))
GetH5ad(obj_squamous, file.path(OUTDIR, "squamous_final_20260104.h5ad"))

assay_use <- DefaultAssay(seurat_obj)

cat("Original layers:\n")
print(Layers(seurat_obj[[assay_use]]))

cat("Basal layers after subset:\n")
print(Layers(obj_basal[[assay_use]]))

cat("Secretory layers after subset:\n")
print(Layers(obj_secretory[[assay_use]]))

cat("Squamous layers after subset:\n")
print(Layers(obj_squamous[[assay_use]]))

seurat_obj_main <- readRDS(
  '/home/h2048/data/R/1228/basal/basal_filtered_20251228.rds'
)


inspect_layer <- function(
  obj,
  layer = c("counts", "data"),
  assay = NULL,
  n_show = 5
) {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }

  get_layer <- function(obj, assay, layer) {
    # Seurat v5 推荐 LayerData；旧接口用 GetAssayData(slot=)
    out <- tryCatch(
      {
        SeuratObject::LayerData(obj[[assay]], layer = layer)
      },
      error = function(e) {
        GetAssayData(obj, assay = assay, slot = layer)
      }
    )
    out
  }

  x <- get_layer(obj, assay, layer)

  cat("\n==============================\n")
  cat("Assay:", assay, " | Layer:", layer, "\n")
  cat("Class:", paste(class(x), collapse = ", "), "\n")
  cat("Dim (features x cells):", paste(dim(x), collapse = " x "), "\n")

  # 稀疏矩阵信息
  nnz <- if ("dgCMatrix" %in% class(x)) length(x@x) else sum(x != 0)
  total <- prod(dim(x))
  cat("Non-zeros:", nnz, sprintf("(%.3f%%)", 100 * nnz / total), "\n")

  # 值域与分位数（只看非零值更有意义）
  vals <- if ("dgCMatrix" %in% class(x)) x@x else as.vector(x)
  if (length(vals) > 0) {
    cat("Value range:", sprintf("[%.4g, %.4g]", min(vals), max(vals)), "\n")
    qs <- quantile(
      vals,
      probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1),
      na.rm = TRUE
    )
    cat("Quantiles (non-zero):\n")
    print(qs)
    # 是否“基本整数”（counts 应该是）
    is_integerish <- mean(abs(vals - round(vals)) < 1e-8) > 0.999
    cat("Mostly integer values?:", is_integerish, "\n")
  } else {
    cat("No non-zero values found.\n")
  }

  # 每细胞总量 / 检测到的基因数（对 counts/data 都能看）
  libsize <- Matrix::colSums(x)
  nfeat <- Matrix::colSums(x > 0)
  cat("Per-cell library size summary:\n")
  print(summary(as.numeric(libsize)))
  cat("Per-cell detected features summary:\n")
  print(summary(as.numeric(nfeat)))

  # 随机挑几个基因/细胞看一个小切片（避免巨大矩阵打印）
  set.seed(1)
  genes <- sample(rownames(x), min(n_show, nrow(x)))
  cells <- sample(colnames(x), min(n_show, ncol(x)))
  cat(
    "\nSmall slice (",
    length(genes),
    "genes x",
    length(cells),
    "cells ):\n",
    sep = ""
  )
  print(as.matrix(x[genes, cells, drop = FALSE]))
}

# 对三个对象分别看 counts/data
objs <- list(
  all = seurat_obj_main
  # ,
  # basal = obj_basal,
  # secretory = obj_secretory,
  # squamous = obj_squamous
)

for (nm in names(objs)) {
  cat("\n\n##########", nm, "##########\n")
  inspect_layer(objs[[nm]], layer = "counts")
  inspect_layer(objs[[nm]], layer = "data")
}


# =========================
# Split from seurat_obj_main by cell IDs and export h5ad
# =========================

OUTDIR <- "h5ad_split_20260104_from_main"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- 0) Assay sanity: 强烈建议确保 DefaultAssay 是 RNA（原始counts通常在RNA）
# 如果你的原始counts在别的assay（例如 "RNA" / "SCT"），这里改成对应的
DefaultAssay(seurat_obj_main) <- "RNA"

# ---- 1) Collect cell IDs from the previously-split objects (even if they lost counts)
cells_basal <- colnames(obj_basal)
cells_secretory <- colnames(obj_secretory)
cells_squamous <- colnames(obj_squamous)

# ---- 2) Safety checks: disjoint & overlap
dup_bs <- intersect(cells_basal, cells_secretory)
dup_bq <- intersect(cells_basal, cells_squamous)
dup_sq <- intersect(cells_secretory, cells_squamous)
dup_all <- unique(c(dup_bs, dup_bq, dup_sq))
stopifnot(length(dup_all) == 0)

# ---- 3) Match cells to seurat_obj_main (in case some cells are missing)
main_cells <- colnames(seurat_obj_main)

cells_basal_in <- intersect(cells_basal, main_cells)
cells_secretory_in <- intersect(cells_secretory, main_cells)
cells_squamous_in <- intersect(cells_squamous, main_cells)

miss_basal <- setdiff(cells_basal, main_cells)
miss_secretory <- setdiff(cells_secretory, main_cells)
miss_squamous <- setdiff(cells_squamous, main_cells)

if (length(miss_basal) > 0) {
  warning("Basal cells not found in seurat_obj_main: ", length(miss_basal))
}
if (length(miss_secretory) > 0) {
  warning(
    "Secretory cells not found in seurat_obj_main: ",
    length(miss_secretory)
  )
}
if (length(miss_squamous) > 0) {
  warning(
    "Squamous cells not found in seurat_obj_main: ",
    length(miss_squamous)
  )
}

# 再次确保导出三份没有重叠（按 main 里实际存在的cells）
dup2 <- unique(c(
  intersect(cells_basal_in, cells_secretory_in),
  intersect(cells_basal_in, cells_squamous_in),
  intersect(cells_secretory_in, cells_squamous_in)
))
stopifnot(length(dup2) == 0)

# ---- 4) Write metadata label into seurat_obj_main for reproducibility
seurat_obj_main$epi_lineage_3way <- "other"
seurat_obj_main$epi_lineage_3way[cells_basal_in] <- "basal"
seurat_obj_main$epi_lineage_3way[cells_secretory_in] <- "secretory"
seurat_obj_main$epi_lineage_3way[cells_squamous_in] <- "squamous"

# ---- 5) Subset from seurat_obj_main by cells (this preserves original layers from main)
obj_basal_main <- subset(seurat_obj_main, cells = cells_basal_in)
obj_secretory_main <- subset(seurat_obj_main, cells = cells_secretory_in)
obj_squamous_main <- subset(seurat_obj_main, cells = cells_squamous_in)

# ---- 6) Quick layer checks (Seurat v5: Layers(); fallback if not available)
check_layers <- function(obj, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  cat("\n---", deparse(substitute(obj)), " Assay:", assay, "---\n")
  out <- tryCatch(
    {
      Layers(obj[[assay]])
    },
    error = function(e) {
      # old-style: just report slots existence
      c("counts", "data", "scale.data")
    }
  )
  print(out)
}

check_layers(obj_basal_main)
check_layers(obj_secretory_main)
check_layers(obj_squamous_main)

# ---- 7) Export to h5ad using your helper
# 注意：不要 DietSeurat（除非你明确要丢 reduction/graph/scale.data 等）
GetH5ad(obj_basal_main, file.path(OUTDIR, "basal_final_20260104.h5ad"))
GetH5ad(obj_secretory_main, file.path(OUTDIR, "secretory_final_20260104.h5ad"))
GetH5ad(obj_squamous_main, file.path(OUTDIR, "squamous_final_20260104.h5ad"))

# ---- 8) (Optional) export full main with lineage column
# GetH5ad(seurat_obj_main, file.path(OUTDIR, "all_with_epi_lineage_20260104.h5ad"))
