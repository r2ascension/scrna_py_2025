# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
output_dir <- '/home/h2048/data/R/0104/B'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
library(reticulate)
library(harmony)
library(ggplot2)
library(patchwork)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
# H5AD_FILE <- "/home/h2048/data/R/1228/B/b_final_20260103.rds"
h5ad_file1 <- "/home/h2048/data/R/1228/B/ciliated_bbknn_integrated.h5ad"
# cat("Reading first h5ad file...\n")
seurat_obj <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)
# seurat_obj <- readRDS(H5AD_FILE)
# GetH5ad(seurat_obj,'b_filtered_20260104.h5ad')
seurat_obj <- NormalizeData(seurat_obj) #归一化
library(Seurat)
library(dplyr)

# ==== 0) 基础设置 ====
# seurat_obj <- readRDS("...")  # 你自己已加载的话可注释
DefaultAssay(seurat_obj) <- "RNA"

ident_col <- "leiden_bbknn_res1.0"
# Idents(seurat_obj) <- ident_col

# ===== 1) 基因面板（按需增减）=====
panel <- list(
  B_core = c("MS4A1", "CD79A", "CD74", "HLA-DRA", "CD37", "CD22"),
  GC = c("BCL6", "AICDA", "RGS13", "MME", "LMO2", "SERPINA9"),
  ABC = c("FCRL4", "ITGAX", "TBX21", "FCRL5", "CXCR3"),
  Plasma = c("PRDM1", "XBP1", "MZB1", "JCHAIN", "SDC1"),
  Activate = c(
    "MIR155HG",
    "CD69",
    "DUSP1",
    "DUSP2",
    "TNFAIP3",
    "NFKBIA",
    "BCL2A1"
  ),
  Contam = c(
    "TRAC",
    "TRBC1",
    "TRBC2",
    "EPCAM",
    "KRT19",
    "KRT8",
    "STATH",
    "BPIFA1",
    "LTF",
    "ALDH3A1",
    "MUC5AC"
  )
)

genes <- unique(unlist(panel))
genes <- intersect(genes, rownames(seurat_obj))
cat(sprintf(
  "Genes kept: %d / %d\n",
  length(genes),
  length(unique(unlist(panel)))
))

# ===== 2) cluster大小 =====
print(table(Idents(seurat_obj)))

# ===== 3) DotPlot（总览）=====
pdf("Bcell_marker_dotplot.pdf", width = 18, height = 7)
p1 <- DotPlot(seurat_obj, features = genes) +
  RotatedAxis() +
  ggtitle("B cell validation markers (DotPlot)") +
  theme(plot.title = element_text(face = "bold"))
print(p1)
dev.off()

# ===== 4) FeaturePlot（重点看污染/状态）=====
fp_genes <- intersect(
  c(
    "MS4A1",
    "CD79A",
    "CD74",
    "BCL6",
    "AICDA",
    "MME",
    "SERPINA9",
    "FCRL4",
    "ITGAX",
    "PRDM1",
    "XBP1",
    "MZB1",
    "JCHAIN",
    "TRBC2",
    "TRAC",
    "EPCAM",
    "KRT19",
    "STATH",
    "ALDH3A1",
    "MUC5AC"
  ),
  rownames(seurat_obj)
)

pdf("Bcell_marker_featureplot.pdf", width = 16, height = 12)
p2 <- FeaturePlot(
  seurat_obj,
  features = fp_genes,
  reduction = "umap",
  ncol = 5,
  raster = TRUE
)
print(p2)
dev.off()

# # ===== 5) 平均表达表（每cluster一列，便于你贴结果）=====
# avg <- AverageExpression(seurat_obj, features = genes, assays = "RNA", slot = "data")$RNA
# write.csv(avg, "Bcell_marker_AverageExpression_byCluster.csv", quote = FALSE)
# write.csv(p1$data,'dotplot.csv')

# sig <- list(
#   epi = c("EPCAM","KRT19","KRT8","ALDH3A1","MUC5AC"),
#   T   = c("TRAC","TRBC1","TRBC2"),
#   plasma = c("PRDM1","XBP1","MZB1","JCHAIN","SDC1"),
#   GC  = c("BCL6","AICDA","RGS13","MME","SERPINA9","LMO2"),
#   ABC = c("FCRL4","FCRL5","ITGAX","TBX21","CXCR3"),
#   act = c("CD69","MIR155HG","DUSP1","DUSP2","NFKBIA","TNFAIP3","BCL2A1")
# )

# for(nm in names(sig)) seurat_obj <- AddModuleScore(seurat_obj, list(intersect(sig[[nm]], rownames(seurat_obj))), name = nm)
# VlnPlot(seurat_obj, features = paste0(names(sig), "1"), group.by = "seurat_clusters", pt.size = 0)

# pdf("qc_umap_harmony.pdf", width = 12, height = 9)
# print(
#   DimPlot(seurat_obj, group.by = 'leiden_bbknn_res1.0', label = TRUE, raster = TRUE) +
#     ggtitle(paste0("UMAP (Harmony) - ", 'leiden_res1.0'))
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

# # ===== Create summary table for all signatures =====
# library(dplyr)

# # Initialize empty list to store results
# summary_tables <- list()

# for (nm in names(sig)) {
#   score_col <- paste0(nm, "1")

#   # Calculate statistics by cluster
#   summary_tables[[nm]] <- seurat_obj@meta.data %>%
#     group_by(leiden_bbknn_res1.0) %>%
#     summarise(
#       n_cells = n(),
#       mean = mean(.data[[score_col]], na.rm = TRUE),
#       median = median(.data[[score_col]], na.rm = TRUE),
#       sd = sd(.data[[score_col]], na.rm = TRUE),
#       q25 = quantile(.data[[score_col]], 0.25, na.rm = TRUE),
#       q75 = quantile(.data[[score_col]], 0.75, na.rm = TRUE)
#     ) %>%
#     mutate(signature = nm) %>%
#     arrange(desc(mean))

#   cat(sprintf("\n===== %s Summary by Cluster =====\n", nm))
#   print(summary_tables[[nm]])
# }

# # Combine all signatures into one table
# combined_summary <- bind_rows(summary_tables)
# print(combined_summary)

# # Optional: Save to file
# write.csv(combined_summary, "signature_summary_by_cluster.csv", row.names = FALSE)

# p_epi <- DotPlot(seurat_obj, features = c("EPCAM","KRT19","KRT8","KRT18","KRT5","KRT14","TP63"),
#         group.by = "leiden_bbknn_res1.0") + RotatedAxis()
# print(p_epi)

# seurat_obj <- subset(
#   seurat_obj,
#   subset = leiden_res1.0 %in% c('10'),
#   invert = TRUE
# )
# saveRDS(seurat_obj, "b_filtered_20260103.rds")

# seurat_obj$seurat_clusters <- NULL

# # 1) 要删除的 metadata 列：leiden_* + (可选) RNA_snn_res.* / SCT_snn_res.*
# pat_drop <- c(
#   "^leiden_Epithelial",     # leiden_Epithelial_res0.2/0.4/... + leiden_Epithelial
#   "^leiden_bbknn",          # leiden_bbknn_res1.2/...
#   "^leiden_harmony",        # leiden_harmony_res1.2/... + leiden_harmony
#   "^RNA_snn_res\\.",        # Seurat FindClusters 生成的 RNA_snn_res.X（如果存在）
#   "^SCT_snn_res\\."         # 如你用过 SCT（如果存在）
# )

# md <- seurat_obj@meta.data
# cols_drop <- unique(unlist(lapply(pat_drop, \(p) grep(p, colnames(md), value = TRUE))))

# cat("Will drop metadata columns (n=", length(cols_drop), "):\n", sep = "")
# print(cols_drop)

# # 2) 防止当前 Idents 依赖被删列：如果你之前 Idents(seurat_obj) <- "leiden_*"，建议先切回一个保留列
# # （按你的对象实际情况改，比如 "seurat_clusters" / "Manual_Annotation" / "celltype"）
# if ("seurat_clusters" %in% colnames(md)) {
#   Idents(seurat_obj) <- "seurat_clusters"
# }

# # 3) 删除 metadata 列
# if (length(cols_drop) > 0) {
#   seurat_obj@meta.data <- md[, setdiff(colnames(md), cols_drop), drop = FALSE]
# }

# cat("Current reductions:\n")
# print(Reductions(seurat_obj))

# # 1) Safety check
# if (!"cnmf_usages" %in% Reductions(seurat_obj)) {
#   stop("Reduction 'cnmf_usages' not found in this Seurat object.")
# }

# # 2) Drop all other reductions
# keep_red <- "cnmf_usages"
# drop_red <- setdiff(Reductions(seurat_obj), keep_red)

# seurat_obj@reductions <- seurat_obj@reductions[keep_red]

# cat("Dropped reductions:\n")
# print(drop_red)
# cat("Remaining reductions:\n")
# print(Reductions(seurat_obj))

# ===== Data Normalization and HVG Selection (Excluding IG genes) =====

# Normalization
seurat_obj <- NormalizeData(seurat_obj, scale.factor = 1e4)

# # --- Identify IG-related genes to exclude ---
# all_genes <- rownames(seurat_obj)

# # Additional genes to exclude
# exclude_patterns <- c(
#   "^IGH",
#   "^IGK",
#   "^IGL", # IG genes
#   "^TRA",
#   "^TRB",
#   "^TRG",
#   "^TRD", # TCR genes (optional)
#   "^MT-", # Mitochondrial genes (optional)
#   "^RP[SL]" # Ribosomal genes (optional)
# )

# # Find all IG genes
# ig_genes <- grep(
#   paste(exclude_patterns, collapse = "|"),
#   all_genes,
#   value = TRUE
# )
# cat(sprintf("Found %d IG-related genes to exclude\n", length(ig_genes)))

# cat("Examples:", head(ig_genes, 10), "\n")

# # --- Method 1: Exclude IG genes before HVG selection ---
# # Create a gene list excluding IG genes
# genes_to_use <- setdiff(all_genes, ig_genes)
# cat(sprintf(
#   "Using %d genes for HVG selection (excluded %d IG genes)\n",
#   length(genes_to_use),
#   length(ig_genes)
# ))

# Find variable features only from non-IG genes
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = TRUE
)

# # Manually filter out any IG genes that might have been selected
# variable_features <- VariableFeatures(seurat_obj)
# variable_features_filtered <- setdiff(variable_features, ig_genes)

# # Update variable features
# VariableFeatures(seurat_obj) <- variable_features_filtered
# cat(sprintf(
#   "Final HVG count: %d (removed %d IG genes from HVG list)\n",
#   length(variable_features_filtered),
#   length(variable_features) - length(variable_features_filtered)
# ))

# # --- Optional: Save excluded gene list for documentation ---
# write.table(
#   data.frame(gene = ig_genes),
#   file = "excluded_ig_genes.txt",
#   row.names = FALSE,
#   quote = FALSE
# )

# --- Continue with standard workflow ---
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# --- Harmony batch correction ---
seurat_obj <- RunHarmony(
  object = seurat_obj,
  group.by.vars = c("sample"),
  theta = c(3),
  lambda = c(1),
  sigma = 0.1,
  nclust = 30,
  reduction.use = "pca",
  max_iter = 20,
  early_stop = TRUE,
  dims = 1:30
)

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 20,
  n.trees = 500,
  min.dist = 0.3,
  metric = "correlation"
)

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  k.param = 20
)

seurat_obj <- FindClusters(
  seurat_obj,
  algorithm = 4,
  group.singletons = TRUE,
  resolution = 3,
  verbose = TRUE
)

# --- Visualization ---
pdf("qc_umap_harmony.pdf", width = 12, height = 9)
print(
  DimPlot(seurat_obj, group.by = 'leiden_res1.0', label = TRUE, raster = TRUE) +
    ggtitle("UMAP (Harmony) - Clustering (IG genes excluded from HVG)")
)
print(
  DimPlot(seurat_obj, group.by = 'dataset', raster = TRUE) +
    ggtitle("UMAP (Harmony) - Dataset")
)
print(
  DimPlot(seurat_obj, group.by = 'tissue', raster = TRUE)
)
print(
  DimPlot(seurat_obj, group.by = 'scanvi_predictions', raster = TRUE) +
    ggtitle("UMAP (Harmony) - scanvi_predictions")
)

print(
  DimPlot(seurat_obj, group.by = 'dominant_gep', raster = TRUE) +
    ggtitle("UMAP (Harmony) - dominant_gep")
)
dev.off()

Idents(seurat_obj) <- "leiden_res1.0"
# saveRDS(seurat_obj, "b_filtered_20260103.rds")

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
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Save results
write.csv(all_markers, "all_markers_wilcox.csv", row.names = FALSE)
write.csv(top_markers, "top10_markers_per_cluster.csv", row.names = FALSE)

cat("Marker detection completed\n")


p0 <- VlnPlot(
  seurat_obj,
  features = c(
    'CD19',
    'MS4A1',
    'CD79A',
    'TRAC',
    'CD3D',
    'IL7R',
    'LYZ',
    'LST1',
    'NKG7'
  ),
  group.by = 'leiden_res1.0',
  pt.size = 0.05,
  ncol = 2
)
print(p0)


library(Seurat)

# ===== settings =====
ASSAY_USE <- "RNA"
SLOT_USE <- "data" # log-normalized; if you want counts use "counts"

b_genes <- c("CD19", "MS4A1", "CD79A")
t_genes <- c("LCK", "CD3D", "CD3E")

# thresholds on normalized expression (adjust if needed)
b_thr <- 1
t_thr <- 1

# strictness: require >=k genes positive in each set
b_k <- 2
t_k <- 2

# ===== fetch expression safely =====
all_genes <- unique(c(b_genes, t_genes))
present <- intersect(all_genes, rownames(seurat_obj))
md <- FetchData(seurat_obj, vars = present, slot = SLOT_USE, assay = ASSAY_USE)

# fill missing genes as 0 (so code won't break)
for (g in setdiff(all_genes, colnames(md))) {
  md[[g]] <- 0
}

# ===== define B-high & T-high (vectorized) =====
b_pos <- rowSums(md[, b_genes, drop = FALSE] > b_thr) >= b_k
t_pos <- rowSums(md[, t_genes, drop = FALSE] > t_thr) >= t_k

bt_doublet_like <- b_pos & t_pos

cat(
  "B&T doublet-like cells:",
  sum(bt_doublet_like),
  " / ",
  nrow(md),
  sprintf(" (%.3f%%)\n", 100 * mean(bt_doublet_like))
)

# ===== tag + remove =====
seurat_obj$flag_BT_doublet_like <- bt_doublet_like

cells_drop <- rownames(md)[bt_doublet_like]
# optional: save list
# writeLines(cells_drop, "BT_doublet_like_cells.txt")

seurat_obj <- subset(
  seurat_obj,
  cells = setdiff(colnames(seurat_obj), cells_drop)
)
seurat_obj <- seurat_obj
# sanity check
table(seurat_obj$flag_BT_doublet_like)

# 假设 Idents(seurat_obj) 已经是你的 cluster
tab <- table(Idents(seurat_obj), seurat_obj$flag_BT_doublet_like)
print(tab)

# 看每个cluster中被flag比例
prop <- prop.table(tab, margin = 1)
print(round(100 * prop, 2))


library(Seurat)

ASSAY_USE <- "RNA"
SLOT_USE <- "data"
CLUSTER_COL <- "seurat_clusters" # 改成你的cluster列

Idents(seurat_obj) <- seurat_obj[[CLUSTER_COL, drop = TRUE]]

keep_present <- function(obj, gs) intersect(gs, rownames(obj))

# B core
B_core <- keep_present(
  seurat_obj,
  c("CD19", "MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA")
)

# Myeloid core (mono/mac)
My_core <- keep_present(
  seurat_obj,
  c(
    "LYZ",
    "LST1",
    "TYROBP",
    "AIF1",
    "FCGR3A",
    "MS4A7",
    "CTSS",
    "CST3",
    "LGALS3",
    "FCN1",
    "C1QA",
    "C1QB",
    "C1QC",
    "APOE",
    "MERTK",
    "SPI1",
    "CSF1R"
  )
)

# DC (CD1C+ etc)
DC_core <- keep_present(
  seurat_obj,
  c("CD1C", "FCER1A", "CLEC10A", "ITGAX", "CCR7")
)

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(B = B_core),
  name = "Score_B"
)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(My = My_core),
  name = "Score_My"
)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(DC = DC_core),
  name = "Score_DC"
)

# quick check
print(VlnPlot(
  seurat_obj,
  features = c("Score_B1", "Score_My1", "Score_DC1"),
  group.by = CLUSTER_COL,
  pt.size = 0,
  ncol = 3
))


md <- FetchData(
  seurat_obj,
  vars = c("Score_B1", "Score_My1", "Score_DC1"),
  slot = SLOT_USE,
  assay = ASSAY_USE
)

# 保守阈值（你也可以改成更激进：0.85/0.5）
my_thr <- as.numeric(quantile(md$Score_My1, 0.90, na.rm = TRUE))
dc_thr <- as.numeric(quantile(md$Score_DC1, 0.90, na.rm = TRUE))
b_low <- as.numeric(quantile(md$Score_B1, 0.40, na.rm = TRUE))

flag_myeloid <- (md$Score_My1 > my_thr | md$Score_DC1 > dc_thr) &
  (md$Score_B1 < b_low)

cat(
  "Myeloid-like removed:",
  sum(flag_myeloid),
  "/",
  nrow(md),
  sprintf("(%.3f%%)\n", 100 * mean(flag_myeloid))
)

seurat_obj$flag_myeloid_like <- flag_myeloid
cells_drop <- rownames(md)[flag_myeloid]

hard_my <- keep_present(
  seurat_obj,
  c("LST1", "TYROBP", "FCGR3A", "MS4A7", "C1QC", "APOE", "MERTK", "FCN1")
)
expr <- FetchData(
  seurat_obj,
  vars = hard_my,
  slot = SLOT_USE,
  assay = ASSAY_USE
)

hard_pos <- rowSums(expr > 1) >= 2 # 阈值1偏保守；可改0.5
flag_myeloid2 <- flag_myeloid & hard_pos

cat(
  "Myeloid-like (hard) removed:",
  sum(flag_myeloid2),
  "/",
  nrow(expr),
  sprintf("(%.3f%%)\n", 100 * mean(flag_myeloid2))
)

cells_drop2 <- rownames(expr)[flag_myeloid2]
seurat_obj_B <- subset(
  seurat_obj,
  cells = setdiff(colnames(seurat_obj), cells_drop2)
)


library(Seurat)
library(Matrix)

ASSAY_USE <- "RNA"
SLOT_USE <- "data"
CLUSTER_COL <- "seurat_clusters"

Idents(seurat_obj) <- seurat_obj[[CLUSTER_COL, drop = TRUE]]

keep_present <- function(obj, gs) intersect(gs, rownames(obj))
pos <- function(x, thr = 1) ifelse(is.na(x), FALSE, x > thr)

# --- hard myeloid markers (more specific) ---
my_hard <- keep_present(
  seurat_obj,
  c("LST1", "TYROBP", "FCER1G", "MS4A7", "FCGR3A", "S100A8", "S100A9", "LGALS3")
)
expr_my <- FetchData(
  seurat_obj,
  vars = my_hard,
  slot = SLOT_USE,
  assay = ASSAY_USE
)
for (g in setdiff(my_hard, colnames(expr_my))) {
  expr_my[[g]] <- 0
}

# require >=2 hard markers positive
flag_my_hard <- rowSums(expr_my > 1) >= 2
cat(
  "Hard-myeloid cells:",
  sum(flag_my_hard),
  "/",
  nrow(expr_my),
  sprintf("(%.3f%%)\n", 100 * mean(flag_my_hard))
)

seurat_obj$flag_my_hard <- flag_my_hard

plasma_genes <- keep_present(
  seurat_obj,
  c("JCHAIN", "MZB1", "XBP1", "TNFRSF17", "SDC1", "DERL3", "FKBP11")
)
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(plasma = plasma_genes),
  name = "Score_Plasma"
)
md2 <- FetchData(
  seurat_obj,
  vars = c("Score_Plasma1"),
  slot = SLOT_USE,
  assay = ASSAY_USE
)
plasma_hi <- md2$Score_Plasma1 >
  as.numeric(quantile(md2$Score_Plasma1, 0.80, na.rm = TRUE)) # 保守保护上20%

# 最终删除：硬髓系阳性 且 不是plasma高
flag_my_final <- seurat_obj$flag_my_hard & (!plasma_hi)
cat(
  "Myeloid removed (hard, plasma-protected):",
  sum(flag_my_final),
  "/",
  nrow(md2),
  sprintf("(%.3f%%)\n", 100 * mean(flag_my_final))
)

seurat_obj$flag_myeloid_like2 <- flag_my_final
seurat_obj_B <- subset(
  seurat_obj,
  cells = setdiff(colnames(seurat_obj), rownames(md2)[flag_my_final])
)
seurat_obj <- seurat_obj_B


t_core <- c("CD3D", "CD3E", "TRAC", "CD247", "LCK")
b_core <- c("CD19", "MS4A1", "CD79A", "CD74", "HLA-DRA")
my_hard <- c(
  "LST1",
  "TYROBP",
  "FCER1G",
  "MS4A7",
  "FCGR3A",
  "S100A8",
  "S100A9",
  "LYZ"
)

t_core <- intersect(t_core, rownames(seurat_obj))
b_core <- intersect(b_core, rownames(seurat_obj))
my_hard <- intersect(my_hard, rownames(seurat_obj))

DotPlot(
  seurat_obj,
  features = list(B = b_core, T = t_core, Myeloid = my_hard),
  group.by = "seurat_clusters"
) +
  RotatedAxis()


mdT <- FetchData(seurat_obj, vars = t_core, slot = "data")
thr <- 1
t_pos_k2 <- rowSums(mdT > thr) >= 2

tabT <- table(seurat_obj$seurat_clusters, t_pos_k2)
propT <- prop.table(tabT, 1)[, "TRUE"]
print(round(100 * propT, 2))

cells_drop_T <- names(t_pos_k2)[t_pos_k2]
seurat_obj_B2 <- subset(
  seurat_obj,
  cells = setdiff(colnames(seurat_obj), cells_drop_T)
)

# 重跑（至少UMAP+neighbors+clusters）
seurat_obj_B2 <- RunPCA(seurat_obj_B2, npcs = 30)
seurat_obj_B2 <- RunUMAP(seurat_obj_B2, dims = 1:30)
seurat_obj_B2 <- FindNeighbors(seurat_obj_B2, dims = 1:30)
seurat_obj_B2 <- FindClusters(seurat_obj_B2, resolution = 0.5)


t_anchor <- intersect(c("CD3D", "CD3E", "TRAC", "CD247"), rownames(seurat_obj))
t_core <- intersect(
  c("CD3D", "CD3E", "TRAC", "CD247", "LCK"),
  rownames(seurat_obj)
)

thr <- 1
mdA <- FetchData(seurat_obj, vars = t_anchor, slot = "data")
mdT <- FetchData(seurat_obj, vars = t_core, slot = "data")

t_anchor_pos <- rowSums(mdA > thr) >= 1
t_k2_pos <- rowSums(mdT > thr) >= 2

flag_T_like <- t_anchor_pos & t_k2_pos

cat(
  "T-like removed:",
  sum(flag_T_like),
  "/",
  nrow(mdT),
  sprintf("(%.3f%%)\n", 100 * mean(flag_T_like))
)

seurat_obj$flag_T_like <- flag_T_like
seurat_obj <- subset(
  seurat_obj,
  cells = setdiff(colnames(seurat_obj), names(flag_T_like)[flag_T_like])
)


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
  theta = c(3), # Higher theta for more diverse clustering
  lambda = c(1), # Higher lambda to reduce overcorrection
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
  n.neighbors = 50,
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
  k.param = 50
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
Idents(seurat_obj) <- "leiden_res1.0"


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


b_core <- c("CD79A", "MS4A1", "CD74", "HLA-DRA", "CD37", "CD72")
t_core <- c("CD3D", "CD3E", "TRAC", "CD247", "LCK")
my_hard <- c(
  "LST1",
  "TYROBP",
  "FCER1G",
  "MS4A7",
  "FCGR3A",
  "S100A8",
  "S100A9",
  "CTSS",
  "LYZ"
)
epi <- c("EPCAM", "KRT19", "KRT8", "KRT18", "MUC1", "KRT5", "KRT14", "DMBT1")

plasma <- c(
  "SDC1",
  "MZB1",
  "XBP1",
  "PRDM1",
  "IRF4",
  "JCHAIN",
  "DERL3",
  "SEC11C",
  "FKBP11",
  "HSP90B1",
  "TNFRSF17",
  "TNFRSF13B",
  "CD27"
)

mk <- unique(c(b_core, plasma, t_core, my_hard, epi))
mk <- intersect(mk, rownames(seurat_obj))

DotPlot(seurat_obj, features = mk, group.by = "seurat_clusters") + RotatedAxis()

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
})

# ====== Config ======
REDUCTION_USE <- "umap"
SLOT_USE <- "data" # 通常用 data；如果你想看counts可改成 "counts"
PT_SIZE <- 0.15
NCOL <- 4

# ====== Marker sets ======
b_core <- c("CD79A", "MS4A1", "CD74", "HLA-DRA", "CD37", "CD72")
plasma <- c(
  "SDC1",
  "MZB1",
  "XBP1",
  "PRDM1",
  "IRF4",
  "JCHAIN",
  "DERL3",
  "SEC11C",
  "FKBP11",
  "HSP90B1",
  "TNFRSF17",
  "TNFRSF13B",
  "CD27"
)
t_core <- c("CD3D", "CD3E", "TRAC", "CD247", "LCK")
my_hard <- c(
  "LST1",
  "TYROBP",
  "FCER1G",
  "MS4A7",
  "FCGR3A",
  "S100A8",
  "S100A9",
  "CTSS",
  "LYZ"
)
epi <- c("EPCAM", "KRT19", "KRT8", "KRT18", "MUC1", "KRT5", "KRT14", "DMBT1")

gene_sets <- list(
  B_core = b_core,
  Plasma = plasma,
  T_core = t_core,
  Myeloid_hard = my_hard,
  Epi = epi
)

# ====== Helper ======
plot_feature_set <- function(
  obj,
  genes,
  title,
  out_pdf,
  reduction = REDUCTION_USE,
  slot = SLOT_USE,
  pt.size = PT_SIZE,
  ncol = NCOL
) {
  genes_ok <- intersect(genes, rownames(obj))
  if (length(genes_ok) == 0) {
    message("[SKIP] ", title, ": no genes found in object.")
    return(invisible(NULL))
  }

  nrow <- ceiling(length(genes_ok) / ncol)
  pdf_w <- max(8, 3.2 * ncol)
  pdf_h <- max(6, 2.8 * nrow)

  p <- FeaturePlot(
    obj,
    features = genes_ok,
    reduction = reduction,
    slot = slot,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    pt.size = pt.size,
    ncol = ncol
  ) +
    plot_annotation(title = title)

  pdf(out_pdf, width = pdf_w, height = pdf_h, onefile = TRUE)
  print(p)
  dev.off()

  message("[OK] ", title, " -> ", out_pdf, " (n=", length(genes_ok), ")")
}

# ====== Run ======
for (nm in names(gene_sets)) {
  plot_feature_set(
    seurat_obj,
    genes = gene_sets[[nm]],
    title = nm,
    out_pdf = paste0("FeaturePlot_", nm, ".pdf")
  )
}


# -----------------------------
# DotPlot marker panel (B / T / Plasma / Myeloid / Epi / Cycling)
# -----------------------------
ASSAY_USE <- DefaultAssay(seurat_obj) # 或者 "RNA"
SLOT_USE <- "data" # 常用：data；若你想看counts可改 "counts"

b_core <- c(
  "CD79A",
  "MS4A1",
  "CD74",
  "HLA-DRA",
  "CD37",
  "CD72",
  "CD19",
  "BANK1"
)
plasma <- c(
  "MZB1",
  "XBP1",
  "JCHAIN",
  "PRDM1",
  "SDC1",
  "TNFRSF17",
  "DERL3",
  "SEC11C",
  "FKBP11",
  "IGHG1",
  "IGHG4",
  "IGHA1",
  "IGHA2"
)
t_core <- c(
  "CD3D",
  "CD3E",
  "TRAC",
  "CD247",
  "LCK",
  "TRBC1",
  "TRBC2",
  "IL7R",
  "CCR7",
  "NKG7"
)
my_hard <- c(
  "LST1",
  "TYROBP",
  "FCER1G",
  "MS4A7",
  "FCGR3A",
  "S100A8",
  "S100A9",
  "CTSS",
  "LYZ",
  "LGALS3",
  "CST3",
  "SPI1",
  "C1QA",
  "C1QB",
  "C1QC",
  "MERTK",
  "ITGAX",
  "FCER1A",
  "CLEC10A"
)
epi <- c(
  "EPCAM",
  "KRT19",
  "KRT8",
  "KRT18",
  "MUC1",
  "KRT5",
  "KRT14",
  "BPIFA1",
  "STATH",
  "PIP",
  "LTF",
  "DMBT1"
)
cycling <- c(
  "MKI67",
  "TOP2A",
  "HMGB2",
  "TUBA1B",
  "TYMS",
  "MCM4",
  "MCM6",
  "PCNA",
  "STMN1",
  "UBE2C",
  "NUSAP1",
  "CENPF"
)

mk <- unique(c(b_core, plasma, t_core, my_hard, epi, cycling))
mk <- intersect(mk, rownames(seurat_obj))

p <- DotPlot(
  seurat_obj,
  features = mk,
  group.by = "seurat_clusters",
  assay = ASSAY_USE,
  # slot  = SLOT_USE,
  dot.scale = 6
) +
  RotatedAxis()

p

# ===== Fix: Check and locate cluster column =====

# First, let's check what columns exist in metadata
cat("=== Available columns in metadata ===\n")
print(colnames(seurat_obj@meta.data))

# Common cluster column names in Seurat
possible_cluster_cols <- c(
  "seurat_clusters",

  "leiden_res1.0"
)

# Find which one exists
existing_cols <- intersect(
  possible_cluster_cols,
  colnames(seurat_obj@meta.data)
)
cat("\n=== Found cluster columns ===\n")
print(existing_cols)

# If you already know the cluster column name, set it here:
# OPTION 1: Manually specify (replace "YOUR_CLUSTER_COLUMN" with actual name)
cluster_col <- "seurat_clusters" # <-- CHANGE THIS to your actual column name

# OPTION 2: Auto-detect (use first found cluster column)
# cluster_col <- existing_cols[1]

# Verify the cluster column
cat("\n=== Using cluster column:", cluster_col, "===\n")
cat("Cluster distribution:\n")
print(table(seurat_obj@meta.data[[cluster_col]]))


# ===== B Cell Annotation - Seurat V5 Direct Assignment =====

library(Seurat)
library(dplyr)

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Current cells:", ncol(seurat_obj), "\n\n")

# Use seurat_clusters column
cluster_col <- "seurat_clusters"
cat("Cluster distribution:\n")
print(table(seurat_obj@meta.data[[cluster_col]]))

# ===== Direct annotation to meta.data =====

# Level 2 mapping
level2_map <- c(
  "1" = "Naive B cells",
  "2" = "Plasma cells",
  "3" = "Activated B cells",
  "4" = "Plasma cells",
  "5" = "Memory B cells",
  "6" = "Memory B cells",
  "7" = "Activated B cells",
  "8" = "Proliferating B cells",
  "9" = "Naive B cells",
  "10" = "Memory B cells",
  "11" = "Plasma cells"
)

# Level 3 mapping
level3_map <- c(
  "1" = "CCR7+ Naive B cells",
  "2" = "Proliferating Plasma cells",
  "3" = "Activated B cells",
  "4" = "IgA+ Secretory Plasma cells",
  "5" = "MARCKSL1+ Memory B cells",
  "6" = "IGKV3+ Memory B cells",
  "7" = "Germinal Center B cells",
  "8" = "Proliferating B cells",
  "9" = "Transitional B cells",
  "10" = "LTA/LTB+ Memory B cells",
  "11" = "Long-lived Plasma cells"
)

# Apply annotations DIRECTLY to meta.data
seurat_obj@meta.data$cell_type_level_2 <- level2_map[as.character(seurat_obj@meta.data[[
  cluster_col
]])]
seurat_obj@meta.data$cell_type_level_3 <- level3_map[as.character(seurat_obj@meta.data[[
  cluster_col
]])]

# Check for NAs
if (any(is.na(seurat_obj@meta.data$cell_type_level_2))) {
  cat("\n⚠️ Found NA values!\n")
  stop("Check annotation maps.")
}

# Set factor levels
level2_order <- c(
  "Naive B cells",
  "Memory B cells",
  "Activated B cells",
  "Proliferating B cells",
  "Plasma cells"
)

level3_order <- c(
  "CCR7+ Naive B cells",
  "Transitional B cells",
  "MARCKSL1+ Memory B cells",
  "IGKV3+ Memory B cells",
  "LTA/LTB+ Memory B cells",
  "Activated B cells",
  "Germinal Center B cells",
  "Proliferating B cells",
  "IgA+ Secretory Plasma cells",
  "Proliferating Plasma cells",
  "Long-lived Plasma cells"
)

seurat_obj@meta.data$cell_type_level_2 <- factor(
  seurat_obj@meta.data$cell_type_level_2,
  levels = level2_order
)

seurat_obj@meta.data$cell_type_level_3 <- factor(
  seurat_obj@meta.data$cell_type_level_3,
  levels = level3_order
)

# Verify
cat("\n=== Level 2 ===\n")
print(table(seurat_obj@meta.data$cell_type_level_2))

cat("\n=== Level 3 ===\n")
print(table(seurat_obj@meta.data$cell_type_level_3))

# Save
saveRDS(seurat_obj, "Bcells_annotated.rds")
cat("\n✓ Saved!\n")

# ===== Plots =====

library(ggplot2)
library(cowplot)

colors_level2 <- c(
  "Naive B cells" = "#4DBBD5",
  "Memory B cells" = "#00A087",
  "Activated B cells" = "#E64B35",
  "Proliferating B cells" = "#F39B7F",
  "Plasma cells" = "#8491B4"
)

colors_level3 <- c(
  "CCR7+ Naive B cells" = "#4DBBD5",
  "Transitional B cells" = "#91D1C2",
  "MARCKSL1+ Memory B cells" = "#00A087",
  "IGKV3+ Memory B cells" = "#3C5488",
  "LTA/LTB+ Memory B cells" = "#7E6148",
  "Activated B cells" = "#E64B35",
  "Germinal Center B cells" = "#DC0000",
  "Proliferating B cells" = "#F39B7F",
  "IgA+ Secretory Plasma cells" = "#8491B4",
  "Proliferating Plasma cells" = "#B09C85",
  "Long-lived Plasma cells" = "#7876B1"
)

p1 <- DimPlot(
  seurat_obj,
  group.by = "cell_type_level_2",
  cols = colors_level2,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("Level 2") +
  theme_classic()

p2 <- DimPlot(
  seurat_obj,
  group.by = "cell_type_level_3",
  cols = colors_level3,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("Level 3") +
  theme_classic()

plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.3))
ggsave("Bcells_UMAP.pdf", width = 18, height = 7)

# Dotplot
markers <- c(
  "CCR7",
  "SELL",
  "YBX3",
  "MARCKSL1",
  "FCRL3",
  "LTA",
  "LTB",
  "CD83",
  "AICDA",
  "MME",
  "MKI67",
  "TOP2A",
  "PLK1",
  "MZB1",
  "XBP1",
  "SDC1",
  "PRDM1",
  "BPIFA1",
  "STATH",
  "LTF",
  "MERTK"
)

DotPlot(
  seurat_obj,
  features = markers,
  group.by = "cell_type_level_3",
  cols = c("lightgrey", "red"),
  dot.scale = 6
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("Bcells_dotplot.pdf", width = 14, height = 10)

# Summary
summary_level2 <- seurat_obj@meta.data %>%
  count(cell_type_level_2) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

summary_level3 <- seurat_obj@meta.data %>%
  count(cell_type_level_2, cell_type_level_3) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

write.csv(summary_level2, "summary_level2.csv", row.names = FALSE)
write.csv(summary_level3, "summary_level3.csv", row.names = FALSE)

cat("\n=== Level 2 Summary ===\n")
print(summary_level2)

cat("\n✓ Done! Files: Bcells_annotated.rds, *.pdf, *.csv\n")
saveRDS(seurat_obj, 'B_final_20260103.rds')

# ===== Comprehensive B Cell Marker Visualization =====
# Date: 2026-01-04
# Purpose: 整合所有11个B细胞亚型的特征marker

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

# seurat_obj <- readRDS("B_final_20260103.rds")

# ===== 完整Marker列表（按细胞类型分层组织）=====

comprehensive_markers <- list(
  # ========== GROUP 1: Pan B Cell Markers ==========
  `Pan B Cell` = c(
    "MS4A1", # CD20
    "CD79A", # B cell receptor component
    "CD19" # B cell co-receptor
  ),

  # ========== GROUP 2: Naive B Cells (Cluster 1) ==========
  `Naive B` = c(
    "CCR7", # Lymph node homing (1.81 log2FC)
    "SELL", # L-selectin/CD62L (1.61 log2FC)
    "CD69", # Early activation (2.10 log2FC)
    "JUNB", # AP-1 TF, early response (1.59 log2FC)
    "DUSP1" # Phosphatase, activation regulation
  ),

  # ========== GROUP 3: Transitional B Cells (Cluster 9) ==========
  `Transitional B` = c(
    "YBX3", # ⭐⭐⭐ Y-box protein (3.54 log2FC!)
    "CD55", # ⭐⭐ Complement regulator (2.14 log2FC)
    "CREM", # ⭐⭐ cAMP-responsive TF (2.21 log2FC)
    "PIK3IP1", # PI3K inhibitor (2.35 log2FC)
    "VPS37B" # ESCRT component (2.39 log2FC)
  ),

  # ========== GROUP 4: Memory B Cells (Clusters 5,6,10) ==========
  `Memory B - General` = c(
    "CD27", # Memory marker
    "MARCKSL1", # ⭐ Memory signature (1.35 log2FC)
    "FCRL3", # Fc receptor-like 3
    "RGS13", # GC/Memory marker (1.11 log2FC)
    "SEMA4A" # Semaphorin
  ),

  `Memory B - Tissue Resident` = c(
    "LTA", # ⭐ Lymphotoxin α (1.33 log2FC, Cl 10)
    "LTB", # ⭐ Lymphotoxin β (1.02 log2FC, Cl 10)
    "CD72", # Co-inhibitory receptor
    "HTR3A" # Serotonin receptor (Cl 10 specific)
  ),

  # ========== GROUP 5: Activated B Cells (Cluster 3) ==========
  `Activated B` = c(
    "CD83", # Activation marker
    "MIR155HG", # ⭐⭐⭐ Inflammatory lncRNA (2.74 log2FC!)
    "MGAT5", # ⭐⭐ N-glycosylation (2.23 log2FC)
    "ADGRE5", # ⭐⭐ Adhesion GPCR (2.45 log2FC)
    "DOCK10" # Rho GTPase regulator
  ),

  # ========== GROUP 6: Germinal Center B Cells (Cluster 7) ==========
  `GC B - Core` = c(
    "AICDA", # ⭐⭐⭐ AID enzyme (2.05 log2FC)
    "MME", # ⭐⭐⭐ CD10 (2.29 log2FC)
    "BCL6" # Master TF
  ),

  `GC B - Metabolism` = c(
    "SUGCT", # ⭐⭐ Succinate metabolism (2.59 log2FC!)
    "SLC2A5" # ⭐ Glucose transporter (2.17 log2FC)
  ),

  `GC B - Ig Editing` = c(
    "IGHV1-69", # ⭐ Heavy chain (2.17 log2FC)
    "IGKV1-33", # ⭐⭐⭐ Light chain (4.70 log2FC!)
    "IGKV1D-33" # ⭐⭐⭐ Light chain (4.66 log2FC!)
  ),

  # ========== GROUP 7: Proliferating B Cells (Cluster 8) ==========
  `Proliferating B - Mitosis` = c(
    "MKI67", # Proliferation marker
    "TOP2A", # DNA topoisomerase
    "KIF20A", # ⭐⭐⭐ Kinesin (6.22 log2FC!)
    "PLK1", # ⭐⭐⭐ Polo-like kinase (5.81 log2FC!)
    "UBE2C" # ⭐⭐⭐ Ubiquitin ligase (5.66 log2FC!)
  ),

  `Proliferating B - DNA Replication` = c(
    "PCNA", # DNA clamp
    "MCM4", # DNA helicase
    "TYMS" # Thymidylate synthase
  ),

  # ========== GROUP 8: Plasma Cells - General (All Clusters 2,4,11) ==========
  `Plasma - Core` = c(
    "MZB1", # ER chaperone
    "XBP1", # UPR transcription factor
    "JCHAIN", # J chain for IgA/IgM secretion
    "SDC1", # CD138, syndecan-1
    "PRDM1" # BLIMP-1, plasma cell master TF
  ),

  # ========== GROUP 9: Proliferating Plasma Cells (Cluster 2) ==========
  `Proliferating Plasma` = c(
    "GINS2", # ⭐⭐⭐ DNA replication (3.64 log2FC!)
    "UHRF1", # ⭐⭐⭐ Epigenetic regulator (3.45 log2FC!)
    "CLSPN", # ⭐⭐⭐ Checkpoint protein (3.41 log2FC!)
    "ASF1B", # Histone chaperone
    "CDT1" # Replication licensing
  ),

  # ========== GROUP 10: IgA+ Secretory Plasma Cells (Cluster 4) ==========
  `IgA Secretory - Antimicrobial` = c(
    "BPIFA1", # ⭐⭐⭐ BPI fold protein (4.49 log2FC!)
    "STATH", # ⭐⭐⭐ Statherin (4.54 log2FC!)
    "LTF", # ⭐⭐⭐ Lactoferrin (4.62 log2FC!)
    "LYZ", # ⭐⭐⭐ Lysozyme (4.60 log2FC!)
    "PIP" # ⭐⭐ Prolactin-induced protein (4.43 log2FC!)
  ),

  `IgA Secretory - Secretory` = c(
    "ZG16B", # ⭐⭐⭐ Zymogen granule (4.70 log2FC!)
    "C6orf58", # ⭐⭐⭐ Secretory protein (4.73 log2FC!)
    "IGHA1", # IgA heavy chain α1
    "IGHA2" # IgA heavy chain α2 (4.20 log2FC)
  ),

  # ========== GROUP 11: Long-lived Plasma Cells (Cluster 11) ==========
  `Long-lived Plasma` = c(
    "MERTK", # ⭐⭐⭐ Survival signal (5.31 log2FC!)
    "DGKI", # ⭐⭐⭐ Lipid kinase (6.32 log2FC!)
    "BMP6", # ⭐⭐⭐ Bone morphogenetic protein (4.85 log2FC!)
    "FAR2" # Fatty acyl-CoA reductase (4.63 log2FC!)
  )
)

# Flatten all markers
all_markers_flat <- unique(unlist(comprehensive_markers))
all_markers_available <- intersect(all_markers_flat, rownames(seurat_obj))

cat("=== Comprehensive Marker Check ===\n")
cat("Total unique markers:", length(all_markers_flat), "\n")
cat("Available in dataset:", length(all_markers_available), "\n")
cat(
  "Missing:",
  length(all_markers_flat) - length(all_markers_available),
  "\n\n"
)

# ===== 1. MASTER DOTPLOT - All Subtypes =====

# 创建分组label（用于dotplot分隔线）
marker_groups <- data.frame(
  Marker = unlist(comprehensive_markers),
  Group = rep(
    names(comprehensive_markers),
    sapply(comprehensive_markers, length)
  ),
  stringsAsFactors = FALSE
)

# Filter available markers
dotplot_markers_ordered <- unlist(comprehensive_markers)
dotplot_markers_ordered <- intersect(
  dotplot_markers_ordered,
  rownames(seurat_obj)
)

cat(
  "Creating master dotplot with",
  length(dotplot_markers_ordered),
  "markers...\n"
)

p_master_dotplot <- DotPlot(
  seurat_obj,
  features = dotplot_markers_ordered,
  group.by = "cell_type_level_3",
  cols = c("lightgrey", "#B2182B"), # Colorbrewer red
  dot.scale = 6
) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90", size = 0.3)
  ) +
  labs(
    title = "B Cell Subtype Comprehensive Marker Panel",
    subtitle = paste0(
      length(dotplot_markers_ordered),
      " markers across 11 subtypes"
    ),
    x = NULL,
    y = NULL
  )

ggsave(
  "Bcells_DotPlot_MASTER_comprehensive.pdf",
  p_master_dotplot,
  width = 14,
  height = 24
)

cat("✓ Saved: Bcells_DotPlot_MASTER_comprehensive.pdf\n\n")

# ===== 2. Feature Plots - Organized by Cell Type =====

# Panel 1: Pan B + Naive + Transitional
cat("Creating FeaturePlot Panel 1: Pan B / Naive / Transitional...\n")

panel1_markers <- c(
  # Pan B
  "MS4A1",
  "CD79A",
  "CD19",
  # Naive B
  "CCR7",
  "SELL",
  "CD69",
  "JUNB",
  # Transitional B
  "YBX3",
  "CD55",
  "CREM",
  "PIK3IP1"
)
panel1_markers <- intersect(panel1_markers, rownames(seurat_obj))

p_feature_panel1 <- FeaturePlot(
  seurat_obj,
  features = panel1_markers,
  ncol = 4,
  pt.size = 0.3,
  order = TRUE,
  cols = c("lightgrey", "#2166AC") # Blue
) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "Bcells_FeaturePlot_Panel1_PanB_Naive_Transitional.pdf",
  p_feature_panel1,
  width = 16,
  height = 12
)

cat("✓ Saved: Panel 1\n")

# Panel 2: Memory B Cells
cat("Creating FeaturePlot Panel 2: Memory B cells...\n")

panel2_markers <- c(
  # General memory
  "CD27",
  "MARCKSL1",
  "FCRL3",
  "RGS13",
  "SEMA4A",
  "SERPINA9",
  "NEIL1",
  # Tissue-resident
  "LTA",
  "LTB",
  "CD72",
  "HTR3A"
)
panel2_markers <- intersect(panel2_markers, rownames(seurat_obj))

p_feature_panel2 <- FeaturePlot(
  seurat_obj,
  features = panel2_markers,
  ncol = 4,
  pt.size = 0.3,
  order = TRUE,
  cols = c("lightgrey", "#1B7837") # Green
) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "Bcells_FeaturePlot_Panel2_Memory.pdf",
  p_feature_panel2,
  width = 16,
  height = 12
)

cat("✓ Saved: Panel 2\n")

# Panel 3: Activated + GC B Cells
cat("Creating FeaturePlot Panel 3: Activated / GC B cells...\n")

panel3_markers <- c(
  # Activated
  "CD83",
  "MIR155HG",
  "MGAT5",
  "ADGRE5",
  # GC core
  "AICDA",
  "MME",
  "BCL6",
  # GC metabolism
  "SUGCT",
  "SLC2A5",
  # GC Ig editing
  "IGHV1-69",
  "IGKV1-33",
  "IGKV1D-33"
)
panel3_markers <- intersect(panel3_markers, rownames(seurat_obj))

p_feature_panel3 <- FeaturePlot(
  seurat_obj,
  features = panel3_markers,
  ncol = 4,
  pt.size = 0.3,
  order = TRUE,
  cols = c("lightgrey", "#D6604D") # Orange-red
) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "Bcells_FeaturePlot_Panel3_Activated_GC.pdf",
  p_feature_panel3,
  width = 16,
  height = 12
)

cat("✓ Saved: Panel 3\n")

# Panel 4: Proliferating Cells (B + Plasma)
cat("Creating FeaturePlot Panel 4: Proliferating cells...\n")

panel4_markers <- c(
  # Proliferating B (Cluster 8)
  "MKI67",
  "TOP2A",
  "KIF20A",
  "PLK1",
  "UBE2C",
  "PCNA",
  # Proliferating Plasma (Cluster 2)
  "GINS2",
  "UHRF1",
  "CLSPN",
  "MCM4",
  "ASF1B",
  "CDT1"
)
panel4_markers <- intersect(panel4_markers, rownames(seurat_obj))

p_feature_panel4 <- FeaturePlot(
  seurat_obj,
  features = panel4_markers,
  ncol = 4,
  pt.size = 0.3,
  order = TRUE,
  cols = c("lightgrey", "#762A83") # Purple
) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "Bcells_FeaturePlot_Panel4_Proliferating.pdf",
  p_feature_panel4,
  width = 16,
  height = 12
)

cat("✓ Saved: Panel 4\n")

# Panel 5: Plasma Cells - General + IgA Secretory
cat("Creating FeaturePlot Panel 5: Plasma cells (General + IgA)...\n")

panel5_markers <- c(
  # General plasma
  "MZB1",
  "XBP1",
  "JCHAIN",
  "SDC1",
  "PRDM1",
  # IgA secretory - antimicrobial
  "BPIFA1",
  "STATH",
  "LTF",
  "LYZ",
  "PIP",
  # IgA secretory - secretory machinery
  "ZG16B",
  "C6orf58",
  "IGHA1",
  "IGHA2"
)
panel5_markers <- intersect(panel5_markers, rownames(seurat_obj))

p_feature_panel5 <- FeaturePlot(
  seurat_obj,
  features = panel5_markers,
  ncol = 4,
  pt.size = 0.3,
  order = TRUE,
  cols = c("lightgrey", "#8C510A") # Brown
) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "Bcells_FeaturePlot_Panel5_Plasma_General_IgA.pdf",
  p_feature_panel5,
  width = 16,
  height = 14
)

cat("✓ Saved: Panel 5\n")

# Panel 6: Long-lived Plasma Cells
cat("Creating FeaturePlot Panel 6: Long-lived Plasma cells...\n")

panel6_markers <- c(
  # Core plasma markers
  "MZB1",
  "XBP1",
  "SDC1",
  "PRDM1",
  # Long-lived specific
  "MERTK",
  "DGKI",
  "BMP6",
  "FAR2"
)
panel6_markers <- intersect(panel6_markers, rownames(seurat_obj))

p_feature_panel6 <- FeaturePlot(
  seurat_obj,
  features = panel6_markers,
  ncol = 4,
  pt.size = 0.3,
  order = TRUE,
  cols = c("lightgrey", "#543005") # Dark brown
) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "Bcells_FeaturePlot_Panel6_LongLived_Plasma.pdf",
  p_feature_panel6,
  width = 16,
  height = 8
)

cat("✓ Saved: Panel 6\n\n")

# ===== 3. Specialized DotPlots =====

# DotPlot 1: Key Markers Only (Simplified for presentation)
cat("Creating simplified key marker dotplot...\n")

key_markers_simple <- c(
  # Pan + Naive
  "MS4A1",
  "CCR7",
  "SELL",
  # Transitional
  "YBX3",
  "CD55",
  "CREM",
  # Memory
  "MARCKSL1",
  "FCRL3",
  "LTA",
  "LTB",
  # Activated
  "CD83",
  "MIR155HG",
  # GC
  "AICDA",
  "MME",
  "SUGCT",
  # Proliferating
  "MKI67",
  "KIF20A",
  "GINS2",
  # Plasma
  "MZB1",
  "XBP1",
  "SDC1",
  # IgA
  "BPIFA1",
  "STATH",
  "LTF",
  # Long-lived
  "MERTK",
  "PRDM1"
)
key_markers_simple <- intersect(key_markers_simple, rownames(seurat_obj))

p_dotplot_simple <- DotPlot(
  seurat_obj,
  features = key_markers_simple,
  group.by = "cell_type_level_3",
  cols = c("lightgrey", "#B2182B"),
  dot.scale = 10
) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "B Cell Subtypes: Key Signature Markers",
    subtitle = "Publication-ready simplified panel",
    x = NULL,
    y = NULL
  )

ggsave(
  "Bcells_DotPlot_KEY_simplified.pdf",
  p_dotplot_simple,
  width = 12,
  height = 10
)

cat("✓ Saved: Simplified key marker dotplot\n")

# DotPlot 2: Transitional vs Other Naive/Memory
cat("Creating Transitional-focused comparison dotplot...\n")

transitional_markers <- c(
  # Transitional specific
  "YBX3",
  "CD55",
  "CREM",
  "PIK3IP1",
  "VPS37B",
  # Should be low in Transitional
  "CCR7",
  "SELL", # Naive
  "CD27",
  "MARCKSL1", # Memory
  "MZB1" # Plasma
)
transitional_markers <- intersect(transitional_markers, rownames(seurat_obj))

p_dotplot_transitional <- DotPlot(
  seurat_obj,
  features = transitional_markers,
  group.by = "cell_type_level_3",
  cols = c("lightgrey", "#4575B4"),
  dot.scale = 12
) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Transitional B Cells: Defining Markers",
    subtitle = "Immature B cells from bone marrow",
    x = NULL,
    y = NULL
  )

ggsave(
  "Bcells_DotPlot_Transitional_focused.pdf",
  p_dotplot_transitional,
  width = 12,
  height = 8
)

cat("✓ Saved: Transitional-focused dotplot\n\n")

# ===== 4. Summary Reference Table =====

marker_reference_complete <- data.frame(
  Marker = unlist(comprehensive_markers),
  Category = rep(
    names(comprehensive_markers),
    sapply(comprehensive_markers, length)
  ),
  Available = unlist(comprehensive_markers) %in% rownames(seurat_obj),
  stringsAsFactors = FALSE
)

# Add cluster information
cluster_info <- data.frame(
  Category = names(comprehensive_markers),
  Primary_Cluster = c(
    "All",
    "1",
    "9",
    "5",
    "10",
    "3",
    "7",
    "7",
    "7",
    "8",
    "8",
    "2,4,11",
    "2",
    "4",
    "4",
    "11"
  ),
  Subtype = c(
    "Pan B",
    "CCR7+ Naive",
    "Transitional",
    "MARCKSL1+ Memory",
    "LTA/LTB+ Memory",
    "Activated",
    "GC-Core",
    "GC-Metabolism",
    "GC-IgEditing",
    "Prolif B-Mitosis",
    "Prolif B-DNA",
    "Plasma-Core",
    "Prolif Plasma",
    "IgA-Antimicrobial",
    "IgA-Secretory",
    "Long-lived Plasma"
  )
)

marker_reference_complete <- merge(
  marker_reference_complete,
  cluster_info,
  by = "Category",
  all.x = TRUE
)

write.csv(
  marker_reference_complete,
  "Bcells_marker_reference_COMPLETE.csv",
  row.names = FALSE
)

cat("✓ Saved: Bcells_marker_reference_COMPLETE.csv\n\n")

# ===== 5. Final Summary =====

cat(rep("=", 70), "\n", sep = "")
cat("✓✓✓ COMPREHENSIVE B CELL VISUALIZATION COMPLETE ✓✓✓\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("📊 Files Generated:\n\n")
cat("MASTER DOTPLOT:\n")
cat("  1. Bcells_DotPlot_MASTER_comprehensive.pdf\n")
cat(
  "     - ",
  length(dotplot_markers_ordered),
  " markers across 11 subtypes\n\n"
)

cat("FEATURE PLOT PANELS (6 panels):\n")
cat(
  "  2. Panel 1: Pan B / Naive / Transitional (",
  length(panel1_markers),
  " markers)\n"
)
cat("  3. Panel 2: Memory B cells (", length(panel2_markers), " markers)\n")
cat(
  "  4. Panel 3: Activated / GC B cells (",
  length(panel3_markers),
  " markers)\n"
)
cat(
  "  5. Panel 4: Proliferating cells (",
  length(panel4_markers),
  " markers)\n"
)
cat(
  "  6. Panel 5: Plasma (General + IgA) (",
  length(panel5_markers),
  " markers)\n"
)
cat(
  "  7. Panel 6: Long-lived Plasma (",
  length(panel6_markers),
  " markers)\n\n"
)

cat("SPECIALIZED DOTPLOTS:\n")
cat("  8. Bcells_DotPlot_KEY_simplified.pdf\n")
cat("     - Publication-ready simplified panel\n")
cat("  9. Bcells_DotPlot_Transitional_focused.pdf\n")
cat("     - Transitional B cell specific markers\n\n")

cat("REFERENCE TABLE:\n")
cat("  10. Bcells_marker_reference_COMPLETE.csv\n")
cat("      - All markers with category and cluster info\n\n")

cat("📈 Marker Statistics:\n")
cat("  - Total unique markers: ", length(all_markers_flat), "\n")
cat("  - Available in dataset: ", length(all_markers_available), "\n")
cat(
  "  - Coverage: ",
  round(length(all_markers_available) / length(all_markers_flat) * 100, 1),
  "%\n\n"
)

cat("🔬 Subtype Coverage:\n")
subtype_summary <- table(marker_reference_complete$Subtype)
for (subtype in names(subtype_summary)) {
  cat("  -", subtype, ":", subtype_summary[subtype], "markers\n")
}

cat("\n✓ All visualizations complete!\n")


# ======================================================================
# Decisive Marker Genes for B Cell Annotation Validation
# ======================================================================
# Usage: Use these markers for DotPlot and FeaturePlot to validate clusters
# Date: 2025-01-04
# ======================================================================

# ======================================================================
# TIER 1: GOLD STANDARD MARKERS (Must check these first!)
# ======================================================================

gold_standard_markers <- list(
  # Pan-B cell markers (should be positive in ALL true B cells)
  Pan_B = c(
    "MS4A1", # CD20 - mature B cells
    "CD19", # Pan-B marker
    "CD79A" # B cell receptor component
  ),

  # Germinal Center (MOST SPECIFIC!)
  Germinal_Center = c(
    "AICDA", # ⭐⭐⭐ GOLD STANDARD for GC B cells
    "BCL6", # ⭐⭐ GC transcription factor (cluster 11)
    "MME" # CD10 - GC marker
  ),

  # Plasma cells (HIGH specificity)
  Plasma = c(
    "JCHAIN", # ⭐⭐⭐ Plasma cell definitive marker
    "MZB1", # ⭐⭐ Plasma cell specific
    "PRDM1", # ⭐⭐ BLIMP-1, master regulator (cluster 21)
    "XBP1" # ⭐ ER stress, antibody production
  ),

  # IgA+ secretory (for clusters 0,1)
  IgA_Secretory = c(
    "IGHA1", # ⭐⭐ IgA heavy chain
    "IGHA2", # ⭐⭐ IgA heavy chain
    "PIGR" # ⭐ Polymeric Ig receptor (secretory)
  ),

  # Naive B cells
  Naive = c(
    "FCMR", # ⭐⭐⭐ Highly specific for naive (clusters 2,7)
    "SELL", # CD62L - lymph node homing
    "FCER2" # CD23 - naive/follicular
  ),

  # Memory B cells
  Memory = c(
    "CD27", # ⭐⭐ Classical memory marker (if available)
    "MARCKSL1" # ⭐⭐ Memory marker (cluster 13)
  ),

  # Atypical Memory (VERY IMPORTANT for CRSwNP!)
  Atypical_Memory = c(
    "FCRL4", # ⭐⭐⭐ GOLD STANDARD for atypical/exhausted (cluster 15)
    "ITGAX", # CD11c - atypical marker
    "TBX21" # T-bet - atypical marker (if available)
  ),

  # Activated B cells
  Activated = c(
    "CD69", # ⭐⭐ Early activation (cluster 6)
    "CD83", # ⭐⭐ Mature activation (cluster 4)
    "CD86" # Co-stimulation (if available)
  ),

  # Proliferating
  Proliferating = c(
    "MKI67", # ⭐⭐⭐ GOLD STANDARD for proliferation
    "TOP2A", # ⭐⭐ Mitosis
    "PCNA" # ⭐ DNA replication
  ),

  # CONTAMINATION markers (should be NEGATIVE!)
  Contamination = c(
    "LCK", # ⭐⭐⭐ T cell specific (cluster 5)
    "CD3E", # T cell marker
    "CD8A" # T cell marker
  )
)

# ======================================================================
# TIER 2: SUPPORTING MARKERS (Use for refined annotation)
# ======================================================================

supporting_markers <- list(
  # Cluster-specific from your data
  Cluster_Specific = c(
    "CCR7", # Cluster 8 (CCR7+ activated)
    "NFKB1", # Cluster 4 (NF-κB activated)
    "REL", # Cluster 4 (NF-κB family)
    "CD1C", # Cluster 16 (antigen-presenting memory)
    "HSPA1A", # Cluster 6 (stress response - potential artifact)
    "FOS", # Cluster 6 (immediate early gene)
    "RGS13" # Germinal center
  ),

  # Cell cycle markers (for clusters 14,17,18)
  Cell_Cycle = c(
    "PCNA", # S phase
    "MCM2", # S phase (if available)
    "CCNB1", # G2/M
    "CDK1" # G2/M
  ),

  # Immunoglobulin isotypes
  Ig_Isotypes = c(
    "IGHG1", # IgG1
    "IGHG4", # IgG4 (cluster 3)
    "IGHM", # IgM
    "IGKC" # Kappa light chain
  )
)

# ======================================================================
# RECOMMENDED DOTPLOT CONFIGURATION
# ======================================================================

# Option 1: Comprehensive DotPlot (all key markers)
all_decisive_markers <- c(
  # Pan-B
  "MS4A1",
  "CD19",
  "CD79A",
  # Naive
  "FCMR",
  "SELL",
  "FCER2",
  # Memory
  "CD27",
  "MARCKSL1",
  # Atypical Memory
  "FCRL4",
  "ITGAX",
  # Germinal Center
  "AICDA",
  "BCL6",
  "MME",
  # Activated
  "CD69",
  "CD83",
  "CCR7",
  "NFKB1",
  # Plasma
  "JCHAIN",
  "MZB1",
  "PRDM1",
  "XBP1",
  # IgA secretory
  "IGHA1",
  "PIGR",
  # Proliferation
  "MKI67",
  "TOP2A",
  # Contamination
  "LCK",
  "CD3E"
)

# Option 2: Minimal Essential Panel (top 15 most decisive)
essential_panel <- c(
  "MS4A1", # Pan-B
  "FCMR", # Naive
  "CD27", # Memory (if available, otherwise use MARCKSL1)
  "FCRL4", # Atypical Memory
  "AICDA", # ⭐ Germinal Center (MUST HAVE!)
  "BCL6", # Germinal Center
  "CD69", # Activated
  "CCR7", # Activated subset
  "JCHAIN", # ⭐ Plasma (MUST HAVE!)
  "MZB1", # Plasma
  "IGHA1", # IgA+ plasma
  "PRDM1", # Long-lived plasma
  "MKI67", # ⭐ Proliferating (MUST HAVE!)
  "LCK", # ⭐ Contamination check (MUST HAVE!)
  "CD3E" # Contamination check
)

# ======================================================================
# FEATUREPLOT RECOMMENDATIONS
# ======================================================================

# Set 1: Core B cell identity (check first!)
featureplot_set1 <- c("MS4A1", "CD19", "JCHAIN", "AICDA")

# Set 2: Differentiation states
featureplot_set2 <- c("FCMR", "CD27", "JCHAIN", "PRDM1")

# Set 3: Germinal center vs Plasma
featureplot_set3 <- c("AICDA", "BCL6", "MZB1", "XBP1")

# Set 4: Activation and proliferation
featureplot_set4 <- c("CD69", "CD83", "MKI67", "TOP2A")

# Set 5: Atypical and contamination check
featureplot_set5 <- c("FCRL4", "ITGAX", "LCK", "CD3E")

# Set 6: IgA secretory markers
featureplot_set6 <- c("IGHA1", "IGHA2", "PIGR", "JCHAIN")

# ======================================================================
# EXAMPLE CODE FOR PLOTTING
# ======================================================================

# Assuming Seurat object named 'seurat_obj'

# --- DotPlot (Recommended) ---
# Check which markers are available first
available_markers <- all_decisive_markers[
  all_decisive_markers %in% rownames(seurat_obj)
]

cat(
  "Available markers:",
  length(available_markers),
  "/",
  length(all_decisive_markers),
  "\n"
)
print(available_markers)

# Create DotPlot
library(Seurat)
library(ggplot2)

# Full panel
DotPlot(
  seurat_obj,
  features = available_markers,
  group.by = "seurat_clusters", # or your cluster column
  cols = c("lightgrey", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ggtitle("B Cell Decisive Markers by Cluster")

# Grouped DotPlot (better visualization)
# You need to manually create groups or use the lists above

# --- FeaturePlot ---
# Set 1: Core identity
FeaturePlot(
  seurat_obj,
  features = featureplot_set1,
  ncol = 2,
  pt.size = 0.5,
  order = TRUE
) # Plot high-expressing cells on top

# Set 2: Differentiation
FeaturePlot(
  seurat_obj,
  features = featureplot_set2,
  ncol = 2,
  pt.size = 0.5,
  order = TRUE
)

# For specific clusters of interest:
# Example: Check if cluster 10/11 are really GC B cells
FeaturePlot(
  seurat_obj,
  features = c("AICDA", "BCL6", "MME", "RGS13"),
  ncol = 2,
  cells = WhichCells(seurat_obj, idents = c("10", "11")),
  pt.size = 1
)

# Example: Check contamination in cluster 5
FeaturePlot(
  seurat_obj,
  features = c("LCK", "CD3E", "MS4A1", "CD19"),
  ncol = 2,
  cells = WhichCells(seurat_obj, idents = "5"),
  pt.size = 1
)

# ======================================================================
# VALIDATION CHECKLIST
# ======================================================================

cat("\n=== VALIDATION CHECKLIST ===\n\n")

cat("□ Step 1: Check Pan-B markers (MS4A1, CD19, CD79A)\n")
cat("   → Should be positive in clusters: 0-4, 6-8, 10-11, 13-19, 21\n")
cat("   → Should be NEGATIVE in: 5, 9, 12, 20\n\n")

cat("□ Step 2: Verify CONTAMINATION (LCK, CD3E)\n")
cat("   → Cluster 5 should be LCK+ (REMOVE)\n")
cat("   → Clusters 9, 12, 20 should lack B cell markers (REMOVE)\n\n")

cat("□ Step 3: Check GERMINAL CENTER (AICDA, BCL6)\n")
cat("   → Clusters 10, 11 should be AICDA+++ (very high)\n")
cat("   → Cluster 11 should be BCL6+ (Light zone)\n\n")

cat("□ Step 4: Verify PLASMA cells (JCHAIN, MZB1, PRDM1)\n")
cat("   → Clusters 0, 1, 3, 21 should be JCHAIN+\n")
cat("   → Cluster 21 should be PRDM1+ (long-lived)\n")
cat("   → Clusters 0, 1 should be IGHA1/2+ (secretory IgA)\n\n")

cat("□ Step 5: Check NAIVE markers (FCMR, SELL)\n")
cat("   → Clusters 2, 7, 19 should be FCMR+\n\n")

cat("□ Step 6: Verify ATYPICAL MEMORY (FCRL4)\n")
cat("   → Cluster 15 should be FCRL4+ (important for CRSwNP!)\n\n")

cat("□ Step 7: Check PROLIFERATION (MKI67, TOP2A)\n")
cat("   → Clusters 14, 17, 18 should be MKI67+++\n")
cat("   → Also check cluster 3 (proliferating plasma)\n\n")

cat("□ Step 8: Validate ACTIVATION (CD69, CD83, CCR7)\n")
cat("   → Cluster 4 should be CD83+ (NF-κB)\n")
cat("   → Cluster 6 should be CD69+ (but check for stress!)\n")
cat("   → Cluster 8 should be CCR7+\n\n")

# ======================================================================
# CLUSTER-SPECIFIC VALIDATION MARKERS
# ======================================================================

cluster_specific_validation <- list(
  "0" = c("IGHA1", "IGHA2", "JCHAIN", "PIGR"), # IgA secretory
  "1" = c("IGHA1", "IGHA2", "JCHAIN", "PIGR"), # IgA secretory
  "2" = c("FCMR", "SELL", "FCER2", "IL4R"), # Naive
  "3" = c("JCHAIN", "MKI67", "PCNA", "IGHG4"), # Proliferating plasma
  "4" = c("CD83", "NFKB1", "REL", "CD40"), # NF-κB activated
  "5" = c("LCK", "CD3E", "CD8A"), # REMOVE - T cell
  "6" = c("CD69", "FOS", "HSPA1A", "JUN"), # Stress/activated
  "7" = c("FCMR", "SELL", "FCER2"), # Naive/transitional
  "8" = c("CCR7", "CD69", "CXCR4"), # CCR7+ activated
  "9" = c("MERTK", "PRKCA"), # REMOVE - non-B
  "10" = c("AICDA", "MME", "BCL6"), # GC dark zone
  "11" = c("BCL6", "AICDA", "MME", "RGS13"), # GC light zone
  "12" = c("DGKI", "ESR1"), # REMOVE - non-B
  "13" = c("MARCKSL1", "CD27", "TNFRSF13B"), # Memory
  "14" = c("MKI67", "TOP2A", "CDK1"), # Cycling G2/M
  "15" = c("FCRL4", "ITGAX", "TBX21", "CCR1"), # Atypical memory
  "16" = c("CD1C", "CD27", "MARCKSL1"), # Antigen-presenting memory
  "17" = c("MKI67", "PCNA", "CCNB1"), # Cycling S/G2/M
  "18" = c("PCNA", "RRM2", "TK1"), # Cycling S phase
  "19" = c("IGKV1-33", "FCMR", "FCER2"), # Naive
  "20" = c(), # REMOVE - low quality
  "21" = c("PRDM1", "MZB1", "XBP1", "IRF4") # Long-lived plasma
)

cat("\n=== Expected Marker Expression by Cluster ===\n")
for (cluster in names(cluster_specific_validation)) {
  markers <- cluster_specific_validation[[cluster]]
  if (length(markers) > 0) {
    cat(sprintf("Cluster %s: %s\n", cluster, paste(markers, collapse = ", ")))
  }
}

# ======================================================================
# FINAL RECOMMENDATION
# ======================================================================

cat("\n=== RECOMMENDED PLOTTING ORDER ===\n\n")
cat("1. Start with DotPlot using 'all_decisive_markers'\n")
cat("   → This gives you overview of all clusters\n\n")

cat("2. Then do FeaturePlots in this order:\n")
cat("   a) featureplot_set1 - Core B cell identity\n")
cat("   b) featureplot_set3 - GC vs Plasma (most important!)\n")
cat("   c) featureplot_set5 - Contamination check\n")
cat("   d) featureplot_set4 - Activation/proliferation\n\n")

cat("3. Focus on these CRITICAL validations:\n")
cat("   ⭐ AICDA expression in clusters 10,11 (MUST be very high)\n")
cat("   ⭐ JCHAIN expression in clusters 0,1,3,21 (MUST be high)\n")
cat("   ⭐ LCK expression in cluster 5 (should be positive → REMOVE)\n")
cat("   ⭐ FCRL4 expression in cluster 15 (CRSwNP relevant!)\n")
cat("   ⭐ MKI67 in clusters 14,17,18 (cycling cells)\n\n")

cat("4. After validation, remove clusters: 5, 9, 12, 20\n")
cat("5. Consider re-clustering after QC\n\n")

# ======================================================================

# ===== B Cell Annotation and Visualization =====
# Purpose: Annotate B cell clusters and generate UMAP + Dotplot
# Author: r2end
# Date: 2025-01-05

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ===== Configuration =====
# INPUT_H5AD <- "/home/h2048/data/R/0103/ciliated/ciliated_bbknn_integrated.h5ad"  # Update path
OUTPUT_DIR <- "/home/h2048/data/R/0104/B"

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===== 1. Load Data =====
# Option A: If you have .rds file from Seurat
# seurat_obj <- readRDS(INPUT_RDS)

# Option B: Convert from h5ad using SeuratDisk (recommended)
library(SeuratDisk)
# First convert h5ad to h5seurat in Python or use Convert function
# For this workflow, assume you have a Seurat object

# Placeholder: assume seurat_obj is loaded
# For demonstration, I'll show the annotation and plotting code
# You need to load your actual Seurat object

cat("Loading Seurat object...\n")
# seurat_obj <- readRDS("path_to_seurat_object.rds")  # Update this

# ===== 2. Define Annotations =====
# Level 2 annotation (Major lineages)
annotation_level2 <- c(
  '0' = 'Naive B',
  '1' = 'GC B',
  '2' = 'GC B',
  '3' = 'Memory B',
  '4' = 'Memory B',
  '5' = 'Memory B',
  '6' = 'Plasma',
  '7' = 'Plasma',
  '8' = 'Follicular B',
  '9' = 'Memory B',
  '10' = 'Plasma'
)

# Level 3 annotation (Detailed subtypes)
annotation_level3 <- c(
  '0' = 'Naive B',
  '1' = 'GC B',
  '2' = 'Proliferating GC B',
  '3' = 'Activated Memory B',
  '4' = 'Atypical Memory B (FCRL4+)',
  '5' = 'Activated Memory B',
  '6' = 'IgA+ Plasma',
  '7' = 'IgA+ Plasma',
  '8' = 'Tissue-resident Follicular B',
  '9' = 'Tissue-resident Memory B',
  '10' = 'IgG+ Plasma'
)

# ✅ Apply annotations - CORRECT WAY
seurat_obj@meta.data$cell_type_level_2 <- annotation_level2[as.character(
  seurat_obj@meta.data$leiden_res1.0
)]
seurat_obj@meta.data$cell_type_level_3 <- annotation_level3[as.character(
  seurat_obj@meta.data$leiden_res1.0
)]

# Set factor levels for proper ordering
level2_order <- c('Naive B', 'GC B', 'Memory B', 'Follicular B', 'Plasma')
level3_order <- c(
  'Naive B',
  'GC B',
  'Proliferating GC B',
  'Activated Memory B',
  'Atypical Memory B (FCRL4+)',
  'Tissue-resident Memory B',
  'Tissue-resident Follicular B',
  'IgA+ Plasma',
  'IgG+ Plasma'
)

seurat_obj@meta.data$cell_type_level_2 <- factor(
  seurat_obj@meta.data$cell_type_level_2,
  levels = level2_order
)
seurat_obj@meta.data$cell_type_level_3 <- factor(
  seurat_obj@meta.data$cell_type_level_3,
  levels = level3_order
)

cat("Annotations applied successfully\n")
cat("Cell type distribution (Level 2):\n")
print(table(seurat_obj@meta.data$cell_type_level_2))
cat("\nCell type distribution (Level 3):\n")
print(table(seurat_obj@meta.data$cell_type_level_3))

# dev.off()

# ===== 3. UMAP Visualization =====
cat("\nGenerating UMAP plots...\n")

# Define color palettes
colors_level2 <- c(
  'Naive B' = '#E64B35',
  'GC B' = '#4DBBD5',
  'Memory B' = '#00A087',
  'Follicular B' = '#3C5488',
  'Plasma' = '#F39B7F'
)

colors_level3 <- c(
  'Naive B' = '#E64B35',
  'GC B' = '#4DBBD5',
  'Proliferating GC B' = '#91D1C2',
  'Activated Memory B' = '#00A087',
  'Atypical Memory B (FCRL4+)' = '#8491B4',
  'Tissue-resident Memory B' = '#7E6148',
  'Tissue-resident Follicular B' = '#3C5488',
  'IgA+ Plasma' = '#F39B7F',
  'IgG+ Plasma' = '#DC0000'
)

# UMAP - Level 2
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type_level_2",
  cols = colors_level2,
  pt.size = 0.5,
  label = TRUE,
  label.size = 4,
  label.box = TRUE,
  repel = TRUE
) +
  ggtitle("B Cell Subtypes (Level 2)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# UMAP - Level 3
p2 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type_level_3",
  cols = colors_level3,
  pt.size = 0.5,
  label = TRUE,
  label.size = 3.5,
  label.box = TRUE,
  repel = TRUE
) +
  ggtitle("B Cell Subtypes (Level 3)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Combined UMAP plot
combined_umap <- p1 / p2 + plot_layout(heights = c(1, 1))

# Save UMAP
ggsave(
  filename = file.path(OUTPUT_DIR, "bcell_umap_annotation.pdf"),
  plot = combined_umap,
  width = 14,
  height = 16,
  dpi = 300
)

cat("✅ UMAP saved: bcell_umap_annotation.pdf\n")

# ===== 4. Dotplot - Key Markers =====
cat("\nGenerating Dotplot...\n")

# Define key markers for each subtype
marker_genes <- c(
  # Naive B
  'S1PR1',
  'SELL',
  'IGHD',
  'IGHM',
  'TCL1A',

  # GC B
  'AICDA',
  'BCL6',
  'RGS13',
  'MME',

  # Proliferation
  'MKI67',
  'TOP2A',

  # Memory B (general)
  'CD27',
  'BCL2',

  # Activated Memory
  'NR4A1',
  'CCR7',
  'GPR183',

  # Atypical Memory
  'FCRL4',
  'FCRL5',
  'TNFRSF13B',

  # Tissue-resident
  'CD69',
  'CXCR4',
  'ITGB2',

  # Plasma
  'PRDM1',
  'XBP1',
  'MZB1',
  'JCHAIN',
  'IGHA1',
  'IGHA2',
  'IGHG1',
  'IGHG3'
)

# Filter markers that exist in the dataset
available_markers <- marker_genes[marker_genes %in% rownames(seurat_obj)]

if (length(available_markers) < length(marker_genes)) {
  missing <- setdiff(marker_genes, available_markers)
  cat("⚠️ Missing markers:", paste(missing, collapse = ", "), "\n")
}

# Dotplot by Level 3
p_dot <- DotPlot(
  seurat_obj,
  features = available_markers,
  group.by = "cell_type_level_3",
  cols = c("lightgrey", "red"),
  dot.scale = 8
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  ggtitle("B Cell Subtype Markers")

# Save Dotplot
ggsave(
  filename = file.path(OUTPUT_DIR, "bcell_dotplot_markers.pdf"),
  plot = p_dot,
  width = 12,
  height = 10,
  dpi = 300
)

cat("✅ Dotplot saved: bcell_dotplot_markers.pdf\n")

# ===== 5. Alternative: Detailed Dotplot by Category =====
cat("\nGenerating categorized Dotplot...\n")

# Organize markers by functional category
marker_list <- list(
  'Naive' = c('S1PR1', 'SELL', 'IGHD', 'IGHM', 'TCL1A', 'FCER2'),
  'GC & Proliferation' = c('AICDA', 'BCL6', 'RGS13', 'MME', 'MKI67', 'TOP2A'),
  'Memory' = c('CD27', 'BCL2', 'NR4A1', 'CCR7', 'GPR183'),
  'Atypical & Resident' = c('FCRL4', 'FCRL5', 'CD69', 'CXCR4', 'ITGB2'),
  'Plasma' = c('PRDM1', 'XBP1', 'MZB1', 'JCHAIN', 'IGHA1', 'IGHG1')
)

# Filter available markers
marker_list_filtered <- lapply(marker_list, function(x) {
  x[x %in% rownames(seurat_obj)]
})
all_markers_ordered <- unlist(marker_list_filtered)

# Create dotplot
p_dot_detailed <- DotPlot(
  seurat_obj,
  features = all_markers_ordered,
  group.by = "cell_type_level_3",
  cols = c("lightgrey", "red"),
  dot.scale = 8
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  ggtitle("B Cell Markers by Functional Category")

# Add category separators (manual annotation)
# This requires adding geom_hline, can be done manually or with annotation

ggsave(
  filename = file.path(OUTPUT_DIR, "bcell_dotplot_detailed.pdf"),
  plot = p_dot_detailed,
  width = 14,
  height = 12,
  dpi = 300
)

cat("✅ Detailed Dotplot saved: bcell_dotplot_detailed.pdf\n")

# ===== 6. Feature Plots for Key Markers =====
cat("\nGenerating feature plots for tissue-resident markers...\n")

resident_markers <- c('CD69', 'CXCR4', 'FCRL4', 'ITGB2')
resident_markers <- resident_markers[resident_markers %in% rownames(seurat_obj)]

if (length(resident_markers) > 0) {
  p_features <- FeaturePlot(
    seurat_obj,
    features = resident_markers,
    reduction = "umap",
    pt.size = 0.3,
    ncol = 2,
    cols = c("lightgrey", "red")
  ) +
    plot_annotation(
      title = "Tissue-Resident B Cell Markers",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )
    )

  ggsave(
    filename = file.path(OUTPUT_DIR, "bcell_resident_markers_featureplot.pdf"),
    plot = p_features,
    width = 12,
    height = 12,
    dpi = 300
  )

  cat("✅ Feature plots saved: bcell_resident_markers_featureplot.pdf\n")
}

# ===== 7. Save Annotated Object =====
cat("\nSaving annotated Seurat object...\n")

saveRDS(
  seurat_obj,
  file = file.path(OUTPUT_DIR, "bcell_annotated.rds")
)

cat("✅ Annotated object saved: bcell_annotated.rds\n")

# ===== 8. Summary Statistics =====
cat("\n" %+% paste(rep("=", 70), collapse = "") %+% "\n")
cat("ANNOTATION SUMMARY\n")
cat(paste(rep("=", 70), collapse = "") %+% "\n")

cat("\nCell counts by Level 2:\n")
print(table(seurat_obj$cell_type_level_2))

cat("\nCell counts by Level 3:\n")
print(table(seurat_obj$cell_type_level_3))

cat("\nCell proportions by Level 2:\n")
print(prop.table(table(seurat_obj$cell_type_level_2)) * 100)

cat("\n" %+% paste(rep("=", 70), collapse = "") %+% "\n")

cat("\n📁 All outputs saved to:", OUTPUT_DIR, "\n")
cat("Files generated:\n")
cat("  1. bcell_umap_annotation.pdf\n")
cat("  2. bcell_dotplot_markers.pdf\n")
cat("  3. bcell_dotplot_detailed.pdf\n")
cat("  4. bcell_resident_markers_featureplot.pdf\n")
cat("  5. bcell_annotated.rds\n")
cat("\n✅ Analysis completed successfully!\n")

# ===== 9. Session Info =====
cat("\nSession Info:\n")
sessionInfo()
GetH5ad(seurat_obj, 'b_final_20260104.h5ad')
