#!/usr/bin/env Rscript
# =============================================================================
# RNA_snn_res.1 Cluster QC + Lineage Validation (Epithelial focus) v1.0
# Output:
#  - cluster_qc_summary.csv
#  - cluster_signature_summary.csv
#  - cluster_problem_report.txt
#  - qc_violin.pdf
#  - marker_dotplot_RNA_snn_res.1.pdf
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


cnt <- get_counts_matrix(seurat_obj)

gene_ncells <- Matrix::rowSums(cnt > 0)
keep_genes <- names(gene_ncells[gene_ncells >= 3])

before_genes <- nrow(seurat_obj)
seurat_obj <- subset(seurat_obj, features = keep_genes)
after_genes <- nrow(seurat_obj)
before_genes
after_genes

seurat_obj <- NormalizeData(seurat_obj) #归一化

# ===== 1) 定义“组蛋白基因”集合（覆盖 HIST* + H1/H2A/H2B/H3/H4 系列）=====
hist_genes <- grep(
  "^(HIST\\d|H1-|H1F|H2A|H2B|H3|H4)",
  rownames(seurat_obj),
  value = TRUE
)

cat(sprintf("Histone-like genes detected: %d\n", length(hist_genes)))
print(head(hist_genes, 30))

# ===== 2) 重新选 HVG，并把组蛋白从 HVG 中剔除 =====
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 3000
)
VariableFeatures(seurat_obj) <- setdiff(
  VariableFeatures(seurat_obj),
  hist_genes
)
seurat_obj <- ScaleData(seurat_obj) #标准化

# 使用高变基因进行主成分分析，降低数据维度
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
seurat_obj <- RunHarmony(
  object = seurat_obj,
  group.by.vars = c("sample"),
  theta = c(6.5), # Higher theta for more diverse clustering
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
  n.neighbors = 30,
  # n.trees = 500,
  min.dist = 0.4,
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
  k.param = 30
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
# saveRDS(seurat_obj, 'ciliated_filtered_20260103.rds')
GetH5ad(seurat_obj, 'ciliated_filtered_20260104.h5ad')
# Find all markers using Wilcoxon test
library(Seurat)
library(future)

# Enable parallel processing
plan("multisession", workers = 16) # Adjust based on your CPU cores
options(future.globals.maxSize = 4000 * 1024^2) # 8GB max object size

# =========================
# 0) 组蛋白基因规则（按你的数据：H2AC*, H2BC*, H4C*, H3-3A 等）
#    - ^HIST   : HIST1H*, HIST2H* 等
#    - ^H[1-4][A-Z] : H2AC18/H2AFZ/H4C3/H1FX 等（digit 后紧跟字母）
#    - ^H3-    : H3-3A/H3-3B
#    - ^H1-    : H1-4 等
#    这样不会误伤 mouse 的 H2-K1 (digit 后是 '-')
# =========================
histone_regex <- "^(HIST|H[1-4][A-Z]|H3-|H1-)"

histone_genes <- grep(histone_regex, rownames(seurat_obj), value = TRUE)
cat(sprintf("[Info] Histone genes detected: %d\n", length(histone_genes)))

features_use <- setdiff(rownames(seurat_obj), histone_genes)
cat(sprintf("[Info] Features used for DE: %d\n", length(features_use)))

# =========================
# 1) FindAllMarkers（Wilcoxon）阶段排除组蛋白
# =========================
all_markers <- FindAllMarkers(
  object = seurat_obj,
  assay = "RNA", # 如果你用 SCT，就改成 "SCT" 并先 PrepSCTFindMarkers
  test.use = "wilcox",
  features = features_use,
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)

# =========================
# 2) 兜底：再过滤一次，保证结果里绝对没有组蛋白
# =========================
all_markers <- subset(all_markers, !(gene %in% histone_genes))

# 可选：快速自检
stopifnot(!any(all_markers$gene %in% histone_genes))


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

# -------------------------
# Config
# -------------------------
CLUSTER_COL <- "RNA_snn_res.1"
OUTDIR <- "RNA_snn_res2_validation_out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Choose assay (prefer SCT if present)
assay_use <- if ("SCT" %in% names(seurat_obj@assays)) "SCT" else "RNA"
DefaultAssay(seurat_obj) <- assay_use

# Choose reduction (prefer umap)
reduction_use <- if ("umap" %in% names(seurat_obj@reductions)) {
  "umap"
} else {
  names(seurat_obj@reductions)[1]
}
message(sprintf(
  "Using assay: %s | reduction: %s | cluster: %s",
  assay_use,
  reduction_use,
  CLUSTER_COL
))

# Ensure clustering column exists
stopifnot(CLUSTER_COL %in% colnames(seurat_obj@meta.data))
Idents(seurat_obj) <- seurat_obj[[CLUSTER_COL, drop = TRUE]]

# Ensure normalized data exists (light check)
# If you know you already have normalized data, you can comment this block.
try(
  {
    dat <- GetAssayData(seurat_obj, slot = "data")
    if (ncol(dat) == 0 || nrow(dat) == 0) stop("Empty data slot")
  },
  silent = TRUE
)

# -------------------------
# Helpers
# -------------------------
present_genes <- function(obj, genes) {
  genes <- unique(genes)
  genes[genes %in% rownames(obj)]
}

add_named_module_score <- function(obj, genes, name, assay = NULL) {
  if (!is.null(assay)) {
    DefaultAssay(obj) <- assay
  }
  genes_use <- present_genes(obj, genes)
  if (length(genes_use) < 3) {
    obj[[paste0("sig_", name)]] <- 0
    return(obj)
  }
  obj <- AddModuleScore(
    obj,
    features = list(genes_use),
    name = paste0("sig_", name, "_")
  )
  score_col <- paste0("sig_", name, "_1")
  obj[[paste0("sig_", name)]] <- obj[[score_col, drop = TRUE]]
  obj[[score_col]] <- NULL
  return(obj)
}

pct_expr <- function(obj, gene, group_col) {
  if (!gene %in% rownames(obj)) {
    return(NULL)
  }
  v <- FetchData(obj, vars = gene)[, 1]
  df <- data.frame(cluster = obj[[group_col, drop = TRUE]], expr = v > 0)
  df %>%
    group_by(cluster) %>%
    summarize(pct = mean(expr) * 100, .groups = "drop") %>%
    rename(!!paste0(gene, "_pct") := pct)
}

robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  madv <- mad(x, na.rm = TRUE)
  if (madv == 0) {
    return(rep(0, length(x)))
  }
  0.6745 * (x - med) / madv
}

# -------------------------
# Marker panels (edit freely)
# -------------------------
markers <- list(
  # Ciliated lineage
  Ciliated_Mature = c(
    "FOXJ1",
    "TPPP3",
    "PIFO",
    "HYDIN",
    "SPEF2",
    "DNAH5",
    "DNAH9",
    "DNAI1",
    "DNAAF1",
    "RSPH1",
    "RSPH10B2",
    "CFAP46",
    "CFAP54",
    "CFAP99",
    "CCDC114",
    "ODAD1",
    "WDR90"
  ),
  Ciliated_Deuterosomal = c("MCIDAS", "DEUP1", "CCNO", "CDC20B", "MYB", "PLK4"),
  Cycling = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "HMGB2",
    "ANLN",
    "DLGAP5",
    "MELK",
    "CDC20"
  ),

  # Epithelial programs
  Basal = c("KRT5", "KRT14", "TP63", "KRT15", "LGALS3", "ITGA6"),
  Secretory_Club = c(
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "KRT8",
    "KRT18",
    "BPIFA1"
  ),
  Goblet = c("MUC5AC", "SPDEF", "AGR2", "CLCA1", "BPIFB1", "TFF3"),
  Ionocyte = c("FOXI1", "CFTR", "ASCL3", "ATP6V1B1", "SLC26A4"),
  AT2 = c("SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NKX2-1", "GPR116"),

  # Non-epithelial / contamination checks
  Immune_Pan = c(
    "PTPRC",
    "LST1",
    "TYROBP",
    "FCER1G",
    "NKG7",
    "TRAC",
    "CD3D",
    "MS4A1",
    "CD79A"
  ),
  Tcell = c("TRAC", "CD3D", "CD3E", "IL7R", "LTB", "CCL5", "GZMA"),
  Bcell = c("MS4A1", "CD79A", "CD74", "CD37", "HLA-DRA"),
  Plasma = c("MZB1", "JCHAIN", "XBP1", "IGHG1"),
  Myeloid = c("LYZ", "S100A8", "S100A9", "FCN1", "C1QA", "C1QB", "C1QC"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),

  # Technical / stress flags
  Stress_AP1 = c(
    "JUN",
    "FOS",
    "FOSB",
    "JUNB",
    "ATF3",
    "DDIT3",
    "HSPA1A",
    "HSPA1B"
  ),
  IFN = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1"),
  Histone = c("H2AC18", "H2AC19", "HIST1H1C", "HIST1H2BD", "HIST2H2BE", "H1FX")
)

# Flatten for dotplot
dot_features <- unique(unlist(markers))
dot_features <- present_genes(seurat_obj, dot_features)

# -------------------------
# QC metrics
# -------------------------
# MT percent (human)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Hemoglobin percent (RBC)
hb_genes <- present_genes(
  seurat_obj,
  c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBM", "HBQ1")
)
if (length(hb_genes) > 0) {
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(
    seurat_obj,
    features = hb_genes
  )
} else {
  seurat_obj[["percent.hb"]] <- 0
}

# Cell cycle scores (optional; harmless if genes absent)
cc_s <- present_genes(seurat_obj, Seurat::cc.genes.updated.2019$s.genes)
cc_g2m <- present_genes(seurat_obj, Seurat::cc.genes.updated.2019$g2m.genes)
if (length(cc_s) >= 10 && length(cc_g2m) >= 10) {
  seurat_obj <- CellCycleScoring(
    seurat_obj,
    s.features = cc_s,
    g2m.features = cc_g2m,
    set.ident = FALSE
  )
} else {
  seurat_obj$S.Score <- 0
  seurat_obj$G2M.Score <- 0
}

# -------------------------
# Signature/module scores
# -------------------------
sig_sets <- list(
  Ciliated_Mature = markers$Ciliated_Mature,
  Ciliated_Deuterosomal = markers$Ciliated_Deuterosomal,
  Cycling = markers$Cycling,
  Basal = markers$Basal,
  Secretory_Club = markers$Secretory_Club,
  Goblet = markers$Goblet,
  Ionocyte = markers$Ionocyte,
  AT2 = markers$AT2,
  Immune = markers$Immune_Pan,
  Tcell = markers$Tcell,
  Bcell = markers$Bcell,
  Plasma = markers$Plasma,
  Myeloid = markers$Myeloid,
  Endothelial = markers$Endothelial,
  Fibroblast = markers$Fibroblast,
  Stress_AP1 = markers$Stress_AP1,
  IFN = markers$IFN,
  Histone = markers$Histone
)

for (nm in names(sig_sets)) {
  seurat_obj <- add_named_module_score(
    seurat_obj,
    sig_sets[[nm]],
    nm,
    assay = assay_use
  )
}

sig_cols <- paste0("sig_", names(sig_sets))
sig_cols <- sig_cols[sig_cols %in% colnames(seurat_obj@meta.data)]

# -------------------------
# Cluster-level summaries
# -------------------------
meta <- seurat_obj@meta.data %>%
  mutate(cluster = .data[[CLUSTER_COL]])

qc_summary <- meta %>%
  group_by(cluster) %>%
  summarize(
    n_cells = n(),
    nCount_median = median(nCount_RNA, na.rm = TRUE),
    nFeature_median = median(nFeature_RNA, na.rm = TRUE),
    percent_mt_median = median(percent.mt, na.rm = TRUE),
    percent_hb_median = median(percent.hb, na.rm = TRUE),
    S_score_median = median(S.Score, na.rm = TRUE),
    G2M_score_median = median(G2M.Score, na.rm = TRUE),
    .groups = "drop"
  )

# Sentinel gene expression % (robust contamination flags)
sentinel_genes <- c(
  "EPCAM",
  "KRT5",
  "FOXJ1",
  "CDC20B",
  "MUC5AC",
  "SCGB1A1",
  "PTPRC",
  "CD3D",
  "MS4A1",
  "LYZ",
  "HBB",
  "PECAM1",
  "COL1A1",
  "SFTPC",
  "BEST4"
)
sentinel_genes <- sentinel_genes[sentinel_genes %in% rownames(seurat_obj)]

sent_tables <- lapply(sentinel_genes, function(g) {
  pct_expr(seurat_obj, g, CLUSTER_COL)
})
sent_tables <- sent_tables[!sapply(sent_tables, is.null)]
if (length(sent_tables) > 0) {
  sent_df <- Reduce(function(x, y) full_join(x, y, by = "cluster"), sent_tables)
  qc_summary <- qc_summary %>% left_join(sent_df, by = "cluster")
}

# Signature summary by cluster (median)
sig_summary <- meta %>%
  group_by(cluster) %>%
  summarize(
    across(all_of(sig_cols), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Save CSVs
write.csv(
  qc_summary,
  file = file.path(OUTDIR, "cluster_qc_summary.csv"),
  row.names = FALSE
)
write.csv(
  sig_summary,
  file = file.path(OUTDIR, "cluster_signature_summary.csv"),
  row.names = FALSE
)

# -------------------------
# Problem flags (robust, cluster-level)
# -------------------------
flag_df <- qc_summary

# metrics to flag (robust z > 3)
flag_metrics <- c(
  "nCount_median",
  "nFeature_median",
  "percent_mt_median",
  "percent_hb_median",
  "PTPRC_pct",
  "HBB_pct",
  "SFTPC_pct",
  "COL1A1_pct",
  "PECAM1_pct"
)
flag_metrics <- flag_metrics[flag_metrics %in% colnames(flag_df)]

for (m in flag_metrics) {
  z <- robust_z(flag_df[[m]])
  flag_df[[paste0(m, "_rz")]] <- z
  flag_df[[paste0("flag_", m)]] <- z > 3
}

# Mixed-lineage heuristic (rank-based): epithelial vs immune co-high
# epithelial signature aggregate
epi_sig_cols <- intersect(
  c(
    "sig_Basal",
    "sig_Secretory_Club",
    "sig_Goblet",
    "sig_Ciliated_Mature",
    "sig_Ciliated_Deuterosomal",
    "sig_AT2"
  ),
  colnames(sig_summary)
)
imm_sig_cols <- intersect(
  c("sig_Immune", "sig_Tcell", "sig_Bcell", "sig_Myeloid", "sig_Plasma"),
  colnames(sig_summary)
)

mix_note <- sig_summary %>%
  mutate(
    epi_score = if (length(epi_sig_cols) > 0) {
      rowMeans(across(all_of(epi_sig_cols)), na.rm = TRUE)
    } else {
      0
    },
    imm_score = if (length(imm_sig_cols) > 0) {
      rowMeans(across(all_of(imm_sig_cols)), na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  mutate(
    epi_rank = rank(-epi_score, ties.method = "min"),
    imm_rank = rank(-imm_score, ties.method = "min"),
    possible_doublet_or_mixed = (epi_rank <= 3 & imm_rank <= 3)
  ) %>%
  select(
    cluster,
    epi_score,
    imm_score,
    epi_rank,
    imm_rank,
    possible_doublet_or_mixed
  )

flag_df <- flag_df %>% left_join(mix_note, by = "cluster")

# Write text report
report_path <- file.path(OUTDIR, "cluster_problem_report.txt")
con <- file(report_path, open = "wt")
writeLines("Cluster Problem Report (RNA_snn_res.1)\n", con)

writeLines("1) Robust outlier flags (robust z > 3)\n", con)
for (i in seq_len(nrow(flag_df))) {
  cl <- as.character(flag_df$cluster[i])
  hits <- c()
  for (m in flag_metrics) {
    if (isTRUE(flag_df[[paste0("flag_", m)]][i])) hits <- c(hits, m)
  }
  if (length(hits) > 0) {
    writeLines(
      sprintf("- Cluster %s: outliers -> %s", cl, paste(hits, collapse = ", ")),
      con
    )
  }
}

writeLines(
  "\n2) Possible mixed-lineage / doublet-like clusters (top3 epithelial AND top3 immune)\n",
  con
)
dbl <- flag_df %>% filter(possible_doublet_or_mixed)
if (nrow(dbl) == 0) {
  writeLines("- None flagged by rank-based rule.", con)
} else {
  apply(dbl, 1, function(r) {
    writeLines(
      sprintf(
        "- Cluster %s: epi_score=%.3f (rank %d), imm_score=%.3f (rank %d)",
        r[["cluster"]],
        as.numeric(r[["epi_score"]]),
        as.integer(r[["epi_rank"]]),
        as.numeric(r[["imm_score"]]),
        as.integer(r[["imm_rank"]])
      ),
      con
    )
  })
}

close(con)

# -------------------------
# Plots
# -------------------------
# QC violin (no points)
p_qc <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"),
  group.by = CLUSTER_COL,
  pt.size = 0
) +
  NoLegend()

pdf(file.path(OUTDIR, "qc_violin.pdf"), width = 14, height = 8)
print(p_qc)
dev.off()

# DotPlot of markers
p_dot <- DotPlot(
  seurat_obj,
  features = dot_features,
  group.by = CLUSTER_COL
) +
  RotatedAxis() +
  ggtitle("Marker DotPlot (RNA_snn_res.1)")

pdf(
  file.path(OUTDIR, "marker_dotplot_RNA_snn_res.1.pdf"),
  width = 18,
  height = 10
)
print(p_dot)
dev.off()

# Signature UMAPs (module scores)
# If too many panels, split automatically
sig_cols_plot <- sig_cols[sig_cols %in% colnames(seurat_obj@meta.data)]
plots <- lapply(sig_cols_plot, function(sc) {
  FeaturePlot(
    seurat_obj,
    features = sc,
    reduction = reduction_use,
    raster = FALSE
  ) +
    ggtitle(sc)
})
p_umap <- wrap_plots(plots, ncol = 4)

pdf(file.path(OUTDIR, "signature_umap.pdf"), width = 18, height = 14)
print(p_umap)
dev.off()

message("Done. Outputs written to: ", OUTDIR)


DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- seurat_obj[["RNA_snn_res.1", drop = TRUE]]

genes <- c(
  "EPCAM",
  "KRT5",
  "FOXJ1",
  "TPPP3",
  "CDC20B",
  "DEUP1",
  "MCIDAS",
  "MUC5AC",
  "SCGB1A1",
  "BPIFA1",
  "SFTPC",
  "BEST4",
  "PTPRC",
  "CD3D",
  "MS4A1",
  "LYZ",
  "HBB",
  "PECAM1",
  "COL1A1"
)
genes <- genes[genes %in% rownames(seurat_obj)]

counts <- GetAssayData(seurat_obj, slot = "counts")[genes, , drop = FALSE]
clusters <- levels(Idents(seurat_obj))

pct_mat <- sapply(
  clusters,
  function(cl) {
    cells <- WhichCells(seurat_obj, idents = cl)
    if (length(cells) == 0) {
      return(rep(NA_real_, length(genes)))
    }
    Matrix::rowMeans(counts[, cells, drop = FALSE] > 0) * 100
  },
  simplify = "matrix"
)

rownames(pct_mat) <- genes

pct_df <- data.frame(
  gene = genes,
  as.data.frame(pct_mat, check.names = FALSE),
  check.names = FALSE
)
write.csv(pct_df, "RNA_snn_res2_true_pct_counts.csv", row.names = FALSE)

# seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% c('21','26','27'), invert =TRUE)

# =========================
# DotPlot after removing 21/26/27 (RNA_snn_res.1)
# =========================

library(Seurat)
library(ggplot2)

# 2) Use RNA for interpretable marker signal (counts/log-normalized)
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- seurat_obj[["RNA_snn_res.1", drop = TRUE]]

# 3) Marker panel (ciliated/deuterosomal/cycling/secretory/basal/stress + contamination checks)
markers <- list(
  Ciliated_Mature = c(
    "FOXJ1",
    "TPPP3",
    "PIFO1",
    "DNAH5",
    "DNAI1",
    "RSPH1",
    "RSPH9",
    "CFAP46",
    "CFAP54",
    "HYDIN"
  ),
  Deuterosomal = c("CDC20B", "DEUP1", "MCIDAS", "CCNO", "MYB", "FOXN4", "PLK4"),
  Cycling = c("MKI67", "TOP2A", "HMGB2", "TYMS", "UBE2C", "CENPF"),
  Secretory_Club = c(
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "BPIFA1",
    "KRT8",
    "KRT18"
  ),
  Goblet = c("MUC5AC", "SPDEF", "AGR2", "CLCA1", "BPIFB1"),
  Basal = c("KRT5", "KRT14", "TP63", "KRT15", "NGFR"),
  Stress_AP1 = c("JUN", "FOS", "FOSB", "ATF3", "JUNB", "HSPA1A", "DUSP1"),
  Ionocyte_Tuft = c("FOXI1", "CFTR", "ASCL3", "POU2F3", "TRPM5", "IL25"),
  QC_Contam = c(
    "PTPRC",
    "CD3D",
    "MS4A1",
    "LYZ",
    "HBB",
    "PECAM1",
    "COL1A1",
    "SFTPC"
  )
)

# keep only genes present
markers <- lapply(markers, function(x) x[x %in% rownames(seurat_obj)])
features_use <- unique(unlist(markers))

# 4) DotPlot
p <- DotPlot(
  seurat_obj,
  features = features_use,
  assay = "RNA",
  cols = c("lightgrey", "red"), # adjust if you want
  dot.scale = 6
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9)
  )

# Optional: group markers with separators (visual)
# You can also split by marker categories by plotting each list separately if preferred.

print(p)

# 5) Save
pdf("DotPlot_RNA_snn_res2_after_rm_21_26_27.pdf", width = 14, height = 7)
print(p)
dev.off()

write.csv(p$data, 'temp.csv')

# seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% c('21','22'), invert =TRUE)

#!/usr/bin/env Rscript
# ============================================================
# RNA_snn_res.1 Marker QC Panel (DotPlot + Key FeaturePlot)
# Version: v1.0
# Purpose:
#   - One-shot marker panel for: histone/technical, ciliated,
#     deuterosomal, secretory/goblet/serous, basal, IEG/stress,
#     DUOX/MHC-II inflammation, immune doublets, neuronal/glial contamination
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# USER CONFIG
# ----------------------------
# seurat_obj <- readRDS("your_seurat.rds")
group_by <- "RNA_snn_res.1"
reduction <- "umap"
assay_use <- DefaultAssay(seurat_obj)

OUT_DOTPLOT <- sprintf("%s_marker_dotplot.pdf", group_by)
OUT_FEATURE <- sprintf("%s_marker_featureplots.pdf", group_by)
OUT_MARKER_TXT <- sprintf("%s_marker_gene_list.txt", group_by)

stopifnot(group_by %in% colnames(seurat_obj[[]]))
stopifnot(reduction %in% Reductions(seurat_obj))

# Make sure grouping exists and is not NULL
seurat_obj[[group_by]] <- as.factor(seurat_obj[[group_by]][, 1])

# ----------------------------
# MARKER SETS (by category)
# ----------------------------
marker_sets <- list(
  Technical_Histone = c(
    "H2AC18",
    "H2AC19",
    "H2AFZ",
    "H2AZ1",
    "H3F3A",
    "H3F3B",
    "H1-0",
    "H1-2"
  ),

  Ciliated_Mature = c(
    "FOXJ1",
    "TPPP3",
    "RSPH1",
    "PIFO",
    "DNAH5",
    "DNAI1",
    "CFAP299",
    "TTC25"
  ),

  Multiciliogenesis_Deuterosomal = c(
    "MCIDAS",
    "DEUP1",
    "CCNO",
    "CDC20B",
    "PLK4",
    "CEP152",
    "POC1A"
  ),

  Secretory_Goblet = c("MUC5AC", "SPDEF", "AGR2", "TFF3"),
  Secretory_Club = c("SCGB1A1", "SCGB3A1", "CYP2F1", "KRT19"),
  Secretory_Serous = c(
    "BPIFB1",
    "BPIFA1",
    "LCN2",
    "SLPI",
    "WFDC2",
    "PI3",
    "PIGR"
  ),

  Basal = c("KRT5", "KRT14", "TP63", "KRT15"),
  Epithelial_General = c("EPCAM", "KRT8", "KRT18", "TACSTD2"),

  IEG = c("FOS", "JUN", "JUNB", "ATF3", "EGR1"),
  Injury_Repair = c("CTGF", "CYR61", "AREG", "GDF15", "KLF4", "SOCS3"),

  Inflammation_DUOX = c("DUOX2", "DUOXA2", "IDO1", "SAA1"),
  Antigen_Presentation_MHCII = c(
    "HLA-DRA",
    "HLA-DPA1",
    "HLA-DPB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DMB"
  ),

  Immune_Doublet_Check = c(
    "PTPRC",
    "CD3D",
    "NKG7",
    "MS4A1",
    "CD79A",
    "LYZ",
    "LST1",
    "S100A8",
    "S100A9"
  ),

  Neuronal = c("RBFOX3", "TUBB3", "MAP2", "SNAP25", "SLC17A7", "GAD1", "GAD2"),
  Glial_Schwann = c("PLP1", "SOX10", "S100B", "MPZ", "GFAP"),

  Neuro_like_LongGenes_Helper = c(
    "SOX5",
    "FOXP2",
    "CNTN3",
    "TENM4",
    "GRM7",
    "NRXN3",
    "MEF2C"
  )
)

# ----------------------------
# FILTER TO PRESENT GENES
# ----------------------------
all_genes <- rownames(seurat_obj)

marker_sets_present <- lapply(marker_sets, function(v) intersect(v, all_genes))
marker_sets_missing <- lapply(marker_sets, function(v) setdiff(v, all_genes))

# Print missing summary (console)
missing_n <- sapply(marker_sets_missing, length)
cat("=== Missing gene counts by category ===\n")
print(missing_n[missing_n > 0])

# Save the final used marker list to txt (for record)
con <- file(OUT_MARKER_TXT, open = "wt")
writeLines(sprintf("group_by: %s", group_by), con)
writeLines(sprintf("assay: %s", assay_use), con)
writeLines(sprintf("reduction: %s", reduction), con)
writeLines("", con)
for (nm in names(marker_sets_present)) {
  writeLines(sprintf("[%s]", nm), con)
  writeLines(paste(marker_sets_present[[nm]], collapse = ","), con)
  writeLines("", con)
}
close(con)
cat(sprintf("Saved marker list: %s\n", OUT_MARKER_TXT))

# ----------------------------
# DOTPLOT (single figure, faceted by category)
# ----------------------------
p_dot <- DotPlot(
  object = seurat_obj,
  features = marker_sets_present, # list => will facet by category
  group.by = group_by,
  assay = assay_use
) +
  RotatedAxis() +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = sprintf("Marker QC Panel (%s)", group_by),
    x = NULL,
    y = group_by
  )

pdf(OUT_DOTPLOT, width = 18, height = 10, onefile = TRUE)
print(p_dot)
dev.off()
cat(sprintf("Saved DotPlot: %s\n", OUT_DOTPLOT))
write.csv(p_dot$data, 'temp2.csv')
# ----------------------------
# KEY FEATUREPLOTS (minimal decisive set)
# ----------------------------
key_genes <- c(
  "H2AC18", # histone/technical
  "FOXJ1",
  "TPPP3", # ciliated mature
  "MCIDAS",
  "DEUP1",
  "CCNO", # deuterosomal
  "MUC5AC",
  "SCGB1A1",
  "BPIFB1", # secretory axes
  "FOS",
  "ATF3", # IEG
  "DUOX2",
  "HLA-DRA", # inflammation/MHC-II
  "PTPRC", # immune doublet
  "TUBB3" # neuronal contamination
)
key_genes <- intersect(key_genes, all_genes)
if (length(key_genes) == 0) {
  stop("No key_genes found in this object.")
}

pdf(OUT_FEATURE, width = 10, height = 8, onefile = TRUE)
for (g in key_genes) {
  p <- FeaturePlot(
    object = seurat_obj,
    features = g,
    reduction = reduction,
    raster = FALSE,
    pt.size = 0.25
  ) +
    ggtitle(g) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 14, face = "bold"))
  print(p)
}
dev.off()
cat(sprintf("Saved FeaturePlots: %s\n", OUT_FEATURE))


suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# =========================
# 1) Gene helpers
# =========================
.present_genes <- function(obj, genes, assay = "RNA") {
  genes <- unique(genes)
  genes[genes %in% rownames(obj[[assay]])]
}

.detect_hb_genes <- function(obj, assay = "RNA") {
  # conservative Hb/RBC genes
  hb_regex <- c(
    "^HBA",
    "^HBB$",
    "^HBD$",
    "^HBE1$",
    "^HBG1$",
    "^HBG2$",
    "^HBM$",
    "^HBQ1$"
  )
  hb <- character(0)
  for (p in hb_regex) {
    hb <- c(hb, grep(p, rownames(obj[[assay]]), value = TRUE))
  }
  unique(hb)
}

# =========================
# 2) Add per-cell contamination metrics
# =========================
add_contam_metrics <- function(
  obj,
  assay = "RNA",
  airway_genes,
  cluster_col = "RNA_snn_res.1",
  verbose = TRUE
) {
  DefaultAssay(obj) <- assay
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop(sprintf("'%s' not in meta.data", cluster_col))
  }

  hb_genes <- .detect_hb_genes(obj, assay = assay)
  airway_genes <- .present_genes(obj, airway_genes, assay = assay)

  # Epithelial proxy to protect true epithelial clusters from being removed as "airway ambient"
  epi_proxy <- .present_genes(
    obj,
    c("EPCAM", "KRT19", "KRT8", "KRT18", "TACSTD2"),
    assay = assay
  )

  if (verbose) {
    message(sprintf("[HB genes] n=%d", length(hb_genes)))
    message(sprintf("[Airway genes] n=%d", length(airway_genes)))
    message(sprintf("[Epi proxy genes] n=%d", length(epi_proxy)))
  }

  obj[["pct_hb"]] <- if (length(hb_genes) > 0) {
    PercentageFeatureSet(obj, hb_genes, assay = assay)
  } else {
    0
  }
  obj[["pct_airway"]] <- if (length(airway_genes) > 0) {
    PercentageFeatureSet(obj, airway_genes, assay = assay)
  } else {
    0
  }
  obj[["pct_epi_proxy"]] <- if (length(epi_proxy) > 0) {
    PercentageFeatureSet(obj, epi_proxy, assay = assay)
  } else {
    0
  }

  obj
}

# =========================
# 3) Cluster-level removal decision
# =========================
remove_contam_clusters <- function(
  obj,
  cluster_col = "RNA_snn_res.1",
  assay = "RNA",
  # ----- thresholds (defaults: conservative) -----
  # RBC/HB cluster rule: remove if cluster median pct_hb high OR large fraction of cells high
  hb_median_cut = 5, # median pct_hb >= 5%
  hb_frac_cut = 0.50, # OR >=50% cells have pct_hb >= hb_cell_cut
  hb_cell_cut = 1, # "hb-high cell" definition: pct_hb >= 1%

  # Airway ambient cluster rule: remove only if "non-epithelial overall" AND airway signal high
  epi_median_max = 5, # non-epithelial cluster: median pct_epi_proxy < 5%
  airway_median_cut = 10, # and median pct_airway >= 10%
  airway_frac_cut = 0.60, # OR >=60% cells have pct_airway >= airway_cell_cut
  airway_cell_cut = 5, # "airway-high cell": pct_airway >= 5%

  min_cells = 30,
  verbose = TRUE
) {
  DefaultAssay(obj) <- assay
  if (
    !all(
      c("pct_hb", "pct_airway", "pct_epi_proxy") %in% colnames(obj@meta.data)
    )
  ) {
    stop(
      "Missing pct_hb / pct_airway / pct_epi_proxy. Run add_contam_metrics() first."
    )
  }
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop(sprintf("'%s' not in meta.data", cluster_col))
  }

  md <- as.data.table(obj@meta.data, keep.rownames = "cell")
  setnames(md, cluster_col, "cluster")

  # cluster stats
  stats <- md[,
    .(
      n = .N,
      med_pct_hb = median(pct_hb, na.rm = TRUE),
      frac_hb_high = mean(pct_hb >= hb_cell_cut, na.rm = TRUE),

      med_pct_airway = median(pct_airway, na.rm = TRUE),
      frac_airway_high = mean(pct_airway >= airway_cell_cut, na.rm = TRUE),

      med_pct_epi = median(pct_epi_proxy, na.rm = TRUE)
    ),
    by = cluster
  ]

  # apply minimum cells
  stats <- stats[n >= min_cells]

  # decision rules
  stats[,
    rm_hb := (med_pct_hb >= hb_median_cut) | (frac_hb_high >= hb_frac_cut)
  ]
  stats[,
    rm_airway_ambient := (med_pct_epi < epi_median_max) &
      ((med_pct_airway >= airway_median_cut) |
        (frac_airway_high >= airway_frac_cut))
  ]

  stats[, rm_cluster := rm_hb | rm_airway_ambient]
  stats[,
    rm_reason := fifelse(
      rm_hb & rm_airway_ambient,
      "HB+AirwayAmbient",
      fifelse(rm_hb, "HB", fifelse(rm_airway_ambient, "AirwayAmbient", "Keep"))
    )
  ]

  rm_clusters <- stats[rm_cluster == TRUE, as.character(cluster)]
  keep_clusters <- stats[rm_cluster == FALSE, as.character(cluster)]

  if (verbose) {
    message("[Cluster removal summary]")
    print(stats[order(-rm_cluster, -med_pct_hb, -med_pct_airway)][,
      .(
        cluster,
        n,
        med_pct_hb,
        frac_hb_high,
        med_pct_airway,
        frac_airway_high,
        med_pct_epi,
        rm_reason
      )
    ])
    message(sprintf(
      "\n[Remove clusters] n=%d: %s",
      length(rm_clusters),
      paste(rm_clusters, collapse = ",")
    ))
  }

  # subset: keep only clusters not removed
  obj$cluster_tmp_for_filter <- obj@meta.data[[cluster_col]]
  obj_clean <- subset(obj, subset = cluster_tmp_for_filter %in% keep_clusters)
  obj_clean$cluster_tmp_for_filter <- NULL

  # return everything for audit trail
  list(
    obj_clean = obj_clean,
    cluster_stats = stats,
    removed_clusters = rm_clusters,
    kept_clusters = keep_clusters
  )
}

# =========================
# 4) Example airway gene set (edit freely)
# =========================
airway_markers <- c(
  "EPCAM",
  "TACSTD2",
  "KRT5",
  "KRT14",
  "KRT19",
  "KRT8",
  "KRT18",
  "SCGB1A1",
  "SCGB3A1",
  "BPIFA1",
  "BPIFB1",
  "SCGB2B2",
  "MUC1",
  "MUC5AC",
  "MUC16",
  "AGR2",
  "WFDC2",
  "SLPI",
  "PI3",
  "LCN2",
  "DUOX2",
  "DUOXA2",
  "ALDH3A1",
  "UGT2A1",
  "CXCL17",
  "BEST4"
)

# =========================
# 5) Run (replace obj with your Seurat object)
# =========================
# obj <- readRDS("your_seurat.rds")

# Step A: compute per-cell metrics based on RNA_snn_res.1
# obj <- add_contam_metrics(obj, airway_genes = airway_markers, cluster_col = "RNA_snn_res.1")

# Step B: remove contaminated clusters (whole clusters)
# res <- remove_contam_clusters(obj, cluster_col = "RNA_snn_res.1")

# Clean object:
# obj_clean <- res$obj_clean

# Save for reproducibility:
# saveRDS(res$cluster_stats, file = "cluster_contam_stats_res2.rds")
# saveRDS(obj_clean, file = "seurat_post_clusterRemoval_res2.rds")

library(dplyr)
library(stringr)

# markers_df: data.frame with columns: cluster, gene, avg_log2FC (and others)

flag_clusters_by_technical_markers <- function(
  markers_df,
  top_n = 50,
  frac_warn = 0.30,
  frac_bad = 0.50
) {
  tech_classify <- function(g) {
    case_when(
      str_detect(g, "^MT-") ~ "mt",
      str_detect(g, "^RPL|^RPS") ~ "ribo",
      str_detect(g, "^H[1-4]") | str_detect(g, "^HIST") ~ "histone",
      str_detect(g, "^HB[AB]") ~ "hemoglobin",
      str_detect(g, "^IGH|^IGK|^IGL") ~ "immunoglobulin",
      str_detect(g, "^FOS|^JUN|^ATF3|^EGR") ~ "IEG_stress",
      str_detect(g, "^HSP") ~ "HSP_stress",
      TRUE ~ "other"
    )
  }

  out <- markers_df %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(tech = tech_classify(gene)) %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      frac_histone = mean(tech == "histone"),
      frac_mt = mean(tech == "mt"),
      frac_ribo = mean(tech == "ribo"),
      frac_hb = mean(tech == "hemoglobin"),
      frac_stress = mean(tech %in% c("IEG_stress", "HSP_stress")),
      frac_tech_total = mean(tech != "other"),
      top10_genes = paste(head(gene, 10), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      flag = case_when(
        frac_tech_total >= frac_bad ~ "BAD",
        frac_tech_total >= frac_warn ~ "WARN",
        TRUE ~ "OK"
      )
    ) %>%
    arrange(desc(frac_tech_total), desc(frac_histone), desc(frac_mt))

  return(out)
}

tech_report <- flag_clusters_by_technical_markers(sig_markers, top_n = 50)
print(tech_report)
print(subset(tech_report, flag != "OK"))
write.csv(tech_report, 'tech_report.csv')


panel <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT5", "KRT14", "TP63"),
  Immune = c("PTPRC", "LST1", "LYZ", "MS4A1", "CD3D", "NKG7"),
  APC_MHCII = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74"),
  Complement = c("C3", "C4A", "C4B", "SERPING1"),
  IFN_ISG = c(
    "ISG15",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFI44L",
    "OAS1",
    "OAS3",
    "MX1"
  ),
  Endothelial = c("PECAM1", "VWF", "KDR"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
  RBC = c("HBB", "HBA1", "HBA2", "ALAS2")
)

p0 <- DotPlot(
  seurat_obj,
  features = unlist(panel),
  group.by = "seurat_clusters"
) +
  RotatedAxis()
p1 <- VizDimLoadings(seurat_obj, dims = 1:5, reduction = "pca")

pdf('panel.pdf', width = 36, height = 24)
p0
p1
dev.off()


# 假设 Idents(seurat_obj) 是 cluster
pcs <- Embeddings(seurat_obj, "pca")[, 1:5]
seurat_obj$PC1 <- pcs[, 1]
seurat_obj$PC2 <- pcs[, 2]
seurat_obj$PC3 <- pcs[, 3]
seurat_obj$PC4 <- pcs[, 4]
seurat_obj$PC5 <- pcs[, 5]

# # QC（如果你已有就跳过）
# seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# seurat_obj$percent.ribo <- PercentageFeatureSet(seurat_obj, pattern = "^RPL|^RPS")

# 1) 看 PC2/PC4 是否和 QC 强相关（强相关=更像技术轴）
qc_mat <- seurat_obj@meta.data[, c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.ribo",
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5"
)]
round(
  cor(qc_mat, use = "pairwise.complete.obs")[
    c("PC2", "PC4"),
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
  ],
  3
)

# 2) 找“极端簇”：哪些 cluster 的 PC2 / PC4 分布显著偏移
VlnPlot(
  seurat_obj,
  features = c("PC2", "PC4"),
  group.by = "seurat_clusters",
  pt.size = 0
)

# 3) 直接在 UMAP 上看 PC2/PC4 的空间结构（是否形成“坏区域/边缘带”）
FeaturePlot(
  seurat_obj,
  features = c("PC2", "PC4"),
  reduction = "umap",
  raster = FALSE
)

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# =========================
# 0) Helpers
# =========================
.present_genes <- function(seurat_obj, genes, assay = "RNA") {
  genes <- unique(genes)
  genes[genes %in% rownames(seurat_obj[[assay]])]
}

# robust z-like outlier using MAD within each cluster
# flag if x > median + k*MAD
.mad_flag <- function(x, k = 4) {
  med <- median(x, na.rm = TRUE)
  m <- mad(x, constant = 1, na.rm = TRUE) # constant=1 => raw MAD
  # If MAD=0 (e.g., all same), never flag
  if (is.na(m) || m == 0) {
    return(rep(FALSE, length(x)))
  }
  x > (med + k * m)
}

# =========================
# 1) Add contamination scores
# =========================
add_contam_scores <- function(
  seurat_obj,
  cluster_col = "RNA_snn_res.1",
  assay = "RNA",
  airway_genes,
  hb_regex = c(
    "^HBA",
    "^HBB$",
    "^HBD$",
    "^HBE1$",
    "^HBG1$",
    "^HBG2$",
    "^HBM$",
    "^HBQ1$"
  ),
  ctrl = 100,
  verbose = TRUE
) {
  DefaultAssay(seurat_obj) <- assay

  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("cluster_col '%s' not found in meta.data", cluster_col))
  }

  # --- HB genes (auto-detect from regex)
  hb_genes <- character(0)
  for (p in hb_regex) {
    hb_genes <- c(
      hb_genes,
      grep(p, rownames(seurat_obj[[assay]]), value = TRUE)
    )
  }
  hb_genes <- unique(hb_genes)

  # --- airway genes (user-defined)
  airway_genes <- .present_genes(seurat_obj, airway_genes, assay = assay)

  if (verbose) {
    message(sprintf("[HB genes] n=%d", length(hb_genes)))
    message(sprintf("[Airway genes] n=%d", length(airway_genes)))
  }

  # Percent expression (more interpretable than module score)
  if (length(hb_genes) > 0) {
    seurat_obj[["pct_hb"]] <- PercentageFeatureSet(
      seurat_obj,
      features = hb_genes,
      assay = assay
    )
  } else {
    seurat_obj[["pct_hb"]] <- 0
  }

  if (length(airway_genes) > 0) {
    seurat_obj[["pct_airway"]] <- PercentageFeatureSet(
      seurat_obj,
      features = airway_genes,
      assay = assay
    )
  } else {
    seurat_obj[["pct_airway"]] <- 0
  }

  # Module scores (captures coordinated weak ambient signal better)
  # AddModuleScore will create columns like "HBScore1", "AirwayScore1"
  if (length(hb_genes) >= 5) {
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(hb_genes),
      assay = assay,
      ctrl = ctrl,
      name = "HBScore"
    )
  } else {
    seurat_obj[["HBScore1"]] <- 0
  }

  if (length(airway_genes) >= 10) {
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(airway_genes),
      assay = assay,
      ctrl = ctrl,
      name = "AirwayScore"
    )
  } else {
    seurat_obj[["AirwayScore1"]] <- 0
  }

  # A simple epithelial identity proxy (to avoid deleting true epithelial)
  epi_proxy <- .present_genes(
    seurat_obj,
    c("EPCAM", "KRT19", "KRT8", "KRT18", "TACSTD2"),
    assay = assay
  )
  if (length(epi_proxy) > 0) {
    seurat_obj[["pct_epi_proxy"]] <- PercentageFeatureSet(
      seurat_obj,
      features = epi_proxy,
      assay = assay
    )
  } else {
    seurat_obj[["pct_epi_proxy"]] <- 0
  }

  seurat_obj
}

# =========================
# 2) Flag contaminant cells (cluster-aware, robust)
# =========================
flag_contaminants <- function(
  seurat_obj,
  cluster_col = "RNA_snn_res.1",
  k_mad_hb = 4,
  k_mad_airway = 4,
  # only remove airway-ambient cells when epithelial proxy is low:
  epi_proxy_min = 5,
  # optional absolute floors (safety nets)
  floor_pct_hb = 1,
  floor_pct_airway = 5,
  verbose = TRUE
) {
  md <- as.data.table(seurat_obj@meta.data, keep.rownames = "cell")
  setnames(md, cluster_col, "cluster")

  # HB outliers per cluster
  md[, hb_outlier_mad := .mad_flag(pct_hb, k = k_mad_hb), by = cluster]
  md[, hb_outlier_floor := pct_hb >= floor_pct_hb]
  md[, hb_flag := hb_outlier_mad | hb_outlier_floor]

  # Airway outliers per cluster
  md[,
    airway_outlier_mad := .mad_flag(pct_airway, k = k_mad_airway),
    by = cluster
  ]
  md[, airway_outlier_floor := pct_airway >= floor_pct_airway]

  # Key logic:
  # - HB: generally safe to remove when very high (RBC carryover)
  # - Airway ambient: remove primarily when the cell is NOT truly epithelial (epi_proxy low)
  md[,
    airway_flag := (airway_outlier_mad | airway_outlier_floor) &
      (pct_epi_proxy < epi_proxy_min)
  ]

  # Combine
  md[, contam_flag := hb_flag | airway_flag]
  md[,
    contam_reason := fifelse(
      hb_flag & airway_flag,
      "HB+Airway",
      fifelse(hb_flag, "HB", fifelse(airway_flag, "Airway", "None"))
    )
  ]

  # write back
  seurat_obj$contam_flag <- md[match(colnames(seurat_obj), cell), contam_flag]
  seurat_obj$contam_reason <- md[
    match(colnames(seurat_obj), cell),
    contam_reason
  ]

  if (verbose) {
    tab <- md[, .N, by = .(contam_reason)][order(-N)]
    message("[Contam summary]")
    print(tab)

    tab2 <- md[contam_flag == TRUE, .N, by = cluster][order(-N)]
    message("[Contam counts by cluster (top 20)]")
    print(head(tab2, 20))
  }

  seurat_obj
}

# =========================
# 3) Filter and re-cluster from RNA_snn_res.1
# =========================
filter_and_recluster <- function(
  seurat_obj,
  cluster_col = "RNA_snn_res.1",
  assay = "RNA",
  keep_contam = FALSE,
  dims = 1:30,
  resolution = 2.0,
  prefix = "post_decontam",
  seed = 1
) {
  if (!"contam_flag" %in% colnames(seurat_obj@meta.data)) {
    stop("contam_flag not found. Run flag_contaminants() first.")
  }

  # filter
  if (!keep_contam) {
    seurat_obj <- subset(
      seurat_obj,
      cells = colnames(seurat_obj)[!seurat_obj$contam_flag]
    )
  } else {
    seurat_obj <- seurat_obj
  }

  DefaultAssay(seurat_obj) <- assay
  set.seed(seed)

  # Standard recluster pipeline
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    nfeatures = 4000,
    verbose = FALSE
  )
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(
    seurat_obj,
    npcs = max(dims),
    verbose = FALSE,
    reduction.name = paste0("pca_", prefix)
  )
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = paste0("pca_", prefix),
    dims = dims,
    graph.name = c(paste0(prefix, "_nn"), paste0(prefix, "_snn")),
    verbose = FALSE
  )
  seurat_obj <- FindClusters(
    seurat_obj,
    graph.name = paste0(prefix, "_snn"),
    resolution = resolution,
    verbose = FALSE
  )
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = paste0("pca_", prefix),
    dims = dims,
    reduction.name = paste0("umap_", prefix),
    seed.use = seed,
    verbose = FALSE
  )

  # cluster column will be:  paste0(prefix, "_snn_res.", resolution)
  seurat_obj
}

# =========================
# 4) Example: plug in your airway markers
# =========================
airway_markers <- c(
  # epithelial/airway typical ambient sources
  "EPCAM",
  "TACSTD2",
  "KRT5",
  "KRT14",
  "KRT19",
  "KRT8",
  "KRT18",
  "SCGB1A1",
  "SCGB3A1",
  "BPIFA1",
  "BPIFB1",
  "SCGB2B2",
  "MUC1",
  "MUC5AC",
  "MUC16",
  "AGR2",
  "WFDC2",
  "SLPI",
  "PI3",
  "LCN2",
  "DUOX2",
  "DUOXA2",
  "ALDH3A1",
  "UGT2A1",
  "CXCL17",
  # if you ALSO consider ciliated ambient (optional;按你策略决定是否加入)
  # "FOXJ1","TPPP3","CDC20B","DEUP1","MCIDAS"
  "BEST4"
)

# ---- RUN (replace 'seurat_obj' with your Seurat seurat_object) ----
# seurat_obj <- readRDS("your_seurat.rds")

# 1) scores
seurat_obj <- add_contam_scores(
  seurat_obj,
  cluster_col = "RNA_snn_res.1",
  airway_genes = airway_markers
)

# 2) flags (默认比较保守：cluster内 MAD outlier + floor，且airway只在“非上皮”里删)
seurat_obj <- flag_contaminants(
  seurat_obj,
  cluster_col = "RNA_snn_res.1",
  k_mad_hb = 4,
  k_mad_airway = 4,
  epi_proxy_min = 5,
  floor_pct_hb = 1,
  floor_pct_airway = 5
)

# 3) filter + recluster at res=2 (你要的“之后的 cluster 源自 RNA_snn_res.1”，这里即在去污染后重跑并产生新的res.2)
seurat_obj_clean <- filter_and_recluster(
  seurat_obj,
  cluster_col = "RNA_snn_res.1",
  resolution = 2.0,
  prefix = ""
)

seurat_obj <- subset(seurat_obj, subset = pct_hb <= 0.01551)
# New cluster col: "clean_r2_snn_res.2"
# Idents(seurat_obj_clean) <- "clean_r2_snn_res.2"
