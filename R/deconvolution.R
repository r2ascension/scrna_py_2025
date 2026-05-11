library(MuSiC)
library(Seurat)
library(SingleCellExperiment)
library(Biobase)
library(org.Hs.eg.db)
library(MuSiC)
library(Biobase)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
source('/home/h2048/script/R/RenameENSG_v2.R')
output_dir <- "/home/h2048/data/R/bulk0121/"
out_dir <- output_dir
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
library(reticulate)
library(data.table)
library(harmony)
library(ggplot2)
library(patchwork)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()

# 1. 读取单细胞参考数据 (RDS文件)
# sc_ref 通常是 SingleCellExperiment 或 ExpressionSet 对象
sc_ref <- GetSeurat(
  "/home/h2048/data/py/The integrated Human Lung Cell Atlas.h5ad"
)
sc_ref <- convert_and_cleanup(
  sc_ref,
  source_layer = "counts",
  new_assay_name = "RNA.symbol",
  handle_duplicates = "sum",
  remove_old_assay = TRUE,
  set_as_default = TRUE
)
sc_ref <- fix_counts_assay_name(
  sc_ref,
  from_assay = "RNA.symbol",
  to_assay = "RNA",
  set_default = TRUE
)
head(sc_ref)
sce_ref <- as.SingleCellExperiment(sc_ref)

# ==============================================================================
# 1. 读取数据（关键：暂时不设置 row.names = 1）
# ==============================================================================
# stringsAsFactors = FALSE 确保基因名是字符型
bulk_raw <- read.csv(
  "/home/h2048/data/bulk/gene_count.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 检查一下第一列的列名，假设是 "gene_name" (根据您的报错信息推测)
# 如果第一列名字不是 gene_name，请在下面代码中替换
colnames(bulk_raw)[1] <- "gene_name"


# ==============================================================================
# 2. 处理重复基因 (去重/聚合)
# ==============================================================================
# 对于做反卷积，如果有重复的 Gene Symbol，最合理的做法是将它们的 Counts 相加
# 这样能保留该基因的总表达量
bulk_agg <- bulk_raw %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum)) %>% # 对所有数值列求和
  as.data.frame()

# 检查是否还有重复（理论上不会有了）
# any(duplicated(bulk_agg$gene_name))

# ==============================================================================
# 3. 转换为 MuSiC 需要的矩阵格式
# ==============================================================================
# 将基因名设为行名
rownames(bulk_agg) <- bulk_agg$gene_name

# 移除原来的基因名列，只保留数值矩阵
bulk_mtx <- as.matrix(bulk_agg[, -1])

# 确保矩阵是数值型（防止意外混入字符）
class(bulk_mtx) # 应该是 "matrix" "array"


# ==============================================================================
# 4. 构建 ExpressionSet 并运行 MuSiC
# ==============================================================================
bulk_eset <- ExpressionSet(assayData = bulk_mtx)


# ==============================================================================
# 4. 执行反卷积
# ==============================================================================
# 请确保 clusters (细胞类型) 和 samples (样本ID) 对应您 sc_ref 的实际列名
music_res <- music_prop(
  bulk.mtx = exprs(bulk_eset), # 或者直接传 bulk_mtx
  sc.sce = sce_ref, # 单细胞参考对象
  clusters = 'ann_level_3', # 替换为您单细胞数据中细胞类型的列名
  samples = 'sample', # 替换为您单细胞数据中供体/样本的列名
  verbose = TRUE
)

# 5. 查看结果
head(music_res$Est.prop.weighted)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 1) Extract results
# -------------------------------
props_w <- as.data.frame(music_res$Est.prop.weighted) # samples x cell_types
props_n <- as.data.frame(music_res$Est.prop.allgene) # NNLS baseline (samples x cell_types)
r2_vec <- music_res$r.squared.full # named vector (per sample) in MuSiC output
var_mat <- music_res$Var.prop # variance matrix (cell_types x samples or similar)

# 如果行名丢了，尝试从 bulk.mtx 的 colnames 或 r2 names 补
if (is.null(rownames(props_w)) || any(rownames(props_w) == "")) {
  if (!is.null(names(r2_vec)) && length(names(r2_vec)) == nrow(props_w)) {
    rownames(props_w) <- names(r2_vec)
    rownames(props_n) <- names(r2_vec)
  } else {
    rownames(props_w) <- paste0("S", seq_len(nrow(props_w)))
    rownames(props_n) <- rownames(props_w)
  }
}

# -------------------------------
# 2) Sanity checks
# -------------------------------
check_df <- data.frame(
  sample = rownames(props_w),
  sum_weighted = rowSums(props_w),
  min_weighted = apply(props_w, 1, min),
  max_weighted = apply(props_w, 1, max)
)

write.csv(
  check_df,
  file.path(out_dir, "QC_prop_weighted_checks.csv"),
  row.names = FALSE
)

# quick console checks
print(check_df)

# 负值检查（一般不应出现；若有，通常是数值问题/异常输入）
if (any(props_w < -1e-8, na.rm = TRUE)) {
  warning(
    "Detected negative proportions in Est.prop.weighted. Please inspect inputs/output."
  )
}

# -------------------------------
# 3) Plot: R^2 per sample
# -------------------------------
qc_r2 <- data.frame(
  sample = names(r2_vec) %||% rownames(props_w),
  r2 = as.numeric(r2_vec)
)

p_r2 <- ggplot(qc_r2, aes(x = sample, y = r2)) +
  geom_col() +
  theme_bw() +
  labs(
    title = "MuSiC goodness-of-fit (R^2) per sample",
    x = "Sample",
    y = "R^2"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "QC_R2_per_sample.pdf"), p_r2, width = 7, height = 4)

# -------------------------------
# 4) Plot: Stacked bar (weighted)
# -------------------------------
df_w_long <- props_w %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "prop") %>%
  group_by(sample) %>%
  mutate(prop = pmax(prop, 0)) %>% # just in case tiny negative numerical noise
  ungroup()

p_stack_w <- ggplot(df_w_long, aes(x = sample, y = prop, fill = cell_type)) +
  geom_col(width = 0.85) +
  theme_bw() +
  labs(
    title = "Estimated cell-type proportions (MuSiC weighted)",
    x = "Sample",
    y = "Proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "Prop_stackedbar_weighted.pdf"),
  p_stack_w,
  width = 9,
  height = 5
)

# -------------------------------
# 5) Plot: Stacked bar (NNLS allgene)
# -------------------------------
df_n_long <- props_n %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "prop") %>%
  group_by(sample) %>%
  mutate(prop = pmax(prop, 0)) %>%
  ungroup()

p_stack_n <- ggplot(df_n_long, aes(x = sample, y = prop, fill = cell_type)) +
  geom_col(width = 0.85) +
  theme_bw() +
  labs(
    title = "Estimated cell-type proportions (NNLS allgene baseline)",
    x = "Sample",
    y = "Proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "Prop_stackedbar_allgene_NNLS.pdf"),
  p_stack_n,
  width = 9,
  height = 5
)

# -------------------------------
# 6) Plot: Heatmap (weighted proportions)
# -------------------------------
mat_w <- as.matrix(props_w)
# 可选：按 cell type 维度做 z-score，突出相对变化（仅 5 样本时更容易看“谁相对高/低”）
# scale = "row" 会对每个 cell type 在不同样本间中心化/标准化
pdf(
  file.path(out_dir, "Prop_heatmap_weighted_scale_row.pdf"),
  width = 7,
  height = 6
)
pheatmap(
  mat_w,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "MuSiC weighted proportions (row-scaled)"
)
dev.off()

# 同时给一张不 scale 的绝对比例热图
pdf(
  file.path(out_dir, "Prop_heatmap_weighted_absolute.pdf"),
  width = 7,
  height = 6
)
pheatmap(
  mat_w,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "MuSiC weighted proportions (absolute)"
)
dev.off()

# -------------------------------
# 7) Plot: PCA on proportions (CLR transform recommended)
#    CLR: log(p) - mean(log(p)) per sample; add small pseudocount to avoid log(0)
# -------------------------------
clr_transform <- function(P, eps = 1e-6) {
  P2 <- as.matrix(P)
  P2[P2 < 0] <- 0
  P2 <- P2 + eps
  P2 <- P2 / rowSums(P2)
  L <- log(P2)
  L - rowMeans(L)
}

mat_clr <- clr_transform(props_w)

pca <- prcomp(mat_clr, center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  sample = rownames(mat_clr),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7) +
  theme_bw() +
  labs(
    title = "PCA on CLR-transformed proportions (MuSiC weighted)",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  file.path(out_dir, "Prop_PCA_CLR_weighted.pdf"),
  p_pca,
  width = 6,
  height = 5
)

# -------------------------------
# 8) Plot: Variability across healthy samples (CV per cell type)
# -------------------------------
cv_df <- data.frame(
  cell_type = colnames(props_w),
  mean_prop = colMeans(props_w),
  sd_prop = apply(props_w, 2, sd),
  cv = apply(props_w, 2, sd) / pmax(colMeans(props_w), 1e-8)
) %>%
  arrange(desc(cv))

write.csv(
  cv_df,
  file.path(out_dir, "Prop_CV_by_celltype.csv"),
  row.names = FALSE
)

p_cv <- ggplot(cv_df, aes(x = reorder(cell_type, cv), y = cv)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Across-sample variability (CV) by cell type (MuSiC weighted)",
    x = "Cell type",
    y = "Coefficient of variation (CV)"
  )

ggsave(
  file.path(out_dir, "Prop_CV_by_celltype_weighted.pdf"),
  p_cv,
  width = 7,
  height = 6
)

# -------------------------------
# 9) Optional: Compare weighted vs NNLS per sample (sanity)
# -------------------------------
# 计算两种方法在每个样本上的相关（在 cell type 维度）
cmp_df <- data.frame(
  sample = rownames(props_w),
  cor_w_vs_nnls = sapply(rownames(props_w), function(s) {
    suppressWarnings(cor(
      as.numeric(props_w[s, ]),
      as.numeric(props_n[s, ]),
      method = "pearson"
    ))
  })
)
write.csv(
  cmp_df,
  file.path(out_dir, "Compare_weighted_vs_NNLS_corr.csv"),
  row.names = FALSE
)

p_cmp <- ggplot(cmp_df, aes(x = sample, y = cor_w_vs_nnls)) +
  geom_point(size = 3) +
  theme_bw() +
  ylim(0, 1) +
  labs(
    title = "Agreement: MuSiC weighted vs NNLS (Pearson r across cell types)",
    x = "Sample",
    y = "Pearson r"
  )

ggsave(
  file.path(out_dir, "Compare_weighted_vs_NNLS_corr.pdf"),
  p_cmp,
  width = 7,
  height = 4
)

message("Done. Plots saved to: ", out_dir)

# helper for base R "%||%" (like rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Paths
# -----------------------------
run_dir <- "/home/h2048/data/R/bulk0121/deconv_run_20260121"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(run_dir, "input"), showWarnings = FALSE)
dir.create(file.path(run_dir, "reference"), showWarnings = FALSE)
dir.create(file.path(run_dir, "mapping"), showWarnings = FALSE)
dir.create(file.path(run_dir, "results"), showWarnings = FALSE)
dir.create(file.path(run_dir, "logs"), showWarnings = FALSE)

# -----------------------------
# Assume you already have:
# bulk_agg (gene x sample counts matrix)
# sce_ref  (SingleCellExperiment aligned to common genes)
# common_genes (character)
# music_res (MuSiC output list)
# -----------------------------

# 1) Bulk input (recommended: RDS for exact type + fast load)
saveRDS(bulk_agg, file.path(run_dir, "input", "bulk_counts_agg.rds"))
write.csv(
  data.frame(sample = colnames(bulk_agg)),
  file.path(run_dir, "input", "bulk_meta.csv"),
  row.names = FALSE
)

# 2) Single-cell reference
# If small in-memory -> saveRDS is fine
# If large / HDF5-backed -> saveHDF5SummarizedExperiment is safer (relocatable)
if (requireNamespace("HDF5Array", quietly = TRUE)) {
  # SingleCellExperiment is a SummarizedExperiment derivative -> compatible
  HDF5Array::saveHDF5SummarizedExperiment(
    sce_ref,
    dir = file.path(run_dir, "reference", "sce_ref_h5se"),
    replace = TRUE
  )
} else {
  saveRDS(sce_ref, file.path(run_dir, "reference", "sce_ref.rds"))
}
write.csv(
  as.data.frame(SummarizedExperiment::colData(sce_ref)),
  file.path(run_dir, "reference", "sce_ref_coldata.csv"),
  row.names = FALSE
)

# 3) Mapping artifacts
writeLines(common_genes, file.path(run_dir, "mapping", "common_genes.txt"))
# 如果您有 gene_map（如 SYMBOL 去重聚合、ENSG->SYMBOL），也建议写出：
# write.table(gene_map, file.path(run_dir,"mapping","gene_map.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# 4) MuSiC outputs: save full + export key tables
saveRDS(music_res, file.path(run_dir, "results", "music_res.rds"))

write.csv(
  music_res$Est.prop.weighted,
  file.path(run_dir, "results", "Est.prop.weighted.csv")
)
write.csv(
  music_res$Est.prop.allgene,
  file.path(run_dir, "results", "Est.prop.allgene.csv")
)
write.csv(
  data.frame(
    sample = names(music_res$r.squared.full),
    r2 = as.numeric(music_res$r.squared.full)
  ),
  file.path(run_dir, "results", "r.squared.full.csv"),
  row.names = FALSE
)
write.csv(music_res$Var.prop, file.path(run_dir, "results", "Var.prop.csv"))

# 5) Manifest (参数固化)
manifest <- list(
  clusters_col = "ann_level_3",
  samples_col = "sample",
  bulk_samples = colnames(bulk_agg),
  n_common_genes = length(common_genes),
  saved_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(run_dir, "logs", "manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# 6) Session info (版本固化)
sink(file.path(run_dir, "logs", "sessionInfo.txt"))
print(sessionInfo())
sink()

# 7) Checksums (防止输入文件悄悄变了)
md5 <- data.frame(
  file = c("input/bulk_counts_agg.rds", "results/music_res.rds"),
  md5 = tools::md5sum(file.path(
    run_dir,
    c("input/bulk_counts_agg.rds", "results/music_res.rds")
  ))
)
write.csv(md5, file.path(run_dir, "logs", "md5sum.tsv"), row.names = FALSE)

wg <- music_res$Weight.gene
wg <- sort(wg, decreasing = TRUE)
head(wg, 50)
music_res$r.squared.full
head(sort(music_res$Weight.gene, decreasing = TRUE), 30)
getwd()
saveRDS(sc_ref, 'LungMap.rds')
