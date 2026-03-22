# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
# setwd("/home/h2048/data/R/1218")
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

# ===== Remove Ribosomal, Mitochondrial, ENSG, and Unannotated Genes =====

# Configuration
remove_mt <- TRUE
remove_ribo <- TRUE
remove_pseudogenes <- TRUE
remove_ensg <- TRUE
remove_unannotated <- TRUE

# Get all gene names
all_genes <- rownames(seurat_obj)

# Initialize genes to remove
genes_to_remove <- c()

# 1. Mitochondrial genes (MT-)
if (remove_mt) {
  mt_genes <- grep('^MT-', all_genes, value = TRUE)
  genes_to_remove <- c(genes_to_remove, mt_genes)
  cat(sprintf('Mitochondrial genes: %d\n', length(mt_genes)))
}

# 2. Ribosomal genes (RPS, RPL, MRPS, MRPL)
if (remove_ribo) {
  ribo_genes <- grep('^RPS|^RPL|^MRPS|^MRPL', all_genes, value = TRUE)
  genes_to_remove <- c(genes_to_remove, ribo_genes)
  cat(sprintf('Ribosomal genes: %d\n', length(ribo_genes)))
}

# 3. Ribosomal pseudogenes (RPS29P1, RPL10P9, MRPS36P1, etc.)
if (remove_pseudogenes) {
  # Match: RPS/RPL/MRPS/MRPL + digits + P + digits
  pseudo_genes <- grep(
    '^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$',
    all_genes,
    value = TRUE
  )
  genes_to_remove <- c(genes_to_remove, pseudo_genes)
  cat(sprintf('Pseudogenes: %d\n', length(pseudo_genes)))
}

# 4. ENSG unannotated genes
if (remove_ensg) {
  ensg_genes <- grep('^ENSG[0-9]+', all_genes, value = TRUE)
  genes_to_remove <- c(genes_to_remove, ensg_genes)
  cat(sprintf('ENSG genes: %d\n', length(ensg_genes)))
}

unannotated_genes <- grep(
  '^(AC|AL|AP|BX|Z)[0-9]+\\.|^RP[0-9]+-|^CTD-|^CTB-|^CTC-|^LINC[0-9]+|-AS[0-9]+$|-OT[0-9]+$|^LOC[0-9]+',
  all_genes,
  value = TRUE
)

# # 5. Unannotated transcripts
# if (remove_unannotated) {
#   unannotated_genes <- grep(
#     '^(AC|AL|AP|BX|Z)[0-9]+\\.|^RP[0-9]+-|^CT[DBCS]-|^LINC[0-9]+|-AS[0-9]+$|-OT[0-9]+$',
#     all_genes,
#     value = TRUE
#   )
#   genes_to_remove <- c(genes_to_remove, unannotated_genes)
#   cat(sprintf('Unannotated transcripts: %d\n', length(unannotated_genes)))
# }

genes_to_remove <- unique(c(unannotated_genes, genes_to_remove))

# Remove duplicates
genes_to_remove <- unique(genes_to_remove)
cat(sprintf('\nTotal genes to remove: %d\n', length(genes_to_remove)))

# Keep genes
genes_to_keep <- setdiff(all_genes, genes_to_remove)
cat(sprintf('Genes to keep: %d\n', length(genes_to_keep)))

# Subset Seurat object
seurat_obj <- subset(seurat_obj, features = genes_to_keep)


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


seurat_obj <- readRDS(
  '/home/h2048/data/R/1221/ciliated/ciliated_20251222.rds'
)
output_dir <- '/home/h2048/data/R/1221/ciliated'
dir.create(output_dir)
setwd(output_dir)
seurat_obj <- NormalizeData(seurat_obj) #归一化
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000
) #寻找变异基因
seurat_obj <- ScaleData(seurat_obj) #标准化

seurat_obj <- RunPCA(seurat_obj, npcs = 30)

anno <- fread("/home/h2048/data/R/1217/Annotation.csv") # 你的CSV：含 Cluster, Size, Percent, Suggestion, Dominant_GEP, Manual_Annotation
anno[, Cluster := as.character(Cluster)]

cluster_col <- "RNA_snn_res.1" # 如果你的cluster列不是这个名，改成对应列名
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

cnt <- get_counts_matrix(seurat_obj)

gene_ncells <- Matrix::rowSums(cnt > 0)
keep_genes <- names(gene_ncells[gene_ncells >= 3])

before_genes <- nrow(seurat_obj)
seurat_obj <- subset(seurat_obj, features = keep_genes)
after_genes <- nrow(seurat_obj)
before_genes
after_genes

seurat_obj <- NormalizeData(seurat_obj) #归一化
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
  theta = c(1), # Higher theta for more diverse clustering
  lambda = c(7), # Higher lambda to reduce overcorrection
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
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  k.param = 45
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
write.csv(top_markers, "all_markers_top20.csv", row.names = FALSE)
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
  group.by = "RNA_snn_res.1",
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

cluster_col <- "RNA_snn_res.1"

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

saveRDS(seurat_obj, 'ciliated_obj_doublet_removed_20251222_2.rds')

# ===== Remove doublets AND pure B cells =====
seurat_obj <- readRDS('/home/h2048/data/R/1221/ciliated/ciliated_20251222.rds')
# Step 1: Calculate safe scores (if not already done)
b_cell_safe <- c('CD79A', 'CD79B', 'MS4A1', 'MZB1', 'JCHAIN', 'IGKC')
epithelial_safe <- c('KRT5', 'TP63', 'EPCAM', 'CDH1')

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(b_cell_safe),
  name = 'B_safe'
)

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(epithelial_safe),
  name = 'Epi_safe'
)

# Step 2: Classify cells into removal categories
seurat_obj$removal_reason <- 'Keep'

# B-Epithelial doublets
seurat_obj$removal_reason[
  seurat_obj$B_safe1 > 0.3 & seurat_obj$Epi_safe1 > 0.3
] <-
  'B-Epithelial_Doublet'

# Pure B cells (contamination)
seurat_obj$removal_reason[
  seurat_obj$B_safe1 > 0.7 & seurat_obj$Epi_safe1 < 0.3
] <-
  'B_cell_Contamination'

# Flag all cells to remove
seurat_obj$to_remove <- seurat_obj$removal_reason != 'Keep'

# Step 3: Extract all cells to be removed
removed_cells <- seurat_obj@meta.data[seurat_obj$to_remove, ]
removed_cells$cell_barcode <- rownames(removed_cells)

# Step 4: Prepare export data
# Check if 'sample' column exists
if ('sample' %in% colnames(removed_cells)) {
  sample_col <- 'sample'
} else if ('Sample' %in% colnames(removed_cells)) {
  sample_col <- 'Sample'
} else {
  sample_col <- 'orig.ident'
  cat('Warning: No "sample" column found, using orig.ident instead\n')
}

removed_export <- removed_cells[, c(
  'cell_barcode',
  sample_col,
  'nCount_RNA',
  'nFeature_RNA',
  'percent.mt',
  'RNA_snn_res.1',
  'B_safe1',
  'Epi_safe1',
  'removal_reason'
)]

# Rename columns
colnames(removed_export) <- c(
  'Cell_Barcode',
  'Sample',
  'Total_UMI',
  'Total_Genes',
  'Percent_MT',
  'Cluster',
  'B_cell_Score',
  'Epithelial_Score',
  'Removal_Reason'
)

# Sort by removal reason for easier viewing
removed_export <- removed_export[order(removed_export$Removal_Reason), ]

# Step 5: Save to CSV
write.csv(
  removed_export,
  'removed_cells_all.csv',
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf('Total cells to remove: %d\n', nrow(removed_export)))
cat(sprintf(
  '  - B-Epithelial doublets: %d\n',
  sum(removed_export$Removal_Reason == 'B-Epithelial_Doublet')
))
cat(sprintf(
  '  - B cell contamination: %d\n',
  sum(removed_export$Removal_Reason == 'B_cell_Contamination')
))

# Step 6: Summary by removal reason and cluster
removal_summary <- table(
  removed_cells$RNA_snn_res.1,
  removed_cells$removal_reason
)

write.csv(
  as.data.frame.matrix(removal_summary),
  'removal_summary_by_cluster.csv',
  row.names = TRUE
)

print(removal_summary)

# Step 7: Summary by sample
removal_by_sample <- table(
  removed_cells[[sample_col]],
  removed_cells$removal_reason
)

write.csv(
  as.data.frame.matrix(removal_by_sample),
  'removal_summary_by_sample.csv',
  row.names = TRUE
)

print(removal_by_sample)

# Step 8: Visualize
library(ggplot2)

# Score scatter plot with removal categories
plot_data <- FetchData(
  seurat_obj,
  vars = c('B_safe1', 'Epi_safe1', 'removal_reason')
)

p1 <- ggplot(
  plot_data,
  aes(x = Epi_safe1, y = B_safe1, color = removal_reason)
) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = 0.3, linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = 0.3, linetype = 'dashed', color = 'black') +
  scale_color_manual(
    values = c(
      'Keep' = 'grey80',
      'B-Epithelial_Doublet' = 'red',
      'B_cell_Contamination' = 'blue'
    )
  ) +
  theme_bw() +
  labs(
    title = 'Cell Removal Strategy',
    x = 'Epithelial Score',
    y = 'B Cell Score',
    color = 'Cell Type'
  ) +
  annotate(
    'text',
    x = 0.15,
    y = 0.6,
    label = 'Pure B cells\n(Remove)',
    color = 'blue',
    size = 3
  ) +
  annotate(
    'text',
    x = 0.6,
    y = 0.6,
    label = 'Doublets\n(Remove)',
    color = 'red',
    size = 3
  ) +
  annotate(
    'text',
    x = 0.6,
    y = 0.15,
    label = 'Epithelial\n(Keep)',
    color = 'grey40',
    size = 3
  )

# UMAP visualization
p2 <- DimPlot(
  seurat_obj,
  group.by = 'removal_reason',
  cols = c(
    'Keep' = 'grey80',
    'B-Epithelial_Doublet' = 'red',
    'B_cell_Contamination' = 'blue'
  )
) +
  ggtitle(sprintf('Cells to Remove (n=%d)', sum(seurat_obj$to_remove)))

pdf('cell_removal_visualization.pdf', width = 14, height = 6)
print(p1 | p2)
dev.off()

# Step 9: Remove flagged cells
cat(sprintf('\nBefore filtering: %d cells\n', ncol(seurat_obj)))

seurat_obj <- subset(seurat_obj, subset = to_remove == FALSE)

cat(sprintf('After filtering: %d cells\n', ncol(seurat_obj)))
cat(sprintf(
  'Removed: %d cells (%.2f%%)\n',
  ncol(seurat_obj) - ncol(seurat_obj),
  100 * (ncol(seurat_obj) - ncol(seurat_obj)) / ncol(seurat_obj)
))

# Step 10: Save cleaned object
saveRDS(seurat_obj, 'ciliated_20251222_2.rds')

cat('\n=== Saved files ===\n')
cat('1. removed_cells_all.csv - All removed cells with reasons\n')
cat('2. removal_summary_by_cluster.csv - Summary by cluster\n')
cat('3. removal_summary_by_sample.csv - Summary by sample\n')
cat('4. cell_removal_visualization.pdf - Visualization\n')
cat('5. seurat_obj_cleaned_final.rds - Cleaned Seurat object\n')

# saveRDS(seurat_obj, 'ciliated_obj_doublet_removed_20251222.rds')
seurat_obj$removal_reason <- 'Keep'
# Pure B cells (contamination)
seurat_obj$removal_reason[seurat_obj$RNA_snn_res.1 == '18'] <-
  'Myeloid_cell_Contamination'

# Flag all cells to remove
seurat_obj$to_remove <- seurat_obj$removal_reason != 'Keep'

# Step 3: Extract all cells to be removed
removed_cells <- seurat_obj@meta.data[seurat_obj$to_remove, ]
removed_cells$cell_barcode <- rownames(removed_cells)

# Step 4: Prepare export data
# Check if 'sample' column exists
if ('sample' %in% colnames(removed_cells)) {
  sample_col <- 'sample'
} else if ('Sample' %in% colnames(removed_cells)) {
  sample_col <- 'Sample'
} else {
  sample_col <- 'orig.ident'
  cat('Warning: No "sample" column found, using orig.ident instead\n')
}

removed_export <- removed_cells[, c(
  'cell_barcode',
  sample_col,
  'nCount_RNA',
  'nFeature_RNA',
  'percent.mt',
  'RNA_snn_res.1',
  'B_safe1',
  'Epi_safe1',
  'removal_reason'
)]

# Rename columns
colnames(removed_export) <- c(
  'Cell_Barcode',
  'Sample',
  'Total_UMI',
  'Total_Genes',
  'Percent_MT',
  'Cluster',
  'B_cell_Score',
  'Epithelial_Score',
  'Removal_Reason'
)

# Sort by removal reason for easier viewing
removed_export <- removed_export[order(removed_export$Removal_Reason), ]

# Step 5: Save to CSV
write.csv(
  removed_export,
  'removed_cells_all.csv',
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf('Total cells to remove: %d\n', nrow(removed_export)))
cat(sprintf(
  '  - B-Epithelial doublets: %d\n',
  sum(removed_export$Removal_Reason == 'B-Epithelial_Doublet')
))
cat(sprintf(
  '  - B cell contamination: %d\n',
  sum(removed_export$Removal_Reason == 'B_cell_Contamination')
))

# Step 6: Summary by removal reason and cluster
removal_summary <- table(
  removed_cells$RNA_snn_res.1,
  removed_cells$removal_reason
)

write.csv(
  as.data.frame.matrix(removal_summary),
  'removal_summary_by_cluster.csv',
  row.names = TRUE
)

print(removal_summary)

# Step 7: Summary by sample
removal_by_sample <- table(
  removed_cells[[sample_col]],
  removed_cells$removal_reason
)

write.csv(
  as.data.frame.matrix(removal_by_sample),
  'removal_summary_by_sample.csv',
  row.names = TRUE
)

print(removal_by_sample)

# Step 9: Remove flagged cells
cat(sprintf('\nBefore filtering: %d cells\n', ncol(seurat_obj)))

seurat_obj <- subset(seurat_obj, subset = to_remove == FALSE)

cat(sprintf('After filtering: %d cells\n', ncol(seurat_obj)))

saveRDS(seurat_obj, 'ciliated_obj_doublet_removed_20251222_3.rds')

# ===== Merge clusters and recalculate markers =====

# ===== Fix: Convert factor to character first =====

# Convert to character to avoid factor level issues
new_clusters <- as.character(seurat_obj$RNA_snn_res.1)

# Merge clusters
new_clusters[new_clusters %in% c('3', '11', '16')] <- 'c1'
new_clusters[new_clusters %in% c('17')] <- 'c2'
new_clusters[new_clusters %in% c('10')] <- 'c3'
new_clusters[new_clusters %in% c('4')] <- 'c4'
new_clusters[new_clusters %in% c('1', '7', '8')] <- 'c5'
new_clusters[new_clusters %in% c('5', '6', '12')] <- 'c6'
new_clusters[new_clusters %in% c('2')] <- 'c7'
new_clusters[new_clusters %in% c('15')] <- 'c8'
new_clusters[new_clusters %in% c('13')] <- 'c9'
new_clusters[new_clusters %in% c('14', '9')] <- 'c10'

seurat_obj$merged_clusters <- factor(new_clusters)
Idents(seurat_obj) <- 'merged_clusters'
table(seurat_obj$merged_clusters)
# Calculate markers
merged_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = 'wilcox'
)

write.csv(merged_markers, 'merged_clusters_markers.csv', row.names = FALSE)

# Top 20
top20 <- merged_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(top20, 'merged_clusters_top20_markers.csv', row.names = FALSE)
# ===== 1.1.1 定义关键marker基因 =====

# Club细胞成熟marker
club_mature_genes <- c("SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFA1")

# 祖细胞marker
progenitor_genes <- c("ALDH1A1", "MSI2", "SOX2", "SOX9", "BMI1")

# 上皮分化转录因子
epithelial_tf_genes <- c("GRHL2", "EHF", "TFCP2L1", "ELF5", "OVOL2")

# 线粒体/代谢marker
mitochondrial_genes <- c("MTRNR2L8", "MTRNR2L12", "MT-CO1", "MT-ND1")

# 杯状细胞分化marker
goblet_genes <- c("MUC5AC", "MUC5B", "SPDEF", "AGR2", "TFF3")

# 纤毛细胞marker
ciliated_genes <- c("FOXJ1", "RSPH1", "DNAH5", "CCDC39")

# 基底细胞marker (应该阴性)
basal_genes <- c("TP63", "KRT5", "KRT14", "KRT15")

# ===== 1.1.2 UMAP可视化关键marker =====

# Cluster 3特异性高亮
p_cluster3_highlight <- DimPlot(
  seurat_obj,
  cells.highlight = WhichCells(seurat_obj, idents = 3),
  cols.highlight = "red",
  cols = "grey90",
  pt.size = 0.5
) +
  ggtitle("Cluster 3 Location") +
  theme_minimal()

# Club细胞marker
p_club <- FeaturePlot(
  seurat_obj,
  features = club_mature_genes,
  ncol = 2,
  pt.size = 0.5,
  order = TRUE,
  cols = c("lightgrey", "red")
)

# 祖细胞marker
p_progenitor <- FeaturePlot(
  seurat_obj,
  features = progenitor_genes,
  ncol = 3,
  pt.size = 0.5,
  order = TRUE,
  cols = c("lightgrey", "blue")
)

# 上皮TF
p_epi_tf <- FeaturePlot(
  seurat_obj,
  features = epithelial_tf_genes,
  ncol = 3,
  pt.size = 0.5,
  order = TRUE,
  cols = c("lightgrey", "darkgreen")
)

# 保存图片
ggsave(
  filename = file.path(output_dir, "Cluster3_01_location.pdf"),
  plot = p_cluster3_highlight,
  width = 6,
  height = 5
)
ggsave(
  filename = file.path(output_dir, "Cluster3_02_club_markers.pdf"),
  plot = p_club,
  width = 8,
  height = 8
)
ggsave(
  filename = file.path(output_dir, "Cluster3_03_progenitor_markers.pdf"),
  plot = p_progenitor,
  width = 12,
  height = 8
)
ggsave(
  filename = file.path(output_dir, "Cluster3_04_epithelial_TF.pdf"),
  plot = p_epi_tf,
  width = 12,
  height = 8
)

# ===== 1.1.3 排除性marker检查 =====

# 检查是否混入其他细胞类型
p_exclusion <- FeaturePlot(
  seurat_obj,
  features = c("TP63", "FOXJ1", "MUC5AC", "CD3D"),
  ncol = 2,
  pt.size = 0.5
)


ggsave(
  filename = file.path(output_dir, "Cluster3_05_exclusion_markers.pdf"),
  plot = p_exclusion,
  width = 8,
  height = 8
)

# ===== 1.1.4 VlnPlot定量对比 =====

# 准备数据
all_markers <- c(
  club_mature_genes,
  progenitor_genes[1:3],
  epithelial_tf_genes[1:3],
  basal_genes[1:2]
)

p_violin <- VlnPlot(
  seurat_obj,
  features = all_markers,
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 4
)


ggsave(
  filename = file.path(output_dir, "Cluster3_06_violin_quantitative.pdf"),
  plot = p_violin,
  width = 16,
  height = 12
)


# ===== 1.2.1 计算多个功能评分 =====

# Club细胞成熟度评分
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(club_mature_genes),
  name = "Club_Maturity"
)

# 祖细胞潜能评分
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(progenitor_genes),
  name = "Progenitor_Score"
)

# 杯状细胞倾向评分
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(goblet_genes),
  name = "Goblet_Tendency"
)

# 纤毛细胞倾向评分
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(ciliated_genes),
  name = "Ciliated_Tendency"
)

# 线粒体/氧化磷酸化评分
oxidative_phos_genes <- c(
  "ATP5A1",
  "ATP5C1",
  "NDUFA4",
  "NDUFB2",
  "UQCRQ",
  "COX5A",
  "COX5B",
  "COX6C",
  "MT-CO1",
  "MT-CO2",
  "MT-ND1",
  "MT-ND2"
)

seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(oxidative_phos_genes),
  name = "OxPhos_Score"
)

# ===== 1.2.2 评分可视化 =====

# UMAP展示各个评分
score_features <- c(
  "Club_Maturity1",
  "Progenitor_Score1",
  "Goblet_Tendency1",
  "Ciliated_Tendency1",
  "OxPhos_Score1"
)

p_scores_umap <- FeaturePlot(
  seurat_obj,
  features = score_features,
  ncol = 3,
  pt.size = 0.5,
  cols = c("lightgrey", "red")
) +
  plot_annotation(title = "Cluster 3 Functional Scores")

ggsave(
  filename = file.path(output_dir, "Cluster3_07_functional_scores_umap.pdf"),
  plot = p_scores_umap,
  width = 15,
  height = 10
)

# VlnPlot对比不同cluster
p_scores_violin <- VlnPlot(
  seurat_obj,
  features = score_features,
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(title = "Functional Scores across Clusters")

ggsave(
  filename = file.path(output_dir, "Cluster3_08_functional_scores_violin.pdf"),
  plot = p_scores_violin,
  width = 15,
  height = 10
)

# ===== 1.2.3 评分统计检验 =====

# 提取Cluster 3的评分
cluster3_cells <- WhichCells(seurat_obj, idents = 3)
other_cells <- setdiff(colnames(seurat_obj), cluster3_cells)

score_comparison <- data.frame(
  Score_Type = character(),
  Cluster3_Mean = numeric(),
  Cluster3_Median = numeric(),
  Others_Mean = numeric(),
  Others_Median = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (score in score_features) {
  cluster3_values <- seurat_obj@meta.data[cluster3_cells, score]
  other_values <- seurat_obj@meta.data[other_cells, score]

  wilcox_test <- wilcox.test(cluster3_values, other_values)

  score_comparison <- rbind(
    score_comparison,
    data.frame(
      Score_Type = score,
      Cluster3_Mean = mean(cluster3_values, na.rm = TRUE),
      Cluster3_Median = median(cluster3_values, na.rm = TRUE),
      Others_Mean = mean(other_values, na.rm = TRUE),
      Others_Median = median(other_values, na.rm = TRUE),
      P_value = wilcox_test$p.value
    )
  )
}

# 保存统计结果
write.csv(
  score_comparison,
  file = file.path(output_dir, "Cluster3_score_statistics.csv"),
  row.names = FALSE
)

print("=== Cluster 3 Functional Score Statistics ===")
print(score_comparison)

# ===== 1.3.1 提取三个分化方向的评分 =====

cluster3_metadata <- seurat_obj@meta.data[cluster3_cells, ]

# 创建三角图数据
differentiation_data <- data.frame(
  Cell = rownames(cluster3_metadata),
  Club_Score = cluster3_metadata$Club_Maturity1,
  Goblet_Score = cluster3_metadata$Goblet_Tendency1,
  Ciliated_Score = cluster3_metadata$Ciliated_Tendency1,
  Progenitor_Score = cluster3_metadata$Progenitor_Score1
)

# 归一化到0-1（用于三角图）
normalize_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

differentiation_data$Club_norm <- normalize_01(differentiation_data$Club_Score)
differentiation_data$Goblet_norm <- normalize_01(
  differentiation_data$Goblet_Score
)
differentiation_data$Ciliated_norm <- normalize_01(
  differentiation_data$Ciliated_Score
)

# ===== 1.3.2 散点图：Club vs Goblet vs Ciliated =====

# Club vs Goblet
p_club_goblet <- ggplot(
  differentiation_data,
  aes(x = Club_Score, y = Goblet_Score, color = Progenitor_Score)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = median(differentiation_data$Progenitor_Score)
  ) +
  theme_minimal() +
  labs(
    title = "Club vs Goblet Differentiation",
    subtitle = "Color = Progenitor Score",
    x = "Club Maturity Score",
    y = "Goblet Tendency Score"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

# Club vs Ciliated
p_club_ciliated <- ggplot(
  differentiation_data,
  aes(x = Club_Score, y = Ciliated_Score, color = Progenitor_Score)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = median(differentiation_data$Progenitor_Score)
  ) +
  theme_minimal() +
  labs(
    title = "Club vs Ciliated Differentiation",
    subtitle = "Color = Progenitor Score",
    x = "Club Maturity Score",
    y = "Ciliated Tendency Score"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

# Goblet vs Ciliated
p_goblet_ciliated <- ggplot(
  differentiation_data,
  aes(x = Goblet_Score, y = Ciliated_Score, color = Progenitor_Score)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = median(differentiation_data$Progenitor_Score)
  ) +
  theme_minimal() +
  labs(
    title = "Goblet vs Ciliated Differentiation",
    subtitle = "Color = Progenitor Score",
    x = "Goblet Tendency Score",
    y = "Ciliated Tendency Score"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

# 合并图
p_diff_combined <- (p_club_goblet | p_club_ciliated | p_goblet_ciliated) +
  plot_annotation(title = "Cluster 3 Differentiation Potential Analysis")

ggsave(
  filename = file.path(output_dir, "Cluster3_09_differentiation_potential.pdf"),
  plot = p_diff_combined,
  width = 18,
  height = 6
)

# ===== 1.3.3 分类细胞的分化倾向 =====

# 定义分类标准（基于评分中位数）
median_club <- median(differentiation_data$Club_Score)
median_goblet <- median(differentiation_data$Goblet_Score)
median_ciliated <- median(differentiation_data$Ciliated_Score)

differentiation_data$Differentiation_Type <- case_when(
  differentiation_data$Club_Score > median_club &
    differentiation_data$Goblet_Score <= median_goblet &
    differentiation_data$Ciliated_Score <= median_ciliated ~ "Club_dominant",

  differentiation_data$Goblet_Score > median_goblet &
    differentiation_data$Club_Score <= median_club ~ "Goblet_biased",

  differentiation_data$Ciliated_Score > median_ciliated &
    differentiation_data$Club_Score <= median_club ~ "Ciliated_biased",

  differentiation_data$Club_Score > median_club &
    differentiation_data$Goblet_Score >
      median_goblet ~ "Club_Goblet_transition",

  TRUE ~ "Undifferentiated_progenitor"
)

# 统计各类型比例
diff_type_summary <- differentiation_data %>%
  group_by(Differentiation_Type) %>%
  summarise(
    Count = n(),
    Percentage = n() / nrow(differentiation_data) * 100,
    Mean_Progenitor_Score = mean(Progenitor_Score, na.rm = TRUE)
  )

print("=== Cluster 3 Differentiation Type Distribution ===")
print(diff_type_summary)

# 保存
write.csv(
  diff_type_summary,
  file = file.path(output_dir, "Cluster3_differentiation_types.csv"),
  row.names = FALSE
)

# 饼图可视化
p_diff_pie <- ggplot(
  diff_type_summary,
  aes(x = "", y = Percentage, fill = Differentiation_Type)
) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(
    title = "Cluster 3 Differentiation State Distribution",
    fill = "Differentiation Type"
  ) +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    position = position_stack(vjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "Cluster3_10_differentiation_pie.pdf"),
  plot = p_diff_pie,
  width = 8,
  height = 6
)


# ===== 2.1.1 定义关键marker基因 =====

# 浆细胞核心marker
plasma_cell_genes <- c("IGKC", "IGLC2", "IGLC3", "IGHA1", "JCHAIN")

# 异常marker (需要解释的)
anomalous_genes <- c("FABP4", "VIM", "CCL14", "IGFBP4")

# Club细胞marker (检查是否混入)
club_check_genes <- c("SCGB1A1", "BPIFA1")

# HLA-II类分子
hla_genes <- c("HLA-DRA", "HLA-DPB1", "HLA-DPA1", "HLA-E")

# ===== 2.1.2 UMAP可视化 =====

# Cluster 4位置
p_cluster4_location <- DimPlot(
  seurat_obj,
  cells.highlight = WhichCells(seurat_obj, idents = 4),
  cols.highlight = "blue",
  cols = "grey90",
  pt.size = 0.5
) +
  ggtitle("Cluster 4 Location") +
  theme_minimal()

# 浆细胞marker
p_plasma_markers <- FeaturePlot(
  seurat_obj,
  features = plasma_cell_genes,
  ncol = 3,
  pt.size = 0.5,
  order = TRUE,
  cols = c("lightgrey", "darkblue")
) +
  plot_annotation(title = "Plasma Cell Markers")

# 异常marker
p_anomalous_markers <- FeaturePlot(
  seurat_obj,
  features = anomalous_genes,
  ncol = 2,
  pt.size = 0.5,
  order = TRUE,
  cols = c("lightgrey", "red")
) +
  plot_annotation(title = "Anomalous Markers (Need Explanation)")

# 保存
ggsave(
  filename = file.path(output_dir, "Cluster4_01_location.pdf"),
  plot = p_cluster4_location,
  width = 6,
  height = 5
)
ggsave(
  filename = file.path(output_dir, "Cluster4_02_plasma_markers.pdf"),
  plot = p_plasma_markers,
  width = 12,
  height = 8
)
ggsave(
  filename = file.path(output_dir, "Cluster4_03_anomalous_markers.pdf"),
  plot = p_anomalous_markers,
  width = 8,
  height = 8
)

# ===== 2.1.3 关键：共表达散点图 =====

# 提取Cluster 4细胞
cluster4_cells <- WhichCells(seurat_obj, idents = 4)

# IGKC vs FABP4
p_igkc_fabp4 <- FeatureScatter(
  seurat_obj,
  cells = cluster4_cells,
  feature1 = "IGKC",
  feature2 = "FABP4",
  pt.size = 1.5
) +
  ggtitle("Cluster 4: IGKC vs FABP4 Co-expression") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "High both = Doublet?\nHigh IGKC only = Plasma\nHigh FABP4 only = Contaminant?",
    hjust = 1.1,
    vjust = 1.1,
    size = 3,
    color = "blue"
  )

# IGKC vs VIM
p_igkc_vim <- FeatureScatter(
  seurat_obj,
  cells = cluster4_cells,
  feature1 = "IGKC",
  feature2 = "VIM",
  pt.size = 1.5
) +
  ggtitle("Cluster 4: IGKC vs VIM Co-expression")

# JCHAIN vs FABP4
p_jchain_fabp4 <- FeatureScatter(
  seurat_obj,
  cells = cluster4_cells,
  feature1 = "JCHAIN",
  feature2 = "FABP4",
  pt.size = 1.5
) +
  ggtitle("Cluster 4: JCHAIN vs FABP4 Co-expression")

# VIM vs FABP4 (检查是否是同一类细胞)
p_vim_fabp4 <- FeatureScatter(
  seurat_obj,
  cells = cluster4_cells,
  feature1 = "VIM",
  feature2 = "FABP4",
  pt.size = 1.5
) +
  ggtitle("Cluster 4: VIM vs FABP4 Co-expression")

# 合并
p_coexpression <- (p_igkc_fabp4 | p_igkc_vim) /
  (p_jchain_fabp4 | p_vim_fabp4) +
  plot_annotation(title = "Cluster 4: Critical Co-expression Analysis")

ggsave(
  filename = file.path(output_dir, "Cluster4_04_coexpression_analysis.pdf"),
  plot = p_coexpression,
  width = 14,
  height = 12
)

# ===== 2.1.4 计算相关性 =====

cluster4_expr <- GetAssayData(seurat_obj, slot = "data")[, cluster4_cells]

correlation_results <- data.frame(
  Gene1 = character(),
  Gene2 = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

gene_pairs <- list(
  c("IGKC", "FABP4"),
  c("IGKC", "VIM"),
  c("IGKC", "CCL14"),
  c("JCHAIN", "FABP4"),
  c("VIM", "FABP4"),
  c("IGKC", "SCGB1A1")
)

for (pair in gene_pairs) {
  gene1 <- pair[1]
  gene2 <- pair[2]

  if (gene1 %in% rownames(cluster4_expr) & gene2 %in% rownames(cluster4_expr)) {
    cor_test <- cor.test(
      as.numeric(cluster4_expr[gene1, ]),
      as.numeric(cluster4_expr[gene2, ]),
      method = "spearman"
    )

    correlation_results <- rbind(
      correlation_results,
      data.frame(
        Gene1 = gene1,
        Gene2 = gene2,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value
      )
    )
  }
}

cat("\n=== Cluster 4 Gene Correlation Analysis ===\n")
print(correlation_results)

write.csv(
  correlation_results,
  file = file.path(output_dir, "Cluster4_gene_correlations.csv"),
  row.names = FALSE
)

# ===== 2.3.1 Cluster 4的子聚类 =====

cluster4_subset <- subset(seurat_obj, idents = 4)

# 重新处理
cluster4_subset <- ScaleData(cluster4_subset)
cluster4_subset <- RunPCA(
  cluster4_subset,
  features = VariableFeatures(cluster4_subset)
)

# Elbow plot
p_elbow <- ElbowPlot(cluster4_subset, ndims = 50)
ggsave(
  filename = file.path(output_dir, "Cluster4_07_elbow_plot.pdf"),
  plot = p_elbow,
  width = 6,
  height = 4
)

# 聚类 (尝试不同resolution)
cluster4_subset <- FindNeighbors(cluster4_subset, dims = 1:30)

for (res in c(0.3, 0.5, 0.8)) {
  cluster4_subset <- FindClusters(cluster4_subset, resolution = res)

  cluster_col <- paste0("RNA_snn_res.", res)

  cat(sprintf(
    "\n=== Resolution %.1f: %d subclusters ===\n",
    res,
    length(unique(cluster4_subset@meta.data[[cluster_col]]))
  ))

  print(table(cluster4_subset@meta.data[[cluster_col]]))
}

# 使用最合适的resolution (这里假设0.5)
Idents(cluster4_subset) <- "RNA_snn_res.0.5"
cluster4_subset <- RunUMAP(cluster4_subset, dims = 1:30)

# 可视化子cluster
p_subcluster <- DimPlot(cluster4_subset, label = TRUE, pt.size = 1.5) +
  ggtitle("Cluster 4 Sub-clustering (Resolution 0.5)")

ggsave(
  filename = file.path(output_dir, "Cluster4_08_subclusters.pdf"),
  plot = p_subcluster,
  width = 8,
  height = 6
)

# ===== 2.3.2 子cluster的marker表达 =====

# 关键marker在子cluster上的表达
all_key_markers <- c(plasma_cell_genes, anomalous_genes)

p_subcluster_markers <- FeaturePlot(
  cluster4_subset,
  features = all_key_markers,
  ncol = 3,
  pt.size = 1,
  order = TRUE
) +
  plot_annotation(title = "Cluster 4 Subcluster: Key Marker Expression")

ggsave(
  filename = file.path(output_dir, "Cluster4_09_subcluster_markers.pdf"),
  plot = p_subcluster_markers,
  width = 15,
  height = 12
)

# DotPlot
p_subcluster_dot <- DotPlot(
  cluster4_subset,
  features = all_key_markers,
  group.by = "RNA_snn_res.0.5"
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cluster 4 Subcluster Marker Expression")

ggsave(
  filename = file.path(output_dir, "Cluster4_10_subcluster_dotplot.pdf"),
  plot = p_subcluster_dot,
  width = 10,
  height = 8
)

# ===== 2.3.3 找每个子cluster的DE genes =====

subcluster_markers <- FindAllMarkers(
  cluster4_subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Top 20 markers per subcluster
top20_subcluster <- subcluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(
  top20_subcluster,
  file = file.path(output_dir, "Cluster4_subcluster_top20_markers.csv"),
  row.names = FALSE
)

# Heatmap
top10_subcluster <- subcluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_subcluster_heatmap <- DoHeatmap(
  cluster4_subset,
  features = top10_subcluster$gene,
  group.by = "RNA_snn_res.0.5"
) +
  theme(axis.text.y = element_text(size = 6))

ggsave(
  filename = file.path(output_dir, "Cluster4_11_subcluster_heatmap.pdf"),
  plot = p_subcluster_heatmap,
  width = 12,
  height = 15
)

# ===== 2.3.4 子cluster的细胞类型预测 =====

# 基于marker判断每个subcluster的身份
subcluster_identity_summary <- data.frame(
  Subcluster = character(),
  Cell_Count = numeric(),
  IGKC_pct = numeric(),
  FABP4_pct = numeric(),
  VIM_pct = numeric(),
  Predicted_Identity = character(),
  stringsAsFactors = FALSE
)

for (sc in unique(Idents(cluster4_subset))) {
  sc_cells <- WhichCells(cluster4_subset, idents = sc)

  igkc_pct <- mean(cluster4_subset[["RNA"]]@layers@data["IGKC", sc_cells] > 0) *
    100
  fabp4_pct <- mean(
    cluster4_subset[["RNA"]]@layers@data["FABP4", sc_cells] > 0
  ) *
    100
  vim_pct <- mean(cluster4_subset[["RNA"]]@layers@adta["VIM", sc_cells] > 0) *
    100

  # 简单规则判断
  predicted_id <- case_when(
    igkc_pct > 70 & fabp4_pct < 30 ~ "Plasma_Cell",
    fabp4_pct > 50 & igkc_pct < 30 ~ "Macrophage_or_Adipocyte",
    vim_pct > 60 & igkc_pct < 30 ~ "Fibroblast",
    igkc_pct > 40 & fabp4_pct > 40 ~ "Likely_Doublet",
    TRUE ~ "Mixed_or_Uncertain"
  )

  subcluster_identity_summary <- rbind(
    subcluster_identity_summary,
    data.frame(
      Subcluster = sc,
      Cell_Count = length(sc_cells),
      IGKC_pct = round(igkc_pct, 1),
      FABP4_pct = round(fabp4_pct, 1),
      VIM_pct = round(vim_pct, 1),
      Predicted_Identity = predicted_id
    )
  )
}

cat("\n=== Cluster 4 Subcluster Identity Prediction ===\n")
print(subcluster_identity_summary)

write.csv(
  subcluster_identity_summary,
  file = file.path(output_dir, "Cluster4_subcluster_identity.csv"),
  row.names = FALSE
)


# ===== 1.5.1 Cluster 3的子聚类 =====

# 提取Cluster 3
cluster3_subset <- subset(seurat_obj, idents = 3)

# 重新进行PCA和聚类
cluster3_subset <- ScaleData(
  cluster3_subset,
  features = rownames(cluster3_subset)
)
cluster3_subset <- RunPCA(
  cluster3_subset,
  features = VariableFeatures(cluster3_subset)
)

# 选择PC数量（查看Elbow plot）
ElbowPlot(cluster3_subset, ndims = 50)

# 聚类
cluster3_subset <- FindNeighbors(cluster3_subset, dims = 1:30)
cluster3_subset <- FindClusters(cluster3_subset, resolution = 0.5)
cluster3_subset <- RunUMAP(cluster3_subset, dims = 1:30)

# 可视化子cluster
p_subcluster_umap <- DimPlot(
  cluster3_subset,
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1
) +
  ggtitle("Cluster 3 Sub-clustering")

ggsave(
  filename = file.path(output_dir, "Cluster3_13_subclusters.pdf"),
  plot = p_subcluster_umap,
  width = 8,
  height = 6
)

# ===== 1.5.2 子cluster的marker表达 =====

# 在子cluster上展示关键marker
p_subcluster_markers <- FeaturePlot(
  cluster3_subset,
  features = c("ALDH1A1", "SCGB1A1", "MUC5B", "GRHL2"),
  ncol = 2,
  pt.size = 1,
  order = TRUE
) +
  plot_annotation(title = "Cluster 3 Subcluster Markers")

ggsave(
  filename = file.path(output_dir, "Cluster3_14_subcluster_markers.pdf"),
  plot = p_subcluster_markers,
  width = 10,
  height = 10
)

# ===== 1.5.3 找每个子cluster的marker =====

Idents(cluster3_subset) <- "seurat_clusters"
subcluster_markers <- FindAllMarkers(
  cluster3_subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 保存top markers
top_subcluster_markers <- subcluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(
  top_subcluster_markers,
  file = file.path(output_dir, "Cluster3_subcluster_top20_markers.csv"),
  row.names = FALSE
)

# Heatmap展示
top10 <- subcluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_subcluster_heatmap <- DoHeatmap(
  cluster3_subset,
  features = top10$gene,
  group.by = "seurat_clusters"
) +
  theme(axis.text.y = element_text(size = 6))

ggsave(
  filename = file.path(output_dir, "Cluster3_15_subcluster_heatmap.pdf"),
  plot = p_subcluster_heatmap,
  width = 12,
  height = 15
)


p10 <- p_subcluster_dot$data
p10


# ===== 1. 重新标注 Cluster 4 =====

# 提取subcluster信息
cluster4_meta <- cluster4_subset@meta.data

# 基于正确的细胞类型重新标注
cluster4_meta$Cell_Type_Corrected <- case_when(
  cluster4_meta$RNA_snn_res.0.5 == 0 ~ "Inflammatory_Myofibroblasts_IL33",
  cluster4_meta$RNA_snn_res.0.5 == 1 ~ "Macrophages_AIF1",
  cluster4_meta$RNA_snn_res.0.5 == 2 ~ "Perivascular_Fibroblasts_FBN1",
  TRUE ~ "Uncertain"
)

# 添加到原始对象
seurat_obj@meta.data[rownames(cluster4_meta), "Cluster4_Corrected"] <-
  cluster4_meta$Cell_Type_Corrected

# ===== 2. 验证关键marker =====

# IL33在Sub0的表达
VlnPlot(
  cluster4_subset,
  features = c("IL33", "TAGLN", "CNN3"),
  group.by = "RNA_snn_res.0.5",
  pt.size = 0
) +
  ggtitle("Myofibroblast Markers in Subcluster 0")

# AIF1在Sub1的表达
VlnPlot(
  cluster4_subset,
  features = c("AIF1", "S100A4", "CD68"),
  group.by = "RNA_snn_res.0.5",
  pt.size = 0
) +
  ggtitle("Macrophage Markers in Subcluster 1")

# FBN1在Sub2的表达
VlnPlot(
  cluster4_subset,
  features = c("FBN1", "ITGA9", "CACNA1C"),
  group.by = "RNA_snn_res.0.5",
  pt.size = 0
) +
  ggtitle("Perivascular Markers in Subcluster 2")

# ===== 3. 验证"吞噬抗体"假说 =====

# 双染：AIF1 vs IGKC
FeatureScatter(
  cluster4_subset,
  feature1 = "AIF1",
  feature2 = "IGKC",
  group.by = "RNA_snn_res.0.5"
) +
  ggtitle("AIF1 vs IGKC: Phagocytosed Antibodies?")

# 如果Sub1是吞噬抗体的巨噬细胞：
# - AIF1高的细胞，IGKC也高
# - 应该看到正相关

# ===== 4. IL33表达的重要性分析 =====

# 计算IL33评分
il33_expressing_cells <- WhichCells(cluster4_subset, expression = IL33 > 0)

cat(sprintf(
  "\nIL33+ cells in Cluster 4: %d (%.1f%%)\n",
  length(il33_expressing_cells),
  100 * length(il33_expressing_cells) / ncol(cluster4_subset)
))

# IL33在各subcluster的比例
il33_by_subcluster <- cluster4_subset@meta.data %>%
  group_by(RNA_snn_res.0.5) %>%
  summarise(
    Total = n(),
    IL33_positive = sum(cluster4_subset[["RNA"]]@layers@data["IL33", ] > 0),
    IL33_pct = 100 * IL33_positive / Total,
    IL33_mean = mean(cluster4_subset[["RNA"]]@layers@data["IL33", ])
  )

print(il33_by_subcluster)

# ===== 5. 更新最终cluster注释 =====

seurat_obj$cluster_final_corrected <- as.character(seurat_obj$seurat_clusters)

# 替换Cluster 4
seurat_obj$cluster_final_corrected[seurat_obj$seurat_clusters == 4] <-
  paste0("C4_", seurat_obj$Cluster4_Corrected[seurat_obj$seurat_clusters == 4])

# 查看
table(seurat_obj$cluster_final_corrected)

# ===== 6. 可视化修正 =====

# UMAP with corrected annotation
p_corrected <- DimPlot(
  seurat_obj,
  group.by = "cluster_final_corrected",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("Corrected Cluster Annotation - Cluster 4 Identified as Stromal")

ggsave("Cluster4_CORRECTED_annotation.pdf", p_corrected, width = 12, height = 8)

# 关键marker的FeaturePlot
key_markers <- c("IL33", "AIF1", "FBN1", "VIM", "IGKC", "FABP4")

p_key_markers <- FeaturePlot(
  seurat_obj,
  features = key_markers,
  ncol = 3,
  order = TRUE,
  pt.size = 0.3
) +
  plot_annotation(title = "Cluster 4 Corrected Identity Markers")

ggsave(
  "Cluster4_corrected_key_markers.pdf",
  p_key_markers,
  width = 15,
  height = 10
)

# ===== 7. 保存修正的对象 =====

saveRDS(seurat_obj, "seurat_obj_cluster4_CORRECTED.rds")


# ===== 1.1 定义需要验证的关键marker =====

# Cluster 3 (Club细胞)
club_markers <- c(
  "SCGB1A1",
  "BPIFA1",
  "GRHL2",
  "EHF",
  "ALDH1A1",
  "MSI2",
  "SOX2",
  "MUC5B"
)

# Cluster 4 Subcluster 0 (假设的肌成纤维细胞)
myofib_markers <- c(
  "VIM",
  "FABP4",
  "IL33",
  "TAGLN",
  "CNN3",
  "CAV1",
  "ADIRF",
  "CCL14"
)

# Cluster 4 Subcluster 1 (假设的巨噬细胞/浆细胞)
macro_plasma_markers <- c(
  "AIF1",
  "S100A4",
  "CD68",
  "FCGR3A", # 巨噬细胞
  "IGKC",
  "JCHAIN",
  "IGLC2",
  "IGHA1"
) # 浆细胞

# Cluster 4 Subcluster 2 (假设的血管周成纤维细胞)
perivascular_markers <- c("FBN1", "ITGA9", "CACNA1C", "EPHA4", "RBMS3")

# 合并所有marker
all_key_markers <- c(
  club_markers,
  myofib_markers,
  macro_plasma_markers,
  perivascular_markers
)

# 去重
all_key_markers <- unique(all_key_markers)

# 检查哪些marker存在于数据中
markers_present <- all_key_markers[all_key_markers %in% rownames(seurat_obj)]

cat("=== Markers Present in Dataset ===\n")
cat(sprintf("Total requested: %d\n", length(all_key_markers)))
cat(sprintf("Actually present: %d\n", length(markers_present)))
cat("\nMissing markers:\n")
print(setdiff(all_key_markers, markers_present))

# ===== 1.2 提取绝对表达量 =====

# 使用normalized data (log1p transformed)
expr_data <- GetAssayData(seurat_obj, slot = "data")

# 创建汇总表
cluster_expr_summary <- data.frame()

for (marker in markers_present) {
  for (cluster_id in sort(unique(seurat_obj$seurat_clusters))) {
    # 获取该cluster的细胞
    cluster_cells <- WhichCells(seurat_obj, idents = cluster_id)

    # 提取该marker在该cluster的表达
    expr_values <- expr_data[marker, cluster_cells]

    # 计算统计量
    cluster_expr_summary <- rbind(
      cluster_expr_summary,
      data.frame(
        Marker = marker,
        Cluster = cluster_id,
        Mean_Expr = mean(expr_values),
        Median_Expr = median(expr_values),
        Max_Expr = max(expr_values),
        Pct_Positive = 100 * mean(expr_values > 0), # 表达>0的细胞比例
        Pct_High = 100 * mean(expr_values > 1), # 表达>1的细胞比例
        SD_Expr = sd(expr_values),
        Cell_Count = length(cluster_cells),
        stringsAsFactors = FALSE
      )
    )
  }
}

# 保存完整表格
write.csv(
  cluster_expr_summary,
  file.path(output_dir, "Cluster_Absolute_Expression_Summary.csv"),
  row.names = FALSE
)

cat("\n=== Expression Summary Saved ===\n")
head(cluster_expr_summary, 20)

# ===== 1.3 为Cluster 4创建subcluster版本 =====

# 如果已经有cluster4_subset和subcluster信息
if (exists("cluster4_subset")) {
  cluster4_expr_data <- GetAssayData(cluster4_subset, slot = "data")

  cluster4_subcluster_summary <- data.frame()

  for (marker in markers_present) {
    for (subcluster_id in sort(unique(cluster4_subset$RNA_snn_res.0.5))) {
      subcluster_cells <- WhichCells(cluster4_subset, idents = subcluster_id)

      expr_values <- cluster4_expr_data[marker, subcluster_cells]

      cluster4_subcluster_summary <- rbind(
        cluster4_subcluster_summary,
        data.frame(
          Marker = marker,
          Subcluster = subcluster_id,
          Mean_Expr = mean(expr_values),
          Median_Expr = median(expr_values),
          Max_Expr = max(expr_values),
          Pct_Positive = 100 * mean(expr_values > 0),
          Pct_High = 100 * mean(expr_values > 1),
          SD_Expr = sd(expr_values),
          Cell_Count = length(subcluster_cells),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  write.csv(
    cluster4_subcluster_summary,
    file.path(output_dir, "Cluster4_Subcluster_Absolute_Expression.csv"),
    row.names = FALSE
  )

  cat("\n=== Cluster 4 Subcluster Expression Summary Saved ===\n")
  head(cluster4_subcluster_summary, 20)
}


# ===== 2.1 VlnPlot - 显示绝对表达分布 =====

# 为每组marker创建VlnPlot

# Club细胞markers
p_club_abs <- VlnPlot(
  seurat_obj,
  features = club_markers[club_markers %in% rownames(seurat_obj)],
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 4,
  log = FALSE
) + # 不对Y轴log转换，显示真实值
  plot_annotation(
    title = "Club Cell Markers - Absolute Expression (log1p normalized)",
    subtitle = "Higher values = higher absolute expression"
  )

ggsave(
  file.path(output_dir, "Absolute_Expr_Club_markers.pdf"),
  p_club_abs,
  width = 16,
  height = 8
)

# 肌成纤维细胞markers
p_myofib_abs <- VlnPlot(
  seurat_obj,
  features = myofib_markers[myofib_markers %in% rownames(seurat_obj)],
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 4
) +
  plot_annotation(title = "Myofibroblast Markers - Absolute Expression")

ggsave(
  file.path(output_dir, "Absolute_Expr_Myofibroblast_markers.pdf"),
  p_myofib_abs,
  width = 16,
  height = 8
)

# 巨噬细胞/浆细胞markers
p_macro_plasma_abs <- VlnPlot(
  seurat_obj,
  features = macro_plasma_markers[
    macro_plasma_markers %in% rownames(seurat_obj)
  ],
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 4
) +
  plot_annotation(
    title = "Macrophage/Plasma Cell Markers - Absolute Expression"
  )

ggsave(
  file.path(output_dir, "Absolute_Expr_Macrophage_Plasma_markers.pdf"),
  p_macro_plasma_abs,
  width = 16,
  height = 8
)

# ===== 2.2 针对Cluster 4的subcluster对比 =====

if (exists("cluster4_subset")) {
  # 关键marker在Cluster 4 subclusters的表达
  key_cluster4_markers <- c(
    "IGKC",
    "JCHAIN",
    "AIF1",
    "S100A4",
    "VIM",
    "FABP4",
    "IL33",
    "TAGLN",
    "FBN1",
    "ITGA9"
  )

  key_cluster4_markers <- key_cluster4_markers[
    key_cluster4_markers %in% rownames(cluster4_subset)
  ]

  p_cluster4_sub_abs <- VlnPlot(
    cluster4_subset,
    features = key_cluster4_markers,
    group.by = "RNA_snn_res.0.5",
    pt.size = 0,
    ncol = 5
  ) +
    plot_annotation(
      title = "Cluster 4 Subclusters - Absolute Expression Comparison"
    )

  ggsave(
    file.path(output_dir, "Absolute_Expr_Cluster4_Subclusters.pdf"),
    p_cluster4_sub_abs,
    width = 20,
    height = 8
  )
}

# ===== 2.3 Heatmap - 所有cluster的关键marker绝对表达 =====

# 创建平均表达矩阵
expr_matrix <- cluster_expr_summary %>%
  select(Marker, Cluster, Mean_Expr) %>%
  pivot_wider(names_from = Cluster, values_from = Mean_Expr) %>%
  column_to_rownames("Marker") %>%
  as.matrix()

# 绘制heatmap
library(pheatmap)

# 分类marker并排序
marker_annotation <- data.frame(
  Type = c(
    rep("Club", length(club_markers)),
    rep("Myofibroblast", length(myofib_markers)),
    rep("Macro_Plasma", length(macro_plasma_markers)),
    rep("Perivascular", length(perivascular_markers))
  ),
  row.names = c(
    club_markers,
    myofib_markers,
    macro_plasma_markers,
    perivascular_markers
  )
)

# 只保留存在的marker
marker_annotation <- marker_annotation[
  rownames(marker_annotation) %in% rownames(expr_matrix),
  ,
  drop = FALSE
]

pdf(
  file.path(output_dir, "Absolute_Expr_Heatmap_AllClusters.pdf"),
  width = 12,
  height = 16
)
pheatmap(
  expr_matrix[rownames(marker_annotation), ],
  scale = "row", # 按行标准化以便比较
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = marker_annotation,
  main = "Marker Absolute Expression Across All Clusters\n(Row-scaled)",
  fontsize = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()

# 不标准化的版本（显示真实绝对值）
pdf(
  file.path(output_dir, "Absolute_Expr_Heatmap_AllClusters_Unscaled.pdf"),
  width = 12,
  height = 16
)
pheatmap(
  expr_matrix[rownames(marker_annotation), ],
  scale = "none", # 不标准化
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = marker_annotation,
  main = "Marker Absolute Expression Across All Clusters\n(Unscaled)",
  fontsize = 8,
  color = colorRampPalette(c("lightgrey", "yellow", "red"))(100)
)
dev.off()


# ===== 提取关键间质/免疫marker =====

key_markers_needed <- c(
  # 间质/成纤维细胞
  "VIM",
  "FABP4",
  "IL33",
  "TAGLN",
  "CNN3",
  "CAV1",
  "CCL14",
  # 巨噬细胞
  "AIF1",
  "S100A4",
  "CD68",
  "FCGR3A",
  "CD14",
  # 浆细胞
  "IGKC",
  "JCHAIN",
  "IGLC2",
  "IGHA1",
  "CD38",
  "SDC1",
  # 血管周
  "FBN1",
  "ITGA9",
  "CACNA1C",
  # 上皮/Club
  "EPCAM",
  "CDH1",
  "KRT8",
  "KRT18"
)

# 检查存在的marker
key_markers_present <- key_markers_needed[
  key_markers_needed %in% rownames(seurat_obj)
]

# 提取表达数据
key_expr_summary <- cluster_expr_summary %>%
  filter(Marker %in% key_markers_present) %>%
  filter(Cluster %in% c(3, 4)) %>%
  select(Marker, Cluster, Mean_Expr, Median_Expr, Pct_Positive, Pct_High) %>%
  arrange(Marker, Cluster)

print(key_expr_summary)

# Cluster 4 subclusters
if (exists("cluster4_subcluster_summary")) {
  key_cluster4_summary <- cluster4_subcluster_summary %>%
    filter(Marker %in% key_markers_present) %>%
    select(
      Marker,
      Subcluster,
      Mean_Expr,
      Median_Expr,
      Pct_Positive,
      Pct_High
    ) %>%
    arrange(Marker, Subcluster)

  print(key_cluster4_summary)
}

# 保存
write.csv(key_expr_summary, "Key_Markers_Cluster3_4_Expression.csv")
if (exists("key_cluster4_summary")) {
  write.csv(key_cluster4_summary, "Key_Markers_Cluster4_Subclusters.csv")
}


# ===== 1.1 基于已有数据识别细胞类型 =====

library(tidyverse)
library(Seurat)

# 创建细胞类型分类
seurat_obj$Cell_Type_Broad <- NA

# 先看看所有cluster的关键marker表达
# 我们需要先确定哪些是成纤维，哪些是上皮

# 提取关键marker的表达（所有cluster）
key_markers_all_clusters <- c(
  "VIM",
  "COL1A1",
  "COL3A1",
  "DCN", # 成纤维marker
  "EPCAM",
  "CDH1",
  "KRT8",
  "KRT18",
  "KRT19", # 上皮marker
  "GRHL2",
  "EHF", # 上皮TF
  "BPIFA1",
  "SCGB1A1",
  "MUC5AC",
  "FOXJ1",
  "TP63"
) # 上皮亚型

# 检查哪些marker存在
key_markers_present <- key_markers_all_clusters[
  key_markers_all_clusters %in% rownames(seurat_obj)
]

# 计算每个cluster的平均表达
cluster_marker_expr <- data.frame()

for (cluster in sort(unique(seurat_obj$seurat_clusters))) {
  cluster_cells <- WhichCells(seurat_obj, idents = cluster)

  for (marker in key_markers_present) {
    expr_values <- GetAssayData(seurat_obj, slot = "data")[
      marker,
      cluster_cells
    ]

    cluster_marker_expr <- rbind(
      cluster_marker_expr,
      data.frame(
        Cluster = cluster,
        Marker = marker,
        Mean_Expr = mean(expr_values),
        Pct_Positive = 100 * mean(expr_values > 0),
        stringsAsFactors = FALSE
      )
    )
  }
}

# 转换为宽表格便于查看
cluster_marker_wide <- cluster_marker_expr %>%
  select(Cluster, Marker, Mean_Expr) %>%
  pivot_wider(names_from = Marker, values_from = Mean_Expr)

print("=== Key Markers Across All Clusters ===")
print(cluster_marker_wide)

# 保存
write.csv(
  cluster_marker_wide,
  file.path(output_dir, "All_Clusters_Key_Markers.csv"),
  row.names = FALSE
)

# ===== 1.2 自动分类细胞类型 =====

# 定义分类规则
classify_cluster <- function(cluster_data) {
  # cluster_data是一个包含marker表达的行

  # 获取关键marker值
  vim <- ifelse("VIM" %in% names(cluster_data), cluster_data$VIM, 0)
  epcam <- ifelse("EPCAM" %in% names(cluster_data), cluster_data$EPCAM, 0)
  cdh1 <- ifelse("CDH1" %in% names(cluster_data), cluster_data$CDH1, 0)
  krt8 <- ifelse("KRT8" %in% names(cluster_data), cluster_data$KRT8, 0)
  grhl2 <- ifelse("GRHL2" %in% names(cluster_data), cluster_data$GRHL2, 0)
  bpifa1 <- ifelse("BPIFA1" %in% names(cluster_data), cluster_data$BPIFA1, 0)

  # 分类逻辑
  if (vim > 2 & epcam < 1 & grhl2 < 0.3) {
    return("Fibroblast")
  } else if (epcam > 1 | grhl2 > 0.3 | (krt8 > 1 & cdh1 > 0.5)) {
    return("Epithelial")
  } else if (bpifa1 > 3 & grhl2 > 0.2) {
    return("Epithelial")
  } else if (vim > 1 & bpifa1 > 2 & grhl2 < 0.2) {
    return("Ambiguous_High_BPIFA1") # 可能是ambient RNA
  } else {
    return("Other")
  }
}

# 应用分类
cluster_classification <- cluster_marker_wide %>%
  rowwise() %>%
  mutate(Cell_Type_Broad = classify_cluster(cur_data()))

cat("\n=== Automatic Cluster Classification ===\n")
print(cluster_classification %>% select(Cluster, Cell_Type_Broad))

# 保存分类结果
write.csv(
  cluster_classification,
  file.path(output_dir, "Cluster_Classification.csv"),
  row.names = FALSE
)

# ===== 1.3 手动确认和调整 =====

# 打印需要手动确认的cluster
cat("\n=== Clusters Need Manual Verification ===\n")
print(
  cluster_classification %>%
    filter(Cell_Type_Broad %in% c("Ambiguous_High_BPIFA1", "Other"))
)

# 基于我们已知的信息手动设置
# 你需要根据实际cluster数量调整

# 示例（请根据你的数据调整）：
manual_classification <- data.frame(
  Cluster = unique(seurat_obj$seurat_clusters),
  Manual_Type = c(
    "Other", # Cluster 0 - 需要确认
    "Other", # Cluster 1 - 需要确认
    "Other", # Cluster 2 - 需要确认
    "Epithelial", # Cluster 3 - Club细胞
    "Fibroblast", # Cluster 4 - 成纤维（我们已确认）
    rep("Other", length(unique(seurat_obj$seurat_clusters)) - 5) # 其余cluster
  )
)

cat(
  "\n⚠️ IMPORTANT: Please review and manually edit the classification above!\n"
)
cat("Based on your DotPlot and UMAP, identify:\n")
cat("1. Which clusters are clearly Epithelial (high EPCAM/KRT/GRHL2)?\n")
cat(
  "2. Which clusters are clearly Fibroblast (high VIM, low epithelial markers)?\n"
)
cat("3. Which are immune cells, endothelial, etc?\n\n")

# 提示用户
cat("Please run the following to see marker expression heatmap:\n")
cat(
  "DoHeatmap(seurat_obj, features = key_markers_present, group.by = 'seurat_clusters')\n\n"
)


# ===== 2.1 DotPlot - 所有cluster的关键marker =====

# 成纤维 vs 上皮分类marker
classification_markers <- c(
  # 成纤维细胞
  "VIM",
  "COL1A1",
  "COL3A1",
  "DCN",
  "PDGFRA",
  "PDGFRB",
  # 上皮细胞
  "EPCAM",
  "CDH1",
  "KRT8",
  "KRT18",
  "KRT19",
  # 上皮TF
  "GRHL2",
  "EHF",
  "TP63",
  # 上皮亚型
  "BPIFA1",
  "SCGB1A1",
  "MUC5AC",
  "FOXJ1"
)

classification_markers <- classification_markers[
  classification_markers %in% rownames(seurat_obj)
]

p_classification_dot <- DotPlot(
  seurat_obj,
  features = classification_markers,
  group.by = "seurat_clusters"
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Fibroblast vs Epithelial Classification Markers")

ggsave(
  file.path(output_dir, "Classification_DotPlot.pdf"),
  p_classification_dot,
  width = 12,
  height = 10
)

# ===== 2.2 Heatmap - 按marker类型分组 =====

library(pheatmap)

# 准备heatmap数据
heatmap_data <- cluster_marker_wide %>%
  column_to_rownames("Cluster") %>%
  select(any_of(classification_markers))

# 转置（行=marker，列=cluster）
heatmap_data_t <- t(heatmap_data)

# Marker分组注释
marker_annotation <- data.frame(
  Type = c(
    rep(
      "Fibroblast",
      sum(
        c("VIM", "COL1A1", "COL3A1", "DCN", "PDGFRA", "PDGFRB") %in%
          rownames(heatmap_data_t)
      )
    ),
    rep(
      "Epithelial_General",
      sum(
        c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19") %in%
          rownames(heatmap_data_t)
      )
    ),
    rep(
      "Epithelial_TF",
      sum(c("GRHL2", "EHF", "TP63") %in% rownames(heatmap_data_t))
    ),
    rep(
      "Epithelial_Subtype",
      sum(
        c("BPIFA1", "SCGB1A1", "MUC5AC", "FOXJ1") %in% rownames(heatmap_data_t)
      )
    )
  ),
  row.names = rownames(heatmap_data_t)
)

pdf(
  file.path(output_dir, "Classification_Heatmap.pdf"),
  width = 10,
  height = 12
)
pheatmap(
  heatmap_data_t,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = marker_annotation,
  main = "Cluster Classification Heatmap",
  fontsize = 8
)
dev.off()

# ===== 2.3 VlnPlot - VIM vs EPCAM 分布 =====

p_vim_epcam <- VlnPlot(
  seurat_obj,
  features = c("VIM", "EPCAM"),
  group.by = "seurat_clusters",
  ncol = 2,
  pt.size = 0
)

ggsave(
  file.path(output_dir, "VIM_vs_EPCAM_Violin.pdf"),
  p_vim_epcam,
  width = 12,
  height = 6
)

# ===== 2.4 Scatter plot - VIM vs GRHL2/EPCAM =====

# 计算每个cluster的平均值
cluster_avg <- cluster_marker_wide

if ("VIM" %in% colnames(cluster_avg) & "EPCAM" %in% colnames(cluster_avg)) {
  p_scatter_vim_epcam <- ggplot(
    cluster_avg,
    aes(x = VIM, y = EPCAM, label = Cluster)
  ) +
    geom_point(size = 4, alpha = 0.7) +
    geom_text(vjust = -1, size = 3) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
    annotate(
      "text",
      x = 4,
      y = 0.5,
      label = "Fibroblast\n(VIM>2, EPCAM<1)",
      color = "red"
    ) +
    annotate(
      "text",
      x = 1,
      y = 3,
      label = "Epithelial\n(EPCAM>1)",
      color = "blue"
    ) +
    theme_minimal() +
    labs(
      title = "Cluster Classification: VIM vs EPCAM",
      subtitle = "Based on mean expression"
    )

  ggsave(
    file.path(output_dir, "Scatter_VIM_vs_EPCAM.pdf"),
    p_scatter_vim_epcam,
    width = 8,
    height = 6
  )
}

if ("VIM" %in% colnames(cluster_avg) & "GRHL2" %in% colnames(cluster_avg)) {
  p_scatter_vim_grhl2 <- ggplot(
    cluster_avg,
    aes(x = VIM, y = GRHL2, label = Cluster)
  ) +
    geom_point(size = 4, alpha = 0.7) +
    geom_text(vjust = -1, size = 3) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "blue") +
    annotate(
      "text",
      x = 4,
      y = 0.1,
      label = "Fibroblast\n(VIM>2, GRHL2<0.3)",
      color = "red"
    ) +
    annotate(
      "text",
      x = 1,
      y = 0.8,
      label = "Epithelial\n(GRHL2>0.3)",
      color = "blue"
    ) +
    theme_minimal() +
    labs(
      title = "Cluster Classification: VIM vs GRHL2",
      subtitle = "Based on mean expression"
    )

  ggsave(
    file.path(output_dir, "Scatter_VIM_vs_GRHL2.pdf"),
    p_scatter_vim_grhl2,
    width = 8,
    height = 6
  )
}
# dev.off()

# ===== 1.1 加载必要的包 =====
library(Seurat)
library(tidyverse)
library(patchwork)

# 设置输出目录
output_dir <- "Final_Annotation_Results"
dir.create(output_dir, showWarnings = FALSE)

# ===== 1.2 提取Cluster 3和4的细胞 =====

# Cluster 3细胞
cluster3_cells <- WhichCells(seurat_obj, idents = 3)

# Cluster 4细胞（如果已经有subcluster信息）
cluster4_cells <- WhichCells(seurat_obj, idents = 4)

# 创建包含Cluster 3和4的subset
clusters_3_4 <- subset(seurat_obj, idents = c(3, 4))

cat(sprintf("Cluster 3: %d cells\n", length(cluster3_cells)))
cat(sprintf("Cluster 4: %d cells\n", length(cluster4_cells)))
cat(sprintf("Total: %d cells\n", ncol(clusters_3_4)))

# ===== 1.3 创建细分的cluster ID =====

# 为Cluster 3和4创建细分ID
clusters_3_4$Cluster_Detailed <- as.character(clusters_3_4$seurat_clusters)

# 如果Cluster 4有subcluster信息，添加进去
if (
  exists("cluster4_subset") &&
    "RNA_snn_res.0.5" %in% colnames(cluster4_subset@meta.data)
) {
  # 为Cluster 4的细胞添加subcluster信息
  cluster4_subcluster <- cluster4_subset$RNA_snn_res.0.5
  names(cluster4_subcluster) <- colnames(cluster4_subset)

  # 更新Cluster_Detailed
  clusters_3_4$Cluster_Detailed[
    colnames(clusters_3_4) %in% names(cluster4_subcluster)
  ] <-
    paste0(
      "C4_Sub",
      cluster4_subcluster[colnames(clusters_3_4)[
        colnames(clusters_3_4) %in% names(cluster4_subcluster)
      ]]
    )

  # Cluster 3保持原样
  clusters_3_4$Cluster_Detailed[clusters_3_4$seurat_clusters == 3] <- "C3"
} else {
  # 如果没有subcluster，就简单标记
  clusters_3_4$Cluster_Detailed <- paste0("C", clusters_3_4$seurat_clusters)
}

# 查看分组
cat("\n=== Cluster Detailed Groups ===\n")
print(table(clusters_3_4$Cluster_Detailed))

# ===== 1.4 定义所有关键Marker =====

all_markers <- list(
  # === Club/Secretory细胞 ===
  Club_Core = c("SCGB1A1", "SCGB3A1", "BPIFA1", "BPIFB1", "BPIFB2"),

  Epithelial_TF = c(
    "GRHL2",
    "GRHL1",
    "EHF",
    "ELF3",
    "ELF5",
    "OVOL2",
    "TFCP2L1"
  ),

  Progenitor = c("ALDH1A1", "MSI2", "SOX2", "SOX9", "BMI1"),

  Epithelial_General = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19"),

  Goblet_Differentiation = c("MUC5B", "MUC5AC", "SPDEF", "AGR2", "TFF3"),

  # === 成纤维细胞 ===
  Fibroblast_Core = c("VIM", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),

  Myofibroblast = c("TAGLN", "CNN3", "ACTA2", "MYH11"),

  Inflammatory_Fibro = c("IL33", "CCL14", "CXCL12", "CXCL14"),

  Lipid_Metabolism = c("FABP4", "FABP5", "ADIRF", "CAV1", "CAVIN2"),

  # === 巨噬细胞 ===
  Macrophage = c("AIF1", "S100A4", "CD68", "CD14", "FCGR3A", "CD163", "MSR1"),

  # === 浆细胞/B细胞 ===
  Plasma_Cell = c(
    "IGKC",
    "IGLC2",
    "IGLC3",
    "IGHA1",
    "IGHA2",
    "IGHG1",
    "JCHAIN",
    "SDC1",
    "CD38"
  ),

  # === 血管周成纤维细胞 ===
  Perivascular = c("FBN1", "ITGA9", "CACNA1C", "EPHA4", "RBMS3"),

  # === 线粒体/代谢 ===
  Mitochondrial = c("MTRNR2L8", "MTRNR2L12", "MT-CO1", "MT-ND1", "MT-ATP6"),

  # === 细胞周期 ===
  Proliferation = c("MKI67", "TOP2A", "PCNA")
)

# 合并所有marker（去重）
all_markers_vector <- unique(unlist(all_markers))

# 检查哪些marker存在于数据中
markers_present <- all_markers_vector[
  all_markers_vector %in% rownames(clusters_3_4)
]

cat("\n=== Marker Summary ===\n")
cat(sprintf("Total markers requested: %d\n", length(all_markers_vector)))
cat(sprintf("Markers present in data: %d\n", length(markers_present)))
cat(sprintf(
  "Missing markers: %d\n",
  length(all_markers_vector) - length(markers_present)
))

# 打印缺失的marker
missing_markers <- setdiff(all_markers_vector, markers_present)
if (length(missing_markers) > 0) {
  cat("\nMissing markers:\n")
  print(missing_markers)
}

# 保存marker列表
marker_categories <- data.frame(
  Marker = all_markers_vector,
  Category = rep(names(all_markers), sapply(all_markers, length)),
  Present = all_markers_vector %in% rownames(clusters_3_4)
)

write.csv(
  marker_categories,
  file.path(output_dir, "Marker_List_with_Categories.csv"),
  row.names = FALSE
)


# ===== 2.1 按类别组织marker（只保留存在的）=====

markers_for_plot <- list()

for (category in names(all_markers)) {
  present_in_category <- all_markers[[category]][
    all_markers[[category]] %in% markers_present
  ]
  if (length(present_in_category) > 0) {
    markers_for_plot[[category]] <- present_in_category
  }
}

# 合并成向量（保持分类顺序）
markers_ordered <- unlist(markers_for_plot)

cat("\n=== Markers for DotPlot ===\n")
cat(sprintf("Total markers to plot: %d\n", length(markers_ordered)))

# ===== 2.2 创建DotPlot =====

# 设置Idents为详细的cluster分组
Idents(clusters_3_4) <- "Cluster_Detailed"

# DotPlot (不缩放)
p_dotplot_unscaled <- DotPlot(
  clusters_3_4,
  features = markers_ordered,
  group.by = "Cluster_Detailed",
  scale = FALSE
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Cluster 3 & 4 Detailed Marker Expression",
    subtitle = "Unscaled - Dot size = % expressing, Color = mean expression",
    x = "Markers",
    y = "Cluster"
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 2,
    name = "Mean\nExpression"
  )

ggsave(
  file.path(output_dir, "DotPlot_Cluster3_4_Detailed_Unscaled.pdf"),
  p_dotplot_unscaled,
  width = 50,
  height = 20,
  limitsize = FALSE
)

# DotPlot (行标准化)
p_dotplot_scaled <- DotPlot(
  clusters_3_4,
  features = markers_ordered,
  group.by = "Cluster_Detailed",
  scale = TRUE
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Cluster 3 & 4 Detailed Marker Expression (Row-scaled)",
    subtitle = "Scaled by row - Easier to compare relative differences",
    x = "Markers",
    y = "Cluster"
  )

ggsave(
  file.path(output_dir, "DotPlot_Cluster3_4_Detailed_Scaled.pdf"),
  p_dotplot_scaled,
  width = 50,
  height = 20,
  limitsize = FALSE
)

# ===== 2.3 按marker类别分面的DotPlot =====

# 创建marker注释
marker_annotation <- data.frame(
  Marker = markers_ordered,
  Category = rep(names(markers_for_plot), sapply(markers_for_plot, length))
)

# 提取DotPlot数据
dotplot_data <- p_dotplot_unscaled$data

# 添加category信息
dotplot_data <- dotplot_data %>%
  left_join(marker_annotation, by = c("features.plot" = "Marker"))

# 按category分面绘图
p_dotplot_faceted <- ggplot(dotplot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_wrap(~Category, scales = "free_y", ncol = 2) +
  scale_color_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 0,
    name = "Scaled\nExpression"
  ) +
  scale_size_continuous(name = "% Expressing") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Cluster 3 & 4: Markers by Category",
    x = "Cluster",
    y = "Markers"
  )

ggsave(
  file.path(output_dir, "DotPlot_Cluster3_4_by_Category.pdf"),
  p_dotplot_faceted,
  width = 16,
  height = 24
)

# ===== 2.4 保存DotPlot数据 =====

# 保存原始数据供后续分析
write.csv(
  dotplot_data,
  file.path(output_dir, "DotPlot_Data_Cluster3_4.csv"),
  row.names = FALSE
)

cat("\n✅ DotPlot generated and saved!\n")


# ===== 3.1 读取用户的annotation.csv =====

cat("\n=== Reading User Annotation File ===\n")
cat("Please make sure your annotation.csv has the following columns:\n")
cat(
  "  - Cluster_ID: matching your cluster names (e.g., 'C3', 'C4_Sub0', 'C4_Sub1', etc.)\n"
)
cat("  - Level1: Broad cell type (e.g., 'Epithelial', 'Stromal', etc.)\n")
cat(
  "  - Level2: Detailed cell type (e.g., 'Club_Progenitor', 'IL33_Myofibroblast', etc.)\n\n"
)

# 示例annotation.csv格式
example_annotation <- data.frame(
  Cluster_ID = c("C3", "C4_Sub0", "C4_Sub1", "C4_Sub2"),
  Level1 = c("Epithelial", "Stromal", "Immune", "Stromal"),
  Level2 = c(
    "Club_Progenitor_Cells",
    "IL33_Myofibroblasts",
    "Antibody_laden_Macrophages",
    "Perivascular_Fibroblasts"
  )
)

# 保存示例文件
write.csv(
  example_annotation,
  file.path(output_dir, "annotation_TEMPLATE.csv"),
  row.names = FALSE
)

cat("Template file saved: annotation_TEMPLATE.csv\n")
cat("Please create your annotation.csv based on this template.\n\n")

# 尝试读取用户的annotation.csv
annotation_file <- "/home/h2048/data/R/1221/serous/Final_Annotation_Results/annotation.csv" # 用户应该把文件放在工作目录

if (file.exists(annotation_file)) {
  cat("✅ Found annotation.csv, reading...\n")
  annotation <- read.csv(annotation_file, stringsAsFactors = FALSE)

  # 验证必需的列
  required_cols <- c("Cluster_ID", "Level1", "Level2")
  missing_cols <- setdiff(required_cols, colnames(annotation))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "ERROR: Missing required columns in annotation.csv: %s\n",
      paste(missing_cols, collapse = ", ")
    ))
  }

  cat("\n=== Annotation File Content ===\n")
  print(annotation)

  # 检查所有cluster是否都有注释
  current_clusters <- unique(clusters_3_4$Cluster_Detailed)
  missing_annotation <- setdiff(current_clusters, annotation$Cluster_ID)

  if (length(missing_annotation) > 0) {
    warning(sprintf(
      "WARNING: The following clusters are missing annotation:\n%s\n",
      paste(missing_annotation, collapse = ", ")
    ))
  }

  extra_annotation <- setdiff(annotation$Cluster_ID, current_clusters)
  if (length(extra_annotation) > 0) {
    warning(sprintf(
      "WARNING: The following Cluster_IDs in annotation.csv don't exist in data:\n%s\n",
      paste(extra_annotation, collapse = ", ")
    ))
  }
} else {
  cat("⚠️ annotation.csv not found. Please create it based on the template.\n")
  cat("After creating annotation.csv, re-run this section.\n\n")

  # 提供一个推荐的注释（基于我们的分析）
  recommended_annotation <- data.frame(
    Cluster_ID = c("C3", "C4_Sub0", "C4_Sub1", "C4_Sub2"),
    Level1 = c("Epithelial", "Stromal", "Immune", "Stromal"),
    Level2 = c(
      "Club_Secretory_Progenitor",
      "Inflammatory_Myofibroblasts_IL33",
      "Macrophages_Antibody_laden",
      "Perivascular_Fibroblasts_FBN1"
    ),
    Notes = c(
      "High BPIFA1, ALDH1A1, GRHL2; progenitor features with goblet differentiation tendency",
      "High VIM, FABP4, IL33, TAGLN; tissue remodeling and type 2 inflammation initiation",
      "High AIF1, S100A4; phagocytosed antibodies (IGKC+ but not plasma cells)",
      "High FBN1, ITGA9, CACNA1C; perivascular support cells"
    )
  )

  write.csv(
    recommended_annotation,
    file.path(output_dir, "annotation_RECOMMENDED.csv"),
    row.names = FALSE
  )

  cat("✅ Recommended annotation saved: annotation_RECOMMENDED.csv\n")
  cat("You can use this as a starting point.\n\n")

  # 使用推荐的注释继续
  annotation <- recommended_annotation %>% select(Cluster_ID, Level1, Level2)
}

# ===== 3.2 应用注释到Seurat对象 =====

if (exists("annotation")) {
  cat("\n=== Applying Annotations ===\n")

  # 创建映射
  cluster_to_level1 <- annotation$Level1
  names(cluster_to_level1) <- annotation$Cluster_ID

  cluster_to_level2 <- annotation$Level2
  names(cluster_to_level2) <- annotation$Cluster_ID

  # 应用到clusters_3_4
  clusters_3_4$Cell_Type_Level1 <- cluster_to_level1[
    clusters_3_4$Cluster_Detailed
  ]
  clusters_3_4$Cell_Type_Level2 <- cluster_to_level2[
    clusters_3_4$Cluster_Detailed
  ]

  # 检查是否有NA（未注释的细胞）
  na_level1 <- sum(is.na(clusters_3_4$Cell_Type_Level1))
  na_level2 <- sum(is.na(clusters_3_4$Cell_Type_Level2))

  if (na_level1 > 0 | na_level2 > 0) {
    warning(sprintf(
      "WARNING: %d cells have NA in Level1, %d cells have NA in Level2\n",
      na_level1,
      na_level2
    ))
  }

  # 统计
  cat("\n=== Annotation Statistics ===\n")
  cat("Level1 Distribution:\n")
  print(table(clusters_3_4$Cell_Type_Level1, useNA = "ifany"))

  cat("\nLevel2 Distribution:\n")
  print(table(clusters_3_4$Cell_Type_Level2, useNA = "ifany"))

  # 也应用到完整的seurat对象
  seurat_obj$Cell_Type_Level1 <- NA
  seurat_obj$Cell_Type_Level2 <- NA

  seurat_obj$Cell_Type_Level1[colnames(
    clusters_3_4
  )] <- clusters_3_4$Cell_Type_Level1
  seurat_obj$Cell_Type_Level2[colnames(
    clusters_3_4
  )] <- clusters_3_4$Cell_Type_Level2

  cat("\n✅ Annotations applied to Seurat object!\n")
}
