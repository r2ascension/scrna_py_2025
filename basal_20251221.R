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

# # ===== Remove Ribosomal, Mitochondrial, ENSG, and Unannotated Genes =====

# # Configuration
# remove_mt <- TRUE
# remove_ribo <- TRUE
# remove_pseudogenes <- TRUE
# remove_ensg <- TRUE
# remove_unannotated <- TRUE

# # Get all gene names
# all_genes <- rownames(seurat_obj)

# # Initialize genes to remove
# genes_to_remove <- c()

# # 1. Mitochondrial genes (MT-)
# if (remove_mt) {
#   mt_genes <- grep('^MT-', all_genes, value = TRUE)
#   genes_to_remove <- c(genes_to_remove, mt_genes)
#   cat(sprintf('Mitochondrial genes: %d\n', length(mt_genes)))
# }

# # 2. Ribosomal genes (RPS, RPL, MRPS, MRPL)
# if (remove_ribo) {
#   ribo_genes <- grep('^RPS|^RPL|^MRPS|^MRPL', all_genes, value = TRUE)
#   genes_to_remove <- c(genes_to_remove, ribo_genes)
#   cat(sprintf('Ribosomal genes: %d\n', length(ribo_genes)))
# }

# # 3. Ribosomal pseudogenes (RPS29P1, RPL10P9, MRPS36P1, etc.)
# if (remove_pseudogenes) {
#   # Match: RPS/RPL/MRPS/MRPL + digits + P + digits
#   pseudo_genes <- grep('^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$', all_genes, value = TRUE)
#   genes_to_remove <- c(genes_to_remove, pseudo_genes)
#   cat(sprintf('Pseudogenes: %d\n', length(pseudo_genes)))
# }

# # 4. ENSG unannotated genes
# if (remove_ensg) {
#   ensg_genes <- grep('^ENSG[0-9]+', all_genes, value = TRUE)
#   genes_to_remove <- c(genes_to_remove, ensg_genes)
#   cat(sprintf('ENSG genes: %d\n', length(ensg_genes)))
# }

# unannotated_genes <- grep(
#   '^(AC|AL|AP|BX|Z)[0-9]+\\.|^RP[0-9]+-|^CTD-|^CTB-|^CTC-|^LINC[0-9]+|-AS[0-9]+$|-OT[0-9]+$|^LOC[0-9]+',
#   all_genes,
#   value = TRUE
# )

# # # 5. Unannotated transcripts
# # if (remove_unannotated) {
# #   unannotated_genes <- grep(
# #     '^(AC|AL|AP|BX|Z)[0-9]+\\.|^RP[0-9]+-|^CT[DBCS]-|^LINC[0-9]+|-AS[0-9]+$|-OT[0-9]+$',
# #     all_genes,
# #     value = TRUE
# #   )
# #   genes_to_remove <- c(genes_to_remove, unannotated_genes)
# #   cat(sprintf('Unannotated transcripts: %d\n', length(unannotated_genes)))
# # }

# genes_to_remove <- unique(c(unannotated_genes,genes_to_remove))

# # Remove duplicates
# genes_to_remove <- unique(genes_to_remove)
# cat(sprintf('\nTotal genes to remove: %d\n', length(genes_to_remove)))

# # Keep genes
# genes_to_keep <- setdiff(all_genes, genes_to_remove)
# cat(sprintf('Genes to keep: %d\n', length(genes_to_keep)))

# # Subset Seurat object
# seurat_obj <- subset(seurat_obj, features = genes_to_keep)

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


# seurat_obj <- readRDS(
#   '/home/h2048/data/R/1221/per_celltype_harmony_rogue/seurat_objects/Basal_harmony.rds'
# )
seurat_obj <- readRDS(
  '/home/h2048/data/R/1221/basal/basal_20251222.rds'
)
output_dir <- '/home/h2048/data/R/1223/basal'
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

pdf("epi_subtype_marker.pdf", width = 48, height = 32, onefile = TRUE)
print(p)
dev.off()


pdf("feature_plots_key_genes_1.pdf", width = 10, height = 8)

for (gene in marker_genes) {
  p <- FeaturePlot(
    seurat_obj,
    features = gene,
    reduction = "umap",
    raster = TRUE,
    pt.size = 0.5
  ) +
    ggtitle(gene) +
    theme(plot.title = element_text(size = 16, face = "bold"))

  cat(gene)
  print(p)
}

dev.off()

pdf(paste0("markers_dotplot.pdf"), width = 64, height = 36)
dot_plot <- DotPlot(
  seurat_obj,
  features = marker_genes,
  group.by = "RNA_snn_res.1",
  cols = c("lightgrey", "red"),
  dot.scale = 8
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(dot_plot)
dev.off()
# 保存点图数据
source_data <- dot_plot$data
write.csv(source_data, paste0("dotplot_data.csv"))
# getwd()
run_tissue_comparison_analysis(
  seurat_obj,
  cell_anno_col = "RNA_snn_res.1",
  tissue_col = "tissue",
  sample_col = "sample",
  min_cell_per_sample = 3,
  min_sample_per_tissue = 3,
  min_pb_libsize = 1000,
  min_pb_detected_genes = 500,
  run_gsva = TRUE,
  run_go = TRUE,
  output_dir = "./analysis_results",
  species = "Homo sapiens",
  padj_thr = 0.05,
  lfc_thr = 1.0,
  min_gene_total_counts = 10,
  top_deg_heatmap = 50,
  gsva_top_var_h = 20,
  gsva_top_var_k = 30,
  n_cores = 1,
  msigdb_gmt = "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt",
  use_lfc_shrink = TRUE
)

# seurat_obj <- readRDS(
#   '/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds'
# )
# seurat_obj <- subset(seurat_obj, subset = Manual_Annotation %in% 'Basal')

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


suppressMessages(library(ROGUE))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
OUTPUT_DIR <- "rogue_results"
GROUPING_COL <- "RNA_snn_res.1" # Change to desired metadata column
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

# seurat_obj <- readRDS(
#   '/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds'
# )
# seurat_obj <- subset(seurat_obj, subset = Manual_Annotation %in% 'Basal')
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
  n.neighbors = 30,
  n.trees = 500,
  min.dist = 0.2,
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
saveRDS(seurat_obj, 'basal_obj_doublet_removed_20251222_2.rds')

# 2. （可选）准备本地GMT
# 下载到固定位置，以后复用
msigdb_gmt <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# 3. 运行（相同调用，结果更准确）
run_one_vs_rest_enrichment(
  seurat_obj,
  group_col = "RNA_snn_res.1",
  msigdb_gmt_file = msigdb_gmt, # ← 新增（可选）
  run_go = TRUE,
  run_gsea = TRUE
)
# ===== Identify, record, and remove doublets =====

# # Step 1: Calculate safe scores
# b_cell_safe <- c('CD79A', 'CD79B', 'MS4A1', 'MZB1', 'JCHAIN', 'IGKC')
# epithelial_safe <- c('KRT5', 'TP63', 'EPCAM', 'CDH1')

# seurat_obj <- AddModuleScore(seurat_obj,
#                              features = list(b_cell_safe),
#                              name = 'B_safe')

# seurat_obj <- AddModuleScore(seurat_obj,
#                              features = list(epithelial_safe),
#                              name = 'Epi_safe')

# # Step 2: Flag doublets
# doublet_threshold_b <- 0.3
# doublet_threshold_epi <- 0.3

# seurat_obj$is_doublet <- (seurat_obj$B_safe1 > doublet_threshold_b) &
#                          (seurat_obj$Epi_safe1 > doublet_threshold_epi)

# seurat_obj$doublet_type <- ifelse(seurat_obj$is_doublet,
#                                   'B-Epithelial_Doublet',
#                                   'Singlet')

# # Step 3: Extract doublet cell information
# doublet_cells <- seurat_obj@meta.data[seurat_obj$is_doublet, ]

# # Add cell barcode as a column
# doublet_cells$cell_barcode <- rownames(doublet_cells)

# # Reorder columns for better readability
# doublet_info <- doublet_cells[, c(
#   'cell_barcode',
#   'orig.ident',
#   'nCount_RNA',
#   'nFeature_RNA',
#   'percent.mt',
#   'RNA_snn_res.1',
#   'B_safe1',
#   'Epi_safe1',
#   'doublet_type'
# )]

# # Rename columns for clarity
# colnames(doublet_info) <- c(
#   'Cell_Barcode',
#   'Sample_ID',
#   'Total_UMI',
#   'Total_Genes',
#   'Percent_MT',
#   'Cluster',
#   'B_cell_Score',
#   'Epithelial_Score',
#   'Doublet_Type'
# )

# # Step 4: Save to CSV
# write.csv(doublet_info,
#           'removed_doublet_cells.csv',
#           row.names = FALSE,
#           quote = FALSE)

# cat(sprintf('Saved %d doublet cells to removed_doublet_cells.csv\n',
#             nrow(doublet_info)))

# # Step 5: Summary statistics by cluster
# doublet_summary <- as.data.frame(table(doublet_cells$RNA_snn_res.1))
# colnames(doublet_summary) <- c('Cluster', 'N_Doublets')
# doublet_summary$Percent <- round(100 * doublet_summary$N_Doublets /
#                                  table(seurat_obj$RNA_snn_res.1)[doublet_summary$Cluster], 2)

# write.csv(doublet_summary,
#           'doublet_summary_by_cluster.csv',
#           row.names = FALSE)

# print(doublet_summary)

# # Step 6: Visualize before removal
# p1 <- DimPlot(seurat_obj, group.by = 'is_doublet',
#               cols = c('grey80', 'red')) +
#       ggtitle(sprintf('Doublets (n=%d)', sum(seurat_obj$is_doublet)))

# p2 <- FeatureScatter(seurat_obj,
#                      feature1 = 'Epi_safe1',
#                      feature2 = 'B_safe1',
#                      group.by = 'is_doublet',
#                      cols = c('grey80', 'red'))

# pdf('doublet_identification.pdf', width = 12, height = 6)
# print(p1 | p2)
# dev.off()

# # Step 7: Remove doublets
# cat(sprintf('\nBefore filtering: %d cells\n', ncol(seurat_obj)))

# seurat_obj <- subset(seurat_obj, subset = is_doublet == FALSE)

# cat(sprintf('After filtering: %d cells\n', ncol(seurat_obj)))
# cat(sprintf('Removed: %d cells (%.2f%%)\n',
#             ncol(seurat_obj) - ncol(seurat_obj),
#             100 * (ncol(seurat_obj) - ncol(seurat_obj)) / ncol(seurat_obj)))

# # Step 8: Save cleaned object

# # ===== Remove doublets AND pure B cells =====
# # seurat_obj <- readRDS('/home/h2048/data/R/1221/basal/basal_20251222.rds')
# # Step 1: Calculate safe scores (if not already done)
# b_cell_safe <- c('CD79A', 'CD79B', 'MS4A1', 'MZB1', 'JCHAIN', 'IGKC')
# epithelial_safe <- c('KRT5', 'TP63', 'EPCAM', 'CDH1')

# seurat_obj <- AddModuleScore(seurat_obj,
#                              features = list(b_cell_safe),
#                              name = 'B_safe')

# seurat_obj <- AddModuleScore(seurat_obj,
#                              features = list(epithelial_safe),
#                              name = 'Epi_safe')

# # Step 2: Classify cells into removal categories
# seurat_obj$removal_reason <- 'Keep'

# # B-Epithelial doublets
# seurat_obj$removal_reason[seurat_obj$B_safe1 > 0.3 & seurat_obj$Epi_safe1 > 0.3] <-
#   'B-Epithelial_Doublet'

# # Pure B cells (contamination)
# seurat_obj$removal_reason[seurat_obj$B_safe1 > 0.7 & seurat_obj$Epi_safe1 < 0.3] <-
#   'B_cell_Contamination'

# # Flag all cells to remove
# seurat_obj$to_remove <- seurat_obj$removal_reason != 'Keep'

# # Step 3: Extract all cells to be removed
# removed_cells <- seurat_obj@meta.data[seurat_obj$to_remove, ]
# removed_cells$cell_barcode <- rownames(removed_cells)

# # Step 4: Prepare export data
# # Check if 'sample' column exists
# if ('sample' %in% colnames(removed_cells)) {
#   sample_col <- 'sample'
# } else if ('Sample' %in% colnames(removed_cells)) {
#   sample_col <- 'Sample'
# } else {
#   sample_col <- 'orig.ident'
#   cat('Warning: No "sample" column found, using orig.ident instead\n')
# }

# removed_export <- removed_cells[, c(
#   'cell_barcode',
#   sample_col,
#   'nCount_RNA',
#   'nFeature_RNA',
#   'percent.mt',
#   'RNA_snn_res.1',
#   'B_safe1',
#   'Epi_safe1',
#   'removal_reason'
# )]

# # Rename columns
# colnames(removed_export) <- c(
#   'Cell_Barcode',
#   'Sample',
#   'Total_UMI',
#   'Total_Genes',
#   'Percent_MT',
#   'Cluster',
#   'B_cell_Score',
#   'Epithelial_Score',
#   'Removal_Reason'
# )

# # Sort by removal reason for easier viewing
# removed_export <- removed_export[order(removed_export$Removal_Reason), ]

# # Step 5: Save to CSV
# write.csv(removed_export,
#           'removed_cells_all.csv',
#           row.names = FALSE,
#           quote = FALSE)

# cat(sprintf('Total cells to remove: %d\n', nrow(removed_export)))
# cat(sprintf('  - B-Epithelial doublets: %d\n',
#             sum(removed_export$Removal_Reason == 'B-Epithelial_Doublet')))
# cat(sprintf('  - B cell contamination: %d\n',
#             sum(removed_export$Removal_Reason == 'B_cell_Contamination')))

# # Step 6: Summary by removal reason and cluster
# removal_summary <- table(removed_cells$RNA_snn_res.1,
#                         removed_cells$removal_reason)

# write.csv(as.data.frame.matrix(removal_summary),
#           'removal_summary_by_cluster.csv',
#           row.names = TRUE)

# print(removal_summary)

# # Step 7: Summary by sample
# removal_by_sample <- table(removed_cells[[sample_col]],
#                           removed_cells$removal_reason)

# write.csv(as.data.frame.matrix(removal_by_sample),
#           'removal_summary_by_sample.csv',
#           row.names = TRUE)

# print(removal_by_sample)

# # Step 8: Visualize
# library(ggplot2)

# # Score scatter plot with removal categories
# plot_data <- FetchData(seurat_obj,
#                        vars = c('B_safe1', 'Epi_safe1', 'removal_reason'))

# p1 <- ggplot(plot_data, aes(x = Epi_safe1, y = B_safe1, color = removal_reason)) +
#   geom_point(alpha = 0.5, size = 0.8) +
#   geom_hline(yintercept = 0.3, linetype = 'dashed', color = 'black') +
#   geom_vline(xintercept = 0.3, linetype = 'dashed', color = 'black') +
#   scale_color_manual(values = c('Keep' = 'grey80',
#                                 'B-Epithelial_Doublet' = 'red',
#                                 'B_cell_Contamination' = 'blue')) +
#   theme_bw() +
#   labs(title = 'Cell Removal Strategy',
#        x = 'Epithelial Score',
#        y = 'B Cell Score',
#        color = 'Cell Type') +
#   annotate('text', x = 0.15, y = 0.6, label = 'Pure B cells\n(Remove)',
#            color = 'blue', size = 3) +
#   annotate('text', x = 0.6, y = 0.6, label = 'Doublets\n(Remove)',
#            color = 'red', size = 3) +
#   annotate('text', x = 0.6, y = 0.15, label = 'Epithelial\n(Keep)',
#            color = 'grey40', size = 3)

# # UMAP visualization
# p2 <- DimPlot(seurat_obj, group.by = 'removal_reason',
#               cols = c('Keep' = 'grey80',
#                       'B-Epithelial_Doublet' = 'red',
#                       'B_cell_Contamination' = 'blue')) +
#       ggtitle(sprintf('Cells to Remove (n=%d)', sum(seurat_obj$to_remove)))

# pdf('cell_removal_visualization.pdf', width = 14, height = 6)
# print(p1 | p2)
# dev.off()

# # Step 9: Remove flagged cells
# cat(sprintf('\nBefore filtering: %d cells\n', ncol(seurat_obj)))

# seurat_obj <- subset(seurat_obj, subset = to_remove == FALSE)

# cat(sprintf('After filtering: %d cells\n', ncol(seurat_obj)))
# cat(sprintf('Removed: %d cells (%.2f%%)\n',
#             ncol(seurat_obj) - ncol(seurat_obj),
#             100 * (ncol(seurat_obj) - ncol(seurat_obj)) / ncol(seurat_obj)))

# # Step 10: Save cleaned object
# # saveRDS(seurat_obj, 'seurat_obj_cleaned_final.rds')

# cat('\n=== Saved files ===\n')
# cat('1. removed_cells_all.csv - All removed cells with reasons\n')
# cat('2. removal_summary_by_cluster.csv - Summary by cluster\n')
# cat('3. removal_summary_by_sample.csv - Summary by sample\n')
# cat('4. cell_removal_visualization.pdf - Visualization\n')
# cat('5. seurat_obj_cleaned_final.rds - Cleaned Seurat object\n')

saveRDS(seurat_obj, 'basal_obj_doublet_removed_20251222.rds')
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

saveRDS(seurat_obj, 'basal_20251222.rds')

# # ===== Merge clusters and recalculate markers =====

# # ===== Fix: Convert factor to character first =====

# # Convert to character to avoid factor level issues
# new_clusters <- as.character(seurat_obj$RNA_snn_res.1)

# # Merge clusters
# new_clusters[new_clusters %in% c('3','11','16')] <- 'c1'
# new_clusters[new_clusters %in% c('17')] <- 'c2'
# new_clusters[new_clusters %in% c('10')] <- 'c3'
# new_clusters[new_clusters %in% c('4')] <- 'c4'
# new_clusters[new_clusters %in% c('1','7','8')] <- 'c5'
# new_clusters[new_clusters %in% c('5','6','12')] <- 'c6'
# new_clusters[new_clusters %in% c('2')] <- 'c7'
# new_clusters[new_clusters %in% c('15')] <- 'c8'
# new_clusters[new_clusters %in% c('13')] <- 'c9'
# new_clusters[new_clusters %in% c('14','9')] <- 'c10'

# seurat_obj$merged_clusters <- factor(new_clusters)
# Idents(seurat_obj) <- 'merged_clusters'
# table(seurat_obj$merged_clusters)
# # Calculate markers
# merged_markers <- FindAllMarkers(
#   seurat_obj,
#   only.pos = TRUE,
#   min.pct = 0.25,
#   logfc.threshold = 0.25,
#   test.use = 'wilcox'
# )

# write.csv(merged_markers, 'merged_clusters_markers.csv', row.names = FALSE)

# # Top 20
# top20 <- merged_markers %>%
#   group_by(cluster) %>%
#   top_n(n = 20, wt = avg_log2FC)

# write.csv(top20, 'merged_clusters_top20_markers.csv', row.names = FALSE)
# # ===== 1.1.1 定义关键marker基因 =====

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
