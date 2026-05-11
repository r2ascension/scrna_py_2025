# ============================================================
# 单细胞RNA测序数据整合分析 - 修复版
# Single-cell RNA-seq Data Integration - Fixed Version
# ============================================================

library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)


setwd("E:/R/1010/")

# ============================================================
# Step 1: 定义数据源配置
# Data Source Configuration
# ============================================================

cat("=== 开始数据整合流程 ===\n\n")

data_sources <- list(
  # Seon_Pyo_Hong_2023 数据集
  list(
    file = "E:/R/Source/20251010/GSE207083_RAW/GSM6276895_female_2_normalized_expression_matrix.csv.gz",
    sample_id = "GSM6276895",
    study = "Seon_Pyo_Hong_2023",
    condition = "Healthy",
    batch = "Seon_Pyo_Hong_2023"
  ),
  list(
    file = "E:/R/Source/20251010/GSE207083_RAW/GSM6276892_male_1_normalized_expression_matrix.csv.gz",
    sample_id = "GSM6276892",
    study = "Seon_Pyo_Hong_2023",
    condition = "Healthy",
    batch = "Seon_Pyo_Hong_2023"
  ),
  list(
    file = "E:/R/Source/20251010/GSE207083_RAW/GSM6276893_male_2_normalized_expression_matrix.csv.gz",
    sample_id = "GSM6276893",
    study = "Seon_Pyo_Hong_2023",
    condition = "Healthy",
    batch = "Seon_Pyo_Hong_2023"
  ),
  list(
    file = "E:/R/Source/20251010/GSE207083_RAW/GSM6276894_female_1_normalized_expression_matrix.csv.gz",
    sample_id = "GSM6276894",
    study = "Seon_Pyo_Hong_2023",
    condition = "Healthy",
    batch = "Seon_Pyo_Hong_2023"
  )
)
# seurat_obj$sample[seurat_obj$sample == "Nasal_Mucosa"] <- 'GSM4695772'
# seurat_obj$sample[seurat_obj$sample == "Healthy_2"] <- 'GSM5226282'
# seurat_obj$sample[seurat_obj$sample == "Healthy_1"] <- 'GSM5226281'
# =# ============================================================
# 使用 fread() 读取所有文件
# ============================================================

matrices_list <- list()

for (i in seq_along(data_sources)) {
  source_info <- data_sources[[i]]
  cat(sprintf(
    "\n【文件 %d/%d】: %s\n",
    i,
    length(data_sources),
    source_info$sample_id
  ))
  cat(sprintf("路径: %s\n", source_info$file))

  tryCatch(
    {
      # 使用 fread() 自动检测分隔符
      mat <- fread(
        source_info$file,
        header = TRUE,
        sep = "auto", # 自动检测是逗号还是制表符
        stringsAsFactors = FALSE,
        data.table = FALSE # 返回data.frame而非data.table
      )

      cat(sprintf("  ✓ 读取成功！维度: %d 行 × %d 列\n", nrow(mat), ncol(mat)))

      # 第一列是基因名
      genes <- mat[[1]]
      expr <- mat[, -1] # 移除第一列

      cat(sprintf("  处理后: %d 基因 × %d 细胞\n", length(genes), ncol(expr)))

      # 检查数据类型
      if (!is.numeric(expr[[1]])) {
        cat("  转换为数值型...\n")
        expr <- as.data.frame(lapply(expr, as.numeric))
      }

      # 验证数值范围
      sample_values <- as.numeric(unlist(expr[
        1:min(100, nrow(expr)),
        1:min(10, ncol(expr))
      ]))
      sample_values <- sample_values[!is.na(sample_values)]

      cat(sprintf(
        "  数值范围: %.2f - %.2f\n",
        min(sample_values),
        max(sample_values)
      ))
      cat(sprintf(
        "  非零值比例: %.2f%%\n",
        sum(sample_values > 0) / length(sample_values) * 100
      ))

      # 设置行名和列名
      rownames(expr) <- genes

      # 为细胞ID添加样本前缀
      original_cell_ids <- colnames(expr)
      new_cell_ids <- paste0(source_info$sample_id, "_", original_cell_ids)
      colnames(expr) <- new_cell_ids

      # 保存到列表
      matrices_list[[i]] <- list(
        genes = genes,
        expr = expr,
        source_info = source_info,
        cell_ids = new_cell_ids,
        original_barcodes = original_cell_ids
      )

      cat("  ✓ 验证通过，已保存\n")
    },
    error = function(e) {
      cat(sprintf("  ❌ 读取失败: %s\n", e$message))
    }
  )
}

# ============================================================
# 检查读取结果
# ============================================================

success_count <- sum(!sapply(matrices_list, is.null))
cat(sprintf(
  "\n=== 成功读取 %d/%d 个文件 ===\n\n",
  success_count,
  length(data_sources)
))

if (success_count == 0) {
  stop("❌ 所有文件读取失败！")
}

# 移除失败的读取
matrices_list <- matrices_list[!sapply(matrices_list, is.null)]

# 显示每个样本的信息
cat("各样本统计:\n")
for (i in seq_along(matrices_list)) {
  mat_data <- matrices_list[[i]]
  cat(sprintf(
    "  %s: %d 基因 × %d 细胞\n",
    mat_data$source_info$sample_id,
    length(mat_data$genes),
    ncol(mat_data$expr)
  ))
}

# ============================================================
# 获取所有基因的并集
# ============================================================

cat("\n标准化基因集...\n")

all_genes <- Reduce(union, lapply(matrices_list, function(x) x$genes))
cat(sprintf("  总基因数(并集): %d\n", length(all_genes)))

# 检查重复基因
dup_genes <- duplicated(all_genes)
if (any(dup_genes)) {
  cat(sprintf("  警告: 发现 %d 个重复基因名\n", sum(dup_genes)))
  all_genes <- unique(all_genes)
}

# ============================================================
# 对齐并合并矩阵
# ============================================================

cat("\n对齐表达矩阵...\n")

aligned_matrices <- lapply(seq_along(matrices_list), function(i) {
  cat(sprintf(
    "  处理样本 %d/%d: %s\n",
    i,
    length(matrices_list),
    matrices_list[[i]]$source_info$sample_id
  ))

  current_genes <- matrices_list[[i]]$genes
  current_expr <- matrices_list[[i]]$expr

  # 创建对齐矩阵
  aligned_matrix <- matrix(
    0,
    nrow = length(all_genes),
    ncol = ncol(current_expr)
  )

  # 找到基因位置
  gene_positions <- match(current_genes, all_genes)

  # 填充数据
  aligned_matrix[gene_positions, ] <- as.matrix(current_expr)

  # 设置列名
  colnames(aligned_matrix) <- colnames(current_expr)

  # 验证
  non_zero <- sum(aligned_matrix > 0)
  cat(sprintf(
    "    非零值数量: %d (%.2f%%)\n",
    non_zero,
    non_zero / length(aligned_matrix) * 100
  ))

  return(aligned_matrix)
})

# ============================================================
# 合并所有矩阵
# ============================================================

cat("\n合并表达矩阵...\n")

merged_expr <- do.call(cbind, aligned_matrices)
rownames(merged_expr) <- all_genes

cat(sprintf(
  "  合并后维度: %d 基因 × %d 细胞\n",
  nrow(merged_expr),
  ncol(merged_expr)
))

# ============================================================
# 创建 Metadata
# ============================================================

cat("\n构建细胞元数据...\n")

metadata_list <- lapply(matrices_list, function(mat_data) {
  n_cells <- length(mat_data$cell_ids)

  metadata_df <- data.frame(
    cell_id = mat_data$cell_ids,
    original_barcode = mat_data$original_barcodes,
    sample_id = rep(mat_data$source_info$sample_id, n_cells),
    study = rep(mat_data$source_info$study, n_cells),
    condition = rep(mat_data$source_info$condition, n_cells),
    batch = rep(mat_data$source_info$batch, n_cells),
    row.names = mat_data$cell_ids,
    stringsAsFactors = FALSE
  )

  return(metadata_df)
})

combined_metadata <- do.call(rbind, metadata_list)

cat("  样本分布:\n")
print(table(combined_metadata$sample_id))
cat("\n  疾病状态分布:\n")
print(table(combined_metadata$condition))

# 验证一致性
if (!all(colnames(merged_expr) == rownames(combined_metadata))) {
  stop("❌ 表达矩阵和metadata顺序不一致!")
}

# ============================================================
# 转换为稀疏矩阵并保存
# ============================================================

cat("\n转换为稀疏矩阵...\n")

sparse_matrix <- Matrix(merged_expr, sparse = TRUE)
rownames(sparse_matrix) <- all_genes
colnames(sparse_matrix) <- colnames(merged_expr)

# 计算稀疏度
total_elements <- prod(dim(sparse_matrix))
non_zero_elements <- length(sparse_matrix@x)
sparsity <- (1 - non_zero_elements / total_elements) * 100

cat(sprintf("  稀疏度: %.2f%% (非零元素: %d)\n", sparsity, non_zero_elements))

# 保存文件
cat("\n保存中间文件...\n")

write.table(
  all_genes,
  "genes_fixed.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
write.table(
  colnames(sparse_matrix),
  "cells_fixed.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
writeMM(sparse_matrix, "expression_matrix_fixed.mtx")
write.csv(combined_metadata, "metadata_fixed.csv", row.names = TRUE)

cat("  已保存:\n")
cat("    - expression_matrix_fixed.mtx\n")
cat("    - genes_fixed.txt\n")
cat("    - cells_fixed.txt\n")
cat("    - metadata_fixed.csv\n")

# ============================================================
# 清理内存
# ============================================================

rm(merged_expr, aligned_matrices)
gc()

# ============================================================
# 创建 Seurat 对象
# ============================================================

cat("\n创建 Seurat 对象...\n")

seurat_obj <- CreateSeuratObject(
  counts = sparse_matrix,
  meta.data = combined_metadata,
  project = "Nasal_Fixed",
  min.cells = 3,
  min.features = 200
)

# 添加质控指标
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^RP[SL]"
)

# ============================================================
# 查看结果
# ============================================================

cat("\n=== Seurat 对象创建完成 ===\n\n")

cat("对象信息:\n")
print(seurat_obj)

cat("\n质控指标统计:\n")
summary_stats <- summary(seurat_obj@meta.data[, c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt"
)])
print(summary_stats)

cat("\n每个样本的细胞数:\n")
print(table(seurat_obj$sample_id))

cat("\n疾病状态分布:\n")
print(table(seurat_obj$condition))

# ============================================================
# 质控可视化
# ============================================================

cat("\n生成质控图...\n")

library(ggplot2)
library(patchwork)

# 按样本分组
p1 <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0.1
)

ggsave("QC_by_sample_fixed.pdf", p1, width = 15, height = 5)

# 按条件分组
p2 <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "condition",
  ncol = 3,
  pt.size = 0.1
)

ggsave("QC_by_condition_fixed.pdf", p2, width = 12, height = 5)

cat("  已保存质控图:\n")
cat("    - QC_by_sample_fixed.pdf\n")
cat("    - QC_by_condition_fixed.pdf\n")

# ============================================================
# 保存 Seurat 对象
# ============================================================

seurat_obj@meta.data$dataset <- 'Seon_Pyo_Hong_2023'
seurat_obj@meta.data$study <- 'Seon_Pyo_Hong_2023'
seurat_obj@meta.data$tissue <- 'Nose'
seurat_obj@meta.data$sample <- seurat_obj$sample_id
seurat_obj@meta.data$sample_id <- NULL
seurat_obj@meta.data$tissue_sampling_method <- 'biopsy'
seurat_obj@meta.data$For_MASC <- FALSE
# seurat_obj$tissue_sampling_method[seurat_obj$sample_id == 'Nasal_Mucosa'] <- 'biopsy'
saveRDS(seurat_obj, "Seon_Pyo_Hong_2023.rds")

cat("\n✓ 所有处理完成！\n")
cat("✓ Seurat 对象已保存: seurat_object_fixed.rds\n\n")

# ============================================================
# 生成分析报告
# ============================================================

cat("=== 数据整合报告 ===\n\n")

report <- data.frame(
  Item = c(
    "总基因数",
    "总细胞数(过滤前)",
    "总细胞数(过滤后)",
    "样本数量",
    "疾病样本细胞数",
    "健康对照细胞数",
    "平均基因数/细胞",
    "平均UMI数/细胞",
    "平均线粒体比例(%)"
  ),
  Value = c(
    length(all_genes),
    ncol(sparse_matrix),
    ncol(seurat_obj),
    length(unique(combined_metadata$sample_id)),
    sum(seurat_obj$condition == "Disease"),
    sum(seurat_obj$condition == "Healthy"),
    round(mean(seurat_obj$nFeature_RNA), 1),
    round(mean(seurat_obj$nCount_RNA), 1),
    round(mean(seurat_obj$percent.mt), 2)
  )
)

print(report, row.names = FALSE)

write.csv(report, "integration_report_fixed.csv", row.names = FALSE)

cat("\n建议的下一步分析:\n")
cat("1. 质量过滤\n")
cat("2. 标准化 (NormalizeData)\n")
cat("3. 高变基因筛选 (FindVariableFeatures)\n")
cat("4. 批次校正 (Harmony)\n")
cat("5. 降维聚类 (PCA, UMAP, FindClusters)\n")
cat("6. 细胞类型注释\n\n")
