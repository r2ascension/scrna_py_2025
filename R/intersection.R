# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(rsvd)
library(harmony)
library(presto)
#!/usr/bin/env Rscript
# ============================================================================
# Sample-Gene Availability Matrix Builder
# 样本-基因可用性矩阵构建工具
# ============================================================================
# Purpose: 从原始RDS文件构建样本×基因可用性矩阵
# Output: CSV文件，记录每个样本原本检测到的基因（1=有，0=填充）
# ============================================================================

library(Seurat)
library(Matrix)

# ============================================================================
# Main Function: 构建样本-基因可用性矩阵
# ============================================================================

#' Build sample-gene availability matrix from original RDS files
#' 从原始RDS文件构建样本-基因可用性矩阵
#'
#' @param rds_paths Character vector of RDS file paths
#' @param batch_column Column name in metadata to use as batch ID (default: "orig.ident")
#'                     If NULL or column not found, use filename as batch
#' @param output_file Output CSV filename
#' @param verbose Print progress messages
#' @return Data frame with samples as rows, genes as columns, batch info included
#' 
build_sample_gene_matrix <- function(rds_paths,
                                     batch_column = "orig.ident",
                                     output_file = "sample_gene_availability.csv",
                                     verbose = TRUE) {
  
  if (verbose) {
    cat("\n")
    cat("╔════════════════════════════════════════════════════════╗\n")
    cat("║  Sample-Gene Availability Matrix Builder              ║\n")
    cat("║  样本-基因可用性矩阵构建                                 ║\n")
    cat("╚════════════════════════════════════════════════════════╝\n")
    cat("\n")
    cat("Number of RDS files:", length(rds_paths), "\n")
    cat("Batch column:", ifelse(is.null(batch_column), "Use filename", batch_column), "\n")
    cat("Output file:", output_file, "\n\n")
  }
  
  # ----- Phase 1: 扫描所有RDS，收集样本和基因信息 -----
  if (verbose) {
    cat("========================================\n")
    cat("Phase 1: Scanning RDS files\n")
    cat("========================================\n")
  }
  
  all_samples <- c()
  all_genes <- c()
  rds_info <- list()
  
  for (i in seq_along(rds_paths)) {
    rds_path <- rds_paths[i]
    
    if (verbose) {
      cat(sprintf("[%d/%d] Scanning: %s\n", i, length(rds_paths), basename(rds_path)))
    }
    
    tryCatch({
      # 读取Seurat对象
      seurat_obj <- readRDS(rds_path)
      
      # 提取样本名和基因名
      samples <- colnames(seurat_obj)
      genes <- rownames(seurat_obj)
      
      # 确定批次ID
      if (!is.null(batch_column) && batch_column %in% colnames(seurat_obj@meta.data)) {
        # 使用metadata中的批次信息
        batch_ids <- seurat_obj@meta.data[[batch_column]]
        names(batch_ids) <- samples
      } else {
        # 使用文件名作为批次ID
        batch_id <- tools::file_path_sans_ext(basename(rds_path))
        batch_ids <- setNames(rep(batch_id, length(samples)), samples)
        
        if (verbose && !is.null(batch_column)) {
          cat(sprintf("  ⚠ Column '%s' not found, using filename as batch\n", batch_column))
        }
      }
      
      # 存储信息
      rds_info[[i]] <- list(
        file = rds_path,
        samples = samples,
        genes = genes,
        batch_ids = batch_ids
      )
      
      # 更新全局样本和基因列表
      all_samples <- c(all_samples, samples)
      all_genes <- unique(c(all_genes, genes))
      
      if (verbose) {
        cat(sprintf("  ✓ Samples: %d, Genes: %d, Batches: %s\n", 
                    length(samples), 
                    length(genes),
                    paste(unique(batch_ids), collapse = ", ")))
      }
      
    }, error = function(e) {
      warning(sprintf("Error reading %s: %s", basename(rds_path), e$message))
    })
  }
  
  # 检查样本名唯一性
  if (any(duplicated(all_samples))) {
    warning("Duplicated sample names detected across RDS files!")
    dup_samples <- all_samples[duplicated(all_samples)]
    cat("Duplicated samples:", paste(head(dup_samples, 10), collapse = ", "), "...\n")
  }
  
  if (verbose) {
    cat("\n")
    cat("Scan Summary:\n")
    cat("-----------------------------------\n")
    cat("Total samples:", length(all_samples), "\n")
    cat("Total unique genes (union):", length(all_genes), "\n")
    cat("\n")
  }
  
  # ----- Phase 2: 构建样本×基因矩阵 -----
  if (verbose) {
    cat("========================================\n")
    cat("Phase 2: Building availability matrix\n")
    cat("========================================\n")
  }
  
  # 初始化矩阵（全部为0）
  # 使用稀疏矩阵以节省内存
  availability_matrix <- Matrix(0, 
                                nrow = length(all_samples),
                                ncol = length(all_genes),
                                sparse = TRUE)
  rownames(availability_matrix) <- all_samples
  colnames(availability_matrix) <- all_genes
  
  # 初始化批次信息向量
  sample_batches <- character(length(all_samples))
  names(sample_batches) <- all_samples
  
  # 填充矩阵
  for (i in seq_along(rds_info)) {
    info <- rds_info[[i]]
    
    if (verbose) {
      cat(sprintf("[%d/%d] Processing: %s\n", 
                  i, length(rds_info), basename(info$file)))
    }
    
    # 为这个RDS中的样本标记其检测到的基因
    for (sample in info$samples) {
      # 标记该样本有的基因为1
      availability_matrix[sample, info$genes] <- 1
      # 记录批次信息
      sample_batches[sample] <- info$batch_ids[sample]
    }
  }
  
  if (verbose) {
    cat("\n")
    cat("Matrix Statistics:\n")
    cat("-----------------------------------\n")
    cat("Matrix dimensions:", nrow(availability_matrix), "×", ncol(availability_matrix), "\n")
    cat("Total elements:", length(availability_matrix), "\n")
    cat("Non-zero elements:", sum(availability_matrix), 
        sprintf("(%.2f%%)\n", sum(availability_matrix)/length(availability_matrix)*100))
    cat("Average genes per sample:", round(rowSums(availability_matrix) %>% mean(), 0), "\n")
    cat("Memory usage:", format(object.size(availability_matrix), units = "MB"), "\n")
    cat("\n")
  }
  
  # ----- Phase 3: 转换为data frame并保存 -----
  if (verbose) {
    cat("========================================\n")
    cat("Phase 3: Converting to CSV format\n")
    cat("========================================\n")
  }
  
  # 转换稀疏矩阵为常规矩阵（CSV需要）
  if (verbose) cat("Converting sparse matrix to dense format...\n")
  availability_df <- as.matrix(availability_matrix)
  
  # 添加批次信息列
  if (verbose) cat("Adding batch information...\n")
  availability_df <- data.frame(
    batch = sample_batches[rownames(availability_df)],
    availability_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # 保存为CSV
  if (verbose) cat(sprintf("Writing to file: %s\n", output_file))
  write.csv(availability_df, 
            output_file, 
            row.names = TRUE,
            quote = FALSE)
  
  if (verbose) {
    cat("\n")
    cat("╔════════════════════════════════════════════════════════╗\n")
    cat("║  Complete!                                             ║\n")
    cat("╚════════════════════════════════════════════════════════╝\n")
    cat("\n")
    cat("Output file created:", output_file, "\n")
    cat("File size:", format(file.size(output_file)/1024^2, digits=2), "MB\n")
    cat("\n")
    cat("CSV structure:\n")
    cat("  - Row names: Sample IDs\n")
    cat("  - First column: 'batch' (batch information)\n")
    cat("  - Remaining columns: Gene names\n")
    cat("  - Values: 1 (gene detected) / 0 (zero-padded)\n")
    cat("\n")
    cat("Quick check:\n")
    cat("  Head of first 5 samples and 5 genes:\n")
    print(availability_df[1:min(5, nrow(availability_df)), 
                          1:min(6, ncol(availability_df))])
    cat("\n")
  }
  
  return(invisible(availability_df))
}


# ============================================================================
# 辅助函数：快速诊断和统计
# ============================================================================

#' Quick diagnosis of availability matrix
#' 快速诊断可用性矩阵
#'
#' @param csv_file Path to the availability CSV file
#' @param output_dir Directory for diagnostic plots
#' 
diagnose_availability_matrix <- function(csv_file = "sample_gene_availability.csv",
                                         output_dir = "diagnostics") {
  
  cat("\n")
  cat("========================================\n")
  cat("Loading and Diagnosing Availability Matrix\n")
  cat("========================================\n\n")
  
  # 读取CSV
  cat("Reading CSV file...\n")
  availability_df <- read.csv(csv_file, row.names = 1, check.names = FALSE)
  
  # 分离批次信息和基因矩阵
  batch_info <- availability_df$batch
  gene_matrix <- as.matrix(availability_df[, -1])
  
  cat("Matrix dimensions:", nrow(gene_matrix), "samples ×", ncol(gene_matrix), "genes\n\n")
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----- 1. 样本层面统计 -----
  cat("1. Sample-level Statistics\n")
  cat("-----------------------------------\n")
  
  genes_per_sample <- rowSums(gene_matrix)
  
  cat("Average genes per sample:", round(mean(genes_per_sample), 0), "\n")
  cat("Range:", min(genes_per_sample), "-", max(genes_per_sample), "\n")
  cat("Median:", median(genes_per_sample), "\n\n")
  
  # 可视化
  pdf(file.path(output_dir, "genes_per_sample_distribution.pdf"), width = 10, height = 6)
  hist(genes_per_sample, 
       breaks = 30,
       col = "#3498db",
       border = "white",
       main = "Distribution of Gene Counts per Sample",
       xlab = "Number of Genes Detected",
       ylab = "Number of Samples")
  abline(v = median(genes_per_sample), col = "red", lwd = 2, lty = 2)
  legend("topright", 
         legend = paste("Median:", median(genes_per_sample)),
         col = "red", lty = 2, lwd = 2)
  dev.off()
  
  # ----- 2. 基因层面统计 -----
  cat("2. Gene-level Statistics\n")
  cat("-----------------------------------\n")
  
  samples_per_gene <- colSums(gene_matrix)
  
  cat("Average samples per gene:", round(mean(samples_per_gene), 1), "\n")
  cat("Genes in all samples:", sum(samples_per_gene == nrow(gene_matrix)), "\n")
  cat("Genes in <10% samples:", sum(samples_per_gene < nrow(gene_matrix)*0.1), "\n\n")
  
  # 可视化
  pdf(file.path(output_dir, "samples_per_gene_distribution.pdf"), width = 10, height = 6)
  hist(samples_per_gene, 
       breaks = 50,
       col = "#e74c3c",
       border = "white",
       main = "Distribution of Sample Coverage per Gene",
       xlab = "Number of Samples Gene Detected In",
       ylab = "Number of Genes")
  abline(v = nrow(gene_matrix)*0.5, col = "blue", lwd = 2, lty = 2)
  legend("topright", 
         legend = paste("50% threshold:", nrow(gene_matrix)*0.5),
         col = "blue", lty = 2, lwd = 2)
  dev.off()
  
  # ----- 3. 批次统计 -----
  cat("3. Batch-level Statistics\n")
  cat("-----------------------------------\n")
  
  batch_summary <- data.frame(
    batch = names(table(batch_info)),
    n_samples = as.vector(table(batch_info)),
    mean_genes = tapply(genes_per_sample, batch_info, mean) %>% round(0),
    stringsAsFactors = FALSE
  )
  
  print(batch_summary)
  cat("\n")
  
  # 保存统计表
  write.csv(batch_summary, 
            file.path(output_dir, "batch_summary.csv"),
            row.names = FALSE)
  
  cat("✓ Diagnostic reports saved to:", output_dir, "\n\n")
  
  return(invisible(list(
    genes_per_sample = genes_per_sample,
    samples_per_gene = samples_per_gene,
    batch_summary = batch_summary
  )))
}


# ============================================================================
# 使用示例 (Usage Examples)
# ============================================================================

# # ----- 基础使用 -----
# 
# # 1. 准备RDS文件路径列表
# rds_files <- c(
#   "data/batch1_samples.rds",
#   "data/batch2_samples.rds",
#   "data/batch3_samples.rds"
# )
# 
# 或者自动扫描文件夹
rds_files <- list.files("E:/R/Source/final/",
                       pattern = "\\.rds$",
                       full.names = TRUE)
# 
# # 2. 构建可用性矩阵
# # 默认使用metadata中的"orig.ident"列作为批次
# build_sample_gene_matrix(
#   rds_paths = rds_files,
#   batch_column = "orig.ident",  # 如果有的话会使用，没有就用文件名
#   output_file = "sample_gene_availability.csv",
#   verbose = TRUE
# )
# 
# # 3. 可选：运行诊断分析
# diagnose_availability_matrix(
#   csv_file = "sample_gene_availability.csv",
#   output_dir = "availability_diagnostics"
# )
# 
# # ----- 读取和使用矩阵 -----
# 
# # 读取生成的CSV文件
# availability <- read.csv("sample_gene_availability.csv", 
#                         row.names = 1, 
#                         check.names = FALSE)
# 
# # 查看结构
# head(availability[, 1:10])
# 
# # 提取批次信息
# sample_batches <- availability$batch
# 
# # 提取基因矩阵（移除batch列）
# gene_matrix <- as.matrix(availability[, -1])
# 
# # 示例：查找在所有样本中都检测到的基因
# universal_genes <- colnames(gene_matrix)[colSums(gene_matrix) == nrow(gene_matrix)]
# cat("Genes in all samples:", length(universal_genes), "\n")
# 
# # 示例：查找某个样本检测到的基因
# sample_id <- rownames(gene_matrix)[1]
# genes_in_sample <- colnames(gene_matrix)[gene_matrix[sample_id, ] == 1]
# cat("Genes in", sample_id, ":", length(genes_in_sample), "\n")
# 
# # 示例：为差异分析过滤基因
# # 保留在至少50%样本中检测到的基因
# min_sample_threshold <- 0.5
# genes_to_keep <- colnames(gene_matrix)[
#   colSums(gene_matrix) >= nrow(gene_matrix) * min_sample_threshold
# ]
# cat("Genes passing 50% threshold:", length(genes_to_keep), "\n")
# 
# # ----- 特殊场景：不同的批次列名 -----
# 
# # 如果您的Seurat对象中批次信息在"sample"或"batch"列
build_sample_gene_matrix(
  rds_paths = rds_files,
  batch_column = "dataset",  # 或 "batch", "library", 等
  output_file = "sample_gene_availability.csv"
)
# 
# # 如果完全不想使用metadata，只用文件名
# build_sample_gene_matrix(
#   rds_paths = rds_files,
#   batch_column = NULL,  # 强制使用文件名
#   output_file = "sample_gene_availability.csv"
# )


