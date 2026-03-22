# ============================================================================================
# 单细胞RNA测序数据处理流程 - 整合版
# Single-cell RNA-seq Analysis Pipeline - Integrated Version
# 
# 功能特性 Features:
# 1. 批次数据读取与质控 Batch data loading and QC
# 2. EmptyDrops空液滴过滤 EmptyDrops filtering
# 3. 回归异常值检测 Regression-based outlier detection
# 4. Checkpoint双细胞检测 Checkpoint-based doublet detection
# 5. Assay5兼容性处理 Assay5 compatibility handling
# 6. 内存优化策略 Memory optimization strategies
# ============================================================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

# ============================================================================================
# 加载必要的包 Load Required Packages
# ============================================================================================
library(DropletUtils)
library(scater)
library(Seurat)
library(DoubletFinder)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
Sys.setenv(RETICULATE_PYTHON = "C:/Users/崔填祎/AppData/Local/Programs/Python/Python311/python.exe")
library(leiden)
library(harmony)
library(reticulate)
library(kBET)
library(pryr) 
library(cluster)
library(presto)

# 设置随机种子和并行计算 Set random seed and parallel computing
set.seed(42)
options(future.globals.maxSize = Inf)

# 设置工作目录 Set working directory
dir <- "E:/R/1010"
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
setwd(dir)

# 创建输出目录 Create output directories
output_dirs <- c(
  "batch_results",     # 批次处理结果 Batch processing results
  "analysis_output",   # 最终分析结果 Final analysis results
  "qaqc_plots",        # QC图 QC plots
  "doublet_results"    # 双细胞检测结果 Doublet detection results
)

for(output_dir in output_dirs) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
}

# ============================================================================================
# 函数定义 Function Definitions
# ============================================================================================

# ============================================================================================
# 函数1: Assay5到Assay转换 (Assay5 to Assay Conversion)
# ============================================================================================
rebuild_assay_from_scratch <- function(seurat_obj, assay = "RNA") {
  #' 从Assay5完全重建v3/v4兼容的Assay
  #' Completely rebuild v3/v4 compatible Assay from Assay5
  #' 
  #' @param seurat_obj Seurat对象 Seurat object
  #' @param assay assay名称 Assay name
  #' @return 重建后的Seurat对象 Reconstructed Seurat object
  
  cat("=== Starting Assay Reconstruction ===\n")
  cat("=== 开始重建Assay ===\n\n")
  
  assay_obj <- seurat_obj[[assay]]
  
  if (!inherits(assay_obj, "Assay5")) {
    cat("Already in Assay format, no reconstruction needed\n")
    cat("当前已是Assay格式，无需重建\n")
    return(seurat_obj)
  }
  
  # Step 1: Extract counts data
  cat("Step 1/3: Extracting counts data...\n")
  cat("步骤1/3: 提取counts数据...\n")
  
  tryCatch({
    counts_data <- LayerData(seurat_obj, assay = assay, layer = "counts")
    cat(sprintf("  Counts matrix: %d genes × %d cells\n", 
                nrow(counts_data), ncol(counts_data)))
  }, error = function(e) {
    stop("Unable to extract counts data: ", e$message)
  })
  
  # Step 2: Extract normalized data (if exists)
  cat("Step 2/3: Extracting normalized data...\n")
  cat("步骤2/3: 提取normalized data...\n")
  
  tryCatch({
    data_data <- LayerData(seurat_obj, assay = assay, layer = "data")
    has_data <- TRUE
    cat(sprintf("  Data matrix: %d genes × %d cells\n", 
                nrow(data_data), ncol(data_data)))
  }, error = function(e) {
    cat("  Data layer not found, will recalculate from counts\n")
    cat("  未找到data层，将从counts重新计算\n")
    data_data <- NULL
    has_data <- FALSE
  })
  
  # Step 3: Create new Assay object
  cat("Step 3/3: Creating new Assay object...\n")
  cat("步骤3/3: 创建新的Assay对象...\n")
  
  new_assay <- CreateAssayObject(counts = counts_data)
  
  if (has_data) {
    new_assay <- SetAssayData(new_assay, slot = "data", new.data = data_data)
    cat("  ✓ Normalized data added\n")
    cat("  ✓ 已添加normalized data\n")
  } else {
    cat("  Performing normalization...\n")
    cat("  正在进行标准化...\n")
    temp_obj <- CreateSeuratObject(counts = counts_data)
    temp_obj <- NormalizeData(temp_obj, verbose = FALSE)
    new_assay <- SetAssayData(new_assay, 
                              slot = "data",
                              new.data = GetAssayData(temp_obj, slot = "data"))
    rm(temp_obj)
    cat("  ✓ Normalization completed\n")
    cat("  ✓ 已完成标准化\n")
  }
  
  # Replace original assay
  seurat_obj[[assay]] <- new_assay
  
  cat("\n✓ Assay reconstruction completed!\n")
  cat("✓ Assay重建完成!\n\n")
  
  # Validation
  cat("Validation | 验证结果:\n")
  cat(sprintf("  - Assay type | Assay类型: %s\n", class(seurat_obj[[assay]])[1]))
  cat(sprintf("  - Counts: %d genes × %d cells\n", 
              nrow(GetAssayData(seurat_obj, slot = "counts")),
              ncol(GetAssayData(seurat_obj, slot = "counts"))))
  cat(sprintf("  - Data: %d genes × %d cells\n", 
              nrow(GetAssayData(seurat_obj, slot = "data")),
              ncol(GetAssayData(seurat_obj, slot = "data"))))
  
  return(seurat_obj)
}

# ============================================================================================
# 函数2: EmptyDrops处理 (EmptyDrops Processing)
# ============================================================================================
process_emptydrops <- function(sce, lower = 100, fdr_threshold = 0.01) {
  #' EmptyDrops空液滴过滤
  #' EmptyDrops empty droplet filtering
  #' 
  #' @param sce Seurat对象或counts矩阵 Seurat object or counts matrix
  #' @param lower 最低UMI阈值 Minimum UMI threshold
  #' @param fdr_threshold FDR阈值 FDR threshold
  #' @return 细胞索引 Cell indices
  
  if (inherits(sce, "Seurat")) {
    mat <- GetAssayData(sce, slot = "counts")
  } else {
    mat <- sce
  }
  
  set.seed(42)
  e.out <- DropletUtils::emptyDrops(
    m = mat,
    lower = lower,
    niters = 10000
  )
  
  is.cell <- e.out$FDR <= fdr_threshold
  is.cell[is.na(is.cell)] <- FALSE
  
  n_cells <- sum(is.cell)
  cat(sprintf("EmptyDrops identified %d cells (FDR < %.3f)\n", 
              n_cells, fdr_threshold))
  
  return(which(is.cell))
}

# ============================================================================================
# 函数3: 回归异常值检测 (Regression-based Outlier Detection)
# ============================================================================================
regression_outliers <- function(seurat_obj, outliers_threshold = 0.999) {
  #' 基于回归模型的异常细胞检测
  #' Regression-based outlier cell detection
  #' 
  #' @param seurat_obj Seurat对象 Seurat object
  #' @param outliers_threshold 异常值置信水平 Outlier confidence level
  #' @return 保留的细胞名称 Names of cells to keep
  
  nCount <- seurat_obj$nCount_RNA
  nFeature <- seurat_obj$nFeature_RNA
  
  log_counts <- log10(nCount + 1)
  log_features <- log10(nFeature + 1)
  
  fit <- lm(log_features ~ log_counts)
  
  pred <- predict(fit, 
                  data.frame(log_counts = log_counts), 
                  interval = "prediction", 
                  level = outliers_threshold)
  
  outliers <- log_features < pred[,"lwr"] | log_features > pred[,"upr"]
  
  # Visualization
  pdf(file.path("qaqc_plots", 
                paste0("regression_outliers_", 
                       unique(seurat_obj$sample), ".pdf")),
      width = 10, height = 8)
  
  plot(log_counts, log_features,
       pch = 16, cex = 0.6,
       col = ifelse(outliers, "red", "grey50"),
       xlab = "log10(nCount_RNA)",
       ylab = "log10(nFeature_RNA)",
       main = "Regression Outliers Detection")
  
  lines(sort(log_counts), pred[order(log_counts),"lwr"], 
        col = "blue", lty = 2)
  lines(sort(log_counts), pred[order(log_counts),"upr"], 
        col = "blue", lty = 2)
  lines(sort(log_counts), pred[order(log_counts),"fit"], 
        col = "blue", lwd = 2)
  
  legend("topleft",
         c("Cells", "Outliers", "Fit", "Prediction Interval"),
         col = c("grey50", "red", "blue", "blue"),
         pch = c(16, 16, NA, NA),
         lty = c(NA, NA, 1, 2),
         lwd = c(NA, NA, 2, 1))
  
  dev.off()
  
  cat(sprintf("Detected %d outliers (%.1f%%)\n", 
              sum(outliers), mean(outliers) * 100))
  
  cells.keep <- colnames(seurat_obj)[!outliers]
  return(cells.keep)
}

# ============================================================================================
# 函数4: 生成QC报告 (Generate QC Report)
# ============================================================================================
generate_qc_report <- function(seurat_obj, stage, output_dir = "qaqc_plots") {
  #' 生成质控报告
  #' Generate quality control report
  #' 
  #' @param seurat_obj Seurat对象 Seurat object
  #' @param stage 分析阶段名称 Analysis stage name
  #' @param output_dir 输出目录 Output directory
  #' @return 统计数据框 Statistics data frame
  
  stats <- list(
    "Number of Cells" = ncol(seurat_obj),
    "Number of Genes" = nrow(seurat_obj),
    "Median Features per Cell" = median(seurat_obj$nFeature_RNA),
    "Median Counts per Cell" = median(seurat_obj$nCount_RNA),
    "Median MT%" = median(seurat_obj$percent.mt)
  )
  
  if("Phase" %in% colnames(seurat_obj@meta.data)) {
    phase_counts <- table(seurat_obj$Phase)
    stats <- c(stats, list(
      "Cells in G1" = as.numeric(phase_counts["G1"]),
      "Cells in G2M" = as.numeric(phase_counts["G2M"]),
      "Cells in S" = as.numeric(phase_counts["S"])
    ))
  } else {
    stats <- c(stats, list(
      "Cells in G1" = NA,
      "Cells in G2M" = NA,
      "Cells in S" = NA
    ))
  }
  
  stats[["Number of Variable Genes"]] <- length(VariableFeatures(seurat_obj))
  
  stats_df <- data.frame(
    Metric = names(stats),
    Value = unlist(stats)
  )
  
  write.csv(stats_df,
            file = file.path(output_dir, paste0("qc_stats_", stage, ".csv")),
            row.names = FALSE)
  
  # Generate QC plots
  pdf(file.path(output_dir, paste0("qc_report_", stage, ".pdf")),
      width = 15, height = 12)
  
  print(VlnPlot(seurat_obj,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3,
                pt.size = 0.1) + 
          ggtitle(paste("QC Metrics -", stage)))
  
  p1 <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
  p2 <- FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt", pt.size = 0.5)
  print(p1 + p2)
  
  if("Phase" %in% colnames(seurat_obj@meta.data)) {
    print(VlnPlot(seurat_obj,
                  features = c("S.Score", "G2M.Score"),
                  ncol = 2,
                  pt.size = 0.1) + 
            ggtitle("Cell Cycle Scores"))
    
    phase_props <- prop.table(table(seurat_obj$Phase))
    barplot(phase_props,
            main = "Cell Cycle Phase Distribution",
            ylab = "Proportion",
            col = c("lightblue", "lightgreen", "pink"))
  }
  
  dev.off()
  
  return(stats_df)
}

# ============================================================================================
# 函数5: 生成基础QC图 (Generate Basic QC Plots)
# ============================================================================================
generate_qc_plots <- function(seurat_obj, sample_name) {
  #' 生成基础QC可视化图表
  #' Generate basic QC visualization plots
  
  p1 <- VlnPlot(seurat_obj, 
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3)
  
  p2 <- FeatureScatter(seurat_obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "nFeature_RNA")
  
  p3 <- FeatureScatter(seurat_obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "percent.mt")
  
  pdf(file.path("qaqc_plots", paste0(sample_name, "_qc.pdf")),
      width = 15, height = 10)
  print(p1)
  print(p2 + p3)
  dev.off()
}

# ============================================================================================
# 函数6: 读取10x数据 (Read 10x Data)
# 注意：移除了旧的DoubletFinder调用，双细胞检测将在所有批次合并后统一进行
# ============================================================================================
read_10x_data <- function(
    data_dir, 
    sample_name,
    tissue_name,
    batch,
    mt_threshold = 20,
    use_emptydrops = FALSE,
    do_regression = TRUE,
    study = NULL,
    dataset = NULL,
    tissue_sampling_method = NULL
) {
  #' 读取10x数据并进行基础处理（不包含双细胞检测）
  #' Read 10x data and perform basic processing (without doublet detection)
  #' 
  #' @note Cell cycle scoring is skipped to avoid JoinLayers operations
  
  cat(sprintf("\n>>> Reading %s (batch: %s) <<<\n", sample_name, batch))
  
  # 1. Read data
  tryCatch({
    sce <- Read10X(data_dir)
    
    # Handle multiple data types (Gene Expression + Antibody Capture, etc.)
    if(is.list(sce)) {
      cat("  - Multiple data types detected:\n")
      for(type_name in names(sce)) {
        cat(sprintf("    * %s: %d features\n", type_name, nrow(sce[[type_name]])))
      }
      
      # Extract Gene Expression data
      if("Gene Expression" %in% names(sce)) {
        cat("  - Using 'Gene Expression' data\n")
        sce <- sce[["Gene Expression"]]
      } else if("RNA" %in% names(sce)) {
        cat("  - Using 'RNA' data\n")
        sce <- sce[["RNA"]]
      } else {
        # Use the first element if standard names not found
        cat(sprintf("  - Using '%s' data (first element)\n", names(sce)[1]))
        sce <- sce[[1]]
      }
    }
    
    if(ncol(sce) == 0) {
      stop("No cells found in the data")
    }
    
    cat("  - Raw cells:", ncol(sce), "| Genes:", nrow(sce), "\n")
    
  }, error = function(e) {
    stop(e)
  })
  
  # 2. EmptyDrops filtering
  if(use_emptydrops) {
    cat("Running EmptyDrops...\n")
    cells.keep <- process_emptydrops(sce)
    sce <- sce[, cells.keep]
    cat(sprintf("  - Kept %d cells after EmptyDrops\n", ncol(sce)))
  }
  
  # 3. Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = sce,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # 4. Add metadata
  seurat_obj$sample <- sample_name
  seurat_obj$tissue <- tissue_name
  seurat_obj$batch <- batch
  
  # Add optional metadata
  if(!is.null(study)) {
    seurat_obj$study <- study
  }
  if(!is.null(dataset)) {
    seurat_obj$dataset <- dataset
  }
  if(!is.null(tissue_sampling_method)) {
    seurat_obj$tissue_sampling_method <- tissue_sampling_method
  }
  
  cat(sprintf("  - Processed: %d cells, %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))
  
  # 5. Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rp"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # 6. Filter based on MT%
  cells_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = percent.mt <= mt_threshold)
  cells_after <- ncol(seurat_obj)
  cat("  - Filtered out", cells_before - cells_after, "cells based on percent.mt\n")
  
  # 7. Regression analysis
  if(do_regression) {
    cat("Running regression analysis...\n") 
    cells.keep <- regression_outliers(seurat_obj)
    seurat_obj <- subset(seurat_obj, cells = cells.keep)
  }
  
  # 8. Cell cycle scoring - SKIPPED
  # Note: Cell cycle scoring is skipped during batch processing to avoid
  # JoinLayers operations. It can be performed after batch integration if needed.
  cat("Cell cycle scoring skipped (will be done after integration if needed)\n")
  
  # 9. Generate QC report
  qc_stats <- generate_qc_report(seurat_obj, 
                                 stage = paste0("initial_", sample_name))
  
  return(seurat_obj)
}

# ============================================================================================
# 函数7: 处理单个批次 (Process Single Batch)
# ============================================================================================
process_batch <- function(
    sample_list, 
    tissue_name, 
    batch,
    mt_threshold = 20,
    use_emptydrops = FALSE,
    study = NULL,
    dataset = NULL,
    tissue_sampling_method = NULL
) {
  #' 处理单个批次的所有样本
  #' Process all samples in a single batch
  
  cat(sprintf("\n=== Processing batch: %s (%s) ===\n", batch, tissue_name))
  if(!is.null(study)) {
    cat(sprintf("    Study: %s\n", study))
  }
  if(!is.null(dataset)) {
    cat(sprintf("    Dataset: %s\n", dataset))
  }
  if(!is.null(tissue_sampling_method)) {
    cat(sprintf("    Sampling method: %s\n", tissue_sampling_method))
  }
  
  tryCatch({
    current_mem <- pryr::mem_used()
    cat(sprintf("Memory usage: %.2f GB\n", current_mem / 1024^3))
  }, error = function(e) {
    cat("Note: Memory tracking unavailable\n")
  })
  
  batch_objs <- lapply(sample_list, function(sample) {
    cat(sprintf("\n=== Processing sample: %s ===\n", sample$name))
    
    sobj <- read_10x_data(
      data_dir = sample$path,
      sample_name = sample$name,
      tissue_name = tissue_name,
      batch = batch,
      mt_threshold = mt_threshold,
      use_emptydrops = use_emptydrops,
      do_regression = TRUE,
      study = study,
      dataset = dataset,
      tissue_sampling_method = tissue_sampling_method
    )
    
    generate_qc_plots(sobj, sample$name)
    
    return(sobj)
  })
  
  batch_objs <- batch_objs[!sapply(batch_objs, is.null)]
  if(length(batch_objs) == 0) {
    stop("No samples were successfully processed in this batch")
  }
  
  cat("\nMerging batch samples...\n")
  if(length(batch_objs) > 1) {
    batch_objs <- lapply(batch_objs, function(x) {
      if(length(Layers(x[["RNA"]])) > 1) {
        x[["RNA"]] <- JoinLayers(x[["RNA"]])
      }
      return(x)
    })
    
    merged_batch <- merge(
      x = batch_objs[[1]], 
      y = batch_objs[-1],
      add.cell.ids = names(batch_objs),
      merge.data = TRUE
    )
  } else {
    merged_batch <- batch_objs[[1]]
  }
  
  if(length(Layers(merged_batch[["RNA"]])) > 1) {
    merged_batch[["RNA"]] <- JoinLayers(merged_batch[["RNA"]])
  }
  
  cat("Splitting RNA assay into layers by batch...\n")
  merged_batch[["RNA"]] <- split(merged_batch[["RNA"]], 
                                 f = merged_batch$batch)
  
  cat("Normalizing data...\n")
  merged_batch <- NormalizeData(merged_batch)
  merged_batch <- FindVariableFeatures(merged_batch)
  
  save_name <- file.path("batch_results", 
                         paste0(tissue_name, "_", batch, ".rds"))
  saveRDS(merged_batch, file = save_name)
  cat(sprintf("  - Batch saved to: %s\n", save_name))
  
  qc_stats <- generate_qc_report(merged_batch,
                                 stage = paste0("batch_", batch))
  
  write.csv(qc_stats,
            file = file.path("batch_results", 
                             paste0("qc_stats_", batch, ".csv")),
            row.names = FALSE)
  
  rm(batch_objs)
  gc()
  
  return(merged_batch)
}

# ============================================================================================
# 函数8: 样本级QC预检查 (Sample-level QC Pre-inspection)
# ============================================================================================
generate_sample_qc_inspection <- function(seurat_obj_main, 
                                          sample_column = "sample",
                                          output_dir = "qaqc_plots") {
  #' 为每个样本生成独立的质控图和统计报告
  #' Generate individual QC plots and statistics for each sample
  #' 
  #' @param seurat_obj_main 主Seurat对象 Main Seurat object
  #' @param sample_column 样本列名 Sample column name
  #' @param output_dir 输出目录 Output directory
  #' @return 样本统计数据框 Sample statistics dataframe
  
  cat("\n============================================================\n")
  cat("   Sample-level QC Inspection\n")
  cat("   样本级质控预检查\n")
  cat("============================================================\n\n")
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Get all sample names
  sample_names <- unique(seurat_obj_main@meta.data[[sample_column]])
  n_samples <- length(sample_names)
  
  cat(sprintf("Total samples to inspect: %d\n\n", n_samples))
  
  # Set identity
  Idents(seurat_obj_main) <- sample_column
  
  # Initialize statistics dataframe
  all_stats <- data.frame()
  
  # Process each sample
  for(i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    
    cat(sprintf("[%d/%d] Processing sample: %s\n", i, n_samples, sample_name))
    
    # Extract cells for this sample
    sample_cells <- WhichCells(seurat_obj_main, 
                               expression = get(sample_column) == sample_name)
    
    if(length(sample_cells) == 0) {
      cat(sprintf("  ⚠️  Warning: No cells found for sample: %s\n\n", sample_name))
      next
    }
    
    # Create sample subset
    sample_obj <- subset(seurat_obj_main, cells = sample_cells)
    n_cells <- ncol(sample_obj)
    
    cat(sprintf("  Cells: %d\n", n_cells))
    
    # Generate QC plots
    tryCatch({
      # Violin plots for QC metrics
      p1 <- VlnPlot(sample_obj, 
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, 
                    pt.size = 0.1) +
        plot_annotation(title = paste("QC Metrics -", sample_name),
                        theme = theme(plot.title = element_text(hjust = 0.5, 
                                                                face = "bold", 
                                                                size = 14)))
      
      # Scatter plot - RNA count vs feature count
      p2 <- FeatureScatter(sample_obj, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA",
                           pt.size = 0.8) +
        ggtitle("UMI Count vs Gene Count")
      
      # Scatter plot - RNA count vs MT%
      p3 <- FeatureScatter(sample_obj, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt",
                           pt.size = 0.8) +
        ggtitle("UMI Count vs MT%")
      
      # Save plots to PDF
      pdf(file.path(output_dir, paste0(sample_name, "_qc_inspection.pdf")),
          width = 15, height = 10)
      print(p1)
      print(p2 + p3)
      dev.off()
      
      cat(sprintf("  ✓ QC plots saved to: %s\n", 
                  paste0(sample_name, "_qc_inspection.pdf")))
      
    }, error = function(e) {
      cat(sprintf("  ✗ Error generating plots: %s\n", e$message))
    })
    
    # Calculate summary statistics
    stats <- data.frame(
      Sample = sample_name,
      Cells = n_cells,
      Median_Genes = median(sample_obj$nFeature_RNA),
      Mean_Genes = round(mean(sample_obj$nFeature_RNA), 1),
      SD_Genes = round(sd(sample_obj$nFeature_RNA), 1),
      Median_UMIs = median(sample_obj$nCount_RNA),
      Mean_UMIs = round(mean(sample_obj$nCount_RNA), 1),
      SD_UMIs = round(sd(sample_obj$nCount_RNA), 1),
      Median_Mt_Percent = round(median(sample_obj$percent.mt), 2),
      Mean_Mt_Percent = round(mean(sample_obj$percent.mt), 2),
      SD_Mt_Percent = round(sd(sample_obj$percent.mt), 2),
      stringsAsFactors = FALSE
    )
    
    # Print statistics
    cat(sprintf("  Median genes/cell: %d | Median UMIs/cell: %d | Median MT%%: %.2f%%\n",
                stats$Median_Genes, stats$Median_UMIs, stats$Median_Mt_Percent))
    
    # Append to all_stats
    all_stats <- rbind(all_stats, stats)
    
    cat("\n")
    
    # Clean up
    rm(sample_obj)
    gc(verbose = FALSE)
  }
  
  # Save summary statistics
  csv_file <- file.path(output_dir, "sample_qc_metrics.csv")
  write.csv(all_stats, csv_file, row.names = FALSE)
  
  cat("============================================================\n")
  cat("✓ Sample QC Inspection Completed\n")
  cat("✓ 样本质控预检查完成\n")
  cat("============================================================\n\n")
  
  cat("Output files:\n")
  cat(sprintf("  📊 %s (summary statistics)\n", csv_file))
  cat(sprintf("  📈 %d × PDF files in %s/\n", n_samples, output_dir))
  
  # Print summary table
  cat("\nSummary Statistics:\n")
  print(all_stats)
  cat("\n")
  
  # Generate a summary comparison plot
  if(nrow(all_stats) > 1) {
    tryCatch({
      # Prepare data for plotting
      plot_data <- all_stats %>%
        select(Sample, Mean_Genes, Mean_UMIs, Mean_Mt_Percent) %>%
        tidyr::pivot_longer(cols = -Sample, 
                            names_to = "Metric", 
                            values_to = "Value")
      
      # Create faceted comparison plot
      p_compare <- ggplot(plot_data, aes(x = Sample, y = Value, fill = Sample)) +
        geom_bar(stat = "identity", color = "black") +
        facet_wrap(~Metric, scales = "free_y", ncol = 1) +
        theme_classic(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              strip.background = element_rect(fill = "grey90"),
              strip.text = element_text(face = "bold")) +
        labs(title = "QC Metrics Comparison Across Samples",
             x = "Sample", y = "Value") +
        scale_fill_brewer(palette = "Set3")
      
      ggsave(file.path(output_dir, "samples_qc_comparison.pdf"),
             p_compare, 
             width = max(10, n_samples * 0.5), 
             height = 12)
      
      cat(sprintf("  📈 samples_qc_comparison.pdf (comparison plot)\n"))
      
    }, error = function(e) {
      cat(sprintf("  Note: Could not generate comparison plot: %s\n", e$message))
    })
  }
  
  cat("\n")
  
  return(all_stats)
}

# ============================================================================================
# 函数9: Checkpoint双细胞检测 (Checkpoint-based Doublet Detection)
# ============================================================================================
run_doublet_detection_checkpoint <- function(
    seurat_obj_main,
    doublet_rate = 0.06,
    sample_column = "sample",
    output_dir = "doublet_results"
) {
  #' 使用checkpoint机制的双细胞检测流程
  #' Doublet detection pipeline with checkpoint mechanism
  #' 
  #' @param seurat_obj_main 主Seurat对象 Main Seurat object
  #' @param doublet_rate 预期双细胞率 Expected doublet rate
  #' @param sample_column 样本列名 Sample column name
  #' @param output_dir 输出目录 Output directory
  #' @return 更新后的Seurat对象 Updated Seurat object
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("\n============================================================\n")
  cat("   Doublet Detection Pipeline - Checkpoint Mode\n")
  cat("   双细胞检测流程 - Checkpoint模式\n")
  cat("============================================================\n\n")
  
  # Step 0: Clean old DoubletFinder columns
  cat("Step 0: Cleaning old DoubletFinder columns...\n")
  old_df_cols <- grep("DF.classifications|pANN", 
                      colnames(seurat_obj_main@meta.data), 
                      value = TRUE)
  
  if(length(old_df_cols) > 0) {
    cat(sprintf("  Found %d old columns, removing...\n", length(old_df_cols)))
    seurat_obj_main@meta.data[, old_df_cols] <- NULL
  }
  
  # Step 1: Initialize checkpoint matrix
  cat("\nStep 1: Initializing checkpoint matrix...\n")
  all_cell_ids <- colnames(seurat_obj_main)
  n_total_cells <- length(all_cell_ids)
  
  checkpoint_df <- data.frame(
    cell_id = all_cell_ids,
    sample = seurat_obj_main@meta.data[[sample_column]],
    doublet_status = NA_character_,
    doublet_score = NA_real_,
    doublet_class = NA_character_,
    processing_time = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(checkpoint_df) <- all_cell_ids
  
  cat(sprintf("  Created checkpoint matrix: %d cells × %d columns\n", 
              nrow(checkpoint_df), ncol(checkpoint_df)))
  
  seurat_obj_main$doublet_status <- NA_character_
  seurat_obj_main$doublet_score <- NA_real_
  seurat_obj_main$doublet_class <- NA_character_
  
  # Step 2: Get sample list
  if(!sample_column %in% colnames(seurat_obj_main@meta.data)) {
    stop(sprintf("Error: Column '%s' not found!", sample_column))
  }
  
  samples <- unique(seurat_obj_main@meta.data[[sample_column]])
  n_samples <- length(samples)
  
  cat(sprintf("\nTotal samples: %d\n", n_samples))
  cat(sprintf("Total cells: %d\n\n", n_total_cells))
  
  stats_list <- list()
  
  # Step 3: Process each sample
  for(i in 1:n_samples) {
    sample_name <- samples[i]
    
    cat("\n============================================================\n")
    cat(sprintf("[%d/%d] Processing Sample: %s\n", i, n_samples, sample_name))
    cat("============================================================\n\n")
    
    start_time <- Sys.time()
    
    cell_idx <- seurat_obj_main@meta.data[[sample_column]] == sample_name
    cell_names <- colnames(seurat_obj_main)[cell_idx]
    n_cells <- length(cell_names)
    
    cat(sprintf("Cells in sample: %d\n", n_cells))
    
    # Check cell count
    if(n_cells < 50) {
      cat("⚠️  WARNING: Insufficient cells (<50), skipping...\n\n")
      
      seurat_obj_main$doublet_status[cell_idx] <- "Insufficient_cells"
      seurat_obj_main$doublet_score[cell_idx] <- NA
      seurat_obj_main$doublet_class[cell_idx] <- "Skipped"
      
      checkpoint_df[cell_names, "doublet_status"] <- "Insufficient_cells"
      checkpoint_df[cell_names, "doublet_score"] <- NA
      checkpoint_df[cell_names, "doublet_class"] <- "Skipped"
      checkpoint_df[cell_names, "processing_time"] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      
      stats_list[[sample_name]] <- data.frame(
        Sample = sample_name,
        Total = n_cells,
        Singlets = 0,
        Doublets = 0,
        Rate = 0,
        Status = "Skipped"
      )
      
      saveRDS(checkpoint_df, file.path(output_dir, "checkpoint_matrix.rds"))
      next
    }
    
    # Extract sample subset
    cat("Extracting sample subset...\n")
    obj_sample <- subset(seurat_obj_main, cells = cell_names)
    
    old_cols_in_subset <- grep("DF.classifications|pANN", 
                               colnames(obj_sample@meta.data), 
                               value = TRUE)
    if(length(old_cols_in_subset) > 0) {
      obj_sample@meta.data[, old_cols_in_subset] <- NULL
    }
    
    nExp <- round(doublet_rate * n_cells)
    cat(sprintf("Expected doublets: %d (%.1f%%)\n\n", nExp, doublet_rate * 100))
    
    # Run DoubletFinder
    tryCatch({
      cat("  [1/4] Normalizing data...\n")
      if(is.null(obj_sample@assays$RNA@data) || 
         nrow(obj_sample@assays$RNA@data) == 0) {
        obj_sample <- NormalizeData(obj_sample, verbose = FALSE)
      }
      
      obj_sample <- FindVariableFeatures(obj_sample, 
                                         selection.method = "vst", 
                                         nfeatures = 2000, 
                                         verbose = FALSE)
      obj_sample <- ScaleData(obj_sample, 
                              features = VariableFeatures(obj_sample), 
                              verbose = FALSE)
      
      cat("  [2/4] Running PCA...\n")
      n_pcs <- min(30, ncol(obj_sample) - 1)
      n_pcs <- max(n_pcs, 10)
      obj_sample <- RunPCA(obj_sample, npcs = n_pcs, verbose = FALSE)
      
      cat("  [3/4] Optimizing pK parameter...\n")
      sweep.res <- paramSweep(obj_sample, PCs = 1:n_pcs, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      
      pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      if(is.na(pK) || pK < 0.01 || pK > 0.3) {
        pK <- 0.09
      }
      cat(sprintf("    Using pK: %.3f\n", pK))
      
      cat("  [4/4] Running DoubletFinder...\n")
      obj_sample <- doubletFinder(obj_sample, 
                                  PCs = 1:n_pcs, 
                                  pN = 0.25, 
                                  pK = pK, 
                                  nExp = nExp, 
                                  sct = FALSE)
      
      # Extract results
      df_class_col <- grep("DF.classifications", 
                           colnames(obj_sample@meta.data), 
                           value = TRUE)
      df_score_col <- grep("pANN", 
                           colnames(obj_sample@meta.data), 
                           value = TRUE)
      
      if(length(df_class_col) > 0) {
        df_class_col <- df_class_col[length(df_class_col)]
        df_score_col <- if(length(df_score_col) > 0) {
          df_score_col[length(df_score_col)]
        } else {
          NULL
        }
      }
      
      if(length(df_class_col) > 0) {
        cat("\n  ✓ Extracting results...\n")
        
        status_values <- obj_sample@meta.data[[df_class_col]]
        score_values <- if(!is.null(df_score_col)) {
          obj_sample@meta.data[[df_score_col]]
        } else {
          rep(NA, n_cells)
        }
        
        # Update main object and checkpoint
        seurat_obj_main$doublet_status[cell_idx] <- status_values
        seurat_obj_main$doublet_score[cell_idx] <- score_values
        seurat_obj_main$doublet_class[cell_idx] <- ifelse(
          status_values == "Singlet", "Singlet", "Doublet"
        )
        
        checkpoint_df[cell_names, "doublet_status"] <- status_values
        checkpoint_df[cell_names, "doublet_score"] <- score_values
        checkpoint_df[cell_names, "doublet_class"] <- ifelse(
          status_values == "Singlet", "Singlet", "Doublet"
        )
        checkpoint_df[cell_names, "processing_time"] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        
        n_singlets <- sum(status_values == "Singlet", na.rm = TRUE)
        n_doublets <- sum(status_values == "Doublet", na.rm = TRUE)
        
        end_time <- Sys.time()
        elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
        
        cat(sprintf("\n  Results:\n"))
        cat(sprintf("    - Singlets: %d (%.1f%%)\n", 
                    n_singlets, n_singlets/n_cells*100))
        cat(sprintf("    - Doublets: %d (%.1f%%)\n", 
                    n_doublets, n_doublets/n_cells*100))
        cat(sprintf("    - Time: %.1f sec\n", elapsed_time))
        
        stats_list[[sample_name]] <- data.frame(
          Sample = sample_name,
          Total = n_cells,
          Singlets = n_singlets,
          Doublets = n_doublets,
          Rate = n_doublets/n_cells*100,
          Status = "Success",
          Time_sec = round(elapsed_time, 1)
        )
        
        # Save checkpoint
        cat(sprintf("\n  💾 Saving checkpoint...\n"))
        saveRDS(checkpoint_df, file.path(output_dir, "checkpoint_matrix.rds"))
        write.csv(checkpoint_df, 
                  file.path(output_dir, "checkpoint_matrix.csv"), 
                  row.names = TRUE)
      }
      
    }, error = function(e) {
      cat(sprintf("\n  ✗ ERROR: %s\n", e$message))
      
      seurat_obj_main$doublet_status[cell_idx] <- "Error"
      seurat_obj_main$doublet_score[cell_idx] <- NA
      seurat_obj_main$doublet_class[cell_idx] <- "Error"
      
      checkpoint_df[cell_names, "doublet_status"] <- "Error"
      checkpoint_df[cell_names, "doublet_score"] <- NA
      checkpoint_df[cell_names, "doublet_class"] <- "Error"
      checkpoint_df[cell_names, "processing_time"] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      
      stats_list[[sample_name]] <- data.frame(
        Sample = sample_name,
        Total = n_cells,
        Singlets = 0,
        Doublets = 0,
        Rate = 0,
        Status = "Error",
        Time_sec = NA
      )
      
      saveRDS(checkpoint_df, file.path(output_dir, "checkpoint_matrix.rds"))
    })
    
    rm(obj_sample)
    gc(verbose = FALSE)
  }
  
  # Step 4: Generate summary and visualizations
  cat("\n\n============================================================\n")
  cat("FINAL SUMMARY | 最终汇总\n")
  cat("============================================================\n\n")
  
  stats <- do.call(rbind, stats_list)
  print(stats)
  
  total_cells <- ncol(seurat_obj_main)
  total_singlets <- sum(seurat_obj_main$doublet_class == "Singlet", na.rm = TRUE)
  total_doublets <- sum(seurat_obj_main$doublet_class == "Doublet", na.rm = TRUE)
  
  cat(sprintf("\nOverall: %d cells | %d singlets (%.1f%%) | %d doublets (%.1f%%)\n",
              total_cells, 
              total_singlets, total_singlets/total_cells*100,
              total_doublets, total_doublets/total_cells*100))
  
  # Save results
  write.csv(stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  saveRDS(seurat_obj_main, file.path(output_dir, "seurat_with_doublets_final.rds"))
  
  # Generate visualizations
  p1 <- ggplot(stats, aes(x = Sample, y = Rate, fill = Status)) +
    geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = doublet_rate * 100, 
               linetype = "dashed", color = "red", size = 1) +
    geom_text(aes(label = sprintf("%.1f%%", Rate)), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Success" = "#3498DB", "Skipped" = "#95A5A6", 
                                 "Failed" = "#E74C3C", "Error" = "#E67E22")) +
    labs(title = "Doublet Detection Rate by Sample", 
         x = "Sample", y = "Doublet Rate (%)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "doublet_rate_barplot.pdf"), 
         p1, width = max(10, n_samples * 0.4), height = 6)
  
  cat("\n✓ Doublet detection completed!\n")
  cat("✓ 双细胞检测完成!\n\n")
  
  return(seurat_obj_main)
}

# ============================================================================================
# 主执行流程 MAIN EXECUTION PIPELINE
# ============================================================================================

cat("\n")
cat("============================================================================================\n")
cat("  Single-cell RNA-seq Analysis Pipeline - Sinus Dataset\n")
cat("  单细胞RNA测序分析流程 - 窦样本数据集\n")
cat("============================================================================================\n")
cat("\n")
cat("Dataset Information | 数据集信息:\n")
cat("  - Source: GSE235715\n")
cat("  - Study: Guanrui Liao 2025\n")
cat("  - Tissue: Sinus (副鼻窦)\n")
cat("  - Sampling method: Biopsy (活检)\n")
cat("  - Number of samples: 5\n")
cat("  - Sample IDs: GSM7508289-GSM7508293\n")
cat("============================================================================================\n\n")

# Define sample lists (定义样本列表)

# Sinus samples (窦样本 - Guanrui Liao 2025 dataset)
# Source: GSE235715
# Study/Dataset/Batch: Guanrui_Liao_2025
# Tissue: sinus
# Sampling method: biopsy
Sinus_samples <- list(
  list(path = "E:/R/Source/20251010/GSE235715_RAW/GSM7508289_Control_1_filtered", 
       name = "GSM7508289"),
  list(path = "E:/R/Source/20251010/GSE235715_RAW/GSM7508290_Control_2_filtered", 
       name = "GSM7508290"),
  list(path = "E:/R/Source/20251010/GSE235715_RAW/GSM7508291_Control_3_filtered", 
       name = "GSM7508291"),
  list(path = "E:/R/Source/20251010/GSE235715_RAW/GSM7508292_Control_4_filtered", 
       name = "GSM7508292"),
  list(path = "E:/R/Source/20251010/GSE235715_RAW/GSM7508293_Control_5_filtered", 
       name = "GSM7508293")
)

# ============================================================================================
# PHASE 1: Batch Processing (批次处理阶段)
# 注意：此阶段不进行双细胞检测
# ============================================================================================

cat("\n")
cat("============================================================================================\n")
cat("PHASE 1: Batch Processing\n")
cat("阶段1: 批次处理\n")
cat("============================================================================================\n\n")

# Process Sinus samples
cat("\nProcessing Sinus samples...\n")
sinus <- process_batch(Sinus_samples, 
                       tissue_name = "sinus", 
                       batch = "Guanrui_Liao_2025", 
                       mt_threshold = 20,
                       study = "Guanrui_Liao_2025",
                       dataset = "Guanrui_Liao_2025",
                       tissue_sampling_method = "biopsy")
rm(sinus); gc()

# ============================================================================================
# PHASE 2: Load Processed Batch (加载处理后的批次)
# ============================================================================================

cat("\n")
cat("============================================================================================\n")
cat("PHASE 2: Load Processed Batch\n")
cat("阶段2: 加载处理后的批次\n")
cat("============================================================================================\n\n")

set.seed(42)

# Load the processed batch file
batch_file <- "batch_results/sinus_Guanrui_Liao_2025.rds"

if(!file.exists(batch_file)) {
  stop("Batch file not found: ", batch_file, "\nPlease ensure PHASE 1 completed successfully.")
}

cat("Loading processed batch file...\n")
merged_obj <- readRDS(batch_file)

cat("\nLoaded data statistics:\n")
cat(sprintf("- Total cells: %d\n", ncol(merged_obj)))
cat(sprintf("- Total features: %d\n", nrow(merged_obj)))
cat(sprintf("- Number of samples: %d\n", length(unique(merged_obj$sample))))

# Display sample distribution
cat("\nSample distribution:\n")
print(table(merged_obj$sample))

# Display metadata summary
cat("\nMetadata summary:\n")
cat(sprintf("- Tissue: %s\n", unique(merged_obj$tissue)))
cat(sprintf("- Batch: %s\n", unique(merged_obj$batch)))
cat(sprintf("- Study: %s\n", unique(merged_obj$study)))
cat(sprintf("- Dataset: %s\n", unique(merged_obj$dataset)))
cat(sprintf("- Sampling method: %s\n", unique(merged_obj$tissue_sampling_method)))

# Save merged object before doublet detection
saveRDS(merged_obj, "analysis_output/merged_before_doublets.rds")
cat("\n✓ Saved to: analysis_output/merged_before_doublets.rds\n")

# ============================================================================================
# PHASE 3: Sample QC Inspection & Doublet Detection (质控检查与双细胞检测阶段)
# ============================================================================================

cat("\n")
cat("============================================================================================\n")
cat("PHASE 3: Sample QC Inspection & Doublet Detection\n")
cat("阶段3: 样本质控检查与双细胞检测\n")
cat("============================================================================================\n\n")

# Step 1: Sample-level QC inspection (样本级质控预检查)
cat("Step 1: Sample-level QC inspection before doublet detection...\n")
cat("步骤1: 双细胞检测前的样本质控预检查...\n\n")

sample_qc_stats <- generate_sample_qc_inspection(
  seurat_obj_main = merged_obj,
  sample_column = "sample",
  output_dir = "qaqc_plots"
)

cat("\n⚠️  CHECKPOINT: Please review the QC plots in 'qaqc_plots/' directory.\n")
cat("⚠️  检查点: 请检查'qaqc_plots/'目录中的质控图。\n")
cat("    - Review sample_qc_metrics.csv for overall statistics\n")
cat("    - Check individual sample PDF files for detailed QC plots\n")
cat("    - Verify samples_qc_comparison.pdf for cross-sample comparison\n\n")

cat("Press [Enter] to continue with doublet detection, or [Ctrl+C] to stop and review...\n")
cat("按 [Enter] 继续双细胞检测，或按 [Ctrl+C] 停止并检查...\n")
readline(prompt = "")

# Step 2: Convert Assay5 to Assay format (关键步骤!)
cat("\nStep 2: Converting Assay5 to Assay format for DoubletFinder compatibility...\n")
cat("步骤2: 转换Assay5格式以兼容DoubletFinder...\n\n")

merged_obj <- rebuild_assay_from_scratch(merged_obj, assay = "RNA")

# Step 3: Run doublet detection
cat("\nStep 3: Running doublet detection with checkpoint...\n")
cat("步骤3: 运行checkpoint双细胞检测...\n\n")

merged_obj <- run_doublet_detection_checkpoint(
  seurat_obj_main = merged_obj,
  doublet_rate = 0.06,
  sample_column = "sample",
  output_dir = "doublet_results"
)

# ============================================================================================
# PHASE 4: Final Processing (最终处理阶段)
# ============================================================================================

cat("\n")
cat("============================================================================================\n")
cat("PHASE 4: Final Processing\n")
cat("阶段4: 最终处理\n")
cat("============================================================================================\n\n")

# Filter out doublets (optional)
cat("Filtering doublets (optional, you can skip this step)...\n")
merged_obj_clean <- subset(merged_obj, subset = doublet_class == "Singlet")

cat(sprintf("\nCells before filtering: %d\n", ncol(merged_obj)))
cat(sprintf("Cells after filtering: %d\n", ncol(merged_obj_clean)))
cat(sprintf("Doublets removed: %d (%.1f%%)\n", 
            ncol(merged_obj) - ncol(merged_obj_clean),
            (ncol(merged_obj) - ncol(merged_obj_clean))/ncol(merged_obj)*100))

# Save final objects
saveRDS(merged_obj, "Guanrui_Liao_2025.rds")
saveRDS(merged_obj_clean, "analysis_output/merged_singlets_only.rds")

# ============================================================================================
# Pipeline Completed!
# ============================================================================================

cat("\n")
cat("============================================================================================\n")
cat("✓ PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("✓ 流程成功完成!\n")
cat("============================================================================================\n\n")

cat("Output files:\n")
cat("输出文件:\n\n")
cat("📁 batch_results/\n")
cat("   - sinus_Guanrui_Liao_2025.rds (processed batch)\n")
cat("   - qc_stats_Guanrui_Liao_2025.csv (batch QC statistics)\n\n")
cat("📁 qaqc_plots/\n")
cat("   - sample_qc_metrics.csv (summary statistics for all 5 samples)\n")
cat("   - samples_qc_comparison.pdf (cross-sample comparison plot)\n")
cat("   - GSM7508289_qc_inspection.pdf (sample 1 detailed QC)\n")
cat("   - GSM7508290_qc_inspection.pdf (sample 2 detailed QC)\n")
cat("   - GSM7508291_qc_inspection.pdf (sample 3 detailed QC)\n")
cat("   - GSM7508292_qc_inspection.pdf (sample 4 detailed QC)\n")
cat("   - GSM7508293_qc_inspection.pdf (sample 5 detailed QC)\n")
cat("   - Regression outlier plots for each sample\n\n")
cat("📁 doublet_results/\n")
cat("   - checkpoint_matrix.rds (lightweight checkpoint)\n")
cat("   - checkpoint_matrix.csv (human-readable checkpoint)\n")
cat("   - seurat_with_doublets_final.rds\n")
cat("   - summary_statistics.csv\n")
cat("   - doublet_rate_barplot.pdf\n")
cat("   - cell_distribution_barplot.pdf\n")
cat("   - score_distribution.pdf\n\n")
cat("📁 analysis_output/\n")
cat("   - merged_before_doublets.rds (before doublet detection)\n")
cat("   - merged_with_doublets.rds (all cells with doublet annotations)\n")
cat("   - merged_singlets_only.rds (singlets only, ready for downstream analysis)\n\n")

cat("Next steps:\n")
cat("下一步:\n")
cat("1. Review sample QC metrics in qaqc_plots/sample_qc_metrics.csv\n")
cat("2. Examine individual sample QC plots in qaqc_plots/\n")
cat("3. Review doublet detection results in doublet_results/\n")
cat("4. Load merged_singlets_only.rds for downstream analysis:\n")
cat("   seurat_obj <- readRDS('analysis_output/merged_singlets_only.rds')\n")
cat("5. Normalization and scaling:\n")
cat("   seurat_obj <- NormalizeData(seurat_obj)\n")
cat("   seurat_obj <- FindVariableFeatures(seurat_obj)\n")
cat("   seurat_obj <- ScaleData(seurat_obj)\n")
cat("6. Cell cycle scoring (optional):\n")
cat("   s.genes <- cc.genes$s.genes\n")
cat("   g2m.genes <- cc.genes$g2m.genes\n")
cat("   seurat_obj <- CellCycleScoring(seurat_obj, s.features=s.genes, g2m.features=g2m.genes)\n")
cat("7. Dimensionality reduction:\n")
cat("   seurat_obj <- RunPCA(seurat_obj)\n")
cat("   seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)\n")
cat("8. Clustering:\n")
cat("   seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)\n")
cat("   seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)\n")
cat("9. Cell type annotation and marker analysis\n")
cat("10. Differential expression analysis\n\n")

cat("============================================================================================\n\n")
seurat_obj<- readRDS("E:/R/0301/seurat_raw_1002.rds")
saveRDS(seurat_obj,'Seon_Pyo_Hong_2023.rds')
# merged_obj <- subset(merged_obj, subset = nFeature_RNA <= 6000 & nFeature_RNA >= 200)
