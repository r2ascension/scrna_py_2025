#' 批量处理RDS文件，使用limma进行批次校正并执行下游分析
#'
#' @param input_dir 输入RDS文件目录
#' @param output_dir 输出目录
#' @param pattern RDS文件匹配模式
#' @param batch_var 批次变量名，可选
#' @param cell_type_markers 细胞类型标记基因列表，可选
#' @param run_downstream 是否进行下游分析，默认TRUE
#' @return 处理完成的文件列表
process_seurat_batch <- function(input_dir,
                                 output_dir,
                                 pattern = "\\.rds$",
                                 batch_var = NULL,
                                 cell_type_markers = NULL,
                                 run_downstream = TRUE) {
  
  # 检查包依赖
  required_packages <- c("Seurat", "limma", "ggplot2", "dplyr")
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("请安装", pkg, "包"))
    }
  }
  
  # 创建输出目录
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 获取RDS文件列表
  rds_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if(length(rds_files) == 0) {
    stop(paste0("在目录 '", input_dir, "' 中未找到匹配的RDS文件"))
  }
  
  # 初始化日志
  log_file <- file.path(output_dir, "batch_processing_log.txt")
  cat(paste0("批处理开始: ", Sys.time(), "\n"), file = log_file)
  
  # 初始化结果追踪
  results_df <- data.frame(
    file_name = character(),
    cell_count = integer(),
    processing_time = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
  
  # 检测Seurat版本
  seurat_v5 <- FALSE
  tryCatch({
    seurat_version <- packageVersion("Seurat")
    seurat_v5 <- seurat_version >= "5.0.0"
    message(paste0("检测到Seurat版本: ", seurat_version, " (V", ifelse(seurat_v5, "5+", "4或更低"), ")"))
  }, error = function(e) {
    message("无法检测Seurat版本，假设为V4或更低")
  })
  
  # 依次处理每个文件
  for(i in seq_along(rds_files)) {
    rds_file <- rds_files[i]
    file_name <- basename(rds_file)
    
    # 设置输出路径
    file_base <- tools::file_path_sans_ext(file_name)
    file_output_dir <- file.path(output_dir, file_base)
    if(!dir.exists(file_output_dir)) {
      dir.create(file_output_dir, recursive = TRUE)
    }
    
    corrected_rds <- file.path(file_output_dir, paste0(file_base, "_limma_corrected.rds"))
    
    # 记录开始时间
    start_time <- Sys.time()
    
    # 处理状态
    status <- "成功"
    
    # 记录日志
    message(paste0("\n处理文件 ", i, "/", length(rds_files), ": ", file_name))
    cat(paste0("处理文件: ", file_name, " (", i, "/", length(rds_files), ")\n"), 
        file = log_file, append = TRUE)
    
    # 尝试处理文件
    tryCatch({
      # 1. 加载RDS文件
      message("  加载Seurat对象...")
      seurat_obj <- readRDS(rds_file)
      
      # 记录细胞数量
      cell_count <- ncol(seurat_obj)
      message(paste0("  细胞数量: ", cell_count))
      
      # 2. 自动检测批次变量(如果未提供)
      if(is.null(batch_var)) {
        potential_batch_vars <- c("sample", "samples", "batch", "orig.ident", "study")
        for(var in potential_batch_vars) {
          if(var %in% colnames(seurat_obj@meta.data) && 
             length(unique(seurat_obj@meta.data[[var]])) > 1) {
            batch_var <- var
            break
          }
        }
        
        if(is.null(batch_var)) {
          message("  未找到有效的批次变量，将使用orig.ident")
          batch_var <- "orig.ident"
        }
      }
      
      message(paste0("  使用批次变量: ", batch_var))
      
      # 3. 批次效应校正(limma)
      message("  执行limma批次校正...")
      
      # 提取表达矩阵 - 只使用GetAssayData (适用于V4及以下)
      expr_matrix <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
      
      # 获取批次信息
      if(!batch_var %in% colnames(seurat_obj@meta.data)) {
        stop(paste0("批次变量 '", batch_var, "' 不存在于元数据中"))
      }
      
      batch_info <- seurat_obj@meta.data[[batch_var]]
      
      # 应用limma批次校正
      corrected_matrix <- limma::removeBatchEffect(expr_matrix, batch = batch_info)
      
      # 创建新对象以保存校正结果
      corrected_obj <- seurat_obj
      
      # 更新校正后的表达矩阵 - 只使用SetAssayData (适用于V4及以下)
      corrected_obj <- SetAssayData(corrected_obj, slot = "data", 
                                    new.data = corrected_matrix, assay = "RNA")
      
      # 添加校正信息到元数据
      corrected_obj@meta.data$batch_corrected <- TRUE
      corrected_obj@misc$correction_info <- list(
        method = "limma_removeBatchEffect",
        batch_var = batch_var,
        date = as.character(Sys.Date())
      )
      
      # 4. 下游分析(如果需要)
      if(run_downstream) {
        # 使用降维和聚类函数
        message("  执行降维和聚类...")
        corrected_obj <- run_dim_reduction(corrected_obj, harmony = FALSE)
        
        # 设置可视化输出前缀
        output_prefix <- file.path(file_output_dir, paste0(file_base, "_"))
        
        # 执行标准分析流程
        message("  执行可视化和分析...")
        corrected_obj <- run_analysis(corrected_obj, 
                                      Markers = cell_type_markers, 
                                      output_prefix = output_prefix)
      }
      
      # 5. 保存校正后的对象
      message("  保存校正后的Seurat对象...")
      saveRDS(corrected_obj, file = corrected_rds)
      
      # 6. 清理内存
      rm(seurat_obj, corrected_obj, expr_matrix, corrected_matrix)
      invisible(gc())
      
    }, error = function(e) {
      # 记录错误
      error_msg <- paste0("处理出错: ", e$message)
      message(paste0("  ", error_msg))
      cat(paste0(error_msg, "\n"), file = log_file, append = TRUE)
      status <- paste0("失败: ", e$message)
    })
    
    # 计算处理时间
    end_time <- Sys.time()
    proc_time <- difftime(end_time, start_time, units = "mins")
    
    # 更新结果表
    new_row <- data.frame(
      file_name = file_name,
      cell_count = ifelse(exists("cell_count"), cell_count, NA),
      processing_time = as.numeric(proc_time),
      status = status,
      stringsAsFactors = FALSE
    )
    
    results_df <- rbind(results_df, new_row)
    
    # 保存处理进度
    write.csv(results_df, file.path(output_dir, "processing_summary.csv"), row.names = FALSE)
    
    # 强制内存回收
    invisible(gc())
  }
  
  # 完成日志
  cat(paste0("\n批处理结束: ", Sys.time(), "\n"), file = log_file, append = TRUE)
  cat(paste0("总计处理 ", length(rds_files), " 个文件\n"), file = log_file, append = TRUE)
  cat(paste0("成功: ", sum(results_df$status == "成功"), "\n"), 
      file = log_file, append = TRUE)
  cat(paste0("失败: ", sum(results_df$status != "成功"), "\n"), 
      file = log_file, append = TRUE)
  
  # 返回处理结果
  return(results_df)
}

#' 运行降维和聚类
#'
#' @param seurat_obj Seurat对象
#' @param harmony 是否使用Harmony进行批次效应校正
#' @return 处理后的Seurat对象
run_dim_reduction <- function(seurat_obj, harmony = TRUE) {
  # 标准化和特征选择
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 4000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 40)
  
  if(harmony) {
    # 识别批次变量 - 自动检测
    batch_vars <- c()
    if("sample" %in% colnames(seurat_obj@meta.data)) batch_vars <- c(batch_vars, "sample")
    if("samples" %in% colnames(seurat_obj@meta.data)) batch_vars <- c(batch_vars, "sample")
    if("study" %in% colnames(seurat_obj@meta.data)) batch_vars <- c(batch_vars, "study")
    if("Tissue" %in% colnames(seurat_obj@meta.data)) batch_vars <- c(batch_vars, "tissue")
    if("tissue" %in% colnames(seurat_obj@meta.data)) batch_vars <- c(batch_vars, "tissue")
    
    # 确保至少有一个批次变量
    if(length(batch_vars) == 0) {
      warning("No batch variables found. Using harmony without batch correction.")
      batch_vars <- c("orig.ident")
    }
    
    cat(paste0("Running Harmony with batch variables: ", paste(batch_vars, collapse=", "), "\n"))
    
    # 直接使用Seurat对象运行Harmony整合
    seurat_obj <- RunHarmony(
      object = seurat_obj,           
      group.by.vars = c("batch_vars"),
      # theta = c(1),                  # 控制整合强度
      # lambda = c(6.7),               # 控制过度校正
      sigma = 0.05,                  # 控制聚类紧密度
      nclust = 60,                   # 聚类数量
      reduction.use = "pca",
      max_iter = 70,                 # 迭代次数 
      cool_down = 10,
      epsilon_harmony = 0.0001,
      monitor = "cluster_assignment",
      factor = 0.5,
      tol = 1e-10,
      min_delta = 0.0001,
      patience = 10,
      epsilon_cluster = 0.00001,
      early_stop = TRUE,
      dims = 1:40                    # 使用的PCA维度
    )
    
    # 计算整合质量指标（仅在有lisi包时）
    if(requireNamespace("lisi", quietly = TRUE)) {
      tryCatch({
        print("Calculating integration quality metrics...")
        
        # 获取降维结果和元数据
        embeddings <- Embeddings(seurat_obj, "harmony")
        metadata <- seurat_obj@meta.data
        
        # 检查有效的注释
        valid_columns <- c()
        if("ann_level_2" %in% colnames(metadata) && sum(!is.na(metadata$ann_level_2)) > 0) {
          valid_columns <- c(valid_columns, "ann_level_2")
        }
        if("CellType" %in% colnames(metadata) && sum(!is.na(metadata$CellType)) > 0) {
          valid_columns <- c(valid_columns, "CellType")
        }
        
        # 检查样本变量
        sample_col <- batch_vars[1]  # 使用第一个批次变量
        
        # 如果有有效的样本列和细胞类型列，计算LISI分数
        if(length(valid_columns) > 0) {
          # 过滤掉NA的细胞
          valid_cells <- !is.na(metadata[[valid_columns[1]]])
          embeddings <- embeddings[valid_cells, ]
          metadata <- metadata[valid_cells, ]
          
          # 检查并处理无效值
          embeddings_valid <- is.finite(rowSums(embeddings))
          if(sum(!embeddings_valid) > 0) {
            cat(sprintf("Removed %d cells with NA/NaN/Inf values\n", sum(!embeddings_valid)))
            embeddings <- embeddings[embeddings_valid, ]
            metadata <- metadata[embeddings_valid, ]
          }
          
          # 设置采样参数
          n_cells_sample <- min(10000, nrow(embeddings))
          n_iterations <- 5
          k <- 30
          
          # 多次抽样计算LISI scores
          ilisi_scores <- c()
          clisi_scores <- c()
          
          for(i in 1:n_iterations) {
            set.seed(i)
            cells_idx <- sample(nrow(embeddings), n_cells_sample)
            
            # 抽取对应的数据
            sampled_embeddings <- embeddings[cells_idx, ]
            sampled_metadata <- metadata[cells_idx, ]
            
            # 尝试计算LISI scores
            tryCatch({
              lisi_res <- lisi::compute_lisi(sampled_embeddings, 
                                             sampled_metadata, 
                                             c(sample_col, valid_columns[1]), 
                                             k)
              
              # 存储结果
              ilisi_scores[i] <- mean(lisi_res[, sample_col])
              clisi_scores[i] <- mean(lisi_res[, valid_columns[1]])
              
              cat(sprintf("Iteration %d completed: iLISI = %.3f, cLISI = %.3f\n", 
                          i, ilisi_scores[i], clisi_scores[i]))
            }, error = function(e) {
              cat(sprintf("Error in iteration %d: %s\n", i, e$message))
            })
          }
          
          # 计算并输出结果
          if(length(ilisi_scores) > 0) {
            cat("\nIntegration Quality Metrics:")
            cat("\n--------------------------")
            cat(sprintf("\nNumber of cells used: %d (after removing NA annotations)", nrow(embeddings)))
            cat(sprintf("\niLISI Score: %.3f ± %.3f", mean(ilisi_scores), sd(ilisi_scores)))
            cat(sprintf("\ncLISI Score: %.3f ± %.3f", mean(clisi_scores), sd(clisi_scores)))
            cat("\n\nInterpretation:")
            cat("\n- iLISI: Higher values (closer to N_batches) indicate better batch mixing")
            cat("\n- cLISI: Lower values (closer to 1) indicate better cell type separation")
            cat("\n--------------------------\n")
          } else {
            cat("\nWarning: Could not compute LISI scores due to errors in all iterations\n")
          }
        }
      }, error = function(e) {
        cat(paste("Error calculating LISI scores:", e$message, "\n"))
        cat("Continuing with analysis...\n")
      })
    }
    
    # 直接使用harmony结果进行UMAP降维
    seurat_obj <- RunUMAP(seurat_obj, 
                          reduction = "harmony", 
                          dims = 1:30, 
                          n.neighbors = 30,
                          min.dist = 0.4,
                          learning.rate = 0.2,    
                          n.epochs = 1400,        
                          spread = 1.2,
                          repulsion.strength = 1.1,
                          metric = "correlation")
    
    # 直接使用harmony结果构建细胞图
    seurat_obj <- FindNeighbors(seurat_obj, 
                                reduction = "harmony", 
                                dims = 1:30,
                                n.trees = 500,
                                k.param = 30) 
  } else {
    # 不使用Harmony时，直接用PCA进行UMAP和聚类
    seurat_obj <- RunUMAP(seurat_obj, 
                          reduction = "pca", 
                          dims = 1:30,
                          n.neighbors = 30,
                          min.dist = 0.3,
                          metric = "cosine")
    
    seurat_obj <- FindNeighbors(seurat_obj, 
                                reduction = "pca", 
                                dims = 1:30,
                                k.param = 20)
  }
  
  # 聚类分析
  seurat_obj <- FindClusters(seurat_obj,
                             algorithm = 4,  # Leiden算法
                             group.singletons = FALSE,
                             resolution = 3)
  
  return(seurat_obj)
}

#' 运行数据分析和可视化
#'
#' @param seurat_obj Seurat对象
#' @param Markers 标记基因列表
#' @param output_prefix 输出文件前缀
#' @return 处理后的Seurat对象
run_analysis <- function(seurat_obj, Markers = NULL, output_prefix = "") {
  # 创建输出目录
  output_dir <- dirname(output_prefix)
  if(!dir.exists(output_dir) && output_dir != "") {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 设置默认idents
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  
  # 首先输出UMAP聚类图
  pdf(paste0(output_prefix, "umap_clusters.pdf"), width = 10, height = 8)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE))
  dev.off()
  
  # Marker基因分析
  if(!is.null(Markers)) {
    # 确保Markers没有重复
    Markers <- unique(Markers)
    
    # 检查Markers是否存在于数据中
    valid_markers <- Markers[Markers %in% rownames(seurat_obj)]
    if(length(valid_markers) < length(Markers)) {
      missing_markers <- setdiff(Markers, valid_markers)
      cat("Warning: The following markers were not found in the data:", 
          paste(missing_markers, collapse = ", "), "\n")
    }
    
    if(length(valid_markers) > 0) {
      # 生成Feature plots
      pdf(paste0(output_prefix, "umap_FeaturePlot.pdf"), width = 8, height = 8)
      for(marker in valid_markers) {
        print(paste("Processing marker:", marker))
        print(FeaturePlot(seurat_obj, features = marker, raster = TRUE))
      }
      dev.off()
      
      # DotPlot
      pdf(paste0(output_prefix, "markers_dotplot.pdf"), width = 32, height = 18)
      dot_plot <- DotPlot(seurat_obj, 
                          features = valid_markers, 
                          group.by = "seurat_clusters",
                          cols = c("lightgrey", "red"),
                          dot.scale = 8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # 保存点图数据
      source_data <- dot_plot$data
      write.csv(source_data, paste0(output_prefix, "dotplot_data.csv"))
      print(dot_plot)
      dev.off()
      
      # 计算每个cluster中每个marker基因的表达统计
      print("Calculating expression statistics...")
      results <- list()
      
      for(cluster in unique(Idents(seurat_obj))) {
        cells_in_cluster <- WhichCells(seurat_obj, idents = cluster)
        
        for(gene in valid_markers) {
          expr_data <- GetAssayData(seurat_obj, slot = "data")[gene, cells_in_cluster]
          
          # 计算表达细胞的平均表达量
          expr_cells <- expr_data[expr_data > 0]
          mean_expr <- if(length(expr_cells) > 0) mean(expr_cells) else 0
          
          # 计算表达比例
          pct_expr <- sum(expr_data > 0) / length(expr_data) * 100
          
          results[[length(results) + 1]] <- data.frame(
            Cluster = cluster,
            Gene = gene,
            Mean_Expression = round(mean_expr, 3),
            Percent_Expressing = round(pct_expr, 2)
          )
        }
      }
      
      results_df <- do.call(rbind, results)
      write.csv(results_df, paste0(output_prefix, "cluster_marker_stats.csv"), row.names = FALSE)
    }
  }
  
  # 可视化不同分组的UMAP分布
  visualize_group_distribution <- function(group_var) {
    if(group_var %in% colnames(seurat_obj@meta.data)) {
      # 确保变量有效，移除NA值创建新列
      var_name <- paste0(group_var, "_no_na")
      seurat_obj@meta.data[[var_name]] <- seurat_obj@meta.data[[group_var]]
      seurat_obj@meta.data[[var_name]][is.na(seurat_obj@meta.data[[var_name]])] <- "Unknown"
      
      # 绘制整体分布
      pdf(paste0(output_prefix, group_var, "_umap.pdf"), width = 15, height = 10)
      print(DimPlot(seurat_obj, reduction = "umap", group.by = group_var, raster = TRUE, 
                    pt.size = 0.5) + ggtitle(group_var))
      dev.off()
      
      # 绘制各组单独高亮
      groups <- unique(seurat_obj@meta.data[[group_var]])
      groups <- groups[!is.na(groups)]  # 移除NA值
      
      pdf(paste0(output_prefix, group_var, "_umap_split.pdf"), width = 12, height = 12)
      Idents(seurat_obj) <- group_var  # 设置identity
      for(group_name in groups) {
        print(DimPlot(seurat_obj, 
                      reduction = "umap",
                      cells.highlight = CellsByIdentities(object = seurat_obj, 
                                                          idents = group_name),
                      cols = "grey",
                      pt.size = 0.5,
                      label = FALSE) +
                ggtitle(paste0(group_var, ": ", group_name)) +
                theme(legend.text = element_text(size = 12)))
      }
      dev.off()
    }
  }
  
  # 可视化不同分组的分布
  meta_columns <- colnames(seurat_obj@meta.data)
  important_vars <- c("study", "sample", "samples", "tissue", "Tissue", 
                      "ann_level_2", "CellType", "cell_type")
  
  for(var in important_vars) {
    if(var %in% meta_columns) {
      visualize_group_distribution(var)
    }
  }
  
  # 寻找cluster标记基因
  print("Finding cluster markers...")
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, paste0(output_prefix, "cluster_markers.csv"))
  
  # 选择每个cluster的top marker并生成热图
  top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
  marker_genes <- unique(top10_markers$gene)
  
  # 如果标记基因太多，分批生成热图
  seurat_obj <- ScaleData(seurat_obj, features = marker_genes)
  batch_size <- 25  # 每批次的基因数
  marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))
  
  # 生成热图
  pdf(paste0(output_prefix, "marker_heatmap_batch.pdf"), width = 32, height = 18)
  for(i in seq_along(marker_batches)) {
    print(paste("Processing marker batch", i, "of", length(marker_batches)))
    print(DoHeatmap(seurat_obj, 
                    features = marker_batches[[i]],
                    size = 24) + 
            NoLegend() +
            ggtitle(paste("Marker Genes Batch", i)) +
            theme(axis.text.y = element_text(size = 36)))
  }
  dev.off()
  
  return(seurat_obj)
}


######################################################################################################
###########################################上面是函数部分#############################################
######################################################################################################


# 主脚本文件: run_batch_processing.R

# 加载必要的包
library(Seurat)
library(limma)
library(dplyr)
library(ggplot2)

# 定义标记基因列表
Marker_Myeloid <- c('CLEC9A','XCR1','CADM1','CLNK','FLT3','ZBTB46',
                    'CLEC10A','CD1E','FCER1A','CD1D','ITGAX','CDIC','FCGR2B','PKIB',
                    'CCR7','CD83','LAMP3','CCL22','CCL17','CCL19','LAD1',
                    'LILRA4','SMPD3','SCT','IRF7','PLD4','CLEC4C',
                    'MARCO','FABP4','CYP27A1','SIGLEC1','ABCG1','PPARG',
                    'C1QA','C1QB','C1QC','HLA-DPA1','SLC40A1',
                    'FOLR2','F13A1',
                    'SPP1','HAMP','VCAN','CCR2','CCR5',
                    'FCN1',
                    'S100A12','RNASE2',
                    'LILRA5','MTSS1',
                    'TPSAB1','MS4A2','TPSB2',
                    'FCGR3B','CSF3R','CXCR1')

Marker_Endothelial <- c('S100B','ALDH1A1',
                        'GJA4','HEY1','DKK2','IGFBP3','EFNB2',
                        'ACKR1','VWF',
                        'RGCC','VWA1','IL7R','FCN3','MT1M',
                        'TFF3','LYVE1','MMRN1','CCL21',
                        'MYL9','ACTA2','TINAGL1','NOTCH3','LAMC3')

Marker_Fibroblast <- c('APOD','FGF7','COL15A1','MFAP5','PI16','CD34',
                       'MMP11','COL10A1','POSTN','LRRC15','HOPX','IGFBP5','TIMP1','MMP1','COL7A1','WNT5A','ISG15','IL7R','SFRP4','SFRP2','COMP','RGS5','PDGFRB','NDUFA4L2','NOTCH3',
                       'CXCL1','CXCL2','IL6','CEBPD','CLU','CTGF','HGF','HSPA6','DNAJB1','MYC','AFT4','PLAU','CHI3L1','MMP3','IL1R1','IL13RA2','TNFSF11','MMP10','OSMR','IL11','STRA6','FAP','WNT2','TWIST1','IL24',
                       'ACTG2','HHIP','CNN1',
                       'MYH11','ACTA2','TAGLN',
                       'KRT18','SLPI','UPK3B','MSLN','CALB2','WT1','KLK11','ITLN1',
                       'WSB1','DDX17','CTNNB1',
                       'RBP1','STAR','STMN1',
                       'CXCL12','CD74','HLA.DRB1','HLA.DRA',
                       'ADAMDEC1','CCL8','APOE','APOC1',
                       'LIMCH1','A2M','ADH1B',
                       'PRG4','CRTAC1',
                       'CXCL14','VSTM2A','SOX6','COL4A5','COL4A6','TSLP','FRZB','BMP5','BMP2','CPM','F3')

Marker_Epithelial <- c('TP63','KRT5',
                       'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1',
                       'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24',
                       'MUC5AC','SPDEF','LYPD2','ITLN1',
                       'ASCL1','GRP',
                       'POU2F3','ASCL2','CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB',
                       'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
                       'SFTPC','LAMP3','MF5D2A','C8orf4','C11orf96',
                       'VIM','SOX9',
                       'KRT14','MYH11','ACTA2',
                       'DMBT1','RNASE1',
                       'MUC5B','SPDEF',
                       'LYZ','LTF',
                       'SFTPB','SCGB3A2','SFTA2')

Marker_T <- c('KLRD1','FCGR3A','GNLY','TYROBP','FCER1G','KLRC1','FGFBP2','SPON2','MYOM2',
              'TRDC','KRT86',
              'GATA3','IL5','AREG','HPGDS',
              'IL23R','RORC','LST1','PCDH9','TNFSF11',
              'TRDC','TRGC1','TRGC2',
              'CD4','CD28','CD40LG','TRAT1','TNFRSF25',
              'CD8A','CD8B','TRGC2',
              'CCR7','TCF7','LEF1','SELL',
              'CD28','IL7R','CCR6',
              'GATA3','IL4','IL13',
              'IL17A','CCL20',
              'CCR7','TCF7','LEF1','SELL',
              'TBX21','EOMES','GZMK','KLRG1',
              'ITGA1','CD8A','CD8B','CCR6','IL7R',
              'TBX21','GZMB','GZMH','FGFBP2',
              'ZNF683','IFNG','CCL4L2','PDCD1',
              'KLRB1','IL7R','NCR3','CEBPD',
              'FOXP3','IL2RA','IKZF2','TNFRSF4',
              'MKI67','TOP2A','TK1','CENPW')

Marker_SMC <- c('RGS5','CD36','NOTCH3',
                'SCIN','EPAS1','HLA.C','IGKC','PTP4A3','FN1','COL18A1','WFDC1','IGHG4','IGHG1',
                'RERGL','PLN','SORBS2','DSTN','TSC22D1','BCAM','C11orf96','FILIP1','WTIP','NRGN',
                'TMEM176B','ANGPTL1','CFH','FHL1','VCAN','GGT5','C1S','COL6A3','STEAP4','ADGRL3')

Marker_B <- c('IGHD','TCL1A',
              'MS4A1','BANK1',
              'MKI67','TOP2A',
              'IGHA1','IGHA2',
              'IGHGP','IGHG1')

# 合并所有标记基因
all_markers <- list(
  Myeloid = Marker_Myeloid,
  Endothelial = Marker_Endothelial,
  Fibroblast = Marker_Fibroblast,
  Epithelial = Marker_Epithelial,
  T = Marker_T,
  SMC = Marker_SMC,
  B = Marker_B
)

# 根据细胞类型选择特定的标记基因
select_markers <- function(cell_type = NULL) {
  if(is.null(cell_type)) {
    # 如果未指定细胞类型，返回所有标记基因的并集
    return(unique(unlist(all_markers)))
  } else if(cell_type %in% names(all_markers)) {
    # 如果指定了具体的细胞类型，返回该类型的标记基因
    return(all_markers[[cell_type]])
  } else {
    warning(paste0("未找到细胞类型: ", cell_type, "，返回所有标记基因"))
    return(unique(unlist(all_markers)))
  }
}

# 设置输入和输出目录
input_dir <- "E:/R/0301/rds"   # RDS文件所在目录
output_dir <- "E:/R/0301/newrds"   # 结果输出目录

# 执行批处理
# 可以根据需要选择特定的细胞类型标记基因
results <- process_seurat_batch(
  input_dir = input_dir,
  output_dir = output_dir,
  pattern = "\\.rds$",                    # RDS文件匹配模式
  batch_var = 'study',                       # 自动检测批次变量
  cell_type_markers = select_markers(),   # 使用所有标记基因
  run_downstream = TRUE                   # 运行下游分析
)

# 查看处理结果
print(results)

# # 如果您想要处理Fibroblast特定的RDS文件，并只使用成纤维细胞的标记基因
# fibroblast_results <- process_seurat_batch(
#   input_dir = input_dir,
#   output_dir = file.path(output_dir, "fibroblast_specific"),
#   pattern = "fibroblast.*\\.rds$",           # 匹配包含fibroblast的RDS文件
#   batch_var = "sample",                      # 指定批次变量
#   cell_type_markers = select_markers("Fibroblast"), # 使用成纤维细胞标记基因
#   run_downstream = TRUE
# )