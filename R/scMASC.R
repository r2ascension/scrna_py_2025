#' 基于MASC的单细胞组织间两两比较分析函数（优化版）
#' 
#' 使用混合效应模型对所有组织对进行两两比较分析
#' 
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名，默认"Annotation"
#' @param sample_col 样本/供体列名，默认"sample"，作为随机效应
#' @param contrast_col 组织/条件列名，默认"tissue" 
#' @param fixed_effects_cols 固定效应列名向量，默认NULL
#' @param exclude_filter 排除条件的表达式字符串，例如："nFeature_RNA < 200 | percent.mt > 20"
#' @param output_dir 输出目录
#' @param p_threshold p值显著性阈值，默认0.05
#' @param fdr_threshold FDR显著性阈值，默认0.1
#' @param min_cells 每种细胞类型的最小细胞数，默认10
#' @param min_prop 考虑分析的最小细胞比例，默认0.001
#' @param color_palette 用于绘图的颜色列表
#' @param save_models 是否保存MASC模型，默认FALSE
#' @param min_samples 每个组织的最小样本数，默认2
#' @param specific_pairs 指定要比较的组织对的列表，默认NULL表示比较所有可能的组织对
#' @param seed 随机数种子
#' @return 分析结果列表
scPairwiseMASCAnalysis <- function(seurat_obj, 
                                   cell_type_col = "Annotation", 
                                   sample_col = "sample", 
                                   contrast_col = "tissue", 
                                   fixed_effects_cols = NULL,
                                   exclude_filter = NULL, 
                                   output_dir = "pairwise_MASC_analysis",
                                   p_threshold = 0.05,
                                   fdr_threshold = 0.1,
                                   min_cells = 10,
                                   min_prop = 0.001,
                                   color_palette = NULL,
                                   save_models = FALSE,
                                   min_samples = 2,
                                   specific_pairs = NULL,
                                   seed = 42) {
  
  # 设置随机数种子
  set.seed(seed)
  
  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 记录函数开始时间
  start_time <- Sys.time()
  message(paste("Starting pairwise MASC cell association analysis:", format(start_time)))
  
  #--------------------------------------------------
  # 1. 数据准备
  #--------------------------------------------------
  message("Preparing data...")
  
  # 确保必要的包已加载
  required_packages <- c("Seurat", "dplyr", "ggplot2", "reshape2", 
                         "tidyr", "lme4", "cowplot", "utils")
  
  # 检查并加载包
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it."))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # 检查pheatmap包
  run_heatmap <- requireNamespace("pheatmap", quietly = TRUE)
  
  # 获取元数据
  meta_data <- seurat_obj@meta.data
  
  # 检查必要的列
  required_cols <- c(cell_type_col, sample_col, contrast_col)
  if (!is.null(fixed_effects_cols)) {
    required_cols <- c(required_cols, fixed_effects_cols)
  }
  
  if (!all(required_cols %in% colnames(meta_data))) {
    stop(paste("Missing required columns:", 
               paste(setdiff(required_cols, colnames(meta_data)), collapse = ", ")))
  }
  
  # 应用筛选条件 - 修复筛选逻辑问题
  original_count <- nrow(meta_data)
  if (!is.null(exclude_filter) && exclude_filter != "") {
    tryCatch({
      # 正确解析exclude_filter为排除条件
      filter_expr <- parse(text = exclude_filter)
      # 评估条件，找出要排除的行
      rows_to_exclude <- eval(filter_expr, meta_data)
      # 保留不符合排除条件的行
      filtered_meta <- meta_data[!rows_to_exclude, ]
      
      excluded_count <- original_count - nrow(filtered_meta)
      message(paste("Cells before filtering:", original_count))
      message(paste("Cells excluded by filter:", excluded_count, 
                    "(", round(excluded_count/original_count * 100, 2), "%)"))
      message(paste("Cells after filtering:", nrow(filtered_meta)))
      
      # 检查是否过滤掉所有细胞
      if (nrow(filtered_meta) == 0) {
        warning("Filter condition excluded all cells! Using original data.")
        filtered_meta <- meta_data
      }
    }, error = function(e) {
      message(paste("Error in filter expression:", e$message))
      message("Using all cells without filtering")
      filtered_meta <- meta_data
    })
  } else {
    filtered_meta <- meta_data
    message(paste("No filter applied. Using all", nrow(filtered_meta), "cells."))
  }
  
  # 移除NA值
  na_cols <- c(cell_type_col, sample_col, contrast_col)
  if (!is.null(fixed_effects_cols)) {
    na_cols <- c(na_cols, fixed_effects_cols)
  }
  
  # 创建逻辑向量指示哪些行有NA值
  has_na <- apply(filtered_meta[, na_cols, drop = FALSE], 1, function(x) any(is.na(x)))
  if (any(has_na)) {
    message(paste("Removing", sum(has_na), "cells with NA values in key columns"))
    filtered_meta <- filtered_meta[!has_na, ]
  }
  
  # 转换因子列
  filtered_meta[[cell_type_col]] <- as.factor(filtered_meta[[cell_type_col]])
  filtered_meta[[sample_col]] <- as.factor(filtered_meta[[sample_col]])
  filtered_meta[[contrast_col]] <- as.factor(filtered_meta[[contrast_col]])
  
  if (!is.null(fixed_effects_cols)) {
    for (col in fixed_effects_cols) {
      filtered_meta[[col]] <- as.factor(filtered_meta[[col]])
    }
  }
  
  # 获取所有组织水平
  all_tissues <- levels(filtered_meta[[contrast_col]])
  n_tissues <- length(all_tissues)
  
  if (n_tissues < 2) {
    stop("Need at least 2 tissues for comparison")
  }
  
  message(paste("Found", n_tissues, "tissues/conditions:", paste(all_tissues, collapse=", ")))
  
  # 统计每个组织的样本数量（修复BUG：使用唯一样本计数）
  samples_by_tissue <- tapply(filtered_meta[[sample_col]], 
                              filtered_meta[[contrast_col]], 
                              function(x) length(unique(x)))
  
  message("Sample count per tissue:")
  for (tissue in names(samples_by_tissue)) {
    message(paste(" -", tissue, ":", samples_by_tissue[tissue], "samples"))
  }
  
  # 检查样本总数
  n_samples <- length(unique(filtered_meta[[sample_col]]))
  n_cell_types <- length(unique(filtered_meta[[cell_type_col]]))
  
  message(paste("Total unique samples/donors:", n_samples))
  message(paste("Number of cell types:", n_cell_types))
  
  # 检查样本数量是否足够进行统计分析
  if (n_samples < 4) {
    warning("Very few samples (<4). Statistical results may not be reliable.")
  }
  
  # 检查每个组织是否有足够样本
  tissues_with_enough_samples <- names(samples_by_tissue)[samples_by_tissue >= min_samples]
  if (length(tissues_with_enough_samples) < 2) {
    stop(paste("Not enough tissues with", min_samples, "or more samples for comparison"))
  }
  
  if (length(tissues_with_enough_samples) < n_tissues) {
    warning(paste("The following tissues have fewer than", min_samples, "samples and will be excluded:"))
    warning(paste(" ", setdiff(names(samples_by_tissue), tissues_with_enough_samples), collapse="\n  "))
    
    # 过滤掉样本不足的组织
    filtered_meta <- filtered_meta[filtered_meta[[contrast_col]] %in% tissues_with_enough_samples, ]
    filtered_meta[[contrast_col]] <- factor(filtered_meta[[contrast_col]])  # 重新设置因子水平
  }
  
  # 更新组织列表
  all_tissues <- levels(filtered_meta[[contrast_col]])
  n_tissues <- length(all_tissues)
  
  #--------------------------------------------------
  # 2. 计算细胞计数和比例（用于可视化和过滤）
  #--------------------------------------------------
  message("Calculating cell counts and proportions...")
  
  # 创建计数表
  counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
  counts_matrix <- as.matrix(counts_table)
  
  # 添加样本的对比信息
  sample_contrast_df <- unique(filtered_meta[, c(sample_col, contrast_col)])
  sample_to_contrast <- setNames(
    as.character(sample_contrast_df[[contrast_col]]), 
    as.character(sample_contrast_df[[sample_col]])
  )
  
  # 检查每种细胞类型的总细胞数
  cell_type_totals <- rowSums(counts_matrix)
  low_count_cell_types <- names(cell_type_totals)[cell_type_totals < min_cells]
  
  if (length(low_count_cell_types) > 0) {
    message(paste("Warning: The following cell types have fewer than", min_cells, "cells:"))
    message(paste(low_count_cell_types, collapse = ", "))
    message("These cell types will be excluded from analysis.")
    
    # 过滤低计数的细胞类型
    filtered_meta <- filtered_meta[!filtered_meta[[cell_type_col]] %in% low_count_cell_types, ]
    # 重建计数矩阵
    counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
    counts_matrix <- as.matrix(counts_table)
  }
  
  # 计算每个样本的细胞比例
  if (nrow(counts_matrix) > 0 && ncol(counts_matrix) > 0) {
    props_by_sample <- prop.table(counts_matrix, margin = 2)
    
    # 过滤低比例的细胞类型
    cell_type_mean_props <- rowMeans(props_by_sample)
    low_prop_cell_types <- names(cell_type_mean_props)[cell_type_mean_props < min_prop]
    
    if (length(low_prop_cell_types) > 0) {
      message(paste("Warning: The following cell types have mean proportion <", min_prop, ":"))
      message(paste(low_prop_cell_types, collapse = ", "))
      message("These cell types will be excluded from analysis.")
      
      # 过滤低比例的细胞类型
      filtered_meta <- filtered_meta[!filtered_meta[[cell_type_col]] %in% low_prop_cell_types, ]
      # 重建计数矩阵和比例矩阵
      counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
      counts_matrix <- as.matrix(counts_table)
      props_by_sample <- prop.table(counts_matrix, margin = 2)
    }
    
    # 转换为长格式用于ggplot绘图
    props_df <- reshape2::melt(props_by_sample)
    colnames(props_df) <- c("CellType", "Sample", "Proportion")
    
    # 添加对比信息到props_df
    props_df$Contrast <- sample_to_contrast[as.character(props_df$Sample)]
    
    # 检查是否有NA值
    if (any(is.na(props_df$Contrast))) {
      warning("Some samples couldn't be matched to contrasts. Check sample_to_contrast mapping.")
      # 移除NA值
      props_df <- props_df[!is.na(props_df$Contrast), ]
    }
    
    # 保存比例数据
    write.csv(props_df, file.path(output_dir, "cell_proportions.csv"), row.names = FALSE)
  } else {
    stop("No data remains after filtering. Please check your filtering criteria.")
  }
  
  #--------------------------------------------------
  # 3. 可视化细胞组成
  #--------------------------------------------------
  message("Generating cell composition visualizations...")
  
  # 定义一致的颜色调色板
  if(is.null(color_palette)) {
    n_cell_types <- length(unique(props_df$CellType))
    
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector <- unlist(mapply(RColorBrewer::brewer.pal, 
                                  qual_col_pals$maxcolors, 
                                  rownames(qual_col_pals)))
      color_palette <- sample(col_vector, n_cell_types)
    } else {
      # 基础R颜色
      color_palette <- sample(colors(), n_cell_types)
    }
    
    names(color_palette) <- unique(props_df$CellType)
  }
  
  # 样本级别的堆叠柱状图
  p1 <- ggplot(props_df, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right") +
    labs(title = "Cell Type Composition by Sample",
         y = "Proportion", 
         x = "Sample", 
         fill = "Cell Type") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_palette)
  
  # 添加样本分组
  p1 <- p1 + facet_grid(. ~ Contrast, scales = "free_x", space = "free")
  
  # 保存图像
  ggsave(file.path(output_dir, "cell_composition_by_sample.pdf"), 
         p1, width = 14, height = 8)
  
  # 对比级别的堆叠柱状图
  props_by_contrast <- props_df %>%
    group_by(CellType, Contrast) %>%
    summarise(Proportion = mean(Proportion), .groups = "drop")
  
  p2 <- ggplot(props_by_contrast, aes(x = Contrast, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right") +
    labs(title = "Average Cell Type Composition by Contrast",
         y = "Proportion", 
         x = "Contrast Group", 
         fill = "Cell Type") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_palette)
  
  # 保存图像
  ggsave(file.path(output_dir, "cell_composition_by_contrast.pdf"), 
         p2, width = 10, height = 8)
  
  # 为每个细胞类型创建箱线图
  props_df %>%
    group_by(CellType) %>%
    group_walk(function(cell_data, key) {
      cell_type <- key$CellType
      p <- ggplot(cell_data, aes(x = Contrast, y = Proportion, fill = Contrast)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
        theme_bw() +
        labs(title = paste0(cell_type, " - Proportion by Condition"),
             y = "Proportion", 
             x = "Contrast Group") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 12, face = "bold"),
              legend.position = "none") +
        scale_y_continuous(labels = scales::percent)
      
      ggsave(file.path(output_dir, paste0(cell_type, "_boxplot.pdf")), 
             p, width = 8, height = 6)
    })
  
  #--------------------------------------------------
  # 4. 定义MASC函数
  #--------------------------------------------------
  
  # 检查是否在环境中找到MASC
  if (!exists("MASC", envir = .GlobalEnv)) {
    MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                     verbose = FALSE, save_models = FALSE, save_model_dir = NULL) {
      # Check inputs
      if (is.factor(dataset[[contrast]]) == FALSE) {
        stop("Specified contrast term is not coded as a factor in dataset")
      }
      
      # 将cluster转换为character
      cluster <- as.character(cluster)
      # 创建设计矩阵
      designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
      dataset <- cbind(designmat, dataset)
      
      # 创建输出列表
      res <- vector(mode = "list", length = length(unique(cluster)))
      names(res) <- attributes(designmat)$dimnames[[2]]
      
      # 创建模型公式
      if (!is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                              paste0("(1|", random_effects, ")", collapse = " + ")),
                            collapse = " + ")
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
        }
      } else if (!is.null(fixed_effects) && is.null(random_effects)) {
        model_rhs <- paste0(fixed_effects, collapse = " + ")
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
          stop("No random effects specified")
        }
      } else if (is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
        }
      } else {
        model_rhs <- "1" # 仅包含截距
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
          stop("No random or fixed effects specified")
        }
      }
      
      # 初始化列表存储每个簇的模型对象
      cluster_models <- vector(mode = "list",
                               length = length(attributes(designmat)$dimnames[[2]]))
      names(cluster_models) <- attributes(designmat)$dimnames[[2]]
      
      # 为每个簇运行嵌套混合效应模型
      for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
        test_cluster <- attributes(designmat)$dimnames[[2]][i]
        if (verbose == TRUE) {
          message(paste("Creating logistic mixed models for", test_cluster))
        }
        
        # 添加检查：确保每个对比水平中该簇有足够的细胞
        cluster_by_contrast <- table(dataset[[contrast]], dataset[[test_cluster]])
        if (nrow(cluster_by_contrast) < 2 || ncol(cluster_by_contrast) < 2) {
          warning(paste("Invalid data structure for cluster", test_cluster, "- skipping"))
          cluster_models[[i]]$error <- "Invalid data structure"
          next
        }
        
        # 检查每个对比水平中是否有足够的0和1（在每个对比水平中都需要有细胞属于和不属于该簇）
        has_enough_data <- TRUE
        for (level in levels(dataset[[contrast]])) {
          contrast_subset <- dataset[dataset[[contrast]] == level, ]
          cluster_count <- sum(contrast_subset[[test_cluster]])
          
          # 检查在该对比水平中，是否有足够的细胞属于和不属于该簇
          if (cluster_count < 1 || cluster_count >= nrow(contrast_subset) || nrow(contrast_subset) < 1) {
            has_enough_data <- FALSE
            break
          }
        }
        
        if (!has_enough_data) {
          warning(paste("Not enough variation for cluster", test_cluster, "in at least one contrast level - skipping"))
          cluster_models[[i]]$error <- "Insufficient variation"
          next
        }
        
        null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                       model_rhs), collapse = ""))
        full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                       model_rhs), collapse = ""))
        # 运行空模型和完整混合效应模型
        tryCatch({
          null_model <- lme4::glmer(formula = null_fm, data = dataset,
                                    family = binomial, nAGQ = 1, verbose = 0,
                                    control = glmerControl(optimizer = "bobyqa"))
          full_model <- lme4::glmer(formula = full_fm, data = dataset,
                                    family = binomial, nAGQ = 1, verbose = 0,
                                    control = glmerControl(optimizer = "bobyqa"))
          model_lrt <- anova(null_model, full_model)
          
          # 计算对比项beta的置信区间
          contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
          contrast_ci <- confint.merMod(full_model, method = "Wald",
                                        parm = contrast_lvl2)
          
          # 保存模型对象到列表
          cluster_models[[i]]$null_model <- null_model
          cluster_models[[i]]$full_model <- full_model
          cluster_models[[i]]$model_lrt <- model_lrt
          cluster_models[[i]]$confint <- contrast_ci
        }, error = function(e) {
          warning(paste("Model for cluster", test_cluster, "failed:", e$message))
          cluster_models[[i]]$error <- e$message
        })
      }
      
      # 组织结果为输出数据框
      output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                           size = colSums(designmat))
      
      # 检查哪些模型成功运行
      successful_models <- sapply(cluster_models, function(x) !is.null(x$model_lrt))
      
      if (any(successful_models)) {
        output$model.pvalue <- NA
        output$model.pvalue[successful_models] <- sapply(
          cluster_models[successful_models], 
          function(x) x$model_lrt[["Pr(>Chisq)"]][2]
        )
        
        # 添加FDR校正
        output$FDR <- p.adjust(output$model.pvalue, method = "BH")
        
        # 添加优势比和置信区间
        contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
        output[[paste(contrast_lvl2, "OR", sep = ".")]] <- NA
        output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- NA
        output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- NA
        
        for (i in which(successful_models)) {
          output[[paste(contrast_lvl2, "OR", sep = ".")]][i] <- 
            exp(fixef(cluster_models[[i]]$full_model)[[contrast_lvl2]])
          output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]][i] <- 
            exp(cluster_models[[i]]$confint[contrast_lvl2, "2.5 %"])
          output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]][i] <- 
            exp(cluster_models[[i]]$confint[contrast_lvl2, "97.5 %"])
        }
      }
      
      # 保存模型对象并返回结果
      if (save_models == TRUE) {
        if (is.null(save_model_dir)) {
          save_dir <- getwd()
        } else {
          save_dir <- save_model_dir
          dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        }
        saveRDS(cluster_models, file = file.path(save_dir, "masc_models.rds"))
        message(paste("Models saved to", file.path(save_dir, "masc_models.rds")))
        return(list(results = output, models = cluster_models))
      } else {
        return(list(results = output, models = cluster_models))
      }
    }
  }
  
  #--------------------------------------------------
  # 5. 构建所有可能的组织对
  #--------------------------------------------------
  message("Building tissue/condition pairs for comparison...")
  
  # 生成所有可能的组织对组合
  if (is.null(specific_pairs)) {
    tissue_pairs <- combn(all_tissues, 2, simplify = FALSE)
    message(paste("Generated", length(tissue_pairs), "tissue pairs for comparison"))
  } else {
    # 验证指定的组织对
    tissue_pairs <- specific_pairs
    valid_pairs <- sapply(tissue_pairs, function(pair) {
      all(pair %in% all_tissues) && length(pair) == 2
    })
    
    if (!all(valid_pairs)) {
      invalid_indices <- which(!valid_pairs)
      warning(paste("Invalid tissue pairs specified at indices:", 
                    paste(invalid_indices, collapse=", ")))
      tissue_pairs <- tissue_pairs[valid_pairs]
    }
    
    message(paste("Using", length(tissue_pairs), "specified tissue pairs for comparison"))
  }
  
  if (length(tissue_pairs) == 0) {
    stop("No valid tissue pairs for comparison")
  }
  
  #--------------------------------------------------
  # 6. 对每个组织对运行MASC分析
  #--------------------------------------------------
  
  # 创建结果存储列表
  pairwise_results <- list()
  
  # 创建进度条
  pb <- utils::txtProgressBar(min = 0, max = length(tissue_pairs), style = 3)
  
  # 对每个组织对进行分析
  for (i in 1:length(tissue_pairs)) {
    pair <- tissue_pairs[[i]]
    tissue1 <- pair[1]
    tissue2 <- pair[2]
    
    comparison_name <- paste(tissue2, "vs", tissue1)
    message(paste("\nComparing", comparison_name, "(",i, "of", length(tissue_pairs),")"))
    
    # 创建只包含这两个组织的数据子集
    pair_meta <- filtered_meta[filtered_meta[[contrast_col]] %in% pair, ]
    
    # 修复: 计算每个组织的唯一样本数 - 使用简单的计数方法避免NA值
    sample_tissue_pairs <- unique(pair_meta[, c(sample_col, contrast_col)])
    samples_per_tissue <- table(sample_tissue_pairs[[contrast_col]])
    
    # 确保只有当前组织对
    samples_per_tissue <- samples_per_tissue[names(samples_per_tissue) %in% pair]
    
    message(paste(" Sample counts:", 
                  paste(names(samples_per_tissue), "=", samples_per_tissue, collapse=", ")))
    
    # 确保有足够样本
    if (length(samples_per_tissue) != 2 || any(samples_per_tissue < min_samples)) {
      skip_reason <- if(length(samples_per_tissue) != 2) 
        "missing data for one tissue" 
      else 
        "insufficient samples"
      warning(paste("Skipping", comparison_name, "-", skip_reason))
      next
    }
    
    # 重新编码对比因子以确保tissue1是参照水平
    pair_meta[[contrast_col]] <- factor(pair_meta[[contrast_col]], levels = pair)
    
    # 修复：检查每个细胞类型在每个组织中的细胞数
    cell_counts_by_tissue <- table(pair_meta[[cell_type_col]], pair_meta[[contrast_col]])
    
    # 检查每个细胞类型在两个组织中是否都有足够的计数
    zero_or_low_count_cells <- rownames(cell_counts_by_tissue)[
      apply(cell_counts_by_tissue, 1, function(x) any(x < min_cells))
    ]
    
    if (length(zero_or_low_count_cells) > 0) {
      message(paste("Warning: The following cell types have insufficient counts (<", min_cells, ") in at least one tissue and will be excluded from this comparison:"))
      message(paste(zero_or_low_count_cells, collapse = ", "))
      
      # 过滤掉这些细胞类型
      pair_meta <- pair_meta[!pair_meta[[cell_type_col]] %in% zero_or_low_count_cells, ]
    }
    
    # 检查过滤后是否还有足够细胞类型进行比较
    if (length(unique(pair_meta[[cell_type_col]])) < 1) {
      warning(paste("Skipping", comparison_name, "- no cell types have sufficient counts in both tissues"))
      next
    }
    
    # 将数据转换为MASC所需格式
    masc_dataset <- pair_meta[, c(cell_type_col, sample_col, contrast_col)]
    
    # 添加固定效应列
    if (!is.null(fixed_effects_cols)) {
      masc_dataset <- cbind(masc_dataset, pair_meta[, fixed_effects_cols, drop = FALSE])
    }
    
    # 运行MASC分析
    pair_results <- tryCatch({
      MASC(dataset = masc_dataset,
           cluster = masc_dataset[[cell_type_col]],
           contrast = contrast_col,
           random_effects = sample_col,
           fixed_effects = fixed_effects_cols,
           verbose = TRUE,
           save_models = save_models,
           save_model_dir = ifelse(save_models, 
                                   file.path(output_dir, paste0("models_", gsub(" ", "_", comparison_name))), 
                                   NULL))
    }, error = function(e) {
      warning(paste("MASC analysis failed for", comparison_name, ":", e$message))
      return(NULL)
    })
    
    # 如果分析成功
    if (!is.null(pair_results)) {
      # 提取结果并添加组织信息
      results_df <- pair_results$results
      
      # 检查结果是否为空
      if (nrow(results_df) > 0) {
        results_df$Tissue1 <- tissue1  # 参照水平
        results_df$Tissue2 <- tissue2  # 对比水平
        results_df$Comparison <- comparison_name
        
        # 将结果添加到列表
        pairwise_results[[comparison_name]] <- list(
          results = results_df,
          models = if(save_models) pair_results$models else NULL
        )
      } else {
        warning(paste("No results returned for", comparison_name))
      }
    }
    
    # 更新进度条
    utils::setTxtProgressBar(pb, i)
  }
  
  # 关闭进度条
  close(pb)
  
  # 检查是否有成功的分析
  if (length(pairwise_results) == 0) {
    stop("No successful comparisons were made")
  }
  
  #--------------------------------------------------
  # 7. 整合所有比较结果
  #--------------------------------------------------
  # 在整合所有比较结果的部分
  # 整合所有比较结果
  message("\nIntegrating all pairwise comparison results...")
  
  # 提取和组合所有结果
  if (length(pairwise_results) > 0) {
    # 确保只处理有效的结果
    valid_results <- sapply(pairwise_results, function(x) {
      !is.null(x$results) && nrow(x$results) > 0
    })
    
    if (any(valid_results)) {
      # 组合所有有效的结果
      all_results_list <- lapply(names(pairwise_results)[valid_results], function(comparison) {
        results_df <- pairwise_results[[comparison]]$results
        
        # 标准化列名
        or_cols <- grep("\\.OR$", colnames(results_df), value = TRUE)
        ci_lower_cols <- grep("\\.OR\\.95pct\\.ci\\.lower$", colnames(results_df), value = TRUE)
        ci_upper_cols <- grep("\\.OR\\.95pct\\.ci\\.upper$", colnames(results_df), value = TRUE)
        
        if (length(or_cols) > 0) {
          colnames(results_df)[colnames(results_df) == or_cols[1]] <- "OR"
          colnames(results_df)[colnames(results_df) == ci_lower_cols[1]] <- "CI_Lower"
          colnames(results_df)[colnames(results_df) == ci_upper_cols[1]] <- "CI_Upper"
        }
        
        # 添加比较名称
        results_df$Comparison <- comparison
        return(results_df)
      })
      
      # 合并所有结果
      all_results <- do.call(rbind, all_results_list)
      
      # 重命名列以匹配细胞类型
      cluster_to_celltype <- function(cluster_name) {
        # 从"clusterCellType"中提取"CellType"
        gsub("^cluster", "", cluster_name)
      }
      
      all_results$CellType <- cluster_to_celltype(all_results$cluster)
      
      # 重要修复：对所有p值进行一次性FDR校正
      all_results$FDR <- p.adjust(all_results$model.pvalue, method = "BH")
      
      # 添加显著性标记
      all_results$Significant_P <- all_results$model.pvalue < p_threshold
      all_results$Significant_FDR <- all_results$FDR < fdr_threshold  # 使用统一阈值
      
      # 添加显著性星号
      all_results$Significance <- ifelse(all_results$FDR < 0.001, "***",
                                         ifelse(all_results$FDR < 0.01, "**",
                                                ifelse(all_results$FDR < 0.05, "*", 
                                                       ifelse(all_results$FDR < 0.1, ".", "ns"))))
      
      # 保存整合结果
      write.csv(all_results, file.path(output_dir, "all_pairwise_masc_results.csv"), row.names = FALSE)
    } else {
      all_results <- empty_results
      warning("No valid results in any comparison")
    }
  } else {
    all_results <- empty_results
    warning("No pairwise results available")
  }
  # 检查是否有结果
  if (nrow(all_results) > 0) {
    # 重命名列以匹配细胞类型
    cluster_to_celltype <- function(cluster_name) {
      # 从"clusterCellType"中提取"CellType"
      gsub("^cluster", "", cluster_name)
    }
    
    all_results$CellType <- cluster_to_celltype(all_results$cluster)
    
    # 添加显著性标记
    all_results$Significant_FDR <- all_results$FDR < fdr_threshold
    all_results$Significant_P <- all_results$model.pvalue < p_threshold
    
    # 添加显著性星号
    all_results$Significance <- ifelse(all_results$FDR < 0.001, "***",
                                       ifelse(all_results$FDR < 0.01, "**",
                                              ifelse(all_results$FDR < 0.05, "*", 
                                                     ifelse(all_results$FDR < 0.1, ".", "ns"))))
    
    # 保存整合结果
    write.csv(all_results, file.path(output_dir, "all_pairwise_masc_results.csv"), row.names = FALSE)
    
    # 拆分结果为每个细胞类型的数据框
    if (length(unique(all_results$CellType)) > 0) {
      cell_type_results <- split(all_results, all_results$CellType)
      
      # 为每个细胞类型保存结果
      for (cell_type in names(cell_type_results)) {
        write.csv(cell_type_results[[cell_type]], 
                  file.path(output_dir, paste0(cell_type, "_pairwise_results.csv")), 
                  row.names = FALSE)
      }
    }
    
    #--------------------------------------------------
    # 8. 显示显著结果摘要
    #--------------------------------------------------
    
    # 找出显著的比较
    sig_results <- all_results[all_results$Significant_FDR, ]
    
    if (nrow(sig_results) > 0) {
      message(paste("Found", nrow(sig_results), "significant cell type differences (FDR <", fdr_threshold, ")"))
      
      # 按细胞类型分组摘要
      sig_by_celltype <- split(sig_results, sig_results$CellType)
      
      for (cell_type in names(sig_by_celltype)) {
        message(paste("\nCell type:", cell_type))
        ct_results <- sig_by_celltype[[cell_type]]
        
        for (i in 1:nrow(ct_results)) {
          row <- ct_results[i, ]
          direction <- ifelse(row$OR > 1, "enriched", "depleted")
          message(paste0("  - ", row$Comparison, ": ", 
                         direction, " (OR = ", round(row$OR, 2), 
                         ", FDR = ", format(row$FDR, digits = 3), ")"))
        }
      }
    } else {
      message("No significant cell type differences found")
    }
    
    #--------------------------------------------------
    # 9. 可视化结果
    #--------------------------------------------------
    message("\nGenerating visualizations for pairwise comparisons...")
    
    # 为每个细胞类型创建森林图
    for (cell_type in unique(all_results$CellType)) {
      # 提取该细胞类型的结果
      ct_results <- all_results[all_results$CellType == cell_type, ]
      
      # 跳过没有OR的结果
      ct_results <- ct_results[!is.na(ct_results$OR), ]
      
      if (nrow(ct_results) > 0) {
        # 计算对数优势比和置信区间
        ct_results$log_OR <- log(ct_results$OR)
        ct_results$log_CI_Lower <- log(ct_results$CI_Lower)
        ct_results$log_CI_Upper <- log(ct_results$CI_Upper)
        
        # 处理无限值和NA值
        ct_results$log_OR[!is.finite(ct_results$log_OR)] <- NA
        ct_results$log_CI_Lower[!is.finite(ct_results$log_CI_Lower)] <- NA
        ct_results$log_CI_Upper[!is.finite(ct_results$log_CI_Upper)] <- NA
        
        # 检查是否有足够的非NA值进行绘图
        if (sum(!is.na(ct_results$log_OR)) > 0) {
          # 按对数优势比排序
          ct_results <- ct_results[order(ct_results$log_OR), ]
          
          # 创建比较因子，保持排序
          ct_results$Comparison <- factor(ct_results$Comparison, levels = ct_results$Comparison)
          
          # 森林图
          p_forest <- ggplot(ct_results, aes(x = log_OR, y = Comparison, color = Significant_FDR)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
            geom_point(size = 3) +
            geom_errorbarh(aes(xmin = log_CI_Lower, xmax = log_CI_Upper), height = 0.2) +
            scale_color_manual(values = c("FALSE" = "darkgray", "TRUE" = "red")) +
            labs(
              title = paste("MASC Odds Ratios for", cell_type),
              subtitle = paste("Red: FDR <", fdr_threshold),
              x = "Log Odds Ratio (95% CI)",
              y = "Comparison"
            ) +
            theme_bw() +
            theme(
              legend.position = "none",
              plot.title = element_text(face = "bold", size = 14)
            )
          
          # 强调显著结果行 (修复向量化问题)
          if (any(ct_results$Significant_FDR)) {
            p_forest <- p_forest + 
              theme(axis.text.y = element_text(face = "bold", 
                                               color = "red", 
                                               size = 10))
          }
          
          ggsave(file.path(output_dir, paste0(cell_type, "_forest_plot.pdf")), p_forest, width = 10, height = 8)
          message(paste("  Created forest plot for", cell_type))
        }
      }
    }
    
    # 检查并确保安装了必要的包
    required_vis_packages <- c("reshape2", "pheatmap")
    for (pkg in required_vis_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste(pkg, "package is required. Please install it with: install.packages('", pkg, "')"))
        if (pkg == "pheatmap") {
          run_heatmap <- FALSE 
        }
      }
    }
    
    # 创建热图显示所有细胞类型和比较
    if (requireNamespace("pheatmap", quietly = TRUE) && nrow(all_results) > 0) {
      message("Generating heatmap of log odds ratios across comparisons...")
      
      tryCatch({
        # 准备热图数据
        # 为缺失的OR值添加占位符
        all_results$log_OR <- log(all_results$OR)
        all_results$log_OR[!is.finite(all_results$log_OR)] <- 0
        
        # 打印调试信息
        message(paste("Total results rows:", nrow(all_results)))
        message(paste("Unique cell types:", length(unique(all_results$CellType))))
        message(paste("Unique comparisons:", length(unique(all_results$Comparison))))
        
        # 过滤掉NA值
        heatmap_data <- all_results[!is.na(all_results$CellType) & 
                                      !is.na(all_results$Comparison) &
                                      !is.na(all_results$log_OR) &
                                      !is.na(all_results$Significant_FDR), ]
        
        # 更可靠的宽格式转换 (使用reshape2)
        if (requireNamespace("reshape2", quietly = TRUE)) {
          library(reshape2)
          
          # 创建日志OR矩阵
          heat_matrix_df <- reshape2::dcast(heatmap_data, 
                                            CellType ~ Comparison, 
                                            value.var = "log_OR", 
                                            fill = 0)
          
          # 创建显著性矩阵
          sig_matrix_df <- reshape2::dcast(heatmap_data, 
                                           CellType ~ Comparison, 
                                           value.var = "Significant_FDR", 
                                           fill = FALSE)
          
          message(paste("Heat matrix dimensions:", nrow(heat_matrix_df), "rows x", 
                        ncol(heat_matrix_df), "columns"))
          
          # 检查矩阵是否有足够的列
          if (ncol(heat_matrix_df) > 1) {  # 至少需要CellType列和一个比较列
            # 转换为矩阵格式
            rownames(heat_matrix_df) <- heat_matrix_df$CellType
            heat_matrix <- as.matrix(heat_matrix_df[, -1, drop = FALSE])
            
            rownames(sig_matrix_df) <- sig_matrix_df$CellType
            sig_matrix <- as.matrix(sig_matrix_df[, -1, drop = FALSE])
            
            # 进行直观截断
            max_abs_log_or <- max(abs(heat_matrix), na.rm = TRUE)
            cap_value <- min(max_abs_log_or, 3)  # 截断极值为±3
            heat_matrix_capped <- pmin(pmax(heat_matrix, -cap_value), cap_value)
            
            # 创建显著性标记矩阵
            sig_symbols <- matrix("", 
                                  nrow = nrow(heat_matrix_capped), 
                                  ncol = ncol(heat_matrix_capped),
                                  dimnames = dimnames(heat_matrix_capped))
            
            # 填充显著性标记
            for(i in 1:nrow(sig_matrix)) {
              for(j in 1:ncol(sig_matrix)) {
                if(!is.na(sig_matrix[i, j]) && sig_matrix[i, j]) {
                  sig_symbols[i, j] <- "*"
                }
              }
            }
            
            # 提取组织信息
            comparison_split <- strsplit(colnames(heat_matrix_capped), " vs ")
            tissue1 <- sapply(comparison_split, function(x) if(length(x) > 1) x[2] else NA)
            tissue2 <- sapply(comparison_split, function(x) if(length(x) > 0) x[1] else NA)
            
            # 提取唯一组织名称
            unique_tissues <- unique(c(tissue1, tissue2))
            unique_tissues <- unique_tissues[!is.na(unique_tissues)]
            
            message(paste("Unique tissues for coloring:", 
                          paste(unique_tissues, collapse=", ")))
            
            # 创建列注释
            col_anno <- data.frame(
              Tissue1 = tissue1,
              Tissue2 = tissue2,
              row.names = colnames(heat_matrix_capped)
            )
            
            # 生成组织的颜色映射
            tissue_colors <- setNames(
              rainbow(length(unique_tissues)),
              unique_tissues
            )
            
            anno_colors <- list(
              Tissue1 = tissue_colors,
              Tissue2 = tissue_colors
            )
            
            # 创建热图
            pdf(file.path(output_dir, "pairwise_log_OR_heatmap.pdf"), width = 12, height = 10)
            pheatmap::pheatmap(
              heat_matrix_capped,
              main = "Log Odds Ratios Across All Comparisons",
              annotation_col = col_anno,
              annotation_colors = anno_colors,
              display_numbers = sig_symbols,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              breaks = seq(-cap_value, cap_value, length.out = 101)
            )
            dev.off()
            
            message("Heatmap created successfully")
          } else {
            message("Insufficient data for heatmap (need multiple comparisons)")
          }
        } else {
          message("reshape2 package is required for heatmap generation.")
        }
      }, error = function(e) {
        warning(paste("Error generating heatmap:", e$message))
        print(traceback())  # 打印完整堆栈跟踪
        message("Continuing with other visualizations...")
      })
    }
    
    # 为细胞类型在各组织中的丰度创建堆叠柱状图
    p_stack <- props_by_contrast %>%
      ggplot(aes(x = Contrast, y = Proportion, fill = CellType)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_bw() +
      labs(
        title = "Cell Type Composition Across Tissues",
        y = "Proportion",
        x = "Tissue"
      ) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = color_palette)
    
    ggsave(file.path(output_dir, "cell_composition_horizontal.pdf"), p_stack, width = 10, height = 8)
    
    # 创建多面板图，显示显著差异的细胞类型比例
    if (nrow(sig_results) > 0) {
      # 获取显著差异的细胞类型
      sig_celltypes <- unique(sig_results$CellType)
      
      # 获取所有组织
      all_tissues <- unique(props_df$Contrast)
      n_tissues <- length(all_tissues)
      
      message(paste("\nGenerating boxplots with adjusted p-values for", 
                    length(sig_celltypes), "significant cell types..."))
      
      # 生成所有可能的组织对比较 (只生成一次)
      global_comparisons <- list()
      global_comparison_names <- character()
      
      # 系统地生成所有可能的组织对
      for (i in 1:(n_tissues-1)) {
        for (j in (i+1):n_tissues) {
          tissue1 <- all_tissues[i]
          tissue2 <- all_tissues[j]
          # 确保比较名称与MASC结果中的格式一致
          comparison_name <- paste(tissue2, "vs", tissue1)
          global_comparisons[[length(global_comparisons) + 1]] <- c(tissue1, tissue2)
          global_comparison_names <- c(global_comparison_names, comparison_name)
        }
      }
      
      message(paste("Generated", length(global_comparison_names), 
                    "possible tissue comparisons"))
      
      # 为每个显著的细胞类型绘制箱线图
      for (cell_type in sig_celltypes) {
        message(paste("  Creating boxplot for", cell_type))
        
        # 获取该细胞类型的比例数据
        ct_data <- props_df[props_df$CellType == cell_type, ]
        
        # 获取该细胞类型的显著对比
        ct_sig_results <- sig_results[sig_results$CellType == cell_type, ]
        
        # 获取该细胞类型的所有比较结果
        ct_all_results <- all_results[all_results$CellType == cell_type, ]
        
        message(paste("  Found", nrow(ct_sig_results), "significant comparisons"))
        message(paste("  Found", nrow(ct_all_results), "total comparisons for this cell type"))
        
        # 确保Contrast是因子类型，保持固定顺序
        ct_data$Contrast <- factor(ct_data$Contrast)
        local_tissues <- levels(ct_data$Contrast)
        
        # 创建箱线图
        p <- ggplot(ct_data, aes(x = Contrast, y = Proportion, fill = Contrast)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) +
          geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
          theme_bw() +
          labs(
            title = paste0(cell_type, " - Tissue Comparisons"),
            subtitle = paste("Number of significant comparisons:", nrow(ct_sig_results)),
            y = "Proportion",
            x = "Tissue"
          ) +
          scale_y_continuous(labels = scales::percent) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"),
            legend.position = "none"
          )
        
        # 计算y轴比例数据范围，用于放置比较标签
        y_max <- max(ct_data$Proportion, na.rm = TRUE)
        y_range <- y_max * 0.15  # 标签间隔
        
        # 为该细胞类型创建本地比较列表
        local_comparisons <- list()
        local_comparison_names <- character()
        
        # 只为当前细胞类型的数据中实际存在的组织创建比较
        for (i in 1:(length(local_tissues)-1)) {
          for (j in (i+1):length(local_tissues)) {
            tissue1 <- local_tissues[i]
            tissue2 <- local_tissues[j]
            # 确保比较名称与MASC结果中的格式一致
            comparison_name <- paste(tissue2, "vs", tissue1)
            local_comparisons[[length(local_comparisons) + 1]] <- c(tissue1, tissue2)
            local_comparison_names <- c(local_comparison_names, comparison_name)
          }
        }
        
        message(paste("  Generated", length(local_comparison_names), 
                      "comparisons for local tissues:", 
                      paste(local_tissues, collapse=", ")))
        
        # 手动添加比较线和标签
        added_comparisons <- 0
        
        # 遍历本地比较
        for (i in 1:length(local_comparison_names)) {
          comp_name <- local_comparison_names[i]
          result_rows <- ct_all_results[ct_all_results$Comparison == comp_name, ]
          
          if (nrow(result_rows) > 0) {
            added_comparisons <- added_comparisons + 1
            pair <- local_comparisons[[i]]
            
            # 提取调整后的p值
            adj_p_value <- result_rows$FDR[1]
            original_p <- result_rows$model.pvalue[1]
            
            # 格式化调整后p值显示
            if (adj_p_value < 0.001) {
              p_text <- "adj.P < 0.001 ***"
            } else if (adj_p_value < 0.01) {
              p_text <- paste0("adj.P = ", sprintf("%.3f", adj_p_value), " **")
            } else if (adj_p_value < 0.05) {
              p_text <- paste0("adj.P = ", sprintf("%.3f", adj_p_value), " *")
            } else if (adj_p_value < 0.1) {
              p_text <- paste0("adj.P = ", sprintf("%.3f", adj_p_value), " .")
            } else {
              p_text <- paste0("adj.P = ", sprintf("%.3f", adj_p_value), " ns")
            }
            
            # 计算位置，确保不重叠
            row_num <- ceiling(added_comparisons / 2)  # 每行最多2个比较
            y_pos <- y_max + (row_num * y_range)
            
            # 获取组织的位置索引
            x1 <- which(local_tissues == as.character(pair[1]))
            x2 <- which(local_tissues == as.character(pair[2]))
            
            message(paste("    Adding comparison:", comp_name, 
                          "(p =", round(original_p, 4), 
                          ", adj.p =", round(adj_p_value, 4), ")"))
            
            if (length(x1) > 0 && length(x2) > 0) {
              # 添加比较线和标签
              p <- p + 
                annotate("segment", x = x1, y = y_pos, xend = x2, yend = y_pos) +
                annotate("segment", x = x1, y = y_pos - 0.05*y_max, xend = x1, yend = y_pos) +
                annotate("segment", x = x2, y = y_pos - 0.05*y_max, xend = x2, yend = y_pos) +
                annotate("text", x = (x1+x2)/2, y = y_pos + 0.02*y_max, 
                         label = p_text, size = 3)
            }
          } else {
            message(paste("    Skipping comparison:", comp_name, "- no results found"))
          }
        }
        
        # 调整Y轴范围，确保所有标签可见
        max_y_needed <- y_max + (ceiling(max(added_comparisons, 1)/2) * y_range) + 0.1*y_max
        p <- p + coord_cartesian(
          ylim = c(0, min(max_y_needed, 1.0)),  # 最大不超过100%
          expand = TRUE
        )
        
        # 保存图像
        safe_celltype <- gsub(" ", "_", cell_type)  # 替换空格为下划线
        ggsave(file.path(output_dir, paste0(safe_celltype, "_significant_boxplot.pdf")), 
               p, width = 12, height = 8)
        
        if (added_comparisons > 0) {
          message(paste("  Successfully added", added_comparisons, "comparisons to boxplot"))
        } else {
          message("  WARNING: No comparisons were added to the boxplot!")
        }
      }
      
      message(paste("Boxplots with adjusted p-values saved to", output_dir))
    }
  } else {
    warning("No results available for visualization")
  }
  
  # 记录分析结束时间
  end_time <- Sys.time()
  message(paste("\nPairwise MASC analysis completed, duration:", format(end_time - start_time)))
  
  # 返回结果
  return(list(
    counts = counts_matrix,
    proportions = props_df,
    all_results = all_results,
    pairwise_results = pairwise_results,
    significant_results = if(exists("sig_results")) sig_results else data.frame()
  ))
}

# 此函数可以替换scPairwiseMASCAnalysis函数中的箱式图部分，也可以单独运行
generate_pvalue_boxplots <- function(props_df, all_results, output_dir, 
                                     fdr_threshold = 0.05,
                                     show_all_celltypes = TRUE) {
  
  # 输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 标记显著结果（使用指定的阈值）
  all_results$Significant_FDR <- all_results$FDR < fdr_threshold
  
  # 仅保留显著的结果用于筛选
  sig_results <- all_results[all_results$Significant_FDR, ]
  
  # 确定要绘制的细胞类型
  if (show_all_celltypes) {
    # 显示所有细胞类型
    celltypes_to_plot <- unique(props_df$CellType)
  } else {
    # 只显示有显著差异的细胞类型
    celltypes_to_plot <- unique(sig_results$CellType)
  }
  
  if (length(celltypes_to_plot) == 0) {
    message("警告: 没有找到细胞类型或没有显著差异的细胞类型")
    return(NULL)
  }
  
  message(paste("将为", length(celltypes_to_plot), "个细胞类型生成箱式图"))
  
  # 获取所有组织
  all_tissues <- unique(props_df$Contrast)
  message(paste("检测到", length(all_tissues), "个组织:", paste(all_tissues, collapse=", ")))
  
  # 生成所有可能的组织对比较
  tissue_pairs <- list()
  comparison_names <- character()
  
  # 生成所有可能的组织对
  for (i in 1:(length(all_tissues)-1)) {
    for (j in (i+1):length(all_tissues)) {
      tissue1 <- as.character(all_tissues[i])
      tissue2 <- as.character(all_tissues[j])
      # 确保比较名称与MASC结果中的格式一致
      comparison_name <- paste(tissue2, "vs", tissue1)
      tissue_pairs[[length(tissue_pairs) + 1]] <- c(tissue1, tissue2)
      comparison_names <- c(comparison_names, comparison_name)
    }
  }
  
  # 检查生成的组织对
  message(paste("生成了", length(tissue_pairs), "个组织对比较:"))
  for (i in 1:length(tissue_pairs)) {
    message(paste("  ", i, ":", comparison_names[i]))
  }
  
  # 检查ggpubr包
  has_ggpubr <- requireNamespace("ggpubr", quietly = TRUE)
  if (!has_ggpubr) {
    message("警告: 建议安装ggpubr包以获得更好的图表效果 (install.packages('ggpubr'))")
  }
  
  # 为每个细胞类型创建箱式图
  for (cell_type in celltypes_to_plot) {
    # 获取该细胞类型的比例数据
    ct_data <- props_df[props_df$CellType == cell_type, ]
    
    # 获取该细胞类型的结果
    ct_results <- all_results[all_results$CellType == cell_type, ]
    
    # 获取显著的比较（用于标题）
    ct_sig_results <- sig_results[sig_results$CellType == cell_type, ]
    
    # 创建基础箱线图
    p <- ggplot(ct_data, aes(x = Contrast, y = Proportion, fill = Contrast)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
      theme_bw() +
      labs(
        title = paste0(cell_type, " - 组织间比例差异"),
        subtitle = if(nrow(ct_sig_results) > 0) {
          paste("显著比较 (FDR <", fdr_threshold, "):", nrow(ct_sig_results), "个")
        } else {
          paste("无显著差异 (FDR <", fdr_threshold, ")")
        },
        y = "比例",
        x = "组织"
      ) +
      scale_y_continuous(labels = scales::percent) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "none"
      )
    
    # 对箱式图添加所有组织对的比较
    # 计算数据的y轴范围
    y_max <- max(ct_data$Proportion, na.rm = TRUE) * 1.1
    y_range <- y_max * 0.15  # 用于分隔不同比较的间距
    
    # 根据组织对数量合理分配比较标签
    n_pairs <- length(tissue_pairs)
    max_pairs_per_row <- min(5, n_pairs)
    
    # 调试输出
    message(paste("处理细胞类型:", cell_type))
    message(paste("  组织对比较数量:", n_pairs))
    
    # 收集比较信息用于调试
    found_comparisons <- character()
    
    if (has_ggpubr) {
      # 使用ggpubr添加比较线和p值
      require(ggpubr)
      
      # 遍历所有组织对并添加比较
      for (i in 1:length(tissue_pairs)) {
        # 获取当前组织对
        pair <- tissue_pairs[[i]]
        comp_name <- comparison_names[i]
        
        # 从结果中查找对应的行
        result_row <- ct_results[ct_results$Comparison == comp_name, ]
        
        # 收集调试信息
        if (nrow(result_row) > 0) {
          found_comparisons <- c(found_comparisons, comp_name)
        } else {
          message(paste("  ⚠️ 未找到比较:", comp_name))
        }
        
        # 无论是否找到结果，都添加比较
        if (nrow(result_row) > 0) {
          # 提取FDR值和显著性
          fdr_value <- result_row$FDR[1]
          is_significant <- result_row$Significant_FDR[1]
          
          # 格式化FDR值文本
          if (fdr_value < 0.001) {
            fdr_text <- "FDR < 0.001"
          } else {
            fdr_text <- sprintf("FDR = %.3f", fdr_value)
          }
          
          # 添加显著性标记
          if (is_significant) {
            if (fdr_value < 0.001) sig_symbol <- "***"
            else if (fdr_value < 0.01) sig_symbol <- "**"
            else if (fdr_value < 0.05) sig_symbol <- "*"
            else sig_symbol <- ""
          } else {
            sig_symbol <- "ns"
          }
          
          # 完整标签文本
          label_text <- paste(fdr_text, sig_symbol)
          
          # 计算标签位置
          row_index <- ceiling(i / max_pairs_per_row)
          y_pos <- y_max + (row_index - 1) * y_range
          
          # 在ggplot图上添加比较
          comparison_list <- list(c(pair[1], pair[2]))
          p <- p + stat_compare_means(
            comparisons = comparison_list,
            label = label_text,
            method = "t.test",
            label.y = y_pos,
            size = 3,
            vjust = -0.5
          )
        }
      }
    } else {
      # 手动添加比较线和p值
      for (i in 1:length(tissue_pairs)) {
        # 获取当前组织对
        pair <- tissue_pairs[[i]]
        comp_name <- comparison_names[i]
        
        # 从结果中查找对应的行
        result_row <- ct_results[ct_results$Comparison == comp_name, ]
        
        # 收集调试信息
        if (nrow(result_row) > 0) {
          found_comparisons <- c(found_comparisons, comp_name)
        } else {
          message(paste("  ⚠️ 未找到比较:", comp_name))
        }
        
        # 无论是否找到结果，都添加比较
        # 只有在存在结果时添加比较
        if (nrow(result_row) > 0) {
          # 提取FDR值和显著性
          fdr_value <- result_row$FDR[1]
          is_significant <- result_row$Significant_FDR[1]
          
          # 格式化FDR值文本
          if (fdr_value < 0.001) {
            fdr_text <- "FDR < 0.001"
          } else {
            fdr_text <- sprintf("FDR = %.3f", fdr_value)
          }
          
          # 添加显著性标记
          if (is_significant) {
            if (fdr_value < 0.001) sig_symbol <- "***"
            else if (fdr_value < 0.01) sig_symbol <- "**"
            else if (fdr_value < 0.05) sig_symbol <- "*"
            else sig_symbol <- ""
          } else {
            sig_symbol <- "ns"
          }
          
          # 完整标签文本
          label_text <- paste(fdr_text, sig_symbol)
          
          # 根据组织名称计算位置
          all_tissue_levels <- levels(factor(all_tissues))
          if (length(all_tissue_levels) == 0) {
            all_tissue_levels <- unique(all_tissues)
          }
          
          x1 <- which(all_tissue_levels == pair[1])
          x2 <- which(all_tissue_levels == pair[2])
          
          # 计算标签位置
          row_index <- ceiling(i / max_pairs_per_row)
          y_pos <- y_max + (row_index - 1) * y_range
          y_line <- y_pos - 0.02 * y_max
          
          # 添加比较线和标签
          p <- p + 
            geom_segment(aes(x = x1, y = y_line, xend = x2, yend = y_line), 
                         inherit.aes = FALSE) +
            geom_segment(aes(x = x1, y = y_line - 0.02*y_max, xend = x1, yend = y_line), 
                         inherit.aes = FALSE) +
            geom_segment(aes(x = x2, y = y_line - 0.02*y_max, xend = x2, yend = y_line), 
                         inherit.aes = FALSE) +
            geom_text(aes(x = (x1+x2)/2, y = y_line + 0.02*y_max, label = label_text), 
                      inherit.aes = FALSE, size = 3)
        }
      }
    }
    
    # 调试输出已找到的比较
    message(paste("  找到的比较:", paste(found_comparisons, collapse=", ")))
    
    # 检查已知对比结果中是否包含所有期望的组织对
    for (comp_name in comparison_names) {
      if (!comp_name %in% ct_results$Comparison) {
        message(paste("  ⚠️ 结果中缺少比较:", comp_name))
      }
    }
    
    # 调整Y轴范围
    needed_height <- y_max + (ceiling(length(tissue_pairs) / max_pairs_per_row) * y_range)
    p <- p + coord_cartesian(
      ylim = c(0, min(needed_height, y_max * 3)),
      expand = TRUE
    )
    
    # 保存图像（使用安全文件名）
    safe_celltype <- gsub(" ", "_", cell_type)
    output_file <- file.path(output_dir, paste0(safe_celltype, "_pvalue_boxplot.pdf"))
    
    message(paste("  保存图表:", output_file))
    ggsave(output_file, p, width = 12, height = 8)
  }
  
  message(paste("所有箱式图已保存到:", output_dir))
}

# MASC收敛问题诊断与修复函数

#' 诊断MASC分析中的收敛问题
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名
#' @param sample_col 样本列名  
#' @param contrast_col 对比列名
#' @param min_cells_per_group 每组最小细胞数
#' @param min_samples_per_group 每组最小样本数
diagnose_masc_convergence <- function(seurat_obj, 
                                      cell_type_col = "Annotation",
                                      sample_col = "sample", 
                                      contrast_col = "tissue",
                                      min_cells_per_group = 20,
                                      min_samples_per_group = 3) {
  
  meta_data <- seurat_obj@meta.data
  
  cat("=== MASC收敛问题诊断报告 ===\n\n")
  
  # 1. 基础数据统计
  cat("1. 基础数据统计:\n")
  cat(sprintf("总细胞数: %d\n", nrow(meta_data)))
  cat(sprintf("总样本数: %d\n", length(unique(meta_data[[sample_col]]))))
  cat(sprintf("细胞类型数: %d\n", length(unique(meta_data[[cell_type_col]]))))
  cat(sprintf("组织/条件数: %d\n", length(unique(meta_data[[contrast_col]]))))
  
  # 2. 检查样本分布平衡性
  cat("\n2. 样本分布检查:\n")
  sample_by_tissue <- table(unique(meta_data[, c(sample_col, contrast_col)])[[contrast_col]])
  print(sample_by_tissue)
  
  unbalanced_tissues <- names(sample_by_tissue)[sample_by_tissue < min_samples_per_group]
  if(length(unbalanced_tissues) > 0) {
    cat("⚠️  样本数不足的组织:", paste(unbalanced_tissues, collapse = ", "), "\n")
  }
  
  # 3. 识别问题细胞类型
  cat("\n3. 细胞类型分布检查:\n")
  problematic_celltypes <- c()
  
  for(ct in unique(meta_data[[cell_type_col]])) {
    ct_data <- meta_data[meta_data[[cell_type_col]] == ct, ]
    
    # 检查每个组织中该细胞类型的计数
    ct_by_tissue <- table(ct_data[[contrast_col]])
    
    # 检查是否存在计数过少的组织
    low_count_tissues <- names(ct_by_tissue)[ct_by_tissue < min_cells_per_group]
    
    if(length(low_count_tissues) > 0) {
      problematic_celltypes <- c(problematic_celltypes, ct)
      cat(sprintf("⚠️  %s: 以下组织细胞数过少 (%s)\n", 
                  ct, paste(paste(low_count_tissues, ct_by_tissue[low_count_tissues], sep="="), 
                            collapse=", ")))
    }
  }
  
  # 4. 计算细胞类型的变异系数
  cat("\n4. 细胞类型分布变异性分析:\n")
  cv_results <- c()
  
  for(ct in unique(meta_data[[cell_type_col]])) {
    ct_counts <- table(meta_data[meta_data[[cell_type_col]] == ct, contrast_col])
    if(length(ct_counts) > 1) {
      cv <- sd(ct_counts) / mean(ct_counts)
      cv_results <- c(cv_results, cv)
      names(cv_results)[length(cv_results)] <- ct
    }
  }
  
  # 排序显示变异系数最高的细胞类型
  cv_sorted <- sort(cv_results, decreasing = TRUE)
  cat("变异系数最高的细胞类型（可能导致收敛问题）:\n")
  print(head(cv_sorted, 10))
  
  # 5. 生成建议
  cat("\n=== 修复建议 ===\n")
  
  if(length(problematic_celltypes) > 0) {
    cat("1. 考虑过滤以下细胞类型:\n")
    cat("   ", paste(problematic_celltypes, collapse = ", "), "\n")
  }
  
  if(length(unbalanced_tissues) > 0) {
    cat("2. 考虑排除样本数不足的组织:\n")
    cat("   ", paste(unbalanced_tissues, collapse = ", "), "\n")
  }
  
  cat("3. 建议的过滤参数:\n")
  cat(sprintf("   min_cells = %d (当前建议)\n", min_cells_per_group))
  cat(sprintf("   min_samples = %d (当前建议)\n", min_samples_per_group))
  
  high_cv_celltypes <- names(cv_sorted)[cv_sorted > 2.0]
  if(length(high_cv_celltypes) > 0) {
    cat("4. 高变异性细胞类型（考虑单独分析）:\n")
    cat("   ", paste(high_cv_celltypes, collapse = ", "), "\n")
  }
  
  return(list(
    problematic_celltypes = problematic_celltypes,
    unbalanced_tissues = unbalanced_tissues,
    cv_results = cv_sorted,
    recommendations = list(
      min_cells = min_cells_per_group,
      min_samples = min_samples_per_group,
      exclude_celltypes = problematic_celltypes,
      exclude_tissues = unbalanced_tissues
    )
  ))
}

#' 改进的MASC函数，增强收敛稳定性
#' @param dataset 数据集
#' @param cluster 细胞类型
#' @param contrast 对比变量
#' @param random_effects 随机效应
#' @param fixed_effects 固定效应
#' @param verbose 是否输出详细信息
#' @param max_iterations 最大迭代次数
#' @param tolerance 收敛容忍度
improved_MASC <- function(dataset, cluster, contrast, 
                          random_effects = NULL, fixed_effects = NULL,
                          verbose = FALSE, max_iterations = 50000,
                          tolerance = 0.01) {
  
  # 检查输入
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }
  
  # 将cluster转换为character
  cluster <- as.character(cluster)
  # 创建设计矩阵
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  
  # 创建输出列表
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # 创建模型公式
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
  } else {
    stop("No random effects specified")
  }
  
  if (verbose) {
    message(paste("Using model:", "cluster ~", model_rhs))
  }
  
  # 初始化列表存储每个簇的模型对象
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # 定义多个优化器选项
  optimizers <- list(
    list(optimizer = "bobyqa", 
         optCtrl = list(maxfun = max_iterations)),
    list(optimizer = "Nelder_Mead", 
         optCtrl = list(maxfun = max_iterations)),
    list(optimizer = c("bobyqa", "Nelder_Mead"), 
         optCtrl = list(maxfun = max_iterations))
  )
  
  # 为每个簇运行嵌套混合效应模型
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    
    if (verbose) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    
    # 数据质量检查
    cluster_by_contrast <- table(dataset[[contrast]], dataset[[test_cluster]])
    
    if (nrow(cluster_by_contrast) < 2 || ncol(cluster_by_contrast) < 2) {
      warning(paste("Invalid data structure for cluster", test_cluster, "- skipping"))
      cluster_models[[i]]$error <- "Invalid data structure"
      next
    }
    
    # 检查每个对比水平中是否有足够的变异
    has_enough_data <- TRUE
    for (level in levels(dataset[[contrast]])) {
      contrast_subset <- dataset[dataset[[contrast]] == level, ]
      cluster_count <- sum(contrast_subset[[test_cluster]])
      
      if (cluster_count < 2 || cluster_count >= (nrow(contrast_subset) - 1)) {
        has_enough_data <- FALSE
        break
      }
    }
    
    if (!has_enough_data) {
      warning(paste("Insufficient variation for cluster", test_cluster, "- skipping"))
      cluster_models[[i]]$error <- "Insufficient variation"
      next
    }
    
    # 构建模型公式
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    
    # 尝试不同的优化器
    model_success <- FALSE
    
    for (opt_config in optimizers) {
      if (model_success) break
      
      tryCatch({
        # 自定义收敛控制
        ctrl <- glmerControl(
          optimizer = opt_config$optimizer,
          optCtrl = opt_config$optCtrl,
          calc.derivs = FALSE,
          check.conv.grad = .makeCC("warning", tol = tolerance),
          check.conv.singular = .makeCC("ignore"),
          check.conv.hess = .makeCC("ignore")
        )
        
        null_model <- glmer(formula = null_fm, data = dataset,
                            family = binomial, nAGQ = 1, 
                            control = ctrl, verbose = 0)
        
        full_model <- glmer(formula = full_fm, data = dataset,
                            family = binomial, nAGQ = 1, 
                            control = ctrl, verbose = 0)
        
        # 检查模型收敛
        null_converged <- is.null(null_model@optinfo$conv$lme4$messages)
        full_converged <- is.null(full_model@optinfo$conv$lme4$messages)
        
        if (null_converged && full_converged) {
          model_lrt <- anova(null_model, full_model)
          
          # 计算置信区间
          contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
          contrast_ci <- confint.merMod(full_model, method = "Wald",
                                        parm = contrast_lvl2)
          
          # 保存成功的模型
          cluster_models[[i]]$null_model <- null_model
          cluster_models[[i]]$full_model <- full_model
          cluster_models[[i]]$model_lrt <- model_lrt
          cluster_models[[i]]$confint <- contrast_ci
          cluster_models[[i]]$optimizer_used <- opt_config$optimizer
          
          model_success <- TRUE
          
          if (verbose) {
            message(paste("Model converged successfully using", 
                          paste(opt_config$optimizer, collapse = "+")))
          }
        }
        
      }, error = function(e) {
        if (verbose) {
          message(paste("Optimizer", paste(opt_config$optimizer, collapse = "+"), 
                        "failed:", e$message))
        }
      })
    }
    
    if (!model_success) {
      warning(paste("All optimizers failed for cluster", test_cluster))
      cluster_models[[i]]$error <- "All optimizers failed"
    }
  }
  
  # 组织结果
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  
  # 检查成功的模型
  successful_models <- sapply(cluster_models, function(x) !is.null(x$model_lrt))
  
  if (any(successful_models)) {
    output$model.pvalue <- NA
    output$model.pvalue[successful_models] <- sapply(
      cluster_models[successful_models], 
      function(x) x$model_lrt[["Pr(>Chisq)"]][2]
    )
    
    # 添加FDR校正
    output$FDR <- p.adjust(output$model.pvalue, method = "BH")
    
    # 添加优势比和置信区间
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    output[[paste(contrast_lvl2, "OR", sep = ".")]] <- NA
    output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- NA
    output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- NA
    
    for (i in which(successful_models)) {
      output[[paste(contrast_lvl2, "OR", sep = ".")]][i] <- 
        exp(fixef(cluster_models[[i]]$full_model)[[contrast_lvl2]])
      output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]][i] <- 
        exp(cluster_models[[i]]$confint[contrast_lvl2, "2.5 %"])
      output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]][i] <- 
        exp(cluster_models[[i]]$confint[contrast_lvl2, "97.5 %"])
    }
    
    # 添加使用的优化器信息
    output$optimizer_used <- NA
    output$optimizer_used[successful_models] <- sapply(
      cluster_models[successful_models],
      function(x) paste(x$optimizer_used, collapse = "+")
    )
  }
  
  return(list(results = output, models = cluster_models))
}



# 使用示例 / Usage Example:
# 
# # 1. 准备数据
# props_data <- prepare_data_for_plotting(
#   seurat_obj = your_seurat_object,
#   cell_type_col = "Annotation",
#   sample_col = "sample",
#   contrast_col = "tissue"
# )
# 
# # 2. 创建图表（确保与您的数据匹配）
# publication_plot <- create_publication_boxplot(
#   props_df = props_data,
#   stats_results = your_masc_results$all_results,  # 来自scMASC分析
#   output_dir = "publication_figures",
#   figure_width = 14,
#   figure_height = 10,
#   cell_types_order = c("Epithelial", "Endothelial", "Myeloid", "SMC", "T", "Fibroblast", "B"),
#   contrast_colors = c(
#     "lung parenchyma" = "#2E8B57",    # 深绿色
#     "sinus" = "#B22222",             # 砖红色  
#     "nose" = "#CD853F",              # 沙褐色
#     "respiratory airway" = "#4682B4"  # 钢蓝色
#   ),
#   show_individual_points = TRUE,
#   add_statistics = TRUE
# )

create_publication_boxplot <- function(props_df, 
                                       stats_results = NULL,
                                       output_dir = "publication_plots",
                                       figure_width = 12,
                                       figure_height = 8,
                                       cell_types_order = NULL,
                                       contrast_colors = NULL,
                                       show_individual_points = TRUE,
                                       font_size = 12,
                                       add_statistics = TRUE) {
  
  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 数据预处理 - Data preprocessing
  message("正在处理数据 / Processing data...")
  
  # 确保数据格式正确
  props_df$CellType <- as.character(props_df$CellType)
  props_df$Contrast <- as.character(props_df$Contrast)
  props_df$Proportion <- as.numeric(props_df$Proportion)
  
  # 设置细胞类型顺序
  if (is.null(cell_types_order)) {
    cell_types_order <- sort(unique(props_df$CellType))
  }
  props_df$CellType <- factor(props_df$CellType, levels = cell_types_order)
  
  # 设置对比组顺序
  contrast_levels <- unique(props_df$Contrast)
  props_df$Contrast <- factor(props_df$Contrast, levels = contrast_levels)
  
  # 设置颜色方案
  if (is.null(contrast_colors)) {
    # 使用临床研究常用的专业配色方案
    n_contrasts <- length(contrast_levels)
    if (n_contrasts <= 4) {
      contrast_colors <- c("#2E8B57", "#4682B4", "#CD853F", "#B22222")  # 深绿、钢蓝、沙褐、砖红
    } else {
      contrast_colors <- rainbow(n_contrasts, start = 0.1, end = 0.9)
    }
    names(contrast_colors) <- contrast_levels
  }
  
  # 主图绘制 - Main plot creation
  message("正在创建主图 / Creating main plot...")
  
  # 设置dodge宽度（这个值必须与add_significance_markers中保持一致）
  dodge_width <- 0.8
  
  p <- ggplot(props_df, aes(x = CellType, y = Proportion, fill = Contrast)) +
    geom_boxplot(position = position_dodge(width = dodge_width), 
                 outlier.shape = NA, 
                 alpha = 0.8,
                 width = 0.6) +
    
    # 添加个体数据点
    {if(show_individual_points) {
      geom_jitter(position = position_jitterdodge(dodge.width = dodge_width, jitter.width = 0.2),
                  size = 1.5, alpha = 0.6, stroke = 0.3)
    }} +
    
    # 主题设置
    theme_bw() +
    theme(
      # 轴文本设置
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                 size = font_size, color = "black"),
      axis.text.y = element_text(size = font_size, color = "black"),
      
      # 轴标题设置
      axis.title.x = element_text(size = font_size + 2, face = "bold", 
                                  margin = margin(t = 10)),
      axis.title.y = element_text(size = font_size + 2, face = "bold",
                                  margin = margin(r = 10)),
      
      # 图例设置
      legend.title = element_text(size = font_size + 1, face = "bold"),
      legend.text = element_text(size = font_size),
      legend.position = "top",
      legend.box = "horizontal",
      
      # 面板设置
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1),
      
      # 背景设置
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # 边距设置
      plot.margin = margin(20, 20, 20, 20)
    ) +
    
    # 标签和颜色设置
    labs(
      x = "Cell Type",
      y = "Fraction of Cell Number",
      fill = "Group"
    ) +
    scale_fill_manual(values = contrast_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.15)))  # 为统计标记留出空间
  
  # 添加统计显著性标记
  if (add_statistics && !is.null(stats_results)) {
    message("正在添加统计显著性标记 / Adding statistical significance markers...")
    
    # 为每个细胞类型添加统计标记
    p <- add_significance_markers(p, props_df, stats_results, contrast_levels)
  }
  
  # 保存高质量图片
  message("正在保存图片 / Saving plots...")
  
  # PDF格式（矢量图，适合期刊发表）
  pdf_file <- file.path(output_dir, "cell_proportions_publication.pdf")
  ggsave(pdf_file, p, width = figure_width, height = figure_height, 
         device = "pdf", dpi = 300, useDingbats = FALSE)
  
  # PNG格式（高分辨率位图，适合演示）
  png_file <- file.path(output_dir, "cell_proportions_publication.png")
  ggsave(png_file, p, width = figure_width, height = figure_height, 
         device = "png", dpi = 300, bg = "white")
  
  # EPS格式（矢量图，某些期刊要求）
  eps_file <- file.path(output_dir, "cell_proportions_publication.eps")
  ggsave(eps_file, p, width = figure_width, height = figure_height, 
         device = "eps", dpi = 300)
  
  message(paste("图片已保存到:", output_dir))
  message("Files saved: PDF, PNG, EPS formats")
  
  return(p)
}

#' 添加统计显著性标记的辅助函数
#' Helper function to add statistical significance markers
add_significance_markers <- function(plot_obj, props_df, stats_results, contrast_levels) {
  
  # 获取Y轴数据范围
  y_max <- max(props_df$Proportion, na.rm = TRUE)
  y_increment <- y_max * 0.12  # 增加标记间距
  
  # 为每个细胞类型处理统计结果
  cell_types <- levels(props_df$CellType)
  n_contrasts <- length(contrast_levels)
  
  # 计算dodge宽度（与ggplot中的position_dodge保持一致）
  dodge_width <- 0.8
  
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    # 获取该细胞类型的统计结果
    ct_stats <- stats_results[stats_results$CellType == cell_type, ]
    
    if (nrow(ct_stats) == 0) next
    
    # 解析比较信息并添加标记
    comparison_count <- 0
    
    for (j in 1:nrow(ct_stats)) {
      comparison <- ct_stats$Comparison[j]
      fdr_value <- ct_stats$FDR[j]
      
      # 解析比较的两个组（假设格式为 "GroupA vs GroupB"）
      groups <- trimws(strsplit(comparison, " vs ")[[1]])
      
      if (length(groups) == 2 && all(groups %in% contrast_levels)) {
        # 确定显著性符号
        if (fdr_value < 0.0001) {
          sig_symbol <- "****"
        } else if (fdr_value < 0.001) {
          sig_symbol <- "***"
        } else if (fdr_value < 0.01) {
          sig_symbol <- "**"
        } else if (fdr_value < 0.05) {
          sig_symbol <- "*"
        } else if (fdr_value < 0.1) {
          sig_symbol <- "."
        } else {
          next  # 不显著，跳过
        }
        
        # 计算标记位置
        comparison_count <- comparison_count + 1
        y_position <- y_max + (comparison_count * y_increment)
        
        # 精确计算箱式图的X轴位置
        # 每个细胞类型中心位置
        center_x <- i
        
        # 计算每个箱子相对中心的偏移
        # position_dodge会将箱子分布在中心周围
        group1_index <- which(contrast_levels == groups[1])
        group2_index <- which(contrast_levels == groups[2])
        
        # 计算偏移量（基于position_dodge的逻辑）
        offset1 <- (group1_index - (n_contrasts + 1) / 2) * (dodge_width / n_contrasts)
        offset2 <- (group2_index - (n_contrasts + 1) / 2) * (dodge_width / n_contrasts)
        
        group1_x <- center_x + offset1
        group2_x <- center_x + offset2
        
        # 添加连接线和显著性标记
        plot_obj <- plot_obj +
          annotate("segment", 
                   x = group1_x, xend = group2_x, 
                   y = y_position, yend = y_position,
                   color = "black", size = 0.5) +
          annotate("segment", 
                   x = group1_x, xend = group1_x, 
                   y = y_position - y_increment * 0.15, yend = y_position,
                   color = "black", size = 0.5) +
          annotate("segment", 
                   x = group2_x, xend = group2_x, 
                   y = y_position - y_increment * 0.15, yend = y_position,
                   color = "black", size = 0.5) +
          annotate("text", 
                   x = (group1_x + group2_x) / 2, 
                   y = y_position + y_increment * 0.2,
                   label = sig_symbol, 
                   size = 3.5, fontface = "bold")
      }
    }
  }
  
  return(plot_obj)
}

#' 从MASC结果中直接提取比例数据（推荐方法）
#' Extract proportion data directly from MASC results
extract_props_from_masc <- function(masc_results) {
  
  message("从MASC结果中提取比例数据...")
  
  # 检查MASC结果中是否有proportions数据
  if ("proportions" %in% names(masc_results)) {
    props_df <- masc_results$proportions
    message("成功从MASC结果中提取比例数据")
    
    # 检查数据格式
    required_cols <- c("CellType", "Sample", "Contrast", "Proportion")
    if (all(required_cols %in% colnames(props_df))) {
      message(paste("比例数据包含", nrow(props_df), "行"))
      message(paste("细胞类型:", length(unique(props_df$CellType)), "种"))
      message(paste("样本数:", length(unique(props_df$Sample)), "个"))
      message(paste("对比组:", length(unique(props_df$Contrast)), "个"))
      
      return(props_df)
    } else {
      warning("MASC比例数据格式不正确，缺少必要列")
      return(NULL)
    }
  } else {
    message("MASC结果中未找到proportions数据，尝试其他方法...")
    return(NULL)
  }
}

#' 与MASC完全一致的数据过滤和比例计算函数
#' MASC-consistent data filtering and proportion calculation
prepare_data_for_plotting_masc_consistent <- function(seurat_obj, 
                                                      cell_type_col = "Annotation",
                                                      sample_col = "sample", 
                                                      contrast_col = "tissue",
                                                      exclude_filter = NULL,
                                                      min_samples = 2,
                                                      min_cells = 10,
                                                      min_prop = 0.01) {
  
  message("使用与MASC一致的方法准备绘图数据...")
  
  # 获取元数据
  meta_data <- seurat_obj@meta.data
  
  # 检查必要的列
  required_cols <- c(cell_type_col, sample_col, contrast_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta_data)]
  if (length(missing_cols) > 0) {
    stop(paste("缺少以下列:", paste(missing_cols, collapse = ", ")))
  }
  
  message(paste("原始细胞数:", nrow(meta_data)))
  
  # 步骤1: 应用exclude_filter（与MASC完全相同的逻辑）
  filtered_meta <- meta_data
  original_count <- nrow(meta_data)
  
  if (!is.null(exclude_filter) && exclude_filter != "") {
    tryCatch({
      filter_expr <- parse(text = exclude_filter)
      rows_to_exclude <- eval(filter_expr, meta_data)
      filtered_meta <- meta_data[!rows_to_exclude, ]
      
      excluded_count <- original_count - nrow(filtered_meta)
      message(paste("应用过滤条件:", exclude_filter))
      message(paste("过滤前:", original_count, "个细胞"))
      message(paste("过滤掉:", excluded_count, "个细胞"))
      message(paste("过滤后:", nrow(filtered_meta), "个细胞"))
      
      if (nrow(filtered_meta) == 0) {
        warning("过滤条件排除了所有细胞！使用原始数据。")
        filtered_meta <- meta_data
      }
    }, error = function(e) {
      message(paste("过滤表达式错误:", e$message))
      message("使用所有细胞，不进行过滤")
      filtered_meta <- meta_data
    })
  }
  
  # 移除NA值
  na_cols <- c(cell_type_col, sample_col, contrast_col)
  has_na <- apply(filtered_meta[, na_cols, drop = FALSE], 1, function(x) any(is.na(x)))
  if (any(has_na)) {
    message(paste("移除", sum(has_na), "个含有NA值的细胞"))
    filtered_meta <- filtered_meta[!has_na, ]
  }
  
  # 转换为因子
  filtered_meta[[cell_type_col]] <- as.factor(filtered_meta[[cell_type_col]])
  filtered_meta[[sample_col]] <- as.factor(filtered_meta[[sample_col]])
  filtered_meta[[contrast_col]] <- as.factor(filtered_meta[[contrast_col]])
  
  # 步骤2: 检查组织样本数量并过滤不足的组织
  all_tissues <- levels(filtered_meta[[contrast_col]])
  n_tissues <- length(all_tissues)
  
  if (n_tissues < 2) {
    stop("需要至少2个组织进行比较")
  }
  
  # 统计每个组织的样本数量
  samples_by_tissue <- tapply(filtered_meta[[sample_col]], 
                              filtered_meta[[contrast_col]], 
                              function(x) length(unique(x)))
  
  message("每个组织的样本数:")
  for (tissue in names(samples_by_tissue)) {
    message(paste(" -", tissue, ":", samples_by_tissue[tissue], "个样本"))
  }
  
  # 过滤样本数量不足的组织
  tissues_with_enough_samples <- names(samples_by_tissue)[samples_by_tissue >= min_samples]
  if (length(tissues_with_enough_samples) < 2) {
    stop(paste("没有足够的组织具有", min_samples, "个或更多样本进行比较"))
  }
  
  if (length(tissues_with_enough_samples) < n_tissues) {
    excluded_tissues <- setdiff(names(samples_by_tissue), tissues_with_enough_samples)
    message(paste("以下组织样本数不足，将被排除:", paste(excluded_tissues, collapse = ", ")))
    
    filtered_meta <- filtered_meta[filtered_meta[[contrast_col]] %in% tissues_with_enough_samples, ]
    filtered_meta[[contrast_col]] <- factor(filtered_meta[[contrast_col]])
  }
  
  # 更新组织列表
  all_tissues <- levels(filtered_meta[[contrast_col]])
  
  # 步骤3: 计算细胞计数并过滤低计数的细胞类型
  message("计算细胞计数和比例...")
  
  counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
  counts_matrix <- as.matrix(counts_table)
  
  # 添加样本的对比信息
  sample_contrast_df <- unique(filtered_meta[, c(sample_col, contrast_col)])
  sample_to_contrast <- setNames(
    as.character(sample_contrast_df[[contrast_col]]), 
    as.character(sample_contrast_df[[sample_col]])
  )
  
  # 过滤总细胞数不足的细胞类型
  cell_type_totals <- rowSums(counts_matrix)
  low_count_cell_types <- names(cell_type_totals)[cell_type_totals < min_cells]
  
  if (length(low_count_cell_types) > 0) {
    message(paste("以下细胞类型总细胞数少于", min_cells, "个，将被排除:"))
    message(paste(low_count_cell_types, collapse = ", "))
    
    filtered_meta <- filtered_meta[!filtered_meta[[cell_type_col]] %in% low_count_cell_types, ]
    counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
    counts_matrix <- as.matrix(counts_table)
  }
  
  # 步骤4: 计算比例并过滤低比例的细胞类型
  if (nrow(counts_matrix) > 0 && ncol(counts_matrix) > 0) {
    props_by_sample <- prop.table(counts_matrix, margin = 2)
    
    # 过滤平均比例太低的细胞类型
    cell_type_mean_props <- rowMeans(props_by_sample)
    low_prop_cell_types <- names(cell_type_mean_props)[cell_type_mean_props < min_prop]
    
    if (length(low_prop_cell_types) > 0) {
      message(paste("以下细胞类型平均比例低于", min_prop, "，将被排除:"))
      message(paste(low_prop_cell_types, collapse = ", "))
      
      filtered_meta <- filtered_meta[!filtered_meta[[cell_type_col]] %in% low_prop_cell_types, ]
      counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
      counts_matrix <- as.matrix(counts_table)
      props_by_sample <- prop.table(counts_matrix, margin = 2)
    }
    
    # 转换为长格式
    props_df <- reshape2::melt(props_by_sample)
    colnames(props_df) <- c("CellType", "Sample", "Proportion")
    
    # 添加对比信息
    props_df$Contrast <- sample_to_contrast[as.character(props_df$Sample)]
    
    # 移除NA值
    if (any(is.na(props_df$Contrast))) {
      warning("某些样本无法匹配到对比组")
      props_df <- props_df[!is.na(props_df$Contrast), ]
    }
    
    # 转换数据类型
    props_df$CellType <- as.character(props_df$CellType)
    props_df$Sample <- as.character(props_df$Sample)
    props_df$Contrast <- as.character(props_df$Contrast)
    props_df$Proportion <- as.numeric(props_df$Proportion)
    
    message("最终数据准备完成:")
    message(paste("- 细胞类型:", length(unique(props_df$CellType)), "种"))
    message(paste("- 对比组:", length(unique(props_df$Contrast)), "个"))
    message(paste("- 样本数:", length(unique(props_df$Sample)), "个"))
    message(paste("- 总行数:", nrow(props_df)))
    
    return(props_df)
    
  } else {
    stop("过滤后没有数据剩余，请检查过滤条件")
  }
}

#' 智能数据准备函数（优先使用MASC结果）
#' Smart data preparation function (prioritize MASC results)
smart_prepare_plotting_data <- function(seurat_obj = NULL, 
                                        masc_results = NULL,
                                        cell_type_col = "Annotation",
                                        sample_col = "sample",
                                        contrast_col = "tissue",
                                        exclude_filter = NULL,
                                        min_samples = 2,
                                        min_cells = 10,
                                        min_prop = 0.01) {
  
  # 方法1: 优先尝试从MASC结果中提取
  if (!is.null(masc_results)) {
    message("尝试从MASC结果中提取比例数据...")
    props_df <- extract_props_from_masc(masc_results)
    
    if (!is.null(props_df)) {
      message("✓ 成功从MASC结果中获取比例数据")
      return(props_df)
    }
  }
  
  # 方法2: 使用与MASC一致的方法重新计算
  if (!is.null(seurat_obj)) {
    message("从Seurat对象重新计算比例数据（与MASC一致的方法）...")
    props_df <- prepare_data_for_plotting_masc_consistent(
      seurat_obj = seurat_obj,
      cell_type_col = cell_type_col,
      sample_col = sample_col,
      contrast_col = contrast_col,
      exclude_filter = exclude_filter,
      min_samples = min_samples,
      min_cells = min_cells,
      min_prop = min_prop
    )
    
    message("✓ 成功使用MASC一致的方法计算比例数据")
    return(props_df)
  }
  
  stop("需要提供MASC结果或Seurat对象")
}

#' 创建基础箱式图（不含显著性标记）
create_base_boxplot <- function(props_df, contrast_colors = NULL, 
                                cell_types_order = NULL, dodge_width = 0.8) {
  
  # 数据预处理
  if (!is.null(cell_types_order)) {
    props_df$CellType <- factor(props_df$CellType, levels = cell_types_order)
  }
  
  contrast_levels <- unique(props_df$Contrast)
  props_df$Contrast <- factor(props_df$Contrast, levels = contrast_levels)
  
  # 设置颜色
  if (is.null(contrast_colors)) {
    contrast_colors <- c("#2E8B57", "#B22222", "#CD853F", "#4682B4")[1:length(contrast_levels)]
    names(contrast_colors) <- contrast_levels
  }
  
  # 创建基础图
  p <- ggplot(props_df, aes(x = CellType, y = Proportion, fill = Contrast)) +
    geom_boxplot(position = position_dodge(width = dodge_width), 
                 outlier.shape = NA, 
                 alpha = 0.8,
                 width = 0.6) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge_width, jitter.width = 0.2),
                size = 1.5, alpha = 0.6, stroke = 0.3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                 size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "top",
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      x = "Cell Type",
      y = "Fraction of Cell Number",
      fill = "Group"
    ) +
    scale_fill_manual(values = contrast_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.30)))  # 为显著性标记预留空间
  
  return(list(plot = p, dodge_width = dodge_width, contrast_levels = contrast_levels))
}

#' 计算精确的箱式图位置
calculate_box_positions <- function(cell_types, contrast_levels, dodge_width = 0.8) {
  
  n_cell_types <- length(cell_types)
  n_contrasts <- length(contrast_levels)
  
  # 存储每个箱子的精确位置
  positions <- data.frame(
    CellType = character(0),
    Contrast = character(0),
    x_position = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    center_x <- i  # 细胞类型的中心位置
    
    # 计算每个对比组在该细胞类型下的位置
    if (n_contrasts == 1) {
      # 只有一个对比组，位置就是中心
      x_pos <- center_x
      positions <- rbind(positions, data.frame(
        CellType = cell_type,
        Contrast = contrast_levels[1],
        x_position = x_pos,
        stringsAsFactors = FALSE
      ))
    } else {
      # 多个对比组，均匀分布在中心周围
      for (j in seq_along(contrast_levels)) {
        contrast <- contrast_levels[j]
        
        # 使用与ggplot2相同的位置计算方法
        offset <- (j - (n_contrasts + 1) / 2) * (dodge_width / n_contrasts)
        x_pos <- center_x + offset
        
        positions <- rbind(positions, data.frame(
          CellType = cell_type,
          Contrast = contrast,
          x_position = x_pos,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(positions)
}

#' 创建显著性标记层
create_significance_layer <- function(stats_results, props_df, positions_df, 
                                      cell_types, contrast_levels) {
  
  # 获取Y轴范围
  y_max <- max(props_df$Proportion, na.rm = TRUE)
  y_increment <- y_max * 0.12
  
  # 存储所有注释
  annotations <- list()
  
  for (cell_type in cell_types) {
    # 获取该细胞类型的统计结果
    ct_stats <- stats_results[stats_results$CellType == cell_type, ]
    
    if (nrow(ct_stats) == 0) next
    
    comparison_count <- 0
    
    for (j in 1:nrow(ct_stats)) {
      comparison <- ct_stats$Comparison[j]
      fdr_value <- ct_stats$FDR[j]
      
      # 解析比较组
      groups <- trimws(strsplit(comparison, " vs ")[[1]])
      
      if (length(groups) == 2 && all(groups %in% contrast_levels)) {
        # 确定显著性符号
        if (fdr_value < 0.0001) {
          sig_symbol <- "****"
        } else if (fdr_value < 0.001) {
          sig_symbol <- "***"
        } else if (fdr_value < 0.01) {
          sig_symbol <- "**"
        } else if (fdr_value < 0.05) {
          sig_symbol <- "*"
        } else if (fdr_value < 0.1) {
          sig_symbol <- "."
        } else {
          next  # 不显著，跳过
        }
        
        comparison_count <- comparison_count + 1
        y_position <- y_max + (comparison_count * y_increment)
        
        # 从预计算的位置表中获取精确位置
        group1_pos <- positions_df[positions_df$CellType == cell_type & 
                                     positions_df$Contrast == groups[1], "x_position"]
        group2_pos <- positions_df[positions_df$CellType == cell_type & 
                                     positions_df$Contrast == groups[2], "x_position"]
        
        if (length(group1_pos) > 0 && length(group2_pos) > 0) {
          # 添加到注释列表
          annotations[[length(annotations) + 1]] <- list(
            type = "horizontal_line",
            x = group1_pos, xend = group2_pos,
            y = y_position, yend = y_position
          )
          
          annotations[[length(annotations) + 1]] <- list(
            type = "vertical_line_1",
            x = group1_pos, xend = group1_pos,
            y = y_position - y_increment * 0.15, yend = y_position
          )
          
          annotations[[length(annotations) + 1]] <- list(
            type = "vertical_line_2", 
            x = group2_pos, xend = group2_pos,
            y = y_position - y_increment * 0.15, yend = y_position
          )
          
          annotations[[length(annotations) + 1]] <- list(
            type = "text",
            x = (group1_pos + group2_pos) / 2,
            y = y_position + y_increment * 0.2,
            label = sig_symbol
          )
          
          message(paste("添加显著性标记:", cell_type, "-", comparison, 
                        "- 位置:", round(group1_pos, 2), "to", round(group2_pos, 2)))
        }
      }
    }
  }
  
  return(annotations)
}

#' 将注释添加到图上
add_annotations_to_plot <- function(base_plot, annotations) {
  
  p <- base_plot
  
  for (annotation in annotations) {
    if (annotation$type == "horizontal_line") {
      p <- p + annotate("segment", 
                        x = annotation$x, xend = annotation$xend,
                        y = annotation$y, yend = annotation$yend,
                        color = "black", size = 0.5)
    } else if (annotation$type %in% c("vertical_line_1", "vertical_line_2")) {
      p <- p + annotate("segment",
                        x = annotation$x, xend = annotation$xend,
                        y = annotation$y, yend = annotation$yend,
                        color = "black", size = 0.5)
    } else if (annotation$type == "text") {
      p <- p + annotate("text",
                        x = annotation$x, y = annotation$y,
                        label = annotation$label,
                        size = 3.5, fontface = "bold")
    }
  }
  
  return(p)
}

#' 主函数：分别生成并合并
create_separate_aligned_plot <- function(props_df, stats_results, 
                                         output_dir = "separate_aligned_plots",
                                         figure_width = 14, figure_height = 10,
                                         cell_types_order = NULL, contrast_colors = NULL,
                                         dodge_width = 0.8) {
  
  # 确保输出目录存在
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Step 1: 创建基础箱式图...")
  
  # 第一步：创建基础箱式图
  base_result <- create_base_boxplot(props_df, contrast_colors, cell_types_order, dodge_width)
  base_plot <- base_result$plot
  contrast_levels <- base_result$contrast_levels
  
  # 获取最终的细胞类型顺序
  if (is.null(cell_types_order)) {
    cell_types <- levels(props_df$CellType)
  } else {
    cell_types <- cell_types_order
  }
  
  # 保存基础图
  ggsave(file.path(output_dir, "base_boxplot.pdf"), base_plot, 
         width = figure_width, height = figure_height)
  
  message("Step 2: 计算精确位置...")
  
  # 第二步：计算所有箱子的精确位置
  positions_df <- calculate_box_positions(cell_types, contrast_levels, dodge_width)
  
  # 保存位置信息用于调试
  write.csv(positions_df, file.path(output_dir, "box_positions.csv"), row.names = FALSE)
  message("箱子位置信息已保存到 box_positions.csv")
  
  if (!is.null(stats_results) && nrow(stats_results) > 0) {
    
    message("Step 3: 创建显著性标记...")
    
    # 第三步：创建显著性标记
    annotations <- create_significance_layer(stats_results, props_df, positions_df, 
                                             cell_types, contrast_levels)
    
    message(paste("创建了", length(annotations), "个注释元素"))
    
    message("Step 4: 合并图层...")
    
    # 第四步：合并
    final_plot <- add_annotations_to_plot(base_plot, annotations)
    
  } else {
    message("没有统计结果，只显示基础箱式图")
    final_plot <- base_plot
  }
  
  # 保存最终图片
  pdf_file <- file.path(output_dir, "final_aligned_plot.pdf")
  ggsave(pdf_file, final_plot, width = figure_width, height = figure_height, 
         device = "pdf", dpi = 300, useDingbats = FALSE)
  
  png_file <- file.path(output_dir, "final_aligned_plot.png")
  ggsave(png_file, final_plot, width = figure_width, height = figure_height, 
         device = "png", dpi = 300, bg = "white")
  
  message(paste("最终对齐图片已保存:", pdf_file))
  
  # 返回结果用于进一步检查
  return(list(
    final_plot = final_plot,
    base_plot = base_plot,
    positions = positions_df,
    annotations = if(exists("annotations")) annotations else NULL,
    contrast_levels = contrast_levels,
    cell_types = cell_types
  ))
}

# 使用示例 / Usage Example:
# 
# # 1. 准备数据
# props_data <- prepare_data_for_plotting(
#   seurat_obj = your_seurat_object,
#   cell_type_col = "Annotation",
#   sample_col = "sample",
#   contrast_col = "tissue"
# )
# 
# # 2. 创建图表
# publication_plot <- create_publication_boxplot(
#   props_df = props_data,
#   stats_results = your_masc_results$all_results,  # 来自scMASC分析
#   output_dir = "publication_figures",
#   figure_width = 14,
#   figure_height = 10,
#   cell_types_order = c("Epithelial cell", "Fibroblast", "Endothelial cell", 
#                        "Smooth muscle cell", "B cell", "T cell", 
#                        "Myeloid cell", "Cycling cell"),
#   show_individual_points = TRUE,
#   add_statistics = TRUE
# )


# 使用示例和建议工作流程
cat("
=== 使用建议工作流程 ===

1. 首先运行诊断:
diagnosis <- diagnose_masc_convergence(seurat_obj)

2. 根据诊断结果调整参数:
results <- scPairwiseMASCAnalysis(
  seurat_obj = seurat_obj,
  min_cells = diagnosis$recommendations$min_cells,
  min_samples = diagnosis$recommendations$min_samples,
  exclude_filter = 'cell_type_col %in% c(\"problematic_types\")'
)

3. 如果仍有收敛问题，可以直接使用improved_MASC函数

注意事项:
- 收敛警告不一定意味着结果不可信，但需要仔细检查
- 对于关键发现，建议使用多种方法验证
- 考虑增加样本量或合并相似的细胞类型
")

cat("
sample:

diagnosis <- diagnose_masc_convergence(seurat_obj_main)

results <- scPairwiseMASCAnalysis(
  seurat_obj = seurat_obj_main,
  cell_type_col = 'Annotation',
  sample_col = 'sample',
  contrast_col = 'tissue',
  exclude_filter = 'tissue_sampling_method == 'scraping'',
  fixed_effects_cols = NULL,  # 添加任何固定效应协变量
  output_dir = 'pairwise_MASC_results_0928',
  # 如果您需要少于2个样本也能进行比较，可以降低此阈值
  min_samples = 2
)
props_data <- smart_prepare_plotting_data(
  masc_results = results,  # 先尝试从MASC结果提取
  seurat_obj = seurat_obj_main,  # 备用重新计算
  cell_type_col = 'Annotation',
  sample_col = 'sample',
  contrast_col = 'tissue',
  exclude_filter = 'tissue_sampling_method != 'scraping'',  # 必须与MASC相同！
  min_samples = 2,
  min_cells = 10,
  min_prop = 0.01
)


# 第三步：验证比例数据
print('比例数据预览:')
head(props_data)

print('与MASC结果中的组织是否一致:')
print('绘图数据中的组织:')
table(props_data$Contrast)

print('MASC结果中的组织:')
if('all_results' %in% names(results)) {
  print(table(results$all_results$Comparison))
}

# 第四步：创建完美对齐的图表
final_result <- create_separate_aligned_plot(
  props_df = props_data,
  stats_results = results$all_results,
  output_dir = 'consistent_plots',
  cell_types_order = c('Epithelial', 'Endothelial', 'Myeloid', 'SMC', 'T', 'Fibroblast', 'B'),
  contrast_colors = c(
    'lung parenchyma' = '#2E8B57',
    'sinus' = '#B22222', 
    'nose' = '#CD853F',
    'respiratory airway' = '#4682B4'
  )
)
")