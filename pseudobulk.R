#' T细胞组织间差异与通路活性分析
#'
#' @description 对T细胞scRNA-seq数据进行组织间差异分析和通路活性评估
#' @param seurat_obj Seurat对象，包含T细胞数据
#' @param cell_anno_col 细胞类型注释列名，默认为"Annotation_2"
#' @param tissue_col 组织类型注释列名，默认为"tissue"
#' @param sample_col 样本ID注释列名，默认为"sample"
#' @param min_cell_per_sample 每个样本中最少细胞数量，默认为3
#' @param min_sample_per_tissue 每个组织中最少样本数量，默认为3
#' @param run_gsva 是否运行GSVA分析，默认为TRUE
#' @param run_go 是否运行GO富集分析，默认为TRUE
#' @param output_dir 输出目录，默认为当前目录
#' @param species 物种，默认为"Homo sapiens"
#' @return 不返回值，结果保存到指定目录
#' @import Seurat DESeq2 GSVA msigdbr pheatmap ggplot2 reshape2 dplyr
#' @export
#'
#' @examples
#' run_tissue_comparison_analysis(seurat_obj = T_object,
#'                               output_dir = "./results")
run_tissue_comparison_analysis <- function(
  seurat_obj,
  cell_anno_col = "Annotation_2",
  tissue_col = "tissue",
  sample_col = "sample",
  min_cell_per_sample = 3,
  min_sample_per_tissue = 3,
  run_gsva = TRUE,
  run_go = TRUE,
  output_dir = "./analysis_results",
  species = "Homo sapiens"
) {
  # 检查依赖包
  required_packages <- c(
    "Seurat",
    "DESeq2",
    "GSVA",
    "msigdbr",
    "pheatmap",
    "ggplot2",
    "reshape2",
    "dplyr"
  )

  # 加载所需包
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("请安装", pkg, "包"))
    }
    library(pkg, character.only = TRUE)
  }

  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 检查GO富集需要的包
  if (run_go) {
    if (
      !requireNamespace("clusterProfiler", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)
    ) {
      message("未安装clusterProfiler或org.Hs.eg.db，将跳过GO富集分析")
      run_go <- FALSE
    } else {
      library(clusterProfiler)
      library(org.Hs.eg.db)
    }
  }

  # 检查输入对象
  if (!inherits(seurat_obj, "Seurat")) {
    stop("输入对象必须是Seurat对象")
  }

  # 检查必要的元数据列
  required_columns <- c(cell_anno_col, tissue_col, sample_col)
  missing_columns <- required_columns[
    !required_columns %in% colnames(seurat_obj@meta.data)
  ]

  if (length(missing_columns) > 0) {
    stop(paste(
      "Seurat对象缺少必要的元数据列:",
      paste(missing_columns, collapse = ", ")
    ))
  }

  # 显示数据集概况
  message("==== 数据集概况 ====")
  message(paste0("共", ncol(seurat_obj), "个细胞"))
  message("组织分布:")
  print(table(seurat_obj@meta.data[[tissue_col]]))
  message("细胞类型分布:")
  print(table(seurat_obj@meta.data[[cell_anno_col]]))

  # 1. Pseudobulk分析
  ########################

  # 创建pseudobulk结果目录
  pseudobulk_dir <- file.path(output_dir, "pseudobulk_analysis")
  dir.create(pseudobulk_dir, showWarnings = FALSE)

  # 设置工作目录
  original_dir <- getwd()
  setwd(pseudobulk_dir)

  # 分析不同组织中的细胞类型分布
  message("分析不同组织中的细胞类型分布...")
  tissue_cell_distribution <- table(
    seurat_obj@meta.data[[tissue_col]],
    seurat_obj@meta.data[[cell_anno_col]]
  )
  write.csv(tissue_cell_distribution, "tissue_celltype_distribution.csv")

  # 按组织、样本和细胞类型聚合表达数据
  message("执行基于组织和细胞类型的pseudobulk聚合...")
  group_by_cols <- c(tissue_col, sample_col, cell_anno_col)
  av_tissue_celltype <- AggregateExpression(
    seurat_obj,
    group.by = group_by_cols,
    # assays = "RNA",
    # layer = "counts",
    return.seurat = FALSE
  )

  # 处理聚合结果
  av_tissue_celltype_df <- as.data.frame(av_tissue_celltype[["RNA"]])
  write.csv(
    av_tissue_celltype_df,
    file = "pseudobulk_tissue_sample_celltype.csv"
  )

  # 从列名中提取元数据
  extract_metadata <- function(column_names) {
    # 按下划线分割列名
    parts <- strsplit(column_names, "_")

    # 创建结果数据框
    result <- data.frame(
      column = column_names,
      tissue = sapply(parts, function(x) x[1]),
      sample = sapply(parts, function(x) x[2]),
      celltype = sapply(parts, function(x) paste(x[-(1:2)], collapse = "_")),
      stringsAsFactors = FALSE
    )

    return(result)
  }

  # 提取元数据
  metadata <- extract_metadata(colnames(av_tissue_celltype_df))
  rownames(metadata) <- metadata$column

  # 为每个细胞类型执行组织间差异分析的函数
  run_pairwise_tissue_comparison <- function(
    counts_matrix,
    metadata,
    cell_type
  ) {
    # 筛选特定细胞类型的样本
    cell_indices <- which(metadata$celltype == cell_type)

    if (length(cell_indices) < (min_sample_per_tissue * 2)) {
      message(paste(
        "细胞类型",
        cell_type,
        "的样本数量不足(",
        length(cell_indices),
        ")，跳过分析"
      ))
      return(NULL)
    }

    # 提取该细胞类型的计数数据和元数据
    cell_counts <- counts_matrix[, cell_indices, drop = FALSE]
    cell_metadata <- metadata[cell_indices, , drop = FALSE]

    # 查看组织分布
    tissue_counts <- table(cell_metadata$tissue)
    message(paste("细胞类型", cell_type, "在各组织中的样本数:"))
    print(tissue_counts)

    # 只保留至少有min_sample_per_tissue个样本的组织
    valid_tissues <- names(tissue_counts[
      tissue_counts >= min_sample_per_tissue
    ])

    if (length(valid_tissues) < 2) {
      message(paste(
        "细胞类型",
        cell_type,
        "中至少有",
        min_sample_per_tissue,
        "个样本的组织不足2种，跳过比较"
      ))
      return(NULL)
    }

    # 过滤掉样本数不足的组织
    valid_indices <- cell_metadata$tissue %in% valid_tissues
    cell_counts <- cell_counts[, valid_indices, drop = FALSE]
    cell_metadata <- cell_metadata[valid_indices, , drop = FALSE]

    # 确保counts_matrix是整数
    cell_counts <- round(cell_counts)

    # 生成所有可能的两两组织比较
    tissue_pairs <- combn(valid_tissues, 2, simplify = FALSE)
    message(paste(
      "为细胞类型",
      cell_type,
      "执行",
      length(tissue_pairs),
      "个组织间两两比较"
    ))

    results_list <- list()

    # 对每对组织进行比较
    for (pair in tissue_pairs) {
      tissue1 <- pair[1]
      tissue2 <- pair[2]

      message(paste("比较:", tissue1, "vs", tissue2))

      # 只选择这两个组织的样本
      pair_indices <- cell_metadata$tissue %in% pair
      pair_counts <- cell_counts[, pair_indices, drop = FALSE]
      pair_metadata <- cell_metadata[pair_indices, , drop = FALSE]

      # 再次检查每个组织的样本数
      tissue_sample_counts <- table(pair_metadata$tissue)
      if (any(tissue_sample_counts < min_sample_per_tissue)) {
        message(paste(
          "跳过比较: 某个组织样本数少于",
          min_sample_per_tissue,
          ":",
          paste(
            names(tissue_sample_counts),
            tissue_sample_counts,
            sep = "=",
            collapse = ", "
          )
        ))
        next
      }

      # 创建DESeq2数据集
      dds <- DESeqDataSetFromMatrix(
        countData = pair_counts,
        colData = pair_metadata,
        design = ~tissue
      )

      # 设置参考水平
      dds$tissue <- relevel(factor(dds$tissue), ref = tissue1)

      # 过滤低表达基因
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep, ]

      # 运行DESeq2
      dds <- DESeq(dds)

      # 获取差异结果
      comparison_name <- paste0(tissue2, "_vs_", tissue1)
      res <- results(dds, contrast = c("tissue", tissue2, tissue1))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df <- res_df[order(res_df$padj), ]
      res_df <- na.omit(res_df)

      # 添加上下调信息
      res_df$regulation <- ifelse(
        res_df$padj > 0.05,
        "stable",
        ifelse(
          abs(res_df$log2FoldChange) < 1,
          "stable",
          ifelse(res_df$log2FoldChange >= 1, "up", "down")
        )
      )

      # 保存结果
      results_list[[comparison_name]] <- res_df

      # 创建输出目录
      result_dir <- paste0(cell_type, "_", tissue1, "_vs_", tissue2)
      result_dir <- gsub("[^a-zA-Z0-9_]", "-", result_dir) # 替换非法字符
      dir.create(result_dir, showWarnings = FALSE)

      # 保存到CSV
      output_file <- file.path(result_dir, "DEGs.csv")
      write.csv(res_df, file = output_file, row.names = FALSE)

      # 生成摘要统计
      summary_stats <- data.frame(
        Comparison = comparison_name,
        Total_DEGs = sum(res_df$regulation != "stable"),
        Up_regulated = sum(res_df$regulation == "up"),
        Down_regulated = sum(res_df$regulation == "down"),
        Total_genes = nrow(res_df)
      )

      write.csv(
        summary_stats,
        file = file.path(result_dir, "summary_stats.csv"),
        row.names = FALSE
      )

      # 创建火山图
      volcano_file <- file.path(result_dir, "volcano_plot.pdf")
      pdf(volcano_file, width = 10, height = 8)
      p <- ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
        geom_point(size = 1.5, alpha = 0.7, aes(color = regulation)) +
        scale_color_manual(
          values = c("down" = "#00468B", "stable" = "gray", "up" = "#E64B35")
        ) +
        labs(
          x = "Log2(fold change)",
          y = "-log10(adjusted p-value)",
          title = paste0(cell_type, ": ", tissue2, " vs ", tissue1)
        ) +
        geom_hline(
          yintercept = -log10(0.05),
          linetype = 2,
          color = 'black',
          linewidth = 0.5
        ) +
        geom_vline(
          xintercept = c(-1, 1),
          linetype = 2,
          color = 'black',
          linewidth = 0.5
        ) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
      print(p)
      dev.off()

      # 统计差异基因数量
      message(paste0(
        cell_type,
        " ",
        comparison_name,
        ": ",
        sum(res_df$regulation == "up"),
        " 上调基因, ",
        sum(res_df$regulation == "down"),
        " 下调基因"
      ))

      # 创建主要差异基因热图
      top_degs <- res_df %>%
        filter(regulation != "stable") %>%
        arrange(padj) %>%
        head(50) # 取前50个差异基因

      if (nrow(top_degs) > 0) {
        # 提取这些基因在所有样本中的表达值
        gene_expr <- pair_counts[top_degs$gene, , drop = FALSE]

        # 标准化
        gene_expr_norm <- t(scale(t(log2(gene_expr + 1))))

        # 创建热图注释
        anno_col <- data.frame(Tissue = pair_metadata$tissue)
        rownames(anno_col) <- colnames(gene_expr_norm)

        # 绘制热图
        heatmap_file <- file.path(result_dir, "DEG_heatmap.pdf")
        pdf(heatmap_file, width = 12, height = 10)
        pheatmap(
          gene_expr_norm,
          annotation_col = anno_col,
          main = paste0(cell_type, " DEGs: ", tissue2, " vs ", tissue1),
          fontsize_row = 8,
          show_colnames = FALSE
        )
        dev.off()

        # 执行GO富集分析
        if (run_go) {
          message("执行GO富集分析...")

          # 提取上调和下调基因
          up_genes <- res_df$gene[res_df$regulation == "up"]
          down_genes <- res_df$gene[res_df$regulation == "down"]

          # GO富集分析 - 上调基因
          if (length(up_genes) >= 10) {
            ego_up <- enrichGO(
              gene = up_genes,
              OrgDb = org.Hs.eg.db,
              keyType = "SYMBOL",
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.2
            )

            if (!is.null(ego_up) && nrow(ego_up) > 0) {
              # 保存结果
              write.csv(
                as.data.frame(ego_up),
                file = file.path(result_dir, "GO_upregulated.csv"),
                row.names = FALSE
              )

              # 绘制GO富集点图
              pdf(
                file.path(result_dir, "GO_upregulated_dotplot.pdf"),
                width = 10,
                height = 8
              )
              print(dotplot(
                ego_up,
                showCategory = 20,
                title = "GO Enrichment: Upregulated Genes"
              ))
              dev.off()
            }
          }

          # GO富集分析 - 下调基因
          if (length(down_genes) >= 10) {
            ego_down <- enrichGO(
              gene = down_genes,
              OrgDb = org.Hs.eg.db,
              keyType = "SYMBOL",
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.2
            )

            if (!is.null(ego_down) && nrow(ego_down) > 0) {
              # 保存结果
              write.csv(
                as.data.frame(ego_down),
                file = file.path(result_dir, "GO_downregulated.csv"),
                row.names = FALSE
              )

              # 绘制GO富集点图
              pdf(
                file.path(result_dir, "GO_downregulated_dotplot.pdf"),
                width = 10,
                height = 8
              )
              print(dotplot(
                ego_down,
                showCategory = 20,
                title = "GO Enrichment: Downregulated Genes"
              ))
              dev.off()
            }
          }
        }
      }
    }

    return(results_list)
  }

  # 对每个细胞类型执行组织间两两比较
  all_cell_types <- unique(metadata$celltype)
  results_by_celltype <- list()

  for (cell_type in all_cell_types) {
    message(paste("分析细胞类型:", cell_type))

    # 跳过NA或空字符串
    if (is.na(cell_type) || cell_type == "") {
      message("跳过无效细胞类型")
      next
    }

    # 创建细胞类型目录
    cell_dir <- gsub("[^a-zA-Z0-9]", "_", cell_type)
    dir.create(cell_dir, showWarnings = FALSE)

    # 设置工作目录
    current_dir <- getwd()
    setwd(cell_dir)

    # 执行组织间两两比较
    tryCatch(
      {
        results <- run_pairwise_tissue_comparison(
          av_tissue_celltype_df,
          metadata,
          cell_type
        )

        if (!is.null(results)) {
          results_by_celltype[[cell_type]] <- results
        }
      },
      error = function(e) {
        message(paste("分析细胞类型", cell_type, "时出错:", e$message))
      }
    )

    # 返回原目录
    setwd(current_dir)
  }

  # 创建组织间比较的总结报告
  message("创建组织间比较的总结报告...")

  # 收集所有比较的结果
  all_comparisons <- data.frame()

  # 遍历目录结构收集结果
  for (cell_dir in list.dirs(recursive = FALSE)) {
    cell_type <- basename(cell_dir)

    # 跳过非细胞类型目录
    if (!dir.exists(cell_dir) || cell_dir == "." || cell_dir == "..") {
      next
    }

    # 查找所有比较目录
    comparison_dirs <- list.dirs(path = cell_dir, recursive = FALSE)

    for (comp_dir in comparison_dirs) {
      # 检查是否有summary_stats.csv
      stats_file <- file.path(comp_dir, "summary_stats.csv")

      if (file.exists(stats_file)) {
        stats <- read.csv(stats_file)
        stats$CellType <- cell_type

        all_comparisons <- rbind(all_comparisons, stats)
      }
    }
  }

  # 保存总结表
  if (nrow(all_comparisons) > 0) {
    write.csv(
      all_comparisons,
      "all_pairwise_comparisons_summary.csv",
      row.names = FALSE
    )

    # 创建总结图
    # 条形图：每个细胞类型的差异基因数
    pdf("diff_gene_counts_by_celltype.pdf", width = 15, height = 10)
    p <- ggplot(
      all_comparisons,
      aes(x = CellType, y = Total_DEGs, fill = Comparison)
    ) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "差异基因数量（按细胞类型）",
        x = "细胞类型",
        y = "差异基因数量"
      )
    print(p)
    dev.off()

    # 堆叠条形图：每个细胞类型的上调/下调基因比例
    all_comparisons_long <- reshape2::melt(
      all_comparisons,
      id.vars = c("CellType", "Comparison", "Total_genes"),
      measure.vars = c("Up_regulated", "Down_regulated"),
      variable.name = "Regulation",
      value.name = "Count"
    )

    pdf("up_down_genes_by_comparison.pdf", width = 15, height = 10)
    p <- ggplot(
      all_comparisons_long,
      aes(x = Comparison, y = Count, fill = Regulation)
    ) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~CellType, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(
        values = c("Up_regulated" = "#E64B35", "Down_regulated" = "#00468B")
      ) +
      labs(
        title = "上调和下调基因数量（按细胞类型和组织比较）",
        x = "组织比较",
        y = "基因数量"
      )
    print(p)
    dev.off()

    # 热图：所有组织比较的差异基因数
    heatmap_data <- reshape2::dcast(
      all_comparisons,
      CellType ~ Comparison,
      value.var = "Total_DEGs"
    )
    rownames(heatmap_data) <- heatmap_data$CellType
    heatmap_data <- heatmap_data[, -1]

    # 替换NA值为0
    heatmap_data[is.na(heatmap_data)] <- 0

    # 绘制热图
    pdf("comparison_heatmap.pdf", width = 12, height = 10)
    pheatmap(
      heatmap_data,
      display_numbers = TRUE,
      main = "组织间比较的差异基因数量热图",
      fontsize_number = 8,
      fontsize = 10
    )
    dev.off()
  }

  message("组织间细胞类型两两比较分析完成！")

  # 2. GSVA分析
  #################
  if (run_gsva) {
    # 创建GSVA目录
    gsva_dir <- file.path(output_dir, "gsva_analysis")
    dir.create(gsva_dir, showWarnings = FALSE)
    setwd(gsva_dir)

    message("开始进行GSVA分析...")

    # 函数定义 - 提高代码重用性
    run_gsva_analysis <- function(
      expr_data,
      gene_sets,
      method = "gsva",
      mx.diff = TRUE,
      output_prefix
    ) {
      # 运行GSVA
      message(paste0("使用", output_prefix, "基因集运行GSVA..."))
      tryCatch(
        {
          gsva_result <- gsva(
            expr_data,
            gene_sets,
            method = method,
            mx.diff = mx.diff
          )
          write.csv(
            gsva_result,
            file = paste0("gsva_", output_prefix, "_results.csv")
          )
          return(gsva_result)
        },
        error = function(e) {
          message(paste0("GSVA分析出错: ", e$message))
          return(NULL)
        }
      )
    }

    # 绘制GSVA热图函数
    plot_gsva_heatmap <- function(
      gsva_result,
      top_n = 20,
      output_prefix,
      title
    ) {
      if (is.null(gsva_result) || nrow(gsva_result) == 0) {
        message(paste0("无法绘制", output_prefix, "热图: 结果为空"))
        return(NULL)
      }

      # 确保有足够行数
      if (nrow(gsva_result) < 2) {
        message(paste0(
          "无法绘制",
          output_prefix,
          "热图: 行数太少(",
          nrow(gsva_result),
          ")"
        ))
        return(NULL)
      }

      # 计算行方差
      row_vars <- apply(gsva_result, 1, var)

      # 检查是否有方差大于零的行
      if (sum(row_vars > 0) == 0) {
        message(paste0("无法绘制", output_prefix, "热图: 所有行方差为零"))
        return(NULL)
      }

      # 选择方差最大的通路
      select_n <- min(top_n, sum(row_vars > 0))
      var_pathways <- names(tail(sort(row_vars), select_n))

      if (length(var_pathways) == 0) {
        message(paste0(
          "无法绘制",
          output_prefix,
          "热图: 无法选择方差最大的通路"
        ))
        return(NULL)
      }

      message(paste0(
        "为",
        output_prefix,
        "热图选择了",
        length(var_pathways),
        "个通路"
      ))

      # 绘制热图
      pdf_path <- paste0(output_prefix, "_gsva_heatmap.pdf")
      pdf(
        pdf_path,
        width = 10,
        height = max(4, ceiling(length(var_pathways) / 3))
      )

      tryCatch(
        {
          # 设置颜色调色板
          col_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

          # 绘制热图并确保打印
          hmap <- pheatmap(
            gsva_result[var_pathways, ],
            scale = "row",
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            fontsize_row = 8,
            fontsize_col = 10,
            angle_col = "45", # 修正为字符串类型
            main = title,
            color = col_palette
          )
          # 强制打印热图
          print(hmap)
          message(paste0("成功绘制热图并保存到", pdf_path))
          invisible(hmap)
        },
        error = function(e) {
          message(paste0("绘制热图出错: ", e$message))
          # 确保在出错时生成有效的PDF
          plot(1, type = "n", axes = FALSE, ann = FALSE)
          text(1, 1, paste0("绘制热图出错: ", e$message), cex = 1.2)
          NULL
        },
        finally = {
          # 确保PDF设备已关闭
          dev.off()
        }
      )
    }

    # 绘制GSVA气泡图函数
    plot_gsva_bubble <- function(
      gsva_result,
      var_pathways,
      output_prefix,
      title
    ) {
      if (is.null(gsva_result) || nrow(gsva_result) == 0) {
        message(paste0("无法绘制", output_prefix, "气泡图: 结果为空"))
        return(NULL)
      }

      # 确保var_pathways不为空且存在于gsva_result中
      if (length(var_pathways) == 0) {
        message(paste0("无法绘制", output_prefix, "气泡图: 通路列表为空"))
        return(NULL)
      }

      # 只保留存在于gsva_result中的通路
      valid_pathways <- intersect(var_pathways, rownames(gsva_result))
      if (length(valid_pathways) == 0) {
        message(paste0("无法绘制", output_prefix, "气泡图: 无有效通路"))
        return(NULL)
      }

      message(paste0(
        "使用",
        length(valid_pathways),
        "个通路绘制",
        output_prefix,
        "气泡图"
      ))

      # 绘制气泡图
      pdf_path <- paste0(output_prefix, "_gsva_bubble.pdf")
      pdf(pdf_path, width = 12, height = 10)

      tryCatch(
        {
          # 准备长格式数据
          gsva_data_long <- reshape2::melt(
            gsva_result[valid_pathways, ],
            varnames = c("Pathway", "CellType"),
            value.name = "GSVA_Score"
          )

          # 绘制气泡图
          p <- ggplot(
            gsva_data_long,
            aes(
              x = CellType,
              y = Pathway,
              size = abs(GSVA_Score),
              color = GSVA_Score
            )
          ) +
            geom_point(alpha = 0.8) +
            scale_size_continuous(range = c(1, 8), name = "绝对分数") +
            scale_color_gradient2(
              low = "blue",
              mid = "white",
              high = "red",
              midpoint = 0,
              name = "GSVA分数"
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(face = "bold"),
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
            ) +
            labs(title = title, x = "细胞类型", y = "通路")
          # 强制打印图形
          print(p)
          message(paste0("成功绘制气泡图并保存到", pdf_path))
          invisible(p)
        },
        error = function(e) {
          message(paste0("绘制气泡图出错: ", e$message))
          # 确保在出错时生成有效的PDF
          plot(1, type = "n", axes = FALSE, ann = FALSE)
          text(1, 1, paste0("绘制气泡图出错: ", e$message), cex = 1.2)
          NULL
        },
        finally = {
          # 确保PDF设备已关闭
          dev.off()
        }
      )
    }

    # 设置细胞类型为标识
    Idents(seurat_obj) <- seurat_obj[[cell_anno_col]]
    message("按细胞类型计算平均表达...")

    # 计算每个细胞类型的平均表达
    expr <- AverageExpression(seurat_obj, assays = "RNA", layer = "data")[[1]]
    expr <- expr[rowSums(expr) > 0, ] # 过滤非表达基因
    expr <- as.matrix(expr)
    message(paste0(
      "共获取",
      nrow(expr),
      "个基因在",
      ncol(expr),
      "个细胞类型中的平均表达"
    ))

    # 从MSigDB获取Hallmark和C2 KEGG基因集
    message("获取基因集数据...")
    hallmark_genesets <- msigdbr(species = species, category = "H")
    hallmark_genesets <- subset(
      hallmark_genesets,
      select = c("gs_name", "gene_symbol")
    ) %>%
      as.data.frame()
    hallmark_genesets <- split(
      hallmark_genesets$gene_symbol,
      hallmark_genesets$gs_name
    )
    message(paste0("成功获取", length(hallmark_genesets), "个Hallmark基因集"))

    # KEGG基因集
    c2_kegg_genesets <- msigdbr(
      species = species,
      category = "C2",
      subcategory = "KEGG"
    )
    c2_kegg_genesets <- subset(
      c2_kegg_genesets,
      select = c("gs_name", "gene_symbol")
    ) %>%
      as.data.frame()
    c2_kegg_genesets <- split(
      c2_kegg_genesets$gene_symbol,
      c2_kegg_genesets$gs_name
    )
    message(paste0("成功获取", length(c2_kegg_genesets), "个KEGG基因集"))

    # 运行Hallmark GSVA
    gsva_result_hallmark <- run_gsva_analysis(
      expr_data = expr,
      gene_sets = hallmark_genesets,
      output_prefix = "hallmark"
    )

    # 运行KEGG GSVA
    gsva_result_kegg <- run_gsva_analysis(
      expr_data = expr,
      gene_sets = c2_kegg_genesets,
      output_prefix = "kegg"
    )

    # 可视化GSVA结果
    # Hallmark基因集热图
    if (!is.null(gsva_result_hallmark)) {
      # 绘制热图
      var_pathways_hallmark <- names(tail(
        sort(apply(gsva_result_hallmark, 1, var)),
        20
      ))
      plot_gsva_heatmap(
        gsva_result = gsva_result_hallmark,
        top_n = 20,
        output_prefix = "hallmark",
        title = "Hallmark基因集GSVA得分"
      )

      # 绘制气泡图
      plot_gsva_bubble(
        gsva_result = gsva_result_hallmark,
        var_pathways = var_pathways_hallmark,
        output_prefix = "hallmark",
        title = "Hallmark通路GSVA得分"
      )
    }

    # KEGG基因集热图
    if (!is.null(gsva_result_kegg)) {
      # 绘制热图
      var_pathways_kegg <- names(tail(
        sort(apply(gsva_result_kegg, 1, var)),
        30
      ))
      plot_gsva_heatmap(
        gsva_result = gsva_result_kegg,
        top_n = 30,
        output_prefix = "kegg",
        title = "KEGG通路GSVA得分"
      )

      # 绘制气泡图
      plot_gsva_bubble(
        gsva_result = gsva_result_kegg,
        var_pathways = var_pathways_kegg,
        output_prefix = "kegg",
        title = "KEGG通路GSVA得分"
      )
    }

    # 按组织和细胞类型进行GSVA分析
    message("按组织和细胞类型进行GSVA分析...")

    # 按组织分组，计算每个细胞类型在每个组织中的平均表达
    tissue_types <- unique(seurat_obj@meta.data[[tissue_col]])
    tissue_gsva_results <- list()

    for (tissue in tissue_types) {
      tissue_cells <- WhichCells(
        seurat_obj,
        expression = get(tissue_col) == tissue
      )
      if (length(tissue_cells) < 50) {
        # 确保有足够细胞进行分析
        message(paste0(
          "组织",
          tissue,
          "的细胞数量不足(",
          length(tissue_cells),
          ")，跳过分析"
        ))
        next
      }

      tissue_obj <- subset(seurat_obj, cells = tissue_cells)
      Idents(tissue_obj) <- tissue_obj[[cell_anno_col]]

      # 计算该组织中每个细胞类型的平均表达
      tissue_expr <- AverageExpression(
        tissue_obj,
        assays = "RNA",
        layer = "data"
      )[[1]]
      tissue_expr <- tissue_expr[rowSums(tissue_expr) > 0, ] # 过滤非表达基因
      tissue_expr <- as.matrix(tissue_expr)

      message(paste0(
        "组织",
        tissue,
        "中共有",
        ncol(tissue_expr),
        "个细胞类型的平均表达"
      ))

      # 运行GSVA
      # Hallmark基因集
      tissue_gsva_hallmark <- run_gsva_analysis(
        expr_data = tissue_expr,
        gene_sets = hallmark_genesets,
        output_prefix = paste0(tissue, "_hallmark")
      )

      # KEGG基因集
      tissue_gsva_kegg <- run_gsva_analysis(
        expr_data = tissue_expr,
        gene_sets = c2_kegg_genesets,
        output_prefix = paste0(tissue, "_kegg")
      )

      # 可视化结果
      if (!is.null(tissue_gsva_hallmark)) {
        # 绘制热图
        plot_gsva_heatmap(
          gsva_result = tissue_gsva_hallmark,
          top_n = 20,
          output_prefix = paste0(tissue, "_hallmark"),
          title = paste0(tissue, " Hallmark基因集GSVA得分")
        )
      }

      if (!is.null(tissue_gsva_kegg)) {
        # 绘制热图
        plot_gsva_heatmap(
          gsva_result = tissue_gsva_kegg,
          top_n = 30,
          output_prefix = paste0(tissue, "_kegg"),
          title = paste0(tissue, " KEGG通路GSVA得分")
        )
      }

      # 保存结果
      tissue_gsva_results[[tissue]] <- list(
        hallmark = tissue_gsva_hallmark,
        kegg = tissue_gsva_kegg
      )
    }

    # 组织间GSVA差异比较
    if (length(tissue_types) >= 2 && length(tissue_gsva_results) >= 2) {
      message("执行组织间通路活性比较...")

      # 获取所有可能的组织对
      tissue_pairs <- combn(names(tissue_gsva_results), 2, simplify = FALSE)

      for (pair in tissue_pairs) {
        tissue_1 <- pair[1]
        tissue_2 <- pair[2]

        message(paste0("比较: ", tissue_1, " vs ", tissue_2))

        # 对Hallmark和KEGG分别进行比较
        for (pathway_type in c("hallmark", "kegg")) {
          gsva_tissue_1 <- tissue_gsva_results[[tissue_1]][[pathway_type]]
          gsva_tissue_2 <- tissue_gsva_results[[tissue_2]][[pathway_type]]

          if (is.null(gsva_tissue_1) || is.null(gsva_tissue_2)) {
            message(paste0("跳过", pathway_type, "比较: 缺少数据"))
            next
          }

          # 查找两个组织之间的共同细胞类型
          common_cell_types <- intersect(
            colnames(gsva_tissue_1),
            colnames(gsva_tissue_2)
          )

          if (length(common_cell_types) == 0) {
            message(paste0(
              "组织'",
              tissue_1,
              "'和'",
              tissue_2,
              "'之间没有共同的细胞类型"
            ))
            next
          }

          message(paste0("找到", length(common_cell_types), "个共同细胞类型"))

          # 查找共同通路
          common_pathways <- intersect(
            rownames(gsva_tissue_1),
            rownames(gsva_tissue_2)
          )

          if (length(common_pathways) == 0) {
            message(paste0(
              "组织'",
              tissue_1,
              "'和'",
              tissue_2,
              "'之间没有共同通路"
            ))
            next
          }

          # 对每个细胞类型计算通路活性差异
          for (cell_type in common_cell_types) {
            tryCatch(
              {
                message(paste0("处理细胞类型: ", cell_type))

                # 提取两个组织中该细胞类型的GSVA得分
                cell_gsva_1 <- gsva_tissue_1[
                  common_pathways,
                  cell_type,
                  drop = FALSE
                ]
                cell_gsva_2 <- gsva_tissue_2[
                  common_pathways,
                  cell_type,
                  drop = FALSE
                ]

                # 计算差异
                diff_gsva <- cell_gsva_1 - cell_gsva_2
                colnames(diff_gsva) <- paste0(cell_type, "_diff")

                # 保存差异结果
                output_csv <- paste0(
                  "gsva_diff_",
                  pathway_type,
                  "_",
                  gsub(" ", "_", cell_type),
                  "_",
                  gsub(" ", "_", tissue_1),
                  "_vs_",
                  gsub(" ", "_", tissue_2),
                  ".csv"
                )
                write.csv(diff_gsva, file = output_csv)
                message(paste0("差异结果保存至: ", output_csv))

                # 获取差异最大的通路
                if (nrow(diff_gsva) == 0) {
                  message(paste0(
                    "细胞类型'",
                    cell_type,
                    "'没有差异数据，跳过可视化"
                  ))
                  next
                }

                # 计算并排序绝对差异
                abs_diffs <- abs(diff_gsva[, 1])
                if (
                  length(abs_diffs) == 0 ||
                    all(is.na(abs_diffs)) ||
                    all(abs_diffs == 0)
                ) {
                  message(paste0(
                    "细胞类型'",
                    cell_type,
                    "'没有有效差异，跳过可视化"
                  ))
                  next
                }

                # 按绝对差异值排序
                sorted_idx <- order(abs_diffs, decreasing = TRUE)
                # 选择前15个或所有通路(如果少于15个)
                n_paths <- min(15, length(sorted_idx))
                top_idx <- sorted_idx[1:n_paths]

                # 确保索引有效
                if (length(top_idx) == 0) {
                  message(paste0(
                    "细胞类型'",
                    cell_type,
                    "'无法选择差异最大的通路"
                  ))
                  next
                }

                # 获取通路名称
                top_diff_pathways <- rownames(diff_gsva)[top_idx]
                message(paste0(
                  "选择了",
                  length(top_diff_pathways),
                  "个差异最大的通路"
                ))

                # 绘制条形图
                pdf_path <- paste0(
                  "gsva_diff_barplot_",
                  pathway_type,
                  "_",
                  gsub(" ", "_", cell_type),
                  "_",
                  gsub(" ", "_", tissue_1),
                  "_vs_",
                  gsub(" ", "_", tissue_2),
                  ".pdf"
                )
                pdf(pdf_path, width = 10, height = 8)

                tryCatch(
                  {
                    # 按差异值排序通路
                    ordered_paths <- top_diff_pathways[order(diff_gsva[
                      top_diff_pathways,
                      1
                    ])]

                    # 创建数据框
                    barplot_data <- data.frame(
                      Pathway = factor(ordered_paths, levels = ordered_paths),
                      Difference = diff_gsva[ordered_paths, 1]
                    )

                    # 绘制条形图
                    p <- ggplot(
                      barplot_data,
                      aes(x = Difference, y = Pathway, fill = Difference > 0)
                    ) +
                      geom_bar(stat = "identity") +
                      scale_fill_manual(
                        values = c("TRUE" = "#E64B35", "FALSE" = "#00468B"),
                        labels = c("TRUE" = tissue_1, "FALSE" = tissue_2),
                        name = "更高活性在"
                      ) +
                      labs(
                        title = paste0(cell_type, " 通路活性差异"),
                        subtitle = paste0(tissue_1, " vs ", tissue_2),
                        x = "GSVA得分差异"
                      ) +
                      theme_bw() +
                      theme(
                        axis.text.y = element_text(size = 9),
                        plot.title = element_text(
                          hjust = 0.5,
                          size = 14,
                          face = "bold"
                        ),
                        plot.subtitle = element_text(hjust = 0.5, size = 12)
                      )

                    # 确保图被打印
                    print(p)
                    message(paste0("成功绘制条形图并保存至", pdf_path))
                  },
                  error = function(e) {
                    message(paste0("绘制条形图出错: ", e$message))
                    # 确保在出错时生成有效的PDF
                    plot(1, type = "n", axes = FALSE, ann = FALSE)
                    text(1, 1, paste0("绘制条形图出错: ", e$message), cex = 1.2)
                  },
                  finally = {
                    # 确保PDF设备已关闭
                    dev.off()
                  }
                )
              },
              error = function(e) {
                message(paste0(
                  "处理细胞类型'",
                  cell_type,
                  "'时出错: ",
                  e$message
                ))
              }
            )
          }

          # 创建组织差异热图
          message("创建组织差异热图...")
          if (length(common_cell_types) >= 2 && length(common_pathways) > 0) {
            # 创建差异矩阵
            diff_matrix <- matrix(
              NA,
              nrow = length(common_pathways),
              ncol = length(common_cell_types)
            )
            rownames(diff_matrix) <- common_pathways
            colnames(diff_matrix) <- common_cell_types

            # 填充差异矩阵
            for (i in 1:length(common_cell_types)) {
              cell_type <- common_cell_types[i]
              if (
                cell_type %in%
                  colnames(gsva_tissue_1) &&
                  cell_type %in% colnames(gsva_tissue_2)
              ) {
                diff_matrix[, i] <- gsva_tissue_1[common_pathways, cell_type] -
                  gsva_tissue_2[common_pathways, cell_type]
              }
            }

            # 移除全部为NA的行和列
            diff_matrix <- diff_matrix[
              rowSums(!is.na(diff_matrix)) > 0,
              colSums(!is.na(diff_matrix)) > 0,
              drop = FALSE
            ]

            if (nrow(diff_matrix) > 0 && ncol(diff_matrix) > 0) {
              # 计算方差最大的通路
              if (nrow(diff_matrix) < 2) {
                message("差异矩阵行数太少，无法绘制热图")
              } else {
                row_vars <- apply(diff_matrix, 1, var, na.rm = TRUE)
                valid_rows <- which(!is.na(row_vars) & row_vars > 0)

                if (length(valid_rows) == 0) {
                  message("没有方差大于零的通路，无法绘制热图")
                } else {
                  # 选择方差最大的通路
                  n_paths <- min(25, length(valid_rows))
                  top_var_idx <- tail(order(row_vars[valid_rows]), n_paths)
                  var_paths <- rownames(diff_matrix)[valid_rows[top_var_idx]]

                  # 绘制热图
                  pdf_path <- paste0(
                    "tissue_diff_heatmap_",
                    pathway_type,
                    "_",
                    gsub(" ", "_", tissue_1),
                    "_vs_",
                    gsub(" ", "_", tissue_2),
                    ".pdf"
                  )
                  pdf(pdf_path, width = 10, height = 12)

                  tryCatch(
                    {
                      hmap <- pheatmap(
                        diff_matrix[var_paths, ],
                        scale = "row",
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        fontsize_row = 8,
                        fontsize_col = 10,
                        angle_col = "45",
                        main = paste0(
                          "通路活性差异: ",
                          tissue_1,
                          " vs ",
                          tissue_2
                        ),
                        color = colorRampPalette(c(
                          "#00468B",
                          "white",
                          "#E64B35"
                        ))(50)
                      )

                      # 确保热图被打印
                      print(hmap)
                      message(paste0("成功绘制组织差异热图并保存至", pdf_path))
                    },
                    error = function(e) {
                      message(paste0("绘制组织差异热图出错: ", e$message))
                      # 确保在出错时生成有效的PDF
                      plot(1, type = "n", axes = FALSE, ann = FALSE)
                      text(
                        1,
                        1,
                        paste0("绘制组织差异热图出错: ", e$message),
                        cex = 1.2
                      )
                    },
                    finally = {
                      # 确保PDF设备已关闭
                      dev.off()
                    }
                  )
                }
              }
            } else {
              message("处理后的差异矩阵为空，无法绘制热图")
            }
          } else {
            message("共同细胞类型或通路数量不足，无法创建组织差异热图")
          }
        }
      }
    }

    message("GSVA分析完成！结果保存在:", normalizePath(gsva_dir))
  }

  # 返回到原始目录
  setwd(original_dir)
  message("分析完成！所有结果保存在:", normalizePath(output_dir))
}

# 加载Seurat对象
# T_object应包含tissue、sample和Annotation_2的元数据列
# ==============================================================================
# Simple Layer-Based ENSG to Symbol Conversion
# Strategy: Extract layer → Convert rownames → Save to new layer
# ==============================================================================

# Step 1: Check layer structure
check_layers <- function(seurat_obj) {
  cat("\n=== LAYER STRUCTURE ===\n\n")

  assay_name <- DefaultAssay(seurat_obj)
  cat(sprintf("Default Assay: %s\n", assay_name))

  # List all layers
  layers <- Layers(seurat_obj, assay = assay_name)
  cat(sprintf("\nAvailable layers: %s\n", paste(layers, collapse = ", ")))

  # Check each layer
  cat("\nLayer details:\n")
  for (layer in layers) {
    mat <- LayerData(seurat_obj, layer = layer, assay = assay_name)
    cat(sprintf("  %s: %d genes × %d cells\n", layer, nrow(mat), ncol(mat)))
    cat(sprintf(
      "    Sample genes: %s\n",
      paste(head(rownames(mat), 3), collapse = ", ")
    ))
  }

  return(layers)
}


# Step 2: Extract layer matrix and convert gene names
convert_layer_genes <- function(
  seurat_obj,
  source_layer = "counts",
  handle_duplicates = "unique"
) {
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install org.Hs.eg.db: BiocManager::install('org.Hs.eg.db')")
  }
  library(AnnotationDbi)

  assay_name <- DefaultAssay(seurat_obj)

  cat("\n=== CONVERTING LAYER GENES ===\n\n")
  cat(sprintf("Source layer: %s\n", source_layer))

  # Extract matrix
  cat("\nStep 1: Extracting matrix...\n")
  mat <- LayerData(seurat_obj, layer = source_layer, assay = assay_name)

  cat(sprintf("  Matrix: %d genes × %d cells\n", nrow(mat), ncol(mat)))
  cat(sprintf(
    "  Current genes: %s\n",
    paste(head(rownames(mat), 5), collapse = ", ")
  ))

  # Check if already symbols
  if (!grepl("^ENSG", rownames(mat)[1])) {
    cat("\n✓ Already using gene symbols\n")
    return(mat)
  }

  # Convert gene names
  cat("\nStep 2: Converting ENSG to symbols...\n")

  ensg_ids <- rownames(mat)
  ensg_clean <- gsub("\\..*$", "", ensg_ids) # Remove version

  symbol_map <- mapIds(
    org.Hs.eg.db,
    keys = ensg_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  new_symbols <- ifelse(is.na(symbol_map), ensg_ids, as.character(symbol_map))

  n_converted <- sum(!is.na(symbol_map))
  cat(sprintf(
    "  Converted: %d / %d (%.1f%%)\n",
    n_converted,
    length(ensg_ids),
    n_converted / length(ensg_ids) * 100
  ))

  # Handle duplicates
  n_dup <- sum(duplicated(new_symbols))
  if (n_dup > 0) {
    cat(sprintf(
      "\n  Found %d duplicates, handling: %s\n",
      n_dup,
      handle_duplicates
    ))

    if (handle_duplicates == "unique") {
      new_symbols <- make.unique(new_symbols, sep = "_")
    } else if (handle_duplicates == "sum") {
      cat("  Aggregating duplicates by summing...\n")
      unique_symbols <- unique(new_symbols)
      new_mat <- do.call(
        rbind,
        lapply(unique_symbols, function(sym) {
          idx <- which(new_symbols == sym)
          if (length(idx) == 1) {
            mat[idx, , drop = FALSE]
          } else {
            Matrix::colSums(mat[idx, , drop = FALSE])
          }
        })
      )
      rownames(new_mat) <- unique_symbols
      mat <- new_mat
      new_symbols <- unique_symbols
      cat(sprintf("  Aggregated to %d unique genes\n", length(unique_symbols)))
    } else if (handle_duplicates == "first") {
      keep <- !duplicated(new_symbols)
      mat <- mat[keep, ]
      new_symbols <- new_symbols[keep]
    }
  }

  # Update rownames
  if (handle_duplicates != "sum") {
    rownames(mat) <- new_symbols
  }

  cat(sprintf("\n✓ Conversion complete\n"))
  cat(sprintf("  Final genes: %d\n", nrow(mat)))
  cat(sprintf("  Sample: %s\n", paste(head(rownames(mat), 5), collapse = ", ")))

  return(mat)
}


# Step 3: Add converted matrix to new layer
add_converted_layer <- function(
  seurat_obj,
  converted_matrix,
  new_layer_name = "counts.symbol"
) {
  cat("\n=== ADDING NEW LAYER ===\n\n")

  assay_name <- DefaultAssay(seurat_obj)

  cat(sprintf("New layer name: %s\n", new_layer_name))
  cat(sprintf(
    "Matrix: %d genes × %d cells\n",
    nrow(converted_matrix),
    ncol(converted_matrix)
  ))

  # Add as new layer
  seurat_obj <- SetAssayData(
    seurat_obj,
    layer = new_layer_name,
    new.data = converted_matrix,
    assay = assay_name
  )

  cat(sprintf("✓ Layer '%s' added\n", new_layer_name))

  # Verify
  cat("\nVerification:\n")
  layers_after <- Layers(seurat_obj, assay = assay_name)
  cat(sprintf("  Available layers: %s\n", paste(layers_after, collapse = ", ")))

  # Check the new layer
  check_mat <- LayerData(seurat_obj, layer = new_layer_name, assay = assay_name)
  cat(sprintf(
    "  New layer genes: %s\n",
    paste(head(rownames(check_mat), 5), collapse = ", ")
  ))

  return(seurat_obj)
}


# ==============================================================================
# All-in-one function
# ==============================================================================

simple_convert <- function(
  seurat_obj,
  source_layer = "counts",
  new_layer_name = "counts.symbol",
  handle_duplicates = "unique"
) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Simple Layer-Based Gene Conversion\n")
  cat(rep("=", 70), "\n", sep = "")

  # Step 1: Check layers
  cat("\n[1/3] Checking layer structure...\n")
  check_layers(seurat_obj)

  # Step 2: Convert
  cat("\n[2/3] Converting gene names...\n")
  converted_mat <- convert_layer_genes(
    seurat_obj,
    source_layer = source_layer,
    handle_duplicates = handle_duplicates
  )

  # Step 3: Add new layer
  cat("\n[3/3] Adding converted layer...\n")
  seurat_obj <- add_converted_layer(
    seurat_obj,
    converted_mat,
    new_layer_name = new_layer_name
  )

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("✓ COMPLETE\n")
  cat(rep("=", 70), "\n\n", sep = "")

  cat("Next steps:\n")
  cat(sprintf(
    "  1. Check new layer: LayerData(seurat_obj, layer = '%s')\n",
    new_layer_name
  ))
  cat(sprintf(
    "  2. View genes: head(rownames(LayerData(seurat_obj, layer = '%s')))\n",
    new_layer_name
  ))
  cat("  3. Use for analysis or set as default\n\n")

  return(seurat_obj)
}


# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# # Step-by-step approach:
#
# # 1. Check what layers exist
# check_layers(seurat_obj1)
#
# # 2. Convert the counts layer
# converted_counts <- convert_layer_genes(seurat_obj1, source_layer = "counts")
#
# # 3. Add as new layer
# seurat_obj1 <- add_converted_layer(seurat_obj1, converted_counts, "counts.symbol")
#
# # 4. Verify
# LayerData(seurat_obj1, layer = "counts.symbol") %>% rownames() %>% head(20)

# # One-command approach:
# seurat_obj1 <- simple_convert(
#   seurat_obj1,
#   source_layer = "counts",
#   new_layer_name = "counts.symbol",
#   handle_duplicates = "sum"  # or "unique" or "first"
# )

# # Quick wrapper with defaults:
# quick_layer_convert <- function(seurat_obj) {
#   simple_convert(
#     seurat_obj,
#     source_layer = "counts",
#     new_layer_name = "counts.symbol",
#     handle_duplicates = "sum"
#   )
# }
#
# # Usage:
# seurat_obj1 <- quick_layer_convert(seurat_obj1)

# ==============================================================================
# Debug: Check what's happening to gene names
# ==============================================================================

debug_gene_names <- function(seurat_obj, converted_matrix) {
  cat("\n=== DEBUGGING GENE NAMES ===\n\n")

  # Check converted matrix
  cat("Converted matrix:\n")
  cat(sprintf(
    "  Dimensions: %d × %d\n",
    nrow(converted_matrix),
    ncol(converted_matrix)
  ))
  cat(sprintf("  Rownames class: %s\n", class(rownames(converted_matrix))))
  cat(sprintf(
    "  First 10 genes: %s\n",
    paste(head(rownames(converted_matrix), 10), collapse = ", ")
  ))
  cat(sprintf(
    "  Gene format: %s\n",
    ifelse(grepl("^ENSG", rownames(converted_matrix)[1]), "ENSG", "Symbols")
  ))

  # Check object's current features
  cat("\nSeurat object features:\n")
  cat(sprintf("  Total features: %d\n", nrow(seurat_obj)))
  cat(sprintf(
    "  First 10: %s\n",
    paste(head(rownames(seurat_obj), 10), collapse = ", ")
  ))

  # Check if they match
  cat("\nComparison:\n")
  cat(sprintf(
    "  Matrix genes == Object genes: %s\n",
    identical(rownames(converted_matrix), rownames(seurat_obj))
  ))
  cat(sprintf(
    "  Matrix dim == Object dim: %s\n",
    nrow(converted_matrix) == nrow(seurat_obj)
  ))
}


# ==============================================================================
# Fixed: Add layer with gene name preservation (v2)
# ==============================================================================

add_converted_layer_v2 <- function(
  seurat_obj,
  converted_matrix,
  new_layer_name = "counts.symbol",
  create_new_assay = FALSE
) {
  cat("\n=== ADDING CONVERTED LAYER (v2) ===\n\n")

  assay_name <- DefaultAssay(seurat_obj)

  cat(sprintf(
    "Strategy: %s\n",
    ifelse(create_new_assay, "Create new assay", "Add as layer")
  ))
  cat(sprintf(
    "Matrix: %d genes × %d cells\n",
    nrow(converted_matrix),
    ncol(converted_matrix)
  ))
  cat(sprintf(
    "Gene format: %s\n",
    paste(head(rownames(converted_matrix), 5), collapse = ", ")
  ))

  if (create_new_assay) {
    # Strategy 1: Create a completely new assay with symbols
    cat("\nCreating new assay with gene symbols...\n")

    new_assay_name <- paste0(assay_name, ".symbol")

    new_assay <- CreateAssayObject(counts = converted_matrix)
    seurat_obj[[new_assay_name]] <- new_assay

    cat(sprintf("✓ New assay '%s' created\n", new_assay_name))
    cat(sprintf(
      "  Genes: %s\n",
      paste(head(rownames(seurat_obj[[new_assay_name]]), 5), collapse = ", ")
    ))

    # Verify
    check_genes <- rownames(seurat_obj[[new_assay_name]])
    if (grepl("^ENSG", check_genes[1])) {
      cat("⚠️  WARNING: Genes still ENSG format\n")
    } else {
      cat("✓ SUCCESS: Genes are symbols\n")
    }
  } else {
    # Strategy 2: Try to add as layer (may lose gene names in v5)
    cat("\nAttempting to add as layer...\n")
    cat("⚠️  Note: Seurat v5 may reset gene names to match object features\n\n")

    # Try direct approach
    seurat_obj <- SetAssayData(
      seurat_obj,
      layer = new_layer_name,
      new.data = converted_matrix,
      assay = assay_name
    )

    # Check if it worked
    check_mat <- LayerData(
      seurat_obj,
      layer = new_layer_name,
      assay = assay_name
    )
    check_genes <- rownames(check_mat)

    if (grepl("^ENSG", check_genes[1])) {
      cat("❌ Layer method failed - genes reset to ENSG\n")
      cat("\nTrying alternative: Create new assay instead...\n")

      new_assay_name <- paste0(assay_name, ".symbol")
      new_assay <- CreateAssayObject(counts = converted_matrix)
      seurat_obj[[new_assay_name]] <- new_assay

      cat(sprintf("✓ Created assay '%s' instead\n", new_assay_name))
      cat(sprintf("  Use: DefaultAssay(seurat_obj) <- '%s'\n", new_assay_name))
    } else {
      cat(sprintf("✓ Layer '%s' added successfully\n", new_layer_name))
      cat(sprintf(
        "  Genes: %s\n",
        paste(head(check_genes, 5), collapse = ", ")
      ))
    }
  }

  # Final verification
  cat("\nFinal status:\n")
  cat(sprintf(
    "  Available assays: %s\n",
    paste(names(seurat_obj@assays), collapse = ", ")
  ))
  cat(sprintf("  Current default: %s\n", DefaultAssay(seurat_obj)))

  return(seurat_obj)
}


# ==============================================================================
# Recommended: Create new assay with symbols (Most reliable)
# ==============================================================================

convert_to_symbol_assay <- function(
  seurat_obj,
  source_layer = "counts",
  new_assay_name = "RNA.symbol",
  handle_duplicates = "sum"
) {
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install: BiocManager::install('org.Hs.eg.db')")
  }
  library(AnnotationDbi)

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Convert to Symbol Assay (Recommended Method)\n")
  cat(rep("=", 70), "\n", sep = "")

  assay_name <- DefaultAssay(seurat_obj)

  # Extract matrix
  cat("\n[1/4] Extracting counts...\n")
  mat <- LayerData(seurat_obj, layer = source_layer, assay = assay_name)
  cat(sprintf("  %d genes × %d cells\n", nrow(mat), ncol(mat)))

  # Convert
  cat("\n[2/4] Converting gene names...\n")
  ensg_clean <- gsub("\\..*$", "", rownames(mat))

  symbol_map <- mapIds(
    org.Hs.eg.db,
    keys = ensg_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  new_symbols <- ifelse(
    is.na(symbol_map),
    rownames(mat),
    as.character(symbol_map)
  )

  n_conv <- sum(!is.na(symbol_map))
  cat(sprintf(
    "  Converted: %d / %d (%.1f%%)\n",
    n_conv,
    length(ensg_clean),
    n_conv / length(ensg_clean) * 100
  ))

  # Handle duplicates
  if (sum(duplicated(new_symbols)) > 0) {
    cat(sprintf("\n[3/4] Handling duplicates (%s)...\n", handle_duplicates))

    if (handle_duplicates == "sum") {
      unique_symbols <- unique(new_symbols)
      new_mat <- do.call(
        rbind,
        lapply(unique_symbols, function(sym) {
          idx <- which(new_symbols == sym)
          if (length(idx) == 1) {
            mat[idx, , drop = FALSE]
          } else {
            Matrix::colSums(mat[idx, , drop = FALSE])
          }
        })
      )
      rownames(new_mat) <- unique_symbols
      mat <- new_mat
      cat(sprintf("  Aggregated to %d unique genes\n", nrow(mat)))
    } else if (handle_duplicates == "unique") {
      rownames(mat) <- make.unique(new_symbols, sep = "_")
    } else {
      keep <- !duplicated(new_symbols)
      mat <- mat[keep, ]
    }
  } else {
    rownames(mat) <- new_symbols
    cat("\n[3/4] No duplicates found\n")
  }

  # Create new assay
  cat(sprintf("\n[4/4] Creating assay '%s'...\n", new_assay_name))
  new_assay <- CreateAssayObject(counts = mat)
  seurat_obj[[new_assay_name]] <- new_assay

  cat(sprintf("✓ Assay created\n"))
  cat(sprintf(
    "  Genes: %s\n",
    paste(head(rownames(seurat_obj[[new_assay_name]]), 5), collapse = ", ")
  ))

  # Verification
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Verification\n")
  cat(rep("=", 70), "\n", sep = "")

  check_genes <- rownames(seurat_obj[[new_assay_name]])
  cat(sprintf("Assay '%s':\n", new_assay_name))
  cat(sprintf("  Total genes: %d\n", length(check_genes)))
  cat(sprintf(
    "  Format: %s\n",
    ifelse(grepl("^ENSG", check_genes[1]), "ENSG (FAILED)", "Symbols (SUCCESS)")
  ))
  cat(sprintf("  Sample: %s\n", paste(head(check_genes, 10), collapse = ", ")))

  cat("\nUsage:\n")
  cat(sprintf("  DefaultAssay(seurat_obj) <- '%s'\n", new_assay_name))
  cat("  Then use normally: FeaturePlot(), DotPlot(), etc.\n")

  cat(rep("=", 70), "\n\n", sep = "")

  return(seurat_obj)
}


# ==============================================================================
# USAGE
# ==============================================================================

# # Method 1: Debug current situation
# debug_gene_names(seurat_obj1, converted_counts)
#
# # Method 2: Try to add as layer (may not preserve names in v5)
# seurat_obj1 <- add_converted_layer_v2(
#   seurat_obj1,
#   converted_counts,
#   create_new_assay = FALSE
# )
#
# # Method 3: Create new assay (RECOMMENDED - always works)
# seurat_obj1 <- add_converted_layer_v2(
#   seurat_obj1,
#   converted_counts,
#   create_new_assay = TRUE
# )
#
# Method 4: All-in-one (EASIEST)
# seurat_obj1 <- convert_to_symbol_assay(
#   seurat_obj1,
#   source_layer = "counts",
#   new_assay_name = "RNA.symbol",
#   handle_duplicates = "unique"
# )
#
# # Then use the symbol assay
# DefaultAssay(seurat_obj1) <- "RNA.symbol"
# FeaturePlot(seurat_obj1, features = c("CD3D", "CD4"))

# ==============================================================================
# Clean Up Old Layers/Assays
# ==============================================================================

# Remove specific layer from an assay
remove_layer <- function(seurat_obj, layer_name, assay_name = NULL) {
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(seurat_obj)
  }

  cat(sprintf(
    "\nRemoving layer '%s' from assay '%s'...\n",
    layer_name,
    assay_name
  ))

  # Check if layer exists
  current_layers <- Layers(seurat_obj, assay = assay_name)

  if (!layer_name %in% current_layers) {
    cat(sprintf("⚠️  Layer '%s' not found\n", layer_name))
    cat(sprintf(
      "Available layers: %s\n",
      paste(current_layers, collapse = ", ")
    ))
    return(seurat_obj)
  }

  # Remove layer by setting to NULL
  seurat_obj[[assay_name]][[layer_name]] <- NULL

  # Verify
  new_layers <- Layers(seurat_obj, assay = assay_name)

  if (layer_name %in% new_layers) {
    cat("❌ Failed to remove layer\n")
  } else {
    cat(sprintf("✓ Layer '%s' removed\n", layer_name))
    cat(sprintf("Remaining layers: %s\n", paste(new_layers, collapse = ", ")))
  }

  return(seurat_obj)
}


# Remove entire assay
remove_assay <- function(seurat_obj, assay_name) {
  cat(sprintf("\nRemoving assay '%s'...\n", assay_name))

  # Check if exists
  if (!assay_name %in% names(seurat_obj@assays)) {
    cat(sprintf("⚠️  Assay '%s' not found\n", assay_name))
    cat(sprintf(
      "Available assays: %s\n",
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
    return(seurat_obj)
  }

  # Don't remove if it's the only assay
  if (length(seurat_obj@assays) == 1) {
    stop("Cannot remove the only assay in object")
  }

  # Don't remove if it's the default
  if (assay_name == DefaultAssay(seurat_obj)) {
    cat("⚠️  This is the default assay\n")
    cat("Switching default to first available assay...\n")
    other_assays <- setdiff(names(seurat_obj@assays), assay_name)
    DefaultAssay(seurat_obj) <- other_assays[1]
    cat(sprintf("New default: %s\n", DefaultAssay(seurat_obj)))
  }

  # Remove
  seurat_obj[[assay_name]] <- NULL

  # Verify
  if (assay_name %in% names(seurat_obj@assays)) {
    cat("❌ Failed to remove assay\n")
  } else {
    cat(sprintf("✓ Assay '%s' removed\n", assay_name))
    cat(sprintf(
      "Remaining assays: %s\n",
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
  }

  return(seurat_obj)
}


# Clean all layers except specified ones
keep_only_layers <- function(
  seurat_obj,
  keep_layers = "counts",
  assay_name = NULL
) {
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(seurat_obj)
  }

  cat(sprintf("\nCleaning assay '%s'...\n", assay_name))
  cat(sprintf("Keeping only: %s\n", paste(keep_layers, collapse = ", ")))

  current_layers <- Layers(seurat_obj, assay = assay_name)
  layers_to_remove <- setdiff(current_layers, keep_layers)

  if (length(layers_to_remove) == 0) {
    cat("✓ Already clean, no layers to remove\n")
    return(seurat_obj)
  }

  cat(sprintf("Removing: %s\n", paste(layers_to_remove, collapse = ", ")))

  for (layer in layers_to_remove) {
    seurat_obj[[assay_name]][[layer]] <- NULL
  }

  # Verify
  final_layers <- Layers(seurat_obj, assay = assay_name)
  cat(sprintf("✓ Final layers: %s\n", paste(final_layers, collapse = ", ")))

  return(seurat_obj)
}


# ==============================================================================
# All-in-one: Convert and clean up in one command
# ==============================================================================

convert_and_cleanup <- function(
  seurat_obj,
  source_layer = "counts",
  new_assay_name = "RNA.symbol",
  handle_duplicates = "sum",
  remove_old_assay = TRUE,
  set_as_default = TRUE
) {
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install: BiocManager::install('org.Hs.eg.db')")
  }
  library(AnnotationDbi)

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Convert to Symbols and Clean Up\n")
  cat(rep("=", 70), "\n", sep = "")

  old_assay_name <- DefaultAssay(seurat_obj)

  # Step 1: Convert
  cat("\n[1/4] Converting to symbols...\n")

  mat <- LayerData(seurat_obj, layer = source_layer, assay = old_assay_name)
  cat(sprintf("  Source: %d genes × %d cells\n", nrow(mat), ncol(mat)))

  ensg_clean <- gsub("\\..*$", "", rownames(mat))

  symbol_map <- mapIds(
    org.Hs.eg.db,
    keys = ensg_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  new_symbols <- ifelse(
    is.na(symbol_map),
    rownames(mat),
    as.character(symbol_map)
  )

  n_conv <- sum(!is.na(symbol_map))
  cat(sprintf(
    "  Converted: %d / %d (%.1f%%)\n",
    n_conv,
    length(ensg_clean),
    n_conv / length(ensg_clean) * 100
  ))

  # Handle duplicates
  if (sum(duplicated(new_symbols)) > 0) {
    cat(sprintf("  Handling duplicates: %s\n", handle_duplicates))

    if (handle_duplicates == "sum") {
      unique_symbols <- unique(new_symbols)
      new_mat <- do.call(
        rbind,
        lapply(unique_symbols, function(sym) {
          idx <- which(new_symbols == sym)
          if (length(idx) == 1) {
            mat[idx, , drop = FALSE]
          } else {
            Matrix::colSums(mat[idx, , drop = FALSE])
          }
        })
      )
      rownames(new_mat) <- unique_symbols
      mat <- new_mat
    } else if (handle_duplicates == "unique") {
      rownames(mat) <- make.unique(new_symbols, sep = "_")
    } else {
      keep <- !duplicated(new_symbols)
      mat <- mat[keep, ]
    }
  } else {
    rownames(mat) <- new_symbols
  }

  # Step 2: Create new assay
  cat("\n[2/4] Creating new assay...\n")
  new_assay <- CreateAssayObject(counts = mat)
  seurat_obj[[new_assay_name]] <- new_assay
  cat(sprintf(
    "  ✓ Assay '%s' created (%d genes)\n",
    new_assay_name,
    nrow(seurat_obj[[new_assay_name]])
  ))

  # Step 3: Set as default
  if (set_as_default) {
    cat("\n[3/4] Setting as default assay...\n")
    DefaultAssay(seurat_obj) <- new_assay_name
    cat(sprintf("  ✓ Default assay: %s\n", DefaultAssay(seurat_obj)))
  } else {
    cat("\n[3/4] Keeping original default assay\n")
  }

  # Step 4: Remove old assay
  if (remove_old_assay) {
    cat("\n[4/4] Removing old assay...\n")
    cat(sprintf("  Removing: %s\n", old_assay_name))
    seurat_obj[[old_assay_name]] <- NULL
    cat("  ✓ Old assay removed\n")
  } else {
    cat("\n[4/4] Keeping old assay\n")
  }

  # Final status
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Final Status\n")
  cat(rep("=", 70), "\n", sep = "")

  cat(sprintf(
    "Object: %d genes × %d cells\n",
    nrow(seurat_obj),
    ncol(seurat_obj)
  ))
  cat(sprintf(
    "Available assays: %s\n",
    paste(names(seurat_obj@assays), collapse = ", ")
  ))
  cat(sprintf("Default assay: %s\n", DefaultAssay(seurat_obj)))
  cat(sprintf(
    "Gene format: %s\n",
    paste(head(rownames(seurat_obj), 5), collapse = ", ")
  ))

  cat(rep("=", 70), "\n\n", sep = "")

  return(seurat_obj)
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# # Example 1: Remove specific layer
# seurat_obj1 <- remove_layer(seurat_obj1, "counts.symbol")
#
# # Example 2: Remove old RNA assay (after creating RNA.symbol)
# seurat_obj1 <- remove_assay(seurat_obj1, "RNA")
#
# # Example 3: Keep only counts layer in RNA assay
# seurat_obj1 <- keep_only_layers(seurat_obj1, keep_layers = "counts", assay_name = "RNA")
#
# # Example 4: ONE COMMAND - Convert and remove old data
# seurat_obj1 <- convert_and_cleanup(
#   seurat_obj1,
#   source_layer = "counts",
#   new_assay_name = "RNA.symbol",
#   handle_duplicates = "sum",
#   remove_old_assay = TRUE,   # Remove old RNA assay
#   set_as_default = TRUE       # Set RNA.symbol as default
# )
#
# # Example 5: Convert but keep both assays
# seurat_obj1 <- convert_and_cleanup(
#   seurat_obj1,
#   remove_old_assay = FALSE    # Keep both RNA and RNA.symbol
# )
