process_cells <- function(
    seurat_object,                   # 输入的Seurat对象
    cell_subset = c("CD4_TM", "CD4_TN"), # 要分析的细胞类型
    annotation_column = "Annotation_2",  # 包含细胞类型注释的列名
    output_prefix = "CD4_Tcells",       # 输出文件的前缀
    output_dir = "results",              # 输出文件夹路径
    pca_dims = 1:40,                    # 用于后续分析的主成分数量
    umap_dims = 1:30,                   # 用于UMAP的主成分数量
    harmony_vars = c("sample"),         # 用于Harmony批次校正的变量
    harmony_theta = 5,                  # Harmony参数theta
    harmony_lambda = 1,                 # Harmony参数lambda
    harmony_sigma = 0.07,               # Harmony参数sigma
    harmony_nclust = 30,                # Harmony聚类数量
    umap_n_neighbors = 20,              # UMAP邻居数量
    umap_min_dist = 0.4,                # UMAP最小距离
    n_hvg = 4500,                       # 高变基因数量
    neighbor_k = 15,                    # FindNeighbors的k参数
    cluster_resolution = 3,             # 聚类分辨率
    cluster_algorithm = 4,              # 聚类算法
    marker_genes = NULL,                # 用于可视化的标记基因列表
    output_plots = TRUE,                # 是否生成并保存可视化结果
    plot_group_vars = c("study", "sample", "ann_level_3", "tissue", "ann_finest_level") # 用于分组可视化的变量
) {
  # 载入必要的包
  require(Seurat)
  require(harmony)
  require(ggplot2)
  require(dplyr)
  
  # 1. 细胞亚群提取
  message("提取", paste(cell_subset, collapse=", "), "细胞亚群...")
  subset_obj <- subset(seurat_object, get(annotation_column) %in% cell_subset)
  
  # 2. 标准数据预处理
  message("进行数据标准化...")
  subset_obj <- NormalizeData(subset_obj)
  
  # 3. 寻找高变异基因
  message("识别高变基因...")
  subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = n_hvg)
  
  # 4. 数据缩放
  message("数据缩放...")
  subset_obj <- ScaleData(subset_obj, features = VariableFeatures(subset_obj))
  
  # 5. 主成分分析
  message("执行主成分分析...")
  subset_obj <- RunPCA(subset_obj)
  
  # 6. Harmony批次效应校正
  message("使用Harmony进行批次效应校正...")
  subset_obj <- RunHarmony(
    object = subset_obj,           
    group.by.vars = harmony_vars,   
    theta = harmony_theta,                    
    lambda = harmony_lambda,                  
    sigma = harmony_sigma,           
    nclust = harmony_nclust,            
    reduction.use = "pca",
    max_iter = 20, 
    early_stop = TRUE,
    dims = pca_dims           
  )
  
  # 7. 运行UMAP降维
  message("执行UMAP降维...")
  subset_obj <- RunUMAP(subset_obj, 
                        reduction = "harmony", 
                        dims = umap_dims, 
                        n.neighbors = umap_n_neighbors,
                        min.dist = umap_min_dist,
                        metric = "correlation")
  
  # 8. 寻找邻居
  message("构建邻居网络...")
  subset_obj <- FindNeighbors(subset_obj, 
                              reduction = "harmony", 
                              dims = umap_dims,
                              k.param = neighbor_k) 
  
  # 9. 细胞聚类
  message("进行细胞聚类分析...")
  subset_obj <- FindClusters(subset_obj,
                             algorithm = cluster_algorithm,
                             group.singletons = FALSE,
                             resolution = cluster_resolution,
                             verbose = TRUE)
  
  # 设置默认的聚类结果
  Idents(subset_obj) <- subset_obj$seurat_clusters
  
  # 10. 如果需要，生成并保存可视化结果
  if(output_plots) {
    # 创建输出目录结构
    # 主输出目录
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message(paste("创建输出目录:", output_dir))
    }
    
    # 创建子目录
    plots_dir <- file.path(output_dir, "plots")
    data_dir <- file.path(output_dir, "data")
    cluster_dir <- file.path(plots_dir, "clusters")
    feature_dir <- file.path(plots_dir, "features")
    dotplot_dir <- file.path(plots_dir, "dotplots")
    annotation_dir <- file.path(plots_dir, "annotations")
    
    # 创建所有子目录
    for(dir_path in c(plots_dir, data_dir, cluster_dir, feature_dir, dotplot_dir, annotation_dir)) {
      if(!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
        message(paste("创建子目录:", dir_path))
      }
    }
    
    # 聚类可视化
    message("生成聚类可视化...")
    cluster_plot <- DimPlot(subset_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
      ggtitle(paste0(output_prefix, "细胞聚类")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    
    # 保存聚类UMAP图
    cluster_file <- file.path(cluster_dir, paste0(output_prefix, "_umap_clusters.pdf"))
    pdf(cluster_file, width = 10, height = 8)
    print(cluster_plot)
    dev.off()
    message(paste("保存聚类图至:", cluster_file))
    
    # 如果提供了marker_genes，则绘制特征图
    if(!is.null(marker_genes)) {
      message("生成标记基因特征图...")
      feature_file <- file.path(feature_dir, paste0(output_prefix, "_FeaturePlot.pdf"))
      pdf(feature_file, width = 8, height = 8)
      for(marker in marker_genes) {
        if(marker %in% rownames(subset_obj)) {
          message(paste("处理标记基因:", marker))
          print(FeaturePlot(subset_obj, features = marker, raster = TRUE))
        } else {
          message(paste("标记基因未找到:", marker))
        }
      }
      dev.off()
      message(paste("保存特征图至:", feature_file))
      
      # 生成点状图
      message("生成点状图...")
      dotplot_file <- file.path(dotplot_dir, paste0(output_prefix, "_markers_dotplot.pdf"))
      pdf(dotplot_file, width = 32, height = 18)
      dot_plot <- DotPlot(subset_obj, 
                          features = marker_genes, 
                          group.by = "seurat_clusters",
                          split.by = NULL,
                          cols = c("lightgrey", "red"),
                          dot.scale = 8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(dot_plot)
      dev.off()
      message(paste("保存点状图至:", dotplot_file))
      
      # 保存点状图数据
      source_data <- dot_plot$data
      dotplot_data_file <- file.path(data_dir, paste0(output_prefix, "_dotplot_data.csv"))
      write.csv(source_data, dotplot_data_file)
      message(paste("保存点状图数据至:", dotplot_data_file))
    }
    
    # 分组可视化
    message("生成分组可视化...")
    annotation_file <- file.path(annotation_dir, paste0(output_prefix, "_umap_annotation.pdf"))
    pdf(annotation_file, width = 15, height = 10)
    for(group_var in plot_group_vars) {
      if(group_var %in% colnames(subset_obj@meta.data)) {
        p <- DimPlot(subset_obj, reduction = "umap", group.by = group_var, raster = TRUE, 
                     pt.size = 0.5) + ggtitle(group_var)
        print(p)
      } else {
        message(paste("分组变量未找到:", group_var))
      }
    }
    dev.off()
    message(paste("保存分组可视化至:", annotation_file))
  }
  
  # 返回处理后的Seurat对象
  return(subset_obj)
}

# 使用示例
# 假设T_object是已准备好的包含所有T细胞的Seurat对象，Marker_T是标记基因列表
# CD4_obj <- process_cells(
#   seurat_object = T_object,
#   cell_subset = c("CD4_TM", "CD4_TN"),
#   annotation_column = "Annotation_2",
#   output_prefix = "CD4_Tcells",
#   output_dir = "results/CD4_analysis", # 指定输出目录
#   marker_genes = Marker_T
# )
#
# # 保存处理后的Seurat对象
# saveRDS(CD4_obj, file.path("results/CD4_analysis", "CD4_Tcells_processed.rds"))