# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(rsvd)
library(harmony)
library(presto)
Sys.setenv(RETICULATE_PYTHON = "C:/Users/崔填祎/AppData/Local/Programs/Python/Python311/python.exe")
library(reticulate)
library(DoubletFinder)
library(DirichletReg)
library(speckle)
library(tidyverse)
library(MAST)
require(DirichletReg)
require(speckle)
require(ggplot2)
require(cowplot)
require(edgeR)
require(reshape2)
require(pheatmap)
require(tidyr)
library(patchwork)
options(future.globals.maxSize = Inf )
library(Seurat)
library(dplyr)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(foreach)
library(doParallel)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(visNetwork)
library(DESeq2)
source("E:/R/scMASC.R")
source("E:/R/addgrids3d.r")
Markers <- c('FXYD3', 'EPCAM', 'ELF3', 'SERPINF1', 'TSPAN1',
             'SCGB1A1', 'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8',#T
             'CD53', 'PTPRC', 'CORO1A', 'CCL5',
             'MS4A1', 'TNFRSF17', 'CD19', 'CD79A', #B
             'CD40LG', 'TNFRSF25', 'CD28', 'CD4', 'CD3E', 'CD8A', 'CD8B', 
             'TRGC2', 'CD2', 'TRBC2', #T
             'FCER1G', 'C1orf162', 'CLEC7A', 'CD1C', 'CD86', 'CD14', 'XCR1', 'HLA-DRA',#Myeloid
             'COL1A2', 'DCN', 'MFAP4', 'LUM', 'COL6A3', 'CFD', 'COL1A1', 'PDGFRA', 
             'MXRA8', 'LEPR',
             'MYH11', 'TINAGL1', 'PLN', 'DES', 'ACTA2', 'CNN1', 'TAGLN', #Fibroblast & SMC
             'CLDN5', 'ECSCR', 'CLEC14A', 'VWF', 'PECAM1', 'DARC', 'PTPRB', 'PDE2A', 
             'PLAT', 'GJA5', 'SPARCL1', 'AQP1', 'MMRN1', 'CCL21', #Endothelial
             'MKI67', 'TOP2A', 'TK1', 'CENPW')#Proliferation

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
              'KLRB1','IL7R','NCR3','CEBPD','SLC4A10','TRAV1-2','NCR3',
              'FOXP3','IL2RA','IKZF2','TNFRSF4',
              'MKI67','TOP2A','TK1','CENPW')

setwd("E:/R/0301/T_analysis/")#设置工作路径
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
T_object <- readRDS("T_analyzed.rds")#读取数据
T_object <- NormalizeData(T_object)#归一化
T_object <- FindVariableFeatures(T_object, selection.method = "vst", nfeatures = 5000)#寻找变异基因
T_object <- ScaleData(T_object)#标准化
T_object <- RunPCA(T_object, npcs = 50)#PCA
T_object <- RunHarmony(
    object = T_object,           
    group.by.vars = c("sample"),   
    theta = c(3),
    lambda = c(0.8),
    sigma = 0.05,           
    # nclust = 60,            
    reduction.use = "pca",
    max_iter = 20,
    # cool_down = 10,
    # epsilon_harmony = 0.0001,
    # monitor = "harmony_score",
    # factor = 0.8,
    # tol = 1e-10,
    # min_delta = 0.0001,
    # patience = 10,
    # epsilon_cluster = 0.00001,
    early_stop = TRUE,
    dims = 1:50           
)
# 获取降维结果和元数据
embeddings <- Embeddings(T_object, "harmony")
metadata <- T_object@meta.data

# 过滤掉ann_level_3为NA的细胞
valid_cells <- !is.na(metadata$ann_level_3)
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
n_iterations <- 10
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
                                   c("sample", "ann_level_3"), 
                                   k)
    
    # 存储结果
    ilisi_scores[i] <- mean(lisi_res[, "sample"])
    clisi_scores[i] <- mean(lisi_res[, "ann_level_3"])
    
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

# # 创建输出目录（如果不存在）
# if(!dir.exists("qc_plots")) {
#   dir.create("qc_plots")
# }

# # 保存所有knee plots到一个PDF文件
# pdf("qc_plots/knee_plots_analysis.pdf", width = 12, height = 10)

# # 1. UMI knee plot
# counts <- GetAssayData(T_object, layer = "counts")
# total_umi <- Matrix::colSums(counts)
# umi_rank <- rank(-total_umi)

# print(ggplot(data.frame(rank = umi_rank, UMI = total_umi), 
#              aes(x = rank, y = UMI)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         scale_x_log10() +
#         scale_y_log10() +
#         theme_bw() +
#         labs(title = "Post-integration UMI Knee Plot",
#              x = "Cell rank",
#              y = "Total UMIs"))

# # 2. 基因数量的knee plot
# genes_per_cell <- Matrix::colSums(counts > 0)
# gene_rank <- rank(-genes_per_cell)

# print(ggplot(data.frame(rank = gene_rank, Genes = genes_per_cell), 
#              aes(x = rank, y = Genes)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         scale_x_log10() +
#         scale_y_log10() +
#         theme_bw() +
#         labs(title = "Post-integration Gene Detection Knee Plot",
#              x = "Cell rank",
#              y = "Number of genes detected"))

# # 3. PC肘部图
# pct_var <- T_object[["pca"]]@stdev / sum(T_object[["pca"]]@stdev) * 100
# cumsum_var <- cumsum(pct_var)

# elbow_data <- data.frame(
#   PC = seq_along(pct_var),
#   variance = pct_var,
#   cumulative = cumsum_var
# )

# print(ggplot(elbow_data, aes(x = PC)) +
#         geom_line(aes(y = variance), color = "blue") +
#         geom_point(aes(y = variance), color = "blue") +
#         geom_line(aes(y = cumulative), color = "red") +
#         geom_point(aes(y = cumulative), color = "red") +
#         theme_bw() +
#         labs(title = "PCA Elbow Plot",
#              x = "Principal Component",
#              y = "Percentage of Variance Explained") +
#         scale_y_continuous(
#           name = "Individual Variance Explained (%)",
#           sec.axis = sec_axis(~., name = "Cumulative Variance Explained (%)")
#         ))

# # 4. Harmony scores knee plot - 使用采样
# set.seed(42)
# n_sample <- min(10000, ncol(T_object))
# sample_idx <- sample(seq_len(ncol(T_object)), n_sample)

# harmony_scores <- T_object[["harmony"]]@cell.embeddings[sample_idx, 1:20]
# harmony_dist <- dist(harmony_scores)
# harmony_scores_ordered <- sort(colMeans(as.matrix(harmony_dist)), decreasing = TRUE)

# print(ggplot(data.frame(rank = seq_along(harmony_scores_ordered),
#                         score = harmony_scores_ordered), 
#              aes(x = rank, y = score)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         theme_bw() +
#         labs(title = "Harmony Scores Distribution (10k cells sampled)",
#              x = "Cell rank",
#              y = "Average Harmony distance"))

# dev.off()

# # 保存统计信息到CSV
# stats_df <- data.frame(
#   Metric = c("Total_cells", "Median_UMI", "Median_genes", 
#              "PC_90percent_var", "Harmony_mean_dist"),
#   Value = c(
#     length(total_umi),
#     median(total_umi),
#     median(genes_per_cell),
#     which(cumsum_var > 90)[1],
#     mean(harmony_scores_ordered)
#   )
# )

# write.csv(stats_df, "qc_plots/knee_plot_statistics.csv", row.names = FALSE)

# 直接运行UMAP
T_object <- RunUMAP(
    T_object, 
    reduction = "harmony", 
    dims = 1:30, 
    n.neighbors = 20,
    min.dist = 0.3,
    # learning.rate = 0.2,
    # n.trees = 500,
    # n.epochs = 800,
    # spread = 1.2,
    # repulsion.strength = 1.1,
    metric = "correlation"
)

# FindNeighbors函数的缩进修正
T_object <- FindNeighbors(
    T_object, 
    reduction = "harmony", 
    dims = 1:30,
        k.param = 15
) 

T_object <- FindClusters(
    T_object, 
    resolution = 2, 
    group.singletons = FALSE,
    algorithm = 4
    )
# 8. 绘制UMAP聚类图
print("Generating plots...")
Idents(T_object) <- T_object$seurat_clusters
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(T_object, reduction = "umap", label = TRUE)
dev.off()
# saveRDS(T_object, "final_analyzed_T_object.rds")
# 9. 保存分析结果
print("Saving analysis results...")

# T_object <- readRDS("E:/R/0224/final_analyzed_T_object.rds")


# 9. 标记基因分析

print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for(marker in Marker_T) {
  if(marker %in% rownames(T_object)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(T_object, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(T_object, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(T_object, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(T_object, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(T_object, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(T_object, reduction = "umap", group.by = "ann_finest_level", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_finest_level")
p6 <- DimPlot(T_object, reduction = "umap", group.by = "Annotation_2", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Annotation_2")
p1
p2
p3
p4
p5
p6
dev.off()

tissues <- unique(T_object$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(T_object) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(T_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = T_object, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(T_object$Annotation_2)
pdf("Annotation_2_umap_split.pdf", width = 12, height = 12)
Idents(T_object) <- "Annotation_2"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(T_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = T_object, 
                                                    idents = sample_name ),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(T_object$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(T_object) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(T_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = T_object, 
                                                    idents = study_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(study_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

# 先设置identities
# 先正确设置identities，排除NA值
T_object$ann_level_3_no_na <- T_object$ann_level_3
T_object$ann_level_3_no_na[is.na(T_object$ann_level_3_no_na)] <- "Unknown"
Idents(T_object) <- "ann_level_3_no_na"

pdf("ann_level_3_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_level_3 <- unique(T_object$ann_level_3)
ann_level_3 <- ann_level_3[!is.na(ann_level_3)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_level_3) {
  print(DimPlot(T_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(T_object, 
                                                    idents = anno),
                cols = "grey",
                pt.size = 0.5,
                raster = TRUE,
                label = FALSE) +
          ggtitle(anno) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()



# 1. 首先检查 Markers 中是否有重复
print(length(Marker_T))
print(length(unique(Marker_T)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Marker_T)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Marker_T <- unique(Marker_T)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot.pdf", width = 32, height = 18)
dot_plot <- DotPlot(T_object, 
                    features = Marker_T, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- dot_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotplot_data.csv")
dot_plot
dev.off()
getwd()
# dev.off()

# 2. 计算每个cluster中每个基因的表达情况
print("Calculating expression statistics...")

# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
Idents(T_object) <- T_object$seurat_clusters
markers <- FindAllMarkers(T_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

# table(subset(T_object, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

# 11. 生成热图
print("Generating heatmap...")
top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
write.csv(top10_markers, "cluster_top10_markers.csv")

# 1. 先缩放marker基因
marker_genes <- unique(top10_markers$gene)
T_object <- ScaleData(T_object, features = marker_genes)

# 1. 将marker基因分成较小的批次
batch_size <- 25  # 每批50个基因
marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))

# 生成文件名
filename <- "marker_heatmap_batch.pdf"
pdf(filename, width = 32, height = 18)
# 2. 为每个批次生成热图
for(i in seq_along(marker_batches)) {
  print(paste("Processing batch", i, "of", length(marker_batches)))
  
  # 绘制热图
  print(DoHeatmap(T_object, 
                  features = marker_batches[[i]],
                  size = 24) + 
          NoLegend() +
          ggtitle(paste("Marker Genes Batch", i))+
          theme(axis.text.y = element_text(size = 36))
  )
  # dev.off()
}
dev.off()

# 1. 读取注释文件（无表头）
annotations <- read.csv("Annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation_2")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- T_object$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
T_object$Annotation_2 <- new_idents
T_object <- SetIdent(T_object, value = new_idents)

# 4. 设置颜色
cell_colors <- c(
  "Epithelial" = "#E41A1C",   
  "Fibroblast" = "#377EB8",    
  "T" = "#4DAF4A",             
  "Myeloid" = "#984EA3",       
  "Endothelial" = "#FF7F00",   
  "SMC" = "#FFFF33",           
  "B" = "#A65628",             
  "Proliferation" = "#F781BF"  
)

# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(T_object, 
        reduction = "umap",
        raster = TRUE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        # cols = cell_colors
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 6. 绘制分割版本
pdf("cell_types_umap_split.pdf", width = 12, height = 12)
# 获取所有细胞类型
cell_types <- unique(T_object$Annotation_2)
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 为每个细胞类型单独绘图
for(cell_type in cell_types) {
  # 创建文件名
  
  
  # 绘制该细胞类型的UMAP图
  
  print(DimPlot(T_object, 
                reduction = "umap",
                cells.highlight = WhichCells(T_object, idents = cell_type),
                cols.highlight = cell_colors[cell_type],
                cols = "grey",  # 其他细胞显示为灰色
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(cell_type) +
          theme(legend.text = element_text(size = 12)))
  
}
dev.off()

# table(T_object@meta.data$seurat_clusters[T_object@meta.data$Annotation=='Epithelial' & T_object@meta.data$ann_level_2 =='Myeloid'])

# 假设我们选择 cell_type 为 "T_cells" 的细胞
# 1. 提取Undefined细胞子集
# 从主Seurat对象中筛选出Annotation标注为"Undefined"的细胞
# 检查T_object是否存在
exists("T_object")

# 检查T_object的结构
str(T_object)

# 检查Annotation_2列的唯一值
unique(T_object$Annotation_2)
setwd('..')
dir.create('./ILCs/')
setwd('./ILCs/')
# 检查是否有拼写错误或空格问题
table(T_object$Annotation_2)
Undefined_obj <- subset(T_object, Annotation_2 %in% c("ILCs"))

# 2. 标准数据预处理
# 数据标准化，默认使用"LogNormalize"方法
Undefined_obj <- NormalizeData(Undefined_obj)

# 3. 寻找高变异基因
# 使用方差稳定化转换(VST)方法选择3000个高变基因，用于后续降维分析
Undefined_obj <- FindVariableFeatures(Undefined_obj, selection.method = "vst", nfeatures = 4500)

# 4. 数据缩放
# 对所有高变基因进行归一化处理，使其均值为0，方差为1
Undefined_obj <- ScaleData(Undefined_obj, features = VariableFeatures(Undefined_obj))

# 5. 主成分分析
# 使用高变基因进行主成分分析，降低数据维度
Undefined_obj <- RunPCA(Undefined_obj)
getwd()
# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
Undefined_obj <- RunHarmony(
  object = Undefined_obj,           
  group.by.vars = c("sample"),   
  theta = c(3),      # Higher theta for more diverse clustering               
  lambda = c(1),   # Higher lambda to reduce overcorrection                
  sigma = 0.05,           # Lower sigma for tighter clusters
  nclust = 30,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20, 
  early_stop = TRUE,
  dims = 1:40           # More iterations for better convergence
)
Undefined_obj <- RunUMAP(Undefined_obj, 
                         reduction = "harmony", 
                         dims = 1:30, 
                         n.neighbors = 20,
                         # n.trees = 500,
                         min.dist = 0.4,
                         # learning.rate = 0.2,    # 相对保守的学习率
                         # n.epochs = 1400,        # 增加迭代次数补偿较小的学习率
                         # spread = 1.2,
                         # repulsion.strength = 1.1,
                         
                         metric = "correlation")
# Find neighbors
Undefined_obj <- FindNeighbors(Undefined_obj, 
                               reduction = "harmony", 
                               dims = 1:30,
                               k.param = 15) 

Undefined_obj <- FindClusters(Undefined_obj,
                              algorithm = 4,
                              group.singletons = FALSE,
                              resolution = 3, # 多个分辨率
                              verbose = TRUE)

# 11. 可视化聚类结果
print("生成聚类可视化...")
# 设置默认使用分辨率为2的聚类结果
Idents(Undefined_obj) <- Undefined_obj$seurat_clusters

# 保存UMAP聚类图
pdf("Undefined_umap_clusters.pdf", width = 10, height = 8)
DimPlot(Undefined_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Undefined细胞聚类") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
dev.off()

print("Generating feature plots...")
pdf("umap_FeaturePlot_Undefined.pdf", width = 8, height = 8)
for(marker in Marker_T) {
  if(marker %in% rownames(Undefined_obj)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Undefined_obj, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 14. 可视化
pdf("umap_anno_integration_Undefined_obj.pdf", width = 15, height = 10)
p1 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "ann_finest_level", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_finest_level")
p6 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "Annotation_2", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Annotation_2")
p1
p2
p3
p4
p5
p6
dev.off()

# 1. 首先检查 Markers 中是否有重复
print(length(Marker_T))
print(length(unique(Marker_T)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Marker_T)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Marker_T <- unique(Marker_T)
# Ident(T_object) <- T_object$Annotation_2
# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot_1.pdf", width = 32, height = 18)
dot_plot <- DotPlot(Undefined_obj, 
                    features = Marker_T, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- dot_plot$data
# 查看或保存数据
# View(source_data)
write.csv(source_data, "dotplot_data_undefined.csv")
dot_plot
dev.off()

# 1. 读取注释文件（无表头）
annotations <- read.csv("Annotation_Undefined_1.csv", header = FALSE)
# 正确命名列名为Annotation_2
colnames(annotations) <- c("Cluster", "Annotation_2")  # 使用Annotation_2作为列名

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- Undefined_obj$seurat_clusters  # 获取当前的cluster标识
# 使用Annotation_2列进行匹配
new_idents <- annotations$Annotation_2[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
# 使用Annotation_2作为列名添加到对象
Undefined_obj$Annotation_2 <- new_idents
# 使用SetIdent函数设置标识
Undefined_obj <- Seurat::SetIdent(Undefined_obj, value = "Annotation_2")  # 明确使用Annotation_2列

# 4. 设置颜色
cell_colors <- c(
  "Epithelial" = "#E41A1C",   
  "Fibroblast" = "#377EB8",    
  "T" = "#4DAF4A",             
  "Myeloid" = "#984EA3",       
  "Endothelial" = "#FF7F00",   
  "SMC" = "#FFFF33",           
  "B" = "#A65628",             
  "Proliferation" = "#F781BF"  
)

# 5. 绘制UMAP图
pdf("Undefined_cell_types_umap.pdf", width = 8, height = 8)
Seurat::DimPlot(Undefined_obj, 
                reduction = "umap",
                raster = TRUE,
                label = TRUE,
                pt.size = 0.5,
                label.size = 4,
                # cols = cell_colors  # 启用颜色设置
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 获取Undefined_obj中细胞的索引
cells_to_update <- rownames(Undefined_obj@meta.data)

# 检查这些细胞是否都存在于主对象中
cells_in_main <- cells_to_update %in% rownames(T_object@meta.data)
valid_cells_to_update <- cells_to_update[cells_in_main]

# 只更新主对象中存在的细胞的Annotation_2
T_object@meta.data[valid_cells_to_update, "Annotation_2"] <- 
  Undefined_obj@meta.data[valid_cells_to_update, "Annotation_2"]

# 验证更新
print(table(T_object@meta.data$Annotation_2))
T_object@meta.data$Annotation <- T_object@meta.data$Annotation_2
# 设置T_object的标识为Annotation_2
Seurat::Idents(T_object) <- "Annotation_2"

# 绘制更新后的UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
Seurat::DimPlot(T_object, 
                reduction = "umap",
                raster = TRUE,
                label = TRUE,
                pt.size = 0.5,
                label.size = 4,
                # cols = cell_colors  # 启用颜色设置
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 释放内存
rm(Undefined_obj)
gc()

# 输出最终细胞注释结果
print(table(T_object$Annotation_2))

# 将"Epthelial"替换为"Epithelial"
T_object$Annotation_2[T_object$Annotation_2 == "NK_cells"] <- "NK"

# 函数：生成高级UMAP可视化(带置信椭圆和坐标轴箭头)
generate_advanced_umap <- function(seurat_obj) {
  # 创建输出目录
  dir.create("10-Advanced_Visualizations", showWarnings = FALSE)
  setwd("10-Advanced_Visualizations")
  message("创建高级UMAP可视化...")
  
  # 整理绘图数据
  mydata <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)  # 取出降维数据
  mydata <- cbind(mydata, as.data.frame(seurat_obj@meta.data))  # 与注释合并
  
  # 模拟坐标轴数据
  umap_range <- apply(mydata[, c("UMAP_1", "UMAP_2")], 2, range)
  min_x <- umap_range[1, 1]
  min_y <- umap_range[1, 2]
  
  # 坐标轴位置
  axis_margin <- 0.5  # 边距
  axis_length <- 2    # 轴长度
  
  line.x.data <- data.frame(
    x = c(min_x - axis_margin, min_x - axis_margin + axis_length),
    y = c(min_y - axis_margin, min_y - axis_margin)
  )
  
  line.y.data <- data.frame(
    x = c(min_x - axis_margin, min_x - axis_margin),
    y = c(min_y - axis_margin, min_y - axis_margin + axis_length)
  )
  
  # 设置配色方案
  mycol <- c(pal_d3()(10), pal_aaas()(7), pal_uchicago()(7), pal_jama()(7))
  
  # 获取分组信息（如果有）
  if("group" %in% colnames(mydata)) {
    group_var <- "group"
  } else if("condition" %in% colnames(mydata)) {
    group_var <- "condition"
  } else {
    group_var <- NULL
    message("未找到group或condition变量，UMAP不会按组分面")
  }
  
  # 绘制基本图片
  p_umap <- ggplot(mydata, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = 0.8, aes(col = Annotation_2)) +
    scale_color_manual(values = mycol) + # 设置点的颜色
    labs(color = "Cell Type")
  
  # 尝试添加置信椭圆（如果每个组有足够的点）
  tryCatch({
    p_umap <- p_umap + 
      stat_ellipse(aes(fill = Annotation_2),
                   geom = "polygon",
                   linetype = 2,        # 置信区间的类型
                   linewidth = 0.7,     # 置信区间的粗细
                   alpha = 0.2) +
      scale_fill_manual(values = mycol)  # 让置信区间与点的颜色相一致
  }, error = function(e) {
    message("无法添加置信椭圆，可能是某些组的点数太少: ", e$message)
  })
  
  # 添加坐标轴箭头
  p_umap <- p_umap +
    geom_line(data = line.x.data,
              aes(x = x, y = y), 
              arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +  # 绘制X轴坐标轴
    geom_line(data = line.y.data,
              aes(x = x, y = y), 
              arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +  # 绘制Y轴坐标轴
    theme(
      panel.border = element_blank(),  # 隐藏边框
      axis.title = element_blank(),    # 隐藏轴标题
      axis.text = element_blank(),     # 隐藏文本
      axis.ticks = element_blank(),    # 隐藏轴线
      panel.background = element_rect(fill = 'white'),  # 背景色
      plot.background = element_rect(fill = "white"),
      panel.grid = element_blank()     # 去除网格线
    )
  
  # 如果有分组变量，按组分面
  if(!is.null(group_var) && length(unique(mydata[[group_var]])) > 1) {
    p_umap <- p_umap + facet_grid(~ get(group_var))
  }
  
  # 保存高级UMAP图
  ggsave("advanced_umap.pdf", plot = p_umap, width = 12, height = 10)
  
  # 计算并添加细胞比例图
  message("创建细胞比例图...")
  
  # 整理细胞比例数据
  if(!is.null(group_var)) {
    cell_props <- as.data.frame(table(mydata$Annotation_2, mydata[[group_var]]))
    colnames(cell_props) <- c("CellType", "Group", "Freq")
    
    # 计算每组内的比例
    cell_props <- cell_props %>%
      group_by(Group) %>%
      mutate(Total = sum(Freq),
             Proportion = Freq / Total) %>%
      ungroup()
    
    # 绘制堆叠柱状图
    p_prop <- ggplot(cell_props, aes(x = Group, y = Proportion, fill = CellType)) +
      geom_bar(stat = "identity", position = "fill") + 
      scale_fill_manual(values = mycol) +
      labs(title = "细胞类型比例", x = "", y = "比例") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
    
    # 保存细胞比例图
    ggsave("cell_proportions.pdf", plot = p_prop, width = 8, height = 6)
  }
  
  # 返回可视化对象
  return(list(
    umap_plot = p_umap,
    prop_plot = if(exists("p_prop")) p_prop else NULL
  ))
}

# 函数：生成标记基因表达的高级可视化
generate_advanced_marker_viz <- function(seurat_obj, top_n = 5) {
  setwd("10-Advanced_Visualizations")
  message("创建标记基因高级可视化...")
  
  # 寻找每个细胞类型的标记基因
  Idents(seurat_obj) <- "Annotation_2"
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 为每个细胞类型选择top N个标记基因
  top_markers <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = top_n, wt = avg_log2FC)
  
  # 保存标记基因结果
  write.csv(markers, file = "all_markers.csv")
  write.csv(top_markers, file = "top_markers.csv")
  
  # 创建高级散点图可视化
  # 准备数据
  markers$gene_rank <- rank(-markers$avg_log2FC)
  
  # 高亮顶部基因
  top_genes <- unique(top_markers$gene)
  markers$highlight <- markers$gene %in% top_genes
  
  # 创建散点图
  p_scatter <- ggplot(markers, 
                      aes(x = pct.1 - pct.2, y = avg_log2FC)) +
    geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
    scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) +
    scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1)) +
    geom_text_repel(data = subset(markers, highlight == TRUE),
                    aes(label = gene),
                    box.padding = 0.5,
                    point.padding = 0.3,
                    force = 10,
                    segment.color = "grey50",
                    size = 3) +
    facet_wrap(~cluster, scales = "free_y") +
    labs(x = "表达细胞百分比差异 (pct.1 - pct.2)",
         y = "平均Log2倍变化",
         title = "各细胞类型的标记基因") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # 保存散点图
  ggsave("marker_genes_scatter.pdf", plot = p_scatter, width = 14, height = 10)
  
  # 绘制热图
  # 从标记基因中提取顶部基因表达数据
  top_genes_union <- unique(top_markers$gene)
  
  # 计算每个细胞类型中的平均表达
  avg_expr <- AverageExpression(seurat_obj, 
                                features = top_genes_union, 
                                assays = "RNA",
                                group.by = "Annotation_2",
                                slot = "data")
  
  # 绘制热图
  heatmap_data <- avg_expr$RNA
  
  # 对数据进行Z分数标准化
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # 为热图准备注释
  gene_cluster <- markers %>%
    filter(gene %in% rownames(heatmap_data_scaled)) %>%
    select(gene, cluster) %>%
    distinct()
  
  gene_anno <- data.frame(
    Cluster = gene_cluster$cluster[match(rownames(heatmap_data_scaled), gene_cluster$gene)]
  )
  rownames(gene_anno) <- rownames(heatmap_data_scaled)
  
  # 绘制热图
  pdf("marker_genes_heatmap.pdf", width = 12, height = 14)
  pheatmap(heatmap_data_scaled,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_row = gene_anno,
           fontsize_row = 8,
           fontsize_col = 10,
           main = "标记基因表达热图")
  dev.off()
  
  return(list(
    markers = markers,
    top_markers = top_markers,
    scatter_plot = p_scatter
  ))
}


saveRDS(T_object,"T_analyzed.RDS")

source("E:/R/enrichment_functions.R")
source('E:/R/scMASC.R')
results <- scPairwiseMASCAnalysis(
  seurat_obj = T_object,
  cell_type_col = "Annotation_2",
  sample_col = "sample",
  contrast_col = "tissue",
  # exclude_filter = "tissue_sampling_method == 'scraping'",
  fixed_effects_cols = NULL,  # 添加任何固定效应协变量
  output_dir = "pairwise_MASC_results",
  # 如果您需要少于2个样本也能进行比较，可以降低此阈值
  min_samples = 2
)

# 加载所需的额外包
library(GSVA)
library(msigdbr)
library(pheatmap)
library(DESeq2)
library(reshape2)
library(ggplot2)
library(tibble)
library(qs)

# 假设T_object已经完成了前面的处理步骤
# 创建目录保存分析结果
dir.create("7-pseudobulk_analysis", showWarnings = FALSE)
dir.create("8-GSVA_analysis", showWarnings = FALSE)

########################
# 第一部分：Pseudobulk分析
########################
# =====================================================================
# Pseudobulk分析: 所有组织间的细胞类型两两比较
# =====================================================================

########################
# 第一部分：改进版Pseudobulk分析 - 带Top基因标注火山图
########################

# 创建目录保存分析结果
dir.create("7-pseudobulk_analysis_enhanced", showWarnings = FALSE)
setwd("7-pseudobulk_analysis_enhanced")
message("开始进行增强版组织间细胞类型两两比较分析...")

# 确保必需包已加载
required_packages <- c("DESeq2", "ggplot2", "ggrepel", "dplyr", "pheatmap")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("安装", pkg, "包..."))
    if(pkg == "DESeq2") {
      if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("DESeq2")
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# 检查T_object必要列
required_columns <- c("Annotation_2", "tissue", "sample")
missing_columns <- required_columns[!required_columns %in% colnames(T_object@meta.data)]
if(length(missing_columns) > 0) {
  stop(paste("T_object缺少必要的元数据列:", paste(missing_columns, collapse=", ")))
}

# 检查组织和细胞类型分布
message("分析不同组织中的细胞类型分布...")
tissue_cell_distribution <- table(T_object$tissue, T_object$Annotation_2)
write.csv(tissue_cell_distribution, "tissue_celltype_distribution.csv")

message("数据集包含以下组织类型:")
print(table(T_object$tissue))
message("数据集包含以下T细胞亚型:")
print(table(T_object$Annotation_2))

# 1. Pseudobulk聚合
message("执行基于组织和细胞类型的pseudobulk聚合...")
av_tissue_celltype <- AggregateExpression(
  T_object,
  group.by = c("tissue", "sample", "Annotation_2"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

av_tissue_celltype_df <- as.data.frame(av_tissue_celltype[[1]])
write.csv(av_tissue_celltype_df, file = "pseudobulk_tissue_sample_celltype.csv")

# 2. 提取元数据
extract_metadata <- function(column_names) {
  parts <- strsplit(column_names, "_")
  result <- data.frame(
    column = column_names,
    tissue = sapply(parts, function(x) x[1]),
    sample = sapply(parts, function(x) x[2]),
    celltype = sapply(parts, function(x) paste(x[-(1:2)], collapse="_")),
    stringsAsFactors = FALSE
  )
  return(result)
}

metadata <- extract_metadata(colnames(av_tissue_celltype_df))
rownames(metadata) <- metadata$column

message("提取的元数据示例(前5行):")
print(head(metadata, 5))

# ============================================================================
# 3. 改进版两两比较函数 - 带Top10基因标注火山图
# ============================================================================
run_enhanced_pairwise_comparison <- function(counts_matrix, metadata, cell_type, 
                                             top_n_genes = 10) {
  
  # 筛选特定细胞类型的样本
  cell_indices <- which(metadata$celltype == cell_type)
  
  if(length(cell_indices) < 6) {
    message(paste("细胞类型", cell_type, "的样本数量不足(", length(cell_indices), ")，跳过分析"))
    return(NULL)
  }
  
  # 提取该细胞类型的数据
  cell_counts <- counts_matrix[, cell_indices, drop = FALSE]
  cell_metadata <- metadata[cell_indices, , drop = FALSE]
  
  # 查看组织分布
  tissue_counts <- table(cell_metadata$tissue)
  message(paste("细胞类型", cell_type, "在各组织中的样本数:"))
  print(tissue_counts)
  
  # 只保留至少有3个样本的组织
  valid_tissues <- names(tissue_counts[tissue_counts >= 3])
  
  if(length(valid_tissues) < 2) {
    message(paste("细胞类型", cell_type, "中至少有3个样本的组织不足2种，跳过比较"))
    return(NULL)
  }
  
  # 过滤数据
  valid_indices <- cell_metadata$tissue %in% valid_tissues
  cell_counts <- cell_counts[, valid_indices, drop = FALSE]
  cell_metadata <- cell_metadata[valid_indices, , drop = FALSE]
  cell_counts <- round(cell_counts)
  
  # 生成所有可能的两两组织比较
  tissue_pairs <- combn(valid_tissues, 2, simplify = FALSE)
  message(paste("为细胞类型", cell_type, "执行", length(tissue_pairs), "个组织间两两比较"))
  
  results_list <- list()
  
  # 对每对组织进行比较
  for(pair in tissue_pairs) {
    tissue1 <- pair[1]
    tissue2 <- pair[2]
    
    message(paste("比较:", tissue1, "vs", tissue2))
    
    # 选择这两个组织的样本
    pair_indices <- cell_metadata$tissue %in% pair
    pair_counts <- cell_counts[, pair_indices, drop = FALSE]
    pair_metadata <- cell_metadata[pair_indices, , drop = FALSE]
    
    # 检查样本数
    tissue_sample_counts <- table(pair_metadata$tissue)
    if(any(tissue_sample_counts < 3)) {
      message(paste("跳过比较: 某个组织样本数少于3:", 
                    paste(names(tissue_sample_counts), tissue_sample_counts, sep="=", collapse=", ")))
      next
    }
    
    # DESeq2分析
    dds <- DESeqDataSetFromMatrix(
      countData = pair_counts,
      colData = pair_metadata,
      design = ~ tissue
    )
    
    dds$tissue <- relevel(factor(dds$tissue), ref = tissue1)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    dds <- DESeq(dds)
    
    # 获取差异结果
    comparison_name <- paste0(tissue2, "_vs_", tissue1)
    res <- results(dds, contrast = c("tissue", tissue2, tissue1))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[order(res_df$padj), ]
    res_df <- na.omit(res_df)
    
    # 添加上下调信息
    res_df$regulation <- ifelse(res_df$padj > 0.05, "stable",
                                ifelse(abs(res_df$log2FoldChange) < 1, "stable",
                                       ifelse(res_df$log2FoldChange >= 1, "up", "down")))
    
    # 保存结果
    results_list[[comparison_name]] <- res_df
    
    # 创建输出目录
    result_dir <- paste0(cell_type, "_", tissue1, "_vs_", tissue2)
    result_dir <- gsub("[^a-zA-Z0-9_]", "-", result_dir)
    dir.create(result_dir, showWarnings = FALSE)
    
    # 保存完整DEG列表
    write.csv(res_df, file = file.path(result_dir, "DEGs_complete.csv"), row.names = FALSE)
    
    # ========================================================================
    # 识别Top差异基因
    # ========================================================================
    # 提取显著差异基因
    sig_genes <- res_df[res_df$regulation != "stable", ]
    
    if(nrow(sig_genes) > 0) {
      # 按p值和log2FC综合排序选择top基因
      # 计算综合得分: -log10(padj) * |log2FC|
      sig_genes$score <- -log10(sig_genes$padj) * abs(sig_genes$log2FoldChange)
      sig_genes <- sig_genes[order(-sig_genes$score), ]
      
      # 分别选择上调和下调的top基因
      up_genes <- sig_genes[sig_genes$regulation == "up", ]
      down_genes <- sig_genes[sig_genes$regulation == "down", ]
      
      # 选择top N/2上调和N/2下调，如果某一类不足则从另一类补充
      n_up <- min(ceiling(top_n_genes/2), nrow(up_genes))
      n_down <- min(ceiling(top_n_genes/2), nrow(down_genes))
      
      # 如果一类不足，从另一类补充
      if(n_up < ceiling(top_n_genes/2) && nrow(down_genes) > n_down) {
        n_down <- min(top_n_genes - n_up, nrow(down_genes))
      } else if(n_down < ceiling(top_n_genes/2) && nrow(up_genes) > n_up) {
        n_up <- min(top_n_genes - n_down, nrow(up_genes))
      }
      
      top_genes_list <- c()
      if(n_up > 0) top_genes_list <- c(top_genes_list, up_genes$gene[1:n_up])
      if(n_down > 0) top_genes_list <- c(top_genes_list, down_genes$gene[1:n_down])
      
      # 保存top基因列表及其统计信息
      if(length(top_genes_list) > 0) {
        top_genes_info <- res_df[res_df$gene %in% top_genes_list, 
                                 c("gene", "baseMean", "log2FoldChange", "padj", "regulation")]
        top_genes_info <- top_genes_info[order(-abs(top_genes_info$log2FoldChange)), ]
        
        write.csv(top_genes_info, 
                  file = file.path(result_dir, "Top_genes_annotated.csv"), 
                  row.names = FALSE)
        
        message(paste("识别到", length(top_genes_list), "个top差异基因"))
      }
    } else {
      top_genes_list <- character(0)
      message("未发现显著差异基因")
    }
    
    # ========================================================================
    # 创建增强版火山图（带Top基因标注）
    # ========================================================================
    volcano_file <- file.path(result_dir, "volcano_plot_annotated.pdf")
    pdf(volcano_file, width = 12, height = 10)
    
    # 准备标注数据
    res_df$label <- ""
    if(length(top_genes_list) > 0) {
      res_df$label[res_df$gene %in% top_genes_list] <- res_df$gene[res_df$gene %in% top_genes_list]
    }
    
    # 创建火山图
    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      # 底层所有点
      geom_point(aes(color = regulation), size = 2, alpha = 0.6) +
      
      # 配色方案
      scale_color_manual(
        values = c("down" = "#00468B", "stable" = "grey70", "up" = "#E64B35"),
        labels = c("down" = paste0("Down (", sum(res_df$regulation == "down"), ")"),
                   "stable" = paste0("Stable (", sum(res_df$regulation == "stable"), ")"),
                   "up" = paste0("Up (", sum(res_df$regulation == "up"), ")"))
      ) +
      
      # 阈值线
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = 'black', linewidth = 0.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = 'black', linewidth = 0.5) +
      
      # 标题和标签
      labs(
        title = paste0(cell_type, ": ", tissue2, " vs ", tissue1),
        subtitle = paste0("Top ", length(top_genes_list), " differential genes annotated"),
        x = "Log2(Fold Change)",
        y = "-Log10(Adjusted P-value)",
        color = "Regulation"
      ) +
      
      # 主题设置
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12)
      )
    
    # 添加基因标注（如果有top基因）
    if(length(top_genes_list) > 0) {
      p <- p + 
        # 高亮top基因点
        geom_point(data = res_df[res_df$label != "", ], 
                   aes(x = log2FoldChange, y = -log10(padj)), 
                   size = 3, shape = 21, fill = NA, color = "black", stroke = 1.5) +
        
        # 添加基因名标注（使用ggrepel避免重叠）
        ggrepel::geom_text_repel(
          data = res_df[res_df$label != "", ],
          aes(label = label),
          size = 3.5,
          fontface = "bold",
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "black",
          segment.size = 0.5,
          max.overlaps = Inf,
          min.segment.length = 0,
          force = 2,
          force_pull = 0.5
        )
    }
    
    print(p)
    dev.off()
    
    message(paste0("增强版火山图已保存: ", volcano_file))
    
    # ========================================================================
    # 生成摘要统计
    # ========================================================================
    summary_stats <- data.frame(
      Comparison = comparison_name,
      CellType = cell_type,
      Total_DEGs = sum(res_df$regulation != "stable"),
      Up_regulated = sum(res_df$regulation == "up"),
      Down_regulated = sum(res_df$regulation == "down"),
      Total_genes_tested = nrow(res_df),
      Top_genes_count = length(top_genes_list),
      stringsAsFactors = FALSE
    )
    
    write.csv(summary_stats, 
              file = file.path(result_dir, "summary_stats.csv"), 
              row.names = FALSE)
    
    # ========================================================================
    # 可选：创建Top基因热图
    # ========================================================================
    if(length(top_genes_list) > 0 && length(top_genes_list) <= 50) {
      tryCatch({
        # 提取top基因的表达值
        top_gene_expr <- pair_counts[top_genes_list, , drop = FALSE]
        
        # 标准化
        top_gene_expr_norm <- t(scale(t(log2(top_gene_expr + 1))))
        
        # 创建注释
        anno_col <- data.frame(
          Tissue = pair_metadata$tissue,
          row.names = colnames(top_gene_expr_norm)
        )
        
        # 添加基因的调控方向注释
        anno_row <- data.frame(
          Regulation = res_df$regulation[match(rownames(top_gene_expr_norm), res_df$gene)],
          row.names = rownames(top_gene_expr_norm)
        )
        
        # 绘制热图
        heatmap_file <- file.path(result_dir, "Top_genes_heatmap.pdf")
        pdf(heatmap_file, width = 10, height = max(8, length(top_genes_list) * 0.3))
        
        # ⭐ 修复：明确保存并打印pheatmap对象
        p_heatmap <- pheatmap(
          top_gene_expr_norm,
          annotation_col = anno_col,
          annotation_row = anno_row,
          annotation_colors = list(
            Tissue = setNames(c("#E64B35", "#00468B"), pair),
            Regulation = c("up" = "#E64B35", "down" = "#00468B")
          ),
          main = paste0("Top ", length(top_genes_list), " DEGs: ", tissue2, " vs ", tissue1),
          fontsize_row = 9,
          fontsize_col = 8,
          show_colnames = TRUE,
          cluster_cols = TRUE,
          cluster_rows = TRUE,
          color = colorRampPalette(c("navy", "white", "firebrick"))(100),
          silent = FALSE  # ⭐ 确保不静默
        )
        
        # ⭐ 关键：明确打印热图对象
        print(p_heatmap)
        
        dev.off()
        message(paste0("Top基因热图已保存: ", heatmap_file))
      }, error = function(e) {
        dev.off()  # 确保即使出错也关闭设备
        message(paste("热图生成失败:", e$message))
      })
    }
    
    # ========================================================================
    # 可选：GO富集分析（如果有clusterProfiler）
    # ========================================================================
    if(requireNamespace("clusterProfiler", quietly = TRUE) && 
       requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      
      message("执行GO富集分析...")
      library(clusterProfiler)
      library(org.Hs.eg.db)
      
      # 上调基因GO富集
      up_genes <- res_df$gene[res_df$regulation == "up"]
      if(length(up_genes) >= 10) {
        ego_up <- enrichGO(
          gene = up_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        if(!is.null(ego_up) && nrow(ego_up) > 0) {
          write.csv(as.data.frame(ego_up), 
                    file = file.path(result_dir, "GO_upregulated.csv"), 
                    row.names = FALSE)
          
          pdf(file.path(result_dir, "GO_upregulated_dotplot.pdf"), width = 10, height = 8)
          print(dotplot(ego_up, showCategory = 20, title = "GO Enrichment: Upregulated Genes"))
          dev.off()
        }
      }
      
      # 下调基因GO富集
      down_genes <- res_df$gene[res_df$regulation == "down"]
      if(length(down_genes) >= 10) {
        ego_down <- enrichGO(
          gene = down_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        if(!is.null(ego_down) && nrow(ego_down) > 0) {
          write.csv(as.data.frame(ego_down), 
                    file = file.path(result_dir, "GO_downregulated.csv"), 
                    row.names = FALSE)
          
          pdf(file.path(result_dir, "GO_downregulated_dotplot.pdf"), width = 10, height = 8)
          print(dotplot(ego_down, showCategory = 20, title = "GO Enrichment: Downregulated Genes"))
          dev.off()
        }
      }
    }
  }
  
  return(results_list)
}

# ============================================================================
# 4. 执行增强版分析
# ============================================================================
all_cell_types <- unique(metadata$celltype)
results_by_celltype <- list()

for(cell_type in all_cell_types) {
  message(paste("\n========================================"))
  message(paste("分析细胞类型:", cell_type))
  message(paste("========================================"))
  
  if(is.na(cell_type) || cell_type == "") {
    message("跳过无效细胞类型")
    next
  }
  
  # 创建细胞类型目录
  cell_dir <- gsub("[^a-zA-Z0-9]", "_", cell_type)
  dir.create(cell_dir, showWarnings = FALSE)
  
  original_dir <- getwd()
  setwd(cell_dir)
  
  # 执行增强版分析
  tryCatch({
    results <- run_enhanced_pairwise_comparison(
      counts_matrix = av_tissue_celltype_df,
      metadata = metadata,
      cell_type = cell_type,
      top_n_genes = 10  # 可调整标注基因数量
    )
    
    if(!is.null(results)) {
      results_by_celltype[[cell_type]] <- results
    }
  }, error = function(e) {
    message(paste("分析细胞类型", cell_type, "时出错:", e$message))
  })
  
  setwd(original_dir)
}

# ============================================================================
# 5. 创建交叉比较总结报告
# ============================================================================
message("\n创建交叉比较总结报告...")

# 收集所有比较的Top基因
all_top_genes <- list()
all_comparisons <- data.frame()

for(cell_dir in list.dirs(recursive = FALSE, full.names = TRUE)) {
  cell_type <- basename(cell_dir)
  
  if(!dir.exists(cell_dir) || cell_dir %in% c(".", "..")) next
  
  comparison_dirs <- list.dirs(path = cell_dir, recursive = FALSE, full.names = TRUE)
  
  for(comp_dir in comparison_dirs) {
    # 读取统计摘要
    stats_file <- file.path(comp_dir, "summary_stats.csv")
    if(file.exists(stats_file)) {
      stats <- read.csv(stats_file)
      all_comparisons <- rbind(all_comparisons, stats)
    }
    
    # 收集Top基因
    top_genes_file <- file.path(comp_dir, "Top_genes_annotated.csv")
    if(file.exists(top_genes_file)) {
      top_genes <- read.csv(top_genes_file)
      comparison_name <- basename(comp_dir)
      all_top_genes[[comparison_name]] <- top_genes
    }
  }
}

# 保存汇总表
if(nrow(all_comparisons) > 0) {
  write.csv(all_comparisons, "all_comparisons_summary.csv", row.names = FALSE)
  
  # 创建Top基因频率统计
  if(length(all_top_genes) > 0) {
    # 统计每个基因在多少个比较中出现
    all_gene_names <- unlist(lapply(all_top_genes, function(x) x$gene))
    gene_frequency <- as.data.frame(table(all_gene_names))
    colnames(gene_frequency) <- c("Gene", "Frequency")
    gene_frequency <- gene_frequency[order(-gene_frequency$Frequency), ]
    
    write.csv(gene_frequency, "top_genes_frequency_across_comparisons.csv", row.names = FALSE)
    
    # 绘制基因频率图
    if(nrow(gene_frequency) > 0) {
      top_freq_genes <- head(gene_frequency, min(30, nrow(gene_frequency)))
      
      pdf("top_genes_frequency_barplot.pdf", width = 12, height = 8)
      p <- ggplot(top_freq_genes, aes(x = reorder(Gene, Frequency), y = Frequency)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
        coord_flip() +
        labs(
          title = "Most Frequently Identified Top DEGs Across All Comparisons",
          subtitle = paste("Based on", length(all_top_genes), "tissue-tissue comparisons"),
          x = "Gene Symbol",
          y = "Number of Comparisons"
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
        )
      print(p)
      dev.off()
    }
  }
  
  # 创建比较结果热图
  if(nrow(all_comparisons) > 1) {
    # DEG数量热图
    heatmap_data <- reshape2::dcast(all_comparisons, 
                                    CellType ~ Comparison, 
                                    value.var = "Total_DEGs")
    rownames(heatmap_data) <- heatmap_data$CellType
    heatmap_data$CellType <- NULL
    heatmap_data[is.na(heatmap_data)] <- 0
    
    if(ncol(heatmap_data) > 0 && nrow(heatmap_data) > 0) {
      pdf("comparison_DEG_counts_heatmap.pdf", width = 12, height = 10)
      pheatmap(
        heatmap_data,
        display_numbers = TRUE,
        number_format = "%.0f",
        fontsize_number = 8,
        main = "Number of DEGs Across Tissue Comparisons by Cell Type",
        color = colorRampPalette(c("white", "orange", "red"))(50),
        cluster_rows = TRUE,
        cluster_cols = TRUE
      )
      dev.off()
    }
  }
}

message("\n========================================")
message("增强版Pseudobulk分析完成！")
message("========================================")
message("主要改进:")
message("1. ✓ 火山图添加Top10差异基因标注")
message("2. ✓ 保存Top基因详细信息（Top_genes_annotated.csv）")
message("3. ✓ 创建Top基因表达热图")
message("4. ✓ 统计Top基因在所有比较中的出现频率")
message("5. ✓ 生成交叉比较可视化报告")

# 返回主目录
setwd("..")
########################
# 第二部分：GSVA分析
########################
getwd()
setwd("./8-GSVA_analysis")
message("开始进行GSVA分析...")
library(msigdbr)

# 设置工作目录
current_dir <- getwd()
output_dir <- "./8-GSVA_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# 加载必要的库包
library(Seurat)
library(msigdbr)
library(GSVA)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)

message("开始进行GSVA分析...")

run_gsva_analysis <- function(expr_data, gene_sets, method = "gsva", mx.diff = TRUE, output_prefix, min.sz = 7) {
  # 运行GSVA
  message(paste0("使用", output_prefix, "基因集运行GSVA..."))
  tryCatch({
    gsva_result <- gsva(expr_data, gene_sets, method = method, mx.diff = mx.diff, min.sz = min.sz)
    write.csv(gsva_result, file = paste0("gsva_", output_prefix, "_results.csv"))
    return(gsva_result)
  }, error = function(e) {
    message(paste0("GSVA分析出错: ", e$message))
    return(NULL)
  })
}

# 改进的热图绘制函数
plot_gsva_heatmap <- function(gsva_result, top_n = 20, output_prefix, title) {
  tryCatch({
    if(is.null(gsva_result) || nrow(gsva_result) == 0) {
      message(paste0("Cannot draw ", output_prefix, " heatmap: Results are empty"))
      return(NULL)
    }
    
    # 确保有足够多的行
    if(nrow(gsva_result) < 2) {
      message(paste0("Cannot draw ", output_prefix, " heatmap: Too few rows (", nrow(gsva_result), ")"))
      return(NULL)
    }
    
    # 计算行方差并检查
    row_vars <- apply(gsva_result, 1, var)
    if(all(is.na(row_vars)) || all(row_vars == 0)) {
      message("All rows have zero or NA variance, cannot create meaningful heatmap")
      # 创建一个简单的占位图以确保PDF不为空
      pdf(paste0(output_prefix, "_gsva_heatmap.pdf"), width = 7, height = 5)
      plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
      text(1.5, 1.5, "Cannot create heatmap: No variance in data", cex=1.2)
      dev.off()
      return(NULL)
    }
    
    # 选择变异最大的通路
    select_n <- min(top_n, sum(!is.na(row_vars) & row_vars > 0))
    if(select_n == 0) {
      message("No pathways with positive variance found")
      # 创建一个简单的占位图
      pdf(paste0(output_prefix, "_gsva_heatmap.pdf"), width = 7, height = 5)
      plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
      text(1.5, 1.5, "Cannot create heatmap: No pathways with variance", cex=1.2)
      dev.off()
      return(NULL)
    }
    
    var_pathways <- names(tail(sort(row_vars), select_n))
    message(paste0("Selected ", length(var_pathways), " pathways for ", output_prefix, " heatmap"))
    
    # 准备热图数据 - 确保数据有效
    heatmap_data <- gsva_result[var_pathways, , drop=FALSE]
    
    # 检查数据是否全为NA
    if(all(is.na(heatmap_data))) {
      message("Heatmap data contains only NA values")
      pdf(paste0(output_prefix, "_gsva_heatmap.pdf"), width = 7, height = 5)
      plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
      text(1.5, 1.5, "Cannot create heatmap: Data contains only NA values", cex=1.2)
      dev.off()
      return(NULL)
    }
    
    # 尝试使用标准R热图函数作为备选
    pdf_path <- paste0(output_prefix, "_gsva_heatmap.pdf")
    pdf(pdf_path, width = 10, height = max(4, ceiling(length(var_pathways)/3)))
    
    # 尝试pheatmap，如果失败则使用基础热图
    success <- FALSE
    
    tryCatch({
      # 设置颜色
      col_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
      
      # 尝试使用pheatmap
      hmap <- pheatmap(heatmap_data,
                       scale = "row",
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       fontsize_row = 8,
                       fontsize_col = 10,
                       angle_col = "45",
                       main = title,
                       color = col_palette)
      print(hmap)
      success <- TRUE
      message("Successfully created heatmap using pheatmap")
    }, error = function(e) {
      message(paste0("pheatmap error: ", e$message, ". Trying base heatmap instead."))
      
      # 如果pheatmap失败，尝试使用基础热图函数
      tryCatch({
        # 数据预处理
        data_for_heatmap <- as.matrix(heatmap_data)
        # 行标准化
        data_scaled <- t(scale(t(data_for_heatmap)))
        # 处理可能的Inf和NA
        data_scaled[is.infinite(data_scaled)] <- NA
        data_scaled[is.na(data_scaled)] <- 0
        
        # 绘制基础热图
        heatmap(data_scaled, 
                main = title,
                col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                scale = "none", # 已经标准化过
                margins = c(8, 8))
        success <- TRUE
        message("Successfully created heatmap using base heatmap function")
      }, error = function(e2) {
        message(paste0("Base heatmap also failed: ", e2$message))
        
        # 如果一切都失败，只画一个简单的图确保PDF不为空
        plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
        text(1.5, 1.5, "Heatmap creation failed. See R console for details.", cex=1.2)
      })
    })
    
    # 确保关闭PDF设备
    dev.off()
    
    # 如果两种方法都失败，至少创建一个包含错误消息的PDF
    if(!success) {
      pdf(paste0(output_prefix, "_gsva_heatmap_error.pdf"), width = 7, height = 5)
      plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
      text(1.5, 1.5, "Failed to create heatmap. Check data and R console.", cex=1.2)
      dev.off()
    }
    
    return(success)
  }, error = function(e) {
    message(paste0("Error in plot_gsva_heatmap function: ", e$message))
    # 确保创建一个包含错误信息的PDF
    pdf(paste0(output_prefix, "_gsva_heatmap_error.pdf"), width = 7, height = 5)
    plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
    text(1.5, 1.5, paste0("Error: ", e$message), cex=1.2)
    dev.off()
    return(FALSE)
  })
}

# Modified plot_gsva_bubble function
plot_gsva_bubble <- function(gsva_result, var_pathways, output_prefix, title) {
  if(is.null(gsva_result) || nrow(gsva_result) == 0) {
    message(paste0("Cannot draw ", output_prefix, " bubble plot: Results are empty"))
    return(NULL)
  }
  
  # Ensure var_pathways is not empty and exists in gsva_result
  if(length(var_pathways) == 0) {
    message(paste0("Cannot draw ", output_prefix, " bubble plot: Pathway list is empty"))
    return(NULL)
  }
  
  # Only keep pathways that exist in gsva_result
  valid_pathways <- intersect(var_pathways, rownames(gsva_result))
  if(length(valid_pathways) == 0) {
    message(paste0("Cannot draw ", output_prefix, " bubble plot: No valid pathways"))
    return(NULL)
  }
  
  message(paste0("Using ", length(valid_pathways), " pathways for ", output_prefix, " bubble plot"))
  
  # Safely draw bubble plot
  pdf_path <- paste0(output_prefix, "_gsva_bubble.pdf")
  pdf(pdf_path, width = 12, height = 10)
  
  tryCatch({
    # Prepare long format data
    gsva_data_long <- reshape2::melt(gsva_result[valid_pathways, ], 
                                     varnames = c("Pathway", "CellType"), 
                                     value.name = "GSVA_Score")
    
    # Draw bubble plot
    p <- ggplot(gsva_data_long, aes(x = CellType, y = Pathway, size = abs(GSVA_Score), color = GSVA_Score)) +
      geom_point(alpha = 0.8) +
      scale_size_continuous(range = c(1, 8), name = "Absolute Score") +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "GSVA Score") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      labs(title = title,
           x = "Cell Type", 
           y = "Pathway")
    # Force print plot
    print(p)
    message(paste0("Successfully drew bubble plot and saved to ", pdf_path))
    invisible(p)
  }, error = function(e) {
    message(paste0("Error drawing bubble plot: ", e$message))
    # Ensure a valid PDF even when error occurs
    plot(1, type="n", axes=FALSE, ann=FALSE)
    text(1, 1, paste0("Error drawing bubble plot: ", e$message), cex=1.2)
    NULL
  }, finally = {
    # Ensure PDF device is closed
    dev.off()
  })
}
# 2. 数据准备
# 检查Seurat对象是否在工作环境中
if(!exists("T_object")) {
  stop("在工作环境中找不到'T_object' Seurat对象。请先加载数据。")
}
# 获取所有基因集
all_gene_sets <- msigdbr(species = "Homo sapiens")

# 查看数据框结构
str(all_gene_sets)
colnames(all_gene_sets)

# 根据实际列名筛选KEGG相关通路
# 假设列名是'gs_name'和'gs_exact_source'或类似名称
# kegg_related <- unique(all_gene_sets[grep("KEGG", all_gene_sets$gs_name), 
#                                      c("gs_collection", "gs_subcat")])
# print(kegg_related)

# 3. 获取基因集
# 从MSigDB获取Hallmark和C2 KEGG基因集
message("获取基因集数据...")
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H") 
hallmark_genesets <- subset(hallmark_genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
hallmark_genesets <- split(hallmark_genesets$gene_symbol, hallmark_genesets$gs_name)
message(paste0("成功获取", length(hallmark_genesets), "个Hallmark基因集"))

# 获取C2集合中的KEGG通路
c2_kegg_genesets <- all_gene_sets[grep("KEGG", all_gene_sets$gs_name) & 
                                    all_gene_sets$gs_collection == "C2", ]
c2_kegg_genesets <- subset(c2_kegg_genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
c2_kegg_genesets <- split(c2_kegg_genesets$gene_symbol, c2_kegg_genesets$gs_name)

# 检查结果
# message(paste0("Successfully obtained ", length(c2_kegg_genesets), " KEGG gene sets"))
message(paste0("成功获取", length(c2_kegg_genesets), "个C2 KEGG基因集"))

# 检查组织注释是否存在
if(!"tissue" %in% colnames(T_object@meta.data)) {
  message("在Seurat对象中找不到'tissue'列。跳过组织特异性分析。")
} else {
  # 获取所有组织类型
  tissue_types <- unique(T_object$tissue)
  tissue_types <- tissue_types[!is.na(tissue_types)]
  message(paste0("共检测到", length(tissue_types), "个组织类型"))
}

# 在GSVA分析开始前
message("DEBUG: Starting GSVA analysis")
message("DEBUG: tissue_types = ", paste(tissue_types, collapse=", "))

# 初始化tissue_gsva_results列表
tissue_gsva_results <- list()
message("DEBUG: tissue_gsva_results initialized as empty list")

# 在每个组织分析前
for(tissue in tissue_types) {
  message("DEBUG: Processing tissue: ", tissue)
  
  # 直接从元数据中筛选细胞
  cells_in_tissue <- rownames(T_object@meta.data[T_object@meta.data$tissue == tissue, ])
  message("DEBUG: Number of cells = ", length(cells_in_tissue))
  
  if(length(cells_in_tissue) < 50) {
    message(paste0("组织'", tissue, "'的细胞数量太少(", length(cells_in_tissue), ")，跳过"))
    next
  }
  
  # 创建该组织的子集
  tissue_subset <- subset(T_object, cells = cells_in_tissue)
  
  # 设置细胞类型为标识
  Idents(tissue_subset) <- tissue_subset$Annotation_2
  
  # 计算平均表达
  expr_tissue <- AverageExpression(tissue_subset, assays = "RNA", slot = "data")[[1]]
  expr_tissue <- expr_tissue[rowSums(expr_tissue) > 0,]
  expr_tissue <- as.matrix(expr_tissue)
  message("DEBUG: Expression matrix dimensions: ", nrow(expr_tissue), " × ", ncol(expr_tissue))
  
  # 运行GSVA (只使用Hallmark基因集)
  message("DEBUG: Running GSVA for tissue ", tissue)
  gsva_tissue <- run_gsva_analysis(
    expr_data = expr_tissue, 
    gene_sets = hallmark_genesets, 
    output_prefix = gsub(" ", "_", tissue),
    method = "gsva",  # 显式指定方法
    mx.diff = TRUE    # 显式指定参数
  )
  
  # 检查GSVA结果
  if(!is.null(gsva_tissue)) {
    message("DEBUG: GSVA result dimensions: ", nrow(gsva_tissue), " × ", ncol(gsva_tissue))
    
    # 保存结果
    tissue_gsva_results[[tissue]] <- gsva_tissue
    message("DEBUG: Added GSVA result for tissue ", tissue, " to tissue_gsva_results")
    
    # 绘制组织特异性热图
    message("DEBUG: Calling plot_gsva_heatmap")
    plot_gsva_heatmap(
      gsva_result = gsva_tissue, 
      top_n = 20, 
      output_prefix = paste0("tissue_", gsub(" ", "_", tissue)), 
      title = paste0("Pathway Activity in ", tissue)
    )
  } else {
    message("DEBUG: GSVA result is NULL for tissue ", tissue)
  }
}

# 检查tissue_gsva_results
message("DEBUG: After processing all tissues, tissue_gsva_results length = ", length(tissue_gsva_results))
# 4. 整体细胞类型GSVA分析
# 确保细胞类型注释存在
if(!"Annotation_2" %in% colnames(T_object@meta.data)) {
  stop("在Seurat对象中找不到'Annotation_2'。请先进行细胞类型注释。")
}

# 设置细胞类型为标识
Idents(T_object) <- T_object$Annotation_2
message("按细胞类型计算平均表达...")

# 计算每个细胞类型的平均表达
expr <- AverageExpression(T_object, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0,]  # 过滤非表达基因
expr <- as.matrix(expr)
message(paste0("共获取", nrow(expr), "个基因在", ncol(expr), "个细胞类型中的平均表达"))

# 5. 运行GSVA分析
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

# 6. 可视化GSVA结果
# Hallmark基因集热图
if(!is.null(gsva_result_hallmark)) {
  # 绘制热图
  var_pathways_hallmark <- names(tail(sort(apply(gsva_result_hallmark, 1, var)), 20))
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
if(!is.null(gsva_result_kegg)) {
  # 绘制热图
  var_pathways_kegg <- names(tail(sort(apply(gsva_result_kegg, 1, var)), 30))
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

# 7. 按组织和细胞类型进行GSVA分析
message("按组织和细胞类型进行GSVA分析...")

# 检查组织注释是否存在
if(!"tissue" %in% colnames(T_object@meta.data)) {
  message("在Seurat对象中找不到'tissue'列。跳过组织特异性分析。")
} else {
  # 获取所有组织类型
  tissue_types <- unique(T_object$tissue)
  tissue_types <- tissue_types[!is.na(tissue_types)]
  message(paste0("共检测到", length(tissue_types), "个组织类型"))
}


# 按组织和细胞类型进行GSVA分析
message("按组织和细胞类型进行GSVA分析...")
# Modified plot_gsva_heatmap function - Fix empty PDF problem
plot_gsva_heatmap <- function(gsva_result, top_n = 20, output_prefix, title) {
  if(is.null(gsva_result) || nrow(gsva_result) == 0) {
    message(paste0("Cannot draw ", output_prefix, " heatmap: Results are empty"))
    return(NULL)
  }
  
  # Ensure there are enough rows
  if(nrow(gsva_result) < 2) {
    message(paste0("Cannot draw ", output_prefix, " heatmap: Too few rows (", nrow(gsva_result), ")"))
    return(NULL)
  }
  
  # Calculate row variance
  row_vars <- apply(gsva_result, 1, var)
  
  # Check if there are rows with variance greater than zero
  if(sum(row_vars > 0) == 0) {
    message(paste0("Cannot draw ", output_prefix, " heatmap: All rows have zero variance"))
    return(NULL)
  }
  
  # Select pathways with the highest variance
  select_n <- min(top_n, sum(row_vars > 0))
  var_pathways <- names(tail(sort(row_vars), select_n))
  
  if(length(var_pathways) == 0) {
    message(paste0("Cannot draw ", output_prefix, " heatmap: Unable to select pathways with highest variance"))
    return(NULL)
  }
  
  message(paste0("Selected ", length(var_pathways), " pathways for ", output_prefix, " heatmap"))
  
  # Safely draw heatmap
  pdf_path <- paste0(output_prefix, "_gsva_heatmap.pdf")
  pdf(pdf_path, width = 10, height = max(4, ceiling(length(var_pathways)/3)))
  
  tryCatch({
    # Set a simple color palette
    col_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
    
    # Draw heatmap and ensure printing
    hmap <- pheatmap(gsva_result[var_pathways, ],
                     scale = "row",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     fontsize_row = 8,
                     fontsize_col = 10,
                     angle_col = "45", # Corrected to string
                     main = title,
                     color = col_palette)
    # Force print heatmap
    print(hmap)
    message(paste0("Successfully drew heatmap and saved to ", pdf_path))
    invisible(hmap)
  }, error = function(e) {
    message(paste0("Error drawing heatmap: ", e$message))
    # Ensure a valid PDF even when error occurs
    plot(1, type="n", axes=FALSE, ann=FALSE)
    text(1, 1, paste0("Error drawing heatmap: ", e$message), cex=1.2)
    NULL
  }, finally = {
    # Ensure PDF device is closed
    dev.off()
  })
}

# Modified plot_gsva_bubble function
plot_gsva_bubble <- function(gsva_result, var_pathways, output_prefix, title) {
  if(is.null(gsva_result) || nrow(gsva_result) == 0) {
    message(paste0("Cannot draw ", output_prefix, " bubble plot: Results are empty"))
    return(NULL)
  }
  
  # Ensure var_pathways is not empty and exists in gsva_result
  if(length(var_pathways) == 0) {
    message(paste0("Cannot draw ", output_prefix, " bubble plot: Pathway list is empty"))
    return(NULL)
  }
  
  # Only keep pathways that exist in gsva_result
  valid_pathways <- intersect(var_pathways, rownames(gsva_result))
  if(length(valid_pathways) == 0) {
    message(paste0("Cannot draw ", output_prefix, " bubble plot: No valid pathways"))
    return(NULL)
  }
  
  message(paste0("Using ", length(valid_pathways), " pathways for ", output_prefix, " bubble plot"))
  
  # Safely draw bubble plot
  pdf_path <- paste0(output_prefix, "_gsva_bubble.pdf")
  pdf(pdf_path, width = 12, height = 10)
  
  tryCatch({
    # Prepare long format data
    gsva_data_long <- reshape2::melt(gsva_result[valid_pathways, ], 
                                     varnames = c("Pathway", "CellType"), 
                                     value.name = "GSVA_Score")
    
    # Draw bubble plot
    p <- ggplot(gsva_data_long, aes(x = CellType, y = Pathway, size = abs(GSVA_Score), color = GSVA_Score)) +
      geom_point(alpha = 0.8) +
      scale_size_continuous(range = c(1, 8), name = "Absolute Score") +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "GSVA Score") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      labs(title = title,
           x = "Cell Type", 
           y = "Pathway")
    # Force print plot
    print(p)
    message(paste0("Successfully drew bubble plot and saved to ", pdf_path))
    invisible(p)
  }, error = function(e) {
    message(paste0("Error drawing bubble plot: ", e$message))
    # Ensure a valid PDF even when error occurs
    plot(1, type="n", axes=FALSE, ann=FALSE)
    text(1, 1, paste0("Error drawing bubble plot: ", e$message), cex=1.2)
    NULL
  }, finally = {
    # Ensure PDF device is closed
    dev.off()
  })
}


# Modified tissue comparison code - Solve empty top_diff_pathways problem
# Tissue comparison section
if(length(tissue_types) >= 2 && length(tissue_gsva_results) >= 2) {
  message("Performing pathway activity comparison between tissues...")
  
  # Get all possible tissue pairs
  tissue_pairs <- combn(names(tissue_gsva_results), 2, simplify = FALSE)
  
  for(pair in tissue_pairs) {
    tissue_1 <- pair[1]
    tissue_2 <- pair[2]
    
    message(paste0("Comparing: ", tissue_1, " vs ", tissue_2))
    
    gsva_tissue_1 <- tissue_gsva_results[[tissue_1]]
    gsva_tissue_2 <- tissue_gsva_results[[tissue_2]]
    
    # Find common cell types between tissues
    common_cell_types <- intersect(colnames(gsva_tissue_1), colnames(gsva_tissue_2))
    
    if(length(common_cell_types) == 0) {
      message(paste0("Tissues '", tissue_1, "' and '", tissue_2, "' have no common cell types"))
      next
    }
    
    message(paste0("Found ", length(common_cell_types), " common cell types"))
    
    # Find common pathways between tissues
    common_pathways <- intersect(rownames(gsva_tissue_1), rownames(gsva_tissue_2))
    
    if(length(common_pathways) == 0) {
      message(paste0("Tissues '", tissue_1, "' and '", tissue_2, "' have no common pathways"))
      next
    }
    
    # Calculate pathway activity differences for each cell type
    for(cell_type in common_cell_types) {
      tryCatch({
        message(paste0("Processing cell type: ", cell_type))
        
        # Extract GSVA scores for this cell type in both tissues
        cell_gsva_1 <- gsva_tissue_1[common_pathways, cell_type, drop = FALSE]
        cell_gsva_2 <- gsva_tissue_2[common_pathways, cell_type, drop = FALSE]
        
        # Calculate differences
        diff_gsva <- cell_gsva_1 - cell_gsva_2
        colnames(diff_gsva) <- paste0(cell_type, "_diff")
        
        # Save difference results
        output_csv <- paste0("gsva_diff_", gsub(" ", "_", cell_type), "_", 
                             gsub(" ", "_", tissue_1), "_vs_", 
                             gsub(" ", "_", tissue_2), ".csv")
        write.csv(diff_gsva, file = output_csv)
        message(paste0("Difference results saved to: ", output_csv))
        
        # Get pathways with largest differences, handle safely
        if(nrow(diff_gsva) == 0) {
          message(paste0("Cell type '", cell_type, "' has no difference data, skipping visualization"))
          next
        }
        
        # Calculate and sort absolute differences
        abs_diffs <- abs(diff_gsva[,1])
        if(length(abs_diffs) == 0 || all(is.na(abs_diffs)) || all(abs_diffs == 0)) {
          message(paste0("Cell type '", cell_type, "' has no valid differences, skipping visualization"))
          next
        }
        
        # Sort by absolute difference value
        sorted_idx <- order(abs_diffs, decreasing = TRUE)
        # Select top 15 or all pathways (if fewer than 15)
        n_paths <- min(15, length(sorted_idx))
        top_idx <- sorted_idx[1:n_paths]
        
        # Ensure indices are valid
        if(length(top_idx) == 0) {
          message(paste0("Cell type '", cell_type, "' cannot select pathways with largest differences"))
          next
        }
        
        # Get pathway names
        top_diff_pathways <- rownames(diff_gsva)[top_idx]
        message(paste0("Selected ", length(top_diff_pathways), " pathways with largest differences"))
        
        # Draw bar plot
        pdf_path <- paste0("gsva_diff_barplot_", gsub(" ", "_", cell_type), "_", 
                           gsub(" ", "_", tissue_1), "_vs_", 
                           gsub(" ", "_", tissue_2), ".pdf")
        pdf(pdf_path, width = 10, height = 8)
        
        tryCatch({
          # Sort pathways by difference value
          ordered_paths <- top_diff_pathways[order(diff_gsva[top_diff_pathways, 1])]
          
          # Create data frame
          barplot_data <- data.frame(
            Pathway = factor(ordered_paths, levels = ordered_paths),
            Difference = diff_gsva[ordered_paths, 1]
          )
          
          # Draw bar plot
          p <- ggplot(barplot_data, aes(x = Difference, y = Pathway, fill = Difference > 0)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#00468B"),
                              labels = c("TRUE" = tissue_1, "FALSE" = tissue_2),
                              name = "Higher activity in") +
            labs(title = paste0("Pathway Activity Differences in ", cell_type),
                 subtitle = paste0(tissue_1, " vs ", tissue_2),
                 x = "GSVA Score Difference") +
            theme_bw() +
            theme(axis.text.y = element_text(size = 9),
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 12))
          
          # Ensure plot is printed
          print(p)
          message(paste0("Successfully drew bar plot and saved to ", pdf_path))
        }, error = function(e) {
          message(paste0("Error drawing bar plot: ", e$message))
          # Ensure a valid PDF even when error occurs
          plot(1, type="n", axes=FALSE, ann=FALSE)
          text(1, 1, paste0("Error drawing bar plot: ", e$message), cex=1.2)
        }, finally = {
          # Ensure PDF device is closed
          dev.off()
        })
        
      }, error = function(e) {
        message(paste0("Error processing cell type '", cell_type, "': ", e$message))
      })
    }
    
    # Create tissue difference heatmap - also add safety checks
    message("Creating tissue difference heatmap...")
    if(length(common_cell_types) >= 2 && length(common_pathways) > 0) {
      # Create difference matrix
      diff_matrix <- matrix(NA, nrow = length(common_pathways), ncol = length(common_cell_types))
      rownames(diff_matrix) <- common_pathways
      colnames(diff_matrix) <- common_cell_types
      
      # Fill difference matrix
      for(i in 1:length(common_cell_types)) {
        cell_type <- common_cell_types[i]
        if(cell_type %in% colnames(gsva_tissue_1) && cell_type %in% colnames(gsva_tissue_2)) {
          diff_matrix[, i] <- gsva_tissue_1[common_pathways, cell_type] - 
            gsva_tissue_2[common_pathways, cell_type]
        }
      }
      
      # Remove rows and columns with all NAs
      diff_matrix <- diff_matrix[rowSums(!is.na(diff_matrix)) > 0, colSums(!is.na(diff_matrix)) > 0, drop = FALSE]
      
      if(nrow(diff_matrix) > 0 && ncol(diff_matrix) > 0) {
        # Calculate pathways with highest variance, handle safely
        if(nrow(diff_matrix) < 2) {
          message("Difference matrix has too few rows, cannot draw heatmap")
        } else {
          row_vars <- apply(diff_matrix, 1, var, na.rm = TRUE)
          valid_rows <- which(!is.na(row_vars) & row_vars > 0)
          
          if(length(valid_rows) == 0) {
            message("No pathways with variance greater than zero, cannot draw heatmap")
          } else {
            # Select pathways with highest variance
            n_paths <- min(25, length(valid_rows))
            top_var_idx <- tail(order(row_vars[valid_rows]), n_paths)
            var_paths <- rownames(diff_matrix)[valid_rows[top_var_idx]]
            
            # 绘制热图
            pdf_path <- paste0("tissue_diff_heatmap_", gsub(" ", "_", tissue_1), "_vs_", 
                               gsub(" ", "_", tissue_2), ".pdf")
            pdf(pdf_path, width = 10, height = 12)
            
            tryCatch({
              # 数据预处理，处理可能的NA值
              plot_data <- diff_matrix[var_paths, ]
              
              # 使用pheatmap尝试绘制
              hmap <- pheatmap(plot_data,
                               scale = "row",
                               cluster_rows = TRUE,
                               cluster_cols = TRUE,
                               fontsize_row = 8,
                               fontsize_col = 10,
                               angle_col = "45",
                               main = paste0("Pathway Activity Differences: ", tissue_1, " vs ", tissue_2),
                               color = colorRampPalette(c("#00468B", "white", "#E64B35"))(50))
              
              # 确保热图被打印
              print(hmap)
              message(paste0("Successfully drew tissue difference heatmap and saved to ", pdf_path))
            }, error = function(e) {
              message(paste0("pheatmap error: ", e$message, ". Trying base heatmap."))
              
              # 如果pheatmap失败，尝试使用基础热图
              tryCatch({
                # 处理数据
                plot_data <- as.matrix(diff_matrix[var_paths, ])
                plot_data[is.na(plot_data)] <- 0
                
                # 使用基础热图
                heatmap(plot_data,
                        main = paste0("Pathway Activity Differences: ", tissue_1, " vs ", tissue_2),
                        col = colorRampPalette(c("#00468B", "white", "#E64B35"))(50),
                        margins = c(8, 8))
                
                message("Successfully drew tissue difference heatmap using base heatmap function")
              }, error = function(e2) {
                message(paste0("Base heatmap also failed: ", e2$message))
                # 如果所有方法都失败，画一个简单的图
                plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
                text(1.5, 1.5, "Heatmap creation failed. See R console for details.", cex=1.2)
              })
            }, finally = {
              # 确保PDF设备关闭
              dev.off()
            })
            
          }
        }
      } else {
        message("Processed difference matrix is empty, cannot draw heatmap")
      }
    } else {
      message("Insufficient number of common cell types or pathways, cannot create tissue difference heatmap")
    }
  }
}

# 恢复原工作目录
setwd(current_dir)
message("GSVA分析完成！结果保存在:", normalizePath(output_dir))

################################################################################################################
################################################################################################################
################################################################################################################
tissue_gsva_kegg_results <- list()
# 在组织分析循环中添加KEGG基因集分析
for(tissue in tissue_types) {
  message(paste0("Processing tissue with KEGG gene sets: ", tissue))
  
  # 获取该组织的表达矩阵（假设已经计算）
  cells_in_tissue <- rownames(T_object@meta.data[T_object@meta.data$tissue == tissue, ])
  
  if(length(cells_in_tissue) < 50) {
    message(paste0("Tissue '", tissue, "' has too few cells (", length(cells_in_tissue), "), skipping"))
    next
  }
  
  # 创建该组织的子集
  tissue_subset <- subset(T_object, cells = cells_in_tissue)
  
  # 设置细胞类型为标识
  Idents(tissue_subset) <- tissue_subset$Annotation_2
  
  # 计算平均表达
  expr_tissue <- AverageExpression(tissue_subset, assays = "RNA", slot = "data")[[1]]
  expr_tissue <- expr_tissue[rowSums(expr_tissue) > 0,]
  expr_tissue <- as.matrix(expr_tissue)
  
  # 运行KEGG GSVA
  message(paste0("Running KEGG GSVA for tissue ", tissue))
  gsva_tissue_kegg <- run_gsva_analysis(
    expr_data = expr_tissue, 
    gene_sets = c2_kegg_genesets, 
    output_prefix = paste0(gsub(" ", "_", tissue), "_kegg"),
    method = "gsva",
    mx.diff = TRUE
  )
  
  # 检查GSVA结果
  if(!is.null(gsva_tissue_kegg)) {
    message(paste0("KEGG GSVA result dimensions: ", nrow(gsva_tissue_kegg), " × ", ncol(gsva_tissue_kegg)))
    
    # 保存结果
    tissue_gsva_kegg_results[[tissue]] <- gsva_tissue_kegg
    
    # 绘制热图
    message("Drawing KEGG heatmap...")
    plot_gsva_heatmap(
      gsva_result = gsva_tissue_kegg, 
      top_n = 30, 
      output_prefix = paste0("tissue_", gsub(" ", "_", tissue), "_kegg"), 
      title = paste0("KEGG Pathway Activity in ", tissue)
    )
  } else {
    message(paste0("KEGG GSVA result is NULL for tissue ", tissue))
  }
}
# KEGG版本的组织间比较分析
if(length(tissue_types) >= 2 && length(tissue_gsva_kegg_results) >= 2) {
  message("Performing KEGG pathway activity comparison between tissues...")
  
  # 获取所有可能的组织对比
  tissue_pairs <- combn(names(tissue_gsva_kegg_results), 2, simplify = FALSE)
  message(paste0("Total KEGG comparison pairs: ", length(tissue_pairs)))
  
  for(pair in tissue_pairs) {
    tissue_1 <- pair[1]
    tissue_2 <- pair[2]
    
    message(paste0("Comparing KEGG pathways: ", tissue_1, " vs ", tissue_2))
    
    gsva_tissue_1 <- tissue_gsva_kegg_results[[tissue_1]]
    gsva_tissue_2 <- tissue_gsva_kegg_results[[tissue_2]]
    
    # 找到两个组织共有的细胞类型
    common_cell_types <- intersect(colnames(gsva_tissue_1), colnames(gsva_tissue_2))
    
    if(length(common_cell_types) == 0) {
      message(paste0("Tissues '", tissue_1, "' and '", tissue_2, "' have no common cell types for KEGG analysis"))
      next
    }
    
    message(paste0("Found ", length(common_cell_types), " common cell types for KEGG analysis: ", 
                   paste(common_cell_types, collapse=", ")))
    
    # 找到两个组织中共有的KEGG通路
    common_pathways <- intersect(rownames(gsva_tissue_1), rownames(gsva_tissue_2))
    message(paste0("Found ", length(common_pathways), " common KEGG pathways"))
    
    if(length(common_pathways) == 0) {
      message(paste0("Tissues '", tissue_1, "' and '", tissue_2, "' have no common KEGG pathways"))
      next
    }
    
    # 创建KEGG组织差异结果目录
    diff_dir <- paste0("kegg_diff_", gsub(" ", "_", tissue_1), "_vs_", gsub(" ", "_", tissue_2))
    dir.create(diff_dir, showWarnings = FALSE, recursive = TRUE)
    
    # 计算每个细胞类型中两个组织间的KEGG通路差异
    cell_processed <- 0
    cell_success <- 0
    cell_heatmaps <- list()
    
    for(cell_type in common_cell_types) {
      tryCatch({
        cell_processed <- cell_processed + 1
        message(paste0("Processing cell type for KEGG comparison (", cell_processed, "/", length(common_cell_types), "): ", cell_type))
        
        # 提取该细胞类型在两个组织中的GSVA分数
        cell_gsva_1 <- gsva_tissue_1[common_pathways, cell_type, drop = FALSE]
        cell_gsva_2 <- gsva_tissue_2[common_pathways, cell_type, drop = FALSE]
        
        # 计算差异
        diff_gsva <- cell_gsva_1 - cell_gsva_2
        colnames(diff_gsva) <- paste0(cell_type, "_diff")
        
        # 保存差异结果
        output_csv <- file.path(diff_dir, paste0("gsva_kegg_diff_", gsub(" ", "_", cell_type), ".csv"))
        write.csv(diff_gsva, file = output_csv)
        message(paste0("KEGG difference results saved to: ", output_csv))
        
        # 安全处理差异最大的通路选择
        if(nrow(diff_gsva) == 0) {
          message(paste0("Cell type '", cell_type, "' has no KEGG difference data, skipping visualization"))
          next
        }
        
        # 计算绝对差异
        abs_diffs <- abs(diff_gsva[,1])
        
        # 处理无效数据
        if(length(abs_diffs) == 0 || all(is.na(abs_diffs)) || all(abs_diffs == 0)) {
          message(paste0("Cell type '", cell_type, "' has no valid KEGG differences, skipping visualization"))
          next
        }
        
        # 排序并选择差异最大的通路
        sorted_idx <- order(abs_diffs, decreasing = TRUE)
        n_paths <- min(15, length(sorted_idx))
        top_diff_pathways <- rownames(diff_gsva)[sorted_idx[1:n_paths]]
        
        message(paste0("Selected ", length(top_diff_pathways), " KEGG pathways with largest differences"))
        
        # 存储选中的通路以便后续组合热图
        cell_heatmaps[[cell_type]] <- list(
          diff_data = diff_gsva[top_diff_pathways, , drop = FALSE],
          top_paths = top_diff_pathways
        )
        
        # 绘制条形图
        pdf_path <- file.path(diff_dir, paste0("gsva_kegg_diff_barplot_", gsub(" ", "_", cell_type), ".pdf"))
        
        pdf(pdf_path, width = 12, height = 10) # 宽度增加，以适应可能较长的KEGG通路名称
        
        tryCatch({
          # 按差异值排序通路
          ordered_paths <- top_diff_pathways[order(diff_gsva[top_diff_pathways, 1])]
          
          # 创建数据框
          barplot_data <- data.frame(
            Pathway = factor(ordered_paths, levels = ordered_paths),
            Difference = diff_gsva[ordered_paths, 1]
          )
          
          # 绘制条形图 - 调整以适应KEGG通路名称
          p <- ggplot(barplot_data, aes(x = Difference, y = Pathway, fill = Difference > 0)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#00468B"),
                              labels = c("TRUE" = tissue_1, "FALSE" = tissue_2),
                              name = "Higher activity in") +
            labs(title = paste0("KEGG Pathway Activity Differences in ", cell_type),
                 subtitle = paste0(tissue_1, " vs ", tissue_2),
                 x = "GSVA Score Difference") +
            theme_bw() +
            theme(axis.text.y = element_text(size = 8), # 字体略小，以适应更多文本
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 12))
          
          # 确保图形被打印
          print(p)
          message(paste0("Successfully drew KEGG bar plot and saved to ", pdf_path))
          cell_success <- cell_success + 1
          
        }, error = function(e) {
          message(paste0("Error drawing KEGG bar plot: ", e$message))
          # 确保在出错时仍然能生成有效的PDF
          plot(1, type="n", axes=FALSE, ann=FALSE)
          text(1, 1, paste0("Error drawing KEGG bar plot: ", e$message), cex=1.2)
        }, finally = {
          # 确保关闭PDF设备
          dev.off()
        })
        
      }, error = function(e) {
        message(paste0("Error processing cell type '", cell_type, "' for KEGG: ", e$message))
      })
    }
    
    message(paste0("Successfully processed ", cell_success, " out of ", length(common_cell_types), " cell types for KEGG comparison"))
    
    # 创建组织间KEGG差异热图
    message("Creating tissue KEGG difference heatmap...")
    if(length(common_cell_types) >= 2 && length(common_pathways) > 0) {
      # 创建差异矩阵
      diff_matrix <- matrix(NA, nrow = length(common_pathways), ncol = length(common_cell_types))
      rownames(diff_matrix) <- common_pathways
      colnames(diff_matrix) <- common_cell_types
      
      # 填充差异矩阵
      for(i in 1:length(common_cell_types)) {
        cell_type <- common_cell_types[i]
        if(cell_type %in% colnames(gsva_tissue_1) && cell_type %in% colnames(gsva_tissue_2)) {
          diff_matrix[, i] <- gsva_tissue_1[common_pathways, cell_type] - 
            gsva_tissue_2[common_pathways, cell_type]
        }
      }
      
      # 移除全NA的行和列
      diff_matrix <- diff_matrix[rowSums(!is.na(diff_matrix)) > 0, colSums(!is.na(diff_matrix)) > 0, drop = FALSE]
      
      message(paste0("KEGG difference matrix dimensions: ", nrow(diff_matrix), " × ", ncol(diff_matrix)))
      
      if(nrow(diff_matrix) > 0 && ncol(diff_matrix) > 0) {
        # 由于KEGG通路数量可能很多，限制选择的通路数量
        if(nrow(diff_matrix) < 2) {
          message("KEGG difference matrix has too few rows, cannot draw heatmap")
        } else {
          # 计算方差
          row_vars <- apply(diff_matrix, 1, var, na.rm = TRUE)
          valid_rows <- which(!is.na(row_vars) & row_vars > 0)
          
          if(length(valid_rows) == 0) {
            message("No KEGG pathways with variance greater than zero, cannot draw heatmap")
          } else {
            # 选择变异最大的通路，但限制数量以避免热图过大
            n_paths <- min(30, length(valid_rows))
            top_var_idx <- tail(order(row_vars[valid_rows]), n_paths)
            var_paths <- rownames(diff_matrix)[valid_rows[top_var_idx]]
            
            message(paste0("Selected ", length(var_paths), " KEGG pathways with highest variance for heatmap"))
            
            # 保存差异矩阵数据
            diff_matrix_csv <- file.path(diff_dir, "tissue_kegg_diff_matrix.csv")
            write.csv(diff_matrix, file = diff_matrix_csv)
            message(paste0("Saved complete KEGG difference matrix to ", diff_matrix_csv))
            
            # 绘制热图
            pdf_path <- file.path(diff_dir, paste0("tissue_kegg_diff_heatmap.pdf"))
            pdf(pdf_path, width = 12, height = 14) # 增加高度以适应更多KEGG通路
            
            # 多种热图绘制方法尝试
            heatmap_success <- FALSE
            
            # 1. 尝试pheatmap
            tryCatch({
              plot_data <- diff_matrix[var_paths, ]
              
              # 数据检查
              message(paste0("KEGG heatmap data dimensions: ", nrow(plot_data), " × ", ncol(plot_data)))
              message(paste0("Any NA values: ", any(is.na(plot_data))))
              message(paste0("Data range: [", min(plot_data, na.rm = TRUE), ", ", max(plot_data, na.rm = TRUE), "]"))
              
              # 调整字体大小以适应KEGG通路名称
              hmap <- pheatmap(plot_data,
                               scale = "row",
                               cluster_rows = TRUE,
                               cluster_cols = TRUE,
                               fontsize_row = 6, # 字体更小以适应更多通路
                               fontsize_col = 10,
                               angle_col = "45",
                               main = paste0("KEGG Pathway Activity Differences: ", tissue_1, " vs ", tissue_2),
                               color = colorRampPalette(c("#00468B", "white", "#E64B35"))(50))
              
              # 确保热图被打印
              print(hmap)
              message("Successfully drew KEGG tissue difference heatmap with pheatmap")
              heatmap_success <- TRUE
              
            }, error = function(e) {
              message(paste0("pheatmap error for KEGG: ", e$message, ". Trying base heatmap."))
              
              # 2. 如果pheatmap失败，尝试基础热图
              tryCatch({
                # 数据处理
                plot_data <- as.matrix(diff_matrix[var_paths, ])
                plot_data[is.na(plot_data)] <- 0
                
                # 基础热图
                heatmap(plot_data,
                        main = paste0("KEGG Pathway Activity Differences: ", tissue_1, " vs ", tissue_2),
                        col = colorRampPalette(c("#00468B", "white", "#E64B35"))(50),
                        margins = c(8, 10),  # 边距更大以适应更长的路径名
                        cexRow = 0.6)        # 行标签字体更小
                
                message("Successfully drew KEGG tissue difference heatmap using base heatmap function")
                heatmap_success <- TRUE
                
              }, error = function(e2) {
                message(paste0("Base heatmap also failed for KEGG: ", e2$message, ". Trying image plot."))
                
                # 3. 如果基础热图也失败，使用最简单的image函数
                tryCatch({
                  plot_data <- as.matrix(diff_matrix[var_paths, ])
                  plot_data[is.na(plot_data)] <- 0
                  
                  # 标准化数据用于可视化
                  z <- t(scale(t(plot_data)))
                  
                  # 使用image函数
                  image(1:ncol(z), 1:nrow(z), z, 
                        main = paste0("KEGG Pathway Activity Differences: ", tissue_1, " vs ", tissue_2),
                        xlab = "Cell Types", ylab = "KEGG Pathways",
                        col = colorRampPalette(c("blue", "white", "red"))(50),
                        axes = FALSE)
                  
                  # 添加轴标签
                  axis(1, at = 1:ncol(z), labels = colnames(z), las = 2, cex.axis = 0.7)
                  axis(2, at = 1:nrow(z), labels = rownames(z), las = 2, cex.axis = 0.5)  # 字体更小
                  
                  message("Successfully drew KEGG tissue difference heatmap using image function")
                  heatmap_success <- TRUE
                  
                }, error = function(e3) {
                  message(paste0("All KEGG heatmap methods failed: ", e3$message))
                  
                  # 4. 如果所有方法都失败，绘制一个简单的错误消息
                  plot(c(1,2), c(1,2), type="n", axes=FALSE, ann=FALSE)
                  text(1, 1, "KEGG heatmap creation failed. See R console for details.", cex=1.2)
                })
              })
            }, finally = {
              # 确保PDF设备关闭
              dev.off()
              message(paste0("KEGG tissue difference heatmap ", 
                             ifelse(heatmap_success, "successfully saved", "attempt failed"), 
                             " to ", pdf_path))
            })
          }
        }
      } else {
        message("Processed KEGG difference matrix is empty, cannot draw heatmap")
      }
    } else {
      message("Insufficient number of common cell types or KEGG pathways, cannot create tissue difference heatmap")
    }
    
    # 创建组合热图展示所有细胞类型的Top差异KEGG通路
    if(length(cell_heatmaps) > 0) {
      message("Creating combined top differential KEGG pathways heatmap...")
      
      # 收集所有细胞类型的Top差异通路
      all_top_paths <- unique(unlist(lapply(cell_heatmaps, function(x) x$top_paths)))
      
      if(length(all_top_paths) > 0 && length(all_top_paths) <= 100) {  # 限制通路数量
        # 创建合并的差异矩阵
        combined_diff <- diff_matrix[all_top_paths, , drop = FALSE]
        
        if(nrow(combined_diff) > 0 && ncol(combined_diff) > 0) {
          pdf_path <- file.path(diff_dir, paste0("combined_top_kegg_diff_heatmap.pdf"))
          pdf(pdf_path, width = 14, height = min(30, length(all_top_paths)/1.5))  # 更大的PDF以适应KEGG通路
          
          tryCatch({
            # 尝试使用pheatmap绘制合并热图
            hmap <- pheatmap(combined_diff,
                             scale = "row",
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             fontsize_row = 5,  # 很小的字体以适应更多通路
                             fontsize_col = 10,
                             angle_col = "45",
                             main = paste0("Top Differential KEGG Pathways: ", tissue_1, " vs ", tissue_2),
                             color = colorRampPalette(c("#00468B", "white", "#E64B35"))(50))
            
            print(hmap)
            message(paste0("Successfully created combined top differential KEGG pathways heatmap"))
            
          }, error = function(e) {
            message(paste0("Error creating combined KEGG heatmap: ", e$message))
            # 尝试基础热图
            heatmap(as.matrix(combined_diff),
                    main = paste0("Top Differential KEGG Pathways: ", tissue_1, " vs ", tissue_2),
                    col = colorRampPalette(c("blue", "white", "red"))(50),
                    margins = c(8, 10),
                    cexRow = 0.5)  # 更小的行标签字体
          }, finally = {
            dev.off()
          })
        }
      } else if(length(all_top_paths) > 100) {
        message("Too many top KEGG pathways (", length(all_top_paths), "), limiting to top 100 by variance")
        
        # 计算所有通路的方差并选择top 100
        sub_diff_matrix <- diff_matrix[all_top_paths, , drop = FALSE]
        path_vars <- apply(sub_diff_matrix, 1, var, na.rm = TRUE)
        top_var_paths <- names(sort(path_vars, decreasing = TRUE)[1:100])
        
        # 创建限制后的合并差异矩阵
        combined_diff <- diff_matrix[top_var_paths, , drop = FALSE]
        
        # 绘制热图
        pdf_path <- file.path(diff_dir, paste0("combined_top_kegg_diff_heatmap.pdf"))
        pdf(pdf_path, width = 14, height = 30)  # 固定高度
        
        tryCatch({
          # 尝试使用pheatmap绘制合并热图
          hmap <- pheatmap(combined_diff,
                           scale = "row",
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           fontsize_row = 5,
                           fontsize_col = 10,
                           angle_col = "45",
                           main = paste0("Top 100 Differential KEGG Pathways by Variance: ", tissue_1, " vs ", tissue_2),
                           color = colorRampPalette(c("#00468B", "white", "#E64B35"))(50))
          
          print(hmap)
          message(paste0("Successfully created limited combined top differential KEGG pathways heatmap"))
          
        }, error = function(e) {
          message(paste0("Error creating limited combined KEGG heatmap: ", e$message))
          # 尝试基础热图
          heatmap(as.matrix(combined_diff),
                  main = paste0("Top 100 Differential KEGG Pathways by Variance: ", tissue_1, " vs ", tissue_2),
                  col = colorRampPalette(c("blue", "white", "red"))(50),
                  margins = c(8, 10),
                  cexRow = 0.4)
        }, finally = {
          dev.off()
        })
      }
    }
  }
  
  message("KEGG tissue comparison analysis completed!")
} else {
  message("Insufficient number of tissues with KEGG results for comparison (need at least 2)")
}


####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

# =============================================================================
# 整合版STARTRAC-dist组织分布分析系统 (Integrated STARTRAC-dist Analysis System)
# 融合稳健性改进、完整统计方法和高质量可视化
# 作者：Claude & 用户协作开发
# =============================================================================

# 加载必需的包
required_packages <- c(
  "Startrac", "Seurat", "ggplot2", "ComplexHeatmap", "RColorBrewer", 
  "circlize", "tidyverse", "ggpubr", "pheatmap", "corrplot", 
  "cowplot", "viridis", "scales", "readr", "qs", 
  "gridExtra", "reshape2", "VennDiagram", "grid", "yaml"
)

# 智能包管理函数
setup_packages <- function(packages) {
  missing_packages <- packages[!sapply(packages, require, character.only = TRUE, quietly = TRUE)]
  
  if(length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    
    # 特殊处理某些包
    for(pkg in missing_packages) {
      tryCatch({
        if(pkg == "Startrac") {
          # STARTRAC可能需要从GitHub安装
          if(!require("devtools", quietly = TRUE)) install.packages("devtools")
          devtools::install_github("Japrin/STARTRAC")
        } else {
          install.packages(pkg)
        }
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        warning(paste("Failed to install", pkg, ":", e$message))
      })
    }
  }
  
  # 设置全局主题
  if(require("ggplot2", quietly = TRUE)) {
    theme_set(theme_bw() + theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ))
  }
}

setup_packages(required_packages)

# =============================================================================
# 核心数据处理与验证函数 (Core Data Processing & Validation)
# =============================================================================

#' 稳健的Ro/e值转换为符号表示 (按STARTRAC标准)
#' @param roe_value Ro/e指数值
#' @param thresholds 阈值向量，默认按STARTRAC标准c(0.5, 0.8, 1.2, 2.0)
#' @return 符号字符串
roe_to_symbol <- function(roe_value, thresholds = c(0.5, 0.8, 1.2, 2.0)) {
  # 处理各种异常值
  if(is.null(roe_value) || length(roe_value) == 0 || is.na(roe_value) || is.nan(roe_value)) {
    return("-")
  }
  
  if(is.infinite(roe_value)) {
    return(ifelse(roe_value > 0, "+++", "---"))
  }
  
  if(!is.numeric(roe_value)) return("-")
  
  # 按照STARTRAC标准分类
  if(roe_value >= thresholds[4]) {        # >= 2.0
    return("+++")
  } else if(roe_value >= thresholds[3]) { # >= 1.2
    return("+")
  } else if(roe_value > thresholds[2]) {  # > 0.8
    return("+/-")
  } else if(roe_value >= thresholds[1]) { # >= 0.5
    return("-")
  } else {                                # < 0.5
    return("---")
  }
}

#' Ro/e值的生物学解释
#' @param roe_value Ro/e指数值
get_roe_interpretation <- function(roe_value) {
  if(is.na(roe_value) || !is.finite(roe_value)) {
    return("Data not available")
  } else if(roe_value >= 2.0) {
    return("Strongly enriched")
  } else if(roe_value >= 1.5) {
    return("Moderately enriched")
  } else if(roe_value >= 1.2) {
    return("Mildly enriched")
  } else if(roe_value >= 0.8) {
    return("Near expected levels")
  } else if(roe_value >= 0.5) {
    return("Mildly depleted")
  } else if(roe_value >= 0.3) {
    return("Moderately depleted")
  } else {
    return("Strongly depleted")
  }
}

#' 数据预处理和验证
#' @param roe_matrix Ro/e指数矩阵
#' @return 清洁的Ro/e矩阵
clean_roe_matrix <- function(roe_matrix) {
  message("Cleaning and validating Ro/e matrix...")
  
  # 检查输入
  if(is.null(roe_matrix) || (!is.matrix(roe_matrix) && !is.data.frame(roe_matrix))) {
    stop("Invalid Ro/e matrix provided")
  }
  
  # 转换为矩阵
  if(is.data.frame(roe_matrix)) {
    roe_matrix <- as.matrix(roe_matrix)
  }
  
  if(nrow(roe_matrix) == 0 || ncol(roe_matrix) == 0) {
    stop("Empty Ro/e matrix")
  }
  
  # 替换无效值
  roe_matrix[is.nan(roe_matrix)] <- NA
  roe_matrix[is.infinite(roe_matrix)] <- NA
  
  # 检查数据范围合理性
  finite_values <- roe_matrix[is.finite(roe_matrix)]
  
  message(sprintf("  Matrix dimensions: %d cell types × %d tissues", 
                  nrow(roe_matrix), ncol(roe_matrix)))
  message(sprintf("  Valid values: %d/%d (%.1f%%)", 
                  sum(is.finite(roe_matrix)), length(roe_matrix),
                  100 * sum(is.finite(roe_matrix)) / length(roe_matrix)))
  
  if(length(finite_values) > 0) {
    message(sprintf("  Value range: %.3f - %.3f", 
                    min(finite_values, na.rm = TRUE), 
                    max(finite_values, na.rm = TRUE)))
  } else {
    warning("No finite values in Ro/e matrix")
  }
  
  return(roe_matrix)
}

#' 显示可用的metadata列信息
#' @param seurat_obj Seurat对象
#' @param show_preview 是否显示每列的预览
show_metadata_info <- function(seurat_obj, show_preview = TRUE) {
  
  meta_data <- seurat_obj@meta.data
  available_cols <- colnames(meta_data)
  
  message("\n=== AVAILABLE METADATA COLUMNS ===")
  message("Total columns: ", length(available_cols))
  
  for(i in seq_along(available_cols)) {
    col_name <- available_cols[i]
    col_data <- meta_data[[col_name]]
    
    # 基本信息
    unique_count <- length(unique(col_data[!is.na(col_data)]))
    na_count <- sum(is.na(col_data))
    completion_rate <- round((1 - na_count/length(col_data)) * 100, 1)
    
    message(sprintf("  %2d: %-30s [%d unique, %.1f%% complete]", 
                    i, col_name, unique_count, completion_rate))
    
    # 显示预览
    if(show_preview) {
      if(is.numeric(col_data)) {
        range_info <- paste0("Range: ", round(min(col_data, na.rm = TRUE), 2), 
                             " - ", round(max(col_data, na.rm = TRUE), 2))
        message(sprintf("      %s", range_info))
      } else {
        unique_vals <- unique(col_data[!is.na(col_data)])
        preview_vals <- if(length(unique_vals) <= 5) {
          paste(unique_vals, collapse = ", ")
        } else {
          paste(c(head(unique_vals, 3), "..."), collapse = ", ")
        }
        message(sprintf("      Values: %s", preview_vals))
      }
    }
  }
  
  return(available_cols)
}

# =============================================================================
# 稳健的STARTRAC-dist可视化系统 (Robust STARTRAC-dist Visualization)
# =============================================================================

#' 创建稳健的STARTRAC-dist风格可视化 (主函数)
#' @param roe_matrix Ro/e指数矩阵
#' @param significance_results 显著性结果（可选）
#' @param output_path 输出文件路径
#' @param title 图像标题
#' @param method 绘图方法 ("grid", "base", "auto")
#' @param cell_type_groups 细胞类型分组（可选）
create_startrac_dist_plot <- function(roe_matrix, 
                                      significance_results = NULL, 
                                      output_path = "startrac_dist_plot.pdf",
                                      title = "STARTRAC-dist Analysis",
                                      method = "auto",
                                      cell_type_groups = NULL) {
  
  message("Creating STARTRAC-dist style visualization...")
  
  # 1. 数据预处理
  roe_matrix <- clean_roe_matrix(roe_matrix)
  
  # 基本参数
  tissues <- colnames(roe_matrix)
  celltypes <- rownames(roe_matrix)
  n_tissues <- length(tissues)
  n_celltypes <- length(celltypes)
  
  if(is.null(tissues) || is.null(celltypes)) {
    stop("Matrix must have row and column names")
  }
  
  # 2. 创建符号矩阵
  symbol_matrix <- matrix("", nrow = n_celltypes, ncol = n_tissues)
  rownames(symbol_matrix) <- celltypes
  colnames(symbol_matrix) <- tissues
  
  for(i in 1:n_celltypes) {
    for(j in 1:n_tissues) {
      symbol_matrix[i, j] <- roe_to_symbol(roe_matrix[i, j])
    }
  }
  
  # 3. 智能选择绘图方法
  if(method == "auto") {
    # 根据数据大小和复杂性选择方法
    if(n_tissues * n_celltypes > 500 || any(nchar(c(tissues, celltypes)) > 20)) {
      method <- "base"  # 大数据或长标签用基础绘图
    } else {
      method <- "grid"  # 默认使用grid系统
    }
  }
  
  # 4. 执行绘图
  if(method == "grid") {
    create_grid_plot(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes)
  } else {
    create_base_plot(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes)
  }
  
  message(paste("STARTRAC-dist plot saved to:", output_path))
}

#' Grid系统绘图方法 (精确布局)
create_grid_plot <- function(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes) {
  
  tryCatch({
    n_tissues <- length(tissues)
    n_celltypes <- length(celltypes)
    
    pdf(output_path, width = max(12, n_tissues * 1.5), height = max(8, n_celltypes * 0.4))
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2, widths = c(1, 1))))
    
    # === 左侧：符号矩阵图 ===
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    pushViewport(viewport(x = 0.15, y = 0.15, width = 0.7, height = 0.7, just = c("left", "bottom")))
    
    cell_width <- 1 / n_tissues
    cell_height <- 1 / n_celltypes
    
    # 颜色映射函数
    get_colors <- function(roe_val) {
      if(!is.finite(roe_val)) return(list(bg = "#F0F0F0", text = "black"))
      
      if(roe_val >= 2.0) return(list(bg = "#8B0000", text = "white"))
      if(roe_val >= 1.2) return(list(bg = "#DC143C", text = "white"))
      if(roe_val > 0.8) return(list(bg = "#FFFACD", text = "black"))
      if(roe_val >= 0.5) return(list(bg = "#87CEEB", text = "black"))
      return(list(bg = "#4682B4", text = "white"))
    }
    
    # 绘制符号网格
    for(i in 1:n_celltypes) {
      for(j in 1:n_tissues) {
        x <- (j - 0.5) / n_tissues
        y <- (n_celltypes - i + 0.5) / n_celltypes
        
        roe_val <- roe_matrix[i, j]
        colors <- get_colors(roe_val)
        symbol <- symbol_matrix[i, j]
        
        grid.rect(x = x, y = y, 
                  width = cell_width * 0.9, height = cell_height * 0.9,
                  gp = gpar(fill = colors$bg, col = "white", lwd = 1))
        grid.text(symbol, x = x, y = y, 
                  gp = gpar(col = colors$text, fontsize = 11, fontface = "bold"))
      }
    }
    
    # 添加标签
    for(i in 1:n_celltypes) {
      y <- (n_celltypes - i + 0.5) / n_celltypes
      grid.text(celltypes[i], x = -0.02, y = y, just = "right", gp = gpar(fontsize = 9))
    }
    
    for(j in 1:n_tissues) {
      x <- (j - 0.5) / n_tissues
      grid.text(tissues[j], x = x, y = -0.02, just = "center", rot = 45, gp = gpar(fontsize = 9))
    }
    
    grid.text("STARTRAC-dist Symbols", x = 0.5, y = 1.05, just = "center",
              gp = gpar(fontsize = 14, fontface = "bold"))
    
    popViewport(2)
    
    # === 右侧：数值热图 ===
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    pushViewport(viewport(x = 0.15, y = 0.15, width = 0.7, height = 0.7, just = c("left", "bottom")))
    
    # 创建颜色渐变
    finite_values <- roe_matrix[is.finite(roe_matrix)]
    if(length(finite_values) > 0) {
      val_min <- min(finite_values, na.rm = TRUE)
      val_max <- max(finite_values, na.rm = TRUE)
      color_func <- colorRamp2(c(val_min, 1, val_max), c("#2166AC", "white", "#B2182B"))
    } else {
      color_func <- colorRamp2(c(0, 1, 2), c("#2166AC", "white", "#B2182B"))
    }
    
    # 绘制数值热图
    for(i in 1:n_celltypes) {
      for(j in 1:n_tissues) {
        x <- (j - 0.5) / n_tissues
        y <- (n_celltypes - i + 0.5) / n_celltypes
        
        roe_val <- roe_matrix[i, j]
        
        if(is.finite(roe_val)) {
          bg_color <- color_func(roe_val)
          text_color <- ifelse(abs(roe_val - 1) > 0.5, "white", "black")
          
          grid.rect(x = x, y = y, 
                    width = cell_width * 0.9, height = cell_height * 0.9,
                    gp = gpar(fill = bg_color, col = "white", lwd = 1))
          grid.text(sprintf("%.2f", roe_val), x = x, y = y, 
                    gp = gpar(col = text_color, fontsize = 9))
        } else {
          grid.rect(x = x, y = y, 
                    width = cell_width * 0.9, height = cell_height * 0.9,
                    gp = gpar(fill = "#F0F0F0", col = "white", lwd = 1))
          grid.text("NA", x = x, y = y, gp = gpar(col = "gray50", fontsize = 8))
        }
      }
    }
    
    # 添加标签
    for(i in 1:n_celltypes) {
      y <- (n_celltypes - i + 0.5) / n_celltypes
      grid.text(celltypes[i], x = 1.02, y = y, just = "left", gp = gpar(fontsize = 9))
    }
    
    for(j in 1:n_tissues) {
      x <- (j - 0.5) / n_tissues
      grid.text(tissues[j], x = x, y = -0.02, just = "center", rot = 45, gp = gpar(fontsize = 9))
    }
    
    grid.text("Ro/e Index Values", x = 0.5, y = 1.05, just = "center",
              gp = gpar(fontsize = 14, fontface = "bold"))
    
    popViewport(2)
    popViewport()
    
    # === 添加总标题和图例 ===
    add_title_and_legend(title)
    
    dev.off()
    
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    message("Grid method failed: ", e$message, ". Trying base graphics...")
    create_base_plot(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes)
  })
}

#' 基础绘图系统方法 (兼容性备用)
create_base_plot <- function(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes) {
  
  tryCatch({
    n_tissues <- length(tissues)
    n_celltypes <- length(celltypes)
    
    pdf(output_path, width = max(12, n_tissues * 1.2), height = max(8, n_celltypes * 0.35))
    
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
    
    # 左侧：符号矩阵
    par(mar = c(8, 8, 4, 2))
    
    plot(1, 1, type = "n", xlim = c(0.5, n_tissues + 0.5), 
         ylim = c(0.5, n_celltypes + 0.5),
         xlab = "", ylab = "", main = "STARTRAC Symbols", 
         axes = FALSE)
    
    for(i in 1:n_celltypes) {
      for(j in 1:n_tissues) {
        symbol <- symbol_matrix[i, j]
        roe_val <- roe_matrix[i, j]
        
        # 设置颜色
        if(is.finite(roe_val)) {
          if(roe_val >= 2.0) { bg_col <- "#8B0000"; text_col <- "white" }
          else if(roe_val >= 1.2) { bg_col <- "#DC143C"; text_col <- "white" }
          else if(roe_val > 0.8) { bg_col <- "#FFFACD"; text_col <- "black" }
          else if(roe_val >= 0.5) { bg_col <- "#87CEEB"; text_col <- "black" }
          else { bg_col <- "#4682B4"; text_col <- "white" }
        } else {
          bg_col <- "#F0F0F0"; text_col <- "black"
        }
        
        rect(j - 0.4, i - 0.4, j + 0.4, i + 0.4, 
             col = bg_col, border = "white", lwd = 1.5)
        text(j, i, symbol, col = text_col, cex = 1.1, font = 2)
      }
    }
    
    axis(1, at = 1:n_tissues, labels = tissues, las = 2, cex.axis = 0.8)
    axis(2, at = 1:n_celltypes, labels = celltypes, las = 2, cex.axis = 0.8)
    
    # 右侧：数值图
    par(mar = c(8, 2, 4, 8))
    
    finite_vals <- roe_matrix[is.finite(roe_matrix)]
    if(length(finite_vals) > 0) {
      val_range <- range(finite_vals, na.rm = TRUE)
      colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(50)
      
      image(1:n_tissues, 1:n_celltypes, t(roe_matrix), 
            col = colors, axes = FALSE, main = "Ro/e Values", zlim = val_range)
      
      # 添加数值标签
      for(i in 1:n_celltypes) {
        for(j in 1:n_tissues) {
          if(is.finite(roe_matrix[i, j])) {
            text_color <- ifelse(abs(roe_matrix[i, j] - 1) > 0.5, "white", "black")
            text(j, i, sprintf("%.2f", roe_matrix[i, j]), col = text_color, cex = 0.8)
          }
        }
      }
      
      axis(1, at = 1:n_tissues, labels = tissues, las = 2, cex.axis = 0.8)
      axis(4, at = 1:n_celltypes, labels = celltypes, las = 2, cex.axis = 0.8)
    }
    
    # 添加总标题
    mtext(title, outer = TRUE, cex = 1.5, font = 2, line = -2)
    
    # 添加图例
    add_base_legend()
    
    dev.off()
    
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    stop("Both grid and base plotting methods failed: ", e$message)
  })
}

#' 添加标题和图例 (Grid版本)
add_title_and_legend <- function(title) {
  pushViewport(viewport())
  
  # 总标题
  grid.text(title, x = 0.5, y = 0.95, just = "center",
            gp = gpar(fontsize = 16, fontface = "bold"))
  
  # 图例
  legend_y <- 0.08
  legend_items <- c(
    "+++: Strong enrichment (Ro/e ≥ 2.0)",
    "+: Moderate enrichment (1.2 ≤ Ro/e < 2.0)", 
    "+/-: Near expected (0.8 < Ro/e < 1.2)",
    "-: Moderate depletion (0.5 ≤ Ro/e ≤ 0.8)",
    "---: Strong depletion (Ro/e < 0.5)"
  )
  
  legend_colors <- c("#8B0000", "#DC143C", "#FFFACD", "#87CEEB", "#4682B4")
  
  grid.text("Symbol Legend (STARTRAC-dist)", x = 0.5, y = legend_y + 0.03, just = "center",
            gp = gpar(fontsize = 12, fontface = "bold"))
  
  item_width <- 0.18
  for(i in 1:length(legend_items)) {
    x_pos <- 0.1 + (i - 1) * item_width
    
    grid.rect(x = x_pos, y = legend_y, width = 0.02, height = 0.02,
              gp = gpar(fill = legend_colors[i], col = "black"))
    grid.text(legend_items[i], x = x_pos + 0.025, y = legend_y, just = "left",
              gp = gpar(fontsize = 9))
  }
  
  grid.text("Method: STARTRAC-dist uses Ro/e = Observed/Expected from chi-square test", 
            x = 0.5, y = 0.02, just = "center", gp = gpar(fontsize = 10, col = "gray30"))
  
  popViewport()
}

#' 添加基础图例 (Base版本)
add_base_legend <- function() {
  par(fig = c(0, 1, 0, 1), oma = c(4, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  
  legend("bottom", 
         legend = c("+++: ≥2.0", "+: 1.2-2.0", "+/-: 0.8-1.2", "-: 0.5-0.8", "---: <0.5"),
         title = "STARTRAC-dist Symbol Legend",
         fill = c("#8B0000", "#DC143C", "#FFFACD", "#87CEEB", "#4682B4"),
         border = "black", horiz = TRUE, cex = 0.9, title.cex = 1.1)
}

# =============================================================================
# 主分析管道 (Main Analysis Pipeline)
# =============================================================================

#' 整合版组织分布分析 (主函数)
#' @param seurat_obj Seurat对象
#' @param celltype_col 细胞类型列名
#' @param tissue_col 组织位置列名  
#' @param patient_col 患者/样本列名
#' @param output_dir 输出目录
#' @param min_cells 最小细胞数阈值
#' @param statistical_test 统计检验方法
#' @param create_plots 是否创建可视化图像
#' @param plot_method 绘图方法选择
#' @return 分析结果列表
integrated_startrac_analysis <- function(
    seurat_obj,
    celltype_col,
    tissue_col, 
    patient_col,
    output_dir = NULL,
    min_cells = 50,
    statistical_test = "chisq",
    create_plots = TRUE,
    plot_method = "auto"
) {
  
  message("=================================================================")
  message("    INTEGRATED STARTRAC-DIST ANALYSIS PIPELINE")
  message("=================================================================")
  
  # 1. 初始化和验证
  if(is.null(output_dir)) {
    output_dir <- paste0("./startrac_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 检查必需列
  meta_data <- seurat_obj@meta.data
  required_cols <- c(celltype_col, tissue_col, patient_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta_data)]
  
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # 2. 数据预处理
  message("Step 1: Data preprocessing and quality control...")
  
  analysis_data <- data.frame(
    clone.id = rownames(meta_data),
    patient = as.character(meta_data[[patient_col]]),
    majorCluster = as.character(meta_data[[celltype_col]]),
    loc = as.character(meta_data[[tissue_col]]),
    Cell_Name = "Cell",
    stringsAsFactors = FALSE
  )
  
  # 数据清理
  original_cells <- nrow(analysis_data)
  analysis_data <- analysis_data[complete.cases(analysis_data), ]
  analysis_data <- analysis_data[
    analysis_data$majorCluster != "" & 
      analysis_data$loc != "" & 
      analysis_data$patient != "" &
      !analysis_data$majorCluster %in% c("NA", "Unknown", "unknown", "Unassigned") &
      !analysis_data$loc %in% c("NA", "Unknown", "unknown"), ]
  
  cells_after_cleaning <- nrow(analysis_data)
  removal_rate <- round((original_cells - cells_after_cleaning) / original_cells * 100, 2)
  
  message(sprintf("Data cleaning: %d → %d cells (%.1f%% removed)", 
                  original_cells, cells_after_cleaning, removal_rate))
  
  # 移除低频细胞类型
  celltype_counts <- table(analysis_data$majorCluster)
  low_count_celltypes <- names(celltype_counts)[celltype_counts < max(5, min_cells/10)]
  if(length(low_count_celltypes) > 0) {
    message("Removing low-count cell types: ", paste(head(low_count_celltypes, 3), collapse = ", "))
    analysis_data <- analysis_data[!analysis_data$majorCluster %in% low_count_celltypes, ]
  }
  
  # 验证数据充足性
  unique_tissues <- unique(analysis_data$loc)
  if(length(unique_tissues) < 2) {
    stop("Need at least 2 tissue locations. Found: ", paste(unique_tissues, collapse = ", "))
  }
  
  message(sprintf("Final dataset: %d cells, %d cell types, %d tissues, %d samples",
                  nrow(analysis_data), 
                  length(unique(analysis_data$majorCluster)),
                  length(unique_tissues),
                  length(unique(analysis_data$patient))))
  
  # 3. 计算Ro/e指数 (STARTRAC核心)
  message("Step 2: Calculating Ro/e indices using STARTRAC method...")
  
  roe_result <- tryCatch({
    calTissueDist(analysis_data,
                  byPatient = FALSE,
                  colname.cluster = "majorCluster",
                  colname.patient = "patient", 
                  colname.tissue = "loc",
                  method = statistical_test,
                  min.rowSum = min_cells)
  }, error = function(e) {
    stop("Failed to calculate Ro/e indices: ", e$message)
  })
  
  write.csv(roe_result, file.path(output_dir, "roe_indices.csv"))
  
  # 4. 统计显著性检验 (符合STARTRAC原文方法)
  message("Step 3: Statistical significance testing...")
  
  contingency_table <- table(analysis_data$majorCluster, analysis_data$loc)
  
  # 全局卡方检验
  global_chisq_result <- tryCatch({
    expected_from_margins <- outer(rowSums(contingency_table), colSums(contingency_table)) / sum(contingency_table)
    min_expected <- min(expected_from_margins)
    
    if(min_expected >= 1 && sum(expected_from_margins < 5) / length(expected_from_margins) <= 0.2) {
      chisq.test(contingency_table)
    } else {
      chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
    }
  }, error = function(e) {
    message("Warning in chi-square test: ", e$message)
    return(NULL)
  })
  
  # 处理统计结果
  significance_results <- NULL
  if(!is.null(global_chisq_result)) {
    global_significant <- global_chisq_result$p.value < 0.05
    
    message(sprintf("Global chi-square test: χ² = %.3f, p = %.6f (%s)", 
                    global_chisq_result$statistic, 
                    global_chisq_result$p.value,
                    ifelse(global_significant, "SIGNIFICANT", "NOT SIGNIFICANT")))
    
    # 细胞类型级别的显著性评估
    if(global_significant) {
      expected_matrix <- global_chisq_result$expected
      observed_matrix <- as.matrix(contingency_table)
      
      chi_contributions <- ((observed_matrix - expected_matrix)^2) / expected_matrix
      celltype_chi_contrib <- rowSums(chi_contributions)
      std_residuals <- (observed_matrix - expected_matrix) / sqrt(expected_matrix)
      max_abs_residuals <- apply(abs(std_residuals), 1, max)
      
      significance_results <- data.frame(
        CellType = rownames(roe_result),
        Chi_Contribution = celltype_chi_contrib,
        Max_Std_Residual = max_abs_residuals,
        Significant = celltype_chi_contrib > qchisq(0.95, df = ncol(contingency_table) - 1),
        stringsAsFactors = FALSE
      )
      
      significance_results$Significance_Level <- ifelse(
        significance_results$Max_Std_Residual > qnorm(0.9995), "***",
        ifelse(significance_results$Max_Std_Residual > qnorm(0.995), "**",
               ifelse(significance_results$Max_Std_Residual > qnorm(0.975), "*", "ns"))
      )
      
      write.csv(significance_results, file.path(output_dir, "significance_results.csv"), row.names = FALSE)
      
      n_significant <- sum(significance_results$Significant)
      message(sprintf("Significant cell types: %d/%d", n_significant, nrow(significance_results)))
    }
  }
  
  # 5. 创建可视化 (集成稳健版本)
  if(create_plots) {
    message("Step 4: Creating visualizations...")
    
    # 主要STARTRAC-dist图像
    create_startrac_dist_plot(
      roe_matrix = roe_result,
      significance_results = significance_results,
      output_path = file.path(output_dir, "startrac_dist_analysis.pdf"),
      title = "STARTRAC-dist Tissue Distribution Analysis",
      method = plot_method
    )
    
    # 补充可视化
    create_supplementary_plots(roe_result, analysis_data, output_dir)
  }
  
  # 6. 生成综合报告
  message("Step 5: Generating comprehensive report...")
  generate_analysis_report(roe_result, significance_results, analysis_data, output_dir)
  
  # 7. 返回结果
  results <- list(
    roe_indices = roe_result,
    significance_results = significance_results,
    raw_data = analysis_data,
    global_test = global_chisq_result,
    output_directory = output_dir,
    parameters = list(
      celltype_col = celltype_col,
      tissue_col = tissue_col,
      patient_col = patient_col,
      min_cells = min_cells,
      statistical_test = statistical_test
    )
  )
  
  message("=== ANALYSIS COMPLETED SUCCESSFULLY ===")
  message("Output directory: ", output_dir)
  message("Files generated: ", length(list.files(output_dir, recursive = TRUE)))
  
  return(results)
}

# =============================================================================
# 辅助分析和可视化函数 (Auxiliary Analysis & Visualization)
# =============================================================================

#' 创建增强版热图 (Enhanced Heatmap with Annotations)
#' @param roe_result Ro/e指数矩阵
#' @param significance_results 显著性结果
#' @param output_dir 输出目录
create_enhanced_roe_heatmap <- function(roe_result, significance_results, output_dir) {
  
  message("Creating enhanced Ro/e heatmap with annotations...")
  
  tryCatch({
    # 检查是否有ComplexHeatmap包
    if(!require("ComplexHeatmap", quietly = TRUE)) {
      message("ComplexHeatmap not available, creating alternative heatmap...")
      create_alternative_heatmap(roe_result, significance_results, output_dir)
      return()
    }
    
    # 准备细胞类型分组（基于名称模式）
    celltype_groups <- sapply(rownames(roe_result), function(x) {
      x_lower <- tolower(x)
      if(grepl("^t[_ ]|^cd[48]|tcr|th[0-9]|treg|cytotoxic", x_lower)) return("T_cells")
      if(grepl("^b[_ ]|plasma|^cd19|^cd20|antibody", x_lower)) return("B_cells")
      if(grepl("macro|monocyte|^m[0-9]|dendritic|dc[_ ]|^cd14|^cd16", x_lower)) return("Myeloid")
      if(grepl("fibroblast|smc|smooth.*muscle|pericyte|myofibro", x_lower)) return("Stromal")
      if(grepl("epithelial|endothelial|acinar|basal|ciliated", x_lower)) return("Epithelial")
      if(grepl("neutrophil|eosinophil|basophil|mast", x_lower)) return("Granulocytes")
      return("Other")
    })
    
    # 颜色函数
    col_fun <- colorRamp2(
      c(min(roe_result, na.rm = TRUE), 1, max(roe_result, na.rm = TRUE)), 
      c("#2166AC", "white", "#B2182B")
    )
    
    # 创建注释
    row_annotation <- NULL
    if(!is.null(significance_results)) {
      # 创建显著性注释向量
      sig_annotation <- rep("ns", nrow(roe_result))
      names(sig_annotation) <- rownames(roe_result)
      
      # 标记显著的细胞类型
      for(i in 1:nrow(significance_results)) {
        celltype <- significance_results$CellType[i]
        if(celltype %in% names(sig_annotation)) {
          sig_annotation[celltype] <- significance_results$Significance_Level[i]
        }
      }
      
      # 颜色方案
      sig_colors <- c(
        "***" = "#8B0000",  # 深红色 - 高显著性
        "**" = "#CD5C5C",   # 中红色 - 中等显著性 
        "*" = "#F0E68C",    # 淡黄色 - 轻度显著
        "ns" = "#FFFFFF"    # 白色 - 无显著性
      )
      
      # 细胞类型分组注释颜色
      group_colors <- c(
        "T_cells" = "#FF6B6B", "B_cells" = "#4ECDC4", "Myeloid" = "#45B7D1",
        "Stromal" = "#96CEB4", "Epithelial" = "#FFEAA7", "Granulocytes" = "#DDA0DD",
        "Other" = "#D3D3D3"
      )
      
      # 创建行注释
      row_annotation <- rowAnnotation(
        `Cell Group` = celltype_groups,
        `Significance` = sig_annotation,
        col = list(
          `Cell Group` = group_colors[names(group_colors) %in% unique(celltype_groups)],
          `Significance` = sig_colors
        ),
        annotation_legend_param = list(
          `Cell Group` = list(title = "Cell Type\nGroup"),
          `Significance` = list(title = "Statistical\nSignificance")
        ),
        width = unit(1.5, "cm")
      )
    }
    
    # 创建PDF文件
    pdf(file.path(output_dir, "enhanced_roe_heatmap_pairwise.pdf"), width = 14, height = max(8, nrow(roe_result) * 0.4))
    
    # 生成复杂热图
    ht <- Heatmap(
      as.matrix(roe_result),
      name = "Ro/e Index",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "right",
      column_names_side = "top",
      right_annotation = row_annotation,
      
      # 自定义单元格内容
      cell_fun = function(j, i, x, y, width, height, fill) {
        value <- roe_result[i, j]
        if(!is.na(value) && is.finite(value)) {
          # 添加Ro/e数值标签
          text_color <- ifelse(abs(value - 1) > 0.3, "white", "black")
          grid.text(sprintf("%.2f", value), x, y, 
                    gp = gpar(fontsize = 9, col = text_color))
          
          # 添加符号标记
          symbol <- roe_to_symbol(value)
          if(symbol %in% c("+++", "---")) {
            # 对于极端值添加小标记
            grid.circle(x + 0.35 * width, y + 0.35 * height, r = unit(1, "mm"), 
                        gp = gpar(fill = ifelse(symbol == "+++", "red", "blue"), 
                                  col = NA))
          }
        }
      },
      
      # 热图图例参数
      heatmap_legend_param = list(
        title = "Ro/e Index",
        at = c(min(roe_result, na.rm = TRUE), 1, max(roe_result, na.rm = TRUE)),
        labels = c("Depleted", "Expected", "Enriched"),
        legend_direction = "vertical",
        legend_width = unit(6, "cm"),
        title_position = "topcenter",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10)
      ),
      
      # 标题
      column_title = "Enhanced Tissue Distribution Analysis with Statistical Annotations",
      column_title_gp = gpar(fontsize = 16, fontface = "bold"),
      row_title = "Cell Types", 
      row_title_gp = gpar(fontsize = 14),
      
      # 添加边框
      border = TRUE,
      border_gp = gpar(col = "black", lty = 1)
    )
    
    # 绘制热图
    draw(ht)
    
    # 添加说明文字
    grid.text("Red/Blue dots: Extreme enrichment/depletion (Ro/e ≥2.0 or ≤0.5)", 
              x = 0.5, y = 0.05, 
              gp = gpar(fontsize = 10, col = "gray40"))
    grid.text("Significance levels: *** P<0.001, ** P<0.01, * P<0.05", 
              x = 0.5, y = 0.02, 
              gp = gpar(fontsize = 10, col = "gray40"))
    
    dev.off()
    
    message("Enhanced Ro/e heatmap saved successfully")
    
  }, error = function(e) {
    message("Failed to create enhanced heatmap: ", e$message)
    message("Creating alternative heatmap...")
    create_alternative_heatmap(roe_result, significance_results, output_dir)
  })
}

#' 创建备用热图 (当ComplexHeatmap不可用时)
create_alternative_heatmap <- function(roe_result, significance_results, output_dir) {
  
  tryCatch({
    pdf(file.path(output_dir, "enhanced_roe_heatmap_pairwise.pdf"), width = 12, height = 10)
    
    # 使用pheatmap作为备用
    if(require("pheatmap", quietly = TRUE)) {
      
      # 创建注释数据框
      annotation_row <- NULL
      if(!is.null(significance_results)) {
        annotation_row <- data.frame(
          Significance = rep("ns", nrow(roe_result)),
          row.names = rownames(roe_result)
        )
        
        for(i in 1:nrow(significance_results)) {
          celltype <- significance_results$CellType[i]
          if(celltype %in% rownames(annotation_row)) {
            annotation_row[celltype, "Significance"] <- significance_results$Significance_Level[i]
          }
        }
        
        # 注释颜色
        ann_colors <- list(
          Significance = c("***" = "#8B0000", "**" = "#CD5C5C", "*" = "#F0E68C", "ns" = "#FFFFFF")
        )
      }
      
      pheatmap(as.matrix(roe_result),
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               display_numbers = TRUE,
               number_format = "%.2f",
               fontsize_number = 8,
               annotation_row = annotation_row,
               annotation_colors = if(exists("ann_colors")) ann_colors else NULL,
               main = "Enhanced Ro/e Index Heatmap",
               cellwidth = 20,
               cellheight = 15)
      
    } else {
      # 最基本的热图
      heatmap(as.matrix(roe_result),
              col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
              main = "Ro/e Index Heatmap",
              margins = c(8, 8))
    }
    
    dev.off()
    message("Alternative heatmap created successfully")
    
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    message("Failed to create alternative heatmap: ", e$message)
  })
}

#' 创建补充可视化图像
create_supplementary_plots <- function(roe_result, analysis_data, output_dir, significance_results = NULL) {
  
  # 1. 创建增强版热图 (重要!)
  create_enhanced_roe_heatmap(roe_result, significance_results, output_dir)
  
  # 2. 柱状图
  tryCatch({
    roe_df <- as.data.frame(roe_result) %>%
      rownames_to_column("CellType") %>%
      pivot_longer(cols = -CellType, names_to = "Tissue", values_to = "Roe_Index") %>%
      filter(is.finite(Roe_Index))
    
    p_bar <- ggplot(roe_df, aes(x = CellType, y = Roe_Index, fill = Tissue)) +
      geom_col(position = "dodge", width = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
      scale_fill_viridis_d() +
      labs(title = "Cell Type Distribution Across Tissues (Ro/e Index)",
           subtitle = "Values > 1: Enriched, Values < 1: Depleted",
           x = "Cell Type", y = "Ro/e Index", fill = "Tissue Location") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, "roe_barplot.pdf"), p_bar, width = 12, height = 8)
  }, error = function(e) message("Warning: Failed to create bar plot - ", e$message))
  
  # 3. 细胞比例堆叠图
  tryCatch({
    prop_data <- analysis_data %>%
      group_by(loc, majorCluster) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(loc) %>%
      mutate(proportion = count / sum(count))
    
    p_stack <- ggplot(prop_data, aes(x = loc, y = proportion, fill = majorCluster)) +
      geom_col() +
      scale_fill_viridis_d() +
      labs(title = "Cell Type Composition by Tissue Location",
           x = "Tissue Location", y = "Proportion", fill = "Cell Type") +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, "cell_composition_stacked.pdf"), p_stack, width = 10, height = 6)
  }, error = function(e) message("Warning: Failed to create stacked plot - ", e$message))
  
  # 4. 相关性热图
  if(ncol(roe_result) >= 2 && nrow(roe_result) >= 2) {
    tryCatch({
      valid_roe_data <- roe_result[complete.cases(roe_result), ]
      if(nrow(valid_roe_data) >= 2) {
        correlation_matrix <- cor(t(valid_roe_data), use = "complete.obs", method = "pearson")
        
        pdf(file.path(output_dir, "celltype_correlation.pdf"), width = 10, height = 10)
        corrplot::corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust",
                           tl.cex = 0.8, tl.col = "black", tl.srt = 45,
                           title = "Cell Type Distribution Pattern Correlation", mar = c(0,0,2,0))
        dev.off()
      }
    }, error = function(e) message("Warning: Failed to create correlation plot - ", e$message))
  }
}

#' 生成综合分析报告
generate_analysis_report <- function(roe_result, significance_results, analysis_data, output_dir) {
  
  report_file <- file.path(output_dir, "comprehensive_analysis_report.txt")
  
  sink(report_file)
  
  cat("=======================================================\n")
  cat("    INTEGRATED STARTRAC-DIST ANALYSIS REPORT\n")
  cat("=======================================================\n\n")
  
  cat("Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total Cells Analyzed: ", nrow(analysis_data), "\n")
  cat("Cell Types: ", length(unique(analysis_data$majorCluster)), "\n")
  cat("Tissue Locations: ", length(unique(analysis_data$loc)), "\n")
  cat("Patients/Samples: ", length(unique(analysis_data$patient)), "\n\n")
  
  cat("STARTRAC-dist Method Summary:\n")
  cat("- Approach: Global chi-square test with cell-type decomposition\n")
  cat("- Ro/e Index: Ratio of Observed to Expected cell frequencies\n")
  cat("- Interpretation: >1 enriched, <1 depleted, =1 expected distribution\n\n")
  
  cat("Ro/e Index Results:\n")
  cat("-------------------\n")
  print(roe_result)
  cat("\n")
  
  if(!is.null(significance_results)) {
    sig_count <- sum(significance_results$Significant)
    cat("Statistical Significance Summary:\n")
    cat("---------------------------------\n")
    cat("Significant cell types: ", sig_count, "/", nrow(significance_results), 
        " (", round(sig_count/nrow(significance_results)*100, 1), "%)\n")
    
    if(sig_count > 0) {
      cat("\nTop significant cell types:\n")
      top_sig <- head(significance_results[significance_results$Significant, ], 5)
      for(i in 1:nrow(top_sig)) {
        cat("  ", i, ". ", top_sig$CellType[i], " (", top_sig$Significance_Level[i], ")\n", sep="")
      }
    }
  }
  
  cat("\nInterpretation Guidelines:\n")
  cat("- +++: Strong enrichment (Ro/e ≥ 2.0)\n")
  cat("- +: Moderate enrichment (1.2 ≤ Ro/e < 2.0)\n")
  cat("- +/-: Near expected (0.8 < Ro/e < 1.2)\n")
  cat("- -: Moderate depletion (0.5 ≤ Ro/e ≤ 0.8)\n")
  cat("- ---: Strong depletion (Ro/e < 0.5)\n\n")
  
  cat("Files Generated:\n")
  cat("- startrac_dist_analysis.pdf: Main STARTRAC-dist visualization\n")
  cat("- roe_indices.csv: Numerical Ro/e index matrix\n")
  cat("- significance_results.csv: Statistical significance testing results\n")
  cat("- roe_barplot.pdf: Bar chart visualization\n")
  cat("- cell_composition_stacked.pdf: Tissue composition analysis\n")
  cat("- celltype_correlation.pdf: Cell type pattern correlation\n")
  
  sink()
  
  message("Comprehensive report saved to: ", report_file)
}

# =============================================================================
# 便捷包装函数和使用示例 (Convenience Wrappers & Usage Examples)
# =============================================================================

#' 快速启动分析 (一键式分析)
#' @param seurat_obj Seurat对象
#' @param celltype_col 细胞类型列名
#' @param tissue_col 组织列名
#' @param patient_col 样本列名
#' @param output_dir 输出目录（可选）
quick_startrac_analysis <- function(seurat_obj, celltype_col, tissue_col, patient_col, output_dir = NULL) {
  
  # 显示数据概览
  message("=== DATA OVERVIEW ===")
  show_metadata_info(seurat_obj, show_preview = FALSE)
  
  # 执行分析
  results <- integrated_startrac_analysis(
    seurat_obj = seurat_obj,
    celltype_col = celltype_col,
    tissue_col = tissue_col,
    patient_col = patient_col,
    output_dir = output_dir,
    min_cells = 30,
    statistical_test = "chisq",
    create_plots = TRUE,
    plot_method = "auto"
  )
  
  # 快速结果解读
  message("\n=== QUICK RESULTS SUMMARY ===")
  interpret_results_summary(results)
  
  return(results)
}

#' 结果快速解读
interpret_results_summary <- function(results) {
  
  roe_matrix <- results$roe_indices
  
  # 找出最显著的模式
  max_enrichment <- apply(roe_matrix, 1, function(x) max(x, na.rm = TRUE))
  min_depletion <- apply(roe_matrix, 1, function(x) min(x, na.rm = TRUE))
  
  # 最富集的模式
  top_enriched_idx <- which.max(max_enrichment)
  if(length(top_enriched_idx) > 0) {
    ct <- names(max_enrichment)[top_enriched_idx]
    tissue <- names(which.max(roe_matrix[ct, ]))
    value <- max_enrichment[top_enriched_idx]
    message("Most enriched pattern: ", ct, " in ", tissue, " (Ro/e = ", round(value, 2), ")")
  }
  
  # 最耗竭的模式  
  top_depleted_idx <- which.min(min_depletion)
  if(length(top_depleted_idx) > 0) {
    ct <- names(min_depletion)[top_depleted_idx]
    tissue <- names(which.min(roe_matrix[ct, ]))
    value <- min_depletion[top_depleted_idx]
    message("Most depleted pattern: ", ct, " in ", tissue, " (Ro/e = ", round(value, 2), ")")
  }
  
  # 统计概览
  n_enriched <- sum(roe_matrix > 1.2, na.rm = TRUE)
  n_depleted <- sum(roe_matrix < 0.8, na.rm = TRUE)
  total_comparisons <- sum(is.finite(roe_matrix))
  
  message("Distribution patterns: ", n_enriched, " enriched, ", n_depleted, " depleted out of ", total_comparisons, " comparisons")
  
  if(!is.null(results$significance_results)) {
    n_sig <- sum(results$significance_results$Significant)
    message("Statistical significance: ", n_sig, "/", nrow(results$significance_results), " cell types show significant tissue distribution patterns")
  }
}

# =============================================================================
# 使用示例 (Usage Examples)
# =============================================================================

# 示例1：基本使用
results <- quick_startrac_analysis(
  seurat_obj = T_object,
  celltype_col = "Annotation",
  tissue_col = "tissue",
  patient_col = "sample"
)

# 示例2：完整参数控制
# results <- integrated_startrac_analysis(
#   seurat_obj = your_seurat_object,
#   celltype_col = "detailed_celltype",
#   tissue_col = "tissue_location",
#   patient_col = "patient_id", 
#   output_dir = "./my_startrac_analysis",
#   min_cells = 100,
#   statistical_test = "chisq",
#   create_plots = TRUE,
#   plot_method = "grid"
# )

# 示例3：仅创建可视化（已有数据）
# create_startrac_dist_plot(
#   roe_matrix = your_roe_matrix,
#   output_path = "custom_startrac_plot.pdf",
#   title = "Custom STARTRAC Analysis",
#   method = "auto"
# )

message("=================================================================")
message("    INTEGRATED STARTRAC-DIST ANALYSIS SYSTEM LOADED")
message("=================================================================")
message("Main functions:")
message("  • quick_startrac_analysis() - One-click analysis")
message("  • integrated_startrac_analysis() - Full parameter control") 
message("  • create_startrac_dist_plot() - Standalone visualization")
message("  • show_metadata_info() - Examine your data columns")
message("")
message("Key improvements:")
message("  ✓ Robust error handling with automatic fallback methods")
message("  ✓ STARTRAC-compliant statistical methodology") 
message("  ✓ High-quality publication-ready visualizations")
message("  ✓ Comprehensive analysis reports")
message("  ✓ Flexible parameter configuration")
message("")
message("Ready to analyze! Use quick_startrac_analysis() to get started.")
message("=================================================================")



######################################################################################################
######################################################################################################
######################################################################################################
#' 使用DimPlot创建多条件栅格化可视化
#' Create multi-condition rasterized visualization using DimPlot
#'
#' @param seurat_obj Seurat对象
#' @param condition_col 条件列名称
#' @param cell_type_col 细胞类型列名称  
#' @param reduction 降维方法
#' @param pt_size 点的大小
#' @param ncol 图形列数
#' @param raster 是否使用栅格化
#' @param alpha 点的透明度
#' @param label 是否显示聚类标签
#' @param label_size 标签大小
#' @param save_pdf 是否保存PDF文件
#' @param output_prefix 输出文件前缀
#' @return ggplot对象列表
create_multiCondition_dimplot <- function(seurat_obj,
                                          condition_col = "tissue",
                                          cell_type_col = "Annotation", 
                                          reduction = "umap",
                                          pt_size = 0.1,
                                          ncol = 3,
                                          raster = TRUE,
                                          alpha = 2,
                                          label = FALSE,
                                          label_size = 4,
                                          save_pdf = TRUE,
                                          output_prefix = "multiCondition_dimplot") {
  
  # 检查必需的列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  # 检查降维结果是否存在
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("降维结果", reduction, "不存在，请先运行", reduction, "分析"))
  }
  
  # 获取唯一的条件
  conditions <- unique(seurat_obj@meta.data[[condition_col]])
  conditions <- conditions[!is.na(conditions)]
  conditions <- sort(conditions)
  
  message(paste("发现", length(conditions), "个条件:", paste(conditions, collapse = ", ")))
  
  # 创建绘图列表
  plot_list <- list()
  
  # 设置细胞类型为idents用于DimPlot
  Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]
  
  # 1. 总体图：所有样本 (Overall plot: All samples)
  p_all <- DimPlot(seurat_obj, 
                   reduction = reduction,
                   raster = raster,
                   pt.size = pt_size,
                   alpha = alpha,
                   label = label,
                   label.size = label_size) +
    labs(title = "All samples") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      panel.grid = element_blank(),
      panel.border = element_blank()  # ✅ 去掉边框
    ) +
    guides(color = guide_legend(title = "Cell types", 
                                override.aes = list(size = 3, alpha = 1)))
  
  plot_list[["All_samples"]] <- p_all
  
  # 2. 各条件的子图 (Subplots for each condition)
  for (cond in conditions) {
    # 使用逻辑索引提取该条件的细胞
    condition_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[condition_col]] == cond]
    
    # 创建该条件的子集
    condition_subset <- subset(seurat_obj, cells = condition_cells)
    
    # 设置细胞类型为idents
    Idents(condition_subset) <- condition_subset@meta.data[[cell_type_col]]
    
    # 创建条件特异性图
    p_cond <- DimPlot(condition_subset, 
                      reduction = reduction,
                      raster = raster,
                      pt.size = pt_size,
                      alpha = alpha,
                      label = label,
                      label.size = label_size) +
      labs(title = cond) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none",  # 在子图中隐藏图例
        panel.grid = element_blank(),
        panel.border = element_blank() 
      )
    
    plot_list[[cond]] <- p_cond
  }
  
  # 3. 组合图形 (Combine plots)
  # 创建图例 (Create legend)
  legend_plot <- get_legend(p_all + 
                              theme(legend.position = "right",
                                    legend.direction = "vertical",
                                    legend.key.size = unit(0.6, "cm")))
  
  # 移除总体图的图例 (Remove legend from overall plot)
  p_all_no_legend <- p_all + theme(legend.position = "none")
  plot_list[["All_samples"]] <- p_all_no_legend
  
  # 计算布局 (Calculate layout)
  n_plots <- length(plot_list)
  
  # 创建组合图 (Create combined plot)
  if (n_plots <= 4) {
    combined_plot <- wrap_plots(plot_list, ncol = 2)
  } else {
    combined_plot <- wrap_plots(plot_list, ncol = ncol)
  }
  
  # 添加图例到右侧 (Add legend to the right)
  final_plot <- combined_plot | legend_plot
  final_plot <- final_plot + plot_layout(widths = c(4, 1))
  
  # 保存图形 (Save plots)
  if (save_pdf) {
    # 计算PDF尺寸 (Calculate PDF dimensions)
    pdf_width <- ifelse(n_plots <= 4, 12, 15)
    pdf_height <- ifelse(n_plots <= 4, 8, 10)
    
    # 保存组合图 (Save combined plot)
    ggsave(paste0(output_prefix, "_combined.pdf"), 
           final_plot, 
           width = pdf_width, 
           height = pdf_height, 
           dpi = 300)
    
    # 保存单独的图 (Save individual plots)
    for (plot_name in names(plot_list)) {
      ggsave(paste0(output_prefix, "_", plot_name, ".pdf"), 
             plot_list[[plot_name]], 
             width = 6, height = 5, dpi = 300)
    }
    
    message(paste("图形已保存:", paste0(output_prefix, "_combined.pdf")))
  }
  
  return(list(
    combined_plot = final_plot,
    individual_plots = plot_list,
    cell_counts = table(seurat_obj@meta.data[[condition_col]], 
                        seurat_obj@meta.data[[cell_type_col]])
  ))
}

# ====================================================================
# 诊断和检查函数 (Diagnostic and Check Functions)
# ====================================================================

#' 检查Seurat对象的数据结构
#' Check Seurat object data structure
#'
#' @param seurat_obj Seurat对象
#' @param condition_col 条件列名称
#' @param cell_type_col 细胞类型列名称
check_seurat_structure <- function(seurat_obj, 
                                   condition_col = "tissue", 
                                   cell_type_col = "Annotation") {
  
  cat("=== Seurat对象数据结构检查 ===\n")
  
  # 1. 基本信息
  cat("\n1. 基本信息:\n")
  cat(paste("   总细胞数:", ncol(seurat_obj), "\n"))
  cat(paste("   基因数:", nrow(seurat_obj), "\n"))
  
  # 2. 元数据列
  cat("\n2. 元数据列:\n")
  meta_cols <- colnames(seurat_obj@meta.data)
  cat(paste("   可用列:", paste(meta_cols, collapse = ", "), "\n"))
  
  # 3. 检查指定列是否存在
  cat("\n3. 指定列检查:\n")
  if (condition_col %in% meta_cols) {
    cat(paste("   ✓ 条件列 '", condition_col, "' 存在\n", sep = ""))
    conditions <- unique(seurat_obj@meta.data[[condition_col]])
    conditions <- conditions[!is.na(conditions)]
    cat(paste("     条件数量:", length(conditions), "\n"))
    cat(paste("     条件列表:", paste(conditions, collapse = ", "), "\n"))
  } else {
    cat(paste("   ✗ 条件列 '", condition_col, "' 不存在！\n", sep = ""))
    cat("     请检查列名或使用以下建议:\n")
    potential_cols <- meta_cols[grepl("tissue|group|condition|sample", meta_cols, ignore.case = TRUE)]
    if(length(potential_cols) > 0) {
      cat(paste("     可能的条件列:", paste(potential_cols, collapse = ", "), "\n"))
    }
  }
  
  if (cell_type_col %in% meta_cols) {
    cat(paste("   ✓ 细胞类型列 '", cell_type_col, "' 存在\n", sep = ""))
    cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])
    cell_types <- cell_types[!is.na(cell_types)]
    cat(paste("     细胞类型数量:", length(cell_types), "\n"))
    cat(paste("     细胞类型列表:", paste(head(cell_types, 10), collapse = ", ")))
    if(length(cell_types) > 10) cat(" ...")
    cat("\n")
  } else {
    cat(paste("   ✗ 细胞类型列 '", cell_type_col, "' 不存在！\n", sep = ""))
    cat("     请检查列名或使用以下建议:\n")
    potential_cols <- meta_cols[grepl("annotation|cell_type|cluster|ident", meta_cols, ignore.case = TRUE)]
    if(length(potential_cols) > 0) {
      cat(paste("     可能的细胞类型列:", paste(potential_cols, collapse = ", "), "\n"))
    }
  }
  
  # 4. 降维结果检查
  cat("\n4. 降维结果检查:\n")
  reductions <- names(seurat_obj@reductions)
  if(length(reductions) > 0) {
    cat(paste("   可用降维方法:", paste(reductions, collapse = ", "), "\n"))
  } else {
    cat("   ✗ 没有发现降维结果，请先运行RunTSNE()或RunUMAP()\n")
  }
  
  # 5. 细胞计数统计
  if (condition_col %in% meta_cols && cell_type_col %in% meta_cols) {
    cat("\n5. 细胞计数统计:\n")
    count_table <- table(seurat_obj@meta.data[[condition_col]], 
                         seurat_obj@meta.data[[cell_type_col]])
    print(count_table)
    
    # 检查是否有空组合
    empty_combinations <- which(count_table == 0, arr.ind = TRUE)
    if(nrow(empty_combinations) > 0) {
      cat("\n   警告: 发现空的组合 (没有细胞):\n")
      for(i in 1:nrow(empty_combinations)) {
        row_name <- rownames(count_table)[empty_combinations[i, 1]]
        col_name <- colnames(count_table)[empty_combinations[i, 2]]
        cat(paste("     ", row_name, " + ", col_name, "\n"))
      }
    }
  }
  
  cat("\n=== 检查完成 ===\n")
}

# ====================================================================
# 使用示例 (Usage Examples)
# ====================================================================

# 使用前先检查数据结构
# check_seurat_structure(T_object, "tissue", "Annotation")

# 基本调用示例
# results <- create_multiCondition_dimplot(
#   seurat_obj = T_object,
#   condition_col = "tissue",
#   cell_type_col = "Annotation", 
#   reduction = "tsne",
#   pt_size = 0.1,
#   raster = TRUE,    # 启用栅格化，类似FeaturePlot的raster=TRUE
#   alpha = 0.8,
#   save_pdf = TRUE,
#   output_prefix = "CRS_dimplot_analysis"
# )

# 高质量发表版本
# results_hq <- create_multiCondition_dimplot(
#   seurat_obj = T_object,
#   condition_col = "tissue",
#   cell_type_col = "Annotation",
#   pt_size = 0.2,     # 较大的点用于发表
#   raster = TRUE,
#   alpha = 0.9,       # 更不透明
#   label = TRUE,      # 显示细胞类型标签
#   label_size = 5,
#   output_prefix = "publication_quality"
# )

# 大数据集优化版本
# results_large <- create_multiCondition_dimplot(
#   seurat_obj = T_object,
#   condition_col = "tissue",
#   cell_type_col = "Annotation",
#   pt_size = 0.05,    # 很小的点避免重叠
#   raster = TRUE,     # 栅格化提高性能
#   alpha = 0.6,       # 较低透明度显示密度
#   output_prefix = "large_dataset_optimized"
# )

# ====================================================================
# 额外的便捷函数 (Additional Convenience Functions)
# ====================================================================

#' 快速预览函数
quick_preview <- function(seurat_obj, condition_col = "tissue", cell_type_col = "Annotation") {
  # 检查列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  results <- create_multiCondition_dimplot(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    cell_type_col = cell_type_col,
    pt_size = 0.05,
    raster = TRUE,
    alpha = 0.6,
    save_pdf = FALSE,  # 不保存文件，只预览
    output_prefix = "preview"
  )
  return(results$combined_plot)
}

#' 发表质量函数
publication_quality <- function(seurat_obj, condition_col = "tissue", cell_type_col = "Annotation") {
  # 检查列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  results <- create_multiCondition_dimplot(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    cell_type_col = cell_type_col,
    pt_size = 0.15,
    raster = TRUE,
    alpha = 0.9,
    label = FALSE,     # 通常发表不需要标签
    save_pdf = TRUE,
    output_prefix = "publication_dimplot"
  )
  return(results)
}

#' 带标签版本函数
labeled_version <- function(seurat_obj, condition_col = "tissue", cell_type_col = "Annotation") {
  # 检查列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  results <- create_multiCondition_dimplot(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    cell_type_col = cell_type_col,
    pt_size = 0.1,
    raster = TRUE,
    alpha = 0.7,
    label = TRUE,      # 显示细胞类型标签
    label_size = 4,
    save_pdf = TRUE,
    output_prefix = "labeled_dimplot"
  )
  return(results)
}

message("基于DimPlot的多条件栅格化可视化函数已加载完成!")
message("DimPlot-based multi-condition rasterized visualization functions loaded!")
message("现在使用Seurat原生的raster支持，无需额外包依赖")
message("Now using Seurat's native raster support, no additional package dependencies required")
message("已修复WhichCells错误，现在可以正常运行")
message("Fixed WhichCells error, should work properly now")
message("使用 check_seurat_structure(seurat_obj, 'tissue', 'Annotation') 检查数据")
message("Use check_seurat_structure(seurat_obj, 'tissue', 'Annotation') to check your data")

# 检查数据结构
check_seurat_structure(T_object, "tissue", "Annotation")

# 基本调用
results <- create_multiCondition_dimplot(
  seurat_obj = T_object,
  condition_col = "tissue",
  cell_type_col = "Annotation",
  reduction = 'umap',
  pt_size = 0.45,
  alpha = 2
)

# 快速预览
preview_plot <- quick_preview(T_object, "tissue", "Annotation")
print(preview_plot)

######################################################################################################
######################################################################################################
######################################################################################################
diagnosis <- diagnose_masc_convergence(T_object)

results <- scPairwiseMASCAnalysis(
  seurat_obj = T_object,
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  exclude_filter = "tissue_sampling_method == 'scraping'",
  fixed_effects_cols = NULL,  # 添加任何固定效应协变量
  output_dir = "pairwise_MASC_results_0928",
  # 如果您需要少于2个样本也能进行比较，可以降低此阈值
  min_samples = 2
)
props_data <- smart_prepare_plotting_data(
  masc_results = results,  # 先尝试从MASC结果提取
  seurat_obj = T_object,  # 备用重新计算
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  exclude_filter = "tissue_sampling_method != 'scraping'",  # 必须与MASC相同！
  min_samples = 2,
  min_cells = 10,
  min_prop = 0.01
)


# 第三步：验证比例数据
print("比例数据预览:")
head(props_data)

print("与MASC结果中的组织是否一致:")
print("绘图数据中的组织:")
table(props_data$Contrast)

print("MASC结果中的组织:")
if("all_results" %in% names(results)) {
  print(table(results$all_results$Comparison))
}

# 第四步：创建完美对齐的图表
final_result <- create_separate_aligned_plot(
  props_df = props_data,
  stats_results = results$all_results,
  output_dir = "consistent_plots",
  cell_types_order = c("Epithelial", "Endothelial", "Myeloid", "SMC", "T", "Fibroblast", "B"),
  contrast_colors = c(
    "lung parenchyma" = "#2E8B57",
    "sinus" = "#B22222", 
    "nose" = "#CD853F",
    "respiratory airway" = "#4682B4"
  )
)

print("箱子位置信息:")
print(final_result$positions)

######################################################################################################
######################################################################################################
######################################################################################################
# 统计每个组织的唯一样本数量
tissue_sample_counts <- sapply(unique(T_object$tissue), function(t) {
  length(unique(T_object$sample[T_object$tissue == t]))
})
names(tissue_sample_counts) <- unique(T_object$tissue)
print(tissue_sample_counts)

# 或者使用table + unique组合
unique_combinations <- unique(T_object@meta.data[, c("tissue", "sample")])
samples_per_tissue <- table(unique_combinations$tissue)
print(samples_per_tissue)


######################################################################################################
######################################################################################################
######################################################################################################

# 大型单细胞表达热图生成函数
# Comprehensive Single-Cell Expression Heatmap Generator

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

#' 创建大型单细胞表达热图（类似发表图片）
#' Create Large-Scale Single-Cell Expression Heatmap
#' 
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名，默认"Annotation"
#' @param genes_to_plot 要绘制的基因列表，默认NULL（自动选择）
#' @param max_cells_per_type 每种细胞类型最大细胞数，用于下采样
#' @param max_genes 最大基因数量
#' @param output_dir 输出目录
#' @param output_prefix 输出文件前缀
#' @param cell_type_colors 细胞类型颜色映射
#' @param expression_color_scheme 表达量颜色方案
#' @param figure_width 图片宽度（英寸）
#' @param figure_height 图片高度（英寸）
#' @param show_gene_names 是否显示基因名称
#' @param show_cell_names 是否显示细胞名称
#' @param cluster_genes 是否对基因进行聚类
#' @param cluster_cells 是否对细胞进行聚类
create_comprehensive_expression_heatmap <- function(
    seurat_obj,
    cell_type_col = "Annotation",
    genes_to_plot = NULL,
    max_cells_per_type = 500,
    max_genes = 200,
    output_dir = "expression_heatmaps",
    output_prefix = "comprehensive_expression",
    cell_type_colors = NULL,
    expression_color_scheme = c("blue", "white", "red"),
    figure_width = 20,
    figure_height = 16,
    show_gene_names = TRUE,
    show_cell_names = FALSE,
    cluster_genes = TRUE,
    cluster_cells = FALSE
) {
  
  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("开始创建大型表达热图 / Starting comprehensive expression heatmap creation...")
  
  # 1. 数据预处理和质控 ====
  message("Step 1: 数据预处理 / Data preprocessing...")
  
  # 检查细胞类型列是否存在
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "不存在于meta.data中"))
  }
  
  # 移除NA细胞类型
  valid_cells <- !is.na(seurat_obj@meta.data[[cell_type_col]])
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_cells])
  
  # 获取细胞类型信息
  cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])
  cell_types <- cell_types[!is.na(cell_types)]
  message(paste("检测到", length(cell_types), "种细胞类型:", paste(cell_types, collapse = ", ")))
  
  # 2. 基因选择策略 ====
  message("Step 2: 基因选择 / Gene selection...")
  
  if (is.null(genes_to_plot)) {
    # 自动选择基因：结合高变基因和marker基因
    
    # 方法1：获取每种细胞类型的top marker基因
    message("寻找细胞类型marker基因 / Finding cell type markers...")
    
    Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]
    
    tryCatch({
      # 寻找所有细胞类型的marker基因
      markers <- FindAllMarkers(
        seurat_obj, 
        only.pos = TRUE, 
        min.pct = 0.25, 
        logfc.threshold = 0.5,
        max.cells.per.ident = 1000  # 限制每种类型的细胞数以加速计算
      )
      
      # 每种细胞类型选择top基因
      top_markers <- markers %>%
        group_by(cluster) %>%
        top_n(n = min(15, max_genes %/% length(cell_types)), wt = avg_log2FC) %>%
        pull(gene) %>%
        unique()
      
      message(paste("找到", length(top_markers), "个marker基因"))
      
    }, error = function(e) {
      message("Marker基因寻找失败，使用高变基因替代")
      top_markers <- character(0)
    })
    
    # 方法2：高变基因
    if (length(VariableFeatures(seurat_obj)) > 0) {
      variable_genes <- head(VariableFeatures(seurat_obj), max_genes %/% 2)
      message(paste("选择", length(variable_genes), "个高变基因"))
    } else {
      message("计算高变基因...")
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = max_genes %/% 2)
      variable_genes <- VariableFeatures(seurat_obj)
    }
    
    # 合并基因列表
    genes_to_plot <- unique(c(top_markers, variable_genes))
    genes_to_plot <- head(genes_to_plot, max_genes)
    
  } else {
    # 使用提供的基因列表，并过滤存在的基因
    genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(seurat_obj)]
  }
  
  message(paste("最终选择", length(genes_to_plot), "个基因用于热图"))
  
  # 3. 细胞下采样 ====
  message("Step 3: 细胞下采样 / Cell subsampling...")
  
  # 按细胞类型进行分层抽样
  set.seed(42)  # 确保结果可重复
  
  sampled_cells <- c()
  for (ct in cell_types) {
    ct_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == ct]
    
    if (length(ct_cells) > max_cells_per_type) {
      # 随机抽样
      sampled_ct_cells <- sample(ct_cells, max_cells_per_type)
    } else {
      sampled_ct_cells <- ct_cells
    }
    
    sampled_cells <- c(sampled_cells, sampled_ct_cells)
    message(paste("细胞类型", ct, ":", length(sampled_ct_cells), "/", length(ct_cells), "个细胞"))
  }
  
  # 创建下采样后的子集
  heatmap_obj <- subset(seurat_obj, cells = sampled_cells)
  message(paste("下采样后总细胞数:", length(sampled_cells)))
  
  # 4. 数据标准化 ====
  message("Step 4: 数据标准化 / Data scaling...")
  
  # 确保选定的基因进行了标准化
  if (!all(genes_to_plot %in% rownames(heatmap_obj[["RNA"]]@scale.data))) {
    message("对选定基因进行标准化...")
    heatmap_obj <- ScaleData(heatmap_obj, features = genes_to_plot)
  }
  
  # 5. 准备热图数据 ====
  message("Step 5: 准备热图数据 / Preparing heatmap data...")
  
  # 提取标准化表达矩阵
  expression_matrix <- heatmap_obj[["RNA"]]@scale.data[genes_to_plot, ]
  
  # 处理无限值和NA值
  expression_matrix[is.infinite(expression_matrix)] <- 0
  expression_matrix[is.na(expression_matrix)] <- 0
  
  # 设置表达值范围（类似Seurat的默认设置）
  expression_matrix[expression_matrix > 2.5] <- 2.5
  expression_matrix[expression_matrix < -2.5] <- -2.5
  
  # 6. 细胞排序和分组 ====
  message("Step 6: 细胞排序 / Cell ordering...")
  
  # 按细胞类型对细胞进行分组和排序
  cell_metadata <- heatmap_obj@meta.data[colnames(expression_matrix), ]
  cell_metadata$CellType <- cell_metadata[[cell_type_col]]
  
  # 按细胞类型排序细胞
  if (!cluster_cells) {
    # 简单按细胞类型分组
    cell_order <- order(cell_metadata$CellType)
  } else {
    # 在每个细胞类型内进行聚类
    cell_order <- c()
    for (ct in sort(cell_types)) {
      ct_cells <- rownames(cell_metadata)[cell_metadata$CellType == ct]
      if (length(ct_cells) > 1) {
        ct_matrix <- expression_matrix[, ct_cells, drop = FALSE]
        if (ncol(ct_matrix) > 2) {
          # 计算细胞间距离并聚类
          cell_dist <- dist(t(ct_matrix))
          cell_hclust <- hclust(cell_dist)
          ct_order <- ct_cells[cell_hclust$order]
        } else {
          ct_order <- ct_cells
        }
      } else {
        ct_order <- ct_cells
      }
      cell_order <- c(cell_order, match(ct_order, colnames(expression_matrix)))
    }
  }
  
  # 重新排序表达矩阵
  expression_matrix <- expression_matrix[, cell_order]
  cell_metadata <- cell_metadata[colnames(expression_matrix), ]
  
  # 7. 设置颜色方案 ====
  message("Step 7: 设置颜色方案 / Setting color schemes...")
  
  # 细胞类型颜色
  if (is.null(cell_type_colors)) {
    if (length(cell_types) <= 12) {
      cell_type_colors <- brewer.pal(min(12, max(3, length(cell_types))), "Set3")
    } else {
      cell_type_colors <- rainbow(length(cell_types))
    }
    names(cell_type_colors) <- sort(cell_types)
  }
  
  # 表达量颜色
  if (length(expression_color_scheme) == 3) {
    expression_colors <- colorRamp2(
      c(-2.5, 0, 2.5), 
      expression_color_scheme
    )
  } else {
    expression_colors <- colorRamp2(
      seq(-2.5, 2.5, length.out = length(expression_color_scheme)),
      expression_color_scheme
    )
  }
  
  # 8. 创建注释 ====
  message("Step 8: 创建注释 / Creating annotations...")
  
  # 顶部细胞类型注释
  top_annotation <- HeatmapAnnotation(
    `Cell Type` = cell_metadata$CellType,
    col = list(`Cell Type` = cell_type_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    simple_anno_size = unit(8, "mm")
  )
  
  # 9. 生成热图 ====
  message("Step 9: 生成热图 / Generating heatmap...")
  
  # 基因聚类（如果需要）
  if (cluster_genes && nrow(expression_matrix) > 2) {
    gene_clustering <- TRUE
  } else {
    gene_clustering <- FALSE
  }
  
  # 创建主热图
  main_heatmap <- Heatmap(
    expression_matrix,
    name = "Scaled Expression",
    
    # 颜色设置
    col = expression_colors,
    
    # 聚类设置
    cluster_rows = gene_clustering,
    cluster_columns = FALSE,  # 列已经按细胞类型排序
    
    # 显示设置
    show_row_names = show_gene_names,
    show_column_names = show_cell_names,
    
    # 注释
    top_annotation = top_annotation,
    
    # 字体大小
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    
    # 热图大小
    width = unit(figure_width - 4, "inch"),
    height = unit(figure_height - 3, "inch"),
    
    # 添加分割线（在细胞类型之间）
    column_split = cell_metadata$CellType,
    column_gap = unit(1, "mm"),
    
    # 图例设置
    heatmap_legend_param = list(
      title = "Scaled\nExpression",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_height = unit(4, "cm")
    )
  )
  
  # 10. 保存热图 ====
  message("Step 10: 保存热图 / Saving heatmap...")
  
  # PDF格式（矢量图，适合期刊）
  pdf_file <- file.path(output_dir, paste0(output_prefix, "_heatmap.pdf"))
  pdf(pdf_file, width = figure_width, height = figure_height)
  draw(main_heatmap, 
       annotation_legend_side = "right",
       heatmap_legend_side = "right")
  dev.off()
  
  # PNG格式（高分辨率位图）
  png_file <- file.path(output_dir, paste0(output_prefix, "_heatmap.png"))
  png(png_file, width = figure_width * 300, height = figure_height * 300, 
      res = 300, bg = "white")
  draw(main_heatmap,
       annotation_legend_side = "right", 
       heatmap_legend_side = "right")
  dev.off()
  
  # 11. 保存相关数据 ====
  message("Step 11: 保存相关数据 / Saving related data...")
  
  # 保存基因列表
  write.csv(
    data.frame(Gene = rownames(expression_matrix)),
    file.path(output_dir, paste0(output_prefix, "_genes.csv")),
    row.names = FALSE
  )
  
  # 保存细胞信息
  write.csv(
    cell_metadata,
    file.path(output_dir, paste0(output_prefix, "_cell_metadata.csv")),
    row.names = TRUE
  )
  
  # 保存表达矩阵（可选，如果数据量不大）
  if (ncol(expression_matrix) * nrow(expression_matrix) < 1e6) {
    write.csv(
      expression_matrix,
      file.path(output_dir, paste0(output_prefix, "_expression_matrix.csv")),
      row.names = TRUE
    )
  }
  
  message("热图创建完成！/ Heatmap creation completed!")
  message(paste("输出文件保存在:", output_dir))
  message(paste("主要文件:", basename(pdf_file), "和", basename(png_file)))
  
  # 返回关键信息
  result <- list(
    expression_matrix = expression_matrix,
    cell_metadata = cell_metadata,
    genes_used = rownames(expression_matrix),
    cell_types = cell_types,
    cell_type_colors = cell_type_colors,
    output_files = c(pdf_file, png_file)
  )
  
  return(result)
}

#' 快速预览热图（小规模版本）
#' Quick heatmap preview with reduced scale
quick_preview_heatmap <- function(seurat_obj,
                                  cell_type_col = "Annotation",
                                  max_cells_total = 1000,
                                  max_genes = 50) {
  
  message("创建快速预览热图...")
  
  # 简化参数创建预览
  result <- create_comprehensive_expression_heatmap(
    seurat_obj = seurat_obj,
    cell_type_col = cell_type_col,
    max_cells_per_type = max_cells_total %/% length(unique(seurat_obj@meta.data[[cell_type_col]])),
    max_genes = max_genes,
    output_prefix = "preview",
    figure_width = 12,
    figure_height = 8,
    show_gene_names = TRUE
  )
  
  return(result)
}

# 使用示例 / Usage Examples:
#
# 1. 基本使用（自动选择基因）
result <- create_comprehensive_expression_heatmap(
  seurat_obj = T_object,
  cell_type_col = "Annotation",
  output_dir = "publication_heatmaps"
)
#
# # 2. 使用特定基因列表
# marker_genes <- c("CD3E", "CD3D", "CD14", "CD68", "EPCAM", "PECAM1", "COL1A1")
# result <- create_comprehensive_expression_heatmap(
#   seurat_obj = T_object,
#   cell_type_col = "Annotation", 
#   genes_to_plot = marker_genes,
#   max_cells_per_type = 300,
#   output_prefix = "marker_expression"
# )
#
# # 3. 快速预览
# preview <- quick_preview_heatmap(T_object)
#
# # 4. 自定义颜色
# custom_colors <- c(
#   "Epithelial" = "#E41A1C",
#   "Fibroblast" = "#377EB8", 
#   "T" = "#4DAF4A",
#   "Myeloid" = "#984EA3",
#   "Endothelial" = "#FF7F00"
# )
# result <- create_comprehensive_expression_heatmap(
#   seurat_obj = T_object,
#   cell_type_col = "Annotation",
#   cell_type_colors = custom_colors,
#   expression_color_scheme = c("#2166AC", "#F7F7F7", "#B2182B")
# )

######################################################################################################
######################################################################################################
######################################################################################################

#########################################
# 组织间细胞通讯对比分析主脚本
#########################################

# 加载所需的依赖包
library(Seurat)
library(CellChat)
library(patchwork)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(future)
library(grid)
library(igraph)

# 设置分析参数
output_dir <- "CellChat_Tissue_Comparison"
group_by <- "Annotation"  # 细胞类型注释列
tissue_column <- "tissue"  # 组织类型列
top_n_pathways <- 10  # 要详细分析的顶级通路数量

# 创建输出目录
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

# 读取预处理好的Seurat对象
# showWarnings = TRUE <- readRDS("final_analyzed_showWarnings = TRUE.rds")

# 获取要比较的组织类型
all_tissues <- unique(T_object[[tissue_column]])
message(paste("将分析以下组织:", paste(all_tissues, collapse = ", ")))
str(all_tissues)
# 将all_tissues从数据框转换为字符向量
all_tissues <- as.character(all_tissues$tissue)
print(all_tissues)  # 检查转换结果
#-------------------------
# 第一步：为每个组织创建和分析CellChat对象
#-------------------------
# 创建CellChat对象列表
cellchat_list <- list()

for(tissue in all_tissues) {
  message(paste("处理组织:", tissue))
  
  # 创建安全的目录名
  safe_tissue_name <- gsub(" ", "_", tissue)
  safe_tissue_name <- gsub("[^a-zA-Z0-9_]", "", safe_tissue_name)
  
  # 构建完整路径
  tissue_dir <- file.path(output_dir, safe_tissue_name)
  dir.create(tissue_dir, recursive = TRUE, showWarnings = TRUE)
  
  # 检查是否已有保存的对象
  saved_file <- file.path(tissue_dir, "cellchat_obj.rds")
  
  if(file.exists(saved_file)) {
    message(paste("加载已有的CellChat对象:", saved_file))
    cellchat_obj <- readRDS(saved_file)
  } else {
    # 提取组织特异的细胞
    message(paste0("提取组织: ", tissue))
    cells_to_keep <- which(T_object@meta.data[[tissue_column]] == tissue)
    if(length(cells_to_keep) == 0) {
      warning(paste0("在组织列中找不到匹配的组织: ", tissue))
      next
    }
    tissue_seurat <- T_object[, cells_to_keep]
    
    # 创建CellChat对象
    message("创建CellChat对象...")
    cellchat_obj <- createCellChat(object = tissue_seurat, 
                                   group.by = group_by, 
                                   assay = "RNA")
    
    # 加载数据库
    message("加载数据库...")
    CellChatDB <- CellChatDB.human  # 人类数据
    # CellChatDB <- CellChatDB.mouse  # 小鼠数据(取消注释选择)
    
    # 使用"Secreted Signaling"子集
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
    cellchat_obj@DB <- CellChatDB.use
    
    # 预处理
    message("数据预处理...")
    cellchat_obj <- subsetData(cellchat_obj)
    cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
    cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
    
    # 计算通讯概率
    message("计算细胞通讯概率...")
    cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean")
    cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
    
    # 计算通路水平通讯
    message("计算信号通路水平通讯...")
    cellchat_obj <- computeCommunProbPathway(cellchat_obj)
    
    # 聚合网络
    message("聚合通讯网络...")
    cellchat_obj <- aggregateNet(cellchat_obj)
    
    # 计算网络中心性
    message("计算网络中心性...")
    cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
    
    # 保存结果
    saveRDS(cellchat_obj, file = saved_file)
  }
  
  # 保存到列表
  cellchat_list[[tissue]] <- cellchat_obj
  
  message(paste("完成", tissue, "的CellChat对象创建和分析"))
}

#-------------------------
# 第二步：合并CellChat对象并进行比较分析
#-------------------------
# 合并对象用于比较分析
message("合并CellChat对象...")
# 确保所有对象都是CellChat类
all_cellchat <- all(sapply(cellchat_list, function(x) inherits(x, "CellChat")))
if(!all_cellchat) {
  stop("输入列表中含有非CellChat对象")
}

# 合并CellChat对象
merged_cellchat <- mergeCellChat(cellchat_list, 
                                 add.names = names(cellchat_list))
# merged_cellchat <- readRDS("E:\\R\\0301\\CellChat_Tissue_Comparison\\CellChat_Tissue_Comparison\\CellChat_Tissue_Comparison\\merged_cellchat.rds")
# 保存合并的对象
saveRDS(merged_cellchat, file.path(output_dir, "merged_cellchat.rds"))

#-------------------------
# 第三步：比较组织间通讯数量和强度
#-------------------------
message("比较组织间通讯数量和强度...")

# 修改为
cell_counts <- sapply(cellchat_list, function(x) {
  # 直接使用meta.data中的细胞数量
  if(length(x@meta$datasets) > 0) {
    return(length(x@meta$datasets))
  } else if(length(x@idents) > 0) {
    # 备选方法：使用idents向量的长度
    return(length(x@idents))
  } else {
    # 备选方法：使用表达矩阵的列数
    return(ncol(x@data))
  }
})

# 创建比较图
# 数量比较图
gg1 <- compareInteractions(merged_cellchat, 
                           show.legend = F, 
                           group = names(cellchat_list),
                           measure = "count",
                           title = "组织间通讯数量比较")

# 强度比较图
gg2 <- compareInteractions(merged_cellchat, 
                           show.legend = F, 
                           group = names(cellchat_list),
                           measure = "weight",
                           title = "组织间通讯强度比较")

# 增加每个数据集的摘要统计
stats_data <- data.frame(
  Tissue = names(cellchat_list),
  Cell_Count = cell_counts,
  Int_Count = sapply(cellchat_list, function(x) sum(x@net$count != 0)),
  Int_Weight = sapply(cellchat_list, function(x) sum(x@net$weight)),
  Unique_Pathways = sapply(cellchat_list, function(x) length(x@netP$pathways))
)

# 计算每千细胞的通讯数
stats_data$Int_Per_1k_Cells <- stats_data$Int_Count / (stats_data$Cell_Count/1000)

# 添加表格到图形
if(requireNamespace("gridExtra", quietly = TRUE)) {
  gt <- gridExtra::tableGrob(stats_data, rows = NULL)
  gg3 <- ggplot() + 
    theme_void() +
    annotation_custom(gt, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    labs(title = "组织通讯统计摘要")
  
  # 组合图形
  if(requireNamespace("patchwork", quietly = TRUE)) {
    p_combined <- (gg1 + gg2) / gg3
  } else {
    p_combined <- gg1 + gg2
  }
} else {
  p_combined <- gg1 + gg2
}

# 保存比较图
ggsave(file.path(output_dir, "interaction_comparison.pdf"), 
       p_combined, width = 12, height = 7)

# 保存统计数据
write.csv(stats_data, 
          file = file.path(output_dir, "interaction_stats.csv"), 
          row.names = FALSE)

#-------------------------
# 第四步：比较细胞类型间通讯差异
#-------------------------
message("比较细胞类型间通讯差异...")

# 创建输出目录
diff_dir <- file.path(output_dir, "Diff_Interactions")
dir.create(diff_dir, showWarnings = TRUE, recursive = TRUE)

# 基于互作数量比较
message("基于互作数量的差异比较...")
pdf(file.path(diff_dir, "diff_interaction_count.pdf"), width = 12, height = 10)
netVisual_diffInteraction(merged_cellchat, weight.scale = TRUE, measure = "count")
dev.off()

# 差异互作热图(数量)
pdf(file.path(diff_dir, "diff_heatmap_count.pdf"), width = 12, height = 10)
gg_count <- netVisual_heatmap(merged_cellchat, measure = "count")
ComplexHeatmap::draw(gg_count)
dev.off()

# 基于互作强度比较
message("基于互作强度的差异比较...")
pdf(file.path(diff_dir, "diff_interaction_weight.pdf"), width = 12, height = 10)
netVisual_diffInteraction(merged_cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

# 差异互作热图(强度)
pdf(file.path(diff_dir, "diff_heatmap_weight.pdf"), width = 12, height = 10)
gg_weight <- netVisual_heatmap(merged_cellchat, measure = "weight")
ComplexHeatmap::draw(gg_weight)
dev.off()

# 各组织网络图比较
message("生成各组织网络图...")
weight_max_count <- getMaxWeight(cellchat_list, attribute = c("idents", "count"))
weight_max_weight <- getMaxWeight(cellchat_list, attribute = c("idents", "weight"))

pdf(file.path(diff_dir, "tissue_count_comparison.pdf"), width = 14, height = 12)
par(mfrow = c(ceiling(length(cellchat_list)/2), 2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(cellchat_list[[i]]@net$count, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.weight.max = weight_max_count[2], 
                   edge.width.max = 12, 
                   title.name = paste0(names(cellchat_list)[i], " - 互作数量"))
}
dev.off()

pdf(file.path(diff_dir, "tissue_weight_comparison.pdf"), width = 14, height = 12)
par(mfrow = c(ceiling(length(cellchat_list)/2), 2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(cellchat_list[[i]]@net$weight, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.weight.max = weight_max_weight[2], 
                   edge.width.max = 12, 
                   title.name = paste0(names(cellchat_list)[i], " - 互作强度"))
}
dev.off()

#-------------------------
# 第五步：比较信号通路强度
#-------------------------
message("比较组织间信号通路强度...")

# 从merged_cellchat对象获取各组织的通路
pathways_all <- list()
for(tissue in names(merged_cellchat@netP)) {
  if("pathways" %in% names(merged_cellchat@netP[[tissue]])) {
    pathways_all[[tissue]] <- names(merged_cellchat@netP[[tissue]][["pathways"]])
  }
}

# 检查是否成功获取
print(lapply(pathways_all, length))

# 获取共有通路
if(length(pathways_all) > 0) {
  common_pathways <- Reduce(intersect, pathways_all)
  message(paste0("共有通路数量: ", length(common_pathways)))
  
  # 获取所有唯一通路
  all_unique_pathways <- unique(unlist(pathways_all))
  message(paste0("唯一通路数量: ", length(all_unique_pathways)))
}

# 获取共有通路
pathways_all <- lapply(cellchat_list, function(x) names(x@netP$pathways))
common_pathways <- Reduce(intersect, pathways_all)

message(paste0("发现", length(common_pathways), "个所有组织共有的信号通路"))

# 1. 两两比较所有组织对的信号通路强度
tissue_pairs <- combn(names(cellchat_list), 2, simplify = FALSE)
message(paste0("创建", length(tissue_pairs), "个组织对的通路强度比较..."))

# 为每对组织创建比较图
for(i in seq_along(tissue_pairs)) {
  pair <- tissue_pairs[[i]]
  pair_name <- paste(pair, collapse = "_vs_")
  message(paste0("比较组织对: ", pair_name))
  
  # 只选择这两个组织的cellchat对象
  cellchat_pair <- list()
  cellchat_pair[[pair[1]]] <- cellchat_list[[pair[1]]]
  cellchat_pair[[pair[2]]] <- cellchat_list[[pair[2]]]
  
  # 合并这两个组织的cellchat对象
  merged_pair <- mergeCellChat(cellchat_pair, add.names = names(cellchat_pair))
  
  # 创建比较图
  pair_dir <- file.path(output_dir, "Pairwise_Comparisons")
  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 尝试使用rankNet进行比较
  tryCatch({
    gg1 <- rankNet(merged_pair, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    gg2 <- rankNet(merged_pair, mode = "comparison", stacked = FALSE, do.stat = TRUE)
    
    # 组合并保存
    p <- gg1 + gg2
    ggsave(file.path(pair_dir, paste0("pathway_comparison_", pair_name, ".pdf")), 
           p, width = 12, height = 8)
  }, error = function(e) {
    message(paste0("为组织对 ", pair_name, " 创建rankNet图时出错: ", e$message))
  })
}

# 2. 创建所有组织的通路强度综合比较
if(length(common_pathways) > 0) {
  # 计算每个组织中每个通路的强度
  pathway_strength <- sapply(common_pathways, function(pathway) {
    sapply(names(cellchat_list), function(tissue) {
      if(pathway %in% names(cellchat_list[[tissue]]@netP$pathways)) {
        sum(cellchat_list[[tissue]]@netP$prob[, pathway])
      } else {
        0
      }
    })
  })
  
  # 转换为数据框
  df <- as.data.frame(pathway_strength)
  rownames(df) <- names(cellchat_list)
  
  # 保存通路强度数据
  write.csv(df, file = file.path(output_dir, "pathway_strength_data.csv"))
  
  # 计算每个通路在所有组织中的总强度，用于排序
  total_strength <- colSums(df)
  top_pathways <- names(sort(total_strength, decreasing = TRUE))
  
  # 限制显示的通路数量（如果太多）
  top_n_pathways <- min(length(top_pathways), 30)  # 最多显示30个通路
  top_pathways <- top_pathways[1:top_n_pathways]
  
  # 热图可视化 - 使用所有组织的数据
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    # 按总强度排序的通路热图
    heatmap_data <- t(df)[top_pathways, ]
    
    # 绘制热图并保存
    pdf(file.path(output_dir, "pathway_strength_heatmap.pdf"), width = 12, height = 10)
    pheatmap::pheatmap(heatmap_data,
                       cluster_rows = FALSE,  # 保持按强度排序
                       cluster_cols = TRUE,   # 组织可以聚类
                       display_numbers = TRUE,
                       fontsize_number = 8,
                       main = "组织间信号通路强度比较")
    dev.off()
    
    # 创建标准化热图（按行标准化，突出组织间差异）
    heatmap_data_scaled <- t(apply(heatmap_data, 1, function(x) {
      if(max(x) > 0) {
        return(x / max(x))  # 每行除以最大值，归一化到0-1
      } else {
        return(x)
      }
    }))
    
    pdf(file.path(output_dir, "pathway_strength_heatmap_scaled.pdf"), width = 12, height = 10)
    pheatmap::pheatmap(heatmap_data_scaled,
                       cluster_rows = FALSE,
                       cluster_cols = TRUE,
                       display_numbers = TRUE,
                       fontsize_number = 8,
                       main = "组织间信号通路强度比较（按行归一化）")
    dev.off()
  }
  
  # 3. 创建自定义条形图展示顶级通路在各组织中的比较
  if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
    # 使用前15个通路创建条形图
    top_15_pathways <- top_pathways[1:min(15, length(top_pathways))]
    plot_data <- df[, top_15_pathways, drop = FALSE]
    
    # 转换为长格式数据
    plot_data_long <- reshape2::melt(as.matrix(plot_data),
                                     varnames = c("Tissue", "Pathway"),
                                     value.name = "Strength")
    
    # 设置因子水平以控制绘图顺序
    plot_data_long$Pathway <- factor(plot_data_long$Pathway, levels = rev(top_15_pathways))
    
    # 创建条形图
    bar_plot <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = Pathway, y = Strength, fill = Tissue)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::coord_flip() +  # 水平条形图
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::labs(
        title = "主要信号通路在不同组织中的强度比较",
        x = "通路",
        y = "通路强度"
      )
    
    # 保存条形图
    ggplot2::ggsave(file.path(output_dir, "top_pathways_barplot.pdf"), 
                    bar_plot, width = 12, height = 10)
    
    # 4. 创建雷达图比较各组织的通路分布
    if(length(names(cellchat_list)) > 2 && requireNamespace("fmsb", quietly = TRUE)) {
      # 准备雷达图数据 - 使用前10个通路
      top_10_pathways <- top_pathways[1:min(10, length(top_pathways))]
      radar_data <- as.data.frame(t(df[, top_10_pathways, drop = FALSE]))
      
      # 添加最大值和最小值行，fmsb包需要
      radar_data <- rbind(
        rep(max(radar_data), ncol(radar_data)),  # 最大值行
        rep(0, ncol(radar_data)),                # 最小值行
        radar_data
      )
      
      # 创建雷达图
      pdf(file.path(output_dir, "pathway_radar_chart.pdf"), width = 10, height = 10)
      par(mar = c(1, 1, 3, 1))  # 调整边距
      fmsb::radarchart(
        radar_data,
        pcol = rainbow(ncol(radar_data)),
        pfcol = sapply(rainbow(ncol(radar_data)), function(x) adjustcolor(x, alpha.f = 0.3)),
        plwd = 2,
        cglcol = "grey",
        cglty = 1,
        axislabcol = "grey30",
        caxislabels = seq(0, round(max(radar_data[-(1:2),])), length.out = 5),
        title = "组织信号通路分布雷达图"
      )
      legend(
        "topright",
        legend = colnames(radar_data),
        col = rainbow(ncol(radar_data)),
        lty = 1,
        lwd = 2,
        bty = "n"
      )
      dev.off()
    }
  }
  # 这部分代码需要修改，处理无共有通路的情况
} else {
  message("未找到所有组织共有的通路，尝试使用组织特异通路进行比较...")
  
  # 获取在任意组织中出现的所有通路
  all_unique_pathways <- unique(unlist(pathways_all))
  
  # 创建存在矩阵 - 显示每个通路在哪些组织中存在
  presence_matrix <- sapply(all_unique_pathways, function(pathway) {
    sapply(names(cellchat_list), function(tissue) {
      if(pathway %in% names(cellchat_list[[tissue]]@netP$pathways)) {
        return(1)
      } else {
        return(0)
      }
    })
  })
  
  # 确保presence_matrix是一个矩阵，而不是向量
  if(!is.matrix(presence_matrix)) {
    # 如果只有一个组织或一个通路，需要手动转换为矩阵
    if(length(all_unique_pathways) == 1) {
      # 只有一个通路
      presence_matrix <- matrix(presence_matrix, 
                                nrow = length(names(cellchat_list)), 
                                ncol = 1,
                                dimnames = list(names(cellchat_list), all_unique_pathways))
    } else if(length(names(cellchat_list)) == 1) {
      # 只有一个组织
      presence_matrix <- matrix(presence_matrix, 
                                nrow = 1, 
                                ncol = length(all_unique_pathways),
                                dimnames = list(names(cellchat_list), all_unique_pathways))
    }
  }
  
  # 计算每个通路在多少组织中存在
  # 使用rowSums如果矩阵被转置了，或者colSums如果未转置
  if(is.matrix(presence_matrix)) {
    pathway_counts <- colSums(presence_matrix)
    
    # 按存在组织数量排序
    sorted_pathways <- names(sort(pathway_counts, decreasing = TRUE))
    
    # 绘制通路存在热图
    if(requireNamespace("pheatmap", quietly = TRUE)) {
      pdf(file.path(output_dir, "pathway_presence_heatmap.pdf"), width = 12, height = 10)
      pheatmap::pheatmap(t(presence_matrix)[sorted_pathways, ],
                         cluster_rows = FALSE,
                         cluster_cols = TRUE,
                         color = c("white", "darkblue"),
                         legend_breaks = c(0, 1),
                         legend_labels = c("缺失", "存在"),
                         display_numbers = TRUE,
                         fontsize_number = 8,
                         main = "信号通路在各组织中的存在情况")
      dev.off()
    }
    
    # 创建通路在各组织中的出现频率条形图
    if(requireNamespace("ggplot2", quietly = TRUE)) {
      # 将数据转换为数据框
      freq_data <- data.frame(
        Pathway = names(pathway_counts),
        Frequency = pathway_counts / length(names(cellchat_list)) * 100  # 转换为百分比
      )
      
      # 按频率排序
      freq_data <- freq_data[order(-freq_data$Frequency), ]
      
      # 限制显示的通路数量
      if(nrow(freq_data) > 30) {
        freq_data <- freq_data[1:30, ]
      }
      
      # 设置因子水平以控制绘图顺序
      freq_data$Pathway <- factor(freq_data$Pathway, levels = rev(freq_data$Pathway))
      
      # 创建条形图
      bar_plot <- ggplot2::ggplot(freq_data, ggplot2::aes(x = Pathway, y = Frequency)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +  # 水平条形图
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 10),
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
          title = "信号通路在组织中的出现频率",
          x = "通路",
          y = "出现频率 (%)"
        )
      
      # 保存条形图
      ggplot2::ggsave(file.path(output_dir, "pathway_frequency_barplot.pdf"), 
                      bar_plot, width = 12, height = 10)
    }
  } else {
    message("无法创建存在矩阵, 跳过热图生成")
    
    # 输出一些调试信息
    message("all_unique_pathways的长度: ", length(all_unique_pathways))
    message("cellchat_list的名称: ", paste(names(cellchat_list), collapse=", "))
    message("presence_matrix的类型: ", class(presence_matrix))
    message("presence_matrix的结构: ")
    str(presence_matrix)
  }
}
message("完成信号通路强度比较分析")
#-------------------------
# Step 6: Compare signaling sending and receiving patterns
#-------------------------
message("Comparing signaling sending and receiving patterns across tissues...")

# Create output directory
pattern_dir <- file.path(output_dir, "Signaling_Patterns")
dir.create(pattern_dir, showWarnings = TRUE, recursive = TRUE)

# Merge all pathways
i <- 1
pathway.union <- union(cellchat_list[[i]]@netP$pathways, 
                       cellchat_list[[i+1]]@netP$pathways)

# 对每个模式生成热图
for(pattern in c("outgoing", "incoming", "all")) {
  message(paste0("Generating heatmap for ", pattern, " pattern..."))
  
  # 设置热图颜色
  color <- switch(pattern,
                  "outgoing" = "Reds",
                  "incoming" = "GnBu",
                  "all" = "OrRd")
  
  # 为每个组织单独创建PDF文件
  for(i in seq_along(cellchat_list)) {
    tissue <- names(cellchat_list)[i]
    message(paste0("Processing tissue: ", tissue))
    
    # 创建单独的文件名
    safe_tissue_name <- gsub(" ", "_", tissue)
    safe_tissue_name <- gsub("[^a-zA-Z0-9_]", "", safe_tissue_name)
    file_path <- file.path(pattern_dir, paste0("signaling_pattern_", pattern, "_", safe_tissue_name, ".pdf"))
    
    # 尝试生成热图
    tryCatch({
      pdf(file_path, width = 10, height = 8)
      ht <- netAnalysis_signalingRole_heatmap(cellchat_list[[i]], 
                                              pattern = pattern, 
                                              signaling = pathway.union, 
                                              title = tissue, 
                                              width = 8, 
                                              height = 10, 
                                              color.heatmap = color)
      
      # 如果成功创建则绘制热图
      if(!is.null(ht)) {
        ComplexHeatmap::draw(ht)
      }
      dev.off()
    }, error = function(e) {
      # 确保即使出错也关闭设备
      if(dev.cur() > 1) dev.off()
      message(paste0("Error generating heatmap for tissue ", tissue, ": ", e$message))
    })
  }
  
  # 也创建一个汇总文件
  file_path_all <- file.path(pattern_dir, paste0("signaling_pattern_", pattern, "_all_tissues.pdf"))
  pdf(file_path_all, width = 14, height = 12)
  
  # 构建一个简单的文本说明页面
  grid::grid.newpage()
  grid::grid.text(paste0("Signaling Pattern: ", pattern, "\n\n",
                         "Individual heatmaps have been saved as separate files for each tissue."),
                  gp = grid::gpar(fontsize = 14))
  
  dev.off()
}
#-------------------------
# Step 7: Analyze pathway similarity (enhanced version)
#-------------------------
message("Analyzing pathway similarity...")

# Create output directory
similarity_dir <- file.path(output_dir, "Pathway_Similarity")
dir.create(similarity_dir, showWarnings = TRUE, recursive = TRUE)

#-------------------------------------------------------
# First check tissue data validity and pathway overlap
#-------------------------------------------------------
# Examine pathways in each tissue
pathway_check_file <- file.path(similarity_dir, "pathway_check.txt")
sink(pathway_check_file)

cat("Pathway Analysis Report\n")
cat("======================\n\n")

# Extract pathways for each tissue
pathways_by_tissue <- list()
for(tissue in names(cellchat_list)) {
  n_pathways <- length(cellchat_list[[tissue]]@netP$pathways)
  cat(paste0("Tissue: ", tissue, ", Number of pathways: ", n_pathways, "\n"))
  
  # Store pathways for this tissue
  if(n_pathways > 0) {
    pathways_by_tissue[[tissue]] <- names(cellchat_list[[tissue]]@netP$pathways)
    top_paths <- pathways_by_tissue[[tissue]][1:min(5, n_pathways)]
    cat(paste0("  Top pathways: ", paste(top_paths, collapse=", "), "\n"))
  } else {
    pathways_by_tissue[[tissue]] <- character(0)
    cat("  No pathways found\n")
  }
  cat("\n")
}

# Identify valid tissues with pathways
valid_tissues <- names(cellchat_list)[sapply(cellchat_list, function(x) {
  length(x@netP$pathways) > 0
})]
cat(paste0("Valid tissues with pathways: ", paste(valid_tissues, collapse=", "), "\n\n"))

# Calculate pathway overlap
if(length(valid_tissues) >= 2) {
  cat("Pathway Overlap Analysis:\n")
  cat("=========================\n\n")
  
  # Find common pathways across all tissues
  common_pathways <- Reduce(intersect, pathways_by_tissue[valid_tissues])
  cat(paste0("Pathways common to all tissues: ", length(common_pathways), "\n"))
  if(length(common_pathways) > 0) {
    cat(paste0("  ", paste(common_pathways, collapse=", "), "\n"))
  }
  cat("\n")
  
  # Find pathways present in at least 2 tissues
  all_pathways <- unique(unlist(pathways_by_tissue))
  pathway_presence <- sapply(all_pathways, function(p) {
    sum(sapply(pathways_by_tissue, function(x) p %in% x))
  })
  
  shared_pathways <- names(pathway_presence[pathway_presence >= 2])
  cat(paste0("Pathways present in at least 2 tissues: ", length(shared_pathways), "\n"))
  if(length(shared_pathways) > 0 && length(shared_pathways) <= 20) {
    cat(paste0("  ", paste(shared_pathways, collapse=", "), "\n"))
  } else if(length(shared_pathways) > 20) {
    cat(paste0("  Top 20: ", paste(shared_pathways[1:20], collapse=", "), "...\n"))
  }
  cat("\n")
  
  # Check unique pathways for each tissue
  cat("Unique pathways by tissue:\n")
  for(tissue in valid_tissues) {
    unique_paths <- setdiff(pathways_by_tissue[[tissue]], 
                            unlist(pathways_by_tissue[setdiff(valid_tissues, tissue)]))
    cat(paste0("  ", tissue, ": ", length(unique_paths), " unique pathways\n"))
    if(length(unique_paths) > 0 && length(unique_paths) <= 10) {
      cat(paste0("    ", paste(unique_paths, collapse=", "), "\n"))
    } else if(length(unique_paths) > 10) {
      cat(paste0("    Top 10: ", paste(unique_paths[1:10], collapse=", "), "...\n"))
    }
  }
}

sink()

#-------------------------------------------------------
# Attempt standard CellChat similarity analysis
#-------------------------------------------------------
standard_similarity_success <- FALSE

if(length(valid_tissues) >= 2) {
  message("Attempting standard CellChat similarity analysis...")
  
  # If we have valid tissues, create a subset of the merged object
  if(length(valid_tissues) < length(cellchat_list)) {
    cellchat_subset <- cellchat_list[valid_tissues]
    merged_subset <- mergeCellChat(cellchat_subset, add.names = names(cellchat_subset))
  } else {
    cellchat_subset <- cellchat_list  # Define cellchat_subset even when using all tissues
    merged_subset <- merged_cellchat
  }
  
  # First, ensure UMAP is available
  if(!requireNamespace("reticulate", quietly = TRUE)) {
    message("Installing reticulate package...")
    install.packages("reticulate")
  }
  
  # tryCatch({
  #   message("Checking UMAP installation...")
  #   reticulate::py_module_available("umap")
  # }, error = function(e) {
  #   message("UMAP not found. Attempting to install...")
  #   reticulate::py_install(packages = "umap-learn")
  # })
  
  # Continue with structural similarity analysis
  tryCatch({
    message("Computing structural similarity...")
    merged_subset <- computeNetSimilarityPairwise(merged_subset, type = "structural")
    
    # Try to use a different dimensionality reduction if UMAP fails
    tryCatch({
      merged_subset <- netEmbedding(merged_subset, type = "structural")
    }, error = function(e) {
      message("Error in netEmbedding with default method. Trying with alternative method...")
      # Try with PCA instead
      merged_subset <- netEmbedding(merged_subset, type = "structural", reduction.method = "pca")
    })
    
    merged_subset <- netClustering(merged_subset, type = "structural")
    
    # Check if computation was successful
    if(!is.null(merged_subset@net$structural.similarity)) {
      standard_similarity_success <- TRUE
      message("Standard structural similarity computation successful")
      
      # Visualization
      pdf(file.path(similarity_dir, "standard_structural_similarity.pdf"), width = 12, height = 10)
      
      # Plot embedding
      netVisual_embeddingPairwise(merged_subset, type = "structural", label.size = 3.5)
      
      # Plot similarity ranking
      p <- rankSimilarity(merged_subset, type = "structural")
      print(p)
      
      # Plot similarity heatmap
      sim_mat <- merged_subset@net$structural.similarity
      if(!is.null(sim_mat) && requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(sim_mat, 
                           color = colorRampPalette(c("blue", "white", "red"))(100),
                           display_numbers = TRUE, 
                           number_format = "%.2f",
                           fontsize_number = 8,
                           main = "Structural Similarity Matrix")
      }
      
      dev.off()
      
      # Save similarity matrix
      write.csv(sim_mat, 
                file = file.path(similarity_dir, "structural_similarity_matrix.csv"))
    }
  }, error = function(e) {
    message(paste0("Standard structural similarity analysis failed: ", e$message))
  })
  
  # Attempt functional similarity analysis if all tissues have the same cell types
  same_celltypes <- all(sapply(cellchat_subset[-1], function(x) {
    all(levels(x@idents) %in% levels(cellchat_subset[[1]]@idents)) &&
      all(levels(cellchat_subset[[1]]@idents) %in% levels(x@idents))
  }))
  
  if(same_celltypes) {
    message("Attempting functional similarity analysis...")
    tryCatch({
      merged_subset <- computeNetSimilarityPairwise(merged_subset, type = "functional")
      merged_subset <- netEmbedding(merged_subset, type = "functional")
      merged_subset <- netClustering(merged_subset, type = "functional")
      
      # Check if computation was successful
      if(!is.null(merged_subset@net$functional.similarity)) {
        message("Standard functional similarity computation successful")
        
        pdf(file.path(similarity_dir, "standard_functional_similarity.pdf"), width = 12, height = 10)
        
        # Plot embedding
        netVisual_embeddingPairwise(merged_subset, type = "functional", label.size = 3.5)
        
        # Plot similarity ranking
        p <- rankSimilarity(merged_subset, type = "functional")
        print(p)
        
        # Plot similarity heatmap
        func_sim_mat <- merged_subset@net$functional.similarity
        if(!is.null(func_sim_mat) && requireNamespace("pheatmap", quietly = TRUE)) {
          pheatmap::pheatmap(func_sim_mat, 
                             color = colorRampPalette(c("blue", "white", "red"))(100),
                             display_numbers = TRUE, 
                             number_format = "%.2f",
                             fontsize_number = 8,
                             main = "Functional Similarity Matrix")
        }
        
        dev.off()
        
        # Save similarity matrix
        write.csv(func_sim_mat, 
                  file = file.path(similarity_dir, "functional_similarity_matrix.csv"))
      }
    }, error = function(e) {
      message(paste0("Standard functional similarity analysis failed: ", e$message))
    })
  } else {
    message("Tissues have different cell type compositions, skipping functional similarity analysis")
  }
}

#-------------------------------------------------------
# Create custom similarity matrices (alternative approach)
#-------------------------------------------------------
message("Creating custom similarity matrices based on pathway overlap...")

# Custom function to compute Jaccard similarity based on pathway overlap
create_pathway_overlap_similarity <- function(cellchat_list) {
  tissues <- names(cellchat_list)
  n_tissues <- length(tissues)
  
  # Extract pathways for each tissue
  pathways_by_tissue <- lapply(cellchat_list, function(x) {
    if(length(x@netP$pathways) > 0) {
      return(names(x@netP$pathways))
    } else {
      return(character(0))
    }
  })
  
  # Create similarity matrix
  sim_mat <- matrix(0, nrow = n_tissues, ncol = n_tissues)
  rownames(sim_mat) <- tissues
  colnames(sim_mat) <- tissues
  
  # Compute Jaccard similarity for each tissue pair
  for(i in 1:n_tissues) {
    for(j in 1:n_tissues) {
      if(i == j) {
        sim_mat[i, j] <- 1  # Self-similarity is 1
      } else {
        paths_i <- pathways_by_tissue[[tissues[i]]]
        paths_j <- pathways_by_tissue[[tissues[j]]]
        
        # Compute Jaccard index only if both tissues have pathways
        if(length(paths_i) > 0 && length(paths_j) > 0) {
          intersection <- length(intersect(paths_i, paths_j))
          union_size <- length(union(paths_i, paths_j))
          sim_mat[i, j] <- intersection / union_size
        }
      }
    }
  }
  
  return(sim_mat)
}

# Custom function to compute similarity based on pathway strength correlation
create_pathway_strength_similarity <- function(cellchat_list) {
  tissues <- names(cellchat_list)
  n_tissues <- length(tissues)
  
  # Get all unique pathways across tissues
  all_pathways <- unique(unlist(lapply(cellchat_list, function(x) {
    if(length(x@netP$pathways) > 0) {
      return(names(x@netP$pathways))
    } else {
      return(character(0))
    }
  })))
  
  # Skip if no pathways found
  if(length(all_pathways) == 0) {
    return(NULL)
  }
  
  # Create a matrix of pathway strengths
  pathway_strength <- matrix(0, nrow = n_tissues, ncol = length(all_pathways))
  rownames(pathway_strength) <- tissues
  colnames(pathway_strength) <- all_pathways
  
  # Fill in pathway strengths
  for(i in 1:n_tissues) {
    tissue <- tissues[i]
    
    if(length(cellchat_list[[tissue]]@netP$pathways) > 0) {
      for(pathway in all_pathways) {
        if(pathway %in% names(cellchat_list[[tissue]]@netP$pathways)) {
          # Sum up the pathway probabilities
          pathway_strength[i, pathway] <- sum(cellchat_list[[tissue]]@netP$prob[, pathway])
        }
      }
    }
  }
  
  # Compute correlation-based similarity
  sim_mat <- matrix(0, nrow = n_tissues, ncol = n_tissues)
  rownames(sim_mat) <- tissues
  colnames(sim_mat) <- tissues
  
  for(i in 1:n_tissues) {
    for(j in 1:n_tissues) {
      if(i == j) {
        sim_mat[i, j] <- 1  # Self-similarity is 1
      } else {
        # Compute correlation only if both tissues have pathways
        if(sum(pathway_strength[i,]) > 0 && sum(pathway_strength[j,]) > 0) {
          sim_mat[i, j] <- cor(pathway_strength[i,], pathway_strength[j,], method = "spearman")
        }
      }
    }
  }
  
  return(list(
    similarity = sim_mat,
    pathway_strength = pathway_strength
  ))
}

# Create and visualize Jaccard similarity matrix
jaccard_sim <- create_pathway_overlap_similarity(cellchat_list)
if(!is.null(jaccard_sim) && requireNamespace("pheatmap", quietly = TRUE)) {
  # Save matrix as CSV
  write.csv(jaccard_sim, file = file.path(similarity_dir, "jaccard_similarity_matrix.csv"))
  
  # Create heatmap
  pdf(file.path(similarity_dir, "jaccard_similarity_heatmap.pdf"), width = 10, height = 8)
  pheatmap::pheatmap(jaccard_sim, 
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     display_numbers = TRUE, 
                     number_format = "%.2f",
                     fontsize_number = 8,
                     main = "Pathway Overlap Similarity (Jaccard Index)")
  dev.off()
  
  # Create hierarchical clustering visualization
  pdf(file.path(similarity_dir, "jaccard_similarity_dendrogram.pdf"), width = 10, height = 8)
  # First create a distance matrix (1 - similarity)
  dist_mat <- as.dist(1 - jaccard_sim)
  # Then perform hierarchical clustering
  hclust_result <- hclust(dist_mat, method = "complete")
  # Plot dendrogram
  plot(hclust_result, main = "Hierarchical Clustering of Tissues by Pathway Overlap",
       xlab = "", sub = "")
  dev.off()
}

# Create and visualize pathway strength similarity
strength_sim_result <- create_pathway_strength_similarity(cellchat_list)
if(!is.null(strength_sim_result) && requireNamespace("pheatmap", quietly = TRUE)) {
  strength_sim <- strength_sim_result$similarity
  pathway_strength <- strength_sim_result$pathway_strength
  
  # Save matrices as CSV
  write.csv(strength_sim, file = file.path(similarity_dir, "strength_correlation_matrix.csv"))
  write.csv(pathway_strength, file = file.path(similarity_dir, "pathway_strength_by_tissue.csv"))
  
  # Create correlation heatmap
  pdf(file.path(similarity_dir, "strength_correlation_heatmap.pdf"), width = 10, height = 8)
  pheatmap::pheatmap(strength_sim, 
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     display_numbers = TRUE, 
                     number_format = "%.2f",
                     fontsize_number = 8,
                     main = "Pathway Strength Correlation Similarity")
  dev.off()
  
  # Create hierarchical clustering visualization
  pdf(file.path(similarity_dir, "strength_correlation_dendrogram.pdf"), width = 10, height = 8)
  # First create a distance matrix (1 - similarity)
  dist_mat <- as.dist(1 - strength_sim)
  # Then perform hierarchical clustering
  hclust_result <- hclust(dist_mat, method = "complete")
  # Plot dendrogram
  plot(hclust_result, main = "Hierarchical Clustering of Tissues by Pathway Strength",
       xlab = "", sub = "")
  dev.off()
  
  # Create heatmap of pathway strengths
  # Select top pathways by total strength
  pathway_totals <- colSums(pathway_strength)
  top_pathways <- names(sort(pathway_totals, decreasing = TRUE))[1:min(30, length(pathway_totals))]
  
  # Create heatmap of top pathways
  pdf(file.path(similarity_dir, "pathway_strength_heatmap.pdf"), width = 12, height = 10)
  pheatmap::pheatmap(pathway_strength[, top_pathways], 
                     color = colorRampPalette(c("white", "navy"))(100),
                     display_numbers = FALSE, 
                     fontsize_row = 10,
                     fontsize_col = 8,
                     main = "Pathway Strength by Tissue")
  dev.off()
}

message("Pathway similarity analysis completed")

#-------------------------
# 第八步：分析信号通路差异
#-------------------------
message("分析组织间信号通路差异...")

# 创建输出目录
diff_path_dir <- file.path(output_dir, "Differential_Pathways")
dir.create(diff_path_dir, showWarnings = TRUE, recursive = TRUE)

# 两两组织间比较
tissues_pairs <- combn(all_tissues, 2, simplify = FALSE)
for(pair in tissues_pairs) {
  message(paste("比较", pair[1], "和", pair[2], "的差异通路"))
  pair_name <- paste(pair, collapse = "_vs_")
  pair_dir <- file.path(diff_path_dir, pair_name)
  dir.create(pair_dir, showWarnings = TRUE, recursive = TRUE)
  
  # 正向分析 (第一个组织 vs 第二个组织)
  features.name1 <- paste0(pair[1], ".vs.", pair[2])
  
  # 差异表达分析
  merged_cellchat <- identifyOverExpressedGenes(merged_cellchat, 
                                                group.dataset = "datasets", 
                                                pos.dataset = pair[1], 
                                                neg.dataset = pair[2],
                                                features.name = features.name1, 
                                                only.pos = FALSE, 
                                                thresh.pc = 0.1, 
                                                thresh.fc = 0.25,
                                                thresh.p = 0.05)
  
  # 映射差异表达结果到通讯网络
  net1 <- netMappingDEG(merged_cellchat, features.name = features.name1)
  
  # 提取上调的配体-受体对
  net1.up <- subsetCommunication(merged_cellchat, net = net1, 
                                 datasets = pair[1],
                                 ligand.logFC = 0.25, 
                                 receptor.logFC = NULL)
  
  # 提取下调的配体-受体对
  net1.down <- subsetCommunication(merged_cellchat, net = net1, 
                                   datasets = pair[2],
                                   ligand.logFC = 0.25, 
                                   receptor.logFC = NULL)
  
  # 保存差异通路数据
  if(nrow(net1.up) > 0) {
    write.csv(net1.up, file = file.path(pair_dir, paste0(pair[1], "_upregulated_LR.csv")))
  }
  if(nrow(net1.down) > 0) {
    write.csv(net1.down, file = file.path(pair_dir, paste0(pair[2], "_upregulated_LR.csv")))
  }
  
  # 可视化上调和下调的信号
  # 1. 气泡图
  if(nrow(net1.up) > 0 || nrow(net1.down) > 0) {
    pdf(file.path(pair_dir, "differential_lr_bubble.pdf"), width = 12, height = 10)
    
    if(nrow(net1.up) > 0) {
      pairLR.up <- net1.up[, "interaction_name", drop = FALSE]
      gg1 <- netVisual_bubble(merged_cellchat, 
                              pairLR.use = pairLR.up, 
                              comparison = c(pair[1], pair[2]), 
                              angle.x = 90, 
                              remove.isolate = TRUE,
                              title.name = paste0("Up-regulated in ", pair[1]))
      print(gg1)
    }
    
    if(nrow(net1.down) > 0) {
      pairLR.down <- net1.down[, "interaction_name", drop = FALSE]
      gg2 <- netVisual_bubble(merged_cellchat, 
                              pairLR.use = pairLR.down, 
                              comparison = c(pair[2], pair[1]), 
                              angle.x = 90, 
                              remove.isolate = TRUE,
                              title.name = paste0("Up-regulated in ", pair[2]))
      print(gg2)
    }
    
    dev.off()
  }
  
  # 2. 和弦图
  pdf(file.path(pair_dir, "differential_lr_chord.pdf"), width = 12, height = 10)
  par(mfrow = c(1,2), xpd=TRUE)
  
  if(nrow(net1.up) > 0) {
    netVisual_chord_gene(cellchat_list[[pair[1]]], 
                         slot.name = 'net', 
                         net = net1.up, 
                         lab.cex = 0.8, 
                         small.gap = 3.5, 
                         title.name = paste0("Up-regulated in ", pair[1]))
  }
  
  if(nrow(net1.down) > 0) {
    netVisual_chord_gene(cellchat_list[[pair[2]]], 
                         slot.name = 'net', 
                         net = net1.down, 
                         lab.cex = 0.8, 
                         small.gap = 3.5, 
                         title.name = paste0("Up-regulated in ", pair[2]))
  }
  
  dev.off()
  
  # 3. 富集分析
  pdf(file.path(pair_dir, "enriched_ligands_wordcloud.pdf"), width = 10, height = 8)
  if(nrow(net1.up) > 0) {
    enriched.up <- computeEnrichmentScore(net1.up, species = 'human', 
                                          variable.both = TRUE)
  }
  if(nrow(net1.down) > 0) {
    enriched.down <- computeEnrichmentScore(net1.down, species = 'human', 
                                            variable.both = TRUE)
  }
  dev.off()
  
  # 反向分析已通过上面的down结果完成，无需再次计算
  
  message(paste("完成", pair[1], "和", pair[2], "的差异通路分析"))
}

#-------------------------
# 第九步：比较特定信号通路
#-------------------------
message("比较特定信号通路...")

# 选择要分析的通路
if(length(common_pathways) > 0) {
  message(paste0("所有组织共有 ", length(common_pathways), " 个信号通路"))
  
  # 计算通路在所有组织中的总强度
  pathway_strengths <- sapply(common_pathways, function(pathway) {
    sum(sapply(cellchat_list, function(x) sum(x@netP$prob[, pathway])))
  })
  
  # 选择前N个强度最高的通路
  selected_pathways <- names(sort(pathway_strengths, decreasing = TRUE)[1:min(top_n_pathways, length(pathway_strengths))])
} else {
  message("未找到所有组织共有的信号通路，将使用出现频率最高的通路")
  
  # 获取并集和出现频率
  all_pathways <- unlist(lapply(cellchat_list, function(x) names(x@netP$pathways)))
  pathway_counts <- table(all_pathways)
  
  # 选择出现频率最高的通路
  selected_pathways <- names(sort(pathway_counts, decreasing = TRUE)[1:min(top_n_pathways, length(pathway_counts))])
}

message(paste0("选择以下通路进行详细比较: ", paste(selected_pathways, collapse = ", ")))

# 比较选择的通路
pathway_dir <- file.path(output_dir, "Pathway_Comparison")
dir.create(pathway_dir, showWarnings = TRUE, recursive = TRUE)

# 处理每个信号通路
for(pathway in selected_pathways) {
  message(paste0("分析通路: ", pathway))
  pathway_dir_specific <- file.path(pathway_dir, pathway)
  dir.create(pathway_dir_specific, showWarnings = TRUE, recursive = TRUE)
  
  # 获取最大权重以统一尺度
  weight.max <- getMaxWeight(cellchat_list, 
                             slot.name = c("netP"), 
                             attribute = pathway)
  
  # 使用不同布局可视化
  for(layout in c("circle", "chord", "hierarchy")) {
    pdf(file.path(pathway_dir_specific, paste0(pathway, "_", layout, ".pdf")), 
        width = 12, height = 10)
    
    par(mfrow = c(1, length(cellchat_list)), xpd=TRUE)
    
    for(i in seq_along(cellchat_list)) {
      tissue <- names(cellchat_list)[i]
      cellchat_obj <- cellchat_list[[i]]
      
      # 检查通路是否存在
      if(!pathway %in% names(cellchat_obj@netP$pathways)) {
        message(paste("通路", pathway, "在", tissue, "中不存在，跳过"))
        next
      }
      
      if(layout == "hierarchy") {
        # 定义层次图的接收者顶点
        n_idents <- length(levels(cellchat_obj@idents))
        vertex.receiver <- seq(1, ceiling(n_idents/2))
        
        netVisual_aggregate(cellchat_obj, 
                            signaling = pathway, 
                            vertex.receiver = vertex.receiver,
                            edge.weight.max = weight.max[1], 
                            edge.width.max = 10, 
                            signaling.name = paste0(pathway, " - ", tissue))
      } else {
        netVisual_aggregate(cellchat_obj, 
                            signaling = pathway, 
                            layout = layout,
                            edge.weight.max = weight.max[1], 
                            edge.width.max = 10, 
                            signaling.name = paste0(pathway, " - ", tissue))
      }
    }
    
    dev.off()
  }
  
  # 热图比较
  pdf(file.path(pathway_dir_specific, paste0(pathway, "_heatmap.pdf")), width = 12, height = 10)
  ht_list <- list()
  
  for(i in seq_along(cellchat_list)) {
    tissue <- names(cellchat_list)[i]
    cellchat_obj <- cellchat_list[[i]]
    
    # 检查通路是否存在
    if(!pathway %in% names(cellchat_obj@netP$pathways)) {
      next
    }
    
    ht_list[[i]] <- netVisual_heatmap(cellchat_obj, 
                                      signaling = pathway, 
                                      color.heatmap = "Reds",
                                      title.name = paste0(pathway, " - ", tissue))
  }
  
  if(length(ht_list) > 1) {
    combined_ht <- do.call(`+`, ht_list)
    ComplexHeatmap::draw(combined_ht, ht_gap = unit(0.5, "cm"))
  } else if(length(ht_list) == 1) {
    ComplexHeatmap::draw(ht_list[[1]])
  }
  dev.off()
  
  # 配体-受体对贡献分析
  pdf(file.path(pathway_dir_specific, paste0(pathway, "_contribution.pdf")), width = 12, height = 10)
  
  for(i in seq_along(cellchat_list)) {
    tissue <- names(cellchat_list)[i]
    cellchat_obj <- cellchat_list[[i]]
    
    # 检查通路是否存在
    if(!pathway %in% names(cellchat_obj@netP$pathways)) {
      next
    }
    
    # 分析配体-受体对对信号通路的贡献
    gg <- netAnalysis_contribution(cellchat_obj, 
                                   signaling = pathway,
                                   title = paste0(pathway, " - ", tissue))
    print(gg)
    
    # 提取主要的配体-受体对进行可视化
    pairLR <- extractEnrichedLR(cellchat_obj, 
                                signaling = pathway, 
                                geneLR.return = FALSE)
    
    # 如果存在配体-受体对，可视化前几个
    if(nrow(pairLR) > 0) {
      for(j in 1:min(3, nrow(pairLR))) {
        LR.show <- pairLR[j, ]
        
        # 使用圆形布局可视化
        netVisual_individual(cellchat_obj, 
                             signaling = pathway, 
                             pairLR.use = LR.show, 
                             layout = "circle",
                             title.name = paste0(pathway, " - ", LR.show, " - ", tissue))
      }
    }
  }
  
  dev.off()
  
  # 表达情况分析
  pdf(file.path(pathway_dir_specific, paste0(pathway, "_expression.pdf")), width = 12, height = 10)
  
  # 设置数据集顺序
  merged_cellchat@meta$datasets <- factor(merged_cellchat@meta$datasets, 
                                          levels = names(cellchat_list))
  
  # 绘制表达分布图
  tryCatch({
    p1 <- plotGeneExpression(merged_cellchat, 
                             signaling = pathway, 
                             split.by = "datasets", 
                             colors.ggplot = TRUE)
    print(p1)
    
    p2 <- plotGeneExpression(merged_cellchat, 
                             signaling = pathway, 
                             split.by = "datasets", 
                             colors.ggplot = TRUE, 
                             type = "dot")
    print(p2)
  }, error = function(e) {
    message(paste0("绘制表达分布图时出错: ", e$message))
  })
  
  dev.off()
}

#-------------------------
# 第十步：比较细胞类型特异通讯
#-------------------------
message("比较细胞类型特异通讯...")

# 获取所有组织中共有的细胞类型
common_celltypes <- Reduce(intersect, lapply(cellchat_list, function(x) levels(x@idents)))

if(length(common_celltypes) > 0) {
  message(paste0("所有组织共有", length(common_celltypes), "个细胞类型: ", 
                 paste(common_celltypes, collapse = ", ")))
  
  # 创建细胞类型比较目录
  celltype_dir <- file.path(output_dir, "CellType_Comparison")
  dir.create(celltype_dir, showWarnings = TRUE, recursive = TRUE)
  
  # 分析每个共有细胞类型
  for(cell_type in common_celltypes) {
    message(paste0("分析细胞类型: ", cell_type))
    
    # 创建细胞类型目录
    cell_dir <- file.path(celltype_dir, paste0("CellType_", cell_type))
    dir.create(cell_dir, showWarnings = TRUE, recursive = TRUE)
    
    # 创建PDF文件
    pdf(file.path(cell_dir, paste0(cell_type, "_across_tissues.pdf")), width = 12, height = 10)
    
    # 准备数据容器
    tissue_stats <- data.frame(
      tissue = character(),
      outgoing_count = numeric(),
      incoming_count = numeric(),
      outgoing_weight = numeric(),
      incoming_weight = numeric(),
      top_outgoing_pathways = character(),
      top_incoming_pathways = character(),
      top_outgoing_targets = character(),
      top_incoming_sources = character(),
      stringsAsFactors = FALSE
    )
    
    # 获取包含该细胞类型的组织
    tissues_with_celltype <- names(cellchat_list)[
      sapply(cellchat_list, function(x) {
        cell_type %in% levels(x@idents)
      })
    ]
    
    # 如果只有一个组织，跳过比较
    if(length(tissues_with_celltype) <= 1) {
      message(paste0("细胞类型", cell_type, "在多个组织中未找到，无法比较"))
      next
    }
    
    # 1. 比较该细胞类型作为信号源的通讯模式
    par(mfrow = c(1, length(tissues_with_celltype)), xpd=TRUE)
    
    for(i in seq_along(tissues_with_celltype)) {
      tissue <- tissues_with_celltype[i]
      cellchat_obj <- cellchat_list[[tissue]]
      mat <- cellchat_obj@net$weight
      
      # 提取该细胞类型发出的信号
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[cell_type, ] <- mat[cell_type, ]
      
      netVisual_circle(mat2, 
                       vertex.weight = as.numeric(table(cellchat_obj@idents)), 
                       weight.scale = TRUE, 
                       edge.weight.max = max(mat), 
                       title.name = paste0(cell_type, " 在 ", tissue, " 中发出的信号"))
      
      # 计算发送统计
      out_count <- sum(mat[cell_type, ] > 0)
      out_weight <- sum(mat[cell_type, ])
      
      # 保存统计数据
      tissue_stats <- rbind(tissue_stats, 
                            data.frame(tissue = tissue,
                                       outgoing_count = out_count,
                                       incoming_count = NA,
                                       outgoing_weight = out_weight,
                                       incoming_weight = NA,
                                       top_outgoing_pathways = NA,
                                       top_incoming_pathways = NA,
                                       top_outgoing_targets = NA,
                                       top_incoming_sources = NA,
                                       stringsAsFactors = FALSE))
    }
    
    # 2. 比较该细胞类型作为信号接收者的通讯模式
    par(mfrow = c(1, length(tissues_with_celltype)), xpd=TRUE)
    
    for(i in seq_along(tissues_with_celltype)) {
      tissue <- tissues_with_celltype[i]
      cellchat_obj <- cellchat_list[[tissue]]
      mat <- cellchat_obj@net$weight
      
      # 提取接收到的信号
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[, cell_type] <- mat[, cell_type]
      
      netVisual_circle(mat2, 
                       vertex.weight = as.numeric(table(cellchat_obj@idents)), 
                       weight.scale = TRUE, 
                       edge.weight.max = max(mat), 
                       title.name = paste0(cell_type, " 在 ", tissue, " 中接收的信号"))
      
      # 计算接收统计
      in_count <- sum(mat[, cell_type] > 0)
      in_weight <- sum(mat[, cell_type])
      
      # 更新统计数据
      row_idx <- which(tissue_stats$tissue == tissue)
      if(length(row_idx) > 0) {
        tissue_stats$incoming_count[row_idx] <- in_count
        tissue_stats$incoming_weight[row_idx] <- in_weight
      }
    }
    
    # 3. 比较该细胞类型的信号通路强度
    pathway_list <- list()
    
    for(i in seq_along(tissues_with_celltype)) {
      tissue <- tissues_with_celltype[i]
      cellchat_obj <- cellchat_list[[tissue]]
      
      # 检查是否存在通路数据
      if(length(cellchat_obj@netP$pathways) == 0) {
        next
      }
      
      # 提取该细胞类型在各通路中的强度
      pathways_strength <- sapply(cellchat_obj@netP$pathways, function(p) {
        # 以该细胞类型作为发送者
        sum_outgoing <- sum(p[cell_type, ])
        # 以该细胞类型作为接收者
        sum_incoming <- sum(p[, cell_type])
        return(c(outgoing = sum_outgoing, incoming = sum_incoming))
      })
      
      pathway_list[[tissue]] <- pathways_strength
      
      # 找出顶级通路
      # 顶级发送通路
      outgoing_paths <- pathways_strength["outgoing", ]
      outgoing_paths <- outgoing_paths[outgoing_paths > 0]
      if(length(outgoing_paths) > 0) {
        top_out_paths <- names(sort(outgoing_paths, decreasing = TRUE)[1:min(3, length(outgoing_paths))])
        tissue_stats$top_outgoing_pathways[tissue_stats$tissue == tissue] <- paste(top_out_paths, collapse = ", ")
      }
      
      # 顶级接收通路
      incoming_paths <- pathways_strength["incoming", ]
      incoming_paths <- incoming_paths[incoming_paths > 0]
      if(length(incoming_paths) > 0) {
        top_in_paths <- names(sort(incoming_paths, decreasing = TRUE)[1:min(3, length(incoming_paths))])
        tissue_stats$top_incoming_pathways[tissue_stats$tissue == tissue] <- paste(top_in_paths, collapse = ", ")
      }
      
      # 找出主要互作细胞类型
      if(nrow(mat) > 1) {
        # 主要目标细胞
        out_weights <- mat[cell_type, ]
        out_weights <- out_weights[out_weights > 0]
        if(length(out_weights) > 0) {
          top_targets <- names(sort(out_weights, decreasing = TRUE)[1:min(3, length(out_weights))])
          tissue_stats$top_outgoing_targets[tissue_stats$tissue == tissue] <- paste(top_targets, collapse = ", ")
        }
        
        # 主要来源细胞
        in_weights <- mat[, cell_type]
        in_weights <- in_weights[in_weights > 0]
        if(length(in_weights) > 0) {
          top_sources <- names(sort(in_weights, decreasing = TRUE)[1:min(3, length(in_weights))])
          tissue_stats$top_incoming_sources[tissue_stats$tissue == tissue] <- paste(top_sources, collapse = ", ")
        }
      }
    }
    
    # 整合数据并可视化
    # 转换为长格式数据用于绘图
    if(length(pathway_list) > 0) {
      pathway_data <- do.call(rbind, lapply(names(pathway_list), function(tissue) {
        df <- as.data.frame(t(pathway_list[[tissue]]))
        df$pathway <- rownames(df)
        df$tissue <- tissue
        return(df)
      }))
      
      if(nrow(pathway_data) > 0) {
        # 创建长格式数据
        pathway_long <- reshape2::melt(pathway_data, 
                                       id.vars = c("pathway", "tissue"),
                                       variable.name = "direction",
                                       value.name = "strength")
        
        # 找出最强的前10个通路
        top_pathways <- pathway_long %>%
          dplyr::group_by(pathway) %>%
          dplyr::summarise(total_strength = sum(strength, na.rm = TRUE)) %>%
          dplyr::arrange(desc(total_strength)) %>%
          dplyr::slice_head(n = 10) %>%
          dplyr::pull(pathway)
        
        # 过滤数据以仅包含顶级通路
        plot_data <- pathway_long %>%
          dplyr::filter(pathway %in% top_pathways)
        
        # 绘制比较图
        p <- ggplot(plot_data, aes(x = pathway, y = strength, fill = tissue)) +
          geom_bar(stat = "identity", position = "dodge") +
          facet_wrap(~ direction) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = paste0(cell_type, " 在不同组织中的信号通路强度"),
               x = "通路", y = "通讯强度")
        
        print(p)
      }
    }
    
    # 4. 气泡图比较不同组织中该细胞类型与其他细胞类型的互作
    for(i in seq_along(tissues_with_celltype)) {
      tissue <- tissues_with_celltype[i]
      cellchat_obj <- cellchat_list[[tissue]]
      
      # 该细胞类型作为信号源
      p1 <- netVisual_bubble(cellchat_obj, 
                             sources.use = cell_type, 
                             remove.isolate = TRUE, 
                             title.name = paste0(cell_type, " 在 ", tissue, " 中作为信号源"))
      
      # 该细胞类型作为信号目标
      p2 <- netVisual_bubble(cellchat_obj, 
                             targets.use = cell_type, 
                             remove.isolate = TRUE, 
                             title.name = paste0(cell_type, " 在 ", tissue, " 中作为信号靶点"))
      
      print(p1)
      print(p2)
    }
    
    dev.off()
    
    # 保存统计数据
    write.csv(tissue_stats, 
              file = file.path(cell_dir, paste0(cell_type, "_stats.csv")),
              row.names = FALSE)
    
    if(exists("pathway_data") && nrow(pathway_data) > 0) {
      write.csv(pathway_data, 
                file = file.path(cell_dir, paste0(cell_type, "_pathway_strength.csv")),
                row.names = FALSE)
    }
    
    message(paste0("完成细胞类型 ", cell_type, " 的通讯模式分析"))
  }
  
  # 创建细胞类型比较摘要报告
  if(length(common_celltypes) > 1) {
    message("创建细胞类型比较摘要报告...")
    
    # 准备综合数据
    all_stats <- data.frame()
    
    # 读取所有细胞类型的统计数据并合并
    for(cell_type in common_celltypes) {
      stats_file <- file.path(celltype_dir, paste0("CellType_", cell_type), 
                              paste0(cell_type, "_stats.csv"))
      
      if(file.exists(stats_file)) {
        stats <- read.csv(stats_file)
        stats$cell_type <- cell_type
        all_stats <- rbind(all_stats, stats)
      }
    }
    
    if(nrow(all_stats) > 0) {
      # 计算额外指标
      all_stats$total_weight <- all_stats$outgoing_weight + all_stats$incoming_weight
      all_stats$send_receive_ratio <- all_stats$outgoing_weight / (all_stats$incoming_weight + 0.1)
      
      # 创建PDF报告
      pdf(file.path(celltype_dir, "CellType_Comparison_Report.pdf"), width = 12, height = 10)
      
      # 各细胞类型在不同组织中的总通讯强度
      # 创建热图数据
      heatmap_data <- reshape2::dcast(all_stats, cell_type ~ tissue, value.var = "total_weight")
      rownames(heatmap_data) <- heatmap_data$cell_type
      heatmap_data$cell_type <- NULL
      
      # 标准化数据以便比较
      heatmap_data_norm <- t(apply(as.matrix(heatmap_data), 1, function(x) {
        if(max(x, na.rm = TRUE) > 0) {
          return(x/max(x, na.rm = TRUE))
        } else {
          return(x)
        }
      }))
      
      # 绘制热图
      if(requireNamespace("pheatmap", quietly = TRUE) && nrow(heatmap_data) > 1) {
        pheatmap::pheatmap(heatmap_data_norm,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           display_numbers = TRUE,
                           fontsize_number = 8,
                           main = "细胞类型在各组织中的相对通讯活跃度")
      }
      
      # 各细胞类型的发送/接收倾向(比率)
      ratio_plot <- ggplot(all_stats, aes(x = cell_type, y = send_receive_ratio, fill = tissue)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_log10() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "细胞类型发送/接收信号比率(对数刻度)",
             subtitle = "值>1表示主要为发送者，值<1表示主要为接收者",
             x = "细胞类型", y = "发送/接收比率(对数)")
      
      print(ratio_plot)
      
      # 细胞类型的通讯多样性(互作细胞类型数量)
      diversity_plot <- ggplot(all_stats, aes(x = cell_type, y = outgoing_count + incoming_count, fill = tissue)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "细胞类型的通讯多样性",
             x = "细胞类型", y = "总互作数量")
      
      print(diversity_plot)
      
      # 创建综合统计表
      summary_table <- all_stats %>%
        dplyr::group_by(cell_type) %>%
        dplyr::summarise(
          avg_outgoing = mean(outgoing_weight, na.rm = TRUE),
          avg_incoming = mean(incoming_weight, na.rm = TRUE),
          avg_ratio = mean(send_receive_ratio, na.rm = TRUE),
          var_across_tissues = sd(total_weight, na.rm = TRUE) / mean(total_weight, na.rm = TRUE)
        ) %>%
        dplyr::arrange(desc(avg_outgoing + avg_incoming))
      
      if(requireNamespace("gridExtra", quietly = TRUE)) {
        gt <- gridExtra::tableGrob(summary_table, rows = NULL)
        grid::grid.newpage()
        grid::grid.draw(gt)
      }
      
      dev.off()
      
      # 保存综合数据
      write.csv(all_stats, 
                file = file.path(celltype_dir, "all_celltype_stats.csv"),
                row.names = FALSE)
      
      write.csv(summary_table,
                file = file.path(celltype_dir, "celltype_summary_table.csv"),
                row.names = FALSE)
    }
  }
} else {
  message("未找到所有组织共有的细胞类型")
}

#-------------------------
# 第十一步：生成综合报告
#-------------------------
message("生成综合报告...")

# 创建输出文件
report_file <- file.path(output_dir, "Tissue_Communication_Summary_Report.pdf")
pdf(report_file, width = 14, height = 10)

# 1. 总体通讯模式比较
gg1 <- compareInteractions(merged_cellchat, show.legend = FALSE, measure = "count")
gg2 <- compareInteractions(merged_cellchat, show.legend = FALSE, measure = "weight")
print(gg1 + gg2)

# 2. 信号通路强度比较
gg3 <- rankNet(merged_cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg4 <- rankNet(merged_cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
print(gg3 + gg4)

# 3. 显示各组织的通讯网络
for(tissue in all_tissues) {
  if(!tissue %in% names(cellchat_list)) next
  
  cellchat_obj <- cellchat_list[[tissue]]
  groupSize <- as.numeric(table(cellchat_obj@idents))
  
  # 数量网络
  netVisual_circle(cellchat_obj@net$count, 
                   vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   label.edge = FALSE,
                   title.name = paste0(tissue, " - 互作数量"))
  
  # 强度网络
  netVisual_circle(cellchat_obj@net$weight, 
                   vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   label.edge = FALSE,
                   title.name = paste0(tissue, " - 互作强度"))
}

# 4. 汇总主要发现
for(tissue in all_tissues) {
  if(!tissue %in% names(cellchat_list)) next
  
  cellchat_obj <- cellchat_list[[tissue]]
  if(nrow(cellchat_obj@net$weight) > 0) {
    # 最活跃的细胞类型
    outgoing <- rowSums(cellchat_obj@net$weight)
    incoming <- colSums(cellchat_obj@net$weight)
    
    active_senders <- sort(outgoing, decreasing = TRUE)
    active_receivers <- sort(incoming, decreasing = TRUE)
    
    # 主要信号通路
    if(length(cellchat_obj@netP$pathways) > 0) {
      pathway_strength <- sapply(cellchat_obj@netP$pathways, function(x) sum(x))
      top_pathways <- sort(pathway_strength, decreasing = TRUE)
      
      # 创建总结表格
      summary_data <- data.frame(
        Tissue = tissue,
        Active_Senders = paste(names(active_senders)[1:min(3, length(active_senders))], collapse = ", "),
        Active_Receivers = paste(names(active_receivers)[1:min(3, length(active_receivers))], collapse = ", "),
        Top_Pathways = paste(names(top_pathways)[1:min(5, length(top_pathways))], collapse = ", ")
      )
      
      # 打印表格
      grid.newpage()
      grid.table(summary_data)
    }
  }
}

# 5. 对比热图
# 合并所有通路
pathway.union <- unique(unlist(lapply(cellchat_list, function(x) names(x@netP$pathways))))

# 对outgoing和incoming模式绘制热图
for(pattern in c("outgoing", "incoming")) {
  color <- ifelse(pattern == "outgoing", "Reds", "GnBu")
  
  for(i in seq_along(cellchat_list)) {
    tissue <- names(cellchat_list)[i]
    cellchat_obj <- cellchat_list[[i]]
    
    if(length(cellchat_obj@netP$pathways) > 0) {
      ht <- netAnalysis_signalingRole_heatmap(cellchat_obj, 
                                              pattern = pattern, 
                                              signaling = pathway.union, 
                                              title = tissue, 
                                              width = 8, 
                                              height = 10, 
                                              color.heatmap = color)
      
      if(class(ht)[1] == "Heatmap") {
        ComplexHeatmap::draw(ht)
      }
    }
  }
}

dev.off()

#-------------------------
# 第十二步：输出总结信息
#-------------------------
message("创建分析总结文件...")

# 创建总结文件
summary_file <- file.path(output_dir, "analysis_summary.txt")
sink(summary_file)

cat("===============================================\n")
cat("组织间细胞通讯分析总结\n")
cat("===============================================\n\n")

cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("分析的组织:\n")
for(tissue in all_tissues) {
  if(!tissue %in% names(cellchat_list)) next
  
  # 修正Cell_Count的计算
  cell_count <- nrow(cellchat_list[[tissue]]@meta)
  if(cell_count == 0) {
    cell_count <- length(cellchat_list[[tissue]]@idents)
  }
  
  int_count <- sum(cellchat_list[[tissue]]@net$count != 0)
  unique_pathways <- length(cellchat_list[[tissue]]@netP$pathways)
  
  cat(paste0("- ", tissue, ": ", cell_count, " 细胞, ", 
             int_count, " 互作, ",
             unique_pathways, " 信号通路\n"))
}
cat("\n")

cat("信号通路分析:\n")
if(length(common_pathways) > 0) {
  cat(paste0("- 共有通路数量: ", length(common_pathways), "\n"))
  cat(paste0("- 前5个共有通路: ", paste(common_pathways[1:min(5, length(common_pathways))], collapse = ", "), "\n"))
} else {
  cat("- 未找到所有组织共有的通路\n")
}

cat(paste0("- 选择分析的通路: ", paste(selected_pathways, collapse = ", "), "\n\n"))

cat("细胞类型分析:\n")
if(length(common_celltypes) > 0) {
  cat(paste0("- 共有细胞类型数量: ", length(common_celltypes), "\n"))
  cat(paste0("- 共有细胞类型: ", paste(common_celltypes, collapse = ", "), "\n"))
} else {
  cat("- 未找到所有组织共有的细胞类型\n")
}
cat("\n")

cat("组织间主要差异:\n")
# 使用pathway_strength比较结果提取差异最大的通路
if(exists("pathway_strength") && ncol(pathway_strength) > 0) {
  # 计算各通路在组织间的差异
  pathway_var <- apply(pathway_strength, 2, function(x) {
    max_val <- max(x, na.rm = TRUE)
    min_val <- min(x, na.rm = TRUE)
    if(is.finite(max_val) && is.finite(min_val)) {
      return(max_val - min_val)
    } else {
      return(NA)
    }
  })
  
  # 过滤掉NA值
  pathway_var <- pathway_var[!is.na(pathway_var)]
  
  if(length(pathway_var) > 0) {
    # 排序并找出差异最大的通路
    pathway_diff <- sort(pathway_var, decreasing = TRUE)
    
    cat("强度差异最大的前5个通路:\n")
    for(i in 1:min(5, length(pathway_diff))) {
      pathway <- names(pathway_diff)[i]
      cat(paste0("- ", pathway, ": 差异值 = ", round(pathway_diff[i], 2), "\n"))
    }
  }
} else if(exists("gg4") && is.ggplot(gg4)) {
  # 尝试从ggplot对象中提取数据
  if(!is.null(gg4$data)) {
    gg4_data <- gg4$data
    if(nrow(gg4_data) > 0) {
      # 计算各通路在组织间的差异
      pathway_diff <- gg4_data %>%
        dplyr::group_by(features.plot) %>%
        dplyr::summarise(max_diff = max(interaction.strength) - min(interaction.strength)) %>%
        dplyr::arrange(desc(max_diff))
      
      cat("强度差异最大的前5个通路:\n")
      for(i in 1:min(5, nrow(pathway_diff))) {
        pathway <- pathway_diff$features.plot[i]
        cat(paste0("- ", pathway, ": 差异值 = ", round(pathway_diff$max_diff[i], 2), "\n"))
      }
    }
  }
}

# 添加组织特异的总结
cat("\n各组织特征总结:\n")
for(tissue in all_tissues) {
  if(!tissue %in% names(cellchat_list)) next
  
  cellchat_obj <- cellchat_list[[tissue]]
  cat(paste0("\n", tissue, ":\n"))
  
  # 最活跃的细胞类型
  if(nrow(cellchat_obj@net$weight) > 0) {
    outgoing <- rowSums(cellchat_obj@net$weight)
    incoming <- colSums(cellchat_obj@net$weight)
    
    # 排序并找出前3个
    active_senders <- sort(outgoing, decreasing = TRUE)
    active_receivers <- sort(incoming, decreasing = TRUE)
    
    cat("- 主要发送者: ", paste(names(active_senders)[1:min(3, length(active_senders))], collapse = ", "), "\n")
    cat("- 主要接收者: ", paste(names(active_receivers)[1:min(3, length(active_receivers))], collapse = ", "), "\n")
    
    # 主要信号通路
    if(length(cellchat_obj@netP$pathways) > 0) {
      pathway_strength <- sapply(cellchat_obj@netP$pathways, function(x) sum(x))
      top_pathways <- sort(pathway_strength, decreasing = TRUE)
      
      cat("- 主要通路: ", paste(names(top_pathways)[1:min(5, length(top_pathways))], collapse = ", "), "\n")
    }
    
    # 中介细胞类型（网络中心性高的细胞）
    if("centrality" %in% names(cellchat_obj@netP)) {
      if(length(cellchat_obj@netP$centrality) > 0) {
        if("information" %in% names(cellchat_obj@netP$centrality)) {
          info_centrality <- cellchat_obj@netP$centrality$information
          if(length(info_centrality) > 0) {
            top_mediators <- sort(info_centrality, decreasing = TRUE)
            cat("- 关键中介细胞: ", paste(names(top_mediators)[1:min(3, length(top_mediators))], collapse = ", "), "\n")
          }
        }
      }
    }
  }
}

cat("\n===============================================\n")
cat("分析完成!\n")
cat("结果保存在: ", output_dir, "\n")
cat("===============================================\n")

sink()

#-------------------------
# 第十三步：创建可视化摘要
#-------------------------
message("创建可视化摘要...")

# 创建可视化摘要目录
vis_dir <- file.path(output_dir, "Visualization_Summary")
dir.create(vis_dir, showWarnings = TRUE, recursive = TRUE)

# 1. 创建互作数量和强度的比较图
if(exists("stats_data") && nrow(stats_data) > 0) {
  # 修正Cell_Count为正确值
  for(i in 1:nrow(stats_data)) {
    tissue <- stats_data$Tissue[i]
    if(tissue %in% names(cellchat_list)) {
      stats_data$Cell_Count[i] <- length(cellchat_list[[tissue]]@idents)
    }
  }
  
  # 重新计算每千细胞互作数
  stats_data$Int_Per_1k_Cells <- stats_data$Int_Count / (stats_data$Cell_Count/1000)
  
  # 创建比较柱状图
  interaction_plot <- ggplot(stats_data, aes(x = Tissue)) +
    geom_bar(aes(y = Int_Count, fill = "互作数量"), stat = "identity", position = "dodge", alpha = 0.7) +
    geom_bar(aes(y = Int_Per_1k_Cells * 10, fill = "每千细胞互作(×10)"), stat = "identity", position = "dodge", alpha = 0.7) +
    scale_fill_manual(values = c("互作数量" = "darkblue", "每千细胞互作(×10)" = "orangered")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
    labs(title = "组织间互作数量比较",
         x = "组织", y = "互作数量")
  
  # 保存图形
  ggsave(file.path(vis_dir, "interaction_counts_comparison.pdf"), 
         interaction_plot, width = 10, height = 6)
  
  # 创建通路数量比较图
  pathway_plot <- ggplot(stats_data, aes(x = Tissue, y = Unique_Pathways, fill = Tissue)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "组织间激活通路数量比较",
         x = "组织", y = "激活通路数量")
  
  # 保存图形
  ggsave(file.path(vis_dir, "pathway_counts_comparison.pdf"), 
         pathway_plot, width = 10, height = 6)
}

# 2. 创建通路强度热图
if(exists("heatmap_data") && ncol(heatmap_data) > 0) {
  # 保存通路热图
  pdf(file.path(vis_dir, "pathway_heatmap_summary.pdf"), width = 12, height = 10)
  
  # 绘制热图
  pheatmap::pheatmap(heatmap_data,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     display_numbers = TRUE,
                     fontsize_number = 8,
                     main = "组织间信号通路强度热图")
  
  dev.off()
}

# 3. 创建细胞类型通讯模式摘要
if(exists("all_stats") && nrow(all_stats) > 0) {
  # 计算各细胞类型的发送/接收总强度
  celltype_summary <- all_stats %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      avg_outgoing = mean(outgoing_weight, na.rm = TRUE),
      avg_incoming = mean(incoming_weight, na.rm = TRUE),
      total_strength = sum(outgoing_weight + incoming_weight, na.rm = TRUE)
    ) %>%
    dplyr::arrange(desc(total_strength))
  
  # 选择前10个最活跃的细胞类型
  top_celltypes <- celltype_summary$cell_type[1:min(10, nrow(celltype_summary))]
  
  # 筛选数据
  plot_data <- all_stats %>%
    dplyr::filter(cell_type %in% top_celltypes)
  
  # 创建通讯强度柱状图
  celltype_plot <- ggplot(plot_data, aes(x = cell_type)) +
    geom_bar(aes(y = outgoing_weight, fill = "发送强度"), stat = "identity", position = "stack") +
    geom_bar(aes(y = -incoming_weight, fill = "接收强度"), stat = "identity", position = "stack") +
    scale_fill_manual(values = c("发送强度" = "firebrick", "接收强度" = "steelblue")) +
    facet_wrap(~tissue) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
    labs(title = "主要细胞类型通讯强度比较",
         x = "细胞类型", y = "通讯强度 (正值:发送, 负值:接收)")
  
  # 保存图形
  ggsave(file.path(vis_dir, "celltype_communication_summary.pdf"), 
         celltype_plot, width = 12, height = 8)
}

# 4. 创建组织差异通路网络图摘要
# 选择一张差异最显著的组织对比图
gg_diff_summary <- netVisual_diffInteraction(merged_cellchat, 
                                             weight.scale = TRUE, 
                                             measure = "weight",
                                             top = 0.25)  # 只显示权重最高的25%

# 保存图形
pdf(file.path(vis_dir, "diff_interaction_summary.pdf"), width = 12, height = 10)
gg_diff_summary
dev.off()

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
getwd()

# =============================================================================
# 整合版STARTRAC-dist组织分布分析系统 (Integrated STARTRAC-dist Analysis System)
# 融合稳健性改进、完整统计方法和高质量可视化
# 作者：Claude & 用户协作开发
# =============================================================================

# 加载必需的包
required_packages <- c(
  "Startrac", "Seurat", "ggplot2", "ComplexHeatmap", "RColorBrewer", 
  "circlize", "tidyverse", "ggpubr", "pheatmap", "corrplot", 
  "cowplot", "viridis", "scales", "readr", "qs", 
  "gridExtra", "reshape2", "VennDiagram", "grid", "yaml"
)

# 智能包管理函数
setup_packages <- function(packages) {
  missing_packages <- packages[!sapply(packages, require, character.only = TRUE, quietly = TRUE)]
  
  if(length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    
    # 特殊处理某些包
    for(pkg in missing_packages) {
      tryCatch({
        if(pkg == "Startrac") {
          # STARTRAC可能需要从GitHub安装
          if(!require("devtools", quietly = TRUE)) install.packages("devtools")
          devtools::install_github("Japrin/STARTRAC")
        } else {
          install.packages(pkg)
        }
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        warning(paste("Failed to install", pkg, ":", e$message))
      })
    }
  }
  
  # 设置全局主题
  if(require("ggplot2", quietly = TRUE)) {
    theme_set(theme_bw() + theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ))
  }
}

setup_packages(required_packages)

# =============================================================================
# 核心数据处理与验证函数 (Core Data Processing & Validation)
# =============================================================================

#' 稳健的Ro/e值转换为符号表示 (按STARTRAC标准)
#' @param roe_value Ro/e指数值
#' @param thresholds 阈值向量，默认按STARTRAC标准c(0.5, 0.8, 1.2, 2.0)
#' @return 符号字符串
roe_to_symbol <- function(roe_value, thresholds = c(0.5, 0.8, 1.2, 2.0)) {
  # 处理各种异常值
  if(is.null(roe_value) || length(roe_value) == 0 || is.na(roe_value) || is.nan(roe_value)) {
    return("-")
  }
  
  if(is.infinite(roe_value)) {
    return(ifelse(roe_value > 0, "+++", "---"))
  }
  
  if(!is.numeric(roe_value)) return("-")
  
  # 按照STARTRAC标准分类
  if(roe_value >= thresholds[4]) {        # >= 2.0
    return("+++")
  } else if(roe_value >= thresholds[3]) { # >= 1.2
    return("+")
  } else if(roe_value > thresholds[2]) {  # > 0.8
    return("+/-")
  } else if(roe_value >= thresholds[1]) { # >= 0.5
    return("-")
  } else {                                # < 0.5
    return("---")
  }
}

#' Ro/e值的生物学解释
#' @param roe_value Ro/e指数值
get_roe_interpretation <- function(roe_value) {
  if(is.na(roe_value) || !is.finite(roe_value)) {
    return("Data not available")
  } else if(roe_value >= 2.0) {
    return("Strongly enriched")
  } else if(roe_value >= 1.5) {
    return("Moderately enriched")
  } else if(roe_value >= 1.2) {
    return("Mildly enriched")
  } else if(roe_value >= 0.8) {
    return("Near expected levels")
  } else if(roe_value >= 0.5) {
    return("Mildly depleted")
  } else if(roe_value >= 0.3) {
    return("Moderately depleted")
  } else {
    return("Strongly depleted")
  }
}

#' 数据预处理和验证
#' @param roe_matrix Ro/e指数矩阵
#' @return 清洁的Ro/e矩阵
clean_roe_matrix <- function(roe_matrix) {
  message("Cleaning and validating Ro/e matrix...")
  
  # 检查输入
  if(is.null(roe_matrix) || (!is.matrix(roe_matrix) && !is.data.frame(roe_matrix))) {
    stop("Invalid Ro/e matrix provided")
  }
  
  # 转换为矩阵
  if(is.data.frame(roe_matrix)) {
    roe_matrix <- as.matrix(roe_matrix)
  }
  
  if(nrow(roe_matrix) == 0 || ncol(roe_matrix) == 0) {
    stop("Empty Ro/e matrix")
  }
  
  # 替换无效值
  roe_matrix[is.nan(roe_matrix)] <- NA
  roe_matrix[is.infinite(roe_matrix)] <- NA
  
  # 检查数据范围合理性
  finite_values <- roe_matrix[is.finite(roe_matrix)]
  
  message(sprintf("  Matrix dimensions: %d cell types × %d tissues", 
                  nrow(roe_matrix), ncol(roe_matrix)))
  message(sprintf("  Valid values: %d/%d (%.1f%%)", 
                  sum(is.finite(roe_matrix)), length(roe_matrix),
                  100 * sum(is.finite(roe_matrix)) / length(roe_matrix)))
  
  if(length(finite_values) > 0) {
    message(sprintf("  Value range: %.3f - %.3f", 
                    min(finite_values, na.rm = TRUE), 
                    max(finite_values, na.rm = TRUE)))
  } else {
    warning("No finite values in Ro/e matrix")
  }
  
  return(roe_matrix)
}

#' 显示可用的metadata列信息
#' @param seurat_obj Seurat对象
#' @param show_preview 是否显示每列的预览
show_metadata_info <- function(seurat_obj, show_preview = TRUE) {
  
  meta_data <- seurat_obj@meta.data
  available_cols <- colnames(meta_data)
  
  message("\n=== AVAILABLE METADATA COLUMNS ===")
  message("Total columns: ", length(available_cols))
  
  for(i in seq_along(available_cols)) {
    col_name <- available_cols[i]
    col_data <- meta_data[[col_name]]
    
    # 基本信息
    unique_count <- length(unique(col_data[!is.na(col_data)]))
    na_count <- sum(is.na(col_data))
    completion_rate <- round((1 - na_count/length(col_data)) * 100, 1)
    
    message(sprintf("  %2d: %-30s [%d unique, %.1f%% complete]", 
                    i, col_name, unique_count, completion_rate))
    
    # 显示预览
    if(show_preview) {
      if(is.numeric(col_data)) {
        range_info <- paste0("Range: ", round(min(col_data, na.rm = TRUE), 2), 
                             " - ", round(max(col_data, na.rm = TRUE), 2))
        message(sprintf("      %s", range_info))
      } else {
        unique_vals <- unique(col_data[!is.na(col_data)])
        preview_vals <- if(length(unique_vals) <= 5) {
          paste(unique_vals, collapse = ", ")
        } else {
          paste(c(head(unique_vals, 3), "..."), collapse = ", ")
        }
        message(sprintf("      Values: %s", preview_vals))
      }
    }
  }
  
  return(available_cols)
}

# =============================================================================
# 稳健的STARTRAC-dist可视化系统 (Robust STARTRAC-dist Visualization)
# =============================================================================

#' 创建稳健的STARTRAC-dist风格可视化 (主函数)
#' @param roe_matrix Ro/e指数矩阵
#' @param significance_results 显著性结果（可选）
#' @param output_path 输出文件路径
#' @param title 图像标题
#' @param method 绘图方法 ("grid", "base", "auto")
#' @param cell_type_groups 细胞类型分组（可选）
create_startrac_dist_plot <- function(roe_matrix, 
                                      significance_results = NULL, 
                                      output_path = "startrac_dist_plot.pdf",
                                      title = "STARTRAC-dist Analysis",
                                      method = "auto",
                                      cell_type_groups = NULL) {
  
  message("Creating STARTRAC-dist style visualization...")
  
  # 1. 数据预处理
  roe_matrix <- clean_roe_matrix(roe_matrix)
  
  # 基本参数
  tissues <- colnames(roe_matrix)
  celltypes <- rownames(roe_matrix)
  n_tissues <- length(tissues)
  n_celltypes <- length(celltypes)
  
  if(is.null(tissues) || is.null(celltypes)) {
    stop("Matrix must have row and column names")
  }
  
  # 2. 创建符号矩阵
  symbol_matrix <- matrix("", nrow = n_celltypes, ncol = n_tissues)
  rownames(symbol_matrix) <- celltypes
  colnames(symbol_matrix) <- tissues
  
  for(i in 1:n_celltypes) {
    for(j in 1:n_tissues) {
      symbol_matrix[i, j] <- roe_to_symbol(roe_matrix[i, j])
    }
  }
  
  # 3. 智能选择绘图方法
  if(method == "auto") {
    # 根据数据大小和复杂性选择方法
    if(n_tissues * n_celltypes > 500 || any(nchar(c(tissues, celltypes)) > 20)) {
      method <- "base"  # 大数据或长标签用基础绘图
    } else {
      method <- "grid"  # 默认使用grid系统
    }
  }
  
  # 4. 执行绘图
  if(method == "grid") {
    create_grid_plot(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes)
  } else {
    create_base_plot(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes)
  }
  
  message(paste("STARTRAC-dist plot saved to:", output_path))
}

#' Grid系统绘图方法 (精确布局)
create_grid_plot <- function(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes) {
  
  tryCatch({
    n_tissues <- length(tissues)
    n_celltypes <- length(celltypes)
    
    pdf(output_path, width = max(12, n_tissues * 1.5), height = max(8, n_celltypes * 0.4))
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2, widths = c(1, 1))))
    
    # === 左侧：符号矩阵图 ===
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    pushViewport(viewport(x = 0.15, y = 0.15, width = 0.7, height = 0.7, just = c("left", "bottom")))
    
    cell_width <- 1 / n_tissues
    cell_height <- 1 / n_celltypes
    
    # 颜色映射函数
    get_colors <- function(roe_val) {
      if(!is.finite(roe_val)) return(list(bg = "#F0F0F0", text = "black"))
      
      if(roe_val >= 2.0) return(list(bg = "#8B0000", text = "white"))
      if(roe_val >= 1.2) return(list(bg = "#DC143C", text = "white"))
      if(roe_val > 0.8) return(list(bg = "#FFFACD", text = "black"))
      if(roe_val >= 0.5) return(list(bg = "#87CEEB", text = "black"))
      return(list(bg = "#4682B4", text = "white"))
    }
    
    # 绘制符号网格
    for(i in 1:n_celltypes) {
      for(j in 1:n_tissues) {
        x <- (j - 0.5) / n_tissues
        y <- (n_celltypes - i + 0.5) / n_celltypes
        
        roe_val <- roe_matrix[i, j]
        colors <- get_colors(roe_val)
        symbol <- symbol_matrix[i, j]
        
        grid.rect(x = x, y = y, 
                  width = cell_width * 0.9, height = cell_height * 0.9,
                  gp = gpar(fill = colors$bg, col = "white", lwd = 1))
        grid.text(symbol, x = x, y = y, 
                  gp = gpar(col = colors$text, fontsize = 11, fontface = "bold"))
      }
    }
    
    # 添加标签
    for(i in 1:n_celltypes) {
      y <- (n_celltypes - i + 0.5) / n_celltypes
      grid.text(celltypes[i], x = -0.02, y = y, just = "right", gp = gpar(fontsize = 9))
    }
    
    for(j in 1:n_tissues) {
      x <- (j - 0.5) / n_tissues
      grid.text(tissues[j], x = x, y = -0.02, just = "center", rot = 45, gp = gpar(fontsize = 9))
    }
    
    grid.text("STARTRAC-dist Symbols", x = 0.5, y = 1.05, just = "center",
              gp = gpar(fontsize = 14, fontface = "bold"))
    
    popViewport(2)
    
    # === 右侧：数值热图 ===
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    pushViewport(viewport(x = 0.15, y = 0.15, width = 0.7, height = 0.7, just = c("left", "bottom")))
    
    # 创建颜色渐变
    finite_values <- roe_matrix[is.finite(roe_matrix)]
    if(length(finite_values) > 0) {
      val_min <- min(finite_values, na.rm = TRUE)
      val_max <- max(finite_values, na.rm = TRUE)
      color_func <- colorRamp2(c(val_min, 1, val_max), c("#2166AC", "white", "#B2182B"))
    } else {
      color_func <- colorRamp2(c(0, 1, 2), c("#2166AC", "white", "#B2182B"))
    }
    
    # 绘制数值热图
    for(i in 1:n_celltypes) {
      for(j in 1:n_tissues) {
        x <- (j - 0.5) / n_tissues
        y <- (n_celltypes - i + 0.5) / n_celltypes
        
        roe_val <- roe_matrix[i, j]
        
        if(is.finite(roe_val)) {
          bg_color <- color_func(roe_val)
          text_color <- ifelse(abs(roe_val - 1) > 0.5, "white", "black")
          
          grid.rect(x = x, y = y, 
                    width = cell_width * 0.9, height = cell_height * 0.9,
                    gp = gpar(fill = bg_color, col = "white", lwd = 1))
          grid.text(sprintf("%.2f", roe_val), x = x, y = y, 
                    gp = gpar(col = text_color, fontsize = 9))
        } else {
          grid.rect(x = x, y = y, 
                    width = cell_width * 0.9, height = cell_height * 0.9,
                    gp = gpar(fill = "#F0F0F0", col = "white", lwd = 1))
          grid.text("NA", x = x, y = y, gp = gpar(col = "gray50", fontsize = 8))
        }
      }
    }
    
    # 添加标签
    for(i in 1:n_celltypes) {
      y <- (n_celltypes - i + 0.5) / n_celltypes
      grid.text(celltypes[i], x = 1.02, y = y, just = "left", gp = gpar(fontsize = 9))
    }
    
    for(j in 1:n_tissues) {
      x <- (j - 0.5) / n_tissues
      grid.text(tissues[j], x = x, y = -0.02, just = "center", rot = 45, gp = gpar(fontsize = 9))
    }
    
    grid.text("Ro/e Index Values", x = 0.5, y = 1.05, just = "center",
              gp = gpar(fontsize = 14, fontface = "bold"))
    
    popViewport(2)
    popViewport()
    
    # === 添加总标题和图例 ===
    add_title_and_legend(title)
    
    dev.off()
    
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    message("Grid method failed: ", e$message, ". Trying base graphics...")
    create_base_plot(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes)
  })
}

#' 基础绘图系统方法 (兼容性备用)
create_base_plot <- function(roe_matrix, symbol_matrix, output_path, title, tissues, celltypes) {
  
  tryCatch({
    n_tissues <- length(tissues)
    n_celltypes <- length(celltypes)
    
    pdf(output_path, width = max(12, n_tissues * 1.2), height = max(8, n_celltypes * 0.35))
    
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
    
    # 左侧：符号矩阵
    par(mar = c(8, 8, 4, 2))
    
    plot(1, 1, type = "n", xlim = c(0.5, n_tissues + 0.5), 
         ylim = c(0.5, n_celltypes + 0.5),
         xlab = "", ylab = "", main = "STARTRAC Symbols", 
         axes = FALSE)
    
    for(i in 1:n_celltypes) {
      for(j in 1:n_tissues) {
        symbol <- symbol_matrix[i, j]
        roe_val <- roe_matrix[i, j]
        
        # 设置颜色
        if(is.finite(roe_val)) {
          if(roe_val >= 2.0) { bg_col <- "#8B0000"; text_col <- "white" }
          else if(roe_val >= 1.2) { bg_col <- "#DC143C"; text_col <- "white" }
          else if(roe_val > 0.8) { bg_col <- "#FFFACD"; text_col <- "black" }
          else if(roe_val >= 0.5) { bg_col <- "#87CEEB"; text_col <- "black" }
          else { bg_col <- "#4682B4"; text_col <- "white" }
        } else {
          bg_col <- "#F0F0F0"; text_col <- "black"
        }
        
        rect(j - 0.4, i - 0.4, j + 0.4, i + 0.4, 
             col = bg_col, border = "white", lwd = 1.5)
        text(j, i, symbol, col = text_col, cex = 1.1, font = 2)
      }
    }
    
    axis(1, at = 1:n_tissues, labels = tissues, las = 2, cex.axis = 0.8)
    axis(2, at = 1:n_celltypes, labels = celltypes, las = 2, cex.axis = 0.8)
    
    # 右侧：数值图
    par(mar = c(8, 2, 4, 8))
    
    finite_vals <- roe_matrix[is.finite(roe_matrix)]
    if(length(finite_vals) > 0) {
      val_range <- range(finite_vals, na.rm = TRUE)
      colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(50)
      
      image(1:n_tissues, 1:n_celltypes, t(roe_matrix), 
            col = colors, axes = FALSE, main = "Ro/e Values", zlim = val_range)
      
      # 添加数值标签
      for(i in 1:n_celltypes) {
        for(j in 1:n_tissues) {
          if(is.finite(roe_matrix[i, j])) {
            text_color <- ifelse(abs(roe_matrix[i, j] - 1) > 0.5, "white", "black")
            text(j, i, sprintf("%.2f", roe_matrix[i, j]), col = text_color, cex = 0.8)
          }
        }
      }
      
      axis(1, at = 1:n_tissues, labels = tissues, las = 2, cex.axis = 0.8)
      axis(4, at = 1:n_celltypes, labels = celltypes, las = 2, cex.axis = 0.8)
    }
    
    # 添加总标题
    mtext(title, outer = TRUE, cex = 1.5, font = 2, line = -2)
    
    # 添加图例
    add_base_legend()
    
    dev.off()
    
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    stop("Both grid and base plotting methods failed: ", e$message)
  })
}

#' 添加标题和图例 (Grid版本)
add_title_and_legend <- function(title) {
  pushViewport(viewport())
  
  # 总标题
  grid.text(title, x = 0.5, y = 0.95, just = "center",
            gp = gpar(fontsize = 16, fontface = "bold"))
  
  # 图例
  legend_y <- 0.08
  legend_items <- c(
    "+++: Strong enrichment (Ro/e ≥ 2.0)",
    "+: Moderate enrichment (1.2 ≤ Ro/e < 2.0)", 
    "+/-: Near expected (0.8 < Ro/e < 1.2)",
    "-: Moderate depletion (0.5 ≤ Ro/e ≤ 0.8)",
    "---: Strong depletion (Ro/e < 0.5)"
  )
  
  legend_colors <- c("#8B0000", "#DC143C", "#FFFACD", "#87CEEB", "#4682B4")
  
  grid.text("Symbol Legend (STARTRAC-dist)", x = 0.5, y = legend_y + 0.03, just = "center",
            gp = gpar(fontsize = 12, fontface = "bold"))
  
  item_width <- 0.18
  for(i in 1:length(legend_items)) {
    x_pos <- 0.1 + (i - 1) * item_width
    
    grid.rect(x = x_pos, y = legend_y, width = 0.02, height = 0.02,
              gp = gpar(fill = legend_colors[i], col = "black"))
    grid.text(legend_items[i], x = x_pos + 0.025, y = legend_y, just = "left",
              gp = gpar(fontsize = 9))
  }
  
  grid.text("Method: STARTRAC-dist uses Ro/e = Observed/Expected from chi-square test", 
            x = 0.5, y = 0.02, just = "center", gp = gpar(fontsize = 10, col = "gray30"))
  
  popViewport()
}

#' 添加基础图例 (Base版本)
add_base_legend <- function() {
  par(fig = c(0, 1, 0, 1), oma = c(4, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  
  legend("bottom", 
         legend = c("+++: ≥2.0", "+: 1.2-2.0", "+/-: 0.8-1.2", "-: 0.5-0.8", "---: <0.5"),
         title = "STARTRAC-dist Symbol Legend",
         fill = c("#8B0000", "#DC143C", "#FFFACD", "#87CEEB", "#4682B4"),
         border = "black", horiz = TRUE, cex = 0.9, title.cex = 1.1)
}

# =============================================================================
# 主分析管道 (Main Analysis Pipeline)
# =============================================================================

#' 整合版组织分布分析 (主函数)
#' @param seurat_obj Seurat对象
#' @param celltype_col 细胞类型列名
#' @param tissue_col 组织位置列名  
#' @param patient_col 患者/样本列名
#' @param output_dir 输出目录
#' @param min_cells 最小细胞数阈值
#' @param statistical_test 统计检验方法
#' @param create_plots 是否创建可视化图像
#' @param plot_method 绘图方法选择
#' @return 分析结果列表
integrated_startrac_analysis <- function(
    seurat_obj,
    celltype_col,
    tissue_col, 
    patient_col,
    output_dir = NULL,
    min_cells = 50,
    statistical_test = "chisq",
    create_plots = TRUE,
    plot_method = "auto"
) {
  
  message("=================================================================")
  message("    INTEGRATED STARTRAC-DIST ANALYSIS PIPELINE")
  message("=================================================================")
  
  # 1. 初始化和验证
  if(is.null(output_dir)) {
    output_dir <- paste0("./startrac_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 检查必需列
  meta_data <- seurat_obj@meta.data
  required_cols <- c(celltype_col, tissue_col, patient_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta_data)]
  
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # 2. 数据预处理
  message("Step 1: Data preprocessing and quality control...")
  
  analysis_data <- data.frame(
    clone.id = rownames(meta_data),
    patient = as.character(meta_data[[patient_col]]),
    majorCluster = as.character(meta_data[[celltype_col]]),
    loc = as.character(meta_data[[tissue_col]]),
    Cell_Name = "Cell",
    stringsAsFactors = FALSE
  )
  
  # 数据清理
  original_cells <- nrow(analysis_data)
  analysis_data <- analysis_data[complete.cases(analysis_data), ]
  analysis_data <- analysis_data[
    analysis_data$majorCluster != "" & 
      analysis_data$loc != "" & 
      analysis_data$patient != "" &
      !analysis_data$majorCluster %in% c("NA", "Unknown", "unknown", "Unassigned") &
      !analysis_data$loc %in% c("NA", "Unknown", "unknown"), ]
  
  cells_after_cleaning <- nrow(analysis_data)
  removal_rate <- round((original_cells - cells_after_cleaning) / original_cells * 100, 2)
  
  message(sprintf("Data cleaning: %d → %d cells (%.1f%% removed)", 
                  original_cells, cells_after_cleaning, removal_rate))
  
  # 移除低频细胞类型
  celltype_counts <- table(analysis_data$majorCluster)
  low_count_celltypes <- names(celltype_counts)[celltype_counts < max(5, min_cells/10)]
  if(length(low_count_celltypes) > 0) {
    message("Removing low-count cell types: ", paste(head(low_count_celltypes, 3), collapse = ", "))
    analysis_data <- analysis_data[!analysis_data$majorCluster %in% low_count_celltypes, ]
  }
  
  # 验证数据充足性
  unique_tissues <- unique(analysis_data$loc)
  if(length(unique_tissues) < 2) {
    stop("Need at least 2 tissue locations. Found: ", paste(unique_tissues, collapse = ", "))
  }
  
  message(sprintf("Final dataset: %d cells, %d cell types, %d tissues, %d samples",
                  nrow(analysis_data), 
                  length(unique(analysis_data$majorCluster)),
                  length(unique_tissues),
                  length(unique(analysis_data$patient))))
  
  # 3. 计算Ro/e指数 (STARTRAC核心)
  message("Step 2: Calculating Ro/e indices using STARTRAC method...")
  
  roe_result <- tryCatch({
    calTissueDist(analysis_data,
                  byPatient = FALSE,
                  colname.cluster = "majorCluster",
                  colname.patient = "patient", 
                  colname.tissue = "loc",
                  method = statistical_test,
                  min.rowSum = min_cells)
  }, error = function(e) {
    stop("Failed to calculate Ro/e indices: ", e$message)
  })
  
  write.csv(roe_result, file.path(output_dir, "roe_indices.csv"))
  
  # 4. 统计显著性检验 (符合STARTRAC原文方法)
  message("Step 3: Statistical significance testing...")
  
  contingency_table <- table(analysis_data$majorCluster, analysis_data$loc)
  
  # 全局卡方检验
  global_chisq_result <- tryCatch({
    expected_from_margins <- outer(rowSums(contingency_table), colSums(contingency_table)) / sum(contingency_table)
    min_expected <- min(expected_from_margins)
    
    if(min_expected >= 1 && sum(expected_from_margins < 5) / length(expected_from_margins) <= 0.2) {
      chisq.test(contingency_table)
    } else {
      chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
    }
  }, error = function(e) {
    message("Warning in chi-square test: ", e$message)
    return(NULL)
  })
  
  # 处理统计结果
  significance_results <- NULL
  if(!is.null(global_chisq_result)) {
    global_significant <- global_chisq_result$p.value < 0.05
    
    message(sprintf("Global chi-square test: χ² = %.3f, p = %.6f (%s)", 
                    global_chisq_result$statistic, 
                    global_chisq_result$p.value,
                    ifelse(global_significant, "SIGNIFICANT", "NOT SIGNIFICANT")))
    
    # 细胞类型级别的显著性评估
    if(global_significant) {
      expected_matrix <- global_chisq_result$expected
      observed_matrix <- as.matrix(contingency_table)
      
      chi_contributions <- ((observed_matrix - expected_matrix)^2) / expected_matrix
      celltype_chi_contrib <- rowSums(chi_contributions)
      std_residuals <- (observed_matrix - expected_matrix) / sqrt(expected_matrix)
      max_abs_residuals <- apply(abs(std_residuals), 1, max)
      
      significance_results <- data.frame(
        CellType = rownames(roe_result),
        Chi_Contribution = celltype_chi_contrib,
        Max_Std_Residual = max_abs_residuals,
        Significant = celltype_chi_contrib > qchisq(0.95, df = ncol(contingency_table) - 1),
        stringsAsFactors = FALSE
      )
      
      significance_results$Significance_Level <- ifelse(
        significance_results$Max_Std_Residual > qnorm(0.9995), "***",
        ifelse(significance_results$Max_Std_Residual > qnorm(0.995), "**",
               ifelse(significance_results$Max_Std_Residual > qnorm(0.975), "*", "ns"))
      )
      
      write.csv(significance_results, file.path(output_dir, "significance_results.csv"), row.names = FALSE)
      
      n_significant <- sum(significance_results$Significant)
      message(sprintf("Significant cell types: %d/%d", n_significant, nrow(significance_results)))
    }
  }
  
  # 5. 创建可视化 (集成稳健版本)
  if(create_plots) {
    message("Step 4: Creating visualizations...")
    
    # 主要STARTRAC-dist图像
    create_startrac_dist_plot(
      roe_matrix = roe_result,
      significance_results = significance_results,
      output_path = file.path(output_dir, "startrac_dist_analysis.pdf"),
      title = "STARTRAC-dist Tissue Distribution Analysis",
      method = plot_method
    )
    
    # 补充可视化
    create_supplementary_plots(roe_result, analysis_data, output_dir)
  }
  
  # 6. 生成综合报告
  message("Step 5: Generating comprehensive report...")
  generate_analysis_report(roe_result, significance_results, analysis_data, output_dir)
  
  # 7. 返回结果
  results <- list(
    roe_indices = roe_result,
    significance_results = significance_results,
    raw_data = analysis_data,
    global_test = global_chisq_result,
    output_directory = output_dir,
    parameters = list(
      celltype_col = celltype_col,
      tissue_col = tissue_col,
      patient_col = patient_col,
      min_cells = min_cells,
      statistical_test = statistical_test
    )
  )
  
  message("=== ANALYSIS COMPLETED SUCCESSFULLY ===")
  message("Output directory: ", output_dir)
  message("Files generated: ", length(list.files(output_dir, recursive = TRUE)))
  
  return(results)
}

# =============================================================================
# 辅助分析和可视化函数 (Auxiliary Analysis & Visualization)
# =============================================================================

#' 创建增强版热图 (Enhanced Heatmap with Annotations)
#' @param roe_result Ro/e指数矩阵
#' @param significance_results 显著性结果
#' @param output_dir 输出目录
create_enhanced_roe_heatmap <- function(roe_result, significance_results, output_dir) {
  
  message("Creating enhanced Ro/e heatmap with annotations...")
  
  tryCatch({
    # 检查是否有ComplexHeatmap包
    if(!require("ComplexHeatmap", quietly = TRUE)) {
      message("ComplexHeatmap not available, creating alternative heatmap...")
      create_alternative_heatmap(roe_result, significance_results, output_dir)
      return()
    }
    
    # 准备细胞类型分组（基于名称模式）
    celltype_groups <- sapply(rownames(roe_result), function(x) {
      x_lower <- tolower(x)
      if(grepl("^t[_ ]|^cd[48]|tcr|th[0-9]|treg|cytotoxic", x_lower)) return("T_cells")
      if(grepl("^b[_ ]|plasma|^cd19|^cd20|antibody", x_lower)) return("B_cells")
      if(grepl("macro|monocyte|^m[0-9]|dendritic|dc[_ ]|^cd14|^cd16", x_lower)) return("Myeloid")
      if(grepl("fibroblast|smc|smooth.*muscle|pericyte|myofibro", x_lower)) return("Stromal")
      if(grepl("epithelial|endothelial|acinar|basal|ciliated", x_lower)) return("Epithelial")
      if(grepl("neutrophil|eosinophil|basophil|mast", x_lower)) return("Granulocytes")
      return("Other")
    })
    
    # 颜色函数
    col_fun <- colorRamp2(
      c(min(roe_result, na.rm = TRUE), 1, max(roe_result, na.rm = TRUE)), 
      c("#2166AC", "white", "#B2182B")
    )
    
    # 创建注释
    row_annotation <- NULL
    if(!is.null(significance_results)) {
      # 创建显著性注释向量
      sig_annotation <- rep("ns", nrow(roe_result))
      names(sig_annotation) <- rownames(roe_result)
      
      # 标记显著的细胞类型
      for(i in 1:nrow(significance_results)) {
        celltype <- significance_results$CellType[i]
        if(celltype %in% names(sig_annotation)) {
          sig_annotation[celltype] <- significance_results$Significance_Level[i]
        }
      }
      
      # 颜色方案
      sig_colors <- c(
        "***" = "#8B0000",  # 深红色 - 高显著性
        "**" = "#CD5C5C",   # 中红色 - 中等显著性 
        "*" = "#F0E68C",    # 淡黄色 - 轻度显著
        "ns" = "#FFFFFF"    # 白色 - 无显著性
      )
      
      # 细胞类型分组注释颜色
      group_colors <- c(
        "T_cells" = "#FF6B6B", "B_cells" = "#4ECDC4", "Myeloid" = "#45B7D1",
        "Stromal" = "#96CEB4", "Epithelial" = "#FFEAA7", "Granulocytes" = "#DDA0DD",
        "Other" = "#D3D3D3"
      )
      
      # 创建行注释
      row_annotation <- rowAnnotation(
        `Cell Group` = celltype_groups,
        `Significance` = sig_annotation,
        col = list(
          `Cell Group` = group_colors[names(group_colors) %in% unique(celltype_groups)],
          `Significance` = sig_colors
        ),
        annotation_legend_param = list(
          `Cell Group` = list(title = "Cell Type\nGroup"),
          `Significance` = list(title = "Statistical\nSignificance")
        ),
        width = unit(1.5, "cm")
      )
    }
    
    # 创建PDF文件
    pdf(file.path(output_dir, "enhanced_roe_heatmap_pairwise.pdf"), width = 14, height = max(8, nrow(roe_result) * 0.4))
    
    # 生成复杂热图
    ht <- Heatmap(
      as.matrix(roe_result),
      name = "Ro/e Index",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "right",
      column_names_side = "top",
      right_annotation = row_annotation,
      
      # 自定义单元格内容
      cell_fun = function(j, i, x, y, width, height, fill) {
        value <- roe_result[i, j]
        if(!is.na(value) && is.finite(value)) {
          # 添加Ro/e数值标签
          text_color <- ifelse(abs(value - 1) > 0.3, "white", "black")
          grid.text(sprintf("%.2f", value), x, y, 
                    gp = gpar(fontsize = 9, col = text_color))
          
          # 添加符号标记
          symbol <- roe_to_symbol(value)
          if(symbol %in% c("+++", "---")) {
            # 对于极端值添加小标记
            grid.circle(x + 0.35 * width, y + 0.35 * height, r = unit(1, "mm"), 
                        gp = gpar(fill = ifelse(symbol == "+++", "red", "blue"), 
                                  col = NA))
          }
        }
      },
      
      # 热图图例参数
      heatmap_legend_param = list(
        title = "Ro/e Index",
        at = c(min(roe_result, na.rm = TRUE), 1, max(roe_result, na.rm = TRUE)),
        labels = c("Depleted", "Expected", "Enriched"),
        legend_direction = "vertical",
        legend_width = unit(6, "cm"),
        title_position = "topcenter",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10)
      ),
      
      # 标题
      column_title = "Enhanced Tissue Distribution Analysis with Statistical Annotations",
      column_title_gp = gpar(fontsize = 16, fontface = "bold"),
      row_title = "Cell Types", 
      row_title_gp = gpar(fontsize = 14),
      
      # 添加边框
      border = TRUE,
      border_gp = gpar(col = "black", lty = 1)
    )
    
    # 绘制热图
    draw(ht)
    
    # 添加说明文字
    grid.text("Red/Blue dots: Extreme enrichment/depletion (Ro/e ≥2.0 or ≤0.5)", 
              x = 0.5, y = 0.05, 
              gp = gpar(fontsize = 10, col = "gray40"))
    grid.text("Significance levels: *** P<0.001, ** P<0.01, * P<0.05", 
              x = 0.5, y = 0.02, 
              gp = gpar(fontsize = 10, col = "gray40"))
    
    dev.off()
    
    message("Enhanced Ro/e heatmap saved successfully")
    
  }, error = function(e) {
    message("Failed to create enhanced heatmap: ", e$message)
    message("Creating alternative heatmap...")
    create_alternative_heatmap(roe_result, significance_results, output_dir)
  })
}

#' 创建备用热图 (当ComplexHeatmap不可用时)
create_alternative_heatmap <- function(roe_result, significance_results, output_dir) {
  
  tryCatch({
    pdf(file.path(output_dir, "enhanced_roe_heatmap_pairwise.pdf"), width = 12, height = 10)
    
    # 使用pheatmap作为备用
    if(require("pheatmap", quietly = TRUE)) {
      
      # 创建注释数据框
      annotation_row <- NULL
      if(!is.null(significance_results)) {
        annotation_row <- data.frame(
          Significance = rep("ns", nrow(roe_result)),
          row.names = rownames(roe_result)
        )
        
        for(i in 1:nrow(significance_results)) {
          celltype <- significance_results$CellType[i]
          if(celltype %in% rownames(annotation_row)) {
            annotation_row[celltype, "Significance"] <- significance_results$Significance_Level[i]
          }
        }
        
        # 注释颜色
        ann_colors <- list(
          Significance = c("***" = "#8B0000", "**" = "#CD5C5C", "*" = "#F0E68C", "ns" = "#FFFFFF")
        )
      }
      
      pheatmap(as.matrix(roe_result),
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               display_numbers = TRUE,
               number_format = "%.2f",
               fontsize_number = 8,
               annotation_row = annotation_row,
               annotation_colors = if(exists("ann_colors")) ann_colors else NULL,
               main = "Enhanced Ro/e Index Heatmap",
               cellwidth = 20,
               cellheight = 15)
      
    } else {
      # 最基本的热图
      heatmap(as.matrix(roe_result),
              col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
              main = "Ro/e Index Heatmap",
              margins = c(8, 8))
    }
    
    dev.off()
    message("Alternative heatmap created successfully")
    
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    message("Failed to create alternative heatmap: ", e$message)
  })
}

#' 创建补充可视化图像
create_supplementary_plots <- function(roe_result, analysis_data, output_dir, significance_results = NULL) {
  
  # 1. 创建增强版热图 (重要!)
  create_enhanced_roe_heatmap(roe_result, significance_results, output_dir)
  
  # 2. 柱状图
  tryCatch({
    roe_df <- as.data.frame(roe_result) %>%
      rownames_to_column("CellType") %>%
      pivot_longer(cols = -CellType, names_to = "Tissue", values_to = "Roe_Index") %>%
      filter(is.finite(Roe_Index))
    
    p_bar <- ggplot(roe_df, aes(x = CellType, y = Roe_Index, fill = Tissue)) +
      geom_col(position = "dodge", width = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
      scale_fill_viridis_d() +
      labs(title = "Cell Type Distribution Across Tissues (Ro/e Index)",
           subtitle = "Values > 1: Enriched, Values < 1: Depleted",
           x = "Cell Type", y = "Ro/e Index", fill = "Tissue Location") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, "roe_barplot.pdf"), p_bar, width = 12, height = 8)
  }, error = function(e) message("Warning: Failed to create bar plot - ", e$message))
  
  # 3. 细胞比例堆叠图
  tryCatch({
    prop_data <- analysis_data %>%
      group_by(loc, majorCluster) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(loc) %>%
      mutate(proportion = count / sum(count))
    
    p_stack <- ggplot(prop_data, aes(x = loc, y = proportion, fill = majorCluster)) +
      geom_col() +
      scale_fill_viridis_d() +
      labs(title = "Cell Type Composition by Tissue Location",
           x = "Tissue Location", y = "Proportion", fill = "Cell Type") +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, "cell_composition_stacked.pdf"), p_stack, width = 10, height = 6)
  }, error = function(e) message("Warning: Failed to create stacked plot - ", e$message))
  
  # 4. 相关性热图
  if(ncol(roe_result) >= 2 && nrow(roe_result) >= 2) {
    tryCatch({
      valid_roe_data <- roe_result[complete.cases(roe_result), ]
      if(nrow(valid_roe_data) >= 2) {
        correlation_matrix <- cor(t(valid_roe_data), use = "complete.obs", method = "pearson")
        
        pdf(file.path(output_dir, "celltype_correlation.pdf"), width = 10, height = 10)
        corrplot::corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust",
                           tl.cex = 0.8, tl.col = "black", tl.srt = 45,
                           title = "Cell Type Distribution Pattern Correlation", mar = c(0,0,2,0))
        dev.off()
      }
    }, error = function(e) message("Warning: Failed to create correlation plot - ", e$message))
  }
}

#' 生成综合分析报告
generate_analysis_report <- function(roe_result, significance_results, analysis_data, output_dir) {
  
  report_file <- file.path(output_dir, "comprehensive_analysis_report.txt")
  
  sink(report_file)
  
  cat("=======================================================\n")
  cat("    INTEGRATED STARTRAC-DIST ANALYSIS REPORT\n")
  cat("=======================================================\n\n")
  
  cat("Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total Cells Analyzed: ", nrow(analysis_data), "\n")
  cat("Cell Types: ", length(unique(analysis_data$majorCluster)), "\n")
  cat("Tissue Locations: ", length(unique(analysis_data$loc)), "\n")
  cat("Patients/Samples: ", length(unique(analysis_data$patient)), "\n\n")
  
  cat("STARTRAC-dist Method Summary:\n")
  cat("- Approach: Global chi-square test with cell-type decomposition\n")
  cat("- Ro/e Index: Ratio of Observed to Expected cell frequencies\n")
  cat("- Interpretation: >1 enriched, <1 depleted, =1 expected distribution\n\n")
  
  cat("Ro/e Index Results:\n")
  cat("-------------------\n")
  print(roe_result)
  cat("\n")
  
  if(!is.null(significance_results)) {
    sig_count <- sum(significance_results$Significant)
    cat("Statistical Significance Summary:\n")
    cat("---------------------------------\n")
    cat("Significant cell types: ", sig_count, "/", nrow(significance_results), 
        " (", round(sig_count/nrow(significance_results)*100, 1), "%)\n")
    
    if(sig_count > 0) {
      cat("\nTop significant cell types:\n")
      top_sig <- head(significance_results[significance_results$Significant, ], 5)
      for(i in 1:nrow(top_sig)) {
        cat("  ", i, ". ", top_sig$CellType[i], " (", top_sig$Significance_Level[i], ")\n", sep="")
      }
    }
  }
  
  cat("\nInterpretation Guidelines:\n")
  cat("- +++: Strong enrichment (Ro/e ≥ 2.0)\n")
  cat("- +: Moderate enrichment (1.2 ≤ Ro/e < 2.0)\n")
  cat("- +/-: Near expected (0.8 < Ro/e < 1.2)\n")
  cat("- -: Moderate depletion (0.5 ≤ Ro/e ≤ 0.8)\n")
  cat("- ---: Strong depletion (Ro/e < 0.5)\n\n")
  
  cat("Files Generated:\n")
  cat("- startrac_dist_analysis.pdf: Main STARTRAC-dist visualization\n")
  cat("- roe_indices.csv: Numerical Ro/e index matrix\n")
  cat("- significance_results.csv: Statistical significance testing results\n")
  cat("- roe_barplot.pdf: Bar chart visualization\n")
  cat("- cell_composition_stacked.pdf: Tissue composition analysis\n")
  cat("- celltype_correlation.pdf: Cell type pattern correlation\n")
  
  sink()
  
  message("Comprehensive report saved to: ", report_file)
}

# =============================================================================
# 便捷包装函数和使用示例 (Convenience Wrappers & Usage Examples)
# =============================================================================

#' 快速启动分析 (一键式分析)
#' @param seurat_obj Seurat对象
#' @param celltype_col 细胞类型列名
#' @param tissue_col 组织列名
#' @param patient_col 样本列名
#' @param output_dir 输出目录（可选）
quick_startrac_analysis <- function(seurat_obj, celltype_col, tissue_col, patient_col, output_dir = NULL) {
  
  # 显示数据概览
  message("=== DATA OVERVIEW ===")
  show_metadata_info(seurat_obj, show_preview = FALSE)
  
  # 执行分析
  results <- integrated_startrac_analysis(
    seurat_obj = seurat_obj,
    celltype_col = celltype_col,
    tissue_col = tissue_col,
    patient_col = patient_col,
    output_dir = output_dir,
    min_cells = 30,
    statistical_test = "chisq",
    create_plots = TRUE,
    plot_method = "auto"
  )
  
  # 快速结果解读
  message("\n=== QUICK RESULTS SUMMARY ===")
  interpret_results_summary(results)
  
  return(results)
}

#' 结果快速解读
interpret_results_summary <- function(results) {
  
  roe_matrix <- results$roe_indices
  
  # 找出最显著的模式
  max_enrichment <- apply(roe_matrix, 1, function(x) max(x, na.rm = TRUE))
  min_depletion <- apply(roe_matrix, 1, function(x) min(x, na.rm = TRUE))
  
  # 最富集的模式
  top_enriched_idx <- which.max(max_enrichment)
  if(length(top_enriched_idx) > 0) {
    ct <- names(max_enrichment)[top_enriched_idx]
    tissue <- names(which.max(roe_matrix[ct, ]))
    value <- max_enrichment[top_enriched_idx]
    message("Most enriched pattern: ", ct, " in ", tissue, " (Ro/e = ", round(value, 2), ")")
  }
  
  # 最耗竭的模式  
  top_depleted_idx <- which.min(min_depletion)
  if(length(top_depleted_idx) > 0) {
    ct <- names(min_depletion)[top_depleted_idx]
    tissue <- names(which.min(roe_matrix[ct, ]))
    value <- min_depletion[top_depleted_idx]
    message("Most depleted pattern: ", ct, " in ", tissue, " (Ro/e = ", round(value, 2), ")")
  }
  
  # 统计概览
  n_enriched <- sum(roe_matrix > 1.2, na.rm = TRUE)
  n_depleted <- sum(roe_matrix < 0.8, na.rm = TRUE)
  total_comparisons <- sum(is.finite(roe_matrix))
  
  message("Distribution patterns: ", n_enriched, " enriched, ", n_depleted, " depleted out of ", total_comparisons, " comparisons")
  
  if(!is.null(results$significance_results)) {
    n_sig <- sum(results$significance_results$Significant)
    message("Statistical significance: ", n_sig, "/", nrow(results$significance_results), " cell types show significant tissue distribution patterns")
  }
}

# =============================================================================
# 使用示例 (Usage Examples)
# =============================================================================

# 示例1：基本使用
results <- quick_startrac_analysis(
  seurat_obj = T_object,
  celltype_col = "Annotation",
  tissue_col = "tissue",
  patient_col = "sample"
)

# 示例2：完整参数控制
# results <- integrated_startrac_analysis(
#   seurat_obj = your_seurat_object,
#   celltype_col = "detailed_celltype",
#   tissue_col = "tissue_location",
#   patient_col = "patient_id", 
#   output_dir = "./my_startrac_analysis",
#   min_cells = 100,
#   statistical_test = "chisq",
#   create_plots = TRUE,
#   plot_method = "grid"
# )

# 示例3：仅创建可视化（已有数据）
# create_startrac_dist_plot(
#   roe_matrix = your_roe_matrix,
#   output_path = "custom_startrac_plot.pdf",
#   title = "Custom STARTRAC Analysis",
#   method = "auto"
# )

message("=================================================================")
message("    INTEGRATED STARTRAC-DIST ANALYSIS SYSTEM LOADED")
message("=================================================================")
message("Main functions:")
message("  • quick_startrac_analysis() - One-click analysis")
message("  • integrated_startrac_analysis() - Full parameter control") 
message("  • create_startrac_dist_plot() - Standalone visualization")
message("  • show_metadata_info() - Examine your data columns")
message("")
message("Key improvements:")
message("  ✓ Robust error handling with automatic fallback methods")
message("  ✓ STARTRAC-compliant statistical methodology") 
message("  ✓ High-quality publication-ready visualizations")
message("  ✓ Comprehensive analysis reports")
message("  ✓ Flexible parameter configuration")
message("")
message("Ready to analyze! Use quick_startrac_analysis() to get started.")
message("=================================================================")



######################################################################################################
######################################################################################################
######################################################################################################
#' 使用DimPlot创建多条件栅格化可视化
#' Create multi-condition rasterized visualization using DimPlot
#'
#' @param seurat_obj Seurat对象
#' @param condition_col 条件列名称
#' @param cell_type_col 细胞类型列名称  
#' @param reduction 降维方法
#' @param pt_size 点的大小
#' @param ncol 图形列数
#' @param raster 是否使用栅格化
#' @param alpha 点的透明度
#' @param label 是否显示聚类标签
#' @param label_size 标签大小
#' @param save_pdf 是否保存PDF文件
#' @param output_prefix 输出文件前缀
#' @return ggplot对象列表
create_multiCondition_dimplot <- function(seurat_obj,
                                          condition_col = "tissue",
                                          cell_type_col = "Annotation", 
                                          reduction = "umap",
                                          pt_size = 0.1,
                                          ncol = 3,
                                          raster = TRUE,
                                          alpha = 2,
                                          label = FALSE,
                                          label_size = 4,
                                          save_pdf = TRUE,
                                          output_prefix = "multiCondition_dimplot") {
  
  # 检查必需的列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  # 检查降维结果是否存在
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("降维结果", reduction, "不存在，请先运行", reduction, "分析"))
  }
  
  # 获取唯一的条件
  conditions <- unique(seurat_obj@meta.data[[condition_col]])
  conditions <- conditions[!is.na(conditions)]
  conditions <- sort(conditions)
  
  message(paste("发现", length(conditions), "个条件:", paste(conditions, collapse = ", ")))
  
  # 创建绘图列表
  plot_list <- list()
  
  # 设置细胞类型为idents用于DimPlot
  Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]
  
  # 1. 总体图：所有样本 (Overall plot: All samples)
  p_all <- DimPlot(seurat_obj, 
                   reduction = reduction,
                   raster = raster,
                   pt.size = pt_size,
                   alpha = alpha,
                   label = label,
                   label.size = label_size) +
    labs(title = "All samples") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      panel.grid = element_blank(),
      panel.border = element_blank()  # ✅ 去掉边框
    ) +
    guides(color = guide_legend(title = "Cell types", 
                                override.aes = list(size = 3, alpha = 1)))
  
  plot_list[["All_samples"]] <- p_all
  
  # 2. 各条件的子图 (Subplots for each condition)
  for (cond in conditions) {
    # 使用逻辑索引提取该条件的细胞
    condition_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[condition_col]] == cond]
    
    # 创建该条件的子集
    condition_subset <- subset(seurat_obj, cells = condition_cells)
    
    # 设置细胞类型为idents
    Idents(condition_subset) <- condition_subset@meta.data[[cell_type_col]]
    
    # 创建条件特异性图
    p_cond <- DimPlot(condition_subset, 
                      reduction = reduction,
                      raster = raster,
                      pt.size = pt_size,
                      alpha = alpha,
                      label = label,
                      label.size = label_size) +
      labs(title = cond) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none",  # 在子图中隐藏图例
        panel.grid = element_blank(),
        panel.border = element_blank() 
      )
    
    plot_list[[cond]] <- p_cond
  }
  
  # 3. 组合图形 (Combine plots)
  # 创建图例 (Create legend)
  legend_plot <- get_legend(p_all + 
                              theme(legend.position = "right",
                                    legend.direction = "vertical",
                                    legend.key.size = unit(0.6, "cm")))
  
  # 移除总体图的图例 (Remove legend from overall plot)
  p_all_no_legend <- p_all + theme(legend.position = "none")
  plot_list[["All_samples"]] <- p_all_no_legend
  
  # 计算布局 (Calculate layout)
  n_plots <- length(plot_list)
  
  # 创建组合图 (Create combined plot)
  if (n_plots <= 4) {
    combined_plot <- wrap_plots(plot_list, ncol = 2)
  } else {
    combined_plot <- wrap_plots(plot_list, ncol = ncol)
  }
  
  # 添加图例到右侧 (Add legend to the right)
  final_plot <- combined_plot | legend_plot
  final_plot <- final_plot + plot_layout(widths = c(4, 1))
  
  # 保存图形 (Save plots)
  if (save_pdf) {
    # 计算PDF尺寸 (Calculate PDF dimensions)
    pdf_width <- ifelse(n_plots <= 4, 12, 15)
    pdf_height <- ifelse(n_plots <= 4, 8, 10)
    
    # 保存组合图 (Save combined plot)
    ggsave(paste0(output_prefix, "_combined.pdf"), 
           final_plot, 
           width = pdf_width, 
           height = pdf_height, 
           dpi = 300)
    
    # 保存单独的图 (Save individual plots)
    for (plot_name in names(plot_list)) {
      ggsave(paste0(output_prefix, "_", plot_name, ".pdf"), 
             plot_list[[plot_name]], 
             width = 6, height = 5, dpi = 300)
    }
    
    message(paste("图形已保存:", paste0(output_prefix, "_combined.pdf")))
  }
  
  return(list(
    combined_plot = final_plot,
    individual_plots = plot_list,
    cell_counts = table(seurat_obj@meta.data[[condition_col]], 
                        seurat_obj@meta.data[[cell_type_col]])
  ))
}

# ====================================================================
# 诊断和检查函数 (Diagnostic and Check Functions)
# ====================================================================

#' 检查Seurat对象的数据结构
#' Check Seurat object data structure
#'
#' @param seurat_obj Seurat对象
#' @param condition_col 条件列名称
#' @param cell_type_col 细胞类型列名称
check_seurat_structure <- function(seurat_obj, 
                                   condition_col = "tissue", 
                                   cell_type_col = "Annotation") {
  
  cat("=== Seurat对象数据结构检查 ===\n")
  
  # 1. 基本信息
  cat("\n1. 基本信息:\n")
  cat(paste("   总细胞数:", ncol(seurat_obj), "\n"))
  cat(paste("   基因数:", nrow(seurat_obj), "\n"))
  
  # 2. 元数据列
  cat("\n2. 元数据列:\n")
  meta_cols <- colnames(seurat_obj@meta.data)
  cat(paste("   可用列:", paste(meta_cols, collapse = ", "), "\n"))
  
  # 3. 检查指定列是否存在
  cat("\n3. 指定列检查:\n")
  if (condition_col %in% meta_cols) {
    cat(paste("   ✓ 条件列 '", condition_col, "' 存在\n", sep = ""))
    conditions <- unique(seurat_obj@meta.data[[condition_col]])
    conditions <- conditions[!is.na(conditions)]
    cat(paste("     条件数量:", length(conditions), "\n"))
    cat(paste("     条件列表:", paste(conditions, collapse = ", "), "\n"))
  } else {
    cat(paste("   ✗ 条件列 '", condition_col, "' 不存在！\n", sep = ""))
    cat("     请检查列名或使用以下建议:\n")
    potential_cols <- meta_cols[grepl("tissue|group|condition|sample", meta_cols, ignore.case = TRUE)]
    if(length(potential_cols) > 0) {
      cat(paste("     可能的条件列:", paste(potential_cols, collapse = ", "), "\n"))
    }
  }
  
  if (cell_type_col %in% meta_cols) {
    cat(paste("   ✓ 细胞类型列 '", cell_type_col, "' 存在\n", sep = ""))
    cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])
    cell_types <- cell_types[!is.na(cell_types)]
    cat(paste("     细胞类型数量:", length(cell_types), "\n"))
    cat(paste("     细胞类型列表:", paste(head(cell_types, 10), collapse = ", ")))
    if(length(cell_types) > 10) cat(" ...")
    cat("\n")
  } else {
    cat(paste("   ✗ 细胞类型列 '", cell_type_col, "' 不存在！\n", sep = ""))
    cat("     请检查列名或使用以下建议:\n")
    potential_cols <- meta_cols[grepl("annotation|cell_type|cluster|ident", meta_cols, ignore.case = TRUE)]
    if(length(potential_cols) > 0) {
      cat(paste("     可能的细胞类型列:", paste(potential_cols, collapse = ", "), "\n"))
    }
  }
  
  # 4. 降维结果检查
  cat("\n4. 降维结果检查:\n")
  reductions <- names(seurat_obj@reductions)
  if(length(reductions) > 0) {
    cat(paste("   可用降维方法:", paste(reductions, collapse = ", "), "\n"))
  } else {
    cat("   ✗ 没有发现降维结果，请先运行RunTSNE()或RunUMAP()\n")
  }
  
  # 5. 细胞计数统计
  if (condition_col %in% meta_cols && cell_type_col %in% meta_cols) {
    cat("\n5. 细胞计数统计:\n")
    count_table <- table(seurat_obj@meta.data[[condition_col]], 
                         seurat_obj@meta.data[[cell_type_col]])
    print(count_table)
    
    # 检查是否有空组合
    empty_combinations <- which(count_table == 0, arr.ind = TRUE)
    if(nrow(empty_combinations) > 0) {
      cat("\n   警告: 发现空的组合 (没有细胞):\n")
      for(i in 1:nrow(empty_combinations)) {
        row_name <- rownames(count_table)[empty_combinations[i, 1]]
        col_name <- colnames(count_table)[empty_combinations[i, 2]]
        cat(paste("     ", row_name, " + ", col_name, "\n"))
      }
    }
  }
  
  cat("\n=== 检查完成 ===\n")
}

# ====================================================================
# 使用示例 (Usage Examples)
# ====================================================================

# 使用前先检查数据结构
# check_seurat_structure(T_object, "tissue", "Annotation")

# 基本调用示例
# results <- create_multiCondition_dimplot(
#   seurat_obj = T_object,
#   condition_col = "tissue",
#   cell_type_col = "Annotation", 
#   reduction = "tsne",
#   pt_size = 0.1,
#   raster = TRUE,    # 启用栅格化，类似FeaturePlot的raster=TRUE
#   alpha = 0.8,
#   save_pdf = TRUE,
#   output_prefix = "CRS_dimplot_analysis"
# )

# 高质量发表版本
# results_hq <- create_multiCondition_dimplot(
#   seurat_obj = T_object,
#   condition_col = "tissue",
#   cell_type_col = "Annotation",
#   pt_size = 0.2,     # 较大的点用于发表
#   raster = TRUE,
#   alpha = 0.9,       # 更不透明
#   label = TRUE,      # 显示细胞类型标签
#   label_size = 5,
#   output_prefix = "publication_quality"
# )

# 大数据集优化版本
# results_large <- create_multiCondition_dimplot(
#   seurat_obj = T_object,
#   condition_col = "tissue",
#   cell_type_col = "Annotation",
#   pt_size = 0.05,    # 很小的点避免重叠
#   raster = TRUE,     # 栅格化提高性能
#   alpha = 0.6,       # 较低透明度显示密度
#   output_prefix = "large_dataset_optimized"
# )

# ====================================================================
# 额外的便捷函数 (Additional Convenience Functions)
# ====================================================================

#' 快速预览函数
quick_preview <- function(seurat_obj, condition_col = "tissue", cell_type_col = "Annotation") {
  # 检查列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  results <- create_multiCondition_dimplot(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    cell_type_col = cell_type_col,
    pt_size = 0.05,
    raster = TRUE,
    alpha = 0.6,
    save_pdf = FALSE,  # 不保存文件，只预览
    output_prefix = "preview"
  )
  return(results$combined_plot)
}

#' 发表质量函数
publication_quality <- function(seurat_obj, condition_col = "tissue", cell_type_col = "Annotation") {
  # 检查列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  results <- create_multiCondition_dimplot(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    cell_type_col = cell_type_col,
    pt_size = 0.15,
    raster = TRUE,
    alpha = 0.9,
    label = FALSE,     # 通常发表不需要标签
    save_pdf = TRUE,
    output_prefix = "publication_dimplot"
  )
  return(results)
}

#' 带标签版本函数
labeled_version <- function(seurat_obj, condition_col = "tissue", cell_type_col = "Annotation") {
  # 检查列是否存在
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("条件列", condition_col, "在元数据中不存在"))
  }
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "在元数据中不存在"))
  }
  
  results <- create_multiCondition_dimplot(
    seurat_obj = seurat_obj,
    condition_col = condition_col,
    cell_type_col = cell_type_col,
    pt_size = 0.1,
    raster = TRUE,
    alpha = 0.7,
    label = TRUE,      # 显示细胞类型标签
    label_size = 4,
    save_pdf = TRUE,
    output_prefix = "labeled_dimplot"
  )
  return(results)
}

message("基于DimPlot的多条件栅格化可视化函数已加载完成!")
message("DimPlot-based multi-condition rasterized visualization functions loaded!")
message("现在使用Seurat原生的raster支持，无需额外包依赖")
message("Now using Seurat's native raster support, no additional package dependencies required")
message("已修复WhichCells错误，现在可以正常运行")
message("Fixed WhichCells error, should work properly now")
message("使用 check_seurat_structure(seurat_obj, 'tissue', 'Annotation') 检查数据")
message("Use check_seurat_structure(seurat_obj, 'tissue', 'Annotation') to check your data")

# 检查数据结构
check_seurat_structure(T_object, "tissue", "Annotation")

# 基本调用
results <- create_multiCondition_dimplot(
  seurat_obj = T_object,
  condition_col = "tissue",
  cell_type_col = "Annotation",
  reduction = 'umap',
  pt_size = 0.45,
  alpha = 2
)

# 快速预览
preview_plot <- quick_preview(T_object, "tissue", "Annotation")
print(preview_plot)

######################################################################################################
######################################################################################################
######################################################################################################
diagnosis <- diagnose_masc_convergence(T_object)

results <- scPairwiseMASCAnalysis(
  seurat_obj = T_object,
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  exclude_filter = "tissue_sampling_method == 'scraping'",
  fixed_effects_cols = NULL,  # 添加任何固定效应协变量
  output_dir = "pairwise_MASC_results_0928",
  # 如果您需要少于2个样本也能进行比较，可以降低此阈值
  min_samples = 2
)
props_data <- smart_prepare_plotting_data(
  masc_results = results,  # 先尝试从MASC结果提取
  seurat_obj = T_object,  # 备用重新计算
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  exclude_filter = "tissue_sampling_method != 'scraping'",  # 必须与MASC相同！
  min_samples = 2,
  min_cells = 10,
  min_prop = 0.01
)


# 第三步：验证比例数据
print("比例数据预览:")
head(props_data)

print("与MASC结果中的组织是否一致:")
print("绘图数据中的组织:")
table(props_data$Contrast)

print("MASC结果中的组织:")
if("all_results" %in% names(results)) {
  print(table(results$all_results$Comparison))
}

# 第四步：创建完美对齐的图表
final_result <- create_separate_aligned_plot(
  props_df = props_data,
  stats_results = results$all_results,
  output_dir = "consistent_plots",
  contrast_colors = c(
    "lung parenchyma" = "#2E8B57",
    "sinus" = "#B22222", 
    "nose" = "#CD853F",
    "respiratory airway" = "#4682B4"
  )
)

print("箱子位置信息:")
print(final_result$positions)

######################################################################################################
######################################################################################################
######################################################################################################
# 统计每个组织的唯一样本数量
tissue_sample_counts <- sapply(unique(T_object$tissue), function(t) {
  length(unique(T_object$sample[T_object$tissue == t]))
})
names(tissue_sample_counts) <- unique(T_object$tissue)
print(tissue_sample_counts)

# 或者使用table + unique组合
unique_combinations <- unique(T_object@meta.data[, c("tissue", "sample")])
samples_per_tissue <- table(unique_combinations$tissue)
print(samples_per_tissue)


######################################################################################################
######################################################################################################
######################################################################################################

# 大型单细胞表达热图生成函数
# Comprehensive Single-Cell Expression Heatmap Generator

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

#' 创建大型单细胞表达热图（类似发表图片）
#' Create Large-Scale Single-Cell Expression Heatmap
#' 
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名，默认"Annotation"
#' @param genes_to_plot 要绘制的基因列表，默认NULL（自动选择）
#' @param max_cells_per_type 每种细胞类型最大细胞数，用于下采样
#' @param max_genes 最大基因数量
#' @param output_dir 输出目录
#' @param output_prefix 输出文件前缀
#' @param cell_type_colors 细胞类型颜色映射
#' @param expression_color_scheme 表达量颜色方案
#' @param figure_width 图片宽度（英寸）
#' @param figure_height 图片高度（英寸）
#' @param show_gene_names 是否显示基因名称
#' @param show_cell_names 是否显示细胞名称
#' @param cluster_genes 是否对基因进行聚类
#' @param cluster_cells 是否对细胞进行聚类
create_comprehensive_expression_heatmap <- function(
    seurat_obj,
    cell_type_col = "Annotation",
    genes_to_plot = NULL,
    max_cells_per_type = 500,
    max_genes = 200,
    output_dir = "expression_heatmaps",
    output_prefix = "comprehensive_expression",
    cell_type_colors = NULL,
    expression_color_scheme = c("blue", "white", "red"),
    figure_width = 20,
    figure_height = 16,
    show_gene_names = TRUE,
    show_cell_names = FALSE,
    cluster_genes = TRUE,
    cluster_cells = FALSE
) {
  
  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("开始创建大型表达热图 / Starting comprehensive expression heatmap creation...")
  
  # 1. 数据预处理和质控 ====
  message("Step 1: 数据预处理 / Data preprocessing...")
  
  # 检查细胞类型列是否存在
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("细胞类型列", cell_type_col, "不存在于meta.data中"))
  }
  
  # 移除NA细胞类型
  valid_cells <- !is.na(seurat_obj@meta.data[[cell_type_col]])
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_cells])
  
  # 获取细胞类型信息
  cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])
  cell_types <- cell_types[!is.na(cell_types)]
  message(paste("检测到", length(cell_types), "种细胞类型:", paste(cell_types, collapse = ", ")))
  
  # 2. 基因选择策略 ====
  message("Step 2: 基因选择 / Gene selection...")
  
  if (is.null(genes_to_plot)) {
    # 自动选择基因：结合高变基因和marker基因
    
    # 方法1：获取每种细胞类型的top marker基因
    message("寻找细胞类型marker基因 / Finding cell type markers...")
    
    Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]
    
    tryCatch({
      # 寻找所有细胞类型的marker基因
      markers <- FindAllMarkers(
        seurat_obj, 
        only.pos = TRUE, 
        min.pct = 0.25, 
        logfc.threshold = 0.5,
        max.cells.per.ident = 1000  # 限制每种类型的细胞数以加速计算
      )
      
      # 每种细胞类型选择top基因
      top_markers <- markers %>%
        group_by(cluster) %>%
        top_n(n = min(15, max_genes %/% length(cell_types)), wt = avg_log2FC) %>%
        pull(gene) %>%
        unique()
      
      message(paste("找到", length(top_markers), "个marker基因"))
      
    }, error = function(e) {
      message("Marker基因寻找失败，使用高变基因替代")
      top_markers <- character(0)
    })
    
    # 方法2：高变基因
    if (length(VariableFeatures(seurat_obj)) > 0) {
      variable_genes <- head(VariableFeatures(seurat_obj), max_genes %/% 2)
      message(paste("选择", length(variable_genes), "个高变基因"))
    } else {
      message("计算高变基因...")
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = max_genes %/% 2)
      variable_genes <- VariableFeatures(seurat_obj)
    }
    
    # 合并基因列表
    genes_to_plot <- unique(c(top_markers, variable_genes))
    genes_to_plot <- head(genes_to_plot, max_genes)
    
  } else {
    # 使用提供的基因列表，并过滤存在的基因
    genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(seurat_obj)]
  }
  
  message(paste("最终选择", length(genes_to_plot), "个基因用于热图"))
  
  # 3. 细胞下采样 ====
  message("Step 3: 细胞下采样 / Cell subsampling...")
  
  # 按细胞类型进行分层抽样
  set.seed(42)  # 确保结果可重复
  
  sampled_cells <- c()
  for (ct in cell_types) {
    ct_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == ct]
    
    if (length(ct_cells) > max_cells_per_type) {
      # 随机抽样
      sampled_ct_cells <- sample(ct_cells, max_cells_per_type)
    } else {
      sampled_ct_cells <- ct_cells
    }
    
    sampled_cells <- c(sampled_cells, sampled_ct_cells)
    message(paste("细胞类型", ct, ":", length(sampled_ct_cells), "/", length(ct_cells), "个细胞"))
  }
  
  # 创建下采样后的子集
  heatmap_obj <- subset(seurat_obj, cells = sampled_cells)
  message(paste("下采样后总细胞数:", length(sampled_cells)))
  
  # 4. 数据标准化 ====
  message("Step 4: 数据标准化 / Data scaling...")
  
  # 确保选定的基因进行了标准化
  if (!all(genes_to_plot %in% rownames(heatmap_obj[["RNA"]]@scale.data))) {
    message("对选定基因进行标准化...")
    heatmap_obj <- ScaleData(heatmap_obj, features = genes_to_plot)
  }
  
  # 5. 准备热图数据 ====
  message("Step 5: 准备热图数据 / Preparing heatmap data...")
  
  # 提取标准化表达矩阵
  expression_matrix <- heatmap_obj[["RNA"]]@scale.data[genes_to_plot, ]
  
  # 处理无限值和NA值
  expression_matrix[is.infinite(expression_matrix)] <- 0
  expression_matrix[is.na(expression_matrix)] <- 0
  
  # 设置表达值范围（类似Seurat的默认设置）
  expression_matrix[expression_matrix > 2.5] <- 2.5
  expression_matrix[expression_matrix < -2.5] <- -2.5
  
  # 6. 细胞排序和分组 ====
  message("Step 6: 细胞排序 / Cell ordering...")
  
  # 按细胞类型对细胞进行分组和排序
  cell_metadata <- heatmap_obj@meta.data[colnames(expression_matrix), ]
  cell_metadata$CellType <- cell_metadata[[cell_type_col]]
  
  # 按细胞类型排序细胞
  if (!cluster_cells) {
    # 简单按细胞类型分组
    cell_order <- order(cell_metadata$CellType)
  } else {
    # 在每个细胞类型内进行聚类
    cell_order <- c()
    for (ct in sort(cell_types)) {
      ct_cells <- rownames(cell_metadata)[cell_metadata$CellType == ct]
      if (length(ct_cells) > 1) {
        ct_matrix <- expression_matrix[, ct_cells, drop = FALSE]
        if (ncol(ct_matrix) > 2) {
          # 计算细胞间距离并聚类
          cell_dist <- dist(t(ct_matrix))
          cell_hclust <- hclust(cell_dist)
          ct_order <- ct_cells[cell_hclust$order]
        } else {
          ct_order <- ct_cells
        }
      } else {
        ct_order <- ct_cells
      }
      cell_order <- c(cell_order, match(ct_order, colnames(expression_matrix)))
    }
  }
  
  # 重新排序表达矩阵
  expression_matrix <- expression_matrix[, cell_order]
  cell_metadata <- cell_metadata[colnames(expression_matrix), ]
  
  # 7. 设置颜色方案 ====
  message("Step 7: 设置颜色方案 / Setting color schemes...")
  
  # 细胞类型颜色
  if (is.null(cell_type_colors)) {
    if (length(cell_types) <= 12) {
      cell_type_colors <- brewer.pal(min(12, max(3, length(cell_types))), "Set3")
    } else {
      cell_type_colors <- rainbow(length(cell_types))
    }
    names(cell_type_colors) <- sort(cell_types)
  }
  
  # 表达量颜色
  if (length(expression_color_scheme) == 3) {
    expression_colors <- colorRamp2(
      c(-2.5, 0, 2.5), 
      expression_color_scheme
    )
  } else {
    expression_colors <- colorRamp2(
      seq(-2.5, 2.5, length.out = length(expression_color_scheme)),
      expression_color_scheme
    )
  }
  
  # 8. 创建注释 ====
  message("Step 8: 创建注释 / Creating annotations...")
  
  # 顶部细胞类型注释
  top_annotation <- HeatmapAnnotation(
    `Cell Type` = cell_metadata$CellType,
    col = list(`Cell Type` = cell_type_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    simple_anno_size = unit(8, "mm")
  )
  
  # 9. 生成热图 ====
  message("Step 9: 生成热图 / Generating heatmap...")
  
  # 基因聚类（如果需要）
  if (cluster_genes && nrow(expression_matrix) > 2) {
    gene_clustering <- TRUE
  } else {
    gene_clustering <- FALSE
  }
  
  # 创建主热图
  main_heatmap <- Heatmap(
    expression_matrix,
    name = "Scaled Expression",
    
    # 颜色设置
    col = expression_colors,
    
    # 聚类设置
    cluster_rows = gene_clustering,
    cluster_columns = FALSE,  # 列已经按细胞类型排序
    
    # 显示设置
    show_row_names = show_gene_names,
    show_column_names = show_cell_names,
    
    # 注释
    top_annotation = top_annotation,
    
    # 字体大小
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    
    # 热图大小
    width = unit(figure_width - 4, "inch"),
    height = unit(figure_height - 3, "inch"),
    
    # 添加分割线（在细胞类型之间）
    column_split = cell_metadata$CellType,
    column_gap = unit(1, "mm"),
    
    # 图例设置
    heatmap_legend_param = list(
      title = "Scaled\nExpression",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_height = unit(4, "cm")
    )
  )
  
  # 10. 保存热图 ====
  message("Step 10: 保存热图 / Saving heatmap...")
  
  # PDF格式（矢量图，适合期刊）
  pdf_file <- file.path(output_dir, paste0(output_prefix, "_heatmap.pdf"))
  pdf(pdf_file, width = figure_width, height = figure_height)
  draw(main_heatmap, 
       annotation_legend_side = "right",
       heatmap_legend_side = "right")
  dev.off()
  
  # PNG格式（高分辨率位图）
  png_file <- file.path(output_dir, paste0(output_prefix, "_heatmap.png"))
  png(png_file, width = figure_width * 300, height = figure_height * 300, 
      res = 300, bg = "white")
  draw(main_heatmap,
       annotation_legend_side = "right", 
       heatmap_legend_side = "right")
  dev.off()
  
  # 11. 保存相关数据 ====
  message("Step 11: 保存相关数据 / Saving related data...")
  
  # 保存基因列表
  write.csv(
    data.frame(Gene = rownames(expression_matrix)),
    file.path(output_dir, paste0(output_prefix, "_genes.csv")),
    row.names = FALSE
  )
  
  # 保存细胞信息
  write.csv(
    cell_metadata,
    file.path(output_dir, paste0(output_prefix, "_cell_metadata.csv")),
    row.names = TRUE
  )
  
  # 保存表达矩阵（可选，如果数据量不大）
  if (ncol(expression_matrix) * nrow(expression_matrix) < 1e6) {
    write.csv(
      expression_matrix,
      file.path(output_dir, paste0(output_prefix, "_expression_matrix.csv")),
      row.names = TRUE
    )
  }
  
  message("热图创建完成！/ Heatmap creation completed!")
  message(paste("输出文件保存在:", output_dir))
  message(paste("主要文件:", basename(pdf_file), "和", basename(png_file)))
  
  # 返回关键信息
  result <- list(
    expression_matrix = expression_matrix,
    cell_metadata = cell_metadata,
    genes_used = rownames(expression_matrix),
    cell_types = cell_types,
    cell_type_colors = cell_type_colors,
    output_files = c(pdf_file, png_file)
  )
  
  return(result)
}

#' 快速预览热图（小规模版本）
#' Quick heatmap preview with reduced scale
quick_preview_heatmap <- function(seurat_obj,
                                  cell_type_col = "Annotation",
                                  max_cells_total = 1000,
                                  max_genes = 50) {
  
  message("创建快速预览热图...")
  
  # 简化参数创建预览
  result <- create_comprehensive_expression_heatmap(
    seurat_obj = seurat_obj,
    cell_type_col = cell_type_col,
    max_cells_per_type = max_cells_total %/% length(unique(seurat_obj@meta.data[[cell_type_col]])),
    max_genes = max_genes,
    output_prefix = "preview",
    figure_width = 12,
    figure_height = 8,
    show_gene_names = TRUE
  )
  
  return(result)
}

# 使用示例 / Usage Examples:
#
# 1. 基本使用（自动选择基因）
result <- create_comprehensive_expression_heatmap(
  seurat_obj = T_object,
  cell_type_col = "Annotation",
  output_dir = "publication_heatmaps"
)
#
# # 2. 使用特定基因列表
# marker_genes <- c("CD3E", "CD3D", "CD14", "CD68", "EPCAM", "PECAM1", "COL1A1")
# result <- create_comprehensive_expression_heatmap(
#   seurat_obj = T_object,
#   cell_type_col = "Annotation", 
#   genes_to_plot = marker_genes,
#   max_cells_per_type = 300,
#   output_prefix = "marker_expression"
# )
#
# # 3. 快速预览
# preview <- quick_preview_heatmap(T_object)
#
# # 4. 自定义颜色
# custom_colors <- c(
#   "Epithelial" = "#E41A1C",
#   "Fibroblast" = "#377EB8", 
#   "T" = "#4DAF4A",
#   "Myeloid" = "#984EA3",
#   "Endothelial" = "#FF7F00"
# )
# result <- create_comprehensive_expression_heatmap(
#   seurat_obj = T_object,
#   cell_type_col = "Annotation",
#   cell_type_colors = custom_colors,
#   expression_color_scheme = c("#2166AC", "#F7F7F7", "#B2182B")
# )

######################################################################################################
######################################################################################################
######################################################################################################
