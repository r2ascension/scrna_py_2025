# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(rsvd)
library(harmony)
library(presto)
Sys.setenv(RETICULATE_PYTHON = "E:/ProgramData/anaconda3/python.exe")
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
source("E:/R/scMASC.R")
Markers <- c('FXYD3', 'EPCAM', 'ELF3', 'SERPINF1', 'TSPAN1',
             'SCGB1A1', 'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8',#B
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

Marker_B <- c("IGHD", "TCL1A", "FCER2", "IL4R",
              "MS4A1", "CD27", "BANK1", "TNFRSF13B",
              "AICDA", "RGS13", "GCSAM",
              "IL10", "CD1D", "CD5", "TGFB1",
              "MKI67", "TOP2A",
              "CD38", "SDC1", "IGHA1", "IGHA2",
              "CD38", "SDC1", "IGHG1", "IGHGP",
              "CD38", "SDC1", "MZB1"
)

setwd("E:/R/0301/B_analysis/")#设置工作路径
B_object <- readRDS("B_analyzed_20250625.rds")#读取数据
B_object <- NormalizeData(B_object)#归一化
B_object <- FindVariableFeatures(B_object, selection.method = "vst", nfeatures = 3000)#寻找变异基因
B_object <- ScaleData(B_object)#标准化
B_object <- RunPCA(B_object, npcs = 40)#PCA
B_object <- RunHarmony(
    object = B_object,           
    group.by.vars = c("sample"),   
    theta = c(3.2),
    # lambda = c(14.5),
    sigma = 0.05,
    nclust = 20,
    reduction.use = "pca",
    max_iter = 70,
    cool_down = 10,
    epsilon_harmony = 0.0001,
    monitor = "harmony_score",
    factor = 0.8,
    tol = 1e-10,
    min_delta = 0.0001,
    patience = 10,
    # epsilon_cluster = 0.00001,
    early_stop = TRUE,
    dims = 1:40           
)
# 获取降维结果和元数据
embeddings <- Embeddings(B_object, "harmony")
metadata <- B_object@meta.data

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
# counts <- GetAssayData(B_object, layer = "counts")
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
# pct_var <- B_object[["pca"]]@stdev / sum(B_object[["pca"]]@stdev) * 100
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
# n_sample <- min(10000, ncol(B_object))
# sample_idx <- sample(seq_len(ncol(B_object)), n_sample)

# harmony_scores <- B_object[["harmony"]]@cell.embeddings[sample_idx, 1:20]
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
B_object <- RunUMAP(
    B_object, 
    reduction = "harmony", 
    dims = 1:30, 
    n.neighbors = 30,
    min.dist = 0.3,
    n.trees = 50,
    # learning.rate = 0.2,    
    # n.epochs = 1400,        
    # spread = 1.2,
    # repulsion.strength = 1.1,
    metric = "correlation"
)

# FindNeighbors函数的缩进修正
B_object <- FindNeighbors(
    B_object, 
    reduction = "harmony", 
    dims = 1:30,
    n.trees = 500,
    k.param = 30
) 

B_object <- FindClusters(
    B_object, 
    resolution = 3, 
    algorithm = 2
    )
# 8. 绘制UMAP聚类图
print("Generating plots...")
Idents(B_object) <- B_object$seurat_clusters
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(B_object, reduction = "umap", label = TRUE)
dev.off()
# saveRDS(B_object, "final_analyzed_B_object.rds")
# 9. 保存分析结果
print("Saving analysis results...")

# B_object <- readRDS("E:/R/0224/final_analyzed_B_object.rds")


# 9. 标记基因分析

print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for(marker in Marker_B) {
  if(marker %in% rownames(B_object)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(B_object, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(B_object, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(B_object, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(B_object, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(B_object, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(B_object, reduction = "umap", group.by = "ann_finest_level", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_finest_level")
# p6 <- DimPlot(B_object, reduction = "umap", group.by = "Annotation", raster = TRUE, 
#               pt.size = 0.5) + ggtitle("Annotation")
p1
p2
p3
p4
p5
# p6
dev.off()

tissues <- unique(B_object$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(B_object) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(B_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = B_object, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(B_object$sample)
pdf("sample_umap_split.pdf", width = 12, height = 12)
Idents(B_object) <- "sample"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(B_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = B_object, 
                                                    idents = sample_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(B_object$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(B_object) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(B_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = B_object, 
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
B_object$ann_finest_level_no_na <- B_object$ann_finest_level
B_object$ann_finest_level_no_na[is.na(B_object$ann_finest_level_no_na)] <- "Unknown"
Idents(B_object) <- "ann_finest_level_no_na"

pdf("ann_level_3_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_finest_level_no_na <- unique(B_object$ann_finest_level_no_na)
ann_finest_level_no_na <- ann_finest_level_no_na[!is.na(ann_finest_level_no_na)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_finest_level_no_na) {
  print(DimPlot(B_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(B_object, 
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
print(length(Markers))
print(length(unique(Markers)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Markers)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Marker_B <- unique(Marker_B)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot.pdf", width = 32, height = 18)
dot_plot <- DotPlot(B_object, 
                    features = Marker_B, 
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
# saveRDS(B_object, "B_analyzed_20250625.rds")

# 1. 读取注释文件（无表头）
annotations <- read.csv("Annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- B_object$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
B_object$Annotation <- new_idents
B_object <- SetIdent(B_object, value = new_idents)

cell_colors <- c(
  "Naive_B" = "#00BFFF",          # 深天蓝色 - 幼稚B细胞
  "Memory_B" = "#228B22",         # 森林绿 - 记忆B细胞
  "Cycling_B" = "#FF8C00",        # 深橙色 - 增殖B细胞
  "IgA_Plasma" = "#FF1493",       # 深粉红 - IgA浆细胞
  "IgG_Plasma" = "#8B0000",       # 暗红色 - IgG浆细胞
  "IgA_G_Plasma" = "#9400D3"      # 紫罗兰 - IgA/G双阳性浆细胞
)
Idents(B_object) <- "Annotation_2"  # 先设置identity
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(B_object, 
        reduction = "umap",
        raster = FALSE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        cols = cell_colors
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 6. 绘制分割版本
pdf("cell_types_umap_split.pdf", width = 12, height = 12)
# 获取所有细胞类型
cell_types <- unique(B_object$Annotation)
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 为每个细胞类型单独绘图
for(cell_type in cell_types) {
  # 创建文件名
  
  
  # 绘制该细胞类型的UMAP图
  
  print(DimPlot(B_object, 
                reduction = "umap",
                cells.highlight = WhichCells(B_object, idents = cell_type),
                cols.highlight = cell_colors[cell_type],
                cols = "grey",  # 其他细胞显示为灰色
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(cell_type) +
          theme(legend.text = element_text(size = 12)))
  
}
dev.off()

# 假设我们选择 cell_type 为 "T_cells" 的细胞
# 1. 提取Undefined细胞子集
# 从主Seurat对象中筛选出Annotation标注为"Undefined"的细胞
Undefined_obj <- subset(B_object, subset = Annotation == "Undefined")

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
Undefined_obj <- RunPCA(Undefined_obj, npcs = 50)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
Undefined_obj <- RunHarmony(
  object = Undefined_obj,           
  group.by.vars = c("sample"),   
  theta = c(2),      # Higher theta for more diverse clustering               
  lambda = c(7),   # Higher lambda to reduce overcorrection                
  sigma = 0.01,           # Lower sigma for tighter clusters
  nclust = 15,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20, 
  early_stop = TRUE,
  dims = 1:40           # More iterations for better convergence
)
Undefined_obj <- RunUMAP(Undefined_obj, 
                         reduction = "harmony", 
                         dims = 1:30, 
                         n.neighbors = 100,
                         n.trees = 500,
                         min.dist = 0.3,
                         learning.rate = 0.15,    # 相对保守的学习率
                         n.epochs = 800,        # 增加迭代次数补偿较小的学习率
                         # spread = 1.2,
                         # repulsion.strength = 1.1,
                         
                         metric = "correlation")
# Find neighbors
Undefined_obj <- FindNeighbors(Undefined_obj, 
                               reduction = "harmony", 
                               dims = 1:30,
                               k.param = 15) 

Undefined_obj <- FindClusters(Undefined_obj,
                              algorithm = 2,
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
for(marker in Marker_B) {
  if(marker %in% rownames(Undefined_obj)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Undefined_obj, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 12. 保存分析结果
# saveRDS(Undefined_obj, "Undefined_analyzed.rds")
print("Undefined细胞分析完成!")

# 16. 生成UMAP可视化图：按样本、组织等元数据
# 按样本(sample)分组绘制UMAP图
pdf("Undefined_umap_by_metadata.pdf", width =
      12, height = 10)
p1 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "sample", pt.size = 0.5) + 
  ggtitle("按样本分组")
p2 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "tissue", pt.size = 0.5) + 
  ggtitle("按组织分组")
p3 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5) + 
  ggtitle("聚类结果")

print(p1)
print(p2)
print(p3)
dev.off()
print("Analysis completed!")

# 1. 首先检查 Markers 中是否有重复
print(length(Markers))
print(length(unique(Markers)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Markers)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Markers <- unique(Markers)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot_Undefined.pdf", width = 32, height = 18)
dot_plot <- DotPlot(Undefined_obj, 
                    features = Marker_B, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- dot_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotplot_data_Undefined.csv")
dot_plot
dev.off()
getwd()

# 14. 可视化
pdf("umap_Undefined_obj_anno_integration.pdf", width = 15, height = 10)
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
p1
p2
p3
p4
p5
dev.off()

# 1. 读取注释文件（无表头）
annotations <- read.csv("Undefined_annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- Undefined_obj$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
Undefined_obj$Annotation <- new_idents
Undefined_obj <- SetIdent(Undefined_obj, value = new_idents)

# 4. 设置颜色
cell_colors <- c(
  "B" = "#E41A1C",   
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
DimPlot(Undefined_obj, 
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

# 获取Undefined_obj中细胞的索引
cells_to_update <- rownames(Undefined_obj@meta.data)

# 检查这些细胞是否都存在于主对象中
cells_in_main <- cells_to_update %in% rownames(B_object@meta.data)
valid_cells_to_update <- cells_to_update[cells_in_main]

# 只更新主对象中存在的细胞的Annotation
B_object@meta.data[valid_cells_to_update, "Annotation"] <- 
  Undefined_obj@meta.data[valid_cells_to_update, "Annotation"]

# 验证更新
table(B_object@meta.data$Annotation)

Idents(B_object) <- "Annotation"
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(B_object, 
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

# B_object$Annotation[B_object@meta.data[["Annotation"]] == "Ciliated secretory cell"] <- "Ciliated cell"
rm(Undefined_obj)
gc()

table(B_object$Annotation)

original_annotations <- B_object@meta.data$Annotation
B_object@meta.data$Annotation_clean <- gsub(" ", "_", original_annotations)

# 使用优化版scPairwiseMASCAnalysis函数
results <- scPairwiseMASCAnalysis(
  seurat_obj = B_object,
  cell_type_col = "Annotation_clean",
  sample_col = "sample",
  contrast_col = "tissue",
  # 在优化版中，此条件正确排除不是'scraping'方法采样的细胞
  exclude_filter = "tissue_sampling_method == 'scraping'",
  fixed_effects_cols = NULL,
  output_dir = "pairwise_MASC_results",
  min_samples = 2,
  # 以下是可选参数，使用默认值
  # p_threshold = 0.05,      # p值显著性阈值
  # fdr_threshold = 0.1,     # FDR显著性阈值
  # min_cells = 10,          # 每种细胞类型的最小细胞数
  # min_prop = 0.001,        # 考虑分析的最小细胞比例
  # save_models = FALSE      # 是否保存MASC模型
)


