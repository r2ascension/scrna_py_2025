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
library(dplyr)
library(scatterplot3d)
library(plotrix)
library(ggsci)
options(future.globals.maxSize = Inf )
library(parallel)
library(ROGUE)
library(RColorBrewer)
library(ggunchull)
source("E:/R/scMASC.R")
source("E:/R/addgrids3d.r")
Markers <- c('FXYD3', 'EPCAM', 'ELF3', 'SERPINF1', 'TSPAN1',
             'SCGB1A1', 'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8',#Epithelial
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

Marker_Epithelial <- c('TP63','KRT5',
                       'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1',
                       'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24',
                       "KRT19","NOTCH3","KRT16","KRT23",
                       'MUC5AC','SPDEF','LYPD2','ITLN1',
                       'ASCL1','GRP',
                       "KRT8",
                       'POU2F3','ASCL2','CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB',"ASCL3","PDE1C",
                       'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
                       'SFTPC','LAMP3','MF5D2A','C8orf4','C11orf96',
                       'VIM','SOX9',
                       'KRT14','MYH11','ACTA2',"MYLK",
                       'DMBT1','RNASE1',
                       'MUC5B','SPDEF',
                       'LYZ','LTF',"PIP","CCL28",
                       "DEUP1","FOXN4","CDC20B","CCNO",
                       'SFTPB','SCGB3A2','SFTA2',
                       'MKI67', 'TOP2A', 'TK1', 'CENPW')

# 设置上皮细胞分析工作目录
Epiwd <- "E:/R/0301/Epithelial_analysis/"

# 检查并创建目录
if (!dir.exists(Epiwd)) {
  dir.create(Epiwd, recursive = TRUE)
  cat("已创建目录:", Epiwd, "\n")
} else {
  cat("目录已存在:", Epiwd, "\n")
}

# 设置工作路径
setwd(Epiwd)
cat("当前工作目录:", getwd(), "\n")

Epithelial_object <- readRDS("Epithelial_analyzed_20250729.rds")#读取数据
Epithelial_object <- NormalizeData(Epithelial_object)#归一化
Epithelial_object <- FindVariableFeatures(Epithelial_object, selection.method = "vst", nfeatures = 4000)#寻找变异基因
Epithelial_object <- ScaleData(Epithelial_object)#标准化
Epithelial_object <- RunPCA(Epithelial_object, npcs = 40)#PCA
Epithelial_object <- RunHarmony(
    object = Epithelial_object,           
    group.by.vars = c("sample"),   
    theta = c(4.5),
    lambda = c(1),
    sigma = 0.015,
    nclust = 30,
    reduction.use = "pca",
    max_iter = 20,
    cool_down = 10,
    # epsilon_harmony = 0.0001,
    # monitor = "harmony_score",
    # factor = 0.8,
    # tol = 1e-10,
    # min_delta = 0.0001,
    patience = 10,
    # epsilon_cluster = 0.00001,
    early_stop = TRUE,
    dims = 1:40           
)
# # 获取降维结果和元数据
# embeddings <- Embeddings(Epithelial_object, "harmony")
# metadata <- Epithelial_object@meta.data
# 
# # 过滤掉ann_level_3为NA的细胞
# valid_cells <- !is.na(metadata$ann_level_3)
# embeddings <- embeddings[valid_cells, ]
# metadata <- metadata[valid_cells, ]
# 
# # 检查并处理无效值
# embeddings_valid <- is.finite(rowSums(embeddings))
# if(sum(!embeddings_valid) > 0) {
#   cat(sprintf("Removed %d cells with NA/NaN/Inf values\n", sum(!embeddings_valid)))
#   embeddings <- embeddings[embeddings_valid, ]
#   metadata <- metadata[embeddings_valid, ]
# }
# 
# # 设置采样参数
# n_cells_sample <- min(10000, nrow(embeddings))
# n_iterations <- 10
# k <- 30
# 
# # 多次抽样计算LISI scores
# ilisi_scores <- c()
# clisi_scores <- c()
# 
# for(i in 1:n_iterations) {
#   set.seed(i)
#   cells_idx <- sample(nrow(embeddings), n_cells_sample)
#   
#   # 抽取对应的数据
#   sampled_embeddings <- embeddings[cells_idx, ]
#   sampled_metadata <- metadata[cells_idx, ]
#   
#   # 尝试计算LISI scores
#   tryCatch({
#     lisi_res <- lisi::compute_lisi(sampled_embeddings, 
#                                    sampled_metadata, 
#                                    c("sample", "ann_level_3"), 
#                                    k)
#     
#     # 存储结果
#     ilisi_scores[i] <- mean(lisi_res[, "sample"])
#     clisi_scores[i] <- mean(lisi_res[, "ann_level_3"])
#     
#     cat(sprintf("Iteration %d completed: iLISI = %.3f, cLISI = %.3f\n", 
#                 i, ilisi_scores[i], clisi_scores[i]))
#   }, error = function(e) {
#     cat(sprintf("Error in iteration %d: %s\n", i, e$message))
#   })
# }
# 
# # 计算并输出结果
# if(length(ilisi_scores) > 0) {
#   cat("\nIntegration Quality Metrics:")
#   cat("\n--------------------------")
#   cat(sprintf("\nNumber of cells used: %d (after removing NA annotations)", nrow(embeddings)))
#   cat(sprintf("\niLISI Score: %.3f ± %.3f", mean(ilisi_scores), sd(ilisi_scores)))
#   cat(sprintf("\ncLISI Score: %.3f ± %.3f", mean(clisi_scores), sd(clisi_scores)))
#   cat("\n\nInterpretation:")
#   cat("\n- iLISI: Higher values (closer to N_batches) indicate better batch mixing")
#   cat("\n- cLISI: Lower values (closer to 1) indicate better cell type separation")
#   cat("\n--------------------------\n")
# } else {
#   cat("\nWarning: Could not compute LISI scores due to errors in all iterations\n")
# }
# Epithelial_object <- process_doublets(Epithelial_object)
# # 创建输出目录（如果不存在）
# if(!dir.exists("qc_plots")) {
#   dir.create("qc_plots")
# }

# # 保存所有knee plots到一个PDF文件
# pdf("qc_plots/knee_plots_analysis.pdf", width = 12, height = 10)

# # 1. UMI knee plot
# counts <- GetAssayData(Epithelial_object, layer = "counts")
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
# pct_var <- Epithelial_object[["pca"]]@stdev / sum(Epithelial_object[["pca"]]@stdev) * 100
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
# n_sample <- min(10000, ncol(Epithelial_object))
# sample_idx <- sample(seq_len(ncol(Epithelial_object)), n_sample)

# harmony_scores <- Epithelial_object[["harmony"]]@cell.embeddings[sample_idx, 1:20]
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
Epithelial_object <- RunUMAP(
    Epithelial_object, 
    reduction = "harmony", 
    dims = 1:40, 
    n.neighbors = 20,
    min.dist = 0.3,
    n_trees = 500,
    # learning.rate = 0.2,    
    # n.epochs = 1400,        
    # spread = 1.2,
    # repulsion.strength = 1.1,
    metric = "correlation"
)

# FindNeighbors函数的缩进修正
Epithelial_object <- FindNeighbors(
    Epithelial_object, 
    reduction = "harmony", 
    dims = 1:40,
    n.trees = 30,
    k.param = 20
) 

Epithelial_object <- FindClusters(
    Epithelial_object, 
    resolution = 3, 
    algorithm = 2
    )
# 8. 绘制UMAP聚类图
print("Generating plots...")
Idents(Epithelial_object) <- Epithelial_object$seurat_clusters
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(Epithelial_object, reduction = "umap", label = TRUE)
dev.off()
# saveRDS(Epithelial_object, "Epithelial_analyzed_20250729.rds")
# 9. 保存分析结果
print("Saving analysis results...")

# Epithelial_object <- readRDS("E:/R/0224/final_analyzed_Epithelial_object.rds")


# 9. 标记基因分析

print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for(marker in Marker_Epithelial) {
  if(marker %in% rownames(Epithelial_object)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Epithelial_object, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(Epithelial_object, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(Epithelial_object, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(Epithelial_object, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(Epithelial_object, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(Epithelial_object, reduction = "umap", group.by = "ann_finest_level", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_finest_level")
p6 <- DimPlot(Epithelial_object, reduction = "umap", group.by = "Annotation", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Annotation")
p1
p2
p3
p4
p5
p6
dev.off()

# 1. 首先检查 Markers 中是否有重复
print(length(Markers))
print(length(unique(Markers)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Markers)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Marker_Epithelial <- unique(Marker_Epithelial)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot.pdf", width = 32, height = 18)
dot_plot <- DotPlot(Epithelial_object, 
                    features = Marker_Epithelial, 
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

tissues <- unique(Epithelial_object$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(Epithelial_object) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(Epithelial_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Epithelial_object, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(Epithelial_object$sample)
pdf("sample_umap_split.pdf", width = 12, height = 12)
Idents(Epithelial_object) <- "sample"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(Epithelial_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Epithelial_object, 
                                                    idents = sample_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(Epithelial_object$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(Epithelial_object) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(Epithelial_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Epithelial_object, 
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
Epithelial_object$Annotation_2 <- Epithelial_object$Annotation
Epithelial_object$ann_finest_level_no_na[is.na(Epithelial_object$ann_finest_level_no_na)] <- "Unknown"
Idents(Epithelial_object) <- "ann_finest_level_no_na"

pdf("ann_level_3_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_finest_level_no_na <- unique(Epithelial_object$ann_finest_level_no_na)
ann_finest_level_no_na <- ann_finest_level_no_na[!is.na(ann_finest_level_no_na)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_finest_level_no_na) {
  print(DimPlot(Epithelial_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(Epithelial_object, 
                                                    idents = anno),
                cols = "grey",
                pt.size = 0.5,
                raster = TRUE,
                label = FALSE) +
          ggtitle(anno) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

# saveRDS(Epithelial_object, "Epithelial_analyzed_20250916.rds")

# 1. 读取注释文件（无表头）
annotations <- read.csv("Annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- Epithelial_object$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
Epithelial_object$Annotation <- new_idents
Epithelial_object <- SetIdent(Epithelial_object, value = new_idents)

# 呼吸道细胞类型配色方案
cell_colors <- c(
  "AT0" = "#E31A1C",          # 鲜红色 - 肺泡上皮祖细胞
  "AT1" = "#FF7F00",          # 橙色 - I型肺泡上皮细胞
  "AT2" = "#1F78B4",          # 蓝色 - II型肺泡上皮细胞
  "Basal_cell" = "#33A02C",   # 绿色 - 基底细胞
  "Ciliated_cell" = "#6A3D9A", # 紫色 - 纤毛细胞
  "Secretory_cell" = "#FF6699",    # 粉红色 - 俱乐部细胞(Clara细胞)
  "Goblet_cell" = "#B15928",  # 棕色 - 杯状细胞
  "Ionocyte" = "#FFFF99",     # 浅黄色 - 离子细胞
  "Mucous_cell" = "#A6CEE3",  # 浅蓝色 - 粘液细胞
  "Serous_cell" = "#FDBF6F",   # 浅橙色 - 浆液细胞
  "PNEC" = "#1F1E33",   # 浅橙色 - 浆液细胞
  "Ciliated_secretory_cell" = "#000000"   # 浅橙色 - 浆液细胞
)
Idents(Epithelial_object) <- "Annotation"  # 先设置identity
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(Epithelial_object, 
        reduction = "umap",
        raster = TRUE,
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
cell_types <- unique(Epithelial_object$Annotation)
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 为每个细胞类型单独绘图
for(cell_type in cell_types) {
  # 创建文件名
  
  
  # 绘制该细胞类型的UMAP图
  
  print(DimPlot(Epithelial_object, 
                reduction = "umap",
                cells.highlight = WhichCells(Epithelial_object, idents = cell_type),
                cols.highlight = cell_colors[cell_type],
                cols = "grey",  # 其他细胞显示为灰色
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(cell_type) +
          theme(legend.text = element_text(size = 12)))
  
}
dev.off()
setwd('..')
getwd()
dir.create('./Basal/')
setwd("./Basal/")
# 假设我们选择 cell_type 为 "T_cells" 的细胞
# 1. 提取Basal细胞子集
# 从主Seurat对象中筛选出Annotation标注为"Basal"的细胞
Basal_obj <- subset(Epithelial_object, subset = Annotation == "Basal_cell")

# 2. 标准数据预处理
# 数据标准化，默认使用"LogNormalize"方法
Basal_obj <- NormalizeData(Basal_obj)

# 3. 寻找高变异基因
# 使用方差稳定化转换(VST)方法选择3000个高变基因，用于后续降维分析
Basal_obj <- FindVariableFeatures(Basal_obj, selection.method = "vst", nfeatures = 3000)

# 4. 数据缩放
# 对所有高变基因进行归一化处理，使其均值为0，方差为1
Basal_obj <- ScaleData(Basal_obj, features = VariableFeatures(Basal_obj))

# 5. 主成分分析
# 使用高变基因进行主成分分析，降低数据维度
Basal_obj <- RunPCA(Basal_obj, npcs = 40)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
Basal_obj <- RunHarmony(
  object = Basal_obj,           
  group.by.vars = c("sample"),   
  theta = c(3),      # Higher theta for more diverse clustering               
  lambda = c(1),   # Higher lambda to reduce overcorrection                
  sigma = 0.03,           # Lower sigma for tighter clusters
  nclust = 15,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20, 
  early_stop = TRUE,
  dims = 1:40           # More iterations for better convergence
)
Basal_obj <- RunUMAP(Basal_obj, 
                         reduction = "harmony", 
                         dims = 1:30, 
                         n.neighbors = 30,
                         # n.trees = 500,
                         # min.dist = 0.3,
                         # learning.rate = 0.15,    # 相对保守的学习率
                         # n.epochs = 800,        # 增加迭代次数补偿较小的学习率
                         # spread = 1.2,
                         # repulsion.strength = 1.1,
                         metric = "correlation")
# Find neighbors
Basal_obj <- FindNeighbors(Basal_obj, 
                               reduction = "harmony", 
                               dims = 1:30,
                               k.param = 15) 

Basal_obj <- FindClusters(Basal_obj,
                              algorithm = 2,
                              group.singletons = FALSE,
                              resolution = 3, # 多个分辨率
                              verbose = TRUE)

# 11. 可视化聚类结果
print("生成聚类可视化...")
# 设置默认使用分辨率为2的聚类结果
Idents(Basal_obj) <- Basal_obj$seurat_clusters

# 保存UMAP聚类图
pdf("Basal_umap_clusters.pdf", width = 10, height = 8)
DimPlot(Basal_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Basal细胞聚类") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
dev.off()

print("Generating feature plots...")
pdf("umap_FeaturePlot_Basal.pdf", width = 8, height = 8)
for(marker in Marker_Epithelial) {
  if(marker %in% rownames(Basal_obj)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Basal_obj, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 12. 保存分析结果
# saveRDS(Basal_obj, "Basal_analyzed.rds")
print("Basal细胞分析完成!")

# 16. 生成UMAP可视化图：按样本、组织等元数据
# 按样本(sample)分组绘制UMAP图

pdf("Basal_umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(Basal_obj, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(Basal_obj, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(Basal_obj, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(Basal_obj, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(Basal_obj, reduction = "umap", group.by = "ann_finest_level", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_finest_level")
p6 <- DimPlot(Basal_obj, reduction = "umap", group.by = "Annotation", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Annotation")
p1
p2
p3
p4
p5
p6
dev.off()
print("Analysis completed!")

# 1. 首先检查 Markers 中是否有重复
print(length(Marker_Epithelial))
print(length(unique(Marker_Epithelial)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Marker_Epithelial[duplicated(Marker_Epithelial)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Marker_Epithelial <- unique(Marker_Epithelial)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot_Basal.pdf", width = 32, height = 18)
dot_plot <- DotPlot(Basal_obj, 
                    features = Marker_Epithelial, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- dot_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotplot_data_Basal.csv")
dot_plot
dev.off()
getwd()

# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
Idents(Basal_obj) <- Basal_obj$seurat_clusters
markers <- FindAllMarkers(Basal_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

# table(subset(Basal_obj, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

# 11. 生成热图
print("Generating heatmap...")
top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)


# 1. 先缩放marker基因
marker_genes <- unique(top10_markers$gene)
Basal_obj <- ScaleData(Basal_obj, features = marker_genes)

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
  print(DoHeatmap(Basal_obj, 
                  features = marker_batches[[i]],
                  size = 24) + 
          NoLegend() +
          ggtitle(paste("Marker Genes Batch", i))+
          theme(axis.text.y = element_text(size = 36))
  )
  # dev.off()
}
dev.off()

# 12. 生成点状图 (Dot Plot)
print("Generating dot plots for marker genes...")

# 12.1 为每个cluster选择top marker基因用于点状图展示
# 可以根据需要调整每个cluster显示的基因数量
top_markers_per_cluster <- markers %>% 
  group_by(cluster) %>% 
  top_n(5, wt = avg_log2FC) %>%  # 每个cluster选择top5基因
  arrange(cluster, desc(avg_log2FC))

# 提取用于点状图的基因列表
dotplot_genes <- unique(top_markers_per_cluster$gene)

# 12.2 生成整体点状图
print("Creating comprehensive dot plot...")
pdf("marker_dotplot_comprehensive.pdf", width = 20, height = 12)
print(DotPlot(Basal_obj, 
              features = dotplot_genes,
              cols = c("lightgrey", "red"),
              dot.scale = 8) + 
        RotatedAxis() +
        ggtitle("Top Marker Genes Expression Across Clusters") +
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              plot.title = element_text(size = 16, hjust = 0.5),
              legend.title = element_text(size = 12)))
dev.off()

# 12.3 按cluster分组生成点状图（如果cluster数量很多）
# 计算cluster数量
n_clusters <- length(unique(Basal_obj$seurat_clusters))

if(n_clusters > 20) {  # 如果cluster数量超过20，按批次展示
  print("Large number of clusters detected. Generating batch dot plots...")
  
  # 将cluster分组
  clusters <- sort(as.numeric(unique(Basal_obj$seurat_clusters)))
  cluster_batches <- split(clusters, ceiling(seq_along(clusters)/10))  # 每批10个cluster
  
  # 为每批cluster生成点状图
  for(i in seq_along(cluster_batches)) {
    batch_clusters <- cluster_batches[[i]]
    
    # 筛选该批次cluster的marker基因
    batch_markers <- markers %>% 
      filter(cluster %in% batch_clusters) %>%
      group_by(cluster) %>% 
      top_n(5, wt = avg_log2FC) %>%
      arrange(cluster, desc(avg_log2FC))
    
    batch_genes <- unique(batch_markers$gene)
    
    # 生成该批次的点状图
    filename <- paste0("marker_dotplot_batch_", i, ".pdf")
    pdf(filename, width = 16, height = 10)
    
    print(DotPlot(Basal_obj, 
                  features = batch_genes,
                  idents = batch_clusters,
                  cols = c("lightgrey", "red"),
                  dot.scale = 8) + 
            RotatedAxis() +
            ggtitle(paste("Marker Genes - Clusters", min(batch_clusters), "to", max(batch_clusters))) +
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 12),
                  plot.title = element_text(size = 14, hjust = 0.5)))
    
    dev.off()
    print(paste("Completed batch", i, "- Clusters:", paste(batch_clusters, collapse = ", ")))
  }
}

# 12.4 生成功能导向的点状图（按基因功能分类）
print("Creating functional category dot plots...")

# 定义一些常见的功能基因类别（可根据研究领域调整）
functional_categories <- list(
  "Cell_Cycle" = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNB2"),
  "Stem_Cell" = c("SOX2", "NANOG", "POU5F1", "KLF4", "MYC"),
  "Differentiation" = c("TP63", "KRT5", "KRT14", "KRT1", "KRT10"),
  "Apoptosis" = c("TP53", "BAX", "BCL2", "CASP3", "CASP7"),
  "Immune_Response" = c("CD3D", "CD4", "CD8A", "CD68", "PTPRC")
)

# 为每个功能类别生成点状图
for(category in names(functional_categories)) {
  genes_in_category <- functional_categories[[category]]
  # 检查哪些基因在数据中存在
  available_genes <- genes_in_category[genes_in_category %in% rownames(Basal_obj)]
  
  if(length(available_genes) > 0) {
    filename <- paste0("dotplot_", category, ".pdf")
    pdf(filename, width = 12, height = 8)
    
    print(DotPlot(Basal_obj, 
                  features = available_genes,
                  cols = c("lightgrey", "red"),
                  dot.scale = 10) + 
            RotatedAxis() +
            ggtitle(paste(category, "Related Genes")) +
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 12),
                  plot.title = element_text(size = 14, hjust = 0.5)))
    
    dev.off()
    print(paste("Generated dot plot for", category, "with", length(available_genes), "genes"))
  }
}

# 12.5 生成split点状图（如果有分组信息，如treatment等）
# 检查tissue分组信息
if("tissue" %in% colnames(Basal_obj@meta.data) && 
   length(unique(Basal_obj$tissue)) > 1) {
  
  print("Creating split dot plot by tissue...")
  
  # 获取分组数量
  n_groups <- length(unique(Basal_obj$tissue))
  print(paste("Number of tissue groups:", n_groups))
  print(paste("Tissue groups:", paste(unique(Basal_obj$tissue), collapse = ", ")))
  
  # 选择top marker基因（减少基因数量以提高可读性）
  top_split_markers <- markers %>% 
    group_by(cluster) %>% 
    top_n(2, wt = avg_log2FC) %>%  # 减少到每个cluster top2基因
    arrange(cluster, desc(avg_log2FC))
  
  split_genes <- unique(top_split_markers$gene)
  print(paste("Number of genes for split plot:", length(split_genes)))
  
  # 生成足够的颜色
  library(RColorBrewer)
  # 根据分组数量生成颜色，每个组需要两种颜色（低表达和高表达）
  colors_needed <- n_groups * 2
  
  # 使用多个调色板组合生成足够的颜色
  if(colors_needed <= 12) {
    colors <- brewer.pal(max(3, colors_needed), "Set3")[1:colors_needed]
  } else {
    # 如果需要更多颜色，组合多个调色板
    colors1 <- brewer.pal(12, "Set3")
    colors2 <- brewer.pal(min(8, colors_needed-12), "Pastel1")
    colors <- c(colors1, colors2)[1:colors_needed]
  }
  
  # 为split dot plot创建颜色渐变
  split_colors <- rep(c("lightgrey", "red"), n_groups)
  
  pdf("marker_dotplot_split_by_tissue.pdf", width = 28, height = 14)
  tryCatch({
    print(DotPlot(Basal_obj, 
                  features = split_genes,
                  cols = split_colors,
                  dot.scale = 6,
                  split.by = "tissue") + 
            RotatedAxis() +
            ggtitle("Top Marker Genes Expression Split by Tissue") +
            theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  plot.title = element_text(size = 16, hjust = 0.5),
                  legend.position = "right"))
  }, error = function(e) {
    print(paste("Error in split dot plot:", e$message))
    # 如果仍然出错，尝试不使用split
    print(DotPlot(Basal_obj, 
                  features = split_genes[1:min(20, length(split_genes))],  # 限制基因数量
                  cols = c("lightgrey", "red"),
                  dot.scale = 6) + 
            RotatedAxis() +
            ggtitle("Top Marker Genes Expression (Simplified)") +
            theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  plot.title = element_text(size = 16, hjust = 0.5)))
  })
  dev.off()
}
# 12.6 保存点状图相关的marker基因信息
print("Saving dot plot gene information...")
write.csv(top_markers_per_cluster, "dotplot_marker_genes.csv", row.names = FALSE)

# 生成基因表达统计摘要 - 改进版本
print("Generating gene expression statistics...")

# 首先检查哪些基因实际存在于对象中
available_genes <- dotplot_genes[dotplot_genes %in% rownames(Basal_obj)]
missing_genes <- dotplot_genes[!dotplot_genes %in% rownames(Basal_obj)]

if(length(missing_genes) > 0) {
  print(paste("Warning: The following genes are not found in the dataset:", 
              paste(missing_genes, collapse = ", ")))
}

print(paste("Processing", length(available_genes), "available genes out of", 
            length(dotplot_genes), "total genes"))

# 使用更安全的方法生成统计信息
gene_stats_list <- list()

for(i in seq_along(available_genes)) {
  gene <- available_genes[i]
  
  if(i %% 10 == 0) {
    print(paste("Processing gene", i, "of", length(available_genes), ":", gene))
  }
  
  tryCatch({
    # 使用AggregateExpression而不是AverageExpression（按Seurat v5建议）
    cluster_expr <- AggregateExpression(Basal_obj, 
                                        features = gene, 
                                        group.by = "seurat_clusters",
                                        return.seurat = FALSE)
    
    # 检查返回的数据结构
    if("RNA" %in% names(cluster_expr) && gene %in% rownames(cluster_expr$RNA)) {
      expr_values <- cluster_expr$RNA[gene, ]
      max_cluster <- names(which.max(expr_values))
      max_expression <- max(expr_values, na.rm = TRUE)
      
      gene_stats_list[[i]] <- data.frame(
        Gene = gene,
        Max_Cluster = max_cluster,
        Max_Expression = max_expression,
        Mean_Expression = mean(expr_values, na.rm = TRUE),
        Median_Expression = median(expr_values, na.rm = TRUE),
        SD_Expression = sd(expr_values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    } else {
      # 如果AggregateExpression失败，尝试备用方法
      print(paste("AggregateExpression failed for gene:", gene, ". Trying alternative method..."))
      
      # 使用手动计算
      Idents(Basal_obj) <- Basal_obj$seurat_clusters
      expr_by_cluster <- sapply(levels(Idents(Basal_obj)), function(cluster) {
        cells_in_cluster <- WhichCells(Basal_obj, idents = cluster)
        if(length(cells_in_cluster) > 0) {
          mean(GetAssayData(Basal_obj, assay = "RNA", slot = "data")[gene, cells_in_cluster], na.rm = TRUE)
        } else {
          0
        }
      })
      
      max_cluster <- names(which.max(expr_by_cluster))
      max_expression <- max(expr_by_cluster, na.rm = TRUE)
      
      gene_stats_list[[i]] <- data.frame(
        Gene = gene,
        Max_Cluster = max_cluster,
        Max_Expression = max_expression,
        Mean_Expression = mean(expr_by_cluster, na.rm = TRUE),
        Median_Expression = median(expr_by_cluster, na.rm = TRUE),
        SD_Expression = sd(expr_by_cluster, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    print(paste("Error processing gene", gene, ":", e$message))
    gene_stats_list[[i]] <<- data.frame(
      Gene = gene,
      Max_Cluster = "Error",
      Max_Expression = NA,
      Mean_Expression = NA,
      Median_Expression = NA,
      SD_Expression = NA,
      stringsAsFactors = FALSE
    )
  })
}

# 合并结果
gene_stats <- do.call(rbind, gene_stats_list[!sapply(gene_stats_list, is.null)])

# 如果仍然为空，创建简化版本
if(nrow(gene_stats) == 0) {
  print("Creating simplified gene statistics...")
  gene_stats <- data.frame(
    Gene = available_genes,
    Max_Cluster = "Not_calculated",
    Max_Expression = NA,
    stringsAsFactors = FALSE
  )
}

# 保存统计结果
write.csv(gene_stats, "dotplot_gene_statistics.csv", row.names = FALSE)

# 生成额外的统计信息
print("Generating additional statistics...")

# 每个cluster的marker基因数量统计
cluster_marker_counts <- markers %>%
  group_by(cluster) %>%
  summarise(
    Total_markers = n(),
    Top5_markers = min(5, n()),
    Avg_logFC = mean(avg_log2FC, na.rm = TRUE),
    Max_logFC = max(avg_log2FC, na.rm = TRUE),
    .groups = 'drop'
  )

write.csv(cluster_marker_counts, "cluster_marker_summary.csv", row.names = FALSE)

# 基因表达分布统计
if(nrow(gene_stats) > 0 && !all(is.na(gene_stats$Max_Expression))) {
  expression_summary <- data.frame(
    Statistic = c("Total_genes", "Mean_max_expression", "Median_max_expression", 
                  "SD_max_expression", "Min_max_expression", "Max_max_expression"),
    Value = c(
      nrow(gene_stats),
      mean(gene_stats$Max_Expression, na.rm = TRUE),
      median(gene_stats$Max_Expression, na.rm = TRUE),
      sd(gene_stats$Max_Expression, na.rm = TRUE),
      min(gene_stats$Max_Expression, na.rm = TRUE),
      max(gene_stats$Max_Expression, na.rm = TRUE)
    )
  )
  
  write.csv(expression_summary, "expression_distribution_summary.csv", row.names = FALSE)
}

print("Dot plot generation completed!")
print(paste("Generated dot plots for", length(available_genes), "available marker genes"))
print(paste("Total clusters analyzed:", length(unique(Basal_obj$seurat_clusters))))
print("Output files created:")
print("- dotplot_marker_genes.csv")
print("- dotplot_gene_statistics.csv")
print("- cluster_marker_summary.csv")
if(exists("expression_summary")) {
  print("- expression_distribution_summary.csv")
}

# 1. 读取注释文件（无表头）
annotations <- read.csv("Basal_annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- Basal_obj$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
Basal_obj$Annotation <- new_idents
Basal_obj <- SetIdent(Basal_obj, value = new_idents)

# 呼吸道细胞类型配色方案
cell_colors <- c(
  "AT0" = "#E31A1C",          # 鲜红色 - 肺泡上皮祖细胞
  "AT1" = "#FF7F00",          # 橙色 - I型肺泡上皮细胞
  "AT2" = "#1F78B4",          # 蓝色 - II型肺泡上皮细胞
  "Basal_cell" = "#33A02C",   # 绿色 - 基底细胞
  "Ciliated_cell" = "#6A3D9A", # 紫色 - 纤毛细胞
  "Secretory_cell" = "#FF6699",    # 粉红色 - 俱乐部细胞(Clara细胞)
  "Goblet_cell" = "#B15928",  # 棕色 - 杯状细胞
  "Ionocyte" = "#FFFF99",     # 浅黄色 - 离子细胞
  "Basal_cell" = "#A6CEE3",  # 浅蓝色 - 粘液细胞
  "Serous_cell" = "#FDBF6F",   # 浅橙色 - 浆液细胞
  "PNEC" = "#1F1E33",   # 浅橙色 - 浆液细胞
  "Basal_cell" = "#000000"   # 浅橙色 - 浆液细胞
)
Idents(Epithelial_object) <- "Annotation"  # 先设置identity
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(Epithelial_object, 
        reduction = "umap",
        raster = TRUE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        cols = cell_colors
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 获取Basal_obj中细胞的索引
cells_to_update <- rownames(Basal_obj@meta.data)

# 检查这些细胞是否都存在于主对象中
cells_in_main <- cells_to_update %in% rownames(Epithelial_object@meta.data)
valid_cells_to_update <- cells_to_update[cells_in_main]

# 只更新主对象中存在的细胞的Annotation
Epithelial_object@meta.data[valid_cells_to_update, "Annotation"] <- 
  Basal_obj@meta.data[valid_cells_to_update, "Annotation"]

# 验证更新
table(Epithelial_object@meta.data$Annotation)

Idents(Epithelial_object) <- "Annotation"
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(Epithelial_object, 
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

Epithelial_object$Annotation[Epithelial_object@meta.data[["Annotation"]] %in% c("Secretory_ciliated_cell")] <- "Secretory_cell"
# rm(Basal_obj)
# gc()

Epithelial_object$Annotation[Epithelial_object@meta.data[["tissue"]] %in% c("sinus",'nose') & Epithelial_object@meta.data[["Annotation"]] %in% c("AT0","AT1","AT2")] <- "Secretory_cell"
# 
table(Epithelial_object$Annotation[Epithelial_object@meta.data[["tissue"]] %in% c('respiratory airway')])

original_annotations <- Epithelial_object@meta.data$Annotation
Epithelial_object@meta.data$Annotation_clean <- gsub(" ", "_", original_annotations)

# 使用优化版scPairwiseMASCAnalysis函数
results <- scPairwiseMASCAnalysis(
  seurat_obj = Epithelial_object,
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  # 在优化版中，此条件正确排除不是'scraping'方法采样的细胞
  exclude_filter = "tissue_sampling_method == 'scraping'",
  fixed_effects_cols = NULL,
  output_dir = "pairwise_MASC_results",
  min_samples = 3,
  # 以下是可选参数，使用默认值
  p_threshold = 0.05,      # p值显著性阈值
  fdr_threshold = 0.1,     # FDR显著性阈值
  min_cells = 20,          # 每种细胞类型的最小细胞数
  min_prop = 0.001,        # 考虑分析的最小细胞比例
  save_models = FALSE      # 是否保存MASC模型
)

################################################################################################################
################################################################################################################

# 函数：生成Harmony整合后的高级UMAP可视化(带置信椭圆和坐标轴箭头)
# 适用于经过Harmony整合的Seurat对象
generate_harmony_umap <- function(seurat_obj, 
                                  reduction = "umap", 
                                  group_by = NULL,
                                  cell_type_col = "Annotation_2",
                                  output_dir = "10-Advanced_Visualizations",
                                  point_size = 0.8,
                                  alpha_level = 0.2,
                                  use_ellipse = TRUE,
                                  show_legend = TRUE,
                                  width = 12, 
                                  height = 10) {
  
  # 载入必要的包
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(ggsci) # 用于配色
  
  # 检查seurat对象
  if(!inherits(seurat_obj, "Seurat")) {
    stop("输入必须是Seurat对象")
  }
  
  # 检查是否已运行Harmony和UMAP
  if(!"harmony" %in% names(seurat_obj@reductions)) {
    warning("未检测到harmony降维结果，请确认对象已运行RunHarmony()函数")
  }
  
  if(!reduction %in% names(seurat_obj@reductions)) {
    stop(paste0(reduction, "降维未在Seurat对象中找到"))
  }
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  original_dir <- getwd()
  setwd(output_dir)
  on.exit(setwd(original_dir)) # 确保函数结束时恢复工作目录
  
  message("创建Harmony整合后的高级UMAP可视化...")
  
  # 整理绘图数据
  mydata <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)  # 取出降维数据
  colnames(mydata) <- c("UMAP_1", "UMAP_2") # 确保列名一致
  mydata <- cbind(mydata, as.data.frame(seurat_obj@meta.data))  # 与注释合并
  
  # 检查细胞类型列是否存在
  if(!cell_type_col %in% colnames(mydata)) {
    stop(paste0("细胞类型列 '", cell_type_col, "' 在meta.data中未找到"))
  }
  
  # 重命名细胞类型列以便于处理
  mydata$CellType <- mydata[[cell_type_col]]
  
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
  
  # 设置配色方案 - 扩展调色板以适应更多细胞类型
  n_cell_types <- length(unique(mydata$CellType))
  mycol <- c(pal_d3()(10), pal_aaas()(10), pal_uchicago()(9), pal_jama()(7))
  
  # 如果cell type数量超过预设颜色数量，扩展颜色
  if(n_cell_types > length(mycol)) {
    mycol <- colorRampPalette(mycol)(n_cell_types)
  }
  
  # 获取分组信息
  if(is.null(group_by)) {
    # 自动检测分组变量
    if("group" %in% colnames(mydata)) {
      group_by <- "group"
    } else if("condition" %in% colnames(mydata)) {
      group_by <- "condition"
    } else if("study" %in% colnames(mydata)) {
      group_by <- "study"
    } else if("tissue" %in% colnames(mydata)) {
      group_by <- "tissue"
    } else {
      group_by <- NULL
      message("未找到分组变量，UMAP不会按组分面")
    }
  } else if(!group_by %in% colnames(mydata)) {
    warning(paste0("指定的分组变量 '", group_by, "' 在meta.data中未找到，不会按组分面"))
    group_by <- NULL
  }
  
  # 绘制基本图片
  p_umap <- ggplot(mydata, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = point_size, aes(col = CellType)) +
    scale_color_manual(values = mycol) + # 设置点的颜色
    labs(color = "Cell Type", 
         title = "UMAP",
         subtitle = if(!is.null(group_by)) paste0("按", group_by, "分组") else NULL)
  
  # 添加置信椭圆（如果需要且每个组有足够的点）
  if(use_ellipse) {
    tryCatch({
      p_umap <- p_umap + 
        stat_ellipse(aes(fill = CellType),
                     geom = "polygon",
                     linetype = 2,        # 置信区间的类型
                     linewidth = 0.7,     # 置信区间的粗细
                     alpha = alpha_level) +
        scale_fill_manual(values = mycol)  # 让置信区间与点的颜色相一致
    }, error = function(e) {
      message("无法添加置信椭圆，可能是某些组的点数太少: ", e$message)
    })
  }
  
  # 添加坐标轴箭头
  p_umap <- p_umap +
    geom_line(data = line.x.data,
              aes(x = x, y = y), 
              arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +  # 绘制X轴坐标轴
    geom_line(data = line.y.data,
              aes(x = x, y = y), 
              arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +  # 绘制Y轴坐标轴
    theme_minimal() +
    theme(
      panel.border = element_blank(),  # 隐藏边框
      axis.title = element_blank(),    # 隐藏轴标题
      axis.text = element_blank(),     # 隐藏文本
      axis.ticks = element_blank(),    # 隐藏轴线
      panel.background = element_rect(fill = 'white'),  # 背景色
      plot.background = element_rect(fill = "white"),
      panel.grid = element_blank(),     # 去除网格线
      legend.position = if(show_legend) "right" else "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # 如果有分组变量，按组分面
  if(!is.null(group_by) && length(unique(mydata[[group_by]])) > 1) {
    # 创建分面公式
    facet_formula <- as.formula(paste("~", group_by))
    p_umap <- p_umap + facet_wrap(facet_formula)
    
    # 调整分面后的图片大小
    n_groups <- length(unique(mydata[[group_by]]))
    if(n_groups > 3) {
      width <- min(16, width + (n_groups - 3) * 2)
    }
  }
  
  # 保存高级UMAP图
  umap_filename <- paste0("harmony_umap_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf")
  ggsave(umap_filename, plot = p_umap, width = width, height = height)
  message(paste0("UMAP图已保存为: ", umap_filename))
  
  # 计算并添加细胞比例图
  if(!is.null(group_by) && length(unique(mydata[[group_by]])) > 1) {
    message("创建细胞比例图...")
    
    # 整理细胞比例数据
    cell_props <- as.data.frame(table(mydata$CellType, mydata[[group_by]]))
    colnames(cell_props) <- c("CellType", "Group", "Freq")
    
    # 计算每组内的比例
    cell_props <- cell_props %>%
      group_by(Group) %>%
      mutate(Total = sum(Freq),
             Proportion = Freq / Total,
             Percentage = Proportion * 100) %>%
      ungroup()
    
    # 添加标签文本 (百分比)
    cell_props$Label <- sprintf("%.1f%%", cell_props$Percentage)
    
    # 计算堆叠位置
    cell_props <- cell_props %>%
      group_by(Group) %>%
      arrange(Group, desc(Freq)) %>%
      mutate(pos = cumsum(Proportion) - 0.5 * Proportion) %>%
      ungroup()
    
    # 绘制堆叠柱状图
    p_prop <- ggplot(cell_props, aes(x = Group, y = Proportion, fill = CellType)) +
      geom_bar(stat = "identity", position = "fill", width = 0.7) + 
      scale_fill_manual(values = mycol) +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      labs(title = "各组中细胞类型比例", x = "", y = "百分比") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            panel.grid.major.x = element_blank())
    
    # 添加百分比标签
    if(length(unique(cell_props$CellType)) <= 10) {  # 只在细胞类型数量不多时添加标签
      p_prop <- p_prop + 
        geom_text(aes(y = pos, label = Label), 
                  size = 3, color = "white", fontface = "bold",
                  data = filter(cell_props, Proportion >= 0.05))  # 只标记占比>=5%的类别
    }
    
    # 创建饼图（每个组一个饼图）
    p_pie_list <- list()
    for(g in unique(cell_props$Group)) {
      # 提取该组数据
      g_data <- cell_props %>% filter(Group == g)
      
      # 创建饼图
      p_pie <- ggplot(g_data, aes(x = "", y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = mycol) +
        labs(title = paste0(g, "组细胞类型比例"), x = NULL, y = NULL) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
      
      # 添加标签（只对主要细胞类型添加）
      if(length(unique(g_data$CellType)) <= 10) {
        p_pie <- p_pie + 
          geom_text(aes(y = pos, label = Label), 
                    size = 3, color = "white", fontface = "bold",
                    data = filter(g_data, Proportion >= 0.05))
      }
      
      p_pie_list[[g]] <- p_pie
    }
    
    # 保存细胞比例图
    prop_filename <- paste0("harmony_cell_proportions_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf")
    ggsave(prop_filename, plot = p_prop, width = 8, height = 6)
    message(paste0("细胞比例图已保存为: ", prop_filename))
    
    # 保存饼图
    for(g in names(p_pie_list)) {
      pie_filename <- paste0("harmony_pie_", g, "_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf")
      ggsave(pie_filename, plot = p_pie_list[[g]], width = 6, height = 6)
      message(paste0(g, "组饼图已保存为: ", pie_filename))
    }
  }
  
  # 返回可视化对象列表
  result_list <- list(
    umap_plot = p_umap,
    prop_plot = if(exists("p_prop")) p_prop else NULL,
    pie_plots = if(exists("p_pie_list")) p_pie_list else NULL
  )
  
  message("可视化函数执行完成")
  return(result_list)
}

# 使用示例
# library(Seurat)
# library(harmony)
# 
# # 加载数据
# seurat_obj <- readRDS("path/to/your/seuraEpithelial_object.rds")
# 
# # 如果需要先运行Harmony和UMAP
# # seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "group", dims.use = 1:30)
# # seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
# 
# 生成高级可视化
# viz_results <- generate_harmony_umap(
#   seurat_obj = Epithelial_object,
#   group_by = "group",          # 按哪个变量分组
#   cell_type_col = "Annotation_2", # 细胞类型注释列
#   output_dir = "10-Harmony_Visualizations", # 输出目录
#   point_size = 0.8,            # 点大小
#   use_ellipse = TRUE,          # 是否添加置信椭圆
#   width = 12,                  # 图像宽度
#   height = 10                  # 图像高度
# )
# # generate_advanced_umap(Epithelial_object)

################################################################################################################
################################################################################################################
################################################################################################################

# 函数：生成标记基因表达的高级可视化（优化层级结构和错误处理）
generate_marker_viz <- function(seurat_obj, 
                                ident_col = "Annotation_2",
                                top_n = 5, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25,
                                only.pos = TRUE,
                                test.use = "wilcox",
                                assay = "RNA",
                                slot = "data",
                                output_dir = "Marker_Visualization",
                                custom_markers = NULL,
                                add_dotplot = TRUE,
                                add_violin = TRUE,
                                save_tables = TRUE,
                                verbose = TRUE) {
  
  # ============= 函数初始化部分 =============
  # 1. 初始化环境和验证输入
  init_result <- tryCatch({
    if(verbose) message("初始化标记基因可视化...")
    
    # 加载必要的包
    suppressPackageStartupMessages({
      require(Seurat)
      require(dplyr)
      require(ggplot2)
      require(pheatmap)
      require(RColorBrewer)
      require(viridis)
    })
    
    # 验证Seurat对象
    if(!inherits(seurat_obj, "Seurat")) {
      stop("输入必须是Seurat对象")
    }
    
    # 创建输出目录
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    original_dir <- getwd()
    setwd(output_dir)
    
    # 初始化结果列表
    result_list <- list()
    
    # 创建时间戳
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
    
    list(
      original_dir = original_dir,
      result_list = result_list,
      timestamp = timestamp
    )
  }, error = function(e) {
    message("初始化失败: ", e$message)
    return(NULL)
  })
  
  # 如果初始化失败，立即返回
  if(is.null(init_result)) {
    return(NULL)
  }
  
  # 提取初始化结果
  original_dir <- init_result$original_dir
  result_list <- init_result$result_list
  timestamp <- init_result$timestamp
  
  # 确保函数结束时恢复工作目录
  on.exit(setwd(original_dir), add = TRUE)
  
  # ============= 标记基因寻找/处理部分 =============
  # 2. 获取或生成标记基因
  markers_result <- tryCatch({
    # 检查标识列
    if(!ident_col %in% colnames(seurat_obj@meta.data)) {
      stop(paste0("标识列 '", ident_col, "' 在meta.data中未找到"))
    }
    
    # 设置标识
    if(verbose) message("设置标识为: ", ident_col)
    Idents(seurat_obj) <- ident_col
    
    # 处理标记基因数据
    if(!is.null(custom_markers)) {
      # 使用自定义标记基因
      if(verbose) message("使用提供的自定义标记基因列表")
      
      if(is.character(custom_markers)) {
        # 如果是字符向量，构建一个虚拟markers数据框
        cluster_labels <- unique(seurat_obj@meta.data[[ident_col]])
        markers <- data.frame()
        
        # 为每个cluster随机分配一些标记基因
        genes_per_cluster <- ceiling(length(custom_markers) / length(cluster_labels))
        
        for(i in seq_along(cluster_labels)) {
          # 选择一部分基因
          start_idx <- (i-1) * genes_per_cluster + 1
          end_idx <- min(i * genes_per_cluster, length(custom_markers))
          cluster_genes <- custom_markers[start_idx:end_idx]
          
          if(length(cluster_genes) > 0) {
            tmp_df <- data.frame(
              p_val = rep(0.01, length(cluster_genes)),
              avg_log2FC = runif(length(cluster_genes), 1, 3),
              pct.1 = runif(length(cluster_genes), 0.5, 0.9),
              pct.2 = runif(length(cluster_genes), 0.1, 0.4),
              p_val_adj = rep(0.05, length(cluster_genes)),
              cluster = rep(cluster_labels[i], length(cluster_genes)),
              gene = cluster_genes,
              stringsAsFactors = FALSE
            )
            markers <- rbind(markers, tmp_df)
          }
        }
      } else if(is.data.frame(custom_markers)) {
        # 如果是数据框，确保有必要的列
        req_cols <- c("gene", "cluster")
        if(!all(req_cols %in% colnames(custom_markers))) {
          stop("自定义标记基因数据框必须包含gene和cluster列")
        }
        markers <- custom_markers
        
        # 如果没有必要的列，添加虚拟数据
        if(!"p_val" %in% colnames(markers)) markers$p_val <- 0.01
        if(!"avg_log2FC" %in% colnames(markers)) markers$avg_log2FC <- 1
        if(!"pct.1" %in% colnames(markers)) markers$pct.1 <- 0.7
        if(!"pct.2" %in% colnames(markers)) markers$pct.2 <- 0.2
        if(!"p_val_adj" %in% colnames(markers)) markers$p_val_adj <- 0.05
      } else {
        stop("自定义标记基因必须是字符向量或数据框")
      }
    } else {
      # 寻找标记基因
      if(verbose) message("正在寻找标记基因，这可能需要一些时间...")
      
      # 寻找所有标记基因
      markers <- FindAllMarkers(seurat_obj, 
                                only.pos = only.pos, 
                                min.pct = min.pct, 
                                logfc.threshold = logfc.threshold,
                                test.use = test.use)
    }
    
    # 如果没有找到标记基因
    if(nrow(markers) == 0) {
      stop("未找到任何标记基因，请尝试降低阈值或提供自定义标记基因")
    }
    
    # 为每个细胞类型选择top N个标记基因
    top_markers <- markers %>% 
      group_by(cluster) %>% 
      top_n(n = min(n(), top_n), wt = avg_log2FC)
    
    # 保存标记基因结果
    if(save_tables) {
      write.csv(markers, file = paste0("all_markers_", timestamp, ".csv"), row.names = FALSE)
      write.csv(top_markers, file = paste0("top_markers_", timestamp, ".csv"), row.names = FALSE)
      if(verbose) message("已保存标记基因表格")
    }
    
    list(
      markers = markers,
      top_markers = top_markers
    )
  }, error = function(e) {
    message("标记基因处理失败: ", e$message)
    # 尝试使用更简单的方法找标记
    tryCatch({
      if(verbose) message("尝试使用简化方法寻找标记基因...")
      markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
      if(nrow(markers) == 0) {
        return(NULL)
      }
      top_markers <- markers %>% 
        group_by(cluster) %>% 
        top_n(n = min(n(), top_n), wt = avg_log2FC)
      
      return(list(
        markers = markers,
        top_markers = top_markers
      ))
    }, error = function(e2) {
      message("简化标记基因寻找也失败: ", e2$message)
      return(NULL)
    })
  })
  
  # 如果标记基因处理失败，立即返回
  if(is.null(markers_result)) {
    return(NULL)
  }
  
  # 提取标记基因结果
  markers <- markers_result$markers
  top_markers <- markers_result$top_markers
  
  # 将标记基因添加到结果列表
  result_list$markers <- markers
  result_list$top_markers <- top_markers
  
  # ============= 可视化部分 =============
  # 3. 创建散点图
  scatter_result <- tryCatch({
    if(verbose) message("创建标记基因散点图...")
    
    # 准备数据
    viz_markers <- markers
    viz_markers$gene_rank <- rank(-viz_markers$avg_log2FC)
    
    # 高亮顶部基因
    top_genes <- unique(top_markers$gene)
    viz_markers$highlight <- viz_markers$gene %in% top_genes
    
    # 计算每个cluster的点数量以动态调整facet大小
    cluster_counts <- table(viz_markers$cluster)
    cluster_order <- names(sort(cluster_counts, decreasing = TRUE))
    viz_markers$cluster <- factor(viz_markers$cluster, levels = cluster_order)
    
    # 为了美观，将标签限制在每个面板的top_n个
    label_genes <- viz_markers %>%
      group_by(cluster) %>%
      top_n(top_n, wt = avg_log2FC) %>%
      pull(gene)
    
    viz_markers$label <- ifelse(viz_markers$gene %in% label_genes, viz_markers$gene, "")
    
    # 创建散点图
    p_scatter <- ggplot(viz_markers, 
                        aes(x = pct.1 - pct.2, y = avg_log2FC)) +
      geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
      scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red3")) +
      scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) +
      scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1)) +
      facet_wrap(~cluster, scales = "free_y", 
                 ncol = min(4, length(unique(viz_markers$cluster)))) +
      labs(x = "表达细胞百分比差异 (pct.1 - pct.2)",
           y = "平均Log2倍变化",
           title = paste0("各", ident_col, "的标记基因")) +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    # 尝试添加标签
    if(requireNamespace("ggrepel", quietly = TRUE)) {
      p_scatter <- p_scatter + 
        ggrepel::geom_text_repel(data = subset(viz_markers, label != ""),
                                 aes(label = label),
                                 box.padding = 0.5,
                                 point.padding = 0.3,
                                 segment.color = "grey50",
                                 size = 3)
    }
    
    # 保存散点图
    scatter_filename <- paste0("marker_genes_scatter_", timestamp, ".pdf")
    ggsave(scatter_filename, plot = p_scatter, 
           width = min(20, 4 + 2 * length(unique(viz_markers$cluster))), 
           height = min(16, 3 + 1.5 * length(unique(viz_markers$cluster))))
    if(verbose) message(paste0("已保存标记基因散点图: ", scatter_filename))
    
    list(
      p_scatter = p_scatter,
      scatter_filename = scatter_filename
    )
  }, error = function(e) {
    message("散点图创建失败: ", e$message)
    return(NULL)
  })
  
  # 保存散点图结果（即使失败也继续）
  if(!is.null(scatter_result)) {
    result_list$scatter_plot <- scatter_result$p_scatter
    scatter_filename <- scatter_result$scatter_filename
  }
  
  # 4. 创建热图
  heatmap_result <- tryCatch({
    if(verbose) message("创建标记基因热图...")
    
    # 获取所有top标记基因的并集
    top_genes_union <- unique(top_markers$gene)
    
    # 检查基因是否在对象中
    genes_present <- top_genes_union[top_genes_union %in% rownames(seurat_obj)]
    if(length(genes_present) == 0) {
      warning("没有标记基因在Seurat对象中找到，跳过热图绘制")
      return(NULL)
    } else if(length(genes_present) < length(top_genes_union)) {
      warning(paste0("只有", length(genes_present), "/", length(top_genes_union), 
                     "个标记基因在Seurat对象中找到"))
      top_genes_union <- genes_present
    }
    
    # 安全计算表达矩阵
    avg_expr_matrix <- tryCatch({
      # 尝试标准方法
      if(verbose) message("  使用标准方法计算平均表达...")
      avg_expr <- AverageExpression(seurat_obj, 
                                    features = top_genes_union, 
                                    assays = assay,
                                    group.by = ident_col,
                                    slot = slot)
      avg_expr[[assay]]
    }, error = function(e) {
      # 如果标准方法失败，尝试手动计算
      if(verbose) message("  标准方法失败，尝试手动计算平均表达...")
      
      metadata <- seurat_obj@meta.data
      expr_data <- GetAssayData(seurat_obj, assay = assay, slot = slot)
      
      # 仅保留存在的基因
      expr_data <- expr_data[rownames(expr_data) %in% top_genes_union, ]
      
      # 为每个细胞类型计算平均表达
      cell_types <- unique(metadata[[ident_col]])
      heatmap_data <- matrix(0, nrow = nrow(expr_data), ncol = length(cell_types))
      rownames(heatmap_data) <- rownames(expr_data)
      colnames(heatmap_data) <- cell_types
      
      for(ct in cell_types) {
        # 获取该细胞类型的细胞
        cells <- rownames(metadata)[metadata[[ident_col]] == ct]
        # 计算平均表达
        if(length(cells) > 0) {
          # 处理SCT对象可能的Rle类型数据
          cell_expr <- expr_data[, cells, drop = FALSE]
          if(is(cell_expr, "dgCMatrix")) {
            heatmap_data[, ct] <- rowMeans(cell_expr)
          } else {
            # 一行一行处理，避免Rle类型错误
            for(g in rownames(cell_expr)) {
              heatmap_data[g, ct] <- mean(as.numeric(cell_expr[g,]), na.rm = TRUE)
            }
          }
        }
      }
      return(heatmap_data)
    })
    
    # 如果表达矩阵获取失败
    if(is.null(avg_expr_matrix) || nrow(avg_expr_matrix) == 0) {
      warning("无法获取有效的表达矩阵，跳过热图绘制")
      return(NULL)
    }
    
    # 安全计算Z-score
    heatmap_data_scaled <- tryCatch({
      if(verbose) message("  计算Z-score标准化...")
      t(scale(t(avg_expr_matrix)))
    }, error = function(e) {
      if(verbose) message("  标准Z-score计算失败，尝试手动计算...")
      # 手动计算z-scores
      result <- avg_expr_matrix
      for(i in 1:nrow(avg_expr_matrix)) {
        row_mean <- mean(as.numeric(avg_expr_matrix[i,]), na.rm = TRUE)
        row_sd <- sd(as.numeric(avg_expr_matrix[i,]), na.rm = TRUE)
        if(row_sd > 0) {
          result[i,] <- (avg_expr_matrix[i,] - row_mean) / row_sd
        }
      }
      return(result)
    })
    
    # 准备热图注释
    gene_cluster <- markers %>%
      filter(gene %in% rownames(heatmap_data_scaled)) %>%
      select(gene, cluster) %>%
      distinct()
    
    # 如果有重复基因，保留第一个出现的cluster
    gene_cluster <- gene_cluster %>%
      group_by(gene) %>%
      slice(1) %>%
      ungroup()
    
    # 创建注释数据框
    gene_anno <- data.frame(
      Cluster = gene_cluster$cluster[match(rownames(heatmap_data_scaled), gene_cluster$gene)]
    )
    rownames(gene_anno) <- rownames(heatmap_data_scaled)
    
    # 设置颜色方案
    n_clusters <- length(unique(gene_anno$Cluster))
    anno_colors <- list(
      Cluster = setNames(colorRampPalette(brewer.pal(min(9, n_clusters), "Set1"))(n_clusters),
                         unique(gene_anno$Cluster))
    )
    heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    
    # 绘制热图
    heatmap_filename <- paste0("marker_genes_heatmap_", timestamp, ".pdf")
    pdf(heatmap_filename, 
        width = min(14, 6 + 0.2 * ncol(heatmap_data_scaled)), 
        height = min(16, 6 + 0.15 * nrow(heatmap_data_scaled)))
    
    # 处理潜在的排序问题
    heat_labels_row <- if(nrow(heatmap_data_scaled) > 50) FALSE else rownames(heatmap_data_scaled)
    
    # 尝试绘制热图
    tryCatch({
      pheatmap(heatmap_data_scaled,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = heat_labels_row,
               show_colnames = TRUE,
               annotation_row = gene_anno,
               annotation_colors = anno_colors,
               fontsize_row = 7,
               fontsize_col = 9,
               color = heatmap_colors,
               main = paste0(ident_col, "标记基因表达热图"),
               angle_col = 45,
               border_color = NA,
               fontsize = 10)
    }, error = function(e) {
      # 如果聚类失败，尝试无聚类版本
      message("  热图聚类出错，使用无聚类版本...")
      pheatmap(heatmap_data_scaled,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               show_rownames = heat_labels_row,
               show_colnames = TRUE,
               annotation_row = gene_anno,
               annotation_colors = anno_colors,
               fontsize_row = 7,
               fontsize_col = 9,
               color = heatmap_colors,
               main = paste0(ident_col, "标记基因表达热图 (无聚类)"),
               angle_col = 45,
               border_color = NA,
               fontsize = 10)
    })
    
    dev.off()
    if(verbose) message(paste0("已保存标记基因热图: ", heatmap_filename))
    
    list(
      heatmap_data = heatmap_data_scaled,
      heatmap_filename = heatmap_filename
    )
  }, error = function(e) {
    message("热图创建失败: ", e$message)
    return(NULL)
  })
  
  # 保存热图结果
  if(!is.null(heatmap_result)) {
    result_list$heatmap_data <- heatmap_result$heatmap_data
    heatmap_filename <- heatmap_result$heatmap_filename
  }
  
  # 5. 创建点图
  dot_result <- NULL
  if(add_dotplot) {
    dot_result <- tryCatch({
      if(verbose) message("创建标记基因点图...")
      
      # 获取top基因
      top_genes_union <- unique(top_markers$gene)
      
      # 过滤存在的基因
      genes_for_dot <- top_genes_union[top_genes_union %in% rownames(seurat_obj)]
      
      if(length(genes_for_dot) == 0) {
        warning("没有找到可用于点图的基因，跳过此步骤")
        return(NULL)
      }
      
      # 创建点图
      p_dot <- DotPlot(seurat_obj, 
                       features = genes_for_dot, 
                       group.by = ident_col,
                       assay = assay,
                       scale = TRUE) +
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
        labs(title = paste0(ident_col, "标记基因点图"))
      
      # 尝试使用viridis配色
      if(requireNamespace("viridis", quietly = TRUE)) {
        p_dot <- p_dot + scale_color_viridis()
      }
      
      # 保存点图
      dot_filename <- paste0("marker_genes_dotplot_", timestamp, ".pdf")
      ggsave(dot_filename, plot = p_dot, 
             width = min(14, 6 + 0.2 * length(genes_for_dot)), 
             height = min(12, 4 + 0.2 * length(unique(seurat_obj@meta.data[[ident_col]]))))
      if(verbose) message(paste0("已保存标记基因点图: ", dot_filename))
      
      list(
        p_dot = p_dot,
        dot_filename = dot_filename
      )
    }, error = function(e) {
      message("点图创建失败: ", e$message)
      return(NULL)
    })
  }
  
  # 保存点图结果
  if(!is.null(dot_result)) {
    result_list$dot_plot <- dot_result$p_dot
    dot_filename <- dot_result$dot_filename
  }
  
  # 6. 创建小提琴图
  violin_result <- NULL
  if(add_violin) {
    violin_result <- tryCatch({
      if(verbose) message("创建标记基因小提琴图...")
      
      # 为每个cluster选择最佳的几个标记
      violin_genes <- top_markers %>%
        group_by(cluster) %>%
        top_n(n = min(n(), 2), wt = avg_log2FC) %>%
        pull(gene)
      
      violin_genes <- unique(violin_genes)
      violin_genes <- violin_genes[violin_genes %in% rownames(seurat_obj)]
      
      if(length(violin_genes) == 0) {
        warning("没有找到可用于小提琴图的基因，跳过此步骤")
        return(NULL)
      }
      
      # 创建小提琴图
      p_violin <- VlnPlot(seurat_obj, 
                          features = violin_genes, 
                          group.by = ident_col,
                          pt.size = 0,  # 隐藏点
                          assay = assay,
                          combine = FALSE)
      
      # 确保p_violin是一个列表
      if(!is.list(p_violin)) {
        p_violin <- list(p_violin)
      }
      
      # 调整小提琴图样式
      for(i in seq_along(p_violin)) {
        p_violin[[i]] <- p_violin[[i]] + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none",
                plot.title = element_text(face = "bold", size = 10))
      }
      
      # 创建多页PDF文件
      multi_violin_filename <- paste0("marker_genes_violin_all_", timestamp, ".pdf")
      pdf(multi_violin_filename, width = 10, height = 6)
      
      # 逐个打印每个小提琴图到同一个PDF文件
      for(i in seq_along(p_violin)) {
        print(p_violin[[i]])
      }
      
      # 关闭PDF设备
      dev.off()
      if(verbose) message(paste0("已将所有小提琴图保存到: ", multi_violin_filename))
      
      # 尝试创建组合图
      if(length(p_violin) > 1 && requireNamespace("patchwork", quietly = TRUE)) {
        combined_filename <- paste0("marker_genes_violin_combined_", timestamp, ".pdf")
        p_combined <- patchwork::wrap_plots(p_violin, ncol = 2)
        ggsave(combined_filename, plot = p_combined, 
               width = min(14, 4 + length(p_violin) * 1.5), 
               height = min(12, 4 + length(p_violin) * 0.8))
        if(verbose) message(paste0("已保存组合小提琴图: ", combined_filename))
        list(
          p_violin = p_violin,
          p_combined = p_combined,
          multi_violin_filename = multi_violin_filename,
          combined_filename = combined_filename
        )
      } else {
        list(
          p_violin = p_violin,
          multi_violin_filename = multi_violin_filename
        )
      }
    }, error = function(e) {
      message("小提琴图创建失败: ", e$message)
      return(NULL)
    })
  }
  
  # ============= 综合报告部分 =============
  # 7. 创建综合可视化
  summary_result <- tryCatch({
    if(verbose) message("创建综合可视化报告...")
    
    # 收集可用的图表
    plot_list <- list()
    
    # 添加散点图
    if(!is.null(scatter_result)) {
      plot_list$scatter <- scatter_result$p_scatter
    }
    
    # 添加点图
    if(!is.null(dot_result)) {
      plot_list$dot <- dot_result$p_dot
    }
    
    # 如果有多个图表，创建组合图
    if(length(plot_list) > 1 && requireNamespace("patchwork", quietly = TRUE)) {
      p_summary <- patchwork::wrap_plots(
        plot_list,
        ncol = 1,
        heights = c(1.5, rep(1, length(plot_list) - 1))
      )
      
      # 保存组合图
      summary_filename <- paste0("marker_genes_summary_", timestamp, ".pdf")
      ggsave(summary_filename, plot = p_summary, width = 12, height = 4 * length(plot_list))
      if(verbose) message(paste0("已保存综合可视化报告: ", summary_filename))
      
      list(
        p_summary = p_summary,
        summary_filename = summary_filename
      )
    } else if(length(plot_list) == 1) {
      # 如果只有一个图表，直接使用
      p_summary <- plot_list[[1]]
      summary_filename <- paste0("marker_genes_summary_", timestamp, ".pdf")
      ggsave(summary_filename, plot = p_summary, width = 12, height = 10)
      if(verbose) message(paste0("已保存单一可视化: ", summary_filename))
      
      list(
        p_summary = p_summary,
        summary_filename = summary_filename
      )
    } else {
      if(verbose) message("没有足够的图表用于创建综合报告")
      NULL
    }
  }, error = function(e) {
    message("综合报告创建失败: ", e$message)
    return(NULL)
  })
  
  # 保存综合报告结果
  if(!is.null(summary_result)) {
    result_list$summary_plot <- summary_result$p_summary
    summary_filename <- summary_result$summary_filename
  }
  
  # 8. 创建HTML报告
  html_result <- tryCatch({
    if(verbose) message("创建HTML报告...")
    
    # 初始化HTML内容
    html_content <- paste0(
      "<!DOCTYPE html>
      <html>
      <head>
        <title>标记基因分析报告</title>
        <style>
          body { font-family: Arial, sans-serif; margin: 20px; }
          h1, h2 { color: #333366; }
          .figure { margin: 20px 0; text-align: center; }
          .figure img { max-width: 100%; border: 1px solid #ddd; }
          table { border-collapse: collapse; width: 100%; }
          th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
          th { background-color: #f2f2f2; }
          tr:nth-child(even) { background-color: #f9f9f9; }
          .note { color: #666; font-style: italic; }
        </style>
      </head>
      <body>
        <h1>单细胞标记基因分析报告</h1>
        <p>分析时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>
        <p>标识变量: ", ident_col, "</p>
        <p>分析参数: min.pct = ", min.pct, ", logfc.threshold = ", logfc.threshold, "</p>
        
        <h2>1. 标记基因总览</h2>
        <p>共找到 ", nrow(markers), " 个标记基因，每个细胞类型展示 ", top_n, " 个表达最高的基因。</p>")
    
    # 添加图表部分
    figure_count <- 1
    
    # 添加散点图
    if(!is.null(scatter_result) && file.exists(scatter_result$scatter_filename)) {
      html_content <- paste0(html_content, "
        <div class='figure'>
          <img src='", scatter_result$scatter_filename, "' alt='标记基因散点图'>
          <p>图", figure_count, ": 标记基因散点图。红色点表示每个类别的顶部标记基因。</p>
        </div>")
      figure_count <- figure_count + 1
    }
    
    # 添加热图
    if(!is.null(heatmap_result) && file.exists(heatmap_result$heatmap_filename)) {
      html_content <- paste0(html_content, "
        <div class='figure'>
          <img src='", heatmap_result$heatmap_filename, "' alt='标记基因热图'>
          <p>图", figure_count, ": 标记基因表达热图。展示了不同细胞类型中标记基因的表达水平。</p>
        </div>")
      figure_count <- figure_count + 1
    }
    
    # 添加点图
    if(!is.null(dot_result) && file.exists(dot_result$dot_filename)) {
      html_content <- paste0(html_content, "
        <div class='figure'>
          <img src='", dot_result$dot_filename, "' alt='标记基因点图'>
          <p>图", figure_count, ": 标记基因点图。点的大小表示表达该基因的细胞比例，颜色表示平均表达水平。</p>
        </div>")
      figure_count <- figure_count + 1
    }
    
    # 添加小提琴图
    if(!is.null(violin_result)) {
      # 优先添加组合图
      if(!is.null(violin_result$combined_filename) && file.exists(violin_result$combined_filename)) {
        html_content <- paste0(html_content, "
          <div class='figure'>
            <img src='", violin_result$combined_filename, "' alt='标记基因小提琴图'>
            <p>图", figure_count, ": 标记基因小提琴图。展示了标记基因在不同细胞类型中的表达分布。</p>
          </div>")
        figure_count <- figure_count + 1
      } else if(length(violin_result$violin_filenames) > 0) {
        # 添加单独的小提琴图
        html_content <- paste0(html_content, "
          <h3>标记基因小提琴图</h3>")
        
        for(i in seq_along(violin_result$violin_filenames)) {
          if(file.exists(violin_result$violin_filenames[i])) {
            gene_name <- gsub("marker_genes_violin_|_.*\\.pdf", "", violin_result$violin_filenames[i])
            html_content <- paste0(html_content, "
              <div class='figure'>
                <img src='", violin_result$violin_filenames[i], "' alt='标记基因小提琴图'>
                <p>图", figure_count, ": 基因 ", gene_name, " 的表达分布。</p>
              </div>")
            figure_count <- figure_count + 1
          }
        }
      }
    }
    
    # 添加标记基因表格
    html_content <- paste0(html_content, "
      <h2>2. 每个细胞类型的顶部标记基因</h2>
      <table>
        <tr>
          <th>细胞类型</th>
          <th>基因</th>
          <th>平均Log2倍变化</th>
          <th>调整后p值</th>
          <th>表达比例差异</th>
        </tr>")
    
    # 填充表格内容
    for(cl in unique(top_markers$cluster)) {
      genes <- top_markers %>% filter(cluster == cl)
      
      if(nrow(genes) > 0) {
        html_content <- paste0(html_content, "<tr><td rowspan='", nrow(genes), "'>", cl, "</td>")
        
        # 添加第一个基因
        html_content <- paste0(html_content, 
                               "<td>", genes$gene[1], "</td>",
                               "<td>", round(genes$avg_log2FC[1], 2), "</td>",
                               "<td>", format(genes$p_val_adj[1], scientific = TRUE, digits = 2), "</td>",
                               "<td>", round(genes$pct.1[1] - genes$pct.2[1], 2), "</td></tr>")
        
        # 添加剩余基因
        if(nrow(genes) > 1) {
          for(i in 2:nrow(genes)) {
            html_content <- paste0(html_content, "<tr>",
                                   "<td>", genes$gene[i], "</td>",
                                   "<td>", round(genes$avg_log2FC[i], 2), "</td>",
                                   "<td>", format(genes$p_val_adj[i], scientific = TRUE, digits = 2), "</td>",
                                   "<td>", round(genes$pct.1[i] - genes$pct.2[i], 2), "</td></tr>")
          }
        }
      }
    }
    
    html_content <- paste0(html_content, "
      </table>")
    
    # 添加下载部分
    html_content <- paste0(html_content, "
      <h2>3. 数据下载</h2>
      <p>完整的标记基因列表和顶部标记基因可以在以下文件中找到：</p>
      <ul>")
    
    if(save_tables) {
      all_markers_file <- paste0("all_markers_", timestamp, ".csv")
      top_markers_file <- paste0("top_markers_", timestamp, ".csv")
      
      if(file.exists(all_markers_file)) {
        html_content <- paste0(html_content, "
          <li><a href='", all_markers_file, "'>所有标记基因</a></li>")
      }
      
      if(file.exists(top_markers_file)) {
        html_content <- paste0(html_content, "
          <li><a href='", top_markers_file, "'>顶部标记基因</a></li>")
      }
    }
    
    html_content <- paste0(html_content, "
      </ul>
      
      <h2>4. 分析方法</h2>
      <p>本分析使用Seurat的FindAllMarkers函数识别每个细胞类型的特异表达基因。标记基因的识别基于以下设置:
      <ul>
        <li>统计检验方法: ", test.use, "</li>
        <li>最小表达细胞比例(min.pct): ", min.pct, "</li>
        <li>最小Log2倍变化阈值: ", logfc.threshold, "</li>
        <li>仅选择上调基因: ", ifelse(only.pos, "是", "否"), "</li>
      </ul>
      </p>
      
      <p class='note'>报告生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>
    </body>
    </html>")
    
    # 保存HTML报告
    html_filename <- paste0("marker_genes_report_", timestamp, ".html")
    writeLines(html_content, html_filename)
    if(verbose) message(paste0("已保存HTML报告: ", html_filename))
    
    list(
      html_content = html_content,
      html_filename = html_filename
    )
  }, error = function(e) {
    message("HTML报告创建失败: ", e$message)
    return(NULL)
  })
  
  # 保存HTML报告结果
  if(!is.null(html_result)) {
    result_list$html_report_file <- html_result$html_filename
  }
  
  # 返回所有结果
  if(verbose) message("标记基因分析完成")
  return(result_list)
}

# 简化版标记基因可视化函数（适用于基本场景或主函数失败时）
simple_marker_viz <- function(seurat_obj, ident_col = "Annotation_2", top_n = 5, output_dir = "Simple_Markers") {
  # 加载必要的包
  suppressPackageStartupMessages({
    require(Seurat)
    require(dplyr)
    require(ggplot2)
  })
  
  message("使用简化版标记基因可视化函数...")
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  original_dir <- getwd()
  setwd(output_dir)
  on.exit(setwd(original_dir))
  
  # 设置分组标识
  if(ident_col %in% colnames(seurat_obj@meta.data)) {
    Idents(seurat_obj) <- ident_col
  } else {
    message("警告: 找不到指定的标识列 ", ident_col, "，使用默认标识")
  }
  
  # 寻找标记基因
  message("寻找标记基因...")
  tryCatch({
    markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    # 为每个类型选择top基因
    top_markers <- markers %>% 
      group_by(cluster) %>% 
      top_n(n = top_n, wt = avg_log2FC)
    
    # 保存标记基因列表
    write.csv(markers, "all_markers_simple.csv", row.names = FALSE)
    write.csv(top_markers, "top_markers_simple.csv", row.names = FALSE)
    
    # 绘制一系列简单可视化
    result_plots <- list()
    
    # 1. 热图
    message("创建热图...")
    tryCatch({
      top_genes <- unique(top_markers$gene)
      top_genes <- top_genes[top_genes %in% rownames(seurat_obj)]
      if(length(top_genes) > 0) {
        p_heatmap <- DoHeatmap(seurat_obj, features = top_genes) + 
          theme(axis.text.y = element_text(size = 8))
        ggsave("simple_heatmap.pdf", plot = p_heatmap, width = 10, height = 12)
        result_plots$heatmap <- p_heatmap
      }
    }, error = function(e) {
      message("热图创建失败: ", e$message)
    })
    
    # 2. 点图
    message("创建点图...")
    tryCatch({
      if(length(top_genes) > 0) {
        p_dot <- DotPlot(seurat_obj, features = top_genes) + 
          coord_flip() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave("simple_dotplot.pdf", plot = p_dot, width = 10, height = 12)
        result_plots$dotplot <- p_dot
      }
    }, error = function(e) {
      message("点图创建失败: ", e$message)
    })
    
    # 3. 小提琴图（仅前几个基因）
    message("创建小提琴图...")
    tryCatch({
      if(length(top_genes) > 0) {
        # 每个类别选择一个top基因
        top_violin_genes <- top_markers %>%
          group_by(cluster) %>%
          slice(1) %>%
          pull(gene)
        
        top_violin_genes <- unique(top_violin_genes)
        top_violin_genes <- top_violin_genes[top_violin_genes %in% rownames(seurat_obj)]
        
        if(length(top_violin_genes) > 0) {
          p_violin <- VlnPlot(seurat_obj, features = head(top_violin_genes, 6), 
                              pt.size = 0, ncol = 2)
          ggsave("simple_violin.pdf", plot = p_violin, width = 12, height = 8)
          result_plots$violin <- p_violin
        }
      }
    }, error = function(e) {
      message("小提琴图创建失败: ", e$message)
    })
    
    # 创建简单的HTML报告
    message("创建简单HTML报告...")
    tryCatch({
      html_content <- paste0(
        "<!DOCTYPE html>
        <html>
        <head>
          <title>简化标记基因报告</title>
          <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1, h2 { color: #333366; }
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            tr:nth-child(even) { background-color: #f9f9f9; }
          </style>
        </head>
        <body>
          <h1>单细胞标记基因分析报告 (简化版)</h1>
          <p>分析时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>
          <p>标识变量: ", ident_col, "</p>
          
          <h2>标记基因列表</h2>
          <p>各细胞类型的顶部标记基因可在以下文件中找到: <a href='top_markers_simple.csv'>top_markers_simple.csv</a></p>
          <p>完整的标记基因列表可在以下文件中找到: <a href='all_markers_simple.csv'>all_markers_simple.csv</a></p>
          
          <h2>可视化结果</h2>
          <p>以下文件包含生成的可视化:</p>
          <ul>")
      
      if(file.exists("simple_heatmap.pdf")) {
        html_content <- paste0(html_content, "
            <li><a href='simple_heatmap.pdf'>标记基因热图</a></li>")
      }
      
      if(file.exists("simple_dotplot.pdf")) {
        html_content <- paste0(html_content, "
            <li><a href='simple_dotplot.pdf'>标记基因点图</a></li>")
      }
      
      if(file.exists("simple_violin.pdf")) {
        html_content <- paste0(html_content, "
            <li><a href='simple_violin.pdf'>标记基因小提琴图</a></li>")
      }
      
      html_content <- paste0(html_content, "
          </ul>
        </body>
        </html>")
      
      # 保存HTML
      writeLines(html_content, "simple_marker_report.html")
    }, error = function(e) {
      message("HTML报告创建失败: ", e$message)
    })
    
    message("简化版标记基因可视化完成")
    return(list(
      markers = markers, 
      top_markers = top_markers,
      plots = result_plots
    ))
  }, error = function(e) {
    message("简化版标记基因分析失败: ", e$message)
    return(NULL)
  })
}

# 使用示例
# library(Seurat)
# 
# # 基本用法
# seurat_obj <- readRDS("path/to/your/seuraEpithelial_object.rds")
results <- generate_marker_viz(
  seurat_obj = Epithelial_object,
  ident_col = "Annotation_2",
  top_n = 5,
  output_dir = "Marker_Results"
)
# 
# # 如果主函数失败，尝试简化版
# # simple_results <- simple_marker_viz(seurat_obj)
################################################################################################################################
################################################################################################################################
################################################################################################################################
# 创建输出目录
dir.create("7-pairwise_tissue_comparison", showWarnings = FALSE)
setwd("7-pairwise_tissue_comparison")
message("开始进行组织间细胞类型两两比较分析...")

# 确保DESeq2已加载
if(!require("DESeq2")) {
  message("安装DESeq2包...")
  if(!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("DESeq2")
  library(DESeq2)
} else {
  library(DESeq2)
}

# 首先检查Epithelial_object是否包含必要的元数据列
required_columns <- c("Annotation_2", "tissue", "sample")
missing_columns <- required_columns[!required_columns %in% colnames(Epithelial_object@meta.data)]

if(length(missing_columns) > 0) {
  stop(paste("Epithelial_object缺少必要的元数据列:", paste(missing_columns, collapse=", ")))
}

# 检查不同组织中的细胞类型分布
message("分析不同组织中的细胞类型分布...")
tissue_cell_distribution <- table(Epithelial_object$tissue, Epithelial_object$Annotation_2)
write.csv(tissue_cell_distribution, "tissue_celltype_distribution.csv")

# 打印组织和细胞类型信息
message("数据集包含以下组织类型:")
print(table(Epithelial_object$tissue))

message("数据集包含以下T细胞亚型:")
print(table(Epithelial_object$Annotation_2))

# 1. 按组织和细胞类型进行pseudobulk聚合
message("执行基于组织和细胞类型的pseudobulk聚合...")
av_tissue_celltype <- AggregateExpression(
  Epithelial_object,
  group.by = c("tissue", "sample", "Annotation_2"),  # 按组织、样本和细胞类型分组
  assays = "RNA",
  slot = "counts",  # 使用原始计数
  return.seurat = FALSE
)

# 将结果转换为数据框
av_tissue_celltype_df <- as.data.frame(av_tissue_celltype[[1]])
write.csv(av_tissue_celltype_df, file = "pseudobulk_tissue_sample_celltype.csv")

# 2. 准备细胞类型间组织比较分析
message("准备细胞类型间组织比较分析...")

# 从聚合矩阵的列名中提取组织、样本和细胞类型信息
extract_metadata <- function(column_names) {
  # 假设列名格式为 "tissue_sample_celltype"
  parts <- strsplit(column_names, "_")
  
  # 创建结果数据框
  result <- data.frame(
    column = column_names,
    tissue = sapply(parts, function(x) x[1]),
    sample = sapply(parts, function(x) x[2]),
    celltype = sapply(parts, function(x) paste(x[-(1:2)], collapse="_")),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# 提取元数据
metadata <- extract_metadata(colnames(av_tissue_celltype_df))
rownames(metadata) <- metadata$column

# 检查提取的元数据
message("提取的元数据示例(前5行):")
print(head(metadata, 5))

# 3. 为每个细胞类型执行所有可能的组织两两比较
message("执行细胞类型的组织间两两比较分析...")

run_pairwise_tissue_comparison <- function(counts_matrix, metadata, cell_type) {
  # 筛选特定细胞类型的样本
  cell_indices <- which(metadata$celltype == cell_type)
  
  if(length(cell_indices) < 6) {  # 至少需要6个样本(每组至少3个)
    message(paste("细胞类型", cell_type, "的样本数量不足(", length(cell_indices), ")，跳过分析"))
    return(NULL)
  }
  
  # 提取该细胞类型的计数数据和元数据
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
  
  # 过滤掉样本数不足的组织
  valid_indices <- cell_metadata$tissue %in% valid_tissues
  cell_counts <- cell_counts[, valid_indices, drop = FALSE]
  cell_metadata <- cell_metadata[valid_indices, , drop = FALSE]
  
  # 确保counts_matrix是整数
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
    
    # 只选择这两个组织的样本
    pair_indices <- cell_metadata$tissue %in% pair
    pair_counts <- cell_counts[, pair_indices, drop = FALSE]
    pair_metadata <- cell_metadata[pair_indices, , drop = FALSE]
    
    # 再次检查每个组织的样本数
    tissue_sample_counts <- table(pair_metadata$tissue)
    if(any(tissue_sample_counts < 3)) {
      message(paste("跳过比较: 某个组织样本数少于3:", 
                    paste(names(tissue_sample_counts), tissue_sample_counts, sep="=", collapse=", ")))
      next
    }
    
    # 创建DESeq2数据集
    dds <- DESeqDataSetFromMatrix(
      countData = pair_counts,
      colData = pair_metadata,
      design = ~ tissue
    )
    
    # 设置参考水平 (任意选择其中一个组织)
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
    res_df$regulation <- ifelse(res_df$padj > 0.05, "stable",
                                ifelse(abs(res_df$log2FoldChange) < 1, "stable",
                                       ifelse(res_df$log2FoldChange >= 1, "up", "down")))
    
    # 保存结果
    results_list[[comparison_name]] <- res_df
    
    # 创建输出目录
    result_dir <- paste0(cell_type, "_", tissue1, "_vs_", tissue2)
    result_dir <- gsub("[^a-zA-Z0-9_]", "-", result_dir)  # 替换非法字符
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
    
    write.csv(summary_stats, file = file.path(result_dir, "summary_stats.csv"), row.names = FALSE)
    
    # 创建火山图
    volcano_file <- file.path(result_dir, "volcano_plot.pdf")
    pdf(volcano_file, width = 10, height = 8)
    p <- ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
      geom_point(size = 1.5, alpha = 0.7, aes(color = regulation)) +
      scale_color_manual(values = c("down" = "#00468B", "stable" = "gray", "up" = "#E64B35")) +
      labs(x = "Log2(fold change)", 
           y = "-log10(adjusted p-value)",
           title = paste0(cell_type, ": ", tissue2, " vs ", tissue1)) +
      geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'black', linewidth = 0.5) + 
      geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'black', linewidth = 0.5) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    print(p)
    dev.off()
    
    # 统计差异基因数量
    message(paste0(cell_type, " ", comparison_name, ": ",
                   sum(res_df$regulation == "up"), " 上调基因, ",
                   sum(res_df$regulation == "down"), " 下调基因"))
    
    # 创建主要差异基因热图 - 修复版本
    top_degs <- res_df %>%
      filter(regulation != "stable") %>%
      arrange(padj) %>%
      head(50)  # 取前50个差异基因
    
    if(nrow(top_degs) > 0) {
      # 提取这些基因在所有样本中的表达值
      gene_expr <- pair_counts[top_degs$gene, , drop = FALSE]
      
      # 检查数据有效性
      message(paste0("热图数据维度: ", nrow(gene_expr), " × ", ncol(gene_expr)))
      message(paste0("差异基因数量: ", nrow(top_degs)))
      
      # 标准化 - 改进方法
      gene_expr_norm <- tryCatch({
        # 先进行log转换
        log_expr <- log2(gene_expr + 1)
        
        # 检查是否有全零的行
        non_zero_rows <- rowSums(log_expr) > 0
        if(sum(non_zero_rows) == 0) {
          message("所有基因的表达值都为0，无法创建热图")
          NULL
        } else {
          # 只对非零行进行标准化
          log_expr_filtered <- log_expr[non_zero_rows, , drop = FALSE]
          
          # 行标准化
          normalized <- t(scale(t(log_expr_filtered)))
          
          # 处理可能的NA值
          normalized[is.na(normalized)] <- 0
          
          normalized
        }
      }, error = function(e) {
        message(paste("标准化出错:", e$message))
        NULL
      })
      
      if(!is.null(gene_expr_norm) && nrow(gene_expr_norm) > 0) {
        # 创建热图注释
        anno_col <- data.frame(Tissue = pair_metadata$tissue)
        rownames(anno_col) <- colnames(gene_expr_norm)
        
        # 绘制热图 - 使用多种方法尝试
        heatmap_file <- file.path(result_dir, "DEG_heatmap.pdf")
        pdf(heatmap_file, width = 12, height = 10)
        
        heatmap_success <- FALSE
        
        tryCatch({
          # 方法1: 使用pheatmap
          message("尝试使用pheatmap绘制DEG热图...")
          
          # 检查数据范围
          data_range <- range(gene_expr_norm, na.rm = TRUE)
          message(paste0("标准化数据范围: [", round(data_range[1], 3), ", ", round(data_range[2], 3), "]"))
          
          hmap <- pheatmap(gene_expr_norm,
                           annotation_col = anno_col,
                           main = paste0(cell_type, " DEGs: ", tissue2, " vs ", tissue1),
                           fontsize_row = 8,
                           fontsize_col = 8,
                           show_colnames = FALSE,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           scale = "none",  # 已经标准化过了
                           color = colorRampPalette(c("blue", "white", "red"))(50))
          
          # 确保打印热图
          print(hmap)
          message("成功使用pheatmap创建DEG热图")
          heatmap_success <- TRUE
          
        }, error = function(e) {
          message(paste0("pheatmap失败: ", e$message, ". 尝试基础热图..."))
          
          # 方法2: 使用基础heatmap函数
          tryCatch({
            # 创建颜色注释
            tissue_colors <- rainbow(length(unique(pair_metadata$tissue)))
            names(tissue_colors) <- unique(pair_metadata$tissue)
            col_colors <- tissue_colors[pair_metadata$tissue]
            
            # 绘制基础热图
            heatmap(as.matrix(gene_expr_norm),
                    main = paste0(cell_type, " DEGs: ", tissue2, " vs ", tissue1),
                    col = colorRampPalette(c("blue", "white", "red"))(50),
                    margins = c(5, 8),
                    cexRow = 0.8,
                    cexCol = 0.8,
                    ColSideColors = col_colors,
                    labCol = rep("", ncol(gene_expr_norm)))  # 隐藏列名
            
            # 添加图例
            legend("topright", 
                   legend = names(tissue_colors),
                   fill = tissue_colors,
                   title = "Tissue",
                   cex = 0.8,
                   bty = "n")
            
            message("成功使用基础heatmap创建DEG热图")
            heatmap_success <- TRUE
            
          }, error = function(e2) {
            message(paste0("基础heatmap也失败: ", e2$message, ". 尝试image函数..."))
            
            # 方法3: 使用image函数
            tryCatch({
              # 准备数据
              plot_data <- as.matrix(gene_expr_norm)
              
              # 绘制image
              image(1:ncol(plot_data), 1:nrow(plot_data), t(plot_data),
                    main = paste0(cell_type, " DEGs: ", tissue2, " vs ", tissue1),
                    xlab = "Samples", ylab = "Genes",
                    col = colorRampPalette(c("blue", "white", "red"))(50),
                    axes = FALSE)
              
              # 添加轴
              axis(1, at = 1:ncol(plot_data), labels = FALSE)
              axis(2, at = 1:min(20, nrow(plot_data)), 
                   labels = rownames(plot_data)[1:min(20, nrow(plot_data))], 
                   las = 2, cex.axis = 0.6)
              
              # 添加组织分组线
              tissue_breaks <- cumsum(table(pair_metadata$tissue))
              abline(v = tissue_breaks[-length(tissue_breaks)] + 0.5, col = "black", lwd = 2)
              
              # 添加组织标签
              tissue_centers <- c(0, tissue_breaks[-length(tissue_breaks)]) + 
                diff(c(0, tissue_breaks)) / 2
              mtext(names(table(pair_metadata$tissue)), side = 1, at = tissue_centers, line = 1)
              
              message("成功使用image函数创建DEG热图")
              heatmap_success <- TRUE
              
            }, error = function(e3) {
              message(paste0("所有热图方法都失败: ", e3$message))
              
              # 方法4: 创建一个简单的错误消息图
              plot(1, type = "n", axes = FALSE, ann = FALSE)
              text(1, 1, paste0("DEG热图创建失败\n", 
                                "基因数: ", nrow(gene_expr_norm), "\n",
                                "样本数: ", ncol(gene_expr_norm)), 
                   cex = 1.2, col = "red")
            })
          })
        }, finally = {
          # 确保关闭PDF设备
          dev.off()
          
          if(heatmap_success) {
            message(paste0("DEG热图已保存到: ", heatmap_file))
          } else {
            message(paste0("DEG热图创建失败，但PDF文件已创建: ", heatmap_file))
          }
        })
      } else {
        message("数据标准化失败或没有有效的差异基因，跳过热图创建")
      }
    } else {
      message("没有找到差异表达基因，跳过热图创建")
    }
    
    # 执行GO富集分析(如果有clusterProfiler) - 整合版本
    if(nrow(top_degs) > 0 && 
       requireNamespace("clusterProfiler", quietly = TRUE) && 
       requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      
      message("执行GO富集分析...")
      library(clusterProfiler)
      library(org.Hs.eg.db)
      
      # 提取上调和下调基因
      up_genes <- res_df$gene[res_df$regulation == "up"]
      down_genes <- res_df$gene[res_df$regulation == "down"]
      
      # 存储GO富集结果
      go_results <- list()
      
      # GO富集分析 - 上调基因
      if(length(up_genes) >= 10) {
        tryCatch({
          ego_up <- enrichGO(gene = up_genes,
                             OrgDb = org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)
          
          if(!is.null(ego_up) && nrow(ego_up) > 0) {
            go_up_df <- as.data.frame(ego_up)
            go_up_df$regulation <- "upregulated"
            go_results[["upregulated"]] <- go_up_df
            message(paste0("上调基因GO富集: ", nrow(go_up_df), "个显著通路"))
          } else {
            message("上调基因GO富集分析无显著结果")
          }
        }, error = function(e) {
          message(paste0("上调基因GO富集分析失败: ", e$message))
        })
      } else {
        message(paste0("上调基因数量不足(", length(up_genes), ")，跳过GO富集分析"))
      }
      
      # GO富集分析 - 下调基因
      if(length(down_genes) >= 10) {
        tryCatch({
          ego_down <- enrichGO(gene = down_genes,
                               OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL",
                               ont = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)
          
          if(!is.null(ego_down) && nrow(ego_down) > 0) {
            go_down_df <- as.data.frame(ego_down)
            go_down_df$regulation <- "downregulated"
            go_results[["downregulated"]] <- go_down_df
            message(paste0("下调基因GO富集: ", nrow(go_down_df), "个显著通路"))
          } else {
            message("下调基因GO富集分析无显著结果")
          }
        }, error = function(e) {
          message(paste0("下调基因GO富集分析失败: ", e$message))
        })
      } else {
        message(paste0("下调基因数量不足(", length(down_genes), ")，跳过GO富集分析"))
      }
      
      # 整合GO富集结果
      if(length(go_results) > 0) {
        # 合并上调和下调的GO结果
        combined_go <- do.call(rbind, go_results)
        rownames(combined_go) <- NULL
        
        # 保存整合的GO结果
        write.csv(combined_go, 
                  file = file.path(result_dir, "GO_enrichment_combined.csv"), 
                  row.names = FALSE)
        
        message(paste0("整合GO富集结果: ", nrow(combined_go), "个显著通路"))
        
        # 创建整合的GO富集可视化
        pdf(file.path(result_dir, "GO_enrichment_combined.pdf"), width = 14, height = 10)
        
        tryCatch({
          # 方法1: 创建对比点图
          if("upregulated" %in% names(go_results) && "downregulated" %in% names(go_results)) {
            # 双向GO富集对比
            
            # 选择每组的前10个通路
            top_up <- go_results[["upregulated"]] %>% 
              arrange(p.adjust) %>% 
              head(10)
            
            top_down <- go_results[["downregulated"]] %>% 
              arrange(p.adjust) %>% 
              head(10)
            
            # 合并用于可视化
            plot_data <- rbind(
              data.frame(
                Description = top_up$Description,
                log10pvalue = -log10(top_up$p.adjust),
                Count = top_up$Count,
                Regulation = "Upregulated",
                GeneRatio = sapply(strsplit(top_up$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
              ),
              data.frame(
                Description = top_down$Description,
                log10pvalue = -log10(top_down$p.adjust) * -1,  # 负值用于向下显示
                Count = top_down$Count,
                Regulation = "Downregulated",
                GeneRatio = sapply(strsplit(top_down$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])) * -1  # 负值
              )
            )
            
            # 创建双向条形图
            p1 <- ggplot(plot_data, aes(x = reorder(Description, abs(log10pvalue)), y = log10pvalue, fill = Regulation)) +
              geom_bar(stat = "identity") +
              scale_fill_manual(values = c("Upregulated" = "#E64B35", "Downregulated" = "#00468B")) +
              coord_flip() +
              labs(title = paste0("GO Enrichment Comparison: ", cell_type, " (", tissue2, " vs ", tissue1, ")"),
                   x = "GO Terms", 
                   y = "-log10(adjusted p-value)",
                   subtitle = "Upregulated (positive) vs Downregulated (negative)") +
              theme_bw() +
              theme(axis.text.y = element_text(size = 8),
                    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 10)) +
              geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
            
            print(p1)
            
            # 添加新页面，绘制气泡图对比
            grid::grid.newpage()
            
            # 重新处理数据用于气泡图
            bubble_data <- rbind(
              data.frame(
                Description = top_up$Description,
                log10pvalue = -log10(top_up$p.adjust),
                Count = top_up$Count,
                Regulation = "Upregulated",
                GeneRatio = sapply(strsplit(top_up$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
              ),
              data.frame(
                Description = top_down$Description,
                log10pvalue = -log10(top_down$p.adjust),
                Count = top_down$Count,
                Regulation = "Downregulated",
                GeneRatio = sapply(strsplit(top_down$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
              )
            )
            
            p2 <- ggplot(bubble_data, aes(x = GeneRatio, y = reorder(Description, log10pvalue), 
                                          size = Count, color = log10pvalue)) +
              geom_point(alpha = 0.8) +
              scale_size_continuous(range = c(2, 8), name = "Gene Count") +
              scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
              facet_wrap(~Regulation, scales = "free_y", ncol = 2) +
              labs(title = paste0("GO Enrichment Bubble Plot: ", cell_type, " (", tissue2, " vs ", tissue1, ")"),
                   x = "Gene Ratio", 
                   y = "GO Terms") +
              theme_bw() +
              theme(axis.text.y = element_text(size = 7),
                    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
                    strip.text = element_text(face = "bold"))
            
            print(p2)
            
          } else {
            # 只有一种类型的基因有GO结果
            single_result <- go_results[[1]]
            regulation_type <- names(go_results)[1]
            
            # 简单点图
            p3 <- ggplot(single_result[1:min(20, nrow(single_result)), ], 
                         aes(x = Count, y = reorder(Description, Count))) +
              geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
              scale_size_continuous(range = c(2, 8), name = "Gene Count") +
              scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
              labs(title = paste0("GO Enrichment: ", regulation_type, " genes"),
                   subtitle = paste0(cell_type, " (", tissue2, " vs ", tissue1, ")"),
                   x = "Gene Count", 
                   y = "GO Terms") +
              theme_bw() +
              theme(axis.text.y = element_text(size = 8),
                    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
            
            print(p3)
          }
          
          message("成功创建整合的GO富集可视化")
          
        }, error = function(e) {
          message(paste0("GO富集可视化绘制失败: ", e$message))
          plot(1, type = "n", axes = FALSE, ann = FALSE)
          text(1, 1, "GO富集可视化绘制失败\n请检查数据格式", cex = 1.2)
        })
        
        dev.off()
        
        # 创建GO富集摘要表
        go_summary <- combined_go %>%
          group_by(regulation) %>%
          summarise(
            Total_pathways = n(),
            Avg_pvalue = mean(p.adjust),
            Max_gene_count = max(Count),
            .groups = "drop"
          )
        
        write.csv(go_summary, 
                  file = file.path(result_dir, "GO_enrichment_summary.csv"), 
                  row.names = FALSE)
        
      } else {
        message("没有显著的GO富集结果")
      }
      
    } else {
      if(nrow(top_degs) == 0) {
        message("没有差异基因，跳过GO富集分析")
      } else {
        message("缺少clusterProfiler或org.Hs.eg.db包，跳过GO富集分析")
      }
    }
  }
  
  return(results_list)
}

# 4. 对每个细胞类型执行组织间两两比较
all_cell_types <- unique(metadata$celltype)
results_by_celltype <- list()

for(cell_type in all_cell_types) {
  message(paste("分析细胞类型:", cell_type))
  
  # 跳过NA或空字符串
  if(is.na(cell_type) || cell_type == "") {
    message("跳过无效细胞类型")
    next
  }
  
  # 创建细胞类型目录
  cell_dir <- gsub("[^a-zA-Z0-9]", "_", cell_type)
  dir.create(cell_dir, showWarnings = FALSE)
  
  # 设置工作目录
  original_dir <- getwd()
  setwd(cell_dir)
  
  # 执行组织间两两比较
  tryCatch({
    results <- run_pairwise_tissue_comparison(
      av_tissue_celltype_df, 
      metadata, 
      cell_type
    )
    
    if(!is.null(results)) {
      results_by_celltype[[cell_type]] <- results
    }
  }, error = function(e) {
    message(paste("分析细胞类型", cell_type, "时出错:", e$message))
  })
  
  # 返回原目录
  setwd(original_dir)
}

# 5. 创建组织间比较的总结报告
message("创建组织间比较的总结报告...")

# 收集所有比较的结果
all_comparisons <- data.frame()

# 遍历目录结构收集结果
for(cell_dir in list.dirs(recursive = FALSE)) {
  cell_type <- basename(cell_dir)
  
  # 跳过非细胞类型目录
  if(!dir.exists(cell_dir) || cell_dir == "." || cell_dir == "..") {
    next
  }
  
  # 查找所有比较目录
  comparison_dirs <- list.dirs(path = cell_dir, recursive = FALSE)
  
  for(comp_dir in comparison_dirs) {
    # 检查是否有summary_stats.csv
    stats_file <- file.path(comp_dir, "summary_stats.csv")
    
    if(file.exists(stats_file)) {
      stats <- read.csv(stats_file)
      stats$CellType <- cell_type
      
      all_comparisons <- rbind(all_comparisons, stats)
    }
  }
}

# 保存总结表
if(nrow(all_comparisons) > 0) {
  write.csv(all_comparisons, "all_pairwise_comparisons_summary.csv", row.names = FALSE)
  
  # 创建总结图
  # 条形图：每个细胞类型的差异基因数
  pdf("diff_gene_counts_by_celltype.pdf", width = 15, height = 10)
  p <- ggplot(all_comparisons, aes(x = CellType, y = Total_DEGs, fill = Comparison)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "差异基因数量（按细胞类型）", 
         x = "细胞类型", 
         y = "差异基因数量")
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
  p <- ggplot(all_comparisons_long, aes(x = Comparison, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ CellType, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = c("Up_regulated" = "#E64B35", "Down_regulated" = "#00468B")) +
    labs(title = "上调和下调基因数量（按细胞类型和组织比较）", 
         x = "组织比较", 
         y = "基因数量")
  print(p)
  dev.off()
  
  # 热图：所有组织比较的差异基因数
  # 准备热图数据
  heatmap_data <- dcast(all_comparisons, CellType ~ Comparison, value.var = "Total_DEGs")
  rownames(heatmap_data) <- heatmap_data$CellType
  heatmap_data <- heatmap_data[, -1]
  
  # 替换NA值为0
  heatmap_data[is.na(heatmap_data)] <- 0
  
  # 绘制热图
  pdf("comparison_heatmap.pdf", width = 12, height = 10)
  pheatmap(heatmap_data,
           display_numbers = TRUE,
           main = "组织间比较的差异基因数量热图",
           fontsize_number = 8,
           fontsize = 10)
  dev.off()
}

message("组织间细胞类型两两比较分析完成！")

# 返回主目录
setwd("..")

################################################################################################################################
################################################################################################################################
################################################################################################################################
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

# 检查组织注释是否存在
if(!"tissue" %in% colnames(Epithelial_object@meta.data)) {
  message("在Seurat对象中找不到'tissue'列。跳过组织特异性分析。")
} else {
  # 获取所有组织类型
  tissue_types <- unique(Epithelial_object$tissue)
  tissue_types <- tissue_types[!is.na(tissue_types)]
  message(paste0("共检测到", length(tissue_types), "个组织类型"))
}

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
if(!exists("Epithelial_object")) {
  stop("在工作环境中找不到'Epithelial_object' Seurat对象。请先加载数据。")
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
  cells_in_tissue <- rownames(Epithelial_object@meta.data[Epithelial_object@meta.data$tissue == tissue, ])
  message("DEBUG: Number of cells = ", length(cells_in_tissue))
  
  if(length(cells_in_tissue) < 50) {
    message(paste0("组织'", tissue, "'的细胞数量太少(", length(cells_in_tissue), ")，跳过"))
    next
  }
  
  # 创建该组织的子集
  tissue_subset <- subset(Epithelial_object, cells = cells_in_tissue)
  
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
if(!"Annotation_2" %in% colnames(Epithelial_object@meta.data)) {
  stop("在Seurat对象中找不到'Annotation_2'。请先进行细胞类型注释。")
}

# 设置细胞类型为标识
Idents(Epithelial_object) <- Epithelial_object$Annotation_2
message("按细胞类型计算平均表达...")

# 计算每个细胞类型的平均表达
expr <- AverageExpression(Epithelial_object, assays = "RNA", slot = "data")[[1]]
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
if(!"tissue" %in% colnames(Epithelial_object@meta.data)) {
  message("在Seurat对象中找不到'tissue'列。跳过组织特异性分析。")
} else {
  # 获取所有组织类型
  tissue_types <- unique(Epithelial_object$tissue)
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
  cells_in_tissue <- rownames(Epithelial_object@meta.data[Epithelial_object@meta.data$tissue == tissue, ])
  
  if(length(cells_in_tissue) < 50) {
    message(paste0("Tissue '", tissue, "' has too few cells (", length(cells_in_tissue), "), skipping"))
    next
  }
  
  # 创建该组织的子集
  tissue_subset <- subset(Epithelial_object, cells = cells_in_tissue)
  
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
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' 综合基因功能可视化分析函数
#' Comprehensive Gene Function Visualization Analysis Function
#' 
#' @description 
#' 该函数接受功能标注的基因列表，结合Seurat对象和pseudoBulk结果，
#' 生成热图(heatmap)、降维图(dimplot)、点图(dotplot)和密度图(densityplot)
#' 
#' @param gene_list 命名列表，包含基因名称和功能标注 (named list with gene names and functional annotations)
#' @param seurat_obj Seurat对象 (Seurat object)
#' @param pseudobulk_data pseudoBulk表达矩阵 (pseudoBulk expression matrix)
#' @param group_by 分组变量名称 (grouping variable name)，默认为"seurat_clusters"
#' @param assay 使用的检测方法 (assay to use)，默认为"RNA" 
#' @param reduction 降维方法 (reduction method)，默认为"umap"
#' @param output_dir 输出目录 (output directory)，默认为"Gene_Function_Visualization"
#' @param prefix 输出文件前缀 (output file prefix)，默认为"Gene_Analysis"
#' @param width PDF图片宽度 (PDF width)，默认为12
#' @param height PDF图片高度 (PDF height)，默认为10
#' @param pt_size 点大小 (point size)，默认为0.5
#' @param max_plot_width 最大图像宽度限制 (maximum plot width)，默认为45英寸
#' @param max_plot_height 最大图像高度限制 (maximum plot height)，默认为30英寸
#' @param verbose 是否显示详细信息 (verbose output)，默认为TRUE
#' 
#' @return 返回包含所有图形对象和分析结果的列表 (list containing all plot objects and analysis results)
#' 
#' @examples
#' # 示例基因列表 (Example gene list)
#' gene_functions <- list(
#'   "T细胞标记" = c("CD3E", "CD3D", "CD5", "TRAC"),
#'   "髓系细胞标记" = c("CD14", "CD68", "ITGAM", "ITGAX"),
#'   "炎症相关" = c("IFNG", "IL10", "TNF")
#' )
#' 
#' # 调用函数 (Function call)
#' results <- comprehensive_gene_visualization(
#'   gene_list = gene_functions,
#'   seurat_obj = your_seurat_object,
#'   pseudobulk_data = your_pseudobulk_matrix,
#'   group_by = "cell_type"
#' )

comprehensive_gene_visualization <- function(
    gene_list,
    seurat_obj,
    pseudobulk_data = NULL,
    group_by = "seurat_clusters",
    assay = "RNA",
    reduction = "umap",
    output_dir = "Gene_Function_Visualization",
    prefix = "Gene_Analysis",
    width = 12,
    height = 10,
    pt_size = 0.5,
    max_plot_width = 45,    # 新增：最大图像宽度限制
    max_plot_height = 30,   # 新增：最大图像高度限制
    verbose = TRUE
) {
  
  # ============================================================================
  # 1. 环境准备和参数验证 (Environment Setup and Parameter Validation)
  # ============================================================================
  
  if(verbose) message("开始综合基因功能可视化分析...")
  
  # 加载必要的R包 (Load required packages)
  required_packages <- c("Seurat", "ggplot2", "dplyr", "pheatmap", 
                         "RColorBrewer", "viridis", "cowplot", "patchwork")
  
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("请安装R包:", pkg, "(Please install package:", pkg, ")"))
    }
  }
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(pheatmap)
    library(RColorBrewer)
    library(viridis)
    library(cowplot)
    library(patchwork)
  })
  
  # 验证输入参数 (Validate input parameters)
  if(missing(gene_list) || is.null(gene_list)) {
    stop("必须提供基因列表 (gene_list is required)")
  }
  
  if(missing(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("必须提供有效的Seurat对象 (Valid Seurat object is required)")
  }
  
  if(!group_by %in% colnames(seurat_obj@meta.data)) {
    warning(paste("分组变量", group_by, "在meta.data中不存在，使用默认聚类"))
    group_by <- "seurat_clusters"
  }
  
  # 验证图像尺寸参数 (Validate plot dimension parameters)
  if(max_plot_width > 50 || max_plot_height > 50) {
    warning("图像尺寸过大可能导致保存失败，建议设置max_plot_width和max_plot_height <= 45")
  }
  
  if(width > max_plot_width || height > max_plot_height) {
    if(verbose) message("调整基础图像尺寸以符合最大限制...")
    width <- min(width, max_plot_width)
    height <- min(height, max_plot_height)
  }
  
  # 创建输出目录 (Create output directory)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # 初始化结果列表 (Initialize result list)
  results <- list(
    parameters = list(
      gene_list = gene_list,
      group_by = group_by,
      assay = assay,
      reduction = reduction,
      timestamp = timestamp
    ),
    plots = list(),
    data = list(),
    files = list()
  )
  
  # ============================================================================
  # 2. 数据预处理 (Data Preprocessing)  
  # ============================================================================
  
  if(verbose) message("进行数据预处理...")
  
  # 提取所有基因名称 (Extract all gene names)
  all_genes <- unique(unlist(gene_list))
  
  # 检查基因在Seurat对象中的存在情况 (Check gene availability)
  genes_in_data <- all_genes[all_genes %in% rownames(seurat_obj)]
  genes_missing <- all_genes[!all_genes %in% rownames(seurat_obj)]
  
  if(verbose) {
    message(paste("总基因数量:", length(all_genes)))
    message(paste("数据中存在的基因:", length(genes_in_data)))
    message(paste("缺失的基因:", length(genes_missing)))
    if(length(genes_missing) > 0) {
      message("缺失基因列表:", paste(genes_missing, collapse = ", "))
    }
  }
  
  # 更新基因列表，只保留存在的基因 (Update gene list with available genes)
  gene_list_filtered <- lapply(gene_list, function(genes) {
    intersect(genes, genes_in_data)
  })
  
  # 移除空的功能类别 (Remove empty functional categories)
  gene_list_filtered <- gene_list_filtered[sapply(gene_list_filtered, length) > 0]
  
  if(length(gene_list_filtered) == 0) {
    stop("没有找到任何有效的基因用于分析 (No valid genes found for analysis)")
  }
  
  # 设置分组标识 (Set grouping identity)
  Idents(seurat_obj) <- group_by
  
  # 记录基因信息 (Record gene information)
  results$data$gene_availability <- list(
    total_genes = length(all_genes),
    available_genes = genes_in_data,
    missing_genes = genes_missing,
    filtered_gene_list = gene_list_filtered
  )
  
  # ============================================================================
  # 3. 生成热图 (Generate Heatmap)
  # ============================================================================
  
  if(verbose) message("生成表达热图...")
  
  tryCatch({
    # 3.1 单细胞数据热图 (Single-cell heatmap)
    genes_for_heatmap <- unique(unlist(gene_list_filtered))
    
    # 确保基因已经标准化 (Ensure genes are scaled)
    genes_to_scale <- genes_for_heatmap[!genes_for_heatmap %in% rownames(seurat_obj[[assay]]@scale.data)]
    if(length(genes_to_scale) > 0) {
      if(verbose) message("对新基因进行标准化...")
      seurat_obj <- ScaleData(seurat_obj, features = genes_to_scale, verbose = FALSE)
    }
    
    # Seurat热图 (Seurat heatmap)
    p_heatmap_seurat <- DoHeatmap(
      seurat_obj, 
      features = genes_for_heatmap,
      group.by = group_by,
      assay = assay,
      raster = TRUE
    ) + 
      scale_fill_viridis_c() +
      theme(axis.text.y = element_text(size = 8)) +
      ggtitle("Gene Expression Heatmap")
    
    results$plots$heatmap_seurat <- p_heatmap_seurat
    
    # 3.2 平均表达热图 (Average expression heatmap)
    avg_expr <- AverageExpression(
      seurat_obj, 
      features = genes_for_heatmap,
      group.by = group_by,
      assays = assay,
      slot = "data"
    )
    
    avg_expr_matrix <- avg_expr[[assay]]
    
    # Z-score标准化 (Z-score normalization)
    avg_expr_scaled <- tryCatch({
      t(scale(t(avg_expr_matrix)))
    }, error = function(e) {
      if(verbose) message("使用手动Z-score计算...")
      result <- avg_expr_matrix
      for(i in 1:nrow(avg_expr_matrix)) {
        row_mean <- mean(as.numeric(avg_expr_matrix[i,]), na.rm = TRUE)
        row_sd <- sd(as.numeric(avg_expr_matrix[i,]), na.rm = TRUE)
        if(row_sd > 0) {
          result[i,] <- (avg_expr_matrix[i,] - row_mean) / row_sd
        }
      }
      return(result)
    })
    
    # 准备功能注释 (Prepare functional annotations)
    gene_annotation <- data.frame(
      Gene = rownames(avg_expr_scaled),
      Function = NA,
      stringsAsFactors = FALSE
    )
    
    for(func_name in names(gene_list_filtered)) {
      func_genes <- gene_list_filtered[[func_name]]
      gene_annotation$Function[gene_annotation$Gene %in% func_genes] <- func_name
    }
    
    rownames(gene_annotation) <- gene_annotation$Gene
    gene_annotation$Gene <- NULL
    
    # 设置颜色 (Set colors)
    n_functions <- length(unique(gene_annotation$Function))
    function_colors <- setNames(
      rainbow(n_functions), 
      unique(gene_annotation$Function)
    )
    
    anno_colors <- list(Function = function_colors)
    
    # 保存pheatmap热图 (Save pheatmap heatmap)
    heatmap_file <- file.path(output_dir, paste0(prefix, "_average_expression_heatmap_", timestamp, ".pdf"))
    
    pdf(heatmap_file, width = width, height = height)
    pheatmap(
      avg_expr_scaled,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      annotation_row = gene_annotation,
      annotation_colors = anno_colors,
      fontsize_row = 8,
      fontsize_col = 10,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      main = "Average Expression Heatmap",
      angle_col = 45
    )
    dev.off()
    
    results$plots$heatmap_average <- avg_expr_scaled
    results$files$heatmap <- heatmap_file
    results$data$average_expression <- avg_expr_matrix
    
    if(verbose) message(paste("热图已保存:", heatmap_file))
    
  }, error = function(e) {
    warning(paste("热图生成失败:", e$message))
  })
  
  # ============================================================================
  # 4. 生成降维图 (Generate Dimension Reduction Plots)
  # ============================================================================
  
  if(verbose) message("生成降维可视化图...")
  
  tryCatch({
    # 4.1 按分组着色的降维图 (Dimplot by grouping)
    p_dimplot_group <- DimPlot(
      seurat_obj,
      reduction = reduction,
      group.by = group_by,
      pt.size = pt_size,
      raster = TRUE
    ) + 
      ggtitle(paste("Cell Grouping", toupper(reduction), "Plot")) +
      theme_minimal()
    
    results$plots$dimplot_group <- p_dimplot_group
    
    # 4.2 特征基因表达降维图 (Feature plots for key genes)
    feature_plots <- list()
    
    # 选择每个功能类别的代表基因 (Select representative genes)
    key_genes <- c()
    for(func_name in names(gene_list_filtered)) {
      func_genes <- gene_list_filtered[[func_name]]
      # 每个功能类别最多选择3个基因 (Max 3 genes per function)
      key_genes <- c(key_genes, head(func_genes, 3))
    }
    key_genes <- unique(key_genes)
    
    if(length(key_genes) > 0) {
      for(gene in key_genes) {
        p_feature <- FeaturePlot(
          seurat_obj,
          features = gene,
          reduction = reduction,
          pt.size = pt_size,
          raster = TRUE
        ) + 
          ggtitle(gene) +
          theme_minimal()
        
        feature_plots[[gene]] <- p_feature
      }
      
      results$plots$feature_plots <- feature_plots
    }
    
  }, error = function(e) {
    warning(paste("降维图生成失败:", e$message))
  })
  
  # ============================================================================
  # 5. 生成点图 (Generate Dot Plot)
  # ============================================================================
  
  if(verbose) message("生成基因表达点图...")
  
  tryCatch({
    genes_for_dotplot <- unique(unlist(gene_list_filtered))
    
    p_dotplot <- DotPlot(
      seurat_obj,
      features = genes_for_dotplot,
      group.by = group_by,
      assay = assay,
      cols = c("lightgrey", "red"),
      dot.scale = 8
    ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      ggtitle("Gene Expression Dot Plot") +
      xlab("Genes") + 
      ylab("Cell Groups")
    
    results$plots$dotplot <- p_dotplot
    
    # 保存点图数据 (Save dot plot data)
    dotplot_data <- p_dotplot$data
    results$data$dotplot_data <- dotplot_data
    
  }, error = function(e) {
    warning(paste("点图生成失败:", e$message))
  })
  
  # ============================================================================
  # 6. 生成密度图 (Generate Density Plots)
  # ============================================================================
  
  if(verbose) message("生成基因表达密度图...")
  
  tryCatch({
    density_plots <- list()
    
    # 为每个功能类别生成密度图 (Generate density plots for each function)
    for(func_name in names(gene_list_filtered)) {
      func_genes <- gene_list_filtered[[func_name]]
      
      if(length(func_genes) > 0) {
        # 计算功能评分 (Calculate function score)
        if(length(func_genes) == 1) {
          func_score <- FetchData(seurat_obj, vars = func_genes[1])[,1]
        } else {
          # 使用AddModuleScore计算模块评分 (Use AddModuleScore for multiple genes)
          temp_obj <- AddModuleScore(
            seurat_obj, 
            features = list(func_genes), 
            name = "temp_score",
            verbose = FALSE
          )
          func_score <- temp_obj$temp_score1
        }
        
        # 获取分组信息 (Get grouping information)
        group_info <- seurat_obj@meta.data[[group_by]]
        
        # 创建数据框 (Create data frame)
        density_data <- data.frame(
          Score = func_score,
          Group = group_info,
          stringsAsFactors = FALSE
        )
        
        # 绘制密度图 (Plot density)
        p_density <- ggplot(density_data, aes(x = Score, fill = Group)) +
          geom_density(alpha = 0.6) +
          scale_fill_viridis_d() +
          theme_minimal() +
          ggtitle(paste(func_name, "Expression Density Distribution")) +
          xlab("Expression Score") +
          ylab("Density")
        
        density_plots[[func_name]] <- p_density
      }
    }
    
    results$plots$density_plots <- density_plots
    
  }, error = function(e) {
    warning(paste("密度图生成失败:", e$message))
  })
  
  # ============================================================================
  # 7. pseudoBulk数据分析 (PseudoBulk Data Analysis)
  # ============================================================================
  
  if(!is.null(pseudobulk_data)) {
    if(verbose) message("分析pseudoBulk数据...")
    
    tryCatch({
      # 检查基因在pseudoBulk数据中的存在情况 (Check gene availability in pseudoBulk)
      pb_genes_available <- intersect(genes_in_data, rownames(pseudobulk_data))
      
      if(length(pb_genes_available) > 0) {
        # 提取pseudoBulk表达数据 (Extract pseudoBulk expression data)
        pb_expr <- pseudobulk_data[pb_genes_available, , drop = FALSE]
        
        # Z-score标准化 (Z-score normalization)
        pb_expr_scaled <- t(scale(t(pb_expr)))
        
        # 生成pseudoBulk热图 (Generate pseudoBulk heatmap)
        pb_heatmap_file <- file.path(output_dir, paste0(prefix, "_pseudobulk_heatmap_", timestamp, ".pdf"))
        
        pdf(pb_heatmap_file, width = width, height = height)
        pheatmap(
          pb_expr_scaled,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          show_rownames = TRUE,
          show_colnames = TRUE,
          color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
          main = "PseudoBulk Expression Heatmap",
          fontsize_row = 8,
          fontsize_col = 10
        )
        dev.off()
        
        results$data$pseudobulk_expression <- pb_expr
        results$files$pseudobulk_heatmap <- pb_heatmap_file
        
        if(verbose) message(paste("PseudoBulk热图已保存:", pb_heatmap_file))
      }
      
    }, error = function(e) {
      warning(paste("PseudoBulk分析失败:", e$message))
    })
  }
  
  # ============================================================================
  # 8. 保存综合结果 (Save Comprehensive Results)
  # ============================================================================
  
  if(verbose) message("保存分析结果...")
  
  # 8.1 保存单独的PDF文件 (Save individual PDF files)
  
  # 保存Seurat热图 (Save Seurat heatmap)
  if(!is.null(results$plots$heatmap_seurat)) {
    seurat_heatmap_file <- file.path(output_dir, paste0(prefix, "_seurat_heatmap_", timestamp, ".pdf"))
    # 智能计算热图尺寸
    heatmap_width <- min(max_plot_width, max(width, 20))
    heatmap_height <- min(max_plot_height, max(height, 16))
    ggsave(seurat_heatmap_file, results$plots$heatmap_seurat, 
           width = heatmap_width, height = heatmap_height, limitsize = FALSE)
    results$files$seurat_heatmap <- seurat_heatmap_file
  }
  
  # 保存降维图 (Save dimension reduction plots)
  if(!is.null(results$plots$dimplot_group)) {
    dimplot_file <- file.path(output_dir, paste0(prefix, "_dimplot_", timestamp, ".pdf"))
    ggsave(dimplot_file, results$plots$dimplot_group, width = width, height = height)
    results$files$dimplot <- dimplot_file
  }
  
  # 保存点图 (Save dot plot)
  if(!is.null(results$plots$dotplot)) {
    dotplot_file <- file.path(output_dir, paste0(prefix, "_dotplot_", timestamp, ".pdf"))
    # 智能计算点图宽度，但不超过最大限制
    calculated_width <- length(unique(unlist(gene_list_filtered))) * 0.25 + 6
    dot_width <- min(max_plot_width, max(width, calculated_width))
    if(verbose) message(paste("点图宽度设置为:", round(dot_width, 1), "英寸"))
    ggsave(dotplot_file, results$plots$dotplot, width = dot_width, height = height, limitsize = FALSE)
    results$files$dotplot <- dotplot_file
  }
  
  # 保存密度图 (Save density plots)
  if(length(results$plots$density_plots) > 0) {
    density_file <- file.path(output_dir, paste0(prefix, "_density_plots_", timestamp, ".pdf"))
    
    # 将所有密度图合并 (Combine all density plots)
    if(length(results$plots$density_plots) > 1) {
      combined_density <- wrap_plots(results$plots$density_plots, ncol = 2)
    } else {
      combined_density <- results$plots$density_plots[[1]]
    }
    
    # 智能计算密度图尺寸
    density_width <- min(max_plot_width, max(width, 16))
    density_height <- min(max_plot_height, height + ceiling(length(results$plots$density_plots)/2) * 4)
    
    ggsave(density_file, combined_density, 
           width = density_width, height = density_height, limitsize = FALSE)
    results$files$density_plots <- density_file
  }
  
  # 保存特征图 (Save feature plots)
  if(length(results$plots$feature_plots) > 0) {
    feature_file <- file.path(output_dir, paste0(prefix, "_feature_plots_", timestamp, ".pdf"))
    
    # 智能计算特征图尺寸 (Smart calculation for feature plot dimensions)
    feature_width <- min(max_plot_width, max(width, 14))
    feature_height <- min(max_plot_height, max(height, 10))
    
    pdf(feature_file, width = feature_width, height = feature_height)
    for(gene_name in names(results$plots$feature_plots)) {
      print(results$plots$feature_plots[[gene_name]])
    }
    dev.off()
    
    results$files$feature_plots <- feature_file
  }
  
  # 8.2 保存数据文件 (Save data files)
  
  # 保存基因可用性信息 (Save gene availability info)
  gene_info_file <- file.path(output_dir, paste0(prefix, "_gene_availability_", timestamp, ".txt"))
  sink(gene_info_file)
  cat("基因可用性分析报告\n")
  cat("==================\n\n")
  cat("总基因数量:", length(all_genes), "\n")
  cat("数据中存在的基因数量:", length(genes_in_data), "\n")
  cat("缺失的基因数量:", length(genes_missing), "\n\n")
  
  if(length(genes_missing) > 0) {
    cat("缺失的基因列表:\n")
    cat(paste(genes_missing, collapse = ", "), "\n\n")
  }
  
  cat("各功能类别的基因数量:\n")
  for(func_name in names(gene_list_filtered)) {
    cat(paste(func_name, ":", length(gene_list_filtered[[func_name]]), "个基因\n"))
  }
  sink()
  
  results$files$gene_availability <- gene_info_file
  
  # 保存点图数据 (Save dot plot data)
  if(!is.null(results$data$dotplot_data)) {
    dotplot_data_file <- file.path(output_dir, paste0(prefix, "_dotplot_data_", timestamp, ".csv"))
    write.csv(results$data$dotplot_data, dotplot_data_file, row.names = FALSE)
    results$files$dotplot_data <- dotplot_data_file
  }
  
  # 保存平均表达数据 (Save average expression data)
  if(!is.null(results$data$average_expression)) {
    avg_expr_file <- file.path(output_dir, paste0(prefix, "_average_expression_", timestamp, ".csv"))
    write.csv(results$data$average_expression, avg_expr_file)
    results$files$average_expression <- avg_expr_file
  }
  
  # ============================================================================
  # 9. 完成分析 (Complete Analysis)
  # ============================================================================
  
  if(verbose) {
    message("分析完成！")
    message(paste("结果保存在目录:", output_dir))
    message("生成的文件:")
    for(file_type in names(results$files)) {
      if(!is.null(results$files[[file_type]])) {
        message(paste("  -", file_type, ":", basename(results$files[[file_type]])))
      }
    }
  }
  
  return(results)
}

# ============================================================================
# 辅助函数 (Helper Functions)
# ============================================================================

#' 快速基因可视化函数 (Quick gene visualization function)
#' 
#' @description 简化版的基因可视化函数，适用于快速分析
#' 
quick_gene_viz <- function(seurat_obj, genes, group_by = "seurat_clusters") {
  # 基本验证
  if(missing(genes) || length(genes) == 0) {
    stop("必须提供基因列表")
  }
  
  # 过滤存在的基因
  available_genes <- intersect(genes, rownames(seurat_obj))
  
  if(length(available_genes) == 0) {
    stop("没有找到有效的基因")
  }
  
  # 设置分组
  Idents(seurat_obj) <- group_by
  
  # 生成图形
  plots <- list()
  
  # 点图
  plots$dotplot <- DotPlot(seurat_obj, features = available_genes) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 降维图
  plots$dimplot <- DimPlot(seurat_obj, group.by = group_by)
  
  # 热图
  if(length(available_genes) <= 50) {
    plots$heatmap <- DoHeatmap(seurat_obj, features = available_genes)
  }
  
  return(plots)
}

#' 基因功能评分函数 (Gene function scoring function)
#' 
calculate_function_scores <- function(seurat_obj, gene_list, method = "AddModuleScore") {
  scores <- list()
  
  for(func_name in names(gene_list)) {
    func_genes <- gene_list[[func_name]]
    func_genes <- intersect(func_genes, rownames(seurat_obj))
    
    if(length(func_genes) > 0) {
      if(method == "AddModuleScore") {
        temp_obj <- AddModuleScore(seurat_obj, features = list(func_genes), 
                                   name = paste0(func_name, "_"), verbose = FALSE)
        scores[[func_name]] <- temp_obj@meta.data[, paste0(func_name, "_1")]
      } else if(method == "mean") {
        expr_data <- FetchData(seurat_obj, vars = func_genes)
        scores[[func_name]] <- rowMeans(expr_data, na.rm = TRUE)
      }
    }
  }
  
  return(scores)
}
#' 多基因UMAP表达图和小提琴图生成函数
#' Multi-Gene UMAP Expression and Violin Plot Generator
#' 
#' @description 
#' 生成类似Nature风格的多基因表达可视化，包括UMAP特征图和小提琴图
#' Generate Nature-style multi-gene expression visualization with UMAP feature plots and violin plots
#' 
#' @param seurat_obj Seurat对象 (Seurat object)
#' @param gene_list 基因向量或命名列表 (vector of genes or named list of gene categories)
#' @param group_by 分组变量 (grouping variable)，默认为"seurat_clusters"
#' @param reduction 降维方法 (reduction method)，默认为"umap"
#' @param output_dir 输出目录 (output directory)，默认为"Multi_Gene_Visualization"
#' @param prefix 文件前缀 (file prefix)，默认为"MultiGene"
#' @param ncol_feature UMAP特征图列数 (number of columns for feature plots)，默认为4
#' @param ncol_violin 小提琴图列数 (number of columns for violin plots)，默认为3
#' @param pt_size 点大小 (point size)，默认为0.1
#' @param feature_width UMAP图宽度 (feature plot width)，默认为16
#' @param feature_height UMAP图高度 (feature plot height)，默认为12
#' @param violin_width 小提琴图宽度 (violin plot width)，默认为18
#' @param violin_height 小提琴图高度 (violin plot height)，默认为12
#' @param color_scheme 颜色方案 (color scheme)，默认为"viridis"
#' @param max_genes_per_plot 每页最大基因数 (max genes per plot)，默认为16
#' @param sample_cells 是否对细胞进行抽样 (whether to sample cells)，默认为TRUE
#' @param max_cells 最大细胞数 (maximum number of cells)，默认为20000
#' @param verbose 是否显示详细信息 (verbose output)，默认为TRUE
#' 
#' @return 返回包含所有图形对象和文件路径的列表 (list containing all plot objects and file paths)
#' 
#' @examples
#' # 基本用法 (Basic usage)
#' epithelial_markers <- c("EPCAM", "KRT5", "TP63", "MUC5AC", "FOXJ1", "PIFO")
#' results <- generate_multi_gene_visualization(
#'   seurat_obj = your_seurat_object,
#'   gene_list = epithelial_markers,
#'   group_by = "cell_type"
#' )
#' 
#' # 分类基因用法 (Categorized genes usage)
#' gene_categories <- list(
#'   "基底细胞标记" = c("KRT5", "TP63", "KRT14"),
#'   "纤毛细胞标记" = c("FOXJ1", "PIFO", "RSPH1"),
#'   "分泌细胞标记" = c("MUC5AC", "SCGB1A1", "SPDEF")
#' )
#' results <- generate_multi_gene_visualization(
#'   seurat_obj = your_seurat_object,
#'   gene_list = gene_categories
#' )

#' 多基因UMAP表达图和小提琴图生成函数
#' Multi-Gene UMAP Expression and Violin Plot Generator
#' 
#' @description 
#' 生成类似Nature风格的多基因表达可视化，包括UMAP特征图和小提琴图
#' Generate Nature-style multi-gene expression visualization with UMAP feature plots and violin plots
#' 
#' @param seurat_obj Seurat对象 (Seurat object)
#' @param gene_list 基因向量或命名列表 (vector of genes or named list of gene categories)
#' @param group_by 分组变量 (grouping variable)，默认为"seurat_clusters"
#' @param reduction 降维方法 (reduction method)，默认为"umap"
#' @param output_dir 输出目录 (output directory)，默认为"Multi_Gene_Visualization"
#' @param prefix 文件前缀 (file prefix)，默认为"MultiGene"
#' @param ncol_feature UMAP特征图列数 (number of columns for feature plots)，默认为4
#' @param ncol_violin 小提琴图列数 (number of columns for violin plots)，默认为3
#' @param pt_size 点大小 (point size)，默认为0.1
#' @param feature_width UMAP图宽度 (feature plot width)，默认为16
#' @param feature_height UMAP图高度 (feature plot height)，默认为12
#' @param violin_width 小提琴图宽度 (violin plot width)，默认为18
#' @param violin_height 小提琴图高度 (violin plot height)，默认为12
#' @param color_scheme 颜色方案 (color scheme)，默认为"viridis"
#' @param max_genes_per_plot 每页最大基因数 (max genes per plot)，默认为16
#' @param sample_cells 是否对细胞进行抽样 (whether to sample cells)，默认为TRUE
#' @param max_cells 最大细胞数 (maximum number of cells)，默认为20000
#' @param combine_all_categories 是否将所有类别合并到一个PDF (combine all categories in one PDF)，默认为FALSE
#' @param verbose 是否显示详细信息 (verbose output)，默认为TRUE
#' 
#' @return 返回包含所有图形对象和文件路径的列表 (list containing all plot objects and file paths)
#' 
#' @examples
#' # 基本用法 (Basic usage)
#' epithelial_markers <- c("EPCAM", "KRT5", "TP63", "MUC5AC", "FOXJ1", "PIFO")
#' results <- generate_multi_gene_visualization(
#'   seurat_obj = your_seurat_object,
#'   gene_list = epithelial_markers,
#'   group_by = "cell_type"
#' )
#' 
#' # 分类基因用法 (Categorized genes usage)
#' gene_categories <- list(
#'   "基底细胞标记" = c("KRT5", "TP63", "KRT14"),
#'   "纤毛细胞标记" = c("FOXJ1", "PIFO", "RSPH1"),
#'   "分泌细胞标记" = c("MUC5AC", "SCGB1A1", "SPDEF")
#' )
#' results <- generate_multi_gene_visualization(
#'   seurat_obj = your_seurat_object,
#'   gene_list = gene_categories,
#'   combine_all_categories = TRUE  # 将所有类别合并到一个PDF中
#' )

generate_multi_gene_visualization <- function(
    seurat_obj,
    gene_list,
    group_by = "seurat_clusters",
    reduction = "umap",
    output_dir = "Multi_Gene_Visualization",
    prefix = "MultiGene",
    ncol_feature = 4,
    ncol_violin = 3,
    pt_size = 0.1,
    feature_width = 16,
    feature_height = 12,
    violin_width = 18,
    violin_height = 12,
    color_scheme = "viridis",
    max_genes_per_plot = 16,
    sample_cells = TRUE,
    max_cells = 20000,
    combine_all_categories = FALSE,  # 新增：是否将所有类别合并到一个PDF
    verbose = TRUE
) {
  
  # ============================================================================
  # 1. 环境准备和参数验证 (Environment Setup and Parameter Validation)
  # ============================================================================
  
  if(verbose) message("开始多基因可视化分析...")
  
  # 加载必要的R包 (Load required packages)
  required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "viridis", "RColorBrewer")
  
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("请安装R包:", pkg))
    }
  }
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(viridis)
    library(RColorBrewer)
  })
  
  # 验证Seurat对象 (Validate Seurat object)
  if(!inherits(seurat_obj, "Seurat")) {
    stop("必须提供有效的Seurat对象")
  }
  
  # 验证分组变量 (Validate grouping variable)
  if(!group_by %in% colnames(seurat_obj@meta.data)) {
    warning(paste("分组变量", group_by, "不存在，使用默认聚类"))
    group_by <- "seurat_clusters"
  }
  
  # 验证降维结果 (Validate reduction)
  if(!reduction %in% names(seurat_obj@reductions)) {
    warning(paste("降维结果", reduction, "不存在，使用可用的第一个降维"))
    reduction <- names(seurat_obj@reductions)[1]
  }
  
  # 创建输出目录 (Create output directory)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # ============================================================================
  # 2. 数据预处理 (Data Preprocessing)
  # ============================================================================
  
  if(verbose) message("处理基因列表...")
  
  # 处理基因列表 (Process gene list)
  if(is.list(gene_list) && !is.null(names(gene_list))) {
    # 命名列表：分类基因 (Named list: categorized genes)
    all_genes <- unique(unlist(gene_list))
    gene_categories <- gene_list
    has_categories <- TRUE
  } else {
    # 简单向量：所有基因为一类 (Simple vector: all genes in one category)
    all_genes <- unique(as.character(gene_list))
    gene_categories <- list("Selected_Genes" = all_genes)
    has_categories <- FALSE
  }
  
  # 检查基因在数据中的存在情况 (Check gene availability)
  genes_available <- all_genes[all_genes %in% rownames(seurat_obj)]
  genes_missing <- all_genes[!all_genes %in% rownames(seurat_obj)]
  
  if(verbose) {
    message(paste("总基因数量:", length(all_genes)))
    message(paste("可用基因数量:", length(genes_available)))
    if(length(genes_missing) > 0) {
      message(paste("缺失基因:", paste(genes_missing, collapse = ", ")))
    }
  }
  
  if(length(genes_available) == 0) {
    stop("没有找到任何有效的基因用于分析")
  }
  
  # 更新基因分类，只保留可用基因 (Update gene categories with available genes)
  gene_categories_filtered <- lapply(gene_categories, function(genes) {
    intersect(genes, genes_available)
  })
  gene_categories_filtered <- gene_categories_filtered[sapply(gene_categories_filtered, length) > 0]
  
  # 细胞抽样 (Cell sampling)
  if(sample_cells && ncol(seurat_obj) > max_cells) {
    if(verbose) message(paste("对细胞进行抽样，从", ncol(seurat_obj), "个细胞中抽取", max_cells, "个"))
    set.seed(42)
    sampled_cells <- sample(colnames(seurat_obj), max_cells)
    seurat_subset <- subset(seurat_obj, cells = sampled_cells)
  } else {
    seurat_subset <- seurat_obj
  }
  
  # 设置分组标识 (Set grouping identity)
  Idents(seurat_subset) <- group_by
  
  # 初始化结果列表 (Initialize result list)
  results <- list(
    parameters = list(
      genes_total = length(all_genes),
      genes_available = genes_available,
      genes_missing = genes_missing,
      gene_categories = gene_categories_filtered,
      group_by = group_by,
      reduction = reduction,
      timestamp = timestamp
    ),
    plots = list(),
    files = list()
  )
  
  # ============================================================================
  # 3. 生成UMAP特征图 (Generate UMAP Feature Plots)
  # ============================================================================
  
  if(verbose) message("生成UMAP特征图...")
  
  # 根据颜色方案设置调色板 (Set color palette based on color scheme)
  if(color_scheme == "viridis") {
    color_palette <- scale_color_viridis_c(option = "plasma", direction = 1)
  } else if(color_scheme == "blues") {
    color_palette <- scale_color_gradient(low = "lightgrey", high = "darkblue")
  } else if(color_scheme == "reds") {
    color_palette <- scale_color_gradient(low = "lightgrey", high = "darkred")
  } else {
    color_palette <- scale_color_viridis_c()
  }
  
  # 按类别生成特征图 (Generate feature plots by category)
  feature_plots_by_category <- list()
  
  for(category_name in names(gene_categories_filtered)) {
    category_genes <- gene_categories_filtered[[category_name]]
    
    if(length(category_genes) == 0) next
    
    if(verbose) message(paste("处理类别:", category_name, "包含", length(category_genes), "个基因"))
    
    # 将基因分批处理 (Process genes in batches)
    gene_batches <- split(category_genes, ceiling(seq_along(category_genes) / max_genes_per_plot))
    
    category_plots <- list()
    
    # 创建该类别的特征图PDF文件 (Create feature plot PDF for this category)
    feature_pdf_file <- file.path(
      output_dir, 
      paste0(prefix, "_", gsub("[^A-Za-z0-9_]", "_", category_name), 
             "_feature_plots_", timestamp, ".pdf")
    )
    
    pdf(feature_pdf_file, width = feature_width, height = feature_height)
    
    for(batch_idx in seq_along(gene_batches)) {
      batch_genes <- gene_batches[[batch_idx]]
      
      if(verbose) message(paste("  处理批次", batch_idx, ":", paste(batch_genes, collapse = ", ")))
      
      # 生成每个基因的特征图 (Generate feature plot for each gene)
      individual_plots <- list()
      
      for(gene in batch_genes) {
        tryCatch({
          p <- FeaturePlot(
            seurat_subset,
            features = gene,
            reduction = reduction,
            pt.size = pt_size,
            raster = TRUE,
            order = TRUE
          ) + 
            color_palette +
            ggtitle(gene) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 9)
            ) +
            guides(color = guide_colorbar(title = "Expression"))
          
          individual_plots[[gene]] <- p
          
        }, error = function(e) {
          warning(paste("基因", gene, "特征图生成失败:", e$message))
        })
      }
      
      if(length(individual_plots) > 0) {
        # 组合图形 (Combine plots)
        combined_plot <- wrap_plots(individual_plots, ncol = ncol_feature)
        
        # 添加总标题 (Add overall title)
        if(has_categories) {
          final_title <- paste(category_name, "- Page", batch_idx)
        } else {
          final_title <- paste("Gene Expression - Page", batch_idx)
        }
        
        combined_plot_with_title <- combined_plot + 
          plot_annotation(
            title = final_title,
            theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
          )
        
        # 打印到PDF页面 (Print to PDF page)
        print(combined_plot_with_title)
        
        category_plots[[paste0("batch_", batch_idx)]] <- combined_plot_with_title
      }
    }
    
    dev.off()
    
    feature_plots_by_category[[category_name]] <- category_plots
    results$files[[paste0(category_name, "_feature_plots")]] <- feature_pdf_file
    
    if(verbose) message(paste("已保存特征图PDF:", basename(feature_pdf_file)))
  }
  
  results$plots$feature_plots <- feature_plots_by_category
  
  # ============================================================================
  # 4. 生成小提琴图 (Generate Violin Plots)
  # ============================================================================
  
  if(verbose) message("生成小提琴图...")
  
  violin_plots_by_category <- list()
  
  for(category_name in names(gene_categories_filtered)) {
    category_genes <- gene_categories_filtered[[category_name]]
    
    if(length(category_genes) == 0) next
    
    if(verbose) message(paste("生成", category_name, "的小提琴图"))
    
    # 将基因分批处理 (Process genes in batches)
    gene_batches <- split(category_genes, ceiling(seq_along(category_genes) / max_genes_per_plot))
    
    category_violin_plots <- list()
    
    # 创建该类别的小提琴图PDF文件 (Create violin plot PDF for this category)
    violin_pdf_file <- file.path(
      output_dir, 
      paste0(prefix, "_", gsub("[^A-Za-z0-9_]", "_", category_name), 
             "_violin_plots_", timestamp, ".pdf")
    )
    
    pdf(violin_pdf_file, width = violin_width, height = violin_height)
    
    for(batch_idx in seq_along(gene_batches)) {
      batch_genes <- gene_batches[[batch_idx]]
      
      tryCatch({
        # 方法1: 传统小提琴图 (Traditional violin plots)
        p_violin_traditional <- VlnPlot(
          seurat_subset,
          features = batch_genes,
          group.by = group_by,
          pt.size = 0,  # 不显示点
          ncol = ncol_violin,
          combine = TRUE
        ) + 
          plot_annotation(
            title = paste(category_name, "Expression Distribution - Page", batch_idx),
            theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
          )
        
        # 打印传统小提琴图到PDF页面 (Print traditional violin plot to PDF page)
        print(p_violin_traditional)
        
        # 方法2: 横向热图风格小提琴图 (Horizontal heatmap-style violin plot)
        if(length(batch_genes) <= 10) {  # 只有基因数量适中时才生成横向图
          
          # 准备数据 (Prepare data)
          plot_data_list <- list()
          
          for(gene in batch_genes) {
            gene_data <- FetchData(seurat_subset, vars = c(gene, group_by))
            colnames(gene_data) <- c("Expression", "Group")
            gene_data$Gene <- gene
            plot_data_list[[gene]] <- gene_data
          }
          
          combined_data <- do.call(rbind, plot_data_list)
          
          # 计算平均表达用于颜色映射 (Calculate average expression for color mapping)
          avg_expr <- combined_data %>%
            group_by(Gene, Group) %>%
            summarise(
              avg_expr = mean(Expression, na.rm = TRUE),
              .groups = 'drop'
            )
          
          # 创建横向小提琴图 (Create horizontal violin plot)
          p_violin_horizontal <- ggplot(combined_data, aes(x = Gene, y = Expression, fill = Group)) +
            geom_violin(scale = "width", trim = FALSE, alpha = 0.8, color = "white", size = 0.3) +
            stat_summary(fun = mean, geom = "point", size = 1, color = "black", alpha = 0.8) +
            scale_fill_viridis_d(option = "plasma", name = "Cell Type") +
            coord_flip() +  # 翻转坐标使基因在Y轴
            theme_minimal() +
            theme(
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 10),
              legend.position = "bottom",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              panel.grid.major = element_line(size = 0.3, alpha = 0.5),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
            ) +
            labs(
              title = paste(category_name, "Expression Profile - Horizontal View"),
              x = "Genes",
              y = "Expression Level"
            )
          
          print(p_violin_horizontal)
          
          # 方法3: 热图风格的表达可视化 (Heatmap-style expression visualization)
          # 创建类似目标图的热图样式
          heatmap_data <- combined_data %>%
            group_by(Gene, Group) %>%
            summarise(
              avg_expr = mean(Expression, na.rm = TRUE),
              pct_expr = sum(Expression > 0) / n() * 100,
              .groups = 'drop'
            )
          
          # 使用ggplot创建热图风格的小提琴图组合
          p_heatmap_violin <- ggplot(combined_data, aes(x = Gene, y = Group)) +
            # 添加小提琴图背景
            geom_violin(aes(fill = Group), alpha = 0.3, scale = "width", trim = FALSE) +
            # 添加热图颜色映射
            geom_tile(data = heatmap_data, aes(fill = avg_expr), alpha = 0.8, height = 0.8) +
            scale_fill_viridis_c(option = "plasma", name = "Average\nExpression") +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 12),
              legend.position = "right",
              legend.title = element_text(size = 10),
              panel.grid = element_blank(),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
            ) +
            labs(
              title = paste(category_name, "Expression Heatmap with Violin Overlay"),
              x = "Genes",
              y = "Cell Types"
            )
          
          print(p_heatmap_violin)
        }
        
        category_violin_plots[[paste0("batch_", batch_idx)]] <- p_violin_traditional
        
        if(verbose) message(paste("  已添加小提琴图页面", batch_idx))
        
      }, error = function(e) {
        warning(paste("类别", category_name, "批次", batch_idx, "小提琴图生成失败:", e$message))
      })
    }
    
    dev.off()
    
    violin_plots_by_category[[category_name]] <- category_violin_plots
    results$files[[paste0(category_name, "_violin_plots")]] <- violin_pdf_file
    
    if(verbose) message(paste("已保存小提琴图PDF:", basename(violin_pdf_file)))
  }
  
  results$plots$violin_plots <- violin_plots_by_category
  
  # ============================================================================
  # 4.5. 生成合并PDF (Generate Combined PDFs if requested)
  # ============================================================================
  
  if(combine_all_categories && length(gene_categories_filtered) > 1) {
    if(verbose) message("生成合并的PDF文件...")
    
    # 生成合并的特征图PDF (Generate combined feature plots PDF)
    combined_feature_pdf <- file.path(output_dir, paste0(prefix, "_ALL_feature_plots_", timestamp, ".pdf"))
    pdf(combined_feature_pdf, width = feature_width, height = feature_height)
    
    for(category_name in names(feature_plots_by_category)) {
      category_plots <- feature_plots_by_category[[category_name]]
      for(batch_name in names(category_plots)) {
        print(category_plots[[batch_name]])
      }
    }
    dev.off()
    
    results$files$combined_feature_plots <- combined_feature_pdf
    if(verbose) message(paste("已保存合并特征图PDF:", basename(combined_feature_pdf)))
    
    # 生成合并的小提琴图PDF (Generate combined violin plots PDF)  
    combined_violin_pdf <- file.path(output_dir, paste0(prefix, "_ALL_violin_plots_", timestamp, ".pdf"))
    pdf(combined_violin_pdf, width = violin_width, height = violin_height)
    
    for(category_name in names(violin_plots_by_category)) {
      category_plots <- violin_plots_by_category[[category_name]]
      for(batch_name in names(category_plots)) {
        print(category_plots[[batch_name]])
      }
    }
    dev.off()
    
    results$files$combined_violin_plots <- combined_violin_pdf
    if(verbose) message(paste("已保存合并小提琴图PDF:", basename(combined_violin_pdf)))
  }
  
  # ============================================================================
  # 5. 生成组合热图 (Generate Combined Heatmap)
  # ============================================================================
  
  if(verbose) message("生成表达热图...")
  
  tryCatch({
    # 使用所有可用基因生成热图 (Generate heatmap with all available genes)
    if(length(genes_available) > 0) {
      
      # 确保基因已标准化 (Ensure genes are scaled)
      genes_to_scale <- genes_available[!genes_available %in% rownames(seurat_subset[["RNA"]]@scale.data)]
      if(length(genes_to_scale) > 0) {
        if(verbose) message("标准化缺失的基因...")
        seurat_subset <- ScaleData(seurat_subset, features = genes_to_scale, verbose = FALSE)
      }
      
      # 生成传统Seurat热图 (Generate traditional Seurat heatmap)
      p_heatmap <- DoHeatmap(
        seurat_subset,
        features = genes_available,
        group.by = group_by,
        raster = TRUE
      ) + 
        scale_fill_viridis_c(option = "plasma") +
        ggtitle("Gene Expression Heatmap") +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 8)
        )
      
      # 保存传统热图 (Save traditional heatmap)
      heatmap_filename <- file.path(
        output_dir, 
        paste0(prefix, "_expression_heatmap_", timestamp, ".pdf")
      )
      
      # 计算热图尺寸 (Calculate heatmap dimensions)
      heatmap_width <- min(20, max(12, length(genes_available) * 0.2 + 8))
      heatmap_height <- min(16, max(10, 12))
      
      ggsave(
        heatmap_filename, 
        p_heatmap, 
        width = heatmap_width, 
        height = heatmap_height,
        device = "pdf"
      )
      
      results$plots$heatmap <- p_heatmap
      results$files$heatmap <- heatmap_filename
      
      if(verbose) message(paste("已保存传统热图:", basename(heatmap_filename)))
      
      # 生成横向热图风格的可视化 (Generate horizontal heatmap-style visualization)
      if(verbose) message("生成横向热图风格可视化...")
      
      # 为每个类别生成横向热图 (Generate horizontal heatmap for each category)
      horizontal_heatmap_files <- list()
      
      for(category_name in names(gene_categories_filtered)) {
        category_genes <- gene_categories_filtered[[category_name]]
        
        if(length(category_genes) > 0 && length(category_genes) <= 15) {  # 限制基因数量以保持可读性
          
          # 生成横向热图风格图 (Generate horizontal heatmap-style plot)
          p_horizontal <- generate_heatmap_violin_plot(
            seurat_obj = seurat_subset,
            genes = category_genes,
            group_by = group_by,
            title = paste(category_name, "Expression Profile"),
            color_scheme = color_scheme,
            show_points = TRUE
          )
          
          # 保存横向热图 (Save horizontal heatmap)
          horizontal_filename <- file.path(
            output_dir, 
            paste0(prefix, "_", gsub("[^A-Za-z0-9_]", "_", category_name), 
                   "_horizontal_heatmap_", timestamp, ".pdf")
          )
          
          ggsave(
            horizontal_filename, 
            p_horizontal, 
            width = max(10, length(category_genes) * 1.2 + 4), 
            height = max(6, length(unique(seurat_subset@meta.data[[group_by]])) * 0.8 + 2),
            device = "pdf"
          )
          
          horizontal_heatmap_files[[category_name]] <- horizontal_filename
          results$files[[paste0(category_name, "_horizontal_heatmap")]] <- horizontal_filename
          
          if(verbose) message(paste("已保存横向热图:", basename(horizontal_filename)))
        }
      }
      
      # 生成综合的平均表达热图 (Generate comprehensive average expression heatmap)
      if(verbose) message("生成平均表达热图...")
      
      # 计算平均表达 (Calculate average expression)
      avg_expr <- AverageExpression(
        seurat_subset, 
        features = genes_available,
        group.by = group_by,
        assays = "RNA",
        slot = "data"
      )
      
      avg_expr_matrix <- avg_expr[["RNA"]]
      
      # 标准化用于热图 (Normalize for heatmap)
      avg_expr_scaled <- t(scale(t(avg_expr_matrix)))
      
      # 使用pheatmap生成高质量热图 (Generate high-quality heatmap with pheatmap)
      library(pheatmap)
      
      # 准备基因注释 (Prepare gene annotations)
      if(has_categories) {
        gene_annotation <- data.frame(
          Category = NA,
          stringsAsFactors = FALSE,
          row.names = rownames(avg_expr_scaled)
        )
        
        for(cat_name in names(gene_categories_filtered)) {
          cat_genes <- gene_categories_filtered[[cat_name]]
          gene_annotation[cat_genes, "Category"] <- cat_name
        }
        
        # 设置类别颜色 (Set category colors)
        n_categories <- length(unique(gene_annotation$Category[!is.na(gene_annotation$Category)]))
        category_colors <- setNames(
          rainbow(n_categories), 
          unique(gene_annotation$Category[!is.na(gene_annotation$Category)])
        )
        
        anno_colors <- list(Category = category_colors)
      } else {
        gene_annotation <- NULL
        anno_colors <- NULL
      }
      
      # 保存平均表达热图 (Save average expression heatmap)
      pheatmap_filename <- file.path(
        output_dir, 
        paste0(prefix, "_average_expression_pheatmap_", timestamp, ".pdf")
      )
      
      pdf(pheatmap_filename, width = max(10, ncol(avg_expr_scaled) * 0.8 + 4), 
          height = max(8, nrow(avg_expr_scaled) * 0.3 + 3))
      
      pheatmap(
        avg_expr_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        annotation_row = gene_annotation,
        annotation_colors = anno_colors,
        fontsize_row = max(6, min(10, 120/nrow(avg_expr_scaled))),
        fontsize_col = 10,
        color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        main = "Average Gene Expression Heatmap",
        angle_col = 45,
        border_color = NA
      )
      
      dev.off()
      
      results$files$pheatmap <- pheatmap_filename
      
      if(verbose) message(paste("已保存平均表达热图:", basename(pheatmap_filename)))
    }
    
  }, error = function(e) {
    warning(paste("热图生成失败:", e$message))
  })
  
  # ============================================================================
  # 6. 生成分析报告 (Generate Analysis Report)
  # ============================================================================
  
  if(verbose) message("生成分析总结...")
  
  # 保存基因信息 (Save gene information)
  gene_info_file <- file.path(output_dir, paste0(prefix, "_gene_info_", timestamp, ".txt"))
  
  sink(gene_info_file)
  cat("多基因可视化分析报告\n")
  cat("===================\n\n")
  cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("细胞数量:", ncol(seurat_subset), "\n")
  cat("分组变量:", group_by, "\n")
  cat("降维方法:", reduction, "\n\n")
  
  cat("基因分析统计:\n")
  cat("总基因数量:", length(all_genes), "\n")
  cat("可用基因数量:", length(genes_available), "\n")
  cat("缺失基因数量:", length(genes_missing), "\n\n")
  
  if(length(genes_missing) > 0) {
    cat("缺失基因列表:\n")
    cat(paste(genes_missing, collapse = ", "), "\n\n")
  }
  
  cat("基因分类信息:\n")
  for(category in names(gene_categories_filtered)) {
    genes_in_category <- gene_categories_filtered[[category]]
    cat(sprintf("%-25s: %d 个基因 (%s)\n", 
                category, 
                length(genes_in_category),
                paste(genes_in_category, collapse = ", ")))
  }
  
  cat("\n生成的文件:\n")
  for(file_type in names(results$files)) {
    cat(sprintf("%-30s: %s\n", file_type, basename(results$files[[file_type]])))
  }
  sink()
  
  results$files$gene_info <- gene_info_file
  
  # ============================================================================
  # 7. 完成分析 (Complete Analysis)
  # ============================================================================
  
  if(verbose) {
    message("分析完成！")
    message(paste("结果保存在:", output_dir))
    message(paste("共生成", length(results$files), "个文件"))
  }
  
  return(results)
}

#' 生成横向热图风格小提琴图 (Generate horizontal heatmap-style violin plot)
#' 
#' @description 专门生成类似目标图形的横向热图小提琴图组合
generate_heatmap_violin_plot <- function(
    seurat_obj,
    genes,
    group_by = "seurat_clusters",
    title = "Gene Expression Profile",
    width = 12,
    height = 6,
    color_scheme = "plasma",
    show_points = TRUE,
    point_size = 0.5,
    violin_alpha = 0.7,
    tile_alpha = 0.8
) {
  
  # 验证基因 (Validate genes)
  available_genes <- intersect(genes, rownames(seurat_obj))
  if(length(available_genes) == 0) {
    stop("没有找到有效的基因")
  }
  
  # 准备数据 (Prepare data)
  plot_data_list <- list()
  
  for(gene in available_genes) {
    gene_data <- FetchData(seurat_obj, vars = c(gene, group_by))
    colnames(gene_data) <- c("Expression", "Group")
    gene_data$Gene <- gene
    plot_data_list[[gene]] <- gene_data
  }
  
  combined_data <- do.call(rbind, plot_data_list)
  combined_data$Gene <- factor(combined_data$Gene, levels = available_genes)
  
  # 计算统计数据 (Calculate statistics)
  stats_data <- combined_data %>%
    group_by(Gene, Group) %>%
    summarise(
      avg_expr = mean(Expression, na.rm = TRUE),
      max_expr = max(Expression, na.rm = TRUE),
      pct_expr = sum(Expression > 0) / n() * 100,
      .groups = 'drop'
    )
  
  # 创建类似目标图的可视化 (Create target-style visualization)
  p <- ggplot(combined_data, aes(x = Gene, y = Group)) +
    # 添加小提琴形状 (Add violin shapes)
    geom_violin(
      aes(fill = Group), 
      alpha = violin_alpha, 
      scale = "width", 
      trim = FALSE,
      color = "white",
      size = 0.3
    ) +
    # 添加平均表达热图 (Add average expression heatmap)
    geom_tile(
      data = stats_data, 
      aes(fill = avg_expr), 
      alpha = tile_alpha,
      height = 0.8,
      color = "white",
      size = 0.5
    ) +
    # 可选：添加平均值点 (Optional: add mean points)
    {if(show_points) 
      geom_point(
        data = stats_data,
        aes(size = pct_expr),
        color = "black",
        alpha = 0.8
      )
    } +
    # 颜色映射 (Color mapping)
    {if(color_scheme == "plasma") 
      scale_fill_viridis_c(option = "plasma", name = "Average\nExpression")
      else if(color_scheme == "reds")
        scale_fill_gradient(low = "white", high = "darkred", name = "Average\nExpression")
      else 
        scale_fill_viridis_c(name = "Average\nExpression")
    } +
    # 可选：点大小映射 (Optional: point size mapping)
    {if(show_points)
      scale_size_continuous(
        name = "% Expressed",
        range = c(0.5, 3),
        guide = guide_legend(override.aes = list(color = "black", alpha = 1))
      )
    } +
    # 主题设置 (Theme settings)
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "italic"),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.key.size = unit(0.8, "cm")
    ) +
    labs(
      title = title,
      x = "",
      y = ""
    )
  
  return(p)
}

# ============================================================================
# 辅助函数 (Helper Functions)
# ============================================================================

#' 快速多基因可视化函数 (Quick multi-gene visualization)
#' 
#' @description 简化版多基因可视化，适用于快速预览
quick_multi_gene_viz <- function(seurat_obj, genes, ncol = 4, pt_size = 0.1) {
  
  # 过滤可用基因 (Filter available genes)
  available_genes <- intersect(genes, rownames(seurat_obj))
  
  if(length(available_genes) == 0) {
    stop("没有找到有效的基因")
  }
  
  # 生成特征图 (Generate feature plots)
  feature_plots <- list()
  for(gene in available_genes) {
    p <- FeaturePlot(seurat_obj, features = gene, pt.size = pt_size, raster = TRUE) +
      ggtitle(gene) + theme_minimal()
    feature_plots[[gene]] <- p
  }
  
  # 组合图形 (Combine plots)
  combined_plot <- wrap_plots(feature_plots, ncol = ncol)
  
  return(combined_plot)
}

#' 基因表达统计函数 (Gene expression statistics)
#' 
calculate_gene_stats <- function(seurat_obj, genes, group_by = "seurat_clusters") {
  
  # 过滤可用基因 (Filter available genes)
  available_genes <- intersect(genes, rownames(seurat_obj))
  
  # 计算统计信息 (Calculate statistics)
  stats_list <- list()
  
  for(gene in available_genes) {
    expr_data <- FetchData(seurat_obj, vars = c(gene, group_by))
    
    stats <- expr_data %>%
      group_by(!!sym(group_by)) %>%
      summarise(
        mean_expr = mean(!!sym(gene), na.rm = TRUE),
        median_expr = median(!!sym(gene), na.rm = TRUE),
        pct_expressed = sum(!!sym(gene) > 0) / n() * 100,
        .groups = 'drop'
      ) %>%
      mutate(gene = gene)
    
    stats_list[[gene]] <- stats
  }
  
  # 合并所有统计结果 (Combine all statistics)
  all_stats <- do.call(rbind, stats_list)
  
  return(all_stats)
}

# ============================================================================
# 使用示例 (Usage Examples)
# ============================================================================

cat("
使用示例 (Usage Examples):
=========================

# 1. 基本用法 - 每个类别单独PDF文件，不同batch为不同页面
epithelial_markers <- c('EPCAM', 'KRT5', 'TP63', 'MUC5AC', 'MUC5B', 
                        'FOXJ1', 'PIFO', 'SCGB1A1', 'SCGB3A1', 'LYZ', 'LTF')

results <- generate_multi_gene_visualization(
  seurat_obj = your_seurat_object,
  gene_list = epithelial_markers,
  group_by = 'cell_type',
  output_dir = 'Epithelial_Multi_Gene_Viz',
  prefix = 'Epithelial'
)

# 2. 生成目标样式的横向热图小提琴图 (推荐)
target_genes <- c('FOXJ1', 'DYNC2H1', 'KIF3A', 'CDFHR3', 'SCGB1A1', 
                  'SCGB3A1', 'CST1', 'SERPINB3', 'SAA1')

gene_categories <- list(
  '纤毛相关基因' = c('FOXJ1', 'DYNC2H1', 'KIF3A'),
  '分泌相关基因' = c('SCGB1A1', 'SCGB3A1'),
  '防御相关基因' = c('CST1', 'SERPINB3', 'SAA1')
)

results <- generate_multi_gene_visualization(
  seurat_obj = your_seurat_object,
  gene_list = gene_categories,
  group_by = 'cell_type',  # 如：'HC', 'CRSsNP', 'neCRSwNP', 'eCRSwNP'
  max_genes_per_plot = 9,  # 每页最多9个基因，便于横向显示
  color_scheme = 'plasma'  # 类似目标图的颜色方案
)

# 3. 直接生成单个横向热图风格图
target_plot <- generate_heatmap_violin_plot(
  seurat_obj = your_seurat_object,
  genes = target_genes,
  group_by = 'cell_type',
  title = 'Gene Expression Profile',
  color_scheme = 'plasma',
  show_points = TRUE,      # 显示表达百分比点
  width = 12,
  height = 6
)

# 保存单个横向图
ggsave('target_style_plot.pdf', target_plot, width = 12, height = 6)

# 4. 分类基因用法 - 每个类别一个PDF，页面按batch分组
gene_categories_detailed <- list(
  '基底细胞标记' = c('KRT5', 'TP63', 'KRT14'),
  '纤毛细胞标记' = c('FOXJ1', 'PIFO', 'RSPH1', 'DYNC2H1', 'KIF3A'),
  '分泌细胞标记' = c('MUC5AC', 'MUC5B', 'SCGB1A1', 'SCGB3A1', 'SPDEF'),
  '抗菌防御' = c('LYZ', 'LTF', 'SLPI', 'DMBT1', 'CST1', 'SERPINB3')
)

results <- generate_multi_gene_visualization(
  seurat_obj = your_seurat_object,
  gene_list = gene_categories_detailed,
  group_by = 'Annotation_2',
  ncol_feature = 3,        # 每行3个特征图
  ncol_violin = 2,         # 每行2个小提琴图
  max_genes_per_plot = 12  # 每页最多12个基因
)

# 5. 合并所有类别到单个PDF - 适合综合展示
results <- generate_multi_gene_visualization(
  seurat_obj = your_seurat_object,
  gene_list = gene_categories_detailed,
  combine_all_categories = TRUE,  # 生成合并PDF
  output_dir = 'Combined_Gene_Viz',
  prefix = 'AllCategories'
)

# 6. 大数据集优化用法
results <- generate_multi_gene_visualization(
  seurat_obj = large_seurat_object,
  gene_list = many_genes,
  sample_cells = TRUE,      # 启用细胞抽样
  max_cells = 15000,       # 最多使用15000个细胞
  pt_size = 0.05,          # 减小点大小
  max_genes_per_plot = 8   # 每页减少基因数量
)

# 7. 生成发表级别的图表组合
publication_genes <- list(
  'Ciliated_Markers' = c('FOXJ1', 'PIFO', 'RSPH1'),
  'Secretory_Markers' = c('SCGB1A1', 'SCGB3A1', 'MUC5AC'),
  'Defense_Markers' = c('LYZ', 'LTF', 'SERPINB3')
)

# 生成完整分析
results <- generate_multi_gene_visualization(
  seurat_obj = epithelial_obj,
  gene_list = publication_genes,
  group_by = 'cell_subtype',
  color_scheme = 'plasma',
  feature_width = 15,
  feature_height = 10,
  violin_width = 18,
  violin_height = 8,
  combine_all_categories = TRUE
)

输出文件说明:
============

## 标准输出文件:
- *_{类别名}_feature_plots_*.pdf: 每个类别的UMAP特征图 (多页PDF)
- *_{类别名}_violin_plots_*.pdf: 每个类别的小提琴图 (多页PDF)
  └── 包含3种风格: 传统小提琴图 + 横向小提琴图 + 热图小提琴图组合
  
## 新增热图输出文件:
- *_expression_heatmap_*.pdf: Seurat风格热图
- *_{类别名}_horizontal_heatmap_*.pdf: 横向热图风格图 (类似目标图)
- *_average_expression_pheatmap_*.pdf: 高质量平均表达热图
- *_gene_info_*.txt: 基因分析统计报告

## 合并模式额外文件 (combine_all_categories = TRUE):
- *_ALL_feature_plots_*.pdf: 所有类别特征图合并
- *_ALL_violin_plots_*.pdf: 所有类别小提琴图合并

目标图形样式特点:
===============
✅ 横向布局: 基因在X轴，细胞类型在Y轴
✅ 小提琴形状: 显示表达分布的形状
✅ 颜色映射: 使用viridis plasma配色显示平均表达
✅ 统计信息: 可选择显示表达百分比点
✅ 发表质量: 高分辨率PDF，适合期刊投稿

使用技巧:
========
1. 对于目标样式图，建议每次不超过10个基因以保持可读性
2. 使用plasma颜色方案获得最佳视觉效果  
3. 横向热图最适合展示功能相关的基因集
4. 可以结合传统分析方法获得全面的可视化
5. 大数据集建议先抽样再生成横向图
")
epithelial_gene_functions <- list(
  
  # 1. 通用上皮标记 (General Epithelial Markers)
  "General Epithelial Markers" = c(
    "EPCAM",      # 上皮细胞粘附分子 (Epithelial Cell Adhesion Molecule)
    "CDH1",       # E-钙粘蛋白 (E-cadherin)
    "CLDN1",      # 紧密连接蛋白1 (Claudin-1)
    "OCLN",       # 咬合蛋白 (Occludin)
    "TJP1",       # 紧密连接蛋白1 (ZO-1)
    "DSP",        # 胞间连接斑蛋白 (Desmoplakin)
    "PKP1",       # 斑块糖蛋白1 (Plakophilin-1)
    "KRT18",      # 角蛋白18 (Keratin 18)
    "KRT19",      # 角蛋白19 (Keratin 19)
    "FXYD3",      # FXYD域包含离子转运调节因子3
    "ELF3",       # E74样因子3 (E74-like factor 3)
    "GRHL1",      # 灰度样同源异型结构域1 (Grainyhead-like 1)
    "OVOL1"       # 卵蛋白相关转录因子1 (Ovo-like 1)
  ),
  
  # 2. 基底细胞标记 (Basal Cell Markers)
  "Basal Cell Markers" = c(
    "KRT5",       # 角蛋白5 (Keratin 5)
    "KRT14",      # 角蛋白14 (Keratin 14)
    "TP63",       # 肿瘤蛋白p63 (Tumor protein p63)
    "ITGA6",      # 整合素α6 (Integrin alpha-6)
    "ITGB4",      # 整合素β4 (Integrin beta-4)
    "COL17A1",    # XVII型胶原α1链 (Collagen XVII alpha-1)
    "LAMB3",      # 层粘连蛋白β3 (Laminin beta-3)
    "SERPINB3",   # 丝氨酸蛋白酶抑制剂B3 (Serpin B3)
    "S100A2",     # S100钙结合蛋白A2 (S100 calcium-binding protein A2)
    "IVL",        # 囊膜蛋白 (Involucrin)
    "SPRR1A",     # 小脯氨酸富集区域1A (Small proline-rich protein 1A)
    "DSC3"        # 胞间连接芯糖蛋白3 (Desmocollin-3)
  ),
  
  # 3. 纤毛上皮标记 (Ciliated Epithelial Markers) 
  "Ciliated Epithelial Markers" = c(
    "FOXJ1",      # 叉头盒蛋白J1 (Forkhead box J1)
    "RSPH1",      # 放射状辐条头蛋白1 (Radial spoke head protein 1)
    "PIFO",       # 纤毛形成蛋白 (Pifo homolog)
    "CCDC78",     # 卷曲螺旋域包含蛋白78 (Coiled-coil domain containing 78)
    "CCDC153",    # 卷曲螺旋域包含蛋白153 (Coiled-coil domain containing 153)
    "DNAH5",      # 动力蛋白重链5 (Dynein heavy chain 5)
    "DNAI1",      # 动力蛋白中间链1 (Dynein intermediate chain 1)
    "C9orf24",    # 9号染色体开放阅读框24
    "C20orf85",   # 20号染色体开放阅读框85
    "TPPP3",      # 微管聚合促进蛋白家族成员3 (Tubulin polymerization promoting protein family member 3)
    "SNTN",       # 辛匹蛋白 (Sinectin)
    "CCDC114"     # 卷曲螺旋域包含蛋白114 (Coiled-coil domain containing 114)
  ),
  
  # 4. 分泌细胞标记 (Secretory Cell Markers)
  "Secretory Cell Markers" = c(
    "MUC5AC",     # 粘蛋白5AC (Mucin 5AC)
    "MUC5B",      # 粘蛋白5B (Mucin 5B)
    "SPDEF",      # SAM尖端域ETS因子 (SAM pointed domain ETS factor)
    "SCGB1A1",    # 分泌球蛋白1A1 (Secretoglobin 1A1)
    "SCGB3A1",    # 分泌球蛋白3A1 (Secretoglobin 3A1)
    "SCGB3A2",    # 分泌球蛋白3A2 (Secretoglobin 3A2)
    "TFF1",       # 三叶因子1 (Trefoil factor 1)
    "TFF3",       # 三叶因子3 (Trefoil factor 3)
    "AGR2",       # 前端梯度蛋白2 (Anterior gradient protein 2)
    "LYPD2",      # LY6/PLAUR域包含蛋白2 (LY6/PLAUR domain containing 2)
    "BPIFA1",     # BPI折叠包含家族A成员1 (BPI fold containing family A member 1)
    "PIGR"        # 多聚免疫球蛋白受体 (Polymeric immunoglobulin receptor)
  ),
  
  # 5. 肺泡上皮标记 (Alveolar Epithelial Markers)
  "Alveolar Epithelial Markers" = c(
    "AGER",       # 晚期糖基化终产物受体 (Advanced glycosylation end-product receptor)
    "PDPN",       # 足细胞糖蛋白 (Podoplanin)
    "CLIC5",      # 氯离子通道蛋白5 (Chloride intracellular channel 5)
    "CAV1",       # 胞膜窖蛋白1 (Caveolin-1)
    "RTKN2",      # Rhotekin 2
    "HOPX",       # HOP同源异型盒蛋白 (HOP homeobox)
    "SFTPC",      # 表面活性蛋白C (Surfactant protein C)
    "SFTPB",      # 表面活性蛋白B (Surfactant protein B)
    "SFTPA1",     # 表面活性蛋白A1 (Surfactant protein A1)
    "SFTPA2",     # 表面活性蛋白A2 (Surfactant protein A2)
    "ABCA3",      # ATP结合盒转运蛋白A3 (ATP binding cassette subfamily A member 3)
    "LAMP3"       # 溶酶体相关膜蛋白3 (Lysosomal associated membrane protein 3)
  ),
  
  # 6. 神经内分泌细胞标记 (Neuroendocrine Cell Markers)
  "Neuroendocrine Cell Markers" = c(
    "ASCL1",      # 酰基螺旋环螺旋样1 (Achaete-scute homolog 1)
    "CHGA",       # 嗜铬粒蛋白A (Chromogranin A)
    "CHGB",       # 嗜铬粒蛋白B (Chromogranin B)
    "SYP",        # 突触素 (Synaptophysin)
    "NCAM1",      # 神经细胞粘附分子1 (Neural cell adhesion molecule 1)
    "GRP",        # 胃泌素释放肽 (Gastrin releasing peptide)
    "CALCA",      # 降钙素相关多肽α (Calcitonin related polypeptide alpha)
    "ENO2",       # 烯醇化酶2 (Enolase 2)
    "PCSK1",      # 前蛋白转化酶枯草杆菌蛋白酶/kexin 1型 (Proprotein convertase subtilisin/kexin type 1)
    "INSM1",      # INSM转录抑制因子1 (INSM transcriptional repressor 1)
    "DDC",        # DOPA脱羧酶 (Dopa decarboxylase)
    "UCHL1"       # 泛素C末端水解酶L1 (Ubiquitin C-terminal hydrolase L1)
  ),
  
  # 7. 感觉细胞标记 (Sensory Cell Markers)
  "Sensory Cell Markers" = c(
    "POU2F3",     # POU类2同源异型盒3 (POU class 2 homeobox 3)
    "TRPM5",      # 瞬时受体电位阳离子通道M5 (Transient receptor potential cation channel subfamily M member 5)
    "AVIL",       # 副肌动蛋白样 (Advillin)
    "GNAT3",      # 鸟苷酸结合蛋白α转导蛋白3 (G protein subunit alpha transducin 3)
    "PLCB2",      # 磷脂酶Cβ2 (Phospholipase C beta 2)
    "ASCL2",      # 酰基螺旋环螺旋样2 (Achaete-scute homolog 2)
    "SOX9",       # SRY-box转录因子9 (SRY-box transcription factor 9)
    "HCK",        # 造血细胞激酶 (Hemopoietic cell kinase)
    "PTGS1",      # 前列腺素G/H合酶1 (Prostaglandin G/H synthase 1)
    "IL25",       # 白细胞介素25 (Interleukin 25)
    "SUCNR1",     # 琥珀酸受体1 (Succinate receptor 1)
    "BMX"         # BMX非受体酪氨酸激酶 (BMX non-receptor tyrosine kinase)
  ),
  
  # 8. 离子转运细胞标记 (Ionocyte Markers)
  "Ionocyte Markers" = c(
    "CFTR",       # 囊性纤维化跨膜传导调节因子 (Cystic fibrosis transmembrane conductance regulator)
    "FOXI1",      # 叉头盒蛋白I1 (Forkhead box I1)
    "ASCL3",      # 酰基螺旋环螺旋样3 (Achaete-scute homolog 3)
    "CLCNKB",     # 氯离子通道Kb (Chloride voltage-gated channel Kb)
    "BSND",       # Bartter综合征、耳聋，常染色体隐性遗传 (Bartter syndrome and deafness)
    "ATP6V1B1",   # ATPase H+转运V型1亚基B1 (ATPase H+ transporting V1 subunit B1)
    "ATP6V0A4",   # ATPase H+转运V0亚基a4 (ATPase H+ transporting V0 subunit a4)
    "CLCNKA",     # 氯离子通道Ka (Chloride voltage-gated channel Ka)
    "SLC26A4",    # 溶质载体家族26成员4 (Solute carrier family 26 member 4)
    "SLC12A2",    # 溶质载体家族12成员2 (Solute carrier family 12 member 2)
    "KCNMA1",     # 钾钙激活通道亚家族M α1 (Potassium calcium-activated channel subfamily M alpha 1)
    "SCNN1A"      # 钠离子通道非电压门控1α亚基 (Sodium channel epithelial 1 alpha subunit)
  ),
  
  # 9. 上皮屏障功能 (Epithelial Barrier Function)
  "Epithelial Barrier Function" = c(
    "CLDN3",      # 紧密连接蛋白3 (Claudin-3)
    "CLDN4",      # 紧密连接蛋白4 (Claudin-4)
    "CLDN7",      # 紧密连接蛋白7 (Claudin-7)
    "CLDN18",     # 紧密连接蛋白18 (Claudin-18)
    "CLDN20",     # 紧密连接蛋白20 (Claudin-20)
    "TJP2",       # 紧密连接蛋白2 (ZO-2)
    "TJP3",       # 紧密连接蛋白3 (ZO-3)
    "MARVELD2",   # MARVEL域包含蛋白2 (MARVEL domain containing 2)
    "JAM2",       # 连接粘附分子2 (Junctional adhesion molecule 2)
    "JAM3",       # 连接粘附分子3 (Junctional adhesion molecule 3)
    "F11R",       # F11受体 (F11 receptor)
    "ESAM"        # 内皮选择性粘附分子 (Endothelial cell selective adhesion molecule)
  ),
  
  # 10. 细胞极性与细胞骨架 (Cell Polarity and Cytoskeleton)
  "Cell Polarity and Cytoskeleton" = c(
    "LLGL1",      # Lethal giant larvae同源蛋白1 (Lethal giant larvae homolog 1)
    "LLGL2",      # Lethal giant larvae同源蛋白2 (Lethal giant larvae homolog 2)
    "SCRIB",      # Scribble细胞极性蛋白 (Scribble planar cell polarity protein)
    "DLG1",       # 盘状大肿瘤抑制因子1 (Discs large MAGUK scaffold protein 1)
    "PARD3",      # 分区缺陷3 (Par-3 family cell polarity regulator)
    "PARD6A",     # 分区缺陷6A (Par-6 family cell polarity regulator alpha)
    "PKC1",       # 蛋白激酶C1 (Protein kinase C)
    "CRB3",       # Crumbs细胞极性复合物成分3 (Crumbs cell polarity complex component 3)
    "MPP5",       # 膜蛋白古鸟苷酸激酶相关蛋白5 (Membrane protein palmitoylated 5)
    "PATJ",       # PALS1相关紧密连接蛋白 (PALS1 associated tight junction protein)
    "KRT8",       # 角蛋白8 (Keratin 8)
    "VIL1"        # 微绒毛蛋白1 (Villin-1)
  ),
  
  # 11. 上皮间质转化 (Epithelial-Mesenchymal Transition)
  "Epithelial-Mesenchymal Transition" = c(
    "CDH2",       # N-钙粘蛋白 (N-cadherin)
    "VIM",        # 波形蛋白 (Vimentin)
    "FN1",        # 纤连蛋白1 (Fibronectin 1)
    "SNAI1",      # 蜗牛家族锌指1 (Snail family zinc finger 1)
    "SNAI2",      # 蜗牛家族锌指2 (Snail family zinc finger 2)
    "TWIST1",     # Twist家族bHLH转录因子1 (Twist family bHLH transcription factor 1)
    "TWIST2",     # Twist家族bHLH转录因子2 (Twist family bHLH transcription factor 2)
    "ZEB1",       # 锌指E-box结合同源异型盒1 (Zinc finger E-box binding homeobox 1)
    "ZEB2",       # 锌指E-box结合同源异型盒2 (Zinc finger E-box binding homeobox 2)
    "MMP2",       # 基质金属蛋白酶2 (Matrix metallopeptidase 2)
    "MMP9",       # 基质金属蛋白酶9 (Matrix metallopeptidase 9)
    "ACTA2"       # 肌动蛋白α2 (Actin alpha 2)
  ),
  
  # 12. 细胞增殖与修复 (Cell Proliferation and Repair)
  "Cell Proliferation and Repair" = c(
    "MKI67",      # 标记物Ki-67 (Marker of proliferation Ki-67)
    "PCNA",       # 增殖细胞核抗原 (Proliferating cell nuclear antigen)
    "TOP2A",      # 拓扑异构酶IIα (Topoisomerase II alpha)
    "TK1",        # 胸苷激酶1 (Thymidine kinase 1)
    "CENPW",      # 着丝粒蛋白W (Centromere protein W)
    "CDK1",       # 周期蛋白依赖性激酶1 (Cyclin dependent kinase 1)
    "CCNB1",      # 周期蛋白B1 (Cyclin B1)
    "CCNB2",      # 周期蛋白B2 (Cyclin B2)
    "MCM2",       # 微染色体维持复合物成分2 (Minichromosome maintenance complex component 2)
    "MCM3",       # 微染色体维持复合物成分3 (Minichromosome maintenance complex component 3)
    "TYMS",       # 胸苷酸合酶 (Thymidylate synthetase)
    "RRM2"        # 核糖核苷酸还原酶调节亚基M2 (Ribonucleoside-diphosphate reductase subunit M2)
  ),
  
  # 13. 抗菌防御功能 (Antimicrobial Defense)
  "Antimicrobial Defense" = c(
    "LYZ",        # 溶菌酶 (Lysozyme)
    "LTF",        # 乳铁蛋白 (Lactoferrin)
    "SLPI",       # 分泌型白细胞蛋白酶抑制剂 (Secretory leukocyte peptidase inhibitor)
    "DMBT1",      # 缺失恶性脑肿瘤1 (Deleted in malignant brain tumors 1)
    "RNASE1",     # 核糖核酸酶A家族成员1 (Ribonuclease A family member 1)
    "DEFB1",      # 防御素β1 (Defensin beta 1)
    "DEFB4A",     # 防御素β4A (Defensin beta 4A)
    "CAMP",       # 抗菌肽 (Cathelicidin antimicrobial peptide)
    "S100A7",     # S100钙结合蛋白A7 (S100 calcium-binding protein A7)
    "S100A8",     # S100钙结合蛋白A8 (S100 calcium-binding protein A8)
    "S100A9",     # S100钙结合蛋白A9 (S100 calcium-binding protein A9)
    "LCN2"        # 脂质载体蛋白2 (Lipocalin 2)
  ),
  
  # 14. 代谢相关功能 (Metabolic Functions)
  "Metabolic Functions" = c(
    "SLC34A2",    # 溶质载体家族34成员2 (Solute carrier family 34 member 2)
    "SLC6A14",    # 溶质载体家族6成员14 (Solute carrier family 6 member 14)
    "ABCC3",      # ATP结合盒转运蛋白C3 (ATP binding cassette subfamily C member 3)
    "CYP2F1",     # 细胞色素P450家族2亚家族F成员1 (Cytochrome P450 family 2 subfamily F member 1)
    "ALDH3A1",    # 醛脱氢酶3家族成员A1 (Aldehyde dehydrogenase 3 family member A1)
    "AKR1C1",     # 醛酮还原酶家族1成员C1 (Aldo-keto reductase family 1 member C1)
    "FASN",       # 脂肪酸合酶 (Fatty acid synthase)
    "ACACA",      # 乙酰辅酶A羧化酶α (Acetyl-CoA carboxylase alpha)
    "G6PD",       # 葡萄糖-6-磷酸脱氢酶 (Glucose-6-phosphate dehydrogenase)
    "PKM",        # 丙酮酸激酶M1/2 (Pyruvate kinase M1/2)
    "LDHA",       # 乳酸脱氢酶A (Lactate dehydrogenase A)
    "HK2"         # 己糖激酶2 (Hexokinase 2)
  )
)

# 使用示例 (Usage Example)
# =======================

# 1. 直接使用完整的基因功能列表
results_epithelial <- comprehensive_gene_visualization(
  gene_list = epithelial_gene_functions,
  seurat_obj = Epithelial_object,
  group_by = "Annotation",
  output_dir = "Epithelial_Function_Analysis",
  prefix = "Epithelial_Study",
  max_plot_width = 35,      # 适度限制宽度
  max_plot_height = 25,     # 适度限制高度
  pt_size = 0.3,           # 减小点大小
  verbose = TRUE
)

results <- generate_multi_gene_visualization(
  seurat_obj = Epithelial_object,
  gene_list = epithelial_gene_functions,
  group_by = 'Annotation',
  max_genes_per_plot = 9,  # 每页最多9个基因，便于横向显示
  color_scheme = 'spectral'  # 类似目标图的颜色方案
)

# # 2. 选择特定功能类别进行分析
# selected_functions <- epithelial_gene_functions[c(
#   "通用上皮标记", 
#   "基底细胞标记", 
#   "纤毛上皮标记", 
#   "分泌细胞标记"
# )]
# 
# results_selected <- comprehensive_gene_visualization(
#   gene_list = selected_functions,
#   seurat_obj = your_epithelial_seurat_object,
#   group_by = "annotation_level_2",
#   output_dir = "Selected_Epithelial_Analysis"
# )
# 
# # 3. 基于组织特异性的分析示例
# lung_epithelial_functions <- epithelial_gene_functions[c(
#   "通用上皮标记",
#   "纤毛上皮标记", 
#   "肺泡上皮标记",
#   "分泌细胞标记",
#   "离子转运细胞标记"
# )]

# 示例使用说明
cat("
使用示例 (Usage Example):
=======================

# 1. 准备基因功能列表
gene_functions <- list(
  'T细胞标记' = c('CD3E', 'CD3D', 'CD5', 'TRAC'),
  '髓系细胞标记' = c('CD14', 'CD68', 'ITGAM', 'ITGAX'),
  '炎症相关' = c('IFNG', 'IL10', 'TNF', 'IL6')
)

# 2. 调用主函数
results <- comprehensive_gene_visualization(
  gene_list = gene_functions,
  seurat_obj = your_seurat_object,
  pseudobulk_data = your_pseudobulk_matrix,  # 可选
  group_by = 'cell_type',
  output_dir = 'My_Gene_Analysis',
  prefix = 'Study1'
)

# 3. 查看结果
names(results$plots)  # 查看生成的图形
results$files  # 查看保存的文件路径

# 4. 快速分析
quick_plots <- quick_gene_viz(your_seurat_object, c('CD3E', 'CD14', 'IFNG'))

注意: 所有PDF图表标题均为英文，便于国际期刊发表使用
")

# ============================================================================
# SCENIC和inferCNV完整分析模块 (Complete SCENIC and inferCNV Analysis Modules)
# 基于单细胞转录组数据的转录调控网络和拷贝数变异分析
# ============================================================================

# 必需的包 (Required packages)
scenic_infercnv_packages <- c(
  "SCENIC", "SCopeLoomR", "AUCell", "RcisTarget", "doParallel",
  "infercnv", "Matrix", "ComplexHeatmap", "circlize", "HiddenMarkov",
  "foreach", "doSNOW", "parallel", "future", "BiocParallel", "arrow"
)

# 12. SCENIC转录调控网络分析模块 =============================

#' 执行完整SCENIC转录调控网络分析
#' Perform complete SCENIC transcriptional regulatory network analysis
#'
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名
#' @param output_dir 输出目录
#' @param organism 物种 ("hgnc"为人类, "mgi"为小鼠)
#' @param n_cores 并行计算核心数
#' @param db_dir 数据库目录路径
#' @param min_genes_per_regulon 每个regulon的最小基因数
#' @param min_regulon_gene_occurrence 基因在regulons中的最小出现次数
#' @return 完整SCENIC分析结果
run_scenic_analysis <- function(seurat_obj,
                                cell_type_col = "Annotation",
                                output_dir = "SCENIC_analysis",
                                organism = "hgnc",
                                n_cores = 8,
                                db_dir = "SCENIC_databases",
                                min_genes_per_regulon = 20,
                                min_regulon_gene_occurrence = 5) {
  
  message("开始SCENIC转录调控网络分析 (Starting SCENIC transcriptional regulatory network analysis)...")
  
  # 检查和安装必要的包
  check_and_install_scenic_packages()
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(db_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 加载必要的包
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SCopeLoomR)
  
  # 准备SCENIC输入数据
  scenic_input <- prepare_scenic_input_data(seurat_obj, cell_type_col, output_dir)
  
  # 下载并设置数据库
  setup_scenic_databases(organism, db_dir)
  
  # 初始化SCENIC选项
  scenicOptions <- initialize_scenic_options(organism, db_dir, output_dir)
  
  # 运行完整SCENIC流程
  scenic_results <- run_complete_scenic_pipeline(
    scenic_input, scenicOptions, n_cores, 
    min_genes_per_regulon, min_regulon_gene_occurrence
  )
  
  # 后处理和可视化
  post_process_scenic_results(scenic_results, seurat_obj, cell_type_col, output_dir)
  
  # 整合结果到Seurat对象
  seurat_obj_updated <- integrate_scenic_to_seurat(seurat_obj, scenic_results)
  
  message("SCENIC分析完成 (SCENIC analysis completed)")
  
  return(list(
    scenic_results = scenic_results,
    seurat_obj = seurat_obj_updated,
    scenicOptions = scenicOptions,
    output_dir = output_dir
  ))
}

#' 检查和安装SCENIC相关包
#' Check and install SCENIC packages
check_and_install_scenic_packages <- function() {
  
  required_packages <- c("SCENIC", "AUCell", "RcisTarget", "SCopeLoomR")
  
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("请先安装", pkg, "包。\n",
                 "安装方法：\n",
                 "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
                 "BiocManager::install(c('SCENIC', 'AUCell', 'RcisTarget'))\n",
                 "devtools::install_github('aertslab/SCopeLoomR')"))
    }
  }
  
  message("SCENIC相关包检查完成 (SCENIC packages check completed)")
}

#' 准备SCENIC输入数据
#' Prepare SCENIC input data
prepare_scenic_input_data <- function(seurat_obj, cell_type_col, output_dir) {
  
  message("准备SCENIC输入数据 (Preparing SCENIC input data)...")
  
  # 提取表达矩阵 (使用counts)
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # 转换为常规矩阵
  if(inherits(expr_matrix, "dgCMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  
  # 过滤低表达基因 (在至少1%的细胞中表达)
  min_cells <- round(ncol(expr_matrix) * 0.01)
  genes_keep <- rowSums(expr_matrix > 0) >= min_cells
  expr_matrix_filtered <- expr_matrix[genes_keep, ]
  
  message(paste("保留基因数:", nrow(expr_matrix_filtered), "/", nrow(expr_matrix)))
  
  # 过滤低质量细胞
  min_genes_per_cell <- 200
  cells_keep <- colSums(expr_matrix_filtered > 0) >= min_genes_per_cell
  expr_matrix_filtered <- expr_matrix_filtered[, cells_keep]
  
  message(paste("保留细胞数:", ncol(expr_matrix_filtered), "/", ncol(expr_matrix)))
  
  # 准备细胞信息
  cell_info <- data.frame(
    Cell = colnames(expr_matrix_filtered),
    CellType = seurat_obj@meta.data[colnames(expr_matrix_filtered), cell_type_col],
    stringsAsFactors = FALSE
  )
  
  # 移除NA细胞类型
  valid_cells <- !is.na(cell_info$CellType)
  cell_info <- cell_info[valid_cells, ]
  expr_matrix_filtered <- expr_matrix_filtered[, cell_info$Cell]
  
  # 保存预处理后的数据
  saveRDS(expr_matrix_filtered, file.path(output_dir, "expr_matrix_filtered.rds"))
  write.csv(cell_info, file.path(output_dir, "cell_info.csv"), row.names = FALSE)
  
  return(list(
    expr_matrix = expr_matrix_filtered,
    cell_info = cell_info
  ))
}

#' 设置SCENIC数据库
#' Setup SCENIC databases
setup_scenic_databases <- function(organism, db_dir) {
  
  message("设置SCENIC数据库 (Setting up SCENIC databases)...")
  
  if(organism == "hgnc") {
    # 人类数据库URLs
    db_urls <- list(
      "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather" = 
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
      "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather" = 
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
    )
  } else if(organism == "mgi") {
    # 小鼠数据库URLs
    db_urls <- list(
      "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather" = 
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
      "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather" = 
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
    )
  } else {
    stop("支持的物种: 'hgnc' (人类) 或 'mgi' (小鼠)")
  }
  
  # 下载数据库文件
  for(db_name in names(db_urls)) {
    db_path <- file.path(db_dir, db_name)
    if(!file.exists(db_path)) {
      message(paste("下载数据库:", db_name))
      tryCatch({
        download.file(db_urls[[db_name]], db_path, mode = "wb")
      }, error = function(e) {
        warning(paste("数据库下载失败:", db_name, "-", e$message))
      })
    } else {
      message(paste("数据库已存在:", db_name))
    }
  }
  
  message("SCENIC数据库设置完成 (SCENIC databases setup completed)")
}

#' 初始化SCENIC选项
#' Initialize SCENIC options
initialize_scenic_options <- function(organism, db_dir, output_dir) {
  
  message("初始化SCENIC选项 (Initializing SCENIC options)...")
  
  # 获取数据库文件
  db_files <- list.files(db_dir, pattern = "\\.feather$", full.names = TRUE)
  
  if(length(db_files) < 2) {
    stop("数据库文件不足，请检查数据库下载是否成功")
  }
  
  # 设置数据库路径
  if(organism == "hgnc") {
    db_10kb <- db_files[grepl("10kb", db_files)][1]
    db_500bp <- db_files[grepl("500bp", db_files)][1]
  } else {
    db_10kb <- db_files[grepl("10kb", db_files)][1]
    db_500bp <- db_files[grepl("500bp", db_files)][1]
  }
  
  # 初始化SCENIC选项
  scenicOptions <- initializeScenic(
    org = organism,
    dbDir = db_dir,
    dbs = c("10kb" = db_10kb, "500bp" = db_500bp),
    datasetTitle = "scRNA_SCENIC_Analysis"
  )
  
  # 设置输出目录
  scenicOptions@settings$outDir <- output_dir
  scenicOptions@settings$verbose <- TRUE
  scenicOptions@settings$nCores <- getOption("mc.cores", 1)
  
  return(scenicOptions)
}

#' 运行完整SCENIC流程
#' Run complete SCENIC pipeline
run_complete_scenic_pipeline <- function(scenic_input, scenicOptions, n_cores, 
                                         min_genes_per_regulon, min_regulon_gene_occurrence) {
  
  message("运行完整SCENIC流程 (Running complete SCENIC pipeline)...")
  
  expr_matrix <- scenic_input$expr_matrix
  
  # 设置并行计算
  library(BiocParallel)
  register(MulticoreParam(n_cores))
  
  # 步骤1: 基因过滤
  message("步骤1: 基因过滤 (Step 1: Gene filtering)...")
  genesKept <- geneFiltering(expr_matrix, scenicOptions)
  expr_matrix_filtered <- expr_matrix[genesKept, ]
  
  # 步骤2: 计算基因间相关性
  message("步骤2: 计算基因间相关性 (Step 2: Gene correlation)...")
  runCorrelation(expr_matrix_filtered, scenicOptions)
  
  # 步骤3: 运行GENIE3推断基因调控网络
  message("步骤3: GENIE3基因调控网络推断 (Step 3: GENIE3 gene regulatory network inference)...")
  # 可以设置更多参数来优化GENIE3
  runGenie3(expr_matrix_filtered, scenicOptions, nParts = n_cores)
  
  # 步骤4: 获取共表达模块
  message("步骤4: 共表达模块识别 (Step 4: Co-expression modules)...")
  runSCENIC_1_coexpressionModules(scenicOptions)
  
  # 步骤5: 获取regulons (TF及其靶基因)
  message("步骤5: Regulons识别 (Step 5: Regulons identification)...")
  runSCENIC_2_createRegulons(scenicOptions, 
                             minGenes = min_genes_per_regulon)
  
  # 步骤6: 评估regulon活性 (AUCell)
  message("步骤6: Regulon活性评估 (Step 6: Regulon activity scoring)...")
  runSCENIC_3_scoreCells(scenicOptions, expr_matrix_filtered)
  
  # 步骤7: 二值化AUC值
  message("步骤7: AUC二值化 (Step 7: AUC binarization)...")
  runSCENIC_4_aucell_binarize(scenicOptions, 
                              minGenes = min_genes_per_regulon)
  
  # 加载结果
  regulons <- readRDS(file.path(scenicOptions@settings$outDir, "int/3.4_regulons_forAUCell.Rds"))
  aucell_rankings <- readRDS(file.path(scenicOptions@settings$outDir, "int/3.1_aucell_rankings.Rds"))
  aucell_matrix <- readRDS(file.path(scenicOptions@settings$outDir, "int/3.4_AUCell_binaryRegulonActivity.Rds"))
  
  message("SCENIC核心流程完成 (SCENIC core pipeline completed)")
  
  return(list(
    regulons = regulons,
    aucell_rankings = aucell_rankings,
    aucell_matrix = aucell_matrix,
    scenicOptions = scenicOptions
  ))
}

# 13. inferCNV拷贝数变异分析模块 ==============================

#' 执行完整inferCNV拷贝数变异分析
#' Perform complete inferCNV copy number variation analysis
#'
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名
#' @param normal_cell_types 正常细胞类型 (作为参考)
#' @param malignant_cell_types 恶性细胞类型 (可选)
#' @param gene_order_file 基因顺序文件路径 (可选)
#' @param output_dir 输出目录
#' @param cutoff 表达量截断值
#' @param denoise 是否去噪
#' @param HMM 是否使用HMM预测CNV
#' @param num_threads 线程数
#' @return 完整inferCNV分析结果
run_infercnv_analysis <- function(seurat_obj,
                                  cell_type_col = "Annotation",
                                  normal_cell_types = c("T_cells", "B_cells", "NK_cells"),
                                  malignant_cell_types = NULL,
                                  gene_order_file = NULL,
                                  output_dir = "inferCNV_analysis",
                                  cutoff = 0.1,
                                  denoise = TRUE,
                                  HMM = TRUE,
                                  num_threads = 8) {
  
  message("开始inferCNV拷贝数变异分析 (Starting inferCNV copy number variation analysis)...")
  
  # 检查和安装必要的包
  check_and_install_infercnv_packages()
  
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  library(infercnv)
  
  # 准备inferCNV输入数据
  infercnv_input <- prepare_infercnv_input_data(seurat_obj, cell_type_col, 
                                                normal_cell_types, malignant_cell_types, 
                                                output_dir)
  
  # 获取或创建基因顺序文件
  if(is.null(gene_order_file)) {
    gene_order_file <- create_gene_order_file_from_biomart(infercnv_input$expr_matrix, output_dir)
  }
  
  # 运行完整inferCNV分析
  infercnv_results <- run_complete_infercnv_pipeline(
    infercnv_input, gene_order_file, output_dir, 
    cutoff, denoise, HMM, num_threads
  )
  
  # 后处理和分析结果
  cnv_analysis <- post_process_infercnv_results(infercnv_results, seurat_obj, 
                                                cell_type_col, output_dir)
  
  message("inferCNV分析完成 (inferCNV analysis completed)")
  
  return(list(
    infercnv_obj = infercnv_results,
    cnv_analysis = cnv_analysis,
    gene_order_file = gene_order_file,
    output_dir = output_dir
  ))
}

#' 检查和安装inferCNV相关包
#' Check and install inferCNV packages
check_and_install_infercnv_packages <- function() {
  
  if(!requireNamespace("infercnv", quietly = TRUE)) {
    stop(paste("请先安装inferCNV包。\n",
               "安装方法：\n",
               "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
               "BiocManager::install('infercnv')\n",
               "# 或者从GitHub安装最新版本：\n",
               "devtools::install_github('broadinstitute/inferCNV')"))
  }
  
  message("inferCNV包检查完成 (inferCNV package check completed)")
}

#' 准备inferCNV输入数据
#' Prepare inferCNV input data
prepare_infercnv_input_data <- function(seurat_obj, cell_type_col, normal_cell_types, 
                                        malignant_cell_types, output_dir) {
  
  message("准备inferCNV输入数据 (Preparing inferCNV input data)...")
  
  # 提取表达矩阵 (使用counts)
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # 提取细胞类型信息
  cell_types <- seurat_obj@meta.data[[cell_type_col]]
  names(cell_types) <- colnames(seurat_obj)
  
  # 移除NA细胞类型
  valid_cells <- !is.na(cell_types)
  cell_types <- cell_types[valid_cells]
  expr_matrix <- expr_matrix[, names(cell_types)]
  
  # 确保有足够的正常细胞作为参考
  normal_cells <- sum(cell_types %in% normal_cell_types)
  if(normal_cells < 50) {
    stop(paste("正常参考细胞数量不足 (", normal_cells, ")，建议至少50个正常细胞"))
  }
  
  message(paste("正常参考细胞数量:", normal_cells))
  
  # 过滤低表达基因 (在至少5%的细胞中表达，且平均表达量>0.1)
  min_cells_expr <- round(ncol(expr_matrix) * 0.05)
  gene_detection_rate <- rowSums(expr_matrix > 0)
  gene_mean_expr <- Matrix::rowMeans(expr_matrix)
  
  genes_keep <- (gene_detection_rate >= min_cells_expr) & (gene_mean_expr >= 0.1)
  expr_matrix_filtered <- expr_matrix[genes_keep, ]
  
  message(paste("保留基因数:", nrow(expr_matrix_filtered), "/", nrow(expr_matrix)))
  
  # 创建注释文件
  annotations_df <- data.frame(
    Cell = names(cell_types),
    CellType = as.character(cell_types),
    stringsAsFactors = FALSE
  )
  
  # 保存输入文件
  # 注意：inferCNV需要特定格式的输入文件
  expr_file <- file.path(output_dir, "expr_matrix.txt")
  annot_file <- file.path(output_dir, "annotations.txt")
  
  # 保存表达矩阵
  write.table(as.matrix(expr_matrix_filtered), file = expr_file, 
              sep = "\t", quote = FALSE, col.names = NA)
  
  # 保存注释文件 (inferCNV要求格式: cell_name\tcell_type)
  write.table(annotations_df, file = annot_file, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(list(
    expr_matrix = expr_matrix_filtered,
    expr_file = expr_file,
    annot_file = annot_file,
    annotations = annotations_df,
    normal_cell_types = normal_cell_types
  ))
}

#' 从BioMart创建基因顺序文件
#' Create gene order file from BioMart
create_gene_order_file_from_biomart <- function(expr_matrix, output_dir) {
  
  message("从BioMart获取基因位置信息 (Getting gene positions from BioMart)...")
  
  # 检查biomaRt包
  if(!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("请安装biomaRt包: BiocManager::install('biomaRt')")
  }
  
  library(biomaRt)
  
  genes <- rownames(expr_matrix)
  
  # 连接到Ensembl数据库
  tryCatch({
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    
    # 获取基因位置信息
    gene_info <- getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
      filters = "hgnc_symbol",
      values = genes,
      mart = ensembl
    )
    
    # 过滤掉非标准染色体
    standard_chrs <- c(1:22, "X", "Y")
    gene_info <- gene_info[gene_info$chromosome_name %in% standard_chrs, ]
    
    # 添加chr前缀
    gene_info$chromosome_name <- paste0("chr", gene_info$chromosome_name)
    
    # 按染色体和位置排序
    gene_info$chromosome_name <- factor(gene_info$chromosome_name, 
                                        levels = paste0("chr", c(1:22, "X", "Y")))
    gene_info <- gene_info[order(gene_info$chromosome_name, gene_info$start_position), ]
    
    # 保存基因顺序文件
    gene_order_file <- file.path(output_dir, "gene_order_file.txt")
    write.table(gene_info, file = gene_order_file, 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    message(paste("基因顺序文件已创建，包含", nrow(gene_info), "个基因"))
    
  }, error = function(e) {
    warning(paste("BioMart查询失败:", e$message, "使用备用方法"))
    gene_order_file <- create_simplified_gene_order(genes, output_dir)
  })
  
  return(gene_order_file)
}

#' 创建简化基因顺序文件 (备用方法)
#' Create simplified gene order file (fallback method)
create_simplified_gene_order <- function(genes, output_dir) {
  
  message("创建简化基因顺序文件 (Creating simplified gene order file)...")
  
  # 使用预定义的常见基因位置信息
  # 这只是一个示例，实际应用中应该使用真实的基因组注释
  set.seed(123)  # 保证可重复性
  
  gene_order <- data.frame(
    gene = genes,
    chr = sample(paste0("chr", c(1:22, "X", "Y")), length(genes), replace = TRUE),
    start = sample(1:200000000, length(genes)),
    end = sample(1:200000000, length(genes)),
    stringsAsFactors = FALSE
  )
  
  # 确保start < end
  temp <- gene_order$start
  gene_order$start <- pmin(gene_order$start, gene_order$end)
  gene_order$end <- pmax(temp, gene_order$end)
  
  # 按染色体和位置排序
  gene_order$chr <- factor(gene_order$chr, levels = paste0("chr", c(1:22, "X", "Y")))
  gene_order <- gene_order[order(gene_order$chr, gene_order$start), ]
  
  # 保存基因顺序文件
  gene_order_file <- file.path(output_dir, "gene_order_file_simplified.txt")
  write.table(gene_order, file = gene_order_file, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  warning("使用了简化的基因位置信息，建议手动提供准确的基因顺序文件")
  return(gene_order_file)
}

#' 运行完整inferCNV流程
#' Run complete inferCNV pipeline
run_complete_infercnv_pipeline <- function(infercnv_input, gene_order_file, output_dir, 
                                           cutoff, denoise, HMM, num_threads) {
  
  message("运行完整inferCNV分析流程 (Running complete inferCNV analysis pipeline)...")
  
  # 创建inferCNV对象
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = infercnv_input$expr_file,
    annotations_file = infercnv_input$annot_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = infercnv_input$normal_cell_types
  )
  
  # 运行inferCNV分析
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = cutoff,                    # 表达量截断值
    out_dir = output_dir,               # 输出目录
    cluster_by_groups = TRUE,           # 按组聚类
    denoise = denoise,                  # 去噪
    HMM = HMM,                         # 使用HMM预测CNV状态
    num_threads = num_threads,          # 线程数
    analysis_mode = "subclusters",      # 分析模式
    tumor_subcluster_partition_method = "leiden",  # 聚类方法
    tumor_subcluster_pval = 0.05,      # 聚类p值阈值
    k_obs_groups = 1,                  # 观察组数
    leiden_resolution = 0.00001,       # Leiden分辨率
    per_chr_hmm_subclusters = TRUE,    # 每条染色体独立HMM
    write_expr_matrix = TRUE,          # 保存表达矩阵
    write_phylo = TRUE,                # 保存系统发育树
    output_format = "pdf"              # 输出格式
  )
  
  message("inferCNV核心分析完成 (inferCNV core analysis completed)")
  
  return(infercnv_obj)
}

# 14. 后处理和可视化函数 ====================================

#' SCENIC结果后处理
#' Post-process SCENIC results
post_process_scenic_results <- function(scenic_results, seurat_obj, cell_type_col, output_dir) {
  
  message("SCENIC结果后处理 (Post-processing SCENIC results)...")
  
  # 创建regulon活性热图
  create_regulon_activity_heatmap(scenic_results, seurat_obj, cell_type_col, output_dir)
  
  # 分析细胞类型特异性regulons
  celltype_regulons <- analyze_celltype_specific_regulons(scenic_results, seurat_obj, cell_type_col, output_dir)
  
  # 创建regulon网络图
  create_regulon_network_plot(scenic_results, output_dir)
  
  # 保存regulon活性分数
  save_regulon_activity_scores(scenic_results, output_dir)
  
  return(celltype_regulons)
}

#' inferCNV结果后处理
#' Post-process inferCNV results
post_process_infercnv_results <- function(infercnv_obj, seurat_obj, cell_type_col, output_dir) {
  
  message("inferCNV结果后处理 (Post-processing inferCNV results)...")
  
  # 提取CNV预测结果
  cnv_predictions <- extract_cnv_predictions(infercnv_obj, output_dir)
  
  # 分析细胞类型特异性CNV
  celltype_cnv <- analyze_celltype_cnv_patterns(cnv_predictions, seurat_obj, cell_type_col, output_dir)
  
  # 创建CNV摘要图
  create_cnv_summary_plots(cnv_predictions, seurat_obj, cell_type_col, output_dir)
  
  # 识别高置信度CNV事件
  high_confidence_cnv <- identify_high_confidence_cnv(cnv_predictions, output_dir)
  
  return(list(
    cnv_predictions = cnv_predictions,
    celltype_cnv = celltype_cnv,
    high_confidence_cnv = high_confidence_cnv
  ))
}

#' 整合SCENIC到Seurat对象
#' Integrate SCENIC results to Seurat object
integrate_scenic_to_seurat <- function(seurat_obj, scenic_results) {
  
  message("整合SCENIC结果到Seurat对象 (Integrating SCENIC results to Seurat object)...")
  
  if("aucell_matrix" %in% names(scenic_results)) {
    # 提取regulon活性分数
    regulon_activity <- scenic_results$aucell_matrix
    
    # 匹配细胞名称
    common_cells <- intersect(colnames(seurat_obj), colnames(regulon_activity))
    
    if(length(common_cells) > 0) {
      # 创建新的assay存储regulon活性
      regulon_assay <- CreateAssayObject(data = regulon_activity[, common_cells])
      seurat_obj[["SCENIC"]] <- regulon_assay
      
      message(paste("已添加", nrow(regulon_activity), "个regulons到SCENIC assay"))
    }
  }
  
  return(seurat_obj)
}

# 15. 综合分析函数 ==========================================

#' 运行SCENIC和inferCNV综合分析
#' Run combined SCENIC and inferCNV analysis
run_scenic_infercnv_comprehensive_analysis <- function(seurat_obj,
                                                       cell_type_col = "Annotation",
                                                       normal_cell_types = c("T_cells", "B_cells", "NK_cells"),
                                                       organism = "hgnc",
                                                       run_scenic = TRUE,
                                                       run_infercnv = TRUE,
                                                       output_dir = "SCENIC_inferCNV_comprehensive",
                                                       n_cores = 8) {
  
  message("=== 开始SCENIC和inferCNV综合分析 (Starting comprehensive SCENIC and inferCNV analysis) ===")
  
  # 创建主输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  results <- list()
  
  # 数据质量检查
  perform_data_quality_check(seurat_obj, cell_type_col, normal_cell_types)
  
  # 1. SCENIC分析
  if(run_scenic) {
    message("执行SCENIC转录调控网络分析 (Running SCENIC transcriptional regulatory network analysis)...")
    scenic_dir <- file.path(output_dir, "SCENIC")
    
    results$scenic <- run_scenic_analysis(
      seurat_obj = seurat_obj,
      cell_type_col = cell_type_col,
      output_dir = scenic_dir,
      organism = organism,
      n_cores = n_cores
    )
  }
  
  # 2. inferCNV分析
  if(run_infercnv) {
    message("执行inferCNV拷贝数变异分析 (Running inferCNV copy number variation analysis)...")
    infercnv_dir <- file.path(output_dir, "inferCNV")
    
    results$infercnv <- run_infercnv_analysis(
      seurat_obj = seurat_obj,
      cell_type_col = cell_type_col,
      normal_cell_types = normal_cell_types,
      output_dir = infercnv_dir,
      num_threads = n_cores
    )
  }
  
  # 3. 整合分析
  if(run_scenic && run_infercnv) {
    message("整合SCENIC和inferCNV结果 (Integrating SCENIC and inferCNV results)...")
    
    integrated_results <- integrate_scenic_infercnv_comprehensive(
      results$scenic, 
      results$infercnv, 
      seurat_obj, 
      cell_type_col,
      output_dir
    )
    results$integrated <- integrated_results
  }
  
  # 4. 生成综合报告
  generate_comprehensive_analysis_report(results, seurat_obj, cell_type_col, output_dir)
  
  message("=== SCENIC和inferCNV综合分析完成 (Comprehensive analysis completed) ===")
  
  return(results)
}

#' 数据质量检查
#' Perform data quality check
perform_data_quality_check <- function(seurat_obj, cell_type_col, normal_cell_types) {
  
  message("执行数据质量检查 (Performing data quality check)...")
  
  # 检查基本信息
  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj)
  
  message(paste("细胞数量:", n_cells))
  message(paste("基因数量:", n_genes))
  
  # 检查细胞类型信息
  cell_types <- seurat_obj@meta.data[[cell_type_col]]
  cell_type_counts <- table(cell_types, useNA = "always")
  
  message("细胞类型分布:")
  print(cell_type_counts)
  
  # 检查正常细胞数量
  normal_cells <- sum(cell_types %in% normal_cell_types, na.rm = TRUE)
  message(paste("正常参考细胞数量:", normal_cells))
  
  # 质量检查警告
  if(n_cells < 1000) {
    warning("细胞数量较少 (<1000)，可能影响分析质量")
  }
  
  if(normal_cells < 100) {
    warning("正常参考细胞数量不足 (<100)，inferCNV分析可能不可靠")
  }
  
  if(sum(is.na(cell_types)) > n_cells * 0.1) {
    warning("超过10%的细胞缺少细胞类型注释")
  }
  
  message("数据质量检查完成 (Data quality check completed)")
}

#' 生成综合分析报告
#' Generate comprehensive analysis report
generate_comprehensive_analysis_report <- function(results, seurat_obj, cell_type_col, output_dir) {
  
  message("生成综合分析报告 (Generating comprehensive analysis report)...")
  
  report_file <- file.path(output_dir, "Comprehensive_Analysis_Report.html")
  
  # 创建HTML报告
  report_content <- c(
    "<!DOCTYPE html>",
    "<html><head><title>SCENIC and inferCNV Comprehensive Analysis Report</title></head>",
    "<body>",
    "<h1>SCENIC and inferCNV Comprehensive Analysis Report</h1>",
    paste("<p>Analysis Date:", Sys.time(), "</p>"),
    "<h2>Dataset Summary</h2>",
    paste("<p>Total Cells:", ncol(seurat_obj), "</p>"),
    paste("<p>Total Genes:", nrow(seurat_obj), "</p>"),
    "<h2>Analysis Results</h2>"
  )
  
  if("scenic" %in% names(results)) {
    report_content <- c(report_content,
                        "<h3>SCENIC Analysis</h3>",
                        "<p>✅ Transcriptional regulatory network analysis completed</p>")
  }
  
  if("infercnv" %in% names(results)) {
    report_content <- c(report_content,
                        "<h3>inferCNV Analysis</h3>",
                        "<p>✅ Copy number variation analysis completed</p>")
  }
  
  if("integrated" %in% names(results)) {
    report_content <- c(report_content,
                        "<h3>Integrated Analysis</h3>",
                        "<p>✅ SCENIC and inferCNV results integration completed</p>")
  }
  
  report_content <- c(report_content, "</body></html>")
  
  # 保存报告
  writeLines(report_content, report_file)
  
  message(paste("综合分析报告已保存:", report_file))
}

# ============================================================================
# 使用示例 (Usage Examples)
# ============================================================================

# 示例1: 运行完整SCENIC分析
scenic_results <- run_scenic_analysis(
  seurat_obj = Epithelial_object,
  cell_type_col = "Annotation",
  output_dir = "SCENIC_results",
  organism = "hgnc",
  n_cores = 8
)

# 示例2: 运行完整inferCNV分析 (癌症研究)
# infercnv_results <- run_infercnv_analysis(
#   seurat_obj = your_seurat_obj,
#   cell_type_col = "Annotation",
#   normal_cell_types = c("T_cells", "B_cells", "NK_cells", "Macrophages"),
#   output_dir = "inferCNV_results",
#   cutoff = 0.1,
#   HMM = TRUE,
#   num_threads = 8
# )

# 示例3: 运行综合分析
# comprehensive_results <- run_scenic_infercnv_comprehensive_analysis(
#   seurat_obj = your_seurat_obj,
#   cell_type_col = "Annotation",
#   normal_cell_types = c("T_cells", "B_cells", "NK_cells"),
#   organism = "hgnc",
#   run_scenic = TRUE,
#   run_infercnv = TRUE,
#   output_dir = "SCENIC_inferCNV_comprehensive",
#   n_cores = 8
# )

message("SCENIC和inferCNV完整分析模块加载完成 (Complete SCENIC and inferCNV analysis modules loaded)")

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

# 提取UMAP降维信息并绘图，当然也可使用tsne的降维信息，这里就不展示了
Epithelial_object <- RunUMAP(Epithelial_object, dims = 1:40, reduction = "pca", n.components = 3)
dim(Epithelial_object@reductions$umap)
tmpumap3<-Embeddings(object = Epithelial_object[["umap"]])

cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", 
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")

cb_palette.use <- cb_palette[1:length(unique(Epithelial_object$Annotation))]
col_match <- data.frame(cluster=unique(Epithelial_object$Annotation),col=cb_palette.use)
col_draw<- col_match[match(Epithelial_object$Annotation,col_match[,1]),2]
library(plotly)
tmpumap3 <- as.data.frame(tmpumap3)
fig <- plot_ly(tmpumap3, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color =Epithelial_object$Annotation, colors = cb_palette.use,size=2)
fig

library(rgl)
plot3d(
  tmpumap3,
  col = col_draw,
  type = 'p', radius = .001,axes=T,box=F)

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

# 修复后的手动ROGUE计算函数
manual_rogue_calculation_fixed <- function(expr_matrix, labels, samples, platform = "UMI") {
  
  # 确保所有输入都是字符向量
  labels <- as.character(labels)
  samples <- as.character(samples)
  
  # 获取唯一值
  unique_cell_types <- unique(labels)
  unique_samples <- unique(samples)
  
  # 移除NA值
  unique_cell_types <- unique_cell_types[!is.na(unique_cell_types)]
  unique_samples <- unique_samples[!is.na(unique_samples)]
  
  cat("开始修复后的ROGUE计算...\n")
  cat("有效细胞类型数:", length(unique_cell_types), "\n")
  cat("有效样本数:", length(unique_samples), "\n")
  
  # 创建结果矩阵
  results_matrix <- matrix(NA, 
                           nrow = length(unique_samples), 
                           ncol = length(unique_cell_types),
                           dimnames = list(unique_samples, unique_cell_types))
  
  for (i in seq_along(unique_cell_types)) {
    cell_type <- unique_cell_types[i]
    cat("处理细胞类型:", cell_type, "\n")
    
    for (j in seq_along(unique_samples)) {
      sample <- unique_samples[j]
      
      # 获取特定细胞类型和样本的细胞索引
      cell_idx <- which(labels == cell_type & samples == sample)
      
      if (length(cell_idx) < 10) {  # 最少需要10个细胞
        cat("  ", sample, ": 细胞数不足 (", length(cell_idx), ")\n")
        results_matrix[j, i] <- NA
        next
      }
      
      # 提取子集数据
      subset_expr <- expr_matrix[, cell_idx, drop = FALSE]
      
      tryCatch({
        # 检查数据质量
        if (ncol(subset_expr) < 10 || nrow(subset_expr) < 100) {
          results_matrix[j, i] <- NA
          next
        }
        
        # 过滤低质量数据
        gene_detection <- rowSums(subset_expr > 0)
        cell_detection <- colSums(subset_expr > 0)
        
        # 保留在至少3个细胞中表达的基因
        keep_genes <- gene_detection >= 3
        # 保留至少表达50个基因的细胞
        keep_cells <- cell_detection >= 50
        
        if (sum(keep_genes) < 100 || sum(keep_cells) < 5) {
          results_matrix[j, i] <- NA
          next
        }
        
        subset_expr_filtered <- subset_expr[keep_genes, keep_cells]
        
        # 计算熵
        ent_res <- SE_fun(subset_expr_filtered, span = 0.6, r = 1, mt.method = "fdr")
        
        # 计算ROGUE值
        rogue_val <- CalculateRogue(ent_res, platform = platform, cutoff = 0.05)
        
        results_matrix[j, i] <- rogue_val
        
        cat("  ", sample, ": ROGUE =", round(rogue_val, 3), 
            " (", ncol(subset_expr_filtered), "细胞,", nrow(subset_expr_filtered), "基因)\n")
        
      }, error = function(e) {
        cat("  ", sample, ": 计算失败 -", e$message, "\n")
        results_matrix[j, i] <- NA
      })
    }
  }
  
  # 转换为数据框
  results_df <- as.data.frame(results_matrix)
  return(results_df)
}

# 现在数据格式正确，运行ROGUE分析
cat("=== 运行ROGUE分析 ===\n")
cat("数据准备完成 - 11种细胞类型，94个样本，158,186个细胞\n")

# 先尝试标准ROGUE函数
tryCatch({
  cat("尝试标准ROGUE函数...\n")
  
  rogue_result_standard <- rogue(
    expr = expr_matrix,
    labels = labels,
    samples = samples,
    platform = "UMI",
    min.cell.n = 30,        # 每个样本中每种细胞类型至少30个细胞
    remove.outlier.n = 2,   # 移除2个异常值细胞
    filter = TRUE,          # 启用过滤
    min.cells = 5,          # 基因至少在5个细胞中表达
    min.genes = 50,         # 细胞至少表达50个基因
    span = 0.6,
    r = 1,
    mt.method = "fdr"
  )
  
  cat("标准ROGUE分析成功！\n")
  print("=== ROGUE分析结果矩阵 ===")
  print(dim(rogue_result_standard))
  print(head(rogue_result_standard))
  
  # 保存结果
  write.csv(rogue_result_standard, "epithelial_rogue_results.csv")
  
  # 转换为长格式便于分析
  rogue_long <- rogue_result_standard %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather(key = "cell_type", value = "rogue", -sample) %>%
    filter(!is.na(rogue))
  
  write.csv(rogue_long, "epithelial_rogue_long_format.csv", row.names = FALSE)
  
  success_standard <- TRUE
  
}, error = function(e) {
  cat("标准ROGUE函数失败:", e$message, "\n")
  success_standard <- FALSE
})

# 如果标准方法失败，使用手动方法
if (!exists("success_standard") || !success_standard) {
  cat("使用手动ROGUE计算方法...\n")
  
  rogue_result_manual <- manual_rogue_calculation_fixed(expr_matrix, labels, samples)
  rogue_long <- rogue_result_manual %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather(key = "cell_type", value = "rogue", -sample) %>%
    filter(!is.na(rogue))
} else {
  rogue_result_manual <- rogue_result_standard
}

# 生成详细的统计报告
cat("=== 生成ROGUE分析报告 ===\n")

# 计算每种细胞类型的统计信息
cell_type_stats <- data.frame(
  cell_type = unique_cell_types,
  total_cells = as.numeric(table(labels)[unique_cell_types]),
  stringsAsFactors = FALSE
)

# 添加样本分布信息
sample_distribution <- table(labels, samples)
cell_type_stats$samples_present <- rowSums(sample_distribution > 0)[unique_cell_types]
cell_type_stats$mean_cells_per_sample <- rowMeans(sample_distribution)[unique_cell_types]

print("=== 细胞类型基本统计 ===")
print(cell_type_stats)

# ROGUE结果统计
if (exists("rogue_long") && nrow(rogue_long) > 0) {
  rogue_summary <- rogue_long %>%
    group_by(cell_type) %>%
    summarise(
      n_samples_with_rogue = n(),
      mean_rogue = mean(rogue, na.rm = TRUE),
      median_rogue = median(rogue, na.rm = TRUE),
      sd_rogue = sd(rogue, na.rm = TRUE),
      min_rogue = min(rogue, na.rm = TRUE),
      max_rogue = max(rogue, na.rm = TRUE),
      high_purity_samples = sum(rogue > 0.8, na.rm = TRUE),
      medium_purity_samples = sum(rogue >= 0.5 & rogue <= 0.8, na.rm = TRUE),
      low_purity_samples = sum(rogue < 0.5, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(cell_type_stats, by = "cell_type") %>%
    arrange(desc(mean_rogue))
  
  print("=== ROGUE纯度统计报告 ===")
  print(rogue_summary)
  write.csv(rogue_summary, "epithelial_rogue_summary.csv", row.names = FALSE)
  
  # 创建可视化
  library(ggplot2)
  library(scales)  # 用于alpha()函数
  
  # 设置字体（解决中文显示问题）
  tryCatch({
    theme_set(theme_bw(base_family = "SimHei"))
  }, error = function(e) {
    theme_set(theme_bw(base_family = ""))
    cat("使用默认字体，如果中文显示异常请安装showtext包\n")
  })
  
  # 1. ROGUE分数箱线图 - 改进版
  p1 <- ggplot(rogue_long, aes(x = reorder(cell_type, rogue, median, na.rm = TRUE), 
                               y = rogue, fill = cell_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    theme_bw() +
    labs(
      title = "Epithelial Cell Type ROGUE Purity Scores",
      subtitle = paste("Based on", length(unique_samples), "samples and", 
                       format(sum(cell_type_stats$total_cells), big.mark = ","), "cells"),
      x = "Cell Type",
      y = "ROGUE Purity Score",
      caption = "Higher ROGUE scores indicate better cell type purity"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, color = "gray50")
    ) +
    geom_hline(yintercept = c(0.5, 0.7, 0.8), 
               linetype = "dashed", 
               color = c("red", "orange", "green"),
               alpha = 0.7) +
    annotate("text", x = 1.5, y = 0.82, label = "High Purity (>0.8)", 
             color = "green", size = 3, hjust = 0) +
    annotate("text", x = 1.5, y = 0.72, label = "Medium Purity (0.5-0.8)", 
             color = "orange", size = 3, hjust = 0) +
    annotate("text", x = 1.5, y = 0.52, label = "Low Purity (<0.5)", 
             color = "red", size = 3, hjust = 0) +
    scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1))
  
  # 2. 纯度分数与细胞数量关系 - 改进版
  p2 <- rogue_summary %>%
    ggplot(aes(x = log10(total_cells), y = mean_rogue, 
               size = n_samples_with_rogue, color = cell_type)) +
    geom_point(alpha = 0.8) +
    geom_text(aes(label = cell_type), hjust = 0, vjust = -0.5, 
              size = 3, nudge_x = 0.02, show.legend = FALSE) +
    theme_bw() +
    labs(
      title = "ROGUE Purity vs Cell Count Relationship",
      x = "Total Cell Count (log10)",
      y = "Mean ROGUE Score",
      size = "Samples with ROGUE",
      color = "Cell Type"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14)
    ) +
    scale_x_continuous(breaks = 2:5, labels = c("100", "1K", "10K", "100K")) +
    scale_y_continuous(limits = c(0.85, 1.0)) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # 3. 按纯度等级分布 - 修复版
  purity_counts <- rogue_summary %>%
    select(cell_type, high_purity_samples, medium_purity_samples, low_purity_samples) %>%
    tidyr::gather(key = "purity_level", value = "count", -cell_type) %>%
    mutate(
      purity_level = factor(purity_level, 
                            levels = c("low_purity_samples", "medium_purity_samples", "high_purity_samples"),
                            labels = c("Low Purity (<0.5)", "Medium Purity (0.5-0.8)", "High Purity (>0.8)"))
    )
  
  p3 <- ggplot(purity_counts, aes(x = reorder(cell_type, count), y = count, fill = purity_level)) +
    geom_col(position = "stack", alpha = 0.8) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c("red", "orange", "green")) +
    labs(
      title = "Purity Level Distribution by Cell Type", 
      x = "Cell Type",
      y = "Number of Samples",
      fill = "Purity Level"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    ) +
    geom_text(data = purity_counts %>% 
                group_by(cell_type) %>% 
                summarise(total = sum(count), .groups = "drop"),
              aes(x = cell_type, y = total + 1, label = total, fill = NULL),
              hjust = 0, size = 3)
  
  # 4. 添加一个新的汇总表格图
  p4 <- rogue_summary %>%
    select(cell_type, mean_rogue, total_cells, n_samples_with_rogue) %>%
    mutate(
      purity_category = case_when(
        mean_rogue > 0.8 ~ "High",
        mean_rogue > 0.5 ~ "Medium", 
        TRUE ~ "Low"
      )
    ) %>%
    ggplot(aes(x = reorder(cell_type, mean_rogue), y = mean_rogue)) +
    geom_col(aes(fill = purity_category), alpha = 0.8) +
    geom_text(aes(label = paste0("n=", format(total_cells, big.mark = ","))), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c("High" = "green", "Medium" = "orange", "Low" = "red")) +
    labs(
      title = "ROGUE Scores Summary by Cell Type",
      x = "Cell Type",
      y = "Mean ROGUE Score",
      fill = "Purity Category"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    ) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2))
  
  # 保存图片
  ggsave("epithelial_rogue_boxplot.pdf", p1, width = 14, height = 10)
  ggsave("epithelial_rogue_vs_cellcount.pdf", p2, width = 12, height = 8)
  ggsave("epithelial_purity_distribution.pdf", p3, width = 12, height = 8)
  ggsave("epithelial_rogue_summary_barplot.pdf", p4, width = 12, height = 8)
  
  # 显示图片
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  # 生成一个综合的多面板图
  library(patchwork)
  combined_plot <- (p1 / p4) | p2
  ggsave("epithelial_rogue_combined.pdf", combined_plot, width = 20, height = 12)
  
} else {
  cat("ROGUE分析未成功，无法生成统计报告\n")
}

cat("=== 分析完成 ===\n")
cat("输出文件:\n")
cat("- epithelial_rogue_results.csv: 完整ROGUE结果矩阵\n")
cat("- epithelial_rogue_long_format.csv: 长格式结果\n")
cat("- epithelial_rogue_summary.csv: 统计汇总报告\n")
cat("- epithelial_rogue_boxplot.pdf: 纯度分数箱线图\n")
cat("- epithelial_rogue_vs_cellcount.pdf: 纯度与细胞数关系图\n")
cat("- epithelial_purity_distribution.pdf: 纯度等级分布图\n")

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

# 使用逻辑索引的可靠方法
split_and_save_by_celltype_reliable <- function(seurat_obj, output_dir = "epithelial_celltype_objects") {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cell_types <- unique(as.character(seurat_obj$Annotation))
  cell_types <- cell_types[!is.na(cell_types)]
  
  cat("=== 开始分割保存 ===\n")
  cat("细胞类型:", length(cell_types), "个\n")
  
  for (cell_type in cell_types) {
    safe_name <- gsub("[^A-Za-z0-9_]", "_", cell_type)
    file_path <- file.path(output_dir, paste0(safe_name, ".rds"))
    
    # 使用逻辑索引直接提取细胞
    cells_to_keep <- seurat_obj$Annotation == cell_type
    cell_subset <- seurat_obj[, cells_to_keep]
    
    saveRDS(cell_subset, file_path)
    cat(sprintf("%-15s: %6d cells -> %s\n", 
                cell_type, ncol(cell_subset), basename(file_path)))
  }
  
  cat("=== 完成 ===\n")
}

# 执行
split_and_save_by_celltype_reliable(Epithelial_object)


####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
# 获取所有RDS文件
rds_files <- list.files("epithelial_celltype_objects", pattern = "\\.rds$", full.names = TRUE)
cat("发现", length(rds_files), "个RDS文件\n")

# 创建结果收集文件
results_file <- "celltype_rogue_individual/all_celltype_results.csv"

# 如果结果文件不存在，创建表头
if (!file.exists(results_file)) {
  header_df <- data.frame(
    cell_type = character(),
    total_cells = integer(),
    total_genes = integer(),
    filtered_genes = integer(),
    significant_genes = integer(),
    rogue_score = numeric(),
    mean_entropy = numeric(),
    analysis_date = character()
  )
  write.csv(header_df, results_file, row.names = FALSE)
}

# 逐个分析每个文件
for (i in seq_along(rds_files)) {
  rds_file <- rds_files[i]
  file_name <- basename(rds_file)
  cell_type <- gsub("\\.rds$", "", file_name)
  
  cat("\n=== 第", i, "/", length(rds_files), "个:", cell_type, "===\n")
  
  # 内存状态
  cat("开始前内存状态:\n")
  print(gc())
  
  tryCatch({
    # 加载对象
    cell_obj <- readRDS(rds_file)
    cat("细胞数:", ncol(cell_obj), "基因数:", nrow(cell_obj), "\n")
    
    # 提取表达矩阵 - 保持稀疏格式减少内存使用
    expr_matrix <- GetAssayData(cell_obj, slot = "counts")
    
    # 数据过滤 - 在稀疏格式下进行
    gene_detection <- Matrix::rowSums(expr_matrix > 0)
    cell_detection <- Matrix::colSums(expr_matrix > 0)
    
    keep_genes <- gene_detection >= 5
    keep_cells <- cell_detection >= 100
    
    cat("过滤: 保留", sum(keep_genes), "基因,", sum(keep_cells), "细胞\n")
    
    if (sum(keep_genes) < 500 || sum(keep_cells) < 50) {
      cat("数据质量不足，跳过\n")
      rm(cell_obj, expr_matrix)
      gc()
      next
    }
    
    # 应用过滤
    expr_matrix <- expr_matrix[keep_genes, keep_cells]
    
    # 现在转换为密集矩阵进行ROGUE分析
    cat("转换为密集矩阵进行分析...\n")
    expr_matrix <- as.matrix(expr_matrix)
    
    # 计算熵
    cat("计算基因表达熵...\n")
    entropy_result <- SE_fun(
      expr = expr_matrix,
      span = 0.6,
      r = 1,
      mt.method = "fdr"
    )
    
    # 计算ROGUE分数
    rogue_score <- CalculateRogue(
      entropy_result,
      platform = "UMI",
      cutoff = 0.05
    )
    
    # 统计信息
    sig_genes <- entropy_result[entropy_result$p.adj < 0.05 & !is.na(entropy_result$p.adj), ]
    
    # 保存详细结果
    write.csv(entropy_result, 
              paste0("celltype_rogue_individual/", cell_type, "_entropy.csv"), 
              row.names = FALSE)
    
    # 创建结果行
    result_row <- data.frame(
      cell_type = cell_type,
      total_cells = ncol(cell_obj),
      total_genes = nrow(cell_obj),
      filtered_genes = nrow(entropy_result),
      significant_genes = nrow(sig_genes),
      rogue_score = round(rogue_score, 4),
      mean_entropy = round(mean(entropy_result$entropy, na.rm = TRUE), 4),
      analysis_date = as.character(Sys.Date())
    )
    
    # 追加到结果文件
    write.table(result_row, results_file, sep = ",", append = TRUE, 
                row.names = FALSE, col.names = FALSE)
    
    cat("完成! ROGUE分数:", round(rogue_score, 4), 
        ", 显著基因:", nrow(sig_genes), "\n")
    
  }, error = function(e) {
    cat("分析失败:", e$message, "\n")
  })
  
  # 强制清理内存
  cat("清理内存...\n")
  rm(list = ls(pattern = "cell_obj|expr_matrix|entropy_result|sig_genes"))
  gc()
  
  cat("清理后内存状态:\n")
  print(gc())
  
  # 短暂暂停让系统释放内存
  Sys.sleep(1)
}

cat("\n=== 所有分析完成 ===\n")

# 读取并显示最终结果
final_results <- read.csv(results_file)
print(final_results)

# 按ROGUE分数排序
final_results_sorted <- final_results[order(-final_results$rogue_score), ]
write.csv(final_results_sorted, "celltype_rogue_individual/results_sorted_by_rogue.csv", row.names = FALSE)

cat("结果按ROGUE分数排序:\n")
print(final_results_sorted[, c("cell_type", "total_cells", "rogue_score")])

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

plotData <- as.data.frame(Epithelial_object[["umap"]]@cell.embeddings)

Anno <-as.vector(unique(Epithelial_object$Annotation))
col <- sample(RColorBrewer::brewer.pal(length(Anno), "Paired"))
# 确保plotData包含注释信息
plotData$Annotation <- Epithelial_object$Annotation

# 方案1：使用stat_unchull（推荐参数设置）
ggplot(plotData, aes(x = umap_1, y = umap_2, fill = Annotation, color = Annotation)) +
  stat_unchull(alpha = 0.5, size = 0.5, lty = 1, 
               delta = 1.0,    # 控制包围线扩展距离，可调整为0.5-2.0
               n = 10,         # kNN算法参数，控制点密度计算
               th = 1.0) +     # 简化曲线的距离阈值
  geom_point(size = 0.5, show.legend = FALSE) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(hjust = 0.05, face = "italic")
  ) +
  guides(color = FALSE) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  labs(title = "Epithelial Cell Subpopulations in UMAP Space",
       x = "UMAP_1", y = "UMAP_2")

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
  seurat_obj = Epithelial_object,
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
# check_seurat_structure(Epithelial_object, "tissue", "Annotation")

# 基本调用示例
# results <- create_multiCondition_dimplot(
#   seurat_obj = Epithelial_object,
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
#   seurat_obj = Epithelial_object,
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
#   seurat_obj = Epithelial_object,
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
check_seurat_structure(Epithelial_object, "tissue", "Annotation")

# 基本调用
results <- create_multiCondition_dimplot(
  seurat_obj = Epithelial_object,
  condition_col = "tissue",
  cell_type_col = "Annotation",
  reduction = 'umap',
  pt_size = 0.45,
  alpha = 2
)

# 快速预览
preview_plot <- quick_preview(Epithelial_object, "tissue", "Annotation")
print(preview_plot)

######################################################################################################
######################################################################################################
######################################################################################################
diagnosis <- diagnose_masc_convergence(Epithelial_object)

results <- scPairwiseMASCAnalysis(
  seurat_obj = Epithelial_object,
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
  seurat_obj = Epithelial_object,  # 备用重新计算
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
tissue_sample_counts <- sapply(unique(Epithelial_object$tissue), function(t) {
  length(unique(Epithelial_object$sample[Epithelial_object$tissue == t]))
})
names(tissue_sample_counts) <- unique(Epithelial_object$tissue)
print(tissue_sample_counts)

# 或者使用table + unique组合
unique_combinations <- unique(Epithelial_object@meta.data[, c("tissue", "sample")])
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
  seurat_obj = Epithelial_object,
  cell_type_col = "Annotation",
  output_dir = "publication_heatmaps"
)
#
# # 2. 使用特定基因列表
# marker_genes <- c("CD3E", "CD3D", "CD14", "CD68", "EPCAM", "PECAM1", "COL1A1")
# result <- create_comprehensive_expression_heatmap(
#   seurat_obj = Epithelial_object,
#   cell_type_col = "Annotation", 
#   genes_to_plot = marker_genes,
#   max_cells_per_type = 300,
#   output_prefix = "marker_expression"
# )
#
# # 3. 快速预览
# preview <- quick_preview_heatmap(Epithelial_object)
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
#   seurat_obj = Epithelial_object,
#   cell_type_col = "Annotation",
#   cell_type_colors = custom_colors,
#   expression_color_scheme = c("#2166AC", "#F7F7F7", "#B2182B")
# )

######################################################################################################
######################################################################################################
######################################################################################################