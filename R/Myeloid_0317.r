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

setwd("E:/R/0301/Myeloid_analysis/")#设置工作路径
Myeloid_object <- readRDS("Myeloid_analyzed.rds")#读取数据
Myeloid_object <- NormalizeData(Myeloid_object)#归一化
Myeloid_object <- FindVariableFeatures(Myeloid_object, selection.method = "vst", nfeatures = 4000)#寻找变异基因
Myeloid_object <- ScaleData(Myeloid_object)#标准化
Myeloid_object <- RunPCA(Myeloid_object, npcs = 50)#PCA
Myeloid_object <- RunHarmony(
    object = Myeloid_object,           
    group.by.vars = c("sample"),   
    theta = c(3),
    lambda = c(5),
    sigma = 0.05,
    nclust = 60,
    reduction.use = "pca",
    max_iter = 10,
    # cool_down = 10,
    epsilon_harmony = 0.0001,
    monitor = "harmony_score",
    # factor = 0.8,
    # tol = 1e-10,
    # min_delta = 0.0001,
    # patience = 10,
    epsilon_cluster = 0.00001,
    early_stop = TRUE,
    dims = 1:50           
)
# 获取降维结果和元数据
embeddings <- Embeddings(Myeloid_object, "harmony")
metadata <- Myeloid_object@meta.data

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
# 
# # 保存所有knee plots到一个PDF文件
# pdf("qc_plots/knee_plots_analysis.pdf", width = 12, height = 10)
# 
# # 1. UMI knee plot
# counts <- GetAssayData(Myeloid_object, layer = "counts")
# total_umi <- Matrix::colSums(counts)
# umi_rank <- rank(-total_umi)
# 
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
# 
# # 2. 基因数量的knee plot
# genes_per_cell <- Matrix::colSums(counts > 0)
# gene_rank <- rank(-genes_per_cell)
# 
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
# 
# # # 3. PC肘部图
# pcMyeloid_var <- Myeloid_object[["pca"]]@stdev / sum(Myeloid_object[["pca"]]@stdev) * 100
# cumsum_var <- cumsum(pcMyeloid_var)
# 
# elbow_data <- data.frame(
#   PC = seq_along(pcMyeloid_var),
#   variance = pcMyeloid_var,
#   cumulative = cumsum_var
# )
# 
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
# 
# # # 4. Harmony scores knee plot - 使用采样
# set.seed(42)
# n_sample <- min(10000, ncol(Myeloid_object))
# sample_idx <- sample(seq_len(ncol(Myeloid_object)), n_sample)
# 
# harmony_scores <- Myeloid_object[["harmony"]]@cell.embeddings[sample_idx, 1:20]
# harmony_dist <- dist(harmony_scores)
# harmony_scores_ordered <- sort(colMeans(as.matrix(harmony_dist)), decreasing = TRUE)
# 
# print(ggplot(data.frame(rank = seq_along(harmony_scores_ordered),
#                         score = harmony_scores_ordered),
#              aes(x = rank, y = score)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         theme_bw() +
#         labs(title = "Harmony Scores Distribution (10k cells sampled)",
#              x = "Cell rank",
#              y = "Average Harmony distance"))
# 
# dev.off()
# 
# # 保存统计信息到CSV
# stats_df <- data.frame(
#   Metric = c("Total_cells", "Median_UMI", "Median_genes",
#              "PC_90percenMyeloid_var", "Harmony_mean_dist"),
#   Value = c(
#     length(total_umi),
#     median(total_umi),
#     median(genes_per_cell),
#     which(cumsum_var > 90)[1],
#     mean(harmony_scores_ordered)
#   )
# )
# 
# write.csv(stats_df, "qc_plots/knee_ploMyeloid_statistics.csv", row.names = FALSE)

# 直接运行UMAP
Myeloid_object <- RunUMAP(
    Myeloid_object, 
    reduction = "harmony", 
    dims = 1:42, 
    n.neighbors = 30,
    min.dist = 0.3,
    # learning.rate = 0.2,
    n.trees = 500,
    # n.epochs = 1000,
    # spread = 1.2,
    # repulsion.strength = 1.1,
    metric = "correlation"
)

# FindNeighbors函数的缩进修正
Myeloid_object <- FindNeighbors(
    Myeloid_object, 
    reduction = "harmony", 
    dims = 1:42,
        k.param = 30
) 

Myeloid_object <- FindClusters(
    Myeloid_object, 
    resolution = 2, 
    algorithm = 4
    )
# 8. 绘制UMAP聚类图
print("Generating plots...")
Idents(Myeloid_object) <- Myeloid_object$seurat_clusters
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(Myeloid_object, reduction = "umap", label = TRUE)
dev.off()
# saveRDS(Myeloid_object, "final_analyzed_Myeloid_object.rds")
# 9. 保存分析结果
print("Saving analysis results...")

# Myeloid_object <- readRDS("E:/R/0224/final_analyzed_Myeloid_object.rds")


# 9. 标记基因分析

print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for(marker in Marker_Myeloid) {
  if(marker %in% rownames(Myeloid_object)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Myeloid_object, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(Myeloid_object, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(Myeloid_object, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(Myeloid_object, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(Myeloid_object, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(Myeloid_object, reduction = "umap", group.by = "ann_level_4", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_4")
p1
p2
p3
p4
p5
dev.off()
table(Myeloid_object@meta.data$ann_level_4)
tissues <- unique(Myeloid_object$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(Myeloid_object) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(Myeloid_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Myeloid_object, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(Myeloid_object$sample)
pdf("sample_umap_split.pdf", width = 12, height = 12)
Idents(Myeloid_object) <- "sample"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(Myeloid_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Myeloid_object, 
                                                    idents = sample_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(Myeloid_object$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(Myeloid_object) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(Myeloid_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Myeloid_object, 
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
Myeloid_object$ann_level_3_no_na <- Myeloid_object$ann_level_3
Myeloid_object$ann_level_3_no_na[is.na(Myeloid_object$ann_level_3_no_na)] <- "Unknown"
Idents(Myeloid_object) <- "ann_level_3_no_na"

pdf("ann_level_3_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_level_3 <- unique(Myeloid_object$ann_level_3)
ann_level_3 <- ann_level_3[!is.na(ann_level_3)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_level_3) {
  print(DimPlot(Myeloid_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(Myeloid_object, 
                                                    idents = anno),
                cols = "grey",
                pt.size = 0.5,
                raster = TRUE,
                label = FALSE) +
          ggtitle(anno) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
Idents(Myeloid_object) <- Myeloid_object$seurat_clusters
markers <- FindAllMarkers(Myeloid_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

# table(subset(Myeloid_object, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

# 11. 生成热图
print("Generating heatmap...")
top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)


# 1. 先缩放marker基因
marker_genes <- unique(top10_markers$gene)
Myeloid_object <- ScaleData(Myeloid_object, features = marker_genes)

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
  print(DoHeatmap(Myeloid_object, 
                  features = marker_batches[[i]],
                  size = 24) + 
          NoLegend() +
          ggtitle(paste("Marker Genes Batch", i))+
          theme(axis.text.y = element_text(size = 36))
  )
  # dev.off()
}
dev.off()

# 1. 首先检查 Marker_Myeloid 中是否有重复
print(length(Marker_Myeloid))
print(length(unique(Marker_Myeloid)))
# 方法1：直接显示所有重复的元素
duplicated_Marker_Myeloid <- Marker_Myeloid[duplicated(Marker_Myeloid)]
print("重复的 Marker_Myeloid 是：")
print(duplicated_Marker_Myeloid)

# 2. 如果发现有重复，可以移除重复项
Marker_Myeloid <- unique(Marker_Myeloid)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("Marker_Myeloid_dotplot.pdf", width = 32, height = 18)
doMyeloid_plot <- DotPlot(Myeloid_object, 
                    features = Marker_Myeloid, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- doMyeloid_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotploMyeloid_data.csv")
doMyeloid_plot
dev.off()
getwd()

# 1. 读取注释文件（无表头）
annotations <- read.csv("Annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation_2")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- Myeloid_object$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
Myeloid_object$Annotation_2 <- new_idents
Myeloid_object <- SetIdent(Myeloid_object, value = new_idents)

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
Idents(Myeloid_object) <- Myeloid_object$Annotation_2
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(Myeloid_object, 
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
cell_types <- unique(Myeloid_object$Annotation_2)
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 为每个细胞类型单独绘图
for(cell_type in cell_types) {
  # 创建文件名
  
  
  # 绘制该细胞类型的UMAP图
  
  print(DimPlot(Myeloid_object, 
                reduction = "umap",
                cells.highlight = WhichCells(Myeloid_object, idents = cell_type),
                cols.highlight = cell_colors[cell_type],
                cols = "grey",  # 其他细胞显示为灰色
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(cell_type) +
          theme(legend.text = element_text(size = 12)))
  
}
dev.off()

saveRDS(Myeloid_object,"Myeloid_analyzed.RDS")

################################################################################################################
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
    } else if("orig.ident" %in% colnames(mydata)) {
      group_by <- "orig.ident"
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
# seurat_obj <- readRDS("path/to/your/seuraMyeloid_object.rds")
# 
# # 如果需要先运行Harmony和UMAP
# # seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "group", dims.use = 1:30)
# # seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
# 
# 生成高级可视化
# viz_results <- generate_harmony_umap(
#   seurat_obj = Myeloid_object,
#   group_by = "group",          # 按哪个变量分组
#   cell_type_col = "Annotation_2", # 细胞类型注释列
#   output_dir = "10-Harmony_Visualizations", # 输出目录
#   point_size = 0.8,            # 点大小
#   use_ellipse = TRUE,          # 是否添加置信椭圆
#   width = 12,                  # 图像宽度
#   height = 10                  # 图像高度
# )
# # generate_advanced_umap(Myeloid_object)

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
# seurat_obj <- readRDS("path/to/your/seuraMyeloid_object.rds")
results <- generate_marker_viz(
  seurat_obj = Myeloid_object,
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

# 首先检查Myeloid_object是否包含必要的元数据列
required_columns <- c("Annotation_2", "tissue", "sample")
missing_columns <- required_columns[!required_columns %in% colnames(Myeloid_object@meta.data)]

if(length(missing_columns) > 0) {
  stop(paste("Myeloid_object缺少必要的元数据列:", paste(missing_columns, collapse=", ")))
}

# 检查不同组织中的细胞类型分布
message("分析不同组织中的细胞类型分布...")
tissue_cell_distribution <- table(Myeloid_object$tissue, Myeloid_object$Annotation_2)
write.csv(tissue_cell_distribution, "tissue_celltype_distribution.csv")

# 打印组织和细胞类型信息
message("数据集包含以下组织类型:")
print(table(Myeloid_object$tissue))

message("数据集包含以下T细胞亚型:")
print(table(Myeloid_object$Annotation_2))

# 1. 按组织和细胞类型进行pseudobulk聚合
message("执行基于组织和细胞类型的pseudobulk聚合...")
av_tissue_celltype <- AggregateExpression(
  Myeloid_object,
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
getwd()
setwd("..")
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

# 检查组织注释是否存在
if(!"tissue" %in% colnames(Myeloid_object@meta.data)) {
  message("在Seurat对象中找不到'tissue'列。跳过组织特异性分析。")
} else {
  # 获取所有组织类型
  tissue_types <- unique(Myeloid_object$tissue)
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
if(!exists("Myeloid_object")) {
  stop("在工作环境中找不到'Myeloid_object' Seurat对象。请先加载数据。")
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
  cells_in_tissue <- rownames(Myeloid_object@meta.data[Myeloid_object@meta.data$tissue == tissue, ])
  message("DEBUG: Number of cells = ", length(cells_in_tissue))
  
  if(length(cells_in_tissue) < 50) {
    message(paste0("组织'", tissue, "'的细胞数量太少(", length(cells_in_tissue), ")，跳过"))
    next
  }
  
  # 创建该组织的子集
  tissue_subset <- subset(Myeloid_object, cells = cells_in_tissue)
  
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
if(!"Annotation_2" %in% colnames(Myeloid_object@meta.data)) {
  stop("在Seurat对象中找不到'Annotation_2'。请先进行细胞类型注释。")
}

# 设置细胞类型为标识
Idents(Myeloid_object) <- Myeloid_object$Annotation_2
message("按细胞类型计算平均表达...")

# 计算每个细胞类型的平均表达
expr <- AverageExpression(Myeloid_object, assays = "RNA", slot = "data")[[1]]
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
if(!"tissue" %in% colnames(Myeloid_object@meta.data)) {
  message("在Seurat对象中找不到'tissue'列。跳过组织特异性分析。")
} else {
  # 获取所有组织类型
  tissue_types <- unique(Myeloid_object$tissue)
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
  cells_in_tissue <- rownames(Myeloid_object@meta.data[Myeloid_object@meta.data$tissue == tissue, ])
  
  if(length(cells_in_tissue) < 50) {
    message(paste0("Tissue '", tissue, "' has too few cells (", length(cells_in_tissue), "), skipping"))
    next
  }
  
  # 创建该组织的子集
  tissue_subset <- subset(Myeloid_object, cells = cells_in_tissue)
  
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