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
             'FCER1G', 'C1orf162', 'CLEC7A', 'CD1C', 'CD86', 'CD14', 'XCR1', 'HLA-DRA',#Fibroblast
             'COL1A2', 'DCN', 'MFAP4', 'LUM', 'COL6A3', 'CFD', 'COL1A1', 'PDGFRA', 
             'MXRA8', 'LEPR',
             'MYH11', 'TINAGL1', 'PLN', 'DES', 'ACTA2', 'CNN1', 'TAGLN', #Fibroblast & SMC
             'CLDN5', 'ECSCR', 'CLEC14A', 'VWF', 'PECAM1', 'DARC', 'PTPRB', 'PDE2A', 
             'PLAT', 'GJA5', 'SPARCL1', 'AQP1', 'MMRN1', 'CCL21', #Endothelial
             'MKI67', 'TOP2A', 'TK1', 'CENPW')#Proliferation
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

Marker_Fibroblast_V2 <- c("ADGRB3", "SFRP2", "CXCL2", "DCN", "PTGDS",
                          "COL4A2", "CCL19", "STEAP4", "RGS5", "CLSTN2",
                          "SFRP4", "COCH", "MIR99AHG", "C2orf40", "PTN",
                          "STATH", "LYZ", "SLPI", "BPIFA1", "ZG16B",
                           "TAGLN", "C11orf96", "CACNB2", "RERGL", "MYH11",
                           "KRT19", "SERPINB3", "KRT17", "AQP3", "S100A2",
                           "POSTN", "ROBO1", "CCL26", "AC079298.3", "PAPPA"
)

setwd("E:/R/0301/Fibroblast_analysis/")#设置工作路径
Fibroblast_object <- readRDS("Fibroblast_analyzed.rds")#读取数据
Fibroblast_object <- NormalizeData(Fibroblast_object)#归一化
Fibroblast_object <- FindVariableFeatures(Fibroblast_object, selection.method = "vst", nfeatures = 2000)#寻找变异基因
Fibroblast_object <- ScaleData(Fibroblast_object)#标准化
Fibroblast_object <- RunPCA(Fibroblast_object, npcs = 40)#PCA
Fibroblast_object <- RunHarmony(
    object = Fibroblast_object,           
    group.by.vars = c("sample"),   
    theta = c(3),
    lambda = c(1),
    sigma = 0.1,           
    # nclust = 4,
    reduction.use = "pca",
    max_iter = 10,
    # cool_down = 10,
    epsilon_harmony = 0.0001,
    monitor = "harmony_score",
    # factor = 0.8,
    # tol = 1e-10,
    # min_delta = 0.0001,
    # patience = 10,
    # epsilon_cluster = 0.00001,
    # early_stop = FALSE,
    dims = 1:40           
)
# # 获取降维结果和元数据
# embeddings <- Embeddings(Fibroblast_object, "harmony")
# metadata <- Fibroblast_object@meta.data
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

# # 创建输出目录（如果不存在）
# if(!dir.exists("qc_plots")) {
#   dir.create("qc_plots")
# }
# 
# # 保存所有knee plots到一个PDF文件
# pdf("qc_plots/knee_plots_analysis.pdf", width = 12, height = 10)
# 
# # 1. UMI knee plot
# counts <- GetAssayData(Fibroblast_object, layer = "counts")
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
# pcFibroblast_var <- Fibroblast_object[["pca"]]@stdev / sum(Fibroblast_object[["pca"]]@stdev) * 100
# cumsum_var <- cumsum(pcFibroblast_var)
# 
# elbow_data <- data.frame(
#   PC = seq_along(pcFibroblast_var),
#   variance = pcFibroblast_var,
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
# n_sample <- min(10000, ncol(Fibroblast_object))
# sample_idx <- sample(seq_len(ncol(Fibroblast_object)), n_sample)
# 
# harmony_scores <- Fibroblast_object[["harmony"]]@cell.embeddings[sample_idx, 1:20]
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
#              "PC_90percenFibroblast_var", "Harmony_mean_dist"),
#   Value = c(
#     length(total_umi),
#     median(total_umi),
#     median(genes_per_cell),
#     which(cumsum_var > 90)[1],
#     mean(harmony_scores_ordered)
#   )
# )
# 
# write.csv(stats_df, "qc_plots/knee_ploFibroblast_statistics.csv", row.names = FALSE)

# 直接运行UMAP
Fibroblast_object <- RunUMAP(
    Fibroblast_object, 
    reduction = "harmony", 
    dims = 1:40, 
    n.neighbors = 20,
    min.dist = 0.3,
    # learning.rate = 0.2,
    n.trees = 500,
    # n.epochs = 1000,
    # spread = 1.2,
    # repulsion.strength = 1.1,
    metric = "correlation"
)

# FindNeighbors函数的缩进修正
Fibroblast_object <- FindNeighbors(
    Fibroblast_object, 
    reduction = "harmony", 
    dims = 1:40,
        k.param = 20
) 

Fibroblast_object <- FindClusters(
    Fibroblast_object, 
    resolution = 2, 
    algorithm = 2
    )
# 8. 绘制UMAP聚类图
print("Generating plots...")
Idents(Fibroblast_object) <- Fibroblast_object$seurat_clusters
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(Fibroblast_object, reduction = "umap", label = TRUE)
dev.off()
# saveRDS(Fibroblast_object, "final_analyzed_Fibroblast_object.rds")
# 9. 保存分析结果
print("Saving analysis results...")

# Fibroblast_object <- readRDS("E:/R/0224/final_analyzed_Fibroblast_object.rds")


# 9. 标记基因分析

print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for(marker in Marker_Fibroblast) {
  if(marker %in% rownames(Fibroblast_object)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Fibroblast_object, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

print("Generating feature plots...")
pdf("umap_FeaturePlot_V2.pdf", width = 8, height = 8)
for(marker in Marker_Fibroblast_V2) {
  if(marker %in% rownames(Fibroblast_object)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Fibroblast_object, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(Fibroblast_object, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(Fibroblast_object, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(Fibroblast_object, reduction = "umap", group.by = "ann_level_3", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_3")
p4 <- DimPlot(Fibroblast_object, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(Fibroblast_object, reduction = "umap", group.by = "ann_level_4", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_4")
p1
p2
p3
p4
p5
dev.off()

tissues <- unique(Fibroblast_object$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(Fibroblast_object) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(Fibroblast_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Fibroblast_object, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(Fibroblast_object$sample)
pdf("sample_umap_split.pdf", width = 12, height = 12)
Idents(Fibroblast_object) <- "sample"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(Fibroblast_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Fibroblast_object, 
                                                    idents = sample_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(Fibroblast_object$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(Fibroblast_object) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(Fibroblast_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = Fibroblast_object, 
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
Fibroblast_object$ann_level_3_no_na <- Fibroblast_object$ann_level_3
Fibroblast_object$ann_level_3_no_na[is.na(Fibroblast_object$ann_level_3_no_na)] <- "Unknown"
Idents(Fibroblast_object) <- "ann_level_3_no_na"

pdf("ann_level_3_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_level_3 <- unique(Fibroblast_object$ann_level_3)
ann_level_3 <- ann_level_3[!is.na(ann_level_3)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_level_3) {
  print(DimPlot(Fibroblast_object, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(Fibroblast_object, 
                                                    idents = anno),
                cols = "grey",
                pt.size = 0.5,
                raster = TRUE,
                label = FALSE) +
          ggtitle(anno) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

Fibroblast_object@meta.data$QC[Fibroblast_object$seurat_clusters %in% c('16','17','33','26')] <- '1'

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
pdf("markers_dotplot.pdf", width = 32, height = 18)
doFibroblast_plot <- DotPlot(Fibroblast_object, 
                    features = Marker_Fibroblast, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- doFibroblast_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotploFibroblast_data.csv")
doFibroblast_plot
dev.off()
getwd()

table(Fibroblast_object$seurat_clusters)



# 2. 计算每个cluster中每个基因的表达情况
print("Calculating expression statistics...")

# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
Idents(Fibroblast_object) <- Fibroblast_object$seurat_clusters
markers <- FindAllMarkers(Fibroblast_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

# table(subset(Fibroblast_object, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

# 11. 生成热图
print("Generating heatmap...")
top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
write.csv(top10_markers, "cluster_top10_markers.csv")

# 1. 先缩放marker基因
marker_genes <- unique(top10_markers$gene)
Fibroblast_object <- ScaleData(Fibroblast_object, features = marker_genes)

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
  print(DoHeatmap(Fibroblast_object, 
                  features = marker_batches[[i]],
                  size = 24) + 
          NoLegend() +
          ggtitle(paste("Marker Genes Batch", i))+
          theme(axis.text.y = element_text(size = 36))
  )
  # dev.off()
}
dev.off()


saveRDS(Fibroblast_object,"Fibroblast_analyzed0928.rds")

#######################################################################################################
#######################################################################################################
#######################################################################################################
# 定义要检查的标记基因
fibroblast_markers <- c("BAMBI", "CLDN1", "COL11A1", "CXCL14", "DPT", "RGS5",'PDGFRB','NDUFA4L2','NOTCH3')

# 检查基因是否存在
available_genes <- intersect(fibroblast_markers, rownames(Fibroblast_object))
cat("可用基因:", paste(available_genes, collapse = ", "), "\n")

# 方法1: 单个基因的violin图（类似原图风格）
violin_plots <- list()
for(gene in available_genes) {
  p <- VlnPlot(Fibroblast_object, 
               features = gene,
               group.by = "seurat_clusters",
               pt.size = 0) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = gene, y = "Expression")
  
  violin_plots[[gene]] <- p
}

# 组合图形（垂直排列，类似原图）
combined_plot <- wrap_plots(violin_plots, ncol = 1)

# 保存图形
ggsave("fibroblast_cluster_markers_violin.pdf", 
       combined_plot, 
       width = 8, 
       height = length(available_genes) * 2.5, 
       dpi = 300)

# 方法2: 所有基因在一个图中（堆叠式）
stacked_plot <- VlnPlot(Fibroblast_object, 
                        features = available_genes,
                        group.by = "seurat_clusters",
                        stack = TRUE,
                        flip = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("fibroblast_markers_stacked_violin.pdf", 
       stacked_plot, 
       width = 10, 
       height = 8)

# 方法3: 网格排列
grid_plot <- VlnPlot(Fibroblast_object, 
                     features = available_genes,
                     group.by = "seurat_clusters",
                     ncol = 2,
                     pt.size = 0) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("fibroblast_markers_grid_violin.pdf", 
       grid_plot, 
       width = 12, 
       height = 8)

# 快速查看cluster信息
cat("Cluster数量:", length(unique(Fibroblast_object$seurat_clusters)), "\n")
print(table(Fibroblast_object$seurat_clusters))

cat("\n图形已保存:\n")
cat("1. fibroblast_cluster_markers_violin.pdf - 垂直排列（类似原图）\n")
cat("2. fibroblast_markers_stacked_violin.pdf - 堆叠式\n") 
cat("3. fibroblast_markers_grid_violin.pdf - 网格排列\n")

#######################################################################################################
#######################################################################################################
#######################################################################################################