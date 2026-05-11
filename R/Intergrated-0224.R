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
options(future.globals.maxSize = Inf)
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

Markers <- c('FXYD3', 'EPCAM', 'ELF3', 'IGFBP2', 'SERPINF1', 'TSPAN1',
             'SCGB1A1', 'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8',
             'CD53', 'PTPRC', 'CORO1A', 'ISG20', 'CCL5',
             'MS4A1', 'TNFRSF17', 'CD19', 'CD79A', 'SDC1',
             'CD40LG', 'TNFRSF25', 'CD28', 'CD4', 'CD3E', 'CD8A', 'CD8B', 
             'TRGC2', 'CD2', 'TRBC2',
             'FCER1G', 'C1orf162', 'CLEC7A', 'CD1C', 'CD86', 'CD14', 'XCR1', 'HLA-DRA',
             'COL1A2', 'DCN', 'MFAP4', 'LUM', 'COL6A3', 'CFD', 'COL1A1', 'PDGFRA', 
             'MXRA8', 'NBL1', 'VCAN', 'LEPR',
             'MYH11', 'TINAGL1', 'PLN', 'DES', 'ACTA2', 'CNN1', 'TAGLN',
             'CLDN5', 'ECSCR', 'CLEC14A', 'VWF', 'PECAM1', 'DARC', 'PTPRB', 'PDE2A', 
             'PLAT', 'GJA5', 'SPARCL1', 'AQP1', 'RNASE1', 'MMRN1', 'CCL21', 'TFF3',
             'MKI67', 'TOP2A', 'TK1', 'CENPW')
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
# 设置工作目录
setwd("E:/R/0225/")

# 1. 读取所有RDS文件
batch_files <- list.files("E:\\R\\0216\\rds", pattern = "\\.rds$", full.names = TRUE)
if(length(batch_files) == 0) {
  stop("No RDS files found!")
}


# 1. 读取第一个数据集作为基础
print(sprintf("Reading first dataset: %s", basename(batch_files[1])))
seurat_obj_main <- readRDS(batch_files[1])
seurat_obj_main <- DietSeurat(seurat_obj_main, 
                              counts = TRUE,
                              data = TRUE,
                              scale.data = FALSE,
                              dimreducs = NULL,
                              graphs = NULL)

# 预处理第一个数据集
print("Processing first dataset...")
seurat_obj_main <- NormalizeData(seurat_obj_main, verbose = FALSE)
seurat_obj_main <- FindVariableFeatures(seurat_obj_main, 
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
seurat_obj_main <- ScaleData(seurat_obj_main, verbose = FALSE)
seurat_obj_main <- RunPCA(seurat_obj_main, verbose = FALSE)

# 2. 逐个处理其他数据集
for(i in 2:length(batch_files)) {
  print(sprintf("\nProcessing dataset %d/%d: %s", i, length(batch_files), basename(batch_files[i])))
  
  # 读取当前数据集
  current_data <- readRDS(batch_files[i])
  current_data <- DietSeurat(current_data,
                             counts = TRUE,
                             data = TRUE,
                             scale.data = FALSE,
                             dimreducs = NULL,
                             graphs = NULL)
  
  # 预处理当前数据集
  print("Normalizing and finding variable features...")
  current_data <- NormalizeData(current_data, verbose = FALSE)
  current_data <- FindVariableFeatures(current_data,
                                       selection.method = "vst",
                                       nfeatures = 2000,
                                       verbose = FALSE)
  
  # 找到共同变异基因
  genes.use <- intersect(VariableFeatures(seurat_obj_main),
                         VariableFeatures(current_data))
  print(sprintf("Number of common variable features: %d", length(genes.use)))
  
  # 准备整合
  print("Finding transfer anchors...")
  integration.anchors <- FindTransferAnchors(
    reference = seurat_obj_main,
    query = current_data,
    dims = 1:30,
    reference.reduction = "pca",
    k.anchor = 5,
    k.filter = 100,
    k.score = 20,
    max.features = 200,
    features = genes.use
  )
  
  # 整合数据
  print("Mapping query data...")
  current_data <- MapQuery(
    anchorset = integration.anchors,
    reference = seurat_obj_main,
    query = current_data,
    refdata = list(
      counts = "counts",
      data = "data"
    ),
    reference.reduction = "pca",
    reduction.model = "pca"
  )
  
  # 合并数据集
  print("Merging with main object...")
  seurat_obj_main <- merge(seurat_obj_main, current_data)
  
  # 清理内存
  rm(current_data, integration.anchors, genes.use)
  gc()
  
  # 打印当前内存使用情况
  print(sprintf("Current object size: %d cells, %d features", 
                ncol(seurat_obj_main), 
                nrow(seurat_obj_main)))
}

# 3. 过滤不需要的样本
print("\nFiltering unwanted samples...")
cells_before <- ncol(seurat_obj_main)
seurat_obj_main <- subset(seurat_obj_main, 
                          subset = tissue_sampling_method != 'brush' & 
                            study != 'Jain_Misharin_2021')
cells_after <- ncol(seurat_obj_main)


# 8. 清理中间对象节省内存
rm(merged_list, anchors)
gc()

# 获取所有unique的sample并进行doublet检测
samples <- unique(seurat_obj_main$sample)
samples <- samples[!is.na(samples)]

# 为每个sample创建一个列表来存储结果
sample_results <- list()

# 循环处理每个sample进行doublet检测
for(sample_name in samples) {
  cat(sprintf("\nProcessing sample: %s\n", sample_name))
  
  # 提取当前sample的数据
  current_obj <- subset(seurat_obj_main, subset = sample == sample_name)
  
  # 获取细胞数量
  n_cells <- ncol(current_obj)
  if(n_cells < 50) {
    cat("Too few cells, skipping...\n")
    next
  }
  
  # 计算预期双细胞率
  doublet_rate <- min(0.075, n_cells / 1000 * 0.008)
  nExp_poi <- round(doublet_rate * n_cells)
  cat(sprintf("Expected doublets: %d (rate: %.3f)\n", nExp_poi, doublet_rate))
  
  # 转换回RNA assay进行doublet检测
  DefaultAssay(current_obj) <- "RNA"
  
  # 数据准备
  current_obj <- NormalizeData(current_obj, verbose = FALSE)
  current_obj <- FindVariableFeatures(current_obj, verbose = FALSE)
  current_obj <- ScaleData(current_obj, verbose = FALSE)
  current_obj <- RunPCA(current_obj, npcs = 30, verbose = FALSE)
  
  # 参数优化和doublet检测
  tryCatch({
    sweep.res.list <- paramSweep(current_obj, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    # 确保pK是数值
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    if(is.na(pK)) pK <- 0.09
    cat(sprintf("Using pK: %.3f\n", pK))
    
    # 运行DoubletFinder
    current_obj <- doubletFinder(
      current_obj,
      PCs = 1:30,
      pN = 0.25,
      pK = pK,
      nExp = nExp_poi,
      sct = FALSE
    )
    
    # 获取classification列名
    DF.name <- colnames(current_obj@meta.data)[grep("DF.classifications", 
                                                    colnames(current_obj@meta.data))]
    
    # 保存doublet结果
    if(length(DF.name) > 0) {
      current_obj$doublet_status <- current_obj@meta.data[[DF.name[1]]]
      doublet_count <- table(current_obj$doublet_status)
      cat("Doublet detection results:\n")
      print(doublet_count)
    } else {
      current_obj$doublet_status <- "Unknown"
      cat("No doublet classification found\n")
    }
    
  }, error = function(e) {
    warning(sprintf("DoubletFinder failed for sample %s: %s", sample_name, e$message))
    current_obj$doublet_status <- "Unknown"
  })
  
  # 保存结果
  sample_results[[sample_name]] <- current_obj$doublet_status
}
# 设置内存优化参数


# 加载必要的包
library(Seurat)
library(harmony)
library(future)
library(future.apply)
library(lisi)
library(ggplot2)

# 设置并行计算
plan("multiprocess", workers = 4)

# 定义回归异常值检测函数
regression_outliers <- function(seurat_obj, outliers_threshold = 0.999) {
  # 获取metrics数据
  nCount <- seurat_obj$nCount_RNA
  nFeature <- seurat_obj$nFeature_RNA
  
  # 对数转换
  log_counts <- log10(nCount + 1)
  log_features <- log10(nFeature + 1)
  
  # 线性回归
  fit <- lm(log_features ~ log_counts)
  
  # 计算预测区间
  pred <- predict(fit, 
                  data.frame(log_counts = log_counts), 
                  interval = "prediction", 
                  level = outliers_threshold
  )
  
  # 标记异常值
  outliers <- log_features < pred[,"lwr"] | log_features > pred[,"upr"]
  
  # 创建QC图目录
  if(!dir.exists("qaqc_plots")) {
    dir.create("qaqc_plots")
  }
  
  # 可视化
  pdf(file.path("qaqc_plots", 
                paste0("regression_outliers_", 
                       paste(unique(seurat_obj$sample), collapse = "_"), 
                       ".pdf")),
      width = 10, height = 8
  )
  
  # 绘制散点图
  plot(log_counts, log_features,
       pch = 16, cex = 0.6,
       col = ifelse(outliers, "red", "grey50"),
       xlab = "log10(nCount_RNA)",
       ylab = "log10(nFeature_RNA)",
       main = "Regression Outliers Detection"
  )
  
  # 添加拟合线和预测区间
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
         lwd = c(NA, NA, 2, 1)
  )
  
  dev.off()
  
  # 输出结果
  cat(sprintf("Detected %d outliers (%.1f%%)\n", 
              sum(outliers), 
              mean(outliers) * 100))
  
  # 返回非异常值的细胞
  cells.keep <- colnames(seurat_obj)[!outliers]
  return(cells.keep)
}

# 1. 整合doublet检测结果
print("Integrating doublet detection results...")
seurat_obj_main$doublet_status <- "Unknown"
for(sample_name in names(sample_results)) {
  cells_in_sample <- WhichCells(seurat_obj_main, expression = sample == sample_name)
  seurat_obj_main$doublet_status[cells_in_sample] <- sample_results[[sample_name]]
}

# 统计和过滤doublets
cells_before_doublet <- ncol(seurat_obj_main)
seurat_obj_main <- subset(seurat_obj_main, subset = doublet_status != "Doublet")
cells_after_doublet <- ncol(seurat_obj_main)
cat(sprintf("\nRemoved %d doublets (%.1f%%)\n", 
            cells_before_doublet - cells_after_doublet, 
            (cells_before_doublet - cells_after_doublet)/cells_before_doublet * 100))

# 2. 回归异常值检测
print("\nPerforming regression outlier detection...")
cells_to_keep <- regression_outliers(seurat_obj_main, outliers_threshold = 0.999)
cells_before_regression <- ncol(seurat_obj_main)
seurat_obj_main <- subset(seurat_obj_main, cells = cells_to_keep)
cells_after_regression <- ncol(seurat_obj_main)

cat(sprintf("\nRegression outlier removal stats:"))
cat(sprintf("\n- Cells before: %d", cells_before_regression))
cat(sprintf("\n- Cells after: %d", cells_after_regression))
cat(sprintf("\n- Removed cells: %d (%.1f%%)\n", 
            cells_before_regression - cells_after_regression,
            (cells_before_regression - cells_after_regression)/cells_before_regression * 100))

# 3. 数据预处理
print("\nPreprocessing data...")
DefaultAssay(seurat_obj_main) <- "RNA"

# 数据标准化和变异基因选择
print("Normalizing data...")
seurat_obj_main <- NormalizeData(seurat_obj_main, verbose = TRUE)

print("Finding variable features...")
seurat_obj_main <- FindVariableFeatures(seurat_obj_main, 
                                        selection.method = "vst", 
                                        nfeatures = 3000,
                                        verbose = TRUE)

print("Scaling data for all cells together...")
seurat_obj_main <- ScaleData(seurat_obj_main,
                             features = VariableFeatures(seurat_obj_main),
                             verbose = TRUE)

# 4. 降维分析
print("Running PCA...")
seurat_obj_main <- RunPCA(seurat_obj_main, npcs = 40, verbose = TRUE)

# 5. Harmony整合
print("Running Harmony integration...")
seurat_obj_main <- RunHarmony(
  object = seurat_obj_main,
  group.by.vars = c("sample", "tissue_sampling_method"),
  theta = c(2, 1),
  lambda = c(5, 5),
  sigma = 0.15,
  nclust = 60,
  max_iter = 30,
  reduction.use = "pca",
  dims.use = 1:40,
  verbose = TRUE
)

# 6. 计算整合质量指标
print("Calculating integration quality metrics...")
set.seed(42)
n_cells_sample <- min(10000, ncol(seurat_obj_main))
n_iterations <- 10
k <- 30

valid_cells <- !is.na(seurat_obj_main$ann_level_2)
valid_embeddings <- Embeddings(seurat_obj_main, "harmony")[valid_cells, ]
valid_metadata <- seurat_obj_main@meta.data[valid_cells, ]

lisi_results <- list()
for(i in 1:n_iterations) {
  set.seed(i)
  sampled_cells <- sample(nrow(valid_embeddings), n_cells_sample)
  
  lisi_res <- lisi::compute_lisi(valid_embeddings[sampled_cells, ],
                                 valid_metadata[sampled_cells, ],
                                 c("sample", "ann_level_2"),
                                 k)
  
  lisi_results[[i]] <- list(
    ilisi = mean(lisi_res[, "sample"]),
    clisi = mean(lisi_res[, "ann_level_2"])
  )
}

ilisi_scores <- sapply(lisi_results, function(x) x$ilisi)
clisi_scores <- sapply(lisi_results, function(x) x$clisi)

cat("\nIntegration Quality Metrics:")
cat(sprintf("\niLISI Score: %.3f ± %.3f", mean(ilisi_scores), sd(ilisi_scores)))
cat(sprintf("\ncLISI Score: %.3f ± %.3f", mean(clisi_scores), sd(clisi_scores)))

# 7. UMAP降维
print("Running UMAP...")
seurat_obj_main <- RunUMAP(seurat_obj_main,
                           reduction = "harmony",
                           dims = 1:40,
                           n.neighbors = 15,
                           min.dist = 0.3,
                           learning.rate = 0.5,
                           n.epochs = 1200,
                           spread = 1.2,
                           metric = "correlation",
                           verbose = TRUE)

# 8. 聚类分析
print("Finding neighbors and clusters...")
seurat_obj_main <- FindNeighbors(seurat_obj_main,
                                 reduction = "harmony",
                                 dims = 1:40,
                                 k.param = 15,
                                 prune.SNN = 1/20,
                                 verbose = TRUE)

seurat_obj_main <- FindClusters(seurat_obj_main,
                                resolution = 3,
                                verbose = TRUE)

# 9. 可视化
print("Generating visualizations...")

# 创建可视化目录
if(!dir.exists("visualization_plots")) {
  dir.create("visualization_plots")
}

# 基本UMAP图
pdf("visualization_plots/umap_plots.pdf", width = 15, height = 15)

# Batch effect相关可视化
p1 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "study", 
              raster = TRUE, pt.size = 0.5) + 
  ggtitle("Studies") + 
  theme_minimal()

p2 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "sample", 
              raster = TRUE, pt.size = 0.5) + 
  ggtitle("Samples") + 
  theme_minimal()

p3 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "ann_level_2", 
              raster = TRUE, pt.size = 0.5) + 
  ggtitle("Cell Types") + 
  theme_minimal()

print(p1)
print(p2)
print(p3)

# 分Study展示
for(study in unique(seurat_obj_main$study)) {
  if(!is.na(study)) {
    p <- DimPlot(seurat_obj_main,
                 reduction = "umap",
                 cells.highlight = WhichCells(seurat_obj_main, 
                                              expression = study == study),
                 cols = c("grey", "red"),
                 pt.size = 0.5,
                 raster = TRUE) +
      ggtitle(study) +
      theme_minimal()
    print(p)
  }
}

# 分细胞类型展示
ann_level_2 <- unique(seurat_obj_main$ann_level_2)
ann_level_2 <- ann_level_2[!is.na(ann_level_2)]

for(anno in ann_level_2) {
  p <- DimPlot(seurat_obj_main,
               reduction = "umap",
               cells.highlight = WhichCells(seurat_obj_main, 
                                            expression = ann_level_2 == anno),
               cols = c("grey", "red"),
               pt.size = 0.5,
               raster = TRUE) +
    ggtitle(anno) +
    theme_minimal()
  print(p)
}

dev.off()

# 10. 保存结果
print("Saving results...")
saveRDS(seurat_obj_main, "final_processed_seurat_object.rds")

integration_metrics <- list(
  cells_before_doublet_removal = cells_before_doublet,
  cells_after_doublet_removal = cells_after_doublet,
  doublet_removal_percentage = (cells_before_doublet - cells_after_doublet)/cells_before_doublet * 100,
  cells_before_regression = cells_before_regression,
  cells_after_regression = cells_after_regression,
  regression_removal_percentage = (cells_before_regression - cells_after_regression)/cells_before_regression * 100,
  ilisi_mean = mean(ilisi_scores),
  ilisi_sd = sd(ilisi_scores),
  clisi_mean = mean(clisi_scores),
  clisi_sd = sd(clisi_scores)
)

saveRDS(integration_metrics, "integration_quality_metrics.rds")

print("Analysis completed successfully!")



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
DotPlot(seurat_obj_main, 
        features = Markers, 
        group.by = "seurat_clusters",
        cols = c("lightgrey", "red"),
        dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")
# 2. 计算每个cluster中每个基因的表达情况
print("Calculating expression statistics...")

# 针对 Seurat V5 获取表达矩阵
data.use <- GetAssayData(seurat_obj_main, layer = "data")  # V5使用layer而不是slot
Idents(seurat_obj_main) <- seurat_obj_main$seurat_clusters
# 初始化结果列表
results <- list()

# 对每个cluster计算
for(cluster in unique(Idents(seurat_obj_main))) {
  # 获取该cluster的细胞
  cells_in_cluster <- WhichCells(seurat_obj_main, idents = cluster)

  # 对每个marker基因计算
  for(gene in Markers) {
    if(gene %in% rownames(seurat_obj_main)) {
      # 获取该基因在该cluster中的表达数据
      expr_data <- data.use[gene, cells_in_cluster]

      # 计算检测率（百分比）
      pct_expressing <- sum(expr_data > 0) / length(expr_data) * 100

      # 计算平均表达量（所有细胞，包括不表达的）
      avg_expr <- mean(expr_data)

      # 添加到结果中
      results[[length(results) + 1]] <- data.frame(
        Cluster = cluster,
        Gene = gene,
        Average_Expression = round(avg_expr, 3),
        Percent_Expressing = round(pct_expressing, 2)
      )
    }
  }
}

# 合并所有结果
results_df <- do.call(rbind, results)

# 重塑数据成矩阵格式
expr_matrix <- reshape2::dcast(results_df, Gene ~ Cluster, value.var = "Average_Expression")
pct_matrix <- reshape2::dcast(results_df, Gene ~ Cluster, value.var = "Percent_Expressing")

# 保存结果
print("Saving results...")
write.csv(results_df, "cluster_marker_stats.csv", row.names = FALSE)
write.csv(expr_matrix, "cluster_marker_expression_matrix.csv", row.names = TRUE)
write.csv(pct_matrix, "cluster_marker_percent_matrix.csv", row.names = TRUE)

saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")
seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
markers <- FindAllMarkers(seurat_obj_main, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

table(subset(seurat_obj_main, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

# 11. 生成热图
print("Generating heatmap...")
top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)


# 1. 先缩放marker基因
marker_genes <- unique(top10_markers$gene)
seurat_obj_main <- ScaleData(seurat_obj_main, features = marker_genes)

# 1. 将marker基因分成较小的批次
batch_size <- 25  # 每批50个基因
marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))

# 生成文件名
filename <- "marker_heatmap_batch.pdf"
pdf(filename, width = 100, height = 750)
# 2. 为每个批次生成热图
for(i in seq_along(marker_batches)) {
  print(paste("Processing batch", i, "of", length(marker_batches)))
  
  # 绘制热图
  print(DoHeatmap(seurat_obj_main, 
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
annotations <- read.csv("annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "cell_ann_me")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- seurat_obj_main$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$cell_ann_me[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
seurat_obj_main$cell_ann_me <- new_idents
seurat_obj_main <- SetIdent(seurat_obj_main, value = new_idents)

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
DimPlot(seurat_obj_main, 
        reduction = "umap",
        raster = TRUE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        cols = cell_colors) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 6. 绘制分割版本
pdf("cell_types_umap_split.pdf", width = 12, height = 12)
# 获取所有细胞类型
cell_types <- unique(seurat_obj_main$cell_ann_me)
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 为每个细胞类型单独绘图
for(cell_type in cell_types) {
  # 创建文件名
  
  
  # 绘制该细胞类型的UMAP图
  
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = WhichCells(seurat_obj_main, idents = cell_type),
                cols.highlight = cell_colors[cell_type],
                cols = "grey",  # 其他细胞显示为灰色
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(cell_type) +
          theme(legend.text = element_text(size = 12)))
  
}
dev.off()
# 方法1
table(seurat_obj_main@meta.data$tissue[seurat_obj_main@meta.data$cell_ann_me == "B"])

# 1. 首先过滤数据，创建用于比例计算的子集
filtered_metadata <- seurat_obj_main@meta.data %>%
  filter(tissue_sampling_method != "scraping")

# 使用过滤后的数据计算细胞数量
cell_counts <- table(filtered_metadata$sample, 
                     filtered_metadata$cell_ann_me) %>%
  as.data.frame() %>%
  tidyr::spread(Var2, Freq) %>%
  setNames(c("Sample", setdiff(colnames(.), "Var1")))

# 2. 计算每个sample的总细胞数
total_cells <- rowSums(cell_counts[,-1])

# 3. 计算比例
proportion_cols <- setdiff(colnames(cell_counts), "Sample")
cell_proportions <- cell_counts
cell_proportions[,proportion_cols] <- sweep(cell_counts[,proportion_cols], 1, total_cells, "/")

# 4. 添加总细胞数列
cell_proportions$Total_Cells <- total_cells

# 5. 获取组织信息 (使用过滤后的数据)
tissue_data <- filtered_metadata %>%
  select(sample, tissue) %>%
  unique() %>%
  as.data.frame()

# 6. 合并组织信息
combined_data <- merge(tissue_data, cell_proportions, by.x = "sample", by.y = "Sample")

# 7. 将原始数量列名修改为带Count的列名
count_cols <- paste0(proportion_cols, "_Count")
combined_data[,count_cols] <- cell_counts[,proportion_cols]

# 8. 将比例转换为百分比
for(col in proportion_cols) {
  combined_data[[col]] <- combined_data[[col]] * 100
}

# 9. 设置列的顺序
cell_type_pairs <- c()
for(cell_type in proportion_cols) {
  cell_type_pairs <- c(cell_type_pairs, paste0(cell_type, "_Count"), cell_type)
}

col_order <- c("sample", "tissue", "Total_Cells", cell_type_pairs)
combined_data <- combined_data[, col_order]

# 10. 按照tissue和Sample排序
combined_data <- combined_data[order(combined_data$tissue, combined_data$sample),]

# 11. 输出到CSV
write.csv(combined_data, "cell_type_counts_and_proportions.csv", row.names = FALSE)

# 12. 为绘图准备数据
tissue_summary <- combined_data %>%
  group_by(tissue) %>%
  summarise(across(proportion_cols, list(
    mean = ~mean(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE)
  )))

tissue_summary_long <- tissue_summary %>%
  tidyr::pivot_longer(
    cols = -tissue,
    names_to = c("cell_ann_me", "Stat"),
    names_pattern = "(.*)_(mean|sd)",
    values_to = "Value"
  ) %>%
  tidyr::pivot_wider(
    names_from = Stat,
    values_from = Value
  )

# 修复 SeuratObject 分析流程

# 1. 首先进行统计检验
stat_results <- list()
cell_types <- unique(seurat_obj_main$cell_ann_me)
cell_types <- cell_types[!is.na(cell_types)]

# 对每个细胞类型进行组织间的比较
for(cell_type in proportion_cols) {
  # 创建用于统计分析的数据框
  test_data <- data.frame(
    Proportion = combined_data[[cell_type]],
    tissue = combined_data$tissue
  )
  
  # 进行 Kruskal-Wallis 检验
  kw_test <- kruskal.test(Proportion ~ tissue, data = test_data)
  
  if(kw_test$p.value < 0.05) {
    # 如果显著，进行两两比较
    tissue_pairs <- combn(unique(test_data$tissue), 2, simplify = FALSE)
    pair_results <- list()
    
    for(pair in tissue_pairs) {
      pair_data <- test_data[test_data$tissue %in% pair,]
      wilcox_test <- wilcox.test(
        Proportion ~ tissue,
        data = pair_data,
        exact = FALSE
      )
      
      pair_results[[length(pair_results) + 1]] <- data.frame(
        cell_ann_me = cell_type,
        tissue1 = pair[1],
        tissue2 = pair[2],
        pvalue = wilcox_test$p.value
      )
    }
    
    # 合并所有配对的结果
    pair_results_df <- do.call(rbind, pair_results)
    
    # 进行多重检验校正
    pair_results_df$qvalue <- p.adjust(pair_results_df$pvalue, method = "BH")
    
    # 存储结果
    stat_results[[cell_type]] <- pair_results_df
  }
}

# 合并所有细胞类型的结果
all_stat_results <- do.call(rbind, stat_results)

# 2. 修复箱线图绘制部分，添加显著性标记
for(cell_type in proportion_cols) {
  # 获取这个细胞类型的显著性结果
  sig_results <- all_stat_results[all_stat_results$cell_ann_me == cell_type,]
  
  # 基础箱线图
  p <- ggplot(combined_data, aes(x = tissue, y = .data[[cell_type]], fill = tissue)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(x = "tissue", 
         y = "Percentage (%)", 
         title = paste(cell_type, "Distribution")) +
    scale_fill_brewer(palette = "Set3")
  
  # 获取y轴范围用于设置标注位置
  y_max <- max(combined_data[,cell_type], na.rm = TRUE)
  y_range <- y_max - min(combined_data[,cell_type], na.rm = TRUE)
  
  # 仅为显著的比较添加标注
  sig_pairs <- sig_results[sig_results$qvalue < 0.05,]
  
  if(nrow(sig_pairs) > 0) {
    for(i in 1:nrow(sig_pairs)) {
      # 获取x轴位置
      x1 <- which(unique(combined_data$tissue) == sig_pairs$tissue1[i])
      x2 <- which(unique(combined_data$tissue) == sig_pairs$tissue2[i])
      
      # 设置标注高度
      y_pos <- y_max + y_range * (0.1 + 0.1 * i)
      
      # 添加连接线
      p <- p + annotate("segment", 
                        x = x1, xend = x2,
                        y = y_pos, yend = y_pos,
                        colour = "black")
      
      # 添加 q 值和星号标注
      q_value <- format(sig_pairs$qvalue[i], digits = 3, scientific = TRUE)
      stars <- ifelse(sig_pairs$qvalue[i] < 0.001, "***",
                      ifelse(sig_pairs$qvalue[i] < 0.01, "**",
                             ifelse(sig_pairs$qvalue[i] < 0.05, "*", "")))
      
      p <- p + annotate("text",
                        x = (x1 + x2)/2,
                        y = y_pos + y_range * 0.02,
                        label = paste0("q = ", q_value, " ", stars),
                        size = 3)
    }
    
    # 调整y轴范围以适应标注
    max_y_with_anno <- y_max + y_range * (0.1 + 0.1 * nrow(sig_pairs) + 0.1)
    p <- p + coord_cartesian(ylim = c(NA, max_y_with_anno))
  }
  
  # 保存图片
  pdf(paste0(cell_type, "_statistical_boxplot.pdf"), width = 10, height = 8)
  print(p)
  dev.off()
}

# 3. 创建详细的统计结果摘要表
stat_summary <- all_stat_results %>%
  # 保留所有结果，不仅仅是显著的
  arrange(cell_ann_me, qvalue) %>%
  mutate(
    # 保持原始组织名称的顺序
    Comparison = paste(tissue1, "vs", tissue2),
    # 格式化 p 值和 q 值
    pvalue_formatted = format(pvalue, digits = 3, scientific = TRUE),
    qvalue_formatted = format(qvalue, digits = 3, scientific = TRUE),
    # 添加显著性标记
    Significance = case_when(
      qvalue < 0.001 ~ "***",
      qvalue < 0.01 ~ "**",
      qvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    # 添加是否显著的标记
    Significant = qvalue < 0.05,
    # 添加效应大小的方向（如果需要的话，可以根据实际数据调整）
    Direction = case_when(
      Significant == FALSE ~ "Not Significant",
      TRUE ~ "Significant"
    )
  ) %>%
  # 重新排列列的顺序
  select(
    cell_ann_me, 
    Comparison, 
    tissue1, 
    tissue2, 
    pvalue,
    pvalue_formatted,
    qvalue,
    qvalue_formatted,
    Significance,
    Significant,
    Direction
  )

# 4. 输出统计结果
# 保存完整的统计结果
write.csv(stat_summary, "tissue_comparison_all_results.csv", row.names = FALSE)

# 保存仅显著性结果
write.csv(
  stat_summary %>% filter(Significant == TRUE),
  "tissue_comparison_significant_results.csv", 
  row.names = FALSE
)

# 5. 保存修改后的 Seurat 对象
saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")

print("Analysis completed successfully!")

# 假设我们选择 cell_type 为 "T_cells" 的细胞
# Epithelial_obj <- subset(seurat_obj_main, subset = cell_ann_me == "Epithelial")
# Epithelial_obj <- NormalizeData(Epithelial_obj)
# Epithelial_obj <- FindVariableFeatures(Epithelial_obj, selection.method = "vst", nfeatures = 3000)
# Epithelial_obj <- ScaleData(Epithelial_obj, features = VariableFeatures(Epithelial_obj))
# Epithelial_obj <- RunPCA(Epithelial_obj)
# 
# pca_embeddings <- Embeddings(Epithelial_obj, 'pca')
# metadata_subset <- data.frame(
#   sample = Epithelial_obj@meta.data$sample,
#   tissue = Epithelial_obj@meta.data$tissue,
#   row.names = rownames(Epithelial_obj@meta.data)
# )
# 
# harmony_out <- RunHarmony(
#   data_mat = pca_embeddings,           
#   meta_data = metadata_subset,
#   vars_use = c("sample", "tissue"),   
#   theta = c(2, 1.5),                   
#   lambda = c(1, 1),                    
#   sigma = 0.2,                         
#   nclust = 50,                         
#   max_iter = 20,                       
#   early_stop = TRUE,                   
#   plot_convergence = TRUE,             
#   verbose = TRUE                       
# )
# 
# # RPCA处理
# harmony_embeddings <- Embeddings(Epithelial_obj, 'harmony')
# rpca_result <- rpca(harmony_embeddings,
#                     k = 50,
#                     center = TRUE,
#                     scale = TRUE,
#                     retx = TRUE,
#                     p = 10,
#                     q = 2,
#                     rand = TRUE)
# 
# # 添加RPCA结果
# colnames(rpca_result$x) <- paste0("RPCA_", 1:ncol(rpca_result$x))
# rownames(rpca_result$x) <- colnames(Epithelial_obj)
# Epithelial_obj[["rpca"]] <- CreateDimReducObject(
#   embeddings = rpca_result$x,
#   key = "RPCA_",
#   assay = DefaultAssay(Epithelial_obj)
# )
# # 7. UMAP降维和聚类
# print("Running UMAP and clustering...")
# Epithelial_obj <- RunUMAP(Epithelial_obj, reduction = "rpca", dims = 1:30, metric = "correlation")
# Epithelial_obj <- FindNeighbors(Epithelial_obj, reduction = "rpca", dims = 1:30)
# Epithelial_obj <- FindClusters(Epithelial_obj, resolution = 2)
# 
# print("Analysis completed!")
# 
# 
# # 9. 保存分析结果
# print("Saving analysis results...")
# saveRDS(Epithelial_obj, "final_analyzed_Epithelial_obj.rds")
# 
# 
# # 8. 绘制UMAP聚类图
# print("Generating plots...")
# pdf("Epithelial_umap_clusters.pdf", width = 10, height = 8)
# DimPlot(Epithelial_obj, reduction = "umap", label = TRUE)

cell_types <- cell_types[! cell_types %in% c('Proliferation')]
for (cell_type in cell_types) {
  new_var_name <- paste(cell_type, 'obj', sep = "_")
  Marker_current <- get(paste0('Marker_',cell_type))
  assign(new_var_name, run_subset_analysis(
    seurat_obj_main,
    cell_type = cell_type,
    markers = Marker_current,
    output_dir = paste(cell_type, "analysis", sep = "_")
  ))
}


