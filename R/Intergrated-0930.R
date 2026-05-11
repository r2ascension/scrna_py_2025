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
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(CellChat)
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = +Inf)

# reticulate::py_install(packages = 'umap-learn')
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
Marker_Stromal <- c('APOD','FGF7','COL15A1','MFAP5','PI16','CD34',
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
                    'CXCL14','VSTM2A','SOX6','COL4A5','COL4A6','TSLP','FRZB','BMP5','BMP2','CPM','F3',
                    'RGS5','CD36','NOTCH3',
                    'SCIN','EPAS1','HLA.C','IGKC','PTP4A3','FN1','COL18A1','WFDC1','IGHG4','IGHG1',
                    'RERGL','PLN','SORBS2','DSTN','TSC22D1','BCAM','C11orf96','FILIP1','WTIP','NRGN',
                    'TMEM176B','ANGPTL1','CFH','FHL1','VCAN','GGT5','C1S','COL6A3','STEAP4','ADGRL3')
Marker_Immune <- c('CLEC9A','XCR1','CADM1','CLNK','FLT3','ZBTB46',
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
                   'FCGR3B','CSF3R','CXCR1',
                   'KLRD1','FCGR3A','GNLY','TYROBP','FCER1G','KLRC1','FGFBP2','SPON2','MYOM2',
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
                   'MKI67','TOP2A','TK1','CENPW',
                   'IGHD','TCL1A',
                   'MS4A1','BANK1',
                   'MKI67','TOP2A',
                   'IGHA1','IGHA2',
                   'IGHGP','IGHG1')
setwd("E:/R/0301/")
source("E:/R/enrichment_functions.R")
source("E:/R/scMASC.R")


###########################################################################################################
###########################################################################################################
###########################################################################################################
# seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main_0630.rds")
# seurat_obj2 <- readRDS("E:/R/Source/final/Weiqing_Wang_2022.rds")
seurat_obj_main <- merge(x = seurat_obj_main, 
                       y = seurat_obj2,
                       project = "Merged_Project")

cat("过滤前细胞数:", ncol(seurat_obj_main), "\n")
# cat("过滤前细胞数:", ncol(seurat_obj2), "\n")
seurat_obj_main <- subset(seurat_obj_main, subset = nFeature_RNA <= 6000 & nFeature_RNA >= 200)
# seurat_obj2 <- subset(seurat_obj2, subset = nFeature_RNA <= 6000 & nFeature_RNA >= 200)
cat("过滤后细胞数:", ncol(seurat_obj_main), "\n")

samples_to_remove <- c("F02609", "F01506", "GRO-02_biopsy", "GRO-01_biopsy")
seurat_obj_main <- subset(
  seurat_obj_main,
  subset = sample %in% samples_to_remove,
  invert = TRUE
)

# =====================================================================
# 方法1: 完全重建Assay（最稳定）
# =====================================================================

rebuild_assay_from_scratch <- function(seurat_obj, assay = "RNA") {
  #' 从Assay5完全重建v3/v4兼容的Assay
  #' Completely rebuild v3/v4 compatible Assay from Assay5
  #' 
  #' @param seurat_obj Seurat对象
  #' @param assay assay名称
  #' @return 重建后的Seurat对象
  
  cat("=== 开始重建Assay ===\n")
  cat("=== Starting Assay Reconstruction ===\n\n")
  
  assay_obj <- seurat_obj[[assay]]
  
  if (!inherits(assay_obj, "Assay5")) {
    cat("当前已是Assay格式，无需重建\n")
    return(seurat_obj)
  }
  
  # 步骤1: 提取counts数据
  cat("步骤1/3: 提取counts数据...\n")
  cat("Step 1/3: Extracting counts data...\n")
  
  tryCatch({
    # 尝试从counts层提取
    counts_data <- LayerData(seurat_obj, assay = assay, layer = "counts")
    cat(sprintf("  Counts矩阵: %d genes × %d cells\n", 
                nrow(counts_data), ncol(counts_data)))
  }, error = function(e) {
    stop("无法提取counts数据: ", e$message)
  })
  
  # 步骤2: 提取normalized data（如果存在）
  cat("步骤2/3: 提取normalized data...\n")
  cat("Step 2/3: Extracting normalized data...\n")
  
  tryCatch({
    data_data <- LayerData(seurat_obj, assay = assay, layer = "data")
    has_data <- TRUE
    cat(sprintf("  Data矩阵: %d genes × %d cells\n", 
                nrow(data_data), ncol(data_data)))
  }, error = function(e) {
    cat("  未找到data层，将从counts重新计算\n")
    data_data <- NULL
    has_data <- FALSE
  })
  
  # 步骤3: 创建新的Assay对象
  cat("步骤3/3: 创建新的Assay对象...\n")
  cat("Step 3/3: Creating new Assay object...\n")
  
  # 创建基础Assay（只包含counts）
  new_assay <- CreateAssayObject(counts = counts_data)
  
  # 如果有normalized data，添加进去
  if (has_data) {
    new_assay <- SetAssayData(new_assay, 
                              slot = "data", 
                              new.data = data_data)
    cat("  ✓ 已添加normalized data\n")
  } else {
    # 如果没有data，进行标准化
    cat("  正在进行标准化...\n")
    temp_obj <- CreateSeuratObject(counts = counts_data)
    temp_obj <- NormalizeData(temp_obj, verbose = FALSE)
    new_assay <- SetAssayData(new_assay, 
                              slot = "data",
                              new.data = GetAssayData(temp_obj, slot = "data"))
    rm(temp_obj)
    cat("  ✓ 已完成标准化\n")
  }
  
  # 替换原assay
  seurat_obj[[assay]] <- new_assay
  
  cat("\n✓ Assay重建完成!\n")
  cat("✓ Assay reconstruction completed!\n\n")
  
  # 验证
  cat("验证结果:\n")
  cat("Validation:\n")
  cat(sprintf("  - Assay类型: %s\n", class(seurat_obj[[assay]])[1]))
  cat(sprintf("  - Counts: %d genes × %d cells\n", 
              nrow(GetAssayData(seurat_obj, slot = "counts")),
              ncol(GetAssayData(seurat_obj, slot = "counts"))))
  cat(sprintf("  - Data: %d genes × %d cells\n", 
              nrow(GetAssayData(seurat_obj, slot = "data")),
              ncol(GetAssayData(seurat_obj, slot = "data"))))
  
  return(seurat_obj)
}

# =====================================================================
# 方法2: 保守的JoinLayers + 转换
# =====================================================================

safe_convert_assay <- function(seurat_obj, assay = "RNA") {
  #' 安全的Assay5到Assay转换
  #' Safe Assay5 to Assay conversion
  
  cat("=== 安全转换模式 ===\n")
  
  assay_obj <- seurat_obj[[assay]]
  
  if (!inherits(assay_obj, "Assay5")) {
    cat("已是Assay格式\n")
    return(seurat_obj)
  }
  
  # 步骤1: 检查并列出所有层
  cat("\n当前Assay5的层:\n")
  print(Layers(seurat_obj, assay = assay))
  
  # 步骤2: 合并层
  cat("\n合并所有层...\n")
  seurat_obj[[assay]] <- JoinLayers(seurat_obj[[assay]])
  
  cat("合并后的层:\n")
  print(Layers(seurat_obj, assay = assay))
  
  # 步骤3: 转换
  cat("\n转换为Assay格式...\n")
  seurat_obj[[assay]] <- as(seurat_obj[[assay]], Class = "Assay")
  
  cat("✓ 转换完成!\n")
  cat(sprintf("最终Assay类型: %s\n", class(seurat_obj[[assay]])[1]))
  
  return(seurat_obj)
}

# =====================================================================
# 方法3: 逐样本处理时临时转换（内存友好）
# =====================================================================

process_sample_with_conversion <- function(seurat_obj, 
                                           sample_name,
                                           sample_col = "sample") {
  #' 提取单个样本并转换为Assay格式
  #' Extract single sample and convert to Assay format
  #' 
  #' @param seurat_obj 主Seurat对象（Assay5格式）
  #' @param sample_name 样本名称
  #' @param sample_col 样本列名
  #' @return 转换后的样本Seurat对象
  
  cat(sprintf("\n提取样本: %s\n", sample_name))
  
  # 提取样本
  cells <- WhichCells(seurat_obj, 
                      expression = get(sample_col) == sample_name)
  
  sample_obj <- subset(seurat_obj, cells = cells)
  
  # 获取数据
  counts <- LayerData(sample_obj, assay = "RNA", layer = "counts")
  
  # 创建新对象（自动为Assay格式）
  sample_obj_new <- CreateSeuratObject(
    counts = counts,
    meta.data = sample_obj@meta.data
  )
  
  cat(sprintf("  细胞数: %d\n", ncol(sample_obj_new)))
  cat(sprintf("  Assay类型: %s\n", class(sample_obj_new[["RNA"]])[1]))
  
  return(sample_obj_new)
}

# =====================================================================
# 使用示例
# =====================================================================

# 选择一种方法使用

# === 方法1: 完全重建（最推荐） ===
# seurat_obj_main <- rebuild_assay_from_scratch(seurat_obj_main)

# # 1. 查看现有的降维结果
# names(seurat_obj_main@reductions)
# # 例如: "pca" "umap" "tsne" "harmony"
# 
# # 2. 删除特定的降维
# seurat_obj_main[["pca"]] <- NULL
# seurat_obj_main[["umap"]] <- NULL
# seurat_obj_main[["harmony"]] <- NULL
# seurat_obj_main@assays$RNA@scale.data <- matrix(0, 0, 0)
# # 3. 删除所有降维（批量）
# for(i in names(seurat_obj_main@reductions)){
#   seurat_obj_main[[i]] <- NULL
# }
# 
# # 4. 验证
# names(seurat_obj_main@reductions)

###########################################################################################################
###########################################################################################################
###########################################################################################################
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

cat("\n")
cat("============================================================\n")
cat("   Single-cell Doublet Detection Pipeline\n")
cat("   Platform: 10X Chromium | Expected rate: 15%\n")
cat("============================================================\n\n")


# ---------- 核心函数 ----------

#' 单样本双细胞检测
process_single_sample <- function(seurat_obj, sample_name, doublet_rate = 0.15) {
  
  cat(sprintf("\n========== Sample: %s ==========\n", sample_name))
  
  n_cells <- ncol(seurat_obj)
  cat(sprintf("Cells: %d\n", n_cells))
  
  # 细胞数检查
  if(n_cells < 50) {
    cat("SKIP: Insufficient cells (<50)\n")
    seurat_obj$doublet_status <- "Insufficient_cells"
    seurat_obj$doublet_score <- NA
    seurat_obj$doublet_class <- "Skipped"
    return(seurat_obj)
  }
  
  # 计算预期双细胞数
  nExp <- round(doublet_rate * n_cells)
  cat(sprintf("Expected doublets: %d (%.1f%%)\n", nExp, doublet_rate * 100))
  
  tryCatch({
    # 标准化
    cat("Step 1/4: Normalizing...\n")
    if(is.null(seurat_obj@assays$RNA@data) || nrow(seurat_obj@assays$RNA@data) == 0) {
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    }
    
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), 
                            verbose = FALSE)
    
    # PCA
    cat("Step 2/4: Running PCA...\n")
    n_pcs <- min(30, ncol(seurat_obj) - 1)
    n_pcs <- max(n_pcs, 10)
    seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
    
    # pK优化
    cat("Step 3/4: Optimizing pK...\n")
    sweep.res <- paramSweep(seurat_obj, PCs = 1:n_pcs, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    if(is.na(pK) || pK < 0.01 || pK > 0.3) {
      pK <- 0.09
      cat("Using default pK: 0.09\n")
    } else {
      cat(sprintf("Optimal pK: %.3f\n", pK))
    }
    
    # DoubletFinder
    cat("Step 4/4: Running DoubletFinder...\n")
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:n_pcs, pN = 0.25, 
                                pK = pK, nExp = nExp, sct = FALSE)
    
    # 提取结果
    df_class <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
    df_score <- grep("pANN", colnames(seurat_obj@meta.data), value = TRUE)
    
    if(length(df_class) > 0) {
      seurat_obj$doublet_status <- seurat_obj@meta.data[[df_class[1]]]
      seurat_obj$doublet_score <- if(length(df_score) > 0) seurat_obj@meta.data[[df_score[1]]] else NA
      seurat_obj$doublet_class <- ifelse(seurat_obj$doublet_status == "Singlet", "Singlet", "Doublet")
      
      # 删除原始列
      seurat_obj@meta.data[, df_class] <- NULL
      seurat_obj@meta.data[, df_score] <- NULL
      
      n_singlets <- sum(seurat_obj$doublet_class == "Singlet")
      n_doublets <- sum(seurat_obj$doublet_class == "Doublet")
      cat(sprintf("Results: %d singlets (%.1f%%), %d doublets (%.1f%%)\n",
                  n_singlets, n_singlets/n_cells*100,
                  n_doublets, n_doublets/n_cells*100))
    } else {
      cat("WARNING: Detection failed\n")
      seurat_obj$doublet_status <- "Failed"
      seurat_obj$doublet_score <- NA
      seurat_obj$doublet_class <- "Failed"
    }
    
    seurat_obj$doublet_sample <- sample_name
    
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    seurat_obj$doublet_status <- "Error"
    seurat_obj$doublet_score <- NA
    seurat_obj$doublet_class <- "Error"
    seurat_obj$doublet_sample <- sample_name
  })
  
  return(seurat_obj)
}


#' 批量处理所有样本
process_all_samples <- function(seurat_obj, sample_column = "orig.ident", 
                                doublet_rate = 0.15, output_dir = "doublet_results") {
  
  # 创建输出目录
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 检查样本列
  if(!sample_column %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found. Available: %s",
                 sample_column, paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  # 样本列表
  samples <- unique(seurat_obj@meta.data[[sample_column]])
  n_samples <- length(samples)
  
  cat(sprintf("Total samples: %d\n", n_samples))
  cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))
  cat("\nSample distribution:\n")
  print(table(seurat_obj@meta.data[[sample_column]]))
  cat("\n")
  
  # 初始化
  processed_list <- list()
  stats <- data.frame()
  
  # 逐样本处理
  for(i in 1:n_samples) {
    sample_name <- samples[i]
    cat(sprintf("\n[%d/%d] Processing: %s\n", i, n_samples, sample_name))
    cat("------------------------------------------------------------\n")
    
    # 提取样本
    cells <- colnames(seurat_obj)[seurat_obj@meta.data[[sample_column]] == sample_name]
    obj_sample <- subset(seurat_obj, cells = cells)
    
    # 处理
    obj_sample <- process_single_sample(obj_sample, sample_name, doublet_rate)
    processed_list[[sample_name]] <- obj_sample
    
    # 统计
    stats <- rbind(stats, data.frame(
      Sample = sample_name,
      Total = ncol(obj_sample),
      Singlets = sum(obj_sample$doublet_class == "Singlet", na.rm = TRUE),
      Doublets = sum(obj_sample$doublet_class == "Doublet", na.rm = TRUE),
      Rate = mean(obj_sample$doublet_class == "Doublet", na.rm = TRUE) * 100,
      Status = ifelse("Singlet" %in% obj_sample$doublet_class, "Success",
                      unique(obj_sample$doublet_class)[1])
    ))
  }
  
  # 合并
  cat("\n============================================================\n")
  cat("Merging all samples...\n")
  merged <- merge(x = processed_list[[1]], 
                  y = processed_list[2:length(processed_list)],
                  merge.data = TRUE)
  
  # 统计摘要
  cat("\n=== Summary Statistics ===\n")
  print(stats)
  cat("\n")
  
  write.csv(stats, file.path(output_dir, "summary.csv"), row.names = FALSE)
  
  # 可视化
  cat("Generating plots...\n")
  
  # 图1: 双细胞率
  p1 <- ggplot(stats, aes(x = Sample, y = Rate, fill = Status)) +
    geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = doublet_rate * 100, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%.1f%%", Rate)), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Success" = "#3498DB", "Skipped" = "#95A5A6", "Error" = "#E74C3C")) +
    labs(title = "Doublet Rate by Sample", 
         subtitle = sprintf("Red line = Expected (%.1f%%)", doublet_rate * 100),
         x = "Sample", y = "Doublet Rate (%)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "doublet_rate.pdf"), p1, width = 10, height = 6)
  
  # 图2: 细胞分布
  stats_long <- stats %>%
    select(Sample, Singlets, Doublets) %>%
    tidyr::pivot_longer(cols = c(Singlets, Doublets), names_to = "Type", values_to = "Count")
  
  p2 <- ggplot(stats_long, aes(x = Sample, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    geom_text(data = stats, aes(x = Sample, y = Total, label = Total),
              vjust = -0.5, inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = c("Singlets" = "#2ECC71", "Doublets" = "#E74C3C")) +
    labs(title = "Cell Distribution", x = "Sample", y = "Cells") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "cell_distribution.pdf"), p2, width = 10, height = 6)
  
  # 图3: Score分布
  valid_data <- merged@meta.data[!is.na(merged$doublet_score), ]
  if(nrow(valid_data) > 0) {
    p3 <- ggplot(valid_data, aes(x = doublet_score, fill = doublet_class)) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(values = c("Singlet" = "#2ECC71", "Doublet" = "#E74C3C")) +
      labs(title = "Doublet Score Distribution", x = "Score", y = "Density") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(file.path(output_dir, "score_distribution.pdf"), p3, width = 8, height = 6)
  }
  
  # 保存对象
  saveRDS(merged, file.path(output_dir, "seurat_annotated.rds"))
  
  # 最终统计
  cat("\n============================================================\n")
  cat("Final Statistics:\n")
  cat(sprintf("Total cells: %d\n", ncol(merged)))
  cat(sprintf("Singlets: %d (%.1f%%)\n",
              sum(merged$doublet_class == "Singlet", na.rm = TRUE),
              mean(merged$doublet_class == "Singlet", na.rm = TRUE) * 100))
  cat(sprintf("Doublets: %d (%.1f%%)\n",
              sum(merged$doublet_class == "Doublet", na.rm = TRUE),
              mean(merged$doublet_class == "Doublet", na.rm = TRUE) * 100))
  cat("============================================================\n\n")
  
  cat(sprintf("Results saved to: %s/\n", output_dir))
  cat("  - summary.csv\n")
  cat("  - seurat_annotated.rds\n")
  cat("  - doublet_rate.pdf\n")
  cat("  - cell_distribution.pdf\n")
  cat("  - score_distribution.pdf\n\n")
  
  return(list(seurat = merged, stats = stats, samples = processed_list))
}

library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

# 假设您已经加载了 process_single_sample 和 process_all_samples 函数
# source("doubletfinder_functions.R")
# seurat_obj_main <- readRDS("E:/R/0301/seurat_raw_with_doublets_1003.rds")
# ------------------------------------------------------------
# 场景 1: 标准分析（无细胞注释）
# 推荐：探索性分析、首次分析
# ------------------------------------------------------------
# 运行双细胞检测
# 加载你的Seurat对象
seurat_obj_main <- readRDS("E:/R/0301/seurat_raw_with_doublets_1003.rds")
# 运行完整流程
result <- process_doublets_by_sample(
  seurat_obj = seurat_obj_main,
  sample_column = "sample",      # 根据实际情况修改
  doublet_rate = 0.1,               # 15%双细胞率
  output_dir = "doublet_qc_results"
)

# 获取处理后的对象
seurat_obj_main <- result$seurat_object

# 移除双细胞
seurat_obj_main <- subset(seurat_obj_main, subset = doublet_class == "Singlet")

# 4. 保存
saveRDS(seurat_obj_main, "seurat_raw_with_doublets_1003.rds")

###########################################################################################################
###########################################################################################################
###########################################################################################################


# 完整的单细胞RNA测序分析流程
# 包含所有函数和批量处理代码
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

#' 分析单个细胞类型子集
#'
#' @param seurat_obj Seurat对象
#' @param cell_type 细胞类型
#' @param cell_type_column 细胞类型列名（可选）
#' @param cell_type_value 细胞类型值（可选）
#' @param markers 标记基因列表
#' @param output_dir 输出目录
#' @param normalize 是否重新标准化
#' @return 处理后的子集Seurat对象
run_subset_analysis <- function(seurat_obj, cell_type = NULL, 
                                cell_type_column = NULL, cell_type_value = NULL,
                                markers = NULL, output_dir = ".", normalize = TRUE) {
  # 兼容两种调用方式
  if(!is.null(cell_type) && (is.null(cell_type_column) || is.null(cell_type_value))) {
    # 自动检测细胞类型列
    possible_columns <- c("Annotation")
    found_column <- NULL
    
    for(col in possible_columns) {
      if(col %in% colnames(seurat_obj@meta.data)) {
        if(cell_type %in% seurat_obj@meta.data[[col]]) {
          found_column <- col
          break
        }
      }
    }
    
    if(is.null(found_column)) {
      stop(paste("Error: Could not find cell type", cell_type, 
                 "in any of the metadata columns:", 
                 paste(possible_columns, collapse=", ")))
    }
    
    cell_type_column <- found_column
    cell_type_value <- cell_type
  } else if(is.null(cell_type_column) || is.null(cell_type_value)) {
    stop("Error: Either provide 'cell_type' or both 'cell_type_column' and 'cell_type_value'")
  }
  
  # 验证cell_type_column参数
  if(!cell_type_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Column", cell_type_column, "not found in Seurat object metadata"))
  }
  
  # 设置输出目录
  if(!is.null(markers)) markers <- unique(markers)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  sanitized_name <- gsub("[^a-zA-Z0-9_]", "_", cell_type_value)  # 安全的文件名
  output_prefix <- file.path(output_dir, paste0(sanitized_name, "_"))
  
  # 提取特定细胞类型的子集
  cat(paste("Extracting cells with", cell_type_column, "=", cell_type_value, "\n"))
  subset_cells <- seurat_obj@meta.data[[cell_type_column]] == cell_type_value
  if(sum(subset_cells) == 0) {
    stop(paste("Error: No cells found with", cell_type_column, "=", cell_type_value))
  }
  cat(paste("Found", sum(subset_cells), "cells\n"))
  
  subset_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[subset_cells])
  
  # 是否需要重新标准化
  # if(normalize) {
  #   cat("Normalizing data...\n")
  #   subset_obj <- NormalizeData(subset_obj)
  # }
  
  # 运行降维聚类
  cat("Running dimensionality reduction and clustering...\n")
  subset_obj <- run_dim_reduction(subset_obj, harmony = TRUE)
  
  # 运行分析
  cat("Running analysis and visualization...\n")
  subset_obj <- run_analysis(subset_obj, markers, output_prefix)
  
  # 保存结果
  cat("Saving results...\n")
  saveRDS(subset_obj, paste0(output_prefix, "analyzed.rds"))
  
  cat(paste("Analysis complete. Results saved to", output_dir, "\n"))
  return(subset_obj)
}

########################
# 批量分析主程序部分
########################
# 加载Seurat对象
# seurat_obj_main <- readRDS("path/to/your/seurat_obj_main.rds")
#' 优化的数据分析函数 - 内存友好版本，仅在可视化时抽样
#'
#' @param seurat_obj Seurat对象
#' @param markers 预设的marker基因列表，默认NULL
#' @param output_prefix 输出文件前缀
#' @param dot_plot_cols 点图的颜色设置，默认c("lightgrey", "red")
#' @param find_markers_params 寻找marker基因的参数列表
#' @param max_cells_for_heatmap 热图最大细胞数，默认3000
#' @return 处理后的Seurat对象
run_optimized_analysis <- function(seurat_obj, 
                                   markers = NULL, 
                                   output_prefix = "",
                                   dot_plot_cols = c("lightgrey", "red"),
                                   find_markers_params = list(only.pos = TRUE, 
                                                              min.pct = 0.25, 
                                                              logfc.threshold = 0.25),
                                   max_cells_for_heatmap = 3000) {
  
  # 确保输出目录存在
  output_dir <- dirname(output_prefix)
  if(output_dir != "" && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 首先输出UMAP聚类图
  message("Generating UMAP cluster plots...")
  tryCatch({
    pdf(paste0(output_prefix, "umap_clusters.pdf"), width = 10, height = 8)
    print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
            ggtitle("UMAP Clustering"))
    
    # 如果数据中有组织或条件信息，也绘制按组织/条件着色的UMAP
    if("Tissue" %in% colnames(seurat_obj@meta.data)) {
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "Tissue") + 
              ggtitle("UMAP by Tissue"))
    }
    
    if("condition" %in% colnames(seurat_obj@meta.data)) {
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "condition") + 
              ggtitle("UMAP by Condition"))
    }
    
    if("tissue" %in% colnames(seurat_obj@meta.data)) {
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "tissue") + 
              ggtitle("UMAP by Tissue"))
    }
    dev.off()
  }, error = function(e) {
    message(paste0("Error generating UMAP plots: ", e$message))
  })
  
  # Marker基因分析 - 使用全部细胞
  if(!is.null(markers)) {
    # 过滤掉不在数据中的marker基因
    valid_markers <- markers[markers %in% rownames(seurat_obj)]
    if(length(valid_markers) > 0) {
      message(paste("Analyzing", length(valid_markers), "marker genes..."))
      
      # 生成Feature plots
      tryCatch({
        feature_plot_width <- min(10, 5 * ceiling(sqrt(length(valid_markers))))
        feature_plot_height <- min(8, 4 * ceiling(length(valid_markers) / ceiling(sqrt(length(valid_markers)))))
        
        pdf(paste0(output_prefix, "umap_FeaturePlot.pdf"), width = feature_plot_width, height = feature_plot_height)
        # 每页最多显示6个基因
        for(i in seq(1, length(valid_markers), by = 6)) {
          end_idx <- min(i + 5, length(valid_markers))
          if(i <= end_idx) {  # 防止索引错误
            current_markers <- valid_markers[i:end_idx]
            p <- FeaturePlot(seurat_obj, 
                             features = current_markers,
                             raster = TRUE,
                             ncol = min(3, length(current_markers)))
            print(p)
          }
        }
        dev.off()
      }, error = function(e) {
        message(paste0("Error generating feature plots: ", e$message))
      })
      
      # DotPlot - 使用全部细胞
      tryCatch({
        dot_plot_width <- max(10, length(valid_markers) * 0.3)
        pdf(paste0(output_prefix, "markers_dotplot.pdf"), width = dot_plot_width, height = 10)
        print(DotPlot(seurat_obj, 
                      features = valid_markers, 
                      group.by = "seurat_clusters",
                      cols = dot_plot_cols,
                      dot.scale = 8) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        dev.off()
      }, error = function(e) {
        message(paste0("Error generating dot plot: ", e$message))
      })
      
      # 计算每个cluster中每个marker基因的表达统计 - 使用全部细胞进行计算
      message("Calculating marker gene expression statistics...")
      tryCatch({
        results <- data.frame()
        
        for(cluster in unique(Idents(seurat_obj))) {
          cells_in_cluster <- WhichCells(seurat_obj, idents = cluster)
          
          for(gene in valid_markers) {
            expr_data <- GetAssayData(seurat_obj, layer = "data")[gene, cells_in_cluster]
            mean_expr <- mean(expr_data)
            mean_expr_nonzero <- if(sum(expr_data > 0) > 0) mean(expr_data[expr_data > 0]) else 0
            
            pct_expr <- sum(expr_data > 0) / length(expr_data) * 100
            
            # 使用rbind而不是列表提高效率
            results <- rbind(results, data.frame(
              Cluster = cluster,
              Gene = gene,
              Mean_Expression = round(mean_expr, 3),
              Mean_Expression_NonZero = round(mean_expr_nonzero, 3),
              Percent_Expressing = round(pct_expr, 2)
            ))
          }
        }
        
        write.csv(results, paste0(output_prefix, "cluster_marker_stats.csv"), row.names = FALSE)
      }, error = function(e) {
        message(paste0("Error calculating marker statistics: ", e$message))
      })
    } else {
      warning("None of the provided markers were found in the dataset.")
    }
  }
  
  # 找marker基因 - 使用全部细胞进行计算
  cluster_markers_file <- paste0(output_prefix, "cluster_markers.csv")
  if(!file.exists(cluster_markers_file)) {
    message("Finding cluster-specific marker genes...")
    tryCatch({
      # 使用do.call允许传递自定义参数列表
      markers_args <- c(list(object = seurat_obj), find_markers_params)
      markers <- do.call(FindAllMarkers, markers_args)
      
      write.csv(markers, cluster_markers_file, row.names = FALSE)
    }, error = function(e) {
      message(paste0("Error finding markers: ", e$message))
      # 创建空的结果文件避免重复尝试
      write.csv(data.frame(), cluster_markers_file, row.names = FALSE)
      return(seurat_obj)
    })
  } else {
    message(paste0("Loading existing marker genes from: ", cluster_markers_file))
    markers <- read.csv(cluster_markers_file)
  }
  
  # 生成热图 - 仅在热图可视化时进行抽样
  if(file.exists(cluster_markers_file) && file.size(cluster_markers_file) > 0) {
    if(!exists("markers")) markers <- read.csv(cluster_markers_file)
    
    if(nrow(markers) > 0) {
      message("Generating marker gene heatmaps...")
      
      # 对每个群集取前5个marker基因(减少基因数量)
      tryCatch({
        top_markers <- markers %>% 
          group_by(cluster) %>% 
          top_n(5, wt = avg_log2FC)  # 减少为每个cluster的前5个
        
        marker_genes <- unique(top_markers$gene)
        
        if(length(marker_genes) > 0) {
          # 仅为热图创建抽样对象
          cell_count <- ncol(seurat_obj)
          
          if(cell_count > max_cells_for_heatmap) {
            message(paste0("Large dataset detected (", cell_count, " cells). Creating sampled object with ", 
                           max_cells_for_heatmap, " cells for heatmap visualization only."))
            
            # 按聚类进行分层抽样
            set.seed(42) # 设置随机种子确保可重复性
            
            # 获取每个聚类的细胞
            cells_by_cluster <- split(colnames(seurat_obj), seurat_obj$seurat_clusters)
            
            # 计算每个聚类应抽取的细胞数量（按比例）
            cells_per_cluster <- lapply(cells_by_cluster, function(cells) {
              n_sample <- ceiling(length(cells) / cell_count * max_cells_for_heatmap)
              if(length(cells) <= n_sample) return(cells)
              return(sample(cells, n_sample))
            })
            
            # 合并所有抽样细胞
            sampled_cells <- unlist(cells_per_cluster)
            
            # 创建子集用于热图
            heatmap_obj <- subset(seurat_obj, cells = sampled_cells)
            message(paste0("Created temporary object with ", ncol(heatmap_obj), " cells for heatmap visualization"))
          } else {
            heatmap_obj <- seurat_obj
          }
          
          # 仅缩放热图所需的基因
          heatmap_obj <- ScaleData(heatmap_obj, features = marker_genes)
          
          # 减小批次大小，每批次最多15个基因
          batch_size <- 15  # 热图中一次显示的基因数
          marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))
          
          for(i in seq_along(marker_batches)) {
            batch_genes <- marker_batches[[i]]
            if(length(batch_genes) > 0) {
              # 设置适当的高度
              pdf_height <- max(8, length(batch_genes) * 0.2)
              
              # 使用tryCatch捕获错误
              tryCatch({
                message(paste0("Creating heatmap batch ", i, " with ", length(batch_genes), " genes"))
                
                heatmap_file <- paste0(output_prefix, "marker_heatmap_batch_", i, ".pdf")
                pdf(heatmap_file, width = 10, height = pdf_height)
                
                # 使用raster=TRUE提高性能，减小点大小
                p <- DoHeatmap(heatmap_obj, 
                               features = batch_genes, 
                               group.by = "seurat_clusters",
                               size = 3,       
                               raster = TRUE)  # 使用光栅化提高性能
                
                print(p)
                dev.off()
                
                # 清理内存
                rm(p)
                gc()
                
              }, error = function(e) {
                message(paste0("Error generating heatmap batch ", i, ": ", e$message))
                # 如果文件已经创建但不完整，删除它
                if(file.exists(paste0(output_prefix, "marker_heatmap_batch_", i, ".pdf"))) {
                  file.remove(paste0(output_prefix, "marker_heatmap_batch_", i, ".pdf"))
                }
              })
            }
          }
          
          # 清理临时对象释放内存
          if(exists("heatmap_obj") && !identical(heatmap_obj, seurat_obj)) {
            rm(heatmap_obj)
            gc()
          }
        }
      }, error = function(e) {
        message(paste0("Error in heatmap preparation: ", e$message))
      })
    }
  }
  
  return(seurat_obj)
}

# 1. 读取所有RDS文件
batch_files <- list.files("E:\\R\\0216\\rds", pattern = "\\.rds$", full.names = TRUE)
if(length(batch_files) == 0) {
  stop("No RDS files found!")
}

# 2. 读取数据集
merged_list <- lapply(batch_files, function(file) {
  cat(sprintf("Processing: %s\n", basename(file)))
  current_batch <- readRDS(file)
  DietSeurat(current_batch, dimreducs = NULL, graphs = NULL)
})

# 3. 获取所有基因
all_genes <- unique(unlist(lapply(merged_list, rownames)))
cat("Total unique genes:", length(all_genes), "\n")

# 4. 对每个样本进行基因填充
filled_list <- list()
for(i in seq_along(merged_list)) {
  obj <- merged_list[[i]]
  missing_genes <- setdiff(all_genes, rownames(obj))
  
  if(length(missing_genes) > 0) {
    cat(sprintf("Batch %d: Adding %d missing genes\n", i, length(missing_genes)))
    
    # 创建填充矩阵
    fill_mat <- Matrix::Matrix(0, 
                               nrow = length(missing_genes), 
                               ncol = ncol(obj), 
                               sparse = TRUE)
    rownames(fill_mat) <- missing_genes
    colnames(fill_mat) <- colnames(obj)
    
    # 添加到原矩阵
    new_counts <- rbind(GetAssayData(obj, layer = "counts"), fill_mat)
    new_counts <- new_counts[all_genes,]  # 确保基因顺序一致
    obj[["RNA"]] <- CreateAssayObject(counts = new_counts)
  }
  filled_list[[i]] <- obj
}

# 5. 处理每个批次
processed_list <- list()
for(i in seq_along(filled_list)) {
  cat(sprintf("\nProcessing batch %d/%d\n", i, length(filled_list)))
  
  obj <- filled_list[[i]]
  cat(sprintf("Cells: %d, Features: %d\n", ncol(obj), nrow(obj)))
  
  # 标准化
  cat("Performing normalization...\n")
  obj <- NormalizeData(obj)
  
  # 找变异基因
  cat("Finding variable features...\n")
  obj <- FindVariableFeatures(obj, 
                              selection.method = "vst",
                              nfeatures = 3000)
  
  processed_list[[i]] <- obj
  cat("Batch processing completed\n")
}

# 6. 合并所有批次
cat("\nMerging all batches...\n")
seurat_obj_main <- merge(x = processed_list[[1]], 
                         y = processed_list[-1],
                         add.cell.ids = paste0("Batch", seq_along(processed_list)))
seurat_obj_main <- subset(seurat_obj_main, subset = tissue_sampling_method !='brush'
                          & study!='Jain_Misharin_2021')
# 8. 清理中间对象节省内存
rm(merged_list, filled_list, processed_list)
gc()

# 获取所有unique的sample
samples <- unique(seurat_obj_main$sample)
samples <- samples[!is.na(samples)]  # 移除NA值

# 为每个sample创建一个列表来存储结果
sample_results <- list()

# 循环处理每个sample
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
  
  # 数据准备
  current_obj <- NormalizeData(current_obj, verbose = TRUE)
  current_obj <- FindVariableFeatures(current_obj, verbose = TRUE)
  current_obj <- ScaleData(current_obj, verbose = TRUE)
  current_obj <- RunPCA(current_obj, npcs = 30, verbose = TRUE)
  
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
      # 统计doublet数量
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

# 将所有结果整合回主对象
seurat_obj_main$doublet_status <- "Unknown"  # 初始化列
for(sample_name in names(sample_results)) {
  cells_in_sample <- WhichCells(seurat_obj_main, expression = sample == sample_name)
  seurat_obj_main$doublet_status[cells_in_sample] <- sample_results[[sample_name]]
}

# 打印总体统计
cat("\nOverall doublet detection results:\n")
print(table(seurat_obj_main$doublet_status))

cells_before <- ncol(seurat_obj_main)
seurat_obj_main <- subset(seurat_obj_main, subset = doublet_status != "Doublet")
cells_after <- ncol(seurat_obj_main)
cat(sprintf("Removed %d doublets (%.1f%%)\n", 
            cells_before - cells_after, 
            (cells_before - cells_after)/cells_before * 100))

#' 优化版本：基于标记基因的双胞体检测函数
#' 
#' 基于细胞类型特异性标记基因的共表达模式来检测潜在双胞体
#' 按照样本分批处理，避免跨样本双胞体假阳性
#' 
#' @param seurat_obj Seurat对象
#' @param marker_genes 各细胞类型的标记基因列表
#' @param sample_col 样本信息列名，默认为"sample"
#' @param min_gene_threshold 共表达检测中的最小基因数阈值，默认为3
#' @param expression_threshold 基因表达阈值，默认为0
#' @param min_cells_per_sample 每个样本的最小细胞数阈值，低于此值的样本将被跳过，默认为50
#' @param adaptive_threshold 是否使用自适应阈值，默认为FALSE
#' @param parallel 是否使用并行计算，默认为FALSE
#' @param cores 并行计算核心数，默认为NULL（自动检测）
#' @param seed 随机数种子，确保结果可重复性
#' @param verbose 是否显示详细日志，默认为TRUE
#' @return 包含分析结果和更新后对象的列表
#' @examples
#' \dontrun{
#' # 基本用法
#' results <- detect_marker_doublets_by_sample(seurat_obj, marker_genes)
#' 
#' # 使用自适应阈值和并行计算
#' results <- detect_marker_doublets_by_sample(
#'   seurat_obj, 
#'   marker_genes,
#'   expr_threshold_method = "adaptive",
#'   parallel = TRUE,
#'   cores = 4
#' )
#' }

#' 整合样本结果
#'
#' @param sample_results 样本处理结果列表
#' @param seurat_obj Seurat对象
#' @param verbose 是否输出详细信息
#' @return 整合后的结果
integrate_sample_results <- function(sample_results, seurat_obj, verbose = TRUE) {
  # 初始化结果
  all_cells <- colnames(seurat_obj)
  is_doublet <- rep(FALSE, length(all_cells))
  names(is_doublet) <- all_cells
  doublet_types <- rep("None", length(all_cells))
  names(doublet_types) <- all_cells
  cell_states <- rep("Normal", length(all_cells))
  names(cell_states) <- all_cells
  
  # 整合各样本的结果
  total_doublets <- 0
  sample_stats <- list()
  
  for (result in sample_results) {
    # 更新双胞体标记
    is_doublet[names(result$is_doublet)] <- result$is_doublet
    doublet_types[names(result$doublet_types)] <- result$doublet_types
    cell_states[names(result$cell_states)] <- result$cell_states
    
    # 累计双胞体数量
    total_doublets <- total_doublets + result$doublet_cells
    
    # 样本统计
    sample_stats[[result$sample]] <- list(
      sample = result$sample,
      total_cells = result$total_cells,
      doublet_cells = result$doublet_cells,
      doublet_percent = result$doublet_percent
    )
  }
  
  # 转换样本统计为数据框
  sample_stats_df <- do.call(rbind, lapply(sample_stats, function(x) {
    data.frame(
      sample = x$sample,
      total_cells = x$total_cells,
      doublet_cells = x$doublet_cells,
      doublet_percent = x$doublet_percent
    )
  }))
  
  # 计算总体双胞体比例
  doublet_percentage <- total_doublets / length(all_cells) * 100
  
  if (verbose) {
    cat(paste("\n所有样本中总共检测到", total_doublets, "个标记基因双胞体 (", 
              round(doublet_percentage, 2), "%)\n"))
  }
  
  return(list(
    is_doublet = is_doublet,
    doublet_types = doublet_types,
    cell_states = cell_states,
    sample_stats = sample_stats_df,
    total_doublets = total_doublets,
    doublet_percentage = doublet_percentage
  ))
}

#' 计算聚类级别的统计信息
#'
#' @param seurat_obj Seurat对象
#' @param verbose 是否输出详细信息
#' @return 聚类统计数据框
calculate_cluster_stats <- function(seurat_obj, verbose = TRUE) {
  # 直接使用factor创建交叉表，避免列名问题
  doublet_factor <- factor(seurat_obj$is_marker_doublet, 
                           levels = c(FALSE, TRUE), 
                           labels = c("Singlet", "Doublet"))
  
  cluster_doublet_table <- table(seurat_obj$seurat_clusters, doublet_factor)
  
  # 计算每个聚类的双胞体比例
  cluster_totals <- rowSums(cluster_doublet_table)
  cluster_doublet_count <- cluster_doublet_table[, "Doublet"]
  cluster_doublet_percent <- (cluster_doublet_count / cluster_totals) * 100
  
  # 保存聚类双胞体分布
  cluster_stats <- data.frame(
    cluster = names(cluster_doublet_count),
    doublet_count = cluster_doublet_count,
    total_cells = cluster_totals,
    doublet_percent = cluster_doublet_percent
  )
  
  # 识别高双胞体比例聚类（通常>25%可能是双胞体富集）
  high_doublet_clusters <- which(cluster_doublet_percent > 25)
  if (length(high_doublet_clusters) > 0 && verbose) {
    cat("\n警告：以下聚类的双胞体比例异常高:\n")
    for (cl in high_doublet_clusters) {
      cat(paste("  聚类", names(cluster_doublet_percent)[cl], ": ", 
                round(cluster_doublet_percent[cl], 2), "%\n"))
    }
  }
  
  return(cluster_stats)
}

#' 创建可视化
#'
#' @param seurat_obj Seurat对象
#' @return 可视化图形列表
create_visualizations <- function(seurat_obj) {
  require(Seurat)
  require(ggplot2)
  
  # 检查是否有UMAP降维结果
  if (!"umap" %in% names(seurat_obj@reductions)) {
    warning("没有找到UMAP降维结果，无法创建可视化")
    return(NULL)
  }
  
  # 创建双胞体分布图
  p1 <- Seurat::DimPlot(seurat_obj, 
                        reduction = "umap", 
                        group.by = "is_marker_doublet", 
                        cols = c("FALSE" = "grey", "TRUE" = "red")) + 
    ggtitle("标记基因双胞体分布")
  
  # 双胞体得分分布图
  if ("doublet_score" %in% colnames(seurat_obj@meta.data)) {
    p2 <- Seurat::FeaturePlot(seurat_obj, 
                              features = "doublet_score", 
                              cols = c("lightgrey", "red")) + 
      ggtitle("双胞体得分分布")
  } else {
    p2 <- NULL
  }
  
  # 双胞体类型分布图
  p3 <- Seurat::DimPlot(seurat_obj, 
                        reduction = "umap", 
                        group.by = "marker_doublet_type", 
                        cols = c("None" = "grey")) + 
    ggtitle("双胞体类型分布")
  
  # 特殊细胞状态分布图
  if (any(seurat_obj$cell_state != "Normal")) {
    p4 <- Seurat::DimPlot(seurat_obj, 
                          reduction = "umap", 
                          group.by = "cell_state") + 
      ggtitle("细胞状态分布")
  } else {
    p4 <- NULL
  }
  
  # 返回可视化列表
  visualizations <- list(
    doublet_distribution = p1,
    doublet_score = p2,
    doublet_type = p3,
    cell_state = p4
  )
  
  # 过滤NULL值
  return(visualizations[!sapply(visualizations, is.null)])
}
detect_marker_doublets_by_sample <- function(seurat_obj, 
                                             marker_genes,
                                             sample_col = "sample",
                                             min_gene_threshold = 3,
                                             expression_threshold = 0,
                                             min_cells_per_sample = 50,
                                             adaptive_threshold = FALSE,
                                             parallel = FALSE,
                                             cores = NULL,
                                             seed = 42,
                                             verbose = TRUE) {
  # 设置随机数种子以确保可重复性
  set.seed(seed)
  
  # 初始化结果和计时
  results <- list()
  start_time <- Sys.time()
  
  # 基础输入验证
  validate_input(seurat_obj, marker_genes, sample_col)
  
  # 处理样本列缺失的情况
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    warning(paste0("样本列'", sample_col, "'不存在，创建单一样本列"))
    seurat_obj@meta.data[[sample_col]] <- "sample1"
  }
  
  # 自适应阈值计算
  if (adaptive_threshold) {
    thresholds <- calculate_adaptive_thresholds(seurat_obj, marker_genes)
    expression_threshold <- thresholds$expression
    min_gene_threshold <- thresholds$min_gene
    if (verbose) {
      cat(paste("使用自适应阈值: 表达阈值 =", round(expression_threshold, 4), 
                "，最小基因阈值 =", min_gene_threshold, "\n"))
    }
  }
  
  # 标记基因准备与验证
  valid_markers <- prepare_markers(marker_genes, seurat_obj, verbose)
  marker_genes <- valid_markers$markers
  valid_cell_types <- valid_markers$valid_types
  
  # 定义双胞体模式
  doublet_patterns <- define_doublet_patterns(
    valid_cell_types, marker_genes, min_gene_threshold, verbose)
  
  # 预计算细胞类型得分（一次性计算所有细胞）
  cell_type_scores <- calculate_cell_type_scores(
    seurat_obj, marker_genes, valid_cell_types, expression_threshold, verbose)
  
  # 添加细胞类型得分到元数据
  seurat_obj <- add_cell_type_scores(seurat_obj, cell_type_scores, verbose)
  
  # 获取样本列表
  samples <- unique(seurat_obj@meta.data[[sample_col]])
  if (verbose) {
    cat(paste("发现", length(samples), "个样本，开始处理...\n"))
  }
  
  # 设置并行处理
  if (parallel) {
    cl <- setup_parallel(cores, samples)
  }
  
  # 按样本处理（并行或串行）
  sample_results <- process_samples(
    seurat_obj, samples, sample_col, doublet_patterns, 
    expression_threshold, min_cells_per_sample, parallel, cl, verbose)
  
  # 关闭并行集群（如果使用）
  if (parallel && !is.null(cl)) {
    parallel::stopCluster(cl)
  }
  
  # 整合结果并更新对象
  integrated_results <- integrate_sample_results(
    sample_results, seurat_obj, verbose)
  
  # 更新Seurat对象的元数据
  seurat_obj$marker_doublet_type <- integrated_results$doublet_types
  seurat_obj$is_marker_doublet <- integrated_results$is_doublet
  seurat_obj$cell_state <- integrated_results$cell_states
  
  # 创建过滤后的对象
  seurat_obj_filtered <- subset(seurat_obj, subset = is_marker_doublet == FALSE)
  
  # 计算聚类级别的统计信息（如果有聚类）
  cluster_stats <- NULL
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_stats <- calculate_cluster_stats(seurat_obj, verbose)
  }
  
  # 计算执行时间
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "mins")
  
  # 汇总结果
  results$summary <- list(
    original_cells = ncol(seurat_obj),
    doublet_cells = sum(seurat_obj$is_marker_doublet, na.rm = TRUE),
    doublet_percent = sum(seurat_obj$is_marker_doublet, na.rm = TRUE) / ncol(seurat_obj) * 100,
    remaining_cells = ncol(seurat_obj_filtered),
    special_state_cells = sum(seurat_obj$cell_state != "Normal", na.rm = TRUE),
    execution_time_mins = as.numeric(exec_time)
  )
  
  # 存储处理后的对象和详细结果
  results$seurat_obj <- seurat_obj
  results$seurat_obj_filtered <- seurat_obj_filtered
  results$marker_genes <- marker_genes
  results$doublet_patterns <- doublet_patterns
  results$sample_stats <- integrated_results$sample_stats
  results$cluster_stats <- cluster_stats
  results$parameters <- list(
    min_gene_threshold = min_gene_threshold,
    expression_threshold = expression_threshold,
    adaptive_threshold = adaptive_threshold
  )
  
  # 添加可视化函数
  results$visualize <- function() {
    return(create_visualizations(seurat_obj))
  }
  
  if (verbose) {
    cat(paste("\n分析完成。原始细胞:", ncol(seurat_obj), 
              "，标记基因双胞体:", results$summary$doublet_cells,
              "，过滤后细胞:", ncol(seurat_obj_filtered),
              "\n执行时间:", round(exec_time, 2), "分钟\n"))
  }
  
  return(results)
}

#' 验证输入参数
#'
#' @param seurat_obj Seurat对象
#' @param marker_genes 标记基因列表
#' @param sample_col 样本列名
#' @return 无返回值，如有问题则停止执行
validate_input <- function(seurat_obj, marker_genes, sample_col) {
  # 验证Seurat对象
  if (!inherits(seurat_obj, "Seurat")) {
    stop("输入对象必须是Seurat对象")
  }
  
  # 验证标记基因
  if (is.null(marker_genes) || length(marker_genes) == 0) {
    stop("必须提供标记基因列表")
  }
  
  # 验证标记基因格式
  if (!is.list(marker_genes)) {
    stop("标记基因必须以列表形式提供，每个元素对应一种细胞类型")
  }
  
  # 更多验证逻辑可以在这里添加
  return(TRUE)
}

#' 计算自适应阈值
#'
#' @param seurat_obj Seurat对象
#' @param marker_genes 标记基因列表
#' @return 包含自适应阈值的列表
calculate_adaptive_thresholds <- function(seurat_obj, marker_genes) {
  # 获取表达矩阵
  expr_matrix <- Seurat::GetAssayData(seurat_obj, slot = "data")
  
  # 提取所有标记基因
  all_markers <- unique(unlist(marker_genes))
  all_markers <- all_markers[all_markers %in% rownames(expr_matrix)]
  
  # 计算非零表达值的分布
  marker_expr <- expr_matrix[all_markers, ]
  non_zero <- marker_expr[marker_expr > 0]
  
  # 使用10%分位数作为表达阈值
  expr_threshold <- quantile(non_zero, 0.1, na.rm = TRUE)
  
  # 计算标记基因列表长度的分布
  marker_lengths <- sapply(marker_genes, length)
  
  # 使用中位数的30%作为最小基因阈值，但至少为2
  min_gene_threshold <- max(2, round(median(marker_lengths) * 0.3))
  
  return(list(
    expression = expr_threshold,
    min_gene = min_gene_threshold
  ))
}

#' 准备和验证标记基因
#'
#' @param marker_genes 标记基因列表
#' @param seurat_obj Seurat对象
#' @param verbose 是否输出详细信息
#' @return 包含验证后标记基因和有效细胞类型的列表
prepare_markers <- function(marker_genes, seurat_obj, verbose = TRUE) {
  if (verbose) {
    cat("准备标记基因...\n")
  }
  
  # 确保标记基因存在于数据集中
  for (cell_type in names(marker_genes)) {
    genes_in_data <- marker_genes[[cell_type]][marker_genes[[cell_type]] %in% rownames(seurat_obj)]
    marker_genes[[cell_type]] <- genes_in_data
    
    # 检查是否有足够的标记基因
    if (length(genes_in_data) == 0) {
      warning(paste("警告:", cell_type, "类型没有可用的标记基因"))
    } else if (verbose) {
      cat(paste(" -", cell_type, "标记基因:", length(genes_in_data), "个\n"))
    }
  }
  
  # 保存有效的细胞类型
  valid_cell_types <- names(marker_genes)[sapply(marker_genes, length) > 0]
  if (length(valid_cell_types) < 2) {
    stop("至少需要两种有效的细胞类型才能检测双胞体")
  }
  
  return(list(
    markers = marker_genes,
    valid_types = valid_cell_types
  ))
}

#' 定义细胞谱系关系
#'
#' @return 细胞谱系关系列表
define_cell_lineages <- function() {
  # 定义细胞谱系关系，用于确定生物学兼容性
  cell_lineages <- list(
    Immune = c("T_cell", "B_cell", "NK_cell", "Myeloid", "Macrophage", "Dendritic", "Neutrophil"),
    Structural = c("Fibroblast", "SMC", "Myofibroblast", "Endothelial", "Pericyte"),
    Epithelial = c("Epithelial", "Basal", "Luminal"),
    Neural = c("Neuron", "Oligodendrocyte", "Astrocyte", "Microglia")
  )
  
  return(cell_lineages)
}

#' 定义双胞体模式
#'
#' @param valid_cell_types 有效的细胞类型列表
#' @param marker_genes 标记基因列表
#' @param min_gene_threshold 最小基因阈值
#' @param verbose 是否输出详细信息
#' @return 双胞体模式列表
define_doublet_patterns <- function(valid_cell_types, marker_genes, 
                                    min_gene_threshold, verbose = TRUE) {
  if (verbose) {
    cat("定义生物学不兼容的细胞类型组合...\n")
  }
  
  # 定义细胞谱系关系
  cell_lineages <- define_cell_lineages()
  
  # 自动生成所有可能的双胞体组合
  doublet_patterns <- list()
  cell_type_pairs <- combn(valid_cell_types, 2)
  
  for (i in 1:ncol(cell_type_pairs)) {
    type1 <- cell_type_pairs[1, i]
    type2 <- cell_type_pairs[2, i]
    
    # 检查是否属于同一谱系
    same_lineage <- FALSE
    for (lineage_name in names(cell_lineages)) {
      lineage <- cell_lineages[[lineage_name]]
      if ((type1 %in% lineage && type2 %in% lineage)) {
        same_lineage <- TRUE
        break
      }
    }
    
    # 特定组合的生物学兼容性检查
    if ((type1 == "T_cell" && type2 == "B_cell") || 
        (type1 == "B_cell" && type2 == "T_cell") ||
        (type1 == "T_cell" && type2 == "Myeloid") || 
        (type1 == "Myeloid" && type2 == "T_cell") ||
        (type1 == "B_cell" && type2 == "Myeloid") || 
        (type1 == "Myeloid" && type2 == "B_cell")) {
      # 免疫细胞之间共存是正常的，跳过这些组合
      same_lineage <- TRUE
    }
    
    # 跳过同一谱系的组合（可能是生物学兼容的）
    if (same_lineage) next
    
    pattern_name <- paste(type1, type2, sep = "_")
    doublet_patterns[[pattern_name]] <- list(
      group1 = marker_genes[[type1]],
      group2 = marker_genes[[type2]],
      min_group1 = min(min_gene_threshold, length(marker_genes[[type1]])),
      min_group2 = min(min_gene_threshold, length(marker_genes[[type2]]))
    )
  }
  
  if (verbose) {
    cat(paste(" - 定义了", length(doublet_patterns), "种生物学不兼容的细胞类型组合\n"))
  }
  
  return(doublet_patterns)
}

#' 计算细胞类型得分
#'
#' @param seurat_obj Seurat对象
#' @param marker_genes 标记基因列表
#' @param valid_cell_types 有效的细胞类型
#' @param expression_threshold 表达阈值
#' @param verbose 是否输出详细信息
#' @return 细胞类型得分数据框
calculate_cell_type_scores <- function(seurat_obj, marker_genes, valid_cell_types, 
                                       expression_threshold, verbose = TRUE) {
  if (verbose) {
    cat("计算所有细胞的标记基因表达特征...\n")
  }
  
  # 获取表达矩阵
  expr_matrix <- Seurat::GetAssayData(seurat_obj, slot = "data")
  
  # 初始化结果数据框
  cell_type_scores <- data.frame(cell_id = colnames(seurat_obj))
  
  # 预先计算每种细胞类型的表达矩阵
  expr_by_type <- list()
  for (cell_type in valid_cell_types) {
    markers <- marker_genes[[cell_type]]
    if (length(markers) == 0) next
    
    # 提取这组标记基因的表达矩阵
    expr_by_type[[cell_type]] <- expr_matrix[markers, , drop = FALSE]
  }
  
  # 计算每种细胞类型的标记基因表达量和比例
  for (cell_type in valid_cell_types) {
    markers <- marker_genes[[cell_type]]
    if (length(markers) == 0) next
    
    # 使用预计算的表达矩阵
    type_expr <- expr_by_type[[cell_type]]
    
    # 使用稀疏矩阵优化的方式计算
    expressed_markers <- Matrix::colSums(type_expr > expression_threshold)
    expression_ratio <- expressed_markers / length(markers)
    
    # 计算表达强度得分 (平均表达值)
    expression_strength <- Matrix::colMeans(type_expr)
    
    # 添加到结果数据框
    cell_type_scores[[paste0(cell_type, "_count")]] <- expressed_markers
    cell_type_scores[[paste0(cell_type, "_ratio")]] <- expression_ratio
    cell_type_scores[[paste0(cell_type, "_score")]] <- expression_strength
  }
  
  return(cell_type_scores)
}

#' 将细胞类型得分添加到Seurat对象
#'
#' @param seurat_obj Seurat对象
#' @param cell_type_scores 细胞类型得分数据框
#' @param verbose 是否输出详细信息
#' @return 更新后的Seurat对象
add_cell_type_scores <- function(seurat_obj, cell_type_scores, verbose = TRUE) {
  # 将分数添加到元数据
  cell_type_scores_subset <- cell_type_scores[, -1, drop = FALSE]
  row.names(cell_type_scores_subset) <- cell_type_scores$cell_id
  
  # 确保行名匹配
  common_cells <- intersect(rownames(cell_type_scores_subset), colnames(seurat_obj))
  seurat_obj@meta.data <- cbind(
    seurat_obj@meta.data, 
    cell_type_scores_subset[common_cells, , drop = FALSE]
  )
  
  # 计算每个细胞的主导细胞类型
  ratio_cols <- grep("_ratio$", colnames(seurat_obj@meta.data), value = TRUE)
  
  if (length(ratio_cols) > 0) {
    # 获取每个细胞的主导细胞类型
    seurat_obj$dominant_cell_type <- apply(seurat_obj@meta.data[, ratio_cols, drop = FALSE], 1, function(x) {
      if (all(is.na(x))) return(NA)
      cell_type <- gsub("_ratio$", "", ratio_cols[which.max(x)])
      return(cell_type)
    })
    
    # 计算主导细胞类型的表达率
    seurat_obj$dominant_ratio <- apply(seurat_obj@meta.data[, ratio_cols, drop = FALSE], 1, function(x) {
      if (all(is.na(x))) return(NA)
      return(max(x, na.rm = TRUE))
    })
    
    # 计算次要细胞类型的表达率
    seurat_obj$secondary_ratio <- apply(seurat_obj@meta.data[, ratio_cols, drop = FALSE], 1, function(x) {
      if (length(x) < 2 || all(is.na(x))) return(NA)
      sorted_x <- sort(x, decreasing = TRUE)
      if (length(sorted_x) >= 2) return(sorted_x[2])
      return(0)
    })
    
    # 计算混杂得分 - 表达多个细胞类型标记的程度
    threshold <- 0.2  # 可调整的阈值
    seurat_obj$marker_mixing_score <- rowSums(seurat_obj@meta.data[, ratio_cols, drop = FALSE] >= threshold, na.rm = TRUE)
    
    # 计算双胞体得分
    seurat_obj$doublet_score <- seurat_obj$secondary_ratio / 
      (seurat_obj$dominant_ratio + 0.01)  # 添加小值避免除零
  }
  
  # 初始化新的元数据列
  seurat_obj$marker_doublet_type <- "None"
  seurat_obj$is_marker_doublet <- FALSE
  seurat_obj$cell_state <- "Normal"
  
  return(seurat_obj)
}

#' 设置并行处理环境
#'
#' @param cores 核心数
#' @param samples 样本列表
#' @return 并行集群对象
setup_parallel <- function(cores, samples) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("parallel包不可用，将使用串行处理")
    return(NULL)
  }
  
  # 设置核心数
  if (is.null(cores)) {
    cores <- min(parallel::detectCores() - 1, length(samples), 8)  # 默认最多8核
  }
  cores <- max(1, min(cores, parallel::detectCores() - 1))  # 安全检查
  
  # 创建集群
  cl <- parallel::makeCluster(cores)
  
  # 如果可用，设置进度条
  if (requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pboptions(use_lb = TRUE)
  }
  
  return(cl)
}

#' 处理样本（并行或串行）
#'
#' @param seurat_obj Seurat对象
#' @param samples 样本列表
#' @param sample_col 样本列名
#' @param doublet_patterns 双胞体模式
#' @param expression_threshold 表达阈值
#' @param min_cells_per_sample 每个样本的最小细胞数
#' @param parallel 是否使用并行处理
#' @param cl 并行集群对象
#' @param verbose 是否输出详细信息
#' @return 样本处理结果列表
# 处理样本（并行或串行）
process_samples <- function(seurat_obj, samples, sample_col, doublet_patterns,
                            expression_threshold, min_cells_per_sample, 
                            parallel = FALSE, cl = NULL, verbose = TRUE) {
  # 获取表达矩阵
  expr_matrix <- Seurat::GetAssayData(seurat_obj, slot = "data")
  
  # 预计算每种双胞体模式的标记基因表达矩阵
  pattern_expr <- list()
  for (pattern_name in names(doublet_patterns)) {
    pattern <- doublet_patterns[[pattern_name]]
    pattern_expr[[pattern_name]] <- list(
      group1 = expr_matrix[pattern$group1, , drop = FALSE],
      group2 = expr_matrix[pattern$group2, , drop = FALSE]
    )
  }
  
  process_one_sample <- function(current_sample) {
    if (verbose) {
      cat(paste("\n处理样本:", current_sample, "\n"))
    }
    
    # 获取当前样本的细胞索引
    sample_cells_idx <- which(seurat_obj@meta.data[[sample_col]] == current_sample)
    
    # 检查样本细胞数量
    if (length(sample_cells_idx) < min_cells_per_sample) {
      if (verbose) {
        cat(paste("  跳过样本", current_sample, "，细胞数量过少(", length(sample_cells_idx), ")\n"))
      }
      return(NULL)
    }
    
    # 获取当前样本的细胞ID
    sample_cells <- rownames(seurat_obj@meta.data)[sample_cells_idx]
    if (verbose) {
      cat(paste("  处理", length(sample_cells), "个细胞...\n"))
    }
    
    # 创建当前样本的双胞体结果存储
    sample_doublets <- list()
    is_doublet <- rep(FALSE, length(sample_cells))
    names(is_doublet) <- sample_cells
    doublet_types <- rep("None", length(sample_cells))
    names(doublet_types) <- sample_cells
    
    # 对每种双胞体模式进行检测
    for (pattern_name in names(doublet_patterns)) {
      pattern <- doublet_patterns[[pattern_name]]
      
      # 使用预计算的表达矩阵 - 仅处理当前样本的细胞
      group1_expr <- pattern_expr[[pattern_name]]$group1[, sample_cells, drop = FALSE]
      group2_expr <- pattern_expr[[pattern_name]]$group2[, sample_cells, drop = FALSE]
      
      # 找出表达足够多第一组标记基因的细胞
      cells_express_group1 <- colnames(group1_expr)[
        Matrix::colSums(group1_expr > expression_threshold) >= pattern$min_group1
      ]
      
      # 找出表达足够多第二组标记基因的细胞
      cells_express_group2 <- colnames(group2_expr)[
        Matrix::colSums(group2_expr > expression_threshold) >= pattern$min_group2
      ]
      
      # 找出同时满足两个条件的细胞
      potential_doublets <- intersect(cells_express_group1, cells_express_group2)
      
      # 保存结果
      sample_doublets[[pattern_name]] <- potential_doublets
      
      # 标记双胞体类型
      if (length(potential_doublets) > 0) {
        is_doublet[potential_doublets] <- TRUE
        doublet_types[potential_doublets] <- pattern_name
      }
      
      if (verbose) {
        cat(paste("    -", pattern_name, ":", length(potential_doublets), "个可能的双胞体\n"))
      }
    }
    
    # 考虑上皮-间质转化 (EMT)
    cell_states <- rep("Normal", length(sample_cells))
    names(cell_states) <- sample_cells
    
    if (all(c("Epithelial_ratio", "Fibroblast_ratio") %in% colnames(seurat_obj@meta.data))) {
      # 获取潜在EMT细胞
      meta_subset <- seurat_obj@meta.data[sample_cells, ]
      potential_emt_cells <- sample_cells[
        meta_subset$Epithelial_ratio >= 0.3 & 
          meta_subset$Fibroblast_ratio >= 0.2 &
          is_doublet == TRUE  # 已被标记为双胞体
      ]
      
      # 将这些潜在的EMT细胞标记为特殊类型而非双胞体
      if (length(potential_emt_cells) > 0) {
        is_doublet[potential_emt_cells] <- FALSE
        cell_states[potential_emt_cells] <- "Potential_EMT"
        
        if (verbose) {
          cat(paste("    - 识别了", length(potential_emt_cells), "个潜在EMT细胞\n"))
        }
      }
    }
    
    # 计算当前样本的双胞体数量
    sample_doublet_count <- sum(is_doublet, na.rm = TRUE)
    sample_doublet_percent <- (sample_doublet_count / length(sample_cells)) * 100
    
    if (verbose) {
      cat(paste("  样本", current_sample, "中检测到", sample_doublet_count, 
                "个双胞体 (", round(sample_doublet_percent, 2), "%)\n"))
    }
    
    # 返回结果
    return(list(
      sample = current_sample,
      total_cells = length(sample_cells),
      doublet_cells = sample_doublet_count,
      doublet_percent = sample_doublet_percent,
      is_doublet = is_doublet,
      doublet_types = doublet_types,
      cell_states = cell_states,
      doublets_by_pattern = sample_doublets
    ))
  }
  
  # 执行样本处理（并行或串行）
  if (parallel && !is.null(cl)) {
    # 修改这里: 导出正确的变量名 - 使用在函数内部定义的变量名
    parallel::clusterExport(cl, varlist = c("seurat_obj", "expr_matrix", 
                                            "doublet_patterns", "pattern_expr",
                                            "expression_threshold", "sample_col",
                                            "min_cells_per_sample", "verbose"))
    
    # 并行处理样本
    if (requireNamespace("pbapply", quietly = TRUE)) {
      sample_results <- pbapply::pblapply(samples, process_one_sample, cl = cl)
    } else {
      sample_results <- parallel::parLapply(cl, samples, process_one_sample)
    }
  } else {
    # 串行处理样本
    sample_results <- lapply(samples, process_one_sample)
  }
  
  # 过滤空结果
  sample_results <- sample_results[!sapply(sample_results, is.null)]
  
  return(sample_results)
}
  
# 最终优化的标记基因列表
marker_genes <- list(
  # 上皮细胞标记基因
  Epithelial = c('EPCAM', 'KRT8', 'FXYD3', 'ELF3', 'CLDN4'),
  
  # T细胞标记基因 - 移除CD4，使用更T细胞特异的标记
  T_cell = c('CD3E', 'CD3D', 'CD2', 'CD8A', 'TRAC'),
  
  # B细胞标记基因
  B_cell = c('MS4A1', 'CD79A', 'CD19', 'TNFRSF17'),
  
  # 髓系细胞标记基因
  Myeloid = c('CD14', 'FCER1G', 'CLEC7A', 'CD86', 'C1orf162'),
  
  # 内皮细胞标记基因
  Endothelial = c('PECAM1', 'VWF', 'CLDN5', 'PTPRB', 'ECSCR'),
  
  # 成纤维细胞标记基因
  Fibroblast = c('COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'),
  
  # 平滑肌细胞标记基因
  SMC = c('ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'DES')
)

message("Running marker-based doublet detection...")
doublet_results <- detect_marker_doublets_by_sample(
  seurat_obj = seurat_obj_main,     # 确保这里使用正确的对象名称
  marker_genes = marker_genes,      
  sample_col = "sample",            
  min_gene_threshold = 3,           
  expression_threshold = 0,         
  min_cells_per_sample = 50,
  adaptive_threshold = TRUE,        
  parallel = FALSE,                  
  cores = NULL,                     
  verbose = TRUE
)
# 查看双胞体检测结果
print(doublet_results$summary)
print(doublet_results$sample_stats)

# 可视化双胞体分布
p1 <- DimPlot(doublet_results$seurat_obj, 
              reduction = "umap", 
              group.by = "is_marker_doublet", 
              cols = c("FALSE" = "grey", "TRUE" = "red")) + 
  ggtitle("标记基因双胞体分布")

# 按样本着色可视化
p2 <- DimPlot(doublet_results$seurat_obj, 
              reduction = "umap", 
              group.by = "sample") + 
  ggtitle("样本分布")

# 共同显示
pdf("doublet_detection_plots.pdf", width = 14, height = 7)
print(p1 + p2)
dev.off()

# 可视化双胞体类型分布
pdf("doublet_types.pdf", width = 10, height = 8)
DimPlot(doublet_results$seurat_obj, 
        reduction = "umap", 
        group.by = "marker_doublet_type", 
        cols = c("None" = "grey")) + 
  ggtitle("双胞体类型分布")
dev.off()

# 更新对象为过滤后的版本
seurat_obj_main <- doublet_results$seurat_obj_filtered
message(paste0("Removed ", doublet_results$summary$doublet_cells, 
               " marker-based doublets (", 
               round(doublet_results$summary$doublet_percent, 2), "%)"))
# 7. 输出最终结果信息
cat(sprintf("\nFinal merged object:\n"))
cat(sprintf("Total cells: %d\n", ncol(seurat_obj_main)))
cat(sprintf("Total features: %d\n", nrow(seurat_obj_main)))

all_marker_genes <- c(
  # T细胞标记基因
  "CD3E",     # CD3ε
  "CD3D",     # CD3δ
  "CD3G",     # CD3γ
  "CD5",      # CD5
  "TRAC",     # T细胞受体α链恒定区
  "TRBC1",    # T细胞受体β链恒定区1
  "TRBC2",    # T细胞受体β链恒定区2
  "TRGC1",    # T细胞受体γ链恒定区1
  "TRGC2",    # T细胞受体γ链恒定区2
  "TRDC",     # T细胞受体δ链恒定区
  "ITGAE",    # CD103，组织驻留T细胞
  "CD69",     # 组织驻留T细胞标记
  
  # 髓系细胞标记基因
  "ITGAX",    # CD11c，树突状细胞
  "ITGAM",    # CD11b，多种髓系细胞
  "CD14",     # 单核/巨噬细胞
  "CD68",     # 巨噬细胞
  "SIGLEC1",  # CD169，某些巨噬细胞亚群
  "MRC1",     # CD206，M2型巨噬细胞
  "MARCO",    # 肺泡巨噬细胞特异性标记
  "FUT4",     # CD15，中性粒细胞
  "CEACAM8",  # CD66b，中性粒细胞
  "MPO",      # 髓过氧化物酶，中性粒细胞
  
  # 功能相关基因
  "IFNG",     # 干扰素γ，T细胞激活
  "IL10",     # 白细胞介素10
  "LYZ",      # 溶菌酶，髓系细胞
  "CST3",     # 胱抑素C，髓系细胞
  "S100A8",   # 钙结合蛋白
  "S100A9",   # 钙结合蛋白
  "LAMP1",    # 溶酶体相关膜蛋白1
  
  # 组织特异性和其他相关基因
  "SPP1",     # 骨桥蛋白
  "EPCAM",    # 上皮细胞标记
  "PECAM1",   # 内皮细胞标记
  "PTPRC",    # CD45，所有白细胞
  "SELL"      # CD62L，L-选择素
)

# 检查基因是否存在于数据集中
genes_present <- all_marker_genes[all_marker_genes %in% rownames(seurat_obj_main)]
genes_missing <- all_marker_genes[!all_marker_genes %in% rownames(seurat_obj_main)]

# 创建输出文件夹（如果不存在）
dir.create("marker_gene_plots", showWarnings = TRUE)

# 打印存在和缺失的基因
cat("存在的基因数量:", length(genes_present), "\n")
cat("缺失的基因数量:", length(genes_missing), "\n")

# 保存基因存在信息到文本文件
sink("marker_gene_plots/gene_status.txt")
cat("存在的基因数量:", length(genes_present), "\n")
cat("存在的基因:", paste(genes_present, collapse=", "), "\n\n")
cat("缺失的基因数量:", length(genes_missing), "\n")
if(length(genes_missing) > 0) {
  cat("缺失的基因:", paste(genes_missing, collapse=", "), "\n")
}
sink()

# 确保基因在scale.data中可用
# 由于很多基因不在scale.data中，我们需要先对这些基因进行标准化
genes_to_scale <- genes_present[!genes_present %in% rownames(seurat_obj_main[["RNA"]]@scale.data)]
if(length(genes_to_scale) > 0) {
  seurat_obj_main <- ScaleData(seurat_obj_main, features = genes_to_scale)
  cat("已对缺失的", length(genes_to_scale), "个基因进行标准化\n")
}

# 创建多页PDF文件
if(length(genes_present) > 0) {
  # 使用英文标题避免中文编码问题
  # 1. 热图 - 所有基因在聚类中的表达情况
  pdf("marker_gene_plots/heatmap_all_markers.pdf", width = 12, height = 20)
  print(
    DoHeatmap(seurat_obj_main, features = genes_present, group.by = "seurat_clusters") +
      ggtitle("Expression of Myeloid and T-cell Markers in Clusters")
  )
  dev.off()
  
  # 2. 关键基因组合的特征图 (一页一个基因)
  key_genes <- c("CD3E", "CD3D", "CD5", "ITGAE", "TRAC", # T细胞标记
                 "ITGAM", "ITGAX", "CD14", "CD68", "MPO") # 髓系标记
  key_genes_present <- key_genes[key_genes %in% genes_present]
  
  pdf("marker_gene_plots/feature_plots_key_genes.pdf", width = 10, height = 8)
  for(gene in key_genes_present) {
    # 使用英文标题，减少样本量
    p <- FeaturePlot(seurat_obj_main, 
                     features = gene, 
                     reduction = "umap",
                     raster = FALSE,
                     pt.size = 0.5,
                     cells = sample(colnames(seurat_obj_main), min(20000, ncol(seurat_obj_main)))) + 
      ggtitle(paste("Expression of", gene)) +
      theme(plot.title = element_text(size = 16, face = "bold"))
    print(p)
  }
  dev.off()
  
  # 3. 使用单独的FeaturePlot
  t_cell_markers <- c("CD3E", "CD3D", "CD5", "ITGAE", "TRAC")
  myeloid_markers <- c("ITGAM", "ITGAX", "CD14", "CD68", "MPO")
  
  t_present <- t_cell_markers[t_cell_markers %in% genes_present]
  myeloid_present <- myeloid_markers[myeloid_markers %in% genes_present]
  
  if(length(t_present) > 0 && length(myeloid_present) > 0) {
    # 尝试使用基本的R绘图功能而不是复杂的patchwork布局
    pdf("marker_gene_plots/paired_feature_plots.pdf", width = 12, height = 8)
    
    for(t_gene in t_present) {
      for(m_gene in myeloid_present) {
        # 绘制T细胞标记
        par(mfrow=c(1,2)) # 设置为1行2列的布局
        
        # 获取表达数据
        t_expr <- GetAssayData(seurat_obj_main, layer = "data")[t_gene, ]
        m_expr <- GetAssayData(seurat_obj_main, layer = "data")[m_gene, ]
        umap_coords <- Embeddings(seurat_obj_main, reduction = "umap")
        
        # 随机抽样20000个细胞
        set.seed(42)
        sample_idx <- sample(1:ncol(seurat_obj_main), min(20000, ncol(seurat_obj_main)))
        
        # T细胞标记图
        plot(umap_coords[sample_idx,1], umap_coords[sample_idx,2], 
             pch=16, cex=0.5, col=colorRampPalette(c("grey", "blue"))(10)[cut(t_expr[sample_idx], breaks=10)],
             main=paste("Expression of", t_gene), xlab="UMAP1", ylab="UMAP2")
        
        # 髓系标记图
        plot(umap_coords[sample_idx,1], umap_coords[sample_idx,2], 
             pch=16, cex=0.5, col=colorRampPalette(c("grey", "red"))(10)[cut(m_expr[sample_idx], breaks=10)],
             main=paste("Expression of", m_gene), xlab="UMAP1", ylab="UMAP2")
        
        # 重置绘图参数
        par(mfrow=c(1,1))
        
        # 添加散点图比较两个基因的表达
        plot(t_expr[sample_idx], m_expr[sample_idx], 
             pch=16, cex=0.5, col=adjustcolor("black", alpha.f=0.5),
             main=paste(t_gene, "vs", m_gene, "Expression"),
             xlab=t_gene, ylab=m_gene)
        abline(h=0, v=0, lty=2, col="grey")
      }
    }
    dev.off()
  }
  
  # 4. 小提琴图展示不同聚类中的标记基因表达
  pdf("marker_gene_plots/violin_plots.pdf", width = 12, height = 8)
  
  # 为避免一次绘制太多小提琴图，分批绘制
  t_batches <- split(t_present, ceiling(seq_along(t_present)/4))
  m_batches <- split(myeloid_present, ceiling(seq_along(myeloid_present)/4))
  
  for(t_batch in t_batches) {
    if(length(t_batch) > 0) {
      p <- VlnPlot(seurat_obj_main, features = t_batch, group.by = "seurat_clusters", pt.size = 0, ncol = 2) +
        plot_annotation(title = "T-cell Markers Expression")
      print(p)
    }
  }
  
  for(m_batch in m_batches) {
    if(length(m_batch) > 0) {
      p <- VlnPlot(seurat_obj_main, features = m_batch, group.by = "seurat_clusters", pt.size = 0, ncol = 2) +
        plot_annotation(title = "Myeloid Markers Expression")
      print(p)
    }
  }
  dev.off()
  
  # 5. 散点图比较潜在双表型标记对的表达 - 使用基本R绘图
  pdf("marker_gene_plots/scatter_plots.pdf", width = 10, height = 10)
  # 选择主要的几对T细胞-髓系标记组合
  key_pairs <- list(
    c("CD3E", "ITGAM"),  # CD3-CD11b
    c("CD3E", "ITGAX"),  # CD3-CD11c
    c("CD5", "CD14"),    # CD5-CD14
    c("CD3E", "CD68"),   # CD3-CD68
    c("ITGAE", "CD68")   # CD103-CD68
  )
  
  for(pair in key_pairs) {
    if(all(pair %in% genes_present)) {
      # 提取基因表达数据
      gene1 <- pair[1]
      gene2 <- pair[2]
      
      # 获取表达数据和聚类信息
      expr_data <- GetAssayData(seurat_obj_main, layer = "data")[c(gene1, gene2), ]
      clusters <- seurat_obj_main$seurat_clusters
      
      # 随机抽样
      set.seed(42)
      sample_idx <- sample(1:ncol(seurat_obj_main), min(10000, ncol(seurat_obj_main)))
      
      # 使用不同颜色表示不同聚类
      cluster_colors <- rainbow(length(unique(clusters)))
      names(cluster_colors) <- unique(clusters)
      
      # 绘制散点图
      plot(expr_data[1, sample_idx], expr_data[2, sample_idx], 
           pch=16, cex=0.5, 
           col=cluster_colors[clusters[sample_idx]],
           main=paste(gene1, "vs", gene2, "Expression"),
           xlab=gene1, ylab=gene2)
      
      # 添加图例
      legend("topright", legend=names(cluster_colors), 
             col=cluster_colors, pch=16, cex=0.8, title="Clusters")
      
      # 添加参考线
      abline(h=0, v=0, lty=2, col="grey")
    }
  }
  dev.off()
}

cat("所有图表已保存到 'marker_gene_plots' 文件夹\n")

# 设置T细胞标记和髓系标记
t_cell_markers <- c("CD3E", "CD5", "CD3D")
myeloid_markers <- c("ITGAM", "ITGAX", "CD14", "CD68")

# 创建输出PDF
pdf("marker_gene_plots/t_myeloid_blend_plots.pdf", width = 24, height = 10)

# 为每个T细胞标记与每个髓系标记创建blend FeaturePlot
for (t_marker in t_cell_markers) {
  if (t_marker %in% rownames(seurat_obj_main)) {
    for (m_marker in myeloid_markers) {
      if (m_marker %in% rownames(seurat_obj_main)) {
        # 检查两个标记是否都存在
        cat(paste("绘制", t_marker, "vs", m_marker, "blend plot\n"))
        
        # 使用blend参数创建共表达图
        p <- FeaturePlot(seurat_obj_main, 
                         features = c(t_marker, m_marker),
                         blend = TRUE)
        
        print(p + plot_annotation(title = paste(t_marker, "vs", m_marker, "Co-expression")))
      }
    }
  }
}

dev.off()

cat("所有blend共表达图表已保存到marker_gene_plots/t_myeloid_blend_plots.pdf文件\n")

# 设置T细胞标记和髓系标记
t_cell_markers <- c("CD3E", "CD5", "CD3D")
myeloid_markers <- c("ITGAM", "ITGAX", "CD14", "CD68")

# 定义函数：识别并提取双表型细胞
identify_dual_phenotype_cells <- function(seurat_obj, 
                                          t_marker, 
                                          m_marker, 
                                          threshold = 0.3,
                                          visualize = TRUE) {
  # 检查标记是否存在
  if(!(t_marker %in% rownames(seurat_obj) && m_marker %in% rownames(seurat_obj))) {
    cat(paste("标记", t_marker, "或", m_marker, "不在数据集中\n"))
    return(NULL)
  }
  
  # 获取表达数据
  expr_data <- GetAssayData(seurat_obj, layer = "data")[c(t_marker, m_marker), ]
  
  # 识别双表型细胞
  dual_cells <- colnames(expr_data)[expr_data[t_marker, ] > threshold & 
                                      expr_data[m_marker, ] > threshold]
  
  # 打印结果
  cat(paste("找到", length(dual_cells), "个同时表达", t_marker, 
            "和", m_marker, "的细胞\n"))
  
  # 如果需要可视化
  if(visualize && length(dual_cells) > 0) {
    # 在UMAP上展示这些细胞
    p1 <- DimPlot(seurat_obj, 
                  cells.highlight = dual_cells,
                  cols.highlight = "purple",
                  sizes.highlight = 1,
                  reduction = "umap") +
      ggtitle(paste(t_marker, "+", m_marker, "+细胞"))
    
    # 在表达散点图上展示
    all_cells <- sample(colnames(seurat_obj), min(10000, ncol(seurat_obj)))
    expr_data_sample <- FetchData(seurat_obj, 
                                  vars = c(t_marker, m_marker, "seurat_clusters"), 
                                  cells = all_cells)
    
    dual_cells_in_sample <- intersect(dual_cells, all_cells)
    expr_data_sample$cell_type <- "Other"
    expr_data_sample$cell_type[rownames(expr_data_sample) %in% dual_cells_in_sample] <- "Dual Phenotype"
    
    p2 <- ggplot(expr_data_sample, aes_string(x = t_marker, y = m_marker)) +
      geom_point(aes(color = cell_type), alpha = 0.7, size = 1) +
      scale_color_manual(values = c("Dual Phenotype" = "purple", "Other" = "grey")) +
      theme_classic() +
      labs(title = paste(t_marker, "vs", m_marker, "表达")) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "red")
    
    # 打印图表
    print(p1)
    print(p2)
  }
  
  # 返回双表型细胞的列表
  return(dual_cells)
}

# 处理所有标记组合并保存结果
pdf("marker_gene_plots/dual_phenotype_cells.pdf", width = 16, height = 8)

# 存储所有结果
all_dual_cells <- list()

# 处理每一对组合
for (t_marker in t_cell_markers) {
  for (m_marker in myeloid_markers) {
    # 定义组合名称
    combo_name <- paste(t_marker, m_marker, sep = "_")
    
    # 识别双表型细胞
    dual_cells <- identify_dual_phenotype_cells(seurat_obj_main, 
                                                t_marker, 
                                                m_marker, 
                                                threshold = 0.3,  # 可以调整阈值
                                                visualize = TRUE)
    
    # 保存结果
    all_dual_cells[[combo_name]] <- dual_cells
  }
}
dev.off()

# 保存所有双表型细胞的统计信息
sink("marker_gene_plots/dual_phenotype_stats.txt")
cat("双表型细胞统计信息:\n\n")

for (combo_name in names(all_dual_cells)) {
  cells <- all_dual_cells[[combo_name]]
  if (!is.null(cells)) {
    cat(paste(combo_name, "组合:", length(cells), "个细胞\n"))
    
    # 如果存在足够的双表型细胞，分析它们的聚类分布
    if (length(cells) >= 10) {
      cluster_counts <- table(seurat_obj_main$seurat_clusters[cells])
      cat("  聚类分布:\n")
      for (cluster in names(cluster_counts)) {
        cat(paste("    聚类", cluster, ":", cluster_counts[cluster], "个细胞\n"))
      }
    }
    cat("\n")
  }
}
sink()

# 创建一个新的元数据列标记双表型细胞
seurat_obj_main$dual_phenotype <- "None"

# 对每种组合的双表型细胞进行标记
for (combo_name in names(all_dual_cells)) {
  cells <- all_dual_cells[[combo_name]]
  if (!is.null(cells) && length(cells) > 0) {
    seurat_obj_main$dual_phenotype[colnames(seurat_obj_main) %in% cells] <- combo_name
  }
}

# 保存标记了双表型细胞的Seurat对象(可选)
# saveRDS(seurat_obj_main, "seurat_obj_with_dual_phenotype.rds")

# 创建一个合并的UMAP图，显示所有双表型细胞
pdf("marker_gene_plots/all_dual_phenotype_cells.pdf", width = 12, height = 10)

# 所有双表型细胞的合并列表
all_dual_cells_combined <- unique(unlist(all_dual_cells))

# 如果有双表型细胞，创建UMAP图
if (length(all_dual_cells_combined) > 0) {
  p <- DimPlot(seurat_obj_main,
               group.by = "dual_phenotype",
               cols = c("None" = "grey", 
                        setNames(rainbow(length(names(all_dual_cells))), 
                                 names(all_dual_cells)))) +
    ggtitle("所有双表型细胞在UMAP上的分布")
  print(p)
}

dev.off()
cat("原始细胞数:", ncol(seurat_obj_main), "\n")
seurat_obj_main <- subset(seurat_obj_main, subset = dual_phenotype == "None")
cat("剩余细胞数:", ncol(seurat_obj_main), "\n")
cat("双表型细胞分析完成，结果已保存到marker_gene_plots文件夹\n")

# saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")




regression_outliers_by_sample <- function(seurat_obj, outliers_threshold = 0.999) {
  # 创建输出目录 [^1]
  dir.create("qaqc_plots", recursive = TRUE, showWarnings = TRUE)
  
  # 存储所有要保留的细胞
  all_cells_keep <- vector()
  
  # 按样本处理
  for(sample_name in unique(seurat_obj$sample)) {
    # 获取当前样本数据
    current_obj <- subset(seurat_obj, subset = sample == sample_name)
    
    # 获取metrics数据
    nCount <- current_obj$nCount_RNA
    nFeature <- current_obj$nFeature_RNA
    
    # 对数转换
    log_counts <- log10(nCount + 1)
    log_features <- log10(nFeature + 1)
    
    # 线性回归
    fit <- lm(log_features ~ log_counts)
    
    # 计算预测区间
    pred <- predict(fit, 
                    data.frame(log_counts = log_counts), 
                    interval = "prediction", 
                    level = outliers_threshold)
    
    # 标记异常值
    outliers <- log_features < pred[,"lwr"] | log_features > pred[,"upr"]
    
    # 可视化
    pdf(file.path("qaqc_plots", 
                  paste0("regression_outliers_", sample_name, ".pdf")),
        width = 10, height = 8)
    
    # 绘制散点图
    plot(log_counts, log_features,
         pch = 16, cex = 0.6,
         col = ifelse(outliers, "red", "grey50"),
         xlab = "log10(nCount_RNA)",
         ylab = "log10(nFeature_RNA)",
         main = paste("Regression Outliers Detection -", sample_name))
    
    # 添加拟合线和预测区间
    lines(sort(log_counts), pred[order(log_counts),"lwr"], col = "blue", lty = 2)
    lines(sort(log_counts), pred[order(log_counts),"upr"], col = "blue", lty = 2)
    lines(sort(log_counts), pred[order(log_counts),"fit"], col = "blue", lwd = 2)
    
    legend("topleft",
           c("Cells", "Outliers", "Fit", "Prediction Interval"),
           col = c("grey50", "red", "blue", "blue"),
           pch = c(16, 16, NA, NA),
           lty = c(NA, NA, 1, 2),
           lwd = c(NA, NA, 2, 1))
    
    dev.off()
    
    # 输出每个样本的结果
    cat(sprintf("Sample %s: Detected %d outliers (%.1f%%)\n", 
                sample_name,
                sum(outliers), 
                mean(outliers) * 100))
    
    # 收集要保留的细胞
    cells_keep <- colnames(current_obj)[!outliers]
    all_cells_keep <- c(all_cells_keep, cells_keep)
  }
  
  # 返回所有要保留的细胞
  return(all_cells_keep)
}

# seurat_obj_main[["RNA"]] <- as(object = seurat_obj_main[["RNA"]], 
#                                Class = "Assay")

# 使用函数
# cells_to_keep <- regression_outliers_by_sample(seurat_obj_main)
# seurat_obj_main <- subset(seurat_obj_main, cells = cells_to_keep)
# seurat_obj_main <- JoinLayers(seurat_obj_main)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_obj_main <- CellCycleScoring(seurat_obj_main, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes)

# 2. 先进行标准化和特征选择
print("Normalizing data and finding variable features...")

# # table(seurat_obj_main@meta.data$sample[seurat_obj_main[["percent.mt"]] > 20])
# seurat_obj_main[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_main, pattern = "^MT-")
# seurat_obj_main <- subset(seurat_obj_main, subset = percent.mt <= 20)
# # 
# # 核糖体基因分析
# seurat_obj_main[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj_main,
#                                                           pattern = "^RP[SL]")
# 
# pdf("ribosomal_ratio_distribution.pdf", width = 8, height = 6)
# 
# ggplot(seurat_obj_main@meta.data, aes(x = percent.ribo)) +
#   geom_histogram(bins = 100, fill = "#4ECDC4", color = "white", alpha = 0.8) +
#   labs(title = "Ribosomal Gene Percentage Distribution",
#        x = "Ribosomal Ratio (%)",
#        y = "Cell Count") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
# 
# dev.off()
# 
# # 核糖体基因统计
# cat("\n=== Ribosomal Gene Statistics ===\n")
# summary(seurat_obj_main$percent.ribo)
# seurat_obj_main <- subset(seurat_obj_main, subset = percent.ribo <= 40)

# rm(doublet_results)
# seurat_obj_main <- readRDS("seurat_clean_0930.rds")
seurat_obj_main <- NormalizeData(seurat_obj_main)
seurat_obj_main <- FindVariableFeatures(seurat_obj_main, selection.method = "vst", nfeatures = 4000)
seurat_obj_main <- ScaleData(seurat_obj_main,vars.to.regress = c("S.Score", "G2M.Score",'percent.mt','percent.ribo'))
# seurat_obj_main <- SCTransform(seurat_obj_main, vars.to.regress = 'percent.mt')
seurat_obj_main <- RunPCA(seurat_obj_main,npcs = 40,verbose = TRUE)
pdf('dimmap.pdf')
DimHeatmap(seurat_obj_main, dims = 1:20, cells = 100, balanced = TRUE)
dev.off()
# #
# seurat_obj_main <- RenameGenesSeurat(
#   obj = seurat_obj_main,
#   newnames_file = "E:/R/sources/10x/featrues/3/merged_all_features.csv",
#   remove_duplicates = TRUE
# )
# saveRDS(seurat_obj_main,"final_analyzed_seurat_obj_main.rds")
# seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
# # 计算每个细胞的线粒体基因占比

seurat_obj_main <- RunHarmony(
  object = seurat_obj_main,           
  group.by.vars = c("sample"),   
  theta = c(1),      # Higher theta for more diverse clustering               
  lambda = c(6.7),   # Higher lambda to reduce overcorrection                
  sigma = 0.1,           # Lower sigma for tighter clusters
  nclust = 60,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 70, 
  # cool_down = 10,
  # epsilon_harmony = 0.0001,
  # monitor = "cluster_assignment",
  # factor = 0.5,
  # tol = 1e-9,
  # min_delta = 0.0001,
  # patience = 10,
  # # block.size = 100000,
  # epsilon_cluster = 0.00001,
  early_stop = TRUE,
  dims = 1:40           # More iterations for better convergence
)

# print("Calculating integration quality metrics...")
# 
# # 获取降维结果和元数据
# embeddings <- Embeddings(seurat_obj_main, "harmony")
# metadata <- seurat_obj_main@meta.data
# 
# # 过滤掉ann_level_2为NA的细胞
# valid_cells <- !is.na(metadata$ann_level_2)
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
#                                    c("sample", "ann_level_2"), 
#                                    k)
#     
#     # 存储结果
#     ilisi_scores[i] <- mean(lisi_res[, "sample"])
#     clisi_scores[i] <- mean(lisi_res[, "ann_level_2"])
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
# counts <- GetAssayData(seurat_obj_main, layer = "counts")
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
# # 3. PC肘部图
# pct_var <- seurat_obj_main[["pca"]]@stdev / sum(seurat_obj_main[["pca"]]@stdev) * 100
# cumsum_var <- cumsum(pct_var)
# 
# elbow_data <- data.frame(
#   PC = 1:length(pct_var),
#   variance = pct_var,
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
# # 4. Harmony scores knee plot - 使用采样
# set.seed(42)
# n_sample <- min(10000, ncol(seurat_obj_main))
# sample_idx <- sample(1:ncol(seurat_obj_main), n_sample)
# 
# harmony_scores <- seurat_obj_main[["harmony"]]@cell.embeddings[sample_idx, 1:20]
# harmony_dist <- dist(harmony_scores)
# harmony_scores_ordered <- sort(colMeans(as.matrix(harmony_dist)), decreasing = TRUE)
# 
# print(ggplot(data.frame(rank = 1:length(harmony_scores_ordered),
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
#              "PC_90percent_var", "Harmony_mean_dist"),
#   Value = c(
#     length(total_umi),
#     median(total_umi),
#     median(genes_per_cell),
#     which(cumsum_var > 90)[1],
#     mean(harmony_scores_ordered)
#   )
# )
# 
# write.csv(stats_df, "qc_plots/knee_plot_statistics.csv", row.names = FALSE)

# 直接运行UMAP
seurat_obj_main <- RunUMAP(seurat_obj_main, 
                           reduction = "harmony", 
                           dims = 1:40, 
                           n.neighbors = 30,
                           min.dist = 0.3,
                           # learning.rate = 0.2,    # 相对保守的学习率
                           # n.epochs = 1400,        # 增加迭代次数补偿较小的学习率
                           # spread = 1.2,
                           # repulsion.strength = 1.1,
                           
                           metric = "correlation")
# Find neighbors
seurat_obj_main <- FindNeighbors(seurat_obj_main, 
                                 reduction = "harmony", 
                                 dims = 1:40,
                                 n.trees = 500,
                                 k.param = 45) 

seurat_obj_main <- FindClusters(seurat_obj_main,
                                algorithm = 4,
                                group.singletons = FALSE,
                                resolution = 3, # 多个分辨率
                                verbose = TRUE)


# 8. 绘制UMAP聚类图
print("Generating plots...")
Idents(seurat_obj_main) <- seurat_obj_main$seurat_clusters
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(seurat_obj_main, reduction = "umap", label = TRUE)
dev.off()
# saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")
# 9. 保存分析结果
print("Saving analysis results...")

# seurat_obj_main <- readRDS("E:/R/0224/final_analyzed_seurat_obj_main.rds")


# 9. 标记基因分析

print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for(marker in Markers) {
  if(marker %in% rownames(seurat_obj_main)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(seurat_obj_main, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()
# gene_test <- c("CD44", "NT5E", "THY1", "ENG")
# pdf("marker_genes_by_tissue.pdf", width = 50, height = 12)
# 
# for(gene in gene_test) {
#   # 获取组织列表
#   tissue_levels <- unique(seurat_obj_main$tissue)
#   
#   p <- FeaturePlot(
#     object = seurat_obj_main,
#     features = gene,
#     reduction = "umap",
#     split.by = "tissue",
#     keep.scale = "feature",
#     cols = c("lightgrey", "blue"),
#     pt.size = 1,
#     order = TRUE,
#     combine = TRUE
#   )
#   
#   # 获取组合图中的子图数量
#   n_plots <- length(p$patches$plots)
#   
#   # 为每个子图添加正确的标题
#   for(i in 1:n_plots) {
#     if(i <= length(tissue_levels)) {
#       tissue_name <- tissue_levels[i]
#       p$patches$plots[[i]] <- p$patches$plots[[i]] + 
#         ggtitle(paste0(gene, " in ", tissue_name))
#     }
#   }
#   
#   print(p)
# }
# 
# dev.off()

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "ann_level_2", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_2")
p4 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p5 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "Annotation", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Annotation")
p1
p2
p3
p4
p5
dev.off()

tissues <- unique(seurat_obj_main$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(seurat_obj_main) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = seurat_obj_main, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(seurat_obj_main$sample)
pdf("sample_umap_split.pdf", width = 12, height = 12)
Idents(seurat_obj_main) <- "sample"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = seurat_obj_main, 
                                                    idents = sample_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(seurat_obj_main$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(seurat_obj_main) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = seurat_obj_main, 
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
seurat_obj_main$ann_level_2_no_na <- seurat_obj_main$ann_level_2
seurat_obj_main$ann_level_2_no_na[is.na(seurat_obj_main$ann_level_2_no_na)] <- "Unknown"
Idents(seurat_obj_main) <- "ann_level_2_no_na"

pdf("ann_level_2_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_level_2 <- unique(seurat_obj_main$ann_level_2)
ann_level_2 <- ann_level_2[!is.na(ann_level_2)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_level_2) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(seurat_obj_main, 
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
Markers <- unique(Markers)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot.pdf", width = 32, height = 18)
dot_plot <- DotPlot(seurat_obj_main, 
                    features = Markers, 
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

# 2. 计算每个cluster中每个基因的表达情况
print("Calculating expression statistics...")

# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
Idents(seurat_obj_main) <- seurat_obj_main$seurat_clusters
markers <- FindAllMarkers(seurat_obj_main, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

# table(subset(seurat_obj_main, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

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
pdf(filename, width = 32, height = 18)
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
# 

# seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
# 1. 读取注释文件（无表头）
annotations <- read.csv("annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- seurat_obj_main$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
seurat_obj_main$Annotation <- new_idents
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
        # cols = cell_colors
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 6. 绘制分割版本
pdf("cell_types_umap_split.pdf", width = 12, height = 12)
# 获取所有细胞类型
cell_types <- unique(seurat_obj_main$Annotation)
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
# table(seurat_obj_main@meta.data$tissue[seurat_obj_main@meta.data$seurat_clusters == "43"])
# table(seurat_obj_main@meta.data$sample[seurat_obj_main@meta.data$percent.mt >= 15])
# 5. 保存修改后的 Seurat 对象
# saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")

# table(seurat_obj_main@meta.data$seurat_clusters[seurat_obj_main@meta.data$Annotation=='Epithelial' & seurat_obj_main@meta.data$ann_level_2 =='Myeloid'])

# 假设我们选择 cell_type 为 "T_cells" 的细胞
# 1. 提取Undefined细胞子集
# 从主Seurat对象中筛选出Annotation标注为"Undefined"的细胞
Undefined_obj <- subset(seurat_obj_main, subset = Annotation == "Undefined")

# 2. 标准数据预处理
# 数据标准化，默认使用"LogNormalize"方法
Undefined_obj <- NormalizeData(Undefined_obj)

# 3. 寻找高变异基因
# 使用方差稳定化转换(VST)方法选择3000个高变基因，用于后续降维分析
Undefined_obj <- FindVariableFeatures(Undefined_obj, selection.method = "vst", nfeatures = 3000)

# 4. 数据缩放
# 对所有高变基因进行归一化处理，使其均值为0，方差为1
Undefined_obj <- ScaleData(Undefined_obj, features = VariableFeatures(Undefined_obj))

# 5. 主成分分析
# 使用高变基因进行主成分分析，降低数据维度
Undefined_obj <- RunPCA(Undefined_obj)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
Undefined_obj <- RunHarmony(
  object = Undefined_obj,           
  group.by.vars = c("sample"),   
  theta = c(1),      # Higher theta for more diverse clustering               
  lambda = c(7),   # Higher lambda to reduce overcorrection                
  sigma = 0.1,           # Lower sigma for tighter clusters
  nclust = 30,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20, 
  early_stop = TRUE,
  dims = 1:40           # More iterations for better convergence
)
Undefined_obj <- RunUMAP(Undefined_obj, 
                           reduction = "harmony", 
                           dims = 1:30, 
                           n.neighbors = 30,
                           n.trees = 500,
                           min.dist = 0.4,
                           learning.rate = 0.2,    # 相对保守的学习率
                           n.epochs = 1400,        # 增加迭代次数补偿较小的学习率
                           spread = 1.2,
                           repulsion.strength = 1.1,
                           
                           metric = "correlation")
# Find neighbors
Undefined_obj <- FindNeighbors(Undefined_obj, 
                                 reduction = "harmony", 
                                 dims = 1:30,
                                 k.param = 45) 

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
for(marker in Markers) {
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

# 15. 批次整合质量评估
print("评估批次整合质量...")
# 获取harmony降维结果和元数据
embeddings <- Embeddings(Undefined_obj, "harmony")
metadata <- Undefined_obj@meta.data

# 抽样计算LISI分数以评估批次整合质量
set.seed(42)
n_cells_sample <- min(5000, nrow(embeddings))
sample_idx <- sample(nrow(embeddings), n_cells_sample)

# 使用lisi包计算整合得分
if(requireNamespace("lisi", quietly = TRUE)) {
  lisi_res <- lisi::compute_lisi(embeddings[sample_idx, 1:30], 
                                 metadata[sample_idx, ], 
                                 c("sample", "tissue"), 
                                 30)
  
  # 保存整合质量指标
  cat(sprintf("\n批次整合质量评估:\n"))
  cat(sprintf("- 样本整合得分(iLISI): %.3f\n", mean(lisi_res[, "sample"])))
  cat(sprintf("- 组织保真度(cLISI): %.3f\n", mean(lisi_res[, "tissue"])))
  
  # 保存到文件
  integration_metrics <- data.frame(
    Metric = c("iLISI_sample", "cLISI_tissue"),
    Value = c(mean(lisi_res[, "sample"]), mean(lisi_res[, "tissue"]))
  )
  write.csv(integration_metrics, "Undefined_integration_metrics.csv", row.names = FALSE)
} else {
  warning("lisi包未安装，跳过批次整合质量评估")
}

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
                    features = Markers, 
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
# 9. 保存分析结果
print("Saving analysis results...")
# saveRDS(Epithelial_obj, "final_analyzed_Epithelial_obj.rds")

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
cells_in_main <- cells_to_update %in% rownames(seurat_obj_main@meta.data)
valid_cells_to_update <- cells_to_update[cells_in_main]

# 只更新主对象中存在的细胞的Annotation
seurat_obj_main@meta.data[valid_cells_to_update, "Annotation"] <- 
  Undefined_obj@meta.data[valid_cells_to_update, "Annotation"]

# 验证更新
table(seurat_obj_main@meta.data$Annotation)

Idents(seurat_obj_main) <- "Annotation"
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(seurat_obj_main, 
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
rm(Undefined_obj)
gc()

table(seurat_obj_main$Annotation)

# 将"Epthelial"替换为"Epithelial"
seurat_obj_main$Annotation[seurat_obj_main$Annotation == "Remove"] <- "B"

# 确认更改已生效
table(seurat_obj_main$Annotation)

seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
# annotations <- read.csv("annotation.csv", header = FALSE)
# colnames(annotations) <- c("Cluster", "Annotation")

print("Analysis completed successfully!")

results <- analyze_cell_proportions_complete(
  seurat_obj = seurat_obj_main,
  exclude_filter = "tissue_sampling_method != 'scraping'",
  output_dir = "proportion_results"
)

# 获取所有细胞类型
cell_column <- "Annotation"  # 根据您的数据调整
cell_types <- unique(seurat_obj_main@meta.data[[cell_column]])
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 批量处理每种细胞类型
for (cell_type in cell_types) {
  new_var_name <- paste(cell_type, 'obj', sep = "_")
  Marker_current <- get(paste0('Marker_',cell_type))
  assign(new_var_name, run_subset_analysis(
    seurat_obj_main,
    cell_type = cell_type,
    markers = Marker_current,
    output_dir = paste(cell_type, "analysis", sep = "_")
  ))
  
  # 清理内存
  gc()
}

results <- scPairwiseMASCAnalysis(
  seurat_obj = seurat_obj_main,
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  exclude_filter = "tissue_sampling_method != 'scraping'",
  fixed_effects_cols = NULL,  # 添加任何固定效应协变量
  output_dir = "pairwise_MASC_results",
  # 如果您需要少于2个样本也能进行比较，可以降低此阈值
  min_samples = 2
)

# 首先创建质控图的输出目录
dir.create("qaqc_plots", showWarnings = TRUE)

# 获取所有样本名称
sample_names <- unique(seurat_obj_main$sample)
Idents(seurat_obj_main) <- "sample"
# 为每个样本生成质控图
for(sample_name in sample_names) {
  # 提取该样本的细胞
  sample_cells <- WhichCells(seurat_obj_main, expression = sample == sample_name)
  
  # 如果该样本有细胞，则继续处理
  if(length(sample_cells) > 0) {
    # 创建该样本的子集
    sample_obj <- subset(seurat_obj_main, cells = sample_cells)
    
    # 生成质控图 - 小提琴图
    p1 <- VlnPlot(sample_obj, 
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, 
                  pt.size = 0.1)  # 减小点的大小以提高可读性
    
    # 散点图 - RNA计数与基因数量
    p2 <- FeatureScatter(sample_obj, 
                         feature1 = "nCount_RNA", 
                         feature2 = "nFeature_RNA")
    
    # 散点图 - RNA计数与线粒体百分比
    p3 <- FeatureScatter(sample_obj, 
                         feature1 = "nCount_RNA", 
                         feature2 = "percent.mt")
    
    # 保存质控图
    pdf(file.path("qaqc_plots", paste0(sample_name, "_qc.pdf")),
        width = 15, height = 10)
    print(p1)
    print(p2 + p3)
    dev.off()
    
    # 计算并输出该样本的汇总统计信息
    stats <- data.frame(
      Sample = sample_name,
      Cells = ncol(sample_obj),
      Median_Genes = median(sample_obj$nFeature_RNA),
      Mean_Genes = mean(sample_obj$nFeature_RNA),
      Median_UMIs = median(sample_obj$nCount_RNA),
      Mean_UMIs = mean(sample_obj$nCount_RNA),
      Median_Mt_Percent = median(sample_obj$percent.mt),
      Mean_Mt_Percent = mean(sample_obj$percent.mt)
    )
    
    # 将统计信息写入文件
    if(!exists("all_stats")) {
      all_stats <- stats
    } else {
      all_stats <- rbind(all_stats, stats)
    }
    
    cat(sprintf("Generated QC plots for sample: %s (Cells: %d)\n", 
                sample_name, ncol(sample_obj)))
  } else {
    cat(sprintf("Warning: No cells found for sample: %s\n", sample_name))
  }
}

# 保存所有样本的统计数据
write.csv(all_stats, "qaqc_plots/sample_qc_metrics.csv", row.names = FALSE)

###############################################################################################
###############################################################################################
###############################################################################################

# Optimized QC Plot Generation - Preprocessing Data to Avoid Repeated Subset Operations
# Author: Bioinformatics Analysis
# Date: 2025
# Performance Optimization: Pre-extract sample data to avoid repeated subset in loops

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Create QC plots output directory
dir.create("qc_plots_by_metric", showWarnings = FALSE)

# ============================================================================
# Data Preprocessing Stage - Extract All Sample Data at Once
# ============================================================================
cat("=== Starting Data Preprocessing ===\n")

# Get all sample names
sample_names <- unique(seurat_obj_main$sample)
cat(sprintf("Detected %d samples: %s\n", 
            length(sample_names), 
            paste(sample_names, collapse = ", ")))

# Pre-extract QC data for each sample
sample_data_list <- list()
sample_stats_list <- list()

cat("Extracting data for each sample...\n")
for(sample_name in sample_names) {
  cat(sprintf("  Processing sample: %s", sample_name))
  
  # Extract cell indices for this sample
  sample_cells <- WhichCells(seurat_obj_main, expression = sample == sample_name)
  
  if(length(sample_cells) > 0) {
    # Extract QC metrics data (avoid creating Seurat subobjects)
    sample_qc_data <- data.frame(
      cell_id = sample_cells,
      sample = sample_name,
      nFeature_RNA = seurat_obj_main$nFeature_RNA[sample_cells],
      nCount_RNA = seurat_obj_main$nCount_RNA[sample_cells],
      percent.mt = seurat_obj_main$percent.mt[sample_cells]
    )
    
    # Store data
    sample_data_list[[sample_name]] <- sample_qc_data
    
    # Pre-calculate statistics
    stats <- data.frame(
      Sample = sample_name,
      Cells = length(sample_cells),
      Median_Genes = median(sample_qc_data$nFeature_RNA, na.rm = TRUE),
      Mean_Genes = round(mean(sample_qc_data$nFeature_RNA, na.rm = TRUE), 1),
      SD_Genes = round(sd(sample_qc_data$nFeature_RNA, na.rm = TRUE), 1),
      Median_UMIs = median(sample_qc_data$nCount_RNA, na.rm = TRUE),
      Mean_UMIs = round(mean(sample_qc_data$nCount_RNA, na.rm = TRUE), 1),
      SD_UMIs = round(sd(sample_qc_data$nCount_RNA, na.rm = TRUE), 1),
      Median_Mt_Percent = round(median(sample_qc_data$percent.mt, na.rm = TRUE), 2),
      Mean_Mt_Percent = round(mean(sample_qc_data$percent.mt, na.rm = TRUE), 2),
      SD_Mt_Percent = round(sd(sample_qc_data$percent.mt, na.rm = TRUE), 2)
    )
    sample_stats_list[[sample_name]] <- stats
    
    cat(sprintf(" ✓ (Cells: %d)\n", nrow(sample_qc_data)))
  } else {
    cat(" ✗ (No cells)\n")
    # Create empty placeholder dataframe
    sample_data_list[[sample_name]] <- data.frame()
    sample_stats_list[[sample_name]] <- data.frame(
      Sample = sample_name, Cells = 0,
      Median_Genes = NA, Mean_Genes = NA, SD_Genes = NA,
      Median_UMIs = NA, Mean_UMIs = NA, SD_UMIs = NA,
      Median_Mt_Percent = NA, Mean_Mt_Percent = NA, SD_Mt_Percent = NA
    )
  }
}

# Combine statistics data
all_stats <- do.call(rbind, sample_stats_list)
rownames(all_stats) <- NULL

cat("Data preprocessing completed!\n\n")

# ============================================================================
# Calculate Layout Parameters and Global Axis Ranges
# ============================================================================
n_samples <- length(sample_names)
n_cols <- min(4, n_samples)  # Maximum 4 columns
n_rows <- ceiling(n_samples / n_cols)

cat(sprintf("Plot layout: %d rows × %d cols (total %d samples)\n\n", n_rows, n_cols, n_samples))

# Calculate global axis ranges for consistent scaling across all plots
all_qc_data_valid <- do.call(rbind, sample_data_list[sapply(sample_data_list, nrow) > 0])

if(nrow(all_qc_data_valid) > 0) {
  # Global ranges with some padding
  nFeature_range <- range(all_qc_data_valid$nFeature_RNA, na.rm = TRUE)
  nFeature_range <- c(nFeature_range[1] * 0.95, nFeature_range[2] * 1.05)
  
  nCount_range <- range(all_qc_data_valid$nCount_RNA, na.rm = TRUE)  
  nCount_range <- c(nCount_range[1] * 0.95, nCount_range[2] * 1.05)
  
  percent_mt_range <- range(all_qc_data_valid$percent.mt, na.rm = TRUE)
  percent_mt_range <- c(max(0, percent_mt_range[1] - 1), percent_mt_range[2] + 1)
  
  cat("Global axis ranges calculated:\n")
  cat(sprintf("  nFeature_RNA: %.0f - %.0f\n", nFeature_range[1], nFeature_range[2]))
  cat(sprintf("  nCount_RNA: %.0f - %.0f\n", nCount_range[1], nCount_range[2]))
  cat(sprintf("  percent.mt: %.1f - %.1f\n", percent_mt_range[1], percent_mt_range[2]))
  cat("\n")
} else {
  # Default ranges if no valid data
  nFeature_range <- c(0, 5000)
  nCount_range <- c(0, 50000)
  percent_mt_range <- c(0, 100)
}

# ============================================================================
# Define Plotting Functions - Unified Plotting Logic
# ============================================================================

# Violin plot function with consistent axes
create_violin_plot <- function(sample_name, feature, ylabel, y_range) {
  if(sample_name %in% names(sample_data_list) && nrow(sample_data_list[[sample_name]]) > 0) {
    data <- sample_data_list[[sample_name]]
    
    p <- ggplot(data, aes_string(x = "'Sample'", y = feature)) +
      geom_violin(fill = rainbow(1), alpha = 0.7, trim = TRUE) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 0.1) +
      stat_summary(fun = median, geom = "point", size = 2, color = "red") +
      ylim(y_range) +  # Apply consistent y-axis range
      ggtitle(paste("Sample:", sample_name)) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      ylab(ylabel)
    return(p)
  } else {
    # Blank plot for samples with no data
    return(ggplot() + 
             ggtitle(paste("Sample:", sample_name, "(No cells)")) +
             ylim(y_range) +  # Keep consistent axes even for empty plots
             theme_void() +
             theme(plot.title = element_text(size = 12, hjust = 0.5)))
  }
}

# Scatter plot function with consistent axes
create_scatter_plot <- function(sample_name, x_feature, y_feature, xlabel, ylabel, x_range, y_range) {
  if(sample_name %in% names(sample_data_list) && nrow(sample_data_list[[sample_name]]) > 0) {
    data <- sample_data_list[[sample_name]]
    
    p <- ggplot(data, aes_string(x = x_feature, y = y_feature)) +
      geom_point(alpha = 0.6, size = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red", size = 0.8) +
      xlim(x_range) + ylim(y_range) +  # Apply consistent axis ranges
      ggtitle(paste("Sample:", sample_name)) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, hjust = 0.5)) +
      xlab(xlabel) + ylab(ylabel)
    return(p)
  } else {
    return(ggplot() + 
             ggtitle(paste("Sample:", sample_name, "(No cells)")) +
             xlim(x_range) + ylim(y_range) +  # Keep consistent axes
             theme_void() +
             theme(plot.title = element_text(size = 12, hjust = 0.5)))
  }
}

# ============================================================================
# Generate All Types of QC Plots
# ============================================================================

# 1. Gene count violin plots
cat("Generating gene count violin plots (nFeature_RNA)...\n")
violin_plots_genes <- lapply(sample_names, create_violin_plot, 
                             feature = "nFeature_RNA", 
                             ylabel = "Number of Genes (nFeature_RNA)",
                             y_range = nFeature_range)
genes_combined <- wrap_plots(violin_plots_genes, ncol = n_cols, nrow = n_rows)

# 2. UMI count violin plots
cat("Generating UMI count violin plots (nCount_RNA)...\n")
violin_plots_counts <- lapply(sample_names, create_violin_plot, 
                              feature = "nCount_RNA", 
                              ylabel = "UMI Counts (nCount_RNA)",
                              y_range = nCount_range)
counts_combined <- wrap_plots(violin_plots_counts, ncol = n_cols, nrow = n_rows)

# 3. Mitochondrial percentage violin plots
cat("Generating mitochondrial percentage violin plots (percent.mt)...\n")
violin_plots_mt <- lapply(sample_names, create_violin_plot, 
                          feature = "percent.mt", 
                          ylabel = "Mitochondrial Gene Percentage (%)",
                          y_range = percent_mt_range)
mt_combined <- wrap_plots(violin_plots_mt, ncol = n_cols, nrow = n_rows)

# 4. UMI vs Gene count scatter plots
cat("Generating UMI vs Gene count scatter plots...\n")
scatter_plots_umi_genes <- lapply(sample_names, create_scatter_plot,
                                  x_feature = "nCount_RNA", 
                                  y_feature = "nFeature_RNA",
                                  xlabel = "UMI Counts (nCount_RNA)",
                                  ylabel = "Number of Genes (nFeature_RNA)",
                                  x_range = nCount_range,
                                  y_range = nFeature_range)
umi_genes_combined <- wrap_plots(scatter_plots_umi_genes, ncol = n_cols, nrow = n_rows)

# 5. UMI vs Mitochondrial percentage scatter plots
cat("Generating UMI vs Mitochondrial percentage scatter plots...\n")
scatter_plots_umi_mt <- lapply(sample_names, create_scatter_plot,
                               x_feature = "nCount_RNA", 
                               y_feature = "percent.mt",
                               xlabel = "UMI Counts (nCount_RNA)",
                               ylabel = "Mitochondrial Gene Percentage (%)",
                               x_range = nCount_range,
                               y_range = percent_mt_range)
umi_mt_combined <- wrap_plots(scatter_plots_umi_mt, ncol = n_cols, nrow = n_rows)

# ============================================================================
# Save All Plots
# ============================================================================
cat("Saving plot files...\n")

# Set plot dimensions
fig_width <- 4 * n_cols
fig_height <- 4 * n_rows

# Save all types of plots with limitsize = FALSE
ggsave("qc_plots_by_metric/nFeature_RNA_all_samples.pdf", 
       plot = genes_combined, 
       width = fig_width, height = fig_height, units = "in", limitsize = FALSE)

ggsave("qc_plots_by_metric/nCount_RNA_all_samples.pdf", 
       plot = counts_combined, 
       width = fig_width, height = fig_height, units = "in", limitsize = FALSE)

ggsave("qc_plots_by_metric/percent_mt_all_samples.pdf", 
       plot = mt_combined, 
       width = fig_width, height = fig_height, units = "in", limitsize = FALSE)

ggsave("qc_plots_by_metric/scatter_UMI_vs_Genes_all_samples.pdf", 
       plot = umi_genes_combined, 
       width = fig_width, height = fig_height, units = "in", limitsize = FALSE)

ggsave("qc_plots_by_metric/scatter_UMI_vs_Mt_all_samples.pdf", 
       plot = umi_mt_combined, 
       width = fig_width, height = fig_height, units = "in", limitsize = FALSE)

# Save statistics data
write.csv(all_stats, "qc_plots_by_metric/sample_qc_metrics_summary.csv", 
          row.names = FALSE)

# Optional: Save raw QC data for further analysis
all_qc_data <- do.call(rbind, sample_data_list[sapply(sample_data_list, nrow) > 0])
write.csv(all_qc_data, "qc_plots_by_metric/all_samples_qc_raw_data.csv", 
          row.names = FALSE)

# ============================================================================
# Output Completion Information
# ============================================================================
cat("\n=== Quality Control Analysis Completed ===\n")
cat(sprintf("Successfully processed %d samples\n", sum(all_stats$Cells > 0)))
cat(sprintf("Total cell count: %s\n", format(sum(all_stats$Cells), big.mark = ",")))

cat("\nGenerated files:\n")
cat("  Plot files:\n")
cat("    - nFeature_RNA_all_samples.pdf (Gene count comparison)\n")
cat("    - nCount_RNA_all_samples.pdf (UMI count comparison)\n") 
cat("    - percent_mt_all_samples.pdf (Mitochondrial percentage comparison)\n")
cat("    - scatter_UMI_vs_Genes_all_samples.pdf (UMI-Gene correlation)\n")
cat("    - scatter_UMI_vs_Mt_all_samples.pdf (UMI-Mitochondrial correlation)\n")
cat("  Data files:\n")
cat("    - sample_qc_metrics_summary.csv (Statistical summary)\n")
cat("    - all_samples_qc_raw_data.csv (Raw QC data)\n")

cat(sprintf("\nAll files saved in: qc_plots_by_metric/\n"))
cat("Plot dimensions:", fig_width, "×", fig_height, "inches\n")

# Display statistical summary
cat("\nSample QC Statistics Summary:\n")
print(all_stats)

###############################################################################################
###############################################################################################
###############################################################################################

# 批量分析多个细胞类型
all_results <- run_all_celltypes_tissue_comparisons(
  seurat_obj = seurat_obj_main,
  cell_types = c( "Fibroblast", "T", "B", "Myeloid","Endothelial"),#"SMC","Epithelial",
  cell_type_col = "Annotation",
  tissue_col = "tissue",
  base_output_dir = "all_tissue_comparisons",
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  min_cells = 50
)

epithelial_results_1 <- run_specific_tissue_comparison(
  seurat_obj = seurat_obj_main,
  cell_type = "Epithelial",
  tissue1 = "lung parenchyma",
  tissue2 = "sinus",
  cell_type_col = "Annotation",
  tissue_col = "tissue",
  output_dir = "Epithelial_lung_vs_sinus",
  log2fc_cutoff = 1.5,
  pval_cutoff = 0.05
)
# seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
epithelial_results <- run_specific_tissue_comparison(
  seurat_obj = seurat_obj_main,
  cell_type = "Epithelial",
  tissue1 = "nose",
  tissue2 = "sinus",
  cell_type_col = "Annotation",
  tissue_col = "tissue",
  output_dir = "Epithelial_nose_vs_sinus",
  log2fc_cutoff = 1.5,
  pval_cutoff = 0.05
)

library(dplyr)
library(ggplot2)
library(ggalluvial)
library(scales)
library(RColorBrewer)
library(colorspace)  # 用于颜色调整
library(tidyr)       # 替代gather函数

# 从元数据创建数据框
meta_data <- as.data.frame(seurat_obj_main@meta.data)

# 以样本为单位汇总数据
sample_summary <- meta_data %>% 
  group_by(sample, study, tissue, tissue_sampling_method) %>%
  summarise(cell_count = n(), .groups = "drop")

# 为ggalluvial准备数据并按研究排序样本
alluv_data <- sample_summary %>%
  select(sample, study, tissue, tissue_sampling_method, cell_count) %>%
  # 为样本创建排序因子
  mutate(
    study_factor = factor(study),
    sample_ordered = factor(sample, 
                            levels = unique(sample[order(study_factor)]))
  )

# 计算每种类别的唯一值数量
n_samples <- length(unique(alluv_data$sample))
n_studies <- length(unique(alluv_data$study))
n_tissues <- length(unique(alluv_data$tissue))
n_methods <- length(unique(alluv_data$tissue_sampling_method))

# 确定唯一的组织类型
unique_tissues <- unique(alluv_data$tissue)
n_unique_tissues <- length(unique_tissues)

# 为每种组织类型创建一个基础颜色
tissue_base_colors <- brewer.pal(min(n_unique_tissues, 8), "Set2")
if(n_unique_tissues > 8) {
  tissue_base_colors <- colorRampPalette(tissue_base_colors)(n_unique_tissues)
}
names(tissue_base_colors) <- unique_tissues

# 创建基于tissue归属的颜色映射
# 样本颜色 - 基于其所属的tissue
sample_tissue_map <- alluv_data %>%
  select(sample_ordered, tissue) %>%
  distinct()

sample_colors <- sapply(sample_tissue_map$tissue, function(t) {
  base <- tissue_base_colors[t]
  lighten(base, amount = 0.3)  # 使用colorspace包的lighten函数
})
names(sample_colors) <- sample_tissue_map$sample_ordered

# 研究颜色 - 基于主要关联的tissue
study_tissue_counts <- alluv_data %>%
  group_by(study, tissue) %>%
  summarise(count = sum(cell_count), .groups = "drop") %>%
  arrange(desc(count))

study_primary_tissue <- study_tissue_counts %>%
  group_by(study) %>%
  slice(1) %>%
  select(study, tissue)

study_colors <- sapply(study_primary_tissue$tissue, function(t) {
  base <- tissue_base_colors[t]
  darken(base, amount = 0.1)  # 微暗化
})
names(study_colors) <- study_primary_tissue$study

# 采样方法颜色 - 基于主要关联的tissue
method_tissue_counts <- alluv_data %>%
  group_by(tissue_sampling_method, tissue) %>%
  summarise(count = sum(cell_count), .groups = "drop") %>%
  arrange(desc(count))

method_primary_tissue <- method_tissue_counts %>%
  group_by(tissue_sampling_method) %>%
  slice(1) %>%
  select(tissue_sampling_method, tissue)

method_colors <- sapply(method_primary_tissue$tissue, function(t) {
  base <- tissue_base_colors[t]
  darken(base, amount = 0.2)  # 进一步暗化
})
names(method_colors) <- method_primary_tissue$tissue_sampling_method

# 创建标注版桑基图
p <- ggplot(alluv_data,
            aes(y = cell_count, axis1 = sample_ordered, axis2 = study, axis3 = tissue, axis4 = tissue_sampling_method)) +
  # 流线 - 半透明，按tissue着色
  geom_alluvium(aes(fill = tissue), alpha = 0.5, width = 1/3) +
  # 节点 - 各层级不同颜色
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/3, color = "grey80") +
  # 标签 - 所有节点，完全移到左侧
  geom_text(stat = "stratum", 
            aes(label = paste0(after_stat(stratum), " (", comma(after_stat(count)), ")")),
            size = 2.3,
            hjust = 1,
            nudge_x = -0.2,  # 增大向左偏移
            check_overlap = TRUE) +
  # 设置填充颜色，包含所有类别的颜色
  scale_fill_manual(values = c(
    sample_colors,
    study_colors,
    tissue_base_colors,
    method_colors
  ), guide = "none") +
  # 坐标轴设置
  scale_x_discrete(limits = c("Sample", "Study", "Tissue", "Sampling Method"),
                   labels = c("Sample\n(by Study)", "Study", "Tissue", "Sampling Method")) +  # 增加左侧空间
  scale_y_continuous(labels = comma) +
  # 主题设置
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.y = element_text(size = 10, family = "Times New Roman"),
    axis.text.x = element_text(size = 12, face = "bold", family = "Times New Roman"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold", family = "Times New Roman"),
    text = element_text(family = "Times New Roman")  # 明确指定为Times New Roman
  ) +
  labs(
    y = "Cell Count"
  )

# 保存图表
figure_height <- max(10, n_samples * 0.15)  # 根据样本数量动态调整高度
figure_width <- 18  # 增加宽度，为左侧标签腾出空间
ggsave("static_sankey_complete.png", p, width = figure_width, height = figure_height, dpi = 300)
ggsave("static_sankey_complete.pdf", p, width = figure_width, height = figure_height)

# 输出节点计数信息
cat("样本数量:", n_samples, "\n")
cat("研究数量:", n_studies, "\n")
cat("组织类型数量:", n_tissues, "\n")
cat("采样方法数量:", n_methods, "\n")
cat("总细胞数:", comma(sum(alluv_data$cell_count)), "\n")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
all_tissues <- unique(seurat_obj_main[[tissue_column]])
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
    cells_to_keep <- which(seurat_obj_main@meta.data[[tissue_column]] == tissue)
    if(length(cells_to_keep) == 0) {
      warning(paste0("在组织列中找不到匹配的组织: ", tissue))
      next
    }
    tissue_seurat <- seurat_obj_main[, cells_to_keep]
    
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
  seurat_obj = seurat_obj_main,
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
# check_seurat_structure(seurat_obj_main, "tissue", "Annotation")

# 基本调用示例
# results <- create_multiCondition_dimplot(
#   seurat_obj = seurat_obj_main,
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
#   seurat_obj = seurat_obj_main,
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
#   seurat_obj = seurat_obj_main,
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
check_seurat_structure(seurat_obj_main, "tissue", "Annotation")

# 基本调用
results <- create_multiCondition_dimplot(
  seurat_obj = seurat_obj_main,
  condition_col = "tissue",
  cell_type_col = "Annotation",
  reduction = 'umap',
  pt_size = 0.45,
  alpha = 2
)

# 快速预览
preview_plot <- quick_preview(seurat_obj_main, "tissue", "Annotation")
print(preview_plot)

######################################################################################################
######################################################################################################
######################################################################################################
diagnosis <- diagnose_masc_convergence(seurat_obj_main)

results <- scPairwiseMASCAnalysis(
  seurat_obj = seurat_obj_main,
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
  seurat_obj = seurat_obj_main,  # 备用重新计算
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
tissue_sample_counts <- sapply(unique(seurat_obj_main$tissue), function(t) {
  length(unique(seurat_obj_main$sample[seurat_obj_main$tissue == t]))
})
names(tissue_sample_counts) <- unique(seurat_obj_main$tissue)
print(tissue_sample_counts)

# 或者使用table + unique组合
unique_combinations <- unique(seurat_obj_main@meta.data[, c("tissue", "sample")])
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
    max_cells_per_type = 2000,
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
  seurat_obj = seurat_obj_main,
  cell_type_col = "Annotation",
  output_dir = "publication_heatmaps"
)
#
# # 2. 使用特定基因列表
# marker_genes <- c("CD3E", "CD3D", "CD14", "CD68", "EPCAM", "PECAM1", "COL1A1")
# result <- create_comprehensive_expression_heatmap(
#   seurat_obj = seurat_obj_main,
#   cell_type_col = "Annotation", 
#   genes_to_plot = marker_genes,
#   max_cells_per_type = 300,
#   output_prefix = "marker_expression"
# )
#
# # 3. 快速预览
# preview <- quick_preview_heatmap(seurat_obj_main)
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
#   seurat_obj = seurat_obj_main,
#   cell_type_col = "Annotation",
#   cell_type_colors = custom_colors,
#   expression_color_scheme = c("#2166AC", "#F7F7F7", "#B2182B")
# )

######################################################################################################
######################################################################################################
######################################################################################################

# ============================================================================
# Seurat对象基因表达质量诊断脚本
# 检测Harmony整合后填充0基因的影响
# ============================================================================

diagnose_seurat_gene_quality <- function(
    seurat_obj,
    celltype_col = "Annotation_2",
    tissue_col = "tissue",
    sample_col = "sample",
    assay = "RNA",
    slot = "counts"
) {
  
  message("╔═══════════════════════════════════════════════════════════╗")
  message("║     Seurat对象基因表达质量诊断                              ║")
  message("║     Detecting Zero-Padded Genes After Harmony Integration ║")
  message("╚═══════════════════════════════════════════════════════════╝\n")
  
  # 检查必需列
  required_cols <- c(celltype_col, tissue_col, sample_col)
  missing_cols <- required_cols[!required_cols %in% colnames(seurat_obj@meta.data)]
  if(length(missing_cols) > 0) {
    stop(paste("缺少必需列:", paste(missing_cols, collapse = ", ")))
  }
  
  # 提取counts矩阵
  if(slot == "counts") {
    counts_matrix <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
  } else {
    counts_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  }
  
  # 转换为密集矩阵（如果是稀疏矩阵）
  if(inherits(counts_matrix, "dgCMatrix")) {
    message("检测到稀疏矩阵，转换为密集矩阵进行诊断...")
    counts_matrix <- as.matrix(counts_matrix)
  }
  
  # 基本信息
  total_genes <- nrow(counts_matrix)
  total_cells <- ncol(counts_matrix)
  
  message(sprintf("数据集概况:"))
  message(sprintf("  基因数: %d", total_genes))
  message(sprintf("  细胞数: %d", total_cells))
  message(sprintf("  细胞类型数: %d", length(unique(seurat_obj@meta.data[[celltype_col]]))))
  message(sprintf("  组织类型数: %d", length(unique(seurat_obj@meta.data[[tissue_col]]))))
  message(sprintf("  样本数: %d\n", length(unique(seurat_obj@meta.data[[sample_col]]))))
  
  # ========================================================================
  # 1. 全局基因质量评估
  # ========================================================================
  message("【1】全局基因表达质量分析")
  message("─────────────────────────────────────")
  
  # 计算基因统计
  gene_stats <- data.frame(
    gene = rownames(counts_matrix),
    total_count = Matrix::rowSums(counts_matrix),
    mean_expr = Matrix::rowMeans(counts_matrix),
    nonzero_cells = Matrix::rowSums(counts_matrix > 0),
    nonzero_prop = Matrix::rowSums(counts_matrix > 0) / total_cells,
    max_expr = apply(counts_matrix, 1, max),
    stringsAsFactors = FALSE
  )
  
  # 零值统计
  total_values <- as.numeric(total_genes) * total_cells
  zero_count <- sum(counts_matrix == 0)
  zero_prop <- zero_count / total_values
  
  message(sprintf("零值比例: %.2f%% (%s / %s)", 
                  zero_prop * 100,
                  format(zero_count, big.mark = ","),
                  format(total_values, big.mark = ",")))
  
  # 基因表达频率分类
  gene_stats$category <- cut(
    gene_stats$nonzero_prop,
    breaks = c(0, 0.01, 0.05, 0.1, 0.3, 0.5, 1),
    labels = c("Extremely Sparse (<1%)", 
               "Very Sparse (1-5%)", 
               "Sparse (5-10%)",
               "Low (10-30%)", 
               "Moderate (30-50%)", 
               "Common (>50%)")
  )
  
  message("\n基因表达频率分布:")
  freq_table <- table(gene_stats$category)
  for(i in seq_along(freq_table)) {
    cat(sprintf("  %-25s: %6d (%.1f%%)\n", 
                names(freq_table)[i], 
                freq_table[i], 
                freq_table[i]/total_genes * 100))
  }
  
  # 识别可能的填充基因
  likely_padded <- gene_stats[
    gene_stats$nonzero_prop < 0.05 & 
      gene_stats$max_expr < 10,
  ]
  
  message(sprintf("\n⚠️  可能的填充0基因: %d (%.1f%%)", 
                  nrow(likely_padded), 
                  nrow(likely_padded)/total_genes * 100))
  
  if(nrow(likely_padded) > 0) {
    message("   示例基因 (前10个):")
    example_genes <- head(likely_padded$gene, 10)
    cat(sprintf("   %s\n", paste(example_genes, collapse = ", ")))
  }
  
  # ========================================================================
  # 2. 按细胞类型诊断
  # ========================================================================
  message("\n【2】按细胞类型的基因表达质量")
  message("─────────────────────────────────────")
  
  cell_types <- unique(seurat_obj@meta.data[[celltype_col]])
  celltype_stats <- list()
  
  for(ct in cell_types) {
    ct_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[celltype_col]] == ct]
    ct_counts <- counts_matrix[, ct_cells, drop = FALSE]
    
    ct_zero_prop <- sum(ct_counts == 0) / length(ct_counts)
    ct_expressed_genes <- sum(Matrix::rowSums(ct_counts > 0) > 0)
    
    celltype_stats[[ct]] <- list(
      n_cells = length(ct_cells),
      zero_prop = ct_zero_prop,
      expressed_genes = ct_expressed_genes
    )
    
    message(sprintf("  %-20s: %5d cells, %.1f%% 零值, %5d 表达基因", 
                    ct, 
                    length(ct_cells),
                    ct_zero_prop * 100,
                    ct_expressed_genes))
  }
  
  # ========================================================================
  # 3. 按组织诊断
  # ========================================================================
  message("\n【3】按组织的基因表达质量")
  message("─────────────────────────────────────")
  
  tissues <- unique(seurat_obj@meta.data[[tissue_col]])
  tissue_stats <- list()
  
  for(tissue in tissues) {
    tissue_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[tissue_col]] == tissue]
    tissue_counts <- counts_matrix[, tissue_cells, drop = FALSE]
    
    tissue_zero_prop <- sum(tissue_counts == 0) / length(tissue_counts)
    tissue_expressed_genes <- sum(Matrix::rowSums(tissue_counts > 0) > 0)
    
    tissue_stats[[tissue]] <- list(
      n_cells = length(tissue_cells),
      zero_prop = tissue_zero_prop,
      expressed_genes = tissue_expressed_genes
    )
    
    message(sprintf("  %-20s: %5d cells, %.1f%% 零值, %5d 表达基因", 
                    tissue, 
                    length(tissue_cells),
                    tissue_zero_prop * 100,
                    tissue_expressed_genes))
  }
  
  # ========================================================================
  # 4. 按组织×细胞类型交叉诊断
  # ========================================================================
  message("\n【4】组织×细胞类型交叉诊断（仅显示异常组合）")
  message("─────────────────────────────────────")
  
  cross_stats <- data.frame()
  
  for(tissue in tissues) {
    for(ct in cell_types) {
      combo_cells <- rownames(seurat_obj@meta.data)[
        seurat_obj@meta.data[[tissue_col]] == tissue & 
          seurat_obj@meta.data[[celltype_col]] == ct
      ]
      
      if(length(combo_cells) > 0) {
        combo_counts <- counts_matrix[, combo_cells, drop = FALSE]
        combo_zero_prop <- sum(combo_counts == 0) / length(combo_counts)
        combo_expressed <- sum(Matrix::rowSums(combo_counts > 0) > 0)
        
        cross_stats <- rbind(cross_stats, data.frame(
          Tissue = tissue,
          CellType = ct,
          N_Cells = length(combo_cells),
          Zero_Prop = combo_zero_prop,
          Expressed_Genes = combo_expressed
        ))
        
        # 只显示零值比例异常高的组合
        if(combo_zero_prop > 0.95) {
          message(sprintf("  ⚠️  %s × %s: %.1f%% 零值 (%d 细胞)", 
                          tissue, ct, combo_zero_prop * 100, length(combo_cells)))
        }
      }
    }
  }
  
  # ========================================================================
  # 5. 推荐过滤阈值
  # ========================================================================
  message("\n【5】推荐的基因过滤阈值")
  message("─────────────────────────────────────")
  
  # 方法1：基于分位数
  q05_cells <- quantile(gene_stats$nonzero_cells, 0.05)
  q05_count <- quantile(gene_stats$total_count, 0.05)
  q05_prop <- quantile(gene_stats$nonzero_prop, 0.05)
  
  # 方法2：基于实际数据分布
  # 找到表达频率的"断崖"位置
  prop_breaks <- seq(0, 0.5, by = 0.05)
  prop_hist <- hist(gene_stats$nonzero_prop[gene_stats$nonzero_prop < 0.5], 
                    breaks = prop_breaks, plot = FALSE)
  
  # 找到第一个大的跳跃（可能是真实基因和填充基因的分界）
  prop_diffs <- diff(prop_hist$counts)
  suggested_cutoff_idx <- which(prop_diffs > median(prop_diffs) * 2)[1]
  
  if(!is.na(suggested_cutoff_idx)) {
    suggested_prop <- prop_breaks[suggested_cutoff_idx + 1]
  } else {
    suggested_prop <- 0.1  # 默认值
  }
  
  message("\n策略A: 保守过滤（保留95%基因）")
  message(sprintf("  • 最小表达细胞数: %.0f (%.1f%%)", 
                  q05_cells, q05_cells/total_cells * 100))
  message(sprintf("  • 最小总表达量: %.0f", q05_count))
  message(sprintf("  • 最小表达比例: %.1f%%", q05_prop * 100))
  
  message("\n策略B: 标准过滤（移除可能的填充基因）")
  message(sprintf("  • 最小表达细胞数: %.0f (%.1f%%)", 
                  total_cells * suggested_prop, suggested_prop * 100))
  message(sprintf("  • 最小总表达量: 50"))
  message(sprintf("  • 最小表达比例: %.1f%%", suggested_prop * 100))
  
  message("\n策略C: 严格过滤（仅保留高质量基因）")
  message(sprintf("  • 最小表达细胞数: %.0f (10%%)", total_cells * 0.1))
  message(sprintf("  • 最小总表达量: 100"))
  message(sprintf("  • 最小表达比例: 10%%"))
  
  # ========================================================================
  # 6. 可视化诊断结果
  # ========================================================================
  message("\n【6】生成诊断可视化图表...")
  message("─────────────────────────────────────")
  
  tryCatch({
    pdf("gene_quality_diagnosis.pdf", width = 14, height = 10)
    
    par(mfrow = c(2, 3), mar = c(5, 4, 4, 2))
    
    # 图1: 基因表达频率分布
    hist(gene_stats$nonzero_prop, breaks = 50, 
         col = "steelblue", border = "white",
         main = "Gene Expression Frequency Distribution",
         xlab = "Proportion of Cells Expressing Gene",
         ylab = "Number of Genes")
    abline(v = suggested_prop, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = "Suggested cutoff", 
           col = "red", lty = 2, lwd = 2, bty = "n")
    
    # 图2: 表达量分布（log scale）
    hist(log10(gene_stats$total_count + 1), breaks = 50,
         col = "coral", border = "white",
         main = "Total Expression Distribution",
         xlab = "log10(Total Count + 1)",
         ylab = "Number of Genes")
    
    # 图3: 散点图（表达频率 vs 平均表达量）
    plot(gene_stats$nonzero_prop, log10(gene_stats$mean_expr + 1),
         pch = 16, col = adjustcolor("blue", alpha = 0.3),
         main = "Expression Frequency vs Mean Expression",
         xlab = "Proportion of Expressing Cells",
         ylab = "log10(Mean Expression + 1)")
    abline(v = suggested_prop, col = "red", lwd = 2, lty = 2)
    abline(h = log10(1), col = "red", lwd = 2, lty = 2)
    
    # 图4: 细胞类型零值比例
    ct_zeros <- sapply(celltype_stats, function(x) x$zero_prop * 100)
    barplot(sort(ct_zeros), las = 2, col = "lightblue",
            main = "Zero Proportion by Cell Type",
            ylab = "Zero Proportion (%)",
            cex.names = 0.7)
    
    # 图5: 组织零值比例
    tissue_zeros <- sapply(tissue_stats, function(x) x$zero_prop * 100)
    barplot(sort(tissue_zeros), las = 2, col = "lightcoral",
            main = "Zero Proportion by Tissue",
            ylab = "Zero Proportion (%)",
            cex.names = 0.7)
    
    # 图6: 表达基因数量比较
    ct_expressed <- sapply(celltype_stats, function(x) x$expressed_genes)
    tissue_expressed <- sapply(tissue_stats, function(x) x$expressed_genes)
    
    boxplot(list(
      "By Cell Type" = ct_expressed,
      "By Tissue" = tissue_expressed
    ), col = c("lightblue", "lightcoral"),
    main = "Number of Expressed Genes",
    ylab = "Number of Genes")
    
    dev.off()
    
    message("✓ 诊断图表已保存: gene_quality_diagnosis.pdf")
  }, error = function(e) {
    if(dev.cur() > 1) dev.off()
    message("✗ 可视化生成失败: ", e$message)
  })
  
  # ========================================================================
  # 7. 保存诊断结果
  # ========================================================================
  message("\n【7】保存诊断数据文件...")
  message("─────────────────────────────────────")
  
  # 保存基因统计
  write.csv(gene_stats, "gene_quality_stats.csv", row.names = FALSE)
  message("✓ 基因统计已保存: gene_quality_stats.csv")
  
  # 保存可能的填充基因列表
  if(nrow(likely_padded) > 0) {
    write.csv(likely_padded, "likely_padded_genes.csv", row.names = FALSE)
    message("✓ 可能填充基因已保存: likely_padded_genes.csv")
  }
  
  # 保存交叉统计
  write.csv(cross_stats, "tissue_celltype_cross_stats.csv", row.names = FALSE)
  message("✓ 交叉统计已保存: tissue_celltype_cross_stats.csv")
  
  # ========================================================================
  # 8. 返回诊断结果
  # ========================================================================
  message("\n╔═══════════════════════════════════════════════════════════╗")
  message("║     诊断完成 Diagnosis Complete                            ║")
  message("╚═══════════════════════════════════════════════════════════╝\n")
  
  results <- list(
    gene_stats = gene_stats,
    likely_padded = likely_padded,
    celltype_stats = celltype_stats,
    tissue_stats = tissue_stats,
    cross_stats = cross_stats,
    thresholds = list(
      conservative = list(
        min_cells = q05_cells,
        min_count = q05_count,
        min_prop = q05_prop
      ),
      standard = list(
        min_cells = total_cells * suggested_prop,
        min_count = 50,
        min_prop = suggested_prop
      ),
      strict = list(
        min_cells = total_cells * 0.1,
        min_count = 100,
        min_prop = 0.1
      )
    ),
    summary = list(
      total_genes = total_genes,
      total_cells = total_cells,
      zero_proportion = zero_prop,
      n_likely_padded = nrow(likely_padded)
    )
  )
  
  return(invisible(results))
}

# ============================================================================
# 便捷包装函数：直接应用推荐过滤
# ============================================================================
filter_seurat_genes <- function(
    seurat_obj,
    diagnosis_result = NULL,
    strategy = c("conservative", "standard", "strict"),
    custom_thresholds = NULL
) {
  
  strategy <- match.arg(strategy)
  
  # 如果没有提供诊断结果，先运行诊断
  if(is.null(diagnosis_result)) {
    message("未提供诊断结果，先运行诊断...")
    diagnosis_result <- diagnose_seurat_gene_quality(seurat_obj)
  }
  
  # 确定阈值
  if(!is.null(custom_thresholds)) {
    thresholds <- custom_thresholds
  } else {
    thresholds <- diagnosis_result$thresholds[[strategy]]
  }
  
  message(sprintf("\n应用过滤策略: %s", strategy))
  message(sprintf("  • 最小表达细胞数: %.0f", thresholds$min_cells))
  message(sprintf("  • 最小总表达量: %.0f", thresholds$min_count))
  message(sprintf("  • 最小表达比例: %.3f", thresholds$min_prop))
  
  # 获取基因统计
  gene_stats <- diagnosis_result$gene_stats
  
  # 应用过滤
  genes_to_keep <- (gene_stats$nonzero_cells >= thresholds$min_cells) &
    (gene_stats$total_count >= thresholds$min_count) &
    (gene_stats$nonzero_prop >= thresholds$min_prop)
  
  message(sprintf("\n过滤前: %d 基因", nrow(seurat_obj)))
  
  # 过滤Seurat对象
  seurat_filtered <- seurat_obj[genes_to_keep, ]
  
  message(sprintf("过滤后: %d 基因", nrow(seurat_filtered)))
  message(sprintf("移除: %d 基因 (%.1f%%)\n", 
                  sum(!genes_to_keep), 
                  sum(!genes_to_keep)/length(genes_to_keep) * 100))
  
  return(seurat_filtered)
}

# ============================================================================
# 完整使用流程
# ============================================================================

# 1. 运行诊断
diagnosis <- diagnose_seurat_gene_quality(
  seurat_obj = seurat_obj_main,
  celltype_col = "Annotation",
  tissue_col = "tissue",
  sample_col = "sample"
)

# 2. 查看诊断结果
# 查看可能的填充基因
head(diagnosis$likely_padded, 20)

# 查看推荐阈值
diagnosis$thresholds

# # 3. 应用过滤（三种策略可选）
# # 策略A: 保守（适合探索性分析）
# seurat_obj_main_filtered_conservative <- filter_seurat_genes(
#   seurat_obj_main, 
#   diagnosis, 
#   strategy = "conservative"
# )
# 
# # 策略B: 标准（推荐用于Pseudobulk分析）
# seurat_obj_main_filtered_standard <- filter_seurat_genes(
#   seurat_obj_main, 
#   diagnosis, 
#   strategy = "standard"
# )

# 策略C: 严格（用于高质量分析）
seurat_obj_main_filtered_strict <- filter_seurat_genes(
  seurat_obj_main, 
  diagnosis, 
  strategy = "strict"
)

# 4. 或者使用自定义阈值
seurat_obj_main_filtered_custom <- filter_seurat_genes(
  seurat_obj_main,
  diagnosis,
  custom_thresholds = list(
    min_cells = 100,
    min_count = 50,
    min_prop = 0.05
  )
)

# 5. 使用过滤后的对象进行后续分析
# 例如Pseudobulk分析
av_filtered <- AggregateExpression(
  seurat_obj_main_filtered_standard,
  group.by = c("tissue", "sample", "Annotation_2"),
  assays = "RNA",
  slot = "counts"
)

######################################################################################################
######################################################################################################
######################################################################################################

# ========================================
# 分泌型上皮细胞标记 (Secretory Epithelial Cell Marker)
# ========================================

# 提取目标基因表达数据
expr <- FetchData(seurat_obj_main, 
                  vars = c("PRR4", "PIP", "SCGB1A1", "SCGB3A2", 
                           "MUC5B", "BPIFB1", "BPIFA1", "STATH", 
                           "ZG16B", "EPCAM"))

# 初始化标记列
seurat_obj_main$secretory_epi_high <- 0

# 定义三个判定条件 (Criteria)
# 条件1: 10个基因中至少2个高表达 (≥2 genes with expression > 1)
condition1 <- rowSums(expr[, 1:10] > 1) >= 2

# 条件2: STATH基因高表达 (STATH expression > 1)
condition2 <- expr[, 8] > 1  # STATH是第8列

# 条件3: SCGB3A2基因高表达 (SCGB3A2 expression > 1)
condition3 <- expr[, 4] > 1  # SCGB3A2是第4列

# 复合标记: 满足任一条件即标记为1 (OR logic)
seurat_obj_main$secretory_epi_high <- as.numeric(
  condition1 | condition2 | condition3
)

# 统计非上皮细胞中的阳性标记分布
# (Distribution of positive cells in non-epithelial populations)
table(seurat_obj_main$Annotation[seurat_obj_main$secretory_epi_high == 1 & 
                                   !(seurat_obj_main$Annotation %in% c('Epithelial', 'Myeloid'))], 
      seurat_obj_main$sample[seurat_obj_main$secretory_epi_high == 1 & 
                               !(seurat_obj_main$Annotation %in% c('Epithelial', 'Myeloid'))])

seurat_obj_main <- subset(
  seurat_obj_main,
  subset = secretory_epi_high == 1 & !(seurat_obj_main$Annotation %in% c('Epithelial', 'Myeloid')),
  invert = TRUE
)

# 综合T细胞标记（包含核心+效应标记）
expr <- FetchData(seurat_obj_main, vars = c("CD3E", "CD247", "TRBC1", "ZAP70",  # 核心T细胞
                                            "CD8A", "GZMA", "PRF1", "GNLY",      # 细胞毒性
                                            "IFNG", "PTPRC"))                     # 效应+泛白
seurat_obj_main$T_NK_high <- 0
seurat_obj_main$T_NK_high <- as.numeric(
  (rowSums(expr[,1:4] > 1) >= 2 | rowSums(expr[,5:9] > 1) >= 3) &  # 核心T标记≥2个 OR 效应标记≥3个
    expr$PTPRC > 1                                                      # 且CD45+
)

# 过滤误注释细胞
seurat_obj_main <- subset(
  seurat_obj_main,
  subset = T_NK_high == 1 & !(seurat_obj_main$Annotation %in% c('T', 'B', 'Myeloid')),
  invert = TRUE
)

# ========================================
# B细胞标记污染检测 (B Cell Marker Contamination Detection)
# ========================================

# 定义B细胞核心标记基因（按特异性分层）
# Tier 1: 高度特异的B细胞标记
b_core_markers <- c(
  "MS4A1",     # CD20 - 成熟B细胞经典标记
  "CD79A",     # B细胞受体信号组件
  "CD79B",     # B细胞受体信号组件
  "CD19"       # B细胞泛标记
)

# Tier 2: 浆细胞特异标记
plasma_markers <- c(
  "JCHAIN",    # 连接链 - 浆细胞极高特异性
  "MZB1",      # 浆细胞特征基因
  "SDC1"       # CD138 - 浆细胞表面标记
)

# Tier 3: 免疫球蛋白基因（谨慎使用，可能有低水平背景表达）
ig_markers <- c(
  "IGHG1",     # IgG1重链
  "IGHG2",     # IgG2重链  
  "IGHG4",     # IgG4重链
  "IGHA1",     # IgA1重链
  "IGLC1",     # Lambda轻链
  "IGLC2",     # Lambda轻链
  "IGKC"       # Kappa轻链
)

# 合并所有B细胞相关标记
all_b_markers <- c(b_core_markers, plasma_markers, ig_markers)

# 检查哪些基因在数据中存在
available_b_markers <- all_b_markers[all_b_markers %in% rownames(seurat_obj_main)]
cat("可用的B细胞标记基因:", length(available_b_markers), "个\n")
print(available_b_markers)

# 提取B细胞标记基因的表达数据
expr <- FetchData(seurat_obj_main, vars = available_b_markers)

# 分层检测策略
# 检查各层标记的可用性
available_core <- b_core_markers[b_core_markers %in% available_b_markers]
available_plasma <- plasma_markers[plasma_markers %in% available_b_markers]
available_ig <- ig_markers[ig_markers %in% available_b_markers]

cat("\n各层标记基因可用情况:\n")
cat("核心B细胞标记:", length(available_core), "个\n")
cat("浆细胞标记:", length(available_plasma), "个\n")
cat("免疫球蛋白标记:", length(available_ig), "个\n")

# 定义表达阈值（根据您的数据调整）
expression_threshold <- 1  # log-normalized expression > 1

# 条件1: 核心B细胞标记 - 至少2个高表达
if(length(available_core) >= 2) {
  condition1 <- rowSums(expr[, available_core, drop = FALSE] > expression_threshold) >= 2
} else if(length(available_core) == 1) {
  condition1 <- expr[, available_core] > expression_threshold
} else {
  condition1 <- rep(FALSE, nrow(expr))
}

# 条件2: 浆细胞标记 - 至少1个高表达（JCHAIN权重最高）
if(length(available_plasma) > 0) {
  # 如果有JCHAIN，单独判断
  if("JCHAIN" %in% available_plasma) {
    condition2 <- expr[, "JCHAIN"] > expression_threshold
  } else {
    condition2 <- rowSums(expr[, available_plasma, drop = FALSE] > expression_threshold) >= 1
  }
} else {
  condition2 <- rep(FALSE, nrow(expr))
}

# 条件3: 免疫球蛋白基因 - 至少3个高表达（更严格，避免假阳性）
if(length(available_ig) >= 3) {
  condition3 <- rowSums(expr[, available_ig, drop = FALSE] > expression_threshold) >= 3
} else if(length(available_ig) >= 1) {
  condition3 <- rowSums(expr[, available_ig, drop = FALSE] > expression_threshold) >= length(available_ig)
} else {
  condition3 <- rep(FALSE, nrow(expr))
}

# 初始化标记列
seurat_obj_main$b_marker_positive <- 0

# 复合判断: 满足任一条件即标记
seurat_obj_main$b_marker_positive <- as.numeric(
  condition1 | condition2 | condition3
)

cat("\n标记为B细胞marker阳性的细胞数:", 
    sum(seurat_obj_main$b_marker_positive), "\n")

table(seurat_obj_main$Annotation[seurat_obj_main$b_marker_positive == 1 & 
                                   !(seurat_obj_main$Annotation %in% c('B'))], 
      seurat_obj_main$sample[seurat_obj_main$b_marker_positive == 1 & 
                               !(seurat_obj_main$Annotation %in% c('B'))])

seurat_obj_main <- subset(
  seurat_obj_main,
  subset = b_marker_positive == 1 & !(seurat_obj_main$Annotation %in% c('B')),
  invert = TRUE
)

table(seurat_obj_main$tissue)
######################################################################################################
######################################################################################################
######################################################################################################
