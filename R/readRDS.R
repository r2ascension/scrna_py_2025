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
setwd("E:/R/1012/")
source("E:/R/enrichment_functions.R")
source("E:/R/scMASC.R")

seurat_obj <- readRDS("E:/R/Source/final/Ordovas_Montanes_2018.rds")
table(seurat_obj@meta.data[["sample"]])
table(seurat_obj@meta.data[["tissue"]])
table(seurat_obj@meta.data[["tissue_sampling_method"]])
table(seurat_obj@meta.data[["dataset"]])
table(seurat_obj@meta.data[["batch"]])
# seurat_obj <- seurat_obj_main
# rm(seurat_obj_main)
# 1. 查看现有的降维结果
names(seurat_obj@reductions)
# 例如: "pca" "umap" "tsne" "harmony"
seurat_obj$tissue <- 'nose'
# 2. 删除特定的降维
seurat_obj[["pca"]] <- NULL
seurat_obj[["umap"]] <- NULL
seurat_obj@graphs[["integrated_snn"]]<- NULL
seurat_obj@assays$RNA@scale.data <- NULL
seurat_obj@assays$SCT <- NULL
seurat_obj@assays[["RNA"]]@layers[["scale.data"]] <- NULL
seurat_obj@assays$integrated <- NULL
# integrated
# 3. 删除所有降维（批量）
for(i in names(seurat_obj@reductions)){
  seurat_obj[[i]] <- NULL
}
# seurat_obj <- subset(seurat_obj, tissue=='nose' & tissue_sampling_method== 'brush')
seurat_obj <- subset(seurat_obj, subset = group %in% c('Control') )
# 4. 验证
names(seurat_obj@reductions)
seurat_obj$dataset <- 'Maya_E_Kotas_2022'
seurat_obj$study <- 'Maya_E_Kotas_2022'
seurat_obj$tissue <- 'sinus'
seurat_obj$tissue_sampling_method <- 'brush'
dir <- "E:/R/1011"
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
setwd(dir)
saveRDS(seurat_obj,'Ordovas_Montanes_2018.rds')
