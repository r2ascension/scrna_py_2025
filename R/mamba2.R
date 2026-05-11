if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
rm(list = ls());gc()#删除所有环境变量整理内存空间
setwd('E:/R/5/')
download.file('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz','my.pbmc.gz')
untar(gunzip("my.pbmc.gz"))
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
## Warning: Feature names cannot have underscores ('_'), replacing with dashes
## ('-')
#初步查看Seurat对象
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
ncol(pbmc)
## [1] 2700
ncol(pbmc.data)
## [1] 2700

lalala <- as.data.frame(pbmc@assays[["RNA"]]@layers[["counts"]])
write.table(lalala,'mycount.txt',sep = '\t')#表达矩阵可以这么存出来

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
#人源的数据为MT，鼠源的需要换成mt
head(pbmc@meta.data,5)
##                  orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759
## AAACATTGAGCTAC-1     pbmc3k       4903         1352  3.7935958
## AAACATTGATCAGC-1     pbmc3k       3147         1129  0.8897363
## AAACCGTGCTTCCG-1     pbmc3k       2639          960  1.7430845
## AAACCGTGTATGCG-1     pbmc3k        980          521  1.2244898
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
CombinePlots(plots = list(plot1, plot2)) 
## Warning: CombinePlots is being deprecated. Plots should now be combined using
## the patchwork system.

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
#CombinePlots这步需要你的绘图窗口足够大
CombinePlots(plots = list(plot1, plot2)) 
## Warning: CombinePlots is being deprecated. Plots should now be combined using
## the patchwork system.

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
ncol(as.data.frame(pbmc@assays[["RNA"]]@layers[["counts"]]))
## [1] 2638

#将表达数据标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#pbmc[["RNA"]]@data

#寻找高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(pbmc), 10)#查看十个程度最强的高变基因

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
## When using repel, set xnudge and ynudge to 0 for optimal results
plot1 + plot2# plot variable features with and without labels
## Warning: Transformation introduced infinite values in continuous x-axis
## Transformation introduced infinite values in continuous x-axis

#scale全部基因
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
## Centering and scaling data matrix
#pbmc <- ScaleData(pbmc) ##默认用高变基因来scale

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
## PC_ 1 
## Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
##     FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
##     PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
## Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
##     CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
##     MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
## PC_ 2 
## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
##     HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
##     BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
## Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
##     CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
##     TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC 
## PC_ 3 
## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA 
##     HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8 
##     PLAC8, BLNK, MALAT1, SMIM14, PLD4, LAT2, IGLL5, P2RX5, SWAP70, FCGR2B 
## Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU 
##     HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1 
##     NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA 
## PC_ 4 
## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HLA-DPB1, HIST1H2AC, PF4, TCL1A 
##     SDPR, HLA-DPA1, HLA-DRB1, HLA-DQA2, HLA-DRA, PPBP, LINC00926, GNG11, HLA-DRB5, SPARC 
##     GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, GZMB, CLU, TUBB1 
## Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL 
##     AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7 
##     LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6 
## PC_ 5 
## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2 
##     GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, CCL5 
##     RBP7, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX 
## Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1 
##     PTGES3, LILRB2, MAL, CD27, HN1, CD2, GDI2, ANXA5, CORO1B, TUBA1B 
##     FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
## PC_ 1 
## Positive:  CST3, TYROBP, LST1, AIF1, FTL 
## Negative:  MALAT1, LTB, IL32, IL7R, CD2 
## PC_ 2 
## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 
## Negative:  NKG7, PRF1, CST7, GZMB, GZMA 
## PC_ 3 
## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1 
## Negative:  PPBP, PF4, SDPR, SPARC, GNG11 
## PC_ 4 
## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 
## Negative:  VIM, IL7R, S100A6, IL32, S100A8 
## PC_ 5 
## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY 
## Negative:  LTB, IL7R, CKB, VIM, MS4A7
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
## Warning: Removed 32525 rows containing missing values (`geom_point()`).

ElbowPlot(pbmc) 


pbmc <- FindNeighbors(pbmc, dims = 1:10)
## Computing nearest neighbor graph
## Computing SNN
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2638
## Number of edges: 95965
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8723
## Number of communities: 9
## Elapsed time: 0 seconds

pbmc <- RunUMAP(pbmc, dims = 1:10)
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
## 14:57:59 UMAP embedding parameters a = 0.9922 b = 1.112
## 14:57:59 Read 2638 rows and found 10 numeric columns
## 14:57:59 Using Annoy for neighbor search, n_neighbors = 30
## 14:57:59 Building Annoy index with metric = cosine, n_trees = 50
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## 14:57:59 Writing NN index file to temp file C:\Users\CPU1302\AppData\Local\Temp\RtmpmSipEx\file5b5c3ac8e77
## 14:57:59 Searching Annoy index using 1 thread, search_k = 3000
## 14:57:59 Annoy recall = 100%
## 14:58:00 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
## 14:58:00 Initializing from normalized Laplacian + noise (using irlba)
## 14:58:00 Commencing optimization for 500 epochs, with 105124 positive edges
## 14:58:06 Optimization finished
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
## Calculating cluster 0
## Calculating cluster 1
## Calculating cluster 2
## Calculating cluster 3
## Calculating cluster 4
## Calculating cluster 5
## Calculating cluster 6
## Calculating cluster 7
## Calculating cluster 8

if(!require(dplyr))install.packages("dplyr")
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
## # A tibble: 18 × 7
## # Groups:   cluster [9]
##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    
##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
##  1 1.74e-109       1.07 0.897 0.593 2.39e-105 0       LDHB    
##  2 1.17e- 83       1.33 0.435 0.108 1.60e- 79 0       CCR7    
##  3 0               5.57 0.996 0.215 0         1       S100A9  
##  4 0               5.48 0.975 0.121 0         1       S100A8  
##  5 7.99e- 87       1.28 0.981 0.644 1.10e- 82 2       LTB     
##  6 2.61e- 59       1.24 0.424 0.111 3.58e- 55 2       AQP3    
##  7 0               4.31 0.936 0.041 0         3       CD79A   
##  8 9.48e-271       3.59 0.622 0.022 1.30e-266 3       TCL1A   
##  9 1.17e-178       2.97 0.957 0.241 1.60e-174 4       CCL5    
## 10 4.93e-169       3.01 0.595 0.056 6.76e-165 4       GZMK    
## 11 3.51e-184       3.31 0.975 0.134 4.82e-180 5       FCGR3A  
## 12 2.03e-125       3.09 1     0.315 2.78e-121 5       LST1    
## 13 1.05e-265       4.89 0.986 0.071 1.44e-261 6       GZMB    
## 14 6.82e-175       4.92 0.958 0.135 9.36e-171 6       GNLY    
## 15 1.48e-220       3.87 0.812 0.011 2.03e-216 7       FCER1A  
## 16 1.67e- 21       2.87 1     0.513 2.28e- 17 7       HLA-DPB1
## 17 7.73e-200       7.24 1     0.01  1.06e-195 8       PF4     
## 18 3.68e-110       8.58 1     0.024 5.05e-106 8       PPBP
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc,'pbmc.rds')
pbmc<- readRDS('pbmc.rds')
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
