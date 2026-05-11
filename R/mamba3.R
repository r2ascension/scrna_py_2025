devtools::install_github('satijalab/seurat-data')
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(tidyverse)

# 如果你的内存不是很充裕，可以进行抽样，降低计算量
# pbmc <- subset(pbmc, downsample = 50)
ifnb <- readRDS('./第四讲.多样本整合/pbmcrenamed.rds')
# 这其实是一个已经整合好的Seurat对象，不要问我哪里来的，看下去，后面你自然会明白
class(ifnb)

# 其中包含了两个分组：
unique(ifnb$group)

ifnb.list <- SplitObject(ifnb, split.by = "group")

class(ifnb.list)

# 分别取出这两个对象
C57 <- ifnb.list$C57
AS1 <- ifnb.list$AS1
# 简单利用merge整合这两个数据
pbmc <- merge(C57, y = c(AS1), add.cell.ids = c("C57", "AS1"), project = "ALL")
# 这时我们就获得了一个新的Seurat对象
pbmc
# 可以发现这是barcode的名称发生了变化：
head(colnames(pbmc))
unique(sapply(X = strsplit(colnames(pbmc), split = "_"), FUN = "[", 1))

table(pbmc$orig.ident)

### 定义一个函数帮助预处理数据
myfunction1 <- function(testA.seu){
  testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)
  return(testA.seu)
}

# 执行对每个数据的NormalizeData（标准化表达矩阵）与FindVariableFeatures（高变基因计算），这两部操作的过程可参考第三讲：https://www.bilibili.com/video/BV1S44y1b76Z?p=4&vd_source=6335356a0d3631ad476b7c7de83892db
# 
C57 <- myfunction1(C57)
AS1 <- myfunction1(AS1)

testAB.anchors <- FindIntegrationAnchors(object.list = list(C57,AS1), dims = 1:20)

testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)

# 将后续计算用默认矩阵由"RNA"改为"integrated"
DefaultAssay(testAB.integrated) <- "integrated"
# 整合后的数据从scale开始运行，后续基本与单样本分析部分无异
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))

testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)

testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.5)

testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)


test.seu <- pbmc
# harmony前需要完成标准化、高变基因计算、scale、PCA等分析
test.seu <-  test.seu%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()

test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)

#按照我们上述的group变量进行Harmony操作
test.seu=test.seu %>% RunHarmony("group", plot_convergence = TRUE)

# 以Harmony的结果进行后续的降维与分群
# UMAP及分群：
test.seu <- test.seu %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

# TSNE
test.seu <- test.seu %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)

# 创建图片对象
p3 <- DimPlot(test.seu, reduction = "tsne", group.by = "group", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p4 <- DimPlot(test.seu, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

# 看一下降维结果展示
p3|p4