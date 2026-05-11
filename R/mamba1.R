Sys.setenv(LANGUAGE = "en")
library(BiocManager)
library(multtest)
library(dplyr)
library(devtools)
library(Seurat)
#devtools::install_github('pzhaonet/mindr')
library('mindr')
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))
if(!require(mindr))install.packages("tidyverse")

#####自动读取cellranger(LINUX)输出的feature barcode matric
setwd('E:/R/4')
rm(list = ls())
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/") 
#自动读取10X的数据，是一些tsv与mtx文件
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "biomamba")

####仅有一个稀疏矩阵时的读取方法#####
matrix_data <- read.table("single_cell_datamatrix.txt", sep="\t", header=T, row.names=1)
dim(matrix_data)
## [1] 13714  2700
seurat_obj <- CreateSeuratObject(counts = matrix_data)

######读取RDS文件############
rm(list = ls())
pbmc <- readRDS("panc8.rds")
saveRDS(pbmc,"pbmc.rds")
out2 <- str(pbmc)  %>%  
  capture.output(.) %>% 
  gsub(
    pattern = "\\.\\. ",
    replacement = "#",
    x = .
  ) %>% 
  gsub(
    pattern = "\\.\\.@",
    replacement = "# ",
    x = .
  ) %>% 
  gsub(pattern="^\\s+#",replace="#")


mindr::mm(
  from = out2,
  input_type = "markdown",
  output_type = "widget",
  root = "Seurat"
)