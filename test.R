library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(celda)
library(viridis)
getwd()
seurat_obj_1 <- readRDS(
    '/home/h2048/data/source/final/Kerstin_B_Meyer_2021.rds'
)
head(seurat_obj_1)
seurat_obj_1$sample <- seurat_obj_1$sample_id
saveRDS(seurat_obj_1, '/home/h2048/data/source/final/Kerstin_B_Meyer_2021.rds')
