# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
# setwd("/home/h2048/data/R/1218")
library(reticulate)
library(harmony)
library(ggplot2)
library(patchwork)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
setwd('/home/h2048/data/R/1215/merge/cleaned_samples_COMPLETE_v4.1/')
seurat_obj <- readRDS(
  '/home/h2048/data/R/1215/merge/cleaned_samples_COMPLETE_v4.1/merged_seurat_standardized.rds'
)
doublet_cols <- grep(
  "^(pANN_|DF\\.classifications_)",
  colnames(seurat_obj@meta.data),
  value = TRUE
)

cat(sprintf("Found %d DoubletFinder columns to remove\n", length(doublet_cols)))
if (length(doublet_cols) > 0) {
  seurat_obj@meta.data <- seurat_obj@meta.data[,
    !colnames(seurat_obj@meta.data) %in% doublet_cols
  ]
  cat("✓ Columns removed\n")
} else {
  cat("No DoubletFinder columns found\n")
}
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
seurat_obj <- subset(
  seurat_obj,
  subset = GEO %in% c('GSE210695'),
  invert = TRUE
)
print(str(seurat_obj))
head(seurat_obj)
saveRDS(seurat_obj, 'merged_seurat_standardized_filterd.rds')
GetH5ad(seurat_obj, 'merged_seurat_standardized_filterd.h5ad')
