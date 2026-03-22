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
  '/home/h2048/data/R/1215/merge/cleaned_samples_COMPLETE_v4.1/merged_seurat_standardized_filterd.rds'
)
# ===== Convert cellLabel to character before export =====
# 将cellLabel转换为字符串类型

# If cellLabel is a factor
if (is.factor(seurat_obj$cellLabel)) {
  seurat_obj$cellLabel <- as.character(seurat_obj$cellLabel)
  cat("✓ Converted cellLabel from factor to character\n")
}

# Handle NA values (replace with "Unknown")
seurat_obj$cellLabel[is.na(seurat_obj$cellLabel)] <- "NA"

# Verify conversion
cat("cellLabel class:", class(seurat_obj$cellLabel), "\n")
cat("Sample values:\n")
print(head(seurat_obj$cellLabel))

# Now try exporting again
GetH5ad(seurat_obj, 'merged_seurat_standardized_filterd.h5ad')
