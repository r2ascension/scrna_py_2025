# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
output_dir <- '/home/h2048/data/R/0118'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
library(reticulate)
library(harmony)
library(ggplot2)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct

seurat_obj <- GetSeurat(
  '/home/h2048/data/py/0115/epithelial_DUAL_SCANVI_v2_7_PRODUCTION/checkpoints/adata_preprocessed.h5ad'
)
table(seurat_obj$disease_level_1, seurat_obj$sample_id)


seurat_obj$sample <- seurat_obj$orig.ident
