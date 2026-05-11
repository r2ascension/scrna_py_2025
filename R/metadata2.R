# Load required libraries
library(Seurat)
library(data.table)
library(Matrix)
setwd("/home/h2048/data/R/1029")
# ============================================================================
# Configuration
# ============================================================================
count_file <- "/home/h2048/data/source/GSE176269/GSE176269_CovidStudy_rawCounts_061721.txt" # Update with your count file path
metadata_file <- "/home/h2048/data/source/GSE176269/GSE176269_CovidStudy_phenotype_061721.txt" # Update with your metadata file path
output_file <- "Jennifer_P_Wang_2021.rds" # Output file name
project_name <- "MyProject" # Project name for Seurat object

# ============================================================================
# Step 1: Read count matrix
# ============================================================================
cat("Reading count matrix...\n")
# Use fread for efficient reading of large files
count_data <- fread(count_file, header = TRUE, data.table = FALSE)

# Extract gene names (first column) and convert to matrix
gene_names <- count_data[, 1]
count_matrix <- as.matrix(count_data[, -1])
rownames(count_matrix) <- gene_names

# Convert to sparse matrix for memory efficiency
count_matrix <- as(count_matrix, "sparseMatrix")

cat(
  "Count matrix dimensions: %d genes x %d cells\n",
  nrow(count_matrix),
  ncol(count_matrix)
)

# ============================================================================
# Step 2: Read metadata
# ============================================================================
cat("Reading metadata...\n")
metadata <- fread(metadata_file, header = TRUE, data.table = FALSE)

# Set cell barcodes as rownames
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1, drop = FALSE]

cat(
  "Metadata dimensions: %d cells x %d features\n",
  nrow(metadata),
  ncol(metadata)
)

# ============================================================================
# Step 3: Quality check - match cells between count and metadata
# ============================================================================
cat("Checking cell matching...\n")
common_cells <- intersect(colnames(count_matrix), rownames(metadata))
cat("Common cells: %d\n", length(common_cells))

# Subset to common cells
count_matrix <- count_matrix[, common_cells]
metadata <- metadata[common_cells, , drop = FALSE]

# ============================================================================
# Step 4: Create Seurat object
# ============================================================================
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = count_matrix,
  project = project_name,
  meta.data = metadata,
  min.cells = 3, # Keep genes expressed in >= 3 cells
  min.features = 200 # Keep cells with >= 200 genes
)

# Add basic QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^RP[SL]"
)

# Print basic statistics
cat("\n=== Seurat Object Summary ===\n")
cat("Cells: %d\n", ncol(seurat_obj))
cat("Genes: %d\n", nrow(seurat_obj))
cat("Metadata columns: %d\n", ncol(seurat_obj@meta.data))
print(head(seurat_obj@meta.data, 3))

table(seurat_obj$status)
seurat_obj <- subset(seurat_obj, subset = status == 'healthy')
seurat_obj$tissue <- 'nose'
seurat_obj$tissue_sampling_method <- 'fluid'
seurat_obj$dataset <- 'Jennifer_P_Wang_2021'
seurat_objsample <- seurat_obj$sampID
# ============================================================================
# Step 5: Save Seurat object
# ============================================================================
cat("\nSaving Seurat object to %s...\n", output_file)
saveRDS(seurat_obj, file = output_file)
cat("Seurat object saved successfully!\n")

# ============================================================================
# Optional: Quick QC plots
# ============================================================================
# Uncomment to generate QC plots
# library(ggplot2)
# pdf("seurat_qc_plots.pdf", width = 12, height = 4)
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# dev.off()
