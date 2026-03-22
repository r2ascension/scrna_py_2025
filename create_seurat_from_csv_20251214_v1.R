#!/usr/bin/env Rscript
# Version 2: For multi-sample data where count matrix needs sample prefixes
# Assumes count matrix and metadata have matching cell order
# Automatically adds "Sample_Barcode" format to ensure uniqueness

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})
setwd("/home/h2048/data/R/1214")
# ===== Configuration Section =====
metadata_path <- "/home/h2048/data/source/1213/GSE210695_SU2C_nasal_single_cell_labels_fine.csv"
count_matrix_path <- "/home/h2048/data/source/1213/GSE210695_SU2C_nasal_single_cell_counts.tsv" # UPDATE THIS
output_rds_path <- "/home/h2048/data/R/1214/GSE210695.rds"

sample_column <- "sample"
celltype_column <- "cellLabel"
barcode_column <- "barcode"

# Set to TRUE if count matrix and metadata are in EXACTLY the same order
assume_same_order <- TRUE # CHANGE to FALSE if unsure

# Helper functions for formatted logging
logi <- function(...) cat(sprintf(...), sep = "")
logl <- function(...) cat(sprintf(...), "\n", sep = "")

logl("Starting Seurat object creation (Multi-sample version)...")

# Cross-platform memory info
if (.Platform$OS.type == "windows") {
  logl("Memory limit (Windows): %.0f MB", utils::memory.limit())
} else {
  logl("Memory limit: not available on this OS (Linux/macOS)")
}

# ===== Step 1: Load Metadata =====
logl("[1/5] Loading metadata from: %s", metadata_path)
if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path)
}

metadata <- fread(metadata_path)
logl(
  "  Metadata dimensions: %d cells x %d columns",
  nrow(metadata),
  ncol(metadata)
)
logl("  Columns: %s", paste(names(metadata), collapse = ", "))

# Check required columns
if (!barcode_column %in% names(metadata)) {
  stop("Barcode column '", barcode_column, "' not found in metadata!")
}
if (!sample_column %in% names(metadata)) {
  stop("Sample column '", sample_column, "' not found in metadata!")
}
if (!celltype_column %in% names(metadata)) {
  stop("Cell type column '", celltype_column, "' not found in metadata!")
}

metadata_df <- as.data.frame(metadata)
rm(metadata)
gc()

# Check for NA in critical columns
if (anyNA(metadata_df[[barcode_column]])) {
  stop("Metadata barcode column contains NA values.")
}
if (anyNA(metadata_df[[sample_column]])) {
  stop("Metadata sample column contains NA values.")
}

metadata_df[[barcode_column]] <- as.character(metadata_df[[barcode_column]])
metadata_df[[sample_column]] <- as.character(metadata_df[[sample_column]])

# Count duplicates in original barcodes
dup_bc <- duplicated(metadata_df[[barcode_column]])
n_samples <- length(unique(metadata_df[[sample_column]]))

if (any(dup_bc)) {
  ndup <- sum(dup_bc)
  logl(
    "  [INFO] Found %d duplicated barcodes across %d samples",
    ndup,
    n_samples
  )
  logl("  This is expected for multi-sample data")
  logl("  Will add sample prefixes to ensure uniqueness")
}

# Create unique cell IDs: Sample_Barcode
metadata_df$cell_id <- paste(
  metadata_df[[sample_column]],
  metadata_df[[barcode_column]],
  sep = "_"
)

# Check uniqueness
if (anyDuplicated(metadata_df$cell_id)) {
  stop(
    "Even after adding sample prefix, cell_id still has duplicates!",
    "\nThis suggests duplicate entries in metadata for same sample+barcode combination."
  )
}

rownames(metadata_df) <- metadata_df$cell_id

# Convert to factor for memory efficiency
metadata_df[[sample_column]] <- as.factor(metadata_df[[sample_column]])
metadata_df[[celltype_column]] <- as.factor(metadata_df[[celltype_column]])

logl("  Unique samples: %d", nlevels(metadata_df[[sample_column]]))
logl("  Unique cell types: %d", nlevels(metadata_df[[celltype_column]]))
logl("  Created unique cell IDs: %d", nrow(metadata_df))

# ===== Step 2: Load Count Matrix =====
logl("[2/5] Loading count matrix from: %s", count_matrix_path)
if (!file.exists(count_matrix_path)) {
  stop("Count matrix file not found: ", count_matrix_path)
}

logl("  [INFO] Reading CSV and converting to sparse matrix...")

count_data <- fread(count_matrix_path, check.names = FALSE)
logl(
  "  CSV loaded: %d genes x %d cells",
  nrow(count_data),
  ncol(count_data) - 1
)

# Extract gene names
gene_col <- names(count_data)[1]
gene_names <- count_data[[gene_col]]
count_data[[gene_col]] <- NULL

# Get original cell barcodes
cell_barcodes <- names(count_data)

# Convert to numeric matrix
logl("  Converting to numeric matrix...")
count_matrix <- as.matrix(count_data)
storage.mode(count_matrix) <- "numeric"
rm(count_data)
gc()

rownames(count_matrix) <- make.unique(as.character(gene_names))
colnames(count_matrix) <- cell_barcodes

logl(
  "  Dense matrix size: %.2f GB",
  as.numeric(object.size(count_matrix)) / 1e9
)

# ===== Step 3: Add Sample Prefixes to Count Matrix =====
logl("[3/5] Adding sample prefixes to count matrix barcodes...")

if (assume_same_order) {
  # Assume count matrix columns are in same order as metadata rows
  if (ncol(count_matrix) != nrow(metadata_df)) {
    stop(
      "Count matrix columns (%d) != metadata rows (%d)",
      "\nCannot assume same order. Set assume_same_order=FALSE and provide barcode matching logic."
    )
  }

  logl("  [INFO] Assuming count matrix and metadata are in the same order")
  logl("  Verifying first 10 barcodes match...")

  # Get first 10 for verification
  n_check <- min(10, ncol(count_matrix))
  count_bc <- colnames(count_matrix)[1:n_check]
  meta_bc <- metadata_df[[barcode_column]][1:n_check]

  if (!all(count_bc == meta_bc)) {
    logl("  [WARNING] First 10 barcodes don't match exactly:")
    logl("    Count: %s", paste(count_bc, collapse = ", "))
    logl("    Metadata: %s", paste(meta_bc, collapse = ", "))
    stop(
      "Barcodes not in same order. Set assume_same_order=FALSE or fix input data."
    )
  }

  logl("  First 10 barcodes match - proceeding with prefix addition")

  # Add prefixes using metadata sample column (in same order)
  new_colnames <- paste(
    metadata_df[[sample_column]],
    colnames(count_matrix),
    sep = "_"
  )
  colnames(count_matrix) <- new_colnames
} else {
  stop(
    "Non-matching order not yet implemented in this version.",
    "\nPlease ensure count matrix and metadata are in the same order,",
    "\nor use Version 1 if count matrix already has unique barcodes."
  )
}

# Convert to sparse
logl("  Converting to sparse matrix...")
counts <- as(count_matrix, "dgCMatrix")
rm(count_matrix)
gc()

logl("  Sparse matrix size: %.2f GB", as.numeric(object.size(counts)) / 1e9)

# Calculate sparsity
nnz <- Matrix::nnzero(counts)
total <- prod(dim(counts))
sparsity <- 100 * (1 - nnz / total)
logl("  Non-zero elements: %d (%.2f%% sparse)", nnz, sparsity)

# ===== Step 4: Align Cells =====
logl("[4/5] Aligning cells with metadata...")

metadata_cells <- rownames(metadata_df) # Sample_Barcode
matrix_cells <- colnames(counts) # Sample_Barcode

logl("  Metadata cells: %d", length(metadata_cells))
logl("  Matrix cells: %d", length(matrix_cells))

# Should be perfect match if same order
common_cells <- intersect(matrix_cells, metadata_cells)
logl("  Common cells: %d", length(common_cells))

if (length(common_cells) != length(matrix_cells)) {
  logl("  [WARNING] Not all matrix cells found in metadata")
}

# Subset (should be no-op if perfect match)
metadata_df <- metadata_df[common_cells, , drop = FALSE]
counts <- counts[, common_cells, drop = FALSE]

logl("  Final dimensions: %d genes x %d cells", nrow(counts), ncol(counts))

# ===== Step 5: Create Seurat Object =====
logl("[5/5] Creating Seurat object...")

seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata_df,
  project = "FullAtlas",
  min.cells = 0,
  min.features = 0
)

logl("  Seurat object created successfully")
logl("  Object size: %.2f GB", as.numeric(object.size(seurat_obj)) / 1e9)
logl("  Number of cells: %d", ncol(seurat_obj))
logl("  Number of features: %d", nrow(seurat_obj))
logl("  Metadata columns: %d", ncol(seurat_obj@meta.data))

logl("\nSample distribution:")
print(table(seurat_obj@meta.data[[sample_column]]))

logl("\nCell type distribution (top 10):")
celltype_counts <- sort(
  table(seurat_obj@meta.data[[celltype_column]]),
  decreasing = TRUE
)
print(head(celltype_counts, 10))

# ===== Save =====
logl("\nSaving Seurat object to: %s", output_rds_path)

output_dir <- dirname(output_rds_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  logl("  Created output directory: %s", output_dir)
}

saveRDS(seurat_obj, file = output_rds_path, compress = "gzip")
logl("  Saved successfully!")
logl("  File size: %.2f GB", as.numeric(file.size(output_rds_path)) / 1e9)

# ===== Quality Check Summary =====
counts_mat <- seurat_obj[["RNA"]]@counts
logl("\n===== Quality Check Summary =====")
logl("Total cells: %d", ncol(seurat_obj))
logl("Total genes: %d", nrow(seurat_obj))
logl("Total UMI counts: %.2e", Matrix::sum(counts_mat))
logl("Mean UMI per cell: %.0f", mean(Matrix::colSums(counts_mat)))
logl("Mean genes per cell: %.0f", mean(Matrix::colSums(counts_mat > 0)))

logl("\n✅ Pipeline completed successfully!")
logl("Output file: %s", output_rds_path)
