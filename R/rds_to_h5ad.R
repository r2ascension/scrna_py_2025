# ===== Batch Convert Seurat RDS to H5AD using SCNT =====
# Purpose: Convert cell-type-specific Harmony-integrated Seurat objects to h5ad format
# Author: r2end
# Date: 2024-12-23

# ----- 1. Setup Environment -----
library(Seurat)
library(SCNT) # Contains GetH5ad() function
library(reticulate)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()

# Input/Output paths
input_dir <- "/home/h2048/data/R/1221/per_celltype_harmony_rogue/seurat_objects"
output_dir <- "/home/h2048/data/R/1223/per_celltype_harmony_rogue/h5ad_objects"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# ----- 2. Get All RDS Files -----
rds_files <- list.files(
  input_dir,
  pattern = "_harmony\\.rds$",
  full.names = TRUE
)
cat(sprintf(
  "\n========== Found %d RDS files to convert ==========\n",
  length(rds_files)
))
print(basename(rds_files))

# ----- 3. Batch Conversion Function -----
convert_rds_to_h5ad <- function(rds_path, output_dir) {
  # Extract cell type name
  file_name <- basename(rds_path)
  cell_type <- gsub("_harmony\\.rds$", "", file_name)

  cat(sprintf("\n[%s] Processing: %s\n", Sys.time(), cell_type))

  # Define output path
  h5ad_path <- file.path(output_dir, paste0(cell_type, "_harmony.h5ad"))

  # Skip if h5ad already exists
  if (file.exists(h5ad_path)) {
    cat(sprintf("  ✓ H5AD already exists, skipping: %s\n", basename(h5ad_path)))
    return(invisible(NULL))
  }

  tryCatch(
    {
      # Load Seurat object
      cat(sprintf("  → Loading RDS...\n"))
      seurat_obj <- readRDS(rds_path)

      # Print basic info
      n_cells <- ncol(seurat_obj)
      n_features <- nrow(seurat_obj)
      assays <- names(seurat_obj@assays)
      reductions <- names(seurat_obj@reductions)

      cat(sprintf("  → Cells: %d, Features: %d\n", n_cells, n_features))
      cat(sprintf("  → Assays: %s\n", paste(assays, collapse = ", ")))
      cat(sprintf("  → Reductions: %s\n", paste(reductions, collapse = ", ")))

      # Convert to h5ad using SCNT::GetH5ad()
      cat(sprintf("  → Converting to h5ad...\n"))
      GetH5ad(
        seurat_obj = seurat_obj,
        output_path = h5ad_path,
        mode = "sc", # Single-cell mode
        assay = "RNA", # Default RNA assay
        debug = FALSE # Set to TRUE if need debugging
      )

      # Verify file creation
      if (file.exists(h5ad_path)) {
        file_size <- file.info(h5ad_path)$size / 1024^2 # MB
        cat(sprintf(
          "  ✓ Successfully created: %s (%.2f MB)\n",
          basename(h5ad_path),
          file_size
        ))
      } else {
        cat(sprintf("  ✗ Failed to create h5ad file\n"))
      }

      # Clean up memory
      rm(seurat_obj)
      gc(verbose = FALSE)
    },
    error = function(e) {
      cat(sprintf("  ✗ Error converting %s: %s\n", cell_type, e$message))
    }
  )
}

# ----- 4. Execute Batch Conversion -----
cat("\n========== Starting Batch Conversion ==========\n")

start_time <- Sys.time()

for (i in seq_along(rds_files)) {
  cat(sprintf("\n--- [%d/%d] ---\n", i, length(rds_files)))
  convert_rds_to_h5ad(rds_files[i], output_dir)
}

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\n========== Conversion Complete ==========\n"))
cat(sprintf("Total time elapsed: %.2f minutes\n", elapsed_time))

# ----- 5. Summary Report -----
h5ad_files <- list.files(output_dir, pattern = "\\.h5ad$", full.names = TRUE)
cat(sprintf("\nSuccessfully created %d h5ad files:\n", length(h5ad_files)))

# Create summary table
summary_df <- data.frame(
  CellType = gsub("_harmony\\.h5ad$", "", basename(h5ad_files)),
  FileName = basename(h5ad_files),
  SizeMB = sapply(h5ad_files, function(x) file.info(x)$size / 1024^2),
  row.names = NULL
)

print(summary_df, row.names = FALSE)

cat(sprintf("\nOutput directory: %s\n", output_dir))
