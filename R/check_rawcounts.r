library(Seurat)
# Load required libraries
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
setwd("/home/h2048/data/R/1124")
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
# source("path/to/SCNT.R")

# seurat_obj <- readRDS("/home/h2048/data/R/1124/merged_object_sc.rds")
counts <- GetAssayData(seurat_obj, slot = "counts")[2000:3000, 2000:3000]
m <- max(counts)
cat(sprintf("Max: %.2f -> ", m))
if (m > 100) {
  cat("✓ RAW COUNTS\n")
} else if (m < 20) {
  cat("✗ LOG-TRANSFORMED\n")
} else {
  cat("⚠️  UNCLEAR\n")
}

#' GetH5ad - Improved Version for Seurat v5
#'
#' Convert a Seurat object into a .h5ad file with robust raw counts extraction.
#' Supports both single-cell (sc) and spatial transcriptomics (st) modes.
#' Compatible with Seurat v5 Assay5 objects.
#'
#' @param seurat_obj A Seurat object to be converted
#' @param output_path Path where the H5AD file will be saved
#' @param mode Export mode: "sc" (single-cell) or "st" (spatial transcriptomics)
#' @param assay Name of the assay to use. Default is "RNA"
#' @param coord_cols Coordinate column names. Default c("x", "y")
#' @param cell_col Column name containing cell identifiers. Default "cell"
#' @param scale Scale for tissue coordinates ("hires" or "lowres"). Default "hires"
#' @param external_coords Optional external coordinates dataframe
#' @param validate_counts Logical. If TRUE, validates that counts are raw. Default TRUE
#' @param debug Logical. If TRUE, prints debugging information. Default FALSE
#'
#' @return No return value. The H5AD file will be saved at the specified location.
#' @export

GetH5ad <- function(
  seurat_obj,
  output_path,
  mode = c("sc", "st"),
  assay = "RNA",
  coord_cols = c("x", "y"),
  cell_col = "cell",
  scale = "hires",
  external_coords = NULL,
  validate_counts = TRUE,
  debug = FALSE
) {
  mode <- match.arg(mode)

  # ============================================================================
  # Check Python Environment
  # ============================================================================
  if (!reticulate::py_available(initialize = TRUE)) {
    stop("Python environment not detected. Please configure Python first.")
  }

  tryCatch(
    {
      anndata <- reticulate::import("anndata", delay_load = FALSE)
      np <- reticulate::import("numpy", delay_load = FALSE)
      pd <- reticulate::import("pandas", delay_load = FALSE)
    },
    error = function(e) {
      stop(
        "Required Python packages not found. Install with:\n",
        "  reticulate::py_install(c('anndata', 'numpy', 'pandas'))\n",
        "Error: ",
        e$message
      )
    }
  )

  # ============================================================================
  # Extract Raw Counts Matrix (Seurat v5 Compatible)
  # ============================================================================
  if (debug) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Step 1: Extracting raw counts matrix\n")
    cat(rep("=", 70), "\n\n", sep = "")
    cat("Assay type:", class(seurat_obj[[assay]])[1], "\n")
  }

  # Try multiple methods to extract counts (Seurat v5 compatible)
  counts_matrix <- tryCatch(
    {
      # Method 1: LayerData (Seurat v5 preferred)
      if (debug) {
        cat("  Trying LayerData(layer='counts')...\n")
      }
      Seurat::LayerData(seurat_obj, assay = assay, layer = "counts")
    },
    error = function(e1) {
      tryCatch(
        {
          # Method 2: GetAssayData with layer parameter
          if (debug) {
            cat("  Trying GetAssayData(layer='counts')...\n")
          }
          Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
        },
        error = function(e2) {
          tryCatch(
            {
              # Method 3: GetAssayData with slot parameter (Seurat v3/v4)
              if (debug) {
                cat("  Trying GetAssayData(slot='counts')...\n")
              }
              Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
            },
            error = function(e3) {
              # Method 4: Direct access
              if (debug) {
                cat("  Trying direct access to @counts...\n")
              }
              seurat_obj[[assay]]@counts
            }
          )
        }
      )
    }
  )

  if (is.null(counts_matrix) || length(counts_matrix) == 0) {
    stop("Failed to extract counts matrix from Seurat object!")
  }

  # ============================================================================
  # Validate Raw Counts (MUST PASS BEFORE PROCEEDING)
  # ============================================================================
  if (validate_counts) {
    if (debug) {
      cat("\n", rep("=", 70), "\n", sep = "")
      cat("CRITICAL: Raw Counts Validation\n")
      cat(rep("=", 70), "\n\n", sep = "")
    }

    # Sample a subset for validation (up to 3000 cells and genes for robust check)
    sample_size_cells <- min(3000, ncol(counts_matrix))
    sample_size_genes <- min(3000, nrow(counts_matrix))

    if (
      inherits(counts_matrix, "dgCMatrix") || inherits(counts_matrix, "Matrix")
    ) {
      sample_counts <- as.matrix(counts_matrix[
        1:sample_size_genes,
        1:sample_size_cells
      ])
    } else {
      sample_counts <- counts_matrix[1:sample_size_genes, 1:sample_size_cells]
    }

    max_val <- max(sample_counts)
    mean_val <- mean(sample_counts)

    # Check if values look like raw counts (should be integers or close to integers)
    is_integer_like <- all(sample_counts == round(sample_counts))

    if (debug) {
      cat(sprintf(
        "Sample statistics (first %d x %d):\n",
        sample_size_genes,
        sample_size_cells
      ))
      cat(sprintf("  Max value: %.2f\n", max_val))
      cat(sprintf("  Mean value: %.4f\n", mean_val))
      cat(sprintf("  Integer-like: %s\n", is_integer_like))
    }

    # STRICT VALIDATION - Stop if data looks log-transformed
    if (max_val < 20 && mean_val < 2 && !is_integer_like) {
      stop(
        "\n",
        "❌ VALIDATION FAILED: Data appears to be log-transformed!\n",
        "   Max value: ",
        round(max_val, 2),
        "\n",
        "   Mean value: ",
        round(mean_val, 4),
        "\n",
        "   Integer-like: ",
        is_integer_like,
        "\n\n",
        "H5AD export requires RAW COUNTS for downstream analysis.\n",
        "Please verify you are extracting the 'counts' layer.\n\n",
        "To bypass validation (not recommended), use: validate_counts=FALSE"
      )
    }

    # Warning if max value is suspiciously low
    if (is_integer_like && max_val < 10) {
      warning(
        "\n",
        "⚠️  WARNING: Max count value is very low (",
        round(max_val, 2),
        ")\n",
        "   This might indicate:\n",
        "   1. Very shallow sequencing depth\n",
        "   2. Wrong data layer extracted\n",
        "   3. Filtered/subset data\n",
        "Proceeding, but please verify this is expected...\n"
      )
    }

    if (is_integer_like && max_val >= 10) {
      if (debug) {
        cat("\n✅ VALIDATION PASSED: Raw counts confirmed\n")
        cat("   Proceeding with export...\n")
      }
    }
  }

  # ============================================================================
  # Prepare Matrix for Export
  # ============================================================================
  if (debug) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Step 2: Preparing matrix for export\n")
    cat(rep("=", 70), "\n\n", sep = "")
  }

  # Ensure matrix is in proper format
  if (
    !inherits(counts_matrix, "matrix") &&
      !inherits(counts_matrix, "dgCMatrix")
  ) {
    counts_matrix <- as.matrix(counts_matrix)
  }

  original_cell_ids <- colnames(counts_matrix)
  gene_names <- rownames(counts_matrix)

  if (debug) {
    cat(sprintf(
      "Original dimensions: %d genes × %d cells\n",
      nrow(counts_matrix),
      ncol(counts_matrix)
    ))
    cat(
      "First 5 cell IDs:",
      paste(head(original_cell_ids, 5), collapse = ", "),
      "\n"
    )
    cat("First 5 genes:", paste(head(gene_names, 5), collapse = ", "), "\n")
  }

  # Transpose: Seurat (genes × cells) → AnnData (cells × genes)
  counts_matrix <- Matrix::t(counts_matrix)
  counts_matrix <- as(counts_matrix, "dgCMatrix")

  if (debug) {
    cat(sprintf(
      "\nTransposed dimensions: %d cells × %d genes\n",
      nrow(counts_matrix),
      ncol(counts_matrix)
    ))
  }

  # ============================================================================
  # Prepare and Clean Metadata
  # ============================================================================
  if (debug) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Step 3: Preparing metadata\n")
    cat(rep("=", 70), "\n\n", sep = "")
  }

  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)

  if (debug) {
    cat(sprintf("Original metadata: %d columns\n", ncol(meta_data)))
  }

  # Clean metadata: convert all columns to appropriate types for H5AD
  meta_data_clean <- as.data.frame(
    lapply(meta_data, function(col) {
      # Convert factors to character
      if (is.factor(col)) {
        col <- as.character(col)
      }

      # Handle different data types
      if (is.numeric(col)) {
        # Keep numeric as is, but replace NaN with NA
        col[is.nan(col)] <- NA
        return(col)
      } else if (is.logical(col)) {
        # Convert logical to character to avoid issues
        return(as.character(col))
      } else if (is.character(col)) {
        # Replace NA with "NA" string for character columns
        col[is.na(col)] <- "NA"
        return(col)
      } else {
        # For any other type, convert to character
        col <- as.character(col)
        col[is.na(col)] <- "NA"
        return(col)
      }
    }),
    stringsAsFactors = FALSE
  )

  # Restore rownames
  rownames(meta_data_clean) <- rownames(meta_data)

  if (debug) {
    cat(sprintf("Cleaned metadata: %d columns\n", ncol(meta_data_clean)))
    cat("\nColumn types after cleaning:\n")
    type_summary <- table(sapply(meta_data_clean, class))
    for (type in names(type_summary)) {
      cat(sprintf("  %s: %d columns\n", type, type_summary[type]))
    }
  }

  # ============================================================================
  # Create AnnData Object
  # ============================================================================
  if (debug) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Step 4: Creating AnnData object\n")
    cat(rep("=", 70), "\n\n", sep = "")
  }

  # Create AnnData object with cleaned metadata
  adata <- anndata$AnnData(
    X = counts_matrix,
    obs = pd$DataFrame(meta_data_clean)
  )
  adata$var_names <- np$array(gene_names)

  if (debug) {
    cat("AnnData object created successfully\n")
    cat(sprintf(
      "  Shape: %d cells × %d genes\n",
      nrow(counts_matrix),
      ncol(counts_matrix)
    ))
  }

  # ============================================================================
  # Add Dimensionality Reductions
  # ============================================================================
  if (length(seurat_obj@reductions) > 0) {
    if (debug) {
      cat("\nAdding dimensionality reductions:\n")
    }

    for (reduction_name in names(seurat_obj@reductions)) {
      embeddings <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)
      obsm_key <- paste0("X_", reduction_name)
      adata$obsm[obsm_key] <- np$array(embeddings)

      if (debug) {
        cat(sprintf(
          "  ✓ Added %s (%d dimensions)\n",
          obsm_key,
          ncol(embeddings)
        ))
      }
    }
  }

  # ============================================================================
  # Add Spatial Information (if mode = "st")
  # ============================================================================
  if (mode == "st") {
    if (debug) {
      cat("\n", rep("=", 70), "\n", sep = "")
      cat("Step 5: Adding spatial information\n")
      cat(rep("=", 70), "\n\n", sep = "")
    }

    # Get coordinates
    if (is.null(external_coords)) {
      tryCatch(
        {
          coords <- Seurat::GetTissueCoordinates(seurat_obj, scale = scale)
          if (debug) {
            cat("Retrieved coordinates using GetTissueCoordinates()\n")
            cat(sprintf("  Dimensions: %d × %d\n", nrow(coords), ncol(coords)))
          }
        },
        error = function(e) {
          stop(
            "Failed to get coordinates: ",
            e$message,
            "\nConsider providing external_coords parameter."
          )
        }
      )
    } else {
      coords <- external_coords
      if (debug) {
        cat("Using provided external coordinates\n")
        cat(sprintf("  Dimensions: %d × %d\n", nrow(coords), ncol(coords)))
      }
    }

    # Prepare spatial coordinates
    if (!is.null(cell_col) && cell_col %in% colnames(coords)) {
      spatial_coords <- coords[, coord_cols, drop = FALSE]
      rownames(spatial_coords) <- coords[[cell_col]]
    } else {
      spatial_coords <- coords[, coord_cols, drop = FALSE]
    }

    # Add to AnnData
    adata$obsm$`__setitem__`("spatial", np$array(as.matrix(spatial_coords)))

    # Add spatial metadata
    if (length(seurat_obj@images) > 0) {
      image_name <- names(seurat_obj@images)[1]
      image_obj <- seurat_obj@images[[image_name]]
      scale_factors <- image_obj@scale.factors

      adata$uns[["spatial"]] <- reticulate::r_to_py(list(
        sample = list(
          images = list(),
          scalefactors = scale_factors,
          metadata = list()
        )
      ))

      if (debug) {
        cat("Added spatial metadata\n")
      }
    }
  }

  # ============================================================================
  # Write H5AD File
  # ============================================================================
  if (debug) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Step 6: Writing H5AD file\n")
    cat(rep("=", 70), "\n\n", sep = "")
  }

  adata$write_h5ad(output_path)

  cat("\n✓ H5AD file successfully created!\n")
  cat(sprintf("  Output: %s\n", output_path))
  cat(sprintf("  Cells: %d\n", nrow(counts_matrix)))
  cat(sprintf("  Genes: %d\n", ncol(counts_matrix)))

  if (mode == "st") {
    cat("\nNote: Spatial image files are not exported.\n")
    cat("Save tissue_lowres_image.png manually to a 'spatial/' folder.\n")
  }

  invisible(NULL)
}

# ==============================================================================
# Usage Examples
# ==============================================================================

# Example 1: Single-cell mode (basic)
# GetH5ad(seurat_obj, "output_sc.h5ad", mode = "sc")

# Example 2: Single-cell mode with debugging
GetH5ad(seurat_obj, "output_sc_1204.h5ad", mode = "sc", debug = TRUE)

# Example 3: Spatial transcriptomics mode
# GetH5ad(seurat_obj, "output_st.h5ad", mode = "st", scale = "hires")

# Example 4: With custom assay and validation disabled
# GetH5ad(seurat_obj, "output.h5ad", mode = "sc",
#         assay = "SCT", validate_counts = FALSE)

# Example 5: Spatial with external coordinates
# GetH5ad(seurat_obj, "output_st.h5ad", mode = "st",
#         external_coords = my_coords_df, debug = TRUE)
