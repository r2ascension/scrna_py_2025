# ==============================================================================
# Clean Up Old Layers/Assays
# ==============================================================================

# Remove specific layer from an assay
remove_layer <- function(seurat_obj, layer_name, assay_name = NULL) {
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(seurat_obj)
  }

  cat(sprintf(
    "\nRemoving layer '%s' from assay '%s'...\n",
    layer_name,
    assay_name
  ))

  # Check if layer exists
  current_layers <- Layers(seurat_obj, assay = assay_name)

  if (!layer_name %in% current_layers) {
    cat(sprintf("⚠️  Layer '%s' not found\n", layer_name))
    cat(sprintf(
      "Available layers: %s\n",
      paste(current_layers, collapse = ", ")
    ))
    return(seurat_obj)
  }

  # Remove layer by setting to NULL
  seurat_obj[[assay_name]][[layer_name]] <- NULL

  # Verify
  new_layers <- Layers(seurat_obj, assay = assay_name)

  if (layer_name %in% new_layers) {
    cat("❌ Failed to remove layer\n")
  } else {
    cat(sprintf("✓ Layer '%s' removed\n", layer_name))
    cat(sprintf("Remaining layers: %s\n", paste(new_layers, collapse = ", ")))
  }

  return(seurat_obj)
}


# Remove entire assay
remove_assay <- function(seurat_obj, assay_name) {
  cat(sprintf("\nRemoving assay '%s'...\n", assay_name))

  # Check if exists
  if (!assay_name %in% names(seurat_obj@assays)) {
    cat(sprintf("⚠️  Assay '%s' not found\n", assay_name))
    cat(sprintf(
      "Available assays: %s\n",
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
    return(seurat_obj)
  }

  # Don't remove if it's the only assay
  if (length(seurat_obj@assays) == 1) {
    stop("Cannot remove the only assay in object")
  }

  # Don't remove if it's the default
  if (assay_name == DefaultAssay(seurat_obj)) {
    cat("⚠️  This is the default assay\n")
    cat("Switching default to first available assay...\n")
    other_assays <- setdiff(names(seurat_obj@assays), assay_name)
    DefaultAssay(seurat_obj) <- other_assays[1]
    cat(sprintf("New default: %s\n", DefaultAssay(seurat_obj)))
  }

  # Remove
  seurat_obj[[assay_name]] <- NULL

  # Verify
  if (assay_name %in% names(seurat_obj@assays)) {
    cat("❌ Failed to remove assay\n")
  } else {
    cat(sprintf("✓ Assay '%s' removed\n", assay_name))
    cat(sprintf(
      "Remaining assays: %s\n",
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
  }

  return(seurat_obj)
}


# Clean all layers except specified ones
keep_only_layers <- function(
  seurat_obj,
  keep_layers = "counts",
  assay_name = NULL
) {
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(seurat_obj)
  }

  cat(sprintf("\nCleaning assay '%s'...\n", assay_name))
  cat(sprintf("Keeping only: %s\n", paste(keep_layers, collapse = ", ")))

  current_layers <- Layers(seurat_obj, assay = assay_name)
  layers_to_remove <- setdiff(current_layers, keep_layers)

  if (length(layers_to_remove) == 0) {
    cat("✓ Already clean, no layers to remove\n")
    return(seurat_obj)
  }

  cat(sprintf("Removing: %s\n", paste(layers_to_remove, collapse = ", ")))

  for (layer in layers_to_remove) {
    seurat_obj[[assay_name]][[layer]] <- NULL
  }

  # Verify
  final_layers <- Layers(seurat_obj, assay = assay_name)
  cat(sprintf("✓ Final layers: %s\n", paste(final_layers, collapse = ", ")))

  return(seurat_obj)
}


# ==============================================================================
# All-in-one: Convert and clean up in one command
# ==============================================================================

convert_and_cleanup <- function(
  seurat_obj,
  source_layer = "counts",
  new_assay_name = "RNA.symbol",
  handle_duplicates = "unique",
  remove_old_assay = TRUE,
  set_as_default = TRUE
) {
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install: BiocManager::install('org.Hs.eg.db')")
  }
  library(AnnotationDbi)

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Convert to Symbols and Clean Up\n")
  cat(rep("=", 70), "\n", sep = "")

  old_assay_name <- DefaultAssay(seurat_obj)

  # Step 1: Convert
  cat("\n[1/4] Converting to symbols...\n")

  mat <- LayerData(seurat_obj, layer = source_layer, assay = old_assay_name)
  cat(sprintf("  Source: %d genes × %d cells\n", nrow(mat), ncol(mat)))

  ensg_clean <- gsub("\\..*$", "", rownames(mat))

  symbol_map <- mapIds(
    org.Hs.eg.db,
    keys = ensg_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  new_symbols <- ifelse(
    is.na(symbol_map),
    rownames(mat),
    as.character(symbol_map)
  )

  n_conv <- sum(!is.na(symbol_map))
  cat(sprintf(
    "  Converted: %d / %d (%.1f%%)\n",
    n_conv,
    length(ensg_clean),
    n_conv / length(ensg_clean) * 100
  ))

  # Handle duplicates
  if (sum(duplicated(new_symbols)) > 0) {
    cat(sprintf("  Handling duplicates: %s\n", handle_duplicates))

    if (handle_duplicates == "sum") {
      unique_symbols <- unique(new_symbols)
      new_mat <- do.call(
        rbind,
        lapply(unique_symbols, function(sym) {
          idx <- which(new_symbols == sym)
          if (length(idx) == 1) {
            mat[idx, , drop = FALSE]
          } else {
            Matrix::colSums(mat[idx, , drop = FALSE])
          }
        })
      )
      rownames(new_mat) <- unique_symbols
      mat <- new_mat
    } else if (handle_duplicates == "unique") {
      rownames(mat) <- make.unique(new_symbols, sep = "_")
    } else {
      keep <- !duplicated(new_symbols)
      mat <- mat[keep, ]
    }
  } else {
    rownames(mat) <- new_symbols
  }

  # Step 2: Create new assay
  cat("\n[2/4] Creating new assay...\n")
  new_assay <- CreateAssayObject(counts = mat)
  seurat_obj[[new_assay_name]] <- new_assay
  cat(sprintf(
    "  ✓ Assay '%s' created (%d genes)\n",
    new_assay_name,
    nrow(seurat_obj[[new_assay_name]])
  ))

  # Step 3: Set as default
  if (set_as_default) {
    cat("\n[3/4] Setting as default assay...\n")
    DefaultAssay(seurat_obj) <- new_assay_name
    cat(sprintf("  ✓ Default assay: %s\n", DefaultAssay(seurat_obj)))
  } else {
    cat("\n[3/4] Keeping original default assay\n")
  }

  # Step 4: Remove old assay
  if (remove_old_assay) {
    cat("\n[4/4] Removing old assay...\n")
    cat(sprintf("  Removing: %s\n", old_assay_name))
    seurat_obj[[old_assay_name]] <- NULL
    cat("  ✓ Old assay removed\n")
  } else {
    cat("\n[4/4] Keeping old assay\n")
  }

  # Final status
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Final Status\n")
  cat(rep("=", 70), "\n", sep = "")

  cat(sprintf(
    "Object: %d genes × %d cells\n",
    nrow(seurat_obj),
    ncol(seurat_obj)
  ))
  cat(sprintf(
    "Available assays: %s\n",
    paste(names(seurat_obj@assays), collapse = ", ")
  ))
  cat(sprintf("Default assay: %s\n", DefaultAssay(seurat_obj)))
  cat(sprintf(
    "Gene format: %s\n",
    paste(head(rownames(seurat_obj), 5), collapse = ", ")
  ))

  cat(rep("=", 70), "\n\n", sep = "")

  return(seurat_obj)
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# # Example 1: Remove specific layer
# seurat_obj1 <- remove_layer(seurat_obj1, "counts.symbol")
#
# # Example 2: Remove old RNA assay (after creating RNA.symbol)
# seurat_obj1 <- remove_assay(seurat_obj1, "RNA")
#
# # Example 3: Keep only counts layer in RNA assay
# seurat_obj1 <- keep_only_layers(seurat_obj1, keep_layers = "counts", assay_name = "RNA")
#
# # Example 4: ONE COMMAND - Convert and remove old data
# seurat_obj1 <- convert_and_cleanup(
#   seurat_obj1,
#   source_layer = "counts",
#   new_assay_name = "RNA.symbol",
#   handle_duplicates = "sum",
#   remove_old_assay = TRUE,   # Remove old RNA assay
#   set_as_default = TRUE       # Set RNA.symbol as default
# )
#
# # Example 5: Convert but keep both assays
# seurat_obj1 <- convert_and_cleanup(
#   seurat_obj1,
#   remove_old_assay = FALSE    # Keep both RNA and RNA.symbol
# )

# sc_ref <- convert_and_cleanup(sc_ref,handle_duplicates = "unique", ,remove_old_assay = TRUE)
