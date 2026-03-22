# ==============================================================================
# Seurat v5 Utility: Layer/Assay Cleanup + ENSG->SYMBOL Conversion (Safe + Fast)
# Author: r2end (reviewed & hardened)
# Notes:
# - ENSG version stripping is applied ONLY to ENSG IDs (prevents AC245033.1 truncation)
# - Duplicate SYMBOL aggregation uses sparse matrix group-sum (fast, memory-safe)
# - Layer removal prefers LayerData<- (Seurat v5), with robust fallbacks
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required. Install.packages('Matrix')")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required. Install.packages('Seurat')")
  }
})

# ------------------------------------------------------------------------------
# Internal helper: remove a layer robustly (Seurat v5 first, then fallbacks)
# ------------------------------------------------------------------------------
.safe_remove_layer <- function(
  seurat_obj,
  assay_name,
  layer_name,
  verbose = TRUE
) {
  ok <- FALSE

  # v5: LayerData setter
  ok <- tryCatch(
    {
      LayerData(seurat_obj, assay = assay_name, layer = layer_name) <- NULL
      TRUE
    },
    error = function(e) FALSE
  )

  # fallback: Assay5 @layers list
  if (!ok) {
    ok <- tryCatch(
      {
        a <- seurat_obj[[assay_name]]
        if (!is.null(a@layers) && layer_name %in% names(a@layers)) {
          a@layers[[layer_name]] <- NULL
          seurat_obj[[assay_name]] <- a
          TRUE
        } else {
          FALSE
        }
      },
      error = function(e) FALSE
    )
  }

  # fallback: legacy slots (counts/data/scale.data) if user passes those names
  if (!ok && layer_name %in% c("counts", "data", "scale.data")) {
    ok <- tryCatch(
      {
        a <- seurat_obj[[assay_name]]
        slot(a, layer_name) <- new("dgCMatrix")
        seurat_obj[[assay_name]] <- a
        TRUE
      },
      error = function(e) FALSE
    )
  }

  if (verbose && !ok) {
    cat(sprintf(
      "⚠️  Could not remove layer '%s' from assay '%s' (unsupported structure)\n",
      layer_name,
      assay_name
    ))
  }

  return(seurat_obj)
}

# ------------------------------------------------------------------------------
# Remove specific layer from an assay
# ------------------------------------------------------------------------------
remove_layer <- function(seurat_obj, layer_name, assay_name = NULL) {
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(seurat_obj)
  }

  cat(sprintf(
    "\nRemoving layer '%s' from assay '%s'...\n",
    layer_name,
    assay_name
  ))

  current_layers <- tryCatch(
    Layers(seurat_obj, assay = assay_name),
    error = function(e) character(0)
  )

  if (length(current_layers) == 0) {
    cat("⚠️  No layers detected (or incompatible Seurat version).\n")
    return(seurat_obj)
  }

  if (!layer_name %in% current_layers) {
    cat(sprintf("⚠️  Layer '%s' not found\n", layer_name))
    cat(sprintf(
      "Available layers: %s\n",
      paste(current_layers, collapse = ", ")
    ))
    return(seurat_obj)
  }

  seurat_obj <- .safe_remove_layer(
    seurat_obj,
    assay_name,
    layer_name,
    verbose = FALSE
  )

  new_layers <- Layers(seurat_obj, assay = assay_name)
  if (layer_name %in% new_layers) {
    cat("❌ Failed to remove layer\n")
  } else {
    cat(sprintf("✓ Layer '%s' removed\n", layer_name))
    cat(sprintf("Remaining layers: %s\n", paste(new_layers, collapse = ", ")))
  }

  return(seurat_obj)
}

# ------------------------------------------------------------------------------
# Remove entire assay
# ------------------------------------------------------------------------------
remove_assay <- function(seurat_obj, assay_name) {
  cat(sprintf("\nRemoving assay '%s'...\n", assay_name))

  if (!assay_name %in% names(seurat_obj@assays)) {
    cat(sprintf("⚠️  Assay '%s' not found\n", assay_name))
    cat(sprintf(
      "Available assays: %s\n",
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
    return(seurat_obj)
  }

  if (length(seurat_obj@assays) == 1) {
    stop("Cannot remove the only assay in object")
  }

  if (assay_name == DefaultAssay(seurat_obj)) {
    cat("⚠️  This is the default assay\n")
    other_assays <- setdiff(names(seurat_obj@assays), assay_name)
    DefaultAssay(seurat_obj) <- other_assays[1]
    cat(sprintf("New default: %s\n", DefaultAssay(seurat_obj)))
  }

  seurat_obj[[assay_name]] <- NULL

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

# ------------------------------------------------------------------------------
# Clean all layers except specified ones
# ------------------------------------------------------------------------------
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

  current_layers <- tryCatch(
    Layers(seurat_obj, assay = assay_name),
    error = function(e) character(0)
  )

  if (length(current_layers) == 0) {
    cat("⚠️  No layers detected (or incompatible Seurat version).\n")
    return(seurat_obj)
  }

  layers_to_remove <- setdiff(current_layers, keep_layers)
  if (length(layers_to_remove) == 0) {
    cat("✓ Already clean, no layers to remove\n")
    return(seurat_obj)
  }

  cat(sprintf("Removing: %s\n", paste(layers_to_remove, collapse = ", ")))
  for (layer in layers_to_remove) {
    seurat_obj <- .safe_remove_layer(
      seurat_obj,
      assay_name,
      layer,
      verbose = FALSE
    )
  }

  final_layers <- Layers(seurat_obj, assay = assay_name)
  cat(sprintf("✓ Final layers: %s\n", paste(final_layers, collapse = ", ")))

  return(seurat_obj)
}

# ------------------------------------------------------------------------------
# Internal helper: fast sparse group-sum by factor labels
# ------------------------------------------------------------------------------
.sparse_group_sum_rows <- function(mat, group_labels) {
  # mat: genes x cells (dgCMatrix recommended)
  # group_labels: length nrow(mat) character
  f <- factor(group_labels, levels = unique(group_labels))
  G <- Matrix::sparseMatrix(
    i = seq_along(f),
    j = as.integer(f),
    x = 1,
    dims = c(length(f), nlevels(f))
  )
  out <- Matrix::t(G) %*% mat
  rownames(out) <- levels(f)
  colnames(out) <- colnames(mat)
  out
}

# ------------------------------------------------------------------------------
# All-in-one: Convert ENSG->SYMBOL and clean up
# ------------------------------------------------------------------------------
convert_and_cleanup <- function(
  seurat_obj,
  source_layer = "counts",
  new_assay_name = "RNA.symbol",
  handle_duplicates = c("sum", "unique", "first"),
  remove_old_assay = TRUE,
  set_as_default = TRUE
) {
  handle_duplicates <- match.arg(handle_duplicates)

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install: BiocManager::install('org.Hs.eg.db')")
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Install: BiocManager::install('AnnotationDbi')")
  }

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Convert to Symbols and Clean Up\n")
  cat(rep("=", 70), "\n", sep = "")

  old_assay_name <- DefaultAssay(seurat_obj)

  # Step 1: Extract matrix from layer
  cat("\n[1/4] Converting to symbols...\n")
  mat <- tryCatch(
    {
      LayerData(seurat_obj, layer = source_layer, assay = old_assay_name)
    },
    error = function(e) {
      stop(sprintf(
        "Failed to read LayerData(assay='%s', layer='%s'): %s",
        old_assay_name,
        source_layer,
        e$message
      ))
    }
  )

  cat(sprintf("  Source: %d genes × %d cells\n", nrow(mat), ncol(mat)))

  # Enforce sparse for memory safety
  if (!inherits(mat, "dgCMatrix")) {
    mat <- Matrix::as(mat, "dgCMatrix")
  }

  rn <- rownames(mat)

  # Only strip version for ENSG IDs; do NOT touch symbols like AC245033.1
  is_ensg <- grepl("^ENSG\\d+(\\.\\d+)?$", rn)
  ensg_base <- rn
  ensg_base[is_ensg] <- sub("\\..*$", "", rn[is_ensg])

  # Map only unique ENSG bases to reduce time
  ensg_unique <- unique(ensg_base[is_ensg])

  sym_unique <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensg_unique,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  # Construct new names with robust fallback
  new_symbols <- rn
  if (length(ensg_unique) > 0) {
    sym_mapped <- unname(sym_unique[ensg_base[is_ensg]])
    new_symbols[is_ensg] <- ifelse(
      !is.na(sym_mapped),
      sym_mapped,
      ensg_base[is_ensg]
    )
  }

  n_conv <- sum(is_ensg) - sum(is.na(sym_unique[ensg_unique]))
  cat(sprintf("  ENSG rows: %d\n", sum(is_ensg)))
  cat(sprintf(
    "  Converted (unique ENSG): %d / %d (%.1f%%)\n",
    sum(!is.na(sym_unique)),
    length(sym_unique),
    ifelse(
      length(sym_unique) == 0,
      0,
      100 * sum(!is.na(sym_unique)) / length(sym_unique)
    )
  ))

  # Handle duplicates
  dup_n <- sum(duplicated(new_symbols))
  if (dup_n > 0) {
    cat(sprintf("  Duplicated gene names after mapping: %d\n", dup_n))
    cat(sprintf("  Handling duplicates: %s\n", handle_duplicates))

    if (handle_duplicates == "sum") {
      mat <- .sparse_group_sum_rows(mat, new_symbols)
    } else if (handle_duplicates == "unique") {
      rownames(mat) <- make.unique(new_symbols, sep = "_")
    } else {
      # "first"
      keep <- !duplicated(new_symbols)
      mat <- mat[keep, , drop = FALSE]
      rownames(mat) <- new_symbols[keep]
    }
  } else {
    rownames(mat) <- new_symbols
  }

  # Step 2: Create new assay
  cat("\n[2/4] Creating new assay...\n")
  new_assay <- Seurat::CreateAssayObject(counts = mat)
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
    seurat_obj <- remove_assay(seurat_obj, old_assay_name)
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
    "Gene head: %s\n",
    paste(head(rownames(seurat_obj), 5), collapse = ", ")
  ))
  cat(rep("=", 70), "\n\n", sep = "")

  return(seurat_obj)
}

# ==============================================================================
# Example usage:
# sc_ref <- convert_and_cleanup(
#   sc_ref,
#   source_layer = "counts",
#   new_assay_name = "RNA",
#   handle_duplicates = "sum",
#   remove_old_assay = TRUE,
#   set_as_default = TRUE
# )
#
# seurat_obj <- keep_only_layers(seurat_obj, keep_layers = c("counts"), assay_name = "RNA.symbol")
# ==============================================================================

# sc_ref <- convert_and_cleanup(sc_ref,handle_duplicates = "unique", ,remove_old_assay = TRUE)

fix_counts_assay_name <- function(
  seurat_obj,
  from_assay = "counts",
  to_assay = "RNA",
  set_default = TRUE,
  overwrite = FALSE
) {
  # ---- checks ----
  if (!from_assay %in% names(seurat_obj@assays)) {
    stop(sprintf(
      "Assay '%s' not found. Available assays: %s",
      from_assay,
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
  }

  if (to_assay %in% names(seurat_obj@assays)) {
    if (!overwrite) {
      # auto-backup existing to_assay
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup <- paste0(to_assay, ".bak_", ts)
      message(sprintf(
        "Target assay '%s' already exists -> backing up to '%s'",
        to_assay,
        backup
      ))
      seurat_obj[[backup]] <- seurat_obj[[to_assay]]
      seurat_obj[[to_assay]] <- NULL
    } else {
      message(sprintf("Overwriting existing assay '%s'", to_assay))
      seurat_obj[[to_assay]] <- NULL
    }
  }

  # ---- move/rename assay ----
  seurat_obj[[to_assay]] <- seurat_obj[[from_assay]]

  if (set_default) {
    DefaultAssay(seurat_obj) <- to_assay
  }

  # ---- update reductions assay.used (important) ----
  if (length(seurat_obj@reductions) > 0) {
    for (rd in names(seurat_obj@reductions)) {
      au <- tryCatch(
        seurat_obj@reductions[[rd]]@assay.used,
        error = function(e) NA
      )
      if (!is.na(au) && identical(au, from_assay)) {
        seurat_obj@reductions[[rd]]@assay.used <- to_assay
      }
    }
  }

  # ---- remove old assay name ----
  seurat_obj[[from_assay]] <- NULL

  message(sprintf(
    "✓ Renamed assay '%s' -> '%s'; DefaultAssay = '%s'",
    from_assay,
    to_assay,
    DefaultAssay(seurat_obj)
  ))
  message(sprintf(
    "Assays now: %s",
    paste(names(seurat_obj@assays), collapse = ", ")
  ))
  message(sprintf(
    "Layers in '%s': %s",
    to_assay,
    paste(Layers(seurat_obj, assay = to_assay), collapse = ", ")
  ))

  return(seurat_obj)
}

# ---- one-click run ----
# sc_ref <- fix_counts_assay_name(sc_ref, from_assay = "counts", to_assay = "RNA", set_default = TRUE)
# 或者你更想显式一点：
# sc_ref <- fix_counts_assay_name(sc_ref, from_assay = "counts", to_assay = "RNA.symbol", set_default = TRUE)
