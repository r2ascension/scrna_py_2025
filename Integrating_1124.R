################################################################################
# Complete Integrated Pipeline: Gene Standardization + QC + Smart Merge
# 完整整合流程：基因标准化 + 质控 + 智能合并
#
# Workflow:
# 1. Read all RDS files from folder
# 2. Gene name standardization (mapping + synonym removal)
# 3. 7-step QC pipeline for each sample
# 4. Smart merge all cleaned samples
# 5. Build batch-level gene availability matrix
#
# Author: r2end
# Date: 2025-11-24
################################################################################

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

#===============================================================================
# Part 1: Gene Standardization Functions (from original workflow)
#===============================================================================

#' Get counts matrix (compatible with all Seurat versions)
get_counts_matrix <- function(seurat_obj, assay = "RNA") {
  counts <- tryCatch(
    {
      LayerData(seurat_obj, assay = assay, layer = "counts")
    },
    error = function(e1) {
      tryCatch(
        {
          GetAssayData(seurat_obj, assay = assay, layer = "counts")
        },
        error = function(e2) {
          tryCatch(
            {
              GetAssayData(seurat_obj, assay = assay, slot = "counts")
            },
            error = function(e3) {
              seurat_obj[[assay]]@counts
            }
          )
        }
      )
    }
  )
  return(counts)
}

#' Create Assay object (compatible with all Seurat versions)
create_assay_object <- function(counts_matrix) {
  assay <- tryCatch(
    {
      CreateAssay5Object(counts = counts_matrix)
    },
    error = function(e) {
      CreateAssayObject(counts = counts_matrix)
    }
  )
  return(assay)
}

#' Map gene names to official symbols
map_gene_names <- function(genes, gene_db = org.Hs.eg.db) {
  gene_mapping <- data.frame(
    original_name = genes,
    official_symbol = NA_character_,
    mapping_source = NA_character_,
    stringsAsFactors = FALSE
  )

  # Layer 1: SYMBOL
  symbol_match <- tryCatch(
    {
      mapIds(
        gene_db,
        keys = genes,
        column = "SYMBOL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    },
    error = function(e) rep(NA, length(genes))
  )

  matched <- !is.na(symbol_match)
  gene_mapping$official_symbol[matched] <- symbol_match[matched]
  gene_mapping$mapping_source[matched] <- "SYMBOL"

  # Layer 2: ALIAS
  unmatched <- is.na(gene_mapping$official_symbol)
  if (sum(unmatched) > 0) {
    alias_match <- tryCatch(
      {
        mapIds(
          gene_db,
          keys = genes[unmatched],
          column = "SYMBOL",
          keytype = "ALIAS",
          multiVals = "first"
        )
      },
      error = function(e) rep(NA, sum(unmatched))
    )

    alias_matched <- !is.na(alias_match)
    gene_mapping$official_symbol[unmatched][alias_matched] <- alias_match[
      alias_matched
    ]
    gene_mapping$mapping_source[unmatched][alias_matched] <- "ALIAS"
  }

  # Layer 3: ENSEMBL
  unmatched <- is.na(gene_mapping$official_symbol)
  if (sum(unmatched) > 0) {
    ensembl_match <- tryCatch(
      {
        mapIds(
          gene_db,
          keys = genes[unmatched],
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      },
      error = function(e) rep(NA, sum(unmatched))
    )

    ensembl_matched <- !is.na(ensembl_match)
    gene_mapping$official_symbol[unmatched][ensembl_matched] <- ensembl_match[
      ensembl_matched
    ]
    gene_mapping$mapping_source[unmatched][ensembl_matched] <- "ENSEMBL"
  }

  # Layer 4: Keep original
  unmatched <- is.na(gene_mapping$official_symbol)
  gene_mapping$official_symbol[unmatched] <- genes[unmatched]
  gene_mapping$mapping_source[unmatched] <- "UNMAPPED"

  return(gene_mapping)
}

#' Remove synonym duplicates
remove_synonym_duplicates <- function(seurat_obj, gene_mapping, assay = "RNA") {
  symbol_counts <- table(gene_mapping$official_symbol)
  duplicate_symbols <- names(symbol_counts)[symbol_counts > 1]

  if (length(duplicate_symbols) == 0) {
    counts_matrix <- get_counts_matrix(seurat_obj, assay)
    rownames(counts_matrix) <- gene_mapping$official_symbol
    seurat_obj[[assay]] <- create_assay_object(counts_matrix)

    return(list(
      cleaned_obj = seurat_obj,
      synonym_report = NULL
    ))
  }

  cat("  Found", length(duplicate_symbols), "synonym groups\n")

  synonym_report <- list()
  genes_to_keep <- logical(nrow(seurat_obj))
  names(genes_to_keep) <- rownames(seurat_obj)
  genes_to_keep[] <- TRUE

  for (symbol in duplicate_symbols) {
    original_names <- gene_mapping$original_name[
      gene_mapping$official_symbol == symbol
    ]
    genes_to_keep[original_names[-1]] <- FALSE

    synonym_report[[symbol]] <- data.frame(
      official_symbol = symbol,
      n_synonyms = length(original_names),
      kept = original_names[1],
      removed = paste(original_names[-1], collapse = " | "),
      stringsAsFactors = FALSE
    )
  }

  seurat_obj <- seurat_obj[genes_to_keep, ]

  gene_mapping_kept <- gene_mapping[
    gene_mapping$original_name %in% rownames(seurat_obj),
  ]
  gene_mapping_kept <- gene_mapping_kept[
    match(rownames(seurat_obj), gene_mapping_kept$original_name),
  ]

  counts_matrix <- get_counts_matrix(seurat_obj, assay)
  rownames(counts_matrix) <- gene_mapping_kept$official_symbol
  seurat_obj[[assay]] <- create_assay_object(counts_matrix)

  synonym_report_df <- if (length(synonym_report) > 0) {
    do.call(rbind, synonym_report)
  } else {
    NULL
  }

  cat("  ✓ Removed", sum(!genes_to_keep), "duplicate genes\n")

  return(list(
    cleaned_obj = seurat_obj,
    synonym_report = synonym_report_df
  ))
}

#===============================================================================
# Part 2: QC Functions (streamlined, minimal saving)
#===============================================================================

#' Run EmptyDrops filtering (by sample)
run_emptydrops_qc <- function(
  seurat_obj,
  sample_column = "orig.ident",
  lower = 100,
  fdr_threshold = 0.01,
  niters = 10000,
  verbose = TRUE
) {
  if (verbose) {
    cat("  [EmptyDrops] Running by sample...\n")
  }

  cells_before <- ncol(seurat_obj)

  # Check if sample column exists
  if (!sample_column %in% colnames(seurat_obj@meta.data)) {
    if (verbose) {
      cat(sprintf(
        "    Warning: column '%s' not found, treating as single sample\n",
        sample_column
      ))
    }
    samples <- rep("all", ncol(seurat_obj))
  } else {
    samples <- seurat_obj@meta.data[[sample_column]]
  }

  unique_samples <- unique(samples)
  all_cell_indices <- c()

  for (sample_id in unique_samples) {
    sample_cells <- which(samples == sample_id)
    sample_obj <- seurat_obj[, sample_cells]

    if (verbose) {
      cat(sprintf(
        "    Sample: %s (%d cells)\n",
        sample_id,
        length(sample_cells)
      ))
    }

    mat <- get_counts_matrix(sample_obj)

    set.seed(42)
    e.out <- DropletUtils::emptyDrops(
      m = mat,
      lower = lower,
      niters = niters,
      test.ambient = TRUE
    )

    is.cell <- e.out$FDR <= fdr_threshold
    is.cell[is.na(is.cell)] <- FALSE

    # Get original indices
    kept_indices <- sample_cells[is.cell]
    all_cell_indices <- c(all_cell_indices, kept_indices)

    if (verbose) {
      cat(sprintf(
        "      Kept: %d/%d (%.1f%%)\n",
        sum(is.cell),
        length(is.cell),
        sum(is.cell) / length(is.cell) * 100
      ))
    }
  }

  cells_after <- length(all_cell_indices)

  if (verbose) {
    cat(sprintf(
      "    Total - Before: %d, After: %d, Removed: %d (%.1f%%)\n",
      cells_before,
      cells_after,
      cells_before - cells_after,
      (cells_before - cells_after) / cells_before * 100
    ))
  }

  return(list(
    cell_indices = all_cell_indices,
    stats = c(before = cells_before, after = cells_after)
  ))
}

#' Run regression-based outlier detection (by sample)
run_outlier_qc <- function(
  seurat_obj,
  sample_column = "orig.ident",
  outlier_threshold = 0.999,
  verbose = TRUE
) {
  if (verbose) {
    cat("  [Outlier Detection] Running by sample...\n")
  }

  cells_before <- ncol(seurat_obj)

  # Check if sample column exists
  if (!sample_column %in% colnames(seurat_obj@meta.data)) {
    if (verbose) {
      cat(sprintf(
        "    Warning: column '%s' not found, treating as single sample\n",
        sample_column
      ))
    }
    samples <- rep("all", ncol(seurat_obj))
  } else {
    samples <- seurat_obj@meta.data[[sample_column]]
  }

  unique_samples <- unique(samples)
  all_cell_names <- c()

  for (sample_id in unique_samples) {
    sample_cells <- which(samples == sample_id)
    sample_obj <- seurat_obj[, sample_cells]

    if (verbose) {
      cat(sprintf(
        "    Sample: %s (%d cells)\n",
        sample_id,
        length(sample_cells)
      ))
    }

    nCount <- sample_obj$nCount_RNA
    nFeature <- sample_obj$nFeature_RNA

    log_counts <- log10(nCount + 1)
    log_features <- log10(nFeature + 1)

    fit <- lm(log_features ~ log_counts)
    pred <- predict(
      fit,
      data.frame(log_counts = log_counts),
      interval = "prediction",
      level = outlier_threshold
    )

    outliers <- log_features < pred[, "lwr"] | log_features > pred[, "upr"]

    kept_cells <- colnames(sample_obj)[!outliers]
    all_cell_names <- c(all_cell_names, kept_cells)

    if (verbose) {
      cat(sprintf(
        "      Kept: %d/%d (%.1f%%)\n",
        sum(!outliers),
        length(outliers),
        sum(!outliers) / length(outliers) * 100
      ))
    }
  }

  cells_after <- length(all_cell_names)

  if (verbose) {
    cat(sprintf(
      "    Total - Before: %d, After: %d, Removed: %d (%.1f%%)\n",
      cells_before,
      cells_after,
      cells_before - cells_after,
      (cells_before - cells_after) / cells_before * 100
    ))
  }

  return(list(
    cell_names = all_cell_names,
    stats = c(before = cells_before, after = cells_after)
  ))
}

#' Run DoubletFinder (by sample, Seurat v5 / Assay5 compatible)
run_doubletfinder_qc <- function(
  seurat_obj,
  sample_column = "sample",
  doublet_rate = 0.06,
  n_pcs = 30,
  verbose = TRUE
) {
  if (verbose) {
    cat("  [DoubletFinder] Running by sample...\n")
  }

  cells_before <- ncol(seurat_obj)

  # 统一从 meta.data 里拿 sample 定义
  if (!sample_column %in% colnames(seurat_obj@meta.data)) {
    if (verbose) {
      cat(sprintf(
        "    Warning: column '%s' not found in meta.data, treating as single sample\n",
        sample_column
      ))
    }
    samples <- rep("all", ncol(seurat_obj))
  } else {
    samples <- seurat_obj@meta.data[[sample_column]]
  }

  unique_samples <- unique(samples)
  all_cell_names <- character(0)

  for (sample_id in unique_samples) {
    sample_cells <- which(samples == sample_id)
    sample_obj <- seurat_obj[, sample_cells, drop = FALSE]
    n_cells <- ncol(sample_obj)

    if (verbose) {
      cat(sprintf("    Sample: %s (%d cells)\n", sample_id, n_cells))
    }

    if (n_cells < 50) {
      if (verbose) {
        cat("      Skipped (too few cells)\n")
      }
      all_cell_names <- c(all_cell_names, colnames(sample_obj))
      next
    }

    nExp <- round(doublet_rate * n_cells)

    result <- tryCatch(
      {
        DefaultAssay(sample_obj) <- DefaultAssay(seurat_obj)

        # ⚠️ 关键修改：不要再用 @data，统一用 GetAssayData()
        data_mat <- tryCatch(
          GetAssayData(sample_obj, slot = "data"),
          error = function(e) NULL
        )

        if (is.null(data_mat) || nrow(data_mat) == 0) {
          sample_obj <- NormalizeData(sample_obj, verbose = FALSE)
        }

        sample_obj <- FindVariableFeatures(
          sample_obj,
          selection.method = "vst",
          nfeatures = 2000,
          verbose = FALSE
        )
        sample_obj <- ScaleData(
          sample_obj,
          features = VariableFeatures(sample_obj),
          verbose = FALSE
        )

        n_pcs_use <- min(n_pcs, max(10, ncol(sample_obj) - 1))
        sample_obj <- RunPCA(sample_obj, npcs = n_pcs_use, verbose = FALSE)

        sweep.res <- paramSweep(sample_obj, PCs = 1:n_pcs_use, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)

        pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        if (is.na(pK) || pK < 0.01 || pK > 0.3) {
          pK <- 0.09
        }

        sample_obj <- doubletFinder(
          sample_obj,
          PCs = 1:n_pcs_use,
          pN = 0.25,
          pK = pK,
          nExp = nExp,
          sct = FALSE
        )

        df_class_col <- grep(
          "DF.classifications",
          colnames(sample_obj@meta.data),
          value = TRUE
        )

        if (length(df_class_col) > 0) {
          df_class_col <- df_class_col[length(df_class_col)]
          status_values <- sample_obj@meta.data[[df_class_col]]

          singlet_cells <- colnames(sample_obj)[status_values == "Singlet"]

          if (verbose) {
            cat(sprintf(
              "      Kept: %d/%d singlets (%.1f%%)\n",
              length(singlet_cells),
              n_cells,
              length(singlet_cells) / n_cells * 100
            ))
          }

          list(status = "Success", cells = singlet_cells)
        } else {
          list(status = "Failed", cells = colnames(sample_obj))
        }
      },
      error = function(e) {
        if (verbose) {
          cat(sprintf("      ✗ ERROR: %s\n", e$message))
          if (!is.null(e$call)) {
            cat(sprintf(
              "      Call: %s\n",
              paste(deparse(e$call), collapse = " ")
            ))
          }
        }
        list(status = "Error", cells = colnames(sample_obj))
      }
    )

    all_cell_names <- c(all_cell_names, result$cells)
  }

  cells_after <- length(all_cell_names)

  if (verbose) {
    removed <- cells_before - cells_after
    rate <- if (cells_before > 0) removed / cells_before * 100 else 0
    cat(sprintf(
      "    Total - Before: %d, After: %d, Removed: %d (%.1f%%)\n",
      cells_before,
      cells_after,
      removed,
      rate
    ))
  }

  return(list(
    status = "Complete",
    cell_names = all_cell_names,
    stats = c(before = cells_before, after = cells_after)
  ))
}

#' Run DecontX (batch-aware)
run_decontx_qc <- function(
  seurat_obj,
  batch_column = NULL,
  max_iter = 500,
  verbose = TRUE
) {
  if (verbose) {
    if (
      !is.null(batch_column) && batch_column %in% colnames(seurat_obj@meta.data)
    ) {
      cat(sprintf(
        "  [DecontX] Running with batch correction (column: %s)...\n",
        batch_column
      ))
    } else {
      cat("  [DecontX] Running without batch correction...\n")
    }
  }

  counts_matrix <- get_counts_matrix(seurat_obj)

  if (!inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix <- as(counts_matrix, "dgCMatrix")
  }

  decontX_params <- list(x = counts_matrix, maxIter = max_iter)

  if (
    !is.null(batch_column) && batch_column %in% colnames(seurat_obj@meta.data)
  ) {
    batch_vector <- seurat_obj@meta.data[[batch_column]]
    decontX_params$batch <- batch_vector

    if (verbose) {
      n_batches <- length(unique(batch_vector))
      cat(sprintf("    Processing %d samples/batches\n", n_batches))
    }
  }

  set.seed(12345)
  decontX_results <- do.call(decontX, decontX_params)

  contamination_vector <- decontX_results$contamination
  seurat_obj$decontX_contamination <- contamination_vector

  if (!is.null(decontX_results$z)) {
    seurat_obj$decontX_clusters <- as.factor(decontX_results$z)
  }

  decontX_counts <- round(decontX_results$decontXcounts)
  seurat_obj[["decontXcounts"]] <- create_assay_object(decontX_counts)

  if (verbose) {
    cat(sprintf(
      "    Mean contamination: %.2f%%\n",
      mean(contamination_vector) * 100
    ))

    # Show per-batch statistics if batch column is provided
    if (
      !is.null(batch_column) && batch_column %in% colnames(seurat_obj@meta.data)
    ) {
      batch_stats <- aggregate(
        contamination_vector,
        by = list(seurat_obj@meta.data[[batch_column]]),
        FUN = function(x) sprintf("%.2f%%", mean(x) * 100)
      )
      colnames(batch_stats) <- c("Sample", "Mean_Contamination")
      cat("    Per-sample contamination:\n")
      for (i in 1:nrow(batch_stats)) {
        cat(sprintf(
          "      %s: %s\n",
          batch_stats$Sample[i],
          batch_stats$Mean_Contamination[i]
        ))
      }
    }
  }

  return(seurat_obj)
}

#===============================================================================
# Part 3: Integrated Pipeline for Single Sample
#===============================================================================

#' Process single RDS file: standardization + QC
process_single_sample <- function(
  rds_path,
  gene_db = org.Hs.eg.db,
  # QC parameters
  sample_column = "orig.ident", # NEW: column name for sample identification
  nfeature_min = 200,
  nfeature_max = 6000,
  mt_pattern = "^MT-",
  mt_max = 20,
  rb_pattern = "^RP[SL]",
  rb_max = 40,
  emptydrops_lower = 100,
  emptydrops_fdr = 0.01,
  outlier_threshold = 0.999,
  doublet_rate = 0.06,
  doublet_pcs = 30,
  decontx_batch_column = NULL,
  decontx_max_iter = 500,
  verbose = TRUE
) {
  file_name <- basename(rds_path)
  file_base <- tools::file_path_sans_ext(file_name)

  if (verbose) {
    cat("\n")
    cat("════════════════════════════════════════════════════════════\n")
    cat(sprintf("Processing: %s\n", file_name))
    cat("════════════════════════════════════════════════════════════\n")
  }

  # Track statistics
  qc_stats <- data.frame(
    step = character(),
    cells_before = integer(),
    cells_after = integer(),
    stringsAsFactors = FALSE
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 0: Load and clean
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[0] Loading and cleaning...\n")
  }

  seurat_obj <- readRDS(rds_path)
  initial_cells <- ncol(seurat_obj)
  initial_genes <- nrow(seurat_obj)

  seurat_obj <- DietSeurat(
    seurat_obj,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL
  )

  if (verbose) {
    cat(sprintf("  Cells: %d, Genes: %d\n", initial_cells, initial_genes))
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Step 1: Gene standardization
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[1] Gene name standardization...\n")
  }

  original_genes <- rownames(seurat_obj)
  gene_mapping <- map_gene_names(original_genes, gene_db)

  n_mapped <- sum(gene_mapping$mapping_source != "UNMAPPED")
  if (verbose) {
    cat(sprintf(
      "  Mapped: %d/%d (%.1f%%)\n",
      n_mapped,
      length(original_genes),
      100 * n_mapped / length(original_genes)
    ))
  }

  result <- remove_synonym_duplicates(seurat_obj, gene_mapping)
  seurat_obj <- result$cleaned_obj

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "1_Gene_Standardization",
      cells_before = initial_cells,
      cells_after = ncol(seurat_obj),
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 2: nFeature filtering
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[2] nFeature filtering...\n")
  }

  cells_before <- ncol(seurat_obj)
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= nfeature_min &
      nFeature_RNA <= nfeature_max
  )
  cells_after <- ncol(seurat_obj)

  if (verbose) {
    cat(sprintf("  Range: %d-%d genes\n", nfeature_min, nfeature_max))
    cat(sprintf(
      "  Before: %d, After: %d, Removed: %d (%.1f%%)\n",
      cells_before,
      cells_after,
      cells_before - cells_after,
      (cells_before - cells_after) / cells_before * 100
    ))
  }

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "2_nFeature",
      cells_before = cells_before,
      cells_after = cells_after,
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 3: Mitochondrial filtering
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[3] Mitochondrial filtering...\n")
  }

  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
      seurat_obj,
      pattern = mt_pattern
    )
  }

  cells_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = percent.mt < mt_max)
  cells_after <- ncol(seurat_obj)

  if (verbose) {
    cat(sprintf("  Threshold: < %d%%\n", mt_max))
    cat(sprintf(
      "  Before: %d, After: %d, Removed: %d (%.1f%%)\n",
      cells_before,
      cells_after,
      cells_before - cells_after,
      (cells_before - cells_after) / cells_before * 100
    ))
  }

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "3_Mitochondrial",
      cells_before = cells_before,
      cells_after = cells_after,
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 4: Ribosomal filtering
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[4] Ribosomal filtering...\n")
  }

  if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
      seurat_obj,
      pattern = rb_pattern
    )
  }

  cells_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = percent.rb < rb_max)
  cells_after <- ncol(seurat_obj)

  if (verbose) {
    cat(sprintf("  Threshold: < %d%%\n", rb_max))
    cat(sprintf(
      "  Before: %d, After: %d, Removed: %d (%.1f%%)\n",
      cells_before,
      cells_after,
      cells_before - cells_after,
      (cells_before - cells_after) / cells_before * 100
    ))
  }

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "4_Ribosomal",
      cells_before = cells_before,
      cells_after = cells_after,
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 5: EmptyDrops (optional, may fail if already filtered)
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[5] EmptyDrops filtering...\n")
  }

  emptydrops_result <- tryCatch(
    {
      run_emptydrops_qc(
        seurat_obj,
        sample_column,
        emptydrops_lower,
        emptydrops_fdr,
        verbose = verbose
      )
    },
    error = function(e) {
      if (verbose) {
        cat("  Skipped (likely already filtered data)\n")
      }
      return(NULL)
    }
  )

  if (!is.null(emptydrops_result)) {
    seurat_obj <- seurat_obj[, emptydrops_result$cell_indices]
    qc_stats <- rbind(
      qc_stats,
      data.frame(
        step = "5_EmptyDrops",
        cells_before = emptydrops_result$stats["before"],
        cells_after = emptydrops_result$stats["after"],
        stringsAsFactors = FALSE
      )
    )
  } else {
    cells_current <- ncol(seurat_obj)
    qc_stats <- rbind(
      qc_stats,
      data.frame(
        step = "5_EmptyDrops",
        cells_before = cells_current,
        cells_after = cells_current,
        stringsAsFactors = FALSE
      )
    )
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Step 6: Outlier detection
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[6] Outlier detection...\n")
  }

  outlier_result <- run_outlier_qc(
    seurat_obj,
    sample_column,
    outlier_threshold,
    verbose
  )
  seurat_obj <- seurat_obj[, outlier_result$cell_names]

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "6_Outlier",
      cells_before = outlier_result$stats["before"],
      cells_after = outlier_result$stats["after"],
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 7: DoubletFinder
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[7] DoubletFinder...\n")
  }

  doublet_result <- run_doubletfinder_qc(
    seurat_obj,
    sample_column,
    doublet_rate,
    doublet_pcs,
    verbose
  )

  if (doublet_result$status == "Complete") {
    seurat_obj <- seurat_obj[, doublet_result$cell_names]
  }

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "7_DoubletFinder",
      cells_before = doublet_result$stats["before"],
      cells_after = doublet_result$stats["after"],
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Step 8: DecontX
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n[8] DecontX...\n")
  }

  cells_before <- ncol(seurat_obj)
  seurat_obj <- run_decontx_qc(
    seurat_obj,
    decontx_batch_column,
    decontx_max_iter,
    verbose
  )

  qc_stats <- rbind(
    qc_stats,
    data.frame(
      step = "8_DecontX",
      cells_before = cells_before,
      cells_after = ncol(seurat_obj),
      stringsAsFactors = FALSE
    )
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Summary
  # ═══════════════════════════════════════════════════════════════════════
  if (verbose) {
    cat("\n")
    cat("────────────────────────────────────────────────────────────\n")
    cat("Summary\n")
    cat("────────────────────────────────────────────────────────────\n")
    cat(sprintf("Initial: %d cells, %d genes\n", initial_cells, initial_genes))
    cat(sprintf(
      "Final: %d cells, %d genes\n",
      ncol(seurat_obj),
      nrow(seurat_obj)
    ))
    cat(sprintf(
      "Total removed: %d cells (%.1f%%)\n",
      initial_cells - ncol(seurat_obj),
      (initial_cells - ncol(seurat_obj)) / initial_cells * 100
    ))
    cat("────────────────────────────────────────────────────────────\n\n")
  }

  return(list(
    seurat_obj = seurat_obj,
    qc_stats = qc_stats,
    file_name = file_base,
    initial_cells = initial_cells,
    final_cells = ncol(seurat_obj),
    initial_genes = initial_genes,
    final_genes = nrow(seurat_obj),
    gene_mapping = gene_mapping,
    synonym_report = result$synonym_report
  ))
}

#===============================================================================
# Part 4: Main Pipeline - Process Folder + Smart Merge
#===============================================================================

#' Complete pipeline: standardization + QC + merge
run_complete_pipeline <- function(
  input_dir,
  output_dir = "pipeline_output",
  species = "human",
  # QC parameters
  sample_column = "orig.ident", # Column name for sample identification (for EmptyDrops, Outlier, DoubletFinder)
  nfeature_min = 200,
  nfeature_max = 6000,
  mt_pattern = "^MT-",
  mt_max = 20,
  rb_pattern = "^RP[SL]",
  rb_max = 40,
  emptydrops_lower = 100,
  emptydrops_fdr = 0.01,
  outlier_threshold = 0.999,
  doublet_rate = 0.06,
  doublet_pcs = 30,
  decontx_batch_column = NULL, # Batch column for DecontX (if NULL, uses sample_column)
  decontx_max_iter = 500,
  # Merge parameters
  add_cell_ids = TRUE,
  # Output parameters
  save_individual_samples = TRUE,
  build_gene_availability = TRUE,
  batch_column = "orig.ident",
  verbose = TRUE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Complete Integrated Pipeline                            ║\n")
    cat("║  Gene Standardization + QC + Smart Merge                 ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }

  start_time <- Sys.time()

  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cleaned_dir <- file.path(output_dir, "cleaned_samples")
  if (save_individual_samples) {
    dir.create(cleaned_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Select gene database
  gene_db <- if (species == "human") {
    org.Hs.eg.db
  } else if (species == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Please install: BiocManager::install('org.Mm.eg.db')")
    }
    library(org.Mm.eg.db)
    org.Mm.eg.db
  } else {
    stop("Unsupported species")
  }

  # Find all RDS files
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

  if (length(rds_files) == 0) {
    stop("No RDS files found in input directory")
  }

  if (verbose) {
    cat(sprintf("Found %d RDS files\n", length(rds_files)))
    cat(sprintf("Species: %s\n", species))
    cat(sprintf("Output directory: %s\n\n", output_dir))
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Process each sample
  # ═══════════════════════════════════════════════════════════════════════

  processed_samples <- list()
  all_qc_stats <- list()

  for (i in seq_along(rds_files)) {
    if (verbose) {
      cat(sprintf(
        "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
      ))
      cat(sprintf("Sample %d/%d\n", i, length(rds_files)))
      cat(sprintf(
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
      ))
    }

    result <- process_single_sample(
      rds_path = rds_files[i],
      gene_db = gene_db,
      sample_column = sample_column, # NEW: pass sample column
      nfeature_min = nfeature_min,
      nfeature_max = nfeature_max,
      mt_pattern = mt_pattern,
      mt_max = mt_max,
      rb_pattern = rb_pattern,
      rb_max = rb_max,
      emptydrops_lower = emptydrops_lower,
      emptydrops_fdr = emptydrops_fdr,
      outlier_threshold = outlier_threshold,
      doublet_rate = doublet_rate,
      doublet_pcs = doublet_pcs,
      decontx_batch_column = if (!is.null(decontx_batch_column)) {
        decontx_batch_column
      } else {
        sample_column
      },
      decontx_max_iter = decontx_max_iter,
      verbose = verbose
    )

    # Save individual sample if requested
    if (save_individual_samples) {
      output_file <- file.path(
        cleaned_dir,
        paste0(result$file_name, "_cleaned.rds")
      )
      saveRDS(result$seurat_obj, output_file)
      if (verbose) cat(sprintf("\nSaved: %s\n", basename(output_file)))
    }

    # Store results
    processed_samples[[result$file_name]] <- result$seurat_obj
    all_qc_stats[[result$file_name]] <- result$qc_stats

    # Save gene mapping
    mapping_dir <- file.path(output_dir, "gene_mappings")
    dir.create(mapping_dir, showWarnings = FALSE)
    fwrite(
      result$gene_mapping,
      file.path(mapping_dir, paste0(result$file_name, "_mapping.csv"))
    )

    # Save synonym report if exists
    if (!is.null(result$synonym_report)) {
      synonym_dir <- file.path(output_dir, "synonym_reports")
      dir.create(synonym_dir, showWarnings = FALSE)
      fwrite(
        result$synonym_report,
        file.path(synonym_dir, paste0(result$file_name, "_synonyms.csv"))
      )
    }

    # Clean up
    gc(verbose = FALSE)
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Smart Merge
  # ═══════════════════════════════════════════════════════════════════════

  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Smart Merge                                             ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }

  sample_names <- names(processed_samples)

  if (verbose) {
    cat(sprintf("Merging %d samples...\n", length(processed_samples)))
  }

  if (add_cell_ids) {
    merged_obj <- merge(
      x = processed_samples[[1]],
      y = processed_samples[-1],
      add.cell.ids = sample_names
    )
  } else {
    merged_obj <- merge(x = processed_samples[[1]], y = processed_samples[-1])
  }

  if (verbose) {
    cat("✓ Merge complete\n")
    cat(sprintf("  Total cells: %d\n", ncol(merged_obj)))
    cat(sprintf("  Total genes: %d\n\n", nrow(merged_obj)))
  }

  # Save merged object
  merged_file <- file.path(output_dir, "merged_seurat.rds")
  saveRDS(merged_obj, merged_file)
  if (verbose) {
    cat(sprintf("Saved: %s\n\n", basename(merged_file)))
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Build gene availability matrix (optional, by metadata batch_column)
  # ═══════════════════════════════════════════════════════════════════════

  if (build_gene_availability) {
    if (verbose) {
      cat("╔══════════════════════════════════════════════════════════╗\n")
      cat("║  Building Gene Availability Matrix                      ║\n")
      cat("╚══════════════════════════════════════════════════════════╝\n\n")
    }

    availability_dir <- file.path(output_dir, "gene_availability")
    dir.create(availability_dir, showWarnings = FALSE)

    # 如果 batch_column 没给，就用 sample_column
    batch_column_use <- if (!is.null(batch_column)) {
      batch_column
    } else {
      sample_column
    }

    if (!batch_column_use %in% colnames(merged_obj@meta.data)) {
      warning(sprintf(
        "batch_column '%s' not found in meta.data; using a single batch 'all'.",
        batch_column_use
      ))
      merged_obj$.__batch_for_availability__ <- "all"
      batch_column_use <- ".__batch_for_availability__"
    }

    counts_mat <- get_counts_matrix(merged_obj)
    all_genes <- rownames(counts_mat)
    batches <- unique(merged_obj@meta.data[[batch_column_use]])

    if (verbose) {
      cat(sprintf("Total unique genes: %d\n", length(all_genes)))
      cat(sprintf(
        "Number of batches (%s): %d\n\n",
        batch_column_use,
        length(batches)
      ))
    }

    batch_summary <- data.frame(
      batch = character(),
      n_genes = integer(),
      n_cells = integer(),
      stringsAsFactors = FALSE
    )

    for (b in batches) {
      cells_idx <- which(merged_obj@meta.data[[batch_column_use]] == b)
      if (length(cells_idx) == 0) {
        next
      }

      present <- Matrix::rowSums(counts_mat[, cells_idx, drop = FALSE] > 0) > 0L

      availability_df <- data.frame(
        gene = all_genes,
        available = as.integer(present),
        stringsAsFactors = FALSE
      )

      csv_file <- file.path(availability_dir, paste0("batch_", b, ".csv"))
      fwrite(availability_df, csv_file)

      batch_summary <- rbind(
        batch_summary,
        data.frame(
          batch = b,
          n_genes = sum(present),
          n_cells = length(cells_idx),
          stringsAsFactors = FALSE
        )
      )
    }

    fwrite(batch_summary, file.path(availability_dir, "batch_summary.csv"))
    writeLines(all_genes, file.path(availability_dir, "all_genes.txt"))

    if (verbose) {
      cat("✓ Gene availability matrices saved\n\n")
    }
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Generate final report
  # ═══════════════════════════════════════════════════════════════════════

  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")

  # Combine all QC stats
  combined_stats <- do.call(
    rbind,
    lapply(names(all_qc_stats), function(name) {
      stats <- all_qc_stats[[name]]
      stats$sample <- name
      stats
    })
  )

  fwrite(combined_stats, file.path(output_dir, "all_samples_qc_stats.csv"))

  # Generate text report
  generate_pipeline_report(
    output_dir = output_dir,
    n_samples = length(processed_samples),
    merged_obj = merged_obj,
    combined_stats = combined_stats,
    elapsed_time = elapsed_time
  )

  if (verbose) {
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Pipeline Complete!                                      ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("Total time: %.2f minutes\n", as.numeric(elapsed_time)))
    cat(sprintf("Output directory: %s\n", output_dir))
    cat(sprintf("Final merged object: %s\n", basename(merged_file)))
    cat("\n")
  }

  return(list(
    merged_obj = merged_obj,
    processed_samples = processed_samples,
    qc_stats = combined_stats,
    elapsed_time = elapsed_time,
    output_dir = output_dir
  ))
}

#===============================================================================
# Part 5: Report Generation
#===============================================================================

generate_pipeline_report <- function(
  output_dir,
  n_samples,
  merged_obj,
  combined_stats,
  elapsed_time
) {
  report_file <- file.path(output_dir, "PIPELINE_REPORT.txt")

  sink(report_file)

  cat("═══════════════════════════════════════════════════════════\n")
  cat("Integrated Pipeline Report\n")
  cat("Gene Standardization + QC + Smart Merge\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  cat(sprintf("Analysis date: %s\n", Sys.Date()))
  cat(sprintf("Processing time: %.2f minutes\n\n", as.numeric(elapsed_time)))

  cat("───────────────────────────────────────────────────────────\n")
  cat("Overview\n")
  cat("───────────────────────────────────────────────────────────\n")
  cat(sprintf("Number of samples processed: %d\n", n_samples))
  cat(sprintf("Final merged cells: %d\n", ncol(merged_obj)))
  cat(sprintf("Final merged genes: %d\n\n", nrow(merged_obj)))

  cat("───────────────────────────────────────────────────────────\n")
  cat("Per-sample Summary\n")
  cat("───────────────────────────────────────────────────────────\n")

  sample_summary <- combined_stats %>%
    group_by(sample) %>%
    summarise(
      initial_cells = first(cells_before),
      final_cells = last(cells_after),
      removal_rate = (first(cells_before) - last(cells_after)) /
        first(cells_before) *
        100,
      .groups = "drop"
    )

  print(sample_summary)
  cat("\n")

  cat("───────────────────────────────────────────────────────────\n")
  cat("QC Steps Across All Samples\n")
  cat("───────────────────────────────────────────────────────────\n")

  step_summary <- combined_stats %>%
    group_by(step) %>%
    summarise(
      total_before = sum(cells_before),
      total_after = sum(cells_after),
      total_removed = sum(cells_before - cells_after),
      removal_rate = (sum(cells_before) - sum(cells_after)) /
        sum(cells_before) *
        100,
      .groups = "drop"
    )

  print(step_summary)
  cat("\n")

  cat("═══════════════════════════════════════════════════════════\n")
  cat("Output Files\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("cleaned_samples/          # Individual cleaned samples\n")
  cat("gene_mappings/            # Gene name mapping tables\n")
  cat("synonym_reports/          # Synonym removal reports\n")
  cat("gene_availability/        # Gene availability matrices\n")
  cat("merged_seurat.rds         # Final merged object\n")
  cat("all_samples_qc_stats.csv  # Combined QC statistics\n")
  cat("PIPELINE_REPORT.txt       # This report\n")
  cat("═══════════════════════════════════════════════════════════\n")

  sink()

  cat(sprintf("Report saved: %s\n", basename(report_file)))
}

#===============================================================================
# Usage Example
#===============================================================================

# Run complete pipeline
results <- run_complete_pipeline(
  input_dir = "/home/h2048/data/source/final",
  output_dir = "/home/h2048/data/R/1124",
  species = "human",

  # QC parameters
  sample_column = "sample", # ✅ 统一使用 meta.data$sample（可按需要改）
  nfeature_min = 200,
  nfeature_max = 6000,
  mt_pattern = "^MT-",
  mt_max = 20,
  rb_pattern = "^RP[SL]",
  rb_max = 40,
  emptydrops_lower = 100,
  emptydrops_fdr = 0.01,
  outlier_threshold = 0.999,
  doublet_rate = 0.06,
  doublet_pcs = 30,
  decontx_batch_column = "sample", # 如果 NULL，就用 sample_column
  decontx_max_iter = 500,
  # Merge parameters
  add_cell_ids = TRUE,
  # Output parameters
  save_individual_samples = TRUE,
  build_gene_availability = TRUE,
  batch_column = "sample", # 如果 NULL，就用 sample_column

  # Merge parameters
  verbose = TRUE
)

# Access results
merged_seurat <- results$merged_obj

merged_seurat <- readRDS('/home/h2048/data/R/1124/merged_seurat.rds')

head(merged_seurat)
merged_seurat <- subset(merged_seurat, subset = decontX_contamination <= 0.25)

library(Seurat)
library(Matrix)

## ------------ 基本设置 ------------
library(Seurat)
library(Matrix)

assay_use <- DefaultAssay(merged_seurat) # 一般是 "RNA"

## 1. 先把多层 counts join 成一个
# 看看现在有哪些 layer（调试用）
print(Layers(merged_seurat[[assay_use]]))
# 例如会打印: "counts.1" "counts.2" "counts.3" "data.1" ...

# join 所有 layer
merged_seurat[[assay_use]] <- JoinLayers(merged_seurat[[assay_use]])

# 再确认一下，现在应该是标准的 counts / data / scale.data 结构
print(Layers(merged_seurat[[assay_use]]))
# 例如: "counts" "data" "scale.data"

## 2. 取 counts 矩阵
cnt <- LayerData(
  object = merged_seurat,
  assay = assay_use,
  layer = "counts"
)

## 3. 按 min.cells = 3 过滤基因
gene_ncells <- Matrix::rowSums(cnt > 0)
keep_genes <- names(gene_ncells[gene_ncells >= 3])

cat("保留基因数: ", length(keep_genes), "\n")

merged_seurat <- subset(merged_seurat, features = keep_genes)
saveRDS(merged_seurat, '/home/h2048/data/source/final/merged_seurat.rds')
