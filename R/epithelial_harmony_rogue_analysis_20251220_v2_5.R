#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

# =========================
# Configuration
# =========================
INPUT_RDS  <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/R/1218/harmony_subset_rogue_v2_5"

MANUAL_COL <- "Manual_Annotation"
BATCH_COL  <- "dataset"

TARGET_MANUAL <- NULL
EXCLUDE_CELLTYPES <- c("Unknown", "Doublet", "LowQuality", "Skip", "Undefined", "Mixed")

# P1-2 FIX: Celltype analyzability thresholds
MIN_CELLS_PER_CELLTYPE <- 200   # Skip if celltype has <200 cells total
MIN_CELLS_PER_BATCH    <- 20    # Skip if any batch has <20 cells

N_VAR_FEATURES <- 4000
N_PCS          <- 50
CLUSTER_RES    <- 2
CLUSTER_ALGO   <- 4

MARKER_MIN_PCT <- 0.25
MARKER_LOGFC   <- 0.25
TOP_N_HEATMAP  <- 10
MAX_CELLS_PER_IDENT <- 5000

RUN_ROGUE <- TRUE
MIN_CELLS_PER_GROUP   <- 30
ROGUE_MIN_CELLS_GENE  <- 10
ROGUE_MIN_GENES_CELL  <- 200
ROGUE_MAX_CELLS_GROUP <- 2000
ROGUE_MAX_GENES_GROUP <- 3000
ROGUE_MAX_COMBOS      <- 500

SAVE_PRE_HARMONY <- TRUE

RANDOM_SEED <- 42

# =========================
# IO
# =========================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_DIR   <- file.path(OUTPUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
TAB_DIR   <- file.path(OUTPUT_DIR, "tables");  dir.create(TAB_DIR, showWarnings = FALSE)
OBJ_DIR   <- file.path(OUTPUT_DIR, "objects"); dir.create(OBJ_DIR, showWarnings = FALSE)
LOG_FILE  <- file.path(OUTPUT_DIR, "analysis.log")

sink(LOG_FILE, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)

cat("========================================\n")
cat("Harmony + ROGUE Analysis (v2.5)\n")
cat("========================================\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("Input :", INPUT_RDS, "\n")
cat("Output:", OUTPUT_DIR, "\n\n")

cat("=== Analysis Design Notes ===\n")
cat("This pipeline analyzes each cell type INDEPENDENTLY:\n")
cat("- Each type has separate HVG, PCA, Harmony embedding\n")
cat("- Results NOT directly comparable across cell types\n")
cat("- Use for intra-celltype substructure discovery\n")
cat("- For cross-celltype comparisons, integrate all types together first\n\n")

# =========================
# Input Validation
# =========================
cat("=== Input Validation ===\n")
if (!file.exists(INPUT_RDS)) {
  stop(sprintf("❌ Input not found: %s", INPUT_RDS))
}

seurat_obj <- readRDS(INPUT_RDS)

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("❌ RNA assay not found")
}

DefaultAssay(seurat_obj) <- "RNA"

if (!MANUAL_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("❌ '%s' not in meta.data", MANUAL_COL))
}

if (!BATCH_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("❌ '%s' not in meta.data", BATCH_COL))
}

counts_available <- FALSE
tryCatch({
  test_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  if (!is.null(test_counts) && length(test_counts) > 0) {
    counts_available <- TRUE
    cat("✓ counts layer detected (v5)\n")
  }
}, error = function(e) {
  tryCatch({
    test_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    if (!is.null(test_counts) && length(test_counts) > 0) {
      counts_available <<- TRUE
      cat("✓ counts slot detected (v4)\n")
    }
  }, error = function(e2) {})
})

if (!counts_available) {
  stop("❌ Cannot access counts")
}

n_batches <- length(unique(seurat_obj@meta.data[[BATCH_COL]]))
if (n_batches < 2) {
  stop(sprintf("❌ Only %d batch. Need ≥2", n_batches))
}

cat(sprintf("✓ Input: %s cells, %s genes\n", 
            format(ncol(seurat_obj), big.mark = ","),
            format(nrow(seurat_obj), big.mark = ",")))
cat(sprintf("✓ %s levels: %d\n", MANUAL_COL, length(unique(seurat_obj[[MANUAL_COL, drop = TRUE]]))))
cat(sprintf("✓ %s levels: %d\n\n", BATCH_COL, n_batches))

# =========================
# Analysis Strategy
# =========================
if (is.null(TARGET_MANUAL)) {
  cat("=== Analysis Strategy: AUTO-DETECT ===\n")
  
  available_types <- unique(seurat_obj[[MANUAL_COL, drop = TRUE]])
  available_types <- available_types[!is.na(available_types) & available_types != ""]
  
  excluded_found <- intersect(available_types, EXCLUDE_CELLTYPES)
  if (length(excluded_found) > 0) {
    cat(sprintf("⚠ Excluding %d problematic types:\n", length(excluded_found)))
    for (ex in excluded_found) {
      n_cells_ex <- sum(seurat_obj[[MANUAL_COL, drop = TRUE]] == ex, na.rm = TRUE)
      cat(sprintf("  - %s: %s cells (excluded)\n", ex, format(n_cells_ex, big.mark = ",")))
    }
    available_types <- setdiff(available_types, EXCLUDE_CELLTYPES)
  }
  
  if (length(available_types) == 0) {
    stop("❌ No valid types after exclusion")
  }
  
  cell_types_to_analyze <- as.list(setNames(available_types, available_types))
  
  cat(sprintf("\nWill analyze %d types:\n", length(cell_types_to_analyze)))
  for (ct in names(cell_types_to_analyze)) {
    n_cells_ct <- sum(seurat_obj[[MANUAL_COL, drop = TRUE]] == ct, na.rm = TRUE)
    cat(sprintf("  - %s: %s cells\n", ct, format(n_cells_ct, big.mark = ",")))
  }
  cat("\n")
  
} else {
  cat("=== Analysis Strategy: USER-SPECIFIED ===\n")
  
  available_types <- unique(seurat_obj[[MANUAL_COL, drop = TRUE]])
  invalid_types <- setdiff(TARGET_MANUAL, available_types)
  if (length(invalid_types) > 0) {
    stop(sprintf("❌ Invalid types: %s", paste(invalid_types, collapse = ", ")))
  }
  
  cell_types_to_analyze <- as.list(setNames(TARGET_MANUAL, TARGET_MANUAL))
  cat(sprintf("Will analyze %d types:\n", length(cell_types_to_analyze)))
  for (ct in names(cell_types_to_analyze)) {
    n_cells_ct <- sum(seurat_obj[[MANUAL_COL, drop = TRUE]] == ct, na.rm = TRUE)
    cat(sprintf("  - %s: %s cells\n", ct, format(n_cells_ct, big.mark = ",")))
  }
  cat("\n")
}

config_list <- list(
  input_rds = INPUT_RDS,
  output_dir = OUTPUT_DIR,
  manual_col = MANUAL_COL,
  batch_col = BATCH_COL,
  target_manual = TARGET_MANUAL,
  excluded_celltypes = EXCLUDE_CELLTYPES,
  auto_detected_types = if(is.null(TARGET_MANUAL)) names(cell_types_to_analyze) else NULL,
  min_cells_per_celltype = MIN_CELLS_PER_CELLTYPE,
  min_cells_per_batch = MIN_CELLS_PER_BATCH,
  n_var_features = N_VAR_FEATURES,
  n_pcs = N_PCS,
  cluster_res = CLUSTER_RES,
  cluster_algo = CLUSTER_ALGO,
  random_seed = RANDOM_SEED,
  save_pre_harmony = SAVE_PRE_HARMONY
)
saveRDS(config_list, file.path(OUTPUT_DIR, "analysis_config.rds"))

# =========================
# Helper Functions
# =========================

#' P0-3 FIX: Strict Harmony with required params fallback
run_harmony_safe <- function(obj, batch_col, n_pcs_use) {
  cat("Running Harmony...\n")
  
  fn <- tryCatch(
    getS3method("RunHarmony", "Seurat"),
    error = function(e) {
      cat("  Cannot get S3 method, using generic\n")
      harmony::RunHarmony
    }
  )
  
  fml <- names(formals(fn))
  cat(sprintf("  Parameters: %s\n", paste(head(fml, 10), collapse = ", ")))
  
  params <- list(object = obj, group.by.vars = batch_col)
  
  if ("reduction.use" %in% fml) {
    params$reduction.use <- "pca"
    cat("  Using: reduction.use = 'pca'\n")
  } else if ("reduction" %in% fml) {
    params$reduction <- "pca"
    cat("  Using: reduction = 'pca'\n")
  } else {
    stop("❌ Cannot determine reduction parameter")
  }
  
  if ("dims.use" %in% fml) {
    params$dims.use <- 1:n_pcs_use
    cat(sprintf("  Using: dims.use = 1:%d\n", n_pcs_use))
  } else if ("dims" %in% fml) {
    params$dims <- 1:n_pcs_use
    cat(sprintf("  Using: dims = 1:%d\n", n_pcs_use))
  } else {
    stop("❌ Cannot determine dims parameter")
  }
  
  if ("max_iter" %in% fml || "max.iter" %in% fml) {
    param_name <- if ("max_iter" %in% fml) "max_iter" else "max.iter"
    params[[param_name]] <- 20
  }
  
  if ("early_stop" %in% fml || "early.stop" %in% fml) {
    param_name <- if ("early_stop" %in% fml) "early_stop" else "early.stop"
    params[[param_name]] <- TRUE
  }
  
  if ("verbose" %in% fml) {
    params$verbose <- FALSE
  }
  
  result <- tryCatch({
    do.call(fn, params)
  }, error = function(e) {
    cat(sprintf("❌ Harmony failed: %s\n", e$message))
    cat("  Retrying with required params only...\n")
    
    required_names <- intersect(names(params), 
                                c("object", "group.by.vars", "reduction", "reduction.use", 
                                  "dims", "dims.use"))
    fallback_params <- params[required_names]
    
    tryCatch({
      do.call(fn, fallback_params)
    }, error = function(e2) {
      stop(sprintf("❌ Harmony failed with required params. Error: %s", e2$message))
    })
  })
  
  cat("✓ Harmony completed\n")
  return(result)
}

#' P1-1 FIX: Improved clustering algorithm test
test_clustering_algorithm <- function(test_obj, requested_algo) {
  cat(sprintf("Testing algorithm %d...\n", requested_algo))
  
  set.seed(42)
  n_test <- min(50, ncol(test_obj))
  test_cells <- sample(colnames(test_obj), n_test)
  test_small <- test_obj[, test_cells]
  
  test_small <- NormalizeData(test_small, verbose = FALSE)
  test_small <- FindVariableFeatures(test_small, nfeatures = 100, verbose = FALSE)
  test_small <- ScaleData(test_small, features = VariableFeatures(test_small), verbose = FALSE)
  test_small <- RunPCA(test_small, npcs = 5, verbose = FALSE)
  test_small <- FindNeighbors(test_small, dims = 1:5, verbose = FALSE)
  
  result <- tryCatch({
    test_cluster <- FindClusters(test_small, resolution = 0.5, 
                                  algorithm = requested_algo, verbose = FALSE)
    cat(sprintf("✓ Algorithm %d available\n", requested_algo))
    requested_algo
  }, error = function(e) {
    cat(sprintf("⚠ Algorithm %d failed: %s\n", requested_algo, e$message))
    cat("  Falling back to algorithm 1 (Louvain)\n")
    1
  })
  
  rm(test_small)
  invisible(gc())
  return(result)
}
#' Harmony sanity check - must run AFTER RunPCA and BEFORE RunHarmony
harmony_sanity_check <- function(seurat_sub, batch_col, celltype_name,
                                 min_cells_per_batch_in_celltype = 20) {
  cat("\n--- Harmony Sanity Check ---\n")

  # 0) Basic existence checks
  if (!batch_col %in% colnames(seurat_sub@meta.data)) {
    stop(sprintf("❌ batch_col '%s' not found in meta.data (celltype=%s)", batch_col, celltype_name))
  }
  if (!"pca" %in% names(seurat_sub@reductions)) {
    stop(sprintf("❌ PCA reduction not found before Harmony (celltype=%s). RunPCA first.", celltype_name))
  }

  # 1) Clean batch labels
  seurat_sub@meta.data[[batch_col]] <- trimws(as.character(seurat_sub@meta.data[[batch_col]]))
  bad_batch <- is.na(seurat_sub@meta.data[[batch_col]]) | seurat_sub@meta.data[[batch_col]] == ""
  if (any(bad_batch)) {
    cat(sprintf("⚠ Removing %d cells with NA/empty '%s'\n", sum(bad_batch), batch_col))
    keep_cells <- rownames(seurat_sub@meta.data)[!bad_batch]
    seurat_sub <- subset(seurat_sub, cells = keep_cells)
  }

  # 2) Drop tiny batches (stability)
  batch_tab <- table(seurat_sub@meta.data[[batch_col]])
  cat("Batch table (post-clean):\n")
  print(batch_tab)

  tiny_batches <- names(batch_tab)[batch_tab < min_cells_per_batch_in_celltype]
  if (length(tiny_batches) > 0) {
    cat(sprintf("⚠ Dropping %d tiny batches (<%d cells): %s\n",
                length(tiny_batches), min_cells_per_batch_in_celltype,
                paste(tiny_batches, collapse = ", ")))
    keep_cells <- rownames(seurat_sub@meta.data)[!(seurat_sub@meta.data[[batch_col]] %in% tiny_batches)]
    seurat_sub <- subset(seurat_sub, cells = keep_cells)
  }

  # Must still have >= 2 batches for Harmony
  nb <- length(unique(seurat_sub@meta.data[[batch_col]]))
  if (nb < 2) {
    stop(sprintf("❌ After cleaning, only %d batch left. Skip Harmony for celltype=%s", nb, celltype_name))
  }

  # 3) Align PCA embedding rows to cell order (prevents non-conformable arguments)
  emb <- Embeddings(seurat_sub, "pca")
  cells <- colnames(seurat_sub)

  if (!all(cells %in% rownames(emb))) {
    stop(sprintf("❌ PCA embeddings missing cells (celltype=%s).", celltype_name))
  }
  seurat_sub[["pca"]]@cell.embeddings <- emb[cells, , drop = FALSE]

  cat(sprintf("✓ Embeddings aligned: %d cells × %d PCs\n",
              nrow(seurat_sub[["pca"]]@cell.embeddings),
              ncol(seurat_sub[["pca"]]@cell.embeddings)))
  cat("✓ Harmony sanity check done.\n")

  return(seurat_sub)
}

#' Clean legacy data
clean_legacy_data <- function(obj) {
  cat("Cleaning legacy reductions/graphs...\n")
  
  obj@graphs <- list()
  obj@neighbors <- list()
  
  reductions_to_remove <- c("pca", "umap", "umap_harmony", "harmony", "umap_pca")
  removed <- intersect(names(obj@reductions), reductions_to_remove)
  if (length(removed) > 0) {
    cat(sprintf("  Removing: %s\n", paste(removed, collapse = ", ")))
    obj@reductions <- obj@reductions[!names(obj@reductions) %in% reductions_to_remove]
  }
  
  tryCatch({
    if ("scale.data" %in% names(obj[["RNA"]]@layers)) {
      obj[["RNA"]]@layers[["scale.data"]] <- NULL
    }
  }, error = function(e) {})
  
  obj
}

#' FindAllMarkers compatibility
find_markers_safe <- function(obj) {
  fml <- names(formals(Seurat::FindAllMarkers))
  
  args <- list(
    object = obj,
    only.pos = TRUE,
    min.pct = MARKER_MIN_PCT,
    logfc.threshold = MARKER_LOGFC,
    test.use = "wilcox",
    max.cells.per.ident = MAX_CELLS_PER_IDENT,
    assay = DefaultAssay(obj),
    verbose = FALSE
  )
  
  if ("layer" %in% fml) {
    args$layer <- "data"
  } else if ("slot" %in% fml) {
    args$slot <- "data"
  }
  
  do.call(Seurat::FindAllMarkers, args)
}

#' Heatmap gene selection
select_heatmap_genes <- function(markers, top_n) {
  MIN_PER_CLUSTER <- max(3, floor(top_n * 0.3))
  
  all_genes <- c()
  clusters <- unique(markers$cluster)
  
  for (clust in clusters) {
    clust_markers <- markers %>%
      filter(cluster == clust) %>%
      arrange(desc(avg_log2FC)) %>%
      head(MIN_PER_CLUSTER) %>%
      pull(gene)
    all_genes <- c(all_genes, clust_markers)
  }
  all_genes <- unique(all_genes)
  
  if (length(all_genes) < top_n * length(clusters)) {
    remaining <- markers %>%
      filter(!gene %in% all_genes) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(avg_log2FC)) %>%
      pull(gene) %>%
      unique()
    
    n_needed <- (top_n * length(clusters)) - length(all_genes)
    all_genes <- c(all_genes, head(remaining, n_needed))
  }
  
  unique(all_genes)
}

#' Stratified sampling
stratified_sample_for_heatmap <- function(obj, max_cells, group_by) {
  meta <- obj@meta.data
  groups <- unique(meta[[group_by]])
  
  n_groups <- length(groups)
  per_group <- ceiling(max_cells / n_groups)
  
  cells_keep <- c()
  for (grp in groups) {
    cells_grp <- rownames(meta)[meta[[group_by]] == grp]
    n_take <- min(length(cells_grp), per_group)
    if (n_take > 0) {
      cells_keep <- c(cells_keep, sample(cells_grp, n_take))
    }
  }
  
  if (length(cells_keep) > max_cells) {
    cells_keep <- sample(cells_keep, max_cells)
  }
  
  cat(sprintf("  Stratified: %d cells from %d groups\n", length(cells_keep), n_groups))
  cells_keep
}

#' P2-2 FIX: Cluster × Sample composition table
compute_cluster_sample_composition <- function(obj, cluster_col, sample_col, output_prefix) {
  cat("Computing cluster × sample composition...\n")
  
  meta <- obj@meta.data
  comp_table <- table(meta[[cluster_col]], meta[[sample_col]])
  
  comp_df <- as.data.frame.matrix(comp_table)
  comp_df <- tibble::rownames_to_column(comp_df, var = "cluster")
  
  write.csv(comp_df, 
            file.path(TAB_DIR, sprintf("%s_cluster_sample_counts.csv", output_prefix)),
            row.names = FALSE)
  
  comp_prop <- prop.table(comp_table, margin = 1)
  comp_prop_df <- as.data.frame.matrix(comp_prop)
  comp_prop_df <- tibble::rownames_to_column(comp_prop_df, var = "cluster")
  
  write.csv(comp_prop_df,
            file.path(TAB_DIR, sprintf("%s_cluster_sample_proportions.csv", output_prefix)),
            row.names = FALSE)
  
  cat(sprintf("✓ Saved cluster × sample tables\n"))
  
  invisible(list(counts = comp_df, proportions = comp_prop_df))
}

#' Process one cell type
process_subset <- function(seurat_obj, celltype_name, celltype_filter, output_prefix, 
                          actual_cluster_algo) {
  cat("\n========================================\n")
  cat(sprintf("Processing: %s\n", celltype_name))
  cat("========================================\n")
  
  cat(sprintf("Subsetting to '%s'...\n", celltype_filter))
  meta <- seurat_obj@meta.data
  cells_keep <- rownames(meta)[meta[[MANUAL_COL]] == celltype_filter]
  
  if (length(cells_keep) == 0) {
    cat("❌ No cells found. Skipping.\n")
    return(NULL)
  }
  
  seurat_sub <- subset(seurat_obj, cells = cells_keep)
  n_cells_total <- ncol(seurat_sub)
  cat(sprintf("Cells: %s\n", format(n_cells_total, big.mark = ",")))
  
  # P1-2 FIX: Check minimum cells
  if (n_cells_total < MIN_CELLS_PER_CELLTYPE) {
    cat(sprintf("❌ Only %d cells (< %d threshold). Skipping.\n", 
                n_cells_total, MIN_CELLS_PER_CELLTYPE))
    return(NULL)
  }
  
  batch_counts <- table(seurat_sub@meta.data[[BATCH_COL]])
  n_batches_sub <- length(batch_counts)
  
  if (n_batches_sub < 2) {
    cat("❌ Only 1 batch. Skipping.\n")
    return(NULL)
  }
  
  # P1-2 FIX: Check minimum cells per batch
  min_batch_size <- min(batch_counts)
  if (min_batch_size < MIN_CELLS_PER_BATCH) {
    cat(sprintf("⚠ WARNING: Smallest batch has only %d cells (< %d threshold)\n",
                min_batch_size, MIN_CELLS_PER_BATCH))
    cat("  This may cause unstable Harmony/clustering results\n")
    cat("  Consider increasing MIN_CELLS_PER_BATCH or filtering small batches\n")
  }
  
  cat(sprintf("Batches: %d\n", n_batches_sub))
  batch_summary <- data.frame(
    batch = names(batch_counts),
    n_cells = as.integer(batch_counts),
    stringsAsFactors = FALSE
  )
  print(batch_summary)
  
  seurat_sub <- clean_legacy_data(seurat_sub)
  
  # Normalize + HVG + Scale + PCA
  cat("\n--- Normalize + HVG + Scale + PCA ---\n")
  seurat_sub <- NormalizeData(seurat_sub, verbose = FALSE)
  seurat_sub <- FindVariableFeatures(seurat_sub, nfeatures = N_VAR_FEATURES, verbose = FALSE)
  cat(sprintf("✓ HVG: %d\n", length(VariableFeatures(seurat_sub))))
  
  seurat_sub <- ScaleData(seurat_sub, features = VariableFeatures(seurat_sub), verbose = FALSE)
  
  set.seed(RANDOM_SEED)
  seurat_sub <- RunPCA(seurat_sub, features = VariableFeatures(seurat_sub), 
                       npcs = N_PCS, verbose = FALSE)
  
  n_pcs_actual <- ncol(Embeddings(seurat_sub, "pca"))
  n_pcs_use <- min(N_PCS, n_pcs_actual)
  cat(sprintf("✓ PCA: %d PCs (use %d)\n", n_pcs_actual, n_pcs_use))
  
  if (n_pcs_use < 10) {
    cat(sprintf("⚠ Only %d PCs available\n", n_pcs_use))
  }
  
  tryCatch({
    if ("scale.data" %in% names(seurat_sub[["RNA"]]@layers)) {
      seurat_sub[["RNA"]]@layers[["scale.data"]] <- NULL
      invisible(gc())
    }
  }, error = function(e) {})
  
  # Pre-Harmony
  if (SAVE_PRE_HARMONY) {
    cat("\n--- Pre-Harmony (PCA) ---\n")
    set.seed(RANDOM_SEED)
    seurat_sub <- RunUMAP(seurat_sub, reduction = "pca", dims = 1:n_pcs_use,
                          reduction.name = "umap_pca", reduction.key = "umappca_",
                          seed.use = RANDOM_SEED, verbose = FALSE)
    
    # P2-1 FIX: Quick pre-Harmony clustering sanity check
    seurat_sub <- tryCatch({
      FindNeighbors(seurat_sub, reduction = "pca", dims = 1:n_pcs_use,
                    graph.name = c("RNA_nn_pca", "RNA_snn_pca"), verbose = FALSE)
    }, error = function(e) {
      FindNeighbors(seurat_sub, reduction = "pca", dims = 1:n_pcs_use, verbose = FALSE)
    })

    set.seed(RANDOM_SEED)
    seurat_sub <- FindClusters(seurat_sub, resolution = 0.5, algorithm = 1, 
                                graph.name = "RNA_snn_pca", verbose = FALSE)
    seurat_sub$clusters_pca <- Idents(seurat_sub)
    
    # P0-1 FIX: Cowplot dependency check
    pdf(file.path(FIG_DIR, sprintf("%s_pre_harmony_umap.pdf", output_prefix)), 
        width = 12, height = 5)
    
    p1 <- DimPlot(seurat_sub, reduction = "umap_pca", group.by = BATCH_COL, raster = TRUE) +
      ggtitle(sprintf("%s - Pre-Harmony | Batch", celltype_name))
    p2 <- DimPlot(seurat_sub, reduction = "umap_pca", group.by = "clusters_pca",
                  label = TRUE, repel = TRUE, raster = TRUE) +
      ggtitle(sprintf("%s - Pre-Harmony | PCA Clusters", celltype_name))
    
    if (requireNamespace("cowplot", quietly = TRUE)) {
      print(cowplot::plot_grid(p1, p2, ncol = 2))
    } else {
      cat("  ⚠ cowplot not available, using separate pages\n")
      print(p1)
      print(p2)
    }
    dev.off()
    
    cat("✓ Pre-Harmony UMAP saved\n")
  }
  seurat_sub <- harmony_sanity_check(
  seurat_sub,
  batch_col = BATCH_COL,
  celltype_name = celltype_name,
  min_cells_per_batch_in_celltype = MIN_CELLS_PER_BATCH  # 你上面配置的阈值
  )

  cat("\n--- Harmony ---\n")
  seurat_sub <- run_harmony_safe(seurat_sub, BATCH_COL, n_pcs_use)
  
  # UMAP + Clustering
  cat("\n--- UMAP + Clustering ---\n")
  set.seed(RANDOM_SEED)
  seurat_sub <- RunUMAP(seurat_sub, reduction = "harmony", dims = 1:n_pcs_use,
                        reduction.name = "umap_harmony", reduction.key = "umaph_",
                        seed.use = RANDOM_SEED, verbose = FALSE)
  
  seurat_sub <- FindNeighbors(seurat_sub, reduction = "harmony", dims = 1:n_pcs_use, 
                               verbose = FALSE)
  
  set.seed(RANDOM_SEED)
  fc_args <- list(
    object = seurat_sub,
    resolution = CLUSTER_RES,
    algorithm = actual_cluster_algo,
    verbose = FALSE
  )
  
  fml_fc <- names(formals(Seurat::FindClusters))
  if ("random.seed" %in% fml_fc) {
    fc_args$random.seed <- RANDOM_SEED
  } else if ("seed.use" %in% fml_fc) {
    fc_args$seed.use <- RANDOM_SEED
  }
  
  seurat_sub <- do.call(Seurat::FindClusters, fc_args)
  seurat_sub$leiden_harmony <- Idents(seurat_sub)
  
  algo_name <- if (actual_cluster_algo == 4) "Leiden" else "Louvain"
  cat(sprintf("✓ %d clusters (%s)\n", length(unique(seurat_sub$leiden_harmony)), algo_name))
  
  # P2-2 FIX: Cluster × Sample composition
  compute_cluster_sample_composition(seurat_sub, "leiden_harmony", BATCH_COL, output_prefix)
  
  # Save
  sub_rds <- file.path(OBJ_DIR, sprintf("%s_harmony_clustered.rds", output_prefix))
  saveRDS(seurat_sub, sub_rds, compress = "xz")
  cat(sprintf("✓ Saved: %s\n", sub_rds))
  
  # UMAP plots
  cat("\n--- UMAP Plots ---\n")
  pdf(file.path(FIG_DIR, sprintf("%s_umap_harmony.pdf", output_prefix)), width = 15, height = 5)
  p1 <- DimPlot(seurat_sub, reduction = "umap_harmony", group.by = "leiden_harmony", 
                label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle(sprintf("%s - Harmony | Clusters", celltype_name))
  p2 <- DimPlot(seurat_sub, reduction = "umap_harmony", group.by = BATCH_COL, raster = TRUE) +
    ggtitle(sprintf("%s - Harmony | Batch", celltype_name))
  p3 <- DimPlot(seurat_sub, reduction = "umap_harmony", group.by = MANUAL_COL, 
                label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle(sprintf("%s - Harmony | Annotation", celltype_name))
  
  if (requireNamespace("cowplot", quietly = TRUE)) {
    print(cowplot::plot_grid(p1, p2, p3, ncol = 3))
  } else {
    print(p1); print(p2); print(p3)
  }
  dev.off()
  
  # Markers
  cat("\n--- Markers ---\n")
  Idents(seurat_sub) <- "leiden_harmony"
  markers <- find_markers_safe(seurat_sub)
  
  write.csv(markers, file.path(TAB_DIR, sprintf("%s_markers.csv", output_prefix)), 
            row.names = FALSE)
  cat(sprintf("✓ %d markers\n", nrow(markers)))
  
  # Heatmap
  cat("\n--- Heatmap ---\n")
  if (nrow(markers) > 0) {
    hm_genes <- select_heatmap_genes(markers, TOP_N_HEATMAP)
    hm_genes <- hm_genes[hm_genes %in% rownames(seurat_sub)]
    
    if (length(hm_genes) >= 2) {
      n_cells_hm <- min(ncol(seurat_sub), 10000)
      if (ncol(seurat_sub) > n_cells_hm) {
        set.seed(RANDOM_SEED)
        cells_hm <- stratified_sample_for_heatmap(seurat_sub, n_cells_hm, "leiden_harmony")
        seurat_sub_hm <- subset(seurat_sub, cells = cells_hm)
        
        if (inherits(seurat_sub_hm[["RNA"]], "Assay5")) {
          cat("  Joining layers (Assay5)\n")
          seurat_sub_hm[["RNA"]] <- JoinLayers(seurat_sub_hm[["RNA"]])
        }
      } else {
        seurat_sub_hm <- seurat_sub
      }
      
      tryCatch({
        if ("scale.data" %in% names(seurat_sub_hm[["RNA"]]@layers)) {
          seurat_sub_hm[["RNA"]]@layers[["scale.data"]] <- NULL
        }
      }, error = function(e) {})
      
      seurat_sub_hm <- ScaleData(seurat_sub_hm, features = hm_genes, verbose = FALSE)
      
      pdf(file.path(FIG_DIR, sprintf("%s_heatmap.pdf", output_prefix)), width = 12, height = 10)
      print(DoHeatmap(seurat_sub_hm, features = hm_genes, group.by = "leiden_harmony", 
                      raster = TRUE) +
              ggtitle(sprintf("%s - Top Markers", celltype_name)))
      dev.off()
      
      rm(seurat_sub_hm)
      invisible(gc())
    }
  }
  
  # ROGUE
  if (RUN_ROGUE) {
    cat("\n--- ROGUE ---\n")
    rogue_cluster <- compute_rogue_by_combo(
      seurat_sub, BATCH_COL, "leiden_harmony",
      MIN_CELLS_PER_GROUP, ROGUE_MIN_CELLS_GENE, ROGUE_MIN_GENES_CELL,
      ROGUE_MAX_CELLS_GROUP, ROGUE_MAX_GENES_GROUP, ROGUE_MAX_COMBOS
    )
    save_rogue_outputs(rogue_cluster, sprintf("%s_clusters", output_prefix),
                       paste0(celltype_name, " | leiden_harmony | "))
  }
  
  cat(sprintf("\n✓ Completed: %s\n", celltype_name))
  invisible(gc())
  return(sub_rds)
}

#' P0-2 FIX: ROGUE with correct skipped combo export
compute_rogue_by_combo <- function(seurat_obj, dataset_col, group_col,
                                    min_cells_per_group, min_cells_gene, min_genes_cell,
                                    max_cells_group, max_genes_group, max_combos) {
  if (!requireNamespace("ROGUE", quietly = TRUE)) {
    stop("ROGUE not installed")
  }
  
  meta <- seurat_obj@meta.data
  
  combo_counts_full <- meta %>%
    dplyr::group_by(.data[[dataset_col]], .data[[group_col]]) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_cells >= min_cells_per_group) %>%
    dplyr::arrange(desc(n_cells))
  
  cat(sprintf("ROGUE combos (%s × %s) n≥%d: %d total\n",
              dataset_col, group_col, min_cells_per_group, nrow(combo_counts_full)))
  
  if (nrow(combo_counts_full) == 0) {
    cat("⚠ No valid combos\n")
    return(data.frame())
  }
  
  # P0-2 FIX: Correct skipped combo export logic
  combo_counts <- combo_counts_full
  if (nrow(combo_counts_full) > max_combos) {
    cat(sprintf("⚠ Limiting to top %d combos (by cell count)\n", max_combos))
    
    combo_keep <- combo_counts_full[1:max_combos, , drop = FALSE]
    combo_skipped <- combo_counts_full[(max_combos + 1):nrow(combo_counts_full), , drop = FALSE]
    
    write.csv(combo_keep,
              file.path(TAB_DIR, sprintf("%s_%s_rogue_combos_kept.csv", dataset_col, group_col)),
              row.names = FALSE)
    write.csv(combo_skipped,
              file.path(TAB_DIR, sprintf("%s_%s_rogue_combos_skipped.csv", dataset_col, group_col)),
              row.names = FALSE)
    
    cat(sprintf("  Kept: %s\n", 
                file.path(TAB_DIR, sprintf("%s_%s_rogue_combos_kept.csv", dataset_col, group_col))))
    cat(sprintf("  Skipped: %s\n",
                file.path(TAB_DIR, sprintf("%s_%s_rogue_combos_skipped.csv", dataset_col, group_col))))
    
    combo_counts <- combo_keep
  }
  
  expr_all <- tryCatch({
    GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts")
  }, error = function(e) {
    GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), slot = "counts")
  })
  
  if (is.null(rownames(expr_all))) rownames(expr_all) <- rownames(seurat_obj)
  if (is.null(colnames(expr_all))) colnames(expr_all) <- colnames(seurat_obj)
  
  rogue_list <- vector("list", nrow(combo_counts))
  pb <- txtProgressBar(min = 0, max = nrow(combo_counts), style = 3, file = stderr())
  
  for (i in seq_len(nrow(combo_counts))) {
    ds  <- combo_counts[[dataset_col]][i]
    grp <- combo_counts[[group_col]][i]
    n0  <- combo_counts$n_cells[i]
    
    cells_keep <- rownames(meta)[meta[[dataset_col]] == ds & meta[[group_col]] == grp]
    
    rogue_val <- NA_real_
    n_genes_used <- NA_integer_
    n_cells_used <- length(cells_keep)
    
    tryCatch({
      if (length(cells_keep) > max_cells_group) {
        set.seed(RANDOM_SEED + i)
        cells_keep <- sample(cells_keep, max_cells_group)
      }
      
      mat <- expr_all[, cells_keep, drop = FALSE]
      
      gene_detect <- Matrix::rowSums(mat > 0)
      keep_genes <- which(gene_detect >= min_cells_gene)
      mat <- mat[keep_genes, , drop = FALSE]
      
      cell_detect <- Matrix::colSums(mat > 0)
      keep_cells <- which(cell_detect >= min_genes_cell)
      mat <- mat[, keep_cells, drop = FALSE]
      
      if (ncol(mat) < min_cells_per_group || nrow(mat) < 100) {
        rogue_val <- NA_real_
      } else {
        if (nrow(mat) > max_genes_group) {
          vf <- VariableFeatures(seurat_obj)
          vf <- vf[vf %in% rownames(mat)]
          if (length(vf) >= 200) {
            mat <- mat[head(vf, max_genes_group), , drop = FALSE]
          } else {
            rs <- Matrix::rowSums(mat)
            topg <- names(sort(rs, decreasing = TRUE))[1:min(max_genes_group, nrow(mat))]
            mat <- mat[topg, , drop = FALSE]
          }
        }
        
        n_genes_used <- nrow(mat)
        n_cells_used <- ncol(mat)
        
        dense <- as.matrix(mat)
        ent_res <- ROGUE::SE_fun(dense)
        
        if (any(!is.finite(ent_res$entropy))) {
          rogue_val <- NA_real_
        } else {
          rogue_val <- ROGUE::CalculateRogue(ent_res, platform = "UMI")
        }
      }
      
    }, error = function(e) {
      rogue_val <<- NA_real_
    })
    
    rogue_list[[i]] <- data.frame(
      dataset = as.character(ds),
      group = as.character(grp),
      n_cells = as.integer(n0),
      n_cells_used = as.integer(n_cells_used),
      n_genes_used = as.integer(n_genes_used),
      rogue_value = as.numeric(rogue_val),
      stringsAsFactors = FALSE
    )
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat("\n⚠ ROGUE Interpretation:\n")
  cat("  - Scores for INTRA-celltype cluster purity comparison\n")
  cat("  - DO NOT compare across different cell types\n")
  cat("  - Each celltype has different HVG/filtering\n\n")
  
  dplyr::bind_rows(rogue_list)
}

#' Save ROGUE outputs
save_rogue_outputs <- function(rogue_df, out_prefix, title_prefix = "") {
  if (nrow(rogue_df) == 0) {
    cat("⚠ Empty ROGUE\n")
    return(invisible(NULL))
  }
  
  write.csv(rogue_df, file.path(TAB_DIR, paste0(out_prefix, "_rogue_long.csv")), 
            row.names = FALSE)
  
  rogue_wide <- rogue_df %>%
    select(dataset, group, rogue_value) %>%
    pivot_wider(names_from = group, values_from = rogue_value)
  write.csv(rogue_wide, file.path(TAB_DIR, paste0(out_prefix, "_rogue_wide.csv")), 
            row.names = FALSE)
  
  by_group <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    group_by(group) %>%
    summarise(n_datasets = n(), total_cells = sum(n_cells),
              median_rogue = median(rogue_value), mean_rogue = mean(rogue_value),
              sd_rogue = sd(rogue_value), .groups = "drop") %>%
    arrange(desc(median_rogue))
  write.csv(by_group, file.path(TAB_DIR, paste0(out_prefix, "_rogue_summary_by_group.csv")), 
            row.names = FALSE)
  
  by_dataset <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    group_by(dataset) %>%
    summarise(n_groups = n(), total_cells = sum(n_cells),
              median_rogue = median(rogue_value), mean_rogue = mean(rogue_value),
              .groups = "drop") %>%
    arrange(desc(median_rogue))
  write.csv(by_dataset, file.path(TAB_DIR, paste0(out_prefix, "_rogue_summary_by_dataset.csv")), 
            row.names = FALSE)
  
  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_boxplot.pdf")), width = 12, height = 6)
  p1 <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    ggplot(aes(x = reorder(group, rogue_value, median), y = rogue_value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.8) +
    geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red") +
    labs(title = paste0(title_prefix, "ROGUE by Group"), x = "Group", y = "ROGUE") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1)
  dev.off()
  
  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_heatmap.pdf")), width = 14, height = 10)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    mat <- rogue_df %>%
      select(dataset, group, rogue_value) %>%
      pivot_wider(names_from = group, values_from = rogue_value) %>%
      tibble::column_to_rownames("dataset") %>%
      as.matrix()
    pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, na_col = "grey90",
                       main = paste0(title_prefix, "ROGUE Heatmap"), fontsize = 8)
  } else {
    plot.new()
    text(0.5, 0.5, "pheatmap not installed", cex = 1.5)
  }
  dev.off()
  
  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_dotplot.pdf")), width = 14, height = 8)
  p3 <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    ggplot(aes(x = group, y = dataset, size = n_cells, color = rogue_value)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 10)) +
    labs(title = paste0(title_prefix, "ROGUE: Group × Dataset"),
         x = "Group", y = "Dataset", size = "Cells", color = "ROGUE") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  tryCatch({
    p3 <- p3 + scale_color_viridis_c(option = "viridis")
  }, error = function(e) {
    p3 <<- p3 + scale_color_gradient(low = "blue", high = "yellow")
  })
  print(p3)
  dev.off()
  
  saveRDS(list(long = rogue_df, wide = rogue_wide, 
               summary_group = by_group, summary_dataset = by_dataset),
          file.path(OBJ_DIR, paste0(out_prefix, "_rogue_results.rds")))
  
  cat(sprintf("✓ ROGUE saved: %s\n", out_prefix))
}

# =========================
# Main Loop
# =========================
cat("\n========================================\n")
cat("Main Analysis Loop\n")
cat("========================================\n\n")

cat("=== Testing Clustering Algorithm ===\n")
actual_cluster_algo <- test_clustering_algorithm(seurat_obj, CLUSTER_ALGO)
if (actual_cluster_algo != CLUSTER_ALGO) {
  cat(sprintf("⚠ Will use algorithm %d instead of %d\n", 
              actual_cluster_algo, CLUSTER_ALGO))
}
cat("\n")

results_paths <- list()

for (celltype_name in names(cell_types_to_analyze)) {
  celltype_filter <- cell_types_to_analyze[[celltype_name]]
  output_prefix <- gsub("[^A-Za-z0-9_-]", "_", celltype_name)
  
  result_path <- tryCatch({
    process_subset(seurat_obj, celltype_name, celltype_filter, output_prefix, actual_cluster_algo)
  }, error = function(e) {
    cat(sprintf("\n❌ ERROR in %s: %s\n", celltype_name, e$message))
    print(e)
    NULL
  })
  
  results_paths[[celltype_name]] <- result_path
  invisible(gc())
}

# =========================
# Finalize
# =========================
cat("\n========================================\n")
cat("Analysis Complete\n")
cat("========================================\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("Output:", OUTPUT_DIR, "\n\n")

cat("=== Results Summary ===\n")
successful <- sum(sapply(results_paths, function(x) !is.null(x)))
cat(sprintf("Completed: %d / %d\n", successful, length(results_paths)))
for (ct_name in names(results_paths)) {
  if (!is.null(results_paths[[ct_name]])) {
    cat(sprintf("  ✓ %s: %s\n", ct_name, results_paths[[ct_name]]))
  } else {
    cat(sprintf("  ✗ %s: FAILED\n", ct_name))
  }
}
cat("\n")

cat("=== Interpretation Notes ===\n")
cat("1. Each cell type analyzed INDEPENDENTLY (separate HVG/PCA/Harmony)\n")
cat("2. Results NOT comparable across cell types\n")
cat("3. Use for intra-celltype substructure only\n")
cat("4. ROGUE compares cluster purity within same celltype only\n")
cat("5. Pre-Harmony UMAPs saved for batch effect assessment\n")
cat("6. Cluster × Sample composition tables exported for diagnostics\n\n")

if (actual_cluster_algo != CLUSTER_ALGO) {
  cat(sprintf("⚠ Clustering used algorithm %d (not %d)\n", 
              actual_cluster_algo, CLUSTER_ALGO))
}

cat("\n=== Session Info ===\n")
print(sessionInfo())

cat("\n✓ Done. Check", LOG_FILE, "\n")
