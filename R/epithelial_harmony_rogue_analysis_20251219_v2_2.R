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
  library(jsonlite)
})

# =========================
# Configuration
# =========================
INPUT_RDS  <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/R/1218/harmony_subset_rogue_v2_2"

MANUAL_COL <- "Manual_Annotation"
BATCH_COL  <- "sample"

TARGET_MANUAL <- NULL
# TARGET_MANUAL <- c("Basal", "Goblet", "Ciliated", "Secretory")

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

RANDOM_SEED <- 42

# =========================
# IO
# =========================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_DIR   <- file.path(OUTPUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
TAB_DIR   <- file.path(OUTPUT_DIR, "tables");  dir.create(TAB_DIR, showWarnings = FALSE)
OBJ_DIR   <- file.path(OUTPUT_DIR, "objects"); dir.create(OBJ_DIR, showWarnings = FALSE)
LOG_FILE  <- file.path(OUTPUT_DIR, "analysis.log")

# P1-2 FIX: Add on.exit protection for sink
sink(LOG_FILE, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)

cat("========================================\n")
cat("Harmony + ROGUE Analysis (v2.2)\n")
cat("========================================\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("Input :", INPUT_RDS, "\n")
cat("Output:", OUTPUT_DIR, "\n\n")

# =========================
# Input Validation
# =========================
cat("=== Input Validation ===\n")
if (!file.exists(INPUT_RDS)) {
  stop(sprintf("❌ Input file not found: %s", INPUT_RDS))
}

seurat_obj <- readRDS(INPUT_RDS)

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("❌ RNA assay not found in Seurat object")
}

DefaultAssay(seurat_obj) <- "RNA"

if (!MANUAL_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("❌ '%s' not found in meta.data.", MANUAL_COL))
}

if (!BATCH_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("❌ '%s' not found in meta.data.", BATCH_COL))
}

# P0-2 FIX: Layer/slot compatible counts check
counts_available <- FALSE
tryCatch({
  test_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  if (!is.null(test_counts) && length(test_counts) > 0) {
    counts_available <- TRUE
    cat("✓ counts layer detected (v5 style)\n")
  }
}, error = function(e) {
  tryCatch({
    test_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    if (!is.null(test_counts) && length(test_counts) > 0) {
      counts_available <<- TRUE
      cat("✓ counts slot detected (v4 style)\n")
    }
  }, error = function(e2) {
    counts_available <<- FALSE
  })
})

if (!counts_available) {
  stop("❌ Cannot access counts data (tried both layer and slot)")
}

n_batches <- length(unique(seurat_obj@meta.data[[BATCH_COL]]))
if (n_batches < 2) {
  stop(sprintf("❌ Only %d batch found. Harmony requires ≥2 batches.", n_batches))
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
  cat("=== Analysis Strategy: ALL CELLS TOGETHER ===\n")
  cell_types_to_analyze <- list("AllCells" = NULL)
} else {
  cat("=== Analysis Strategy: PER CELL TYPE SEPARATELY ===\n")
  available_types <- unique(seurat_obj[[MANUAL_COL, drop = TRUE]])
  invalid_types <- setdiff(TARGET_MANUAL, available_types)
  if (length(invalid_types) > 0) {
    stop(sprintf("❌ Invalid cell types in TARGET_MANUAL: %s\nAvailable: %s",
                 paste(invalid_types, collapse = ", "),
                 paste(head(available_types, 20), collapse = ", ")))
  }
  
  cell_types_to_analyze <- as.list(setNames(TARGET_MANUAL, TARGET_MANUAL))
  cat(sprintf("Will analyze %d cell types separately:\n", length(cell_types_to_analyze)))
  cat(paste("-", names(cell_types_to_analyze), collapse = "\n"), "\n\n")
}

config_list <- list(
  input_rds = INPUT_RDS,
  output_dir = OUTPUT_DIR,
  manual_col = MANUAL_COL,
  batch_col = BATCH_COL,
  target_manual = TARGET_MANUAL,
  n_var_features = N_VAR_FEATURES,
  n_pcs = N_PCS,
  cluster_res = CLUSTER_RES,
  random_seed = RANDOM_SEED
)
saveRDS(config_list, file.path(OUTPUT_DIR, "analysis_config.rds"))

# =========================
# Helper Functions
# =========================

#' Harmony wrapper (Seurat) - strict args, no projection to avoid non-conformable
run_harmony_safe <- function(obj, batch_col, n_pcs_use) {
  # #region agent log
  log_entry <- jsonlite::toJSON(list(
    sessionId = "debug-session",
    runId = "post-fix",
    hypothesisId = "A",
    location = "run_harmony_safe:entry",
    message = "Harmony function entry",
    data = list(batch_col = batch_col, n_pcs_use = n_pcs_use, has_pca = "pca" %in% names(obj@reductions)),
    timestamp = as.numeric(Sys.time()) * 1000
  ), auto_unbox = TRUE)
  cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion
  
  cat("Running Harmony (strict Seurat args, project.dim=FALSE)...\n")
  
  # hard checks (fail-fast, prevents silent dimension mismatch)
  stopifnot("pca" %in% names(obj@reductions))
  pca_emb <- Seurat::Embeddings(obj, "pca")
  if (ncol(pca_emb) < 2) stop("❌ PCA embeddings have <2 PCs; cannot run Harmony")
  n_pcs_use <- min(n_pcs_use, ncol(pca_emb))
  
  # #region agent log
  pca_loadings <- tryCatch(Seurat::Loadings(obj, "pca"), error = function(e) NULL)
  log_entry <- jsonlite::toJSON(list(
    sessionId = "debug-session",
    runId = "post-fix",
    hypothesisId = "B",
    location = "run_harmony_safe:before_harmony",
    message = "PCA dimensions before Harmony",
    data = list(
      pca_emb_dim = paste(dim(pca_emb), collapse="x"),
      pca_loadings_dim = if(!is.null(pca_loadings)) paste(dim(pca_loadings), collapse="x") else "NULL",
      n_pcs_use = n_pcs_use
    ),
    timestamp = as.numeric(Sys.time()) * 1000
  ), auto_unbox = TRUE)
  cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion
  
  # IMPORTANT: project.dim=FALSE avoids post-harmony projection step that often triggers
  #            'non-conformable arguments' on some Seurat/Harmony combinations
  obj <- harmony::RunHarmony(
    object        = obj,
    group.by.vars = batch_col,
    reduction.use = "pca",
    dims.use      = 1:n_pcs_use,
    reduction.save= "harmony",
    project.dim   = FALSE,
    max_iter      = 20,
    early_stop    = TRUE,
    verbose       = TRUE
  )
  
  # #region agent log
  log_entry <- jsonlite::toJSON(list(
    sessionId = "debug-session",
    runId = "post-fix",
    hypothesisId = "B",
    location = "run_harmony_safe:after_harmony",
    message = "Harmony completed",
    data = list(
      has_harmony = "harmony" %in% names(obj@reductions),
      harmony_dim = if("harmony" %in% names(obj@reductions)) {
        paste(dim(Seurat::Embeddings(obj, "harmony")), collapse="x")
      } else "NULL"
    ),
    timestamp = as.numeric(Sys.time()) * 1000
  ), auto_unbox = TRUE)
  cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion
  
  # sanity check
  if (!"harmony" %in% names(obj@reductions)) {
    stop("❌ Harmony finished but 'harmony' reduction not found in object@reductions")
  }
  cat("✓ Harmony completed successfully\n")
  return(obj)
}

#' Clean legacy data
clean_legacy_data <- function(obj) {
  cat("Cleaning legacy reductions/graphs/neighbors...\n")
  
  obj@graphs <- list()
  obj@neighbors <- list()
  
  reductions_to_remove <- c("pca", "umap", "umap_harmony", "harmony")
  removed <- intersect(names(obj@reductions), reductions_to_remove)
  if (length(removed) > 0) {
    cat(sprintf("  Removing reductions: %s\n", paste(removed, collapse = ", ")))
    obj@reductions <- obj@reductions[!names(obj@reductions) %in% reductions_to_remove]
  }
  
  # Clear scale.data if exists
  tryCatch({
    if ("scale.data" %in% names(obj[["RNA"]]@layers)) {
      cat("  Removing old scale.data layer\n")
      obj[["RNA"]]@layers[["scale.data"]] <- NULL
    }
  }, error = function(e) {
    # Silently continue if layers access fails (v4 object)
  })
  
  obj
}

#' P0-4 FIX: FindAllMarkers with layer/slot compatibility
find_markers_safe <- function(obj, ...) {
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
  
  # Try layer parameter first (v5), fallback to slot (v4)
  if ("layer" %in% fml) {
    args$layer <- "data"
    cat("  Using layer='data' for markers\n")
  } else if ("slot" %in% fml) {
    args$slot <- "data"
    cat("  Using slot='data' for markers\n")
  }
  
  do.call(Seurat::FindAllMarkers, args)
}

#' Stratified sampling for heatmap
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
  
  cat(sprintf("  Stratified sampling: %d cells from %d groups\n", 
              length(cells_keep), n_groups))
  
  cells_keep
}

#' Process one cell type
process_subset <- function(seurat_obj, celltype_name, celltype_filter, output_prefix) {
  cat("\n========================================\n")
  cat(sprintf("Processing: %s\n", celltype_name))
  cat("========================================\n")
  
  # Subset cells
  if (is.null(celltype_filter)) {
    cat("Using all cells...\n")
    seurat_sub <- seurat_obj
  } else {
    cat(sprintf("Subsetting to '%s' cells...\n", celltype_filter))
    meta <- seurat_obj@meta.data
    cells_keep <- rownames(meta)[meta[[MANUAL_COL]] == celltype_filter]
    
    if (length(cells_keep) == 0) {
      cat(sprintf("❌ No cells found for '%s'. Skipping.\n", celltype_filter))
      return(NULL)
    }
    
    seurat_sub <- subset(seurat_obj, cells = cells_keep)
  }
  
  cat(sprintf("Cells after subset: %s\n", format(ncol(seurat_sub), big.mark = ",")))
  
  # Check batch diversity
  batch_counts <- table(seurat_sub@meta.data[[BATCH_COL]])
  n_batches_sub <- length(batch_counts)
  
  if (n_batches_sub < 2) {
    cat(sprintf("❌ Only %d batch. Cannot run Harmony. Skipping.\n", n_batches_sub))
    return(NULL)
  }
  
  cat(sprintf("Batch distribution: %d batches\n", n_batches_sub))
  batch_summary <- data.frame(
    batch = names(batch_counts),
    n_cells = as.integer(batch_counts),
    stringsAsFactors = FALSE
  )
  print(batch_summary)
  
  # Clean legacy data
  seurat_sub <- clean_legacy_data(seurat_sub)
  
  # Normalize + HVG + Scale + PCA
  cat("\n--- Normalize + HVG + Scale + PCA ---\n")
  seurat_sub <- NormalizeData(seurat_sub, verbose = FALSE)
  
  seurat_sub <- FindVariableFeatures(
    seurat_sub, 
    selection.method = "vst", 
    nfeatures = N_VAR_FEATURES, 
    verbose = FALSE
  )
  cat(sprintf("✓ HVG computed: %d\n", length(VariableFeatures(seurat_sub))))
  
  seurat_sub <- ScaleData(
    seurat_sub, 
    features = VariableFeatures(seurat_sub), 
    verbose = FALSE
  )
  
  set.seed(RANDOM_SEED)
  seurat_sub <- RunPCA(
    seurat_sub, 
    features = VariableFeatures(seurat_sub), 
    npcs = N_PCS, 
    verbose = FALSE
  )
  
  # P0-1 FIX: Get actual PC dimensions and bound n_pcs_use
  n_pcs_actual <- ncol(Embeddings(seurat_sub, "pca"))
  n_pcs_use <- min(N_PCS, n_pcs_actual)
  
  cat(sprintf("✓ PCA computed: %d PCs (will use %d)\n", n_pcs_actual, n_pcs_use))
  
  if (n_pcs_use < 10) {
    cat(sprintf("⚠ Warning: Only %d PCs available. Results may be unstable.\n", n_pcs_use))
  }
  
  # Clear scale.data after PCA
  tryCatch({
    if ("scale.data" %in% names(seurat_sub[["RNA"]]@layers)) {
      cat("Clearing scale.data after PCA...\n")
      seurat_sub[["RNA"]]@layers[["scale.data"]] <- NULL
      invisible(gc())
    }
  }, error = function(e) {
    # Silently continue if layers access fails
  })
  
  # Harmony
  cat("\n--- Harmony Integration ---\n")
  seurat_sub <- run_harmony_safe(seurat_sub, BATCH_COL, n_pcs_use)
  
  # UMAP + Clustering
  cat("\n--- UMAP + Clustering ---\n")
  set.seed(RANDOM_SEED)
  seurat_sub <- RunUMAP(
    seurat_sub,
    reduction = "harmony",
    dims = 1:n_pcs_use,  # P0-1 FIX: Use bounded n_pcs_use
    reduction.name = "umap_harmony",
    reduction.key = "umaph_",
    seed.use = RANDOM_SEED,
    verbose = FALSE
  )
  
  seurat_sub <- FindNeighbors(
    seurat_sub, 
    reduction = "harmony", 
    dims = 1:n_pcs_use,  # P0-1 FIX
    verbose = FALSE
  )
  
  # FindClusters with seed compatibility
  set.seed(RANDOM_SEED)
  fc_args <- list(
    object = seurat_sub,
    resolution = CLUSTER_RES,
    algorithm = CLUSTER_ALGO,
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
  cat(sprintf("✓ Found %d clusters\n", length(unique(seurat_sub$leiden_harmony))))
  
  # Save object
  sub_rds <- file.path(OBJ_DIR, sprintf("%s_harmony_clustered.rds", output_prefix))
  saveRDS(seurat_sub, sub_rds, compress = "xz")
  cat(sprintf("✓ Saved: %s\n", sub_rds))
  
  # UMAP plots
  cat("\n--- Generating UMAP plots ---\n")
  pdf(file.path(FIG_DIR, sprintf("%s_umap_overview.pdf", output_prefix)), width = 12, height = 10)
  
  print(DimPlot(seurat_sub, reduction = "umap_harmony", group.by = "leiden_harmony", 
                label = TRUE, repel = TRUE, raster = TRUE) +
          ggtitle(sprintf("%s - Leiden Clusters", celltype_name)))
  
  if (!is.null(celltype_filter)) {
    print(DimPlot(seurat_sub, reduction = "umap_harmony", group.by = MANUAL_COL, 
                  label = TRUE, repel = TRUE, raster = TRUE) +
            ggtitle(sprintf("%s - %s", celltype_name, MANUAL_COL)))
  }
  
  print(DimPlot(seurat_sub, reduction = "umap_harmony", group.by = BATCH_COL, 
                raster = TRUE) +
          ggtitle(sprintf("%s - %s", celltype_name, BATCH_COL)))
  
  dev.off()
  
  # FindAllMarkers
  cat("\n--- Finding Markers ---\n")
  Idents(seurat_sub) <- "leiden_harmony"
  
  markers <- find_markers_safe(seurat_sub)
  
  marker_csv <- file.path(TAB_DIR, sprintf("%s_markers.csv", output_prefix))
  write.csv(markers, marker_csv, row.names = FALSE)
  cat(sprintf("✓ Saved %d markers: %s\n", nrow(markers), marker_csv))
  
  # Heatmap
  cat("\n--- Generating Heatmap ---\n")
  if (nrow(markers) > 0) {
    top_markers <- markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = TOP_N_HEATMAP, with_ties = FALSE) %>%
      ungroup()
    
    hm_genes <- unique(top_markers$gene)
    hm_genes <- hm_genes[hm_genes %in% rownames(seurat_sub)]
    
    if (length(hm_genes) >= 2) {
      n_cells_hm <- min(ncol(seurat_sub), 10000)
      if (ncol(seurat_sub) > n_cells_hm) {
        cat(sprintf("Downsampling for heatmap (stratified)...\n"))
        set.seed(RANDOM_SEED)
        cells_hm <- stratified_sample_for_heatmap(seurat_sub, n_cells_hm, "leiden_harmony")
        seurat_sub_hm <- subset(seurat_sub, cells = cells_hm)
      } else {
        seurat_sub_hm <- seurat_sub
      }
      
      # Clear old scale.data before scaling
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
    } else {
      cat("⚠ Not enough valid genes for heatmap\n")
    }
  }
  
  # ROGUE
  if (RUN_ROGUE) {
    cat("\n--- Computing ROGUE ---\n")
    
    if (!is.null(celltype_filter)) {
      cat(sprintf("Skipping ROGUE by %s (single cell type)\n", MANUAL_COL))
    } else {
      rogue_manual <- compute_rogue_by_combo(
        seurat_sub,
        dataset_col = BATCH_COL,
        group_col   = MANUAL_COL,
        min_cells_per_group = MIN_CELLS_PER_GROUP,
        min_cells_gene = ROGUE_MIN_CELLS_GENE,
        min_genes_cell = ROGUE_MIN_GENES_CELL,
        max_cells_group = ROGUE_MAX_CELLS_GROUP,
        max_genes_group = ROGUE_MAX_GENES_GROUP
      )
      save_rogue_outputs(rogue_manual, sprintf("%s_manual", output_prefix), 
                         paste0(celltype_name, " | ", MANUAL_COL, " | "))
    }
    
    rogue_cluster <- compute_rogue_by_combo(
      seurat_sub,
      dataset_col = BATCH_COL,
      group_col   = "leiden_harmony",
      min_cells_per_group = MIN_CELLS_PER_GROUP,
      min_cells_gene = ROGUE_MIN_CELLS_GENE,
      min_genes_cell = ROGUE_MIN_GENES_CELL,
      max_cells_group = ROGUE_MAX_CELLS_GROUP,
      max_genes_group = ROGUE_MAX_GENES_GROUP
    )
    save_rogue_outputs(rogue_cluster, sprintf("%s_clusters", output_prefix),
                       paste0(celltype_name, " | leiden_harmony | "))
  }
  
  cat(sprintf("\n✓ Completed: %s\n", celltype_name))
  invisible(gc())
  
  return(sub_rds)
}

#' Compute ROGUE
compute_rogue_by_combo <- function(seurat_obj,
                                    dataset_col,
                                    group_col,
                                    min_cells_per_group = 30,
                                    min_cells_gene = 10,
                                    min_genes_cell = 200,
                                    max_cells_group = 2000,
                                    max_genes_group = 3000) {
  
  if (!requireNamespace("ROGUE", quietly = TRUE)) {
    stop("Package 'ROGUE' not installed.")
  }
  
  meta <- seurat_obj@meta.data
  stopifnot(dataset_col %in% colnames(meta), group_col %in% colnames(meta))
  
  combo_counts <- meta %>%
    dplyr::group_by(.data[[dataset_col]], .data[[group_col]]) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_cells >= min_cells_per_group)
  
  cat(sprintf("ROGUE combos (%s × %s) with n≥%d: %d\n",
              dataset_col, group_col, min_cells_per_group, nrow(combo_counts)))
  
  if (nrow(combo_counts) == 0) {
    cat("⚠ No valid combos. Skipping.\n")
    return(data.frame())
  }
  
  # P0-2 FIX: Layer/slot compatible counts retrieval
  expr_all <- tryCatch({
    GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts")
  }, error = function(e) {
    GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), slot = "counts")
  })
  
  if (is.null(rownames(expr_all))) rownames(expr_all) <- rownames(seurat_obj)
  if (is.null(colnames(expr_all))) colnames(expr_all) <- colnames(seurat_obj)
  
  rogue_list <- vector("list", nrow(combo_counts))
  pb <- txtProgressBar(min = 0, max = nrow(combo_counts), style = 3)
  
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
            vf <- head(vf, max_genes_group)
            mat <- mat[vf, , drop = FALSE]
          } else {
            rs <- Matrix::rowSums(mat)
            topg <- names(sort(rs, decreasing = TRUE))[seq_len(max_genes_group)]
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
      n_genes_used <<- NA_integer_
    })
    
    rogue_list[[i]] <- data.frame(
      dataset = as.character(ds),
      group   = as.character(grp),
      n_cells = as.integer(n0),
      n_cells_used = as.integer(n_cells_used),
      n_genes_used = as.integer(n_genes_used),
      rogue_value  = as.numeric(rogue_val),
      stringsAsFactors = FALSE
    )
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  dplyr::bind_rows(rogue_list)
}

#' Save ROGUE outputs
save_rogue_outputs <- function(rogue_df, out_prefix, title_prefix = "") {
  if (nrow(rogue_df) == 0) {
    cat("⚠ Empty ROGUE results. Skipping.\n")
    return(invisible(NULL))
  }
  
  write.csv(rogue_df, file.path(TAB_DIR, paste0(out_prefix, "_rogue_long.csv")), 
            row.names = FALSE)
  
  rogue_wide <- rogue_df %>%
    select(dataset, group, rogue_value) %>%
    tidyr::pivot_wider(names_from = group, values_from = rogue_value)
  write.csv(rogue_wide, file.path(TAB_DIR, paste0(out_prefix, "_rogue_wide.csv")), 
            row.names = FALSE)
  
  by_group <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    group_by(group) %>%
    summarise(
      n_datasets = n(),
      total_cells = sum(n_cells),
      median_rogue = median(rogue_value),
      mean_rogue = mean(rogue_value),
      sd_rogue = sd(rogue_value),
      .groups = "drop"
    ) %>% arrange(desc(median_rogue))
  write.csv(by_group, file.path(TAB_DIR, paste0(out_prefix, "_rogue_summary_by_group.csv")), 
            row.names = FALSE)
  
  by_dataset <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    group_by(dataset) %>%
    summarise(
      n_groups = n(),
      total_cells = sum(n_cells),
      median_rogue = median(rogue_value),
      mean_rogue = mean(rogue_value),
      .groups = "drop"
    ) %>% arrange(desc(median_rogue))
  write.csv(by_dataset, file.path(TAB_DIR, paste0(out_prefix, "_rogue_summary_by_dataset.csv")), 
            row.names = FALSE)
  
  # Plots
  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_boxplot.pdf")), width = 12, height = 6)
  p1 <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    ggplot(aes(x = reorder(group, rogue_value, median), y = rogue_value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.8) +
    geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red") +
    labs(title = paste0(title_prefix, "ROGUE by Group"),
         x = "Group", y = "ROGUE Score") +
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
    
    pheatmap::pheatmap(
      mat, cluster_rows = TRUE, cluster_cols = TRUE, na_col = "grey90",
      main = paste0(title_prefix, "ROGUE Heatmap: Dataset × Group"),
      fontsize = 8
    )
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
  
  saveRDS(
    list(long = rogue_df, wide = rogue_wide, summary_group = by_group, summary_dataset = by_dataset),
    file.path(OBJ_DIR, paste0(out_prefix, "_rogue_results.rds"))
  )
  
  cat(sprintf("✓ ROGUE outputs saved: %s\n", out_prefix))
}

# =========================
# Main Loop
# =========================
cat("\n========================================\n")
cat("Starting Main Analysis Loop\n")
cat("========================================\n\n")

results_paths <- list()

for (celltype_name in names(cell_types_to_analyze)) {
  celltype_filter <- cell_types_to_analyze[[celltype_name]]
  output_prefix <- gsub("[^A-Za-z0-9_-]", "_", celltype_name)
  
  result_path <- tryCatch({
    res <- process_subset(seurat_obj, celltype_name, celltype_filter, output_prefix)
    # #region agent log
    log_entry <- jsonlite::toJSON(list(
      sessionId = "debug-session",
      runId = "pre-fix",
      hypothesisId = "E",
      location = "main_loop:process_subset_success",
      message = "process_subset returned value",
      data = list(
        celltype_name = celltype_name,
        result_type = typeof(res),
        result_class = class(res)[1],
        result_is_null = is.null(res),
        result_is_char = is.character(res),
        result_length = if(!is.null(res)) length(res) else 0
      ),
      timestamp = as.numeric(Sys.time()) * 1000
    ), auto_unbox = TRUE)
    cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
    # #endregion
    res
  }, error = function(e) {
    cat(sprintf("\n❌ ERROR in %s: %s\n", celltype_name, e$message))
    cat("Traceback:\n")
    print(e)
    # #region agent log
    log_entry <- jsonlite::toJSON(list(
      sessionId = "debug-session",
      runId = "pre-fix",
      hypothesisId = "E",
      location = "main_loop:process_subset_error",
      message = "process_subset error",
      data = list(
        celltype_name = celltype_name,
        error_msg = e$message
      ),
      timestamp = as.numeric(Sys.time()) * 1000
    ), auto_unbox = TRUE)
    cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
    # #endregion
    NULL
  })
  
  # #region agent log
  log_entry <- jsonlite::toJSON(list(
    sessionId = "debug-session",
    runId = "pre-fix",
    hypothesisId = "C",
    location = "main_loop:before_assign",
    message = "Before assigning to results_paths",
    data = list(
      celltype_name = celltype_name,
      result_type = typeof(result_path),
      result_class = class(result_path)[1],
      result_is_null = is.null(result_path),
      result_is_char = is.character(result_path),
      result_is_list = is.list(result_path)
    ),
    timestamp = as.numeric(Sys.time()) * 1000
  ), auto_unbox = TRUE)
  cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
  # #endregion
  
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
# #region agent log
log_entry <- jsonlite::toJSON(list(
  sessionId = "debug-session",
  runId = "post-fix",
  hypothesisId = "C",
  location = "finalize:before_vapply",
  message = "Before vapply call",
  data = list(
    results_paths_length = length(results_paths),
    results_paths_names = names(results_paths),
    results_paths_types = vapply(results_paths, typeof, character(1)),
    results_paths_classes = vapply(results_paths, function(x) class(x)[1], character(1)),
    results_paths_is_null = vapply(results_paths, is.null, logical(1))
  ),
  timestamp = as.numeric(Sys.time()) * 1000
), auto_unbox = TRUE)
cat(log_entry, "\n", file = "/home/h2048/.cursor/debug.log", append = TRUE)
# #endregion

# Use vapply instead of sapply to ensure logical(1) return type
successful <- sum(vapply(results_paths, function(x) !is.null(x), logical(1)))
cat(sprintf("Successfully completed: %d / %d cell types\n", successful, length(results_paths)))
for (ct_name in names(results_paths)) {
  if (!is.null(results_paths[[ct_name]])) {
    cat(sprintf("  ✓ %s: %s\n", ct_name, results_paths[[ct_name]]))
  } else {
    cat(sprintf("  ✗ %s: FAILED\n", ct_name))
  }
}
cat("\n")

cat("=== Session Info ===\n")
print(sessionInfo())

# Sink will be closed by on.exit()
cat("\n✓ All done. Check", LOG_FILE, "for details.\n")
