#!/usr/bin/env Rscript
# =============================================================================
# Harmony-based Subset Analysis + Markers + ROGUE (Seurat v5 / Assay5)
# VERSION: v1.2 (2025-12-19)
#
# Must-fix implemented:
# - P0: Force fresh HVG->Scale->PCA from uncorrected data for Harmony input
# - P1: Manual_Annotation must exist (no Idents fallback)
# - P1: Harmony API compatibility (reduction vs reduction.use, dims vs dims.use)
# - P1: Input validation (file exists, columns exist, >=2 batches, subset size)
#
# Should-fix implemented:
# - Unified deterministic downsampling utilities
# - QC plots (UMAP + key markers)
# - Params recorded to seurat@misc
# - saveRDS with xz compression
# - ROGUE sparse-stage filter once (no redundant filtering by default)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =========================
# Configuration
# =========================
INPUT_RDS  <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/R/1218/harmony_subset_rogue"

MANUAL_COL <- "Manual_Annotation"   # MUST exist
BATCH_COL  <- "dataset"             # Harmony batch key

# subset selection: NULL keeps all
TARGET_MANUAL <- NULL
# TARGET_MANUAL <- c("Basal", "Goblet")

# Harmony input PCA (fresh)
PCA_NFEATURES <- 2000
N_PCS         <- 50

# UMAP / clustering
CLUSTER_RES   <- 0.8
CLUSTER_ALGO  <- 4   # 4=Leiden (fallback to 1 if not available)

# Markers / heatmap
MARKER_MIN_PCT <- 0.25
MARKER_LOGFC   <- 0.25
MARKER_TEST    <- "wilcox"
TOP_N_HEATMAP  <- 10

# Unified downsampling (deterministic)
SEED <- 42
MAX_CELLS_PER_CLUSTER_FOR_MARKERS <- 5000
MAX_CELLS_PER_CLUSTER_FOR_HEATMAP <- 500
ROGUE_MAX_CELLS_PER_COMBO         <- 2000

# ROGUE
RUN_ROGUE <- TRUE
MIN_CELLS_PER_GROUP  <- 30
ROGUE_MIN_CELLS_GENE <- 10
ROGUE_MIN_GENES_CELL <- 200
ROGUE_MAX_GENES      <- 3000
ROGUE_USE_MATR_FILTER <- FALSE  # default: avoid redundant filtering

# QC marker panel (edit as needed)
QC_MARKERS <- c("TP63", "KRT5", "KRT14", "KRT13", "MUC5AC", "SCGB1A1", "FOXJ1", "AGER", "SFTPC")

# saveRDS compression
RDS_COMPRESS <- "xz"  # "xz" | TRUE | FALSE

# =========================
# Setup output + logging
# =========================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_DIR <- file.path(OUTPUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
TAB_DIR <- file.path(OUTPUT_DIR, "tables");  dir.create(TAB_DIR, showWarnings = FALSE)
OBJ_DIR <- file.path(OUTPUT_DIR, "objects"); dir.create(OBJ_DIR, showWarnings = FALSE)
LOG_DIR <- file.path(OUTPUT_DIR, "logs");    dir.create(LOG_DIR, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, paste0("run_v1.2_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(LOG_FILE, split = TRUE)

section <- function(title) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat(title, "\n")
  cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("=", 80), "\n\n", sep = "")
}

cat("Input : ", INPUT_RDS, "\n", sep = "")
cat("Output: ", OUTPUT_DIR, "\n", sep = "")
cat("Log   : ", LOG_FILE, "\n\n", sep = "")

# =========================
# Package checks
# =========================
if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("Package 'harmony' not installed.")
}
if (RUN_ROGUE && !requireNamespace("ROGUE", quietly = TRUE)) {
  stop("Package 'ROGUE' not installed.")
}

cat("Seurat       :", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject :", as.character(packageVersion("SeuratObject")), "\n")
cat("harmony      :", as.character(packageVersion("harmony")), "\n")
if (RUN_ROGUE) cat("ROGUE        :", as.character(packageVersion("ROGUE")), "\n")
cat("\n")

# =============================================================================
# Helpers
# =============================================================================
stopf <- function(...) stop(sprintf(...), call. = FALSE)

validate_inputs <- function(path, manual_col, batch_col) {
  if (!file.exists(path)) stopf("❌ Input file not found: %s", path)
  sz_gb <- file.size(path) / 1e9
  cat(sprintf("📁 Loading %s (%.2f GB)\n", basename(path), sz_gb))
  obj <- readRDS(path)

  if (!inherits(obj, "Seurat")) stop("❌ Input is not a Seurat object.")
  if (!"RNA" %in% names(obj@assays)) stop("❌ RNA assay not found.")

  md <- obj@meta.data
  if (!manual_col %in% colnames(md)) {
    stopf("❌ '%s' not found in meta.data.\nAvailable columns (head): %s",
          manual_col, paste(head(colnames(md), 30), collapse = ", "))
  }
  if (!batch_col %in% colnames(md)) {
    stopf("❌ '%s' not found in meta.data.\nAvailable columns (head): %s",
          batch_col, paste(head(colnames(md), 30), collapse = ", "))
  }

  n_batches <- length(unique(md[[batch_col]]))
  if (n_batches < 2) stopf("❌ Only %d batch found in '%s'. Harmony requires ≥2.", n_batches, batch_col)
  cat(sprintf("✓ %d batches detected (%s)\n", n_batches, batch_col))

  # layers presence
  layers <- names(obj[["RNA"]]@layers)
  if (!"counts" %in% layers) stop("❌ RNA@layers[['counts']] not found.")
  if (!"data" %in% layers) {
    cat("⚠️ RNA@layers[['data']] not found -> will NormalizeData()\n")
  }

  obj
}

subset_by_manual <- function(obj, manual_col, target_manual) {
  if (is.null(target_manual)) {
    cat("TARGET_MANUAL = NULL -> keeping all cells\n")
    return(obj)
  }
  n_before <- ncol(obj)
  obj2 <- subset(obj, subset = .data[[manual_col]] %in% target_manual)
  n_after <- ncol(obj2)
  if (n_after < 100) stopf("❌ Only %d cells after subset (from %d). Too few.", n_after, n_before)
  cat(sprintf("✓ Subset: %s → %s cells (%.1f%%)\n",
              format(n_before, big.mark = ","),
              format(n_after, big.mark = ","),
              100 * n_after / n_before))
  obj2
}

clear_graphs_neighbors <- function(obj) {
  obj@graphs <- list()
  obj@neighbors <- list()
  obj
}

drop_scale_layer_if_exists <- function(obj, assay = "RNA") {
  lyr <- names(obj[[assay]]@layers)
  if ("scale.data" %in% lyr) {
    obj[[assay]]@layers[["scale.data"]] <- NULL
  }
  obj
}

# Unified downsampling: per group (Idents or meta column)
downsample_cells_by_group <- function(meta, group_col, max_per_group, seed = 42) {
  set.seed(seed)
  groups <- unique(meta[[group_col]])
  keep <- character(0)
  for (g in groups) {
    cells <- rownames(meta)[meta[[group_col]] == g]
    if (length(cells) > max_per_group) {
      keep <- c(keep, sample(cells, max_per_group))
    } else {
      keep <- c(keep, cells)
    }
  }
  keep
}

# Harmony compatibility wrapper
run_harmony_compat <- function(obj, group.by.vars, reduction_in, dims_vec, reduction.save = "harmony", verbose = FALSE) {
  # Prefer method formals to avoid partial matching surprises
  fn <- NULL
  fn <- tryCatch(getS3method("RunHarmony", "Seurat"), error = function(e) NULL)

  if (is.null(fn)) {
    # fallback: generic (still should dispatch)
    fn <- harmony::RunHarmony
  }

  fml <- names(formals(fn))
  cat("Harmony call formals:\n")
  cat(paste(fml, collapse = ", "), "\n\n")

  args <- list(object = obj, group.by.vars = group.by.vars, reduction.save = reduction.save, verbose = verbose)

  # reduction argument name
  if ("reduction.use" %in% fml) {
    args$reduction.use <- reduction_in
  } else if ("reduction" %in% fml) {
    args$reduction <- reduction_in
  } else {
    stop("❌ Cannot find reduction argument in harmony::RunHarmony.Seurat formals.")
  }

  # dims argument name
  if ("dims.use" %in% fml) {
    args$dims.use <- dims_vec
  } else if ("dims" %in% fml) {
    args$dims <- dims_vec
  } else {
    # some versions may not expose dims; safest to pass dims.use if exists; otherwise stop
    stop("❌ Cannot find dims argument (dims.use/dims) in harmony::RunHarmony.Seurat formals.")
  }

  # assay argument (optional)
  if ("assay.use" %in% fml) args$assay.use <- DefaultAssay(obj)

  do.call(fn, args)
}

safe_findclusters <- function(obj, res, algo_pref = 4) {
  out <- tryCatch({
    FindClusters(obj, resolution = res, algorithm = algo_pref, verbose = FALSE)
  }, error = function(e) {
    cat("⚠️ FindClusters(algo=4 Leiden) failed; fallback to Louvain (algo=1)\n")
    FindClusters(obj, resolution = res, algorithm = 1, verbose = FALSE)
  })
  out
}

# =============================================================================
# 1) Load + validate + subset
# =============================================================================
section("SECTION 1: Load + Validate + Subset")

seurat_obj <- validate_inputs(INPUT_RDS, MANUAL_COL, BATCH_COL)
DefaultAssay(seurat_obj) <- "RNA"

# Normalize if needed
if (!"data" %in% names(seurat_obj[["RNA"]]@layers)) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
}

seurat_sub <- subset_by_manual(seurat_obj, MANUAL_COL, TARGET_MANUAL)
rm(seurat_obj); invisible(gc())

cat(sprintf("Cells: %s | Genes: %s\n",
            format(ncol(seurat_sub), big.mark = ","),
            format(nrow(seurat_sub), big.mark = ",")))

cat(sprintf("Manual levels (n=%d): %s\n\n",
            length(unique(seurat_sub[[MANUAL_COL, drop = TRUE]])),
            paste(sort(unique(seurat_sub[[MANUAL_COL, drop = TRUE]])), collapse = ", ")))

# keep pipeline params in misc
seurat_sub@misc$pipeline_params <- list(
  version = "v1.2",
  input_rds = INPUT_RDS,
  output_dir = OUTPUT_DIR,
  manual_col = MANUAL_COL,
  batch_col = BATCH_COL,
  target_manual = TARGET_MANUAL,
  pca_nfeatures = PCA_NFEATURES,
  n_pcs = N_PCS,
  cluster_res = CLUSTER_RES,
  cluster_algo_pref = CLUSTER_ALGO,
  seed = SEED,
  downsample_markers = MAX_CELLS_PER_CLUSTER_FOR_MARKERS,
  downsample_heatmap = MAX_CELLS_PER_CLUSTER_FOR_HEATMAP,
  rogue_max_cells_combo = ROGUE_MAX_CELLS_PER_COMBO,
  timestamp = as.character(Sys.time())
)

# Clear graphs/neighbors to avoid confusion with legacy integrations
seurat_sub <- clear_graphs_neighbors(seurat_sub)

# =============================================================================
# 2) Force fresh HVG -> Scale -> PCA for Harmony input
# =============================================================================
section("SECTION 2: Fresh HVG -> Scale -> PCA (Harmony input)")

# Drop old scale layer to avoid feature-mismatch warning
seurat_sub <- drop_scale_layer_if_exists(seurat_sub, assay = "RNA")

# Try to reuse existing vf_vst_* if present and meaningful; otherwise recompute HVG
assay_md <- seurat_sub[["RNA"]]@meta.data
use_precomputed_vf <- FALSE
if ("vf_vst_counts_variable" %in% colnames(assay_md)) {
  n_vf <- sum(assay_md$vf_vst_counts_variable %in% TRUE, na.rm = TRUE)
  if (is.finite(n_vf) && n_vf >= PCA_NFEATURES) {
    use_precomputed_vf <- TRUE
    # prefer rank if present
    if ("vf_vst_counts_rank" %in% colnames(assay_md)) {
      ord <- order(assay_md$vf_vst_counts_rank, na.last = NA)
      vf <- rownames(assay_md)[ord]
      vf <- vf[seq_len(min(PCA_NFEATURES, length(vf)))]
    } else {
      vf <- rownames(assay_md)[which(assay_md$vf_vst_counts_variable %in% TRUE)]
      vf <- head(vf, PCA_NFEATURES)
    }
    VariableFeatures(seurat_sub) <- vf
    cat(sprintf("Using precomputed vf_vst_counts_variable HVGs: %d\n", length(vf)))
  }
}
if (!use_precomputed_vf) {
  cat("Computing HVGs with FindVariableFeatures(vst)...\n")
  seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = PCA_NFEATURES, verbose = FALSE)
  cat(sprintf("HVGs computed: %d\n", length(VariableFeatures(seurat_sub))))
}

cat("Scaling HVGs (temporary scale.data layer)...\n")
seurat_sub <- ScaleData(seurat_sub, features = VariableFeatures(seurat_sub), verbose = FALSE)

# Put PCA into a dedicated reduction name to avoid legacy 'pca' confusion
PCA_RED_NAME <- "pca_harmony_input"
if (PCA_RED_NAME %in% names(seurat_sub@reductions)) {
  seurat_sub@reductions[[PCA_RED_NAME]] <- NULL
}

cat("Running PCA (Harmony input)...\n")
seurat_sub <- RunPCA(
  object = seurat_sub,
  features = VariableFeatures(seurat_sub),
  npcs = N_PCS,
  reduction.name = PCA_RED_NAME,
  reduction.key = "pcain_",
  verbose = FALSE
)

# Drop scale.data immediately to avoid huge memory
seurat_sub <- drop_scale_layer_if_exists(seurat_sub, assay = "RNA")
invisible(gc())

# =============================================================================
# 3) Harmony -> UMAP -> Neighbors -> Clustering
# =============================================================================
section("SECTION 3: Harmony -> UMAP -> Clustering")

cat("Running Harmony (batch = ", BATCH_COL, ")...\n", sep = "")
seurat_sub <- run_harmony_compat(
  obj = seurat_sub,
  group.by.vars = BATCH_COL,
  reduction_in = PCA_RED_NAME,
  dims_vec = 1:N_PCS,
  reduction.save = "harmony",
  verbose = FALSE
)

# Name UMAP clearly
UMAP_RED_NAME <- "umap_harmony"
if (UMAP_RED_NAME %in% names(seurat_sub@reductions)) {
  # avoid overwriting silently
  UMAP_RED_NAME <- paste0("umap_harmony_res", gsub("\\.", "p", as.character(CLUSTER_RES)))
}
cat("Running UMAP (reduction = harmony) -> ", UMAP_RED_NAME, "\n", sep = "")
seurat_sub <- RunUMAP(
  object = seurat_sub,
  reduction = "harmony",
  dims = 1:N_PCS,
  reduction.name = UMAP_RED_NAME,
  reduction.key  = "umaph_",
  verbose = FALSE
)

cat("Neighbors + clustering...\n")
seurat_sub <- FindNeighbors(seurat_sub, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
seurat_sub <- safe_findclusters(seurat_sub, res = CLUSTER_RES, algo_pref = CLUSTER_ALGO)

# Save cluster id
CLUSTER_COL <- paste0("leiden_harmony_res", gsub("\\.", "p", as.character(CLUSTER_RES)))
seurat_sub[[CLUSTER_COL]] <- Idents(seurat_sub)

cat("Cluster column: ", CLUSTER_COL, "\n", sep = "")
cat("Clusters (n): ", length(levels(seurat_sub[[CLUSTER_COL]][, 1])), "\n\n", sep = "")

# Save object (post clustering)
SUB_RDS <- file.path(OBJ_DIR, "seurat_harmony_clustered.rds")
saveRDS(seurat_sub, SUB_RDS, compress = RDS_COMPRESS)
cat("Saved clustered object: ", SUB_RDS, "\n\n", sep = "")

# =============================================================================
# 4) QC Plots
# =============================================================================
section("SECTION 4: QC Plots")

pdf(file.path(FIG_DIR, "qc_umap_harmony.pdf"), width = 12, height = 9)
print(DimPlot(seurat_sub, reduction = UMAP_RED_NAME, group.by = CLUSTER_COL, label = TRUE, repel = TRUE, raster = TRUE) +
        ggtitle(paste0("UMAP (Harmony) - ", CLUSTER_COL)))
print(DimPlot(seurat_sub, reduction = UMAP_RED_NAME, group.by = MANUAL_COL, label = TRUE, repel = TRUE, raster = TRUE) +
        ggtitle(paste0("UMAP (Harmony) - ", MANUAL_COL)))
print(DimPlot(seurat_sub, reduction = UMAP_RED_NAME, group.by = BATCH_COL, raster = TRUE) +
        ggtitle(paste0("UMAP (Harmony) - ", BATCH_COL)))
dev.off()

# Marker preservation (FeaturePlot)
genes_present <- QC_MARKERS[QC_MARKERS %in% rownames(seurat_sub)]
if (length(genes_present) > 0) {
  pdf(file.path(FIG_DIR, "qc_featureplots_key_markers.pdf"), width = 12, height = 9)
  for (g in genes_present) {
    print(FeaturePlot(seurat_sub, features = g, reduction = UMAP_RED_NAME, raster = TRUE) +
            ggtitle(paste0("FeaturePlot: ", g)))
  }
  dev.off()
} else {
  cat("QC_MARKERS not found in rownames; skipping FeaturePlot QC.\n")
}

# =============================================================================
# 5) FindAllMarkers (with unified downsampling) + Heatmap
# =============================================================================
section("SECTION 5: FindAllMarkers + Heatmap")

Idents(seurat_sub) <- seurat_sub[[CLUSTER_COL, drop = TRUE]]

cells_for_markers <- downsample_cells_by_group(
  meta = seurat_sub@meta.data,
  group_col = CLUSTER_COL,
  max_per_group = MAX_CELLS_PER_CLUSTER_FOR_MARKERS,
  seed = SEED
)
cat(sprintf("Cells for markers (downsampled): %s\n", format(length(cells_for_markers), big.mark = ",")))

obj_mrk <- subset(seurat_sub, cells = cells_for_markers)

markers <- FindAllMarkers(
  object = obj_mrk,
  only.pos = TRUE,
  min.pct = MARKER_MIN_PCT,
  logfc.threshold = MARKER_LOGFC,
  test.use = MARKER_TEST,
  assay = DefaultAssay(obj_mrk),
  slot = "data"
)

MARKER_CSV <- file.path(TAB_DIR, "markers_FindAllMarkers_harmony_clusters.csv")
write.csv(markers, MARKER_CSV, row.names = FALSE)
cat("Saved markers: ", MARKER_CSV, "\n", sep = "")

# Top markers per cluster -> heatmap genes
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N_HEATMAP, with_ties = FALSE) %>%
  ungroup()

hm_genes <- unique(top_markers$gene)
hm_genes <- hm_genes[hm_genes %in% rownames(seurat_sub)]
cat(sprintf("Heatmap genes: %d\n", length(hm_genes)))

# Heatmap downsample
cells_for_heatmap <- downsample_cells_by_group(
  meta = seurat_sub@meta.data,
  group_col = CLUSTER_COL,
  max_per_group = MAX_CELLS_PER_CLUSTER_FOR_HEATMAP,
  seed = SEED
)
cat(sprintf("Cells for heatmap (downsampled): %s\n", format(length(cells_for_heatmap), big.mark = ",")))

# Scale only heatmap genes (temporary)
seurat_sub <- drop_scale_layer_if_exists(seurat_sub, assay = "RNA")
seurat_sub <- ScaleData(seurat_sub, features = hm_genes, verbose = FALSE)

pdf(file.path(FIG_DIR, "heatmap_top_markers_harmony_clusters.pdf"), width = 14, height = 10)
print(DoHeatmap(
  object = subset(seurat_sub, cells = cells_for_heatmap),
  features = hm_genes,
  group.by = CLUSTER_COL,
  raster = TRUE
) + ggtitle(paste0("Top markers per ", CLUSTER_COL, " (scaled on selected genes)")))
dev.off()

# Drop scale layer again
seurat_sub <- drop_scale_layer_if_exists(seurat_sub, assay = "RNA")
invisible(gc())

# Save post-marker object
SUB_RDS2 <- file.path(OBJ_DIR, "seurat_harmony_clustered_postmarkers.rds")
saveRDS(seurat_sub, SUB_RDS2, compress = RDS_COMPRESS)
cat("Saved post-marker object: ", SUB_RDS2, "\n", sep = "")

# =============================================================================
# 6) ROGUE (dataset × Manual_Annotation) and (dataset × Harmony clusters)
# =============================================================================
compute_rogue_by_combo <- function(seurat_obj,
                                  dataset_col,
                                  group_col,
                                  min_cells_per_group,
                                  min_cells_gene,
                                  min_genes_cell,
                                  max_cells_combo,
                                  max_genes,
                                  use_matr_filter = FALSE,
                                  counts_layer = "counts") {
  stopifnot(requireNamespace("ROGUE", quietly = TRUE))

  meta <- seurat_obj@meta.data
  stopifnot(dataset_col %in% colnames(meta), group_col %in% colnames(meta))

  combo_counts <- meta %>%
    group_by(.data[[dataset_col]], .data[[group_col]]) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    filter(n_cells >= min_cells_per_group)

  cat(sprintf("Valid combos (%s × %s), n>=%d: %d\n",
              dataset_col, group_col, min_cells_per_group, nrow(combo_counts)))

  expr_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = counts_layer)

  rogue_list <- vector("list", nrow(combo_counts))
  pb <- txtProgressBar(min = 0, max = nrow(combo_counts), style = 3)

  for (i in seq_len(nrow(combo_counts))) {
    ds  <- combo_counts[[dataset_col]][i]
    grp <- combo_counts[[group_col]][i]
    n0  <- combo_counts$n_cells[i]

    cells_keep <- rownames(meta)[meta[[dataset_col]] == ds & meta[[group_col]] == grp]

    # deterministic downsample per combo
    if (length(cells_keep) > max_cells_combo) {
      set.seed(SEED)
      cells_keep <- sample(cells_keep, max_cells_combo)
    }

    rogue_val <- NA_real_
    n_cells_used <- length(cells_keep)
    n_genes_used <- NA_integer_

    tryCatch({
      mat <- expr_all[, cells_keep, drop = FALSE]

      # Sparse-stage single-pass filtering (avoid redundancy)
      gene_detect <- Matrix::rowSums(mat > 0)
      mat <- mat[gene_detect >= min_cells_gene, , drop = FALSE]

      cell_detect <- Matrix::colSums(mat > 0)
      mat <- mat[, cell_detect >= min_genes_cell, drop = FALSE]

      if (ncol(mat) < min_cells_per_group || nrow(mat) < 100) {
        rogue_val <- NA_real_
      } else {
        # cap genes by expression
        if (nrow(mat) > max_genes) {
          rs <- Matrix::rowSums(mat)
          topg <- names(sort(rs, decreasing = TRUE))[seq_len(max_genes)]
          mat <- mat[topg, , drop = FALSE]
        }

        n_genes_used <- nrow(mat)
        n_cells_used <- ncol(mat)

        dense <- as.matrix(mat)

        if (use_matr_filter) {
          dense <- ROGUE::matr.filter(
            dense,
            min.cells = max(3, floor(min_cells_gene / 2)),
            min.genes = max(50, floor(min_genes_cell / 2))
          )
        }

        if (ncol(dense) < min_cells_per_group || nrow(dense) < 100) {
          rogue_val <- NA_real_
        } else {
          ent_res <- ROGUE::SE_fun(dense)
          if (any(!is.finite(ent_res$entropy))) {
            rogue_val <- NA_real_
          } else {
            rogue_val <- ROGUE::CalculateRogue(ent_res, platform = "UMI")
          }
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

  bind_rows(rogue_list)
}

save_rogue_outputs <- function(rogue_df, out_prefix, title_prefix = "") {
  write.csv(rogue_df, file.path(TAB_DIR, paste0(out_prefix, "_rogue_long.csv")), row.names = FALSE)

  rogue_wide <- rogue_df %>%
    select(dataset, group, rogue_value) %>%
    pivot_wider(names_from = group, values_from = rogue_value)
  write.csv(rogue_wide, file.path(TAB_DIR, paste0(out_prefix, "_rogue_wide.csv")), row.names = FALSE)

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
  write.csv(by_group, file.path(TAB_DIR, paste0(out_prefix, "_rogue_summary_by_group.csv")), row.names = FALSE)

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
  write.csv(by_dataset, file.path(TAB_DIR, paste0(out_prefix, "_rogue_summary_by_dataset.csv")), row.names = FALSE)

  # plots
  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_boxplot_by_group.pdf")), width = 12, height = 6)
  p1 <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    ggplot(aes(x = reorder(group, rogue_value, median), y = rogue_value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.6) +
    geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed") +
    labs(title = paste0(title_prefix, "ROGUE by Group"), x = "Group", y = "ROGUE") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1); dev.off()

  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_dotplot.pdf")), width = 14, height = 8)
  p2 <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    ggplot(aes(x = group, y = dataset, size = n_cells, color = rogue_value)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 10)) +
    labs(title = paste0(title_prefix, "ROGUE: Group × Dataset"),
         x = "Group", y = "Dataset", size = "Cells", color = "ROGUE") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p2); dev.off()

  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_heatmap.pdf")), width = 14, height = 10)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    mat <- rogue_df %>%
      select(dataset, group, rogue_value) %>%
      pivot_wider(names_from = group, values_from = rogue_value) %>%
      tibble::column_to_rownames("dataset") %>%
      as.matrix()
    pheatmap::pheatmap(
      mat,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      na_col = "grey90",
      main = paste0(title_prefix, "ROGUE Heatmap: Dataset × Group"),
      fontsize = 8
    )
  } else {
    plot.new(); text(0.5, 0.5, "pheatmap not installed")
  }
  dev.off()

  saveRDS(
    list(long = rogue_df, wide = rogue_wide, summary_group = by_group, summary_dataset = by_dataset),
    file.path(OBJ_DIR, paste0(out_prefix, "_rogue_results.rds")),
    compress = RDS_COMPRESS
  )
}

if (RUN_ROGUE) {
  section("SECTION 6: ROGUE")

  suppressPackageStartupMessages(library(ROGUE))

  cat("=== ROGUE: dataset × Manual_Annotation ===\n")
  rogue_manual <- compute_rogue_by_combo(
    seurat_obj = seurat_sub,
    dataset_col = BATCH_COL,
    group_col   = MANUAL_COL,
    min_cells_per_group = MIN_CELLS_PER_GROUP,
    min_cells_gene = ROGUE_MIN_CELLS_GENE,
    min_genes_cell = ROGUE_MIN_GENES_CELL,
    max_cells_combo = ROGUE_MAX_CELLS_PER_COMBO,
    max_genes = ROGUE_MAX_GENES,
    use_matr_filter = ROGUE_USE_MATR_FILTER,
    counts_layer = "counts"
  )
  save_rogue_outputs(rogue_manual, out_prefix = "manual_annotation", title_prefix = paste0(MANUAL_COL, " | "))

  cat("\n=== ROGUE: dataset × Harmony clusters ===\n")
  rogue_cluster <- compute_rogue_by_combo(
    seurat_obj = seurat_sub,
    dataset_col = BATCH_COL,
    group_col   = CLUSTER_COL,
    min_cells_per_group = MIN_CELLS_PER_GROUP,
    min_cells_gene = ROGUE_MIN_CELLS_GENE,
    min_genes_cell = ROGUE_MIN_GENES_CELL,
    max_cells_combo = ROGUE_MAX_CELLS_PER_COMBO,
    max_genes = ROGUE_MAX_GENES,
    use_matr_filter = ROGUE_USE_MATR_FILTER,
    counts_layer = "counts"
  )
  save_rogue_outputs(rogue_cluster, out_prefix = "harmony_clusters", title_prefix = paste0(CLUSTER_COL, " | "))

  cat("\nROGUE done.\n")
  cat("Tables : ", TAB_DIR, "\n", sep = "")
  cat("Figures: ", FIG_DIR, "\n", sep = "")
  cat("Objects: ", OBJ_DIR, "\n", sep = "")
}

# =============================================================================
# Finish
# =============================================================================
section("SECTION 7: SessionInfo + Done")
print(sessionInfo())
cat("\nAll done.\n")

sink()

