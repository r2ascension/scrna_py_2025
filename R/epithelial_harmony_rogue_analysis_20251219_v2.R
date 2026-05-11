#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# =========================
# Configuration
# =========================
INPUT_RDS  <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/R/1218/harmony_subset_rogue"

MANUAL_COL <- "Manual_Annotation"   # If not exists, will create from Idents()
BATCH_COL  <- "sample"             # Harmony batch key

# Manual_Annotation categories to keep; set to NULL to keep all
TARGET_MANUAL <- NULL
# Example:
# TARGET_MANUAL <- c("Basal", "Goblet")

# Dimensionality reduction / clustering
N_VAR_FEATURES <- 4000
N_PCS          <- 50
CLUSTER_RES    <- 2
CLUSTER_ALGO   <- 4    # 4 = Leiden

# FindAllMarkers
MARKER_MIN_PCT <- 0.25
MARKER_LOGFC   <- 0.25
TOP_N_HEATMAP  <- 10

# ROGUE
RUN_ROGUE <- TRUE
MIN_CELLS_PER_GROUP   <- 30
ROGUE_MIN_CELLS_GENE  <- 10
ROGUE_MIN_GENES_CELL  <- 200
ROGUE_MAX_CELLS_GROUP <- 2000   # Downsample if too many cells per group (speed + memory)
ROGUE_MAX_GENES_GROUP <- 3000   # Gene cap per group (sparse filter then truncate)

# =========================
# IO
# =========================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_DIR   <- file.path(OUTPUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
TAB_DIR   <- file.path(OUTPUT_DIR, "tables");  dir.create(TAB_DIR, showWarnings = FALSE)
OBJ_DIR   <- file.path(OUTPUT_DIR, "objects"); dir.create(OBJ_DIR, showWarnings = FALSE)

cat("Input :", INPUT_RDS, "\n")
cat("Output:", OUTPUT_DIR, "\n\n")

# =========================
# Load Seurat
# =========================
seurat_obj <- readRDS(INPUT_RDS)
DefaultAssay(seurat_obj) <- "RNA"

# Ensure MANUAL_COL exists
if (!MANUAL_COL %in% colnames(seurat_obj@meta.data)) {
  cat(sprintf("'%s' not found in meta.data -> creating from Idents(seurat_obj)\n", MANUAL_COL))
  seurat_obj[[MANUAL_COL]] <- Idents(seurat_obj)
}
stopifnot(BATCH_COL %in% colnames(seurat_obj@meta.data))

# =========================
# Subset by Manual_Annotation
# =========================
if (!is.null(TARGET_MANUAL)) {
  cat("Subsetting by Manual_Annotation in TARGET_MANUAL...\n")
  seurat_sub <- subset(seurat_obj, subset = .data[[MANUAL_COL]] %in% TARGET_MANUAL)
} else {
  cat("TARGET_MANUAL = NULL -> keeping all cells\n")
  seurat_sub <- seurat_obj
}
rm(seurat_obj); invisible(gc())

cat(sprintf("Cells after subset: %s\n", format(ncol(seurat_sub), big.mark = ",")))
cat(sprintf("Manual levels (n=%d): %s\n\n",
            length(unique(seurat_sub[[MANUAL_COL, drop = TRUE]])),
            paste(head(sort(unique(seurat_sub[[MANUAL_COL, drop = TRUE]])), 20), collapse = ", ")))

# =========================
# Recompute HVG -> PCA -> Harmony -> UMAP -> Clustering
# =========================
# Normalize (write to data layer). If you already trust data layer, this is safe (overwrites deterministically).
seurat_sub <- NormalizeData(seurat_sub, verbose = FALSE)

seurat_sub <- FindVariableFeatures(
  seurat_sub, selection.method = "vst", nfeatures = N_VAR_FEATURES, verbose = FALSE
)

seurat_sub <- ScaleData(
  seurat_sub, features = VariableFeatures(seurat_sub), verbose = FALSE
)

seurat_sub <- RunPCA(
  seurat_sub, features = VariableFeatures(seurat_sub), npcs = N_PCS, verbose = FALSE
)

# ===== FIX: Use correct Harmony parameters =====
# Based on successful examples from project:
# - reduction.use (not reduction)
# - dims (not dims.use)
# - no assay.use parameter needed
seurat_sub <- RunHarmony(
  object = seurat_sub,
  group.by.vars = BATCH_COL,
  reduction.use = "pca",       # ✓ Correct parameter name
  dims = 1:N_PCS,              # ✓ Correct parameter name
  theta = NULL,                 # Use default
  lambda = NULL,                # Use default
  sigma = 0.1,
  nclust = NULL,                # Auto-determine
  max_iter = 20,
  early_stop = TRUE,
  verbose = FALSE
)

seurat_sub <- RunUMAP(
  seurat_sub,
  reduction = "harmony",
  dims = 1:N_PCS,
  reduction.name = "umap_harmony",
  reduction.key = "umaph_",
  verbose = FALSE
)

seurat_sub <- FindNeighbors(seurat_sub, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
seurat_sub <- FindClusters(seurat_sub, resolution = CLUSTER_RES, algorithm = CLUSTER_ALGO, verbose = FALSE)

# Save cluster id into a dedicated column
seurat_sub$leiden_harmony <- Idents(seurat_sub)

# Save object
SUB_RDS <- file.path(OBJ_DIR, "seurat_subset_harmony_clustered.rds")
saveRDS(seurat_sub, SUB_RDS)
cat("Saved clustered object:", SUB_RDS, "\n\n")

# =========================
# Quick plots (UMAP)
# =========================
pdf(file.path(FIG_DIR, "umap_harmony_overview.pdf"), width = 10, height = 8)
print(DimPlot(seurat_sub, reduction = "umap_harmony", group.by = "leiden_harmony", label = TRUE, repel = TRUE) +
        ggtitle("UMAP (Harmony) - leiden_harmony"))
print(DimPlot(seurat_sub, reduction = "umap_harmony", group.by = MANUAL_COL, label = TRUE, repel = TRUE) +
        ggtitle(paste0("UMAP (Harmony) - ", MANUAL_COL)))
print(DimPlot(seurat_sub, reduction = "umap_harmony", group.by = BATCH_COL, raster = TRUE) +
        ggtitle(paste0("UMAP (Harmony) - ", BATCH_COL)))
dev.off()

# =========================
# FindAllMarkers (by harmony clusters)
# =========================
Idents(seurat_sub) <- "leiden_harmony"

markers <- FindAllMarkers(
  seurat_sub,
  only.pos = TRUE,
  min.pct = MARKER_MIN_PCT,
  logfc.threshold = MARKER_LOGFC,
  test.use = "wilcox",
  assay = DefaultAssay(seurat_sub),
  slot = "data"   # Use log-normalized data layer
)

MARKER_CSV <- file.path(TAB_DIR, "markers_FindAllMarkers_leiden_harmony.csv")
write.csv(markers, MARKER_CSV, row.names = FALSE)
cat("Saved markers:", MARKER_CSV, "\n")

# Heatmap of top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N_HEATMAP, with_ties = FALSE) %>%
  ungroup()

hm_genes <- unique(top_markers$gene)
# Only scale heatmap genes to avoid creating huge scale.data matrix
seurat_sub <- ScaleData(seurat_sub, features = hm_genes, verbose = FALSE)

pdf(file.path(FIG_DIR, "heatmap_top_markers_leiden_harmony.pdf"), width = 12, height = 10)
print(DoHeatmap(seurat_sub, features = hm_genes, group.by = "leiden_harmony", raster = TRUE) +
        ggtitle("Top markers per leiden_harmony (scaled on selected genes only)"))
dev.off()

# =========================
# ROGUE helpers
# =========================
compute_rogue_by_combo <- function(seurat_obj,
                                  dataset_col,
                                  group_col,
                                  min_cells_per_group = 30,
                                  min_cells_gene = 10,
                                  min_genes_cell = 200,
                                  max_cells_group = 2000,
                                  max_genes_group = 3000,
                                  counts_layer = "counts") {

  if (!requireNamespace("ROGUE", quietly = TRUE)) {
    stop("Package 'ROGUE' not installed.")
  }

  meta <- seurat_obj@meta.data
  stopifnot(dataset_col %in% colnames(meta), group_col %in% colnames(meta))

  combo_counts <- meta %>%
    dplyr::group_by(.data[[dataset_col]], .data[[group_col]]) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_cells >= min_cells_per_group)

  cat(sprintf("ROGUE combos (%s x %s) with n>=%d: %d\n",
              dataset_col, group_col, min_cells_per_group, nrow(combo_counts)))

  # counts matrix (sparse)
  expr_all <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = counts_layer)
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
      # optional downsample cells to control runtime
      if (length(cells_keep) > max_cells_group) {
        set.seed(42)
        cells_keep <- sample(cells_keep, max_cells_group)
      }

      mat <- expr_all[, cells_keep, drop = FALSE]

      # filter genes by detected cells (sparse-safe)
      gene_detect <- Matrix::rowSums(mat > 0)
      keep_genes <- which(gene_detect >= min_cells_gene)
      mat <- mat[keep_genes, , drop = FALSE]

      # filter cells by detected genes (sparse-safe)
      cell_detect <- Matrix::colSums(mat > 0)
      keep_cells <- which(cell_detect >= min_genes_cell)
      mat <- mat[, keep_cells, drop = FALSE]

      if (ncol(mat) < min_cells_per_group || nrow(mat) < 100) {
        rogue_val <- NA_real_
      } else {
        # cap genes (prefer VariableFeatures if possible)
        if (nrow(mat) > max_genes_group) {
          vf <- VariableFeatures(seurat_obj)
          vf <- vf[vf %in% rownames(mat)]
          if (length(vf) >= 200) {
            vf <- head(vf, max_genes_group)
            mat <- mat[vf, , drop = FALSE]
          } else {
            # fallback: top expressed genes
            rs <- Matrix::rowSums(mat)
            topg <- names(sort(rs, decreasing = TRUE))[seq_len(max_genes_group)]
            mat <- mat[topg, , drop = FALSE]
          }
        }

        n_genes_used <- nrow(mat)
        n_cells_used <- ncol(mat)

        # ROGUE expects dense matrix in most common usage
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

save_rogue_outputs <- function(rogue_df, out_prefix, title_prefix = "") {
  # long / wide
  write.csv(rogue_df, file.path(TAB_DIR, paste0(out_prefix, "_rogue_long.csv")), row.names = FALSE)

  rogue_wide <- rogue_df %>%
    select(dataset, group, rogue_value) %>%
    tidyr::pivot_wider(names_from = group, values_from = rogue_value)
  write.csv(rogue_wide, file.path(TAB_DIR, paste0(out_prefix, "_rogue_wide.csv")), row.names = FALSE)

  # summaries
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
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.8) +
    geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed") +
    labs(
      title = paste0(title_prefix, "ROGUE by Group (Across Datasets)"),
      x = "Group",
      y = "ROGUE"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1); dev.off()

  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_heatmap.pdf")), width = 14, height = 10)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    mat <- rogue_df %>%
      select(dataset, group, rogue_value) %>%
      pivot_wider(names_from = group, values_from = rogue_value) %>%
      column_to_rownames("dataset") %>%
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

  pdf(file.path(FIG_DIR, paste0(out_prefix, "_rogue_dotplot.pdf")), width = 14, height = 8)
  p3 <- rogue_df %>%
    filter(!is.na(rogue_value)) %>%
    ggplot(aes(x = group, y = dataset, size = n_cells, color = rogue_value)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 10)) +
    labs(
      title = paste0(title_prefix, "ROGUE: Group × Dataset"),
      x = "Group", y = "Dataset", size = "Cells", color = "ROGUE"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p3); dev.off()

  saveRDS(
    list(long = rogue_df, wide = rogue_wide, summary_group = by_group, summary_dataset = by_dataset),
    file.path(OBJ_DIR, paste0(out_prefix, "_rogue_results.rds"))
  )
}

# =========================
# Run ROGUE
# =========================
if (RUN_ROGUE) {
  suppressPackageStartupMessages(library(ROGUE))

  # 1) dataset × Manual_Annotation
  rogue_manual <- compute_rogue_by_combo(
    seurat_sub,
    dataset_col = BATCH_COL,
    group_col   = MANUAL_COL,
    min_cells_per_group = MIN_CELLS_PER_GROUP,
    min_cells_gene = ROGUE_MIN_CELLS_GENE,
    min_genes_cell = ROGUE_MIN_GENES_CELL,
    max_cells_group = ROGUE_MAX_CELLS_GROUP,
    max_genes_group = ROGUE_MAX_GENES_GROUP,
    counts_layer = "counts"
  )
  save_rogue_outputs(rogue_manual, out_prefix = "manual_annotation", title_prefix = paste0(MANUAL_COL, " | "))

  # 2) dataset × harmony clusters
  rogue_cluster <- compute_rogue_by_combo(
    seurat_sub,
    dataset_col = BATCH_COL,
    group_col   = "leiden_harmony",
    min_cells_per_group = MIN_CELLS_PER_GROUP,
    min_cells_gene = ROGUE_MIN_CELLS_GENE,
    min_genes_cell = ROGUE_MIN_GENES_CELL,
    max_cells_group = ROGUE_MAX_CELLS_GROUP,
    max_genes_group = ROGUE_MAX_GENES_GROUP,
    counts_layer = "counts"
  )
  save_rogue_outputs(rogue_cluster, out_prefix = "leiden_harmony", title_prefix = "leiden_harmony | ")

  cat("\nROGUE done. Outputs in:\n", TAB_DIR, "\n", FIG_DIR, "\n", OBJ_DIR, "\n")
}

cat("\nAll done.\n")
