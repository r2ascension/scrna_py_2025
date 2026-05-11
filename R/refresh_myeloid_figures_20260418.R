#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

source("/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414_v3.R")

OBJECT_PATH <- "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.rds"
OUTPUT_DIR <- dirname(OBJECT_PATH)
FIG_DIR <- file.path(OUTPUT_DIR, "figures")
MANIFEST_PATH <- file.path(FIG_DIR, "figure_refresh_manifest.tsv")

TISSUE_COL <- "tissue"
SAMPLE_COL <- "sample"
CELLTYPE_L2_COL <- "cell_type_L2"
CELLTYPE_L3_COL <- "cell_type_L3"
LABEL_COL <- CELLTYPE_L3_COL
UMAP_REDUCTION_PREFERRED <- c("umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap")
CHOIR_REDUCTION_CANDIDATES <- c("scanvi", "scvi", "pca")
UMAP_PT_SIZE <- 0.35
HEATMAP_CELLS_PER_TYPE <- 100L
MARKER_PANEL_GROUP_COL <- CELLTYPE_L3_COL
MARKER_PANEL_FIG_SUBDIR <- "marker_panels"
MARKER_PANEL_HEATMAP_CELLS_PER_GROUP <- 80L
PLOT_SHUFFLE <- FALSE
REGENERATE_HEAVY_PLOTS <- FALSE

UMAP_TISSUE_COLORS <- c(
  "lung parenchyma"    = "#D55E00",
  "nose"               = "#009E73",
  "respiratory airway" = "#0072B2",
  "sinus"              = "#CC79A7"
)

pick_reduction <- function(obj, preferred) {
  red <- Reductions(obj)
  hit <- preferred[preferred %in% red]
  if (length(hit) > 0) hit[[1]] else NULL
}

build_umap_plot <- function(obj, reduction_name, group_col, title,
                            split_col = NULL, label = FALSE,
                            width = 10, height = 8,
                            cols = NULL,
                            shuffle = PLOT_SHUFFLE) {
  p <- DimPlot(
    obj,
    reduction = reduction_name,
    group.by = group_col,
    split.by = split_col,
    pt.size = UMAP_PT_SIZE,
    shuffle = shuffle,
    label = label,
    repel = label,
    cols = cols
  ) +
    ggtitle(title) +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  list(plot = p, width = width, height = height)
}

save_plot <- function(p, path_no_ext, width = 10, height = 8) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
}

stratified_downsample <- function(obj, group_col, n_per = HEATMAP_CELLS_PER_TYPE) {
  set.seed(42)
  meta <- obj@meta.data
  groups <- sort(unique(as.character(meta[[group_col]])))
  groups <- groups[!is.na(groups) & nzchar(trimws(groups))]
  cells <- unlist(lapply(groups, function(group_name) {
    group_cells <- rownames(meta)[as.character(meta[[group_col]]) == group_name]
    sample(group_cells, min(length(group_cells), n_per))
  }), use.names = FALSE)
  subset(obj, cells = cells)
}

record_manifest_row <- function(manifest_rows,
                                figure,
                                rel_path,
                                display_reduction = "",
                                cluster_source_reduction = "",
                                shuffle = PLOT_SHUFFLE,
                                note = "") {
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    figure = figure,
    path = rel_path,
    display_reduction = ifelse(is.null(display_reduction), "", display_reduction),
    cluster_source_reduction = ifelse(is.null(cluster_source_reduction), "", cluster_source_reduction),
    shuffle = as.character(shuffle),
    refreshed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    note = note,
    stringsAsFactors = FALSE
  )
  manifest_rows
}

cat("=== Refresh Myeloid Figures (2026-04-18) ===\n")
if (!file.exists(OBJECT_PATH)) {
  stop(sprintf("Final object not found: %s", OBJECT_PATH))
}
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cfg <- tc_myeloid_wrapper_preset()
known_markers <- cfg$KNOWN_MARKERS
marker_panels <- cfg$MARKER_PANELS

obj <- readRDS(OBJECT_PATH)
if (!"RNA" %in% names(obj@assays)) {
  stop("RNA assay not found in myeloid final object")
}
if (!"data" %in% Layers(obj[["RNA"]])) {
  cat("[INFO] data layer missing; running NormalizeData() before visualization\n")
  obj <- NormalizeData(obj, verbose = FALSE)
}

obj <- obj[, sort(colnames(obj))]
meta <- obj@meta.data

required_cols <- c(TISSUE_COL, SAMPLE_COL, CELLTYPE_L2_COL, CELLTYPE_L3_COL)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing metadata columns: %s", paste(missing_cols, collapse = ", ")))
}

umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
if (is.null(umap_reduction)) {
  stop(sprintf(
    "No display UMAP reduction found. Candidates: %s | available: %s",
    paste(UMAP_REDUCTION_PREFERRED, collapse = ", "),
    paste(Reductions(obj), collapse = ", ")
  ))
}
choir_reduction <- pick_reduction(obj, CHOIR_REDUCTION_CANDIDATES)
if (is.null(choir_reduction)) {
  stop(sprintf(
    "No CHOIR latent reduction found. Candidates: %s | available: %s",
    paste(CHOIR_REDUCTION_CANDIDATES, collapse = ", "),
    paste(Reductions(obj), collapse = ", ")
  ))
}
choir_col_hits <- grep("^CHOIR_clusters", colnames(meta), value = TRUE)
if (length(choir_col_hits) == 0L) {
  stop("No CHOIR cluster column found in metadata")
}
choir_col <- choir_col_hits[[1]]

cat(sprintf("[OK] display reduction: %s\n", umap_reduction))
cat(sprintf("[OK] CHOIR latent reduction: %s\n", choir_reduction))
cat(sprintf("[OK] CHOIR metadata column: %s\n", choir_col))
cat(sprintf("[OK] deterministic rendering: shuffle=%s | sorted cells by barcode\n", as.character(PLOT_SHUFFLE)))

manifest_rows <- list()
tissues <- sort(unique(as.character(meta[[TISSUE_COL]])))
tissues <- tissues[!is.na(tissues) & nzchar(trimws(tissues))]
tissue_cols_use <- UMAP_TISSUE_COLORS[names(UMAP_TISSUE_COLORS) %in% tissues]

# UMAP family ---------------------------------------------------------------
cat("[STEP] Redrawing UMAP family...\n")
p_tissue <- build_umap_plot(
  obj,
  umap_reduction,
  TISSUE_COL,
  title = sprintf("Myeloid - Tissue (%s)", umap_reduction),
  cols = tissue_cols_use,
  width = 11,
  height = 8
)
save_plot(p_tissue$plot, file.path(FIG_DIR, "umap_tissue"), p_tissue$width, p_tissue$height)
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "umap_tissue",
  rel_path = "figures/umap_tissue.png",
  display_reduction = umap_reduction,
  note = "Deterministic redraw on shared display UMAP"
)

p_l2 <- build_umap_plot(
  obj,
  umap_reduction,
  CELLTYPE_L2_COL,
  title = sprintf("Myeloid - Cell Type L2 (%s)", umap_reduction),
  label = TRUE,
  width = 14,
  height = 10
)
save_plot(p_l2$plot, file.path(FIG_DIR, "umap_celltype_L2"), p_l2$width, p_l2$height)
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "umap_celltype_L2",
  rel_path = "figures/umap_celltype_L2.png",
  display_reduction = umap_reduction,
  note = "Deterministic redraw on shared display UMAP"
)

p_l2_split <- build_umap_plot(
  obj,
  umap_reduction,
  CELLTYPE_L2_COL,
  title = sprintf("Myeloid - L2 by Tissue (%s)", umap_reduction),
  split_col = TISSUE_COL,
  label = TRUE,
  width = max(12, 4.5 * length(tissues)),
  height = 8
)
save_plot(p_l2_split$plot, file.path(FIG_DIR, "umap_L2_split_tissue"), p_l2_split$width, p_l2_split$height)
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "umap_L2_split_tissue",
  rel_path = "figures/umap_L2_split_tissue.png",
  display_reduction = umap_reduction,
  note = "Deterministic redraw on shared display UMAP"
)

p_l3 <- build_umap_plot(
  obj,
  umap_reduction,
  LABEL_COL,
  title = sprintf("Myeloid - Cell Type L3 (%s)", umap_reduction),
  label = TRUE,
  width = 18,
  height = 12
)
save_plot(p_l3$plot, file.path(FIG_DIR, "umap_celltype_L3"), p_l3$width, p_l3$height)
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "umap_celltype_L3",
  rel_path = "figures/umap_celltype_L3.png",
  display_reduction = umap_reduction,
  note = "Deterministic redraw on shared display UMAP"
)

p_choir <- build_umap_plot(
  obj,
  umap_reduction,
  choir_col,
  title = sprintf(
    "Myeloid - CHOIR Clusters (display: %s; clustered on: %s)",
    umap_reduction,
    choir_reduction
  ),
  label = TRUE,
  width = 14,
  height = 10
)
save_plot(p_choir$plot, file.path(FIG_DIR, "choir_umap"), p_choir$width, p_choir$height)
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "choir_umap",
  rel_path = "figures/choir_umap.png",
  display_reduction = umap_reduction,
  cluster_source_reduction = choir_reduction,
  note = "CHOIR clusters projected onto the same deterministic display UMAP as L3"
)

p_choir_split <- build_umap_plot(
  obj,
  umap_reduction,
  choir_col,
  title = sprintf("Myeloid - CHOIR by Tissue (%s)", umap_reduction),
  split_col = TISSUE_COL,
  label = TRUE,
  width = max(12, 4.5 * length(tissues)),
  height = 8
)
save_plot(p_choir_split$plot, file.path(FIG_DIR, "choir_umap_split_tissue"), p_choir_split$width, p_choir_split$height)
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "choir_umap_split_tissue",
  rel_path = "figures/choir_umap_split_tissue.png",
  display_reduction = umap_reduction,
  cluster_source_reduction = choir_reduction,
  note = "CHOIR clusters projected onto the same deterministic display UMAP as L3"
)

# Marker dotplot -------------------------------------------------------------
cat("[STEP] Redrawing marker dotplot...\n")
markers_present <- intersect(known_markers, rownames(obj))
if (length(markers_present) >= 3L) {
  Idents(obj) <- LABEL_COL
  p_dot <- DotPlot(obj, features = markers_present) +
    RotatedAxis() +
    ggtitle("Myeloid - Known Markers (L3)") +
    theme(axis.text.x = element_text(size = 7))
  save_plot(
    p_dot,
    file.path(FIG_DIR, "dotplot_markers"),
    width = max(12, length(markers_present) * 0.45),
    height = max(6, length(unique(meta[[LABEL_COL]])) * 0.4)
  )
  manifest_rows <- record_manifest_row(
    manifest_rows,
    figure = "dotplot_markers",
    rel_path = "figures/dotplot_markers.png",
    note = sprintf("Known marker dotplot across %d present markers", length(markers_present))
  )
}

# Composition plots ----------------------------------------------------------
cat("[STEP] Redrawing composition plots...\n")
comp_df <- meta %>%
  filter(!is.na(.data[[TISSUE_COL]]), !is.na(.data[[CELLTYPE_L2_COL]])) %>%
  count(.data[[TISSUE_COL]], .data[[CELLTYPE_L2_COL]]) %>%
  group_by(.data[[TISSUE_COL]]) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp <- ggplot(comp_df, aes(x = .data[[TISSUE_COL]], y = pct, fill = .data[[CELLTYPE_L2_COL]])) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Tissue",
    y = "Percentage (%)",
    title = "Myeloid - L2 Composition by Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_comp, file.path(FIG_DIR, "composition_tissue_L2"))
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "composition_tissue_L2",
  rel_path = "figures/composition_tissue_L2.png",
  note = "Regenerated from final object metadata"
)

p_count <- ggplot(comp_df, aes(x = .data[[TISSUE_COL]], y = n, fill = .data[[CELLTYPE_L2_COL]])) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Tissue",
    y = "Cell Count",
    title = "Myeloid - Absolute Count by Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_count, file.path(FIG_DIR, "count_tissue_L2"))
manifest_rows <- record_manifest_row(
  manifest_rows,
  figure = "count_tissue_L2",
  rel_path = "figures/count_tissue_L2.png",
  note = "Regenerated from final object metadata"
)

sample_comp <- meta %>%
  filter(
    !is.na(.data[[TISSUE_COL]]),
    !is.na(.data[[CELLTYPE_L2_COL]]),
    !is.na(.data[[SAMPLE_COL]])
  ) %>%
  count(.data[[SAMPLE_COL]], .data[[TISSUE_COL]], .data[[CELLTYPE_L2_COL]]) %>%
  group_by(.data[[SAMPLE_COL]]) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

if (nrow(sample_comp) > 0L) {
  l2_types <- sort(unique(as.character(meta[[CELLTYPE_L2_COL]])))
  facet_ncol <- min(4L, length(l2_types))
  facet_nrow <- ceiling(length(l2_types) / facet_ncol)
  p_sample <- ggplot(sample_comp, aes(x = .data[[TISSUE_COL]], y = pct, fill = .data[[TISSUE_COL]])) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(as.formula(paste("~", CELLTYPE_L2_COL)), scales = "free_y", ncol = facet_ncol) +
    labs(
      x = "Tissue",
      y = "Proportion per Sample (%)",
      title = "Myeloid - Sample-level Composition by Tissue"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  save_plot(
    p_sample,
    file.path(FIG_DIR, "composition_sample_level"),
    width = max(10, facet_ncol * 4),
    height = max(6, facet_nrow * 3.5)
  )
  manifest_rows <- record_manifest_row(
    manifest_rows,
    figure = "composition_sample_level",
    rel_path = "figures/composition_sample_level.png",
    note = "Regenerated from final object metadata"
  )
}

if (isTRUE(REGENERATE_HEAVY_PLOTS)) {
  cat("[STEP] Regenerating heavy plots (marker panels + heatmap)...\n")
  if (exists("tc_generate_marker_panel_visualizations", mode = "function") && length(marker_panels) > 0L) {
    marker_panel_summary <- tc_generate_marker_panel_visualizations(
      obj = obj,
      marker_panels = marker_panels,
      fig_dir = file.path(FIG_DIR, MARKER_PANEL_FIG_SUBDIR),
      group_col = MARKER_PANEL_GROUP_COL,
      lineage_display = "Myeloid",
      reduction_name = umap_reduction,
      report_dir = file.path(OUTPUT_DIR, "reports"),
      panel_prefix = "marker_panel",
      downsample_n = MARKER_PANEL_HEATMAP_CELLS_PER_GROUP,
      verbose = TRUE
    )
    manifest_rows <- record_manifest_row(
      manifest_rows,
      figure = "marker_panels",
      rel_path = "figures/marker_panels",
      display_reduction = umap_reduction,
      note = sprintf(
        "Marker panels regenerated (%d summary rows)",
        ifelse(is.data.frame(marker_panel_summary$summary), nrow(marker_panel_summary$summary), 0L)
      )
    )
  }

  Idents(obj) <- CELLTYPE_L2_COL
  top_mk <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    test.use = "wilcox"
  )
  if (!is.null(top_mk) && nrow(top_mk) > 0L) {
    top10 <- top_mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
    obj_ds <- stratified_downsample(obj, CELLTYPE_L2_COL, HEATMAP_CELLS_PER_TYPE)
    obj_ds <- ScaleData(obj_ds, features = unique(top10$gene), verbose = FALSE)
    p_heatmap <- DoHeatmap(obj_ds, features = unique(top10$gene), size = 3) +
      ggtitle("Myeloid - Top Markers per L2 (downsampled, deterministic redraw)")
    save_plot(p_heatmap, file.path(FIG_DIR, "heatmap_top_markers"), width = 14, height = 10)
    manifest_rows <- record_manifest_row(
      manifest_rows,
      figure = "heatmap_top_markers",
      rel_path = "figures/heatmap_top_markers.png",
      note = sprintf("Top-marker heatmap from %d marker rows", nrow(top_mk))
    )
  }
} else {
  manifest_rows <- record_manifest_row(
    manifest_rows,
    figure = "retained_existing_heavy_plots",
    rel_path = "figures/marker_panels | figures/heatmap_top_markers.png",
    display_reduction = umap_reduction,
    note = "Existing heavy plots were retained to keep the refresh fast; all main summary/UMAP figures were redrawn deterministically."
  )
}

manifest_df <- rbindlist(manifest_rows, fill = TRUE)
fwrite(manifest_df, MANIFEST_PATH, sep = "\t")
cat(sprintf("[OK] Manifest written: %s\n", MANIFEST_PATH))
cat(sprintf("[OK] Refreshed %d figure groups\n", nrow(manifest_df)))
cat("=== Myeloid figure refresh complete ===\n")
