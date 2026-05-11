#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("/home/h2048/script/R/tissue_comparison_llm_outlier_visual_batch_20260507.R")

REMOVAL_OUTPUT_STEM <- "posthoc_cluster_removal_20260507"
FIGURE_OUTPUT_STEM <- "figures"
FIGURE_DOTPLOT_MAX_GENES <- 18L
FIGURE_FEATUREPLOT_MAX_GENES <- 8L

TARGET_OUTPUT_DIRS <- c(
  "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
  "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun",
  "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414"
)

split_gene_field <- function(x) {
  x <- sanitize_text(x)
  if (length(x) == 0L || all(!nzchar(x))) return(character())
  unique(trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)))
}

read_requested_review_rows <- function(output_dir) {
  review_table <- file.path(
    output_dir,
    "figures",
    "llm_outlier_review_20260507",
    "tables",
    "all_outlier_clusters.tsv"
  )
  if (!file.exists(review_table)) {
    stop(sprintf("Requested-review table not found: %s", review_table))
  }
  dt <- fread(review_table, sep = "\t")
  out <- dt[sanitize_text(user_remove_request) == "yes"]
  if (nrow(out) == 0L) {
    stop(sprintf("No rows with user_remove_request=yes in %s", review_table))
  }
  out
}

build_removal_status <- function(cluster_vec, removed_clusters_chr) {
  vals <- ifelse(
    as.character(cluster_vec) %in% removed_clusters_chr,
    paste0("removed_c", as.character(cluster_vec)),
    "retained"
  )
  level_removed <- paste0("removed_c", removed_clusters_chr)
  factor(vals, levels = c("retained", level_removed))
}

build_removal_group_labels <- function(cluster_ids_chr) {
  stats::setNames(paste0("removed c", cluster_ids_chr), cluster_ids_chr)
}

build_removal_dotplot_group <- function(cluster_vec, removed_clusters_chr) {
  removed_map <- build_removal_group_labels(removed_clusters_chr)
  vals <- ifelse(
    as.character(cluster_vec) %in% removed_clusters_chr,
    unname(removed_map[as.character(cluster_vec)]),
    "retained"
  )
  factor(vals, levels = c(unname(removed_map), "retained"))
}

build_removal_colors <- function(removed_clusters_chr) {
  removed_levels <- paste0("removed_c", removed_clusters_chr)
  removed_cols <- grDevices::hcl.colors(length(removed_levels), palette = "Dark 3")
  c(retained = "#d1d5db", stats::setNames(removed_cols, removed_levels))
}

draw_removed_status_umap <- function(obj, reduction, status_col, output_path, title, subtitle = "") {
  removed_levels <- as.character(levels(obj@meta.data[[status_col]]))
  removed_levels <- removed_levels[removed_levels != "retained"]
  colors <- build_removal_colors(gsub("^removed_c", "", removed_levels))
  p <- DimPlot(
    obj,
    reduction = reduction,
    group.by = status_col,
    cols = colors,
    pt.size = point_size_for_obj(obj),
    shuffle = FALSE
  ) +
    labs(title = title, subtitle = subtitle) +
    coord_equal() +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  save_plot_png(p, output_path, width = 9, height = 7.5)
  p
}

draw_filtered_umap <- function(obj, reduction, output_path, title, subtitle = "") {
  obj$llm_filtered_state <- factor("retained_after_filter", levels = "retained_after_filter")
  p <- DimPlot(
    obj,
    reduction = reduction,
    group.by = "llm_filtered_state",
    cols = c(retained_after_filter = "#2563eb"),
    pt.size = point_size_for_obj(obj),
    shuffle = FALSE
  ) +
    labs(title = title, subtitle = subtitle) +
    coord_equal() +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "none"
    )
  save_plot_png(p, output_path, width = 9, height = 7.5)
  p
}

draw_before_after_umap <- function(before_plot, after_plot, output_path, title) {
  combo <- patchwork::wrap_plots(
    before_plot + theme(legend.position = "none"),
    after_plot + theme(legend.position = "none"),
    ncol = 2
  ) + patchwork::plot_annotation(title = title)
  save_plot_png(combo, output_path, width = 16, height = 7.8)
}

draw_removed_dotplot <- function(obj, group_col, genes, output_path, title, subtitle = "") {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(obj)]
  if (length(genes) == 0L) return(invisible(NULL))

  p <- DotPlot(
    obj,
    features = genes,
    group.by = group_col,
    cols = c("#f3f4f6", "#2563eb"),
    dot.scale = 6,
    scale = FALSE
  ) +
    RotatedAxis() +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      color = "Avg.Exp",
      size = "Pct.Exp"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )

  width <- max(10, 0.42 * length(genes) + 4)
  height <- max(5, 0.5 * length(levels(obj@meta.data[[group_col]])) + 2.2)
  save_plot_png(p, output_path, width = width, height = height)
}

process_one_output_dir <- function(output_dir) {
  cat(sprintf("\n=== posthoc figures | %s ===\n", output_dir))
  removal_dir <- file.path(output_dir, REMOVAL_OUTPUT_STEM)
  summary_path <- file.path(removal_dir, "removal_summary.tsv")
  if (!file.exists(summary_path)) {
    stop(sprintf("Removal summary not found: %s", summary_path))
  }

  removal_summary <- fread(summary_path, sep = "\t")
  if (nrow(removal_summary) != 1L) {
    stop(sprintf("Expected exactly one row in %s", summary_path))
  }

  review_rows <- read_requested_review_rows(output_dir)
  review_rows[, cluster_id_chr := as.character(cluster_id)]
  removed_clusters_chr <- unique(as.character(unlist(strsplit(as.character(removal_summary$remove_clusters[[1]]), ",", fixed = TRUE))))
  review_rows <- review_rows[cluster_id_chr %in% removed_clusters_chr]
  if (nrow(review_rows) == 0L) {
    stop(sprintf("No requested review rows matched removed clusters for %s", output_dir))
  }
  review_rows[, cluster_id_int := suppressWarnings(as.integer(cluster_id_chr))]
  setorder(review_rows, screen_key, cluster_id_int)
  review_rows[, cluster_id_int := NULL]

  original_rds <- find_final_object_path(output_dir)
  filtered_rds <- removal_summary$filtered_rds[[1]]
  cluster_kind <- sanitize_text(removal_summary$cluster_kind[[1]])
  lineage_tag <- sanitize_text(removal_summary$lineage_tag[[1]])

  original_obj <- readRDS(original_rds)
  filtered_obj <- readRDS(filtered_rds)
  original_obj <- original_obj[, sort(colnames(original_obj))]
  filtered_obj <- filtered_obj[, sort(colnames(filtered_obj))]
  original_obj <- ensure_expression_ready(original_obj)
  filtered_obj <- ensure_expression_ready(filtered_obj)

  assign_df <- read_cluster_assignments(output_dir, cluster_kind)
  original_obj <- add_plot_cluster_column(original_obj, assign_df, plot_col = PLOT_CLUSTER_COL)
  filtered_obj <- add_plot_cluster_column(filtered_obj, assign_df, plot_col = PLOT_CLUSTER_COL)
  umap_reduction <- pick_umap_reduction(original_obj, lineage_tag)

  match_original <- match_cluster_vector(colnames(original_obj), assign_df)
  original_cluster_vec <- as.integer(match_original$cluster_vec)
  match_filtered <- match_cluster_vector(colnames(filtered_obj), assign_df)
  original_obj$llm_removal_status <- build_removal_status(original_cluster_vec, removed_clusters_chr)
  original_obj$llm_removed_dotplot_group <- build_removal_dotplot_group(original_cluster_vec, removed_clusters_chr)
  filtered_obj$llm_removal_status <- factor("retained_after_filter", levels = "retained_after_filter")

  figure_dir <- file.path(removal_dir, FIGURE_OUTPUT_STEM)
  umap_dir <- file.path(figure_dir, "umap")
  dotplot_dir <- file.path(figure_dir, "dotplots")
  featureplot_dir <- file.path(figure_dir, "featureplots")
  tables_dir <- file.path(figure_dir, "tables")
  dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dotplot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(featureplot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  top_marker_genes <- unique(unlist(lapply(review_rows$top_marker_genes, split_gene_field), use.names = FALSE))
  contamination_genes <- unique(unlist(lapply(review_rows$contamination_focus_genes, split_gene_field), use.names = FALSE))
  top_marker_genes <- tc_match_features_to_available(top_marker_genes, rownames(original_obj))
  contamination_genes <- tc_match_features_to_available(contamination_genes, rownames(original_obj))
  top_marker_genes <- head(top_marker_genes, FIGURE_DOTPLOT_MAX_GENES)
  contamination_genes <- head(contamination_genes, FIGURE_DOTPLOT_MAX_GENES)
  contamination_feature_genes <- head(contamination_genes, FIGURE_FEATUREPLOT_MAX_GENES)

  before_umap_path <- file.path(umap_dir, "removed_clusters_highlight_before.png")
  after_umap_path <- file.path(umap_dir, "filtered_object_umap.png")
  combined_umap_path <- file.path(umap_dir, "before_after_removal_umap.png")
  top_dotplot_path <- file.path(dotplot_dir, "removed_clusters_top_marker_dotplot.png")
  contamination_dotplot_path <- file.path(dotplot_dir, "removed_clusters_contamination_dotplot.png")
  contamination_featureplot_path <- file.path(featureplot_dir, "removed_clusters_contamination_featureplot.png")
  review_subset_path <- file.path(tables_dir, "requested_removed_outliers.tsv")

  before_plot <- draw_removed_status_umap(
    obj = original_obj,
    reduction = umap_reduction,
    status_col = "llm_removal_status",
    output_path = before_umap_path,
    title = sprintf("%s before removal: clusters requested for removal", lineage_tag),
    subtitle = paste(sprintf("c%s", removed_clusters_chr), collapse = ", ")
  )

  after_plot <- draw_filtered_umap(
    obj = filtered_obj,
    reduction = umap_reduction,
    output_path = after_umap_path,
    title = sprintf("%s after posthoc removal", lineage_tag),
    subtitle = sprintf("retained cells: %d", ncol(filtered_obj))
  )

  draw_before_after_umap(
    before_plot = before_plot,
    after_plot = after_plot,
    output_path = combined_umap_path,
    title = sprintf("%s posthoc cluster removal: before vs after", lineage_tag)
  )

  if (length(top_marker_genes) > 0L) {
    draw_removed_dotplot(
      obj = original_obj,
      group_col = "llm_removed_dotplot_group",
      genes = top_marker_genes,
      output_path = top_dotplot_path,
      title = sprintf("%s removed clusters: top markers vs retained", lineage_tag),
      subtitle = paste(sprintf("c%s", removed_clusters_chr), collapse = ", ")
    )
  }

  if (length(contamination_genes) > 0L) {
    draw_removed_dotplot(
      obj = original_obj,
      group_col = "llm_removed_dotplot_group",
      genes = contamination_genes,
      output_path = contamination_dotplot_path,
      title = sprintf("%s removed clusters: contamination markers vs retained", lineage_tag),
      subtitle = paste(sprintf("c%s", removed_clusters_chr), collapse = ", ")
    )
  }

  if (length(contamination_feature_genes) > 0L) {
    draw_featureplot_grid(
      obj = original_obj,
      genes = contamination_feature_genes,
      reduction = umap_reduction,
      output_path = contamination_featureplot_path,
      title = sprintf("%s removed clusters: contamination FeaturePlot", lineage_tag),
      subtitle = paste(contamination_feature_genes, collapse = ", ")
    )
  }

  fwrite(review_rows, review_subset_path, sep = "\t")

  readme_lines <- c(
    sprintf("# %s posthoc removal figures", lineage_tag),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    sprintf("- Source output: `%s`", output_dir),
    sprintf("- Source filtered object: `%s`", basename(filtered_rds)),
    sprintf("- Cluster assignment used for removal: `%s`", cluster_kind),
    sprintf("- Removed clusters: `%s`", paste(sprintf("c%s", removed_clusters_chr), collapse = ", ")),
    sprintf("- UMAP reduction: `%s`", umap_reduction),
    "",
    "## UMAP",
    "",
    sprintf("![before vs after](umap/%s)", basename(combined_umap_path)),
    "",
    sprintf("![before highlight](umap/%s)", basename(before_umap_path)),
    "",
    sprintf("![after filtered](umap/%s)", basename(after_umap_path)),
    "",
    "## Dotplots",
    "",
    if (file.exists(top_dotplot_path)) sprintf("![top-marker dotplot](dotplots/%s)", basename(top_dotplot_path)) else NULL,
    "",
    if (file.exists(contamination_dotplot_path)) sprintf("![contamination dotplot](dotplots/%s)", basename(contamination_dotplot_path)) else NULL,
    "",
    "## Featureplot",
    "",
    if (file.exists(contamination_featureplot_path)) sprintf("![contamination featureplot](featureplots/%s)", basename(contamination_featureplot_path)) else NULL,
    "",
    "## Tables",
    "",
    sprintf("- [`%s`](tables/%s)", basename(review_subset_path), basename(review_subset_path)),
    ""
  )
  writeLines(Filter(Negate(is.null), readme_lines), file.path(figure_dir, "README.md"))

  data.table(
    lineage_tag = lineage_tag,
    output_dir = output_dir,
    figure_dir = figure_dir,
    removed_clusters = paste(removed_clusters_chr, collapse = ","),
    n_original_cells = ncol(original_obj),
    n_filtered_cells = ncol(filtered_obj),
    n_review_rows = nrow(review_rows),
    before_after_umap = combined_umap_path,
    contamination_dotplot = if (file.exists(contamination_dotplot_path)) contamination_dotplot_path else NA_character_,
    contamination_featureplot = if (file.exists(contamination_featureplot_path)) contamination_featureplot_path else NA_character_
  )
}

cat("=== Tissue Comparison Posthoc Cluster Removal Figures (2026-05-07) ===\n")
rows <- lapply(TARGET_OUTPUT_DIRS, function(output_dir) {
  tryCatch(
    process_one_output_dir(output_dir),
    error = function(e) {
      data.table(
        lineage_tag = basename(output_dir),
        output_dir = output_dir,
        figure_dir = file.path(output_dir, REMOVAL_OUTPUT_STEM, FIGURE_OUTPUT_STEM),
        removed_clusters = NA_character_,
        n_original_cells = NA_integer_,
        n_filtered_cells = NA_integer_,
        n_review_rows = NA_integer_,
        before_after_umap = NA_character_,
        contamination_dotplot = NA_character_,
        contamination_featureplot = NA_character_,
        error = conditionMessage(e)
      )
    }
  )
})
summary_dt <- rbindlist(rows, fill = TRUE)
print(summary_dt)
if ("error" %in% colnames(summary_dt) && any(!is.na(summary_dt$error) & nzchar(summary_dt$error))) {
  stop("One or more posthoc figure targets failed")
}
cat("[DONE] Posthoc removal figures generated successfully.\n")
