#!/usr/bin/env Rscript
# ==============================================================================
# Tissue-comparison pathway enrichment barplots
# ==============================================================================
# Date: 2026-05-08
# Purpose:
#   Postprocess a tissue-comparison output directory and add compact barplots for
#   ssGSEA top pathways and OFA enrichment tables.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  args <- c(
    "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
    "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
    "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
    "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508"
  )
}

sanitize_file <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  substr(x, 1L, 180L)
}

save_plot <- function(p, path, width = 12, height = 8) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, p, width = width, height = height, dpi = 300)
  invisible(path)
}

read_csv_safe <- function(path) {
  tryCatch(fread(path), error = function(e) data.table())
}

plot_ssgsea_bar <- function(path, out_dir, top_n = 40L) {
  dt <- read_csv_safe(path)
  if (nrow(dt) == 0L || !("pathway" %in% colnames(dt))) return(NULL)
  metric_col <- if ("z_score" %in% colnames(dt)) "z_score" else if ("score" %in% colnames(dt)) "score" else NULL
  if (is.null(metric_col)) return(NULL)
  dt[, metric := suppressWarnings(as.numeric(get(metric_col)))]
  dt <- dt[is.finite(metric)]
  if (nrow(dt) == 0L) return(NULL)
  dt[, abs_metric := abs(metric)]
  dt <- dt[order(-abs_metric)][seq_len(min(.N, top_n))]
  dt[, pathway_short := str_trunc(as.character(pathway), 58)]
  group_col <- if ("group" %in% colnames(dt)) "group" else if ("cluster" %in% colnames(dt)) "cluster" else NULL
  if (!is.null(group_col)) {
    dt[, plot_label := str_trunc(paste0(get(group_col), " | ", pathway_short), 90)]
  } else {
    dt[, plot_label := pathway_short]
  }
  dt[, plot_label := factor(plot_label, levels = rev(unique(plot_label[order(metric)])))]
  if (!("direction" %in% colnames(dt))) dt[, direction := ifelse(metric >= 0, "positive", "negative")]
  title <- sprintf("ssGSEA top pathways: %s", basename(path))
  p <- ggplot(dt, aes(x = metric, y = plot_label, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey40") +
    scale_fill_manual(values = c(positive = "#2563eb", negative = "#dc2626", UP = "#2563eb", DOWN = "#dc2626"), na.value = "#6b7280") +
    labs(title = title, x = metric_col, y = NULL, fill = "direction") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  rel <- gsub("^reports/", "", sub(paste0("^", normalizePath(dirname(dirname(path)), mustWork = FALSE), "/"), "", normalizePath(path, mustWork = FALSE)))
  out_path <- file.path(out_dir, paste0("ssgsea_", sanitize_file(rel), ".png"))
  save_plot(p, out_path, width = 13, height = max(7, 0.22 * nrow(dt) + 2.2))
}

plot_ofa_enrichment_bar <- function(path, out_dir, top_n = 30L) {
  dt <- read_csv_safe(path)
  if (nrow(dt) == 0L) return(NULL)
  term_col <- intersect(c("Description", "term", "pathway"), colnames(dt))[1]
  if (is.na(term_col)) return(NULL)
  if ("neg_log10_padj" %in% colnames(dt)) {
    dt[, metric := suppressWarnings(as.numeric(neg_log10_padj))]
    metric_label <- "-log10(adj.P)"
  } else if ("p.adjust" %in% colnames(dt)) {
    dt[, metric := -log10(pmax(suppressWarnings(as.numeric(`p.adjust`)), 1e-300))]
    metric_label <- "-log10(adj.P)"
  } else {
    return(NULL)
  }
  dt <- dt[is.finite(metric)]
  if (nrow(dt) == 0L) return(NULL)
  dt <- dt[order(-metric)][seq_len(min(.N, top_n))]
  dt[, term_short := str_trunc(as.character(get(term_col)), 72)]
  if ("direction" %in% colnames(dt)) {
    dt[, plot_label := str_trunc(paste0(direction, " | ", term_short), 90)]
  } else {
    dt[, direction := "enriched"]
    dt[, plot_label := term_short]
  }
  dt[, plot_label := factor(plot_label, levels = rev(unique(plot_label[order(metric)])))]
  p <- ggplot(dt, aes(x = metric, y = plot_label, fill = direction)) +
    geom_col(width = 0.72) +
    labs(title = sprintf("OFA enrichment: %s", basename(dirname(path))), x = metric_label, y = NULL, fill = "direction") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  rel <- sub("^reports/", "", sub(paste0("^", normalizePath(dirname(dirname(path)), mustWork = FALSE), "/"), "", normalizePath(path, mustWork = FALSE)))
  out_path <- file.path(out_dir, paste0("ofa_", sanitize_file(rel), ".png"))
  save_plot(p, out_path, width = 12, height = max(6.5, 0.24 * nrow(dt) + 2.2))
}

process_output_dir <- function(output_dir) {
  reports_dir <- file.path(output_dir, "reports")
  if (!dir.exists(reports_dir)) {
    warning(sprintf("reports dir not found: %s", reports_dir))
    return(data.table(output_dir = output_dir, status = "missing_reports", n_plots = 0L))
  }
  out_dir <- file.path(output_dir, "figures", "pathway_enrichment_barplots_20260508")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ssgsea_files <- list.files(
    reports_dir,
    pattern = "ssgsea(_l3)?_top_pathways_.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  ofa_files <- list.files(
    reports_dir,
    pattern = "enrichment_bubble_data\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  plot_paths <- c(
    unlist(lapply(ssgsea_files, plot_ssgsea_bar, out_dir = out_dir), use.names = FALSE),
    unlist(lapply(ofa_files, plot_ofa_enrichment_bar, out_dir = out_dir), use.names = FALSE)
  )
  plot_paths <- plot_paths[!is.na(plot_paths) & nzchar(plot_paths)]

  manifest <- data.table(plot_path = plot_paths)
  manifest[, plot_file := basename(plot_path)]
  fwrite(manifest, file.path(out_dir, "pathway_barplot_manifest.tsv"), sep = "\t")

  readme <- c(
    sprintf("# Pathway enrichment barplots — %s", basename(output_dir)),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    sprintf("- Source reports: `%s`", reports_dir),
    sprintf("- Plots generated: `%d`", length(plot_paths)),
    "",
    "These barplots summarize ssGSEA top-pathway CSV files and OFA enrichment bubble-data CSV files when present.",
    "",
    "## Manifest",
    "",
    "- `pathway_barplot_manifest.tsv`",
    ""
  )
  writeLines(readme, file.path(out_dir, "README.md"))

  data.table(output_dir = output_dir, figure_dir = out_dir, status = "ok", n_plots = length(plot_paths))
}

cat("=== Tissue comparison pathway enrichment barplots (2026-05-08) ===\n")
rows <- rbindlist(lapply(args, process_output_dir), fill = TRUE)
print(rows)
if (any(rows$status != "ok")) quit(save = "no", status = 1L)
