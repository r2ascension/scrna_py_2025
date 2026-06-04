#!/usr/bin/env Rscript
# ==============================================================================
# Tissue-comparison ORA enrichment barplots (GO_BP + KEGG from per-contrast results)
# ==============================================================================
# Date: 2026-06-03
# Purpose:
#   Generate bar plots from per-contrast ORA enrichment results (GO_BP.csv, KEGG.csv
#   in enrichment_up/ and enrichment_down/) for lineages that are missing them.
#   Uses traditional ORA (clusterProfiler), not ssGSEA.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ---- config ----
TOP_N_PER_CONTRAST  <- 5L    # top terms per contrast (up) + per contrast (down)
TOP_N_PER_CELLTYPE  <- 15L   # top terms per cell type summary
TOP_N_PER_LINEAGE   <- 25L   # top terms across whole lineage
MAX_TERM_CHARS      <- 65L

# ---- helpers ----

sanitize_file <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  substr(x, 1L, 180L)
}

save_plot <- function(p, path, width = 14, height = 9) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, p, width = width, height = height, dpi = 300)
  message("  -> ", path)
  invisible(path)
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(data.table())
  tryCatch(fread(path), error = function(e) data.table())
}

# ---- barplot generators ----

#' Per-contrast barplot: top GO_BP + KEGG up/down for a single contrast
plot_contrast_bar <- function(enrich_up_dir, enrich_down_dir, label, out_dir) {
  plots <- list()

  for (db in c("GO_BP", "KEGG")) {
    for (dir_side in c("up", "down")) {
      direction <- if (dir_side == "up") "UP" else "DOWN"
      base_dir  <- if (dir_side == "up") enrich_up_dir else enrich_down_dir
      csv_path  <- file.path(base_dir, paste0(db, ".csv"))

      dt <- read_csv_safe(csv_path)
      if (nrow(dt) == 0L) next

      # ensure required columns
      if (!("Description" %in% colnames(dt))) next
      if (!("p.adjust" %in% colnames(dt)) && !("pvalue" %in% colnames(dt))) next

      pval_col <- if ("p.adjust" %in% colnames(dt)) "p.adjust" else "pvalue"
      dt[, neg_log10_p := -log10(pmax(suppressWarnings(as.numeric(get(pval_col))), 1e-300))]
      dt <- dt[is.finite(neg_log10_p)][order(-neg_log10_p)]
      if (nrow(dt) == 0L) next

      dt <- dt[seq_len(min(.N, TOP_N_PER_CONTRAST))]
      dt[, term_short := str_trunc(as.character(Description), MAX_TERM_CHARS)]
      dt[, term_short := factor(term_short, levels = rev(unique(term_short[order(neg_log10_p)])))]

      title <- sprintf("%s | %s %s | %s", label, db, direction, basename(dirname(base_dir)))

      p <- ggplot(dt, aes(x = neg_log10_p, y = term_short)) +
        geom_col(fill = if (direction == "UP") "#2563eb" else "#dc2626", width = 0.7) +
        labs(title = title, x = "-log10(p.adjust)", y = NULL) +
        theme_bw(base_size = 10) +
        theme(
          plot.title    = element_text(face = "bold", size = 10),
          panel.grid.major.y = element_blank(),
          axis.text.y   = element_text(size = 8)
        )

      fname <- paste0("contrast_", sanitize_file(label), "_", db, "_", direction, ".png")
      save_plot(p, file.path(out_dir, "per_contrast", fname),
                width = 14, height = max(3, 0.28 * nrow(dt) + 2.5))
      plots[[length(plots) + 1L]] <- fname
    }
  }
  return(plots)
}

#' Per-celltype summary barplot: aggregate top terms across all contrasts for a cell type
plot_celltype_summary <- function(contrast_results, out_dir) {
  if (nrow(contrast_results) == 0L) return(NULL)

  # group by celltype × level × db × direction, pick top across all contrasts
  summary <- contrast_results[, .(
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    n_contrasts     = uniqueN(contrast)
  ), by = .(celltype, level, db, direction, Description)]

  summary <- summary[order(-max_neg_log10_p)]
  summary <- summary[, head(.SD, TOP_N_PER_CELLTYPE), by = .(celltype, level, db, direction)]
  if (nrow(summary) == 0L) return(NULL)

  for (ct in unique(summary$celltype)) {
    sub <- summary[celltype == ct]
    sub[, term_short := str_trunc(as.character(Description), MAX_TERM_CHARS)]
    sub[, panel_label := paste0(db, " | ", direction)]

    sub[, term_short := factor(term_short,
      levels = rev(unique(term_short[order(max_neg_log10_p)])))]

    title <- sprintf("ORA top terms: %s (L2+L3)", ct)

    p <- ggplot(sub, aes(x = max_neg_log10_p, y = term_short, fill = direction)) +
      geom_col(width = 0.7) +
      facet_wrap(~ panel_label, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = c(UP = "#2563eb", DOWN = "#dc2626")) +
      labs(title = title, x = "-log10(p.adjust)", y = NULL,
           subtitle = sprintf("Max across %d contrasts", uniqueN(sub$n_contrasts))) +
      theme_bw(base_size = 10) +
      theme(
        plot.title         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 7.5),
        strip.text         = element_text(face = "bold", size = 9)
      )

    fname <- paste0("celltype_summary_", sanitize_file(ct), ".png")
    save_plot(p, file.path(out_dir, "per_celltype", fname),
              width = 16, height = max(6, 0.24 * nrow(sub) + 3.5))
  }
  invisible(NULL)
}

#' Lineage-level summary: top GO_BP + KEGG terms across the whole lineage
plot_lineage_summary <- function(contrast_results, out_dir, lineage_name) {
  if (nrow(contrast_results) == 0L) return(NULL)

  for (db in c("GO_BP", "KEGG")) {
    sub <- contrast_results[db == db]
    if (nrow(sub) == 0L) next

    # Aggregate: for each term, take the best -log10(padj) across all contrasts
    lineage_summary <- sub[, .(
      max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
      n_celltypes     = uniqueN(celltype)
    ), by = .(direction, Description)]

    lineage_summary <- lineage_summary[order(-max_neg_log10_p)]
    lineage_summary <- lineage_summary[, head(.SD, TOP_N_PER_LINEAGE), by = direction]
    if (nrow(lineage_summary) == 0L) next

    lineage_summary[, term_short := str_trunc(as.character(Description), MAX_TERM_CHARS)]
    lineage_summary[, term_short := factor(term_short,
      levels = rev(unique(term_short[order(max_neg_log10_p)])))]

    title <- sprintf("%s — %s ORA Enrichment (L2+L3 combined)", lineage_name, db)

    p <- ggplot(lineage_summary, aes(x = max_neg_log10_p, y = term_short, fill = direction)) +
      geom_col(width = 0.72) +
      geom_vline(xintercept = -log10(0.05), linewidth = 0.3, color = "grey50", linetype = "dashed") +
      scale_fill_manual(values = c(UP = "#2563eb", DOWN = "#dc2626")) +
      labs(title = title, x = "-log10(p.adjust)", y = NULL) +
      theme_bw(base_size = 11) +
      theme(
        plot.title         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 8)
      )

    fname <- paste0("lineage_summary_", db, ".png")
    save_plot(p, file.path(out_dir, fname),
              width = 15, height = max(7, 0.26 * nrow(lineage_summary) + 2.5))
  }
  invisible(NULL)
}

# ---- main processing ----

process_lineage <- function(lineage_dir) {
  lineage_name <- basename(lineage_dir)
  short_name   <- sub("_tissue_comparison.*", "", lineage_name)

  message("\n========================================")
  message("Processing: ", short_name)
  message("Dir: ", lineage_dir)

  out_dir <- file.path(lineage_dir, "figures", "ora_enrichment_barplots_20260603")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "per_contrast"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "per_celltype"), recursive = TRUE, showWarnings = FALSE)

  # Collect all contrasts with enrichment results
  all_results <- data.table()

  for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
    level_path <- file.path(lineage_dir, level_dir_name)
    if (!dir.exists(level_path)) next

    level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

    # find all DESeq2_results.csv
    deseq_files <- list.files(level_path, pattern = "DESeq2_results\\.csv$",
                              recursive = TRUE, full.names = TRUE)
    message("  ", level_label, ": ", length(deseq_files), " contrasts found")

    for (df in deseq_files) {
      contrast_dir <- dirname(df)
      rel_parts    <- strsplit(sub(paste0("^", level_path, "/"), "", contrast_dir), "/")[[1]]
      celltype     <- rel_parts[1]
      contrast     <- if (length(rel_parts) >= 2) rel_parts[2] else "unknown"

      enrich_up   <- file.path(contrast_dir, "enrichment_up")
      enrich_down <- file.path(contrast_dir, "enrichment_down")

      if (!dir.exists(enrich_up) && !dir.exists(enrich_down)) next

      label <- paste0(short_name, "_", level_label, "_", celltype, "_", contrast)

      # Generate per-contrast barplots
      plot_contrast_bar(enrich_up, enrich_down, label, out_dir)

      # Collect data for summary plots
      for (db in c("GO_BP", "KEGG")) {
        for (dir_side in c("up", "down")) {
          direction <- if (dir_side == "up") "UP" else "DOWN"
          base_dir  <- if (dir_side == "up") enrich_up else enrich_down
          csv_path  <- file.path(base_dir, paste0(db, ".csv"))

          dt <- read_csv_safe(csv_path)
          if (nrow(dt) == 0L) next
          if (!("Description" %in% colnames(dt))) next

          pval_col <- if ("p.adjust" %in% colnames(dt)) "p.adjust" else "pvalue"
          dt[, neg_log10_p := -log10(pmax(suppressWarnings(as.numeric(get(pval_col))), 1e-300))]
          dt <- dt[is.finite(neg_log10_p)]
          if (nrow(dt) == 0L) next

          dt[, `:=`(
            celltype  = celltype,
            level     = level_label,
            db        = db,
            direction = direction,
            contrast  = contrast,
            lineage   = short_name
          )]

          keep_cols <- intersect(c("celltype", "level", "db", "direction", "contrast",
                                   "lineage", "Description", "neg_log10_p",
                                   "Count", "GeneRatio", "p.adjust"),
                                 colnames(dt))
          all_results <- rbind(all_results, dt[, ..keep_cols], fill = TRUE)
        }
      }
    }
  }

  if (nrow(all_results) == 0L) {
    message("  WARNING: No enrichment results found for ", short_name)
    return(data.table(lineage = short_name, status = "no_results", n_plots = 0L))
  }

  # Generate celltype-level summary barplots
  message("  Generating per-celltype summary barplots...")
  plot_celltype_summary(all_results, out_dir)

  # Generate lineage-level summary barplots
  message("  Generating lineage-level summary barplots...")
  plot_lineage_summary(all_results, out_dir, short_name)

  # Count total plots
  n_plots <- length(list.files(out_dir, pattern = "\\.png$", recursive = TRUE))

  # Write manifest
  manifest <- data.table(
    lineage   = short_name,
    n_contrast_plots = length(list.files(file.path(out_dir, "per_contrast"), pattern = "\\.png$")),
    n_celltype_plots = length(list.files(file.path(out_dir, "per_celltype"), pattern = "\\.png$")),
    n_lineage_plots  = length(list.files(out_dir, pattern = "lineage_summary.*\\.png$")),
    total_plots      = n_plots
  )

  fwrite(all_results[, .(
    lineage, level, celltype, contrast, db, direction,
    n_terms = .N,
    top_term = Description[which.max(neg_log10_p)],
    top_neg_log10_p = max(neg_log10_p, na.rm = TRUE)
  ), by = .(lineage, level, celltype, contrast, db, direction)],
    file.path(out_dir, "ora_enrichment_summary.tsv"), sep = "\t")

  fwrite(manifest, file.path(out_dir, "ora_barplot_manifest.tsv"), sep = "\t")

  # Write README
  readme <- c(
    sprintf("# ORA Enrichment Barplots — %s", short_name),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    "## Method",
    "- **ORA (Over-Representation Analysis)** via clusterProfiler",
    "- Databases: GO_BP, KEGG",
    "- Performed separately for UP and DOWN DEGs in each tissue contrast",
    "- Bar plots show -log10(p.adjust) for top enriched terms",
    "",
    "## Output Structure",
    "- `per_contrast/` — Bar plots for each individual DEG contrast",
    "- `per_celltype/` — Summary bar plots aggregating top terms per cell type",
    "- `lineage_summary_GO_BP.png` / `lineage_summary_KEGG.png` — Lineage-level combined",
    "",
    sprintf("## Summary"),
    sprintf("- Total plots: %d", n_plots),
    sprintf("- Contrast-level: %d plots", manifest$n_contrast_plots),
    sprintf("- Celltype-level: %d plots", manifest$n_celltype_plots),
    sprintf("- Lineage-level: %d plots", manifest$n_lineage_plots),
    "",
    "## Data Files",
    "- `ora_enrichment_summary.tsv` — Summary statistics per contrast/database/direction",
    "- `ora_barplot_manifest.tsv` — Plot manifest",
    ""
  )
  writeLines(readme, file.path(out_dir, "README.md"))

  message("  Done: ", n_plots, " plots generated")
  data.table(lineage = short_name, status = "ok", n_plots = n_plots, out_dir = out_dir)
}

# ---- run ----

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  args <- c(
    "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
    "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414",
    "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
  )
}

cat("=== Tissue comparison ORA enrichment barplots (2026-06-03) ===\n")
cat("Input lineages:", length(args), "\n\n")

results <- rbindlist(lapply(args, process_lineage), fill = TRUE)
cat("\n=== Summary ===\n")
print(results)
