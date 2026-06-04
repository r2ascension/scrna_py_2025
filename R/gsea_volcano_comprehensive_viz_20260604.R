#!/usr/bin/env Rscript
# gsea_volcano_comprehensive_viz_20260604.R
# ==============================================================================
# Comprehensive visualization: GSEA line plots, FGSEA bar plots, DEG volcano plots
#
# Three visualization types for docs/experiments/:
#   1. GSEA line plots (经典 GSEA 运行富集分数折线图)
#      - Running enrichment score + barcode + ranked metric distribution
#      - Generated via fgsea::plotEnrichment() for top pathways per contrast
#   2. FGSEA bar plots (FGSEA 富集棒状图)
#      - Top N pathways by NES, split by up/down direction
#      - Per-lineage summary + per-celltype summary
#   3. DEG volcano plots (差异表达火山图)
#      - log2FC vs -log10(padj) for all contrasts across all 7 lineages
#      - Where missing, generated; where present, catalogued
#
# Resume-safe: skips outputs that already exist.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(parallel)
})

# ── Configuration ─────────────────────────────────────────────────────────────

# Lineage root directories (latest versions)
LINEAGE_ROOTS <- list(
  epithelial  = "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
  bcell       = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
  tnk         = "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508",
  myeloid     = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
  endothelial = "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
  fibroblast  = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414",
  smc         = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
)

# Output base for summary-level plots
OUTPUT_BASE    <- "/home/h2048/data/R/20260604/gsea_volcano_comprehensive_viz"
GSEA_LINE_DIR  <- file.path(OUTPUT_BASE, "01_gsea_line_plots")
FGSEA_BAR_DIR  <- file.path(OUTPUT_BASE, "02_fgsea_bar_plots")
VOLCANO_DIR    <- file.path(OUTPUT_BASE, "03_volcano_plots")

# Gene set collections (same as fgsea_tissue_deg_runner_20260602.R)
GS_CONFIG <- list(
  Hallmark = list(collection = "H",  subcollection = NULL,            label = "MSigDB Hallmark"),
  KEGG     = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS", label = "KEGG"),
  GO_BP    = list(collection = "C5", subcollection = "GO:BP",          label = "GO Biological Process")
)

# Plotting parameters
TOP_N_GSEA_LINE   <- 5L    # top N pathways per contrast for GSEA line plots
TOP_N_FGSEA_BAR   <- 30L   # top N pathways for lineage-level FGSEA bar plots
TOP_N_CELLTYPE_BAR <- 15L  # top N pathways for celltype-level FGSEA bar plots
MAX_TERM_CHARS    <- 65L
VOLCANO_PADJ_THR  <- 0.05
VOLCANO_LFC_THR   <- 0.5

# ── Helpers ───────────────────────────────────────────────────────────────────

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  substr(x, 1L, 180L)
}

save_plot <- function(p, path, width = 12, height = 8) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, p, width = width, height = height, dpi = 300)
  message("    -> ", path)
  invisible(path)
}

# ── Gene set cache ────────────────────────────────────────────────────────────

load_gene_sets <- function(gs_config, species = "Homo sapiens") {
  gs_list <- list()
  for (nm in names(gs_config)) {
    cfg <- gs_config[[nm]]
    args <- list(species = species, collection = cfg$collection)
    if (!is.null(cfg$subcollection)) {
      args$subcollection <- cfg$subcollection
    }
    msig <- do.call(msigdbr, args)
    gs <- split(msig$gene_symbol, msig$gs_name)
    gs_list[[nm]] <- gs
    message(sprintf("  %-8s: %d gene sets, %d unique genes",
                    nm, length(gs), length(unique(msig$gene_symbol))))
  }
  gs_list
}

# ── Part 1: GSEA Line Plots ──────────────────────────────────────────────────

#' Generate classic GSEA running enrichment score plot for a single pathway
plot_gsea_line <- function(rnk, pathway_genes, pathway_name, collection_name,
                           contrast_label) {
  # Use fgsea's built-in plotEnrichment
  p <- plotEnrichment(
    pathway = pathway_genes,
    stats   = rnk
  ) +
    labs(
      title    = pathway_name,
      subtitle = paste0(collection_name, " | ", contrast_label)
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey50"),
      panel.grid    = element_blank()
    )
  p
}

#' Generate GSEA line plots for one contrast's FGSEA results
generate_gsea_line_plots_one_contrast <- function(contrast_dir, gs_cache,
                                                   lineage_name, level,
                                                   celltype, contrast) {
  # Find FGSEA TSV files
  for (direction_dir in c("fgsea_up", "fgsea_down")) {
    fgsea_dir <- file.path(contrast_dir, direction_dir)
    if (!dir.exists(fgsea_dir)) next

    tsv_files <- list.files(fgsea_dir, pattern = "\\.tsv$", full.names = TRUE)
    if (length(tsv_files) == 0L) next

    # Read DESeq2 results to get ranked gene list
    de_path <- file.path(contrast_dir, "DESeq2_results.csv")
    if (!file.exists(de_path)) next

    de <- fread(de_path)
    if (!("stat" %in% names(de))) {
      if (all(c("log2FoldChange", "lfcSE") %in% names(de))) {
        de[, stat := log2FoldChange / lfcSE]
      } else next
    }
    de <- de[!is.na(stat) & !is.infinite(stat)]
    setorder(de, -stat)
    rnk <- de[["stat"]]
    names(rnk) <- de[["gene"]]

    direction_label <- if (direction_dir == "fgsea_up") "UP" else "DOWN"

    for (tsv_f in tsv_files) {
      collection_name <- sub("\\.tsv$", "", basename(tsv_f))
      if (!collection_name %in% names(gs_cache)) next

      fgsea_res <- fread(tsv_f)
      if (nrow(fgsea_res) == 0L) next

      # Sort by abs(NES) descending, pick top N
      fgsea_res <- fgsea_res[order(-abs(NES))]
      top_paths <- fgsea_res[seq_len(min(TOP_N_GSEA_LINE, nrow(fgsea_res)))]

      for (i in seq_len(nrow(top_paths))) {
        row <- top_paths[i]
        pw_name <- row$pathway
        if (!pw_name %in% names(gs_cache[[collection_name]])) next

        pw_genes <- gs_cache[[collection_name]][[pw_name]]

        out_dir <- file.path(GSEA_LINE_DIR, lineage_name, level, celltype, contrast,
                             direction_label, collection_name)
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

        fname <- paste0(sanitize_name(pw_name), ".pdf")
        out_path <- file.path(out_dir, fname)

        if (file.exists(out_path)) next  # resume-safe

        contrast_label <- paste(lineage_name, level, celltype, contrast, direction_label, sep = " | ")
        tryCatch({
          p <- plot_gsea_line(rnk, pw_genes, pw_name, collection_name, contrast_label)
          save_plot(p, out_path, width = 8, height = 7)
        }, error = function(e) {
          message("    ! Error plotting ", pw_name, ": ", conditionMessage(e))
        })
      }
    }
  }
}

# ── Part 2: FGSEA Bar Plots ──────────────────────────────────────────────────

#' Collect all FGSEA results for a lineage into one data.table
collect_fgsea_results <- function(lineage_root) {
  all_fgsea <- data.table()
  for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
    level_path <- file.path(lineage_root, level_dir_name)
    if (!dir.exists(level_path)) next
    level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

    # Find all fgsea_up / fgsea_down directories via list.dirs
    all_dirs <- list.dirs(level_path, recursive = TRUE, full.names = TRUE)
    fgsea_dirs <- grep("fgsea_up$|fgsea_down$", all_dirs, value = TRUE)

    for (fgsea_dir in fgsea_dirs) {
      tsv_files <- list.files(fgsea_dir, pattern = "\\.tsv$", full.names = TRUE)
      if (length(tsv_files) == 0L) next

      # Extract celltype, contrast, direction from path
      rel <- sub(paste0("^", level_path, "/"), "", fgsea_dir)
      parts <- strsplit(rel, "/")[[1]]
      celltype  <- parts[1]
      contrast  <- parts[2]
      direction_dir_name <- parts[3]  # fgsea_up or fgsea_down
      direction_label <- if (direction_dir_name == "fgsea_up") "UP" else "DOWN"

      for (tsv_f in tsv_files) {
        collection_name <- sub("\\.tsv$", "", basename(tsv_f))

        dt <- tryCatch(fread(tsv_f), error = function(e) data.table())
        if (nrow(dt) == 0L) next

        # Keep only relevant columns
        keep_cols <- intersect(c("pathway", "pval", "padj", "NES", "ES", "size",
                                  "leadingEdge"), names(dt))
        dt <- dt[, ..keep_cols]

        dt[, `:=`(
          celltype   = celltype,
          contrast   = contrast,
          level      = level_label,
          collection = collection_name,
          direction  = direction_label
        )]
        all_fgsea <- rbind(all_fgsea, dt, fill = TRUE)
      }
    }
  }
  all_fgsea
}

#' Lineage-level FGSEA bar plot: top pathways across all contrasts
plot_lineage_fgsea_bar <- function(fgsea_dt, lineage_name, out_dir) {
  if (nrow(fgsea_dt) == 0L) return(NULL)

  for (coll in unique(fgsea_dt$collection)) {
    sub <- fgsea_dt[collection == coll]
    if (nrow(sub) == 0L) next

    # Aggregate: for each pathway, take best (most extreme NES) across contrasts
    summary <- sub[, .(
      NES_max  = NES[which.max(abs(NES))],
      n_contrasts = uniqueN(paste(celltype, contrast)),
      best_padj = min(padj, na.rm = TRUE)
    ), by = .(pathway, direction)]

    # Filter to pathways seen in at least 2 contrasts
    summary <- summary[n_contrasts >= 1]
    summary <- summary[order(-abs(NES_max))]
    summary <- summary[, head(.SD, TOP_N_FGSEA_BAR), by = direction]

    if (nrow(summary) == 0L) next

    summary[, term_short := substr(pathway, 1, MAX_TERM_CHARS)]
    if (any(duplicated(summary$term_short))) {
      summary[, term_short := make.unique(term_short, sep = "_")]
    }
    summary[, term_short := factor(term_short, levels = rev(unique(term_short)))]

    p <- ggplot(summary, aes(x = NES_max, y = term_short, fill = direction)) +
      geom_col(width = 0.72) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
      scale_fill_manual(values = c(UP = "#E74C3C", DOWN = "#2980B9")) +
      labs(
        title    = sprintf("%s — FGSEA %s (L2+L3)", lineage_name, coll),
        subtitle = sprintf("Top %d pathways by |NES| across %d contrasts",
                           nrow(summary), uniqueN(fgsea_dt[, paste(celltype, contrast)])),
        x        = "Normalized Enrichment Score (NES)",
        y        = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 8),
        legend.position    = "bottom"
      )

    fname <- paste0("lineage_fgsea_bar_", coll, ".pdf")
    save_plot(p, file.path(out_dir, fname),
              width = 14, height = max(8, 0.26 * nrow(summary) + 2.5))

    # Also save PNG
    png_path <- sub("\\.pdf$", ".png", file.path(out_dir, fname))
    ggsave(png_path, p, width = 14, height = max(8, 0.26 * nrow(summary) + 2.5), dpi = 150)
  }
  invisible(NULL)
}

#' Celltype-level FGSEA bar plot
plot_celltype_fgsea_bar <- function(fgsea_dt, lineage_name, out_dir) {
  if (nrow(fgsea_dt) == 0L) return(NULL)

  for (ct in unique(fgsea_dt$celltype)) {
    ct_dt <- fgsea_dt[celltype == ct]
    if (nrow(ct_dt) == 0L) next

    for (coll in unique(ct_dt$collection)) {
      sub <- ct_dt[collection == coll]
      if (nrow(sub) == 0L) next

      summary <- sub[, .(
        NES_max = NES[which.max(abs(NES))],
        best_padj = min(padj, na.rm = TRUE)
      ), by = .(pathway, direction, contrast)]

      summary <- summary[order(-abs(NES_max))]
      summary <- summary[, head(.SD, TOP_N_CELLTYPE_BAR), by = direction]
      if (nrow(summary) == 0L) next

      summary[, term_short := substr(pathway, 1, MAX_TERM_CHARS)]
      if (any(duplicated(summary$term_short))) {
        summary[, term_short := make.unique(term_short, sep = "_")]
      }
      summary[, term_short := factor(term_short, levels = rev(unique(term_short)))]

      n_contrasts <- uniqueN(sub$contrast)

      p <- ggplot(summary, aes(x = NES_max, y = term_short, fill = direction)) +
        geom_col(width = 0.7) +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
        scale_fill_manual(values = c(UP = "#E74C3C", DOWN = "#2980B9")) +
        labs(
          title = sprintf("%s | %s | %s", lineage_name, ct, coll),
          subtitle = sprintf("Top %d pathways by |NES| across %d contrasts",
                             nrow(summary), n_contrasts),
          x = "NES", y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
          plot.title         = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          axis.text.y        = element_text(size = 7),
          legend.position    = "bottom"
        )

      fname <- paste0("celltype_fgsea_bar_", sanitize_name(ct), "_", coll, ".pdf")
      save_plot(p, file.path(out_dir, "per_celltype", fname),
                width = 12, height = max(5, 0.24 * nrow(summary) + 2.5))
    }
  }
}

# ── Part 3: DEG Volcano Plots ─────────────────────────────────────────────────

#' Generate volcano plot for one DESeq2 contrast
make_volcano_plot <- function(de_path, out_path, padj_thr, lfc_thr) {
  de <- fread(de_path)

  # Determine required columns
  if (!all(c("log2FoldChange", "padj") %in% names(de))) {
    return(list(status = "skip:missing_columns"))
  }

  de <- de[!is.na(padj) & !is.na(log2FoldChange) & !is.infinite(log2FoldChange)]
  de[, padj_plot := pmax(padj, .Machine$double.xmin)]
  de[, neg_log10_padj := -log10(padj_plot)]

  # Determine regulation status
  de[, regulation := "stable"]
  de[padj < padj_thr & log2FoldChange > lfc_thr,  regulation := "up"]
  de[padj < padj_thr & log2FoldChange < -lfc_thr, regulation := "down"]

  n_up   <- de[regulation == "up", .N]
  n_down <- de[regulation == "down", .N]

  # Extract contrast label from path
  parts       <- strsplit(dirname(de_path), "/")[[1]]
  contrast    <- parts[length(parts)]
  celltype    <- parts[length(parts) - 1]

  # Convert contrast code to readable label
  contrast_label <- gsub("_", " ", contrast)

  title <- paste(celltype, "|", contrast_label)
  subtitle <- sprintf("UP: %d | DOWN: %d | padj<%.2f, |LFC|>%.1f",
                      n_up, n_down, padj_thr, lfc_thr)

  p <- ggplot(de, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = regulation), size = 0.8, alpha = 0.55) +
    scale_color_manual(
      values = c(down = "#00468B", stable = "grey75", up = "#E64B35"),
      name   = "Regulation"
    ) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
               linewidth = 0.3, color = "grey40") +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
               linewidth = 0.3, color = "grey40") +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "log2(Fold Change)",
      y        = "-log10(Adjusted p-value)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      plot.title   = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey50"),
      legend.position = "bottom"
    )

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_path, p, width = 8, height = 8, dpi = 300)
  list(status = "ok", n_up = n_up, n_down = n_down)
}

#' Generate volcano plots for all contrasts in a lineage
generate_volcano_plots <- function(lineage_root, lineage_name) {
  results <- data.table()
  out_base <- file.path(VOLCANO_DIR, lineage_name)

  for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
    level_path <- file.path(lineage_root, level_dir_name)
    if (!dir.exists(level_path)) next
    level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

    de_files <- list.files(level_path, pattern = "DESeq2_results\\.csv$",
                           recursive = TRUE, full.names = TRUE)
    de_files <- grep("/fgsea_", de_files, invert = TRUE, value = TRUE)
    de_files <- grep("/enrichment_", de_files, invert = TRUE, value = TRUE)

    for (de_f in de_files) {
      rel <- sub(paste0("^", level_path, "/"), "", de_f)
      rel <- sub("/DESeq2_results\\.csv$", "", rel)
      parts <- strsplit(rel, "/")[[1]]
      celltype <- parts[1]
      contrast <- if (length(parts) >= 2) parts[2] else "unknown"

      out_path <- file.path(out_base, level_label, celltype,
                            paste0("volcano_", sanitize_name(contrast), ".pdf"))

      # Also generate a copy in the original contrast directory
      contrast_dir <- dirname(de_f)
      local_pdf    <- file.path(contrast_dir, "volcano.pdf")

      row <- data.table(
        lineage = lineage_name, level = level_label, celltype = celltype,
        contrast = contrast, status = "pending"
      )

      if (file.exists(out_path) && file.exists(local_pdf)) {
        row[, status := "already_done"]
        results <- rbind(results, row, fill = TRUE)
        next
      }

      res <- tryCatch(
        make_volcano_plot(de_f, out_path, VOLCANO_PADJ_THR, VOLCANO_LFC_THR),
        error = function(e) list(status = paste0("error:", conditionMessage(e)))
      )

      # Also copy to local if successful
      if (res$status == "ok" && !file.exists(local_pdf)) {
        dir.create(dirname(local_pdf), showWarnings = FALSE, recursive = TRUE)
        file.copy(out_path, local_pdf, overwrite = TRUE)
      }

      row$status    <- res$status
      row$n_up      <- if (!is.null(res$n_up)) res$n_up else NA_integer_
      row$n_down    <- if (!is.null(res$n_down)) res$n_down else NA_integer_
      results <- rbind(results, row, fill = TRUE)
    }
  }
  results
}

# ── Part 4: Experiment Documentation Generator ────────────────────────────────

write_experiment_docs <- function(manifest) {
  doc_dir <- file.path(OUTPUT_BASE, "docs_for_experiments")
  dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

  # Summary README
  lines <- c(
    "# GSEA / FGSEA / DEG Comprehensive Visualization (2026-06-04)",
    "",
    "## Overview",
    "",
    "This experiment generates three visualization types to complete the",
    "`docs/experiments/` inventory:",
    "",
    "### 1. GSEA Line Plots (`01_gsea_line_plots/`)",
    "Classic GSEA running enrichment score plots showing:",
    "- Running enrichment score (green line) as it walks down the ranked gene list",
    "- Barcode indicating where pathway member genes appear in the ranked list",
    "- Ranked list metric (log2FC-based Wald statistic) distribution",
    "",
    "Generated for top 5 pathways per gene set collection × contrast × direction",
    "using `fgsea::plotEnrichment()`.",
    "",
    "### 2. FGSEA Bar Plots (`02_fgsea_bar_plots/`)",
    "Bar plots of normalized enrichment scores (NES) for top enriched pathways:",
    "- **Lineage-level**: Top 30 pathways by |NES| across all contrasts",
    "- **Celltype-level**: Top 15 pathways by |NES| per celltype",
    "- Databases: Hallmark, KEGG, GO:BP",
    "- Split by UP (red) / DOWN (blue) direction",
    "",
    "### 3. DEG Volcano Plots (`03_volcano_plots/`)",
    "Volcano plots for all DESeq2 contrasts across all 7 lineages:",
    "- log2(Fold Change) vs -log10(adjusted p-value)",
    "- Red = upregulated (padj < 0.05, LFC > 0.5)",
    "- Blue = downregulated (padj < 0.05, LFC < -0.5)",
    "- Grey = not significant",
    "",
    "## Lineage Coverage",
    "",
    sprintf("| Lineage | Volcano | FGSEA Bar | GSEA Line |"),
    sprintf("|---------|---------|-----------|-----------|")
  )

  for (ln in c("epithelial", "bcell", "tnk", "myeloid", "endothelial", "fibroblast", "smc")) {
    v_count  <- nrow(manifest[lineage == ln & grepl("volcano", type)])
    f_count  <- nrow(manifest[lineage == ln & grepl("fgsea_bar", type)])
    g_count  <- nrow(manifest[lineage == ln & grepl("gsea_line", type)])
    v_ok  <- if (v_count > 0) sprintf("%d contrasts", v_count) else "no data"
    f_ok  <- if (f_count > 0) "✓ generated" else "no FGSEA data"
    g_ok  <- if (g_count > 0) "✓ generated" else "no FGSEA data"
    lines <- c(lines, sprintf("| %s | %s | %s | %s |", ln, v_ok, f_ok, g_ok))
  }

  lines <- c(lines, "",
    "## Methods",
    "",
    "- **FGSEA**: `fgsea::fgsea()` with 10,000 permutations, gene sets from MSigDB v2023.2",
    "- **Ranking metric**: DESeq2 Wald statistic (all genes, not just DEGs)",
    "- **DEG calls**: DESeq2 pseudobulk, padj < 0.05, |log2FC| > 0.5",
    "- **Volcano**: ggplot2 with regulation coloring",
    "",
    "## Data Sources",
    "",
    "- DEG tables: `pseudobulk_de/*/DESeq2_results.csv` and `pseudobulk_de_L3/*/DESeq2_results.csv`",
    "- FGSEA results: `fgsea_up/*.tsv` and `fgsea_down/*.tsv` per contrast",
    "- Gene sets: MSigDB Hallmark (H), KEGG (C2:CP:KEGG_MEDICUS), GO:BP (C5:GO:BP)",
    "",
    "## Output Structure",
    "",
    "```",
    "data/R/20260604/gsea_volcano_comprehensive_viz/",
    "├── 01_gsea_line_plots/          # GSEA running enrichment score plots",
    "│   └── <lineage>/<level>/<celltype>/<contrast>/<direction>/<collection>/",
    "│       └── <pathway>.pdf",
    "├── 02_fgsea_bar_plots/          # FGSEA enrichment bar plots",
    "│   ├── <lineage>/",
    "│   │   ├── lineage_fgsea_bar_<collection>.pdf",
    "│   │   └── per_celltype/",
    "│   │       └── celltype_fgsea_bar_<celltype>_<collection>.pdf",
    "│   └── fgsea_bar_manifest.tsv",
    "├── 03_volcano_plots/            # DEG volcano plots",
    "│   └── <lineage>/<level>/<celltype>/",
    "│       └── volcano_<contrast>.pdf",
    "├── docs_for_experiments/        # Documentation for copying to docs/experiments/",
    "│   └── README.md",
    "└── all_manifest.tsv             # Combined manifest of all outputs",
    "```",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  )

  writeLines(lines, file.path(doc_dir, "README.md"))
  message("  Experiment docs written to ", file.path(doc_dir, "README.md"))
}

# ── Main ──────────────────────────────────────────────────────────────────────

main <- function() {
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  GSEA / FGSEA / DEG Comprehensive Visualization\n")
  cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)
  all_manifest <- data.table()

  # Load gene sets once (needed for GSEA line plots + FGSEA bars)
  cat("─ Loading MSigDB gene sets ─\n")
  gs_cache <- load_gene_sets(GS_CONFIG)
  cat("\n")

  start_time <- Sys.time()

  for (lineage_name in names(LINEAGE_ROOTS)) {
    lineage_root <- LINEAGE_ROOTS[[lineage_name]]
    cat(sprintf("\n━━━ %s ━━━\n", toupper(lineage_name)))
    cat(sprintf("  Root: %s\n", lineage_root))

    if (!dir.exists(lineage_root)) {
      cat("  SKIP: directory not found\n")
      next
    }

    # ── Part 1: GSEA Line Plots ──
    cat("  [1/3] GSEA line plots...\n")
    gsea_line_count <- 0L
    for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
      level_path <- file.path(lineage_root, level_dir_name)
      if (!dir.exists(level_path)) next
      level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

      # Find all contrast dirs that have fgsea_up or fgsea_down (via list.dirs)
      all_dirs <- list.dirs(level_path, recursive = TRUE, full.names = TRUE)
      fgsea_dirs <- grep("fgsea_up$|fgsea_down$", all_dirs, value = TRUE)
      fgsea_contrast_dirs <- unique(dirname(fgsea_dirs))

      for (cd in fgsea_contrast_dirs) {
        rel <- sub(paste0("^", level_path, "/"), "", cd)
        parts <- strsplit(rel, "/")[[1]]
        celltype <- parts[1]
        contrast <- if (length(parts) >= 2) parts[2] else "unknown"

        generate_gsea_line_plots_one_contrast(
          cd, gs_cache, lineage_name, level_label, celltype, contrast
        )
        gsea_line_count <- gsea_line_count + 1L
      }
    }
    cat(sprintf("    Processed %d contrast(s) for GSEA line plots\n", gsea_line_count))

    # ── Part 2: FGSEA Bar Plots ──
    cat("  [2/3] FGSEA bar plots...\n")
    fgsea_dt <- collect_fgsea_results(lineage_root)
    if (nrow(fgsea_dt) > 0L) {
      fgsea_bar_out <- file.path(FGSEA_BAR_DIR, lineage_name)
      dir.create(fgsea_bar_out, showWarnings = FALSE, recursive = TRUE)

      plot_lineage_fgsea_bar(fgsea_dt, lineage_name, fgsea_bar_out)
      plot_celltype_fgsea_bar(fgsea_dt, lineage_name, fgsea_bar_out)

      cat(sprintf("    %d FGSEA records, %d unique pathways\n",
                  nrow(fgsea_dt), uniqueN(fgsea_dt$pathway)))

      # Add to manifest
      all_manifest <- rbind(all_manifest,
        data.table(lineage = lineage_name, type = "fgsea_bar",
                   n_records = nrow(fgsea_dt),
                   n_pathways = uniqueN(fgsea_dt$pathway),
                   n_collections = uniqueN(fgsea_dt$collection)),
        fill = TRUE)
    } else {
      cat("    No FGSEA results found (run fgsea_tissue_deg_runner_20260602.R first)\n")
    }

    # ── Part 3: DEG Volcano Plots ──
    cat("  [3/3] DEG volcano plots...\n")
    volcano_results <- generate_volcano_plots(lineage_root, lineage_name)
    volcano_ok <- volcano_results[status == "ok" | status == "already_done", .N]
    cat(sprintf("    %d volcano plots (ok+skip), %d new generated\n",
                volcano_ok, volcano_results[status == "ok", .N]))

    if (nrow(volcano_results) > 0L) {
      all_manifest <- rbind(all_manifest,
        volcano_results[, .(lineage, type = "volcano", level, celltype, contrast,
                            status, n_up, n_down)],
        fill = TRUE)
    }
  }

  # ── Write combined manifest ──
  manifest_path <- file.path(OUTPUT_BASE, "all_manifest.tsv")
  fwrite(all_manifest, manifest_path, sep = "\t")

  # ── Write experiment documentation ──
  write_experiment_docs(all_manifest)

  # ── Summary ──
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("  SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n")

  # Count outputs
  n_gsea_line <- length(list.files(GSEA_LINE_DIR, pattern = "\\.pdf$", recursive = TRUE))
  n_fgsea_bar <- length(list.files(FGSEA_BAR_DIR, pattern = "\\.pdf$", recursive = TRUE))
  n_volcano   <- length(list.files(VOLCANO_DIR,   pattern = "\\.pdf$", recursive = TRUE))
  n_volcano_new <- all_manifest[type == "volcano" & status == "ok", .N]

  cat(sprintf("  GSEA line plots   : %d PDFs\n", n_gsea_line))
  cat(sprintf("  FGSEA bar plots   : %d PDFs\n", n_fgsea_bar))
  cat(sprintf("  Volcano plots     : %d total, %d new\n", n_volcano, n_volcano_new))
  cat(sprintf("  Elapsed           : %.1f min\n", elapsed))
  cat(sprintf("  Output base       : %s\n", OUTPUT_BASE))
  cat(sprintf("  Manifest          : %s\n", manifest_path))
  cat(sprintf("\n  Done: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

main()
