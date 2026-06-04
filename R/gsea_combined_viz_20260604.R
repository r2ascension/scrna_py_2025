#!/usr/bin/env Rscript
# gsea_combined_viz_20260604.R
# ==============================================================================
# Combined visualizations per user request:
#   1. GSEA line plots: UP+DOWN+methods combined into one figure per contrast
#   2. FGSEA bar plots: all levels (per-contrast, celltype, lineage) with methods combined
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

# ── Configuration ─────────────────────────────────────────────────────────────

LINEAGE_ROOTS <- list(
  epithelial  = "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
  bcell       = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
  tnk         = "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508",
  myeloid     = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416"
)

OUTPUT_BASE       <- "/home/h2048/data/R/20260604/gsea_volcano_comprehensive_viz"
GSEA_COMBINED_DIR <- file.path(OUTPUT_BASE, "01b_gsea_combined_figures")
FGSEA_BAR_ALL_DIR <- file.path(OUTPUT_BASE, "02b_fgsea_bar_all_levels")

GS_COLLECTIONS <- c("Hallmark", "KEGG", "GO_BP")  # primary 3 for combined figures

TOP_N_GSEA_PER_DIR  <- 3L   # top N pathways per direction per collection for GSEA panels
TOP_N_BAR_CONTRAST  <- 8L   # top N per contrast bar
TOP_N_BAR_CELLTYPE  <- 15L  # top N per celltype bar
TOP_N_BAR_LINEAGE   <- 25L  # top N per lineage combined bar
MAX_TERM_CHARS      <- 55L

# ── Helpers ───────────────────────────────────────────────────────────────────

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  substr(x, 1L, 180L)
}

# ── Gene set cache ────────────────────────────────────────────────────────────

load_gene_sets <- function() {
  gs_list <- list()
  gs_config <- list(
    Hallmark = list(collection = "H",  subcollection = NULL),
    KEGG     = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
    GO_BP    = list(collection = "C5", subcollection = "GO:BP"),
    GO_CC    = list(collection = "C5", subcollection = "GO:CC"),
    GO_MF    = list(collection = "C5", subcollection = "GO:MF")
  )
  for (nm in names(gs_config)) {
    cfg <- gs_config[[nm]]
    args <- list(species = "Homo sapiens", collection = cfg$collection)
    if (!is.null(cfg$subcollection)) args$subcollection <- cfg$subcollection
    msig <- do.call(msigdbr, args)
    gs_list[[nm]] <- split(msig$gene_symbol, msig$gs_name)
    message(sprintf("  %-8s: %d gene sets", nm, length(gs_list[[nm]])))
  }
  gs_list
}

# ── Part 1: Combined GSEA Line Figure ────────────────────────────────────────

#' Create a combined GSEA figure for one contrast
#' Layout: 2 rows (UP, DOWN) × 3 columns (Hallmark, KEGG, GO_BP) × top N per cell
#' Total panels = 2 * 3 * TOP_N = 18 panels max, arranged via patchwork
generate_combined_gsea_figure <- function(contrast_dir, gs_cache, lineage_name,
                                           level, celltype, contrast) {
  de_path <- file.path(contrast_dir, "DESeq2_results.csv")
  if (!file.exists(de_path)) return(NULL)

  de <- fread(de_path)
  if (!("stat" %in% names(de))) {
    if (all(c("log2FoldChange", "lfcSE") %in% names(de))) {
      de[, stat := log2FoldChange / lfcSE]
    } else return(NULL)
  }
  de <- de[!is.na(stat) & !is.infinite(stat)]
  setorder(de, -stat)
  rnk <- de[["stat"]]
  names(rnk) <- de[["gene"]]

  # Collect all top pathways for this contrast
  all_panels <- list()  # list of ggplots
  panel_labels <- character(0)

  for (direction_dir in c("fgsea_up", "fgsea_down")) {
    fgsea_dir <- file.path(contrast_dir, direction_dir)
    if (!dir.exists(fgsea_dir)) next
    dir_label <- if (direction_dir == "fgsea_up") "UP" else "DOWN"

    for (coll in GS_COLLECTIONS) {
      tsv_path <- file.path(fgsea_dir, paste0(coll, ".tsv"))
      if (!file.exists(tsv_path)) next

      res <- tryCatch(fread(tsv_path), error = function(e) data.table())
      if (nrow(res) == 0L) next

      res <- res[order(-abs(NES))]
      top <- res[seq_len(min(TOP_N_GSEA_PER_DIR, nrow(res)))]

      for (i in seq_len(nrow(top))) {
        row <- top[i]
        pw_name <- row$pathway
        if (!pw_name %in% names(gs_cache[[coll]])) next

        pw_genes <- gs_cache[[coll]][[pw_name]]
        nes_val <- round(row$NES, 2)
        padj_val <- formatC(row$padj, format = "e", digits = 1)

        label <- sprintf("%s %s\n%s\nNES=%.2f padj=%s",
                         dir_label, coll, pw_name, nes_val, padj_val)

        # Truncate long pathway names for panel title
        short_name <- if (nchar(pw_name) > 50) paste0(substr(pw_name, 1, 47), "...") else pw_name

        p <- tryCatch({
          plotEnrichment(pathway = pw_genes, stats = rnk) +
            labs(title = short_name,
                 subtitle = sprintf("%s | %s | NES=%.2f", dir_label, coll, nes_val)) +
            theme_minimal(base_size = 7) +
            theme(
              plot.title    = element_text(face = "bold", size = 7),
              plot.subtitle = element_text(size = 6, color = "grey50"),
              axis.title    = element_text(size = 6),
              axis.text     = element_text(size = 5),
              panel.grid    = element_blank(),
              plot.margin   = margin(2, 2, 2, 2)
            )
        }, error = function(e) NULL)

        if (!is.null(p)) {
          all_panels[[length(all_panels) + 1L]] <- p
          panel_labels <- c(panel_labels, label)
        }
      }
    }
  }

  if (length(all_panels) == 0L) return(NULL)

  # Arrange panels in grid: columns = 3 (Hallmark, KEGG, GO_BP), rows determined by content
  # Use patchwork wrap_plots
  n_cols <- min(3L, length(all_panels))
  combined <- wrap_plots(all_panels, ncol = n_cols) +
    plot_annotation(
      title = paste(celltype, "|", gsub("_", " ", contrast)),
      subtitle = paste(lineage_name, level, "— Top", TOP_N_GSEA_PER_DIR,
                       "pathways per direction × collection"),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8, color = "grey50", hjust = 0.5)
      )
    )

  # Save
  out_dir <- file.path(GSEA_COMBINED_DIR, lineage_name, level, celltype)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir, paste0("gsea_combined_", sanitize_name(contrast), ".pdf"))

  # Dynamic sizing
  n_rows <- ceiling(length(all_panels) / n_cols)
  pdf_w <- n_cols * 4.5
  pdf_h <- n_rows * 3.5 + 1

  ggsave(out_path, combined, width = pdf_w, height = pdf_h, dpi = 200, limitsize = FALSE)
  message("    -> ", out_path, sprintf(" (%d panels)", length(all_panels)))

  length(all_panels)
}

# ── Part 2: FGSEA Bar Plots at All Levels ─────────────────────────────────────

#' Collect all FGSEA results into a data.table
collect_all_fgsea <- function(lineage_root) {
  all_fgsea <- data.table()
  for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
    level_path <- file.path(lineage_root, level_dir_name)
    if (!dir.exists(level_path)) next
    level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

    all_dirs <- list.dirs(level_path, recursive = TRUE, full.names = TRUE)
    fgsea_dirs <- grep("fgsea_up$|fgsea_down$", all_dirs, value = TRUE)

    for (fgsea_dir in fgsea_dirs) {
      tsv_files <- list.files(fgsea_dir, pattern = "\\.tsv$", full.names = TRUE)
      if (length(tsv_files) == 0L) next

      rel <- sub(paste0("^", level_path, "/"), "", fgsea_dir)
      parts <- strsplit(rel, "/")[[1]]
      celltype  <- parts[1]
      contrast  <- parts[2]
      dir_label <- if (grepl("fgsea_up$", fgsea_dir)) "UP" else "DOWN"

      for (tsv_f in tsv_files) {
        coll <- sub("\\.tsv$", "", basename(tsv_f))
        if (!coll %in% c("Hallmark", "KEGG", "GO_BP")) next  # focus on main 3

        dt <- tryCatch(fread(tsv_f), error = function(e) data.table())
        if (nrow(dt) == 0L) next

        keep_cols <- intersect(c("pathway", "pval", "padj", "NES", "ES", "size"),
                               names(dt))
        if (length(keep_cols) == 0L) next
        dt <- dt[, keep_cols, with = FALSE]

        dt[, `:=`(
          celltype   = celltype,
          contrast   = contrast,
          level      = level_label,
          collection = coll,
          direction  = dir_label
        )]
        all_fgsea <- rbind(all_fgsea, dt, fill = TRUE)
      }
    }
  }
  all_fgsea
}

#' Per-contrast combined bar plot: all 3 methods in one figure
plot_contrast_fgsea_bar <- function(fgsea_dt, lineage_name, out_dir) {
  if (nrow(fgsea_dt) == 0L) return(0L)
  n_generated <- 0L

  contrast_groups <- fgsea_dt[, .N, by = .(level, celltype, contrast)]
  for (i in seq_len(nrow(contrast_groups))) {
    grp <- contrast_groups[i]
    sub <- fgsea_dt[level == grp$level & celltype == grp$celltype & contrast == grp$contrast]
    if (nrow(sub) == 0L) next

    # Top N per direction, combining all 3 collections
    summary <- sub[, .(
      NES_max = NES[which.max(abs(NES))],
      best_padj = min(padj, na.rm = TRUE),
      collections = paste(sort(unique(collection)), collapse = "+")
    ), by = .(pathway, direction, collection)]

    # Pick top per direction
    summary <- summary[order(-abs(NES_max))]
    summary <- summary[, head(.SD, TOP_N_BAR_CONTRAST), by = direction]
    if (nrow(summary) == 0L) next

    summary[, term_short := substr(pathway, 1, MAX_TERM_CHARS)]
    if (any(duplicated(summary$term_short))) {
      summary[, term_short := make.unique(term_short, sep = "_")]
    }
    # Create unique ordering term to handle duplicates across collections
    summary[, term_unique := paste0(term_short, " [", collection, "]")]
    summary[, term_unique := factor(term_unique, levels = rev(unique(term_unique[order(NES_max)])))]

    p <- ggplot(summary, aes(x = NES_max, y = term_unique, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
      scale_fill_manual(values = c(UP = "#E74C3C", DOWN = "#2980B9")) +
      facet_wrap(~ collection, scales = "free_y", nrow = 1) +
      labs(
        title = sprintf("%s | %s | %s", grp$celltype, gsub("_", " ", grp$contrast), grp$level),
        subtitle = sprintf("%s — Top %d pathways per direction", lineage_name, TOP_N_BAR_CONTRAST),
        x = "NES", y = NULL
      ) +
      theme_bw(base_size = 8) +
      theme(
        plot.title         = element_text(face = "bold", size = 9),
        plot.subtitle      = element_text(size = 7, color = "grey50"),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 6),
        strip.text         = element_text(face = "bold", size = 8),
        legend.position    = "bottom"
      )

    fname <- paste0("contrast_bar_", sanitize_name(grp$celltype), "_",
                    sanitize_name(grp$contrast), "_", grp$level, ".pdf")
    n_terms <- nrow(summary)
    ggsave(file.path(out_dir, "per_contrast", fname), p,
           width = 16, height = max(5, 0.22 * n_terms + 3), dpi = 200, limitsize = FALSE)
    n_generated <- n_generated + 1L
  }
  n_generated
}

#' Celltype-level combined bar plot: top pathways per celltype across all contrasts
plot_celltype_fgsea_bar <- function(fgsea_dt, lineage_name, out_dir) {
  if (nrow(fgsea_dt) == 0L) return(0L)
  n_generated <- 0L

  for (ct in unique(fgsea_dt$celltype)) {
    ct_dt <- fgsea_dt[celltype == ct]
    if (nrow(ct_dt) == 0L) next

    # Aggregate across contrasts within this celltype
    summary <- ct_dt[, .(
      NES_max = NES[which.max(abs(NES))],
      n_contrasts = uniqueN(contrast),
      best_padj = min(padj, na.rm = TRUE)
    ), by = .(pathway, direction, collection)]

    summary <- summary[order(-abs(NES_max))]
    summary <- summary[, head(.SD, TOP_N_BAR_CELLTYPE), by = .(direction, collection)]
    if (nrow(summary) == 0L) next

    summary[, term_short := substr(pathway, 1, MAX_TERM_CHARS)]
    if (any(duplicated(summary$term_short))) {
      summary[, term_short := make.unique(term_short, sep = "_")]
    }
    summary[, term_unique := paste0(term_short, " [", collection, "]")]
    summary[, term_unique := factor(term_unique, levels = rev(unique(term_unique[order(NES_max)])))]

    p <- ggplot(summary, aes(x = NES_max, y = term_unique, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
      scale_fill_manual(values = c(UP = "#E74C3C", DOWN = "#2980B9")) +
      facet_wrap(~ collection, scales = "free_y", nrow = 1) +
      labs(
        title = sprintf("%s | %s", lineage_name, ct),
        subtitle = sprintf("Top %d pathways per direction×collection across %d contrasts",
                           TOP_N_BAR_CELLTYPE, uniqueN(ct_dt$contrast)),
        x = "NES", y = NULL
      ) +
      theme_bw(base_size = 9) +
      theme(
        plot.title         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 7),
        strip.text         = element_text(face = "bold", size = 8),
        legend.position    = "bottom"
      )

    fname <- paste0("celltype_bar_", sanitize_name(ct), ".pdf")
    n_terms <- nrow(summary)
    ggsave(file.path(out_dir, "per_celltype", fname), p,
           width = 16, height = max(6, 0.22 * n_terms + 3), dpi = 200, limitsize = FALSE)
    n_generated <- n_generated + 1L
  }
  n_generated
}

#' Lineage-level combined bar plot: all 3 methods in ONE figure
plot_lineage_combined_bar <- function(fgsea_dt, lineage_name, out_dir) {
  if (nrow(fgsea_dt) == 0L) return(NULL)

  # Aggregate: best NES per pathway per direction per collection
  summary <- fgsea_dt[, .(
    NES_max = NES[which.max(abs(NES))],
    n_celltypes = uniqueN(celltype),
    n_contrasts = uniqueN(paste(celltype, contrast, sep = "::")),
    best_padj = min(padj, na.rm = TRUE)
  ), by = .(pathway, direction, collection)]

  summary <- summary[order(-abs(NES_max))]
  summary <- summary[, head(.SD, TOP_N_BAR_LINEAGE), by = .(direction, collection)]
  if (nrow(summary) == 0L) return(NULL)

  summary[, term_short := substr(pathway, 1, MAX_TERM_CHARS)]
  if (any(duplicated(summary$term_short))) {
    summary[, term_short := make.unique(term_short, sep = "_")]
  }
  summary[, term_unique := paste0(term_short, " [", collection, "]")]
  summary[, term_unique := factor(term_unique, levels = rev(unique(term_unique[order(NES_max)])))]

  p <- ggplot(summary, aes(x = NES_max, y = term_unique, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    scale_fill_manual(values = c(UP = "#E74C3C", DOWN = "#2980B9")) +
    facet_wrap(~ collection, scales = "free_y", nrow = 1) +
    labs(
      title = sprintf("%s — FGSEA Combined (L2+L3)", lineage_name),
      subtitle = sprintf("Top %d pathways per direction×collection | %d cell types, %d contrast combinations",
                         TOP_N_BAR_LINEAGE,
                         uniqueN(fgsea_dt$celltype),
                         uniqueN(fgsea_dt[, paste(celltype, contrast, sep = "::")])),
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 7),
      strip.text         = element_text(face = "bold", size = 9),
      legend.position    = "bottom"
    )

  fname <- "lineage_combined_bar_all_methods.pdf"
  n_terms <- nrow(summary)
  ggsave(file.path(out_dir, fname), p,
         width = 18, height = max(9, 0.2 * n_terms + 3), dpi = 200, limitsize = FALSE)

  # Also PNG
  png_path <- sub("\\.pdf$", ".png", file.path(out_dir, fname))
  ggsave(png_path, p, width = 18, height = max(9, 0.2 * n_terms + 3), dpi = 150, limitsize = FALSE)

  n_terms
}

# ── Main ──────────────────────────────────────────────────────────────────────

main <- function() {
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Combined GSEA + All-level FGSEA Bar Plots\n")
  cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  dir.create(GSEA_COMBINED_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(FGSEA_BAR_ALL_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("─ Loading MSigDB gene sets ─\n")
  gs_cache <- load_gene_sets()
  cat("\n")

  start_time <- Sys.time()
  total_gsea_combined <- 0L
  total_bar_contrast  <- 0L
  total_bar_celltype  <- 0L
  all_lineage_manifests <- list()

  for (lineage_name in names(LINEAGE_ROOTS)) {
    lineage_root <- LINEAGE_ROOTS[[lineage_name]]
    cat(sprintf("\n━━━ %s ━━━\n", toupper(lineage_name)))

    # ── Part 1: Combined GSEA Figures ──
    cat("  [1/3] Combined GSEA figures...\n")
    lineage_gsea <- 0L
    lineage_gsea_panels <- 0L

    for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
      level_path <- file.path(lineage_root, level_dir_name)
      if (!dir.exists(level_path)) next
      level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

      all_dirs <- list.dirs(level_path, recursive = TRUE, full.names = TRUE)
      fgsea_dirs <- grep("fgsea_up$|fgsea_down$", all_dirs, value = TRUE)
      contrast_dirs <- unique(dirname(fgsea_dirs))

      for (cd in contrast_dirs) {
        rel <- sub(paste0("^", level_path, "/"), "", cd)
        parts <- strsplit(rel, "/")[[1]]
        celltype <- parts[1]
        contrast <- if (length(parts) >= 2) parts[2] else "unknown"

        n_panels <- generate_combined_gsea_figure(
          cd, gs_cache, lineage_name, level_label, celltype, contrast
        )
        if (!is.null(n_panels) && n_panels > 0L) {
          lineage_gsea <- lineage_gsea + 1L
          lineage_gsea_panels <- lineage_gsea_panels + n_panels
        }
      }
    }
    cat(sprintf("    %d combined figures (%d panels total)\n", lineage_gsea, lineage_gsea_panels))
    total_gsea_combined <- total_gsea_combined + lineage_gsea

    # ── Part 2: All-level FGSEA Bar Plots ──
    cat("  [2/3] Collecting FGSEA results...\n")
    fgsea_dt <- collect_all_fgsea(lineage_root)
    cat(sprintf("    %d records collected\n", nrow(fgsea_dt)))

    if (nrow(fgsea_dt) == 0L) {
      cat("    SKIP: No FGSEA results\n")
      next
    }

    bar_out <- file.path(FGSEA_BAR_ALL_DIR, lineage_name)
    dir.create(file.path(bar_out, "per_contrast"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(bar_out, "per_celltype"), showWarnings = FALSE, recursive = TRUE)

    # Per-contrast
    cat("  [2a] Per-contrast bar plots...\n")
    n_contrast <- plot_contrast_fgsea_bar(fgsea_dt, lineage_name, bar_out)
    cat(sprintf("    %d per-contrast bar plots\n", n_contrast))
    total_bar_contrast <- total_bar_contrast + n_contrast

    # Per-celltype
    cat("  [2b] Per-celltype bar plots...\n")
    n_celltype <- plot_celltype_fgsea_bar(fgsea_dt, lineage_name, bar_out)
    cat(sprintf("    %d per-celltype bar plots\n", n_celltype))
    total_bar_celltype <- total_bar_celltype + n_celltype

    # Lineage combined
    cat("  [2c] Lineage combined bar plot...\n")
    n_lineage_terms <- plot_lineage_combined_bar(fgsea_dt, lineage_name, bar_out)
    cat(sprintf("    Combined lineage bar: %d terms\n",
                if (is.null(n_lineage_terms)) 0L else n_lineage_terms))

    # Manifest
    manifest <- data.table(
      lineage = lineage_name,
      n_fgsea_records   = nrow(fgsea_dt),
      n_unique_pathways = uniqueN(fgsea_dt$pathway),
      n_celltypes       = uniqueN(fgsea_dt$celltype),
      n_contrasts       = uniqueN(fgsea_dt[, paste(celltype, contrast, sep = "::")]),
      n_contrast_bars   = n_contrast,
      n_celltype_bars   = n_celltype,
      n_gsea_combined   = lineage_gsea
    )
    all_lineage_manifests[[length(all_lineage_manifests) + 1L]] <- manifest
  }

  # ── Write master manifest ──
  master <- rbindlist(all_lineage_manifests)
  fwrite(master, file.path(OUTPUT_BASE, "combined_viz_manifest.tsv"), sep = "\t")

  # ── Summary ──
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  n_gsea_combined <- length(list.files(GSEA_COMBINED_DIR, pattern = "\\.pdf$", recursive = TRUE))
  n_bar_total <- length(list.files(FGSEA_BAR_ALL_DIR, pattern = "\\.pdf$", recursive = TRUE))

  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("  SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("  Combined GSEA figures : %d\n", n_gsea_combined))
  cat(sprintf("  Bar plots total       : %d (contrast=%d, celltype=%d, lineage=%d)\n",
              n_bar_total, total_bar_contrast, total_bar_celltype, length(LINEAGE_ROOTS)))
  cat(sprintf("  Elapsed               : %.1f min\n", elapsed))
  cat(sprintf("\n  Outputs:\n"))
  cat(sprintf("    %s\n", GSEA_COMBINED_DIR))
  cat(sprintf("    %s\n", FGSEA_BAR_ALL_DIR))
  cat(sprintf("\n  Done: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

main()
