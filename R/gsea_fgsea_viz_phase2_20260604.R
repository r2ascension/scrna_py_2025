#!/usr/bin/env Rscript
# gsea_fgsea_viz_phase2_20260604.R
# Phase 2: GSEA line plots + FGSEA bar plots (volcano done in phase 1)
# Fixed: use list.dirs() instead of list.files(include.dirs=TRUE) for directory search

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

# ── Configuration ─────────────────────────────────────────────────────────────

LINEAGE_ROOTS <- list(
  epithelial  = "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
  bcell       = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
  tnk         = "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508",
  myeloid     = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
  endothelial = "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
  fibroblast  = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414",
  smc         = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
)

OUTPUT_BASE    <- "/home/h2048/data/R/20260604/gsea_volcano_comprehensive_viz"
GSEA_LINE_DIR  <- file.path(OUTPUT_BASE, "01_gsea_line_plots")
FGSEA_BAR_DIR  <- file.path(OUTPUT_BASE, "02_fgsea_bar_plots")

GS_CONFIG <- list(
  Hallmark = list(collection = "H",  subcollection = NULL,            label = "MSigDB Hallmark"),
  KEGG     = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS", label = "KEGG"),
  GO_BP    = list(collection = "C5", subcollection = "GO:BP",          label = "GO Biological Process")
)

TOP_N_GSEA_LINE    <- 5L
TOP_N_FGSEA_BAR    <- 30L
TOP_N_CELLTYPE_BAR  <- 15L
MAX_TERM_CHARS     <- 65L

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

# ── Gene sets ─────────────────────────────────────────────────────────────────

load_gene_sets <- function(gs_config, species = "Homo sapiens") {
  gs_list <- list()
  for (nm in names(gs_config)) {
    cfg <- gs_config[[nm]]
    args <- list(species = species, collection = cfg$collection)
    if (!is.null(cfg$subcollection)) args$subcollection <- cfg$subcollection
    msig <- do.call(msigdbr, args)
    gs <- split(msig$gene_symbol, msig$gs_name)
    gs_list[[nm]] <- gs
    message(sprintf("  %-8s: %d gene sets, %d unique genes",
                    nm, length(gs), length(unique(msig$gene_symbol))))
  }
  gs_list
}

# ── Part 1: GSEA Line Plots ──────────────────────────────────────────────────

generate_gsea_line_plots_one_contrast <- function(contrast_dir, gs_cache,
                                                   lineage_name, level,
                                                   celltype, contrast) {
  n_generated <- 0L
  for (direction_dir in c("fgsea_up", "fgsea_down")) {
    fgsea_dir <- file.path(contrast_dir, direction_dir)
    if (!dir.exists(fgsea_dir)) next

    tsv_files <- list.files(fgsea_dir, pattern = "\\.tsv$", full.names = TRUE)
    if (length(tsv_files) == 0L) next

    # Read DESeq2 results for ranked gene list
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

      fgsea_res <- tryCatch(fread(tsv_f), error = function(e) data.table())
      if (nrow(fgsea_res) == 0L) next

      fgsea_res <- fgsea_res[order(-abs(NES))]
      top_paths <- fgsea_res[seq_len(min(TOP_N_GSEA_LINE, nrow(fgsea_res)))]

      for (i in seq_len(nrow(top_paths))) {
        row <- top_paths[i]
        pw_name <- row$pathway
        if (!pw_name %in% names(gs_cache[[collection_name]])) next

        pw_genes <- gs_cache[[collection_name]][[pw_name]]

        out_dir <- file.path(GSEA_LINE_DIR, lineage_name, level, celltype, contrast,
                             direction_label, collection_name)
        fname <- paste0(sanitize_name(pw_name), ".pdf")
        out_path <- file.path(out_dir, fname)

        if (file.exists(out_path)) next  # resume-safe

        contrast_label <- paste(lineage_name, level, celltype, contrast, direction_label, sep = " | ")
        tryCatch({
          p <- plotEnrichment(pathway = pw_genes, stats = rnk) +
            labs(
              title    = pw_name,
              subtitle = paste0(collection_name, " | ", contrast_label)
            ) +
            theme_bw(base_size = 11) +
            theme(
              plot.title    = element_text(face = "bold", size = 11),
              plot.subtitle = element_text(size = 8, color = "grey50"),
              panel.grid    = element_blank()
            )
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          ggsave(out_path, p, width = 8, height = 7, dpi = 200)
          n_generated <- n_generated + 1L
        }, error = function(e) {
          message("    ! Error: ", pw_name, " - ", conditionMessage(e))
        })
      }
    }
  }
  n_generated
}

# ── Part 2: FGSEA Bar Plots ──────────────────────────────────────────────────

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

        keep_cols <- intersect(c("pathway", "pval", "padj", "NES", "ES", "size",
                                  "leadingEdge"), names(dt))
        if (length(keep_cols) == 0L) next
        dt <- dt[, keep_cols, with = FALSE]

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

plot_lineage_fgsea_bar <- function(fgsea_dt, lineage_name, out_dir) {
  if (nrow(fgsea_dt) == 0L) return(NULL)

  for (coll in unique(fgsea_dt$collection)) {
    sub <- fgsea_dt[collection == coll]
    if (nrow(sub) == 0L) next

    summary <- sub[, .(
      NES_max  = NES[which.max(abs(NES))],
      n_contrasts = uniqueN(paste(celltype, contrast, sep = "::")),
      best_padj = min(padj, na.rm = TRUE)
    ), by = .(pathway, direction)]

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
        subtitle = sprintf("Top %d pathways by |NES| across %d contrast×celltype combinations",
                           nrow(summary), uniqueN(sub[, paste(celltype, contrast, sep = "::")])),
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

    # Also PNG for easy viewing
    png_path <- sub("\\.pdf$", ".png", file.path(out_dir, fname))
    ggsave(png_path, p, width = 14, height = max(8, 0.26 * nrow(summary) + 2.5), dpi = 150)
  }
  invisible(NULL)
}

# ── Main ──────────────────────────────────────────────────────────────────────

main <- function() {
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Phase 2: GSEA line + FGSEA bar plots\n")
  cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)
  start_time <- Sys.time()

  cat("─ Loading MSigDB gene sets ─\n")
  gs_cache <- load_gene_sets(GS_CONFIG)
  cat("\n")

  total_gsea_lines <- 0L
  total_fgsea_records <- 0L

  for (lineage_name in names(LINEAGE_ROOTS)) {
    lineage_root <- LINEAGE_ROOTS[[lineage_name]]
    cat(sprintf("\n━━━ %s ━━━\n", toupper(lineage_name)))

    if (!dir.exists(lineage_root)) {
      cat("  SKIP: directory not found\n")
      next
    }

    # ── Part 1: GSEA Line Plots ──
    cat("  [1/2] GSEA line plots...\n")
    lineage_gsea_count <- 0L
    for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
      level_path <- file.path(lineage_root, level_dir_name)
      if (!dir.exists(level_path)) next
      level_label <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"

      all_dirs <- list.dirs(level_path, recursive = TRUE, full.names = TRUE)
      fgsea_dirs <- grep("fgsea_up$|fgsea_down$", all_dirs, value = TRUE)
      fgsea_contrast_dirs <- unique(dirname(fgsea_dirs))

      for (cd in fgsea_contrast_dirs) {
        rel <- sub(paste0("^", level_path, "/"), "", cd)
        parts <- strsplit(rel, "/")[[1]]
        celltype <- parts[1]
        contrast <- if (length(parts) >= 2) parts[2] else "unknown"

        n <- generate_gsea_line_plots_one_contrast(
          cd, gs_cache, lineage_name, level_label, celltype, contrast
        )
        lineage_gsea_count <- lineage_gsea_count + n
      }
    }
    cat(sprintf("    %d GSEA line plots generated\n", lineage_gsea_count))
    total_gsea_lines <- total_gsea_lines + lineage_gsea_count

    # ── Part 2: FGSEA Bar Plots ──
    cat("  [2/2] FGSEA bar plots...\n")
    fgsea_dt <- collect_fgsea_results(lineage_root)
    if (nrow(fgsea_dt) > 0L) {
      fgsea_bar_out <- file.path(FGSEA_BAR_DIR, lineage_name)
      plot_lineage_fgsea_bar(fgsea_dt, lineage_name, fgsea_bar_out)

      cat(sprintf("    %d FGSEA records, %d unique pathways\n",
                  nrow(fgsea_dt), uniqueN(fgsea_dt$pathway)))
      total_fgsea_records <- total_fgsea_records + nrow(fgsea_dt)
    } else {
      cat("    No FGSEA results found\n")
    }
  }

  # ── Summary ──
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  n_gsea_line <- length(list.files(GSEA_LINE_DIR, pattern = "\\.pdf$", recursive = TRUE))
  n_fgsea_bar <- length(list.files(FGSEA_BAR_DIR, pattern = "\\.(pdf|png)$", recursive = TRUE))

  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("  SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("  GSEA line plots   : %d PDFs\n", n_gsea_line))
  cat(sprintf("  FGSEA bar plots   : %d files (PDF+PNG)\n", n_fgsea_bar))
  cat(sprintf("  Elapsed           : %.1f min\n", elapsed))
  cat(sprintf("  Output base       : %s\n", OUTPUT_BASE))
  cat(sprintf("\n  Done: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

main()
