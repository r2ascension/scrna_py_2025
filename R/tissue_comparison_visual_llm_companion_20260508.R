#!/usr/bin/env Rscript
# ==============================================================================
# Tissue-comparison visual + LLM companion postprocessor
# ==============================================================================
# Date: 2026-05-08
# Purpose:
#   Add method-level visualizations and standardized LLM companion artifacts to
#   tissue-comparison outputs that already contain raw result tables/figures but
#   do not yet expose paired *_data.csv, *_LLM_PROMPT.md, *_LLM_ANALYSIS.md,
#   *_LLM_STATUS.json, and *_figure_companion_manifest.tsv files.
#
# This is intentionally a light-weight postprocessor: it does not rerun Seurat,
# DESeq2, ssGSEA, OFA, CHOIR, or any expensive analysis.
# ============================================================================== 

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

HELPER <- "/home/h2048/script/R/smc_anno_visual_llm_helper_20260506.R"
if (!file.exists(HELPER)) stop(sprintf("Missing visual LLM helper: %s", HELPER), call. = FALSE)
source(HELPER)

DEFAULT_OUTPUT_DIRS <- c(
  "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
  "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
  "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
  "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) args <- DEFAULT_OUTPUT_DIRS
args <- normalizePath(args, winslash = "/", mustWork = FALSE)

sanitize_file <- function(x, max_len = 160L) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  substr(x, 1L, max_len)
}

shorten <- function(x, width = 80L) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  too_long <- nchar(x) > width
  x[too_long] <- paste0(substr(x[too_long], 1L, max(1L, width - 1L)), "…")
  x
}

factor_unique <- function(labels, reverse = FALSE) {
  labels <- as.character(labels)
  lev <- unique(labels)
  if (isTRUE(reverse)) lev <- rev(lev)
  factor(labels, levels = lev)
}

read_dt <- function(path) {
  tryCatch(data.table::fread(path, data.table = TRUE, showProgress = FALSE), error = function(e) data.table())
}

safe_write_dt <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(dt), path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

save_plot <- function(p, path, width = 12, height = 7) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 300, limitsize = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

existing_plot_path <- function(paths) {
  paths <- as.character(paths)
  paths[file.exists(paths)][1]
}

first_col <- function(dt, candidates) {
  hit <- intersect(candidates, colnames(dt))
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

as_num <- function(x) suppressWarnings(as.numeric(x))

register_safe <- function(figure_path,
                          data,
                          method,
                          figure_type,
                          title,
                          manifest_path,
                          extra_context = NULL,
                          max_rows_for_prompt = 5000L) {
  if (is.null(figure_path) || !length(figure_path) || is.na(figure_path) || !nzchar(figure_path) || !file.exists(figure_path)) {
    return(list(status = "missing_figure", figure_path = figure_path))
  }
  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) > max_rows_for_prompt) {
    df <- df[seq_len(max_rows_for_prompt), , drop = FALSE]
  }
  tryCatch({
    smcanno_viz_register_figure(
      figure_path = figure_path,
      data = df,
      method = method,
      figure_type = figure_type,
      title = title,
      extra_context = extra_context,
      llm_config = list(enabled = TRUE, mode = "queued", model = "deepseek-reasoner", timeout_sec = 240),
      manifest_path = manifest_path
    )
  }, error = function(e) {
    warning(sprintf("Failed to register companion for %s: %s", figure_path, conditionMessage(e)), call. = FALSE)
    list(status = "error", error = conditionMessage(e), figure_path = figure_path)
  })
}

parse_de_path <- function(path, output_dir) {
  rel <- sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  level <- if (length(parts) >= 1L) parts[[1]] else NA_character_
  cluster <- if (length(parts) >= 2L) parts[[2]] else NA_character_
  contrast <- if (length(parts) >= 3L) parts[[3]] else NA_character_
  list(level = level, cluster = cluster, contrast = contrast, rel = rel)
}

summarize_de_file <- function(path, output_dir) {
  dt <- read_dt(path)
  meta <- parse_de_path(path, output_dir)
  if (nrow(dt) == 0L) {
    return(data.table(
      source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
      level = meta$level, cluster = meta$cluster, contrast = meta$contrast,
      n_genes = 0L, n_sig = 0L, n_up = 0L, n_down = 0L,
      top_up_gene = NA_character_, top_down_gene = NA_character_, min_padj = NA_real_
    ))
  }
  gene_col <- first_col(dt, c("gene", "gene_symbol", "symbol", "Gene", "feature", "rowname"))
  if (is.na(gene_col)) {
    gene_col <- "gene"
    dt[, gene := if (!is.null(rownames(dt))) rownames(dt) else as.character(seq_len(.N))]
  }
  padj_col <- first_col(dt, c("padj", "p_val_adj", "p.adjust", "adj.P.Val", "qvalue", "FDR"))
  p_col <- first_col(dt, c("pvalue", "p_val", "P.Value", "p.value"))
  fc_col <- first_col(dt, c("log2FoldChange", "avg_log2FC", "avg_logFC", "logFC", "log2FC"))
  padj <- if (!is.na(padj_col)) as_num(dt[[padj_col]]) else if (!is.na(p_col)) p.adjust(as_num(dt[[p_col]]), method = "BH") else rep(NA_real_, nrow(dt))
  fc <- if (!is.na(fc_col)) as_num(dt[[fc_col]]) else rep(NA_real_, nrow(dt))
  sig <- is.finite(padj) & padj < 0.05 & is.finite(fc) & abs(fc) >= 0.25
  up <- sig & fc > 0
  down <- sig & fc < 0
  top_up_gene <- NA_character_
  top_down_gene <- NA_character_
  if (any(up, na.rm = TRUE)) {
    idx <- which(up)[order(padj[up], -fc[up], na.last = NA)][1]
    top_up_gene <- as.character(dt[[gene_col]][idx])
  }
  if (any(down, na.rm = TRUE)) {
    idx <- which(down)[order(padj[down], fc[down], na.last = NA)][1]
    top_down_gene <- as.character(dt[[gene_col]][idx])
  }
  data.table(
    source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
    volcano_pdf = normalizePath(file.path(dirname(path), "volcano.pdf"), winslash = "/", mustWork = FALSE),
    volcano_png = normalizePath(file.path(dirname(path), "volcano.png"), winslash = "/", mustWork = FALSE),
    level = meta$level,
    cluster = meta$cluster,
    contrast = meta$contrast,
    n_genes = nrow(dt),
    n_sig = sum(sig, na.rm = TRUE),
    n_up = sum(up, na.rm = TRUE),
    n_down = sum(down, na.rm = TRUE),
    top_up_gene = top_up_gene,
    top_down_gene = top_down_gene,
    min_padj = suppressWarnings(min(padj, na.rm = TRUE))
  )
}

summarize_enrichment_file <- function(path, output_dir) {
  dt <- read_dt(path)
  if (nrow(dt) == 0L) return(NULL)
  rel <- sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  term_col <- first_col(dt, c("Description", "term", "pathway", "ID", "GeneSet", "gs_name"))
  padj_col <- first_col(dt, c("p.adjust", "padj", "qvalue", "FDR", "adj.P.Val"))
  p_col <- first_col(dt, c("pvalue", "p.value", "P.Value"))
  if (is.na(term_col)) return(NULL)
  padj <- if (!is.na(padj_col)) as_num(dt[[padj_col]]) else if (!is.na(p_col)) p.adjust(as_num(dt[[p_col]]), method = "BH") else rep(NA_real_, nrow(dt))
  metric <- -log10(pmax(padj, 1e-300))
  keep <- is.finite(metric)
  if (!any(keep)) return(NULL)
  db <- sub("\\.csv$", "", basename(path))
  direction <- if (any(parts %in% c("enrichment_up", "enrichment_down"))) parts[parts %in% c("enrichment_up", "enrichment_down")][1] else NA_character_
  contrast <- if (length(parts) >= 3L) parts[[3]] else NA_character_
  cluster <- if (length(parts) >= 2L) parts[[2]] else NA_character_
  level <- if (length(parts) >= 1L) parts[[1]] else NA_character_
  out <- data.table(
    source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
    level = level,
    cluster = cluster,
    contrast = contrast,
    direction = direction,
    database = db,
    term = as.character(dt[[term_col]]),
    padj = padj,
    neg_log10_padj = metric
  )
  out[keep][order(-neg_log10_padj)]
}

summarize_ssgsea_file <- function(path, level_label, db_label) {
  dt <- read_dt(path)
  if (nrow(dt) == 0L) return(data.table())
  metric_col <- first_col(dt, c("z_score", "score", "mean_z", "delta", "NES", "enrichment_score"))
  pathway_col <- first_col(dt, c("pathway", "term", "Description", "gene_set"))
  if (is.na(metric_col) || is.na(pathway_col)) return(data.table())
  group_col <- first_col(dt, c("group", "cluster", "cell_type", "celltype", "label"))
  if (is.na(group_col)) {
    dt[, group := "all"]
    group_col <- "group"
  }
  dt[, metric := as_num(get(metric_col))]
  dt <- dt[is.finite(metric)]
  if (nrow(dt) == 0L) return(data.table())
  dt[, `:=`(
    source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
    method_level = level_label,
    database = db_label,
    pathway = as.character(get(pathway_col)),
    group = as.character(get(group_col)),
    abs_metric = abs(metric),
    direction = ifelse(metric >= 0, "positive", "negative")
  )]
  dt[order(-abs_metric)]
}

collect_ssgsea_llm_index <- function(output_dir) {
  files <- list.files(file.path(output_dir, "reports"), pattern = "interpretation\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0L) return(data.table())
  rows <- lapply(files, function(path) {
    rel <- sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
    parent <- basename(dirname(path))
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    cls <- paste(class(obj), collapse = ";")
    txt <- NA_character_
    if (is.character(obj)) txt <- paste(head(obj, 3), collapse = "\n")
    if (is.list(obj)) {
      flat <- unlist(obj, recursive = TRUE, use.names = TRUE)
      flat <- flat[nzchar(as.character(flat))]
      if (length(flat) > 0L) txt <- paste(head(paste(names(flat), flat, sep = ": "), 5), collapse = "\n")
    }
    data.table(source_file = normalizePath(path, winslash = "/", mustWork = FALSE), rel_path = rel, group_id = parent, object_class = cls, preview = txt)
  })
  rbindlist(rows, fill = TRUE)
}

plot_count_bar <- function(dt, x_col, y_col, fill_col = NULL, title, xlab = NULL, ylab = NULL) {
  p <- ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]]))
  if (!is.null(fill_col) && fill_col %in% colnames(dt)) {
    p <- p + geom_col(aes(fill = .data[[fill_col]]), width = 0.72)
  } else {
    p <- p + geom_col(fill = "#2563eb", width = 0.72)
  }
  p + coord_flip() + theme_bw(base_size = 11) + labs(title = title, x = xlab, y = ylab, fill = fill_col) +
    theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())
}

make_core_figure_inventory <- function(output_dir, comp_dir, manifest_path) {
  fig_dir <- file.path(output_dir, "figures")
  files <- if (dir.exists(fig_dir)) list.files(fig_dir, pattern = "\\.(png|pdf)$", recursive = TRUE, full.names = TRUE) else character()
  if (length(files) == 0L) return(data.table(method = "core_figures", status = "missing", n_records = 0L, n_figures = 0L))
  dt <- data.table(
    figure_path = normalizePath(files, winslash = "/", mustWork = FALSE),
    rel_path = sub(paste0("^", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "/"), "", normalizePath(files, winslash = "/", mustWork = FALSE)),
    file_ext = tools::file_ext(files),
    bytes = file.info(files)$size
  )
  dt[, figure_group := fifelse(grepl("umap", rel_path, ignore.case = TRUE), "UMAP",
    fifelse(grepl("composition|count", rel_path, ignore.case = TRUE), "composition/count",
      fifelse(grepl("marker_panel", rel_path, ignore.case = TRUE), "marker panels",
        fifelse(grepl("marker|dotplot|heatmap", rel_path, ignore.case = TRUE), "markers", "other"))))]
  safe_write_dt(dt, file.path(comp_dir, "core_figure_inventory.tsv"))
  plot_dt <- dt[, .(n_figures = .N, total_mb = sum(bytes, na.rm = TRUE) / 1024^2), by = figure_group][order(n_figures)]
  plot_dt[, figure_group := factor(figure_group, levels = figure_group)]
  p <- plot_count_bar(plot_dt, "figure_group", "n_figures", title = "Core tissue-comparison figure inventory", xlab = NULL, ylab = "Number of files")
  fig <- save_plot(p, file.path(comp_dir, "core_figure_inventory_barplot.png"), width = 9, height = 5.5)
  register_safe(fig, dt, "core tissue-comparison figures", "figure_inventory", sprintf("Core figure inventory: %s", basename(output_dir)), manifest_path,
    extra_context = "This companion summarizes existing UMAP, composition, marker, and ssGSEA heatmap files that previously lacked paired LLM CSV/prompt/status artifacts.")
  data.table(method = "core_figures", status = "ok", n_records = nrow(dt), n_figures = length(files))
}

make_ssgsea_companions <- function(output_dir, comp_dir, manifest_path) {
  reports <- file.path(output_dir, "reports")
  specs <- data.table(
    level = c("L2", "L2", "L3", "L3"),
    database = c("hallmark", "go_bp", "hallmark", "go_bp"),
    report_file = file.path(reports, c("ssgsea_top_pathways_hallmark.csv", "ssgsea_top_pathways_go_bp.csv", "ssgsea_l3_top_pathways_hallmark.csv", "ssgsea_l3_top_pathways_go_bp.csv")),
    heatmap_stem = file.path(output_dir, "figures", c("ssgsea_heatmap_hallmark", "ssgsea_heatmap_go_bp", "ssgsea_l3_heatmap_hallmark", "ssgsea_l3_heatmap_go_bp"))
  )
  rows <- list()
  all_top <- list()
  for (i in seq_len(nrow(specs))) {
    if (!file.exists(specs$report_file[[i]])) next
    dt <- summarize_ssgsea_file(specs$report_file[[i]], specs$level[[i]], specs$database[[i]])
    if (nrow(dt) == 0L) next
    all_top[[length(all_top) + 1L]] <- dt
    heatmap <- existing_plot_path(c(paste0(specs$heatmap_stem[[i]], ".png"), paste0(specs$heatmap_stem[[i]], ".pdf")))
    if (!is.na(heatmap) && length(heatmap)) {
      register_safe(heatmap, head(dt, 300L), sprintf("ssGSEA %s %s", specs$level[[i]], specs$database[[i]]), "ssGSEA_heatmap", sprintf("ssGSEA %s %s heatmap", specs$level[[i]], specs$database[[i]]), manifest_path,
        extra_context = "Existing ssGSEA heatmap retrofitted with top-pathway data CSV and queued LLM prompt/status.")
    }
    top_dt <- copy(dt[order(-abs_metric)][seq_len(min(.N, 40L))])
    top_dt[, plot_label := factor_unique(shorten(paste(group, pathway, sep = " | "), 100), reverse = TRUE)]
    p <- ggplot(top_dt, aes(x = metric, y = plot_label, fill = direction)) +
      geom_col(width = 0.72) + geom_vline(xintercept = 0, linewidth = 0.25, color = "grey40") +
      scale_fill_manual(values = c(positive = "#2563eb", negative = "#dc2626"), na.value = "#6b7280") +
      theme_bw(base_size = 10) + labs(title = sprintf("Top ssGSEA pathways — %s %s", specs$level[[i]], specs$database[[i]]), x = "z/score", y = NULL, fill = "direction") +
      theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 7), panel.grid.major.y = element_blank())
    fig <- save_plot(p, file.path(comp_dir, sprintf("ssgsea_%s_%s_top_pathway_barplot.png", tolower(specs$level[[i]]), specs$database[[i]])), width = 13, height = max(7, 0.22 * nrow(top_dt) + 2))
    register_safe(fig, top_dt, sprintf("ssGSEA %s %s", specs$level[[i]], specs$database[[i]]), "top_pathway_barplot", sprintf("Top ssGSEA pathways — %s %s", specs$level[[i]], specs$database[[i]]), manifest_path)
    rows[[length(rows) + 1L]] <- data.table(method = sprintf("ssGSEA_%s_%s", specs$level[[i]], specs$database[[i]]), status = "ok", n_records = nrow(dt), n_figures = 2L)
  }
  if (length(all_top) > 0L) safe_write_dt(rbindlist(all_top, fill = TRUE), file.path(comp_dir, "ssgsea_top_pathway_companion_data.tsv"))
  idx <- collect_ssgsea_llm_index(output_dir)
  if (nrow(idx) > 0L) {
    safe_write_dt(idx, file.path(comp_dir, "ssgsea_existing_llm_interpretation_index.tsv"))
    idx_plot <- idx[, .(n_interpretations = .N), by = .(llm_level = fifelse(grepl("ssgsea_llm_l3", rel_path), "L3", "L2"))]
    idx_plot[, llm_level := factor(llm_level, levels = llm_level)]
    p <- plot_count_bar(idx_plot, "llm_level", "n_interpretations", title = "Existing ssGSEA LLM interpretation RDS index", xlab = NULL, ylab = "RDS files")
    fig <- save_plot(p, file.path(comp_dir, "ssgsea_existing_llm_interpretation_index.png"), width = 7, height = 4.8)
    register_safe(fig, idx, "ssGSEA existing LLM interpretations", "llm_interpretation_index", sprintf("Existing ssGSEA LLM interpretation index: %s", basename(output_dir)), manifest_path,
      extra_context = "The upstream helper wrote interpretation.rds objects; this companion exports an index CSV and queues a method-level LLM audit.")
    rows[[length(rows) + 1L]] <- data.table(method = "ssGSEA_existing_LLM_RDS_index", status = "ok", n_records = nrow(idx), n_figures = 1L)
  }
  if (length(rows) == 0L) data.table(method = "ssGSEA", status = "missing", n_records = 0L, n_figures = 0L) else rbindlist(rows, fill = TRUE)
}

make_pseudobulk_companions <- function(output_dir, comp_dir, manifest_path) {
  de_files <- list.files(output_dir, pattern = "DESeq2_results\\.csv$", recursive = TRUE, full.names = TRUE)
  de_files <- de_files[grepl("/pseudobulk_de", de_files)]
  if (length(de_files) == 0L) return(data.table(method = "pseudobulk_DE", status = "missing", n_records = 0L, n_figures = 0L))
  summary_dt <- rbindlist(lapply(de_files, summarize_de_file, output_dir = output_dir), fill = TRUE)
  safe_write_dt(summary_dt, file.path(comp_dir, "pseudobulk_de_method_summary.tsv"))
  plot_dt <- summary_dt[order(-n_sig)][seq_len(min(.N, 45L))]
  plot_dt[, label := factor_unique(shorten(paste(level, cluster, contrast, sep = " | "), 110), reverse = TRUE)]
  p <- ggplot(plot_dt, aes(x = n_sig, y = label, fill = level)) +
    geom_col(width = 0.72) + theme_bw(base_size = 10) +
    labs(title = "Pseudobulk DE significant-gene burden", x = "Significant genes (padj < 0.05, |logFC| >= 0.25)", y = NULL, fill = "level") +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 7), panel.grid.major.y = element_blank())
  fig <- save_plot(p, file.path(comp_dir, "pseudobulk_de_significant_gene_burden.png"), width = 13.5, height = max(7, 0.23 * nrow(plot_dt) + 2))
  register_safe(fig, summary_dt, "pseudobulk DE", "significant_gene_burden", sprintf("Pseudobulk DE significant-gene burden: %s", basename(output_dir)), manifest_path,
    extra_context = "Method-level companion for DESeq2_results.csv files and existing volcano plots across L2/L3 pseudobulk comparisons.")

  volcano_dt <- summary_dt[file.exists(volcano_png) | file.exists(volcano_pdf), .(source_file, volcano_png, volcano_pdf, level, cluster, contrast, n_sig, n_up, n_down, top_up_gene, top_down_gene)]
  if (nrow(volcano_dt) > 0L) safe_write_dt(volcano_dt, file.path(comp_dir, "pseudobulk_existing_volcano_index.tsv"))

  evidence_files <- list.files(output_dir, pattern = "llm_multi_db_evidence_integrated_up_down\\.csv$", recursive = TRUE, full.names = TRUE)
  evidence_files <- evidence_files[grepl("/pseudobulk_de", evidence_files)]
  rows <- list(data.table(method = "pseudobulk_DE", status = "ok", n_records = nrow(summary_dt), n_figures = 1L))
  if (length(evidence_files) > 0L) {
    ev_rows <- lapply(evidence_files, function(path) {
      dt <- read_dt(path)
      meta <- parse_de_path(path, output_dir)
      if (nrow(dt) == 0L) return(NULL)
      dt[, `:=`(
        source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
        level = meta$level,
        cluster = meta$cluster,
        contrast = meta$contrast
      )]
      dt
    })
    ev <- rbindlist(ev_rows, fill = TRUE)
    if (nrow(ev) > 0L) {
      safe_write_dt(ev, file.path(comp_dir, "pseudobulk_integrated_llm_evidence_companion_data.tsv"))
      class_col <- first_col(ev, c("biological_signal_class", "signal_class", "class", "short_call", "discovery_flag"))
      if (is.na(class_col)) {
        ev[, biological_signal_class := "record"]
        class_col <- "biological_signal_class"
      }
      ev_plot <- ev[, .(n_records = .N), by = .(signal_class = as.character(get(class_col)))][order(n_records)]
      ev_plot[, signal_class := factor(shorten(signal_class, 90), levels = shorten(signal_class, 90))]
      p2 <- plot_count_bar(ev_plot, "signal_class", "n_records", title = "Integrated LLM evidence record classes", xlab = NULL, ylab = "Records")
      fig2 <- save_plot(p2, file.path(comp_dir, "pseudobulk_integrated_llm_evidence_classes.png"), width = 10, height = max(5, 0.25 * nrow(ev_plot) + 2))
      register_safe(fig2, ev, "pseudobulk integrated evidence", "llm_evidence_class_summary", sprintf("Integrated LLM evidence classes: %s", basename(output_dir)), manifest_path,
        extra_context = "Upstream interpret_agent outputs existed as per-comparison CSV/RDS; this companion makes them visible to the standard visual LLM batch runner.")
      rows[[length(rows) + 1L]] <- data.table(method = "pseudobulk_integrated_LLM_evidence", status = "ok", n_records = nrow(ev), n_figures = 1L)
    }
  }
  rbindlist(rows, fill = TRUE)
}

make_enrichment_companions <- function(output_dir, comp_dir, manifest_path) {
  enrich_files <- list.files(output_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  enrich_files <- enrich_files[grepl("/pseudobulk_de", enrich_files) & grepl("/enrichment_(up|down)/", enrich_files)]
  enrich_files <- enrich_files[!grepl("all_enrichment|README|manifest", enrich_files, ignore.case = TRUE)]
  if (length(enrich_files) == 0L) return(data.table(method = "pathway_enrichment", status = "missing", n_records = 0L, n_figures = 0L))
  rows <- lapply(enrich_files, summarize_enrichment_file, output_dir = output_dir)
  enrich <- rbindlist(rows[!vapply(rows, is.null, logical(1))], fill = TRUE)
  if (nrow(enrich) == 0L) return(data.table(method = "pathway_enrichment", status = "empty", n_records = 0L, n_figures = 0L))
  safe_write_dt(enrich, file.path(comp_dir, "pathway_enrichment_companion_data.tsv"))
  top_dt <- enrich[order(-neg_log10_padj)][seq_len(min(.N, 60L))]
  top_dt[, label := factor_unique(shorten(paste(level, cluster, contrast, direction, database, term, sep = " | "), 125), reverse = TRUE)]
  p <- ggplot(top_dt, aes(x = neg_log10_padj, y = label, fill = database)) +
    geom_col(width = 0.72) + theme_bw(base_size = 9.5) +
    labs(title = "Top pseudobulk pathway enrichment terms", x = "-log10(adj.P)", y = NULL, fill = "database") +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 6.5), panel.grid.major.y = element_blank())
  fig <- save_plot(p, file.path(comp_dir, "pathway_enrichment_top_terms.png"), width = 15, height = max(8, 0.22 * nrow(top_dt) + 2))
  register_safe(fig, top_dt, "pathway enrichment", "top_enrichment_terms", sprintf("Top pseudobulk pathway enrichment terms: %s", basename(output_dir)), manifest_path,
    extra_context = "Aggregates existing enrichment_up/enrichment_down CSV files across pseudobulk DE comparisons, including databases such as Hallmark, GO, KEGG, PanglaoDB, CellMarker, or custom sets when present.")

  count_dt <- enrich[, .(n_terms = .N, max_neg_log10_padj = max(neg_log10_padj, na.rm = TRUE)), by = .(level, direction, database)][order(n_terms)]
  count_dt[, label := factor_unique(paste(level, direction, database, sep = " | "))]
  p2 <- ggplot(count_dt, aes(x = n_terms, y = label, fill = database)) +
    geom_col(width = 0.72) + theme_bw(base_size = 10) +
    labs(title = "Pathway enrichment coverage by level/direction/database", x = "Terms", y = NULL, fill = "database") +
    theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())
  fig2 <- save_plot(p2, file.path(comp_dir, "pathway_enrichment_coverage.png"), width = 11, height = max(5.5, 0.25 * nrow(count_dt) + 2))
  register_safe(fig2, count_dt, "pathway enrichment", "coverage_summary", sprintf("Pathway enrichment coverage: %s", basename(output_dir)), manifest_path)
  data.table(method = "pathway_enrichment", status = "ok", n_records = nrow(enrich), n_figures = 2L)
}

make_method_coverage_dashboard <- function(method_rows, output_dir, comp_dir, manifest_path) {
  dt <- rbindlist(method_rows, fill = TRUE)
  dt[is.na(status), status := "unknown"]
  dt[, has_records := n_records > 0]
  dt[, has_figures := n_figures > 0]
  dt[, coverage_score := as.integer(status == "ok") + as.integer(has_records) + as.integer(has_figures)]
  safe_write_dt(dt, file.path(comp_dir, "method_companion_coverage_summary.tsv"))
  plot_dt <- dt[, .(status = status[1], n_records = sum(n_records, na.rm = TRUE), n_figures = sum(n_figures, na.rm = TRUE), coverage_score = max(coverage_score, na.rm = TRUE)), by = method][order(coverage_score, n_records)]
  plot_dt[, method := factor(method, levels = method)]
  p <- ggplot(plot_dt, aes(x = coverage_score, y = method, fill = status)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sprintf("records=%s | figs=%s", n_records, n_figures)), hjust = -0.02, size = 3) +
    scale_x_continuous(limits = c(0, max(3.5, max(plot_dt$coverage_score, na.rm = TRUE) + 1.4))) +
    theme_bw(base_size = 11) +
    labs(title = "Method visual/CSV/LLM companion coverage", x = "Coverage score", y = NULL, fill = "status") +
    theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())
  fig <- save_plot(p, file.path(comp_dir, "method_companion_coverage_dashboard.png"), width = 12, height = max(5.5, 0.35 * nrow(plot_dt) + 2))
  register_safe(fig, dt, "tissue-comparison companion audit", "method_coverage_dashboard", sprintf("Method companion coverage: %s", basename(output_dir)), manifest_path,
    extra_context = "Final audit showing which tissue-comparison methods now have method-level visualizations, LLM data CSV/TSV, queued prompts/status, and manifest rows.")
  dt
}

process_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    warning(sprintf("Output directory not found: %s", output_dir), call. = FALSE)
    return(data.table(output_dir = output_dir, status = "missing_output", n_manifest_rows = 0L))
  }
  comp_dir <- file.path(output_dir, "figures", "method_llm_companions_20260508")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(comp_dir, "method_llm_companion_figure_companion_manifest.tsv")
  if (file.exists(manifest_path)) unlink(manifest_path, force = TRUE)

  cat(sprintf("\n=== Companion postprocess: %s ===\n", output_dir))
  method_rows <- list(
    make_core_figure_inventory(output_dir, comp_dir, manifest_path),
    make_ssgsea_companions(output_dir, comp_dir, manifest_path),
    make_pseudobulk_companions(output_dir, comp_dir, manifest_path),
    make_enrichment_companions(output_dir, comp_dir, manifest_path)
  )
  coverage <- make_method_coverage_dashboard(method_rows, output_dir, comp_dir, manifest_path)
  manifest <- if (file.exists(manifest_path)) read_dt(manifest_path) else data.table()
  readme <- c(
    sprintf("# Tissue-comparison method LLM companions — %s", basename(output_dir)),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    "This directory retrofits methods that had raw outputs but lacked standardized visual + CSV + LLM companion artifacts.",
    "",
    "## Main files",
    "",
    "- `method_llm_companion_figure_companion_manifest.tsv`: manifest consumed by the decoupled LLM batch runner.",
    "- `method_companion_coverage_summary.tsv`: method-level coverage audit.",
    "- `*_data.csv`: paired data exported by `smcanno_viz_register_figure()` next to each registered figure.",
    "- `*_LLM_PROMPT.md`, `*_LLM_ANALYSIS.md`, `*_LLM_STATUS.json`: queued or completed LLM artifacts.",
    "",
    "## Covered method families",
    "",
    paste0("- ", coverage$method, ": ", coverage$status, " (records=", coverage$n_records, ", figures=", coverage$n_figures, ")"),
    ""
  )
  writeLines(readme, file.path(comp_dir, "README.md"), useBytes = TRUE)
  data.table(
    output_dir = output_dir,
    companion_dir = normalizePath(comp_dir, winslash = "/", mustWork = FALSE),
    manifest_path = normalizePath(manifest_path, winslash = "/", mustWork = FALSE),
    status = "ok",
    n_manifest_rows = nrow(manifest),
    n_methods = length(unique(coverage$method))
  )
}

cat("=== Tissue-comparison visual + LLM companion postprocessor (2026-05-08) ===\n")
summary <- rbindlist(lapply(args, process_output_dir), fill = TRUE)
print(summary)
summary_path <- file.path(dirname(args[[1]]), "tissue_companion_postprocess_summary_20260508.tsv")
safe_write_dt(summary, summary_path)
cat(sprintf("\nSummary written: %s\n", summary_path))
if (any(summary$status != "ok")) quit(save = "no", status = 1L)
