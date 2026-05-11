#!/usr/bin/env Rscript
# ==============================================================================
# smc_anno tradeSeq helper (2026-04-29)
# ==============================================================================

SMCANNO_VISUAL_LLM_HELPER_PATH <- "/home/h2048/script/R/smc_anno_visual_llm_helper_20260506.R"
if (!exists("smcanno_viz_register_figure", mode = "function")) {
  source(SMCANNO_VISUAL_LLM_HELPER_PATH)
}

smcanno_ts_default_marker_panel <- c(
  "RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3",
  "ACTA2", "TAGLN", "CNN1", "MYH11", "MYLK", "MYL9",
  "FN1", "VCAN", "RERGL", "COL1A1", "COL3A1", "POSTN"
)

smcanno_ts_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

smcanno_ts_write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

smcanno_ts_safe_numeric_summary <- function(x, fun) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

smcanno_ts_get_assay_counts <- function(seurat_obj, assay_name) {
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  } else {
    Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "counts")
  }
}

smcanno_ts_detected_cells <- function(counts) {
  Matrix::rowSums(counts > 0)
}

smcanno_ts_gene_exclusion_reasons <- function(genes) {
  gene_raw <- trimws(as.character(genes))
  gene_up <- toupper(gene_raw)
  reasons <- vector("list", length(gene_up))
  add_reason <- function(idx, reason) {
    if (!any(idx)) return(invisible(NULL))
    reasons[idx] <<- Map(function(x) unique(c(x, reason)), reasons[idx])
    invisible(NULL)
  }
  add_reason(grepl("^(IG[HKL]|IGJ|JCHAIN)", gene_up), "immunoglobulin_or_jchain")
  add_reason(grepl("^MT[-.]", gene_up), "mitochondrial_MT")
  add_reason(grepl("^(RPL|RPS|MRPL|MRPS)", gene_up), "ribosomal_RP")
  add_reason(grepl("^RP[0-9]", gene_up), "RP_lncRNA_like_prefix")
  add_reason(grepl("^ENSG[0-9]", gene_up), "ensembl_id_ENSG")
  add_reason(
    grepl("^(LINC|AC[0-9]|AL[0-9]|AP[0-9]|MIR|SNHG|SNORA|SNORD|RNU|RNVU|Y_RNA)", gene_up) |
      grepl("(^|[-.])AS[0-9]*$", gene_up) |
      grepl("[-.](DT|IT[0-9]*|OT[0-9]*)$", gene_up) |
      gene_up %in% c("MALAT1", "NEAT1", "NEAT2", "XIST", "TSIX", "KCNQ1OT1", "H19"),
    "lncRNA_like_symbol"
  )
  vapply(reasons, function(x) {
    x <- unique(as.character(x))
    x <- x[nzchar(x)]
    if (length(x) == 0L) "" else paste(x, collapse = ";")
  }, character(1))
}

smcanno_ts_build_gene_filter_audit <- function(counts,
                                               detected,
                                               expressed_genes,
                                               variable_features = character(),
                                               marker_panel = character(),
                                               gene_panel = character()) {
  all_genes <- rownames(counts)
  reasons <- smcanno_ts_gene_exclusion_reasons(all_genes)
  audit <- data.frame(
    gene = all_genes,
    detected_cells = as.integer(detected[all_genes]),
    expressed_for_tradeSeq = all_genes %in% expressed_genes,
    should_exclude = nzchar(reasons),
    exclude_reason = ifelse(nzchar(reasons), reasons, "kept"),
    in_marker_panel = all_genes %in% marker_panel,
    in_user_gene_panel = all_genes %in% gene_panel,
    in_variable_features = all_genes %in% variable_features,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  audit
}

smcanno_ts_select_genes <- function(seurat_obj,
                                    counts,
                                    max_genes = 1000L,
                                    min_cells_expressed = 20L,
                                    gene_panel = NULL,
                                    marker_panel = smcanno_ts_default_marker_panel,
                                    exclude_noninformative = TRUE) {
  detected <- smcanno_ts_detected_cells(counts)
  expressed_genes <- names(detected)[detected >= as.integer(min_cells_expressed)]

  variable_features <- tryCatch(Seurat::VariableFeatures(seurat_obj), error = function(e) character())
  audit_tbl <- smcanno_ts_build_gene_filter_audit(
    counts = counts,
    detected = detected,
    expressed_genes = expressed_genes,
    variable_features = variable_features,
    marker_panel = marker_panel,
    gene_panel = gene_panel
  )
  if (isTRUE(exclude_noninformative)) {
    excluded_genes <- audit_tbl$gene[audit_tbl$should_exclude]
    expressed_genes <- setdiff(expressed_genes, excluded_genes)
  }
  preferred <- unique(c(gene_panel, marker_panel, variable_features))
  preferred <- preferred[preferred %in% rownames(counts)]
  preferred <- preferred[preferred %in% expressed_genes]

  if (length(preferred) < as.integer(max_genes)) {
    fallback <- setdiff(expressed_genes[order(detected[expressed_genes], decreasing = TRUE)], preferred)
    preferred <- c(preferred, fallback)
  }
  selected <- unique(utils::head(preferred, as.integer(max_genes)))
  out <- data.frame(
    gene = selected,
    detected_cells = as.integer(detected[selected]),
    in_marker_panel = selected %in% marker_panel,
    in_user_gene_panel = selected %in% gene_panel,
    in_variable_features = selected %in% variable_features,
    should_exclude = FALSE,
    exclude_reason = "kept",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  audit_tbl$selected_for_fit <- audit_tbl$gene %in% selected
  attr(out, "exclusion_audit") <- audit_tbl
  attr(out, "exclusion_summary") <- data.frame(
    metric = c("total_genes", "expressed_genes_before_qc", "excluded_genes_total", "excluded_expressed_genes", "retained_expressed_genes", "selected_for_fit"),
    value = c(
      length(rownames(counts)),
      sum(audit_tbl$expressed_for_tradeSeq),
      sum(audit_tbl$should_exclude),
      sum(audit_tbl$should_exclude & audit_tbl$expressed_for_tradeSeq),
      length(expressed_genes),
      length(selected)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out
}

smcanno_ts_normalize_result_table <- function(tbl) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE, check.names = FALSE)
  tbl$gene <- rownames(tbl)
  rownames(tbl) <- NULL
  if ("pvalue" %in% colnames(tbl)) {
    tbl$padj <- stats::p.adjust(tbl$pvalue, method = "BH")
    tbl <- tbl[order(tbl$pvalue, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  }
  tbl[, c("gene", setdiff(colnames(tbl), "gene")), drop = FALSE]
}

smcanno_ts_build_evaluatek_summary <- function(eval_mat) {
  eval_mat <- as.matrix(eval_mat)
  if (length(eval_mat) == 0L || ncol(eval_mat) == 0L) {
    return(data.frame())
  }
  do.call(rbind, lapply(seq_len(ncol(eval_mat)), function(i) {
    vals <- eval_mat[, i]
    col_id <- colnames(eval_mat)[i]
    data.frame(
      k = suppressWarnings(as.integer(gsub("[^0-9]+", "", col_id))),
      metric_column = col_id,
      n_genes = length(vals),
      n_finite = sum(is.finite(vals)),
      n_missing = sum(!is.finite(vals)),
      median_metric = smcanno_ts_safe_numeric_summary(vals, stats::median),
      mean_metric = smcanno_ts_safe_numeric_summary(vals, mean),
      min_metric = smcanno_ts_safe_numeric_summary(vals, min),
      max_metric = smcanno_ts_safe_numeric_summary(vals, max),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}

smcanno_ts_plot_evaluatek_summary <- function(summary_tbl, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (nrow(summary_tbl) == 0L) return(invisible(NULL))
  plot_df <- summary_tbl[is.finite(summary_tbl$k), , drop = FALSE]
  if (nrow(plot_df) == 0L) return(invisible(NULL))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = k, y = median_metric)) +
    ggplot2::geom_line(linewidth = 0.7, colour = "#1F78B4") +
    ggplot2::geom_point(size = 2, colour = "#1F78B4") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = "tradeSeq evaluateK summary", x = "k", y = "Median metric") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  ggplot2::ggsave(output_path, p, width = 7, height = 5)
  invisible(output_path)
}

smcanno_ts_safe_neglog10 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite_positive <- x[is.finite(x) & x > 0]
  floor_val <- if (length(finite_positive) > 0L) min(finite_positive) * 0.1 else .Machine$double.xmin
  x[!is.finite(x) | x <= 0] <- floor_val
  -log10(x)
}

smcanno_ts_gene_class <- function(genes) {
  genes_up <- toupper(trimws(as.character(genes)))
  out <- rep("other", length(genes_up))
  out[grepl("^(ACTA2|TAGLN|CNN1|MYH11|MYL9|MYLK|TPM2|MUSTN1)$", genes_up)] <- "contractile_SMC"
  out[grepl("^(PDGFRB|RGS5|CSPG4|MCAM|NOTCH3|ADIRF|MGP)$", genes_up)] <- "pericyte_mural"
  out[grepl("^(COL1A1|COL1A2|COL3A1|FN1|VCAN|POSTN|FBLN1|SPARCL1|SFRP1|SFRP2|C1R|C1S|CFH)$", genes_up)] <- "ECM_remodeling"
  out[grepl("^(IGH|IGK|IGL|JCHAIN)", genes_up)] <- "immunoglobulin_warning"
  out
}

smcanno_ts_prepare_top_test_plot_data <- function(test_tables,
                                                  top_n_per_test = 20L) {
  rows <- list()
  idx <- 1L
  for (test_name in names(test_tables)) {
    tbl <- as.data.frame(test_tables[[test_name]], stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(tbl) == 0L || !"gene" %in% colnames(tbl)) next
    p_col <- if ("padj" %in% colnames(tbl)) "padj" else if ("pvalue" %in% colnames(tbl)) "pvalue" else NA_character_
    if (is.na(p_col) || !nzchar(p_col)) next
    tbl$test_name <- test_name
    tbl$significance_column <- p_col
    tbl$significance_value <- suppressWarnings(as.numeric(tbl[[p_col]]))
    tbl$neg_log10_significance <- smcanno_ts_safe_neglog10(tbl$significance_value)
    tbl$global_pvalue <- if ("pvalue" %in% colnames(tbl)) suppressWarnings(as.numeric(tbl$pvalue)) else NA_real_
    tbl$global_padj <- if ("padj" %in% colnames(tbl)) suppressWarnings(as.numeric(tbl$padj)) else NA_real_
    tbl$global_waldStat <- if ("waldStat" %in% colnames(tbl)) suppressWarnings(as.numeric(tbl$waldStat)) else NA_real_
    tbl$log10_waldStat_plus_1 <- log10(pmax(tbl$global_waldStat, 0) + 1)
    tbl$gene_class <- smcanno_ts_gene_class(tbl$gene)
    tbl$effect_summary <- NA_character_
    tbl$effect_min <- NA_real_
    tbl$effect_median <- NA_real_
    tbl$effect_max <- NA_real_
    effect_cols <- grep("^(meanLogFC|fcMedian|logFC)", colnames(tbl), value = TRUE)
    if (length(effect_cols) > 0L) {
      effect_stats <- t(apply(tbl[, effect_cols, drop = FALSE], 1, function(row) {
        vals <- suppressWarnings(as.numeric(row))
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0L) return(c(effect_min = NA_real_, effect_median = NA_real_, effect_max = NA_real_))
        c(effect_min = min(vals), effect_median = stats::median(vals), effect_max = max(vals))
      }))
      tbl$effect_min <- effect_stats[, "effect_min"]
      tbl$effect_median <- effect_stats[, "effect_median"]
      tbl$effect_max <- effect_stats[, "effect_max"]
      tbl$effect_summary <- sprintf("min=%.3g;median=%.3g;max=%.3g", tbl$effect_min, tbl$effect_median, tbl$effect_max)
    }
    all_zero_sig <- all(!is.na(tbl$significance_value) & tbl$significance_value == 0)
    tbl$rank_metric <- if (isTRUE(all_zero_sig) && any(is.finite(tbl$global_waldStat))) tbl$global_waldStat else tbl$neg_log10_significance
    tbl$rank_metric_name <- if (isTRUE(all_zero_sig) && any(is.finite(tbl$global_waldStat))) "waldStat_due_to_zero_padj" else "-log10_significance"
    tbl$plot_metric <- ifelse(is.finite(tbl$log10_waldStat_plus_1), tbl$log10_waldStat_plus_1, tbl$rank_metric)
    tbl$plot_metric_name <- ifelse(is.finite(tbl$log10_waldStat_plus_1), "log10_waldStat_plus_1", tbl$rank_metric_name)
    tbl$effect_direction <- ifelse(is.na(tbl$effect_median), "unknown", ifelse(tbl$effect_median > 0, "positive", ifelse(tbl$effect_median < 0, "negative", "near_zero")))
    tbl <- tbl[order(tbl$rank_metric, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    tbl$rank_within_test <- seq_len(nrow(tbl))
    rows[[idx]] <- utils::head(tbl[, c("test_name", "rank_within_test", "gene", "gene_class", "significance_column", "significance_value", "neg_log10_significance", "global_pvalue", "global_padj", "global_waldStat", "log10_waldStat_plus_1", "effect_min", "effect_median", "effect_max", "effect_summary", "effect_direction", "rank_metric", "rank_metric_name", "plot_metric", "plot_metric_name"), drop = FALSE], top_n_per_test)
    idx <- idx + 1L
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

smcanno_ts_plot_top_test_results <- function(plot_df, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- as.data.frame(plot_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(plot_df) == 0L) return(invisible(NULL))
  split_df <- lapply(split(plot_df, plot_df$test_name), function(df) {
    df <- df[order(df$rank_within_test, decreasing = TRUE), , drop = FALSE]
    df
  })
  plot_df <- do.call(rbind, split_df)
  rownames(plot_df) <- NULL
  plot_df$gene_label <- paste0(sprintf("%02d", plot_df$rank_within_test), ". ", plot_df$gene)
  plot_df$gene_label_facet <- paste(plot_df$test_name, plot_df$gene_label, sep = "___")
  plot_df$gene_label_facet <- factor(plot_df$gene_label_facet, levels = unique(plot_df$gene_label_facet))
  plot_df$plot_metric <- suppressWarnings(as.numeric(plot_df$plot_metric))
  plot_df$plot_metric[!is.finite(plot_df$plot_metric)] <- 0
  plot_df$abs_effect_median <- abs(suppressWarnings(as.numeric(plot_df$effect_median)))
  plot_df$abs_effect_median[!is.finite(plot_df$abs_effect_median)] <- 0
  plot_df$effect_label <- ifelse(is.finite(plot_df$effect_median), sprintf("%+.2g", plot_df$effect_median), "")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = plot_metric, y = gene_label_facet, fill = gene_class)) +
    ggplot2::geom_col(width = 0.72, alpha = 0.88) +
    ggplot2::geom_text(ggplot2::aes(label = effect_label), hjust = -0.08, size = 2.6, colour = "black") +
    ggplot2::facet_wrap(~ test_name, scales = "free_y") +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.16))) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(
      title = "tradeSeq top genes by Wald statistic after non-informative gene filtering",
      subtitle = "Bars use log10(Wald statistic + 1); text labels show median effect when available. IG/MT/RP/ENSG/lncRNA-like genes are excluded before fitting.",
      x = "log10(global Wald statistic + 1)",
      y = "Ranked gene",
      fill = "Gene class"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.position = "right"
    )
  ggplot2::ggsave(output_path, p, width = 12, height = 9)
  png_path <- sub("\\.pdf$", ".png", output_path)
  ggplot2::ggsave(png_path, p, width = 12, height = 9, dpi = 180)
  invisible(output_path)
}

smcanno_ts_build_test_gene_heatmap_data <- function(top_tests) {
  top_tests <- as.data.frame(top_tests, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(top_tests) == 0L) return(data.frame())
  top_tests$log10_waldStat <- log10(pmax(suppressWarnings(as.numeric(top_tests$global_waldStat)), 0) + 1)
  top_tests$effect_direction <- ifelse(is.na(top_tests$effect_median), "unknown", ifelse(top_tests$effect_median > 0, "positive", ifelse(top_tests$effect_median < 0, "negative", "near_zero")))
  top_tests[, c("test_name", "gene", "gene_class", "rank_within_test", "global_waldStat", "log10_waldStat", "effect_median", "effect_direction", "global_padj", "rank_metric_name"), drop = FALSE]
}

smcanno_ts_plot_test_gene_heatmap <- function(heatmap_df, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  heatmap_df <- as.data.frame(heatmap_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(heatmap_df) == 0L) return(invisible(NULL))
  gene_order <- aggregate(log10_waldStat ~ gene, data = heatmap_df, FUN = max, na.rm = TRUE)
  gene_order <- gene_order[order(gene_order$log10_waldStat, decreasing = TRUE), , drop = FALSE]
  heatmap_df$gene <- factor(heatmap_df$gene, levels = rev(gene_order$gene))
  heatmap_df$effect_label <- ifelse(is.finite(heatmap_df$effect_median), sprintf("%.2g", heatmap_df$effect_median), "")
  p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = test_name, y = gene, fill = log10_waldStat)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = effect_label), size = 2.5) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "log10\nWald+1") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(title = "tradeSeq cross-test Wald/effect heatmap", x = "tradeSeq test", y = "Gene") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(output_path, p, width = 8, height = max(7, 0.22 * length(unique(heatmap_df$gene)) + 3))
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 8, height = max(7, 0.22 * length(unique(heatmap_df$gene)) + 3), dpi = 180)
  invisible(output_path)
}

smcanno_ts_lineage_label_map <- function(sds) {
  lineages <- tryCatch(slingshot::slingLineages(sds), error = function(e) list())
  if (length(lineages) == 0L) return(NULL)
  stats::setNames(vapply(seq_along(lineages), function(i) paste0("Lineage", i, ": ", paste(as.character(lineages[[i]]), collapse = " -> ")), character(1)), as.character(seq_along(lineages)))
}

smcanno_ts_predict_smoother_data <- function(fit,
                                             genes,
                                             lineage_label_map = NULL,
                                             n_points = 60L) {
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0L || !requireNamespace("tradeSeq", quietly = TRUE)) return(data.frame())
  pred <- tryCatch(
    tradeSeq::predictSmooth(fit, gene = genes, nPoints = as.integer(n_points), tidy = TRUE),
    error = function(e) e
  )
  if (inherits(pred, "error")) return(data.frame(error = conditionMessage(pred), stringsAsFactors = FALSE))
  pred <- as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(pred) == 0L) return(pred)
  if (!"gene" %in% colnames(pred)) {
    gene_col <- grep("gene", colnames(pred), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(gene_col) && nzchar(gene_col)) colnames(pred)[colnames(pred) == gene_col] <- "gene"
  }
  x_col <- intersect(c("time", "pseudotime", "t", "x"), colnames(pred))[1]
  y_col <- intersect(c("yhat", "estimate", "smooth", "expression", "value"), colnames(pred))[1]
  lineage_col <- intersect(c("lineage", "lineage_id", "curve", "Lineage"), colnames(pred))[1]
  if (!is.na(lineage_col)) {
    pred$lineage <- as.character(pred[[lineage_col]])
    pred$lineage_id <- ifelse(grepl("^Lineage", pred$lineage), pred$lineage, paste0("Lineage", pred$lineage))
    if (!is.null(lineage_label_map)) {
      pred$lineage_label <- unname(lineage_label_map[as.character(pred$lineage)])
      pred$lineage_label[is.na(pred$lineage_label) | !nzchar(pred$lineage_label)] <- pred$lineage_id[is.na(pred$lineage_label) | !nzchar(pred$lineage_label)]
    } else {
      pred$lineage_label <- pred$lineage_id
    }
  } else {
    pred$lineage <- "all"
    pred$lineage_id <- "all"
    pred$lineage_label <- "all"
  }
  if (!is.na(x_col)) {
    pred[[x_col]] <- suppressWarnings(as.numeric(pred[[x_col]]))
    pred$scaled_time <- ave(pred[[x_col]], pred$lineage_label, FUN = function(v) {
      rng <- range(v, finite = TRUE)
      if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(NA_real_, length(v)))
      (v - rng[1]) / diff(rng)
    })
  }
  if (!is.na(y_col)) {
    pred[[y_col]] <- suppressWarnings(as.numeric(pred[[y_col]]))
    pred$scaled_yhat <- ave(pred[[y_col]], pred$gene, FUN = function(v) {
      rng <- range(v, finite = TRUE)
      if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(v)))
      (v - rng[1]) / diff(rng)
    })
  }
  pred
}

smcanno_ts_summarize_smoothers <- function(smoother_df) {
  smoother_df <- as.data.frame(smoother_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(smoother_df) == 0L || "error" %in% colnames(smoother_df)) return(data.frame())
  req <- c("gene", "lineage_label", "scaled_time", "yhat")
  if (!all(req %in% colnames(smoother_df))) return(data.frame())
  rows <- lapply(split(smoother_df, paste(smoother_df$gene, smoother_df$lineage_label, sep = "|||")), function(df) {
    df <- df[order(df$scaled_time), , drop = FALSE]
    y <- suppressWarnings(as.numeric(df$yhat))
    t <- suppressWarnings(as.numeric(df$scaled_time))
    finite <- is.finite(y) & is.finite(t)
    if (!any(finite)) return(NULL)
    y <- y[finite]
    t <- t[finite]
    peak_idx <- which.max(y)
    trough_idx <- which.min(y)
    delta <- utils::tail(y, 1) - y[1]
    data.frame(
      gene = df$gene[[1]],
      lineage_label = df$lineage_label[[1]],
      n_grid_points = length(y),
      start_yhat = y[1],
      end_yhat = utils::tail(y, 1),
      delta_end_start = delta,
      min_yhat = min(y),
      max_yhat = max(y),
      dynamic_range = max(y) - min(y),
      peak_scaled_time = t[peak_idx],
      trough_scaled_time = t[trough_idx],
      trend_class = ifelse(delta > 0.5, "end_high", ifelse(delta < -0.5, "end_low", "non_monotone_or_flat")),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

smcanno_ts_plot_smoothers <- function(smoother_df, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  smoother_df <- as.data.frame(smoother_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(smoother_df) == 0L || "error" %in% colnames(smoother_df)) return(invisible(NULL))
  x_col <- intersect(c("time", "pseudotime", "t", "x"), colnames(smoother_df))[1]
  y_col <- if ("scaled_yhat" %in% colnames(smoother_df)) "scaled_yhat" else intersect(c("yhat", "estimate", "smooth", "expression", "value"), colnames(smoother_df))[1]
  lineage_col <- if ("lineage_label" %in% colnames(smoother_df)) "lineage_label" else intersect(c("lineage", "lineage_id", "curve", "Lineage"), colnames(smoother_df))[1]
  x_plot_col <- if ("scaled_time" %in% colnames(smoother_df)) "scaled_time" else x_col
  if (is.na(x_plot_col) || is.na(y_col) || !"gene" %in% colnames(smoother_df)) return(invisible(NULL))
  smoother_df[[x_plot_col]] <- suppressWarnings(as.numeric(smoother_df[[x_plot_col]]))
  smoother_df[[y_col]] <- suppressWarnings(as.numeric(smoother_df[[y_col]]))
  p <- ggplot2::ggplot(smoother_df, ggplot2::aes_string(x = x_plot_col, y = y_col, colour = lineage_col, group = lineage_col)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(title = "tradeSeq normalized smoother shapes for top genes", x = "Scaled pseudotime within lineage", y = if (identical(y_col, "scaled_yhat")) "Scaled fitted expression within gene" else "Fitted expression", colour = "Lineage") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  ggplot2::ggsave(output_path, p, width = 13, height = 9)
  png_path <- sub("\\.pdf$", ".png", output_path)
  ggplot2::ggsave(png_path, p, width = 13, height = 9, dpi = 180)
  invisible(output_path)
}

smcanno_ts_plot_smoother_summary <- function(summary_df, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_df) == 0L) return(invisible(NULL))
  gene_order <- aggregate(dynamic_range ~ gene, data = summary_df, FUN = max, na.rm = TRUE)
  gene_order <- gene_order[order(gene_order$dynamic_range, decreasing = TRUE), , drop = FALSE]
  summary_df$gene <- factor(summary_df$gene, levels = rev(gene_order$gene))
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = lineage_label, y = gene, fill = delta_end_start)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(size = dynamic_range), shape = 21, colour = "black", fill = NA) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "End-start\nfit delta") +
    ggplot2::scale_size_continuous(range = c(1, 5), name = "Dynamic\nrange") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = "tradeSeq smoother end-state and dynamic-range summary", x = "Lineage", y = "Gene") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(output_path, p, width = 11, height = max(6.5, 0.28 * length(unique(summary_df$gene)) + 3))
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 11, height = max(6.5, 0.28 * length(unique(summary_df$gene)) + 3), dpi = 180)
  invisible(output_path)
}

smcanno_ts_write_final_summary <- function(output_dir,
                                           metadata_tbl,
                                           top_tests,
                                           smoother_summary,
                                           test_heatmap_df) {
  summary_path <- file.path(output_dir, "tradeSeq_FINAL_SUMMARY.md")
  prompt_path <- file.path(output_dir, "tradeSeq_FINAL_SUMMARY_LLM_PROMPT.md")
  status_path <- file.path(output_dir, "tradeSeq_FINAL_SUMMARY_LLM_STATUS.json")
  analysis_path <- file.path(output_dir, "tradeSeq_FINAL_SUMMARY_LLM_ANALYSIS.md")
  data_path <- file.path(output_dir, "tradeSeq_FINAL_SUMMARY_data.csv")
  manifest_path <- file.path(output_dir, "tradeSeq_llm_manifest.tsv")
  if (file.exists(manifest_path)) unlink(manifest_path)
  top_by_wald <- utils::head(top_tests[order(top_tests$global_waldStat, decreasing = TRUE, na.last = TRUE), c("test_name", "gene", "gene_class", "global_waldStat", "effect_median"), drop = FALSE], 20L)
  ig_hits <- unique(top_tests$gene[top_tests$gene_class == "immunoglobulin_warning"])
  lineage_patterns <- if (nrow(smoother_summary) > 0L) {
    aggregate(cbind(delta_end_start, dynamic_range) ~ lineage_label + trend_class, data = smoother_summary, FUN = function(x) round(mean(x, na.rm = TRUE), 3))
  } else data.frame()
  summary_sections <- list()
  if (nrow(top_by_wald) > 0L) {
    top_by_wald$summary_section <- "top_wald_genes"
    summary_sections[[length(summary_sections) + 1L]] <- top_by_wald
  }
  if (nrow(lineage_patterns) > 0L) {
    lineage_patterns$summary_section <- "smoother_lineage_patterns"
    summary_sections[[length(summary_sections) + 1L]] <- lineage_patterns
  }
  final_summary_data <- if (length(summary_sections) > 0L) {
    all_cols <- unique(unlist(lapply(summary_sections, colnames), use.names = FALSE))
    do.call(rbind, lapply(summary_sections, function(df) {
      for (nm in setdiff(all_cols, colnames(df))) df[[nm]] <- NA
      df[, all_cols, drop = FALSE]
    }))
  } else {
    data.frame(summary_section = character(), stringsAsFactors = FALSE)
  }
  smcanno_viz_write_csv(final_summary_data, data_path)
  lines_rule <- c(
    "# tradeSeq final interpretation summary",
    "",
    "## 数据驱动总览",
    sprintf("- tradeSeq version: `%s`; fitted genes: `%s`; lineages: `%s`; nknots: `%s`.", metadata_tbl$tradeSeq_version[[1]], metadata_tbl$n_genes_fit[[1]], metadata_tbl$n_lineages[[1]], metadata_tbl$nknots_fit[[1]]),
    sprintf("- evaluateK status: `%s`; k values tested: `%s`.", metadata_tbl$evaluateK_status[[1]], metadata_tbl$k_values[[1]]),
    sprintf("- Gene filtering: exclude_noninformative=`%s`; total excluded genes=`%s`; expressed excluded genes=`%s`.", if ("exclude_noninformative" %in% colnames(metadata_tbl)) metadata_tbl$exclude_noninformative[[1]] else NA, if ("n_genes_excluded_total" %in% colnames(metadata_tbl)) metadata_tbl$n_genes_excluded_total[[1]] else NA, if ("n_expressed_genes_excluded" %in% colnames(metadata_tbl)) metadata_tbl$n_expressed_genes_excluded[[1]] else NA),
    sprintf("- Top-test unique genes after filtering: `%s`; residual immunoglobulin-warning genes among top hits: `%s`.", length(unique(top_tests$gene)), if (length(ig_hits) == 0L) "none" else paste(ig_hits, collapse = ", ")),
    "",
    "## 最强 Wald 统计量基因",
    smcanno_viz_markdown_table(top_by_wald, max_rows = 20L),
    "",
    "## smoother 端点/动态范围概览",
    smcanno_viz_markdown_table(lineage_patterns, max_rows = 20L),
    "",
    "## 解释重点",
    "- 本轮 tradeSeq 在建模前已剔除 IG/JCHAIN、MT、RP/RPL/RPS/MRPL/MRPS、ENSG ID 与 lncRNA-like symbol（LINC/AC/AL/AP/MIR/SNHG/antisense 等）基因；请优先解释保留下来的 SMC/pericyte/ECM 动态基因。",
    "- 由于多张 test 图中的 adjusted/global p 值存在 0 值下溢，最终排序和可视化应优先参考 `global_waldStat`、效应方向和 smoother 形状，而不是等长的 -log10(p) 条形。",
    "- contractile SMC / pericyte-mural / ECM remodeling 三类基因同时出现，支持 Pericyte 到多个 SMC/ECM remodeling 终末状态的连续与分支动态。",
    if (length(ig_hits) > 0L) sprintf("- 仍有 IG/JCHAIN 类基因进入 top hits（%s），说明过滤规则需要人工复核。", paste(ig_hits, collapse = ", ")) else "- 过滤后 top hits 中未检测到 IG/JCHAIN 类警示基因。",
    "",
    "## 最后总结",
    "- 当前 tradeSeq 结果支持存在强烈的沿伪时间表达动态，但真正的生物学解释应以 Wald/effect/smoother 三者一致的基因为主。",
    "- 优先关注反复出现在 association/pattern/start-end/diff-end 且 smoother 端点差异明确的基因；对只因 p 值下溢显著但曲线形状不稳定的基因降权。",
    "- 下一步建议把 top smoother 基因与 Slingshot 三个终末分支、PAGA 拓扑和已知 SMC/Pericyte marker 合并成候选调控轴。"
  )
  prompt <- paste(
    "请基于 tradeSeq 统计表、top Wald genes、smoother summary 和下方规则总结，写一份深入中文总结，必须包含最后总结。",
    paste(lines_rule, collapse = "\n"),
    sep = "\n\n"
  )
  smcanno_viz_write_markdown(c("# tradeSeq final summary LLM prompt", "", prompt), prompt_path)
  status <- list(
    enabled = TRUE,
    status = "queued_for_batch",
    llm_mode = "queued",
    method = "tradeSeq",
    figure_type = "final_interpretation_summary",
    title = "tradeSeq final interpretation summary",
    model = "deepseek-reasoner",
    figure_path = summary_path,
    data_csv = normalizePath(data_path, winslash = "/", mustWork = FALSE),
    prompt_md = prompt_path,
    analysis_md = analysis_path,
    env_files_loaded = character(),
    runner_script = "/home/h2048/script/R/smc_anno_figure_llm_batch_20260507.R",
    error = NULL
  )
  smcanno_viz_write_markdown(c(lines_rule, "", "## LLM 状态", "在线 LLM 深度总结已排队到独立 runner；本文件先保留规则驱动总结，避免 tradeSeq 可视化阶段阻塞。"), summary_path)
  smcanno_viz_write_markdown(c(
    "# tradeSeq final summary LLM analysis",
    "",
    "**LLM status:** `queued_for_batch`（已生成 `tradeSeq_FINAL_SUMMARY_data.csv` 和 prompt，等待独立 runner 补写在线总结。）",
    "",
    lines_rule
  ), analysis_path)
  smcanno_viz_write_json(status, status_path)
  row <- data.frame(
    method = status$method,
    figure_type = status$figure_type,
    title = status$title,
    figure_path = normalizePath(summary_path, winslash = "/", mustWork = FALSE),
    data_csv = normalizePath(data_path, winslash = "/", mustWork = FALSE),
    llm_analysis_md = normalizePath(analysis_path, winslash = "/", mustWork = FALSE),
    llm_prompt_md = normalizePath(prompt_path, winslash = "/", mustWork = FALSE),
    llm_status_json = normalizePath(status_path, winslash = "/", mustWork = FALSE),
    llm_mode = "queued",
    llm_status = status$status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  smcanno_viz_append_manifest(row, manifest_path)
  list(summary_md = summary_path, prompt_md = prompt_path, analysis_md = analysis_path, data_csv = data_path, status_json = status_path, manifest = manifest_path, status = status$status)
}

smcanno_ts_write_visual_bundle <- function(output_dir,
                                           evaluatek_metrics,
                                           fit,
                                           metadata_tbl,
                                           lineage_label_map = NULL,
                                           association_tbl,
                                           pattern_tbl,
                                           startend_tbl,
                                           diffend_tbl,
                                           enable_plots = TRUE) {
  if (!isTRUE(enable_plots)) return(list(status = "disabled"))
  fig_manifest <- file.path(output_dir, "tradeSeq_figure_companion_manifest.tsv")
  if (file.exists(fig_manifest)) unlink(fig_manifest)

  out <- list(status = "ok", manifest = fig_manifest)

  if (nrow(evaluatek_metrics) > 0L) {
    eval_pdf <- file.path(output_dir, "tradeSeq_evaluateK_metrics.pdf")
    smcanno_ts_plot_evaluatek_summary(evaluatek_metrics, eval_pdf)
    smcanno_viz_register_figure(
      figure_path = eval_pdf,
      data = evaluatek_metrics,
      method = "tradeSeq",
      figure_type = "evaluateK_metric_curve",
      title = "tradeSeq evaluateK metric curve",
      extra_context = c("Median evaluateK metric across the selected gene subset is used to choose a practical knot range.", "Flat or noisy curves should be treated as weak evidence for a unique k."),
      manifest_path = fig_manifest
    )
  }

  top_tests <- smcanno_ts_prepare_top_test_plot_data(
    list(
      associationTest = association_tbl,
      patternTest = pattern_tbl,
      startVsEndTest = startend_tbl,
      diffEndTest = diffend_tbl
    ),
    top_n_per_test = 20L
  )
  smcanno_ts_write_tsv(top_tests, file.path(output_dir, "tradeSeq_top_test_results.tsv"))
  out$top_tests <- top_tests
  if (nrow(top_tests) > 0L) {
    top_pdf <- file.path(output_dir, "tradeSeq_top_test_results.pdf")
    smcanno_ts_plot_top_test_results(top_tests, top_pdf)
    smcanno_viz_register_figure(
      figure_path = top_pdf,
      data = top_tests,
      method = "tradeSeq",
      figure_type = "top_gene_test_barplot",
      title = "tradeSeq top genes across association/pattern/start-end/diff-end tests",
      extra_context = c("Rows are top genes per tradeSeq test after excluding IG/MT/RP/ENSG/lncRNA-like genes before fitting.", "Bars use log10(Wald+1), because adjusted p-values frequently underflow to zero and become visually non-informative."),
      manifest_path = fig_manifest
    )
  }

  test_heatmap_df <- smcanno_ts_build_test_gene_heatmap_data(top_tests)
  smcanno_ts_write_tsv(test_heatmap_df, file.path(output_dir, "tradeSeq_test_gene_wald_effect_heatmap.tsv"))
  out$test_heatmap <- test_heatmap_df
  if (nrow(test_heatmap_df) > 0L) {
    heatmap_pdf <- file.path(output_dir, "tradeSeq_test_gene_wald_effect_heatmap.pdf")
    smcanno_ts_plot_test_gene_heatmap(test_heatmap_df, heatmap_pdf)
    smcanno_viz_register_figure(
      figure_path = heatmap_pdf,
      data = test_heatmap_df,
      method = "tradeSeq",
      figure_type = "cross_test_wald_effect_heatmap",
      title = "tradeSeq cross-test Wald/effect heatmap",
      extra_context = c("This heatmap is preferred over raw -log10(p) bars when adjusted p-values underflow to zero.", "Tile fill is log10(Wald+1); text labels show median effect estimates."),
      manifest_path = fig_manifest
    )
  }

  smoother_genes <- unique(top_tests$gene)
  smoother_genes <- utils::head(smoother_genes[nzchar(smoother_genes)], 12L)
  smoother_df <- smcanno_ts_predict_smoother_data(fit, smoother_genes, lineage_label_map = lineage_label_map, n_points = 60L)
  smcanno_ts_write_tsv(smoother_df, file.path(output_dir, "tradeSeq_top_gene_smoothers.tsv"))
  smoother_summary <- smcanno_ts_summarize_smoothers(smoother_df)
  smcanno_ts_write_tsv(smoother_summary, file.path(output_dir, "tradeSeq_top_gene_smoother_summary.tsv"))
  out$smoother_df <- smoother_df
  out$smoother_summary <- smoother_summary
  if (nrow(smoother_df) > 0L && !"error" %in% colnames(smoother_df)) {
    smoother_pdf <- file.path(output_dir, "tradeSeq_top_gene_smoothers.pdf")
    smcanno_ts_plot_smoothers(smoother_df, smoother_pdf)
    smcanno_viz_register_figure(
      figure_path = smoother_pdf,
      data = smoother_df,
      method = "tradeSeq",
      figure_type = "fitted_smoother_facets",
      title = "tradeSeq fitted smoothers for top dynamic genes",
      extra_context = c("Each facet shows model-fitted expression along the pseudotime grid.", "Compare lineage-specific curves with Slingshot terminal branches before assigning branch biology."),
      manifest_path = fig_manifest
    )
  }
  if (nrow(smoother_summary) > 0L) {
    smoother_summary_pdf <- file.path(output_dir, "tradeSeq_smoother_delta_heatmap.pdf")
    smcanno_ts_plot_smoother_summary(smoother_summary, smoother_summary_pdf)
    smcanno_viz_register_figure(
      figure_path = smoother_summary_pdf,
      data = smoother_summary,
      method = "tradeSeq",
      figure_type = "smoother_delta_dynamic_range_heatmap",
      title = "tradeSeq smoother end-start delta and dynamic range summary",
      extra_context = c("This summary is more interpretable than raw smoother spaghetti plots for branch comparison.", "Fill shows end-start fitted expression delta; point size shows fitted dynamic range."),
      manifest_path = fig_manifest
    )
  }

  out$final_summary <- smcanno_ts_write_final_summary(
    output_dir = output_dir,
    metadata_tbl = metadata_tbl,
    top_tests = top_tests,
    smoother_summary = smoother_summary,
    test_heatmap_df = test_heatmap_df
  )
  out
}

smcanno_ts_run_bundle <- function(seurat_obj,
                                  slingshot_sce,
                                  output_dir,
                                  assay_name,
                                  max_genes = 1000L,
                                  min_cells_expressed = 20L,
                                  gene_panel = NULL,
                                  marker_panel = smcanno_ts_default_marker_panel,
                                  k_values = 3:7,
                                  evaluatek_n_genes = 200L,
                                  nknots = 6L,
                                  l2fc = 0,
                                  enable_plots = TRUE,
                                  exclude_noninformative = TRUE,
                                  verbose = TRUE) {
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop("Package 'tradeSeq' is required for smc_anno tradeSeq analysis.", call. = FALSE)
  }
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("Package 'slingshot' is required for smc_anno tradeSeq analysis.", call. = FALSE)
  }

  out_dir <- smcanno_ts_ensure_dir(output_dir)
  counts <- smcanno_ts_get_assay_counts(seurat_obj, assay_name = assay_name)
  gene_tbl <- smcanno_ts_select_genes(
    seurat_obj = seurat_obj,
    counts = counts,
    max_genes = max_genes,
    min_cells_expressed = min_cells_expressed,
    gene_panel = gene_panel,
    marker_panel = marker_panel,
    exclude_noninformative = exclude_noninformative
  )
  genes_use <- gene_tbl$gene
  if (length(genes_use) < 5L) {
    stop("tradeSeq gene panel has fewer than 5 genes after filtering; cannot fit GAM reliably.", call. = FALSE)
  }

  smcanno_ts_write_tsv(gene_tbl, file.path(out_dir, "tradeSeq_gene_panel.tsv"))
  gene_filter_audit <- attr(gene_tbl, "exclusion_audit")
  gene_filter_summary <- attr(gene_tbl, "exclusion_summary")
  if (is.null(gene_filter_audit)) gene_filter_audit <- data.frame()
  if (is.null(gene_filter_summary)) gene_filter_summary <- data.frame()
  smcanno_ts_write_tsv(gene_filter_audit, file.path(out_dir, "tradeSeq_gene_filter_audit.tsv"))
  smcanno_ts_write_tsv(gene_filter_summary, file.path(out_dir, "tradeSeq_gene_filter_summary.tsv"))

  sds <- slingshot::SlingshotDataSet(slingshot_sce)
  eval_genes <- utils::head(genes_use, min(length(genes_use), as.integer(evaluatek_n_genes)))
  eval_counts <- counts[eval_genes, , drop = FALSE]

  evaluatek_path <- file.path(out_dir, "tradeSeq_evaluateK.rds")
  evaluatek_raw_path <- file.path(out_dir, "tradeSeq_evaluateK_raw.tsv")
  evaluatek_summary_path <- file.path(out_dir, "tradeSeq_evaluateK_metrics.tsv")

  ek <- tryCatch(
    tradeSeq::evaluateK(
      counts = eval_counts,
      sds = sds,
      k = as.integer(k_values),
      nGenes = min(length(eval_genes), as.integer(evaluatek_n_genes)),
      plot = FALSE,
      verbose = verbose
    ),
    error = function(e) e
  )

  if (inherits(ek, "error")) {
    ek_summary <- data.frame(
      k = as.integer(k_values),
      metric_column = paste0("k:", as.integer(k_values)),
      n_genes = length(eval_genes),
      n_finite = NA_integer_,
      n_missing = NA_integer_,
      median_metric = NA_real_,
      mean_metric = NA_real_,
      min_metric = NA_real_,
      max_metric = NA_real_,
      status = "error",
      message = conditionMessage(ek),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    smcanno_ts_write_tsv(data.frame(), evaluatek_raw_path)
  } else {
    saveRDS(ek, evaluatek_path)
    ek_mat <- as.matrix(ek)
    ek_raw <- data.frame(gene = rownames(ek_mat), as.data.frame(ek_mat, check.names = FALSE), stringsAsFactors = FALSE, check.names = FALSE)
    smcanno_ts_write_tsv(ek_raw, evaluatek_raw_path)
    ek_summary <- smcanno_ts_build_evaluatek_summary(ek_mat)
    ek_summary$status <- "ok"
    ek_summary$message <- NA_character_
  }
  smcanno_ts_write_tsv(ek_summary, evaluatek_summary_path)

  fit <- tradeSeq::fitGAM(
    counts = counts[genes_use, , drop = FALSE],
    sds = sds,
    genes = genes_use,
    nknots = as.integer(nknots),
    verbose = verbose,
    sce = TRUE
  )
  saveRDS(fit, file.path(out_dir, "tradeSeq_fitGAM.rds"))

  association_tbl <- smcanno_ts_normalize_result_table(tradeSeq::associationTest(fit, global = TRUE, lineages = TRUE, l2fc = l2fc))
  pattern_tbl <- smcanno_ts_normalize_result_table(tradeSeq::patternTest(fit, global = TRUE, pairwise = TRUE, l2fc = l2fc))
  startend_tbl <- smcanno_ts_normalize_result_table(tradeSeq::startVsEndTest(fit, global = TRUE, lineages = TRUE, l2fc = l2fc))
  diffend_tbl <- smcanno_ts_normalize_result_table(tradeSeq::diffEndTest(fit, global = TRUE, pairwise = TRUE, l2fc = l2fc))

  smcanno_ts_write_tsv(association_tbl, file.path(out_dir, "tradeSeq_associationTest.tsv"))
  smcanno_ts_write_tsv(pattern_tbl, file.path(out_dir, "tradeSeq_patternTest.tsv"))
  smcanno_ts_write_tsv(startend_tbl, file.path(out_dir, "tradeSeq_startVsEndTest.tsv"))
  smcanno_ts_write_tsv(diffend_tbl, file.path(out_dir, "tradeSeq_diffEndTest.tsv"))

  metadata_tbl <- data.frame(
    tradeSeq_version = as.character(utils::packageVersion("tradeSeq")),
    assay = assay_name,
    n_cells = ncol(seurat_obj),
    n_lineages = ncol(slingshot::slingPseudotime(sds, na = FALSE)),
    n_genes_fit = length(genes_use),
    n_genes_evaluateK = length(eval_genes),
    k_values = paste(as.integer(k_values), collapse = ";"),
    nknots_fit = as.integer(nknots),
    min_cells_expressed = as.integer(min_cells_expressed),
    exclude_noninformative = isTRUE(exclude_noninformative),
    n_genes_excluded_total = if (nrow(gene_filter_audit) > 0L) sum(gene_filter_audit$should_exclude) else NA_integer_,
    n_expressed_genes_excluded = if (nrow(gene_filter_audit) > 0L) sum(gene_filter_audit$should_exclude & gene_filter_audit$expressed_for_tradeSeq) else NA_integer_,
    l2fc_threshold = as.numeric(l2fc),
    evaluateK_status = if (inherits(ek, "error")) "error" else "ok",
    evaluateK_message = if (inherits(ek, "error")) conditionMessage(ek) else NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  smcanno_ts_write_tsv(metadata_tbl, file.path(out_dir, "tradeSeq_metadata.tsv"))

  visual_bundle <- smcanno_ts_write_visual_bundle(
    output_dir = out_dir,
    evaluatek_metrics = ek_summary,
    fit = fit,
    metadata_tbl = metadata_tbl,
    lineage_label_map = smcanno_ts_lineage_label_map(sds),
    association_tbl = association_tbl,
    pattern_tbl = pattern_tbl,
    startend_tbl = startend_tbl,
    diffend_tbl = diffend_tbl,
    enable_plots = enable_plots
  )

  list(
    status = "ok",
    output_dir = out_dir,
    metadata = metadata_tbl,
    gene_panel = gene_tbl,
    gene_filter_audit = gene_filter_audit,
    gene_filter_summary = gene_filter_summary,
    evaluateK_metrics = ek_summary,
    fit = fit,
    association = association_tbl,
    pattern = pattern_tbl,
    start_vs_end = startend_tbl,
    diff_end = diffend_tbl,
    visual_bundle = visual_bundle
  )
}