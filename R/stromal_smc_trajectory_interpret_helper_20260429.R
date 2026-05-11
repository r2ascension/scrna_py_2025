#!/usr/bin/env Rscript
# ==============================================================================
# Stromal SMC trajectory interpretation helper (2026-04-29)
# ==============================================================================

SMCANNO_VISUAL_LLM_HELPER_PATH <- "/home/h2048/script/R/smc_anno_visual_llm_helper_20260506.R"
if (!exists("smcanno_viz_register_figure", mode = "function")) {
  source(SMCANNO_VISUAL_LLM_HELPER_PATH)
}

smcti_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

smcti_write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

smcti_safe_numeric_summary <- function(x, fun) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

smcti_safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, na.rm = TRUE, type = 8))
}

smcti_sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))

smcti_get_assay_matrix <- function(seurat_obj, assay_name, slot_name = c("counts", "data")) {
  slot_name <- match.arg(slot_name)
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = slot_name)
  } else {
    Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = slot_name)
  }
}

smcti_ensure_named_matrix <- function(x, prefix) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0(prefix, seq_len(ncol(x)))
  x
}

smcti_build_primary_lineage <- function(weight_mat) {
  apply(weight_mat, 1, function(x) {
    finite_idx <- which(is.finite(x))
    if (length(finite_idx) == 0L) return(NA_character_)
    best_idx <- finite_idx[which.max(x[finite_idx])]
    if (!is.finite(x[best_idx]) || x[best_idx] <= 0) return(NA_character_)
    colnames(weight_mat)[best_idx]
  })
}

smcti_lineage_terminal_map <- function(lineage_summary) {
  setNames(as.character(lineage_summary$terminal_state), as.character(lineage_summary$lineage_id))
}

smcti_primary_terminal_state <- function(pseudotime_table, lineage_summary) {
  terminal_map <- smcti_lineage_terminal_map(lineage_summary)
  out <- terminal_map[as.character(pseudotime_table$primary_lineage)]
  unname(ifelse(is.na(pseudotime_table$primary_lineage), NA_character_, out))
}

smcti_run_slingshot_analysis <- function(seurat_obj,
                                         assay_name,
                                         reduction_name,
                                         cluster_col,
                                         start_cluster = NULL,
                                         end_clusters = NULL,
                                         meta_cols = NULL,
                                         weight_threshold = 0.5) {
  cluster_values <- as.character(seurat_obj@meta.data[[cluster_col]])
  cluster_values[is.na(cluster_values)] <- ""
  if (any(!nzchar(trimws(cluster_values)))) {
    stop(sprintf("Metadata column '%s' contains empty or NA labels", cluster_col), call. = FALSE)
  }
  cluster_labels <- factor(cluster_values)
  if (!is.null(start_cluster) && !start_cluster %in% levels(cluster_labels)) {
    stop(sprintf("Start cluster not found for reduction '%s': %s", reduction_name, start_cluster), call. = FALSE)
  }
  if (!is.null(end_clusters)) {
    missing_end <- setdiff(as.character(end_clusters), levels(cluster_labels))
    if (length(missing_end) > 0L) {
      stop(sprintf("End clusters not found for reduction '%s': %s", reduction_name, paste(missing_end, collapse = ", ")), call. = FALSE)
    }
  }

  sce <- suppressWarnings(as.SingleCellExperiment(seurat_obj, assay = assay_name))
  trajectory_coords <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)[colnames(seurat_obj), , drop = FALSE]
  SingleCellExperiment::reducedDim(sce, "TRAJECTORY") <- trajectory_coords
  SummarizedExperiment::colData(sce)[[cluster_col]] <- cluster_labels

  sling_args <- list(
    data = sce,
    clusterLabels = cluster_col,
    reducedDim = "TRAJECTORY"
  )
  if (!is.null(start_cluster) && nzchar(start_cluster)) sling_args$start.clus <- start_cluster
  if (!is.null(end_clusters) && length(end_clusters) > 0L) sling_args$end.clus <- as.character(end_clusters)
  sce <- do.call(slingshot::slingshot, sling_args)

  pseudotime_mat <- smcti_ensure_named_matrix(slingshot::slingPseudotime(sce), prefix = "Lineage")
  curve_weight_mat <- smcti_ensure_named_matrix(slingshot::slingCurveWeights(sce), prefix = "Lineage")
  lineages <- slingshot::slingLineages(sce)

  pseudotime_table <- data.frame(
    cell_id = colnames(seurat_obj),
    seurat_obj@meta.data[colnames(seurat_obj), meta_cols, drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (lineage_name in colnames(pseudotime_mat)) {
    pseudotime_table[[paste0(lineage_name, "_pseudotime")]] <- pseudotime_mat[, lineage_name]
    pseudotime_table[[paste0(lineage_name, "_curve_weight")]] <- curve_weight_mat[, lineage_name]
  }
  pseudotime_table$primary_lineage <- smcti_build_primary_lineage(curve_weight_mat)

  lineage_summary <- do.call(rbind, lapply(seq_along(lineages), function(i) {
    lineage_name <- names(lineages)[i]
    lineage_path <- as.character(lineages[[i]])
    data.frame(
      lineage_id = lineage_name,
      root_state = lineage_path[[1]],
      terminal_state = utils::tail(lineage_path, 1),
      lineage_path = paste(lineage_path, collapse = " -> "),
      n_cells_with_pseudotime = sum(is.finite(pseudotime_mat[, lineage_name])),
      n_cells_curve_weight_ge_0_5 = sum(curve_weight_mat[, lineage_name] >= weight_threshold, na.rm = TRUE),
      median_pseudotime = smcti_safe_numeric_summary(pseudotime_mat[, lineage_name], stats::median),
      max_pseudotime = smcti_safe_numeric_summary(pseudotime_mat[, lineage_name], max),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  cluster_summary <- do.call(rbind, lapply(levels(cluster_labels), function(cluster_name) {
    idx <- which(cluster_labels == cluster_name)
    out <- data.frame(cluster_label = cluster_name, n_cells = length(idx), stringsAsFactors = FALSE, check.names = FALSE)
    for (lineage_name in colnames(pseudotime_mat)) {
      out[[paste0(lineage_name, "_median_pseudotime")]] <- smcti_safe_numeric_summary(pseudotime_mat[idx, lineage_name], stats::median)
      out[[paste0(lineage_name, "_mean_curve_weight")]] <- smcti_safe_numeric_summary(curve_weight_mat[idx, lineage_name], mean)
      out[[paste0(lineage_name, "_cells_weight_ge_", gsub("\\.", "_", sprintf("%.2f", weight_threshold)))]] <- sum(curve_weight_mat[idx, lineage_name] >= weight_threshold, na.rm = TRUE)
    }
    out
  }))

  list(
    sce = sce,
    cluster_labels = cluster_labels,
    pseudotime_mat = pseudotime_mat,
    curve_weight_mat = curve_weight_mat,
    lineages = lineages,
    pseudotime_table = pseudotime_table,
    lineage_summary = lineage_summary,
    cluster_summary = cluster_summary,
    primary_terminal_state = smcti_primary_terminal_state(pseudotime_table, lineage_summary),
    reduction_name = reduction_name,
    cluster_col = cluster_col,
    start_cluster = start_cluster,
    end_clusters = if (is.null(end_clusters)) character() else as.character(end_clusters)
  )
}

smcti_plot_slingshot_clusters <- function(sce, labels, output_path, title_text, axis_prefix = "UMAP") {
  coords <- SingleCellExperiment::reducedDim(sce, "TRAJECTORY")
  palette_vals <- setNames(grDevices::hcl.colors(length(levels(labels)), palette = "Dark 3"), levels(labels))
  grDevices::pdf(output_path, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(coords, col = palette_vals[as.character(labels)], pch = 16, asp = 1,
       xlab = paste0(axis_prefix, "_1"), ylab = paste0(axis_prefix, "_2"), main = title_text)
  lines(slingshot::SlingshotDataSet(sce), lwd = 2, col = "black")
  legend("topright", legend = levels(labels), col = palette_vals, pch = 16, cex = 0.8, bty = "n")
}

smcti_plot_slingshot_pseudotime <- function(sce, pseudotime_vec, output_path, title_text, axis_prefix = "UMAP") {
  coords <- SingleCellExperiment::reducedDim(sce, "TRAJECTORY")
  point_cols <- rep("grey85", length(pseudotime_vec))
  valid <- is.finite(pseudotime_vec)
  if (any(valid)) {
    bins <- cut(pseudotime_vec[valid], breaks = 100, include.lowest = TRUE)
    pal <- grDevices::hcl.colors(100, palette = "viridis")
    point_cols[valid] <- pal[as.integer(bins)]
  }
  grDevices::pdf(output_path, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(coords, col = point_cols, pch = 16, asp = 1,
       xlab = paste0(axis_prefix, "_1"), ylab = paste0(axis_prefix, "_2"), main = title_text)
  lines(slingshot::SlingshotDataSet(sce), lwd = 2, col = "black")
}

smcti_build_lineage_qc_table <- function(res, weight_threshold = 0.5) {
  do.call(rbind, lapply(colnames(res$pseudotime_mat), function(lineage_name) {
    pt <- res$pseudotime_mat[, lineage_name]
    wt <- res$curve_weight_mat[, lineage_name]
    data.frame(
      lineage_id = lineage_name,
      reduction = res$reduction_name,
      n_total_cells = length(pt),
      n_cells_with_pseudotime = sum(is.finite(pt)),
      frac_na_pseudotime = mean(!is.finite(pt)),
      n_cells_weight_ge_threshold = sum(wt >= weight_threshold, na.rm = TRUE),
      mean_curve_weight = smcti_safe_numeric_summary(wt, mean),
      median_curve_weight = smcti_safe_numeric_summary(wt, stats::median),
      curve_weight_q05 = smcti_safe_quantile(wt, 0.05),
      curve_weight_q95 = smcti_safe_quantile(wt, 0.95),
      median_pseudotime = smcti_safe_numeric_summary(pt, stats::median),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}

smcti_build_branch_overlap_cell_table <- function(res, ambiguity_threshold = 0.2) {
  wt <- as.matrix(res$curve_weight_mat)
  wt_ord <- t(apply(wt, 1, function(x) sort(x[is.finite(x)], decreasing = TRUE)))
  top1 <- rep(NA_real_, nrow(wt))
  top2 <- rep(NA_real_, nrow(wt))
  if (ncol(wt_ord) >= 1) top1 <- wt_ord[, 1]
  if (ncol(wt_ord) >= 2) top2 <- wt_ord[, 2]
  n_above <- rowSums(wt >= ambiguity_threshold, na.rm = TRUE)
  data.frame(
    cell_id = rownames(wt),
    n_lineages_weight_ge_threshold = n_above,
    top1_curve_weight = top1,
    top2_curve_weight = top2,
    ambiguity_margin = top1 - top2,
    is_multi_branch = n_above >= 2,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

smcti_build_branch_overlap_summary <- function(res, ambiguity_threshold = 0.2) {
  cell_tbl <- smcti_build_branch_overlap_cell_table(res, ambiguity_threshold = ambiguity_threshold)
  wt <- as.matrix(res$curve_weight_mat)
  pair_tbl <- do.call(rbind, lapply(utils::combn(colnames(wt), 2, simplify = FALSE), function(pair) {
    idx <- wt[, pair[1]] >= ambiguity_threshold & wt[, pair[2]] >= ambiguity_threshold
    data.frame(
      summary_level = "lineage_pair",
      lineage_a = pair[1],
      lineage_b = pair[2],
      n_cells = sum(idx, na.rm = TRUE),
      frac_cells = mean(idx, na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  overall <- data.frame(
    summary_level = "overall",
    lineage_a = NA_character_,
    lineage_b = NA_character_,
    n_cells = nrow(cell_tbl),
    frac_cells = NA_real_,
    n_cells_unassigned = sum(cell_tbl$n_lineages_weight_ge_threshold == 0, na.rm = TRUE),
    n_cells_unique_branch = sum(cell_tbl$n_lineages_weight_ge_threshold == 1, na.rm = TRUE),
    n_cells_multi_branch = sum(cell_tbl$is_multi_branch, na.rm = TRUE),
    frac_cells_multi_branch = mean(cell_tbl$is_multi_branch, na.rm = TRUE),
    mean_top1_curve_weight = smcti_safe_numeric_summary(cell_tbl$top1_curve_weight, mean),
    mean_top2_curve_weight = smcti_safe_numeric_summary(cell_tbl$top2_curve_weight, mean),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pair_tbl$n_cells_unassigned <- NA_real_
  pair_tbl$n_cells_unique_branch <- NA_real_
  pair_tbl$n_cells_multi_branch <- NA_real_
  pair_tbl$frac_cells_multi_branch <- NA_real_
  pair_tbl$mean_top1_curve_weight <- NA_real_
  pair_tbl$mean_top2_curve_weight <- NA_real_
  rbind(overall, pair_tbl)
}

smcti_build_lineage_context_table <- function(res, context_cols = c("sample", "tissue"), weight_threshold = 0.5) {
  keep_cols <- intersect(context_cols, colnames(res$pseudotime_table))
  if (length(keep_cols) == 0L) return(data.frame())
  out_rows <- list()
  idx <- 1L
  for (context_col in keep_cols) {
    values <- unique(as.character(res$pseudotime_table[[context_col]]))
    values <- values[nzchar(values) & !is.na(values)]
    for (context_value in values) {
      row_idx <- which(as.character(res$pseudotime_table[[context_col]]) == context_value)
      for (lineage_name in res$lineage_summary$lineage_id) {
        pt <- res$pseudotime_table[[paste0(lineage_name, "_pseudotime")]][row_idx]
        wt <- res$pseudotime_table[[paste0(lineage_name, "_curve_weight")]][row_idx]
        out_rows[[idx]] <- data.frame(
          context_type = context_col,
          context_value = context_value,
          lineage_id = lineage_name,
          n_cells = length(row_idx),
          n_cells_with_pseudotime = sum(is.finite(pt)),
          frac_na_pseudotime = mean(!is.finite(pt)),
          median_pseudotime = smcti_safe_numeric_summary(pt, stats::median),
          mean_curve_weight = smcti_safe_numeric_summary(wt, mean),
          n_cells_weight_ge_threshold = sum(wt >= weight_threshold, na.rm = TRUE),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, out_rows)
}

smcti_plot_density_table <- function(plot_df, x_col, color_col, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = x_col, color = color_col, fill = color_col)) +
    ggplot2::geom_density(alpha = 0.20, linewidth = 0.8) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = title_text, x = "Pseudotime", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  ggplot2::ggsave(output_path, p, width = 10, height = 6)
  invisible(output_path)
}

smcti_write_qc_bundle <- function(res,
                                  output_dir,
                                  fig_dir,
                                  context_cols = c("sample", "tissue"),
                                  ambiguity_threshold = 0.2,
                                  weight_threshold = 0.5) {
  qc_lineage <- smcti_build_lineage_qc_table(res, weight_threshold = weight_threshold)
  overlap_summary <- smcti_build_branch_overlap_summary(res, ambiguity_threshold = ambiguity_threshold)
  overlap_cells <- smcti_build_branch_overlap_cell_table(res, ambiguity_threshold = ambiguity_threshold)
  context_tbl <- smcti_build_lineage_context_table(res, context_cols = context_cols, weight_threshold = weight_threshold)

  smcti_write_tsv(qc_lineage, file.path(output_dir, "trajectory_qc_lineage.tsv"))
  smcti_write_tsv(overlap_summary, file.path(output_dir, "trajectory_qc_branch_overlap.tsv"))
  smcti_write_tsv(overlap_cells, file.path(output_dir, "trajectory_qc_branch_overlap_cells.tsv"))
  if (nrow(context_tbl) > 0L) {
    smcti_write_tsv(context_tbl, file.path(output_dir, "trajectory_qc_lineage_context.tsv"))
  } else {
    smcti_write_tsv(data.frame(), file.path(output_dir, "trajectory_qc_lineage_context.tsv"))
  }

  lineage_density <- do.call(rbind, lapply(res$lineage_summary$lineage_id, function(lineage_name) {
    data.frame(
      lineage_id = lineage_name,
      pseudotime = res$pseudotime_table[[paste0(lineage_name, "_pseudotime")]],
      stringsAsFactors = FALSE
    )
  }))
  lineage_density <- lineage_density[is.finite(lineage_density$pseudotime), , drop = FALSE]
  if (nrow(lineage_density) > 0L) {
    smcti_plot_density_table(lineage_density, x_col = "pseudotime", color_col = "lineage_id",
                             output_path = file.path(fig_dir, "trajectory_qc_pseudotime_density_by_lineage.pdf"),
                             title_text = "Pseudotime density by lineage")
  }

  if ("tissue" %in% colnames(res$pseudotime_table)) {
    tissue_density <- data.frame(
      tissue = as.character(res$pseudotime_table$tissue),
      primary_pseudotime = NA_real_,
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(res$pseudotime_table))) {
      lineage_name <- res$pseudotime_table$primary_lineage[i]
      if (!is.na(lineage_name) && nzchar(lineage_name)) {
        tissue_density$primary_pseudotime[i] <- res$pseudotime_table[[paste0(lineage_name, "_pseudotime")]][i]
      }
    }
    tissue_density <- tissue_density[is.finite(tissue_density$primary_pseudotime), , drop = FALSE]
    if (nrow(tissue_density) > 0L) {
      smcti_plot_density_table(tissue_density, x_col = "primary_pseudotime", color_col = "tissue",
                               output_path = file.path(fig_dir, "trajectory_qc_pseudotime_density_by_tissue.pdf"),
                               title_text = "Primary-lineage pseudotime density by tissue")
    }
  }

  list(lineage_qc = qc_lineage, overlap_summary = overlap_summary, overlap_cells = overlap_cells, context = context_tbl)
}

smcti_build_sample_lineage_fraction <- function(res, sample_col = "sample") {
  if (!sample_col %in% colnames(res$pseudotime_table)) return(data.frame())
  tbl <- table(as.character(res$pseudotime_table[[sample_col]]), as.character(res$pseudotime_table$primary_lineage), useNA = "no")
  sample_totals <- rowSums(tbl)
  out <- do.call(rbind, lapply(rownames(tbl), function(sample_id) {
    do.call(rbind, lapply(colnames(tbl), function(lineage_id) {
      data.frame(
        sample = sample_id,
        primary_lineage = lineage_id,
        n_cells = as.integer(tbl[sample_id, lineage_id]),
        sample_total_cells = as.integer(sample_totals[[sample_id]]),
        fraction_of_sample = if (sample_totals[[sample_id]] == 0) NA_real_ else as.integer(tbl[sample_id, lineage_id]) / sample_totals[[sample_id]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
  }))
  out
}

smcti_build_sample_pseudotime_summary <- function(res, sample_col = "sample") {
  if (!sample_col %in% colnames(res$pseudotime_table)) return(data.frame())
  out_rows <- list()
  idx <- 1L
  for (sample_id in unique(as.character(res$pseudotime_table[[sample_col]]))) {
    sample_idx <- which(as.character(res$pseudotime_table[[sample_col]]) == sample_id)
    for (lineage_name in res$lineage_summary$lineage_id) {
      pt <- res$pseudotime_table[[paste0(lineage_name, "_pseudotime")]][sample_idx]
      out_rows[[idx]] <- data.frame(
        sample = sample_id,
        lineage_id = lineage_name,
        n_cells = length(sample_idx),
        n_cells_with_pseudotime = sum(is.finite(pt)),
        mean_pseudotime = smcti_safe_numeric_summary(pt, mean),
        median_pseudotime = smcti_safe_numeric_summary(pt, stats::median),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, out_rows)
}

smcti_build_tissue_branch_bias <- function(res, tissue_col = "tissue") {
  if (!tissue_col %in% colnames(res$pseudotime_table)) return(data.frame())
  tbl <- table(as.character(res$pseudotime_table[[tissue_col]]), as.character(res$pseudotime_table$primary_lineage), useNA = "no")
  totals <- rowSums(tbl)
  do.call(rbind, lapply(rownames(tbl), function(tissue_id) {
    do.call(rbind, lapply(colnames(tbl), function(lineage_id) {
      data.frame(
        tissue = tissue_id,
        primary_lineage = lineage_id,
        n_cells = as.integer(tbl[tissue_id, lineage_id]),
        tissue_total_cells = as.integer(totals[[tissue_id]]),
        fraction_of_tissue = if (totals[[tissue_id]] == 0) NA_real_ else as.integer(tbl[tissue_id, lineage_id]) / totals[[tissue_id]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
  }))
}

smcti_plot_lineage_fraction_bar <- function(df, group_col, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(df) == 0L) return(invisible(NULL))
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = group_col, y = "fraction_of_sample", fill = "primary_lineage")) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = title_text, x = group_col, y = "Fraction") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(output_path, p, width = 9, height = 6)
  invisible(output_path)
}

smcti_write_association_bundle <- function(res, output_dir, fig_dir, sample_col = "sample", tissue_col = "tissue") {
  assoc_dir <- smcti_ensure_dir(output_dir)
  sample_lineage_fraction <- smcti_build_sample_lineage_fraction(res, sample_col = sample_col)
  sample_pseudotime_summary <- smcti_build_sample_pseudotime_summary(res, sample_col = sample_col)
  tissue_branch_bias <- smcti_build_tissue_branch_bias(res, tissue_col = tissue_col)
  smcti_write_tsv(sample_lineage_fraction, file.path(assoc_dir, "sample_lineage_fraction.tsv"))
  smcti_write_tsv(sample_pseudotime_summary, file.path(assoc_dir, "sample_pseudotime_summary.tsv"))
  smcti_write_tsv(tissue_branch_bias, file.path(assoc_dir, "tissue_branch_bias.tsv"))
  if (nrow(sample_lineage_fraction) > 0L) {
    smcti_plot_lineage_fraction_bar(sample_lineage_fraction, group_col = "sample", output_path = file.path(fig_dir, "sample_lineage_fraction_barplot.pdf"), title_text = "Sample-level branch composition")
  }
  list(sample_lineage_fraction = sample_lineage_fraction, sample_pseudotime_summary = sample_pseudotime_summary, tissue_branch_bias = tissue_branch_bias)
}

smcti_compare_terminal_assignments <- function(reference_states, other_states) {
  keep <- !is.na(reference_states) & !is.na(other_states) & nzchar(reference_states) & nzchar(other_states)
  if (!any(keep)) return(NA_real_)
  mean(reference_states[keep] == other_states[keep])
}

smcti_compare_slingshot_results <- function(reference_res,
                                            comparison_res,
                                            scenario_id,
                                            scenario_type,
                                            reduction,
                                            start_cluster,
                                            end_clusters,
                                            status = "ok",
                                            notes = NA_character_) {
  ref_states <- smcti_primary_terminal_state(reference_res$pseudotime_table, reference_res$lineage_summary)
  comp_states <- smcti_primary_terminal_state(comparison_res$pseudotime_table, comparison_res$lineage_summary)
  agreement <- smcti_compare_terminal_assignments(ref_states, comp_states)
  data.frame(
    scenario_id = scenario_id,
    scenario_type = scenario_type,
    reduction = reduction,
    start_cluster = ifelse(is.null(start_cluster), NA_character_, start_cluster),
    end_clusters = ifelse(length(end_clusters) == 0L, "AUTO", paste(end_clusters, collapse = ";")),
    status = status,
    n_lineages = nrow(comparison_res$lineage_summary),
    terminal_states = paste(sort(unique(comparison_res$lineage_summary$terminal_state)), collapse = ";"),
    finite_pseudotime_cells = sum(apply(comparison_res$pseudotime_mat, 1, function(x) any(is.finite(x)))),
    topology_match_reference = identical(sort(unique(reference_res$lineage_summary$terminal_state)), sort(unique(comparison_res$lineage_summary$terminal_state))),
    primary_terminal_agreement = agreement,
    notes = notes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

smcti_write_sensitivity_bundle <- function(seurat_obj,
                                           reference_res,
                                           output_dir,
                                           assay_name,
                                           cluster_col,
                                           start_cluster,
                                           end_clusters,
                                           reductions = c("umap", "pca", "harmony"),
                                           meta_cols = NULL,
                                           weight_threshold = 0.5) {
  sens_dir <- smcti_ensure_dir(output_dir)
  detail_dir <- smcti_ensure_dir(file.path(sens_dir, "details"))
  reductions <- unique(reductions[reductions %in% names(seurat_obj@reductions)])
  rows <- list()
  idx <- 1L
  for (red in reductions) {
    if (identical(red, reference_res$reduction_name)) next
    run_res <- tryCatch(
      smcti_run_slingshot_analysis(seurat_obj, assay_name = assay_name, reduction_name = red, cluster_col = cluster_col,
                                   start_cluster = start_cluster, end_clusters = end_clusters, meta_cols = meta_cols,
                                   weight_threshold = weight_threshold),
      error = function(e) e
    )
    if (inherits(run_res, "error")) {
      rows[[idx]] <- data.frame(scenario_id = paste0("reduction_", red), scenario_type = "reduction", reduction = red,
                                start_cluster = start_cluster, end_clusters = paste(end_clusters, collapse = ";"),
                                status = "error", n_lineages = NA_integer_, terminal_states = NA_character_,
                                finite_pseudotime_cells = NA_integer_, topology_match_reference = NA,
                                primary_terminal_agreement = NA_real_, notes = conditionMessage(run_res),
                                stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      smcti_write_tsv(run_res$lineage_summary, file.path(detail_dir, sprintf("%s_lineage_summary.tsv", smcti_sanitize_name(red))))
      rows[[idx]] <- smcti_compare_slingshot_results(reference_res, run_res,
                                                     scenario_id = paste0("reduction_", red), scenario_type = "reduction",
                                                     reduction = red, start_cluster = start_cluster, end_clusters = end_clusters)
    }
    idx <- idx + 1L
  }

  auto_res <- tryCatch(
    smcti_run_slingshot_analysis(seurat_obj, assay_name = assay_name, reduction_name = reference_res$reduction_name,
                                 cluster_col = cluster_col, start_cluster = start_cluster, end_clusters = NULL,
                                 meta_cols = meta_cols, weight_threshold = weight_threshold),
    error = function(e) e
  )
  if (inherits(auto_res, "error")) {
    rows[[idx]] <- data.frame(scenario_id = "auto_terminal", scenario_type = "terminal", reduction = reference_res$reduction_name,
                              start_cluster = start_cluster, end_clusters = "AUTO", status = "error", n_lineages = NA_integer_,
                              terminal_states = NA_character_, finite_pseudotime_cells = NA_integer_, topology_match_reference = NA,
                              primary_terminal_agreement = NA_real_, notes = conditionMessage(auto_res),
                              stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    smcti_write_tsv(auto_res$lineage_summary, file.path(detail_dir, "auto_terminal_lineage_summary.tsv"))
    rows[[idx]] <- smcti_compare_slingshot_results(reference_res, auto_res,
                                                   scenario_id = "auto_terminal", scenario_type = "terminal",
                                                   reduction = reference_res$reduction_name, start_cluster = start_cluster,
                                                   end_clusters = character())
  }
  idx <- idx + 1L

  alt_start_candidates <- grep("pericyte", levels(reference_res$cluster_labels), value = TRUE, ignore.case = TRUE)
  alt_start_candidates <- setdiff(alt_start_candidates, start_cluster)
  if (length(alt_start_candidates) > 0L) {
    alt_start <- alt_start_candidates[[1]]
    alt_res <- tryCatch(
      smcti_run_slingshot_analysis(seurat_obj, assay_name = assay_name, reduction_name = reference_res$reduction_name,
                                   cluster_col = cluster_col, start_cluster = alt_start, end_clusters = end_clusters,
                                   meta_cols = meta_cols, weight_threshold = weight_threshold),
      error = function(e) e
    )
    if (inherits(alt_res, "error")) {
      rows[[idx]] <- data.frame(scenario_id = paste0("alt_start_", smcti_sanitize_name(alt_start)), scenario_type = "start_cluster",
                                reduction = reference_res$reduction_name, start_cluster = alt_start, end_clusters = paste(end_clusters, collapse = ";"),
                                status = "error", n_lineages = NA_integer_, terminal_states = NA_character_, finite_pseudotime_cells = NA_integer_,
                                topology_match_reference = NA, primary_terminal_agreement = NA_real_, notes = conditionMessage(alt_res),
                                stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      smcti_write_tsv(alt_res$lineage_summary, file.path(detail_dir, sprintf("alt_start_%s_lineage_summary.tsv", smcti_sanitize_name(alt_start))))
      rows[[idx]] <- smcti_compare_slingshot_results(reference_res, alt_res,
                                                     scenario_id = paste0("alt_start_", smcti_sanitize_name(alt_start)),
                                                     scenario_type = "start_cluster", reduction = reference_res$reduction_name,
                                                     start_cluster = alt_start, end_clusters = end_clusters)
    }
  }

  summary_tbl <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  smcti_write_tsv(summary_tbl, file.path(sens_dir, "sensitivity_summary.tsv"))
  summary_tbl
}

smcti_get_monocle3_root_node <- function(cds, cluster_col, start_cluster) {
  cell_ids <- rownames(SummarizedExperiment::colData(cds))[as.character(SummarizedExperiment::colData(cds)[[cluster_col]]) == start_cluster]
  if (length(cell_ids) == 0L) stop(sprintf("No cells found for root cluster '%s'", start_cluster), call. = FALSE)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex)
  closest_vertex <- closest_vertex[cell_ids, , drop = FALSE]
  vertex_counts <- table(closest_vertex[, 1])
  if (length(vertex_counts) == 0L) stop("No principal graph vertices available for root selection", call. = FALSE)
  graph_nodes <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name
  graph_nodes[as.numeric(names(which.max(vertex_counts)))]
}

smcti_run_monocle3_bundle <- function(seurat_obj,
                                      output_dir,
                                      assay_name,
                                      cluster_col,
                                      start_cluster,
                                      end_clusters = NULL,
                                      reduction_name = "umap",
                                      num_dim = 30L) {
  out_dir <- smcti_ensure_dir(output_dir)
  figure_manifest <- file.path(out_dir, "monocle3_figure_companion_manifest.tsv")
  if (file.exists(figure_manifest)) unlink(figure_manifest)
  res <- tryCatch({
    counts <- smcti_get_assay_matrix(seurat_obj, assay_name, slot_name = "counts")
    cell_metadata <- seurat_obj@meta.data[colnames(seurat_obj), , drop = FALSE]
    gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE)
    cds <- monocle3::new_cell_data_set(expression_data = counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
    preprocess_dim <- max(5L, min(as.integer(num_dim), ncol(seurat_obj) - 1L, 30L))
    cds <- monocle3::preprocess_cds(cds, num_dim = preprocess_dim)
    if (reduction_name %in% names(seurat_obj@reductions)) {
      umap_coords <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)[colnames(seurat_obj), , drop = FALSE]
      if (ncol(umap_coords) >= 2L) {
        colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
        SingleCellExperiment::reducedDims(cds)$UMAP <- umap_coords[, seq_len(min(2L, ncol(umap_coords))), drop = FALSE]
      }
    }
    if (!"UMAP" %in% names(SingleCellExperiment::reducedDims(cds))) {
      cds <- monocle3::reduce_dimension(cds, preprocess_method = "PCA")
    }
    cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
    cds <- monocle3::learn_graph(cds, use_partition = FALSE)
    root_pr_node <- smcti_get_monocle3_root_node(cds, cluster_col = cluster_col, start_cluster = start_cluster)
    cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_node)

    pt <- monocle3::pseudotime(cds)
    pt_tbl <- data.frame(
      cell_id = names(pt),
      cluster_label = as.character(SummarizedExperiment::colData(cds)[names(pt), cluster_col]),
      pseudotime = as.numeric(pt),
      partition = as.character(monocle3::partitions(cds))[match(names(pt), names(monocle3::partitions(cds)))],
      monocle3_cluster = as.character(monocle3::clusters(cds))[match(names(pt), names(monocle3::clusters(cds)))],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    for (extra_col in intersect(c("sample", "tissue", "condition"), colnames(seurat_obj@meta.data))) {
      pt_tbl[[extra_col]] <- seurat_obj@meta.data[pt_tbl$cell_id, extra_col]
    }

    cluster_summary <- do.call(rbind, lapply(unique(pt_tbl$cluster_label), function(cluster_name) {
      idx <- which(pt_tbl$cluster_label == cluster_name)
      data.frame(
        cluster_label = cluster_name,
        n_cells = length(idx),
        n_cells_with_pseudotime = sum(is.finite(pt_tbl$pseudotime[idx])),
        median_pseudotime = smcti_safe_numeric_summary(pt_tbl$pseudotime[idx], stats::median),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
    partition_summary <- as.data.frame.matrix(table(pt_tbl$partition, pt_tbl$cluster_label))
    partition_summary$partition <- rownames(partition_summary)
    rownames(partition_summary) <- NULL

    smcti_write_tsv(pt_tbl, file.path(out_dir, "monocle3_pseudotime.tsv"))
    smcti_write_tsv(cluster_summary, file.path(out_dir, "monocle3_cluster_summary.tsv"))
    smcti_write_tsv(partition_summary, file.path(out_dir, "monocle3_partition_summary.tsv"))
    saveRDS(cds, file.path(out_dir, "monocle3_cds.rds"))

    umap_coords <- as.data.frame(SingleCellExperiment::reducedDims(cds)$UMAP, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
    umap_coords$cell_id <- rownames(umap_coords)
    m3_plot_df <- merge(umap_coords, pt_tbl, by = "cell_id", all.x = TRUE, sort = FALSE)
    smcti_write_tsv(m3_plot_df, file.path(out_dir, "monocle3_plot_cells_data.tsv"))

    p_celltype <- monocle3::plot_cells(cds, color_cells_by = cluster_col, label_cell_groups = FALSE,
                                       label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3)
    m3_celltype_pdf <- file.path(out_dir, "monocle3_trajectory_celltypes.pdf")
    ggplot2::ggsave(m3_celltype_pdf, p_celltype, width = 10, height = 8)
    ggplot2::ggsave(sub("\\.pdf$", ".png", m3_celltype_pdf), p_celltype, width = 10, height = 8, dpi = 180)
    smcanno_viz_register_figure(
      figure_path = m3_celltype_pdf,
      data = m3_plot_df,
      method = "Monocle3",
      figure_type = "trajectory_celltype_embedding",
      title = "Monocle3 trajectory colored by cell type",
      extra_context = c("The paired CSV includes UMAP coordinates, Monocle3 pseudotime, partitions, and clusters.", sprintf("Root cluster requested: %s", start_cluster)),
      manifest_path = figure_manifest
    )
    p_pt <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                                 label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 3)
    m3_pt_pdf <- file.path(out_dir, "monocle3_trajectory_pseudotime.pdf")
    ggplot2::ggsave(m3_pt_pdf, p_pt, width = 10, height = 8)
    ggplot2::ggsave(sub("\\.pdf$", ".png", m3_pt_pdf), p_pt, width = 10, height = 8, dpi = 180)
    smcanno_viz_register_figure(
      figure_path = m3_pt_pdf,
      data = m3_plot_df,
      method = "Monocle3",
      figure_type = "trajectory_pseudotime_embedding",
      title = "Monocle3 trajectory colored by pseudotime",
      extra_context = c("Use the paired CSV to compare Monocle3 pseudotime by cluster rather than relying only on the color gradient.", sprintf("Root graph node selected from `%s`: %s", start_cluster, root_pr_node)),
      manifest_path = figure_manifest
    )

    packet <- pa_build_monocle3_trajectory_packet(
      pseudotime_table = pt_tbl,
      lineage_summary = cluster_summary,
      branch_summary = cluster_summary,
      partition_summary = partition_summary,
      source_path = out_dir,
      root_state = start_cluster,
      terminal_states = if (is.null(end_clusters)) character() else as.character(end_clusters)
    )
    saveRDS(packet, file.path(out_dir, "monocle3_trajectory_packet.rds"))
    list(status = "ok", pseudotime_table = pt_tbl, cluster_summary = cluster_summary, partition_summary = partition_summary,
         packet = packet, root_pr_node = root_pr_node, output_dir = out_dir)
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(msg, con = file.path(out_dir, "monocle3_error.txt"), useBytes = TRUE)
    list(status = "error", error = msg, output_dir = out_dir)
  })
  res
}

smcti_patch_igraph_for_monocle2 <- function() {
  if (!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("monocle", quietly = TRUE)) return(invisible(FALSE))
  patch_targets <- c(
    "assign_cell_lineage",
    "buildBranchCellDataSet",
    "count_leaf_descendents",
    "cth_classifier_cds",
    "cth_classifier_cell",
    "extract_ddrtree_ordering",
    "extract_good_branched_ordering",
    "extract_good_ordering",
    "make_canonical",
    "measure_diameter_path",
    "project2MST"
  )
  ns <- asNamespace("monocle")
  for (fn_name in patch_targets) {
    if (!exists(fn_name, envir = ns, inherits = FALSE)) next
    fn <- get(fn_name, envir = ns, inherits = FALSE)
    if (!is.function(fn)) next
    body_text <- paste(deparse(body(fn)), collapse = "\n")
    body_text <- gsub(
      "V\\(([^)]+)\\)\\[suppressWarnings\\(nei\\(([^,\\)]+),\\s*mode = ([^)]+)\\)\\)\\]",
      "suppressWarnings(igraph::neighbors(\\1, \\2, mode = \\3))",
      body_text,
      perl = TRUE
    )
    body_text <- gsub(
      "V\\(([^)]+)\\)\\[suppressWarnings\\(nei\\(([^)]+)\\)\\)\\]",
      "suppressWarnings(igraph::neighbors(\\1, \\2, mode = 'all'))",
      body_text,
      perl = TRUE
    )
    body_text <- gsub("neimode\\s*=", "mode =", body_text, perl = TRUE)
    body(fn) <- parse(text = body_text)[[1]]
    if (bindingIsLocked(fn_name, ns)) unlockBinding(fn_name, ns)
    assign(fn_name, fn, envir = ns)
    lockBinding(fn_name, ns)
  }
  invisible(TRUE)
}

smcti_patch_dplyr_for_monocle2 <- function() {
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("rlang", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("dplyr")
  normalize_underscore_dots <- function(..., .dots = list()) {
    dots <- c(list(...), if (is.null(.dots)) list() else as.list(.dots))
    dots <- dots[!vapply(dots, is.null, logical(1))]
    unlist(lapply(dots, function(x) {
      if (inherits(x, "formula")) {
        return(list(rlang::f_rhs(x)))
      }
      if (inherits(x, "quosure") || is.symbol(x) || is.call(x)) {
        return(list(x))
      }
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      as.list(x)
    }), recursive = FALSE)
  }
  compat_select_ <- function(.data, ..., .dots = list()) {
    dots <- normalize_underscore_dots(..., .dots = .dots)
    if (length(dots) == 0L) return(.data)
    exprs <- lapply(dots, function(x) if (is.character(x)) rlang::sym(x) else x)
    dplyr::select(.data, !!!exprs)
  }
  compat_group_by_ <- function(.data, ..., .dots = list(), add = FALSE) {
    dots <- normalize_underscore_dots(..., .dots = .dots)
    if (length(dots) == 0L) return(dplyr::group_by(.data, .add = isTRUE(add)))
    exprs <- lapply(dots, function(x) if (is.character(x)) rlang::sym(x) else x)
    dplyr::group_by(.data, !!!exprs, .add = isTRUE(add))
  }
  for (fn_name in c("select_", "group_by_")) {
    if (!exists(fn_name, envir = ns, inherits = FALSE)) next
    if (bindingIsLocked(fn_name, ns)) unlockBinding(fn_name, ns)
    assign(fn_name, if (identical(fn_name, "select_")) compat_select_ else compat_group_by_, envir = ns)
    lockBinding(fn_name, ns)
  }
  if (requireNamespace("monocle", quietly = TRUE)) {
    monocle_ns <- asNamespace("monocle")
    for (fn_name in c("plot_cell_trajectory")) {
      if (!exists(fn_name, envir = monocle_ns, inherits = FALSE)) next
      fn <- get(fn_name, envir = monocle_ns, inherits = FALSE)
      if (!is.function(fn)) next
      body_text <- paste(deparse(body(fn)), collapse = "\n")
      body_text <- gsub("(?<!::)\\bselect_\\s*\\(", "dplyr::select(", body_text, perl = TRUE)
      body(fn) <- parse(text = body_text)[[1]]
      if (bindingIsLocked(fn_name, monocle_ns)) unlockBinding(fn_name, monocle_ns)
      assign(fn_name, fn, envir = monocle_ns)
      lockBinding(fn_name, monocle_ns)
    }
  }
  invisible(TRUE)
}

smcti_monocle2_ordering_genes <- function(seurat_obj, max_genes = 2000L, marker_panel = NULL) {
  ordering <- tryCatch(Seurat::VariableFeatures(seurat_obj), error = function(e) character())
  ordering <- unique(c(ordering, marker_panel))
  ordering <- ordering[ordering %in% rownames(seurat_obj)]
  if (length(ordering) == 0L) ordering <- rownames(seurat_obj)
  unique(utils::head(ordering, max_genes))
}

smcti_run_monocle2_bundle <- function(seurat_obj,
                                      output_dir,
                                      assay_name,
                                      cluster_col,
                                      start_cluster,
                                      end_clusters = NULL,
                                      max_ordering_genes = 2000L,
                                      marker_panel = NULL) {
  out_dir <- smcti_ensure_dir(output_dir)
  figure_manifest <- file.path(out_dir, "monocle2_figure_companion_manifest.tsv")
  if (file.exists(figure_manifest)) unlink(figure_manifest)
  res <- tryCatch({
    if (!requireNamespace("Biobase", quietly = TRUE)) {
      stop("Package 'Biobase' is required for Monocle2 cross-check", call. = FALSE)
    }
    if (!requireNamespace("DDRTree", quietly = TRUE)) {
      stop("Package 'DDRTree' is required for Monocle2 cross-check", call. = FALSE)
    }
    suppressPackageStartupMessages(require(DDRTree))
    smcti_patch_igraph_for_monocle2()
    smcti_patch_dplyr_for_monocle2()
    counts <- smcti_get_assay_matrix(seurat_obj, assay_name, slot_name = "counts")
    ordering_genes <- smcti_monocle2_ordering_genes(seurat_obj, max_genes = max_ordering_genes, marker_panel = marker_panel)
    counts_use <- counts[ordering_genes, , drop = FALSE]
    pheno_data <- Biobase::AnnotatedDataFrame(seurat_obj@meta.data[colnames(seurat_obj), , drop = FALSE])
    feature_data <- Biobase::AnnotatedDataFrame(data.frame(gene_short_name = rownames(counts_use), row.names = rownames(counts_use), stringsAsFactors = FALSE))
    cds <- monocle::newCellDataSet(counts_use, phenoData = pheno_data, featureData = feature_data, expressionFamily = VGAM::negbinomial.size())
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- suppressWarnings(BiocGenerics::estimateDispersions(cds))
    cds <- monocle::setOrderingFilter(cds, ordering_genes)
    cds <- monocle::reduceDimension(cds, max_components = 2, method = "DDRTree", verbose = FALSE)
    cds <- monocle::orderCells(cds)
    pd <- Biobase::pData(cds)
    root_state <- names(sort(table(pd$State[pd[[cluster_col]] == start_cluster]), decreasing = TRUE))[1]
    if (!is.na(root_state) && nzchar(root_state)) cds <- monocle::orderCells(cds, root_state = root_state)
    pd <- Biobase::pData(cds)
    pt_tbl <- data.frame(
      cell_id = rownames(pd),
      cluster_label = as.character(pd[[cluster_col]]),
      pseudotime = as.numeric(pd$Pseudotime),
      state = as.character(pd$State),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    for (extra_col in intersect(c("sample", "tissue", "condition"), colnames(seurat_obj@meta.data))) {
      pt_tbl[[extra_col]] <- seurat_obj@meta.data[pt_tbl$cell_id, extra_col]
    }
    cluster_summary <- do.call(rbind, lapply(unique(pt_tbl$cluster_label), function(cluster_name) {
      idx <- which(pt_tbl$cluster_label == cluster_name)
      state_tab <- sort(table(pt_tbl$state[idx]), decreasing = TRUE)
      data.frame(
        cluster_label = cluster_name,
        n_cells = length(idx),
        n_cells_with_pseudotime = sum(is.finite(pt_tbl$pseudotime[idx])),
        median_pseudotime = smcti_safe_numeric_summary(pt_tbl$pseudotime[idx], stats::median),
        dominant_state = if (length(state_tab) == 0L) NA_character_ else names(state_tab)[1],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
    state_summary <- as.data.frame.matrix(table(pt_tbl$state, pt_tbl$cluster_label))
    state_summary$state <- rownames(state_summary)
    rownames(state_summary) <- NULL
    smcti_write_tsv(pt_tbl, file.path(out_dir, "monocle2_pseudotime.tsv"))
    smcti_write_tsv(cluster_summary, file.path(out_dir, "monocle2_cluster_summary.tsv"))
    smcti_write_tsv(state_summary, file.path(out_dir, "monocle2_state_summary.tsv"))
    saveRDS(cds, file.path(out_dir, "monocle2_cds.rds"))

    m2_coords <- tryCatch(t(monocle::reducedDimS(cds)), error = function(e) NULL)
    if (!is.null(m2_coords) && nrow(m2_coords) == nrow(pt_tbl)) {
      m2_coords <- as.data.frame(m2_coords, stringsAsFactors = FALSE, check.names = FALSE)
      colnames(m2_coords) <- paste0("DDRTree_", seq_len(ncol(m2_coords)))
      m2_coords$cell_id <- rownames(pd)
      m2_plot_df <- merge(m2_coords, pt_tbl, by = "cell_id", all.x = TRUE, sort = FALSE)
    } else {
      m2_plot_df <- pt_tbl
    }
    smcti_write_tsv(m2_plot_df, file.path(out_dir, "monocle2_plot_cells_data.tsv"))

    m2_celltype_pdf <- file.path(out_dir, "monocle2_trajectory_celltypes.pdf")
    grDevices::pdf(m2_celltype_pdf, width = 10, height = 8)
    print(monocle::plot_cell_trajectory(cds, color_by = cluster_col))
    grDevices::dev.off()
    grDevices::png(sub("\\.pdf$", ".png", m2_celltype_pdf), width = 1800, height = 1440, res = 180)
    print(monocle::plot_cell_trajectory(cds, color_by = cluster_col))
    grDevices::dev.off()
    smcanno_viz_register_figure(
      figure_path = m2_celltype_pdf,
      data = m2_plot_df,
      method = "Monocle2",
      figure_type = "trajectory_celltype_embedding",
      title = "Monocle2 DDRTree trajectory colored by cell type",
      extra_context = c("The paired CSV includes DDRTree coordinates when available, Monocle2 pseudotime, state, and cell type.", sprintf("Root state inferred from `%s`: %s", start_cluster, root_state)),
      manifest_path = figure_manifest
    )
    m2_pt_pdf <- file.path(out_dir, "monocle2_trajectory_pseudotime.pdf")
    grDevices::pdf(m2_pt_pdf, width = 10, height = 8)
    print(monocle::plot_cell_trajectory(cds, color_by = "Pseudotime"))
    grDevices::dev.off()
    grDevices::png(sub("\\.pdf$", ".png", m2_pt_pdf), width = 1800, height = 1440, res = 180)
    print(monocle::plot_cell_trajectory(cds, color_by = "Pseudotime"))
    grDevices::dev.off()
    smcanno_viz_register_figure(
      figure_path = m2_pt_pdf,
      data = m2_plot_df,
      method = "Monocle2",
      figure_type = "trajectory_pseudotime_embedding",
      title = "Monocle2 DDRTree trajectory colored by pseudotime",
      extra_context = c("Use the paired CSV to compare Monocle2 state/pseudotime against Slingshot and Monocle3.", sprintf("Ordering genes used: %s", length(ordering_genes))),
      manifest_path = figure_manifest
    )

    list(status = "ok", pseudotime_table = pt_tbl, cluster_summary = cluster_summary, state_summary = state_summary,
         root_state = root_state, output_dir = out_dir)
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(msg, con = file.path(out_dir, "monocle2_error.txt"), useBytes = TRUE)
    list(status = "error", error = msg, output_dir = out_dir)
  })
  res
}

smcti_write_minimal_h5ad_for_paga <- function(seurat_obj,
                                              output_dir,
                                              groupby_key,
                                              reduction_name,
                                              assay_name,
                                              python_bin = "/home/h2048/miniconda3/bin/python") {
  out_dir <- smcti_ensure_dir(output_dir)
  obs_cols <- unique(c(groupby_key, intersect(c("sample", "tissue", "condition"), colnames(seurat_obj@meta.data))))
  obs <- data.frame(cell_id = colnames(seurat_obj), seurat_obj@meta.data[colnames(seurat_obj), obs_cols, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in colnames(obs)) {
    if (is.factor(obs[[nm]])) obs[[nm]] <- as.character(obs[[nm]])
  }
  obs_path <- file.path(out_dir, "paga_obs.tsv")
  emb_path <- file.path(out_dir, "paga_embedding.tsv")
  h5ad_path <- file.path(out_dir, "paga_input.h5ad")
  exporter_path <- file.path(out_dir, "build_minimal_paga_input.py")
  obsm_key <- paste0("X_", smcti_sanitize_name(reduction_name))
  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)[colnames(seurat_obj), , drop = FALSE]
  emb_df <- data.frame(cell_id = rownames(emb), emb, check.names = FALSE, stringsAsFactors = FALSE)
  smcti_write_tsv(obs, obs_path)
  smcti_write_tsv(emb_df, emb_path)
  writeLines(c(
    "#!/usr/bin/env python3",
    "import sys",
    "import pandas as pd",
    "import numpy as np",
    "import anndata as ad",
    sprintf("obs = pd.read_csv(r'%s', sep='\\t')", obs_path),
    sprintf("emb = pd.read_csv(r'%s', sep='\\t')", emb_path),
    "obs = obs.set_index('cell_id')",
    "emb = emb.set_index('cell_id')",
    "obs = obs.loc[emb.index]",
    "for col in obs.columns:",
    "    if obs[col].dtype == object:",
    "        obs[col] = obs[col].fillna('').astype(str)",
    "X = np.zeros((obs.shape[0], 1), dtype=np.float32)",
    "var = pd.DataFrame(index=['dummy_feature'])",
    "adata = ad.AnnData(X=X, obs=obs, var=var)",
    sprintf("adata.obsm['%s'] = emb.values.astype('float32')", obsm_key),
    sprintf("adata.write_h5ad(r'%s', compression='gzip')", h5ad_path)
  ), con = exporter_path, useBytes = TRUE)
  Sys.chmod(exporter_path, mode = "0755")
  status <- system2(python_bin, args = exporter_path)
  if (!identical(status, 0L) || !file.exists(h5ad_path)) {
    stop("Failed to build minimal h5ad input for PAGA", call. = FALSE)
  }
  list(h5ad_path = h5ad_path, obsm_key = obsm_key, obs_path = obs_path, emb_path = emb_path, exporter_path = exporter_path)
}

smcti_build_paga_plot_tables <- function(out_dir,
                                         topology_packet,
                                         groupby_key) {
  mat <- as.matrix(topology_packet$connectivity_matrix)
  heatmap_df <- expand.grid(
    from = rownames(mat),
    to = colnames(mat),
    stringsAsFactors = FALSE
  )
  heatmap_df$weight <- as.numeric(mat[cbind(match(heatmap_df$from, rownames(mat)), match(heatmap_df$to, colnames(mat)))])

  obs <- utils::read.delim(file.path(out_dir, "paga_obs.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  emb <- utils::read.delim(file.path(out_dir, "paga_embedding.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  coord_cols <- setdiff(colnames(emb), "cell_id")
  if (length(coord_cols) < 2L) stop("PAGA embedding needs at least two coordinates for plotting", call. = FALSE)
  cell_df <- merge(obs, emb[, c("cell_id", coord_cols[1:2]), drop = FALSE], by = "cell_id", all.x = TRUE, sort = FALSE)
  names(cell_df)[names(cell_df) == coord_cols[1]] <- "dim1"
  names(cell_df)[names(cell_df) == coord_cols[2]] <- "dim2"
  cell_df[[groupby_key]] <- as.character(cell_df[[groupby_key]])
  cell_df$dim1 <- suppressWarnings(as.numeric(cell_df$dim1))
  cell_df$dim2 <- suppressWarnings(as.numeric(cell_df$dim2))

  node_xy <- stats::aggregate(cbind(dim1, dim2) ~ group_label, data = data.frame(group_label = cell_df[[groupby_key]], dim1 = cell_df$dim1, dim2 = cell_df$dim2), FUN = mean, na.rm = TRUE)
  node_counts <- as.data.frame(table(cell_df[[groupby_key]]), stringsAsFactors = FALSE)
  colnames(node_counts) <- c("group_label", "n_cells")
  node_df <- merge(node_xy, node_counts, by = "group_label", all.x = TRUE, sort = FALSE)
  node_df$fraction_cells <- node_df$n_cells / sum(node_df$n_cells)

  edge_df <- topology_packet$edge_table
  if (nrow(edge_df) > 0L) {
    edge_plot <- merge(edge_df, node_df[, c("group_label", "dim1", "dim2", "n_cells"), drop = FALSE], by.x = "from", by.y = "group_label", all.x = TRUE, sort = FALSE)
    names(edge_plot)[names(edge_plot) %in% c("dim1", "dim2", "n_cells")] <- c("from_dim1", "from_dim2", "from_n_cells")
    edge_plot <- merge(edge_plot, node_df[, c("group_label", "dim1", "dim2", "n_cells"), drop = FALSE], by.x = "to", by.y = "group_label", all.x = TRUE, sort = FALSE)
    names(edge_plot)[names(edge_plot) %in% c("dim1", "dim2", "n_cells")] <- c("to_dim1", "to_dim2", "to_n_cells")
  } else {
    edge_plot <- data.frame(from = character(), to = character(), weight = numeric(), from_dim1 = numeric(), from_dim2 = numeric(), from_n_cells = numeric(), to_dim1 = numeric(), to_dim2 = numeric(), to_n_cells = numeric(), stringsAsFactors = FALSE)
  }

  list(heatmap = heatmap_df, cells = cell_df, nodes = node_df, edges = edge_plot)
}

smcti_plot_paga_heatmap <- function(heatmap_df, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = to, y = from, fill = weight)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", weight)), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", name = "PAGA\nconnectivity") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = "PAGA cluster connectivity heatmap", x = "To cluster", y = "From cluster") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(output_path, p, width = 8, height = 7)
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 8, height = 7, dpi = 180)
  invisible(output_path)
}

smcti_plot_paga_network <- function(node_df, edge_df, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  p <- ggplot2::ggplot()
  if (nrow(edge_df) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, linewidth = weight, alpha = weight),
      colour = "#525252",
      lineend = "round"
    ) + ggplot2::scale_linewidth(range = c(0.6, 3.5)) + ggplot2::scale_alpha(range = c(0.35, 0.9))
  }
  p <- p +
    ggplot2::geom_point(data = node_df, ggplot2::aes(x = dim1, y = dim2, size = n_cells, fill = group_label), shape = 21, colour = "black", alpha = 0.95) +
    ggplot2::geom_text(data = node_df, ggplot2::aes(x = dim1, y = dim2, label = group_label), vjust = -1.1, size = 3.5) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = "PAGA topology projected onto mean embedding", x = "Mean embedding dim1", y = "Mean embedding dim2", size = "Cells") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "right")
  ggplot2::ggsave(output_path, p, width = 9, height = 7)
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 9, height = 7, dpi = 180)
  invisible(output_path)
}

smcti_plot_paga_embedding <- function(cell_df, groupby_key, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  p <- ggplot2::ggplot(cell_df, ggplot2::aes_string(x = "dim1", y = "dim2", colour = groupby_key)) +
    ggplot2::geom_point(size = 0.65, alpha = 0.75) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = "PAGA input embedding colored by group", x = "Embedding dim1", y = "Embedding dim2", colour = groupby_key) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  ggplot2::ggsave(output_path, p, width = 9, height = 7)
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 9, height = 7, dpi = 180)
  invisible(output_path)
}

smcti_write_paga_visual_bundle <- function(out_dir,
                                           topology_packet,
                                           groupby_key) {
  figure_manifest <- file.path(out_dir, "paga_figure_companion_manifest.tsv")
  if (file.exists(figure_manifest)) unlink(figure_manifest)
  plot_tables <- smcti_build_paga_plot_tables(out_dir, topology_packet = topology_packet, groupby_key = groupby_key)
  smcti_write_tsv(plot_tables$heatmap, file.path(out_dir, "paga_connectivity_heatmap_data.tsv"))
  smcti_write_tsv(plot_tables$cells, file.path(out_dir, "paga_embedding_group_plot_data.tsv"))
  smcti_write_tsv(plot_tables$nodes, file.path(out_dir, "paga_network_nodes.tsv"))
  smcti_write_tsv(plot_tables$edges, file.path(out_dir, "paga_network_edges_for_plot.tsv"))

  heatmap_pdf <- file.path(out_dir, "paga_connectivity_heatmap.pdf")
  smcti_plot_paga_heatmap(plot_tables$heatmap, heatmap_pdf)
  smcanno_viz_register_figure(
    figure_path = heatmap_pdf,
    data = plot_tables$heatmap,
    method = "PAGA",
    figure_type = "connectivity_heatmap",
    title = "PAGA connectivity matrix heatmap",
    extra_context = c(sprintf("Group key: %s", groupby_key), sprintf("PAGA groups: %s; retained topology edges: %s", topology_packet$n_groups, topology_packet$n_edges), "Connectivity values are undirected cluster-level graph connectivities from Scanpy PAGA."),
    manifest_path = figure_manifest
  )

  network_pdf <- file.path(out_dir, "paga_network_graph.pdf")
  smcti_plot_paga_network(plot_tables$nodes, plot_tables$edges, network_pdf)
  network_data <- rbind(
    data.frame(row_type = "node", group = plot_tables$nodes$group_label, partner = NA_character_, weight = NA_real_, dim1 = plot_tables$nodes$dim1, dim2 = plot_tables$nodes$dim2, n_cells = plot_tables$nodes$n_cells, stringsAsFactors = FALSE),
    data.frame(row_type = "edge", group = plot_tables$edges$from, partner = plot_tables$edges$to, weight = plot_tables$edges$weight, dim1 = plot_tables$edges$from_dim1, dim2 = plot_tables$edges$from_dim2, n_cells = plot_tables$edges$from_n_cells, stringsAsFactors = FALSE)
  )
  smcanno_viz_register_figure(
    figure_path = network_pdf,
    data = network_data,
    method = "PAGA",
    figure_type = "projected_network_graph",
    title = "PAGA projected network graph",
    extra_context = c("Node positions are mean embedding coordinates of cells in each group; edge width/alpha reflects PAGA connectivity.", "Use this as topology support, not as pseudotime direction by itself."),
    manifest_path = figure_manifest
  )

  embedding_pdf <- file.path(out_dir, "paga_embedding_groups.pdf")
  smcti_plot_paga_embedding(plot_tables$cells, groupby_key = groupby_key, output_path = embedding_pdf)
  smcanno_viz_register_figure(
    figure_path = embedding_pdf,
    data = plot_tables$cells,
    method = "PAGA",
    figure_type = "embedding_colored_by_group",
    title = "PAGA input embedding colored by cell type",
    extra_context = c("This plot verifies whether PAGA groups occupy coherent regions in the input embedding.", "Separated groups with weak edges should be interpreted differently from adjacent groups with strong edges."),
    manifest_path = figure_manifest
  )

  list(status = "ok", plot_tables = plot_tables, manifest = figure_manifest)
}

smcti_run_paga_bundle <- function(seurat_obj,
                                  output_dir,
                                  groupby_key,
                                  reduction_name,
                                  assay_name,
                                  python_bin = "/home/h2048/miniconda3/bin/python") {
  out_dir <- smcti_ensure_dir(output_dir)
  res <- tryCatch({
    export_info <- smcti_write_minimal_h5ad_for_paga(seurat_obj, output_dir = out_dir, groupby_key = groupby_key,
                                                     reduction_name = reduction_name, assay_name = assay_name,
                                                     python_bin = python_bin)
    paga_res <- pa_run_paga_runner(adata_h5ad_path = export_info$h5ad_path, output_dir = out_dir,
                                   groupby_key = groupby_key, python_cmd = python_bin,
                                   use_rep = export_info$obsm_key, n_neighbors = 20L, random_seed = 42L)
    if (!identical(paga_res$status, "ok") || is.null(paga_res$topology_packet)) {
      stop("PAGA runner did not return a topology packet", call. = FALSE)
    }
    edge_tbl <- paga_res$topology_packet$edge_table
    smcti_write_tsv(edge_tbl, file.path(out_dir, "paga_edge_table.tsv"))
    saveRDS(paga_res$topology_packet, file.path(out_dir, "paga_topology_packet.rds"))
        visual_bundle <- smcti_write_paga_visual_bundle(out_dir = out_dir, topology_packet = paga_res$topology_packet, groupby_key = groupby_key)
    list(status = "ok", topology_packet = paga_res$topology_packet, manifest = paga_res$manifest,
          output_dir = out_dir, h5ad_path = export_info$h5ad_path, visual_bundle = visual_bundle)
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(msg, con = file.path(out_dir, "paga_error.txt"), useBytes = TRUE)
    list(status = "error", error = msg, output_dir = out_dir)
  })
  res
}

smcti_dynamic_gene_marker_panel <- c(
  "RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3",
  "ACTA2", "TAGLN", "CNN1", "MYH11", "MYL9", "TPM2",
  "COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2",
  "CXCL12", "IL6", "CCL2", "POSTN", "LRRC15"
)

smcti_compute_gene_pseudotime_stats <- function(expr_mat, pseudotime_vec) {
  stats_list <- lapply(seq_len(nrow(expr_mat)), function(i) {
    gene_expr <- as.numeric(expr_mat[i, ])
    keep <- is.finite(gene_expr) & is.finite(pseudotime_vec)
    if (sum(keep) < 5L) return(c(rho = NA_real_, p_value = NA_real_))
    ct <- suppressWarnings(stats::cor.test(gene_expr[keep], pseudotime_vec[keep], method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p_value = ct$p.value)
  })
  out <- do.call(rbind, stats_list)
  rownames(out) <- rownames(expr_mat)
  out
}

smcti_run_dynamic_gene_bundle <- function(seurat_obj,
                                          slingshot_res,
                                          output_dir,
                                          assay_name,
                                          run_pathway = TRUE,
                                          max_genes = 2500L) {
  out_dir <- smcti_ensure_dir(output_dir)
  expr <- smcti_get_assay_matrix(seurat_obj, assay_name, slot_name = "data")
  genes_use <- unique(c(Seurat::VariableFeatures(seurat_obj), smcti_dynamic_gene_marker_panel))
  genes_use <- genes_use[genes_use %in% rownames(expr)]
  if (length(genes_use) == 0L) genes_use <- rownames(expr)
  genes_use <- utils::head(genes_use, max_genes)
  expr <- expr[genes_use, , drop = FALSE]

  method_used <- if (requireNamespace("tradeSeq", quietly = TRUE)) "tradeSeq" else "spearman_fallback"
  writeLines(sprintf("dynamic_gene_method=%s", method_used), con = file.path(out_dir, "dynamic_gene_method.txt"), useBytes = TRUE)

  dynamic_rows <- list()
  idx <- 1L
  enrichment_rows <- list()
  eidx <- 1L
  hallmark_sets <- NULL
  if (isTRUE(run_pathway) && requireNamespace("msigdbr", quietly = TRUE) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }

  for (lineage_name in slingshot_res$lineage_summary$lineage_id) {
    pt <- slingshot_res$pseudotime_table[[paste0(lineage_name, "_pseudotime")]]
    wt <- slingshot_res$pseudotime_table[[paste0(lineage_name, "_curve_weight")]]
    keep_cells <- which(is.finite(pt) & is.finite(wt) & wt >= 0.5)
    if (length(keep_cells) < 8L) next
    expr_sub <- as.matrix(expr[, slingshot_res$pseudotime_table$cell_id[keep_cells], drop = FALSE])
    stat_mat <- smcti_compute_gene_pseudotime_stats(expr_sub, pt[keep_cells])
    res_tbl <- data.frame(
      lineage_id = lineage_name,
      gene = rownames(stat_mat),
      spearman_rho = stat_mat[, "rho"],
      p_value = stat_mat[, "p_value"],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    res_tbl$padj <- stats::p.adjust(res_tbl$p_value, method = "BH")
    res_tbl <- res_tbl[order(abs(res_tbl$spearman_rho), decreasing = TRUE, na.last = TRUE), ]
    dynamic_rows[[idx]] <- res_tbl
    idx <- idx + 1L
    smcti_write_tsv(res_tbl, file.path(out_dir, sprintf("dynamic_genes_%s.tsv", smcti_sanitize_name(lineage_name))))

    if (requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)) {
      top_pos <- head(res_tbl$gene[which(is.finite(res_tbl$spearman_rho) & res_tbl$spearman_rho > 0)], 15)
      top_neg <- head(res_tbl$gene[which(is.finite(res_tbl$spearman_rho) & res_tbl$spearman_rho < 0)], 15)
      heatmap_genes <- unique(c(rev(top_neg), top_pos))
      heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(expr_sub)]
      if (length(heatmap_genes) >= 4L) {
        ord <- order(pt[keep_cells])
        mat <- expr_sub[heatmap_genes, ord, drop = FALSE]
        mat <- t(scale(t(mat)))
        mat[!is.finite(mat)] <- 0
        grDevices::pdf(file.path(out_dir, sprintf("dynamic_genes_%s_heatmap.pdf", smcti_sanitize_name(lineage_name))), width = 12, height = 8)
        ComplexHeatmap::draw(ComplexHeatmap::Heatmap(mat, name = "zscore", cluster_rows = FALSE, cluster_columns = FALSE,
                                                     show_column_names = FALSE,
                                                     col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
                                                     column_title = paste(lineage_name, "ordered by pseudotime")))
        grDevices::dev.off()
      }
    }

    if (!is.null(hallmark_sets)) {
      pos_genes <- res_tbl$gene[which(is.finite(res_tbl$padj) & res_tbl$padj < 0.05 & res_tbl$spearman_rho > 0.20)]
      neg_genes <- res_tbl$gene[which(is.finite(res_tbl$padj) & res_tbl$padj < 0.05 & res_tbl$spearman_rho < -0.20)]
      for (direction in c("positive", "negative")) {
        gene_vec <- if (identical(direction, "positive")) pos_genes else neg_genes
        if (length(gene_vec) < 5L) next
        enr <- suppressWarnings(clusterProfiler::enricher(
          gene = unique(gene_vec),
          TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
          TERM2NAME = unique(hallmark_sets[, c("gs_name", "gs_name")]),
          pAdjustMethod = "BH"
        ))
        if (!is.null(enr) && nrow(as.data.frame(enr)) > 0L) {
          enr_df <- as.data.frame(enr)
          enr_df$lineage_id <- lineage_name
          enr_df$direction <- direction
          enrichment_rows[[eidx]] <- enr_df
          eidx <- eidx + 1L
          smcti_write_tsv(enr_df, file.path(out_dir, sprintf("pathway_dynamics_%s_%s.tsv", smcti_sanitize_name(lineage_name), direction)))
        }
      }
    }
  }

  combined_dynamic <- if (length(dynamic_rows) > 0L) do.call(rbind, dynamic_rows) else data.frame()
  smcti_write_tsv(combined_dynamic, file.path(out_dir, "dynamic_genes_combined.tsv"))
  if (length(enrichment_rows) > 0L) {
    combined_enrichment <- do.call(rbind, enrichment_rows)
    smcti_write_tsv(combined_enrichment, file.path(out_dir, "pathway_dynamics_combined.tsv"))
  } else {
    combined_enrichment <- data.frame()
    smcti_write_tsv(combined_enrichment, file.path(out_dir, "pathway_dynamics_combined.tsv"))
  }
  list(status = "ok", method = method_used, dynamic_genes = combined_dynamic, pathway_dynamics = combined_enrichment, output_dir = out_dir)
}

smcti_infer_phenotype_col <- function(meta_data,
                                      candidates = c("disease_status", "disease", "phenotype", "status", "group", "Group", "condition", "diagnosis", "case_control", "Disease", "Condition"),
                                      min_levels = 2L,
                                      max_levels = 12L) {
  candidates <- unique(candidates[candidates %in% colnames(meta_data)])
  if (length(candidates) == 0L) return(list(column = NA_character_, candidates = data.frame()))
  rows <- lapply(candidates, function(col) {
    vals <- as.character(meta_data[[col]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    n_unique <- length(unique(vals))
    disease_like <- any(grepl("control|healthy|normal|hc|case|disease|patient|crs|polyp|tumou?r|cancer|inflam|non", vals, ignore.case = TRUE))
    sample_like <- n_unique > max_levels || n_unique > max(20L, floor(nrow(meta_data) * 0.05))
    score <- as.integer(n_unique >= min_levels && n_unique <= max_levels) * 10L + as.integer(disease_like) * 5L - as.integer(sample_like) * 8L
    data.frame(
      column = col,
      n_non_missing = length(vals),
      n_unique = n_unique,
      example_values = paste(utils::head(unique(vals), 8L), collapse = ";"),
      disease_like = disease_like,
      sample_like = sample_like,
      score = score,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  cand_tbl <- do.call(rbind, rows)
  cand_tbl <- cand_tbl[order(cand_tbl$score, -cand_tbl$n_unique, decreasing = TRUE), , drop = FALSE]
  best <- cand_tbl$column[[1]]
  if (!is.finite(cand_tbl$score[[1]]) || cand_tbl$score[[1]] <= 0) best <- NA_character_
  list(column = best, candidates = cand_tbl)
}

smcti_filter_gene_panel_for_trends <- function(genes,
                                               seurat_obj,
                                               max_genes = 16L,
                                               exclude_noninformative = TRUE) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  genes <- genes[genes %in% rownames(seurat_obj)]
  if (isTRUE(exclude_noninformative)) {
    reasons <- if (exists("smcanno_ts_gene_exclusion_reasons", mode = "function")) {
      smcanno_ts_gene_exclusion_reasons(genes)
    } else {
      gene_up <- toupper(genes)
      ifelse(
        grepl("^(IG[HKL]|IGJ|JCHAIN|MT[-.]|RPL|RPS|MRPL|MRPS|ENSG|LINC|AC[0-9]|AL[0-9]|AP[0-9]|MIR|SNHG|SNORA|SNORD|RNU|RNVU|Y_RNA)", gene_up) |
          grepl("(^|[-.])AS[0-9]*$", gene_up) |
          gene_up %in% c("MALAT1", "NEAT1", "NEAT2", "XIST", "TSIX", "KCNQ1OT1", "H19"),
        "excluded", ""
      )
    }
    genes <- genes[!nzchar(reasons)]
  }
  utils::head(genes, as.integer(max_genes))
}

smcti_build_method_pseudotime_long <- function(pseudotime_table,
                                               method_name,
                                               seurat_obj = NULL,
                                               cluster_col = NULL,
                                               phenotype_col = NULL,
                                               weight_threshold = 0.5) {
  pseudotime_table <- as.data.frame(pseudotime_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(pseudotime_table) == 0L || !"cell_id" %in% colnames(pseudotime_table)) return(data.frame())
  if (!is.null(seurat_obj)) {
    add_cols <- unique(c(cluster_col, phenotype_col))
    add_cols <- add_cols[!is.na(add_cols) & nzchar(add_cols) & add_cols %in% colnames(seurat_obj@meta.data)]
    for (col in add_cols) {
      if (!col %in% colnames(pseudotime_table)) {
        pseudotime_table[[col]] <- as.character(seurat_obj@meta.data[pseudotime_table$cell_id, col, drop = TRUE])
      }
    }
  }
  pt_cols <- grep("(^pseudotime$|_pseudotime$)", colnames(pseudotime_table), value = TRUE)
  if (length(pt_cols) == 0L) return(data.frame())
  rows <- lapply(pt_cols, function(pt_col) {
    trajectory_id <- sub("_pseudotime$", "", pt_col)
    if (identical(trajectory_id, pt_col)) trajectory_id <- method_name
    wt_col <- paste0(trajectory_id, "_curve_weight")
    curve_weight <- if (wt_col %in% colnames(pseudotime_table)) suppressWarnings(as.numeric(pseudotime_table[[wt_col]])) else rep(1, nrow(pseudotime_table))
    df <- data.frame(
      method = method_name,
      trajectory_id = trajectory_id,
      cell_id = as.character(pseudotime_table$cell_id),
      pseudotime = suppressWarnings(as.numeric(pseudotime_table[[pt_col]])),
      curve_weight = curve_weight,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!is.null(cluster_col) && !is.na(cluster_col) && cluster_col %in% colnames(pseudotime_table)) df$cluster_label <- as.character(pseudotime_table[[cluster_col]])
    if (!is.null(phenotype_col) && !is.na(phenotype_col) && phenotype_col %in% colnames(pseudotime_table)) df$phenotype <- as.character(pseudotime_table[[phenotype_col]])
    df <- df[is.finite(df$pseudotime), , drop = FALSE]
    if (wt_col %in% colnames(pseudotime_table)) df <- df[is.finite(df$curve_weight) & df$curve_weight >= weight_threshold, , drop = FALSE]
    if (nrow(df) == 0L) return(df)
    rng <- range(df$pseudotime, finite = TRUE)
    df$scaled_pseudotime <- if (all(is.finite(rng)) && diff(rng) > 0) (df$pseudotime - rng[1]) / diff(rng) else NA_real_
    df
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

smcti_build_gene_trend_long <- function(seurat_obj,
                                        assay_name,
                                        pt_long,
                                        genes,
                                        slot_name = "data") {
  genes <- smcti_filter_gene_panel_for_trends(genes, seurat_obj = seurat_obj, max_genes = length(genes), exclude_noninformative = TRUE)
  if (length(genes) == 0L || nrow(pt_long) == 0L) return(data.frame())
  expr <- tryCatch(smcti_get_assay_matrix(seurat_obj, assay_name = assay_name, slot_name = slot_name), error = function(e) NULL)
  if (is.null(expr) || nrow(expr) == 0L) expr <- smcti_get_assay_matrix(seurat_obj, assay_name = assay_name, slot_name = "counts")
  genes <- genes[genes %in% rownames(expr)]
  cells <- unique(pt_long$cell_id[pt_long$cell_id %in% colnames(expr)])
  if (length(genes) == 0L || length(cells) == 0L) return(data.frame())
  expr_mat <- as.matrix(expr[genes, cells, drop = FALSE])
  rows <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    df <- pt_long[pt_long$cell_id %in% cells, , drop = FALSE]
    df$gene <- genes[[i]]
    df$expression <- as.numeric(expr_mat[genes[[i]], df$cell_id])
    rows[[i]] <- df
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

smcti_build_gene_trend_summary <- function(trend_long, n_bins = 30L, by_phenotype = FALSE) {
  trend_long <- as.data.frame(trend_long, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(trend_long) == 0L) return(data.frame())
  trend_long$bin_id <- pmax(1L, pmin(as.integer(n_bins), floor(suppressWarnings(as.numeric(trend_long$scaled_pseudotime)) * as.integer(n_bins)) + 1L))
  trend_long$bin_mid <- (trend_long$bin_id - 0.5) / as.integer(n_bins)
  group_cols <- c("method", "trajectory_id", "gene", "bin_id", "bin_mid")
  if (isTRUE(by_phenotype) && "phenotype" %in% colnames(trend_long)) group_cols <- c(group_cols, "phenotype")
  key <- do.call(interaction, c(trend_long[group_cols], list(drop = TRUE, sep = "|||")))
  rows <- lapply(split(trend_long, key), function(df) {
    first_vals <- df[1, group_cols, drop = FALSE]
    data.frame(
      first_vals,
      n_cells = nrow(df),
      mean_pseudotime = smcti_safe_numeric_summary(df$pseudotime, mean),
      mean_scaled_pseudotime = smcti_safe_numeric_summary(df$scaled_pseudotime, mean),
      mean_curve_weight = smcti_safe_numeric_summary(df$curve_weight, mean),
      mean_expression = smcti_safe_numeric_summary(df$expression, mean),
      median_expression = smcti_safe_numeric_summary(df$expression, stats::median),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$method, out$trajectory_id, out$gene, out$bin_id), , drop = FALSE]
}

smcti_add_smoothed_expression_to_summary <- function(summary_df,
                                                     smoother_method = "smooth_spline",
                                                     smoother_min_points = 5L) {
  summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_df) == 0L) {
    summary_df$smoothed_expression <- numeric(0)
    summary_df$smoothing_applied <- logical(0)
    summary_df$smoother_method <- character(0)
    return(summary_df)
  }
  if (!all(c("mean_scaled_pseudotime", "mean_expression") %in% colnames(summary_df))) {
    stop("summary_df must contain mean_scaled_pseudotime and mean_expression", call. = FALSE)
  }

  smoother_method <- trimws(as.character(smoother_method[[1]]))
  if (!nzchar(smoother_method)) smoother_method <- "none"
  smoother_min_points <- max(3L, as.integer(smoother_min_points))

  summary_df$smoothed_expression <- suppressWarnings(as.numeric(summary_df$mean_expression))
  summary_df$smoothing_applied <- FALSE
  summary_df$smoother_method <- "none"

  group_cols <- intersect(c("method", "trajectory_id", "gene", "phenotype"), colnames(summary_df))
  if (length(group_cols) == 0L) group_cols <- c("gene")
  split_key <- do.call(interaction, c(summary_df[group_cols], list(drop = TRUE, sep = "|||")))

  for (group_idx in split(seq_len(nrow(summary_df)), split_key)) {
    df <- summary_df[group_idx, , drop = FALSE]
    ord <- order(df$mean_scaled_pseudotime, df$bin_id, na.last = TRUE)
    df <- df[ord, , drop = FALSE]
    row_idx <- group_idx[ord]
    x <- suppressWarnings(as.numeric(df$mean_scaled_pseudotime))
    y <- suppressWarnings(as.numeric(df$mean_expression))
    keep <- is.finite(x) & is.finite(y)
    if (!any(keep)) next

    x_keep <- x[keep]
    y_keep <- y[keep]
    row_keep <- row_idx[keep]
    applied <- FALSE
    y_smooth <- y_keep

    if (identical(smoother_method, "smooth_spline") && length(y_keep) >= smoother_min_points && length(unique(x_keep)) >= smoother_min_points) {
      fit <- tryCatch(stats::smooth.spline(x = x_keep, y = y_keep), error = function(e) NULL)
      pred <- if (!is.null(fit)) tryCatch(stats::predict(fit, x = x_keep)$y, error = function(e) NULL) else NULL
      if (!is.null(pred) && length(pred) == length(y_keep) && all(is.finite(pred))) {
        y_smooth <- as.numeric(pred)
        applied <- TRUE
      }
    }

    summary_df$smoothed_expression[row_keep] <- y_smooth
    summary_df$smoothing_applied[row_keep] <- applied
    summary_df$smoother_method[row_keep] <- if (applied) smoother_method else "none"
  }

  summary_df$smoothed_expression[!is.finite(summary_df$smoothed_expression) & is.finite(summary_df$mean_expression)] <-
    suppressWarnings(as.numeric(summary_df$mean_expression[!is.finite(summary_df$smoothed_expression) & is.finite(summary_df$mean_expression)]))
  summary_df$smoother_method[is.na(summary_df$smoother_method) | !nzchar(summary_df$smoother_method)] <- "none"
  summary_df$smoothing_applied[is.na(summary_df$smoothing_applied)] <- FALSE
  summary_df
}

smcti_build_gene_trend_plot_data <- function(summary_df,
                                             smoother_method = "smooth_spline",
                                             smoother_grid_n = 201L,
                                             smoother_min_points = 5L) {
  summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_df) == 0L) return(data.frame())
  if (!"smoothed_expression" %in% colnames(summary_df)) {
    summary_df <- smcti_add_smoothed_expression_to_summary(summary_df, smoother_method = smoother_method, smoother_min_points = smoother_min_points)
  }

  observed_df <- summary_df
  observed_df$plot_layer <- "observed_bin"
  observed_df$plot_expression <- suppressWarnings(as.numeric(observed_df$mean_expression))
  observed_df$line_group <- interaction(observed_df$trajectory_id, observed_df$gene,
                                        if ("phenotype" %in% colnames(observed_df)) observed_df$phenotype else "all",
                                        drop = TRUE, sep = "||")

  curve_rows <- list()
  idx <- 1L
  group_cols <- intersect(c("method", "trajectory_id", "gene", "phenotype"), colnames(summary_df))
  if (length(group_cols) == 0L) group_cols <- c("gene")
  split_key <- do.call(interaction, c(summary_df[group_cols], list(drop = TRUE, sep = "|||")))
  smoother_grid_n <- max(11L, as.integer(smoother_grid_n))
  smoother_min_points <- max(3L, as.integer(smoother_min_points))

  for (group_idx in split(seq_len(nrow(summary_df)), split_key)) {
    df <- summary_df[group_idx, , drop = FALSE]
    df <- df[order(df$mean_scaled_pseudotime, df$bin_id, na.last = TRUE), , drop = FALSE]
    x <- suppressWarnings(as.numeric(df$mean_scaled_pseudotime))
    y <- suppressWarnings(as.numeric(df$mean_expression))
    keep <- is.finite(x) & is.finite(y)
    if (!any(keep)) next
    x_keep <- x[keep]
    y_keep <- y[keep]
    fit_ok <- FALSE
    x_plot <- x_keep
    y_plot <- suppressWarnings(as.numeric(df$smoothed_expression[keep]))
    if (!all(is.finite(y_plot))) y_plot <- y_keep

    if (identical(smoother_method, "smooth_spline") && length(y_keep) >= smoother_min_points && length(unique(x_keep)) >= smoother_min_points) {
      fit <- tryCatch(stats::smooth.spline(x = x_keep, y = y_keep), error = function(e) NULL)
      if (!is.null(fit)) {
        x_grid <- seq(min(x_keep), max(x_keep), length.out = smoother_grid_n)
        pred <- tryCatch(stats::predict(fit, x = x_grid)$y, error = function(e) NULL)
        if (!is.null(pred) && length(pred) == length(x_grid) && all(is.finite(pred))) {
          x_plot <- x_grid
          y_plot <- as.numeric(pred)
          fit_ok <- TRUE
        }
      }
    }

    curve_df <- df[rep(which(keep)[1], length(x_plot)), , drop = FALSE]
    curve_df$bin_id <- NA_integer_
    curve_df$bin_mid <- NA_real_
    curve_df$n_cells <- NA_integer_
    curve_df$mean_pseudotime <- NA_real_
    curve_df$mean_curve_weight <- NA_real_
    curve_df$mean_expression <- NA_real_
    curve_df$median_expression <- NA_real_
    curve_df$smoothed_expression <- y_plot
    curve_df$plot_layer <- "smoothed_curve"
    curve_df$plot_expression <- y_plot
    curve_df$mean_scaled_pseudotime <- x_plot
    curve_df$smoothing_applied <- fit_ok
    curve_df$smoother_method <- if (fit_ok) smoother_method else "none"
    curve_df$line_group <- interaction(curve_df$trajectory_id, curve_df$gene,
                                       if ("phenotype" %in% colnames(curve_df)) curve_df$phenotype else "all",
                                       drop = TRUE, sep = "||")
    curve_rows[[idx]] <- curve_df
    idx <- idx + 1L
  }

  curve_df <- if (length(curve_rows) > 0L) do.call(rbind, curve_rows) else observed_df[0, , drop = FALSE]
  plot_df <- rbind(observed_df, curve_df)
  rownames(plot_df) <- NULL
  plot_df
}

smcti_plot_gene_trend_summary <- function(summary_df,
                                          output_path,
                                          title_text,
                                          color_col = "trajectory_id",
                                          smoother_method = "smooth_spline",
                                          smoother_grid_n = 201L,
                                          smoother_min_points = 5L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_df) == 0L || !color_col %in% colnames(summary_df)) return(invisible(NULL))
  plot_df <- smcti_build_gene_trend_plot_data(summary_df, smoother_method = smoother_method,
                                              smoother_grid_n = smoother_grid_n,
                                              smoother_min_points = smoother_min_points)
  if (nrow(plot_df) == 0L) return(invisible(NULL))
  observed_df <- plot_df[plot_df$plot_layer == "observed_bin", , drop = FALSE]
  curve_df <- plot_df[plot_df$plot_layer == "smoothed_curve", , drop = FALSE]
  multi_traj <- length(unique(summary_df$trajectory_id)) > 1L && !identical(color_col, "trajectory_id")

  p <- ggplot2::ggplot()
  if (nrow(curve_df) > 0L) {
    if (isTRUE(multi_traj)) {
      p <- p + ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes_string(x = "mean_scaled_pseudotime", y = "plot_expression", colour = color_col, group = "line_group", linetype = "trajectory_id"),
        linewidth = 0.9,
        alpha = 0.95
      )
    } else {
      p <- p + ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes_string(x = "mean_scaled_pseudotime", y = "plot_expression", colour = color_col, group = "line_group"),
        linewidth = 0.9,
        alpha = 0.95
      )
    }
  }
  p <- p + ggplot2::geom_point(
    data = observed_df,
    ggplot2::aes_string(x = "mean_scaled_pseudotime", y = "mean_expression", colour = color_col, group = "line_group"),
    size = 0.9,
    alpha = 0.75
  ) +
    ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = title_text, x = "Scaled pseudotime", y = "Mean expression", colour = color_col) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "right")
  if (isTRUE(multi_traj)) {
    p <- p + ggplot2::labs(linetype = "Trajectory")
  }
  ggplot2::ggsave(output_path, p, width = 13, height = 9)
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 13, height = 9, dpi = 180)
  invisible(plot_df)
}

smcti_write_gene_pseudotime_trend_bundle <- function(seurat_obj,
                                                     pseudotime_table,
                                                     output_dir,
                                                     assay_name,
                                                     method_name,
                                                     gene_panel,
                                                     cluster_col = NULL,
                                                     phenotype_col = NULL,
                                                     max_genes = 16L,
                                                     weight_threshold = 0.5,
                                                     trend_smoother_method = "smooth_spline",
                                                     trend_smoother_grid_n = 201L,
                                                     trend_smoother_min_points = 5L) {
  out_dir <- smcti_ensure_dir(output_dir)
  manifest_path <- file.path(out_dir, sprintf("%s_gene_trend_figure_companion_manifest.tsv", smcti_sanitize_name(method_name)))
  if (file.exists(manifest_path)) unlink(manifest_path)
  genes_use <- smcti_filter_gene_panel_for_trends(gene_panel, seurat_obj = seurat_obj, max_genes = max_genes, exclude_noninformative = TRUE)
  if (length(genes_use) == 0L) {
    smcti_write_tsv(data.frame(), file.path(out_dir, "gene_trend_selected_genes.tsv"))
    return(list(status = "no_genes", output_dir = out_dir, manifest = manifest_path))
  }
  smcti_write_tsv(data.frame(gene = genes_use, stringsAsFactors = FALSE), file.path(out_dir, "gene_trend_selected_genes.tsv"))
  pt_long <- smcti_build_method_pseudotime_long(pseudotime_table = pseudotime_table, method_name = method_name, seurat_obj = seurat_obj, cluster_col = cluster_col, phenotype_col = phenotype_col, weight_threshold = weight_threshold)
  trend_long <- smcti_build_gene_trend_long(seurat_obj = seurat_obj, assay_name = assay_name, pt_long = pt_long, genes = genes_use, slot_name = "data")
  smcti_write_tsv(trend_long, file.path(out_dir, "gene_pseudotime_trend_long.tsv"))
  trend_summary <- smcti_build_gene_trend_summary(trend_long, n_bins = 30L, by_phenotype = FALSE)
  trend_summary <- smcti_add_smoothed_expression_to_summary(trend_summary,
                                                            smoother_method = trend_smoother_method,
                                                            smoother_min_points = trend_smoother_min_points)
  smcti_write_tsv(trend_summary, file.path(out_dir, "gene_pseudotime_trend_summary.tsv"))
  if (nrow(trend_summary) > 0L) {
    trend_plot_data <- smcti_build_gene_trend_plot_data(trend_summary,
                                                        smoother_method = trend_smoother_method,
                                                        smoother_grid_n = trend_smoother_grid_n,
                                                        smoother_min_points = trend_smoother_min_points)
    smcti_write_tsv(trend_plot_data, file.path(out_dir, "gene_pseudotime_trend_plot_data.tsv"))
    trend_pdf <- file.path(out_dir, sprintf("%s_gene_pseudotime_trends.pdf", smcti_sanitize_name(method_name)))
    smcti_plot_gene_trend_summary(trend_summary, trend_pdf,
                                  sprintf("%s gene expression trends along pseudotime", method_name),
                                  color_col = "trajectory_id",
                                  smoother_method = trend_smoother_method,
                                  smoother_grid_n = trend_smoother_grid_n,
                                  smoother_min_points = trend_smoother_min_points)
    smcanno_viz_register_figure(
      figure_path = trend_pdf,
      data = trend_plot_data,
      method = method_name,
      figure_type = "gene_expression_pseudotime_trends",
      title = sprintf("%s gene expression trends along pseudotime", method_name),
      extra_context = c(
        "Observed points show 30-bin mean expression summaries against scaled pseudotime for selected SMC/pericyte/ECM dynamic genes.",
        sprintf("A display-layer `%s` smoother is overlaid for readability; this is not equivalent to a tradeSeq fitted GAM.", trend_smoother_method),
        "IG/MT/RP/ENSG/lncRNA-like genes are excluded from this trend panel.",
        if (!is.null(phenotype_col) && !is.na(phenotype_col)) sprintf("Disease phenotype column detected for companion comparisons: %s", phenotype_col) else "No phenotype column was provided."
      ),
      manifest_path = manifest_path
    )
  }
  phenotype_summary <- data.frame()
  if (!is.null(phenotype_col) && !is.na(phenotype_col) && "phenotype" %in% colnames(trend_long)) {
    phenotype_summary <- smcti_build_gene_trend_summary(trend_long, n_bins = 30L, by_phenotype = TRUE)
    phenotype_summary <- smcti_add_smoothed_expression_to_summary(phenotype_summary,
                                                                  smoother_method = trend_smoother_method,
                                                                  smoother_min_points = trend_smoother_min_points)
    smcti_write_tsv(phenotype_summary, file.path(out_dir, "gene_pseudotime_trend_by_phenotype_summary.tsv"))
    if (nrow(phenotype_summary) > 0L && length(unique(phenotype_summary$phenotype)) >= 2L) {
      phenotype_plot_data <- smcti_build_gene_trend_plot_data(phenotype_summary,
                                                              smoother_method = trend_smoother_method,
                                                              smoother_grid_n = trend_smoother_grid_n,
                                                              smoother_min_points = trend_smoother_min_points)
      smcti_write_tsv(phenotype_plot_data, file.path(out_dir, "gene_pseudotime_trend_by_phenotype_plot_data.tsv"))
      phenotype_pdf <- file.path(out_dir, sprintf("%s_gene_pseudotime_trends_by_%s.pdf", smcti_sanitize_name(method_name), smcti_sanitize_name(phenotype_col)))
      smcti_plot_gene_trend_summary(phenotype_summary, phenotype_pdf,
                                    sprintf("%s gene trends by %s", method_name, phenotype_col),
                                    color_col = "phenotype",
                                    smoother_method = trend_smoother_method,
                                    smoother_grid_n = trend_smoother_grid_n,
                                    smoother_min_points = trend_smoother_min_points)
      smcanno_viz_register_figure(
        figure_path = phenotype_pdf,
        data = phenotype_plot_data,
        method = method_name,
        figure_type = "gene_expression_pseudotime_trends_by_phenotype",
        title = sprintf("%s gene expression trends by %s", method_name, phenotype_col),
        extra_context = c(
          "Observed points compare binned mean expression along scaled pseudotime across disease phenotype groups.",
          sprintf("The overlaid `%s` curve is a display smoother for readability, not a tradeSeq-like fitted model.", trend_smoother_method),
          "Use the paired CSV to identify phenotype-specific divergence rather than judging only from line overlap.",
          sprintf("Phenotype column: %s", phenotype_col)
        ),
        manifest_path = manifest_path
      )
    }
  } else {
    smcti_write_tsv(phenotype_summary, file.path(out_dir, "gene_pseudotime_trend_by_phenotype_summary.tsv"))
  }
  list(status = "ok", output_dir = out_dir, manifest = manifest_path, genes = genes_use, trend_summary = trend_summary, phenotype_summary = phenotype_summary)
}

smcti_build_phenotype_pseudotime_summary <- function(method_tables,
                                                     seurat_obj,
                                                     phenotype_col,
                                                     cluster_col = NULL,
                                                     weight_threshold = 0.5) {
  if (is.null(phenotype_col) || is.na(phenotype_col) || !phenotype_col %in% colnames(seurat_obj@meta.data)) return(data.frame())
  rows <- list()
  idx <- 1L
  for (method_name in names(method_tables)) {
    pt_long <- smcti_build_method_pseudotime_long(method_tables[[method_name]], method_name = method_name, seurat_obj = seurat_obj, cluster_col = cluster_col, phenotype_col = phenotype_col, weight_threshold = weight_threshold)
    if (nrow(pt_long) == 0L || !"phenotype" %in% colnames(pt_long)) next
    key <- interaction(pt_long$method, pt_long$trajectory_id, pt_long$phenotype, drop = TRUE, sep = "|||")
    for (df in split(pt_long, key)) {
      rows[[idx]] <- data.frame(
        method = df$method[[1]],
        trajectory_id = df$trajectory_id[[1]],
        phenotype = df$phenotype[[1]],
        n_cells = nrow(df),
        mean_pseudotime = smcti_safe_numeric_summary(df$pseudotime, mean),
        median_pseudotime = smcti_safe_numeric_summary(df$pseudotime, stats::median),
        q25_pseudotime = smcti_safe_quantile(df$pseudotime, 0.25),
        q75_pseudotime = smcti_safe_quantile(df$pseudotime, 0.75),
        mean_scaled_pseudotime = smcti_safe_numeric_summary(df$scaled_pseudotime, mean),
        median_scaled_pseudotime = smcti_safe_numeric_summary(df$scaled_pseudotime, stats::median),
        q25_scaled_pseudotime = smcti_safe_quantile(df$scaled_pseudotime, 0.25),
        q75_scaled_pseudotime = smcti_safe_quantile(df$scaled_pseudotime, 0.75),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      idx <- idx + 1L
    }
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$method, out$trajectory_id, out$phenotype), , drop = FALSE]
}

smcti_plot_phenotype_pseudotime_summary <- function(summary_df, output_path, phenotype_col) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_df) == 0L) return(invisible(NULL))
  summary_df$method_trajectory <- paste(summary_df$method, summary_df$trajectory_id, sep = ": ")
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = phenotype, y = median_scaled_pseudotime, colour = phenotype)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = q25_scaled_pseudotime, ymax = q75_scaled_pseudotime), size = 0.4) +
    ggplot2::facet_wrap(~ method_trajectory, scales = "free_y") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = sprintf("Pseudotime distribution summary by %s", phenotype_col), x = phenotype_col, y = "Median scaled pseudotime") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
  ggplot2::ggsave(output_path, p, width = 12, height = 8)
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 12, height = 8, dpi = 180)
  invisible(output_path)
}

smcti_write_phenotype_pseudotime_comparison_bundle <- function(method_tables,
                                                               seurat_obj,
                                                               output_dir,
                                                               phenotype_col,
                                                               cluster_col = NULL,
                                                               weight_threshold = 0.5) {
  out_dir <- smcti_ensure_dir(output_dir)
  manifest_path <- file.path(out_dir, "phenotype_pseudotime_figure_companion_manifest.tsv")
  if (file.exists(manifest_path)) unlink(manifest_path)
  summary_df <- smcti_build_phenotype_pseudotime_summary(method_tables, seurat_obj = seurat_obj, phenotype_col = phenotype_col, cluster_col = cluster_col, weight_threshold = weight_threshold)
  smcti_write_tsv(summary_df, file.path(out_dir, "phenotype_pseudotime_summary.tsv"))
  if (nrow(summary_df) > 0L) {
    pdf_path <- file.path(out_dir, sprintf("phenotype_pseudotime_comparison_by_%s.pdf", smcti_sanitize_name(phenotype_col)))
    smcti_plot_phenotype_pseudotime_summary(summary_df, pdf_path, phenotype_col = phenotype_col)
    smcanno_viz_register_figure(
      figure_path = pdf_path,
      data = summary_df,
      method = "Cross-method trajectory comparison",
      figure_type = "phenotype_pseudotime_distribution_comparison",
      title = sprintf("Pseudotime distribution by %s across methods", phenotype_col),
      extra_context = c("Slingshot branch-specific pseudotime uses cells with curve weight >= 0.5; Monocle2/3 use their scalar pseudotime.", "PAGA is not included here because it is a topology graph method, not a single-cell pseudotime estimator.", sprintf("Disease phenotype column: %s", phenotype_col)),
      manifest_path = manifest_path
    )
  }
  list(status = if (nrow(summary_df) > 0L) "ok" else "empty", summary = summary_df, output_dir = out_dir, manifest = manifest_path)
}

smcti_build_cross_method_summary <- function(slingshot_res,
                                             monocle3_res = NULL,
                                             monocle2_res = NULL,
                                             paga_res = NULL) {
  rows <- list()
  rows[[1]] <- data.frame(
    method = "Slingshot",
    status = "ok",
    n_cells = nrow(slingshot_res$pseudotime_table),
    n_cells_with_pseudotime = sum(apply(slingshot_res$pseudotime_mat, 1, function(x) any(is.finite(x)))),
    topology_or_state = paste(sort(unique(slingshot_res$lineage_summary$terminal_state)), collapse = ";"),
    cluster_median_rank_cor_vs_slingshot = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  add_method <- function(method_name, method_res, pt_col = "pseudotime") {
    if (is.null(method_res)) return(NULL)
    if (!identical(method_res$status, "ok")) {
      return(data.frame(method = method_name, status = method_res$status, n_cells = NA_integer_, n_cells_with_pseudotime = NA_integer_,
                        topology_or_state = NA_character_, cluster_median_rank_cor_vs_slingshot = NA_real_, stringsAsFactors = FALSE, check.names = FALSE))
    }
    ref_cluster_pt <- aggregate(slingshot_res$pseudotime_table[[paste0(slingshot_res$lineage_summary$lineage_id[[1]], "_pseudotime")]],
                                by = list(cluster_label = slingshot_res$pseudotime_table[[slingshot_res$cluster_col]]), FUN = function(x) smcti_safe_numeric_summary(x, stats::median))
    colnames(ref_cluster_pt)[2] <- "ref_median_pseudotime"
    method_cluster_pt <- aggregate(method_res$pseudotime_table[[pt_col]], by = list(cluster_label = method_res$pseudotime_table$cluster_label),
                                   FUN = function(x) smcti_safe_numeric_summary(x, stats::median))
    colnames(method_cluster_pt)[2] <- "method_median_pseudotime"
    merged <- merge(ref_cluster_pt, method_cluster_pt, by = "cluster_label")
    rho <- if (nrow(merged) >= 2L) suppressWarnings(stats::cor(merged$ref_median_pseudotime, merged$method_median_pseudotime, method = "spearman", use = "pairwise.complete.obs")) else NA_real_
    data.frame(
      method = method_name,
      status = "ok",
      n_cells = nrow(method_res$pseudotime_table),
      n_cells_with_pseudotime = sum(is.finite(method_res$pseudotime_table[[pt_col]])),
      topology_or_state = if (identical(method_name, "Monocle3")) paste(sort(unique(method_res$pseudotime_table$partition)), collapse = ";") else paste(sort(unique(method_res$pseudotime_table$state)), collapse = ";"),
      cluster_median_rank_cor_vs_slingshot = rho,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (!is.null(monocle3_res)) rows[[length(rows) + 1L]] <- add_method("Monocle3", monocle3_res, pt_col = "pseudotime")
  if (!is.null(monocle2_res)) rows[[length(rows) + 1L]] <- add_method("Monocle2", monocle2_res, pt_col = "pseudotime")
  if (!is.null(paga_res)) {
    if (!identical(paga_res$status, "ok")) {
      rows[[length(rows) + 1L]] <- data.frame(method = "PAGA", status = paga_res$status, n_cells = NA_integer_, n_cells_with_pseudotime = NA_integer_, topology_or_state = NA_character_, cluster_median_rank_cor_vs_slingshot = NA_real_, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      rows[[length(rows) + 1L]] <- data.frame(method = "PAGA", status = "ok", n_cells = length(slingshot_res$pseudotime_table$cell_id), n_cells_with_pseudotime = NA_integer_,
                                              topology_or_state = sprintf("n_groups=%s;n_edges=%s", paga_res$topology_packet$n_groups, paga_res$topology_packet$n_edges),
                                              cluster_median_rank_cor_vs_slingshot = NA_real_, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  do.call(rbind, rows)
}
