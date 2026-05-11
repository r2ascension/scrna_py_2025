#!/usr/bin/env Rscript
# ==============================================================================
# Program Support Helper (2026-04-29 v1)
# ============================================================================== 

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)

PA_SUPPORT_HELPER_VERSION_20260429_V1 <- "20260429_v1"

pa_roe_to_symbol <- function(roe_value, thresholds = c(0.5, 0.8, 1.2, 2.0)) {
  if (is.null(roe_value) || length(roe_value) == 0L || is.na(roe_value) || is.nan(roe_value)) return("-")
  if (is.infinite(roe_value)) return(ifelse(roe_value > 0, "+++", "---"))
  if (!is.numeric(roe_value)) return("-")
  if (roe_value >= thresholds[[4]]) return("+++")
  if (roe_value >= thresholds[[3]]) return("+")
  if (roe_value > thresholds[[2]]) return("+/-")
  if (roe_value >= thresholds[[1]]) return("-")
  "---"
}

pa_roe_interpretation <- function(roe_value) {
  if (is.na(roe_value) || !is.finite(roe_value)) {
    "Data not available"
  } else if (roe_value >= 2.0) {
    "Strongly enriched"
  } else if (roe_value >= 1.5) {
    "Moderately enriched"
  } else if (roe_value >= 1.2) {
    "Mildly enriched"
  } else if (roe_value >= 0.8) {
    "Near expected levels"
  } else if (roe_value >= 0.5) {
    "Mildly depleted"
  } else if (roe_value >= 0.3) {
    "Moderately depleted"
  } else {
    "Strongly depleted"
  }
}

pa_extract_metadata_df <- function(input_obj) {
  if (is.data.frame(input_obj)) {
    meta <- input_obj
    if (!"cell_id" %in% colnames(meta)) {
      meta$cell_id <- if (!is.null(rownames(meta))) rownames(meta) else sprintf("cell_%d", seq_len(nrow(meta)))
    }
    return(meta)
  }
  if (inherits(input_obj, "Seurat")) {
    meta <- input_obj@meta.data
    meta$cell_id <- rownames(meta)
    return(meta)
  }
  stop("input_obj must be either a data.frame or a Seurat object", call. = FALSE)
}

pa_get_seurat_assay_matrix <- function(seurat_obj,
                                       assay_name = NULL,
                                       preferred_layers = c("counts", "data")) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object", call. = FALSE)
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required", call. = FALSE)
  }
  assay_name <- pa_null_coalesce(assay_name, Seurat::DefaultAssay(seurat_obj))
  layer_names <- tryCatch({
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      SeuratObject::Layers(seurat_obj[[assay_name]])
    } else {
      character()
    }
  }, error = function(e) character())

  chosen_layer <- preferred_layers[preferred_layers %in% layer_names][1]
  if (!is.na(chosen_layer) && nzchar(chosen_layer)) {
    return(Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = chosen_layer))
  }

  for (slot_name in preferred_layers) {
    mat <- tryCatch(
      Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = slot_name),
      error = function(e) NULL
    )
    if (!is.null(mat)) return(mat)
  }

  stop(sprintf("Could not resolve a usable assay layer/slot for assay '%s'", assay_name), call. = FALSE)
}

pa_fill_numeric_matrix <- function(mat, fill_value = NULL) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (all(is.finite(mat))) return(mat)
  finite_vals <- as.numeric(mat[is.finite(mat)])
  replacement <- if (!is.null(fill_value) && is.finite(fill_value)) {
    as.numeric(fill_value)
  } else if (length(finite_vals) > 0L) {
    stats::median(finite_vals)
  } else {
    0
  }
  mat[!is.finite(mat)] <- replacement
  mat
}

pa_plot_matrix_heatmap <- function(mat,
                                   output_prefix,
                                   title,
                                   display_numbers = NULL,
                                   color_palette = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
                                   na_color = "#F0F0F0") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    return(invisible(NULL))
  }
  safe_mat <- pa_fill_numeric_matrix(mat)
  finite_vals <- as.numeric(safe_mat[is.finite(safe_mat)])
  if (length(finite_vals) == 0L) finite_vals <- 0
  val_range <- range(finite_vals)
  breaks <- if (diff(val_range) < .Machine$double.eps) {
    seq(val_range[[1]] - 0.5, val_range[[2]] + 0.5, length.out = length(color_palette) + 1L)
  } else {
    seq(val_range[[1]], val_range[[2]], length.out = length(color_palette) + 1L)
  }
  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)

  grDevices::pdf(pdf_path, width = max(8, ncol(safe_mat) * 1.2), height = max(6, nrow(safe_mat) * 0.45))
  pheatmap::pheatmap(
    safe_mat,
    color = color_palette,
    breaks = breaks,
    display_numbers = display_numbers,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    main = title,
    na_col = na_color,
    fontsize = 10,
    fontsize_number = 8
  )
  grDevices::dev.off()

  grDevices::png(png_path, width = 1800, height = 1200, res = 180)
  pheatmap::pheatmap(
    safe_mat,
    color = color_palette,
    breaks = breaks,
    display_numbers = display_numbers,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    main = title,
    na_col = na_color,
    fontsize = 10,
    fontsize_number = 8
  )
  grDevices::dev.off()

  invisible(list(pdf = pdf_path, png = png_path))
}

pa_build_roe_packet <- function(observed_table,
                                expected_matrix,
                                roe_matrix,
                                roe_long_table,
                                row_summary,
                                global_test_summary,
                                row_key,
                                col_key) {
  list(
    row_key = pa_scalar_chr(row_key, "row_key"),
    col_key = pa_scalar_chr(col_key, "col_key"),
    observed_table = observed_table,
    expected_matrix = expected_matrix,
    roe_matrix = roe_matrix,
    roe_long_table = roe_long_table,
    row_summary = row_summary,
    global_test_summary = global_test_summary,
    helper_version = PA_SUPPORT_HELPER_VERSION_20260429_V1
  )
}

pa_run_roe_runner <- function(input_obj,
                              output_dir,
                              row_key,
                              col_key,
                              min_row_sum = 10L,
                              pseudocount = 0,
                              create_heatmap = TRUE,
                              heatmap_prefix = "roe") {
  meta <- pa_extract_metadata_df(input_obj)
  row_key <- pa_scalar_chr(row_key, "row_key")
  col_key <- pa_scalar_chr(col_key, "col_key")
  output_dir <- pa_prepare_output_dir(output_dir)

  pa_validate_required_columns(meta, c(row_key, col_key), "Ro/e metadata")
  row_vals <- pa_safe_trim(meta[[row_key]])
  col_vals <- pa_safe_trim(meta[[col_key]])
  keep <- nzchar(row_vals) & nzchar(col_vals)
  if (!any(keep)) {
    stop("No non-empty observations available for Ro/e calculation", call. = FALSE)
  }

  observed_table <- table(row_vals[keep], col_vals[keep], useNA = "no")
  observed_table <- observed_table[rowSums(observed_table) >= as.integer(min_row_sum), , drop = FALSE]
  observed_table <- observed_table[, colSums(observed_table) > 0, drop = FALSE]
  if (nrow(observed_table) == 0L || ncol(observed_table) == 0L) {
    stop("Observed table is empty after filtering", call. = FALSE)
  }

  total_n <- sum(observed_table)
  expected_matrix <- outer(rowSums(observed_table), colSums(observed_table)) / total_n
  rownames(expected_matrix) <- rownames(observed_table)
  colnames(expected_matrix) <- colnames(observed_table)
  roe_matrix <- (observed_table + as.numeric(pseudocount)) / (expected_matrix + as.numeric(pseudocount))
  std_resid <- (observed_table - expected_matrix) / sqrt(pmax(expected_matrix, .Machine$double.eps))

  use_simulated_p <- min(expected_matrix) < 1 || mean(expected_matrix < 5) > 0.2
  global_test <- suppressWarnings(stats::chisq.test(observed_table, simulate.p.value = use_simulated_p, B = if (use_simulated_p) 10000 else 2000))
  global_test_summary <- list(
    statistic = unname(global_test$statistic),
    parameter = if (!is.null(global_test$parameter)) unname(global_test$parameter) else NA_real_,
    p_value = unname(global_test$p.value),
    method = global_test$method,
    simulated_p_value = isTRUE(use_simulated_p)
  )

  roe_rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(observed_table))) {
    for (j in seq_len(ncol(observed_table))) {
      current_roe <- as.numeric(roe_matrix[i, j])
      roe_rows[[idx]] <- data.frame(
        row_level = rownames(observed_table)[[i]],
        col_level = colnames(observed_table)[[j]],
        observed = as.numeric(observed_table[i, j]),
        expected = as.numeric(expected_matrix[i, j]),
        roe = current_roe,
        symbol = pa_roe_to_symbol(current_roe),
        interpretation = pa_roe_interpretation(current_roe),
        std_residual = as.numeric(std_resid[i, j]),
        row_key = row_key,
        col_key = col_key,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  roe_long_table <- do.call(rbind, roe_rows)

  row_summary <- do.call(rbind, lapply(seq_len(nrow(roe_matrix)), function(i) {
    row_name <- rownames(roe_matrix)[[i]]
    row_values <- as.numeric(roe_matrix[i, ])
    max_idx <- which.max(row_values)
    min_idx <- which.min(row_values)
    data.frame(
      row_level = row_name,
      n_cells = as.numeric(rowSums(observed_table)[[i]]),
      dominant_col = colnames(roe_matrix)[[max_idx]],
      dominant_roe = row_values[[max_idx]],
      depleted_col = colnames(roe_matrix)[[min_idx]],
      depleted_roe = row_values[[min_idx]],
      stringsAsFactors = FALSE
    )
  }))

  symbol_matrix <- matrix(
    vapply(as.numeric(roe_matrix), pa_roe_to_symbol, character(1)),
    nrow = nrow(roe_matrix),
    ncol = ncol(roe_matrix),
    dimnames = dimnames(roe_matrix)
  )

  observed_path <- file.path(output_dir, "roe_observed_counts.csv")
  expected_path <- file.path(output_dir, "roe_expected_counts.csv")
  roe_path <- file.path(output_dir, "roe_indices.csv")
  long_path <- file.path(output_dir, "roe_long.tsv")
  row_summary_path <- file.path(output_dir, "roe_row_summary.tsv")
  manifest_path <- file.path(output_dir, "roe_manifest.json")

  utils::write.csv(as.data.frame.matrix(observed_table), observed_path, quote = FALSE)
  utils::write.csv(as.data.frame(expected_matrix), expected_path, quote = FALSE)
  utils::write.csv(as.data.frame(roe_matrix), roe_path, quote = FALSE)
  pa_write_tsv(roe_long_table, long_path)
  pa_write_tsv(row_summary, row_summary_path)

  heatmap_paths <- NULL
  if (isTRUE(create_heatmap)) {
    heatmap_paths <- pa_plot_matrix_heatmap(
      mat = roe_matrix,
      output_prefix = file.path(output_dir, heatmap_prefix),
      title = sprintf("Ro/e heatmap: %s × %s", row_key, col_key),
      display_numbers = symbol_matrix
    )
  }

  manifest <- list(
    output_dir = output_dir,
    row_key = row_key,
    col_key = col_key,
    observed_path = observed_path,
    expected_path = expected_path,
    roe_path = roe_path,
    long_path = long_path,
    row_summary_path = row_summary_path,
    heatmap_paths = heatmap_paths,
    global_test_summary = global_test_summary
  )
  pa_write_json(manifest, manifest_path)

  list(
    status = "ok",
    roe_packet = pa_build_roe_packet(
      observed_table = observed_table,
      expected_matrix = expected_matrix,
      roe_matrix = roe_matrix,
      roe_long_table = roe_long_table,
      row_summary = row_summary,
      global_test_summary = global_test_summary,
      row_key = row_key,
      col_key = col_key
    ),
    manifest = manifest
  )
}

pa_run_ggtree_runner <- function(input_matrix,
                                 output_dir,
                                 tree_prefix = "matrix_tree",
                                 margins = c("rows", "cols"),
                                 fill_value = NULL,
                                 layout = "rectangular") {
  if (!requireNamespace("ggtree", quietly = TRUE)) stop("Package 'ggtree' is required", call. = FALSE)
  if (!requireNamespace("ape", quietly = TRUE)) stop("Package 'ape' is required", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required", call. = FALSE)

  output_dir <- pa_prepare_output_dir(output_dir)
  tree_prefix <- pa_scalar_chr(tree_prefix, "tree_prefix")
  mat <- pa_fill_numeric_matrix(input_matrix, fill_value = fill_value)
  margins <- unique(tolower(pa_safe_trim(margins)))
  margins[margins %in% c("row", "rows")] <- "rows"
  margins[margins %in% c("col", "cols", "columns")] <- "cols"
  margins <- intersect(margins, c("rows", "cols"))
  if (length(margins) == 0L) stop("margins must include at least one of 'rows' or 'cols'", call. = FALSE)

  summary_rows <- list()
  idx <- 1L
  for (margin in margins) {
    cluster_mat <- if (identical(margin, "rows")) mat else t(mat)
    tip_labels <- rownames(cluster_mat)
    if (is.null(tip_labels) || !length(tip_labels)) {
      tip_labels <- sprintf("%s_%d", margin, seq_len(nrow(cluster_mat)))
      rownames(cluster_mat) <- tip_labels
    }
    if (nrow(cluster_mat) < 2L) {
      summary_rows[[idx]] <- data.frame(
        margin = margin,
        n_tips = nrow(cluster_mat),
        status = "skipped",
        pdf_path = NA_character_,
        png_path = NA_character_,
        newick_path = NA_character_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
      next
    }

    hc <- stats::hclust(stats::dist(cluster_mat))
    phy <- ape::as.phylo(hc)
    pdf_path <- file.path(output_dir, sprintf("%s_%s_tree.pdf", tree_prefix, margin))
    png_path <- file.path(output_dir, sprintf("%s_%s_tree.png", tree_prefix, margin))
    newick_path <- file.path(output_dir, sprintf("%s_%s_tree.newick", tree_prefix, margin))

    tree_plot <- ggtree::ggtree(phy, layout = layout) +
      ggtree::geom_tiplab(size = 3) +
      ggplot2::ggtitle(sprintf("ggtree %s clustering: %s", margin, tree_prefix))

    ggplot2::ggsave(pdf_path, plot = tree_plot, width = 8, height = max(6, length(phy$tip.label) * 0.3))
    ggplot2::ggsave(png_path, plot = tree_plot, width = 8, height = max(6, length(phy$tip.label) * 0.3), dpi = 200)
    ape::write.tree(phy, file = newick_path)

    summary_rows[[idx]] <- data.frame(
      margin = margin,
      n_tips = length(phy$tip.label),
      status = "ok",
      pdf_path = pdf_path,
      png_path = png_path,
      newick_path = newick_path,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  tree_summary <- do.call(rbind, summary_rows)
  pa_write_tsv(tree_summary, file.path(output_dir, sprintf("%s_tree_summary.tsv", tree_prefix)))
  pa_write_json(
    list(
      tree_prefix = tree_prefix,
      output_dir = output_dir,
      margins = margins,
      tree_summary = tree_summary
    ),
    file.path(output_dir, sprintf("%s_tree_manifest.json", tree_prefix))
  )

  list(status = if (any(tree_summary$status == "ok")) "ok" else "empty", tree_summary = tree_summary)
}

pa_with_tibble_attached <- function(expr) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for ROGUE compatibility in this environment", call. = FALSE)
  }
  attached_before <- "package:tibble" %in% search()
  if (!attached_before) {
    suppressPackageStartupMessages(library(tibble))
    on.exit({
      if ("package:tibble" %in% search()) {
        try(detach("package:tibble", unload = TRUE, character.only = TRUE), silent = TRUE)
      }
    }, add = TRUE)
  }
  eval.parent(substitute(expr))
}

pa_calc_rogue_safe <- function(mat,
                               min.cells = 10,
                               min.genes = 200,
                               platform = "UMI") {
  if (!requireNamespace("ROGUE", quietly = TRUE)) {
    stop("Package 'ROGUE' is required", call. = FALSE)
  }
  rogue_ns <- asNamespace("ROGUE")
  matr_filter <- get("matr.filter", envir = rogue_ns)
  se_fun <- get("SE_fun", envir = rogue_ns)
  calc_rogue <- get("CalculateRogue", envir = rogue_ns)

  compute_rogue <- function(expr_mat, method_label) {
    pa_with_tibble_attached({
      mat_filt <- matr_filter(expr_mat, min.cells = min.cells, min.genes = min.genes)
      if (is.null(dim(mat_filt)) || length(dim(mat_filt)) != 2L) stop("Filtered expression matrix lost two-dimensional structure")
      if (ncol(mat_filt) < min.cells || nrow(mat_filt) < max(2L, floor(min.genes / 2))) stop("Insufficient data after filtering")
      ent <- se_fun(mat_filt)
      if (any(is.na(ent$entropy)) || any(is.infinite(ent$entropy))) stop("Invalid entropy values")
      rogue_val <- calc_rogue(ent, platform = platform)
      if (!is.finite(rogue_val) || rogue_val < 0 || rogue_val > 1) stop(sprintf("Invalid ROGUE value: %.4f", rogue_val))
      list(rogue_value = as.numeric(rogue_val), n_genes = nrow(mat_filt), n_cells = ncol(mat_filt), method = method_label)
    })
  }

  sparse_result <- try({
    compute_rogue(mat, "sparse")
  }, silent = TRUE)
  if (!inherits(sparse_result, "try-error")) return(sparse_result)

  dense_result <- try({
    mat_dense <- as.matrix(mat)
    storage.mode(mat_dense) <- "numeric"
    compute_rogue(mat_dense, "dense")
  }, silent = TRUE)
  if (inherits(dense_result, "try-error")) return(NULL)
  dense_result
}

pa_run_rogue_runner <- function(seurat_obj,
                                output_dir,
                                group_col,
                                split_col = NULL,
                                assay_name = NULL,
                                preferred_layers = c("counts", "data"),
                                min_cells_per_group = 10L,
                                max_cells_per_group = 2000L,
                                rogue_min_cells = 10L,
                                rogue_min_genes = 200L,
                                platform = "UMI",
                                seed = 42L,
                                create_heatmap = TRUE) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object", call. = FALSE)
  output_dir <- pa_prepare_output_dir(output_dir)
  meta <- seurat_obj@meta.data
  required_cols <- c(group_col, if (!is.null(split_col)) split_col else character())
  pa_validate_required_columns(meta, required_cols, "ROGUE metadata")
  expr_all <- pa_get_seurat_assay_matrix(seurat_obj, assay_name = assay_name, preferred_layers = preferred_layers)

  group_vals <- pa_safe_trim(meta[[group_col]])
  split_vals <- if (!is.null(split_col)) pa_safe_trim(meta[[split_col]]) else rep("ALL", nrow(meta))
  cell_ids <- rownames(meta)
  groups <- sort(unique(group_vals[nzchar(group_vals)]))

  result_rows <- list()
  skipped_rows <- list()
  result_idx <- 1L
  skipped_idx <- 1L

  compute_one_subset <- function(current_cells, group_value, split_value, scope_label) {
    n_cells_total <- length(current_cells)
    if (n_cells_total < as.integer(min_cells_per_group)) {
      return(list(
        result = NULL,
        skipped = data.frame(
          group = group_value,
          split = split_value,
          scope = scope_label,
          n_cells_total = n_cells_total,
          reason = "too_few_cells",
          stringsAsFactors = FALSE
        )
      ))
    }

    sampled_cells <- current_cells
    if (n_cells_total > as.integer(max_cells_per_group)) {
      set.seed(as.integer(seed))
      sampled_cells <- sample(current_cells, as.integer(max_cells_per_group))
    }
    expr_subset <- expr_all[, sampled_cells, drop = FALSE]
    rogue_result <- pa_calc_rogue_safe(
      expr_subset,
      min.cells = as.integer(rogue_min_cells),
      min.genes = as.integer(rogue_min_genes),
      platform = platform
    )
    if (is.null(rogue_result)) {
      return(list(
        result = NULL,
        skipped = data.frame(
          group = group_value,
          split = split_value,
          scope = scope_label,
          n_cells_total = n_cells_total,
          reason = "rogue_failed",
          stringsAsFactors = FALSE
        )
      ))
    }

    list(
      result = data.frame(
        group = group_value,
        split = split_value,
        scope = scope_label,
        n_cells_total = n_cells_total,
        n_cells_used = length(sampled_cells),
        n_cells_after_filter = rogue_result$n_cells,
        n_genes_after_filter = rogue_result$n_genes,
        filter_fraction = rogue_result$n_cells / max(1, length(sampled_cells)),
        rogue_value = rogue_result$rogue_value,
        calc_method = rogue_result$method,
        stringsAsFactors = FALSE
      ),
      skipped = NULL
    )
  }

  for (group_value in groups) {
    group_cells <- cell_ids[group_vals == group_value]
    overall_res <- compute_one_subset(group_cells, group_value, "ALL", "overall")
    if (is.data.frame(overall_res$result)) {
      result_rows[[result_idx]] <- overall_res$result
      result_idx <- result_idx + 1L
    }
    if (is.data.frame(overall_res$skipped)) {
      skipped_rows[[skipped_idx]] <- overall_res$skipped
      skipped_idx <- skipped_idx + 1L
    }

    if (!is.null(split_col)) {
      current_splits <- sort(unique(split_vals[group_vals == group_value & nzchar(split_vals)]))
      for (split_value in current_splits) {
        split_cells <- cell_ids[group_vals == group_value & split_vals == split_value]
        split_res <- compute_one_subset(split_cells, group_value, split_value, "split")
        if (is.data.frame(split_res$result)) {
          result_rows[[result_idx]] <- split_res$result
          result_idx <- result_idx + 1L
        }
        if (is.data.frame(split_res$skipped)) {
          skipped_rows[[skipped_idx]] <- split_res$skipped
          skipped_idx <- skipped_idx + 1L
        }
      }
    }
  }

  rogue_table <- if (length(result_rows) > 0L) do.call(rbind, result_rows) else data.frame()
  skipped_table <- if (length(skipped_rows) > 0L) do.call(rbind, skipped_rows) else data.frame()

  rogue_summary <- if (nrow(rogue_table) > 0L) {
    do.call(rbind, lapply(split(rogue_table, rogue_table$group), function(df) {
      split_df <- df[df$scope == "split", , drop = FALSE]
      overall_df <- df[df$scope == "overall", , drop = FALSE]
      data.frame(
        group = df$group[[1]],
        overall_rogue = if (nrow(overall_df) > 0L) overall_df$rogue_value[[1]] else NA_real_,
        n_successful_scopes = nrow(df),
        n_split_success = nrow(split_df),
        median_split_rogue = if (nrow(split_df) > 0L) stats::median(split_df$rogue_value) else NA_real_,
        max_split_rogue = if (nrow(split_df) > 0L) max(split_df$rogue_value) else NA_real_,
        min_split_rogue = if (nrow(split_df) > 0L) min(split_df$rogue_value) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame()
  }

  rogue_wide <- data.frame()
  if (!is.null(split_col) && nrow(rogue_table) > 0L) {
    split_df <- rogue_table[rogue_table$scope == "split", c("group", "split", "rogue_value"), drop = FALSE]
    if (nrow(split_df) > 0L) {
      split_groups <- sort(unique(split_df$group))
      split_levels <- sort(unique(split_df$split))
      wide_mat <- matrix(NA_real_, nrow = length(split_groups), ncol = length(split_levels), dimnames = list(split_groups, split_levels))
      for (i in seq_len(nrow(split_df))) {
        wide_mat[split_df$group[[i]], split_df$split[[i]]] <- split_df$rogue_value[[i]]
      }
      rogue_wide <- data.frame(group = rownames(wide_mat), as.data.frame(wide_mat), check.names = FALSE, stringsAsFactors = FALSE)
      if (isTRUE(create_heatmap)) {
        pa_plot_matrix_heatmap(
          mat = wide_mat,
          output_prefix = file.path(output_dir, "rogue_heatmap"),
          title = sprintf("ROGUE heatmap: %s × %s", group_col, split_col),
          display_numbers = matrix(sprintf("%.2f", pa_fill_numeric_matrix(wide_mat)), nrow = nrow(wide_mat), dimnames = dimnames(wide_mat))
        )
      }
    }
  }

  long_path <- file.path(output_dir, "rogue_long.tsv")
  summary_path <- file.path(output_dir, "rogue_summary.tsv")
  skipped_path <- file.path(output_dir, "rogue_skipped.tsv")
  wide_path <- file.path(output_dir, "rogue_wide.csv")
  manifest_path <- file.path(output_dir, "rogue_manifest.json")
  pa_write_tsv(rogue_table, long_path)
  pa_write_tsv(rogue_summary, summary_path)
  pa_write_tsv(skipped_table, skipped_path)
  utils::write.csv(rogue_wide, wide_path, row.names = FALSE, quote = FALSE)
  pa_write_json(
    list(
      output_dir = output_dir,
      group_col = group_col,
      split_col = pa_null_coalesce(split_col, NA_character_),
      n_successful_rows = nrow(rogue_table),
      n_skipped_rows = nrow(skipped_table),
      long_path = long_path,
      summary_path = summary_path,
      skipped_path = skipped_path,
      wide_path = wide_path
    ),
    manifest_path
  )

  list(
    status = if (nrow(rogue_table) > 0L) "ok" else "empty",
    rogue_table = rogue_table,
    rogue_summary = rogue_summary,
    rogue_wide = rogue_wide,
    skipped_table = skipped_table,
    manifest = list(
      output_dir = output_dir,
      long_path = long_path,
      summary_path = summary_path,
      skipped_path = skipped_path,
      wide_path = wide_path
    )
  )
}

if (sys.nframe() == 0) {
  cat("Program Support Helper (2026-04-29 v1) loaded.\n")
}
