
# ============================================================================== 
# 8.3 Enhanced visualization helpers (v1.2)
# ============================================================================== 
# Added figures:
#   hdWGCNA:
#     - module size barplot
#     - module eigengene heatmap by metadata group
#     - module eigengene violin/boxplot by metadata group
#     - optional module dendrogram wrapper
#     - hub-gene rank barplot
#     - one-shot hdWGCNA visualization summary wrapper
#   CoVarNet:
#     - network size summary
#     - edge direction summary
#     - degree distribution
#     - per-celltype hub-rank barplot
#     - per-celltype signed adjacency heatmap
#     - per-celltype raw correlation heatmap
#     - one-shot CoVarNet visualization summary wrapper

coexpr_safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "unknown"
  x
}

coexpr_save_ggplot <- function(plot,
                               filename,
                               width  = FIG_WIDTH,
                               height = FIG_HEIGHT,
                               out_dir = OUTPUT_DIR) {
  ensure_dir(out_dir)
  pdf_path <- file.path(out_dir, filename)
  grDevices::pdf(pdf_path, width = width, height = height)
  print(plot)
  grDevices::dev.off()
  cat(sprintf("[OK] Plot saved: %s\n", pdf_path))
  invisible(pdf_path)
}

coexpr_zscore_rows <- function(mat) {
  z <- t(apply(mat, 1, function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sx
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

coexpr_get_hub_gene_col <- function(hub_genes_df) {
  candidates <- c("gene_name", "gene", "name")
  hit <- candidates[candidates %in% colnames(hub_genes_df)][1]
  if (is.na(hit)) stop("[ERROR] Hub table has no gene column: gene_name/gene/name")
  hit
}

coexpr_get_hub_score_col <- function(hub_genes_df) {
  candidates <- c("kME", "kME_rank", "hub_score", "degree", "connectivity")
  hit <- candidates[candidates %in% colnames(hub_genes_df)][1]
  if (is.na(hit)) {
    numeric_cols <- colnames(hub_genes_df)[vapply(hub_genes_df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, c("rank"))
    hit <- numeric_cols[1]
  }
  if (is.na(hit)) stop("[ERROR] Hub table has no numeric score column.")
  hit
}

# ----- 8.3.1 hdWGCNA module size barplot -----
hdwgcna_plot_module_sizes <- function(seurat_obj,
                                      out_dir    = OUTPUT_DIR,
                                      wgcna_name = WGCNA_NAME) {
  cat("\n--- hdWGCNA Visualization: Module Sizes ---\n")
  ensure_dir(out_dir)

  mods <- tryCatch(GetModules(seurat_obj, wgcna_name = wgcna_name),
                   error = function(e) NULL)
  if (is.null(mods) || !"module" %in% colnames(mods)) {
    cat("[WARN] Module table unavailable -- skipping module-size plot.\n")
    return(invisible(NULL))
  }

  size_df <- mods %>%
    dplyr::filter(module != "grey") %>%
    dplyr::count(module, name = "n_genes") %>%
    dplyr::arrange(dplyr::desc(n_genes))

  if (nrow(size_df) == 0L) {
    cat("[WARN] No non-grey modules -- skipping module-size plot.\n")
    return(invisible(NULL))
  }

  utils::write.csv(size_df,
                   file.path(out_dir, "hdwgcna_module_sizes.csv"),
                   row.names = FALSE)

  p <- ggplot(size_df,
              aes(x = stats::reorder(module, n_genes), y = n_genes)) +
    geom_col(width = 0.8, fill = "grey35") +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      title = "hdWGCNA Module Sizes",
      x = "Module", y = "Number of genes"
    )

  coexpr_save_ggplot(
    p,
    filename = "hdwgcna_module_sizes.pdf",
    width = FIG_WIDTH,
    height = max(5, 0.28 * nrow(size_df) + 2),
    out_dir = out_dir
  )
  invisible(size_df)
}

# ----- 8.3.2 hdWGCNA ME heatmap by group -----
hdwgcna_plot_me_group_heatmap <- function(seurat_obj,
                                           group_col  = CONDITION_COL,
                                           out_dir    = OUTPUT_DIR,
                                           wgcna_name = WGCNA_NAME,
                                           min_cells_per_group = 10L,
                                           max_groups = 40L,
                                           scale_rows = TRUE) {
  cat(sprintf("\n--- hdWGCNA Visualization: ME Heatmap by %s ---\n", group_col))
  ensure_dir(out_dir)

  meta <- seurat_obj@meta.data
  if (!group_col %in% colnames(meta)) {
    cat(sprintf("[WARN] Column '%s' not found -- skipping ME heatmap.\n", group_col))
    return(invisible(NULL))
  }

  me <- tryCatch(GetMEs(seurat_obj, wgcna_name = wgcna_name),
                 error = function(e) NULL)
  if (is.null(me) || ncol(me) == 0L) {
    cat("[WARN] Module eigengenes unavailable -- skipping ME heatmap.\n")
    return(invisible(NULL))
  }

  shared_cells <- intersect(rownames(me), rownames(meta))
  if (length(shared_cells) < 20L) {
    cat("[WARN] Too few shared cells -- skipping ME heatmap.\n")
    return(invisible(NULL))
  }

  plot_df <- data.frame(
    group = as.character(meta[shared_cells, group_col, drop = TRUE]),
    me[shared_cells, , drop = FALSE],
    check.names = FALSE
  )
  plot_df <- plot_df[!is.na(plot_df$group) & nzchar(plot_df$group), , drop = FALSE]

  group_counts <- sort(table(plot_df$group), decreasing = TRUE)
  keep_groups <- names(group_counts[group_counts >= min_cells_per_group])
  keep_groups <- head(keep_groups, max_groups)
  plot_df <- plot_df[plot_df$group %in% keep_groups, , drop = FALSE]

  if (length(unique(plot_df$group)) < 2L) {
    cat("[WARN] Need at least two valid groups -- skipping ME heatmap.\n")
    return(invisible(NULL))
  }

  mean_df <- stats::aggregate(. ~ group, data = plot_df, FUN = mean, na.rm = TRUE)
  rownames(mean_df) <- mean_df$group
  mean_df$group <- NULL
  mat <- t(as.matrix(mean_df))
  mat <- mat[apply(mat, 1, function(x) any(is.finite(x))), , drop = FALSE]
  plot_mat <- if (scale_rows) coexpr_zscore_rows(mat) else mat

  csv_path <- file.path(out_dir,
                        sprintf("hdwgcna_ME_mean_by_%s.csv", coexpr_safe_filename(group_col)))
  utils::write.csv(mat, csv_path)

  pdf_path <- file.path(out_dir,
                        sprintf("hdwgcna_ME_heatmap_by_%s.pdf", coexpr_safe_filename(group_col)))
  grDevices::pdf(pdf_path,
                 width = max(7, 0.35 * ncol(plot_mat) + 4),
                 height = max(6, 0.25 * nrow(plot_mat) + 3))
  pheatmap::pheatmap(
    plot_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(101),
    fontsize_row = 8,
    fontsize_col = 9,
    main = sprintf("Module Eigengenes by %s%s",
                   group_col,
                   ifelse(scale_rows, " (row z-score)", ""))
  )
  grDevices::dev.off()
  cat(sprintf("[OK] ME heatmap: %s\n", pdf_path))

  invisible(list(mean_matrix = mat, plot_matrix = plot_mat, csv = csv_path))
}

# ----- 8.3.3 hdWGCNA ME violin / boxplot by group -----
hdwgcna_plot_me_violin <- function(seurat_obj,
                                    group_col  = CONDITION_COL,
                                    out_dir    = OUTPUT_DIR,
                                    wgcna_name = WGCNA_NAME,
                                    max_modules = 16L,
                                    max_cells = 50000L,
                                    seed = 42L) {
  cat(sprintf("\n--- hdWGCNA Visualization: ME Violin by %s ---\n", group_col))
  ensure_dir(out_dir)

  meta <- seurat_obj@meta.data
  if (!group_col %in% colnames(meta)) {
    cat(sprintf("[WARN] Column '%s' not found -- skipping ME violin.\n", group_col))
    return(invisible(NULL))
  }

  me <- tryCatch(GetMEs(seurat_obj, wgcna_name = wgcna_name),
                 error = function(e) NULL)
  if (is.null(me) || ncol(me) == 0L) {
    cat("[WARN] Module eigengenes unavailable -- skipping ME violin.\n")
    return(invisible(NULL))
  }

  shared_cells <- intersect(rownames(me), rownames(meta))
  if (length(shared_cells) < 20L) {
    cat("[WARN] Too few shared cells -- skipping ME violin.\n")
    return(invisible(NULL))
  }

  me_sub <- me[shared_cells, , drop = FALSE]
  # Select high-variance MEs for readability.
  mod_vars <- apply(me_sub, 2, stats::var, na.rm = TRUE)
  keep_mods <- names(sort(mod_vars, decreasing = TRUE))[seq_len(min(max_modules, length(mod_vars)))]

  plot_df <- data.frame(
    cell  = shared_cells,
    group = as.character(meta[shared_cells, group_col, drop = TRUE]),
    me_sub[, keep_mods, drop = FALSE],
    check.names = FALSE
  )
  plot_df <- plot_df[!is.na(plot_df$group) & nzchar(plot_df$group), , drop = FALSE]

  if (!is.null(max_cells) && nrow(plot_df) > max_cells) {
    set.seed(seed)
    plot_df <- plot_df[sample(seq_len(nrow(plot_df)), max_cells), , drop = FALSE]
  }

  long_df <- plot_df %>%
    tidyr::pivot_longer(
      cols = -c(cell, group),
      names_to = "module",
      values_to = "ME"
    )

  p <- ggplot(long_df, aes(x = group, y = ME)) +
    geom_violin(scale = "width", trim = TRUE, fill = "grey85", color = "grey40") +
    geom_boxplot(width = 0.12, outlier.size = 0.1, alpha = 0.75) +
    facet_wrap(~ module, scales = "free_y", ncol = 4) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = sprintf("Module Eigengene Distribution by %s", group_col),
      x = group_col, y = "Module eigengene"
    )

  coexpr_save_ggplot(
    p,
    filename = sprintf("hdwgcna_ME_violin_by_%s.pdf", coexpr_safe_filename(group_col)),
    width = max(10, FIG_WIDTH * 1.5),
    height = max(8, ceiling(length(keep_mods) / 4) * 2.5 + 2),
    out_dir = out_dir
  )
  invisible(long_df)
}

# ----- 8.3.4 Optional hdWGCNA dendrogram wrapper -----
hdwgcna_plot_dendrogram <- function(seurat_obj,
                                     out_dir    = OUTPUT_DIR,
                                     wgcna_name = WGCNA_NAME) {
  cat("\n--- hdWGCNA Visualization: Module Dendrogram ---\n")
  ensure_dir(out_dir)

  if (!exists("PlotDendrogram", mode = "function")) {
    cat("[WARN] PlotDendrogram() not found in this hdWGCNA version -- skipping.\n")
    return(invisible(NULL))
  }

  p <- tryCatch(
    PlotDendrogram(seurat_obj, wgcna_name = wgcna_name),
    error = function(e) {
      cat(sprintf("[WARN] PlotDendrogram failed: %s\n", e$message))
      NULL
    }
  )
  if (is.null(p)) return(invisible(NULL))

  coexpr_save_ggplot(
    p,
    filename = "hdwgcna_module_dendrogram.pdf",
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    out_dir = out_dir
  )
  invisible(p)
}

# ----- 8.3.5 hdWGCNA hub-gene rank barplot -----
hdwgcna_plot_hub_rank <- function(hub_genes_df,
                                  out_dir = OUTPUT_DIR,
                                  n_top_per_module = 8L) {
  cat("\n--- hdWGCNA Visualization: Hub Gene Rank ---\n")
  ensure_dir(out_dir)

  if (is.null(hub_genes_df) || nrow(hub_genes_df) == 0L) {
    cat("[WARN] Empty hub gene table -- skipping hub-rank plot.\n")
    return(invisible(NULL))
  }
  if (!"module" %in% colnames(hub_genes_df)) {
    cat("[WARN] Hub gene table has no module column -- skipping hub-rank plot.\n")
    return(invisible(NULL))
  }

  gene_col  <- coexpr_get_hub_gene_col(hub_genes_df)
  score_col <- coexpr_get_hub_score_col(hub_genes_df)

  plot_df <- hub_genes_df %>%
    dplyr::filter(module != "grey") %>%
    dplyr::group_by(module) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]), .by_group = TRUE) %>%
    dplyr::slice_head(n = n_top_per_module) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene_label = .data[[gene_col]])

  if (nrow(plot_df) == 0L) {
    cat("[WARN] No non-grey hub genes -- skipping hub-rank plot.\n")
    return(invisible(NULL))
  }

  utils::write.csv(plot_df,
                   file.path(out_dir, "hdwgcna_top_hubs_for_plot.csv"),
                   row.names = FALSE)

  p <- ggplot(plot_df,
              aes(x = stats::reorder(gene_label, .data[[score_col]]),
                  y = .data[[score_col]])) +
    geom_col(fill = "grey35", width = 0.8) +
    coord_flip() +
    facet_wrap(~ module, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(
      title = sprintf("Top Hub Genes per Module (%s)", score_col),
      x = "Gene", y = score_col
    )

  coexpr_save_ggplot(
    p,
    filename = "hdwgcna_hub_rank_barplot.pdf",
    width = max(10, FIG_WIDTH * 1.4),
    height = max(8, 0.22 * nrow(plot_df) + 3),
    out_dir = out_dir
  )
  invisible(plot_df)
}

# ----- 8.3.6 hdWGCNA visualization summary wrapper -----
hdwgcna_plot_summary <- function(seurat_obj,
                                  hub_genes_df = NULL,
                                  group_cols = c(CONDITION_COL, TISSUE_COL),
                                  out_dir    = OUTPUT_DIR,
                                  wgcna_name = WGCNA_NAME,
                                  make_umap = TRUE,
                                  make_violin = TRUE) {
  cat("\n"); cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat("hdWGCNA: Enhanced visualization summary\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  ensure_dir(out_dir)

  res <- list()
  res$module_sizes <- hdwgcna_plot_module_sizes(seurat_obj, out_dir, wgcna_name)
  res$dendrogram   <- hdwgcna_plot_dendrogram(seurat_obj, out_dir, wgcna_name)

  valid_groups <- unique(group_cols[group_cols %in% colnames(seurat_obj@meta.data)])
  for (gc in valid_groups) {
    n_levels <- length(unique(na.omit(seurat_obj@meta.data[[gc]])))
    if (n_levels < 2L) next
    res[[paste0("ME_heatmap_", gc)]] <- hdwgcna_plot_me_group_heatmap(
      seurat_obj, group_col = gc, out_dir = out_dir, wgcna_name = wgcna_name
    )
    if (make_violin && n_levels <= 20L) {
      res[[paste0("ME_violin_", gc)]] <- hdwgcna_plot_me_violin(
        seurat_obj, group_col = gc, out_dir = out_dir, wgcna_name = wgcna_name
      )
    } else if (make_violin) {
      cat(sprintf("[WARN] %s has %d levels; violin skipped for readability.\n", gc, n_levels))
    }
  }

  if (make_umap) {
    res$module_umap <- hdwgcna_plot_module_umap(seurat_obj, out_dir, wgcna_name)
  }

  if (!is.null(hub_genes_df)) {
    res$hub_rank <- hdwgcna_plot_hub_rank(hub_genes_df, out_dir = out_dir)
    if (CONDITION_COL %in% colnames(seurat_obj@meta.data)) {
      res$hub_dotplot <- hdwgcna_hub_dotplot(
        seurat_obj,
        hub_genes_df = hub_genes_df,
        group_col = CONDITION_COL,
        out_dir = out_dir,
        wgcna_name = wgcna_name
      )
    }
  }

  invisible(res)
}

# ----- 8.4.1 CoVarNet network summary barplot -----
plot_covarnet_summary <- function(covarnet_results,
                                  out_dir = OUTPUT_DIR) {
  cat("\n--- CoVarNet Visualization: Network Summary ---\n")
  ensure_dir(out_dir)

  if (is.null(covarnet_results) || length(covarnet_results) == 0L) {
    cat("[WARN] Empty CoVarNet result list -- skipping summary.\n")
    return(invisible(NULL))
  }

  summary_df <- do.call(rbind, lapply(names(covarnet_results), function(ct) {
    res <- covarnet_results[[ct]]
    g <- res$graph
    if (is.null(g)) {
      return(data.frame(celltype = ct, n_nodes = 0L, n_edges = 0L,
                        mean_degree = NA_real_, n_communities = NA_integer_,
                        stringsAsFactors = FALSE))
    }
    data.frame(
      celltype = ct,
      n_nodes = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      mean_degree = mean(igraph::degree(g)),
      n_communities = length(unique(igraph::V(g)$community)),
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- summary_df[order(summary_df$n_edges, decreasing = TRUE), ]

  utils::write.csv(summary_df,
                   file.path(out_dir, "covarnet_network_summary.csv"),
                   row.names = FALSE)

  long_df <- summary_df %>%
    tidyr::pivot_longer(cols = c(n_nodes, n_edges),
                        names_to = "metric", values_to = "value")

  p <- ggplot(long_df,
              aes(x = stats::reorder(celltype, value), y = value, fill = metric)) +
    geom_col(position = "dodge", width = 0.75) +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(
      title = "CoVarNet Network Size Summary",
      x = "Cell type", y = "Count", fill = "Metric"
    )

  coexpr_save_ggplot(
    p,
    filename = "covarnet_network_summary.pdf",
    width = FIG_WIDTH,
    height = max(5, 0.35 * nrow(summary_df) + 2),
    out_dir = out_dir
  )
  invisible(summary_df)
}

# ----- 8.4.2 CoVarNet edge direction summary -----
plot_covarnet_edge_direction_summary <- function(covarnet_results,
                                                 out_dir = OUTPUT_DIR) {
  cat("\n--- CoVarNet Visualization: Edge Direction Summary ---\n")
  ensure_dir(out_dir)

  edge_df <- do.call(rbind, lapply(names(covarnet_results), function(ct) {
    x <- covarnet_results[[ct]]$edges
    if (is.null(x) || nrow(x) == 0L) return(NULL)
    x
  }))

  if (is.null(edge_df) || nrow(edge_df) == 0L) {
    cat("[WARN] No edges available -- skipping edge-direction summary.\n")
    return(invisible(NULL))
  }

  sum_df <- edge_df %>%
    dplyr::count(celltype, direction, name = "n_edges")

  utils::write.csv(sum_df,
                   file.path(out_dir, "covarnet_edge_direction_summary.csv"),
                   row.names = FALSE)

  p <- ggplot(sum_df, aes(x = stats::reorder(celltype, n_edges),
                          y = n_edges, fill = direction)) +
    geom_col(width = 0.8) +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(
      title = "CoVarNet Edge Direction Summary",
      x = "Cell type", y = "Number of edges", fill = "Direction"
    )

  coexpr_save_ggplot(
    p,
    filename = "covarnet_edge_direction_summary.pdf",
    width = FIG_WIDTH,
    height = max(5, 0.35 * length(unique(sum_df$celltype)) + 2),
    out_dir = out_dir
  )
  invisible(sum_df)
}

# ----- 8.4.3 CoVarNet degree distribution across cell types -----
plot_covarnet_degree_distribution <- function(covarnet_results,
                                              out_dir = OUTPUT_DIR) {
  cat("\n--- CoVarNet Visualization: Degree Distribution ---\n")
  ensure_dir(out_dir)

  deg_df <- do.call(rbind, lapply(names(covarnet_results), function(ct) {
    g <- covarnet_results[[ct]]$graph
    if (is.null(g) || igraph::vcount(g) == 0L) return(NULL)
    data.frame(
      celltype = ct,
      gene = igraph::V(g)$name,
      degree = igraph::degree(g),
      hub_score = igraph::V(g)$hub_score,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(deg_df) || nrow(deg_df) == 0L) {
    cat("[WARN] No graph nodes available -- skipping degree distribution.\n")
    return(invisible(NULL))
  }

  utils::write.csv(deg_df,
                   file.path(out_dir, "covarnet_node_degree_all.csv"),
                   row.names = FALSE)

  p <- ggplot(deg_df, aes(x = degree)) +
    geom_histogram(bins = 40, fill = "grey45", color = "white") +
    facet_wrap(~ celltype, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(
      title = "CoVarNet Degree Distribution",
      x = "Node degree", y = "Number of genes"
    )

  coexpr_save_ggplot(
    p,
    filename = "covarnet_degree_distribution.pdf",
    width = max(10, FIG_WIDTH * 1.4),
    height = max(8, ceiling(length(unique(deg_df$celltype)) / 3) * 2.5 + 2),
    out_dir = out_dir
  )
  invisible(deg_df)
}

# ----- 8.4.4 CoVarNet hub-rank barplot per cell type -----
plot_covarnet_hub_rank <- function(covarnet_result,
                                   celltype_name = "Unknown",
                                   n_top = N_HUB_GENES,
                                   out_dir = OUTPUT_DIR) {
  if (is.null(covarnet_result) || is.null(covarnet_result$hubs)) {
    return(invisible(NULL))
  }
  ensure_dir(out_dir)

  node_df <- covarnet_result$hubs$all_nodes
  if (is.null(node_df) || nrow(node_df) == 0L) return(invisible(NULL))

  plot_df <- node_df %>%
    dplyr::arrange(dplyr::desc(hub_score)) %>%
    dplyr::slice_head(n = n_top)

  p <- ggplot(plot_df, aes(x = stats::reorder(gene, hub_score), y = hub_score)) +
    geom_col(fill = "grey35", width = 0.8) +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(
      title = sprintf("Top CoVarNet Hubs: %s", celltype_name),
      x = "Gene", y = "Hub score"
    )

  coexpr_save_ggplot(
    p,
    filename = sprintf("covarnet_%s_hub_rank.pdf", coexpr_safe_filename(celltype_name)),
    width = FIG_WIDTH,
    height = max(5, 0.28 * nrow(plot_df) + 2),
    out_dir = out_dir
  )
  invisible(plot_df)
}

# ----- 8.4.5 CoVarNet signed adjacency heatmap per cell type -----
plot_covarnet_adjacency_heatmap <- function(graph,
                                            celltype_name = "Unknown",
                                            max_nodes = 80L,
                                            out_dir = OUTPUT_DIR) {
  if (is.null(graph) || igraph::vcount(graph) == 0L || igraph::ecount(graph) == 0L) {
    return(invisible(NULL))
  }
  ensure_dir(out_dir)

  deg <- igraph::degree(graph)
  top_nodes <- names(sort(deg, decreasing = TRUE))[seq_len(min(max_nodes, length(deg)))]
  g_sub <- igraph::induced_subgraph(graph, top_nodes)

  adj <- as.matrix(igraph::as_adjacency_matrix(g_sub, attr = "r", sparse = FALSE))
  adj <- adj[top_nodes[top_nodes %in% rownames(adj)], top_nodes[top_nodes %in% colnames(adj)], drop = FALSE]

  pdf_path <- file.path(out_dir,
                        sprintf("covarnet_%s_signed_adjacency_heatmap.pdf",
                                coexpr_safe_filename(celltype_name)))
  grDevices::pdf(pdf_path,
                 width = max(8, 0.14 * nrow(adj) + 3),
                 height = max(8, 0.14 * nrow(adj) + 3))
  pheatmap::pheatmap(
    adj,
    color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(101),
    breaks = seq(-1, 1, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 6,
    fontsize_col = 6,
    main = sprintf("Signed Co-variation Adjacency: %s", celltype_name)
  )
  grDevices::dev.off()
  cat(sprintf("[OK] Signed adjacency heatmap: %s\n", pdf_path))
  invisible(adj)
}

# ----- 8.4.6 CoVarNet raw correlation heatmap per cell type -----
plot_covarnet_correlation_heatmap <- function(cor_result,
                                              celltype_name = NULL,
                                              max_genes = 80L,
                                              out_dir = OUTPUT_DIR) {
  if (is.null(cor_result) || is.null(cor_result$cor_mat)) return(invisible(NULL))
  ensure_dir(out_dir)

  cor_mat <- cor_result$cor_mat
  ct <- if (!is.null(celltype_name)) celltype_name else cor_result$celltype

  if (nrow(cor_mat) > max_genes) {
    score <- rowSums(abs(cor_mat), na.rm = TRUE)
    keep <- names(sort(score, decreasing = TRUE))[seq_len(max_genes)]
    cor_mat <- cor_mat[keep, keep, drop = FALSE]
  }

  pdf_path <- file.path(out_dir,
                        sprintf("covarnet_%s_correlation_heatmap.pdf",
                                coexpr_safe_filename(ct)))
  grDevices::pdf(pdf_path,
                 width = max(8, 0.14 * nrow(cor_mat) + 3),
                 height = max(8, 0.14 * nrow(cor_mat) + 3))
  pheatmap::pheatmap(
    cor_mat,
    color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(101),
    breaks = seq(-1, 1, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 6,
    fontsize_col = 6,
    main = sprintf("Gene-Gene Correlation: %s", ct)
  )
  grDevices::dev.off()
  cat(sprintf("[OK] Correlation heatmap: %s\n", pdf_path))
  invisible(cor_mat)
}

# ----- 8.4.7 CoVarNet visualization summary wrapper -----
covarnet_plot_all <- function(covarnet_results,
                              out_dir = OUTPUT_DIR,
                              n_label_top = 20L,
                              n_top_hubs = N_HUB_GENES,
                              max_heatmap_nodes = 80L,
                              plot_raw_correlation = TRUE) {
  cat("\n"); cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat("CoVarNet: Enhanced visualization summary\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  ensure_dir(out_dir)

  res <- list()
  res$summary <- plot_covarnet_summary(covarnet_results, out_dir = out_dir)
  res$edge_direction <- plot_covarnet_edge_direction_summary(covarnet_results, out_dir = out_dir)
  res$degree_distribution <- plot_covarnet_degree_distribution(covarnet_results, out_dir = out_dir)

  for (ct in names(covarnet_results)) {
    ct_res <- covarnet_results[[ct]]
    plot_covarnet_graph(
      graph = ct_res$graph,
      celltype_name = ct,
      n_label_top = n_label_top,
      out_dir = out_dir
    )
    plot_covarnet_hub_rank(
      covarnet_result = ct_res,
      celltype_name = ct,
      n_top = n_top_hubs,
      out_dir = out_dir
    )
    plot_covarnet_adjacency_heatmap(
      graph = ct_res$graph,
      celltype_name = ct,
      max_nodes = max_heatmap_nodes,
      out_dir = out_dir
    )
    if (isTRUE(plot_raw_correlation)) {
      plot_covarnet_correlation_heatmap(
        cor_result = ct_res$cor_result,
        celltype_name = ct,
        max_genes = max_heatmap_nodes,
        out_dir = out_dir
      )
    }
  }

  res$hub_overlap <- plot_hub_overlap(covarnet_results, n_top = n_top_hubs, out_dir = out_dir)
  invisible(res)
}
