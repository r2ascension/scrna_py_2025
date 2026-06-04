#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(Seurat)
  library(scop)
  library(SeuratExtend)
  library(ggplot2)
  library(cowplot)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.character(x) && !nzchar(x))) y else x
}

pick_first_seurat <- function(env) {
  nms <- ls(env, all.names = TRUE)
  for (nm in nms) {
    obj <- get(nm, envir = env)
    if (inherits(obj, "Seurat")) {
      attr(obj, "loaded_name") <- nm
      return(obj)
    }
  }
  NULL
}

load_demo_or_user_object <- function(input_rdata = "") {
  if (nzchar(input_rdata) && file.exists(input_rdata)) {
    e <- new.env(parent = emptyenv())
    load(input_rdata, envir = e)
    seu <- pick_first_seurat(e)
    if (is.null(seu)) {
      stop(sprintf("No Seurat object found in: %s", input_rdata))
    }
    attr(seu, "source_tag") <- sprintf(
      "user_rdata:%s::%s",
      input_rdata,
      attr(seu, "loaded_name") %||% "unknown"
    )
    return(seu)
  }

  data("pbmc", package = "SeuratExtend")
  attr(pbmc, "source_tag") <- "SeuratExtend::pbmc"
  pbmc
}

ensure_cluster_column <- function(seu) {
  if (!"cluster" %in% colnames(seu@meta.data)) {
    seu$cluster <- as.character(Idents(seu))
  }
  seu$cluster <- as.character(seu$cluster)
  seu
}

ensure_idents_from_cluster <- function(seu) {
  cluster_levels <- names(sort(table(seu$cluster), decreasing = TRUE))
  seu$cluster <- factor(seu$cluster, levels = cluster_levels)
  Idents(seu) <- seu$cluster
  seu
}

ensure_variable_features <- function(seu, nfeatures = 2000) {
  feats <- tryCatch(VariableFeatures(seu), error = function(e) character(0))
  if (length(feats) == 0) {
    seu <- FindVariableFeatures(seu, nfeatures = nfeatures, verbose = FALSE)
  }
  seu
}

ensure_umap <- function(seu) {
  if (!"pca" %in% names(seu@reductions)) {
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, verbose = FALSE)
  }
  if (!"umap" %in% names(seu@reductions)) {
    seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
    seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
  }
  seu
}

pick_groups <- function(seu, group1 = "", group2 = "", preferred = c("B cell", "Mono CD14")) {
  available <- levels(seu$cluster)
  if (nzchar(group1) && nzchar(group2) && all(c(group1, group2) %in% available)) {
    return(c(group1, group2))
  }
  if (all(preferred %in% available)) {
    return(preferred)
  }
  if (length(available) < 2) {
    stop("Need at least two groups in `cluster` metadata for comparison")
  }
  available[1:2]
}

pick_sample_column <- function(seu) {
  meta <- seu@meta.data
  preferred <- c("orig.ident", "sample", "Sample", "batch", "Batch", "dataset", "Dataset")
  for (nm in preferred) {
    if (nm %in% colnames(meta)) {
      x <- meta[[nm]]
      if (length(unique(x[!is.na(x)])) > 1) {
        return(nm)
      }
    }
  }
  for (nm in colnames(meta)) {
    if (nm == "cluster") {
      next
    }
    x <- meta[[nm]]
    nunique <- length(unique(x[!is.na(x)]))
    if (nunique >= 2 && nunique <= 12) {
      return(nm)
    }
  }
  ""
}

short_hallmark <- function(x) {
  gsub("_", " ", sub("^HALLMARK_", "", x))
}

theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 4),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

save_plot_bundle <- function(plot, stem, output_dir, width, height, dpi = 320, bg = "white") {
  png_file <- file.path(output_dir, paste0(stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = bg)
  ggsave(pdf_file, plot, width = width, height = height, bg = bg, device = grDevices::pdf, useDingbats = FALSE)
  c(png = png_file, pdf = pdf_file)
}

build_umap_composition_plot <- function(seu, sample_key = "") {
  p_umap <- SeuratExtend::DimPlot2(
    seu,
    cols = "pro_blue",
    label = TRUE,
    repel = TRUE,
    raster = TRUE,
    theme = Seurat::NoAxes() +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        text = element_text(colour = "black"),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 20)
      )
  ) +
    labs(title = "Cell states")

  if (!nzchar(sample_key)) {
    return(p_umap)
  }

  p_bar <- SeuratExtend::ClusterDistrBar(
    origin = seu@meta.data[[sample_key]],
    cluster = seu$cluster,
    cols = "pro_blue",
    flip = FALSE,
    border = "black"
  ) +
    labs(
      title = sprintf("Composition by %s", sample_key),
      x = NULL,
      y = "Percentage of cell counts"
    ) +
    theme_pub(14) +
    theme(
      plot.title = element_text(size = 16),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

  cowplot::plot_grid(
    p_umap,
    p_bar,
    ncol = 2,
    rel_widths = c(1.15, 0.85)
  )
}

run_pathway_analysis <- function(seu) {
  seu_gsa <- SeuratExtend::GeneSetAnalysis(
    seu,
    genesets = hall50$human,
    nCores = 1,
    verbose = FALSE
  )
  matr <- seu_gsa@misc$AUCell$genesets
  list(seu = seu_gsa, matr = matr)
}

build_pathway_contrast_table <- function(matr, clusters, group1, group2) {
  idx1 <- which(as.character(clusters) == group1)
  idx2 <- which(as.character(clusters) == group2)
  if (length(idx1) == 0 || length(idx2) == 0) {
    stop("One of the requested groups is absent in pathway matrix")
  }
  mean1 <- rowMeans(matr[, idx1, drop = FALSE])
  mean2 <- rowMeans(matr[, idx2, drop = FALSE])
  out <- data.frame(
    pathway = rownames(matr),
    pathway_label = short_hallmark(rownames(matr)),
    mean_group1 = mean1,
    mean_group2 = mean2,
    delta_group1_minus_group2 = mean1 - mean2,
    abs_delta = abs(mean1 - mean2),
    stringsAsFactors = FALSE
  )
  out[order(out$abs_delta, decreasing = TRUE), ]
}

build_pathway_heatmap_plot <- function(matr, clusters, per_cluster = 2) {
  zscore_tbl <- SeuratExtend::CalcStats(matr, f = clusters, order = "p", n = per_cluster)
  rownames(zscore_tbl) <- short_hallmark(rownames(zscore_tbl))
  SeuratExtend::Heatmap(
    zscore_tbl,
    lab_fill = "zscore",
    color_scheme = "A",
    text.size = 8,
    legend_position = "right"
  ) +
    labs(title = "Cluster-level hallmark programs") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    )
}

build_pathway_waterfall_plot <- function(matr, clusters, group1, group2, top_n = 12, len_threshold = 0.15) {
  matr_plot <- matr
  rownames(matr_plot) <- short_hallmark(rownames(matr_plot))
  SeuratExtend::WaterfallPlot(
    matr_plot,
    f = clusters,
    ident.1 = group1,
    ident.2 = group2,
    style = "segment",
    color_theme = "D",
    top.n = top_n,
    len.threshold = len_threshold,
    title = sprintf("Top pathway shifts: %s vs %s", group1, group2)
  ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      legend.title = element_text(face = "bold")
    )
}

safe_marker_table <- function(seu, group1, group2) {
  tryCatch({
    old_ident <- Idents(seu)
    on.exit(Idents(seu) <- old_ident, add = TRUE)
    Idents(seu) <- seu$cluster
    tbl <- Seurat::FindMarkers(
      seu,
      ident.1 = group1,
      ident.2 = group2,
      logfc.threshold = 0,
      min.pct = 0.05,
      verbose = FALSE
    )
    tbl$gene <- rownames(tbl)
    rownames(tbl) <- NULL
    fc_col <- grep("^avg_log.*FC$|^log2FC$", colnames(tbl), value = TRUE)[1]
    p_col <- c("p_val_adj", "p_adj", "pvalue_adj", "p_val")[c("p_val_adj", "p_adj", "pvalue_adj", "p_val") %in% colnames(tbl)][1]
    list(table = tbl, fc_col = fc_col, p_col = p_col)
  }, error = function(e) {
    message("[WARN] FindMarkers export skipped: ", conditionMessage(e))
    NULL
  })
}

select_marker_labels <- function(tbl, fc_col, p_col, n_each = 6) {
  if (is.null(fc_col) || is.null(p_col)) {
    return(tbl[0, , drop = FALSE])
  }
  ord_fun <- function(x) x[order(x[[p_col]], -abs(x[[fc_col]])), , drop = FALSE]
  pos <- ord_fun(tbl[tbl[[fc_col]] > 0, , drop = FALSE])
  neg <- ord_fun(tbl[tbl[[fc_col]] < 0, , drop = FALSE])
  out <- unique(rbind(utils::head(pos, n_each), utils::head(neg, n_each)))
  out[order(out[[p_col]], -abs(out[[fc_col]])), , drop = FALSE]
}

build_seuratextend_volcano_plot <- function(seu, group1, group2, top_n = 8) {
  SeuratExtend::VolcanoPlot(
    seu,
    ident.1 = group1,
    ident.2 = group2,
    x.quantile = 0.98,
    y.quantile = 0.95,
    top.n = top_n,
    log.base = "2",
    color = c("grey85", "steelblue3", "firebrick3")
  ) +
    theme_pub(14)
}

build_scop_volcano <- function(seu, group1, group2) {
  tryCatch({
    seu_de <- scop::RunDEtest(
      seu,
      group.by = "cluster",
      group1 = group1,
      group2 = group2,
      cores = 1,
      verbose = FALSE
    )
    de_res <- seu_de@tools$DEtest_custom$AllMarkers_wilcox
    p <- scop::VolcanoPlot(
      seu_de,
      group.by = "cluster",
      group_use = c(group1, group2),
      res = de_res,
      nlabel = 10,
      pt.size = 1.8,
      pt.alpha = 0.9,
      cols.background = "grey85",
      label.bg = "white",
      label.size = 4.2,
      theme_use = "theme_classic",
      xlab = sprintf("%s <- log2FC -> %s", group2, group1),
      ylab = "-log10(adj. p)"
    )
    list(plot = p, de_res = de_res)
  }, error = function(e) {
    message("[WARN] scop volcano skipped: ", conditionMessage(e))
    NULL
  })
}

main <- function() {
  input_rdata <- Sys.getenv("SEURAT_RDATA", "")
  output_dir <- Sys.getenv("OUTPUT_DIR", "/home/h2048/temp/seuratextend_plot_optimization_20260522")
  group1_env <- Sys.getenv("GROUP_1", "")
  group2_env <- Sys.getenv("GROUP_2", "")
  pathway_top_n <- as.integer(Sys.getenv("PATHWAY_TOP_N", "12"))
  heatmap_per_cluster <- as.integer(Sys.getenv("HEATMAP_PER_CLUSTER", "2"))
  volcano_top_n <- as.integer(Sys.getenv("VOLCANO_TOP_N", "8"))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  figure_contract <- data.frame(
    figure = c(
      "01_umap_composition",
      "02_pathway_heatmap",
      "03_pathway_waterfall",
      "04_seuratextend_volcano",
      "05_scop_volcano"
    ),
    conclusion = c(
      "Embedding should show cluster separation clearly, while the companion bar plot shows whether sample composition is imbalanced.",
      "Cluster-level Hallmark programs should be compared with limited rows and shortened names to reduce cognitive load.",
      "The contrast figure should only show the strongest pathway shifts instead of the full pathway universe.",
      "The main DE figure should focus attention on the most interpretable markers with restrained labels and thresholds.",
      "scop should remain usable as an alternative volcano workflow when explicit res is supplied."
    ),
    stringsAsFactors = FALSE
  )
  write.table(
    figure_contract,
    file.path(output_dir, "figure_contract.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  seu <- load_demo_or_user_object(input_rdata)
  source_tag <- attr(seu, "source_tag") %||% "unknown"
  seu <- ensure_cluster_column(seu)
  seu <- ensure_variable_features(seu)
  seu <- ensure_umap(seu)
  seu <- ensure_idents_from_cluster(seu)

  sample_key <- pick_sample_column(seu)
  groups <- pick_groups(seu, group1 = group1_env, group2 = group2_env)
  group1 <- groups[[1]]
  group2 <- groups[[2]]

  manifest <- list()

  p_umap <- build_umap_composition_plot(seu, sample_key = sample_key)
  manifest[["01_umap_composition"]] <- save_plot_bundle(
    p_umap,
    "01_umap_composition",
    output_dir,
    width = if (nzchar(sample_key)) 11 else 7.2,
    height = if (nzchar(sample_key)) 5.4 else 6.0
  )

  pathway_res <- run_pathway_analysis(seu)
  seu_gsa <- pathway_res$seu
  matr <- pathway_res$matr

  pathway_tbl <- build_pathway_contrast_table(matr, seu_gsa$cluster, group1, group2)
  write.table(
    pathway_tbl,
    file.path(output_dir, "pathway_contrast_table.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    utils::head(pathway_tbl, pathway_top_n),
    file.path(output_dir, "pathway_contrast_top.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  p_heatmap <- build_pathway_heatmap_plot(matr, seu_gsa$cluster, per_cluster = heatmap_per_cluster)
  manifest[["02_pathway_heatmap"]] <- save_plot_bundle(
    p_heatmap,
    "02_pathway_heatmap",
    output_dir,
    width = 7.4,
    height = 6.8
  )

  p_waterfall <- build_pathway_waterfall_plot(
    matr,
    seu_gsa$cluster,
    group1,
    group2,
    top_n = pathway_top_n,
    len_threshold = 0.15
  )
  manifest[["03_pathway_waterfall"]] <- save_plot_bundle(
    p_waterfall,
    "03_pathway_waterfall",
    output_dir,
    width = 7.2,
    height = 5.8
  )

  p_se <- build_seuratextend_volcano_plot(seu_gsa, group1, group2, top_n = volcano_top_n)
  manifest[["04_seuratextend_volcano"]] <- save_plot_bundle(
    p_se,
    "04_seuratextend_volcano",
    output_dir,
    width = 7.2,
    height = 6.2
  )

  marker_res <- safe_marker_table(seu_gsa, group1, group2)
  if (!is.null(marker_res)) {
    write.table(
      marker_res$table,
      file.path(output_dir, "markers_seurat_findmarkers.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    marker_labels <- select_marker_labels(marker_res$table, marker_res$fc_col, marker_res$p_col, n_each = ceiling(volcano_top_n / 2))
    write.table(
      marker_labels,
      file.path(output_dir, "markers_seurat_labels.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  scop_res <- build_scop_volcano(seu_gsa, group1, group2)
  if (!is.null(scop_res)) {
    manifest[["05_scop_volcano"]] <- save_plot_bundle(
      scop_res$plot,
      "05_scop_volcano",
      output_dir,
      width = 7.2,
      height = 6.2
    )
    write.table(
      data.frame(scop_res$de_res),
      file.path(output_dir, "markers_scop_detest.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  overlap <- sort(intersect(getNamespaceExports("scop"), getNamespaceExports("SeuratExtend")))

  manifest_tbl <- do.call(
    rbind,
    lapply(names(manifest), function(nm) {
      data.frame(
        figure = nm,
        png = unname(manifest[[nm]][["png"]]),
        pdf = unname(manifest[[nm]][["pdf"]]),
        stringsAsFactors = FALSE
      )
    })
  )
  write.table(
    manifest_tbl,
    file.path(output_dir, "figure_manifest.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  summary_lines <- c(
    sprintf("source\t%s", source_tag),
    sprintf("group1\t%s", group1),
    sprintf("group2\t%s", group2),
    sprintf("sample_key\t%s", if (nzchar(sample_key)) sample_key else "<none>"),
    sprintf("hall50_dim\t%d x %d", nrow(matr), ncol(matr)),
    sprintf("namespace_overlaps\t%s", paste(overlap, collapse = ",")),
    "note\tOptimization strategy: white background, direct labels, reduced clutter, shortened pathway names, PNG/PDF dual export.",
    "note\tIf both packages are attached, prefer explicit namespaces for VolcanoPlot and RunSlingshot-related calls.",
    "note\tFor scop::VolcanoPlot on custom comparisons, pass res = seu_de@tools$DEtest_custom$AllMarkers_wilcox after scop::RunDEtest()."
  )
  writeLines(summary_lines, file.path(output_dir, "optimization_summary.tsv"))

  cat("[OK] Plot optimization suite complete\n")
  cat(sprintf("[OK] Source: %s\n", source_tag))
  cat(sprintf("[OK] Groups: %s vs %s\n", group1, group2))
  cat(sprintf("[OK] Sample key: %s\n", if (nzchar(sample_key)) sample_key else "<none>"))
  cat(sprintf("[OK] Output dir: %s\n", output_dir))
  cat(sprintf("[OK] Figures: %s\n", paste(manifest_tbl$figure, collapse = ", ")))
}

main()
