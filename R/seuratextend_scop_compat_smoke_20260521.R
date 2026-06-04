#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(Seurat)
  library(scop)
  library(SeuratExtend)
  library(ggplot2)
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
    source_tag <- sprintf("user_rdata:%s::%s", input_rdata, attr(seu, "loaded_name") %||% "unknown")
    attr(seu, "source_tag") <- source_tag
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

pick_groups <- function(seu, preferred = c("B cell", "Mono CD14")) {
  tab <- sort(table(seu$cluster), decreasing = TRUE)
  available <- names(tab)
  if (all(preferred %in% available)) {
    return(preferred)
  }
  if (length(available) < 2) {
    stop("Need at least two groups in `cluster` metadata for comparison")
  }
  available[1:2]
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

main <- function() {
  input_rdata <- Sys.getenv("SEURAT_RDATA", "")
  output_dir <- Sys.getenv("OUTPUT_DIR", "/home/h2048/temp/seuratextend_scop_compat_20260521")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  seu <- load_demo_or_user_object(input_rdata)
  seu <- ensure_cluster_column(seu)
  seu <- ensure_variable_features(seu)
  seu <- ensure_umap(seu)

  groups <- pick_groups(seu)
  g1 <- groups[[1]]
  g2 <- groups[[2]]

  overlap <- sort(intersect(getNamespaceExports("scop"), getNamespaceExports("SeuratExtend")))

  p1 <- SeuratExtend::DimPlot2(
    seu,
    features = intersect(c("orig.ident", "cluster"), colnames(seu@meta.data)),
    cols = "bright",
    ncol = 2,
    theme = Seurat::NoAxes()
  )
  ggsave(file.path(output_dir, "01_seuratextend_dimplot2.png"), p1, width = 10, height = 5, dpi = 180)

  seu_gsa <- SeuratExtend::GeneSetAnalysis(
    seu,
    genesets = hall50$human,
    nCores = 1,
    verbose = FALSE
  )
  matr <- seu_gsa@misc$AUCell$genesets
  write.table(
    data.frame(pathway = rownames(matr)),
    file.path(output_dir, "hall50_pathways.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  p2 <- SeuratExtend::WaterfallPlot(
    matr,
    f = seu_gsa$cluster,
    ident.1 = g1,
    ident.2 = g2,
    style = "segment",
    color_theme = "D"
  )
  ggsave(file.path(output_dir, "02_seuratextend_waterfall_pathway.png"), p2, width = 8, height = 7, dpi = 180)

  p3 <- SeuratExtend::VolcanoPlot(
    seu_gsa,
    ident.1 = g1,
    ident.2 = g2,
    log.base = "2",
    top.n = 5
  )
  ggsave(file.path(output_dir, "03_seuratextend_volcano.png"), p3, width = 7, height = 6, dpi = 180)

  p4 <- scop::CellDimPlot(
    seu_gsa,
    group.by = "cluster",
    reduction = "umap",
    label = TRUE,
    theme_use = "theme_blank"
  )
  ggsave(file.path(output_dir, "04_scop_celldimplot.png"), p4, width = 7, height = 6, dpi = 180)

  seu_de <- scop::RunDEtest(
    seu_gsa,
    group.by = "cluster",
    group1 = g1,
    group2 = g2,
    cores = 1,
    verbose = FALSE
  )
  de_res <- seu_de@tools$DEtest_custom$AllMarkers_wilcox

  p5 <- scop::VolcanoPlot(
    seu_de,
    group.by = "cluster",
    group_use = c(g1, g2),
    res = de_res
  )
  ggsave(file.path(output_dir, "05_scop_volcano.png"), p5, width = 7, height = 6, dpi = 180)

  summary_lines <- c(
    sprintf("source\t%s", attr(seu, "source_tag") %||% "unknown"),
    sprintf("group1\t%s", g1),
    sprintf("group2\t%s", g2),
    sprintf("hall50_dim\t%d x %d", nrow(matr), ncol(matr)),
    sprintf("namespace_overlaps\t%s", paste(overlap, collapse = ",")),
    "note\tIf both packages are attached, prefer scop::VolcanoPlot / SeuratExtend::VolcanoPlot explicitly.",
    "note\tFor scop::VolcanoPlot on custom comparisons, pass res = seu_de@tools$DEtest_custom$AllMarkers_wilcox after scop::RunDEtest()."
  )
  writeLines(summary_lines, file.path(output_dir, "compat_summary.tsv"))

  cat("[OK] Compatibility smoke test complete\n")
  cat(sprintf("[OK] Source: %s\n", attr(seu, "source_tag") %||% "unknown"))
  cat(sprintf("[OK] Groups: %s vs %s\n", g1, g2))
  cat(sprintf("[OK] Output dir: %s\n", output_dir))
  cat(sprintf("[OK] Overlaps: %s\n", paste(overlap, collapse = ", ")))
  cat(sprintf("[OK] Files: %s\n", paste(sort(list.files(output_dir)), collapse = ", ")))
}

main()
