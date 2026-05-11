#!/usr/bin/env Rscript
# ==============================================================================
# smc_anno first-tier expression-sufficient trajectory methods (2026-05-08)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

TRADESEQ_HELPER_PATH <- "/home/h2048/script/R/smc_anno_tradeseq_helper_20260429.R"
TRAJECTORY_INTERPRET_HELPER_PATH <- "/home/h2048/script/R/stromal_smc_trajectory_interpret_helper_20260429.R"
VISUAL_LLM_HELPER_PATH <- "/home/h2048/script/R/smc_anno_visual_llm_helper_20260506.R"
if (!exists("smcanno_ts_select_genes", mode = "function")) source(TRADESEQ_HELPER_PATH)
if (!exists("smcti_write_gene_pseudotime_trend_bundle", mode = "function")) source(TRAJECTORY_INTERPRET_HELPER_PATH)
if (!exists("smcanno_viz_register_figure", mode = "function")) source(VISUAL_LLM_HELPER_PATH)

if (!exists("INPUT_RDATA")) INPUT_RDATA <- "/home/h2048/temp/smc_anno.Rdata"
if (!exists("OBJECT_NAME")) OBJECT_NAME <- "smc_anno"
if (!exists("BASE_OUTPUT_DIR")) BASE_OUTPUT_DIR <- "/home/h2048/temp/smc_anno_slingshot_20260428"
if (!exists("OUTPUT_DIR")) OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, "first_tier_methods")
if (!exists("ASSAY_NAME")) ASSAY_NAME <- "RNA"
if (!exists("CLUSTER_COL")) CLUSTER_COL <- "celltype"
if (!exists("PHENOTYPE_COL")) PHENOTYPE_COL <- "group"
if (!exists("ROOT_CLUSTER")) ROOT_CLUSTER <- "Pericyte"
if (!exists("PYTHON_BIN")) PYTHON_BIN <- "/home/h2048/miniconda3/envs/bbknn_env/bin/python"
if (!exists("PYTHON_RUNNER")) PYTHON_RUNNER <- "/home/h2048/script/py/smc_anno_first_tier_python_methods_20260508.py"
if (!exists("MAX_H5AD_GENES")) MAX_H5AD_GENES <- 2000L
if (!exists("GENE_TREND_MAX_GENES")) GENE_TREND_MAX_GENES <- 16L
if (!exists("RUN_CYTOTRACE2")) RUN_CYTOTRACE2 <- TRUE
if (!exists("CYTOTRACE2_NCORES")) CYTOTRACE2_NCORES <- 4L

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

safe_numeric_summary <- function(x, fun) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, na.rm = TRUE, type = 8))
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))

pick_seurat_object <- function(env, object_name = NULL) {
  env_names <- ls(env, all.names = TRUE)
  if (!is.null(object_name) && nzchar(object_name) && object_name %in% env_names) {
    obj <- get(object_name, envir = env)
    if (inherits(obj, "Seurat")) return(obj)
  }
  seurat_names <- env_names[vapply(env_names, function(nm) inherits(get(nm, envir = env), "Seurat"), logical(1))]
  if (length(seurat_names) != 1L) stop("Could not uniquely identify Seurat object in INPUT_RDATA", call. = FALSE)
  get(seurat_names[[1]], envir = env)
}

orient_root_low <- function(pt, root_mask) {
  pt <- suppressWarnings(as.numeric(pt))
  finite <- is.finite(pt)
  out <- rep(NA_real_, length(pt))
  if (!any(finite)) return(out)
  rng <- range(pt[finite], finite = TRUE)
  if (diff(rng) > 0) out[finite] <- (pt[finite] - rng[1]) / diff(rng) else out[finite] <- 0
  if (any(root_mask, na.rm = TRUE)) {
    root_mean <- mean(out[root_mask & finite], na.rm = TRUE)
    other_mean <- mean(out[!root_mask & finite], na.rm = TRUE)
    if (is.finite(root_mean) && is.finite(other_mean) && root_mean > other_mean) out[finite] <- 1 - out[finite]
  }
  out
}

plot_umap_pseudotime <- function(obj, pt_tbl, pseudotime_col, output_path, title_text, reduction_name = "umap") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!reduction_name %in% names(obj@reductions)) return(invisible(NULL))
  emb <- Seurat::Embeddings(obj, reduction = reduction_name)[pt_tbl$cell_id, 1:2, drop = FALSE]
  df <- data.frame(
    cell_id = pt_tbl$cell_id,
    dim1 = as.numeric(emb[, 1]),
    dim2 = as.numeric(emb[, 2]),
    pseudotime = suppressWarnings(as.numeric(pt_tbl[[pseudotime_col]])),
    cluster_label = if (CLUSTER_COL %in% colnames(pt_tbl)) as.character(pt_tbl[[CLUSTER_COL]]) else NA_character_,
    phenotype = if (PHENOTYPE_COL %in% colnames(pt_tbl)) as.character(pt_tbl[[PHENOTYPE_COL]]) else NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = pseudotime)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.85) +
    ggplot2::scale_colour_viridis_c(option = "viridis", na.value = "grey85") +
    ggplot2::coord_equal() +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(title = title_text, x = "UMAP_1", y = "UMAP_2", colour = "pseudotime") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  ggplot2::ggsave(output_path, p, width = 8, height = 6.5)
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 8, height = 6.5, dpi = 180)
  invisible(df)
}

plot_method_status <- function(inventory, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(inventory) == 0L) return(invisible(NULL))
  inventory$method <- factor(inventory$method, levels = rev(inventory$method))
  p <- ggplot2::ggplot(inventory, ggplot2::aes(x = status, y = method, fill = status)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(title = "First-tier expression-sufficient trajectory method status", x = "Status", y = "Method") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")
  ggplot2::ggsave(output_path, p, width = 9, height = max(5.5, 0.35 * nrow(inventory) + 2))
  ggplot2::ggsave(sub("\\.pdf$", ".png", output_path), p, width = 9, height = max(5.5, 0.35 * nrow(inventory) + 2), dpi = 180)
  invisible(output_path)
}

export_expression_h5ad <- function(obj, out_dir, genes, python_bin) {
  export_dir <- ensure_dir(file.path(out_dir, "expression_h5ad_export"))
  cells <- colnames(obj)
  genes <- unique(genes[genes %in% rownames(obj)])
  counts <- smcanno_ts_get_assay_counts(obj, ASSAY_NAME)
  data_mat <- tryCatch(smcti_get_assay_matrix(obj, ASSAY_NAME, slot_name = "data"), error = function(e) NULL)
  if (is.null(data_mat) || !all(genes %in% rownames(data_mat))) {
    data_mat <- counts
    data_mat@x <- log1p(data_mat@x)
  }
  counts_sub <- counts[genes, cells, drop = FALSE]
  data_sub <- data_mat[genes, cells, drop = FALSE]

  obs_cols <- intersect(unique(c(CLUSTER_COL, PHENOTYPE_COL, "sample", "orig.ident", "seurat_clusters", "group")), colnames(obj@meta.data))
  obs <- data.frame(cell_id = cells, obj@meta.data[cells, obs_cols, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in colnames(obs)) if (is.factor(obs[[nm]])) obs[[nm]] <- as.character(obs[[nm]])
  var <- data.frame(gene = genes, stringsAsFactors = FALSE, check.names = FALSE)
  obs_path <- file.path(export_dir, "obs.tsv")
  var_path <- file.path(export_dir, "var.tsv")
  counts_path <- file.path(export_dir, "counts_cells_by_genes.mtx")
  data_path <- file.path(export_dir, "data_cells_by_genes.mtx")
  write_tsv(obs, obs_path)
  write_tsv(var, var_path)
  Matrix::writeMM(Matrix::t(counts_sub), counts_path)
  Matrix::writeMM(Matrix::t(data_sub), data_path)

  emb_paths <- list()
  for (red in intersect(c("umap", "pca", "harmony", "scanvi", "scvi"), names(obj@reductions))) {
    emb <- Seurat::Embeddings(obj, reduction = red)[cells, , drop = FALSE]
    emb_df <- data.frame(cell_id = rownames(emb), emb, stringsAsFactors = FALSE, check.names = FALSE)
    emb_path <- file.path(export_dir, sprintf("embedding_%s.tsv", sanitize_name(red)))
    write_tsv(emb_df, emb_path)
    emb_paths[[red]] <- emb_path
  }
  h5ad_path <- file.path(export_dir, "smc_anno_expression_hvg.h5ad")
  builder_path <- file.path(export_dir, "build_expression_h5ad.py")
  emb_lines <- unlist(lapply(names(emb_paths), function(red) {
    key <- paste0("X_", sanitize_name(red))
    c(
      sprintf("emb = pd.read_csv(r'%s', sep='\\t').set_index('cell_id')", emb_paths[[red]]),
      "emb = emb.loc[obs.index]",
      sprintf("adata.obsm['%s'] = emb.values.astype('float32')", key)
    )
  }))
  writeLines(c(
    "#!/usr/bin/env python3",
    "import pandas as pd",
    "import anndata as ad",
    "from scipy.io import mmread",
    sprintf("obs = pd.read_csv(r'%s', sep='\\t').set_index('cell_id')", obs_path),
    sprintf("var = pd.read_csv(r'%s', sep='\\t').set_index('gene')", var_path),
    sprintf("X = mmread(r'%s').tocsr().astype('float32')", data_path),
    sprintf("counts = mmread(r'%s').tocsr().astype('float32')", counts_path),
    "for col in obs.columns:",
    "    if obs[col].dtype == object:",
    "        obs[col] = obs[col].fillna('').astype(str)",
    "adata = ad.AnnData(X=X, obs=obs, var=var)",
    "adata.layers['counts'] = counts",
    emb_lines,
    sprintf("adata.write_h5ad(r'%s', compression='gzip')", h5ad_path)
  ), con = builder_path, useBytes = TRUE)
  Sys.chmod(builder_path, "0755")
  status <- system2(python_bin, args = builder_path)
  if (!identical(status, 0L) || !file.exists(h5ad_path)) stop("Failed to export expression h5ad", call. = FALSE)
  list(h5ad_path = h5ad_path, genes = genes, obs_path = obs_path, var_path = var_path, export_dir = export_dir)
}

register_method_figure <- function(method, figure_path, data_path, figure_type, title, manifest_path, extra_context = character()) {
  if (!file.exists(figure_path) || !file.exists(data_path)) return(NULL)
  df <- utils::read.delim(data_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  smcanno_viz_register_figure(
    figure_path = figure_path,
    data = df,
    method = method,
    figure_type = figure_type,
    title = title,
    extra_context = extra_context,
    manifest_path = manifest_path
  )
}

register_python_outputs <- function(py_root, manifest_path) {
  specs <- list(
    DPT = list(fig = "DPT/DPT_umap_pseudotime.pdf", data = "DPT/DPT_pseudotime.tsv", type = "dpt_pseudotime_embedding", title = "DPT pseudotime on UMAP"),
    PHATE = list(fig = "PHATE/PHATE_embedding_pseudotime.pdf", data = "PHATE/PHATE_pseudotime.tsv", type = "phate_embedding_root_order", title = "PHATE embedding ordered from root"),
    Palantir = list(fig = "Palantir/Palantir_umap_pseudotime.pdf", data = "Palantir/Palantir_pseudotime.tsv", type = "palantir_pseudotime_embedding", title = "Palantir pseudotime on UMAP"),
    VIA = list(fig = "VIA/VIA_umap_pseudotime.pdf", data = "VIA/VIA_pseudotime.tsv", type = "via_pseudotime_embedding", title = "VIA pseudotime on UMAP"),
    CellRank = list(fig = "CellRank/CellRank_kernel_summary.pdf", data = "CellRank/CellRank_kernel_summary.tsv", type = "cellrank_kernel_transition_summary", title = "CellRank PseudotimeKernel and ConnectivityKernel summary")
  )
  for (method in names(specs)) {
    sp <- specs[[method]]
    register_method_figure(
      method = method,
      figure_path = file.path(py_root, sp$fig),
      data_path = file.path(py_root, sp$data),
      figure_type = sp$type,
      title = sp$title,
      manifest_path = manifest_path,
      extra_context = c("Supplemental first-tier expression-sufficient trajectory run on smc_anno.Rdata.", sprintf("Root cluster: %s; phenotype column: %s", ROOT_CLUSTER, PHENOTYPE_COL))
    )
  }
}

read_existing_method_result <- function(method, out_root) {
  out_dir <- file.path(out_root, method)
  pt_path <- file.path(out_dir, sprintf("%s_pseudotime.tsv", method))
  if (!file.exists(pt_path)) return(NULL)
  pt_tbl <- utils::read.delim(pt_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  list(status = "ok", output_dir = out_dir, pseudotime_table = pt_tbl, n_cells = sum(is.finite(suppressWarnings(as.numeric(pt_tbl$pseudotime)))), reused = TRUE)
}

run_tscan <- function(obj, expr_mat, out_root, trend_genes) {
  method <- "TSCAN"
  out_dir <- ensure_dir(file.path(out_root, method))
  status <- tryCatch({
    if (!requireNamespace("TSCAN", quietly = TRUE)) stop("TSCAN package is not installed", call. = FALSE)
    mat <- as.matrix(expr_mat)
    mclust <- TSCAN::exprmclust(mat, clusternum = 2:7, reduce = TRUE)
    ord <- TSCAN::TSCANorder(mclust, orderonly = FALSE)
    pt_tbl <- data.frame(cell_id = colnames(obj), pseudotime = NA_real_, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(pt_tbl) <- pt_tbl$cell_id
    ord$sample_name <- as.character(ord$sample_name)
    pt_tbl[ord$sample_name, "pseudotime"] <- suppressWarnings(as.numeric(ord$Pseudotime))
    pt_tbl$pseudotime <- orient_root_low(pt_tbl$pseudotime, as.character(obj@meta.data[pt_tbl$cell_id, CLUSTER_COL]) == ROOT_CLUSTER)
    pt_tbl[[CLUSTER_COL]] <- as.character(obj@meta.data[pt_tbl$cell_id, CLUSTER_COL])
    pt_tbl[[PHENOTYPE_COL]] <- as.character(obj@meta.data[pt_tbl$cell_id, PHENOTYPE_COL])
    if ("sample" %in% colnames(obj@meta.data)) pt_tbl$sample <- as.character(obj@meta.data[pt_tbl$cell_id, "sample"])
    write_tsv(pt_tbl, file.path(out_dir, "TSCAN_pseudotime.tsv"))
    plot_df <- plot_umap_pseudotime(obj, pt_tbl, "pseudotime", file.path(out_dir, "TSCAN_umap_pseudotime.pdf"), "TSCAN pseudotime")
    if (!is.null(plot_df)) write_tsv(plot_df, file.path(out_dir, "TSCAN_umap_pseudotime_plot_data.tsv"))
    saveRDS(mclust, file.path(out_dir, "TSCAN_mclust_mst.rds"))
    smcti_write_gene_pseudotime_trend_bundle(obj, pt_tbl, file.path(out_dir, "gene_trends"), ASSAY_NAME, method, trend_genes, cluster_col = CLUSTER_COL, phenotype_col = PHENOTYPE_COL, max_genes = GENE_TREND_MAX_GENES)
    list(status = "ok", output_dir = out_dir, pseudotime_table = pt_tbl, n_cells = sum(is.finite(pt_tbl$pseudotime)))
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(out_dir, "TSCAN_error.txt"), useBytes = TRUE)
    list(status = "error", output_dir = out_dir, error = conditionMessage(e))
  })
  status_json <- status
  status_json$pseudotime_table <- NULL
  smcanno_viz_write_json(status_json, file.path(out_dir, "TSCAN_status.json"))
  status
}

run_scorpius <- function(obj, expr_mat, out_root, trend_genes) {
  method <- "SCORPIUS"
  out_dir <- ensure_dir(file.path(out_root, method))
  status <- tryCatch({
    if (!requireNamespace("SCORPIUS", quietly = TRUE)) stop("SCORPIUS package is not installed", call. = FALSE)
    mat <- t(as.matrix(expr_mat))
    space <- SCORPIUS::reduce_dimensionality(mat, dist = "spearman", ndim = 3)
    traj <- SCORPIUS::infer_trajectory(space)
    pt <- orient_root_low(traj$time, as.character(obj@meta.data[colnames(obj), CLUSTER_COL]) == ROOT_CLUSTER)
    pt_tbl <- data.frame(cell_id = colnames(obj), pseudotime = pt, stringsAsFactors = FALSE, check.names = FALSE)
    pt_tbl[[CLUSTER_COL]] <- as.character(obj@meta.data[pt_tbl$cell_id, CLUSTER_COL])
    pt_tbl[[PHENOTYPE_COL]] <- as.character(obj@meta.data[pt_tbl$cell_id, PHENOTYPE_COL])
    if ("sample" %in% colnames(obj@meta.data)) pt_tbl$sample <- as.character(obj@meta.data[pt_tbl$cell_id, "sample"])
    write_tsv(pt_tbl, file.path(out_dir, "SCORPIUS_pseudotime.tsv"))
    space_df <- data.frame(cell_id = colnames(obj), space1 = space[, 1], space2 = space[, 2], space3 = if (ncol(space) >= 3) space[, 3] else NA_real_, pseudotime = pt, stringsAsFactors = FALSE, check.names = FALSE)
    write_tsv(space_df, file.path(out_dir, "SCORPIUS_space.tsv"))
    plot_df <- plot_umap_pseudotime(obj, pt_tbl, "pseudotime", file.path(out_dir, "SCORPIUS_umap_pseudotime.pdf"), "SCORPIUS pseudotime")
    if (!is.null(plot_df)) write_tsv(plot_df, file.path(out_dir, "SCORPIUS_umap_pseudotime_plot_data.tsv"))
    saveRDS(traj, file.path(out_dir, "SCORPIUS_trajectory.rds"))
    smcti_write_gene_pseudotime_trend_bundle(obj, pt_tbl, file.path(out_dir, "gene_trends"), ASSAY_NAME, method, trend_genes, cluster_col = CLUSTER_COL, phenotype_col = PHENOTYPE_COL, max_genes = GENE_TREND_MAX_GENES)
    list(status = "ok", output_dir = out_dir, pseudotime_table = pt_tbl, n_cells = sum(is.finite(pt_tbl$pseudotime)))
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(out_dir, "SCORPIUS_error.txt"), useBytes = TRUE)
    list(status = "error", output_dir = out_dir, error = conditionMessage(e))
  })
  status_json <- status
  status_json$pseudotime_table <- NULL
  smcanno_viz_write_json(status_json, file.path(out_dir, "SCORPIUS_status.json"))
  status
}

run_cytotrace2 <- function(obj, out_root, trend_genes) {
  method <- "CytoTRACE2"
  out_dir <- ensure_dir(file.path(out_root, method))
  status <- tryCatch({
    if (!requireNamespace("CytoTRACE2", quietly = TRUE)) stop("CytoTRACE2 package is not installed", call. = FALSE)
    cyto_obj <- CytoTRACE2::cytotrace2(obj, species = "human", is_seurat = TRUE, slot_type = "counts", ncores = CYTOTRACE2_NCORES, parallelize_models = FALSE, parallelize_smoothing = FALSE)
    score <- suppressWarnings(as.numeric(cyto_obj@meta.data[colnames(obj), "CytoTRACE2_Score"]))
    rel <- suppressWarnings(as.numeric(cyto_obj@meta.data[colnames(obj), "CytoTRACE2_Relative"]))
    pt <- orient_root_low(1 - score, as.character(obj@meta.data[colnames(obj), CLUSTER_COL]) == ROOT_CLUSTER)
    pt_tbl <- data.frame(cell_id = colnames(obj), pseudotime = pt, CytoTRACE2_Score = score, CytoTRACE2_Relative = rel, CytoTRACE2_Potency = as.character(cyto_obj@meta.data[colnames(obj), "CytoTRACE2_Potency"]), stringsAsFactors = FALSE, check.names = FALSE)
    pt_tbl[[CLUSTER_COL]] <- as.character(obj@meta.data[pt_tbl$cell_id, CLUSTER_COL])
    pt_tbl[[PHENOTYPE_COL]] <- as.character(obj@meta.data[pt_tbl$cell_id, PHENOTYPE_COL])
    if ("sample" %in% colnames(obj@meta.data)) pt_tbl$sample <- as.character(obj@meta.data[pt_tbl$cell_id, "sample"])
    write_tsv(pt_tbl, file.path(out_dir, "CytoTRACE2_pseudotime.tsv"))
    plot_df <- plot_umap_pseudotime(obj, pt_tbl, "pseudotime", file.path(out_dir, "CytoTRACE2_umap_pseudotime.pdf"), "CytoTRACE2 differentiation-oriented pseudotime")
    if (!is.null(plot_df)) write_tsv(plot_df, file.path(out_dir, "CytoTRACE2_umap_pseudotime_plot_data.tsv"))
    potency_summary <- aggregate(CytoTRACE2_Score ~ celltype + group + CytoTRACE2_Potency, data = transform(pt_tbl, celltype = pt_tbl[[CLUSTER_COL]], group = pt_tbl[[PHENOTYPE_COL]]), FUN = function(x) round(mean(x, na.rm = TRUE), 4))
    write_tsv(potency_summary, file.path(out_dir, "CytoTRACE2_score_summary_by_celltype_group.tsv"))
    smcti_write_gene_pseudotime_trend_bundle(obj, pt_tbl, file.path(out_dir, "gene_trends"), ASSAY_NAME, method, trend_genes, cluster_col = CLUSTER_COL, phenotype_col = PHENOTYPE_COL, max_genes = GENE_TREND_MAX_GENES)
    list(status = "ok", output_dir = out_dir, pseudotime_table = pt_tbl, n_cells = sum(is.finite(pt_tbl$pseudotime)))
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(out_dir, "CytoTRACE2_error.txt"), useBytes = TRUE)
    list(status = "error", output_dir = out_dir, error = conditionMessage(e))
  })
  status_json <- status
  status_json$pseudotime_table <- NULL
  smcanno_viz_write_json(status_json, file.path(out_dir, "CytoTRACE2_status.json"))
  status
}

collect_python_method_tables <- function(py_root) {
  methods <- c("DPT", "PHATE", "Palantir", "VIA")
  out <- list()
  for (m in methods) {
    f <- file.path(py_root, m, sprintf("%s_pseudotime.tsv", m))
    if (file.exists(f)) out[[m]] <- utils::read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
  out
}

cat("\n", strrep("=", 80), "\n", sep = "")
cat("SMC_ANNO FIRST-TIER EXPRESSION-SUFFICIENT TRAJECTORY METHODS\n")
cat(strrep("=", 80), "\n", sep = "")
cat("Input : ", INPUT_RDATA, "\n", sep = "")
cat("Output: ", OUTPUT_DIR, "\n", sep = "")
cat("Python: ", PYTHON_BIN, "\n\n", sep = "")

OUTPUT_DIR <- ensure_dir(OUTPUT_DIR)
MANIFEST <- file.path(OUTPUT_DIR, "first_tier_figure_companion_manifest.tsv")
if (file.exists(MANIFEST)) unlink(MANIFEST)

load_env <- new.env(parent = emptyenv())
load(INPUT_RDATA, envir = load_env)
obj <- pick_seurat_object(load_env, OBJECT_NAME)
if (!ASSAY_NAME %in% names(obj@assays)) ASSAY_NAME <- Seurat::DefaultAssay(obj)
if (!PHENOTYPE_COL %in% colnames(obj@meta.data)) {
  detected <- smcti_infer_phenotype_col(obj@meta.data)
  PHENOTYPE_COL <- detected$column
}
if (is.na(PHENOTYPE_COL) || !nzchar(PHENOTYPE_COL)) stop("Disease phenotype column could not be detected", call. = FALSE)

counts <- smcanno_ts_get_assay_counts(obj, ASSAY_NAME)
gene_tbl <- smcanno_ts_select_genes(obj, counts, max_genes = MAX_H5AD_GENES, min_cells_expressed = 20L, exclude_noninformative = TRUE)
trend_panel_path <- file.path(BASE_OUTPUT_DIR, "gene_pseudotime_trends", "gene_trend_panel.tsv")
trend_genes <- if (file.exists(trend_panel_path)) utils::read.delim(trend_panel_path, sep = "\t", stringsAsFactors = FALSE)$gene else smcanno_ts_default_marker_panel
trend_genes <- smcti_filter_gene_panel_for_trends(unique(c(trend_genes, head(gene_tbl$gene, 40L))), obj, max_genes = GENE_TREND_MAX_GENES, exclude_noninformative = TRUE)
write_tsv(gene_tbl, file.path(OUTPUT_DIR, "first_tier_expression_gene_panel.tsv"))
write_tsv(data.frame(gene = trend_genes, stringsAsFactors = FALSE), file.path(OUTPUT_DIR, "first_tier_gene_trend_panel.tsv"))

export_info <- export_expression_h5ad(obj, OUTPUT_DIR, genes = gene_tbl$gene, python_bin = PYTHON_BIN)
cat("Expression h5ad: ", export_info$h5ad_path, "\n", sep = "")

python_out <- ensure_dir(file.path(OUTPUT_DIR, "python_methods"))
python_status_path <- file.path(python_out, "python_first_tier_method_status.tsv")
if (file.exists(python_status_path)) {
  cat("Reusing existing Python first-tier outputs: ", python_out, "\n", sep = "")
  py_status <- 0L
} else {
  py_status <- system2(
    PYTHON_BIN,
    args = c(PYTHON_RUNNER, "--h5ad", export_info$h5ad_path, "--output-dir", python_out, "--cluster-col", CLUSTER_COL, "--phenotype-col", PHENOTYPE_COL, "--root-cluster", ROOT_CLUSTER)
  )
}
cat("Python methods exit status: ", py_status, "\n", sep = "")
register_python_outputs(python_out, MANIFEST)

python_tables <- collect_python_method_tables(python_out)
for (m in names(python_tables)) {
  smcti_write_gene_pseudotime_trend_bundle(
    seurat_obj = obj,
    pseudotime_table = python_tables[[m]],
    output_dir = ensure_dir(file.path(python_out, m, "gene_trends")),
    assay_name = ASSAY_NAME,
    method_name = m,
    gene_panel = trend_genes,
    cluster_col = CLUSTER_COL,
    phenotype_col = PHENOTYPE_COL,
    max_genes = GENE_TREND_MAX_GENES
  )
}

expr_data <- smcti_get_assay_matrix(obj, ASSAY_NAME, slot_name = "data")
expr_genes <- head(gene_tbl$gene[gene_tbl$gene %in% rownames(expr_data)], 1000L)
expr_mat <- expr_data[expr_genes, colnames(obj), drop = FALSE]

tscan_res <- read_existing_method_result("TSCAN", OUTPUT_DIR)
if (is.null(tscan_res)) tscan_res <- run_tscan(obj, expr_mat, OUTPUT_DIR, trend_genes)
scorpius_res <- read_existing_method_result("SCORPIUS", OUTPUT_DIR)
if (is.null(scorpius_res)) scorpius_res <- run_scorpius(obj, expr_mat, OUTPUT_DIR, trend_genes)
cyto_res <- read_existing_method_result("CytoTRACE2", OUTPUT_DIR)
if (is.null(cyto_res)) cyto_res <- if (isTRUE(RUN_CYTOTRACE2)) run_cytotrace2(obj, OUTPUT_DIR, trend_genes) else list(status = "disabled", output_dir = file.path(OUTPUT_DIR, "CytoTRACE2"))

for (spec in list(
  list(method = "TSCAN", fig = "TSCAN/TSCAN_umap_pseudotime.pdf", data = "TSCAN/TSCAN_pseudotime.tsv", title = "TSCAN pseudotime on UMAP"),
  list(method = "SCORPIUS", fig = "SCORPIUS/SCORPIUS_umap_pseudotime.pdf", data = "SCORPIUS/SCORPIUS_pseudotime.tsv", title = "SCORPIUS pseudotime on UMAP"),
  list(method = "CytoTRACE2", fig = "CytoTRACE2/CytoTRACE2_umap_pseudotime.pdf", data = "CytoTRACE2/CytoTRACE2_pseudotime.tsv", title = "CytoTRACE2 differentiation-oriented pseudotime on UMAP")
)) {
  register_method_figure(spec$method, file.path(OUTPUT_DIR, spec$fig), file.path(OUTPUT_DIR, spec$data), "supplemental_pseudotime_embedding", spec$title, MANIFEST, extra_context = c("Supplemental first-tier expression-sufficient trajectory run on smc_anno.Rdata.", sprintf("Root cluster: %s; phenotype column: %s", ROOT_CLUSTER, PHENOTYPE_COL)))
}

method_tables <- list()
existing_method_paths <- c(Slingshot = file.path(BASE_OUTPUT_DIR, "slingshot_pseudotime.tsv"), Monocle3 = file.path(BASE_OUTPUT_DIR, "monocle3", "monocle3_pseudotime.tsv"), Monocle2 = file.path(BASE_OUTPUT_DIR, "monocle2", "monocle2_pseudotime.tsv"))
for (method_name in names(existing_method_paths)) {
  p <- existing_method_paths[[method_name]]
  if (file.exists(p)) method_tables[[method_name]] <- utils::read.delim(p, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
method_tables <- c(method_tables, python_tables)
if (identical(tscan_res$status, "ok")) method_tables$TSCAN <- tscan_res$pseudotime_table
if (identical(scorpius_res$status, "ok")) method_tables$SCORPIUS <- scorpius_res$pseudotime_table
if (identical(cyto_res$status, "ok")) method_tables$CytoTRACE2 <- cyto_res$pseudotime_table

first_tier_pheno <- smcti_write_phenotype_pseudotime_comparison_bundle(
  method_tables = method_tables,
  seurat_obj = obj,
  output_dir = ensure_dir(file.path(OUTPUT_DIR, "phenotype_pseudotime_comparison")),
  phenotype_col = PHENOTYPE_COL,
  cluster_col = CLUSTER_COL
)

status_rows <- list(
  data.frame(method = "Monocle3", tier = "expression_sufficient", status = if (file.exists(file.path(BASE_OUTPUT_DIR, "monocle3", "monocle3_pseudotime.tsv"))) "already_completed" else "missing", output_dir = file.path(BASE_OUTPUT_DIR, "monocle3"), notes = "Main run output", stringsAsFactors = FALSE),
  data.frame(method = "Slingshot", tier = "expression_sufficient", status = if (file.exists(file.path(BASE_OUTPUT_DIR, "slingshot_pseudotime.tsv"))) "already_completed" else "missing", output_dir = BASE_OUTPUT_DIR, notes = "Main run output", stringsAsFactors = FALSE),
  data.frame(method = "PAGA", tier = "expression_sufficient_topology", status = if (file.exists(file.path(BASE_OUTPUT_DIR, "paga", "paga_edge_table.tsv"))) "already_completed" else "missing", output_dir = file.path(BASE_OUTPUT_DIR, "paga"), notes = "Topology method; no single-cell pseudotime", stringsAsFactors = FALSE),
  data.frame(method = "DPT", tier = "expression_sufficient", status = if ("DPT" %in% names(python_tables)) "ok" else "error_or_missing", output_dir = file.path(python_out, "DPT"), notes = "Scanpy DPT", stringsAsFactors = FALSE),
  data.frame(method = "Palantir", tier = "expression_sufficient", status = if ("Palantir" %in% names(python_tables)) "ok" else "error_or_missing", output_dir = file.path(python_out, "Palantir"), notes = "Python Palantir", stringsAsFactors = FALSE),
  data.frame(method = "PHATE", tier = "expression_sufficient_embedding", status = if ("PHATE" %in% names(python_tables)) "ok" else "error_or_missing", output_dir = file.path(python_out, "PHATE"), notes = "PHATE embedding; PHATE1 root-oriented ordering", stringsAsFactors = FALSE),
  data.frame(method = "VIA", tier = "expression_sufficient", status = if ("VIA" %in% names(python_tables)) "ok" else "error_or_missing", output_dir = file.path(python_out, "VIA"), notes = "pyVIA", stringsAsFactors = FALSE),
  data.frame(method = "CellRank PseudotimeKernel/ConnectivityKernel", tier = "expression_sufficient_kernel", status = if (file.exists(file.path(python_out, "CellRank", "CellRank_kernel_summary.tsv"))) "ok" else "error_or_missing", output_dir = file.path(python_out, "CellRank"), notes = "Kernel computation uses DPT pseudotime and neighbor connectivity", stringsAsFactors = FALSE),
  data.frame(method = "TSCAN", tier = "expression_sufficient", status = tscan_res$status, output_dir = tscan_res$output_dir, notes = if (!is.null(tscan_res$error)) tscan_res$error else "R TSCAN MST ordering", stringsAsFactors = FALSE),
  data.frame(method = "SCORPIUS", tier = "expression_sufficient", status = scorpius_res$status, output_dir = scorpius_res$output_dir, notes = if (!is.null(scorpius_res$error)) scorpius_res$error else "R SCORPIUS principal curve", stringsAsFactors = FALSE),
  data.frame(method = "CytoTRACE2", tier = "expression_sufficient_maturity", status = cyto_res$status, output_dir = cyto_res$output_dir, notes = if (!is.null(cyto_res$error)) cyto_res$error else "Human CytoTRACE2 maturity score; pseudotime = root-oriented 1-score", stringsAsFactors = FALSE),
  data.frame(method = "STREAM", tier = "expression_sufficient", status = "not_configured", output_dir = NA_character_, notes = "STREAM package/runtime not configured in this workspace; skipped rather than installing unverified package", stringsAsFactors = FALSE)
)
inventory <- do.call(rbind, status_rows)
write_tsv(inventory, file.path(OUTPUT_DIR, "first_tier_method_inventory.tsv"))
status_pdf <- file.path(OUTPUT_DIR, "first_tier_method_status_barplot.pdf")
plot_method_status(inventory, status_pdf)
register_method_figure("First-tier trajectory inventory", status_pdf, file.path(OUTPUT_DIR, "first_tier_method_inventory.tsv"), "method_availability_and_run_status", "First-tier expression-sufficient trajectory method status", MANIFEST, extra_context = c("This inventory tracks which user-requested first-tier methods were already run, newly run, or skipped due to missing/unconfigured runtime.", "STREAM remains not configured; PAGA is topology-only and is not forced into single-cell pseudotime."))

manifests <- smcanno_viz_collect_companion_manifests(OUTPUT_DIR)
smcanno_viz_write_manifest_index(manifests, file.path(OUTPUT_DIR, "first_tier_llm_manifest_index.csv"))

readme <- c(
  "# smc_anno first-tier expression-sufficient trajectory methods",
  "",
  sprintf("- Run time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Input: `%s`", INPUT_RDATA),
  sprintf("- Output: `%s`", OUTPUT_DIR),
  sprintf("- Cluster column: `%s`; root cluster: `%s`; phenotype column: `%s`", CLUSTER_COL, ROOT_CLUSTER, PHENOTYPE_COL),
  sprintf("- Expression h5ad: `%s`", export_info$h5ad_path),
  "",
  "## Method status",
  smcanno_viz_markdown_table(inventory, max_rows = nrow(inventory)),
  "",
  "## Notes",
  "- Slingshot, Monocle3, and PAGA were reused from the full main trajectory run.",
  "- DPT, PHATE, Palantir, VIA, CellRank kernels, TSCAN, SCORPIUS, and CytoTRACE2 were attempted in this supplemental batch.",
  "- PAGA and CellRank kernels are graph/topology/kernel outputs; PAGA is not treated as a scalar pseudotime method.",
  "- Non-tradeSeq gene trend plots keep observed binned mean-expression points and overlay a display-layer smoother for readability; the smoother is not a tradeSeq-like fitted model.",
  "- Every generated figure has paired data CSV/TSV and queued LLM companion files via the manifest system."
)
writeLines(readme, file.path(OUTPUT_DIR, "README.md"), useBytes = TRUE)

cat("\nFirst-tier method inventory:\n")
print(inventory)
cat("\nSaved outputs to: ", OUTPUT_DIR, "\n", sep = "")
cat(strrep("=", 80), "\n", sep = "")
