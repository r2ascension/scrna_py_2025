#!/usr/bin/env Rscript
# ==============================================================================
# smc_anno.Rdata Slingshot Trajectory Analysis (2026-04-28)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
})

TRAJECTORY_BRANCH_HELPER_PATH <- "/home/h2048/script/R/trajectory_branch_helper_20260428_v1.R"
TRADESEQ_HELPER_PATH <- "/home/h2048/script/R/smc_anno_tradeseq_helper_20260429.R"
TRAJECTORY_INTERPRET_HELPER_PATH <- "/home/h2048/script/R/stromal_smc_trajectory_interpret_helper_20260429.R"
VISUAL_LLM_HELPER_PATH <- "/home/h2048/script/R/smc_anno_visual_llm_helper_20260506.R"
if (!exists("pa_build_slingshot_trajectory_packet", mode = "function")) {
  source(TRAJECTORY_BRANCH_HELPER_PATH)
}
if (!exists("smcanno_ts_run_bundle", mode = "function")) {
  source(TRADESEQ_HELPER_PATH)
}
if (!exists("smcti_run_monocle3_bundle", mode = "function")) {
  source(TRAJECTORY_INTERPRET_HELPER_PATH)
}
if (!exists("smcanno_viz_register_figure", mode = "function")) {
  source(VISUAL_LLM_HELPER_PATH)
}

if (!exists("INPUT_RDATA")) {
  INPUT_RDATA <- "/home/h2048/temp/smc_anno.Rdata"
}
if (!exists("OBJECT_NAME")) {
  OBJECT_NAME <- "smc_anno"
}
if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- "/home/h2048/temp/smc_anno_slingshot_20260428"
}
if (!exists("ASSAY_NAME")) {
  ASSAY_NAME <- "RNA"
}
if (!exists("REDUCTION_NAME")) {
  REDUCTION_NAME <- "umap"
}
if (!exists("CLUSTER_COL")) {
  CLUSTER_COL <- "celltype"
}
if (!exists("START_CLUSTER")) {
  START_CLUSTER <- "Pericyte"
}
if (!exists("END_CLUSTERS")) {
  END_CLUSTERS <- c("FN1+ SMC", "RERGL+ SMC", "VCAN+ SMC")
}
if (!exists("SAVE_SCE")) {
  SAVE_SCE <- TRUE
}
if (!exists("RUN_TRADESEQ")) {
  RUN_TRADESEQ <- TRUE
}
if (!exists("TRADESEQ_SUBDIR")) {
  TRADESEQ_SUBDIR <- "tradeSeq"
}
if (!exists("TRADESEQ_MAX_GENES")) {
  TRADESEQ_MAX_GENES <- 1000L
}
if (!exists("TRADESEQ_MIN_CELLS_EXPRESSED")) {
  TRADESEQ_MIN_CELLS_EXPRESSED <- 20L
}
if (!exists("TRADESEQ_GENE_PANEL")) {
  TRADESEQ_GENE_PANEL <- NULL
}
if (!exists("TRADESEQ_MARKER_PANEL")) {
  TRADESEQ_MARKER_PANEL <- smcanno_ts_default_marker_panel
}
if (!exists("TRADESEQ_K_VALUES")) {
  TRADESEQ_K_VALUES <- 3:7
}
if (!exists("TRADESEQ_EVALUATEK_N_GENES")) {
  TRADESEQ_EVALUATEK_N_GENES <- 200L
}
if (!exists("TRADESEQ_NK")) {
  TRADESEQ_NK <- 6L
}
if (!exists("TRADESEQ_L2FC")) {
  TRADESEQ_L2FC <- 0
}
if (!exists("TRADESEQ_ENABLE_PLOTS")) {
  TRADESEQ_ENABLE_PLOTS <- TRUE
}
if (!exists("TRADESEQ_EXCLUDE_NONINFORMATIVE")) {
  TRADESEQ_EXCLUDE_NONINFORMATIVE <- TRUE
}
if (!exists("TRADESEQ_VERBOSE")) {
  TRADESEQ_VERBOSE <- FALSE
}
if (!exists("RUN_MONOCLE3")) {
  RUN_MONOCLE3 <- TRUE
}
if (!exists("MONOCLE3_SUBDIR")) {
  MONOCLE3_SUBDIR <- "monocle3"
}
if (!exists("MONOCLE3_NUM_DIM")) {
  MONOCLE3_NUM_DIM <- 30L
}
if (!exists("RUN_MONOCLE2")) {
  RUN_MONOCLE2 <- TRUE
}
if (!exists("MONOCLE2_SUBDIR")) {
  MONOCLE2_SUBDIR <- "monocle2"
}
if (!exists("MONOCLE2_MAX_ORDERING_GENES")) {
  MONOCLE2_MAX_ORDERING_GENES <- 2000L
}
if (!exists("RUN_PAGA")) {
  RUN_PAGA <- TRUE
}
if (!exists("PAGA_SUBDIR")) {
  PAGA_SUBDIR <- "paga"
}
if (!exists("PAGA_PYTHON_BIN")) {
  PAGA_PYTHON_BIN <- "/home/h2048/miniconda3/envs/bbknn_env/bin/python"
}
if (!exists("RUN_CROSS_METHOD_SUMMARY")) {
  RUN_CROSS_METHOD_SUMMARY <- TRUE
}
if (!exists("PHENOTYPE_COL")) {
  PHENOTYPE_COL <- NULL
}
if (!exists("RUN_GENE_PSEUDOTIME_TRENDS")) {
  RUN_GENE_PSEUDOTIME_TRENDS <- TRUE
}
if (!exists("GENE_TREND_MAX_GENES")) {
  GENE_TREND_MAX_GENES <- 16L
}
if (!exists("RUN_PHENOTYPE_PSEUDOTIME_COMPARISON")) {
  RUN_PHENOTYPE_PSEUDOTIME_COMPARISON <- TRUE
}
if (!exists("RUN_FIGURE_LLM_BATCH")) {
  RUN_FIGURE_LLM_BATCH <- smcanno_viz_env_flag("SMCANNO_RUN_FIGURE_LLM_BATCH", TRUE)
}
if (!exists("FIGURE_LLM_BATCH_ASYNC")) {
  FIGURE_LLM_BATCH_ASYNC <- TRUE
}
if (!exists("FIGURE_LLM_BATCH_MODEL")) {
  FIGURE_LLM_BATCH_MODEL <- "deepseek-reasoner"
}
if (!exists("FIGURE_LLM_BATCH_TIMEOUT_SEC")) {
  FIGURE_LLM_BATCH_TIMEOUT_SEC <- 180
}
if (!exists("FIGURE_LLM_RUNNER_SCRIPT")) {
  FIGURE_LLM_RUNNER_SCRIPT <- "/home/h2048/script/R/smc_anno_figure_llm_batch_20260507.R"
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", as.character(x))
}

write_tsv <- function(df, path) {
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

safe_numeric_summary <- function(x, fun) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  fun(x)
}

ensure_named_matrix <- function(x, prefix) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0(prefix, seq_len(ncol(x)))
  }
  x
}

pick_seurat_object <- function(env, object_name = NULL) {
  env_names <- ls(env, all.names = TRUE)
  if (!is.null(object_name) && nzchar(object_name) && object_name %in% env_names) {
    obj <- get(object_name, envir = env)
    if (inherits(obj, "Seurat")) {
      return(list(name = object_name, object = obj))
    }
    stop(sprintf("Object '%s' exists but is not a Seurat object", object_name), call. = FALSE)
  }
  seurat_names <- env_names[vapply(env_names, function(nm) inherits(get(nm, envir = env), "Seurat"), logical(1))]
  if (length(seurat_names) == 0L) {
    stop("No Seurat object found in INPUT_RDATA", call. = FALSE)
  }
  if (length(seurat_names) > 1L) {
    stop(sprintf("Multiple Seurat objects found: %s", paste(seurat_names, collapse = ", ")), call. = FALSE)
  }
  list(name = seurat_names[[1]], object = get(seurat_names[[1]], envir = env))
}

build_primary_lineage <- function(weight_mat) {
  apply(weight_mat, 1, function(x) {
    finite_idx <- which(is.finite(x))
    if (length(finite_idx) == 0L) {
      return(NA_character_)
    }
    best_idx <- finite_idx[which.max(x[finite_idx])]
    if (!is.finite(x[best_idx]) || x[best_idx] <= 0) {
      return(NA_character_)
    }
    colnames(weight_mat)[best_idx]
  })
}

plot_slingshot_clusters <- function(sce, labels, output_path, title_text) {
  coords <- SingleCellExperiment::reducedDim(sce, "TRAJECTORY")
  palette_vals <- setNames(
    grDevices::hcl.colors(length(levels(labels)), palette = "Dark 3"),
    levels(labels)
  )
  grDevices::pdf(output_path, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(
    coords,
    col = palette_vals[as.character(labels)],
    pch = 16,
    asp = 1,
    xlab = paste0(REDUCTION_NAME, "_1"),
    ylab = paste0(REDUCTION_NAME, "_2"),
    main = title_text
  )
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  legend(
    "topright",
    legend = levels(labels),
    col = palette_vals,
    pch = 16,
    cex = 0.8,
    bty = "n"
  )
}

plot_slingshot_pseudotime <- function(sce, pseudotime_vec, output_path, title_text) {
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
  plot(
    coords,
    col = point_cols,
    pch = 16,
    asp = 1,
    xlab = paste0(REDUCTION_NAME, "_1"),
    ylab = paste0(REDUCTION_NAME, "_2"),
    main = title_text
  )
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
}

cat("\n", strrep("=", 80), "\n", sep = "")
cat("SMC_ANNO RDATA SLINGSHOT TRAJECTORY ANALYSIS\n")
cat(strrep("=", 80), "\n", sep = "")
cat("Started    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Input Rdata: ", INPUT_RDATA, "\n", sep = "")
cat("Output dir : ", OUTPUT_DIR, "\n", sep = "")
cat("Reduction  : ", REDUCTION_NAME, "\n", sep = "")
cat("Cluster col: ", CLUSTER_COL, "\n", sep = "")
cat("Start      : ", START_CLUSTER, "\n", sep = "")
cat("End        : ", paste(END_CLUSTERS, collapse = ", "), "\n\n", sep = "")
cat("tradeSeq   : ", if (isTRUE(RUN_TRADESEQ)) "enabled" else "disabled", "\n", sep = "")
if (isTRUE(RUN_TRADESEQ)) {
  cat("  max genes      : ", TRADESEQ_MAX_GENES, "\n", sep = "")
  cat("  min cells expr : ", TRADESEQ_MIN_CELLS_EXPRESSED, "\n", sep = "")
  cat("  evaluateK k    : ", paste(TRADESEQ_K_VALUES, collapse = ", "), "\n", sep = "")
  cat("  fitGAM nknots  : ", TRADESEQ_NK, "\n\n", sep = "")
}
cat("Monocle3   : ", if (isTRUE(RUN_MONOCLE3)) "enabled" else "disabled", "\n", sep = "")
cat("Monocle2   : ", if (isTRUE(RUN_MONOCLE2)) "enabled" else "disabled", "\n", sep = "")
cat("PAGA       : ", if (isTRUE(RUN_PAGA)) "enabled" else "disabled", "\n\n", sep = "")
cat("Figure LLM : ", if (isTRUE(RUN_FIGURE_LLM_BATCH)) if (isTRUE(FIGURE_LLM_BATCH_ASYNC)) "queued + async runner" else "queued + sync runner" else "queued only", "\n\n", sep = "")

if (!file.exists(INPUT_RDATA)) {
  stop(sprintf("Input Rdata does not exist: %s", INPUT_RDATA), call. = FALSE)
}

load_env <- new.env(parent = emptyenv())
loaded_names <- load(INPUT_RDATA, envir = load_env)
cat("Loaded objects: ", paste(loaded_names, collapse = ", "), "\n", sep = "")

picked <- pick_seurat_object(load_env, OBJECT_NAME)
object_name_used <- picked$name
obj <- picked$object

if (!ASSAY_NAME %in% names(obj@assays)) {
  ASSAY_NAME <- Seurat::DefaultAssay(obj)
}
if (!REDUCTION_NAME %in% names(obj@reductions)) {
  stop(sprintf("Reduction '%s' not found in Seurat object", REDUCTION_NAME), call. = FALSE)
}
if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop(sprintf("Metadata column '%s' not found in Seurat object", CLUSTER_COL), call. = FALSE)
}

cluster_values <- as.character(obj@meta.data[[CLUSTER_COL]])
cluster_values[is.na(cluster_values)] <- ""
if (any(!nzchar(trimws(cluster_values)))) {
  stop(sprintf("Metadata column '%s' contains NA/empty labels", CLUSTER_COL), call. = FALSE)
}
cluster_labels <- factor(cluster_values)

missing_start <- setdiff(START_CLUSTER, levels(cluster_labels))
missing_end <- setdiff(as.character(END_CLUSTERS), levels(cluster_labels))
if (length(missing_start) > 0L) {
  stop(sprintf("Start cluster not found: %s", paste(missing_start, collapse = ", ")), call. = FALSE)
}
if (length(missing_end) > 0L) {
  stop(sprintf("End clusters not found: %s", paste(missing_end, collapse = ", ")), call. = FALSE)
}

OUTPUT_DIR <- ensure_dir(OUTPUT_DIR)
FIG_DIR <- ensure_dir(file.path(OUTPUT_DIR, "figures"))
SLINGSHOT_FIGURE_MANIFEST <- file.path(FIG_DIR, "slingshot_figure_companion_manifest.tsv")
if (file.exists(SLINGSHOT_FIGURE_MANIFEST)) {
  unlink(SLINGSHOT_FIGURE_MANIFEST)
}

cat("Using object : ", object_name_used, "\n", sep = "")
cat("Cells/features: ", ncol(obj), " / ", nrow(obj), "\n", sep = "")
cat("Assay        : ", ASSAY_NAME, "\n", sep = "")
cat("Label counts:\n")
print(sort(table(cluster_labels), decreasing = TRUE))

phenotype_detection <- smcti_infer_phenotype_col(obj@meta.data)
if (is.null(PHENOTYPE_COL) || length(PHENOTYPE_COL) == 0L || is.na(PHENOTYPE_COL) || !nzchar(PHENOTYPE_COL)) {
  PHENOTYPE_COL <- phenotype_detection$column
}
write_tsv(phenotype_detection$candidates, file.path(OUTPUT_DIR, "phenotype_column_candidates.tsv"))
if (!is.na(PHENOTYPE_COL) && nzchar(PHENOTYPE_COL) && PHENOTYPE_COL %in% colnames(obj@meta.data)) {
  cat("Phenotype col: ", PHENOTYPE_COL, "\n", sep = "")
  cat("Phenotype counts:\n")
  print(sort(table(as.character(obj@meta.data[[PHENOTYPE_COL]]), useNA = "ifany"), decreasing = TRUE))
} else {
  PHENOTYPE_COL <- NA_character_
  cat("Phenotype col: not detected\n")
}

sce <- suppressWarnings(as.SingleCellExperiment(obj, assay = ASSAY_NAME))
trajectory_coords <- Seurat::Embeddings(obj, reduction = REDUCTION_NAME)[colnames(obj), , drop = FALSE]
SingleCellExperiment::reducedDim(sce, "TRAJECTORY") <- trajectory_coords
SummarizedExperiment::colData(sce)[[CLUSTER_COL]] <- cluster_labels

sce <- slingshot(
  sce,
  clusterLabels = CLUSTER_COL,
  reducedDim = "TRAJECTORY",
  start.clus = START_CLUSTER,
  end.clus = as.character(END_CLUSTERS)
)

pseudotime_mat <- ensure_named_matrix(slingPseudotime(sce), prefix = "Lineage")
curve_weight_mat <- ensure_named_matrix(slingCurveWeights(sce), prefix = "Lineage")
lineages <- slingLineages(sce)

meta_cols <- intersect(
  unique(c("celltype", "group", "seurat_clusters", "Contractile1", "Synthetic1", "orig.ident", "sample", PHENOTYPE_COL)),
  colnames(obj@meta.data)
)

pseudotime_table <- data.frame(
  cell_id = colnames(obj),
  obj@meta.data[colnames(obj), meta_cols, drop = FALSE],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (lineage_name in colnames(pseudotime_mat)) {
  pseudotime_table[[paste0(lineage_name, "_pseudotime")]] <- pseudotime_mat[, lineage_name]
}
for (lineage_name in colnames(curve_weight_mat)) {
  pseudotime_table[[paste0(lineage_name, "_curve_weight")]] <- curve_weight_mat[, lineage_name]
}
pseudotime_table$primary_lineage <- build_primary_lineage(curve_weight_mat)

plot_coord_df <- data.frame(
  cell_id = rownames(trajectory_coords),
  dim1 = as.numeric(trajectory_coords[, 1]),
  dim2 = as.numeric(trajectory_coords[, 2]),
  cluster_label = as.character(cluster_labels),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
slingshot_plot_df <- merge(plot_coord_df, pseudotime_table, by = "cell_id", all.x = TRUE, sort = FALSE)
write_tsv(slingshot_plot_df, file.path(FIG_DIR, "slingshot_plot_cells_data.tsv"))

lineage_summary <- do.call(rbind, lapply(seq_along(lineages), function(i) {
  lineage_name <- names(lineages)[i]
  lineage_path <- lineages[[i]]
  data.frame(
    lineage_id = lineage_name,
    root_state = lineage_path[[1]],
    terminal_state = utils::tail(lineage_path, 1),
    lineage_path = paste(lineage_path, collapse = " -> "),
    n_cells_with_pseudotime = sum(is.finite(pseudotime_mat[, lineage_name])),
    n_cells_curve_weight_ge_0_5 = sum(curve_weight_mat[, lineage_name] >= 0.5, na.rm = TRUE),
    median_pseudotime = safe_numeric_summary(pseudotime_mat[, lineage_name], stats::median),
    max_pseudotime = safe_numeric_summary(pseudotime_mat[, lineage_name], max),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))

cluster_summary_list <- lapply(levels(cluster_labels), function(cluster_name) {
  idx <- which(cluster_labels == cluster_name)
  out <- data.frame(
    cluster_label = cluster_name,
    n_cells = length(idx),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (lineage_name in colnames(pseudotime_mat)) {
    out[[paste0(lineage_name, "_median_pseudotime")]] <- safe_numeric_summary(pseudotime_mat[idx, lineage_name], stats::median)
    out[[paste0(lineage_name, "_mean_curve_weight")]] <- safe_numeric_summary(curve_weight_mat[idx, lineage_name], mean)
    out[[paste0(lineage_name, "_cells_weight_ge_0_5")]] <- sum(curve_weight_mat[idx, lineage_name] >= 0.5, na.rm = TRUE)
  }
  out
})
cluster_summary <- do.call(rbind, cluster_summary_list)

trajectory_packet <- pa_build_slingshot_trajectory_packet(
  pseudotime_table = pseudotime_table,
  lineage_summary = lineage_summary,
  branch_summary = cluster_summary,
  source_path = normalizePath(INPUT_RDATA, winslash = "/", mustWork = FALSE),
  root_state = START_CLUSTER,
  terminal_states = as.character(END_CLUSTERS)
)
trajectory_branch_packet <- pa_build_trajectory_branch_packet(
  primary_trajectory = trajectory_packet,
  preferred_engine = "Slingshot"
)

write_tsv(pseudotime_table, file.path(OUTPUT_DIR, "slingshot_pseudotime.tsv"))
write_tsv(lineage_summary, file.path(OUTPUT_DIR, "slingshot_lineage_summary.tsv"))
write_tsv(cluster_summary, file.path(OUTPUT_DIR, "slingshot_cluster_summary.tsv"))
saveRDS(trajectory_packet, file.path(OUTPUT_DIR, "slingshot_trajectory_packet.rds"))
saveRDS(trajectory_branch_packet, file.path(OUTPUT_DIR, "trajectory_branch_packet.rds"))
if (isTRUE(SAVE_SCE)) {
  saveRDS(sce, file.path(OUTPUT_DIR, "slingshot_sce.rds"))
}

slingshot_cross_method_res <- list(
  pseudotime_table = pseudotime_table,
  lineage_summary = lineage_summary,
  pseudotime_mat = pseudotime_mat,
  curve_weight_mat = curve_weight_mat,
  cluster_col = CLUSTER_COL
)

tradeSeq_result <- NULL
if (isTRUE(RUN_TRADESEQ)) {
  cat("\nRunning tradeSeq on smc_anno trajectory...\n")
  tradeSeq_dir <- ensure_dir(file.path(OUTPUT_DIR, TRADESEQ_SUBDIR))
  tradeSeq_result <- smcanno_ts_run_bundle(
    seurat_obj = obj,
    slingshot_sce = sce,
    output_dir = tradeSeq_dir,
    assay_name = ASSAY_NAME,
    max_genes = TRADESEQ_MAX_GENES,
    min_cells_expressed = TRADESEQ_MIN_CELLS_EXPRESSED,
    gene_panel = TRADESEQ_GENE_PANEL,
    marker_panel = TRADESEQ_MARKER_PANEL,
    k_values = TRADESEQ_K_VALUES,
    evaluatek_n_genes = TRADESEQ_EVALUATEK_N_GENES,
    nknots = TRADESEQ_NK,
    l2fc = TRADESEQ_L2FC,
    enable_plots = TRADESEQ_ENABLE_PLOTS,
    exclude_noninformative = TRADESEQ_EXCLUDE_NONINFORMATIVE,
    verbose = TRADESEQ_VERBOSE
  )
  cat("tradeSeq outputs: ", tradeSeq_dir, "\n", sep = "")
}

monocle3_result <- NULL
if (isTRUE(RUN_MONOCLE3)) {
  cat("\nRunning Monocle3 on smc_anno trajectory...\n")
  monocle3_result <- smcti_run_monocle3_bundle(
    seurat_obj = obj,
    output_dir = ensure_dir(file.path(OUTPUT_DIR, MONOCLE3_SUBDIR)),
    assay_name = ASSAY_NAME,
    cluster_col = CLUSTER_COL,
    start_cluster = START_CLUSTER,
    end_clusters = END_CLUSTERS,
    reduction_name = REDUCTION_NAME,
    num_dim = MONOCLE3_NUM_DIM
  )
  cat("Monocle3 status: ", monocle3_result$status, "\n", sep = "")
}

monocle2_result <- NULL
if (isTRUE(RUN_MONOCLE2)) {
  cat("\nRunning Monocle2 on smc_anno trajectory...\n")
  monocle2_result <- smcti_run_monocle2_bundle(
    seurat_obj = obj,
    output_dir = ensure_dir(file.path(OUTPUT_DIR, MONOCLE2_SUBDIR)),
    assay_name = ASSAY_NAME,
    cluster_col = CLUSTER_COL,
    start_cluster = START_CLUSTER,
    end_clusters = END_CLUSTERS,
    max_ordering_genes = MONOCLE2_MAX_ORDERING_GENES,
    marker_panel = TRADESEQ_MARKER_PANEL
  )
  cat("Monocle2 status: ", monocle2_result$status, "\n", sep = "")
}

paga_result <- NULL
if (isTRUE(RUN_PAGA)) {
  cat("\nRunning PAGA on smc_anno trajectory...\n")
  paga_result <- smcti_run_paga_bundle(
    seurat_obj = obj,
    output_dir = ensure_dir(file.path(OUTPUT_DIR, PAGA_SUBDIR)),
    groupby_key = CLUSTER_COL,
    reduction_name = REDUCTION_NAME,
    assay_name = ASSAY_NAME,
    python_bin = PAGA_PYTHON_BIN
  )
  cat("PAGA status: ", paga_result$status, "\n", sep = "")
}

cross_method_summary <- NULL
if (isTRUE(RUN_CROSS_METHOD_SUMMARY)) {
  cross_method_summary <- smcti_build_cross_method_summary(
    slingshot_res = slingshot_cross_method_res,
    monocle3_res = monocle3_result,
    monocle2_res = monocle2_result,
    paga_res = paga_result
  )
  write_tsv(cross_method_summary, file.path(OUTPUT_DIR, "cross_method_trajectory_summary.tsv"))
}

gene_pseudotime_trend_results <- list()
phenotype_pseudotime_comparison <- NULL
if (isTRUE(RUN_GENE_PSEUDOTIME_TRENDS)) {
  trend_root <- ensure_dir(file.path(OUTPUT_DIR, "gene_pseudotime_trends"))
  trend_genes <- TRADESEQ_MARKER_PANEL
  if (!is.null(tradeSeq_result) && !is.null(tradeSeq_result$visual_bundle$top_tests) && nrow(tradeSeq_result$visual_bundle$top_tests) > 0L) {
    trend_genes <- unique(c(trend_genes, tradeSeq_result$visual_bundle$top_tests$gene))
  }
  trend_genes <- smcti_filter_gene_panel_for_trends(trend_genes, seurat_obj = obj, max_genes = GENE_TREND_MAX_GENES, exclude_noninformative = TRUE)
  write_tsv(data.frame(gene = trend_genes, stringsAsFactors = FALSE), file.path(trend_root, "gene_trend_panel.tsv"))
  if (length(trend_genes) > 0L) {
    cat("\nWriting gene-over-pseudotime trend bundles for Slingshot/Monocle methods...\n")
    gene_pseudotime_trend_results$Slingshot <- smcti_write_gene_pseudotime_trend_bundle(
      seurat_obj = obj,
      pseudotime_table = pseudotime_table,
      output_dir = ensure_dir(file.path(trend_root, "Slingshot")),
      assay_name = ASSAY_NAME,
      method_name = "Slingshot",
      gene_panel = trend_genes,
      cluster_col = CLUSTER_COL,
      phenotype_col = PHENOTYPE_COL,
      max_genes = GENE_TREND_MAX_GENES
    )
    if (!is.null(monocle3_result) && identical(monocle3_result$status, "ok")) {
      gene_pseudotime_trend_results$Monocle3 <- smcti_write_gene_pseudotime_trend_bundle(
        seurat_obj = obj,
        pseudotime_table = monocle3_result$pseudotime_table,
        output_dir = ensure_dir(file.path(trend_root, "Monocle3")),
        assay_name = ASSAY_NAME,
        method_name = "Monocle3",
        gene_panel = trend_genes,
        cluster_col = CLUSTER_COL,
        phenotype_col = PHENOTYPE_COL,
        max_genes = GENE_TREND_MAX_GENES
      )
    }
    if (!is.null(monocle2_result) && identical(monocle2_result$status, "ok")) {
      gene_pseudotime_trend_results$Monocle2 <- smcti_write_gene_pseudotime_trend_bundle(
        seurat_obj = obj,
        pseudotime_table = monocle2_result$pseudotime_table,
        output_dir = ensure_dir(file.path(trend_root, "Monocle2")),
        assay_name = ASSAY_NAME,
        method_name = "Monocle2",
        gene_panel = trend_genes,
        cluster_col = CLUSTER_COL,
        phenotype_col = PHENOTYPE_COL,
        max_genes = GENE_TREND_MAX_GENES
      )
    }
  }
}

if (isTRUE(RUN_PHENOTYPE_PSEUDOTIME_COMPARISON) && !is.na(PHENOTYPE_COL) && nzchar(PHENOTYPE_COL)) {
  method_tables <- list(Slingshot = pseudotime_table)
  if (!is.null(monocle3_result) && identical(monocle3_result$status, "ok")) method_tables$Monocle3 <- monocle3_result$pseudotime_table
  if (!is.null(monocle2_result) && identical(monocle2_result$status, "ok")) method_tables$Monocle2 <- monocle2_result$pseudotime_table
  phenotype_pseudotime_comparison <- smcti_write_phenotype_pseudotime_comparison_bundle(
    method_tables = method_tables,
    seurat_obj = obj,
    output_dir = ensure_dir(file.path(OUTPUT_DIR, "phenotype_pseudotime_comparison")),
    phenotype_col = PHENOTYPE_COL,
    cluster_col = CLUSTER_COL
  )
}

slingshot_cluster_pdf <- file.path(FIG_DIR, "slingshot_umap_clusters.pdf")
plot_slingshot_clusters(
  sce = sce,
  labels = cluster_labels,
  output_path = slingshot_cluster_pdf,
  title_text = sprintf("smc_anno Slingshot (%s on %s)", CLUSTER_COL, REDUCTION_NAME)
)
smcanno_viz_register_figure(
  figure_path = slingshot_cluster_pdf,
  data = slingshot_plot_df,
  method = "Slingshot",
  figure_type = "embedding_clusters_with_curves",
  title = sprintf("smc_anno Slingshot clusters on %s", REDUCTION_NAME),
  extra_context = c(
    sprintf("Start cluster: %s", START_CLUSTER),
    sprintf("Terminal clusters: %s", paste(as.character(END_CLUSTERS), collapse = "; ")),
    "The paired CSV contains embedding coordinates, cluster labels, pseudotime, and curve weights for each cell."
  ),
  manifest_path = SLINGSHOT_FIGURE_MANIFEST
)

for (lineage_name in colnames(pseudotime_mat)) {
  lineage_row <- lineage_summary[lineage_summary$lineage_id == lineage_name, , drop = FALSE]
  lineage_pdf <- file.path(FIG_DIR, sprintf("slingshot_%s_pseudotime.pdf", sanitize_name(lineage_name)))
  lineage_plot_df <- slingshot_plot_df[, c("cell_id", "dim1", "dim2", "cluster_label", lineage_name_col <- paste0(lineage_name, "_pseudotime"), paste0(lineage_name, "_curve_weight"), "primary_lineage"), drop = FALSE]
  colnames(lineage_plot_df)[colnames(lineage_plot_df) == lineage_name_col] <- "lineage_pseudotime"
  colnames(lineage_plot_df)[colnames(lineage_plot_df) == paste0(lineage_name, "_curve_weight")] <- "lineage_curve_weight"
  lineage_plot_df$lineage_id <- lineage_name
  lineage_plot_df$lineage_path <- lineage_row$lineage_path[[1]]
  plot_slingshot_pseudotime(
    sce = sce,
    pseudotime_vec = pseudotime_mat[, lineage_name],
    output_path = lineage_pdf,
    title_text = sprintf("%s: %s", lineage_name, lineage_row$lineage_path[[1]])
  )
  smcanno_viz_register_figure(
    figure_path = lineage_pdf,
    data = lineage_plot_df,
    method = "Slingshot",
    figure_type = "lineage_pseudotime_embedding",
    title = sprintf("Slingshot %s pseudotime: %s", lineage_name, lineage_row$lineage_path[[1]]),
    extra_context = c(
      sprintf("Lineage path: %s", lineage_row$lineage_path[[1]]),
      sprintf("Cells with finite pseudotime: %s", lineage_row$n_cells_with_pseudotime[[1]]),
      "Cells with low or zero curve weight should not be over-interpreted as part of this branch."
    ),
    manifest_path = SLINGSHOT_FIGURE_MANIFEST
  )
}

figure_llm_manifests <- smcanno_viz_collect_companion_manifests(OUTPUT_DIR)
figure_llm_index_path <- file.path(OUTPUT_DIR, "figure_llm_manifest_index.csv")
smcanno_viz_write_manifest_index(figure_llm_manifests, figure_llm_index_path)
figure_llm_launch <- list(status = "skipped_disabled", manifest_index = figure_llm_index_path, manifest_count = length(figure_llm_manifests))
if (isTRUE(RUN_FIGURE_LLM_BATCH) && length(figure_llm_manifests) > 0L) {
  cat("\nStarting decoupled figure LLM runner", if (isTRUE(FIGURE_LLM_BATCH_ASYNC)) " asynchronously" else " synchronously", "...\n", sep = "")
  figure_llm_launch <- tryCatch(
    smcanno_viz_launch_llm_batch(
      manifest_paths = figure_llm_manifests,
      output_dir = OUTPUT_DIR,
      runner_script = FIGURE_LLM_RUNNER_SCRIPT,
      async = FIGURE_LLM_BATCH_ASYNC,
      model = FIGURE_LLM_BATCH_MODEL,
      timeout_sec = FIGURE_LLM_BATCH_TIMEOUT_SEC,
      log_path = file.path(OUTPUT_DIR, "figure_llm_batch.log")
    ),
    error = function(e) list(status = "error", error = conditionMessage(e), manifest_index = figure_llm_index_path, manifest_count = length(figure_llm_manifests))
  )
  cat("Figure LLM runner status: ", figure_llm_launch$status, "\n", sep = "")
  if (!is.null(figure_llm_launch$log_path)) cat("Figure LLM runner log   : ", figure_llm_launch$log_path, "\n", sep = "")
}
smcanno_viz_write_json(figure_llm_launch, file.path(OUTPUT_DIR, "figure_llm_batch_launch_status.json"))

readme_lines <- c(
  "# smc_anno.Rdata Slingshot trajectory run",
  "",
  sprintf("- Run date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Input Rdata: `%s`", INPUT_RDATA),
  sprintf("- Loaded Seurat object: `%s`", object_name_used),
  sprintf("- Output dir: `%s`", OUTPUT_DIR),
  sprintf("- Assay: `%s`", ASSAY_NAME),
  sprintf("- Reduction: `%s`", REDUCTION_NAME),
  sprintf("- Cluster column: `%s`", CLUSTER_COL),
  sprintf("- Disease phenotype column: `%s`", if (!is.na(PHENOTYPE_COL) && nzchar(PHENOTYPE_COL)) PHENOTYPE_COL else "not detected"),
  sprintf("- Start cluster: `%s`", START_CLUSTER),
  sprintf("- End clusters: `%s`", paste(as.character(END_CLUSTERS), collapse = "`, `")),
  "",
  "## Lineages",
  apply(lineage_summary, 1, function(x) sprintf(
    "- %s: %s (n pseudotime cells=%s)",
    x[["lineage_id"]], x[["lineage_path"]], x[["n_cells_with_pseudotime"]]
  )),
  "",
  "## tradeSeq",
  if (isTRUE(RUN_TRADESEQ) && !is.null(tradeSeq_result)) sprintf(
    "- Enabled: yes (n genes=%s; nknots=%s; evaluateK status=%s; non-informative gene exclusion=%s)",
    tradeSeq_result$metadata$n_genes_fit[[1]],
    tradeSeq_result$metadata$nknots_fit[[1]],
    tradeSeq_result$metadata$evaluateK_status[[1]],
    tradeSeq_result$metadata$exclude_noninformative[[1]]
  ) else "- Enabled: no",
  "",
  "## Cross-method checks",
  sprintf("- Monocle3: %s", if (!is.null(monocle3_result)) monocle3_result$status else "disabled"),
  sprintf("- Monocle2: %s", if (!is.null(monocle2_result)) monocle2_result$status else "disabled"),
  sprintf("- PAGA: %s", if (!is.null(paga_result)) paga_result$status else "disabled"),
  if (isTRUE(RUN_GENE_PSEUDOTIME_TRENDS)) "- Non-tradeSeq gene trend plots keep observed 30-bin mean-expression points and overlay a display-layer smoother for readability; this smoother is not equivalent to a tradeSeq fitted GAM." else NULL,
  "",
  "## Files",
  "- `slingshot_pseudotime.tsv`: per-cell pseudotime + curve weights",
  "- `slingshot_lineage_summary.tsv`: lineage-level summary",
  "- `slingshot_cluster_summary.tsv`: cluster-level summary",
  "- `slingshot_trajectory_packet.rds`: primary Slingshot packet",
  "- `trajectory_branch_packet.rds`: trajectory branch wrapper",
  "- `slingshot_sce.rds`: SingleCellExperiment with fitted curves",
  "- `figures/`: Slingshot PDFs plus paired `*_data.csv`, `*_LLM_PROMPT.md`, queued/rule `*_LLM_ANALYSIS.md`, and `*_LLM_STATUS.json`",
  if (isTRUE(RUN_TRADESEQ)) "- `tradeSeq/`: evaluateK, fitGAM, tradeSeq tests, enhanced plots, paired CSVs, queued prompts, and batch-updatable LLM analyses" else NULL,
  if (isTRUE(RUN_MONOCLE3)) "- `monocle3/`: Monocle3 pseudotime, partitions, trajectory plots, paired CSVs, and queued LLM companion files" else NULL,
  if (isTRUE(RUN_MONOCLE2)) "- `monocle2/`: Monocle2 pseudotime, states, trajectory plots, paired CSVs, and queued LLM companion files" else NULL,
  if (isTRUE(RUN_PAGA)) "- `paga/`: PAGA connectivity matrix, edge table, topology packet, heatmap/network/embedding plots, paired CSVs, and queued LLM companion files" else NULL,
  if (isTRUE(RUN_CROSS_METHOD_SUMMARY)) "- `cross_method_trajectory_summary.tsv`: Slingshot vs Monocle/PAGA summary" else NULL,
  if (isTRUE(RUN_GENE_PSEUDOTIME_TRENDS)) "- `gene_pseudotime_trends/`: per-method gene expression vs pseudotime trend plots with observed bin points plus a display-layer smoother, phenotype-stratified trend summaries, paired CSV/TSV bundles, prompts, and queued LLM analyses" else NULL,
  if (isTRUE(RUN_PHENOTYPE_PSEUDOTIME_COMPARISON)) "- `phenotype_pseudotime_comparison/`: disease phenotype pseudotime summaries across Slingshot/Monocle methods with paired CSV and LLM companion" else NULL,
  "- `phenotype_column_candidates.tsv`: metadata columns scored for disease phenotype detection",
  "- `figure_llm_manifest_index.csv`: all companion manifests collected for the decoupled LLM runner",
  "- `figure_llm_batch_launch_status.json`: whether the independent LLM runner was started and where its log is"
)
writeLines(readme_lines, con = file.path(OUTPUT_DIR, "README.md"), useBytes = TRUE)

cat("\nLineage summary:\n")
print(lineage_summary)
cat("\nSaved outputs to: ", OUTPUT_DIR, "\n", sep = "")
cat("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("=", 80), "\n", sep = "")