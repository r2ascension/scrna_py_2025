#!/usr/bin/env Rscript
# ==============================================================================
# Stromal SMC Slingshot Trajectory Analysis (2026-04-28)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
})

TRAJECTORY_BRANCH_HELPER_PATH <- "/home/h2048/script/R/trajectory_branch_helper_20260428_v1.R"
if (!exists("pa_build_slingshot_trajectory_packet", mode = "function")) {
  source(TRAJECTORY_BRANCH_HELPER_PATH)
}
TRAJECTORY_INTERPRET_HELPER_PATH <- "/home/h2048/script/R/stromal_smc_trajectory_interpret_helper_20260429.R"
if (!exists("smcti_run_slingshot_analysis", mode = "function")) {
  source(TRAJECTORY_INTERPRET_HELPER_PATH)
}

if (!exists("INPUT_RDS")) {
  INPUT_RDS <- "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.rds"
}
if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- "/home/h2048/data/R/0428/stromal_smc_slingshot_20260428"
}
if (!exists("ASSAY_NAME")) {
  ASSAY_NAME <- "RNA"
}
if (!exists("REDUCTION_NAME")) {
  REDUCTION_NAME <- "umap"
}
if (!exists("CLUSTER_COL")) {
  CLUSTER_COL <- "cell_type_L3"
}
if (!exists("START_CLUSTER")) {
  START_CLUSTER <- "Muscle_pericyte_pulmonary"
}
if (!exists("END_CLUSTERS")) {
  END_CLUSTERS <- c("Muscle_smooth_pulmonary", "Muscle_smooth_arterial_systemic")
}
if (!exists("SAVE_SCE")) {
  SAVE_SCE <- TRUE
}
if (!exists("RUN_QC")) {
  RUN_QC <- TRUE
}
if (!exists("RUN_SENSITIVITY")) {
  RUN_SENSITIVITY <- TRUE
}
if (!exists("SENSITIVITY_REDUCTIONS")) {
  SENSITIVITY_REDUCTIONS <- c("umap", "pca", "harmony")
}
if (!exists("RUN_DYNAMIC_GENES")) {
  RUN_DYNAMIC_GENES <- TRUE
}
if (!exists("RUN_PATHWAY_DYNAMICS")) {
  RUN_PATHWAY_DYNAMICS <- TRUE
}
if (!exists("RUN_BRANCH_ASSOCIATION")) {
  RUN_BRANCH_ASSOCIATION <- TRUE
}
if (!exists("RUN_MONOCLE3")) {
  RUN_MONOCLE3 <- TRUE
}
if (!exists("RUN_MONOCLE2")) {
  RUN_MONOCLE2 <- TRUE
}
if (!exists("RUN_PAGA")) {
  RUN_PAGA <- TRUE
}
if (!exists("RUN_CROSS_METHOD_SUMMARY")) {
  RUN_CROSS_METHOD_SUMMARY <- TRUE
}
if (!exists("SENSITIVITY_WEIGHT_THRESHOLD")) {
  SENSITIVITY_WEIGHT_THRESHOLD <- 0.5
}
if (!exists("AMBIGUITY_WEIGHT_THRESHOLD")) {
  AMBIGUITY_WEIGHT_THRESHOLD <- 0.2
}
if (!exists("PAGA_GROUP_KEY")) {
  PAGA_GROUP_KEY <- CLUSTER_COL
}
if (!exists("PAGA_PYTHON_BIN")) {
  PAGA_PYTHON_BIN <- "/home/h2048/miniconda3/envs/bbknn_env/bin/python"
}
if (!exists("MONOCLE3_NUM_DIM")) {
  MONOCLE3_NUM_DIM <- 30L
}
if (!exists("MONOCLE2_MAX_ORDERING_GENES")) {
  MONOCLE2_MAX_ORDERING_GENES <- 2000L
}
if (!exists("DYNAMIC_MAX_GENES")) {
  DYNAMIC_MAX_GENES <- 2500L
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", as.character(x))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
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
cat("STROMAL SMC SLINGSHOT TRAJECTORY ANALYSIS\n")
cat(strrep("=", 80), "\n", sep = "")
cat("Started   : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Input RDS : ", INPUT_RDS, "\n", sep = "")
cat("Output dir: ", OUTPUT_DIR, "\n", sep = "")
cat("Reduction : ", REDUCTION_NAME, "\n", sep = "")
cat("Clusters  : ", CLUSTER_COL, "\n", sep = "")
cat("Start     : ", START_CLUSTER, "\n", sep = "")
cat("End       : ", paste(END_CLUSTERS, collapse = ", "), "\n\n", sep = "")

if (!file.exists(INPUT_RDS)) {
  stop(sprintf("Input RDS does not exist: %s", INPUT_RDS), call. = FALSE)
}

obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) {
  stop("INPUT_RDS is not a Seurat object", call. = FALSE)
}

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
  stop(sprintf("Metadata column '%s' contains empty or NA labels", CLUSTER_COL), call. = FALSE)
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

cat("Loaded Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features\n", sep = "")
cat("Using assay: ", ASSAY_NAME, "\n", sep = "")
cat("Cluster counts:\n")
print(sort(table(cluster_labels), decreasing = TRUE))

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
  c("cell_type_L2", "cell_type_L3", "tissue", "sample", "group", "condition", "seurat_clusters"),
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
    out[[paste0(lineage_name, "_median_pseudotime")]] <- safe_numeric_summary(
      pseudotime_mat[idx, lineage_name],
      stats::median
    )
    out[[paste0(lineage_name, "_mean_curve_weight")]] <- safe_numeric_summary(
      curve_weight_mat[idx, lineage_name],
      mean
    )
    out[[paste0(lineage_name, "_cells_weight_ge_0_5")]] <- sum(
      curve_weight_mat[idx, lineage_name] >= 0.5,
      na.rm = TRUE
    )
  }
  out
})
cluster_summary <- do.call(rbind, cluster_summary_list)

trajectory_packet <- pa_build_slingshot_trajectory_packet(
  pseudotime_table = pseudotime_table,
  lineage_summary = lineage_summary,
  branch_summary = cluster_summary,
  source_path = normalizePath(INPUT_RDS, winslash = "/", mustWork = FALSE),
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

plot_slingshot_clusters(
  sce = sce,
  labels = cluster_labels,
  output_path = file.path(FIG_DIR, "slingshot_umap_clusters.pdf"),
  title_text = sprintf("Stromal SMC Slingshot (%s on %s)", CLUSTER_COL, REDUCTION_NAME)
)

for (lineage_name in colnames(pseudotime_mat)) {
  lineage_row <- lineage_summary[lineage_summary$lineage_id == lineage_name, , drop = FALSE]
  plot_slingshot_pseudotime(
    sce = sce,
    pseudotime_vec = pseudotime_mat[, lineage_name],
    output_path = file.path(FIG_DIR, sprintf("slingshot_%s_pseudotime.pdf", sanitize_name(lineage_name))),
    title_text = sprintf(
      "%s: %s",
      lineage_name,
      lineage_row$lineage_path[[1]]
    )
  )
}

main_res <- list(
  sce = sce,
  cluster_labels = cluster_labels,
  pseudotime_mat = pseudotime_mat,
  curve_weight_mat = curve_weight_mat,
  lineages = lineages,
  pseudotime_table = pseudotime_table,
  lineage_summary = lineage_summary,
  cluster_summary = cluster_summary,
  primary_terminal_state = smcti_primary_terminal_state(pseudotime_table, lineage_summary),
  reduction_name = REDUCTION_NAME,
  cluster_col = CLUSTER_COL,
  start_cluster = START_CLUSTER,
  end_clusters = as.character(END_CLUSTERS)
)

qc_bundle <- NULL
association_bundle <- NULL
sensitivity_summary <- data.frame()
dynamic_bundle <- NULL
monocle3_res <- list(status = "skipped", error = "RUN_MONOCLE3=FALSE")
monocle2_res <- list(status = "skipped", error = "RUN_MONOCLE2=FALSE")
paga_res <- list(status = "skipped", error = "RUN_PAGA=FALSE")
cross_method_summary <- data.frame()

if (!exists("PAGA_USE_REDUCTION") || is.null(PAGA_USE_REDUCTION) || !nzchar(PAGA_USE_REDUCTION) || !PAGA_USE_REDUCTION %in% names(obj@reductions)) {
  PAGA_USE_REDUCTION <- if ("harmony" %in% names(obj@reductions)) {
    "harmony"
  } else if ("pca" %in% names(obj@reductions)) {
    "pca"
  } else {
    REDUCTION_NAME
  }
}

if (isTRUE(RUN_QC)) {
  cat("\n[QC] Writing trajectory reliability summaries...\n")
  qc_bundle <- smcti_write_qc_bundle(
    res = main_res,
    output_dir = OUTPUT_DIR,
    fig_dir = FIG_DIR,
    context_cols = c("sample", "tissue"),
    ambiguity_threshold = AMBIGUITY_WEIGHT_THRESHOLD,
    weight_threshold = SENSITIVITY_WEIGHT_THRESHOLD
  )
}

if (isTRUE(RUN_BRANCH_ASSOCIATION)) {
  cat("\n[ASSOCIATION] Summarizing sample/tissue branch composition...\n")
  association_bundle <- smcti_write_association_bundle(
    res = main_res,
    output_dir = file.path(OUTPUT_DIR, "association"),
    fig_dir = FIG_DIR,
    sample_col = "sample",
    tissue_col = "tissue"
  )
}

if (isTRUE(RUN_SENSITIVITY)) {
  cat("\n[SENSITIVITY] Comparing alternative reductions and terminal settings...\n")
  sensitivity_summary <- smcti_write_sensitivity_bundle(
    seurat_obj = obj,
    reference_res = main_res,
    output_dir = file.path(OUTPUT_DIR, "sensitivity"),
    assay_name = ASSAY_NAME,
    cluster_col = CLUSTER_COL,
    start_cluster = START_CLUSTER,
    end_clusters = as.character(END_CLUSTERS),
    reductions = SENSITIVITY_REDUCTIONS,
    meta_cols = meta_cols,
    weight_threshold = SENSITIVITY_WEIGHT_THRESHOLD
  )
}

if (isTRUE(RUN_DYNAMIC_GENES)) {
  cat("\n[DYNAMIC] Running fallback dynamic gene / pathway interpretation...\n")
  dynamic_bundle <- smcti_run_dynamic_gene_bundle(
    seurat_obj = obj,
    slingshot_res = main_res,
    output_dir = file.path(OUTPUT_DIR, "dynamic_genes"),
    assay_name = ASSAY_NAME,
    run_pathway = RUN_PATHWAY_DYNAMICS,
    max_genes = DYNAMIC_MAX_GENES
  )
}

if (isTRUE(RUN_MONOCLE3)) {
  cat("\n[MONOCLE3] Running monocle3 cross-check...\n")
  if (requireNamespace("monocle3", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("igraph", quietly = TRUE)) {
    monocle3_res <- smcti_run_monocle3_bundle(
      seurat_obj = obj,
      output_dir = file.path(OUTPUT_DIR, "monocle3"),
      assay_name = ASSAY_NAME,
      cluster_col = CLUSTER_COL,
      start_cluster = START_CLUSTER,
      end_clusters = as.character(END_CLUSTERS),
      reduction_name = REDUCTION_NAME,
      num_dim = MONOCLE3_NUM_DIM
    )
  } else {
    monocle3_res <- list(status = "skipped", error = "monocle3/ggplot2/igraph unavailable")
    writeLines(monocle3_res$error, con = file.path(ensure_dir(file.path(OUTPUT_DIR, "monocle3")), "monocle3_error.txt"), useBytes = TRUE)
  }
}

if (isTRUE(RUN_MONOCLE2)) {
  cat("\n[MONOCLE2] Running monocle2 cross-check...\n")
  if (requireNamespace("monocle", quietly = TRUE)) {
    monocle2_res <- smcti_run_monocle2_bundle(
      seurat_obj = obj,
      output_dir = file.path(OUTPUT_DIR, "monocle2"),
      assay_name = ASSAY_NAME,
      cluster_col = CLUSTER_COL,
      start_cluster = START_CLUSTER,
      end_clusters = as.character(END_CLUSTERS),
      max_ordering_genes = MONOCLE2_MAX_ORDERING_GENES,
      marker_panel = smcti_dynamic_gene_marker_panel
    )
  } else {
    monocle2_res <- list(status = "skipped", error = "monocle unavailable")
    writeLines(monocle2_res$error, con = file.path(ensure_dir(file.path(OUTPUT_DIR, "monocle2")), "monocle2_error.txt"), useBytes = TRUE)
  }
}

if (isTRUE(RUN_PAGA)) {
  cat("\n[PAGA] Running graph abstraction cross-check...\n")
  if (nzchar(Sys.which(PAGA_PYTHON_BIN)) || file.exists(PAGA_PYTHON_BIN)) {
    paga_res <- smcti_run_paga_bundle(
      seurat_obj = obj,
      output_dir = file.path(OUTPUT_DIR, "paga"),
      groupby_key = PAGA_GROUP_KEY,
      reduction_name = PAGA_USE_REDUCTION,
      assay_name = ASSAY_NAME,
      python_bin = PAGA_PYTHON_BIN
    )
    if (identical(paga_res$status, "ok")) {
      trajectory_branch_packet <- pa_build_trajectory_branch_packet(
        topology_screen = paga_res$topology_packet,
        primary_trajectory = trajectory_packet,
        preferred_engine = "Slingshot"
      )
      saveRDS(trajectory_branch_packet, file.path(OUTPUT_DIR, "trajectory_branch_packet.rds"))
    }
  } else {
    paga_res <- list(status = "skipped", error = sprintf("PAGA python not found: %s", PAGA_PYTHON_BIN))
    writeLines(paga_res$error, con = file.path(ensure_dir(file.path(OUTPUT_DIR, "paga")), "paga_error.txt"), useBytes = TRUE)
  }
}

if (isTRUE(RUN_CROSS_METHOD_SUMMARY)) {
  cross_method_summary <- smcti_build_cross_method_summary(
    slingshot_res = main_res,
    monocle3_res = monocle3_res,
    monocle2_res = monocle2_res,
    paga_res = paga_res
  )
  write_tsv(cross_method_summary, file.path(OUTPUT_DIR, "cross_method_trajectory_summary.tsv"))
}

module_status_lines <- c(
  "## Interpretation modules",
  sprintf("- QC: `%s`", if (isTRUE(RUN_QC)) "completed" else "skipped"),
  sprintf("- Sensitivity: `%s`", if (isTRUE(RUN_SENSITIVITY)) "completed" else "skipped"),
  sprintf("- Dynamic genes: `%s`", if (isTRUE(RUN_DYNAMIC_GENES)) ifelse(is.null(dynamic_bundle), "skipped", dynamic_bundle$method) else "skipped"),
  sprintf("- Association: `%s`", if (isTRUE(RUN_BRANCH_ASSOCIATION)) "completed" else "skipped"),
  sprintf("- Monocle3: `%s`", monocle3_res$status),
  sprintf("- Monocle2: `%s`", monocle2_res$status),
  sprintf("- PAGA: `%s`", paga_res$status)
)
extra_file_lines <- c(
  "- `trajectory_qc_lineage.tsv`: lineage-level QC summary",
  "- `trajectory_qc_branch_overlap.tsv`: branch ambiguity / overlap summary",
  "- `trajectory_qc_lineage_context.tsv`: sample/tissue distribution per lineage",
  "- `sensitivity/sensitivity_summary.tsv`: reduction / terminal sensitivity checks",
  "- `association/sample_lineage_fraction.tsv`: sample-level primary-lineage composition",
  "- `dynamic_genes/`: pseudotime-correlated genes + pathway summaries",
  "- `monocle3/`: monocle3 pseudotime cross-check outputs",
  "- `monocle2/`: monocle2 pseudotime cross-check outputs",
  "- `paga/`: PAGA connectivity graph outputs",
  "- `cross_method_trajectory_summary.tsv`: Slingshot/Monocle/PAGA concordance overview"
)

readme_lines <- c(
  "# Stromal SMC Slingshot trajectory run",
  "",
  sprintf("- Run date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Input RDS: `%s`", INPUT_RDS),
  sprintf("- Output dir: `%s`", OUTPUT_DIR),
  sprintf("- Assay: `%s`", ASSAY_NAME),
  sprintf("- Reduction: `%s`", REDUCTION_NAME),
  sprintf("- Cluster column: `%s`", CLUSTER_COL),
  sprintf("- Start cluster: `%s`", START_CLUSTER),
  sprintf("- End clusters: `%s`", paste(as.character(END_CLUSTERS), collapse = "`, `")),
  "",
  "## Lineages",
  apply(lineage_summary, 1, function(x) sprintf(
    "- %s: %s (n pseudotime cells=%s)",
    x[["lineage_id"]],
    x[["lineage_path"]],
    x[["n_cells_with_pseudotime"]]
  )),
  "",
  module_status_lines,
  "",
  "## Files",
  "- `slingshot_pseudotime.tsv`: per-cell pseudotime + curve weights",
  "- `slingshot_lineage_summary.tsv`: lineage-level summary",
  "- `slingshot_cluster_summary.tsv`: cluster-level pseudotime/weight summary",
  "- `slingshot_trajectory_packet.rds`: primary Slingshot packet",
  "- `trajectory_branch_packet.rds`: branch packet wrapper",
  "- `slingshot_sce.rds`: SingleCellExperiment with Slingshot fit",
  "- `figures/`: UMAP lineage, QC density, and association PDFs",
  extra_file_lines
)
writeLines(readme_lines, con = file.path(OUTPUT_DIR, "README.md"), useBytes = TRUE)

cat("\nLineage summary:\n")
print(lineage_summary)
cat("\nSaved outputs to: ", OUTPUT_DIR, "\n", sep = "")
cat("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("=", 80), "\n", sep = "")