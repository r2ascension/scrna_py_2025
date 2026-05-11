#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

input_rdata <- "/home/h2048/temp/smc_anno.Rdata"
assert_true(file.exists(input_rdata), sprintf("Missing test input: %s", input_rdata))

tmp_dir <- file.path(tempdir(), "smc_anno_cross_methods_20260430")
output_dir <- file.path(tmp_dir, "output")
unlink(output_dir, recursive = TRUE, force = TRUE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

script_env <- new.env(parent = globalenv())
script_env$INPUT_RDATA <- input_rdata
script_env$OBJECT_NAME <- "smc_anno"
script_env$OUTPUT_DIR <- output_dir
script_env$RUN_TRADESEQ <- FALSE
script_env$RUN_MONOCLE3 <- TRUE
script_env$RUN_MONOCLE2 <- TRUE
script_env$RUN_PAGA <- TRUE
script_env$RUN_CROSS_METHOD_SUMMARY <- TRUE
script_env$RUN_FIGURE_LLM_BATCH <- FALSE
script_env$PAGA_PYTHON_BIN <- "/home/h2048/miniconda3/envs/bbknn_env/bin/python"
script_env$MONOCLE2_MAX_ORDERING_GENES <- 300L
script_env$MONOCLE3_NUM_DIM <- 20L

sys.source("/home/h2048/script/R/smc_anno_slingshot_20260428.R", envir = script_env)

expected_files <- c(
  file.path("monocle3", "monocle3_pseudotime.tsv"),
  file.path("monocle3", "monocle3_trajectory_celltypes_data.csv"),
  file.path("monocle3", "monocle3_trajectory_celltypes_LLM_ANALYSIS.md"),
  file.path("monocle3", "monocle3_trajectory_pseudotime_data.csv"),
  file.path("monocle3", "monocle3_trajectory_pseudotime_LLM_ANALYSIS.md"),
  file.path("monocle2", "monocle2_pseudotime.tsv"),
  file.path("monocle2", "monocle2_trajectory_celltypes_data.csv"),
  file.path("monocle2", "monocle2_trajectory_celltypes_LLM_ANALYSIS.md"),
  file.path("monocle2", "monocle2_trajectory_pseudotime_data.csv"),
  file.path("monocle2", "monocle2_trajectory_pseudotime_LLM_ANALYSIS.md"),
  file.path("paga", "paga_connectivities.csv"),
  file.path("paga", "paga_connectivity_heatmap.pdf"),
  file.path("paga", "paga_connectivity_heatmap_data.csv"),
  file.path("paga", "paga_connectivity_heatmap_LLM_ANALYSIS.md"),
  file.path("paga", "paga_network_graph.pdf"),
  file.path("paga", "paga_network_graph_data.csv"),
  file.path("paga", "paga_network_graph_LLM_ANALYSIS.md"),
  file.path("paga", "paga_embedding_groups.pdf"),
  file.path("paga", "paga_embedding_groups_data.csv"),
  file.path("paga", "paga_embedding_groups_LLM_ANALYSIS.md"),
  file.path("figures", "slingshot_umap_clusters_data.csv"),
  file.path("figures", "slingshot_umap_clusters_LLM_ANALYSIS.md"),
  file.path("gene_pseudotime_trends", "gene_trend_panel.tsv"),
  file.path("gene_pseudotime_trends", "Slingshot", "gene_pseudotime_trend_plot_data.tsv"),
  file.path("gene_pseudotime_trends", "Monocle3", "gene_pseudotime_trend_plot_data.tsv"),
  file.path("gene_pseudotime_trends", "Monocle2", "gene_pseudotime_trend_plot_data.tsv"),
  "figure_llm_manifest_index.csv",
  "figure_llm_batch_launch_status.json",
  "cross_method_trajectory_summary.tsv"
)

missing_files <- expected_files[!file.exists(file.path(output_dir, expected_files))]
assert_true(length(missing_files) == 0L, paste("Expected cross-method outputs are missing:", paste(missing_files, collapse = ", ")))

status <- jsonlite::read_json(file.path(output_dir, "paga", "paga_network_graph_LLM_STATUS.json"), simplifyVector = TRUE)
assert_true(identical(status$status, "queued_for_batch"), "PAGA figure LLM should be queued, not run synchronously during plotting")

trend_plot_data <- read.delim(
  file.path(output_dir, "gene_pseudotime_trends", "Slingshot", "gene_pseudotime_trend_plot_data.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert_true("plot_layer" %in% colnames(trend_plot_data), "Trend plot data must contain plot_layer")
assert_true("plot_expression" %in% colnames(trend_plot_data), "Trend plot data must contain plot_expression")
assert_true(any(trend_plot_data$plot_layer == "observed_bin"), "Trend plot data must preserve observed bin rows")
assert_true(any(trend_plot_data$plot_layer == "smoothed_curve"), "Trend plot data must contain smoothed curve rows")

readme_lines <- readLines(file.path(output_dir, "README.md"), warn = FALSE)
assert_true(
  any(grepl("display-layer smoother", readme_lines, fixed = TRUE)),
  "README must explain that non-tradeSeq trend smoothing is display-layer only"
)

cat("smc_anno cross-method trajectory smoke test passed.\n")