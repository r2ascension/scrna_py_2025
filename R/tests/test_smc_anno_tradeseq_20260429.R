#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

input_rdata <- "/home/h2048/temp/smc_anno.Rdata"
assert_true(file.exists(input_rdata), sprintf("Missing test input: %s", input_rdata))

tmp_dir <- file.path(tempdir(), "smc_anno_tradeseq_20260429")
output_dir <- file.path(tmp_dir, "output")
unlink(output_dir, recursive = TRUE, force = TRUE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

script_env <- new.env(parent = globalenv())
script_env$INPUT_RDATA <- input_rdata
script_env$OBJECT_NAME <- "smc_anno"
script_env$OUTPUT_DIR <- output_dir
script_env$RUN_TRADESEQ <- TRUE
script_env$RUN_MONOCLE3 <- FALSE
script_env$RUN_MONOCLE2 <- FALSE
script_env$RUN_PAGA <- FALSE
script_env$RUN_CROSS_METHOD_SUMMARY <- FALSE
script_env$RUN_FIGURE_LLM_BATCH <- FALSE
script_env$TRADESEQ_MAX_GENES <- 50L
script_env$TRADESEQ_EVALUATEK_N_GENES <- 20L
script_env$TRADESEQ_K_VALUES <- 3:5
script_env$TRADESEQ_NK <- 4L
script_env$TRADESEQ_MIN_CELLS_EXPRESSED <- 10L
script_env$TRADESEQ_ENABLE_PLOTS <- TRUE

sys.source("/home/h2048/script/R/smc_anno_slingshot_20260428.R", envir = script_env)

expected_files <- c(
  "slingshot_pseudotime.tsv",
  file.path("tradeSeq", "tradeSeq_metadata.tsv"),
  file.path("tradeSeq", "tradeSeq_gene_panel.tsv"),
  file.path("tradeSeq", "tradeSeq_evaluateK_metrics.tsv"),
  file.path("tradeSeq", "tradeSeq_associationTest.tsv"),
  file.path("tradeSeq", "tradeSeq_patternTest.tsv"),
  file.path("tradeSeq", "tradeSeq_startVsEndTest.tsv"),
  file.path("tradeSeq", "tradeSeq_diffEndTest.tsv"),
  file.path("tradeSeq", "tradeSeq_evaluateK_metrics.pdf"),
  file.path("tradeSeq", "tradeSeq_evaluateK_metrics_data.csv"),
  file.path("tradeSeq", "tradeSeq_evaluateK_metrics_LLM_ANALYSIS.md"),
  file.path("tradeSeq", "tradeSeq_top_test_results.tsv"),
  file.path("tradeSeq", "tradeSeq_top_test_results.pdf"),
  file.path("tradeSeq", "tradeSeq_top_test_results_data.csv"),
  file.path("tradeSeq", "tradeSeq_top_test_results_LLM_ANALYSIS.md"),
  file.path("tradeSeq", "tradeSeq_test_gene_wald_effect_heatmap.pdf"),
  file.path("tradeSeq", "tradeSeq_test_gene_wald_effect_heatmap_data.csv"),
  file.path("tradeSeq", "tradeSeq_test_gene_wald_effect_heatmap_LLM_ANALYSIS.md"),
  file.path("tradeSeq", "tradeSeq_top_gene_smoother_summary.tsv"),
  file.path("tradeSeq", "tradeSeq_smoother_delta_heatmap.pdf"),
  file.path("tradeSeq", "tradeSeq_smoother_delta_heatmap_data.csv"),
  file.path("tradeSeq", "tradeSeq_smoother_delta_heatmap_LLM_ANALYSIS.md"),
  file.path("tradeSeq", "tradeSeq_FINAL_SUMMARY.md"),
  file.path("tradeSeq", "tradeSeq_FINAL_SUMMARY_data.csv"),
  file.path("tradeSeq", "tradeSeq_FINAL_SUMMARY_LLM_ANALYSIS.md"),
  file.path("figures", "slingshot_umap_clusters_data.csv"),
  file.path("figures", "slingshot_umap_clusters_LLM_ANALYSIS.md"),
  "figure_llm_manifest_index.csv",
  "figure_llm_batch_launch_status.json"
)

missing_files <- expected_files[!file.exists(file.path(output_dir, expected_files))]
assert_true(length(missing_files) == 0L, paste("Expected tradeSeq outputs are missing:", paste(missing_files, collapse = ", ")))

top_plot_data <- read.csv(file.path(output_dir, "tradeSeq", "tradeSeq_top_test_results_data.csv"), stringsAsFactors = FALSE, check.names = FALSE)
assert_true("rank_metric_name" %in% colnames(top_plot_data), "Top test plot data should include rank metric provenance")
assert_true("global_waldStat" %in% colnames(top_plot_data), "Top test plot data should include global Wald statistic")
assert_true(all(top_plot_data$rank_metric_name %in% c("waldStat_due_to_zero_padj", "-log10_significance")), "Top test plot should explicitly record whether ranking used Wald fallback or p-value significance")

final_summary <- paste(readLines(file.path(output_dir, "tradeSeq", "tradeSeq_FINAL_SUMMARY.md"), warn = FALSE), collapse = "\n")
assert_true(grepl("最后总结", final_summary, fixed = TRUE), "tradeSeq final summary should include a final summary section")

status <- jsonlite::read_json(file.path(output_dir, "tradeSeq", "tradeSeq_top_test_results_LLM_STATUS.json"), simplifyVector = TRUE)
assert_true(identical(status$status, "queued_for_batch"), "Figure LLM should be queued, not run synchronously during plotting")

cat("smc_anno tradeSeq smoke test passed.\n")