#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison v2.6.10 — lightweight CHOIR recovery
# ==============================================================================
# Date: 2026-05-13
# Purpose:
#   Restart the long-running 20260511 B-cell CHOIR recovery with lower-cost CHOIR
#   settings. Reuse cached pairwise DE/enrichment/LLM/ssGSEA artifacts from the
#   current 0508 output directory and only re-run CHOIR/OFA/finalization.
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)
source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R")
source("/home/h2048/script/R/tissue_comparison_recovery_lock_20260512.R")

output_dir <- "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508"

cleanup_cluster_artifacts <- function(output_dir, backend = "choir") {
  backend <- tolower(backend)
  report_dir <- file.path(output_dir, "reports")
  figure_dir <- file.path(output_dir, "figures")
  unlink(c(
    file.path(report_dir, backend),
    file.path(report_dir, sprintf("%s_llm_structured.tsv", backend)),
    file.path(report_dir, sprintf("%s_llm_structured_all.rds", backend)),
    file.path(report_dir, sprintf("llm_%s_discovery_screen.tsv", backend)),
    file.path(output_dir, sprintf("LLM_%s_INTERPRETATION.md", toupper(backend))),
    file.path(output_dir, sprintf("LLM_%s_DISCOVERY_REVIEW.md", toupper(backend)))
  ), recursive = TRUE, force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, sprintf("%s_umap*", backend))), force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, sprintf("ssgsea_%s_heatmap*", backend))), force = TRUE)
  unlink(c(
    file.path(output_dir, "bcell_tissue_comparison_final.rds"),
    file.path(output_dir, "bcell_tissue_comparison_final.h5ad"),
    file.path(output_dir, "REPORT.md")
  ), force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.recovery_lock_release <- tc_recovery_register_lock(
  output_dir = output_dir,
  lock_name = ".bcell_recovery_20260513_light.lock",
  completion_paths = c(
    file.path(output_dir, "bcell_tissue_comparison_final.rds"),
    file.path(output_dir, "bcell_tissue_comparison_final.h5ad"),
    file.path(output_dir, "REPORT.md"),
    file.path(output_dir, "reports", "choir_llm_structured.tsv")
  ),
  label = "B-cell lightweight CHOIR recovery"
)
on.exit(.recovery_lock_release(), add = TRUE)
cleanup_cluster_artifacts(output_dir, backend = "choir")

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.10-c22-c13-c25-c14drop-l3-freshscvi-choir-light",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.10",
      "(20260513 lightweight CHOIR restart from cached 0508 pairwise/ssGSEA artifacts)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_10_c22_c13_c25_c14drop_choir_light_20260513.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.10 lightweight CHOIR recovery)",
    H5AD_PATH = "/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad",
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    REQUIRE_LLM = TRUE,
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CLUSTERING = TRUE,
    RUN_CHOIR = TRUE,
    CLUSTER_BACKEND = "CHOIR",
    RUN_CLUSTER_SSGSEA = TRUE,
    RUN_OFA = TRUE,
    RUN_L3_OFA = TRUE,
    RUN_L3_OFA_VS_REST = TRUE,
    RUN_L3_OFA_SAME_L2 = FALSE,
    RUN_L3_OFA_INTER_TISSUE = TRUE,
    N_CORES = 2L,
    CHOIR_N_CORES = 2L,
    CHOIR_N_ITERATIONS = 25L,
    CHOIR_N_TREES = 20L,
    CHOIR_VAR_FEATURES_MAX = 2000L,
    CHOIR_SAMPLE_MAX = 5000L,
    CHOIR_DOWNSAMPLING_RATE = 0.20,
    CHOIR_MAX_CLUSTERS = 20L,
    CHOIR_MIN_CLUSTER_DEPTH = 2500L,
    CHOIR_MAX_REPEAT_ERRORS = 2L,
    CHOIR_SUBTREE_REDUCTIONS = FALSE,
    CHOIR_DISTANCE_AWARENESS = 1L,
    OFA_MAX_CELLS_PER_IDENT = 2000L,
    L3_OFA_MAX_CELLS_PER_IDENT = 2000L,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "2")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [LIGHT-1] Restart CHOIR with cached pairwise DE/enrichment/LLM/ssGSEA artifacts; no re-run of upstream pairwise stages.",
      "  [LIGHT-2] Reduce CHOIR permutations 100->25, RF trees 50->20, variable features 4000->2000, and cores 4->2.",
      "  [LIGHT-3] Limit CHOIR sampled cells/downsampling and cap max clusters to 20 to avoid prolonged underclustering pruning.",
      "  [LIGHT-4] Disable subtree dimensionality reductions and increase minimum cluster depth for a shallower tree.",
      "  [LLM-1] Keep the 20260508 batch-parallel LLM screening overlay."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
