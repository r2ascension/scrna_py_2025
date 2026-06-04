#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison v2.6.11 — mid-weight CHOIR probe
# ==============================================================================
# Date: 2026-05-14
# Purpose:
#   Probe whether a modestly stronger CHOIR configuration yields more than the
#   overly coarse 2-cluster lightweight result, without re-running upstream
#   pairwise DE / ssGSEA / OFA / final export.
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)
source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R")

output_dir <- "/home/h2048/data/R/0514/bcell_tissue_comparison_v2_6_11_c22_c13_c25_c14drop_choir_probe_20260514"
if (dir.exists(output_dir)) {
  unlink(list.files(output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.11-c22-c13-c25-c14drop-choir-probe-mid",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.11",
      "(20260514 mid-weight CHOIR probe with slightly stronger root/subtree settings)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_11_c22_c13_c25_c14drop_choir_probe_20260514.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON PROBE COMPLETE (v2.6.11 mid-weight CHOIR probe)",
    H5AD_PATH = "/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad",
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = FALSE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    REQUIRE_LLM = FALSE,
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CLUSTERING = TRUE,
    RUN_CHOIR = TRUE,
    CLUSTER_BACKEND = "CHOIR",
    RUN_CLUSTER_SSGSEA = FALSE,
    RUN_OFA = FALSE,
    RUN_L3_OFA = FALSE,
    RUN_L3_OFA_VS_REST = FALSE,
    RUN_L3_OFA_SAME_L2 = FALSE,
    RUN_L3_OFA_INTER_TISSUE = FALSE,
    EXIT_AFTER_L3_OFA = TRUE,
    RUN_MILOPY = FALSE,
    N_CORES = 2L,
    CHOIR_N_CORES = 2L,
    CHOIR_ALPHA = 0.24,
    CHOIR_N_ITERATIONS = 35L,
    CHOIR_N_TREES = 30L,
    CHOIR_VAR_FEATURES_MAX = 3000L,
    CHOIR_SAMPLE_MAX = 7000L,
    CHOIR_DOWNSAMPLING_RATE = 0.08,
    CHOIR_MAX_CLUSTERS = "auto",
    CHOIR_MIN_CLUSTER_DEPTH = 1800L,
    CHOIR_MAX_REPEAT_ERRORS = 3L,
    CHOIR_SUBTREE_REDUCTIONS = FALSE,
    CHOIR_DISTANCE_AWARENESS = 1L,
    PIPELINE_CHANGELOG_LINES = c(
      "  [PROBE-1] Use a separate 20260514 probe output directory; do not overwrite completed 0508 results.",
      "  [PROBE-2] Disable pairwise DE / ssGSEA / OFA / final export and stop after clustering-only evaluation.",
      "  [PROBE-3] Keep CHOIR on max_clusters='auto' root/subtree path while modestly increasing alpha, iterations, trees, var_features, sample_max, and downsampling.",
      "  [PROBE-4] Lower min_cluster_depth from the ultra-light setting to 1800 for a slightly finer subtree split."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
