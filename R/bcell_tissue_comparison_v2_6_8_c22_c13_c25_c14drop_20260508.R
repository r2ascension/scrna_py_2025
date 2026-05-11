#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison v2.6.8 — c22/c13/c25/c14 full rerun
# ==============================================================================
# Date: 2026-05-08
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.8-c22-c13-c25-c14drop-l3-freshscvi",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.8",
      "(full rerun after removing CHOIR c22 + c13 + c25 + c14; fresh scVI/scANVI)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.8 c22/c13/c25/c14drop)",
    H5AD_PATH = "/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    REQUIRE_LLM = TRUE,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "3")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [DROP-1] Preserve prior B-cell c22/c13/c25 removals and additionally remove user-curated CHOIR c14.",
      "  [DROP-2] Start from the upstream B-cell reference h5ad and retrain fresh scVI/scANVI before all downstream analyses.",
      "  [LLM-1] Use 2026-05-08 batch-parallel LLM screening overlay for discovery/outlier interpretation.",
      "  [PLOT-1] Pathway enrichment barplots are generated as a post-processing artifact under figures/pathway_enrichment_barplots_20260508."
    )
  )
)

if (dir.exists(PIPELINE_CONFIG$OUTPUT_DIR)) {
  unlink(list.files(PIPELINE_CONFIG$OUTPUT_DIR, all.files = TRUE, no.. = TRUE, full.names = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(PIPELINE_CONFIG$OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
