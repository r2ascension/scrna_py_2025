#!/usr/bin/env Rscript
# ==============================================================================
# Stromal Endothelial Tissue Comparison v1.1.2 — rm CHOIR/OFA 6/52 full rerun
# ==============================================================================
# Date: 2026-05-08
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_ENDOTHELIAL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.1.2-Endothelial-rmCHOIR6-52-fullscvi",
    PIPELINE_SUBTITLE = paste(
      "Stromal Endothelial Tissue Comparison v1.1.2",
      "(remove user-curated OFA/CHOIR clusters 6/52 and rerun branch-specific scVI/scANVI)"
    ),
    GENERATED_BY_LABEL = "stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508.R",
    LINEAGE_COMPLETION_BANNER = "STROMAL ENDOTHELIAL TISSUE COMPARISON COMPLETE (v1.1.2 rmCHOIR6/52 full rerun)",
    H5AD_PATH = "/home/h2048/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508/endothelial/adata_endothelial_reference_v1_5_branchwise.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    REQUIRE_LLM = TRUE,
    CHOIR_N_CORES = 1L,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "3")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ENDO-RM-1] Preserve prior endothelial c5/c33 cleanup and additionally remove user-curated OFA/CHOIR clusters 6 and 52 at h5ad level.",
      "  [ENDO-RM-2] Retrain branch-specific endothelial scVI/scANVI from the filtered stromal input h5ad.",
      "  [LLM-1] Use 2026-05-08 batch-parallel LLM screening overlay.",
      "  [PLOT-1] Pathway enrichment barplots are generated post hoc from ssGSEA/OFA CSV reports."
    )
  )
)

if (dir.exists(PIPELINE_CONFIG$OUTPUT_DIR)) {
  unlink(list.files(PIPELINE_CONFIG$OUTPUT_DIR, all.files = TRUE, no.. = TRUE, full.names = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(PIPELINE_CONFIG$OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
