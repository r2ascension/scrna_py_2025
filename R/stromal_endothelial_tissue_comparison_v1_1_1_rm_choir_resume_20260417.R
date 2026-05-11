#!/usr/bin/env Rscript
# ==============================================================================
# Stromal Endothelial Tissue Comparison v1.1.1-rmCHOIR Resume Wrapper
# ==============================================================================
#
# Purpose:
#   - resume the stalled endothelial rerun from cached stage artifacts already
#     written into the 0414 v1.1.1 output directory
#   - skip visualization / pairwise DE / grouped ssGSEA, which already finished
#   - rerun CHOIR/OFA/report/export with single-core CHOIR for stability
#
# Date: 2026-04-17
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_ENDOTHELIAL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.1.1-Endothelial-rmCHOIR",
    PIPELINE_SUBTITLE = paste(
      "Stromal Endothelial Tissue Comparison v1.1.1",
      "(resume cached stages and rerun CHOIR/OFA after multi-core stall)"
    ),
    GENERATED_BY_LABEL = "stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_resume_20260417.R",
    LINEAGE_COMPLETION_BANNER = "STROMAL ENDOTHELIAL TISSUE COMPARISON COMPLETE (v1.1.1-Endothelial-rmCHOIR)",
    H5AD_PATH = "/home/h2048/data/py/0414/stromal_branch_rerun_20260414/endothelial/adata_endothelial_reference_rm_choir_5_33_20260414_scanvi_umap_refresh_20260416.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CHOIR = TRUE,
    RUN_OFA = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    CHOIR_N_CORES = 1L,
    PIPELINE_CHANGELOG_LINES = c(
      "  [ENDO-RM-1] Remove CHOIR clusters 5 and 33 from the 0414 endothelial tissue-comparison result, then rebuild the branch reference.",
      "  [ENDO-RM-2] Retrain endothelial scVI + scANVI from the filtered reference; downstream L3 still comes from `cell_type_scanvi_pred`.",
      "  [ENDO-RM-3] 0414 v1.1.0 output remains the historical baseline only; this wrapper forces a fresh h5ad load from the rerun artifact.",
      "  [ENDO-RM-4] Wrapper still requires a valid DEEPSEEK key for full LLM completion.",
      "  [ENDO-RM-5] Resume from cached v1.1.1 stage artifacts and rerun CHOIR/OFA with CHOIR_N_CORES=1 after the observed multi-core prune stall."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
