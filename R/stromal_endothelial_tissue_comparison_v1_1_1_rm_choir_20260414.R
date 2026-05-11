#!/usr/bin/env Rscript
# ==============================================================================
# Stromal Endothelial Tissue Comparison v1.1.1-rmCHOIR Wrapper
# ==============================================================================
#
# Purpose:
#   - rerun stromal endothelial tissue comparison after removing CHOIR clusters 5/33
#   - consume the 2026-04-14 endothelial branch rerun reference h5ad
#   - keep the generic stromal wrapper architecture unchanged
#
# Date: 2026-04-14
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_ENDOTHELIAL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.1.1-Endothelial-rmCHOIR",
    PIPELINE_SUBTITLE = paste(
      "Stromal Endothelial Tissue Comparison v1.1.1",
      "(remove CHOIR clusters 5/33 and rerun branch-specific scVI/scANVI)"
    ),
    GENERATED_BY_LABEL = "stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414.R",
    LINEAGE_COMPLETION_BANNER = "STROMAL ENDOTHELIAL TISSUE COMPARISON COMPLETE (v1.1.1-Endothelial-rmCHOIR)",
    H5AD_PATH = "/home/h2048/data/py/0414/stromal_branch_rerun_20260414/endothelial/adata_endothelial_reference_rm_choir_5_33_20260414_scanvi_umap_refresh_20260416.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    CHOIR_N_CORES = 1L,
    PIPELINE_CHANGELOG_LINES = c(
      "  [ENDO-RM-1] Remove CHOIR clusters 5 and 33 from the 0414 endothelial tissue-comparison result, then rebuild the branch reference.",
      "  [ENDO-RM-2] Retrain endothelial scVI + scANVI from the filtered reference; downstream L3 still comes from `cell_type_scanvi_pred`.",
      "  [ENDO-RM-3] 0414 v1.1.0 output remains the historical baseline only; this wrapper forces a fresh h5ad load from the rerun artifact.",
      "  [ENDO-RM-4] Wrapper still requires a valid DEEPSEEK key for full LLM completion."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
