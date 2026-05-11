#!/usr/bin/env Rscript
# ==============================================================================
# Stromal SMC/Pericyte Tissue Comparison v1.1.1-neuronlike Wrapper
# ==============================================================================
#
# Purpose:
#   - rerun stromal SMC/pericyte tissue comparison after relabeling CHOIR cluster 6
#   - conservatively rename the neuron-like outlier state to `Peripheral_neuron_like`
#   - consume the 2026-04-14 SMC branch rerun reference h5ad
#
# Date: 2026-04-14
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_SMC",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.1.1-SMC-neuronlike",
    PIPELINE_SUBTITLE = paste(
      "Stromal SMC/Pericyte Tissue Comparison v1.1.1",
      "(CHOIR cluster 6 relabeled to Peripheral_neuron_like and scANVI rerun)"
    ),
    GENERATED_BY_LABEL = "stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414.R",
    LINEAGE_COMPLETION_BANNER = "STROMAL SMC/PERICYTE TISSUE COMPARISON COMPLETE (v1.1.1-SMC-neuronlike)",
    H5AD_PATH = "/home/h2048/data/py/0414/stromal_branch_rerun_20260414/smc/adata_smc_reference_neuronlike_c6_20260414_scanvi_umap_refresh_20260416.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_0_20260414",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    CHOIR_N_CORES = 1L,
    USE_EXISTING_L2 = FALSE,
    ANALYSIS_L3_DESCRIPTION = "branch-specific scANVI predictions with conservative Peripheral_neuron_like relabel for former CHOIR cluster 6",
    L3_TO_L2_REMAP = c(
      "Muscle_pericyte_pulmonary" = "Pericyte",
      "Muscle_pericyte_systemic" = "Pericyte",
      "Muscle_perivascular_immune_recruiting" = "Smooth_Muscle",
      "Muscle_smooth_pulmonary" = "Smooth_Muscle",
      "Muscle_smooth_arterial_systemic" = "Smooth_Muscle",
      "Peripheral_neuron_like" = "Peripheral_neuron_like"
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [SMC-NL-1] Reinterpret CHOIR cluster 6 as a conservative neuron-like mural state and relabel it to `Peripheral_neuron_like`.",
      "  [SMC-NL-2] Keep the naming intentionally conservative: evidence supports synaptic / peripheral-neuron-like features, but not a confident autonomic, sensory, or Schwann identity.",
      "  [SMC-NL-3] Rerun SMC scANVI on the relabeled reference and derive L2 explicitly from rerun L3 labels instead of reusing stale input `cell_type_L2`.",
      "  [SMC-NL-4] 0414 v1.1.0 output remains the historical baseline only; this wrapper forces a fresh h5ad load from the rerun artifact.",
      "  [SMC-NL-5] Wrapper still requires a valid DEEPSEEK key for full LLM completion."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
