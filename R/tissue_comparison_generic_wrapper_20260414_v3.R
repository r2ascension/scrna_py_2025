#!/usr/bin/env Rscript
# ==============================================================================
# Generic Tissue Comparison Wrapper Interface (Myeloid CHOIR Tune Extension)
# ==============================================================================
#
# Purpose:
#   - keep older wrapper/helper versions immutable
#   - extend the 2026-04-14 myeloid wrapper with safer CHOIR controls
#   - point myeloid runs at the v2.6.2 shared engine that restores CHOIR tuning
#
# Date: 2026-04-14
# ==============================================================================

BASE_GENERIC_WRAPPER_PATH_V3 <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414_v2.R"
if (!file.exists(BASE_GENERIC_WRAPPER_PATH_V3)) {
  stop(sprintf("Base myeloid wrapper extension not found: %s", BASE_GENERIC_WRAPPER_PATH_V3))
}
source(BASE_GENERIC_WRAPPER_PATH_V3)

tc_myeloid_wrapper_preset_base_20260414_v3 <- tc_myeloid_wrapper_preset

tc_myeloid_wrapper_preset <- function() {
  base_cfg <- tc_myeloid_wrapper_preset_base_20260414_v3()
  tuned_overrides <- utils::modifyList(
    base_cfg$SHARED_OVERRIDES,
    list(
      CHOIR_N_CORES = 1L,
      CHOIR_SAMPLE_MAX = 5000L,
      CHOIR_DOWNSAMPLING_RATE = 0.05,
      CHOIR_SUBTREE_REDUCTIONS = FALSE,
      LLM_EXTRA_RULES = c(
        base_cfg$SHARED_OVERRIDES$LLM_EXTRA_RULES,
        "CHOIR tuning is intentionally conservative here; prefer stable, reproducible cluster state calls over ultra-fine over-segmentation."
      )
    )
  )

  utils::modifyList(
    base_cfg,
    list(
      PIPELINE_VERSION_LABEL = "v1.2.2-MYELOID",
      PIPELINE_SUBTITLE = paste(
        "Myeloid Tissue Comparison v1.2.2-MYELOID",
        "(generic wrapper interface, tuned CHOIR controls)"
      ),
      GENERATED_BY_LABEL = "myeloid_tissue_comparison_v1_2_2_20260414.R",
      LINEAGE_COMPLETION_BANNER = "MYELOID TISSUE COMPARISON COMPLETE (v1.2.2-MYELOID)",
      OUTPUT_DIR = "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414",
      PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_1_20260414",
      SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_2_20260414.R",
      SHARED_OVERRIDES = tuned_overrides,
      PIPELINE_CHANGELOG_LINES = c(
        base_cfg$PIPELINE_CHANGELOG_LINES,
        "  [MYELOID-5] Switch shared engine to v2.6.2 so wrapper-based runs can pass explicit CHOIR tuning arguments.",
        "  [MYELOID-6] Use conservative CHOIR controls (`sample_max=5000`, `downsampling_rate=0.05`, `subtree_reductions=FALSE`, `CHOIR_N_CORES=1`) to avoid the pathological heavy-default branch observed in v1.2.1."
      )
    )
  )
}

if (sys.nframe() == 0) {
  if (!exists("PIPELINE_CONFIG", inherits = FALSE)) {
    cat("Myeloid CHOIR-tuned wrapper extension loaded.\n")
    cat("Supported extra lineage: MYELOID\n")
  } else {
    tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
  }
}
