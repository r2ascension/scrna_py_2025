#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Tissue Comparison Pipeline v2.6.1-TNK Wrapper
# ==============================================================================
#
# Generic-interface entrypoint:
#   - uses tissue_comparison_generic_wrapper_20260412.R
#   - keeps the T/NK-specific preset thin and declarative
#   - reuses helper-synchronized ssGSEA / CHOIR / discovery-screening logic
#
# Date: 2026-04-13
# ============================================================================== 

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "TNK",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.1-TNK",
    PIPELINE_SUBTITLE = paste(
      "T/NK Tissue Comparison v2.6.1-TNK",
      "(generic wrapper interface)"
    ),
    GENERATED_BY_LABEL = "tnk_tissue_comparison_v2_6_1_20260413.R",
    LINEAGE_COMPLETION_BANNER = "T/NK TISSUE COMPARISON COMPLETE (v2.6.1-TNK)",
    OUTPUT_DIR = "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413",
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] TNK wrapper is a thin preset over the generic wrapper instead of a standalone fork.",
      "  [TNK-1] Shared helper now owns TNK marker panels plus ssGSEA/discovery screening helpers.",
      "  [TNK-2] TNK run inherits additive LLM outputs for ssGSEA review, CHOIR review, and family-level discovery/outlier screening.",
      "  [TNK-3] Existing 0407 TNK outputs are used as the previous-run baseline while writing new outputs to an isolated 0413 directory."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
