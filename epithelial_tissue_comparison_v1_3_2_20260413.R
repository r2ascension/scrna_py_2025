#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison Pipeline v1.3.2-EPI Wrapper
# ==============================================================================
#
# Generic-interface entrypoint:
#   - uses tissue_comparison_generic_wrapper_20260412.R
#   - keeps the epithelial-specific preset thin and declarative
#   - reuses helper-synchronized marker panels, ssGSEA / CHOIR / discovery logic
#
# Date: 2026-04-13
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "EPITHELIAL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.3.2-EPI",
    PIPELINE_SUBTITLE = paste(
      "Epithelial Tissue Comparison v1.3.2-EPI",
      "(generic wrapper interface)"
    ),
    GENERATED_BY_LABEL = "epithelial_tissue_comparison_v1_3_2_20260413.R",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL TISSUE COMPARISON COMPLETE (v1.3.2-EPI)",
    OUTPUT_DIR = "/home/h2048/data/R/0413/epithelial_tissue_comparison_v1_3_2_20260413",
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] Epithelial wrapper is a thin preset over the generic wrapper instead of a standalone fork.",
      "  [EPI-1] Shared helper now owns epithelial marker panels plus common ssGSEA / discovery screening helpers.",
      "  [EPI-2] Shared engine now supports pseudobulk sample+tissue disambiguation for epithelial-safe aggregation.",
      "  [EPI-3] New run writes into an isolated 0413 output directory while using 0407 epithelial outputs as the historical baseline."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
