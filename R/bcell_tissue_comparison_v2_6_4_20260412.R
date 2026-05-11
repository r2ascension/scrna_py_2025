#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.4 Wrapper
# ==============================================================================
#
# Generic-interface entrypoint:
#   - uses tissue_comparison_generic_wrapper_20260412.R
#   - keeps the B-cell-specific preset thin and declarative
#   - preserves explicit preflight / env / output-dir separation
#
# Date: 2026-04-12
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.4",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.4",
      "(generic wrapper interface)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_4_20260412.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.4)",
    OUTPUT_DIR = "/home/h2048/data/R/0412/bcell_tissue_comparison_v2_6_4_20260412",
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] B-cell wrapper is a thin preset over the generic wrapper instead of a one-off script.",
      "  [ARCH-3] OUTPUT_DIR and PREVIOUS_OUTPUT_DIR remain separated by default; in-place resume requires explicit opt-in.",
      "  [FINAL-1] Compact integrated LLM outputs only; no split final conclusions by up/down or per database.",
      "  [FINAL-2] Non-informative pathways (energy / mitochondrial / ribosomal / ENSG-like) do not occupy top-pathway slots.",
      "  [FINAL-3] Non-epithelial cilia-related genes/pathways are suppressed from enrichment and LLM emphasis.",
      "  [FINAL-4] Previous output reading remains compatible with legacy table/RDS layouts."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
