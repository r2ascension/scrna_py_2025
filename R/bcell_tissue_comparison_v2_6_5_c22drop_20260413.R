#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.5-c22drop Wrapper
# ==============================================================================
#
# Purpose:
#   - rerun B-cell tissue comparison after removing CHOIR cluster 22 contamination
#   - force a fresh load from the updated reference h5ad (no previous final-object reuse)
#   - keep the generic wrapper / shared engine architecture unchanged
#
# Date: 2026-04-13
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.5-c22drop",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.5-c22drop",
      "(fresh rerun after removing CHOIR cluster 22)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_5_c22drop_20260413.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.5-c22drop)",
    H5AD_PATH = "/home/h2048/data/py/0413/bcell_scvi_scanvi_ref_c22drop_20260413/bcell_reference_c22drop_scanvi_ref_20260413.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0413/bcell_tissue_comparison_v2_6_5_c22drop_20260413",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0412/bcell_tissue_comparison_v2_6_4_20260412",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    PIPELINE_CHANGELOG_LINES = c(
      "  [C22-1] Remove CHOIR cluster 22 fibroblast-like contamination before rerunning reference scVI/scANVI.",
      "  [C22-2] Force fresh h5ad load from the updated c22-drop reference; do not reuse the previous final Seurat object.",
      "  [C22-3] Updated reference h5ad preserves raw.X plus layers['counts'] for safer downstream reuse.",
      "  [ARCH-1] Generic wrapper interface still owns config validation, env loading, and previous-run preflight.",
      "  [FINAL-1] Compact integrated LLM outputs only; no split final conclusions by up/down or per database.",
      "  [FINAL-2] Non-informative pathways (energy / mitochondrial / ribosomal / ENSG-like) do not occupy top-pathway slots.",
      "  [FINAL-3] Non-epithelial cilia-related genes/pathways are suppressed from enrichment and LLM emphasis.",
      "  [FINAL-4] Previous output reading remains compatible with legacy table/RDS layouts."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
