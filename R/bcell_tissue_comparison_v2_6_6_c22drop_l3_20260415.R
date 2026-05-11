#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.6-c22drop-l3 Wrapper
# ==============================================================================
#
# Purpose:
#   - rerun B-cell tissue comparison after switching c22-drop scANVI supervision
#     from L2 to L3 labels
#   - collapse all non-nose GC-B cells into Memory_B from the scANVI rerun onward
#   - force a fresh load from the updated L3-rerun reference h5ad
#   - keep the generic wrapper / shared engine architecture unchanged
#
# Date: 2026-04-15
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.6-c22drop-l3-gcb2mem",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.6-c22drop-l3-gcb2mem",
      "(fresh rerun after collapsing all non-nose GC-B labels into Memory_B from scANVI onward)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.6-c22drop-l3-gcb2mem)",
    H5AD_PATH = "/home/h2048/data/py/0415/bcell_scvi_scanvi_ref_c22drop_l3_20260415/bcell_reference_c22drop_scanvi_L3_ref_20260415.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0413/bcell_tissue_comparison_v2_6_5_c22drop_20260413",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    PIPELINE_CHANGELOG_LINES = c(
      "  [L3-1] Remove CHOIR cluster 22 contamination as before, but rerun scANVI with L3 supervision (`cell_type_expert`) instead of L2.",
      "  [L3-2] Fresh scANVI L3 predictions are written back to `cell_type_scanvi_pred` for downstream tissue comparison and LLM review.",
      "  [L3-3] L2 labels are now derived from rerun L3 predictions, not reused from the previous reference metadata.",
      "  [L3-4] All non-nose GC-B labels are collapsed into `Memory_B` before scANVI training; final non-nose GC-B predictions are guard-railed to `Memory_B` as well.",
      "  [L3-5] This rerun intentionally overwrites the previous 0415 outputs in place after pre-cleaning stale artifacts.",
      "  [ARCH-1] Generic wrapper interface still owns config validation, env loading, and previous-run preflight.",
      "  [FINAL-1] Compact integrated LLM outputs only; no split final conclusions by up/down or per database.",
      "  [FINAL-2] Non-informative pathways (energy / mitochondrial / ribosomal / ENSG-like) do not occupy top-pathway slots.",
      "  [FINAL-3] Non-epithelial cilia-related genes/pathways are suppressed from enrichment and LLM emphasis.",
      "  [FINAL-4] Previous output reading remains compatible with legacy table/RDS layouts."
    )
  )
)

if (dir.exists(PIPELINE_CONFIG$OUTPUT_DIR)) {
  stale_paths <- list.files(
    PIPELINE_CONFIG$OUTPUT_DIR,
    all.files = TRUE,
    no.. = TRUE,
    full.names = TRUE
  )
  if (length(stale_paths) > 0) {
    unlink(stale_paths, recursive = TRUE, force = TRUE)
  }
}
dir.create(PIPELINE_CONFIG$OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
