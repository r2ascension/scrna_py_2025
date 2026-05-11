#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.7-c13c25drop Wrapper
# ==============================================================================
#
# Purpose:
#   - preserve the previous c22 removal, additionally remove CHOIR clusters 13 and 25
#   - force a completely fresh rerun from scVI to scANVI to the final LLM outputs
#   - force a fresh load from the updated L3-rerun reference h5ad
#   - overwrite the previous 0415 R outputs in place after cleaning stale artifacts
#   - keep the generic wrapper / shared engine architecture unchanged
#
# Date: 2026-04-17
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v2.6.7-c22-c13-c25drop-l3-gcb2mem-freshscvi",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.7-c22-c13-c25drop-l3-gcb2mem-freshscvi",
      "(full rerun from original reference after removing CHOIR c22 + c13 + c25 and retraining fresh scVI/scANVI)"
    ),
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_7_c13c25drop_20260417.R",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.7-c22-c13-c25drop-l3-gcb2mem-freshscvi)",
    H5AD_PATH = "/home/h2048/data/py/0415/bcell_scvi_scanvi_ref_c22drop_l3_20260415/bcell_reference_c22_c13_c25drop_scanvi_L3_ref_20260417.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0413/bcell_tissue_comparison_v2_6_5_c22drop_20260413",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    PIPELINE_CHANGELOG_LINES = c(
      "  [DROP-1] Preserve the validated c22 removal and additionally remove CHOIR clusters 13 and 25 from the current 0415 B-cell run.",
      "  [DROP-2] Start again from the original full-gene reference h5ad so HVG selection, scVI, scANVI, UMAP, and downstream reports are all regenerated fresh.",
      "  [DROP-3] Fresh scANVI L3 predictions are written back to `cell_type_scanvi_pred` for downstream tissue comparison and LLM review.",
      "  [DROP-4] L2 labels are derived from fresh rerun L3 predictions, not reused from previous reference metadata.",
      "  [DROP-5] All non-nose GC-B labels are collapsed into `Memory_B` before scANVI training; final non-nose GC-B predictions are guard-railed to `Memory_B` as well.",
      "  [DROP-6] This rerun intentionally overwrites the previous 0415 outputs in place after pre-cleaning stale artifacts.",
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
