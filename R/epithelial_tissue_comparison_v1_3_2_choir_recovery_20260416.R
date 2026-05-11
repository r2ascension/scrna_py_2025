#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison v1.3.2-EPI - CHOIR Recovery Runner
# ==============================================================================
#
# Purpose:
#   Recover the interrupted 2026-04-15 epithelial full rerun by:
#     - reusing the validated 0413 final Seurat object as execution base
#     - reusing existing/cached 0415 grouped-ssGSEA artifacts when present
#     - falling back to 0413 report RDS artifacts for pairwise DE / enrichment / LLM
#     - skipping the expensive already-finished front half of the pipeline
#     - re-running CHOIR + cluster ssGSEA + OFA + CHOIR LLM + final REPORT/export
#
# Notes:
#   - N_CORES is intentionally reduced to 2 to lower peak memory pressure during
#     CHOIR subtree processing after the previous 4-core run exited mid-CHOIR.
#   - OUTPUT_DIR remains the original 0415 rerun directory so missing artifacts
#     are backfilled in place.
#
# Date: 2026-04-16
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

previous_dir <- "/home/h2048/data/R/0413/epithelial_tissue_comparison_v1_3_2_20260413"
output_dir <- "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun"

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "EPITHELIAL",
  overrides = list(
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = previous_dir,
    PREVIOUS_FINAL_OBJECT_RDS = file.path(previous_dir, "epithelial_tissue_comparison_final.rds"),
    REUSE_PREVIOUS_FINAL_OBJECT = TRUE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    REQUIRE_LLM = TRUE,
    PIPELINE_VERSION_LABEL = "v1.3.2-EPI-CHOIR-RECOVERY",
    PIPELINE_SUBTITLE = paste(
      "Epithelial Tissue Comparison v1.3.2-EPI",
      "(2026-04-16 CHOIR recovery from interrupted 0415 rerun)"
    ),
    GENERATED_BY_LABEL = "epithelial_tissue_comparison_v1_3_2_choir_recovery_20260416.R",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL TISSUE COMPARISON COMPLETE (v1.3.2-EPI CHOIR recovery 2026-04-16)",
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CHOIR = TRUE,
    RUN_OFA = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    N_CORES = 2L,
    PIPELINE_CHANGELOG_LINES = c(
      "  [RECOVERY-1] 2026-04-16 recovery run resumes the interrupted 0415 rerun in-place.",
      "  [RECOVERY-2] Pairwise DE / enrichment / grouped ssGSEA / wilcox are skipped and backfilled from existing report RDS artifacts.",
      "  [RECOVERY-3] Current 0415 grouped-ssGSEA artifacts are preferred; missing front-half report RDS objects fall back to 0413.",
      "  [RECOVERY-4] Visualization is skipped to reuse existing 0415 figure outputs and shorten restart latency.",
      "  [RECOVERY-5] CHOIR recovery runs with N_CORES=2 to reduce peak memory and avoid another mid-CHOIR exit."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
