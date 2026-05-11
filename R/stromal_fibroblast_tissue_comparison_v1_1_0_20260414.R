#!/usr/bin/env Rscript
# ==============================================================================
# Stromal Fibroblast Tissue Comparison Pipeline v1.1.0
# ==============================================================================
#
# Generic-interface entrypoint:
#   - uses tissue_comparison_generic_wrapper_20260414.R
#   - keeps the stromal fibroblast branch thin and declarative
#   - forces a fresh full rerun with wrapper-managed LLM requirement
#
# Smoke test:
#   PIPELINE_TEST_MODE=true Rscript stromal_fibroblast_tissue_comparison_v1_1_0_20260414.R
#
# Date: 2026-04-14
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_FIBROBLAST"
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)