#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Tissue Comparison Pipeline v1.2.1-MYELOID Wrapper
# ==============================================================================
#
# Generic-interface entrypoint:
#   - uses tissue_comparison_generic_wrapper_20260414_v2.R
#   - keeps the myeloid branch thin and declarative
#   - forces a fresh full rerun with wrapper-managed LLM requirement
#
# Smoke test:
#   PIPELINE_TEST_MODE=true Rscript myeloid_tissue_comparison_v1_2_1_20260414.R
#
# Date: 2026-04-14
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414_v2.R"
source(GENERIC_WRAPPER_PATH)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "MYELOID"
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
