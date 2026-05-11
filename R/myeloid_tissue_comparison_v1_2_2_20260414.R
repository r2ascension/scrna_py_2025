#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Tissue Comparison v1.2.2 (Generic Wrapper, Tuned CHOIR)
# ==============================================================================
#
# Purpose:
#   Thin myeloid entrypoint over the generic wrapper, with versioned CHOIR tuning
#   restored in the shared engine to avoid the pathological heavy-default branch
#   seen in v1.2.1.
#
# Smoke test:
#   PIPELINE_TEST_MODE=true Rscript myeloid_tissue_comparison_v1_2_2_20260414.R
#
# Date: 2026-04-14
# ==============================================================================

source("/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414_v3.R")

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(lineage = "MYELOID")
tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
