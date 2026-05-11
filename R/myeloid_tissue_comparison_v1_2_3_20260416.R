#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Tissue Comparison Pipeline v1.2.3 (2026-04-16 rerun)
# ==============================================================================
#
# Purpose:
#   - rerun the myeloid tissue-comparison pipeline on top of a tissue-aware
#     patched `cell_type_L3_refined` h5ad
#   - keep previous validated 2026-04-14 scripts immutable
#   - use a 2026-04-16 helper overlay that makes ssGSEA-only outlier calls more
#     conservative
#
# Date: 2026-04-16
# ==============================================================================

source("/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414_v3.R")

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "MYELOID",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.2.3-MYELOID",
    PIPELINE_SUBTITLE = paste(
      "Myeloid Tissue Comparison v1.2.3-MYELOID",
      "(tissue-aware alveolar relabel + conservative ssGSEA review)"
    ),
    GENERATED_BY_LABEL = "myeloid_tissue_comparison_v1_2_3_20260416.R",
    LINEAGE_COMPLETION_BANNER = "MYELOID TISSUE COMPARISON COMPLETE (v1.2.3-MYELOID)",
    OUTPUT_DIR = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414",
    H5AD_PATH = "/home/h2048/data/py/0416/adata_myeloid_L3refined_tissueaware_patched_v1.h5ad",
    REQUIRE_LLM = FALSE
  )
)

PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS <- file.path(
  PIPELINE_CONFIG$PREVIOUS_OUTPUT_DIR,
  paste0(PIPELINE_CONFIG$FINAL_FILE_PREFIX, ".rds")
)

PIPELINE_CONFIG$SHARED_OVERRIDES <- utils::modifyList(
  PIPELINE_CONFIG$SHARED_OVERRIDES,
  list(
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260416_myeloid_ssgsea_v1.R",
    ADVANCED_HELPER_ALREADY_LOADED = FALSE,
    LLM_EXTRA_RULES = c(
      PIPELINE_CONFIG$SHARED_OVERRIDES$LLM_EXTRA_RULES,
      "Outside lung parenchyma and respiratory airway, alveolar macrophage labels were tissue-aware remapped upstream to Interstitial macrophages before R analysis.",
      "For ssGSEA-only judgments, absence of canonical macrophage pathways is not by itself evidence of contamination; prefer mixed_or_uncertain unless explicit alternative-lineage or artifact evidence is present."
    )
  )
)

PIPELINE_CONFIG$PIPELINE_CHANGELOG_LINES <- c(
  PIPELINE_CONFIG$PIPELINE_CHANGELOG_LINES,
  "  [MYELOID-7] Input h5ad now comes from a 2026-04-16 tissue-aware patch that remaps non-respiratory alveolar macrophage labels to Interstitial macrophages.",
  "  [MYELOID-8] ssGSEA grouped review and ssGSEA discovery screening use a conservative 2026-04-16 helper overlay so non-canonical pathway profiles are downgraded to mixed/uncertain unless contamination is directly supported.",
  "  [MYELOID-9] Wrapper defaults to REQUIRE_LLM=FALSE so the 2026-04-16 rerun can complete in degraded mode when only placeholder DEEPSEEK keys are available; set PIPELINE_REQUIRE_LLM=true with a real key to restore full LLM outputs."
)

PIPELINE_CONFIG$L3_TO_L2_REMAP <- c(
  PIPELINE_CONFIG$L3_TO_L2_REMAP,
  "Non-classical monocytes_c0" = "Monocyte",
  "Non-classical monocytes_c1" = "Monocyte"
)

PIPELINE_CONFIG$PIPELINE_CHANGELOG_LINES <- c(
  PIPELINE_CONFIG$PIPELINE_CHANGELOG_LINES,
  "  [MYELOID-10] Accept retained hierarchical non-classical monocyte labels (`Non-classical monocytes_c0/c1`) and remap both to Monocyte during L2 validation."
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
