#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Tissue Comparison Pipeline v2.6.3-TNK-relabel
# ==============================================================================
#
# Purpose:
#   - consume the 2026-04-14 corrected-label TNK scanvi rerun
#   - upgrade TNK wrapper execution to the 2026-04-14 generic helper path
#   - use the 2026-04-14 shared engine (v2.6.2) so LLM evidence includes
#     DEG expression proportions when available
#
# Date: 2026-04-14
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
source(GENERIC_WRAPPER_PATH)

base_preset <- tc_get_lineage_preset("TNK")

# Synchronize TNK label vocabulary with corrected scanvi labels.
tnk_custom_markers_db <- base_preset$CUSTOM_MARKERS_DB
tnk_custom_markers_db$subtype[tnk_custom_markers_db$subtype == "CD4 Naive/TCM"] <- "CD4 Naive"
tnk_custom_markers_db$subtype[tnk_custom_markers_db$subtype == "CD4 Tfr"] <- "CD4 Tfh"

tnk_l3_to_l2_remap <- c(
  "CD4 Naive" = "CD4 T cells",
  "CD4 Tcm" = "CD4 T cells",
  "CD4 Tfh" = "CD4 T cells",
  "CD4 Th1" = "CD4 T cells",
  "CD4 Th17" = "CD4 T cells",
  "CD4 Treg" = "CD4 T cells",
  "CD4 Trm" = "CD4 T cells",
  "CD8 Naive" = "CD8 T cells",
  "CD8 Teff" = "CD8 T cells",
  "CD8 Tem" = "CD8 T cells",
  "CD8 Temra" = "CD8 T cells",
  "CD8 Trm" = "CD8 T cells",
  "gdT" = "CD8 T cells",
  "MAIT" = "CD8 T cells",
  "ILC3" = "NK cells",
  "NK" = "NK cells",
  "NK Exhausted" = "NK cells"
)

tnk_shared_overrides <- utils::modifyList(
  base_preset$SHARED_OVERRIDES,
  list(
    LLM_INCLUDE_TOP_DEG = TRUE,
    LLM_TOP_DEG_N = 10L,
    LLM_DEG_PADJ_THR = 0.10,
    LLM_DEG_LFC_THR = 0.15,
    LLM_REQUIRE_INTEGRATED_UP_DOWN = TRUE,
    LLM_EXTRA_RULES = c(
      base_preset$SHARED_OVERRIDES$LLM_EXTRA_RULES,
      "When DEG expression proportion evidence is present, explicitly use it to judge whether the claimed TNK subtype is supported by broad within-group expression or only by sparse marker leakage."
    )
  )
)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "TNK",
  overrides = list(
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    PIPELINE_VERSION_LABEL = "v2.6.3-TNK-relabel",
    PIPELINE_SUBTITLE = paste(
      "T/NK Tissue Comparison v2.6.3-TNK",
      "(corrected scanvi labels + 2026-04-14 helper/engine)"
    ),
    GENERATED_BY_LABEL = "tnk_tissue_comparison_v2_6_3_20260414_relabel_helper.R",
    LINEAGE_COMPLETION_BANNER = "T/NK TISSUE COMPARISON COMPLETE (v2.6.3-TNK-relabel)",
    H5AD_PATH = paste0(
      "/home/h2048/data/py/0414/tnk_scanvi_relabel_rerun_20260414/",
      "adata_tnk_scanvi_refined_relabel_rerun_20260414.h5ad"
    ),
    OUTPUT_DIR = "/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_2_20260413_rm_choir",
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_2_20260414.R",
    L3_TO_L2_REMAP = tnk_l3_to_l2_remap,
    CUSTOM_MARKERS_DB = tnk_custom_markers_db,
    SHARED_OVERRIDES = tnk_shared_overrides,
    PIPELINE_CHANGELOG_LINES = c(
      base_preset$PIPELINE_CHANGELOG_LINES,
      "  [TNK-5] Reuses the 0413 filtered scVI state but reruns scANVI only after relabeling `CD4 Naive/TCM` -> `CD4 Naive` and `CD4 Tfr` -> `CD4 Tfh`.",
      "  [TNK-6] TNK wrapper now runs through the 2026-04-14 generic helper path and the shared engine v2.6.2 (2026-04-14) so LLM evidence can incorporate DEG expression proportions when available.",
      "  [TNK-7] TNK-specific L3 remap table and custom marker DB were synchronized to the corrected label vocabulary to prevent old names from leaking into reports."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
