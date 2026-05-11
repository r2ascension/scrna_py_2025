#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Tissue Comparison v2.6.4 — rm CHOIR/OFA outliers full rerun
# ==============================================================================
# Date: 2026-05-08
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)

base_preset <- tc_get_lineage_preset("TNK")

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
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "3")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
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

output_dir <- "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508"
if (dir.exists(output_dir)) {
  unlink(list.files(output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "TNK",
  overrides = list(
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    REQUIRE_LLM = TRUE,
    PIPELINE_VERSION_LABEL = "v2.6.4-TNK-rmCHOIR23-28-31-41-OFA41-66",
    PIPELINE_SUBTITLE = paste(
      "T/NK Tissue Comparison v2.6.4",
      "(full rerun after removing CHOIR c23/c28/c31/c41 and OFA c41/c66 before scVI/scANVI)"
    ),
    GENERATED_BY_LABEL = "tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508.R",
    LINEAGE_COMPLETION_BANNER = "T/NK TISSUE COMPARISON COMPLETE (v2.6.4 rmCHOIR/OFA full rerun)",
    H5AD_PATH = paste0(
      "/home/h2048/data/py/0508/tnk_scvi_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508/",
      "adata_tnk_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508.h5ad"
    ),
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper",
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_2_20260414.R",
    L3_TO_L2_REMAP = tnk_l3_to_l2_remap,
    CUSTOM_MARKERS_DB = tnk_custom_markers_db,
    SHARED_OVERRIDES = tnk_shared_overrides,
    PIPELINE_CHANGELOG_LINES = c(
      base_preset$PIPELINE_CHANGELOG_LINES,
      "  [TNK-8] Remove user-requested CHOIR clusters c23/c28/c31/c41 and OFA clusters c41/c66 using the current CHOIR backend cluster assignment.",
      "  [TNK-9] Retrain scVI and scANVI from the corrected 0414 refined-label TNK h5ad before rerunning the complete tissue-comparison workflow.",
      "  [LLM-1] Use 2026-05-08 batch-parallel LLM screening overlay.",
      "  [PLOT-1] Pathway enrichment barplots are generated post hoc from ssGSEA/OFA CSV reports."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
