#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison v1.3.3 — rm Leiden/OFA 14/17 full rerun
# ==============================================================================
# Date: 2026-05-08
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)

output_dir <- "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508"
h5ad_path <- "/home/h2048/data/py/0508/epithelial_scanvi_rm_leiden14_17_20260508/epithelial_scanvi_rm_leiden14_17_SELF_for_R.h5ad"
previous_dir <- "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun"

if (dir.exists(output_dir)) {
  unlink(list.files(output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "EPITHELIAL",
  overrides = list(
    H5AD_PATH = h5ad_path,
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = previous_dir,
    REQUIRE_LLM = TRUE,
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS = FALSE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = FALSE,
    PIPELINE_VERSION_LABEL = "v1.3.3-EPI-rmLeiden14-17-fullscvi",
    PIPELINE_SUBTITLE = paste(
      "Epithelial Tissue Comparison v1.3.3",
      "(full rerun after removing OFA/Leiden c14 and c17 before scVI/scANVI)"
    ),
    GENERATED_BY_LABEL = "epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508.R",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL TISSUE COMPARISON COMPLETE (v1.3.3 rmLeiden14/17 full rerun)",
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    ANALYSIS_L3_DESCRIPTION = "scANVI-predicted epithelial L3 labels",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    RUN_VISUALIZATION = TRUE,
    RUN_PAIRWISE_DE = TRUE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = TRUE,
    RUN_CLUSTERING = TRUE,
    RUN_CHOIR = FALSE,
    RUN_OFA = TRUE,
    CLUSTER_BACKEND = "LEIDEN",
    LEIDEN_RESOLUTION = 0.8,
    LEIDEN_N_DIMS = 30L,
    LEIDEN_K_PARAM = 30L,
    LEIDEN_ALGORITHM = 4L,
    RUN_CLUSTER_SSGSEA = TRUE,
    RUN_L3_OFA = TRUE,
    RUN_L3_OFA_VS_REST = TRUE,
    RUN_L3_OFA_INTER_TISSUE = TRUE,
    N_CORES = 2L,
    CHOIR_N_CORES = 1L,
    SSGSEA_CHOIR_METHODS = c("hallmark"),
    OFA_MAX_CELLS_PER_IDENT = 3000L,
    L3_OFA_MAX_CELLS_PER_IDENT = 3000L,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "3")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [DROP-1] Remove the user-curated epithelial OFA/Leiden outlier clusters 14 and 17 at h5ad level before training.",
      "  [FULL-1] Regenerate grouped ssGSEA, pairwise DE, Leiden cluster, OFA, L3 OFA, reports, and LLM outputs from fresh scANVI annotations.",
      "  [LLM-1] Use 2026-05-08 batch-parallel LLM screening overlay.",
      "  [PLOT-1] Pathway enrichment barplots are generated post hoc from ssGSEA/OFA CSV reports."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
