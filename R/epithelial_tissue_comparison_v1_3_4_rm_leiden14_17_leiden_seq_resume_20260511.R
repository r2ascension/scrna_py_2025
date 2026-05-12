#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison v1.3.4 — sequential Leiden recovery for 20260508 rerun
# ==============================================================================
# Date: 2026-05-11
# Purpose:
#   Resume the unfinished 2026-05-08 epithelial full rerun, but use the Leiden
#   backend instead of CHOIR. The input h5ad already contains
#   `leiden_Epithelial_res0.8`, which is reused directly as the cluster backend.
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)
source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R")
source("/home/h2048/script/R/tissue_comparison_recovery_lock_20260512.R")

output_dir <- "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508"
h5ad_path <- "/home/h2048/data/py/0508/epithelial_scanvi_rm_leiden14_17_20260508/epithelial_scanvi_rm_leiden14_17_SELF_for_R.h5ad"
previous_dir <- "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun"

cleanup_cluster_artifacts <- function(output_dir, backend = "leiden") {
  backend <- tolower(backend)
  report_dir <- file.path(output_dir, "reports")
  figure_dir <- file.path(output_dir, "figures")
  unlink(c(
    file.path(report_dir, backend),
    file.path(report_dir, sprintf("%s_llm_structured.tsv", backend)),
    file.path(report_dir, sprintf("%s_llm_structured_all.rds", backend)),
    file.path(report_dir, sprintf("llm_%s_discovery_screen.tsv", backend)),
    file.path(output_dir, sprintf("LLM_%s_INTERPRETATION.md", toupper(backend))),
    file.path(output_dir, sprintf("LLM_%s_DISCOVERY_REVIEW.md", toupper(backend)))
  ), recursive = TRUE, force = TRUE)
  unlink(file.path(report_dir, "choir"), recursive = TRUE, force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, sprintf("%s_umap*", backend))), force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, sprintf("ssgsea_%s_heatmap*", backend))), force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, "choir_umap*")), force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, "ssgsea_choir_heatmap*")), force = TRUE)
  unlink(c(
    file.path(output_dir, "epithelial_tissue_comparison_final.rds"),
    file.path(output_dir, "epithelial_tissue_comparison_final.h5ad"),
    file.path(output_dir, "REPORT.md")
  ), force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.recovery_lock_release <- tc_recovery_register_lock(
  output_dir = output_dir,
  lock_name = ".epithelial_recovery_20260511.lock",
  completion_paths = c(
    file.path(output_dir, "epithelial_tissue_comparison_final.rds"),
    file.path(output_dir, "epithelial_tissue_comparison_final.h5ad"),
    file.path(output_dir, "REPORT.md"),
    file.path(output_dir, "reports", "leiden_llm_structured.tsv")
  ),
  label = "epithelial Leiden recovery"
)
on.exit(.recovery_lock_release(), add = TRUE)
cleanup_cluster_artifacts(output_dir, backend = "leiden")

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "EPITHELIAL",
  overrides = list(
    H5AD_PATH = h5ad_path,
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = previous_dir,
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    REQUIRE_LLM = TRUE,
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    PIPELINE_VERSION_LABEL = "v1.3.4-EPI-rmLeiden14-17-fullscvi-leiden-seq-resume",
    PIPELINE_SUBTITLE = paste(
      "Epithelial Tissue Comparison v1.3.4",
      "(20260511 sequential recovery; Leiden backend from existing leiden_Epithelial_res0.8)"
    ),
    GENERATED_BY_LABEL = "epithelial_tissue_comparison_v1_3_4_rm_leiden14_17_leiden_seq_resume_20260511.R",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL TISSUE COMPARISON COMPLETE (v1.3.4 rmLeiden14/17 Leiden recovery)",
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    ANALYSIS_L3_DESCRIPTION = "scANVI-predicted epithelial L3 labels",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = TRUE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CLUSTERING = TRUE,
    RUN_CHOIR = FALSE,
    RUN_OFA = TRUE,
    CLUSTER_BACKEND = "LEIDEN",
    CLUSTER_EXISTING_COL = "leiden_Epithelial_res0.8",
    LEIDEN_RESOLUTION = 0.8,
    LEIDEN_N_DIMS = 30L,
    LEIDEN_K_PARAM = 30L,
    LEIDEN_ALGORITHM = 4L,
    RUN_CLUSTER_SSGSEA = TRUE,
    RUN_L3_OFA = TRUE,
    RUN_L3_OFA_VS_REST = TRUE,
    RUN_L3_OFA_SAME_L2 = FALSE,
    RUN_L3_OFA_INTER_TISSUE = TRUE,
    N_CORES = 2L,
    CHOIR_N_CORES = 1L,
    CHOIR_FIND_VAR_FEATURES_IF_MISSING = FALSE,
    CHOIR_VAR_FEATURES_MAX = 10L,
    SSGSEA_CHOIR_METHODS = c("hallmark"),
    OFA_MAX_CELLS_PER_IDENT = 3000L,
    L3_OFA_MAX_CELLS_PER_IDENT = 3000L,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "3")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [RECOVERY-1] Reuse the completed 20260508 epithelial scVI/scANVI h5ad and rerun unfinished downstream stages sequentially.",
      "  [RECOVERY-2] Use Leiden, not CHOIR, for epithelial cluster-level OFA/LLM as requested.",
      "  [RECOVERY-3] Reuse the existing h5ad column `leiden_Epithelial_res0.8` directly; do not recompute the 272k-cell clustering graph.",
      "  [LLM-1] Keep the 20260508 batch-parallel LLM screening overlay.",
      "  [PLOT-1] Pathway enrichment barplots are regenerated after successful completion."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
