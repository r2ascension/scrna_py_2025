#!/usr/bin/env Rscript
# ==============================================================================
# Stromal Endothelial Tissue Comparison v1.1.3 — sequential CHOIR recovery
# ==============================================================================
# Date: 2026-05-11
# Purpose:
#   Resume the unfinished 2026-05-08 endothelial full rerun in-place, rerunning
#   the unfinished CHOIR/OFA/LLM/finalization stages sequentially. Python
#   branchwise scVI/scANVI h5ad is reused.
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)
source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R")
source("/home/h2048/script/R/tissue_comparison_recovery_lock_20260512.R")

output_dir <- "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508"

cleanup_cluster_artifacts <- function(output_dir, backend = "choir") {
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
  unlink(Sys.glob(file.path(figure_dir, sprintf("%s_umap*", backend))), force = TRUE)
  unlink(Sys.glob(file.path(figure_dir, sprintf("ssgsea_%s_heatmap*", backend))), force = TRUE)
  unlink(c(
    file.path(output_dir, "stromal_endothelial_tissue_comparison_final.rds"),
    file.path(output_dir, "stromal_endothelial_tissue_comparison_final.h5ad"),
    file.path(output_dir, "REPORT.md")
  ), force = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.recovery_lock_release <- tc_recovery_register_lock(
  output_dir = output_dir,
  lock_name = ".stromal_endothelial_recovery_20260511.lock",
  completion_paths = c(
    file.path(output_dir, "stromal_endothelial_tissue_comparison_final.rds"),
    file.path(output_dir, "stromal_endothelial_tissue_comparison_final.h5ad"),
    file.path(output_dir, "REPORT.md"),
    file.path(output_dir, "reports", "choir_llm_structured.tsv")
  ),
  label = "stromal endothelial CHOIR recovery"
)
on.exit(.recovery_lock_release(), add = TRUE)
cleanup_cluster_artifacts(output_dir, backend = "choir")

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_ENDOTHELIAL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.1.3-Endothelial-rmCHOIR6-52-fullscvi-choir-seq-resume",
    PIPELINE_SUBTITLE = paste(
      "Stromal Endothelial Tissue Comparison v1.1.3",
      "(20260511 sequential CHOIR recovery from the 20260508 branchwise scVI/scANVI rerun)"
    ),
    GENERATED_BY_LABEL = "stromal_endothelial_tissue_comparison_v1_1_3_rm_choir6_52_choir_seq_resume_20260511.R",
    LINEAGE_COMPLETION_BANNER = "STROMAL ENDOTHELIAL TISSUE COMPARISON COMPLETE (v1.1.3 sequential CHOIR recovery)",
    H5AD_PATH = "/home/h2048/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508/endothelial/adata_endothelial_reference_v1_5_branchwise.h5ad",
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    REQUIRE_LLM = TRUE,
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CLUSTERING = TRUE,
    RUN_CHOIR = TRUE,
    CLUSTER_BACKEND = "CHOIR",
    RUN_CLUSTER_SSGSEA = TRUE,
    RUN_OFA = TRUE,
    RUN_L3_OFA = TRUE,
    RUN_L3_OFA_VS_REST = TRUE,
    RUN_L3_OFA_SAME_L2 = FALSE,
    RUN_L3_OFA_INTER_TISSUE = TRUE,
    N_CORES = 4L,
    CHOIR_N_CORES = 4L,
    CHOIR_MAX_REPEAT_ERRORS = 5L,
    CHOIR_SUBTREE_REDUCTIONS = TRUE,
    OFA_MAX_CELLS_PER_IDENT = 3000L,
    L3_OFA_MAX_CELLS_PER_IDENT = 3000L,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "3")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [RECOVERY-1] Reuse the completed 20260508 endothelial branchwise scVI/scANVI h5ad and rerun unfinished downstream stages sequentially.",
      "  [RECOVERY-2] Use the 20260410 shared engine with explicit CHOIR/Leiden backend support and immediate stage checkpoints.",
      "  [RECOVERY-3] Load cached pairwise DE and grouped ssGSEA artifacts from the current 0508 directory, then re-run CHOIR/OFA/LLM/finalization.",
      "  [LLM-1] Keep the 20260508 batch-parallel LLM screening overlay.",
      "  [PLOT-1] Pathway enrichment barplots are regenerated after successful completion."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
