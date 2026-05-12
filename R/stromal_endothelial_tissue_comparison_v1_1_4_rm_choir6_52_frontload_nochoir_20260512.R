#!/usr/bin/env Rscript
# ==============================================================================
# Stromal Endothelial Tissue Comparison v1.1.4 — front-load no-CHOIR recovery
# ==============================================================================
# Date: 2026-05-12
# Purpose:
#   While another lineage is still in a long CHOIR step, pre-load the 20260508
#   endothelial h5ad, metadata harmonization, cached pairwise/ssGSEA artifacts,
#   and report/final-object serialization without starting a second CHOIR job.
#   The full CHOIR wrapper can later reuse the same output directory.
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")
source(GENERIC_WRAPPER_PATH)
source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R")
source("/home/h2048/script/R/tissue_comparison_recovery_lock_20260512.R")

output_dir <- "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508"
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
  label = "stromal endothelial front-load no-CHOIR recovery"
)
on.exit(.recovery_lock_release(), add = TRUE)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "STROMAL_ENDOTHELIAL",
  overrides = list(
    PIPELINE_VERSION_LABEL = "v1.1.4-Endothelial-rmCHOIR6-52-frontload-noCHOIR",
    PIPELINE_SUBTITLE = paste(
      "Stromal Endothelial Tissue Comparison v1.1.4",
      "(20260512 front-load run: no CHOIR; cached upstream stages only)"
    ),
    GENERATED_BY_LABEL = "stromal_endothelial_tissue_comparison_v1_1_4_rm_choir6_52_frontload_nochoir_20260512.R",
    LINEAGE_COMPLETION_BANNER = "STROMAL ENDOTHELIAL FRONT-LOAD COMPLETE (v1.1.4 no CHOIR)",
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
    RUN_CLUSTERING = FALSE,
    RUN_CHOIR = FALSE,
    CLUSTER_BACKEND = "CHOIR",
    RUN_CLUSTER_SSGSEA = FALSE,
    RUN_OFA = FALSE,
    RUN_L3_OFA = FALSE,
    RUN_L3_OFA_VS_REST = FALSE,
    RUN_L3_OFA_SAME_L2 = FALSE,
    RUN_L3_OFA_INTER_TISSUE = FALSE,
    RUN_MILOPY = FALSE,
    N_CORES = 2L,
    CHOIR_N_CORES = 1L,
    ADVANCED_HELPER_PATH = "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "2")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),
    PIPELINE_CHANGELOG_LINES = c(
      "  [FRONTLOAD-1] Reuse the completed 20260508 endothelial branchwise scVI/scANVI h5ad while B cell CHOIR is still running.",
      "  [FRONTLOAD-2] Load cached pairwise DE, enrichment, LLM, and grouped ssGSEA artifacts from the current 0508 directory.",
      "  [FRONTLOAD-3] Explicitly skip CHOIR/OFA/L3-OFA so no second long CHOIR job is started in parallel.",
      "  [FRONTLOAD-4] Serialize a no-CHOIR front-load final object/report; the v1.1.3 CHOIR wrapper will overwrite final outputs after full completion."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
