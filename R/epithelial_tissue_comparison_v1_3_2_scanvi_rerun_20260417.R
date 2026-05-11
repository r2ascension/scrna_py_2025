#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison v1.3.2-EPI — scANVI rerun from Python GPU output
# ==============================================================================
# Purpose:
#   - consume fresh epithelial scANVI predictions exported from scvi_env (GPU)
#   - resume the 0415 epithelial rerun from refreshed annotations without redoing finished grouped ssGSEA
#   - replace the old CHOIR cluster branch with standard Leiden clustering
#   - add L3-vs-rest and same-L3 cross-tissue OFA outputs
#   - regenerate downstream pairwise / cluster / report / LLM outputs on the refreshed annotations
#
# Date: 2026-04-17
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

previous_dir <- "/home/h2048/data/R/0413/epithelial_tissue_comparison_v1_3_2_20260413"
output_dir <- "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun"
h5ad_path <- "/home/h2048/data/py/0417/epithelial_scanvi_v2_8_gpu_rerun/epithelial_scanvi_v2_8_gpu_SELF_for_R.h5ad"

paths_to_clean <- c(
  file.path(output_dir, "figures"),
  file.path(output_dir, "pertpy_milo"),
  file.path(output_dir, "pseudobulk_de"),
  file.path(output_dir, "pseudobulk_de_L3"),
  file.path(output_dir, "wilcox_exploratory"),
  file.path(output_dir, "wilcox_exploratory_L3"),
  file.path(output_dir, "reports", "milopy"),
  file.path(output_dir, "reports", "choir"),
  file.path(output_dir, "reports", "leiden"),
  file.path(output_dir, "reports", "l3_ofa"),
  file.path(output_dir, "reports", "pseudobulk_de_all.rds"),
  file.path(output_dir, "reports", "pseudobulk_de_L3_all.rds"),
  file.path(output_dir, "reports", "enrichment_all.rds"),
  file.path(output_dir, "reports", "enrichment_L3_all.rds"),
  file.path(output_dir, "reports", "interpret_agent_all.rds"),
  file.path(output_dir, "reports", "interpret_agent_L3_all.rds"),
  file.path(output_dir, "reports", "interpret_agent_structured_all.rds"),
  file.path(output_dir, "reports", "interpret_agent_structured_L3_all.rds"),
  file.path(output_dir, "reports", "interpret_agent_structured.tsv"),
  file.path(output_dir, "reports", "interpret_agent_structured_L3.tsv"),
  file.path(output_dir, "reports", "choir_llm_structured.tsv"),
  file.path(output_dir, "reports", "choir_llm_structured_all.rds"),
  file.path(output_dir, "reports", "leiden_llm_structured.tsv"),
  file.path(output_dir, "reports", "leiden_llm_structured_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_vs_rest_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_same_l2_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_inter_tissue_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_vs_rest_summary.tsv"),
  file.path(output_dir, "reports", "l3_ofa_same_l2_summary.tsv"),
  file.path(output_dir, "reports", "l3_ofa_inter_tissue_summary.tsv"),
  file.path(output_dir, "reports", "milopy_run.rds"),
  file.path(output_dir, "reports", "milopy_level_summary.tsv"),
  file.path(output_dir, "reports", "milopy_pairwise_summary.tsv"),
  file.path(output_dir, "reports", "milopy_contrast_summary.tsv"),
  file.path(output_dir, "REPORT.md"),
  file.path(output_dir, "LLM_INTERPRETATION.md"),
  file.path(output_dir, "LLM_SSGSEA_INTERPRETATION.md"),
  file.path(output_dir, "LLM_CHOIR_INTERPRETATION.md"),
  file.path(output_dir, "LLM_LEIDEN_INTERPRETATION.md"),
  file.path(output_dir, "session_info.txt"),
  file.path(output_dir, "epithelial_tissue_comparison_final.rds"),
  file.path(output_dir, "epithelial_tissue_comparison_final.h5ad")
)

for (path in paths_to_clean) {
  if (file.exists(path) || dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
    cat(sprintf("[CLEAN] Removed stale path: %s\n", path))
  }
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
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    PIPELINE_VERSION_LABEL = "v1.3.2-EPI-SCANVI-LEIDEN-RESUME",
    PIPELINE_SUBTITLE = paste(
      "Epithelial Tissue Comparison v1.3.2-EPI",
      "(2026-04-17 Leiden resume from fresh GPU scANVI predictions)"
    ),
    GENERATED_BY_LABEL = "epithelial_tissue_comparison_v1_3_2_scanvi_rerun_20260417.R",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL TISSUE COMPARISON COMPLETE (v1.3.2-EPI scanvi Leiden resume 2026-04-17)",
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    ANALYSIS_L3_DESCRIPTION = "scANVI-predicted epithelial L3 labels",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    RUN_VISUALIZATION = TRUE,
    RUN_PAIRWISE_DE = TRUE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
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
    CHOIR_VAR_FEATURES_MAX = 2000L,
    CHOIR_N_ITERATIONS = 50L,
    CHOIR_N_TREES = 25L,
    CHOIR_SAMPLE_MAX = 8000L,
    CHOIR_MIN_CLUSTER_DEPTH = 3000L,
    SSGSEA_CHOIR_METHODS = c("hallmark"),
    OFA_MAX_CELLS_PER_IDENT = 3000L,
    L3_OFA_MAX_CELLS_PER_IDENT = 3000L,
    PIPELINE_CHANGELOG_LINES = c(
      "  [SCANVI-1] Refresh epithelial tissue comparison from fresh GPU scANVI predictions exported by scvi_env.",
      "  [SCANVI-2] Downstream R analysis consumes `cell_type_scanvi_pred` as the primary L3 annotation source.",
      "  [RESUME-1] Preserve completed grouped ssGSEA artifacts in the 0415 output directory and regenerate only downstream pairwise / cluster / report outputs.",
      "  [LEIDEN-1] Replace the previous CHOIR cluster branch with standard Seurat Leiden clustering on the transferred latent space.",
      "  [LEIDEN-2] Cluster-level ssGSEA remains limited to hallmark and one-vs-rest OFA downsampling is capped at 3000 cells per identity.",
      "  [L3-OFA-1] Add per-L3 vs rest OFA and same-L3 cross-tissue OFA summaries under reports/l3_ofa.",
      "  [WX-1] Exploratory Wilcoxon is skipped to prioritize pairwise, cached grouped ssGSEA, cluster, and L3 OFA outputs."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)

if (sys.nframe() == 0L) {
  quit(save = "no", status = 0, runLast = FALSE)
}
