#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial L3 OFA only runner
# ==============================================================================
# Purpose:
#   - load the refreshed epithelial scANVI h5ad
#   - skip heavyweight visualization / pairwise DE / ssGSEA / clustering branches
#   - run L3 OFA vs rest, same-L2 sibling OFA, and same-L3 cross-tissue OFA
#   - materialize the requested enrichment dotplots under reports/l3_ofa/
#
# Date: 2026-04-19
# ==============================================================================

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
source(GENERIC_WRAPPER_PATH)

args <- commandArgs(trailingOnly = TRUE)
append_mode <- FALSE
if ("--append" %in% args) {
  append_mode <- TRUE
  args <- setdiff(args, "--append")
}
if (length(args) > 0) {
  L3_TYPES_INCLUDE <- unique(trimws(args[nzchar(trimws(args))]))
  cat(sprintf("[INFO] Requested L3 subset: %s\n", paste(L3_TYPES_INCLUDE, collapse = ", ")))
}

previous_dir <- "/home/h2048/data/R/0413/epithelial_tissue_comparison_v1_3_2_20260413"
output_dir <- "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun"
h5ad_path <- "/home/h2048/data/py/0417/epithelial_scanvi_v2_8_gpu_rerun/epithelial_scanvi_v2_8_gpu_SELF_for_R.h5ad"

paths_to_clean <- c(
  file.path(output_dir, "reports", "l3_ofa"),
  file.path(output_dir, "reports", "l3_ofa_vs_rest_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_same_l2_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_inter_tissue_all.rds"),
  file.path(output_dir, "reports", "l3_ofa_vs_rest_summary.tsv"),
  file.path(output_dir, "reports", "l3_ofa_same_l2_summary.tsv"),
  file.path(output_dir, "reports", "l3_ofa_inter_tissue_summary.tsv")
)

if (!append_mode) {
  for (path in paths_to_clean) {
    if (file.exists(path) || dir.exists(path)) {
      unlink(path, recursive = TRUE, force = TRUE)
      cat(sprintf("[CLEAN] Removed stale path: %s\n", path))
    }
  }
} else {
  cat("[INFO] Append mode enabled; keeping existing L3 OFA outputs and adding new batch results\n")
}

dir.create(file.path(output_dir, "reports"), recursive = TRUE, showWarnings = FALSE)

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "EPITHELIAL",
  overrides = list(
    H5AD_PATH = h5ad_path,
    OUTPUT_DIR = output_dir,
    PREVIOUS_OUTPUT_DIR = previous_dir,
    REQUIRE_LLM = FALSE,
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = FALSE,
    LOAD_EXISTING_STAGE_ARTIFACTS = TRUE,
    LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT = TRUE,
    SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE = TRUE,
    PIPELINE_VERSION_LABEL = "v1.3.2-EPI-L3-OFA-ONLY",
    PIPELINE_SUBTITLE = paste(
      "Epithelial L3 OFA only",
      "(2026-04-19 slim rerun from fresh GPU scANVI predictions)"
    ),
    GENERATED_BY_LABEL = "epithelial_l3_ofa_only_20260419.R",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL L3 OFA ONLY COMPLETE (2026-04-19)",
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    ANALYSIS_L3_DESCRIPTION = "scANVI-predicted epithelial L3 labels",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CLUSTERING = FALSE,
    RUN_CHOIR = FALSE,
    RUN_OFA = FALSE,
    RUN_CLUSTER_SSGSEA = FALSE,
    RUN_L3_OFA = TRUE,
    RUN_L3_OFA_VS_REST = TRUE,
    RUN_L3_OFA_SAME_L2 = TRUE,
    RUN_L3_OFA_INTER_TISSUE = TRUE,
    RUN_MILOPY = FALSE,
    EXIT_AFTER_L3_OFA = TRUE,
    L3_TYPES_INCLUDE = if (exists("L3_TYPES_INCLUDE")) L3_TYPES_INCLUDE else NULL,
    N_CORES = 2L,
    OFA_MAX_CELLS_PER_IDENT = 3000L,
    L3_OFA_MAX_CELLS_PER_IDENT = 3000L,
    PIPELINE_CHANGELOG_LINES = c(
      "  [L3-OFA-ONLY-1] Slim rerun that skips visualization, pairwise DE, ssGSEA, clustering, and Milo.",
      "  [L3-OFA-ONLY-2] Runs per-L3 vs rest OFA, same-L2 sibling OFA, and same-L3 cross-tissue OFA outputs.",
      "  [L3-OFA-ONLY-3] Keeps enrichment dotplots under reports/l3_ofa for requested bubble-style summaries."
    )
  )
)

tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
