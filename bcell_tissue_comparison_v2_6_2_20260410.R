#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.2 Wrapper
# ==============================================================================
#
# Final consolidated wrapper:
#   - reads / summarizes previous outputs before rerun
#   - requires LLM execution
#   - keeps pairwise / ssGSEA / CHOIR LLM outputs integrated and compact
#   - excludes non-informative top pathways (energy / MT / RP / ENSG-like)
#   - excludes cilia-related genes for non-epithelial interpretation contexts
#   - keeps IG evidence only for B-cell lineage; non-B lineages should filter IG
#   - remains compatible with historical output layouts via advanced helper fallbacks
#
# Date: 2026-04-10
# ==============================================================================

LINEAGE_TAG           <- "BCELL"
LINEAGE_DISPLAY       <- "B Cell"
LINEAGE_CONTEXT_LABEL <- "B cells and plasma cells"
LINEAGE_CONTEXT_LOWER <- "B-cell"

PIPELINE_VERSION_LABEL <- "v2.6.2"
PIPELINE_SUBTITLE <- paste(
  "B Cell Tissue Comparison v2.6.2",
  "(final compact integrated LLM, legacy-output compatible)"
)
REPORT_TITLE <- "# B Cell Tissue Comparison Report (Normal Respiratory Tract)"
GENERATED_BY_LABEL <- "bcell_tissue_comparison_v2_6_2_20260410.R"
FINAL_FILE_PREFIX  <- "bcell_tissue_comparison_final"
LINEAGE_COMPLETION_BANNER <- "B CELL TISSUE COMPARISON COMPLETE (v2.6.2)"

H5AD_PATH  <- "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0408/bcell_tissue_comparison_v2_6_20260408"
PREVIOUS_OUTPUT_DIR <- OUTPUT_DIR
REUSE_PREVIOUS_FINAL_OBJECT <- TRUE
REUSE_PREVIOUS_OUTPUT_SUMMARY <- TRUE

ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
source(ADVANCED_HELPER_PATH)

invisible(lapply(c("/home/h2048/.env", "/home/h2048/script/.env"), function(path) {
  if (file.exists(path)) tc_load_env_file(path)
}))

if (nchar(Sys.getenv("DEEPSEEK_API_KEY", unset = "")) < 10) {
  if (!file.exists("/home/h2048/.env")) {
    tc_ensure_env_placeholder("/home/h2048/.env", "DEEPSEEK_API_KEY")
  }
  stop("DEEPSEEK_API_KEY not found in environment or loaded .env; final v2.6.2 pipeline requires LLM and will not continue without it.")
}

SSGSEA_METHODS       <- c("hallmark", "go_bp")
SSGSEA_CHOIR_METHODS <- c("hallmark", "go_bp")

FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT <- FALSE
FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
DEPRIORITIZE_ENERGY_PATHWAYS <- TRUE
DEPRIORITIZE_TECHNICAL_PATHWAYS <- TRUE
DEPRIORITIZE_CILIA_PATHWAYS <- TRUE
APPEND_ENERGY_PATHWAYS_AFTER_TOP <- FALSE

PIPELINE_CHANGELOG_TITLE <- sprintf("%s Changes:", PIPELINE_VERSION_LABEL)
PIPELINE_CHANGELOG_LINES <- c(
  "  [FINAL-1] Compact LLM outputs only: no split final conclusions by up/down or by database.",
  "  [FINAL-2] Non-informative pathways (energy / mitochondrial / ribosomal / ENSG-like) do not occupy top-pathway slots.",
  "  [FINAL-3] Non-epithelial cilia-related genes/pathways are suppressed from enrichment and LLM emphasis.",
  "  [FINAL-4] Previous output reading remains compatible with legacy table/RDS layouts.",
  "  [FINAL-5] Previous output summary is collected before rerun and written back to reports/previous_run_summary.tsv."
)

invisible(tc_apply_advanced_shared_overrides(
  overrides = list(
    INTERPRET_MULTI_DB_MIN_TERMS = 1L,
    INTERPRET_MULTI_DB_MAX_DBS = 5L,
    INTERPRET_MULTI_DB_TERMS_PER_DB = 4L,
    INTERPRET_SSGSEA_TERMS_PER_DB = 6L,
    LLM_INCLUDE_TOP_DEG = TRUE,
    LLM_TOP_DEG_N = 8L,
    LLM_DEG_PADJ_THR = 0.10,
    LLM_DEG_LFC_THR = 0.15,
    LLM_SSGSEA_TERMS_PER_DIRECTION = 5L,
    LLM_REQUIRE_INTEGRATED_UP_DOWN = TRUE,
    LLM_SSGSEA_PRIMARY_USE = "cell_type_judgment",
    LLM_ALLOW_COMPARATIVE_HYPOTHESIS = TRUE,
    FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    DEPRIORITIZE_ENERGY_PATHWAYS = TRUE,
    LLM_EXTRA_RULES = c(
      "Keep each field compact and synthesis-first; never split the final answer into separate up/down or per-database mini-conclusions.",
      "Do not let energy metabolism, mitochondrial, ribosomal/translation, or ENSG-like placeholder pathways occupy top-pathway slots.",
      "For non-epithelial interpretation contexts, exclude cilia/ciliogenesis-related genes and pathways from top-gene emphasis and pathway reasoning.",
      "For non-B-cell interpretation contexts, immunoglobulin genes should not be used as primary evidence.",
      "Prefer lineage-informative biology over housekeeping-like programs when summarizing top pathways."
    )
  )
))

source("/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R")
