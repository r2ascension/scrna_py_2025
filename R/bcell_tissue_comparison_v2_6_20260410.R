#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.1 Wrapper
# ==============================================================================
#
# Purpose:
#   Re-run the B-cell tissue comparison pipeline with:
#   - previous-output inspection from the existing 0408 directory
#   - mandatory LLM execution
#   - integrated pairwise up/down evidence for comparison LLM
#   - multi-database ssGSEA aggregation
#   - energy-pathway exclusion from top-pathway selection
#   - MT/RP/ENSG/LINC technical-gene exclusion from enrichment + LLM top genes
#
# Date: 2026-04-10
# ==============================================================================

LINEAGE_TAG           <- "BCELL"
LINEAGE_DISPLAY       <- "B Cell"
LINEAGE_CONTEXT_LABEL <- "B cells and plasma cells"
LINEAGE_CONTEXT_LOWER <- "B-cell"

PIPELINE_VERSION_LABEL <- "v2.6.1"
PIPELINE_SUBTITLE <- paste(
  "B Cell Tissue Comparison v2.6.1",
  "(integrated pairwise LLM, grouped ssGSEA LLM, CHOIR cluster LLM)"
)
REPORT_TITLE <- "# B Cell Tissue Comparison Report (Normal Respiratory Tract)"
GENERATED_BY_LABEL <- "bcell_tissue_comparison_v2_6_20260410.R"
FINAL_FILE_PREFIX  <- "bcell_tissue_comparison_final"
LINEAGE_COMPLETION_BANNER <- "B CELL TISSUE COMPARISON COMPLETE (v2.6.1)"

H5AD_PATH  <- "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0408/bcell_tissue_comparison_v2_6_20260408"
PREVIOUS_OUTPUT_DIR <- OUTPUT_DIR

PIPELINE_CHANGELOG_TITLE <- sprintf("%s Changes:", PIPELINE_VERSION_LABEL)
PIPELINE_CHANGELOG_LINES <- c(
  "  [BCELL-1] Reads the existing 0408 output directory before rerun and records a summary.",
  "  [LLM-6] Pairwise tissue-comparison LLM now integrates up/down DEG and multi-database enrichment into one record.",
  "  [PATH-1] Energy metabolism pathways do not occupy top-pathway slots when more lineage-informative pathways are available.",
  "  [GENE-1] MT/RP/ENSG/LINC technical genes are excluded from enrichment input and LLM top-gene emphasis.",
  "  [SSGSEA-1] ssGSEA LLM remains one record per tissue x cell-type group using combined database evidence."
)

source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R")

invisible(lapply(c("/home/h2048/.env", "/home/h2048/script/.env"), function(path) {
  if (file.exists(path)) tc_load_env_file(path)
}))

if (dir.exists(PREVIOUS_OUTPUT_DIR)) {
  cat("=== Previous Output Summary ===\n")
  previous_run <- tc_read_previous_run(
    PREVIOUS_OUTPUT_DIR,
    include_rds = FALSE,
    include_markdown = FALSE
  )
  tc_print_previous_run_summary(previous_run)
  dir.create(file.path(OUTPUT_DIR, "reports"), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(previous_run$summary) && nrow(previous_run$summary) > 0) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fwrite(
        previous_run$summary,
        file.path(OUTPUT_DIR, "reports", "previous_run_summary.tsv"),
        sep = "\t"
      )
    } else {
      utils::write.table(
        previous_run$summary,
        file.path(OUTPUT_DIR, "reports", "previous_run_summary.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = TRUE
      )
    }
  }
  cat("\n")
}

if (nchar(Sys.getenv("DEEPSEEK_API_KEY", unset = "")) < 10) {
  stop("DEEPSEEK_API_KEY not found in environment or loaded .env; LLM pipeline must not be skipped.")
}

SSGSEA_METHODS       <- c("hallmark", "go_bp")
SSGSEA_CHOIR_METHODS <- c("hallmark", "go_bp")

FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
DEPRIORITIZE_ENERGY_PATHWAYS <- TRUE
APPEND_ENERGY_PATHWAYS_AFTER_TOP <- FALSE

invisible(tc_apply_advanced_shared_overrides(
  overrides = list(
    INTERPRET_MULTI_DB_MIN_TERMS = 1L,
    INTERPRET_MULTI_DB_MAX_DBS = 5L,
    INTERPRET_MULTI_DB_TERMS_PER_DB = 4L,
    INTERPRET_SSGSEA_TERMS_PER_DB = 8L,
    LLM_INCLUDE_TOP_DEG = TRUE,
    LLM_TOP_DEG_N = 10L,
    LLM_DEG_PADJ_THR = 0.10,
    LLM_DEG_LFC_THR = 0.15,
    LLM_SSGSEA_TERMS_PER_DIRECTION = 6L,
    LLM_REQUIRE_INTEGRATED_UP_DOWN = TRUE,
    LLM_SSGSEA_PRIMARY_USE = "cell_type_judgment",
    LLM_ALLOW_COMPARATIVE_HYPOTHESIS = TRUE,
    FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    DEPRIORITIZE_ENERGY_PATHWAYS = TRUE,
    APPEND_ENERGY_PATHWAYS_AFTER_TOP = FALSE,
    LLM_EXTRA_RULES = c(
      "For pairwise tissue comparisons, integrate evidence from both directions into one comparative judgment.",
      "Do not let housekeeping-like energy metabolism pathways occupy top-pathway slots when more lineage-informative pathways are available.",
      "Exclude MT-, RPS, RPL, MRPS, MRPL, ENSG, and LINC genes from top-gene emphasis and pathway reasoning.",
      "When ssGSEA evidence comes from multiple databases, summarize them together instead of interpreting each database separately.",
      "For B-cell comparisons, keep lineage-relevant immunoglobulin evidence when present, but do not over-weight technical mitochondrial or ribosomal genes."
    )
  )
))

source("/home/h2048/script/R/bcell_tissue_comparison_v2_6_20260406.R")
