#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.12
# ==============================================================================
#
# New versioned entrypoint for the ranked-GSEA-aligned B-cell tissue comparison
# pipeline.
#
# Design principles:
#   - old versioned files remain untouched
#   - this file is the new preferred entrypoint
#   - grouped ssGSEA is retained for state/cell-type judgment
#   - pairwise pseudobulk comparisons add true ranked GSEA and expose the
#     resulting evidence to the integrated LLM summaries
#
# Date: 2026-05-23
# ==============================================================================

ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
if (!file.exists(ADVANCED_HELPER_PATH)) {
  stop(sprintf("Advanced helper not found: %s", ADVANCED_HELPER_PATH))
}
source(ADVANCED_HELPER_PATH)

env_flag <- function(name, default = FALSE) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) return(isTRUE(default))
  tolower(raw) %in% c("1", "true", "t", "yes", "y", "on")
}

is_placeholder_secret <- function(x) {
  x <- trimws(as.character(x))
  if (!length(x) || !nzchar(x)) return(TRUE)
  lowered <- tolower(x[[1]])
  lowered %in% c(
    "your-deepseek-api-key",
    "your_deepseek_api_key_here",
    "your_deepseek_api_key",
    "your-key",
    "replace_me",
    "changeme"
  ) || grepl("^your[-_a-z]*api[-_a-z]*key", lowered)
}

PIPELINE_CONFIG <- list(
  LINEAGE_TAG = "BCELL",
  LINEAGE_DISPLAY = "B Cell",
  LINEAGE_CONTEXT_LABEL = "B cells and plasma cells",
  LINEAGE_CONTEXT_LOWER = "B-cell",
  BIOLOGICAL_QUESTION_FRAGMENT = "B cell biology, mucosal immunity, or tissue microenvironment",
  PIPELINE_VERSION_LABEL = "v2.6.12",
  PIPELINE_SUBTITLE = paste(
    "B Cell Tissue Comparison v2.6.12",
    "(new versioned entrypoint with ranked GSEA pairwise evidence)"
  ),
  REPORT_TITLE = "# B Cell Tissue Comparison Report (Normal Respiratory Tract)",
  GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_12_20260523.R",
  FINAL_FILE_PREFIX = "bcell_tissue_comparison_final",
  LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.12)",
  H5AD_PATH = "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad",
  OUTPUT_DIR = "/home/h2048/data/R/0523/bcell_tissue_comparison_v2_6_12_20260523",
  PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0411/bcell_tissue_comparison_v2_6_3_20260411",
  REUSE_PREVIOUS_FINAL_OBJECT = TRUE,
  REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
  REQUIRE_PREVIOUS_OUTPUT_DIR = TRUE,
  ALLOW_IN_PLACE_RESUME = FALSE,
  REQUIRE_LLM = TRUE,
  CREATE_ENV_PLACEHOLDER_IF_MISSING = TRUE,
  ENV_FILE_CANDIDATES = c("/home/h2048/.env", "/home/h2048/script/.env"),
  SSGSEA_METHODS = c("hallmark", "go_bp"),
  SSGSEA_CHOIR_METHODS = c("hallmark", "go_bp"),
  FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT = FALSE,
  FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
  FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
  DEPRIORITIZE_ENERGY_PATHWAYS = TRUE,
  DEPRIORITIZE_TECHNICAL_PATHWAYS = TRUE,
  DEPRIORITIZE_CILIA_PATHWAYS = TRUE,
  APPEND_ENERGY_PATHWAYS_AFTER_TOP = FALSE,
  ADVANCED_HELPER_PATH = ADVANCED_HELPER_PATH,
  ADVANCED_HELPER_ALREADY_LOADED = TRUE,
  SKIP_ENV_AUTOLOAD = TRUE,
  SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE = TRUE,
  SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
  SHARED_OVERRIDES = list(
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
    RUN_RANKED_GSEA = TRUE,
    RANKED_GSEA_DATABASES = c("Hallmark", "GO_BP"),
    RANKED_GSEA_MAX_DBS = 2L,
    RANKED_GSEA_TERMS_PER_DB = 4L,
    RANKED_GSEA_PVALUE_CUTOFF = 0.20,
    RANKED_GSEA_BACKEND = if (requireNamespace("fgsea", quietly = TRUE)) "fgsea" else "DOSE",
    LLM_EXTRA_RULES = c(
      "Keep each field compact and synthesis-first; never split the final answer into separate up/down or per-database mini-conclusions.",
      "Do not let energy metabolism, mitochondrial, ribosomal/translation, or ENSG-like placeholder pathways occupy top-pathway slots.",
      "For non-epithelial interpretation contexts, exclude cilia/ciliogenesis-related genes and pathways from top-gene emphasis and pathway reasoning.",
      "For non-B-cell interpretation contexts, immunoglobulin genes should not be used as primary evidence.",
      "Prefer lineage-informative biology over housekeeping-like programs when summarizing top pathways.",
      "For pairwise tissue comparisons, incorporate ranked GSEA evidence together with ORA instead of treating grouped ssGSEA as a mechanistic substitute."
    )
  )
)

PIPELINE_CONFIG$REQUIRE_LLM <- env_flag("PIPELINE_REQUIRE_LLM", PIPELINE_CONFIG$REQUIRE_LLM)
PIPELINE_CONFIG$ALLOW_IN_PLACE_RESUME <- env_flag("PIPELINE_ALLOW_IN_PLACE_RESUME", PIPELINE_CONFIG$ALLOW_IN_PLACE_RESUME)

PIPELINE_CONFIG$PIPELINE_CHANGELOG_TITLE <- sprintf("%s Changes:", PIPELINE_CONFIG$PIPELINE_VERSION_LABEL)
PIPELINE_CONFIG$PIPELINE_CHANGELOG_LINES <- c(
  "  [ARCH-1] New standalone versioned entrypoint; older versioned files remain unchanged.",
  "  [ARCH-2] Keeps explicit preflight, separated OUTPUT_DIR/PREVIOUS_OUTPUT_DIR, and wrapper-owned env loading.",
  "  [GSEA-1] Pairwise pseudobulk comparisons now add true ranked GSEA in addition to ORA.",
  "  [GSEA-2] Integrated pairwise LLM summaries now absorb ranked GSEA NES and leading-edge evidence.",
  "  [GSEA-3] Grouped ssGSEA remains available for cell-type/state judgment and is not repurposed as fgsea-style evidence.",
  "  [FINAL-1] Compact LLM outputs only: no split final conclusions by up/down or by database.",
  "  [FINAL-2] Non-informative pathways (energy / mitochondrial / ribosomal / ENSG-like) do not occupy top-pathway slots.",
  "  [FINAL-3] Non-epithelial cilia-related genes/pathways are suppressed from enrichment and LLM emphasis."
)

tc_assert_named_list_keys(
  PIPELINE_CONFIG,
  required_keys = c(
    "OUTPUT_DIR", "PREVIOUS_OUTPUT_DIR", "H5AD_PATH", "ENV_FILE_CANDIDATES",
    "REUSE_PREVIOUS_FINAL_OBJECT", "REUSE_PREVIOUS_OUTPUT_SUMMARY",
    "REQUIRE_PREVIOUS_OUTPUT_DIR", "ALLOW_IN_PLACE_RESUME",
    "REQUIRE_LLM", "SHARED_ENGINE_PATH", "SHARED_OVERRIDES"
  ),
  object_name = "PIPELINE_CONFIG"
)

PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS <- file.path(
  PIPELINE_CONFIG$PREVIOUS_OUTPUT_DIR,
  paste0(PIPELINE_CONFIG$FINAL_FILE_PREFIX, ".rds")
)

if (!file.exists(PIPELINE_CONFIG$SHARED_ENGINE_PATH)) {
  stop(sprintf("Shared engine not found: %s", PIPELINE_CONFIG$SHARED_ENGINE_PATH))
}

if (!file.exists(PIPELINE_CONFIG$H5AD_PATH) &&
    !(isTRUE(PIPELINE_CONFIG$REUSE_PREVIOUS_FINAL_OBJECT) && file.exists(PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS))) {
  stop(paste(
    sprintf("H5AD input not found: %s", PIPELINE_CONFIG$H5AD_PATH),
    "No reusable previous final object was found either, so the wrapper cannot continue."
  ))
}

if (isTRUE(PIPELINE_CONFIG$REUSE_PREVIOUS_FINAL_OBJECT) && !file.exists(PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS)) {
  cat(sprintf("[WARN] Previous final object requested but not found: %s\n", PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS))
}

loaded_env_files <- tc_load_env_candidates(PIPELINE_CONFIG$ENV_FILE_CANDIDATES)
PIPELINE_CONFIG$PRELOADED_ENV_FILES <- loaded_env_files

llm_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
if (nchar(llm_key) < 10 && isTRUE(PIPELINE_CONFIG$CREATE_ENV_PLACEHOLDER_IF_MISSING)) {
  tc_ensure_env_placeholder(PIPELINE_CONFIG$ENV_FILE_CANDIDATES[[1]], "DEEPSEEK_API_KEY")
}
llm_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
has_llm_key <- nchar(llm_key) >= 10 && !is_placeholder_secret(llm_key)
if (!has_llm_key) {
  if (nchar(llm_key) >= 10 && is_placeholder_secret(llm_key)) {
    Sys.unsetenv("DEEPSEEK_API_KEY")
  }
  llm_message <- paste(
    "DEEPSEEK_API_KEY not found or is still a placeholder value in loaded .env files.",
    "Set REQUIRE_LLM=FALSE if you want a degraded non-LLM run; otherwise provide a valid key."
  )
  if (isTRUE(PIPELINE_CONFIG$REQUIRE_LLM)) stop(llm_message)
  cat(sprintf("[WARN] %s\n", llm_message))
} else if (length(loaded_env_files) > 0) {
  cat(sprintf("[OK] Loaded environment file(s): %s\n", paste(loaded_env_files, collapse = ", ")))
}

summary_output_path <- file.path(PIPELINE_CONFIG$OUTPUT_DIR, "reports", "previous_run_summary.tsv")
preflight <- tc_preflight_previous_run(
  previous_output_dir = PIPELINE_CONFIG$PREVIOUS_OUTPUT_DIR,
  current_output_dir = PIPELINE_CONFIG$OUTPUT_DIR,
  allow_in_place_resume = PIPELINE_CONFIG$ALLOW_IN_PLACE_RESUME,
  require_previous_dir = PIPELINE_CONFIG$REQUIRE_PREVIOUS_OUTPUT_DIR,
  reuse_previous_output_summary = PIPELINE_CONFIG$REUSE_PREVIOUS_OUTPUT_SUMMARY,
  summary_output_path = summary_output_path,
  include_rds = FALSE,
  include_markdown = FALSE,
  verbose = TRUE
)
PIPELINE_CONFIG$PREVIOUS_RUN_SUMMARY_OUTPUT_PATH <- summary_output_path
PIPELINE_CONFIG$PREVIOUS_RUN_SUMMARY_DF <- preflight$summary

dir.create(file.path(PIPELINE_CONFIG$OUTPUT_DIR, "reports"), recursive = TRUE, showWarnings = FALSE)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(
    list(
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      pipeline_version = PIPELINE_CONFIG$PIPELINE_VERSION_LABEL,
      generated_by = PIPELINE_CONFIG$GENERATED_BY_LABEL,
      output_dir = PIPELINE_CONFIG$OUTPUT_DIR,
      previous_output_dir = PIPELINE_CONFIG$PREVIOUS_OUTPUT_DIR,
      previous_output_exists = isTRUE(preflight$exists),
      same_dir = isTRUE(preflight$same_dir),
      require_previous_output_dir = isTRUE(PIPELINE_CONFIG$REQUIRE_PREVIOUS_OUTPUT_DIR),
      require_llm = isTRUE(PIPELINE_CONFIG$REQUIRE_LLM),
      h5ad_exists = file.exists(PIPELINE_CONFIG$H5AD_PATH),
      previous_final_object_rds = PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS,
      previous_final_object_exists = file.exists(PIPELINE_CONFIG$PREVIOUS_FINAL_OBJECT_RDS),
      loaded_env_files = loaded_env_files,
      shared_engine_path = PIPELINE_CONFIG$SHARED_ENGINE_PATH,
      shared_override_names = names(PIPELINE_CONFIG$SHARED_OVERRIDES)
    ),
    path = file.path(PIPELINE_CONFIG$OUTPUT_DIR, "reports", "wrapper_preflight.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
}

tc_apply_named_list(PIPELINE_CONFIG, envir = environment())
invisible(tc_apply_advanced_shared_overrides(overrides = PIPELINE_CONFIG$SHARED_OVERRIDES))

source(SHARED_ENGINE_PATH)
