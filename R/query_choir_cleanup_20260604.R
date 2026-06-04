#!/usr/bin/env Rscript
# ==============================================================================
# Query CHOIR/Leiden Clustering + OFA + LLM Cleanup — Unified Wrapper
# ==============================================================================
# Date: 2026-06-04
#
# Minimal R pipeline for processed query h5ad files:
#   Phase 1: CHOIR clustering (with Leiden fallback)
#   Phase 2: Cluster-level OFA (one-vs-rest DE + enrichment)
#   Phase 3: LLM batch screening → discovery_screen.tsv
#
# Accepts command-line arguments so a single wrapper serves all lineages.
# Reuses the generic tissue-comparison wrapper + per-lineage shared engines.
#
# Usage (via Bash orchestrator):
#   Rscript query_choir_cleanup_20260604.R \
#     --lineage-tag BCELL \
#     --h5ad-path /path/to/query.h5ad \
#     --output-dir /path/to/output \
#     --shared-engine /path/to/shared_engine.R
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ── CLI arguments ────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--lineage-tag", type = "character", default = "BCELL",
              help = "Lineage tag matching tc_get_lineage_preset(): BCELL, EPITHELIAL, TNK, MYELOID, ENDOTHELIAL, FIBROBLAST, SMC"),
  make_option("--h5ad-path", type = "character", default = "",
              help = "Path to input query h5ad file"),
  make_option("--output-dir", type = "character", default = "",
              help = "Output directory for R pipeline artifacts"),
  make_option("--shared-engine", type = "character", default = "",
              help = "Path to lineage-specific shared engine R script"),
  make_option("--cluster-backend", type = "character", default = "CHOIR",
              help = "Clustering backend: CHOIR or LEIDEN"),
  make_option("--leiden-resolution", type = "double", default = 0.8,
              help = "Leiden resolution (only used if cluster-backend=LEIDEN)"),
  make_option("--n-cores", type = "integer", default = 2L,
              help = "Number of cores for parallel steps"),
  make_option("--force-rerun", action = "store_true", default = FALSE,
              help = "Force re-run even if outputs exist")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

LINEAGE_TAG      <- args[["lineage-tag"]]
H5AD_PATH        <- args[["h5ad-path"]]
OUTPUT_DIR       <- args[["output-dir"]]
SHARED_ENGINE    <- args[["shared-engine"]]
CLUSTER_BACKEND  <- args[["cluster-backend"]]
LEIDEN_RES       <- args[["leiden-resolution"]]
N_CORES          <- args[["n-cores"]]
FORCE_RERUN      <- args[["force-rerun"]]

# Derive file prefix from lineage tag
lineage_lower <- tolower(LINEAGE_TAG)
FINAL_FILE_PREFIX <- paste0(lineage_lower, "_tissue_comparison_final")

# ── Validation ───────────────────────────────────────────────────────────────

if (!nzchar(H5AD_PATH)) stop("--h5ad-path is required")
if (!nzchar(OUTPUT_DIR)) stop("--output-dir is required")
if (!nzchar(SHARED_ENGINE)) stop("--shared-engine is required")
if (!file.exists(H5AD_PATH)) stop(sprintf("H5AD not found: %s", H5AD_PATH))
if (!file.exists(SHARED_ENGINE)) stop(sprintf("Shared engine not found: %s", SHARED_ENGINE))

# ── Source infrastructure ────────────────────────────────────────────────────

Sys.setenv(RETICULATE_PYTHON = "/home/h2048/miniconda3/envs/scvi_env/bin/python")

GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
if (!file.exists(GENERIC_WRAPPER_PATH)) {
  stop(sprintf("Generic wrapper not found: %s", GENERIC_WRAPPER_PATH))
}
source(GENERIC_WRAPPER_PATH)

PARALLEL_LLM_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260508_parallel_llm.R"
if (file.exists(PARALLEL_LLM_PATH)) {
  source(PARALLEL_LLM_PATH)
} else {
  source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R")
}

RECOVERY_LOCK_PATH <- "/home/h2048/script/R/tissue_comparison_recovery_lock_20260512.R"
if (file.exists(RECOVERY_LOCK_PATH)) {
  source(RECOVERY_LOCK_PATH)
}

# ── Resume safety ────────────────────────────────────────────────────────────

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

FINAL_RDS <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".rds"))
FINAL_H5AD <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".h5ad"))
LLM_SCREEN_TSV <- file.path(OUTPUT_DIR, "reports", "llm_choir_discovery_screen.tsv")

if (!FORCE_RERUN && file.exists(FINAL_RDS) && file.exists(FINAL_H5AD) && file.exists(LLM_SCREEN_TSV)) {
  cat(sprintf("[query_choir] Outputs already exist in %s — skipping (use --force-rerun to override)\n", OUTPUT_DIR))
  quit(save = "no", status = 0)
}

# ── Cleanup stale cluster artifacts for fresh run ────────────────────────────

backend_lower <- tolower(CLUSTER_BACKEND)
report_dir <- file.path(OUTPUT_DIR, "reports")
figure_dir <- file.path(OUTPUT_DIR, "figures")

stale_paths <- c(
  file.path(report_dir, backend_lower),
  file.path(report_dir, sprintf("%s_llm_structured.tsv", backend_lower)),
  file.path(report_dir, sprintf("%s_llm_structured_all.rds", backend_lower)),
  file.path(report_dir, sprintf("llm_%s_discovery_screen.tsv", backend_lower)),
  file.path(OUTPUT_DIR, sprintf("LLM_%s_INTERPRETATION.md", toupper(backend_lower))),
  file.path(OUTPUT_DIR, sprintf("LLM_%s_DISCOVERY_REVIEW.md", toupper(backend_lower)))
)
unlink(stale_paths, recursive = TRUE, force = TRUE)

stale_glob <- c(
  Sys.glob(file.path(figure_dir, sprintf("%s_umap*", backend_lower))),
  Sys.glob(file.path(figure_dir, sprintf("ssgsea_%s_heatmap*", backend_lower)))
)
unlink(stale_glob, force = TRUE)

# Final outputs
unlink(c(FINAL_RDS, FINAL_H5AD, file.path(OUTPUT_DIR, "REPORT.md")), force = TRUE)

# ── Recovery lock ────────────────────────────────────────────────────────────

lock_name <- sprintf(".%s_query_cleanup_20260604.lock", lineage_lower)
completion_paths <- c(FINAL_RDS, FINAL_H5AD, LLM_SCREEN_TSV)

if (exists("tc_recovery_register_lock")) {
  .recovery_lock_release <- tc_recovery_register_lock(
    output_dir = OUTPUT_DIR,
    lock_name = lock_name,
    completion_paths = completion_paths,
    label = sprintf("%s query CHOIR cleanup", LINEAGE_TAG)
  )
  on.exit(.recovery_lock_release(), add = TRUE)
}

# ── Build pipeline config ────────────────────────────────────────────────────

# Determine which presets are available for this lineage.
# MYELOID and stromal lineages don't have their own presets in the base
# wrapper; fall back to BCELL which has the most general engine interface.
lineage_for_preset <- switch(
  LINEAGE_TAG,
  ENDOTHELIAL = "STROMAL_ENDOTHELIAL",
  FIBROBLAST  = "STROMAL_FIBROBLAST",
  SMC         = "STROMAL_SMC",
  MYELOID     = "BCELL",       # no myeloid preset; BCELL provides fallback
  "BCELL"                      # default safe fallback
)

cat(sprintf("[query_choir] Building config for lineage=%s (preset key=%s)\n",
            LINEAGE_TAG, lineage_for_preset))

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = lineage_for_preset,
  overrides = list(
    # ── Paths ──
    H5AD_PATH = H5AD_PATH,
    OUTPUT_DIR = OUTPUT_DIR,
    PREVIOUS_OUTPUT_DIR = OUTPUT_DIR,   # same-dir resume-safe
    SHARED_ENGINE_PATH = SHARED_ENGINE,

    # ── Version labels ──
    PIPELINE_VERSION_LABEL = "v1.0-query-cleanup",
    GENERATED_BY_LABEL = "query_choir_cleanup_20260604.R",
    FINAL_FILE_PREFIX = FINAL_FILE_PREFIX,
    LINEAGE_COMPLETION_BANNER = sprintf("%s QUERY CLEANUP COMPLETE", toupper(LINEAGE_TAG)),

    # ── Stage control: minimal clustering+OFA+LLM, skip expensive DE ──
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = FALSE,
    REQUIRE_PREVIOUS_OUTPUT_DIR = FALSE,
    ALLOW_IN_PLACE_RESUME = TRUE,
    REQUIRE_LLM = TRUE,

    RUN_VISUALIZATION = FALSE,
    RUN_PAIRWISE_DE = FALSE,
    RUN_WILCOX = FALSE,
    RUN_SSGSEA = FALSE,
    RUN_CLUSTERING = TRUE,
    RUN_CHOIR = (CLUSTER_BACKEND == "CHOIR"),
    RUN_LEIDEN = (CLUSTER_BACKEND == "LEIDEN"),
    CLUSTER_BACKEND = CLUSTER_BACKEND,
    RUN_CLUSTER_SSGSEA = TRUE,
    RUN_OFA = TRUE,
    RUN_L3_OFA = FALSE,
    RUN_L3_OFA_VS_REST = FALSE,
    RUN_L3_OFA_SAME_L2 = FALSE,
    RUN_L3_OFA_INTER_TISSUE = FALSE,
    ENABLE_LLM = TRUE,

    # ── CHOIR lightweight parameters ──
    N_CORES = as.integer(N_CORES),
    CHOIR_N_CORES = as.integer(N_CORES),
    CHOIR_N_ITERATIONS = 25L,
    CHOIR_N_TREES = 20L,
    CHOIR_VAR_FEATURES_MAX = 2000L,
    CHOIR_SAMPLE_MAX = 5000L,
    CHOIR_DOWNSAMPLING_RATE = 0.05,
    CHOIR_MAX_CLUSTERS = "auto",
    CHOIR_MIN_CLUSTER_DEPTH = 2500L,
    CHOIR_MAX_REPEAT_ERRORS = 2L,
    CHOIR_SUBTREE_REDUCTIONS = FALSE,
    CHOIR_DISTANCE_AWARENESS = 1L,
    OFA_MAX_CELLS_PER_IDENT = 2000L,

    # ── Leiden fallback parameters ──
    LEIDEN_RESOLUTION = as.numeric(LEIDEN_RES),
    LEIDEN_N_DIMS = 30L,
    LEIDEN_K_PARAM = 30L,

    # ── LLM batch screening ──
    ADVANCED_HELPER_PATH = if (file.exists(PARALLEL_LLM_PATH)) PARALLEL_LLM_PATH
                          else "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R",
    LLM_SCREEN_PARALLEL_WORKERS = as.integer(Sys.getenv("LLM_SCREEN_PARALLEL_WORKERS", "2")),
    LLM_SCREEN_PARALLEL_STAGGER_SEC = as.numeric(Sys.getenv("LLM_SCREEN_PARALLEL_STAGGER_SEC", "0.5")),
    STANDARDIZE_LLM_RETRY_SLEEP_SEC = as.numeric(Sys.getenv("STANDARDIZE_LLM_RETRY_SLEEP_SEC", "1")),

    # ── Lineage-specific display metadata (safe defaults for non-preset lineages) ──
    LINEAGE_DISPLAY = switch(LINEAGE_TAG,
      BCELL="B Cell", TNK="T/NK", EPITHELIAL="Epithelial", MYELOID="Myeloid",
      ENDOTHELIAL="Endothelial", FIBROBLAST="Fibroblast", SMC="SMC", LINEAGE_TAG),
    LINEAGE_CONTEXT_LABEL = switch(LINEAGE_TAG,
      BCELL="B cells and plasma cells", TNK="T, NK, and NKT cells",
      EPITHELIAL="epithelial cells", MYELOID="myeloid cells",
      ENDOTHELIAL="endothelial cells", FIBROBLAST="fibroblast cells",
      SMC="smooth muscle cells and pericytes", LINEAGE_TAG),
    BIOLOGICAL_QUESTION_FRAGMENT = switch(LINEAGE_TAG,
      BCELL="B cell biology, mucosal immunity, or tissue microenvironment",
      TNK="T cell biology, NK surveillance, or tissue adaptation",
      EPITHELIAL="epithelial differentiation, barrier function, or mucosal defense",
      MYELOID="myeloid cell states, tissue macrophage programs, or innate immunity",
      ENDOTHELIAL="vascular zonation, lymphatic drainage, or endothelial activation",
      FIBROBLAST="fibroblast specialization, ECM remodeling, or tissue support",
      SMC="smooth muscle contractility, pericyte coverage, or vascular tone",
      "cellular biology"),
    SSGSEA_METHODS = c("hallmark", "go_bp"),
    SSGSEA_CHOIR_METHODS = c("hallmark", "go_bp")
  )
)

# ── Execute ──────────────────────────────────────────────────────────────────

cat(sprintf("[query_choir] Launching pipeline for %s\n", LINEAGE_TAG))
cat(sprintf("  h5ad:      %s\n", H5AD_PATH))
cat(sprintf("  output:    %s\n", OUTPUT_DIR))
cat(sprintf("  engine:    %s\n", SHARED_ENGINE))
cat(sprintf("  backend:   %s\n", CLUSTER_BACKEND))

result <- tc_run_generic_tissue_comparison(PIPELINE_CONFIG)

cat(sprintf("[query_choir] Pipeline complete for %s\n", LINEAGE_TAG))
cat(sprintf("  final RDS: %s\n", FINAL_RDS))
cat(sprintf("  final H5AD: %s\n", FINAL_H5AD))
cat(sprintf("  LLM screen: %s\n", LLM_SCREEN_TSV))
