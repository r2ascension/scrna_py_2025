#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.6.1 - PRODUCTION
# ==============================================================================
#
# Purpose:
#   B cell lineage tissue comparison across normal respiratory tract sites:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. L3 -> L2 remapping (merge Memory/GC/Plasma subtypes)
#     3. Visualization (UMAP, marker dotplot/heatmap, composition)
#     4. Pseudobulk DESeq2 per L2/L3 subtype x tissue pair
#     5. Exploratory cell-level wilcox per L2/L3 subtype (marker discovery only)
#     6. Multi-database enrichment (GO BP/MF/CC, KEGG, Hallmark,
#        CellMarker, PanglaoDB, custom B cell markers) - symbol direct
#     7. interpret_agent (DeepSeek) with tissue-pair-specific context
#     8. Structured REPORT.md (png embeds)
#
# v2.6.1 Changes (vs v2.4.0):
#   [LLM-6] Pairwise comparison LLM now integrates up/down evidence together
#           (including multi-database evidence) instead of splitting by direction
#   [RESUME-1] Can reuse previous B-cell outputs / final object to accelerate reruns
#   [ENV-1] Auto-load DEEPSEEK_API_KEY from common .env files
#
# v2.4.0 Changes (vs v2.3.0):
#   [LLM-4] run_grouped_ssgsea: per tissue x L2/L3 group LLM interpretation
#            after top pathways computed; results in ssgsea_llm_structured.tsv
#            and LLM_SSGSEA_INTERPRETATION.md
#   [LLM-5] CHOIR OFA per-cluster interpret_agent (up + down); results in
#            choir_llm_structured.tsv and LLM_CHOIR_INTERPRETATION.md
#
# v2.3.0 Changes (vs v2.2.2):
#   [PB-1] Global pseudobulk cache per level (L2/L3), avoiding repeated subset()
#          and delaying dense conversion until DESeq2 input construction
#   [LLM-3] Added validated gene->pathway map from real enrichment hits and passed
#          it into both interpret context and LLM standardization evidence
#   [WX-1] Exploratory wilcox outputs now carry explicit evidence-tier metadata and
#          use exploratory_* filenames to reduce downstream misread risk
#   [P2-1] stratified_downsample() now ignores NA groups explicitly
#   [P2-2] h5ad sparse export avoids COO triplet data.frame materialization
#   [P2-3] normalization provenance recorded under obj@misc$normalization_for_report
#
# v2.2.1 Fixes (vs v2.2):
#   [P1-1] UMAP output without raster / alpha; png dpi 300
#   [P1-2] aggregate_pseudobulk: sparse-friendly list-cbind (no dense pre-alloc)
#   [P1-3] Unmapped L3 values -> stop() fail-fast (not silent drop)
#   [P1-4] CellMarker/PanglaoDB: strict B cell / plasma / GC term filter
#   [P1-5] interpret_agent retry + atomic-result-tolerant normalization
#   [P1-6] Secondary DeepSeek standardization to fixed JSON schema
#   [P2-6] Heatmap title annotated "visualization only"
#   [P2-7] set.seed(42) inside stratified_downsample for cross-session reproducibility
#
# Author: r2end + GitHub Copilot
# Date:   2026-04-10
# ==============================================================================

# ==============================================================================
# 0. Thread Control
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS  = "1"
)

# ==============================================================================
# 1. Configuration
# ==============================================================================

if (!exists("LINEAGE_TAG"))               LINEAGE_TAG <- "BCELL"
if (!exists("LINEAGE_DISPLAY"))           LINEAGE_DISPLAY <- "B Cell"
if (!exists("LINEAGE_CONTEXT_LABEL"))     LINEAGE_CONTEXT_LABEL <- "B cells and plasma cells"
if (!exists("LINEAGE_CONTEXT_LOWER"))     LINEAGE_CONTEXT_LOWER <- "B-cell"
if (!exists("BIOLOGICAL_QUESTION_FRAGMENT"))
  BIOLOGICAL_QUESTION_FRAGMENT <- "B cell biology, mucosal immunity, or tissue microenvironment"
if (!exists("PIPELINE_VERSION_LABEL"))    PIPELINE_VERSION_LABEL <- "v2.6.1"
if (!exists("PIPELINE_SUBTITLE"))
  PIPELINE_SUBTITLE <- sprintf("%s Tissue Comparison %s (pseudobulk DESeq2, multi-database enrichment)",
                               LINEAGE_DISPLAY, PIPELINE_VERSION_LABEL)
if (!exists("REPORT_TITLE"))
  REPORT_TITLE <- sprintf("# %s Tissue Comparison Report (Normal Respiratory Tract)", LINEAGE_DISPLAY)
if (!exists("GENERATED_BY_LABEL"))        GENERATED_BY_LABEL <- "bcell_tissue_comparison_v2_6_1_20260410.R"
if (!exists("FINAL_FILE_PREFIX"))         FINAL_FILE_PREFIX <- "bcell_tissue_comparison_final"
if (!exists("LINEAGE_COMPLETION_BANNER"))
  LINEAGE_COMPLETION_BANNER <- sprintf("%s TISSUE COMPARISON COMPLETE (%s)", toupper(LINEAGE_DISPLAY), PIPELINE_VERSION_LABEL)
if (!exists("INITIALIZE_ONLY"))          INITIALIZE_ONLY <- FALSE
if (!exists("USE_EXISTING_L2"))           USE_EXISTING_L2 <- FALSE
if (!exists("L2_SOURCE_COL"))             L2_SOURCE_COL <- "cell_type_L2"
if (!exists("PIPELINE_TEST_MODE")) {
  PIPELINE_TEST_MODE <- identical(tolower(Sys.getenv("PIPELINE_TEST_MODE", "false")), "true")
}
if (!exists("ADVANCED_HELPER_PATH")) {
  ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
}
if (!exists("ADVANCED_HELPER_ALREADY_LOADED")) ADVANCED_HELPER_ALREADY_LOADED <- FALSE
if (!exists("SKIP_ENV_AUTOLOAD")) SKIP_ENV_AUTOLOAD <- FALSE
if (!exists("SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE")) SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE <- FALSE
if (!isTRUE(ADVANCED_HELPER_ALREADY_LOADED) && file.exists(ADVANCED_HELPER_PATH)) {
  source(ADVANCED_HELPER_PATH)
  ADVANCED_HELPER_ALREADY_LOADED <- TRUE
}
if (!exists("PREVIOUS_OUTPUT_DIR")) {
  PREVIOUS_OUTPUT_DIR <- "/home/h2048/data/R/0408/bcell_tissue_comparison_v2_6_20260408"
}
if (!exists("REUSE_PREVIOUS_FINAL_OBJECT")) REUSE_PREVIOUS_FINAL_OBJECT <- TRUE
if (!exists("REUSE_PREVIOUS_OUTPUT_SUMMARY")) REUSE_PREVIOUS_OUTPUT_SUMMARY <- TRUE
if (!exists("CUSTOM_DB_NAME"))            CUSTOM_DB_NAME <- "B_cell_custom"
if (!exists("CUSTOM_DB_LABEL"))           CUSTOM_DB_LABEL <- "B cell custom markers"
if (!exists("CUSTOM_ENRICHMENT_SIZE_RULE"))
  CUSTOM_ENRICHMENT_SIZE_RULE <- c(min = 2L, max = 100L)
if (!exists("CELLMARKER_CELLTYPE_PATTERN")) {
  CELLMARKER_CELLTYPE_PATTERN <- paste(
    "B cell|B-cell|B lymphocyte|Plasma cell|Plasmablast|Germinal center|",
    "Memory B|Naive B|Follicular B|Marginal zone B|Breg|Age-associated B|",
    "Atypical B|Pre-B|Pro-B|Transitional B"
  )
}
if (!exists("PANGLAODB_CELLTYPE_PATTERN")) {
  PANGLAODB_CELLTYPE_PATTERN <- "B cell|Plasma cell|Plasmablast|Germinal center|Memory B|Naive B|Follicular"
}
if (!exists("FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT")) {
  FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT <- !identical(toupper(trimws(as.character(LINEAGE_TAG))), "BCELL")
}
if (!exists("IG_RELATED_GENE_REGEX")) {
  IG_RELATED_GENE_REGEX <- "^IG[HKL]"
}
if (!exists("FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT")) {
  lineage_hint <- paste(
    toupper(trimws(as.character(LINEAGE_TAG))),
    toupper(trimws(as.character(LINEAGE_DISPLAY))),
    toupper(trimws(as.character(LINEAGE_CONTEXT_LABEL))),
    toupper(trimws(as.character(LINEAGE_CONTEXT_LOWER)))
  )
  FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT <- !grepl("EPITHEL|CILIATED", lineage_hint, perl = TRUE)
}
if (!exists("CILIA_RELATED_GENE_REGEX")) {
  CILIA_RELATED_GENE_REGEX <- paste(
    c(
      "^FOXJ1$", "^TPPP3$", "^PIFO$", "^RSPH", "^DNAH", "^DNAI", "^CFAP",
      "^CCNO$", "^MCIDAS$", "^DEUP1$", "^FOXN4$", "^CDC20B$", "^HYDIN$",
      "^DNAAF", "^DYX1C1$", "^LRRC6$", "^RFX2$", "^RFX3$", "^IFT", "^BBS",
      "^ARMC4$", "^SPEF2$", "^CCDC", "^CEP"
    ),
    collapse = "|"
  )
}
if (!exists("FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT")) {
  FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
}
if (!exists("FILTER_NONCODING_LOCUS_GENES_FOR_LLM_AND_ENRICHMENT")) {
  FILTER_NONCODING_LOCUS_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
}
if (!exists("NONCODING_LOCUS_GENE_REGEX")) {
  NONCODING_LOCUS_GENE_REGEX <- if (exists("tc_placeholder_noncoding_gene_regex", mode = "function")) {
    tc_placeholder_noncoding_gene_regex()
  } else {
    paste(
      c(
        "^(AC|AL|AP|BX|Z)[0-9]+\\.",
        "^CT[BCD]-",
        "-OT[0-9]+$",
        "^LOC[0-9]+"
      ),
      collapse = "|"
    )
  }
}
if (!exists("TECHNICAL_GENE_REGEX")) {
  TECHNICAL_GENE_REGEX <- if (exists("tc_default_technical_gene_regex", mode = "function")) {
    tc_default_technical_gene_regex()
  } else {
    paste(
      c(
        "^ENSG",
        "^LINC",
        "^MT-",
        "^MT\\.",
        "^MTRNR",
        "^RPS",
        "^RPL",
        "^MRPS",
        "^MRPL",
        "^RP[0-9]+-",
        "^RP[0-9]+$"
      ),
      collapse = "|"
    )
  }
}
if (!exists("DEPRIORITIZE_ENERGY_PATHWAYS")) {
  DEPRIORITIZE_ENERGY_PATHWAYS <- TRUE
}
if (!exists("DEPRIORITIZE_TECHNICAL_PATHWAYS")) {
  DEPRIORITIZE_TECHNICAL_PATHWAYS <- TRUE
}
if (!exists("DEPRIORITIZE_CILIA_PATHWAYS")) {
  DEPRIORITIZE_CILIA_PATHWAYS <- isTRUE(FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT)
}
if (!exists("APPEND_ENERGY_PATHWAYS_AFTER_TOP")) {
  APPEND_ENERGY_PATHWAYS_AFTER_TOP <- FALSE
}
if (!exists("ENERGY_PATHWAY_REGEX")) {
  ENERGY_PATHWAY_REGEX <- paste(
    c(
      "oxidative phosphorylation",
      "electron transport",
      "respiratory chain",
      "cellular respiration",
      "aerobic respiration",
      "atp synthesis",
      "atp metabolic",
      "mitochondrial atp",
      "generation of precursor metabolites and energy",
      "tricarboxylic acid cycle",
      "tca cycle",
      "citric acid cycle",
      "glycolysis",
      "gluconeogenesis",
      "fatty acid beta[- ]oxidation",
      "fatty acid oxidation"
    ),
    collapse = "|"
  )
}
if (!exists("TECHNICAL_PATHWAY_REGEX")) {
  TECHNICAL_PATHWAY_REGEX <- paste(
    c(
      "mitochond", "mitophag", "mitochondrial gene expression", "mitochondrial transl",
      "mt[- _.]", "mtrnr", "cytosolic ribosome", "ribosom", "ribosomal",
      "translation", "translational", "rRNA processing", "rRNA metabolic",
      "peptide biosynthetic", "structural constituent of ribosome",
      "SRP-dependent cotranslational protein targeting", "\\bENSG"
    ),
    collapse = "|"
  )
}
if (!exists("CILIA_PATHWAY_REGEX")) {
  CILIA_PATHWAY_REGEX <- paste(
    c(
      "cilium", "ciliary", "ciliogenesis", "motile cilium", "axoneme",
      "dynein complex", "intraflagellar transport", "flagell", "microtubule-based movement"
    ),
    collapse = "|"
  )
}

# ----- Input / Output -----
if (!exists("H5AD_PATH")) {
  H5AD_PATH <- "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
}
if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- PREVIOUS_OUTPUT_DIR
}
if (!exists("PREVIOUS_FINAL_OBJECT_RDS")) {
  PREVIOUS_FINAL_OBJECT_RDS <- file.path(PREVIOUS_OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".rds"))
}

# ----- Environment Files -----
ENV_FILE_CANDIDATES <- c("/home/h2048/.env", "/home/h2048/script/.env")

load_env_file <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  lines <- readLines(path, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- trimws(gsub("^export\\s+", "", sub("=.*$", "", line)))
    val <- trimws(gsub("^['\"]|['\"]$", "", sub("^[^=]*=", "", line)))
    if (nzchar(key) && !nzchar(Sys.getenv(key, unset = ""))) {
      Sys.setenv(structure(val, names = key))
    }
  }
  invisible(TRUE)
}

env_files_loaded <- ENV_FILE_CANDIDATES[file.exists(ENV_FILE_CANDIDATES)]
if (!isTRUE(SKIP_ENV_AUTOLOAD)) {
  invisible(lapply(env_files_loaded, load_env_file))
} else if (exists("PRELOADED_ENV_FILES")) {
  env_files_loaded <- PRELOADED_ENV_FILES[file.exists(PRELOADED_ENV_FILES)]
} else {
  env_files_loaded <- character()
}

# ----- L3 -> L2 Remapping -----
if (!exists("L3_SOURCE_COL")) L3_SOURCE_COL <- "cell_type_scanvi_pred"
if (!exists("ANALYSIS_L3_DESCRIPTION")) ANALYSIS_L3_DESCRIPTION <- "fine scanvi labels"
if (!exists("L3_TO_L2_TABLE_HEADER_LEFT")) {
  L3_TO_L2_TABLE_HEADER_LEFT <- sprintf("L3 (`%s`)", L3_SOURCE_COL)
}

if (!exists("L3_TO_L2_REMAP")) {
  L3_TO_L2_REMAP <- c(
    "Atypical_Memory_B"                  = "Memory_B",
    "IGHEplus_Atypical_Memory_B"         = "Memory_B",
    "Memory_B"                           = "Memory_B",
    "GC_B_Dark_Zone_Centroblast_Cycling" = "GC_B",
    "GC_B_Light_Zone_Centrocyte"         = "GC_B",
    "GC_B_Transitional"                  = "GC_B",
    "Plasma_IgA"                         = "Plasma",
    "Plasma_IgG"                         = "Plasma",
    "Naive_B"                            = "Naive_B"
  )
}

# ----- Reference Database Paths -----
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL      <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH  <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"

# ----- DeepSeek API -----
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
ENABLE_LLM <- nchar(DEEPSEEK_API_KEY) >= 10
if (!ENABLE_LLM) {
  cat("[WARN] DEEPSEEK_API_KEY not set. interpret_agent will be SKIPPED.\n")
  cat("       Set via: export DEEPSEEK_API_KEY='your-key' to enable LLM.\n")
} else if (length(env_files_loaded) > 0) {
  cat(sprintf("[OK] Loaded environment file(s): %s\n", paste(env_files_loaded, collapse = ", ")))
}
INTERPRET_AGENT_MODEL           <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS      <- 50
INTERPRET_AGENT_ADD_PPI         <- FALSE
INTERPRET_AGENT_MAX_RETRIES     <- 3L
INTERPRET_AGENT_RETRY_SLEEP_SEC <- 2
STANDARDIZE_LLM_OUTPUT          <- TRUE
STANDARDIZE_LLM_MODEL           <- "deepseek-chat"
STANDARDIZE_LLM_MAX_RETRIES     <- 3L
STANDARDIZE_LLM_RETRY_SLEEP_SEC <- 2
STANDARDIZE_LLM_MAX_INPUT_CHARS <- 12000L
if (!exists("INTERPRET_MULTI_DB_PREFERRED")) {
  INTERPRET_MULTI_DB_PREFERRED <- c(
    "GO_BP", "Hallmark", "KEGG", "CellMarker",
    CUSTOM_DB_NAME, "PanglaoDB", "GO_MF", "GO_CC"
  )
}
if (!exists("INTERPRET_MULTI_DB_MAX_DBS"))       INTERPRET_MULTI_DB_MAX_DBS <- 5L
if (!exists("INTERPRET_MULTI_DB_TERMS_PER_DB")) INTERPRET_MULTI_DB_TERMS_PER_DB <- 4L
if (!exists("INTERPRET_MULTI_DB_MIN_TERMS"))    INTERPRET_MULTI_DB_MIN_TERMS <- 2L
if (!exists("INTERPRET_SSGSEA_TERMS_PER_DB"))   INTERPRET_SSGSEA_TERMS_PER_DB <- 8L
if (!exists("LLM_INCLUDE_TOP_DEG"))              LLM_INCLUDE_TOP_DEG <- TRUE
if (!exists("LLM_TOP_DEG_N"))                    LLM_TOP_DEG_N <- 10L
if (!exists("LLM_DEG_PADJ_THR"))                 LLM_DEG_PADJ_THR <- 0.10
if (!exists("LLM_DEG_LFC_THR"))                  LLM_DEG_LFC_THR <- 0.15
if (!exists("LLM_SSGSEA_TERMS_PER_DIRECTION"))   LLM_SSGSEA_TERMS_PER_DIRECTION <- 4L
if (!exists("LLM_REQUIRE_INTEGRATED_UP_DOWN"))   LLM_REQUIRE_INTEGRATED_UP_DOWN <- TRUE
if (!exists("LLM_SSGSEA_PRIMARY_USE"))           LLM_SSGSEA_PRIMARY_USE <- "cell_type_judgment"
if (!exists("LLM_ALLOW_COMPARATIVE_HYPOTHESIS")) LLM_ALLOW_COMPARATIVE_HYPOTHESIS <- TRUE
if (!exists("LLM_EXTRA_RULES"))                  LLM_EXTRA_RULES <- character()
if (!exists("CHOIR_LLM_TOP_DEG_N"))               CHOIR_LLM_TOP_DEG_N <- 20L
if (!exists("CHOIR_LLM_TOP_ENRICH_TERMS_PER_DB")) CHOIR_LLM_TOP_ENRICH_TERMS_PER_DB <- 6L
if (!exists("CHOIR_LLM_TOP_SSGSEA_TERMS_PER_DIRECTION")) CHOIR_LLM_TOP_SSGSEA_TERMS_PER_DIRECTION <- 8L
if (!exists("CHOIR_LLM_MAX_DBS"))                 CHOIR_LLM_MAX_DBS <- 8L

# ----- Python / reticulate -----
PYTHON_CONDA_ENV <- "bbknn_env"

# ----- Metadata Column Names -----
TISSUE_COL      <- "tissue"
SAMPLE_COL      <- "sample"
CELLTYPE_L2_COL <- "cell_type_L2"
CELLTYPE_L3_COL <- "cell_type_L3"
LABEL_COL       <- CELLTYPE_L3_COL

UMAP_REDUCTION_PREFERRED <- c("umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap")
UMAP_PT_SIZE <- 0.35
UMAP_TISSUE_COLORS <- c(
  "lung parenchyma"    = "#D55E00",
  "nose"               = "#009E73",
  "respiratory airway" = "#0072B2",
  "sinus"              = "#CC79A7"
)

# ----- Pseudobulk DE Parameters -----
MIN_CELLS_PER_PSEUDOBULK <- 10
MIN_SAMPLES_PER_TISSUE   <- 3
PADJ_THR                 <- 0.05
LFC_THR                  <- 1.0
DESIGN_COVARIATE         <- NULL
if (!exists("PSEUDOBULK_DISAMBIGUATION_CANDIDATES")) {
  PSEUDOBULK_DISAMBIGUATION_CANDIDATES <- c("dataset", "source", "batch", "orig.ident")
}

# ----- Reproducibility -----
set.seed(42)

current_future_maxsize <- getOption("future.globals.maxSize")
if (is.null(current_future_maxsize) || !is.numeric(current_future_maxsize) || !is.finite(current_future_maxsize)) {
  current_future_maxsize <- 0
}

# ----- Exploratory Wilcox Parameters -----
WILCOX_MIN_CELLS <- 50
WILCOX_LFC_THR   <- 0.25
WILCOX_TOP_N     <- 200

# ----- Enrichment -----
TOP_N_DEG_ENRICHMENT <- 200

if (!exists("ENRICHMENT_GS_SIZE_RULES")) {
  ENRICHMENT_GS_SIZE_RULES <- list(
    default    = c(min = 10L, max = 500L),
    CellMarker = c(min = 3L,  max = 200L),
    PanglaoDB  = c(min = 3L,  max = 200L)
  )
  ENRICHMENT_GS_SIZE_RULES[[CUSTOM_DB_NAME]] <- CUSTOM_ENRICHMENT_SIZE_RULE
}

# ----- Visualization -----
HEATMAP_CELLS_PER_TYPE <- 100
if (!exists("RUN_VISUALIZATION")) RUN_VISUALIZATION <- TRUE

# ----- Stage gating / recovery -----
if (!exists("RUN_PAIRWISE_DE")) RUN_PAIRWISE_DE <- TRUE
if (!exists("RUN_WILCOX")) RUN_WILCOX <- TRUE
if (!exists("LOAD_EXISTING_STAGE_ARTIFACTS")) LOAD_EXISTING_STAGE_ARTIFACTS <- TRUE
if (!exists("LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT")) LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT <- TRUE

# ----- ssGSEA -----
if (!exists("RUN_SSGSEA")) RUN_SSGSEA <- TRUE
if (!exists("SSGSEA_METHODS")) SSGSEA_METHODS <- c("hallmark", "go_bp")
if (!exists("SSGSEA_N_TOP")) SSGSEA_N_TOP <- 20
if (!exists("SSGSEA_N_HEATMAP")) SSGSEA_N_HEATMAP <- 30
if (!exists("SSGSEA_CUSTOM_GMT")) SSGSEA_CUSTOM_GMT <- NULL
if (!exists("N_CORES")) N_CORES <- 4

# ----- Clustering backend -----
if (!exists("RUN_CLUSTERING")) {
  if (exists("RUN_CHOIR")) {
    RUN_CLUSTERING <- isTRUE(RUN_CHOIR)
  } else {
    RUN_CLUSTERING <- TRUE
  }
}
RUN_CLUSTERING <- isTRUE(RUN_CLUSTERING)
if (!exists("RUN_CHOIR")) RUN_CHOIR <- RUN_CLUSTERING
RUN_CHOIR <- isTRUE(RUN_CHOIR)
if (!exists("CHOIR_N_CORES")) CHOIR_N_CORES <- N_CORES
if (!exists("CLUSTER_BACKEND")) CLUSTER_BACKEND <- "CHOIR"
CLUSTER_BACKEND <- toupper(trimws(as.character(CLUSTER_BACKEND)))
if (!CLUSTER_BACKEND %in% c("CHOIR", "LEIDEN")) {
  stop(sprintf("Unsupported CLUSTER_BACKEND: %s (expected CHOIR or LEIDEN)", CLUSTER_BACKEND))
}
if (!exists("LEIDEN_RESOLUTION")) LEIDEN_RESOLUTION <- 0.8
if (!exists("LEIDEN_N_DIMS")) LEIDEN_N_DIMS <- 30L
if (!exists("LEIDEN_K_PARAM")) LEIDEN_K_PARAM <- 30L
if (!exists("LEIDEN_ALGORITHM")) LEIDEN_ALGORITHM <- 4L
if (!exists("CHOIR_REDUCTION_CANDIDATES")) CHOIR_REDUCTION_CANDIDATES <- c("scanvi", "scvi", "harmony", "pca")
if (!exists("CHOIR_ALPHA")) CHOIR_ALPHA <- 0.20
if (!exists("CHOIR_VAR_FEATURES_MAX")) CHOIR_VAR_FEATURES_MAX <- 4000L
if (!exists("CHOIR_FIND_VAR_FEATURES_IF_MISSING")) CHOIR_FIND_VAR_FEATURES_IF_MISSING <- TRUE
if (!exists("CHOIR_VAR_FEATURES_METHOD")) CHOIR_VAR_FEATURES_METHOD <- "vst"
if (!exists("CHOIR_N_ITERATIONS")) CHOIR_N_ITERATIONS <- 100L
if (!exists("CHOIR_N_TREES")) CHOIR_N_TREES <- 50L
if (!exists("CHOIR_USE_VARIANCE")) CHOIR_USE_VARIANCE <- TRUE
if (!exists("CHOIR_MIN_ACCURACY")) CHOIR_MIN_ACCURACY <- 0.5
if (!exists("CHOIR_MIN_CONNECTIONS")) CHOIR_MIN_CONNECTIONS <- 1L
if (!exists("CHOIR_MAX_REPEAT_ERRORS")) CHOIR_MAX_REPEAT_ERRORS <- 20L
if (!exists("CHOIR_DISTANCE_APPROX")) CHOIR_DISTANCE_APPROX <- TRUE
if (!exists("CHOIR_DISTANCE_AWARENESS")) CHOIR_DISTANCE_AWARENESS <- 2L
if (!exists("CHOIR_COLLECT_ALL_METRICS")) CHOIR_COLLECT_ALL_METRICS <- FALSE
if (!exists("CHOIR_SAMPLE_MAX")) CHOIR_SAMPLE_MAX <- Inf
if (!exists("CHOIR_DOWNSAMPLING_RATE")) CHOIR_DOWNSAMPLING_RATE <- "auto"
if (!exists("CHOIR_MIN_READS")) CHOIR_MIN_READS <- NULL
if (!exists("CHOIR_MAX_CLUSTERS")) CHOIR_MAX_CLUSTERS <- "auto"
if (!exists("CHOIR_MIN_CLUSTER_DEPTH")) CHOIR_MIN_CLUSTER_DEPTH <- 2000L
if (!exists("CHOIR_NORMALIZATION_METHOD")) CHOIR_NORMALIZATION_METHOD <- "none"
if (!exists("CHOIR_SUBTREE_REDUCTIONS")) CHOIR_SUBTREE_REDUCTIONS <- TRUE
if (!exists("CHOIR_REDUCTION_METHOD")) CHOIR_REDUCTION_METHOD <- NULL
if (!exists("CHOIR_REDUCTION_PARAMS")) CHOIR_REDUCTION_PARAMS <- list()
if (!exists("CHOIR_BATCH_CORRECTION_METHOD")) CHOIR_BATCH_CORRECTION_METHOD <- "none"
if (!exists("CHOIR_BATCH_CORRECTION_PARAMS")) CHOIR_BATCH_CORRECTION_PARAMS <- list()
if (!exists("CHOIR_BATCH_LABELS")) CHOIR_BATCH_LABELS <- NULL
if (!exists("CHOIR_NEIGHBOR_PARAMS")) CHOIR_NEIGHBOR_PARAMS <- list(verbose = FALSE)
if (!exists("CHOIR_CLUSTER_PARAMS")) CHOIR_CLUSTER_PARAMS <- list(algorithm = 1, group.singletons = TRUE, verbose = FALSE)
if (!exists("CHOIR_FUTURE_GLOBALS_MAXSIZE")) CHOIR_FUTURE_GLOBALS_MAXSIZE <- 8 * 1024^3
options(future.globals.maxSize = max(current_future_maxsize, CHOIR_FUTURE_GLOBALS_MAXSIZE))

# ----- OFA -----
if (!exists("RUN_OFA")) RUN_OFA <- TRUE
if (!exists("OFA_MIN_CELLS_FOCAL")) OFA_MIN_CELLS_FOCAL <- 50
if (!exists("OFA_MIN_CELLS_REST")) OFA_MIN_CELLS_REST <- 100
if (!exists("OFA_PADJ_THR")) OFA_PADJ_THR <- 0.05
if (!exists("OFA_LFC_THR")) OFA_LFC_THR <- 0.25
if (!exists("OFA_TOP_N")) OFA_TOP_N <- 200
if (!exists("OFA_MAX_CELLS_PER_IDENT")) OFA_MAX_CELLS_PER_IDENT <- 5000

# ----- ssGSEA for CHOIR clusters -----
if (!exists("SSGSEA_CHOIR_METHODS")) SSGSEA_CHOIR_METHODS <- c("hallmark", "go_bp")
if (!exists("SSGSEA_CHOIR_N_TOP")) SSGSEA_CHOIR_N_TOP <- 20
if (!exists("SSGSEA_CHOIR_N_HEATMAP")) SSGSEA_CHOIR_N_HEATMAP <- 30
if (!exists("RUN_CLUSTER_SSGSEA")) RUN_CLUSTER_SSGSEA <- RUN_SSGSEA

# ----- L3 OFA -----
if (!exists("RUN_L3_OFA")) RUN_L3_OFA <- RUN_OFA
if (!exists("RUN_L3_OFA_VS_REST")) RUN_L3_OFA_VS_REST <- RUN_L3_OFA
if (!exists("RUN_L3_OFA_SAME_L2")) RUN_L3_OFA_SAME_L2 <- RUN_L3_OFA
if (!exists("RUN_L3_OFA_INTER_TISSUE")) RUN_L3_OFA_INTER_TISSUE <- RUN_L3_OFA
if (!exists("L3_OFA_MIN_CELLS_FOCAL")) L3_OFA_MIN_CELLS_FOCAL <- OFA_MIN_CELLS_FOCAL
if (!exists("L3_OFA_MIN_CELLS_REST")) L3_OFA_MIN_CELLS_REST <- OFA_MIN_CELLS_REST
if (!exists("L3_OFA_PADJ_THR")) L3_OFA_PADJ_THR <- OFA_PADJ_THR
if (!exists("L3_OFA_LFC_THR")) L3_OFA_LFC_THR <- OFA_LFC_THR
if (!exists("L3_OFA_TOP_N")) L3_OFA_TOP_N <- OFA_TOP_N
if (!exists("L3_OFA_MAX_CELLS_PER_IDENT")) L3_OFA_MAX_CELLS_PER_IDENT <- OFA_MAX_CELLS_PER_IDENT
if (!exists("EXIT_AFTER_L3_OFA")) EXIT_AFTER_L3_OFA <- FALSE

# ----- MiloPy tissue DA -----
if (!exists("RUN_MILOPY")) RUN_MILOPY <- TRUE
if (!exists("MILOPY_SCRIPT_PATH")) MILOPY_SCRIPT_PATH <- "/home/h2048/script/py/tissue_subcluster_milopy_20260418_v1.py"
if (!exists("MILOPY_CONDA_EXE")) MILOPY_CONDA_EXE <- "/home/h2048/miniconda3/bin/conda"
if (!exists("MILOPY_ENV_PREFIX")) MILOPY_ENV_PREFIX <- "/home/h2048/miniconda3/envs/scarches_stable_pertpy"
if (!exists("MILOPY_LEVEL_COLS")) MILOPY_LEVEL_COLS <- c(CELLTYPE_L2_COL, CELLTYPE_L3_COL)
if (!exists("MILOPY_LATENT_KEY_CANDIDATES")) MILOPY_LATENT_KEY_CANDIDATES <- c(
  "X_scanvi", "X_scANVI", "X_scanvi_corrected", "X_scvi", "X_harmony", "X_pca"
)
if (!exists("MILOPY_UMAP_KEY_CANDIDATES")) MILOPY_UMAP_KEY_CANDIDATES <- c(
  "X_umap_scanvi", "X_umap_scANVI", "X_umap_scanvi_corrected", "X_umap_scvi", "X_umap_harmony", "X_umap"
)
if (!exists("MILOPY_MIN_CELLS_PER_CELLTYPE")) MILOPY_MIN_CELLS_PER_CELLTYPE <- 200L
if (!exists("MILOPY_MIN_CELLS_PER_SAMPLE")) MILOPY_MIN_CELLS_PER_SAMPLE <- 10L
if (!exists("MILOPY_MIN_SAMPLES_PER_TISSUE")) MILOPY_MIN_SAMPLES_PER_TISSUE <- 3L
if (!exists("MILOPY_N_NEIGHBORS")) MILOPY_N_NEIGHBORS <- 30L
if (!exists("MILOPY_NHOOD_PROP")) MILOPY_NHOOD_PROP <- 0.1
if (!exists("MILOPY_ALPHA")) MILOPY_ALPHA <- 0.1
if (!exists("MILOPY_RANDOM_SEED")) MILOPY_RANDOM_SEED <- 0L
if (!exists("MILOPY_MAX_CELLTYPES_PER_LEVEL")) MILOPY_MAX_CELLTYPES_PER_LEVEL <- NULL
if (!exists("MILOPY_MAKE_PLOTS")) MILOPY_MAKE_PLOTS <- TRUE
if (!exists("MILOPY_WRITE_MILO_H5AD")) MILOPY_WRITE_MILO_H5AD <- FALSE
if (!exists("MILOPY_LOG_FILENAME")) MILOPY_LOG_FILENAME <- "milopy_run.log"
if (!exists("MILOPY_STOP_ON_ERROR")) MILOPY_STOP_ON_ERROR <- FALSE
if (!exists("MILOPY_OUTPUT_DIRNAME")) MILOPY_OUTPUT_DIRNAME <- "pertpy_milo"
if (!exists("FINAL_SCHEMA_VERSION")) FINAL_SCHEMA_VERSION <- "tc_final_v1"
if (!exists("FINAL_SCHEMA_OBJECT_ROLE")) FINAL_SCHEMA_OBJECT_ROLE <- "final"
if (!exists("FINAL_SCHEMA_LATENT_REDUCTION_CANDIDATES")) FINAL_SCHEMA_LATENT_REDUCTION_CANDIDATES <- NULL
if (!exists("FINAL_SCHEMA_UMAP_REDUCTION_CANDIDATES")) FINAL_SCHEMA_UMAP_REDUCTION_CANDIDATES <- NULL

if (!exists("KNOWN_MARKERS")) {
  KNOWN_MARKERS <- c(
    "CD79A", "CD79B", "MS4A1", "CD19", "PAX5",
    "IGHD", "IGHM", "TCL1A", "FCER2", "IL4R",
    "CD27", "TNFRSF13B", "AIM2",
    "BCL6", "AICDA", "RGS13", "MEF2B", "MME", "MKI67",
    "PRDM1", "XBP1", "MZB1", "JCHAIN", "SDC1", "DERL3", "SSR4",
    "IGHA1", "IGHG1", "IGHE",
    "ITGAX", "TBX21", "FCRL4", "FCRL5",
    "CD69", "CD86", "ISG15",
    "LAPTM5", "CD74"
  )
}

if (!exists("MARKER_PANEL_GROUP_COL")) MARKER_PANEL_GROUP_COL <- CELLTYPE_L3_COL
if (!exists("MARKER_PANEL_FIG_SUBDIR")) MARKER_PANEL_FIG_SUBDIR <- "marker_panels"
if (!exists("MARKER_PANEL_HEATMAP_CELLS_PER_GROUP")) MARKER_PANEL_HEATMAP_CELLS_PER_GROUP <- 80L

if (!exists("CUSTOM_MARKERS_DB")) {
  CUSTOM_MARKERS_DB <- data.frame(
    subtype = c(
      "Naive_B", "Transitional_B",
      "Memory_B_Unswitched", "Memory_B_Switched",
      "Atypical_Memory_B",
      "GC_B_Dark_Zone", "GC_B_Light_Zone",
      "Plasmablast", "Plasma_IgA", "Plasma_IgG", "Plasma_IgE",
      "Breg"
    ),
    markers = c(
      "IGHD,IGHM,TCL1A,FCER2,IL4R,CD38,CD24",
      "CD24,CD38,IGHM,IGHD,MME,SOX4",
      "CD27,IGHM,IGHD,TNFRSF13B",
      "CD27,IGHA1,IGHG1,IGHG2,AIM2,TNFRSF13B",
      "ITGAX,TBX21,FCRL4,FCRL5,ZEB2,CXCR3,FGR",
      "BCL6,AICDA,MKI67,TOP2A,CXCR4,FOXO1",
      "BCL6,LMO2,RGS13,MEF2B,CD83,CXCR5",
      "PRDM1,XBP1,MZB1,JCHAIN,MKI67,IRF4",
      "PRDM1,XBP1,MZB1,JCHAIN,SDC1,IGHA1,IGHA2,DERL3,SSR4",
      "PRDM1,XBP1,MZB1,JCHAIN,SDC1,IGHG1,IGHG2,IGHG3,DERL3,SSR4",
      "PRDM1,XBP1,MZB1,JCHAIN,SDC1,IGHE,DERL3,SSR4",
      "IL10,CD24,CD27,GZMB,TGFB1"
    ),
    stringsAsFactors = FALSE
  )
}

if (!exists("MARKER_PANELS")) {
  lineage_tag_upper <- toupper(trimws(as.character(LINEAGE_TAG)))
  MARKER_PANELS <- if (exists("tc_default_bcell_marker_panels", mode = "function") &&
                        identical(lineage_tag_upper, "BCELL")) {
    tc_default_bcell_marker_panels(
      known_markers = KNOWN_MARKERS,
      custom_markers_db = CUSTOM_MARKERS_DB
    )
  } else if (exists("tc_default_tnk_marker_panels", mode = "function") &&
             identical(lineage_tag_upper, "TNK")) {
    tc_default_tnk_marker_panels(
      known_markers = KNOWN_MARKERS,
      custom_markers_db = CUSTOM_MARKERS_DB
    )
  } else if (exists("tc_default_epithelial_marker_panels", mode = "function") &&
             identical(lineage_tag_upper, "EPITHELIAL")) {
    tc_default_epithelial_marker_panels(
      known_markers = KNOWN_MARKERS,
      custom_markers_db = CUSTOM_MARKERS_DB
    )
  } else {
    list()
  }
}

# ----- Tissue-Pair-Specific Context for LLM -----
if (!exists("LINEAGE_BASE_CONTEXT")) {
  LINEAGE_BASE_CONTEXT <- paste(
    "B cells and plasma cells from NORMAL (non-diseased) human respiratory tract tissues.",
    "This is a cross-site anatomical comparison, NOT a disease vs healthy comparison.",
    "Key B cell subtypes: Naive B, Memory B (unswitched/switched/atypical),",
    "GC B (dark zone centroblast / light zone centrocyte / transitional),",
    "Plasma cell (IgA/IgG).",
    "Focus: regional variation in B cell composition, differentiation state,",
    "and mucosal humoral immunity along the respiratory tract."
  )
}

if (!exists("TISSUE_CONTEXT")) {
  TISSUE_CONTEXT <- list(
    "nose" = paste(
      "Nasal cavity: first-line mucosal barrier, high antigen exposure,",
      "NALT (nasal-associated lymphoid tissue) supports local B cell responses,",
      "IgA-dominant secretory immunity, resident memory B cells."
    ),
    "sinus" = paste(
      "Paranasal sinus: semi-enclosed mucosal cavity, lower antigen load than nose,",
      "dependent on drainage for antigen clearance,",
      "potential site for ectopic lymphoid structures in inflammation,",
      "IgA secretion for mucosal defense."
    ),
    "respiratory airway" = paste(
      "Conducting airways (trachea/bronchi): ciliated epithelium with mucociliary clearance,",
      "BALT (bronchus-associated lymphoid tissue) in some individuals,",
      "IgA/IgG transport across epithelium, tissue-resident memory B cells,",
      "interface between upper and lower respiratory immunity."
    ),
    "lung parenchyma" = paste(
      "Lung parenchyma (alveolar region): gas exchange surface,",
      "thin epithelial barrier, alveolar macrophage-dominated immunity,",
      "lower B cell density than conducting airways,",
      "IgG-dominant (vs IgA in upper airways), systemic-like immune features,",
      "interstitial B cells and plasma cells near bronchovascular bundles."
    )
  )
}

if (!exists("PIPELINE_CHANGELOG_TITLE")) PIPELINE_CHANGELOG_TITLE <- sprintf("%s Changes:", PIPELINE_VERSION_LABEL)
if (!exists("PIPELINE_CHANGELOG_LINES")) {
  PIPELINE_CHANGELOG_LINES <- c(
    "  [LLM-6] Pairwise DE LLM integrates up/down evidence together across multiple databases",
    "  [RESUME-1] Prefer previous final object and previous-run summary when available",
    "  [ENV-1] Auto-load DEEPSEEK_API_KEY from common .env files",
    "  [PB-1] One-shot pseudobulk cache by sample x tissue x group; no repeated subset() in DE loops",
    "  [LLM-3] Validated gene->pathway map derived from enrichment geneID hits and passed to LLM evidence",
    "  [LLM-4] ssGSEA per-group LLM: each tissue x L2/L3 group gets an independent interpretation",
    "  [LLM-5] CHOIR per-cluster LLM: each CHOIR cluster gets one integrated up/down interpretation",
    "  [WX-1] Exploratory wilcox outputs labeled with analysis_tier/statistical_unit/replicate_model",
    "  [P2-1] stratified_downsample() ignores NA groups explicitly",
    "  [P2-2] h5ad sparse export builds CSR directly from dgCMatrix slots",
    "  [P2-3] NormalizeData provenance recorded in obj@misc$normalization_for_report"
  )
}
if (!exists("PIPELINE_INHERITED_FIXES_TITLE")) PIPELINE_INHERITED_FIXES_TITLE <- "v2.2.1 Fixes (inherited):"
if (!exists("PIPELINE_INHERITED_FIX_LINES")) {
  PIPELINE_INHERITED_FIX_LINES <- c(
    "  [P1-1] UMAP output without raster / alpha; png dpi 300",
    "  [P1-2] aggregate_pseudobulk sparse-friendly (no dense pre-alloc)",
    "  [P1-3] Unmapped L3 -> stop() fail-fast",
    sprintf("  [P1-4] CellMarker/PanglaoDB strict %s term filter", LINEAGE_DISPLAY),
    "  [P1-5] interpret_agent retry + atomic-result-tolerant normalization",
    "  [P1-6] Secondary DeepSeek standardization to fixed JSON schema",
    "  [P2-6] Heatmap title: 'visualization only'",
    "  [P2-7] set.seed(42) in stratified_downsample"
  )
}

# ==============================================================================
# [LLM-1] build_tissue_pair_context
# ==============================================================================
build_tissue_pair_context <- function(t1, t2, celltype_label, dir_name,
                                      celltype_level = "L2",
                                      gene_pathway_map_text = "") {
  ctx1 <- TISSUE_CONTEXT[[t1]]
  ctx2 <- TISSUE_CONTEXT[[t2]]
  if (is.null(ctx1)) ctx1 <- paste(t1, ": no specific context available.")
  if (is.null(ctx2)) ctx2 <- paste(t2, ": no specific context available.")

  validated_map_block <- if (nzchar(safe_trim(gene_pathway_map_text))) {
    paste(
      "\n\nValidated gene-pathway hits (use ONLY these for explicit gene->pathway claims):",
      gene_pathway_map_text,
      sep = "\n"
    )
  } else {
    "\n\nValidated gene-pathway hits: none available from current evidence. Do not invent gene->pathway links."
  }

  direction_text <- if (identical(dir_name, "integrated_up_down")) {
    sprintf("integrated up/down programs distinguishing %s from %s", t2, t1)
  } else {
    sprintf("genes %s-regulated in %s relative to %s", dir_name, t2, t1)
  }

  sprintf(
    paste(
      "%s",
      "\nComparing %s vs %s for %s-level group '%s'.",
      "Direction: %s.",
      "\n--- %s ---\n%s",
      "\n--- %s ---\n%s",
      "\nBiological question: What regional differences in %s explain the observed",
      "transcriptional divergence between these two anatomical sites?",
      "\n\nOutput requirements:",
      "1. Write ALL interpretive content in Chinese (Chinese characters).",
      "2. Only state a gene->pathway relationship when it is explicitly supported by the",
      "   validated gene-pathway hits provided below.",
      "3. If no validated map is available for a highlighted gene, describe the gene and",
      "   the enriched programs separately instead of forcing a direct mapping.",
      "4. When a validated mapping is available, use the format:",
      "   gene(log2FC=X.X) -> pathway name.",
      "%s"
    ),
    LINEAGE_BASE_CONTEXT, t2, t1, celltype_level, celltype_label,
    direction_text,
    t2, ctx2, t1, ctx1,
    BIOLOGICAL_QUESTION_FRAGMENT,
    validated_map_block
  )
}

build_integrated_tissue_pair_context <- function(t1, t2, celltype_label,
                                                 celltype_level = "L2",
                                                 gene_pathway_map_text = "") {
  ctx1 <- TISSUE_CONTEXT[[t1]]
  ctx2 <- TISSUE_CONTEXT[[t2]]
  if (is.null(ctx1)) ctx1 <- paste(t1, ": no specific context available.")
  if (is.null(ctx2)) ctx2 <- paste(t2, ": no specific context available.")

  validated_map_block <- if (nzchar(safe_trim(gene_pathway_map_text))) {
    paste(
      "\n\nValidated gene-pathway hits (use ONLY these for explicit gene->pathway claims):",
      gene_pathway_map_text,
      sep = "\n"
    )
  } else {
    "\n\nValidated gene-pathway hits: none available from current evidence. Do not invent gene->pathway links."
  }

  sprintf(
    paste(
      "%s",
      "\nComparing %s vs %s for %s-level group '%s'.",
      "This is an integrated bidirectional comparison: jointly consider features higher in %s relative to %s",
      "and features higher in %s relative to %s.",
      "\n--- %s ---\n%s",
      "\n--- %s ---\n%s",
      "\nBiological question: When evidence from both directions is integrated together, what regional differences in %s best explain the observed transcriptional divergence between these two anatomical sites?",
      "\n\nOutput requirements:",
      "1. Write ALL interpretive content in Chinese (Chinese characters).",
      "2. Only state a gene->pathway relationship when it is explicitly supported by the",
      "   validated gene-pathway hits provided below.",
      "3. If no validated map is available for a highlighted gene, describe the gene and",
      "   the enriched programs separately instead of forcing a direct mapping.",
      "4. Integrate evidence from both directions into ONE unified interpretation; do not",
      "   write separate disconnected up/down conclusions.",
      "5. If a validated mapping is available, use the format:",
      "   gene(log2FC=X.X) -> pathway name.",
      "6. Start overview with a concise cell-type/state judgment, then summarize the",
      "   comparative biological state in 1-2 sentences.",
      "%s"
    ),
    LINEAGE_BASE_CONTEXT,
    t2, t1, celltype_level, celltype_label,
    t2, t1, t1, t2,
    t2, ctx2, t1, ctx1,
    BIOLOGICAL_QUESTION_FRAGMENT,
    validated_map_block
  )
}

# ==============================================================================
# 2. Load Libraries
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("%s Tissue Comparison Pipeline %s (PRODUCTION)\n",
            LINEAGE_DISPLAY, PIPELINE_VERSION_LABEL))
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

suppressPackageStartupMessages({
  library(reticulate)
  library(SCNT)
  library(Seurat)
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(fanyi)
  library(GSVA)
  library(BiocParallel)
  library(CHOIR)
})

if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
  py_bin <- Sys.getenv("RETICULATE_PYTHON")
} else {
  use_condaenv(PYTHON_CONDA_ENV, required = TRUE)
  py_bin <- tryCatch(py_config()$python, error = function(e) paste0("conda:", PYTHON_CONDA_ENV))
}

LOCAL_GETSEURAT <- "/home/h2048/script/R/GetSeurat.R"
if (file.exists(LOCAL_GETSEURAT)) source(LOCAL_GETSEURAT)

if (ENABLE_LLM) {
  fanyi::set_translate_option(key = DEEPSEEK_API_KEY, source = "deepseek")
  cat("[OK] LLM enabled (DeepSeek)\n")
} else {
  cat("[INFO] LLM disabled (no API key)\n")
}
cat(sprintf("[OK] Python configured: %s\n", py_bin))
if (file.exists(LOCAL_GETSEURAT)) cat(sprintf("[OK] Local GetSeurat loaded: %s\n", LOCAL_GETSEURAT))
cat("[OK] Libraries loaded\n\n")

tc_should_filter_ig_genes <- function() {
  isTRUE(FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT)
}

tc_should_filter_technical_genes <- function() {
  isTRUE(FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT)
}

tc_should_filter_noncoding_locus_genes <- function() {
  isTRUE(FILTER_NONCODING_LOCUS_GENES_FOR_LLM_AND_ENRICHMENT)
}

tc_should_filter_cilia_genes <- function() {
  isTRUE(FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT)
}

tc_is_ig_related_gene <- function(genes) {
  genes <- toupper(trimws(as.character(genes)))
  genes[is.na(genes)] <- ""
  nzchar(genes) & grepl(IG_RELATED_GENE_REGEX, genes, perl = TRUE)
}

tc_is_cilia_related_gene <- function(genes) {
  genes <- toupper(trimws(as.character(genes)))
  genes[is.na(genes)] <- ""
  nzchar(genes) & grepl(CILIA_RELATED_GENE_REGEX, genes, perl = TRUE)
}

tc_is_technical_gene <- function(genes) {
  genes <- toupper(trimws(as.character(genes)))
  genes[is.na(genes)] <- ""
  nzchar(genes) & grepl(TECHNICAL_GENE_REGEX, genes, perl = TRUE)
}

tc_is_noncoding_locus_gene <- function(genes) {
  genes <- toupper(trimws(as.character(genes)))
  genes[is.na(genes)] <- ""
  nzchar(genes) & grepl(NONCODING_LOCUS_GENE_REGEX, genes, perl = TRUE)
}

tc_is_excluded_gene <- function(genes) {
  exclude <- rep(FALSE, length(genes))
  if (tc_should_filter_technical_genes()) exclude <- exclude | tc_is_technical_gene(genes)
  if (tc_should_filter_noncoding_locus_genes()) exclude <- exclude | tc_is_noncoding_locus_gene(genes)
  if (tc_should_filter_ig_genes()) exclude <- exclude | tc_is_ig_related_gene(genes)
  if (tc_should_filter_cilia_genes()) exclude <- exclude | tc_is_cilia_related_gene(genes)
  exclude
}

tc_filter_gene_symbols <- function(genes, unique_only = TRUE) {
  genes <- toupper(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes <- genes[!tc_is_excluded_gene(genes)]
  if (unique_only) unique(genes) else genes
}

tc_filter_named_gene_fc <- function(gene_fc) {
  if (is.null(gene_fc) || length(gene_fc) == 0) return(numeric())
  keep <- !is.na(names(gene_fc)) & nzchar(trimws(names(gene_fc)))
  gene_fc <- gene_fc[keep]
  if (length(gene_fc) == 0) return(numeric())
  gene_fc <- gene_fc[!tc_is_excluded_gene(names(gene_fc))]
  stats::setNames(as.numeric(gene_fc), toupper(names(gene_fc)))
}

tc_filter_df_by_gene_col <- function(df, gene_col) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col) || !gene_col %in% colnames(df)) {
    return(df)
  }
  keep <- !is.na(df[[gene_col]]) & nzchar(trimws(as.character(df[[gene_col]])))
  keep <- keep & !tc_is_excluded_gene(df[[gene_col]])
  df[keep, , drop = FALSE]
}

tc_filter_term2gene_df <- function(t2g_df) {
  if (is.null(t2g_df) || !is.data.frame(t2g_df) || nrow(t2g_df) == 0 || !"gene" %in% colnames(t2g_df)) {
    return(t2g_df)
  }
  genes_norm <- toupper(trimws(as.character(t2g_df$gene)))
  keep <- !is.na(genes_norm) & nzchar(genes_norm)
  keep <- keep & !tc_is_excluded_gene(genes_norm)
  t2g_df <- t2g_df[keep, , drop = FALSE]
  if (nrow(t2g_df) == 0) return(t2g_df)
  t2g_df$gene <- toupper(trimws(as.character(t2g_df$gene)))
  unique(t2g_df)
}

tc_filter_enrichment_gene_ids <- function(x) {
  genes <- tc_filter_gene_symbols(unlist(strsplit(as.character(x), "/", fixed = TRUE)), unique_only = TRUE)
  paste(genes, collapse = "/")
}

tc_sanitize_enrichment_df <- function(df, pathway_col = "Description", geneid_col = "geneID") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)
  if (!is.null(geneid_col) && geneid_col %in% colnames(df)) {
    df[[geneid_col]] <- vapply(df[[geneid_col]], tc_filter_enrichment_gene_ids, character(1))
    df <- df[is.na(df[[geneid_col]]) | nzchar(trimws(df[[geneid_col]])), , drop = FALSE]
    if ("Count" %in% colnames(df) && nrow(df) > 0) {
      gene_counts <- vapply(df[[geneid_col]], function(x) {
        length(tc_filter_gene_symbols(unlist(strsplit(as.character(x), "/", fixed = TRUE)), unique_only = TRUE))
      }, integer(1))
      df$Count <- ifelse(gene_counts > 0L, gene_counts, df$Count)
    }
  }
  if (!is.null(pathway_col) && pathway_col %in% colnames(df)) {
    df[[pathway_col]] <- trimws(as.character(df[[pathway_col]]))
  }
  df
}

tc_sanitize_enrichment_obj <- function(enrich_obj) {
  if (is.null(enrich_obj)) return(NULL)
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(enrich_obj)
  enrich_df <- tc_sanitize_enrichment_df(enrich_df)
  has_result_slot <- tryCatch("result" %in% methods::slotNames(enrich_obj), error = function(e) FALSE)
  if (inherits(enrich_obj, "enrichResult") && isTRUE(has_result_slot)) {
    enrich_obj@result <- enrich_df
    return(enrich_obj)
  }
  enrich_df
}

tc_is_energy_pathway <- function(pathways) {
  if (!isTRUE(DEPRIORITIZE_ENERGY_PATHWAYS)) return(rep(FALSE, length(pathways)))
  pathways <- tolower(trimws(as.character(pathways)))
  pathways[is.na(pathways)] <- ""
  nzchar(pathways) & grepl(ENERGY_PATHWAY_REGEX, pathways, ignore.case = TRUE, perl = TRUE)
}

tc_is_technical_pathway <- function(pathways) {
  if (!isTRUE(DEPRIORITIZE_TECHNICAL_PATHWAYS)) return(rep(FALSE, length(pathways)))
  pathways <- tolower(trimws(as.character(pathways)))
  pathways[is.na(pathways)] <- ""
  nzchar(pathways) & grepl(TECHNICAL_PATHWAY_REGEX, pathways, ignore.case = TRUE, perl = TRUE)
}

tc_is_cilia_pathway <- function(pathways) {
  if (!isTRUE(DEPRIORITIZE_CILIA_PATHWAYS)) return(rep(FALSE, length(pathways)))
  pathways <- tolower(trimws(as.character(pathways)))
  pathways[is.na(pathways)] <- ""
  nzchar(pathways) & grepl(CILIA_PATHWAY_REGEX, pathways, ignore.case = TRUE, perl = TRUE)
}

tc_is_excluded_top_pathway <- function(pathways) {
  tc_is_energy_pathway(pathways) | tc_is_technical_pathway(pathways) | tc_is_cilia_pathway(pathways)
}

tc_select_top_pathway_rows <- function(df, n_top, pathway_col = "Description") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)
  n_top <- max(1L, as.integer(n_top))
  if (!pathway_col %in% colnames(df)) {
    return(utils::head(df, n_top))
  }
  exclude_flag <- tc_is_excluded_top_pathway(df[[pathway_col]])
  keep_idx <- integer()
  exclude_idx <- integer()
  for (idx in seq_len(nrow(df))) {
    if (isTRUE(exclude_flag[idx])) {
      exclude_idx <- c(exclude_idx, idx)
    } else {
      keep_idx <- c(keep_idx, idx)
    }
    if (length(keep_idx) >= n_top) break
  }
  if (length(keep_idx) == 0) {
    if (isTRUE(APPEND_ENERGY_PATHWAYS_AFTER_TOP)) return(utils::head(df, n_top))
    return(df[0, , drop = FALSE])
  }
  selected_idx <- if (isTRUE(APPEND_ENERGY_PATHWAYS_AFTER_TOP)) c(keep_idx, exclude_idx) else keep_idx
  df[selected_idx, , drop = FALSE]
}

# ==============================================================================
# 3. Load GMT + Reference Databases
# ==============================================================================

cat("=== Loading Reference Databases ===\n")

gmt_all <- tryCatch({
  g <- read.gmt(MSIGDB_GMT_PATH); g$gene <- toupper(g$gene)
  cat(sprintf("[OK] MSigDB: %d gene sets\n", length(unique(g$term)))); g
}, error = function(e) { cat("[WARN] MSigDB failed\n"); NULL })

hallmark_t2g <- if (!is.null(gmt_all)) {
  h <- gmt_all %>% filter(grepl("^HALLMARK_", term))
  if (nrow(h) > 0) { cat(sprintf("[OK] Hallmark: %d\n", length(unique(h$term)))); h } else NULL
} else NULL

kegg_t2g <- if (!is.null(gmt_all)) {
  k <- gmt_all %>% filter(grepl("^KEGG_", term))
  if (nrow(k) > 0) { cat(sprintf("[OK] KEGG: %d\n", length(unique(k$term)))); k } else NULL
} else NULL

go_bp_t2g <- NULL; go_mf_t2g <- NULL; go_cc_t2g <- NULL
tryCatch({
  gmt_go <- read.gmt(GMT_GO_ALL); gmt_go$gene <- toupper(gmt_go$gene)
  go_bp_t2g <- gmt_go %>% filter(grepl("^GOBP_", term))
  go_mf_t2g <- gmt_go %>% filter(grepl("^GOMF_", term))
  go_cc_t2g <- gmt_go %>% filter(grepl("^GOCC_", term))
  cat(sprintf("[OK] GO: BP=%d, MF=%d, CC=%d\n",
              length(unique(go_bp_t2g$term)),
              length(unique(go_mf_t2g$term)),
              length(unique(go_cc_t2g$term))))
}, error = function(e) cat("[WARN] GO GMT failed\n"))

cellmarker_t2g <- NULL
tryCatch({
  cm_db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
  cm_db <- cm_db[grepl("Human", species, ignore.case = TRUE)]
  cat(sprintf("[OK] CellMarker raw: %d human entries\n", nrow(cm_db)))
  cm_bcell <- cm_db[grepl(CELLMARKER_CELLTYPE_PATTERN, cell_name, ignore.case = TRUE)]
  cat(sprintf("[OK] CellMarker lineage-specific: %d entries\n", nrow(cm_bcell)))
  if (nrow(cm_bcell) < 10) {
    cat(sprintf("[WARN] Too few %s entries in CellMarker; skipping\n", LINEAGE_DISPLAY))
  } else {
    t2g_list <- list()
    for (i in seq_len(nrow(cm_bcell))) {
      markers_raw <- cm_bcell$marker[i]; cell_type <- cm_bcell$cell_name[i]
      if (is.na(markers_raw) || markers_raw == "") next
      markers <- unlist(strsplit(markers_raw, "[,;\\s]+"))
      markers <- gsub('["\r\n\\[\\]]', '', markers)
      markers <- trimws(toupper(markers))
      markers <- unique(markers[markers != "" & !is.na(markers)])
      if (length(markers) > 0)
        t2g_list[[length(t2g_list) + 1]] <- data.frame(term = cell_type, gene = markers, stringsAsFactors = FALSE)
    }
    cellmarker_t2g <- tc_filter_term2gene_df(dplyr::bind_rows(t2g_list) %>% distinct())
    cat(sprintf("[OK] CellMarker TERM2GENE: %d pairs (%d terms)\n",
                nrow(cellmarker_t2g), length(unique(cellmarker_t2g$term))))
  }
}, error = function(e) cat(sprintf("[WARN] CellMarker failed: %s\n", e$message)))

panglaodb_t2g <- NULL
tryCatch({
  pdb <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)
  setnames(pdb, old = c("official gene symbol", "cell type"),
           new = c("gene_symbol", "cell_type"), skip_absent = TRUE)
  pdb <- pdb[grepl("Hs", species, fixed = TRUE)]
  cat(sprintf("[OK] PanglaoDB raw: %d human markers\n", nrow(pdb)))
  pdb_bcell <- pdb[grepl(PANGLAODB_CELLTYPE_PATTERN, cell_type, ignore.case = TRUE)]
  cat(sprintf("[OK] PanglaoDB lineage-specific: %d markers (%d types)\n",
              nrow(pdb_bcell), length(unique(pdb_bcell$cell_type))))
  if (nrow(pdb_bcell) < 5) {
    cat(sprintf("[WARN] Too few %s markers in PanglaoDB; skipping\n", LINEAGE_DISPLAY))
  } else {
    panglaodb_t2g <- pdb_bcell %>%
      dplyr::select(cell_type, gene_symbol) %>%
      mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
      filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
      distinct() %>%
      dplyr::rename(term = cell_type, gene = gene_symbol)
    panglaodb_t2g <- tc_filter_term2gene_df(panglaodb_t2g)
    cat(sprintf("[OK] PanglaoDB TERM2GENE: %d pairs (%d terms)\n",
                nrow(panglaodb_t2g), length(unique(panglaodb_t2g$term))))
  }
}, error = function(e) cat(sprintf("[WARN] PanglaoDB failed: %s\n", e$message)))

custom_t2g <- CUSTOM_MARKERS_DB %>%
  tidyr::separate_rows(markers, sep = ",") %>%
  dplyr::mutate(markers = trimws(toupper(markers))) %>%
  dplyr::filter(markers != "") %>%
  dplyr::select(term = subtype, gene = markers) %>%
  distinct()
custom_t2g <- tc_filter_term2gene_df(custom_t2g)
cat(sprintf("[OK] Custom %s: %d subtypes, %d pairs\n",
            CUSTOM_DB_LABEL, length(unique(custom_t2g$term)), nrow(custom_t2g)))
cat("\n")

# ==============================================================================
# 4. Helper Functions
# ==============================================================================

safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

level_file_stem <- function(base_stem, level_name, keep_legacy_l2 = TRUE) {
  if (keep_legacy_l2 && identical(level_name, "L2")) return(base_stem)
  sprintf("%s_%s", base_stem, tolower(level_name))
}

pick_reduction <- function(obj, preferred) {
  red <- Reductions(obj)
  hit <- preferred[preferred %in% red]
  if (length(hit) > 0) hit[[1]] else NULL
}

build_umap_plot <- function(obj, reduction_name, group_col, title,
                            split_col = NULL, label = FALSE, width = 10, height = 8,
                            cols = NULL) {
  p <- DimPlot(
    obj, reduction = reduction_name, group.by = group_col, split.by = split_col,
    pt.size = UMAP_PT_SIZE, shuffle = TRUE, label = label, repel = label, cols = cols
  ) + ggtitle(title) + coord_equal() + theme_classic(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
  list(plot = p, width = width, height = height)
}

# ----- 4.1 GMT enrichment -----
get_enrichment_size_rule <- function(db_name) {
  rule <- ENRICHMENT_GS_SIZE_RULES[[db_name]]
  if (is.null(rule)) rule <- ENRICHMENT_GS_SIZE_RULES[["default"]]
  list(min = as.integer(rule[["min"]]), max = as.integer(rule[["max"]]))
}

run_gmt_enrichment <- function(gene_list, t2g, db_name, tested_genes = NULL) {
  if (is.null(t2g) || nrow(t2g) == 0 || length(gene_list) < 5) return(NULL)
  t2g_use <- t2g %>%
    mutate(term = as.character(term), gene = toupper(as.character(gene))) %>%
    filter(!is.na(term), !is.na(gene), term != "", gene != "") %>%
    distinct(term, gene)
  t2g_use <- tc_filter_term2gene_df(t2g_use)
  if (nrow(t2g_use) == 0) return(NULL)
  universe_use <- if (!is.null(tested_genes)) {
    unique(intersect(tc_filter_gene_symbols(tested_genes), unique(t2g_use$gene)))
  } else unique(t2g_use$gene)
  if (length(universe_use) < 5) {
    cat(sprintf("    [INFO] %s: too few universe genes (%d)\n", db_name, length(universe_use)))
    return(NULL)
  }
  t2g_use <- t2g_use %>% filter(gene %in% universe_use)
  size_rule <- get_enrichment_size_rule(db_name)
  gs_sizes  <- t2g_use %>% count(term, name = "gs_size")
  max_gs    <- suppressWarnings(max(gs_sizes$gs_size, na.rm = TRUE))
  if (!is.finite(max_gs) || max_gs < 2) return(NULL)
  min_gs    <- max(2L, min(size_rule$min, as.integer(max_gs)))
  max_gs_ok <- max(min_gs, size_rule$max)
  valid_terms <- gs_sizes %>% filter(gs_size >= min_gs, gs_size <= max_gs_ok) %>% pull(term)
  if (length(valid_terms) == 0) return(NULL)
  t2g_use   <- t2g_use %>% filter(term %in% valid_terms)
  genes_use <- intersect(tc_filter_gene_symbols(gene_list), unique(t2g_use$gene))
  if (length(genes_use) < 5) return(NULL)
  tryCatch(
    suppressMessages(enricher(
      gene = genes_use, TERM2GENE = t2g_use, universe = unique(t2g_use$gene),
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH",
      minGSSize = min_gs, maxGSSize = max_gs_ok
    )),
    error = function(e) { cat(sprintf("    [WARN] %s: %s\n", db_name, e$message)); NULL }
  )
}

# ----- 4.2 Pseudobulk aggregation -----
build_pseudobulk_cache <- function(obj, group_col) {
  meta <- obj@meta.data
  required <- c(SAMPLE_COL, TISSUE_COL, group_col)
  valid <- complete.cases(meta[, required, drop = FALSE])
  valid <- valid & trimws(as.character(meta[[SAMPLE_COL]])) != ""
  valid <- valid & trimws(as.character(meta[[TISSUE_COL]])) != ""
  valid <- valid & trimws(as.character(meta[[group_col]])) != ""
  if (!any(valid)) return(NULL)
  meta_sub <- meta[valid, , drop = FALSE]
  if (exists("tc_resolve_pseudobulk_group_columns", mode = "function")) {
    resolved <- tc_resolve_pseudobulk_group_columns(
      meta_sub,
      sample_col = SAMPLE_COL,
      tissue_col = TISSUE_COL,
      candidate_cols = PSEUDOBULK_DISAMBIGUATION_CANDIDATES
    )
  } else {
    resolved <- list(
      group_cols = c(SAMPLE_COL, TISSUE_COL),
      disambiguation_cols = character(),
      ambiguous_keys_remaining = character()
    )
  }
  if (length(resolved$ambiguous_keys_remaining) > 0) {
    stop(sprintf(
      "Non-unique sample+tissue keys remain after disambiguation for group column '%s'. Ambiguous: %s",
      group_col,
      paste(utils::head(resolved$ambiguous_keys_remaining, 10), collapse = ", ")
    ))
  }
  group_cols <- unique(c(resolved$group_cols, group_col))
  meta_sub <- meta_sub[, group_cols, drop = FALSE]
  meta_sub$cell <- rownames(meta[valid, , drop = FALSE])
  meta_sub$group_value <- as.character(meta_sub[[group_col]])
  meta_sub$pb_group <- apply(meta_sub[, c(resolved$group_cols, "group_value"), drop = FALSE], 1, function(v) {
    if (exists("tc_normalize_meta_values", mode = "function")) {
      paste(tc_normalize_meta_values(v), collapse = "__")
    } else {
      paste(trimws(as.character(v)), collapse = "__")
    }
  })
  pb_meta <- meta_sub %>%
    group_by(pb_group, sample = .data[[SAMPLE_COL]], tissue = .data[[TISSUE_COL]], group_value) %>%
    summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    filter(n_cells >= MIN_CELLS_PER_PSEUDOBULK) %>%
    arrange(group_value, tissue, sample)
  if (nrow(pb_meta) == 0) return(NULL)
  pb_meta$pb_group <- as.character(pb_meta$pb_group)
  if (anyDuplicated(pb_meta$pb_group)) {
    dup_groups <- unique(pb_meta$pb_group[duplicated(pb_meta$pb_group)])
    cat(sprintf(
      "[WARN] build_pseudobulk_cache(): deduplicating %d repeated pseudobulk keys for '%s'\n",
      length(dup_groups), group_col
    ))
    pb_meta <- pb_meta[!duplicated(pb_meta$pb_group), , drop = FALSE]
  }
  pb_levels <- pb_meta$pb_group
  dup_level_idx <- which(duplicated(pb_levels) | duplicated(pb_levels, fromLast = TRUE))
  if (length(dup_level_idx) > 0) {
    cat(sprintf(
      "[DEBUG] build_pseudobulk_cache(): duplicated factor levels detected for '%s' after dedup step (%d rows)\n",
      group_col, length(dup_level_idx)
    ))
    debug_cols <- intersect(c("pb_group", "sample", "tissue", "group_value"), colnames(pb_meta))
    print(utils::head(pb_meta[dup_level_idx, debug_cols, drop = FALSE], 20))
  }
  keep_cells  <- meta_sub$pb_group %in% pb_meta$pb_group
  meta_keep   <- meta_sub[keep_cells, , drop = FALSE]
  group_factor <- factor(as.character(meta_keep$pb_group), levels = pb_levels)
  counts_mat  <- GetAssayData(obj, layer = "counts")[, meta_keep$cell, drop = FALSE]
  cell_to_pb  <- Matrix::sparseMatrix(
    i = seq_len(nrow(meta_keep)), j = as.integer(group_factor), x = 1,
    dims = c(nrow(meta_keep), nrow(pb_meta))
  )
  pb_counts <- as(counts_mat %*% cell_to_pb, "dgCMatrix")
  pb_meta   <- as.data.frame(pb_meta, stringsAsFactors = FALSE)
  colnames(pb_counts) <- pb_meta$pb_group
  rownames(pb_meta)   <- pb_meta$pb_group
  list(
    counts = pb_counts,
    meta = pb_meta,
    group_col = group_col,
    group_cols = group_cols,
    disambiguation_cols = resolved$disambiguation_cols
  )
}

aggregate_pseudobulk <- function(pseudobulk_cache, group_value, min_total_cells = 30L) {
  if (is.null(pseudobulk_cache)) return(NULL)
  pb_meta <- pseudobulk_cache$meta
  keep    <- pb_meta$group_value == group_value
  if (!any(keep)) return(NULL)
  if (sum(pb_meta$n_cells[keep], na.rm = TRUE) < min_total_cells) return(NULL)
  pb_counts   <- pseudobulk_cache$counts[, keep, drop = FALSE]
  pb_meta_sub <- pb_meta[keep, , drop = FALSE]
  valid <- colnames(pb_counts)[Matrix::colSums(pb_counts) > 0 & !is.na(pb_meta_sub[colnames(pb_counts), "tissue"])]
  if (length(valid) < 4) return(NULL)
  list(
    counts = pb_counts[, valid, drop = FALSE],
    meta = pb_meta_sub[valid, , drop = FALSE],
    group_cols = pseudobulk_cache$group_cols,
    disambiguation_cols = pseudobulk_cache$disambiguation_cols
  )
}

# ----- 4.3 DESeq2 pairwise -----
run_deseq2_pairwise <- function(pb, t1, t2) {
  keep <- pb$meta$tissue %in% c(t1, t2)
  if (sum(keep) < 4) return(NULL)
  counts_sub <- pb$counts[, keep, drop = FALSE]
  meta_sub   <- pb$meta[keep, , drop = FALSE]
  tissue_levels_safe <- make.names(c(t1, t2), unique = TRUE)
  meta_sub$tissue      <- droplevels(factor(meta_sub$tissue, levels = c(t1, t2)))
  meta_sub$tissue_safe <- factor(
    ifelse(meta_sub$tissue == t1, tissue_levels_safe[1], tissue_levels_safe[2]),
    levels = tissue_levels_safe
  )
  n1 <- sum(meta_sub$tissue == t1); n2 <- sum(meta_sub$tissue == t2)
  if (n1 < MIN_SAMPLES_PER_TISSUE || n2 < MIN_SAMPLES_PER_TISSUE) return(NULL)
  keep_genes <- rowSums(counts_sub >= 1) >= max(3, ncol(counts_sub) * 0.2)
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  if (nrow(counts_sub) < 100) return(NULL)
  tested_genes <- rownames(counts_sub)
  counts_int   <- round(as.matrix(counts_sub))
  storage.mode(counts_int) <- "integer"
  design_formula <- ~ tissue_safe
  design_used    <- design_formula
  if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
    n_levels <- length(unique(meta_sub[[DESIGN_COVARIATE]]))
    if (n_levels >= 2 && n_levels < nrow(meta_sub)) {
      meta_sub[[DESIGN_COVARIATE]] <- factor(meta_sub[[DESIGN_COVARIATE]])
      design_formula <- as.formula(paste("~", DESIGN_COVARIATE, "+ tissue_safe"))
    }
  }
  dds <- tryCatch({
    design_used <- design_formula
    dds <- DESeqDataSetFromMatrix(counts_int, meta_sub, design_formula)
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    tryCatch({
      design_used <<- ~ tissue_safe
      dds2 <- DESeqDataSetFromMatrix(counts_int, meta_sub, ~ tissue_safe)
      DESeq(dds2, quiet = TRUE)
    }, error = function(e2) { cat(sprintf("    [ERROR] DESeq2: %s\n", e2$message)); NULL })
  })
  if (is.null(dds)) return(NULL)
  res <- results(dds, contrast = c("tissue_safe", tissue_levels_safe[2], tissue_levels_safe[1]),
                 alpha = PADJ_THR)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(sig = ifelse(padj < PADJ_THR & abs(log2FoldChange) > LFC_THR, "sig", "ns"),
           direction = ifelse(log2FoldChange > 0, "up", "down"))
  list(de_table = res_df, tested_genes = tested_genes,
       n_up = sum(res_df$sig == "sig" & res_df$direction == "up"),
       n_down = sum(res_df$sig == "sig" & res_df$direction == "down"),
       n_samples_1 = n1, n_samples_2 = n2, tissue_1 = t1, tissue_2 = t2,
      design = deparse(design_used))
}

# ----- 4.4 Exploratory wilcox -----
run_wilcox_exploratory <- function(obj, group_value, group_col = CELLTYPE_L2_COL) {
  cells <- colnames(obj)[obj@meta.data[[group_col]] == group_value]
  if (length(cells) < WILCOX_MIN_CELLS) return(NULL)
  sub     <- subset(obj, cells = cells)
  tissues <- sort(unique(na.omit(sub@meta.data[[TISSUE_COL]])))
  if (length(tissues) < 2) return(NULL)
  Idents(sub) <- TISSUE_COL
  results <- list()
  for (pair in combn(tissues, 2, simplify = FALSE)) {
    t1 <- pair[1]; t2 <- pair[2]
    n1 <- sum(sub@meta.data[[TISSUE_COL]] == t1, na.rm = TRUE)
    n2 <- sum(sub@meta.data[[TISSUE_COL]] == t2, na.rm = TRUE)
    if (n1 < WILCOX_MIN_CELLS || n2 < WILCOX_MIN_CELLS) next
    de <- tryCatch(
      FindMarkers(sub, ident.1 = t2, ident.2 = t1, test.use = "wilcox",
                  min.pct = 0.1, logfc.threshold = 0, max.cells.per.ident = 5000),
      error = function(e) NULL
    )
    if (!is.null(de) && nrow(de) > 0) {
      de$gene <- rownames(de)
      de$padj <- p.adjust(de$p_val, method = "BH")
      results[[paste0(t2, "_vs_", t1)]] <- de %>%
        filter(abs(avg_log2FC) >= WILCOX_LFC_THR) %>%
        arrange(padj) %>% head(WILCOX_TOP_N) %>%
        mutate(analysis_tier = "exploratory_cell_level", statistical_unit = "cell",
               replicate_model = "none",
               celltype_level = ifelse(identical(group_col, CELLTYPE_L3_COL), "L3", "L2"),
               celltype_group = group_value, comparison = paste0(t2, "_vs_", t1),
               tissue_case = t2, tissue_reference = t1,
               n_cells_case = n2, n_cells_reference = n1, grouping_column = group_col,
               .before = 1)
    }
  }
  results
}

# ----- 4.5 interpret_agent + standardize -----

sanitize_utf8_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- as.character(x)
  x <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = ""))
  x[is.na(x)] <- ""; x
}

safe_trim <- function(x) {
  x <- to_scalar(x); x <- sanitize_utf8_text(x); trimws(x)
}

capture_object_text <- function(x) {
  if (is.null(x)) return("")
  txt <- tryCatch(paste(capture.output(print(x)), collapse = "\n"), error = function(e) "")
  if (!nzchar(trimws(txt)))
    txt <- tryCatch(paste(capture.output(str(x, max.level = 3)), collapse = "\n"), error = function(e) "")
  safe_trim(txt)
}

unwrap_interpret_agent_result <- function(x) {
  if (inherits(x, "interpretation_list") && length(x) >= 1) return(x[[1]])
  if (is.list(x) && length(x) == 1 && (is.list(x[[1]]) || !is.null(names(x[[1]])))) return(x[[1]])
  x
}

has_named_entries <- function(x) !is.null(names(x)) && any(nzchar(names(x)))

extract_named_text <- function(x, candidates) {
  for (nm in candidates) {
    val <- NULL
    if (is.list(x) && !is.null(x[[nm]])) val <- x[[nm]]
    else if (!is.list(x) && has_named_entries(x) && nm %in% names(x)) val <- x[[nm]]
    if (!is.null(val)) { txt <- safe_trim(val); if (nzchar(txt)) return(txt) }
  }
  ""
}

is_retryable_interpret_agent_issue <- function(warnings = character(), error_message = NULL) {
  signals <- c(warnings, error_message)
  if (length(signals) == 0) return(FALSE)
  any(grepl("Failed to parse JSON response|premature EOF|invalid for atomic vectors|429|5[0-9]{2}|timeout|temporar",
            signals, ignore.case = TRUE, perl = TRUE))
}

looks_structured_interpret_agent_result <- function(x) {
  core <- unwrap_interpret_agent_result(x)
  any(nzchar(c(
    extract_named_text(core, c("overview", "summary", "interpretation", "narrative")),
    extract_named_text(core, c("key_mechanisms", "mechanisms", "keyMechanisms")),
    extract_named_text(core, c("hypothesis", "model", "working_hypothesis")),
    extract_named_text(core, c("key_drivers", "drivers", "genes", "gene_drivers"))
  )))
}

placeholder_text <- function(x, default = "Not available from current evidence.") {
  x <- safe_trim(x); if (nzchar(x)) x else default
}

collapse_driver_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- trimws(as.character(x)); x <- x[nzchar(x)]
  paste(unique(x), collapse = ", ")
}

truncate_text <- function(x, max_chars = STANDARDIZE_LLM_MAX_INPUT_CHARS) {
  x <- safe_trim(x); if (!nzchar(x)) return("")
  if (nchar(x, type = "chars") <= max_chars) return(x)
  paste0(substr(x, 1, max_chars), "\n...[truncated]")
}

build_enrichment_evidence_text <- function(enrich_obj, gene_fc = NULL, n_terms = 8L, n_genes = 15L) {
  if (is.null(enrich_obj)) return("")
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return("")
  enrich_df <- tc_sanitize_enrichment_df(enrich_df)
  enrich_df <- enrich_df[order(enrich_df$p.adjust), , drop = FALSE]
  top_terms <- tc_select_top_pathway_rows(enrich_df, n_terms)
  term_lines <- vapply(seq_len(nrow(top_terms)), function(i) {
    sprintf("- %s | padj=%s | Count=%s | Genes=%s",
            safe_trim(top_terms$Description[i]),
            format(top_terms$p.adjust[i], scientific = TRUE, digits = 3),
            safe_trim(top_terms$Count[i]),
            truncate_text(gsub("/", ", ", safe_trim(top_terms$geneID[i])), 200))
  }, character(1))
  gene_lines <- character()
  if (!is.null(gene_fc) && length(gene_fc) > 0) {
    gene_fc <- tc_filter_named_gene_fc(gene_fc)
    gene_fc <- sort(gene_fc[!is.na(gene_fc)], decreasing = TRUE)
    top_up   <- utils::head(gene_fc, n_genes)
    top_down <- utils::head(sort(gene_fc, decreasing = FALSE), n_genes)
    gene_lines <- c(
      sprintf("Top up genes: %s",   paste(sprintf("%s(%.2f)", names(top_up),   top_up),   collapse = ", ")),
      sprintf("Top down genes: %s", paste(sprintf("%s(%.2f)", names(top_down), top_down), collapse = ", "))
    )
  }
  gene_pathway_map  <- build_gene_pathway_map(enrich_obj, gene_fc = gene_fc, n_terms = n_terms)
  gene_pathway_lines <- format_gene_pathway_map_text(gene_pathway_map)
  paste(c("Top enrichment evidence:", term_lines,
          if (length(gene_pathway_lines) > 0) "Validated gene-pathway hits:" else NULL,
          gene_pathway_lines, gene_lines), collapse = "\n")
}

split_gene_ids <- function(x) {
  x <- safe_trim(x); if (!nzchar(x)) return(character())
  tc_filter_gene_symbols(unlist(strsplit(x, "/", fixed = TRUE)), unique_only = TRUE)
}

build_gene_pathway_map <- function(enrich_obj, gene_fc = NULL, n_terms = 8L,
                                   n_top_genes = 10L, max_pathways_per_gene = 3L) {
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
  if (is.null(enrich_df) || nrow(enrich_df) == 0 || is.null(gene_fc) || length(gene_fc) == 0)
    return(data.frame())
  enrich_df <- tc_sanitize_enrichment_df(enrich_df)
  gene_fc   <- tc_filter_named_gene_fc(gene_fc)
  gene_fc   <- gene_fc[!is.na(gene_fc)]
  if (length(gene_fc) == 0) return(data.frame())
  enrich_df <- enrich_df[order(enrich_df$p.adjust), , drop = FALSE]
  top_terms <- tc_select_top_pathway_rows(enrich_df, n_terms)
  top_genes <- names(sort(abs(gene_fc), decreasing = TRUE))[seq_len(min(n_top_genes, length(gene_fc)))]
  top_genes <- toupper(top_genes)
  map_rows  <- lapply(seq_len(nrow(top_terms)), function(i) {
    genes_i <- intersect(split_gene_ids(top_terms$geneID[i]), top_genes)
    if (length(genes_i) == 0) return(NULL)
    data.frame(gene = genes_i, log2FC = unname(gene_fc[genes_i]),
               pathway = safe_trim(top_terms$Description[i]),
               pathway_padj = top_terms$p.adjust[i], stringsAsFactors = FALSE)
  })
  map_df <- dplyr::bind_rows(map_rows)
  if (nrow(map_df) == 0) return(data.frame())
  map_df %>%
    mutate(abs_log2FC = abs(log2FC)) %>%
    arrange(desc(abs_log2FC), pathway_padj, pathway) %>%
    group_by(gene) %>% slice_head(n = max_pathways_per_gene) %>% ungroup() %>%
    select(gene, log2FC, pathway, pathway_padj)
}

format_gene_pathway_map_text <- function(gene_pathway_map) {
  if (is.null(gene_pathway_map) || nrow(gene_pathway_map) == 0) return(character())
  apply(gene_pathway_map, 1, function(row) {
    sprintf("- %s(log2FC=%.2f) -> %s | pathway_padj=%s",
            row[["gene"]], as.numeric(row[["log2FC"]]), row[["pathway"]],
            format(as.numeric(row[["pathway_padj"]]), scientific = TRUE, digits = 3))
  })
}

extract_json_string <- function(text) {
  text      <- safe_trim(text); if (!nzchar(text)) return("")
  text      <- gsub("^```(?:json)?\\s*", "", text, perl = TRUE)
  text      <- gsub("\\s*```$",           "", text, perl = TRUE)
  start_idx <- regexpr("\\{", text, perl = TRUE)[1]
  end_positions <- gregexpr("\\}", text, perl = TRUE)[[1]]
  if (start_idx < 1 || length(end_positions) == 0 || end_positions[1] < 1) return("")
  end_idx <- end_positions[length(end_positions)]
  substr(text, start_idx, end_idx)
}

parse_standardized_json <- function(text) {
  candidates <- unique(c(safe_trim(text), extract_json_string(text)))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) == 0) return(NULL)
  for (candidate in candidates) {
    parsed <- tryCatch(jsonlite::fromJSON(candidate, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(parsed) && is.list(parsed)) return(parsed)
  }
  NULL
}

# [LLM-2] standardize_result_with_llm
standardize_result_with_llm <- function(raw_text, context_str, evidence_text,
                                        gene_pathway_map_text,
                                        celltype_label, celltype_level,
                                        comp_name, dir_name, source_db,
                                        warnings = character(), error_message = NULL) {
  if (!ENABLE_LLM || !STANDARDIZE_LLM_OUTPUT)
    return(list(result = NULL, warnings = character(), error = "standardization disabled"))

  noninformative_gene_rule_text <- if (exists("tc_noninformative_gene_rule_text", mode = "function")) {
    tc_noninformative_gene_rule_text()
  } else {
    paste(
      "For non-B-cell analyses, do not use IG genes as primary evidence; for all analyses,",
      "treat MT-, ribosomal (RPS/RPL/MRPS/MRPL), ENSG/LINC, and placeholder non-coding loci",
      "(for example AC123456.1, AL123456.1, AP000000.1, CTB-/CTC-/CTD- clone-style loci, LOC genes, and *-OT transcripts)",
      "as non-informative unless no stronger lineage-relevant evidence exists."
    )
  }

  prompt <- paste(
    sprintf("You are standardizing a biological interpretation for single-cell %s tissue comparison.",
            LINEAGE_CONTEXT_LOWER),
    "Return valid JSON only. No markdown, no code fences, no extra commentary.",
    "Use exactly these keys:",
    "overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "",
    "Rules:",
    "1. overview/key_mechanisms/hypothesis/narrative/evidence/limitations must be strings.",
    "2. key_drivers must be an array of short gene or regulator names.",
    "3. Every field must be present.",
    "4. If evidence is weak, keep the content conservative and say so in limitations.",
    "5. If a field cannot be recovered, use 'Not available from current evidence.'.",
    "6. Write ALL narrative content (overview, key_mechanisms, hypothesis, narrative,",
    "   evidence, limitations) in Chinese (Simplified Chinese characters).",
    "   key_drivers should remain as gene symbol strings (e.g. [\"BCL6\", \"AICDA\"]).",
    "7. Gene->pathway links must be conservative: only use the validated gene-pathway",
    "   map provided below. If a highlighted gene lacks a validated mapping, describe",
    "   the gene and enriched programs separately instead of inventing a direct link.",
    "8. In key_mechanisms, prioritize genes with the largest absolute log2FoldChange",
    "   values. When a validated mapping exists, use the format:",
    "   gene(log2FC=X.X) -> pathway: explanation.",
    "9. overview must start with a concise cell-type/state judgment.",
    "10. If both up and down evidence are provided, integrate them into one judgment",
    "    instead of discussing them independently.",
    "11. When ssGSEA evidence is provided, use it mainly for cell-type/state judgment,",
    "    and avoid over-interpreting broad pathways as direct mechanisms.",
    paste0("12. ", noninformative_gene_rule_text),
    "13. Keep the output concise: overview <= 2 sentences, key_mechanisms <= 3 short clauses,",
    "    hypothesis <= 1 sentence, narrative <= 3 sentences, and do not split the final",
    "    answer into separate mini-conclusions by database or by up/down direction.",
    "14. hypothesis should be a brief comparative scientific hypothesis only when the",
    "    comparison is explicit; otherwise use 'Not applicable.'.",
    "15. In evidence, explicitly cite top up/down DEG and main positive/negative ssGSEA",
    "    signals when they are present in the evidence block.",
    if (length(LLM_EXTRA_RULES) > 0) paste("16. Additional rules:", paste(LLM_EXTRA_RULES, collapse = " | ")) else NULL,
    "",
    sprintf("Cell type level: %s", celltype_level),
    sprintf("Cell type group: %s", celltype_label),
    sprintf("Comparison: %s", comp_name),
    sprintf("Direction: %s", dir_name),
    sprintf("Source DB: %s", source_db),
    "Context:", truncate_text(context_str),
    "Original interpret_agent output:", truncate_text(raw_text),
    "Validated gene-pathway map:", truncate_text(gene_pathway_map_text),
    "Enrichment evidence:", truncate_text(evidence_text),
    "Warnings and errors:", truncate_text(paste(c(warnings, error_message), collapse = "\n")),
    sep = "\n"
  )

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(STANDARDIZE_LLM_MAX_RETRIES)) {
    if (attempt > 1) cat(sprintf("    [INFO] standardize_llm retry %d/%d\n", attempt, STANDARDIZE_LLM_MAX_RETRIES))
    else cat(sprintf("    [INFO] standardize_llm (%s)\n", STANDARDIZE_LLM_MODEL))
    response_text <- tryCatch(
      safe_trim(fanyi::chat_request(prompt, model = STANDARDIZE_LLM_MODEL, api_key = DEEPSEEK_API_KEY)),
      error = function(e) { last_error <<- conditionMessage(e); "" }
    )
    parsed <- parse_standardized_json(response_text)
    if (!is.null(parsed))
      return(list(result = parsed, warnings = unique(collected_warnings), error = NULL))
    if (nzchar(response_text))
      collected_warnings <- c(collected_warnings, sprintf("standardize attempt %d returned non-JSON text", attempt))
    if (attempt < STANDARDIZE_LLM_MAX_RETRIES) Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
  }
  list(result = NULL, warnings = unique(collected_warnings),
       error = ifelse(is.null(last_error) || !nzchar(last_error),
                      "failed to standardize interpret_agent output", last_error))
}

normalize_interpret_agent_result <- function(raw_result, celltype_label, celltype_level,
                                             comp_name, dir_name,
                                             celltype_l2 = NA_character_,
                                             source_db = NA_character_,
                                             warnings = character(),
                                             error_message = NULL,
                                             enrich_obj = NULL,
                                             context_str = NULL,
                                             gene_fc = NULL,
                                             gene_pathway_map = NULL) {
  core     <- unwrap_interpret_agent_result(raw_result)
  raw_text <- capture_object_text(raw_result)
  overview       <- extract_named_text(core, c("overview", "summary", "interpretation", "narrative"))
  key_mechanisms <- extract_named_text(core, c("key_mechanisms", "mechanisms", "keyMechanisms"))
  hypothesis     <- extract_named_text(core, c("hypothesis", "model", "working_hypothesis"))
  narrative      <- extract_named_text(core, c("narrative", "story", "details"))
  key_drivers    <- extract_named_text(core, c("key_drivers", "drivers", "genes", "gene_drivers"))
  evidence       <- extract_named_text(core, c("evidence", "supporting_evidence", "rationale"))
  limitations    <- extract_named_text(core, c("limitations", "caveats", "notes"))
  if (!nzchar(overview) && nzchar(raw_text)) overview <- raw_text
  major_fields_present <- c(overview, key_mechanisms, hypothesis, narrative, evidence)
  needs_standardization <- isTRUE(STANDARDIZE_LLM_OUTPUT) && (
    !all(nzchar(major_fields_present)) || !nzchar(key_drivers) ||
    !is.null(error_message) || !is.list(core)
  )
  if (needs_standardization) {
    evidence_text <- build_enrichment_evidence_text(enrich_obj, gene_fc = gene_fc)
    gene_pathway_map_text <- paste(format_gene_pathway_map_text(gene_pathway_map), collapse = "\n")
    std_payload <- standardize_result_with_llm(
      raw_text = raw_text, context_str = context_str,
      evidence_text = evidence_text, gene_pathway_map_text = gene_pathway_map_text,
      celltype_label = celltype_label, celltype_level = celltype_level,
      comp_name = comp_name, dir_name = dir_name, source_db = source_db,
      warnings = warnings, error_message = error_message
    )
    if (!is.null(std_payload$result)) {
      std            <- std_payload$result
      overview       <- placeholder_text(std$overview)
      key_mechanisms <- placeholder_text(std$key_mechanisms)
      hypothesis     <- placeholder_text(std$hypothesis)
      narrative      <- placeholder_text(std$narrative)
      key_drivers    <- placeholder_text(collapse_driver_field(std$key_drivers))
      evidence       <- placeholder_text(std$evidence)
      limitations    <- placeholder_text(std$limitations)
      warnings <- unique(c(warnings, std_payload$warnings,
                           if (!is.null(error_message) && nzchar(error_message))
                             sprintf("initial interpret_agent error: %s", error_message)))
      error_message <- NULL
    } else {
      warnings <- unique(c(warnings, std_payload$warnings))
      limitations <- placeholder_text(limitations)
      if (!is.null(std_payload$error) && nzchar(std_payload$error))
        error_message <- paste(unique(c(error_message, std_payload$error)), collapse = " | ")
    }
  }
  overview       <- placeholder_text(overview)
  key_mechanisms <- placeholder_text(key_mechanisms)
  hypothesis     <- placeholder_text(hypothesis)
  narrative      <- placeholder_text(narrative)
  key_drivers    <- placeholder_text(key_drivers)
  evidence       <- placeholder_text(evidence)
  limitations    <- placeholder_text(limitations)
  status <- dplyr::case_when(
    !is.null(error_message) ~ "error",
    is.null(raw_result)     ~ "missing",
    any(nzchar(c(overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations))) ~ "structured",
    nzchar(raw_text) ~ "raw_text",
    TRUE ~ "empty"
  )
  list(celltype_level = celltype_level, celltype_label = celltype_label,
       celltype_l2 = if (is.null(celltype_l2)) NA_character_ else celltype_l2,
       comparison = comp_name, direction = dir_name, source_db = source_db, status = status,
       warnings = unique(warnings), error = if (is.null(error_message)) "" else error_message,
       overview = overview, key_mechanisms = key_mechanisms, hypothesis = hypothesis,
       narrative = narrative, key_drivers = key_drivers, evidence = evidence, limitations = limitations,
       raw_text = raw_text, raw_result = raw_result)
}

resolve_group_l2_label <- function(obj, group_col, group_value) {
  if (identical(group_col, CELLTYPE_L2_COL)) return(as.character(group_value))
  hits <- unique(as.character(obj@meta.data[obj@meta.data[[group_col]] == group_value, CELLTYPE_L2_COL]))
  hits <- hits[!is.na(hits) & nzchar(trimws(hits))]
  if (length(hits) == 1) return(hits)
  if (length(hits) > 1) {
    cat(sprintf("[WARN] Non-unique L2 mapping for %s '%s'; structured celltype_l2 set to NA\n",
                group_col, group_value))
  }
  NA_character_
}

resolve_dominant_label <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(trimws(values))]
  if (length(values) == 0) return(NA_character_)
  freq <- sort(table(values), decreasing = TRUE)
  names(freq)[1]
}

llm_pick_existing_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

llm_num <- function(x) suppressWarnings(as.numeric(x))

append_nonempty_blocks <- function(..., sep = "\n\n") {
  blocks <- unlist(list(...), use.names = FALSE)
  blocks <- blocks[!is.na(blocks) & nzchar(trimws(blocks))]
  paste(blocks, collapse = sep)
}

format_top_deg_entries <- function(df, gene_col, logfc_col = NULL, padj_col = NULL,
                                   top_n = LLM_TOP_DEG_N, decreasing = TRUE,
                                   pct_lookup = NULL,
                                   pct1_col = NULL,
                                   pct2_col = NULL,
                                   pct1_label = "group_1",
                                   pct2_label = "group_2") {
  if (exists("tc_format_top_deg_entries", mode = "function")) {
    return(tc_format_top_deg_entries(
      df = df,
      gene_col = gene_col,
      logfc_col = logfc_col,
      padj_col = padj_col,
      top_n = top_n,
      decreasing = decreasing,
      pct_lookup = pct_lookup,
      pct1_col = pct1_col,
      pct2_col = pct2_col,
      pct1_label = pct1_label,
      pct2_label = pct2_label
    ))
  }
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col)) return("None")
  lfc <- if (!is.null(logfc_col)) llm_num(df[[logfc_col]]) else rep(NA_real_, nrow(df))
  padj <- if (!is.null(padj_col)) llm_num(df[[padj_col]]) else rep(Inf, nrow(df))
  ord <- if (!all(is.na(lfc))) {
    order(if (decreasing) -lfc else lfc, padj, na.last = TRUE)
  } else {
    order(padj, na.last = TRUE)
  }
  df <- df[ord, , drop = FALSE]
  df <- utils::head(df, max(1L, as.integer(top_n)))
  genes <- as.character(df[[gene_col]])
  genes[is.na(genes) | !nzchar(trimws(genes))] <- "NA"
  if (!is.null(logfc_col) && !is.null(padj_col)) {
    return(paste0(
      genes, "(", round(llm_num(df[[logfc_col]]), 2), ", padj=",
      signif(llm_num(df[[padj_col]]), 2), ")",
      collapse = ", "
    ))
  }
  if (!is.null(logfc_col)) {
    return(paste0(genes, "(", round(llm_num(df[[logfc_col]]), 2), ")", collapse = ", "))
  }
  paste(genes, collapse = ", ")
}

format_integrated_deg_entries <- function(df, gene_col, logfc_col, padj_col = NULL,
                                          top_n = LLM_TOP_DEG_N,
                                          pct_lookup = NULL,
                                          pct1_col = NULL,
                                          pct2_col = NULL,
                                          pct1_label = "group_1",
                                          pct2_label = "group_2") {
  if (exists("tc_format_integrated_deg_entries", mode = "function")) {
    return(tc_format_integrated_deg_entries(
      df = df,
      gene_col = gene_col,
      logfc_col = logfc_col,
      padj_col = padj_col,
      top_n = top_n,
      pct_lookup = pct_lookup,
      pct1_col = pct1_col,
      pct2_col = pct2_col,
      pct1_label = pct1_label,
      pct2_label = pct2_label
    ))
  }
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col) || is.null(logfc_col)) return("None")
  lfc <- llm_num(df[[logfc_col]])
  padj <- if (!is.null(padj_col)) llm_num(df[[padj_col]]) else rep(Inf, nrow(df))
  ord <- order(-abs(lfc), padj, na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- utils::head(df, max(1L, as.integer(top_n)))
  genes <- as.character(df[[gene_col]])
  genes[is.na(genes) | !nzchar(trimws(genes))] <- "NA"
  paste0(
    genes,
    "(log2FC=", round(llm_num(df[[logfc_col]]), 2),
    if (!is.null(padj_col)) paste0(", padj=", signif(llm_num(df[[padj_col]]), 2)) else "",
    ")",
    collapse = ", "
  )
}

build_top_deg_context <- function(de_df,
                                  top_n = LLM_TOP_DEG_N,
                                  padj_thr = LLM_DEG_PADJ_THR,
                                  lfc_thr = LLM_DEG_LFC_THR,
                                  pct_lookup_builder = NULL,
                                  pct1_label = "group_1",
                                  pct2_label = "group_2") {
  if (!isTRUE(LLM_INCLUDE_TOP_DEG) || is.null(de_df) || !is.data.frame(de_df) || nrow(de_df) == 0) return("")
  gene_col <- llm_pick_existing_col(de_df, c("gene", "symbol", "feature", "features", "genes"))
  logfc_col <- llm_pick_existing_col(de_df, c("avg_log2FC", "log2FoldChange", "log2FC", "avg_logFC"))
  padj_col <- llm_pick_existing_col(de_df, c("padj", "p_val_adj", "FDR", "qvalue", "p_adj"))
  if (is.null(gene_col) || is.null(logfc_col)) return("")
  de_df <- tc_filter_df_by_gene_col(de_df, gene_col)
  if (nrow(de_df) == 0) return("")

  lfc <- llm_num(de_df[[logfc_col]])
  keep <- !is.na(lfc)
  if (!is.null(padj_col)) {
    padj <- llm_num(de_df[[padj_col]])
    keep <- keep & !is.na(padj) & padj <= padj_thr
  }
  up_df <- de_df[keep & lfc >= lfc_thr, , drop = FALSE]
  down_df <- de_df[keep & lfc <= -lfc_thr, , drop = FALSE]
  if (nrow(up_df) == 0 && nrow(down_df) == 0) {
    up_df <- de_df[!is.na(lfc) & lfc > 0, , drop = FALSE]
    down_df <- de_df[!is.na(lfc) & lfc < 0, , drop = FALSE]
  }
  integrated_df <- dplyr::bind_rows(up_df, down_df)
  if (nrow(integrated_df) == 0) integrated_df <- de_df[!is.na(lfc), , drop = FALSE]
  pct_lookup <- if (is.function(pct_lookup_builder)) {
    tryCatch(
      pct_lookup_builder(unique(as.character(integrated_df[[gene_col]]))),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  paste(
    sprintf(
      "Integrated top DEG summary (combine both directions; use directly for cell-type/state judgment; expression proportions shown as pct[%s] and pct[%s]):",
      pct1_label,
      pct2_label
    ),
    paste0("- Top DEG across both directions: ",
           format_integrated_deg_entries(
             integrated_df,
             gene_col,
             logfc_col,
             padj_col,
             top_n,
             pct_lookup = pct_lookup,
             pct1_label = pct1_label,
             pct2_label = pct2_label
           )),
    sep = "\n"
  )
}

build_choir_deg_evidence_text <- function(de_df,
                                          cluster_id = NULL,
                                          top_n_abs = CHOIR_LLM_TOP_DEG_N,
                                          top_n_up = 12L,
                                          top_n_down = 12L,
                                          padj_thr = OFA_PADJ_THR,
                                          lfc_thr = OFA_LFC_THR,
                                          pct1_label = "cluster",
                                          pct2_label = "rest") {
  if (!exists("tc_build_choir_deg_evidence_text", mode = "function")) {
    stop("tc_build_choir_deg_evidence_text() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_build_choir_deg_evidence_text(
    de_df = de_df,
    cluster_id = cluster_id,
    top_n_abs = top_n_abs,
    top_n_up = top_n_up,
    top_n_down = top_n_down,
    padj_thr = padj_thr,
    lfc_thr = lfc_thr,
    pct1_label = pct1_label,
    pct2_label = pct2_label,
    caller_env = parent.frame()
  )
}

build_integrated_directional_enrichment_text <- function(evidence_df, direction_map = NULL) {
  if (is.null(evidence_df) || !is.data.frame(evidence_df) || nrow(evidence_df) == 0) return("")
  if (!"direction" %in% colnames(evidence_df)) return(build_multi_db_context_text(evidence_df))
  df <- evidence_df
  if (!"db" %in% colnames(df)) df$db <- "NA"
  if (!"Description" %in% colnames(df)) df$Description <- rownames(df)
  if (!"geneID" %in% colnames(df)) df$geneID <- ""
  if (!"Count" %in% colnames(df)) df$Count <- NA_integer_
  if (!"p.adjust" %in% colnames(df)) df$p.adjust <- NA_real_
  df <- tc_sanitize_enrichment_df(df)
  df$direction <- tolower(trimws(as.character(df$direction)))
  df <- df[order(match(df$direction, c("up", "down")), df$db, df$p.adjust, na.last = TRUE), , drop = FALSE]
  lines <- vapply(seq_len(nrow(df)), function(i) {
    dir_key <- df$direction[i]
    dir_label <- if (!is.null(direction_map) && dir_key %in% names(direction_map)) direction_map[[dir_key]] else toupper(dir_key)
    sprintf(
      "- [%s][%s] %s | padj=%s | Count=%s | Genes=%s",
      dir_label,
      safe_trim(df$db[i]),
      safe_trim(gsub("^\\[[^]]+\\]\\s*", "", safe_trim(df$Description[i]))),
      format(df$p.adjust[i], scientific = TRUE, digits = 3),
      safe_trim(df$Count[i]),
      truncate_text(gsub("/", ", ", safe_trim(df$geneID[i])), 220)
    )
  }, character(1))
  paste(c(
    "Integrated cross-database enrichment summary (all comparison directions are provided together; synthesize them in one judgment):",
    lines
  ), collapse = "\n")
}

flatten_interpretation_records <- function(x) {
  rows <- list(); idx <- 1L
  for (celltype_label in names(x)) {
    for (comp_name in names(x[[celltype_label]])) {
      for (dir_name in names(x[[celltype_label]][[comp_name]])) {
        rec <- x[[celltype_label]][[comp_name]][[dir_name]]
        if (is.null(rec)) next
        rows[[idx]] <- data.frame(
          celltype_level = safe_trim(rec$celltype_level),
          celltype_label = safe_trim(rec$celltype_label),
          celltype_l2    = safe_trim(rec$celltype_l2),
          comparison     = safe_trim(rec$comparison),
          direction      = safe_trim(rec$direction),
          source_db      = safe_trim(rec$source_db),
          status         = safe_trim(rec$status),
          warnings       = paste(rec$warnings, collapse = " | "),
          error          = safe_trim(rec$error),
          overview       = safe_trim(rec$overview),
          key_mechanisms = safe_trim(rec$key_mechanisms),
          hypothesis     = safe_trim(rec$hypothesis),
          narrative      = safe_trim(rec$narrative),
          key_drivers    = safe_trim(rec$key_drivers),
          evidence       = safe_trim(rec$evidence),
          limitations    = safe_trim(rec$limitations),
          raw_text       = safe_trim(rec$raw_text),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  dplyr::bind_rows(rows)
}

write_interpretation_markdown <- function(records, path, title, include_raw = FALSE) {
  md <- c(title, "", sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
  if (nrow(records) == 0) {
    md <- c(md, "No interpret_agent records available.")
    writeLines(md, path); return(invisible(path))
  }
  for (i in seq_len(nrow(records))) {
    rec <- records[i, , drop = FALSE]
    md  <- c(md, sprintf("## %s | %s | %s | %s",
                         rec$celltype_level, rec$celltype_label, rec$comparison, rec$direction), "")
    md  <- c(md, sprintf("- **Status:** %s", rec$status),
             sprintf("- **Cell Type Level:** %s", rec$celltype_level),
             sprintf("- **Source DB:** %s", ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
    if ("confidence" %in% colnames(records) && nzchar(safe_trim(rec$confidence))) {
      md <- c(md, sprintf("- **Confidence:** %s", rec$confidence))
    }
    if (nzchar(rec$warnings)) md <- c(md, sprintf("- **Warnings:** %s", rec$warnings))
    if (nzchar(rec$error))    md <- c(md, sprintf("- **Error:** %s",    rec$error))
    md <- c(md, "")
    if ("cell_type_judgment" %in% colnames(records) && nzchar(safe_trim(rec$cell_type_judgment))) {
      md <- c(md, "### Cell Type Judgment", "", placeholder_text(rec$cell_type_judgment), "")
    }
    if (any(c("annotated_l3_correspondence", "annotation_match_degree", "annotation_count_summary",
              "annotated_l3_dominant", "annotated_l3_dominant_pct", "annotated_l3_secondary", "annotated_l3_secondary_pct") %in% colnames(records))) {
      annotation_lines <- character()
      if ("annotation_match_degree" %in% colnames(records) && nzchar(safe_trim(rec$annotation_match_degree))) {
        annotation_lines <- c(annotation_lines, sprintf("- **Match Degree:** %s", rec$annotation_match_degree))
      }
      if ("annotated_l3_dominant" %in% colnames(records) && nzchar(safe_trim(rec$annotated_l3_dominant))) {
        annotation_lines <- c(
          annotation_lines,
          sprintf(
            "- **Dominant annotated L3:** %s%s",
            rec$annotated_l3_dominant,
            if ("annotated_l3_dominant_pct" %in% colnames(records) && !is.na(suppressWarnings(as.numeric(rec$annotated_l3_dominant_pct)))) {
              sprintf(" (%s)", tc_pct_label(rec$annotated_l3_dominant_pct))
            } else ""
          )
        )
      }
      if ("annotated_l3_secondary" %in% colnames(records) && nzchar(safe_trim(rec$annotated_l3_secondary))) {
        annotation_lines <- c(
          annotation_lines,
          sprintf(
            "- **Secondary annotated L3:** %s%s",
            rec$annotated_l3_secondary,
            if ("annotated_l3_secondary_pct" %in% colnames(records) && !is.na(suppressWarnings(as.numeric(rec$annotated_l3_secondary_pct)))) {
              sprintf(" (%s)", tc_pct_label(rec$annotated_l3_secondary_pct))
            } else ""
          )
        )
      }
      if ("annotated_l3_correspondence" %in% colnames(records) && nzchar(safe_trim(rec$annotated_l3_correspondence))) {
        annotation_lines <- c(annotation_lines, "", placeholder_text(rec$annotated_l3_correspondence))
      }
      if ("annotation_count_summary" %in% colnames(records) && nzchar(safe_trim(rec$annotation_count_summary))) {
        annotation_lines <- c(annotation_lines, "", placeholder_text(rec$annotation_count_summary))
      }
      if (length(annotation_lines) > 0) {
        md <- c(md, "### Annotation Alignment", "", annotation_lines, "")
      }
    }
    if ("outlier_assessment" %in% colnames(records) && nzchar(safe_trim(rec$outlier_assessment))) {
      md <- c(md, "### Outlier Assessment", "", placeholder_text(rec$outlier_assessment), "")
    }
    if ("discovery_assessment" %in% colnames(records) && nzchar(safe_trim(rec$discovery_assessment))) {
      md <- c(md, "### Discovery Assessment", "", placeholder_text(rec$discovery_assessment), "")
    }
    if ("integrated_diagnostic_comment" %in% colnames(records) && nzchar(safe_trim(rec$integrated_diagnostic_comment))) {
      md <- c(md, "### Integrated Diagnostic Comment", "", placeholder_text(rec$integrated_diagnostic_comment), "")
    }
    for (nm in c("Overview", "Key Mechanisms", "Hypothesis", "Narrative", "Key Drivers", "Evidence", "Limitations")) {
      field_key <- gsub(" ", "_", tolower(nm))
      field_key <- gsub("key_mechanisms", "key_mechanisms", field_key)
      val <- switch(nm,
        "Overview"       = rec$overview,       "Key Mechanisms" = rec$key_mechanisms,
        "Hypothesis"     = rec$hypothesis,     "Narrative"      = rec$narrative,
        "Key Drivers"    = rec$key_drivers,    "Evidence"       = rec$evidence,
        "Limitations"    = rec$limitations
      )
      md <- c(md, sprintf("### %s", nm), "", placeholder_text(val), "")
    }
    if (include_raw && nzchar(rec$raw_text)) md <- c(md, "### Raw Output", "", rec$raw_text, "")
  }
  writeLines(md, path); invisible(path)
}

choose_best_enrichment_db <- function(enr_list, preferred = NULL) {
  if (is.null(preferred))
    preferred <- c("GO_BP", "Hallmark", "KEGG", "GO_MF", "CellMarker", CUSTOM_DB_NAME)
  for (db in preferred) {
    er_df <- tryCatch(as.data.frame(enr_list[[db]]), error = function(e) data.frame())
    er_df <- tc_sanitize_enrichment_df(er_df)
    er_df <- tc_select_top_pathway_rows(er_df, 3)
    if (!is.null(enr_list[[db]]) && nrow(er_df) >= 1)
      return(list(enrich = enr_list[[db]], db = db))
  }
  list(enrich = NULL, db = NA_character_)
}

collapse_source_db_label <- function(db_names, representative_db = NULL) {
  db_names <- unique(db_names[nzchar(db_names)])
  if (length(db_names) == 0) {
    if (!is.null(representative_db) && nzchar(representative_db)) return(representative_db)
    return(NA_character_)
  }
  label <- paste(db_names, collapse = " + ")
  if (!is.null(representative_db) && nzchar(representative_db))
    return(sprintf("multi_db[%s] | seed=%s", label, representative_db))
  sprintf("multi_db[%s]", label)
}

select_enrichment_databases_for_llm <- function(enr_list, preferred = NULL,
                                                max_dbs = INTERPRET_MULTI_DB_MAX_DBS,
                                                min_terms = INTERPRET_MULTI_DB_MIN_TERMS) {
  if (is.null(preferred)) preferred <- INTERPRET_MULTI_DB_PREFERRED
  available <- names(enr_list)[vapply(names(enr_list), function(db) {
    er <- enr_list[[db]]
    er_df <- tryCatch(as.data.frame(er), error = function(e) data.frame())
    er_df <- tc_sanitize_enrichment_df(er_df)
    er_df <- tc_select_top_pathway_rows(er_df, max(1L, as.integer(min_terms)))
    !is.null(er) && nrow(er_df) > 0
  }, logical(1))]
  if (length(available) == 0) return(character())
  ordered <- unique(c(intersect(preferred, available), setdiff(available, preferred)))
  selected <- ordered[vapply(ordered, function(db) {
    er_df <- tryCatch(as.data.frame(enr_list[[db]]), error = function(e) data.frame())
    er_df <- tc_sanitize_enrichment_df(er_df)
    nrow(tc_select_top_pathway_rows(er_df, max(1L, as.integer(min_terms)))) >= 1
  }, logical(1))]
  if (length(selected) == 0) selected <- ordered
  utils::head(selected, max(1L, as.integer(max_dbs)))
}

combine_enrichment_results_for_llm <- function(enr_list, db_names,
                                               n_terms_per_db = INTERPRET_MULTI_DB_TERMS_PER_DB) {
  if (length(db_names) == 0) return(data.frame())
  rows <- lapply(db_names, function(db) {
    er <- enr_list[[db]]
    df <- tryCatch(as.data.frame(er), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!"Description" %in% colnames(df)) df$Description <- rownames(df)
    if (!"geneID" %in% colnames(df))      df$geneID <- ""
    if (!"Count" %in% colnames(df))       df$Count <- NA_integer_
    if (!"p.adjust" %in% colnames(df))    df$p.adjust <- NA_real_
    df <- tc_sanitize_enrichment_df(df)
    df <- df[order(df$p.adjust, na.last = TRUE), , drop = FALSE]
    df <- tc_select_top_pathway_rows(df, max(1L, as.integer(n_terms_per_db)))
    df$db <- db
    df$Description <- sprintf("[%s] %s", db, safe_trim(df$Description))
    df
  })
  dplyr::bind_rows(rows)
}

build_multi_db_context_text <- function(enrich_df) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return("")
  db_blocks <- unlist(lapply(split(enrich_df, enrich_df$db), function(df) {
    db_label <- unique(df$db)[1]
    lines <- vapply(seq_len(nrow(df)), function(i) {
      sprintf("- %s | padj=%s | Count=%s | Genes=%s",
              safe_trim(gsub("^\\[[^]]+\\]\\s*", "", safe_trim(df$Description[i]))),
              format(df$p.adjust[i], scientific = TRUE, digits = 3),
              safe_trim(df$Count[i]),
              truncate_text(gsub("/", ", ", safe_trim(df$geneID[i])), 220))
    }, character(1))
    c(sprintf("[%s]", db_label), lines, "")
  }), use.names = FALSE)
  paste(c(
    "Cross-database enrichment summary (prioritize signals supported across multiple databases; do not overfit to a single database):",
    db_blocks
  ), collapse = "\n")
}

prepare_llm_enrichment_bundle <- function(enr_list, gene_fc = NULL, preferred = NULL,
                                          max_dbs = INTERPRET_MULTI_DB_MAX_DBS,
                                          n_terms_per_db = INTERPRET_MULTI_DB_TERMS_PER_DB) {
  selected_dbs <- select_enrichment_databases_for_llm(enr_list, preferred = preferred, max_dbs = max_dbs)
  representative <- choose_best_enrichment_db(enr_list, preferred = selected_dbs)
  combined_df <- combine_enrichment_results_for_llm(enr_list, selected_dbs, n_terms_per_db = n_terms_per_db)
  evidence_obj <- if (nrow(combined_df) > 0) combined_df else representative$enrich
  gene_pathway_map <- build_gene_pathway_map(
    evidence_obj,
    gene_fc = gene_fc,
    n_terms = if (is.data.frame(evidence_obj)) max(8L, nrow(evidence_obj)) else 8L
  )
  list(
    interpret_enrich = representative$enrich,
    interpret_db = representative$db,
    selected_dbs = selected_dbs,
    source_db = collapse_source_db_label(selected_dbs, representative$db),
    evidence_df = combined_df,
    evidence_obj = evidence_obj,
    evidence_text = build_multi_db_context_text(combined_df),
    gene_pathway_map = gene_pathway_map
  )
}

# ----- 4.6 Stratified downsampling -----
stratified_downsample <- function(obj, group_col, n_per = HEATMAP_CELLS_PER_TYPE) {
  set.seed(42)
  meta   <- obj@meta.data
  groups <- na.omit(unique(as.character(meta[[group_col]])))
  groups <- groups[nzchar(trimws(groups))]
  cells  <- unlist(lapply(groups, function(g) {
    gc <- rownames(meta)[as.character(meta[[group_col]]) == g]
    sample(gc, min(length(gc), n_per))
  }))
  subset(obj, cells = cells)
}

ensure_dir_exists <- function(path, label = NULL, max_attempts = 3L, sleep_sec = 0.2) {
  path <- trimws(as.character(path))
  if (!nzchar(path)) stop("Directory path is empty.")
  if (file.exists(path) && !dir.exists(path)) {
    stop(sprintf("Path exists but is not a directory: %s", path))
  }
  max_attempts <- max(1L, as.integer(max_attempts))
  if (!dir.exists(path)) {
    for (attempt in seq_len(max_attempts)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(path)) break
      if (attempt < max_attempts) Sys.sleep(sleep_sec)
    }
  }
  if (!dir.exists(path)) {
    label_suffix <- if (!is.null(label) && nzchar(trimws(as.character(label)))) {
      sprintf(" for %s", trimws(as.character(label)))
    } else ""
    stop(sprintf("Failed to create directory%s: %s", label_suffix, path))
  }
  invisible(path)
}

# ----- 4.7 Save pdf + png -----
save_plot <- function(p, path_no_ext, width = 10, height = 8) {
  ensure_dir_exists(dirname(path_no_ext), label = "plot output")
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
}

sanitize_obs_for_h5ad <- function(df) {
  out <- df
  for (col in colnames(out)) {
    vec <- out[[col]]
    if      (is.factor(vec))                         out[[col]] <- as.character(vec)
    else if (is.logical(vec))                        out[[col]] <- ifelse(is.na(vec), NA_integer_, as.integer(vec))
    else if (inherits(vec, c("POSIXct", "POSIXt"))) out[[col]] <- format(vec, tz = "UTC", usetz = TRUE)
    else if (inherits(vec, "Date"))                  out[[col]] <- as.character(vec)
    else if (is.list(vec))                           out[[col]] <- vapply(vec, function(x) {
      if (length(x) == 0 || all(is.na(x))) return("")
      paste(as.character(x), collapse = "; ")
    }, character(1))
    else if (is.character(vec))                      out[[col]] <- trimws(vec)
    else if (!(is.integer(vec) || is.numeric(vec)))  out[[col]] <- as.character(vec)
    if (is.character(out[[col]])) {
      out[[col]][is.na(out[[col]])] <- ""
      out[[col]] <- trimws(out[[col]])
    }
  }
  out
}

matrix_to_scipy_csr <- function(mat, scipy_sparse, np) {
  mat_csc <- as(mat, "dgCMatrix")
  scipy_sparse$csc_matrix(
    reticulate::tuple(
      np$array(as.numeric(mat_csc@x), dtype = np$float32),
      np$array(as.integer(mat_csc@i), dtype = np$int32),
      np$array(as.integer(mat_csc@p), dtype = np$int32)
    ),
    shape = reticulate::tuple(as.integer(nrow(mat_csc)), as.integer(ncol(mat_csc)))
  )$transpose()$tocsr()
}

# ----- 4.8 ssGSEA helpers -----
term2gene_to_list <- function(term2gene_df) split(toupper(term2gene_df$gene), term2gene_df$term)

map_gene_sets_to_features <- function(gene_sets, features) {
  feats_upper <- toupper(features)
  keep <- !duplicated(feats_upper)
  m    <- setNames(features[keep], feats_upper[keep])
  gs2  <- lapply(gene_sets, function(g) unique(na.omit(m[toupper(g)])))
  gs2[lengths(gs2) > 0]
}

filter_gs_size <- function(gene_sets, minSize = 10, maxSize = 500) {
  gene_sets[lengths(gene_sets) >= minSize & lengths(gene_sets) <= maxSize]
}

extract_directional_ssgsea_top_rows <- function(score_mat, z_mat, top_n,
                                                group_field = "group",
                                                level_name = NULL,
                                                method = NULL) {
  rows <- lapply(colnames(score_mat), function(group_id) {
    sc <- score_mat[, group_id]
    z  <- z_mat[, group_id]
    pos_ord <- order(z, sc, decreasing = TRUE, na.last = TRUE)
    neg_ord <- order(z, sc, decreasing = FALSE, na.last = TRUE)
    select_idx <- function(ord) {
      ord <- ord[!is.na(z[ord])]
      if (length(ord) == 0) return(integer())
      ord_df <- data.frame(pathway = rownames(score_mat)[ord], idx = ord, stringsAsFactors = FALSE)
      ord_df <- tc_select_top_pathway_rows(ord_df, max(1L, as.integer(top_n)), pathway_col = "pathway")
      ord_df$idx
    }
    pos_idx <- select_idx(pos_ord)
    neg_idx <- select_idx(neg_ord)
    make_df <- function(idx, direction_label) {
      if (length(idx) == 0) return(NULL)
      df <- data.frame(
        pathway = rownames(score_mat)[idx],
        score = sc[idx],
        z_score = z[idx],
        rank = seq_along(idx),
        direction = direction_label,
        stringsAsFactors = FALSE
      )
      if (!is.null(level_name)) df$level <- level_name
      if (!is.null(method)) df$method <- method
      df[[group_field]] <- group_id
      df
    }
    dplyr::bind_rows(
      make_df(pos_idx, "positive"),
      make_df(neg_idx, "negative")
    )
  })
  dplyr::bind_rows(rows)
}

# ----- 4.9 ssGSEA group LLM interpretation helpers -----

build_ssgsea_evidence_text <- function(top_df_group, n_top = 20L,
                                       n_top_per_method = INTERPRET_SSGSEA_TERMS_PER_DB) {
  if (is.null(top_df_group) || nrow(top_df_group) == 0) return("")
  if (!"direction" %in% colnames(top_df_group)) {
    top_df_group$direction <- ifelse(top_df_group$z_score < 0, "negative", "positive")
  }
  top_df_group$direction <- ifelse(
    tolower(as.character(top_df_group$direction)) %in% c("negative", "down", "neg", "suppressed"),
    "negative", "positive"
  )
  method_vals <- if ("method" %in% colnames(top_df_group)) unique(as.character(top_df_group$method)) else character()
  method_vals <- method_vals[!is.na(method_vals) & nzchar(trimws(method_vals))]
  format_ssgsea_lines <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    label <- if ("method" %in% colnames(df) && length(method_vals) > 1) {
      sprintf("[%s] %s", as.character(df$method), safe_trim(df$pathway))
    } else {
      safe_trim(df$pathway)
    }
    vapply(seq_len(nrow(df)), function(i) {
      sprintf("- %s | z_score=%.3f | score=%.3f | rank=%d",
              label[i], df$z_score[i], df$score[i], df$rank[i])
    }, character(1))
  }
  pos_df <- top_df_group[top_df_group$direction == "positive", , drop = FALSE]
  neg_df <- top_df_group[top_df_group$direction == "negative", , drop = FALSE]
  pos_df <- pos_df[order(pos_df$z_score, pos_df$score, decreasing = TRUE), , drop = FALSE]
  neg_df <- neg_df[order(neg_df$z_score, neg_df$score, decreasing = FALSE), , drop = FALSE]
  n_select <- if (length(method_vals) > 1) max(as.integer(n_top), as.integer(n_top_per_method)) else as.integer(n_top)
  pos_df <- tc_select_top_pathway_rows(pos_df, n_select, pathway_col = "pathway")
  neg_df <- tc_select_top_pathway_rows(neg_df, n_select, pathway_col = "pathway")
  pos_lines <- if (nrow(pos_df) > 0) c(
    "Integrated positive / identity-supporting pathways:",
    format_ssgsea_lines(pos_df)
  ) else NULL
  neg_lines <- if (nrow(neg_df) > 0) c(
    "Integrated negative / suppressed pathways:",
    format_ssgsea_lines(neg_df)
  ) else NULL
  paste(c(
    "Integrated ssGSEA pathway summary across all selected gene set databases (do not interpret database-by-database; integrate positive and negative signals together):",
    pos_lines,
    neg_lines
  ), collapse = "\n")
}

build_ssgsea_group_context <- function(group_id, methods, level_name) {
  parts        <- strsplit(group_id, "__", fixed = TRUE)[[1]]
  tissue_label <- if (length(parts) >= 1) parts[1] else group_id
  ct_label     <- if (length(parts) >= 2) paste(parts[-1], collapse = " / ") else group_id
  tissue_ctx   <- TISSUE_CONTEXT[[tissue_label]]
  if (is.null(tissue_ctx)) tissue_ctx <- paste(tissue_label, ": no specific context available.")
  methods <- unique(as.character(methods))
  methods <- methods[!is.na(methods) & nzchar(trimws(methods))]
  method_label <- if (length(methods) == 0) "ssGSEA" else paste(methods, collapse = ", ")
  sprintf(
    paste(
      "%s",
      "\nssGSEA pathway activity profile for %s-level group: [%s].",
      "\n  Tissue: %s | Cell type: %s | Gene set database: %s",
      "\n--- Tissue context ---\n%s",
      "\nBiological question: What does the pathway activity pattern in this specific",
      "tissue x cell type combination reveal about %s biology at this anatomical site?",
      "\n\nOutput requirements:",
      "1. Write ALL interpretive content in Chinese (Simplified Chinese characters).",
      "2. key_drivers must remain as English gene symbols (e.g. [\"BCL6\", \"AICDA\"]).",
      "3. Synthesize concordant evidence across multiple ssGSEA databases/methods instead of interpreting each database separately.",
      "4. Use ssGSEA mainly to judge the most likely cell type/state, then briefly summarize the biological state.",
      "5. Integrate positive and negative pathways together; do not write separate disconnected conclusions.",
      "6. Keep the interpretation concise.",
      "7. Do not invent gene->pathway links not supported by the evidence."
    ),
    LINEAGE_BASE_CONTEXT,
    level_name, group_id,
    tissue_label, ct_label, method_label,
    tissue_ctx,
    BIOLOGICAL_QUESTION_FRAGMENT
  )
}

run_ssgsea_group_llm <- function(group_id, top_df_group, methods,
                                 level_name, out_dir = NULL,
                                 celltype_l2 = NA_character_) {
  if (!ENABLE_LLM || is.null(top_df_group) || nrow(top_df_group) == 0) return(NULL)
  evidence_text <- build_ssgsea_evidence_text(top_df_group)
  if (!nzchar(safe_trim(evidence_text))) return(NULL)
  methods  <- unique(as.character(methods))
  methods  <- methods[!is.na(methods) & nzchar(trimws(methods))]
  source_db_label <- collapse_source_db_label(methods)
  dir_tag  <- if (length(methods) > 1) "ssgsea_multi_db" else paste0("ssgsea_", methods[1])
  parts    <- strsplit(group_id, "__", fixed = TRUE)[[1]]
  ct_label <- if (length(parts) >= 2) paste(parts[-1], collapse = "__") else group_id
  ctx_str  <- build_ssgsea_group_context(group_id, methods, level_name)
  std_payload <- run_ssgsea_group_review_llm(
    group_id = group_id,
    ss_ctx = ctx_str,
    evidence_text = evidence_text,
    source_db_label = source_db_label,
    annotated_label = ct_label,
    annotated_level = level_name,
    warnings = character(),
    error_message = NULL
  )
  if (!is.null(std_payload$result)) {
    std <- std_payload$result
    rec <- list(
      celltype_level = level_name, celltype_label = ct_label, celltype_l2 = celltype_l2,
      comparison = group_id, direction = dir_tag, source_db = source_db_label,
      status = "structured", warnings = unique(std_payload$warnings), error = "",
      confidence = placeholder_text(std$confidence, default = "medium"),
      cell_type_judgment = placeholder_text(std$cell_type_judgment),
      annotation_match_degree = placeholder_text(std$annotation_match_degree),
      annotated_l3_correspondence = placeholder_text(std$annotated_l3_correspondence),
      outlier_assessment = placeholder_text(std$outlier_assessment),
      discovery_assessment = placeholder_text(std$discovery_assessment),
      integrated_diagnostic_comment = placeholder_text(std$integrated_diagnostic_comment),
      annotated_l2_dominant = ifelse(is.na(celltype_l2), "", safe_trim(celltype_l2)),
      annotated_l2_dominant_pct = NA_real_,
      annotated_l2_purity_class = "group_defined",
      annotated_l3_dominant = ct_label,
      annotated_l3_dominant_pct = 100,
      annotated_l3_secondary = "",
      annotated_l3_secondary_pct = NA_real_,
      annotated_l3_purity_class = "group_defined",
      annotation_count_summary = sprintf("ssGSEA grouped profile is directly defined by annotated %s label '%s' within group %s.", level_name, ct_label, group_id),
      overview       = placeholder_text(std$overview),
      key_mechanisms = placeholder_text(std$key_mechanisms),
      hypothesis     = placeholder_text(std$hypothesis),
      narrative      = placeholder_text(std$narrative),
      key_drivers    = placeholder_text(collapse_driver_field(std$key_drivers)),
      evidence       = placeholder_text(std$evidence),
      limitations    = placeholder_text(std$limitations),
      raw_text       = evidence_text, raw_result = std_payload$result
    )
  } else {
    rec <- list(
      celltype_level = level_name, celltype_label = ct_label, celltype_l2 = celltype_l2,
      comparison = group_id, direction = dir_tag, source_db = source_db_label,
      status = "error", warnings = unique(std_payload$warnings),
      error = if (!is.null(std_payload$error)) std_payload$error else "standardization failed",
      confidence = "low",
      cell_type_judgment = "Not available from current evidence.",
      annotation_match_degree = "mixed",
      annotated_l3_correspondence = sprintf("当前无法完成对 %s 注释标签 '%s' 的可靠对照判断。", level_name, ct_label),
      outlier_assessment = "当前因 LLM 结构化失败，无法可靠判断是否为异常值或误标。",
      discovery_assessment = "当前因 LLM 结构化失败，无法可靠判断是否为潜在生物学发现。",
      integrated_diagnostic_comment = "当前仅保留 ssGSEA 原始证据，需要重新运行 LLM 才能完成 annotation/discovery/outlier 审阅。",
      annotated_l2_dominant = ifelse(is.na(celltype_l2), "", safe_trim(celltype_l2)),
      annotated_l2_dominant_pct = NA_real_,
      annotated_l2_purity_class = "group_defined",
      annotated_l3_dominant = ct_label,
      annotated_l3_dominant_pct = 100,
      annotated_l3_secondary = "",
      annotated_l3_secondary_pct = NA_real_,
      annotated_l3_purity_class = "group_defined",
      annotation_count_summary = sprintf("ssGSEA grouped profile is directly defined by annotated %s label '%s' within group %s.", level_name, ct_label, group_id),
      overview = "Not available from current evidence.",
      key_mechanisms = "Not available from current evidence.",
      hypothesis = "Not available from current evidence.",
      narrative = "Not available from current evidence.",
      key_drivers = "Not available from current evidence.",
      evidence = safe_trim(evidence_text),
      limitations = "LLM standardization failed.",
      raw_text = evidence_text, raw_result = NULL
    )
  }
  if (!is.null(out_dir)) {
    group_safe <- safe_name(group_id)
    meth_dir   <- file.path(out_dir, sprintf("grp_%s_%s", dir_tag, group_safe))
    dir.create(meth_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rec, file.path(meth_dir, "interpretation.rds"))
    writeLines(
      c(sprintf("# ssGSEA LLM: %s | %s | %s", level_name, group_id, source_db_label), "",
        sprintf("**Overview:** %s",       rec$overview),       "",
        sprintf("**Key Mechanisms:** %s", rec$key_mechanisms), "",
        sprintf("**Hypothesis:** %s",     rec$hypothesis),     "",
        sprintf("**Key Drivers:** %s",    rec$key_drivers)),
      file.path(meth_dir, "interpretation_summary.txt")
    )
  }
  rec
}

get_ssgsea_method_names <- function(ssgsea_results) {
  names(ssgsea_results)[vapply(ssgsea_results, function(x) {
    is.list(x) && !is.null(x$scores) && !is.null(x$z_scores)
  }, logical(1))]
}

# Helper: collect ssGSEA LLM records from results list into flat data.frame
collect_ssgsea_llm_records <- function(ssgsea_results) {
  rows <- list()
  if (!is.null(ssgsea_results[["llm_combined"]]) && length(ssgsea_results[["llm_combined"]]) > 0) {
    llm_sources <- list(ssgsea_multi_db = ssgsea_results[["llm_combined"]])
  } else {
    method_names <- get_ssgsea_method_names(ssgsea_results)
    llm_sources <- stats::setNames(lapply(method_names, function(method) {
      ssgsea_results[[method]][["llm"]]
    }), method_names)
  }
  for (method in names(llm_sources)) {
    llm_list <- llm_sources[[method]]
    if (is.null(llm_list) || length(llm_list) == 0) next
    for (gid in names(llm_list)) {
      rec <- llm_list[[gid]]
      if (is.null(rec)) next
      row_df <- data.frame(
        celltype_level = safe_trim(rec$celltype_level),
        celltype_label = safe_trim(rec$celltype_label),
        celltype_l2    = safe_trim(rec$celltype_l2),
        comparison     = safe_trim(rec$comparison),
        direction      = safe_trim(rec$direction),
        source_db      = safe_trim(rec$source_db),
        status         = safe_trim(rec$status),
        warnings       = paste(rec$warnings, collapse = " | "),
        error          = safe_trim(rec$error),
        overview       = safe_trim(rec$overview),
        key_mechanisms = safe_trim(rec$key_mechanisms),
        hypothesis     = safe_trim(rec$hypothesis),
        narrative      = safe_trim(rec$narrative),
        key_drivers    = safe_trim(rec$key_drivers),
        evidence       = safe_trim(rec$evidence),
        limitations    = safe_trim(rec$limitations),
        raw_text       = safe_trim(rec$raw_text),
        stringsAsFactors = FALSE
      )
      optional_char_fields <- c(
        "confidence", "cell_type_judgment", "annotation_match_degree", "annotated_l3_correspondence",
        "outlier_assessment", "discovery_assessment", "integrated_diagnostic_comment",
        "annotated_l2_dominant", "annotated_l2_purity_class",
        "annotated_l3_dominant", "annotated_l3_secondary", "annotated_l3_purity_class",
        "annotation_count_summary"
      )
      optional_num_fields <- c("annotated_l2_dominant_pct", "annotated_l3_dominant_pct", "annotated_l3_secondary_pct")
      for (nm in optional_char_fields) {
        row_df[[nm]] <- if (!is.null(rec[[nm]])) safe_trim(rec[[nm]]) else ""
      }
      for (nm in optional_num_fields) {
        row_df[[nm]] <- if (!is.null(rec[[nm]])) suppressWarnings(as.numeric(rec[[nm]])) else NA_real_
      }
      rows[[length(rows) + 1]] <- row_df
    }
  }
  if (length(rows) == 0) {
    if (exists("tc_empty_interpretation_records_df", mode = "function")) {
      return(tc_empty_interpretation_records_df())
    }
    return(data.frame())
  }
  dplyr::bind_rows(rows)
}

# Helper: collect CHOIR cluster LLM records
collect_choir_llm_records <- function(choir_ofa) {
  if (!exists("tc_collect_choir_llm_records", mode = "function")) {
    stop("tc_collect_choir_llm_records() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_collect_choir_llm_records(choir_ofa)
}

build_choir_annotation_count_tables <- function(meta_df, choir_col,
                                                annotation_cols = c(L2 = CELLTYPE_L2_COL, L3 = CELLTYPE_L3_COL)) {
  if (exists("tc_build_choir_annotation_count_tables", mode = "function")) {
    return(tc_build_choir_annotation_count_tables(
      meta_df = meta_df,
      choir_col = choir_col,
      annotation_cols = annotation_cols
    ))
  }
  if (is.null(meta_df) || !is.data.frame(meta_df) || nrow(meta_df) == 0) {
    return(list(long = data.frame(), wide = list()))
  }
  if (is.null(choir_col) || !choir_col %in% colnames(meta_df)) {
    return(list(long = data.frame(), wide = list()))
  }

  meta_use <- meta_df %>%
    dplyr::mutate(choir_cluster = as.character(.data[[choir_col]])) %>%
    dplyr::filter(!is.na(choir_cluster), nzchar(trimws(choir_cluster)))
  if (nrow(meta_use) == 0) return(list(long = data.frame(), wide = list()))

  cluster_sizes <- meta_use %>% dplyr::count(choir_cluster, name = "cluster_size")
  long_rows <- list()
  wide_tables <- list()

  for (level_name in names(annotation_cols)) {
    annotation_col <- annotation_cols[[level_name]]
    if (is.null(annotation_col) || !annotation_col %in% colnames(meta_use)) next
    count_df <- meta_use %>%
      dplyr::mutate(annotation = as.character(.data[[annotation_col]])) %>%
      dplyr::filter(!is.na(annotation), nzchar(trimws(annotation))) %>%
      dplyr::count(choir_cluster, annotation, name = "n_cells") %>%
      dplyr::left_join(cluster_sizes, by = "choir_cluster") %>%
      dplyr::mutate(
        annotation_level = level_name,
        pct_in_cluster = ifelse(cluster_size > 0, 100 * n_cells / cluster_size, NA_real_)
      ) %>%
      dplyr::arrange(suppressWarnings(as.numeric(choir_cluster)), choir_cluster, dplyr::desc(n_cells), annotation)
    if (nrow(count_df) == 0) next

    long_rows[[length(long_rows) + 1L]] <- count_df %>%
      dplyr::select(choir_cluster, annotation_level, annotation, n_cells, cluster_size, pct_in_cluster)

    wide_tables[[level_name]] <- count_df %>%
      dplyr::select(choir_cluster, annotation, n_cells) %>%
      tidyr::pivot_wider(names_from = annotation, values_from = n_cells, values_fill = 0) %>%
      dplyr::left_join(cluster_sizes, by = "choir_cluster") %>%
      dplyr::relocate(cluster_size, .after = choir_cluster)
  }

  list(
    long = if (length(long_rows) > 0) dplyr::bind_rows(long_rows) else data.frame(),
    wide = wide_tables
  )
}

write_cluster_annotation_count_tables <- function(meta_df, cluster_col, cluster_dir,
                                                  annotation_cols = c(L2 = CELLTYPE_L2_COL, L3 = CELLTYPE_L3_COL)) {
  tables <- build_choir_annotation_count_tables(
    meta_df = meta_df,
    choir_col = cluster_col,
    annotation_cols = annotation_cols
  )
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  long_path <- file.path(cluster_dir, "cluster_L2_L3_counts_long.csv")
  if (is.data.frame(tables$long) && nrow(tables$long) > 0) {
    long_df <- tables$long %>% dplyr::rename(cluster_id = choir_cluster)
    fwrite(long_df, long_path)
  }
  wide_paths <- list()
  for (level_name in names(tables$wide)) {
    wide_df <- tables$wide[[level_name]]
    if (!is.data.frame(wide_df) || nrow(wide_df) == 0) next
    wide_df <- wide_df %>% dplyr::rename(cluster_id = choir_cluster)
    out_path <- file.path(cluster_dir, sprintf("cluster_%s_counts.csv", level_name))
    fwrite(wide_df, out_path)
    wide_paths[[level_name]] <- out_path
  }
  list(tables = tables, paths = c(wide_paths, long = long_path))
}

run_ofa_marker_enrichment <- function(obj,
                                      focal_cells,
                                      rest_cells,
                                      focal_label,
                                      rest_label = "rest",
                                      out_dir,
                                      plot_title,
                                      group_col_name = "ofa_group",
                                      min_cells_focal = OFA_MIN_CELLS_FOCAL,
                                      min_cells_rest = OFA_MIN_CELLS_REST,
                                      padj_thr = OFA_PADJ_THR,
                                      lfc_thr = OFA_LFC_THR,
                                      top_n = OFA_TOP_N,
                                      max_cells_per_ident = OFA_MAX_CELLS_PER_IDENT) {
  n_focal <- length(focal_cells)
  n_rest <- length(rest_cells)
  if (n_focal < min_cells_focal || n_rest < min_cells_rest) {
    return(list(
      skipped = TRUE,
      reason = sprintf("focal=%d rest=%d (min %d/%d)", n_focal, n_rest, min_cells_focal, min_cells_rest),
      n_focal = n_focal,
      n_rest = n_rest,
      de_table = NULL,
      enrich = list(),
      n_sig = 0L
    ))
  }

  sub_ofa <- subset(obj, cells = c(focal_cells, rest_cells))
  sub_ofa@meta.data[[group_col_name]] <- ifelse(colnames(sub_ofa) %in% focal_cells, as.character(focal_label), as.character(rest_label))
  Idents(sub_ofa) <- group_col_name

  if (n_focal > max_cells_per_ident || n_rest > max_cells_per_ident) {
    set.seed(42)
    keep_focal <- sample(focal_cells, min(n_focal, max_cells_per_ident))
    keep_rest <- sample(rest_cells, min(n_rest, max_cells_per_ident))
    sub_ofa <- subset(sub_ofa, cells = c(keep_focal, keep_rest))
    Idents(sub_ofa) <- group_col_name
  }

  de_ofa <- tryCatch(
    FindMarkers(
      sub_ofa,
      ident.1 = as.character(focal_label),
      ident.2 = as.character(rest_label),
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0
    ),
    error = function(e) {
      cat(sprintf("    [ERROR] FindMarkers: %s\n", e$message))
      NULL
    }
  )
  rm(sub_ofa)
  gc()

  if (is.null(de_ofa) || nrow(de_ofa) == 0) {
    return(list(
      skipped = TRUE,
      reason = "No markers returned",
      n_focal = n_focal,
      n_rest = n_rest,
      de_table = NULL,
      enrich = list(),
      n_sig = 0L
    ))
  }

  de_ofa$gene <- rownames(de_ofa)
  de_ofa$padj <- p.adjust(de_ofa$p_val, method = "BH")
  de_ofa <- de_ofa %>% filter(abs(avg_log2FC) >= lfc_thr) %>% arrange(padj)
  n_sig_ofa <- sum(de_ofa$padj < padj_thr, na.rm = TRUE)

  ensure_dir_exists(out_dir, label = sprintf("OFA output for %s", focal_label))
  fwrite(de_ofa, file.path(out_dir, "markers.csv"))

  de_v <- de_ofa %>%
    mutate(
      sig_flag = padj < padj_thr & abs(avg_log2FC) >= lfc_thr,
      label = ifelse(sig_flag & dplyr::row_number() <= 20, gene, "")
    )
  pv_ofa <- ggplot(de_v, aes(x = avg_log2FC, y = -log10(padj), color = sig_flag)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), guide = "none") +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "blue") +
    labs(title = plot_title, x = "avg_log2FC", y = "-log10(padj)") +
    theme_minimal()
  save_plot(pv_ofa, file.path(out_dir, "volcano"), width = 8, height = 6)

  all_tested_ofa <- unique(tc_filter_gene_symbols(de_ofa$gene))
  ofa_enrich <- list()
  for (dir_name in c("up", "down")) {
    ofa_genes <- if (dir_name == "up") {
      de_ofa %>% filter(padj < padj_thr, avg_log2FC > 0) %>%
        arrange(desc(avg_log2FC)) %>% head(top_n) %>% pull(gene)
    } else {
      de_ofa %>% filter(padj < padj_thr, avg_log2FC < 0) %>%
        arrange(avg_log2FC) %>% head(top_n) %>% pull(gene)
    }
    ofa_genes <- tc_filter_gene_symbols(ofa_genes)
    if (length(ofa_genes) < 5) next
    enr_ofa <- list(
      GO_BP = run_gmt_enrichment(ofa_genes, go_bp_t2g, "GO_BP", all_tested_ofa),
      GO_MF = run_gmt_enrichment(ofa_genes, go_mf_t2g, "GO_MF", all_tested_ofa),
      GO_CC = run_gmt_enrichment(ofa_genes, go_cc_t2g, "GO_CC", all_tested_ofa),
      KEGG = run_gmt_enrichment(ofa_genes, kegg_t2g, "KEGG", all_tested_ofa),
      Hallmark = run_gmt_enrichment(ofa_genes, hallmark_t2g, "Hallmark", all_tested_ofa),
      CellMarker = run_gmt_enrichment(ofa_genes, cellmarker_t2g, "CellMarker", all_tested_ofa),
      PanglaoDB = run_gmt_enrichment(ofa_genes, panglaodb_t2g, "PanglaoDB", all_tested_ofa)
    )
    enr_ofa[[CUSTOM_DB_NAME]] <- run_gmt_enrichment(ofa_genes, custom_t2g, CUSTOM_DB_NAME, all_tested_ofa)
    enr_ofa <- enr_ofa[!sapply(enr_ofa, is.null)]
    ofa_enrich[[dir_name]] <- enr_ofa

    enr_out_ofa <- file.path(out_dir, paste0("enrichment_", dir_name))
    ensure_dir_exists(enr_out_ofa, label = sprintf("OFA enrichment %s %s", focal_label, dir_name))
    for (db in names(enr_ofa)) {
      er <- enr_ofa[[db]]
      if (nrow(as.data.frame(er)) > 0) {
        ensure_dir_exists(enr_out_ofa, label = sprintf("OFA enrichment %s %s", focal_label, dir_name))
        fwrite(as.data.frame(er), file.path(enr_out_ofa, paste0(db, ".csv")))
        tryCatch({
          pdf(file.path(enr_out_ofa, paste0(db, "_dotplot.pdf")), width = 10, height = 8)
          print(dotplot(er, showCategory = 15, title = sprintf("%s | %s", db, plot_title)))
          dev.off()
        }, error = function(e) NULL)
      }
    }
    saveRDS(enr_ofa, file.path(enr_out_ofa, "all_enrichment.rds"))
  }

  bubble_res <- tryCatch(
    if (exists("tc_save_ofa_bubble_plot", mode = "function")) {
      tc_save_ofa_bubble_plot(
        enrich_by_direction = ofa_enrich,
        out_dir = out_dir,
        stem = "enrichment_bubble",
        title = plot_title,
        top_n_per_db = 6L,
        max_term_chars = 80L
      )
    } else {
      list(data = data.frame(), data_csv = NA_character_, plots = list(pdf = NA_character_, png = NA_character_))
    },
    error = function(e) {
      cat(sprintf("    [WARN] OFA bubble plot failed: %s\n", e$message))
      list(data = data.frame(), data_csv = NA_character_, plots = list(pdf = NA_character_, png = NA_character_))
    }
  )

  bubble_overview_base <- file.path(out_dir, "bubbleplot_overview")
  bubble_overview_tsv <- paste0(bubble_overview_base, ".tsv")
  bubble_overview_png <- paste0(bubble_overview_base, ".png")
  bubble_overview_pdf <- paste0(bubble_overview_base, ".pdf")
  bubble_overview_df <- if (is.data.frame(bubble_res$data)) bubble_res$data else data.frame()
  fwrite(bubble_overview_df, bubble_overview_tsv, sep = "\t")
  if (!is.null(bubble_res$plots$png) && nzchar(bubble_res$plots$png) && file.exists(bubble_res$plots$png)) {
    invisible(file.copy(bubble_res$plots$png, bubble_overview_png, overwrite = TRUE))
  }
  if (!is.null(bubble_res$plots$pdf) && nzchar(bubble_res$plots$pdf) && file.exists(bubble_res$plots$pdf)) {
    invisible(file.copy(bubble_res$plots$pdf, bubble_overview_pdf, overwrite = TRUE))
  }
  if (!file.exists(bubble_overview_png) || !file.exists(bubble_overview_pdf)) {
    bubble_placeholder <- ggplot() +
      annotate("text", x = 1, y = 1, label = "No enriched pathways passed current thresholds", size = 5) +
      xlim(0, 2) + ylim(0, 2) +
      theme_void() +
      ggtitle(plot_title)
    save_plot(bubble_placeholder, bubble_overview_base, width = 10, height = 4)
  }

  list(
    skipped = FALSE,
    reason = NULL,
    n_focal = n_focal,
    n_rest = n_rest,
    n_sig = n_sig_ofa,
    de_table = de_ofa,
    enrich = ofa_enrich,
    bubble_data = bubble_res$data,
    bubble_data_csv = bubble_res$data_csv,
    bubble_paths = bubble_res$plots
  )
}

build_choir_cluster_ssgsea_evidence <- function(ssgsea_choir_all, cluster_id,
                                                n_top_per_method = CHOIR_LLM_TOP_SSGSEA_TERMS_PER_DIRECTION) {
  if (!exists("tc_build_choir_cluster_ssgsea_evidence", mode = "function")) {
    stop("tc_build_choir_cluster_ssgsea_evidence() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_build_choir_cluster_ssgsea_evidence(
    ssgsea_choir_all = ssgsea_choir_all,
    cluster_id = cluster_id,
    n_top_per_method = n_top_per_method,
    caller_env = parent.frame()
  )
}

run_choir_cluster_review_llm <- function(cluster_id,
                                         choir_ctx,
                                         evidence_text,
                                         source_db_label,
                                         warnings = character(),
                                         error_message = NULL) {
  if (!exists("tc_run_choir_cluster_review_llm", mode = "function")) {
    stop("tc_run_choir_cluster_review_llm() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_run_choir_cluster_review_llm(
    cluster_id = cluster_id,
    choir_ctx = choir_ctx,
    evidence_text = evidence_text,
    source_db_label = source_db_label,
    warnings = warnings,
    error_message = error_message,
    caller_env = parent.frame()
  )
}

run_ssgsea_group_review_llm <- function(group_id,
                                        ss_ctx,
                                        evidence_text,
                                        source_db_label,
                                        annotated_label,
                                        annotated_level,
                                        warnings = character(),
                                        error_message = NULL) {
  if (!exists("tc_run_ssgsea_group_review_llm", mode = "function")) {
    stop("tc_run_ssgsea_group_review_llm() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_run_ssgsea_group_review_llm(
    group_id = group_id,
    ss_ctx = ss_ctx,
    evidence_text = evidence_text,
    source_db_label = source_db_label,
    annotated_label = annotated_label,
    annotated_level = annotated_level,
    warnings = warnings,
    error_message = error_message,
    caller_env = parent.frame()
  )
}

build_record_screening_prompt <- function(batch_df, family_label) {
  if (!exists("tc_build_record_screening_prompt", mode = "function")) {
    stop("tc_build_record_screening_prompt() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_build_record_screening_prompt(batch_df, family_label)
}

run_record_screening_llm_batch <- function(batch_df, family_label) {
  if (!exists("tc_run_record_screening_llm_batch", mode = "function")) {
    stop("tc_run_record_screening_llm_batch() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_run_record_screening_llm_batch(
    batch_df = batch_df,
    family_label = family_label,
    caller_env = parent.frame()
  )
}

run_record_screening_llm <- function(records_df,
                                     family_label,
                                     batch_size = 10L) {
  if (!exists("tc_run_record_screening_llm", mode = "function")) {
    stop("tc_run_record_screening_llm() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_run_record_screening_llm(
    records_df = records_df,
    family_label = family_label,
    batch_size = batch_size,
    caller_env = parent.frame()
  )
}

write_discovery_screen_markdown <- function(screen_df, path, title, family_label) {
  if (!exists("tc_write_discovery_screen_markdown", mode = "function")) {
    stop("tc_write_discovery_screen_markdown() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_write_discovery_screen_markdown(
    screen_df = screen_df,
    path = path,
    title = title,
    family_label = family_label
  )
}

build_choir_ofa_screening_table <- function(choir_ofa_all) {
  if (!exists("tc_build_choir_ofa_screening_table", mode = "function")) {
    stop("tc_build_choir_ofa_screening_table() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_build_choir_ofa_screening_table(
    choir_ofa_all = choir_ofa_all,
    caller_env = parent.frame()
  )
}

run_and_write_discovery_screen <- function(records_df,
                                           family_label,
                                           output_dir,
                                           prefix,
                                           title,
                                           batch_size = 10L) {
  if (!exists("tc_run_and_write_discovery_screen", mode = "function")) {
    stop("tc_run_and_write_discovery_screen() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_run_and_write_discovery_screen(
    records_df = records_df,
    family_label = family_label,
    output_dir = output_dir,
    prefix = prefix,
    title = title,
    batch_size = batch_size,
    caller_env = parent.frame()
  )
}

merge_directional_gene_pathway_maps <- function(bundle_list) {
  rows <- lapply(names(bundle_list), function(dir_name) {
    bundle <- bundle_list[[dir_name]]
    if (is.null(bundle) || is.null(bundle$gene_pathway_map) || nrow(bundle$gene_pathway_map) == 0) return(NULL)
    df <- bundle$gene_pathway_map
    df$pathway <- sprintf("[%s] %s", toupper(dir_name), df$pathway)
    df
  })
  dplyr::bind_rows(rows)
}

run_pairwise_integrated_llm <- function(group_name, level_name, group_l2_label,
                                        t1, t2, de_df, enrich_by_direction,
                                        comp_dir) {
  if (!ENABLE_LLM || is.null(de_df) || !is.data.frame(de_df) || nrow(de_df) == 0) return(NULL)

  gene_col  <- llm_pick_existing_col(de_df, c("gene", "symbol", "feature", "features", "genes"))
  logfc_col <- llm_pick_existing_col(de_df, c("avg_log2FC", "log2FoldChange", "log2FC", "avg_logFC"))
  gene_fc_all <- if (!is.null(gene_col) && !is.null(logfc_col)) {
    tc_filter_named_gene_fc(stats::setNames(llm_num(de_df[[logfc_col]]), toupper(as.character(de_df[[gene_col]]))))
  } else numeric()

  bundle_list <- list()
  for (dir_name in c("up", "down")) {
    dir_enrich <- enrich_by_direction[[dir_name]]
    if (is.null(dir_enrich) || length(dir_enrich) == 0) next
    gene_fc_dir <- if (dir_name == "up") gene_fc_all[gene_fc_all > 0] else gene_fc_all[gene_fc_all < 0]
    bundle_list[[dir_name]] <- prepare_llm_enrichment_bundle(dir_enrich, gene_fc = gene_fc_dir)
  }

  gene_pathway_map <- merge_directional_gene_pathway_maps(bundle_list)
  gene_pathway_map_text <- paste(format_gene_pathway_map_text(gene_pathway_map), collapse = "\n")

  evidence_rows <- lapply(names(bundle_list), function(dir_name) {
    df <- bundle_list[[dir_name]]$evidence_df
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$direction <- dir_name
    df
  })
  evidence_df <- dplyr::bind_rows(evidence_rows)
  pct_lookup_builder <- NULL
  if (exists("tc_compute_deg_expression_pct_lookup", mode = "function") && exists("obj", inherits = TRUE)) {
    group_col <- if (identical(toupper(level_name), "L3")) CELLTYPE_L3_COL else CELLTYPE_L2_COL
    meta_obj <- obj@meta.data
    cells_t2 <- rownames(meta_obj)[
      as.character(meta_obj[[group_col]]) == as.character(group_name) &
        as.character(meta_obj[[TISSUE_COL]]) == as.character(t2)
    ]
    cells_t1 <- rownames(meta_obj)[
      as.character(meta_obj[[group_col]]) == as.character(group_name) &
        as.character(meta_obj[[TISSUE_COL]]) == as.character(t1)
    ]
    if (length(cells_t2) > 0 || length(cells_t1) > 0) {
      pct_lookup_builder <- function(genes) {
        tc_compute_deg_expression_pct_lookup(
          obj = obj,
          genes = genes,
          cells_1 = cells_t2,
          cells_2 = cells_t1
        )
      }
    }
  }
  evidence_text <- append_nonempty_blocks(
    build_top_deg_context(
      de_df,
      pct_lookup_builder = pct_lookup_builder,
      pct1_label = t2,
      pct2_label = t1
    ),
    build_integrated_directional_enrichment_text(
      evidence_df,
      direction_map = c(
        up = sprintf("higher_in_%s_vs_%s", safe_name(t2), safe_name(t1)),
        down = sprintf("higher_in_%s_vs_%s", safe_name(t1), safe_name(t2))
      )
    )
  )
  if (!nzchar(safe_trim(evidence_text))) return(NULL)

  if (nrow(evidence_df) > 0) {
    fwrite(evidence_df, file.path(comp_dir, "llm_multi_db_evidence_integrated_up_down.csv"))
  }
  if (nrow(gene_pathway_map) > 0) {
    fwrite(gene_pathway_map, file.path(comp_dir, "validated_gene_pathway_map_integrated_up_down.csv"))
  }

  selected_dbs_all <- unique(unlist(lapply(bundle_list, function(bundle) bundle$selected_dbs), use.names = FALSE))
  source_db_label <- if (length(selected_dbs_all) > 0) {
    paste0("integrated_", collapse_source_db_label(selected_dbs_all))
  } else {
    "integrated_multi_db"
  }

  context_str <- build_integrated_tissue_pair_context(
    t1 = t1, t2 = t2,
    celltype_label = group_name,
    celltype_level = level_name,
    gene_pathway_map_text = gene_pathway_map_text
  )
  context_str <- append_nonempty_blocks(context_str, evidence_text)

  std_payload <- standardize_result_with_llm(
    raw_text = "",
    context_str = context_str,
    evidence_text = evidence_text,
    gene_pathway_map_text = gene_pathway_map_text,
    celltype_label = group_name,
    celltype_level = level_name,
    comp_name = paste0(t2, "_vs_", t1),
    dir_name = "integrated_up_down",
    source_db = source_db_label,
    warnings = character(),
    error_message = NULL
  )

  if (!is.null(std_payload$result)) {
    std <- std_payload$result
    rec <- list(
      celltype_level = level_name,
      celltype_label = group_name,
      celltype_l2 = group_l2_label,
      comparison = paste0(t2, "_vs_", t1),
      direction = "integrated_up_down",
      source_db = source_db_label,
      status = "structured",
      warnings = unique(std_payload$warnings),
      error = "",
      overview = placeholder_text(std$overview),
      key_mechanisms = placeholder_text(std$key_mechanisms),
      hypothesis = placeholder_text(std$hypothesis),
      narrative = placeholder_text(std$narrative),
      key_drivers = placeholder_text(collapse_driver_field(std$key_drivers)),
      evidence = placeholder_text(std$evidence),
      limitations = placeholder_text(std$limitations),
      raw_text = evidence_text,
      raw_result = std_payload$result
    )
  } else {
    rec <- list(
      celltype_level = level_name,
      celltype_label = group_name,
      celltype_l2 = group_l2_label,
      comparison = paste0(t2, "_vs_", t1),
      direction = "integrated_up_down",
      source_db = source_db_label,
      status = "error",
      warnings = unique(std_payload$warnings),
      error = ifelse(is.null(std_payload$error), "pairwise integrated LLM failed", std_payload$error),
      overview = placeholder_text(""),
      key_mechanisms = placeholder_text(""),
      hypothesis = placeholder_text(if (isTRUE(LLM_ALLOW_COMPARATIVE_HYPOTHESIS)) "" else "Not applicable."),
      narrative = placeholder_text(""),
      key_drivers = placeholder_text(""),
      evidence = placeholder_text(evidence_text),
      limitations = placeholder_text("Evidence integration failed during LLM standardization."),
      raw_text = evidence_text,
      raw_result = NULL
    )
  }

  saveRDS(rec, file.path(comp_dir, "interpret_agent_integrated_up_down_structured.rds"))
  if (nzchar(rec$raw_text)) {
    writeLines(rec$raw_text, file.path(comp_dir, "interpret_agent_integrated_up_down_raw.txt"))
  }
  rec
}

run_choir_cluster_llm <- function(cluster_id, de_df, ofa_enrich,
                                  ssgsea_choir_all = NULL,
                                  choir_cluster_l2 = NA_character_,
                                  choir_annotation_tables = NULL) {
  if (!exists("tc_run_choir_cluster_llm", mode = "function")) {
    stop("tc_run_choir_cluster_llm() not found. Please load tissue_comparison_advanced_helper_20260408.R first.")
  }
  tc_run_choir_cluster_llm(
    cluster_id = cluster_id,
    de_df = de_df,
    ofa_enrich = ofa_enrich,
    ssgsea_choir_all = ssgsea_choir_all,
    choir_cluster_l2 = choir_cluster_l2,
    choir_annotation_tables = choir_annotation_tables,
    caller_env = parent.frame()
  )
}

# interpret_agent wrappers
run_interpret_agent_once <- function(enrich_obj, context_str, gene_fc = NULL) {
  enrich_obj <- tc_sanitize_enrichment_obj(enrich_obj)
  gene_fc <- tc_filter_named_gene_fc(gene_fc)
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0)
    return(list(result = NULL, warnings = character(), error = "empty enrichment"))
  warn_msgs <- character()
  res <- tryCatch(
    withCallingHandlers(
      clusterProfiler::interpret_agent(
        x = enrich_obj, context = context_str,
        n_pathways = INTERPRET_AGENT_N_PATHWAYS, model = INTERPRET_AGENT_MODEL,
        api_key = DEEPSEEK_API_KEY, add_ppi = INTERPRET_AGENT_ADD_PPI,
        gene_fold_change = gene_fc
      ),
      warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
    ),
    error = function(e) {
      cat(sprintf("    [ERROR] agent: %s\n", e$message))
      structure(list(message = e$message), class = "interpret_agent_error")
    }
  )
  if (inherits(res, "interpret_agent_error"))
    return(list(result = NULL, warnings = unique(warn_msgs), error = res$message))
  list(result = res, warnings = unique(warn_msgs), error = NULL)
}

run_interpret_agent_safe <- function(enrich_obj, context_str, gene_fc = NULL) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0)
    return(list(result = NULL, warnings = character(), error = "empty enrichment"))
  combined_warnings <- character()
  last_payload <- list(result = NULL, warnings = character(), error = "interpret_agent did not run")
  for (attempt in seq_len(INTERPRET_AGENT_MAX_RETRIES)) {
    if (attempt > 1) cat(sprintf("    [INFO] interpret_agent retry %d/%d\n", attempt, INTERPRET_AGENT_MAX_RETRIES))
    payload      <- run_interpret_agent_once(enrich_obj, context_str, gene_fc)
    last_payload <- payload
    if (length(payload$warnings) > 0)
      combined_warnings <- c(combined_warnings, sprintf("attempt %d: %s", attempt, payload$warnings))
    retryable      <- is_retryable_interpret_agent_issue(payload$warnings, payload$error)
    structured     <- !is.null(payload$result) && looks_structured_interpret_agent_result(payload$result)
    raw_text_valid <- !is.null(payload$result) && nzchar(capture_object_text(payload$result))
    if (is.null(payload$error) && (structured || (raw_text_valid && !retryable))) {
      payload$warnings <- unique(c(combined_warnings, payload$warnings))
      return(payload)
    }
    if (attempt < INTERPRET_AGENT_MAX_RETRIES && (retryable || is.null(payload$result))) {
      Sys.sleep(INTERPRET_AGENT_RETRY_SLEEP_SEC); next
    }
    break
  }
  last_payload$warnings <- unique(c(combined_warnings, last_payload$warnings))
  if (is.null(last_payload$error) && !is.null(last_payload$result)) return(last_payload)
  if (is.null(last_payload$error) || !nzchar(last_payload$error))
    last_payload$error <- sprintf("interpret_agent failed after %d attempts", INTERPRET_AGENT_MAX_RETRIES)
  last_payload
}

to_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.character(x) && length(x) == 1) return(sanitize_utf8_text(x))
  if (is.character(x)) return(paste(sanitize_utf8_text(x), collapse = "; "))
  if (is.list(x))
    return(trimws(sanitize_utf8_text(paste(capture.output(str(x, max.level = 2)), collapse = "\n"))))
  trimws(sanitize_utf8_text(as.character(x)))
}

safe_read_rds <- function(path) {
  if (is.null(path) || length(path) != 1 || !nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

if (isTRUE(INITIALIZE_ONLY)) {
  cat("[INIT ONLY] Reference databases and helper functions loaded; skipping pipeline execution.\n")
  stop(structure(
    list(message = "bcell pipeline initialized without execution"),
    class = c("bcell_pipeline_init_only", "error", "condition")
  ))
}

# ==============================================================================
# 5. Load Data
# ==============================================================================

cat(sprintf("=== Loading %s Data ===\n", LINEAGE_DISPLAY))
if (isTRUE(REUSE_PREVIOUS_OUTPUT_SUMMARY) && !isTRUE(SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE) && exists("tc_read_previous_run") && dir.exists(PREVIOUS_OUTPUT_DIR)) {
  cat(sprintf("[INFO] Previous output summary: %s\n", PREVIOUS_OUTPUT_DIR))
  tryCatch(
    tc_print_previous_run_summary(tc_read_previous_run(PREVIOUS_OUTPUT_DIR, include_rds = FALSE, include_markdown = FALSE)),
    error = function(e) cat(sprintf("[WARN] Previous output summary failed: %s\n", e$message))
  )
}

obj <- NULL
if (isTRUE(REUSE_PREVIOUS_FINAL_OBJECT) && file.exists(PREVIOUS_FINAL_OBJECT_RDS)) {
  obj <- safe_read_rds(PREVIOUS_FINAL_OBJECT_RDS)
  if (!is.null(obj)) cat(sprintf("[OK] Loaded previous final object: %s\n", PREVIOUS_FINAL_OBJECT_RDS))
}
if (is.null(obj)) {
  if (!file.exists(H5AD_PATH)) stop(sprintf("File not found: %s", H5AD_PATH))
  obj <- GetSeurat(h5ad_path = H5AD_PATH, prefer_raw = FALSE,
                   prefer_layer_counts = TRUE, validate_counts = TRUE, debug = TRUE)
}
cat(sprintf("[OK] %d cells x %d genes\n", ncol(obj), nrow(obj)))
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA assay missing 'counts' layer.")
cat("[OK] counts layer verified\n\n")

fig_dir   <- file.path(OUTPUT_DIR, "figures")
rpt_dir   <- file.path(OUTPUT_DIR, "reports")
de_dir    <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir    <- file.path(OUTPUT_DIR, "wilcox_exploratory")
de_l3_dir <- file.path(OUTPUT_DIR, "pseudobulk_de_L3")
wx_l3_dir <- file.path(OUTPUT_DIR, "wilcox_exploratory_L3")
for (d in c(fig_dir, rpt_dir, de_dir, wx_dir, de_l3_dir, wx_l3_dir)) {
  ensure_dir_exists(d, label = "pipeline output")
}

previous_rpt_dir <- file.path(PREVIOUS_OUTPUT_DIR, "reports")

load_cached_report_rds <- function(filename, default = NULL, label = NULL) {
  if (!isTRUE(LOAD_EXISTING_STAGE_ARTIFACTS)) return(default)
  label <- if (!is.null(label) && nzchar(label)) label else filename
  current_path  <- file.path(rpt_dir, filename)
  previous_path <- file.path(previous_rpt_dir, filename)
  candidate_paths <- if (isTRUE(LOAD_EXISTING_STAGE_ARTIFACTS_PREFER_CURRENT)) {
    c(current_path, previous_path)
  } else {
    c(previous_path, current_path)
  }
  candidate_paths <- unique(candidate_paths[file.exists(candidate_paths)])
  if (length(candidate_paths) == 0) return(default)
  for (path in candidate_paths) {
    obj <- safe_read_rds(path)
    if (!is.null(obj)) {
      cat(sprintf("[OK] Loaded cached %s: %s\n", label, path))
      return(obj)
    }
  }
  default
}

# ==============================================================================
# 6. Validate Metadata & Apply L2 Remapping
# ==============================================================================

cat("=== Validating Metadata & L2/L3 Annotation ===\n")
meta <- obj@meta.data
required_cols <- c(TISSUE_COL, SAMPLE_COL, L3_SOURCE_COL)
if (USE_EXISTING_L2) required_cols <- c(required_cols, L2_SOURCE_COL)
for (col in required_cols) {
  if (!col %in% colnames(meta)) stop(sprintf("Missing required column: %s", col))
}
l3_vals <- as.character(meta[[L3_SOURCE_COL]])
l2_vals <- NULL
if (USE_EXISTING_L2) {
  l2_vals <- as.character(meta[[L2_SOURCE_COL]])
  if (all(is.na(l2_vals)) || all(trimws(l2_vals) == ""))
    stop(sprintf("Existing L2 column '%s' is empty.", L2_SOURCE_COL))
  cat(sprintf("[OK] Using existing L2 from '%s'\n", L2_SOURCE_COL))
} else {
  l2_vals <- unname(L3_TO_L2_REMAP[l3_vals])
  unmapped <- unique(l3_vals[is.na(l2_vals)])
  if (length(unmapped) > 0)
    stop(sprintf("Unmapped L3 values in '%s': %s\nUpdate L3_TO_L2_REMAP.", L3_SOURCE_COL, paste(unmapped, collapse = ", ")))
  original_l2_backup_col <- paste0(CELLTYPE_L2_COL, "_input")
  if (CELLTYPE_L2_COL %in% colnames(obj@meta.data) && !original_l2_backup_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[original_l2_backup_col]] <- obj@meta.data[[CELLTYPE_L2_COL]]
    cat(sprintf("[OK] Preserved original '%s' as '%s'\n", CELLTYPE_L2_COL, original_l2_backup_col))
  }
}
obj@meta.data[[CELLTYPE_L2_COL]] <- l2_vals
obj@meta.data[[CELLTYPE_L3_COL]] <- l3_vals
l3_l2_mapping_summary <- as.data.frame(table(
  L3 = obj@meta.data[[CELLTYPE_L3_COL]],
  L2 = obj@meta.data[[CELLTYPE_L2_COL]],
  useNA = "no"
), stringsAsFactors = FALSE) %>%
  dplyr::filter(Freq > 0) %>%
  dplyr::arrange(L3, dplyr::desc(Freq), L2)
cat("\nL3 / L2 Summary:\n")
print(table(L3 = l3_vals, L2 = l2_vals, useNA = "ifany"))
cat("\nL2 Distribution:\n"); print(table(obj@meta.data[[CELLTYPE_L2_COL]], useNA = "ifany"))
cat("\nL3 Distribution:\n"); print(table(obj@meta.data[[CELLTYPE_L3_COL]], useNA = "ifany"))
cat("\n")
bad_idx <- is.na(meta[[TISSUE_COL]]) | trimws(as.character(meta[[TISSUE_COL]])) == "" |
           is.na(meta[[SAMPLE_COL]]) | trimws(as.character(meta[[SAMPLE_COL]])) == "" |
           is.na(obj@meta.data[[CELLTYPE_L2_COL]]) | is.na(obj@meta.data[[CELLTYPE_L3_COL]]) |
           trimws(as.character(obj@meta.data[[CELLTYPE_L3_COL]])) == ""
if (sum(bad_idx) > 0) {
  cat(sprintf("[INFO] Dropping %d cells with NA tissue/sample/L2/L3\n", sum(bad_idx)))
  obj <- subset(obj, cells = colnames(obj)[!bad_idx])
}
n_samples <- length(unique(obj@meta.data[[SAMPLE_COL]]))
if (n_samples < 2) stop(sprintf("Only %d unique sample(s). Need >= 2.", n_samples))
meta     <- obj@meta.data
tissues  <- sort(unique(na.omit(meta[[TISSUE_COL]])))
l2_types <- sort(unique(na.omit(meta[[CELLTYPE_L2_COL]])))
l3_types <- sort(unique(na.omit(meta[[CELLTYPE_L3_COL]])))
if (exists("L3_TYPES_INCLUDE") && length(L3_TYPES_INCLUDE) > 0) {
  l3_keep <- unique(as.character(L3_TYPES_INCLUDE))
  l3_types <- intersect(l3_types, l3_keep)
  if (length(l3_types) == 0) {
    stop("L3_TYPES_INCLUDE provided but none of the requested labels are present in the object.")
  }
  cat(sprintf("[INFO] Restricting L3 analysis to %d requested label(s): %s\n",
              length(l3_types), paste(l3_types, collapse = ", ")))
}
cat(sprintf("[OK] Cells=%d | Tissues=%s\n", ncol(obj), paste(tissues, collapse = ", ")))
cat(sprintf("[OK] L2 types (%d): %s\n", length(l2_types), paste(l2_types, collapse = ", ")))
cat(sprintf("[OK] L3 types (%d): %s\n", length(l3_types), paste(l3_types, collapse = ", ")))
cat("\n")
input_had_data_layer <- "data" %in% Layers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
obj@misc$normalization_for_report <- list(
  rerun = TRUE, method = "LogNormalize",
  input_had_data_layer = input_had_data_layer,
  purpose = c("visualization", "marker_analysis"), date = as.character(Sys.time())
)
cat("[OK] data layer ready\n\n")
if (PIPELINE_TEST_MODE) {
  cat("[TEST MODE] Preflight complete.\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# 7. Visualization
# ==============================================================================

marker_panel_viz <- NULL
cat("=== Visualization ===\n")
umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
if (isTRUE(RUN_VISUALIZATION)) {
  if (!is.null(umap_reduction)) {
    tissue_cols_use <- UMAP_TISSUE_COLORS[names(UMAP_TISSUE_COLORS) %in% tissues]
    p1 <- build_umap_plot(obj, umap_reduction, TISSUE_COL,
      title = sprintf("%s - Tissue (%s)", LINEAGE_DISPLAY, umap_reduction),
      cols = tissue_cols_use, width = 11, height = 8)
    save_plot(p1$plot, file.path(fig_dir, "umap_tissue"), width = p1$width, height = p1$height)
    p2 <- build_umap_plot(obj, umap_reduction, CELLTYPE_L2_COL,
      title = sprintf("%s - Cell Type L2 (%s)", LINEAGE_DISPLAY, umap_reduction),
      label = TRUE, width = 14, height = 10)
    save_plot(p2$plot, file.path(fig_dir, "umap_celltype_L2"), width = p2$width, height = p2$height)
    p2s <- build_umap_plot(obj, umap_reduction, CELLTYPE_L2_COL,
      title = sprintf("%s - L2 by Tissue (%s)", LINEAGE_DISPLAY, umap_reduction),
      split_col = TISSUE_COL, label = TRUE,
      width = max(12, 4.5 * length(tissues)), height = 8)
    save_plot(p2s$plot, file.path(fig_dir, "umap_L2_split_tissue"), width = p2s$width, height = p2s$height)
    p2_l3 <- build_umap_plot(obj, umap_reduction, LABEL_COL,
      title = sprintf("%s - Cell Type L3 (%s)", LINEAGE_DISPLAY, umap_reduction),
      label = TRUE, width = 18, height = 12)
    save_plot(p2_l3$plot, file.path(fig_dir, "umap_celltype_L3"), width = p2_l3$width, height = p2_l3$height)
    cat(sprintf("[OK] UMAP saved: %s\n", umap_reduction))
  } else cat("[WARN] No UMAP reduction\n")

  markers_present <- intersect(KNOWN_MARKERS, rownames(obj))
  if (length(markers_present) >= 3) {
    Idents(obj) <- LABEL_COL
    p3 <- DotPlot(obj, features = markers_present) + RotatedAxis() +
      ggtitle(sprintf("%s - Known Markers (L3)", LINEAGE_DISPLAY)) +
      theme(axis.text.x = element_text(size = 7))
    save_plot(p3, file.path(fig_dir, "dotplot_markers"),
              width  = max(12, length(markers_present) * 0.45),
              height = max(6, length(unique(meta[[LABEL_COL]])) * 0.4))
    cat("[OK] Dotplot saved\n")
  }

  if (exists("tc_generate_marker_panel_visualizations", mode = "function") && length(MARKER_PANELS) > 0) {
    marker_panel_viz <- tryCatch(
      tc_generate_marker_panel_visualizations(
        obj = obj,
        marker_panels = MARKER_PANELS,
        fig_dir = file.path(fig_dir, MARKER_PANEL_FIG_SUBDIR),
        group_col = MARKER_PANEL_GROUP_COL,
        lineage_display = LINEAGE_DISPLAY,
        reduction_name = umap_reduction,
        report_dir = rpt_dir,
        panel_prefix = "marker_panel",
        downsample_n = MARKER_PANEL_HEATMAP_CELLS_PER_GROUP,
        verbose = TRUE
      ),
      error = function(e) {
        cat(sprintf("[WARN] Marker-panel visualization failed: %s\n", e$message))
        NULL
      }
    )
  }

  comp_df <- meta %>%
    filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL))) %>%
    count(!!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
    group_by(!!sym(TISSUE_COL)) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
  p4 <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(CELLTYPE_L2_COL))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Tissue", y = "Percentage (%)",
         title = sprintf("%s - L2 Composition by Tissue", LINEAGE_DISPLAY)) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p4, file.path(fig_dir, "composition_tissue_L2"))
  p4b <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = n, fill = !!sym(CELLTYPE_L2_COL))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Tissue", y = "Cell Count",
         title = sprintf("%s - Absolute Count by Tissue", LINEAGE_DISPLAY)) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p4b, file.path(fig_dir, "count_tissue_L2"))

  sample_comp <- meta %>%
    filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL)), !is.na(!!sym(SAMPLE_COL))) %>%
    count(!!sym(SAMPLE_COL), !!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
    group_by(!!sym(SAMPLE_COL)) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
  if (nrow(sample_comp) > 0) {
    facet_ncol <- min(4, length(l2_types)); facet_nrow <- ceiling(length(l2_types) / facet_ncol)
    p4c <- ggplot(sample_comp, aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(TISSUE_COL))) +
      geom_boxplot(outlier.size = 0.5) + geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
      facet_wrap(as.formula(paste("~", CELLTYPE_L2_COL)), scales = "free_y", ncol = facet_ncol) +
      labs(x = "Tissue", y = "Proportion per Sample (%)",
           title = sprintf("%s - Sample-level Composition by Tissue", LINEAGE_DISPLAY)) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    save_plot(p4c, file.path(fig_dir, "composition_sample_level"),
              width = max(10, facet_ncol * 4), height = max(6, facet_nrow * 3.5))
    cat("[OK] Sample-level composition saved\n")
  }

  Idents(obj) <- CELLTYPE_L2_COL
  top_mk <- tryCatch(
    FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                   max.cells.per.ident = 500, test.use = "wilcox"),
    error = function(e) { cat("[WARN] FindAllMarkers failed\n"); NULL }
  )
  if (!is.null(top_mk) && nrow(top_mk) > 0) {
    fwrite(top_mk, file.path(rpt_dir, "all_markers_per_L2.csv"))
    top10  <- top_mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
    obj_ds <- stratified_downsample(obj, CELLTYPE_L2_COL, HEATMAP_CELLS_PER_TYPE)
    obj_ds <- ScaleData(obj_ds, features = unique(top10$gene), verbose = FALSE)
    p5 <- DoHeatmap(obj_ds, features = unique(top10$gene), size = 3) +
      ggtitle(sprintf("%s - Top Markers per L2 (downsampled, visualization only)", LINEAGE_DISPLAY))
    save_plot(p5, file.path(fig_dir, "heatmap_top_markers"), width = 14, height = 10)
    rm(obj_ds); gc()
    cat("[OK] Heatmap saved\n")
  }
} else {
  cat("[INFO] Visualization skipped; existing figures (if any) will be reused.\n")
}
cat("\n")

# ==============================================================================
# 8. Pseudobulk DESeq2 + Multi-Database Enrichment + interpret_agent
# ==============================================================================

run_pairwise_tissue_comparison <- function(obj, group_col, group_values, level_name, out_dir) {
  pb_de_level            <- list()
  enrich_level           <- list()
  agent_level            <- list()
  agent_structured_level <- list()
  pseudobulk_cache       <- build_pseudobulk_cache(obj, group_col = group_col)
  cat(sprintf("=== Pseudobulk DESeq2 Tissue Comparison [%s] ===\n", level_name))
  cat(sprintf("Thresholds: padj < %s, |log2FC| > %s\n", PADJ_THR, LFC_THR))
  if (!is.null(pseudobulk_cache))
    cat(sprintf("[INFO] Pseudobulk cache [%s]: %d columns\n\n", level_name, ncol(pseudobulk_cache$counts)))

  for (group_name in group_values) {
    cat(sprintf("\n>> %s: %s\n", level_name, group_name))
    pb <- aggregate_pseudobulk(pseudobulk_cache, group_name)
    if (is.null(pb)) { cat("  [SKIP] Insufficient pseudobulk\n"); next }
    group_l2_label <- resolve_group_l2_label(obj, group_col, group_name)
    pb_tissues <- unique(na.omit(pb$meta$tissue))
    if (length(pb_tissues) < 2) { cat("  [SKIP] < 2 tissues\n"); next }
    cat(sprintf("  Pseudobulk samples: %d (%s)\n", nrow(pb$meta), paste(pb_tissues, collapse = ", ")))
    pb_de_level[[group_name]]            <- list()
    enrich_level[[group_name]]           <- list()
    agent_level[[group_name]]            <- list()
    agent_structured_level[[group_name]] <- list()

    for (pair in combn(as.character(pb_tissues), 2, simplify = FALSE)) {
      t1 <- pair[1]; t2 <- pair[2]
      comp_name <- paste0(t2, "_vs_", t1)
      cat(sprintf("  DESeq2: %s\n", comp_name))
      res <- run_deseq2_pairwise(pb, t1, t2)
      if (is.null(res)) { cat("    [SKIP] Too few samples per tissue\n"); next }
      pb_de_level[[group_name]][[comp_name]] <- res
      cat(sprintf("    DEGs: %d up, %d down (samples: %d vs %d)\n",
                  res$n_up, res$n_down, res$n_samples_1, res$n_samples_2))
      comp_dir <- file.path(out_dir, safe_name(group_name), safe_name(comp_name))
      ensure_dir_exists(comp_dir, label = sprintf("pairwise comparison %s %s", group_name, comp_name))
      fwrite(res$de_table, file.path(comp_dir, "DESeq2_results.csv"))

      vp <- res$de_table %>% mutate(label = ifelse(sig == "sig" & rank(padj) <= 20, gene, ""))
      pv <- ggplot(vp, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
        geom_point(alpha = 0.5, size = 0.8) +
        scale_color_manual(values = c("sig" = "red", "ns" = "grey70")) +
        geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
        geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "blue") +
        labs(title = sprintf("DESeq2: %s (%s %s)", comp_name, level_name, group_name),
             x = "log2 Fold Change", y = "-log10(padj)") + theme_minimal()
      save_plot(pv, file.path(comp_dir, "volcano"), width = 8, height = 6)

      enrich_level[[group_name]][[comp_name]] <- list()
      agent_level[[group_name]][[comp_name]]  <- list()
      agent_structured_level[[group_name]][[comp_name]] <- list()

      for (dir_name in c("up", "down")) {
        genes <- if (dir_name == "up") {
          res$de_table %>% filter(sig == "sig", log2FoldChange > 0) %>%
            arrange(desc(log2FoldChange)) %>% head(TOP_N_DEG_ENRICHMENT) %>% pull(gene)
        } else {
          res$de_table %>% filter(sig == "sig", log2FoldChange < 0) %>%
            arrange(log2FoldChange) %>% head(TOP_N_DEG_ENRICHMENT) %>% pull(gene)
        }
        genes <- tc_filter_gene_symbols(genes)
        if (length(genes) < 5) { enrich_level[[group_name]][[comp_name]][[dir_name]] <- NULL; next }
        cat(sprintf("    Enrichment [%s]: %d genes\n", dir_name, length(genes)))
        tested  <- res$tested_genes
        enr_list <- list(
          GO_BP = run_gmt_enrichment(genes, go_bp_t2g, "GO_BP", tested),
          GO_MF = run_gmt_enrichment(genes, go_mf_t2g, "GO_MF", tested),
          GO_CC = run_gmt_enrichment(genes, go_cc_t2g, "GO_CC", tested),
          KEGG  = run_gmt_enrichment(genes, kegg_t2g,  "KEGG",  tested),
          Hallmark   = run_gmt_enrichment(genes, hallmark_t2g,   "Hallmark",   tested),
          CellMarker = run_gmt_enrichment(genes, cellmarker_t2g, "CellMarker", tested),
          PanglaoDB  = run_gmt_enrichment(genes, panglaodb_t2g,  "PanglaoDB",  tested)
        )
        enr_list[[CUSTOM_DB_NAME]] <- run_gmt_enrichment(genes, custom_t2g, CUSTOM_DB_NAME, tested)
        enr_list <- enr_list[!sapply(enr_list, is.null)]
        enrich_level[[group_name]][[comp_name]][[dir_name]] <- enr_list

        enr_out <- file.path(comp_dir, paste0("enrichment_", dir_name))
        ensure_dir_exists(enr_out, label = sprintf("pairwise enrichment %s %s", comp_name, dir_name))
        for (db in names(enr_list)) {
          er <- enr_list[[db]]
          if (nrow(as.data.frame(er)) > 0) {
            fwrite(as.data.frame(er), file.path(enr_out, paste0(db, ".csv")))
            tryCatch({
              pdf(file.path(enr_out, paste0(db, "_dotplot.pdf")), width = 10, height = 8)
              print(dotplot(er, showCategory = 15,
                            title = sprintf("%s %s (%s %s %s)", db, dir_name, level_name, group_name, comp_name)))
              dev.off()
            }, error = function(e) NULL)
          }
        }
        saveRDS(enr_list, file.path(enr_out, "all_enrichment.rds"))
      }

      if (ENABLE_LLM) {
        cat("    LLM [integrated_up_down] (combine DE + multi-db enrichment across both directions)\n")
        pairwise_llm_rec <- tryCatch(
          run_pairwise_integrated_llm(
            group_name = group_name,
            level_name = level_name,
            group_l2_label = group_l2_label,
            t1 = t1,
            t2 = t2,
            de_df = res$de_table,
            enrich_by_direction = enrich_level[[group_name]][[comp_name]],
            comp_dir = comp_dir
          ),
          error = function(e) {
            cat(sprintf("    [WARN] pairwise integrated LLM failed: %s\n", e$message))
            NULL
          }
        )
        if (!is.null(pairwise_llm_rec)) {
          agent_level[[group_name]][[comp_name]][["integrated_up_down"]] <- pairwise_llm_rec$raw_result
          agent_structured_level[[group_name]][[comp_name]][["integrated_up_down"]] <- pairwise_llm_rec
          cat(sprintf("    [OK] pairwise integrated LLM status=%s\n", pairwise_llm_rec$status))
          Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
        }
      }
    }
  }
  list(pb_de = pb_de_level, enrich = enrich_level,
       agent = agent_level, agent_structured = agent_structured_level)
}

if (isTRUE(RUN_PAIRWISE_DE)) {
  comparison_l2 <- run_pairwise_tissue_comparison(
    obj = obj, group_col = CELLTYPE_L2_COL, group_values = l2_types, level_name = "L2", out_dir = de_dir)
  pb_de_all            <- comparison_l2$pb_de
  enrich_all           <- comparison_l2$enrich
  agent_all            <- comparison_l2$agent
  agent_structured_all <- comparison_l2$agent_structured

  comparison_l3 <- run_pairwise_tissue_comparison(
    obj = obj, group_col = CELLTYPE_L3_COL, group_values = l3_types, level_name = "L3", out_dir = de_l3_dir)
  pb_de_l3_all            <- comparison_l3$pb_de
  enrich_l3_all           <- comparison_l3$enrich
  agent_l3_all            <- comparison_l3$agent
  agent_structured_l3_all <- comparison_l3$agent_structured
  saveRDS(pb_de_all,               file.path(rpt_dir, "pseudobulk_de_all.rds"))
  saveRDS(pb_de_l3_all,            file.path(rpt_dir, "pseudobulk_de_L3_all.rds"))
  saveRDS(enrich_all,              file.path(rpt_dir, "enrichment_all.rds"))
  saveRDS(enrich_l3_all,           file.path(rpt_dir, "enrichment_L3_all.rds"))
  saveRDS(agent_all,               file.path(rpt_dir, "interpret_agent_all.rds"))
  saveRDS(agent_l3_all,            file.path(rpt_dir, "interpret_agent_L3_all.rds"))
  saveRDS(agent_structured_all,    file.path(rpt_dir, "interpret_agent_structured_all.rds"))
  saveRDS(agent_structured_l3_all, file.path(rpt_dir, "interpret_agent_structured_L3_all.rds"))
  cat("[OK] Pairwise stage checkpoints saved\n")
} else {
  cat("=== Pseudobulk DESeq2 Tissue Comparison [SKIPPED; loading cached artifacts] ===\n")
  pb_de_all <- load_cached_report_rds("pseudobulk_de_all.rds", default = list(), label = "L2 pseudobulk DE")
  pb_de_l3_all <- load_cached_report_rds("pseudobulk_de_L3_all.rds", default = list(), label = "L3 pseudobulk DE")
  enrich_all <- load_cached_report_rds("enrichment_all.rds", default = list(), label = "L2 enrichment")
  enrich_l3_all <- load_cached_report_rds("enrichment_L3_all.rds", default = list(), label = "L3 enrichment")
  agent_all <- load_cached_report_rds("interpret_agent_all.rds", default = list(), label = "L2 interpret_agent")
  agent_l3_all <- load_cached_report_rds("interpret_agent_L3_all.rds", default = list(), label = "L3 interpret_agent")
  agent_structured_all <- load_cached_report_rds("interpret_agent_structured_all.rds", default = list(), label = "L2 structured interpret_agent")
  agent_structured_l3_all <- load_cached_report_rds("interpret_agent_structured_L3_all.rds", default = list(), label = "L3 structured interpret_agent")
  comparison_l2 <- NULL
  comparison_l3 <- NULL
}

# ==============================================================================
# 8.5 Grouped ssGSEA (average-expression) + [LLM-4] Per-group LLM
# ==============================================================================

run_grouped_ssgsea <- function(obj, group_col, level_name, methods, n_top, n_heatmap) {
  results_all <- list()
  group_key   <- paste0("ssgsea_group_", tolower(level_name))
  file_stem   <- level_file_stem("ssgsea", level_name)
  obj@meta.data[[group_key]] <- paste0(obj@meta.data[[TISSUE_COL]], "__", obj@meta.data[[group_col]])
  group_l2_map <- obj@meta.data %>%
    dplyr::transmute(
      group_id = .data[[group_key]],
      celltype_label = as.character(.data[[group_col]]),
      celltype_l2 = as.character(.data[[CELLTYPE_L2_COL]])
    ) %>%
    dplyr::group_by(group_id, celltype_label) %>%
    dplyr::summarise(
      celltype_l2 = {
        vals <- unique(celltype_l2[!is.na(celltype_l2) & nzchar(trimws(celltype_l2))])
        if (length(vals) == 1) vals else NA_character_
      },
      .groups = "drop"
    )
  ssgsea_groups <- sort(unique(obj@meta.data[[group_key]]))
  cat(sprintf("\n=== Grouped ssGSEA (average-expression; tissue x %s) ===\n", level_name))
  cat(sprintf("Methods: %s\n\n", paste(methods, collapse = ", ")))
  cat(sprintf("[INFO] ssGSEA groups (%d): %s\n\n", length(ssgsea_groups), paste(ssgsea_groups, collapse = ", ")))

  for (method in methods) {
    cat(sprintf("--- ssGSEA method [%s]: %s ---\n", level_name, method))
    selected_t2g <- switch(
      method,
      hallmark = hallmark_t2g, go_bp = go_bp_t2g, go_mf = go_mf_t2g,
      go_cc = go_cc_t2g, kegg = kegg_t2g,
      custom = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT)) {
        tryCatch(read.gmt(SSGSEA_CUSTOM_GMT) %>% mutate(gene = toupper(gene)),
                 error = function(e) { cat("[WARN] custom GMT failed\n"); NULL })
      } else NULL, NULL
    )
    if (is.null(selected_t2g) || nrow(selected_t2g) == 0) {
      cat(sprintf("[WARN] Method '%s': no gene sets, skipping\n", method)); next
    }
    features <- rownames(obj)
    gs_list  <- term2gene_to_list(selected_t2g)
    gs_list  <- map_gene_sets_to_features(gs_list, features)
    gs_list  <- filter_gs_size(gs_list)
    if (length(gs_list) == 0) {
      cat(sprintf("[WARN] Method '%s': no valid gene sets after filtering\n", method)); next
    }
    cat(sprintf("[OK] %d gene sets prepared\n", length(gs_list)))
    needed_genes <- unique(unlist(gs_list, use.names = FALSE))
    cat(sprintf("[INFO] Computing group-average expression for %d genes...\n", length(needed_genes)))
    avg_expr <- tryCatch(
      AverageExpression(obj, assays = "RNA", slot = "data",
                        group.by = group_key, features = needed_genes, verbose = FALSE)[["RNA"]],
      error = function(e) { cat(sprintf("[ERROR] AverageExpression failed: %s\n", e$message)); NULL }
    )
    if (is.null(avg_expr)) { cat("[WARN] Skipping this method\n"); next }
    cat("[INFO] Running ssGSEA...\n")
    bp <- BiocParallel::SnowParam(workers = N_CORES, type = "SOCK", progressbar = FALSE)
      ssgsea_scores <- tryCatch({
        if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
          param <- GSVA::ssgseaParam(exprData = as.matrix(avg_expr), geneSets = gs_list,
                                     alpha = 0.25, normalize = TRUE, minSize = 10, maxSize = 500)
          GSVA::gsva(param, BPPARAM = bp, verbose = FALSE)
        } else {
          GSVA::gsva(as.matrix(avg_expr), gs_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
        }
      }, error = function(e) { cat(sprintf("[ERROR] ssGSEA failed: %s\n", e$message)); NULL },
      finally = {
        try(BiocParallel::bpstop(bp), silent = TRUE)
      })
    if (is.null(ssgsea_scores)) { cat("[WARN] Skipping this method\n"); next }
    ssgsea_z <- t(scale(t(ssgsea_scores)))
    cat(sprintf("[OK] ssGSEA done: %d pathways x %d groups\n", nrow(ssgsea_scores), ncol(ssgsea_scores)))

    results_all[[method]] <- list(scores = ssgsea_scores, z_scores = ssgsea_z, gene_sets = gs_list)

    saveRDS(ssgsea_scores, file.path(rpt_dir, sprintf("%s_scores_%s.rds", file_stem, method)))
    saveRDS(ssgsea_z,      file.path(rpt_dir, sprintf("%s_z_%s.rds",      file_stem, method)))

    top_df <- extract_directional_ssgsea_top_rows(
      ssgsea_scores, ssgsea_z,
      top_n = n_top,
      group_field = "group",
      level_name = level_name,
      method = method
    )
    fwrite(top_df, file.path(rpt_dir, sprintf("%s_top_pathways_%s.csv", file_stem, method)))

    # Heatmap
    mean_abs_z <- rowMeans(abs(ssgsea_z), na.rm = TRUE)
    top_paths  <- names(sort(mean_abs_z, decreasing = TRUE))[seq_len(min(n_heatmap, length(mean_abs_z)))]
    heat_mat   <- ssgsea_z[top_paths, , drop = FALSE]
    heat_mat[is.nan(heat_mat) | is.infinite(heat_mat)] <- 0
    if (nrow(heat_mat) >= 3 && ncol(heat_mat) >= 2) {
      tryCatch({
        heat_path <- file.path(fig_dir, sprintf("%s_heatmap_%s", file_stem, method))
        ht <- pheatmap::pheatmap(
          heat_mat, cluster_rows = TRUE, cluster_cols = TRUE,
          color  = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
          breaks = seq(-3, 3, length.out = 101),
          main   = sprintf("%s ssGSEA Z-score (%s) — top %d pathways by tissue x %s",
                           LINEAGE_DISPLAY, method, nrow(heat_mat), level_name),
          fontsize_row = 7, fontsize_col = 9, cellwidth = 22, cellheight = 10,
          filename = paste0(heat_path, ".pdf"),
          width  = max(10, ncol(heat_mat) * 1.8),
          height = max(8,  nrow(heat_mat) * 0.35 + 3)
        )
        grDevices::png(paste0(heat_path, ".png"),
                       width = max(10, ncol(heat_mat) * 1.8),
                       height = max(8, nrow(heat_mat) * 0.35 + 3),
                       units = "in", res = 300)
        grid::grid.newpage(); grid::grid.draw(ht$gtable)
        grDevices::dev.off()
        cat(sprintf("[OK] Heatmap saved: %s_heatmap_%s\n", file_stem, method))
      }, error = function(e) cat(sprintf("[WARN] Heatmap failed for %s: %s\n", method, e$message)))
    }
    results_all[[method]][["top_df"]] <- top_df
    results_all[[method]][["llm"]] <- list()

    cat(sprintf("[OK] Method '%s' complete\n\n", method))
  }

  ssgsea_llm_combined <- list()
  ssgsea_method_names <- get_ssgsea_method_names(results_all)
  if (ENABLE_LLM && length(ssgsea_method_names) > 0) {
    cat(sprintf("[LLM-4] Interpreting ssGSEA groups with combined multi-database evidence (level=%s)...\n", level_name))
    ssgsea_llm_out <- file.path(rpt_dir, sprintf("ssgsea_llm_%s", tolower(level_name)))
    dir.create(ssgsea_llm_out, recursive = TRUE, showWarnings = FALSE)
    top_df_all <- dplyr::bind_rows(lapply(ssgsea_method_names, function(method) {
      results_all[[method]][["top_df"]]
    }))
    all_group_ids <- unique(top_df_all$group)
    for (gid in all_group_ids) {
      gid_df   <- top_df_all[top_df_all$group == gid, , drop = FALSE]
      gid_meta <- group_l2_map[group_l2_map$group_id == gid, , drop = FALSE]
      gid_methods <- unique(as.character(gid_df$method))
      rec_gid <- tryCatch(
        run_ssgsea_group_llm(
          gid, gid_df, gid_methods, level_name, ssgsea_llm_out,
          celltype_l2 = if (nrow(gid_meta) >= 1) gid_meta$celltype_l2[1] else NA_character_
        ),
        error = function(e) {
          cat(sprintf("    [WARN] ssGSEA multi-db LLM failed for %s: %s\n", gid, e$message)); NULL
        }
      )
      if (!is.null(rec_gid)) {
        ssgsea_llm_combined[[gid]] <- rec_gid
        cat(sprintf("    [OK] %s | dbs=%s | status=%s\n",
                    gid, paste(gid_methods, collapse = ", "), rec_gid$status))
      }
      Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
    }
    cat(sprintf("[OK] Combined ssGSEA LLM interpretations: %d / %d groups\n",
                length(ssgsea_llm_combined), length(all_group_ids)))
  }
  results_all[["llm_combined"]] <- ssgsea_llm_combined

  cat(sprintf("[SUMMARY] Completed %d ssGSEA methods for %s\n", length(ssgsea_method_names), level_name))
  for (mn in ssgsea_method_names) {
    cat(sprintf("  - %s: %d pathways x %d groups\n",
                mn, nrow(results_all[[mn]]$scores), ncol(results_all[[mn]]$scores)))
  }
  if (ENABLE_LLM)
    cat(sprintf("  - combined multi-db ssGSEA LLM: %d records\n", length(ssgsea_llm_combined)))
  cat("\n")
  results_all
}

ssgsea_results_all    <- list()
ssgsea_results_l3_all <- list()

if (RUN_SSGSEA) {
  ssgsea_results_all <- run_grouped_ssgsea(
    obj = obj, group_col = CELLTYPE_L2_COL, level_name = "L2",
    methods = SSGSEA_METHODS, n_top = SSGSEA_N_TOP, n_heatmap = SSGSEA_N_HEATMAP)
  ssgsea_results_l3_all <- run_grouped_ssgsea(
    obj = obj, group_col = CELLTYPE_L3_COL, level_name = "L3",
    methods = SSGSEA_METHODS, n_top = SSGSEA_N_TOP, n_heatmap = SSGSEA_N_HEATMAP)
  saveRDS(ssgsea_results_all,    file.path(rpt_dir, "ssgsea_results_all.rds"))
  saveRDS(ssgsea_results_l3_all, file.path(rpt_dir, "ssgsea_results_L3_all.rds"))
  ssgsea_primary_method <- if ("hallmark" %in% names(ssgsea_results_all)) "hallmark"
                           else if (length(ssgsea_results_all) > 0) names(ssgsea_results_all)[1]
                           else NULL
  cat(sprintf("[INFO] Primary ssGSEA method: %s\n\n",
              if (is.null(ssgsea_primary_method)) "none" else ssgsea_primary_method))
} else {
  cat("[INFO] RUN_SSGSEA=FALSE; loading cached grouped ssGSEA artifacts when available\n")
  ssgsea_results_all <- load_cached_report_rds("ssgsea_results_all.rds", default = list(), label = "L2 ssGSEA")
  ssgsea_results_l3_all <- load_cached_report_rds("ssgsea_results_L3_all.rds", default = list(), label = "L3 ssGSEA")
  ssgsea_primary_method <- if ("hallmark" %in% names(ssgsea_results_all)) "hallmark"
                           else if (length(ssgsea_results_all) > 0) names(ssgsea_results_all)[1]
                           else NULL
  cat(sprintf("[INFO] Primary ssGSEA method: %s\n\n",
              if (is.null(ssgsea_primary_method)) "none" else ssgsea_primary_method))
}

# ==============================================================================
# 8.6 Clustering backend (CHOIR / Leiden) + Per-cluster ssGSEA + OFA + [LLM-5]
# ============================================================================== 

cluster_backend <- toupper(trimws(as.character(CLUSTER_BACKEND)))
cluster_backend_label <- if (identical(cluster_backend, "LEIDEN")) "Leiden" else "CHOIR"
cluster_output_prefix <- tolower(cluster_backend_label)
cluster_llm_tsv_filename <- sprintf("%s_llm_structured.tsv", cluster_output_prefix)
cluster_llm_rds_filename <- sprintf("%s_llm_structured_all.rds", cluster_output_prefix)
cluster_llm_markdown_filename <- sprintf("LLM_%s_INTERPRETATION.md", toupper(cluster_output_prefix))
cluster_llm_markdown_path <- file.path(OUTPUT_DIR, cluster_llm_markdown_filename)
cluster_results_rds_name <- sprintf("%s_results_all.rds", cluster_output_prefix)
cluster_ofa_rds_name <- sprintf("%s_ofa_all.rds", cluster_output_prefix)
cluster_umap_stem <- paste0(cluster_output_prefix, "_umap")
cluster_umap_split_stem <- paste0(cluster_output_prefix, "_umap_split_tissue")

choir_results_all <- list()
choir_ofa_all <- list()
empty_l3_ofa_summary_df <- function() {
  data.frame(
    family = character(),
    celltype_label = character(),
    celltype_l2 = character(),
    comparison = character(),
    n_focal = integer(),
    n_rest = integer(),
    n_sig = integer(),
    bubble_png = character(),
    bubble_pdf = character(),
    output_dir = character(),
    stringsAsFactors = FALSE
  )
}
l3_ofa_vs_rest_all <- list()
l3_ofa_same_l2_all <- list()
l3_ofa_inter_tissue_all <- list()
l3_ofa_vs_rest_summary_df <- empty_l3_ofa_summary_df()
l3_ofa_same_l2_summary_df <- empty_l3_ofa_summary_df()
l3_ofa_inter_tissue_summary_df <- empty_l3_ofa_summary_df()

if (RUN_CLUSTERING) {
  cat(sprintf("\n=== %s Clustering (scANVI/scVI latent) ===\n", cluster_backend_label))
  choir_reduction <- pick_reduction(obj, CHOIR_REDUCTION_CANDIDATES)
  if (is.null(choir_reduction)) {
    cat(sprintf("[INFO] No scANVI/scVI/harmony reduction; running PCA as %s fallback\n", cluster_backend_label))
    obj <- tryCatch(RunPCA(obj, verbose = FALSE),
                    error = function(e) { cat(sprintf("[ERROR] PCA: %s\n", e$message)); obj })
    choir_reduction <- if ("pca" %in% Reductions(obj)) "pca" else NULL
  }

  if (is.null(choir_reduction)) {
    cat(sprintf("[WARN] No suitable reduction for %s clustering; skipping Section 8.6\n", cluster_backend_label))
    RUN_CLUSTERING <- FALSE
  } else {
    cat(sprintf("[INFO] %s reduction: %s\n", cluster_backend_label, choir_reduction))
    choir_embedding <- tryCatch(
      as.matrix(Embeddings(obj, reduction = choir_reduction)),
      error = function(e) {
        cat(sprintf("[ERROR] Failed to extract embeddings for '%s': %s\n", choir_reduction, e$message)); NULL
      }
    )
    if (is.null(choir_embedding)) {
      cat(sprintf("[WARN] %s embedding unavailable; skipping Section 8.6\n", cluster_backend_label))
      RUN_CLUSTERING <- FALSE
    }

    choir_var_features <- VariableFeatures(obj)
    choir_var_cap <- min(as.integer(CHOIR_VAR_FEATURES_MAX), nrow(obj))
    if (length(choir_var_features) == 0 && isTRUE(CHOIR_FIND_VAR_FEATURES_IF_MISSING)) {
      cat(sprintf("[INFO] %s var_features missing; computing %d features via FindVariableFeatures(%s)\n",
                  cluster_backend_label, choir_var_cap, CHOIR_VAR_FEATURES_METHOD))
      obj <- tryCatch(
        FindVariableFeatures(obj, assay = DefaultAssay(obj), selection.method = CHOIR_VAR_FEATURES_METHOD,
                             nfeatures = choir_var_cap, verbose = FALSE),
        error = function(e) {
          cat(sprintf("[WARN] FindVariableFeatures for %s failed: %s\n", cluster_backend_label, e$message))
          obj
        }
      )
      choir_var_features <- VariableFeatures(obj)
    }
    if (length(choir_var_features) == 0) {
      choir_var_features <- head(rownames(obj), choir_var_cap)
      cat(sprintf("[WARN] %s var_features still unavailable; fallback to first %d genes\n",
                  cluster_backend_label, length(choir_var_features)))
    } else if (length(choir_var_features) > choir_var_cap) {
      choir_var_features <- choir_var_features[seq_len(choir_var_cap)]
      cat(sprintf("[INFO] %s var_features capped at %d features\n",
                  cluster_backend_label, length(choir_var_features)))
    } else {
      cat(sprintf("[INFO] %s var_features: %d features\n", cluster_backend_label, length(choir_var_features)))
    }

    cluster_col <- NULL
    cluster_dir <- file.path(rpt_dir, cluster_output_prefix)
  ensure_dir_exists(cluster_dir, label = sprintf("%s cluster output", cluster_backend_label))
    obj_clustered <- NULL

    if (RUN_CLUSTERING) {
      if (identical(cluster_backend, "LEIDEN")) {
        cat(sprintf("[INFO] Running Leiden (resolution=%.3f, dims=%d, k=%d)...\n",
                    as.numeric(LEIDEN_RESOLUTION), as.integer(LEIDEN_N_DIMS), as.integer(LEIDEN_K_PARAM)))
        cluster_tag <- gsub("_+", "_", safe_name(format(as.numeric(LEIDEN_RESOLUTION), scientific = FALSE, trim = TRUE)))
        cluster_col <- paste0("LEIDEN_clusters_", cluster_tag)
        dims_use <- seq_len(min(as.integer(LEIDEN_N_DIMS), ncol(choir_embedding)))
        obj_clustered <- tryCatch({
          obj_tmp <- FindNeighbors(
            obj,
            reduction = choir_reduction,
            dims = dims_use,
            k.param = as.integer(LEIDEN_K_PARAM),
            verbose = FALSE
          )
          obj_tmp <- FindClusters(
            obj_tmp,
            resolution = as.numeric(LEIDEN_RESOLUTION),
            algorithm = as.integer(LEIDEN_ALGORITHM),
            verbose = FALSE
          )
          if ("seurat_clusters" %in% colnames(obj_tmp@meta.data)) {
            obj_tmp@meta.data[[cluster_col]] <- as.character(obj_tmp@meta.data[["seurat_clusters"]])
          }
          obj_tmp
        }, error = function(e) {
          cat(sprintf("[ERROR] Leiden failed: %s\n", e$message))
          NULL
        })
      } else {
        cat(sprintf("[INFO] Running CHOIR (alpha=%.3f, n_cores=%d)...\n", CHOIR_ALPHA, CHOIR_N_CORES))
        cluster_col <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
        choir_args <- list(
          object = obj,
          use_assay = "RNA",
          reduction = choir_embedding,
          var_features = choir_var_features,
          n_cores = CHOIR_N_CORES,
          alpha = CHOIR_ALPHA,
          n_iterations = as.integer(CHOIR_N_ITERATIONS),
          n_trees = as.integer(CHOIR_N_TREES),
          use_variance = isTRUE(CHOIR_USE_VARIANCE),
          min_accuracy = as.numeric(CHOIR_MIN_ACCURACY),
          min_connections = as.integer(CHOIR_MIN_CONNECTIONS),
          max_repeat_errors = as.integer(CHOIR_MAX_REPEAT_ERRORS),
          distance_approx = isTRUE(CHOIR_DISTANCE_APPROX),
          distance_awareness = as.numeric(CHOIR_DISTANCE_AWARENESS),
          collect_all_metrics = isTRUE(CHOIR_COLLECT_ALL_METRICS),
          sample_max = CHOIR_SAMPLE_MAX,
          downsampling_rate = CHOIR_DOWNSAMPLING_RATE,
          min_reads = CHOIR_MIN_READS,
          max_clusters = CHOIR_MAX_CLUSTERS,
          min_cluster_depth = as.integer(CHOIR_MIN_CLUSTER_DEPTH),
          normalization_method = CHOIR_NORMALIZATION_METHOD,
          subtree_reductions = isTRUE(CHOIR_SUBTREE_REDUCTIONS),
          reduction_method = CHOIR_REDUCTION_METHOD,
          reduction_params = CHOIR_REDUCTION_PARAMS,
          batch_correction_method = CHOIR_BATCH_CORRECTION_METHOD,
          batch_correction_params = CHOIR_BATCH_CORRECTION_PARAMS,
          batch_labels = CHOIR_BATCH_LABELS,
          neighbor_params = CHOIR_NEIGHBOR_PARAMS,
          cluster_params = CHOIR_CLUSTER_PARAMS,
          random_seed = 42
        )
        choir_args <- choir_args[!vapply(choir_args, is.null, logical(1))]
        obj_clustered <- tryCatch(
          do.call(CHOIR::CHOIR, choir_args),
          error = function(e) { cat(sprintf("[ERROR] CHOIR failed: %s\n", e$message)); NULL }
        )
      }
    }

    if (!is.null(obj_clustered)) obj <- obj_clustered
    rm(obj_clustered, choir_embedding, choir_var_features)
    gc()

    if (!is.null(cluster_col) && !cluster_col %in% colnames(obj@meta.data)) {
      fallback <- if (identical(cluster_backend, "LEIDEN")) {
        grep("^LEIDEN_clusters", colnames(obj@meta.data), value = TRUE)
      } else {
        grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
      }
      cluster_col <- if (length(fallback) > 0) fallback[1] else NULL
      if (!is.null(cluster_col)) cat(sprintf("[INFO] Using %s cluster column: %s\n", cluster_backend_label, cluster_col))
    }

    if (is.null(cluster_col)) {
      cat(sprintf("[WARN] %s cluster column not found; skipping Section 8.6\n", cluster_backend_label))
    } else {
      choir_clusters <- sort(unique(na.omit(as.character(obj@meta.data[[cluster_col]]))))
      cat(sprintf("[OK] %s: %d clusters\n", cluster_backend_label, length(choir_clusters)))
      fwrite(
        data.frame(cell = colnames(obj), cluster_id = as.character(obj@meta.data[[cluster_col]])),
        file.path(cluster_dir, sprintf("%s_clusters.csv", cluster_output_prefix))
      )
      choir_annotation_tables <- write_cluster_annotation_count_tables(
        meta_df = obj@meta.data,
        cluster_col = cluster_col,
        cluster_dir = cluster_dir,
        annotation_cols = c(L2 = CELLTYPE_L2_COL, L3 = CELLTYPE_L3_COL)
      )$tables

      if (!is.null(umap_reduction)) {
        cluster_param_text <- if (identical(cluster_backend, "LEIDEN")) {
          sprintf("resolution=%.3f", as.numeric(LEIDEN_RESOLUTION))
        } else {
          sprintf("alpha=%.3f", CHOIR_ALPHA)
        }
        p_cluster <- build_umap_plot(
          obj,
          umap_reduction,
          cluster_col,
          title = sprintf("%s - %s Clusters (%s, %s)", LINEAGE_DISPLAY, cluster_backend_label, cluster_param_text, choir_reduction),
          label = TRUE,
          width = 14,
          height = 10
        )
        save_plot(p_cluster$plot, file.path(fig_dir, cluster_umap_stem), p_cluster$width, p_cluster$height)
        p_cluster_s <- build_umap_plot(
          obj,
          umap_reduction,
          cluster_col,
          title = sprintf("%s - %s by Tissue (%s)", LINEAGE_DISPLAY, cluster_backend_label, umap_reduction),
          split_col = TISSUE_COL,
          label = TRUE,
          width = max(12, 4.5 * length(tissues)),
          height = 8
        )
        save_plot(p_cluster_s$plot, file.path(fig_dir, cluster_umap_split_stem), p_cluster_s$width, p_cluster_s$height)
        cat(sprintf("[OK] %s UMAPs saved\n", cluster_backend_label))
      }

      # 8.6a ssGSEA per cluster
      ssgsea_choir_all <- list()
      if (RUN_CLUSTER_SSGSEA && length(SSGSEA_CHOIR_METHODS) > 0) {
        cat(sprintf("\n--- ssGSEA per %s cluster ---\n", cluster_backend_label))
        for (method in SSGSEA_CHOIR_METHODS) {
          cat(sprintf("  method: %s\n", method))
          selected_t2g_choir <- switch(
            method, hallmark = hallmark_t2g, go_bp = go_bp_t2g, go_mf = go_mf_t2g,
            go_cc = go_cc_t2g, kegg = kegg_t2g,
            custom = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT))
              tryCatch(read.gmt(SSGSEA_CUSTOM_GMT) %>% mutate(gene = toupper(gene)), error = function(e) NULL)
            else NULL, NULL
          )
          if (is.null(selected_t2g_choir) || nrow(selected_t2g_choir) == 0) {
            cat(sprintf("  [WARN] '%s': no gene sets\n", method)); next
          }
          gs_choir <- term2gene_to_list(selected_t2g_choir)
          gs_choir <- map_gene_sets_to_features(gs_choir, rownames(obj))
          gs_choir <- filter_gs_size(gs_choir)
          if (length(gs_choir) == 0) { cat(sprintf("  [WARN] '%s': no valid gene sets\n", method)); next }
          needed_genes_choir <- unique(unlist(gs_choir, use.names = FALSE))
          avg_choir <- tryCatch(
            AverageExpression(obj, assays = "RNA", slot = "data",
                              group.by = cluster_col, features = needed_genes_choir, verbose = FALSE)[["RNA"]],
            error = function(e) { cat(sprintf("  [ERROR] AverageExpression: %s\n", e$message)); NULL }
          )
          if (is.null(avg_choir)) next
          bp_choir <- BiocParallel::SnowParam(workers = N_CORES, type = "SOCK", progressbar = FALSE)
            ss_choir <- tryCatch({
              if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
                GSVA::gsva(
                  GSVA::ssgseaParam(
                    exprData = as.matrix(avg_choir),
                    geneSets = gs_choir,
                    alpha = 0.25,
                    normalize = TRUE,
                    minSize = 10,
                    maxSize = 500
                  ),
                  BPPARAM = bp_choir,
                  verbose = FALSE
                )
              } else {
                GSVA::gsva(as.matrix(avg_choir), gs_choir, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
              }
            }, error = function(e) { cat(sprintf("  [ERROR] ssGSEA: %s\n", e$message)); NULL },
            finally = {
              try(BiocParallel::bpstop(bp_choir), silent = TRUE)
            })
          if (is.null(ss_choir)) next
          sz_choir <- t(scale(t(ss_choir)))
          top_choir_df <- extract_directional_ssgsea_top_rows(
            ss_choir,
            sz_choir,
            top_n = SSGSEA_CHOIR_N_TOP,
            group_field = "cluster",
            method = method
          )
          ssgsea_choir_all[[method]] <- list(scores = ss_choir, z_scores = sz_choir, top_df = top_choir_df)
          saveRDS(ss_choir, file.path(cluster_dir, sprintf("ssgsea_scores_%s.rds", method)))
          saveRDS(sz_choir, file.path(cluster_dir, sprintf("ssgsea_z_%s.rds", method)))
          fwrite(top_choir_df, file.path(cluster_dir, sprintf("ssgsea_top_pathways_%s.csv", method)))
          maz <- rowMeans(abs(sz_choir), na.rm = TRUE)
          top_pc <- names(sort(maz, decreasing = TRUE))[seq_len(min(SSGSEA_CHOIR_N_HEATMAP, length(maz)))]
          hm_c <- sz_choir[top_pc, , drop = FALSE]
          hm_c[is.nan(hm_c) | is.infinite(hm_c)] <- 0
          if (nrow(hm_c) >= 3 && ncol(hm_c) >= 2) {
            tryCatch({
              hp <- file.path(fig_dir, sprintf("ssgsea_%s_heatmap_%s", cluster_output_prefix, method))
              ht_c <- pheatmap::pheatmap(
                hm_c,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
                breaks = seq(-3, 3, length.out = 101),
                main = sprintf("%s %s ssGSEA Z-score (%s) — top %d pathways", LINEAGE_DISPLAY, cluster_backend_label, method, nrow(hm_c)),
                fontsize_row = 7,
                fontsize_col = 9,
                cellwidth = 22,
                cellheight = 10,
                filename = paste0(hp, ".pdf"),
                width = max(10, ncol(hm_c) * 1.8),
                height = max(8, nrow(hm_c) * 0.35 + 3)
              )
              grDevices::png(paste0(hp, ".png"), width = max(10, ncol(hm_c) * 1.8),
                             height = max(8, nrow(hm_c) * 0.35 + 3), units = "in", res = 300)
              grid::grid.newpage(); grid::grid.draw(ht_c$gtable); grDevices::dev.off()
              cat(sprintf("  [OK] Heatmap: ssgsea_%s_heatmap_%s\n", cluster_output_prefix, method))
            }, error = function(e) cat(sprintf("  [WARN] Heatmap %s: %s\n", method, e$message)))
          }
          cat(sprintf("  [OK] %s: %d pathways x %d clusters\n", method, nrow(ss_choir), ncol(ss_choir)))
        }
        choir_results_all[["ssgsea"]] <- ssgsea_choir_all
        cat(sprintf("[OK] %s cluster ssGSEA complete\n", cluster_backend_label))
      }

      # 8.6b OFA: one-vs-rest + enrichment + [LLM-5] per-cluster LLM
      if (RUN_OFA) {
        cat(sprintf("\n--- OFA (one-vs-rest) per %s cluster ---\n", cluster_backend_label))
        for (cl in choir_clusters) {
          cl_name <- as.character(cl)
          focal_cells <- colnames(obj)[as.character(obj@meta.data[[cluster_col]]) == cl_name]
          rest_cells <- colnames(obj)[as.character(obj@meta.data[[cluster_col]]) != cl_name]
          choir_cluster_l2 <- resolve_dominant_label(obj@meta.data[focal_cells, CELLTYPE_L2_COL])
          cat(sprintf("  OFA: %s cluster %s vs rest\n", cluster_backend_label, cl_name))

          ofa_dir <- file.path(cluster_dir, sprintf("ofa_%s", safe_name(cl_name)))
          ofa_res <- run_ofa_marker_enrichment(
            obj = obj,
            focal_cells = focal_cells,
            rest_cells = rest_cells,
            focal_label = cl_name,
            rest_label = "rest",
            out_dir = ofa_dir,
            plot_title = sprintf("OFA: %s %s cluster %s vs rest", LINEAGE_DISPLAY, cluster_backend_label, cl_name),
            group_col_name = "ofa_group",
            min_cells_focal = OFA_MIN_CELLS_FOCAL,
            min_cells_rest = OFA_MIN_CELLS_REST,
            padj_thr = OFA_PADJ_THR,
            lfc_thr = OFA_LFC_THR,
            top_n = OFA_TOP_N,
            max_cells_per_ident = OFA_MAX_CELLS_PER_IDENT
          )
          if (isTRUE(ofa_res$skipped)) {
            cat(sprintf("    [SKIP] %s\n", ofa_res$reason))
            next
          }
          cat(sprintf("    %d sig markers (padj<%.2f |log2FC|>%.2f)\n", ofa_res$n_sig, OFA_PADJ_THR, OFA_LFC_THR))

          choir_ofa_llm <- list()
          if (ENABLE_LLM) {
            choir_llm_rec <- tryCatch(
              run_choir_cluster_llm(
                cluster_id = cl_name,
                de_df = ofa_res$de_table,
                ofa_enrich = ofa_res$enrich,
                ssgsea_choir_all = choir_results_all[["ssgsea"]],
                choir_cluster_l2 = choir_cluster_l2,
                choir_annotation_tables = choir_annotation_tables
              ),
              error = function(e) {
                cat(sprintf("    [WARN] %s integrated LLM: %s\n", cluster_backend_label, e$message))
                NULL
              }
            )
            if (!is.null(choir_llm_rec)) {
              if (!identical(cluster_backend, "CHOIR")) {
                choir_llm_rec$celltype_level <- cluster_backend_label
                choir_llm_rec$celltype_label <- paste0(cluster_backend_label, "_", cl_name)
              }
              choir_ofa_llm[["integrated"]] <- choir_llm_rec
              saveRDS(choir_llm_rec, file.path(ofa_dir, "llm_integrated_structured.rds"))
              if (nzchar(choir_llm_rec$raw_text)) {
                writeLines(choir_llm_rec$raw_text, file.path(ofa_dir, "llm_integrated_raw.txt"))
              }
              cat(sprintf("    [OK] %s cluster %s integrated LLM status=%s\n", cluster_backend_label, cl_name, choir_llm_rec$status))
              Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
            }
          }

          choir_ofa_all[[cl_name]] <- list(
            cluster_backend = cluster_backend_label,
            de_table = ofa_res$de_table,
            enrich = ofa_res$enrich,
            llm = choir_ofa_llm,
            n_focal = ofa_res$n_focal,
            n_rest = ofa_res$n_rest,
            n_sig = ofa_res$n_sig
          )
        }
        choir_results_all[["ofa"]] <- choir_ofa_all
        cat(sprintf("[OK] %s OFA complete: %d clusters processed\n", cluster_backend_label, length(choir_ofa_all)))
      }

      saveRDS(choir_results_all, file.path(cluster_dir, cluster_results_rds_name))
      saveRDS(choir_ofa_all, file.path(cluster_dir, cluster_ofa_rds_name))
      cat(sprintf("[OK] %s results saved: %s\n", cluster_backend_label, cluster_dir))
    }
  }
}

# ==============================================================================
# 8.7 L3 OFA: vs other cell types + same L2 siblings + same L3 across tissues
# ============================================================================== 

if (RUN_L3_OFA && (RUN_L3_OFA_VS_REST || RUN_L3_OFA_SAME_L2 || RUN_L3_OFA_INTER_TISSUE)) {
  cat("\n=== L3 OFA (cell-level one-vs-rest) ===\n")
  meta_l3_ofa <- obj@meta.data
  l3_ofa_dir <- file.path(rpt_dir, "l3_ofa")
  l3_ofa_vs_rest_dir <- file.path(l3_ofa_dir, "vs_other_cell_types")
  l3_ofa_same_l2_dir <- file.path(l3_ofa_dir, "vs_same_l2_other_l3")
  l3_ofa_inter_tissue_dir <- file.path(l3_ofa_dir, "same_l3_other_tissues")
  dir.create(l3_ofa_vs_rest_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(l3_ofa_same_l2_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(l3_ofa_inter_tissue_dir, recursive = TRUE, showWarnings = FALSE)

  l3_vs_rows <- list()
  l3_same_l2_rows <- list()
  l3_tissue_rows <- list()

  for (l3_name in l3_types) {
    focal_cells <- rownames(meta_l3_ofa)[as.character(meta_l3_ofa[[CELLTYPE_L3_COL]]) == as.character(l3_name)]
    if (length(focal_cells) == 0) next
    focal_l2 <- resolve_dominant_label(meta_l3_ofa[focal_cells, CELLTYPE_L2_COL, drop = TRUE])

    if (RUN_L3_OFA_VS_REST) {
      cat(sprintf("  OFA: L3 %s vs other cell types\n", l3_name))
      rest_cells <- setdiff(colnames(obj), focal_cells)
      ofa_dir <- file.path(l3_ofa_vs_rest_dir, safe_name(l3_name))
      ofa_res <- run_ofa_marker_enrichment(
        obj = obj,
        focal_cells = focal_cells,
        rest_cells = rest_cells,
        focal_label = as.character(l3_name),
        rest_label = "other_cell_types",
        out_dir = ofa_dir,
        plot_title = sprintf("OFA: %s L3 %s vs other cell types", LINEAGE_DISPLAY, l3_name),
        group_col_name = "l3_ofa_group",
        min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
        min_cells_rest = L3_OFA_MIN_CELLS_REST,
        padj_thr = L3_OFA_PADJ_THR,
        lfc_thr = L3_OFA_LFC_THR,
        top_n = L3_OFA_TOP_N,
        max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
      )
      if (isTRUE(ofa_res$skipped)) {
        cat(sprintf("    [SKIP] %s\n", ofa_res$reason))
      } else {
        l3_ofa_vs_rest_all[[as.character(l3_name)]] <- c(
          ofa_res,
          list(comparison = "vs_other_cell_types", celltype_l2 = focal_l2)
        )
        l3_vs_rows[[length(l3_vs_rows) + 1L]] <- data.frame(
          family = "L3_vs_other_cell_types",
          celltype_label = as.character(l3_name),
          celltype_l2 = ifelse(is.na(focal_l2), "", as.character(focal_l2)),
          comparison = "vs_other_cell_types",
          n_focal = ofa_res$n_focal,
          n_rest = ofa_res$n_rest,
          n_sig = ofa_res$n_sig,
          bubble_png = if (!is.null(ofa_res$bubble_paths$png)) ofa_res$bubble_paths$png else NA_character_,
          bubble_pdf = if (!is.null(ofa_res$bubble_paths$pdf)) ofa_res$bubble_paths$pdf else NA_character_,
          output_dir = ofa_dir,
          stringsAsFactors = FALSE
        )
        cat(sprintf("    [OK] %d sig markers\n", ofa_res$n_sig))
      }
    }

    if (RUN_L3_OFA_SAME_L2) {
      if (is.na(focal_l2) || !nzchar(trimws(as.character(focal_l2)))) {
        cat(sprintf("  [SKIP] L3 %s same-L2 siblings: no dominant L2\n", l3_name))
      } else {
        same_l2_cells <- rownames(meta_l3_ofa)[
          as.character(meta_l3_ofa[[CELLTYPE_L2_COL]]) == as.character(focal_l2) &
            as.character(meta_l3_ofa[[CELLTYPE_L3_COL]]) != as.character(l3_name)
        ]
        if (length(same_l2_cells) == 0) {
          cat(sprintf("  [SKIP] L3 %s same-L2 siblings: no sibling L3 cells within %s\n", l3_name, focal_l2))
        } else {
          cat(sprintf("  OFA: L3 %s vs same-L2 siblings (%s)\n", l3_name, focal_l2))
          ofa_dir <- file.path(l3_ofa_same_l2_dir, safe_name(l3_name))
          ofa_res <- run_ofa_marker_enrichment(
            obj = obj,
            focal_cells = focal_cells,
            rest_cells = same_l2_cells,
            focal_label = as.character(l3_name),
            rest_label = paste0("other_", safe_name(focal_l2), "_L3"),
            out_dir = ofa_dir,
            plot_title = sprintf("OFA: %s L3 %s vs same-L2 siblings (%s)", LINEAGE_DISPLAY, l3_name, focal_l2),
            group_col_name = "l3_same_l2_ofa_group",
            min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
            min_cells_rest = L3_OFA_MIN_CELLS_REST,
            padj_thr = L3_OFA_PADJ_THR,
            lfc_thr = L3_OFA_LFC_THR,
            top_n = L3_OFA_TOP_N,
            max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
          )
          if (isTRUE(ofa_res$skipped)) {
            cat(sprintf("    [SKIP] %s\n", ofa_res$reason))
          } else {
            l3_ofa_same_l2_all[[as.character(l3_name)]] <- c(
              ofa_res,
              list(comparison = "vs_same_l2_other_l3", celltype_l2 = focal_l2)
            )
            l3_same_l2_rows[[length(l3_same_l2_rows) + 1L]] <- data.frame(
              family = "L3_vs_same_L2_other_L3",
              celltype_label = as.character(l3_name),
              celltype_l2 = as.character(focal_l2),
              comparison = "vs_same_l2_other_l3",
              n_focal = ofa_res$n_focal,
              n_rest = ofa_res$n_rest,
              n_sig = ofa_res$n_sig,
              bubble_png = if (!is.null(ofa_res$bubble_paths$png)) ofa_res$bubble_paths$png else NA_character_,
              bubble_pdf = if (!is.null(ofa_res$bubble_paths$pdf)) ofa_res$bubble_paths$pdf else NA_character_,
              output_dir = ofa_dir,
              stringsAsFactors = FALSE
            )
            cat(sprintf("    [OK] %d sig markers\n", ofa_res$n_sig))
          }
        }
      }
    }

    if (RUN_L3_OFA_INTER_TISSUE) {
      l3_meta <- meta_l3_ofa[focal_cells, , drop = FALSE]
      tissue_values <- sort(unique(na.omit(as.character(l3_meta[[TISSUE_COL]]))))
      if (length(tissue_values) < 2) next
      for (tissue_name in tissue_values) {
        focal_cells <- rownames(l3_meta)[as.character(l3_meta[[TISSUE_COL]]) == tissue_name]
        rest_cells <- rownames(l3_meta)[as.character(l3_meta[[TISSUE_COL]]) != tissue_name]
        comparison_label <- paste0(tissue_name, "_vs_other_tissues")
        cat(sprintf("  OFA: L3 %s | %s vs other tissues\n", l3_name, tissue_name))
        ofa_dir <- file.path(l3_ofa_inter_tissue_dir, safe_name(l3_name), safe_name(comparison_label))
        ofa_res <- run_ofa_marker_enrichment(
          obj = obj,
          focal_cells = focal_cells,
          rest_cells = rest_cells,
          focal_label = as.character(tissue_name),
          rest_label = "other_tissues",
          out_dir = ofa_dir,
          plot_title = sprintf("OFA: %s L3 %s | %s vs same L3 in other tissues", LINEAGE_DISPLAY, l3_name, tissue_name),
          group_col_name = "l3_tissue_ofa_group",
          min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
          min_cells_rest = L3_OFA_MIN_CELLS_REST,
          padj_thr = L3_OFA_PADJ_THR,
          lfc_thr = L3_OFA_LFC_THR,
          top_n = L3_OFA_TOP_N,
          max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
        )
        if (isTRUE(ofa_res$skipped)) {
          cat(sprintf("    [SKIP] %s\n", ofa_res$reason))
          next
        }
        if (is.null(l3_ofa_inter_tissue_all[[as.character(l3_name)]])) {
          l3_ofa_inter_tissue_all[[as.character(l3_name)]] <- list()
        }
        l3_ofa_inter_tissue_all[[as.character(l3_name)]][[comparison_label]] <- c(
          ofa_res,
          list(comparison = comparison_label, focal_tissue = tissue_name)
        )
        l3_tissue_rows[[length(l3_tissue_rows) + 1L]] <- data.frame(
          family = "L3_same_celltype_across_tissues",
          celltype_label = as.character(l3_name),
          celltype_l2 = ifelse(is.na(focal_l2), "", as.character(focal_l2)),
          comparison = comparison_label,
          n_focal = ofa_res$n_focal,
          n_rest = ofa_res$n_rest,
          n_sig = ofa_res$n_sig,
          bubble_png = if (!is.null(ofa_res$bubble_paths$png)) ofa_res$bubble_paths$png else NA_character_,
          bubble_pdf = if (!is.null(ofa_res$bubble_paths$pdf)) ofa_res$bubble_paths$pdf else NA_character_,
          output_dir = ofa_dir,
          stringsAsFactors = FALSE
        )
        cat(sprintf("    [OK] %d sig markers\n", ofa_res$n_sig))
      }
    }
  }

  l3_ofa_vs_rest_summary_df <- if (length(l3_vs_rows) > 0) dplyr::bind_rows(l3_vs_rows) else empty_l3_ofa_summary_df()
  l3_ofa_same_l2_summary_df <- if (length(l3_same_l2_rows) > 0) dplyr::bind_rows(l3_same_l2_rows) else empty_l3_ofa_summary_df()
  l3_ofa_inter_tissue_summary_df <- if (length(l3_tissue_rows) > 0) dplyr::bind_rows(l3_tissue_rows) else empty_l3_ofa_summary_df()

  saveRDS(l3_ofa_vs_rest_all, file.path(l3_ofa_dir, "l3_ofa_vs_rest_all.rds"))
  saveRDS(l3_ofa_same_l2_all, file.path(l3_ofa_dir, "l3_ofa_same_l2_all.rds"))
  saveRDS(l3_ofa_inter_tissue_all, file.path(l3_ofa_dir, "l3_ofa_inter_tissue_all.rds"))
  fwrite(l3_ofa_vs_rest_summary_df, file.path(l3_ofa_dir, "l3_ofa_vs_rest_summary.tsv"), sep = "\t")
  fwrite(l3_ofa_same_l2_summary_df, file.path(l3_ofa_dir, "l3_ofa_same_l2_summary.tsv"), sep = "\t")
  fwrite(l3_ofa_inter_tissue_summary_df, file.path(l3_ofa_dir, "l3_ofa_inter_tissue_summary.tsv"), sep = "\t")
  cat(sprintf("[OK] L3 OFA complete: vs-rest=%d | same-L2=%d | inter-tissue=%d\n",
              ifelse(nrow(l3_ofa_vs_rest_summary_df) > 0, nrow(l3_ofa_vs_rest_summary_df), 0),
              ifelse(nrow(l3_ofa_same_l2_summary_df) > 0, nrow(l3_ofa_same_l2_summary_df), 0),
              ifelse(nrow(l3_ofa_inter_tissue_summary_df) > 0, nrow(l3_ofa_inter_tissue_summary_df), 0)))
}

if (isTRUE(EXIT_AFTER_L3_OFA)) {
  cat("[INFO] EXIT_AFTER_L3_OFA=TRUE; stopping pipeline immediately after L3 OFA stage.\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

# ==============================================================================
# 9. Exploratory Wilcoxon
# ==============================================================================

cat("\n=== Exploratory Cell-level Wilcoxon ===\n")
cat("NOTE: P-values are inflated (pseudoreplication). For marker discovery ONLY.\n\n")

wilcox_all <- list()
wilcox_l3_all <- list()
if (isTRUE(RUN_WILCOX)) {
  for (ct_l2 in l2_types) {
    wx <- run_wilcox_exploratory(obj, ct_l2, group_col = CELLTYPE_L2_COL)
    if (!is.null(wx) && length(wx) > 0) {
      wilcox_all[[ct_l2]] <- wx
      for (comp_name in names(wx)) {
        wx_out <- file.path(wx_dir, safe_name(ct_l2)); dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
        fwrite(wx[[comp_name]], file.path(wx_out, paste0("exploratory_", safe_name(comp_name), "_wilcox.csv")))
      }
      cat(sprintf("  [OK] L2 %s: %d comparisons\n", ct_l2, length(wx)))
    }
  }
  for (ct_l3 in l3_types) {
    wx <- run_wilcox_exploratory(obj, ct_l3, group_col = CELLTYPE_L3_COL)
    if (!is.null(wx) && length(wx) > 0) {
      wilcox_l3_all[[ct_l3]] <- wx
      for (comp_name in names(wx)) {
        wx_out <- file.path(wx_l3_dir, safe_name(ct_l3)); dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
        fwrite(wx[[comp_name]], file.path(wx_out, paste0("exploratory_", safe_name(comp_name), "_wilcox.csv")))
      }
      cat(sprintf("  [OK] L3 %s: %d comparisons\n", ct_l3, length(wx)))
    }
  }
} else {
  cat("[INFO] Exploratory Wilcoxon skipped; loading cached artifacts when available\n")
  wilcox_all <- load_cached_report_rds("wilcox_exploratory_all.rds", default = list(), label = "L2 exploratory wilcox")
  wilcox_l3_all <- load_cached_report_rds("wilcox_exploratory_L3_all.rds", default = list(), label = "L3 exploratory wilcox")
}

# ==============================================================================
# 10. Save Summary RDS
# ==============================================================================

cat("\n=== Saving Summary ===\n")
saveRDS(pb_de_all,               file.path(rpt_dir, "pseudobulk_de_all.rds"))
saveRDS(pb_de_l3_all,            file.path(rpt_dir, "pseudobulk_de_L3_all.rds"))
saveRDS(wilcox_all,              file.path(rpt_dir, "wilcox_exploratory_all.rds"))
saveRDS(wilcox_l3_all,           file.path(rpt_dir, "wilcox_exploratory_L3_all.rds"))
saveRDS(enrich_all,              file.path(rpt_dir, "enrichment_all.rds"))
saveRDS(enrich_l3_all,           file.path(rpt_dir, "enrichment_L3_all.rds"))
saveRDS(agent_all,               file.path(rpt_dir, "interpret_agent_all.rds"))
saveRDS(agent_l3_all,            file.path(rpt_dir, "interpret_agent_L3_all.rds"))
saveRDS(agent_structured_all,    file.path(rpt_dir, "interpret_agent_structured_all.rds"))
saveRDS(agent_structured_l3_all, file.path(rpt_dir, "interpret_agent_structured_L3_all.rds"))
saveRDS(l3_ofa_vs_rest_all,      file.path(rpt_dir, "l3_ofa_vs_rest_all.rds"))
saveRDS(l3_ofa_same_l2_all,      file.path(rpt_dir, "l3_ofa_same_l2_all.rds"))
saveRDS(l3_ofa_inter_tissue_all, file.path(rpt_dir, "l3_ofa_inter_tissue_all.rds"))
if (is.data.frame(l3_ofa_vs_rest_summary_df)) {
  fwrite(l3_ofa_vs_rest_summary_df, file.path(rpt_dir, "l3_ofa_vs_rest_summary.tsv"), sep = "\t")
}
if (is.data.frame(l3_ofa_same_l2_summary_df)) {
  fwrite(l3_ofa_same_l2_summary_df, file.path(rpt_dir, "l3_ofa_same_l2_summary.tsv"), sep = "\t")
}
if (is.data.frame(l3_ofa_inter_tissue_summary_df)) {
  fwrite(l3_ofa_inter_tissue_summary_df, file.path(rpt_dir, "l3_ofa_inter_tissue_summary.tsv"), sep = "\t")
}
fwrite(l3_l2_mapping_summary,    file.path(rpt_dir, "L3_to_L2_mapping_summary.csv"))
saveRDS(if (isTRUE(USE_EXISTING_L2)) l3_l2_mapping_summary else L3_TO_L2_REMAP,
  file.path(rpt_dir, "L3_to_L2_remap.rds"))

agent_structured_df_l2 <- flatten_interpretation_records(agent_structured_all)
agent_structured_df_l3 <- flatten_interpretation_records(agent_structured_l3_all)
agent_structured_df    <- dplyr::bind_rows(agent_structured_df_l2, agent_structured_df_l3)
if (nrow(agent_structured_df) > 0) {
  fwrite(agent_structured_df, file.path(rpt_dir, "interpret_agent_structured.tsv"), sep = "\t")
}
if (nrow(agent_structured_df_l3) > 0) {
  fwrite(agent_structured_df_l3, file.path(rpt_dir, "interpret_agent_structured_L3.tsv"), sep = "\t")
}
write_interpretation_markdown(
  agent_structured_df, file.path(OUTPUT_DIR, "LLM_INTERPRETATION.md"),
  "# LLM Interpretation Summary (Normal Tissue Comparison)", include_raw = FALSE
)
write_interpretation_markdown(
  agent_structured_df, file.path(rpt_dir, "LLM_INTERPRETATION_FOR_LLM.md"),
  "# LLM Interpretation Structured Input", include_raw = TRUE
)

# --- [LLM-4] ssGSEA LLM structured records ---
ssgsea_llm_df_l2 <- collect_ssgsea_llm_records(ssgsea_results_all)
ssgsea_llm_df_l3 <- collect_ssgsea_llm_records(ssgsea_results_l3_all)
ssgsea_llm_df    <- dplyr::bind_rows(ssgsea_llm_df_l2, ssgsea_llm_df_l3)
saveRDS(ssgsea_llm_df_l2, file.path(rpt_dir, "ssgsea_llm_structured_L2_all.rds"))
saveRDS(ssgsea_llm_df_l3, file.path(rpt_dir, "ssgsea_llm_structured_L3_all.rds"))
saveRDS(ssgsea_llm_df,    file.path(rpt_dir, "ssgsea_llm_structured_all.rds"))
if (nrow(ssgsea_llm_df) > 0) {
  fwrite(ssgsea_llm_df, file.path(rpt_dir, "ssgsea_llm_structured.tsv"), sep = "\t")
  write_interpretation_markdown(
    ssgsea_llm_df, file.path(OUTPUT_DIR, "LLM_SSGSEA_INTERPRETATION.md"),
    "# LLM ssGSEA Interpretation (Tissue x Cell Type Groups)", include_raw = FALSE
  )
  cat(sprintf("[OK] ssGSEA LLM records: %d rows -> ssgsea_llm_structured.tsv\n", nrow(ssgsea_llm_df)))
}

# --- [LLM-5] Cluster LLM structured records ---
choir_llm_df <- collect_choir_llm_records(choir_ofa_all)
saveRDS(choir_llm_df, file.path(rpt_dir, cluster_llm_rds_filename))
if (nrow(choir_llm_df) > 0) {
  fwrite(choir_llm_df, file.path(rpt_dir, cluster_llm_tsv_filename), sep = "\t")
  write_interpretation_markdown(
    choir_llm_df,
    cluster_llm_markdown_path,
    sprintf("# LLM %s Cluster Interpretation (One-vs-Rest)", cluster_backend_label),
    include_raw = FALSE
  )
  cat(sprintf("[OK] %s LLM records: %d rows -> %s\n", cluster_backend_label, nrow(choir_llm_df), cluster_llm_tsv_filename))
}

# --- Additional full-record LLM screening for biological discoveries / outliers ---
ssgsea_discovery_screen_df <- data.frame()
choir_discovery_screen_df <- data.frame()
ofa_discovery_screen_df <- data.frame()
if (ENABLE_LLM) {
  if (nrow(ssgsea_llm_df) > 0) {
    ssgsea_screen_input <- ssgsea_llm_df %>%
      dplyr::mutate(
        record_id = sprintf("ssgsea_%04d", dplyr::row_number()),
        record_label = paste(celltype_level, comparison, direction, sep = " | "),
        annotation_label = ifelse(nzchar(celltype_label), celltype_label, celltype_l2),
        primary_text = ifelse(nzchar(cell_type_judgment), cell_type_judgment, overview),
        supporting_text = paste(
          ifelse(nzchar(annotated_l3_correspondence), annotated_l3_correspondence, ""),
          ifelse(nzchar(discovery_assessment), discovery_assessment, ""),
          ifelse(nzchar(outlier_assessment), outlier_assessment, ""),
          ifelse(nzchar(integrated_diagnostic_comment), integrated_diagnostic_comment, evidence),
          sep = "\n"
        )
      ) %>%
      dplyr::select(record_id, record_label, annotation_label, primary_text, supporting_text)
    ssgsea_discovery_screen_df <- run_and_write_discovery_screen(
      records_df = ssgsea_screen_input,
      family_label = "ssGSEA grouped interpretations",
      output_dir = OUTPUT_DIR,
      prefix = "llm_ssgsea",
      title = "# LLM ssGSEA Discovery / Outlier Review",
      batch_size = 12L
    )
    cat(sprintf("[OK] ssGSEA discovery screen: %d rows\n", nrow(ssgsea_discovery_screen_df)))
  }
  if (nrow(choir_llm_df) > 0) {
    choir_screen_input <- choir_llm_df %>%
      dplyr::mutate(
        record_id = sprintf("%s_%04d", cluster_output_prefix, dplyr::row_number()),
        record_label = paste(celltype_level, comparison, direction, sep = " | "),
        annotation_label = ifelse(nzchar(annotated_l3_dominant), annotated_l3_dominant, celltype_l2),
        primary_text = ifelse(nzchar(cell_type_judgment), cell_type_judgment, overview),
        supporting_text = paste(
          ifelse(nzchar(annotated_l3_correspondence), annotated_l3_correspondence, ""),
          ifelse(nzchar(discovery_assessment), discovery_assessment, ""),
          ifelse(nzchar(outlier_assessment), outlier_assessment, ""),
          ifelse(nzchar(integrated_diagnostic_comment), integrated_diagnostic_comment, evidence),
          sep = "\n"
        )
      ) %>%
      dplyr::select(record_id, record_label, annotation_label, primary_text, supporting_text)
    choir_discovery_screen_df <- run_and_write_discovery_screen(
      records_df = choir_screen_input,
      family_label = sprintf("%s cluster interpretations", cluster_backend_label),
      output_dir = OUTPUT_DIR,
      prefix = paste0("llm_", cluster_output_prefix),
      title = sprintf("# LLM %s Discovery / Outlier Review", cluster_backend_label),
      batch_size = 10L
    )
    cat(sprintf("[OK] %s discovery screen: %d rows\n", cluster_backend_label, nrow(choir_discovery_screen_df)))
  }
  if (length(choir_ofa_all) > 0) {
    ofa_screen_input <- build_choir_ofa_screening_table(choir_ofa_all)
    if (nrow(ofa_screen_input) > 0) {
      ofa_discovery_screen_df <- run_and_write_discovery_screen(
        records_df = ofa_screen_input,
        family_label = sprintf("%s OFA one-vs-rest summaries", cluster_backend_label),
        output_dir = OUTPUT_DIR,
        prefix = "llm_ofa",
        title = sprintf("# LLM %s OFA Discovery / Outlier Review", cluster_backend_label),
        batch_size = 10L
      )
      cat(sprintf("[OK] OFA discovery screen: %d rows\n", nrow(ofa_discovery_screen_df)))
    }
  }
}

# ==============================================================================
# 11. Generate REPORT.md
# ==============================================================================

cat("\n=== Generating REPORT.md ===\n")
md <- character()
add <- function(...) md <<- c(md, paste0(...))

add(REPORT_TITLE); add("")
add("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M")); add("")
add("**Pipeline:** ", PIPELINE_SUBTITLE); add("")
add("**Note:** Cross-site anatomical comparison of NORMAL tissues, NOT disease vs healthy."); add("")
add("---"); add("")

add("## 1. Data Overview"); add("")
add(sprintf("- **Input:** `%s`", basename(H5AD_PATH)))
add(sprintf("- **Total cells:** %s", format(ncol(obj), big.mark = ",")))
add(sprintf("- **Tissues:** %s", paste(tissues, collapse = ", ")))
add(sprintf("- **L2 subtypes (%d):** %s", length(l2_types), paste(l2_types, collapse = ", ")))
add(sprintf("- **L3 subtypes (%d):** %s", length(l3_types), paste(l3_types, collapse = ", ")))
add(sprintf("- **L3 source column:** `%s` -> standardized `%s`", L3_SOURCE_COL, CELLTYPE_L3_COL))
if (isTRUE(USE_EXISTING_L2)) {
  add(sprintf("- **L2 source column:** `%s` (reused from input metadata)", L2_SOURCE_COL))
  add(""); add("### Observed L3 / L2 Mapping"); add("")
  add(sprintf("| %s | Existing `%s` | Cells |", L3_TO_L2_TABLE_HEADER_LEFT, L2_SOURCE_COL)); add("|---|---|---|")
  for (i in seq_len(nrow(l3_l2_mapping_summary))) {
    add(sprintf("| %s | %s | %d |",
                l3_l2_mapping_summary$L3[i],
                l3_l2_mapping_summary$L2[i],
                l3_l2_mapping_summary$Freq[i]))
  }
} else {
  add(""); add("### L3 -> L2 Remapping"); add("")
  add(sprintf("| %s | L2 (merged) |", L3_TO_L2_TABLE_HEADER_LEFT)); add("|---|---|")
  for (i in seq_along(L3_TO_L2_REMAP)) add(sprintf("| %s | %s |", names(L3_TO_L2_REMAP)[i], L3_TO_L2_REMAP[i]))
}
add("")

add("## 2. Visualization"); add("")
add("![UMAP tissue](figures/umap_tissue.png)"); add("")
add("![UMAP L2](figures/umap_celltype_L2.png)"); add("")
add("![UMAP L3](figures/umap_celltype_L3.png)"); add("")
add("![UMAP split](figures/umap_L2_split_tissue.png)"); add("")
add("![Dotplot](figures/dotplot_markers.png)"); add("")
add("![Composition](figures/composition_tissue_L2.png)"); add("")
add("![Sample composition](figures/composition_sample_level.png)"); add("")
add("![Heatmap](figures/heatmap_top_markers.png)"); add("")
if (!is.null(marker_panel_viz) && is.data.frame(marker_panel_viz$summary) && nrow(marker_panel_viz$summary) > 0 &&
    exists("tc_marker_panel_report_lines", mode = "function")) {
  for (line in tc_marker_panel_report_lines(marker_panel_viz$summary, figure_dir_rel = sprintf("figures/%s", MARKER_PANEL_FIG_SUBDIR))) {
    add(line)
  }
}

add("## 3. Pseudobulk DESeq2 (Primary Inference)"); add("")
add(sprintf("padj < %s, |log2FC| > %s", PADJ_THR, LFC_THR)); add("")
add("### 3.1 L2"); add("")
add("| L2 Subtype | Comparison | Samples (ref/case) | Up | Down | Total |"); add("|---|---|---|---|---|---|")
for (ct in names(pb_de_all)) for (comp in names(pb_de_all[[ct]])) {
  r <- pb_de_all[[ct]][[comp]]; if (is.null(r)) next
  add(sprintf("| %s | %s | %d / %d | %d | %d | %d |", ct, comp, r$n_samples_1, r$n_samples_2, r$n_up, r$n_down, r$n_up + r$n_down))
}
add(""); add("### 3.2 L3"); add("")
add("| L3 Subtype | Comparison | Samples (ref/case) | Up | Down | Total |"); add("|---|---|---|---|---|---|")
for (ct in names(pb_de_l3_all)) for (comp in names(pb_de_l3_all[[ct]])) {
  r <- pb_de_l3_all[[ct]][[comp]]; if (is.null(r)) next
  add(sprintf("| %s | %s | %d / %d | %d | %d | %d |", ct, comp, r$n_samples_1, r$n_samples_2, r$n_up, r$n_down, r$n_up + r$n_down))
}
add("")

add("## 4. Multi-Database Enrichment"); add("")
add(sprintf("Databases: GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, %s", CUSTOM_DB_LABEL)); add("")
for (ct in names(enrich_all)) for (comp in names(enrich_all[[ct]])) for (dir in names(enrich_all[[ct]][[comp]])) {
  enr_l <- enrich_all[[ct]][[comp]][[dir]]; if (is.null(enr_l) || length(enr_l) == 0) next
  add(sprintf("### %s | %s | %s", ct, comp, dir)); add("")
  for (db in names(enr_l)) {
    er <- enr_l[[db]]; if (is.null(er) || nrow(as.data.frame(er)) == 0) next
    top5 <- head(as.data.frame(er), 5)
    add(sprintf("**%s (top 5):**", db)); add("")
    add("| Term | p.adjust | Count |"); add("|---|---|---|")
    for (j in seq_len(nrow(top5))) add(sprintf("| %s | %.2e | %s |", top5$Description[j], top5$p.adjust[j], top5$Count[j]))
    add("")
  }
}

add("## 5. Grouped ssGSEA (Average-Expression)"); add("")
if (length(ssgsea_results_all) == 0 && length(ssgsea_results_l3_all) == 0) {
  add("ssGSEA not run or no methods succeeded."); add("")
} else {
  ssgsea_methods_report <- unique(c(get_ssgsea_method_names(ssgsea_results_all), get_ssgsea_method_names(ssgsea_results_l3_all)))
  add(sprintf("Methods: %s", paste(ssgsea_methods_report, collapse = ", "))); add("")
  add("### 5.1 Tissue x L2"); add("")
  add("| Method | Pathways | Groups | Output |"); add("|---|---|---|---|")
  for (mn in get_ssgsea_method_names(ssgsea_results_all)) {
    add(sprintf("| %s | %d | %d | `reports/ssgsea_scores_%s.rds` |",
                mn, nrow(ssgsea_results_all[[mn]]$scores), ncol(ssgsea_results_all[[mn]]$scores), mn))
  }
  add("")
  for (mn in get_ssgsea_method_names(ssgsea_results_all)) { add(sprintf("![ssGSEA heatmap %s](figures/ssgsea_heatmap_%s.png)", mn, mn)); add("") }
  add("### 5.2 Tissue x L3"); add("")
  add("| Method | Pathways | Groups | Output |"); add("|---|---|---|---|")
  for (mn in get_ssgsea_method_names(ssgsea_results_l3_all)) {
    add(sprintf("| %s | %d | %d | `reports/ssgsea_l3_scores_%s.rds` |",
                mn, nrow(ssgsea_results_l3_all[[mn]]$scores), ncol(ssgsea_results_l3_all[[mn]]$scores), mn))
  }
  add("")
  for (mn in get_ssgsea_method_names(ssgsea_results_l3_all)) { add(sprintf("![L3 ssGSEA heatmap %s](figures/ssgsea_l3_heatmap_%s.png)", mn, mn)); add("") }
  if (exists("ssgsea_llm_df") && nrow(ssgsea_llm_df) > 0) {
    add(sprintf("### 5.3 ssGSEA LLM Summary (%d records)", nrow(ssgsea_llm_df))); add("")
    add("Structured outputs: `reports/ssgsea_llm_structured.tsv` | `LLM_SSGSEA_INTERPRETATION.md`")
    add("- LLM is run once per tissue x cell type group using combined ssGSEA evidence across multiple databases/methods, rather than once per database."); add("")
  }
}

add(sprintf("## 6. %s Clustering + Per-cluster ssGSEA + OFA", cluster_backend_label)); add("")
if (!RUN_CLUSTERING || length(choir_ofa_all) == 0) {
  add(sprintf("%s clustering not run or no clusters processed.", cluster_backend_label)); add("")
} else {
  cluster_col_report <- grep(sprintf("^%s_clusters", toupper(cluster_backend_label)), colnames(obj@meta.data), value = TRUE)
  cluster_col_report <- if (length(cluster_col_report) > 0) cluster_col_report[1] else NULL
  n_clusters_report <- if (!is.null(cluster_col_report) && cluster_col_report %in% colnames(obj@meta.data)) {
    length(unique(na.omit(as.character(obj@meta.data[[cluster_col_report]]))))
  } else {
    NA
  }
  cluster_param_report <- if (identical(cluster_backend, "LEIDEN")) {
    sprintf("resolution=%.3f", as.numeric(LEIDEN_RESOLUTION))
  } else {
    sprintf("alpha=%.3f", CHOIR_ALPHA)
  }
  add(sprintf("- **Backend:** %s | **Parameter:** %s | **Clusters:** %s | **Reduction:** `%s`",
              cluster_backend_label,
              cluster_param_report,
              ifelse(is.na(n_clusters_report), "unknown", n_clusters_report),
              if (exists("choir_reduction")) choir_reduction else "N/A")); add("")
  add(sprintf("### 6.1 %s UMAP", cluster_backend_label))
  add(sprintf("![%s UMAP](figures/%s.png)", cluster_backend_label, cluster_umap_stem)); add("")
  add(sprintf("![%s UMAP by Tissue](figures/%s.png)", cluster_backend_label, cluster_umap_split_stem)); add("")
  if (length(choir_results_all[["ssgsea"]]) > 0) {
    add(sprintf("### 6.2 ssGSEA per %s cluster", cluster_backend_label)); add("")
    for (mn in names(choir_results_all[["ssgsea"]])) {
      add(sprintf("![%s ssGSEA %s](figures/ssgsea_%s_heatmap_%s.png)", cluster_backend_label, mn, cluster_output_prefix, mn)); add("")
    }
  }
  add(sprintf("Cluster cell-count tables: `reports/%s/cluster_L2_counts.csv` | `reports/%s/cluster_L3_counts.csv` | `reports/%s/cluster_L2_L3_counts_long.csv`",
              cluster_output_prefix, cluster_output_prefix, cluster_output_prefix))
  add("")
  add("### 6.3 OFA One-vs-Rest Markers"); add("")
  add("| Cluster | n_focal | n_rest | Sig markers | Integrated LLM |"); add("|---|---|---|---|---|")
  for (cl_name in names(choir_ofa_all)) {
    r <- choir_ofa_all[[cl_name]]
    n_sig <- if (!is.null(r$de_table)) sum(r$de_table$padj < OFA_PADJ_THR, na.rm = TRUE) else 0
    has_integrated <- if (!is.null(r$llm$integrated)) r$llm$integrated$status else "skipped"
    add(sprintf("| %s | %d | %d | %d | %s |", cl_name, r$n_focal, r$n_rest, n_sig, has_integrated))
  }
  add("")
  if (exists("choir_llm_df") && nrow(choir_llm_df) > 0) {
    add(sprintf("### 6.4 %s Cluster LLM Summary (%d records)", cluster_backend_label, nrow(choir_llm_df))); add("")
    add(sprintf("Structured outputs: `reports/%s` | `%s`", cluster_llm_tsv_filename, cluster_llm_markdown_filename))
    add("- Language: Chinese (Simplified). key_drivers: English gene symbols.")
    add(sprintf("- Each %s cluster record includes annotation composition, annotation match degree, outlier assessment, discovery assessment, and integrated DEG + ssGSEA + OFA evidence.", cluster_backend_label)); add("")
  }
}

add("## 7. L3 OFA (Cell-level One-vs-Rest)"); add("")
if ((!is.data.frame(l3_ofa_vs_rest_summary_df) || nrow(l3_ofa_vs_rest_summary_df) == 0) &&
    (!is.data.frame(l3_ofa_same_l2_summary_df) || nrow(l3_ofa_same_l2_summary_df) == 0) &&
    (!is.data.frame(l3_ofa_inter_tissue_summary_df) || nrow(l3_ofa_inter_tissue_summary_df) == 0)) {
  add("L3 OFA not run or no valid comparisons passed minimum-cell thresholds."); add("")
} else {
  if (is.data.frame(l3_ofa_vs_rest_summary_df) && nrow(l3_ofa_vs_rest_summary_df) > 0) {
    add(sprintf("### 7.1 L3 vs Other Cell Types (%d comparisons)", nrow(l3_ofa_vs_rest_summary_df))); add("")
    add("Structured outputs: `reports/l3_ofa_vs_rest_summary.tsv` and per-L3 folders under `reports/l3_ofa/vs_other_cell_types/` (each folder includes `bubbleplot_overview.png/.pdf`)."); add("")
    add("| L3 | L2 | n_focal | n_rest | Sig markers |"); add("|---|---|---|---|---|")
    for (i in seq_len(nrow(l3_ofa_vs_rest_summary_df))) {
      rec <- l3_ofa_vs_rest_summary_df[i, , drop = FALSE]
      add(sprintf("| %s | %s | %d | %d | %d |", rec$celltype_label, rec$celltype_l2, rec$n_focal, rec$n_rest, rec$n_sig))
    }
    add("")
  }
  if (is.data.frame(l3_ofa_same_l2_summary_df) && nrow(l3_ofa_same_l2_summary_df) > 0) {
    add(sprintf("### 7.2 L3 vs Same-L2 Other L3 (%d comparisons)", nrow(l3_ofa_same_l2_summary_df))); add("")
    add("Structured outputs: `reports/l3_ofa_same_l2_summary.tsv` and per-L3 folders under `reports/l3_ofa/vs_same_l2_other_l3/` (each folder includes `bubbleplot_overview.png/.pdf`)."); add("")
    add("| L3 | L2 | n_focal | n_rest | Sig markers |"); add("|---|---|---|---|---|")
    for (i in seq_len(nrow(l3_ofa_same_l2_summary_df))) {
      rec <- l3_ofa_same_l2_summary_df[i, , drop = FALSE]
      add(sprintf("| %s | %s | %d | %d | %d |", rec$celltype_label, rec$celltype_l2, rec$n_focal, rec$n_rest, rec$n_sig))
    }
    add("")
  }
  if (is.data.frame(l3_ofa_inter_tissue_summary_df) && nrow(l3_ofa_inter_tissue_summary_df) > 0) {
    add(sprintf("### 7.3 Same L3 Across Tissues (%d comparisons)", nrow(l3_ofa_inter_tissue_summary_df))); add("")
    add("Structured outputs: `reports/l3_ofa_inter_tissue_summary.tsv` and per-comparison folders under `reports/l3_ofa/same_l3_other_tissues/` (each folder includes `bubbleplot_overview.png/.pdf`)."); add("")
    add("| L3 | L2 | Comparison | n_focal | n_rest | Sig markers |"); add("|---|---|---|---|---|---|")
    for (i in seq_len(nrow(l3_ofa_inter_tissue_summary_df))) {
      rec <- l3_ofa_inter_tissue_summary_df[i, , drop = FALSE]
      add(sprintf("| %s | %s | %s | %d | %d | %d |", rec$celltype_label, rec$celltype_l2, rec$comparison, rec$n_focal, rec$n_rest, rec$n_sig))
    }
    add("")
  }
}

add("## 8. LLM Interpretation (interpret_agent)"); add("")
if (!ENABLE_LLM) {
  add("**SKIPPED:** DEEPSEEK_API_KEY not set. Re-run with API key to enable."); add("")
} else {
  add("### 8.1 Tissue-Pair DESeq2 LLM"); add("")
  if (nrow(agent_structured_df) == 0) {
    add("No interpret_agent results generated."); add("")
  } else {
    add("Structured outputs: `reports/interpret_agent_structured.tsv` | `LLM_INTERPRETATION.md`")
    add("- Language: Chinese (Simplified). Validated gene->pathway map used.")
    add("- Pairwise comparisons are interpreted once per comparison using integrated up/down evidence across multiple databases."); add("")
    for (i in seq_len(nrow(agent_structured_df))) {
      rec <- agent_structured_df[i, , drop = FALSE]
      add(sprintf("#### %s | %s | %s | %s", rec$celltype_level, rec$celltype_label, rec$comparison, rec$direction)); add("")
      add(sprintf("- **Status:** %s | **Source DB:** %s", rec$status, ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
      if (nzchar(rec$overview))   add(sprintf("- **Overview:** %s", rec$overview))
      if (nzchar(rec$hypothesis)) add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
      add("")
    }
  }
  add("### 8.2 ssGSEA Group LLM (Tissue x Cell Type)"); add("")
  if (!exists("ssgsea_llm_df") || nrow(ssgsea_llm_df) == 0) {
    add("No ssGSEA LLM results generated."); add("")
  } else {
    add(sprintf("**%d records** across %d unique tissue x cell type groups using combined ssGSEA database evidence.",
                nrow(ssgsea_llm_df), length(unique(ssgsea_llm_df$comparison))))
    add("Full details: `LLM_SSGSEA_INTERPRETATION.md`"); add("")
    for (i in seq_len(nrow(ssgsea_llm_df))) {
      rec <- ssgsea_llm_df[i, , drop = FALSE]
      add(sprintf("#### %s | %s | %s", rec$celltype_level, rec$comparison, rec$direction)); add("")
      add(sprintf("- **Status:** %s | **DB:** %s", rec$status, ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
      if (nzchar(rec$overview))   add(sprintf("- **Overview:** %s", rec$overview))
      if (nzchar(rec$hypothesis)) add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
      add("")
    }
  }
  add(sprintf("### 8.3 %s Cluster LLM (One-vs-Rest)", cluster_backend_label)); add("")
  if (!exists("choir_llm_df") || nrow(choir_llm_df) == 0) {
    add(if (!RUN_CLUSTERING) sprintf("%s clustering was not run; cluster LLM skipped.", cluster_backend_label) else sprintf("No %s cluster LLM results generated.", cluster_backend_label)); add("")
  } else {
    add(sprintf("**%d records** (%d %s clusters; integrated up/down evidence per cluster).", nrow(choir_llm_df), length(unique(choir_llm_df$comparison)), cluster_backend_label))
    add(sprintf("Full details: `%s`", cluster_llm_markdown_filename)); add("")
    for (i in seq_len(nrow(choir_llm_df))) {
      rec <- choir_llm_df[i, , drop = FALSE]
      add(sprintf("#### %s | %s | %s", rec$celltype_level, rec$comparison, rec$direction)); add("")
      add(sprintf("- **Status:** %s | **DB:** %s", rec$status, ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
      if (nzchar(rec$overview))       add(sprintf("- **Overview:** %s", rec$overview))
      if (nzchar(rec$key_mechanisms)) add(sprintf("- **Key Mechanisms:** %s", rec$key_mechanisms))
      add("")
    }
  }
}

add("## 9. Methods"); add("")
add("- **DE:** Pseudobulk DESeq2 (Squair et al. 2021 Nat Commun)")
add("- **Exploratory:** Cell-level Wilcoxon (marker discovery ONLY, NOT inference)")
add(sprintf("- **Enrichment:** clusterProfiler::enricher() + 8 databases (GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, %s)", CUSTOM_DB_LABEL))
add(sprintf("- **ssGSEA:** GSVA::ssgseaParam on tissue x L2/L3 grouped average expression (not sample-level pseudobulk); methods: %s", paste(SSGSEA_METHODS, collapse = ", ")))
add("- **ssGSEA LLM [LLM-4]:** standardize_result_with_llm once per tissue x cell-type group using combined top pathways across multiple ssGSEA databases")
add("- **Pairwise comparison LLM:** integrated up/down DEG plus combined multi-database enrichment are provided together in one record per comparison")
add(sprintf("- **%s cluster LLM [LLM-5]:** standardize_result_with_llm once per %s cluster using integrated one-vs-rest DEG, cluster-level ssGSEA, and multi-database enrichment", cluster_backend_label, cluster_backend_label))
add("- **L3 OFA:** cell-level one-vs-rest Wilcoxon marker analysis for (i) each L3 vs all other cell types, (ii) each L3 vs same-L2 sibling L3 states, and (iii) one tissue vs all remaining tissues within the same L3")
if (ENABLE_LLM) {
  add("- **LLM:** DeepSeek-based structured interpretation with integrated evidence prompts")
  add("- **LLM language:** Chinese (Simplified); gene symbols retained in English")
  add("- **Gene emphasis:** top-|log2FC| genes explicitly linked to enriched pathways")
} else {
  add("- **LLM:** SKIPPED (no API key)")
}
add("")
add("### Output Objects"); add("")
add(sprintf("- `%s.rds` -- Seurat object (`lineage_family` + `lineage_branch` + `cell_type_L2` + `cell_type_L3`)", FINAL_FILE_PREFIX))
add(sprintf("- `%s.h5ad` -- AnnData object (`lineage_family` + `lineage_branch` + `cell_type_L2` + `cell_type_L3`, `layers['counts']`, schema tags in `uns`)", FINAL_FILE_PREFIX))
add(sprintf("- `reports/%s` -- structured cluster-level LLM records", cluster_llm_tsv_filename))
add("- `reports/l3_ofa_vs_rest_summary.tsv` -- per-L3 vs all other cell types OFA summary")
add("- `reports/l3_ofa_same_l2_summary.tsv` -- per-L3 vs same-L2 sibling L3 OFA summary")
add("- `reports/l3_ofa_inter_tissue_summary.tsv` -- same-L3 across tissues OFA summary")
add(sprintf("- `%s/` -- pertpy Milo tissue DA outputs for L2/L3 subclusters", MILOPY_OUTPUT_DIRNAME))
add("")
add("---")
add(sprintf("*Generated by %s*", GENERATED_BY_LABEL))

writeLines(md, file.path(OUTPUT_DIR, "REPORT.md"))
cat("[OK] REPORT.md written\n")

# ==============================================================================
# 12. Save Final Object (RDS + h5ad)
# ==============================================================================

cat("\n=== Saving Final Object (RDS + h5ad) ===\n")
if (!"cell_type_L3" %in% colnames(obj@meta.data)) {
  obj@meta.data[["cell_type_L3"]] <- as.character(obj@meta.data[[L3_SOURCE_COL]])
  cat(sprintf("[OK] Created cell_type_L3 from '%s'\n", L3_SOURCE_COL))
}
stopifnot("cell_type_L2" %in% colnames(obj@meta.data))
stopifnot("cell_type_L3" %in% colnames(obj@meta.data))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L2"]])))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L3"]])))
cat(sprintf("[OK] cell_type_L2: %d unique\n", length(unique(obj@meta.data[["cell_type_L2"]]))))
cat(sprintf("[OK] cell_type_L3: %d unique\n", length(unique(obj@meta.data[["cell_type_L3"]]))))

final_schema <- NULL
if (exists("tc_apply_final_schema_to_seurat", mode = "function")) {
  obj <- tc_apply_final_schema_to_seurat(
    obj,
    lineage_tag = LINEAGE_TAG,
    generated_by = GENERATED_BY_LABEL,
    final_file_prefix = FINAL_FILE_PREFIX,
    default_latent_reduction_candidates = FINAL_SCHEMA_LATENT_REDUCTION_CANDIDATES,
    default_umap_reduction_candidates = FINAL_SCHEMA_UMAP_REDUCTION_CANDIDATES,
    schema_version = FINAL_SCHEMA_VERSION,
    object_role = FINAL_SCHEMA_OBJECT_ROLE
  )
  final_schema <- obj@misc$tc_schema
  cat(sprintf(
    "[OK] final schema prepared: family=%s branch=%s latent=%s umap=%s\n",
    final_schema$lineage_family,
    final_schema$lineage_branch,
    ifelse(nzchar(final_schema$default_latent_obsm), final_schema$default_latent_obsm, "NA"),
    ifelse(nzchar(final_schema$default_umap_obsm), final_schema$default_umap_obsm, "NA")
  ))
}

rds_path <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".rds"))
saveRDS(obj, rds_path)
cat(sprintf("[OK] RDS saved: %s (%.1f MB)\n", basename(rds_path), file.size(rds_path) / 1e6))

h5ad_out_path <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".h5ad"))
tryCatch({
  anndata      <- reticulate::import("anndata",       convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse",  convert = FALSE)
  np           <- reticulate::import("numpy",          convert = FALSE)
  builtins     <- reticulate::import("builtins",       convert = FALSE)
  counts_mat      <- GetAssayData(obj, layer = "counts")
  norm_layer_name <- if ("data" %in% Layers(obj[["RNA"]])) "data" else "counts"
  expr_mat        <- GetAssayData(obj, layer = norm_layer_name)
  counts_scipy <- matrix_to_scipy_csr(counts_mat, scipy_sparse, np)
  expr_scipy   <- if (identical(norm_layer_name, "counts")) counts_scipy else matrix_to_scipy_csr(expr_mat, scipy_sparse, np)
  meta_export  <- sanitize_obs_for_h5ad(obj@meta.data)
  obs_df       <- reticulate::r_to_py(meta_export)
  var_df       <- data.frame(gene_symbol = rownames(obj), row.names = rownames(obj), stringsAsFactors = FALSE)
  var_py       <- reticulate::r_to_py(var_df)
  adata <- anndata$AnnData(X = expr_scipy, obs = obs_df, var = var_py)
  adata$layers$`__setitem__`("counts", counts_scipy)
  adata$uns$`__setitem__`("X_layer", norm_layer_name)
  for (red_name in Reductions(obj)) {
    emb <- Embeddings(obj, reduction = red_name)
    adata$obsm$`__setitem__`(paste0("X_", red_name), np$array(emb, dtype = np$float32))
  }
  if (!is.null(final_schema) && exists("tc_write_final_schema_to_anndata", mode = "function")) {
    tc_write_final_schema_to_anndata(adata, final_schema)
  }
  adata$write_h5ad(h5ad_out_path, compression = "gzip")
  cat(sprintf("[OK] h5ad saved: %s (%.1f MB)\n", basename(h5ad_out_path), file.size(h5ad_out_path) / 1e6))
  adata_check <- anndata$read_h5ad(h5ad_out_path)
  obs_cols    <- reticulate::py_to_r(builtins$list(adata_check$obs$columns))
  layer_keys  <- reticulate::py_to_r(builtins$list(adata_check$layers$keys()))
  stopifnot("cell_type_L2" %in% obs_cols, "cell_type_L3" %in% obs_cols, "counts" %in% layer_keys)
  if (!is.null(final_schema) && exists("tc_validate_final_h5ad_contract", mode = "function")) {
    tc_validate_final_h5ad_contract(adata_check, builtins, schema = final_schema)
  }
  cat(sprintf("[OK] h5ad verified: %d cells x %d genes, final schema + cell_type_L2/L3 + counts layer present\n",
              reticulate::py_to_r(adata_check$n_obs), reticulate::py_to_r(adata_check$n_vars)))
  rm(adata, adata_check, builtins, counts_mat, expr_mat, counts_scipy, expr_scipy, meta_export, var_df); gc()
}, error = function(e) {
  cat(sprintf("[ERROR] h5ad export failed: %s\n", e$message))
  cat("[INFO] RDS saved successfully.\n")
})

milopy_run <- list(status = "not_run", reason = "RUN_MILOPY=FALSE")
if (isTRUE(RUN_MILOPY)) {
  cat("\n=== MiloPy tissue differential abundance (pertpy) ===\n")
  if (file.exists(h5ad_out_path)) {
    milopy_dir <- file.path(OUTPUT_DIR, MILOPY_OUTPUT_DIRNAME)
    milopy_run <- tc_run_tissue_milopy(
      input_h5ad = h5ad_out_path,
      output_dir = milopy_dir,
      celltype_cols = MILOPY_LEVEL_COLS,
      sample_col = SAMPLE_COL,
      tissue_col = TISSUE_COL,
      conda_exe = MILOPY_CONDA_EXE,
      env_prefix = MILOPY_ENV_PREFIX,
      script_path = MILOPY_SCRIPT_PATH,
      latent_key_candidates = MILOPY_LATENT_KEY_CANDIDATES,
      umap_key_candidates = MILOPY_UMAP_KEY_CANDIDATES,
      min_cells_per_celltype = MILOPY_MIN_CELLS_PER_CELLTYPE,
      min_cells_per_sample = MILOPY_MIN_CELLS_PER_SAMPLE,
      min_samples_per_tissue = MILOPY_MIN_SAMPLES_PER_TISSUE,
      n_neighbors = MILOPY_N_NEIGHBORS,
      nhood_prop = MILOPY_NHOOD_PROP,
      alpha = MILOPY_ALPHA,
      random_seed = MILOPY_RANDOM_SEED,
      max_celltypes_per_level = MILOPY_MAX_CELLTYPES_PER_LEVEL,
      make_plots = MILOPY_MAKE_PLOTS,
      write_milo_h5ad = MILOPY_WRITE_MILO_H5AD,
      log_path = file.path(milopy_dir, MILOPY_LOG_FILENAME),
      stop_on_error = MILOPY_STOP_ON_ERROR
    )
    saveRDS(milopy_run, file.path(rpt_dir, "milopy_run.rds"))
    if (is.data.frame(milopy_run$level_summary) && nrow(milopy_run$level_summary) > 0) {
      fwrite(milopy_run$level_summary, file.path(rpt_dir, "milopy_level_summary.tsv"), sep = "\t")
    }
    if (is.data.frame(milopy_run$pairwise_summary) && nrow(milopy_run$pairwise_summary) > 0) {
      fwrite(milopy_run$pairwise_summary, file.path(rpt_dir, "milopy_pairwise_summary.tsv"), sep = "\t")
    }
    if (is.data.frame(milopy_run$contrast_summary) && nrow(milopy_run$contrast_summary) > 0) {
      fwrite(milopy_run$contrast_summary, file.path(rpt_dir, "milopy_contrast_summary.tsv"), sep = "\t")
    }
    if (identical(milopy_run$status, "completed")) {
      n_levels_done <- if (is.data.frame(milopy_run$level_summary)) nrow(milopy_run$level_summary) else 0L
      cat(sprintf("[OK] MiloPy complete: %d level summaries | outputs: %s\n", n_levels_done, milopy_dir))
    } else {
      cat(sprintf("[WARN] MiloPy failed (exit=%s). See log: %s\n", milopy_run$exit_status, milopy_run$log_path))
    }
  } else {
    milopy_run <- list(status = "skipped", reason = "final h5ad unavailable")
    cat("[WARN] Skipping MiloPy because final h5ad output is unavailable.\n")
  }
}

# ==============================================================================
# 13. Final Summary
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("%s\n", LINEAGE_COMPLETION_BANNER))
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Key outputs:\n")
cat("  REPORT.md\n")
cat("  LLM_INTERPRETATION.md          <- tissue-pair DESeq2 LLM (Chinese)\n")
cat("  LLM_SSGSEA_INTERPRETATION.md   <- ssGSEA group LLM (Chinese) [LLM-4]\n")
cat(sprintf("  %s    <- %s cluster LLM (Chinese) [LLM-5]\n", cluster_llm_markdown_filename, cluster_backend_label))
cat(sprintf("  %s.rds / .h5ad\n", FINAL_FILE_PREFIX))
cat("  figures/                        <- UMAP, dotplot, heatmap, composition\n")
cat("  pseudobulk_de/                  <- DESeq2 L2 + enrichment + volcano\n")
cat("  pseudobulk_de_L3/               <- DESeq2 L3 + enrichment + volcano\n")
cat("  wilcox_exploratory/             <- L2 cell-level wilcox (marker discovery)\n")
cat("  wilcox_exploratory_L3/          <- L3 cell-level wilcox (marker discovery)\n")
cat("  reports/ssgsea_scores_{method}.rds       <- L2 ssGSEA scores\n")
cat("  reports/ssgsea_z_{method}.rds            <- L2 ssGSEA Z-scores\n")
cat("  reports/ssgsea_l3_scores_{method}.rds    <- L3 ssGSEA scores\n")
cat("  reports/ssgsea_results_all.rds           <- complete L2 ssGSEA result structure\n")
cat("  reports/ssgsea_results_L3_all.rds        <- complete L3 ssGSEA result structure\n")
cat("  reports/ssgsea_llm_{level}/              <- per-group LLM interpretations [LLM-4]\n")
cat("  reports/ssgsea_llm_structured.tsv        <- all ssGSEA LLM records\n")
cat("  reports/ssgsea_llm_structured_*_all.rds  <- ssGSEA LLM structured R objects\n")
cat(sprintf("  reports/%s/%s_clusters.csv\n", cluster_output_prefix, cluster_output_prefix))
cat(sprintf("  reports/%s/cluster_L2_counts.csv         <- %s cluster x L2 cell counts\n", cluster_output_prefix, cluster_backend_label))
cat(sprintf("  reports/%s/cluster_L3_counts.csv         <- %s cluster x L3 cell counts\n", cluster_output_prefix, cluster_backend_label))
cat(sprintf("  reports/%s/cluster_L2_L3_counts_long.csv <- long-format %s cluster annotation counts\n", cluster_output_prefix, cluster_backend_label))
cat(sprintf("  reports/%s/ssgsea_{scores|z}_{method}.rds\n", cluster_output_prefix))
cat(sprintf("  reports/%s/ofa_{cluster}/markers.csv + enrichment + integrated LLM [LLM-5]\n", cluster_output_prefix))
cat(sprintf("  reports/%s         <- all %s cluster LLM records\n", cluster_llm_tsv_filename, cluster_backend_label))
cat(sprintf("  reports/%s         <- %s cluster LLM structured R object\n", cluster_llm_rds_filename, cluster_backend_label))
cat("  reports/l3_ofa_vs_rest_summary.tsv       <- L3 vs all other cell types OFA summary\n")
cat("  reports/l3_ofa_same_l2_summary.tsv       <- L3 vs same-L2 sibling L3 OFA summary\n")
cat("  reports/l3_ofa_inter_tissue_summary.tsv  <- same-L3 across tissues OFA summary\n")
cat("  reports/l3_ofa/                          <- detailed L3 OFA marker/enrichment outputs\n")
cat(sprintf("  %s/                          <- pertpy Milo tissue DA outputs for L2/L3\n", MILOPY_OUTPUT_DIRNAME))
cat("  reports/milopy_level_summary.tsv         <- Milo per-level execution summary\n")
cat("  reports/milopy_pairwise_summary.tsv      <- Milo per-comparison execution summary\n")
cat("  reports/milopy_contrast_summary.tsv      <- Milo significant-neighborhood summary\n")
cat("  reports/interpret_agent_structured.tsv   <- tissue-pair DESeq2 LLM records\n")
cat("\n")
cat(sprintf("%s\n", PIPELINE_CHANGELOG_TITLE))
for (line in PIPELINE_CHANGELOG_LINES) cat(sprintf("%s\n", line))
cat("\n")
cat(sprintf("%s\n", PIPELINE_INHERITED_FIXES_TITLE))
for (line in PIPELINE_INHERITED_FIX_LINES) cat(sprintf("%s\n", line))
cat("\n")
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
cat("[OK] Done\n")