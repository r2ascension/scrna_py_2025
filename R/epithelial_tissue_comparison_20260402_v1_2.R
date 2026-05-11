#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison Pipeline v1.2.0 - PRODUCTION
# ==============================================================================
#
# Purpose:
#   Respiratory epithelial tissue comparison across anatomical sites:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. L3 -> L2 remapping (fine epithelial states -> major lineages)
#     3. Visualization (UMAP, marker dotplot/heatmap, composition)
#     4. Pseudobulk DESeq2 per L2/L3 lineage x tissue pair
#     5. Exploratory cell-level wilcox per L2/L3 lineage (marker discovery only)
#     6. Multi-database enrichment (GO BP/MF/CC, KEGG, Hallmark,
#        CellMarker, PanglaoDB, custom epithelial markers)
#     7. interpret_agent (DeepSeek) with epithelial tissue-pair-specific context
#     8. Structured REPORT.md (png embeds)
#     [NEW] 8.5. Pseudobulk ssGSEA by tissue x L2/L3
#     [NEW] 8.6. CHOIR clustering (scANVI latent) + per-cluster ssGSEA + OFA
#
# v1.2.0 Changes (vs v1.1.1):
#   [NEW-12] Create a standalone v1.2 release script with independent output root
#   [FIX-1] Unify version / script-name / output-path strings via constants
#   [FIX-2] Correct REPORT.md footer and summary labels to avoid stale version text
#
# v1.1.1 Changes (vs v1.1.0):
#   [NEW-10] Add full L3 tissue-comparison workflow alongside L2
#            (DE/enrichment/interpret/wilcox/ssGSEA/reporting)
#   [NEW-11] Remove IG-related genes from DE / Wilcoxon / enrichment / ssGSEA / OFA
#
# v1.1.0 Changes (vs v1.0.0):
#   [NEW-1] ENRICHMENT_GS_SIZE_RULES: adaptive min/maxGSSize per database
#   [NEW-2] get_enrichment_size_rule(): helper for per-DB size lookup
#   [NEW-3] run_gmt_enrichment(): upgraded to use adaptive size rules (matches B cell v2.2)
#   [NEW-4] ssGSEA helpers: term2gene_to_list(), map_gene_sets_to_features(), filter_gs_size()
#   [NEW-5] Section 8.5: pseudobulk ssGSEA by tissue x L2 (hallmark + go_bp default)
#   [NEW-6] Section 8.6: CHOIR clustering on scANVI latent + per-cluster ssGSEA (8.6a) + OFA 8-DB (8.6b)
#   [NEW-7] Library loading: GSVA, BiocParallel, CHOIR
#   [NEW-8] Config params: RUN_SSGSEA, SSGSEA_*, RUN_CHOIR, CHOIR_*, RUN_OFA, OFA_*
#   [NEW-9] REPORT.md sections 4.5 (ssGSEA) + 4.6 (CHOIR)
#
# Notes:
#   - This is intended for cross-site anatomical comparison of respiratory
#     epithelial cells. Unless disease covariates are explicitly modeled,
#     do NOT over-interpret site-associated differences as disease effects.
#   - Uses the finalized epithelial scANVI SELF reference with curated L3 labels.
#
# Author: r2end + GitHub Copilot
# Date:   2026-04-02
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

# ----- Input / Output -----
PIPELINE_VERSION <- "v1.2.0"
SCRIPT_BASENAME  <- "epithelial_tissue_comparison_20260402_v1_2.R"
OUTPUT_TAG       <- "epithelial_tissue_comparison_v1_2"

H5AD_PATH  <- "/home/h2048/data/core20260322/epithelial_scanvi_v2_7_HOTFIX_SELF_final.h5ad"
OUTPUT_DIR <- file.path("/home/h2048/data/R/0402", OUTPUT_TAG)

# ----- L3 -> L2 Remapping -----
L3_SOURCE_COL <- "cell_type_L3"

L3_TO_L2_REMAP <- c(
  "AT1_Canonical"               = "Alveolar",
  "AT1_MatrixRemodeling"        = "Alveolar",
  "AT2"                         = "Alveolar",
  "AT2_Cycling"                 = "Alveolar",
  "Basal_Progenitor"            = "Basal_Lineage",
  "Basal_Cycling"               = "Basal_Lineage",
  "Basal_Inflammatory"          = "Basal_Lineage",
  "Basal_EMT_ECM"               = "Basal_Lineage",
  "Suprabasal_Progenitor"       = "Basal_Lineage",
  "Suprabasal_Cycling"          = "Basal_Lineage",
  "Ciliated_Mature"             = "Ciliated_Lineage",
  "Ciliogenesis_Deuterosomal"   = "Ciliated_Lineage",
  "Ciliated_Cycling_Immature"   = "Ciliated_Lineage",
  "Goblet"                      = "Secretory_Lineage",
  "Club"                        = "Secretory_Lineage",
  "SMG_Mucous"                  = "Secretory_Lineage",
  "Goblet_Defense_DUOX2"        = "Secretory_Lineage",
  "SMG_Serous"                  = "SMG",
  "SMG_Duct_Secretory_Defense"  = "SMG",
  "Squamous_Metaplasia"         = "Rare_Specialized",
  "Ionocyte_Brush"              = "Rare_Specialized"
)

# ----- Reference Database Paths -----
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL      <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH  <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"

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

ensure_env_placeholder <- function(path, key) {
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    if (any(grepl(sprintf("^%s\\s*=", key), trimws(lines)))) return(invisible(FALSE))
    write(sprintf("%s=YOUR_%s_HERE", key, key), file = path, append = TRUE)
    return(invisible(TRUE))
  }
  writeLines(sprintf("%s=YOUR_%s_HERE", key, key), path)
  invisible(TRUE)
}

env_files_loaded <- ENV_FILE_CANDIDATES[file.exists(ENV_FILE_CANDIDATES)]
invisible(lapply(env_files_loaded, load_env_file))

# ----- DeepSeek API (OPTIONAL) -----
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
ENABLE_LLM <- nchar(DEEPSEEK_API_KEY) >= 10
if (!ENABLE_LLM) {
  ensure_env_placeholder("/home/h2048/.env", "DEEPSEEK_API_KEY")
  cat("[WARN] DEEPSEEK_API_KEY not set. interpret_agent will be SKIPPED.\n")
} else if (length(env_files_loaded) > 0) {
  cat(sprintf("[OK] Loaded environment file(s): %s\n", paste(env_files_loaded, collapse = ", ")))
}
INTERPRET_AGENT_MODEL                           <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS                      <- 50
INTERPRET_AGENT_ADD_PPI                         <- TRUE
INTERPRET_AGENT_MAX_RETRIES                     <- 3L
INTERPRET_AGENT_RETRY_SLEEP_SEC                 <- 2
INTERPRET_AGENT_MIN_TERMS                       <- 3L
INTERPRET_AGENT_DISABLE_PPI_ON_STRINGDB_TIMEOUT <- TRUE
INTERPRET_AGENT_DB_PRIORITY <- c(
  "GO_BP", "Hallmark", "KEGG", "GO_MF", "GO_CC",
  "CellMarker", "PanglaoDB", "Epithelial_custom"
)
STANDARDIZE_LLM_OUTPUT          <- TRUE
STANDARDIZE_LLM_MODEL           <- "deepseek-chat"
STANDARDIZE_LLM_MAX_RETRIES     <- 3L
STANDARDIZE_LLM_RETRY_SLEEP_SEC <- 2L
STANDARDIZE_LLM_MAX_INPUT_CHARS <- 12000L

# ----- Python / reticulate -----
PYTHON_CONDA_ENV <- "bbknn_env"

# ----- Metadata Column Names -----
TISSUE_COL      <- "tissue"
SAMPLE_COL      <- "sample"
CELLTYPE_L2_COL <- "cell_type_L2"
CELLTYPE_L3_COL <- "cell_type_L3"
LABEL_COL       <- "cell_type_L3"

UMAP_REDUCTION_PREFERRED <- c(
  "umap_fine", "umap_major",
  "umap_scanvi_fine", "umap_scanvi_major",
  "umap_scanvi", "umap_harmony", "umap"
)
UMAP_PT_SIZE <- 0.35
UMAP_TISSUE_COLORS <- c(
  "lung parenchyma"    = "#D55E00",
  "nose"               = "#009E73",
  "respiratory airway" = "#0072B2",
  "sinus"              = "#CC79A7"
)

# ----- Pseudobulk DE Parameters -----
MIN_CELLS_PER_PSEUDOBULK              <- 15
MIN_SAMPLES_PER_TISSUE                <- 3
PADJ_THR                              <- 0.05
LFC_THR                               <- 1.0
DESIGN_COVARIATE                      <- NULL
PSEUDOBULK_DISAMBIGUATION_CANDIDATES  <- c("dataset", "source", "batch", "orig.ident")

# ----- Reproducibility -----
set.seed(42)

current_future_maxsize <- getOption("future.globals.maxSize")
if (is.null(current_future_maxsize) || !is.numeric(current_future_maxsize) || !is.finite(current_future_maxsize)) {
  current_future_maxsize <- 0
}

# ----- Exploratory Wilcox Parameters -----
WILCOX_MIN_CELLS            <- 100
WILCOX_LFC_THR              <- 0.25
WILCOX_TOP_N_PER_DIRECTION  <- 100
WILCOX_MAX_CELLS_PER_IDENT  <- 10000

# ----- Analysis gene filtering -----
FILTER_IG_GENES             <- TRUE
IG_GENE_REGEX               <- "^(IGH|IGK|IGL)"

# ----- Enrichment -----
TOP_N_DEG_ENRICHMENT <- 200

# [NEW-1] Adaptive gene-set size rules per database
ENRICHMENT_GS_SIZE_RULES <- list(
  default            = c(min = 10L, max = 500L),
  CellMarker         = c(min = 3L,  max = 200L),
  PanglaoDB          = c(min = 3L,  max = 200L),
  Epithelial_custom  = c(min = 2L,  max = 100L)
)

# ----- Visualization -----
HEATMAP_CELLS_PER_TYPE <- 120

# ----- [NEW] ssGSEA (pseudobulk, tissue x L2/L3) -----
RUN_SSGSEA         <- TRUE
SSGSEA_METHODS     <- c("hallmark", "go_bp")   # primary + mechanistic depth
SSGSEA_N_TOP       <- 20                        # top pathways per group per method (CSV)
SSGSEA_N_HEATMAP   <- 30                        # pathways in heatmap (by mean |Z|)
SSGSEA_CUSTOM_GMT  <- NULL                      # path to custom GMT; NULL = skip
N_CORES            <- 4                         # BiocParallel SnowParam workers

# ----- [NEW] CHOIR clustering (scANVI latent space) -----
RUN_CHOIR                  <- TRUE
CHOIR_REDUCTION_CANDIDATES <- c("scanvi_fine", "scanvi_major", "scanvi", "scvi", "harmony", "pca")
CHOIR_ALPHA                <- 0.05
CHOIR_FUTURE_GLOBALS_MAXSIZE <- 8 * 1024^3

options(future.globals.maxSize = max(current_future_maxsize, CHOIR_FUTURE_GLOBALS_MAXSIZE))

# ----- [NEW] OFA: one-vs-rest FindMarkers on CHOIR clusters -----
RUN_OFA                 <- TRUE
OFA_MIN_CELLS_FOCAL     <- 50
OFA_MIN_CELLS_REST      <- 100
OFA_PADJ_THR            <- 0.05
OFA_LFC_THR             <- 0.25
OFA_TOP_N               <- 200
OFA_MAX_CELLS_PER_IDENT <- 5000

# ----- [NEW] ssGSEA for CHOIR clusters -----
SSGSEA_CHOIR_METHODS   <- c("hallmark", "go_bp")
SSGSEA_CHOIR_N_TOP     <- 20
SSGSEA_CHOIR_N_HEATMAP <- 30

# ----- Epithelial Known Markers -----
KNOWN_MARKERS <- c(
  "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19",
  "AGER", "HOPX", "CAV1", "SFTPC", "SFTPA1", "SFTPA2", "SFTPB",
  "ABCA3", "NAPSA", "SLC34A2", "CHI3L1",
  "TP63", "KRT5", "KRT14", "KRT15", "ITGA6", "NGFR", "KRT4", "KRT13",
  "KRT17", "FN1", "COL17A1", "VIM",
  "FOXJ1", "TPPP3", "DNAH5", "DNAH9", "RSPH1", "DEUP1", "CCNO", "MCIDAS",
  "SCGB1A1", "SCGB3A1", "SCGB3A2", "SPDEF", "FOXA3", "MUC5AC", "MUC5B",
  "FCGBP", "TFF3", "DUOX2", "DUOXA2", "LCN2", "BPIFA2", "CEACAM5",
  "LTF", "LYZ", "SLPI", "DMBT1", "BPIFA1", "AZGP1", "WFDC2", "PIGR", "TCN1",
  "SPRR1A", "SPRR2A", "IVL", "KRT6A", "S100A7", "FOXI1", "ASCL3", "CFTR",
  "ATP6V0D2", "CLCNKA", "CLCNKB", "BSND",
  "MKI67", "TOP2A", "UBE2C", "BIRC5", "AURKB", "CENPA", "CCNB1"
)

# ----- Custom Epithelial Marker Database (for enrichment) -----
EPITHELIAL_MARKERS_DB <- data.frame(
  subtype = c(
    "AT1_Canonical", "AT1_MatrixRemodeling", "AT2", "AT2_Cycling",
    "Basal_Progenitor", "Basal_Cycling", "Basal_Inflammatory", "Basal_EMT_ECM",
    "Suprabasal_Progenitor", "Suprabasal_Cycling",
    "Ciliated_Mature", "Ciliogenesis_Deuterosomal", "Ciliated_Cycling_Immature",
    "Goblet", "Club", "SMG_Mucous", "Goblet_Defense_DUOX2",
    "SMG_Serous", "SMG_Duct_Secretory_Defense",
    "Squamous_Metaplasia", "Ionocyte_Brush"
  ),
  markers = c(
    "AGER,HOPX,CAV1,AQP4,RTKN2,CLDN18,EMP2",
    "AGER,CAV1,SPARC,COL4A1,COL4A2,SPOCK2",
    "SFTPC,SFTPA1,SFTPA2,SFTPB,ABCA3,NAPSA,SLC34A2,CHI3L1,CXCL8,SAA1",
    "MKI67,TOP2A,UBE2C,BIRC5,AURKB,CENPA,CCNB1",
    "KRT5,KRT14,TP63,KRT15,KRT19,ITGA6,NGFR",
    "KRT14,KRT5,TP63,MKI67,TOP2A,BIRC5",
    "KRT17,CXCL8,CXCL1,CXCL2,TNFAIP3,FOS,JUN",
    "KRT14,TP63,NGFR,FN1,COL17A1,MMP2,VIM",
    "KRT4,KRT13,KRT19,NOTCH1,NOTCH3,TP63,KRT5",
    "KRT4,KRT13,MKI67,TOP2A,BIRC5,TP63,KRT5",
    "FOXJ1,TPPP3,DNAH5,DNAH9,RSPH1,RFX2,RFX3",
    "DEUP1,CCNO,FOXN4,MCIDAS,CDC20B,E2F7,PLK4",
    "TPPP3,RSPH1,MKI67,TOP2A,FOXN4",
    "SCGB1A1,SCGB3A1,SCGB3A2,AGR2,AGR3,CYP2F1",
    "SCGB1A1,SCGB3A1,SFTPB,NAPSA,GPR116,CLDN18",
    "SPDEF,FOXA3,MUC5AC,MUC5B,FCGBP,TFF3,BPIFB2,AZGP1",
    "DUOX2,DUOXA2,LCN2,BPIFA2,CEACAM5",
    "LTF,LYZ,SLPI,DMBT1,BPIFA1,AZGP1,WFDC2",
    "PIGR,SCGB3A1,TCN1,WFDC2,DMBT1,SLPI",
    "SPRR1A,SPRR2A,SPRR2E,IVL,KRT6A,KLK7,S100A7",
    "FOXI1,ASCL3,CFTR,ATP6V0D2,CLCNKA,CLCNKB,BSND"
  ),
  stringsAsFactors = FALSE
)

# ----- Tissue-Pair Context for LLM -----
EPITHELIAL_BASE_CONTEXT <- paste(
  "Respiratory epithelial cells from human airway and alveolar tissues.",
  "This is primarily a cross-site anatomical comparison rather than a disease-vs-healthy design.",
  "Key epithelial programs include basal stem/progenitor maintenance, suprabasal transition,",
  "mucociliary differentiation, club/goblet/SMG secretory specialization, ion transport,",
  "alveolar surfactant biology, and barrier defense.",
  "Focus on regional variation in epithelial composition, differentiation state, barrier function,",
  "mucus production, glandular secretion, and mucociliary clearance."
)

TISSUE_CONTEXT <- list(
  "nose" = paste(
    "Nasal cavity: first-line barrier with intense environmental exposure,",
    "strong innate defense, mucus production, goblet/club programs, and surface epithelial remodeling."
  ),
  "sinus" = paste(
    "Paranasal sinus: semi-enclosed mucosa with dependence on drainage and mucus clearance,",
    "secretory defense, and potential gland-associated epithelial specialization."
  ),
  "respiratory airway" = paste(
    "Conducting airways: mucociliary clearance, club-cell secretory biology, ciliated differentiation,",
    "submucosal gland contribution, and epithelial barrier maintenance."
  ),
  "lung parenchyma" = paste(
    "Lung parenchyma: alveolar gas-exchange surface with AT1/AT2 specialization,",
    "surfactant production, alveolar repair programs, and relatively less mucus-dominant biology."
  )
)

build_tissue_pair_context <- function(t1, t2, celltype_label, dir_name,
                                      celltype_level = "L2") {
  ctx1 <- TISSUE_CONTEXT[[t1]]
  ctx2 <- TISSUE_CONTEXT[[t2]]
  if (is.null(ctx1)) ctx1 <- paste(t1, ": no specific epithelial context available.")
  if (is.null(ctx2)) ctx2 <- paste(t2, ": no specific epithelial context available.")

  sprintf(
    paste(
      "%s",
      "\nComparing %s vs %s for %s-level group '%s'.",
      "Direction: genes %s-regulated in %s relative to %s.",
      "\n--- %s ---\n%s",
      "\n--- %s ---\n%s",
      "\nBiological question: What regional epithelial programs explain the transcriptional divergence?",
      "Consider lineage balance (basal, ciliated, secretory, SMG, alveolar, rare cells),",
      "differentiation trajectories, barrier defense, mucus biology, ion transport,",
      "glandular secretion, and alveolar specialization.",
      "\n\nOutput requirements:",
      "1. Write ALL interpretive content in Chinese (Simplified Chinese).",
      "2. Explicitly identify the genes with the largest absolute log2FoldChange and explain WHICH enriched pathway(s)",
      "   they belong to and WHY those pathways matter for epithelial site biology.",
      "3. For each highlighted gene, state its log2FC and top pathway(s) using:",
      "   gene(log2FC=X.X) -> pathway name."
    ),
    EPITHELIAL_BASE_CONTEXT,
    t2, t1, celltype_level, celltype_label, dir_name, t2, t1,
    t2, ctx2, t1, ctx1
  )
}

# ==============================================================================
# 2. Load Libraries
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("Epithelial Tissue Comparison Pipeline %s (PRODUCTION)\n", PIPELINE_VERSION))
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
  library(jsonlite)
  library(fanyi)
  library(GSVA)           # [NEW]
  library(BiocParallel)   # [NEW]
  library(CHOIR)          # [NEW]
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

# ==============================================================================
# 3. Load GMT + Reference Databases
# ==============================================================================

cat("=== Loading Reference Databases ===\n")

gmt_all <- tryCatch({
  g <- read.gmt(MSIGDB_GMT_PATH); g$gene <- toupper(g$gene)
  cat(sprintf("[OK] MSigDB: %d gene sets\n", length(unique(g$term)))); g
}, error = function(e) { cat(sprintf("[WARN] MSigDB failed: %s\n", e$message)); NULL })

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
}, error = function(e) cat(sprintf("[WARN] GO GMT failed: %s\n", e$message)))

cellmarker_t2g <- NULL
tryCatch({
  cm_db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
  cm_db <- cm_db[grepl("Human", species, ignore.case = TRUE)]
  cat(sprintf("[OK] CellMarker raw: %d human entries\n", nrow(cm_db)))

  epi_pattern     <- paste(c("Epithelial", "Basal", "Goblet", "Club", "Secretory", "Ciliated",
                              "Ionocyte", "Brush", "Tuft", "Deuterosomal", "Serous", "Duct",
                              "Airway", "Respiratory", "Alveolar", "AT1", "AT2", "Suprabasal",
                              "Squamous"), collapse = "|")
  airway_pattern  <- "Lung|Airway|Respiratory|Nasal|Sinus|Bronch|Trachea|Alveol"
  cm_epi <- cm_db[grepl(epi_pattern, cell_name, ignore.case = TRUE) &
                    grepl(airway_pattern, tissue_type, ignore.case = TRUE)]
  cat(sprintf("[OK] CellMarker epithelial-related: %d entries\n", nrow(cm_epi)))

  if (nrow(cm_epi) < 10) {
    cat("[WARN] Too few epithelial entries in CellMarker; skipping\n")
  } else {
    t2g_list <- list()
    for (i in seq_len(nrow(cm_epi))) {
      markers_raw <- cm_epi$marker[i]; cell_type <- cm_epi$cell_name[i]
      if (is.na(markers_raw) || markers_raw == "") next
      markers <- trimws(toupper(gsub('["\r\n\\[\\]]', '',
                                      unlist(strsplit(markers_raw, "[,;\\s]+")))))
      markers <- unique(markers[markers != "" & !is.na(markers)])
      if (length(markers) > 0)
        t2g_list[[length(t2g_list) + 1]] <- data.frame(term = cell_type, gene = markers,
                                                        stringsAsFactors = FALSE)
    }
    cellmarker_t2g <- bind_rows(t2g_list) %>% distinct()
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

  pdb_epi <- pdb[grepl(paste(c("Epithelial", "Basal", "Goblet", "Club", "Secretory", "Ciliated",
                                "Ionocyte", "Brush", "Tuft", "Serous", "Duct", "Airway",
                                "Alveolar", "AT1", "AT2", "Suprabasal", "Squamous"), collapse = "|"),
                         cell_type, ignore.case = TRUE)]
  cat(sprintf("[OK] PanglaoDB epithelial: %d markers (%d types)\n",
              nrow(pdb_epi), length(unique(pdb_epi$cell_type))))

  if (nrow(pdb_epi) < 5) {
    cat("[WARN] Too few epithelial markers in PanglaoDB; skipping\n")
  } else {
    panglaodb_t2g <- pdb_epi %>%
      select(cell_type, gene_symbol) %>%
      mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
      filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
      distinct() %>%
      rename(term = cell_type, gene = gene_symbol)
    cat(sprintf("[OK] PanglaoDB TERM2GENE: %d pairs (%d terms)\n",
                nrow(panglaodb_t2g), length(unique(panglaodb_t2g$term))))
  }
}, error = function(e) cat(sprintf("[WARN] PanglaoDB failed: %s\n", e$message)))

epithelial_custom_t2g <- EPITHELIAL_MARKERS_DB %>%
  separate_rows(markers, sep = ",") %>%
  mutate(markers = trimws(toupper(markers))) %>%
  filter(markers != "") %>%
  select(term = subtype, gene = markers) %>%
  distinct()
cat(sprintf("[OK] Custom epithelial markers: %d subtypes, %d pairs\n",
            length(unique(epithelial_custom_t2g$term)), nrow(epithelial_custom_t2g)))
cat("\n")

# ==============================================================================
# 4. Helper Functions
# ==============================================================================

safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

level_file_stem <- function(base_stem, level_name, keep_legacy_l2 = TRUE) {
  if (keep_legacy_l2 && identical(level_name, "L2")) return(base_stem)
  sprintf("%s_%s", base_stem, tolower(level_name))
}

normalize_meta_values <- function(x) {
  x <- trimws(as.character(x)); x[is.na(x) | x == ""] <- "<NA>"; x
}

is_ig_gene <- function(gene_vec) {
  gene_vec <- toupper(trimws(as.character(gene_vec)))
  !is.na(gene_vec) & gene_vec != "" & grepl(IG_GENE_REGEX, gene_vec, perl = TRUE)
}

filter_analysis_gene_vector <- function(gene_vec, context = NULL, verbose = FALSE) {
  gene_vec <- as.character(gene_vec)
  if (!FILTER_IG_GENES || length(gene_vec) == 0) return(gene_vec)
  keep <- !is_ig_gene(gene_vec)
  removed_n <- sum(!keep, na.rm = TRUE)
  if (isTRUE(verbose) && removed_n > 0) {
    cat(sprintf("[INFO] %s removed %d IG-related genes\n",
                ifelse(is.null(context), "Analysis gene filter:", paste0(context, ":")),
                removed_n))
  }
  unique(gene_vec[keep])
}

filter_marker_table_by_gene <- function(df, gene_col = "gene", context = NULL, verbose = FALSE) {
  if (is.null(df) || nrow(df) == 0 || !FILTER_IG_GENES || !gene_col %in% colnames(df)) return(df)
  keep <- !is_ig_gene(df[[gene_col]])
  removed_n <- sum(!keep, na.rm = TRUE)
  if (isTRUE(verbose) && removed_n > 0) {
    cat(sprintf("[INFO] %s removed %d IG-related rows\n",
                ifelse(is.null(context), "Marker table filter:", paste0(context, ":")),
                removed_n))
  }
  df[keep, , drop = FALSE]
}

filter_term2gene_for_analysis <- function(t2g, context = NULL) {
  if (is.null(t2g) || nrow(t2g) == 0 || !FILTER_IG_GENES) return(t2g)
  t2g %>% filter(!is_ig_gene(gene))
}

filter_gene_sets_for_analysis <- function(gene_sets) {
  if (!FILTER_IG_GENES || length(gene_sets) == 0) return(gene_sets)
  gene_sets <- lapply(gene_sets, function(gs) unique(gs[!is_ig_gene(gs)]))
  gene_sets[lengths(gene_sets) > 0]
}

pick_reduction <- function(obj, preferred) {
  red <- Reductions(obj); hit <- preferred[preferred %in% red]
  if (length(hit) > 0) hit[[1]] else NULL
}

build_umap_plot <- function(obj, reduction_name, group_col, title,
                            split_col = NULL, label = FALSE,
                            width = 10, height = 8, cols = NULL) {
  p <- DimPlot(obj, reduction = reduction_name, group.by = group_col,
               split.by = split_col, pt.size = UMAP_PT_SIZE, shuffle = TRUE,
               label = label, repel = label, cols = cols) +
    ggtitle(title) + coord_equal() + theme_classic(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
  list(plot = p, width = width, height = height)
}

# ----- [NEW-2] Adaptive gene-set size rules -----
get_enrichment_size_rule <- function(db_name) {
  rule <- ENRICHMENT_GS_SIZE_RULES[[db_name]]
  if (is.null(rule)) rule <- ENRICHMENT_GS_SIZE_RULES[["default"]]
  list(min = as.integer(rule[["min"]]), max = as.integer(rule[["max"]]))
}

# ----- [NEW-3] GMT enrichment with adaptive size rules -----
run_gmt_enrichment <- function(gene_list, t2g, db_name, tested_genes = NULL) {
  if (is.null(t2g) || nrow(t2g) == 0 || length(gene_list) < 5) return(NULL)

  gene_list <- filter_analysis_gene_vector(gene_list)
  tested_genes <- filter_analysis_gene_vector(tested_genes)
  t2g <- filter_term2gene_for_analysis(t2g, context = db_name)
  if (length(gene_list) < 5 || is.null(t2g) || nrow(t2g) == 0) return(NULL)

  t2g_use <- t2g %>%
    mutate(term = as.character(term), gene = toupper(as.character(gene))) %>%
    filter(!is.na(term), !is.na(gene), term != "", gene != "") %>%
    distinct(term, gene)
  if (nrow(t2g_use) == 0) return(NULL)

  universe_use <- if (!is.null(tested_genes)) {
    intersect(toupper(tested_genes), unique(t2g_use$gene))
  } else {
    unique(t2g_use$gene)
  }
  universe_use <- unique(universe_use)
  if (length(universe_use) < 5) {
    cat(sprintf("    [INFO] %s: too few genes in tested universe overlap (%d)\n", db_name, length(universe_use)))
    return(NULL)
  }

  t2g_use <- t2g_use %>% filter(gene %in% universe_use)
  if (nrow(t2g_use) == 0) { cat(sprintf("    [INFO] %s: no overlap with tested gene universe\n", db_name)); return(NULL) }

  size_rule  <- get_enrichment_size_rule(db_name)
  gs_sizes   <- t2g_use %>% count(term, name = "gs_size")
  max_gs_obs <- suppressWarnings(max(gs_sizes$gs_size, na.rm = TRUE))
  if (!is.finite(max_gs_obs) || max_gs_obs < 2) {
    cat(sprintf("    [INFO] %s: no gene sets with >=2 overlapping genes\n", db_name)); return(NULL)
  }

  min_gs  <- max(2L, min(size_rule$min, as.integer(max_gs_obs)))
  max_gs  <- max(min_gs, size_rule$max)
  valid_t <- gs_sizes %>% filter(gs_size >= min_gs, gs_size <= max_gs) %>% pull(term)
  if (length(valid_t) == 0) {
    cat(sprintf("    [INFO] %s: no gene sets remain after adaptive size filtering (max_obs=%d)\n",
                db_name, as.integer(max_gs_obs)))
    return(NULL)
  }

  t2g_use    <- t2g_use %>% filter(term %in% valid_t)
  genes_use  <- intersect(toupper(unique(gene_list)), unique(t2g_use$gene))
  if (length(genes_use) < 5) {
    cat(sprintf("    [INFO] %s: too few input genes after overlap filtering (%d)\n", db_name, length(genes_use)))
    return(NULL)
  }

  tryCatch(
    suppressMessages(
      enricher(gene = genes_use, TERM2GENE = t2g_use,
               universe = unique(t2g_use$gene),
               pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               pAdjustMethod = "BH", minGSSize = min_gs, maxGSSize = max_gs)
    ),
    error = function(e) { cat(sprintf("    [WARN] %s: %s\n", db_name, e$message)); NULL }
  )
}

# ----- [NEW-4] ssGSEA helpers -----
term2gene_to_list <- function(term2gene_df) {
  split(toupper(term2gene_df$gene), term2gene_df$term)
}

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

# ----- Pseudobulk disambiguation helpers -----
find_ambiguous_sample_tissue_keys <- function(meta_sub, resolved_cols = character(),
                                               candidate_cols = character()) {
  compare_cols <- unique(c(resolved_cols, candidate_cols))
  if (length(compare_cols) == 0) return(character())
  key_df <- unique(meta_sub[, unique(c(SAMPLE_COL, TISSUE_COL, compare_cols)), drop = FALSE])
  key_df[] <- lapply(key_df, normalize_meta_values)
  resolved_key <- do.call(paste, c(key_df[, unique(c(SAMPLE_COL, TISSUE_COL, resolved_cols)), drop = FALSE], sep = "__"))
  base_key     <- paste(key_df[[SAMPLE_COL]], key_df[[TISSUE_COL]], sep = "__")
  unique(base_key[duplicated(resolved_key) | duplicated(resolved_key, fromLast = TRUE)])
}

resolve_pseudobulk_group_columns <- function(meta_sub) {
  group_cols <- c(SAMPLE_COL, TISSUE_COL)
  candidate_cols <- setdiff(unique(c(if (!is.null(DESIGN_COVARIATE)) DESIGN_COVARIATE,
                                      PSEUDOBULK_DISAMBIGUATION_CANDIDATES)),
                             c(SAMPLE_COL, TISSUE_COL))
  candidate_cols <- candidate_cols[candidate_cols %in% colnames(meta_sub)]
  ambiguous_keys <- find_ambiguous_sample_tissue_keys(meta_sub, candidate_cols = candidate_cols)
  disambiguation_cols <- character()

  while (length(ambiguous_keys) > 0 && length(candidate_cols) > 0) {
    informative_cols <- candidate_cols[vapply(candidate_cols, function(col) {
      probe_df  <- unique(meta_sub[, c(SAMPLE_COL, TISSUE_COL, col), drop = FALSE])
      probe_df[] <- lapply(probe_df, normalize_meta_values)
      probe_key <- paste(probe_df[[SAMPLE_COL]], probe_df[[TISSUE_COL]], sep = "__")
      probe_df  <- probe_df[probe_key %in% ambiguous_keys, , drop = FALSE]
      probe_key <- probe_key[probe_key %in% ambiguous_keys]
      if (nrow(probe_df) == 0) return(FALSE)
      any(vapply(split(probe_df[[col]], probe_key), function(v) length(unique(v)) > 1, logical(1)))
    }, logical(1))]
    if (length(informative_cols) == 0) break
    chosen_col         <- informative_cols[[1]]
    group_cols         <- c(group_cols, chosen_col)
    disambiguation_cols <- c(disambiguation_cols, chosen_col)
    candidate_cols     <- setdiff(candidate_cols, chosen_col)
    ambiguous_keys     <- find_ambiguous_sample_tissue_keys(
      meta_sub, resolved_cols = disambiguation_cols,
      candidate_cols = setdiff(candidate_cols, disambiguation_cols))
  }
  list(group_cols = unique(group_cols), disambiguation_cols = unique(disambiguation_cols),
       ambiguous_keys_remaining = ambiguous_keys)
}

# ----- Pseudobulk aggregation -----
aggregate_pseudobulk <- function(obj, group_value, group_col = CELLTYPE_L2_COL,
                                 min_total_cells = 50L) {
  cells <- colnames(obj)[obj@meta.data[[group_col]] == group_value]
  if (length(cells) < min_total_cells) return(NULL)

  sub        <- subset(obj, cells = cells)
  meta_sub   <- sub@meta.data
  counts_mat <- GetAssayData(sub, layer = "counts")

  resolved  <- resolve_pseudobulk_group_columns(meta_sub)
  if (length(resolved$ambiguous_keys_remaining) > 0) {
    stop(sprintf(
      "Non-unique sample+tissue keys for '%s'. Ambiguous: %s",
      group_value, paste(utils::head(resolved$ambiguous_keys_remaining, 10), collapse = ", ")))
  }

  group_cols     <- resolved$group_cols
  extra_cols     <- setdiff(group_cols, c(SAMPLE_COL, TISSUE_COL))
  meta_sub$pb_group <- apply(meta_sub[, group_cols, drop = FALSE], 1,
                              function(v) paste(normalize_meta_values(v), collapse = "__"))
  groups <- unique(meta_sub$pb_group)

  pb_cols <- list(); pb_meta <- data.frame(row.names = character(), stringsAsFactors = FALSE)

  for (g in groups) {
    cells_g <- rownames(meta_sub)[meta_sub$pb_group == g]
    if (length(cells_g) < MIN_CELLS_PER_PSEUDOBULK) next
    col_sum <- if (length(cells_g) == 1) counts_mat[, cells_g, drop = FALSE] else
      Matrix::rowSums(counts_mat[, cells_g, drop = FALSE])
    if (is.numeric(col_sum) && is.null(dim(col_sum))) {
      col_sum <- Matrix::Matrix(col_sum, ncol = 1, sparse = TRUE)
      rownames(col_sum) <- rownames(counts_mat)
    }
    colnames(col_sum) <- g; pb_cols[[g]] <- col_sum
    idx <- which(meta_sub$pb_group == g)[1]
    pb_meta[g, "tissue"]  <- meta_sub[[TISSUE_COL]][idx]
    pb_meta[g, "sample"]  <- meta_sub[[SAMPLE_COL]][idx]
    pb_meta[g, "n_cells"] <- length(cells_g)
    if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub))
      pb_meta[g, DESIGN_COVARIATE] <- meta_sub[[DESIGN_COVARIATE]][idx]
    for (ec in extra_cols) pb_meta[g, ec] <- meta_sub[[ec]][idx]
  }

  if (length(pb_cols) < 4) return(NULL)
  pb_counts <- as.matrix(do.call(cbind, pb_cols))
  valid <- colnames(pb_counts)[colSums(pb_counts) > 0 &
                                  !is.na(pb_meta[colnames(pb_counts), "tissue"])]
  if (length(valid) < 4) return(NULL)
  list(counts = pb_counts[, valid, drop = FALSE], meta = pb_meta[valid, , drop = FALSE],
       group_cols = group_cols, disambiguation_cols = resolved$disambiguation_cols)
}

# ----- DESeq2 pairwise -----
run_deseq2_pairwise <- function(pb, t1, t2) {
  keep <- pb$meta$tissue %in% c(t1, t2)
  if (sum(keep) < 4) return(NULL)
  counts_sub <- pb$counts[, keep, drop = FALSE]
  meta_sub   <- pb$meta[keep, , drop = FALSE]

  if (FILTER_IG_GENES) {
    keep_ig <- !is_ig_gene(rownames(counts_sub))
    if (!any(keep_ig)) {
      cat("    [WARN] All genes removed by IG filter before DESeq2\n")
      return(NULL)
    }
    counts_sub <- counts_sub[keep_ig, , drop = FALSE]
  }

  tls        <- make.names(c(t1, t2), unique = TRUE)
  meta_sub$tissue      <- droplevels(factor(meta_sub$tissue, levels = c(t1, t2)))
  meta_sub$tissue_safe <- factor(ifelse(meta_sub$tissue == t1, tls[1], tls[2]), levels = tls)
  n1 <- sum(meta_sub$tissue == t1); n2 <- sum(meta_sub$tissue == t2)
  if (n1 < MIN_SAMPLES_PER_TISSUE || n2 < MIN_SAMPLES_PER_TISSUE) return(NULL)
  keep_g <- rowSums(counts_sub >= 1) >= max(3, ncol(counts_sub) * 0.2)
  counts_sub <- counts_sub[keep_g, , drop = FALSE]
  if (nrow(counts_sub) < 100) return(NULL)
  tested_genes <- rownames(counts_sub)
  counts_int   <- round(counts_sub); storage.mode(counts_int) <- "integer"

  design_formula <- ~tissue_safe
  if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
    nlev <- length(unique(meta_sub[[DESIGN_COVARIATE]]))
    if (nlev >= 2 && nlev < nrow(meta_sub)) {
      meta_sub[[DESIGN_COVARIATE]] <- factor(meta_sub[[DESIGN_COVARIATE]])
      design_formula <- as.formula(paste("~", DESIGN_COVARIATE, "+ tissue_safe"))
      cat(sprintf("    Design: %s\n", deparse(design_formula)))
    }
  }

  dds <- tryCatch({
    d <- DESeqDataSetFromMatrix(counts_int, meta_sub, design_formula)
    DESeq(d, quiet = TRUE)
  }, error = function(e) {
    if (!is.null(DESIGN_COVARIATE)) {
      cat(sprintf("    [WARN] Design failed (%s), fallback ~ tissue_safe\n", e$message))
      tryCatch({ d2 <- DESeqDataSetFromMatrix(counts_int, meta_sub, ~tissue_safe); DESeq(d2, quiet = TRUE) },
               error = function(e2) { cat(sprintf("    [ERROR] DESeq2: %s\n", e2$message)); NULL })
    } else { cat(sprintf("    [ERROR] DESeq2: %s\n", e$message)); NULL }
  })
  if (is.null(dds)) return(NULL)

  res <- results(dds, contrast = c("tissue_safe", tls[2], tls[1]), alpha = PADJ_THR)
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>% arrange(padj) %>%
    mutate(sig       = ifelse(padj < PADJ_THR & abs(log2FoldChange) > LFC_THR, "sig", "ns"),
           direction = ifelse(log2FoldChange > 0, "up", "down"))
  list(de_table = res_df, tested_genes = tested_genes,
       n_up = sum(res_df$sig == "sig" & res_df$direction == "up"),
       n_down = sum(res_df$sig == "sig" & res_df$direction == "down"),
       n_samples_1 = n1, n_samples_2 = n2, tissue_1 = t1, tissue_2 = t2,
       design = deparse(design_formula))
}

# ----- Exploratory Wilcoxon -----
run_wilcox_exploratory <- function(obj, group_value, group_col = CELLTYPE_L2_COL) {
  cells <- colnames(obj)[obj@meta.data[[group_col]] == group_value]
  if (length(cells) < WILCOX_MIN_CELLS) return(NULL)
  sub     <- subset(obj, cells = cells)
  tissue_vec <- as.character(sub@meta.data[[TISSUE_COL]])
  tissues <- sort(unique(stats::na.omit(tissue_vec)))
  if (length(tissues) < 2) return(NULL)
  Idents(sub) <- TISSUE_COL; results <- list()
  for (pair in combn(tissues, 2, simplify = FALSE)) {
    t1 <- pair[1]; t2 <- pair[2]
    n1 <- sum(tissue_vec == t1, na.rm = TRUE)
    n2 <- sum(tissue_vec == t2, na.rm = TRUE)
    if (n1 < WILCOX_MIN_CELLS || n2 < WILCOX_MIN_CELLS) next
    cells_pair <- rownames(sub@meta.data)[tissue_vec %in% c(t1, t2)]
    sub_pair   <- subset(sub, cells = cells_pair); Idents(sub_pair) <- TISSUE_COL
    tissue_vec_pair <- as.character(sub_pair@meta.data[[TISSUE_COL]])
    if (n1 > WILCOX_MAX_CELLS_PER_IDENT || n2 > WILCOX_MAX_CELLS_PER_IDENT) {
      set.seed(42)
      k1 <- sample(rownames(sub_pair@meta.data)[tissue_vec_pair == t1],
                   min(n1, WILCOX_MAX_CELLS_PER_IDENT))
      k2 <- sample(rownames(sub_pair@meta.data)[tissue_vec_pair == t2],
                   min(n2, WILCOX_MAX_CELLS_PER_IDENT))
      sub_pair <- subset(sub_pair, cells = c(k1, k2)); Idents(sub_pair) <- TISSUE_COL
    }
    de <- tryCatch(FindMarkers(sub_pair, ident.1 = t2, ident.2 = t1, test.use = "wilcox",
                               min.pct = 0.1, logfc.threshold = 0), error = function(e) NULL)
    if (!is.null(de) && nrow(de) > 0) {
      de$gene <- rownames(de); de$padj <- p.adjust(de$p_val, method = "BH")
      de <- filter_marker_table_by_gene(
        de,
        gene_col = "gene",
        context = sprintf("Wilcoxon %s_vs_%s (%s)", t2, t1, group_value)
      )
      if (nrow(de) == 0) next
      results[[paste0(t2, "_vs_", t1)]] <- de %>%
        filter(abs(avg_log2FC) >= WILCOX_LFC_THR) %>%
        mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
        arrange(direction, padj, desc(abs(avg_log2FC))) %>%
        group_by(direction) %>% slice_head(n = WILCOX_TOP_N_PER_DIRECTION) %>% ungroup()
    }
  }
  results
}

# ----- LLM helper functions -----
sanitize_utf8_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = "byte"))
  bad <- is.na(out)
  if (any(bad)) out[bad] <- enc2utf8(x[bad])
  out
}

safe_trim <- function(x) {
  x <- to_scalar(x)
  x <- sanitize_utf8_text(x)
  trimws(x)
}

capture_object_text <- function(x) {
  if (is.null(x)) return("")
  txt <- tryCatch(paste(capture.output(print(x)), collapse = "\n"), error = function(e) "")
  if (!nzchar(trimws(txt)))
    txt <- tryCatch(paste(capture.output(str(x, max.level = 3)), collapse = "\n"), error = function(e) "")
  trimws(txt)
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
  signals <- sanitize_utf8_text(c(warnings, error_message))
  if (length(signals) == 0) return(FALSE)
  any(grepl("Failed to parse JSON response|premature EOF|invalid for atomic vectors|429|5[0-9]{2}|timeout|temporar",
            signals, ignore.case = TRUE, perl = TRUE))
}

has_stringdb_timeout_issue <- function(warnings = character(), error_message = NULL) {
  signals <- sanitize_utf8_text(c(warnings, error_message))
  if (length(signals) == 0) return(FALSE)
  any(grepl("stringdb-static\\.org|species\\.v12\\.txt", signals, ignore.case = TRUE, perl = TRUE))
}

looks_structured_interpret_agent_result <- function(x) {
  core <- unwrap_interpret_agent_result(x)
  any(nzchar(c(extract_named_text(core, c("overview", "summary", "interpretation", "narrative")),
               extract_named_text(core, c("key_mechanisms", "mechanisms", "keyMechanisms")),
               extract_named_text(core, c("hypothesis", "model", "working_hypothesis")),
               extract_named_text(core, c("key_drivers", "drivers", "genes", "gene_drivers")))))
}

placeholder_text <- function(x, default = "Not available from current evidence.") {
  x <- safe_trim(x); if (nzchar(x)) x else default
}

collapse_driver_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- trimws(as.character(x)); x <- x[nzchar(x)]; paste(unique(x), collapse = ", ")
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
  top_t <- utils::head(enrich_df[order(enrich_df$p.adjust), , drop = FALSE], n_terms)
  term_lines <- vapply(seq_len(nrow(top_t)), function(i)
    sprintf("- %s | padj=%s | Count=%s | Genes=%s",
            safe_trim(top_t$Description[i]),
            format(top_t$p.adjust[i], scientific = TRUE, digits = 3),
            safe_trim(top_t$Count[i]),
            truncate_text(gsub("/", ", ", safe_trim(top_t$geneID[i])), 200)),
    character(1))
  gene_lines <- character()
  if (!is.null(gene_fc) && length(gene_fc) > 0) {
    gene_fc <- sort(gene_fc[!is.na(gene_fc)], decreasing = TRUE)
    top_up   <- utils::head(gene_fc, n_genes)
    top_down <- utils::head(sort(gene_fc, decreasing = FALSE), n_genes)
    gene_lines <- c(sprintf("Top up genes: %s",
                             paste(sprintf("%s(%.2f)", names(top_up), top_up), collapse = ", ")),
                    sprintf("Top down genes: %s",
                             paste(sprintf("%s(%.2f)", names(top_down), top_down), collapse = ", ")))
  }
  paste(c("Top enrichment evidence:", term_lines, gene_lines), collapse = "\n")
}

select_interpret_agent_databases <- function(enr_list) {
  if (is.null(enr_list) || length(enr_list) == 0) return(character())
  db_order <- c(INTERPRET_AGENT_DB_PRIORITY, setdiff(names(enr_list), INTERPRET_AGENT_DB_PRIORITY))
  db_order <- db_order[db_order %in% names(enr_list)]
  db_order[vapply(db_order, function(db) {
    er <- enr_list[[db]]; if (is.null(er)) return(FALSE)
    er_df <- tryCatch(as.data.frame(er), error = function(e) NULL)
    !is.null(er_df) && nrow(er_df) >= INTERPRET_AGENT_MIN_TERMS
  }, logical(1))]
}

build_multi_database_evidence_text <- function(enrich_list, gene_fc = NULL, n_terms_per_db = 5L) {
  if (is.null(enrich_list) || length(enrich_list) == 0) return("")
  sections <- lapply(names(enrich_list), function(db) {
    er_df <- tryCatch(as.data.frame(enrich_list[[db]]), error = function(e) NULL)
    if (is.null(er_df) || nrow(er_df) == 0) return(NULL)
    top_t <- utils::head(er_df[order(er_df$p.adjust), , drop = FALSE], n_terms_per_db)
    tl    <- vapply(seq_len(nrow(top_t)), function(i)
      sprintf("- %s | padj=%s | Count=%s | Genes=%s",
              safe_trim(top_t$Description[i]),
              format(top_t$p.adjust[i], scientific = TRUE, digits = 3),
              safe_trim(top_t$Count[i]),
              truncate_text(gsub("/", ", ", safe_trim(top_t$geneID[i])), 180)),
      character(1))
    c(sprintf("[%s]", db), tl)
  })
  sections <- Filter(Negate(is.null), sections)
  gene_lines <- character()
  if (!is.null(gene_fc) && length(gene_fc) > 0) {
    gene_fc <- gene_fc[!is.na(gene_fc)]
    top_abs <- names(utils::head(sort(abs(gene_fc), decreasing = TRUE), 20))
    gene_lines <- sprintf("Top |log2FC| genes: %s",
                           paste(sprintf("%s(%.2f)", top_abs, gene_fc[top_abs]), collapse = ", "))
  }
  paste(c("Integrated multi-database evidence:", unlist(sections), gene_lines), collapse = "\n")
}

build_multi_database_raw_text <- function(per_db_payloads) {
  if (length(per_db_payloads) == 0) return("")
  blocks <- lapply(names(per_db_payloads), function(db) {
    payload  <- per_db_payloads[[db]]
    raw_text <- capture_object_text(payload$result)
    c(sprintf("=== %s ===", db),
      sprintf("error: %s", ifelse(is.null(payload$error) || !nzchar(payload$error), "", payload$error)),
      sprintf("warnings: %s", paste(payload$warnings, collapse = " | ")),
      if (nzchar(raw_text)) raw_text else "<empty interpret_agent output>")
  })
  paste(unlist(blocks), collapse = "\n")
}

run_interpret_agent_multi <- function(enr_list, context_str, dir_name, gene_fc = NULL) {
  selected_dbs <- select_interpret_agent_databases(enr_list)
  if (length(selected_dbs) == 0) {
    return(list(selected_dbs = character(), per_db_payloads = list(), combined_raw_text = "",
                combined_warnings = "", combined_error = "no enrichment database met interpret_agent criteria"))
  }
  per_db_payloads <- list(); combined_warnings <- character(); combined_errors <- character()
  for (db in selected_dbs) {
    cat(sprintf("    interpret_agent [%s] (db: %s)\n", dir_name, db))
    db_ctx  <- paste(context_str,
                     sprintf("Primary enrichment database for this pass: %s.", db),
                     sprintf("Integrate with broader multi-database evidence, but prioritize terms from %s.", db),
                     sep = "\n\n")
    payload <- run_interpret_agent_safe(enr_list[[db]], db_ctx, gene_fc)
    per_db_payloads[[db]] <- payload
    if (length(payload$warnings) > 0)
      combined_warnings <- c(combined_warnings, sprintf("%s: %s", db, paste(payload$warnings, collapse = " | ")))
    if (!is.null(payload$error) && nzchar(payload$error))
      combined_errors <- c(combined_errors, sprintf("%s: %s", db, payload$error))
  }
  list(selected_dbs = selected_dbs, per_db_payloads = per_db_payloads,
       combined_raw_text = build_multi_database_raw_text(per_db_payloads),
       combined_warnings = unique(combined_warnings),
       combined_error    = paste(unique(combined_errors), collapse = " | "))
}

build_multi_database_interpretation_record <- function(multi_payload, enrich_list,
                                                        celltype_label, celltype_level,
                                                        comp_name, dir_name, context_str,
                                                        gene_fc = NULL) {
  selected_dbs     <- multi_payload$selected_dbs
  combined_raw_text <- safe_trim(multi_payload$combined_raw_text)
  combined_warnings <- unique(multi_payload$combined_warnings)
  combined_error    <- safe_trim(multi_payload$combined_error)
  source_db         <- paste(selected_dbs, collapse = "; ")
  evidence_text     <- build_multi_database_evidence_text(enrich_list[selected_dbs], gene_fc = gene_fc)

  std_payload <- standardize_result_with_llm(
    raw_text = combined_raw_text,
    context_str = paste(context_str,
                        sprintf("Integrate evidence across ALL eligible databases: %s.", source_db),
                        "When databases agree, state the consensus clearly. When they differ, explain what each contributes.",
                        sep = "\n\n"),
    evidence_text = evidence_text,
    celltype_label = celltype_label,
    celltype_level = celltype_level,
    comp_name = comp_name,
    dir_name = dir_name, source_db = source_db, warnings = combined_warnings,
    error_message = if (nzchar(combined_error)) combined_error else NULL
  )

  std    <- std_payload$result
  status <- if (!is.null(std)) "structured" else if (nzchar(combined_raw_text)) "raw_text" else
    if (nzchar(combined_error)) "error" else "empty"

    list(celltype_level = celltype_level,
      celltype_label = celltype_label,
      celltype_l2 = celltype_label,
      comparison = comp_name, direction = dir_name,
       source_db = source_db, status = status,
       warnings = unique(c(combined_warnings, std_payload$warnings)),
       error = ifelse(!is.null(std_payload$error) && nzchar(std_payload$error),
                      std_payload$error, combined_error),
       overview       = placeholder_text(if (!is.null(std)) std$overview       else ""),
       key_mechanisms = placeholder_text(if (!is.null(std)) std$key_mechanisms else ""),
       hypothesis     = placeholder_text(if (!is.null(std)) std$hypothesis     else ""),
       narrative      = placeholder_text(if (!is.null(std)) std$narrative      else combined_raw_text),
       key_drivers    = placeholder_text(if (!is.null(std)) collapse_driver_field(std$key_drivers) else ""),
       evidence       = placeholder_text(if (!is.null(std)) std$evidence       else evidence_text),
       limitations    = placeholder_text(if (!is.null(std)) std$limitations    else
         "Interpretation could not be standardized across databases."),
       raw_text = combined_raw_text, raw_result = multi_payload)
}

extract_json_string <- function(text) {
  text <- safe_trim(text); if (!nzchar(text)) return("")
  text <- gsub("^```(?:json)?\\s*", "", text, perl = TRUE)
  text <- gsub("\\s*```$", "", text, perl = TRUE)
  si   <- regexpr("\\{", text, perl = TRUE)[1]
  ep   <- gregexpr("\\}", text, perl = TRUE)[[1]]
  if (si < 1 || length(ep) == 0 || ep[1] < 1) return("")
  substr(text, si, ep[length(ep)])
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

standardize_result_with_llm <- function(raw_text, context_str, evidence_text,
                                         celltype_label, celltype_level,
                                         comp_name, dir_name, source_db,
                                         warnings = character(), error_message = NULL) {
  if (!ENABLE_LLM || !STANDARDIZE_LLM_OUTPUT)
    return(list(result = NULL, warnings = character(), error = "standardization disabled"))

  prompt <- paste(
    "You are standardizing a biological interpretation for single-cell epithelial tissue comparison.",
    "Return valid JSON only. No markdown, no code fences, no extra commentary.",
    "Use exactly these keys: overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "",
    "Rules:",
    "1. overview/key_mechanisms/hypothesis/narrative/evidence/limitations must be strings.",
    "2. key_drivers must be an array of short gene or regulator names.",
    "3. Every field must be present.",
    "4. If evidence is weak, stay conservative and say so in limitations.",
    "5. If a field cannot be recovered, use 'Not available from current evidence.'.",
    "6. Write ALL narrative content in Simplified Chinese.",
    "7. In key_mechanisms and key_drivers, prioritize genes with the largest absolute log2FoldChange values.",
    "   For each top gene, explicitly state: gene, log2FC, pathway(s), and why the relation matters biologically.",
    "   Use the format: gene(log2FC=X.X) -> pathway: explanation.",
    "",
    sprintf("Cell type level: %s", celltype_level),
    sprintf("Cell type group: %s", celltype_label),
    sprintf("Comparison: %s", comp_name),
    sprintf("Direction: %s", dir_name),  sprintf("Source DB: %s", source_db),
    "Context:", truncate_text(context_str),
    "Original interpret_agent output:", truncate_text(raw_text),
    "Enrichment evidence:", truncate_text(evidence_text),
    "Warnings and errors:", truncate_text(paste(c(warnings, error_message), collapse = "\n")),
    sep = "\n"
  )

  collected_warnings <- character(); last_error <- NULL
  for (attempt in seq_len(STANDARDIZE_LLM_MAX_RETRIES)) {
    response_text <- tryCatch(
      safe_trim(fanyi::chat_request(prompt, model = STANDARDIZE_LLM_MODEL, api_key = DEEPSEEK_API_KEY)),
      error = function(e) { last_error <<- conditionMessage(e); "" })
    parsed <- parse_standardized_json(response_text)
    if (!is.null(parsed))
      return(list(result = parsed, warnings = unique(collected_warnings), error = NULL))
    if (nzchar(response_text))
      collected_warnings <- c(collected_warnings, sprintf("standardize attempt %d returned non-JSON", attempt))
    if (attempt < STANDARDIZE_LLM_MAX_RETRIES) Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
  }
  list(result = NULL, warnings = unique(collected_warnings),
       error = ifelse(is.null(last_error) || !nzchar(last_error),
                      "failed to standardize interpret_agent output", last_error))
}

flatten_interpretation_records <- function(x) {
  rows <- list(); idx <- 1L
  for (ct_l2 in names(x)) for (comp_name in names(x[[ct_l2]])) for (dir_name in names(x[[ct_l2]][[comp_name]])) {
    rec <- x[[ct_l2]][[comp_name]][[dir_name]]; if (is.null(rec)) next
    rows[[idx]] <- data.frame(
      celltype_level = safe_trim(rec$celltype_level),
      celltype_label = safe_trim(rec$celltype_label),
      celltype_l2 = safe_trim(rec$celltype_l2),
      comparison = safe_trim(rec$comparison),
      direction = safe_trim(rec$direction), source_db = safe_trim(rec$source_db),
      status = safe_trim(rec$status), warnings = paste(rec$warnings, collapse = " | "),
      error = safe_trim(rec$error), overview = safe_trim(rec$overview),
      key_mechanisms = safe_trim(rec$key_mechanisms), hypothesis = safe_trim(rec$hypothesis),
      narrative = safe_trim(rec$narrative), key_drivers = safe_trim(rec$key_drivers),
      evidence = safe_trim(rec$evidence), limitations = safe_trim(rec$limitations),
      raw_text = safe_trim(rec$raw_text), stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}

write_interpretation_markdown <- function(records, path, title, include_raw = FALSE) {
  md <- c(title, "", sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
  if (nrow(records) == 0) { md <- c(md, "No interpret_agent records available."); writeLines(md, path); return(invisible(path)) }
  for (i in seq_len(nrow(records))) {
    rec <- records[i, , drop = FALSE]
    md  <- c(md, sprintf("## %s | %s | %s | %s",
                         rec$celltype_level, rec$celltype_label,
                         rec$comparison, rec$direction), "")
    md  <- c(md, sprintf("- **Status:** %s", rec$status),
             sprintf("- **Cell Type Level:** %s", rec$celltype_level),
             sprintf("- **Source DB:** %s", ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
    if (nzchar(rec$warnings)) md <- c(md, sprintf("- **Warnings:** %s", rec$warnings))
    if (nzchar(rec$error))    md <- c(md, sprintf("- **Error:** %s",    rec$error))
    md <- c(md, "")
    for (nm in c("Overview", "Key Mechanisms", "Hypothesis", "Narrative", "Key Drivers", "Evidence", "Limitations")) {
      field_name <- tolower(gsub(" ", "_", nm))
      field_name <- gsub("key_mechanisms", "key_mechanisms", gsub("key_drivers", "key_drivers", field_name))
      val <- placeholder_text(rec[[field_name]])
      md  <- c(md, sprintf("### %s", nm), "", val, "")
    }
    if (include_raw && nzchar(rec$raw_text)) md <- c(md, "### Raw Output", "", rec$raw_text, "")
  }
  writeLines(md, path); invisible(path)
}

run_pairwise_tissue_comparison <- function(obj, group_col, group_values, level_name,
                                           out_dir) {
  pb_de_level            <- list()
  enrich_level           <- list()
  agent_level            <- list()
  agent_structured_level <- list()

  cat(sprintf("=== Pseudobulk DESeq2 Tissue Comparison [%s] ===\n", level_name))
  cat(sprintf("Thresholds: padj < %s, |log2FC| > %s\n", PADJ_THR, LFC_THR))
  cat("Enrichment databases: GO_BP, GO_MF, GO_CC, KEGG, Hallmark, CellMarker, PanglaoDB, Epithelial_custom\n\n")

  for (group_name in group_values) {
    cat(sprintf("\n>> %s: %s\n", level_name, group_name))

    pb <- aggregate_pseudobulk(obj, group_name, group_col = group_col)
    if (is.null(pb)) { cat("  [SKIP] Insufficient pseudobulk\n"); next }

    pb_tissues <- unique(na.omit(pb$meta$tissue))
    if (length(pb_tissues) < 2) { cat("  [SKIP] < 2 tissues\n"); next }

    cat(sprintf("  Pseudobulk samples: %d (%s)\n",
                nrow(pb$meta), paste(pb_tissues, collapse = ", ")))
    cat(sprintf("  Grouping columns: %s\n", paste(pb$group_cols, collapse = ", ")))

    pb_de_level[[group_name]]            <- list()
    enrich_level[[group_name]]           <- list()
    agent_level[[group_name]]            <- list()
    agent_structured_level[[group_name]] <- list()

    for (pair in combn(as.character(pb_tissues), 2, simplify = FALSE)) {
      t1 <- pair[1]
      t2 <- pair[2]
      comp_name <- paste0(t2, "_vs_", t1)
      cat(sprintf("  DESeq2: %s\n", comp_name))

      res <- run_deseq2_pairwise(pb, t1, t2)
      if (is.null(res)) { cat("    [SKIP] Too few samples per tissue\n"); next }

      pb_de_level[[group_name]][[comp_name]] <- res
      cat(sprintf("    DEGs: %d up, %d down (samples: %d vs %d)\n",
                  res$n_up, res$n_down, res$n_samples_1, res$n_samples_2))

      comp_dir <- file.path(out_dir, safe_name(group_name), safe_name(comp_name))
      dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
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

        if (length(genes) < 5) {
          enrich_level[[group_name]][[comp_name]][[dir_name]] <- NULL
          next
        }
        cat(sprintf("    Enrichment [%s]: %d genes\n", dir_name, length(genes)))

        tested <- res$tested_genes
        enr_list <- list(
          GO_BP             = run_gmt_enrichment(genes, go_bp_t2g,             "GO_BP",             tested),
          GO_MF             = run_gmt_enrichment(genes, go_mf_t2g,             "GO_MF",             tested),
          GO_CC             = run_gmt_enrichment(genes, go_cc_t2g,             "GO_CC",             tested),
          KEGG              = run_gmt_enrichment(genes, kegg_t2g,              "KEGG",              tested),
          Hallmark          = run_gmt_enrichment(genes, hallmark_t2g,          "Hallmark",          tested),
          CellMarker        = run_gmt_enrichment(genes, cellmarker_t2g,        "CellMarker",        tested),
          PanglaoDB         = run_gmt_enrichment(genes, panglaodb_t2g,         "PanglaoDB",         tested),
          Epithelial_custom = run_gmt_enrichment(genes, epithelial_custom_t2g, "Epithelial_custom", tested)
        )
        enr_list <- enr_list[!sapply(enr_list, is.null)]
        enrich_level[[group_name]][[comp_name]][[dir_name]] <- enr_list

        enr_out <- file.path(comp_dir, paste0("enrichment_", dir_name))
        dir.create(enr_out, showWarnings = FALSE)
        for (db in names(enr_list)) {
          er <- enr_list[[db]]
          if (nrow(as.data.frame(er)) > 0) {
            fwrite(as.data.frame(er), file.path(enr_out, paste0(db, ".csv")))
            tryCatch({
              pdf(file.path(enr_out, paste0(db, "_dotplot.pdf")), width = 10, height = 8)
              print(dotplot(er, showCategory = 15,
                            title = sprintf("%s %s (%s %s %s)",
                                            db, dir_name, level_name, group_name, comp_name)))
              dev.off()
            }, error = function(e) NULL)
          }
        }
        saveRDS(enr_list, file.path(enr_out, "all_enrichment.rds"))

        if (ENABLE_LLM) {
          tissue_ctx   <- build_tissue_pair_context(
            t1, t2, group_name, dir_name, celltype_level = level_name)
          sig_de       <- res$de_table %>% filter(sig == "sig")
          gene_fc      <- setNames(sig_de$log2FoldChange, toupper(sig_de$gene))
          selected_dbs <- select_interpret_agent_databases(enr_list)

          if (length(selected_dbs) > 0) {
            cat(sprintf("    interpret_agent [%s] integrating %d db(s): %s\n",
                        dir_name, length(selected_dbs), paste(selected_dbs, collapse = ", ")))
            multi_payload <- run_interpret_agent_multi(enr_list[selected_dbs], tissue_ctx, dir_name, gene_fc)
            agent_level[[group_name]][[comp_name]][[dir_name]] <- multi_payload
            ia_structured <- build_multi_database_interpretation_record(
              multi_payload = multi_payload,
              enrich_list = enr_list,
              celltype_label = group_name,
              celltype_level = level_name,
              comp_name = comp_name,
              dir_name = dir_name,
              context_str = tissue_ctx,
              gene_fc = gene_fc
            )
            agent_structured_level[[group_name]][[comp_name]][[dir_name]] <- ia_structured
            saveRDS(multi_payload, file.path(comp_dir, paste0("interpret_agent_", dir_name, "_payload.rds")))
            saveRDS(ia_structured, file.path(comp_dir, paste0("interpret_agent_", dir_name, "_structured.rds")))
            for (db in names(multi_payload$per_db_payloads)) {
              saveRDS(multi_payload$per_db_payloads[[db]],
                      file.path(comp_dir, paste0("interpret_agent_", dir_name, "_", safe_name(db), "_payload.rds")))
            }
            if (nzchar(ia_structured$raw_text)) {
              writeLines(ia_structured$raw_text,
                         file.path(comp_dir, paste0("interpret_agent_", dir_name, "_raw.txt")))
            }
          } else {
            cat(sprintf("    [SKIP] interpret_agent [%s]: no database met minimum term count\n", dir_name))
          }
        }
      }
    }
  }

  list(
    pb_de = pb_de_level,
    enrich = enrich_level,
    agent = agent_level,
    agent_structured = agent_structured_level
  )
}

run_grouped_ssgsea <- function(obj, group_col, level_name, methods,
                               n_top, n_heatmap) {
  results_all <- list()
  group_key   <- paste0("ssgsea_group_", tolower(level_name))
  file_stem   <- level_file_stem("ssgsea", level_name)

  obj@meta.data[[group_key]] <- paste0(
    obj@meta.data[[TISSUE_COL]], "__", obj@meta.data[[group_col]]
  )
  ssgsea_groups <- sort(unique(obj@meta.data[[group_key]]))

  cat(sprintf("\n=== Pseudobulk ssGSEA (tissue x %s) ===\n", level_name))
  cat(sprintf("Methods: %s\n\n", paste(methods, collapse = ", ")))
  cat(sprintf("[INFO] ssGSEA groups (%d): %s\n\n",
              length(ssgsea_groups), paste(ssgsea_groups, collapse = ", ")))

  for (method in methods) {
    cat(sprintf("--- ssGSEA method [%s]: %s ---\n", level_name, method))

    selected_t2g <- switch(
      method,
      hallmark = hallmark_t2g, go_bp = go_bp_t2g, go_mf = go_mf_t2g,
      go_cc = go_cc_t2g, kegg = kegg_t2g,
      custom = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT)) {
        tryCatch(read.gmt(SSGSEA_CUSTOM_GMT) %>% mutate(gene = toupper(gene)),
                 error = function(e) { cat("[WARN] custom GMT failed\n"); NULL })
      } else NULL,
      NULL
    )

    if (is.null(selected_t2g) || nrow(selected_t2g) == 0) {
      cat(sprintf("[WARN] Method '%s': no gene sets, skipping\n", method)); next
    }

    gs_list <- term2gene_to_list(selected_t2g)
    gs_list <- map_gene_sets_to_features(gs_list, rownames(obj))
    gs_list <- filter_gene_sets_for_analysis(gs_list)
    gs_list <- filter_gs_size(gs_list)
    if (length(gs_list) == 0) {
      cat(sprintf("[WARN] Method '%s': no valid gene sets after filtering\n", method)); next
    }
    cat(sprintf("[OK] %d gene sets prepared\n", length(gs_list)))

    needed_genes <- unique(unlist(gs_list, use.names = FALSE))
    cat(sprintf("[INFO] Computing pseudobulk for %d genes...\n", length(needed_genes)))
    avg_expr <- tryCatch(
      AverageExpression(obj, assays = "RNA", slot = "data",
                        group.by = group_key, features = needed_genes,
                        verbose = FALSE)[["RNA"]],
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
        GSVA::gsva(as.matrix(avg_expr), gs_list, method = "ssgsea",
                   ssgsea.norm = TRUE, verbose = FALSE)
      }
    }, error = function(e) { cat(sprintf("[ERROR] ssGSEA failed: %s\n", e$message)); NULL })
    if (is.null(ssgsea_scores)) { cat("[WARN] Skipping this method\n"); next }

    ssgsea_z <- t(scale(t(ssgsea_scores)))
    cat(sprintf("[OK] ssGSEA done: %d pathways x %d groups\n",
                nrow(ssgsea_scores), ncol(ssgsea_scores)))

    results_all[[method]] <- list(scores = ssgsea_scores, z_scores = ssgsea_z, gene_sets = gs_list)
    saveRDS(ssgsea_scores, file.path(rpt_dir, sprintf("%s_scores_%s.rds", file_stem, method)))
    saveRDS(ssgsea_z,      file.path(rpt_dir, sprintf("%s_z_%s.rds",     file_stem, method)))

    top_rows <- lapply(colnames(ssgsea_scores), function(gid) {
      sc  <- ssgsea_scores[, gid]
      z   <- ssgsea_z[, gid]
      ord <- order(z, sc, decreasing = TRUE, na.last = TRUE)
      ti  <- ord[seq_len(min(n_top, length(ord)))]
      data.frame(
        level = level_name,
        method = method,
        group = gid,
        pathway = rownames(ssgsea_scores)[ti],
        score = sc[ti],
        z_score = z[ti],
        rank = seq_along(ti),
        stringsAsFactors = FALSE
      )
    })
    fwrite(bind_rows(top_rows), file.path(rpt_dir, sprintf("%s_top_pathways_%s.csv", file_stem, method)))

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
          main   = sprintf("Epithelial ssGSEA Z-score (%s) — top %d pathways by tissue x %s",
                           method, nrow(heat_mat), level_name),
          fontsize_row = 7, fontsize_col = 9, cellwidth = 22, cellheight = 10,
          filename = paste0(heat_path, ".pdf"),
          width  = max(10, ncol(heat_mat) * 1.8),
          height = max(8,  nrow(heat_mat) * 0.35 + 3)
        )
        grDevices::png(paste0(heat_path, ".png"),
                       width  = max(10, ncol(heat_mat) * 1.8),
                       height = max(8,  nrow(heat_mat) * 0.35 + 3),
                       units = "in", res = 300)
        grid::grid.newpage(); grid::grid.draw(ht$gtable); grDevices::dev.off()
        cat(sprintf("[OK] Heatmap saved: %s_heatmap_%s\n", file_stem, method))
      }, error = function(e) cat(sprintf("[WARN] Heatmap failed for %s: %s\n", method, e$message)))
    }

    cat(sprintf("[OK] Method '%s' complete\n\n", method))
  }

  cat(sprintf("[SUMMARY] Completed %d ssGSEA methods for %s\n", length(results_all), level_name))
  for (mn in names(results_all)) {
    cat(sprintf("  - %s: %d pathways x %d groups\n",
                mn, nrow(results_all[[mn]]$scores), ncol(results_all[[mn]]$scores)))
  }
  cat("\n")

  results_all
}

interpret_agent_runtime <- new.env(parent = emptyenv())
interpret_agent_runtime$add_ppi <- INTERPRET_AGENT_ADD_PPI

run_interpret_agent_once <- function(enrich_obj, context_str, gene_fc = NULL,
                                      add_ppi = interpret_agent_runtime$add_ppi) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0)
    return(list(result = NULL, warnings = character(), error = "empty enrichment"))
  warn_msgs <- character()
  res <- tryCatch(
    withCallingHandlers(
      clusterProfiler::interpret_agent(x = enrich_obj, context = context_str,
                                        n_pathways = INTERPRET_AGENT_N_PATHWAYS,
                                        model = INTERPRET_AGENT_MODEL, api_key = DEEPSEEK_API_KEY,
                                        add_ppi = add_ppi, gene_fold_change = gene_fc),
      warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
    ),
    error = function(e) { cat(sprintf("    [ERROR] agent: %s\n", e$message)); structure(list(message = e$message), class = "interpret_agent_error") }
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
    payload      <- run_interpret_agent_once(enrich_obj, context_str, gene_fc, interpret_agent_runtime$add_ppi)
    last_payload <- payload
    if (isTRUE(INTERPRET_AGENT_DISABLE_PPI_ON_STRINGDB_TIMEOUT) &&
        isTRUE(interpret_agent_runtime$add_ppi) &&
        has_stringdb_timeout_issue(payload$warnings, payload$error)) {
      interpret_agent_runtime$add_ppi <- FALSE
      combined_warnings <- c(combined_warnings, "STRINGdb timed out; disabling add_ppi for remaining calls")
      cat("    [WARN] STRINGdb timed out; disabling add_ppi\n")
    }
    if (length(payload$warnings) > 0)
      combined_warnings <- c(combined_warnings, sprintf("attempt %d: %s", attempt, payload$warnings))
    retryable      <- is_retryable_interpret_agent_issue(payload$warnings, payload$error)
    structured     <- !is.null(payload$result) && looks_structured_interpret_agent_result(payload$result)
    raw_text_valid <- !is.null(payload$result) && nzchar(capture_object_text(payload$result))
    if (is.null(payload$error) && (structured || (raw_text_valid && !retryable))) {
      payload$warnings <- unique(c(combined_warnings, payload$warnings)); return(payload)
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
  if (is.list(x)) return(paste(sanitize_utf8_text(capture.output(str(x, max.level = 2))), collapse = "\n"))
  sanitize_utf8_text(as.character(x))
}

stratified_downsample <- function(obj, group_col, n_per = HEATMAP_CELLS_PER_TYPE) {
  set.seed(42)
  meta  <- obj@meta.data
  cells <- unlist(lapply(unique(meta[[group_col]]), function(g) {
    gc <- rownames(meta)[meta[[group_col]] == g]; sample(gc, min(length(gc), n_per))
  }))
  subset(obj, cells = cells)
}

save_plot <- function(p, path_no_ext, width = 10, height = 8) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
}

add_md_image_if_exists <- function(add_fun, output_dir, rel_path, alt_text) {
  if (file.exists(file.path(output_dir, rel_path)))
    add_fun(sprintf("![%s](%s)", alt_text, rel_path))
  else
    add_fun(sprintf("_Not generated: `%s`_", rel_path))
  add_fun("")
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
      if (length(x) == 0 || all(is.na(x))) return(NA_character_)
      paste(as.character(x), collapse = "; ")
    }, character(1))
    else if (is.character(vec))                      out[[col]] <- trimws(vec)
    else if (!(is.integer(vec) || is.numeric(vec)))  out[[col]] <- as.character(vec)
  }
  out
}

matrix_to_scipy_csr <- function(mat, scipy_sparse, np) {
  mat_t   <- Matrix::t(as(mat, "dgCMatrix"))
  mat_coo <- as.data.frame(Matrix::summary(mat_t))
  scipy_sparse$coo_matrix(
    reticulate::tuple(np$array(mat_coo$x, dtype = np$float32),
                      reticulate::tuple(np$array(mat_coo$i - 1L, dtype = np$int32),
                                        np$array(mat_coo$j - 1L, dtype = np$int32))),
    shape = reticulate::tuple(as.integer(nrow(mat_t)), as.integer(ncol(mat_t)))
  )$tocsr()
}

# ==============================================================================
# 5. Load Data
# ==============================================================================

cat("=== Loading Epithelial Data ===\n")
if (!file.exists(H5AD_PATH)) stop(sprintf("File not found: %s", H5AD_PATH))
obj <- GetSeurat(h5ad_path = H5AD_PATH, prefer_raw = FALSE,
                 prefer_layer_counts = TRUE, validate_counts = TRUE, debug = TRUE)
cat(sprintf("[OK] %d cells x %d genes\n", ncol(obj), nrow(obj)))
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA assay missing 'counts' layer.")
cat("[OK] counts layer verified\n\n")
if (FILTER_IG_GENES) {
  ig_genes_present <- rownames(obj)[is_ig_gene(rownames(obj))]
  cat(sprintf("[INFO] IG filter enabled: %d IG-related genes detected in feature space\n",
              length(ig_genes_present)))
}

fig_dir <- file.path(OUTPUT_DIR, "figures")
rpt_dir <- file.path(OUTPUT_DIR, "reports")
de_dir  <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir  <- file.path(OUTPUT_DIR, "wilcox_exploratory")
de_l3_dir <- file.path(OUTPUT_DIR, "pseudobulk_de_L3")
wx_l3_dir <- file.path(OUTPUT_DIR, "wilcox_exploratory_L3")
for (d in c(fig_dir, rpt_dir, de_dir, wx_dir, de_l3_dir, wx_l3_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. Validate Metadata & Apply L2 Remapping
# ==============================================================================

cat("=== Validating Metadata & L2 Remapping ===\n")
meta <- obj@meta.data
for (col in c(TISSUE_COL, SAMPLE_COL, L3_SOURCE_COL))
  if (!col %in% colnames(meta)) stop(sprintf("Missing required column: %s", col))

l3_vals  <- as.character(meta[[L3_SOURCE_COL]])
l2_vals  <- L3_TO_L2_REMAP[l3_vals]
unmapped <- unique(l3_vals[is.na(l2_vals)])
if (length(unmapped) > 0)
  stop(sprintf("Unmapped L3 values in '%s': %s\nUpdate L3_TO_L2_REMAP before running.",
               L3_SOURCE_COL, paste(unmapped, collapse = ", ")))

obj@meta.data[[CELLTYPE_L2_COL]] <- l2_vals
obj@meta.data[[CELLTYPE_L3_COL]] <- l3_vals

cat("\nL3 -> L2 Remapping Summary:\n")
print(table(L3 = l3_vals, L2 = l2_vals, useNA = "ifany"))
cat("\nL2 Distribution:\n"); print(table(obj@meta.data[[CELLTYPE_L2_COL]], useNA = "ifany")); cat("\n")
cat("L3 Distribution:\n"); print(table(obj@meta.data[[CELLTYPE_L3_COL]], useNA = "ifany")); cat("\n")

bad_idx <- is.na(meta[[TISSUE_COL]]) | trimws(as.character(meta[[TISSUE_COL]])) == "" |
           is.na(meta[[SAMPLE_COL]]) | trimws(as.character(meta[[SAMPLE_COL]])) == "" |
           is.na(obj@meta.data[[CELLTYPE_L2_COL]]) |
           is.na(meta[[LABEL_COL]])  | trimws(as.character(meta[[LABEL_COL]])) == ""
if (sum(bad_idx) > 0) {
  cat(sprintf("[INFO] Dropping %d cells with NA tissue/sample/L2/label\n", sum(bad_idx)))
  obj <- subset(obj, cells = colnames(obj)[!bad_idx])
}

n_samples <- length(unique(obj@meta.data[[SAMPLE_COL]]))
if (n_samples < 2) stop(sprintf("Only %d sample(s). Need >= 2 for pseudobulk DE.", n_samples))
cat(sprintf("[OK] %d unique samples\n", n_samples))

meta     <- obj@meta.data
tissues  <- sort(unique(na.omit(meta[[TISSUE_COL]])))
l2_types <- sort(unique(na.omit(meta[[CELLTYPE_L2_COL]])))
l3_types <- sort(unique(na.omit(meta[[CELLTYPE_L3_COL]])))

cat(sprintf("[OK] Cells: %d\n", ncol(obj)))
cat(sprintf("[OK] Tissues: %s\n", paste(tissues, collapse = ", ")))
cat(sprintf("[OK] L2 types (%d): %s\n", length(l2_types), paste(l2_types, collapse = ", ")))
cat(sprintf("[OK] L3 types (%d): %s\n", length(l3_types), paste(l3_types, collapse = ", ")))
cat(sprintf("[OK] L3 column: %s (%d unique)\n", LABEL_COL, length(unique(meta[[LABEL_COL]]))))
cat("\n[INFO] Running NormalizeData() for visualization and marker analyses\n")
obj <- NormalizeData(obj, verbose = FALSE)
cat("[OK] data layer ready\n\n")

# ==============================================================================
# 7. Visualization
# ==============================================================================

cat("=== Visualization ===\n")
umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
if (!is.null(umap_reduction)) {
  tissue_cols_use <- UMAP_TISSUE_COLORS[names(UMAP_TISSUE_COLORS) %in% tissues]
  p1 <- build_umap_plot(obj, umap_reduction, TISSUE_COL,
                         title = sprintf("Epithelial - Tissue (%s)", umap_reduction),
                         cols = tissue_cols_use, width = 11, height = 8)
  save_plot(p1$plot, file.path(fig_dir, "umap_tissue"), p1$width, p1$height)
  p2 <- build_umap_plot(obj, umap_reduction, CELLTYPE_L2_COL,
                         title = sprintf("Epithelial - Cell Type L2 (%s)", umap_reduction),
                         label = TRUE, width = 14, height = 10)
  save_plot(p2$plot, file.path(fig_dir, "umap_celltype_L2"), p2$width, p2$height)
  p2s <- build_umap_plot(obj, umap_reduction, CELLTYPE_L2_COL,
                          title = sprintf("Epithelial - L2 by Tissue (%s)", umap_reduction),
                          split_col = TISSUE_COL, label = TRUE,
                          width = max(12, 4.5 * length(tissues)), height = 8)
  save_plot(p2s$plot, file.path(fig_dir, "umap_L2_split_tissue"), p2s$width, p2s$height)
  p2_l3 <- build_umap_plot(obj, umap_reduction, LABEL_COL,
                             title = sprintf("Epithelial - Cell Type L3 (%s)", umap_reduction),
                             label = TRUE, width = 18, height = 12)
  save_plot(p2_l3$plot, file.path(fig_dir, "umap_celltype_L3"), p2_l3$width, p2_l3$height)
  cat(sprintf("[OK] UMAP saved using reduction: %s\n", umap_reduction))
} else { cat("[WARN] No UMAP reduction found\n") }

markers_present <- intersect(KNOWN_MARKERS, rownames(obj))
if (length(markers_present) >= 3) {
  Idents(obj) <- LABEL_COL
  p3 <- DotPlot(obj, features = markers_present) + RotatedAxis() +
    ggtitle("Epithelial - Known Markers (L3)") + theme(axis.text.x = element_text(size = 7))
  save_plot(p3, file.path(fig_dir, "dotplot_markers"),
            width = max(14, length(markers_present) * 0.35),
            height = max(7, length(unique(meta[[LABEL_COL]])) * 0.4))
  cat("[OK] Dotplot saved\n")
}

comp_df <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL))) %>%
  count(!!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(TISSUE_COL)) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

p4  <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(CELLTYPE_L2_COL))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Tissue", y = "Percentage (%)", title = "Epithelial - L2 Composition by Tissue") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4, file.path(fig_dir, "composition_tissue_L2"))

p4b <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = n, fill = !!sym(CELLTYPE_L2_COL))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Tissue", y = "Cell Count", title = "Epithelial - Absolute Count by Tissue") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4b, file.path(fig_dir, "count_tissue_L2"))

sample_comp <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL)), !is.na(!!sym(SAMPLE_COL))) %>%
  count(!!sym(SAMPLE_COL), !!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(SAMPLE_COL)) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
if (nrow(sample_comp) > 0) {
  facet_ncol <- min(3, length(l2_types)); facet_nrow <- ceiling(length(l2_types) / facet_ncol)
  p4c <- ggplot(sample_comp, aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(TISSUE_COL))) +
    geom_boxplot(outlier.size = 0.5) + geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(as.formula(paste("~", CELLTYPE_L2_COL)), scales = "free_y", ncol = facet_ncol) +
    labs(x = "Tissue", y = "Proportion per Sample (%)",
         title = "Epithelial - Sample-level Composition by Tissue") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  save_plot(p4c, file.path(fig_dir, "composition_sample_level"),
            width = max(10, facet_ncol * 4), height = max(6, facet_nrow * 3.5))
  cat("[OK] Sample-level composition saved\n")
}

Idents(obj) <- CELLTYPE_L2_COL
top_mk <- tryCatch(
  FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25,
                 max.cells.per.ident = 1000, test.use = "wilcox"),
  error = function(e) { cat(sprintf("[WARN] FindAllMarkers failed: %s\n", e$message)); NULL }
)
if (!is.null(top_mk) && nrow(top_mk) > 0) {
  if (!"gene" %in% colnames(top_mk)) top_mk$gene <- rownames(top_mk)
  top_mk <- filter_marker_table_by_gene(top_mk, gene_col = "gene", context = "FindAllMarkers L2", verbose = TRUE)
  fwrite(top_mk, file.path(rpt_dir, "all_markers_per_L2.csv"))
  if (nrow(top_mk) > 0) {
    top10  <- top_mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
    obj_ds <- stratified_downsample(obj, CELLTYPE_L2_COL, HEATMAP_CELLS_PER_TYPE)
    obj_ds <- ScaleData(obj_ds, features = unique(top10$gene), verbose = FALSE)
    p5 <- DoHeatmap(obj_ds, features = unique(top10$gene), size = 3) +
      ggtitle("Epithelial - Top Markers per L2 (downsampled, visualization only)")
    save_plot(p5, file.path(fig_dir, "heatmap_top_markers"), width = 14, height = 10)
    rm(obj_ds); gc()
    cat("[OK] Heatmap saved\n")
  } else {
    cat("[WARN] All L2 markers removed by IG filter; skipping heatmap\n")
  }
}
cat("\n")

# ==============================================================================
# 8. Pseudobulk DESeq2 + Multi-Database Enrichment + interpret_agent
# ==============================================================================

cat("=== Pseudobulk DESeq2 Tissue Comparison ===\n")
cat(sprintf("Thresholds: padj < %s, |log2FC| > %s\n", PADJ_THR, LFC_THR))
cat("Enrichment databases: GO_BP, GO_MF, GO_CC, KEGG, Hallmark, CellMarker, PanglaoDB, Epithelial_custom\n\n")

comparison_l2 <- run_pairwise_tissue_comparison(
  obj = obj,
  group_col = CELLTYPE_L2_COL,
  group_values = l2_types,
  level_name = "L2",
  out_dir = de_dir
)

pb_de_all            <- comparison_l2$pb_de
enrich_all           <- comparison_l2$enrich
agent_all            <- comparison_l2$agent
agent_structured_all <- comparison_l2$agent_structured

comparison_l3 <- run_pairwise_tissue_comparison(
  obj = obj,
  group_col = CELLTYPE_L3_COL,
  group_values = l3_types,
  level_name = "L3",
  out_dir = de_l3_dir
)

pb_de_l3_all            <- comparison_l3$pb_de
enrich_l3_all           <- comparison_l3$enrich
agent_l3_all            <- comparison_l3$agent
agent_structured_l3_all <- comparison_l3$agent_structured

# ==============================================================================
# 8.5 Pseudobulk ssGSEA by Tissue x L2/L3
# ==============================================================================
# Groups: tissue x L2 and tissue x L3 (e.g. "nose__Basal_Lineage", "nose__AT2")
# Output (same convention as OFA subcluster scripts):
#   reports/ssgsea_scores_{method}.rds
#   reports/ssgsea_z_{method}.rds
#   reports/ssgsea_top_pathways_{method}.csv
#   figures/ssgsea_heatmap_{method}.pdf/.png
# ==============================================================================

ssgsea_results_all <- list()
ssgsea_results_l3_all <- list()

if (RUN_SSGSEA) {
  ssgsea_results_all <- run_grouped_ssgsea(
    obj = obj,
    group_col = CELLTYPE_L2_COL,
    level_name = "L2",
    methods = SSGSEA_METHODS,
    n_top = SSGSEA_N_TOP,
    n_heatmap = SSGSEA_N_HEATMAP
  )

  ssgsea_results_l3_all <- run_grouped_ssgsea(
    obj = obj,
    group_col = CELLTYPE_L3_COL,
    level_name = "L3",
    methods = SSGSEA_METHODS,
    n_top = SSGSEA_N_TOP,
    n_heatmap = SSGSEA_N_HEATMAP
  )

  ssgsea_primary_method <- if ("hallmark" %in% names(ssgsea_results_all)) "hallmark" else
    if (length(ssgsea_results_all) > 0) names(ssgsea_results_all)[1] else NULL
  cat(sprintf("[INFO] Primary ssGSEA method: %s\n\n",
              if (is.null(ssgsea_primary_method)) "none" else ssgsea_primary_method))
}

# ==============================================================================
# 8.6 CHOIR Clustering (scANVI latent) + Per-cluster ssGSEA + OFA
# ==============================================================================
# CHOIR uses the scANVI latent space (not UMAP) for statistically grounded
# hierarchical clustering. Per-cluster OFA and ssGSEA use the same 8-DB stack.
#
# Output (mirrors OFA subcluster interpret scripts):
#   reports/choir/choir_clusters.csv
#   reports/choir/ofa_{cluster}/markers.csv, volcano.pdf/png
#   reports/choir/ofa_{cluster}/enrichment_{dir}/{db}.csv, dotplot.pdf
#   reports/choir/ssgsea_scores_{method}.rds, ssgsea_z_{method}.rds
#   reports/choir/ssgsea_top_pathways_{method}.csv
#   figures/choir_umap.pdf/.png, choir_umap_split_tissue.pdf/.png
#   figures/ssgsea_choir_heatmap_{method}.pdf/.png
# ==============================================================================

choir_results_all <- list()
choir_ofa_all     <- list()

if (RUN_CHOIR) {
  cat("\n=== CHOIR Clustering (scANVI latent) ===\n")

  choir_reduction <- pick_reduction(obj, CHOIR_REDUCTION_CANDIDATES)
  if (is.null(choir_reduction)) {
    cat("[INFO] No preferred reduction found; running PCA as CHOIR fallback\n")
    obj <- tryCatch(RunPCA(obj, verbose = FALSE),
                    error = function(e) { cat(sprintf("[ERROR] PCA: %s\n", e$message)); obj })
    choir_reduction <- if ("pca" %in% Reductions(obj)) "pca" else NULL
  }

  if (is.null(choir_reduction)) {
    cat("[WARN] No suitable reduction for CHOIR; skipping Section 8.6\n")
    RUN_CHOIR <- FALSE
  } else {
    cat(sprintf("[INFO] CHOIR reduction: %s\n", choir_reduction))
    choir_embedding <- tryCatch(
      as.matrix(Embeddings(obj, reduction = choir_reduction)),
      error = function(e) {
        cat(sprintf("[ERROR] Failed to extract CHOIR embeddings '%s': %s\n", choir_reduction, e$message))
        NULL
      }
    )
    if (is.null(choir_embedding)) {
      cat("[WARN] CHOIR embedding unavailable; skipping Section 8.6\n")
      RUN_CHOIR <- FALSE
    }

    choir_var_features <- VariableFeatures(obj)
    if (length(choir_var_features) == 0) {
      choir_var_features <- rownames(obj)
      cat(sprintf("[INFO] CHOIR var_features not set; using all %d genes\n", length(choir_var_features)))
    } else {
      cat(sprintf("[INFO] CHOIR var_features: %d features\n", length(choir_var_features)))
    }
  }

  if (RUN_CHOIR) {
    cat(sprintf("[INFO] Running CHOIR (alpha=%.3f, n_cores=%d)...\n", CHOIR_ALPHA, N_CORES))
    obj_choir <- tryCatch(
      CHOIR::CHOIR(obj, use_assay = "RNA", reduction = choir_embedding,
                   var_features = choir_var_features, n_cores = N_CORES,
                   alpha = CHOIR_ALPHA, random_seed = 42),
      error = function(e) { cat(sprintf("[ERROR] CHOIR failed: %s\n", e$message)); NULL }
    )
    if (!is.null(obj_choir)) obj <- obj_choir
    rm(obj_choir, choir_embedding, choir_var_features); gc()

    choir_col <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
    if (!choir_col %in% colnames(obj@meta.data)) {
      fallback  <- grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
      choir_col <- if (length(fallback) > 0) fallback[1] else NULL
      if (!is.null(choir_col)) cat(sprintf("[INFO] Using CHOIR column: %s\n", choir_col))
    }

    if (is.null(choir_col)) {
      cat("[WARN] CHOIR cluster column not found; skipping Section 8.6\n")
    } else {
      choir_clusters <- sort(unique(na.omit(obj@meta.data[[choir_col]])))
      cat(sprintf("[OK] CHOIR: %d clusters\n", length(choir_clusters)))

      choir_dir <- file.path(rpt_dir, "choir")
      dir.create(choir_dir, recursive = TRUE, showWarnings = FALSE)
      fwrite(data.frame(cell = colnames(obj), choir_cluster = obj@meta.data[[choir_col]]),
             file.path(choir_dir, "choir_clusters.csv"))
      cat("[OK] Cluster assignments saved: choir_clusters.csv\n")

      if (!is.null(umap_reduction)) {
        p_choir <- build_umap_plot(obj, umap_reduction, choir_col,
                                    title = sprintf("Epithelial - CHOIR Clusters (alpha=%.3f, %s)",
                                                    CHOIR_ALPHA, choir_reduction),
                                    label = TRUE, width = 14, height = 10)
        save_plot(p_choir$plot, file.path(fig_dir, "choir_umap"), p_choir$width, p_choir$height)
        p_choir_s <- build_umap_plot(obj, umap_reduction, choir_col,
                                      title = sprintf("Epithelial - CHOIR by Tissue (%s)", umap_reduction),
                                      split_col = TISSUE_COL, label = TRUE,
                                      width = max(12, 4.5 * length(tissues)), height = 8)
        save_plot(p_choir_s$plot, file.path(fig_dir, "choir_umap_split_tissue"),
                  p_choir_s$width, p_choir_s$height)
        cat("[OK] CHOIR UMAPs saved\n")
      }

      # ------------------------------------------------------------------
      # 8.6a ssGSEA per CHOIR cluster
      # ------------------------------------------------------------------
      if (RUN_SSGSEA && length(SSGSEA_CHOIR_METHODS) > 0) {
        cat("\n--- ssGSEA per CHOIR cluster ---\n")
        ssgsea_choir_all <- list()

        for (method in SSGSEA_CHOIR_METHODS) {
          cat(sprintf("  method: %s\n", method))
          selected_t2g_choir <- switch(
            method,
            hallmark = hallmark_t2g, go_bp = go_bp_t2g, go_mf = go_mf_t2g,
            go_cc = go_cc_t2g, kegg = kegg_t2g,
            custom = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT)) {
              tryCatch(read.gmt(SSGSEA_CUSTOM_GMT) %>% mutate(gene = toupper(gene)),
                       error = function(e) NULL)
            } else NULL,
            NULL
          )
          if (is.null(selected_t2g_choir) || nrow(selected_t2g_choir) == 0) {
            cat(sprintf("  [WARN] '%s': no gene sets\n", method)); next
          }
          gs_choir <- term2gene_to_list(selected_t2g_choir)
          gs_choir <- map_gene_sets_to_features(gs_choir, rownames(obj))
          gs_choir <- filter_gene_sets_for_analysis(gs_choir)
          gs_choir <- filter_gs_size(gs_choir)
          if (length(gs_choir) == 0) { cat(sprintf("  [WARN] '%s': no valid gene sets\n", method)); next }

          needed_genes_choir <- unique(unlist(gs_choir, use.names = FALSE))
          avg_choir <- tryCatch(
            AverageExpression(obj, assays = "RNA", slot = "data", group.by = choir_col,
                              features = needed_genes_choir, verbose = FALSE)[["RNA"]],
            error = function(e) { cat(sprintf("  [ERROR] AverageExpression: %s\n", e$message)); NULL }
          )
          if (is.null(avg_choir)) next

          bp_choir <- BiocParallel::SnowParam(workers = N_CORES, type = "SOCK", progressbar = FALSE)
          ss_choir <- tryCatch({
            if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
              GSVA::gsva(GSVA::ssgseaParam(exprData = as.matrix(avg_choir), geneSets = gs_choir,
                                            alpha = 0.25, normalize = TRUE, minSize = 10, maxSize = 500),
                         BPPARAM = bp_choir, verbose = FALSE)
            } else {
              GSVA::gsva(as.matrix(avg_choir), gs_choir, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
            }
          }, error = function(e) { cat(sprintf("  [ERROR] ssGSEA: %s\n", e$message)); NULL })
          if (is.null(ss_choir)) next

          sz_choir <- t(scale(t(ss_choir)))
          ssgsea_choir_all[[method]] <- list(scores = ss_choir, z_scores = sz_choir)
          saveRDS(ss_choir, file.path(choir_dir, sprintf("ssgsea_scores_%s.rds", method)))
          saveRDS(sz_choir, file.path(choir_dir, sprintf("ssgsea_z_%s.rds",     method)))

          fwrite(bind_rows(lapply(colnames(ss_choir), function(cid) {
            sc  <- ss_choir[, cid]; z <- sz_choir[, cid]
            ord <- order(z, sc, decreasing = TRUE, na.last = TRUE)
            ti  <- ord[seq_len(min(SSGSEA_CHOIR_N_TOP, length(ord)))]
            data.frame(method = method, cluster = cid, pathway = rownames(ss_choir)[ti],
                       score = sc[ti], z_score = z[ti], rank = seq_along(ti), stringsAsFactors = FALSE)
          })), file.path(choir_dir, sprintf("ssgsea_top_pathways_%s.csv", method)))

          maz     <- rowMeans(abs(sz_choir), na.rm = TRUE)
          top_p_c <- names(sort(maz, decreasing = TRUE))[seq_len(min(SSGSEA_CHOIR_N_HEATMAP, length(maz)))]
          hm_c    <- sz_choir[top_p_c, , drop = FALSE]
          hm_c[is.nan(hm_c) | is.infinite(hm_c)] <- 0

          if (nrow(hm_c) >= 3 && ncol(hm_c) >= 2) {
            tryCatch({
              hp <- file.path(fig_dir, sprintf("ssgsea_choir_heatmap_%s", method))
              ht_c <- pheatmap::pheatmap(
                hm_c, cluster_rows = TRUE, cluster_cols = TRUE,
                color  = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
                breaks = seq(-3, 3, length.out = 101),
                main   = sprintf("Epithelial CHOIR ssGSEA Z-score (%s) — top %d pathways", method, nrow(hm_c)),
                fontsize_row = 7, fontsize_col = 9, cellwidth = 22, cellheight = 10,
                filename = paste0(hp, ".pdf"),
                width  = max(10, ncol(hm_c) * 1.8),
                height = max(8,  nrow(hm_c) * 0.35 + 3)
              )
              grDevices::png(paste0(hp, ".png"),
                             width  = max(10, ncol(hm_c) * 1.8),
                             height = max(8,  nrow(hm_c) * 0.35 + 3),
                             units = "in", res = 300)
              grid::grid.newpage(); grid::grid.draw(ht_c$gtable); grDevices::dev.off()
              cat(sprintf("  [OK] Heatmap: ssgsea_choir_heatmap_%s\n", method))
            }, error = function(e) cat(sprintf("  [WARN] Heatmap %s: %s\n", method, e$message)))
          }
          cat(sprintf("  [OK] %s: %d pathways x %d clusters\n", method, nrow(ss_choir), ncol(ss_choir)))
        }
        choir_results_all[["ssgsea"]] <- ssgsea_choir_all
        cat("[OK] CHOIR ssGSEA complete\n")
      }

      # ------------------------------------------------------------------
      # 8.6b OFA: one-vs-rest FindMarkers + 8-DB enrichment per cluster
      # ------------------------------------------------------------------
      if (RUN_OFA) {
        cat("\n--- OFA (one-vs-rest) per CHOIR cluster ---\n")

        for (cl in choir_clusters) {
          cl_name     <- as.character(cl)
          focal_cells <- colnames(obj)[obj@meta.data[[choir_col]] == cl]
          rest_cells  <- colnames(obj)[obj@meta.data[[choir_col]] != cl]
          n_focal     <- length(focal_cells); n_rest <- length(rest_cells)
          cat(sprintf("  OFA: cluster %s vs rest\n", cl_name))

          if (n_focal < OFA_MIN_CELLS_FOCAL || n_rest < OFA_MIN_CELLS_REST) {
            cat(sprintf("    [SKIP] focal=%d rest=%d (min %d/%d)\n",
                        n_focal, n_rest, OFA_MIN_CELLS_FOCAL, OFA_MIN_CELLS_REST)); next
          }

          sub_ofa <- subset(obj, cells = c(focal_cells, rest_cells))
          sub_ofa@meta.data$ofa_group <- ifelse(sub_ofa@meta.data[[choir_col]] == cl, cl_name, "rest")
          Idents(sub_ofa) <- "ofa_group"

          if (n_focal > OFA_MAX_CELLS_PER_IDENT || n_rest > OFA_MAX_CELLS_PER_IDENT) {
            set.seed(42)
            kf <- sample(focal_cells, min(n_focal, OFA_MAX_CELLS_PER_IDENT))
            kr <- sample(rest_cells,  min(n_rest,  OFA_MAX_CELLS_PER_IDENT))
            sub_ofa <- subset(sub_ofa, cells = c(kf, kr)); Idents(sub_ofa) <- "ofa_group"
            cat(sprintf("    [INFO] Downsampled to %d focal + %d rest\n", length(kf), length(kr)))
          }

          de_ofa <- tryCatch(
            FindMarkers(sub_ofa, ident.1 = cl_name, ident.2 = "rest", test.use = "wilcox",
                        min.pct = 0.1, logfc.threshold = 0),
            error = function(e) { cat(sprintf("    [ERROR] FindMarkers: %s\n", e$message)); NULL }
          )
          rm(sub_ofa); gc()
          if (is.null(de_ofa) || nrow(de_ofa) == 0) { cat("    [SKIP] No markers returned\n"); next }

          de_ofa$gene <- rownames(de_ofa); de_ofa$padj <- p.adjust(de_ofa$p_val, method = "BH")
          de_ofa <- filter_marker_table_by_gene(
            de_ofa,
            gene_col = "gene",
            context = sprintf("CHOIR OFA cluster %s", cl_name),
            verbose = TRUE
          )
          if (nrow(de_ofa) == 0) { cat("    [SKIP] All markers removed by IG filter\n"); next }
          de_ofa <- de_ofa %>% filter(abs(avg_log2FC) >= OFA_LFC_THR) %>% arrange(padj)
          n_sig_ofa <- sum(de_ofa$padj < OFA_PADJ_THR, na.rm = TRUE)
          cat(sprintf("    %d sig markers (padj<%.2f |log2FC|>%.2f)\n",
                      n_sig_ofa, OFA_PADJ_THR, OFA_LFC_THR))

          ofa_dir <- file.path(choir_dir, sprintf("ofa_%s", safe_name(cl_name)))
          dir.create(ofa_dir, recursive = TRUE, showWarnings = FALSE)
          fwrite(de_ofa, file.path(ofa_dir, "markers.csv"))

          de_v <- de_ofa %>%
            mutate(sig_flag = padj < OFA_PADJ_THR & abs(avg_log2FC) >= OFA_LFC_THR,
                   label    = ifelse(sig_flag & dplyr::row_number() <= 20, gene, ""))
          pv_ofa <- ggplot(de_v, aes(x = avg_log2FC, y = -log10(padj), color = sig_flag)) +
            geom_point(alpha = 0.5, size = 0.8) +
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), guide = "none") +
            geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
            geom_vline(xintercept = c(-OFA_LFC_THR, OFA_LFC_THR), linetype = "dashed", color = "blue") +
            labs(title = sprintf("OFA: Epithelial CHOIR cluster %s vs rest", cl_name),
                 x = "avg_log2FC", y = "-log10(padj)") + theme_minimal()
          save_plot(pv_ofa, file.path(ofa_dir, "volcano"), width = 8, height = 6)

          ofa_enrich <- list(); all_tested_ofa <- de_ofa$gene
          for (dir_name in c("up", "down")) {
            ofa_genes <- if (dir_name == "up") {
              de_ofa %>% filter(padj < OFA_PADJ_THR, avg_log2FC > 0) %>%
                arrange(desc(avg_log2FC)) %>% head(OFA_TOP_N) %>% pull(gene)
            } else {
              de_ofa %>% filter(padj < OFA_PADJ_THR, avg_log2FC < 0) %>%
                arrange(avg_log2FC) %>% head(OFA_TOP_N) %>% pull(gene)
            }
            if (length(ofa_genes) < 5) next
            cat(sprintf("    Enrichment [%s]: %d genes\n", dir_name, length(ofa_genes)))

            enr_ofa <- list(
              GO_BP             = run_gmt_enrichment(ofa_genes, go_bp_t2g,            "GO_BP",            all_tested_ofa),
              GO_MF             = run_gmt_enrichment(ofa_genes, go_mf_t2g,            "GO_MF",            all_tested_ofa),
              GO_CC             = run_gmt_enrichment(ofa_genes, go_cc_t2g,            "GO_CC",            all_tested_ofa),
              KEGG              = run_gmt_enrichment(ofa_genes, kegg_t2g,             "KEGG",             all_tested_ofa),
              Hallmark          = run_gmt_enrichment(ofa_genes, hallmark_t2g,         "Hallmark",         all_tested_ofa),
              CellMarker        = run_gmt_enrichment(ofa_genes, cellmarker_t2g,       "CellMarker",       all_tested_ofa),
              PanglaoDB         = run_gmt_enrichment(ofa_genes, panglaodb_t2g,        "PanglaoDB",        all_tested_ofa),
              Epithelial_custom = run_gmt_enrichment(ofa_genes, epithelial_custom_t2g,"Epithelial_custom",all_tested_ofa)
            )
            enr_ofa <- enr_ofa[!sapply(enr_ofa, is.null)]
            ofa_enrich[[dir_name]] <- enr_ofa

            enr_out_ofa <- file.path(ofa_dir, paste0("enrichment_", dir_name))
            dir.create(enr_out_ofa, showWarnings = FALSE)
            for (db in names(enr_ofa)) {
              er <- enr_ofa[[db]]
              if (nrow(as.data.frame(er)) > 0) {
                fwrite(as.data.frame(er), file.path(enr_out_ofa, paste0(db, ".csv")))
                tryCatch({
                  pdf(file.path(enr_out_ofa, paste0(db, "_dotplot.pdf")), width = 10, height = 8)
                  print(dotplot(er, showCategory = 15,
                                title = sprintf("%s %s (CHOIR %s vs rest)", db, dir_name, cl_name)))
                  dev.off()
                }, error = function(e) NULL)
              }
            }
            saveRDS(enr_ofa, file.path(enr_out_ofa, "all_enrichment.rds"))
          }

          choir_ofa_all[[cl_name]] <- list(de_table = de_ofa, enrich = ofa_enrich,
                                            n_focal = n_focal, n_rest = n_rest)
        }

        choir_results_all[["ofa"]] <- choir_ofa_all
        cat(sprintf("[OK] OFA complete: %d clusters processed\n", length(choir_ofa_all)))
      }

      saveRDS(choir_results_all, file.path(choir_dir, "choir_results_all.rds"))
      saveRDS(choir_ofa_all,     file.path(choir_dir, "choir_ofa_all.rds"))
      cat(sprintf("[OK] CHOIR results saved: %s\n", choir_dir))
    }
  }
}

# ==============================================================================
# 9. Exploratory Wilcoxon (marker discovery only)
# ==============================================================================

cat("\n=== Exploratory Cell-level Wilcoxon ===\n")
cat("NOTE: P-values are inflated (pseudoreplication). For marker discovery ONLY.\n\n")

wilcox_all <- list()
for (ct_l2 in l2_types) {
  wx <- run_wilcox_exploratory(obj, ct_l2, group_col = CELLTYPE_L2_COL)
  if (!is.null(wx) && length(wx) > 0) {
    wilcox_all[[ct_l2]] <- wx
    for (comp_name in names(wx)) {
      wx_out <- file.path(wx_dir, safe_name(ct_l2))
      dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(wx[[comp_name]], file.path(wx_out, paste0(safe_name(comp_name), "_wilcox.csv")))
    }
    cat(sprintf("  [OK] L2 %s: %d comparisons\n", ct_l2, length(wx)))
  }
}

wilcox_l3_all <- list()
for (ct_l3 in l3_types) {
  wx <- run_wilcox_exploratory(obj, ct_l3, group_col = CELLTYPE_L3_COL)
  if (!is.null(wx) && length(wx) > 0) {
    wilcox_l3_all[[ct_l3]] <- wx
    for (comp_name in names(wx)) {
      wx_out <- file.path(wx_l3_dir, safe_name(ct_l3))
      dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(wx[[comp_name]], file.path(wx_out, paste0(safe_name(comp_name), "_wilcox.csv")))
    }
    cat(sprintf("  [OK] L3 %s: %d comparisons\n", ct_l3, length(wx)))
  }
}

# ==============================================================================
# 10. Save Summary RDS
# ==============================================================================

cat("\n=== Saving Summary ===\n")
saveRDS(pb_de_all,            file.path(rpt_dir, "pseudobulk_de_all.rds"))
saveRDS(pb_de_l3_all,         file.path(rpt_dir, "pseudobulk_de_L3_all.rds"))
saveRDS(wilcox_all,           file.path(rpt_dir, "wilcox_exploratory_all.rds"))
saveRDS(wilcox_l3_all,        file.path(rpt_dir, "wilcox_exploratory_L3_all.rds"))
saveRDS(enrich_all,           file.path(rpt_dir, "enrichment_all.rds"))
saveRDS(enrich_l3_all,        file.path(rpt_dir, "enrichment_L3_all.rds"))
saveRDS(agent_all,            file.path(rpt_dir, "interpret_agent_all.rds"))
saveRDS(agent_l3_all,         file.path(rpt_dir, "interpret_agent_L3_all.rds"))
saveRDS(agent_structured_all, file.path(rpt_dir, "interpret_agent_structured_all.rds"))
saveRDS(agent_structured_l3_all, file.path(rpt_dir, "interpret_agent_structured_L3_all.rds"))
saveRDS(L3_TO_L2_REMAP,       file.path(rpt_dir, "L3_to_L2_remap.rds"))

agent_structured_df_l2 <- flatten_interpretation_records(agent_structured_all)
agent_structured_df_l3 <- flatten_interpretation_records(agent_structured_l3_all)
agent_structured_df <- bind_rows(agent_structured_df_l2, agent_structured_df_l3)
if (nrow(agent_structured_df) > 0)
  fwrite(agent_structured_df, file.path(rpt_dir, "interpret_agent_structured.tsv"), sep = "\t")
if (nrow(agent_structured_df_l3) > 0)
  fwrite(agent_structured_df_l3, file.path(rpt_dir, "interpret_agent_structured_L3.tsv"), sep = "\t")

write_interpretation_markdown(
  agent_structured_df, file.path(OUTPUT_DIR, "LLM_INTERPRETATION.md"),
  "# LLM Interpretation Summary (Epithelial Tissue Comparison)", include_raw = FALSE)
write_interpretation_markdown(
  agent_structured_df, file.path(rpt_dir, "LLM_INTERPRETATION_FOR_LLM.md"),
  "# LLM Interpretation Structured Input", include_raw = TRUE)

# ==============================================================================
# 11. Generate REPORT.md
# ==============================================================================

cat("\n=== Generating REPORT.md ===\n")
md <- character(); add <- function(...) md <<- c(md, paste0(...))

add("# Epithelial Tissue Comparison Report")
add(""); add("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M")); add("")
add(sprintf("**Pipeline:** Epithelial Tissue Comparison %s (pseudobulk DESeq2, ssGSEA, CHOIR, multi-database enrichment)", PIPELINE_VERSION))
add("")
add("**Note:** Cross-site anatomical comparison of respiratory epithelial cells.")
add("**Interpretation caution:** Unless disease covariates are explicitly modeled, do not over-interpret differences as disease effects.")
add(""); add("---"); add("")

add("## 1. Data Overview"); add("")
add(sprintf("- **Input:** `%s`", basename(H5AD_PATH)))
add(sprintf("- **Total cells:** %s", format(ncol(obj), big.mark = ",")))
add(sprintf("- **Tissues:** %s", paste(tissues, collapse = ", ")))
add(sprintf("- **L2 lineages (%d):** %s", length(l2_types), paste(l2_types, collapse = ", ")))
add(sprintf("- **L3 states (%d):** %s", length(l3_types), paste(l3_types, collapse = ", ")))
add(sprintf("- **L3 source column:** `%s` (%d unique)", L3_SOURCE_COL, length(unique(meta[[LABEL_COL]]))))
if (FILTER_IG_GENES) add(sprintf("- **Analysis gene filter:** IG-related genes removed using regex `%s`", IG_GENE_REGEX))
add(""); add("### L3 -> L2 Remapping"); add(""); add("| L3 | L2 |"); add("|---|---|")
for (i in seq_along(L3_TO_L2_REMAP)) add(sprintf("| %s | %s |", names(L3_TO_L2_REMAP)[i], L3_TO_L2_REMAP[i]))
add("")

add("## 2. Visualization"); add("")
add("### 2.1 UMAP")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_tissue.png"),      "UMAP tissue")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_celltype_L2.png"), "UMAP L2")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_celltype_L3.png"), "UMAP L3")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_L2_split_tissue.png"), "UMAP split")
add("### 2.2 Marker Dotplot")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "dotplot_markers.png"), "Dotplot")
add("### 2.3 Cell Composition")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "composition_tissue_L2.png"), "Composition")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "count_tissue_L2.png"),       "Counts")
add("### 2.4 Sample-level Composition")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "composition_sample_level.png"), "Sample composition")
add("### 2.5 Top Marker Heatmap")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "heatmap_top_markers.png"), "Heatmap")

add("## 3. Pseudobulk DESeq2 (Primary Inference)"); add("")
add(sprintf("Thresholds: padj < %s, |log2FC| > %s.", PADJ_THR, LFC_THR)); add("")
add("### 3.1 L2 comparisons"); add("")
add("| L2 Lineage | Comparison | Samples (ref/case) | Design | Up | Down | Total |")
add("|---|---|---|---|---|---|---|")
for (ct_l2 in names(pb_de_all)) for (comp_name in names(pb_de_all[[ct_l2]])) {
  r <- pb_de_all[[ct_l2]][[comp_name]]; if (is.null(r)) next
  add(sprintf("| %s | %s | %d / %d | `%s` | %d | %d | %d |",
              ct_l2, comp_name, r$n_samples_1, r$n_samples_2, r$design,
              r$n_up, r$n_down, r$n_up + r$n_down))
}
add("")
add("### 3.2 L3 comparisons"); add("")
add("| L3 State | Comparison | Samples (ref/case) | Design | Up | Down | Total |")
add("|---|---|---|---|---|---|---|")
for (ct_l3 in names(pb_de_l3_all)) for (comp_name in names(pb_de_l3_all[[ct_l3]])) {
  r <- pb_de_l3_all[[ct_l3]][[comp_name]]; if (is.null(r)) next
  add(sprintf("| %s | %s | %d / %d | `%s` | %d | %d | %d |",
              ct_l3, comp_name, r$n_samples_1, r$n_samples_2, r$design,
              r$n_up, r$n_down, r$n_up + r$n_down))
}
add("")

add("## 4. Multi-Database Enrichment"); add("")
add("Databases: GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, custom epithelial markers"); add("")
add("### 4.1 L2 tissue-pair enrichment"); add("")
for (ct_l2 in names(enrich_all)) for (comp_name in names(enrich_all[[ct_l2]])) for (dir_name in names(enrich_all[[ct_l2]][[comp_name]])) {
  enr_l <- enrich_all[[ct_l2]][[comp_name]][[dir_name]]
  if (is.null(enr_l) || length(enr_l) == 0) next
  add(sprintf("### %s | %s | %s", ct_l2, comp_name, dir_name)); add("")
  for (db in names(enr_l)) {
    er <- enr_l[[db]]; if (is.null(er) || nrow(as.data.frame(er)) == 0) next
    top5 <- head(as.data.frame(er), 5)
    add(sprintf("**%s (top 5):**", db)); add(""); add("| Term | p.adjust | Count |"); add("|---|---|---|")
    for (j in seq_len(nrow(top5))) add(sprintf("| %s | %.2e | %s |", top5$Description[j], top5$p.adjust[j], top5$Count[j]))
    add("")
  }
}
add("### 4.2 L3 tissue-pair enrichment"); add("")
for (ct_l3 in names(enrich_l3_all)) for (comp_name in names(enrich_l3_all[[ct_l3]])) for (dir_name in names(enrich_l3_all[[ct_l3]][[comp_name]])) {
  enr_l <- enrich_l3_all[[ct_l3]][[comp_name]][[dir_name]]
  if (is.null(enr_l) || length(enr_l) == 0) next
  add(sprintf("### %s | %s | %s", ct_l3, comp_name, dir_name)); add("")
  for (db in names(enr_l)) {
    er <- enr_l[[db]]; if (is.null(er) || nrow(as.data.frame(er)) == 0) next
    top5 <- head(as.data.frame(er), 5)
    add(sprintf("**%s (top 5):**", db)); add(""); add("| Term | p.adjust | Count |"); add("|---|---|---|")
    for (j in seq_len(nrow(top5))) add(sprintf("| %s | %.2e | %s |", top5$Description[j], top5$p.adjust[j], top5$Count[j]))
    add("")
  }
}

# Section 4.5 ssGSEA
add("## 4.5 Pseudobulk ssGSEA"); add("")
if (!RUN_SSGSEA || (length(ssgsea_results_all) == 0 && length(ssgsea_results_l3_all) == 0)) {
  add("ssGSEA not run or no methods succeeded."); add("")
} else {
  add(sprintf("Methods: %s | Groups: tissue x L2/L3 combinations",
              paste(unique(c(names(ssgsea_results_all), names(ssgsea_results_l3_all))), collapse = ", "))); add("")
  add("### 4.5.1 Tissue x L2"); add("")
  add("| Method | Pathways | Groups | Output |"); add("|---|---|---|---|")
  for (mn in names(ssgsea_results_all))
    add(sprintf("| %s | %d | %d | `reports/ssgsea_scores_%s.rds` |",
                mn, nrow(ssgsea_results_all[[mn]]$scores), ncol(ssgsea_results_all[[mn]]$scores), mn))
  add("")
  if (length(ssgsea_results_all) == 0) add("No L2 ssGSEA results generated.")
  for (mn in names(ssgsea_results_all)) {
    add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", sprintf("ssgsea_heatmap_%s.png", mn)),
                           sprintf("ssGSEA heatmap %s", mn))
  }

  add("### 4.5.2 Tissue x L3"); add("")
  add("| Method | Pathways | Groups | Output |"); add("|---|---|---|---|")
  for (mn in names(ssgsea_results_l3_all))
    add(sprintf("| %s | %d | %d | `reports/ssgsea_l3_scores_%s.rds` |",
                mn, nrow(ssgsea_results_l3_all[[mn]]$scores), ncol(ssgsea_results_l3_all[[mn]]$scores), mn))
  add("")
  if (length(ssgsea_results_l3_all) == 0) add("No L3 ssGSEA results generated.")
  for (mn in names(ssgsea_results_l3_all)) {
    add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", sprintf("ssgsea_l3_heatmap_%s.png", mn)),
                           sprintf("L3 ssGSEA heatmap %s", mn))
  }
}

# Section 4.6 CHOIR
add("## 4.6 CHOIR Clustering + Per-cluster ssGSEA + OFA"); add("")
if (!RUN_CHOIR || length(choir_ofa_all) == 0) {
  add("CHOIR not run or no clusters processed."); add("")
} else {
  choir_col_report <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
  if (!choir_col_report %in% colnames(obj@meta.data)) {
    fb <- grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
    choir_col_report <- if (length(fb) > 0) fb[1] else NULL
  }
  n_cl_report <- if (!is.null(choir_col_report) && choir_col_report %in% colnames(obj@meta.data))
    length(unique(na.omit(obj@meta.data[[choir_col_report]]))) else NA
  add(sprintf("- **CHOIR alpha:** %.3f | **Clusters:** %s | **Reduction:** `%s`",
              CHOIR_ALPHA, ifelse(is.na(n_cl_report), "unknown", n_cl_report),
              if (exists("choir_reduction")) choir_reduction else "N/A"))
  add("")
  add("### 4.6.1 CHOIR UMAP")
  add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "choir_umap.png"),              "CHOIR UMAP")
  add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "choir_umap_split_tissue.png"), "CHOIR UMAP by Tissue")
  if (length(choir_results_all[["ssgsea"]]) > 0) {
    add("### 4.6.2 ssGSEA per CHOIR cluster"); add("")
    for (mn in names(choir_results_all[["ssgsea"]]))
      add_md_image_if_exists(add, OUTPUT_DIR,
                             file.path("figures", sprintf("ssgsea_choir_heatmap_%s.png", mn)),
                             sprintf("CHOIR ssGSEA %s", mn))
  }
  add("### 4.6.3 OFA One-vs-Rest Markers"); add("")
  add("| Cluster | n_focal | n_rest | Sig markers |"); add("|---|---|---|---|")
  for (cl_name in names(choir_ofa_all)) {
    r     <- choir_ofa_all[[cl_name]]
    n_sig <- if (!is.null(r$de_table)) sum(r$de_table$padj < OFA_PADJ_THR, na.rm = TRUE) else 0
    add(sprintf("| %s | %d | %d | %d |", cl_name, r$n_focal, r$n_rest, n_sig))
  }
  add(""); add("`Full OFA results: reports/choir/ofa_*/markers.csv, enrichment_{up|down}/{db}.csv`"); add("")
}

add("## 5. LLM Interpretation (interpret_agent)"); add("")
if (!ENABLE_LLM) {
  add("**SKIPPED:** DEEPSEEK_API_KEY not set."); add("")
} else if (nrow(agent_structured_df) == 0) {
  add("No interpret_agent results generated."); add("")
} else {
  add("Structured outputs:")
  add(sprintf("- `%s`", file.path("reports", "interpret_agent_structured.tsv")))
  add(sprintf("- `%s`", file.path("reports", "interpret_agent_structured_L3.tsv")))
  add(sprintf("- `%s`", "LLM_INTERPRETATION.md"))
  add("- LLM output language: Chinese (Simplified).")
  add("- LLM emphasis: top-|log2FC| genes linked to site-relevant epithelial pathways.")
  add("- LLM evidence integration: all eligible databases merged into one consensus summary.")
  add("- Structured rows merge both L2 and L3 results via `celltype_level` + `celltype_label`.")
  add("")
  for (i in seq_len(nrow(agent_structured_df))) {
    rec <- agent_structured_df[i, , drop = FALSE]
    add(sprintf("### %s | %s | %s | %s",
                rec$celltype_level, rec$celltype_label, rec$comparison, rec$direction)); add("")
    add(sprintf("- **Status:** %s", rec$status))
    add(sprintf("- **Cell Type Level:** %s", rec$celltype_level))
    add(sprintf("- **Source DB(s):** %s", ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
    if (nzchar(rec$overview))       add(sprintf("- **Overview:** %s", rec$overview))
    if (nzchar(rec$key_mechanisms)) add(sprintf("- **Key Mechanisms:** %s", rec$key_mechanisms))
    if (nzchar(rec$hypothesis))     add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
    add("")
  }
}

add("## 6. Methods"); add("")
add("- **DE:** Pseudobulk DESeq2 (sample-level aggregation; Squair et al. 2021 Nat Commun)")
add("- **Exploratory:** Cell-level Wilcoxon rank-sum (marker discovery only, NOT inference)")
add("- **Enrichment:** clusterProfiler::enricher() + MSigDB GMT + CellMarker + PanglaoDB + custom epithelial markers (symbol-based)")
add("- **Enrichment universe:** DESeq2-tested genes intersected with TERM2GENE")
add("- **Enrichment databases (8):** GO BP, GO MF, GO CC, KEGG, Hallmark, CellMarker, PanglaoDB, epithelial custom")
add("- **Enrichment size rules:** adaptive per database (ENRICHMENT_GS_SIZE_RULES)")
if (FILTER_IG_GENES) add(sprintf("- **Gene filtering:** IG-related genes removed from DE / Wilcoxon / enrichment / ssGSEA / OFA using regex `%s`", IG_GENE_REGEX))
add(sprintf("- **ssGSEA:** GSVA::ssgseaParam / gsva() on AverageExpression per tissue x L2/L3; methods: %s",
            paste(SSGSEA_METHODS, collapse = ", ")))
add("- **ssGSEA grouping:** tissue x L2 and tissue x L3 (e.g. 'nose__Basal_Lineage', 'nose__AT2'); Z-score row-wise across groups")
if (ENABLE_LLM) {
  add("- **LLM:** interpret_agent with DeepSeek deepseek-reasoner + multi-database integration")
  add("- **LLM post-processing:** secondary deepseek-chat pass standardizes outputs to fixed JSON schema")
  add("- **LLM output language:** Chinese (Simplified); gene symbols retained in English")
} else {
  add("- **LLM:** SKIPPED (no API key)")
}
if (!is.null(DESIGN_COVARIATE)) {
  add(sprintf("- **Requested covariate:** `%s`", DESIGN_COVARIATE))
} else {
  add("- **DESeq2 design:** ~ tissue (no covariate)")
}
add("- **Analysis levels:** both L2 (major lineages) and L3 (fine epithelial states)")
add("- **L2 remapping:** 21 L3 states merged into 6 major lineages")
add(sprintf("- **Context:** Respiratory epithelial site comparison across %s", paste(tissues, collapse = ", ")))
add(""); add("### Output Objects"); add("")
add("- `epithelial_tissue_comparison_final.rds` -- Seurat object (cell_type_L2 + cell_type_L3)")
add("- `epithelial_tissue_comparison_final.h5ad` -- AnnData object (cell_type_L2 + cell_type_L3)")
add(""); add("---"); add(sprintf("*Generated by %s*", SCRIPT_BASENAME))

writeLines(md, file.path(OUTPUT_DIR, "REPORT.md"))
cat("[OK] REPORT.md written\n")

# ==============================================================================
# 12. Save Final Object (RDS + h5ad)
# ==============================================================================

cat("\n=== Saving Final Object (RDS + h5ad) ===\n")

stopifnot("cell_type_L2" %in% colnames(obj@meta.data))
stopifnot("cell_type_L3" %in% colnames(obj@meta.data))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L2"]])))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L3"]])))

cat(sprintf("[OK] cell_type_L2: %d unique (%s)\n",
            length(unique(obj@meta.data[["cell_type_L2"]])),
            paste(sort(unique(obj@meta.data[["cell_type_L2"]])), collapse = ", ")))
cat(sprintf("[OK] cell_type_L3: %d unique (%s)\n",
            length(unique(obj@meta.data[["cell_type_L3"]])),
            paste(sort(unique(obj@meta.data[["cell_type_L3"]])), collapse = ", ")))

rds_path <- file.path(OUTPUT_DIR, "epithelial_tissue_comparison_final.rds")
saveRDS(obj, rds_path)
cat(sprintf("[OK] RDS saved: %s (%.1f MB)\n", basename(rds_path), file.size(rds_path) / 1e6))

h5ad_out_path <- file.path(OUTPUT_DIR, "epithelial_tissue_comparison_final.h5ad")
tryCatch({
  anndata      <- reticulate::import("anndata", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  np           <- reticulate::import("numpy", convert = FALSE)

  counts_mat      <- GetAssayData(obj, layer = "counts")
  norm_layer_name <- if ("data" %in% Layers(obj[["RNA"]])) "data" else "counts"
  expr_mat        <- GetAssayData(obj, layer = norm_layer_name)

  counts_scipy <- matrix_to_scipy_csr(counts_mat, scipy_sparse, np)
  expr_scipy   <- if (identical(norm_layer_name, "counts")) counts_scipy else
    matrix_to_scipy_csr(expr_mat, scipy_sparse, np)

  meta_export <- sanitize_obs_for_h5ad(obj@meta.data)
  obs_df      <- reticulate::r_to_py(meta_export)
  var_df      <- data.frame(gene_symbol = rownames(obj), row.names = rownames(obj), stringsAsFactors = FALSE)
  var_py      <- reticulate::r_to_py(var_df)

  adata <- anndata$AnnData(X = expr_scipy, obs = obs_df, var = var_py)
  adata$layers$`__setitem__`("counts", counts_scipy)
  adata$uns$`__setitem__`("X_layer", norm_layer_name)
  for (red_name in Reductions(obj))
    adata$obsm$`__setitem__`(paste0("X_", red_name), np$array(Embeddings(obj, reduction = red_name), dtype = np$float32))

  adata$write_h5ad(h5ad_out_path, compression = "gzip")
  cat(sprintf("[OK] h5ad saved: %s (%.1f MB)\n", basename(h5ad_out_path), file.size(h5ad_out_path) / 1e6))

  adata_check <- anndata$read_h5ad(h5ad_out_path)
  builtins     <- reticulate::import("builtins", convert = FALSE)
  obs_cols     <- reticulate::py_to_r(builtins$list(adata_check$obs$columns))
  layer_keys   <- reticulate::py_to_r(builtins$list(adata_check$layers$keys()))
  stopifnot("cell_type_L2" %in% obs_cols)
  stopifnot("cell_type_L3" %in% obs_cols)
  stopifnot("counts" %in% layer_keys)
  cat(sprintf("[OK] h5ad verified: %d cells x %d genes, cell_type_L2 + cell_type_L3 present, layers['counts'] present\n",
              reticulate::py_to_r(adata_check$n_obs), reticulate::py_to_r(adata_check$n_vars)))
  rm(adata, adata_check, builtins, obs_cols, layer_keys,
     counts_mat, expr_mat, counts_scipy, expr_scipy, meta_export, var_df); gc()
}, error = function(e) {
  cat(sprintf("[ERROR] h5ad export failed: %s\n", e$message))
  cat("[INFO] RDS saved successfully; convert manually if h5ad needed.\n")
})

# ==============================================================================
# 13. Final Summary
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("EPITHELIAL TISSUE COMPARISON COMPLETE (%s)\n", PIPELINE_VERSION))
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Directory structure:\n")
cat("  REPORT.md                                <- structured report (png embeds)\n")
cat("  LLM_INTERPRETATION.md                    <- LLM interpretation summary (Chinese)\n")
cat("  epithelial_tissue_comparison_final.rds   <- Seurat object (cell_type_L2 + cell_type_L3)\n")
cat("  epithelial_tissue_comparison_final.h5ad  <- AnnData object (cell_type_L2 + cell_type_L3)\n")
cat("  figures/                                 <- UMAP, dotplot, heatmap, composition (pdf+png)\n")
cat("  figures/ssgsea_heatmap_{method}.*        <- ssGSEA Z-score heatmap per method (pdf+png)\n")
cat("  figures/ssgsea_l3_heatmap_{method}.*     <- L3 ssGSEA Z-score heatmap per method (pdf+png)\n")
cat("  figures/choir_umap.*                     <- CHOIR cluster UMAP (pdf+png)\n")
cat("  figures/ssgsea_choir_heatmap_{method}.*  <- CHOIR ssGSEA heatmap (pdf+png)\n")
cat("  pseudobulk_de/                           <- DESeq2 per L2 x tissue pair + volcano + 8-DB enrichment\n")
cat("  pseudobulk_de_L3/                        <- DESeq2 per L3 x tissue pair + volcano + 8-DB enrichment\n")
cat("  wilcox_exploratory/                      <- cell-level wilcox (marker discovery only)\n")
cat("  wilcox_exploratory_L3/                   <- L3 cell-level wilcox (marker discovery only)\n")
cat("  reports/                                 <- RDS summary objects + structured LLM files\n")
cat("  reports/ssgsea_scores_{method}.rds       <- raw ssGSEA scores (pathways x groups)\n")
cat("  reports/ssgsea_z_{method}.rds            <- Z-score matrix\n")
cat("  reports/ssgsea_top_pathways_{method}.csv <- top pathways per tissue x L2 group\n")
cat("  reports/ssgsea_l3_scores_{method}.rds    <- L3 raw ssGSEA scores (pathways x groups)\n")
cat("  reports/ssgsea_l3_z_{method}.rds         <- L3 Z-score matrix\n")
cat("  reports/ssgsea_l3_top_pathways_{method}.csv <- top pathways per tissue x L3 group\n")
cat("  reports/interpret_agent_structured.tsv   <- merged L2+L3 structured interpret table\n")
cat("  reports/interpret_agent_structured_L3.tsv <- L3-only structured interpret table\n")
cat("  reports/choir/                           <- CHOIR clusters + per-cluster OFA + ssGSEA\n")
cat("  reports/choir/choir_clusters.csv         <- cell -> cluster assignment\n")
cat("  reports/choir/ssgsea_scores_{method}.rds / ssgsea_z_{method}.rds\n")
cat("  reports/choir/ssgsea_top_pathways_{method}.csv\n")
cat("  reports/choir/ofa_{cluster}/             <- per-cluster markers + 8-DB enrichment\n")
cat("\n")

cat(sprintf("%s Changes (vs v1.1.1):\n", PIPELINE_VERSION))
cat("  [NEW-12] Create a standalone v1.2 release script with independent output root\n")
cat("  [FIX-1] Unify version / script-name / output-path strings via constants\n")
cat("  [FIX-2] Correct REPORT.md footer and summary labels to avoid stale version text\n")
cat("\n")
cat("v1.1.1 Changes (vs v1.1.0):\n")
cat("  [NEW-10] Add full L3 tissue-comparison workflow alongside L2\n")
cat("  [NEW-11] Remove IG-related genes from DE / Wilcoxon / enrichment / ssGSEA / OFA\n")
cat("\n")
cat("v1.1.0 Changes (vs v1.0.0):\n")
cat("  [NEW-1] ENRICHMENT_GS_SIZE_RULES: adaptive min/maxGSSize per database\n")
cat("  [NEW-2] get_enrichment_size_rule(): helper for per-DB size lookup\n")
cat("  [NEW-3] run_gmt_enrichment(): upgraded to adaptive size rules\n")
cat("  [NEW-4] ssGSEA helpers: term2gene_to_list(), map_gene_sets_to_features(), filter_gs_size()\n")
cat("  [NEW-5] Section 8.5: pseudobulk ssGSEA by tissue x L2\n")
cat("  [NEW-6] Section 8.6: CHOIR clustering + per-cluster ssGSEA (8.6a) + OFA 8-DB (8.6b)\n")
cat("  [NEW-7] Library loading: GSVA, BiocParallel, CHOIR\n")
cat("  [NEW-8] Config: RUN_SSGSEA/CHOIR/OFA + all related parameters\n")
cat("  [NEW-9] REPORT.md sections 4.5 (ssGSEA) + 4.6 (CHOIR)\n")
cat("\n")

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
cat("[OK] Done\n")
