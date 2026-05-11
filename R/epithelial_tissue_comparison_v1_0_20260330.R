#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Tissue Comparison Pipeline v1.0.0 - PRODUCTION
# ==============================================================================
#
# Purpose:
#   Respiratory epithelial tissue comparison across anatomical sites:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. L3 -> L2 remapping (fine epithelial states -> major lineages)
#     3. Visualization (UMAP, marker dotplot/heatmap, composition)
#     4. Pseudobulk DESeq2 per L2 lineage x tissue pair
#     5. Exploratory cell-level wilcox (marker discovery only)
#     6. Multi-database enrichment (GO BP/MF/CC, KEGG, Hallmark,
#        CellMarker, PanglaoDB, custom epithelial markers)
#     7. interpret_agent (DeepSeek) with epithelial tissue-pair-specific context
#     8. Structured REPORT.md (png embeds)
#
# Notes:
#   - This is intended for cross-site anatomical comparison of respiratory
#     epithelial cells. Unless disease covariates are explicitly modeled,
#     do NOT over-interpret site-associated differences as disease effects.
#   - Uses the finalized epithelial scANVI SELF reference with curated L3 labels.
#
# Author: r2end + GitHub Copilot
# Date:   2026-03-30
# ==============================================================================

# ==============================================================================
# 0. Thread Control
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ==============================================================================
# 1. Configuration
# ==============================================================================

# ----- Input / Output -----
H5AD_PATH <- "/home/h2048/data/core20260322/epithelial_scanvi_v2_7_HOTFIX_SELF_final.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0330/epithelial_tissue_comparison_v1_0"

# ----- L3 -> L2 Remapping -----
L3_SOURCE_COL <- "cell_type_L3"

L3_TO_L2_REMAP <- c(
  "AT1_Canonical" = "Alveolar",
  "AT1_MatrixRemodeling" = "Alveolar",
  "AT2" = "Alveolar",
  "AT2_Cycling" = "Alveolar",
  "Basal_Progenitor" = "Basal_Lineage",
  "Basal_Cycling" = "Basal_Lineage",
  "Basal_Inflammatory" = "Basal_Lineage",
  "Basal_EMT_ECM" = "Basal_Lineage",
  "Suprabasal_Progenitor" = "Basal_Lineage",
  "Suprabasal_Cycling" = "Basal_Lineage",
  "Ciliated_Mature" = "Ciliated_Lineage",
  "Ciliogenesis_Deuterosomal" = "Ciliated_Lineage",
  "Ciliated_Cycling_Immature" = "Ciliated_Lineage",
  "Goblet" = "Secretory_Lineage",
  "Club" = "Secretory_Lineage",
  "SMG_Mucous" = "Secretory_Lineage",
  "Goblet_Defense_DUOX2" = "Secretory_Lineage",
  "SMG_Serous" = "SMG",
  "SMG_Duct_Secretory_Defense" = "SMG",
  "Squamous_Metaplasia" = "Rare_Specialized",
  "Ionocyte_Brush" = "Rare_Specialized"
)

# ----- Reference Database Paths -----
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"

# ----- Environment Files -----
ENV_FILE_CANDIDATES <- c(
  "/home/h2048/.env",
  "/home/h2048/script/.env"
)

load_env_file <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))

  lines <- readLines(path, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- sub("=.*$", "", line)
    val <- sub("^[^=]*=", "", line)
    key <- trimws(gsub("^export\\s+", "", key))
    val <- trimws(val)
    val <- gsub("^['\"]|['\"]$", "", val)
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
    write(
      sprintf("%s=YOUR_%s_HERE", key, key),
      file = path,
      append = TRUE
    )
    return(invisible(TRUE))
  }

  writeLines(sprintf("%s=YOUR_%s_HERE", key, key), path)
  invisible(TRUE)
}

env_files_loaded <- ENV_FILE_CANDIDATES[file.exists(ENV_FILE_CANDIDATES)]
invisible(lapply(env_files_loaded, load_env_file))

# ----- DeepSeek API (OPTIONAL - pipeline runs without it) -----
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
ENABLE_LLM <- nchar(DEEPSEEK_API_KEY) >= 10
if (!ENABLE_LLM) {
  ensure_env_placeholder("/home/h2048/.env", "DEEPSEEK_API_KEY")
  cat("[WARN] DEEPSEEK_API_KEY not set. interpret_agent will be SKIPPED.\n")
  cat("       Placeholder written to /home/h2048/.env for convenience.\n")
} else if (length(env_files_loaded) > 0) {
  cat(sprintf("[OK] Loaded environment file(s): %s\n", paste(env_files_loaded, collapse = ", ")))
}
INTERPRET_AGENT_MODEL <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS <- 50
INTERPRET_AGENT_ADD_PPI <- TRUE
INTERPRET_AGENT_MAX_RETRIES <- 3L
INTERPRET_AGENT_RETRY_SLEEP_SEC <- 2
INTERPRET_AGENT_MIN_TERMS <- 3L
INTERPRET_AGENT_DISABLE_PPI_ON_STRINGDB_TIMEOUT <- TRUE
INTERPRET_AGENT_DB_PRIORITY <- c(
  "GO_BP",
  "Hallmark",
  "KEGG",
  "GO_MF",
  "GO_CC",
  "CellMarker",
  "PanglaoDB",
  "Epithelial_custom"
)
STANDARDIZE_LLM_OUTPUT <- TRUE
STANDARDIZE_LLM_MODEL <- "deepseek-chat"
STANDARDIZE_LLM_MAX_RETRIES <- 3L
STANDARDIZE_LLM_RETRY_SLEEP_SEC <- 2L
STANDARDIZE_LLM_MAX_INPUT_CHARS <- 12000L

# ----- Python / reticulate -----
PYTHON_CONDA_ENV <- "bbknn_env"

# ----- Metadata Column Names -----
TISSUE_COL <- "tissue"
SAMPLE_COL <- "sample"
CELLTYPE_L2_COL <- "cell_type_L2"
CELLTYPE_L3_COL <- "cell_type_L3"
LABEL_COL <- "cell_type_L3"

UMAP_REDUCTION_PREFERRED <- c(
  "umap_fine",
  "umap_major",
  "umap_scanvi_fine",
  "umap_scanvi_major",
  "umap_scanvi",
  "umap_harmony",
  "umap"
)
UMAP_PT_SIZE <- 0.35
UMAP_TISSUE_COLORS <- c(
  "lung parenchyma" = "#D55E00",
  "nose" = "#009E73",
  "respiratory airway" = "#0072B2",
  "sinus" = "#CC79A7"
)

# ----- Pseudobulk DE Parameters -----
MIN_CELLS_PER_PSEUDOBULK <- 15
MIN_SAMPLES_PER_TISSUE <- 3
PADJ_THR <- 0.05
LFC_THR <- 1.0
DESIGN_COVARIATE <- NULL # e.g. "dataset"
PSEUDOBULK_DISAMBIGUATION_CANDIDATES <- c("dataset", "source", "batch", "orig.ident")

# ----- Reproducibility -----
set.seed(42)

# ----- Exploratory Wilcox Parameters -----
WILCOX_MIN_CELLS <- 100
WILCOX_LFC_THR <- 0.25
WILCOX_TOP_N_PER_DIRECTION <- 100
WILCOX_MAX_CELLS_PER_IDENT <- 10000

# ----- Enrichment -----
TOP_N_DEG_ENRICHMENT <- 200

# ----- Visualization -----
HEATMAP_CELLS_PER_TYPE <- 120

# ----- Epithelial Known Markers -----
KNOWN_MARKERS <- c(
  # Epithelial core
  "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19",
  # Alveolar
  "AGER", "HOPX", "CAV1", "SFTPC", "SFTPA1", "SFTPA2", "SFTPB",
  "ABCA3", "NAPSA", "SLC34A2", "CHI3L1",
  # Basal / suprabasal
  "TP63", "KRT5", "KRT14", "KRT15", "ITGA6", "NGFR", "KRT4", "KRT13",
  "KRT17", "FN1", "COL17A1", "VIM",
  # Ciliated
  "FOXJ1", "TPPP3", "DNAH5", "DNAH9", "RSPH1", "DEUP1", "CCNO", "MCIDAS",
  # Secretory / goblet / club
  "SCGB1A1", "SCGB3A1", "SCGB3A2", "SPDEF", "FOXA3", "MUC5AC", "MUC5B",
  "FCGBP", "TFF3", "DUOX2", "DUOXA2", "LCN2", "BPIFA2", "CEACAM5",
  # SMG
  "LTF", "LYZ", "SLPI", "DMBT1", "BPIFA1", "AZGP1", "WFDC2", "PIGR", "TCN1",
  # Rare / special
  "SPRR1A", "SPRR2A", "IVL", "KRT6A", "S100A7", "FOXI1", "ASCL3", "CFTR",
  "ATP6V0D2", "CLCNKA", "CLCNKB", "BSND",
  # Cycling / activation
  "MKI67", "TOP2A", "UBE2C", "BIRC5", "AURKB", "CENPA", "CCNB1"
)

# ----- Custom Epithelial Marker Database (for enrichment) -----
EPITHELIAL_MARKERS_DB <- data.frame(
  subtype = c(
    "AT1_Canonical",
    "AT1_MatrixRemodeling",
    "AT2",
    "AT2_Cycling",
    "Basal_Progenitor",
    "Basal_Cycling",
    "Basal_Inflammatory",
    "Basal_EMT_ECM",
    "Suprabasal_Progenitor",
    "Suprabasal_Cycling",
    "Ciliated_Mature",
    "Ciliogenesis_Deuterosomal",
    "Ciliated_Cycling_Immature",
    "Goblet",
    "Club",
    "SMG_Mucous",
    "Goblet_Defense_DUOX2",
    "SMG_Serous",
    "SMG_Duct_Secretory_Defense",
    "Squamous_Metaplasia",
    "Ionocyte_Brush"
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

# ----- Tissue-Pair-Specific Context for LLM -----
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

build_tissue_pair_context <- function(t1, t2, ct_l2, dir_name) {
  ctx1 <- TISSUE_CONTEXT[[t1]]
  ctx2 <- TISSUE_CONTEXT[[t2]]
  if (is.null(ctx1)) ctx1 <- paste(t1, ": no specific epithelial context available.")
  if (is.null(ctx2)) ctx2 <- paste(t2, ": no specific epithelial context available.")

  sprintf(
    paste(
      "%s",
      "\nComparing %s vs %s for %s cells.",
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
    t2,
    t1,
    ct_l2,
    dir_name,
    t2,
    t1,
    t2,
    ctx2,
    t1,
    ctx1
  )
}

# ==============================================================================
# 2. Load Libraries
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Epithelial Tissue Comparison Pipeline v1.0.0 (PRODUCTION)\n")
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
})

if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
  py_bin <- Sys.getenv("RETICULATE_PYTHON")
} else {
  use_condaenv(PYTHON_CONDA_ENV, required = TRUE)
  py_bin <- tryCatch(py_config()$python, error = function(e) {
    paste0("conda:", PYTHON_CONDA_ENV)
  })
}

LOCAL_GETSEURAT <- "/home/h2048/script/R/GetSeurat.R"
if (file.exists(LOCAL_GETSEURAT)) {
  source(LOCAL_GETSEURAT)
}

if (ENABLE_LLM) {
  fanyi::set_translate_option(key = DEEPSEEK_API_KEY, source = "deepseek")
  cat("[OK] LLM enabled (DeepSeek)\n")
} else {
  cat("[INFO] LLM disabled (no API key)\n")
}
cat(sprintf("[OK] Python configured: %s\n", py_bin))
if (file.exists(LOCAL_GETSEURAT)) {
  cat(sprintf("[OK] Local GetSeurat loaded: %s\n", LOCAL_GETSEURAT))
}
cat("[OK] Libraries loaded\n\n")

# ==============================================================================
# 3. Load GMT + Reference Databases
# ==============================================================================

cat("=== Loading Reference Databases ===\n")

gmt_all <- tryCatch(
  {
    g <- read.gmt(MSIGDB_GMT_PATH)
    g$gene <- toupper(g$gene)
    cat(sprintf("[OK] MSigDB: %d gene sets\n", length(unique(g$term))))
    g
  },
  error = function(e) {
    cat(sprintf("[WARN] MSigDB failed: %s\n", e$message))
    NULL
  }
)

hallmark_t2g <- if (!is.null(gmt_all)) {
  h <- gmt_all %>% filter(grepl("^HALLMARK_", term))
  if (nrow(h) > 0) h else NULL
} else {
  NULL
}
if (!is.null(hallmark_t2g)) {
  cat(sprintf("[OK] Hallmark: %d\n", length(unique(hallmark_t2g$term))))
}

kegg_t2g <- if (!is.null(gmt_all)) {
  k <- gmt_all %>% filter(grepl("^KEGG_", term))
  if (nrow(k) > 0) k else NULL
} else {
  NULL
}
if (!is.null(kegg_t2g)) {
  cat(sprintf("[OK] KEGG: %d\n", length(unique(kegg_t2g$term))))
}

go_bp_t2g <- NULL
go_mf_t2g <- NULL
go_cc_t2g <- NULL
tryCatch(
  {
    gmt_go <- read.gmt(GMT_GO_ALL)
    gmt_go$gene <- toupper(gmt_go$gene)
    go_bp_t2g <- gmt_go %>% filter(grepl("^GOBP_", term))
    go_mf_t2g <- gmt_go %>% filter(grepl("^GOMF_", term))
    go_cc_t2g <- gmt_go %>% filter(grepl("^GOCC_", term))
    cat(sprintf(
      "[OK] GO: BP=%d, MF=%d, CC=%d\n",
      length(unique(go_bp_t2g$term)),
      length(unique(go_mf_t2g$term)),
      length(unique(go_cc_t2g$term))
    ))
  },
  error = function(e) cat(sprintf("[WARN] GO GMT failed: %s\n", e$message))
)

cellmarker_t2g <- NULL
tryCatch(
  {
    cm_db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
    cm_db <- cm_db[grepl("Human", species, ignore.case = TRUE)]
    cat(sprintf("[OK] CellMarker raw: %d human entries\n", nrow(cm_db)))

    epi_pattern <- paste(
      c(
        "Epithelial", "Basal", "Goblet", "Club", "Secretory", "Ciliated",
        "Ionocyte", "Brush", "Tuft", "Deuterosomal", "Serous", "Duct",
        "Airway", "Respiratory", "Alveolar", "AT1", "AT2", "Suprabasal", "Squamous"
      ),
      collapse = "|"
    )
    airway_pattern <- "Lung|Airway|Respiratory|Nasal|Sinus|Bronch|Trachea|Alveol"

    cm_epi <- cm_db[
      grepl(
        epi_pattern,
        cell_name,
        ignore.case = TRUE
      ) &
      grepl(
        airway_pattern,
        tissue_type,
        ignore.case = TRUE
      )
    ]
    cat(sprintf("[OK] CellMarker epithelial-related: %d entries\n", nrow(cm_epi)))

    if (nrow(cm_epi) < 10) {
      cat("[WARN] Too few epithelial entries in CellMarker; skipping this database\n")
    } else {
      t2g_list <- list()
      for (i in seq_len(nrow(cm_epi))) {
        markers_raw <- cm_epi$marker[i]
        cell_type <- cm_epi$cell_name[i]
        if (is.na(markers_raw) || markers_raw == "") next
        markers <- unlist(strsplit(markers_raw, "[,;\\s]+"))
        markers <- gsub('["\r\n\\[\\]]', '', markers)
        markers <- trimws(toupper(markers))
        markers <- unique(markers[markers != "" & !is.na(markers)])
        if (length(markers) > 0) {
          t2g_list[[length(t2g_list) + 1]] <- data.frame(
            term = cell_type,
            gene = markers,
            stringsAsFactors = FALSE
          )
        }
      }
      cellmarker_t2g <- bind_rows(t2g_list) %>% distinct()
      cat(sprintf(
        "[OK] CellMarker TERM2GENE: %d pairs (%d epithelial terms)\n",
        nrow(cellmarker_t2g),
        length(unique(cellmarker_t2g$term))
      ))
    }
  },
  error = function(e) cat(sprintf("[WARN] CellMarker failed: %s\n", e$message))
)

panglaodb_t2g <- NULL
tryCatch(
  {
    pdb <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)
    setnames(
      pdb,
      old = c("official gene symbol", "cell type"),
      new = c("gene_symbol", "cell_type"),
      skip_absent = TRUE
    )
    pdb <- pdb[grepl("Hs", species, fixed = TRUE)]
    cat(sprintf("[OK] PanglaoDB raw: %d human markers\n", nrow(pdb)))

    pdb_epi <- pdb[
      grepl(
        paste(
          c(
            "Epithelial", "Basal", "Goblet", "Club", "Secretory", "Ciliated",
            "Ionocyte", "Brush", "Tuft", "Serous", "Duct", "Airway", "Alveolar",
            "AT1", "AT2", "Suprabasal", "Squamous"
          ),
          collapse = "|"
        ),
        cell_type,
        ignore.case = TRUE
      )
    ]
    cat(sprintf(
      "[OK] PanglaoDB epithelial-related: %d markers (%d cell types)\n",
      nrow(pdb_epi),
      length(unique(pdb_epi$cell_type))
    ))

    if (nrow(pdb_epi) < 5) {
      cat("[WARN] Too few epithelial markers in PanglaoDB; skipping this database\n")
    } else {
      panglaodb_t2g <- pdb_epi %>%
        select(cell_type, gene_symbol) %>%
        mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
        filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
        distinct() %>%
        rename(term = cell_type, gene = gene_symbol)
      cat(sprintf(
        "[OK] PanglaoDB TERM2GENE: %d pairs (%d epithelial terms)\n",
        nrow(panglaodb_t2g),
        length(unique(panglaodb_t2g$term))
      ))
    }
  },
  error = function(e) cat(sprintf("[WARN] PanglaoDB failed: %s\n", e$message))
)

epithelial_custom_t2g <- EPITHELIAL_MARKERS_DB %>%
  separate_rows(markers, sep = ",") %>%
  mutate(markers = trimws(toupper(markers))) %>%
  filter(markers != "") %>%
  select(term = subtype, gene = markers) %>%
  distinct()
cat(sprintf(
  "[OK] Custom epithelial markers: %d subtypes, %d pairs\n",
  length(unique(epithelial_custom_t2g$term)),
  nrow(epithelial_custom_t2g)
))
cat("\n")

# ==============================================================================
# 4. Helper Functions
# ==============================================================================

safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

normalize_meta_values <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "<NA>"
  x
}

pick_reduction <- function(obj, preferred) {
  red <- Reductions(obj)
  hit <- preferred[preferred %in% red]
  if (length(hit) > 0) hit[[1]] else NULL
}

build_umap_plot <- function(
  obj,
  reduction_name,
  group_col,
  title,
  split_col = NULL,
  label = FALSE,
  width = 10,
  height = 8,
  cols = NULL
) {
  p <- DimPlot(
    obj,
    reduction = reduction_name,
    group.by = group_col,
    split.by = split_col,
    pt.size = UMAP_PT_SIZE,
    shuffle = TRUE,
    label = label,
    repel = label,
    cols = cols
  ) +
    ggtitle(title) +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  list(plot = p, width = width, height = height)
}

run_gmt_enrichment <- function(gene_list, t2g, db_name, tested_genes = NULL) {
  if (is.null(t2g) || nrow(t2g) == 0 || length(gene_list) < 5) return(NULL)
  if (!is.null(tested_genes)) {
    universe_use <- intersect(toupper(tested_genes), unique(t2g$gene))
  } else {
    universe_use <- unique(t2g$gene)
  }
  tryCatch(
    enricher(
      gene = toupper(unique(gene_list)),
      TERM2GENE = t2g,
      universe = universe_use,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    ),
    error = function(e) {
      cat(sprintf("    [WARN] %s: %s\n", db_name, e$message))
      NULL
    }
  )
}

find_ambiguous_sample_tissue_keys <- function(
  meta_sub,
  resolved_cols = character(),
  candidate_cols = character()
) {
  compare_cols <- unique(c(resolved_cols, candidate_cols))
  if (length(compare_cols) == 0) return(character())

  key_cols <- unique(c(SAMPLE_COL, TISSUE_COL, compare_cols))
  key_df <- unique(meta_sub[, key_cols, drop = FALSE])
  key_df[] <- lapply(key_df, normalize_meta_values)

  resolved_key <- do.call(paste, c(key_df[, unique(c(SAMPLE_COL, TISSUE_COL, resolved_cols)), drop = FALSE], sep = "__"))
  base_key <- paste(key_df[[SAMPLE_COL]], key_df[[TISSUE_COL]], sep = "__")
  unique(base_key[duplicated(resolved_key) | duplicated(resolved_key, fromLast = TRUE)])
}

resolve_pseudobulk_group_columns <- function(meta_sub) {
  group_cols <- c(SAMPLE_COL, TISSUE_COL)
  candidate_cols <- character()

  if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
    candidate_cols <- c(candidate_cols, DESIGN_COVARIATE)
  }

  candidate_cols <- unique(c(candidate_cols, PSEUDOBULK_DISAMBIGUATION_CANDIDATES))
  candidate_cols <- setdiff(candidate_cols, c(SAMPLE_COL, TISSUE_COL))
  candidate_cols <- candidate_cols[candidate_cols %in% colnames(meta_sub)]

  ambiguous_keys <- find_ambiguous_sample_tissue_keys(meta_sub, candidate_cols = candidate_cols)
  disambiguation_cols <- character()

  while (length(ambiguous_keys) > 0 && length(candidate_cols) > 0) {
    informative_cols <- candidate_cols[vapply(candidate_cols, function(col) {
      probe_df <- unique(meta_sub[, c(SAMPLE_COL, TISSUE_COL, col), drop = FALSE])
      probe_df[] <- lapply(probe_df, normalize_meta_values)
      probe_key <- paste(probe_df[[SAMPLE_COL]], probe_df[[TISSUE_COL]], sep = "__")
      probe_df <- probe_df[probe_key %in% ambiguous_keys, , drop = FALSE]
      probe_key <- probe_key[probe_key %in% ambiguous_keys]
      if (nrow(probe_df) == 0) return(FALSE)
      any(vapply(split(probe_df[[col]], probe_key), function(vals) length(unique(vals)) > 1, logical(1)))
    }, logical(1))]

    if (length(informative_cols) == 0) break

    chosen_col <- informative_cols[[1]]
    group_cols <- c(group_cols, chosen_col)
    disambiguation_cols <- c(disambiguation_cols, chosen_col)
    candidate_cols <- setdiff(candidate_cols, chosen_col)
    ambiguous_keys <- find_ambiguous_sample_tissue_keys(
      meta_sub,
      resolved_cols = disambiguation_cols,
      candidate_cols = setdiff(candidate_cols, disambiguation_cols)
    )
  }

  list(
    group_cols = unique(group_cols),
    disambiguation_cols = unique(disambiguation_cols),
    ambiguous_keys_remaining = ambiguous_keys
  )
}

aggregate_pseudobulk <- function(obj, celltype_l2) {
  cells <- colnames(obj)[obj@meta.data[[CELLTYPE_L2_COL]] == celltype_l2]
  if (length(cells) < 50) return(NULL)

  sub <- subset(obj, cells = cells)
  meta_sub <- sub@meta.data
  counts_mat <- GetAssayData(sub, layer = "counts")

  resolved_grouping <- resolve_pseudobulk_group_columns(meta_sub)
  if (length(resolved_grouping$ambiguous_keys_remaining) > 0) {
    stop(sprintf(
      paste(
        "Non-unique sample+tissue pseudobulk keys remain for cell type '%s'.",
        "Please provide a globally unique sample identifier or set DESIGN_COVARIATE / metadata disambiguation columns.",
        "Ambiguous keys: %s"
      ),
      celltype_l2,
      paste(utils::head(resolved_grouping$ambiguous_keys_remaining, 10), collapse = ", ")
    ))
  }

  group_cols <- resolved_grouping$group_cols
  extra_meta_cols <- setdiff(group_cols, c(SAMPLE_COL, TISSUE_COL))
  meta_sub$pb_group <- apply(meta_sub[, group_cols, drop = FALSE], 1, function(row_vals) {
    paste(normalize_meta_values(row_vals), collapse = "__")
  })
  groups <- unique(meta_sub$pb_group)

  pb_cols <- list()
  pb_meta <- data.frame(row.names = character(), stringsAsFactors = FALSE)

  for (g in groups) {
    cells_g <- rownames(meta_sub)[meta_sub$pb_group == g]
    if (length(cells_g) < MIN_CELLS_PER_PSEUDOBULK) next

    col_sum <- if (length(cells_g) == 1) {
      counts_mat[, cells_g, drop = FALSE]
    } else {
      Matrix::rowSums(counts_mat[, cells_g, drop = FALSE])
    }

    if (is.numeric(col_sum) && is.null(dim(col_sum))) {
      col_sum <- Matrix::Matrix(col_sum, ncol = 1, sparse = TRUE)
      rownames(col_sum) <- rownames(counts_mat)
    }
    colnames(col_sum) <- g
    pb_cols[[g]] <- col_sum

    first_idx <- which(meta_sub$pb_group == g)[1]
    pb_meta[g, "tissue"] <- meta_sub[[TISSUE_COL]][first_idx]
    pb_meta[g, "sample"] <- meta_sub[[SAMPLE_COL]][first_idx]
    pb_meta[g, "n_cells"] <- length(cells_g)
    if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
      pb_meta[g, DESIGN_COVARIATE] <- meta_sub[[DESIGN_COVARIATE]][first_idx]
    }
    if (length(extra_meta_cols) > 0) {
      for (extra_col in extra_meta_cols) {
        pb_meta[g, extra_col] <- meta_sub[[extra_col]][first_idx]
      }
    }
  }

  if (length(pb_cols) < 4) return(NULL)

  pb_counts <- as.matrix(do.call(cbind, pb_cols))
  valid <- colnames(pb_counts)[
    colSums(pb_counts) > 0 & !is.na(pb_meta[colnames(pb_counts), "tissue"])
  ]
  if (length(valid) < 4) return(NULL)

  list(
    counts = pb_counts[, valid, drop = FALSE],
    meta = pb_meta[valid, , drop = FALSE],
    group_cols = group_cols,
    disambiguation_cols = resolved_grouping$disambiguation_cols
  )
}

run_deseq2_pairwise <- function(pb, t1, t2) {
  keep <- pb$meta$tissue %in% c(t1, t2)
  if (sum(keep) < 4) return(NULL)

  counts_sub <- pb$counts[, keep, drop = FALSE]
  meta_sub <- pb$meta[keep, , drop = FALSE]
  tissue_levels_safe <- make.names(c(t1, t2), unique = TRUE)
  meta_sub$tissue <- droplevels(factor(meta_sub$tissue, levels = c(t1, t2)))
  meta_sub$tissue_safe <- factor(
    ifelse(meta_sub$tissue == t1, tissue_levels_safe[1], tissue_levels_safe[2]),
    levels = tissue_levels_safe
  )

  n1 <- sum(meta_sub$tissue == t1)
  n2 <- sum(meta_sub$tissue == t2)
  if (n1 < MIN_SAMPLES_PER_TISSUE || n2 < MIN_SAMPLES_PER_TISSUE) return(NULL)

  keep_genes <- rowSums(counts_sub >= 1) >= max(3, ncol(counts_sub) * 0.2)
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  if (nrow(counts_sub) < 100) return(NULL)

  tested_genes <- rownames(counts_sub)
  counts_int <- round(counts_sub)
  storage.mode(counts_int) <- "integer"

  design_formula <- ~tissue_safe
  if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
    n_levels <- length(unique(meta_sub[[DESIGN_COVARIATE]]))
    if (n_levels >= 2 && n_levels < nrow(meta_sub)) {
      meta_sub[[DESIGN_COVARIATE]] <- factor(meta_sub[[DESIGN_COVARIATE]])
      design_formula <- as.formula(paste("~", DESIGN_COVARIATE, "+ tissue_safe"))
      cat(sprintf("    Design: %s\n", deparse(design_formula)))
    }
  }

  dds <- tryCatch(
    {
      dds <- DESeqDataSetFromMatrix(counts_int, meta_sub, design_formula)
      DESeq(dds, quiet = TRUE)
    },
    error = function(e) {
      if (!is.null(DESIGN_COVARIATE)) {
        cat(sprintf("    [WARN] Design failed (%s), fallback to ~ tissue_safe\n", e$message))
        tryCatch(
          {
            dds2 <- DESeqDataSetFromMatrix(counts_int, meta_sub, ~tissue_safe)
            DESeq(dds2, quiet = TRUE)
          },
          error = function(e2) {
            cat(sprintf("    [ERROR] DESeq2: %s\n", e2$message))
            NULL
          }
        )
      } else {
        cat(sprintf("    [ERROR] DESeq2: %s\n", e$message))
        NULL
      }
    }
  )
  if (is.null(dds)) return(NULL)

  res <- results(
    dds,
    contrast = c("tissue_safe", tissue_levels_safe[2], tissue_levels_safe[1]),
    alpha = PADJ_THR
  )
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(
      sig = ifelse(padj < PADJ_THR & abs(log2FoldChange) > LFC_THR, "sig", "ns"),
      direction = ifelse(log2FoldChange > 0, "up", "down")
    )

  list(
    de_table = res_df,
    tested_genes = tested_genes,
    n_up = sum(res_df$sig == "sig" & res_df$direction == "up"),
    n_down = sum(res_df$sig == "sig" & res_df$direction == "down"),
    n_samples_1 = n1,
    n_samples_2 = n2,
    tissue_1 = t1,
    tissue_2 = t2,
    design = deparse(design_formula)
  )
}

run_wilcox_exploratory <- function(obj, celltype_l2) {
  cells <- colnames(obj)[obj@meta.data[[CELLTYPE_L2_COL]] == celltype_l2]
  if (length(cells) < WILCOX_MIN_CELLS) return(NULL)

  sub <- subset(obj, cells = cells)
  tissues <- sort(unique(na.omit(sub@meta.data[[TISSUE_COL]])))
  if (length(tissues) < 2) return(NULL)

  Idents(sub) <- TISSUE_COL
  results <- list()
  for (pair in combn(tissues, 2, simplify = FALSE)) {
    t1 <- pair[1]
    t2 <- pair[2]
    n1 <- sum(sub@meta.data[[TISSUE_COL]] == t1, na.rm = TRUE)
    n2 <- sum(sub@meta.data[[TISSUE_COL]] == t2, na.rm = TRUE)
    if (n1 < WILCOX_MIN_CELLS || n2 < WILCOX_MIN_CELLS) next

    cells_pair <- rownames(sub@meta.data)[sub@meta.data[[TISSUE_COL]] %in% c(t1, t2)]
    if (length(cells_pair) == 0) next
    sub_pair <- subset(sub, cells = cells_pair)
    Idents(sub_pair) <- TISSUE_COL

    if (n1 > WILCOX_MAX_CELLS_PER_IDENT || n2 > WILCOX_MAX_CELLS_PER_IDENT) {
      set.seed(42)
      keep_t1 <- rownames(sub_pair@meta.data)[sub_pair@meta.data[[TISSUE_COL]] == t1]
      keep_t2 <- rownames(sub_pair@meta.data)[sub_pair@meta.data[[TISSUE_COL]] == t2]
      keep_cells <- c(
        sample(keep_t1, min(length(keep_t1), WILCOX_MAX_CELLS_PER_IDENT)),
        sample(keep_t2, min(length(keep_t2), WILCOX_MAX_CELLS_PER_IDENT))
      )
      sub_pair <- subset(sub_pair, cells = keep_cells)
      Idents(sub_pair) <- TISSUE_COL
      cat(sprintf(
        "  [INFO] Wilcox downsampled %s: %s=%d, %s=%d\n",
        paste0(t2, "_vs_", t1),
        t1,
        sum(sub_pair@meta.data[[TISSUE_COL]] == t1),
        t2,
        sum(sub_pair@meta.data[[TISSUE_COL]] == t2)
      ))
    }

    de <- tryCatch(
      FindMarkers(
        sub_pair,
        ident.1 = t2,
        ident.2 = t1,
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0
      ),
      error = function(e) NULL
    )
    if (!is.null(de) && nrow(de) > 0) {
      de$gene <- rownames(de)
      de$padj <- p.adjust(de$p_val, method = "BH")
      results[[paste0(t2, "_vs_", t1)]] <- de %>%
        filter(abs(avg_log2FC) >= WILCOX_LFC_THR) %>%
        mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
        arrange(direction, padj, desc(abs(avg_log2FC))) %>%
        group_by(direction) %>%
        slice_head(n = WILCOX_TOP_N_PER_DIRECTION) %>%
        ungroup()
    }
  }
  results
}

safe_trim <- function(x) {
  x <- to_scalar(x)
  trimws(x)
}

capture_object_text <- function(x) {
  if (is.null(x)) return("")
  txt <- tryCatch(paste(capture.output(print(x)), collapse = "\n"), error = function(e) "")
  if (!nzchar(trimws(txt))) {
    txt <- tryCatch(paste(capture.output(str(x, max.level = 3)), collapse = "\n"), error = function(e) "")
  }
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
    if (is.list(x) && !is.null(x[[nm]])) {
      val <- x[[nm]]
    } else if (!is.list(x) && has_named_entries(x) && nm %in% names(x)) {
      val <- x[[nm]]
    }
    if (!is.null(val)) {
      txt <- safe_trim(val)
      if (nzchar(txt)) return(txt)
    }
  }
  ""
}

is_retryable_interpret_agent_issue <- function(warnings = character(), error_message = NULL) {
  signals <- c(warnings, error_message)
  if (length(signals) == 0) return(FALSE)
  any(grepl(
    "Failed to parse JSON response|premature EOF|invalid for atomic vectors|429|5[0-9]{2}|timeout|temporar",
    signals,
    ignore.case = TRUE,
    perl = TRUE
  ))
}

has_stringdb_timeout_issue <- function(warnings = character(), error_message = NULL) {
  signals <- c(warnings, error_message)
  if (length(signals) == 0) return(FALSE)
  any(grepl(
    "stringdb-static\\.org|species\\.v12\\.txt",
    signals,
    ignore.case = TRUE,
    perl = TRUE
  ))
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
  x <- safe_trim(x)
  if (nzchar(x)) x else default
}

collapse_driver_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  paste(unique(x), collapse = ", ")
}

truncate_text <- function(x, max_chars = STANDARDIZE_LLM_MAX_INPUT_CHARS) {
  x <- safe_trim(x)
  if (!nzchar(x)) return("")
  if (nchar(x, type = "chars") <= max_chars) return(x)
  paste0(substr(x, 1, max_chars), "\n...[truncated]")
}

build_enrichment_evidence_text <- function(enrich_obj, gene_fc = NULL, n_terms = 8L, n_genes = 15L) {
  if (is.null(enrich_obj)) return("")
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return("")

  top_terms <- utils::head(enrich_df[order(enrich_df$p.adjust), , drop = FALSE], n_terms)
  term_lines <- vapply(
    seq_len(nrow(top_terms)),
    function(i) {
      sprintf(
        "- %s | padj=%s | Count=%s | Genes=%s",
        safe_trim(top_terms$Description[i]),
        format(top_terms$p.adjust[i], scientific = TRUE, digits = 3),
        safe_trim(top_terms$Count[i]),
        truncate_text(gsub("/", ", ", safe_trim(top_terms$geneID[i])), 200)
      )
    },
    character(1)
  )

  gene_lines <- character()
  if (!is.null(gene_fc) && length(gene_fc) > 0) {
    gene_fc <- sort(gene_fc[!is.na(gene_fc)], decreasing = TRUE)
    top_up <- utils::head(gene_fc, n_genes)
    top_down <- utils::head(sort(gene_fc, decreasing = FALSE), n_genes)
    gene_lines <- c(
      sprintf("Top up genes: %s", paste(sprintf("%s(%.2f)", names(top_up), top_up), collapse = ", ")),
      sprintf("Top down genes: %s", paste(sprintf("%s(%.2f)", names(top_down), top_down), collapse = ", "))
    )
  }

  paste(c("Top enrichment evidence:", term_lines, gene_lines), collapse = "\n")
}

select_interpret_agent_databases <- function(enr_list) {
  if (is.null(enr_list) || length(enr_list) == 0) return(character())

  db_order <- c(INTERPRET_AGENT_DB_PRIORITY, setdiff(names(enr_list), INTERPRET_AGENT_DB_PRIORITY))
  db_order <- db_order[db_order %in% names(enr_list)]

  db_order[vapply(db_order, function(db) {
    er <- enr_list[[db]]
    if (is.null(er)) return(FALSE)
    er_df <- tryCatch(as.data.frame(er), error = function(e) NULL)
    !is.null(er_df) && nrow(er_df) >= INTERPRET_AGENT_MIN_TERMS
  }, logical(1))]
}

build_multi_database_evidence_text <- function(enrich_list, gene_fc = NULL, n_terms_per_db = 5L) {
  if (is.null(enrich_list) || length(enrich_list) == 0) return("")

  sections <- lapply(names(enrich_list), function(db) {
    er <- enrich_list[[db]]
    er_df <- tryCatch(as.data.frame(er), error = function(e) NULL)
    if (is.null(er_df) || nrow(er_df) == 0) return(NULL)

    top_terms <- utils::head(er_df[order(er_df$p.adjust), , drop = FALSE], n_terms_per_db)
    term_lines <- vapply(
      seq_len(nrow(top_terms)),
      function(i) {
        sprintf(
          "- %s | padj=%s | Count=%s | Genes=%s",
          safe_trim(top_terms$Description[i]),
          format(top_terms$p.adjust[i], scientific = TRUE, digits = 3),
          safe_trim(top_terms$Count[i]),
          truncate_text(gsub("/", ", ", safe_trim(top_terms$geneID[i])), 180)
        )
      },
      character(1)
    )
    c(sprintf("[%s]", db), term_lines)
  })
  sections <- Filter(Negate(is.null), sections)

  gene_lines <- character()
  if (!is.null(gene_fc) && length(gene_fc) > 0) {
    gene_fc <- gene_fc[!is.na(gene_fc)]
    abs_rank <- sort(abs(gene_fc), decreasing = TRUE)
    top_abs_genes <- names(utils::head(abs_rank, 20))
    gene_lines <- sprintf(
      "Top |log2FC| genes: %s",
      paste(sprintf("%s(%.2f)", top_abs_genes, gene_fc[top_abs_genes]), collapse = ", ")
    )
  }

  paste(c("Integrated multi-database evidence:", unlist(sections), gene_lines), collapse = "\n")
}

build_multi_database_raw_text <- function(per_db_payloads) {
  if (length(per_db_payloads) == 0) return("")

  blocks <- lapply(names(per_db_payloads), function(db) {
    payload <- per_db_payloads[[db]]
    raw_text <- capture_object_text(payload$result)
    c(
      sprintf("=== %s ===", db),
      sprintf("error: %s", ifelse(is.null(payload$error) || !nzchar(payload$error), "", payload$error)),
      sprintf("warnings: %s", paste(payload$warnings, collapse = " | ")),
      if (nzchar(raw_text)) raw_text else "<empty interpret_agent output>"
    )
  })

  paste(unlist(blocks), collapse = "\n")
}

run_interpret_agent_multi <- function(enr_list, context_str, dir_name, gene_fc = NULL) {
  selected_dbs <- select_interpret_agent_databases(enr_list)
  if (length(selected_dbs) == 0) {
    return(list(
      selected_dbs = character(),
      per_db_payloads = list(),
      combined_raw_text = "",
      combined_warnings = "",
      combined_error = "no enrichment database met interpret_agent criteria"
    ))
  }

  per_db_payloads <- list()
  combined_warnings <- character()
  combined_errors <- character()
  for (db in selected_dbs) {
    cat(sprintf("    interpret_agent [%s] (db: %s)\n", dir_name, db))
    db_context <- paste(
      context_str,
      sprintf("Primary enrichment database for this pass: %s.", db),
      sprintf("Integrate with the broader multi-database evidence, but prioritize terms from %s in this pass.", db),
      sep = "\n\n"
    )
    payload <- run_interpret_agent_safe(enr_list[[db]], db_context, gene_fc)
    per_db_payloads[[db]] <- payload
    if (length(payload$warnings) > 0) {
      combined_warnings <- c(combined_warnings, sprintf("%s: %s", db, paste(payload$warnings, collapse = " | ")))
    }
    if (!is.null(payload$error) && nzchar(payload$error)) {
      combined_errors <- c(combined_errors, sprintf("%s: %s", db, payload$error))
    }
  }

  list(
    selected_dbs = selected_dbs,
    per_db_payloads = per_db_payloads,
    combined_raw_text = build_multi_database_raw_text(per_db_payloads),
    combined_warnings = unique(combined_warnings),
    combined_error = paste(unique(combined_errors), collapse = " | ")
  )
}

build_multi_database_interpretation_record <- function(
  multi_payload,
  enrich_list,
  ct_l2,
  comp_name,
  dir_name,
  context_str,
  gene_fc = NULL
) {
  selected_dbs <- multi_payload$selected_dbs
  selected_enrich <- enrich_list[selected_dbs]
  combined_raw_text <- safe_trim(multi_payload$combined_raw_text)
  combined_warnings <- unique(multi_payload$combined_warnings)
  combined_error <- safe_trim(multi_payload$combined_error)
  source_db <- paste(selected_dbs, collapse = "; ")
  evidence_text <- build_multi_database_evidence_text(selected_enrich, gene_fc = gene_fc)

  standardized_payload <- standardize_result_with_llm(
    raw_text = combined_raw_text,
    context_str = paste(
      context_str,
      sprintf("Integrate evidence across ALL eligible databases: %s.", source_db),
      "When databases agree, state the consensus clearly. When they differ, explain what each database contributes and keep the conclusion conservative.",
      sep = "\n\n"
    ),
    evidence_text = evidence_text,
    ct_l2 = ct_l2,
    comp_name = comp_name,
    dir_name = dir_name,
    source_db = source_db,
    warnings = combined_warnings,
    error_message = if (nzchar(combined_error)) combined_error else NULL
  )

  status <- if (!is.null(standardized_payload$result)) {
    "structured"
  } else if (nzchar(combined_raw_text)) {
    "raw_text"
  } else if (nzchar(combined_error)) {
    "error"
  } else {
    "empty"
  }

  std <- standardized_payload$result
  list(
    celltype_l2 = ct_l2,
    comparison = comp_name,
    direction = dir_name,
    source_db = source_db,
    status = status,
    warnings = unique(c(combined_warnings, standardized_payload$warnings)),
    error = ifelse(
      !is.null(standardized_payload$error) && nzchar(standardized_payload$error),
      standardized_payload$error,
      combined_error
    ),
    overview = placeholder_text(if (!is.null(std)) std$overview else ""),
    key_mechanisms = placeholder_text(if (!is.null(std)) std$key_mechanisms else ""),
    hypothesis = placeholder_text(if (!is.null(std)) std$hypothesis else ""),
    narrative = placeholder_text(if (!is.null(std)) std$narrative else combined_raw_text),
    key_drivers = placeholder_text(if (!is.null(std)) collapse_driver_field(std$key_drivers) else ""),
    evidence = placeholder_text(if (!is.null(std)) std$evidence else evidence_text),
    limitations = placeholder_text(if (!is.null(std)) std$limitations else "Interpretation could not be standardized across databases."),
    raw_text = combined_raw_text,
    raw_result = multi_payload
  )
}

extract_json_string <- function(text) {
  text <- safe_trim(text)
  if (!nzchar(text)) return("")
  text <- gsub("^```(?:json)?\\s*", "", text, perl = TRUE)
  text <- gsub("\\s*```$", "", text, perl = TRUE)
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

standardize_result_with_llm <- function(
  raw_text,
  context_str,
  evidence_text,
  ct_l2,
  comp_name,
  dir_name,
  source_db,
  warnings = character(),
  error_message = NULL
) {
  if (!ENABLE_LLM || !STANDARDIZE_LLM_OUTPUT) {
    return(list(result = NULL, warnings = character(), error = "standardization disabled"))
  }

  prompt <- paste(
    "You are standardizing a biological interpretation for single-cell epithelial tissue comparison.",
    "Return valid JSON only. No markdown, no code fences, no extra commentary.",
    "Use exactly these keys:",
    "overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
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
    sprintf("Cell type L2: %s", ct_l2),
    sprintf("Comparison: %s", comp_name),
    sprintf("Direction: %s", dir_name),
    sprintf("Source DB: %s", source_db),
    "Context:",
    truncate_text(context_str),
    "Original interpret_agent output:",
    truncate_text(raw_text),
    "Enrichment evidence:",
    truncate_text(evidence_text),
    "Warnings and errors:",
    truncate_text(paste(c(warnings, error_message), collapse = "\n")),
    sep = "\n"
  )

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(STANDARDIZE_LLM_MAX_RETRIES)) {
    response_text <- tryCatch(
      safe_trim(fanyi::chat_request(
        prompt,
        model = STANDARDIZE_LLM_MODEL,
        api_key = DEEPSEEK_API_KEY
      )),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )

    parsed <- parse_standardized_json(response_text)
    if (!is.null(parsed)) {
      return(list(result = parsed, warnings = unique(collected_warnings), error = NULL))
    }

    if (nzchar(response_text)) {
      collected_warnings <- c(collected_warnings, sprintf("standardize attempt %d returned non-JSON text", attempt))
    }
    if (attempt < STANDARDIZE_LLM_MAX_RETRIES) Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
  }

  list(
    result = NULL,
    warnings = unique(collected_warnings),
    error = ifelse(is.null(last_error) || !nzchar(last_error), "failed to standardize interpret_agent output", last_error)
  )
}

normalize_interpret_agent_result <- function(
  raw_result,
  ct_l2,
  comp_name,
  dir_name,
  source_db = NA_character_,
  warnings = character(),
  error_message = NULL,
  enrich_obj = NULL,
  context_str = NULL,
  gene_fc = NULL
) {
  core <- unwrap_interpret_agent_result(raw_result)
  raw_text <- capture_object_text(raw_result)

  overview <- extract_named_text(core, c("overview", "summary", "interpretation", "narrative"))
  key_mechanisms <- extract_named_text(core, c("key_mechanisms", "mechanisms", "keyMechanisms"))
  hypothesis <- extract_named_text(core, c("hypothesis", "model", "working_hypothesis"))
  narrative <- extract_named_text(core, c("narrative", "story", "details"))
  key_drivers <- extract_named_text(core, c("key_drivers", "drivers", "genes", "gene_drivers"))
  evidence <- extract_named_text(core, c("evidence", "supporting_evidence", "rationale"))
  limitations <- extract_named_text(core, c("limitations", "caveats", "notes"))

  if (!nzchar(overview) && nzchar(raw_text)) overview <- raw_text

  overview_raw <- overview
  key_mechanisms_raw <- key_mechanisms
  hypothesis_raw <- hypothesis
  narrative_raw <- narrative
  key_drivers_raw <- key_drivers
  evidence_raw <- evidence
  limitations_raw <- limitations

  major_fields_present <- c(overview, key_mechanisms, hypothesis, narrative, evidence)
  needs_standardization <- isTRUE(STANDARDIZE_LLM_OUTPUT) &&
    (!all(nzchar(major_fields_present)) || !nzchar(key_drivers) || !is.null(error_message) || !is.list(core))

  if (needs_standardization) {
    evidence_text <- build_enrichment_evidence_text(enrich_obj, gene_fc = gene_fc)
    standardized_payload <- standardize_result_with_llm(
      raw_text = raw_text,
      context_str = context_str,
      evidence_text = evidence_text,
      ct_l2 = ct_l2,
      comp_name = comp_name,
      dir_name = dir_name,
      source_db = source_db,
      warnings = warnings,
      error_message = error_message
    )

    if (!is.null(standardized_payload$result)) {
      std <- standardized_payload$result
      overview <- placeholder_text(std$overview)
      key_mechanisms <- placeholder_text(std$key_mechanisms)
      hypothesis <- placeholder_text(std$hypothesis)
      narrative <- placeholder_text(std$narrative)
      key_drivers <- placeholder_text(collapse_driver_field(std$key_drivers))
      evidence <- placeholder_text(std$evidence)
      limitations <- placeholder_text(std$limitations)
      warnings <- unique(c(
        warnings,
        standardized_payload$warnings,
        if (!is.null(error_message) && nzchar(error_message)) sprintf("initial interpret_agent error: %s", error_message)
      ))
      error_message <- NULL
    } else {
      warnings <- unique(c(warnings, standardized_payload$warnings))
      limitations <- placeholder_text(limitations)
      if (!is.null(standardized_payload$error) && nzchar(standardized_payload$error)) {
        error_message <- paste(unique(c(error_message, standardized_payload$error)), collapse = " | ")
      }
    }
  }

  status <- dplyr::case_when(
    !is.null(error_message) ~ "error",
    is.null(raw_result) ~ "missing",
    any(nzchar(c(
      overview_raw,
      key_mechanisms_raw,
      hypothesis_raw,
      narrative_raw,
      key_drivers_raw,
      evidence_raw,
      limitations_raw
    ))) ~ "structured",
    nzchar(raw_text) ~ "raw_text",
    TRUE ~ "empty"
  )

  overview <- placeholder_text(overview)
  key_mechanisms <- placeholder_text(key_mechanisms)
  hypothesis <- placeholder_text(hypothesis)
  narrative <- placeholder_text(narrative)
  key_drivers <- placeholder_text(key_drivers)
  evidence <- placeholder_text(evidence)
  limitations <- placeholder_text(limitations)

  list(
    celltype_l2 = ct_l2,
    comparison = comp_name,
    direction = dir_name,
    source_db = source_db,
    status = status,
    warnings = unique(warnings),
    error = if (is.null(error_message)) "" else error_message,
    overview = overview,
    key_mechanisms = key_mechanisms,
    hypothesis = hypothesis,
    narrative = narrative,
    key_drivers = key_drivers,
    evidence = evidence,
    limitations = limitations,
    raw_text = raw_text,
    raw_result = raw_result
  )
}

flatten_interpretation_records <- function(x) {
  rows <- list()
  idx <- 1L
  for (ct_l2 in names(x)) {
    for (comp_name in names(x[[ct_l2]])) {
      for (dir_name in names(x[[ct_l2]][[comp_name]])) {
        rec <- x[[ct_l2]][[comp_name]][[dir_name]]
        if (is.null(rec)) next
        rows[[idx]] <- data.frame(
          celltype_l2 = safe_trim(rec$celltype_l2),
          comparison = safe_trim(rec$comparison),
          direction = safe_trim(rec$direction),
          source_db = safe_trim(rec$source_db),
          status = safe_trim(rec$status),
          warnings = paste(rec$warnings, collapse = " | "),
          error = safe_trim(rec$error),
          overview = safe_trim(rec$overview),
          key_mechanisms = safe_trim(rec$key_mechanisms),
          hypothesis = safe_trim(rec$hypothesis),
          narrative = safe_trim(rec$narrative),
          key_drivers = safe_trim(rec$key_drivers),
          evidence = safe_trim(rec$evidence),
          limitations = safe_trim(rec$limitations),
          raw_text = safe_trim(rec$raw_text),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}

write_interpretation_markdown <- function(records, path, title, include_raw = FALSE) {
  md <- c(title, "")
  md <- c(md, sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
  if (nrow(records) == 0) {
    md <- c(md, "No interpret_agent records available.")
    writeLines(md, path)
    return(invisible(path))
  }
  for (i in seq_len(nrow(records))) {
    rec <- records[i, , drop = FALSE]
    md <- c(
      md,
      sprintf("## %s | %s | %s", rec$celltype_l2, rec$comparison, rec$direction),
      ""
    )
    md <- c(md, sprintf("- **Status:** %s", rec$status))
    md <- c(md, sprintf("- **Source DB:** %s", ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
    if (nzchar(rec$warnings)) md <- c(md, sprintf("- **Warnings:** %s", rec$warnings))
    if (nzchar(rec$error)) md <- c(md, sprintf("- **Error:** %s", rec$error))
    md <- c(md, "")
    field_block <- list(
      Overview = rec$overview,
      `Key Mechanisms` = rec$key_mechanisms,
      Hypothesis = rec$hypothesis,
      Narrative = rec$narrative,
      `Key Drivers` = rec$key_drivers,
      Evidence = rec$evidence,
      Limitations = rec$limitations
    )
    for (nm in names(field_block)) {
      val <- placeholder_text(field_block[[nm]])
      md <- c(md, sprintf("### %s", nm), "", val, "")
    }
    if (include_raw && nzchar(rec$raw_text)) {
      md <- c(md, "### Raw Output", "", rec$raw_text, "")
    }
  }
  writeLines(md, path)
  invisible(path)
}

interpret_agent_runtime <- new.env(parent = emptyenv())
interpret_agent_runtime$add_ppi <- INTERPRET_AGENT_ADD_PPI

run_interpret_agent_once <- function(enrich_obj, context_str, gene_fc = NULL, add_ppi = interpret_agent_runtime$add_ppi) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(list(result = NULL, warnings = character(), error = "empty enrichment"))
  }
  warn_msgs <- character()
  res <- tryCatch(
    withCallingHandlers(
      clusterProfiler::interpret_agent(
        x = enrich_obj,
        context = context_str,
        n_pathways = INTERPRET_AGENT_N_PATHWAYS,
        model = INTERPRET_AGENT_MODEL,
        api_key = DEEPSEEK_API_KEY,
        add_ppi = add_ppi,
        gene_fold_change = gene_fc
      ),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat(sprintf("    [ERROR] agent: %s\n", e$message))
      structure(list(message = e$message), class = "interpret_agent_error")
    }
  )
  if (inherits(res, "interpret_agent_error")) {
    return(list(result = NULL, warnings = unique(warn_msgs), error = res$message))
  }
  list(result = res, warnings = unique(warn_msgs), error = NULL)
}

run_interpret_agent_safe <- function(enrich_obj, context_str, gene_fc = NULL) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(list(result = NULL, warnings = character(), error = "empty enrichment"))
  }

  combined_warnings <- character()
  last_payload <- list(result = NULL, warnings = character(), error = "interpret_agent did not run")

  for (attempt in seq_len(INTERPRET_AGENT_MAX_RETRIES)) {
    if (attempt > 1) cat(sprintf("    [INFO] interpret_agent retry %d/%d\n", attempt, INTERPRET_AGENT_MAX_RETRIES))
    payload <- run_interpret_agent_once(
      enrich_obj,
      context_str,
      gene_fc,
      add_ppi = interpret_agent_runtime$add_ppi
    )
    last_payload <- payload

    if (
      isTRUE(INTERPRET_AGENT_DISABLE_PPI_ON_STRINGDB_TIMEOUT) &&
      isTRUE(interpret_agent_runtime$add_ppi) &&
      has_stringdb_timeout_issue(payload$warnings, payload$error)
    ) {
      interpret_agent_runtime$add_ppi <- FALSE
      combined_warnings <- c(
        combined_warnings,
        "STRINGdb request timed out; disabling add_ppi for remaining interpret_agent calls"
      )
      cat("    [WARN] STRINGdb request timed out; disabling add_ppi for remaining interpret_agent calls\n")
    }

    if (length(payload$warnings) > 0) {
      combined_warnings <- c(combined_warnings, sprintf("attempt %d: %s", attempt, payload$warnings))
    }

    retryable <- is_retryable_interpret_agent_issue(payload$warnings, payload$error)
    structured <- !is.null(payload$result) && looks_structured_interpret_agent_result(payload$result)
    raw_text_valid <- !is.null(payload$result) && nzchar(capture_object_text(payload$result))

    if (is.null(payload$error) && (structured || (raw_text_valid && !retryable))) {
      payload$warnings <- unique(c(combined_warnings, payload$warnings))
      return(payload)
    }
    if (attempt < INTERPRET_AGENT_MAX_RETRIES && (retryable || is.null(payload$result))) {
      Sys.sleep(INTERPRET_AGENT_RETRY_SLEEP_SEC)
      next
    }
    break
  }

  last_payload$warnings <- unique(c(combined_warnings, last_payload$warnings))
  if (is.null(last_payload$error) && !is.null(last_payload$result)) return(last_payload)
  if (is.null(last_payload$error) || !nzchar(last_payload$error)) {
    last_payload$error <- sprintf("interpret_agent failed after %d attempts", INTERPRET_AGENT_MAX_RETRIES)
  }
  last_payload
}

to_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.character(x) && length(x) == 1) return(x)
  if (is.character(x)) return(paste(x, collapse = "; "))
  if (is.list(x)) return(paste(capture.output(str(x, max.level = 2)), collapse = "\n"))
  as.character(x)
}

stratified_downsample <- function(obj, group_col, n_per = HEATMAP_CELLS_PER_TYPE) {
  set.seed(42)
  meta <- obj@meta.data
  cells <- unlist(lapply(unique(meta[[group_col]]), function(g) {
    gc <- rownames(meta)[meta[[group_col]] == g]
    sample(gc, min(length(gc), n_per))
  }))
  subset(obj, cells = cells)
}

save_plot <- function(p, path_no_ext, width = 10, height = 8) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
}

add_md_image_if_exists <- function(add_fun, output_dir, rel_path, alt_text) {
  full_path <- file.path(output_dir, rel_path)
  if (file.exists(full_path)) {
    add_fun(sprintf("![%s](%s)", alt_text, rel_path))
  } else {
    add_fun(sprintf("_Not generated: `%s`_", rel_path))
  }
  add_fun("")
}

sanitize_obs_for_h5ad <- function(df) {
  out <- df
  for (col in colnames(out)) {
    vec <- out[[col]]
    if (is.factor(vec)) {
      out[[col]] <- as.character(vec)
    } else if (is.logical(vec)) {
      out[[col]] <- ifelse(is.na(vec), NA_integer_, as.integer(vec))
    } else if (inherits(vec, c("POSIXct", "POSIXt"))) {
      out[[col]] <- format(vec, tz = "UTC", usetz = TRUE)
    } else if (inherits(vec, "Date")) {
      out[[col]] <- as.character(vec)
    } else if (is.list(vec)) {
      out[[col]] <- vapply(vec, function(x) {
        if (length(x) == 0 || all(is.na(x))) return(NA_character_)
        paste(as.character(x), collapse = "; ")
      }, character(1))
    } else if (is.character(vec)) {
      out[[col]] <- trimws(vec)
    } else if (!(is.integer(vec) || is.numeric(vec))) {
      out[[col]] <- as.character(vec)
    }
  }
  out
}

matrix_to_scipy_csr <- function(mat, scipy_sparse, np) {
  mat_coo <- as(Matrix::t(as(mat, "dgCMatrix")), "dgTMatrix")
  scipy_sparse$coo_matrix(
    reticulate::tuple(
      np$array(mat_coo@x, dtype = np$float32),
      reticulate::tuple(
        np$array(mat_coo@i, dtype = np$int32),
        np$array(mat_coo@j, dtype = np$int32)
      )
    ),
    shape = reticulate::tuple(as.integer(nrow(mat_coo)), as.integer(ncol(mat_coo)))
  )$tocsr()
}

# ==============================================================================
# 5. Load Data
# ==============================================================================

cat("=== Loading Epithelial Data ===\n")

if (!file.exists(H5AD_PATH)) stop(sprintf("File not found: %s", H5AD_PATH))
obj <- GetSeurat(
  h5ad_path = H5AD_PATH,
  prefer_raw = FALSE,
  prefer_layer_counts = TRUE,
  validate_counts = TRUE,
  debug = TRUE
)
cat(sprintf("[OK] %d cells x %d genes\n", ncol(obj), nrow(obj)))

if (!"counts" %in% Layers(obj[["RNA"]])) {
  stop("RNA assay missing 'counts' layer. Check GetSeurat output.")
}
cat("[OK] counts layer verified\n\n")

fig_dir <- file.path(OUTPUT_DIR, "figures")
rpt_dir <- file.path(OUTPUT_DIR, "reports")
de_dir <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir <- file.path(OUTPUT_DIR, "wilcox_exploratory")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rpt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wx_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. Validate Metadata & Apply L2 Remapping
# ==============================================================================

cat("=== Validating Metadata & L2 Remapping ===\n")

meta <- obj@meta.data
for (col in c(TISSUE_COL, SAMPLE_COL, L3_SOURCE_COL)) {
  if (!col %in% colnames(meta)) stop(sprintf("Missing required column: %s", col))
}

l3_vals <- as.character(meta[[L3_SOURCE_COL]])
l2_vals <- L3_TO_L2_REMAP[l3_vals]

unmapped <- unique(l3_vals[is.na(l2_vals)])
if (length(unmapped) > 0) {
  stop(sprintf(
    "Unmapped L3 values found in '%s': %s\nUpdate L3_TO_L2_REMAP before running.",
    L3_SOURCE_COL,
    paste(unmapped, collapse = ", ")
  ))
}

obj@meta.data[[CELLTYPE_L2_COL]] <- l2_vals
obj@meta.data[[CELLTYPE_L3_COL]] <- l3_vals

cat("\nL3 -> L2 Remapping Summary:\n")
print(table(L3 = l3_vals, L2 = l2_vals, useNA = "ifany"))
cat("\nL2 Distribution:\n")
print(table(obj@meta.data[[CELLTYPE_L2_COL]], useNA = "ifany"))
cat("\n")

bad_idx <- is.na(meta[[TISSUE_COL]]) |
  trimws(as.character(meta[[TISSUE_COL]])) == "" |
  is.na(meta[[SAMPLE_COL]]) |
  trimws(as.character(meta[[SAMPLE_COL]])) == "" |
  is.na(obj@meta.data[[CELLTYPE_L2_COL]]) |
  is.na(meta[[LABEL_COL]]) |
  trimws(as.character(meta[[LABEL_COL]])) == ""
if (sum(bad_idx) > 0) {
  cat(sprintf("[INFO] Dropping %d cells with NA tissue/sample/L2/label\n", sum(bad_idx)))
  obj <- subset(obj, cells = colnames(obj)[!bad_idx])
}

n_samples <- length(unique(obj@meta.data[[SAMPLE_COL]]))
if (n_samples < 2) stop(sprintf("Only %d unique sample(s). Need >= 2 for pseudobulk DE.", n_samples))
cat(sprintf("[OK] %d unique samples\n", n_samples))

meta <- obj@meta.data
tissues <- sort(unique(na.omit(meta[[TISSUE_COL]])))
l2_types <- sort(unique(na.omit(meta[[CELLTYPE_L2_COL]])))

cat(sprintf("[OK] Cells: %d\n", ncol(obj)))
cat(sprintf("[OK] Tissues: %s\n", paste(tissues, collapse = ", ")))
cat(sprintf("[OK] L2 types (%d): %s\n", length(l2_types), paste(l2_types, collapse = ", ")))
cat(sprintf("[OK] L3 column (label): %s (%d unique)\n", LABEL_COL, length(unique(meta[[LABEL_COL]]))))
cat("\n")

cat("[INFO] Running NormalizeData() for visualization and marker analyses\n")
obj <- NormalizeData(obj, verbose = FALSE)
cat("[OK] data layer ready\n\n")

# ==============================================================================
# 7. Visualization
# ==============================================================================

cat("=== Visualization ===\n")

umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
if (!is.null(umap_reduction)) {
  tissue_cols_use <- UMAP_TISSUE_COLORS[names(UMAP_TISSUE_COLORS) %in% tissues]

  p1 <- build_umap_plot(
    obj,
    umap_reduction,
    TISSUE_COL,
    title = sprintf("Epithelial - Tissue (%s)", umap_reduction),
    cols = tissue_cols_use,
    width = 11,
    height = 8
  )
  save_plot(p1$plot, file.path(fig_dir, "umap_tissue"), p1$width, p1$height)

  p2 <- build_umap_plot(
    obj,
    umap_reduction,
    CELLTYPE_L2_COL,
    title = sprintf("Epithelial - Cell Type L2 (%s)", umap_reduction),
    label = TRUE,
    width = 14,
    height = 10
  )
  save_plot(p2$plot, file.path(fig_dir, "umap_celltype_L2"), p2$width, p2$height)

  p2s <- build_umap_plot(
    obj,
    umap_reduction,
    CELLTYPE_L2_COL,
    title = sprintf("Epithelial - L2 by Tissue (%s)", umap_reduction),
    split_col = TISSUE_COL,
    label = TRUE,
    width = max(12, 4.5 * length(tissues)),
    height = 8
  )
  save_plot(p2s$plot, file.path(fig_dir, "umap_L2_split_tissue"), p2s$width, p2s$height)

  p2_l3 <- build_umap_plot(
    obj,
    umap_reduction,
    LABEL_COL,
    title = sprintf("Epithelial - Cell Type L3 (%s)", umap_reduction),
    label = TRUE,
    width = 18,
    height = 12
  )
  save_plot(p2_l3$plot, file.path(fig_dir, "umap_celltype_L3"), p2_l3$width, p2_l3$height)

  cat(sprintf("[OK] UMAP saved using reduction: %s\n", umap_reduction))
} else {
  cat("[WARN] No UMAP reduction found\n")
}

markers_present <- intersect(KNOWN_MARKERS, rownames(obj))
if (length(markers_present) >= 3) {
  Idents(obj) <- LABEL_COL
  p3 <- DotPlot(obj, features = markers_present) +
    RotatedAxis() +
    ggtitle("Epithelial - Known Markers (L3)") +
    theme(axis.text.x = element_text(size = 7))
  save_plot(
    p3,
    file.path(fig_dir, "dotplot_markers"),
    width = max(14, length(markers_present) * 0.35),
    height = max(7, length(unique(meta[[LABEL_COL]])) * 0.4)
  )
  cat("[OK] Dotplot saved\n")
}

comp_df <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL))) %>%
  count(!!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(TISSUE_COL)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p4 <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(CELLTYPE_L2_COL))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Tissue", y = "Percentage (%)", title = "Epithelial - L2 Composition by Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4, file.path(fig_dir, "composition_tissue_L2"))

p4b <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = n, fill = !!sym(CELLTYPE_L2_COL))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Tissue", y = "Cell Count", title = "Epithelial - Absolute Count by Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4b, file.path(fig_dir, "count_tissue_L2"))
cat("[OK] Composition plots saved\n")

sample_comp <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL)), !is.na(!!sym(SAMPLE_COL))) %>%
  count(!!sym(SAMPLE_COL), !!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(SAMPLE_COL)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

if (nrow(sample_comp) > 0) {
  facet_ncol <- min(3, length(l2_types))
  facet_nrow <- ceiling(length(l2_types) / facet_ncol)
  p4c <- ggplot(sample_comp, aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(TISSUE_COL))) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(as.formula(paste("~", CELLTYPE_L2_COL)), scales = "free_y", ncol = facet_ncol) +
    labs(x = "Tissue", y = "Proportion per Sample (%)", title = "Epithelial - Sample-level Composition by Tissue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  save_plot(
    p4c,
    file.path(fig_dir, "composition_sample_level"),
    width = max(10, facet_ncol * 4),
    height = max(6, facet_nrow * 3.5)
  )
  cat("[OK] Sample-level composition saved\n")
}

Idents(obj) <- CELLTYPE_L2_COL
top_mk <- tryCatch(
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.2,
    logfc.threshold = 0.25,
    max.cells.per.ident = 1000,
    test.use = "wilcox"
  ),
  error = function(e) {
    cat(sprintf("[WARN] FindAllMarkers failed: %s\n", e$message))
    NULL
  }
)

if (!is.null(top_mk) && nrow(top_mk) > 0) {
  fwrite(top_mk, file.path(rpt_dir, "all_markers_per_L2.csv"))
  top10 <- top_mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
  obj_ds <- stratified_downsample(obj, CELLTYPE_L2_COL, HEATMAP_CELLS_PER_TYPE)
  obj_ds <- ScaleData(obj_ds, features = unique(top10$gene), verbose = FALSE)
  p5 <- DoHeatmap(obj_ds, features = unique(top10$gene), size = 3) +
    ggtitle("Epithelial - Top Markers per L2 (downsampled, visualization only)")
  save_plot(p5, file.path(fig_dir, "heatmap_top_markers"), width = 14, height = 10)
  rm(obj_ds)
  gc()
  cat("[OK] Heatmap saved\n")
}
cat("\n")

# ==============================================================================
# 8. Pseudobulk DESeq2 + Enrichment + interpret_agent
# ==============================================================================

cat("=== Pseudobulk DESeq2 Tissue Comparison ===\n")
cat(sprintf("Thresholds: padj < %s, |log2FC| > %s\n", PADJ_THR, LFC_THR))
cat("Enrichment databases: GO_BP, GO_MF, GO_CC, KEGG, Hallmark, CellMarker, PanglaoDB, Epithelial_custom\n\n")

pb_de_all <- list()
enrich_all <- list()
agent_all <- list()
agent_structured_all <- list()

for (ct_l2 in l2_types) {
  cat(sprintf("\n>> L2: %s\n", ct_l2))
  pb <- aggregate_pseudobulk(obj, ct_l2)
  if (is.null(pb)) {
    cat("  [SKIP] Insufficient pseudobulk\n")
    next
  }

  pb_tissues <- unique(na.omit(pb$meta$tissue))
  if (length(pb_tissues) < 2) {
    cat("  [SKIP] < 2 tissues\n")
    next
  }

  cat(sprintf("  Pseudobulk samples: %d (%s)\n", nrow(pb$meta), paste(pb_tissues, collapse = ", ")))
  cat(sprintf("  Grouping columns: %s\n", paste(pb$group_cols, collapse = ", ")))
  pb_de_all[[ct_l2]] <- list()
  enrich_all[[ct_l2]] <- list()
  agent_all[[ct_l2]] <- list()
  agent_structured_all[[ct_l2]] <- list()

  for (pair in combn(as.character(pb_tissues), 2, simplify = FALSE)) {
    t1 <- pair[1]
    t2 <- pair[2]
    comp_name <- paste0(t2, "_vs_", t1)
    cat(sprintf("  DESeq2: %s\n", comp_name))

    res <- run_deseq2_pairwise(pb, t1, t2)
    if (is.null(res)) {
      cat("    [SKIP] Too few samples per tissue\n")
      next
    }

    pb_de_all[[ct_l2]][[comp_name]] <- res
    cat(sprintf(
      "    DEGs: %d up, %d down (samples: %d vs %d)\n",
      res$n_up, res$n_down, res$n_samples_1, res$n_samples_2
    ))

    comp_dir <- file.path(de_dir, safe_name(ct_l2), safe_name(comp_name))
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(res$de_table, file.path(comp_dir, "DESeq2_results.csv"))

    vp <- res$de_table %>% mutate(label = ifelse(sig == "sig" & rank(padj) <= 20, gene, ""))
    pv <- ggplot(vp, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("sig" = "red", "ns" = "grey70")) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
      geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "blue") +
      labs(
        title = sprintf("DESeq2: %s (%s)", comp_name, ct_l2),
        x = "log2 Fold Change",
        y = "-log10(padj)"
      ) +
      theme_minimal()
    save_plot(pv, file.path(comp_dir, "volcano"), width = 8, height = 6)

    enrich_all[[ct_l2]][[comp_name]] <- list()
    agent_all[[ct_l2]][[comp_name]] <- list()
    agent_structured_all[[ct_l2]][[comp_name]] <- list()

    for (dir_name in c("up", "down")) {
      genes <- if (dir_name == "up") {
        res$de_table %>% filter(sig == "sig", log2FoldChange > 0) %>%
          arrange(desc(log2FoldChange)) %>% head(TOP_N_DEG_ENRICHMENT) %>% pull(gene)
      } else {
        res$de_table %>% filter(sig == "sig", log2FoldChange < 0) %>%
          arrange(log2FoldChange) %>% head(TOP_N_DEG_ENRICHMENT) %>% pull(gene)
      }

      if (length(genes) < 5) {
        enrich_all[[ct_l2]][[comp_name]][[dir_name]] <- NULL
        next
      }
      cat(sprintf("    Enrichment [%s]: %d genes\n", dir_name, length(genes)))

      tested <- res$tested_genes
      enr_list <- list(
        GO_BP = run_gmt_enrichment(genes, go_bp_t2g, "GO_BP", tested),
        GO_MF = run_gmt_enrichment(genes, go_mf_t2g, "GO_MF", tested),
        GO_CC = run_gmt_enrichment(genes, go_cc_t2g, "GO_CC", tested),
        KEGG = run_gmt_enrichment(genes, kegg_t2g, "KEGG", tested),
        Hallmark = run_gmt_enrichment(genes, hallmark_t2g, "Hallmark", tested),
        CellMarker = run_gmt_enrichment(genes, cellmarker_t2g, "CellMarker", tested),
        PanglaoDB = run_gmt_enrichment(genes, panglaodb_t2g, "PanglaoDB", tested),
        Epithelial_custom = run_gmt_enrichment(genes, epithelial_custom_t2g, "Epithelial_custom", tested)
      )
      enr_list <- enr_list[!sapply(enr_list, is.null)]
      enrich_all[[ct_l2]][[comp_name]][[dir_name]] <- enr_list

      enr_out <- file.path(comp_dir, paste0("enrichment_", dir_name))
      dir.create(enr_out, showWarnings = FALSE)
      for (db in names(enr_list)) {
        er <- enr_list[[db]]
        if (nrow(as.data.frame(er)) > 0) {
          fwrite(as.data.frame(er), file.path(enr_out, paste0(db, ".csv")))
          tryCatch(
            {
              pdf(file.path(enr_out, paste0(db, "_dotplot.pdf")), width = 10, height = 8)
              print(dotplot(er, showCategory = 15, title = sprintf("%s %s (%s %s)", db, dir_name, ct_l2, comp_name)))
              dev.off()
            },
            error = function(e) NULL
          )
        }
      }
      saveRDS(enr_list, file.path(enr_out, "all_enrichment.rds"))

      if (ENABLE_LLM) {
        tissue_ctx <- build_tissue_pair_context(t1, t2, ct_l2, dir_name)
        sig_de <- res$de_table %>% filter(sig == "sig")
        gene_fc <- setNames(sig_de$log2FoldChange, toupper(sig_de$gene))
        selected_dbs <- select_interpret_agent_databases(enr_list)

        if (length(selected_dbs) > 0) {
          cat(sprintf(
            "    interpret_agent [%s] integrating %d database(s): %s\n",
            dir_name,
            length(selected_dbs),
            paste(selected_dbs, collapse = ", ")
          ))

          multi_payload <- run_interpret_agent_multi(enr_list[selected_dbs], tissue_ctx, dir_name, gene_fc)
          agent_all[[ct_l2]][[comp_name]][[dir_name]] <- multi_payload

          ia_structured <- build_multi_database_interpretation_record(
            multi_payload = multi_payload,
            enrich_list = enr_list,
            ct_l2 = ct_l2,
            comp_name = comp_name,
            dir_name = dir_name,
            context_str = tissue_ctx,
            gene_fc = gene_fc
          )
          agent_structured_all[[ct_l2]][[comp_name]][[dir_name]] <- ia_structured

          saveRDS(multi_payload, file.path(comp_dir, paste0("interpret_agent_", dir_name, "_payload.rds")))
          saveRDS(ia_structured, file.path(comp_dir, paste0("interpret_agent_", dir_name, "_structured.rds")))
          if (length(multi_payload$per_db_payloads) > 0) {
            for (db in names(multi_payload$per_db_payloads)) {
              saveRDS(
                multi_payload$per_db_payloads[[db]],
                file.path(comp_dir, paste0("interpret_agent_", dir_name, "_", safe_name(db), "_payload.rds"))
              )
            }
          }
          if (nzchar(ia_structured$raw_text)) {
            writeLines(ia_structured$raw_text, file.path(comp_dir, paste0("interpret_agent_", dir_name, "_raw.txt")))
          }
        } else {
          cat(sprintf("    [SKIP] interpret_agent [%s]: no enrichment database met minimum term count\n", dir_name))
        }
      }
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
  wx <- run_wilcox_exploratory(obj, ct_l2)
  if (!is.null(wx) && length(wx) > 0) {
    wilcox_all[[ct_l2]] <- wx
    for (comp_name in names(wx)) {
      wx_out <- file.path(wx_dir, safe_name(ct_l2))
      dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(wx[[comp_name]], file.path(wx_out, paste0(safe_name(comp_name), "_wilcox.csv")))
    }
    cat(sprintf("  [OK] %s: %d comparisons\n", ct_l2, length(wx)))
  }
}

# ==============================================================================
# 10. Save Summary RDS
# ==============================================================================

cat("\n=== Saving Summary ===\n")
saveRDS(pb_de_all, file.path(rpt_dir, "pseudobulk_de_all.rds"))
saveRDS(wilcox_all, file.path(rpt_dir, "wilcox_exploratory_all.rds"))
saveRDS(enrich_all, file.path(rpt_dir, "enrichment_all.rds"))
saveRDS(agent_all, file.path(rpt_dir, "interpret_agent_all.rds"))
saveRDS(agent_structured_all, file.path(rpt_dir, "interpret_agent_structured_all.rds"))
saveRDS(L3_TO_L2_REMAP, file.path(rpt_dir, "L3_to_L2_remap.rds"))

agent_structured_df <- flatten_interpretation_records(agent_structured_all)
if (nrow(agent_structured_df) > 0) {
  fwrite(agent_structured_df, file.path(rpt_dir, "interpret_agent_structured.tsv"), sep = "\t")
}

write_interpretation_markdown(
  agent_structured_df,
  file.path(OUTPUT_DIR, "LLM_INTERPRETATION.md"),
  "# LLM Interpretation Summary (Epithelial Tissue Comparison)",
  include_raw = FALSE
)
write_interpretation_markdown(
  agent_structured_df,
  file.path(rpt_dir, "LLM_INTERPRETATION_FOR_LLM.md"),
  "# LLM Interpretation Structured Input",
  include_raw = TRUE
)

# ==============================================================================
# 11. Generate REPORT.md
# ==============================================================================

cat("\n=== Generating REPORT.md ===\n")

md <- character()
add <- function(...) md <<- c(md, paste0(...))

add("# Epithelial Tissue Comparison Report")
add("")
add("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M"))
add("")
add("**Pipeline:** Epithelial Tissue Comparison v1.0.0 (pseudobulk DESeq2, multi-database enrichment)")
add("")
add("**Note:** This is intended as a cross-site anatomical comparison of respiratory epithelial cells.")
add("**Interpretation caution:** Unless disease covariates are explicitly modeled, do not over-interpret differences as disease effects.")
add("")
add("---")
add("")

add("## 1. Data Overview")
add("")
add(sprintf("- **Input:** `%s`", basename(H5AD_PATH)))
add(sprintf("- **Total cells:** %s", format(ncol(obj), big.mark = ",")))
add(sprintf("- **Tissues:** %s", paste(tissues, collapse = ", ")))
add(sprintf("- **L2 lineages (%d):** %s", length(l2_types), paste(l2_types, collapse = ", ")))
add(sprintf("- **L3 source column:** `%s` (%d unique)", L3_SOURCE_COL, length(unique(meta[[LABEL_COL]]))))
add("")
add("### L3 -> L2 Remapping")
add("")
add("| L3 | L2 |")
add("|---|---|")
for (i in seq_along(L3_TO_L2_REMAP)) {
  add(sprintf("| %s | %s |", names(L3_TO_L2_REMAP)[i], L3_TO_L2_REMAP[i]))
}
add("")

add("## 2. Visualization")
add("")
add("### 2.1 UMAP")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_tissue.png"), "UMAP tissue")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_celltype_L2.png"), "UMAP L2")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_celltype_L3.png"), "UMAP L3")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "umap_L2_split_tissue.png"), "UMAP split")
add("")
add("### 2.2 Marker Dotplot")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "dotplot_markers.png"), "Dotplot")
add("")
add("### 2.3 Cell Composition")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "composition_tissue_L2.png"), "Composition")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "count_tissue_L2.png"), "Counts")
add("")
add("### 2.4 Sample-level Composition")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "composition_sample_level.png"), "Sample composition")
add("")
add("### 2.5 Top Marker Heatmap")
add_md_image_if_exists(add, OUTPUT_DIR, file.path("figures", "heatmap_top_markers.png"), "Heatmap")
add("")

add("## 3. Pseudobulk DESeq2 (Primary Inference)")
add("")
add(sprintf(
  "Statistical unit: pseudobulk (sum counts per sample x tissue x L2). Thresholds: padj < %s, |log2FC| > %s.",
  PADJ_THR,
  LFC_THR
))
add("")
add("| L2 Lineage | Comparison | Samples (ref/case) | Design | Up | Down | Total |")
add("|---|---|---|---|---|---|---|")
for (ct_l2 in names(pb_de_all)) {
  for (comp_name in names(pb_de_all[[ct_l2]])) {
    r <- pb_de_all[[ct_l2]][[comp_name]]
    if (is.null(r)) next
    add(sprintf(
      "| %s | %s | %d / %d | `%s` | %d | %d | %d |",
      ct_l2,
      comp_name,
      r$n_samples_1,
      r$n_samples_2,
      r$design,
      r$n_up,
      r$n_down,
      r$n_up + r$n_down
    ))
  }
}
add("")

add("## 4. Multi-Database Enrichment")
add("")
add("Databases: GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, custom epithelial markers")
add("")
for (ct_l2 in names(enrich_all)) {
  for (comp_name in names(enrich_all[[ct_l2]])) {
    for (dir_name in names(enrich_all[[ct_l2]][[comp_name]])) {
      enr_l <- enrich_all[[ct_l2]][[comp_name]][[dir_name]]
      if (is.null(enr_l) || length(enr_l) == 0) next
      add(sprintf("### %s | %s | %s", ct_l2, comp_name, dir_name))
      add("")
      for (db in names(enr_l)) {
        er <- enr_l[[db]]
        if (is.null(er) || nrow(as.data.frame(er)) == 0) next
        top5 <- head(as.data.frame(er), 5)
        add(sprintf("**%s (top 5):**", db))
        add("")
        add("| Term | p.adjust | Count |")
        add("|---|---|---|")
        for (j in seq_len(nrow(top5))) {
          add(sprintf("| %s | %.2e | %s |", top5$Description[j], top5$p.adjust[j], top5$Count[j]))
        }
        add("")
      }
    }
  }
}

add("## 5. LLM Interpretation (interpret_agent)")
add("")
if (!ENABLE_LLM) {
  add("**SKIPPED:** DEEPSEEK_API_KEY not set. Re-run with API key to enable.")
  add("")
} else if (nrow(agent_structured_df) == 0) {
  add("No interpret_agent results generated (insufficient enrichment or all calls failed).")
  add("")
} else {
  add("Structured outputs:")
  add(sprintf("- `%s`", file.path("reports", "interpret_agent_structured.tsv")))
  add(sprintf("- `%s`", "LLM_INTERPRETATION.md"))
  add("- LLM output language: Chinese (Simplified).")
  add("- LLM emphasis: top-|log2FC| genes linked to site-relevant epithelial pathways.")
  add("- LLM evidence integration: all eligible enrichment databases are interpreted and merged into a single consensus summary.")
  add("")
  for (i in seq_len(nrow(agent_structured_df))) {
    rec <- agent_structured_df[i, , drop = FALSE]
    add(sprintf("### %s | %s | %s", rec$celltype_l2, rec$comparison, rec$direction))
    add("")
    add(sprintf("- **Status:** %s", rec$status))
    add(sprintf("- **Source DB(s):** %s", ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
    if (nzchar(rec$overview)) add(sprintf("- **Overview:** %s", rec$overview))
    if (nzchar(rec$key_mechanisms)) add(sprintf("- **Key Mechanisms:** %s", rec$key_mechanisms))
    if (nzchar(rec$hypothesis)) add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
    add("")
  }
}

add("## 6. Methods")
add("")
add("- **DE:** Pseudobulk DESeq2 (sample-level aggregation)")
add("- **Exploratory:** Cell-level Wilcoxon rank-sum (marker discovery only, NOT inference; only top markers per direction are exported)")
add("- **Enrichment:** clusterProfiler::enricher() + MSigDB GMT + CellMarker + PanglaoDB + custom epithelial markers")
add("- **Enrichment universe:** DESeq2-tested genes intersected with TERM2GENE")
add("- **Enrichment databases (8):** GO BP, GO MF, GO CC, KEGG, Hallmark, CellMarker, PanglaoDB, epithelial custom")
if (ENABLE_LLM) {
  add("- **LLM:** interpret_agent with DeepSeek deepseek-reasoner (tissue-pair-specific epithelial context)")
  add("- **LLM integration:** all eligible enrichment databases are interpreted separately, then merged into one consensus summary")
  add("- **LLM post-processing:** secondary deepseek-chat pass standardizes outputs to fixed JSON schema when needed")
  add("- **LLM output language:** Chinese (Simplified); gene symbols retained in English")
} else {
  add("- **LLM:** SKIPPED (no API key)")
}
if (!is.null(DESIGN_COVARIATE)) {
  add(sprintf("- **Requested DESeq2 covariate:** `%s` (propagated into pseudobulk metadata when available)", DESIGN_COVARIATE))
} else {
  add("- **Requested DESeq2 covariate:** none")
}
add("- **DESeq2 design:** actual formula used for each comparison is recorded in Section 3 and in each DESeq2 result object")
add("- **L2 remapping:** L3 epithelial states merged into 6 major lineages")
add(sprintf("- **Context:** Respiratory epithelial site comparison across %s", paste(tissues, collapse = ", ")))
add("")
add("### Output Objects")
add("")
add("- `epithelial_tissue_comparison_final.rds` -- Seurat object with `cell_type_L2` + `cell_type_L3`")
add("- `epithelial_tissue_comparison_final.h5ad` -- AnnData object with `cell_type_L2` + `cell_type_L3`")
add("")
add("---")
add("*Generated by epithelial_tissue_comparison_v1_0_20260330.R*")

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

cat(sprintf(
  "[OK] cell_type_L2: %d unique (%s)\n",
  length(unique(obj@meta.data[["cell_type_L2"]])),
  paste(sort(unique(obj@meta.data[["cell_type_L2"]])), collapse = ", ")
))
cat(sprintf(
  "[OK] cell_type_L3: %d unique (%s)\n",
  length(unique(obj@meta.data[["cell_type_L3"]])),
  paste(sort(unique(obj@meta.data[["cell_type_L3"]])), collapse = ", ")
))

rds_path <- file.path(OUTPUT_DIR, "epithelial_tissue_comparison_final.rds")
saveRDS(obj, rds_path)
cat(sprintf("[OK] RDS saved: %s (%.1f MB)\n", basename(rds_path), file.size(rds_path) / 1e6))

h5ad_out_path <- file.path(OUTPUT_DIR, "epithelial_tissue_comparison_final.h5ad")

tryCatch(
  {
    anndata <- reticulate::import("anndata")
    scipy_sparse <- reticulate::import("scipy.sparse")
    np <- reticulate::import("numpy")

    counts_mat <- GetAssayData(obj, layer = "counts")
    norm_layer_name <- if ("data" %in% Layers(obj[["RNA"]])) "data" else "counts"
    expr_mat <- GetAssayData(obj, layer = norm_layer_name)

    counts_scipy <- matrix_to_scipy_csr(counts_mat, scipy_sparse, np)
    expr_scipy <- if (identical(norm_layer_name, "counts")) counts_scipy else matrix_to_scipy_csr(expr_mat, scipy_sparse, np)

    meta_export <- sanitize_obs_for_h5ad(obj@meta.data)
    obs_df <- reticulate::r_to_py(meta_export)
    var_df <- data.frame(gene_symbol = rownames(obj), row.names = rownames(obj), stringsAsFactors = FALSE)
    var_py <- reticulate::r_to_py(var_df)

    adata <- anndata$AnnData(X = expr_scipy, obs = obs_df, var = var_py)
    adata$layers[["counts"]] <- counts_scipy
    adata$uns[["X_layer"]] <- norm_layer_name
    for (red_name in Reductions(obj)) {
      emb <- Embeddings(obj, reduction = red_name)
      adata$obsm[[paste0("X_", red_name)]] <- np$array(emb, dtype = np$float32)
    }

    adata$write_h5ad(h5ad_out_path, compression = "gzip")
    cat(sprintf("[OK] h5ad saved: %s (%.1f MB)\n", basename(h5ad_out_path), file.size(h5ad_out_path) / 1e6))

    adata_check <- anndata$read_h5ad(h5ad_out_path)
    stopifnot("cell_type_L2" %in% reticulate::py_to_r(adata_check$obs$columns$tolist()))
    stopifnot("cell_type_L3" %in% reticulate::py_to_r(adata_check$obs$columns$tolist()))
    stopifnot("counts" %in% reticulate::py_to_r(adata_check$layers$keys()))
    cat(sprintf(
      "[OK] h5ad verified: %d cells x %d genes, cell_type_L2 + cell_type_L3 present, layers['counts'] present\n",
      reticulate::py_to_r(adata_check$n_obs),
      reticulate::py_to_r(adata_check$n_vars)
    ))
    rm(adata, adata_check, counts_mat, expr_mat, counts_scipy, expr_scipy, meta_export, var_df)
    gc()
  },
  error = function(e) {
    cat(sprintf("[ERROR] h5ad export failed: %s\n", e$message))
    cat("[INFO] RDS saved successfully; convert manually if h5ad needed.\n")
  }
)

# ==============================================================================
# 13. Final Summary
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("EPITHELIAL TISSUE COMPARISON COMPLETE (v1.0.0)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Directory structure:\n")
cat("  REPORT.md                               <- structured report (png embeds)\n")
cat("  LLM_INTERPRETATION.md                   <- structured LLM interpretation summary\n")
cat("  epithelial_tissue_comparison_final.rds  <- Seurat object (cell_type_L2 + cell_type_L3)\n")
cat("  epithelial_tissue_comparison_final.h5ad <- AnnData object (cell_type_L2 + cell_type_L3)\n")
cat("  figures/                                <- UMAP, dotplot, heatmap, composition (pdf+png)\n")
cat("  pseudobulk_de/                          <- DESeq2 per L2 x tissue pair + volcano + enrichment\n")
cat("  wilcox_exploratory/                     <- cell-level wilcox (marker discovery only)\n")
cat("  reports/                                <- RDS summary objects + structured LLM files\n")
cat("\n")

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
cat("[OK] Done\n")
