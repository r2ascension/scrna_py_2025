#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Tissue Comparison Pipeline v1.2.0
# ==============================================================================
#
# Purpose:
#   Adapt the current validated tissue-comparison workflow used for B cells
#   (shared engine synced to 2026-04-06 template) to the myeloid reference
#   object, including DESeq2/ssGSEA/CHOIR and all LLM layers.
#
# Data structure inferred from the myeloid reference / patch workflow:
#   - input h5ad: /home/h2048/data/py/0329/adata_myeloid_L3refined_patched_v1.h5ad
#   - key metadata: sample, tissue, cell_type_L3_refined, ann_finest_level
#   - tissues: nose, sinus, respiratory airway, lung parenchyma
#   - refined L3 states: 13 myeloid states
#   - legacy cell_type_L2 exists but is inconsistent with refined L3 for part of
#     the dataset (e.g. neutrophils nested under classical monocytes; cDC2 split
#     across DC/DC2; some interstitial macrophages under intestinal macrophages)
#
# Strategy:
#   - use cell_type_L3_refined as the standardized L3 source
#   - rebuild a cleaner lineage-level L2 from refined L3 labels
#   - reuse the shared tissue-comparison engine from the latest B-cell script by
#     overriding lineage-specific configuration only
#
# Author: GitHub Copilot
# Date:   2026-04-07
# ==============================================================================

# ----- Lineage identity -----
LINEAGE_TAG           <- "MYELOID"
LINEAGE_DISPLAY       <- "Myeloid"
LINEAGE_CONTEXT_LABEL <- "myeloid cells"
LINEAGE_CONTEXT_LOWER <- "myeloid"
BIOLOGICAL_QUESTION_FRAGMENT <- paste(
  "innate immune surveillance, monocyte/macrophage tissue adaptation,",
  "antigen presentation, interferon tone, inflammatory recruitment,",
  "or tissue-resident phagocyte specialization"
)

# ----- Versioning / outputs -----
PIPELINE_VERSION_LABEL <- "v1.2.0"
PIPELINE_SUBTITLE <- sprintf(
  "%s Tissue Comparison %s (pseudobulk DESeq2, multi-database enrichment, LLM-extended)",
  LINEAGE_DISPLAY, PIPELINE_VERSION_LABEL
)
REPORT_TITLE <- "# Myeloid Tissue Comparison Report (Normal Respiratory Tract)"
GENERATED_BY_LABEL <- "myeloid_tissue_comparison_20260407_v1_2.R"
FINAL_FILE_PREFIX  <- "myeloid_tissue_comparison_final"
LINEAGE_COMPLETION_BANNER <- sprintf(
  "%s TISSUE COMPARISON COMPLETE (%s)",
  toupper(LINEAGE_DISPLAY),
  PIPELINE_VERSION_LABEL
)

# ----- Input / output -----
H5AD_PATH  <- "/home/h2048/data/py/0329/adata_myeloid_L3refined_patched_v1.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0407/myeloid_tissue_comparison_v1_2"

# ----- Interpretation / CHOIR overrides -----
INTERPRET_MULTI_DB_MAX_DBS       <- 5L
INTERPRET_MULTI_DB_TERMS_PER_DB  <- 4L
CHOIR_ALPHA                      <- 0.20
CHOIR_VAR_FEATURES_MAX           <- 4000L

# ----- Standardized label construction -----
# We intentionally do NOT trust the legacy cell_type_L2 column for downstream
# comparison because it is partly inconsistent with cell_type_L3_refined.
USE_EXISTING_L2 <- FALSE
L2_SOURCE_COL   <- "cell_type_L2"
L3_SOURCE_COL   <- "cell_type_L3_refined"
ANALYSIS_L3_DESCRIPTION <- "fine refined myeloid states from `cell_type_L3_refined`"
L3_TO_L2_TABLE_HEADER_LEFT <- "L3 refined (`cell_type_L3_refined`)"

L3_TO_L2_REMAP <- c(
  "Resident Alveolar macrophages"         = "Alveolar_Macrophage",
  "Resting Alveolar macrophages"          = "Alveolar_Macrophage",
  "Interstitial macrophages"              = "Interstitial_Macrophage",
  "Inflammatory Interstitial macrophages" = "Interstitial_Macrophage",
  "CD163L1+ Interstitial macrophages"     = "Interstitial_Macrophage",
  "Typical Classical monocytes"           = "Monocyte",
  "Inflammatory Classical monocytes"      = "Monocyte",
  "Non-classical monocytes"               = "Monocyte",
  "Conventional cDC2"                     = "DC",
  "Langerhans-like cDC2"                  = "DC",
  "pDC"                                   = "DC",
  "Mast cells"                            = "Mast_cell",
  "Neutrophils"                           = "Neutrophil"
)

# ----- Reference DB filters -----
CUSTOM_DB_NAME  <- "Myeloid_custom"
CUSTOM_DB_LABEL <- "Myeloid custom markers"
CUSTOM_ENRICHMENT_SIZE_RULE <- c(min = 2L, max = 120L)
CELLMARKER_CELLTYPE_PATTERN <- paste(
  "Myeloid|Monocyte|Macrophage|Alveolar macrophage|Interstitial macrophage|",
  "Dendritic|DC|cDC|pDC|Plasmacytoid DC|Mast cell|Neutrophil|Langerhans"
)
PANGLAODB_CELLTYPE_PATTERN <- paste(
  "Myeloid|Monocyte|Macrophage|Alveolar macrophage|Interstitial macrophage|",
  "Dendritic|DC|cDC|pDC|Mast cell|Neutrophil|Langerhans"
)

# ----- Visualization markers -----
KNOWN_MARKERS <- c(
  "LYZ", "LST1", "TYMP", "CTSS", "FCER1G", "HLA-DRA",
  "FCN1", "VCAN", "S100A8", "S100A9", "CTSD", "SAT1", "FCGR3A", "IFITM3",
  "FABP4", "PPARG", "MARCO", "C1QA", "C1QB", "C1QC", "INHBA",
  "APOE", "CD163", "CD163L1", "MRC1", "FOLR2",
  "FCER1A", "CD1C", "CLEC10A", "CLEC4C", "LILRA4", "GZMB", "IRF7",
  "KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2",
  "FCGR3B", "CXCL8", "CSF3R", "NAMPT"
)

# ----- Custom lineage marker DB -----
CUSTOM_MARKERS_DB <- data.frame(
  subtype = c(
    "Resident Alveolar macrophages",
    "Resting Alveolar macrophages",
    "Interstitial macrophages",
    "Inflammatory Interstitial macrophages",
    "CD163L1+ Interstitial macrophages",
    "Typical Classical monocytes",
    "Inflammatory Classical monocytes",
    "Non-classical monocytes",
    "Conventional cDC2",
    "Langerhans-like cDC2",
    "pDC",
    "Mast cells",
    "Neutrophils"
  ),
  markers = c(
    "FABP4,PPARG,MARCO,C1QA,C1QB,C1QC,INHBA",
    "FABP4,PPARG,MARCO,C1QA,C1QB,EAR1,ABCA1",
    "APOE,C1QC,CTSB,CTSD,LYZ,FCER1G,TYMP",
    "IL1B,CXCL8,CCL3,NFKBIA,SOD2,CTSB,TYMP",
    "CD163L1,CD163,MRC1,FOLR2,C1QC,APOE,MSR1",
    "FCN1,VCAN,CTSS,LYZ,S100A8,S100A9,LILRB1",
    "FCN1,S100A8,S100A9,IL1B,TNF,CXCL8,NFKBIA",
    "FCGR3A,IFITM3,SAT1,TYMP,LST1,IFI30,MS4A7",
    "CD1C,CLEC10A,FCER1A,HLA-DRA,HLA-DPA1,CST3",
    "CD1A,CD207,FCER1A,CLEC10A,HLA-DRA,CST3",
    "GZMB,IRF7,LILRA4,CLEC4C,TCF4,JCHAIN,IFITM1",
    "KIT,TPSAB1,TPSB2,CPA3,HPGDS,MS4A2,HDC",
    "FCGR3B,CXCL8,CSF3R,NAMPT,S100A8,S100A9,FPR1"
  ),
  stringsAsFactors = FALSE
)

# ----- Myeloid-specific tissue context for LLM -----
LINEAGE_BASE_CONTEXT <- paste(
  "Myeloid cells from NORMAL (non-diseased) human respiratory tract tissues.",
  "This is a cross-site anatomical comparison, NOT a disease vs healthy comparison.",
  "Key populations include alveolar macrophages, interstitial macrophages,",
  "classical and non-classical monocytes, dendritic cells, pDC, mast cells,",
  "and neutrophils.",
  "Focus: regional variation in innate immune surveillance, tissue residency,",
  "antigen presentation, inflammatory recruitment, and mucosal barrier defense."
)

TISSUE_CONTEXT <- list(
  "nose" = paste(
    "Nasal cavity: first-line mucosal barrier with high microbial and particulate exposure,",
    "often enriched for recruited monocytes, neutrophils, mast cells, and antigen-sampling DCs.",
    "Myeloid programs here may emphasize rapid sensing, cytokine response, and epithelial interaction."
  ),
  "sinus" = paste(
    "Paranasal sinus: semi-enclosed mucosal compartment with drainage-dependent clearance.",
    "Myeloid cells here may reflect persistent barrier surveillance, antigen retention,",
    "and macrophage/DC states adapted to low-flow mucosal niches."
  ),
  "respiratory airway" = paste(
    "Conducting airways (trachea/bronchi): mucociliary epithelium exposed to inhaled antigens.",
    "Myeloid compartments here may balance monocyte recruitment, dendritic-cell antigen presentation,",
    "and epithelial-supportive macrophage programs."
  ),
  "lung parenchyma" = paste(
    "Lung parenchyma (alveolar region): gas-exchange surface with resident alveolar macrophages as a dominant",
    "innate immune population, plus interstitial macrophages and DCs for tissue homeostasis and antigen sampling.",
    "Programs here often reflect tissue residency, lipid handling, phagocytosis, and restrained inflammation."
  )
)

# ----- Myeloid changelog text -----
PIPELINE_CHANGELOG_TITLE <- sprintf("%s Changes:", PIPELINE_VERSION_LABEL)
PIPELINE_CHANGELOG_LINES <- c(
  "  [MYELOID-1] Uses `cell_type_L3_refined` as the standardized L3 source column.",
  "  [MYELOID-2] Rebuilds cleaner L2 labels from refined L3 to fix inconsistencies in legacy `cell_type_L2`.",
  "  [MYELOID-3] Adds myeloid-specific marker panels, enrichment filters, and custom marker database.",
  "  [MYELOID-4] Syncs to shared v2.6.0 engine with pseudobulk cache + validated gene->pathway evidence.",
  "  [MYELOID-5] Enables ssGSEA group LLM outputs with combined multi-database evidence (LLM_SSGSEA_INTERPRETATION.md).",
  "  [MYELOID-6] Enables CHOIR per-cluster LLM outputs (LLM_CHOIR_INTERPRETATION.md).",
  "  [MYELOID-7] LLM interpretation now synthesizes multiple enrichment databases in a single call.",
  "  [MYELOID-8] CHOIR alpha increased to 0.20 and variable features capped to reduce over-clustering and runtime."
)
PIPELINE_INHERITED_FIXES_TITLE <- "Shared engine fixes inherited from the validated template:"
PIPELINE_INHERITED_FIX_LINES <- c(
  "  [ENGINE-1] Sparse-friendly pseudobulk aggregation with one-shot cache per level.",
  "  [ENGINE-2] Unmapped L3 labels fail fast instead of being silently dropped.",
  "  [ENGINE-3] interpret_agent retry + malformed-output normalization to fixed JSON schema.",
  "  [ENGINE-4] UMAP/heatmap export tuned for publication-friendly PDF+PNG outputs.",
  "  [ENGINE-5] Exploratory wilcox outputs explicitly labeled as cell-level exploratory evidence.",
  "  [ENGINE-6] ssGSEA and CHOIR layers now support dedicated LLM interpretation outputs."
)

# ----- Run the shared engine -----
source("/home/h2048/script/R/bcell_tissue_comparison_v2_6_20260406.R")
