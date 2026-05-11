#!/usr/bin/env Rscript
# ==============================================================================
# Generic Tissue Comparison Wrapper Interface (Myeloid Extension)
# ==============================================================================
#
# Purpose:
#   - keep the validated 2026-04-14 wrapper immutable
#   - extend preset coverage to the myeloid branch without touching older files
#   - point wrapper-based myeloid runs at a versioned helper extension
#
# Date: 2026-04-14
# ==============================================================================

BASE_GENERIC_WRAPPER_PATH_V2 <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
if (!file.exists(BASE_GENERIC_WRAPPER_PATH_V2)) {
  stop(sprintf("Base generic wrapper not found: %s", BASE_GENERIC_WRAPPER_PATH_V2))
}
source(BASE_GENERIC_WRAPPER_PATH_V2)

ADVANCED_HELPER_PATH_V2 <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260414_v2.R"
if (!file.exists(ADVANCED_HELPER_PATH_V2)) {
  stop(sprintf("Extended advanced helper not found: %s", ADVANCED_HELPER_PATH_V2))
}
source(ADVANCED_HELPER_PATH_V2)

tc_get_lineage_preset_base_20260414_v2 <- tc_get_lineage_preset

tc_myeloid_shared_overrides <- function(extra_rules = character()) {
  list(
    INTERPRET_MULTI_DB_MIN_TERMS = 1L,
    INTERPRET_MULTI_DB_MAX_DBS = 5L,
    INTERPRET_MULTI_DB_TERMS_PER_DB = 4L,
    INTERPRET_SSGSEA_TERMS_PER_DB = 6L,
    LLM_INCLUDE_TOP_DEG = TRUE,
    LLM_TOP_DEG_N = 10L,
    LLM_DEG_PADJ_THR = 0.10,
    LLM_DEG_LFC_THR = 0.15,
    LLM_SSGSEA_TERMS_PER_DIRECTION = 6L,
    LLM_REQUIRE_INTEGRATED_UP_DOWN = TRUE,
    LLM_SSGSEA_PRIMARY_USE = "cell_type_judgment",
    LLM_ALLOW_COMPARATIVE_HYPOTHESIS = TRUE,
    FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    DEPRIORITIZE_ENERGY_PATHWAYS = TRUE,
    DEPRIORITIZE_TECHNICAL_PATHWAYS = TRUE,
    DEPRIORITIZE_CILIA_PATHWAYS = TRUE,
    APPEND_ENERGY_PATHWAYS_AFTER_TOP = FALSE,
    LLM_EXTRA_RULES = c(
      "Start overview with a concise myeloid state judgment.",
      "Use ssGSEA mainly to support cell-state judgment; avoid converting broad pathway shifts into unsupported direct mechanisms.",
      "Prefer resident macrophage, inflammatory monocyte, dendritic antigen-presentation, mast-cell, and neutrophil programs over generic stress or metabolic summaries when ranking discoveries.",
      "For CHOIR clusters, integrate annotation composition, up/down DEG, ssGSEA, and OFA evidence into one compact myeloid-state interpretation.",
      extra_rules
    )
  )
}

tc_myeloid_wrapper_preset <- function() {
  known_markers <- unique(c(
    "LYZ", "LST1", "TYMP", "CTSS", "FCER1G", "HLA-DRA",
    "FCN1", "VCAN", "S100A8", "S100A9", "CTSD", "SAT1", "FCGR3A", "IFITM3",
    "FABP4", "PPARG", "MARCO", "C1QA", "C1QB", "C1QC", "INHBA",
    "APOE", "CD163", "CD163L1", "MRC1", "FOLR2",
    "FCER1A", "CD1C", "CLEC10A", "CLEC4C", "LILRA4", "GZMB", "IRF7",
    "KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2",
    "FCGR3B", "CXCL8", "CSF3R", "NAMPT"
  ))

  custom_markers_db <- data.frame(
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

  marker_panels <- tc_default_myeloid_marker_panels(
    known_markers = known_markers,
    custom_markers_db = custom_markers_db
  )

  list(
    LINEAGE_TAG = "MYELOID",
    LINEAGE_DISPLAY = "Myeloid",
    LINEAGE_CONTEXT_LABEL = "myeloid cells",
    LINEAGE_CONTEXT_LOWER = "myeloid",
    BIOLOGICAL_QUESTION_FRAGMENT = paste(
      "innate immune surveillance, monocyte/macrophage tissue adaptation,",
      "antigen presentation, interferon tone, inflammatory recruitment,",
      "or tissue-resident phagocyte specialization"
    ),
    PIPELINE_VERSION_LABEL = "v1.2.1-MYELOID",
    PIPELINE_SUBTITLE = paste(
      "Myeloid Tissue Comparison v1.2.1-MYELOID",
      "(generic wrapper interface, refined L3 remap)"
    ),
    REPORT_TITLE = "# Myeloid Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "myeloid_tissue_comparison_v1_2_1_20260414.R",
    FINAL_FILE_PREFIX = "myeloid_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "MYELOID TISSUE COMPARISON COMPLETE (v1.2.1-MYELOID)",
    H5AD_PATH = "/home/h2048/data/py/0329/adata_myeloid_L3refined_patched_v1.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_1_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0407/myeloid_tissue_comparison_v1_2",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    L3_SOURCE_COL = "cell_type_L3_refined",
    L2_SOURCE_COL = "cell_type_L2",
    USE_EXISTING_L2 = FALSE,
    ANALYSIS_L3_DESCRIPTION = "fine refined myeloid states from `cell_type_L3_refined`",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 refined (`cell_type_L3_refined`)",
    L3_TO_L2_REMAP = c(
      "Resident Alveolar macrophages" = "Alveolar_Macrophage",
      "Resting Alveolar macrophages" = "Alveolar_Macrophage",
      "Interstitial macrophages" = "Interstitial_Macrophage",
      "Inflammatory Interstitial macrophages" = "Interstitial_Macrophage",
      "CD163L1+ Interstitial macrophages" = "Interstitial_Macrophage",
      "Typical Classical monocytes" = "Monocyte",
      "Inflammatory Classical monocytes" = "Monocyte",
      "Non-classical monocytes" = "Monocyte",
      "Conventional cDC2" = "DC",
      "Langerhans-like cDC2" = "DC",
      "pDC" = "DC",
      "Mast cells" = "Mast_cell",
      "Neutrophils" = "Neutrophil"
    ),
    CHOIR_ALPHA = 0.20,
    CHOIR_VAR_FEATURES_MAX = 4000L,
    PSEUDOBULK_DISAMBIGUATION_CANDIDATES = c("dataset", "source", "batch", "orig.ident"),
    CUSTOM_DB_NAME = "Myeloid_custom",
    CUSTOM_DB_LABEL = "Myeloid custom markers",
    CUSTOM_ENRICHMENT_SIZE_RULE = c(min = 2L, max = 120L),
    CELLMARKER_CELLTYPE_PATTERN = paste(
      "Myeloid|Monocyte|Macrophage|Alveolar macrophage|Interstitial macrophage|",
      "Dendritic|DC|cDC|pDC|Plasmacytoid DC|Mast cell|Neutrophil|Langerhans"
    ),
    PANGLAODB_CELLTYPE_PATTERN = paste(
      "Myeloid|Monocyte|Macrophage|Alveolar macrophage|Interstitial macrophage|",
      "Dendritic|DC|cDC|pDC|Mast cell|Neutrophil|Langerhans"
    ),
    KNOWN_MARKERS = known_markers,
    CUSTOM_MARKERS_DB = custom_markers_db,
    MARKER_PANELS = marker_panels,
    LINEAGE_BASE_CONTEXT = paste(
      "Myeloid cells from NORMAL (non-diseased) human respiratory tract tissues.",
      "This is a cross-site anatomical comparison, NOT a disease-vs-healthy comparison.",
      "Key populations include alveolar macrophages, interstitial macrophages,",
      "classical and non-classical monocytes, dendritic cells, pDC, mast cells,",
      "and neutrophils.",
      "Focus: regional variation in innate immune surveillance, tissue residency,",
      "antigen presentation, inflammatory recruitment, and mucosal barrier defense."
    ),
    TISSUE_CONTEXT = list(
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
    ),
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    SHARED_OVERRIDES = tc_myeloid_shared_overrides(
      extra_rules = c(
        "Prefer alveolar-resident, interstitial-resident, inflammatory-monocyte, and dendritic antigen-presentation programs over generic interferon-only summaries when ranking findings.",
        "If neutrophil or mast-cell signals dominate a cluster, state that explicitly instead of forcing macrophage or DC language."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Myeloid wrapper is now a thin preset over a versioned generic wrapper extension instead of a standalone sourced-engine script.",
      "  [ARCH-2] New helper extension owns myeloid marker panels without modifying older helper versions.",
      "  [MYELOID-1] Uses `cell_type_L3_refined` as the standardized L3 source column.",
      "  [MYELOID-2] Rebuilds cleaner L2 labels from refined L3 to fix inconsistencies in legacy `cell_type_L2`.",
      "  [MYELOID-3] Forces fresh H5AD loading and full rerun; 0407 myeloid v1.2 is used only as the previous-run baseline summary.",
      "  [MYELOID-4] Wrapper requires a valid DEEPSEEK key, so the rerun cannot silently skip LLM steps."
    ),
    PIPELINE_INHERITED_FIXES_TITLE = "Shared engine fixes inherited from the validated template:",
    PIPELINE_INHERITED_FIX_LINES = c(
      "  [ENGINE-1] Sparse-friendly pseudobulk aggregation with sample+tissue disambiguation when needed.",
      "  [ENGINE-2] Unmapped L3 labels fail fast instead of being silently dropped.",
      "  [ENGINE-3] interpret_agent retry + malformed-output normalization to fixed JSON schema.",
      "  [ENGINE-4] ssGSEA and CHOIR layers support annotation-match / outlier / discovery review outputs.",
      "  [ENGINE-5] Previous output reading remains compatible with legacy table/RDS layouts."
    )
  )
}

tc_get_lineage_preset <- function(lineage = c(
  "BCELL", "TNK", "EPITHELIAL", "EPI",
  "STROMAL_ENDOTHELIAL", "STROMAL_FIBROBLAST", "STROMAL_SMC",
  "MYELOID"
)) {
  lineage <- toupper(tc_safe_trim(lineage %||% "BCELL"))
  if (identical(lineage, "MYELOID")) return(tc_myeloid_wrapper_preset())
  tc_get_lineage_preset_base_20260414_v2(lineage)
}

if (sys.nframe() == 0) {
  if (!exists("PIPELINE_CONFIG", inherits = FALSE)) {
    cat("Generic tissue-comparison wrapper extension (2026-04-14 v2) loaded.\n")
    cat("Supported extra lineage: MYELOID\n")
  } else {
    tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
  }
}
