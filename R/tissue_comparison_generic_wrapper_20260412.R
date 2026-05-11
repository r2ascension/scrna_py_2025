#!/usr/bin/env Rscript
# ==============================================================================
# Generic Tissue Comparison Wrapper Interface
# ==============================================================================
#
# Purpose:
#   - Provide an explicit, reusable wrapper interface around the shared
#     tissue-comparison engine(s)
#   - Centralize wrapper-owned env loading, previous-run preflight, and
#     config validation
#   - Make lineage-specific entrypoints thin and declarative
#
# Date: 2026-04-12
# ==============================================================================

ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
if (!file.exists(ADVANCED_HELPER_PATH)) {
  stop(sprintf("Advanced helper not found: %s", ADVANCED_HELPER_PATH))
}
source(ADVANCED_HELPER_PATH)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) return(y)
  x
}

tc_wrapper_env_flag <- function(name, default = FALSE) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) return(isTRUE(default))
  tolower(raw) %in% c("1", "true", "t", "yes", "y", "on")
}

tc_wrapper_scalar_env <- function(name, default = NULL) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) return(default)
  raw
}

tc_wrapper_is_placeholder_secret <- function(x) {
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

tc_generic_wrapper_base_config <- function() {
  list(
    ENV_FILE_CANDIDATES = c("/home/h2048/.env", "/home/h2048/script/.env"),
    REUSE_PREVIOUS_FINAL_OBJECT = TRUE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    REQUIRE_PREVIOUS_OUTPUT_DIR = TRUE,
    ALLOW_IN_PLACE_RESUME = FALSE,
    REQUIRE_LLM = TRUE,
    CREATE_ENV_PLACEHOLDER_IF_MISSING = TRUE,
    ADVANCED_HELPER_PATH = ADVANCED_HELPER_PATH,
    ADVANCED_HELPER_ALREADY_LOADED = TRUE,
    SKIP_ENV_AUTOLOAD = TRUE,
    SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE = TRUE,
    RUN_INITIAL_SUPPORT = TRUE,
    INITIAL_SUPPORT_STRICT = FALSE,
    INITIAL_SUPPORT_SUBDIR = "support",
    INITIAL_SUPPORT_ROGUE_MAX_CELLS = 1000L,
    SHARED_OVERRIDES = list()
  )
}

tc_bcell_wrapper_preset <- function() {
  list(
    LINEAGE_TAG = "BCELL",
    LINEAGE_DISPLAY = "B Cell",
    LINEAGE_CONTEXT_LABEL = "B cells and plasma cells",
    LINEAGE_CONTEXT_LOWER = "B-cell",
    BIOLOGICAL_QUESTION_FRAGMENT = "B cell biology, mucosal immunity, or tissue microenvironment",
    PIPELINE_VERSION_LABEL = "v2.6.4",
    PIPELINE_SUBTITLE = paste(
      "B Cell Tissue Comparison v2.6.4",
      "(generic wrapper interface)"
    ),
    REPORT_TITLE = "# B Cell Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "bcell_tissue_comparison_v2_6_4_20260412.R",
    FINAL_FILE_PREFIX = "bcell_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "B CELL TISSUE COMPARISON COMPLETE (v2.6.4)",
    H5AD_PATH = "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0412/bcell_tissue_comparison_v2_6_4_20260412",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0408/bcell_tissue_comparison_v2_6_20260408",
    SSGSEA_METHODS = c("hallmark", "go_bp"),
    SSGSEA_CHOIR_METHODS = c("hallmark", "go_bp"),
    FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT = FALSE,
    FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    DEPRIORITIZE_ENERGY_PATHWAYS = TRUE,
    DEPRIORITIZE_TECHNICAL_PATHWAYS = TRUE,
    DEPRIORITIZE_CILIA_PATHWAYS = TRUE,
    APPEND_ENERGY_PATHWAYS_AFTER_TOP = FALSE,
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_2_20260414.R",
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
      LLM_EXTRA_RULES = c(
        "Keep each field compact and synthesis-first; never split the final answer into separate up/down or per-database mini-conclusions.",
        "Do not let energy metabolism, mitochondrial, ribosomal/translation, or ENSG-like placeholder pathways occupy top-pathway slots.",
        "For non-epithelial interpretation contexts, exclude cilia/ciliogenesis-related genes and pathways from top-gene emphasis and pathway reasoning.",
        "For non-B-cell interpretation contexts, immunoglobulin genes should not be used as primary evidence.",
        "Prefer lineage-informative biology over housekeeping-like programs when summarizing top pathways."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface provides explicit config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] B-cell entrypoint is now a thin preset over the shared generic wrapper.",
      "  [ARCH-3] OUTPUT_DIR and PREVIOUS_OUTPUT_DIR remain separated by default; in-place resume requires explicit opt-in.",
      "  [SUPPORT-1] Upfront Ro/e + ggtree + ROGUE diagnostics now run before the main analysis and write to reports/support.",
      "  [FINAL-1] Compact integrated LLM outputs only; no split final conclusions by up/down or per database.",
      "  [FINAL-2] Non-informative pathways (energy / mitochondrial / ribosomal / ENSG-like) do not occupy top-pathway slots.",
      "  [FINAL-3] Non-epithelial cilia-related genes/pathways are suppressed from enrichment and LLM emphasis.",
      "  [FINAL-4] Previous output reading remains compatible with legacy table/RDS layouts."
    )
  )
}

tc_tnk_wrapper_preset <- function() {
  list(
    LINEAGE_TAG = "TNK",
    LINEAGE_DISPLAY = "T/NK",
    LINEAGE_CONTEXT_LABEL = "T cells, NK cells, and innate lymphoid cells",
    LINEAGE_CONTEXT_LOWER = "T/NK-cell",
    BIOLOGICAL_QUESTION_FRAGMENT = paste(
      "T/NK cell residency, helper-vs-cytotoxic polarization, innate-like lymphocyte programs,",
      "or tissue microenvironment"
    ),
    PIPELINE_VERSION_LABEL = "v2.6.1-TNK",
    PIPELINE_SUBTITLE = paste(
      "T/NK Tissue Comparison v2.6.1-TNK",
      "(generic wrapper interface)"
    ),
    REPORT_TITLE = "# T/NK Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "tnk_tissue_comparison_v2_6_1_20260413.R",
    FINAL_FILE_PREFIX = "tnk_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "T/NK TISSUE COMPARISON COMPLETE (v2.6.1-TNK)",
    H5AD_PATH = "/home/h2048/data/py/0318/tnk_subcluster_retrain/adata_tnk_scanvi_ref_retrain_v1_2.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0",
    L3_SOURCE_COL = "scanvi_label_refined",
    L2_SOURCE_COL = "cell_type_L2",
    USE_EXISTING_L2 = FALSE,
    ANALYSIS_L3_DESCRIPTION = "refined scANVI labels",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`scanvi_label_refined`)",
    L3_TO_L2_REMAP = c(
      "CD4 Naive/TCM" = "CD4 T cells",
      "CD4 Tcm" = "CD4 T cells",
      "CD4 Tfh" = "CD4 T cells",
      "CD4 Tfr" = "CD4 T cells",
      "CD4 Th1" = "CD4 T cells",
      "CD4 Th17" = "CD4 T cells",
      "CD4 Treg" = "CD4 T cells",
      "CD4 Trm" = "CD4 T cells",
      "CD8 Naive" = "CD8 T cells",
      "CD8 Teff" = "CD8 T cells",
      "CD8 Tem" = "CD8 T cells",
      "CD8 Temra" = "CD8 T cells",
      "CD8 Trm" = "CD8 T cells",
      "gdT" = "CD8 T cells",
      "MAIT" = "CD8 T cells",
      "ILC3" = "NK cells",
      "NK" = "NK cells",
      "NK Exhausted" = "NK cells"
    ),
    UMAP_REDUCTION_PREFERRED = c("umap_refined", "umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap"),
    CHOIR_REDUCTION_CANDIDATES = c("scanvi_refined", "scanvi", "scvi", "harmony", "pca"),
    CUSTOM_DB_NAME = "TNK_custom",
    CUSTOM_DB_LABEL = "T/NK custom markers",
    CUSTOM_ENRICHMENT_SIZE_RULE = c(min = 2L, max = 120L),
    CELLMARKER_CELLTYPE_PATTERN = paste(
      "T cell|T-cell|CD4 T|CD8 T|Helper T|Cytotoxic T|Regulatory T|Treg|",
      "Th1|Th17|Tfh|MAIT|gamma delta|γδ T|Innate lymphoid|ILC|NK cell|Natural killer"
    ),
    PANGLAODB_CELLTYPE_PATTERN = paste(
      "T cell|CD4 T|CD8 T|cytotoxic T|helper T|regulatory T|Treg|",
      "Th1|Th17|Tfh|MAIT|gamma delta|NK cell|natural killer|ILC"
    ),
    KNOWN_MARKERS = unique(c(
      "CD3D", "CD3E", "CD247", "TRAC", "IL7R",
      "CCR7", "SELL", "TCF7", "LEF1", "LTB",
      "CXCR5", "PDCD1", "ICOS", "BCL6", "MAF",
      "TBX21", "IFNG", "CXCR3", "RORC", "IL17A", "CCR6",
      "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
      "CD8A", "CD8B", "KLRG1", "FGFBP2", "PRF1", "GZMB", "GZMK",
      "NKG7", "GNLY", "CX3CR1",
      "CD69", "ITGAE", "CXCR6", "ZNF683", "XCL1",
      "TRDC", "TRGC1", "TRGC2", "TRAV1-2", "KLRB1", "SLC4A10", "ZBTB16",
      "FCGR3A", "KLRD1", "NCR1", "XCL2", "CCL3", "CCL4",
      "HAVCR2", "LAG3", "TOX", "LAYN",
      "KIT", "AHR", "IL23R", "NCR2"
    )),
    CUSTOM_MARKERS_DB = data.frame(
      subtype = c(
        "CD4 Naive/TCM", "CD4 Tcm", "CD4 Tfh", "CD4 Tfr", "CD4 Th1", "CD4 Th17", "CD4 Treg", "CD4 Trm",
        "CD8 Naive", "CD8 Teff", "CD8 Tem", "CD8 Temra", "CD8 Trm", "ILC3", "MAIT", "NK", "NK Exhausted", "gdT"
      ),
      markers = c(
        "CCR7,SELL,TCF7,LEF1,IL7R,LTB,MAL",
        "CCR7,SELL,TCF7,LEF1,LTB,MAL,IL7R,MHC2TA",
        "CXCR5,PDCD1,ICOS,BCL6,MAF,IL21R,TOX2",
        "CXCR5,FOXP3,IL2RA,CTLA4,TIGIT,PDCD1,MAF",
        "TBX21,IFNG,CXCR3,STAT1,CCL5,GZMK,IL7R",
        "RORC,IL17A,IL17F,CCR6,KLRB1,IL23R,AHR",
        "FOXP3,IL2RA,CTLA4,TIGIT,IKZF2,TNFRSF4,LAYN",
        "CD69,ITGAE,CXCR6,ZNF683,IL7R,XCL1,RGS1",
        "CCR7,SELL,TCF7,LEF1,LTB,IL7R,MAL",
        "PRF1,GZMB,NKG7,GNLY,CCL5,FGFBP2,CX3CR1",
        "GZMK,CCL5,NKG7,IL7R,CXCR3,ANXA1",
        "KLRG1,FGFBP2,PRF1,GZMB,CX3CR1,NKG7,FCGR3A",
        "CD69,ITGAE,CXCR6,ZNF683,XCL1,IFNG,NKG7",
        "RORC,IL7R,KIT,AHR,IL23R,NCR2,LTB",
        "TRAV1-2,KLRB1,SLC4A10,ZBTB16,RORA,IL7R,NKG7",
        "NKG7,GNLY,PRF1,GZMB,FCGR3A,KLRD1,NCR1,XCL2",
        "TIGIT,LAG3,HAVCR2,TOX,LAYN,PDCD1,CTLA4",
        "TRDC,TRGC1,TRGC2,NKG7,KLRD1,CCL5,GNLY"
      ),
      stringsAsFactors = FALSE
    ),
    LINEAGE_BASE_CONTEXT = paste(
      "T cells, NK cells, and innate lymphoid cells from NORMAL (non-diseased) human respiratory tract tissues.",
      "This is a cross-site anatomical comparison, NOT a disease vs healthy comparison.",
      "Key lymphocyte programs include CD4 helper / Treg states, CD8 cytotoxic / resident memory states,",
      "innate-like T cells (MAIT, gamma-delta), NK cell effector states, and ILC3 programs.",
      "Focus: regional variation in tissue residency, effector polarization, cytotoxic surveillance,",
      "and mucosal immune adaptation along the respiratory tract."
    ),
    TISSUE_CONTEXT = list(
      "nose" = paste(
        "Nasal cavity: first-line mucosal barrier with intense environmental antigen exposure,",
        "enriched for tissue-resident memory T cells, rapid innate lymphocyte surveillance,",
        "and local helper / regulatory programs that tune mucosal defense."
      ),
      "sinus" = paste(
        "Paranasal sinus: semi-enclosed mucosal cavity with drainage-dependent antigen clearance,",
        "supports compartmentalized resident lymphocyte niches and barrier-immune surveillance,",
        "with balanced helper, regulatory, and cytotoxic programs in normal tissue."
      ),
      "respiratory airway" = paste(
        "Conducting airways (trachea/bronchi): interface between upper and lower airway immunity,",
        "supports airway TRM cells, NK surveillance, and rapid epithelial-adjacent effector responses,",
        "while maintaining tolerance to inhaled exposures."
      ),
      "lung parenchyma" = paste(
        "Lung parenchyma (alveolar region): gas exchange surface with tightly regulated immune tone,",
        "contains interstitial T cells, alveolar-adjacent resident memory populations, and cytotoxic surveillance,",
        "with stronger systemic-like immune constraints than upper airway mucosa."
      )
    ),
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_2_20260414.R",
    SHARED_OVERRIDES = list(
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
      LLM_EXTRA_RULES = c(
        "Start overview with a concise cell-type/state judgment.",
        "Use ssGSEA mainly to support cell-type judgment; avoid over-interpreting pathways as direct mechanisms.",
        "For CHOIR clusters, integrate up/down DEG and cluster-level ssGSEA into one compact interpretation.",
        "Prefer lineage-informative lymphocyte programs over broad stress/metabolic signals when ranking discoveries."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] TNK entrypoint is now a thin preset over the shared generic wrapper.",
      "  [SUPPORT-1] Upfront Ro/e + ggtree + ROGUE diagnostics now run before the main analysis and write to reports/support.",
      "  [TNK-1] Uses `scanvi_label_refined` as L3 and rebuilds deterministic L2 labels from refined TNK states.",
      "  [TNK-2] Reuses helper-owned TNK marker panels and shared ssGSEA / CHOIR / discovery-screening logic.",
      "  [TNK-3] Enables additive LLM outputs for ssGSEA review plus family-level discovery/outlier screening.",
      "  [TNK-4] OUTPUT_DIR and PREVIOUS_OUTPUT_DIR remain separated by default; in-place resume requires explicit opt-in."
    ),
    PIPELINE_INHERITED_FIXES_TITLE = "Shared engine fixes inherited from the validated template:",
    PIPELINE_INHERITED_FIX_LINES = c(
      "  [ENGINE-1] Sparse-friendly pseudobulk aggregation with one-shot cache per level.",
      "  [ENGINE-2] Unmapped L3 labels fail fast instead of being silently dropped.",
      "  [ENGINE-3] interpret_agent retry + malformed-output normalization to fixed JSON schema.",
      "  [ENGINE-4] ssGSEA and CHOIR layers support annotation-match / outlier / discovery review outputs."
    )
  )
}

tc_epithelial_wrapper_preset <- function() {
  list(
    LINEAGE_TAG = "EPITHELIAL",
    LINEAGE_DISPLAY = "Epithelial",
    LINEAGE_CONTEXT_LABEL = "respiratory epithelial cells",
    LINEAGE_CONTEXT_LOWER = "epithelial",
    BIOLOGICAL_QUESTION_FRAGMENT = paste(
      "epithelial differentiation, barrier defense, mucus biology, ion transport,",
      "glandular secretion, or alveolar specialization"
    ),
    PIPELINE_VERSION_LABEL = "v1.3.2-EPI",
    PIPELINE_SUBTITLE = paste(
      "Epithelial Tissue Comparison v1.3.2-EPI",
      "(generic wrapper interface)"
    ),
    REPORT_TITLE = "# Epithelial Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "epithelial_tissue_comparison_v1_3_2_20260413.R",
    FINAL_FILE_PREFIX = "epithelial_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "EPITHELIAL TISSUE COMPARISON COMPLETE (v1.3.2-EPI)",
    H5AD_PATH = "/home/h2048/data/core20260322/epithelial_scanvi_v2_7_HOTFIX_SELF_final.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0413/epithelial_tissue_comparison_v1_3_2_20260413",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0407/epithelial_tissue_comparison_v1_3_1",
    L3_SOURCE_COL = "cell_type_L3",
    ANALYSIS_L3_DESCRIPTION = "curated epithelial L3 labels",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_L3`)",
    L3_TO_L2_REMAP = c(
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
    ),
    UMAP_REDUCTION_PREFERRED = c(
      "umap_fine", "umap_major",
      "umap_scanvi_fine", "umap_scanvi_major",
      "umap_scanvi", "umap_harmony", "umap"
    ),
    CHOIR_REDUCTION_CANDIDATES = c("scanvi_major", "scanvi", "scvi", "harmony", "pca", "scanvi_fine"),
    CHOIR_ALPHA = 0.20,
    CHOIR_VAR_FEATURES_MAX = 4000L,
    CHOIR_FIND_VAR_FEATURES_IF_MISSING = TRUE,
    CHOIR_VAR_FEATURES_METHOD = "vst",
    MIN_CELLS_PER_PSEUDOBULK = 15L,
    WILCOX_MIN_CELLS = 100L,
    HEATMAP_CELLS_PER_TYPE = 120L,
    PSEUDOBULK_DISAMBIGUATION_CANDIDATES = c("dataset", "source", "batch", "orig.ident"),
    CUSTOM_DB_NAME = "Epithelial_custom",
    CUSTOM_DB_LABEL = "Epithelial custom markers",
    CUSTOM_ENRICHMENT_SIZE_RULE = c(min = 2L, max = 100L),
    FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT = TRUE,
    FILTER_CILIA_GENES_FOR_LLM_AND_ENRICHMENT = FALSE,
    DEPRIORITIZE_ENERGY_PATHWAYS = TRUE,
    DEPRIORITIZE_TECHNICAL_PATHWAYS = TRUE,
    DEPRIORITIZE_CILIA_PATHWAYS = FALSE,
    APPEND_ENERGY_PATHWAYS_AFTER_TOP = FALSE,
    CELLMARKER_CELLTYPE_PATTERN = paste(
      "Epithelial|Basal|Goblet|Club|Secretory|Ciliated|Ionocyte|Brush|Tuft|Deuterosomal|Serous|Duct|",
      "Airway|Respiratory|Alveolar|AT1|AT2|Suprabasal|Squamous"
    ),
    PANGLAODB_CELLTYPE_PATTERN = paste(
      "Epithelial|Basal|Goblet|Club|Secretory|Ciliated|Ionocyte|Brush|Tuft|Serous|Duct|Airway|Alveolar|AT1|AT2|Suprabasal|Squamous"
    ),
    KNOWN_MARKERS = c(
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
    ),
    CUSTOM_MARKERS_DB = data.frame(
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
    ),
    LINEAGE_BASE_CONTEXT = paste(
      "Respiratory epithelial cells from NORMAL (non-diseased) human airway and alveolar tissues.",
      "This is a cross-site anatomical comparison, NOT a disease-vs-healthy design.",
      "Key epithelial programs include basal stem/progenitor maintenance, suprabasal transition,",
      "mucociliary differentiation, club/goblet/SMG secretory specialization, ion transport,",
      "alveolar surfactant biology, and barrier defense.",
      "Focus on regional variation in epithelial composition, differentiation state, barrier function,",
      "mucus production, glandular secretion, and alveolar specialization."
    ),
    TISSUE_CONTEXT = list(
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
    ),
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_2_20260414.R",
    SHARED_OVERRIDES = list(
      INTERPRET_MULTI_DB_MIN_TERMS = 3L,
      INTERPRET_MULTI_DB_MAX_DBS = 6L,
      INTERPRET_MULTI_DB_TERMS_PER_DB = 5L,
      INTERPRET_SSGSEA_TERMS_PER_DB = 8L,
      LLM_INCLUDE_TOP_DEG = TRUE,
      LLM_TOP_DEG_N = 10L,
      LLM_DEG_PADJ_THR = 0.10,
      LLM_DEG_LFC_THR = 0.15,
      LLM_SSGSEA_TERMS_PER_DIRECTION = 6L,
      LLM_REQUIRE_INTEGRATED_UP_DOWN = TRUE,
      LLM_SSGSEA_PRIMARY_USE = "cell_type_judgment",
      LLM_ALLOW_COMPARATIVE_HYPOTHESIS = TRUE,
      LLM_EXTRA_RULES = c(
        "Start overview with a concise epithelial lineage/state judgment.",
        "Treat cilia/ciliogenesis programs as biologically meaningful for epithelial interpretation; do not suppress them.",
        "Prefer alveolar, basal, ciliated, secretory, SMG, ionocyte, and squamous differentiation signals over generic stress summaries.",
        "For CHOIR clusters, integrate annotation composition, up/down DEG, ssGSEA, and OFA into one compact epithelial-state interpretation."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] Epithelial entrypoint is now a thin preset over the shared generic wrapper.",
      "  [SUPPORT-1] Upfront Ro/e + ggtree + ROGUE diagnostics now run before the main analysis and write to reports/support.",
      "  [EPI-1] Reuses helper-owned epithelial marker panels and shared ssGSEA / CHOIR / discovery-screening logic.",
      "  [EPI-2] Preserves epithelial-specific cilia biology while still suppressing generic technical/energy noise.",
      "  [EPI-3] Enables pseudobulk sample+tissue disambiguation columns in the shared engine for safer epithelial aggregation.",
      "  [EPI-4] OUTPUT_DIR and PREVIOUS_OUTPUT_DIR remain separated by default; in-place resume requires explicit opt-in."
    ),
    PIPELINE_INHERITED_FIXES_TITLE = "Shared engine fixes inherited from the validated template:",
    PIPELINE_INHERITED_FIX_LINES = c(
      "  [ENGINE-1] Sparse-friendly pseudobulk aggregation with one-shot cache per level.",
      "  [ENGINE-2] Unmapped L3 labels fail fast instead of being silently dropped.",
      "  [ENGINE-3] interpret_agent retry + malformed-output normalization to fixed JSON schema.",
      "  [ENGINE-4] ssGSEA and CHOIR layers support annotation-match / outlier / discovery review outputs."
    )
  )
}

tc_get_lineage_preset <- function(lineage = c("BCELL", "TNK", "EPITHELIAL", "EPI")) {
  lineage <- toupper(tc_safe_trim(lineage %||% "BCELL"))
  switch(
    lineage,
    BCELL = tc_bcell_wrapper_preset(),
    TNK = tc_tnk_wrapper_preset(),
    EPITHELIAL = tc_epithelial_wrapper_preset(),
    EPI = tc_epithelial_wrapper_preset(),
    stop(sprintf("Unsupported lineage preset: %s", lineage))
  )
}

tc_build_generic_tissue_comparison_config <- function(lineage = "BCELL", overrides = list()) {
  cfg <- utils::modifyList(tc_generic_wrapper_base_config(), tc_get_lineage_preset(lineage))
  if (!is.null(overrides) && length(overrides) > 0) {
    cfg <- utils::modifyList(cfg, overrides)
  }

  cfg$REQUIRE_LLM <- tc_wrapper_env_flag("PIPELINE_REQUIRE_LLM", cfg$REQUIRE_LLM %||% TRUE)
  cfg$ALLOW_IN_PLACE_RESUME <- tc_wrapper_env_flag(
    "PIPELINE_ALLOW_IN_PLACE_RESUME",
    cfg$ALLOW_IN_PLACE_RESUME %||% FALSE
  )
  cfg$OUTPUT_DIR <- tc_wrapper_scalar_env("PIPELINE_OUTPUT_DIR", cfg$OUTPUT_DIR)
  cfg$PREVIOUS_OUTPUT_DIR <- tc_wrapper_scalar_env("PIPELINE_PREVIOUS_OUTPUT_DIR", cfg$PREVIOUS_OUTPUT_DIR)
  cfg$H5AD_PATH <- tc_wrapper_scalar_env("PIPELINE_H5AD_PATH", cfg$H5AD_PATH)
  cfg$SHARED_ENGINE_PATH <- tc_wrapper_scalar_env("PIPELINE_SHARED_ENGINE_PATH", cfg$SHARED_ENGINE_PATH)

  cfg$PIPELINE_CHANGELOG_TITLE <- cfg$PIPELINE_CHANGELOG_TITLE %||% sprintf("%s Changes:", cfg$PIPELINE_VERSION_LABEL)
  cfg$PREVIOUS_FINAL_OBJECT_RDS <- file.path(
    cfg$PREVIOUS_OUTPUT_DIR,
    paste0(cfg$FINAL_FILE_PREFIX, ".rds")
  )

  tc_assert_named_list_keys(
    cfg,
    required_keys = c(
      "OUTPUT_DIR", "PREVIOUS_OUTPUT_DIR", "H5AD_PATH", "ENV_FILE_CANDIDATES",
      "REUSE_PREVIOUS_FINAL_OBJECT", "REUSE_PREVIOUS_OUTPUT_SUMMARY",
      "REQUIRE_PREVIOUS_OUTPUT_DIR", "ALLOW_IN_PLACE_RESUME",
      "REQUIRE_LLM", "SHARED_ENGINE_PATH", "SHARED_OVERRIDES",
      "FINAL_FILE_PREFIX", "GENERATED_BY_LABEL", "PIPELINE_VERSION_LABEL"
    ),
    object_name = "PIPELINE_CONFIG"
  )

  cfg
}

tc_run_generic_tissue_comparison <- function(config) {
  tc_assert_named_list_keys(
    config,
    required_keys = c(
      "OUTPUT_DIR", "PREVIOUS_OUTPUT_DIR", "H5AD_PATH", "ENV_FILE_CANDIDATES",
      "REUSE_PREVIOUS_FINAL_OBJECT", "REUSE_PREVIOUS_OUTPUT_SUMMARY",
      "REQUIRE_PREVIOUS_OUTPUT_DIR", "ALLOW_IN_PLACE_RESUME",
      "REQUIRE_LLM", "SHARED_ENGINE_PATH", "SHARED_OVERRIDES",
      "FINAL_FILE_PREFIX", "GENERATED_BY_LABEL", "PIPELINE_VERSION_LABEL",
      "PREVIOUS_FINAL_OBJECT_RDS"
    ),
    object_name = "PIPELINE_CONFIG"
  )

  if (!file.exists(config$SHARED_ENGINE_PATH)) {
    stop(sprintf("Shared engine not found: %s", config$SHARED_ENGINE_PATH))
  }

  if (!file.exists(config$H5AD_PATH) &&
      !(isTRUE(config$REUSE_PREVIOUS_FINAL_OBJECT) && file.exists(config$PREVIOUS_FINAL_OBJECT_RDS))) {
    stop(paste(
      sprintf("H5AD input not found: %s", config$H5AD_PATH),
      "No reusable previous final object was found either, so the wrapper cannot continue."
    ))
  }

  if (isTRUE(config$REUSE_PREVIOUS_FINAL_OBJECT) && !file.exists(config$PREVIOUS_FINAL_OBJECT_RDS)) {
    cat(sprintf("[WARN] Previous final object requested but not found: %s\n", config$PREVIOUS_FINAL_OBJECT_RDS))
  }

  loaded_env_files <- tc_load_env_candidates(config$ENV_FILE_CANDIDATES)
  config$PRELOADED_ENV_FILES <- loaded_env_files

  llm_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  if (nchar(llm_key) < 10 && isTRUE(config$CREATE_ENV_PLACEHOLDER_IF_MISSING)) {
    tc_ensure_env_placeholder(config$ENV_FILE_CANDIDATES[[1]], "DEEPSEEK_API_KEY")
  }
  llm_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  has_llm_key <- nchar(llm_key) >= 10 && !tc_wrapper_is_placeholder_secret(llm_key)
  if (!has_llm_key) {
    if (nchar(llm_key) >= 10 && tc_wrapper_is_placeholder_secret(llm_key)) {
      Sys.unsetenv("DEEPSEEK_API_KEY")
    }
    llm_message <- paste(
      "DEEPSEEK_API_KEY not found or is still a placeholder value in loaded .env files.",
      "Set PIPELINE_REQUIRE_LLM=false if you want a degraded non-LLM run; otherwise provide a valid key."
    )
    if (isTRUE(config$REQUIRE_LLM)) stop(llm_message)
    cat(sprintf("[WARN] %s\n", llm_message))
  } else if (length(loaded_env_files) > 0) {
    cat(sprintf("[OK] Loaded environment file(s): %s\n", paste(loaded_env_files, collapse = ", ")))
  }

  summary_output_path <- file.path(config$OUTPUT_DIR, "reports", "previous_run_summary.tsv")
  preflight <- tc_preflight_previous_run(
    previous_output_dir = config$PREVIOUS_OUTPUT_DIR,
    current_output_dir = config$OUTPUT_DIR,
    allow_in_place_resume = config$ALLOW_IN_PLACE_RESUME,
    require_previous_dir = config$REQUIRE_PREVIOUS_OUTPUT_DIR,
    reuse_previous_output_summary = config$REUSE_PREVIOUS_OUTPUT_SUMMARY,
    summary_output_path = summary_output_path,
    include_rds = FALSE,
    include_markdown = FALSE,
    verbose = TRUE
  )
  config$PREVIOUS_RUN_SUMMARY_OUTPUT_PATH <- summary_output_path
  config$PREVIOUS_RUN_SUMMARY_DF <- preflight$summary

  dir.create(file.path(config$OUTPUT_DIR, "reports"), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      list(
        generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        pipeline_version = config$PIPELINE_VERSION_LABEL,
        generated_by = config$GENERATED_BY_LABEL,
        output_dir = config$OUTPUT_DIR,
        previous_output_dir = config$PREVIOUS_OUTPUT_DIR,
        previous_output_exists = isTRUE(preflight$exists),
        same_dir = isTRUE(preflight$same_dir),
        require_previous_output_dir = isTRUE(config$REQUIRE_PREVIOUS_OUTPUT_DIR),
        require_llm = isTRUE(config$REQUIRE_LLM),
        h5ad_exists = file.exists(config$H5AD_PATH),
        previous_final_object_rds = config$PREVIOUS_FINAL_OBJECT_RDS,
        previous_final_object_exists = file.exists(config$PREVIOUS_FINAL_OBJECT_RDS),
        loaded_env_files = loaded_env_files,
        shared_engine_path = config$SHARED_ENGINE_PATH,
        shared_override_names = names(config$SHARED_OVERRIDES)
      ),
      path = file.path(config$OUTPUT_DIR, "reports", "wrapper_preflight.json"),
      pretty = TRUE,
      auto_unbox = TRUE,
      null = "null"
    )
  }

  execution_env <- new.env(parent = environment())
  tc_apply_named_list(config, envir = execution_env)
  invisible(tc_apply_advanced_shared_overrides(envir = execution_env, overrides = config$SHARED_OVERRIDES))
  source(config$SHARED_ENGINE_PATH, local = execution_env)

  invisible(list(config = config, preflight = preflight, execution_env = execution_env))
}

if (sys.nframe() == 0) {
  if (!exists("PIPELINE_CONFIG", inherits = FALSE)) {
    cat("Generic tissue-comparison wrapper interface loaded.\n")
    cat("Source this file, build PIPELINE_CONFIG with tc_build_generic_tissue_comparison_config(), then call tc_run_generic_tissue_comparison().\n")
  } else {
    tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
  }
}
