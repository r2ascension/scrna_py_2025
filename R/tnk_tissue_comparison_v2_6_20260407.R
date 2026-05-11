#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Tissue Comparison Pipeline - v2.6.0-TNK
# ==============================================================================
#
# Purpose:
#   T/NK lineage tissue comparison across normal respiratory tract sites.
#   This script reuses the validated v2.6 B-cell tissue-comparison backbone,
#   while overriding lineage-specific configuration for TNK data:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. Use refined scanvi_label_refined as L3 and rebuild deterministic L2
#     3. Visualization (UMAP, marker dotplot/heatmap, composition)
#     4. Pseudobulk DESeq2 per L2/L3 subtype x tissue pair
#     5. Exploratory cell-level wilcox per L2/L3 subtype
#     6. Multi-database enrichment + interpret_agent
#     7. Pseudobulk ssGSEA + per-group LLM
#     8. CHOIR + OFA + per-cluster LLM
#     9. Structured REPORT.md + RDS/h5ad export
#
# Notes:
#   - L3 is set to `scanvi_label_refined` (18 refined TNK states).
#   - L2 is rebuilt from refined L3 labels to enforce a deterministic,
#     lineage-consistent coarse grouping for L2-level analysis.
#
# Smoke test:
#   PIPELINE_TEST_MODE=true Rscript tnk_tissue_comparison_v2_6_20260407.R
#   -> exits after load + metadata validation + NormalizeData()
#
# Date: 2026-04-07
# ==============================================================================

LINEAGE_TAG           <- "TNK"
LINEAGE_DISPLAY       <- "T/NK"
LINEAGE_CONTEXT_LABEL <- "T cells, NK cells, and innate lymphoid cells"
LINEAGE_CONTEXT_LOWER <- "T/NK-cell"
BIOLOGICAL_QUESTION_FRAGMENT <- paste(
  "T/NK cell residency, helper-vs-cytotoxic polarization, innate-like lymphocyte programs,",
  "or tissue microenvironment"
)
PIPELINE_VERSION_LABEL <- "v2.6.0-TNK"
PIPELINE_SUBTITLE      <- "T/NK Tissue Comparison v2.6.0-TNK (pseudobulk DESeq2, multi-database enrichment, ssGSEA, CHOIR, LLM)"
REPORT_TITLE           <- "# T/NK Tissue Comparison Report (Normal Respiratory Tract)"
GENERATED_BY_LABEL     <- "tnk_tissue_comparison_v2_6_20260407.R"
FINAL_FILE_PREFIX      <- "tnk_tissue_comparison_final"
LINEAGE_COMPLETION_BANNER <- "T/NK TISSUE COMPARISON COMPLETE (v2.6.0-TNK)"

# ----- Input / Output -----
H5AD_PATH  <- "/home/h2048/data/py/0318/tnk_subcluster_retrain/adata_tnk_scanvi_ref_retrain_v1_2.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0"

# ----- Label columns -----
L3_SOURCE_COL <- "scanvi_label_refined"
L2_SOURCE_COL <- "cell_type_L2"
USE_EXISTING_L2 <- FALSE
ANALYSIS_L3_DESCRIPTION <- "refined scANVI labels"
L3_TO_L2_TABLE_HEADER_LEFT <- sprintf("L3 (`%s`)", L3_SOURCE_COL)
L3_TO_L2_REMAP <- c(
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
)

# ----- Reduction preferences -----
UMAP_REDUCTION_PREFERRED <- c("umap_refined", "umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap")
CHOIR_REDUCTION_CANDIDATES <- c("scanvi_refined", "scanvi", "scvi", "harmony", "pca")

# ----- Database naming / filtering -----
CUSTOM_DB_NAME  <- "TNK_custom"
CUSTOM_DB_LABEL <- "T/NK custom markers"
CUSTOM_ENRICHMENT_SIZE_RULE <- c(min = 2L, max = 120L)
CELLMARKER_CELLTYPE_PATTERN <- paste(
  "T cell|T-cell|CD4 T|CD8 T|Helper T|Cytotoxic T|Regulatory T|Treg|",
  "Th1|Th17|Tfh|MAIT|gamma delta|γδ T|Innate lymphoid|ILC|NK cell|Natural killer"
)
PANGLAODB_CELLTYPE_PATTERN <- paste(
  "T cell|CD4 T|CD8 T|cytotoxic T|helper T|regulatory T|Treg|",
  "Th1|Th17|Tfh|MAIT|gamma delta|NK cell|natural killer|ILC"
)

# ----- Known markers for visualization -----
KNOWN_MARKERS <- unique(c(
  # Pan-T
  "CD3D", "CD3E", "CD247", "TRAC", "IL7R",
  # CD4 naive / memory
  "CCR7", "SELL", "TCF7", "LEF1", "LTB",
  # CD4 effector / helper
  "CXCR5", "PDCD1", "ICOS", "BCL6", "MAF",
  "TBX21", "IFNG", "CXCR3", "RORC", "IL17A", "CCR6",
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
  # CD8 / cytotoxic
  "CD8A", "CD8B", "KLRG1", "FGFBP2", "PRF1", "GZMB", "GZMK",
  "NKG7", "GNLY", "CX3CR1",
  # TRM / tissue residency
  "CD69", "ITGAE", "CXCR6", "ZNF683", "XCL1",
  # gamma-delta / MAIT / innate-like
  "TRDC", "TRGC1", "TRGC2", "TRAV1-2", "KLRB1", "SLC4A10", "ZBTB16",
  # NK
  "FCGR3A", "KLRD1", "NCR1", "XCL2", "CCL3", "CCL4",
  # exhaustion / dysfunction
  "HAVCR2", "LAG3", "TOX", "LAYN",
  # ILC3
  "KIT", "AHR", "IL23R", "NCR2"
))

# ----- Custom T/NK marker database -----
CUSTOM_MARKERS_DB <- data.frame(
  subtype = c(
    "CD4 Naive/TCM",
    "CD4 Tcm",
    "CD4 Tfh",
    "CD4 Tfr",
    "CD4 Th1",
    "CD4 Th17",
    "CD4 Treg",
    "CD4 Trm",
    "CD8 Naive",
    "CD8 Teff",
    "CD8 Tem",
    "CD8 Temra",
    "CD8 Trm",
    "ILC3",
    "MAIT",
    "NK",
    "NK Exhausted",
    "gdT"
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
)

# ----- Tissue-pair-specific context for LLM -----
LINEAGE_BASE_CONTEXT <- paste(
  "T cells, NK cells, and innate lymphoid cells from NORMAL (non-diseased) human respiratory tract tissues.",
  "This is a cross-site anatomical comparison, NOT a disease vs healthy comparison.",
  "Key lymphocyte programs include CD4 helper / Treg states, CD8 cytotoxic / resident memory states,",
  "innate-like T cells (MAIT, gamma-delta), NK cell effector states, and ILC3 programs.",
  "Focus: regional variation in tissue residency, effector polarization, cytotoxic surveillance,",
  "and mucosal immune adaptation along the respiratory tract."
)

TISSUE_CONTEXT <- list(
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
)

# ----- TNK changelog text -----
PIPELINE_CHANGELOG_TITLE <- sprintf("%s Changes:", PIPELINE_VERSION_LABEL)
PIPELINE_CHANGELOG_LINES <- c(
  "  [TNK-1] Uses `scanvi_label_refined` as the standardized L3 source column.",
  "  [TNK-2] Rebuilds deterministic L2 labels from refined L3 states instead of reusing inconsistent input `cell_type_L2`.",
  "  [TNK-3] Adds TNK-specific marker panels, enrichment filters, and custom marker database.",
  "  [TNK-4] Syncs to the shared v2.6 engine with pseudobulk cache + validated gene->pathway evidence.",
  "  [TNK-5] Enables grouped ssGSEA multi-database LLM outputs (`LLM_SSGSEA_INTERPRETATION.md`).",
  "  [TNK-6] Enables CHOIR per-cluster LLM outputs (`LLM_CHOIR_INTERPRETATION.md`).",
  "  [TNK-7] Keeps the TNK wrapper lightweight by overriding lineage-specific configuration only."
)
PIPELINE_INHERITED_FIXES_TITLE <- "Shared engine fixes inherited from the validated template:"
PIPELINE_INHERITED_FIX_LINES <- c(
  "  [ENGINE-1] Sparse-friendly pseudobulk aggregation with one-shot cache per level.",
  "  [ENGINE-2] Unmapped L3 labels fail fast instead of being silently dropped.",
  "  [ENGINE-3] interpret_agent retry + malformed-output normalization to fixed JSON schema.",
  "  [ENGINE-4] UMAP/heatmap export tuned for publication-friendly PDF+PNG outputs.",
  "  [ENGINE-5] Exploratory wilcox outputs explicitly labeled as cell-level exploratory evidence.",
  "  [ENGINE-6] ssGSEA and CHOIR layers support dedicated LLM interpretation outputs."
)

# ----- TNK-specific LLM interpretation tuning -----
ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
if (file.exists(ADVANCED_HELPER_PATH)) source(ADVANCED_HELPER_PATH)
if (exists("tc_apply_advanced_shared_overrides")) {
  tc_apply_advanced_shared_overrides(overrides = list(
    INTERPRET_MULTI_DB_MIN_TERMS = 1L,
    LLM_TOP_DEG_N = 10L,
    LLM_DEG_PADJ_THR = 0.10,
    LLM_DEG_LFC_THR = 0.15,
    LLM_SSGSEA_TERMS_PER_DIRECTION = 6L,
    LLM_EXTRA_RULES = c(
      "Start overview with a concise cell-type/state judgment.",
      "Use ssGSEA mainly to support cell-type judgment; avoid over-interpreting pathways as direct mechanisms.",
      "For CHOIR clusters, integrate up/down DEG and cluster-level ssGSEA into one compact interpretation."
    )
  ))
}

source("/home/h2048/script/R/bcell_tissue_comparison_v2_6_20260406.R")
