#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Tissue Comparison Pipeline - v2.4.0-TNK
# ==============================================================================
#
# Purpose:
#   T/NK lineage tissue comparison across normal respiratory tract sites.
#   This script reuses the validated B-cell tissue-comparison backbone, while
#   overriding lineage-specific configuration for TNK data:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. Use existing cell_type_L2 and refined scanvi_label_refined as L3
#     3. Visualization (UMAP, marker dotplot/heatmap, composition)
#     4. Pseudobulk DESeq2 per L2/L3 subtype x tissue pair
#     5. Exploratory cell-level wilcox per L2/L3 subtype
#     6. Multi-database enrichment + interpret_agent
#     7. Pseudobulk ssGSEA
#     8. CHOIR + OFA
#     9. Structured REPORT.md + RDS/h5ad export
#
# Notes:
#   - L2 is kept from the existing `cell_type_L2` column.
#   - L3 is overridden to `scanvi_label_refined` (18 refined TNK states).
#
# Smoke test:
#   PIPELINE_TEST_MODE=true Rscript tnk_tissue_comparison_v2_4_20260401.R
#   -> exits after load + metadata validation + NormalizeData()
#
# ==============================================================================

LINEAGE_TAG           <- "TNK"
LINEAGE_DISPLAY       <- "T/NK"
LINEAGE_CONTEXT_LABEL <- "T cells, NK cells, and innate lymphoid cells"
LINEAGE_CONTEXT_LOWER <- "T/NK-cell"
BIOLOGICAL_QUESTION_FRAGMENT <- paste(
  "T/NK cell residency, helper-vs-cytotoxic polarization, innate-like lymphocyte programs,",
  "or tissue microenvironment"
)
PIPELINE_VERSION_LABEL <- "v2.4.0-TNK"
PIPELINE_SUBTITLE      <- "T/NK Tissue Comparison v2.4.0-TNK (pseudobulk DESeq2, multi-database enrichment)"
REPORT_TITLE           <- "# T/NK Tissue Comparison Report (Normal Respiratory Tract)"
GENERATED_BY_LABEL     <- "tnk_tissue_comparison_v2_4_20260401.R"
FINAL_FILE_PREFIX      <- "tnk_tissue_comparison_final"
LINEAGE_COMPLETION_BANNER <- "T/NK TISSUE COMPARISON COMPLETE (v2.4.0-TNK)"

# ----- Input / Output -----
H5AD_PATH  <- "/home/h2048/data/py/0318/tnk_subcluster_retrain/adata_tnk_scanvi_ref_retrain_v1_2.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0401/tnk_tissue_comparison_v2_4_0"

# ----- Label columns -----
L3_SOURCE_COL <- "scanvi_label_refined"
L2_SOURCE_COL <- "cell_type_L2"
USE_EXISTING_L2 <- TRUE

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

source("/home/h2048/script/R/bcell_tissue_comparison_v2_4_20260390.R")
