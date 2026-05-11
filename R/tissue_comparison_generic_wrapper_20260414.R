#!/usr/bin/env Rscript
# ==============================================================================
# Generic Tissue Comparison Wrapper Interface (Stromal Extension)
# ==============================================================================
#
# Purpose:
#   - keep the validated 2026-04-12 wrapper immutable
#   - extend preset coverage to stromal endothelial / fibroblast / SMC branches
#   - point stromal wrappers at a versioned helper extension without changing
#     the base helper or base wrapper
#
# Date: 2026-04-14
# ==============================================================================

BASE_GENERIC_WRAPPER_PATH <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
if (!file.exists(BASE_GENERIC_WRAPPER_PATH)) {
  stop(sprintf("Base generic wrapper not found: %s", BASE_GENERIC_WRAPPER_PATH))
}
source(BASE_GENERIC_WRAPPER_PATH)

ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260414.R"
if (!file.exists(ADVANCED_HELPER_PATH)) {
  stop(sprintf("Extended advanced helper not found: %s", ADVANCED_HELPER_PATH))
}
source(ADVANCED_HELPER_PATH)

tc_get_lineage_preset_base_20260412 <- tc_get_lineage_preset

tc_stromal_shared_overrides <- function(state_judgment_label,
                                        extra_rules = character()) {
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
      sprintf("Start overview with a concise %s state judgment.", state_judgment_label),
      "Use ssGSEA mainly to support cell-state judgment; avoid converting broad pathway shifts into unsupported direct mechanisms.",
      "Prefer lineage-informative stromal programs over broad stress, mitochondrial, ribosomal, interferon-only, or leukocyte-contamination summaries when ranking discoveries.",
      "For CHOIR clusters, integrate annotation composition, up/down DEG, ssGSEA, and OFA evidence into one compact interpretation.",
      extra_rules
    )
  )
}

tc_stromal_endothelial_wrapper_preset <- function() {
  known_markers <- unique(c(
    "PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EMCN",
    "PROX1", "PDPN", "LYVE1", "CCL21", "FLT4", "MMRN1",
    "CA4", "EDNRB", "RGCC", "GPIHBP1", "BTNL9", "AQP1",
    "GJA5", "EFNB2", "SOX17", "SEMA3G", "BMX",
    "ACKR1", "SELE", "SELP", "VCAM1", "CXCL12", "CXCR4",
    "ISG15", "IFIT1", "IFIT3", "STAT1",
    "PTPRC", "LST1"
  ))
  custom_markers_db <- data.frame(
    subtype = c(
      "Endothelia_Lymphatic",
      "Endothelia_vascular_Cap_a",
      "Endothelia_vascular_Cap_g",
      "Endothelia_vascular_arterial_pulmonary",
      "Endothelia_vascular_arterial_systemic",
      "Endothelia_vascular_venous_pulmonary",
      "Endothelia_vascular_venous_systemic"
    ),
    markers = c(
      "PROX1,PDPN,LYVE1,CCL21,FLT4,MMRN1",
      "CA4,EDNRB,RGCC,GPIHBP1,CD36,AQP1",
      "KDR,EMCN,BTNL9,ADGRL4,EPAS1,PLVAP",
      "GJA5,EFNB2,SOX17,CXCL12,KCNJ8,FBLN5",
      "GJA5,EFNB2,SEMA3G,BMX,HEY1,SOX17",
      "ACKR1,SELP,VWF,NR2F2,PLVAP,EMCN",
      "ACKR1,SELE,VWF,VCAM1,CXCR4,ISG15"
    ),
    stringsAsFactors = FALSE
  )
  marker_panels <- tc_default_stromal_endothelial_marker_panels(
    known_markers = known_markers,
    custom_markers_db = custom_markers_db
  )

  list(
    LINEAGE_TAG = "STROMAL_ENDOTHELIAL",
    LINEAGE_DISPLAY = "Stromal Endothelial",
    LINEAGE_CONTEXT_LABEL = "stromal endothelial cells",
    LINEAGE_CONTEXT_LOWER = "stromal-endothelial",
    BIOLOGICAL_QUESTION_FRAGMENT = paste(
      "vascular zonation, lymphatic drainage, capillary specialization, arterial-versus-venous programs,",
      "barrier adaptation, endothelial activation, or tissue microenvironment"
    ),
    PIPELINE_VERSION_LABEL = "v1.1.0-Endothelial",
    PIPELINE_SUBTITLE = paste(
      "Stromal Endothelial Tissue Comparison v1.1.0",
      "(generic wrapper interface, branch-specific scANVI reference)"
    ),
    REPORT_TITLE = "# Stromal Endothelial Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "stromal_endothelial_tissue_comparison_v1_1_0_20260414.R",
    FINAL_FILE_PREFIX = "stromal_endothelial_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "STROMAL ENDOTHELIAL TISSUE COMPARISON COMPLETE (v1.1.0-Endothelial)",
    H5AD_PATH = "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/endothelial/adata_endothelial_reference_v1_5_branchwise.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0408/stromal_endothelial_tissue_comparison_v1_0",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    L2_SOURCE_COL = "cell_type_L2",
    USE_EXISTING_L2 = TRUE,
    ANALYSIS_L3_DESCRIPTION = "branch-specific endothelial scANVI predictions",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    UMAP_REDUCTION_PREFERRED = c("umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap"),
    CHOIR_REDUCTION_CANDIDATES = c("scanvi", "scvi", "pca"),
    CHOIR_ALPHA = 0.20,
    CHOIR_VAR_FEATURES_MAX = 4000L,
    PSEUDOBULK_DISAMBIGUATION_CANDIDATES = c("dataset", "source", "batch", "orig.ident"),
    CUSTOM_DB_NAME = "Stromal_Endothelial_custom",
    CUSTOM_DB_LABEL = "stromal endothelial custom markers",
    CUSTOM_ENRICHMENT_SIZE_RULE = c(min = 2L, max = 120L),
    CELLMARKER_CELLTYPE_PATTERN = paste(
      "Endothelial|Endothelia|vascular endothelial|capillary endothelial|arterial endothelial|",
      "venous endothelial|lymphatic endothelial|microvascular endothelial"
    ),
    PANGLAODB_CELLTYPE_PATTERN = paste(
      "Endothelial|vascular endothelial|capillary endothelial|arterial endothelial|",
      "venous endothelial|lymphatic endothelial"
    ),
    KNOWN_MARKERS = known_markers,
    CUSTOM_MARKERS_DB = custom_markers_db,
    MARKER_PANELS = marker_panels,
    LINEAGE_BASE_CONTEXT = paste(
      "Stromal endothelial cells from NORMAL (non-diseased) human respiratory tract tissues.",
      "This is a cross-site anatomical comparison, NOT a disease-vs-healthy comparison.",
      "Key endothelial programs include lymphatic drainage, capillary specialization, arterial/venous zonation,",
      "vascular barrier maintenance, leukocyte trafficking, and inflammatory activation states.",
      "Focus: regional variation in endothelial specialization along nose, sinus, conducting airway, and lung parenchyma."
    ),
    TISSUE_CONTEXT = list(
      "nose" = paste(
        "Nasal mucosa: first-line vascular interface under intense environmental exposure,",
        "with active leukocyte trafficking, vascular barrier tuning, and lymphatic drainage."
      ),
      "sinus" = paste(
        "Paranasal sinus: semi-enclosed mucosal compartment where endothelial programs may reflect drainage,",
        "fluid handling, and low-flow barrier surveillance."
      ),
      "respiratory airway" = paste(
        "Conducting airway: endothelial cells support mucosal immune cell recruitment, epithelial-stromal crosstalk,",
        "and vessel adaptation near bronchi and glands."
      ),
      "lung parenchyma" = paste(
        "Lung parenchyma: endothelial specialization is tightly linked to alveolar gas exchange,",
        "capillary zonation, vascular permeability control, and pulmonary circulation biology."
      )
    ),
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    SHARED_OVERRIDES = tc_stromal_shared_overrides(
      state_judgment_label = "endothelial subtype",
      extra_rules = c(
        "Prioritize lymphatic, capillary, arterial, venous, barrier, and leukocyte-trafficking programs over fibroblast or immune contamination.",
        "Endothelial interferon/activation signals should be described as state modifiers, not as replacement identities, unless the broader endothelial program is lost."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] Stromal endothelial entrypoint is now a thin preset over a versioned stromal wrapper extension.",
      "  [ENDO-1] Force fresh h5ad load from the branchwise endothelial reference; do not reuse the 0408 final object.",
      "  [ENDO-2] Add endothelial-focused marker panels while keeping existing L2 and branch-specific scANVI L3 labels.",
      "  [ENDO-3] New 0414 output directory uses 0408 endothelial results as historical baseline only.",
      "  [ENDO-4] Wrapper requires a valid DEEPSEEK key, so the full rerun cannot silently skip LLM steps."
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

tc_stromal_fibroblast_wrapper_preset <- function() {
  known_markers <- unique(c(
    "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "MFAP4", "PDGFRA",
    "PI16", "MFAP5", "CD34", "TCF21", "C7", "DPT",
    "POSTN", "CTHRC1", "FN1", "ACTA2", "TAGLN", "COL11A1",
    "APOE", "FABP4", "INMT", "CFD", "CXCL14",
    "ISG15", "IFIT1", "IFIT3", "CXCL10", "STAT1", "BST2",
    "PTPRC", "MS4A1"
  ))
  custom_markers_db <- data.frame(
    subtype = c(
      "Fibro_adventitial",
      "Fibro_peribronchial",
      "Fibro_alveolar",
      "Fibro_stress_activated",
      "Fibro_myofibroblast"
    ),
    markers = c(
      "PI16,MFAP5,CD34,C7,DPT,TCF21",
      "COL1A1,COL1A2,COL3A1,LUM,DCN,MFAP4",
      "APOE,FABP4,INMT,NPNT,CFD,CXCL14",
      "ISG15,IFIT1,IFIT3,CXCL10,STAT1,BST2",
      "POSTN,CTHRC1,ACTA2,TAGLN,COL11A1,FN1"
    ),
    stringsAsFactors = FALSE
  )
  marker_panels <- tc_default_stromal_fibroblast_marker_panels(
    known_markers = known_markers,
    custom_markers_db = custom_markers_db
  )

  list(
    LINEAGE_TAG = "STROMAL_FIBROBLAST",
    LINEAGE_DISPLAY = "Stromal Fibroblast",
    LINEAGE_CONTEXT_LABEL = "stromal fibroblasts",
    LINEAGE_CONTEXT_LOWER = "stromal-fibroblast",
    BIOLOGICAL_QUESTION_FRAGMENT = paste(
      "ECM remodeling, niche support, adventitial-versus-alveolar specialization,",
      "peribronchial support programs, stress activation, wound-response, or tissue microenvironment"
    ),
    PIPELINE_VERSION_LABEL = "v1.1.0-Fibroblast",
    PIPELINE_SUBTITLE = paste(
      "Stromal Fibroblast Tissue Comparison v1.1.0",
      "(generic wrapper interface, branch-specific scANVI reference)"
    ),
    REPORT_TITLE = "# Stromal Fibroblast Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "stromal_fibroblast_tissue_comparison_v1_1_0_20260414.R",
    FINAL_FILE_PREFIX = "stromal_fibroblast_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "STROMAL FIBROBLAST TISSUE COMPARISON COMPLETE (v1.1.0-Fibroblast)",
    H5AD_PATH = "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/fibroblast/adata_fibroblast_reference_v1_5_branchwise.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_0_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0408/stromal_fibroblast_tissue_comparison_v1_0",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    L2_SOURCE_COL = "cell_type_L2",
    USE_EXISTING_L2 = TRUE,
    ANALYSIS_L3_DESCRIPTION = "branch-specific fibroblast scANVI predictions",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    UMAP_REDUCTION_PREFERRED = c("umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap"),
    CHOIR_REDUCTION_CANDIDATES = c("scanvi", "scvi", "pca"),
    CHOIR_ALPHA = 0.20,
    CHOIR_VAR_FEATURES_MAX = 4000L,
    PSEUDOBULK_DISAMBIGUATION_CANDIDATES = c("dataset", "source", "batch", "orig.ident"),
    CUSTOM_DB_NAME = "Stromal_Fibroblast_custom",
    CUSTOM_DB_LABEL = "stromal fibroblast custom markers",
    CUSTOM_ENRICHMENT_SIZE_RULE = c(min = 2L, max = 120L),
    CELLMARKER_CELLTYPE_PATTERN = paste(
      "Fibroblast|stromal fibroblast|adventitial fibroblast|alveolar fibroblast|",
      "peribronchial fibroblast|myofibroblast|mesenchymal stromal"
    ),
    PANGLAODB_CELLTYPE_PATTERN = paste(
      "Fibroblast|stromal fibroblast|adventitial fibroblast|alveolar fibroblast|",
      "peribronchial fibroblast|myofibroblast"
    ),
    KNOWN_MARKERS = known_markers,
    CUSTOM_MARKERS_DB = custom_markers_db,
    MARKER_PANELS = marker_panels,
    LINEAGE_BASE_CONTEXT = paste(
      "Stromal fibroblasts from NORMAL (non-diseased) human respiratory tract tissues.",
      "This is a cross-site anatomical comparison, NOT a disease-vs-healthy comparison.",
      "Key fibroblast programs include ECM deposition, matrix remodeling, airway wall support,",
      "adventitial niche formation, alveolar stromal support, and stress / repair activation.",
      "Focus: regional variation in fibroblast specialization across upper airway, conducting airway, and parenchyma."
    ),
    TISSUE_CONTEXT = list(
      "nose" = paste(
        "Nasal mucosa: fibroblasts contribute to barrier support, lamina propria structure,",
        "and rapid remodeling under constant environmental exposure."
      ),
      "sinus" = paste(
        "Paranasal sinus: fibroblast programs may reflect matrix maintenance in a semi-enclosed cavity,",
        "with drainage-dependent remodeling and stromal support."
      ),
      "respiratory airway" = paste(
        "Conducting airway: fibroblasts support peribronchial ECM, airway wall mechanics, gland-associated stroma,",
        "and epithelial-stromal signaling."
      ),
      "lung parenchyma" = paste(
        "Lung parenchyma: fibroblasts support alveolar architecture, interstitial matrix balance,",
        "and niche signaling tied to gas-exchange tissue homeostasis."
      )
    ),
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    SHARED_OVERRIDES = tc_stromal_shared_overrides(
      state_judgment_label = "fibroblast subtype",
      extra_rules = c(
        "Prioritize adventitial, alveolar-support, peribronchial matrix, ECM-remodeling, and myofibroblast programs over endothelial or immune contamination.",
        "Stress or interferon signals should be interpreted as fibroblast state modifiers unless the broader fibroblast matrix program is lost."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] Stromal fibroblast entrypoint is now a thin preset over a versioned stromal wrapper extension.",
      "  [FIB-1] Force fresh h5ad load from the branchwise fibroblast reference; do not reuse the 0408 final object.",
      "  [FIB-2] Add fibroblast-focused marker panels while keeping existing L2 and branch-specific scANVI L3 labels.",
      "  [FIB-3] New 0414 output directory uses 0408 fibroblast results as historical baseline only.",
      "  [FIB-4] Wrapper requires a valid DEEPSEEK key, so the full rerun cannot silently skip LLM steps."
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

tc_stromal_smc_wrapper_preset <- function() {
  known_markers <- unique(c(
    "RGS5", "PDGFRB", "CSPG4", "MCAM", "ABCC9", "KCNJ8", "NOTCH3",
    "ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "PRKG1", "CALD1", "DES",
    "CXCL12", "CCL2", "SEMA3G", "CASQ2",
    "SOX10", "S100B",
    "PTPRC", "LST1"
  ))
  custom_markers_db <- data.frame(
    subtype = c(
      "Muscle_pericyte_pulmonary",
      "Muscle_pericyte_systemic",
      "Muscle_perivascular_immune_recruiting",
      "Muscle_smooth_pulmonary",
      "Muscle_smooth_arterial_systemic"
    ),
    markers = c(
      "RGS5,PDGFRB,CSPG4,MCAM,KCNJ8,ABCC9",
      "RGS5,PDGFRB,NOTCH3,MCAM,DES,CSPG4",
      "RGS5,PDGFRB,CXCL12,CCL2,SEMA3G,NOTCH3",
      "ACTA2,TAGLN,MYH11,CNN1,MYLK,PRKG1",
      "ACTA2,MYH11,TAGLN,CALD1,CNN1,MYLK"
    ),
    stringsAsFactors = FALSE
  )
  marker_panels <- tc_default_stromal_smc_marker_panels(
    known_markers = known_markers,
    custom_markers_db = custom_markers_db
  )

  list(
    LINEAGE_TAG = "STROMAL_SMC",
    LINEAGE_DISPLAY = "Stromal SMC/Pericyte",
    LINEAGE_CONTEXT_LABEL = "smooth muscle cells and pericytes",
    LINEAGE_CONTEXT_LOWER = "stromal-smc-pericyte",
    BIOLOGICAL_QUESTION_FRAGMENT = paste(
      "vascular tone, mural-cell specialization, pulmonary-versus-systemic pericyte programs,",
      "contractility, perivascular immune recruitment, or tissue microenvironment"
    ),
    PIPELINE_VERSION_LABEL = "v1.1.0-SMC",
    PIPELINE_SUBTITLE = paste(
      "Stromal SMC/Pericyte Tissue Comparison v1.1.0",
      "(generic wrapper interface, branch-specific scANVI reference)"
    ),
    REPORT_TITLE = "# Stromal SMC / Pericyte Tissue Comparison Report (Normal Respiratory Tract)",
    GENERATED_BY_LABEL = "stromal_smc_tissue_comparison_v1_1_0_20260414.R",
    FINAL_FILE_PREFIX = "stromal_smc_tissue_comparison_final",
    LINEAGE_COMPLETION_BANNER = "STROMAL SMC/PERICYTE TISSUE COMPARISON COMPLETE (v1.1.0-SMC)",
    H5AD_PATH = "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/smc/adata_smc_reference_v1_5_branchwise.h5ad",
    OUTPUT_DIR = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_0_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0408/stromal_smc_tissue_comparison_v1_0",
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE,
    L3_SOURCE_COL = "cell_type_scanvi_pred",
    L2_SOURCE_COL = "cell_type_L2",
    USE_EXISTING_L2 = TRUE,
    ANALYSIS_L3_DESCRIPTION = "branch-specific SMC/pericyte scANVI predictions",
    L3_TO_L2_TABLE_HEADER_LEFT = "L3 (`cell_type_scanvi_pred`)",
    UMAP_REDUCTION_PREFERRED = c("umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap"),
    CHOIR_REDUCTION_CANDIDATES = c("scanvi", "scvi", "pca"),
    CHOIR_ALPHA = 0.20,
    CHOIR_VAR_FEATURES_MAX = 4000L,
    PSEUDOBULK_DISAMBIGUATION_CANDIDATES = c("dataset", "source", "batch", "orig.ident"),
    CUSTOM_DB_NAME = "Stromal_SMC_custom",
    CUSTOM_DB_LABEL = "stromal SMC/pericyte custom markers",
    CUSTOM_ENRICHMENT_SIZE_RULE = c(min = 2L, max = 120L),
    CELLMARKER_CELLTYPE_PATTERN = paste(
      "Smooth muscle|vascular smooth muscle|pericyte|mural cell|perivascular|",
      "vascular mural|perivascular smooth muscle"
    ),
    PANGLAODB_CELLTYPE_PATTERN = paste(
      "Smooth muscle|vascular smooth muscle|pericyte|mural cell|perivascular"
    ),
    KNOWN_MARKERS = known_markers,
    CUSTOM_MARKERS_DB = custom_markers_db,
    MARKER_PANELS = marker_panels,
    LINEAGE_BASE_CONTEXT = paste(
      "Smooth muscle cells and pericytes from NORMAL (non-diseased) human respiratory tract tissues.",
      "This is a cross-site anatomical comparison, NOT a disease-vs-healthy comparison.",
      "Key mural-cell programs include vascular tone control, contractility, vessel stabilization,",
      "pericyte support, and perivascular immune recruitment.",
      "Focus: regional variation in pulmonary versus systemic mural-cell programs across respiratory tissues."
    ),
    TISSUE_CONTEXT = list(
      "nose" = paste(
        "Nasal mucosa: mural cells support vascular reactivity, microvascular stability,",
        "and rapid regulation of leukocyte entry in an exposure-heavy environment."
      ),
      "sinus" = paste(
        "Paranasal sinus: mural-cell programs may reflect drainage-linked vascular tone,",
        "low-flow vessel maintenance, and perivascular immune positioning."
      ),
      "respiratory airway" = paste(
        "Conducting airway: smooth muscle and pericytes contribute to airway-adjacent vessel behavior,",
        "bronchovascular niche support, and contractile / inflammatory crosstalk."
      ),
      "lung parenchyma" = paste(
        "Lung parenchyma: mural-cell specialization is linked to pulmonary circulation,",
        "alveolar microvascular stability, and local control of perfusion and immune trafficking."
      )
    ),
    SHARED_ENGINE_PATH = "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R",
    SHARED_OVERRIDES = tc_stromal_shared_overrides(
      state_judgment_label = "mural-cell subtype",
      extra_rules = c(
        "Prioritize pericyte, contractile smooth-muscle, pulmonary/systemic mural, and immune-recruiting perivascular programs over fibroblast or leukocyte contamination.",
        "Do not let isolated neural-crest-like markers dominate the interpretation if the broader mural program indicates pericyte or smooth-muscle identity."
      )
    ),
    PIPELINE_CHANGELOG_LINES = c(
      "  [ARCH-1] Generic wrapper interface now owns config validation, env loading, and previous-run preflight.",
      "  [ARCH-2] Stromal SMC/pericyte entrypoint is now a thin preset over a versioned stromal wrapper extension.",
      "  [SMC-1] Force fresh h5ad load from the branchwise SMC/pericyte reference; do not reuse the 0408 final object.",
      "  [SMC-2] Add mural-cell-focused marker panels while keeping existing L2 and branch-specific scANVI L3 labels.",
      "  [SMC-3] New 0414 output directory uses 0408 SMC/pericyte results as historical baseline only.",
      "  [SMC-4] Wrapper requires a valid DEEPSEEK key, so the full rerun cannot silently skip LLM steps."
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
  "STROMAL_ENDOTHELIAL", "STROMAL_FIBROBLAST", "STROMAL_SMC"
)) {
  lineage <- toupper(tc_safe_trim(lineage %||% "BCELL"))
  if (identical(lineage, "STROMAL_ENDOTHELIAL")) return(tc_stromal_endothelial_wrapper_preset())
  if (identical(lineage, "STROMAL_FIBROBLAST")) return(tc_stromal_fibroblast_wrapper_preset())
  if (identical(lineage, "STROMAL_SMC")) return(tc_stromal_smc_wrapper_preset())
  tc_get_lineage_preset_base_20260412(lineage)
}

if (sys.nframe() == 0) {
  if (!exists("PIPELINE_CONFIG", inherits = FALSE)) {
    cat("Stromal tissue-comparison wrapper extension loaded.\n")
    cat("Supported extra lineages: STROMAL_ENDOTHELIAL, STROMAL_FIBROBLAST, STROMAL_SMC\n")
  } else {
    tc_run_generic_tissue_comparison(PIPELINE_CONFIG)
  }
}