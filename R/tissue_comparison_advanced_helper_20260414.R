#!/usr/bin/env Rscript
# ==============================================================================
# Tissue Comparison Advanced Helper Extensions
# ==============================================================================
#
# Purpose:
#   - keep the 2026-04-08 helper immutable
#   - layer stromal lineage marker-panel defaults on top of the validated helper
#   - provide a new versioned helper path for wrapper-based stromal pipelines
#
# Date: 2026-04-14
# ==============================================================================

BASE_ADVANCED_HELPER_PATH <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
if (!file.exists(BASE_ADVANCED_HELPER_PATH)) {
  stop(sprintf("Base advanced helper not found: %s", BASE_ADVANCED_HELPER_PATH))
}
source(BASE_ADVANCED_HELPER_PATH)

tc_panel_genes_with_custom <- function(base_genes,
                                       custom_markers_db = NULL,
                                       subtype_pattern = NULL) {
  base_genes <- toupper(trimws(as.character(base_genes)))
  base_genes <- base_genes[!is.na(base_genes) & nzchar(base_genes)]
  unique(c(
    base_genes,
    tc_extract_custom_panel_markers(custom_markers_db, subtype_pattern = subtype_pattern)
  ))
}

tc_default_stromal_endothelial_marker_panels <- function(known_markers = NULL,
                                                         custom_markers_db = NULL) {
  list(
    pan_endothelial = list(
      label = "Pan-endothelial",
      genes = tc_panel_genes_with_custom(
        c("PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EMCN", "PLVAP"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "ENDOTHE"
      )
    ),
    lymphatic = list(
      label = "Lymphatic",
      genes = tc_panel_genes_with_custom(
        c("PROX1", "PDPN", "LYVE1", "CCL21", "FLT4", "MMRN1"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "LYMPH"
      )
    ),
    capillary = list(
      label = "Capillary",
      genes = tc_panel_genes_with_custom(
        c("CA4", "EDNRB", "RGCC", "GPIHBP1", "BTNL9", "AQP1", "ADGRL4", "EPAS1"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "CAP|CAPILL"
      )
    ),
    arterial = list(
      label = "Arterial",
      genes = tc_panel_genes_with_custom(
        c("GJA5", "EFNB2", "SOX17", "SEMA3G", "BMX", "HEY1", "CXCL12", "FBLN5"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "ARTERIAL"
      )
    ),
    venous_inflammatory = list(
      label = "Venous / activated",
      genes = tc_panel_genes_with_custom(
        c("ACKR1", "SELE", "SELP", "VCAM1", "NR2F2", "PLVAP", "CXCR4", "ISG15", "IFIT3"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "VENOUS|INFLAM|ACKR1"
      )
    )
  )
}

tc_default_stromal_fibroblast_marker_panels <- function(known_markers = NULL,
                                                        custom_markers_db = NULL) {
  list(
    pan_fibroblast = list(
      label = "Pan-fibroblast",
      genes = tc_panel_genes_with_custom(
        c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "MFAP4", "PDGFRA"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "FIBRO"
      )
    ),
    adventitial = list(
      label = "Adventitial",
      genes = tc_panel_genes_with_custom(
        c("PI16", "MFAP5", "CD34", "C7", "DPT", "TCF21"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "ADVENTITIAL"
      )
    ),
    airway_matrix = list(
      label = "Airway matrix / peribronchial",
      genes = tc_panel_genes_with_custom(
        c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "MFAP4", "FN1", "POSTN"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "PERIBRONCHIAL|AIRWAY"
      )
    ),
    alveolar = list(
      label = "Alveolar-support fibroblast",
      genes = tc_panel_genes_with_custom(
        c("APOE", "FABP4", "INMT", "NPNT", "CFD", "CXCL14"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "ALVEOLAR"
      )
    ),
    myofibroblast = list(
      label = "Myofibroblast / activated ECM",
      genes = tc_panel_genes_with_custom(
        c("POSTN", "CTHRC1", "ACTA2", "TAGLN", "COL11A1", "FN1"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "MYOFIBROBLAST|ECM|ACTIVAT"
      )
    ),
    stress_interferon = list(
      label = "Stress / interferon",
      genes = tc_panel_genes_with_custom(
        c("ISG15", "IFIT1", "IFIT3", "CXCL10", "STAT1", "BST2"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "STRESS|IFN|INTERFERON"
      )
    )
  )
}

tc_default_stromal_smc_marker_panels <- function(known_markers = NULL,
                                                 custom_markers_db = NULL) {
  list(
    pan_mural = list(
      label = "Pan-mural",
      genes = tc_panel_genes_with_custom(
        c("RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2", "TAGLN", "MYH11"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "MUSCLE|PERICYTE|MURAL"
      )
    ),
    pericyte = list(
      label = "Pericyte",
      genes = tc_panel_genes_with_custom(
        c("RGS5", "PDGFRB", "CSPG4", "MCAM", "ABCC9", "KCNJ8", "NOTCH3", "DES"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "PERICYTE"
      )
    ),
    smooth_muscle = list(
      label = "Contractile smooth muscle",
      genes = tc_panel_genes_with_custom(
        c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "PRKG1", "CALD1", "DES"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "SMOOTH"
      )
    ),
    immune_recruiting = list(
      label = "Perivascular immune-recruiting",
      genes = tc_panel_genes_with_custom(
        c("CXCL12", "CCL2", "SEMA3G", "RGS5", "PDGFRB", "NOTCH3"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "IMMUNE|RECRUIT|CXCL12|CCL2"
      )
    ),
    pulmonary_systemic_specialization = list(
      label = "Pulmonary / systemic specialization",
      genes = tc_panel_genes_with_custom(
        c("RGS5", "PDGFRB", "KCNJ8", "ABCC9", "ACTA2", "MYH11", "SEMA3G", "CXCL12"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "PULMONARY|SYSTEMIC"
      )
    )
  )
}

if (sys.nframe() == 0) {
  cat("Tissue Comparison Advanced Helper extensions (2026-04-14) loaded.\n")
  cat(sprintf("Base helper: %s\n", BASE_ADVANCED_HELPER_PATH))
}