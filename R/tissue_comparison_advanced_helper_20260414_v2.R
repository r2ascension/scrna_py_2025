#!/usr/bin/env Rscript
# ==============================================================================
# Tissue Comparison Advanced Helper Extensions (Myeloid Extension)
# ==============================================================================
#
# Purpose:
#   - keep the 2026-04-14 helper immutable
#   - layer myeloid lineage marker-panel defaults on top of the stromal helper
#   - provide a new versioned helper path for wrapper-based myeloid pipelines
#
# Date: 2026-04-14
# ==============================================================================

BASE_ADVANCED_HELPER_PATH_V2 <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260414.R"
if (!file.exists(BASE_ADVANCED_HELPER_PATH_V2)) {
  stop(sprintf("Base advanced helper not found: %s", BASE_ADVANCED_HELPER_PATH_V2))
}
source(BASE_ADVANCED_HELPER_PATH_V2)

tc_load_env_file <- function(path, overwrite_placeholder = TRUE) {
  if (!tc_path_exists(path)) return(invisible(FALSE))
  lines <- readLines(path, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- trimws(gsub("^export\\s+", "", sub("=.*$", "", line)))
    val <- trimws(gsub("^['\"]|['\"]$", "", sub("^[^=]*=", "", line)))
    if (!nzchar(key)) next
    cur <- Sys.getenv(key, unset = "")
    should_set <- !nzchar(cur)
    if (!should_set && isTRUE(overwrite_placeholder)) {
      should_set <- tc_is_placeholder_secret(cur)
    }
    if (should_set) {
      do.call(Sys.setenv, stats::setNames(list(val), key))
    }
  }
  invisible(TRUE)
}

tc_default_myeloid_marker_panels <- function(known_markers = NULL,
                                             custom_markers_db = NULL) {
  list(
    pan_myeloid = list(
      label = "Pan-myeloid / APC",
      genes = tc_panel_genes_with_custom(
        c("LYZ", "LST1", "TYMP", "CTSS", "FCER1G", "HLA-DRA", "CST3", "IFI30"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "MACROPH|MONOCYTE|DC|PDC|MAST|NEUTROPHIL|MYELOID"
      )
    ),
    monocyte = list(
      label = "Monocyte",
      genes = tc_panel_genes_with_custom(
        c("FCN1", "VCAN", "S100A8", "S100A9", "CTSD", "SAT1", "FCGR3A", "IFITM3"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "MONOCYTE"
      )
    ),
    alveolar_macrophage = list(
      label = "Alveolar macrophage",
      genes = tc_panel_genes_with_custom(
        c("FABP4", "PPARG", "MARCO", "C1QA", "C1QB", "C1QC", "INHBA", "ABCA1"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "ALVEOLAR"
      )
    ),
    interstitial_macrophage = list(
      label = "Interstitial macrophage",
      genes = tc_panel_genes_with_custom(
        c("APOE", "CD163", "CD163L1", "MRC1", "FOLR2", "MSR1", "C1QC", "CTSB"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "INTERSTITIAL|CD163L1"
      )
    ),
    dendritic_pdc = list(
      label = "Dendritic / pDC",
      genes = tc_panel_genes_with_custom(
        c("FCER1A", "CD1C", "CLEC10A", "HLA-DPA1", "CLEC4C", "LILRA4", "GZMB", "IRF7", "TCF4"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "CDC|PDC|DC|LANGERHANS"
      )
    ),
    mast_cell = list(
      label = "Mast cell",
      genes = tc_panel_genes_with_custom(
        c("KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2", "HPGDS", "HDC"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "MAST"
      )
    ),
    neutrophil = list(
      label = "Neutrophil",
      genes = tc_panel_genes_with_custom(
        c("FCGR3B", "CXCL8", "CSF3R", "NAMPT", "FPR1", "S100A8", "S100A9"),
        custom_markers_db = custom_markers_db,
        subtype_pattern = "NEUTROPHIL"
      )
    )
  )
}

if (sys.nframe() == 0) {
  cat("Tissue Comparison Advanced Helper extensions (2026-04-14 v2) loaded.\n")
  cat(sprintf("Base helper: %s\n", BASE_ADVANCED_HELPER_PATH_V2))
  cat("Added marker panels: MYELOID\n")
}
