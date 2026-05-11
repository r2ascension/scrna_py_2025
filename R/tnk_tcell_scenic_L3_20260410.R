#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# T/NK -> T-cell L3 SCENIC wrapper
# ==============================================================================
# Purpose:
#   1. Load the finalized T/NK tissue-comparison Seurat object
#   2. Keep T-cell compartments only (CD4/CD8 axis; includes MAIT/gdT via L2)
#   3. Create a SCENIC-specific L3 column with naive/tcm collapsed to naive
#   4. Reuse scenic_core_20260410.R for regulon inference and AUCell scoring
#
# Notes:
#   - Original `cell_type_L3` is preserved as `cell_type_L3_original`
#   - SCENIC grouping uses `cell_type_L3_scenic`
#   - `CD4 Naive/TCM` and `CD4 Tcm` are merged into `CD4 Naive`
#   - Generic support is also kept for future `CD8 Tcm` / `CD8 Naive/TCM`
#
# Preflight only:
#   SCENIC_PREFLIGHT_ONLY=true Rscript tnk_tcell_scenic_L3_20260410.R
#   -> load object, subset T cells, write label summary, then exit
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS  = "1"
)

INPUT_RDS  <- "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.rds"
OUTPUT_DIR <- "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/SCENIC_tcell_L3_naive_merged"
DATABASE_DIR <- "/home/h2048/data/index_genome/cisTarget_databases_rscenic"
SCENIC_DB_10KB <- file.path(DATABASE_DIR, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
SCENIC_DB_500BP <- file.path(DATABASE_DIR, "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
ORGANISM   <- "hgnc"
ASSAY_USE  <- "RNA"
CELL_TYPE_COL <- "cell_type_L3_scenic"
SAMPLE_COL <- "sample"
N_CORES    <- 8
SCENIC_DATASET_TITLE <- "TNK_Tcell_L3_SCENIC"
REDUCTION_CANDIDATES <- c("umap_refined", "umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap", "harmony", "pca")

TCELL_L2_KEEP <- c("CD4 T cells", "CD8 T cells")
SOURCE_L3_COL <- "cell_type_L3"
TARGET_L3_COL <- "cell_type_L3_scenic"
NAIVE_REMAP <- c(
  "CD4 Naive/TCM" = "CD4 Naive",
  "CD4 Tcm" = "CD4 Naive",
  "CD8 Naive/TCM" = "CD8 Naive",
  "CD8 Tcm" = "CD8 Naive"
)

SCENIC_PREFLIGHT_ONLY <- tolower(Sys.getenv("SCENIC_PREFLIGHT_ONLY", "false")) %in% c("1", "true", "yes")
SCENIC_SOURCE_ONLY <- TRUE
STATUS_LOG <- file.path(OUTPUT_DIR, "wrapper_status.log")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

status_msg <- function(...) {
  txt <- sprintf(...)
  cat(txt)
  cat(txt, file = STATUS_LOG, append = TRUE)
  try(flush(stdout()), silent = TRUE)
  try(flush.console(), silent = TRUE)
  invisible(txt)
}

CORE_SCRIPT <- "/home/h2048/script/R/scenic_core_20260410.R"
status_msg("[BOOT] wrapper started at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
status_msg("[BOOT] SCENIC_PREFLIGHT_ONLY=%s\n", SCENIC_PREFLIGHT_ONLY)
status_msg("[BOOT] CORE_SCRIPT=%s\n", CORE_SCRIPT)
source(CORE_SCRIPT)
status_msg("[BOOT] sourced scenic_core successfully\n")

normalize_tcell_l3_for_scenic <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  matched <- x %in% names(NAIVE_REMAP)
  x[matched] <- unname(NAIVE_REMAP[x[matched]])

  x <- trimws(x)
  x[nzchar(x) == FALSE] <- NA_character_
  x
}

prepare_tcell_object_for_scenic <- function(obj) {
  stopifnot(inherits(obj, "Seurat"))

  meta <- obj@meta.data
  required_cols <- c("cell_type_L2", SOURCE_L3_COL, SAMPLE_COL)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required metadata columns: %s", paste(missing_cols, collapse = ", ")))
  }

  keep_cells <- rownames(meta)[as.character(meta$cell_type_L2) %in% TCELL_L2_KEEP]
  if (length(keep_cells) == 0) {
    stop("No T-cell cells found after filtering by cell_type_L2.")
  }

  obj_t <- subset(obj, cells = keep_cells)
  obj_t@meta.data[["cell_type_L3_original"]] <- as.character(obj_t@meta.data[[SOURCE_L3_COL]])
  obj_t@meta.data[[TARGET_L3_COL]] <- normalize_tcell_l3_for_scenic(obj_t@meta.data[[SOURCE_L3_COL]])

  bad_cells <- rownames(obj_t@meta.data)[is.na(obj_t@meta.data[[TARGET_L3_COL]])]
  if (length(bad_cells) > 0) {
    stop(sprintf("Found %d cells with NA %s after remapping.", length(bad_cells), TARGET_L3_COL))
  }

  obj_t
}

write_tcell_l3_summary <- function(obj_t, output_dir = OUTPUT_DIR) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

  summary_dt <- data.table::data.table(
    cell = colnames(obj_t),
    cell_type_L2 = as.character(obj_t@meta.data[["cell_type_L2"]]),
    cell_type_L3_original = as.character(obj_t@meta.data[["cell_type_L3_original"]]),
    cell_type_L3_scenic = as.character(obj_t@meta.data[[TARGET_L3_COL]]),
    sample = as.character(obj_t@meta.data[[SAMPLE_COL]])
  )
  data.table::fwrite(summary_dt, file.path(output_dir, "tables", "tcell_l3_scenic_cell_metadata.csv"))

  remap_dt <- summary_dt[, .N, by = .(cell_type_L3_original, cell_type_L3_scenic)][order(cell_type_L3_original, -N)]
  data.table::fwrite(remap_dt, file.path(output_dir, "tables", "tcell_l3_scenic_remap_summary.csv"))

  invisible(remap_dt)
}

main <- function() {
  status_msg("[LOAD] %s\n", INPUT_RDS)
  obj_input <- readRDS(INPUT_RDS)
  status_msg("[INFO] Input object loaded: %d cells\n", ncol(obj_input))

  obj_tcell <- prepare_tcell_object_for_scenic(obj_input)
  status_msg("[INFO] T-cell subset prepared: %d cells\n", ncol(obj_tcell))

  remap_dt <- write_tcell_l3_summary(obj_tcell)
  status_msg("[INFO] SCENIC L3 groups: %s\n", paste(sort(unique(obj_tcell@meta.data[[TARGET_L3_COL]])), collapse = ", "))
  status_msg("[INFO] Remap rows written: %d\n", nrow(remap_dt))

  if (SCENIC_PREFLIGHT_ONLY) {
    status_msg("[DONE] SCENIC_PREFLIGHT_ONLY=true; preprocessing finished without running SCENIC.\n")
    return(invisible(NULL))
  }

  status_msg("[RUN] Starting run_scenic_module()\n")
  result <- run_scenic_module(obj_tcell)
  status_msg("[DONE] run_scenic_module() completed successfully\n")
  invisible(result)
}

tryCatch(
  main(),
  error = function(e) {
    status_msg("[ERROR] %s\n", conditionMessage(e))
    quit(save = "no", status = 1)
  }
)
