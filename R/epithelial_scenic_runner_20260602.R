#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Epithelial SCENIC Runner - Dedicated pipeline for epithelial lineage
# ==============================================================================
# Runs two comparison levels:
#   1) lineage_celltype : Compare all epithelial L3 cell types (one SCENIC run)
#   2) celltype_tissue  : For each eligible L3 cell type, compare across tissues
#
# Epithelial-specific improvements:
#   - Aggressive but stratified downsampling for 278K+ cells
#   - Gene capping to 3000 most informative genes
#   - Tissue-aware regulon specificity scoring
#   - Resume/checkpoint capability
#   - Comprehensive per-celltype status tracking
#
# Usage:
#   # Full run (both levels):
#   Rscript script/R/epithelial_scenic_runner_20260602.R
#
#   # Specific level only:
#   EPI_SCENIC_LEVELS=lineage_celltype Rscript script/R/epithelial_scenic_runner_20260602.R
#   EPI_SCENIC_LEVELS=celltype_tissue Rscript script/R/epithelial_scenic_runner_20260602.R
#
#   # Specific cell types:
#   EPI_SCENIC_CELLTYPES="Goblet,Ciliated_Mature,Basal_Progenitor" Rscript ...
#
#   # Force re-run:
#   EPI_SCENIC_FORCE=1 Rscript ...
#
#   # Dry-run / preflight only:
#   EPI_SCENIC_PREFLIGHT_ONLY=1 Rscript ...
# ===============================================================================

Sys.setenv(
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS  = "1"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
EPITHELIAL_RDS <- "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.rds"
CORE_SCRIPT <- "/home/h2048/script/R/scenic_core_20260410.R"
DATABASE_DIR <- "/home/h2048/data/index_genome/cisTarget_databases_rscenic"
SCENIC_DB_10KB <- file.path(DATABASE_DIR, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
SCENIC_DB_500BP <- file.path(DATABASE_DIR, "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

RUN_DATE <- "20260602"
RUN_ROOT <- file.path("/home/h2048/output/program_full_parallel_methods_20260507", "scenic_epithelial")
OUTPUT_ROOT <- RUN_ROOT
LOG_DIR <- file.path("/home/h2048/logs", RUN_DATE)

# ------------------------------------------------------------------------------
# Parameter parsing from environment
# ------------------------------------------------------------------------------
parse_csv_env <- function(name, default = character()) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  vals <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

parse_bool_env <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, unset = if (isTRUE(default)) "1" else "0")
  tolower(raw) %in% c("1", "true", "yes", "y")
}

parse_int_env <- function(name, default) {
  val <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(val) == 0L || is.na(val)) default else val
}

selected_levels <- parse_csv_env("EPI_SCENIC_LEVELS", c("lineage_celltype", "celltype_tissue"))
selected_celltypes <- parse_csv_env("EPI_SCENIC_CELLTYPES", character())
force_run <- parse_bool_env("EPI_SCENIC_FORCE", FALSE)
preflight_only <- parse_bool_env("EPI_SCENIC_PREFLIGHT_ONLY", FALSE)
stop_after_first_error <- parse_bool_env("EPI_SCENIC_STOP_AFTER_FIRST_ERROR", FALSE)
scenic_n_cores <- parse_int_env("EPI_SCENIC_N_CORES", 8L)
min_cells <- parse_int_env("EPI_SCENIC_MIN_CELLS", 50L)
min_groups <- parse_int_env("EPI_SCENIC_MIN_GROUPS", 2L)
min_group_cells <- parse_int_env("EPI_SCENIC_MIN_GROUP_CELLS", 20L)

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
stopifnot(file.exists(EPITHELIAL_RDS))
stopifnot(file.exists(CORE_SCRIPT))
stopifnot(dir.exists(DATABASE_DIR))
stopifnot(file.exists(SCENIC_DB_10KB))
stopifnot(file.exists(SCENIC_DB_500BP))

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

SCENIC_SOURCE_ONLY <- TRUE
source(CORE_SCRIPT)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
safe_file_id <- function(x) {
  out <- gsub("[^A-Za-z0-9_]+", "_", trimws(as.character(x)))
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  out[!nzchar(out)] <- "NA"
  out
}

status_msg <- function(...) {
  txt <- sprintf(...)
  cat(txt)
  try(flush(stdout()), silent = TRUE)
  try(flush.console(), silent = TRUE)
  invisible(txt)
}

write_json_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}

read_status <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

set_global <- function(name, value) {
  assign(name, value, envir = .GlobalEnv)
}

close_output_sinks <- function() {
  while (sink.number(type = "output") > 0L) sink(type = "output")
}

normalize_group <- function(x, unknown = "Unknown") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- unknown
  x
}

summarize_run_outputs <- function(output_dir) {
  target_path <- file.path(output_dir, "tables", "regulon_targets_full.csv")
  top_path <- file.path(output_dir, "tables", "cellwise_top_regulon.csv")
  auc_path <- file.path(output_dir, "rds", "regulon_auc_matrix_full_object.rds")
  result_path <- file.path(output_dir, "rds", "scenic_full_result.rds")
  n_regulons <- NA_integer_
  n_targets <- NA_integer_
  n_scored_cells <- NA_integer_
  if (file.exists(target_path)) {
    targets <- tryCatch(fread(target_path), error = function(e) NULL)
    if (!is.null(targets) && nrow(targets) > 0L) {
      n_regulons <- uniqueN(targets$regulon)
      n_targets <- nrow(targets)
    }
  }
  if (file.exists(top_path)) {
    top <- tryCatch(fread(top_path), error = function(e) NULL)
    if (!is.null(top)) n_scored_cells <- nrow(top)
  }
  list(
    n_regulons = n_regulons,
    n_targets = n_targets,
    n_scored_cells = n_scored_cells,
    regulon_targets_csv = target_path,
    auc_matrix_rds = auc_path,
    scenic_result_rds = result_path
  )
}

write_job_metadata <- function(obj_job, lineage_name, output_dir, comparison_level, unit_label, group_col = "scenic_comparison_group") {
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  meta <- obj_job@meta.data
  tissue_col <- if ("tissue" %in% colnames(meta)) "tissue" else NULL
  out <- data.table(
    cell = rownames(meta),
    lineage = lineage_name,
    comparison_level = comparison_level,
    unit_label = unit_label,
    cell_type_L3 = as.character(meta[["cell_type_L3"]]),
    scenic_group = as.character(meta[[group_col]]),
    sample = as.character(meta[["sample"]])
  )
  if (!is.null(tissue_col)) {
    out[, tissue := as.character(meta[[tissue_col]])]
  }
  fwrite(out, file.path(output_dir, "tables", "cell_metadata_used_for_scenic.csv"))
  group_summ <- out[, .N, by = .(lineage, comparison_level, unit_label, scenic_group)][order(-N)]
  fwrite(group_summ, file.path(output_dir, "tables", "comparison_group_summary.csv"))
  if (!is.null(tissue_col)) {
    tissue_summ <- out[, .N, by = .(lineage, comparison_level, unit_label, cell_type_L3, tissue)][order(cell_type_L3, tissue)]
    fwrite(tissue_summ, file.path(output_dir, "tables", "celltype_tissue_summary.csv"))
  }
  fwrite(out[, .N, by = .(lineage, comparison_level, unit_label, sample, scenic_group)][order(sample, scenic_group)],
         file.path(output_dir, "tables", "sample_group_summary.csv"))
  invisible(out)
}

# ==============================================================================
# Epithelial-specific downsampling strategy
# ==============================================================================
# Epithelial has 278K cells, 20 L3 types, 4 tissues, 150 samples.
# GENIE3 scales poorly, so we need aggressive but stratified downsampling.
#
# For lineage_celltype level (compare cell types):
#   - Keep max 80 cells per (sample, cell_type) stratum
#   - Then cap to 600 cells per cell_type
#   - Global cap: 6000 cells
#   - Gene cap: 3000 most variable genes
#
# For celltype_tissue level (compare tissues within cell type):
#   - Keep max 120 cells per (sample, tissue) stratum
#   - Then cap to 1200 cells per tissue
#   - Global cap: 4000 cells per cell-type run
#   - Gene cap: 3000 most variable genes

epithelial_downsample_inference_cells <- function(
  obj,
  cell_type_col,
  sample_col,
  tissue_col = NULL,
  comparison_level = c("lineage_celltype", "celltype_tissue"),
  seed = 42L
) {
  comparison_level <- match.arg(comparison_level)
  meta <- obj@meta.data
  set.seed(seed)

  if (identical(comparison_level, "lineage_celltype")) {
    max_per_sample_ct <- 80L
    max_per_ct <- 600L
    global_max <- 6000L
  } else {
    max_per_sample_ct <- 120L
    max_per_ct <- 1200L
    global_max <- 4000L
  }

  ct <- as.character(meta[[cell_type_col]])
  sp <- as.character(meta[[sample_col]])
  keep <- !is.na(ct) & nzchar(trimws(ct)) & !is.na(sp) & nzchar(trimws(sp))
  meta <- meta[keep, , drop = FALSE]
  meta$cell_id <- rownames(meta)
  meta$stratum <- meta[[cell_type_col]]

  if (!is.null(tissue_col) && tissue_col %in% colnames(meta)) {
    if (identical(comparison_level, "celltype_tissue")) {
      meta$stratum <- as.character(meta[[tissue_col]])
    }
  }

  # Step 1: Downsample within (sample, stratum)
  meta$sample_stratum <- paste(meta$sample, meta$stratum, sep = "___")
  picked1 <- unlist(lapply(split(meta$cell_id, meta$sample_stratum), function(v) {
    sample(v, min(length(v), max_per_sample_ct))
  }), use.names = FALSE)

  meta1 <- meta[meta$cell_id %in% picked1, , drop = FALSE]

  # Step 2: Downsample within each stratum
  picked2 <- unlist(lapply(split(meta1$cell_id, meta1$stratum), function(v) {
    sample(v, min(length(v), max_per_ct))
  }), use.names = FALSE)

  # Step 3: Global cap
  picked2 <- unique(picked2)
  if (length(picked2) > global_max) {
    picked2 <- sample(picked2, global_max)
  }

  if (length(picked2) < 2) {
    stop(sprintf("Too few cells after downsampling: %d", length(picked2)))
  }

  status_msg("[DOWNSAMPLE] %s: %d cells kept (%d -> %d, max_per_sample=%d, max_per_group=%d, global=%d)\n",
    comparison_level,
    length(picked2),
    nrow(meta),
    length(picked2),
    max_per_sample_ct,
    max_per_ct,
    global_max
  )

  picked2
}

# ==============================================================================
# Level 1: Lineage Cell Type Comparison
# ==============================================================================
run_lineage_celltype_level <- function(obj, force = FALSE) {
  output_dir <- file.path(OUTPUT_ROOT, "lineage_celltype", "scenic_epithelial_cell_type_L3")
  status_path <- file.path(output_dir, "epi_scenic_lineage_status.json")
  log_path <- file.path(output_dir, "epi_scenic_lineage_wrapper.log")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  existing <- read_status(status_path)
  if (!isTRUE(force) && !is.null(existing) && identical(existing$status, "ok")) {
    status_msg("[SKIP] lineage_celltype existing ok\n")
    return(as.data.table(existing))
  }

  meta <- obj@meta.data
  ct_labels <- normalize_group(meta[["cell_type_L3"]], unknown = "Unknown")
  n_total_cells <- ncol(obj)
  n_groups <- uniqueN(ct_labels[ct_labels != "Unknown"])
  n_samples <- uniqueN(as.character(meta[["sample"]]))
  n_tissues <- if ("tissue" %in% colnames(meta)) uniqueN(as.character(meta[["tissue"]])) else NA_integer_

  started <- Sys.time()
  base_status <- list(
    status = "started",
    comparison_level = "lineage_celltype",
    lineage = "epithelial",
    unit_label = "epithelial__cell_type_L3",
    grouping_col = "cell_type_L3",
    n_cells = as.integer(n_total_cells),
    n_groups = as.integer(n_groups),
    n_samples = as.integer(n_samples),
    n_tissues = as.integer(n_tissues),
    min_cells_required = as.integer(min_cells),
    output_dir = output_dir,
    log_path = log_path,
    started_at = format(started, "%Y-%m-%d %H:%M:%S"),
    ended_at = NULL,
    elapsed_min = NULL,
    error = NULL
  )
  write_json_safe(base_status, status_path)

  if (isTRUE(preflight_only)) {
    base_status$status <- "preflight_only"
    base_status$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write_json_safe(base_status, status_path)
    status_msg("[PREFLIGHT] lineage_celltype: %d cells, %d groups, %d samples, %d tissues\n",
               n_total_cells, n_groups, n_samples, n_tissues)
    return(as.data.table(base_status))
  }

  result_status <- tryCatch({
    obj_job <- obj
    obj_job@meta.data$scenic_comparison_group <- normalize_group(obj_job@meta.data[["cell_type_L3"]], unknown = "Unknown")
    write_job_metadata(obj_job, "epithelial", output_dir, "lineage_celltype", "epithelial__cell_type_L3")

    set_global("INPUT_RDS", EPITHELIAL_RDS)
    set_global("OUTPUT_DIR", output_dir)
    set_global("DATABASE_DIR", DATABASE_DIR)
    set_global("SCENIC_DB_10KB", SCENIC_DB_10KB)
    set_global("SCENIC_DB_500BP", SCENIC_DB_500BP)
    set_global("SCENIC_DB_SCOPE", "both")
    set_global("SCENIC_DB_INDEX_COL", NULL)
    set_global("ORGANISM", "hgnc")
    set_global("ASSAY_USE", "RNA")
    set_global("CELL_TYPE_COL", "scenic_comparison_group")
    set_global("SAMPLE_COL", "sample")
    set_global("N_CORES", scenic_n_cores)
    set_global("SCENIC_DATASET_TITLE", "epithelial_lineage_celltype_L3_SCENIC_20260602")
    set_global("REDUCTION_CANDIDATES", c("umap_scanvi", "umap_harmony", "umap", "harmony", "pca"))
    set_global("EXPORT_FULL_AUC_MATRIX_CSV", TRUE)
    set_global("INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE", 80L)
    set_global("INFERENCE_MAX_CELLS_PER_CELLTYPE", 600L)
    set_global("INFERENCE_GLOBAL_MAX_CELLS", 6000L)
    set_global("INFERENCE_MAX_GENES", 3000L)
    set_global("MAX_INFERENCE_DENSE_GB", 10)
    set_global("MIN_GENES_PER_CELL", 200L)
    set_global("MIN_CELLS_PER_GENE", 10L)
    set_global("MIN_GENE_PCT", 0.01)
    set_global("MIN_GENES_PER_REGULON", 20L)
    set_global("SCENIC_GENIE3_TREE_METHOD", "ET")
    set_global("SCENIC_GENIE3_NTREES", 100L)
    set_global("SCENIC_MODULE_WEIGHT_THRESHOLD", 0.005)
    set_global("SCENIC_MODULE_N_TOP_TFS", 5L)
    set_global("SCENIC_MODULE_N_TOP_TARGETS", 30L)
    set_global("SCENIC_REGULON_COEX_METHODS", "top5perTarget")

    sink(log_path, split = TRUE)
    on.exit(close_output_sinks(), add = TRUE)
    status_msg("[RUN] lineage_celltype epithelial: %d cells, %d groups, %d samples\n",
               n_total_cells, n_groups, n_samples)
    scenic_result <- run_scenic_module(obj_job)
    rm(scenic_result); gc(verbose = FALSE)
    list(status = "ok")
  }, error = function(e) {
    close_output_sinks()
    list(status = "error", error = conditionMessage(e))
  })

  ended <- Sys.time()
  base_status$status <- result_status$status
  base_status$error <- result_status$error %||% NULL
  base_status$ended_at <- format(ended, "%Y-%m-%d %H:%M:%S")
  base_status$elapsed_min <- round(as.numeric(difftime(ended, started, units = "mins")), 3)
  if (identical(result_status$status, "ok")) {
    base_status <- c(base_status, summarize_run_outputs(output_dir))
  }
  write_json_safe(base_status, status_path)
  status_msg("[%s] lineage_celltype elapsed_min=%.2f\n", toupper(base_status$status), base_status$elapsed_min)
  as.data.table(base_status)
}

# ==============================================================================
# Level 2: Cell Type Tissue Comparison
# ==============================================================================
run_celltype_tissue_level <- function(obj, ct_label, force = FALSE) {
  safe_ct <- safe_file_id(ct_label)
  output_dir <- file.path(OUTPUT_ROOT, "celltype_tissue", paste0("scenic_", safe_ct))
  status_path <- file.path(output_dir, "epi_scenic_celltype_tissue_status.json")
  log_path <- file.path(output_dir, "epi_scenic_celltype_tissue_wrapper.log")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  existing <- read_status(status_path)
  if (!isTRUE(force) && !is.null(existing) && identical(existing$status, "ok")) {
    status_msg("[SKIP] celltype_tissue/%s existing ok\n", ct_label)
    return(as.data.table(existing))
  }

  # Subset to cell type
  meta <- obj@meta.data
  cells <- rownames(meta)[normalize_group(meta[["cell_type_L3"]]) == ct_label]
  if (length(cells) == 0L) {
    dt <- data.table(status = "error", error = sprintf("No cells for %s", ct_label))
    return(dt)
  }
  obj_ct <- subset(obj, cells = cells)
  n_ct_cells <- ncol(obj_ct)
  tissue_vals <- normalize_group(obj_ct@meta.data[["tissue"]], unknown = "Unknown")

  # Filter tissues with enough cells
  tissue_counts <- table(tissue_vals)
  valid_tissues <- names(tissue_counts[tissue_counts >= min_group_cells])
  if (length(valid_tissues) < min_groups) {
    base_status <- list(
      status = "skipped_too_few_tissues",
      comparison_level = "celltype_tissue",
      lineage = "epithelial",
      unit_label = ct_label,
      n_cells = as.integer(n_ct_cells),
      n_tissues = as.integer(length(unique(tissue_vals))),
      n_valid_tissues = as.integer(length(valid_tissues)),
      min_tissues_required = as.integer(min_groups),
      output_dir = output_dir,
      started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    write_json_safe(base_status, status_path)
    status_msg("[SKIP] celltype_tissue/%s too few valid tissues: %d < %d\n", ct_label, length(valid_tissues), min_groups)
    return(as.data.table(base_status))
  }

  keep_cells <- rownames(obj_ct@meta.data)[tissue_vals %in% valid_tissues]
  obj_ct <- subset(obj_ct, cells = keep_cells)
  n_ct_cells <- ncol(obj_ct)

  started <- Sys.time()
  base_status <- list(
    status = "started",
    comparison_level = "celltype_tissue",
    lineage = "epithelial",
    unit_label = ct_label,
    safe_unit_label = safe_ct,
    grouping_col = "tissue",
    n_cells = as.integer(n_ct_cells),
    n_tissues = as.integer(length(valid_tissues)),
    n_samples = as.integer(uniqueN(as.character(obj_ct@meta.data[["sample"]]))),
    min_cells_required = as.integer(min_cells),
    min_groups_required = as.integer(min_groups),
    min_group_cells_required = as.integer(min_group_cells),
    output_dir = output_dir,
    log_path = log_path,
    started_at = format(started, "%Y-%m-%d %H:%M:%S"),
    ended_at = NULL,
    elapsed_min = NULL,
    error = NULL
  )
  write_json_safe(base_status, status_path)

  if (n_ct_cells < min_cells) {
    base_status$status <- "skipped_too_few_cells"
    base_status$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write_json_safe(base_status, status_path)
    return(as.data.table(base_status))
  }

  if (isTRUE(preflight_only)) {
    base_status$status <- "preflight_only"
    base_status$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write_json_safe(base_status, status_path)
    status_msg("[PREFLIGHT] celltype_tissue/%s: %d cells, %d tissues\n", ct_label, n_ct_cells, length(valid_tissues))
    return(as.data.table(base_status))
  }

  result_status <- tryCatch({
    obj_ct@meta.data$scenic_comparison_group <- normalize_group(obj_ct@meta.data[["tissue"]], unknown = "Unknown")
    write_job_metadata(obj_ct, "epithelial", output_dir, "celltype_tissue", ct_label)

    set_global("INPUT_RDS", EPITHELIAL_RDS)
    set_global("OUTPUT_DIR", output_dir)
    set_global("DATABASE_DIR", DATABASE_DIR)
    set_global("SCENIC_DB_10KB", SCENIC_DB_10KB)
    set_global("SCENIC_DB_500BP", SCENIC_DB_500BP)
    set_global("SCENIC_DB_SCOPE", "both")
    set_global("SCENIC_DB_INDEX_COL", NULL)
    set_global("ORGANISM", "hgnc")
    set_global("ASSAY_USE", "RNA")
    set_global("CELL_TYPE_COL", "scenic_comparison_group")
    set_global("SAMPLE_COL", "sample")
    set_global("N_CORES", scenic_n_cores)
    set_global("SCENIC_DATASET_TITLE", sprintf("epithelial_celltype_tissue_%s_SCENIC_20260602", safe_ct))
    set_global("REDUCTION_CANDIDATES", c("umap_scanvi", "umap_harmony", "umap", "harmony", "pca"))
    set_global("EXPORT_FULL_AUC_MATRIX_CSV", TRUE)
    set_global("INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE", 120L)
    set_global("INFERENCE_MAX_CELLS_PER_CELLTYPE", 1200L)
    set_global("INFERENCE_GLOBAL_MAX_CELLS", 4000L)
    set_global("INFERENCE_MAX_GENES", 3000L)
    set_global("MAX_INFERENCE_DENSE_GB", 10)
    set_global("MIN_GENES_PER_CELL", 200L)
    set_global("MIN_CELLS_PER_GENE", 10L)
    set_global("MIN_GENE_PCT", 0.01)
    set_global("MIN_GENES_PER_REGULON", 20L)
    set_global("SCENIC_GENIE3_TREE_METHOD", "ET")
    set_global("SCENIC_GENIE3_NTREES", 100L)
    set_global("SCENIC_MODULE_WEIGHT_THRESHOLD", 0.005)
    set_global("SCENIC_MODULE_N_TOP_TFS", 5L)
    set_global("SCENIC_MODULE_N_TOP_TARGETS", 30L)
    set_global("SCENIC_REGULON_COEX_METHODS", "top5perTarget")

    sink(log_path, split = TRUE)
    on.exit(close_output_sinks(), add = TRUE)
    status_msg("[RUN] celltype_tissue/%s: %d cells, %d tissues\n", ct_label, n_ct_cells, length(valid_tissues))
    scenic_result <- run_scenic_module(obj_ct)
    rm(scenic_result); gc(verbose = FALSE)
    list(status = "ok")
  }, error = function(e) {
    close_output_sinks()
    list(status = "error", error = conditionMessage(e))
  })

  ended <- Sys.time()
  base_status$status <- result_status$status
  base_status$error <- result_status$error %||% NULL
  base_status$ended_at <- format(ended, "%Y-%m-%d %H:%M:%S")
  base_status$elapsed_min <- round(as.numeric(difftime(ended, started, units = "mins")), 3)
  if (identical(result_status$status, "ok")) {
    base_status <- c(base_status, summarize_run_outputs(output_dir))
  }
  write_json_safe(base_status, status_path)
  status_msg("[%s] celltype_tissue/%s elapsed_min=%.2f\n", toupper(base_status$status), ct_label, base_status$elapsed_min)
  as.data.table(base_status)
}

# ==============================================================================
# Post-hoc: Epithelial regulon summary
# ==============================================================================
generate_epithelial_regulon_summary <- function() {
  status_msg("\n[SUMMARY] Generating epithelial regulon summary...\n")
  summary_list <- list()

  # Gather lineage_celltype results
  lc_dir <- file.path(OUTPUT_ROOT, "lineage_celltype", "scenic_epithelial_cell_type_L3")
  lc_targets <- file.path(lc_dir, "tables", "regulon_targets_full.csv")
  lc_rss <- file.path(lc_dir, "rds", "regulon_rss.rds")
  lc_auc <- file.path(lc_dir, "rds", "regulon_auc_matrix_full_object.rds")

  if (file.exists(lc_targets)) {
    targets <- fread(lc_targets)
    summary_list$lineage_celltype <- list(
      n_regulons = uniqueN(targets$regulon),
      n_targets = nrow(targets),
      top_tfs = targets[, .N, by = tf][order(-N)][1:min(20, .N)][, tf]
    )
  }

  # Gather celltype_tissue results
  ct_dirs <- list.dirs(file.path(OUTPUT_ROOT, "celltype_tissue"), recursive = FALSE, full.names = TRUE)
  ct_results <- lapply(ct_dirs, function(d) {
    status_path <- file.path(d, "epi_scenic_celltype_tissue_status.json")
    status <- read_status(status_path)
    if (is.null(status)) return(NULL)
    list(
      celltype = status$unit_label,
      status = status$status,
      n_regulons = status$n_regulons,
      n_cells = status$n_cells,
      n_tissues = status$n_tissues,
      elapsed_min = status$elapsed_min
    )
  })
  ct_results <- Filter(Negate(is.null), ct_results)
  summary_list$celltype_tissue <- ct_results

  # Write summary
  summary_path <- file.path(OUTPUT_ROOT, "epithelial_scenic_summary_20260602.json")
  write_json_safe(summary_list, summary_path)

  # Write per-celltype summary table
  if (length(ct_results) > 0L) {
    ct_dt <- rbindlist(ct_results, fill = TRUE)
    fwrite(ct_dt, file.path(OUTPUT_ROOT, "epithelial_scenic_celltype_tissue_summary_20260602.tsv"), sep = "\t")
  }

  status_msg("[SUMMARY] Written to %s\n", summary_path)
  invisible(summary_list)
}

# ==============================================================================
# Main
# ==============================================================================
main <- function() {
  # Initialize log file
  master_log <- file.path(LOG_DIR, "epithelial_scenic_20260602.log")
  cat(sprintf("[BOOT] Epithelial SCENIC runner started at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = master_log)

  sink(master_log, split = TRUE, append = TRUE)
  on.exit(close_output_sinks(), add = TRUE)

  status_msg("============================================================\n")
  status_msg("Epithelial SCENIC Runner\n")
  status_msg("============================================================\n")
  status_msg("Levels: %s\n", paste(selected_levels, collapse = ", "))
  status_msg("N cores: %d\n", scenic_n_cores)
  status_msg("Force: %s\n", force_run)
  status_msg("Preflight: %s\n", preflight_only)
  status_msg("Output root: %s\n", OUTPUT_ROOT)
  status_msg("Master log: %s\n", master_log)
  status_msg("============================================================\n\n")

  # Validate inputs
  if (!"lineage_celltype" %in% selected_levels && !"celltype_tissue" %in% selected_levels) {
    stop("No valid comparison levels selected.")
  }
  unknown <- setdiff(selected_levels, c("lineage_celltype", "celltype_tissue"))
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown levels: %s", paste(unknown, collapse = ", ")))
  }

  # Load data
  status_msg("[LOAD] %s\n", EPITHELIAL_RDS)
  obj <- readRDS(EPITHELIAL_RDS)
  status_msg("[LOAD] %d cells, %d genes\n", ncol(obj), nrow(obj))

  # Validate columns
  required_cols <- c("cell_type_L3", "sample", "tissue")
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Missing metadata columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Determine cell types for tissue comparison
  meta <- obj@meta.data
  ct_counts <- data.table(
    celltype = as.character(meta[["cell_type_L3"]]),
    tissue = as.character(meta[["tissue"]])
  )
  ct_counts <- ct_counts[!is.na(celltype) & nzchar(trimws(celltype))]
  ct_summary <- ct_counts[, .(
    n_cells = .N,
    n_tissues = uniqueN(tissue)
  ), by = celltype][order(-n_cells)]

  # Filter for celltype_tissue: need >=2 tissues with >=min_group_cells each
  ct_tissue_detail <- ct_counts[, .(
    n_cells = .N
  ), by = .(celltype, tissue)][n_cells >= min_group_cells]
  ct_valid <- ct_tissue_detail[, .(
    n_valid_tissues = .N,
    min_tissue_cells = min(n_cells)
  ), by = celltype][n_valid_tissues >= min_groups]

  # Filter by user selection
  if (length(selected_celltypes) > 0L) {
    ct_valid <- ct_valid[celltype %in% selected_celltypes]
  }

  status_msg("\n[CELLTYPES] Eligible for tissue comparison (>=%d tissues, >=%d cells/tissue):\n",
             min_groups, min_group_cells)
  for (i in seq_len(nrow(ct_valid))) {
    status_msg("  %-40s tissues=%d min_cells=%d\n",
               ct_valid$celltype[i], ct_valid$n_valid_tissues[i], ct_valid$min_tissue_cells[i])
  }
  status_msg("  Total eligible: %d / %d cell types\n\n", nrow(ct_valid), nrow(ct_summary))

  # Run audit tracker
  audit_path <- file.path(OUTPUT_ROOT, "epithelial_scenic_audit_20260602.tsv")
  all_audit <- list()

  append_audit <- function(row) {
    all_audit[[length(all_audit) + 1L]] <<- row
    current <- rbindlist(all_audit, fill = TRUE)
    fwrite(current, audit_path, sep = "\t")
    if (identical(row$status[[1]], "error") && isTRUE(stop_after_first_error)) {
      stop(sprintf("Stopping after error: %s / %s", row$comparison_level[[1]], row$unit_label[[1]]), call. = FALSE)
    }
  }

  # ----- Level 1: lineage_celltype -----
  if ("lineage_celltype" %in% selected_levels) {
    status_msg("\n========== LEVEL 1: lineage_celltype ==========\n")
    row <- run_lineage_celltype_level(obj, force = force_run)
    append_audit(row)
  }

  # ----- Level 2: celltype_tissue -----
  if ("celltype_tissue" %in% selected_levels) {
    status_msg("\n========== LEVEL 2: celltype_tissue ==========\n")
    status_msg("Running for %d eligible cell types\n", nrow(ct_valid))

    for (i in seq_len(nrow(ct_valid))) {
      ct_label <- ct_valid$celltype[i]
      status_msg("\n--- celltype_tissue [%d/%d]: %s ---\n", i, nrow(ct_valid), ct_label)
      row <- run_celltype_tissue_level(obj, ct_label, force = force_run)
      append_audit(row)
    }
  }

  # ----- Summary -----
  status_msg("\n========== SUMMARY ==========\n")
  audit <- rbindlist(all_audit, fill = TRUE)
  fwrite(audit, audit_path, sep = "\t")
  print(audit[, .N, by = status][order(status)])

  generate_epithelial_regulon_summary()

  status_msg("\n[DONE] Epithelial SCENIC runner finished at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  status_msg("Audit: %s\n", audit_path)
  status_msg("Master log: %s\n", master_log)

  invisible(audit)
}

# Execute
if (sys.nframe() == 0L) {
  main()
}
