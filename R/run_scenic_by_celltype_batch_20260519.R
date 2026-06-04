#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Run R-SCENIC for every final L3 cell type across all healthy tissue lineages
# ==============================================================================
# Output layout:
#   /home/h2048/output/program_full_parallel_methods_20260507/scenic_by_celltype/
#     <lineage>/<celltype_safe>/scenic_full/
#       scenic_celltype_status.json
#       scenic_celltype_wrapper.log
#       tables/ figures/ rds/ int/
#
# Resume controls:
#   SCENIC_FORCE=1                 Re-run even if status is ok
#   SCENIC_LINEAGES=bcell,tnk      Restrict lineages
#   SCENIC_CELLTYPES=CD4_Tcm,...   Restrict exact L3 labels across selected lineages
#   SCENIC_MAX_CELLTYPES=1         Smoke-test limit after filtering
#   SCENIC_PREFLIGHT_ONLY=1        Write metadata/status only; do not run SCENIC
#   SCENIC_N_CORES=4               Cores used inside each SCENIC run
#   SCENIC_GENIE3_NTREES=100       Fast-mode GENIE3 tree count (override as needed)
#   SCENIC_GENIE3_TREE_METHOD=ET   Fast-mode tree method; RF is the original default
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

RUN_ROOT <- "/home/h2048/output/program_full_parallel_methods_20260507"
SCENIC_ROOT <- file.path(RUN_ROOT, "scenic_by_celltype")
SCENIC_METHOD_DIRNAME <- "scenic_full"
AUDIT_TSV <- file.path(SCENIC_ROOT, "scenic_celltype_audit_20260519.tsv")
CORE_SCRIPT <- "/home/h2048/script/R/scenic_core_20260410.R"
DATABASE_DIR_DEFAULT <- "/home/h2048/data/index_genome/cisTarget_databases_rscenic"
SCENIC_DB_10KB_DEFAULT <- file.path(DATABASE_DIR_DEFAULT, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
SCENIC_DB_500BP_DEFAULT <- file.path(DATABASE_DIR_DEFAULT, "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

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

parse_numeric_vector_env <- function(name, default = NULL) {
	raw <- Sys.getenv(name, unset = "")
	if (!nzchar(raw)) return(default)
	vals <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
	vals <- vals[nzchar(vals)]
	out <- suppressWarnings(as.numeric(vals))
	if (any(is.na(out))) stop(sprintf("%s must be a comma-separated numeric vector.", name), call. = FALSE)
	out
}

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

celltype_method_output_dir <- function(root_dir, lineage, safe_celltype, method_dirname = SCENIC_METHOD_DIRNAME) {
	file.path(root_dir, lineage, safe_celltype, method_dirname)
}

read_status <- function(path) {
	if (!file.exists(path)) return(NULL)
	tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

set_global <- function(name, value) {
	assign(name, value, envir = .GlobalEnv)
}

lineage_configs <- list(
	bcell = list(
		lineage = "bcell",
		rds_path = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/bcell_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/bcell_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_refined", "umap_scanvi", "umap", "harmony", "pca")
	),
	stromal_smc = list(
		lineage = "stromal_smc",
		rds_path = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_scanvi", "umap", "harmony", "pca")
	),
	stromal_fibroblast = list(
		lineage = "stromal_fibroblast",
		rds_path = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_scanvi", "umap", "harmony", "pca")
	),
	stromal_endothelial = list(
		lineage = "stromal_endothelial",
		rds_path = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_endothelial_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_endothelial_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_scanvi", "umap", "harmony", "pca")
	),
	tnk = list(
		lineage = "tnk",
		rds_path = "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_refined", "umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap", "harmony", "pca")
	),
	myeloid = list(
		lineage = "myeloid",
		rds_path = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_refined", "umap_scanvi", "umap", "harmony", "pca")
	),
	epithelial = list(
		lineage = "epithelial",
		rds_path = "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.rds",
		h5ad_path = "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.h5ad",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		reductions = c("umap_refined", "umap_scanvi", "umap", "harmony", "pca")
	)
)

selected_lineages <- parse_csv_env("SCENIC_LINEAGES", names(lineage_configs))
unknown_lineages <- setdiff(selected_lineages, names(lineage_configs))
if (length(unknown_lineages) > 0L) {
	stop(sprintf("Unknown SCENIC_LINEAGES: %s", paste(unknown_lineages, collapse = ", ")), call. = FALSE)
}
selected_celltypes <- parse_csv_env("SCENIC_CELLTYPES", character())
force_run <- parse_bool_env("SCENIC_FORCE", FALSE)
preflight_only <- parse_bool_env("SCENIC_PREFLIGHT_ONLY", FALSE)
stop_after_first_error <- parse_bool_env("SCENIC_STOP_AFTER_FIRST_ERROR", FALSE)
fast_modules <- parse_bool_env("SCENIC_FAST_MODULES", TRUE)
min_cells <- parse_int_env("SCENIC_MIN_CELLS", 50L)
max_celltypes <- parse_int_env("SCENIC_MAX_CELLTYPES", 0L)
scenic_n_cores <- parse_int_env("SCENIC_N_CORES", 4L)

dir.create(SCENIC_ROOT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(CORE_SCRIPT)) stop(sprintf("Missing SCENIC core script: %s", CORE_SCRIPT), call. = FALSE)
if (!dir.exists(DATABASE_DIR_DEFAULT)) stop(sprintf("Missing SCENIC database dir: %s", DATABASE_DIR_DEFAULT), call. = FALSE)
if (!file.exists(SCENIC_DB_10KB_DEFAULT) || !file.exists(SCENIC_DB_500BP_DEFAULT)) {
	stop("Missing compatible R-SCENIC feather database files.", call. = FALSE)
}

SCENIC_SOURCE_ONLY <- TRUE
source(CORE_SCRIPT)

list_celltypes <- function(obj, cfg) {
	meta <- obj@meta.data
	ct_col <- cfg$celltype_col
	sample_col <- cfg$sample_col
	missing_cols <- setdiff(c(ct_col, sample_col), colnames(meta))
	if (length(missing_cols) > 0L) stop(sprintf("%s missing metadata columns: %s", cfg$lineage, paste(missing_cols, collapse = ", ")), call. = FALSE)
	dt <- data.table(
		cell = rownames(meta),
		celltype = as.character(meta[[ct_col]]),
		sample = as.character(meta[[sample_col]])
	)
	dt <- dt[!is.na(celltype) & nzchar(trimws(celltype))]
	dt[, .(n_cells = .N, n_samples = uniqueN(sample[!is.na(sample) & nzchar(sample)])), by = celltype][order(-n_cells, celltype)]
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

write_cell_metadata <- function(obj_ct, cfg, output_dir) {
	dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
	meta <- obj_ct@meta.data
	out <- data.table(
		cell = rownames(meta),
		lineage = cfg$lineage,
		cell_type_L3 = as.character(meta[[cfg$celltype_col]]),
		scenic_celltype_run = as.character(meta$scenic_celltype_run),
		sample = as.character(meta[[cfg$sample_col]]),
		tissue = if ("tissue" %in% colnames(meta)) as.character(meta$tissue) else NA_character_
	)
	fwrite(out, file.path(output_dir, "tables", "cell_metadata_used_for_scenic_celltype_run.csv"))
	tissue_summary <- out[, .N, by = .(lineage, cell_type_L3, tissue)][order(lineage, cell_type_L3, tissue)]
	fwrite(tissue_summary, file.path(output_dir, "tables", "celltype_tissue_summary.csv"))
	sample_summary <- out[, .N, by = .(lineage, cell_type_L3, sample, tissue)][order(lineage, cell_type_L3, tissue, sample)]
	fwrite(sample_summary, file.path(output_dir, "tables", "celltype_sample_summary.csv"))
	invisible(out)
}

close_output_sinks <- function() {
	while (sink.number(type = "output") > 0L) sink(type = "output")
}

run_one_celltype <- function(obj, cfg, ct_label, n_cells, n_samples) {
	safe_ct <- safe_file_id(ct_label)
	output_dir <- celltype_method_output_dir(SCENIC_ROOT, cfg$lineage, safe_ct)
	status_path <- file.path(output_dir, "scenic_celltype_status.json")
	log_path <- file.path(output_dir, "scenic_celltype_wrapper.log")
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

	existing <- read_status(status_path)
	if (!isTRUE(force_run) && !is.null(existing) && identical(existing$status, "ok")) {
		status_msg("[SKIP] %s / %s existing ok\n", cfg$lineage, ct_label)
		return(as.data.table(existing))
	}

	started <- Sys.time()
	base_status <- list(
		status = "started",
		lineage = cfg$lineage,
		celltype = ct_label,
		safe_celltype = safe_ct,
		path_layout = "celltype_method",
		n_cells = as.integer(n_cells),
		n_samples = as.integer(n_samples),
		min_cells_required = as.integer(min_cells),
		output_dir = output_dir,
		log_path = log_path,
		started_at = format(started, "%Y-%m-%d %H:%M:%S"),
		ended_at = NULL,
		elapsed_min = NULL,
		error = NULL
	)
	write_json_safe(base_status, status_path)

	if (is.finite(min_cells) && n_cells < min_cells) {
		base_status$status <- "skipped_too_few_cells"
		base_status$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
		write_json_safe(base_status, status_path)
		status_msg("[SKIP] %s / %s too few cells: %d < %d\n", cfg$lineage, ct_label, n_cells, min_cells)
		return(as.data.table(base_status))
	}

	result_status <- tryCatch({
		cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[cfg$celltype_col]]) == ct_label]
		if (length(cells) == 0L) stop("No cells matched cell type label after object load.")
		obj_ct <- subset(obj, cells = cells)
		obj_ct@meta.data[["scenic_celltype_run"]] <- ct_label
		write_cell_metadata(obj_ct, cfg, output_dir)

		set_global("INPUT_RDS", cfg$rds_path)
		set_global("OUTPUT_DIR", output_dir)
		set_global("DATABASE_DIR", DATABASE_DIR_DEFAULT)
		set_global("SCENIC_DB_10KB", SCENIC_DB_10KB_DEFAULT)
		set_global("SCENIC_DB_500BP", SCENIC_DB_500BP_DEFAULT)
		set_global("SCENIC_DB_INDEX_COL", NULL)
		set_global("ORGANISM", "hgnc")
		set_global("ASSAY_USE", "RNA")
		set_global("CELL_TYPE_COL", "scenic_celltype_run")
		set_global("SAMPLE_COL", cfg$sample_col)
		set_global("N_CORES", scenic_n_cores)
		set_global("SCENIC_DATASET_TITLE", sprintf("%s_%s_SCENIC_by_celltype_20260519", cfg$lineage, safe_ct))
		set_global("REDUCTION_CANDIDATES", cfg$reductions)
		set_global("EXPORT_FULL_AUC_MATRIX_CSV", parse_bool_env("SCENIC_EXPORT_FULL_AUC_MATRIX_CSV", FALSE))
		set_global("INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE", parse_int_env("SCENIC_INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE", 200L))
		set_global("INFERENCE_MAX_CELLS_PER_CELLTYPE", parse_int_env("SCENIC_INFERENCE_MAX_CELLS_PER_CELLTYPE", 3000L))
		set_global("INFERENCE_GLOBAL_MAX_CELLS", parse_int_env("SCENIC_INFERENCE_GLOBAL_MAX_CELLS", 30000L))
		set_global("MAX_INFERENCE_DENSE_GB", as.numeric(Sys.getenv("SCENIC_MAX_INFERENCE_DENSE_GB", unset = "8")))
		set_global("MIN_GENES_PER_CELL", parse_int_env("SCENIC_MIN_GENES_PER_CELL", 200L))
		set_global("MIN_CELLS_PER_GENE", parse_int_env("SCENIC_MIN_CELLS_PER_GENE", 10L))
		set_global("MIN_GENES_PER_REGULON", parse_int_env("SCENIC_MIN_GENES_PER_REGULON", 20L))
		set_global("SCENIC_GENIE3_TREE_METHOD", Sys.getenv("SCENIC_GENIE3_TREE_METHOD", unset = if (isTRUE(fast_modules)) "ET" else "RF"))
		set_global("SCENIC_GENIE3_NTREES", parse_int_env("SCENIC_GENIE3_NTREES", if (isTRUE(fast_modules)) 100L else 1000L))
		set_global("SCENIC_MODULE_WEIGHT_THRESHOLD", parse_numeric_vector_env("SCENIC_MODULE_WEIGHT_THRESHOLD", if (isTRUE(fast_modules)) 0.005 else NULL))
		set_global("SCENIC_MODULE_TOP_THR", parse_numeric_vector_env("SCENIC_MODULE_TOP_THR", if (isTRUE(fast_modules)) 0.005 else NULL))
		set_global("SCENIC_MODULE_N_TOP_TFS", parse_numeric_vector_env("SCENIC_MODULE_N_TOP_TFS", if (isTRUE(fast_modules)) 5 else NULL))
		set_global("SCENIC_MODULE_N_TOP_TARGETS", parse_numeric_vector_env("SCENIC_MODULE_N_TOP_TARGETS", if (isTRUE(fast_modules)) 50 else NULL))
		set_global("SCENIC_MODULE_CORR_THR", parse_numeric_vector_env("SCENIC_MODULE_CORR_THR", NULL))
		coex_methods <- parse_csv_env("SCENIC_REGULON_COEX_METHODS", if (isTRUE(fast_modules)) "top5perTarget" else character())
		set_global("SCENIC_REGULON_COEX_METHODS", if (length(coex_methods) > 0L) coex_methods else NULL)

		if (isTRUE(preflight_only)) {
			rm(obj_ct)
			gc(verbose = FALSE)
			list(status = "preflight_only")
		} else {
			sink(log_path, split = TRUE)
			on.exit(close_output_sinks(), add = TRUE)
			status_msg("[RUN] %s / %s cells=%d samples=%d output=%s\n", cfg$lineage, ct_label, n_cells, n_samples, output_dir)
			scenic_result <- run_scenic_module(obj_ct)
			rm(scenic_result, obj_ct)
			gc(verbose = FALSE)
			list(status = "ok")
		}
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
	status_msg("[%s] %s / %s elapsed_min=%.2f\n", toupper(base_status$status), cfg$lineage, ct_label, base_status$elapsed_min)
	as.data.table(base_status)
}

audit_rows <- list()
task_counter <- 0L

status_msg("[BOOT] SCENIC by-celltype runner started at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
status_msg("[BOOT] selected lineages: %s | n_cores=%d | force=%s | preflight_only=%s | fast_modules=%s\n", paste(selected_lineages, collapse = ","), scenic_n_cores, force_run, preflight_only, fast_modules)

for (lineage_name in selected_lineages) {
	cfg <- lineage_configs[[lineage_name]]
	status_msg("\n[LINEAGE] %s\n", lineage_name)
	status_msg("[LOAD] %s\n", cfg$rds_path)
	obj <- readRDS(cfg$rds_path)
	ct_summary <- list_celltypes(obj, cfg)
	if (length(selected_celltypes) > 0L) {
		ct_summary <- ct_summary[celltype %in% selected_celltypes]
	}
	if (nrow(ct_summary) == 0L) {
		status_msg("[WARN] no selected cell types for lineage %s\n", lineage_name)
		rm(obj); gc(verbose = FALSE)
		next
	}

	for (i in seq_len(nrow(ct_summary))) {
		if (max_celltypes > 0L && task_counter >= max_celltypes) break
		task_counter <- task_counter + 1L
		row <- ct_summary[i]
		audit_rows[[length(audit_rows) + 1L]] <- run_one_celltype(
			obj = obj,
			cfg = cfg,
			ct_label = row$celltype[[1]],
			n_cells = row$n_cells[[1]],
			n_samples = row$n_samples[[1]]
		)
		current <- rbindlist(audit_rows, fill = TRUE)
		fwrite(current, AUDIT_TSV, sep = "\t")
		if (identical(audit_rows[[length(audit_rows)]]$status[[1]], "error") && isTRUE(stop_after_first_error)) {
			stop(sprintf("Stopping after first SCENIC error: %s / %s", lineage_name, row$celltype[[1]]), call. = FALSE)
		}
	}
	rm(obj)
	gc(verbose = FALSE)
	if (max_celltypes > 0L && task_counter >= max_celltypes) break
}

audit <- if (length(audit_rows) > 0L) rbindlist(audit_rows, fill = TRUE) else data.table()
fwrite(audit, AUDIT_TSV, sep = "\t")
status_msg("\n[DONE] SCENIC by-celltype runner finished. tasks=%d audit=%s\n", nrow(audit), AUDIT_TSV)
if (nrow(audit) > 0L) print(audit[, .N, by = status][order(status)])
