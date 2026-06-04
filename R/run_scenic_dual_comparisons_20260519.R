#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Run R-SCENIC at two comparison levels:
#   1) lineage_celltype : within each major lineage, compare cell_type_L3 groups
#   2) celltype_tissue  : within each cell_type_L3, compare tissue-origin groups
# ==============================================================================

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
SCENIC_ROOT <- file.path(RUN_ROOT, "scenic_dual_comparisons")
AUDIT_TSV <- Sys.getenv(
	"SCENIC_AUDIT_TSV",
	unset = file.path(SCENIC_ROOT, "scenic_dual_comparison_audit_20260519.tsv")
)
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

lineage_configs <- list(
	bcell = list(
		lineage = "bcell",
		rds_path = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/bcell_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_refined", "umap_scanvi", "umap", "harmony", "pca")
	),
	stromal_smc = list(
		lineage = "stromal_smc",
		rds_path = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_scanvi", "umap", "harmony", "pca")
	),
	stromal_fibroblast = list(
		lineage = "stromal_fibroblast",
		rds_path = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_scanvi", "umap", "harmony", "pca")
	),
	stromal_endothelial = list(
		lineage = "stromal_endothelial",
		rds_path = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_endothelial_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_scanvi", "umap", "harmony", "pca")
	),
	tnk = list(
		lineage = "tnk",
		rds_path = "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_refined", "umap_scanvi", "umap_scanvi_corrected", "umap_scvi", "umap", "harmony", "pca")
	),
	myeloid = list(
		lineage = "myeloid",
		rds_path = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_refined", "umap_scanvi", "umap", "harmony", "pca")
	),
	epithelial = list(
		lineage = "epithelial",
		rds_path = "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.rds",
		celltype_col = "cell_type_L3",
		sample_col = "sample",
		tissue_col = "tissue",
		reductions = c("umap_refined", "umap_scanvi", "umap", "harmony", "pca")
	)
)

selected_lineages <- parse_csv_env("SCENIC_LINEAGES", names(lineage_configs))
unknown_lineages <- setdiff(selected_lineages, names(lineage_configs))
if (length(unknown_lineages) > 0L) {
	stop(sprintf("Unknown SCENIC_LINEAGES: %s", paste(unknown_lineages, collapse = ", ")), call. = FALSE)
}
selected_celltypes <- parse_csv_env("SCENIC_CELLTYPES", character())
comparison_levels <- parse_csv_env("SCENIC_COMPARISON_LEVELS", c("lineage_celltype", "celltype_tissue"))
unknown_levels <- setdiff(comparison_levels, c("lineage_celltype", "celltype_tissue"))
if (length(unknown_levels) > 0L) {
	stop(sprintf("Unknown SCENIC_COMPARISON_LEVELS: %s", paste(unknown_levels, collapse = ", ")), call. = FALSE)
}

force_run <- parse_bool_env("SCENIC_FORCE", FALSE)
preflight_only <- parse_bool_env("SCENIC_PREFLIGHT_ONLY", FALSE)
stop_after_first_error <- parse_bool_env("SCENIC_STOP_AFTER_FIRST_ERROR", FALSE)
fast_modules <- parse_bool_env("SCENIC_FAST_MODULES", TRUE)
min_cells <- parse_int_env("SCENIC_MIN_CELLS", 50L)
min_groups <- parse_int_env("SCENIC_MIN_GROUPS", 2L)
min_group_cells <- parse_int_env("SCENIC_MIN_GROUP_CELLS", 20L)
max_jobs <- parse_int_env("SCENIC_MAX_JOBS", 0L)
scenic_n_cores <- parse_int_env("SCENIC_N_CORES", 4L)
heavy_lineages <- parse_csv_env("SCENIC_HEAVY_LINEAGES", c("myeloid", "epithelial"))

dir.create(SCENIC_ROOT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(CORE_SCRIPT)) stop(sprintf("Missing SCENIC core script: %s", CORE_SCRIPT), call. = FALSE)
if (!dir.exists(DATABASE_DIR_DEFAULT)) stop(sprintf("Missing SCENIC database dir: %s", DATABASE_DIR_DEFAULT), call. = FALSE)
if (!file.exists(SCENIC_DB_10KB_DEFAULT) || !file.exists(SCENIC_DB_500BP_DEFAULT)) {
	stop("Missing compatible R-SCENIC feather database files.", call. = FALSE)
}

SCENIC_SOURCE_ONLY <- TRUE
source(CORE_SCRIPT)

validate_cfg_columns <- function(obj, cfg) {
	missing_cols <- setdiff(c(cfg$celltype_col, cfg$sample_col, cfg$tissue_col), colnames(obj@meta.data))
	if (length(missing_cols) > 0L) {
		stop(sprintf("%s missing metadata columns: %s", cfg$lineage, paste(missing_cols, collapse = ", ")), call. = FALSE)
	}
}

normalize_group <- function(x, unknown = "Unknown") {
	x <- as.character(x)
	x[is.na(x) | !nzchar(trimws(x))] <- unknown
	x
}

list_celltypes <- function(obj, cfg) {
	meta <- obj@meta.data
	dt <- data.table(
		cell = rownames(meta),
		celltype = normalize_group(meta[[cfg$celltype_col]], unknown = "Unknown_celltype"),
		sample = normalize_group(meta[[cfg$sample_col]], unknown = "Unknown_sample"),
		tissue = normalize_group(meta[[cfg$tissue_col]], unknown = "Unknown_tissue")
	)
	dt[celltype != "Unknown_celltype", .(
		n_cells = .N,
		n_samples = uniqueN(sample),
		n_tissues = uniqueN(tissue)
	), by = celltype][order(-n_cells, celltype)]
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

write_job_metadata <- function(obj_job, cfg, output_dir, comparison_level, unit_label, group_col = "scenic_comparison_group") {
	dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
	meta <- obj_job@meta.data
	out <- data.table(
		cell = rownames(meta),
		lineage = cfg$lineage,
		comparison_level = comparison_level,
		unit_label = unit_label,
		cell_type_L3 = normalize_group(meta[[cfg$celltype_col]], unknown = "Unknown_celltype"),
		comparison_group = normalize_group(meta[[group_col]], unknown = "Unknown_group"),
		sample = normalize_group(meta[[cfg$sample_col]], unknown = "Unknown_sample"),
		tissue = normalize_group(meta[[cfg$tissue_col]], unknown = "Unknown_tissue")
	)
	fwrite(out, file.path(output_dir, "tables", "cell_metadata_used_for_scenic_dual_comparison.csv"))
	fwrite(out[, .N, by = .(lineage, comparison_level, unit_label, comparison_group)][order(-N)], file.path(output_dir, "tables", "comparison_group_summary.csv"))
	fwrite(out[, .N, by = .(lineage, comparison_level, unit_label, cell_type_L3, tissue)][order(cell_type_L3, tissue)], file.path(output_dir, "tables", "celltype_tissue_summary.csv"))
	fwrite(out[, .N, by = .(lineage, comparison_level, unit_label, sample, tissue, comparison_group)][order(tissue, sample, comparison_group)], file.path(output_dir, "tables", "sample_tissue_group_summary.csv"))
	invisible(out)
}

infer_downsampling_defaults <- function(cfg, comparison_level) {
	is_heavy <- cfg$lineage %in% heavy_lineages
	if (isTRUE(is_heavy) && identical(comparison_level, "lineage_celltype")) {
		return(list(per_sample_group = 40L, per_group = 400L, global = 4000L, max_genes = 3500L))
	}
	if (isTRUE(is_heavy) && identical(comparison_level, "celltype_tissue")) {
		return(list(per_sample_group = 60L, per_group = 800L, global = 4000L, max_genes = 3000L))
	}
	list(per_sample_group = 200L, per_group = 3000L, global = 30000L, max_genes = 0L)
}

configure_core_globals <- function(cfg, output_dir, dataset_title, comparison_level) {
	downsample_defaults <- infer_downsampling_defaults(cfg, comparison_level)
	set_global("INPUT_RDS", cfg$rds_path)
	set_global("OUTPUT_DIR", output_dir)
	set_global("DATABASE_DIR", DATABASE_DIR_DEFAULT)
	set_global("SCENIC_DB_10KB", SCENIC_DB_10KB_DEFAULT)
	set_global("SCENIC_DB_500BP", SCENIC_DB_500BP_DEFAULT)
	set_global("SCENIC_DB_SCOPE", Sys.getenv("SCENIC_DB_SCOPE", unset = "500bp"))
	set_global("SCENIC_DB_INDEX_COL", NULL)
	set_global("ORGANISM", "hgnc")
	set_global("ASSAY_USE", "RNA")
	set_global("CELL_TYPE_COL", "scenic_comparison_group")
	set_global("SAMPLE_COL", cfg$sample_col)
	set_global("N_CORES", scenic_n_cores)
	set_global("SCENIC_DATASET_TITLE", dataset_title)
	set_global("REDUCTION_CANDIDATES", cfg$reductions)
	set_global("EXPORT_FULL_AUC_MATRIX_CSV", parse_bool_env("SCENIC_EXPORT_FULL_AUC_MATRIX_CSV", FALSE))
	set_global("INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE", parse_int_env("SCENIC_INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE", downsample_defaults$per_sample_group))
	set_global("INFERENCE_MAX_CELLS_PER_CELLTYPE", parse_int_env("SCENIC_INFERENCE_MAX_CELLS_PER_CELLTYPE", downsample_defaults$per_group))
	set_global("INFERENCE_GLOBAL_MAX_CELLS", parse_int_env("SCENIC_INFERENCE_GLOBAL_MAX_CELLS", downsample_defaults$global))
	set_global("INFERENCE_MAX_GENES", parse_int_env("SCENIC_INFERENCE_MAX_GENES", downsample_defaults$max_genes))
	set_global("MAX_INFERENCE_DENSE_GB", as.numeric(Sys.getenv("SCENIC_MAX_INFERENCE_DENSE_GB", unset = "8")))
	set_global("MIN_GENES_PER_CELL", parse_int_env("SCENIC_MIN_GENES_PER_CELL", 200L))
	set_global("MIN_CELLS_PER_GENE", parse_int_env("SCENIC_MIN_CELLS_PER_GENE", 10L))
	set_global("MIN_GENE_PCT", as.numeric(Sys.getenv("SCENIC_MIN_GENE_PCT", unset = "0.01")))
	set_global("MIN_GENES_PER_REGULON", parse_int_env("SCENIC_MIN_GENES_PER_REGULON", 20L))
	set_global("SCENIC_GENIE3_TREE_METHOD", Sys.getenv("SCENIC_GENIE3_TREE_METHOD", unset = if (isTRUE(fast_modules)) "ET" else "RF"))
	set_global("SCENIC_GENIE3_NTREES", parse_int_env("SCENIC_GENIE3_NTREES", if (isTRUE(fast_modules)) 100L else 1000L))
	set_global("SCENIC_MODULE_WEIGHT_THRESHOLD", parse_numeric_vector_env("SCENIC_MODULE_WEIGHT_THRESHOLD", if (isTRUE(fast_modules)) 0.01 else NULL))
	set_global("SCENIC_MODULE_TOP_THR", parse_numeric_vector_env("SCENIC_MODULE_TOP_THR", if (isTRUE(fast_modules)) 0.002 else NULL))
	set_global("SCENIC_MODULE_N_TOP_TFS", parse_numeric_vector_env("SCENIC_MODULE_N_TOP_TFS", if (isTRUE(fast_modules)) 3 else NULL))
	set_global("SCENIC_MODULE_N_TOP_TARGETS", parse_numeric_vector_env("SCENIC_MODULE_N_TOP_TARGETS", if (isTRUE(fast_modules)) 30 else NULL))
	set_global("SCENIC_MODULE_CORR_THR", parse_numeric_vector_env("SCENIC_MODULE_CORR_THR", NULL))
	coex_methods <- parse_csv_env("SCENIC_REGULON_COEX_METHODS", if (isTRUE(fast_modules)) "top3perTarget" else character())
	set_global("SCENIC_REGULON_COEX_METHODS", if (length(coex_methods) > 0L) coex_methods else NULL)
	status_msg(
		"[CONFIG] %s / %s inference caps: per_sample_group=%d per_group=%d global=%d max_genes=%d heavy_lineage=%s\n",
		cfg$lineage,
		comparison_level,
		INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE,
		INFERENCE_MAX_CELLS_PER_CELLTYPE,
		INFERENCE_GLOBAL_MAX_CELLS,
		INFERENCE_MAX_GENES,
		cfg$lineage %in% heavy_lineages
	)
}

run_one_job <- function(obj_job, cfg, comparison_level, unit_label, output_dir) {
	status_path <- file.path(output_dir, "scenic_dual_status.json")
	log_path <- file.path(output_dir, "scenic_dual_wrapper.log")
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

	n_cells <- ncol(obj_job)
	groups <- normalize_group(obj_job@meta.data$scenic_comparison_group, unknown = "Unknown_group")
	samples <- normalize_group(obj_job@meta.data[[cfg$sample_col]], unknown = "Unknown_sample")
	tissues <- normalize_group(obj_job@meta.data[[cfg$tissue_col]], unknown = "Unknown_tissue")
	group_sizes <- table(groups)
	n_groups <- uniqueN(groups)
	n_groups_min_cells <- sum(group_sizes >= min_group_cells)
	n_samples <- uniqueN(samples)
	n_tissues <- uniqueN(tissues)

	existing <- read_status(status_path)
	if (!isTRUE(force_run) && !is.null(existing) && identical(existing$status, "ok")) {
		status_msg("[SKIP] %s / %s / %s existing ok\n", comparison_level, cfg$lineage, unit_label)
		return(as.data.table(existing))
	}

	started <- Sys.time()
	base_status <- list(
		status = "started",
		comparison_level = comparison_level,
		lineage = cfg$lineage,
		unit_label = unit_label,
		safe_unit_label = safe_file_id(unit_label),
		grouping_col = "scenic_comparison_group",
		n_cells = as.integer(n_cells),
		n_groups = as.integer(n_groups),
		n_groups_min_cells = as.integer(n_groups_min_cells),
		n_samples = as.integer(n_samples),
		n_tissues = as.integer(n_tissues),
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

	if (n_cells < min_cells) {
		base_status$status <- "skipped_too_few_cells"
		base_status$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
		write_json_safe(base_status, status_path)
		status_msg("[SKIP] %s / %s / %s too few cells: %d < %d\n", comparison_level, cfg$lineage, unit_label, n_cells, min_cells)
		return(as.data.table(base_status))
	}
	if (n_groups_min_cells < min_groups) {
		base_status$status <- "skipped_too_few_groups"
		base_status$ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
		write_json_safe(base_status, status_path)
		status_msg("[SKIP] %s / %s / %s too few groups with >=%d cells: %d < %d\n", comparison_level, cfg$lineage, unit_label, min_group_cells, n_groups_min_cells, min_groups)
		return(as.data.table(base_status))
	}

	write_job_metadata(obj_job, cfg, output_dir, comparison_level, unit_label)

	result_status <- tryCatch({
		configure_core_globals(
			cfg = cfg,
			output_dir = output_dir,
			dataset_title = sprintf("%s_%s_%s_SCENIC_20260519", comparison_level, cfg$lineage, safe_file_id(unit_label)),
			comparison_level = comparison_level
		)
		if (isTRUE(preflight_only)) {
			list(status = "preflight_only")
		} else {
			sink(log_path, split = TRUE)
			on.exit(close_output_sinks(), add = TRUE)
			status_msg("[RUN] %s / %s / %s cells=%d groups=%d tissues=%d output=%s\n", comparison_level, cfg$lineage, unit_label, n_cells, n_groups, n_tissues, output_dir)
			scenic_result <- run_scenic_module(obj_job)
			rm(scenic_result)
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
	status_msg("[%s] %s / %s / %s elapsed_min=%.2f\n", toupper(base_status$status), comparison_level, cfg$lineage, unit_label, base_status$elapsed_min)
	as.data.table(base_status)
}

make_lineage_celltype_job <- function(obj, cfg) {
	obj_job <- obj
	obj_job@meta.data$scenic_comparison_group <- normalize_group(obj_job@meta.data[[cfg$celltype_col]], unknown = "Unknown_celltype")
	if (length(selected_celltypes) > 0L) {
		cells <- rownames(obj_job@meta.data)[obj_job@meta.data$scenic_comparison_group %in% selected_celltypes]
		obj_job <- subset(obj_job, cells = cells)
	}
	obj_job
}

make_celltype_tissue_job <- function(obj, cfg, ct_label) {
	meta <- obj@meta.data
	cells <- rownames(meta)[normalize_group(meta[[cfg$celltype_col]], unknown = "Unknown_celltype") == ct_label]
	if (length(cells) == 0L) return(NULL)
	obj_ct <- subset(obj, cells = cells)
	tissue_vals <- normalize_group(obj_ct@meta.data[[cfg$tissue_col]], unknown = "Unknown_tissue")
	group_dt <- data.table(cell = colnames(obj_ct), tissue = tissue_vals)
	keep_groups <- group_dt[, .N, by = tissue][N >= min_group_cells, tissue]
	if (length(keep_groups) > 0L) {
		keep_cells <- group_dt[tissue %in% keep_groups, cell]
		obj_ct <- subset(obj_ct, cells = keep_cells)
		tissue_vals <- normalize_group(obj_ct@meta.data[[cfg$tissue_col]], unknown = "Unknown_tissue")
	}
	obj_ct@meta.data$scenic_comparison_group <- tissue_vals
	obj_ct
}

audit_rows <- list()
job_counter <- 0L

status_msg("[BOOT] SCENIC dual-comparison runner started at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
status_msg("[BOOT] selected lineages: %s | levels=%s | n_cores=%d | force=%s | preflight_only=%s | fast_modules=%s | min_group_cells=%d\n",
	paste(selected_lineages, collapse = ","), paste(comparison_levels, collapse = ","), scenic_n_cores, force_run, preflight_only, fast_modules, min_group_cells)

append_audit <- function(row) {
	audit_rows[[length(audit_rows) + 1L]] <<- row
	current <- rbindlist(audit_rows, fill = TRUE)
	fwrite(current, AUDIT_TSV, sep = "\t")
	if (identical(row$status[[1]], "error") && isTRUE(stop_after_first_error)) {
		stop(sprintf("Stopping after first SCENIC error: %s / %s / %s", row$comparison_level[[1]], row$lineage[[1]], row$unit_label[[1]]), call. = FALSE)
	}
}

for (lineage_name in selected_lineages) {
	cfg <- lineage_configs[[lineage_name]]
	status_msg("\n[LINEAGE] %s\n", lineage_name)
	status_msg("[LOAD] %s\n", cfg$rds_path)
	obj <- readRDS(cfg$rds_path)
	validate_cfg_columns(obj, cfg)

	if ("lineage_celltype" %in% comparison_levels) {
		obj_lineage <- make_lineage_celltype_job(obj, cfg)
		output_dir <- file.path(SCENIC_ROOT, "lineage_celltype", cfg$lineage, paste0("scenic_", safe_file_id(cfg$lineage), "_cell_type_L3"))
		job_counter <- job_counter + 1L
		append_audit(run_one_job(obj_lineage, cfg, "lineage_celltype", paste0(cfg$lineage, "__cell_type_L3"), output_dir))
		rm(obj_lineage); gc(verbose = FALSE)
		if (max_jobs > 0L && job_counter >= max_jobs) break
	}

	if ("celltype_tissue" %in% comparison_levels && !(max_jobs > 0L && job_counter >= max_jobs)) {
		ct_summary <- list_celltypes(obj, cfg)
		if (length(selected_celltypes) > 0L) {
			ct_summary <- ct_summary[celltype %in% selected_celltypes]
		}
		if (nrow(ct_summary) == 0L) {
			status_msg("[WARN] no selected cell types for tissue comparison in lineage %s\n", lineage_name)
		} else {
			for (i in seq_len(nrow(ct_summary))) {
				if (max_jobs > 0L && job_counter >= max_jobs) break
				ct_label <- ct_summary$celltype[[i]]
				obj_ct <- make_celltype_tissue_job(obj, cfg, ct_label)
				if (is.null(obj_ct)) next
				output_dir <- file.path(SCENIC_ROOT, "celltype_tissue", cfg$lineage, paste0("scenic_", safe_file_id(ct_label)))
				job_counter <- job_counter + 1L
				append_audit(run_one_job(obj_ct, cfg, "celltype_tissue", ct_label, output_dir))
				rm(obj_ct); gc(verbose = FALSE)
			}
		}
	}

	rm(obj)
	gc(verbose = FALSE)
	if (max_jobs > 0L && job_counter >= max_jobs) break
}

audit <- if (length(audit_rows) > 0L) rbindlist(audit_rows, fill = TRUE) else data.table()
fwrite(audit, AUDIT_TSV, sep = "\t")
status_msg("\n[DONE] SCENIC dual-comparison runner finished. tasks=%d audit=%s\n", nrow(audit), AUDIT_TSV)
if (nrow(audit) > 0L) print(audit[, .N, by = status][order(status)])
