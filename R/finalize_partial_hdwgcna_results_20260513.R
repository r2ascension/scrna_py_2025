#!/usr/bin/env Rscript

source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

RUN_ROOT <- Sys.getenv('PROGRAM_RUN_ROOT', unset = '/home/h2048/output/program_full_parallel_methods_20260507')
STATUS_TSV <- file.path(RUN_ROOT, 'run_status.tsv')
LINEAGES <- {
  x <- Sys.getenv('HDWGCNA_FINALIZE_LINEAGES', unset = 'bcell,stromal_fibroblast')
  trimws(strsplit(x, ',', fixed = TRUE)[[1]])
}

safe_id <- function(x) gsub('[^A-Za-z0-9_]', '_', as.character(x))
now_chr <- function() format(Sys.time(), '%Y-%m-%d %H:%M:%S')

write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = sub('\\.json$', '.rds', path))
  if (requireNamespace('jsonlite', quietly = TRUE)) {
    jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = 'null')
  }
  invisible(path)
}

read_json <- function(path) {
  if (!file.exists(path) || !requireNamespace('jsonlite', quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

append_status <- function(lineage, method, status, output_dir, note = NA_character_) {
  row <- data.frame(
    timestamp = now_chr(), lineage = lineage, method = method, status = status,
    output_dir = output_dir, note = as.character(note), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(row, file = STATUS_TSV, sep = '\t', row.names = FALSE,
                     col.names = !file.exists(STATUS_TSV), quote = FALSE,
                     append = file.exists(STATUS_TSV), na = '')
}

status_path <- function(ct_dir) file.path(ct_dir, 'hdwgcna_celltype_status.json')
write_ct_status <- function(ct_dir, rec) {
  rec$updated_at <- now_chr()
  rec$output_dir <- ct_dir
  write_json(rec, status_path(ct_dir))
}

infer_no_modules_from_log <- function(log_lines) {
  out <- character()
  m <- regexec('\\[WARN\\] ([^:]+): no non-grey modules detected', log_lines)
  hits <- regmatches(log_lines, m)
  for (h in hits) if (length(h) >= 2) out <- c(out, h[[2]])
  unique(out)
}

infer_errors_from_log <- function(log_lines) {
  out <- list()
  m <- regexec('\\[ERROR\\] ([^ ]+) failed: (.*) -- skipping\\.', log_lines)
  hits <- regmatches(log_lines, m)
  for (h in hits) {
    if (length(h) >= 3) out[[h[[2]]]] <- h[[3]]
  }
  out
}

infer_celltypes_from_worker_config <- function(lineage_dir) {
  cfg_path <- file.path(lineage_dir, 'worker_configs', 'hdwgcna_worker_config.json')
  cfg <- read_json(cfg_path)
  if (is.null(cfg) || is.null(cfg$rds_path) || is.null(cfg$celltype_col) || !file.exists(cfg$rds_path)) {
    return(character())
  }
  out <- tryCatch({
    suppressPackageStartupMessages(library(Seurat))
    obj <- readRDS(cfg$rds_path)
    vals <- sort(unique(as.character(obj@meta.data[[cfg$celltype_col]])))
    rm(obj)
    vals[nzchar(vals) & !is.na(vals)]
  }, error = function(e) {
    cat(sprintf('[WARN] Could not infer celltypes from %s: %s\n', cfg_path, conditionMessage(e)))
    character()
  })
  out
}

finalize_lineage <- function(lineage) {
  lineage_dir <- file.path(RUN_ROOT, lineage)
  out_dir <- file.path(lineage_dir, 'hdwgcna_full')
  if (!dir.exists(out_dir)) {
    cat(sprintf('[SKIP] %s: no hdwgcna_full dir\n', lineage))
    return(invisible(NULL))
  }
  log_path <- file.path(lineage_dir, 'worker_logs', 'hdwgcna_worker.log')
  log_lines <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
  no_module_ct <- infer_no_modules_from_log(log_lines)
  error_map <- infer_errors_from_log(log_lines)
  all_ct <- infer_celltypes_from_worker_config(lineage_dir)

  # Existing celltype dirs, plus log-only celltypes.
  ct_dirs <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)
  ct_dirs <- ct_dirs[grepl('/hdwgcna_', ct_dirs)]
  ct_names <- sub('^hdwgcna_', '', basename(ct_dirs))
  names(ct_dirs) <- ct_names
  for (ct in unique(c(all_ct, no_module_ct, names(error_map)))) {
    ct_safe <- safe_id(ct)
    if (!ct_safe %in% names(ct_dirs)) {
      ct_dirs[[ct_safe]] <- file.path(out_dir, paste0('hdwgcna_', ct_safe))
      dir.create(ct_dirs[[ct_safe]], recursive = TRUE, showWarnings = FALSE)
    }
  }

  results <- list()
  for (ct_safe in names(ct_dirs)) {
    ct_dir <- ct_dirs[[ct_safe]]
    ct <- ct_safe
    rec <- read_json(status_path(ct_dir))
    membership <- file.path(ct_dir, 'hdwgcna_module_membership.csv')
    soft_power <- file.path(ct_dir, 'hdwgcna_soft_power.pdf')

    if (!is.null(rec) && !is.null(rec$status)) {
      status <- as.character(rec$status)[1]
      if (identical(status, 'started')) {
        rec$status <- 'timeout'
        rec$error <- 'stale started status from interrupted worker; finalized as timeout_skipped'
        rec$finalized_from <- 'stale_started'
        write_ct_status(ct_dir, rec)
      }
      results[[ct]] <- rec
      next
    }

    if (file.exists(membership)) {
      mods <- tryCatch(utils::read.csv(membership, stringsAsFactors = FALSE), error = function(e) NULL)
      module_ids <- if (!is.null(mods) && 'module' %in% colnames(mods)) setdiff(unique(as.character(mods$module)), 'grey') else character()
      rec <- list(celltype = ct, status = 'ok', module_ids = module_ids, n_modules = length(module_ids), membership_csv = membership, finalized_from = 'existing_membership')
    } else if (ct %in% safe_id(no_module_ct) || ct %in% no_module_ct) {
      rec <- list(celltype = ct, status = 'no_modules', module_ids = character(), n_modules = 0L, finalized_from = 'worker_log_no_modules')
    } else if (ct %in% names(error_map) || ct_safe %in% safe_id(names(error_map))) {
      err_ct <- names(error_map)[safe_id(names(error_map)) == ct_safe][1]
      rec <- list(celltype = ct, status = 'error', error = unname(error_map[[err_ct]]), finalized_from = 'worker_log_error')
    } else if (file.exists(soft_power)) {
      rec <- list(celltype = ct, status = 'timeout', error = 'soft-power exists but no final status/membership; finalized as interrupted_or_timeout', finalized_from = 'soft_power_only')
    } else {
      rec <- list(celltype = ct, status = 'not_started', finalized_from = 'directory_only')
    }
    write_ct_status(ct_dir, rec)
    results[[ct]] <- rec
  }

  imported <- pa_import_hdwgcna_results(out_dir, unit_id = paste0(lineage, '_hdwgcna_partial_20260513'))
  status_vec <- vapply(results, function(x) as.character(x$status)[1], character(1))
  status_counts <- as.list(table(status_vec))
  partial <- any(!(status_vec %in% c('ok', 'no_modules')))
  runner <- c(imported, list(
    status = 'ok',
    method_status = if (partial) 'partial_with_timeouts_or_errors' else 'ok',
    partial = partial,
    status_counts = status_counts,
    hdwgcna_results = results,
    finalized_at = now_chr(),
    finalized_by = normalizePath(sys.frame(1)$ofile %||% '/home/h2048/script/R/finalize_partial_hdwgcna_results_20260513.R', mustWork = FALSE)
  ))
  saveRDS(runner, file.path(out_dir, 'hdwgcna_runner_result.rds'))
  worker_result <- list(
    status = 'ok',
    method_status = runner$method_status,
    partial = partial,
    lineage = lineage,
    method = 'hdwgcna',
    output_dir = out_dir,
    status_counts = status_counts,
    finalized_at = runner$finalized_at,
    result_rds = file.path(out_dir, 'hdwgcna_runner_result.rds')
  )
  write_json(worker_result, file.path(out_dir, 'method_worker_result.json'))
  append_status(lineage, 'hdwgcna', 'partial_finalized', out_dir, sprintf('method_status=%s; status_counts=%s', runner$method_status, paste(names(status_counts), unlist(status_counts), sep=':', collapse=';')))
  cat(sprintf('[OK] %s finalized: %s\n', lineage, paste(names(status_counts), unlist(status_counts), sep='=', collapse=', ')))
  invisible(worker_result)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
for (lineage in LINEAGES) finalize_lineage(lineage)
