#!/usr/bin/env Rscript

source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

RUN_ROOT <- pa_prepare_output_dir('/home/h2048/output/program_full_parallel_methods_20260507')
STATUS_TSV <- file.path(RUN_ROOT, 'run_status.tsv')
MANIFEST_JSON <- file.path(RUN_ROOT, 'lineage_manifest.json')
WORKER_R <- '/home/h2048/script/R/program_lineage_method_worker_20260507.R'
MONITOR_SH <- '/home/h2048/script/bash/monitor_process_tree_20260507.sh'

CNMF_PYTHON <- '/home/h2048/miniconda3/envs/bbknn_env/bin/python'
HARMONY2_PYTHON <- '/home/h2048/miniconda3/envs/bbknn_env/bin/python'
HARMONY2_HELPER <- '/home/h2048/script/py/harmony2_helper_20260511_v1.py'
PYCOGAPS_PYTHON <- '/home/h2048/miniconda3/envs/scarches_stable_pertpy/bin/python'
PYCOGAPS_HELPER <- '/home/h2048/script/py/pycogaps_helper_20260505_v1.py'
RUN_STAMP <- '20260507_parallel_methods'

COMPUTE_METHODS <- c('hdwgcna', 'covarnet', 'cnmf', 'pycogaps', 'harmony2')
LLM_METHODS <- paste0('llm_', COMPUTE_METHODS)
POLL_SECONDS <- 60L
MONITOR_INTERVAL_DEFAULT <- 60L
MONITOR_INTERVAL_CNMF <- 30L
SKIP_COMPLETED_METHODS <- Sys.getenv('SKIP_COMPLETED_METHODS', unset = '0') %in% c('1', 'true', 'TRUE', 'yes', 'YES')
HDWGCNA_RESUME_CELLTYPES <- !(Sys.getenv('HDWGCNA_RESUME_CELLTYPES', unset = '1') %in% c('0', 'false', 'FALSE', 'no', 'NO'))
HDWGCNA_CELLTYPE_TIMEOUT_SEC <- suppressWarnings(as.integer(Sys.getenv('HDWGCNA_CELLTYPE_TIMEOUT_SEC', unset = '21600')))
if (length(HDWGCNA_CELLTYPE_TIMEOUT_SEC) == 0L || is.na(HDWGCNA_CELLTYPE_TIMEOUT_SEC) || HDWGCNA_CELLTYPE_TIMEOUT_SEC < 0L) {
  HDWGCNA_CELLTYPE_TIMEOUT_SEC <- 21600L
}

lineage_configs <- list(
  list(
    lineage = 'bcell',
    rds_path = '/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/bcell_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/bcell_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  ),
  list(
    lineage = 'stromal_smc',
    rds_path = '/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  ),
  list(
    lineage = 'stromal_fibroblast',
    rds_path = '/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  ),
  list(
    lineage = 'stromal_endothelial',
    rds_path = '/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_endothelial_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_endothelial_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  ),
  list(
    lineage = 'tnk',
    rds_path = '/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  ),
  list(
    lineage = 'myeloid',
    rds_path = '/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  ),
  list(
    lineage = 'epithelial',
    rds_path = '/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.rds',
    h5ad_path = '/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.h5ad',
    celltype_col = 'cell_type_L3',
    sample_col = 'sample',
    batch_col = 'dataset',
    condition_col = 'condition',
    tissue_col = 'tissue'
  )
)

select_lineage_window <- function(configs) {
  start_lineage <- Sys.getenv('START_LINEAGE', unset = '')
  end_lineage <- Sys.getenv('END_LINEAGE', unset = '')
  lineage_names <- vapply(configs, function(cfg) cfg$lineage, character(1))

  start_idx <- 1L
  end_idx <- length(configs)
  if (nzchar(start_lineage)) {
    start_idx <- match(start_lineage, lineage_names)
    if (is.na(start_idx)) {
      stop(sprintf('START_LINEAGE=%s not found. Available lineages: %s', start_lineage, paste(lineage_names, collapse = ', ')), call. = FALSE)
    }
  }
  if (nzchar(end_lineage)) {
    end_idx <- match(end_lineage, lineage_names)
    if (is.na(end_idx)) {
      stop(sprintf('END_LINEAGE=%s not found. Available lineages: %s', end_lineage, paste(lineage_names, collapse = ', ')), call. = FALSE)
    }
  }
  if (start_idx > end_idx) {
    stop(sprintf('START_LINEAGE (%s) occurs after END_LINEAGE (%s)', lineage_names[[start_idx]], lineage_names[[end_idx]]), call. = FALSE)
  }
  selected <- configs[start_idx:end_idx]
  attr(selected, 'start_lineage') <- if (nzchar(start_lineage)) start_lineage else lineage_names[[start_idx]]
  attr(selected, 'end_lineage') <- if (nzchar(end_lineage)) end_lineage else lineage_names[[end_idx]]
  selected
}

lineage_configs <- select_lineage_window(lineage_configs)

status_lock_eval <- function(expr) {
  lock_dir <- paste0(STATUS_TSV, '.lock')
  acquired <- FALSE
  for (i in seq_len(1200L)) {
    acquired <- dir.create(lock_dir, showWarnings = FALSE)
    if (isTRUE(acquired)) break
    Sys.sleep(0.1)
  }
  if (isTRUE(acquired)) {
    on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  } else {
    warning(sprintf('Could not acquire status lock: %s', lock_dir), call. = FALSE)
  }
  force(expr)
}

append_status <- function(lineage, method, status, output_dir, note = NA_character_) {
  row <- data.frame(
    timestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
    lineage = as.character(lineage),
    method = as.character(method),
    status = as.character(status),
    output_dir = as.character(output_dir),
    note = as.character(note),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  status_lock_eval({
    utils::write.table(
      row,
      file = STATUS_TSV,
      sep = '\t',
      row.names = FALSE,
      col.names = !file.exists(STATUS_TSV),
      quote = FALSE,
      append = file.exists(STATUS_TSV),
      na = ''
    )
  })
}

write_banner <- function(...) {
  cat('\n', paste0(rep('=', 88), collapse = ''), '\n', sep = '')
  cat(sprintf(...), '\n')
  cat(paste0(rep('=', 88), collapse = ''), '\n', sep = '')
}

is_pid_running <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (length(pid) == 0L || is.na(pid) || pid <= 0L) return(FALSE)
  status <- suppressWarnings(system2('kill', args = c('-0', as.character(pid)), stdout = FALSE, stderr = FALSE))
  identical(status, 0L)
}

method_output_subdir <- function(method) {
  switch(
    method,
    hdwgcna = 'hdwgcna_full',
    covarnet = 'covarnet_full',
    cnmf = 'cnmf_full',
    pycogaps = 'pycogaps_full',
    harmony2 = 'harmony2_full',
    if (grepl('^llm_', method)) file.path('llm_parallel', sub('^llm_', '', method)) else method
  )
}

read_monitor_peak <- function(path) {
  if (!file.exists(path)) return(list(peak_rss_kb = NA_real_, min_available_kb = NA_real_))
  df <- tryCatch(utils::read.delim(path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(list(peak_rss_kb = NA_real_, min_available_kb = NA_real_))
  list(
    peak_rss_kb = suppressWarnings(max(as.numeric(df$rss_kb_sum), na.rm = TRUE)),
    min_available_kb = suppressWarnings(min(as.numeric(df$mem_available_kb), na.rm = TRUE))
  )
}

launch_monitor <- function(job) {
  interval <- if (identical(job$method, 'cnmf')) MONITOR_INTERVAL_CNMF else MONITOR_INTERVAL_DEFAULT
  cmd <- sprintf(
    '%s %s %s %s %s > %s 2>&1 & echo $!',
    shQuote(MONITOR_SH),
    shQuote(as.character(job$pid)),
    shQuote(job$monitor_tsv),
    shQuote(as.character(interval)),
    shQuote(sprintf('%s/%s', job$lineage, job$method)),
    shQuote(job$monitor_log)
  )
  monitor_pid <- suppressWarnings(as.integer(system(cmd, intern = TRUE)[1]))
  job$monitor_pid <- monitor_pid
  job
}

launch_method_worker <- function(cfg, method, lineage_dir) {
  output_dir <- pa_prepare_output_dir(file.path(lineage_dir, method_output_subdir(method)))
  config_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'worker_configs'))
  logs_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'worker_logs'))
  monitor_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'resource_monitor'))
  worker_config_path <- file.path(config_dir, sprintf('%s_worker_config.json', method))
  worker_log_path <- file.path(logs_dir, sprintf('%s_worker.log', method))
  monitor_tsv <- file.path(monitor_dir, sprintf('%s_memory_monitor.tsv', method))
  monitor_log <- file.path(monitor_dir, sprintf('%s_memory_monitor.log', method))

  worker_config <- modifyList(cfg, list(
    method = method,
    lineage_dir = lineage_dir,
    output_dir = output_dir,
    status_tsv = STATUS_TSV,
    run_root = RUN_ROOT,
    run_stamp = RUN_STAMP,
    cnmf_python = CNMF_PYTHON,
    harmony2_python = HARMONY2_PYTHON,
    harmony2_helper = HARMONY2_HELPER,
    pycogaps_python = PYCOGAPS_PYTHON,
    pycogaps_helper = PYCOGAPS_HELPER,
    hdwgcna_runner_args = list(
      resume_celltypes = HDWGCNA_RESUME_CELLTYPES,
      resume_skip_statuses = c('ok', 'no_modules'),
      celltype_timeout_sec = HDWGCNA_CELLTYPE_TIMEOUT_SEC
    ),
    llm_model = 'deepseek-reasoner',
    llm_timeout_sec = 240L
  ))
  pa_write_json(worker_config, worker_config_path)

  cmd <- sprintf(
    'env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 /usr/bin/Rscript %s %s > %s 2>&1 & echo $!',
    shQuote(WORKER_R),
    shQuote(worker_config_path),
    shQuote(worker_log_path)
  )
  pid <- suppressWarnings(as.integer(system(cmd, intern = TRUE)[1]))
  job <- list(
    lineage = cfg$lineage,
    method = method,
    pid = pid,
    output_dir = output_dir,
    config_path = worker_config_path,
    log_path = worker_log_path,
    result_path = file.path(output_dir, 'method_worker_result.json'),
    monitor_tsv = monitor_tsv,
    monitor_log = monitor_log,
    monitor_pid = NA_integer_
  )
  job <- launch_monitor(job)
  append_status(
    cfg$lineage,
    method,
    'launched_async',
    output_dir,
    sprintf('pid=%s; monitor_pid=%s; log=%s; monitor=%s', job$pid, job$monitor_pid, worker_log_path, monitor_tsv)
  )
  cat(sprintf('[LAUNCH] %s / %s pid=%s monitor=%s log=%s\n', cfg$lineage, method, job$pid, job$monitor_pid, worker_log_path))
  job
}

launch_or_skip_method_worker <- function(cfg, method, lineage_dir) {
  output_dir <- file.path(lineage_dir, method_output_subdir(method))
  result_path <- file.path(output_dir, 'method_worker_result.json')
  if (isTRUE(SKIP_COMPLETED_METHODS) && file.exists(result_path)) {
    result <- pa_json_read(result_path, default = list(status = 'missing_result'))
    result_status <- pa_null_coalesce(result$status, 'missing_result')
    if (identical(result_status, 'ok')) {
      append_status(
        cfg$lineage,
        method,
        'skipped_existing_result',
        output_dir,
        sprintf('result_status=ok; result=%s', result_path)
      )
      cat(sprintf('[SKIP] %s / %s existing ok result=%s\n', cfg$lineage, method, result_path))
      return(NULL)
    }
  }
  launch_method_worker(cfg, method, lineage_dir)
}

wait_for_jobs <- function(jobs, phase_label, poll_seconds = POLL_SECONDS) {
  if (length(jobs) == 0L) return(invisible(jobs))
  write_banner('Waiting for %s jobs (%d)', phase_label, length(jobs))

  remaining <- seq_along(jobs)
  while (length(remaining) > 0L) {
    still_running <- logical(length(remaining))
    for (i in seq_along(remaining)) {
      idx <- remaining[[i]]
      job <- jobs[[idx]]
      running <- is_pid_running(job$pid)
      still_running[[i]] <- running
      if (!running) {
        peak <- read_monitor_peak(job$monitor_tsv)
        result <- pa_json_read(job$result_path, default = list(status = 'missing_result'))
        result_status <- pa_null_coalesce(result$status, 'missing_result')
        append_status(
          job$lineage,
          job$method,
          if (identical(result_status, 'ok')) 'joined_completed' else 'joined_check_result',
          job$output_dir,
          sprintf('pid=%s; result_status=%s; peak_rss_gb=%s; min_available_gb=%s; log=%s',
            job$pid,
            result_status,
            ifelse(is.finite(peak$peak_rss_kb), sprintf('%.2f', peak$peak_rss_kb / 1024^2), 'NA'),
            ifelse(is.finite(peak$min_available_kb), sprintf('%.2f', peak$min_available_kb / 1024^2), 'NA'),
            job$log_path
          )
        )
        cat(sprintf('[JOIN] %s / %s pid=%s status=%s peak_rss_gb=%s\n',
          job$lineage,
          job$method,
          job$pid,
          result_status,
          ifelse(is.finite(peak$peak_rss_kb), sprintf('%.2f', peak$peak_rss_kb / 1024^2), 'NA')
        ))
      }
    }
    remaining <- remaining[still_running]
    if (length(remaining) > 0L) {
      labels <- vapply(jobs[remaining], function(job) sprintf('%s/%s(pid=%s)', job$lineage, job$method, job$pid), character(1))
      cat(sprintf('[WAIT] %s still running: %s\n', phase_label, paste(labels, collapse = ', ')))
      Sys.sleep(as.numeric(poll_seconds))
    }
  }
  invisible(jobs)
}

write_resource_summary <- function(lineage_dir) {
  monitor_files <- list.files(file.path(lineage_dir, 'resource_monitor'), pattern = '_memory_monitor\\.tsv$', full.names = TRUE)
  if (length(monitor_files) == 0L) return(NULL)
  rows <- lapply(monitor_files, function(path) {
    peak <- read_monitor_peak(path)
    method <- sub('_memory_monitor\\.tsv$', '', basename(path))
    data.frame(
      method = method,
      monitor_tsv = path,
      peak_rss_gb = ifelse(is.finite(peak$peak_rss_kb), peak$peak_rss_kb / 1024^2, NA_real_),
      min_available_gb = ifelse(is.finite(peak$min_available_kb), peak$min_available_kb / 1024^2, NA_real_),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  summary_df <- do.call(rbind, rows)
  out_path <- file.path(lineage_dir, 'lineage_resource_summary.tsv')
  pa_write_tsv(summary_df, out_path)
  out_path
}

validate_inputs <- function() {
  missing <- character()
  for (cfg in lineage_configs) {
    if (!file.exists(cfg$rds_path)) missing <- c(missing, cfg$rds_path)
    if (!file.exists(cfg$h5ad_path)) missing <- c(missing, cfg$h5ad_path)
  }
  for (path in c(WORKER_R, MONITOR_SH, PYCOGAPS_HELPER, HARMONY2_HELPER, CNMF_PYTHON, HARMONY2_PYTHON, PYCOGAPS_PYTHON)) {
    if (!file.exists(path)) missing <- c(missing, path)
  }
  if (length(missing) > 0L) {
    stop(sprintf('Missing required inputs/executables:\n%s', paste(unique(missing), collapse = '\n')), call. = FALSE)
  }
  TRUE
}

validate_inputs()

manifest <- list(
  run_root = RUN_ROOT,
  status_tsv = STATUS_TSV,
  run_stamp = RUN_STAMP,
  scheduler = 'lineage_sequential__methods_parallel__llm_parallel_second_phase',
  compute_methods = COMPUTE_METHODS,
  llm_methods = LLM_METHODS,
  monitor = list(default_interval_sec = MONITOR_INTERVAL_DEFAULT, cnmf_interval_sec = MONITOR_INTERVAL_CNMF),
  cnmf_python = CNMF_PYTHON,
  harmony2_python = HARMONY2_PYTHON,
  harmony2_helper = HARMONY2_HELPER,
  pycogaps_python = PYCOGAPS_PYTHON,
  lineages = lineage_configs
)
pa_write_json(manifest, MANIFEST_JSON)
append_status('all', 'scheduler', 'started', RUN_ROOT, sprintf('manifest=%s', MANIFEST_JSON))

for (cfg in lineage_configs) {
  lineage <- cfg$lineage
  lineage_dir <- pa_prepare_output_dir(file.path(RUN_ROOT, lineage))
  write_banner('Lineage sequential run: %s | methods parallel', lineage)
  cat(sprintf('[INFO] RDS  : %s\n', cfg$rds_path))
  cat(sprintf('[INFO] h5ad : %s\n', cfg$h5ad_path))
  append_status(lineage, 'lineage', 'started', lineage_dir, 'launching compute methods in parallel')

  compute_jobs <- Filter(Negate(is.null), lapply(COMPUTE_METHODS, function(method) launch_or_skip_method_worker(cfg, method, lineage_dir)))
  wait_for_jobs(compute_jobs, sprintf('%s compute', lineage), poll_seconds = POLL_SECONDS)

  append_status(lineage, 'lineage', 'compute_phase_completed', lineage_dir, 'launching downstream LLM workers in parallel')
  llm_jobs <- Filter(Negate(is.null), lapply(LLM_METHODS, function(method) launch_or_skip_method_worker(cfg, method, lineage_dir)))
  wait_for_jobs(llm_jobs, sprintf('%s LLM', lineage), poll_seconds = POLL_SECONDS)

  resource_summary <- write_resource_summary(lineage_dir)
  append_status(lineage, 'lineage', 'completed', lineage_dir, sprintf('resource_summary=%s', pa_null_coalesce(resource_summary, NA_character_)))
  gc(verbose = FALSE)
}

append_status('all', 'scheduler', 'completed', RUN_ROOT, sprintf('status=%s', STATUS_TSV))
write_banner('Parallel-method full run finished: %s', RUN_ROOT)
cat(sprintf('[INFO] Status TSV   : %s\n', STATUS_TSV))
cat(sprintf('[INFO] Manifest JSON: %s\n', MANIFEST_JSON))
