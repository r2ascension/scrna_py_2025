#!/usr/bin/env Rscript

source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

RUN_ROOT <- Sys.getenv('PROGRAM_RUN_ROOT', unset = '/home/h2048/output/program_full_parallel_methods_20260507')
STATUS_TSV <- file.path(RUN_ROOT, 'run_status.tsv')
WORKER_R <- Sys.getenv('PROGRAM_METHOD_WORKER_R', unset = '/home/h2048/script/R/program_lineage_method_worker_20260507.R')
RUN_STAMP <- Sys.getenv('PROGRAM_LLM_RUN_STAMP', unset = '20260512_llm_existing_methods')

split_env <- function(name, default) {
  x <- Sys.getenv(name, unset = '')
  if (!nzchar(x)) return(default)
  trimws(strsplit(x, ',', fixed = TRUE)[[1]])
}

SOURCE_METHODS <- split_env('PROGRAM_LLM_SOURCE_METHODS', c('covarnet', 'cnmf', 'pycogaps', 'hdwgcna'))
LINEAGES <- split_env('PROGRAM_LLM_LINEAGES', character())
MAX_PARALLEL <- suppressWarnings(as.integer(Sys.getenv('PROGRAM_LLM_MAX_PARALLEL', unset = '8')))
if (length(MAX_PARALLEL) == 0L || is.na(MAX_PARALLEL) || MAX_PARALLEL <= 0L) MAX_PARALLEL <- 8L
FORCE <- Sys.getenv('PROGRAM_LLM_FORCE', unset = '0') %in% c('1', 'true', 'TRUE', 'yes', 'YES')
LLM_ENABLE_LIVE <- !(Sys.getenv('PROGRAM_LLM_ENABLE_LIVE', unset = '1') %in% c('0', 'false', 'FALSE', 'no', 'NO'))
LLM_MODEL <- Sys.getenv('PROGRAM_LLM_MODEL', unset = 'deepseek-reasoner')
LLM_TIMEOUT_SEC <- suppressWarnings(as.integer(Sys.getenv('PROGRAM_LLM_TIMEOUT_SEC', unset = '240')))
if (length(LLM_TIMEOUT_SEC) == 0L || is.na(LLM_TIMEOUT_SEC) || LLM_TIMEOUT_SEC <= 0L) LLM_TIMEOUT_SEC <- 240L

method_output_subdir <- function(method) {
  switch(method,
    hdwgcna = 'hdwgcna_full',
    covarnet = 'covarnet_full',
    cnmf = 'cnmf_full',
    pycogaps = 'pycogaps_full',
    method
  )
}

append_status <- function(lineage, method, status, output_dir, note = NA_character_) {
  row <- data.frame(
    timestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
    lineage = lineage,
    method = method,
    status = status,
    output_dir = output_dir,
    note = as.character(note),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(row, file = STATUS_TSV, sep = '\t', row.names = FALSE,
                     col.names = !file.exists(STATUS_TSV), quote = FALSE,
                     append = file.exists(STATUS_TSV), na = '')
}

is_pid_running <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (length(pid) == 0L || is.na(pid) || pid <= 0L) return(FALSE)
  identical(suppressWarnings(system2('kill', args = c('-0', as.character(pid)), stdout = FALSE, stderr = FALSE)), 0L)
}

existing_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  x <- pa_json_read(path, default = NULL)
  is.list(x) && identical(as.character(x$status)[1], 'ok')
}

if (!dir.exists(RUN_ROOT)) stop(sprintf('RUN_ROOT does not exist: %s', RUN_ROOT), call. = FALSE)
if (length(LINEAGES) == 0L) {
  dirs <- list.dirs(RUN_ROOT, recursive = FALSE, full.names = FALSE)
  LINEAGES <- sort(dirs[nzchar(dirs)])
}

launch_one <- function(lineage, source_method) {
  lineage_dir <- file.path(RUN_ROOT, lineage)
  source_dir <- file.path(lineage_dir, method_output_subdir(source_method))
  if (!dir.exists(source_dir)) {
    append_status(lineage, paste0('llm_', source_method), 'skipped_missing_source', source_dir, 'source method dir missing')
    return(NULL)
  }
  source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  source_files <- source_files[file.exists(source_files) & !dir.exists(source_files)]
  has_evidence <- any(grepl('\\.(csv|tsv|txt|json|md|png|pdf|svg)$', source_files, ignore.case = TRUE))
  if (!has_evidence) {
    append_status(lineage, paste0('llm_', source_method), 'skipped_no_evidence', source_dir, 'no prompt-ready evidence files')
    return(NULL)
  }

  out_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'llm_parallel', source_method))
  result_path <- file.path(out_dir, 'method_worker_result.json')
  if (!FORCE && existing_ok(result_path)) {
    append_status(lineage, paste0('llm_', source_method), 'skipped_existing_result', out_dir, sprintf('result=%s', result_path))
    return(NULL)
  }

  config_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'worker_configs'))
  logs_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'worker_logs'))
  config_path <- file.path(config_dir, sprintf('llm_%s_worker_config.json', source_method))
  log_path <- file.path(logs_dir, sprintf('llm_%s_worker.log', source_method))
  cfg <- list(
    lineage = lineage,
    method = paste0('llm_', source_method),
    lineage_dir = lineage_dir,
    output_dir = out_dir,
    status_tsv = STATUS_TSV,
    run_root = RUN_ROOT,
    run_stamp = RUN_STAMP,
    llm_model = LLM_MODEL,
    llm_timeout_sec = LLM_TIMEOUT_SEC,
    llm_enable_live = LLM_ENABLE_LIVE
  )
  pa_write_json(cfg, config_path)
  cmd <- sprintf(
    'env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 /usr/bin/Rscript %s %s > %s 2>&1 & echo $!',
    shQuote(WORKER_R), shQuote(config_path), shQuote(log_path)
  )
  pid <- suppressWarnings(as.integer(system(cmd, intern = TRUE)[1]))
  append_status(lineage, paste0('llm_', source_method), 'launched_async', out_dir, sprintf('pid=%s; log=%s; source_dir=%s', pid, log_path, source_dir))
  list(lineage = lineage, method = paste0('llm_', source_method), source_method = source_method,
       pid = pid, output_dir = out_dir, result_path = result_path, log_path = log_path)
}

jobs <- list()
for (lineage in LINEAGES) {
  for (source_method in SOURCE_METHODS) {
    while (length(jobs) >= MAX_PARALLEL) {
      alive <- vapply(jobs, function(job) is_pid_running(job$pid), logical(1))
      jobs <- jobs[alive]
      if (length(jobs) >= MAX_PARALLEL) Sys.sleep(5)
    }
    job <- launch_one(lineage, source_method)
    if (!is.null(job)) jobs[[length(jobs) + 1L]] <- job
  }
}

cat(sprintf('[INFO] Launched %d LLM workers; waiting for completion.\n', length(jobs)))
repeat {
  if (length(jobs) == 0L) break
  alive <- vapply(jobs, function(job) is_pid_running(job$pid), logical(1))
  finished <- jobs[!alive]
  if (length(finished) > 0L) {
    for (job in finished) {
      status <- if (existing_ok(job$result_path)) 'completed' else 'check_result'
      append_status(job$lineage, job$method, status, job$output_dir, sprintf('pid=%s; result=%s; log=%s', job$pid, job$result_path, job$log_path))
      cat(sprintf('[%s] %s/%s %s\n', format(Sys.time(), '%F %T'), job$lineage, job$method, status))
    }
  }
  jobs <- jobs[alive]
  if (length(jobs) > 0L) {
    cat(sprintf('[WAIT] %d LLM workers still running: %s\n', length(jobs), paste(vapply(jobs, function(job) sprintf('%s/%s(pid=%s)', job$lineage, job$method, job$pid), character(1)), collapse = ', ')))
    Sys.sleep(15)
  }
}
cat('[OK] LLM existing-method launcher finished.\n')
