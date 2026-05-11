#!/usr/bin/env Rscript

source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop('Usage: program_python_tail_worker_20260507.R <tail_config.json>', call. = FALSE)
}

cfg_path <- normalizePath(args[[1]], winslash = '/', mustWork = TRUE)
cfg <- pa_json_read(cfg_path, default = NULL)
if (is.null(cfg) || !is.list(cfg)) {
  stop(sprintf('Failed to read tail config: %s', cfg_path), call. = FALSE)
}

STATUS_TSV <- pa_scalar_chr(cfg$status_tsv, 'status_tsv')
PYCOGAPS_HELPER <- pa_scalar_chr(cfg$pycogaps_helper, 'pycogaps_helper')

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
}

run_pycogaps_runner_local <- function(adata_h5ad_path,
                                      output_dir,
                                      run_name,
                                      python_cmd,
                                      celltype_col,
                                      helper_py = PYCOGAPS_HELPER,
                                      n_patterns = 6L,
                                      n_iterations = 50L,
                                      seed = 42L,
                                      n_threads = 1L,
                                      layer = 'counts',
                                      gene_exclusion_config = list()) {
  adata_h5ad_path <- pa_scalar_chr(adata_h5ad_path, 'adata_h5ad_path')
  helper_py <- pa_scalar_chr(helper_py, 'helper_py')
  output_dir <- pa_prepare_output_dir(output_dir)
  py_exec <- pa_resolve_python_executable(python_cmd)

  config <- list(
    helper_py = normalizePath(helper_py, winslash = '/', mustWork = FALSE),
    adata_h5ad_path = normalizePath(adata_h5ad_path, winslash = '/', mustWork = FALSE),
    output_dir = output_dir,
    run_name = pa_scalar_chr(run_name, 'run_name'),
    n_patterns = as.integer(n_patterns),
    n_iterations = as.integer(n_iterations),
    seed = as.integer(seed),
    n_threads = as.integer(n_threads),
    layer = pa_scalar_chr(layer, 'layer'),
    celltype_col = pa_null_coalesce(celltype_col, NULL),
    gene_exclusion_config = pa_null_coalesce(gene_exclusion_config, list())
  )

  config_path <- file.path(output_dir, 'pycogaps_runner_config.json')
  launcher_path <- file.path(output_dir, 'run_pycogaps_runner_launcher.py')

  launcher_lines <- c(
    '#!/usr/bin/env python3',
    'import importlib.util',
    'import json',
    'import pathlib',
    'import scanpy as sc',
    'import sys',
    '',
    'cfg_path = pathlib.Path(sys.argv[1])',
    'cfg = json.loads(cfg_path.read_text(encoding="utf-8"))',
    'helper_path = pathlib.Path(cfg["helper_py"])',
    'spec = importlib.util.spec_from_file_location("pa_pycogaps_helper", helper_path)',
    'module = importlib.util.module_from_spec(spec)',
    'assert spec.loader is not None',
    'spec.loader.exec_module(module)',
    'adata = sc.read_h5ad(cfg["adata_h5ad_path"])',
    'res = module.run_pycogaps_full(',
    '    adata=adata,',
    '    output_dir=pathlib.Path(cfg["output_dir"]),',
    '    run_name=cfg["run_name"],',
    '    n_patterns=cfg.get("n_patterns", 6),',
    '    n_iterations=cfg.get("n_iterations", 50),',
    '    seed=cfg.get("seed", 42),',
    '    n_threads=cfg.get("n_threads", 1),',
    '    layer=cfg.get("layer", "counts"),',
    '    celltype_col=cfg.get("celltype_col"),',
    '    gene_exclusion_config=cfg.get("gene_exclusion_config") or None,',
    ')',
    'out_path = pathlib.Path(cfg["output_dir"]) / "pycogaps_runner_result.json"',
    'out_path.write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str), encoding="utf-8")',
    'print(json.dumps({"success": bool(res.get("success")), "result_json": str(out_path)}, ensure_ascii=False))'
  )

  pa_write_json(config, config_path)
  pa_write_markdown(launcher_lines, launcher_path)
  Sys.chmod(launcher_path, mode = '0755')

  cmd_res <- pa_run_system_command(py_exec, args = c(launcher_path, config_path), fail_on_error = FALSE)
  runner_result <- pa_json_read(file.path(output_dir, 'pycogaps_runner_result.json'), default = list(success = FALSE, output = cmd_res$output))

  list(
    status = if (isTRUE(runner_result$success)) 'ok' else 'error',
    plan = list(
      python = py_exec,
      config_path = config_path,
      launcher_path = launcher_path,
      output_dir = output_dir
    ),
    command_result = cmd_res,
    runner_result = runner_result
  )
}

is_result_error <- function(res) {
  inherits(res, 'error') || (is.list(res) && identical(res$status, 'error'))
}

result_note <- function(res) {
  if (inherits(res, 'error')) return(conditionMessage(res))
  if (is.list(res) && identical(res$status, 'error')) return('runner returned status=error')
  'ok'
}

run_method_safe <- function(lineage, method, output_dir, expr) {
  append_status(lineage, method, 'started', output_dir, NA_character_)
  started <- Sys.time()
  res <- tryCatch(force(expr), error = function(e) e)
  elapsed_min <- round(as.numeric(difftime(Sys.time(), started, units = 'mins')), 2)

  if (is_result_error(res)) {
    append_status(lineage, method, 'error', output_dir, sprintf('%s | elapsed_min=%.2f', result_note(res), elapsed_min))
    cat(sprintf('[ERROR] %s / %s failed after %.2f min: %s\n', lineage, method, elapsed_min, result_note(res)))
    return(list(ok = FALSE, result = res, elapsed_min = elapsed_min))
  }

  append_status(lineage, method, 'completed', output_dir, sprintf('elapsed_min=%.2f', elapsed_min))
  cat(sprintf('[OK] %s / %s completed in %.2f min\n', lineage, method, elapsed_min))
  list(ok = TRUE, result = res, elapsed_min = elapsed_min)
}

lineage <- pa_scalar_chr(cfg$lineage, 'lineage')
lineage_dir <- pa_prepare_output_dir(cfg$lineage_dir)
tail_result_path <- file.path(lineage_dir, 'python_tail_result.json')

append_status(lineage, 'python_tail', 'worker_started', lineage_dir, sprintf('config=%s', cfg_path))

cnmf_dir <- pa_prepare_output_dir(cfg$cnmf_dir)
cnmf_run <- run_method_safe(
  lineage,
  'cnmf',
  cnmf_dir,
  {
    res <- pa_run_cnmf_runner(
      adata_h5ad_path = cfg$h5ad_path,
      output_dir = cnmf_dir,
      run_name = sprintf('%s_cnmf_full_%s', lineage, cfg$run_stamp),
      python_cmd = cfg$cnmf_python,
      dry_run = FALSE,
      celltype_col = cfg$celltype_col,
      batch_col = cfg$batch_col,
      use_batch_hvg = TRUE,
      cnmf_config = list(),
      viz_config = list()
    )
    saveRDS(res, file.path(cnmf_dir, 'cnmf_runner_wrapper_result.rds'))
    res
  }
)

gc(verbose = FALSE)

pycogaps_dir <- pa_prepare_output_dir(cfg$pycogaps_dir)
pycogaps_run <- run_method_safe(
  lineage,
  'pycogaps',
  pycogaps_dir,
  {
    res <- run_pycogaps_runner_local(
      adata_h5ad_path = cfg$h5ad_path,
      output_dir = pycogaps_dir,
      run_name = sprintf('%s_pycogaps_full_%s', lineage, cfg$run_stamp),
      python_cmd = cfg$pycogaps_python,
      celltype_col = cfg$celltype_col,
      n_patterns = 6L,
      n_iterations = 50L,
      seed = 42L,
      n_threads = 1L,
      layer = 'counts',
      gene_exclusion_config = list()
    )
    saveRDS(res, file.path(pycogaps_dir, 'pycogaps_runner_wrapper_result.rds'))
    res
  }
)

ok <- isTRUE(cnmf_run$ok) && isTRUE(pycogaps_run$ok)
summary <- list(
  lineage = lineage,
  status = if (ok) 'ok' else 'error',
  cnmf_ok = isTRUE(cnmf_run$ok),
  pycogaps_ok = isTRUE(pycogaps_run$ok),
  cnmf_elapsed_min = cnmf_run$elapsed_min,
  pycogaps_elapsed_min = pycogaps_run$elapsed_min,
  finished_at = format(Sys.time(), '%Y-%m-%d %H:%M:%S')
)
pa_write_json(summary, tail_result_path)
append_status(lineage, 'python_tail', if (ok) 'completed' else 'error', lineage_dir, sprintf('result=%s', tail_result_path))

if (!ok) quit(status = 1L, save = 'no')
quit(status = 0L, save = 'no')
