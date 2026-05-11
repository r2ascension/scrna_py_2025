#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

RUN_ROOT <- pa_prepare_output_dir('/home/h2048/output/program_full_serial_20260506')
STATUS_TSV <- file.path(RUN_ROOT, 'run_status.tsv')
MANIFEST_JSON <- file.path(RUN_ROOT, 'lineage_manifest.json')

CNMF_PYTHON <- '/home/h2048/miniconda3/envs/bbknn_env/bin/python'
PYCOGAPS_PYTHON <- '/home/h2048/miniconda3/envs/scarches_stable_pertpy/bin/python'
PYCOGAPS_HELPER <- '/home/h2048/script/py/pycogaps_helper_20260505_v1.py'
PYTHON_TAIL_WORKER <- '/home/h2048/script/R/program_python_tail_worker_20260507.R'
RUN_STAMP <- '20260506'

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

pa_write_json(lineage_configs, MANIFEST_JSON)

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

write_banner <- function(...) {
  cat('\n', paste0(rep('=', 88), collapse = ''), '\n', sep = '')
  cat(sprintf(...), '\n')
  cat(paste0(rep('=', 88), collapse = ''), '\n', sep = '')
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

run_method_safe <- function(lineage, method, output_dir, expr) {
  append_status(lineage, method, 'started', output_dir, NA_character_)
  started <- Sys.time()
  res <- tryCatch(
    force(expr),
    error = function(e) e
  )
  elapsed_min <- round(as.numeric(difftime(Sys.time(), started, units = 'mins')), 2)

  if (inherits(res, 'error')) {
    append_status(lineage, method, 'error', output_dir, sprintf('%s | elapsed_min=%.2f', conditionMessage(res), elapsed_min))
    cat(sprintf('[ERROR] %s / %s failed after %.2f min: %s\n', lineage, method, elapsed_min, conditionMessage(res)))
    return(NULL)
  }

  append_status(lineage, method, 'completed', output_dir, sprintf('elapsed_min=%.2f', elapsed_min))
  cat(sprintf('[OK] %s / %s completed in %.2f min\n', lineage, method, elapsed_min))
  res
}

is_pid_running <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (length(pid) == 0L || is.na(pid) || pid <= 0L) return(FALSE)
  status <- suppressWarnings(system2('kill', args = c('-0', as.character(pid)), stdout = FALSE, stderr = FALSE))
  identical(status, 0L)
}

launch_python_tail_async <- function(cfg, lineage_dir) {
  lineage <- cfg$lineage
  cnmf_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'cnmf_full'))
  pycogaps_dir <- pa_prepare_output_dir(file.path(lineage_dir, 'pycogaps_full'))
  tail_config_path <- file.path(lineage_dir, 'python_tail_config.json')
  tail_log_path <- file.path(lineage_dir, 'python_tail_async.log')
  tail_result_path <- file.path(lineage_dir, 'python_tail_result.json')

  tail_config <- list(
    lineage = lineage,
    lineage_dir = lineage_dir,
    h5ad_path = cfg$h5ad_path,
    celltype_col = cfg$celltype_col,
    batch_col = cfg$batch_col,
    cnmf_dir = cnmf_dir,
    pycogaps_dir = pycogaps_dir,
    cnmf_python = CNMF_PYTHON,
    pycogaps_python = PYCOGAPS_PYTHON,
    pycogaps_helper = PYCOGAPS_HELPER,
    status_tsv = STATUS_TSV,
    run_stamp = RUN_STAMP
  )
  pa_write_json(tail_config, tail_config_path)

  cmd <- sprintf(
    'env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 /usr/bin/Rscript %s %s > %s 2>&1 & echo $!',
    shQuote(PYTHON_TAIL_WORKER),
    shQuote(tail_config_path),
    shQuote(tail_log_path)
  )
  pid <- suppressWarnings(as.integer(system(cmd, intern = TRUE)[1]))
  append_status(
    lineage,
    'python_tail',
    'started_async',
    lineage_dir,
    sprintf('pid=%s; log=%s; result=%s', pid, tail_log_path, tail_result_path)
  )
  cat(sprintf('[ASYNC] %s / cNMF+PyCoGAPS tail started pid=%s log=%s\n', lineage, pid, tail_log_path))

  list(
    lineage = lineage,
    pid = pid,
    config_path = tail_config_path,
    log_path = tail_log_path,
    result_path = tail_result_path,
    lineage_dir = lineage_dir
  )
}

wait_for_python_tail_jobs <- function(jobs, poll_seconds = 60L) {
  if (length(jobs) == 0L) return(invisible(jobs))
  write_banner('Waiting for asynchronous cNMF/PyCoGAPS tails (%d jobs)', length(jobs))

  remaining <- seq_along(jobs)
  while (length(remaining) > 0L) {
    still_running <- logical(length(remaining))
    for (i in seq_along(remaining)) {
      idx <- remaining[[i]]
      job <- jobs[[idx]]
      still_running[[i]] <- is_pid_running(job$pid)
      if (!still_running[[i]]) {
        result <- pa_json_read(job$result_path, default = list(status = 'missing_result'))
        append_status(
          job$lineage,
          'python_tail',
          if (identical(result$status, 'ok')) 'joined_completed' else 'joined_check_result',
          job$lineage_dir,
          sprintf('pid=%s; result_status=%s; log=%s', job$pid, pa_null_coalesce(result$status, 'missing'), job$log_path)
        )
        cat(sprintf('[JOIN] %s python tail finished pid=%s status=%s\n', job$lineage, job$pid, pa_null_coalesce(result$status, 'missing')))
      }
    }
    remaining <- remaining[still_running]
    if (length(remaining) > 0L) {
      cat(sprintf('[WAIT] %d python tails still running: %s\n', length(remaining), paste(vapply(jobs[remaining], `[[`, character(1), 'lineage'), collapse = ', ')))
      Sys.sleep(as.numeric(poll_seconds))
    }
  }

  invisible(jobs)
}

python_tail_jobs <- list()

for (cfg in lineage_configs) {
  lineage <- cfg$lineage
  lineage_dir <- pa_prepare_output_dir(file.path(RUN_ROOT, lineage))
  write_banner('Serial full run: %s', lineage)

  cat(sprintf('[INFO] RDS  : %s\n', cfg$rds_path))
  cat(sprintf('[INFO] h5ad : %s\n', cfg$h5ad_path))

  seurat_obj <- readRDS(cfg$rds_path)

  hd_dir <- file.path(lineage_dir, 'hdwgcna_full')
  hd_res <- run_method_safe(
    lineage,
    'hdwgcna',
    hd_dir,
    {
      res <- pa_run_hdwgcna_runner(
        seurat_obj = seurat_obj,
        output_dir = hd_dir,
        celltypes = NULL,
        celltype_col = cfg$celltype_col,
        sample_col = cfg$sample_col,
        condition_col = cfg$condition_col,
        tissue_col = cfg$tissue_col,
        dry_run = FALSE
      )
      saveRDS(res, file.path(hd_dir, 'hdwgcna_runner_result.rds'))
      res
    }
  )
  rm(hd_res)
  gc(verbose = FALSE)

  covar_dir <- file.path(lineage_dir, 'covarnet_full')
  covar_res <- run_method_safe(
    lineage,
    'covarnet',
    covar_dir,
    {
      res <- pa_run_covarnet_runner(
        seurat_obj = seurat_obj,
        output_dir = covar_dir,
        celltypes = NULL,
        celltype_col = cfg$celltype_col,
        dry_run = FALSE
      )
      saveRDS(res, file.path(covar_dir, 'covarnet_runner_result.rds'))
      res
    }
  )
  rm(covar_res, seurat_obj)
  gc(verbose = FALSE)

  python_tail_jobs[[length(python_tail_jobs) + 1L]] <- launch_python_tail_async(cfg, lineage_dir)
  gc(verbose = FALSE)
}

wait_for_python_tail_jobs(python_tail_jobs)

write_banner('Serial full run finished: %s', RUN_ROOT)
cat(sprintf('[INFO] Status TSV   : %s\n', STATUS_TSV))
cat(sprintf('[INFO] Manifest JSON: %s\n', MANIFEST_JSON))
