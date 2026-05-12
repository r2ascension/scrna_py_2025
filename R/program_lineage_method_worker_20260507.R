#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

options(future.globals.maxSize = 16 * 1024^3)

PYCOGAPS_HELPER_DEFAULT <- '/home/h2048/script/py/pycogaps_helper_20260505_v1.py'
GENE_EXCLUSION_HELPER_R <- '/home/h2048/script/R/program_gene_exclusion_helper_20260505_v1.R'
DEEPSEEK_HELPER_R <- '/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R'

status_lock_eval <- function(status_tsv, expr) {
  lock_dir <- paste0(status_tsv, '.lock')
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

append_status_locked <- function(status_tsv, lineage, method, status, output_dir, note = NA_character_) {
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
  status_lock_eval(status_tsv, {
    utils::write.table(
      row,
      file = status_tsv,
      sep = '\t',
      row.names = FALSE,
      col.names = !file.exists(status_tsv),
      quote = FALSE,
      append = file.exists(status_tsv),
      na = ''
    )
  })
}

write_worker_result <- function(result, path) {
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file = sub('\\.json$', '.rds', path))
  if (requireNamespace('jsonlite', quietly = TRUE)) {
    jsonlite::write_json(result, path = path, pretty = TRUE, auto_unbox = TRUE, null = 'null')
  }
  invisible(path)
}

run_pycogaps_runner_local <- function(adata_h5ad_path,
                                      output_dir,
                                      run_name,
                                      python_cmd,
                                      celltype_col,
                                      helper_py = PYCOGAPS_HELPER_DEFAULT,
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

load_env_file_safely <- function(path) {
  if (!file.exists(path)) return(FALSE)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl('^\\s*(#|$)', lines)]
  for (ln in lines) {
    if (!grepl('=', ln, fixed = TRUE)) next
    key <- trimws(sub('=.*$', '', ln))
    val <- sub('^[^=]*=', '', ln)
    val <- trimws(val)
    val <- sub('^export\\s+', '', val)
    val <- sub('^(["\'])(.*)\\1$', '\\2', val)
    if (nzchar(key) && !nzchar(Sys.getenv(key, unset = ''))) {
      do.call(Sys.setenv, stats::setNames(list(val), key))
    }
  }
  TRUE
}

has_live_deepseek_key <- function() {
  for (env_path in c('/home/h2048/.env', '/home/h2048/script/.env')) load_env_file_safely(env_path)
  key <- Sys.getenv('DEEPSEEK_API_KEY', unset = '')
  if (!nzchar(key) || nchar(key) < 20L) return(FALSE)
  !grepl('your|placeholder|api[_-]?key[_-]?here|dummy|example', key, ignore.case = TRUE)
}

read_table_preview <- function(path, n = 8L, max_chars = 1800L) {
  out <- tryCatch({
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c('md')) {
      txt <- paste(utils::head(readLines(path, warn = FALSE), n), collapse = '\n')
      if (nchar(txt) > max_chars) paste0(substr(txt, 1L, max_chars), '\n...[truncated]') else txt
    } else {
    df <- pa_read_table_auto(path)
    if (!is.data.frame(df) || nrow(df) == 0L) return('[empty table]')
    df <- utils::head(df, n)
    txt <- paste(capture.output(utils::write.table(df, sep = '\t', row.names = FALSE, quote = FALSE, na = '')), collapse = '\n')
    if (nchar(txt) > max_chars) paste0(substr(txt, 1L, max_chars), '\n...[truncated]') else txt
    }
  }, error = function(e) sprintf('[preview unavailable: %s]', conditionMessage(e)))
  out
}

llm_plot_exts <- function() c('png', 'pdf', 'svg', 'jpg', 'jpeg', 'tif', 'tiff')

llm_table_exts <- function() c('csv', 'tsv', 'txt', 'md', 'json')

llm_source_role <- function(base, rel_path, ext) {
  base_l <- tolower(base)
  rel_l <- tolower(rel_path)
  if (grepl('gene_exclusion.*(llm|prompt|interpret)', base_l)) return('exclude')
  if (ext %in% llm_plot_exts()) return('visualization_evidence')
  if (grepl('visual|plot|figure|network|heatmap|clustergram|dotplot|umap|stability|trait|hub_overlap|direction_summary|node_degree|network_summary', base_l) ||
      grepl('visual|plot|figure|network|heatmap|clustergram|dotplot|umap|stability|trait|hub_overlap', rel_l)) return('visualization_evidence')
  if (grepl('gene_exclusion|audit', base_l)) return('qc_context')
  if (grepl('manifest|summary', base_l)) return('run_context')
  if (grepl('LLM|prompt|interpret', base, ignore.case = TRUE)) return('llm_sidecar')
  if (grepl('gep|program|module|edge|node|usage|score|factor|pattern', base_l)) return('program_evidence')
  'other'
}

find_candidate_data_for_plot <- function(plot_rel, table_rel) {
  if (length(table_rel) == 0L) return(character())
  plot_l <- tolower(plot_rel)
  table_l <- tolower(table_rel)
  candidates <- character()
  add_matches <- function(pattern) {
    hits <- table_rel[grepl(pattern, table_l, perl = TRUE)]
    candidates <<- unique(c(candidates, hits))
  }

  if (grepl('gene_exclusion', plot_l)) add_matches('gene_exclusion')
  if (grepl('covarnet', plot_l)) {
    add_matches('covarnet_.*_(edges|nodes)\\.(csv|tsv|txt)$')
    add_matches('covarnet_.*_(network|edge_direction|node_degree).*\\.(csv|tsv|json)$')
  }
  if (grepl('hdwgcna', plot_l)) {
    add_matches('hdwgcna_module_membership\\.csv$')
    add_matches('hdwgcna_celltype_status\\.json$')
  }
  if (grepl('cnmf|clustergram|gep|stability', plot_l)) {
    add_matches('gep_(top_genes|gene_scores)_k[0-9]+\\.(csv|tsv)$')
    add_matches('(k_stability|k_selection|run_summary).*\\.json$')
    add_matches('(usage|spectra|score).*\\.(txt|tsv|csv)$')
  }
  if (grepl('pycogaps|pattern', plot_l)) {
    add_matches('pycogaps_(gene_patterns|cell_patterns)\\.tsv$')
    add_matches('pycogaps_(top_genes_by_pattern|summary)\\.json$')
  }

  same_dir <- dirname(plot_rel)
  if (!identical(same_dir, '.')) {
    candidates <- unique(c(candidates, table_rel[dirname(table_rel) == same_dir]))
  }
  head(unique(candidates), 12L)
}

write_program_visualization_manifest <- function(source_dir) {
  if (!dir.exists(source_dir)) return(NULL)
  all_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  all_files <- all_files[file.exists(all_files) & !dir.exists(all_files)]
  if (length(all_files) == 0L) return(NULL)
  source_norm <- normalizePath(source_dir, winslash = '/', mustWork = FALSE)
  file_norm <- normalizePath(all_files, winslash = '/', mustWork = FALSE)
  rel <- ifelse(
    startsWith(file_norm, paste0(source_norm, '/')),
    substr(file_norm, nchar(source_norm) + 2L, nchar(file_norm)),
    basename(file_norm)
  )
  ext <- tolower(tools::file_ext(all_files))
  plot_idx <- ext %in% llm_plot_exts()
  if (!any(plot_idx)) return(NULL)

  table_rel <- rel[ext %in% llm_table_exts()]
  rows <- lapply(which(plot_idx), function(i) {
    candidates <- find_candidate_data_for_plot(rel[[i]], table_rel)
    data.frame(
      figure_path = normalizePath(all_files[[i]], winslash = '/', mustWork = FALSE),
      figure_rel_path = rel[[i]],
      figure_ext = ext[[i]],
      figure_role = if (grepl('gene_exclusion', rel[[i]], ignore.case = TRUE)) 'qc_visualization' else 'method_visualization',
      paired_data_rel_paths = paste(candidates, collapse = ';'),
      n_paired_data_candidates = length(candidates),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  manifest <- do.call(rbind, rows)
  manifest <- manifest[order(manifest$figure_role, manifest$figure_rel_path), , drop = FALSE]
  tsv_path <- file.path(source_dir, 'program_visualization_manifest.tsv')
  utils::write.table(manifest, file = tsv_path, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
  pa_write_json(as.data.frame(manifest, stringsAsFactors = FALSE), file.path(source_dir, 'program_visualization_manifest.json'))
  tsv_path
}

discover_prior_llm_context <- function(output_dir, source_method) {
  prior_path <- file.path(output_dir, sprintf('%s_LLM_interpretation.md', source_method))
  if (!file.exists(prior_path)) return(data.frame())
  data.frame(
    role = 'prior_llm_context',
    ext = 'md',
    path = normalizePath(prior_path, winslash = '/', mustWork = FALSE),
    rel_path = basename(prior_path),
    size_bytes = file.info(prior_path)$size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

discover_llm_source_files <- function(source_dir) {
  if (!dir.exists(source_dir)) return(data.frame())
  write_program_visualization_manifest(source_dir)
  all_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  all_files <- all_files[file.exists(all_files) & !dir.exists(all_files)]
  if (length(all_files) == 0L) return(data.frame())

  source_norm <- normalizePath(source_dir, winslash = '/', mustWork = FALSE)
  file_norm <- normalizePath(all_files, winslash = '/', mustWork = FALSE)
  rel <- ifelse(
    startsWith(file_norm, paste0(source_norm, '/')),
    substr(file_norm, nchar(source_norm) + 2L, nchar(file_norm)),
    basename(file_norm)
  )
  ext <- tolower(tools::file_ext(all_files))
  base <- basename(all_files)
  role <- mapply(llm_source_role, base = base, rel_path = rel, ext = ext, USE.NAMES = FALSE)
  keep <- ext %in% c(llm_table_exts(), llm_plot_exts()) & !(role %in% c('other', 'exclude'))
  df <- data.frame(
    role = role[keep],
    ext = ext[keep],
    path = normalizePath(all_files[keep], winslash = '/', mustWork = FALSE),
    rel_path = rel[keep],
    size_bytes = file.info(all_files[keep])$size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(df) == 0L) return(df)
  role_order <- c('prior_llm_context', 'visualization_evidence', 'program_evidence', 'run_context', 'qc_context', 'llm_sidecar')
  df <- df[order(match(df$role, role_order), df$ext %in% llm_plot_exts(), df$size_bytes, df$rel_path), , drop = FALSE]
  utils::head(df, 80L)
}

run_llm_summary_worker <- function(cfg, output_dir) {
  source_method <- sub('^llm_', '', cfg$method)
  source_dir <- file.path(cfg$lineage_dir, switch(
    source_method,
    hdwgcna = 'hdwgcna_full',
    covarnet = 'covarnet_full',
    cnmf = 'cnmf_full',
    pycogaps = 'pycogaps_full',
    source_method
  ))
  output_dir <- pa_prepare_output_dir(output_dir)
  index_df <- discover_llm_source_files(source_dir)
  prior_idx <- discover_prior_llm_context(output_dir, source_method)
  if (nrow(prior_idx) > 0L) {
    index_df <- rbind(prior_idx, index_df)
  }
  index_path <- file.path(output_dir, sprintf('%s_llm_input_index.tsv', source_method))
  if (nrow(index_df) > 0L) {
    pa_write_tsv(index_df, index_path)
  } else {
    pa_write_tsv(data.frame(role = character(), ext = character(), path = character(), rel_path = character(), size_bytes = numeric()), index_path)
  }

  preview_files <- if (nrow(index_df) > 0L) index_df[index_df$ext %in% llm_table_exts(), , drop = FALSE] else index_df
  preview_role_order <- c('visualization_evidence', 'program_evidence', 'prior_llm_context', 'run_context', 'qc_context', 'llm_sidecar')
  if (nrow(preview_files) > 0L) {
    preview_files <- preview_files[order(match(preview_files$role, preview_role_order), preview_files$size_bytes, preview_files$rel_path), , drop = FALSE]
  }
  preview_files <- utils::head(preview_files, 16L)
  preview_blocks <- character()
  if (nrow(preview_files) > 0L) {
    for (i in seq_len(nrow(preview_files))) {
      preview_blocks <- c(
        preview_blocks,
        sprintf('### %s | %s', preview_files$role[[i]], preview_files$rel_path[[i]]),
        '```',
        read_table_preview(preview_files$path[[i]], n = 8L, max_chars = 1800L),
        '```',
        ''
      )
    }
  }

  figure_files <- if (nrow(index_df) > 0L) index_df[index_df$ext %in% llm_plot_exts(), , drop = FALSE] else index_df
  figure_lines <- if (nrow(figure_files) > 0L) {
    paste0('- ', utils::head(figure_files$rel_path, 30L))
  } else {
    '[No figure files indexed.]'
  }

  prompt <- c(
    sprintf('# %s / %s program-analysis LLM review', cfg$lineage, source_method),
    '',
    '你是单细胞转录组和功能程序分析专家。请基于下面索引、图表清单和表格预览，给出简体中文的初步解释。',
    '注意：当前在线 LLM 主要读取文本预览，不能直接渲染 PDF/PNG；若索引中有 `program_visualization_manifest.tsv`，请优先用它理解“图表 ↔ 支撑数据表”的对应关系，并指出最值得人工查看的图表。',
    '',
    '请输出：',
    '1. 主要可解释的生物学程序/模块；',
    '2. 可能的技术或组成混杂（尤其 IG/MT/RPS/RPL/lncRNA 等已排除或仍需警惕的信号）；',
    '3. 最值得人工复核的程序、文件或图表；',
    '4. 后续验证/统计建议。',
    '',
    sprintf('- lineage: `%s`', cfg$lineage),
    sprintf('- method: `%s`', source_method),
    sprintf('- source_dir: `%s`', source_dir),
    sprintf('- input_index: `%s`', index_path),
    '',
    '## Indexed figure files',
    figure_lines,
    '',
    '## Evidence previews',
    if (length(preview_blocks) > 0L) preview_blocks else '[No prompt-ready source files found yet.]'
  )
  prompt_path <- file.path(output_dir, sprintf('%s_LLM_prompt.md', source_method))
  interpretation_path <- file.path(output_dir, sprintf('%s_LLM_interpretation.md', source_method))
  raw_path <- file.path(output_dir, sprintf('%s_LLM_raw_response.txt', source_method))
  status_path <- file.path(output_dir, sprintf('%s_LLM_status.json', source_method))
  pa_write_markdown(prompt, prompt_path)

  status <- list(
    status = 'queued_no_live_key',
    lineage = cfg$lineage,
    method = source_method,
    prompt_md = prompt_path,
    interpretation_md = interpretation_path,
    raw_response_txt = raw_path,
    input_index_tsv = index_path,
    n_indexed_files = nrow(index_df),
    model = pa_null_coalesce(cfg$llm_model, 'deepseek-reasoner'),
    error = NULL
  )

  llm_enable_live <- isTRUE(pa_null_coalesce(cfg$llm_enable_live, TRUE))
  if (!llm_enable_live || !has_live_deepseek_key()) {
    status$status <- if (llm_enable_live) 'queued_no_live_key' else 'queued_live_disabled'
    pa_write_markdown(c(
      sprintf('# %s / %s LLM interpretation', cfg$lineage, source_method),
      '',
      if (llm_enable_live) '当前没有检测到可用的 `DEEPSEEK_API_KEY`，因此已生成 LLM prompt 和索引，但未调用在线模型。' else '当前 worker 配置禁用了在线 LLM 调用，因此已生成 LLM prompt 和索引，但未调用在线模型。',
      '',
      sprintf('- Prompt: `%s`', prompt_path),
      sprintf('- Input index: `%s`', index_path),
      sprintf('- Indexed files: %d', nrow(index_df)),
      '',
      '可在设置真实 API key 后重跑对应 `llm_*` worker。'
    ), interpretation_path)
    pa_write_json(status, status_path)
    status$status_json <- status_path
    return(list(status = 'ok', llm = status, index = index_df))
  }

  response <- tryCatch({
    if (!exists('tc_deepseek_chat_request', mode = 'function')) source(DEEPSEEK_HELPER_R)
    tc_deepseek_chat_request(
      prompt = paste(prompt, collapse = '\n'),
      model = pa_null_coalesce(cfg$llm_model, 'deepseek-reasoner'),
      api_key = Sys.getenv('DEEPSEEK_API_KEY', unset = ''),
      timeout_sec = as.numeric(pa_null_coalesce(cfg$llm_timeout_sec, 240))
    )
  }, error = function(e) structure(conditionMessage(e), class = 'pa_llm_worker_error'))

  if (inherits(response, 'pa_llm_worker_error')) {
    status$status <- 'error'
    status$error <- as.character(response[[1]])
    pa_write_markdown(c(
      sprintf('# %s / %s LLM interpretation failed', cfg$lineage, source_method),
      '',
      sprintf('Error: `%s`', status$error),
      '',
      sprintf('Prompt retained at `%s`.', prompt_path)
    ), interpretation_path)
  } else {
    status$status <- 'ok'
    pa_write_markdown(c(sprintf('# %s / %s LLM interpretation', cfg$lineage, source_method), '', response), interpretation_path)
    writeLines(as.character(response), raw_path, useBytes = TRUE)
  }
  pa_write_json(status, status_path)
  status$status_json <- status_path
  list(status = 'ok', llm = status, index = index_df)
}

run_compute_method <- function(cfg, output_dir) {
  method <- cfg$method
  lineage <- cfg$lineage
  if (identical(method, 'hdwgcna')) {
    seurat_obj <- readRDS(cfg$rds_path)
    on.exit(rm(seurat_obj), add = TRUE)
    hdwgcna_runner_args <- pa_null_coalesce(cfg$hdwgcna_runner_args, list())
    if (is.null(hdwgcna_runner_args$resume_celltypes)) {
      hdwgcna_runner_args$resume_celltypes <- TRUE
    }
    if (is.null(hdwgcna_runner_args$resume_skip_statuses)) {
      skip_env <- Sys.getenv('HDWGCNA_RESUME_SKIP_STATUSES', unset = 'ok,no_modules,timeout')
      hdwgcna_runner_args$resume_skip_statuses <- trimws(strsplit(skip_env, ',', fixed = TRUE)[[1]])
      hdwgcna_runner_args$resume_skip_statuses <- hdwgcna_runner_args$resume_skip_statuses[nzchar(hdwgcna_runner_args$resume_skip_statuses)]
      if (length(hdwgcna_runner_args$resume_skip_statuses) == 0L) {
        hdwgcna_runner_args$resume_skip_statuses <- c('ok', 'no_modules', 'timeout')
      }
    }
    if (is.null(hdwgcna_runner_args$celltype_timeout_sec)) {
      timeout_env <- Sys.getenv('HDWGCNA_CELLTYPE_TIMEOUT_SEC', unset = '21600')
      hdwgcna_runner_args$celltype_timeout_sec <- suppressWarnings(as.integer(timeout_env))
      if (length(hdwgcna_runner_args$celltype_timeout_sec) == 0L || is.na(hdwgcna_runner_args$celltype_timeout_sec)) {
        hdwgcna_runner_args$celltype_timeout_sec <- 21600L
      }
    }
    hdwgcna_celltypes <- pa_null_coalesce(cfg$hdwgcna_celltypes, NULL)
    if (is.null(hdwgcna_celltypes)) {
      celltypes_env <- Sys.getenv('HDWGCNA_CELLTYPES', unset = '')
      if (nzchar(celltypes_env)) {
        hdwgcna_celltypes <- trimws(strsplit(celltypes_env, ',', fixed = TRUE)[[1]])
        hdwgcna_celltypes <- hdwgcna_celltypes[nzchar(hdwgcna_celltypes)]
      }
    }
    res <- pa_run_hdwgcna_runner(
      seurat_obj = seurat_obj,
      output_dir = output_dir,
      celltypes = hdwgcna_celltypes,
      celltype_col = cfg$celltype_col,
      sample_col = cfg$sample_col,
      condition_col = cfg$condition_col,
      tissue_col = cfg$tissue_col,
      runner_args = hdwgcna_runner_args,
      dry_run = FALSE
    )
    saveRDS(res, file.path(output_dir, 'hdwgcna_runner_result.rds'))
    return(res)
  }

  if (identical(method, 'covarnet')) {
    seurat_obj <- readRDS(cfg$rds_path)
    on.exit(rm(seurat_obj), add = TRUE)
    res <- pa_run_covarnet_runner(
      seurat_obj = seurat_obj,
      output_dir = output_dir,
      celltypes = NULL,
      celltype_col = cfg$celltype_col,
      plot_networks = TRUE,
      plot_overlap = TRUE,
      enhanced_visualizations = TRUE,
      dry_run = FALSE
    )
    saveRDS(res, file.path(output_dir, 'covarnet_runner_result.rds'))
    return(res)
  }

  if (identical(method, 'cnmf')) {
    res <- pa_run_cnmf_runner(
      adata_h5ad_path = cfg$h5ad_path,
      output_dir = output_dir,
      run_name = paste0(lineage, '_cnmf_full_', cfg$run_stamp),
      python_cmd = cfg$cnmf_python,
      unit_id = paste0(lineage, '__cnmf_full'),
      dry_run = FALSE,
      k_range = NULL,
      celltype_col = cfg$celltype_col,
      batch_col = cfg$batch_col,
      use_batch_hvg = TRUE,
      cnmf_config = list(n_workers = 1L),
      viz_config = list()
    )
    saveRDS(res, file.path(output_dir, 'cnmf_runner_result.rds'))
    if (identical(res$status, 'error')) stop('cNMF runner returned status=error', call. = FALSE)
    return(res)
  }

  if (identical(method, 'pycogaps')) {
    res <- run_pycogaps_runner_local(
      adata_h5ad_path = cfg$h5ad_path,
      output_dir = output_dir,
      run_name = paste0(lineage, '_pycogaps_full_', cfg$run_stamp),
      python_cmd = cfg$pycogaps_python,
      celltype_col = cfg$celltype_col,
      helper_py = pa_null_coalesce(cfg$pycogaps_helper, PYCOGAPS_HELPER_DEFAULT),
      n_patterns = as.integer(pa_null_coalesce(cfg$pycogaps_n_patterns, 6L)),
      n_iterations = as.integer(pa_null_coalesce(cfg$pycogaps_n_iterations, 50L)),
      seed = as.integer(pa_null_coalesce(cfg$seed, 42L)),
      n_threads = as.integer(pa_null_coalesce(cfg$pycogaps_n_threads, 1L)),
      layer = 'counts',
      gene_exclusion_config = list(lineage_context = lineage)
    )
    saveRDS(res, file.path(output_dir, 'pycogaps_runner_result.rds'))
    if (identical(res$status, 'error')) stop('PyCoGAPS runner returned status=error', call. = FALSE)
    return(res)
  }

  if (grepl('^llm_', method)) {
    return(run_llm_summary_worker(cfg, output_dir))
  }

  stop(sprintf('Unknown method: %s', method), call. = FALSE)
}

program_lineage_method_worker_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop('Usage: program_lineage_method_worker_20260507.R <worker_config.json>', call. = FALSE)
  config_path <- args[[1]]
  cfg <- pa_json_read(config_path, default = NULL)
  if (is.null(cfg)) stop(sprintf('Could not read worker config JSON: %s', config_path), call. = FALSE)

  required <- c('lineage', 'method', 'lineage_dir', 'output_dir', 'status_tsv', 'run_stamp')
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0L) stop(sprintf('Worker config missing fields: %s', paste(missing, collapse = ', ')), call. = FALSE)

  lineage <- cfg$lineage
  method <- cfg$method
  output_dir <- pa_prepare_output_dir(cfg$output_dir)
  result_path <- file.path(output_dir, 'method_worker_result.json')

  append_status_locked(cfg$status_tsv, lineage, method, 'started', output_dir, sprintf('pid=%s; config=%s', Sys.getpid(), config_path))
  started <- Sys.time()
  res <- tryCatch(
    run_compute_method(cfg, output_dir),
    error = function(e) e
  )
  elapsed_min <- round(as.numeric(difftime(Sys.time(), started, units = 'mins')), 2)

  if (inherits(res, 'error')) {
    result <- list(
      status = 'error',
      lineage = lineage,
      method = method,
      output_dir = output_dir,
      config_path = config_path,
      elapsed_min = elapsed_min,
      error = conditionMessage(res),
      ended_at = format(Sys.time(), '%Y-%m-%d %H:%M:%S')
    )
    write_worker_result(result, result_path)
    append_status_locked(cfg$status_tsv, lineage, method, 'error', output_dir, sprintf('%s | elapsed_min=%.2f', conditionMessage(res), elapsed_min))
    return(invisible(1L))
  }

  method_status <- if (!is.null(res$status)) as.character(res$status)[[1]] else 'ok'
  result <- list(
    status = if (identical(method_status, 'error')) 'error' else 'ok',
    method_status = method_status,
    lineage = lineage,
    method = method,
    output_dir = output_dir,
    config_path = config_path,
    elapsed_min = elapsed_min,
    ended_at = format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  )
  write_worker_result(result, result_path)
  append_status_locked(cfg$status_tsv, lineage, method, 'completed', output_dir, sprintf('method_status=%s; elapsed_min=%.2f; result=%s', method_status, elapsed_min, result_path))
  invisible(if (identical(result$status, 'error')) 1L else 0L)
}

if (identical(sys.nframe(), 0L)) {
  status_code <- program_lineage_method_worker_main()
  quit(save = 'no', status = status_code)
}
