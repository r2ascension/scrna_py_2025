#!/usr/bin/env Rscript

WORKER_R <- Sys.getenv('PROGRAM_METHOD_WORKER_R', unset = '/home/h2048/script/R/program_lineage_method_worker_20260507.R')
if (file.exists(WORKER_R)) source(WORKER_R)
if (!exists('pa_prepare_output_dir', mode = 'function')) source('/home/h2048/script/R/program_architecture_bundle_20260428_v1.R')

RUN_ROOT_DEFAULT <- '/home/h2048/output/program_full_parallel_methods_20260507'

split_env <- function(name, default = character()) {
  x <- Sys.getenv(name, unset = '')
  if (!nzchar(x)) return(default)
  out <- trimws(strsplit(x, ',', fixed = TRUE)[[1]])
  out[nzchar(out)]
}

safe_unit_id <- function(x) {
  y <- gsub('[^A-Za-z0-9_]+', '_', as.character(x))
  y <- gsub('_+', '_', y)
  y <- gsub('^_|_$', '', y)
  ifelse(nzchar(y), y, 'unit')
}

method_output_subdir_unit_llm <- function(method) {
  switch(method,
    cnmf = 'cnmf_full',
    covarnet = 'covarnet_full',
    pycogaps = 'pycogaps_full',
    hdwgcna = 'hdwgcna_full',
    method
  )
}

read_json_unit_llm <- function(path, default = list()) {
  if (!file.exists(path)) return(default)
  if (exists('pa_json_read', mode = 'function')) return(pa_json_read(path, default = default))
  if (!requireNamespace('jsonlite', quietly = TRUE)) return(default)
  tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) default)
}

write_json_unit_llm <- function(x, path) {
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (exists('pa_write_json', mode = 'function')) return(pa_write_json(x, path))
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = 'null')
}

write_lines_unit_llm <- function(lines, path) {
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(lines), con = path, useBytes = TRUE)
  invisible(path)
}

read_table_unit_llm <- function(path) {
  if (exists('pa_read_table_auto', mode = 'function')) return(pa_read_table_auto(path))
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, 'csv')) return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  utils::read.delim(path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
}

preview_any_file <- function(path, n = 50L, max_chars = 4500L) {
  out <- tryCatch({
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c('png', 'pdf', 'svg', 'jpg', 'jpeg', 'tif', 'tiff')) {
      sprintf('[figure file: %s]', normalizePath(path, winslash = '/', mustWork = FALSE))
    } else if (ext %in% c('md', 'txt', 'json')) {
      txt <- paste(utils::head(readLines(path, warn = FALSE), n), collapse = '\n')
      if (nchar(txt) > max_chars) paste0(substr(txt, 1L, max_chars), '\n...[truncated]') else txt
    } else {
      df <- read_table_unit_llm(path)
      if (!is.data.frame(df) || nrow(df) == 0L) return('[empty table]')
      df <- utils::head(df, n)
      txt <- paste(capture.output(utils::write.table(df, sep = '\t', row.names = FALSE, quote = FALSE, na = '')), collapse = '\n')
      if (nchar(txt) > max_chars) paste0(substr(txt, 1L, max_chars), '\n...[truncated]') else txt
    }
  }, error = function(e) sprintf('[preview unavailable: %s]', conditionMessage(e)))
  out
}

program_low_information_gene_pattern <- function() {
  paste(
    c(
      '^MT-',
      '^RPS[0-9A-Z]*$', '^RPL[0-9A-Z]*$', '^MRPS[0-9A-Z]*$', '^MRPL[0-9A-Z]*$',
      '^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$',
      '^(AC|AL|AP|BX|Z)[0-9]+[.]',
      '^RP[0-9]+-',
      '^CTD-', '^CTB-', '^CTC-',
      '-OT[0-9]+$',
      '^LOC[0-9]+'
    ),
    collapse = '|'
  )
}

filter_program_low_information_genes <- function(tbl, gene_col = 'gene') {
  if (!is.data.frame(tbl) || nrow(tbl) == 0L || !gene_col %in% colnames(tbl)) return(tbl)
  genes <- as.character(tbl[[gene_col]])
  keep <- !grepl(program_low_information_gene_pattern(), genes)
  keep[is.na(keep)] <- FALSE
  out <- tbl[keep, , drop = FALSE]
  if ('rank' %in% colnames(out)) {
    out$source_rank <- out$rank
    out$rank <- seq_len(nrow(out))
  }
  rownames(out) <- NULL
  out
}

select_hdwgcna_tissue_evidence <- function(unit_dir) {
  files <- list.files(unit_dir, pattern = '[.](json|png|pdf|svg|csv)$', full.names = TRUE)
  if (length(files) == 0L) return(character())
  base <- basename(files)
  is_condition_or_disease <- grepl('condition|disease|crswnp|healthy|case|control', base, ignore.case = TRUE)
  is_generic_condition_dotplot <- identical(base, 'hdwgcna_hub_dotplot.pdf')
  keep <- !is_condition_or_disease & !is_generic_condition_dotplot & (
    grepl('tissue', base, ignore.case = TRUE) |
      grepl('hub_rank|top_hubs|module_sizes|soft_power|celltype_status|module_umap', base, ignore.case = TRUE)
  )
  sort(files[keep])
}

selected_cnmf_score_files <- function(source_dir, all_cnmf_k = TRUE) {
  gep_dir <- file.path(source_dir, 'gep_gene_tables')
  score_files <- list.files(gep_dir, pattern = '^gep_gene_scores_k[0-9]+\\.tsv$', full.names = TRUE)
  if (length(score_files) == 0L) return(character())
  if (isTRUE(all_cnmf_k)) return(sort(score_files))
  rec <- read_json_unit_llm(file.path(source_dir, 'k_selection_recommendation.json'), default = list())
  run_summary <- read_json_unit_llm(file.path(source_dir, 'run_summary.json'), default = list())
  selected_k <- suppressWarnings(as.integer(unlist(c(run_summary$recommendation$recommended_k, rec$recommended_k), use.names = FALSE)[1]))
  if (length(selected_k) == 0L || is.na(selected_k)) {
    parsed <- suppressWarnings(as.integer(sub('^.*_k([0-9]+)\\.tsv$', '\\1', score_files)))
    selected_k <- parsed[which.max(parsed)]
  }
  selected_path <- file.path(gep_dir, sprintf('gep_gene_scores_k%d.tsv', selected_k))
  if (file.exists(selected_path)) selected_path else sort(score_files)[1]
}

new_task_row <- function(lineage, source_method, unit_type, unit_id, unit_label, source_dir, output_dir, evidence_paths, extra = list()) {
  payload_path <- file.path(output_dir, paste0(unit_id, '_payload.json'))
  prompt_path <- file.path(output_dir, paste0(unit_id, '_LLM_prompt.md'))
  interpretation_path <- file.path(output_dir, paste0(unit_id, '_LLM_interpretation.md'))
  status_path <- file.path(output_dir, paste0(unit_id, '_LLM_status.json'))
  data.frame(
    lineage = lineage,
    source_method = source_method,
    unit_type = unit_type,
    unit_id = unit_id,
    unit_label = unit_label,
    source_dir = normalizePath(source_dir, winslash = '/', mustWork = FALSE),
    output_dir = normalizePath(output_dir, winslash = '/', mustWork = FALSE),
    evidence_paths = paste(normalizePath(evidence_paths[file.exists(evidence_paths)], winslash = '/', mustWork = FALSE), collapse = ';'),
    payload_json = normalizePath(payload_path, winslash = '/', mustWork = FALSE),
    prompt_md = normalizePath(prompt_path, winslash = '/', mustWork = FALSE),
    interpretation_md = normalizePath(interpretation_path, winslash = '/', mustWork = FALSE),
    status_json = normalizePath(status_path, winslash = '/', mustWork = FALSE),
    extra_json = if (length(extra)) jsonlite::toJSON(extra, auto_unbox = TRUE) else '{}',
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_cnmf_unit_rows <- function(lineage, lineage_dir, source_dir, all_cnmf_k = TRUE) {
  score_files <- selected_cnmf_score_files(source_dir, all_cnmf_k = all_cnmf_k)
  if (length(score_files) == 0L) return(data.frame())
  rows <- list()
  for (score_path in score_files) {
    score_tbl <- tryCatch(read_table_unit_llm(score_path), error = function(e) data.frame())
    if (!is.data.frame(score_tbl) || nrow(score_tbl) == 0L || !all(c('gep', 'gene') %in% colnames(score_tbl))) next
    k <- suppressWarnings(as.integer(sub('^.*_k([0-9]+)\\.tsv$', '\\1', score_path)))
    for (gep in unique(as.character(score_tbl$gep))) {
      unit_id <- safe_unit_id(sprintf('cnmf_k%s_%s', k, gep))
      out_dir <- file.path(lineage_dir, 'llm_parallel', 'cnmf', 'units', unit_id)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      gep_gene_path <- file.path(out_dir, paste0(unit_id, '_gene_scores.csv'))
      gep_tbl <- score_tbl[as.character(score_tbl$gep) == gep, , drop = FALSE]
      gep_tbl <- filter_program_low_information_genes(gep_tbl, gene_col = 'gene')
      tryCatch(utils::write.csv(gep_tbl, gep_gene_path, row.names = FALSE), error = function(e) NULL)
      fig_paths <- list.files(file.path(source_dir, 'visualizations'), pattern = sprintf('k%s\\.', k), full.names = TRUE)
      rows[[length(rows) + 1L]] <- new_task_row(
        lineage = lineage,
        source_method = 'cnmf',
        unit_type = 'gep',
        unit_id = unit_id,
        unit_label = sprintf('k%s / %s', k, gep),
        source_dir = source_dir,
        output_dir = out_dir,
        evidence_paths = c(gep_gene_path, file.path(source_dir, 'k_selection_recommendation.json'), fig_paths),
        extra = list(k = k, gep = gep, gep_gene_table = normalizePath(gep_gene_path, winslash = '/', mustWork = FALSE), score_table = normalizePath(score_path, winslash = '/', mustWork = FALSE))
      )
    }
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

build_pycogaps_unit_rows <- function(lineage, lineage_dir, source_dir) {
  top_path <- file.path(source_dir, 'pycogaps_top_genes_by_pattern.json')
  if (!file.exists(top_path)) return(data.frame())
  top_genes <- read_json_unit_llm(top_path, default = list())
  patterns <- names(top_genes)
  if (length(patterns) == 0L) return(data.frame())
  rows <- lapply(patterns, function(pattern) {
    unit_id <- safe_unit_id(paste0('pycogaps_', pattern))
    out_dir <- file.path(lineage_dir, 'llm_parallel', 'pycogaps', 'units', unit_id)
    new_task_row(
      lineage = lineage,
      source_method = 'pycogaps',
      unit_type = 'pattern',
      unit_id = unit_id,
      unit_label = pattern,
      source_dir = source_dir,
      output_dir = out_dir,
      evidence_paths = c(top_path, file.path(source_dir, 'pycogaps_gene_patterns.tsv'), file.path(source_dir, 'pycogaps_cell_patterns.tsv'), list.files(source_dir, pattern = '^pycogaps_.*[.](png|pdf|svg)$', full.names = TRUE)),
      extra = list(pattern = pattern)
    )
  })
  do.call(rbind, rows)
}

build_covarnet_unit_rows <- function(lineage, lineage_dir, source_dir) {
  node_files <- list.files(source_dir, pattern = '^covarnet_.*_nodes[.]csv$', full.names = TRUE)
  if (length(node_files) == 0L) return(data.frame())
  rows <- lapply(node_files, function(node_path) {
    celltype <- sub('^covarnet_', '', sub('_nodes[.]csv$', '', basename(node_path)))
    unit_id <- safe_unit_id(paste0('covarnet_', celltype))
    out_dir <- file.path(lineage_dir, 'llm_parallel', 'covarnet', 'units', unit_id)
    edge_path <- file.path(source_dir, sprintf('covarnet_%s_edges.csv', celltype))
    network_pdf <- file.path(source_dir, sprintf('covarnet_%s_network.pdf', celltype))
    new_task_row(
      lineage = lineage,
      source_method = 'covarnet',
      unit_type = 'celltype_network',
      unit_id = unit_id,
      unit_label = celltype,
      source_dir = source_dir,
      output_dir = out_dir,
      evidence_paths = c(node_path, edge_path, network_pdf, file.path(source_dir, 'covarnet_hub_overlap_heatmap.pdf')),
      extra = list(celltype = celltype)
    )
  })
  do.call(rbind, rows)
}

build_hdwgcna_unit_rows <- function(lineage, lineage_dir, source_dir) {
  unit_dirs <- list.dirs(source_dir, recursive = FALSE, full.names = TRUE)
  unit_dirs <- unit_dirs[grepl('/hdwgcna_', unit_dirs)]
  if (length(unit_dirs) == 0L) return(data.frame())
  rows <- list()
  for (unit_dir in unit_dirs) {
    celltype <- sub('^hdwgcna_', '', basename(unit_dir))
    membership_path <- file.path(unit_dir, 'hdwgcna_module_membership.csv')
    if (!file.exists(membership_path)) next
    module_tbl <- tryCatch(read_table_unit_llm(membership_path), error = function(e) data.frame())
    if (!is.data.frame(module_tbl) || nrow(module_tbl) == 0L || !'module' %in% colnames(module_tbl)) next
    modules <- sort(setdiff(unique(as.character(module_tbl$module)), 'grey'))
    if (length(modules) == 0L) next
    shared_evidence <- select_hdwgcna_tissue_evidence(unit_dir)
    for (module in modules) {
      unit_id <- safe_unit_id(paste('hdwgcna', celltype, module, sep = '_'))
      out_dir <- file.path(lineage_dir, 'llm_parallel', 'hdwgcna', 'units', unit_id)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      module_gene_path <- file.path(out_dir, paste0(unit_id, '_module_genes.csv'))
      module_gene_tbl <- module_tbl[as.character(module_tbl$module) == module, , drop = FALSE]
      tryCatch(utils::write.csv(module_gene_tbl, module_gene_path, row.names = FALSE), error = function(e) NULL)
      rows[[length(rows) + 1L]] <- new_task_row(
        lineage = lineage,
        source_method = 'hdwgcna',
        unit_type = 'hdwgcna_module',
        unit_id = unit_id,
        unit_label = paste(celltype, module, sep = ' / '),
        source_dir = source_dir,
        output_dir = out_dir,
        evidence_paths = c(module_gene_path, membership_path, shared_evidence),
        extra = list(celltype = celltype, module = module, unit_dir = normalizePath(unit_dir, winslash = '/', mustWork = FALSE))
      )
    }
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

build_unit_llm_task_index <- function(run_root = RUN_ROOT_DEFAULT,
                                      lineages = character(),
                                      source_methods = c('cnmf', 'pycogaps', 'covarnet'),
                                      all_cnmf_k = TRUE) {
  run_root <- normalizePath(run_root, winslash = '/', mustWork = FALSE)
  if (length(lineages) == 0L) {
    dirs <- list.dirs(run_root, recursive = FALSE, full.names = FALSE)
    lineages <- sort(dirs[nzchar(dirs)])
  }
  rows <- list()
  for (lineage in lineages) {
    lineage_dir <- file.path(run_root, lineage)
    if (!dir.exists(lineage_dir)) next
    for (method in source_methods) {
      source_dir <- file.path(lineage_dir, method_output_subdir_unit_llm(method))
      if (!dir.exists(source_dir)) next
      part <- switch(method,
        cnmf = build_cnmf_unit_rows(lineage, lineage_dir, source_dir, all_cnmf_k = all_cnmf_k),
        pycogaps = build_pycogaps_unit_rows(lineage, lineage_dir, source_dir),
        covarnet = build_covarnet_unit_rows(lineage, lineage_dir, source_dir),
        hdwgcna = build_hdwgcna_unit_rows(lineage, lineage_dir, source_dir),
        data.frame()
      )
      if (is.data.frame(part) && nrow(part) > 0L) rows[[length(rows) + 1L]] <- part
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      lineage = character(), source_method = character(), unit_type = character(), unit_id = character(),
      unit_label = character(), source_dir = character(), output_dir = character(), evidence_paths = character(),
      payload_json = character(), prompt_md = character(), interpretation_md = character(), status_json = character(),
      extra_json = character(), stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

parse_task_extra_unit_llm <- function(task) {
  raw <- as.character(task$extra_json)[1]
  if (!nzchar(raw) || identical(raw, '{}')) return(list())
  if (!requireNamespace('jsonlite', quietly = TRUE)) return(list())
  out <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = TRUE), error = function(e) list())
  if (is.null(out)) list() else out
}

build_task_metadata_lines_unit_llm <- function(task, extra = NULL) {
  if (is.null(extra)) extra <- parse_task_extra_unit_llm(task)
  lines <- character()
  if (identical(as.character(task$source_method)[1], 'cnmf')) {
    hierarchy <- character()
    if (length(extra$celltype_l2) > 0L && nzchar(as.character(extra$celltype_l2)[1])) {
      hierarchy <- c(hierarchy, sprintf('L2=%s', as.character(extra$celltype_l2)[1]))
    }
    if (length(extra$celltype_l3) > 0L && nzchar(as.character(extra$celltype_l3)[1])) {
      hierarchy <- c(hierarchy, sprintf('L3=%s', as.character(extra$celltype_l3)[1]))
    }
    if (length(hierarchy) > 0L) {
      lines <- c(lines, sprintf('- celltype_hierarchy: `%s`', paste(hierarchy, collapse = ' | ')))
    }
    if (length(extra$celltype_l2) > 0L && nzchar(as.character(extra$celltype_l2)[1])) {
      lines <- c(lines, sprintf('- celltype_l2: `%s`', as.character(extra$celltype_l2)[1]))
    }
    if (length(extra$celltype_l3) > 0L && nzchar(as.character(extra$celltype_l3)[1])) {
      lines <- c(lines, sprintf('- celltype_l3: `%s`', as.character(extra$celltype_l3)[1]))
    }
    if (length(extra$k) > 0L && !is.na(extra$k[1])) {
      lines <- c(lines, sprintf('- k: `%s`', as.character(extra$k[1])))
    }
    if (length(extra$gep) > 0L && nzchar(as.character(extra$gep)[1])) {
      lines <- c(lines, sprintf('- gep: `%s`', as.character(extra$gep)[1]))
    }
  }
  lines
}

build_unit_prompt <- function(task) {
  evidence <- unlist(strsplit(task$evidence_paths, ';', fixed = TRUE), use.names = FALSE)
  evidence <- evidence[nzchar(evidence) & file.exists(evidence)]
  extra <- parse_task_extra_unit_llm(task)
  metadata_lines <- build_task_metadata_lines_unit_llm(task, extra = extra)
  blocks <- character()
  for (path in evidence) {
    blocks <- c(blocks, sprintf('### %s', basename(path)), '```', preview_any_file(path), '```', '')
  }
  if (length(blocks) == 0L) blocks <- '[No unit evidence files found.]'
  method_notes <- character()
  if (identical(as.character(task$source_method), 'cnmf')) {
    method_notes <- c(
      '',
      '## Method-specific guardrails',
      '- 当前 unit 是一个明确的 `k × GEP` 单元；只解释这个 k 下这个 GEP，不要把同一 k 的其它 GEP 或其它 k 的 GEP 混入结论。',
      '- 核心证据以当前 unit 的 `_gene_scores.csv` 为准；该表已过滤 AC/AL/AP/RP/CTD/CTB/CTC/LOC、MT、核糖体等低可解释性或技术性基因。',
      '- 如果剩余 top genes 仍缺少清晰生物学一致性，应降低置信度，而不是强行命名。'
    )
  } else if (identical(as.character(task$source_method), 'hdwgcna')) {
    method_notes <- c(
      '',
      '## Method-specific guardrails',
      '- 所有输入细胞均按健康样本/健康组织背景解释；当前 hdWGCNA unit 应用于正常组织间比较，而不是疾病组、CRSwNP、case/control 或 Healthy-vs-disease 对比。',
      '- 只允许围绕模块基因、hub/kME、模块大小以及 tissue-level eigengene evidence 解释组织生态位、稳态功能或组织适应性。',
      '- 如果 unit label 或基因名出现 inflammatory/immune 等字样，也不要自动推断为疾病炎症；除非基因集合本身支持，只能表述为健康组织中的免疫/应激基线或细胞状态。',
      '- 输出正文不要复述排除性说明；如果需要说明局限，只写“仅限健康组织稳态解释，不能做组间方向性结论”。',
      '- 输出正文不得出现这些字面字符串或等价措辞：CRSwNP、Healthy、disease、condition、case/control、疾病组、病例组、对照组、疾病升高、疾病降低、病例对照差异。'
    )
  }
  c(
    sprintf('# %s / %s / %s LLM unit interpretation', task$lineage, task$source_method, task$unit_label),
    '',
    '你是单细胞转录组功能程序分析专家。请只解释当前这个 unit，不要泛泛总结整个方法。',
    '',
    sprintf('- lineage: `%s`', task$lineage),
    sprintf('- method: `%s`', task$source_method),
    sprintf('- unit_type: `%s`', task$unit_type),
    sprintf('- unit_id: `%s`', task$unit_id),
    sprintf('- unit_label: `%s`', task$unit_label),
    metadata_lines,
    '',
    '请输出简体中文 Markdown，包含：',
    '1. 一个不超过 8 个字的功能/状态命名；',
    '2. 关键支持基因或网络 hub；',
    '3. 可能的技术混杂或不应过度解释的信号；',
    '4. 置信度（high/medium/low）和一句理由；',
    '5. 下一步人工复核建议。',
    method_notes,
    '',
    '## Unit evidence',
    blocks
  )
}

materialize_unit_task <- function(task) {
  dir.create(task$output_dir, recursive = TRUE, showWarnings = FALSE)
  evidence <- unlist(strsplit(task$evidence_paths, ';', fixed = TRUE), use.names = FALSE)
  payload <- as.list(task)
  payload$evidence_paths <- evidence[nzchar(evidence)]
  payload$extra <- parse_task_extra_unit_llm(task)
  payload$created_at <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  write_json_unit_llm(payload, task$payload_json)
  prompt <- build_unit_prompt(task)
  write_lines_unit_llm(prompt, task$prompt_md)
  invisible(prompt)
}

unit_status_ok <- function(path) {
  x <- read_json_unit_llm(path, default = list())
  is.list(x) && identical(as.character(x$status)[1], 'ok')
}

call_unit_llm <- function(task, enable_live = TRUE, force = FALSE, model = 'deepseek-reasoner', timeout_sec = 240L) {
  if (!isTRUE(force) && file.exists(task$status_json) && unit_status_ok(task$status_json)) {
    return(list(status = 'skipped_existing_ok', unit_id = task$unit_id, status_json = task$status_json))
  }
  prompt <- materialize_unit_task(task)
  status <- list(
    status = 'queued_no_live_key',
    lineage = task$lineage,
    method = task$source_method,
    unit_type = task$unit_type,
    unit_id = task$unit_id,
    unit_label = task$unit_label,
    prompt_md = task$prompt_md,
    interpretation_md = task$interpretation_md,
    payload_json = task$payload_json,
    status_json = task$status_json,
    model = model,
    error = NULL,
    ended_at = format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  )
  if (!isTRUE(enable_live) || !has_live_deepseek_key()) {
    status$status <- if (isTRUE(enable_live)) 'queued_no_live_key' else 'queued_live_disabled'
    write_lines_unit_llm(c(
      sprintf('# %s / %s / %s unit LLM interpretation', task$lineage, task$source_method, task$unit_label),
      '',
      if (isTRUE(enable_live)) '未检测到可用 DEEPSEEK_API_KEY，已生成 per-unit prompt/payload，未调用在线模型。' else '在线 LLM 调用被禁用，已生成 per-unit prompt/payload。',
      '',
      sprintf('- Prompt: `%s`', task$prompt_md),
      sprintf('- Payload: `%s`', task$payload_json)
    ), task$interpretation_md)
    write_json_unit_llm(status, task$status_json)
    return(status)
  }
  response <- tryCatch({
    if (!exists('pa_call_deepseek_chat', mode = 'function')) source(WORKER_R)
    pa_call_deepseek_chat(
      prompt = paste(prompt, collapse = '\n'),
      model = model,
      api_key = Sys.getenv('DEEPSEEK_API_KEY', unset = ''),
      timeout_sec = as.numeric(timeout_sec)
    )
  }, error = function(e) structure(conditionMessage(e), class = 'program_unit_llm_error'))
  if (inherits(response, 'program_unit_llm_error')) {
    status$status <- 'error'
    status$error <- as.character(response[[1]])
    write_lines_unit_llm(c(sprintf('# %s failed', task$unit_id), '', sprintf('Error: `%s`', status$error), '', sprintf('Prompt retained at `%s`.', task$prompt_md)), task$interpretation_md)
  } else {
    status$status <- 'ok'
    write_lines_unit_llm(c(sprintf('# %s / %s / %s unit LLM interpretation', task$lineage, task$source_method, task$unit_label), '', as.character(response)), task$interpretation_md)
  }
  status$ended_at <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  write_json_unit_llm(status, task$status_json)
  status
}

run_unit_llm_tasks <- function(index_df, max_parallel = 8L, enable_live = TRUE, force = FALSE, model = 'deepseek-reasoner', timeout_sec = 240L) {
  if (!is.data.frame(index_df) || nrow(index_df) == 0L) return(list())
  max_parallel <- max(1L, min(as.integer(max_parallel), nrow(index_df)))
  tasks <- split(index_df, seq_len(nrow(index_df)))
  worker <- function(df) call_unit_llm(df[1, , drop = FALSE], enable_live = enable_live, force = force, model = model, timeout_sec = timeout_sec)
  if (.Platform$OS.type == 'unix' && max_parallel > 1L) {
    parallel::mclapply(tasks, worker, mc.cores = max_parallel)
  } else {
    lapply(tasks, worker)
  }
}

write_tsv_unit_llm <- function(x, path) {
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (exists('pa_write_tsv', mode = 'function')) return(pa_write_tsv(x, path))
  utils::write.table(x, file = path, sep = '\t', row.names = FALSE, quote = FALSE)
}

sanitize_text_unit_llm <- function(x) {
  if (is.null(x)) return(NA_character_)
  value <- as.character(x)
  out <- suppressWarnings(iconv(value, from = '', to = 'UTF-8', sub = 'byte'))
  replace <- is.na(out) & !is.na(value)
  if (any(replace)) out[replace] <- enc2utf8(value[replace])
  out
}

scalar_extra_unit_llm <- function(extra, name, default = NA_character_) {
  if (is.null(extra) || is.null(extra[[name]]) || length(extra[[name]]) == 0L) return(default)
  value <- sanitize_text_unit_llm(extra[[name]][1])
  if (is.null(value) || length(value) == 0L || is.na(value) || !nzchar(trimws(value))) default else value
}

numeric_suffix_unit_llm <- function(x) {
  x <- sanitize_text_unit_llm(x)
  out <- suppressWarnings(as.integer(sub('^.*?([0-9]+)$', '\\1', x)))
  out[is.na(out)] <- .Machine$integer.max
  out
}

read_unit_llm_markdown_body <- function(path) {
  if (!file.exists(path)) return('*[missing interpretation markdown]*')
  lines <- tryCatch(sanitize_text_unit_llm(readLines(path, warn = FALSE)), error = function(e) character())
  if (length(lines) == 0L) return('*[empty interpretation markdown]*')
  first_nonempty <- which(nzchar(trimws(lines)))[1]
  if (!is.na(first_nonempty) && grepl('^#\\s+', trimws(lines[first_nonempty]))) {
    lines <- lines[-seq_len(first_nonempty)]
    while (length(lines) > 0L && !nzchar(trimws(lines[1]))) lines <- lines[-1]
  }
  if (length(lines) == 0L) '*[empty interpretation markdown]*' else lines
}

task_status_record_unit_llm <- function(task) {
  status_path <- as.character(task$status_json)[1]
  status <- read_json_unit_llm(status_path, default = list())
  status_value <- sanitize_text_unit_llm(status$status %||% NA_character_)[1]
  if (is.na(status_value) || !nzchar(status_value)) {
    interp_path <- as.character(task$interpretation_md)[1]
    status_value <- if (file.exists(interp_path)) 'present_no_status' else 'missing'
  }
  list(
    status = status_value,
    error = sanitize_text_unit_llm(status$error %||% NA_character_)[1],
    started_at = sanitize_text_unit_llm(status$started_at %||% NA_character_)[1],
    ended_at = sanitize_text_unit_llm(status$ended_at %||% NA_character_)[1]
  )
}

unit_llm_index_records <- function(index_df) {
  if (!is.data.frame(index_df) || nrow(index_df) == 0L) return(data.frame())
  rows <- lapply(seq_len(nrow(index_df)), function(i) {
    task <- index_df[i, , drop = FALSE]
    extra <- parse_task_extra_unit_llm(task)
    status <- task_status_record_unit_llm(task)
    celltype_l1_value <- scalar_extra_unit_llm(extra, 'celltype_l1')
    celltype_l2_value <- scalar_extra_unit_llm(extra, 'celltype_l2')
    celltype_l3_value <- scalar_extra_unit_llm(extra, 'celltype_l3')
    celltype_value <- scalar_extra_unit_llm(extra, 'celltype')
    safe_celltype_value <- scalar_extra_unit_llm(
      extra,
      'safe_celltype',
      if (!is.na(celltype_l3_value) && nzchar(celltype_l3_value)) {
        safe_unit_id(celltype_l3_value)
      } else if (!is.na(celltype_value) && nzchar(celltype_value)) {
        safe_unit_id(celltype_value)
      } else if (!is.na(celltype_l2_value) && nzchar(celltype_l2_value)) {
        safe_unit_id(celltype_l2_value)
      } else {
        safe_unit_id(as.character(task$unit_label)[1])
      }
    )
    data.frame(
      lineage = sanitize_text_unit_llm(task$lineage)[1],
      source_method = sanitize_text_unit_llm(task$source_method)[1],
      unit_type = sanitize_text_unit_llm(task$unit_type)[1],
      unit_id = sanitize_text_unit_llm(task$unit_id)[1],
      unit_label = sanitize_text_unit_llm(task$unit_label)[1],
      source_dir = sanitize_text_unit_llm(task$source_dir)[1],
      unit_output_dir = sanitize_text_unit_llm(task$output_dir)[1],
      payload_json = sanitize_text_unit_llm(task$payload_json)[1],
      prompt_md = sanitize_text_unit_llm(task$prompt_md)[1],
      interpretation_md = sanitize_text_unit_llm(task$interpretation_md)[1],
      status_json = sanitize_text_unit_llm(task$status_json)[1],
      celltype = celltype_value,
      celltype_l1 = celltype_l1_value,
      celltype_l2 = celltype_l2_value,
      celltype_l3 = celltype_l3_value,
      celltype_hierarchy = scalar_extra_unit_llm(extra, 'celltype_hierarchy'),
      safe_celltype = safe_celltype_value,
      module = scalar_extra_unit_llm(extra, 'module'),
      pattern = scalar_extra_unit_llm(extra, 'pattern'),
      unit_dir = scalar_extra_unit_llm(extra, 'unit_dir'),
      k = suppressWarnings(as.integer(extra$k %||% NA_integer_)),
      gep = scalar_extra_unit_llm(extra, 'gep'),
      status = as.character(status$status %||% NA_character_)[1],
      error = as.character(status$error %||% NA_character_)[1],
      started_at = as.character(status$started_at %||% NA_character_)[1],
      ended_at = as.character(status$ended_at %||% NA_character_)[1],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$k_sort <- out$k
  out$k_sort[is.na(out$k_sort)] <- .Machine$integer.max
  out$gep_sort <- numeric_suffix_unit_llm(out$gep)
  out
}

markdown_table_lines_unit_llm <- function(df, cols) {
  if (!is.data.frame(df) || nrow(df) == 0L) return(c('| unit_id | status |', '|---|---|'))
  labels <- gsub('_', '\\_', sanitize_text_unit_llm(cols), fixed = TRUE)
  header <- paste0('| ', paste(labels, collapse = ' | '), ' |')
  divider <- paste0('| ', paste(rep('---', length(cols)), collapse = ' | '), ' |')
  body <- apply(df[, cols, drop = FALSE], 1, function(row) {
    values <- vapply(row, function(x) {
      value <- sanitize_text_unit_llm(x)
      if (is.na(value) || !nzchar(value)) '' else gsub('\\|', '\\\\|', value)
    }, character(1))
    paste0('| ', paste(values, collapse = ' | '), ' |')
  })
  c(header, divider, body)
}

collect_existing_cnmf_celltype_master_rows <- function(run_root,
                                                       output_subdir,
                                                       level_tag) {
  json_paths <- Sys.glob(file.path(
    run_root,
    '*',
    output_subdir,
    'llm_parallel',
    'cnmf',
    'celltypes',
    '*',
    sprintf('cnmf_%s_*_LLM_combined_manifest.json', level_tag)
  ))
  if (length(json_paths) == 0L) {
    return(data.frame(stringsAsFactors = FALSE, check.names = FALSE))
  }
  rows <- lapply(sort(unique(json_paths)), function(path) {
    manifest <- read_json_unit_llm(path, default = list())
    lineage <- sanitize_text_unit_llm(manifest$lineage %||% NA_character_)[1]
    safe_celltype <- sanitize_text_unit_llm(manifest$safe_celltype %||% basename(dirname(path)))[1]
    celltype_field <- sanitize_text_unit_llm(manifest$celltype_field %||% sprintf('celltype_%s', level_tag))[1]
    celltype_value <- sanitize_text_unit_llm(
      manifest$celltype_value %||% manifest[[celltype_field]] %||% safe_celltype %||% NA_character_
    )[1]
    if (is.na(lineage) || !nzchar(lineage) || is.na(safe_celltype) || !nzchar(safe_celltype)) {
      return(NULL)
    }
    data.frame(
      lineage = lineage,
      level_tag = level_tag,
      celltype_field = if (!is.na(celltype_field) && nzchar(celltype_field)) celltype_field else sprintf('celltype_%s', level_tag),
      celltype_value = celltype_value,
      safe_celltype = safe_celltype,
      total_units = suppressWarnings(as.integer(manifest$total_units %||% NA_integer_))[1],
      ok_units = suppressWarnings(as.integer(manifest$ok_units %||% NA_integer_))[1],
      error_units = suppressWarnings(as.integer(manifest$error_units %||% NA_integer_))[1],
      combined_md = sanitize_text_unit_llm(manifest$combined_md %||% NA_character_)[1],
      combined_manifest_tsv = sanitize_text_unit_llm(manifest$combined_manifest_tsv %||% sub('[.]json$', '.tsv', path))[1],
      combined_manifest_json = normalizePath(path, winslash = '/', mustWork = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(stringsAsFactors = FALSE, check.names = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- out[order(out$lineage, out$safe_celltype, out$celltype_value), , drop = FALSE]
  rownames(out) <- NULL
  out
}

combine_cnmf_unit_llm_by_celltype <- function(index_df,
                                              run_root,
                                              output_subdir,
                                              level_tag = c('l2', 'l3'),
                                              celltype_field = c('celltype_l2', 'celltype_l3')) {
  level_tag <- tolower(as.character(level_tag)[1])
  celltype_field <- as.character(celltype_field)[1]
  empty_result <- list(
    master_manifest_df = data.frame(stringsAsFactors = FALSE, check.names = FALSE),
    master_manifest_tsv = NA_character_,
    group_results = list()
  )
  records <- unit_llm_index_records(index_df)
  if (!is.data.frame(records) || nrow(records) == 0L) return(empty_result)
  keep <- records$source_method == 'cnmf' & !is.na(records[[celltype_field]]) & nzchar(records[[celltype_field]])
  records <- records[keep, , drop = FALSE]
  if (nrow(records) == 0L) return(empty_result)

  split_key <- paste(records$lineage, records$safe_celltype, records[[celltype_field]], sep = '||')
  groups <- split(records, split_key)
  group_results <- list()
  master_rows <- list()

  for (group_name in names(groups)) {
    group_df <- groups[[group_name]]
    if (!is.data.frame(group_df) || nrow(group_df) == 0L) next
    group_df <- group_df[order(group_df$k_sort, group_df$gep_sort, group_df$unit_id), , drop = FALSE]
    lineage <- as.character(group_df$lineage[1])
    safe_celltype <- as.character(group_df$safe_celltype[1])
    celltype_value <- as.character(group_df[[celltype_field]][1])
    group_dir <- file.path(run_root, lineage, output_subdir, 'llm_parallel', 'cnmf', 'celltypes', safe_celltype)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    prefix <- sprintf('cnmf_%s_%s', level_tag, safe_celltype)
    combined_md <- file.path(group_dir, sprintf('%s_LLM_combined.md', prefix))
    combined_manifest_tsv <- file.path(group_dir, sprintf('%s_LLM_combined_manifest.tsv', prefix))
    combined_manifest_json <- file.path(group_dir, sprintf('%s_LLM_combined_manifest.json', prefix))

    manifest_df <- group_df[, c(
      'lineage', 'source_method', 'unit_type', 'unit_id', 'unit_label',
      'celltype_l1', 'celltype_l2', 'celltype_l3', 'celltype_hierarchy', 'safe_celltype',
      'k', 'gep', 'status', 'error', 'source_dir', 'unit_output_dir',
      'payload_json', 'prompt_md', 'interpretation_md', 'status_json', 'started_at', 'ended_at'
    ), drop = FALSE]
    write_tsv_unit_llm(manifest_df, combined_manifest_tsv)

    ok_units <- sum(group_df$status == 'ok', na.rm = TRUE)
    error_units <- sum(group_df$status == 'error', na.rm = TRUE)
    combined_manifest <- list(
      lineage = lineage,
      source_method = 'cnmf',
      level_tag = level_tag,
      celltype_field = celltype_field,
      safe_celltype = safe_celltype,
      celltype_value = celltype_value,
      celltype_l1 = if (!is.na(group_df$celltype_l1[1]) && nzchar(group_df$celltype_l1[1])) group_df$celltype_l1[1] else NULL,
      celltype_l2 = if (!is.na(group_df$celltype_l2[1]) && nzchar(group_df$celltype_l2[1])) group_df$celltype_l2[1] else NULL,
      celltype_l3 = if (!is.na(group_df$celltype_l3[1]) && nzchar(group_df$celltype_l3[1])) group_df$celltype_l3[1] else NULL,
      total_units = nrow(group_df),
      ok_units = ok_units,
      error_units = error_units,
      combined_md = normalizePath(combined_md, winslash = '/', mustWork = FALSE),
      combined_manifest_tsv = normalizePath(combined_manifest_tsv, winslash = '/', mustWork = FALSE),
      units = split(manifest_df, seq_len(nrow(manifest_df)))
    )
    write_json_unit_llm(combined_manifest, combined_manifest_json)

    status_summary <- table(factor(group_df$status, levels = sort(unique(group_df$status))))
    status_summary <- status_summary[status_summary > 0]
    metadata_lines <- c(
      sprintf('- lineage: `%s`', lineage),
      sprintf('- level: `%s`', toupper(level_tag)),
      sprintf('- source_method: `%s`', 'cnmf'),
      sprintf('- %s: `%s`', celltype_field, celltype_value),
      sprintf('- safe_celltype: `%s`', safe_celltype),
      sprintf('- total_units: `%d`', nrow(group_df)),
      sprintf('- ok_units: `%d`', ok_units),
      sprintf('- error_units: `%d`', error_units),
      sprintf('- combined_manifest_tsv: `%s`', normalizePath(combined_manifest_tsv, winslash = '/', mustWork = FALSE)),
      sprintf('- combined_manifest_json: `%s`', normalizePath(combined_manifest_json, winslash = '/', mustWork = FALSE))
    )
    if (!is.na(group_df$celltype_l1[1]) && nzchar(group_df$celltype_l1[1])) metadata_lines <- c(metadata_lines, sprintf('- celltype_l1: `%s`', group_df$celltype_l1[1]))
    if (!is.na(group_df$celltype_l2[1]) && nzchar(group_df$celltype_l2[1])) metadata_lines <- c(metadata_lines, sprintf('- celltype_l2: `%s`', group_df$celltype_l2[1]))
    if (!is.na(group_df$celltype_l3[1]) && nzchar(group_df$celltype_l3[1])) metadata_lines <- c(metadata_lines, sprintf('- celltype_l3: `%s`', group_df$celltype_l3[1]))

    summary_df <- data.frame(
      unit_id = group_df$unit_id,
      unit_label = group_df$unit_label,
      k = ifelse(is.na(group_df$k), '', as.character(group_df$k)),
      gep = group_df$gep,
      status = group_df$status,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    lines <- c(
      sprintf('# %s / %s / %s cNMF LLM 合并汇总', lineage, toupper(level_tag), celltype_value),
      '',
      metadata_lines,
      '',
      '## 状态概览',
      ''
    )
    if (length(status_summary) == 0L) {
      lines <- c(lines, '- 无可用 unit 状态。', '')
    } else {
      lines <- c(lines, sprintf('- `%s`: %d', names(status_summary), as.integer(status_summary)), '')
    }
    lines <- c(lines, '## Unit 清单', '', markdown_table_lines_unit_llm(summary_df, c('unit_id', 'unit_label', 'k', 'gep', 'status')), '')

    for (i in seq_len(nrow(group_df))) {
      row <- group_df[i, , drop = FALSE]
      unit_header <- sprintf('## %s (`%s`)', row$unit_label[[1]], row$unit_id[[1]])
      unit_meta <- c(
        sprintf('- status: `%s`', row$status[[1]]),
        sprintf('- interpretation_md: `%s`', row$interpretation_md[[1]]),
        sprintf('- status_json: `%s`', row$status_json[[1]])
      )
      if (!is.na(row$error[[1]]) && nzchar(row$error[[1]])) unit_meta <- c(unit_meta, sprintf('- error: `%s`', row$error[[1]]))
      if (!is.na(row$k[[1]])) unit_meta <- c(unit_meta, sprintf('- k: `%s`', row$k[[1]]))
      if (!is.na(row$gep[[1]]) && nzchar(row$gep[[1]])) unit_meta <- c(unit_meta, sprintf('- gep: `%s`', row$gep[[1]]))
      unit_body <- read_unit_llm_markdown_body(row$interpretation_md[[1]])
      lines <- c(lines, unit_header, '', unit_meta, '', unit_body, '')
    }
    write_lines_unit_llm(lines, combined_md)

    group_results[[group_name]] <- list(
      manifest_df = manifest_df,
      combined_md = normalizePath(combined_md, winslash = '/', mustWork = FALSE),
      combined_manifest_tsv = normalizePath(combined_manifest_tsv, winslash = '/', mustWork = FALSE),
      combined_manifest_json = normalizePath(combined_manifest_json, winslash = '/', mustWork = FALSE)
    )
    master_rows[[length(master_rows) + 1L]] <- data.frame(
      lineage = lineage,
      level_tag = level_tag,
      celltype_field = celltype_field,
      celltype_value = celltype_value,
      safe_celltype = safe_celltype,
      total_units = nrow(group_df),
      ok_units = ok_units,
      error_units = error_units,
      combined_md = normalizePath(combined_md, winslash = '/', mustWork = FALSE),
      combined_manifest_tsv = normalizePath(combined_manifest_tsv, winslash = '/', mustWork = FALSE),
      combined_manifest_json = normalizePath(combined_manifest_json, winslash = '/', mustWork = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  master_manifest_df <- if (length(master_rows) == 0L) {
    data.frame(stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    do.call(rbind, master_rows)
  }
  scanned_master_manifest_df <- collect_existing_cnmf_celltype_master_rows(
    run_root = run_root,
    output_subdir = output_subdir,
    level_tag = level_tag
  )
  if (is.data.frame(scanned_master_manifest_df) && nrow(scanned_master_manifest_df) > 0L) {
    master_manifest_df <- scanned_master_manifest_df
  }
  master_manifest_tsv <- if (nrow(master_manifest_df) > 0L) file.path(run_root, sprintf('cnmf_%s_celltype_LLM_combined_manifest.tsv', level_tag)) else NA_character_
  if (!is.na(master_manifest_tsv) && nrow(master_manifest_df) > 0L) write_tsv_unit_llm(master_manifest_df, master_manifest_tsv)

  list(
    master_manifest_df = master_manifest_df,
    master_manifest_tsv = if (!is.na(master_manifest_tsv)) normalizePath(master_manifest_tsv, winslash = '/', mustWork = FALSE) else NA_character_,
    group_results = group_results
  )
}

program_unit_llm_main <- function() {
  run_root <- Sys.getenv('PROGRAM_RUN_ROOT', unset = RUN_ROOT_DEFAULT)
  lineages <- split_env('PROGRAM_UNIT_LLM_LINEAGES', character())
  source_methods <- split_env('PROGRAM_UNIT_LLM_SOURCE_METHODS', c('cnmf', 'pycogaps', 'covarnet'))
  all_cnmf_k <- !(Sys.getenv('PROGRAM_UNIT_LLM_ALL_CNMF_K', unset = '1') %in% c('0', 'false', 'FALSE', 'no', 'NO'))
  max_parallel <- suppressWarnings(as.integer(Sys.getenv('PROGRAM_UNIT_LLM_MAX_PARALLEL', unset = '8')))
  if (length(max_parallel) == 0L || is.na(max_parallel) || max_parallel <= 0L) max_parallel <- 8L
  force <- Sys.getenv('PROGRAM_UNIT_LLM_FORCE', unset = '0') %in% c('1', 'true', 'TRUE', 'yes', 'YES')
  enable_live <- !(Sys.getenv('PROGRAM_UNIT_LLM_ENABLE_LIVE', unset = '1') %in% c('0', 'false', 'FALSE', 'no', 'NO'))
  model <- Sys.getenv('PROGRAM_UNIT_LLM_MODEL', unset = 'deepseek-reasoner')
  timeout_sec <- suppressWarnings(as.integer(Sys.getenv('PROGRAM_UNIT_LLM_TIMEOUT_SEC', unset = '240')))
  if (length(timeout_sec) == 0L || is.na(timeout_sec) || timeout_sec <= 0L) timeout_sec <- 240L

  idx <- build_unit_llm_task_index(run_root = run_root, lineages = lineages, source_methods = source_methods, all_cnmf_k = all_cnmf_k)
  index_path <- file.path(run_root, sprintf('program_unit_llm_task_index_%s.tsv', format(Sys.time(), '%Y%m%d_%H%M%S')))
  if (exists('pa_write_tsv', mode = 'function')) pa_write_tsv(idx, index_path) else utils::write.table(idx, file = index_path, sep = '\t', row.names = FALSE, quote = FALSE)
  cat(sprintf('[INFO] Unit LLM task index: %s (%d tasks)\n', index_path, nrow(idx)))
  results <- run_unit_llm_tasks(idx, max_parallel = max_parallel, enable_live = enable_live, force = force, model = model, timeout_sec = timeout_sec)
  summary <- data.frame(
    unit_id = vapply(results, function(x) as.character(x$unit_id %||% NA_character_), character(1)),
    status = vapply(results, function(x) as.character(x$status %||% NA_character_), character(1)),
    status_json = vapply(results, function(x) as.character(x$status_json %||% NA_character_), character(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  summary_path <- file.path(run_root, sprintf('program_unit_llm_summary_%s.tsv', format(Sys.time(), '%Y%m%d_%H%M%S')))
  if (exists('pa_write_tsv', mode = 'function')) pa_write_tsv(summary, summary_path) else utils::write.table(summary, file = summary_path, sep = '\t', row.names = FALSE, quote = FALSE)
  cat(sprintf('[OK] Unit LLM launcher finished: %s\n', summary_path))
  invisible(list(index = idx, results = results, summary = summary, index_path = index_path, summary_path = summary_path))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

if (identical(sys.nframe(), 0L)) {
  program_unit_llm_main()
}