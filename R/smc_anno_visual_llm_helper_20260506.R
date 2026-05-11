#!/usr/bin/env Rscript
# ==============================================================================
# smc_anno visual + LLM companion helper (2026-05-06)
# ==============================================================================

smcanno_viz_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

smcanno_viz_safe_trim <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

smcanno_viz_null_coalesce <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

smcanno_viz_env_flag <- function(name, default = FALSE) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(trimws(val))) return(isTRUE(default))
  tolower(trimws(val)) %in% c("1", "true", "yes", "y", "on")
}

smcanno_viz_resolve_llm_mode <- function(llm_config = list()) {
  if (is.null(llm_config) || !is.list(llm_config)) llm_config <- list()
  mode <- smcanno_viz_to_scalar(smcanno_viz_null_coalesce(llm_config$mode, NULL))
  if (!nzchar(mode)) {
    mode <- Sys.getenv("SMCANNO_FIGURE_LLM_MODE", unset = "queued")
  }
  mode <- tolower(trimws(mode))
  aliases <- c(
    queue = "queued",
    batch = "queued",
    deferred = "queued",
    none = "offline",
    rule = "offline",
    rule_based = "offline",
    live = "online",
    sync = "online"
  )
  if (mode %in% names(aliases)) mode <- aliases[[mode]]
  if (!mode %in% c("queued", "offline", "online")) mode <- "queued"
  mode
}

smcanno_viz_to_scalar <- function(x, collapse = "; ") {
  if (is.null(x) || length(x) == 0L) return("")
  if (is.list(x) && !is.data.frame(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  x <- smcanno_viz_safe_trim(x)
  x <- x[nzchar(x)]
  if (length(x) == 0L) return("")
  paste(x, collapse = collapse)
}

smcanno_viz_stem <- function(path) {
  sub("\\.[^.]+$", "", path)
}

smcanno_viz_write_markdown <- function(lines, path) {
  smcanno_viz_ensure_dir(dirname(path))
  writeLines(as.character(lines), con = path, useBytes = TRUE)
  invisible(path)
}

smcanno_viz_write_csv <- function(df, path) {
  smcanno_viz_ensure_dir(dirname(path))
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.csv(df, file = path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

smcanno_viz_write_json <- function(x, path) {
  smcanno_viz_ensure_dir(dirname(path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    smcanno_viz_write_markdown(utils::capture.output(str(x)), path)
  }
  invisible(path)
}

smcanno_viz_read_json <- function(path) {
  if (!file.exists(path) || !requireNamespace("jsonlite", quietly = TRUE)) return(list())
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) list())
}

smcanno_viz_markdown_table <- function(df, max_rows = 10L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) return("No tabular rows available.")
  df <- utils::head(df, max_rows)
  df[] <- lapply(df, function(col) gsub("\\|", "/", smcanno_viz_safe_trim(col)))
  header <- paste(c("", colnames(df), ""), collapse = "|")
  sep <- paste(c("", rep("---", ncol(df)), ""), collapse = "|")
  rows <- apply(df, 1, function(row) paste(c("", row, ""), collapse = "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

smcanno_viz_numeric_summary_lines <- function(df, max_numeric_cols = 8L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  numeric_cols <- colnames(df)[vapply(df, function(col) {
    vals <- suppressWarnings(as.numeric(col))
    sum(is.finite(vals)) > 0L
  }, logical(1))]
  numeric_cols <- utils::head(numeric_cols, max_numeric_cols)
  if (length(numeric_cols) == 0L) return("- No numeric columns detected for summary.")
  unlist(lapply(numeric_cols, function(nm) {
    vals <- suppressWarnings(as.numeric(df[[nm]]))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) return(sprintf("- `%s`: no finite values", nm))
    sprintf(
      "- `%s`: n=%s, min=%.4g, median=%.4g, mean=%.4g, max=%.4g",
      nm, length(vals), min(vals), stats::median(vals), mean(vals), max(vals)
    )
  }), use.names = FALSE)
}

smcanno_viz_categorical_summary_lines <- function(df, max_cols = 6L, max_levels = 8L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L) return("- No categorical rows detected for summary.")
  non_numeric_cols <- colnames(df)[!vapply(df, function(col) {
    vals <- suppressWarnings(as.numeric(col))
    sum(is.finite(vals)) > 0L && sum(is.finite(vals)) >= max(1L, floor(0.5 * length(vals)))
  }, logical(1))]
  non_numeric_cols <- utils::head(non_numeric_cols, max_cols)
  if (length(non_numeric_cols) == 0L) return("- No categorical columns detected for summary.")
  unlist(lapply(non_numeric_cols, function(nm) {
    vals <- smcanno_viz_safe_trim(df[[nm]])
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0L) return(sprintf("- `%s`: no non-empty values", nm))
    tab <- sort(table(vals), decreasing = TRUE)
    top <- utils::head(tab, max_levels)
    sprintf("- `%s`: %s", nm, paste(sprintf("%s=%s", names(top), as.integer(top)), collapse = "; "))
  }), use.names = FALSE)
}

smcanno_viz_is_placeholder_secret <- function(x) {
  x <- smcanno_viz_safe_trim(x)
  if (length(x) == 0L || !nzchar(x[[1]])) return(TRUE)
  lowered <- tolower(x[[1]])
  lowered %in% c(
    "your-deepseek-api-key",
    "your_deepseek_api_key_here",
    "your_deepseek_api_key",
    "your-key",
    "replace_me",
    "changeme"
  ) || grepl("^your[-_a-z]*api[-_a-z]*key", lowered) || grepl("placeholder|replace|changeme", lowered)
}

smcanno_viz_load_env_file <- function(path, overwrite_placeholder = TRUE) {
  if (is.null(path) || length(path) == 0L || !file.exists(path[[1]])) return(invisible(FALSE))
  lines <- readLines(path[[1]], warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- trimws(gsub("^export\\s+", "", sub("=.*$", "", line)))
    val <- trimws(gsub("^['\"]|['\"]$", "", sub("^[^=]*=", "", line)))
    if (!nzchar(key)) next
    cur <- Sys.getenv(key, unset = "")
    should_set <- !nzchar(cur)
    if (!should_set && isTRUE(overwrite_placeholder)) should_set <- smcanno_viz_is_placeholder_secret(cur)
    if (should_set) do.call(Sys.setenv, stats::setNames(list(val), key))
  }
  invisible(TRUE)
}

smcanno_viz_resolve_deepseek_key <- function(api_key = NULL,
                                             env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env")) {
  env_candidates <- as.character(env_candidates)
  env_candidates <- env_candidates[file.exists(env_candidates)]
  invisible(lapply(env_candidates, smcanno_viz_load_env_file))
  key <- if (is.null(api_key) || length(api_key) == 0L || !nzchar(smcanno_viz_safe_trim(api_key)[[1]])) {
    Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  } else {
    api_key[[1]]
  }
  key <- smcanno_viz_safe_trim(key)[[1]]
  has_live_key <- nchar(key) >= 10L && !smcanno_viz_is_placeholder_secret(key)
  list(api_key = if (has_live_key) key else "", has_live_key = has_live_key, env_files_loaded = env_candidates)
}

smcanno_viz_default_rule_analysis <- function(method,
                                              figure_type,
                                              title,
                                              data,
                                              extra_context = NULL) {
  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  context <- smcanno_viz_safe_trim(extra_context)
  context <- context[nzchar(context)]
  lines <- c(
    "## 自动规则解读",
    sprintf("- 图表：`%s`", title),
    sprintf("- 方法：`%s`; 图类型：`%s`", method, figure_type),
    sprintf("- 配套数据：`%s` 行 × `%s` 列。", nrow(df), ncol(df)),
    "",
    "### 数值字段概览",
    smcanno_viz_numeric_summary_lines(df),
    "",
    "### 分组/类别字段概览",
    smcanno_viz_categorical_summary_lines(df),
    "",
    "## 解读要点",
    if (length(context) > 0L) paste0("- ", context) else "- 请结合同目录 CSV 的具体数值、显著性和细胞数分布解读；不要只看图形外观。",
    "",
    "## 人工复核建议",
    "- 检查图中最强信号是否由少数细胞、单一样本或极端值驱动。",
    "- 对 trajectory / PAGA / tradeSeq 结果，优先复核与 Slingshot 终末分支和已知 SMC/Pericyte 标记是否一致。",
    "- 若后续要写入报告，请优先引用本图对应 CSV 中的数值，而不是仅凭视觉判断。"
  )
  lines
}

smcanno_viz_tradeSeq_specific_lines <- function(df, figure_type) {
  if (!grepl("tradeseq", tolower(figure_type))) return(character())
  character()
}

smcanno_viz_final_summary_lines <- function(method,
                                            figure_type,
                                            title,
                                            data,
                                            extra_context = NULL) {
  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  context <- smcanno_viz_safe_trim(extra_context)
  context <- context[nzchar(context)]
  key_numeric <- colnames(df)[vapply(df, function(col) {
    vals <- suppressWarnings(as.numeric(col))
    sum(is.finite(vals)) > 0L
  }, logical(1))]
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  dominant_group <- ""
  for (candidate in c("test_name", "lineage_label", "lineage", "cluster_label", "group", "gene")) {
    if (!candidate %in% colnames(df)) next
    vals <- smcanno_viz_safe_trim(df[[candidate]])
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0L) next
    tab <- sort(table(vals), decreasing = TRUE)
    dominant_group <- sprintf("`%s` 中最常见的是 `%s` (n=%s)。", candidate, names(tab)[1], as.integer(tab[[1]]))
    break
  }
  c(
    "## 最后总结",
    sprintf("- 这张图应和其配套 CSV 一起解读：当前数据规模为 `%s` 行 × `%s` 列。", n_rows, n_cols),
    if (length(key_numeric) > 0L) sprintf("- 主要数值字段包括：`%s`。", paste(utils::head(key_numeric, 8L), collapse = "`, `")) else "- 未检测到可量化数值字段，主要用于类别/拓扑核查。",
    if (nzchar(dominant_group)) paste0("- ", dominant_group) else "- 没有明显单一主导分组字段，需回到 CSV 查看逐行证据。",
    if (length(context) > 0L) sprintf("- 关键上下文：%s", paste(context, collapse = " ")) else "- 未提供额外上下文；建议结合上游 trajectory、branch 和 marker 结果综合判断。",
    "- 结论应优先用于提出可复核假设，而不是单独作为最终生物学判定。"
  )
}

smcanno_viz_build_prompt <- function(method,
                                     figure_type,
                                     title,
                                     data,
                                     extra_context = NULL) {
  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  list(
    system_prompt = paste(
      "You are a single-cell trajectory and gene-program analysis expert.",
      "Interpret the figure using the paired CSV data, not visual impression alone.",
      "Respond in concise Chinese Markdown with caveats."
    ),
    user_prompt = paste(
      "请根据下面的图表配套数据和上下文，生成与该图对应的 LLM 分析。",
      "",
      sprintf("## Figure title\n%s", title),
      sprintf("## Method\n%s", method),
      sprintf("## Figure type\n%s", figure_type),
      "",
      "## Data dimensions",
      sprintf("- rows: %s", nrow(df)),
      sprintf("- columns: %s", paste(colnames(df), collapse = ", ")),
      "",
      "## Numeric summaries",
      paste(smcanno_viz_numeric_summary_lines(df), collapse = "\n"),
      "",
      "## Categorical summaries",
      paste(smcanno_viz_categorical_summary_lines(df), collapse = "\n"),
      "",
      "## Preview rows",
      smcanno_viz_markdown_table(df, max_rows = 12L),
      "",
      if (!is.null(extra_context) && length(extra_context) > 0L) paste("## Extra context\n", paste(smcanno_viz_safe_trim(extra_context), collapse = "\n"), sep = "") else "",
      "",
      "请输出深入分析：1) 主要发现；2) 哪些结论由 CSV 数值直接支持；3) 是否支持轨迹/分支解释；4) 风险和混杂；5) 最值得复核的基因/边/分支；6) 最后总结。不要只复述预览行，要综合 numeric summaries、categorical summaries 和上下文。",
      sep = "\n"
    )
  )
}

smcanno_viz_read_prompt_sections <- function(prompt_path) {
  lines <- readLines(prompt_path, warn = FALSE)
  system_idx <- which(trimws(lines) == "## System")
  user_idx <- which(trimws(lines) == "## User")
  if (length(system_idx) == 0L || length(user_idx) == 0L || user_idx[[1]] <= system_idx[[1]]) {
    return(list(system_prompt = "You are a rigorous single-cell analysis expert. Respond in Chinese Markdown.", user_prompt = paste(lines, collapse = "\n")))
  }
  s0 <- system_idx[[1]] + 1L
  s1 <- user_idx[[1]] - 1L
  u0 <- user_idx[[1]] + 1L
  list(
    system_prompt = paste(lines[s0:s1], collapse = "\n"),
    user_prompt = paste(lines[u0:length(lines)], collapse = "\n")
  )
}

smcanno_viz_deepseek_chat <- function(prompt,
                                      system_prompt,
                                      model,
                                      api_key,
                                      timeout_sec = 180,
                                      base_url = NULL,
                                      thinking = NULL,
                                      reasoning_effort = NULL) {
  direct_http <- function() {
    if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Packages 'httr' and 'jsonlite' are required for direct DeepSeek HTTP calls.", call. = FALSE)
    }
    base_url_use <- smcanno_viz_to_scalar(smcanno_viz_null_coalesce(base_url, "https://api.deepseek.com"))
    base_url_use <- sub("/+$", "", base_url_use)
    endpoint <- if (grepl("/chat/completions$", base_url_use)) base_url_use else paste0(base_url_use, "/chat/completions")
    messages <- list()
    if (nzchar(smcanno_viz_to_scalar(system_prompt))) {
      messages[[length(messages) + 1L]] <- list(role = "system", content = smcanno_viz_to_scalar(system_prompt))
    }
    messages[[length(messages) + 1L]] <- list(role = "user", content = smcanno_viz_to_scalar(prompt))
    payload <- list(model = model, messages = messages)
    if (!is.null(thinking)) payload$thinking <- thinking
    if (!is.null(reasoning_effort) && nzchar(smcanno_viz_to_scalar(reasoning_effort))) payload$reasoning_effort <- smcanno_viz_to_scalar(reasoning_effort)
    resp <- httr::POST(
      url = endpoint,
      httr::add_headers(Authorization = paste("Bearer", api_key), `Content-Type` = "application/json"),
      body = payload,
      encode = "json",
      httr::timeout(as.numeric(timeout_sec))
    )
    body_txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    code <- httr::status_code(resp)
    if (code < 200L || code >= 300L) {
      stop(sprintf("DeepSeek HTTP %s: %s", code, substr(body_txt, 1L, 1000L)), call. = FALSE)
    }
    parsed <- jsonlite::fromJSON(body_txt, simplifyVector = FALSE)
    choice <- parsed$choices[[1]]
    content <- smcanno_viz_to_scalar(choice$message$content)
    if (!nzchar(content) && !is.null(choice$message$reasoning_content)) {
      content <- smcanno_viz_to_scalar(choice$message$reasoning_content)
    }
    if (!nzchar(content)) stop("DeepSeek response did not contain message content.", call. = FALSE)
    content
  }
  call_tc <- function(fun) {
    args <- list(
      prompt = prompt,
      system_prompt = system_prompt,
      model = model,
      api_key = api_key,
      timeout_sec = timeout_sec,
      base_url = base_url,
      thinking = thinking,
      reasoning_effort = reasoning_effort
    )
    fml <- names(formals(fun))
    if (!"..." %in% fml) args <- args[names(args) %in% fml]
    do.call(fun, args)
  }
  if (exists("tc_deepseek_chat_request", mode = "function")) {
    return(call_tc(tc_deepseek_chat_request))
  }
  helper_path <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
  if (file.exists(helper_path)) {
    source(helper_path)
    if (exists("tc_deepseek_chat_request", mode = "function")) {
      return(call_tc(tc_deepseek_chat_request))
    }
  }
  http_response <- tryCatch(direct_http(), error = function(e) structure(conditionMessage(e), class = "smcanno_viz_direct_http_error"))
  if (!inherits(http_response, "smcanno_viz_direct_http_error")) return(http_response)
  if (!requireNamespace("fanyi", quietly = TRUE)) {
    stop(sprintf("Direct DeepSeek HTTP failed and fanyi is unavailable: %s", as.character(http_response[[1]])), call. = FALSE)
  }
  try(fanyi::set_translate_option(key = api_key, source = "deepseek"), silent = TRUE)
  fanyi::chat_request(
    paste(system_prompt, prompt, sep = "\n\n"),
    model = model,
    api_key = api_key
  )
}

smcanno_viz_append_manifest <- function(row, manifest_path) {
  smcanno_viz_ensure_dir(dirname(manifest_path))
  row <- as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.table(
    row,
    file = manifest_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = !file.exists(manifest_path),
    append = file.exists(manifest_path),
    na = ""
  )
  invisible(manifest_path)
}

smcanno_viz_collect_companion_manifests <- function(root_dir) {
  if (is.null(root_dir) || length(root_dir) == 0L || !dir.exists(root_dir[[1]])) return(character())
  files <- list.files(root_dir[[1]], pattern = "(figure_companion_manifest|llm_manifest)\\.tsv$", recursive = TRUE, full.names = TRUE)
  normalizePath(files[file.exists(files)], winslash = "/", mustWork = FALSE)
}

smcanno_viz_write_manifest_index <- function(manifest_paths, output_path) {
  manifest_paths <- normalizePath(unique(as.character(manifest_paths)), winslash = "/", mustWork = FALSE)
  manifest_paths <- manifest_paths[file.exists(manifest_paths)]
  df <- data.frame(
    manifest_path = manifest_paths,
    n_rows = vapply(manifest_paths, function(path) {
      tryCatch(nrow(utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)), error = function(e) NA_integer_)
    }, integer(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  smcanno_viz_write_csv(df, output_path)
  invisible(output_path)
}

smcanno_viz_launch_llm_batch <- function(manifest_paths = character(),
                                         output_dir = NULL,
                                         runner_script = "/home/h2048/script/R/smc_anno_figure_llm_batch_20260507.R",
                                         rscript_bin = Sys.which("Rscript"),
                                         async = TRUE,
                                         model = NULL,
                                         timeout_sec = 180,
                                         log_path = NULL,
                                         max_items = NULL,
                                         force = FALSE) {
  if (!file.exists(runner_script)) stop(sprintf("LLM runner script not found: %s", runner_script), call. = FALSE)
  if (is.null(rscript_bin) || !nzchar(rscript_bin)) rscript_bin <- "/usr/bin/Rscript"
  args <- c(runner_script)
  manifest_paths <- unique(as.character(manifest_paths))
  manifest_paths <- manifest_paths[file.exists(manifest_paths)]
  for (manifest_path in manifest_paths) args <- c(args, "--manifest", normalizePath(manifest_path, winslash = "/", mustWork = FALSE))
  if (!is.null(output_dir) && nzchar(smcanno_viz_to_scalar(output_dir))) args <- c(args, "--output-dir", output_dir[[1]])
  if (!is.null(model) && nzchar(smcanno_viz_to_scalar(model))) args <- c(args, "--model", smcanno_viz_to_scalar(model))
  if (!is.null(timeout_sec) && is.finite(as.numeric(timeout_sec))) args <- c(args, "--timeout-sec", as.character(as.numeric(timeout_sec)))
  if (!is.null(max_items) && is.finite(as.numeric(max_items))) args <- c(args, "--max-items", as.character(as.integer(max_items)))
  if (isTRUE(force)) args <- c(args, "--force")
  if (is.null(log_path) || !nzchar(smcanno_viz_to_scalar(log_path))) {
    log_path <- if (!is.null(output_dir) && nzchar(smcanno_viz_to_scalar(output_dir))) file.path(output_dir[[1]], "figure_llm_batch.log") else tempfile("figure_llm_batch_", fileext = ".log")
  }
  smcanno_viz_ensure_dir(dirname(log_path))
  status <- system2(rscript_bin, args = args, stdout = log_path, stderr = log_path, wait = !isTRUE(async))
  list(
    runner_script = normalizePath(runner_script, winslash = "/", mustWork = FALSE),
    rscript_bin = rscript_bin,
    args = args,
    async = isTRUE(async),
    log_path = normalizePath(log_path, winslash = "/", mustWork = FALSE),
    status = if (isTRUE(async)) "started_async" else as.character(status)
  )
}

smcanno_viz_register_figure <- function(figure_path,
                                        data,
                                        method,
                                        figure_type,
                                        title,
                                        extra_context = NULL,
                                        llm_config = list(),
                                        manifest_path = NULL,
                                        data_path = NULL) {
  figure_path <- normalizePath(figure_path, winslash = "/", mustWork = FALSE)
  stem <- smcanno_viz_stem(figure_path)
  if (is.null(data_path) || !nzchar(smcanno_viz_to_scalar(data_path))) data_path <- paste0(stem, "_data.csv")
  prompt_path <- paste0(stem, "_LLM_PROMPT.md")
  analysis_path <- paste0(stem, "_LLM_ANALYSIS.md")
  status_path <- paste0(stem, "_LLM_STATUS.json")
  if (is.null(manifest_path) || !nzchar(smcanno_viz_to_scalar(manifest_path))) {
    manifest_path <- file.path(dirname(figure_path), "figure_companion_manifest.tsv")
  }

  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  smcanno_viz_write_csv(df, data_path)
  prompt <- smcanno_viz_build_prompt(method = method, figure_type = figure_type, title = title, data = df, extra_context = extra_context)
  smcanno_viz_write_markdown(c("# Figure LLM prompt", "", "## System", prompt$system_prompt, "", "## User", prompt$user_prompt), prompt_path)
  rule_lines <- smcanno_viz_default_rule_analysis(method, figure_type, title, df, extra_context = extra_context)
  final_lines <- smcanno_viz_final_summary_lines(method, figure_type, title, df, extra_context = extra_context)

  if (is.null(llm_config) || !is.list(llm_config)) llm_config <- list()
  enabled <- if (is.null(llm_config$enabled)) TRUE else isTRUE(llm_config$enabled)
  llm_mode <- smcanno_viz_resolve_llm_mode(llm_config)
  model <- smcanno_viz_to_scalar(if (is.null(llm_config$model)) "deepseek-reasoner" else llm_config$model)
  if (!nzchar(model)) model <- "deepseek-reasoner"
  key_info <- if (identical(llm_mode, "online") && isTRUE(enabled)) {
    smcanno_viz_resolve_deepseek_key(api_key = llm_config$api_key)
  } else {
    list(api_key = "", has_live_key = FALSE, env_files_loaded = character())
  }
  status <- list(
    enabled = enabled,
    status = "skipped_disabled",
    llm_mode = llm_mode,
    method = method,
    figure_type = figure_type,
    title = title,
    model = model,
    figure_path = figure_path,
    data_csv = normalizePath(data_path, winslash = "/", mustWork = FALSE),
    prompt_md = prompt_path,
    analysis_md = analysis_path,
    env_files_loaded = key_info$env_files_loaded,
    error = NULL
  )

  if (!enabled) {
    status$status <- "skipped_disabled"
    smcanno_viz_write_markdown(c("# Figure LLM analysis", "", "**LLM status:** `skipped_disabled`", "", rule_lines, "", final_lines), analysis_path)
  } else if (identical(llm_mode, "queued")) {
    status$status <- "queued_for_batch"
    status$runner_script <- "/home/h2048/script/R/smc_anno_figure_llm_batch_20260507.R"
    smcanno_viz_write_markdown(c(
      "# Figure LLM analysis",
      "",
      "**LLM status:** `queued_for_batch`（已生成配套 CSV 与 prompt；在线 LLM 已从画图流程解耦，请由批处理脚本异步/单独补写。）",
      "",
      "## 数据驱动审计",
      "",
      rule_lines,
      "",
      final_lines,
      "",
      "## 在线 LLM 补写方式",
      sprintf("- Runner: `%s`", status$runner_script),
      sprintf("- Manifest: `%s`", normalizePath(manifest_path, winslash = "/", mustWork = FALSE)),
      "- 该 runner 会读取本图 `*_data.csv` 与 `*_LLM_PROMPT.md`，成功后覆盖本文件中的在线解读段并更新 `*_LLM_STATUS.json`。"
    ), analysis_path)
  } else if (identical(llm_mode, "offline")) {
    status$status <- "skipped_offline"
    smcanno_viz_write_markdown(c("# Figure LLM analysis", "", "**LLM status:** `skipped_offline`（仅生成规则解读；未排队在线 LLM。）", "", rule_lines, "", final_lines), analysis_path)
  } else if (!isTRUE(key_info$has_live_key)) {
    status$status <- "skipped_no_live_key"
    smcanno_viz_write_markdown(c(
      "# Figure LLM analysis",
      "",
      "**LLM status:** `skipped_no_live_key`（未检测到真实 `DEEPSEEK_API_KEY`，已生成可复用 prompt 和规则解读。）",
      "",
      rule_lines,
      "",
      final_lines,
      "",
      "## 如何生成在线 LLM 版本",
      "在 `/home/h2048/.env` 写入真实 `DEEPSEEK_API_KEY` 后重跑对应脚本。"
    ), analysis_path)
  } else {
    response <- tryCatch({
      smcanno_viz_deepseek_chat(
        prompt = prompt$user_prompt,
        system_prompt = prompt$system_prompt,
        model = model,
        api_key = key_info$api_key,
        timeout_sec = if (is.null(llm_config$timeout_sec)) 180 else as.numeric(llm_config$timeout_sec),
        base_url = llm_config$base_url,
        thinking = llm_config$thinking,
        reasoning_effort = llm_config$reasoning_effort
      )
    }, error = function(e) structure(conditionMessage(e), class = "smcanno_viz_llm_error"))
    if (inherits(response, "smcanno_viz_llm_error")) {
      status$status <- "error"
      status$error <- as.character(response[[1]])
      smcanno_viz_write_markdown(c("# Figure LLM analysis", "", "**LLM status:** `error`", "", sprintf("LLM error: `%s`", status$error), "", rule_lines, "", final_lines), analysis_path)
    } else {
      status$status <- "ok"
      smcanno_viz_write_markdown(c("# Figure LLM analysis", "", "**LLM status:** `ok`", "", "## 数据驱动审计", "", rule_lines, "", "## 在线 LLM 深度解读", "", response, "", final_lines), analysis_path)
    }
  }

  smcanno_viz_write_json(status, status_path)
  row <- data.frame(
    method = method,
    figure_type = figure_type,
    title = title,
    figure_path = figure_path,
    data_csv = normalizePath(data_path, winslash = "/", mustWork = FALSE),
    llm_analysis_md = normalizePath(analysis_path, winslash = "/", mustWork = FALSE),
    llm_prompt_md = normalizePath(prompt_path, winslash = "/", mustWork = FALSE),
    llm_status_json = normalizePath(status_path, winslash = "/", mustWork = FALSE),
    llm_mode = llm_mode,
    llm_status = status$status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  smcanno_viz_append_manifest(row, manifest_path)
  invisible(c(status, list(manifest_path = manifest_path)))
}

if (sys.nframe() == 0) {
  cat("smc_anno visual + LLM companion helper (2026-05-06) loaded.\n")
}