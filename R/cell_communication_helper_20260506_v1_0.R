#!/usr/bin/env Rscript
# ==============================================================================
# Cell-Cell Communication Helper (cci_*)
# cell_communication_helper_20260506_v1_0.R
# ==============================================================================
#
# Purpose:
#   Helper functions for two cell communication analysis scenarios:
#
#   Scenario A - scRNA only (no spatial, no protein validation):
#     CellPhoneDB / CellChat  ->  candidate LR screening
#     LIANA+                  ->  multi-method consensus
#     NicheNet                ->  receiver-side downstream response
#     Language rule: putative / candidate / potential - never "direct interaction"
#
#   Scenario B - Disease vs Control with biological replicates:
#     MultiNicheNet           ->  sample-level DE-based CCI (primary)
#     LIANA+                  ->  auxiliary per-sample aggregation
#     CellChat                ->  global network visualization
#     Core rule: samples are the statistical unit, NOT cells
#
# Output convention:
#   - Every exported interpretation table gets a paired CSV for LLM review.
#   - Every generated figure gets *_data.csv, *_LLM_PROMPT.md,
#     *_LLM_ANALYSIS.md, *_LLM_STATUS.json, and a companion manifest row.
#
# Prefix:  cci_
# Version: v1.1  (2026-05-06)
# ==============================================================================

# ---------------------------------------------------------------------------- #
#  Section 0  Utilities                                                         #
# ---------------------------------------------------------------------------- #

cci_trim <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x[[1]])
}

cci_safe_mkdir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

cci_sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}

cci_to_scalar <- function(x, collapse = "; ") {
  if (is.null(x) || length(x) == 0L) return("")
  if (is.list(x) && !is.data.frame(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0L) return("")
  paste(x, collapse = collapse)
}

cci_null_coalesce <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

cci_path_normalize <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

cci_stem <- function(path) {
  sub("\\.[^.]+$", "", path)
}

cci_write_markdown <- function(lines, path) {
  cci_safe_mkdir(dirname(path))
  writeLines(as.character(lines), con = path, useBytes = TRUE)
  invisible(path)
}

cci_write_csv <- function(df, path) {
  cci_safe_mkdir(dirname(path))
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.csv(df, file = path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

cci_write_tsv <- function(df, path) {
  cci_safe_mkdir(dirname(path))
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  invisible(path)
}

cci_write_json <- function(x, path) {
  cci_safe_mkdir(dirname(path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    cci_write_markdown(utils::capture.output(str(x)), path)
  }
  invisible(path)
}

cci_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

cci_pick_existing_col <- function(df, candidates) {
  if (is.null(df) || !is.data.frame(df) || ncol(df) == 0L) return(NULL)
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) return(NULL)
  hit[[1]]
}

cci_markdown_table <- function(df, max_rows = 12L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) return("No tabular rows available.")
  df <- utils::head(df, max_rows)
  df[] <- lapply(df, function(col) gsub("\\|", "/", cci_to_scalar_list(col)))
  header <- paste(c("", colnames(df), ""), collapse = "|")
  sep <- paste(c("", rep("---", ncol(df)), ""), collapse = "|")
  rows <- apply(df, 1, function(row) paste(c("", row, ""), collapse = "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

cci_to_scalar_list <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character())
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

cci_numeric_summary_lines <- function(df, max_numeric_cols = 10L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  numeric_cols <- colnames(df)[vapply(df, function(col) {
    vals <- cci_num(col)
    sum(is.finite(vals)) > 0L
  }, logical(1))]
  numeric_cols <- utils::head(numeric_cols, max_numeric_cols)
  if (length(numeric_cols) == 0L) return("- No numeric columns detected.")
  unlist(lapply(numeric_cols, function(nm) {
    vals <- cci_num(df[[nm]])
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) return(sprintf("- `%s`: no finite values", nm))
    sprintf(
      "- `%s`: n=%s, min=%.4g, median=%.4g, mean=%.4g, max=%.4g",
      nm, length(vals), min(vals), stats::median(vals), mean(vals), max(vals)
    )
  }), use.names = FALSE)
}

cci_categorical_summary_lines <- function(df, max_cols = 8L, max_levels = 10L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L) return("- No categorical rows detected.")
  non_numeric_cols <- colnames(df)[!vapply(df, function(col) {
    vals <- cci_num(col)
    sum(is.finite(vals)) > 0L && sum(is.finite(vals)) >= max(1L, floor(0.5 * length(vals)))
  }, logical(1))]
  non_numeric_cols <- utils::head(non_numeric_cols, max_cols)
  if (length(non_numeric_cols) == 0L) return("- No categorical columns detected.")
  unlist(lapply(non_numeric_cols, function(nm) {
    vals <- cci_to_scalar_list(df[[nm]])
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0L) return(sprintf("- `%s`: no non-empty values", nm))
    tab <- sort(table(vals), decreasing = TRUE)
    top <- utils::head(tab, max_levels)
    sprintf("- `%s`: %s", nm, paste(sprintf("%s=%s", names(top), as.integer(top)), collapse = "; "))
  }), use.names = FALSE)
}

cci_is_placeholder_secret <- function(x) {
  x <- cci_to_scalar(x)
  if (!nzchar(x)) return(TRUE)
  lowered <- tolower(x)
  lowered %in% c(
    "your-deepseek-api-key",
    "your_deepseek_api_key_here",
    "your_deepseek_api_key",
    "your-key",
    "replace_me",
    "changeme"
  ) || grepl("^your[-_a-z]*api[-_a-z]*key", lowered) || grepl("placeholder|replace|changeme", lowered)
}

cci_load_env_file <- function(path, overwrite_placeholder = TRUE) {
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
    if (!should_set && isTRUE(overwrite_placeholder)) should_set <- cci_is_placeholder_secret(cur)
    if (should_set) do.call(Sys.setenv, stats::setNames(list(val), key))
  }
  invisible(TRUE)
}

cci_ensure_env_placeholder <- function(path = "/home/h2048/.env", key = "DEEPSEEK_API_KEY") {
  cci_safe_mkdir(dirname(path))
  line <- sprintf("%s=your_deepseek_api_key_here", key)
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    if (any(grepl(sprintf("^%s\\s*=", key), trimws(lines)))) return(invisible(FALSE))
    write(line, file = path, append = TRUE)
    return(invisible(TRUE))
  }
  writeLines(line, path)
  invisible(TRUE)
}

cci_resolve_deepseek_key <- function(api_key = NULL,
                                     env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env"),
                                     write_env_placeholder = FALSE) {
  env_candidates <- as.character(env_candidates)
  env_candidates <- env_candidates[file.exists(env_candidates)]
  invisible(lapply(env_candidates, cci_load_env_file))
  key <- if (is.null(api_key) || length(api_key) == 0L || !nzchar(cci_to_scalar(api_key))) {
    Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  } else {
    api_key[[1]]
  }
  key <- cci_to_scalar(key)
  has_live_key <- nchar(key) >= 10L && !cci_is_placeholder_secret(key)
  if (!has_live_key && isTRUE(write_env_placeholder)) {
    cci_ensure_env_placeholder("/home/h2048/.env", "DEEPSEEK_API_KEY")
  }
  list(
    api_key = if (has_live_key) key else "",
    has_live_key = has_live_key,
    env_files_loaded = env_candidates,
    status = if (has_live_key) "ready" else "missing_or_placeholder_key"
  )
}

cci_deepseek_chat <- function(prompt,
                              system_prompt,
                              model,
                              api_key,
                              timeout_sec = 180,
                              base_url = NULL,
                              thinking = NULL,
                              reasoning_effort = NULL) {
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
  if (!requireNamespace("fanyi", quietly = TRUE)) {
    stop("Neither tc_deepseek_chat_request() nor fanyi is available for live LLM calls.", call. = FALSE)
  }
  fanyi::chat_request(
    paste(system_prompt, prompt, sep = "\n\n"),
    model = model,
    api_key = api_key
  )
}

cci_rule_based_llm_analysis <- function(df,
                                        title,
                                        scenario,
                                        extra_context = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  context <- cci_to_scalar_list(extra_context)
  context <- context[nzchar(context)]
  c(
    "## 自动规则解读",
    sprintf("- 分析对象：`%s`", title),
    sprintf("- 场景：`%s`", scenario),
    sprintf("- 配套 CSV：`%s` 行 × `%s` 列。", nrow(df), ncol(df)),
    "",
    "### 数值字段概览",
    cci_numeric_summary_lines(df),
    "",
    "### 分组/类别字段概览",
    cci_categorical_summary_lines(df),
    "",
    "### 预览行",
    cci_markdown_table(df, max_rows = 10L),
    "",
    "## CCI 解读护栏",
    "- 普通解离 scRNA-seq 只能支持 computationally inferred / putative / candidate communication。",
    "- 没有空间共定位或蛋白层验证时，不写 confirmed/direct signaling。",
    "- 疾病 vs 对照场景必须把样本/供体作为统计重复；细胞数不能替代生物学重复。",
    "- LIANA+/CellChat/CellPhoneDB 更适合候选与共识筛选；MultiNicheNet 的样本级 DE 证据优先用于组间差异 CCI。",
    "",
    "## 人工复核建议",
    if (length(context) > 0L) paste0("- ", context) else "- 回到配套 CSV 查看 top interaction 的 sender/receiver 细胞数、样本覆盖、排序分数和方法一致性。",
    "- 检查候选 LR 是否由少数样本、低细胞数群体、环境 RNA 或双细胞驱动。",
    "- 对 B 细胞/浆细胞相关结果，额外复核 IG/JCHAIN/MZB1/XBP1 等高表达是否代表真实状态还是 ambient/plasma spillover。"
  )
}

cci_build_llm_prompt <- function(df,
                                 title,
                                 scenario,
                                 extra_context = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  list(
    system_prompt = paste(
      "You are a rigorous single-cell cell-cell communication analysis expert.",
      "Interpret only what is supported by the paired CSV data.",
      "Use conservative language: potential, putative, candidate, computationally inferred.",
      "Do not claim direct interaction or confirmed signaling without spatial/protein validation.",
      "Respond in concise Simplified Chinese Markdown."
    ),
    user_prompt = paste(
      "请根据下面的 CCI 配套 CSV 摘要和上下文，生成可写入报告的保守解读。",
      "",
      sprintf("## Analysis title\n%s", title),
      sprintf("## Scenario\n%s", scenario),
      "",
      "## Data dimensions",
      sprintf("- rows: %s", nrow(df)),
      sprintf("- columns: %s", paste(colnames(df), collapse = ", ")),
      "",
      "## Numeric summaries",
      paste(cci_numeric_summary_lines(df), collapse = "\n"),
      "",
      "## Categorical summaries",
      paste(cci_categorical_summary_lines(df), collapse = "\n"),
      "",
      "## Preview rows",
      cci_markdown_table(df, max_rows = 12L),
      "",
      if (!is.null(extra_context) && length(extra_context) > 0L) paste("## Extra context\n", paste(cci_to_scalar_list(extra_context), collapse = "\n"), sep = "") else "",
      "",
      "请输出：1) 主要候选通信模式；2) 哪些结论直接由 CSV 支持；3) 若为疾病 vs 对照，明确样本是统计单位；4) 风险和混杂；5) 值得验证的 LR / sender / receiver；6) 最后总结。务必使用 potential/candidate/putative，不要写 direct interaction 或 confirmed signaling。",
      sep = "\n"
    )
  )
}

cci_run_table_llm_interpretation <- function(df,
                                             output_dir,
                                             prefix,
                                             title,
                                             scenario,
                                             extra_context = NULL,
                                             llm_config = list()) {
  cci_safe_mkdir(output_dir)
  if (is.null(llm_config) || !is.list(llm_config)) llm_config <- list()
  enabled <- if (is.null(llm_config$enabled)) TRUE else isTRUE(llm_config$enabled)
  model <- cci_to_scalar(cci_null_coalesce(llm_config$model, "deepseek-reasoner"))
  if (!nzchar(model)) model <- "deepseek-reasoner"
  timeout_sec <- suppressWarnings(as.numeric(cci_null_coalesce(llm_config$timeout_sec, 180)))
  if (is.na(timeout_sec) || timeout_sec <= 0) timeout_sec <- 180

  prompt <- cci_build_llm_prompt(df, title = title, scenario = scenario, extra_context = extra_context)
  prompt_path <- file.path(output_dir, sprintf("%s_LLM_PROMPT.md", prefix))
  analysis_path <- file.path(output_dir, sprintf("%s_LLM_ANALYSIS.md", prefix))
  raw_path <- file.path(output_dir, sprintf("%s_LLM_RAW_RESPONSE.txt", prefix))
  status_path <- file.path(output_dir, sprintf("%s_LLM_STATUS.json", prefix))

  cci_write_markdown(c("# CCI LLM prompt", "", "## System", prompt$system_prompt, "", "## User", prompt$user_prompt), prompt_path)

  key_info <- cci_resolve_deepseek_key(
    api_key = llm_config$api_key,
    env_candidates = cci_null_coalesce(llm_config$env_candidates, c("/home/h2048/.env", "/home/h2048/script/.env")),
    write_env_placeholder = isTRUE(cci_null_coalesce(llm_config$write_env_placeholder, FALSE))
  )
  status <- list(
    enabled = enabled,
    status = "skipped_disabled",
    title = title,
    scenario = scenario,
    model = model,
    prompt_md = cci_path_normalize(prompt_path),
    analysis_md = cci_path_normalize(analysis_path),
    raw_response_txt = cci_path_normalize(raw_path),
    env_files_loaded = key_info$env_files_loaded,
    error = NULL
  )

  rule_lines <- cci_rule_based_llm_analysis(df, title = title, scenario = scenario, extra_context = extra_context)
  if (!enabled) {
    status$status <- "skipped_disabled"
    cci_write_markdown(c("# CCI LLM analysis", "", "**LLM status:** `skipped_disabled`", "", rule_lines), analysis_path)
  } else if (!isTRUE(key_info$has_live_key)) {
    status$status <- "skipped_no_live_key"
    cci_write_markdown(c(
      "# CCI LLM analysis",
      "",
      "**LLM status:** `skipped_no_live_key`（未检测到真实 `DEEPSEEK_API_KEY`，已生成可复用 prompt 和规则解读。）",
      "",
      rule_lines,
      "",
      "## 如何生成在线 LLM 版本",
      "在 `/home/h2048/.env` 写入真实 `DEEPSEEK_API_KEY` 后重跑对应脚本，或传入 `llm_config = list(api_key = ...)`。"
    ), analysis_path)
  } else {
    response <- tryCatch({
      cci_deepseek_chat(
        prompt = prompt$user_prompt,
        system_prompt = prompt$system_prompt,
        model = model,
        api_key = key_info$api_key,
        timeout_sec = timeout_sec,
        base_url = llm_config$base_url,
        thinking = llm_config$thinking,
        reasoning_effort = llm_config$reasoning_effort
      )
    }, error = function(e) structure(conditionMessage(e), class = "cci_llm_error"))
    if (inherits(response, "cci_llm_error")) {
      status$status <- "error"
      status$error <- as.character(response[[1]])
      cci_write_markdown(c("# CCI LLM analysis", "", "**LLM status:** `error`", "", sprintf("LLM error: `%s`", status$error), "", rule_lines), analysis_path)
    } else {
      status$status <- "ok"
      cci_write_markdown(c("# CCI LLM analysis", "", "**LLM status:** `ok`", "", "## 数据驱动审计", "", rule_lines, "", "## 在线 LLM 深度解读", "", response), analysis_path)
      writeLines(response, raw_path, useBytes = TRUE)
    }
  }
  cci_write_json(status, status_path)
  status$status_json <- cci_path_normalize(status_path)
  status
}

cci_append_manifest_row <- function(row, manifest_path) {
  cci_safe_mkdir(dirname(manifest_path))
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

cci_register_figure <- function(figure_path,
                                data,
                                method,
                                figure_type,
                                title,
                                extra_context = NULL,
                                llm_config = list(),
                                manifest_path = NULL,
                                data_path = NULL,
                                scenario = "CCI") {
  figure_path <- cci_path_normalize(figure_path)
  stem <- cci_stem(figure_path)
  if (is.null(data_path) || !nzchar(cci_to_scalar(data_path))) data_path <- paste0(stem, "_data.csv")
  if (is.null(manifest_path) || !nzchar(cci_to_scalar(manifest_path))) {
    manifest_path <- file.path(dirname(figure_path), "figure_companion_manifest.tsv")
  }
  cci_write_csv(data, data_path)
  llm_status <- cci_run_table_llm_interpretation(
    df = data,
    output_dir = dirname(figure_path),
    prefix = basename(stem),
    title = title,
    scenario = paste(scenario, method, figure_type, sep = " | "),
    extra_context = extra_context,
    llm_config = llm_config
  )
  row <- data.frame(
    method = method,
    figure_type = figure_type,
    title = title,
    figure_path = figure_path,
    data_csv = cci_path_normalize(data_path),
    llm_analysis_md = cci_path_normalize(file.path(dirname(figure_path), sprintf("%s_LLM_ANALYSIS.md", basename(stem)))),
    llm_prompt_md = cci_path_normalize(file.path(dirname(figure_path), sprintf("%s_LLM_PROMPT.md", basename(stem)))),
    llm_status_json = cci_path_normalize(file.path(dirname(figure_path), sprintf("%s_LLM_STATUS.json", basename(stem)))),
    llm_status = llm_status$status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cci_append_manifest_row(row, manifest_path)
  invisible(c(llm_status, list(manifest_path = cci_path_normalize(manifest_path))))
}

cci_bind_rows_flexible <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && is.data.frame(x), rows)
  if (length(rows) == 0L) return(data.frame())
  all_names <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    for (nm in setdiff(all_names, colnames(df))) df[[nm]] <- NA
    df[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

cci_empty_llm_interaction_df <- function() {
  data.frame(
    record_id = character(),
    scenario = character(),
    method = character(),
    evidence_scope = character(),
    statistical_unit = character(),
    sender = character(),
    receiver = character(),
    ligand = character(),
    receptor = character(),
    interaction = character(),
    score_name = character(),
    score_value = numeric(),
    p_value = numeric(),
    rank_value = numeric(),
    support_count = numeric(),
    confidence_label = character(),
    language_guardrail = character(),
    interpretation_hint = character(),
    caveat = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

cci_choose_score_col <- function(df) {
  cci_pick_existing_col(df, c(
    "prob",
    "prioritization_score",
    "prioritization_score_ligand_receptor",
    "lr_prioritization_score",
    "ligand_receptor_score",
    "activity_score",
    "aupr_corrected",
    "detection_fraction",
    "mean_aggregate_rank",
    "aggregate_rank",
    "score",
    "rank"
  ))
}

cci_normalize_interactions_for_llm <- function(df,
                                               scenario,
                                               method,
                                               evidence_scope,
                                               statistical_unit,
                                               confidence_default = "putative",
                                               top_n = 80L) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(cci_empty_llm_interaction_df())
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  source_col <- cci_pick_existing_col(df, c("source", "sender", "sender_celltype", "sender_cell_type", "from"))
  target_col <- cci_pick_existing_col(df, c("target", "receiver", "receiver_celltype", "receiver_cell_type", "to"))
  ligand_col <- cci_pick_existing_col(df, c("ligand", "ligand.complex", "ligand_complex", "ligand_symbol", "test_ligand", "ligand_oi"))
  receptor_col <- cci_pick_existing_col(df, c("receptor", "receptor.complex", "receptor_complex", "receptor_symbol", "receptor_oi"))
  interaction_col <- cci_pick_existing_col(df, c("interaction", "interaction_name", "interaction_string", "lr_pair", "ligand_receptor"))
  score_col <- cci_choose_score_col(df)
  p_col <- cci_pick_existing_col(df, c("pval", "p_val", "p.value", "pvalue", "padj", "p_adj", "qvalue", "q_value", "p_val_adj"))
  rank_col <- cci_pick_existing_col(df, c("aggregate_rank", "mean_aggregate_rank", "rank", "prioritization_rank"))
  support_col <- cci_pick_existing_col(df, c("n_methods_supporting", "n_samples_detected", "support_count", "n_supporting_methods"))
  confidence_col <- cci_pick_existing_col(df, c("confidence_label", "evidence_level", "confidence"))

  score <- if (!is.null(score_col)) cci_num(df[[score_col]]) else rep(NA_real_, nrow(df))
  p_value <- if (!is.null(p_col)) cci_num(df[[p_col]]) else rep(NA_real_, nrow(df))
  rank_value <- if (!is.null(rank_col)) cci_num(df[[rank_col]]) else rep(NA_real_, nrow(df))
  support_count <- if (!is.null(support_col)) cci_num(df[[support_col]]) else rep(NA_real_, nrow(df))
  confidence <- if (!is.null(confidence_col)) cci_to_scalar_list(df[[confidence_col]]) else rep(confidence_default, nrow(df))
  confidence[!nzchar(confidence)] <- confidence_default

  sender <- if (!is.null(source_col)) cci_to_scalar_list(df[[source_col]]) else rep("", nrow(df))
  receiver <- if (!is.null(target_col)) cci_to_scalar_list(df[[target_col]]) else rep("", nrow(df))
  ligand <- if (!is.null(ligand_col)) cci_to_scalar_list(df[[ligand_col]]) else rep("", nrow(df))
  receptor <- if (!is.null(receptor_col)) cci_to_scalar_list(df[[receptor_col]]) else rep("", nrow(df))
  interaction <- if (!is.null(interaction_col)) cci_to_scalar_list(df[[interaction_col]]) else rep("", nrow(df))
  empty_interaction <- !nzchar(interaction)
  interaction[empty_interaction] <- paste(
    ifelse(nzchar(sender[empty_interaction]), sender[empty_interaction], "sender?"),
    "->",
    ifelse(nzchar(receiver[empty_interaction]), receiver[empty_interaction], "receiver?"),
    sprintf("[%s-%s]", ifelse(nzchar(ligand[empty_interaction]), ligand[empty_interaction], "ligand?"), ifelse(nzchar(receptor[empty_interaction]), receptor[empty_interaction], "receptor?"))
  )
  ligand_only <- !nzchar(sender) & !nzchar(receiver) & nzchar(ligand)
  interaction[ligand_only] <- ligand[ligand_only]

  ord <- seq_len(nrow(df))
  if (!is.null(score_col) && any(is.finite(score))) {
    if (grepl("rank", score_col, ignore.case = TRUE)) {
      ord <- order(score, p_value, na.last = TRUE)
    } else {
      ord <- order(-score, p_value, na.last = TRUE)
    }
  } else if (!is.null(p_col) && any(is.finite(p_value))) {
    ord <- order(p_value, na.last = TRUE)
  }
  ord <- utils::head(ord, max(1L, as.integer(top_n)))

  out <- data.frame(
    record_id = sprintf("%s_%03d", cci_sanitize_name(tolower(method)), seq_along(ord)),
    scenario = scenario,
    method = method,
    evidence_scope = evidence_scope,
    statistical_unit = statistical_unit,
    sender = sender[ord],
    receiver = receiver[ord],
    ligand = ligand[ord],
    receptor = receptor[ord],
    interaction = interaction[ord],
    score_name = ifelse(is.null(score_col), "", score_col),
    score_value = score[ord],
    p_value = p_value[ord],
    rank_value = rank_value[ord],
    support_count = support_count[ord],
    confidence_label = confidence[ord],
    language_guardrail = "Use potential/candidate/putative; do not claim direct interaction without spatial/protein validation.",
    interpretation_hint = if (grepl("disease|control|MultiNicheNet", paste(scenario, method), ignore.case = TRUE)) {
      "Prioritize sample-level replicate evidence; cells are not biological replicates."
    } else {
      "Use as scRNA-seq candidate screening / multi-method consensus evidence only."
    },
    caveat = if (grepl("CellChat|CellPhoneDB|LIANA", method, ignore.case = TRUE)) {
      "Exploratory LR inference from dissociated scRNA-seq; validate with orthogonal assays where possible."
    } else if (grepl("NicheNet", method, ignore.case = TRUE)) {
      "Receiver downstream response is inferred from expression and prior networks; not protein-level proof."
    } else {
      "Interpret as computationally inferred CCI evidence."
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  out
}

cci_build_llm_interaction_table <- function(tables,
                                            scenario,
                                            default_statistical_unit,
                                            top_n_per_method = 80L) {
  if (is.null(tables) || length(tables) == 0L) return(cci_empty_llm_interaction_df())
  rows <- lapply(names(tables), function(nm) {
    spec <- tables[[nm]]
    if (is.null(spec)) return(NULL)
    if (is.data.frame(spec)) spec <- list(data = spec)
    if (!is.list(spec) || is.null(spec$data)) return(NULL)
    cci_normalize_interactions_for_llm(
      df = spec$data,
      scenario = scenario,
      method = cci_to_scalar(cci_null_coalesce(spec$method, nm)),
      evidence_scope = cci_to_scalar(cci_null_coalesce(spec$evidence_scope, nm)),
      statistical_unit = cci_to_scalar(cci_null_coalesce(spec$statistical_unit, default_statistical_unit)),
      confidence_default = cci_to_scalar(cci_null_coalesce(spec$confidence_default, "putative")),
      top_n = cci_null_coalesce(spec$top_n, top_n_per_method)
    )
  })
  out <- cci_bind_rows_flexible(rows)
  if (nrow(out) == 0L) return(cci_empty_llm_interaction_df())
  out$record_id <- sprintf("cci_%03d", seq_len(nrow(out)))
  out
}

cci_plot_bar <- function(plot_df,
                         label_col,
                         value_col,
                         output_prefix,
                         title,
                         xlab = "Evidence metric",
                         top_n = 20L,
                         fill = "#4E79A7") {
  if (is.null(plot_df) || !is.data.frame(plot_df) || nrow(plot_df) == 0L) return(NULL)
  plot_df <- as.data.frame(plot_df, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c(label_col, value_col) %in% colnames(plot_df))) return(NULL)
  vals <- cci_num(plot_df[[value_col]])
  keep <- is.finite(vals)
  plot_df <- plot_df[keep, , drop = FALSE]
  vals <- vals[keep]
  if (nrow(plot_df) == 0L) return(NULL)
  ord <- order(vals, decreasing = TRUE, na.last = TRUE)
  plot_df <- plot_df[utils::head(ord, max(1L, as.integer(top_n))), , drop = FALSE]
  vals <- cci_num(plot_df[[value_col]])
  labels <- cci_to_scalar_list(plot_df[[label_col]])
  labels[!nzchar(labels)] <- "NA"
  labels <- ifelse(nchar(labels) > 65L, paste0(substr(labels, 1L, 62L), "..."), labels)

  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)
  draw_once <- function(device_fun) {
    device_fun()
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mar = c(5, 14, 4, 2) + 0.1)
    vals_rev <- rev(vals)
    labels_rev <- rev(labels)
    bp <- barplot(
      vals_rev,
      names.arg = labels_rev,
      horiz = TRUE,
      las = 1,
      cex.names = 0.72,
      col = fill,
      border = NA,
      main = title,
      xlab = xlab,
      xlim = c(0, max(vals_rev, na.rm = TRUE) * 1.18 + 1e-9)
    )
    text(vals_rev, bp, labels = signif(vals_rev, 3), pos = 4, cex = 0.7)
    grDevices::dev.off()
  }
  cci_safe_mkdir(dirname(pdf_path))
  draw_once(function() grDevices::pdf(pdf_path, width = 12, height = 8))
  draw_once(function() grDevices::png(png_path, width = 2160, height = 1440, res = 180))
  invisible(list(pdf = cci_path_normalize(pdf_path), png = cci_path_normalize(png_path), data = plot_df))
}

cci_plot_ranked_interactions <- function(df,
                                         output_prefix,
                                         title,
                                         scenario,
                                         method,
                                         evidence_scope,
                                         statistical_unit,
                                         top_n = 20L,
                                         fill = "#4E79A7") {
  norm <- cci_normalize_interactions_for_llm(
    df = df,
    scenario = scenario,
    method = method,
    evidence_scope = evidence_scope,
    statistical_unit = statistical_unit,
    top_n = top_n
  )
  if (nrow(norm) == 0L || !any(is.finite(norm$score_value))) return(NULL)
  norm$plot_value <- norm$score_value
  rank_like <- grepl("rank", norm$score_name, ignore.case = TRUE)
  norm$plot_value[rank_like] <- -log10(pmax(norm$score_value[rank_like], .Machine$double.xmin))
  norm$plot_label <- norm$interaction
  metric_label <- if (any(rank_like)) "Evidence metric (rank transformed as -log10(rank); higher is stronger)" else unique(norm$score_name)[1]
  cci_plot_bar(
    norm,
    label_col = "plot_label",
    value_col = "plot_value",
    output_prefix = output_prefix,
    title = title,
    xlab = metric_label,
    top_n = top_n,
    fill = fill
  )
}

cci_summarise_sender_receiver <- function(df,
                                          scenario,
                                          method,
                                          evidence_scope,
                                          statistical_unit,
                                          top_n = 25L) {
  norm <- cci_normalize_interactions_for_llm(
    df = df,
    scenario = scenario,
    method = method,
    evidence_scope = evidence_scope,
    statistical_unit = statistical_unit,
    top_n = max(top_n * 5L, 100L)
  )
  if (nrow(norm) == 0L || !all(c("sender", "receiver") %in% colnames(norm))) return(data.frame())
  keep <- nzchar(norm$sender) & nzchar(norm$receiver)
  norm <- norm[keep, , drop = FALSE]
  if (nrow(norm) == 0L) return(data.frame())
  norm$sender_receiver <- paste(norm$sender, "->", norm$receiver)
  split_df <- split(norm, norm$sender_receiver)
  rows <- lapply(split_df, function(x) {
    score <- x$score_value
    rank_like <- grepl("rank", x$score_name, ignore.case = TRUE)
    score[rank_like] <- -log10(pmax(score[rank_like], .Machine$double.xmin))
    data.frame(
      sender_receiver = x$sender_receiver[1],
      sender = x$sender[1],
      receiver = x$receiver[1],
      n_interactions = nrow(x),
      mean_evidence_metric = mean(score, na.rm = TRUE),
      top_interaction = x$interaction[which.max(ifelse(is.finite(score), score, -Inf))][1],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- cci_bind_rows_flexible(rows)
  if (nrow(out) == 0L) return(out)
  out <- out[order(-out$n_interactions, -out$mean_evidence_metric), , drop = FALSE]
  utils::head(out, top_n)
}

cci_plot_sender_receiver_summary <- function(df,
                                             output_prefix,
                                             title,
                                             scenario,
                                             method,
                                             evidence_scope,
                                             statistical_unit,
                                             top_n = 20L,
                                             fill = "#59A14F") {
  pair_df <- cci_summarise_sender_receiver(
    df = df,
    scenario = scenario,
    method = method,
    evidence_scope = evidence_scope,
    statistical_unit = statistical_unit,
    top_n = top_n
  )
  if (nrow(pair_df) == 0L) return(NULL)
  cci_plot_bar(
    pair_df,
    label_col = "sender_receiver",
    value_col = "n_interactions",
    output_prefix = output_prefix,
    title = title,
    xlab = "Number of candidate LR records",
    top_n = top_n,
    fill = fill
  )
}

cci_replicate_check_to_df <- function(replicate_check) {
  if (is.null(replicate_check) || is.null(replicate_check$sample_table)) return(data.frame())
  tab <- replicate_check$sample_table
  df <- as.data.frame(as.table(tab), stringsAsFactors = FALSE)
  colnames(df) <- c("celltype", "group", "n_samples")
  inadequate <- rep(FALSE, nrow(df))
  if (!is.null(replicate_check$inadequate_pairs) && nrow(replicate_check$inadequate_pairs) > 0L) {
    bad <- apply(replicate_check$inadequate_pairs, 1, function(idx) {
      paste(rownames(tab)[idx[1]], colnames(tab)[idx[2]], sep = "||")
    })
    inadequate <- paste(df$celltype, df$group, sep = "||") %in% bad
  }
  df$replicate_status <- ifelse(inadequate, "below_min_samples", "ok")
  df
}

cci_extract_multinichenet_prioritized_df <- function(mn_result) {
  if (is.null(mn_result)) return(data.frame())
  if (is.data.frame(mn_result)) return(mn_result)
  if (is.list(mn_result) && is.data.frame(mn_result$group_prioritization_tbl)) return(mn_result$group_prioritization_tbl)
  if (is.list(mn_result) && is.list(mn_result$prioritized_tbl) && is.data.frame(mn_result$prioritized_tbl$group_prioritization_tbl)) {
    return(mn_result$prioritized_tbl$group_prioritization_tbl)
  }
  data.frame()
}

cci_extract_cellchat_net_summary <- function(cellchat_viz) {
  cc_obj <- NULL
  if (is.null(cellchat_viz)) return(data.frame())
  if (is.list(cellchat_viz) && !is.null(cellchat_viz$cellchat_net_summary) && is.data.frame(cellchat_viz$cellchat_net_summary)) {
    return(cellchat_viz$cellchat_net_summary)
  }
  if (is.list(cellchat_viz) && !is.null(cellchat_viz$cc_merged)) cc_obj <- cellchat_viz$cc_merged else cc_obj <- cellchat_viz
  net <- tryCatch(cc_obj@net, error = function(e) NULL)
  if (is.null(net) && is.list(cc_obj) && !is.null(cc_obj$net)) net <- cc_obj$net
  if (is.null(net) || (!is.matrix(net$count) && !is.matrix(net$weight))) return(data.frame())
  count <- if (is.matrix(net$count)) net$count else matrix(NA_real_, nrow = nrow(net$weight), ncol = ncol(net$weight), dimnames = dimnames(net$weight))
  weight <- if (is.matrix(net$weight)) net$weight else matrix(NA_real_, nrow = nrow(net$count), ncol = ncol(net$count), dimnames = dimnames(net$count))
  src <- rownames(count)
  tgt <- colnames(count)
  if (is.null(src)) src <- sprintf("source_%s", seq_len(nrow(count)))
  if (is.null(tgt)) tgt <- sprintf("target_%s", seq_len(ncol(count)))
  grid <- expand.grid(source = src, target = tgt, stringsAsFactors = FALSE)
  grid$count <- as.numeric(count[cbind(match(grid$source, src), match(grid$target, tgt))])
  grid$weight <- as.numeric(weight[cbind(match(grid$source, src), match(grid$target, tgt))])
  grid$sender_receiver <- paste(grid$source, "->", grid$target)
  grid <- grid[(is.finite(grid$count) & grid$count > 0) | (is.finite(grid$weight) & grid$weight > 0), , drop = FALSE]
  grid[order(-grid$count, -grid$weight), , drop = FALSE]
}

cci_plot_scrna_visualizations <- function(cellchat_candidates,
                                          liana_consensus,
                                          nichenet_result,
                                          output_dir,
                                          prefix = "cci_scrna",
                                          llm_config = list()) {
  fig_dir <- file.path(output_dir, "figures")
  cci_safe_mkdir(fig_dir)
  manifest_path <- file.path(fig_dir, paste0(prefix, "_figure_companion_manifest.tsv"))
  if (file.exists(manifest_path)) unlink(manifest_path)
  registered <- list()

  register_if <- function(paths, method, figure_type, title, context) {
    if (is.null(paths) || is.null(paths$png) || is.null(paths$data)) return(invisible(NULL))
    status <- cci_register_figure(
      figure_path = paths$png,
      data = paths$data,
      method = method,
      figure_type = figure_type,
      title = title,
      extra_context = context,
      llm_config = llm_config,
      manifest_path = manifest_path,
      scenario = "A_scrna"
    )
    registered[[length(registered) + 1L]] <<- list(paths = paths, llm = status)
    invisible(NULL)
  }

  cc_rank <- cci_plot_ranked_interactions(
    cellchat_candidates,
    output_prefix = file.path(fig_dir, paste0(prefix, "_cellchat_top_candidates")),
    title = "CellChat top putative/candidate LR interactions",
    scenario = "A_scrna",
    method = "CellChat candidates",
    evidence_scope = "candidate LR screening",
    statistical_unit = "cells for screening only",
    fill = "#4E79A7"
  )
  register_if(cc_rank, "CellChat", "top_lr_barplot", "CellChat top putative/candidate LR interactions", "Use as candidate screening only; no spatial/protein validation is implied.")

  cc_pairs <- cci_plot_sender_receiver_summary(
    cellchat_candidates,
    output_prefix = file.path(fig_dir, paste0(prefix, "_cellchat_sender_receiver_summary")),
    title = "CellChat sender-receiver candidate counts",
    scenario = "A_scrna",
    method = "CellChat candidates",
    evidence_scope = "sender-receiver candidate count",
    statistical_unit = "cells for screening only",
    fill = "#59A14F"
  )
  register_if(cc_pairs, "CellChat", "sender_receiver_barplot", "CellChat sender-receiver candidate counts", "Counts summarize inferred candidate LR records, not validated physical contacts.")

  li_rank <- cci_plot_ranked_interactions(
    liana_consensus,
    output_prefix = file.path(fig_dir, paste0(prefix, "_liana_consensus_top")),
    title = "LIANA+ multi-method consensus top interactions",
    scenario = "A_scrna",
    method = "LIANA+ consensus",
    evidence_scope = "multi-method consensus ranking",
    statistical_unit = "cells for screening only",
    fill = "#F28E2B"
  )
  register_if(li_rank, "LIANA+", "consensus_rank_barplot", "LIANA+ multi-method consensus top interactions", "Lower LIANA aggregate ranks are transformed to stronger visual evidence; interpret conservatively.")

  nn_df <- if (!is.null(nichenet_result) && is.data.frame(nichenet_result$ligand_activity)) nichenet_result$ligand_activity else data.frame()
  nn_rank <- cci_plot_ranked_interactions(
    nn_df,
    output_prefix = file.path(fig_dir, paste0(prefix, "_nichenet_top_ligands")),
    title = "NicheNet top receiver-response ligands",
    scenario = "A_scrna",
    method = "NicheNet ligand activity",
    evidence_scope = "receiver downstream response",
    statistical_unit = "cells for receiver DEG / geneset inference",
    fill = "#E15759"
  )
  register_if(nn_rank, "NicheNet", "ligand_activity_barplot", "NicheNet top receiver-response ligands", "NicheNet supports potential upstream ligands explaining receiver response, not direct signaling proof.")

  list(
    figure_dir = cci_path_normalize(fig_dir),
    manifest_path = cci_path_normalize(manifest_path),
    registered = registered
  )
}

cci_plot_disease_control_visualizations <- function(mn_result,
                                                    liana_per_sample,
                                                    cellchat_viz,
                                                    output_dir,
                                                    prefix = "cci_disease_ctrl",
                                                    llm_config = list()) {
  fig_dir <- file.path(output_dir, "figures")
  cci_safe_mkdir(fig_dir)
  manifest_path <- file.path(fig_dir, paste0(prefix, "_figure_companion_manifest.tsv"))
  if (file.exists(manifest_path)) unlink(manifest_path)
  registered <- list()

  register_if <- function(paths, method, figure_type, title, context) {
    if (is.null(paths) || is.null(paths$png) || is.null(paths$data)) return(invisible(NULL))
    status <- cci_register_figure(
      figure_path = paths$png,
      data = paths$data,
      method = method,
      figure_type = figure_type,
      title = title,
      extra_context = context,
      llm_config = llm_config,
      manifest_path = manifest_path,
      scenario = "B_disease_control"
    )
    registered[[length(registered) + 1L]] <<- list(paths = paths, llm = status)
    invisible(NULL)
  }

  mn_df <- cci_extract_multinichenet_prioritized_df(mn_result)
  mn_rank <- cci_plot_ranked_interactions(
    mn_df,
    output_prefix = file.path(fig_dir, paste0(prefix, "_multinichenet_top_prioritized")),
    title = "MultiNicheNet top prioritized disease-control CCI candidates",
    scenario = "B_disease_control",
    method = "MultiNicheNet prioritized",
    evidence_scope = "sample-level DE-based CCI prioritization",
    statistical_unit = "sample / donor",
    fill = "#4E79A7"
  )
  register_if(mn_rank, "MultiNicheNet", "prioritized_lr_barplot", "MultiNicheNet top prioritized disease-control CCI candidates", "Primary evidence should come from sample-level replicated DE; cells are not replicates.")

  li_df <- if (!is.null(liana_per_sample) && is.data.frame(liana_per_sample$summary)) liana_per_sample$summary else data.frame()
  li_rank <- cci_plot_ranked_interactions(
    li_df,
    output_prefix = file.path(fig_dir, paste0(prefix, "_liana_per_sample_support")),
    title = "LIANA+ per-sample support for candidate LR interactions",
    scenario = "B_disease_control",
    method = "LIANA+ per-sample",
    evidence_scope = "auxiliary per-sample aggregation",
    statistical_unit = "sample / donor",
    fill = "#F28E2B"
  )
  register_if(li_rank, "LIANA+", "per_sample_support_barplot", "LIANA+ per-sample support for candidate LR interactions", "Auxiliary evidence: prioritize interactions recurring across biological samples.")

  rep_df <- data.frame()
  if (is.list(mn_result) && !is.null(mn_result$replicate_check)) rep_df <- cci_replicate_check_to_df(mn_result$replicate_check)
  if (nrow(rep_df) > 0L) {
    rep_df$celltype_group <- paste(rep_df$celltype, rep_df$group, sep = " | ")
    rep_paths <- cci_plot_bar(
      rep_df,
      label_col = "celltype_group",
      value_col = "n_samples",
      output_prefix = file.path(fig_dir, paste0(prefix, "_replicate_adequacy")),
      title = "Replicate adequacy by cell type and group",
      xlab = "Number of biological samples",
      fill = "#59A14F"
    )
    register_if(rep_paths, "Replicate QC", "sample_count_barplot", "Replicate adequacy by cell type and group", "This plot checks the core rule: samples/donors are replicates, not cells.")
  }

  cc_net <- cci_extract_cellchat_net_summary(cellchat_viz)
  if (nrow(cc_net) > 0L) {
    value_col <- if ("weight" %in% colnames(cc_net) && any(is.finite(cci_num(cc_net$weight)))) "weight" else "count"
    cc_paths <- cci_plot_bar(
      cc_net,
      label_col = "sender_receiver",
      value_col = value_col,
      output_prefix = file.path(fig_dir, paste0(prefix, "_cellchat_global_network_summary")),
      title = "CellChat global network summary (visualization only)",
      xlab = paste("CellChat", value_col),
      fill = "#E15759"
    )
    register_if(cc_paths, "CellChat", "global_network_barplot", "CellChat global network summary (visualization only)", "Exploratory network display only; statistical claims should use MultiNicheNet sample-level evidence.")
  }

  list(
    figure_dir = cci_path_normalize(fig_dir),
    manifest_path = cci_path_normalize(manifest_path),
    registered = registered
  )
}

cci_pkg_check <- function(pkgs, scenario = "") {
  # Check required packages; stop with informative message on missing
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "[cci] Missing packages for %s: %s\nInstall via BiocManager::install() or remotes::install_github().",
      scenario,
      paste(missing_pkgs, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

# Validate that a Seurat object or data.frame metadata has required columns
cci_check_required_cols <- function(meta, required, context = "") {
  missing_cols <- setdiff(required, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(sprintf("[cci%s] Missing required metadata columns: %s",
                 if (nzchar(context)) paste0(":", context) else "",
                 paste(missing_cols, collapse = ", ")))
  }
  invisible(TRUE)
}

# Compute per-(sender, receiver) cell counts; needed before any CCI run
cci_count_cells_per_pair <- function(meta, celltype_col, min_cells = 10L) {
  ct_counts <- table(meta[[celltype_col]])
  small <- names(ct_counts)[ct_counts < min_cells]
  if (length(small) > 0) {
    message(sprintf(
      "[cci] Warning: %d cell type(s) below min_cells=%d and will be flagged: %s",
      length(small), min_cells, paste(small, collapse = ", ")
    ))
  }
  list(counts = ct_counts, below_threshold = small)
}

# Validate sample-level adequacy for Scenario B
# Ensures >= min_samples donors per group before any DE-based CCI
cci_validate_replicates <- function(meta,
                                    sample_col    = "sample",
                                    group_col     = "disease_status",
                                    celltype_col  = "cell_type_L2",
                                    min_samples   = 3L) {
  cci_check_required_cols(meta, c(sample_col, group_col, celltype_col))

  # One row per (celltype, group, sample) combination
  per_cg <- unique(meta[, c(celltype_col, group_col, sample_col)])
  sample_tab <- table(per_cg[[celltype_col]], per_cg[[group_col]])
  
  # Flag (celltype, group) combinations with insufficient replicates
  inadequate <- which(sample_tab > 0 & sample_tab < min_samples, arr.ind = TRUE)
  
  if (nrow(inadequate) > 0) {
    msgs <- apply(inadequate, 1, function(idx) {
      sprintf("  celltype='%s' group='%s'  n_samples=%d",
              rownames(sample_tab)[idx[1]],
              colnames(sample_tab)[idx[2]],
              sample_tab[idx[1], idx[2]])
    })
    message(sprintf(
      "[cci] Warning: %d (celltype x group) pairs have < %d samples.\n%s",
      nrow(inadequate), min_samples, paste(msgs, collapse = "\n")
    ))
  } else {
    message(sprintf("[cci] Replicate check passed: all (celltype x group) pairs >= %d samples.", min_samples))
  }

  list(sample_table = sample_tab, inadequate_pairs = inadequate)
}

# ---------------------------------------------------------------------------- #
#  Section 1  Scenario A  CellChat screening                                   #
# ---------------------------------------------------------------------------- #

# Prepare a CellChat object from a Seurat object (v5-compatible)
# Returns a CellChat object ready for computeCommunProb()
cci_prepare_cellchat <- function(seurat_obj,
                                 celltype_col    = "cell_type_L2",
                                 species         = c("human", "mouse"),
                                 assay           = "RNA",
                                 slot_use        = "data",
                                 min_cells       = 10L) {
  cci_pkg_check(c("CellChat", "SeuratObject"), scenario = "CellChat screening")
  species <- match.arg(species)

  # Load species-specific ligand-receptor database
  db <- if (species == "human") CellChat::CellChatDB.human else CellChat::CellChatDB.mouse
  message(sprintf("[cci] CellChat database: %s  (%d interactions)", species, nrow(db$interaction)))

  # Extract normalized expression matrix
  expr_mat <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot_use),
    error = function(e) stop(sprintf("[cci] Cannot extract expression from Seurat: %s", e$message))
  )

  meta <- seurat_obj@meta.data
  cci_check_required_cols(meta, celltype_col)
  cell_counts <- cci_count_cells_per_pair(meta, celltype_col, min_cells)

  # Remove cell types below threshold to prevent unstable estimates
  keep_cells <- !(meta[[celltype_col]] %in% cell_counts$below_threshold)
  if (sum(!keep_cells) > 0) {
    message(sprintf("[cci] Removing %d cells from %d low-count cell types.",
                    sum(!keep_cells), length(cell_counts$below_threshold)))
    expr_mat <- expr_mat[, keep_cells, drop = FALSE]
    meta <- meta[keep_cells, , drop = FALSE]
  }

  cellchat_obj <- tryCatch({
    cc <- CellChat::createCellChat(
      object   = expr_mat,
      meta     = meta,
      group.by = celltype_col
    )
    cc@DB <- db
    cc
  }, error = function(e) {
    stop(sprintf("[cci] CellChat object creation failed: %s", e$message))
  })

  message("[cci] CellChat object created. Run cci_run_cellchat_screening() next.")
  cellchat_obj
}

# Run the core CellChat inference pipeline
# Returns cellchat object with computed probabilities + netP
cci_run_cellchat_screening <- function(cellchat_obj,
                                       population_size = TRUE,
                                       nboot           = 100L,
                                       workers         = 4L) {
  cci_pkg_check(c("CellChat", "future"), scenario = "CellChat screening")

  message("[cci] Subsetting CellChat database to expressed interactions...")
  cellchat_obj <- CellChat::subsetData(cellchat_obj)

  message(sprintf("[cci] Computing communication probabilities (nboot=%d, workers=%d)...",
                  nboot, workers))
  future::plan("multisession", workers = workers)
  on.exit(future::plan("sequential"), add = TRUE)
  cellchat_obj <- tryCatch(
    CellChat::computeCommunProb(
      cellchat_obj,
      type              = "triMean",
      population.size   = population_size,
      nboot             = nboot
    ),
    error = function(e) stop(sprintf("[cci] computeCommunProb failed: %s", e$message))
  )
  # Filter: require >= 10 cells in both sender and receiver
  cellchat_obj <- CellChat::filterCommunication(cellchat_obj, min.cells = 10)

  message("[cci] Computing pathway-level communication...")
  cellchat_obj <- CellChat::computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- CellChat::aggregateNet(cellchat_obj)

  message("[cci] CellChat screening complete.")
  cellchat_obj
}

# Extract a tidy candidate LR table from CellChat with conservative labeling
# All interactions labeled as "putative" or "candidate" per language convention
cci_extract_cellchat_candidates <- function(cellchat_obj,
                                            prob_thresh    = 0.0,
                                            pval_thresh    = 0.05,
                                            top_n          = NULL) {
  lr_df <- tryCatch(
    CellChat::subsetCommunication(cellchat_obj),
    error = function(e) {
      stop(sprintf("[cci] Failed to extract LR table from CellChat: %s", e$message))
    }
  )

  # Apply thresholds
  lr_df <- lr_df[lr_df$prob   >= prob_thresh, , drop = FALSE]
  lr_df <- lr_df[lr_df$pval   <= pval_thresh, , drop = FALSE]

  # Sort by probability descending
  lr_df <- lr_df[order(lr_df$prob, decreasing = TRUE), , drop = FALSE]

  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
    lr_df <- head(lr_df, as.integer(top_n))
  }

  # Language rule: annotate confidence level conservatively
  # Never write "confirmed" or "direct"; use "putative" / "candidate"
  lr_df$confidence_label <- character(nrow(lr_df))
  if (nrow(lr_df) > 0L) {
    lr_df$confidence_label <- "putative"
    prob_q75 <- suppressWarnings(stats::quantile(lr_df$prob, 0.75, na.rm = TRUE))
    if (length(prob_q75) == 1L && is.finite(prob_q75)) {
      lr_df$confidence_label[lr_df$pval <= 0.05 & lr_df$prob >= prob_q75] <- "candidate"
    }
  }

  # Add human-readable interaction string
  lr_df$interaction_string <- paste0(
    lr_df$source, " -> ", lr_df$target,
    " [", lr_df$interaction_name, "] (", lr_df$confidence_label, ")"
  )

  message(sprintf("[cci] Extracted %d putative/candidate LR interactions from CellChat.", nrow(lr_df)))
  lr_df
}

# ---------------------------------------------------------------------------- #
#  Section 2  Scenario A  LIANA+ multi-method consensus                        #
# ---------------------------------------------------------------------------- #

# Run LIANA+ with multiple methods and return consensus ranking
# Input: Seurat object or SingleCellExperiment
# Returns: tidy data.frame with aggregate_rank and method vote counts
cci_run_liana_consensus <- function(seurat_obj,
                                    celltype_col = "cell_type_L2",
                                    methods      = c("natmi", "connectome", "logfc",
                                                     "sca", "cellphonedb"),
                                    species      = c("human", "mouse"),
                                    workers      = 4L,
                                    verbose      = TRUE) {
  cci_pkg_check(c("liana", "OmnipathR"), scenario = "LIANA+ consensus")
  species <- match.arg(species)

  # Set identity for LIANA
  if (inherits(seurat_obj, "Seurat")) {
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[celltype_col]]
  }

  resource_name <- if (species == "human") "Consensus" else "mouseCons"
  message(sprintf("[cci] LIANA+ resource: %s | methods: %s",
                  resource_name, paste(methods, collapse = ", ")))

  liana_result <- tryCatch(
    liana::liana_wrap(
      seurat_obj,
      method      = methods,
      resource    = resource_name,
      verbose     = verbose,
      parallelize = workers > 1L,
      workers     = workers
    ),
    error = function(e) stop(sprintf("[cci] LIANA+ wrap failed: %s", e$message))
  )

  # Aggregate to consensus ranking (lower rank = stronger evidence)
  consensus_df <- tryCatch(
    liana::liana_aggregate(liana_result),
    error = function(e) stop(sprintf("[cci] LIANA+ aggregate failed: %s", e$message))
  )

  # Count how many methods support each interaction
  method_cols <- intersect(
    paste0(methods, ".rank"),
    colnames(consensus_df)
  )
  if (length(method_cols) > 0) {
    consensus_df$n_methods_supporting <- rowSums(
      !is.na(consensus_df[, method_cols, drop = FALSE])
    )
  } else {
    consensus_df$n_methods_supporting <- NA_integer_
  }

  # Language rule: confidence label based on method consensus
  consensus_df$confidence_label <- "putative"
  consensus_df$confidence_label[
    !is.na(consensus_df$n_methods_supporting) &
    consensus_df$n_methods_supporting >= ceiling(length(methods) * 0.6)
  ] <- "candidate"

  message(sprintf(
    "[cci] LIANA+ consensus: %d interactions; %d candidate (>= %.0f%% method agreement).",
    nrow(consensus_df),
    sum(consensus_df$confidence_label == "candidate", na.rm = TRUE),
    ceiling(length(methods) * 0.6) / length(methods) * 100
  ))

  consensus_df
}

# Cross-validate CellChat candidates against LIANA+ consensus
# Returns only interactions supported by both approaches
cci_intersect_cellchat_liana <- function(cellchat_candidates,
                                          liana_consensus,
                                          rank_thresh = 0.05) {
  # LIANA consensus rank in [0,1]; lower = better supported
  liana_top <- liana_consensus[
    !is.na(liana_consensus$aggregate_rank) &
    liana_consensus$aggregate_rank <= rank_thresh, , drop = FALSE
  ]

  # Build a join key: source|target|ligand.complex|receptor.complex
  cellchat_key <- paste(
    cellchat_candidates$source,
    cellchat_candidates$target,
    cellchat_candidates$ligand,
    cellchat_candidates$receptor,
    sep = "|"
  )
  liana_key <- paste(
    liana_top$source,
    liana_top$target,
    liana_top$ligand.complex,
    liana_top$receptor.complex,
    sep = "|"
  )

  shared_idx <- cellchat_candidates[cellchat_key %in% liana_key, , drop = FALSE]
  shared_idx$evidence_level <- "multi-method consensus (putative)"

  message(sprintf(
    "[cci] Intersection: %d of %d CellChat candidates supported by LIANA+ (top %.0f%% rank).",
    nrow(shared_idx), nrow(cellchat_candidates), rank_thresh * 100
  ))
  shared_idx
}

# ---------------------------------------------------------------------------- #
#  Section 3  Scenario A  NicheNet receiver downstream analysis                #
# ---------------------------------------------------------------------------- #

# Prepare NicheNet inputs: define sender/receiver, potential ligands, geneset
cci_prepare_nichenet_inputs <- function(seurat_obj,
                                        receiver_celltype,
                                        sender_celltypes    = "all",
                                        celltype_col        = "cell_type_L2",
                                        condition_col       = "disease_status",
                                        condition_test      = NULL,
                                        condition_reference = NULL,
                                        top_n_ligands       = 30L,
                                        geneset_method      = c("DE", "custom"),
                                        custom_geneset      = NULL) {
  cci_pkg_check(c("nichenetr", "Matrix", "SeuratObject"), scenario = "NicheNet")
  geneset_method <- match.arg(geneset_method)

  meta <- seurat_obj@meta.data
  required_cols <- celltype_col
  if (identical(geneset_method, "DE")) {
    if (is.null(condition_col) || !nzchar(cci_trim(condition_col))) {
      stop("[cci:NicheNet] condition_col is required when geneset_method='DE'.")
    }
    required_cols <- c(required_cols, condition_col)
  }
  cci_check_required_cols(meta, required_cols)

  # Receiver cells
  receiver_cells <- rownames(meta)[meta[[celltype_col]] == receiver_celltype]
  if (length(receiver_cells) == 0) {
    stop(sprintf("[cci:NicheNet] receiver_celltype='%s' not found in %s.", receiver_celltype, celltype_col))
  }

  # Sender cells
  if (identical(sender_celltypes, "all")) {
    sender_celltypes <- setdiff(unique(meta[[celltype_col]]), receiver_celltype)
  }

  # Background gene set: all expressed genes in receiver
  receiver_expr <- SeuratObject::GetAssayData(
    seurat_obj[, receiver_cells],
    slot = "data"
  )
  background_genes <- rownames(receiver_expr)[
    Matrix::rowMeans(receiver_expr > 0) > 0.1
  ]

  # Geneset of interest: DEGs in receiver between conditions
  geneset_oi <- if (geneset_method == "DE" && !is.null(condition_test)) {
    message(sprintf("[cci:NicheNet] Computing DE genes in receiver '%s': %s vs %s",
                    receiver_celltype, condition_test, condition_reference))
    sub_obj <- seurat_obj[, receiver_cells]
    Seurat::Idents(sub_obj) <- sub_obj@meta.data[[condition_col]]
    de_markers <- tryCatch(
      Seurat::FindMarkers(
        sub_obj,
        ident.1   = condition_test,
        ident.2   = condition_reference,
        min.pct   = 0.1,
        logfc.threshold = 0.25
      ),
      error = function(e) stop(sprintf("[cci:NicheNet] FindMarkers failed: %s", e$message))
    )
    rownames(de_markers)[
      de_markers$p_val_adj < 0.05 & abs(de_markers$avg_log2FC) > 0.25
    ]
  } else if (geneset_method == "custom" && !is.null(custom_geneset)) {
    intersect(custom_geneset, background_genes)
  } else {
    stop("[cci:NicheNet] Provide condition_test/reference (geneset_method='DE') or custom_geneset.")
  }

  message(sprintf("[cci:NicheNet] geneset_oi: %d genes | background: %d genes",
                  length(geneset_oi), length(background_genes)))

  list(
    receiver_celltype  = receiver_celltype,
    sender_celltypes   = sender_celltypes,
    celltype_col       = celltype_col,
    background_genes   = background_genes,
    geneset_oi         = geneset_oi,
    receiver_cells     = receiver_cells,
    top_n_ligands      = top_n_ligands
  )
}

# Run NicheNet ligand activity + ligand-target matrix extraction
# Returns: list(ligand_activity, ligand_target_matrix, top_ligands)
cci_run_nichenet <- function(seurat_obj,
                              nichenet_inputs,
                              nichenet_networks_path = NULL,
                              species               = c("human", "mouse")) {
  cci_pkg_check(c("nichenetr", "Matrix", "SeuratObject"), scenario = "NicheNet")
  species <- match.arg(species)

  # Load NicheNet networks
  if (!is.null(nichenet_networks_path) && dir.exists(nichenet_networks_path)) {
    lr_network    <- readRDS(file.path(nichenet_networks_path, "lr_network.rds"))
    weighted_nets <- readRDS(file.path(nichenet_networks_path, "weighted_networks.rds"))
    ligand_target_matrix <- readRDS(file.path(nichenet_networks_path, "ligand_target_matrix.rds"))
  } else {
    message("[cci:NicheNet] Loading networks from nichenetr package (slow on first run)...")
    lr_network    <- nichenetr::lr_network
    weighted_nets <- nichenetr::weighted_networks
    ligand_target_matrix <- nichenetr::ligand_target_matrix
    if (species == "mouse") {
      lr_network$from <- nichenetr::convert_human_to_mouse_symbols(lr_network$from)
      lr_network$to   <- nichenetr::convert_human_to_mouse_symbols(lr_network$to)
      lr_network <- lr_network[!is.na(lr_network$from) & !is.na(lr_network$to), , drop = FALSE]
    }
  }

  # Collect expressed ligands in sender cells
  meta <- seurat_obj@meta.data
  celltype_col <- cci_null_coalesce(nichenet_inputs$celltype_col, "cell_type_L2")
  cci_check_required_cols(meta, celltype_col, context = "NicheNet sender selection")
  sender_cells <- rownames(meta)[
    as.character(meta[[celltype_col]]) %in% as.character(nichenet_inputs$sender_celltypes)
  ]
  if (length(sender_cells) == 0L) {
    stop(sprintf(
      "[cci:NicheNet] No sender cells found in '%s' for sender_celltypes: %s",
      celltype_col,
      paste(nichenet_inputs$sender_celltypes, collapse = ", ")
    ))
  }
  all_ligands_in_network <- unique(lr_network$from)

  sender_expr <- SeuratObject::GetAssayData(seurat_obj[, sender_cells], slot = "data")
  expressed_ligands <- all_ligands_in_network[
    all_ligands_in_network %in% rownames(sender_expr)
  ]
  expressed_ligands <- expressed_ligands[
    Matrix::rowMeans(sender_expr[expressed_ligands, , drop = FALSE] > 0) > 0.1
  ]
  if (length(expressed_ligands) == 0L) {
    stop("[cci:NicheNet] No expressed ligands found in selected sender cells.")
  }

  message(sprintf("[cci:NicheNet] Expressed ligands from senders: %d", length(expressed_ligands)))

  # Ligand activity analysis
  ligand_activity <- tryCatch(
    nichenetr::predict_ligand_activities(
      geneset               = nichenet_inputs$geneset_oi,
      background_expressed_genes = nichenet_inputs$background_genes,
      ligand_target_matrix  = ligand_target_matrix,
      potential_ligands     = expressed_ligands
    ),
    error = function(e) stop(sprintf("[cci:NicheNet] predict_ligand_activities failed: %s", e$message))
  )

  ligand_activity <- ligand_activity[order(ligand_activity$aupr_corrected, decreasing = TRUE), ]
  top_ligands <- head(ligand_activity$test_ligand, nichenet_inputs$top_n_ligands)

  # Extract ligand-target matrix for top ligands
  lt_subset <- ligand_target_matrix[
    top_ligands[top_ligands %in% rownames(ligand_target_matrix)],
    nichenet_inputs$geneset_oi[nichenet_inputs$geneset_oi %in% colnames(ligand_target_matrix)],
    drop = FALSE
  ]

  message(sprintf(
    "[cci:NicheNet] Top %d ligands identified. Ligand-target matrix: %d x %d.",
    length(top_ligands), nrow(lt_subset), ncol(lt_subset)
  ))

  list(
    ligand_activity       = ligand_activity,
    ligand_target_matrix  = lt_subset,
    top_ligands           = top_ligands,
    # Language rule reminder embedded in output
    interpretation_note   = paste(
      "All interactions are putative/candidate based on computational inference from scRNA-seq.",
      "Conclusions should be framed as: 'X may potentially signal to Y via Z'.",
      "Do NOT write 'direct interaction' without orthogonal protein-level evidence."
    )
  )
}

# ---------------------------------------------------------------------------- #
#  Section 4  Scenario B  MultiNicheNet  disease vs control                    #
# ---------------------------------------------------------------------------- #

# Validate and assemble MultiNicheNet input SummarizedExperiment
# Core rule: sample_col defines the statistical unit
cci_prepare_multinichenet <- function(seurat_obj,
                                       celltype_col  = "cell_type_L2",
                                       sample_col    = "sample",
                                       group_col     = "disease_status",
                                       case_group    = NULL,
                                       control_group = NULL,
                                       min_cells_per_sample = 10L,
                                       species       = c("human", "mouse")) {
  cci_pkg_check(c("multinichenetr", "SummarizedExperiment", "Seurat"), scenario = "MultiNicheNet")
  species <- match.arg(species)

  meta <- seurat_obj@meta.data
  cci_check_required_cols(meta, c(celltype_col, sample_col, group_col))

  # Infer case/control labels if not provided
  groups <- unique(as.character(meta[[group_col]]))
  if (is.null(case_group) || is.null(control_group)) {
    if (length(groups) != 2) {
      stop(sprintf(
        "[cci:MultiNicheNet] Exactly 2 groups required. Found: %s. Specify case_group and control_group.",
        paste(groups, collapse = ", ")
      ))
    }
    case_group    <- groups[1]
    control_group <- groups[2]
    message(sprintf("[cci:MultiNicheNet] Auto-assigned: case='%s'  control='%s'", case_group, control_group))
  }

  # ====  CRITICAL RULE ====
  # Validate replicate structure BEFORE building input object.
  # Samples are the statistical unit. CCI comparisons without >= 3 samples per group
  # per cell type are underpowered and should be excluded, not downsampled at cell level.
  replicate_check <- cci_validate_replicates(
    meta,
    sample_col   = sample_col,
    group_col    = group_col,
    celltype_col = celltype_col,
    min_samples  = 3L
  )

  group_sample_counts <- tapply(
    as.character(meta[[sample_col]]),
    as.character(meta[[group_col]]),
    function(x) length(unique(x))
  )

  # Remove (celltype, group) pairs with < 3 samples from analysis
  if (nrow(replicate_check$inadequate_pairs) > 0) {
    inadequate_ct <- rownames(replicate_check$sample_table)[
      replicate_check$inadequate_pairs[, 1]
    ]
    message(sprintf(
      "[cci:MultiNicheNet] Excluding %d cell types due to insufficient replicates: %s",
      length(unique(inadequate_ct)), paste(unique(inadequate_ct), collapse = ", ")
    ))
    meta <- meta[!(meta[[celltype_col]] %in% unique(inadequate_ct)), , drop = FALSE]
  }

  # Per-celltype, per-sample minimum cell threshold
  cs_counts <- table(meta[[celltype_col]], meta[[sample_col]])
  low_pairs <- which(cs_counts > 0 & cs_counts < min_cells_per_sample, arr.ind = TRUE)
  if (nrow(low_pairs) > 0) {
    message(sprintf(
      "[cci:MultiNicheNet] %d (celltype x sample) pairs have < %d cells; MultiNicheNet will handle internally.",
      nrow(low_pairs), min_cells_per_sample
    ))
  }

  # Load species-appropriate NicheNet networks
  lr_network <- if (species == "human") {
    multinichenetr::lr_network_human
  } else {
    multinichenetr::lr_network_mouse
  }
  ligand_target_matrix <- if (species == "human") {
    multinichenetr::ligand_target_matrix_human
  } else {
    multinichenetr::ligand_target_matrix_mouse
  }

  message(sprintf(
    "[cci:MultiNicheNet] Building input SE: case='%s' control='%s' | %d cell types | %d samples",
    case_group, control_group,
    length(unique(meta[[celltype_col]])),
    length(unique(meta[[sample_col]]))
  ))

  # Build the SummarizedExperiment input
  mn_input <- tryCatch(
    multinichenetr::get_abundance_expression_info(
      sce                    = Seurat::as.SingleCellExperiment(seurat_obj[, rownames(meta)]),
      sample_id              = sample_col,
      group_id               = group_col,
      celltype_id            = celltype_col,
      min_cells              = min_cells_per_sample,
      lr_network             = lr_network,
      senders_oi             = NULL,
      receivers_oi           = NULL
    ),
    error = function(e) stop(sprintf("[cci:MultiNicheNet] get_abundance_expression_info failed: %s", e$message))
  )

  list(
    mn_input              = mn_input,
    case_group            = case_group,
    control_group         = control_group,
    lr_network            = lr_network,
    ligand_target_matrix  = ligand_target_matrix,
    replicate_check       = replicate_check,
    group_sample_counts   = group_sample_counts,
    sample_col            = sample_col,
    group_col             = group_col,
    celltype_col          = celltype_col
  )
}

# Run MultiNicheNet DE-based CCI pipeline
# Returns: prioritized LR interactions with sample-level DE support
cci_run_multinichenet <- function(mn_prep,
                                   contrast_tbl   = NULL,
                                   top_n_ligands  = 25L,
                                   top_n_targets  = 25L,
                                   workers        = 4L) {
  cci_pkg_check("multinichenetr", scenario = "MultiNicheNet")

  # Build default two-group contrast if not supplied
  if (is.null(contrast_tbl)) {
    contrast_tbl <- data.frame(
      contrast      = paste0(mn_prep$case_group, "-", mn_prep$control_group),
      group_oi      = mn_prep$case_group,
      group_ref     = mn_prep$control_group,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    message(sprintf("[cci:MultiNicheNet] Using contrast: %s", contrast_tbl$contrast))
  }

  # Step 1: DE analysis per (sender, receiver, sample) — samples are the unit
  message("[cci:MultiNicheNet] Step 1: pseudo-bulk DE per sender x receiver x sample...")
  de_info <- tryCatch(
    multinichenetr::get_DE_info(
      abundances_expression_info = mn_prep$mn_input,
      sample_id                  = mn_prep$sample_col,
      group_id                   = mn_prep$group_col,
      celltype_id                = mn_prep$celltype_col,
      contrasts_oi               = contrast_tbl$contrast,
      min_cells                  = 10L
    ),
    error = function(e) stop(sprintf("[cci:MultiNicheNet] get_DE_info failed: %s", e$message))
  )

  # Step 2: Filter expressed LR pairs
  message("[cci:MultiNicheNet] Step 2: filtering expressed LR pairs...")
  expressed_info <- tryCatch(
    multinichenetr::get_expressed_ligands_receptors(
      abundances_expression_info = mn_prep$mn_input,
      sample_id                  = mn_prep$sample_col,
      group_id                   = mn_prep$group_col,
      celltype_id                = mn_prep$celltype_col,
      contrasts_oi               = contrast_tbl$contrast,
      fraction_cutoff            = 0.05
    ),
    error = function(e) stop(sprintf("[cci:MultiNicheNet] get_expressed_ligands_receptors failed: %s", e$message))
  )

  # Step 3: Prioritize LR interactions
  message("[cci:MultiNicheNet] Step 3: prioritizing CCI interactions...")
  prioritized_tbl <- tryCatch(
    multinichenetr::generate_prioritization_tables(
      sender_receiver_info       = expressed_info$sender_receiver_info,
      sender_receiver_de         = de_info$celltype_de,
      ligand_activities_targets  = de_info$ligand_activities_targets,
      contrast_tbl               = contrast_tbl,
      sender_receiver_tbl        = expressed_info$sender_receiver_tbl,
      ligand_target_matrix       = mn_prep$ligand_target_matrix,
      top_n_target               = top_n_targets
    ),
    error = function(e) stop(sprintf("[cci:MultiNicheNet] prioritization failed: %s", e$message))
  )

  n_prioritized <- nrow(prioritized_tbl$group_prioritization_tbl)
  message(sprintf("[cci:MultiNicheNet] Prioritization complete: %d LR interactions scored.", n_prioritized))

  n_case <- if (!is.null(mn_prep$group_sample_counts) && mn_prep$case_group %in% names(mn_prep$group_sample_counts)) {
    mn_prep$group_sample_counts[[mn_prep$case_group]]
  } else {
    NA_integer_
  }
  n_ctrl <- if (!is.null(mn_prep$group_sample_counts) && mn_prep$control_group %in% names(mn_prep$group_sample_counts)) {
    mn_prep$group_sample_counts[[mn_prep$control_group]]
  } else {
    NA_integer_
  }

  list(
    prioritized_tbl   = prioritized_tbl,
    de_info           = de_info,
    expressed_info    = expressed_info,
    contrast_tbl      = contrast_tbl,
    replicate_check   = mn_prep$replicate_check,
    group_sample_counts = mn_prep$group_sample_counts,
    # Embedded statistical note for write-up
    statistical_note  = paste(
      sprintf("MultiNicheNet analysis used %s as the statistical unit (case='%s', n_case=%s; control='%s', n_ctrl=%s).",
              mn_prep$sample_col, mn_prep$case_group, n_case, mn_prep$control_group, n_ctrl),
      "Pseudo-bulk DE was computed per (sender, receiver, sample) group.",
      "Cell-level resampling was NOT used as the primary statistical test."
    )
  )
}

# ---------------------------------------------------------------------------- #
#  Section 5  Scenario B  LIANA+ per-sample aggregation                       #
# ---------------------------------------------------------------------------- #

# Run LIANA+ on each sample separately and aggregate across samples
# Returns: aggregated data.frame with per-sample support statistics
cci_run_liana_per_sample <- function(seurat_obj,
                                      celltype_col = "cell_type_L2",
                                      sample_col   = "sample",
                                      group_col    = "disease_status",
                                      methods      = c("natmi", "connectome", "logfc", "sca"),
                                      min_cells    = 10L,
                                      workers      = 4L) {
  cci_pkg_check(c("liana", "dplyr", "Seurat"), scenario = "LIANA per-sample")

  meta <- seurat_obj@meta.data
  cci_check_required_cols(meta, c(celltype_col, sample_col, group_col))

  samples <- sort(unique(as.character(meta[[sample_col]])))
  message(sprintf("[cci:LIANA] Running per-sample LIANA+ across %d samples...", length(samples)))

  per_sample_results <- lapply(samples, function(samp) {
    cells_in_sample <- rownames(meta)[meta[[sample_col]] == samp]
    sub_obj <- seurat_obj[, cells_in_sample]

    # Skip samples with insufficient cell types (< 2 types with >= min_cells)
    ct_counts <- table(sub_obj@meta.data[[celltype_col]])
    valid_cts  <- names(ct_counts[ct_counts >= min_cells])
    if (length(valid_cts) < 2) {
      message(sprintf("[cci:LIANA] Skipping sample '%s': < 2 cell types with >= %d cells.", samp, min_cells))
      return(NULL)
    }

    keep_cells <- cells_in_sample[sub_obj@meta.data[[celltype_col]] %in% valid_cts]
    sub_obj    <- seurat_obj[, keep_cells]
    Seurat::Idents(sub_obj) <- sub_obj@meta.data[[celltype_col]]

    result <- tryCatch(
      liana::liana_wrap(sub_obj, method = methods, resource = "Consensus", verbose = FALSE,
                        parallelize = FALSE),
      error = function(e) {
        message(sprintf("[cci:LIANA] Sample '%s' failed: %s", samp, e$message))
        return(NULL)
      }
    )
    if (is.null(result)) return(NULL)

    agg <- tryCatch(liana::liana_aggregate(result), error = function(e) NULL)
    if (is.null(agg)) return(NULL)

    agg[[sample_col]] <- samp
    agg[[group_col]]  <- unique(meta[[group_col]][meta[[sample_col]] == samp])[1]
    agg
  })

  per_sample_results <- Filter(Negate(is.null), per_sample_results)
  if (length(per_sample_results) == 0) {
    stop("[cci:LIANA] All samples failed. Check cell type counts and package installation.")
  }

  combined <- dplyr::bind_rows(per_sample_results)

  # Compute per-interaction frequency and mean rank across samples
  interaction_key <- paste(combined$source, combined$target,
                           combined$ligand.complex, combined$receptor.complex, sep = "|")
  combined$interaction_key <- interaction_key

  split_df <- split(combined, combined$interaction_key)
  summary_rows <- lapply(split_df, function(x) {
    data.frame(
      interaction_key = x$interaction_key[1],
      source = if ("source" %in% colnames(x)) as.character(x$source[1]) else "",
      target = if ("target" %in% colnames(x)) as.character(x$target[1]) else "",
      ligand.complex = if ("ligand.complex" %in% colnames(x)) as.character(x$ligand.complex[1]) else "",
      receptor.complex = if ("receptor.complex" %in% colnames(x)) as.character(x$receptor.complex[1]) else "",
      n_samples_detected = nrow(x),
      mean_aggregate_rank = if ("aggregate_rank" %in% colnames(x)) mean(cci_num(x$aggregate_rank), na.rm = TRUE) else NA_real_,
      groups_detected = if (group_col %in% colnames(x)) paste(sort(unique(as.character(x[[group_col]]))), collapse = "/") else "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  summary_df <- cci_bind_rows_flexible(summary_rows)
  if (nrow(summary_df) > 0L) {
    summary_df <- summary_df[order(summary_df$mean_aggregate_rank, na.last = TRUE), , drop = FALSE]
  }

  # Fraction of samples with detection
  n_total_samples <- length(samples)
  summary_df$detection_fraction <- summary_df$n_samples_detected / n_total_samples

  message(sprintf(
    "[cci:LIANA] Per-sample aggregation: %d unique LR pairs across %d samples.",
    nrow(summary_df), length(per_sample_results)
  ))

  list(
    per_sample_raw = combined,
    summary        = summary_df
  )
}

# ---------------------------------------------------------------------------- #
#  Section 6  Scenario B  CellChat global network visualization                #
# ---------------------------------------------------------------------------- #

# Build and compare CellChat objects for disease vs control
# For visualization only — statistical claims must come from MultiNicheNet
cci_run_cellchat_global_viz <- function(seurat_obj,
                                         celltype_col  = "cell_type_L2",
                                         group_col     = "disease_status",
                                         case_group    = NULL,
                                         control_group = NULL,
                                         species       = c("human", "mouse"),
                                         nboot         = 100L,
                                         workers       = 4L) {
  cci_pkg_check("CellChat", scenario = "CellChat global viz")
  species <- match.arg(species)

  meta   <- seurat_obj@meta.data
  groups <- unique(as.character(meta[[group_col]]))

  if (is.null(case_group) || is.null(control_group)) {
    if (length(groups) != 2) {
      stop("[cci:CellChat-viz] Exactly 2 groups required. Specify case_group and control_group.")
    }
    case_group    <- groups[1]
    control_group <- groups[2]
  }

  run_one_group <- function(grp_label) {
    cells <- rownames(meta)[meta[[group_col]] == grp_label]
    sub_obj <- seurat_obj[, cells]
    cc <- cci_prepare_cellchat(
      sub_obj, celltype_col = celltype_col, species = species
    )
    cc <- cci_run_cellchat_screening(cc, nboot = nboot, workers = workers)
    cc@meta$group <- grp_label
    cc
  }

  message(sprintf("[cci:CellChat-viz] Building CellChat for '%s'...", case_group))
  cc_case <- run_one_group(case_group)

  message(sprintf("[cci:CellChat-viz] Building CellChat for '%s'...", control_group))
  cc_ctrl <- run_one_group(control_group)

  # Lift-over to shared cell types for comparison
  cc_list <- stats::setNames(list(cc_case, cc_ctrl), c(case_group, control_group))
  cellchat_merged <- tryCatch(
    CellChat::mergeCellChat(
      cc_list,
      add.names = names(cc_list)
    ),
    error = function(e) stop(sprintf("[cci:CellChat-viz] mergeCellChat failed: %s", e$message))
  )

  message("[cci:CellChat-viz] CellChat objects merged for comparative visualization.")
  message("[cci:CellChat-viz] NOTE: Use CellChat outputs for visualization only.")
  message("[cci:CellChat-viz] Statistical claims require MultiNicheNet (sample-level DE).")

  list(
    cc_case    = cc_case,
    cc_control = cc_ctrl,
    cc_merged  = cellchat_merged,
    case_label = case_group,
    ctrl_label = control_group,
    viz_note   = paste(
      "CellChat comparison figures show global network structure differences.",
      "These are exploratory and cannot substitute for sample-replicated DE analysis (MultiNicheNet).",
      "Do not compute p-values by treating cells as independent replicates."
    )
  )
}

# ---------------------------------------------------------------------------- #
#  Section 7  Output and Reporting                                              #
# ---------------------------------------------------------------------------- #

# Export Scenario A results: CellChat + LIANA consensus + NicheNet summary
cci_export_scrna_results <- function(cellchat_candidates,
                                      liana_consensus,
                                      nichenet_result,
                                      output_dir,
                                      prefix = "cci_scrna",
                                      shared_candidates = NULL,
                                      write_visuals = TRUE,
                                      llm_config = list(),
                                      top_n_llm = 80L) {
  cci_safe_mkdir(output_dir)
  if (is.null(cellchat_candidates)) cellchat_candidates <- data.frame()
  if (is.null(liana_consensus)) liana_consensus <- data.frame()
  if (is.null(shared_candidates)) shared_candidates <- data.frame()
  if (is.null(nichenet_result) || !is.list(nichenet_result)) {
    nichenet_result <- list(ligand_activity = data.frame(), ligand_target_matrix = matrix(), interpretation_note = "NicheNet not available.")
  }
  nn_ligand_activity <- if (is.data.frame(nichenet_result$ligand_activity)) nichenet_result$ligand_activity else data.frame()

  # CellChat candidates
  cc_path <- file.path(output_dir, paste0(prefix, "_cellchat_candidates.tsv"))
  cc_csv_path <- file.path(output_dir, paste0(prefix, "_cellchat_candidates.csv"))
  cci_write_tsv(cellchat_candidates, cc_path)
  cci_write_csv(cellchat_candidates, cc_csv_path)

  # LIANA consensus
  li_path <- file.path(output_dir, paste0(prefix, "_liana_consensus.tsv"))
  li_csv_path <- file.path(output_dir, paste0(prefix, "_liana_consensus.csv"))
  cci_write_tsv(liana_consensus, li_path)
  cci_write_csv(liana_consensus, li_csv_path)

  shared_path <- file.path(output_dir, paste0(prefix, "_cellchat_liana_shared_candidates.tsv"))
  shared_csv_path <- file.path(output_dir, paste0(prefix, "_cellchat_liana_shared_candidates.csv"))
  cci_write_tsv(shared_candidates, shared_path)
  cci_write_csv(shared_candidates, shared_csv_path)

  # NicheNet ligand activity
  nn_act_path <- file.path(output_dir, paste0(prefix, "_nichenet_ligand_activity.tsv"))
  nn_act_csv_path <- file.path(output_dir, paste0(prefix, "_nichenet_ligand_activity.csv"))
  cci_write_tsv(nn_ligand_activity, nn_act_path)
  cci_write_csv(nn_ligand_activity, nn_act_csv_path)

  # NicheNet ligand-target matrix
  nn_lt_path <- file.path(output_dir, paste0(prefix, "_nichenet_ligand_target_matrix.rds"))
  saveRDS(nichenet_result$ligand_target_matrix, nn_lt_path)

  # LLM-ready interaction CSV and paired LLM interpretation files
  llm_tables <- list(
    CellChat_candidates = list(
      data = cellchat_candidates,
      method = "CellChat candidates",
      evidence_scope = "candidate LR screening",
      statistical_unit = "cells for screening only",
      confidence_default = "putative"
    ),
    LIANA_consensus = list(
      data = liana_consensus,
      method = "LIANA+ consensus",
      evidence_scope = "multi-method consensus",
      statistical_unit = "cells for screening only",
      confidence_default = "putative"
    ),
    CellChat_LIANA_shared = list(
      data = shared_candidates,
      method = "CellChat ∩ LIANA+ shared candidates",
      evidence_scope = "cross-method support",
      statistical_unit = "cells for screening only",
      confidence_default = "candidate"
    ),
    NicheNet_ligand_activity = list(
      data = nn_ligand_activity,
      method = "NicheNet ligand activity",
      evidence_scope = "receiver downstream response",
      statistical_unit = "cells for receiver-response inference",
      confidence_default = "putative"
    )
  )
  llm_df <- cci_build_llm_interaction_table(
    llm_tables,
    scenario = "A_scrna_only_no_spatial_no_protein",
    default_statistical_unit = "cells for screening only",
    top_n_per_method = top_n_llm
  )
  llm_csv_path <- file.path(output_dir, paste0(prefix, "_LLM_interactions.csv"))
  cci_write_csv(llm_df, llm_csv_path)
  llm_status <- cci_run_table_llm_interpretation(
    df = llm_df,
    output_dir = output_dir,
    prefix = paste0(prefix, "_interactions"),
    title = "Scenario A scRNA-only CCI candidate review",
    scenario = "A_scrna_only_no_spatial_no_protein",
    extra_context = c(
      "Recommended route: CellPhoneDB/CellChat candidate screening, LIANA+ multi-method consensus, NicheNet receiver downstream response.",
      "Required language: potential / candidate / putative; no direct interaction claims without spatial or protein validation."
    ),
    llm_config = llm_config
  )

  visualization_manifest <- NULL
  if (isTRUE(write_visuals)) {
    visualization_manifest <- cci_plot_scrna_visualizations(
      cellchat_candidates = cellchat_candidates,
      liana_consensus = liana_consensus,
      nichenet_result = nichenet_result,
      output_dir = output_dir,
      prefix = prefix,
      llm_config = llm_config
    )
  }

  # Language note
  note_path <- file.path(output_dir, paste0(prefix, "_INTERPRETATION_NOTE.txt"))
  writeLines(c(
    "CELL COMMUNICATION INTERPRETATION GUIDELINES",
    "============================================",
    "",
    "Source: scRNA-seq only (no spatial or protein validation)",
    "",
    "Required language for manuscript writing:",
    "  - Use: 'putative', 'candidate', 'potential', 'computationally predicted'",
    "  - Avoid: 'direct interaction', 'confirmed signaling', 'demonstrated communication'",
    "",
    nichenet_result$interpretation_note,
    "",
    sprintf("Generated: %s", Sys.time())
  ), note_path)

  message(sprintf("[cci] Scenario A results exported to: %s", output_dir))
  manifest_path <- file.path(output_dir, paste0(prefix, "_manifest.json"))
  manifest <- list(
    scenario = "A_scrna_only_no_spatial_no_protein",
    output_dir = cci_path_normalize(output_dir),
    cellchat_path = cci_path_normalize(cc_path),
    cellchat_csv = cci_path_normalize(cc_csv_path),
    liana_path = cci_path_normalize(li_path),
    liana_csv = cci_path_normalize(li_csv_path),
    shared_path = cci_path_normalize(shared_path),
    shared_csv = cci_path_normalize(shared_csv_path),
    nichenet_paths = list(activity = cci_path_normalize(nn_act_path), activity_csv = cci_path_normalize(nn_act_csv_path), lt_matrix = cci_path_normalize(nn_lt_path)),
    llm_interactions_csv = cci_path_normalize(llm_csv_path),
    llm = llm_status,
    visualizations = visualization_manifest,
    note_path = cci_path_normalize(note_path),
    helper_version = "20260506_v1.1"
  )
  cci_write_json(manifest, manifest_path)
  manifest$manifest_json <- cci_path_normalize(manifest_path)
  invisible(list(
    cellchat_path  = cc_path,
    cellchat_csv   = cc_csv_path,
    liana_path     = li_path,
    liana_csv      = li_csv_path,
    shared_path    = shared_path,
    shared_csv     = shared_csv_path,
    nichenet_paths = list(activity = nn_act_path, activity_csv = nn_act_csv_path, lt_matrix = nn_lt_path),
    llm_interactions_csv = llm_csv_path,
    llm_status     = llm_status,
    visualizations = visualization_manifest,
    note_path      = note_path,
    manifest       = manifest
  ))
}

# Export Scenario B results: MultiNicheNet + LIANA per-sample + CellChat viz
cci_export_disease_control_results <- function(mn_result,
                                                liana_per_sample,
                                                cellchat_viz,
                                                output_dir,
                                                prefix = "cci_disease_ctrl",
                                                write_visuals = TRUE,
                                                llm_config = list(),
                                                top_n_llm = 80L) {
  cci_safe_mkdir(output_dir)
  mn_df <- cci_extract_multinichenet_prioritized_df(mn_result)
  liana_summary <- if (!is.null(liana_per_sample) && is.data.frame(liana_per_sample$summary)) liana_per_sample$summary else data.frame()

  # MultiNicheNet prioritized table
  mn_path <- file.path(output_dir, paste0(prefix, "_multinichenet_prioritized.tsv"))
  mn_csv_path <- file.path(output_dir, paste0(prefix, "_multinichenet_prioritized.csv"))
  cci_write_tsv(mn_df, mn_path)
  cci_write_csv(mn_df, mn_csv_path)

  # LIANA per-sample summary
  li_path <- file.path(output_dir, paste0(prefix, "_liana_per_sample_summary.tsv"))
  li_csv_path <- file.path(output_dir, paste0(prefix, "_liana_per_sample_summary.csv"))
  cci_write_tsv(liana_summary, li_path)
  cci_write_csv(liana_summary, li_csv_path)

  replicate_df <- if (is.list(mn_result) && !is.null(mn_result$replicate_check)) cci_replicate_check_to_df(mn_result$replicate_check) else data.frame()
  replicate_csv_path <- file.path(output_dir, paste0(prefix, "_replicate_adequacy.csv"))
  cci_write_csv(replicate_df, replicate_csv_path)

  # Save CellChat objects for visualization
  cc_path <- file.path(output_dir, paste0(prefix, "_cellchat_merged.rds"))
  cc_obj_to_save <- if (is.list(cellchat_viz) && !is.null(cellchat_viz$cc_merged)) cellchat_viz$cc_merged else cellchat_viz
  saveRDS(cc_obj_to_save, cc_path)

  # LLM-ready interaction CSV and paired LLM interpretation files
  llm_df <- cci_build_llm_interaction_table(
    list(
      MultiNicheNet_prioritized = list(
        data = mn_df,
        method = "MultiNicheNet prioritized",
        evidence_scope = "sample-level DE-based CCI prioritization",
        statistical_unit = "sample / donor",
        confidence_default = "candidate"
      ),
      LIANA_per_sample = list(
        data = liana_summary,
        method = "LIANA+ per-sample support",
        evidence_scope = "auxiliary per-sample aggregation",
        statistical_unit = "sample / donor",
        confidence_default = "putative"
      )
    ),
    scenario = "B_disease_vs_control_with_biological_replicates",
    default_statistical_unit = "sample / donor",
    top_n_per_method = top_n_llm
  )
  llm_csv_path <- file.path(output_dir, paste0(prefix, "_LLM_interactions.csv"))
  cci_write_csv(llm_df, llm_csv_path)
  llm_status <- cci_run_table_llm_interpretation(
    df = llm_df,
    output_dir = output_dir,
    prefix = paste0(prefix, "_interactions"),
    title = "Scenario B disease-control CCI review",
    scenario = "B_disease_vs_control_with_biological_replicates",
    extra_context = c(
      "Primary route: MultiNicheNet sample-level DE-based CCI; LIANA+ per-sample support is auxiliary; CellChat is global network visualization only.",
      "Core statistical rule: samples/donors are biological replicates; cells are not independent replicates."
    ),
    llm_config = llm_config
  )

  visualization_manifest <- NULL
  if (isTRUE(write_visuals)) {
    visualization_manifest <- cci_plot_disease_control_visualizations(
      mn_result = mn_result,
      liana_per_sample = liana_per_sample,
      cellchat_viz = cellchat_viz,
      output_dir = output_dir,
      prefix = prefix,
      llm_config = llm_config
    )
  }

  # Statistical and methods note for Methods section
  methods_path <- file.path(output_dir, paste0(prefix, "_METHODS_NOTE.txt"))
  writeLines(c(
    "STATISTICAL DESIGN NOTE",
    "=======================",
    "",
    mn_result$statistical_note,
    "",
    "VISUALIZATION NOTE",
    "==================",
    "",
    cellchat_viz$viz_note,
    "",
    "RECOMMENDED METHODS LANGUAGE",
    "============================",
    "",
    "Primary: MultiNicheNet v[X.X] was used to identify differentially active",
    "  cell-cell communication events between [case] and [control] groups.",
    "  Pseudo-bulk DE analysis was performed using samples as the statistical unit",
    "  (n = [N] per group).",
    "",
    "Auxiliary: LIANA+ (multi-method consensus: natmi, connectome, logfc, sca)",
    "  was run per sample and interactions detected in >= 50% of samples per group",
    "  were retained for cross-validation.",
    "",
    "Visualization: CellChat v2 was used to construct global interaction networks",
    "  for visual comparison. Network-level statistics (number/weight of interactions)",
    "  are shown as exploratory summaries only.",
    "",
    sprintf("Generated: %s", Sys.time())
  ), methods_path)

  message(sprintf("[cci] Scenario B results exported to: %s", output_dir))
  manifest_path <- file.path(output_dir, paste0(prefix, "_manifest.json"))
  manifest <- list(
    scenario = "B_disease_vs_control_with_biological_replicates",
    output_dir = cci_path_normalize(output_dir),
    mn_path = cci_path_normalize(mn_path),
    mn_csv = cci_path_normalize(mn_csv_path),
    liana_path = cci_path_normalize(li_path),
    liana_csv = cci_path_normalize(li_csv_path),
    replicate_adequacy_csv = cci_path_normalize(replicate_csv_path),
    cc_path = cci_path_normalize(cc_path),
    llm_interactions_csv = cci_path_normalize(llm_csv_path),
    llm = llm_status,
    visualizations = visualization_manifest,
    methods_path = cci_path_normalize(methods_path),
    helper_version = "20260506_v1.1"
  )
  cci_write_json(manifest, manifest_path)
  manifest$manifest_json <- cci_path_normalize(manifest_path)
  invisible(list(
    mn_path      = mn_path,
    mn_csv       = mn_csv_path,
    liana_path   = li_path,
    liana_csv    = li_csv_path,
    replicate_csv = replicate_csv_path,
    cc_path      = cc_path,
    llm_interactions_csv = llm_csv_path,
    llm_status   = llm_status,
    visualizations = visualization_manifest,
    methods_path = methods_path,
    manifest     = manifest
  ))
}

# Write a markdown summary report for either scenario
cci_write_markdown_report <- function(scenario    = c("A_scrna", "B_disease_control"),
                                       results_dir,
                                       prefix      = "cci",
                                       celltype    = "all",
                                       extra_notes = character()) {
  scenario <- match.arg(scenario)

  report_lines <- c(
    "# Cell-Cell Communication Analysis Report",
    sprintf("**Scenario:** %s", if (scenario == "A_scrna") "A — scRNA-seq only" else "B — Disease vs Control"),
    sprintf("**Date:** %s", Sys.Date()),
    sprintf("**Cell type(s):** %s", celltype),
    "",
    "## Methods Summary",
    "",
    if (scenario == "A_scrna") c(
      "1. **CellChat** (candidate LR screening, n_boot = 100)",
      "2. **LIANA+** (multi-method consensus: natmi, connectome, logfc, sca, cellphonedb)",
      "3. **NicheNet** (receiver-side downstream target prediction)",
      "",
      "> Language note: All interactions are described as *putative* or *candidate*.",
      "> The term 'direct interaction' is not used without orthogonal evidence."
    ) else c(
      "1. **MultiNicheNet** (primary; DE-based CCI using samples as statistical unit)",
      "2. **LIANA+** (auxiliary; per-sample aggregation and cross-validation)",
      "3. **CellChat** (global network visualization only; not used for statistical claims)",
      "",
      "> Statistical unit: **samples** (biological replicates), not cells.",
      "> Cell-level pseudoreplication was explicitly avoided."
    ),
    "",
    "## Output Files",
    sprintf("- Results directory: `%s`", results_dir),
        sprintf("- LLM interaction CSV: `%s`", file.path(results_dir, paste0(prefix, "_LLM_interactions.csv"))),
        sprintf("- LLM prompt/status: `%s`, `%s`",
          file.path(results_dir, paste0(prefix, "_interactions_LLM_PROMPT.md")),
          file.path(results_dir, paste0(prefix, "_interactions_LLM_STATUS.json"))),
        sprintf("- Figure companion manifest: `%s`", file.path(results_dir, "figures", paste0(prefix, "_figure_companion_manifest.tsv"))),
    "",
    if (length(extra_notes) > 0) c("## Notes", extra_notes) else character()
  )

  report_path <- file.path(results_dir, paste0(prefix, "_REPORT.md"))
  writeLines(report_lines, report_path)
  message(sprintf("[cci] Markdown report written to: %s", report_path))
  invisible(report_path)
}

# ---------------------------------------------------------------------------- #
#  Section 8  Quick-start wrappers                                             #
# ---------------------------------------------------------------------------- #

# Scenario A one-call wrapper: CellChat -> LIANA+ -> NicheNet -> export
cci_run_scrna_pipeline <- function(seurat_obj,
                                    output_dir,
                                    celltype_col        = "cell_type_L2",
                                    receiver_celltype   = NULL,
                                    condition_col       = NULL,
                                    condition_test      = NULL,
                                    condition_reference = NULL,
                                    species             = "human",
                                    workers             = 4L,
                                    prefix              = "cci_scrna",
                                    write_visuals       = TRUE,
                                    llm_config          = list()) {
  cci_safe_mkdir(output_dir)
  message("[cci] === Scenario A: scRNA-only CCI pipeline ===")

  # Step 1: CellChat
  cc_obj   <- cci_prepare_cellchat(seurat_obj, celltype_col = celltype_col, species = species)
  cc_obj   <- cci_run_cellchat_screening(cc_obj, workers = workers)
  cc_cands <- cci_extract_cellchat_candidates(cc_obj)

  # Step 2: LIANA+ consensus
  liana_res <- cci_run_liana_consensus(seurat_obj, celltype_col = celltype_col,
                                       species = species, workers = workers)

  # Step 3: Intersection
  shared <- cci_intersect_cellchat_liana(cc_cands, liana_res)

  # Step 4: NicheNet (optional, requires receiver specification)
  nn_result <- NULL
  if (!is.null(receiver_celltype) && !is.null(condition_test)) {
    nn_inputs <- cci_prepare_nichenet_inputs(
      seurat_obj,
      receiver_celltype   = receiver_celltype,
      celltype_col        = celltype_col,
      condition_col       = condition_col,
      condition_test      = condition_test,
      condition_reference = condition_reference
    )
    nn_result <- cci_run_nichenet(seurat_obj, nn_inputs, species = species)
  } else {
    message("[cci] Skipping NicheNet: receiver_celltype or condition_test not specified.")
    nn_result <- list(
      ligand_activity      = data.frame(),
      ligand_target_matrix = matrix(),
      top_ligands          = character(),
      interpretation_note  = "NicheNet not run (no receiver or condition specified)."
    )
  }

  # Step 5: Export
  paths <- cci_export_scrna_results(
    cc_cands,
    liana_res,
    nn_result,
    output_dir,
    prefix = prefix,
    shared_candidates = shared,
    write_visuals = write_visuals,
    llm_config = llm_config
  )
  cci_write_markdown_report("A_scrna", output_dir, prefix = prefix)

  message("[cci] === Scenario A complete ===")
  invisible(list(
    cellchat_obj        = cc_obj,
    cellchat_candidates = cc_cands,
    liana_consensus     = liana_res,
    shared_candidates   = shared,
    nichenet            = nn_result,
    output_paths        = paths
  ))
}

# Scenario B one-call wrapper: MultiNicheNet -> LIANA per-sample -> CellChat viz -> export
cci_run_disease_control_pipeline <- function(seurat_obj,
                                              output_dir,
                                              celltype_col  = "cell_type_L2",
                                              sample_col    = "sample",
                                              group_col     = "disease_status",
                                              case_group    = NULL,
                                              control_group = NULL,
                                              species       = "human",
                                              workers       = 4L,
                                              prefix        = "cci_disease_ctrl",
                                              write_visuals = TRUE,
                                              llm_config    = list()) {
  cci_safe_mkdir(output_dir)
  message("[cci] === Scenario B: Disease vs Control CCI pipeline ===")
  message("[cci] Statistical unit: SAMPLES (biological replicates), NOT cells.")

  # Step 1: MultiNicheNet
  mn_prep   <- cci_prepare_multinichenet(seurat_obj, celltype_col = celltype_col,
                                         sample_col = sample_col, group_col = group_col,
                                         case_group = case_group, control_group = control_group,
                                         species = species)
  mn_result <- cci_run_multinichenet(mn_prep, workers = workers)

  # Step 2: LIANA per-sample (auxiliary)
  liana_ps <- cci_run_liana_per_sample(seurat_obj, celltype_col = celltype_col,
                                        sample_col = sample_col, group_col = group_col,
                                        workers = workers)

  # Step 3: CellChat global visualization
  cc_viz <- cci_run_cellchat_global_viz(seurat_obj, celltype_col = celltype_col,
                                         group_col = group_col,
                                         case_group    = mn_prep$case_group,
                                         control_group = mn_prep$control_group,
                                         species = species, workers = workers)

  # Step 4: Export
  paths <- cci_export_disease_control_results(
    mn_result,
    liana_ps,
    cc_viz,
    output_dir,
    prefix = prefix,
    write_visuals = write_visuals,
    llm_config = llm_config
  )
  cci_write_markdown_report("B_disease_control", output_dir, prefix = prefix)

  message("[cci] === Scenario B complete ===")
  invisible(list(
    multinichenet     = mn_result,
    liana_per_sample  = liana_ps,
    cellchat_viz      = cc_viz,
    output_paths      = paths
  ))
}

message("[cci] cell_communication_helper_20260506_v1_0.R loaded (helper version 20260506_v1.1).")
message("[cci] Scenario A: cci_run_scrna_pipeline()")
message("[cci] Scenario B: cci_run_disease_control_pipeline()")
