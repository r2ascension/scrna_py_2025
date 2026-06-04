#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

DEFAULTS <- list(
  run_root = "/home/h2048/data/R/20260531/epithelial_airway_deg_ppi_workflow_20260531",
  output_dir = NA_character_,
  target_pairs = NULL,
  top_seed_n = 20L,
  top_hub_n = 12L,
  top_edge_n = 12L,
  enable_live = TRUE,
  model = "deepseek-reasoner",
  timeout_sec = 240,
  write_env_placeholder = TRUE
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) return(y)
  x
}

safe_trim <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x[[1L]])
}

safe_name <- function(x) {
  x <- safe_trim(x)
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

path_normalize <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_markdown <- function(lines, path) {
  dir_create(dirname(path))
  writeLines(as.character(lines), con = path, useBytes = TRUE)
  invisible(path)
}

write_json <- function(x, path) {
  dir_create(dirname(path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    dput(x, file = path)
  }
  invisible(path)
}

markdown_table <- function(df, max_rows = 12L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) return("No rows available.")
  df <- utils::head(df, max_rows)
  df[] <- lapply(df, function(col) {
    col <- as.character(col)
    col[is.na(col)] <- ""
    gsub("\\|", "/", trimws(col))
  })
  header <- paste(c("", colnames(df), ""), collapse = "|")
  sep <- paste(c("", rep("---", ncol(df)), ""), collapse = "|")
  rows <- apply(df, 1, function(row) paste(c("", row, ""), collapse = "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

print_usage <- function(defaults) {
  cat(
    paste(
      "Usage:",
      "  Rscript epithelial_airway_deg_ppi_llm_bridge_20260531.R [--run-root PATH] [--output-dir PATH]",
      "         [--target-pairs PAIR1,PAIR2] [--top-seed-n INT] [--top-hub-n INT] [--top-edge-n INT]",
      "         [--enable-live true|false] [--model NAME] [--timeout-sec INT] [--write-env-placeholder true|false]",
      "",
      "Defaults:",
      paste0("  --run-root ", defaults$run_root),
      paste0("  --output-dir ", ifelse(is.na(defaults$output_dir), "<run-root>/llm", defaults$output_dir)),
      paste0("  --top-seed-n ", defaults$top_seed_n),
      paste0("  --top-hub-n ", defaults$top_hub_n),
      paste0("  --top-edge-n ", defaults$top_edge_n),
      paste0("  --enable-live ", tolower(as.character(defaults$enable_live))),
      paste0("  --model ", defaults$model),
      paste0("  --timeout-sec ", defaults$timeout_sec),
      paste0("  --write-env-placeholder ", tolower(as.character(defaults$write_env_placeholder))),
      sep = "\n"
    ),
    "\n"
  )
}

parse_bool <- function(x, default = FALSE) {
  x <- tolower(safe_trim(x))
  if (!nzchar(x)) return(isTRUE(default))
  x %in% c("1", "true", "t", "yes", "y", "on")
}

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- defaults
  if (!length(args)) {
    opts$output_dir <- file.path(opts$run_root, "llm")
    return(opts)
  }
  if (any(args %in% c("-h", "--help"))) {
    print_usage(defaults)
    quit(save = "no", status = 0)
  }

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    value <- args[[i + 1L]]
    key <- sub("^--", "", key)
    switch(
      key,
      "run-root" = opts$run_root <- value,
      "output-dir" = opts$output_dir <- value,
      "target-pairs" = {
        value <- trimws(value)
        opts$target_pairs <- if (!nzchar(value)) NULL else trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
      },
      "top-seed-n" = opts$top_seed_n <- as.integer(value),
      "top-hub-n" = opts$top_hub_n <- as.integer(value),
      "top-edge-n" = opts$top_edge_n <- as.integer(value),
      "enable-live" = opts$enable_live <- parse_bool(value, default = opts$enable_live),
      "model" = opts$model <- value,
      "timeout-sec" = opts$timeout_sec <- as.numeric(value),
      "write-env-placeholder" = opts$write_env_placeholder <- parse_bool(value, default = opts$write_env_placeholder),
      stop("Unknown argument: --", key, call. = FALSE)
    )
    i <- i + 2L
  }

  if (is.na(opts$output_dir) || !nzchar(safe_trim(opts$output_dir))) {
    opts$output_dir <- file.path(opts$run_root, "llm")
  }
  opts
}

validate_config <- function(config) {
  config$run_root <- path_normalize(config$run_root)
  config$output_dir <- path_normalize(config$output_dir)
  if (!dir.exists(config$run_root)) stop("run_root does not exist: ", config$run_root, call. = FALSE)
  for (nm in c("top_seed_n", "top_hub_n", "top_edge_n")) {
    value <- as.integer(config[[nm]])
    if (is.na(value) || value < 1L) stop("--", gsub("_", "-", nm), " must be a positive integer", call. = FALSE)
    config[[nm]] <- value
  }
  config$timeout_sec <- suppressWarnings(as.numeric(config$timeout_sec))
  if (is.na(config$timeout_sec) || config$timeout_sec <= 0) config$timeout_sec <- 240
  config$model <- safe_trim(config$model)
  if (!nzchar(config$model)) config$model <- "deepseek-reasoner"
  if (!is.null(config$target_pairs) && length(config$target_pairs) > 0L) {
    config$target_pairs <- unique(trimws(as.character(config$target_pairs)))
    config$target_pairs <- config$target_pairs[nzchar(config$target_pairs)]
  } else {
    config$target_pairs <- NULL
  }
  config
}

load_env_file_safely <- function(path) {
  if (!file.exists(path)) return(FALSE)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*(#|$)", lines)]
  for (ln in lines) {
    if (!grepl("=", ln, fixed = TRUE)) next
    key <- trimws(sub("=.*$", "", ln))
    val <- sub("^[^=]*=", "", ln)
    val <- trimws(val)
    val <- sub("^export\\s+", "", val)
    val <- sub("^([\"'])(.*)\\1$", "\\2", val)
    if (nzchar(key) && !nzchar(Sys.getenv(key, unset = ""))) {
      do.call(Sys.setenv, stats::setNames(list(val), key))
    }
  }
  TRUE
}

is_placeholder_secret <- function(x) {
  x <- tolower(safe_trim(x))
  if (!nzchar(x)) return(TRUE)
  x %in% c(
    "your-deepseek-api-key",
    "your_deepseek_api_key_here",
    "your_deepseek_api_key",
    "your-key",
    "replace_me",
    "changeme"
  ) || grepl("your|placeholder|api[_-]?key[_-]?here|dummy|example", x, ignore.case = TRUE)
}

ensure_env_placeholder_if_missing <- function(env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env"), key = "DEEPSEEK_API_KEY") {
  existing <- env_candidates[file.exists(env_candidates)]
  if (length(existing) > 0L) return(invisible(FALSE))
  target <- env_candidates[[1L]]
  dir_create(dirname(target))
  writeLines(sprintf("%s=your_deepseek_api_key_here", key), con = target, useBytes = TRUE)
  invisible(TRUE)
}

resolve_deepseek_key <- function(api_key = NULL,
                                 env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env"),
                                 write_env_placeholder = FALSE) {
  env_candidates <- as.character(env_candidates)
  existing_envs <- env_candidates[file.exists(env_candidates)]
  invisible(lapply(existing_envs, load_env_file_safely))
  placeholder_written <- FALSE
  if (length(existing_envs) == 0L && isTRUE(write_env_placeholder)) {
    placeholder_written <- ensure_env_placeholder_if_missing(env_candidates = env_candidates)
  }
  key <- if (is.null(api_key) || length(api_key) == 0L || !nzchar(safe_trim(api_key))) {
    Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  } else {
    api_key[[1L]]
  }
  key <- safe_trim(key)
  has_live_key <- nchar(key) >= 20L && !is_placeholder_secret(key)
  list(
    api_key = if (has_live_key) key else "",
    has_live_key = has_live_key,
    env_files_loaded = existing_envs,
    placeholder_written = placeholder_written,
    placeholder_path = if (placeholder_written) env_candidates[[1L]] else NA_character_,
    status = if (has_live_key) "ready" else "missing_or_placeholder_key"
  )
}

call_deepseek_chat <- function(prompt,
                               system_prompt = NULL,
                               model = "deepseek-reasoner",
                               api_key,
                               timeout_sec = 240,
                               base_url = "https://api.deepseek.com/v1/chat/completions") {
  api_key <- safe_trim(api_key)
  if (!nzchar(api_key)) stop("DEEPSEEK_API_KEY not set", call. = FALSE)

  if (requireNamespace("httr2", quietly = TRUE)) {
    messages <- list()
    if (!is.null(system_prompt) && nzchar(safe_trim(system_prompt))) {
      messages[[length(messages) + 1L]] <- list(role = "system", content = as.character(system_prompt))
    }
    messages[[length(messages) + 1L]] <- list(role = "user", content = as.character(prompt))
    req <- httr2::request(base_url) |>
      httr2::req_headers(
        "Content-Type" = "application/json",
        "Authorization" = paste("Bearer", api_key)
      ) |>
      httr2::req_body_json(list(
        model = model,
        messages = messages,
        stream = FALSE
      ))
    if (is.finite(as.numeric(timeout_sec)) && as.numeric(timeout_sec) > 0) {
      req <- httr2::req_timeout(req, seconds = as.numeric(timeout_sec))
    }
    resp <- httr2::req_perform(req)
    body <- httr2::resp_body_json(resp, simplifyVector = FALSE)
    if (httr2::resp_status(resp) >= 300L) {
      msg <- body$error$message %||% httr2::resp_status_desc(resp)
      stop(sprintf("DeepSeek API HTTP %s: %s", httr2::resp_status(resp), msg), call. = FALSE)
    }
    text <- body$choices[[1L]]$message$content %||% ""
    text <- safe_trim(text)
    if (!nzchar(text)) stop("DeepSeek API returned an empty response", call. = FALSE)
    return(text)
  }

  if (requireNamespace("fanyi", quietly = TRUE)) {
    return(safe_trim(fanyi::chat_request(prompt, model = model, api_key = api_key)))
  }

  stop("Neither httr2 nor fanyi is available for live LLM calls", call. = FALSE)
}

find_seed_summary_path <- function(run_root) {
  candidates <- list.files(
    run_root,
    pattern = "^seed_summary_by_pair_top[0-9]+\\.tsv$",
    full.names = TRUE
  )
  if (!length(candidates)) stop("No seed_summary_by_pair_top*.tsv found under run_root", call. = FALSE)
  top_n <- suppressWarnings(as.integer(sub("^.*top([0-9]+)\\.tsv$", "\\1", candidates)))
  candidates[[order(top_n, decreasing = TRUE)[[1L]]]]
}

read_required_tsv <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
  fread(path, sep = "\t", showProgress = FALSE)
}

read_optional_tsv <- function(path) {
  if (!file.exists(path)) return(data.table())
  fread(path, sep = "\t", showProgress = FALSE)
}

split_semicolon_values <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vals <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  vals[nzchar(vals)]
}

summarize_seed_context <- function(pair_seed_dt, top_n = 20L) {
  pair_seed_dt <- copy(pair_seed_dt)
  if (!nrow(pair_seed_dt)) {
    return(list(
      direction_summary = data.table(direction = character(), n = integer()),
      celltype_summary = data.table(cell_type = character(), n = integer()),
      seed_preview = data.table(),
      seed_text = "No seed genes available."
    ))
  }
  preview_dt <- pair_seed_dt[order(seed_rank)][seq_len(min(top_n, .N)), .(
    seed_rank,
    gene_symbol,
    dominant_direction,
    n_cell_types,
    n_occurrences,
    cell_types,
    max_abs_log2FoldChange,
    min_padj,
    ppi_priority_score
  )]
  direction_summary <- pair_seed_dt[, .N, by = dominant_direction][order(-N, dominant_direction)]
  setnames(direction_summary, "N", "n")
  celltype_values <- split_semicolon_values(pair_seed_dt$cell_types)
  celltype_summary <- if (length(celltype_values)) {
    as.data.table(sort(table(celltype_values), decreasing = TRUE), keep.rownames = "cell_type")
  } else {
    data.table(cell_type = character(), N = integer())
  }
  if (!"cell_type" %in% colnames(celltype_summary) && ncol(celltype_summary) >= 1L) {
    setnames(celltype_summary, colnames(celltype_summary)[[1L]], "cell_type")
  }
  if (!"n" %in% colnames(celltype_summary) && ncol(celltype_summary) >= 2L) {
    setnames(celltype_summary, colnames(celltype_summary)[[2L]], "n")
  }
  seed_text <- paste0(
    "Top seeds are dominated by ",
    paste(sprintf("%s=%s", direction_summary$dominant_direction, direction_summary$n), collapse = "; "),
    ". Recurrent contributing cell types among top seeds: ",
    if (nrow(celltype_summary)) paste(sprintf("%s=%s", head(celltype_summary$cell_type, 6L), head(celltype_summary$n, 6L)), collapse = "; ") else "NA",
    "."
  )
  list(
    direction_summary = direction_summary,
    celltype_summary = celltype_summary,
    seed_preview = preview_dt,
    seed_text = seed_text
  )
}

network_strength_label <- function(n_edges) {
  n_edges <- suppressWarnings(as.numeric(n_edges))
  if (is.na(n_edges) || n_edges <= 0) return("none")
  if (n_edges <= 2) return("very_sparse")
  if (n_edges <= 10) return("sparse")
  if (n_edges <= 40) return("moderate")
  "dense"
}

consensus_interpretation_line <- function(consensus_row, string_row, biogrid_row, hint_row) {
  pieces <- c(
    sprintf("Consensus edges=%s (%s)", consensus_row$n_edges, network_strength_label(consensus_row$n_edges)),
    sprintf("STRING edges=%s", string_row$n_edges),
    sprintf("BioGRID edges=%s", biogrid_row$n_edges),
    sprintf("HINT edges=%s", hint_row$n_edges)
  )
  paste(pieces, collapse = "; ")
}

extract_top_hubs <- function(run_root, pair_label, db_name, top_n = 12L) {
  path <- file.path(run_root, "networks", pair_label, db_name, "hub_genes.tsv")
  dt <- read_optional_tsv(path)
  if (!nrow(dt)) return(dt)
  dt[order(-degree, seed_rank)][seq_len(min(top_n, .N))]
}

extract_top_edges <- function(run_root, pair_label, db_name, top_n = 12L) {
  path <- file.path(run_root, "networks", pair_label, db_name, "edges.tsv")
  dt <- read_optional_tsv(path)
  if (!nrow(dt)) return(dt)
  if ("db_support" %in% colnames(dt)) {
    return(dt[order(-db_support, -weight, gene_a, gene_b)][seq_len(min(top_n, .N))])
  }
  dt[order(-weight, gene_a, gene_b)][seq_len(min(top_n, .N))]
}

build_pair_packet <- function(run_root,
                              pair_label,
                              db_manifest_dt,
                              master_manifest_dt,
                              network_summary_dt,
                              pair_seed_dt,
                              top_seed_n = 20L,
                              top_hub_n = 12L,
                              top_edge_n = 12L) {
	pair_value <- pair_label
	pair_network_dt <- copy(network_summary_dt[pair_label == pair_value])
  if (!nrow(pair_network_dt)) stop("No network_summary rows for pair: ", pair_label, call. = FALSE)
  pair_network_dt[, db_order := match(database, c("consensus", "STRING", "BioGRID", "HINT"))]
  setorder(pair_network_dt, db_order)
  pair_network_dt[, db_order := NULL]
  seed_context <- summarize_seed_context(pair_seed_dt, top_n = top_seed_n)
  hub_tables <- list(
    consensus = extract_top_hubs(run_root, pair_label, "consensus", top_hub_n),
    STRING = extract_top_hubs(run_root, pair_label, "STRING", top_hub_n),
    BioGRID = extract_top_hubs(run_root, pair_label, "BioGRID", top_hub_n),
    HINT = extract_top_hubs(run_root, pair_label, "HINT", top_hub_n)
  )
  edge_tables <- list(
    consensus = extract_top_edges(run_root, pair_label, "consensus", top_edge_n),
    STRING = extract_top_edges(run_root, pair_label, "STRING", top_edge_n),
    BioGRID = extract_top_edges(run_root, pair_label, "BioGRID", top_edge_n),
    HINT = extract_top_edges(run_root, pair_label, "HINT", top_edge_n)
  )

  row_or_empty <- function(db) {
    dt <- pair_network_dt[database == db]
    if (!nrow(dt)) data.table(database = db, n_edges = 0, n_isolated_seeds = NA_integer_, top_hub_gene = NA_character_, top_hub_degree = NA_real_) else dt[1L]
  }
  consensus_row <- row_or_empty("consensus")
  string_row <- row_or_empty("STRING")
  biogrid_row <- row_or_empty("BioGRID")
  hint_row <- row_or_empty("HINT")

  guardrails <- c(
    "当前 PPI 结果是把 pseudobulk DEG 编码蛋白映射到已知数据库互作网络，不是从 RNA 表达直接推断蛋白物理互作。",
    "边是否存在高度依赖数据库覆盖、基因别名映射和文献偏倚；absence of edge ≠ absence of biology。",
    "STRING/BioGRID/HINT 的 score 或 evidence_count 不是表达效应量，不能与 log2FC 直接等价比较。",
    "跨数据库 consensus 边更适合优先讨论；单库边只能作为候选支持。",
    "高 degree hub 可能反映数据库先验中心性或研究热度，不自动等于因果 driver。",
    "解释必须回到 seed gene 的表达方向、涉及的 epithelial cell type、以及后续实验/文献验证。"
  )
  hint_note <- db_manifest_dt[database == "HINT_human_binary_hq", status_note]
  if (length(hint_note) > 0L && nzchar(safe_trim(hint_note[[1L]]))) {
    guardrails <- c(guardrails, paste0("本次 HINT 接入状态：", hint_note[[1L]], "。因此 HINT 证据应视为不完整补充，而不是完整负证据。"))
  }

  packet <- list(
    packet_type = "epithelial_airway_deg_ppi_llm_pair_packet",
    packet_version = "20260531_v1",
    pair_label = pair_label,
    run_root = path_normalize(run_root),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    guardrails = guardrails,
    database_access = as.data.frame(db_manifest_dt, stringsAsFactors = FALSE),
    master_edge_manifest = as.data.frame(master_manifest_dt, stringsAsFactors = FALSE),
    pair_network_summary = as.data.frame(pair_network_dt, stringsAsFactors = FALSE),
    pair_seed_preview = as.data.frame(seed_context$seed_preview, stringsAsFactors = FALSE),
    pair_seed_direction_summary = as.data.frame(seed_context$direction_summary, stringsAsFactors = FALSE),
    pair_seed_celltype_summary = as.data.frame(seed_context$celltype_summary, stringsAsFactors = FALSE),
    seed_context_text = seed_context$seed_text,
    consensus_snapshot = consensus_interpretation_line(consensus_row, string_row, biogrid_row, hint_row),
    top_hubs = lapply(hub_tables, as.data.frame, stringsAsFactors = FALSE),
    top_edges = lapply(edge_tables, as.data.frame, stringsAsFactors = FALSE),
    source_paths = list(
      db_manifest_tsv = file.path(run_root, "db_manifest.tsv"),
      master_edge_manifest_tsv = file.path(run_root, "master_edge_manifest.tsv"),
      network_summary_tsv = file.path(run_root, "network_summary.tsv"),
      seed_summary_tsv = NA_character_,
      pair_network_dir = file.path(run_root, "networks", pair_label)
    )
  )
  packet
}

build_pair_prompt <- function(packet) {
  db_summary_df <- packet$database_access[, c("database", "access_mode", "downloaded", "status_note"), drop = FALSE]
  network_df <- packet$pair_network_summary[, c("database", "n_seed_genes", "n_edges", "n_isolated_seeds", "largest_component_nodes", "top_hub_gene", "top_hub_degree"), drop = FALSE]
  seed_df <- packet$pair_seed_preview[, c("seed_rank", "gene_symbol", "dominant_direction", "n_cell_types", "cell_types", "max_abs_log2FoldChange", "min_padj", "ppi_priority_score"), drop = FALSE]
  consensus_edges_df <- packet$top_edges$consensus
  string_edges_df <- packet$top_edges$STRING
  consensus_hubs_df <- packet$top_hubs$consensus
  string_hubs_df <- packet$top_hubs$STRING

  system_prompt <- paste(
    "You are a rigorous bioinformatics interpreter for DEG-to-PPI bridge analyses.",
    "The network is a knowledge-base overlay on DEG-encoded proteins, not a direct measurement of physical interaction in this experiment.",
    "Use cautious language in Simplified Chinese Markdown.",
    "Do not claim context-specific direct binding, confirmed rewiring, or causal driver status from network centrality alone.",
    "Prioritize cross-database consensus edges over single-database edges, and explicitly mention database coverage / prior-knowledge bias when relevant."
  )

  user_prompt <- paste(
    sprintf("请解读 PPI pair `%s` 的 DEG→PPI bridge 结果。", packet$pair_label),
    "",
    "## 必须遵守的 PPI 护栏",
    paste0("- ", packet$guardrails, collapse = "\n"),
    "",
    "## 当前 pair 的数据库网络概况",
    markdown_table(network_df, max_rows = 10L),
    "",
    "## 数据库接入说明",
    markdown_table(db_summary_df, max_rows = 10L),
    "",
    "## seed context",
    packet$seed_context_text,
    "",
    "## Top seed genes",
    markdown_table(seed_df, max_rows = 20L),
    "",
    "## Consensus top hubs",
    markdown_table(consensus_hubs_df, max_rows = 12L),
    "",
    "## STRING top hubs",
    markdown_table(string_hubs_df, max_rows = 12L),
    "",
    "## Consensus top edges",
    markdown_table(consensus_edges_df, max_rows = 12L),
    "",
    "## STRING top edges",
    markdown_table(string_edges_df, max_rows = 12L),
    "",
    "请输出：",
    "1. 这个 pair 在 PPI 层面最稳健的 candidate core 是什么；",
    "2. 哪些结论主要来自跨库 consensus，哪些只能算单库候选；",
    "3. 如何把网络结果与 seed gene 的表达方向/涉及 cell types 联系起来；",
    "4. 当前结果最重要的局限与偏差来源；",
    "5. 最值得优先文献核对或实验验证的 hub / edge / module。",
    "务必避免把 PPI 数据库网络写成当前样本里已被直接测到的物理互作。",
    sep = "\n"
  )

  list(system_prompt = system_prompt, user_prompt = user_prompt)
}

build_rule_based_interpretation <- function(packet, status = "rule_only") {
  network_dt <- as.data.table(packet$pair_network_summary)
  consensus_row <- network_dt[database == "consensus"]
  string_row <- network_dt[database == "STRING"]
  biogrid_row <- network_dt[database == "BioGRID"]
  hint_row <- network_dt[database == "HINT"]
  if (!nrow(consensus_row)) consensus_row <- data.table(n_edges = 0, n_isolated_seeds = NA_integer_, top_hub_gene = NA_character_, top_hub_degree = NA_real_)
  if (!nrow(string_row)) string_row <- data.table(n_edges = 0, n_isolated_seeds = NA_integer_, top_hub_gene = NA_character_, top_hub_degree = NA_real_)
  if (!nrow(biogrid_row)) biogrid_row <- data.table(n_edges = 0)
  if (!nrow(hint_row)) hint_row <- data.table(n_edges = 0)

  consensus_edges <- as.numeric(consensus_row$n_edges[[1L]])
  consensus_label <- switch(
    network_strength_label(consensus_edges),
    none = "几乎没有跨库共识边",
    very_sparse = "只有极少量跨库共识边",
    sparse = "跨库共识核心较小",
    moderate = "已有可讨论的跨库共识核心",
    dense = "跨库共识网络相对丰富",
    "跨库共识信号未明"
  )

  consensus_hubs_df <- as.data.table(packet$top_hubs$consensus)
  string_hubs_df <- as.data.table(packet$top_hubs$STRING)
  seed_df <- as.data.table(packet$pair_seed_preview)
  top_consensus_hubs <- if (nrow(consensus_hubs_df)) paste(sprintf("%s(deg=%s)", head(consensus_hubs_df$gene_symbol, 5L), head(consensus_hubs_df$degree, 5L)), collapse = ", ") else "NA"
  top_string_hubs <- if (nrow(string_hubs_df)) paste(sprintf("%s(deg=%s)", head(string_hubs_df$gene_symbol, 5L), head(string_hubs_df$degree, 5L)), collapse = ", ") else "NA"
  top_seed_text <- if (nrow(seed_df)) paste(sprintf("%s[%s]", head(seed_df$gene_symbol, 8L), head(seed_df$dominant_direction, 8L)), collapse = ", ") else "NA"
  celltype_summary <- as.data.table(packet$pair_seed_celltype_summary)
  top_celltypes <- if (nrow(celltype_summary)) paste(sprintf("%s=%s", head(celltype_summary$cell_type, 6L), head(celltype_summary$n, 6L)), collapse = "; ") else "NA"
  db_note <- packet$database_access[packet$database_access$database == "HINT_human_binary_hq", "status_note"]
  db_note <- if (length(db_note) > 0L) safe_trim(db_note[[1L]]) else ""

  c(
    sprintf("# PPI LLM interpretation | %s", packet$pair_label),
    "",
    sprintf("**LLM status:** `%s`", status),
    "",
    "## Rule-based overview",
    sprintf("- `%s`：STRING=%s edges，BioGRID=%s edges，HINT=%s edges，consensus=%s edges。", packet$pair_label, string_row$n_edges[[1L]], biogrid_row$n_edges[[1L]], hint_row$n_edges[[1L]], consensus_row$n_edges[[1L]]),
    sprintf("- 当前判断：%s；因此更适合把结果当作 DEG 进入已知 PPI 空间后的 `%s`，而不是当作稳定完成的 context-specific rewiring 图谱。", consensus_label, if (consensus_edges <= 2) "candidate scaffold" else "candidate interaction core"),
    sprintf("- seed 侧最值得优先联系到网络的基因包括：%s。", top_seed_text),
    sprintf("- 这些 seed 主要反复来自的 cell types：%s。", top_celltypes),
    sprintf("- cross-db 优先看的 hub：%s。", top_consensus_hubs),
    sprintf("- 若需要扩大候选范围，可再看 STRING 内部 hub：%s。", top_string_hubs),
    sprintf("- isolated seeds 仍然较多（consensus isolated=%s, STRING isolated=%s），这更像数据库覆盖/别名映射限制，而不应被解读为这些 DEG 没有生物学意义。", consensus_row$n_isolated_seeds[[1L]], string_row$n_isolated_seeds[[1L]]),
    if (nzchar(db_note)) sprintf("- HINT 额外说明：%s。", db_note) else "- HINT 本次没有额外接入告警。",
    "",
    "## PPI-specific guardrails",
    paste0("- ", packet$guardrails, collapse = "\n"),
    "",
    "## Suggested next checks",
    "- 先优先核对 consensus top edge / hub 的文献背景，而不是从单库边直接下结论。",
    "- 回到 seed gene 的表达方向和涉及 epithelial cell types，判断网络是否真的对应同一生物学状态。",
    "- 若某个 hub 主要由 STRING 支撑但跨库没有复制，应把它当作候选扩展节点，而不是主结论。",
    "- 对 degree 很高但 seed_rank 并不靠前的节点，要警惕数据库中心性偏差。"
  )
}

run_pair_llm <- function(packet,
                         output_dir,
                         llm_config = list()) {
  dir_create(output_dir)
  prompt <- build_pair_prompt(packet)
  prompt_path <- file.path(output_dir, "PPI_LLM_PROMPT.md")
  interpretation_path <- file.path(output_dir, "PPI_LLM_INTERPRETATION.md")
  raw_path <- file.path(output_dir, "PPI_LLM_RAW_RESPONSE.txt")
  status_path <- file.path(output_dir, "PPI_LLM_STATUS.json")
  packet_path <- file.path(output_dir, "PPI_LLM_PACKET.json")
  write_json(packet, packet_path)
  write_markdown(c("# PPI LLM prompt", "", "## System", prompt$system_prompt, "", "## User", prompt$user_prompt), prompt_path)

  key_info <- resolve_deepseek_key(
    api_key = llm_config$api_key,
    env_candidates = llm_config$env_candidates %||% c("/home/h2048/.env", "/home/h2048/script/.env"),
    write_env_placeholder = isTRUE(llm_config$write_env_placeholder)
  )

  status <- list(
    pair_label = packet$pair_label,
    status = "skipped_disabled",
    enabled = isTRUE(llm_config$enabled),
    model = llm_config$model,
    prompt_md = path_normalize(prompt_path),
    interpretation_md = path_normalize(interpretation_path),
    raw_response_txt = path_normalize(raw_path),
    packet_json = path_normalize(packet_path),
    env_files_loaded = key_info$env_files_loaded,
    placeholder_written = key_info$placeholder_written,
    placeholder_path = key_info$placeholder_path,
    error = NULL
  )

  rule_lines <- build_rule_based_interpretation(packet, status = if (!isTRUE(llm_config$enabled)) "skipped_disabled" else if (!isTRUE(key_info$has_live_key)) "skipped_no_live_key" else "ok")

  if (!isTRUE(llm_config$enabled)) {
    status$status <- "skipped_disabled"
    write_markdown(rule_lines, interpretation_path)
  } else if (!isTRUE(key_info$has_live_key)) {
    status$status <- "skipped_no_live_key"
    extra_note <- if (isTRUE(key_info$placeholder_written) && nzchar(safe_trim(key_info$placeholder_path))) {
      c("", "## How to enable live LLM", sprintf("已自动写入占位 `.env`：`%s`；填入真实 `DEEPSEEK_API_KEY` 后重跑即可。", key_info$placeholder_path))
    } else {
      c("", "## How to enable live LLM", "在 `/home/h2048/.env` 中提供真实 `DEEPSEEK_API_KEY` 后重跑即可。")
    }
    write_markdown(c(rule_lines, extra_note), interpretation_path)
  } else {
    response <- tryCatch(
      call_deepseek_chat(
        prompt = prompt$user_prompt,
        system_prompt = prompt$system_prompt,
        model = llm_config$model,
        api_key = key_info$api_key,
        timeout_sec = llm_config$timeout_sec,
        base_url = llm_config$base_url %||% "https://api.deepseek.com/v1/chat/completions"
      ),
      error = function(e) structure(conditionMessage(e), class = "ppi_llm_error")
    )
    if (inherits(response, "ppi_llm_error")) {
      status$status <- "error"
      status$error <- as.character(response[[1L]])
      write_markdown(c(rule_lines, "", "## LLM error", sprintf("`%s`", status$error)), interpretation_path)
    } else {
      status$status <- "ok"
      write_markdown(c(rule_lines, "", "## Online LLM interpretation", "", response), interpretation_path)
      writeLines(as.character(response), raw_path, useBytes = TRUE)
    }
  }

  write_json(status, status_path)
  status$status_json <- path_normalize(status_path)
  status
}

build_summary_markdown <- function(config, manifest_dt, db_manifest_dt) {
  lines <- c(
    "# PPI → LLM bridge summary",
    "",
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    sprintf("- Run root: `%s`", config$run_root),
    sprintf("- Output dir: `%s`", config$output_dir),
    sprintf("- Live LLM enabled: `%s`", tolower(as.character(config$enable_live))),
    sprintf("- Model: `%s`", config$model),
    "",
    "## Database access notes",
    "",
    markdown_table(db_manifest_dt[, c("database", "access_mode", "downloaded", "status_note"), drop = FALSE], max_rows = 10L),
    "",
    "## Pair interpretation manifest",
    "",
    markdown_table(manifest_dt[, c("pair_label", "llm_status", "packet_json", "prompt_md", "interpretation_md"), drop = FALSE], max_rows = 20L),
    "",
    "## PPI-specific reminder",
    "",
    "- 这些解释建立在 DEG → known-PPI mapping 上，不是样本内直接测到的 protein interaction。",
    "- 优先阅读跨库 consensus，再把单库 STRING / BioGRID / HINT 视为扩展候选。",
    "- isolated seeds、hub degree、以及库之间边数差异都必须结合数据库覆盖和文献偏倚一起解释。"
  )
  lines
}

main <- function() {
  config <- validate_config(parse_args(DEFAULTS))
  dir_create(config$output_dir)

  db_manifest_path <- file.path(config$run_root, "db_manifest.tsv")
  master_manifest_path <- file.path(config$run_root, "master_edge_manifest.tsv")
  network_summary_path <- file.path(config$run_root, "network_summary.tsv")
  seed_summary_path <- find_seed_summary_path(config$run_root)

  cat("[1/5] Reading workflow outputs...\n")
  db_manifest_dt <- read_required_tsv(db_manifest_path)
  master_manifest_dt <- read_required_tsv(master_manifest_path)
  network_summary_dt <- read_required_tsv(network_summary_path)
  seed_dt <- read_required_tsv(seed_summary_path)

  pair_labels <- unique(as.character(network_summary_dt$pair_label))
  pair_labels <- pair_labels[nzchar(pair_labels)]
  if (!is.null(config$target_pairs) && length(config$target_pairs) > 0L) {
    pair_labels <- intersect(pair_labels, config$target_pairs)
    if (!length(pair_labels)) stop("No matching pair labels found after --target-pairs filter", call. = FALSE)
  }

  cat("[2/5] Building pair packets and prompts...\n")
  manifest_rows <- list()
  for (pair_label in pair_labels) {
    cat("  - Pair:", pair_label, "\n")
    pair_dir <- file.path(config$output_dir, "pairs", pair_label)
	    pair_value <- pair_label
	    pair_seed_dt <- seed_dt[pair_label == pair_value]
    packet <- build_pair_packet(
      run_root = config$run_root,
      pair_label = pair_label,
      db_manifest_dt = db_manifest_dt,
      master_manifest_dt = master_manifest_dt,
      network_summary_dt = network_summary_dt,
      pair_seed_dt = pair_seed_dt,
      top_seed_n = config$top_seed_n,
      top_hub_n = config$top_hub_n,
      top_edge_n = config$top_edge_n
    )
    packet$source_paths$seed_summary_tsv <- seed_summary_path
    status <- run_pair_llm(
      packet = packet,
      output_dir = pair_dir,
      llm_config = list(
        enabled = config$enable_live,
        model = config$model,
        timeout_sec = config$timeout_sec,
        write_env_placeholder = config$write_env_placeholder
      )
    )
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
      pair_label = pair_label,
      llm_status = status$status,
      packet_json = status$packet_json,
      prompt_md = status$prompt_md,
      interpretation_md = status$interpretation_md,
      raw_response_txt = status$raw_response_txt,
      status_json = status$status_json,
      placeholder_written = isTRUE(status$placeholder_written),
      placeholder_path = status$placeholder_path %||% NA_character_,
      error = status$error %||% NA_character_
    )
  }

  cat("[3/5] Writing root manifest...\n")
  manifest_dt <- rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
  fwrite(manifest_dt, file.path(config$output_dir, "PPI_LLM_manifest.tsv"), sep = "\t")

  cat("[4/5] Writing bridge summary...\n")
  write_markdown(
    build_summary_markdown(config, manifest_dt, db_manifest_dt),
    file.path(config$output_dir, "PPI_LLM_BRIDGE_SUMMARY.md")
  )

  cat("[5/5] Done. Key outputs:\n")
  cat("  -", file.path(config$output_dir, "PPI_LLM_manifest.tsv"), "\n")
  cat("  -", file.path(config$output_dir, "PPI_LLM_BRIDGE_SUMMARY.md"), "\n")
  cat("  -", file.path(config$output_dir, "pairs"), "\n")
}

main()
