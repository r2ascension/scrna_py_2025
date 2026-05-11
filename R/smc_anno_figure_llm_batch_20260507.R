#!/usr/bin/env Rscript
# ==============================================================================
# smc_anno figure/table LLM batch runner (2026-05-07)
# ------------------------------------------------------------------------------
# Reads queued companion manifest rows produced by smcanno_viz_register_figure(),
# calls the online LLM outside the plotting/analysis process, and updates the
# paired *_LLM_ANALYSIS.md, *_LLM_STATUS.json, and manifest status columns.
# ==============================================================================

HELPER_PATH <- "/home/h2048/script/R/smc_anno_visual_llm_helper_20260506.R"
if (!file.exists(HELPER_PATH)) stop(sprintf("Missing helper: %s", HELPER_PATH), call. = FALSE)
source(HELPER_PATH)

smcanno_batch_usage <- function() {
  paste(c(
    "Usage:",
    "  Rscript smc_anno_figure_llm_batch_20260507.R --output-dir <dir>",
    "  Rscript smc_anno_figure_llm_batch_20260507.R --manifest <manifest.tsv> [--manifest <manifest2.tsv>]",
    "",
    "Options:",
    "  --output-dir <dir>       Recursively find *_figure_companion_manifest.tsv / *_llm_manifest.tsv under dir.",
    "  --manifest <path>        Manifest to process; may be repeated.",
    "  --model <name>           Override model, default uses row/status model then deepseek-reasoner.",
    "  --timeout-sec <seconds>  Per-record LLM timeout; default 180.",
    "  --max-items <n>          Process at most n queued rows.",
    "  --statuses <csv>         Statuses to process; default queued_for_batch,error.",
    "  --force                  Re-run rows even when status is ok.",
    "  --dry-run                Print queued rows without calling the LLM."
  ), collapse = "\n")
}

smcanno_batch_parse_args <- function(args) {
  out <- list(manifests = character(), output_dirs = character(), model = "", timeout_sec = 180,
              max_items = Inf, statuses = c("queued_for_batch", "error"), force = FALSE, dry_run = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--manifest")) {
      i <- i + 1L; out$manifests <- c(out$manifests, args[[i]])
    } else if (identical(arg, "--output-dir")) {
      i <- i + 1L; out$output_dirs <- c(out$output_dirs, args[[i]])
    } else if (identical(arg, "--model")) {
      i <- i + 1L; out$model <- args[[i]]
    } else if (identical(arg, "--timeout-sec")) {
      i <- i + 1L; out$timeout_sec <- suppressWarnings(as.numeric(args[[i]]))
    } else if (identical(arg, "--max-items")) {
      i <- i + 1L; out$max_items <- suppressWarnings(as.integer(args[[i]]))
    } else if (identical(arg, "--statuses")) {
      i <- i + 1L; out$statuses <- trimws(strsplit(args[[i]], ",", fixed = TRUE)[[1]])
      out$statuses <- out$statuses[nzchar(out$statuses)]
    } else if (identical(arg, "--force")) {
      out$force <- TRUE
    } else if (identical(arg, "--dry-run")) {
      out$dry_run <- TRUE
    } else if (identical(arg, "--help") || identical(arg, "-h")) {
      cat(smcanno_batch_usage(), "\n")
      quit(status = 0)
    } else {
      stop(sprintf("Unknown argument: %s\n\n%s", arg, smcanno_batch_usage()), call. = FALSE)
    }
    i <- i + 1L
  }
  if (!is.finite(out$timeout_sec) || out$timeout_sec <= 0) out$timeout_sec <- 180
  if (!is.finite(out$max_items) || out$max_items <= 0) out$max_items <- Inf
  out
}

smcanno_batch_find_manifests <- function(output_dirs) {
  unique(unlist(lapply(output_dirs, function(root) {
    if (!dir.exists(root)) return(character())
    files <- list.files(root, pattern = "(figure_companion_manifest|llm_manifest)\\.tsv$", recursive = TRUE, full.names = TRUE)
    files[file.exists(files)]
  }), use.names = FALSE))
}

smcanno_batch_update_manifest_row <- function(manifest_df, idx, status, llm_mode = "online_batch") {
  for (nm in names(status)) {
    if (length(status[[nm]]) != 1L || is.list(status[[nm]])) next
    col <- if (identical(nm, "status")) "llm_status" else nm
    if (!col %in% colnames(manifest_df)) manifest_df[[col]] <- ""
    manifest_df[[col]][idx] <- as.character(status[[nm]])
  }
  if (!"llm_mode" %in% colnames(manifest_df)) manifest_df$llm_mode <- ""
  manifest_df$llm_mode[idx] <- llm_mode
  manifest_df
}

smcanno_batch_write_analysis <- function(status, data_df, response = NULL, error = NULL) {
  method <- smcanno_viz_to_scalar(status$method)
  figure_type <- smcanno_viz_to_scalar(status$figure_type)
  title <- smcanno_viz_to_scalar(status$title)
  if (!nzchar(method)) method <- "unknown"
  if (!nzchar(figure_type)) figure_type <- "unknown"
  if (!nzchar(title)) title <- basename(smcanno_viz_to_scalar(status$analysis_md))
  rule_lines <- smcanno_viz_default_rule_analysis(method, figure_type, title, data_df, extra_context = NULL)
  final_lines <- smcanno_viz_final_summary_lines(method, figure_type, title, data_df, extra_context = NULL)
  if (!is.null(response) && nzchar(smcanno_viz_to_scalar(response))) {
    return(c("# Figure LLM analysis", "", "**LLM status:** `ok`", "", "## 数据驱动审计", "", rule_lines, "", "## 在线 LLM 深度解读", "", as.character(response), "", final_lines))
  }
  if (!is.null(error) && nzchar(smcanno_viz_to_scalar(error))) {
    return(c("# Figure LLM analysis", "", "**LLM status:** `error`", "", sprintf("LLM error: `%s`", smcanno_viz_to_scalar(error)), "", rule_lines, "", final_lines))
  }
  c("# Figure LLM analysis", "", "**LLM status:** `skipped_no_live_key`", "", "未检测到真实 `DEEPSEEK_API_KEY`，保留规则解读；稍后可重新运行 batch runner。", "", rule_lines, "", final_lines)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

smcanno_batch_process_row <- function(row, timeout_sec, model_override = "", dry_run = FALSE) {
  status_path <- smcanno_viz_to_scalar(row$llm_status_json)
  prompt_path <- smcanno_viz_to_scalar(row$llm_prompt_md)
  analysis_path <- smcanno_viz_to_scalar(row$llm_analysis_md)
  data_path <- smcanno_viz_to_scalar(row$data_csv)
  status <- smcanno_viz_read_json(status_path)
  if (length(status) == 0L) status <- as.list(row)
  status$status <- smcanno_viz_to_scalar(smcanno_viz_null_coalesce(status$status, row$llm_status))
  status$prompt_md <- prompt_path
  status$analysis_md <- analysis_path
  status$data_csv <- data_path
  status$runner_script <- "/home/h2048/script/R/smc_anno_figure_llm_batch_20260507.R"
  status$started_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  if (isTRUE(dry_run)) {
    cat(sprintf("[DRY-RUN] %s | %s | %s\n", smcanno_viz_to_scalar(status$method), smcanno_viz_to_scalar(status$figure_type), data_path))
    status$status <- "dry_run"
    return(status)
  }

  if (!file.exists(prompt_path)) stop(sprintf("Missing prompt: %s", prompt_path), call. = FALSE)
  if (!file.exists(data_path)) stop(sprintf("Missing data CSV: %s", data_path), call. = FALSE)
  data_df <- utils::read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
  key_info <- smcanno_viz_resolve_deepseek_key()
  status$env_files_loaded <- key_info$env_files_loaded
  if (!isTRUE(key_info$has_live_key)) {
    status$status <- "skipped_no_live_key"
    status$error <- NULL
    status$completed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    smcanno_viz_write_markdown(smcanno_batch_write_analysis(status, data_df), analysis_path)
    smcanno_viz_write_json(status, status_path)
    return(status)
  }

  prompt <- smcanno_viz_read_prompt_sections(prompt_path)
  model <- smcanno_viz_to_scalar(model_override)
  if (!nzchar(model)) model <- smcanno_viz_to_scalar(status$model)
  if (!nzchar(model)) model <- "deepseek-reasoner"
  status$model <- model
  t0 <- proc.time()[["elapsed"]]
  response <- tryCatch(
    smcanno_viz_deepseek_chat(
      prompt = prompt$user_prompt,
      system_prompt = prompt$system_prompt,
      model = model,
      api_key = key_info$api_key,
      timeout_sec = timeout_sec
    ),
    error = function(e) structure(conditionMessage(e), class = "smcanno_batch_llm_error")
  )
  status$elapsed_sec <- round(proc.time()[["elapsed"]] - t0, 3)
  status$completed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (inherits(response, "smcanno_batch_llm_error")) {
    status$status <- "error"
    status$error <- as.character(response[[1]])
    smcanno_viz_write_markdown(smcanno_batch_write_analysis(status, data_df, error = status$error), analysis_path)
  } else {
    status$status <- "ok"
    status$error <- NULL
    smcanno_viz_write_markdown(smcanno_batch_write_analysis(status, data_df, response = response), analysis_path)
  }
  smcanno_viz_write_json(status, status_path)
  status
}

args <- smcanno_batch_parse_args(commandArgs(trailingOnly = TRUE))
manifest_paths <- unique(c(args$manifests, smcanno_batch_find_manifests(args$output_dirs)))
manifest_paths <- normalizePath(manifest_paths[file.exists(manifest_paths)], winslash = "/", mustWork = FALSE)
if (length(manifest_paths) == 0L) stop(sprintf("No manifests found.\n\n%s", smcanno_batch_usage()), call. = FALSE)

cat(sprintf("[LLM batch] started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("[LLM batch] manifests: %s\n", length(manifest_paths)))
processed <- 0L
for (manifest_path in manifest_paths) {
  cat(sprintf("[LLM batch] reading manifest: %s\n", manifest_path))
  manifest <- utils::read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("data_csv", "llm_prompt_md", "llm_analysis_md", "llm_status_json")
  if (!all(required %in% colnames(manifest))) {
    cat(sprintf("[LLM batch] skip manifest without required columns: %s\n", manifest_path))
    next
  }
  if (!"llm_status" %in% colnames(manifest)) manifest$llm_status <- ""
  row_status <- manifest$llm_status
  todo <- seq_len(nrow(manifest))[isTRUE(args$force) | row_status %in% args$statuses]
  if (length(todo) == 0L) {
    cat("[LLM batch] no queued rows in this manifest.\n")
    next
  }
  for (idx in todo) {
    if (processed >= args$max_items) break
    row <- manifest[idx, , drop = FALSE]
    label <- paste(c(smcanno_viz_to_scalar(row$method), smcanno_viz_to_scalar(row$figure_type), basename(smcanno_viz_to_scalar(row$data_csv))), collapse = " | ")
    cat(sprintf("[LLM batch] processing %s\n", label))
    status <- tryCatch(
      smcanno_batch_process_row(row, timeout_sec = args$timeout_sec, model_override = args$model, dry_run = args$dry_run),
      error = function(e) {
        msg <- conditionMessage(e)
        cat(sprintf("[LLM batch] row error: %s\n", msg))
        list(status = "error", error = msg)
      }
    )
    if (!isTRUE(args$dry_run)) {
      manifest <- smcanno_batch_update_manifest_row(manifest, idx, status, llm_mode = "online_batch")
    }
    processed <- processed + 1L
    cat(sprintf("[LLM batch] status=%s\n", smcanno_viz_to_scalar(status$status)))
  }
  if (!isTRUE(args$dry_run)) {
    utils::write.table(manifest, file = manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  }
  if (processed >= args$max_items) break
}
cat(sprintf("[LLM batch] finished: %s | processed=%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), processed))
