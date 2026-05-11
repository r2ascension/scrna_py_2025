#!/usr/bin/env Rscript
# ==============================================================================
# Tissue-comparison helper overlay: parallel LLM discovery screening
# ==============================================================================
# Date: 2026-05-08
# Purpose:
#   Source after the 2026-04 tissue-comparison helpers to override only
#   tc_run_record_screening_llm() with batch-level parallel execution.
#   This keeps the validated prompt/schema logic intact while allowing LLM
#   screening batches to run concurrently when LLM_SCREEN_PARALLEL_WORKERS > 1.
# ==============================================================================

if (!exists("tc_run_record_screening_llm_batch", mode = "function")) {
  BASE_HELPER <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260414_v2.R"
  if (!file.exists(BASE_HELPER)) {
    stop(sprintf("Base tissue-comparison helper not found: %s", BASE_HELPER))
  }
  source(BASE_HELPER)
}

.tc_parallel_llm_scalar <- function(name, default, caller_env = parent.frame()) {
  val <- tryCatch(tc_lookup_in_caller(name, NULL, caller_env), error = function(e) NULL)
  if (is.null(val) || length(val) == 0L || (length(val) == 1L && is.na(val))) {
    val <- Sys.getenv(name, unset = as.character(default))
  }
  val
}

.tc_record_screening_rows_from_batch <- function(batch_out, idx, batch_df) {
  rows <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    rows[[j]] <- data.frame(
      .row_order = idx[[j]],
      record_id = batch_out[[j]]$record_id,
      biological_signal_class = batch_out[[j]]$biological_signal_class,
      confidence = batch_out[[j]]$confidence,
      short_call = batch_out[[j]]$short_call,
      discovery_flag = batch_out[[j]]$discovery_flag,
      outlier_flag = batch_out[[j]]$outlier_flag,
      evidence_summary = batch_out[[j]]$evidence_summary,
      followup = batch_out[[j]]$followup,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

tc_run_record_screening_llm <- function(records_df,
                                        family_label,
                                        batch_size = 10L,
                                        caller_env = parent.frame()) {
  if (is.null(records_df) || !is.data.frame(records_df) || nrow(records_df) == 0) return(data.frame())

  batch_idx <- split(
    seq_len(nrow(records_df)),
    ceiling(seq_len(nrow(records_df)) / max(1L, as.integer(batch_size)))
  )
  retry_sleep <- as.numeric(.tc_parallel_llm_scalar("STANDARDIZE_LLM_RETRY_SLEEP_SEC", 2, caller_env))
  workers <- as.integer(.tc_parallel_llm_scalar("LLM_SCREEN_PARALLEL_WORKERS", 1L, caller_env))
  stagger_sec <- as.numeric(.tc_parallel_llm_scalar("LLM_SCREEN_PARALLEL_STAGGER_SEC", 0.35, caller_env))
  workers <- max(1L, min(workers, length(batch_idx)))

  run_one_batch <- function(i) {
    idx <- batch_idx[[i]]
    if (workers > 1L && stagger_sec > 0) {
      Sys.sleep(((i - 1L) %% workers) * stagger_sec)
    }
    cat(sprintf(
      "[LLM screen%s] %s batch %d/%d (%d records)\n",
      if (workers > 1L) sprintf(" parallel:%d", workers) else "",
      family_label,
      i,
      length(batch_idx),
      length(idx)
    ))
    batch_out <- tc_run_record_screening_llm_batch(
      records_df[idx, , drop = FALSE],
      family_label = family_label,
      caller_env = caller_env
    )
    .tc_record_screening_rows_from_batch(batch_out, idx, records_df[idx, , drop = FALSE])
  }

  if (workers <= 1L || .Platform$OS.type == "windows") {
    rows <- list()
    for (i in seq_along(batch_idx)) {
      rows[[length(rows) + 1L]] <- run_one_batch(i)
      Sys.sleep(retry_sleep)
    }
  } else {
    rows <- parallel::mclapply(
      seq_along(batch_idx),
      run_one_batch,
      mc.cores = workers,
      mc.preschedule = FALSE
    )
  }

  out <- dplyr::bind_rows(rows)
  out <- out[order(out$.row_order), , drop = FALSE]
  out$.row_order <- NULL
  dplyr::left_join(records_df, out, by = "record_id")
}

cat("[parallel-llm] tc_run_record_screening_llm override loaded (set LLM_SCREEN_PARALLEL_WORKERS > 1 to enable).\n")
