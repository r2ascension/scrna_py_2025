#!/usr/bin/env Rscript

UNIT_LLM_SCRIPT <- "/home/h2048/script/R/launch_program_unit_llm_20260513.R"
CNMF_L2_SCRIPT <- "/home/h2048/script/R/launch_cnmf_l2_unit_llm_20260521.R"
CNMF_L3_SCRIPT <- "/home/h2048/script/R/launch_cnmf_l3_unit_llm_20260519.R"

source(UNIT_LLM_SCRIPT)

RUN_ROOT <- Sys.getenv("CROSS_METHOD_LLM_RUN_ROOT", unset = "/home/h2048/output/program_full_parallel_methods_20260507")
OUTPUT_SUBDIR <- Sys.getenv("CROSS_METHOD_LLM_OUTPUT_SUBDIR", unset = "llm_cross_method")

split_env_cross_method <- function(name, default = character()) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  out <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

load_isolated_script_env_cross_method <- function(script_path) {
  env <- new.env(parent = globalenv())
  sys.source(script_path, envir = env)
  env
}

cnmf_launcher_envs_cross_method <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- list(
        l2 = load_isolated_script_env_cross_method(CNMF_L2_SCRIPT),
        l3 = load_isolated_script_env_cross_method(CNMF_L3_SCRIPT)
      )
    }
    cache
  }
})

empty_df_cross_method <- function(cols) {
  out <- as.data.frame(setNames(vector("list", length(cols)), cols), stringsAsFactors = FALSE, check.names = FALSE)
  out[0, , drop = FALSE]
}

normalize_scalar_cross_method <- function(x) {
  value <- sanitize_text_unit_llm(x)[1]
  if (is.na(value)) return(NA_character_)
  value <- trimws(value)
  if (!nzchar(value)) NA_character_ else value
}

normalize_vector_cross_method <- function(x) {
  value <- sanitize_text_unit_llm(x)
  value <- trimws(value)
  value[!nzchar(value)] <- NA_character_
  value
}

safe_celltype_cross_method <- function(x) {
  value <- normalize_scalar_cross_method(x)
  if (is.na(value)) NA_character_ else safe_unit_id(value)
}

logical_yes_no_cross_method <- function(x) if (isTRUE(x)) "yes" else "no"

status_from_counts_cross_method <- function(ok_units, error_units, total_units) {
  ok_units <- suppressWarnings(as.integer(ok_units)[1])
  error_units <- suppressWarnings(as.integer(error_units)[1])
  total_units <- suppressWarnings(as.integer(total_units)[1])
  if (!is.na(ok_units) && ok_units > 0L) return("ok")
  if (!is.na(error_units) && error_units > 0L) return("error_only")
  if (!is.na(total_units) && total_units > 0L) return("present_no_ok")
  "missing"
}

match_celltype_cross_method <- function(value, safe_value, target_value, target_safe) {
  value <- normalize_scalar_cross_method(value)
  safe_value <- normalize_scalar_cross_method(safe_value)
  target_value <- normalize_scalar_cross_method(target_value)
  target_safe <- normalize_scalar_cross_method(target_safe)
  if (!is.na(value) && !is.na(target_value) && identical(value, target_value)) return(TRUE)
  if (!is.na(safe_value) && !is.na(target_safe) && identical(safe_value, target_safe)) return(TRUE)
  FALSE
}

filter_records_by_lineage_cross_method <- function(df, lineage) {
  if (!is.data.frame(df) || nrow(df) == 0L) return(df)
  df[normalize_vector_cross_method(df$lineage) == normalize_scalar_cross_method(lineage), , drop = FALSE]
}

filter_records_by_celltype_cross_method <- function(df, field, target_value, target_safe) {
  if (!is.data.frame(df) || nrow(df) == 0L) return(df)
  values <- df[[field]]
  keep <- vapply(seq_len(nrow(df)), function(i) {
    match_celltype_cross_method(values[[i]], df$safe_celltype[[i]], target_value = target_value, target_safe = target_safe)
  }, logical(1))
  df[keep, , drop = FALSE]
}

combined_group_records_cross_method <- function(combined, level_tag, celltype_field) {
  cols <- c(
    "lineage", "level_tag", "celltype_field", "celltype_value", "safe_celltype",
    "celltype_l1", "celltype_l2", "celltype_l3",
    "total_units", "ok_units", "error_units", "group_status",
    "combined_md", "combined_manifest_tsv", "combined_manifest_json"
  )
  if (!is.list(combined) || length(combined$group_results) == 0L) return(empty_df_cross_method(cols))
  rows <- lapply(names(combined$group_results), function(group_name) {
    group <- combined$group_results[[group_name]]
    manifest_df <- group$manifest_df
    if (!is.data.frame(manifest_df) || nrow(manifest_df) == 0L) return(NULL)
    first <- manifest_df[1, , drop = FALSE]
    total_units <- nrow(manifest_df)
    ok_units <- sum(manifest_df$status == "ok", na.rm = TRUE)
    error_units <- sum(manifest_df$status == "error", na.rm = TRUE)
    data.frame(
      lineage = normalize_scalar_cross_method(first$lineage[[1]]),
      level_tag = level_tag,
      celltype_field = celltype_field,
      celltype_value = normalize_scalar_cross_method(first[[celltype_field]][[1]]),
      safe_celltype = normalize_scalar_cross_method(first$safe_celltype[[1]]),
      celltype_l1 = normalize_scalar_cross_method(first$celltype_l1[[1]]),
      celltype_l2 = normalize_scalar_cross_method(first$celltype_l2[[1]]),
      celltype_l3 = normalize_scalar_cross_method(first$celltype_l3[[1]]),
      total_units = total_units,
      ok_units = ok_units,
      error_units = error_units,
      group_status = status_from_counts_cross_method(ok_units, error_units, total_units),
      combined_md = normalizePath(group$combined_md, winslash = "/", mustWork = FALSE),
      combined_manifest_tsv = normalizePath(group$combined_manifest_tsv, winslash = "/", mustWork = FALSE),
      combined_manifest_json = normalizePath(group$combined_manifest_json, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(empty_df_cross_method(cols))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

collect_cnmf_combined_context_cross_method <- function(run_root = RUN_ROOT,
                                                       lineages = character(),
                                                       celltypes = character(),
                                                       all_cnmf_k = TRUE) {
  envs <- cnmf_launcher_envs_cross_method()
  l2_idx <- envs$l2$build_l2_cnmf_unit_rows(run_root = run_root, lineages = lineages, celltypes = celltypes, all_cnmf_k = all_cnmf_k)
  l2_combined <- combine_cnmf_unit_llm_by_celltype(
    l2_idx,
    run_root = run_root,
    output_subdir = envs$l2$OUTPUT_SUBDIR,
    level_tag = "l2",
    celltype_field = "celltype_l2"
  )
  l3_idx <- envs$l3$build_l3_cnmf_unit_rows(run_root = run_root, lineages = lineages, celltypes = celltypes, all_cnmf_k = all_cnmf_k)
  l3_combined <- combine_cnmf_unit_llm_by_celltype(
    l3_idx,
    run_root = run_root,
    output_subdir = envs$l3$OUTPUT_SUBDIR,
    level_tag = "l3",
    celltype_field = "celltype_l3"
  )
  list(
    l2_index = l2_idx,
    l3_index = l3_idx,
    l2_combined = l2_combined,
    l3_combined = l3_combined,
    l2_records = combined_group_records_cross_method(l2_combined, level_tag = "l2", celltype_field = "celltype_l2"),
    l3_records = combined_group_records_cross_method(l3_combined, level_tag = "l3", celltype_field = "celltype_l3")
  )
}

collect_root_unit_records_cross_method <- function(run_root = RUN_ROOT,
                                                   lineages = character(),
                                                   source_methods = c("covarnet", "hdwgcna")) {
  idx <- build_unit_llm_task_index(run_root = run_root, lineages = lineages, source_methods = source_methods, all_cnmf_k = TRUE)
  records <- unit_llm_index_records(idx)
  if (!is.data.frame(records) || nrow(records) == 0L) return(records)
  records <- records[records$source_method %in% source_methods, , drop = FALSE]
  records <- records[!is.na(records$celltype) & nzchar(records$celltype), , drop = FALSE]
  records$report_md <- records$interpretation_md
  records
}

collect_global_method_reports_cross_method <- function(run_root = RUN_ROOT,
                                                       lineages = character(),
                                                       source_methods = c("pycogaps")) {
  cols <- c("lineage", "source_method", "report_md", "status_json", "status", "input_index_tsv")
  if (length(lineages) == 0L) {
    dirs <- list.dirs(run_root, recursive = FALSE, full.names = FALSE)
    lineages <- sort(dirs[nzchar(dirs)])
  }
  rows <- list()
  for (lineage in lineages) {
    lineage_dir <- file.path(run_root, lineage, "llm_parallel")
    for (method in source_methods) {
      method_dir <- file.path(lineage_dir, method)
      if (!dir.exists(method_dir)) next
      report_md <- file.path(method_dir, sprintf("%s_LLM_interpretation.md", method))
      status_json <- file.path(method_dir, sprintf("%s_LLM_status.json", method))
      input_index_tsv <- file.path(method_dir, sprintf("%s_llm_input_index.tsv", method))
      if (!file.exists(report_md) && !file.exists(status_json)) next
      status <- read_json_unit_llm(status_json, default = list())
      status_value <- normalize_scalar_cross_method(status$status %||% if (file.exists(report_md)) "present_no_status" else NA_character_)
      rows[[length(rows) + 1L]] <- data.frame(
        lineage = normalize_scalar_cross_method(lineage),
        source_method = method,
        report_md = normalizePath(report_md, winslash = "/", mustWork = FALSE),
        status_json = normalizePath(status_json, winslash = "/", mustWork = FALSE),
        status = status_value,
        input_index_tsv = normalizePath(input_index_tsv, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(empty_df_cross_method(cols))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

build_cross_method_targets <- function(l2_records, l3_records, root_unit_records) {
  cols <- c("lineage", "target_level", "target_celltype", "safe_celltype")
  exact_rows <- list()
  if (is.data.frame(l3_records) && nrow(l3_records) > 0L) {
    exact_rows[[length(exact_rows) + 1L]] <- data.frame(
      lineage = l3_records$lineage,
      target_level = "exact",
      target_celltype = l3_records$celltype_l3,
      safe_celltype = l3_records$safe_celltype,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (is.data.frame(root_unit_records) && nrow(root_unit_records) > 0L) {
    exact_rows[[length(exact_rows) + 1L]] <- data.frame(
      lineage = root_unit_records$lineage,
      target_level = "exact",
      target_celltype = root_unit_records$celltype,
      safe_celltype = root_unit_records$safe_celltype,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  exact_df <- if (length(exact_rows) == 0L) empty_df_cross_method(cols) else do.call(rbind, exact_rows)
  if (nrow(exact_df) > 0L) {
    exact_df <- exact_df[!is.na(exact_df$target_celltype) & nzchar(exact_df$target_celltype), , drop = FALSE]
    exact_df <- exact_df[!duplicated(paste(exact_df$lineage, exact_df$safe_celltype, sep = "||")), , drop = FALSE]
  }

  l2_only_df <- empty_df_cross_method(cols)
  if (is.data.frame(l2_records) && nrow(l2_records) > 0L) {
    l2_only_df <- data.frame(
      lineage = l2_records$lineage,
      target_level = "l2_only",
      target_celltype = l2_records$celltype_l2,
      safe_celltype = l2_records$safe_celltype,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    l2_only_df <- l2_only_df[!is.na(l2_only_df$target_celltype) & nzchar(l2_only_df$target_celltype), , drop = FALSE]
    if (nrow(exact_df) > 0L) {
      exact_key <- paste(exact_df$lineage, exact_df$safe_celltype, sep = "||")
      l2_only_df <- l2_only_df[!paste(l2_only_df$lineage, l2_only_df$safe_celltype, sep = "||") %in% exact_key, , drop = FALSE]
    }
    l2_only_df <- l2_only_df[!duplicated(paste(l2_only_df$lineage, l2_only_df$safe_celltype, sep = "||")), , drop = FALSE]
  }
  out <- rbind(exact_df, l2_only_df)
  if (nrow(out) == 0L) return(out)
  out <- out[order(out$lineage, out$target_level, out$target_celltype), , drop = FALSE]
  rownames(out) <- NULL
  out
}

new_dossier_manifest_row_cross_method <- function(lineage,
                                                  target_level,
                                                  target_celltype,
                                                  safe_celltype,
                                                  section_scope,
                                                  source_method,
                                                  source_kind,
                                                  item_id,
                                                  item_label,
                                                  status,
                                                  report_md,
                                                  manifest_tsv = NA_character_,
                                                  manifest_json = NA_character_,
                                                  status_json = NA_character_,
                                                  inline_in_dossier = TRUE,
                                                  note = NA_character_) {
  data.frame(
    lineage = normalize_scalar_cross_method(lineage),
    target_level = normalize_scalar_cross_method(target_level),
    target_celltype = normalize_scalar_cross_method(target_celltype),
    safe_celltype = normalize_scalar_cross_method(safe_celltype),
    section_scope = normalize_scalar_cross_method(section_scope),
    source_method = normalize_scalar_cross_method(source_method),
    source_kind = normalize_scalar_cross_method(source_kind),
    item_id = normalize_scalar_cross_method(item_id),
    item_label = normalize_scalar_cross_method(item_label),
    status = normalize_scalar_cross_method(status),
    report_md = normalize_scalar_cross_method(report_md),
    manifest_tsv = normalize_scalar_cross_method(manifest_tsv),
    manifest_json = normalize_scalar_cross_method(manifest_json),
    status_json = normalize_scalar_cross_method(status_json),
    inline_in_dossier = logical_yes_no_cross_method(inline_in_dossier),
    note = normalize_scalar_cross_method(note),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

append_report_section_cross_method <- function(lines,
                                               title,
                                               report_md,
                                               meta_lines = character()) {
  body <- read_unit_llm_markdown_body(report_md)
  c(lines, title, "", meta_lines, if (length(meta_lines) > 0L) "" else NULL, body, "")
}

write_cross_method_celltype_dossiers <- function(run_root = RUN_ROOT,
                                                 l2_records,
                                                 l3_records,
                                                 root_unit_records,
                                                 global_reports,
                                                 output_subdir = OUTPUT_SUBDIR) {
  targets <- build_cross_method_targets(l2_records, l3_records, root_unit_records)
  empty_result <- list(
    master_manifest_df = empty_df_cross_method(c(
      "lineage", "target_level", "target_celltype", "safe_celltype", "parent_l2",
      "child_l3_count", "has_cnmf_l2", "has_cnmf_l3", "covarnet_units", "hdwgcna_units",
      "has_pycogaps_global", "dossier_md", "dossier_manifest_tsv", "dossier_manifest_json"
    )),
    master_manifest_tsv = NA_character_,
    master_manifest_json = NA_character_,
    dossier_results = list()
  )
  if (!is.data.frame(targets) || nrow(targets) == 0L) return(empty_result)

  master_rows <- list()
  dossier_results <- list()
  for (i in seq_len(nrow(targets))) {
    target <- targets[i, , drop = FALSE]
    lineage <- target$lineage[[1]]
    target_level <- target$target_level[[1]]
    target_celltype <- target$target_celltype[[1]]
    target_safe <- target$safe_celltype[[1]]

    lineage_l2 <- filter_records_by_lineage_cross_method(l2_records, lineage)
    lineage_l3 <- filter_records_by_lineage_cross_method(l3_records, lineage)
    lineage_root <- filter_records_by_lineage_cross_method(root_unit_records, lineage)
    lineage_global <- filter_records_by_lineage_cross_method(global_reports, lineage)

    l3_exact <- filter_records_by_celltype_cross_method(lineage_l3, "celltype_l3", target_celltype, target_safe)
    l2_exact <- filter_records_by_celltype_cross_method(lineage_l2, "celltype_l2", target_celltype, target_safe)
    covarnet_exact <- lineage_root[lineage_root$source_method == "covarnet", , drop = FALSE]
    covarnet_exact <- filter_records_by_celltype_cross_method(covarnet_exact, "celltype", target_celltype, target_safe)
    hdwgcna_exact <- lineage_root[lineage_root$source_method == "hdwgcna", , drop = FALSE]
    hdwgcna_exact <- filter_records_by_celltype_cross_method(hdwgcna_exact, "celltype", target_celltype, target_safe)

    parent_l2_value <- NA_character_
    if (nrow(l3_exact) > 0L) {
      parent_l2_value <- normalize_scalar_cross_method(l3_exact$celltype_l2[[1]])
    } else if (nrow(l2_exact) > 0L) {
      parent_l2_value <- normalize_scalar_cross_method(l2_exact$celltype_l2[[1]])
    }
    parent_l2_safe <- safe_celltype_cross_method(parent_l2_value)
    l2_context <- if (!is.na(parent_l2_value)) {
      filter_records_by_celltype_cross_method(lineage_l2, "celltype_l2", parent_l2_value, parent_l2_safe)
    } else {
      empty_df_cross_method(colnames(lineage_l2))
    }
    if (target_level == "l2_only" && nrow(l2_exact) > 0L) l2_context <- l2_exact

    child_l3 <- empty_df_cross_method(colnames(lineage_l3))
    if (nrow(lineage_l3) > 0L) {
      child_l3 <- lineage_l3[normalize_vector_cross_method(lineage_l3$celltype_l2) == normalize_scalar_cross_method(target_celltype), , drop = FALSE]
      if (nrow(child_l3) > 0L) child_l3 <- child_l3[order(child_l3$celltype_l3), , drop = FALSE]
    }

    pycogaps_global <- lineage_global[lineage_global$source_method == "pycogaps", , drop = FALSE]

    target_dir <- file.path(run_root, lineage, output_subdir, "celltypes", target_safe)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    prefix <- sprintf("cross_method_%s", target_safe)
    dossier_md <- file.path(target_dir, sprintf("%s_LLM_dossier.md", prefix))
    dossier_manifest_tsv <- file.path(target_dir, sprintf("%s_LLM_dossier_manifest.tsv", prefix))
    dossier_manifest_json <- file.path(target_dir, sprintf("%s_LLM_dossier_manifest.json", prefix))

    manifest_rows <- list()
    if (nrow(l3_exact) > 0L) {
      for (j in seq_len(nrow(l3_exact))) {
        row <- l3_exact[j, , drop = FALSE]
        manifest_rows[[length(manifest_rows) + 1L]] <- new_dossier_manifest_row_cross_method(
          lineage = lineage,
          target_level = target_level,
          target_celltype = target_celltype,
          safe_celltype = target_safe,
          section_scope = "cnmf_l3_exact",
          source_method = "cnmf",
          source_kind = "combined_bundle",
          item_id = sprintf("cnmf_l3_%s", row$safe_celltype[[1]]),
          item_label = sprintf("cNMF L3 bundle: %s", row$celltype_l3[[1]]),
          status = row$group_status[[1]],
          report_md = row$combined_md[[1]],
          manifest_tsv = row$combined_manifest_tsv[[1]],
          manifest_json = row$combined_manifest_json[[1]],
          status_json = NA_character_,
          inline_in_dossier = TRUE,
          note = "exact L3 celltype"
        )
      }
    }
    if (nrow(l2_context) > 0L) {
      for (j in seq_len(nrow(l2_context))) {
        row <- l2_context[j, , drop = FALSE]
        l2_note <- if (target_level == "l2_only") "exact L2 celltype" else "parent L2 context"
        manifest_rows[[length(manifest_rows) + 1L]] <- new_dossier_manifest_row_cross_method(
          lineage = lineage,
          target_level = target_level,
          target_celltype = target_celltype,
          safe_celltype = target_safe,
          section_scope = if (target_level == "l2_only") "cnmf_l2_exact" else "cnmf_l2_parent",
          source_method = "cnmf",
          source_kind = "combined_bundle",
          item_id = sprintf("cnmf_l2_%s", row$safe_celltype[[1]]),
          item_label = sprintf("cNMF L2 bundle: %s", row$celltype_l2[[1]]),
          status = row$group_status[[1]],
          report_md = row$combined_md[[1]],
          manifest_tsv = row$combined_manifest_tsv[[1]],
          manifest_json = row$combined_manifest_json[[1]],
          status_json = NA_character_,
          inline_in_dossier = TRUE,
          note = l2_note
        )
      }
    }
    if (nrow(covarnet_exact) > 0L) {
      for (j in seq_len(nrow(covarnet_exact))) {
        row <- covarnet_exact[j, , drop = FALSE]
        manifest_rows[[length(manifest_rows) + 1L]] <- new_dossier_manifest_row_cross_method(
          lineage = lineage,
          target_level = target_level,
          target_celltype = target_celltype,
          safe_celltype = target_safe,
          section_scope = "covarnet_exact",
          source_method = "covarnet",
          source_kind = "unit",
          item_id = row$unit_id[[1]],
          item_label = row$unit_label[[1]],
          status = row$status[[1]],
          report_md = row$report_md[[1]],
          status_json = row$status_json[[1]],
          inline_in_dossier = TRUE,
          note = "exact celltype covariation profile"
        )
      }
    }
    if (nrow(hdwgcna_exact) > 0L) {
      for (j in seq_len(nrow(hdwgcna_exact))) {
        row <- hdwgcna_exact[j, , drop = FALSE]
        manifest_rows[[length(manifest_rows) + 1L]] <- new_dossier_manifest_row_cross_method(
          lineage = lineage,
          target_level = target_level,
          target_celltype = target_celltype,
          safe_celltype = target_safe,
          section_scope = "hdwgcna_exact",
          source_method = "hdwgcna",
          source_kind = "unit",
          item_id = row$unit_id[[1]],
          item_label = row$unit_label[[1]],
          status = row$status[[1]],
          report_md = row$report_md[[1]],
          status_json = row$status_json[[1]],
          inline_in_dossier = TRUE,
          note = sprintf("module=%s", normalize_scalar_cross_method(row$module[[1]]))
        )
      }
    }
    if (nrow(pycogaps_global) > 0L) {
      for (j in seq_len(nrow(pycogaps_global))) {
        row <- pycogaps_global[j, , drop = FALSE]
        manifest_rows[[length(manifest_rows) + 1L]] <- new_dossier_manifest_row_cross_method(
          lineage = lineage,
          target_level = target_level,
          target_celltype = target_celltype,
          safe_celltype = target_safe,
          section_scope = "lineage_global_context",
          source_method = "pycogaps",
          source_kind = "method_summary",
          item_id = sprintf("%s_global_summary", row$source_method[[1]]),
          item_label = sprintf("%s lineage-global summary", row$source_method[[1]]),
          status = row$status[[1]],
          report_md = row$report_md[[1]],
          status_json = row$status_json[[1]],
          inline_in_dossier = FALSE,
          note = "lineage-global context only; not inlined to avoid repetition"
        )
      }
    }

    manifest_df <- if (length(manifest_rows) == 0L) empty_df_cross_method(c(
      "lineage", "target_level", "target_celltype", "safe_celltype", "section_scope",
      "source_method", "source_kind", "item_id", "item_label", "status", "report_md",
      "manifest_tsv", "manifest_json", "status_json", "inline_in_dossier", "note"
    )) else do.call(rbind, manifest_rows)
    write_tsv_unit_llm(manifest_df, dossier_manifest_tsv)
    write_json_unit_llm(list(
      lineage = lineage,
      target_level = target_level,
      target_celltype = target_celltype,
      safe_celltype = target_safe,
      parent_l2 = parent_l2_value,
      child_l3 = if (nrow(child_l3) > 0L) unique(as.character(child_l3$celltype_l3)) else character(),
      sources = split(manifest_df, seq_len(nrow(manifest_df)))
    ), dossier_manifest_json)

    header_lines <- c(
      sprintf("# %s / %s 跨方法 LLM 总整合", lineage, target_celltype),
      "",
      sprintf("- target_level: `%s`", target_level),
      sprintf("- safe_celltype: `%s`", target_safe)
    )
    if (!is.na(parent_l2_value) && nzchar(parent_l2_value)) header_lines <- c(header_lines, sprintf("- parent_l2: `%s`", parent_l2_value))
    if (nrow(child_l3) > 0L) header_lines <- c(header_lines, sprintf("- child_l3_count: `%d`", nrow(child_l3)))
    header_lines <- c(
      header_lines,
      sprintf("- dossier_manifest_tsv: `%s`", normalizePath(dossier_manifest_tsv, winslash = "/", mustWork = FALSE)),
      sprintf("- dossier_manifest_json: `%s`", normalizePath(dossier_manifest_json, winslash = "/", mustWork = FALSE)),
      ""
    )

    source_table_df <- manifest_df[, c("section_scope", "source_method", "item_label", "status", "inline_in_dossier"), drop = FALSE]
    lines <- c(
      header_lines,
      "## 来源清单",
      "",
      markdown_table_lines_unit_llm(source_table_df, c("section_scope", "source_method", "item_label", "status", "inline_in_dossier")),
      ""
    )

    if (target_level == "l2_only" && nrow(child_l3) > 0L) {
      child_df <- data.frame(
        celltype_l3 = child_l3$celltype_l3,
        safe_celltype = child_l3$safe_celltype,
        combined_md = child_l3$combined_md,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      lines <- c(lines, "## 关联 L3 子类型", "", markdown_table_lines_unit_llm(child_df, c("celltype_l3", "safe_celltype", "combined_md")), "")
    }

    if (nrow(l3_exact) > 0L) {
      for (j in seq_len(nrow(l3_exact))) {
        row <- l3_exact[j, , drop = FALSE]
        lines <- append_report_section_cross_method(
          lines,
          title = sprintf("## cNMF L3：%s", row$celltype_l3[[1]]),
          report_md = row$combined_md[[1]],
          meta_lines = c(
            sprintf("- total_units: `%d`", row$total_units[[1]]),
            sprintf("- ok_units: `%d`", row$ok_units[[1]]),
            sprintf("- error_units: `%d`", row$error_units[[1]]),
            sprintf("- combined_manifest_tsv: `%s`", row$combined_manifest_tsv[[1]])
          )
        )
      }
    }

    if (nrow(l2_context) > 0L) {
      for (j in seq_len(nrow(l2_context))) {
        row <- l2_context[j, , drop = FALSE]
        section_title <- if (target_level == "l2_only") {
          sprintf("## cNMF L2：%s", row$celltype_l2[[1]])
        } else {
          sprintf("## cNMF L2 背景：%s", row$celltype_l2[[1]])
        }
        lines <- append_report_section_cross_method(
          lines,
          title = section_title,
          report_md = row$combined_md[[1]],
          meta_lines = c(
            sprintf("- total_units: `%d`", row$total_units[[1]]),
            sprintf("- ok_units: `%d`", row$ok_units[[1]]),
            sprintf("- error_units: `%d`", row$error_units[[1]]),
            sprintf("- combined_manifest_tsv: `%s`", row$combined_manifest_tsv[[1]])
          )
        )
      }
    }

    if (nrow(covarnet_exact) > 0L) {
      covarnet_exact <- covarnet_exact[order(covarnet_exact$unit_id), , drop = FALSE]
      for (j in seq_len(nrow(covarnet_exact))) {
        row <- covarnet_exact[j, , drop = FALSE]
        lines <- append_report_section_cross_method(
          lines,
          title = sprintf("## CoVarNet：%s", row$unit_label[[1]]),
          report_md = row$report_md[[1]],
          meta_lines = c(
            sprintf("- unit_id: `%s`", row$unit_id[[1]]),
            sprintf("- status: `%s`", row$status[[1]]),
            sprintf("- status_json: `%s`", row$status_json[[1]])
          )
        )
      }
    }

    if (nrow(hdwgcna_exact) > 0L) {
      hdwgcna_exact <- hdwgcna_exact[order(hdwgcna_exact$module, hdwgcna_exact$unit_id), , drop = FALSE]
      for (j in seq_len(nrow(hdwgcna_exact))) {
        row <- hdwgcna_exact[j, , drop = FALSE]
        lines <- append_report_section_cross_method(
          lines,
          title = sprintf("## hdWGCNA：%s", row$unit_label[[1]]),
          report_md = row$report_md[[1]],
          meta_lines = c(
            sprintf("- unit_id: `%s`", row$unit_id[[1]]),
            sprintf("- module: `%s`", row$module[[1]]),
            sprintf("- status: `%s`", row$status[[1]]),
            sprintf("- status_json: `%s`", row$status_json[[1]])
          )
        )
      }
    }

    if (nrow(pycogaps_global) > 0L) {
      pycogaps_df <- pycogaps_global[, c("source_method", "status", "report_md", "input_index_tsv"), drop = FALSE]
      lines <- c(
        lines,
        "## Lineage-global 背景",
        "",
        "- `pycogaps` 只提供 lineage-global context，这里不把正文重复内嵌到每个 celltype dossier 里，避免同一份全局模式在多个细胞类型中反复膨胀。",
        "",
        markdown_table_lines_unit_llm(pycogaps_df, c("source_method", "status", "report_md", "input_index_tsv")),
        ""
      )
    }

    write_lines_unit_llm(lines, dossier_md)

    master_rows[[length(master_rows) + 1L]] <- data.frame(
      lineage = lineage,
      target_level = target_level,
      target_celltype = target_celltype,
      safe_celltype = target_safe,
      parent_l2 = parent_l2_value,
      child_l3_count = nrow(child_l3),
      has_cnmf_l2 = logical_yes_no_cross_method(nrow(l2_context) > 0L),
      has_cnmf_l3 = logical_yes_no_cross_method(nrow(l3_exact) > 0L),
      covarnet_units = nrow(covarnet_exact),
      hdwgcna_units = nrow(hdwgcna_exact),
      has_pycogaps_global = logical_yes_no_cross_method(nrow(pycogaps_global) > 0L),
      dossier_md = normalizePath(dossier_md, winslash = "/", mustWork = FALSE),
      dossier_manifest_tsv = normalizePath(dossier_manifest_tsv, winslash = "/", mustWork = FALSE),
      dossier_manifest_json = normalizePath(dossier_manifest_json, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    dossier_results[[paste(lineage, target_safe, sep = "||")]] <- list(
      target = target,
      manifest_df = manifest_df,
      dossier_md = normalizePath(dossier_md, winslash = "/", mustWork = FALSE),
      dossier_manifest_tsv = normalizePath(dossier_manifest_tsv, winslash = "/", mustWork = FALSE),
      dossier_manifest_json = normalizePath(dossier_manifest_json, winslash = "/", mustWork = FALSE)
    )
  }

  master_manifest_df <- if (length(master_rows) == 0L) empty_result$master_manifest_df else do.call(rbind, master_rows)
  master_manifest_tsv <- if (nrow(master_manifest_df) > 0L) file.path(run_root, "cross_method_celltype_LLM_dossier_manifest.tsv") else NA_character_
  master_manifest_json <- if (nrow(master_manifest_df) > 0L) file.path(run_root, "cross_method_celltype_LLM_dossier_manifest.json") else NA_character_
  if (!is.na(master_manifest_tsv)) write_tsv_unit_llm(master_manifest_df, master_manifest_tsv)
  if (!is.na(master_manifest_json)) write_json_unit_llm(list(dossiers = split(master_manifest_df, seq_len(nrow(master_manifest_df)))), master_manifest_json)

  list(
    master_manifest_df = master_manifest_df,
    master_manifest_tsv = if (!is.na(master_manifest_tsv)) normalizePath(master_manifest_tsv, winslash = "/", mustWork = FALSE) else NA_character_,
    master_manifest_json = if (!is.na(master_manifest_json)) normalizePath(master_manifest_json, winslash = "/", mustWork = FALSE) else NA_character_,
    dossier_results = dossier_results
  )
}

launch_cross_method_celltype_llm_main <- function() {
  run_root <- Sys.getenv("CROSS_METHOD_LLM_RUN_ROOT", unset = RUN_ROOT)
  lineages <- split_env_cross_method("CROSS_METHOD_LLM_LINEAGES", character())
  celltypes <- split_env_cross_method("CROSS_METHOD_LLM_CELLTYPES", character())
  root_methods <- split_env_cross_method("CROSS_METHOD_LLM_ROOT_METHODS", c("covarnet", "hdwgcna"))
  global_methods <- split_env_cross_method("CROSS_METHOD_LLM_GLOBAL_METHODS", c("pycogaps"))
  all_cnmf_k <- !(Sys.getenv("CROSS_METHOD_LLM_ALL_CNMF_K", unset = "1") %in% c("0", "false", "FALSE", "no", "NO"))

  cnmf_context <- collect_cnmf_combined_context_cross_method(
    run_root = run_root,
    lineages = lineages,
    celltypes = celltypes,
    all_cnmf_k = all_cnmf_k
  )
  root_unit_records <- collect_root_unit_records_cross_method(run_root = run_root, lineages = lineages, source_methods = root_methods)
  global_reports <- collect_global_method_reports_cross_method(run_root = run_root, lineages = lineages, source_methods = global_methods)
  combined <- write_cross_method_celltype_dossiers(
    run_root = run_root,
    l2_records = cnmf_context$l2_records,
    l3_records = cnmf_context$l3_records,
    root_unit_records = root_unit_records,
    global_reports = global_reports,
    output_subdir = OUTPUT_SUBDIR
  )

  cat(sprintf("[INFO] cross-method celltype dossiers: %d targets\n", nrow(combined$master_manifest_df)))
  if (!is.na(combined$master_manifest_tsv)) {
    cat(sprintf("[OK] cross-method master manifest: %s\n", combined$master_manifest_tsv))
  }
  invisible(list(cnmf_context = cnmf_context, root_unit_records = root_unit_records, global_reports = global_reports, combined = combined))
}

if (identical(sys.nframe(), 0L)) {
  launch_cross_method_celltype_llm_main()
}