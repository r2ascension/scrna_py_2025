#!/usr/bin/env Rscript

# Launch per-unit LLM interpretation for L2 cell-type cNMF outputs created by
# script/py/run_cnmf_by_l2_celltype_20260521.py.
# Preferred layout is `<run_root>/<lineage>/cnmf_by_celltype_l2/<safe_L2>/cnmf_full/`.
# Legacy `cnmf_<safe_L2>/cnmf_full/` directories are still discovered.

UNIT_LLM_SCRIPT <- "/home/h2048/script/R/launch_program_unit_llm_20260513.R"
source(UNIT_LLM_SCRIPT)

RUN_ROOT <- Sys.getenv("CNMF_L2_RUN_ROOT", unset = "/home/h2048/output/program_full_parallel_methods_20260507")
OUTPUT_SUBDIR <- Sys.getenv("CNMF_L2_OUTPUT_SUBDIR", unset = "cnmf_by_celltype_l2")

split_env_l2 <- function(name, default = character()) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  out <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

discover_l2_cnmf_sources <- function(l2_root) {
  candidate_dirs <- sort(unique(c(
    Sys.glob(file.path(l2_root, "*", "cnmf_full")),
    Sys.glob(file.path(l2_root, "cnmf_*", "cnmf_full"))
  )))
  if (length(candidate_dirs) == 0L) return(list())
  records <- lapply(candidate_dirs, function(source_dir) {
    source_dir <- normalizePath(source_dir, winslash = "/", mustWork = FALSE)
    status <- read_json_unit_llm(file.path(source_dir, "cnmf_l2_status.json"), default = list())
    ct_dir <- dirname(source_dir)
    ct_basename <- basename(ct_dir)
    status_ok <- identical(as.character(status$status)[1], "ok")
    path_layout <- as.character(status$path_layout %||% NA_character_)[1]
    is_legacy <- !identical(path_layout, "celltype_method") && grepl("^cnmf_", ct_basename)
    safe_celltype <- safe_unit_id(status$safe_celltype %||% if (grepl("^cnmf_", ct_basename)) sub("^cnmf_", "", ct_basename) else ct_basename)
    celltype_l2 <- as.character(status$celltype_l2 %||% status$cell_type_L2 %||% NA_character_)[1]
    list(
      key = if (!is.na(celltype_l2) && nzchar(celltype_l2)) paste0("celltype:", celltype_l2) else paste0("safe:", safe_celltype),
      status_ok = status_ok,
      is_legacy = is_legacy,
      source_dir = source_dir,
      status = status
    )
  })
  ord <- order(
    vapply(records, function(x) x$key, character(1)),
    !vapply(records, function(x) x$status_ok, logical(1)),
    vapply(records, function(x) x$is_legacy, logical(1)),
    vapply(records, function(x) x$source_dir, character(1))
  )
  records <- records[ord]
  records[!duplicated(vapply(records, function(x) x$key, character(1)))]
}

pick_primary_value_l2 <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x) & !is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  counts <- sort(table(x), decreasing = TRUE)
  names(counts)[1]
}

infer_l2_hierarchy <- function(source_dir, status = list(), fallback_l2 = NA_character_) {
  meta_path <- file.path(source_dir, "cell_metadata_for_l2_cnmf.csv")
  meta <- tryCatch(read_table_unit_llm(meta_path), error = function(e) data.frame())
  l1 <- as.character(status$celltype_l1 %||% status$cell_type_L1 %||% NA_character_)[1]
  l2 <- as.character(status$celltype_l2 %||% status$cell_type_L2 %||% fallback_l2 %||% NA_character_)[1]
  if (is.data.frame(meta) && nrow(meta) > 0L) {
    if ((!nzchar(l1) || is.na(l1)) && "cell_type_L1" %in% colnames(meta)) l1 <- pick_primary_value_l2(meta$cell_type_L1)
    if ((!nzchar(l2) || is.na(l2)) && "cell_type_L2" %in% colnames(meta)) l2 <- pick_primary_value_l2(meta$cell_type_L2)
  }
  list(
    celltype_l1 = if (!is.na(l1) && nzchar(l1)) l1 else NA_character_,
    celltype_l2 = if (!is.na(l2) && nzchar(l2)) l2 else fallback_l2,
    hierarchy_label = if (!is.na(l2) && nzchar(l2)) l2 else fallback_l2,
    metadata_path = if (file.exists(meta_path)) normalizePath(meta_path, winslash = "/", mustWork = FALSE) else NA_character_
  )
}

build_l2_cnmf_unit_rows <- function(run_root = RUN_ROOT,
                                    lineages = character(),
                                    celltypes = character(),
                                    all_cnmf_k = TRUE) {
  run_root <- normalizePath(run_root, winslash = "/", mustWork = FALSE)
  if (length(lineages) == 0L) {
    dirs <- list.dirs(run_root, recursive = FALSE, full.names = FALSE)
    lineages <- sort(dirs[nzchar(dirs)])
  }
  rows <- list()
  for (lineage in lineages) {
    l2_root <- file.path(run_root, lineage, OUTPUT_SUBDIR)
    if (!dir.exists(l2_root)) next
    source_records <- discover_l2_cnmf_sources(l2_root)
    for (record in source_records) {
      source_dir <- record$source_dir
      status <- record$status
      if (!isTRUE(record$status_ok)) next
      ct_basename <- basename(dirname(source_dir))
      celltype <- as.character(status$celltype_l2 %||% if (grepl("^cnmf_", ct_basename)) sub("^cnmf_", "", ct_basename) else ct_basename)
      hierarchy <- infer_l2_hierarchy(source_dir, status = status, fallback_l2 = celltype)
      celltype_l2 <- as.character(hierarchy$celltype_l2 %||% celltype)
      celltype_l1 <- as.character(hierarchy$celltype_l1 %||% NA_character_)
      hierarchy_label <- as.character(hierarchy$hierarchy_label %||% celltype)
      safe_celltype <- safe_unit_id(status$safe_celltype %||% celltype)
      if (length(celltypes) > 0L && !celltype %in% celltypes && !safe_celltype %in% celltypes) next
      score_files <- selected_cnmf_score_files(source_dir, all_cnmf_k = all_cnmf_k)
      if (length(score_files) == 0L) next
      for (score_path in score_files) {
        score_tbl <- tryCatch(read_table_unit_llm(score_path), error = function(e) data.frame())
        if (!is.data.frame(score_tbl) || nrow(score_tbl) == 0L || !all(c("gep", "gene") %in% colnames(score_tbl))) next
        k <- suppressWarnings(as.integer(sub("^.*_k([0-9]+)\\.tsv$", "\\1", score_path)))
        for (gep in unique(as.character(score_tbl$gep))) {
          unit_id <- safe_unit_id(sprintf("cnmf_l2_%s_k%s_%s", safe_celltype, k, gep))
          out_dir <- file.path(l2_root, "llm_parallel", "cnmf", "units", unit_id)
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
          gep_gene_path <- file.path(out_dir, paste0(unit_id, "_gene_scores.csv"))
          gep_tbl <- score_tbl[as.character(score_tbl$gep) == gep, , drop = FALSE]
          gep_tbl <- filter_program_low_information_genes(gep_tbl, gene_col = "gene")
          tryCatch(utils::write.csv(gep_tbl, gep_gene_path, row.names = FALSE), error = function(e) NULL)
          fig_paths <- list.files(file.path(source_dir, "visualizations"), pattern = sprintf("k%s\\.", k), full.names = TRUE)
          evidence_paths <- c(gep_gene_path, file.path(source_dir, "k_selection_recommendation.json"), hierarchy$metadata_path, fig_paths)
          rows[[length(rows) + 1L]] <- new_task_row(
            lineage = lineage,
            source_method = "cnmf",
            unit_type = "l2_gep",
            unit_id = unit_id,
            unit_label = sprintf("%s | k%s / %s", hierarchy_label, k, gep),
            source_dir = source_dir,
            output_dir = out_dir,
            evidence_paths = evidence_paths,
            extra = list(
              celltype_l1 = celltype_l1,
              celltype_l2 = celltype_l2,
              safe_celltype = safe_celltype,
              k = k,
              gep = gep,
              gep_gene_table = normalizePath(gep_gene_path, winslash = "/", mustWork = FALSE),
              score_table = normalizePath(score_path, winslash = "/", mustWork = FALSE),
              cnmf_l2_status_json = normalizePath(file.path(source_dir, "cnmf_l2_status.json"), winslash = "/", mustWork = FALSE),
              cell_metadata_csv = hierarchy$metadata_path
            )
          )
        }
      }
    }
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

launch_cnmf_l2_unit_llm_main <- function() {
  lineages <- split_env_l2("CNMF_L2_LLM_LINEAGES", split_env_l2("CNMF_L2_LINEAGES", character()))
  celltypes <- split_env_l2("CNMF_L2_LLM_CELLTYPES", split_env_l2("CNMF_L2_CELLTYPES", character()))
  all_cnmf_k <- !(Sys.getenv("CNMF_L2_LLM_ALL_K", unset = "1") %in% c("0", "false", "FALSE", "no", "NO"))
  max_parallel <- suppressWarnings(as.integer(Sys.getenv("CNMF_L2_LLM_MAX_PARALLEL", unset = "8")))
  if (length(max_parallel) == 0L || is.na(max_parallel) || max_parallel <= 0L) max_parallel <- 8L
  force <- Sys.getenv("CNMF_L2_LLM_FORCE", unset = "0") %in% c("1", "true", "TRUE", "yes", "YES")
  enable_live <- !(Sys.getenv("CNMF_L2_LLM_ENABLE_LIVE", unset = "1") %in% c("0", "false", "FALSE", "no", "NO"))
  model <- Sys.getenv("CNMF_L2_LLM_MODEL", unset = "deepseek-reasoner")
  timeout_sec <- suppressWarnings(as.integer(Sys.getenv("CNMF_L2_LLM_TIMEOUT_SEC", unset = "240")))
  if (length(timeout_sec) == 0L || is.na(timeout_sec) || timeout_sec <= 0L) timeout_sec <- 240L

  idx <- build_l2_cnmf_unit_rows(RUN_ROOT, lineages = lineages, celltypes = celltypes, all_cnmf_k = all_cnmf_k)
  index_path <- file.path(RUN_ROOT, sprintf("cnmf_l2_unit_llm_task_index_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (exists("pa_write_tsv", mode = "function")) pa_write_tsv(idx, index_path) else utils::write.table(idx, file = index_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("[INFO] L2 cNMF unit LLM task index: %s (%d tasks)\n", index_path, nrow(idx)))
  results <- run_unit_llm_tasks(idx, max_parallel = max_parallel, enable_live = enable_live, force = force, model = model, timeout_sec = timeout_sec)
  summary <- data.frame(
    unit_id = vapply(results, function(x) as.character(x$unit_id %||% NA_character_), character(1)),
    status = vapply(results, function(x) as.character(x$status %||% NA_character_), character(1)),
    status_json = vapply(results, function(x) as.character(x$status_json %||% NA_character_), character(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  summary_path <- file.path(RUN_ROOT, sprintf("cnmf_l2_unit_llm_summary_%s.tsv", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (exists("pa_write_tsv", mode = "function")) pa_write_tsv(summary, summary_path) else utils::write.table(summary, file = summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
  combined <- combine_cnmf_unit_llm_by_celltype(
    idx,
    run_root = RUN_ROOT,
    output_subdir = OUTPUT_SUBDIR,
    level_tag = "l2",
    celltype_field = "celltype_l2"
  )
  if (is.data.frame(combined$master_manifest_df) && nrow(combined$master_manifest_df) > 0L) {
    cat(sprintf("[OK] L2 cNMF celltype combined LLM manifest: %s (%d celltypes)\n", combined$master_manifest_tsv, nrow(combined$master_manifest_df)))
  }
  cat(sprintf("[OK] L2 cNMF unit LLM launcher finished: %s\n", summary_path))
  invisible(list(index = idx, results = results, summary = summary, index_path = index_path, summary_path = summary_path, combined = combined))
}

if (identical(sys.nframe(), 0L)) {
  launch_cnmf_l2_unit_llm_main()
}
