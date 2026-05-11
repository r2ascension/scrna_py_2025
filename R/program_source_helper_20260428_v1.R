#!/usr/bin/env Rscript
# ==============================================================================
# Program Source Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
PA_PROGRAM_REGISTRY_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_registry_helper_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)
if (!exists("pa_register_programs", mode = "function")) source(PA_PROGRAM_REGISTRY_HELPER_PATH_20260428_V1)

PA_HDWGCNA_COVARNET_HELPER_PATH_20260423_V1_1 <- "/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R"
PA_COEXPR_VISUALIZATION_PATCH_PATH_20260508_V1_2 <- "/home/h2048/script/R/coexpr_visualization_patch_v1_2.R"
PA_CNMF_HELPER_PY_PATH_20260419_V1_1 <- "/home/h2048/script/py/cnmf_helper_20260419_v1_1.py"

pa_safe_file_id <- function(x) {
  out <- gsub("[^A-Za-z0-9_]+", "_", pa_safe_trim(x))
  out[!nzchar(out)] <- "NA"
  out
}

pa_source_coexpr_helper <- function(helper_path = PA_HDWGCNA_COVARNET_HELPER_PATH_20260423_V1_1) {
  helper_path <- pa_scalar_chr(helper_path, "helper_path")
  if (!file.exists(helper_path)) {
    stop(sprintf("Coexpression helper file does not exist: %s", helper_path), call. = FALSE)
  }
  source(helper_path, local = .GlobalEnv)
  if (exists("load_coexpr_libs", mode = "function", inherits = TRUE)) {
    load_coexpr_libs()
  }
  if (file.exists(PA_COEXPR_VISUALIZATION_PATCH_PATH_20260508_V1_2)) {
    source(PA_COEXPR_VISUALIZATION_PATCH_PATH_20260508_V1_2, local = .GlobalEnv)
  }
  invisible(helper_path)
}

pa_prepare_program_source_tbl <- function(tbl, id_col, gene_col, object_name = "program_tbl") {
  pa_validate_required_columns(tbl, c(id_col, gene_col), object_name)
  out <- data.frame(
    program_id = pa_safe_trim(tbl[[id_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$gene_vector <- I(lapply(tbl[[gene_col]], pa_normalize_gene_vector))
  if ("optional_weight" %in% colnames(tbl)) out$optional_weight <- tbl$optional_weight
  out
}

pa_register_cnmf_programs <- function(
  unit_id,
  program_tbl,
  program_id_col = "program_id",
  gene_col = "gene_vector",
  subtype = "GEP",
  score_level = "cell",
  score_object_path = NA_character_,
  validation_status = "raw",
  usage_col = "usage_column"
) {
  reg_input <- pa_prepare_program_source_tbl(program_tbl, program_id_col, gene_col, "cNMF program_tbl")
  out <- pa_register_programs(
    unit_id = unit_id,
    source_type = "cNMF",
    source_subtype = subtype,
    program_tbl = reg_input,
    score_level = score_level,
    score_object_path = score_object_path,
    validation_status = validation_status
  )
  if (usage_col %in% colnames(program_tbl)) out$usage_column <- pa_safe_trim(program_tbl[[usage_col]])
  out
}

pa_register_hdwgcna_programs <- function(
  unit_id,
  module_tbl,
  module_id_col = "module_id",
  gene_col = "gene_vector",
  subtype = "module",
  score_level = "metacell",
  score_object_path = NA_character_,
  validation_status = "raw"
) {
  reg_input <- pa_prepare_program_source_tbl(module_tbl, module_id_col, gene_col, "hdWGCNA module_tbl")
  reg_input$program_id <- paste0("hdwgcna_", reg_input$program_id)
  out <- pa_register_programs(
    unit_id = unit_id,
    source_type = "hdWGCNA",
    source_subtype = subtype,
    program_tbl = reg_input,
    score_level = score_level,
    score_object_path = score_object_path,
    validation_status = validation_status
  )
  out$module_id <- pa_safe_trim(module_tbl[[module_id_col]])
  out
}

pa_register_covarnet_programs <- function(
  unit_id,
  program_tbl,
  program_id_col = "program_id",
  gene_col = "gene_vector",
  subtype = "hub_gene_set",
  score_level = "celltype",
  score_object_path = NA_character_,
  validation_status = "raw",
  hub_genes_col = "hub_genes"
) {
  reg_input <- pa_prepare_program_source_tbl(program_tbl, program_id_col, gene_col, "CoVarNet program_tbl")
  out <- pa_register_programs(
    unit_id = unit_id,
    source_type = "CoVarNet",
    source_subtype = subtype,
    program_tbl = reg_input,
    score_level = score_level,
    score_object_path = score_object_path,
    validation_status = validation_status
  )
  if (hub_genes_col %in% colnames(program_tbl)) {
    out$hub_genes <- I(lapply(program_tbl[[hub_genes_col]], pa_normalize_gene_vector))
  }
  out
}

pa_bind_program_registries <- function(...) {
  regs <- list(...)
  if (length(regs) == 1L && is.list(regs[[1]]) && !is.data.frame(regs[[1]])) regs <- regs[[1]]
  regs <- regs[!vapply(regs, is.null, logical(1))]
  if (length(regs) == 0L) return(data.frame())
  if (!all(vapply(regs, is.data.frame, logical(1)))) {
    stop("All program registries must be data.frame objects", call. = FALSE)
  }

  all_cols <- unique(unlist(lapply(regs, colnames), use.names = FALSE))
  list_cols <- unique(unlist(lapply(regs, function(df) names(df)[vapply(df, is.list, logical(1))]), use.names = FALSE))

  regs_aligned <- lapply(regs, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    for (nm in missing_cols) {
      if (nm %in% list_cols) {
        df[[nm]] <- vector("list", nrow(df))
      } else {
        df[[nm]] <- rep(NA, nrow(df))
      }
    }
    df[, all_cols, drop = FALSE]
  })

  out <- do.call(rbind, regs_aligned)
  rownames(out) <- NULL
  out
}

pa_import_hdwgcna_results <- function(output_dir, unit_id = NULL) {
  output_dir <- pa_prepare_output_dir(output_dir)
  membership_files <- list.files(
    output_dir,
    pattern = "hdwgcna_module_membership\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(membership_files) == 0L) {
    return(list(
      program_tbl = data.frame(),
      registry = NULL,
      manifest = list(output_dir = output_dir, module_membership_files = character(), n_programs = 0L)
    ))
  }

  tbls <- lapply(membership_files, function(path) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    module_col <- c("module", "module_id")[c("module", "module_id") %in% colnames(df)][1]
    gene_col <- c("gene_name", "gene", "feature")[c("gene_name", "gene", "feature") %in% colnames(df)][1]
    if (is.na(module_col) || is.na(gene_col)) return(NULL)
    df <- df[!is.na(df[[module_col]]) & nzchar(trimws(df[[module_col]])), , drop = FALSE]
    df <- df[pa_safe_trim(df[[module_col]]) != "grey", , drop = FALSE]
    if (nrow(df) == 0L) return(NULL)

    celltype_label <- sub("^hdwgcna_", "", basename(dirname(path)))
    split_genes <- split(df[[gene_col]], df[[module_col]])
    out <- data.frame(
      module_id = paste(pa_safe_file_id(celltype_label), names(split_genes), sep = "__"),
      celltype = rep(celltype_label, length(split_genes)),
      module_label = names(split_genes),
      module_membership_path = rep(normalizePath(path, winslash = "/", mustWork = FALSE), length(split_genes)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    out$gene_vector <- I(lapply(split_genes, pa_normalize_gene_vector))
    out
  })

  tbls <- Filter(Negate(is.null), tbls)
  if (length(tbls) == 0L) {
    return(list(
      program_tbl = data.frame(),
      registry = NULL,
      manifest = list(output_dir = output_dir, module_membership_files = membership_files, n_programs = 0L)
    ))
  }

  program_tbl <- do.call(rbind, tbls)
  rownames(program_tbl) <- NULL

  registry <- NULL
  if (!is.null(unit_id) && nzchar(pa_safe_trim(unit_id))) {
    registry <- pa_register_hdwgcna_programs(
      unit_id = unit_id,
      module_tbl = program_tbl,
      module_id_col = "module_id",
      gene_col = "gene_vector",
      score_object_path = output_dir
    )
    registry$celltype <- program_tbl$celltype
    registry$module_label <- program_tbl$module_label
    registry$module_membership_path <- program_tbl$module_membership_path
  }

  list(
    program_tbl = program_tbl,
    registry = registry,
    manifest = list(
      output_dir = output_dir,
      module_membership_files = membership_files,
      n_programs = nrow(program_tbl)
    )
  )
}

pa_run_hdwgcna_runner <- function(seurat_obj,
                                  output_dir,
                                  celltypes = NULL,
                                  helper_path = PA_HDWGCNA_COVARNET_HELPER_PATH_20260423_V1_1,
                                  unit_id = NULL,
                                  celltype_col = NULL,
                                  sample_col = NULL,
                                  condition_col = NULL,
                                  tissue_col = NULL,
                                  dry_run = FALSE,
                                  global_setup_args = list(),
                                  runner_args = list()) {
  output_dir <- pa_prepare_output_dir(output_dir)
  plan <- list(
    status = if (isTRUE(dry_run)) "dry_run" else "ready",
    helper_path = helper_path,
    output_dir = output_dir,
    celltypes = pa_null_coalesce(celltypes, character())
  )
  if (isTRUE(dry_run)) return(plan)

  pa_source_coexpr_helper(helper_path)
  if (!is.null(celltype_col) && nzchar(pa_safe_trim(celltype_col))) assign("CELLTYPE_COL", celltype_col, envir = .GlobalEnv)
  if (!is.null(sample_col) && nzchar(pa_safe_trim(sample_col))) assign("SAMPLE_COL", sample_col, envir = .GlobalEnv)
  if (!is.null(condition_col) && nzchar(pa_safe_trim(condition_col))) assign("CONDITION_COL", condition_col, envir = .GlobalEnv)
  if (!is.null(tissue_col) && nzchar(pa_safe_trim(tissue_col))) assign("TISSUE_COL", tissue_col, envir = .GlobalEnv)
  global_setup_args <- utils::modifyList(list(out_dir = output_dir), pa_null_coalesce(global_setup_args, list()))
  setup_obj <- do.call(hdwgcna_global_setup, c(list(seurat_obj = seurat_obj), global_setup_args))
  hdwgcna_results <- do.call(
    hdwgcna_run_celltypes,
    c(list(seurat_obj = setup_obj, celltypes = celltypes, out_dir = output_dir), runner_args)
  )
  imported <- pa_import_hdwgcna_results(output_dir = output_dir, unit_id = unit_id)
  c(imported, list(status = "ok", hdwgcna_results = hdwgcna_results, plan = plan))
}

pa_covarnet_results_to_program_tbl <- function(results,
                                               output_dir,
                                               n_hubs = 25L) {
  if (is.null(results) || !is.list(results) || length(results) == 0L) {
    return(data.frame())
  }

  rows <- list()
  idx <- 1L
  for (celltype_label in names(results)) {
    res <- results[[celltype_label]]
    if (is.null(res$hubs) || is.null(res$hubs$hub_genes) || !is.data.frame(res$hubs$hub_genes) || nrow(res$hubs$hub_genes) == 0L) next
    hub_df <- res$hubs$hub_genes
    gene_col <- c("gene", "gene_name")[c("gene", "gene_name") %in% colnames(hub_df)][1]
    if (is.na(gene_col)) next
    hub_genes <- head(pa_normalize_gene_vector(hub_df[[gene_col]]), max(1L, as.integer(n_hubs)))
    edge_path <- file.path(output_dir, sprintf("covarnet_%s_edges.csv", pa_safe_file_id(celltype_label)))
    node_path <- file.path(output_dir, sprintf("covarnet_%s_nodes.csv", pa_safe_file_id(celltype_label)))
    row <- data.frame(
      program_id = sprintf("covarnet_%s", pa_safe_file_id(celltype_label)),
      celltype = celltype_label,
      edge_path = normalizePath(edge_path, winslash = "/", mustWork = FALSE),
      node_path = normalizePath(node_path, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    row$gene_vector <- I(list(hub_genes))
    row$hub_genes <- I(list(hub_genes))
    rows[[idx]] <- row
    idx <- idx + 1L
  }

  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

pa_run_covarnet_runner <- function(seurat_obj,
                                   output_dir,
                                   celltypes = NULL,
                                   helper_path = PA_HDWGCNA_COVARNET_HELPER_PATH_20260423_V1_1,
                                   unit_id = NULL,
                                   celltype_col = NULL,
                                   dry_run = FALSE,
                                   plot_networks = FALSE,
                                   plot_overlap = FALSE,
                                   enhanced_visualizations = TRUE,
                                   plot_raw_correlation = FALSE,
                                   n_hubs = 25L,
                                   cor_thr = NULL,
                                   pval_thr = NULL,
                                   ...) {
  output_dir <- pa_prepare_output_dir(output_dir)
  plan <- list(
    status = if (isTRUE(dry_run)) "dry_run" else "ready",
    helper_path = helper_path,
    output_dir = output_dir,
    celltypes = pa_null_coalesce(celltypes, character())
  )
  if (isTRUE(dry_run)) return(plan)

  pa_source_coexpr_helper(helper_path)
  if (!is.null(celltype_col) && nzchar(pa_safe_trim(celltype_col))) assign("CELLTYPE_COL", celltype_col, envir = .GlobalEnv)

  old_cor_thr <- if (exists("COVAR_COR_THR", envir = .GlobalEnv, inherits = FALSE)) get("COVAR_COR_THR", envir = .GlobalEnv) else NULL
  old_pval_thr <- if (exists("COVAR_PVAL_THR", envir = .GlobalEnv, inherits = FALSE)) get("COVAR_PVAL_THR", envir = .GlobalEnv) else NULL
  if (!is.null(cor_thr)) assign("COVAR_COR_THR", cor_thr, envir = .GlobalEnv)
  if (!is.null(pval_thr)) assign("COVAR_PVAL_THR", pval_thr, envir = .GlobalEnv)
  on.exit({
    if (!is.null(old_cor_thr)) assign("COVAR_COR_THR", old_cor_thr, envir = .GlobalEnv)
    if (!is.null(old_pval_thr)) assign("COVAR_PVAL_THR", old_pval_thr, envir = .GlobalEnv)
  }, add = TRUE)

  covarnet_results <- covarnet_run_all(seurat_obj = seurat_obj, celltypes = celltypes, out_dir = output_dir, ...)

  if (isTRUE(enhanced_visualizations) && exists("covarnet_plot_all", mode = "function", inherits = TRUE)) {
    tryCatch(
      covarnet_plot_all(
        covarnet_results,
        out_dir = output_dir,
        n_top_hubs = n_hubs,
        plot_raw_correlation = isTRUE(plot_raw_correlation)
      ),
      error = function(e) cat(sprintf("[WARN] Enhanced CoVarNet visualizations failed: %s\n", e$message))
    )
  }

  if (isTRUE(plot_networks)) {
    for (ct in names(covarnet_results)) {
      plot_covarnet_graph(covarnet_results[[ct]]$graph, celltype_name = ct, out_dir = output_dir)
    }
  }
  if (isTRUE(plot_overlap)) {
    plot_hub_overlap(covarnet_results, out_dir = output_dir)
  }

  program_tbl <- pa_covarnet_results_to_program_tbl(covarnet_results, output_dir = output_dir, n_hubs = n_hubs)
  registry <- NULL
  if (nrow(program_tbl) > 0L && !is.null(unit_id) && nzchar(pa_safe_trim(unit_id))) {
    registry <- pa_register_covarnet_programs(
      unit_id = unit_id,
      program_tbl = program_tbl,
      program_id_col = "program_id",
      gene_col = "gene_vector",
      hub_genes_col = "hub_genes",
      score_object_path = output_dir
    )
    registry$celltype <- program_tbl$celltype
    registry$edge_path <- program_tbl$edge_path
    registry$node_path <- program_tbl$node_path
  }

  list(
    status = "ok",
    plan = plan,
    covarnet_results = covarnet_results,
    program_tbl = program_tbl,
    registry = registry,
    manifest = list(
      output_dir = output_dir,
      edge_files = list.files(output_dir, pattern = "_edges\\.csv$", full.names = TRUE),
      node_files = list.files(output_dir, pattern = "_nodes\\.csv$", full.names = TRUE),
      n_programs = nrow(program_tbl)
    )
  )
}

pa_write_cnmf_launcher <- function(path) {
  lines <- c(
    "#!/usr/bin/env python3",
    "import importlib.util",
    "import json",
    "import pathlib",
    "import scanpy as sc",
    "import sys",
    "",
    "cfg_path = pathlib.Path(sys.argv[1])",
    "cfg = json.loads(cfg_path.read_text())",
    "helper_path = pathlib.Path(cfg['helper_py'])",
    "spec = importlib.util.spec_from_file_location('pa_cnmf_helper', helper_path)",
    "module = importlib.util.module_from_spec(spec)",
    "spec.loader.exec_module(module)",
    "adata = sc.read_h5ad(cfg['adata_h5ad_path'])",
    "res = module.run_cnmf_full(",
    "    adata=adata,",
    "    output_dir=pathlib.Path(cfg['output_dir']),",
    "    run_name=cfg['run_name'],",
    "    k_range=cfg.get('k_range'),",
    "    celltype_col=cfg.get('celltype_col'),",
    "    batch_col=cfg.get('batch_col'),",
    "    use_batch_hvg=cfg.get('use_batch_hvg', True),",
    "    cnmf_config=cfg.get('cnmf_config') or {},",
    "    viz_config=cfg.get('viz_config') or {}",
    ")",
    "out_path = pathlib.Path(cfg['output_dir']) / 'cnmf_runner_result.json'",
    "out_path.write_text(json.dumps(res, indent=2, default=str), encoding='utf-8')",
    "print(json.dumps({'success': bool(res.get('success')), 'result_json': str(out_path)}, ensure_ascii=False))"
  )
  pa_write_markdown(lines, path)
  Sys.chmod(path, mode = "0755")
  invisible(path)
}

pa_run_cnmf_runner <- function(adata_h5ad_path,
                               output_dir,
                               run_name,
                               helper_py = PA_CNMF_HELPER_PY_PATH_20260419_V1_1,
                               python_cmd = NULL,
                               unit_id = NULL,
                               dry_run = FALSE,
                               k_range = NULL,
                               celltype_col = NULL,
                               batch_col = NULL,
                               use_batch_hvg = TRUE,
                               cnmf_config = list(),
                               viz_config = list()) {
  adata_h5ad_path <- pa_scalar_chr(adata_h5ad_path, "adata_h5ad_path")
  if (!file.exists(helper_py)) stop(sprintf("cNMF helper Python file does not exist: %s", helper_py), call. = FALSE)
  output_dir <- pa_prepare_output_dir(output_dir)
  py_exec <- pa_resolve_python_executable(python_cmd)

  config <- list(
    helper_py = normalizePath(helper_py, winslash = "/", mustWork = FALSE),
    adata_h5ad_path = normalizePath(adata_h5ad_path, winslash = "/", mustWork = FALSE),
    output_dir = output_dir,
    run_name = pa_scalar_chr(run_name, "run_name"),
    k_range = k_range,
    celltype_col = pa_null_coalesce(celltype_col, NULL),
    batch_col = pa_null_coalesce(batch_col, NULL),
    use_batch_hvg = isTRUE(use_batch_hvg),
    cnmf_config = pa_null_coalesce(cnmf_config, list()),
    viz_config = pa_null_coalesce(viz_config, list())
  )
  config_path <- file.path(output_dir, "cnmf_runner_config.json")
  launcher_path <- file.path(output_dir, "run_cnmf_runner_launcher.py")
  pa_write_json(config, config_path)
  pa_write_cnmf_launcher(launcher_path)

  plan <- list(
    status = if (isTRUE(dry_run)) "dry_run" else "ready",
    python = py_exec,
    config_path = config_path,
    launcher_path = launcher_path,
    output_dir = output_dir
  )
  if (isTRUE(dry_run)) return(plan)

  cmd_res <- pa_run_system_command(py_exec, args = c(launcher_path, config_path), fail_on_error = FALSE)
  runner_result <- pa_json_read(file.path(output_dir, "cnmf_runner_result.json"), default = list(success = FALSE, output = cmd_res$output))
  imported <- NULL
  if (isTRUE(runner_result$success)) {
    imported <- pa_import_cnmf_results(output_dir = output_dir, unit_id = unit_id)
  }

  list(
    status = if (isTRUE(runner_result$success)) "ok" else "error",
    plan = plan,
    command_result = cmd_res,
    runner_result = runner_result,
    program_tbl = if (!is.null(imported)) imported$program_tbl else data.frame(),
    registry = if (!is.null(imported)) imported$registry else NULL,
    import = imported
  )
}

pa_import_cnmf_results <- function(output_dir,
                                   unit_id = NULL,
                                   selected_k = NULL) {
  output_dir <- pa_prepare_output_dir(output_dir)
  recommendation <- pa_json_read(file.path(output_dir, "k_selection_recommendation.json"), default = list())
  run_summary <- pa_json_read(file.path(output_dir, "run_summary.json"), default = list())

  if (is.null(selected_k)) {
    selected_k <- pa_null_coalesce(run_summary$recommendation$recommended_k, recommendation$recommended_k)
  }
  if (!is.null(selected_k)) selected_k <- suppressWarnings(as.integer(unlist(selected_k, use.names = FALSE)[1]))

  gep_dir <- file.path(output_dir, "gep_gene_tables")
  score_files <- list.files(gep_dir, pattern = "^gep_gene_scores_k[0-9]+\\.tsv$", full.names = TRUE)
  if (length(score_files) == 0L) {
    return(list(
      program_tbl = data.frame(),
      registry = NULL,
      manifest = list(output_dir = output_dir, score_files = character(), selected_k = selected_k),
      recommendation = recommendation,
      run_summary = run_summary,
      selected_k = selected_k
    ))
  }

  if (length(selected_k) == 0L || is.na(selected_k)) {
    parsed_k <- suppressWarnings(as.integer(sub("^.*_k([0-9]+)\\.tsv$", "\\1", score_files)))
    selected_k <- parsed_k[which.max(parsed_k)]
  }

  score_path <- file.path(gep_dir, sprintf("gep_gene_scores_k%d.tsv", selected_k))
  if (!file.exists(score_path)) {
    score_path <- score_files[1]
  }
  score_tbl <- pa_read_table_auto(score_path)
  pa_validate_required_columns(score_tbl, c("gep", "gene"), "cNMF gene score table")
  if (!"rank" %in% colnames(score_tbl)) score_tbl$rank <- ave(seq_len(nrow(score_tbl)), score_tbl$gep, FUN = seq_along)

  split_genes <- lapply(split(score_tbl, score_tbl$gep), function(df) {
    df <- df[order(df$rank), , drop = FALSE]
    pa_normalize_gene_vector(df$gene)
  })
  program_tbl <- data.frame(
    program_id = names(split_genes),
    usage_column = names(split_genes),
    selected_k = rep(selected_k, length(split_genes)),
    score_table_path = rep(normalizePath(score_path, winslash = "/", mustWork = FALSE), length(split_genes)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  program_tbl$gene_vector <- I(unname(split_genes))

  registry <- NULL
  if (!is.null(unit_id) && nzchar(pa_safe_trim(unit_id))) {
    registry <- pa_register_cnmf_programs(
      unit_id = unit_id,
      program_tbl = program_tbl,
      program_id_col = "program_id",
      gene_col = "gene_vector",
      usage_col = "usage_column",
      score_object_path = score_path
    )
    registry$selected_k <- program_tbl$selected_k
    registry$score_table_path <- program_tbl$score_table_path
  }

  list(
    program_tbl = program_tbl,
    registry = registry,
    manifest = list(
      output_dir = output_dir,
      score_files = score_files,
      selected_k = selected_k,
      score_table_path = normalizePath(score_path, winslash = "/", mustWork = FALSE)
    ),
    recommendation = recommendation,
    run_summary = run_summary,
    selected_k = selected_k
  )
}

if (sys.nframe() == 0) {
  cat("Program Source Helper (2026-04-28 v1) loaded.\n")
}
