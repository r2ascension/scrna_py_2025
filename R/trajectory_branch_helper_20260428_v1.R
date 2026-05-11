#!/usr/bin/env Rscript
# ==============================================================================
# Trajectory Branch Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)

pa_write_paga_launcher <- function(path) {
  lines <- c(
    "#!/usr/bin/env python3",
    "import json",
    "import pathlib",
    "import scanpy as sc",
    "import pandas as pd",
    "import numpy as np",
    "import scipy.sparse as sp",
    "import sys",
    "",
    "# Compatibility patch for scanpy/igraph with newer SciPy sparse constructors.",
    "try:",
    "    from scanpy import _utils as sc_utils",
    "    def _smcanno_get_sparse_from_igraph(graph, weight_attr=None):",
    "        edges = graph.get_edgelist()",
    "        shape = (graph.vcount(), graph.vcount())",
    "        if len(edges) == 0:",
    "            return sp.csr_matrix(shape, dtype=np.float64)",
    "        rows = np.asarray([edge[0] for edge in edges], dtype=np.int64)",
    "        cols = np.asarray([edge[1] for edge in edges], dtype=np.int64)",
    "        if weight_attr is None:",
    "            weights = np.ones(len(edges), dtype=np.float64)",
    "        else:",
    "            weights = np.asarray(graph.es[weight_attr], dtype=np.float64)",
    "        return sp.csr_matrix((weights, (rows, cols)), shape=shape)",
    "    sc_utils.get_sparse_from_igraph = _smcanno_get_sparse_from_igraph",
    "except Exception:",
    "    pass",
    "",
    "cfg = json.loads(pathlib.Path(sys.argv[1]).read_text())",
    "out_dir = pathlib.Path(cfg['output_dir'])",
    "out_dir.mkdir(parents=True, exist_ok=True)",
    "adata = sc.read_h5ad(cfg['adata_h5ad_path'])",
    "groupby_key = cfg['groupby_key']",
    "if groupby_key not in adata.obs.columns:",
    "    raise KeyError(f\"Missing groupby key: {groupby_key}\")",
    "adata.obs[groupby_key] = pd.Categorical(adata.obs[groupby_key])",
    "use_rep = cfg.get('use_rep')",
    "if use_rep:",
    "    if use_rep not in adata.obsm:",
    "        raise KeyError(f\"Missing embedding: {use_rep}\")",
    "else:",
    "    if 'X_pca' not in adata.obsm:",
    "        sc.pp.pca(adata)",
    "    use_rep = 'X_pca'",
    "sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=int(cfg.get('n_neighbors', 30)), random_state=int(cfg.get('random_seed', 42)))",
    "sc.tl.paga(adata, groups=groupby_key)",
    "conn = adata.uns['paga']['connectivities'].toarray()",
    "cats = list(adata.obs[groupby_key].cat.categories)",
    "pd.DataFrame(conn, index=cats, columns=cats).to_csv(out_dir / 'paga_connectivities.csv')",
    "manifest = {",
    "    'groupby_key': groupby_key,",
    "    'use_rep': use_rep,",
    "    'n_groups': len(cats),",
    "    'connectivity_csv': str(out_dir / 'paga_connectivities.csv')",
    "}",
    "(out_dir / 'paga_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')",
    "print(json.dumps(manifest, ensure_ascii=False))"
  )
  pa_write_markdown(lines, path)
  Sys.chmod(path, mode = "0755")
  invisible(path)
}

pa_slingshot_seurat_to_sce <- function(seurat_obj, cluster_col, reduced_dim) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required for Slingshot runner.", call. = FALSE)
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Package 'SeuratObject' is required for Slingshot runner.", call. = FALSE)
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("Package 'SingleCellExperiment' is required for Slingshot runner.", call. = FALSE)
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("cluster_col not found in Seurat metadata: %s", cluster_col), call. = FALSE)
  }

  assay_name <- Seurat::DefaultAssay(seurat_obj)
  layer_name <- if ("data" %in% SeuratObject::Layers(seurat_obj[[assay_name]])) "data" else "counts"
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = layer_name)
  embed <- Seurat::Embeddings(seurat_obj, reduction = reduced_dim)
  common_cells <- intersect(colnames(seurat_obj), rownames(embed))
  if (length(common_cells) == 0L) stop("No shared cells between Seurat object and selected embedding", call. = FALSE)

  expr <- expr[, common_cells, drop = FALSE]
  embed <- embed[common_cells, , drop = FALSE]
  meta <- seurat_obj@meta.data[common_cells, , drop = FALSE]

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = expr))
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(meta)
  SingleCellExperiment::reducedDims(sce)$PA_SLINGSHOT <- as.matrix(embed)
  sce
}

pa_branch_summary_from_lineages <- function(lineages) {
  if (length(lineages) == 0L) return(data.frame(branch_id = character(), from = character(), to = character(), stringsAsFactors = FALSE))
  rows <- list()
  idx <- 1L
  for (i in seq_along(lineages)) {
    path <- as.character(lineages[[i]])
    if (length(path) < 2L) next
    for (j in seq_len(length(path) - 1L)) {
      rows[[idx]] <- data.frame(
        branch_id = sprintf("L%d_B%d", i, j),
        from = path[[j]],
        to = path[[j + 1L]],
        lineage = sprintf("Lineage_%d", i),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  if (length(rows) == 0L) return(data.frame(branch_id = character(), from = character(), to = character(), lineage = character(), stringsAsFactors = FALSE))
  unique(do.call(rbind, rows))
}

pa_summarize_cytotrace2_scores <- function(score_table,
                                           score_col,
                                           cell_id_col = "cell_id",
                                           pseudotime_table = NULL,
                                           min_abs_rho = 0.3) {
  pa_validate_required_columns(score_table, c(cell_id_col, score_col), "CytoTRACE2 score_table")
  score_df <- score_table[, c(cell_id_col, score_col), drop = FALSE]
  colnames(score_df) <- c("cell_id", "cytotrace2_score")
  score_df$cytotrace2_score <- suppressWarnings(as.numeric(score_df$cytotrace2_score))
  score_df <- score_df[is.finite(score_df$cytotrace2_score), , drop = FALSE]

  if (is.null(pseudotime_table) || !is.data.frame(pseudotime_table) || nrow(pseudotime_table) == 0L) {
    summary_tbl <- data.frame(
      lineage = "overall",
      n_cells = nrow(score_df),
      mean_score = mean(score_df$cytotrace2_score),
      median_score = stats::median(score_df$cytotrace2_score),
      rho_with_pseudotime = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(summary_table = summary_tbl, direction_consistency = NA_character_))
  }

  pt_df <- pseudotime_table
  pa_validate_required_columns(pt_df, c("cell_id", "pseudotime"), "pseudotime_table")
  if (!"lineage" %in% colnames(pt_df)) pt_df$lineage <- "overall"
  pt_df$pseudotime <- suppressWarnings(as.numeric(pt_df$pseudotime))
  merged <- merge(score_df, pt_df[, c("cell_id", "lineage", "pseudotime")], by = "cell_id")
  merged <- merged[is.finite(merged$pseudotime), , drop = FALSE]

  if (nrow(merged) == 0L) {
    summary_tbl <- data.frame(
      lineage = "overall",
      n_cells = nrow(score_df),
      mean_score = mean(score_df$cytotrace2_score),
      median_score = stats::median(score_df$cytotrace2_score),
      rho_with_pseudotime = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(summary_table = summary_tbl, direction_consistency = "insufficient_overlap"))
  }

  split_tbl <- split(merged, merged$lineage)
  summary_rows <- lapply(names(split_tbl), function(lineage_id) {
    df <- split_tbl[[lineage_id]]
    rho <- if (nrow(df) >= 3L) suppressWarnings(stats::cor(df$cytotrace2_score, df$pseudotime, method = "spearman", use = "pairwise.complete.obs")) else NA_real_
    data.frame(
      lineage = lineage_id,
      n_cells = nrow(df),
      mean_score = mean(df$cytotrace2_score),
      median_score = stats::median(df$cytotrace2_score),
      rho_with_pseudotime = rho,
      stringsAsFactors = FALSE
    )
  })
  summary_tbl <- do.call(rbind, summary_rows)
  rho_vals <- summary_tbl$rho_with_pseudotime[is.finite(summary_tbl$rho_with_pseudotime)]
  direction_consistency <- if (length(rho_vals) == 0L) {
    "insufficient_overlap"
  } else if (all(rho_vals <= -abs(min_abs_rho))) {
    "consistent"
  } else if (any(rho_vals >= abs(min_abs_rho))) {
    "inconsistent"
  } else {
    "weak"
  }
  list(summary_table = summary_tbl, direction_consistency = direction_consistency)
}

pa_import_paga_connectivities_csv <- function(path) {
  path <- pa_scalar_chr(path, "path")
  if (!file.exists(path)) stop(sprintf("PAGA connectivity file does not exist: %s", path), call. = FALSE)
  df <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  if (nrow(mat) != ncol(mat)) stop("PAGA connectivity matrix must be square", call. = FALSE)
  list(
    connectivity_matrix = mat,
    group_labels = rownames(mat),
    source_path = normalizePath(path, winslash = "/", mustWork = FALSE)
  )
}

pa_paga_edge_table <- function(connectivity_matrix, group_labels, min_edge_weight = 0) {
  idx <- which(upper.tri(connectivity_matrix) & connectivity_matrix > min_edge_weight, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE))
  }
  data.frame(
    from = group_labels[idx[, 1]],
    to = group_labels[idx[, 2]],
    weight = as.numeric(connectivity_matrix[idx]),
    stringsAsFactors = FALSE
  )
}

pa_build_paga_topology_packet <- function(
  connectivity_matrix,
  group_labels = NULL,
  source_path = NA_character_,
  group_key = NA_character_,
  min_edge_weight = 0.05
) {
  mat <- as.matrix(connectivity_matrix)
  storage.mode(mat) <- "numeric"
  if (nrow(mat) != ncol(mat)) stop("connectivity_matrix must be square", call. = FALSE)
  group_labels <- pa_null_coalesce(group_labels, rownames(mat))
  if (is.null(group_labels) || length(group_labels) != nrow(mat)) {
    stop("group_labels must have same length as connectivity_matrix dimensions", call. = FALSE)
  }
  rownames(mat) <- group_labels
  colnames(mat) <- group_labels
  edge_table <- pa_paga_edge_table(mat, group_labels, min_edge_weight = min_edge_weight)
  list(
    screen_type = "PAGA",
    group_key = pa_null_coalesce(group_key, NA_character_),
    group_labels = group_labels,
    n_groups = length(group_labels),
    n_edges = nrow(edge_table),
    connectivity_matrix = mat,
    edge_table = edge_table,
    source_path = pa_null_coalesce(source_path, NA_character_),
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_build_slingshot_trajectory_packet <- function(
  pseudotime_table = NULL,
  lineage_summary = NULL,
  branch_summary = NULL,
  source_path = NA_character_,
  root_state = NA_character_,
  terminal_states = NULL
) {
  list(
    engine = "Slingshot",
    pseudotime_table = pseudotime_table,
    lineage_summary = lineage_summary,
    branch_summary = branch_summary,
    root_state = pa_null_coalesce(root_state, NA_character_),
    terminal_states = pa_null_coalesce(terminal_states, character()),
    source_path = pa_null_coalesce(source_path, NA_character_),
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_build_monocle3_trajectory_packet <- function(
  pseudotime_table = NULL,
  lineage_summary = NULL,
  branch_summary = NULL,
  partition_summary = NULL,
  source_path = NA_character_,
  root_state = NA_character_,
  terminal_states = NULL
) {
  list(
    engine = "Monocle3",
    pseudotime_table = pseudotime_table,
    lineage_summary = lineage_summary,
    branch_summary = branch_summary,
    partition_summary = partition_summary,
    root_state = pa_null_coalesce(root_state, NA_character_),
    terminal_states = pa_null_coalesce(terminal_states, character()),
    source_path = pa_null_coalesce(source_path, NA_character_),
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_import_cytotrace2_scores_csv <- function(path, score_col = NULL) {
  path <- pa_scalar_chr(path, "path")
  if (!file.exists(path)) stop(sprintf("CytoTRACE2 score file does not exist: %s", path), call. = FALSE)
  df <- utils::read.csv(path, check.names = FALSE)
  guess_cols <- c("cytotrace2_score", "CytoTRACE2", "cytotrace2", "score")
  score_col <- pa_null_coalesce(score_col, guess_cols[guess_cols %in% colnames(df)][1])
  if (is.na(score_col) || !nzchar(score_col)) stop("Could not determine CytoTRACE2 score column", call. = FALSE)
  list(score_table = df, score_col = score_col, source_path = normalizePath(path, winslash = "/", mustWork = FALSE))
}

pa_build_cytotrace2_validation_packet <- function(
  score_table = NULL,
  maturity_summary = NULL,
  source_path = NA_character_,
  direction_consistency = NA_character_,
  score_col = "cytotrace2_score"
) {
  list(
    validator = "CytoTRACE2",
    score_table = score_table,
    score_col = pa_null_coalesce(score_col, "cytotrace2_score"),
    maturity_summary = maturity_summary,
    direction_consistency = pa_null_coalesce(direction_consistency, NA_character_),
    source_path = pa_null_coalesce(source_path, NA_character_),
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_build_trajectory_branch_packet <- function(
  topology_screen = NULL,
  primary_trajectory = NULL,
  maturity_validation = NULL,
  preferred_engine = "Slingshot"
) {
  consensus_status <- "insufficient"
  if (is.list(primary_trajectory) && !is.null(primary_trajectory$engine)) {
    consensus_status <- if (is.list(maturity_validation) && nzchar(pa_safe_trim(maturity_validation$direction_consistency))) {
      paste("primary+validation", pa_safe_trim(maturity_validation$direction_consistency), sep = ":")
    } else {
      "primary_only"
    }
  }
  list(
    topology_screen = topology_screen,
    primary_trajectory = primary_trajectory,
    maturity_validation = maturity_validation,
    preferred_engine = pa_scalar_chr(preferred_engine, "preferred_engine"),
    consensus_status = consensus_status,
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_trajectory_packet_summary_lines <- function(trajectory_branch_packet) {
  if (is.null(trajectory_branch_packet) || !is.list(trajectory_branch_packet)) {
    return("- 未提供 trajectory branch packet。")
  }
  top_type <- if (is.list(trajectory_branch_packet$topology_screen)) trajectory_branch_packet$topology_screen$screen_type else "无"
  primary_engine <- if (is.list(trajectory_branch_packet$primary_trajectory)) trajectory_branch_packet$primary_trajectory$engine else "无"
  validator <- if (is.list(trajectory_branch_packet$maturity_validation)) trajectory_branch_packet$maturity_validation$validator else "无"
  c(
    sprintf("- Topology screen：%s", pa_null_coalesce(top_type, "无")),
    sprintf("- Primary trajectory：%s", pa_null_coalesce(primary_engine, "无")),
    sprintf("- Maturity validation：%s", pa_null_coalesce(validator, "无")),
    sprintf("- Consensus status：%s", pa_null_coalesce(trajectory_branch_packet$consensus_status, "无"))
  )
}

pa_run_paga_runner <- function(adata_h5ad_path,
                               output_dir,
                               groupby_key,
                               python_cmd = NULL,
                               use_rep = NULL,
                               n_neighbors = 30L,
                               random_seed = 42L,
                               dry_run = FALSE) {
  adata_h5ad_path <- pa_scalar_chr(adata_h5ad_path, "adata_h5ad_path")
  output_dir <- pa_prepare_output_dir(output_dir)
  py_exec <- pa_resolve_python_executable(python_cmd)

  cfg <- list(
    adata_h5ad_path = normalizePath(adata_h5ad_path, winslash = "/", mustWork = FALSE),
    output_dir = output_dir,
    groupby_key = pa_scalar_chr(groupby_key, "groupby_key"),
    use_rep = pa_null_coalesce(use_rep, NULL),
    n_neighbors = as.integer(n_neighbors),
    random_seed = as.integer(random_seed)
  )
  cfg_path <- file.path(output_dir, "paga_runner_config.json")
  launcher_path <- file.path(output_dir, "run_paga_runner_launcher.py")
  pa_write_json(cfg, cfg_path)
  pa_write_paga_launcher(launcher_path)

  plan <- list(status = if (isTRUE(dry_run)) "dry_run" else "ready", python = py_exec, config_path = cfg_path, launcher_path = launcher_path)
  if (isTRUE(dry_run)) return(plan)

  cmd_res <- pa_run_system_command(py_exec, args = c(launcher_path, cfg_path), fail_on_error = FALSE)
  manifest <- pa_json_read(file.path(output_dir, "paga_manifest.json"), default = list())
  if (!file.exists(file.path(output_dir, "paga_connectivities.csv"))) {
    return(list(status = "error", command_result = cmd_res, manifest = manifest, topology_packet = NULL, plan = plan))
  }
  paga_import <- pa_import_paga_connectivities_csv(file.path(output_dir, "paga_connectivities.csv"))
  topology_packet <- pa_build_paga_topology_packet(
    connectivity_matrix = paga_import$connectivity_matrix,
    group_labels = paga_import$group_labels,
    source_path = paga_import$source_path,
    group_key = groupby_key
  )

  list(
    status = if (identical(cmd_res$status, 0L)) "ok" else "error",
    command_result = cmd_res,
    manifest = manifest,
    topology_packet = topology_packet,
    plan = plan
  )
}

pa_run_slingshot_runner <- function(seurat_obj,
                                    output_dir,
                                    cluster_col,
                                    reduced_dim = "umap",
                                    start_cluster = NULL,
                                    end_clusters = NULL) {
  if (!requireNamespace("slingshot", quietly = TRUE)) stop("Package 'slingshot' is required for Slingshot runner.", call. = FALSE)
  output_dir <- pa_prepare_output_dir(output_dir)

  sce <- pa_slingshot_seurat_to_sce(seurat_obj, cluster_col = cluster_col, reduced_dim = reduced_dim)
  SummarizedExperiment::colData(sce)[[cluster_col]] <- as.factor(as.character(SummarizedExperiment::colData(sce)[[cluster_col]]))

  sling_obj <- tryCatch(
    slingshot::slingshot(
      sce,
      clusterLabels = cluster_col,
      reducedDim = "PA_SLINGSHOT",
      start.clus = start_cluster,
      end.clus = end_clusters
    ),
    error = function(e) {
      if (!grepl("singular", conditionMessage(e), ignore.case = TRUE)) stop(e)
      rd <- SingleCellExperiment::reducedDims(sce)$PA_SLINGSHOT
      set.seed(42)
      rd <- rd + matrix(stats::rnorm(length(rd), sd = 1e-6), nrow = nrow(rd), ncol = ncol(rd))
      SingleCellExperiment::reducedDims(sce)$PA_SLINGSHOT <- rd
      slingshot::slingshot(
        sce,
        clusterLabels = cluster_col,
        reducedDim = "PA_SLINGSHOT",
        start.clus = start_cluster,
        end.clus = end_clusters
      )
    }
  )

  pseudotime_mat <- slingshot::slingPseudotime(sling_obj)
  weight_mat <- slingshot::slingCurveWeights(sling_obj)
  lineages <- slingshot::slingLineages(sling_obj)
  lineage_ids <- colnames(pseudotime_mat)
  if (is.null(lineage_ids)) lineage_ids <- sprintf("Lineage_%d", seq_len(ncol(pseudotime_mat)))

  pt_rows <- lapply(seq_len(ncol(pseudotime_mat)), function(i) {
    data.frame(
      cell_id = rownames(pseudotime_mat),
      lineage = lineage_ids[[i]],
      pseudotime = as.numeric(pseudotime_mat[, i]),
      curve_weight = as.numeric(weight_mat[, i]),
      stringsAsFactors = FALSE
    )
  })
  pseudotime_table <- do.call(rbind, pt_rows)
  pseudotime_table <- pseudotime_table[is.finite(pseudotime_table$pseudotime), , drop = FALSE]

  lineage_summary <- do.call(rbind, lapply(seq_along(lineages), function(i) {
    path <- as.character(lineages[[i]])
    data.frame(
      lineage = lineage_ids[[i]],
      n_cells = sum(is.finite(pseudotime_mat[, i])),
      root_state = path[[1]],
      terminal_state = path[[length(path)]],
      cluster_path = paste(path, collapse = " -> "),
      stringsAsFactors = FALSE
    )
  }))
  branch_summary <- pa_branch_summary_from_lineages(lineages)

  pseudotime_path <- file.path(output_dir, "slingshot_pseudotime.tsv")
  lineage_path <- file.path(output_dir, "slingshot_lineage_summary.tsv")
  branch_path <- file.path(output_dir, "slingshot_branch_summary.tsv")
  pa_write_tsv(pseudotime_table, pseudotime_path)
  pa_write_tsv(lineage_summary, lineage_path)
  pa_write_tsv(branch_summary, branch_path)
  pa_write_rds(sling_obj, file.path(output_dir, "slingshot_object.rds"))

  trajectory_packet <- pa_build_slingshot_trajectory_packet(
    pseudotime_table = pseudotime_table,
    lineage_summary = lineage_summary,
    branch_summary = branch_summary,
    source_path = output_dir,
    root_state = lineage_summary$root_state[[1]],
    terminal_states = lineage_summary$terminal_state
  )

  list(
    status = "ok",
    trajectory_packet = trajectory_packet,
    manifest = list(
      output_dir = output_dir,
      pseudotime_path = pseudotime_path,
      lineage_path = lineage_path,
      branch_path = branch_path
    )
  )
}

pa_run_cytotrace2_runner <- function(output_dir,
                                     score_csv = NULL,
                                     pseudotime_table = NULL,
                                     score_col = NULL,
                                     cell_id_col = "cell_id",
                                     command = NULL,
                                     command_args = character(),
                                     expected_score_csv = NULL,
                                     min_abs_rho = 0.3,
                                     dry_run = FALSE) {
  output_dir <- pa_prepare_output_dir(output_dir)
  if (is.null(expected_score_csv)) expected_score_csv <- file.path(output_dir, "cytotrace2_scores.csv")

  plan <- list(
    status = if (isTRUE(dry_run)) "dry_run" else "ready",
    output_dir = output_dir,
    score_csv = pa_null_coalesce(score_csv, NA_character_),
    command = pa_null_coalesce(command, NA_character_),
    expected_score_csv = expected_score_csv
  )
  if (isTRUE(dry_run)) return(plan)

  if (is.null(score_csv) || !nzchar(pa_safe_trim(score_csv))) {
    if (is.null(command) || !nzchar(pa_safe_trim(command))) {
      stop("CytoTRACE2 runner needs either score_csv or an external command to generate one.", call. = FALSE)
    }
    pa_run_system_command(command, args = command_args)
    score_csv <- expected_score_csv
  }

  imported <- pa_import_cytotrace2_scores_csv(score_csv, score_col = score_col)
  summary_info <- pa_summarize_cytotrace2_scores(
    score_table = imported$score_table,
    score_col = imported$score_col,
    cell_id_col = cell_id_col,
    pseudotime_table = pseudotime_table,
    min_abs_rho = min_abs_rho
  )

  summary_path <- file.path(output_dir, "cytotrace2_maturity_summary.tsv")
  pa_write_tsv(summary_info$summary_table, summary_path)
  validation_packet <- pa_build_cytotrace2_validation_packet(
    score_table = imported$score_table,
    maturity_summary = summary_info$summary_table,
    source_path = imported$source_path,
    direction_consistency = summary_info$direction_consistency,
    score_col = imported$score_col
  )

  list(
    status = "ok",
    validation_packet = validation_packet,
    manifest = list(
      output_dir = output_dir,
      score_csv = imported$source_path,
      maturity_summary_path = summary_path,
      direction_consistency = summary_info$direction_consistency
    )
  )
}

if (sys.nframe() == 0) {
  cat("Trajectory Branch Helper (2026-04-28 v1) loaded.\n")
}
