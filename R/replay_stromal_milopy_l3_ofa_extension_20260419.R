#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

options(warn = 1)

helper_path <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
wrapper_path <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414.R"
engine_path <- "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R"
master_manifest_path <- "/home/h2048/data/R/0414/stromal_milopy_l3_ofa_extension_manifest_20260419.tsv"

jobs <- data.frame(
  lineage = c("STROMAL_ENDOTHELIAL", "STROMAL_FIBROBLAST", "STROMAL_SMC"),
  output_dir = c(
    "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
    "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414",
    "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
  ),
  stringsAsFactors = FALSE
)

if (!file.exists(helper_path)) stop(sprintf("Helper not found: %s", helper_path))
if (!file.exists(wrapper_path)) stop(sprintf("Wrapper not found: %s", wrapper_path))
if (!file.exists(engine_path)) stop(sprintf("Engine not found: %s", engine_path))

source(helper_path)
source(wrapper_path)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) return(y)
  x
}

safe_unlink <- function(paths) {
  for (path in unique(paths)) {
    if (dir.exists(path) || file.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  }
}

empty_l3_ofa_summary_df <- function() {
  data.frame(
    family = character(),
    celltype_label = character(),
    celltype_l2 = character(),
    comparison = character(),
    n_focal = integer(),
    n_rest = integer(),
    n_sig = integer(),
    bubble_png = character(),
    bubble_pdf = character(),
    output_dir = character(),
    stringsAsFactors = FALSE
  )
}

bootstrap_exec_env <- function(lineage, output_dir, generated_by_label) {
  cfg <- tc_build_generic_tissue_comparison_config(
    lineage = lineage,
    overrides = list(
      OUTPUT_DIR = output_dir,
      REUSE_PREVIOUS_FINAL_OBJECT = TRUE,
      REUSE_PREVIOUS_OUTPUT_SUMMARY = FALSE,
      GENERATED_BY_LABEL = generated_by_label
    )
  )

  loaded_env_files <- tc_load_env_candidates(cfg$ENV_FILE_CANDIDATES)
  if (length(loaded_env_files) > 0) {
    message(sprintf("[OK] %s loaded env files: %s", lineage, paste(loaded_env_files, collapse = ", ")))
  }

  exec_env <- new.env(parent = globalenv())
  tc_apply_named_list(cfg, envir = exec_env)
  invisible(tc_apply_advanced_shared_overrides(envir = exec_env, overrides = cfg$SHARED_OVERRIDES))
  assign("INITIALIZE_ONLY", TRUE, envir = exec_env)
  assign("SKIP_PREVIOUS_RUN_SUMMARY_IN_ENGINE", TRUE, envir = exec_env)
  assign("ADVANCED_HELPER_ALREADY_LOADED", TRUE, envir = exec_env)
  assign("SKIP_ENV_AUTOLOAD", TRUE, envir = exec_env)
  assign("REUSE_PREVIOUS_OUTPUT_SUMMARY", FALSE, envir = exec_env)
  assign("REUSE_PREVIOUS_FINAL_OBJECT", TRUE, envir = exec_env)
  assign("PIPELINE_TEST_MODE", FALSE, envir = exec_env)

  tryCatch(
    source(engine_path, local = exec_env),
    error = function(e) {
      if (!inherits(e, "bcell_pipeline_init_only")) stop(e)
      message(sprintf("[OK] %s init-only bootstrap loaded", lineage))
    }
  )

  list(cfg = cfg, exec_env = exec_env, loaded_env_files = loaded_env_files)
}

install_replay_safe_enrichment <- function(exec_env) {
  engine_fn <- function(name) get(name, envir = exec_env, inherits = TRUE)
  tc_filter_gene_symbols <- engine_fn("tc_filter_gene_symbols")
  tc_filter_term2gene_df <- engine_fn("tc_filter_term2gene_df")
  get_enrichment_size_rule <- engine_fn("get_enrichment_size_rule")

  run_gmt_enrichment_replay_safe <- function(gene_list, t2g, db_name, tested_genes = NULL) {
    if (is.null(t2g) || !is.data.frame(t2g) || nrow(t2g) == 0 || length(gene_list) < 5) return(NULL)
    t2g_use <- data.frame(
      term = as.character(t2g$term),
      gene = toupper(as.character(t2g$gene)),
      stringsAsFactors = FALSE
    )
    t2g_use <- t2g_use[!is.na(t2g_use$term) & !is.na(t2g_use$gene), , drop = FALSE]
    t2g_use <- t2g_use[trimws(t2g_use$term) != "" & trimws(t2g_use$gene) != "", , drop = FALSE]
    t2g_use <- tc_filter_term2gene_df(t2g_use)
    if (!is.data.frame(t2g_use) || nrow(t2g_use) == 0) return(NULL)

    universe_use <- if (!is.null(tested_genes)) {
      unique(intersect(tc_filter_gene_symbols(tested_genes), unique(t2g_use$gene)))
    } else {
      unique(t2g_use$gene)
    }
    if (length(universe_use) < 5) {
      cat(sprintf("    [INFO] %s: too few universe genes (%d)\n", db_name, length(universe_use)))
      return(NULL)
    }

    t2g_use <- t2g_use[t2g_use$gene %in% universe_use, , drop = FALSE]
    if (nrow(t2g_use) == 0) return(NULL)

    size_rule <- get_enrichment_size_rule(db_name)
    gs_sizes <- as.data.frame(table(t2g_use$term), stringsAsFactors = FALSE)
    colnames(gs_sizes) <- c("term", "gs_size")
    gs_sizes$gs_size <- suppressWarnings(as.integer(gs_sizes$gs_size))
    max_gs <- suppressWarnings(max(gs_sizes$gs_size, na.rm = TRUE))
    if (!is.finite(max_gs) || max_gs < 2) return(NULL)

    min_gs <- max(2L, min(size_rule$min, as.integer(max_gs)))
    max_gs_ok <- max(min_gs, size_rule$max)
    valid_terms <- gs_sizes$term[gs_sizes$gs_size >= min_gs & gs_sizes$gs_size <= max_gs_ok]
    if (length(valid_terms) == 0) return(NULL)

    t2g_use <- t2g_use[t2g_use$term %in% valid_terms, , drop = FALSE]
    genes_use <- intersect(tc_filter_gene_symbols(gene_list), unique(t2g_use$gene))
    if (length(genes_use) < 5) return(NULL)

    tryCatch(
      suppressMessages(enricher(
        gene = genes_use,
        TERM2GENE = t2g_use,
        universe = unique(t2g_use$gene),
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        pAdjustMethod = "BH",
        minGSSize = min_gs,
        maxGSSize = max_gs_ok
      )),
      error = function(e) {
        cat(sprintf("    [WARN] %s: %s\n", db_name, e$message))
        NULL
      }
    )
  }

  assign("run_gmt_enrichment", run_gmt_enrichment_replay_safe, envir = exec_env)
}

run_one_job <- function(lineage, output_dir) {
  generated_by_label <- "replay_stromal_milopy_l3_ofa_extension_20260419.R"
  boot <- bootstrap_exec_env(lineage = lineage, output_dir = output_dir, generated_by_label = generated_by_label)
  exec_env <- boot$exec_env
  install_replay_safe_enrichment(exec_env)

  engine_fn <- function(name) get(name, envir = exec_env, inherits = TRUE)
  safe_name <- engine_fn("safe_name")
  resolve_dominant_label <- engine_fn("resolve_dominant_label")
  run_ofa_marker_enrichment <- engine_fn("run_ofa_marker_enrichment")
  tc_run_tissue_milopy <- engine_fn("tc_run_tissue_milopy")

  FINAL_FILE_PREFIX <- get("FINAL_FILE_PREFIX", envir = exec_env, inherits = TRUE)
  LINEAGE_DISPLAY <- get("LINEAGE_DISPLAY", envir = exec_env, inherits = TRUE)
  CELLTYPE_L2_COL <- get("CELLTYPE_L2_COL", envir = exec_env, inherits = TRUE)
  CELLTYPE_L3_COL <- get("CELLTYPE_L3_COL", envir = exec_env, inherits = TRUE)
  SAMPLE_COL <- get("SAMPLE_COL", envir = exec_env, inherits = TRUE)
  TISSUE_COL <- get("TISSUE_COL", envir = exec_env, inherits = TRUE)

  L3_OFA_MIN_CELLS_FOCAL <- get("L3_OFA_MIN_CELLS_FOCAL", envir = exec_env, inherits = TRUE)
  L3_OFA_MIN_CELLS_REST <- get("L3_OFA_MIN_CELLS_REST", envir = exec_env, inherits = TRUE)
  L3_OFA_PADJ_THR <- get("L3_OFA_PADJ_THR", envir = exec_env, inherits = TRUE)
  L3_OFA_LFC_THR <- get("L3_OFA_LFC_THR", envir = exec_env, inherits = TRUE)
  L3_OFA_TOP_N <- get("L3_OFA_TOP_N", envir = exec_env, inherits = TRUE)
  L3_OFA_MAX_CELLS_PER_IDENT <- get("L3_OFA_MAX_CELLS_PER_IDENT", envir = exec_env, inherits = TRUE)

  MILOPY_OUTPUT_DIRNAME <- get("MILOPY_OUTPUT_DIRNAME", envir = exec_env, inherits = TRUE)
  MILOPY_LEVEL_COLS <- get("MILOPY_LEVEL_COLS", envir = exec_env, inherits = TRUE)
  MILOPY_CONDA_EXE <- get("MILOPY_CONDA_EXE", envir = exec_env, inherits = TRUE)
  MILOPY_ENV_PREFIX <- get("MILOPY_ENV_PREFIX", envir = exec_env, inherits = TRUE)
  MILOPY_SCRIPT_PATH <- get("MILOPY_SCRIPT_PATH", envir = exec_env, inherits = TRUE)
  MILOPY_LATENT_KEY_CANDIDATES <- get("MILOPY_LATENT_KEY_CANDIDATES", envir = exec_env, inherits = TRUE)
  MILOPY_UMAP_KEY_CANDIDATES <- get("MILOPY_UMAP_KEY_CANDIDATES", envir = exec_env, inherits = TRUE)
  MILOPY_MIN_CELLS_PER_CELLTYPE <- get("MILOPY_MIN_CELLS_PER_CELLTYPE", envir = exec_env, inherits = TRUE)
  MILOPY_MIN_CELLS_PER_SAMPLE <- get("MILOPY_MIN_CELLS_PER_SAMPLE", envir = exec_env, inherits = TRUE)
  MILOPY_MIN_SAMPLES_PER_TISSUE <- get("MILOPY_MIN_SAMPLES_PER_TISSUE", envir = exec_env, inherits = TRUE)
  MILOPY_N_NEIGHBORS <- get("MILOPY_N_NEIGHBORS", envir = exec_env, inherits = TRUE)
  MILOPY_NHOOD_PROP <- get("MILOPY_NHOOD_PROP", envir = exec_env, inherits = TRUE)
  MILOPY_ALPHA <- get("MILOPY_ALPHA", envir = exec_env, inherits = TRUE)
  MILOPY_RANDOM_SEED <- get("MILOPY_RANDOM_SEED", envir = exec_env, inherits = TRUE)
  MILOPY_MAX_CELLTYPES_PER_LEVEL <- get("MILOPY_MAX_CELLTYPES_PER_LEVEL", envir = exec_env, inherits = TRUE)
  MILOPY_MAKE_PLOTS <- get("MILOPY_MAKE_PLOTS", envir = exec_env, inherits = TRUE)
  MILOPY_WRITE_MILO_H5AD <- get("MILOPY_WRITE_MILO_H5AD", envir = exec_env, inherits = TRUE)
  MILOPY_LOG_FILENAME <- get("MILOPY_LOG_FILENAME", envir = exec_env, inherits = TRUE)
  MILOPY_STOP_ON_ERROR <- get("MILOPY_STOP_ON_ERROR", envir = exec_env, inherits = TRUE)

  obj_path <- file.path(output_dir, paste0(FINAL_FILE_PREFIX, ".rds"))
  h5ad_path <- file.path(output_dir, paste0(FINAL_FILE_PREFIX, ".h5ad"))
  rpt_dir <- file.path(output_dir, "reports")
  milopy_dir <- file.path(output_dir, MILOPY_OUTPUT_DIRNAME)
  l3_ofa_dir <- file.path(rpt_dir, "l3_ofa")
  vs_rest_dir <- file.path(l3_ofa_dir, "vs_other_cell_types")
  same_l2_dir <- file.path(l3_ofa_dir, "vs_same_l2_other_l3")
  summary_md_path <- file.path(output_dir, "POSTHOC_MILOPY_L3_OFA_EXTENSION_20260419.md")
  summary_tsv_path <- file.path(output_dir, "POSTHOC_MILOPY_L3_OFA_EXTENSION_20260419.tsv")

  if (!file.exists(obj_path)) stop(sprintf("[%s] Missing final object: %s", lineage, obj_path))
  if (!file.exists(h5ad_path)) stop(sprintf("[%s] Missing final h5ad: %s", lineage, h5ad_path))

  dir.create(rpt_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(obj_path)
  meta <- obj@meta.data
  if (!all(c(CELLTYPE_L2_COL, CELLTYPE_L3_COL, SAMPLE_COL, TISSUE_COL) %in% colnames(meta))) {
    missing_cols <- setdiff(c(CELLTYPE_L2_COL, CELLTYPE_L3_COL, SAMPLE_COL, TISSUE_COL), colnames(meta))
    stop(sprintf("[%s] Missing required metadata columns: %s", lineage, paste(missing_cols, collapse = ", ")))
  }

  message(sprintf("\n=== [%s] MiloPy L2/L3 tissue comparison ===", lineage))
  safe_unlink(c(
    milopy_dir,
    file.path(rpt_dir, "milopy_run.rds"),
    file.path(rpt_dir, "milopy_level_summary.tsv"),
    file.path(rpt_dir, "milopy_pairwise_summary.tsv"),
    file.path(rpt_dir, "milopy_contrast_summary.tsv")
  ))

  milopy_run <- tc_run_tissue_milopy(
    input_h5ad = h5ad_path,
    output_dir = milopy_dir,
    celltype_cols = MILOPY_LEVEL_COLS,
    sample_col = SAMPLE_COL,
    tissue_col = TISSUE_COL,
    conda_exe = MILOPY_CONDA_EXE,
    env_prefix = MILOPY_ENV_PREFIX,
    script_path = MILOPY_SCRIPT_PATH,
    latent_key_candidates = MILOPY_LATENT_KEY_CANDIDATES,
    umap_key_candidates = MILOPY_UMAP_KEY_CANDIDATES,
    min_cells_per_celltype = MILOPY_MIN_CELLS_PER_CELLTYPE,
    min_cells_per_sample = MILOPY_MIN_CELLS_PER_SAMPLE,
    min_samples_per_tissue = MILOPY_MIN_SAMPLES_PER_TISSUE,
    n_neighbors = MILOPY_N_NEIGHBORS,
    nhood_prop = MILOPY_NHOOD_PROP,
    alpha = MILOPY_ALPHA,
    random_seed = MILOPY_RANDOM_SEED,
    max_celltypes_per_level = MILOPY_MAX_CELLTYPES_PER_LEVEL,
    make_plots = MILOPY_MAKE_PLOTS,
    write_milo_h5ad = MILOPY_WRITE_MILO_H5AD,
    log_path = file.path(milopy_dir, MILOPY_LOG_FILENAME),
    stop_on_error = MILOPY_STOP_ON_ERROR
  )
  saveRDS(milopy_run, file.path(rpt_dir, "milopy_run.rds"))
  if (is.data.frame(milopy_run$level_summary) && nrow(milopy_run$level_summary) > 0) {
    fwrite(milopy_run$level_summary, file.path(rpt_dir, "milopy_level_summary.tsv"), sep = "\t")
  }
  if (is.data.frame(milopy_run$pairwise_summary) && nrow(milopy_run$pairwise_summary) > 0) {
    fwrite(milopy_run$pairwise_summary, file.path(rpt_dir, "milopy_pairwise_summary.tsv"), sep = "\t")
  }
  if (is.data.frame(milopy_run$contrast_summary) && nrow(milopy_run$contrast_summary) > 0) {
    fwrite(milopy_run$contrast_summary, file.path(rpt_dir, "milopy_contrast_summary.tsv"), sep = "\t")
  }

  message(sprintf("\n=== [%s] L3 OFA replay (vs rest + same-L2 siblings) ===", lineage))
  safe_unlink(c(
    vs_rest_dir,
    same_l2_dir,
    file.path(l3_ofa_dir, "l3_ofa_vs_rest_all.rds"),
    file.path(l3_ofa_dir, "l3_ofa_same_l2_all.rds"),
    file.path(l3_ofa_dir, "l3_ofa_vs_rest_summary.tsv"),
    file.path(l3_ofa_dir, "l3_ofa_same_l2_summary.tsv"),
    file.path(l3_ofa_dir, "l3_ofa_extension_manifest.tsv"),
    file.path(rpt_dir, "l3_ofa_vs_rest_all.rds"),
    file.path(rpt_dir, "l3_ofa_same_l2_all.rds"),
    file.path(rpt_dir, "l3_ofa_vs_rest_summary.tsv"),
    file.path(rpt_dir, "l3_ofa_same_l2_summary.tsv")
  ))

  dir.create(vs_rest_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(same_l2_dir, recursive = TRUE, showWarnings = FALSE)

  l3_types <- sort(unique(na.omit(as.character(meta[[CELLTYPE_L3_COL]]))))
  l3_types <- l3_types[nzchar(trimws(l3_types))]

  l3_ofa_vs_rest_all <- list()
  l3_ofa_same_l2_all <- list()
  vs_rest_rows <- list()
  same_l2_rows <- list()

  for (l3_name in l3_types) {
    focal_cells <- rownames(meta)[as.character(meta[[CELLTYPE_L3_COL]]) == l3_name]
    if (length(focal_cells) == 0) next
    focal_l2 <- resolve_dominant_label(meta[focal_cells, CELLTYPE_L2_COL, drop = TRUE])

    message(sprintf("[RUN] [%s] L3 vs rest: %s", lineage, l3_name))
    rest_cells <- setdiff(colnames(obj), focal_cells)
    rest_out_dir <- file.path(vs_rest_dir, safe_name(l3_name))
    rest_res <- run_ofa_marker_enrichment(
      obj = obj,
      focal_cells = focal_cells,
      rest_cells = rest_cells,
      focal_label = l3_name,
      rest_label = "other_cell_types",
      out_dir = rest_out_dir,
      plot_title = sprintf("OFA: %s L3 %s vs other cell types", LINEAGE_DISPLAY, l3_name),
      group_col_name = "l3_ofa_group_replay",
      min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
      min_cells_rest = L3_OFA_MIN_CELLS_REST,
      padj_thr = L3_OFA_PADJ_THR,
      lfc_thr = L3_OFA_LFC_THR,
      top_n = L3_OFA_TOP_N,
      max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
    )
    if (!isTRUE(rest_res$skipped)) {
      l3_ofa_vs_rest_all[[l3_name]] <- c(rest_res, list(comparison = "vs_other_cell_types", celltype_l2 = focal_l2))
      vs_rest_rows[[length(vs_rest_rows) + 1L]] <- data.frame(
        family = "L3_vs_other_cell_types",
        celltype_label = l3_name,
        celltype_l2 = focal_l2 %||% "",
        comparison = "vs_other_cell_types",
        n_focal = rest_res$n_focal,
        n_rest = rest_res$n_rest,
        n_sig = rest_res$n_sig,
        bubble_png = if (!is.null(rest_res$bubble_paths$png)) rest_res$bubble_paths$png else NA_character_,
        bubble_pdf = if (!is.null(rest_res$bubble_paths$pdf)) rest_res$bubble_paths$pdf else NA_character_,
        output_dir = rest_out_dir,
        stringsAsFactors = FALSE
      )
      message(sprintf("[OK] [%s] L3 vs rest: %s | sig=%d", lineage, l3_name, rest_res$n_sig))
    } else {
      message(sprintf("[SKIP] [%s] L3 vs rest: %s | %s", lineage, l3_name, rest_res$reason %||% "unknown"))
    }

    if (is.na(focal_l2) || !nzchar(trimws(as.character(focal_l2)))) {
      message(sprintf("[SKIP] [%s] L3 same-L2: %s | no dominant L2", lineage, l3_name))
      next
    }

    same_l2_cells <- rownames(meta)[
      as.character(meta[[CELLTYPE_L2_COL]]) == as.character(focal_l2) &
        as.character(meta[[CELLTYPE_L3_COL]]) != l3_name
    ]
    if (length(same_l2_cells) == 0) {
      message(sprintf("[SKIP] [%s] L3 same-L2: %s | no sibling L3 cells within %s", lineage, l3_name, focal_l2))
      next
    }

    message(sprintf("[RUN] [%s] L3 vs same-L2 siblings: %s | L2=%s", lineage, l3_name, focal_l2))
    same_l2_out_dir <- file.path(same_l2_dir, safe_name(l3_name))
    same_l2_res <- run_ofa_marker_enrichment(
      obj = obj,
      focal_cells = focal_cells,
      rest_cells = same_l2_cells,
      focal_label = l3_name,
      rest_label = paste0("other_", safe_name(focal_l2), "_L3"),
      out_dir = same_l2_out_dir,
      plot_title = sprintf("OFA: %s L3 %s vs same-L2 siblings (%s)", LINEAGE_DISPLAY, l3_name, focal_l2),
      group_col_name = "l3_same_l2_ofa_group_replay",
      min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
      min_cells_rest = L3_OFA_MIN_CELLS_REST,
      padj_thr = L3_OFA_PADJ_THR,
      lfc_thr = L3_OFA_LFC_THR,
      top_n = L3_OFA_TOP_N,
      max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
    )
    if (!isTRUE(same_l2_res$skipped)) {
      l3_ofa_same_l2_all[[l3_name]] <- c(same_l2_res, list(comparison = "vs_same_l2_other_l3", celltype_l2 = focal_l2))
      same_l2_rows[[length(same_l2_rows) + 1L]] <- data.frame(
        family = "L3_vs_same_L2_other_L3",
        celltype_label = l3_name,
        celltype_l2 = focal_l2 %||% "",
        comparison = "vs_same_l2_other_l3",
        n_focal = same_l2_res$n_focal,
        n_rest = same_l2_res$n_rest,
        n_sig = same_l2_res$n_sig,
        bubble_png = if (!is.null(same_l2_res$bubble_paths$png)) same_l2_res$bubble_paths$png else NA_character_,
        bubble_pdf = if (!is.null(same_l2_res$bubble_paths$pdf)) same_l2_res$bubble_paths$pdf else NA_character_,
        output_dir = same_l2_out_dir,
        stringsAsFactors = FALSE
      )
      message(sprintf("[OK] [%s] L3 same-L2: %s | sig=%d", lineage, l3_name, same_l2_res$n_sig))
    } else {
      message(sprintf("[SKIP] [%s] L3 same-L2: %s | %s", lineage, l3_name, same_l2_res$reason %||% "unknown"))
    }

    gc()
  }

  vs_rest_summary_df <- if (length(vs_rest_rows) > 0) bind_rows(vs_rest_rows) else empty_l3_ofa_summary_df()
  same_l2_summary_df <- if (length(same_l2_rows) > 0) bind_rows(same_l2_rows) else empty_l3_ofa_summary_df()
  l3_manifest_df <- data.frame(
    family = c("L3_vs_other_cell_types", "L3_vs_same_L2_other_L3"),
    n_comparisons = c(nrow(vs_rest_summary_df), nrow(same_l2_summary_df)),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )

  saveRDS(l3_ofa_vs_rest_all, file.path(l3_ofa_dir, "l3_ofa_vs_rest_all.rds"))
  saveRDS(l3_ofa_same_l2_all, file.path(l3_ofa_dir, "l3_ofa_same_l2_all.rds"))
  fwrite(vs_rest_summary_df, file.path(l3_ofa_dir, "l3_ofa_vs_rest_summary.tsv"), sep = "\t")
  fwrite(same_l2_summary_df, file.path(l3_ofa_dir, "l3_ofa_same_l2_summary.tsv"), sep = "\t")
  fwrite(l3_manifest_df, file.path(l3_ofa_dir, "l3_ofa_extension_manifest.tsv"), sep = "\t")

  saveRDS(l3_ofa_vs_rest_all, file.path(rpt_dir, "l3_ofa_vs_rest_all.rds"))
  saveRDS(l3_ofa_same_l2_all, file.path(rpt_dir, "l3_ofa_same_l2_all.rds"))
  fwrite(vs_rest_summary_df, file.path(rpt_dir, "l3_ofa_vs_rest_summary.tsv"), sep = "\t")
  fwrite(same_l2_summary_df, file.path(rpt_dir, "l3_ofa_same_l2_summary.tsv"), sep = "\t")

  milopy_levels_done <- if (is.data.frame(milopy_run$level_summary)) nrow(milopy_run$level_summary) else 0L
  milopy_pairwise_done <- if (is.data.frame(milopy_run$pairwise_summary)) nrow(milopy_run$pairwise_summary) else 0L
  milopy_contrast_done <- if (is.data.frame(milopy_run$contrast_summary)) nrow(milopy_run$contrast_summary) else 0L

  per_job_summary <- data.frame(
    lineage = lineage,
    lineage_display = LINEAGE_DISPLAY,
    output_dir = output_dir,
    obj_path = obj_path,
    h5ad_path = h5ad_path,
    milopy_status = as.character(milopy_run$status %||% "unknown"),
    milopy_levels = milopy_levels_done,
    milopy_pairwise = milopy_pairwise_done,
    milopy_contrast = milopy_contrast_done,
    l3_total = length(l3_types),
    l3_ofa_vs_rest = nrow(vs_rest_summary_df),
    l3_ofa_same_l2 = nrow(same_l2_summary_df),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  fwrite(per_job_summary, summary_tsv_path, sep = "\t")

  writeLines(
    c(
      sprintf("# %s post-hoc MiloPy + L3 OFA extension (2026-04-19)", LINEAGE_DISPLAY),
      "",
      sprintf("- Lineage preset: `%s`", lineage),
      sprintf("- Output root: `%s`", output_dir),
      sprintf("- Final object: `%s`", obj_path),
      sprintf("- Final h5ad: `%s`", h5ad_path),
      sprintf("- MiloPy status: `%s`", milopy_run$status %||% "unknown"),
      sprintf("- MiloPy level rows: %d", milopy_levels_done),
      sprintf("- MiloPy pairwise rows: %d", milopy_pairwise_done),
      sprintf("- MiloPy contrast rows: %d", milopy_contrast_done),
      sprintf("- L3 vs rest comparisons: %d", nrow(vs_rest_summary_df)),
      sprintf("- L3 vs same-L2 sibling comparisons: %d", nrow(same_l2_summary_df)),
      "",
      "## Structured outputs",
      "",
      sprintf("- `%s/`", file.path(output_dir, MILOPY_OUTPUT_DIRNAME)),
      sprintf("- `%s`", file.path(rpt_dir, "milopy_level_summary.tsv")),
      sprintf("- `%s`", file.path(rpt_dir, "milopy_pairwise_summary.tsv")),
      sprintf("- `%s`", file.path(rpt_dir, "milopy_contrast_summary.tsv")),
      sprintf("- `%s`", file.path(rpt_dir, "l3_ofa_vs_rest_summary.tsv")),
      sprintf("- `%s`", file.path(rpt_dir, "l3_ofa_same_l2_summary.tsv")),
      sprintf("- `%s/<L3>/bubbleplot_overview.(png|pdf)`", vs_rest_dir),
      sprintf("- `%s/<L3>/bubbleplot_overview.(png|pdf)`", same_l2_dir),
      ""
    ),
    summary_md_path
  )

  message(sprintf(
    "[OK] [%s] extension complete | Milo=%s | vs_rest=%d | same_L2=%d",
    lineage,
    milopy_run$status %||% "unknown",
    nrow(vs_rest_summary_df),
    nrow(same_l2_summary_df)
  ))

  per_job_summary
}

all_summaries <- bind_rows(lapply(seq_len(nrow(jobs)), function(i) {
  run_one_job(jobs$lineage[[i]], jobs$output_dir[[i]])
}))

dir.create(dirname(master_manifest_path), recursive = TRUE, showWarnings = FALSE)
fwrite(all_summaries, master_manifest_path, sep = "\t")

message("\n[OK] Stromal MiloPy + L3 OFA extension finished for all jobs")
print(all_summaries)
