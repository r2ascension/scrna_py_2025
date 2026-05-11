#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

options(warn = 1)

current_base <- "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415"
helper_path  <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"
wrapper_path <- "/home/h2048/script/R/tissue_comparison_generic_wrapper_20260412.R"
engine_path  <- "/home/h2048/script/R/bcell_tissue_comparison_v2_6_1_20260410.R"

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
    output_dir = character(),
    stringsAsFactors = FALSE
  )
}

parse_gene_ratio <- function(x) {
  vals <- suppressWarnings(as.character(x))
  out <- rep(NA_real_, length(vals))
  has_ratio <- grepl("/", vals, fixed = TRUE)
  if (any(has_ratio)) {
    parts <- strsplit(vals[has_ratio], "/", fixed = TRUE)
    out[has_ratio] <- vapply(parts, function(p) {
      if (length(p) != 2) return(NA_real_)
      num <- suppressWarnings(as.numeric(p[[1]]))
      den <- suppressWarnings(as.numeric(p[[2]]))
      if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
      num / den
    }, numeric(1))
  }
  out
}

truncate_text_local <- function(x, width = 90L) {
  x <- trimws(as.character(x))
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3L), "..."), x)
}

collect_bubbleplot_df <- function(enrich_bundle, top_n_per_db = 4L, max_rows = 24L) {
  rows <- list()
  if (is.null(enrich_bundle) || length(enrich_bundle) == 0) return(data.frame())
  for (direction in names(enrich_bundle)) {
    db_list <- enrich_bundle[[direction]]
    if (is.null(db_list) || length(db_list) == 0) next
    for (db_name in names(db_list)) {
      er <- db_list[[db_name]]
      er_df <- tryCatch(as.data.frame(er), error = function(e) NULL)
      if (is.null(er_df) || nrow(er_df) == 0) next
      if (!"Description" %in% colnames(er_df)) er_df$Description <- rownames(er_df)
      if (!"Count" %in% colnames(er_df)) er_df$Count <- NA_integer_
      if (!"p.adjust" %in% colnames(er_df)) er_df$p.adjust <- NA_real_
      if (!"GeneRatio" %in% colnames(er_df)) er_df$GeneRatio <- NA_character_
      er_df <- er_df %>%
        mutate(
          direction = as.character(direction),
          db = as.character(db_name),
          Count = suppressWarnings(as.numeric(Count)),
          p.adjust = suppressWarnings(as.numeric(p.adjust)),
          gene_ratio_num = parse_gene_ratio(GeneRatio),
          Description = trimws(as.character(Description))
        ) %>%
        arrange(p.adjust, desc(Count), desc(gene_ratio_num)) %>%
        slice_head(n = top_n_per_db)
      rows[[length(rows) + 1L]] <- er_df
    }
  }
  if (length(rows) == 0) return(data.frame())
  out <- bind_rows(rows) %>%
    mutate(
      score = -log10(pmax(p.adjust, 1e-300)),
      bubble_size = ifelse(!is.na(Count) & Count > 0, Count, ifelse(!is.na(gene_ratio_num), pmax(1, round(100 * gene_ratio_num)), 1)),
      term_label = sprintf("[%s] %s", db, truncate_text_local(Description, width = 90L))
    ) %>%
    arrange(p.adjust, desc(bubble_size), term_label) %>%
    slice_head(n = max_rows)
  out
}

write_bubbleplot_overview <- function(enrich_bundle, out_dir, plot_title, save_plot_fn) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bubble_df <- collect_bubbleplot_df(enrich_bundle)
  fwrite(bubble_df, file.path(out_dir, "bubbleplot_overview.tsv"), sep = "\t")
  if (nrow(bubble_df) == 0) {
    p <- ggplot() +
      annotate("text", x = 1, y = 1, label = "No enriched pathways passed current thresholds", size = 5) +
      xlim(0, 2) + ylim(0, 2) +
      theme_void() +
      ggtitle(plot_title)
    save_plot_fn(p, file.path(out_dir, "bubbleplot_overview"), width = 10, height = 4)
    return(invisible(bubble_df))
  }
  bubble_df$term_label <- factor(bubble_df$term_label, levels = rev(unique(bubble_df$term_label)))
  p <- ggplot(bubble_df, aes(x = direction, y = term_label, size = bubble_size, color = score)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(name = "Gene count", range = c(2.5, 10)) +
    scale_color_gradient(name = "-log10(padj)", low = "#4C78A8", high = "#D62728") +
    labs(title = plot_title, x = "Direction", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  save_plot_fn(p, file.path(out_dir, "bubbleplot_overview"), width = 11, height = max(5, 0.28 * nrow(bubble_df) + 2.5))
  invisible(bubble_df)
}

cfg <- tc_build_generic_tissue_comparison_config(
  lineage = "BCELL",
  overrides = list(
    OUTPUT_DIR = current_base,
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0413/bcell_tissue_comparison_v2_6_5_c22drop_20260413",
    REUSE_PREVIOUS_FINAL_OBJECT = TRUE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = FALSE,
    GENERATED_BY_LABEL = "replay_bcell_l3_ofa_extension_20260419.R"
  )
)

loaded_env_files <- tc_load_env_candidates(cfg$ENV_FILE_CANDIDATES)
if (length(loaded_env_files) > 0) {
  message(sprintf("[OK] Loaded env files: %s", paste(loaded_env_files, collapse = ", ")))
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
    message("[OK] init-only B-cell bootstrap loaded")
  }
)

engine_fn <- function(name) get(name, envir = exec_env, inherits = TRUE)
safe_name <- engine_fn("safe_name")
save_plot <- engine_fn("save_plot")
resolve_dominant_label <- engine_fn("resolve_dominant_label")
tc_filter_gene_symbols <- engine_fn("tc_filter_gene_symbols")
tc_filter_term2gene_df <- engine_fn("tc_filter_term2gene_df")
get_enrichment_size_rule <- engine_fn("get_enrichment_size_rule")

CELLTYPE_L2_COL <- get("CELLTYPE_L2_COL", envir = exec_env, inherits = TRUE)
CELLTYPE_L3_COL <- get("CELLTYPE_L3_COL", envir = exec_env, inherits = TRUE)
LINEAGE_DISPLAY <- get("LINEAGE_DISPLAY", envir = exec_env, inherits = TRUE)
L3_OFA_MIN_CELLS_FOCAL <- get("L3_OFA_MIN_CELLS_FOCAL", envir = exec_env, inherits = TRUE)
L3_OFA_MIN_CELLS_REST <- get("L3_OFA_MIN_CELLS_REST", envir = exec_env, inherits = TRUE)
L3_OFA_PADJ_THR <- get("L3_OFA_PADJ_THR", envir = exec_env, inherits = TRUE)
L3_OFA_LFC_THR <- get("L3_OFA_LFC_THR", envir = exec_env, inherits = TRUE)
L3_OFA_TOP_N <- get("L3_OFA_TOP_N", envir = exec_env, inherits = TRUE)
L3_OFA_MAX_CELLS_PER_IDENT <- get("L3_OFA_MAX_CELLS_PER_IDENT", envir = exec_env, inherits = TRUE)

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
run_ofa_marker_enrichment <- engine_fn("run_ofa_marker_enrichment")

obj_path <- file.path(current_base, "bcell_tissue_comparison_final.rds")
if (!file.exists(obj_path)) stop(sprintf("Missing final object: %s", obj_path))
obj <- readRDS(obj_path)

meta <- obj@meta.data
if (!all(c(CELLTYPE_L2_COL, CELLTYPE_L3_COL) %in% colnames(meta))) {
  stop(sprintf("Missing required columns in final object metadata: %s", paste(setdiff(c(CELLTYPE_L2_COL, CELLTYPE_L3_COL), colnames(meta)), collapse = ", ")))
}

rpt_dir <- file.path(current_base, "reports")
l3_ofa_dir <- file.path(rpt_dir, "l3_ofa")
vs_rest_dir <- file.path(l3_ofa_dir, "vs_other_cell_types")
same_l2_dir <- file.path(l3_ofa_dir, "vs_same_l2_other_l3")
manifest_path <- file.path(l3_ofa_dir, "l3_ofa_extension_manifest.tsv")
summary_md_path <- file.path(current_base, "L3_OFA_EXTENSION_20260419.md")

safe_unlink(c(
  vs_rest_dir,
  same_l2_dir,
  file.path(l3_ofa_dir, "l3_ofa_vs_rest_all.rds"),
  file.path(l3_ofa_dir, "l3_ofa_same_l2_all.rds"),
  file.path(l3_ofa_dir, "l3_ofa_vs_rest_summary.tsv"),
  file.path(l3_ofa_dir, "l3_ofa_same_l2_summary.tsv"),
  manifest_path,
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

  message(sprintf("[RUN] L3 vs rest: %s", l3_name))
  rest_cells <- setdiff(colnames(obj), focal_cells)
  rest_out_dir <- file.path(vs_rest_dir, safe_name(l3_name))
  rest_res <- run_ofa_marker_enrichment(
    obj = obj,
    focal_cells = focal_cells,
    rest_cells = rest_cells,
    focal_label = l3_name,
    rest_label = "other_cell_types",
    out_dir = rest_out_dir,
    plot_title = sprintf("OFA: %s L3 %s vs rest", LINEAGE_DISPLAY, l3_name),
    group_col_name = "l3_ofa_group_ext",
    min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
    min_cells_rest = L3_OFA_MIN_CELLS_REST,
    padj_thr = L3_OFA_PADJ_THR,
    lfc_thr = L3_OFA_LFC_THR,
    top_n = L3_OFA_TOP_N,
    max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
  )
  if (!isTRUE(rest_res$skipped)) {
    write_bubbleplot_overview(
      enrich_bundle = rest_res$enrich,
      out_dir = rest_out_dir,
      plot_title = sprintf("Bubble plot: %s vs rest", l3_name),
      save_plot_fn = save_plot
    )
    l3_ofa_vs_rest_all[[l3_name]] <- c(rest_res, list(comparison = "vs_other_cell_types", celltype_l2 = focal_l2))
    vs_rest_rows[[length(vs_rest_rows) + 1L]] <- data.frame(
      family = "L3_vs_rest",
      celltype_label = l3_name,
      celltype_l2 = focal_l2 %||% "",
      comparison = "vs_other_cell_types",
      n_focal = rest_res$n_focal,
      n_rest = rest_res$n_rest,
      n_sig = rest_res$n_sig,
      output_dir = rest_out_dir,
      stringsAsFactors = FALSE
    )
    message(sprintf("[OK] L3 vs rest: %s | sig=%d", l3_name, rest_res$n_sig))
  } else {
    message(sprintf("[SKIP] L3 vs rest: %s | %s", l3_name, rest_res$reason))
  }

  if (is.na(focal_l2) || !nzchar(trimws(focal_l2))) {
    message(sprintf("[SKIP] L3 same-L2: %s | no dominant L2", l3_name))
    next
  }

  same_l2_cells <- rownames(meta)[
    as.character(meta[[CELLTYPE_L2_COL]]) == focal_l2 &
      as.character(meta[[CELLTYPE_L3_COL]]) != l3_name
  ]
  if (length(same_l2_cells) == 0) {
    message(sprintf("[SKIP] L3 same-L2: %s | no sibling L3 cells within %s", l3_name, focal_l2))
    next
  }

  message(sprintf("[RUN] L3 vs same-L2 siblings: %s | L2=%s", l3_name, focal_l2))
  same_l2_out_dir <- file.path(same_l2_dir, safe_name(l3_name))
  same_l2_res <- run_ofa_marker_enrichment(
    obj = obj,
    focal_cells = focal_cells,
    rest_cells = same_l2_cells,
    focal_label = l3_name,
    rest_label = paste0("other_", safe_name(focal_l2), "_L3"),
    out_dir = same_l2_out_dir,
    plot_title = sprintf("OFA: %s L3 %s vs same-L2 siblings (%s)", LINEAGE_DISPLAY, l3_name, focal_l2),
    group_col_name = "l3_same_l2_ofa_group_ext",
    min_cells_focal = L3_OFA_MIN_CELLS_FOCAL,
    min_cells_rest = L3_OFA_MIN_CELLS_REST,
    padj_thr = L3_OFA_PADJ_THR,
    lfc_thr = L3_OFA_LFC_THR,
    top_n = L3_OFA_TOP_N,
    max_cells_per_ident = L3_OFA_MAX_CELLS_PER_IDENT
  )
  if (!isTRUE(same_l2_res$skipped)) {
    write_bubbleplot_overview(
      enrich_bundle = same_l2_res$enrich,
      out_dir = same_l2_out_dir,
      plot_title = sprintf("Bubble plot: %s vs same-L2 siblings", l3_name),
      save_plot_fn = save_plot
    )
    l3_ofa_same_l2_all[[l3_name]] <- c(same_l2_res, list(comparison = "vs_same_l2_other_l3", celltype_l2 = focal_l2))
    same_l2_rows[[length(same_l2_rows) + 1L]] <- data.frame(
      family = "L3_vs_same_L2_other_L3",
      celltype_label = l3_name,
      celltype_l2 = focal_l2 %||% "",
      comparison = "vs_same_l2_other_l3",
      n_focal = same_l2_res$n_focal,
      n_rest = same_l2_res$n_rest,
      n_sig = same_l2_res$n_sig,
      output_dir = same_l2_out_dir,
      stringsAsFactors = FALSE
    )
    message(sprintf("[OK] L3 same-L2: %s | sig=%d", l3_name, same_l2_res$n_sig))
  } else {
    message(sprintf("[SKIP] L3 same-L2: %s | %s", l3_name, same_l2_res$reason))
  }

  gc()
}

vs_rest_summary_df <- if (length(vs_rest_rows) > 0) bind_rows(vs_rest_rows) else empty_l3_ofa_summary_df()
same_l2_summary_df <- if (length(same_l2_rows) > 0) bind_rows(same_l2_rows) else empty_l3_ofa_summary_df()
manifest_df <- data.frame(
  family = c("L3_vs_rest", "L3_vs_same_L2_other_L3"),
  n_comparisons = c(nrow(vs_rest_summary_df), nrow(same_l2_summary_df)),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)

saveRDS(l3_ofa_vs_rest_all, file.path(l3_ofa_dir, "l3_ofa_vs_rest_all.rds"))
saveRDS(l3_ofa_same_l2_all, file.path(l3_ofa_dir, "l3_ofa_same_l2_all.rds"))
fwrite(vs_rest_summary_df, file.path(l3_ofa_dir, "l3_ofa_vs_rest_summary.tsv"), sep = "\t")
fwrite(same_l2_summary_df, file.path(l3_ofa_dir, "l3_ofa_same_l2_summary.tsv"), sep = "\t")
fwrite(manifest_df, manifest_path, sep = "\t")

saveRDS(l3_ofa_vs_rest_all, file.path(rpt_dir, "l3_ofa_vs_rest_all.rds"))
saveRDS(l3_ofa_same_l2_all, file.path(rpt_dir, "l3_ofa_same_l2_all.rds"))
fwrite(vs_rest_summary_df, file.path(rpt_dir, "l3_ofa_vs_rest_summary.tsv"), sep = "\t")
fwrite(same_l2_summary_df, file.path(rpt_dir, "l3_ofa_same_l2_summary.tsv"), sep = "\t")

writeLines(
  c(
    "# B-cell L3 OFA Extension (2026-04-19)",
    "",
    sprintf("- Output root: `%s`", current_base),
    sprintf("- Final object: `%s`", obj_path),
    sprintf("- L3 vs rest comparisons: %d", nrow(vs_rest_summary_df)),
    sprintf("- L3 vs same-L2 sibling comparisons: %d", nrow(same_l2_summary_df)),
    "",
    "## Structured outputs",
    "",
    "- `reports/l3_ofa/l3_ofa_vs_rest_summary.tsv`",
    "- `reports/l3_ofa/l3_ofa_same_l2_summary.tsv`",
    "- `reports/l3_ofa/vs_other_cell_types/<L3>/bubbleplot_overview.(pdf|png)`",
    "- `reports/l3_ofa/vs_same_l2_other_l3/<L3>/bubbleplot_overview.(pdf|png)`",
    ""
  ),
  summary_md_path
)

message(sprintf("[OK] B-cell L3 OFA extension complete | vs_rest=%d | same_L2=%d", nrow(vs_rest_summary_df), nrow(same_l2_summary_df)))
print(manifest_df)
