#!/usr/bin/env Rscript

source('/home/h2048/script/R/modules/hdwgcna/programs/hdwgcna_tissue_preservation_runner_20260531.R')

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

set.seed(1)
mat <- matrix(rpois(8 * 12, lambda = 5), nrow = 8, dimnames = list(paste0('g', 1:8), paste0('c', 1:12)))
obj <- CreateSeuratObject(counts = mat)
obj$cell_type_final_l3 <- rep(c('Naive_B', 'Memory_B'), each = 6)
obj$tissue <- rep(c('nose', 'sinus', 'bronchial'), each = 4)
obj$sample <- rep(c('S1', 'S2'), times = 6)

grid <- prepare_tissue_grid(
  seurat_obj = obj,
  celltype_col = 'cell_type_final_l3',
  tissue_col = 'tissue',
  sample_col = 'sample',
  min_cells_per_tissue = 1L,
  min_samples_per_tissue = 1L
)

assert_true(is.data.frame(grid), 'prepare_tissue_grid should return a data.frame')
assert_true('cells' %in% colnames(grid), 'prepare_tissue_grid should keep the cells list-column')
assert_true(length(grid$cells[[1]]) >= 1L, 'prepare_tissue_grid should populate cell IDs in the cells list-column')
assert_true(any(grid$celltype == 'Naive_B' & grid$tissue == 'nose' & grid$n_cells == 4L), 'prepare_tissue_grid should count Naive_B nose cells correctly')

pres_df <- data.frame(
  lineage = 'bcell',
  celltype = c('Naive_B', 'Naive_B', 'Naive_B', 'Naive_B'),
  reference_network = c('nose', 'nose', 'sinus', 'sinus'),
  query_network = c('sinus', 'sinus', 'nose', 'nose'),
  module = c('blue', 'turquoise', 'blue', 'turquoise'),
  Zsummary.pres = c(6.0, 1.5, 3.0, 0.8),
  medianRank.pres = c(1.0, 2.0, 1.5, 2.5),
  n_common_genes = c(120L, 120L, 120L, 120L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pair_summary <- summarize_preservation_pairs(preservation_df = pres_df)
assert_true(is.data.frame(pair_summary), 'summarize_preservation_pairs should return a data.frame')
assert_true(nrow(pair_summary) == 2L, 'summarize_preservation_pairs should collapse module rows into pair-level summaries')
assert_true(any(pair_summary$reference_network == 'nose' & pair_summary$query_network == 'sinus' & pair_summary$modules_z_gt_2 == 1L), 'pair summary should count modules with Zsummary > 2')
assert_true(all(pair_summary$n_common_genes == 120L), 'pair summary should preserve n_common_genes per pair')

parsed <- parse_args(c('--rds-path', '/tmp/example.rds'))
assert_true(identical(parsed$metacell_min_cells, parsed$metacell_k), 'parse_args should default metacell_min_cells to metacell_k when not specified')

.captured_group_by <- NULL
.captured_metacell_min_cells <- NULL
.call_order <- character()

assign(
  'hdwgcna_global_setup',
  function(seurat_obj,
           out_dir,
           group_by,
           gene_select_mode,
           gene_fraction,
           gene_n_top,
           metacell_k,
           metacell_target,
           max_shared,
           target_metacells,
           metacell_reduction,
           metacell_dims,
           metacell_min_cells,
           wgcna_name) {
    .captured_group_by <<- group_by
    .captured_metacell_min_cells <<- metacell_min_cells
    seurat_obj
  },
  envir = .GlobalEnv
)
assign('hdwgcna_set_datexpr', function(seurat_obj, group_name, group_by, wgcna_name) seurat_obj, envir = .GlobalEnv)
assign('hdwgcna_test_soft_power', function(seurat_obj, out_dir, wgcna_name, r2_cutoff, fallback_power) list(obj = seurat_obj, recommended_power = 8L), envir = .GlobalEnv)
assign('hdwgcna_construct_network', function(seurat_obj, soft_power, net_type, tom_type, cor_type, deep_split, min_module_size, merge_cut_height, detect_cut_height, wgcna_name) seurat_obj, envir = .GlobalEnv)
assign('GetModules', function(seurat_obj, wgcna_name) data.frame(module = c('blue', 'grey'), gene_name = c('g1', 'g2'), stringsAsFactors = FALSE), envir = .GlobalEnv)
assign('hdwgcna_export_modules', function(seurat_obj, out_dir, wgcna_name) invisible(NULL), envir = .GlobalEnv)
assign('hdwgcna_module_eigengenes', function(seurat_obj, group_by, harmonize_by, wgcna_name) {
  .call_order <<- c(.call_order, 'module_eigengenes')
  seurat_obj
}, envir = .GlobalEnv)
assign('hdwgcna_hub_genes', function(seurat_obj, wgcna_name) {
  .call_order <<- c(.call_order, 'hub_genes')
  list(obj = seurat_obj, hub_genes = data.frame(module = 'blue', gene_name = 'g1', stringsAsFactors = FALSE))
}, envir = .GlobalEnv)
assign('GetDatExpr', function(seurat_obj, wgcna_name) matrix(1, nrow = 2, ncol = 2), envir = .GlobalEnv)

on.exit({
  rm('hdwgcna_global_setup', envir = .GlobalEnv)
  rm('hdwgcna_set_datexpr', envir = .GlobalEnv)
  rm('hdwgcna_test_soft_power', envir = .GlobalEnv)
  rm('hdwgcna_construct_network', envir = .GlobalEnv)
  rm('GetModules', envir = .GlobalEnv)
  rm('hdwgcna_export_modules', envir = .GlobalEnv)
  rm('hdwgcna_module_eigengenes', envir = .GlobalEnv)
  rm('hdwgcna_hub_genes', envir = .GlobalEnv)
  rm('GetDatExpr', envir = .GlobalEnv)
}, add = TRUE)

tmp_out <- tempfile('runner_out_')
dir.create(tmp_out, recursive = TRUE, showWarnings = FALSE)
run_args <- list(
  out_dir = tmp_out,
  lineage_name = 'lineage',
  sample_col = 'sample',
  tissue_col = 'tissue',
  gene_select_mode = 'fraction',
  gene_fraction = 0.05,
  gene_n_top = 3000L,
  metacell_k = 1L,
  metacell_target_use = 50000,
  max_shared = 0L,
  target_metacells = 250L,
  metacell_reduction = NA_character_,
  metacell_dims = integer(),
  metacell_min_cells = 1L,
  wgcna_name = 'hdWGCNA_tissue',
  soft_power_r2_cutoff = 0.85,
  soft_power_fallback = 12L,
  soft_power = NA_integer_,
  network_type = 'signed hybrid',
  tom_type = 'signed',
  cor_type = 'pearson',
  deep_split = 4L,
  min_module_size = 20L,
  merge_cut_height = 0.25,
  detect_cut_height = NA_real_
)

res <- run_single_tissue_network(
  seurat_obj = obj,
  celltype_label = 'Naive_B',
  tissue_label = 'nose',
  cells = grid$cells[[which(grid$celltype == 'Naive_B' & grid$tissue == 'nose')[1]]],
  args = run_args
)

assert_true(identical(.captured_group_by, c('sample', 'tissue')), 'run_single_tissue_network should preserve tissue metadata by grouping metacells with sample + tissue')
assert_true(identical(.captured_metacell_min_cells, 2L), 'run_single_tissue_network should enforce a two-metacell-capable minimum cell threshold')
assert_true(identical(.call_order, c('module_eigengenes', 'hub_genes')), 'run_single_tissue_network should compute module eigengenes before hub genes')
assert_true(identical(res$status, 'ok'), 'run_single_tissue_network stubbed smoke should succeed')

skip_args <- run_args
skip_args$metacell_min_cells <- 3L
skip_res <- run_single_tissue_network(
  seurat_obj = obj,
  celltype_label = 'Naive_B',
  tissue_label = 'nose',
  cells = grid$cells[[which(grid$celltype == 'Naive_B' & grid$tissue == 'nose')[1]]],
  args = skip_args
)

assert_true(identical(skip_res$status, 'skipped'), 'run_single_tissue_network should skip tissues with fewer than two metacell-ready sample groups')

assign('hdwgcna_set_datexpr', function(seurat_obj, group_name, group_by, wgcna_name) {
  stop('Too few genes with valid expression levels in the required number of samples.', call. = FALSE)
}, envir = .GlobalEnv)

setdatexpr_skip_res <- run_single_tissue_network(
  seurat_obj = obj,
  celltype_label = 'Naive_B',
  tissue_label = 'nose',
  cells = grid$cells[[which(grid$celltype == 'Naive_B' & grid$tissue == 'nose')[1]]],
  args = run_args
)

assert_true(identical(setdatexpr_skip_res$status, 'skipped'), 'run_single_tissue_network should downgrade sparse SetDatExpr failures to skipped')
assert_true(grepl('SetDatExpr skipped due to insufficient valid metacell gene coverage', setdatexpr_skip_res$record$reason), 'SetDatExpr sparse-data skips should record an explanatory reason')

cat('hdwgcna tissue-preservation runner regression test passed.\n')
