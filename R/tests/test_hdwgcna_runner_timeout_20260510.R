#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

make_tiny_seurat <- function() {
  set.seed(42)
  genes <- c('MS4A1','CD79A','CD79B','BANK1','HLA-DRA','CD74','MKI67','TOP2A')
  cells <- paste0('cell_', seq_len(12))
  counts <- matrix(
    rpois(length(genes) * length(cells), lambda = 5),
    nrow = length(genes),
    dimnames = list(genes, cells)
  )
  seu <- CreateSeuratObject(counts = counts)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu$cell_type_final_l3 <- rep(c('Naive_B', 'Memory_B'), each = 6)
  seu
}

seu <- make_tiny_seurat()
CELLTYPE_COL <- 'cell_type_final_l3'
SAMPLE_COL <- 'sample'

assign('hdwgcna_set_datexpr', function(seurat_obj, group_name, group_by, wgcna_name) {
  e <- structure(
    list(message = 'simulated celltype timeout', call = NULL),
    class = c('TimeoutException', 'error', 'condition')
  )
  stop(e)
}, envir = .GlobalEnv)
assign('hdwgcna_test_soft_power', function(seurat_obj, out_dir, wgcna_name) list(obj = seurat_obj, recommended_power = 6L), envir = .GlobalEnv)
assign('hdwgcna_construct_network', function(seurat_obj, soft_power, wgcna_name) seurat_obj, envir = .GlobalEnv)
assign('GetModules', function(seurat_obj, wgcna_name) data.frame(module = 'turquoise', gene_name = 'MS4A1'), envir = .GlobalEnv)
assign('ModuleEigengenes', function(...) stop('ModuleEigengenes should not be reached after timeout'), envir = .GlobalEnv)
assign('GetMEs', function(...) data.frame(ME1 = numeric(0)), envir = .GlobalEnv)

on.exit({
  rm('hdwgcna_set_datexpr', envir = .GlobalEnv)
  rm('hdwgcna_test_soft_power', envir = .GlobalEnv)
  rm('hdwgcna_construct_network', envir = .GlobalEnv)
  rm('GetModules', envir = .GlobalEnv)
  rm('ModuleEigengenes', envir = .GlobalEnv)
  rm('GetMEs', envir = .GlobalEnv)
}, add = TRUE)

out_dir <- file.path(tempdir(), 'hdwgcna_timeout_test')
res <- hdwgcna_run_celltypes(
  seurat_obj = seu,
  celltypes = 'Naive_B',
  out_dir = out_dir,
  group_by = 'cell_type_final_l3',
  wgcna_name = 'hdWGCNA',
  celltype_timeout_sec = 1
)

assert_true(is.list(res), 'hdwgcna_run_celltypes should return a list after timeout')
assert_true('Naive_B' %in% names(res), 'timeout result should contain the tested celltype')
assert_true(identical(res[['Naive_B']]$status, 'timeout'), 'timed-out celltype should be marked timeout')
assert_true(file.exists(file.path(out_dir, 'hdwgcna_Naive_B', 'hdwgcna_celltype_status.json')), 'timeout should write per-celltype status JSON')

cat('hdwgcna celltype-timeout regression test passed.\n')
