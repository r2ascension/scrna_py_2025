#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

write_json <- function(x, path) {
  if (!requireNamespace('jsonlite', quietly = TRUE)) stop('jsonlite is required for this test', call. = FALSE)
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = 'null')
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
called <- FALSE

assign('hdwgcna_set_datexpr', function(...) {
  called <<- TRUE
  stop('resume should skip completed celltype before SetDatExpr')
}, envir = .GlobalEnv)

on.exit({
  rm('hdwgcna_set_datexpr', envir = .GlobalEnv)
}, add = TRUE)

out_dir <- file.path(tempdir(), 'hdwgcna_resume_celltypes_test')
ct_dir <- file.path(out_dir, 'hdwgcna_Naive_B')
dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
write_json(
  list(
    celltype = 'Naive_B',
    status = 'ok',
    output_dir = ct_dir,
    module_ids = c('turquoise'),
    n_modules = 1L
  ),
  file.path(ct_dir, 'hdwgcna_celltype_status.json')
)

res <- hdwgcna_run_celltypes(
  seurat_obj = seu,
  celltypes = 'Naive_B',
  out_dir = out_dir,
  group_by = 'cell_type_final_l3',
  wgcna_name = 'hdWGCNA',
  resume_celltypes = TRUE
)

assert_true(is.list(res), 'hdwgcna_run_celltypes should return a list')
assert_true('Naive_B' %in% names(res), 'resumed result should contain the skipped celltype')
assert_true(identical(res[['Naive_B']]$status, 'ok'), 'existing ok status should be preserved')
assert_true(identical(res[['Naive_B']]$skipped_existing, TRUE), 'resumed celltype should be marked skipped_existing')
assert_true(!called, 'completed celltype should not call hdwgcna_set_datexpr when resume_celltypes=TRUE')

cat('hdwgcna resume-celltypes regression test passed.\n')
