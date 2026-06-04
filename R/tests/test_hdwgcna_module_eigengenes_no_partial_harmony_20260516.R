#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

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

.captured_group_by_vars <- 'unset'
.captured_getmes_harmonized <- 'unset'
helper_env <- environment(hdwgcna_module_eigengenes)

assign(
  'ModuleEigengenes',
  function(seurat_obj, group.by.vars = NULL, wgcna_name = NULL, ...) {
    .captured_group_by_vars <<- group.by.vars
    seurat_obj
  },
  envir = helper_env
)
assign(
  'GetMEs',
  function(seurat_obj, harmonized = TRUE, wgcna_name = NULL) {
    .captured_getmes_harmonized <<- harmonized
    as.data.frame(
      matrix(
        0,
        nrow = ncol(seurat_obj),
        ncol = 2,
        dimnames = list(colnames(seurat_obj), c('ME1', 'ME2'))
      )
    )
  },
  envir = helper_env
)

on.exit({
  rm('ModuleEigengenes', envir = helper_env)
  rm('GetMEs', envir = helper_env)
}, add = TRUE)

hdwgcna_module_eigengenes(
  seurat_obj = seu,
  group_by = 'cell_type_final_l3',
  harmonize_by = NULL,
  wgcna_name = 'hdWGCNA'
)

assert_true(is.null(.captured_group_by_vars), 'group_by should not be partially matched to group.by.vars when harmonize_by is NULL')
assert_true(isFALSE(.captured_getmes_harmonized), 'GetMEs should use harmonized=FALSE when Harmony is not requested')

cat('hdwgcna ModuleEigengenes no-partial-Harmony regression test passed.\n')