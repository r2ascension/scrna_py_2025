#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
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

assign(
  'ModuleEigengenes',
  function(seurat_obj, group.by, wgcna_name) {
    layers <- SeuratObject::Layers(seurat_obj[['RNA']])
    if (!'scale.data' %in% layers) {
      stop('Need to run ScaleData before running ModuleEigengenes with group.by.vars option.')
    }
    seurat_obj
  },
  envir = .GlobalEnv
)

assign(
  'GetMEs',
  function(seurat_obj, wgcna_name) {
    as.data.frame(
      matrix(
        0,
        nrow = ncol(seurat_obj),
        ncol = 2,
        dimnames = list(colnames(seurat_obj), c('ME1', 'ME2'))
      )
    )
  },
  envir = .GlobalEnv
)

on.exit({
  rm('ModuleEigengenes', envir = .GlobalEnv)
  rm('GetMEs', envir = .GlobalEnv)
}, add = TRUE)

res <- hdwgcna_module_eigengenes(
  seurat_obj = seu,
  group_by = 'cell_type_final_l3',
  wgcna_name = 'hdWGCNA'
)

assert_true(
  'scale.data' %in% SeuratObject::Layers(res[['RNA']]),
  'hdwgcna_module_eigengenes should create scale.data before calling ModuleEigengenes'
)

cat('hdwgcna helper scaling regression test passed.\n')
