#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_no_error <- function(expr, msg) {
  tryCatch(
    force(expr),
    error = function(e) stop(sprintf('%s: %s', msg, conditionMessage(e)), call. = FALSE)
  )
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
  CreateSeuratObject(counts = counts)
}

seu <- make_tiny_seurat()
shim_active <- hdwgcna_install_seurat_v5_slot_compat(verbose = FALSE)

if (packageVersion('SeuratObject') >= '5.0.0') {
  assert_true(isTRUE(shim_active), 'slot/layer compatibility shim should activate for SeuratObject >= 5')

  counts_from_seurat <- assert_no_error(
    Seurat::GetAssayData(seu, slot = 'counts'),
    'Seurat::GetAssayData(slot=) should be translated to layer='
  )
  counts_from_seuratobject <- assert_no_error(
    SeuratObject::GetAssayData(seu, slot = 'counts'),
    'SeuratObject::GetAssayData(slot=) should be translated to layer='
  )

  assert_true(
    identical(dim(counts_from_seurat), dim(counts_from_seuratobject)),
    'Seurat and SeuratObject slot-compatible calls should return matrices with matching dimensions'
  )

  seu2 <- assert_no_error(
    Seurat::SetAssayData(seu, slot = 'data', new.data = counts_from_seurat),
    'Seurat::SetAssayData(slot=) should be translated to layer='
  )
  assert_true(
    'data' %in% SeuratObject::Layers(seu2[['RNA']]),
    'SetAssayData(slot="data") compatibility call should create/update the data layer'
  )
}

cat('hdWGCNA Seurat v5 slot/layer compatibility regression test passed.\n')