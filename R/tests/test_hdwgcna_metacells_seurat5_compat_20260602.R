#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

set.seed(20260602)
genes <- paste0('gene', seq_len(300))
cells <- paste0('cell', seq_len(120))
counts <- matrix(
  rpois(length(genes) * length(cells), lambda = 5),
  nrow = length(genes),
  dimnames = list(genes, cells)
)
counts[seq_len(40), ] <- counts[seq_len(40), ] + matrix(
  rpois(40 * length(cells), lambda = 10),
  nrow = 40
)

seu <- CreateSeuratObject(counts = counts)
seu$sample <- rep(c('S1', 'S2'), each = 60)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 10, verbose = FALSE)

setup_obj <- suppressWarnings(hdwgcna_setup(
  seu,
  wgcna_name = 'hdWGCNA_metacell_test',
  gene_select_mode = 'fraction',
  gene_fraction = 0.05
))

res <- hdwgcna_metacells(
  seurat_obj = setup_obj,
  group_by = 'sample',
  metacell_k = 10,
  metacell_target = 50000,
  max_shared = 5,
  target_metacells = 20,
  metacell_reduction = 'pca',
  metacell_dims = 1:10,
  metacell_min_cells = 20,
  wgcna_name = 'hdWGCNA_metacell_test'
)

assert_true(inherits(res, 'Seurat'), 'hdwgcna_metacells should succeed on Seurat v5 objects with the compatibility shim enabled')

cat('hdwgcna metacells Seurat v5 compatibility smoke test passed.\n')
