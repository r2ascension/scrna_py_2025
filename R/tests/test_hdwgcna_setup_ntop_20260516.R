#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

set.seed(20260516)
genes <- paste0('gene_', seq_len(80))
cells <- paste0('cell_', seq_len(30))
counts <- matrix(
  rpois(length(genes) * length(cells), lambda = 5),
  nrow = length(genes),
  dimnames = list(genes, cells)
)
counts[seq_len(12), ] <- counts[seq_len(12), ] + matrix(
  rpois(12 * length(cells), lambda = 15),
  nrow = 12
)

seu <- CreateSeuratObject(counts = counts)
seu <- NormalizeData(seu, verbose = FALSE)

res <- suppressWarnings(hdwgcna_setup(
  seu,
  wgcna_name = 'hdWGCNA_ntop_test',
  gene_select_mode = 'n_top',
  gene_n_top = 25L
))

selected <- GetWGCNAGenes(res, wgcna_name = 'hdWGCNA_ntop_test')
assert_true(length(selected) == 25L, 'hdwgcna_setup(n_top) should select exactly gene_n_top features when enough genes exist')
assert_true(all(selected %in% rownames(seu)), 'selected n_top features should be present in the Seurat object')

cat('hdwgcna n_top setup regression test passed.\n')