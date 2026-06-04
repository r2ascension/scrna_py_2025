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
genes <- paste0('gene', seq_len(120))
cells <- paste0('cell_', seq_len(40))
counts <- matrix(0L, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))

high_genes <- genes[1:10]
threshold_genes <- genes[11:15]
below_threshold_genes <- genes[16:20]

counts[high_genes, 1:20] <- matrix(
  rpois(length(high_genes) * 20, lambda = 8) + 1L,
  nrow = length(high_genes)
)
counts[threshold_genes, 1:4] <- matrix(
  rpois(length(threshold_genes) * 4, lambda = 6) + 1L,
  nrow = length(threshold_genes)
)
counts[below_threshold_genes, 1:3] <- matrix(
  rpois(length(below_threshold_genes) * 3, lambda = 6) + 1L,
  nrow = length(below_threshold_genes)
)

seu <- CreateSeuratObject(counts = counts)
seu <- NormalizeData(seu, verbose = FALSE)

res <- suppressWarnings(hdwgcna_setup(
  seu,
  wgcna_name = 'hdWGCNA_fraction_test',
  gene_select_mode = 'fraction',
  gene_fraction = 0.10
))

selected <- GetWGCNAGenes(res, wgcna_name = 'hdWGCNA_fraction_test')
assert_true(length(selected) > 0L, 'hdwgcna_setup(fraction) should select at least one gene')
assert_true(all(selected %in% rownames(seu)), 'selected fraction-mode genes should exist in the Seurat object')
assert_true(length(unique(selected)) == length(selected), 'selected fraction-mode genes should be unique')
assert_true(all(high_genes %in% selected), 'highly detected genes should be retained in fraction mode')
assert_true(all(threshold_genes %in% selected), 'genes detected in exactly round(fraction * n_cells) cells should be retained')
assert_true(!any(below_threshold_genes %in% selected), 'genes detected below the fraction threshold should be excluded')

cat('hdwgcna fraction setup regression test passed.\n')