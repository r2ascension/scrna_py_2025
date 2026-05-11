#!/usr/bin/env Rscript

bundle_path <- "/home/h2048/script/R/program_architecture_bundle_20260428_v1.R"
source(bundle_path)

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_identical <- function(x, y, msg) {
  if (!identical(x, y)) stop(msg, call. = FALSE)
}

required_fns <- c(
  "pa_run_roe_runner",
  "pa_run_ggtree_runner",
  "pa_run_rogue_runner"
)

for (fn in required_fns) {
  assert_true(exists(fn, mode = "function"), sprintf("%s should be exported by the bundle", fn))
}

make_tiny_seurat <- function() {
  suppressPackageStartupMessages({
    library(Seurat)
  })
  set.seed(42)
  genes <- c("MS4A1","CD79A","CD79B","BANK1","HLA-DRA","CD74","MKI67","TOP2A", paste0("GENE_", seq_len(192)))
  cells <- paste0("cell_", seq_len(24))
  counts <- matrix(rpois(length(genes) * length(cells), lambda = 3), nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
  counts[1:4, 1:12] <- counts[1:4, 1:12] + 6
  counts[5:8, 13:24] <- counts[5:8, 13:24] + 6
  counts[9:40, seq(1, 24, by = 2)] <- counts[9:40, seq(1, 24, by = 2)] + 2
  seu <- CreateSeuratObject(counts = counts)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu$tissue <- rep(rep(c("nose", "lung parenchyma"), each = 6), times = 2)
  seu$sample <- rep(paste0("S", seq_len(6)), each = 4)
  seu$cell_type_L3 <- rep(c("Naive_B", "Memory_B"), each = 12)
  seu
}

tmp_dir <- file.path(tempdir(), "pa_support_runners")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

seu <- make_tiny_seurat()

roe_res <- pa_run_roe_runner(
  input_obj = seu,
  output_dir = file.path(tmp_dir, "roe"),
  row_key = "cell_type_L3",
  col_key = "tissue",
  min_row_sum = 1
)
assert_true(is.list(roe_res$roe_packet), "Ro/e runner should return roe_packet")
assert_true(file.exists(file.path(tmp_dir, "roe", "roe_indices.csv")), "Ro/e runner should write roe_indices.csv")

ggtree_res <- pa_run_ggtree_runner(
  input_matrix = roe_res$roe_packet$roe_matrix,
  output_dir = file.path(tmp_dir, "ggtree"),
  tree_prefix = "toy_roe",
  margins = c("rows", "cols")
)
assert_true(is.list(ggtree_res$tree_summary), "ggtree runner should return tree_summary")
assert_true(file.exists(file.path(tmp_dir, "ggtree", "toy_roe_rows_tree.pdf")), "ggtree runner should write row tree pdf")

rogue_res <- pa_run_rogue_runner(
  seurat_obj = seu,
  output_dir = file.path(tmp_dir, "rogue"),
  group_col = "cell_type_L3",
  split_col = "tissue",
  min_cells_per_group = 4,
  max_cells_per_group = 100,
  rogue_min_cells = 4,
  rogue_min_genes = 20
)
assert_true(is.data.frame(rogue_res$rogue_table), "ROGUE runner should return rogue_table")
assert_true(nrow(rogue_res$rogue_table) >= 2, "ROGUE runner should emit at least two rows")
assert_true(file.exists(file.path(tmp_dir, "rogue", "rogue_long.tsv")), "ROGUE runner should write rogue_long.tsv")

validation_packet <- pa_build_program_validation_packet(
  roe_table = roe_res$roe_packet$roe_long_table,
  roe_matrix = roe_res$roe_packet$roe_matrix,
  rogue_table = rogue_res$rogue_table,
  rogue_summary = rogue_res$rogue_summary,
  tree_summary = ggtree_res$tree_summary
)
summary_lines <- pa_validation_packet_summary_lines(validation_packet)
assert_true(any(grepl("Ro/e", summary_lines, fixed = TRUE)), "Validation summary should mention Ro/e support")
assert_true(any(grepl("ROGUE", summary_lines, fixed = TRUE)), "Validation summary should mention ROGUE support")
assert_true(any(grepl("ggtree", summary_lines, fixed = TRUE)), "Validation summary should mention ggtree support")

cat("All support runner tests passed.\n")
