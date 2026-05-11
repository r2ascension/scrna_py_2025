#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

make_toy_stromal_smc <- function() {
  set.seed(123)
  genes <- c(
    "RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3",
    "ACTA2", "TAGLN", "CNN1", "MYH11", "MYL9",
    "COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2",
    "CXCL12", "IL6", "CCL2", "POSTN", "LRRC15"
  )
  cluster_levels <- c("Pericyte_like", "Intermediate_SMC", "SMC_pulmonary", "SMC_systemic")
  n_per_cluster <- c(12, 12, 12, 12)
  cells <- paste0("toy_cell_", seq_len(sum(n_per_cluster)))
  cluster_assign <- rep(cluster_levels, times = n_per_cluster)

  base_lambda <- matrix(2, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
  base_lambda[match(c("RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3"), genes), cluster_assign == "Pericyte_like"] <- 8
  base_lambda[match(c("ACTA2", "TAGLN", "CNN1"), genes), cluster_assign == "Intermediate_SMC"] <- 6
  base_lambda[match(c("ACTA2", "TAGLN", "CNN1", "MYH11", "MYL9"), genes), cluster_assign == "SMC_pulmonary"] <- 10
  base_lambda[match(c("ACTA2", "TAGLN", "CNN1", "MYH11", "MYL9"), genes), cluster_assign == "SMC_systemic"] <- 9
  base_lambda[match(c("COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2"), genes), cluster_assign == "SMC_pulmonary"] <- 7
  base_lambda[match(c("POSTN", "LRRC15", "CXCL12"), genes), cluster_assign == "SMC_systemic"] <- 7

  counts <- matrix(rpois(length(base_lambda), lambda = as.vector(base_lambda)), nrow = nrow(base_lambda), dimnames = dimnames(base_lambda))
  seu <- CreateSeuratObject(counts = counts)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu$cell_type_L3 <- cluster_assign
  seu$tissue <- rep(c("pulmonary", "systemic"), each = 24)
  seu$sample <- rep(c("S1", "S2", "S3", "S4"), each = 12)
  seu$condition <- rep(c("control", "disease"), length.out = length(cells))
  seu$seurat_clusters <- rep(0:3, times = n_per_cluster)

  centers_umap <- rbind(
    c(0, 0),
    c(1, 0),
    c(2, 1),
    c(2, -1)
  )
  centers_pca <- rbind(
    c(0, 0, 0, 0, 0),
    c(2, 0, 0, 0, 0),
    c(4, 1, 0, 0, 0),
    c(4, -1, 0, 0, 0)
  )
  umap <- do.call(rbind, lapply(seq_along(cluster_levels), function(i) {
    cbind(
      rnorm(n_per_cluster[i], mean = centers_umap[i, 1], sd = 0.08),
      rnorm(n_per_cluster[i], mean = centers_umap[i, 2], sd = 0.08)
    )
  }))
  pca <- do.call(rbind, lapply(seq_along(cluster_levels), function(i) {
    mat <- matrix(rnorm(n_per_cluster[i] * 5, mean = 0, sd = 0.1), ncol = 5)
    sweep(mat, 2, centers_pca[i, ], "+")
  }))
  harmony <- pca + matrix(rnorm(length(pca), sd = 0.03), ncol = 5)

  rownames(umap) <- rownames(pca) <- rownames(harmony) <- colnames(seu)
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))
  colnames(harmony) <- paste0("harmony_", seq_len(ncol(harmony)))

  seu[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(seu))
  seu[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_", assay = DefaultAssay(seu))
  seu[["harmony"]] <- CreateDimReducObject(embeddings = harmony, key = "harmony_", assay = DefaultAssay(seu))
  seu
}

tmp_dir <- file.path(tempdir(), "stromal_smc_trajectory_interpretation_20260429")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
input_rds <- file.path(tmp_dir, "toy_stromal_smc.rds")
output_dir <- file.path(tmp_dir, "output")

obj <- make_toy_stromal_smc()
saveRDS(obj, input_rds)

script_env <- new.env(parent = globalenv())
script_env$INPUT_RDS <- input_rds
script_env$OUTPUT_DIR <- output_dir
script_env$ASSAY_NAME <- "RNA"
script_env$REDUCTION_NAME <- "umap"
script_env$CLUSTER_COL <- "cell_type_L3"
script_env$START_CLUSTER <- "Pericyte_like"
script_env$END_CLUSTERS <- c("SMC_pulmonary", "SMC_systemic")
script_env$RUN_QC <- TRUE
script_env$RUN_SENSITIVITY <- TRUE
script_env$RUN_DYNAMIC_GENES <- FALSE
script_env$RUN_PATHWAY_DYNAMICS <- FALSE
script_env$RUN_BRANCH_ASSOCIATION <- TRUE
script_env$RUN_MONOCLE3 <- TRUE
script_env$RUN_MONOCLE2 <- TRUE
script_env$RUN_PAGA <- TRUE
script_env$RUN_CROSS_METHOD_SUMMARY <- TRUE

sys.source("/home/h2048/script/R/stromal_smc_slingshot_20260428.R", envir = script_env)

expected_files <- c(
  "slingshot_pseudotime.tsv",
  "trajectory_qc_lineage.tsv",
  "trajectory_qc_branch_overlap.tsv",
  "trajectory_qc_lineage_context.tsv",
  file.path("sensitivity", "sensitivity_summary.tsv"),
  file.path("association", "sample_lineage_fraction.tsv"),
  file.path("monocle3", "monocle3_pseudotime.tsv"),
  file.path("monocle2", "monocle2_pseudotime.tsv"),
  file.path("paga", "paga_connectivities.csv"),
  "cross_method_trajectory_summary.tsv"
)

missing_files <- expected_files[!file.exists(file.path(output_dir, expected_files))]
assert_true(length(missing_files) == 0L, paste("Expected new interpretation outputs are missing:", paste(missing_files, collapse = ", ")))

cat("Stromal SMC trajectory interpretation smoke test passed.\n")
