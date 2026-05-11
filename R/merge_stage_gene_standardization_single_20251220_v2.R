#!/usr/bin/env Rscript
################################################################################
# Ultra-Simple Batch Fix: ENSG Unmapped Issues (All-in-One)
# 超简洁批量修复：ENSG未映射问题（单文件版）
#
# Author: r2end
# Date: 2024-12-20
################################################################################

library(Seurat)
library(Matrix)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)

#===============================================================================
# ⭐ CONFIGURATION - ONLY CHANGE THESE
# ⭐ 配置 - 只需改这里
#===============================================================================

CLEANED_DIR <- "/home/h2048/data/R/1215/merge/cleaned_samples"
OUTPUT_DIR  <- "/home/h2048/data/R/1215/merge/cleaned_samples_FIXED"

# Optional: specify samples (NULL = all)
SAMPLES_TO_PROCESS <- NULL

#===============================================================================
# FUNCTIONS (DO NOT MODIFY)
#===============================================================================

clean_ensembl_versions <- function(gene_names) {
  sub("\\.\\d+$", "", gene_names)
}

get_genes_from_assay <- function(seurat_obj, assay = "RNA") {
  tryCatch(
    Features(seurat_obj, assay = assay),
    error = function(e) rownames(seurat_obj[[assay]])
  )
}

build_global_gene_mapping_fixed <- function(
  all_genes_unique,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Building Global Gene Mapping (ENSG-FIRST)               ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }

  base_keys <- if (clean_ensembl) clean_ensembl_versions(all_genes_unique) else all_genes_unique
  is_ensg   <- grepl("^ENSG[0-9]+$", base_keys)

  gene_mapping <- data.frame(
    original_name   = all_genes_unique,
    official_symbol = NA_character_,
    mapping_source  = NA_character_,
    stringsAsFactors = FALSE
  )

  # Layer 0: ENSEMBL FIRST
  if (verbose) cat("Layer 0: ENSEMBL-first for ENSG...\n")
  idx0 <- which(is_ensg)
  if (length(idx0) > 0) {
    m0 <- tryCatch(
      mapIds(gene_db, keys = base_keys[idx0], column = "SYMBOL", 
             keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx0))
    )
    ok0 <- !is.na(m0) & nzchar(m0)
    gene_mapping$official_symbol[idx0[ok0]] <- as.character(m0[ok0])
    gene_mapping$mapping_source[idx0[ok0]]  <- "ENSEMBL"
    if (verbose) cat(sprintf("  ✓ %d ENSG mapped\n", sum(ok0)))
  }

  # Layer 1: SYMBOL
  if (verbose) cat("Layer 1: SYMBOL...\n")
  idx1 <- which(is.na(gene_mapping$official_symbol) & !is_ensg)
  if (length(idx1) > 0) {
    m1 <- tryCatch(
      mapIds(gene_db, keys = all_genes_unique[idx1], column = "SYMBOL",
             keytype = "SYMBOL", multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx1))
    )
    ok1 <- !is.na(m1) & nzchar(m1)
    gene_mapping$official_symbol[idx1[ok1]] <- as.character(m1[ok1])
    gene_mapping$mapping_source[idx1[ok1]]  <- "SYMBOL"
    if (verbose) cat(sprintf("  ✓ %d mapped\n", sum(ok1)))
  }

  # Layer 2: ALIAS
  if (verbose) cat("Layer 2: ALIAS...\n")
  idx2 <- which(is.na(gene_mapping$official_symbol))
  if (length(idx2) > 0) {
    m2 <- tryCatch(
      mapIds(gene_db, keys = all_genes_unique[idx2], column = "SYMBOL",
             keytype = "ALIAS", multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx2))
    )
    ok2 <- !is.na(m2) & nzchar(m2)
    gene_mapping$official_symbol[idx2[ok2]] <- as.character(m2[ok2])
    gene_mapping$mapping_source[idx2[ok2]]  <- "ALIAS"
    if (verbose) cat(sprintf("  ✓ %d mapped\n", sum(ok2)))
  }

  # Layer 3: ENSEMBL fallback
  if (verbose) cat("Layer 3: ENSEMBL fallback...\n")
  idx3 <- which(is.na(gene_mapping$official_symbol))
  if (length(idx3) > 0) {
    m3 <- tryCatch(
      mapIds(gene_db, keys = base_keys[idx3], column = "SYMBOL",
             keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx3))
    )
    ok3 <- !is.na(m3) & nzchar(m3)
    gene_mapping$official_symbol[idx3[ok3]] <- as.character(m3[ok3])
    gene_mapping$mapping_source[idx3[ok3]]  <- "ENSEMBL"
    if (verbose) cat(sprintf("  ✓ %d mapped\n", sum(ok3)))
  }

  # Layer 4: UNMAPPED
  idx4 <- which(is.na(gene_mapping$official_symbol) | !nzchar(gene_mapping$official_symbol))
  gene_mapping$official_symbol[idx4] <- gene_mapping$original_name[idx4]
  gene_mapping$mapping_source[idx4]  <- "UNMAPPED"

  if (verbose) {
    cat("\n📊 Summary:\n")
    print(table(gene_mapping$mapping_source))
    cat("\n")
  }

  gene_mapping
}

apply_gene_mapping_fixed <- function(
  seurat_obj,
  global_mapping,
  assay = "RNA",
  aggregate_method = "sum",
  verbose = TRUE
) {
  current_genes <- get_genes_from_assay(seurat_obj, assay = assay)

  counts_mat <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = "counts")
  )
  if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- as(counts_mat, "dgCMatrix")

  msub <- global_mapping[match(current_genes, global_mapping$original_name), ]
  if (any(is.na(msub$original_name))) stop("Some genes missing from global mapping!")

  official <- msub$official_symbol
  
  # ENSEMBL patch
  base_keys <- clean_ensembl_versions(current_genes)
  is_ensg   <- grepl("^ENSG[0-9]+$", base_keys)
  needs_patch <- is_ensg & (official == current_genes)
  
  if (any(needs_patch)) {
    if (verbose) cat("  🔧 ENSEMBL patch...\n")
    patched <- mapIds(org.Hs.eg.db, keys = base_keys[needs_patch],
                      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    patched <- as.character(patched)
    ok <- !is.na(patched) & nzchar(patched)
    idx <- which(needs_patch)
    official[idx[ok]] <- patched[ok]
    if (verbose) cat(sprintf("    ✓ %d genes patched\n", sum(ok)))
  }

  # Aggregate
  f <- factor(official, levels = unique(official))
  A <- sparseMatrix(i = as.integer(f), j = seq_along(f), x = 1,
                    dims = c(nlevels(f), length(f)))
  agg <- A %*% counts_mat
  rownames(agg) <- levels(f)
  colnames(agg) <- colnames(counts_mat)

  if (aggregate_method == "mean") {
    gs <- as.numeric(table(f))
    agg <- Diagonal(x = 1/gs) %*% agg
  }

  new_assay <- tryCatch(
    CreateAssay5Object(counts = agg),
    error = function(e) CreateAssayObject(counts = agg)
  )
  seurat_obj[[assay]] <- new_assay

  if (verbose) {
    cat(sprintf("  ✓ Final: %d genes\n", nrow(agg)))
  }

  seurat_obj
}

#===============================================================================
# MAIN WORKFLOW
#===============================================================================

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  Batch Fix ENSG - Ultra Simple (All-in-One)             ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Scan directory
cat("📁 Scanning: ", CLEANED_DIR, "\n")
if (!dir.exists(CLEANED_DIR)) stop("Directory not found!")

all_files <- list.files(CLEANED_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(all_files) == 0) stop("No RDS files found!")

sample_names <- tools::file_path_sans_ext(basename(all_files))
sample_paths <- setNames(all_files, sample_names)

if (!is.null(SAMPLES_TO_PROCESS)) {
  sample_paths <- sample_paths[names(sample_paths) %in% SAMPLES_TO_PROCESS]
}

cat(sprintf("✓ Found %d samples\n\n", length(sample_paths)))

# Quick diagnostic
cat("🔍 Checking ENSG genes...\n")
for (nm in names(sample_paths)) {
  obj <- readRDS(sample_paths[nm])
  n_ensg <- sum(grepl("^ENSG", rownames(obj)))
  cat(sprintf("  %s: %d ENSG\n", nm, n_ensg))
  rm(obj); gc(verbose = FALSE)
}
cat("\n")

# Build mapping
cat("Building global mapping...\n")
all_genes <- unique(unlist(lapply(sample_paths, function(p) {
  obj <- readRDS(p); genes <- rownames(obj); rm(obj); gc(verbose = FALSE); genes
})))

mapping <- build_global_gene_mapping_fixed(all_genes, verbose = TRUE)

# Create output
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(mapping, file.path(OUTPUT_DIR, "global_gene_mapping_FIXED.csv"))

# Process samples
cat("\n🔧 Processing samples...\n\n")
report <- list()

for (nm in names(sample_paths)) {
  cat("─────────────────────────────────────────────────────────\n")
  cat(sprintf("Processing: %s\n", nm))
  cat("─────────────────────────────────────────────────────────\n")
  
  obj <- readRDS(sample_paths[nm])
  genes_before <- nrow(obj)
  ensg_before <- sum(grepl("^ENSG", rownames(obj)))
  
  obj <- apply_gene_mapping_fixed(obj, mapping, verbose = TRUE)
  
  genes_after <- nrow(obj)
  ensg_after <- sum(grepl("^ENSG", rownames(obj)))
  
  output_file <- file.path(OUTPUT_DIR, paste0(nm, "_FIXED.rds"))
  saveRDS(obj, output_file)
  
  report[[nm]] <- data.frame(
    sample = nm,
    genes_before = genes_before,
    genes_after = genes_after,
    ensg_before = ensg_before,
    ensg_after = ensg_after,
    ensg_fixed = ensg_before - ensg_after
  )
  
  cat(sprintf("✓ Saved: %s\n\n", basename(output_file)))
  rm(obj); gc(verbose = FALSE)
}

# Summary
report_df <- do.call(rbind, report)
fwrite(report_df, file.path(OUTPUT_DIR, "batch_fix_report.csv"))

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  ✓ Complete!                                             ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

print(report_df, row.names = FALSE)
cat("\n💾 Output: ", OUTPUT_DIR, "\n\n")