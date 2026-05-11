#!/usr/bin/env Rscript
################################################################################
# Ultra-Simple Quick Start: Auto-Scan All RDS Files
# 超简洁快速开始：自动扫描所有RDS文件
#
# 只需设置一个路径，自动处理目录中的所有RDS文件
# Just set one path, automatically process all RDS files in directory
#
# Author: r2end
# Date: 2024-12-18
################################################################################
library(Seurat)
library(org.Hs.eg.db)
#===============================================================================
# ⭐ ONLY SET THIS: Your cleaned RDS directory
# ⭐ 只需设置这个：你的cleaned RDS文件目录
#===============================================================================

CLEANED_DIR <- "/home/h2048/data/R/1215/merge/cleaned_samples"  # 👈 改这里！
OUTPUT_DIR <- paste0(CLEANED_DIR, "_standardized_output")  # 自动生成输出目录

#===============================================================================
# Auto-scan and process
# 自动扫描和处理
#===============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  Auto Gene Standardization and Merge\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Step 1: Scan directory
cat("Step 1: Scanning directory...\n")
cat(sprintf("  Directory: %s\n", CLEANED_DIR))

if (!dir.exists(CLEANED_DIR)) {
  stop(sprintf("Directory not found: %s", CLEANED_DIR))
}

# Find all .rds files
sample_files <- list.files(
  path = CLEANED_DIR,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = FALSE  # Don't search subdirectories
)

if (length(sample_files) == 0) {
  stop(sprintf("No .rds files found in: %s", CLEANED_DIR))
}

cat(sprintf("  Found %d RDS files\n\n", length(sample_files)))

# Step 2: Generate sample names (remove .rds extension)
sample_names <- tools::file_path_sans_ext(basename(sample_files))
sample_paths <- setNames(sample_files, sample_names)

cat("Step 2: Sample list:\n")
for (i in seq_along(sample_paths)) {
  cat(sprintf("  [%d] %s\n", i, names(sample_paths)[i]))
}
cat("\n")

# Step 3: Load module
cat("Step 3: Loading gene standardization module...\n")
source("/home/h2048/script/R/merge_stage_gene_standardization_20251218.R")
cat("  ✓ Module loaded\n\n")

# Step 4: Run workflow
cat("Step 4: Running standardization and merge...\n")
cat("  (This may take several minutes)\n\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

result <- standardize_and_merge_workflow(
  sample_paths = sample_paths,
  output_dir = OUTPUT_DIR,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  aggregate_method = "sum",
  alignment_strategy = "union",
  merge_data = FALSE,
  verbose = TRUE
)

# Step 5: Summary
cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  ✓ Complete!\n")
cat("═══════════════════════════════════════════════════════════\n\n")

if (!result$validation$all_passed) {
  stop("Validation failed! Please check output above.")
}

merged_obj <- result$merged_obj

cat(sprintf("Final object: %d genes × %d cells (%d samples)\n\n",
            nrow(merged_obj), ncol(merged_obj), length(sample_paths)))

cat("Output files saved in:\n")
cat(sprintf("  %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat(sprintf("  - merged_seurat_standardized.rds (main output)\n"))
cat(sprintf("  - global_gene_mapping/ (mapping details)\n"))
cat(sprintf("  - standardized_samples/ (individual samples)\n\n"))

# Show any issues
unmapped_rate <- sum(result$global_mapping$mapping_source == "UNMAPPED") / 
                 nrow(result$global_mapping) * 100

if (unmapped_rate > 5) {
  cat(sprintf("⚠️  Note: %.1f%% genes are UNMAPPED\n", unmapped_rate))
  cat("   (May indicate non-standard gene names)\n\n")
}

if (!is.null(result$conflicts)) {
  cat(sprintf("ℹ️  %d synonym conflicts were resolved by aggregation\n\n", 
              nrow(result$conflicts)))
}

cat("Next: Load merged object with:\n")
cat(sprintf("  merged <- readRDS('%s/merged_seurat_standardized.rds')\n\n",
            OUTPUT_DIR))
