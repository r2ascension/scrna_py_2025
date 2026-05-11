#!/usr/bin/env Rscript
# =============================================================================
# MINIMAL TEST: Isolate 'match' error
# Run this to see EXACTLY where it fails
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

test_file <- "/home/h2048/data/source/1210/rds/GSM4695772_seurat.rds"

cat("=" * 60, "\n")
cat("MINIMAL TEST: Finding 'match' error\n")
cat("=" * 60, "\n")

# === STEP 1: Read ===
cat("STEP 1: Reading RDS...\n")
obj <- readRDS(test_file)
cat(sprintf("  ✓ Read OK: %d cells, %d features\n", ncol(obj), nrow(obj)))

# === STEP 2: Check metadata ===
cat("\nSTEP 2: Checking metadata types...\n")
for (col in colnames(obj@meta.data)) {
  cl <- class(obj@meta.data[[col]])[1]
  cat(sprintf("  %s: %s", col, cl))
  if (cl == "factor") {
    cat(" [FACTOR - needs conversion!]")
  }
  cat("\n")
}

# === STEP 3: Test direct assignment (OLD WAY - SHOULD FAIL) ===
cat("\nSTEP 3: Testing DIRECT assignment (old way)...\n")
tryCatch(
  {
    obj_test <- obj # Make a copy
    obj_test$sample <- obj_test$orig.ident # Direct assignment
    cat("  ✓ Direct assignment worked (unexpected!)\n")
  },
  error = function(e) {
    cat(sprintf("  ✗ Direct assignment FAILED: %s\n", e$message))
    cat("  This confirms the factor issue!\n")
  }
)

# === STEP 4: Test with conversion (NEW WAY - SHOULD WORK) ===
cat("\nSTEP 4: Testing with as.character() conversion...\n")
tryCatch(
  {
    obj$sample <- as.character(obj$orig.ident)
    cat("  ✓ Conversion worked!\n")
    cat(sprintf("  sample class: %s\n", class(obj$sample)[1]))
  },
  error = function(e) {
    cat(sprintf("  ✗ Conversion FAILED: %s\n", e$message))
    stop("Even conversion failed - different issue!")
  }
)

# === STEP 5: Test metadata operations ===
cat("\nSTEP 5: Testing metadata operations...\n")

# Test CreateSeuratObject with factor metadata
cat("  5a. CreateSeuratObject with factor metadata...\n")
tryCatch(
  {
    counts_small <- GetAssayData(obj, slot = "counts")[1:100, 1:50]
    md_factor <- obj@meta.data[colnames(counts_small), ]

    cat(sprintf(
      "     Metadata has %d factor columns\n",
      sum(sapply(md_factor, is.factor))
    ))

    tmp1 <- CreateSeuratObject(counts = counts_small)
    tmp1@meta.data <- md_factor

    cat("  ✗ FAILED at line above!\n")
  },
  error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n", e$message))
    cat(
      "  This is the problem! CreateSeuratObject doesn't like factor metadata!\n"
    )
  }
)

# Test with converted metadata
cat("  5b. CreateSeuratObject with character metadata...\n")
tryCatch(
  {
    counts_small <- GetAssayData(obj, slot = "counts")[1:100, 1:50]
    md_char <- obj@meta.data[colnames(counts_small), ]

    # Convert all factors
    for (col in colnames(md_char)) {
      if (is.factor(md_char[[col]])) {
        md_char[[col]] <- as.character(md_char[[col]])
      }
    }

    cat(sprintf(
      "     Metadata has %d factor columns (should be 0)\n",
      sum(sapply(md_char, is.factor))
    ))

    tmp2 <- CreateSeuratObject(counts = counts_small)
    tmp2@meta.data <- md_char

    cat("  ✓ SUCCESS with character conversion!\n")
  },
  error = function(e) {
    cat(sprintf("  ✗ Still failed: %s\n", e$message))
  }
)

cat("\n" + "=" * 60 + "\n")
cat("TEST COMPLETE\n")
cat("=" * 60 + "\n")
