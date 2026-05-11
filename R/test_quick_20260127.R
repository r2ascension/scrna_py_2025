#!/usr/bin/env Rscript
# ==============================================================================
# Quick Test: User's Exact Code
# ==============================================================================

library(clusterProfiler)
library(fanyi)

# Load data
cat("Loading GO BP enrichment...\n")
go_bp_enrich <- readRDS("/home/h2048/data/R/0127/bcell_interpret_v2_4_1/reports/go_bp_enrich.rds")
cat(sprintf("Loaded: %d clusters\n\n", 
           length(unique(go_bp_enrich@compareClusterResult$Cluster))))

# Setup
ctx <- paste(
  "Human normal tissue B cells (non-diseased baseline).",
  "Goal: annotate B-cell subclusters and summarize baseline functional states.",
  "Focus on: naive vs memory vs GC-related programs, class-switching, plasma cell differentiation,",
  "antigen presentation, proliferation/cell cycle, interferon-like signatures (if present), and tissue residency/trafficking."
)

Sys.setenv(DEEPSEEK_API_KEY = "sk-ed1879cf6fa14b04aac9cb6c078a3d05")
fanyi::set_translate_option(
  key = Sys.getenv("DEEPSEEK_API_KEY"),
  source = "deepseek"
)

MODEL <- "deepseek-chat"

# Test 1: List input (your code)
cat("=== Test 1: list(go_bp_enrich) - annotation ===\n")
y_anno <- tryCatch(
  {
    clusterProfiler::interpret(
      list(go_bp_enrich),
      task = "annotation",
      context = ctx,
      model = MODEL,
      n_pathways = 10
    )
  },
  error = function(e) {
    cat("[ERROR]:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(y_anno)) {
  cat(sprintf("[OK] Success! Got %d clusters\n", length(y_anno)))
  cat(sprintf("First cluster: %s\n", names(y_anno)[1]))
  cat(sprintf("Fields: %s\n\n", paste(names(y_anno[[1]]), collapse = ", ")))
  
  # Save for inspection
  saveRDS(y_anno, "test_anno_output.rds")
  cat("Saved to: test_anno_output.rds\n\n")
  
  # Check for cell_type
  if (!is.null(y_anno[[1]]$cell_type)) {
    cat("✓ Has cell_type field\n")
  } else if (!is.null(y_anno[[1]]$phenotype)) {
    cat("⚠ Has phenotype field (not cell_type)\n")
  }
}

cat("\n")

# Test 2: List input - phenotype
cat("=== Test 2: list(go_bp_enrich) - phenotype ===\n")
y_pheno <- tryCatch(
  {
    clusterProfiler::interpret(
      list(go_bp_enrich),
      task = "phenotype",
      context = ctx,
      model = MODEL,
      n_pathways = 10
    )
  },
  error = function(e) {
    cat("[ERROR]:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(y_pheno)) {
  cat(sprintf("[OK] Success! Got %d clusters\n", length(y_pheno)))
  cat(sprintf("First cluster: %s\n", names(y_pheno)[1]))
  cat(sprintf("Fields: %s\n\n", paste(names(y_pheno[[1]]), collapse = ", ")))
  
  saveRDS(y_pheno, "test_pheno_output.rds")
  cat("Saved to: test_pheno_output.rds\n\n")
  
  # Check for phenotype
  if (!is.null(y_pheno[[1]]$phenotype)) {
    cat("✓ Has phenotype field\n")
    cat(sprintf("  Example: %s\n", substr(y_pheno[[1]]$phenotype, 1, 80)))
  }
}

cat("\n")

# Test 3: Single object - interpretation
cat("=== Test 3: go_bp_enrich (no list) - interpretation ===\n")
y_mech <- tryCatch(
  {
    clusterProfiler::interpret(
      go_bp_enrich,
      task = "interpretation",
      context = ctx,
      model = MODEL,
      n_pathways = 10
    )
  },
  error = function(e) {
    cat("[ERROR]:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(y_mech)) {
  cat(sprintf("[OK] Success! Got %d clusters\n", length(y_mech)))
  cat(sprintf("First cluster: %s\n", names(y_mech)[1]))
  cat(sprintf("Fields: %s\n\n", paste(names(y_mech[[1]]), collapse = ", ")))
  
  saveRDS(y_mech, "test_mech_output.rds")
  cat("Saved to: test_mech_output.rds\n\n")
}

# Summary
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("list(go_bp_enrich) + annotation:     %s\n", 
           ifelse(!is.null(y_anno), "✓", "✗")))
cat(sprintf("list(go_bp_enrich) + phenotype:      %s\n", 
           ifelse(!is.null(y_pheno), "✓", "✗")))
cat(sprintf("go_bp_enrich (no list) + interpretation: %s\n", 
           ifelse(!is.null(y_mech), "✓", "✗")))

cat("\nOutput files created:\n")
if (!is.null(y_anno)) cat("  - test_anno_output.rds\n")
if (!is.null(y_pheno)) cat("  - test_pheno_output.rds\n")
if (!is.null(y_mech)) cat("  - test_mech_output.rds\n")

cat("\nNext: Inspect with str(readRDS('test_anno_output.rds'))\n")
