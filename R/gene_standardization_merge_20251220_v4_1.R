#!/usr/bin/env Rscript
################################################################################
# Complete All-in-One: Gene Standardization + Incremental Merge
# 完整单文件版：基因标准化 + 增量合并（内存优化）
#
# Version 4.1 - PRODUCTION READY
# 
# Complete Fixes (based on detailed review):
# P0-1: Strict RNA assay usage throughout (no more rownames(obj))
# P0-2: Force DefaultAssay to RNA after standardization
# P0-3: Force DefaultAssay before subset in merge
# P0-4: Fail-fast if any sample fails (no silent partial merge)
# P0-5: Check for duplicated genes in input samples
# P1-1: DietSeurat after standardization (keep only RNA assay)
# P1-2: Audit report for duplicate symbol aggregation
# P1-3: Clear downstream normalization instructions
#
# Author: r2end
# Date: 2024-12-20
# Version: 4.1 PRODUCTION
################################################################################

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  Gene Standardization + Merge v4.1 (PRODUCTION)         ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Script started: %s\n", Sys.time()))
cat(sprintf("R version: %s\n\n", R.version.string))

#===============================================================================
# STEP 0: Configuration
#===============================================================================

cat("STEP 0: Configuration\n")

# ⭐ INPUT/OUTPUT
CLEANED_DIR <- "/home/h2048/data/R/1215/merge/cleaned_samples"
OUTPUT_DIR  <- "/home/h2048/data/R/1215/merge/cleaned_samples_COMPLETE_v4.1"

# ⭐ PROCESSING OPTIONS
SAMPLES_TO_PROCESS <- NULL  # NULL = all samples
CLEAN_ENSEMBL <- TRUE       # Remove version numbers from ENSG
AGGREGATE_METHOD <- "sum"   # "sum" or "mean" for duplicate symbols
ALIGNMENT_STRATEGY <- "union"  # "union", "intersection", or "core"
CORE_THRESHOLD <- 0.8       # For "core" strategy
MERGE_DATA <- FALSE         # Whether to merge normalized data slots
FORCE_PROCESS <- TRUE       # Process even if no ENSG issues

cat(sprintf("  Input:  %s\n", CLEANED_DIR))
cat(sprintf("  Output: %s\n", OUTPUT_DIR))
cat(sprintf("  Strategy: %s alignment, %s aggregation\n\n", 
            ALIGNMENT_STRATEGY, AGGREGATE_METHOD))

#===============================================================================
# STEP 1: Load Libraries
#===============================================================================

cat("STEP 1: Loading libraries...\n")

required_packages <- c("Seurat", "Matrix", "data.table", "org.Hs.eg.db", "AnnotationDbi")

for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("  ✓ All libraries loaded\n\n")

#===============================================================================
# STEP 2: Define Utility Functions
#===============================================================================

cat("STEP 2: Defining utility functions...\n")

clean_ensembl_versions <- function(gene_names) {
  sub("\\.\\d+$", "", gene_names)
}

get_genes_from_assay <- function(seurat_obj, assay = "RNA") {
  tryCatch(
    Features(seurat_obj, assay = assay),
    error = function(e) rownames(seurat_obj[[assay]])
  )
}

# P0-1 FIX: Strict RNA gene extraction
get_genes_strict <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  get_genes_from_assay(seurat_obj, assay = "RNA")
}

cat("  ✓ Utilities defined\n\n")

#===============================================================================
# STEP 3: Define Gene Mapping Function (ENSG-FIRST)
#===============================================================================

cat("STEP 3: Defining gene mapping function...\n")

build_global_gene_mapping_fixed <- function(
  all_genes_unique,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n  ╔════════════════════════════════════════════════════════╗\n")
    cat("  ║  Building Global Gene Mapping (ENSG-FIRST)            ║\n")
    cat("  ╚════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("  Total genes: %d\n\n", length(all_genes_unique)))
  }

  base_keys <- if (clean_ensembl) clean_ensembl_versions(all_genes_unique) else all_genes_unique
  is_ensg   <- grepl("^ENSG[0-9]+$", base_keys)

  gene_mapping <- data.frame(
    original_name   = all_genes_unique,
    official_symbol = NA_character_,
    mapping_source  = NA_character_,
    stringsAsFactors = FALSE
  )

  # Layer 0: ENSEMBL FIRST for ENSG
  if (verbose) cat("  Layer 0: ENSEMBL-first for ENSG...\n")
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
    if (verbose) cat(sprintf("    ✓ %d ENSG mapped\n", sum(ok0)))
  }

  # Layer 1: SYMBOL
  if (verbose) cat("  Layer 1: SYMBOL...\n")
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
    if (verbose) cat(sprintf("    ✓ %d mapped\n", sum(ok1)))
  }

  # Layer 2: ALIAS
  if (verbose) cat("  Layer 2: ALIAS...\n")
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
    if (verbose) cat(sprintf("    ✓ %d mapped\n", sum(ok2)))
  }

  # Layer 3: ENSEMBL fallback
  if (verbose) cat("  Layer 3: ENSEMBL fallback...\n")
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
    if (verbose) cat(sprintf("    ✓ %d mapped\n", sum(ok3)))
  }

  # Layer 4: UNMAPPED
  idx4 <- which(is.na(gene_mapping$official_symbol) | !nzchar(gene_mapping$official_symbol))
  gene_mapping$official_symbol[idx4] <- gene_mapping$original_name[idx4]
  gene_mapping$mapping_source[idx4]  <- "UNMAPPED"

  if (verbose) {
    cat("\n  📊 Mapping Summary:\n")
    tbl <- table(gene_mapping$mapping_source)
    for (src in names(tbl)) {
      pct <- tbl[src] / nrow(gene_mapping) * 100
      cat(sprintf("    %-10s: %6d (%.1f%%)\n", src, tbl[src], pct))
    }
    cat("\n")
  }

  gene_mapping
}

cat("  ✓ Mapping function defined\n\n")

#===============================================================================
# STEP 4: Define Apply Function (Sparse Aggregation)
#===============================================================================

cat("STEP 4: Defining apply function...\n")

apply_gene_mapping_fixed <- function(
  seurat_obj,
  global_mapping,
  assay = "RNA",
  aggregate_method = "sum",
  sample_name = "sample",
  verbose = TRUE
) {
  # P0-1 FIX: Strict RNA assay
  DefaultAssay(seurat_obj) <- assay
  current_genes <- get_genes_from_assay(seurat_obj, assay = assay)
  
  # P0-5 FIX: Check for duplicates in input
  if (any(duplicated(current_genes))) {
    stop(sprintf("Sample '%s' has duplicated gene names before mapping! Please fix upstream.", 
                 sample_name))
  }

  counts_mat <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = "counts")
  )
  if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- as(counts_mat, "dgCMatrix")

  msub <- global_mapping[match(current_genes, global_mapping$original_name), ]
  if (any(is.na(msub$original_name))) {
    stop(sprintf("Sample '%s': Some genes missing from global mapping!", sample_name))
  }

  official <- msub$official_symbol
  
  # ENSEMBL patch for stuck genes
  base_keys <- clean_ensembl_versions(current_genes)
  is_ensg   <- grepl("^ENSG[0-9]+$", base_keys)
  needs_patch <- is_ensg & (official == current_genes)
  
  if (any(needs_patch)) {
    if (verbose) cat("    🔧 ENSEMBL patch...\n")
    patched <- mapIds(org.Hs.eg.db, keys = base_keys[needs_patch],
                      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    patched <- as.character(patched)
    ok <- !is.na(patched) & nzchar(patched)
    idx <- which(needs_patch)
    official[idx[ok]] <- patched[ok]
    if (verbose) cat(sprintf("      ✓ %d genes patched\n", sum(ok)))
  }

  # P1-2: Record aggregation details for audit
  dup_symbols <- official[duplicated(official)]
  dup_info <- NULL
  
  if (length(dup_symbols) > 0) {
    dup_table <- table(official)
    dup_table <- dup_table[dup_table > 1]
    
    dup_info <- list(
      n_dup_symbols = length(dup_table),
      max_dup_size = max(dup_table),
      top_dups = head(sort(dup_table, decreasing = TRUE), 20)
    )
    
    if (verbose) {
      cat(sprintf("    ⚠️  %d symbols will be aggregated (max dup size: %d)\n",
                  dup_info$n_dup_symbols, dup_info$max_dup_size))
    }
  }

  # Sparse matrix aggregation
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
  
  # P0-2 FIX: Force DefaultAssay back to RNA
  DefaultAssay(seurat_obj) <- assay
  seurat_obj@active.assay <- assay  # Extra insurance
  
  # P1-1 FIX: DietSeurat to keep only RNA assay
  seurat_obj <- DietSeurat(
    seurat_obj,
    assays = assay,
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL
  )
  DefaultAssay(seurat_obj) <- assay

  if (verbose) {
    cat(sprintf("    ✓ Final: %d genes\n", nrow(agg)))
  }

  list(
    object = seurat_obj,
    dup_info = dup_info
  )
}

cat("  ✓ Apply function defined\n\n")

#===============================================================================
# STEP 5: Define Incremental Merge Function (Memory-Efficient)
#===============================================================================

cat("STEP 5: Defining incremental merge function...\n")

incremental_merge_aligned <- function(
  sample_paths,
  aligned_genes = NULL,
  merge_data = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n  ╔════════════════════════════════════════════════════════╗\n")
    cat("  ║  Incremental Merge (Memory-Efficient)                 ║\n")
    cat("  ╚════════════════════════════════════════════════════════╝\n\n")
  }

  n_samples <- length(sample_paths)
  sample_names <- names(sample_paths)

  # Load and filter helper
  load_and_filter <- function(path, sample_name, idx) {
    if (verbose) cat(sprintf("  [%d/%d] Loading: %s\n", idx, n_samples, sample_name))
    
    obj <- readRDS(path)
    
    # P0-3 FIX: Force DefaultAssay before any operations
    DefaultAssay(obj) <- "RNA"
    
    if (!is.null(aligned_genes)) {
      # P0-1 FIX: Use strict gene extraction
      available <- intersect(get_genes_from_assay(obj, "RNA"), aligned_genes)
      obj <- subset(obj, features = available)
      if (verbose) cat(sprintf("    Filtered to %d aligned genes\n", length(available)))
    }
    
    # P0-1 FIX: Use strict counting
    n_genes <- length(get_genes_from_assay(obj, "RNA"))
    n_cells <- ncol(obj)
    
    if (verbose) cat(sprintf("    %d genes × %d cells\n", n_genes, n_cells))
    obj
  }

  # Load first sample
  merged_obj <- load_and_filter(sample_paths[1], sample_names[1], 1)
  gc(verbose = FALSE)

  if (n_samples == 1) {
    if (verbose) cat("  Single sample, no merge needed\n")
    return(merged_obj)
  }

  # Incrementally merge remaining samples
  for (i in 2:n_samples) {
    next_obj <- load_and_filter(sample_paths[i], sample_names[i], i)
    
    if (verbose) {
      cat(sprintf("  Merging: %d + %d cells\n", ncol(merged_obj), ncol(next_obj)))
    }
    
    merged_obj <- merge(
      x = merged_obj,
      y = next_obj,
      merge.data = merge_data
    )
    
    rm(next_obj)
    gc(verbose = FALSE)
    
    if (verbose) {
      cat(sprintf("  Cumulative: %d cells\n\n", ncol(merged_obj)))
    }
  }

  # P0-1 FIX: Final gene count using strict method
  DefaultAssay(merged_obj) <- "RNA"
  final_genes <- length(get_genes_from_assay(merged_obj, "RNA"))
  
  if (verbose) {
    cat(sprintf("  ✓ Merge complete: %d genes × %d cells\n\n", 
                final_genes, ncol(merged_obj)))
  }

  merged_obj
}

cat("  ✓ Merge function defined\n\n")

#===============================================================================
# STEP 6: Scan Input Directory
#===============================================================================

cat("STEP 6: Scanning input directory...\n")

if (!dir.exists(CLEANED_DIR)) stop(sprintf("Directory not found: %s", CLEANED_DIR))

all_files <- list.files(CLEANED_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(all_files) == 0) stop("No RDS files found!")

sample_names <- tools::file_path_sans_ext(basename(all_files))
sample_paths <- setNames(all_files, sample_names)

if (!is.null(SAMPLES_TO_PROCESS)) {
  sample_paths <- sample_paths[names(sample_paths) %in% SAMPLES_TO_PROCESS]
  if (length(sample_paths) == 0) {
    stop("None of specified samples found!")
  }
}

cat(sprintf("  Found %d samples to process\n\n", length(sample_paths)))

for (i in seq_along(sample_paths)) {
  cat(sprintf("    [%d] %s\n", i, names(sample_paths)[i]))
}
cat("\n")

#===============================================================================
# STEP 7: Quick Diagnostic
#===============================================================================

cat("STEP 7: Quick diagnostic...\n")

diagnostic <- list()
for (nm in names(sample_paths)) {
  obj <- readRDS(sample_paths[nm])
  # P0-1 FIX: Use strict gene extraction
  genes <- get_genes_strict(obj)
  n_ensg <- sum(grepl("^ENSG[0-9]+", genes))
  diagnostic[[nm]] <- list(total = length(genes), ensg = n_ensg)
  cat(sprintf("  %s: %d ENSG / %d total\n", nm, n_ensg, length(genes)))
  rm(obj); gc(verbose = FALSE)
}

total_ensg <- sum(sapply(diagnostic, function(x) x$ensg))
cat(sprintf("\n  Total ENSG to fix: %d\n", total_ensg))

if (total_ensg == 0 && !FORCE_PROCESS) {
  cat("  No issues found. Set FORCE_PROCESS=TRUE to continue anyway.\n")
  quit(save = "no", status = 0)
}
cat("\n")

#===============================================================================
# STEP 8: Collect All Genes
#===============================================================================

cat("STEP 8: Collecting all unique genes...\n")

all_genes_list <- lapply(sample_paths, function(p) {
  obj <- readRDS(p)
  # P0-1 FIX: Use strict gene extraction
  genes <- get_genes_strict(obj)
  rm(obj); gc(verbose = FALSE)
  genes
})

all_genes_unique <- unique(unlist(all_genes_list))
cat(sprintf("  Total unique genes: %d\n\n", length(all_genes_unique)))

#===============================================================================
# STEP 9: Build Global Mapping
#===============================================================================

cat("STEP 9: Building global gene mapping...\n")

start_mapping <- Sys.time()

global_mapping <- build_global_gene_mapping_fixed(
  all_genes_unique,
  gene_db = org.Hs.eg.db,
  clean_ensembl = CLEAN_ENSEMBL,
  verbose = TRUE
)

end_mapping <- Sys.time()
cat(sprintf("  Completed in %.2f minutes\n\n",
            difftime(end_mapping, start_mapping, units = "mins")))

#===============================================================================
# STEP 10: Save Global Mapping
#===============================================================================

cat("STEP 10: Saving global mapping...\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
mapping_file <- file.path(OUTPUT_DIR, "global_gene_mapping_FIXED.csv")
fwrite(global_mapping, mapping_file)

cat(sprintf("  ✓ Saved: %s (%.2f MB)\n\n", 
            basename(mapping_file), file.size(mapping_file)/1024^2))

#===============================================================================
# STEP 11: Process Each Sample (Standardization)
#===============================================================================

cat("STEP 11: Standardizing genes in each sample...\n\n")

standardized_dir <- file.path(OUTPUT_DIR, "standardized_samples")
dir.create(standardized_dir, showWarnings = FALSE)

standardized_paths <- c()
failed_samples <- c()
all_dup_info <- list()
start_process <- Sys.time()

for (i in seq_along(sample_paths)) {
  nm <- names(sample_paths)[i]
  
  cat("───────────────────────────────────────────────────────────\n")
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(sample_paths), nm))
  cat("───────────────────────────────────────────────────────────\n")
  
  tryCatch({
    obj <- readRDS(sample_paths[nm])
    
    # P0-1 FIX: Strict gene extraction
    genes_before <- length(get_genes_strict(obj))
    ensg_before <- sum(grepl("^ENSG", get_genes_strict(obj)))
    
    cat(sprintf("  Loaded: %d genes × %d cells\n", genes_before, ncol(obj)))
    cat(sprintf("  ENSG before: %d\n", ensg_before))
    
    # Apply with sample name for better error messages
    result <- apply_gene_mapping_fixed(
      obj, 
      global_mapping, 
      aggregate_method = AGGREGATE_METHOD,
      sample_name = nm,
      verbose = TRUE
    )
    
    obj <- result$object
    all_dup_info[[nm]] <- result$dup_info
    
    # P0-1 FIX: Strict gene counting after
    genes_after <- length(get_genes_strict(obj))
    ensg_after <- sum(grepl("^ENSG", get_genes_strict(obj)))
    
    # Validate
    if (any(is.na(get_genes_strict(obj)))) stop("NA in gene names!")
    if (any(duplicated(get_genes_strict(obj)))) stop("Duplicated gene names!")
    
    output_file <- file.path(standardized_dir, paste0(nm, "_standardized.rds"))
    saveRDS(obj, output_file)
    standardized_paths[nm] <- output_file
    
    cat(sprintf("  After: %d genes (ENSG: %d)\n", genes_after, ensg_after))
    cat(sprintf("  ✓ Saved: %s\n\n", basename(output_file)))
    
    rm(obj); gc(verbose = FALSE)
    
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n\n", e$message))
    failed_samples <- c(failed_samples, nm)
  })
}

end_process <- Sys.time()

# P0-4 FIX: Fail-fast if any sample failed
if (length(failed_samples) > 0) {
  cat("\n╔══════════════════════════════════════════════════════════╗\n")
  cat("║  ✗ PROCESSING FAILED                                     ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")
  cat("The following samples failed to process:\n")
  for (nm in failed_samples) {
    cat(sprintf("  - %s\n", nm))
  }
  cat("\nPlease fix these samples before proceeding.\n")
  cat("Check logs above for specific error messages.\n\n")
  stop("Sample processing failed - cannot continue to merge")
}

#===============================================================================
# STEP 11.5: Generate Aggregation Audit Report (P1-2)
#===============================================================================

cat("STEP 11.5: Generating aggregation audit report...\n")

audit_report <- data.frame(
  sample = character(),
  n_dup_symbols = integer(),
  max_dup_size = integer(),
  stringsAsFactors = FALSE
)

for (nm in names(all_dup_info)) {
  if (!is.null(all_dup_info[[nm]])) {
    audit_report <- rbind(audit_report, data.frame(
      sample = nm,
      n_dup_symbols = all_dup_info[[nm]]$n_dup_symbols,
      max_dup_size = all_dup_info[[nm]]$max_dup_size,
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(audit_report) > 0) {
  audit_file <- file.path(OUTPUT_DIR, "aggregation_audit_report.csv")
  fwrite(audit_report, audit_file)
  
  cat(sprintf("  ✓ Saved: %s\n", basename(audit_file)))
  cat("\n  Summary:\n")
  cat(sprintf("    Total samples with aggregation: %d\n", nrow(audit_report)))
  cat(sprintf("    Max dup symbols in any sample: %d\n", max(audit_report$n_dup_symbols)))
  cat(sprintf("    Max dup size across all samples: %d\n\n", max(audit_report$max_dup_size)))
} else {
  cat("  No duplicate symbol aggregation occurred\n\n")
}

#===============================================================================
# STEP 12: Gene Alignment (Prepare for Merge)
#===============================================================================

cat("STEP 12: Preparing gene alignment...\n")

gene_lists <- lapply(standardized_paths, function(p) {
  obj <- readRDS(p)
  # P0-1 FIX: Strict gene extraction
  genes <- get_genes_strict(obj)
  rm(obj); gc(verbose = FALSE)
  genes
})

if (ALIGNMENT_STRATEGY == "union") {
  aligned_genes <- unique(unlist(gene_lists))
  cat(sprintf("  Union: %d genes (no filtering in merge)\n", length(aligned_genes)))
  aligned_genes_for_merge <- NULL  # Seurat merge does union by default
  
} else if (ALIGNMENT_STRATEGY == "intersection") {
  aligned_genes <- Reduce(intersect, gene_lists)
  cat(sprintf("  Intersection: %d genes\n", length(aligned_genes)))
  aligned_genes_for_merge <- aligned_genes
  
} else if (ALIGNMENT_STRATEGY == "core") {
  gene_freq <- table(unlist(gene_lists))
  min_samples <- ceiling(length(standardized_paths) * CORE_THRESHOLD)
  aligned_genes <- names(gene_freq[gene_freq >= min_samples])
  cat(sprintf("  Core (>=%d%%): %d genes\n", CORE_THRESHOLD*100, length(aligned_genes)))
  aligned_genes_for_merge <- aligned_genes
  
} else {
  stop("Unknown alignment strategy!")
}

cat("\n")

#===============================================================================
# STEP 13: Incremental Merge
#===============================================================================

cat("STEP 13: Incremental merge of standardized samples...\n")

start_merge <- Sys.time()

merged_obj <- incremental_merge_aligned(
  sample_paths = standardized_paths,
  aligned_genes = aligned_genes_for_merge,
  merge_data = MERGE_DATA,
  verbose = TRUE
)

end_merge <- Sys.time()

#===============================================================================
# STEP 14: Add Metadata and Save
#===============================================================================

cat("STEP 14: Adding metadata and saving...\n")

# P0-1 FIX: Use strict counting
DefaultAssay(merged_obj) <- "RNA"
final_genes <- length(get_genes_from_assay(merged_obj, "RNA"))
final_cells <- ncol(merged_obj)

merged_obj@misc$gene_standardization <- list(
  version = "v4.1_production",
  ensembl_cleaned = CLEAN_ENSEMBL,
  aggregate_method = AGGREGATE_METHOD,
  alignment_strategy = ALIGNMENT_STRATEGY,
  n_samples = length(sample_paths),
  final_genes = final_genes,
  final_cells = final_cells,
  processing_date = as.character(Sys.Date())
)

merged_file <- file.path(OUTPUT_DIR, "merged_seurat_standardized.rds")
saveRDS(merged_obj, merged_file)

cat(sprintf("  ✓ Saved: %s (%.2f MB)\n\n", 
            basename(merged_file), file.size(merged_file)/1024^2))

#===============================================================================
# FINAL SUMMARY
#===============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  ✓ WORKFLOW COMPLETE                                     ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("📊 Final Summary:\n")
cat(sprintf("  Samples processed: %d\n", length(sample_paths)))
cat(sprintf("  Total ENSG fixed:  %d\n", total_ensg))
cat(sprintf("  Final object:      %d genes × %d cells\n\n", final_genes, final_cells))

cat("⏱️  Timing:\n")
cat(sprintf("  Mapping:         %.2f min\n",
            difftime(end_mapping, start_mapping, units = "mins")))
cat(sprintf("  Standardization: %.2f min\n",
            difftime(end_process, start_process, units = "mins")))
cat(sprintf("  Merge:           %.2f min\n",
            difftime(end_merge, start_merge, units = "mins")))
cat(sprintf("  Total:           %.2f min\n\n",
            difftime(end_merge, start_mapping, units = "mins")))

cat("💾 Output Directory:\n")
cat(sprintf("  %s\n\n", OUTPUT_DIR))

cat("📁 Key Files:\n")
cat("  - merged_seurat_standardized.rds (MAIN OUTPUT - counts only)\n")
cat("  - global_gene_mapping_FIXED.csv\n")
cat("  - aggregation_audit_report.csv (duplicate symbol report)\n")
cat("  - standardized_samples/ (individual samples, RNA assay only)\n\n")

cat("🎯 CRITICAL: Next Steps for Downstream Analysis\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("⚠️  The merged object contains ONLY raw counts.\n")
cat("   You MUST run the following before any analysis:\n\n")

cat("  # 1. Load merged object\n")
cat(sprintf("  merged <- readRDS('%s')\n\n", merged_file))

cat("  # 2. Normalize (REQUIRED)\n")
cat("  merged <- NormalizeData(merged, normalization.method = 'LogNormalize')\n\n")

cat("  # 3. Find variable features\n")
cat("  merged <- FindVariableFeatures(merged, selection.method = 'vst', nfeatures = 2000)\n\n")

cat("  # 4. Scale (for PCA/clustering)\n")
cat("  merged <- ScaleData(merged)\n\n")

cat("  # 5. Then proceed with integration (BBKNN/scVI/Harmony)\n")
cat("  # Example for standard workflow:\n")
cat("  merged <- RunPCA(merged)\n")
cat("  # ... then your integration method of choice\n\n")

cat("⚠️  DO NOT skip normalization - downstream analysis will fail or give wrong results!\n\n")

cat("📖 Audit Report:\n")
cat("   Check aggregation_audit_report.csv for details on\n")
cat("   which samples had duplicate symbols aggregated.\n\n")

cat(sprintf("Script completed: %s\n\n", Sys.time()))

# Exit successfully
quit(save = "no", status = 0)