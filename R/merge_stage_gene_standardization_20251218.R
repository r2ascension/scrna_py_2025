################################################################################
# Gene Standardization and Alignment Module for Merge Stage - FIXED VERSION
# 合并阶段的基因标准化与对齐模块 - 修复版
#
# Key Fixes:
# - P0: Synonym conflicts now resolved per-sample with count aggregation
# - P1: Ensembl version number handling
# - P1: Strict assay consistency for gene names and counts
# - Added comprehensive validation checks
#
# Author: r2end
# Date: 2024-12-18
################################################################################

library(Seurat)
library(Matrix)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)

#===============================================================================
# Part 0: Utility Functions
#===============================================================================

#' Clean Ensembl IDs by removing version numbers
#' @param gene_names Character vector of gene names
#' @return Character vector with version numbers removed
clean_ensembl_versions <- function(gene_names) {
  # Remove version numbers from Ensembl IDs (e.g., ENSG00000123456.12 -> ENSG00000123456)
  cleaned <- sub("\\.\\d+$", "", gene_names)
  return(cleaned)
}

#' Get gene names from specific assay (strict)
#' @param seurat_obj Seurat object
#' @param assay Assay name
#' @return Character vector of gene names
get_genes_from_assay <- function(seurat_obj, assay = "RNA") {
  # Strict: get features from the specified assay, not default rownames
  genes <- tryCatch(
    Features(seurat_obj, assay = assay),
    error = function(e) {
      # Fallback for older Seurat versions
      rownames(seurat_obj[[assay]])
    }
  )
  return(genes)
}

#===============================================================================
# Part 1: Pre-merge Gene Analysis
#===============================================================================

#' Analyze gene names across all samples before merge
analyze_gene_landscape <- function(sample_paths, verbose = TRUE) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Pre-Merge Gene Landscape Analysis                      ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }
  
  n_samples <- length(sample_paths)
  sample_names <- names(sample_paths)
  
  if (verbose) cat("Step 1: Collecting gene lists from all samples...\n")
  
  gene_info <- list()
  for (i in seq_along(sample_paths)) {
    if (verbose) cat(sprintf("  [%d/%d] Loading: %s\n", i, n_samples, sample_names[i]))
    
    obj <- readRDS(sample_paths[i])
    
    # Strict: get genes from RNA assay
    genes <- get_genes_from_assay(obj, assay = "RNA")
    
    gene_info[[sample_names[i]]] <- list(
      genes = genes,
      n_genes = length(genes),
      n_cells = ncol(obj)
    )
  }
  
  if (verbose) cat("\nStep 2: Building global gene universe...\n")
  
  all_genes_list <- lapply(gene_info, function(x) x$genes)
  all_genes_unique <- unique(unlist(all_genes_list))
  
  gene_frequency <- table(unlist(all_genes_list))
  
  if (verbose) {
    cat(sprintf("  Total unique genes across all samples: %d\n", length(all_genes_unique)))
    cat(sprintf("  Genes present in all %d samples: %d\n", 
                n_samples, sum(gene_frequency == n_samples)))
    cat(sprintf("  Genes present in only 1 sample: %d\n", 
                sum(gene_frequency == 1)))
  }
  
  if (verbose) cat("\nStep 3: Analyzing pairwise gene overlap...\n")
  
  overlap_rates <- c()
  for (i in 1:(n_samples-1)) {
    for (j in (i+1):n_samples) {
      genes_i <- gene_info[[i]]$genes
      genes_j <- gene_info[[j]]$genes
      
      overlap <- length(intersect(genes_i, genes_j))
      union_size <- length(union(genes_i, genes_j))
      overlap_rate <- overlap / union_size
      
      overlap_rates <- c(overlap_rates, overlap_rate)
    }
  }
  
  if (verbose) {
    cat(sprintf("  Mean pairwise overlap: %.2f%%\n", mean(overlap_rates) * 100))
    cat(sprintf("  Min pairwise overlap: %.2f%%\n", min(overlap_rates) * 100))
    cat(sprintf("  Max pairwise overlap: %.2f%%\n", max(overlap_rates) * 100))
  }
  
  return(list(
    gene_info = gene_info,
    all_genes_unique = all_genes_unique,
    gene_frequency = gene_frequency,
    overlap_rates = overlap_rates,
    n_samples = n_samples,
    sample_names = sample_names
  ))
}

#===============================================================================
# Part 2: Global Gene Mapping (NO deletion of synonyms)
#===============================================================================

#' Build global gene mapping for all samples
#' @param all_genes_unique Vector of all unique genes across samples
#' @param gene_db Gene annotation database
#' @param clean_ensembl Whether to clean Ensembl version numbers
#' @param verbose Print progress
#' @return Data frame with complete gene mapping (ALL original_names preserved)
build_global_gene_mapping <- function(
  all_genes_unique,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Building Global Gene Mapping                            ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("Total genes to map: %d\n", length(all_genes_unique)))
    if (clean_ensembl) {
      cat("Ensembl version cleaning: ENABLED\n")
    }
    cat("\n")
  }
  
  gene_mapping <- data.frame(
    original_name = all_genes_unique,
    official_symbol = NA_character_,
    mapping_source = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Prepare cleaned versions for Ensembl mapping
  genes_for_ensembl <- if (clean_ensembl) {
    clean_ensembl_versions(all_genes_unique)
  } else {
    all_genes_unique
  }
  
  # Layer 1: SYMBOL
  if (verbose) cat("Layer 1: Mapping through SYMBOL...\n")
  symbol_match <- tryCatch(
    {
      mapIds(
        gene_db,
        keys = all_genes_unique,
        column = "SYMBOL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    },
    error = function(e) rep(NA, length(all_genes_unique))
  )
  
  matched <- !is.na(symbol_match)
  gene_mapping$official_symbol[matched] <- symbol_match[matched]
  gene_mapping$mapping_source[matched] <- "SYMBOL"
  
  if (verbose) cat(sprintf("  Mapped: %d genes\n", sum(matched)))
  
  # Layer 2: ALIAS
  unmatched <- is.na(gene_mapping$official_symbol)
  if (sum(unmatched) > 0) {
    if (verbose) cat("Layer 2: Mapping through ALIAS...\n")
    
    alias_match <- tryCatch(
      {
        mapIds(
          gene_db,
          keys = all_genes_unique[unmatched],
          column = "SYMBOL",
          keytype = "ALIAS",
          multiVals = "first"
        )
      },
      error = function(e) rep(NA, sum(unmatched))
    )
    
    alias_matched <- !is.na(alias_match)
    gene_mapping$official_symbol[unmatched][alias_matched] <- alias_match[alias_matched]
    gene_mapping$mapping_source[unmatched][alias_matched] <- "ALIAS"
    
    if (verbose) cat(sprintf("  Mapped: %d additional genes\n", sum(alias_matched)))
  }
  
  # Layer 3: ENSEMBL (with cleaned version)
  unmatched <- is.na(gene_mapping$official_symbol)
  if (sum(unmatched) > 0) {
    if (verbose) cat("Layer 3: Mapping through ENSEMBL...\n")
    
    ensembl_match <- tryCatch(
      {
        mapIds(
          gene_db,
          keys = genes_for_ensembl[unmatched],  # Use cleaned versions
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      },
      error = function(e) rep(NA, sum(unmatched))
    )
    
    ensembl_matched <- !is.na(ensembl_match)
    gene_mapping$official_symbol[unmatched][ensembl_matched] <- ensembl_match[ensembl_matched]
    gene_mapping$mapping_source[unmatched][ensembl_matched] <- "ENSEMBL"
    
    if (verbose) cat(sprintf("  Mapped: %d additional genes\n", sum(ensembl_matched)))
  }
  
  # Layer 4: Keep original (UNMAPPED)
  unmatched <- is.na(gene_mapping$official_symbol)
  gene_mapping$official_symbol[unmatched] <- all_genes_unique[unmatched]
  gene_mapping$mapping_source[unmatched] <- "UNMAPPED"
  
  if (verbose) {
    cat(sprintf("\nUnmapped genes: %d (%.2f%%)\n", 
                sum(unmatched), sum(unmatched)/length(all_genes_unique)*100))
  }
  
  # Summary
  if (verbose) {
    cat("\nMapping Summary:\n")
    print(table(gene_mapping$mapping_source))
    cat("\n")
  }
  
  return(gene_mapping)
}

#' Identify synonym conflicts (report only, NO deletion)
#' 
#' @param gene_mapping Global gene mapping data frame
#' @param verbose Print information
#' @return List with mapping (unchanged) and conflict report
analyze_synonym_conflicts <- function(gene_mapping, verbose = TRUE) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Analyzing Synonym Conflicts (No Deletion)               ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }
  
  # Find duplicate official symbols
  symbol_counts <- table(gene_mapping$official_symbol)
  duplicate_symbols <- names(symbol_counts)[symbol_counts > 1]
  
  if (length(duplicate_symbols) == 0) {
    if (verbose) cat("No synonym conflicts detected.\n\n")
    return(list(
      mapping = gene_mapping,  # UNCHANGED
      conflicts = NULL
    ))
  }
  
  if (verbose) {
    cat(sprintf("Found %d official symbols with multiple original names\n", 
                length(duplicate_symbols)))
    cat("These will be aggregated (counts summed) when applied to samples\n\n")
  }
  
  # Generate conflict report
  conflict_report <- list()
  
  for (symbol in duplicate_symbols) {
    rows <- gene_mapping[gene_mapping$official_symbol == symbol, ]
    original_names <- rows$original_name
    sources <- rows$mapping_source
    
    # Priority ranking for reporting
    priority_order <- c("SYMBOL", "ALIAS", "ENSEMBL", "UNMAPPED")
    best_source_idx <- which.min(match(sources, priority_order))
    
    conflict_report[[symbol]] <- data.frame(
      official_symbol = symbol,
      n_synonyms = length(original_names),
      all_original_names = paste(original_names, collapse = " | "),
      all_sources = paste(sources, collapse = " | "),
      recommended = original_names[best_source_idx],
      recommended_source = sources[best_source_idx],
      stringsAsFactors = FALSE
    )
    
    if (verbose) {
      cat(sprintf("  %s: %d synonyms\n", symbol, length(original_names)))
      cat(sprintf("    Original names: %s\n", paste(original_names, collapse = ", ")))
      cat(sprintf("    Sources: %s\n", paste(sources, collapse = ", ")))
      cat(sprintf("    → Will aggregate all when applying to samples\n"))
    }
  }
  
  conflict_df <- do.call(rbind, conflict_report)
  
  if (verbose) {
    cat(sprintf("\n✓ Identified %d synonym groups\n", nrow(conflict_df)))
    cat("✓ All mappings preserved for per-sample aggregation\n\n")
  }
  
  return(list(
    mapping = gene_mapping,  # UNCHANGED - all mappings preserved
    conflicts = conflict_df
  ))
}

#===============================================================================
# Part 3: Apply Mapping to Individual Samples (with aggregation)
#===============================================================================

#' Apply global gene mapping to a single sample with duplicate aggregation
#' 
#' @param seurat_obj Seurat object
#' @param global_mapping Global gene mapping data frame (complete, no deletions)
#' @param assay Assay name
#' @param aggregate_method Method for aggregating duplicates ("sum" or "mean")
#' @param verbose Print detailed info
#' @return Modified Seurat object with standardized gene names
apply_gene_mapping_to_sample <- function(
  seurat_obj,
  global_mapping,
  assay = "RNA",
  aggregate_method = "sum",
  verbose = TRUE
) {
  # Step 1: Get current gene names from the specified assay (strict)
  current_genes <- get_genes_from_assay(seurat_obj, assay = assay)
  
  # Step 2: Extract counts from the same assay
  counts_mat <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay, slot = "counts")
    }
  )
  
  # Ensure it's sparse
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
  }
  
  # Step 3: Match current genes to global mapping
  mapping_subset <- global_mapping[global_mapping$original_name %in% current_genes, ]
  
  if (nrow(mapping_subset) == 0) {
    stop("No genes matched between sample and global mapping!")
  }
  
  # Ensure order matches current_genes
  mapping_subset <- mapping_subset[match(current_genes, mapping_subset$original_name), ]
  
  # Check for NA matches
  if (any(is.na(mapping_subset$original_name))) {
    stop("Some genes in sample are missing from global mapping!")
  }
  
  # Step 4: Identify genes that will map to the same official symbol
  official_symbols <- mapping_subset$official_symbol
  
  duplicate_symbols <- unique(official_symbols[duplicated(official_symbols)])
  
  if (length(duplicate_symbols) > 0 && verbose) {
    cat(sprintf("  Found %d official symbols with multiple original names in this sample\n", 
                length(duplicate_symbols)))
    cat(sprintf("  These will be aggregated using method: %s\n", aggregate_method))
  }
  
  # Step 5: Aggregate counts for duplicate official symbols
  if (length(duplicate_symbols) > 0) {
    # Create aggregated matrix
    unique_symbols <- unique(official_symbols)
    aggregated_counts <- Matrix(0, 
                                nrow = length(unique_symbols), 
                                ncol = ncol(counts_mat),
                                sparse = TRUE)
    rownames(aggregated_counts) <- unique_symbols
    colnames(aggregated_counts) <- colnames(counts_mat)
    
    for (symbol in unique_symbols) {
      indices <- which(official_symbols == symbol)
      
      if (length(indices) == 1) {
        # No aggregation needed
        aggregated_counts[symbol, ] <- counts_mat[indices, ]
      } else {
        # Aggregate multiple rows
        if (aggregate_method == "sum") {
          aggregated_counts[symbol, ] <- Matrix::colSums(counts_mat[indices, , drop = FALSE])
        } else if (aggregate_method == "mean") {
          aggregated_counts[symbol, ] <- Matrix::colMeans(counts_mat[indices, , drop = FALSE])
        }
        
        if (verbose) {
          original_names <- mapping_subset$original_name[indices]
          cat(sprintf("    %s: aggregated %d rows (%s)\n", 
                      symbol, length(indices), paste(original_names, collapse = ", ")))
        }
      }
    }
    
    counts_mat <- aggregated_counts
  } else {
    # No duplicates, just rename
    rownames(counts_mat) <- official_symbols
  }
  
  # Step 6: Create new assay with standardized counts
  new_assay <- tryCatch(
    CreateAssay5Object(counts = counts_mat),
    error = function(e) {
      CreateAssayObject(counts = counts_mat)
    }
  )
  
  seurat_obj[[assay]] <- new_assay
  
  if (verbose) {
    cat(sprintf("  Final genes after standardization: %d\n", nrow(counts_mat)))
  }
  
  return(seurat_obj)
}

#' Standardize genes for all samples using global mapping
standardize_all_samples <- function(
  sample_paths,
  global_mapping,
  output_dir,
  aggregate_method = "sum",
  verbose = TRUE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Applying Global Gene Standardization to All Samples    ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }
  
  n_samples <- length(sample_paths)
  sample_names <- names(sample_paths)
  
  # Create output directory
  standardized_dir <- file.path(output_dir, "standardized_samples")
  dir.create(standardized_dir, recursive = TRUE, showWarnings = FALSE)
  
  standardized_paths <- c()
  
  for (i in seq_along(sample_paths)) {
    if (verbose) {
      cat(sprintf("[%d/%d] Processing: %s\n", i, n_samples, sample_names[i]))
    }
    
    # Load sample
    obj <- readRDS(sample_paths[i])
    genes_before <- length(get_genes_from_assay(obj, assay = "RNA"))
    
    # Apply standardization with aggregation
    obj <- apply_gene_mapping_to_sample(
      obj, 
      global_mapping, 
      assay = "RNA",
      aggregate_method = aggregate_method,
      verbose = verbose
    )
    
    genes_after <- length(get_genes_from_assay(obj, assay = "RNA"))
    
    # Validation checks (gene names)
    genes_standardized <- get_genes_from_assay(obj, assay = "RNA")
    if (any(is.na(genes_standardized))) {
      stop(sprintf("Sample %s has NA gene names after standardization!", sample_names[i]))
    }
    
    if (any(duplicated(genes_standardized))) {
      dup_genes <- unique(genes_standardized[duplicated(genes_standardized)])
      stop(sprintf(
        "Sample %s has duplicated gene names after standardization! Examples: %s",
        sample_names[i],
        paste(head(dup_genes), collapse = ", ")
      ))
    }
    
    # Save standardized sample
    output_path <- file.path(
      standardized_dir,
      paste0(sample_names[i], "_standardized.rds")
    )
    saveRDS(obj, output_path)
    
    standardized_paths[sample_names[i]] <- output_path
    
    if (verbose) {
      cat(sprintf("  Genes: %d → %d\n", genes_before, genes_after))
      cat(sprintf("  ✓ Validation passed\n"))
      cat(sprintf("  Saved: %s\n\n", basename(output_path)))
    }
    
    # Clean up
    rm(obj)
    gc(verbose = FALSE)
  }
  
  if (verbose) {
    cat("✓ All samples standardized and validated\n\n")
  }
  
  return(standardized_paths)
}

#===============================================================================
# Part 4: Gene Alignment and Smart Merge (unchanged from original)
#===============================================================================

prepare_gene_alignment <- function(
  sample_paths,
  alignment_strategy = "union",
  core_threshold = 0.8,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Gene Alignment Preparation                              ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("Strategy: %s\n", alignment_strategy))
  }
  
  # Collect gene lists
  gene_lists <- lapply(sample_paths, function(path) {
    obj <- readRDS(path)
    get_genes_from_assay(obj, assay = "RNA")  # Strict
  })
  
  # Calculate alignment
  if (alignment_strategy == "union") {
    aligned_genes <- unique(unlist(gene_lists))
    
    if (verbose) {
      cat(sprintf("  Union strategy: %d total genes\n", length(aligned_genes)))
    }
    
  } else if (alignment_strategy == "intersection") {
    aligned_genes <- Reduce(intersect, gene_lists)
    
    if (verbose) {
      cat(sprintf("  Intersection strategy: %d shared genes\n", length(aligned_genes)))
    }
    
  } else if (alignment_strategy == "core") {
    gene_frequency <- table(unlist(gene_lists))
    min_samples <- ceiling(length(sample_paths) * core_threshold)
    aligned_genes <- names(gene_frequency[gene_frequency >= min_samples])
    
    if (verbose) {
      cat(sprintf("  Core strategy (>= %.0f%% samples): %d genes\n", 
                  core_threshold * 100, length(aligned_genes)))
    }
    
  } else {
    stop("Unknown alignment strategy: ", alignment_strategy)
  }
  
  # Gene availability per sample
  availability <- data.frame(
    sample = names(sample_paths),
    total_genes = sapply(gene_lists, length),
    aligned_genes = sapply(gene_lists, function(genes) {
      sum(genes %in% aligned_genes)
    }),
    stringsAsFactors = FALSE
  )
  availability$coverage <- availability$aligned_genes / availability$total_genes
  
  if (verbose) {
    cat("\nGene coverage per sample:\n")
    print(availability)
    cat("\n")
  }
  
  return(list(
    aligned_genes = aligned_genes,
    strategy = alignment_strategy,
    availability = availability,
    gene_lists = gene_lists
  ))
}

incremental_merge_aligned <- function(
  sample_paths,
  aligned_genes = NULL,
  verbose = TRUE,
  merge_data = FALSE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Incremental Merge with Gene Alignment                  ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }
  
  n_samples <- length(sample_paths)
  sample_names <- names(sample_paths)
  
  # Helper function
  load_and_filter <- function(path, sample_name, idx, total) {
    if (verbose) {
      cat(sprintf("[%d/%d] Loading: %s\n", idx, total, sample_name))
    }
    
    obj <- readRDS(path)
    
    # Filter to aligned genes if specified
    if (!is.null(aligned_genes)) {
      available_genes <- intersect(get_genes_from_assay(obj, "RNA"), aligned_genes)
      obj <- subset(obj, features = available_genes)
      
      if (verbose) {
        cat(sprintf("  Filtered to %d aligned genes\n", length(available_genes)))
      }
    }
    
    if (verbose) {
      cat(sprintf("  Final: %d genes × %d cells\n", nrow(obj), ncol(obj)))
    }
    
    return(obj)
  }
  
  # Load first sample
  merged_obj <- load_and_filter(sample_paths[1], sample_names[1], 1, n_samples)
  gc(verbose = FALSE)
  
  if (n_samples == 1) {
    if (verbose) cat("\nOnly 1 sample, no merge needed\n")
    return(merged_obj)
  }
  
  # Incrementally merge
  for (i in 2:n_samples) {
    next_obj <- load_and_filter(sample_paths[i], sample_names[i], i, n_samples)
    
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
  
  if (verbose) {
    cat("✓ Merge complete\n")
    cat(sprintf("  Final: %d genes × %d cells\n\n", nrow(merged_obj), ncol(merged_obj)))
  }
  
  return(merged_obj)
}

#===============================================================================
# Part 5: Validation Functions
#===============================================================================

#' Comprehensive validation of standardized samples
#' 
#' @param standardized_paths Named vector of standardized sample paths
#' @param global_mapping Global gene mapping
#' @param conflicts Conflict report
#' @param verbose Print detailed checks
#' @return List with validation results
validate_standardization_quality <- function(
  standardized_paths,
  global_mapping,
  conflicts = NULL,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Validation: Gene Standardization Quality               ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
  }
  
  validation_results <- list()
  all_passed <- TRUE
  
  # Check 1: NA in rownames
  if (verbose) cat("Check 1: Verifying no NA in gene names...\n")
  
  for (sample_name in names(standardized_paths)) {
    obj <- readRDS(standardized_paths[sample_name])
    genes <- get_genes_from_assay(obj, assay = "RNA")
    
    if (any(is.na(genes))) {
      cat(sprintf("  ✗ FAIL: Sample %s has %d NA gene names\n", 
                  sample_name, sum(is.na(genes))))
      all_passed <- FALSE
    }
  }
  
  if (verbose && all_passed) {
    cat("  ✓ PASS: No NA gene names detected\n\n")
  }
  
  # Check 2: Duplicate rownames
  if (verbose) cat("Check 2: Verifying no duplicate gene names...\n")
  
  for (sample_name in names(standardized_paths)) {
    obj <- readRDS(standardized_paths[sample_name])
    genes <- get_genes_from_assay(obj, assay = "RNA")
    
    if (any(duplicated(genes))) {
      dup_genes <- genes[duplicated(genes)]
      cat(sprintf("  ✗ FAIL: Sample %s has %d duplicated genes\n", 
                  sample_name, length(unique(dup_genes))))
      cat(sprintf("    Examples: %s\n", paste(head(unique(dup_genes)), collapse = ", ")))
      all_passed <- FALSE
    }
  }
  
  if (verbose && all_passed) {
    cat("  ✓ PASS: No duplicate gene names detected\n\n")
  }
  
  # Check 3: Critical marker genes preserved
  if (verbose) cat("Check 3: Verifying critical marker genes...\n")
  
  critical_markers <- c(
    # Immune
    "CD3D", "CD3E", "CD4", "CD8A", "CD8B",
    # Cytokines - check both names
    "IL8", "CXCL8",  # Should be unified to one
    "TNF", "TNFA", "TNFAIP3",  # TNF family
    # Epithelial
    "EPCAM", "KRT18", "KRT19", "TP63",
    # Common markers
    "GAPDH", "ACTB"
  )
  
  marker_presence <- matrix(FALSE, 
                           nrow = length(standardized_paths),
                           ncol = length(critical_markers),
                           dimnames = list(names(standardized_paths), critical_markers))
  
  for (sample_name in names(standardized_paths)) {
    obj <- readRDS(standardized_paths[sample_name])
    genes <- get_genes_from_assay(obj, assay = "RNA")
    marker_presence[sample_name, ] <- critical_markers %in% genes
  }
  
  # Report markers present in all samples
  markers_in_all <- colSums(marker_presence) == nrow(marker_presence)
  
  if (verbose) {
    cat(sprintf("  Critical markers present in all samples: %d / %d\n",
                sum(markers_in_all), length(critical_markers)))
    
    if (sum(!markers_in_all) > 0) {
      cat("\n  Missing or inconsistent markers:\n")
      for (marker in critical_markers[!markers_in_all]) {
        n_present <- sum(marker_presence[, marker])
        cat(sprintf("    %s: present in %d/%d samples\n", 
                    marker, n_present, length(standardized_paths)))
      }
    }
  }
  
  # Check for synonym unification (e.g., IL8 vs CXCL8)
  if (!is.null(conflicts)) {
    if (verbose) cat("\nCheck 4: Verifying synonym unification...\n")
    
    # Check a few key synonyms
    key_synonyms <- list(
      c("IL8", "CXCL8"),
      c("TNF", "TNFA"),
      c("SEPT1", "SEP1", "SEPTIN1")
    )
    
    for (syn_group in key_synonyms) {
      for (sample_name in names(standardized_paths)) {
        obj <- readRDS(standardized_paths[sample_name])
        genes <- get_genes_from_assay(obj, assay = "RNA")
        
        present <- syn_group[syn_group %in% genes]
        
        if (length(present) > 1) {
          cat(sprintf("  ⚠ WARNING: Sample %s has multiple synonyms: %s\n",
                      sample_name, paste(present, collapse = ", ")))
          cat(sprintf("    This suggests aggregation may have failed\n"))
        }
      }
    }
    
    if (verbose) {
      cat("  ✓ Synonym unification check complete\n\n")
    }
  }
  
  validation_results$all_passed <- all_passed
  validation_results$marker_presence <- marker_presence
  
  if (verbose) {
    if (all_passed) {
      cat("╔══════════════════════════════════════════════════════════╗\n")
      cat("║  ✓ ALL VALIDATION CHECKS PASSED                         ║\n")
      cat("╚══════════════════════════════════════════════════════════╝\n\n")
    } else {
      cat("╔══════════════════════════════════════════════════════════╗\n")
      cat("║  ✗ SOME VALIDATION CHECKS FAILED - REVIEW ABOVE         ║\n")
      cat("╚══════════════════════════════════════════════════════════╝\n\n")
    }
  }
  
  return(validation_results)
}

#===============================================================================
# Part 6: Master Workflow Function
#===============================================================================

#' Complete gene standardization and merge workflow (FIXED VERSION)
standardize_and_merge_workflow <- function(
  sample_paths,
  output_dir,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  aggregate_method = "sum",
  alignment_strategy = "union",
  core_threshold = 0.8,
  merge_data = FALSE,
  verbose = TRUE
) {
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Gene Standardization and Merge Workflow (FIXED)         ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n")
    cat(sprintf("\nSamples: %d\n", length(sample_paths)))
    cat(sprintf("Output: %s\n", output_dir))
    cat(sprintf("Ensembl cleaning: %s\n", clean_ensembl))
    cat(sprintf("Aggregation method: %s\n", aggregate_method))
    cat(sprintf("Alignment strategy: %s\n", alignment_strategy))
  }
  
  # Step 1: Analyze gene landscape
  landscape <- analyze_gene_landscape(sample_paths, verbose = verbose)
  
  # Step 2: Build global gene mapping (with Ensembl cleaning)
  global_mapping <- build_global_gene_mapping(
    landscape$all_genes_unique,
    gene_db = gene_db,
    clean_ensembl = clean_ensembl,
    verbose = verbose
  )
  
  # Step 3: Analyze synonyms (NO deletion)
  analysis <- analyze_synonym_conflicts(global_mapping, verbose = verbose)
  
  # Save mapping
  mapping_dir <- file.path(output_dir, "global_gene_mapping")
  dir.create(mapping_dir, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(analysis$mapping, file.path(mapping_dir, "global_gene_mapping_complete.csv"))
  
  if (!is.null(analysis$conflicts)) {
    fwrite(analysis$conflicts, file.path(mapping_dir, "synonym_conflicts_report.csv"))
  }
  
  if (verbose) {
    cat("✓ Saved global gene mapping (all original names preserved)\n\n")
  }
  
  # Step 4: Apply standardization with aggregation
  standardized_paths <- standardize_all_samples(
    sample_paths = sample_paths,
    global_mapping = analysis$mapping,
    output_dir = output_dir,
    aggregate_method = aggregate_method,
    verbose = verbose
  )
  
  # Step 5: Validation
  validation <- validate_standardization_quality(
    standardized_paths = standardized_paths,
    global_mapping = analysis$mapping,
    conflicts = analysis$conflicts,
    verbose = verbose
  )
  
  if (!validation$all_passed) {
    stop("Validation failed! Please review errors above.")
  }
  
  # Step 6: Prepare gene alignment
  alignment <- prepare_gene_alignment(
    sample_paths = standardized_paths,
    alignment_strategy = alignment_strategy,
    core_threshold = core_threshold,
    verbose = verbose
  )
  
  # Save alignment info
  fwrite(alignment$availability, 
         file.path(output_dir, "gene_alignment_availability.csv"))
  writeLines(alignment$aligned_genes,
             file.path(output_dir, "aligned_genes.txt"))
  
  # Step 7: Incremental merge
  merged_obj <- incremental_merge_aligned(
    sample_paths = standardized_paths,
    aligned_genes = if (alignment_strategy != "union") alignment$aligned_genes else NULL,
    verbose = verbose,
    merge_data = merge_data
  )
  
  # Add metadata
  merged_obj@misc$gene_standardization <- list(
    method = "merge_stage_fixed_v1",
    ensembl_cleaned = clean_ensembl,
    aggregate_method = aggregate_method,
    global_mapping_genes = nrow(analysis$mapping),
    synonym_conflicts = if (!is.null(analysis$conflicts)) nrow(analysis$conflicts) else 0,
    alignment_strategy = alignment_strategy,
    final_genes = nrow(merged_obj),
    processing_date = as.character(Sys.Date())
  )
  
  # Save merged object
  merged_file <- file.path(output_dir, "merged_seurat_standardized.rds")
  saveRDS(merged_obj, merged_file)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Workflow Complete                                       ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("Total time: %.2f minutes\n", as.numeric(elapsed)))
    cat(sprintf("Final object: %d genes × %d cells\n", nrow(merged_obj), ncol(merged_obj)))
    cat(sprintf("Saved: %s\n\n", basename(merged_file)))
  }
  
  return(list(
    merged_obj = merged_obj,
    global_mapping = analysis$mapping,
    conflicts = analysis$conflicts,
    alignment = alignment,
    standardized_paths = standardized_paths,
    validation = validation,
    elapsed_time = elapsed
  ))
}

#===============================================================================
# Usage Example
#===============================================================================

# # After QC, you have cleaned RDS files without gene standardization
# sample_paths <- c(
#   "sample1" = "/path/to/sample1_cleaned.rds",
#   "sample2" = "/path/to/sample2_cleaned.rds",
#   "sample3" = "/path/to/sample3_cleaned.rds"
# )
# 
# # Run FIXED workflow
# result <- standardize_and_merge_workflow(
#   sample_paths = sample_paths,
#   output_dir = "/path/to/output",
#   gene_db = org.Hs.eg.db,
#   clean_ensembl = TRUE,  # Clean version numbers
#   aggregate_method = "sum",  # Sum counts for synonyms
#   alignment_strategy = "union",
#   core_threshold = 0.8,
#   merge_data = FALSE,
#   verbose = TRUE
# )
# 
# # Results with validation
# merged_seurat <- result$merged_obj
# gene_mapping <- result$global_mapping
# validation <- result$validation
