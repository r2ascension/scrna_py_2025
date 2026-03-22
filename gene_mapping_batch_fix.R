################################################################################
# Batch Fix for All Samples with Unmapped ENSG Issues
# 批量修复所有样本的 unmapped ENSG 问题
#
# Date: 2024-12-20
################################################################################

library(Seurat)
library(Matrix)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)

#===============================================================================
# Part 0: Utility Functions (SIMPLIFIED)
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

#===============================================================================
# Part 1: FIXED Global Gene Mapping (ENSG-FIRST)
#===============================================================================

build_global_gene_mapping_fixed <- function(
  all_genes_unique,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  verbose = TRUE
) {
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Building Global Gene Mapping (FIXED + ENSG-FIRST)       ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("Total genes to map: %d\n\n", length(all_genes_unique)))
  }

  # Normalize keys
  base_keys <- if (clean_ensembl) clean_ensembl_versions(all_genes_unique) else all_genes_unique
  is_ensg   <- grepl("^ENSG[0-9]+$", base_keys)

  gene_mapping <- data.frame(
    original_name   = all_genes_unique,
    official_symbol = NA_character_,
    mapping_source  = NA_character_,
    stringsAsFactors = FALSE
  )

  # ---- Layer 0: ENSEMBL FIRST for ENSG ----
  if (verbose) cat("Layer 0: ENSEMBL-first for ENSG patterns...\n")
  idx0 <- which(is_ensg)
  if (length(idx0) > 0) {
    m0 <- tryCatch(
      mapIds(gene_db,
             keys = base_keys[idx0],
             column = "SYMBOL",
             keytype = "ENSEMBL",
             multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx0))
    )
    ok0 <- !is.na(m0) & nzchar(m0)
    gene_mapping$official_symbol[idx0[ok0]] <- as.character(m0[ok0])
    gene_mapping$mapping_source[idx0[ok0]]  <- "ENSEMBL"
    if (verbose) cat(sprintf("  ✓ Mapped: %d ENSG genes\n", sum(ok0)))
  }

  # ---- Layer 1: SYMBOL (only for non-ENSG unmapped) ----
  if (verbose) cat("Layer 1: SYMBOL...\n")
  idx1 <- which(is.na(gene_mapping$official_symbol) & !is_ensg)
  if (length(idx1) > 0) {
    m1 <- tryCatch(
      mapIds(gene_db,
             keys = all_genes_unique[idx1],
             column = "SYMBOL",
             keytype = "SYMBOL",
             multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx1))
    )
    ok1 <- !is.na(m1) & nzchar(m1)
    gene_mapping$official_symbol[idx1[ok1]] <- as.character(m1[ok1])
    gene_mapping$mapping_source[idx1[ok1]]  <- "SYMBOL"
    if (verbose) cat(sprintf("  ✓ Mapped: %d genes\n", sum(ok1)))
  }

  # ---- Layer 2: ALIAS ----
  if (verbose) cat("Layer 2: ALIAS...\n")
  idx2 <- which(is.na(gene_mapping$official_symbol))
  if (length(idx2) > 0) {
    m2 <- tryCatch(
      mapIds(gene_db,
             keys = all_genes_unique[idx2],
             column = "SYMBOL",
             keytype = "ALIAS",
             multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx2))
    )
    ok2 <- !is.na(m2) & nzchar(m2)
    gene_mapping$official_symbol[idx2[ok2]] <- as.character(m2[ok2])
    gene_mapping$mapping_source[idx2[ok2]]  <- "ALIAS"
    if (verbose) cat(sprintf("  ✓ Mapped: %d genes\n", sum(ok2)))
  }

  # ---- Layer 3: ENSEMBL fallback ----
  if (verbose) cat("Layer 3: ENSEMBL (fallback for remaining)...\n")
  idx3 <- which(is.na(gene_mapping$official_symbol))
  if (length(idx3) > 0) {
    m3 <- tryCatch(
      mapIds(gene_db,
             keys = base_keys[idx3],
             column = "SYMBOL",
             keytype = "ENSEMBL",
             multiVals = "first"),
      error = function(e) rep(NA_character_, length(idx3))
    )
    ok3 <- !is.na(m3) & nzchar(m3)
    gene_mapping$official_symbol[idx3[ok3]] <- as.character(m3[ok3])
    gene_mapping$mapping_source[idx3[ok3]]  <- "ENSEMBL"
    if (verbose) cat(sprintf("  ✓ Mapped: %d genes\n", sum(ok3)))
  }

  # ---- Layer 4: UNMAPPED keep original ----
  idx4 <- which(is.na(gene_mapping$official_symbol) | !nzchar(gene_mapping$official_symbol))
  gene_mapping$official_symbol[idx4] <- gene_mapping$original_name[idx4]
  gene_mapping$mapping_source[idx4]  <- "UNMAPPED"

  if (verbose) {
    cat("\n📊 Mapping Summary:\n")
    print(table(gene_mapping$mapping_source))
    cat(sprintf("\nUnmapped rate: %.2f%%\n", 
                sum(gene_mapping$mapping_source == "UNMAPPED") / nrow(gene_mapping) * 100))
    cat("\n")
  }

  gene_mapping
}

#===============================================================================
# Part 2: FIXED Apply Mapping (Sparse Matrix Aggregation)
#===============================================================================

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

  # Subset mapping
  msub <- global_mapping[match(current_genes, global_mapping$original_name), ]
  if (any(is.na(msub$original_name))) {
    stop("Some genes in sample are missing from global mapping!")
  }

  official <- msub$official_symbol
  source_type <- msub$mapping_source

  # ⭐ CRITICAL PATCH: Re-map stuck ENSG genes
  base_keys <- clean_ensembl_versions(current_genes)
  is_ensg   <- grepl("^ENSG[0-9]+$", base_keys)
  
  needs_patch <- is_ensg & 
                 (official == current_genes) & 
                 (source_type %in% c("SYMBOL", "ALIAS", "UNMAPPED") | is.na(source_type))
  
  if (any(needs_patch)) {
    if (verbose) cat("  🔧 Applying ENSEMBL patch for stuck genes...\n")
    
    patched <- mapIds(org.Hs.eg.db,
                      keys = base_keys[needs_patch],
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
    patched <- as.character(patched)
    ok <- !is.na(patched) & nzchar(patched)
    idx <- which(needs_patch)
    official[idx[ok]] <- patched[ok]
    
    if (verbose) cat(sprintf("    ✓ Patched: %d ENSG genes\n", sum(ok)))
  }

  # ⭐ Aggregate duplicates using sparse group matrix (FAST)
  f <- factor(official, levels = unique(official))
  A <- sparseMatrix(
    i = as.integer(f),
    j = seq_along(f),
    x = 1,
    dims = c(nlevels(f), length(f))
  )

  agg <- A %*% counts_mat
  rownames(agg) <- levels(f)
  colnames(agg) <- colnames(counts_mat)

  if (aggregate_method == "mean") {
    gs <- as.numeric(table(f))
    agg <- Diagonal(x = 1 / gs) %*% agg
  } else if (aggregate_method != "sum") {
    stop("aggregate_method must be 'sum' or 'mean'")
  }

  # Replace assay
  new_assay <- tryCatch(
    CreateAssay5Object(counts = agg),
    error = function(e) CreateAssayObject(counts = agg)
  )
  seurat_obj[[assay]] <- new_assay

  if (verbose) {
    n_dup <- sum(table(official) > 1)
    cat(sprintf("  ✓ Final: %d genes (aggregated %d duplicate symbols)\n", 
                nrow(agg), n_dup))
  }

  seurat_obj
}

#===============================================================================
# Part 3: Batch Processing Workflow
#===============================================================================

batch_fix_all_samples <- function(
  sample_paths,
  output_dir,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  aggregate_method = "sum",
  save_mapping = TRUE,
  verbose = TRUE
) {
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Batch Gene Mapping Fix - All Samples                    ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    cat(sprintf("📁 Samples to process: %d\n", length(sample_paths)))
    cat(sprintf("💾 Output directory: %s\n", output_dir))
    cat(sprintf("🔧 Aggregate method: %s\n\n", aggregate_method))
  }

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ===== STEP 1: Build global gene mapping =====
  if (verbose) cat("STEP 1: Collecting all genes and building global mapping...\n\n")

  all_genes_list <- list()
  for (i in seq_along(sample_paths)) {
    sample_name <- names(sample_paths)[i]
    if (verbose) cat(sprintf("  [%d/%d] Loading: %s\n", i, length(sample_paths), sample_name))
    
    obj <- readRDS(sample_paths[i])
    all_genes_list[[sample_name]] <- get_genes_from_assay(obj, "RNA")
    rm(obj)
    gc(verbose = FALSE)
  }

  all_genes_unique <- unique(unlist(all_genes_list))
  if (verbose) cat(sprintf("\n  Total unique genes: %d\n", length(all_genes_unique)))

  # Build fixed mapping
  global_mapping <- build_global_gene_mapping_fixed(
    all_genes_unique,
    gene_db = gene_db,
    clean_ensembl = clean_ensembl,
    verbose = verbose
  )

  # Save global mapping
  if (save_mapping) {
    mapping_file <- file.path(output_dir, "global_gene_mapping_FIXED.csv")
    fwrite(global_mapping, mapping_file)
    if (verbose) cat(sprintf("💾 Saved global mapping: %s\n\n", basename(mapping_file)))
  }

  # ===== STEP 2: Apply to each sample =====
  if (verbose) {
    cat("STEP 2: Applying standardization to each sample...\n\n")
  }

  fixed_paths <- c()
  fix_report <- list()

  for (i in seq_along(sample_paths)) {
    sample_name <- names(sample_paths)[i]
    
    if (verbose) {
      cat("─────────────────────────────────────────────────────────\n")
      cat(sprintf("[%d/%d] Processing: %s\n", i, length(sample_paths), sample_name))
      cat("─────────────────────────────────────────────────────────\n")
    }

    # Load
    obj <- readRDS(sample_paths[i])
    genes_before <- nrow(obj)
    cells <- ncol(obj)

    # Check ENSG unmapped issue
    current_genes <- rownames(obj)
    is_ensg_pattern <- grepl("^ENSG[0-9]+", current_genes)
    n_ensg_before <- sum(is_ensg_pattern)

    # Apply fixed mapping
    obj_fixed <- apply_gene_mapping_fixed(
      obj,
      global_mapping,
      assay = "RNA",
      aggregate_method = aggregate_method,
      verbose = verbose
    )

    genes_after <- nrow(obj_fixed)
    
    # Check ENSG after fix
    fixed_genes <- rownames(obj_fixed)
    is_ensg_after <- grepl("^ENSG[0-9]+", fixed_genes)
    n_ensg_after <- sum(is_ensg_after)

    # Validation
    if (any(is.na(rownames(obj_fixed)))) {
      stop(sprintf("✗ FAIL: Sample %s has NA in rownames!", sample_name))
    }
    if (any(duplicated(rownames(obj_fixed)))) {
      stop(sprintf("✗ FAIL: Sample %s has duplicated rownames!", sample_name))
    }

    # Save
    output_file <- file.path(output_dir, paste0(sample_name, "_FIXED.rds"))
    saveRDS(obj_fixed, output_file)
    fixed_paths[sample_name] <- output_file

    # Report
    fix_report[[sample_name]] <- data.frame(
      sample = sample_name,
      cells = cells,
      genes_before = genes_before,
      genes_after = genes_after,
      ensg_before = n_ensg_before,
      ensg_after = n_ensg_after,
      ensg_fixed = n_ensg_before - n_ensg_after,
      output_file = basename(output_file),
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat(sprintf("  📊 Summary:\n"))
      cat(sprintf("    Genes: %d → %d\n", genes_before, genes_after))
      cat(sprintf("    ENSG unmapped: %d → %d (fixed: %d)\n", 
                  n_ensg_before, n_ensg_after, n_ensg_before - n_ensg_after))
      cat(sprintf("  ✓ Saved: %s\n\n", basename(output_file)))
    }

    rm(obj, obj_fixed)
    gc(verbose = FALSE)
  }

  # ===== STEP 3: Summary Report =====
  report_df <- do.call(rbind, fix_report)
  report_file <- file.path(output_dir, "batch_fix_report.csv")
  fwrite(report_df, report_file)

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════╗\n")
    cat("║  Batch Fix Complete                                      ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n\n")
    
    cat("📊 Overall Statistics:\n")
    cat(sprintf("  Total samples processed: %d\n", nrow(report_df)))
    cat(sprintf("  Total ENSG fixed: %d\n", sum(report_df$ensg_fixed)))
    cat(sprintf("  Processing time: %.2f minutes\n\n", as.numeric(elapsed)))
    
    cat("📋 Summary Report:\n")
    print(report_df, row.names = FALSE)
    cat("\n")
    
    cat(sprintf("💾 Report saved: %s\n", basename(report_file)))
    cat(sprintf("💾 Output directory: %s\n\n", output_dir))
  }

  return(list(
    fixed_paths = fixed_paths,
    global_mapping = global_mapping,
    report = report_df,
    elapsed_time = elapsed
  ))
}

#===============================================================================
# Usage Example: Batch Process All Samples
#===============================================================================

# 准备所有需要修复的样本路径
sample_paths <- c(
  Kerstin_B_Meyer_2021_covid_cleaned = "/home/h2048/data/R/1215/merge/cleaned_samples/Kerstin_B_Meyer_2021_covid_cleaned.rds",
  # 添加其他有 unmapped ENSG 问题的样本
  # sample2 = "/path/to/sample2.rds",
  # sample3 = "/path/to/sample3.rds"
)

# 设置输出目录
output_dir <- "/home/h2048/data/R/1215/merge/cleaned_samples_FIXED"

# 批量修复
result <- batch_fix_all_samples(
  sample_paths = sample_paths,
  output_dir = output_dir,
  gene_db = org.Hs.eg.db,
  clean_ensembl = TRUE,
  aggregate_method = "sum",
  save_mapping = TRUE,
  verbose = TRUE
)

# 查看结果
print(result$report)

# 修复后的文件路径
fixed_files <- result$fixed_paths