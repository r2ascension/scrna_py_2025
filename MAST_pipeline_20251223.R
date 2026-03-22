# ===== per_celltype_nmf_batch_analysis_v1_0_PRODUCTION.R =====
# Purpose: Batch processing for per-celltype NMF meta-program analysis
# Author: r2end
# Date: 2024-12-24
# 
# Key Features:
# - Sparse-optimized NMF (handles 80k+ cells)
# - Safe QC metrics handling
# - Robust error handling for enrichment analysis
# - Memory-efficient batch processing
#
# Memory Usage: <5GB peak for entire batch
# Runtime: ~2-3 hours for 9 cell types
# ================================================================

# ===== Load Libraries =====
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(RcppML)

# ===== Load Custom Functions =====
source('/home/h2048/script/R/tissue_comparison_analysis_20251222.R')

# ===== Configuration =====
INPUT_DIR <- "/home/h2048/data/R/1221/per_celltype_harmony_rogue/seurat_objects"
OUTPUT_BASE_DIR <- "/home/h2048/data/R/1223/per_celltype_nmf_analysis"
MSIGDB_GMT <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# NMF Parameters
NMF_RANK <- 8
NMF_N_HVG <- 2000
NMF_SEED <- 42

# Clustering Parameters
HARMONY_DIMS <- 30
CLUSTER_RESOLUTION <- 1.0

# File list
rds_files <- list.files(INPUT_DIR, pattern = "_harmony.rds$", full.names = TRUE)
cell_type_names <- gsub("_harmony.rds$", "", basename(rds_files))

cat(sprintf("Found %d RDS files to process:\n", length(rds_files)))
print(data.frame(CellType = cell_type_names, File = basename(rds_files)))

# ===== Helper Functions =====

# 1. Safe QC Metrics Calculation
ensure_qc_metrics <- function(seurat_obj) {
  # Ensure QC metrics exist and are valid BEFORE gene filtering
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Check and calculate percent.mt
  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    cat("  Creating percent.mt column...\n")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  } else {
    mt_summary <- summary(seurat_obj$percent.mt)
    if (mt_summary["Max."] == 0 || all(is.na(seurat_obj$percent.mt))) {
      cat("  ⚠️ percent.mt exists but invalid, recalculating...\n")
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    } else {
      cat(sprintf("  ✓ percent.mt exists (mean=%.2f%%)\n", mt_summary["Mean"]))
    }
  }
  
  # Check and calculate percent.rb
  if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
    cat("  Creating percent.rb column...\n")
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  } else {
    rb_summary <- summary(seurat_obj$percent.rb)
    if (rb_summary["Max."] == 0 || all(is.na(seurat_obj$percent.rb))) {
      cat("  ⚠️ percent.rb exists but invalid, recalculating...\n")
      seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
    } else {
      cat(sprintf("  ✓ percent.rb exists (mean=%.2f%%)\n", rb_summary["Mean"]))
    }
  }
  
  return(seurat_obj)
}

# 2. Recalculate Basic QC After Gene Filtering
recalc_basic_qc <- function(seurat_obj) {
  # Recalculate nCount_RNA and nFeature_RNA after subsetting features
  # This ensures QC metrics match the current gene set
  DefaultAssay(seurat_obj) <- "RNA"
  
  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )
  
  seurat_obj$nCount_RNA <- Matrix::colSums(counts_mat)
  seurat_obj$nFeature_RNA <- Matrix::colSums(counts_mat > 0)
  
  cat(sprintf("  Updated QC: mean nCount=%.0f, mean nFeature=%.0f\n",
              mean(seurat_obj$nCount_RNA), mean(seurat_obj$nFeature_RNA)))
  
  return(seurat_obj)
}

# 3. Gene Filtering Function
filter_low_quality_genes <- function(seurat_obj) {
  all_genes <- rownames(seurat_obj)
  
  # Define patterns to remove
  # Note: ENSG removal commented out - only remove if genes are truly unannotated
  # LINC/lncRNA kept for potential biological relevance in enrichment
  patterns_to_remove <- c(
    '^MT-',                               # Mitochondrial
    '^RPS|^RPL|^MRPS|^MRPL',             # Ribosomal
    '^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$', # Ribosomal pseudogenes
    # '^ENSG[0-9]+',                      # Commented: may remove valid genes
    '^(AC|AL|AP|BX|Z)[0-9]+\\.',         # Unannotated contigs
    '^RP[0-9]+-',                         # RP transcripts
    '^CTD-|^CTB-|^CTC-',                 # Clone-based IDs
    # '^LINC[0-9]+',                      # Commented: lncRNA may be relevant
    # '-AS[0-9]+$',                       # Commented: antisense may be relevant
    '-OT[0-9]+$',                         # Overlapping transcripts
    '^LOC[0-9]+'                          # LOC IDs
  )
  
  genes_to_remove <- c()
  for (pattern in patterns_to_remove) {
    pattern_genes <- grep(pattern, all_genes, value = TRUE)
    if (length(pattern_genes) > 0) {
      cat(sprintf("  Pattern '%s': %d genes\n", pattern, length(pattern_genes)))
      genes_to_remove <- c(genes_to_remove, pattern_genes)
    }
  }
  genes_to_remove <- unique(genes_to_remove)
  
  cat(sprintf("  Total genes to remove: %d\n", length(genes_to_remove)))
  
  genes_to_keep <- setdiff(all_genes, genes_to_remove)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)
  
  # Safe LayerData extraction with fallback
  DefaultAssay(seurat_obj) <- "RNA"
  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) {
      cat("  LayerData failed, using GetAssayData fallback\n")
      GetAssayData(seurat_obj, slot = "counts")
    }
  )
  
  # Remove genes expressed in < 3 cells
  gene_ncells <- Matrix::rowSums(counts_mat > 0)
  keep_genes <- names(gene_ncells[gene_ncells >= 3])
  
  before_n <- nrow(seurat_obj)
  seurat_obj <- subset(seurat_obj, features = keep_genes)
  after_n <- nrow(seurat_obj)
  
  cat(sprintf("  Gene filtering: %d → %d (removed %d)\n", 
              before_n, after_n, before_n - after_n))
  
  return(seurat_obj)
}

# 4. NMF Analysis Function (SPARSE-OPTIMIZED)
run_nmf_analysis <- function(seurat_obj, rank = 8, n_hvg = 2000, seed = 42) {
  cat("  Running NMF analysis (sparse-optimized)...\n")
  
  # Get HVG list and validate count
  all_hvg <- VariableFeatures(seurat_obj)
  n_available_hvg <- length(all_hvg)
  
  if (n_available_hvg < n_hvg) {
    cat(sprintf("    ⚠️ Only %d HVGs available (requested %d), using all\n", 
                n_available_hvg, n_hvg))
    hvg_genes <- all_hvg
  } else {
    hvg_genes <- all_hvg[1:n_hvg]
  }
  
  hvg_genes <- hvg_genes[hvg_genes %in% rownames(seurat_obj)]
  cat(sprintf("    Using %d HVGs for NMF\n", length(hvg_genes)))
  
  # Safe data extraction with fallback
  DefaultAssay(seurat_obj) <- "RNA"
  data_mat <- tryCatch(
    LayerData(seurat_obj, layer = "data")[hvg_genes, , drop = FALSE],
    error = function(e) {
      cat("    LayerData failed, using GetAssayData fallback\n")
      GetAssayData(seurat_obj, slot = "data")[hvg_genes, , drop = FALSE]
    }
  )
  
  # ⭐ CRITICAL: Keep sparse throughout expm1 transformation
  # Direct expm1(dgCMatrix) may densify in some R environments
  # Operate on @x slot to guarantee sparsity
  data_for_nmf <- data_mat
  
  if (inherits(data_for_nmf, "dgCMatrix")) {
    cat("    ✓ Input is sparse, applying expm1 to non-zero entries only\n")
    data_for_nmf@x <- expm1(data_for_nmf@x)
    data_for_nmf@x[data_for_nmf@x < 0] <- 0
  } else {
    cat("    ⚠️ Input is dense, converting to sparse...\n")
    data_for_nmf <- as(data_for_nmf, "dgCMatrix")
    data_for_nmf@x <- expm1(data_for_nmf@x)
    data_for_nmf@x[data_for_nmf@x < 0] <- 0
  }
  
  cat(sprintf("    NMF input: %d genes × %d cells (%.2f MB)\n", 
              nrow(data_for_nmf), ncol(data_for_nmf),
              object.size(data_for_nmf) / 1024^2))
  
  # Run NMF (RcppML accepts sparse matrices)
  set.seed(seed)
  nmf_result <- RcppML::nmf(
    data_for_nmf,
    k = rank,
    tol = 1e-4,
    maxit = 100,
    verbose = FALSE,
    seed = seed
  )
  
  # Extract results (RcppML returns list)
  w_matrix <- nmf_result$w
  h_matrix <- nmf_result$h
  
  cat(sprintf("    W matrix: %d genes × %d factors\n", nrow(w_matrix), ncol(w_matrix)))
  cat(sprintf("    H matrix: %d factors × %d cells\n", nrow(h_matrix), ncol(h_matrix)))
  
  # Add NMF scores to metadata using AddMetaData (safer than cbind)
  h_matrix_t <- t(h_matrix)
  colnames(h_matrix_t) <- paste0("NMF_", 1:rank)
  rownames(h_matrix_t) <- colnames(seurat_obj)
  
  seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(h_matrix_t))
  
  # Store NMF results in misc slot
  seurat_obj@misc$nmf_result <- list(
    W = w_matrix,
    H = h_matrix,
    rank = rank,
    hvg_genes = hvg_genes,
    seed = seed
  )
  
  cat(sprintf("    ✓ NMF completed: %d meta-programs identified\n", rank))
  
  return(seurat_obj)
}

# 5. Extract Top Genes per NMF Factor
extract_nmf_top_genes <- function(seurat_obj, top_n = 50) {
  nmf_result <- seurat_obj@misc$nmf_result
  w_matrix <- nmf_result$W
  hvg_genes <- nmf_result$hvg_genes
  
  rownames(w_matrix) <- hvg_genes
  
  top_genes_list <- list()
  for (i in 1:ncol(w_matrix)) {
    factor_loadings <- w_matrix[, i]
    top_genes <- names(sort(factor_loadings, decreasing = TRUE)[1:top_n])
    top_genes_list[[paste0("NMF_", i)]] <- top_genes
  }
  
  return(top_genes_list)
}

# 6. Safe Enrichment Wrapper
run_enrichment_safe <- function(seurat_obj, cluster_col, msigdb_gmt) {
  cat("  Running enrichment analysis (with error handling)...\n")
  
  tryCatch({
    run_one_vs_rest_enrichment(
      seurat_obj,
      group_col = cluster_col,
      msigdb_gmt_file = msigdb_gmt,
      run_go = TRUE,
      run_gsea = TRUE
    )
    cat("  ✓ Enrichment analysis completed\n")
  }, error = function(e) {
    cat(sprintf("  ⚠️ run_one_vs_rest_enrichment failed: %s\n", e$message))
    cat("  Trying GO-only enrichment...\n")
    
    tryCatch({
      run_one_vs_rest_enrichment(
        seurat_obj,
        group_col = cluster_col,
        msigdb_gmt_file = msigdb_gmt,
        run_go = TRUE,
        run_gsea = FALSE
      )
      cat("  ✓ GO enrichment completed (GSEA skipped)\n")
    }, error = function(e2) {
      cat(sprintf("  ⚠️ All enrichment methods failed: %s\n", e2$message))
      cat("  Continuing without enrichment analysis...\n")
    })
  })
}

# 7. NMF Factor GO Enrichment
run_nmf_factor_enrichment <- function(nmf_top_genes, output_dir) {
  cat("  Running NMF factor enrichment...\n")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cat("  ⚠️ clusterProfiler or org.Hs.eg.db not available, skipping\n")
    return(NULL)
  }
  
  for (i in seq_along(nmf_top_genes)) {
    factor_name <- names(nmf_top_genes)[i]
    factor_genes <- nmf_top_genes[[i]]
    
    cat(sprintf("    Analyzing %s (%d genes)...\n", factor_name, length(factor_genes)))
    
    tryCatch({
      ego <- clusterProfiler::enrichGO(
        gene = factor_genes,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(
          as.data.frame(ego),
          file.path(output_dir, paste0(factor_name, "_GO_BP.csv")),
          row.names = FALSE
        )
        cat(sprintf("      ✓ Found %d enriched GO terms\n", nrow(as.data.frame(ego))))
      } else {
        cat(sprintf("      No significant GO terms for %s\n", factor_name))
      }
      
    }, error = function(e) {
      cat(sprintf("      ⚠️ GO enrichment failed for %s: %s\n", factor_name, e$message))
    })
  }
}

# ===== Main Processing Loop =====
for (i in seq_along(rds_files)) {
  cell_type <- cell_type_names[i]
  rds_file <- rds_files[i]
  
  cat(sprintf("\n========== Processing %s (%d/%d) ==========\n", 
              cell_type, i, length(rds_files)))
  
  # Create output directory
  output_dir <- file.path(OUTPUT_BASE_DIR, cell_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===== 1. Load Data =====
  cat("Step 1: Loading RDS...\n")
  seurat_obj <- readRDS(rds_file)
  cat(sprintf("  Loaded: %d cells × %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))
  
  # ===== 2. Ensure QC Metrics (BEFORE gene filtering!) =====
  cat("Step 2: Ensuring QC metrics...\n")
  seurat_obj <- ensure_qc_metrics(seurat_obj)
  
  # ===== 3. Gene Filtering =====
  cat("Step 3: Gene filtering...\n")
  seurat_obj <- filter_low_quality_genes(seurat_obj)
  
  # ===== 3b. Recalculate QC After Filtering =====
  seurat_obj <- recalc_basic_qc(seurat_obj)
  
  # ===== 4. Standard Processing =====
  cat("Step 4: Normalization and PCA...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
  
  # ===== 5. Harmony Batch Correction =====
  cat("Step 5: Harmony batch correction...\n")
  seurat_obj <- RunHarmony(
    seurat_obj,
    group.by.vars = "sample",
    theta = 2,
    lambda = 1,
    sigma = 0.1,
    nclust = 30,
    max_iter = 20,
    dims = 1:HARMONY_DIMS,
    verbose = FALSE
  )
  
  # ===== 6. UMAP and Clustering =====
  cat("Step 6: UMAP and clustering...\n")
  set.seed(NMF_SEED)  # Set seed for reproducibility
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "harmony",
    dims = 1:HARMONY_DIMS,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "correlation",
    seed.use = NMF_SEED,  # Explicit seed for UMAP
    verbose = FALSE
  )
  
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = "harmony",
    dims = 1:HARMONY_DIMS,
    k.param = 30,
    verbose = FALSE
  )
  
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = CLUSTER_RESOLUTION,
    algorithm = 4,
    random.seed = NMF_SEED,  # Explicit seed for clustering
    verbose = FALSE
  )
  
  cluster_col <- paste0("RNA_snn_res.", CLUSTER_RESOLUTION)
  n_clusters <- length(unique(seurat_obj@meta.data[[cluster_col]]))
  cat(sprintf("  Identified %d clusters\n", n_clusters))
  
  # ===== 7. Find Markers =====
  # Note: Markers are computed on log-normalized RNA (not Harmony-corrected)
  # This is standard practice - clusters from Harmony space, markers from original expression
  cat("Step 7: Finding cluster markers...\n")
  Idents(seurat_obj) <- cluster_col
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  write.csv(markers, file.path(output_dir, "all_markers_complete.csv"), row.names = FALSE)
  
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
  write.csv(top_markers, file.path(output_dir, "all_markers_top10.csv"), row.names = FALSE)
  
  # ===== 8. Enrichment Analysis =====
  cat("Step 8: Running enrichment analysis...\n")
  run_enrichment_safe(seurat_obj, cluster_col, MSIGDB_GMT)
  
  # ===== 9. NMF Analysis (SPARSE-OPTIMIZED) =====
  cat("Step 9: NMF meta-program identification...\n")
  seurat_obj <- run_nmf_analysis(
    seurat_obj,
    rank = NMF_RANK,
    n_hvg = NMF_N_HVG,
    seed = NMF_SEED
  )
  
  # Extract top genes per NMF factor
  nmf_top_genes <- extract_nmf_top_genes(seurat_obj, top_n = 50)
  saveRDS(nmf_top_genes, file.path(output_dir, "nmf_top_genes.rds"))
  
  nmf_genes_df <- data.frame(
    Factor = rep(names(nmf_top_genes), each = 50),
    Rank = rep(1:50, length(nmf_top_genes)),
    Gene = unlist(nmf_top_genes)
  )
  write.csv(nmf_genes_df, file.path(output_dir, "nmf_top_genes.csv"), row.names = FALSE)
  
  # ===== 10. NMF Factor Enrichment =====
  cat("Step 10: NMF factor enrichment analysis...\n")
  run_nmf_factor_enrichment(nmf_top_genes, file.path(output_dir, "nmf_enrichment"))
  
  # ===== 11. Visualization =====
  cat("Step 11: Generating visualizations...\n")
  
  # 11.1 Basic UMAP
  pdf(file.path(output_dir, "01_umap_overview.pdf"), width = 15, height = 5)
  p1 <- DimPlot(seurat_obj, group.by = cluster_col, label = TRUE, raster = TRUE) +
    ggtitle(paste0(cell_type, " - Clusters"))
  p2 <- DimPlot(seurat_obj, group.by = "dataset", raster = TRUE) +
    ggtitle("Dataset")
  p3 <- DimPlot(seurat_obj, group.by = "tissue", raster = TRUE) +
    ggtitle("Tissue")
  print(p1 + p2 + p3)
  dev.off()
  
  # 11.2 NMF Factor UMAP
  pdf(file.path(output_dir, "02_nmf_factors_umap.pdf"), width = 15, height = 10)
  nmf_cols <- paste0("NMF_", 1:NMF_RANK)
  p_nmf <- FeaturePlot(
    seurat_obj,
    features = nmf_cols,
    ncol = 4,
    raster = TRUE
  )
  print(p_nmf)
  dev.off()
  
  # 11.3 NMF Factor Violin Plot
  pdf(file.path(output_dir, "03_nmf_factors_by_cluster.pdf"), width = 15, height = 12)
  for (nmf_idx in 1:NMF_RANK) {
    p <- VlnPlot(
      seurat_obj,
      features = paste0("NMF_", nmf_idx),
      group.by = cluster_col,
      pt.size = 0
    ) + ggtitle(paste0("NMF Factor ", nmf_idx))
    print(p)
  }
  dev.off()
  
  # 11.4 Marker Heatmap
  marker_genes <- unique(top_markers$gene)
  marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]
  
  if (length(marker_genes) > 0) {
    seurat_obj <- ScaleData(seurat_obj, features = marker_genes, verbose = FALSE)
    
    pdf(file.path(output_dir, "04_marker_heatmap.pdf"), width = 12, height = 10)
    p <- DoHeatmap(
      seurat_obj,
      features = marker_genes,
      group.by = cluster_col,
      raster = TRUE
    ) + NoLegend()
    print(p)
    dev.off()
  }
  
  # 11.5 QC Metrics
  pdf(file.path(output_dir, "05_qc_metrics.pdf"), width = 12, height = 8)
  p_qc <- VlnPlot(
    seurat_obj,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
    group.by = cluster_col,
    pt.size = 0,
    ncol = 2
  )
  print(p_qc)
  dev.off()
  
  # 11.6 QC Detailed
  pdf(file.path(output_dir, "06_qc_by_cluster_detailed.pdf"), width = 12, height = 6, onefile = TRUE)
  qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
  for (f in qc_feats) {
    p <- VlnPlot(seurat_obj, features = f, group.by = cluster_col, pt.size = 0)
    print(p)
  }
  dev.off()
  
  # 11.7 NMF Heatmap by Cluster
  pdf(file.path(output_dir, "07_nmf_heatmap_by_cluster.pdf"), width = 10, height = 8)
  
  nmf_data <- seurat_obj@meta.data[, c(cluster_col, paste0("NMF_", 1:NMF_RANK))]
  colnames(nmf_data)[1] <- "cluster"
  
  nmf_avg <- nmf_data %>%
    group_by(cluster) %>%
    summarise(across(starts_with("NMF_"), mean))
  
  nmf_matrix <- as.matrix(nmf_avg[, -1])
  rownames(nmf_matrix) <- nmf_avg$cluster
  
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
      t(nmf_matrix),
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      scale = "row",
      main = paste0(cell_type, " - NMF Factor Enrichment by Cluster")
    )
  }
  dev.off()
  
  # ===== 12. Save Final Object =====
  cat("Step 12: Saving results...\n")
  output_rds <- file.path(output_dir, paste0(cell_type, "_analyzed_nmf.rds"))
  saveRDS(seurat_obj, output_rds)
  cat(sprintf("  Saved: %s\n", basename(output_rds)))
  
  # Save analysis summary
  summary_info <- list(
    cell_type = cell_type,
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    n_clusters = n_clusters,
    nmf_rank = NMF_RANK,
    processing_date = Sys.time(),
    cluster_column = cluster_col,
    notes = "Clusters from Harmony space; markers from log-normalized RNA"
  )
  saveRDS(summary_info, file.path(output_dir, "analysis_summary.rds"))
  
  # Clean up memory
  rm(seurat_obj, markers, top_markers, nmf_top_genes)
  gc()
  
  cat(sprintf("✓ Completed %s\n", cell_type))
}

cat("\n========== All Analyses Complete! ==========\n")
cat(sprintf("Results saved in: %s\n", OUTPUT_BASE_DIR))

# ===== Generate Summary Report =====
cat("\nGenerating summary report...\n")

summary_df <- data.frame()
for (ct in cell_type_names) {
  summary_file <- file.path(OUTPUT_BASE_DIR, ct, "analysis_summary.rds")
  if (file.exists(summary_file)) {
    info <- readRDS(summary_file)
    summary_df <- rbind(summary_df, data.frame(
      CellType = info$cell_type,
      N_Cells = info$n_cells,
      N_Genes = info$n_genes,
      N_Clusters = info$n_clusters,
      NMF_Rank = info$nmf_rank,
      ProcessingDate = as.character(info$processing_date)
    ))
  }
}

write.csv(summary_df, file.path(OUTPUT_BASE_DIR, "analysis_summary_all.csv"), row.names = FALSE)
cat("\n✓ Summary report saved: analysis_summary_all.csv\n")
cat("\n=== Pipeline Complete ===\n")