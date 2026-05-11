#!/usr/bin/env Rscript
# ===== h5ad_mast_go_kegg_single_v1_1_OPTIMIZED.R =====
# Purpose: Optimized version with parallel processing
# Changes from v1.0:
#   - Parallel FindAllMarkers
#   - Optimized gene filtering (compiled regex)
#   - Sparse matrix operations optimization
#   - Pre-filtered low-quality data before heavy computation
#   - Option to skip MAST for speed (use faster methods)
# Author: r2end
# Date: 2024-12-28
# =====================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(SCNT)
  library(data.table)
  library(harmony)
  library(future) # For parallel processing
  library(future.apply)
})

# ===== Parallel Setup =====
# Set number of cores (留1-2个核给系统)
N_CORES <- 8
cat(sprintf("Setting up parallel processing with %d cores\n", N_CORES))

plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 10 * 1024^3) # 10GB limit per worker

# ===== Python / Conda Environment =====
use_condaenv("bbknn_env", required = TRUE)
py_config()

# ===== Load Custom Functions =====
source("/home/h2048/script/R/tissue_comparison_analysis_20251222.R")

# ===== Configuration =====
H5AD_FILE <- "/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2/Ciliated/batch_aware/cnmf_analysis_k40/Ciliated_with_cnmf_k40.h5ad"

OUTPUT_DIR <- "/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2/Ciliated/h5ad_single_mast_go_kegg_optimized"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTPUT_DIR)

MSIGDB_GMT <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

CLUSTER_COL <- "RNA_snn_res.1"

# ===== Performance Options =====
RUN_FIND_ALL_MARKERS <- TRUE

# Speed vs Accuracy trade-off for DE
# Options: "wilcox" (fastest), "t" (fast), "MAST" (slow but accurate)
DE_METHOD <- "wilcox" # ⚡ Changed from MAST for 10-20x speedup

# For very large datasets, consider subsampling for marker finding
USE_SUBSAMPLING_FOR_MARKERS <- FALSE
MAX_CELLS_PER_CLUSTER <- 500 # Only if USE_SUBSAMPLING_FOR_MARKERS = TRUE

# Gene filtering
MIN_CELLS_PER_GENE <- 3
MIN_GENES_PER_CELL <- 200 # Pre-filter low-quality cells

# =====================================================
# Optimized Helper Functions
# =====================================================

ensure_qc_metrics <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"

  # Get counts matrix once
  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  # Calculate all QC metrics at once to avoid repeated matrix access
  if (
    !"percent.mt" %in% colnames(seurat_obj@meta.data) ||
      max(seurat_obj$percent.mt, na.rm = TRUE) == 0
  ) {
    cat("  Creating percent.mt...\n")
    mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
    if (length(mt_genes) > 0) {
      seurat_obj[["percent.mt"]] <- (Matrix::colSums(counts_mat[
        mt_genes,
        ,
        drop = FALSE
      ]) /
        Matrix::colSums(counts_mat)) *
        100
    } else {
      seurat_obj[["percent.mt"]] <- 0
    }
  }

  if (
    !"percent.rb" %in% colnames(seurat_obj@meta.data) ||
      max(seurat_obj$percent.rb, na.rm = TRUE) == 0
  ) {
    cat("  Creating percent.rb...\n")
    rb_genes <- grep("^RP[SL]", rownames(seurat_obj), value = TRUE)
    if (length(rb_genes) > 0) {
      seurat_obj[["percent.rb"]] <- (Matrix::colSums(counts_mat[
        rb_genes,
        ,
        drop = FALSE
      ]) /
        Matrix::colSums(counts_mat)) *
        100
    } else {
      seurat_obj[["percent.rb"]] <- 0
    }
  }

  seurat_obj
}

recalc_basic_qc <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"

  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  # Use sparse matrix operations
  seurat_obj$nCount_RNA <- Matrix::colSums(counts_mat)
  seurat_obj$nFeature_RNA <- Matrix::colSums(counts_mat > 0)

  cat(sprintf(
    "  Updated QC: mean nCount=%.0f, mean nFeature=%.0f\n",
    mean(seurat_obj$nCount_RNA),
    mean(seurat_obj$nFeature_RNA)
  ))
  seurat_obj
}

# ⚡ Optimized gene filtering with compiled regex
filter_low_quality_genes_fast <- function(seurat_obj, min_cells = 3) {
  DefaultAssay(seurat_obj) <- "RNA"
  all_genes <- rownames(seurat_obj)

  # Combine patterns into single regex for efficiency
  combined_pattern <- paste(
    c(
      "^MT-",
      "^RPS|^RPL|^MRPS|^MRPL",
      "^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$",
      "^(AC|AL|AP|BX|Z)[0-9]+\\.",
      "^RP[0-9]+-",
      "^CTD-|^CTB-|^CTC-",
      "-OT[0-9]+$",
      "^LOC[0-9]+"
    ),
    collapse = "|"
  )

  # Single grep call instead of multiple
  genes_to_remove <- grep(combined_pattern, all_genes, value = TRUE)
  cat(sprintf(
    "  Genes to remove by pattern: %d\n",
    length(genes_to_remove)
  ))

  genes_to_keep <- setdiff(all_genes, genes_to_remove)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)

  # Optimize min_cells filtering
  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  # Use sparse matrix operations
  gene_ncells <- Matrix::rowSums(counts_mat > 0)
  keep_genes <- names(gene_ncells[gene_ncells >= min_cells])

  before_n <- nrow(seurat_obj)
  seurat_obj <- subset(seurat_obj, features = keep_genes)
  after_n <- nrow(seurat_obj)

  cat(sprintf(
    "  Gene filtering: %d → %d (removed %d)\n",
    before_n,
    after_n,
    before_n - after_n
  ))
  seurat_obj
}

pick_cluster_col <- function(seurat_obj, cluster_col) {
  md <- colnames(seurat_obj@meta.data)
  if (cluster_col %in% md) {
    return(cluster_col)
  }

  candidates <- c(
    "leiden",
    "louvain",
    "clusters",
    "scanpy_leiden",
    "scanpy_louvain",
    "seurat_clusters",
    "RNA_snn_res.0.5",
    "RNA_snn_res.0.8",
    "RNA_snn_res.1",
    "scanvi_predictions",
    "celltypist_pred",
    "manual_annotations"
  )
  hit <- candidates[candidates %in% md]
  if (length(hit) > 0) {
    cat(sprintf(
      "  ⚠️ CLUSTER_COL '%s' not found. Auto-using '%s'\n",
      cluster_col,
      hit[1]
    ))
    return(hit[1])
  }

  stop(sprintf(
    "Cluster column '%s' not found. Available columns:\n%s",
    cluster_col,
    paste(head(md, 60), collapse = ", ")
  ))
}

ensure_umap <- function(seurat_obj, reduction_name = "umap", seed = 42) {
  if (reduction_name %in% names(seurat_obj@reductions)) {
    return(seurat_obj)
  }

  cat(sprintf(
    "  ⚠️ Reduction '%s' not found. Computing PCA+UMAP...\n",
    reduction_name
  ))
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    nfeatures = 2000,
    verbose = FALSE
  )
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

  set.seed(seed)
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca",
    dims = 1:30,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "correlation",
    seed.use = seed,
    verbose = FALSE
  )
  seurat_obj
}

# ⚡ Parallel FindAllMarkers wrapper
find_all_markers_parallel <- function(
  seurat_obj,
  test_use = "wilcox",
  only_pos = TRUE,
  min_pct = 0.25,
  logfc_threshold = 0.25,
  latent_vars = NULL
) {
  clusters <- levels(Idents(seurat_obj))
  cat(sprintf(
    "Finding markers for %d clusters in parallel...\n",
    length(clusters)
  ))

  # Split by cluster and run in parallel
  marker_list <- future_lapply(
    clusters,
    function(cluster_id) {
      tryCatch(
        {
          if (test_use == "MAST" && !is.null(latent_vars)) {
            FindMarkers(
              seurat_obj,
              ident.1 = cluster_id,
              only.pos = only_pos,
              min.pct = min_pct,
              logfc.threshold = logfc_threshold,
              test.use = "MAST",
              latent.vars = latent_vars,
              verbose = FALSE
            )
          } else {
            FindMarkers(
              seurat_obj,
              ident.1 = cluster_id,
              only.pos = only_pos,
              min.pct = min_pct,
              logfc.threshold = logfc_threshold,
              test.use = test_use,
              verbose = FALSE
            )
          }
        },
        error = function(e) {
          cat(sprintf("  ⚠️ Error in cluster %s: %s\n", cluster_id, e$message))
          return(NULL)
        }
      )
    },
    future.seed = TRUE
  )

  # Combine results
  names(marker_list) <- clusters
  marker_list <- marker_list[!sapply(marker_list, is.null)]

  markers_df <- bind_rows(lapply(names(marker_list), function(cid) {
    df <- marker_list[[cid]]
    df$cluster <- cid
    df$gene <- rownames(df)
    df
  }))

  return(markers_df)
}

run_enrichment_safe <- function(seurat_obj, cluster_col, msigdb_gmt) {
  cat("Running enrichment analysis (with error handling)...\n")

  tryCatch(
    {
      run_one_vs_rest_enrichment(
        seurat_obj,
        group_col = cluster_col,
        msigdb_gmt_file = msigdb_gmt,
        run_go = TRUE,
        run_gsea = TRUE
      )
      cat("✓ Enrichment analysis completed\n")
    },
    error = function(e) {
      cat(sprintf("⚠️ run_one_vs_rest_enrichment failed: %s\n", e$message))
      cat("Trying GO-only enrichment...\n")

      tryCatch(
        {
          run_one_vs_rest_enrichment(
            seurat_obj,
            group_col = cluster_col,
            msigdb_gmt_file = msigdb_gmt,
            run_go = TRUE,
            run_gsea = FALSE
          )
          cat("✓ GO enrichment completed (GSEA skipped)\n")
        },
        error = function(e2) {
          cat(sprintf("⚠️ All enrichment methods failed: %s\n", e2$message))
          cat("Continuing without enrichment analysis...\n")
        }
      )
    }
  )
}

# =====================================================
# Main Pipeline
# =====================================================

cat("=== OPTIMIZED PIPELINE START ===\n")
start_time <- Sys.time()

cat("\n[1/9] Reading h5ad...\n")
t1 <- Sys.time()
seurat_obj <- GetSeurat(h5ad_path = H5AD_FILE, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"
cat(sprintf(
  "  Loaded: %d cells × %d genes (%.1f sec)\n",
  ncol(seurat_obj),
  nrow(seurat_obj),
  as.numeric(difftime(Sys.time(), t1, units = "secs"))
))

cat("\n[2/9] Ensuring QC metrics...\n")
t2 <- Sys.time()
seurat_obj <- ensure_qc_metrics(seurat_obj)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t2, units = "secs"))
))

# ⚡ Pre-filter low-quality cells BEFORE heavy computation
cat("\n[2.5/9] Pre-filtering low-quality cells...\n")
t25 <- Sys.time()
before_cells <- ncol(seurat_obj)
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= MIN_GENES_PER_CELL &
    percent.mt < 20 &
    nCount_RNA > 0
)
after_cells <- ncol(seurat_obj)
cat(sprintf(
  "  Cells: %d → %d (removed %d low-quality, %.1f sec)\n",
  before_cells,
  after_cells,
  before_cells - after_cells,
  as.numeric(difftime(Sys.time(), t25, units = "secs"))
))

cat("\n[3/9] Gene filtering (optimized)...\n")
t3 <- Sys.time()
seurat_obj <- filter_low_quality_genes_fast(
  seurat_obj,
  min_cells = MIN_CELLS_PER_GENE
)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t3, units = "secs"))
))

cat("\n[4/9] Recalculate basic QC...\n")
t4 <- Sys.time()
seurat_obj <- recalc_basic_qc(seurat_obj)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t4, units = "secs"))
))

cat("\n[5/9] Normalize data...\n")
t5 <- Sys.time()
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  nfeatures = 2000,
  verbose = FALSE
)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t5, units = "secs"))
))

cat("\n[6/9] Determine cluster column...\n")
cluster_col <- pick_cluster_col(seurat_obj, CLUSTER_COL)
cat(sprintf("  Using: %s\n", cluster_col))
cluster_table <- sort(
  table(seurat_obj@meta.data[[cluster_col]]),
  decreasing = TRUE
)
print(cluster_table)

# ⚡ Optional: Downsample for marker finding
if (USE_SUBSAMPLING_FOR_MARKERS && any(cluster_table > MAX_CELLS_PER_CLUSTER)) {
  cat(sprintf(
    "\n  ⚡ Downsampling large clusters to max %d cells for marker finding...\n",
    MAX_CELLS_PER_CLUSTER
  ))

  cells_to_keep <- unlist(lapply(names(cluster_table), function(cid) {
    cells_in_cluster <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[cluster_col]] == cid
    ]
    if (length(cells_in_cluster) > MAX_CELLS_PER_CLUSTER) {
      sample(cells_in_cluster, MAX_CELLS_PER_CLUSTER)
    } else {
      cells_in_cluster
    }
  }))

  seurat_obj_subset <- subset(seurat_obj, cells = cells_to_keep)
  cat(sprintf("  Subset: %d cells\n", ncol(seurat_obj_subset)))
} else {
  seurat_obj_subset <- seurat_obj
}

# ⚡ Parallel FindAllMarkers
if (isTRUE(RUN_FIND_ALL_MARKERS)) {
  cat("\n[7/9] FindAllMarkers (PARALLEL)...\n")
  t7 <- Sys.time()
  Idents(seurat_obj_subset) <- cluster_col

  if (DE_METHOD == "MAST" && requireNamespace("MAST", quietly = TRUE)) {
    cat("  Using MAST (slow but accurate)...\n")
    markers <- find_all_markers_parallel(
      seurat_obj_subset,
      test_use = "MAST",
      latent_vars = intersect(
        c("nCount_RNA", "percent.mt"),
        colnames(seurat_obj_subset@meta.data)
      )
    )
  } else {
    if (DE_METHOD == "MAST") {
      cat("  ⚠️ MAST not available; using Wilcoxon.\n")
    } else {
      cat(sprintf("  Using %s (fast)...\n", DE_METHOD))
    }
    markers <- find_all_markers_parallel(
      seurat_obj_subset,
      test_use = DE_METHOD
    )
  }

  t7_elapsed <- as.numeric(difftime(Sys.time(), t7, units = "secs"))
  cat(sprintf(
    "  Found %d markers (%.1f sec, %.1f markers/sec)\n",
    nrow(markers),
    t7_elapsed,
    nrow(markers) / t7_elapsed
  ))

  write.csv(markers, "all_markers_complete.csv", row.names = FALSE)

  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
  write.csv(top_markers, "all_markers_top10.csv", row.names = FALSE)
}

# Use full object for enrichment
cat("\n[8/9] Enrichment analysis...\n")
t8 <- Sys.time()
run_enrichment_safe(seurat_obj, cluster_col, MSIGDB_GMT)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t8, units = "secs"))
))

# Plots
cat("\n[9/9] Generating plots...\n")
t9 <- Sys.time()
seurat_obj <- ensure_umap(seurat_obj, reduction_name = "umap", seed = 42)

tryCatch(
  {
    pdf("01_umap_overview.pdf", width = 14, height = 5)
    p1 <- DimPlot(
      seurat_obj,
      group.by = cluster_col,
      label = TRUE,
      raster = TRUE
    ) +
      ggtitle(cluster_col)
    p_list <- list(p1)
    if ("dataset" %in% colnames(seurat_obj@meta.data)) {
      p_list[[length(p_list) + 1]] <- DimPlot(
        seurat_obj,
        group.by = "dataset",
        raster = TRUE
      ) +
        ggtitle("dataset")
    }
    if ("tissue" %in% colnames(seurat_obj@meta.data)) {
      p_list[[length(p_list) + 1]] <- DimPlot(
        seurat_obj,
        group.by = "tissue",
        raster = TRUE
      ) +
        ggtitle("tissue")
    }
    print(wrap_plots(p_list, nrow = 1))
    dev.off()
    cat("  ✓ UMAP overview saved\n")
  },
  error = function(e) {
    cat(sprintf("  ⚠️ UMAP plot failed: %s\n", e$message))
    tryCatch(dev.off(), error = function(e) NULL)
  }
)

tryCatch(
  {
    pdf("02_qc_violin.pdf", width = 12, height = 8)
    print(VlnPlot(
      seurat_obj,
      features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
      group.by = cluster_col,
      pt.size = 0,
      ncol = 2
    ))
    dev.off()
    cat("  ✓ QC violin saved\n")
  },
  error = function(e) {
    cat(sprintf("  ⚠️ QC violin failed: %s\n", e$message))
    tryCatch(dev.off(), error = function(e) NULL)
  }
)

cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t9, units = "secs"))
))

cat("\n[Final] Saving object...\n")
saveRDS(seurat_obj, "seurat_from_h5ad_for_enrichment.rds")

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n=== PIPELINE COMPLETE ===\n")
cat(sprintf("Total time: %.1f minutes\n", total_time))
cat(sprintf("Output dir: %s\n", OUTPUT_DIR))
