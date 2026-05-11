#!/usr/bin/env Rscript
# ===== h5ad_mast_go_kegg_single_v1_3_PRODUCTION.R =====
# Purpose: Production-ready version with code review fixes
# Major fixes from v1.2:
#   1. QC metrics consistency: separate raw vs filtered metrics
#   2. Remove MT/RP genes from matrix (standard practice)
#   3. Gene ID sanity check before enrichment
#   4. DietSeurat for parallel marker finding (memory safety)
#   5. ScaleData only on HVG
#   6. Empty result protection
# Author: r2end
# Date: 2024-12-30
# ===================================================================

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
  library(future)
  library(future.apply)
})

# ===== Parallel Setup =====
N_CORES <- 8
cat(sprintf("Parallel cores available: %d\n", N_CORES))
cat("Default mode: SEQUENTIAL (parallel only for FindAllMarkers)\n")

# ===== Python / Conda Environment =====
use_condaenv("bbknn_env", required = TRUE)
py_config()

# ===== Load Custom Functions =====
source("/home/h2048/script/R/tissue_comparison_analysis_20251222.R")

plan("sequential")
options(future.globals.maxSize = 20 * 1024^3)


# ===== Configuration =====
H5AD_FILE <- "/home/h2048/data/R/1228/basal/basal_filtered_20251228_2.rds"

OUTPUT_DIR <- "/home/h2048/data/R/1228/basal/h5ad_single_mast_go_kegg_v1_3"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTPUT_DIR)

MSIGDB_GMT <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

CLUSTER_COL <- "RNA_snn_res.1"

# ===== Performance Options =====
RUN_FIND_ALL_MARKERS <- TRUE

DE_METHOD <- "wilcox" # Options: "wilcox" (fastest), "t", "MAST"

USE_SUBSAMPLING_FOR_MARKERS <- FALSE
MAX_CELLS_PER_CLUSTER <- 500

# Gene filtering
MIN_CELLS_PER_GENE <- 3
MIN_GENES_PER_CELL <- 200

# =====================================================
# Optimized Helper Functions
# =====================================================

ensure_qc_metrics <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"

  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  # ⭐ Calculate nCount/nFeature ONCE on original matrix (before gene filtering)
  if (!"nCount_RNA_raw" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nCount_RNA_raw <- Matrix::colSums(counts_mat)
    seurat_obj$nFeature_RNA_raw <- Matrix::colSums(counts_mat > 0)
    cat("  Created nCount_RNA_raw / nFeature_RNA_raw\n")
  }

  # Copy to standard names if missing
  if (!"nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nCount_RNA <- seurat_obj$nCount_RNA_raw
  }
  if (!"nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nFeature_RNA <- seurat_obj$nFeature_RNA_raw
  }

  # ⭐ Calculate percent.mt/rb on original matrix (before gene filtering)
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

# ⭐ Remove technical genes (MT/RP/pseudogenes/LINC)
filter_low_quality_genes_fast <- function(seurat_obj, min_cells = 3) {
  DefaultAssay(seurat_obj) <- "RNA"
  all_genes <- rownames(seurat_obj)

  # Combine patterns into single regex for efficiency
  # ⭐ Include MT/RP removal (standard practice for scRNA-seq)
  combined_pattern <- paste(
    c(
      "^MT-", # Mitochondrial genes
      "^RPS|^RPL|^MRPS|^MRPL", # Ribosomal proteins
      "^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$", # Ribosomal pseudogenes
      "^(AC|AL|AP|BX|Z)[0-9]+\\.", # Non-coding LINC
      "^RP[0-9]+-", # Readthrough transcripts
      "^CTD-|^CTB-|^CTC-", # Clone-based annotations
      "-OT[0-9]+$", # Overlapping transcripts
      "^LOC[0-9]+" # Generic loci
    ),
    collapse = "|"
  )

  genes_to_remove <- grep(combined_pattern, all_genes, value = TRUE)
  cat(sprintf(
    "  Genes to remove by pattern: %d\n",
    length(genes_to_remove)
  ))

  # Count by category for reporting
  mt_count <- length(grep("^MT-", genes_to_remove, value = TRUE))
  rp_count <- length(grep(
    "^RPS|^RPL|^MRPS|^MRPL",
    genes_to_remove,
    value = TRUE
  ))
  cat(sprintf(
    "    MT: %d, RP: %d, Other: %d\n",
    mt_count,
    rp_count,
    length(genes_to_remove) - mt_count - rp_count
  ))

  genes_to_keep <- setdiff(all_genes, genes_to_remove)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)

  # Min cells filtering
  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  gene_ncells <- Matrix::rowSums(counts_mat > 0)
  keep_genes <- names(gene_ncells[gene_ncells >= min_cells])

  before_n <- nrow(seurat_obj)
  seurat_obj <- subset(seurat_obj, features = keep_genes)
  after_n <- nrow(seurat_obj)

  cat(sprintf(
    "  Min cells filtering: %d → %d (removed %d)\n",
    before_n,
    after_n,
    before_n - after_n
  ))

  # ⭐ Update filtered QC metrics (but keep raw versions)
  counts_mat_filtered <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  seurat_obj$nCount_RNA_filtered <- Matrix::colSums(counts_mat_filtered)
  seurat_obj$nFeature_RNA_filtered <- Matrix::colSums(counts_mat_filtered > 0)

  cat(sprintf(
    "  QC after filtering: mean nCount_filtered=%.0f, mean nFeature_filtered=%.0f\n",
    mean(seurat_obj$nCount_RNA_filtered),
    mean(seurat_obj$nFeature_RNA_filtered)
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

  # ⭐ ScaleData only on HVG
  seurat_obj <- ScaleData(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    verbose = FALSE
  )

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

# ⭐ Use DietSeurat to reduce memory for parallel marker finding
find_all_markers_parallel_v2 <- function(
  seurat_obj,
  test_use = "wilcox",
  only_pos = TRUE,
  min_pct = 0.25,
  logfc_threshold = 0.25,
  latent_vars = NULL
) {
  cat("Preparing lean object for parallel marker finding...\n")

  # Create a diet version: only keep counts/data, remove scale.data and reductions
  seurat_lean <- DietSeurat(
    seurat_obj,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    assays = "RNA"
  )

  # Keep essential metadata
  essential_cols <- intersect(
    c(
      "nCount_RNA_raw",
      "nFeature_RNA_raw",
      "percent.mt",
      "percent.rb"
    ),
    colnames(seurat_obj@meta.data)
  )

  seurat_lean@meta.data <- seurat_obj@meta.data[, essential_cols, drop = FALSE]

  # Copy identities
  Idents(seurat_lean) <- Idents(seurat_obj)

  cat(sprintf(
    "  Lean object size: %.1f MB (vs %.1f MB original)\n",
    object.size(seurat_lean) / 1024^2,
    object.size(seurat_obj) / 1024^2
  ))

  clusters <- levels(Idents(seurat_lean))
  cat(sprintf(
    "Finding markers for %d clusters in parallel...\n",
    length(clusters)
  ))

  # Parallel execution
  marker_list <- future_lapply(
    clusters,
    function(cluster_id) {
      tryCatch(
        {
          if (test_use == "MAST" && !is.null(latent_vars)) {
            # Use raw QC metrics for MAST
            available_latent <- intersect(
              latent_vars,
              colnames(seurat_lean@meta.data)
            )
            FindMarkers(
              seurat_lean,
              ident.1 = cluster_id,
              only.pos = only_pos,
              min.pct = min_pct,
              logfc.threshold = logfc_threshold,
              test.use = "MAST",
              latent.vars = available_latent,
              verbose = FALSE
            )
          } else {
            FindMarkers(
              seurat_lean,
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

  # ⭐ Check for empty results
  if (length(marker_list) == 0) {
    warning("No markers found in any cluster!")
    return(data.frame())
  }

  markers_df <- bind_rows(lapply(names(marker_list), function(cid) {
    df <- marker_list[[cid]]
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$cluster <- cid
    df$gene <- rownames(df)
    df
  }))

  # Remove NULL entries
  markers_df <- markers_df[!is.na(markers_df$gene), ]

  return(markers_df)
}

# ⭐ Gene ID sanity check before enrichment
check_gene_id_consistency <- function(seurat_obj, msigdb_gmt) {
  cat("\n=== Gene ID Consistency Check ===\n")

  all_genes <- rownames(seurat_obj)

  # Check ENSG pattern
  ensg_pattern <- "^ENSG[0-9]+"
  ensg_genes <- grep(ensg_pattern, all_genes, value = TRUE)
  ensg_ratio <- length(ensg_genes) / length(all_genes)

  cat(sprintf(
    "  Total genes: %d\n",
    length(all_genes)
  ))
  cat(sprintf(
    "  ENSG genes: %d (%.1f%%)\n",
    length(ensg_genes),
    ensg_ratio * 100
  ))

  # Check MSigDB GMT format
  if (file.exists(msigdb_gmt)) {
    gmt_sample <- readLines(msigdb_gmt, n = 5)
    cat("  MSigDB GMT sample:\n")
    cat(paste0("    ", substr(gmt_sample[1], 1, 100), "...\n"))
  }

  # Warning if high ENSG ratio
  if (ensg_ratio > 0.1) {
    cat("\n  ⚠️ WARNING: >10% genes are ENSG IDs!\n")
    cat("  MSigDB uses HGNC symbols. Enrichment hit rate may be low.\n")
    cat("  Recommendation: Convert ENSG to SYMBOL before enrichment.\n")
  } else {
    cat("  ✓ Gene IDs look good (mostly SYMBOL format)\n")
  }

  # Sample common basal markers
  common_genes <- intersect(
    all_genes,
    c(
      "TP63",
      "KRT5",
      "KRT14", # Basal markers
      "MKI67",
      "PCNA",
      "TOP2A", # Proliferation
      "KRT13",
      "KRT4",
      "IVL", # Squamous differentiation
      "IL6",
      "CXCL8",
      "TNF" # Inflammation
    )
  )
  cat(sprintf(
    "  Common basal/inflammation markers found: %d/12\n",
    length(common_genes)
  ))
  if (length(common_genes) > 0) {
    cat(sprintf("    %s\n", paste(common_genes, collapse = ", ")))
  }

  cat("=================================\n\n")
}

run_enrichment_safe <- function(seurat_obj, cluster_col, msigdb_gmt) {
  cat("Running enrichment analysis (with error handling)...\n")

  # ⭐ Gene ID check before enrichment
  check_gene_id_consistency(seurat_obj, msigdb_gmt)

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

cat("=== PRODUCTION PIPELINE v1.3 START ===\n")
start_time <- Sys.time()

cat("\n[1/9] Reading data...\n")
t1 <- Sys.time()

file_ext <- tools::file_ext(H5AD_FILE)

if (file_ext == "h5ad") {
  cat("Loading h5ad file...\n")
  seurat_obj <- GetSeurat(h5ad_path = H5AD_FILE, debug = TRUE)
} else if (file_ext == "rds") {
  cat("Loading RDS file...\n")
  seurat_obj <- readRDS(H5AD_FILE)
} else {
  stop(sprintf(
    "Unsupported file format: %s. Only .h5ad and .rds are supported.",
    file_ext
  ))
}

cat(sprintf("Successfully loaded: %s\n", basename(H5AD_FILE)))
cat(sprintf("Cells: %d, Genes: %d\n", ncol(seurat_obj), nrow(seurat_obj)))
DefaultAssay(seurat_obj) <- "RNA"
cat(sprintf(
  "  Loaded: %d cells × %d genes (%.1f sec)\n",
  ncol(seurat_obj),
  nrow(seurat_obj),
  as.numeric(difftime(Sys.time(), t1, units = "secs"))
))

cat("\n[2/9] Ensuring QC metrics (on original matrix)...\n")
t2 <- Sys.time()
seurat_obj <- ensure_qc_metrics(seurat_obj)
cat(sprintf(
  "  QC metrics (raw): mean nCount=%.0f, mean nFeature=%.0f, mean percent.mt=%.2f%%\n",
  mean(seurat_obj$nCount_RNA_raw),
  mean(seurat_obj$nFeature_RNA_raw),
  mean(seurat_obj$percent.mt)
))
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t2, units = "secs"))
))

cat("\n[2.5/9] Pre-filtering low-quality cells...\n")
t25 <- Sys.time()
before_cells <- ncol(seurat_obj)

# ⭐ Use raw QC metrics for filtering
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA_raw >= MIN_GENES_PER_CELL &
    percent.mt < 20 &
    nCount_RNA_raw > 0
)

after_cells <- ncol(seurat_obj)
cat(sprintf(
  "  Cells: %d → %d (removed %d low-quality, %.1f sec)\n",
  before_cells,
  after_cells,
  before_cells - after_cells,
  as.numeric(difftime(Sys.time(), t25, units = "secs"))
))

cat("\n[3/9] Gene filtering (remove MT/RP/pseudogenes)...\n")
t3 <- Sys.time()
seurat_obj <- filter_low_quality_genes_fast(
  seurat_obj,
  min_cells = MIN_CELLS_PER_GENE
)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t3, units = "secs"))
))

cat("\n[4/9] Normalize data (sequential)...\n")
t4 <- Sys.time()
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  nfeatures = 2000,
  verbose = FALSE
)

cat(sprintf(
  "  Variable features: %d\n",
  length(VariableFeatures(seurat_obj))
))

# ⭐ ScaleData only on HVG
cat("  Scaling HVG only (memory efficient)...\n")
seurat_obj <- ScaleData(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  verbose = FALSE
)

cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t4, units = "secs"))
))

cat("\n[5/9] Determine cluster column...\n")
cluster_col <- pick_cluster_col(seurat_obj, CLUSTER_COL)
cat(sprintf("  Using: %s\n", cluster_col))
cluster_table <- sort(
  table(seurat_obj@meta.data[[cluster_col]]),
  decreasing = TRUE
)
print(cluster_table)

# Optional downsampling
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

# Parallel FindAllMarkers
if (isTRUE(RUN_FIND_ALL_MARKERS)) {
  cat("\n[6/9] FindAllMarkers (PARALLEL MODE)...\n")
  t6 <- Sys.time()

  # ⭐ Platform compatibility
  if (future::supportsMulticore()) {
    cat(sprintf("  Switching to multicore (%d cores)\n", N_CORES))
    plan("multicore", workers = N_CORES)
  } else {
    cat(sprintf("  Switching to multisession (%d workers)\n", N_CORES))
    plan("multisession", workers = N_CORES)
  }

  Idents(seurat_obj_subset) <- cluster_col

  if (DE_METHOD == "MAST" && requireNamespace("MAST", quietly = TRUE)) {
    cat("  Using MAST (accurate, uses raw QC covariates)...\n")
    markers <- find_all_markers_parallel_v2(
      seurat_obj_subset,
      test_use = "MAST",
      latent_vars = c("nCount_RNA_raw", "percent.mt")
    )
  } else {
    if (DE_METHOD == "MAST") {
      cat("  ⚠️ MAST not available; using Wilcoxon.\n")
    } else {
      cat(sprintf("  Using %s (fast)...\n", DE_METHOD))
    }
    markers <- find_all_markers_parallel_v2(
      seurat_obj_subset,
      test_use = DE_METHOD
    )
  }

  plan("sequential")
  cat("  ⭐ Switched back to sequential mode\n")

  t6_elapsed <- as.numeric(difftime(Sys.time(), t6, units = "secs"))

  # ⭐ Check for empty results
  if (is.data.frame(markers) && nrow(markers) > 0) {
    cat(sprintf(
      "  Found %d markers (%.1f sec, %.1f markers/sec)\n",
      nrow(markers),
      t6_elapsed,
      nrow(markers) / t6_elapsed
    ))

    write.csv(markers, "all_markers_complete.csv", row.names = FALSE)

    top_markers <- markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
    write.csv(top_markers, "all_markers_top10.csv", row.names = FALSE)

    cat("  ✓ Marker files saved\n")
  } else {
    cat("  ⚠️ No markers found (results empty)\n")
  }
}

# Enrichment analysis
cat("\n[7/9] Enrichment analysis...\n")
t7 <- Sys.time()
run_enrichment_safe(seurat_obj, cluster_col, MSIGDB_GMT)
cat(sprintf(
  "  Done (%.1f sec)\n",
  as.numeric(difftime(Sys.time(), t7, units = "secs"))
))

# Plots
cat("\n[8/9] Generating plots...\n")
t8 <- Sys.time()
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
    # ⭐ Use raw QC metrics for plotting
    qc_features <- intersect(
      c(
        "nCount_RNA_raw",
        "nFeature_RNA_raw",
        "nCount_RNA_filtered",
        "nFeature_RNA_filtered",
        "percent.mt",
        "percent.rb"
      ),
      colnames(seurat_obj@meta.data)
    )
    print(VlnPlot(
      seurat_obj,
      features = qc_features[1:min(6, length(qc_features))],
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
  as.numeric(difftime(Sys.time(), t8, units = "secs"))
))

cat("\n[9/9] Saving object...\n")
saveRDS(seurat_obj, "seurat_final_with_markers.rds")

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n=== PIPELINE COMPLETE ===\n")
cat(sprintf("Total time: %.1f minutes\n", total_time))
cat(sprintf("Output dir: %s\n", OUTPUT_DIR))
cat("\nKey improvements in v1.3:\n")
cat("  ✓ QC metrics: raw vs filtered separation\n")
cat("  ✓ MT/RP genes: removed from matrix (standard practice)\n")
cat("  ✓ Gene ID check before enrichment\n")
cat("  ✓ DietSeurat for memory-safe parallel markers\n")
cat("  ✓ Platform-aware parallel strategy\n")
