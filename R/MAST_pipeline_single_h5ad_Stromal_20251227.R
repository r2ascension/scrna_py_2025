#!/usr/bin/env Rscript
# ===== h5ad_mast_go_kegg_single_v1_0_PRODUCTION.R =====
# Purpose: Single h5ad → use existing cluster column → MAST/GO/GSEA/KEGG enrichment
# Author: r2end
# Date: 2024-12-24
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
})

# ===== Python / Conda Environment =====
use_condaenv("bbknn_env", required = TRUE)
py_config()

# ===== Load Custom Functions =====
# Must provide:
# - GetSeurat(h5ad_path, debug=TRUE/FALSE)
# - run_one_vs_rest_enrichment(seurat_obj, group_col, msigdb_gmt_file, run_go, run_gsea, ...)
source("/home/h2048/script/R/tissue_comparison_analysis_20251222.R")

# ===== Configuration =====
H5AD_FILE <- "/home/h2048/data/py/1217/cnmf_batch_production_v1_1_1/Stromal_Vascular/batch_aware/cnmf_analysis_k50_1/Stromal_Vascular_with_cnmf_k50.h5ad"

OUTPUT_DIR <- "/home/h2048/data/py/1217/cnmf_batch_production_v1_1_1/Stromal_Vascular/h5ad_single_mast_go_kegg"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTPUT_DIR)

MSIGDB_GMT <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# 你要使用的“分群列”（来自 h5ad obs → Seurat meta.data）
# 例如： "leiden" / "louvain" / "scanvi_predictions" / "manual_annotations" 等
CLUSTER_COL <- "leiden_bbknn_res1.0" # <<< 改成你h5ad里实际的cluster列名

# Marker / enrichment options
RUN_FIND_ALL_MARKERS <- TRUE
DE_METHOD <- "MAST" # "wilcox" or "MAST"（仅影响 FindAllMarkers；enrichment 仍由 run_one_vs_rest_enrichment 控制）

# Gene filtering
MIN_CELLS_PER_GENE <- 3

# =====================================================
# Helper Functions
# =====================================================

ensure_qc_metrics <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"

  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    cat("  Creating percent.mt...\n")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
      seurat_obj,
      pattern = "^MT-"
    )
  } else {
    mt_summary <- summary(seurat_obj$percent.mt)
    if (mt_summary["Max."] == 0 || all(is.na(seurat_obj$percent.mt))) {
      cat("  ⚠️ percent.mt invalid, recalculating...\n")
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
        seurat_obj,
        pattern = "^MT-"
      )
    } else {
      cat(sprintf("  ✓ percent.mt exists (mean=%.2f%%)\n", mt_summary["Mean"]))
    }
  }

  if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
    cat("  Creating percent.rb...\n")
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
      seurat_obj,
      pattern = "^RP[SL]"
    )
  } else {
    rb_summary <- summary(seurat_obj$percent.rb)
    if (rb_summary["Max."] == 0 || all(is.na(seurat_obj$percent.rb))) {
      cat("  ⚠️ percent.rb invalid, recalculating...\n")
      seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
        seurat_obj,
        pattern = "^RP[SL]"
      )
    } else {
      cat(sprintf("  ✓ percent.rb exists (mean=%.2f%%)\n", rb_summary["Mean"]))
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

  seurat_obj$nCount_RNA <- Matrix::colSums(counts_mat)
  seurat_obj$nFeature_RNA <- Matrix::colSums(counts_mat > 0)

  cat(sprintf(
    "  Updated QC: mean nCount=%.0f, mean nFeature=%.0f\n",
    mean(seurat_obj$nCount_RNA),
    mean(seurat_obj$nFeature_RNA)
  ))
  seurat_obj
}

filter_low_quality_genes <- function(seurat_obj, min_cells = 3) {
  DefaultAssay(seurat_obj) <- "RNA"
  all_genes <- rownames(seurat_obj)

  patterns_to_remove <- c(
    "^MT-",
    "^RPS|^RPL|^MRPS|^MRPL",
    "^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$",
    "^(AC|AL|AP|BX|Z)[0-9]+\\.",
    "^RP[0-9]+-",
    "^CTD-|^CTB-|^CTC-",
    "-OT[0-9]+$",
    "^LOC[0-9]+"
  )

  genes_to_remove <- unique(unlist(lapply(patterns_to_remove, function(p) {
    grep(p, all_genes, value = TRUE)
  })))
  cat(sprintf(
    "  Total genes to remove by pattern: %d\n",
    length(genes_to_remove)
  ))

  genes_to_keep <- setdiff(all_genes, genes_to_remove)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)

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
    "Cluster column '%s' not found. Example available columns:\n%s",
    cluster_col,
    paste(head(md, 60), collapse = ", ")
  ))
}

ensure_umap <- function(seurat_obj, reduction_name = "umap", seed = 42) {
  if (reduction_name %in% names(seurat_obj@reductions)) {
    return(seurat_obj)
  }

  cat(sprintf(
    "  ⚠️ Reduction '%s' not found. Computing PCA+UMAP fallback for plotting...\n",
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
# Main
# =====================================================

cat("Step 1: Reading h5ad...\n")
seurat_obj <- GetSeurat(h5ad_path = H5AD_FILE, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"
cat(sprintf(
  "  Loaded: %d cells × %d genes\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
))

cat("Step 2: Ensuring QC metrics...\n")
seurat_obj <- ensure_qc_metrics(seurat_obj)

cat("Step 3: Gene filtering...\n")
seurat_obj <- filter_low_quality_genes(
  seurat_obj,
  min_cells = MIN_CELLS_PER_GENE
)

cat("Step 3b: Recalculate basic QC after filtering...\n")
seurat_obj <- recalc_basic_qc(seurat_obj)

cat("Step 4: Normalize data (for markers/enrichment downstream)...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  nfeatures = 2000,
  verbose = FALSE
)
# ScaleData 不是必须，但很多下游可视化更稳
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

cat("Step 5: Determine cluster/group column from h5ad obs...\n")
cluster_col <- pick_cluster_col(seurat_obj, CLUSTER_COL)
cat(sprintf("  Using cluster column: %s\n", cluster_col))
cat("  Group sizes:\n")
print(sort(table(seurat_obj@meta.data[[cluster_col]]), decreasing = TRUE))

# ===== Optional: Markers (not required for enrichment if your function computes its own DE) =====
if (isTRUE(RUN_FIND_ALL_MARKERS)) {
  cat("Step 6: FindAllMarkers...\n")
  Idents(seurat_obj) <- cluster_col

  if (DE_METHOD == "MAST" && requireNamespace("MAST", quietly = TRUE)) {
    markers <- FindAllMarkers(
      seurat_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "MAST",
      latent.vars = intersect(
        c("nCount_RNA", "percent.mt"),
        colnames(seurat_obj@meta.data)
      ),
      verbose = FALSE
    )
  } else {
    if (DE_METHOD == "MAST") {
      cat("  ⚠️ MAST not available; falling back to Wilcoxon.\n")
    }
    markers <- FindAllMarkers(
      seurat_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      verbose = FALSE
    )
  }

  write.csv(markers, "all_markers_complete.csv", row.names = FALSE)

  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
  write.csv(top_markers, "all_markers_top10.csv", row.names = FALSE)
}

# ===== Enrichment (MAST/GO/GSEA/KEGG handled in your custom function) =====
cat("Step 7: Enrichment (MAST→GO/GSEA/KEGG)...\n")
run_enrichment_safe(seurat_obj, cluster_col, MSIGDB_GMT)

# ===== Basic plots (optional but useful) =====
cat("Step 8: Plots...\n")
seurat_obj <- ensure_umap(seurat_obj, reduction_name = "umap", seed = 42)

# Plot 1: UMAP overview
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
    cat("✓ UMAP overview plot saved\n")
  },
  error = function(e) {
    cat(sprintf("⚠️ UMAP overview plot failed: %s\n", e$message))
    tryCatch(dev.off(), error = function(e) NULL)
  }
)

# Plot 2: QC violin plot
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
    cat("✓ QC violin plot saved\n")
  },
  error = function(e) {
    cat(sprintf("⚠️ QC violin plot failed: %s\n", e$message))
    tryCatch(dev.off(), error = function(e) NULL)
  }
)

# ===== Save object =====
cat("Step 9: Save analyzed object...\n")
saveRDS(seurat_obj, "seurat_from_h5ad_for_enrichment.rds")

cat("\n=== DONE ===\n")
cat(sprintf("Output dir: %s\n", OUTPUT_DIR))
