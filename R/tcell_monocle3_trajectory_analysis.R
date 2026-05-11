#!/usr/bin/env Rscript
################################################################################
# T Cell Trajectory Analysis with Monocle3
################################################################################
#
# Objective: Analyze T cell differentiation trajectories using Monocle3
#
# Workflow:
# 1. Load scANVI-integrated data from h5ad
# 2. Convert to Monocle3 CDS format
# 3. Analyze CD4 trajectory (Naive -> CM -> EM, with Treg branch)
# 4. Analyze CD8 trajectory (Naive -> CM -> EM -> TEMRA)
# 5. Identify pseudotime-dependent genes
# 6. Analyze branch points
#
# Input: adata_tcell_scANVI_annotated.h5ad (with full gene set)
# Output: Trajectory plots, pseudotime analysis, branch genes
#
# Author: r2end
# Date: 2024-12-04
################################################################################

# ==============================================================================
# SETUP
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
})

cat("\n", rep("=", 80), "\n", sep = "")
cat("T CELL TRAJECTORY ANALYSIS WITH MONOCLE3\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Monocle3 version:", as.character(packageVersion("monocle3")), "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
# Load required libraries
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
# source("path/to/SCNT.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

cat(rep("=", 80), "\n")
cat("CONFIGURATION\n")
cat(rep("=", 80), "\n\n")

# Input/Output paths
INPUT_H5AD <- '/home/h2048/data/py/1204/Tcell_scANVI/adata_tcell_scANVI_annotated_fullRaw.h5ad'
OUTPUT_DIR <- "/home/h2048/data/py/1204/Tcell_monocle3"

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "cd4_trajectory"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "cd8_trajectory"), showWarnings = FALSE)

# Data keys
LABELS_KEY <- "Multinomial_Label"
BATCH_KEY <- "dataset"

# Cell type definitions
CD4_CELLTYPES <- c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg")
CD8_CELLTYPES <- c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA")

# Root cell types for trajectory
CD4_ROOT <- "CD4_Naive"
CD8_ROOT <- "CD8_Naive"

# Monocle3 parameters
N_DIM_PREPROCESS <- 50
N_NEIGHBORS_CLUSTER <- 15

# Visualization
FIGURE_DPI <- 300
FIGURE_WIDTH <- 10
FIGURE_HEIGHT <- 8

cat("Configuration loaded:\n")
cat("  Input:", INPUT_H5AD, "\n")
cat("  Output:", OUTPUT_DIR, "\n")
cat("  CD4 trajectory:", paste(CD4_CELLTYPES, collapse = ", "), "\n")
cat("  CD8 trajectory:", paste(CD8_CELLTYPES, collapse = ", "), "\n\n")

# ==============================================================================
# STEP 1: LOAD DATA FROM H5AD
# ==============================================================================

# cat(rep("=", 80), "\n")
# cat("STEP 1: Loading Data from H5AD\n")
# cat(rep("=", 80), "\n\n")

# cat("Converting h5ad to h5seurat...\n")

# # Convert h5ad to h5seurat
# h5seurat_file <- file.path(OUTPUT_DIR, "temp_adata.h5seurat")

# if (!file.exists(h5seurat_file)) {
#   Convert(INPUT_H5AD, dest = "h5seurat", overwrite = TRUE,
#           assay = "RNA")
#   file.rename("temp_adata.h5seurat", h5seurat_file)
#   cat("  Conversion complete\n")
# } else {
#   cat("  Using existing h5seurat file\n")
# }

# # Load as Seurat object
# cat("\nLoading Seurat object...\n")
# seurat_obj <- LoadH5Seurat(h5seurat_file)
seurat_obj <- SCNT::GetSeurat(h5ad_path = INPUT_H5AD)
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")
cat("  Cell types:", length(unique(seurat_obj@meta.data[[LABELS_KEY]])), "\n\n")

# Display cell type distribution
cat("Cell type distribution:\n")
cell_counts <- table(seurat_obj@meta.data[[LABELS_KEY]])
for (ct in names(cell_counts)) {
  pct <- cell_counts[ct] / ncol(seurat_obj) * 100
  cat(sprintf("  %s: %d cells (%.1f%%)\n", ct, cell_counts[ct], pct))
}

# ==============================================================================
# HELPER FUNCTION: CONVERT SEURAT TO CDS
# ==============================================================================

seurat_to_cds <- function(seurat_obj, assay = "RNA") {
  cat("\nConverting Seurat to CDS format...\n")

  # Extract counts matrix
  counts <- GetAssayData(seurat_obj, slot = "counts", assay = assay)

  # Create cell metadata
  cell_metadata <- seurat_obj@meta.data

  # Create gene metadata
  gene_metadata <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts)
  )

  # Create CDS object
  cds <- new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Transfer UMAP if available
  if ("X_umap_scANVI" %in% names(seurat_obj@reductions)) {
    cat("  Transferring scANVI UMAP coordinates...\n")
    umap_coords <- Embeddings(seurat_obj, reduction = "X_umap_scANVI")
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

    # Add to CDS reducedDims
    reducedDims(cds)[["UMAP"]] <- umap_coords
  }

  # Transfer scANVI latent if available
  if ("X_scANVI" %in% names(seurat_obj@reductions)) {
    cat("  Transferring scANVI latent space...\n")
    latent <- Embeddings(seurat_obj, reduction = "X_scANVI")
    reducedDims(cds)[["scANVI"]] <- latent
  }

  cat("  CDS created successfully\n")
  return(cds)
}

# ==============================================================================
# HELPER FUNCTION: PROGRAMMATIC ROOT SELECTION
# ==============================================================================

get_root_principal_node <- function(cds, cell_type_key, root_cell_type) {
  cat("\nSelecting root node programmatically...\n")
  cat("  Root cell type:", root_cell_type, "\n")

  # Get cells of root cell type
  root_cells <- colData(cds)[[cell_type_key]] == root_cell_type
  cell_ids <- which(root_cells)

  cat("  Root cells:", length(cell_ids), "\n")

  if (length(cell_ids) == 0) {
    stop("No cells found for root cell type: ", root_cell_type)
  }

  # Get closest vertex for each cell
  closest_vertex <- cds@principal_graph_aux[[
    "UMAP"
  ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])

  # Find node with most root cells
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]

  cat("  Selected root node:", root_pr_nodes, "\n")

  return(root_pr_nodes)
}

# ==============================================================================
# HELPER FUNCTION: COMPLETE MONOCLE3 TRAJECTORY ANALYSIS
# ==============================================================================

run_monocle3_trajectory <- function(
  seurat_obj,
  cell_types,
  root_cell_type,
  output_prefix,
  fig_dir
) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("MONOCLE3 TRAJECTORY ANALYSIS:", output_prefix, "\n")
  cat(rep("=", 80), "\n\n")

  # ========== 1. Subset cells ==========
  cat("Step 1: Subsetting cells...\n")

  cells_to_keep <- seurat_obj@meta.data[[LABELS_KEY]] %in% cell_types
  seurat_subset <- subset(
    seurat_obj,
    cells = colnames(seurat_obj)[cells_to_keep]
  )

  cat("  Retained:", ncol(seurat_subset), "cells\n")
  for (ct in cell_types) {
    n <- sum(seurat_subset@meta.data[[LABELS_KEY]] == ct)
    cat(sprintf("    %s: %d cells\n", ct, n))
  }

  # ========== 2. Convert to CDS ==========
  cat("\nStep 2: Converting to CDS...\n")
  cds <- seurat_to_cds(seurat_subset)

  # ========== 3. Preprocess ==========
  cat("\nStep 3: Preprocessing...\n")

  # Use existing scANVI latent space if available
  if ("scANVI" %in% names(reducedDims(cds))) {
    cat("  Using existing scANVI latent space\n")
    # Monocle3 expects PCA slot for preprocessing
    reducedDims(cds)[["PCA"]] <- reducedDims(cds)[["scANVI"]]
  } else {
    cat("  Running PCA...\n")
    cds <- preprocess_cds(cds, num_dim = N_DIM_PREPROCESS)
  }

  # ========== 4. Reduce dimension ==========
  cat("\nStep 4: Dimensionality reduction...\n")

  if ("UMAP" %in% names(reducedDims(cds))) {
    cat("  Using existing scANVI UMAP\n")
  } else {
    cat("  Computing UMAP...\n")
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
  }

  # ========== 5. Cluster cells ==========
  cat("\nStep 5: Clustering cells...\n")
  cds <- cluster_cells(cds, k = N_NEIGHBORS_CLUSTER)

  n_clusters <- length(unique(clusters(cds)))
  n_partitions <- length(unique(partitions(cds)))
  cat("  Clusters:", n_clusters, "\n")
  cat("  Partitions:", n_partitions, "\n")

  # ========== 6. Learn trajectory graph ==========
  cat("\nStep 6: Learning trajectory graph...\n")
  cds <- learn_graph(cds, use_partition = FALSE)
  cat("  Trajectory graph learned\n")

  # ========== 7. Order cells in pseudotime ==========
  cat("\nStep 7: Ordering cells in pseudotime...\n")

  # Get root node programmatically
  root_node <- get_root_principal_node(cds, LABELS_KEY, root_cell_type)

  # Order cells
  cds <- order_cells(cds, root_pr_nodes = root_node)

  # Check pseudotime distribution
  pseudotime_vals <- pseudotime(cds)
  finite_pseudotime <- pseudotime_vals[is.finite(pseudotime_vals)]

  cat("  Pseudotime statistics:\n")
  cat(sprintf(
    "    Cells with finite pseudotime: %d (%.1f%%)\n",
    length(finite_pseudotime),
    length(finite_pseudotime) / length(pseudotime_vals) * 100
  ))
  cat(sprintf(
    "    Pseudotime range: %.2f - %.2f\n",
    min(finite_pseudotime),
    max(finite_pseudotime)
  ))

  # ========== 8. Visualization ==========
  cat("\nStep 8: Generating visualizations...\n")

  # Plot 1: Cell types on trajectory
  p1 <- plot_cells(
    cds,
    color_cells_by = LABELS_KEY,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = TRUE,
    graph_label_size = 3
  ) +
    ggtitle(paste(output_prefix, "- Cell Types on Trajectory")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_trajectory_celltypes.pdf")),
    p1,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    dpi = FIGURE_DPI
  )
  cat("  Saved: trajectory_celltypes.pdf\n")

  # Plot 2: Pseudotime
  p2 <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = TRUE,
    graph_label_size = 3
  ) +
    scale_color_viridis(option = "viridis", na.value = "grey80") +
    ggtitle(paste(output_prefix, "- Pseudotime")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_trajectory_pseudotime.pdf")),
    p2,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    dpi = FIGURE_DPI
  )
  cat("  Saved: trajectory_pseudotime.pdf\n")

  # Plot 3: Partitions
  p3 <- plot_cells(
    cds,
    color_cells_by = "partition",
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3
  ) +
    ggtitle(paste(output_prefix, "- Partitions and Branch Points")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_trajectory_partitions.pdf")),
    p3,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    dpi = FIGURE_DPI
  )
  cat("  Saved: trajectory_partitions.pdf\n")

  # ========== 9. Export pseudotime ==========
  cat("\nStep 9: Exporting pseudotime values...\n")

  pseudotime_df <- data.frame(
    cell_id = colnames(cds),
    cell_type = colData(cds)[[LABELS_KEY]],
    pseudotime = pseudotime(cds),
    partition = partitions(cds),
    cluster = clusters(cds)
  )

  output_file <- file.path(fig_dir, paste0(output_prefix, "_pseudotime.csv"))
  write.csv(pseudotime_df, output_file, row.names = FALSE)
  cat("  Saved: pseudotime.csv\n")

  # ========== 10. Summary statistics ==========
  cat("\nSummary Statistics:\n")
  cat(rep("-", 80), "\n", sep = "")

  for (ct in cell_types) {
    ct_pseudotime <- pseudotime_df$pseudotime[pseudotime_df$cell_type == ct]
    ct_pseudotime <- ct_pseudotime[is.finite(ct_pseudotime)]

    if (length(ct_pseudotime) > 0) {
      cat(sprintf("  %s:\n", ct))
      cat(sprintf("    Cells: %d\n", length(ct_pseudotime)))
      cat(sprintf(
        "    Pseudotime range: %.2f - %.2f\n",
        min(ct_pseudotime),
        max(ct_pseudotime)
      ))
      cat(sprintf("    Mean pseudotime: %.2f\n", mean(ct_pseudotime)))
    }
  }

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Trajectory analysis complete:", output_prefix, "\n")
  cat(rep("=", 80), "\n\n", sep = "")

  return(cds)
}

# ==============================================================================
# STEP 2: ANALYZE CD4 TRAJECTORY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("STEP 2: CD4 T Cell Trajectory Analysis\n")
cat(rep("=", 80), "\n\n")

cds_cd4 <- run_monocle3_trajectory(
  seurat_obj = seurat_obj,
  cell_types = CD4_CELLTYPES,
  root_cell_type = CD4_ROOT,
  output_prefix = "CD4",
  fig_dir = file.path(OUTPUT_DIR, "cd4_trajectory")
)

# ==============================================================================
# STEP 3: ANALYZE CD8 TRAJECTORY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("STEP 3: CD8 T Cell Trajectory Analysis\n")
cat(rep("=", 80), "\n\n")

cds_cd8 <- run_monocle3_trajectory(
  seurat_obj = seurat_obj,
  cell_types = CD8_CELLTYPES,
  root_cell_type = CD8_ROOT,
  output_prefix = "CD8",
  fig_dir = file.path(OUTPUT_DIR, "cd8_trajectory")
)

# ==============================================================================
# STEP 4: IDENTIFY PSEUDOTIME-DEPENDENT GENES
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("STEP 4: Identifying Pseudotime-Dependent Genes\n")
cat(rep("=", 80), "\n\n")

# Function to find pseudotime-dependent genes
find_pseudotime_genes <- function(
  cds,
  output_prefix,
  fig_dir,
  num_genes = 100
) {
  cat("Finding genes that change along pseudotime...\n")
  cat("  (This may take 5-10 minutes)\n\n")

  # Test genes for differential expression along pseudotime
  # Using graph_test which is optimized for trajectory analysis
  tryCatch(
    {
      gene_fits <- graph_test(
        cds,
        neighbor_graph = "principal_graph",
        cores = 4
      )
      gene_fits <- gene_fits[order(gene_fits$q_value), ]

      # Filter significant genes
      sig_genes <- gene_fits[gene_fits$q_value < 0.05, ]

      cat("  Significant genes (q < 0.05):", nrow(sig_genes), "\n")

      # Save results
      output_file <- file.path(
        fig_dir,
        paste0(output_prefix, "_pseudotime_genes.csv")
      )
      write.csv(gene_fits, output_file, row.names = FALSE)
      cat("  Saved:", paste0(output_prefix, "_pseudotime_genes.csv\n"))

      # Plot top genes
      if (nrow(sig_genes) > 0) {
        top_genes <- head(rownames(sig_genes), num_genes)

        cat("\n  Top 10 pseudotime-dependent genes:\n")
        for (i in 1:min(10, length(top_genes))) {
          gene <- top_genes[i]
          q_val <- sig_genes[gene, "q_value"]
          cat(sprintf("    %2d. %s (q = %.2e)\n", i, gene, q_val))
        }

        # Visualize top genes
        cat("\n  Plotting top genes...\n")

        # Select representative genes for visualization
        genes_to_plot <- head(top_genes, 12)

        p <- plot_cells(
          cds,
          genes = genes_to_plot,
          show_trajectory_graph = FALSE,
          label_cell_groups = FALSE,
          label_leaves = FALSE
        )

        ggsave(
          file.path(fig_dir, paste0(output_prefix, "_top_genes.pdf")),
          p,
          width = 14,
          height = 10,
          dpi = FIGURE_DPI
        )
        cat("  Saved: top_genes.pdf\n")
      }

      return(gene_fits)
    },
    error = function(e) {
      cat("  Error in graph_test:", conditionMessage(e), "\n")
      cat("  Skipping pseudotime gene analysis\n")
      return(NULL)
    }
  )
}

# Run for CD4
cat("\nAnalyzing CD4 pseudotime-dependent genes...\n")
cd4_genes <- find_pseudotime_genes(
  cds_cd4,
  output_prefix = "CD4",
  fig_dir = file.path(OUTPUT_DIR, "cd4_trajectory")
)

# Run for CD8
cat("\nAnalyzing CD8 pseudotime-dependent genes...\n")
cd8_genes <- find_pseudotime_genes(
  cds_cd8,
  output_prefix = "CD8",
  fig_dir = file.path(OUTPUT_DIR, "cd8_trajectory")
)

# ==============================================================================
# STEP 5: ANALYZE KEY MARKER GENES ALONG TRAJECTORY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("STEP 5: Analyzing Key Marker Genes\n")
cat(rep("=", 80), "\n\n")

# Define key marker genes for T cell differentiation
cd4_markers <- c("CCR7", "SELL", "IL7R", "GZMK", "FOXP3", "IL2RA")
cd8_markers <- c("CCR7", "SELL", "IL7R", "GZMK", "GZMB", "PRF1")

# Function to plot marker genes
plot_marker_genes <- function(cds, markers, output_prefix, fig_dir) {
  cat("Plotting marker genes for", output_prefix, "...\n")

  # Check which markers are available
  available_markers <- markers[markers %in% rownames(cds)]

  if (length(available_markers) == 0) {
    cat("  No marker genes found in dataset\n")
    return(NULL)
  }

  cat("  Available markers:", paste(available_markers, collapse = ", "), "\n")

  # Plot genes on trajectory
  p <- plot_cells(
    cds,
    genes = available_markers,
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_marker_genes.pdf")),
    p,
    width = 12,
    height = 10,
    dpi = FIGURE_DPI
  )
  cat("  Saved: marker_genes.pdf\n\n")
}

# Plot CD4 markers
plot_marker_genes(
  cds_cd4,
  cd4_markers,
  "CD4",
  file.path(OUTPUT_DIR, "cd4_trajectory")
)

# Plot CD8 markers
plot_marker_genes(
  cds_cd8,
  cd8_markers,
  "CD8",
  file.path(OUTPUT_DIR, "cd8_trajectory")
)

# ==============================================================================
# STEP 6: SAVE CDS OBJECTS
# ==============================================================================

cat(rep("=", 80), "\n")
cat("STEP 6: Saving CDS Objects\n")
cat(rep("=", 80), "\n\n")

# Save CD4 CDS
cd4_file <- file.path(OUTPUT_DIR, "cd4_trajectory", "cds_cd4.rds")
saveRDS(cds_cd4, cd4_file)
cat("Saved CD4 CDS:", cd4_file, "\n")

# Save CD8 CDS
cd8_file <- file.path(OUTPUT_DIR, "cd8_trajectory", "cds_cd8.rds")
saveRDS(cds_cd8, cd8_file)
cat("Saved CD8 CDS:", cd8_file, "\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat(rep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n\n")

cat("Summary:\n")
cat(rep("-", 80), "\n", sep = "")

cat("\nCD4 Trajectory:\n")
cat("  Cell types:", paste(CD4_CELLTYPES, collapse = ", "), "\n")
cat("  Root:", CD4_ROOT, "\n")
cat("  Cells:", ncol(cds_cd4), "\n")
cat("  Output:", file.path(OUTPUT_DIR, "cd4_trajectory"), "\n")

cat("\nCD8 Trajectory:\n")
cat("  Cell types:", paste(CD8_CELLTYPES, collapse = ", "), "\n")
cat("  Root:", CD8_ROOT, "\n")
cat("  Cells:", ncol(cds_cd8), "\n")
cat("  Output:", file.path(OUTPUT_DIR, "cd8_trajectory"), "\n")

cat("\nGenerated files:\n")
cat("  - trajectory_celltypes.pdf\n")
cat("  - trajectory_pseudotime.pdf\n")
cat("  - trajectory_partitions.pdf\n")
cat("  - pseudotime.csv\n")
cat("  - pseudotime_genes.csv\n")
cat("  - top_genes.pdf\n")
cat("  - marker_genes.pdf\n")
cat("  - cds.rds (saved CDS objects)\n")

cat("\nNext steps:\n")
cat("  1. Review trajectory plots to verify biological plausibility\n")
cat("  2. Examine pseudotime_genes.csv for differentially expressed genes\n")
cat("  3. Validate key marker genes along pseudotime\n")
cat("  4. Analyze branch points (if present) for fate decisions\n")
cat("  5. Compare CD4 vs CD8 differentiation dynamics\n")

cat("\n", rep("=", 80), "\n", sep = "")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 80), "\n\n")

cat("Happy analyzing! 🔬\n\n")

################################################################################
# END OF SCRIPT
################################################################################
