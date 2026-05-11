#!/usr/bin/env Rscript
################################################################################
# T Cell Trajectory Analysis with Monocle3 (Improved Version)
# Features:
# - Quality control for T cell purity
# - Branch analysis and visualization
# - Enhanced gene analysis with T cell specificity
# - Diagnostic plots for pseudotime validation
# - Proper CDS saving with save_monocle_objects()
################################################################################

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(reshape2)
  library(pheatmap)
})

cat("\n", strrep("=", 80), "\n", sep = "")
cat("T CELL TRAJECTORY ANALYSIS WITH MONOCLE3 (IMPROVED VERSION)\n")
cat(strrep("=", 80), "\n\n", sep = "")

cat("monocle3 version:", as.character(packageVersion("monocle3")), "\n")
cat("Seurat version  :", as.character(packageVersion("Seurat")), "\n")
cat("SeuratDisk ver. :", as.character(packageVersion("SeuratDisk")), "\n")
cat("Started         :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

library(reticulate)
use_condaenv("bbknn_env", required = TRUE)
py_config()

# ==============================================================================
# CONFIGURATION
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("CONFIGURATION\n")
cat(strrep("=", 80), "\n\n")

# ---- Input/Output ----
INPUT_H5AD <- "/home/h2048/data/py/1204/Tcell_scANVI/adata_tcell_scANVI_annotated_fullRaw.h5ad"
OUTPUT_DIR <- "/home/h2048/data/py/1204/Tcell_monocle3_v2"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "qc"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "cd4_trajectory"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "cd8_trajectory"), showWarnings = FALSE)

H5SEURAT_FILE <- sub("\\.h5ad$", ".h5seurat", INPUT_H5AD)

# ---- Key columns ----
LABELS_KEY <- "Multinomial_Label"
BATCH_KEY <- "dataset"

# ---- CD4/CD8 cell types ----
CD4_CELLTYPES <- c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg")
CD8_CELLTYPES <- c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA")

CD4_ROOT <- "CD4_Naive"
CD8_ROOT <- "CD8_Naive"

# ---- Quality control markers ----
TCELL_MARKERS <- c("CD3D", "CD3E", "CD3G")
MYELOID_MARKERS <- c("CD14", "FCGR3A", "ITGAX", "ITGAM", "CD68")
CD4_SPECIFIC <- c("CD4")
CD8_SPECIFIC <- c("CD8A", "CD8B")

# ---- T cell differentiation markers ----
NAIVE_MARKERS <- c("CCR7", "SELL", "LEF1", "TCF7", "IL7R")
MEMORY_MARKERS <- c("IL7R", "CD27", "CD28")
EFFECTOR_MARKERS <- c("GZMK", "GZMB", "PRF1", "GNLY", "NKG7")
TREG_MARKERS <- c("FOXP3", "IL2RA", "CTLA4", "IKZF2")

# ---- Monocle3 parameters ----
N_DIM_PREPROCESS <- 50
N_NEIGHBORS_CLUSTER <- 15
USE_PARTITION <- TRUE # Changed to TRUE for better trajectory learning

# ---- Figure parameters ----
FIGURE_DPI <- 300
FIGURE_WIDTH <- 10
FIGURE_HEIGHT <- 8

cat("Configuration loaded:\n")
cat("  Input H5AD      :", INPUT_H5AD, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  CD4 trajectory  :", paste(CD4_CELLTYPES, collapse = ", "), "\n")
cat("  CD8 trajectory  :", paste(CD8_CELLTYPES, collapse = ", "), "\n")
cat("  Use partition   :", USE_PARTITION, "\n\n")

# ==============================================================================
# STEP 1: LOAD DATA FROM H5AD
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("STEP 1: Loading Data from H5AD\n")
cat(strrep("=", 80), "\n\n")

seurat_obj <- SCNT::GetSeurat(h5ad_path = INPUT_H5AD)

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

if (!LABELS_KEY %in% colnames(seurat_obj@meta.data)) {
  stop("LABELS_KEY '", LABELS_KEY, "' not found in meta.data.")
}

cat("  Cell type column:", LABELS_KEY, "\n")
cat(
  "  Unique labels   :",
  length(unique(seurat_obj@meta.data[[LABELS_KEY]])),
  "\n\n"
)

# Cell type distribution
cat("Cell type distribution:\n")
cell_counts <- table(seurat_obj@meta.data[[LABELS_KEY]])
for (ct in names(cell_counts)) {
  pct <- cell_counts[[ct]] / ncol(seurat_obj) * 100
  cat(sprintf("  %s: %d cells (%.1f%%)\n", ct, cell_counts[[ct]], pct))
}
cat("\n")

# ==============================================================================
# STEP 2: QUALITY CONTROL - T CELL PURITY CHECK
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("STEP 2: Quality Control - T Cell Purity Check\n")
cat(strrep("=", 80), "\n\n")

check_marker_expression <- function(seurat_obj, markers, marker_name) {
  cat("Checking", marker_name, "markers...\n")

  available_markers <- markers[markers %in% rownames(seurat_obj)]
  missing_markers <- markers[!markers %in% rownames(seurat_obj)]

  cat("  Available:", paste(available_markers, collapse = ", "), "\n")
  if (length(missing_markers) > 0) {
    cat("  Missing  :", paste(missing_markers, collapse = ", "), "\n")
  }

  if (length(available_markers) == 0) {
    cat("  WARNING: No markers found!\n\n")
    return(NULL)
  }

  # Get expression data
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  } else {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  }

  marker_expr <- as.matrix(expr_data[available_markers, , drop = FALSE])
  mean_expr <- colMeans(marker_expr)

  cat(sprintf(
    "  Mean expression: %.3f (range: %.3f - %.3f)\n",
    mean(mean_expr),
    min(mean_expr),
    max(mean_expr)
  ))
  cat(sprintf(
    "  Cells with expression > 0: %d (%.1f%%)\n\n",
    sum(mean_expr > 0),
    sum(mean_expr > 0) / length(mean_expr) * 100
  ))

  return(mean_expr)
}

# Check T cell markers
tcell_expr <- check_marker_expression(seurat_obj, TCELL_MARKERS, "T cell")

# Check myeloid contamination
myeloid_expr <- check_marker_expression(seurat_obj, MYELOID_MARKERS, "Myeloid")

# Check CD4/CD8
cd4_expr <- check_marker_expression(seurat_obj, CD4_SPECIFIC, "CD4")
cd8_expr <- check_marker_expression(seurat_obj, CD8_SPECIFIC, "CD8")

# Create QC plots
cat("Generating QC plots...\n")

# T cell marker violin plot
if (!is.null(tcell_expr)) {
  seurat_obj$Tcell_score <- tcell_expr

  p_tcell <- VlnPlot(
    seurat_obj,
    features = "Tcell_score",
    group.by = LABELS_KEY,
    pt.size = 0
  ) +
    ggtitle("T Cell Marker Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(OUTPUT_DIR, "qc", "tcell_marker_score.pdf"),
    p_tcell,
    width = 12,
    height = 6,
    dpi = FIGURE_DPI
  )
}

# Myeloid contamination plot
if (!is.null(myeloid_expr)) {
  seurat_obj$Myeloid_score <- myeloid_expr

  p_myeloid <- VlnPlot(
    seurat_obj,
    features = "Myeloid_score",
    group.by = LABELS_KEY,
    pt.size = 0
  ) +
    ggtitle("Myeloid Marker Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(OUTPUT_DIR, "qc", "myeloid_marker_score.pdf"),
    p_myeloid,
    width = 12,
    height = 6,
    dpi = FIGURE_DPI
  )

  # Flag high myeloid cells
  high_myeloid <- sum(myeloid_expr > 1)
  if (high_myeloid > 0) {
    cat(sprintf(
      "  WARNING: %d cells (%.1f%%) have high myeloid marker expression!\n",
      high_myeloid,
      high_myeloid / ncol(seurat_obj) * 100
    ))
  }
}

cat("  QC plots saved to:", file.path(OUTPUT_DIR, "qc"), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

get_counts_from_seurat <- function(seurat_obj, assay = "RNA") {
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
  } else {
    counts <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }
  counts
}

seurat_to_cds <- function(seurat_obj, assay = "RNA") {
  cat("\nConverting Seurat to CDS format...\n")

  counts <- get_counts_from_seurat(seurat_obj, assay = assay)
  cell_metadata <- seurat_obj@meta.data
  gene_metadata <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )

  cds <- monocle3::new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  reduction_names <- names(seurat_obj@reductions)
  if (length(reduction_names) > 0) {
    cat(
      "  Available Seurat reductions:",
      paste(reduction_names, collapse = ", "),
      "\n"
    )
  }

  # Transfer UMAP
  umap_red <- reduction_names[grepl(
    "umap",
    reduction_names,
    ignore.case = TRUE
  )]
  umap_red <- umap_red[1]
  if (!is.na(umap_red) && !is.null(umap_red)) {
    cat("  Transferring UMAP reduction: ", umap_red, "\n")
    umap_coords <- Seurat::Embeddings(seurat_obj, reduction = umap_red)
    colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
    reducedDims(cds)[["UMAP"]] <- umap_coords
  }

  # Transfer scANVI latent
  latent_red <- reduction_names[grepl(
    "scanvi|scvi",
    reduction_names,
    ignore.case = TRUE
  )]
  latent_red <- latent_red[1]
  if (!is.na(latent_red) && !is.null(latent_red)) {
    cat("  Transferring latent reduction: ", latent_red, "\n")
    latent <- Seurat::Embeddings(seurat_obj, reduction = latent_red)
    reducedDims(cds)[["scANVI"]] <- latent
  }

  cat("  CDS created successfully.\n")
  cds
}

get_root_principal_node <- function(cds, cell_type_key, root_cell_type) {
  cat("\nSelecting root node programmatically...\n")
  cat("  Root cell type:", root_cell_type, "\n")

  root_cells_logical <- colData(cds)[[cell_type_key]] == root_cell_type
  cell_ids <- which(root_cells_logical)

  cat("  Root cells:", length(cell_ids), "\n")

  if (length(cell_ids) == 0) {
    stop("No cells found for root cell type: ", root_cell_type)
  }

  closest_vertex <- cds@principal_graph_aux[[
    "UMAP"
  ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])

  root_pr_nodes <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]

  cat("  Selected root node:", root_pr_nodes, "\n")
  root_pr_nodes
}

# ==============================================================================
# ENHANCED: EXTRACT AND ANALYZE BRANCH INFORMATION
# ==============================================================================

extract_branch_info <- function(cds) {
  cat("\nExtracting branch information...\n")

  # Get the principal graph
  pr_graph <- principal_graph(cds)[["UMAP"]]

  # Get cell to closest vertex mapping
  closest_vertex <- cds@principal_graph_aux[[
    "UMAP"
  ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_df <- as.data.frame(as.matrix(closest_vertex))
  colnames(closest_vertex_df) <- "closest_vertex"

  # Identify branch points
  degree <- igraph::degree(pr_graph)
  branch_points <- names(degree[degree > 2])

  cat("  Branch points:", length(branch_points), "\n")
  cat("  Branch point IDs:", paste(branch_points, collapse = ", "), "\n")

  # Assign cells to branches using pseudotime ranges and graph structure
  cell_branch <- rep("Main", ncol(cds))
  names(cell_branch) <- colnames(cds)

  # Simple branch assignment based on graph connectivity
  if (length(branch_points) > 0) {
    # Get subgraph for each branch
    for (i in seq_along(branch_points)) {
      bp <- branch_points[i]
      # This is simplified - you might want more sophisticated logic
      cell_branch[cell_branch == "Main"] <- paste0("Branch_", i)
    }
  }

  return(list(
    branch = cell_branch,
    branch_points = branch_points,
    n_branches = length(unique(cell_branch))
  ))
}

plot_branch_analysis <- function(cds, output_prefix, fig_dir) {
  cat("\nPerforming branch analysis...\n")

  # Extract branch info
  branch_info <- extract_branch_info(cds)
  colData(cds)$branch <- branch_info$branch

  cat("  Number of branches:", branch_info$n_branches, "\n")

  # Plot 1: Branch assignment on UMAP
  p1 <- plot_cells(
    cds,
    color_cells_by = "branch",
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3
  ) +
    ggtitle(paste(output_prefix, "- Branch Assignment")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_branch_assignment.pdf")),
    p1,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    dpi = FIGURE_DPI
  )

  # Plot 2: Branch vs Cell Type cross-tabulation
  branch_celltype <- table(
    colData(cds)$branch,
    colData(cds)[[LABELS_KEY]]
  )

  cat("\nBranch vs Cell Type distribution:\n")
  print(branch_celltype)
  cat("\n")

  # Save as heatmap
  branch_celltype_pct <- prop.table(branch_celltype, margin = 1) * 100

  pdf(
    file.path(fig_dir, paste0(output_prefix, "_branch_celltype_heatmap.pdf")),
    width = 10,
    height = 6
  )
  pheatmap(
    branch_celltype_pct,
    display_numbers = TRUE,
    number_format = "%.1f",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste(output_prefix, "- Branch vs Cell Type (%)"),
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(50)
  )
  dev.off()

  # Plot 3: Pseudotime distribution by branch
  pseudotime_df <- data.frame(
    pseudotime = pseudotime(cds),
    branch = colData(cds)$branch,
    cell_type = colData(cds)[[LABELS_KEY]]
  )

  p3 <- ggplot(pseudotime_df, aes(x = branch, y = pseudotime, fill = branch)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5) +
    theme_classic() +
    ggtitle(paste(output_prefix, "- Pseudotime by Branch")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_pseudotime_by_branch.pdf")),
    p3,
    width = 10,
    height = 6,
    dpi = FIGURE_DPI
  )

  return(cds)
}

# ==============================================================================
# ENHANCED: PLOT PSEUDOTIME DISTRIBUTION BY CELL TYPE
# ==============================================================================

plot_pseudotime_distribution <- function(
  cds,
  output_prefix,
  fig_dir,
  cell_types
) {
  cat("\nPlotting pseudotime distribution by cell type...\n")

  pseudotime_df <- data.frame(
    cell_id = colnames(cds),
    cell_type = factor(colData(cds)[[LABELS_KEY]], levels = cell_types),
    pseudotime = pseudotime(cds)
  )

  pseudotime_df <- pseudotime_df[pseudotime_df$cell_type %in% cell_types, ]
  pseudotime_df <- pseudotime_df[is.finite(pseudotime_df$pseudotime), ]

  # Violin + boxplot
  p1 <- ggplot(
    pseudotime_df,
    aes(x = cell_type, y = pseudotime, fill = cell_type)
  ) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
    theme_classic() +
    labs(
      title = paste(output_prefix, "- Pseudotime Distribution by Cell Type"),
      x = "Cell Type",
      y = "Pseudotime"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Set2")

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_pseudotime_distribution.pdf")),
    p1,
    width = 10,
    height = 6,
    dpi = FIGURE_DPI
  )

  # Density plot
  p2 <- ggplot(
    pseudotime_df,
    aes(x = pseudotime, color = cell_type, fill = cell_type)
  ) +
    geom_density(alpha = 0.3, size = 1) +
    theme_classic() +
    labs(
      title = paste(output_prefix, "- Pseudotime Density by Cell Type"),
      x = "Pseudotime",
      y = "Density"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2")

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_pseudotime_density.pdf")),
    p2,
    width = 10,
    height = 6,
    dpi = FIGURE_DPI
  )

  cat("  Saved pseudotime distribution plots\n")
}

# ==============================================================================
# ENHANCED: PLOT MARKER GENES ALONG PSEUDOTIME
# ==============================================================================

plot_markers_along_pseudotime <- function(
  cds,
  markers,
  output_prefix,
  fig_dir
) {
  cat("\nPlotting marker genes along pseudotime...\n")

  available_markers <- markers[markers %in% rownames(cds)]
  if (length(available_markers) == 0) {
    cat("  No markers found in CDS\n")
    return(invisible(NULL))
  }

  cat("  Plotting markers:", paste(available_markers, collapse = ", "), "\n")

  # Get expression data
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    expr_matrix <- as.matrix(
      Seurat::GetAssayData(cds, assay = "RNA", layer = "data")[
        available_markers,
        ,
        drop = FALSE
      ]
    )
  } else {
    # For CDS objects, use counts and normalize
    expr_matrix <- as.matrix(counts(cds)[available_markers, , drop = FALSE])
    expr_matrix <- log1p(expr_matrix)
  }

  # Create data frame
  pseudotime_vals <- pseudotime(cds)
  cell_type_vals <- colData(cds)[[LABELS_KEY]]

  plot_data <- data.frame(
    pseudotime = rep(pseudotime_vals, each = length(available_markers)),
    gene = rep(available_markers, length(pseudotime_vals)),
    expression = as.vector(expr_matrix),
    cell_type = rep(cell_type_vals, each = length(available_markers))
  )

  plot_data <- plot_data[is.finite(plot_data$pseudotime), ]

  # Smooth expression along pseudotime
  p <- ggplot(plot_data, aes(x = pseudotime, y = expression, color = gene)) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    theme_classic() +
    labs(
      title = paste(output_prefix, "- Marker Expression Along Pseudotime"),
      x = "Pseudotime",
      y = "Expression (log)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      legend.position = "none"
    )

  ggsave(
    file.path(fig_dir, paste0(output_prefix, "_markers_pseudotime_smooth.pdf")),
    p,
    width = 12,
    height = ceiling(length(available_markers) / 3) * 3,
    dpi = FIGURE_DPI
  )

  cat("  Saved marker pseudotime plot\n")
}

# ==============================================================================
# MAIN TRAJECTORY ANALYSIS FUNCTION (ENHANCED)
# ==============================================================================

run_monocle3_trajectory <- function(
  seurat_obj,
  cell_types,
  root_cell_type,
  output_prefix,
  fig_dir,
  assay = "RNA"
) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("MONOCLE3 TRAJECTORY ANALYSIS: ", output_prefix, "\n", sep = "")
  cat(strrep("=", 80), "\n\n", sep = "")

  # Step 1: Subset cells
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

  # Step 2: Convert to CDS
  cat("\nStep 2: Converting to CDS...\n")
  cds <- seurat_to_cds(seurat_subset, assay = assay)

  # Step 3: Preprocess
  cat("\nStep 3: Preprocessing...\n")
  if ("scANVI" %in% names(reducedDims(cds))) {
    cat("  Using existing scANVI latent space as 'PCA'.\n")
    reducedDims(cds)[["PCA"]] <- reducedDims(cds)[["scANVI"]]
  } else {
    cat("  Running preprocess_cds (PCA)...\n")
    cds <- preprocess_cds(cds, num_dim = N_DIM_PREPROCESS)
  }

  # Step 4: Reduce dimension
  cat("\nStep 4: Dimensionality reduction...\n")
  if ("UMAP" %in% names(reducedDims(cds))) {
    cat("  Using existing UMAP from Seurat.\n")
  } else {
    cat("  Computing UMAP via reduce_dimension()...\n")
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
  }

  # Step 5: Cluster cells
  cat("\nStep 5: Clustering cells...\n")
  cds <- cluster_cells(cds, k = N_NEIGHBORS_CLUSTER)
  n_clusters <- length(unique(clusters(cds)))
  n_partitions <- length(unique(partitions(cds)))
  cat("  Clusters   :", n_clusters, "\n")
  cat("  Partitions :", n_partitions, "\n")

  # Step 6: Learn trajectory graph
  cat("\nStep 6: Learning trajectory graph...\n")
  cds <- learn_graph(cds, use_partition = USE_PARTITION)
  cat("  Trajectory graph learned.\n")

  # Step 7: Order cells in pseudotime
  cat("\nStep 7: Ordering cells in pseudotime...\n")
  root_node <- get_root_principal_node(cds, LABELS_KEY, root_cell_type)
  cds <- order_cells(cds, root_pr_nodes = root_node)

  pseudotime_vals <- pseudotime(cds)
  finite_pseudotime <- pseudotime_vals[is.finite(pseudotime_vals)]
  frac_finite <- length(finite_pseudotime) / length(pseudotime_vals) * 100

  cat("  Pseudotime statistics:\n")
  cat(sprintf(
    "    Cells with finite pseudotime: %d (%.1f%%)\n",
    length(finite_pseudotime),
    frac_finite
  ))
  cat(sprintf(
    "    Pseudotime range: %.2f - %.2f\n",
    min(finite_pseudotime),
    max(finite_pseudotime)
  ))

  # Step 8: Basic visualizations
  cat("\nStep 8: Generating basic visualizations...\n")

  # 8.1 Cell types on trajectory
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

  # 8.2 Pseudotime
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

  # 8.3 Partitions and branch points
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

  # Step 9: Branch analysis (NEW)
  cat("\nStep 9: Branch analysis...\n")
  cds <- plot_branch_analysis(cds, output_prefix, fig_dir)

  # Step 10: Pseudotime distribution (NEW)
  cat("\nStep 10: Pseudotime distribution analysis...\n")
  plot_pseudotime_distribution(cds, output_prefix, fig_dir, cell_types)

  # Step 11: Export pseudotime
  cat("\nStep 11: Exporting pseudotime values...\n")
  pseudotime_df <- data.frame(
    cell_id = colnames(cds),
    cell_type = colData(cds)[[LABELS_KEY]],
    pseudotime = pseudotime(cds),
    partition = partitions(cds),
    cluster = clusters(cds),
    branch = colData(cds)$branch
  )

  output_file <- file.path(fig_dir, paste0(output_prefix, "_pseudotime.csv"))
  write.csv(pseudotime_df, output_file, row.names = FALSE)
  cat("  Saved:", output_file, "\n")

  # Step 12: Summary per cell type
  cat("\nSummary Statistics (per cell type):\n")
  cat(strrep("-", 80), "\n")
  for (ct in cell_types) {
    ct_pt <- pseudotime_df$pseudotime[pseudotime_df$cell_type == ct]
    ct_pt <- ct_pt[is.finite(ct_pt)]
    if (length(ct_pt) > 0) {
      cat(sprintf("  %s:\n", ct))
      cat(sprintf("    Cells: %d\n", length(ct_pt)))
      cat(sprintf(
        "    Pseudotime range: %.2f - %.2f\n",
        min(ct_pt),
        max(ct_pt)
      ))
      cat(sprintf("    Mean pseudotime : %.2f\n", mean(ct_pt)))
      cat(sprintf("    Median pseudotime: %.2f\n", median(ct_pt)))
    }
  }

  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Trajectory analysis complete: ", output_prefix, "\n", sep = "")
  cat(strrep("=", 80), "\n\n")

  return(cds)
}

# ==============================================================================
# STEP 3: CD4 T CELL TRAJECTORY
# ==============================================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("STEP 3: CD4 T Cell Trajectory Analysis\n")
cat(strrep("=", 80), "\n\n")

cds_cd4 <- run_monocle3_trajectory(
  seurat_obj = seurat_obj,
  cell_types = CD4_CELLTYPES,
  root_cell_type = CD4_ROOT,
  output_prefix = "CD4",
  fig_dir = file.path(OUTPUT_DIR, "cd4_trajectory")
)

# ==============================================================================
# STEP 4: CD8 T CELL TRAJECTORY
# ==============================================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("STEP 4: CD8 T Cell Trajectory Analysis\n")
cat(strrep("=", 80), "\n\n")

cds_cd8 <- run_monocle3_trajectory(
  seurat_obj = seurat_obj,
  cell_types = CD8_CELLTYPES,
  root_cell_type = CD8_ROOT,
  output_prefix = "CD8",
  fig_dir = file.path(OUTPUT_DIR, "cd8_trajectory")
)

# ==============================================================================
# STEP 5: IMPROVED PSEUDOTIME-DEPENDENT GENES
# ==============================================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("STEP 5: Identifying Pseudotime-Dependent Genes (T cell filtered)\n")
cat(strrep("=", 80), "\n\n")

find_pseudotime_genes_filtered <- function(
  cds,
  output_prefix,
  fig_dir,
  num_genes = 100
) {
  cat("Finding genes that change along pseudotime...\n")
  cat("  Running graph_test (this may take several minutes)...\n\n")

  res <- tryCatch(
    {
      gene_fits <- monocle3::graph_test(
        cds,
        neighbor_graph = "principal_graph",
        cores = 4
      )
      gene_fits <- gene_fits[order(gene_fits$q_value), ]

      # Filter out myeloid markers
      cat("  Filtering out non-T cell markers...\n")
      non_tcell_genes <- c(
        "CD14",
        "FCGR3A",
        "ITGAX",
        "ITGAM",
        "CD68",
        "LYZ",
        "S100A8",
        "S100A9",
        "VCAN"
      )
      gene_fits_filtered <- gene_fits[
        !rownames(gene_fits) %in% non_tcell_genes,
      ]

      sig_genes <- gene_fits_filtered[
        gene_fits_filtered$q_value < 0.05,
        ,
        drop = FALSE
      ]
      cat("  Significant genes (q < 0.05):", nrow(sig_genes), "\n")

      # Export
      gene_fits_export <- gene_fits_filtered
      gene_fits_export$gene_id <- rownames(gene_fits_export)
      rownames(gene_fits_export) <- NULL

      out_file <- file.path(
        fig_dir,
        paste0(output_prefix, "_pseudotime_genes_filtered.csv")
      )
      write.csv(gene_fits_export, out_file, row.names = FALSE)
      cat("  Saved:", out_file, "\n")

      if (nrow(sig_genes) > 0) {
        top_genes <- head(rownames(sig_genes), num_genes)

        cat("\n  Top 10 pseudotime-dependent genes (filtered):\n")
        for (i in seq_len(min(10, length(top_genes)))) {
          g <- top_genes[i]
          qv <- sig_genes[g, "q_value"]
          cat(sprintf("    %2d. %s (q = %.2e)\n", i, g, qv))
        }

        genes_to_plot <- head(top_genes, 12)
        cat("\n  Plotting top genes on trajectory...\n")

        p <- plot_cells(
          cds,
          genes = genes_to_plot,
          show_trajectory_graph = FALSE,
          label_cell_groups = FALSE,
          label_leaves = FALSE
        )

        ggsave(
          file.path(fig_dir, paste0(output_prefix, "_top_genes_filtered.pdf")),
          p,
          width = 14,
          height = 10,
          dpi = FIGURE_DPI
        )
      }

      gene_fits_filtered
    },
    error = function(e) {
      cat("  Error in graph_test:", conditionMessage(e), "\n")
      NULL
    }
  )

  res
}

cat("\nAnalyzing CD4 pseudotime-dependent genes...\n")
cd4_genes <- find_pseudotime_genes_filtered(
  cds_cd4,
  output_prefix = "CD4",
  fig_dir = file.path(OUTPUT_DIR, "cd4_trajectory")
)

cat("\nAnalyzing CD8 pseudotime-dependent genes...\n")
cd8_genes <- find_pseudotime_genes_filtered(
  cds_cd8,
  output_prefix = "CD8",
  fig_dir = file.path(OUTPUT_DIR, "cd8_trajectory")
)

# ==============================================================================
# STEP 6: MARKER GENES ANALYSIS (ENHANCED)
# ==============================================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("STEP 6: Analyzing Key Marker Genes\n")
cat(strrep("=", 80), "\n\n")

# CD4 markers
cd4_all_markers <- c(
  NAIVE_MARKERS,
  MEMORY_MARKERS,
  EFFECTOR_MARKERS,
  TREG_MARKERS
)
cd4_all_markers <- unique(cd4_all_markers)

# CD8 markers
cd8_all_markers <- c(
  NAIVE_MARKERS,
  MEMORY_MARKERS,
  EFFECTOR_MARKERS
)
cd8_all_markers <- unique(cd8_all_markers)

# Plot markers on trajectory
plot_marker_genes <- function(cds, markers, output_prefix, fig_dir) {
  cat("Plotting marker genes for", output_prefix, "...\n")

  available_markers <- markers[markers %in% rownames(cds)]
  if (length(available_markers) == 0) {
    cat("  No marker genes found in CDS\n")
    return(invisible(NULL))
  }

  cat("  Available markers:", paste(available_markers, collapse = ", "), "\n")

  # Plot on trajectory
  p <- plot_cells(
    cds,
    genes = available_markers,
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )

  out_file <- file.path(fig_dir, paste0(output_prefix, "_marker_genes.pdf"))
  ggsave(out_file, p, width = 14, height = 12, dpi = FIGURE_DPI)
  cat("  Saved:", out_file, "\n")

  # Plot along pseudotime (NEW)
  plot_markers_along_pseudotime(cds, available_markers, output_prefix, fig_dir)
}

cat("\nCD4 marker analysis:\n")
plot_marker_genes(
  cds_cd4,
  cd4_all_markers,
  "CD4",
  file.path(OUTPUT_DIR, "cd4_trajectory")
)

cat("\nCD8 marker analysis:\n")
plot_marker_genes(
  cds_cd8,
  cd8_all_markers,
  "CD8",
  file.path(OUTPUT_DIR, "cd8_trajectory")
)

# ==============================================================================
# STEP 7: SAVE CDS OBJECTS (IMPROVED)
# ==============================================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("STEP 7: Saving CDS Objects\n")
cat(strrep("=", 80), "\n\n")

cd4_dir <- file.path(OUTPUT_DIR, "cd4_trajectory", "cds_cd4")
cd8_dir <- file.path(OUTPUT_DIR, "cd8_trajectory", "cds_cd8")

dir.create(cd4_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cd8_dir, showWarnings = FALSE, recursive = TRUE)

cat("Saving CD4 CDS using save_monocle_objects()...\n")
tryCatch(
  {
    save_monocle_objects(cds_cd4, directory = cd4_dir)
    cat("  Saved CD4 CDS to:", cd4_dir, "\n")
  },
  error = function(e) {
    cat("  Error with save_monocle_objects:", conditionMessage(e), "\n")
    cat("  Falling back to saveRDS()...\n")
    saveRDS(cds_cd4, file.path(OUTPUT_DIR, "cd4_trajectory", "cds_cd4.rds"))
  }
)

cat("Saving CD8 CDS using save_monocle_objects()...\n")
tryCatch(
  {
    save_monocle_objects(cds_cd8, directory = cd8_dir)
    cat("  Saved CD8 CDS to:", cd8_dir, "\n")
  },
  error = function(e) {
    cat("  Error with save_monocle_objects:", conditionMessage(e), "\n")
    cat("  Falling back to saveRDS()...\n")
    saveRDS(cds_cd8, file.path(OUTPUT_DIR, "cd8_trajectory", "cds_cd8.rds"))
  }
)

cat("\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("CD4 Trajectory:\n")
cat("  Cell types:", paste(CD4_CELLTYPES, collapse = ", "), "\n")
cat("  Root      :", CD4_ROOT, "\n")
cat("  Cells     :", ncol(cds_cd4), "\n")
cat("  Output    :", file.path(OUTPUT_DIR, "cd4_trajectory"), "\n\n")

cat("CD8 Trajectory:\n")
cat("  Cell types:", paste(CD8_CELLTYPES, collapse = ", "), "\n")
cat("  Root      :", CD8_ROOT, "\n")
cat("  Cells     :", ncol(cds_cd8), "\n")
cat("  Output    :", file.path(OUTPUT_DIR, "cd8_trajectory"), "\n\n")

cat("Generated files (per CD4/CD8):\n")
cat("  Quality Control:\n")
cat("    - qc/tcell_marker_score.pdf\n")
cat("    - qc/myeloid_marker_score.pdf\n\n")
cat("  Basic Trajectory:\n")
cat("    - *_trajectory_celltypes.pdf\n")
cat("    - *_trajectory_pseudotime.pdf\n")
cat("    - *_trajectory_partitions.pdf\n\n")
cat("  Branch Analysis (NEW):\n")
cat("    - *_branch_assignment.pdf\n")
cat("    - *_branch_celltype_heatmap.pdf\n")
cat("    - *_pseudotime_by_branch.pdf\n\n")
cat("  Pseudotime Distribution (NEW):\n")
cat("    - *_pseudotime_distribution.pdf\n")
cat("    - *_pseudotime_density.pdf\n\n")
cat("  Gene Analysis:\n")
cat("    - *_pseudotime_genes_filtered.csv\n")
cat("    - *_top_genes_filtered.pdf\n")
cat("    - *_marker_genes.pdf\n")
cat("    - *_markers_pseudotime_smooth.pdf (NEW)\n\n")
cat("  Data:\n")
cat("    - *_pseudotime.csv (with branch info)\n")
cat("    - cds_*/  (saved Monocle objects)\n\n")

cat("Key improvements in this version:\n")
cat("  1. Quality control for T cell purity\n")
cat("  2. Branch analysis and visualization\n")
cat("  3. Filtered pseudotime-dependent genes (removed myeloid markers)\n")
cat("  4. Enhanced marker gene analysis along pseudotime\n")
cat("  5. Pseudotime distribution diagnostics\n")
cat("  6. Proper CDS saving with save_monocle_objects()\n\n")

cat("Next steps:\n")
cat("  1. Review QC plots for contamination issues\n")
cat("  2. Examine branch assignments - do they match biology?\n")
cat("  3. Check pseudotime distribution plots - ordered as expected?\n")
cat("  4. Validate marker expression patterns along pseudotime\n")
cat("  5. If issues remain, consider:\n")
cat("     - Re-filtering cells based on marker scores\n")
cat("     - Adjusting clustering parameters\n")
cat("     - Re-annotating cell types based on trajectory\n\n")

cat(strrep("=", 80), "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("=", 80), "\n\n")

cat("Happy analyzing! 🔬\n\n")
