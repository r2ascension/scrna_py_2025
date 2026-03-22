#!/usr/bin/env Rscript
################################################################################
# Basal Cell Trajectory Analysis - Multi-Tissue with Hierarchical Annotations
# Version: v2.8.3.1 HOTFIX - Seurat v5 Compatible
#
# HOTFIX IN v2.8.3.1:
# - FIXED Seurat v5 Assay5 compatibility in check_root_validity()
# - Error: "no slot of name 'data' for this object of class 'Assay5'"
# - Solution: Use GetAssayData() for Seurat v5, fallback to @data for v4
#
# CHANGES IN v2.8.3 (Conservative improvements based on code review):
# 1. ROOT CHANGED: cluster 9 (Cycling) → cluster 1 (Basal_Core)
# 2. ADDED root validity check: outputs marker enrichment for each tissue
# 3. ADDED unmapped cluster warnings: prevents silent NA annotations
# 4. ADDED skipped tissues log: transparency for failed analyses
# 5. OPTIMIZED branch assignment: vertex-level computation (faster)
# 6. DISABLED custom tree plots: coordinate system bug (use standard monocle3)
# 7. CLUSTER 11 RETAINED: verified as bona fide basal (SAA1+/S100A9+ wound response)
#
# Based on: basal_monocle3_trajectory_v2_7_PRODUCTION.R
################################################################################

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(pheatmap)
  library(igraph)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(grid)
  library(patchwork)
  library(cowplot)
  library(RColorBrewer)
})

cat("\n", strrep("=", 80), "\n", sep = "")
cat("BASAL CELL MULTI-TISSUE TRAJECTORY ANALYSIS - v2.8.3.1 HOTFIX\n")
cat("(Seurat v5 compatible, root at Basal_Core, production ready)\n")
cat(strrep("=", 80), "\n\n", sep = "")
cat("monocle3:", as.character(packageVersion("monocle3")), "\n")
cat("Seurat  :", as.character(packageVersion("Seurat")), "\n")
cat("igraph  :", as.character(packageVersion("igraph")), "\n")
cat("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# CONFIG - v2.8.3 CONSERVATIVE PRODUCTION VERSION
# ==============================================================================
# ROOT CHANGED: cluster 9 (Cycling) → cluster 1 (Basal_Core)
# Rationale: Classic stem cell differentiation model
#   - Basal_Core (TP63+/KRT15+/SELM+) represents stem-like basal cells
#   - Trajectory: Core → Primed/Squamous/EMT/Wound response
#   - Better for publication interpretation vs cycling-centric view
#
# CLUSTER 16 REMOVED: IGHG1+/IGJ+/LTF+ immune signatures
# CLUSTER 11 RETAINED: SAA1+/S100A9+ verified as basal wound response
#
# NEW FEATURES:
# - Root validity check: outputs marker enrichment analysis
# - Unmapped cluster warnings: prevents silent NA annotations  
# - Skipped tissues log: tracks failed analyses
# - Optimized branch assignment: vertex-level (faster)
# - Standard monocle3 trajectory graphs (custom tree disabled)
#
# Level 3 annotations: 11 functional groups (optimized from marker analysis)
# Level 4 annotations: 15 detailed types
# ==============================================================================

INPUT_RDS <- "/home/h2048/data/R/1228/basal/basal_filtered_20251228.rds"
OUTPUT_DIR <- "/home/h2048/data/R/1228/basal/Basal_MultiTissue_Trajectory_v2_8_20250103"

# Workflow control
SKIP_PART1 <- FALSE
LOAD_EXISTING_CDS <- FALSE

# ⭐ KEY: Root cluster and tissue handling
ROOT_LABEL <- "1"  # Basal_Core_SELM+DLK2+ (stem-like basal, changed from cycling)
CLUSTER_COL <- "RNA_snn_res.1"
TISSUE_KEY <- "tissue"  # Set to actual tissue column name, or NULL for single dataset

# ⭐ HIERARCHICAL CELL TYPE ANNOTATIONS - OPTIMIZED based on marker genes
# Level 4: 15 cluster types (cluster 16 removed due to immune contamination)
# Level 3: Functional grouping (11 categories, optimized from marker gene analysis)
#
# CLUSTER 16 REMOVED: IGHG1+/IGJ+/LTF+ signatures indicate immune contamination
# rather than bona fide basal cell program
#
# L3 Grouping Logic (based on top marker genes from FindAllMarkers):
# - Basal_Core (1,3): SELM+/DLK2+/KRT15+ stem-like basal cells
# - Basal_Primed (2,14,15): SCGB1A1+/PIGR+/MUC2+ secretory/mucous-primed
# - Basal_Squamous (5,6): SERPINB3/B4+/KRT13+/KRT6A+ squamous differentiation
# - Basal_Inflammatory (4): CXCL2+/IL8+/IEG+ inflammatory response
# - Basal_ECM_Signaling (7): TGFB2+/WNT4+/FN1+ ECM signaling pathways
# - Basal_Metabolic (8): ATP synthases, translation machinery
# - Basal_Cycling (9): MKI67+/TOP2A+ proliferating cells (ROOT)
# - Basal_Antimicrobial (10): WFDC2+/SLPI+ antimicrobial defense
# - Basal_Wound_Response (11): SAA1+/S100A9+/COL17A1+ acute wound
# - Basal_ECM_Adhesion (12): FN1+/COL4A5+/ANK3+ ECM remodeling
# - Basal_EMT_Repair (13): VIM+/CTGF+/SOX9+ EMT/wound repair

# Level 4: Full 16 cluster names (most detailed) - Cluster 16 removed
CLUSTER_TO_CELLTYPE_L4 <- c(
  "1"  = "Basal_Core_SELM+DLK2+",
  "2"  = "Basal_SecretoryPrimed_SCGB1A1_PIGR",
  "3"  = "Basal_Core_KRT15+",
  "4"  = "Basal_Inflammatory_IEG",
  "5"  = "Basal_Squamous_SERPINB3B4",
  "6"  = "Basal_Squamous_KRT13KRT6A",
  "7"  = "Basal_TGFB_WNT_ECM",
  "8"  = "Basal_Metabolic_TranslationHigh",
  "9"  = "Basal_Cycling",
  "10" = "Basal_Antimicrobial_WFDC2_SLPI",
  "11" = "Basal_AcuteStress_Wound",
  "12" = "Basal_ECM_Adhesion",
  "13" = "Basal_WoundRepair_EMT",
  "14" = "Basal_MucousPrimed_MUC2",
  "15" = "Basal_SCGB2B2_LncRNA"
)

# Level 3: Functional grouping (intermediate) - OPTIMIZED based on marker genes
# NOTE: Cluster 16 removed due to suspected immune contamination (IGHG1+/IGJ+)
CLUSTER_TO_CELLTYPE_L3 <- c(
  "1"  = "Basal_Core",              # SELM+, DLK2+ core basal
  "2"  = "Basal_Primed",             # SCGB1A1+, PIGR+ secretory-primed
  "3"  = "Basal_Core",               # KRT15+ core basal
  "4"  = "Basal_Inflammatory",       # CXCL2+, IL8+, IEG+ inflammatory
  "5"  = "Basal_Squamous",           # SERPINB3/B4+ squamous
  "6"  = "Basal_Squamous",           # KRT13+, KRT6A+ squamous
  "7"  = "Basal_ECM_Signaling",      # TGFB2+, WNT4+, FN1+ ECM/signaling
  "8"  = "Basal_Metabolic",          # ATP synthases, translation-high
  "9"  = "Basal_Cycling",            # MKI67+, TOP2A+ proliferating (ROOT)
  "10" = "Basal_Antimicrobial",      # WFDC2+, SLPI+ antimicrobial
  "11" = "Basal_Wound_Response",     # SAA1+, S100A9+, COL17A1+ wound
  "12" = "Basal_ECM_Adhesion",       # FN1+, COL4A5+, ANK3+ ECM/adhesion
  "13" = "Basal_EMT_Repair",         # VIM+, CTGF+, SOX9+ EMT/wound repair
  "14" = "Basal_Primed",             # MUC2+ mucous-primed
  "15" = "Basal_Primed"              # SCGB2B2+, lncRNA-enriched
)

# Color schemes (match Science paper style)
# Color schemes (match Science paper style) - 11 L3 categories (cluster 16 removed)
CELLTYPE_L3_COLORS <- c(
  "Basal_Core"              = "#0bb1da",  # Cyan - stem-like basal
  "Basal_Primed"            = "#32b34e",  # Green - secretory/mucous primed
  "Basal_Squamous"          = "#ff527d",  # Pink - squamous differentiation
  "Basal_Inflammatory"      = "#ff6600",  # Orange - IEG/inflammatory
  "Basal_ECM_Signaling"     = "#9966cc",  # Purple - TGFB/WNT/ECM signaling
  "Basal_Metabolic"         = "#7a00cc",  # Deep purple - metabolic/translation
  "Basal_Cycling"           = "#3333ff",  # Blue - proliferating (ROOT)
  "Basal_Antimicrobial"     = "#00cc99",  # Teal - WFDC2/SLPI antimicrobial
  "Basal_Wound_Response"    = "#ffcc00",  # Yellow - SAA1/acute stress
  "Basal_ECM_Adhesion"      = "#cc0099",  # Magenta - ECM/adhesion
  "Basal_EMT_Repair"        = "#cc6600"   # Brown - VIM/EMT/wound repair
)

CELLTYPE_L4_COLORS <- colorRampPalette(brewer.pal(12, "Set3"))(15)  # 15 clusters after removing 16
names(CELLTYPE_L4_COLORS) <- unname(CLUSTER_TO_CELLTYPE_L4)  # FIX: use celltype names, not cluster IDs

# Monocle3 parameters
N_DIM_PREPROCESS <- 50
N_NEIGHBORS_CLUSTER <- 15
USE_PARTITION <- FALSE
RUN_GRAPH_TEST <- TRUE
REDUCTION_METHOD <- "UMAP"
USE_SEURAT_UMAP <- FALSE

# Tree visualization parameters
GENERATE_TREE_PLOTS <- FALSE  # DISABLED: Custom tree has coordinate system bug (FR layout vs UMAP)
                               # Use standard monocle3::plot_cells(show_trajectory_graph=TRUE) instead
TREE_CELL_SIZE <- 0.5
TREE_LINK_SIZE <- 0.8

# Analysis parameters
TOP_N_MARKERS <- 10
MIN_PSEUDOBULK_CELLS <- 10
MAX_CELLS_MARKER_TRENDS <- 20000
GENERATE_BRANCH_HEATMAPS <- TRUE
N_GENES_HEATMAP <- 50
MIN_CELLS_PER_BRANCH <- 50
N_GENE_CLUSTERS <- 3
USE_BINNED_HEATMAP <- TRUE
N_PSEUDOTIME_BINS <- 100

# Lineage markers - UPDATED to match optimized L3 groups
LINEAGE_MARKERS <- list(
  Basal_Core = c("KRT5", "TP63", "KRT15", "DLK2", "SELM"),
  Basal_Cycling = c("MKI67", "TOP2A", "PCNA", "BIRC5", "UBE2C"),
  Basal_Squamous = c("KRT13", "KRT6A", "SERPINB3", "SERPINB4"),
  Basal_Primed = c("SCGB1A1", "PIGR", "MUC2", "SCGB2B2"),
  Basal_Inflammatory = c("CXCL2", "IL8", "AREG", "EGR2", "FOSL1"),
  Basal_Wound_EMT = c("SAA1", "VIM", "CTGF", "SOX9", "FN1"),
  Basal_ECM = c("FN1", "COL4A5", "TGFB2", "WNT4", "ANK3"),
  Basal_Antimicrobial = c("WFDC2", "SLPI", "LTF", "DEFB1")
)

# ==============================================================================
# Root validity check
# ==============================================================================
check_root_validity <- function(seurat_obj, root_cluster, cluster_col) {
  cat("  Validating root cluster selection...\n")
  
  # Get cells in root cluster
  root_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cluster_col]] == root_cluster]
  other_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cluster_col]] != root_cluster]
  
  if (length(root_cells) == 0) {
    warning("Root cluster ", root_cluster, " has no cells!")
    return(NULL)
  }
  
  # Get expression data (Seurat v5 compatible)
  if (inherits(seurat_obj@assays$RNA, "Assay5")) {
    # Seurat v5: use LayerData or GetAssayData
    expr_data <- tryCatch({
      Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    }, error = function(e) {
      # Fallback to counts if data layer doesn't exist
      Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    })
  } else {
    # Seurat v4: old style
    expr_data <- seurat_obj@assays$RNA@data
  }
  
  # Calculate module scores for different programs
  results <- list()
  
  for (program in names(LINEAGE_MARKERS)) {
    markers <- LINEAGE_MARKERS[[program]]
    markers_present <- markers[markers %in% rownames(expr_data)]
    
    if (length(markers_present) >= 2) {
      # Calculate mean expression for root vs others
      root_expr <- colMeans(as.matrix(expr_data[markers_present, root_cells, drop = FALSE]))
      other_expr <- colMeans(as.matrix(expr_data[markers_present, other_cells, drop = FALSE]))
      
      results[[program]] <- data.frame(
        Program = program,
        Root_Mean = mean(root_expr),
        Other_Mean = mean(other_expr),
        Fold_Enrichment = mean(root_expr) / (mean(other_expr) + 0.01)
      )
    }
  }
  
  result_df <- do.call(rbind, results)
  result_df <- result_df[order(-result_df$Fold_Enrichment), ]
  
  cat("  Root cluster ", root_cluster, " marker enrichment:\n", sep = "")
  print(result_df, row.names = FALSE)
  
  # Check if root has expected signature
  top_program <- result_df$Program[1]
  expected_root_programs <- c("Basal_Core", "Basal_Cycling")
  
  if (top_program %in% expected_root_programs) {
    cat("  ✓ Root shows expected signature: ", top_program, "\n", sep = "")
  } else {
    warning("  ⚠ Root enriched for: ", top_program, " (expected Core or Cycling)")
  }
  
  cat("\n")
  return(result_df)
}

FIG_W <- 10
FIG_H <- 8

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Helpers
# ==============================================================================
sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

safe_write_log <- function(line, path) {
  tryCatch({
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    write(line, file = path, append = TRUE)
    write("\n", file = path, append = TRUE)
  }, error = function(e) NULL)
}

DEBUG_LOG_PATH <- "/home/h2048/.cursor/debug.log"

save_plot_pdf <- function(p, file, width = FIG_W, height = FIG_H) {
  pdf(file, width = width, height = height, onefile = TRUE)
  print(p)
  dev.off()
}

# ==============================================================================
# Add hierarchical cell type annotations
# ==============================================================================
# ==============================================================================
# Add hierarchical cell type annotations
# ==============================================================================
add_celltype_annotations <- function(seurat_obj, cluster_col, map_l3, map_l4) {
  clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  
  # Level 3 (functional groups)
  celltype_l3 <- map_l3[clusters]
  
  # Check for unmapped clusters
  unmapped_l3 <- setdiff(unique(clusters), names(map_l3))
  if (length(unmapped_l3) > 0) {
    warning("⚠ Unmapped clusters in Level 3: ", paste(unmapped_l3, collapse=", "))
    warning("  These clusters will have NA cell type annotations.")
  }
  
  seurat_obj@meta.data$cell_type_level_3 <- celltype_l3
  
  # Level 4 (detailed)
  celltype_l4 <- map_l4[clusters]
  
  # Check for unmapped clusters
  unmapped_l4 <- setdiff(unique(clusters), names(map_l4))
  if (length(unmapped_l4) > 0) {
    warning("⚠ Unmapped clusters in Level 4: ", paste(unmapped_l4, collapse=", "))
    warning("  These clusters will have NA cell type annotations.")
  }
  
  seurat_obj@meta.data$cell_type_level_4 <- celltype_l4
  
  cat("  Added hierarchical annotations:\n")
  cat("    Level 3 (functional): ", length(unique(na.omit(celltype_l3))), " groups\n", sep = "")
  cat("    Level 4 (detailed)  : ", length(unique(na.omit(celltype_l4))), " types\n\n", sep = "")
  
  return(seurat_obj)
}

# ==============================================================================
# Data conversion
# ==============================================================================
get_counts_from_seurat <- function(seurat_obj, assay = "RNA") {
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    x <- tryCatch(
      Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts"),
      error = function(e) NULL
    )
    if (!is.null(x)) return(x)
  }
  x <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts"),
    error = function(e) NULL
  )
  if (is.null(x)) {
    stop("Cannot find raw counts in assay=", assay)
  }
  x
}

seurat_to_cds <- function(seurat_obj, assay = "RNA", use_seurat_umap = FALSE) {
  counts <- get_counts_from_seurat(seurat_obj, assay = assay)
  
  cell_md <- seurat_obj@meta.data
  gene_md <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )
  
  cds <- monocle3::new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_md,
    gene_metadata = gene_md
  )
  
  if (use_seurat_umap) {
    red_names <- names(seurat_obj@reductions)
    umap_red <- red_names[grepl("umap", red_names, ignore.case = TRUE)][1]
    if (!is.na(umap_red) && !is.null(umap_red)) {
      umap_coords <- Seurat::Embeddings(seurat_obj, reduction = umap_red)
      colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
      reducedDims(cds)[["UMAP"]] <- umap_coords
      message("  [INFO] Using Seurat UMAP: ", umap_red)
    }
  } else {
    message("  [INFO] Monocle3 will compute its own UMAP")
  }
  
  cds
}

pick_root_label <- function(labels_vec, root_label = "9") {
  labs <- unique(as.character(labels_vec))
  if (root_label %in% labs) {
    return(root_label)
  }
  message("  [ERROR] Root cluster '", root_label, "' not found")
  message("  Available clusters: ", paste(labs, collapse = ", "))
  return(NULL)
}

get_root_principal_node <- function(cds, cell_type_key, root_cell_type) {
  root_cells <- which(colData(cds)[[cell_type_key]] == root_cell_type)
  if (length(root_cells) == 0) {
    stop("No cells found for root cell type: ", root_cell_type)
  }
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])[, 1]
  
  pr_graph <- monocle3::principal_graph(cds)[["UMAP"]]
  vertex_names <- igraph::V(pr_graph)$name
  
  if (is.numeric(closest_vertex)) {
    closest_vertex_names <- vertex_names[closest_vertex]
  } else {
    closest_vertex_names <- as.character(closest_vertex)
  }
  
  root_vertex_name <- names(which.max(table(closest_vertex_names[root_cells])))
  return(root_vertex_name)
}

assign_branch_by_leaf_fast <- function(cds) {
  g <- monocle3::principal_graph(cds)[["UMAP"]]
  deg <- igraph::degree(g)
  leaves <- names(deg[deg == 1])
  
  if (length(leaves) < 2) {
    colData(cds)$branch <- "Main"
    return(cds)
  }
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])[, 1]
  
  pr_graph <- monocle3::principal_graph(cds)[["UMAP"]]
  vertex_names <- igraph::V(pr_graph)$name
  
  if (is.numeric(closest_vertex)) {
    closest_vertex_names <- vertex_names[closest_vertex]
  } else {
    closest_vertex_names <- as.character(closest_vertex)
  }
  
  # OPTIMIZATION: Compute distances at vertex level, then map to cells
  cat("    Computing vertex-to-leaf distances...\n")
  unique_vertices <- unique(closest_vertex_names)
  vertex_to_leaf_dist <- igraph::distances(g, v = unique_vertices, to = leaves)
  
  # Assign each vertex to its nearest leaf
  vertex_branch <- apply(vertex_to_leaf_dist, 1, which.min)
  names(vertex_branch) <- unique_vertices
  
  # Map cells to branches via their closest vertex
  cell_branch <- vertex_branch[closest_vertex_names]
  colData(cds)$branch <- paste0("Leaf_", cell_branch)
  
  cat("    Identified ", length(leaves), " trajectory branches\n", sep = "")
  cat("    Branch composition:\n")
  print(table(colData(cds)$branch))
  
  cds
}

# ==============================================================================
# ⭐ NEW: Tree structure visualization (inspired by Science paper)
# ==============================================================================
plot_trajectory_tree <- function(cds, 
                                color_by = "cell_type_level_3",
                                color_palette = NULL,
                                out_file,
                                title = "Trajectory Tree Structure",
                                cell_size = 0.5,
                                link_size = 0.8) {
  
  cat("  Generating tree structure plot...\n")
  
  # Extract principal graph
  pr_graph <- monocle3::principal_graph(cds)[["UMAP"]]
  
  # Get graph layout (vertex coordinates)
  graph_layout <- igraph::layout_with_fr(pr_graph)
  vertex_names <- igraph::V(pr_graph)$name
  
  # Create vertex dataframe
  vertex_df <- data.frame(
    vertex = vertex_names,
    x = graph_layout[, 1],
    y = graph_layout[, 2]
  )
  
  # Get edges
  edge_list <- igraph::as_edgelist(pr_graph)
  edge_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2]
  )
  
  # Merge with coordinates
  edge_df <- merge(edge_df, vertex_df, by.x = "from", by.y = "vertex")
  names(edge_df)[3:4] <- c("x_from", "y_from")
  edge_df <- merge(edge_df, vertex_df, by.x = "to", by.y = "vertex")
  names(edge_df)[5:6] <- c("x_to", "y_to")
  
  # Identify branch points and leaves
  deg <- igraph::degree(pr_graph)
  branch_points <- names(deg[deg > 2])
  leaves <- names(deg[deg == 1])
  
  vertex_df$type <- "intermediate"
  vertex_df$type[vertex_df$vertex %in% branch_points] <- "branch_point"
  vertex_df$type[vertex_df$vertex %in% leaves] <- "leaf"
  
  # Get cell coordinates and colors
  umap_coords <- reducedDims(cds)[["UMAP"]]
  cell_df <- data.frame(
    cell = colnames(cds),
    x = umap_coords[, 1],
    y = umap_coords[, 2],
    color_var = colData(cds)[[color_by]]
  )
  
  # Plot
  p <- ggplot() +
    # Cells (background)
    geom_point(data = cell_df, aes(x = x, y = y, color = color_var), 
               size = cell_size, alpha = 0.4) +
    # Graph edges
    geom_segment(data = edge_df, 
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 size = link_size, color = "black", alpha = 0.8) +
    # Graph vertices
    geom_point(data = vertex_df, aes(x = x, y = y, shape = type),
               size = 3, color = "black", fill = "white", stroke = 1.5) +
    scale_shape_manual(
      values = c("intermediate" = 21, "branch_point" = 23, "leaf" = 24),
      labels = c("Intermediate", "Branch Point", "Leaf")
    ) +
    labs(
      title = title,
      x = "Component 1",
      y = "Component 2",
      color = gsub("_", " ", color_by),
      shape = "Node Type"
    ) +
    theme_cowplot(font_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Add color palette if provided
  if (!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  }
  
  # Save
  pdf(out_file, width = 12, height = 10)
  print(p)
  dev.off()
  
  cat("  Saved: ", out_file, "\n", sep = "")
  
  return(p)
}

# Two-panel trajectory plot
save_two_panel_pdf <- function(p_left, p_right, file, title, width = 14, height = 6) {
  pdf(file, width = width, height = height, onefile = TRUE)
  grid::grid.newpage()
  
  lay <- grid::grid.layout(
    nrow = 2, ncol = 2,
    heights = grid::unit.c(grid::unit(0.8, "in"), grid::unit(1, "null")),
    widths  = grid::unit.c(grid::unit(1, "null"), grid::unit(1, "null"))
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  grid::grid.text(
    title,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2),
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )
  
  print(p_left, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  dev.off()
}

make_two_panel_trajectory <- function(cds, cell_type_key, out_file,
                                      reduction_method = "UMAP",
                                      title = "Basal Cell Trajectory",
                                      color_palette = NULL) {
  
  p_left <- monocle3::plot_cells(
    cds,
    reduction_method = reduction_method,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = TRUE,
    label_roots = FALSE
  ) +
    ggtitle(NULL) +
    labs(x = "Component 1", y = "Component 2", color = "Pseudotime") +
    scale_color_viridis_c(option = "plasma", na.value = "grey80") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
  
  p_right <- monocle3::plot_cells(
    cds,
    reduction_method = reduction_method,
    color_cells_by = cell_type_key,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE
  ) +
    ggtitle(NULL) +
    labs(x = "Component 1", y = "Component 2", color = "Cell type") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
  
  if (!is.null(color_palette)) {
    p_right <- p_right + scale_color_manual(values = color_palette)
  }
  
  save_two_panel_pdf(p_left, p_right, out_file, title = title)
  cat("  Saved: ", out_file, "\n", sep = "")
}

# [Include other helper functions from original script]
# ... (plot_pseudotime_markers, generate_branch_heatmap, etc.)

# ==============================================================================
# Main per-tissue trajectory analysis
# ==============================================================================
# ==============================================================================
# Main per-tissue trajectory analysis
# ==============================================================================
run_monocle3_for_tissue <- function(seurat_obj, tissue_value, out_dir, skipped_log = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE)
  
  # Subset by tissue
  md <- seurat_obj@meta.data
  cells <- rownames(md)[as.character(md[[TISSUE_KEY]]) == tissue_value]
  
  if (length(cells) < 200) {
    skip_msg <- paste0("tissue=", tissue_value, " | cells=", length(cells), " | reason=too_few_cells (<200)")
    cat("  [SKIP] ", skip_msg, "\n", sep = "")
    
    # Log to skipped file
    if (!is.null(skipped_log)) {
      write(skip_msg, file = skipped_log, append = TRUE)
    }
    
    return(NULL)
  }
  
  seu <- subset(seurat_obj, cells = cells)
  
  cat("  tissue=", tissue_value, " | cells=", ncol(seu), "\n", sep = "")
  
  # Validate root cluster
  root_label <- pick_root_label(seu@meta.data[[CLUSTER_COL]], ROOT_LABEL)
  if (is.null(root_label)) {
    skip_msg <- paste0("tissue=", tissue_value, " | cells=", ncol(seu), " | reason=root_cluster_not_found")
    cat("  [SKIP] ", skip_msg, "\n", sep = "")
    
    if (!is.null(skipped_log)) {
      write(skip_msg, file = skipped_log, append = TRUE)
    }
    
    return(NULL)
  }
  
  # ===== ROOT VALIDITY CHECK =====
  root_validity <- check_root_validity(seu, root_label, CLUSTER_COL)
  if (!is.null(root_validity)) {
    write.csv(root_validity, file.path(out_dir, "root_validity_check.csv"), row.names = FALSE)
  }
  
  # Convert to CDS
  cds <- seurat_to_cds(seu, use_seurat_umap = USE_SEURAT_UMAP)
  
  # Monocle3 workflow
  cat("  Running preprocess_cds...\n")
  cds <- preprocess_cds(cds, num_dim = N_DIM_PREPROCESS, method = "PCA")
  
  if (!("UMAP" %in% names(reducedDims(cds)))) {
    cat("  Running reduce_dimension...\n")
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
  }
  
  cat("  Running cluster_cells...\n")
  cds <- cluster_cells(cds, reduction_method = REDUCTION_METHOD, k = N_NEIGHBORS_CLUSTER)
  
  cat("  Running learn_graph...\n")
  cds <- learn_graph(cds, use_partition = USE_PARTITION)
  
  # Find root node
  cat("  Finding root node for cluster ", root_label, "...\n", sep = "")
  root_node <- get_root_principal_node(cds, CLUSTER_COL, root_label)
  cat("  Root node: ", root_node, "\n", sep = "")
  
  # Order cells
  cds <- order_cells(cds, root_pr_nodes = root_node)
  
  # Assign branches
  cat("  Assigning branches...\n")
  cds <- assign_branch_by_leaf_fast(cds)
  
  # Save metadata
  metadata(cds)$tissue <- tissue_value
  metadata(cds)$root_label <- root_label
  metadata(cds)$root_node <- root_node
  
  # ========== PLOTS ==========
  
  # Basic trajectory plots
  p1 <- plot_cells(cds, color_cells_by = "cell_type_level_3", 
                   label_cell_groups = FALSE) +
    scale_color_manual(values = CELLTYPE_L3_COLORS) +
    ggtitle(paste0(tissue_value, " - Level 3 Cell Types"))
  save_plot_pdf(p1, file.path(fig_dir, "trajectory_celltype_L3.pdf"))
  
  p2 <- plot_cells(cds, color_cells_by = "pseudotime", 
                   label_cell_groups = FALSE) +
    scale_color_viridis(na.value = "grey80") +
    ggtitle(paste0(tissue_value, " - Pseudotime"))
  save_plot_pdf(p2, file.path(fig_dir, "trajectory_pseudotime.pdf"))
  
  # Two-panel trajectory
  two_panel_file <- file.path(fig_dir, "trajectory_two_panel.pdf")
  make_two_panel_trajectory(
    cds, 
    "cell_type_level_3", 
    two_panel_file,
    title = paste0("Basal Cell Trajectory - ", tissue_value),
    color_palette = CELLTYPE_L3_COLORS
  )
  
  # ⭐ Tree structure plot
  if (GENERATE_TREE_PLOTS) {
    tree_file_L3 <- file.path(fig_dir, "trajectory_tree_level3.pdf")
    plot_trajectory_tree(
      cds,
      color_by = "cell_type_level_3",
      color_palette = CELLTYPE_L3_COLORS,
      out_file = tree_file_L3,
      title = paste0("Basal Trajectory Tree - ", tissue_value, " (Level 3)"),
      cell_size = TREE_CELL_SIZE,
      link_size = TREE_LINK_SIZE
    )
    
    tree_file_L4 <- file.path(fig_dir, "trajectory_tree_level4.pdf")
    plot_trajectory_tree(
      cds,
      color_by = "cell_type_level_4",
      color_palette = CELLTYPE_L4_COLORS,
      out_file = tree_file_L4,
      title = paste0("Basal Trajectory Tree - ", tissue_value, " (Level 4)"),
      cell_size = TREE_CELL_SIZE,
      link_size = TREE_LINK_SIZE
    )
  } else {
    # Use standard monocle3 trajectory visualization instead
    cat("  Generating standard trajectory graphs (custom tree plots disabled)...\n")
    
    # Level 3 with trajectory graph
    p_traj_L3 <- plot_cells(cds, 
                            color_cells_by = "cell_type_level_3",
                            label_cell_groups = FALSE,
                            label_leaves = TRUE,
                            label_branch_points = TRUE,
                            label_roots = TRUE,
                            show_trajectory_graph = TRUE) +
      scale_color_manual(values = CELLTYPE_L3_COLORS) +
      ggtitle(paste0(tissue_value, " - Trajectory Graph (Level 3)"))
    save_plot_pdf(p_traj_L3, file.path(fig_dir, "trajectory_graph_level3.pdf"))
    
    # Level 4 with trajectory graph
    p_traj_L4 <- plot_cells(cds,
                            color_cells_by = "cell_type_level_4",
                            label_cell_groups = FALSE,
                            label_leaves = FALSE,
                            label_branch_points = FALSE,
                            label_roots = FALSE,
                            show_trajectory_graph = TRUE) +
      scale_color_manual(values = CELLTYPE_L4_COLORS) +
      ggtitle(paste0(tissue_value, " - Trajectory Graph (Level 4)"))
    save_plot_pdf(p_traj_L4, file.path(fig_dir, "trajectory_graph_level4.pdf"))
  }
  
  # Export data
  pt <- data.frame(
    cell_id = colnames(cds),
    tissue = tissue_value,
    cluster = colData(cds)[[CLUSTER_COL]],
    cell_type_L3 = colData(cds)$cell_type_level_3,
    cell_type_L4 = colData(cds)$cell_type_level_4,
    pseudotime = pseudotime(cds),
    branch = colData(cds)$branch
  )
  write.csv(pt, file.path(out_dir, "pseudotime.csv"), row.names = FALSE)
  
  # Graph test
  if (RUN_GRAPH_TEST) {
    cat("  Running graph_test...\n")
    gt_mat <- tryCatch(
      monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4),
      error = function(e) NULL
    )
    if (!is.null(gt_mat)) {
      gt_mat <- gt_mat[order(gt_mat$q_value), , drop = FALSE]
      gt_df <- data.frame(gene_id = rownames(gt_mat), gt_mat, row.names = NULL)
      write.csv(gt_df, file.path(out_dir, "pseudotime_de_genes.csv"), row.names = FALSE)
    }
  }
  
  # Save CDS
  cds_dir <- file.path(out_dir, "cds")
  dir.create(cds_dir, showWarnings = FALSE)
  saveRDS(cds, file.path(cds_dir, "cds.rds"))
  
  cat("  ✓ Trajectory analysis complete for ", tissue_value, "\n\n", sep = "")
  
  return(cds)
}

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("LOADING DATA\n")
cat(strrep("=", 80), "\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("RDS not found: ", INPUT_RDS)
}

seurat_obj <- readRDS(INPUT_RDS)
cat("Loaded Seurat: cells=", ncol(seurat_obj), " genes=", nrow(seurat_obj), "\n", sep = "")

# ========== FILTER CLUSTER 16 (IMMUNE CONTAMINATION) ==========
if (CLUSTER_COL %in% colnames(seurat_obj@meta.data)) {
  cluster_16_cells <- sum(seurat_obj@meta.data[[CLUSTER_COL]] == "16", na.rm = TRUE)
  if (cluster_16_cells > 0) {
    cat("\n⚠️  Removing cluster 16 (n=", cluster_16_cells, " cells)\n", sep = "")
    cat("    Reason: Suspected immune contamination (IGHG1+/IGJ+/LTF+)\n")
    cat("    This cluster shows strong Ig signatures inconsistent with basal identity.\n\n")
    
    seurat_obj <- subset(seurat_obj, subset = !!sym(CLUSTER_COL) != "16")
    cat("After filtering: cells=", ncol(seurat_obj), "\n\n", sep = "")
  } else {
    cat("\nℹ️  Cluster 16 not found in data (already filtered or absent).\n\n")
  }
} else {
  warning("CLUSTER_COL '", CLUSTER_COL, "' not found in metadata. Skipping cluster 16 filter.")
}

# ========== ADD HIERARCHICAL ANNOTATIONS ==========
cat(strrep("=", 80), "\n")
cat("ADDING HIERARCHICAL CELL TYPE ANNOTATIONS\n")
cat(strrep("=", 80), "\n\n")

seurat_obj <- add_celltype_annotations(
  seurat_obj, 
  CLUSTER_COL, 
  CLUSTER_TO_CELLTYPE_L3, 
  CLUSTER_TO_CELLTYPE_L4
)

# Check annotations
cat("Level 3 cell types:\n")
print(table(seurat_obj@meta.data$cell_type_level_3))
cat("\nLevel 4 cell types:\n")
print(table(seurat_obj@meta.data$cell_type_level_4))

# ========== TISSUE-SPECIFIC ANALYSIS ==========
if (!is.null(TISSUE_KEY) && TISSUE_KEY %in% colnames(seurat_obj@meta.data)) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("MULTI-TISSUE TRAJECTORY ANALYSIS\n")
  cat(strrep("=", 80), "\n\n", sep = "")
  
  # Get top tissues
  tissue_counts <- sort(table(seurat_obj@meta.data[[TISSUE_KEY]]), decreasing = TRUE)
  tissues_to_analyze <- names(tissue_counts)[seq_len(min(4, length(tissue_counts)))]
  
  cat("Tissues to analyze:\n")
  for (tt in tissues_to_analyze) {
    cat("  ", tt, ": ", tissue_counts[[tt]], " cells\n", sep = "")
  }
  cat("\n")
  
  # Create skipped tissues log
  skipped_log_path <- file.path(OUTPUT_DIR, "skipped_tissues.txt")
  if (file.exists(skipped_log_path)) file.remove(skipped_log_path)
  write("# Skipped tissues log", file = skipped_log_path)
  write(paste0("# Generated: ", Sys.time()), file = skipped_log_path, append = TRUE)
  write("# Format: tissue | cells | reason", file = skipped_log_path, append = TRUE)
  write("", file = skipped_log_path, append = TRUE)
  
  cds_list <- list()
  for (tt in tissues_to_analyze) {
    cat(strrep("-", 80), "\n")
    cat("Processing tissue: ", tt, "\n", sep = "")
    out_dir <- file.path(OUTPUT_DIR, paste0("tissue_", sanitize_name(tt)))
    cds_list[[tt]] <- run_monocle3_for_tissue(seurat_obj, tt, out_dir, skipped_log = skipped_log_path)
  }
  
  # Summary of skipped tissues
  if (file.exists(skipped_log_path)) {
    skipped_lines <- readLines(skipped_log_path)
    skipped_count <- sum(!grepl("^#", skipped_lines) & nchar(skipped_lines) > 0)
    if (skipped_count > 0) {
      cat("\n⚠ ", skipped_count, " tissue(s) skipped. See: ", skipped_log_path, "\n", sep = "")
    } else {
      cat("\n✓ All tissues processed successfully.\n")
    }
  }
  
} else {
  # Single dataset analysis
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("SINGLE DATASET TRAJECTORY ANALYSIS\n")
  cat(strrep("=", 80), "\n\n", sep = "")
  
  # Create skipped log even for single dataset
  skipped_log_path <- file.path(OUTPUT_DIR, "skipped_tissues.txt")
  if (file.exists(skipped_log_path)) file.remove(skipped_log_path)
  write("# Single dataset mode - no tissues skipped", file = skipped_log_path)
  
  cds <- run_monocle3_for_tissue(seurat_obj, "all_cells", OUTPUT_DIR, skipped_log = skipped_log_path)
}

# ========== SUMMARY ==========
cat("\n", strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat("Output directory: ", OUTPUT_DIR, "\n", sep = "")
cat("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("=", 80), "\n")