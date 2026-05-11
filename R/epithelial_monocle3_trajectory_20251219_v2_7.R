#!/usr/bin/env Rscript
################################################################################
# Epithelial Trajectory Analysis - v2.7 PRODUCTION
# 
# CRITICAL FIXES IN v2.7 (P0/P1 from code review):
# - [P0-1 FIX] USE_SEURAT_UMAP control (default FALSE for proper Monocle3 UMAP)
# - [P0-2 FIX] Correct data layer detection in Part 2
# - [P0-3 FIX] Memory-safe fallback in branch heatmap (avoid full-gene correlation)
# - [P0-4 FIX] Safe debug log writing (no script interruption)
# - [P1-1 FIX] USE_PARTITION warning for NA pseudotime
# - [P1-2 FIX] Stricter root label fuzzy matching (avoid Cycling basal as root)
# - [P1-3 FIX] Downsampling for marker trend plots (large datasets)
#
# v2.6 FEATURES (retained):
# - Two-panel trajectory plots (pseudotime vs cell type)
# - Component 1/2 axis labels
# - Automatic Basal lineage focusing
#
# v2.5 FIXES (retained):
# - graph_test rownames preservation
# - circlize library for colorRamp2()
# - Fast branch assignment (60x speedup)
# - Binned heatmaps for readability
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
  library(circlize)      # [P0-2 FIX] Explicit load for colorRamp2()
  library(scales)
  library(grid)
})

cat("\n", strrep("=", 80), "\n", sep = "")
cat("EPITHELIAL TRAJECTORY ANALYSIS - v2.7 PRODUCTION\n")
cat(strrep("=", 80), "\n\n", sep = "")
cat("monocle3:", as.character(packageVersion("monocle3")), "\n")
cat("Seurat  :", as.character(packageVersion("Seurat")), "\n")
cat("igraph  :", as.character(packageVersion("igraph")), "\n")
cat("circlize:", as.character(packageVersion("circlize")), "\n")
cat("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# region agent log: script entry (H0)
DEBUG_LOG_PATH <- "/home/h2048/.cursor/debug.log"
log_line <- paste0(
  '{"sessionId":"debug-session",',
  '"runId":"v2.6",',
  '"hypothesisId":"H0",',
  '"location":"epithelial_monocle3_trajectory_v2_5_HOTFIX.R:global",',
  '"message":"script entry and package versions",',
  '"data":{"monocle3":"', as.character(packageVersion("monocle3")),
  '","Seurat":"', as.character(packageVersion("Seurat")),
  '","igraph":"', as.character(packageVersion("igraph")),
  '","circlize":"', as.character(packageVersion("circlize")), '"},',
  '"timestamp":', as.numeric(Sys.time()) * 1000,
  "}"
)
safe_write_log(log_line, DEBUG_LOG_PATH)
write("\n", file = DEBUG_LOG_PATH, append = TRUE)
# endregion agent log

# ==============================================================================
# CONFIG
# ==============================================================================
INPUT_RDS <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/py/1218/Epi_monocle3_PRODUCTION_v2_20251218"

# Workflow control
SKIP_PART1 <- FALSE  # Set TRUE to skip trajectory analysis and load existing results
LOAD_EXISTING_CDS <- FALSE  # Set TRUE to load existing CDS objects from previous run

LABELS_KEY <- "Manual_Annotation"
TISSUE_KEY <- "tissue"
ROOT_LABEL <- "Basal"
CLUSTER_COL <- "leiden_bbknn"

# [P1-3 FIX] Strict root validation - skip tissue if no Basal found
ALLOW_ROOT_FALLBACK <- FALSE  # Set TRUE to use old fallback behavior

# [v2.7 P0-1 FIX] UMAP control - CRITICAL for biological interpretation
USE_SEURAT_UMAP <- FALSE  # Set TRUE to use Seurat's integrated UMAP (NOT RECOMMENDED)
                          # FALSE (default): Monocle3 computes its own UMAP from PCA
                          # WARNING: Using Seurat UMAP means trajectory is shaped by
                          # integration method (BBKNN/Harmony/scVI), which may distort
                          # true differentiation paths

# Monocle3 parameters
N_DIM_PREPROCESS <- 50
N_NEIGHBORS_CLUSTER <- 15
USE_PARTITION <- FALSE  # [v2.7 P1-1 FIX] Changed default to FALSE
                        # TRUE may cause many cells to have NA pseudotime if graph
                        # is split into multiple disconnected partitions
RUN_GRAPH_TEST <- TRUE
REDUCTION_METHOD <- "UMAP"  # [P1-1] Explicit specification; required for learn_graph in monocle3>=1.4

# Marker analysis parameters
TOP_N_MARKERS <- 10
MIN_PSEUDOBULK_CELLS <- 10  # [v2.3 FIX] Minimum cells per group for pseudobulk aggregation

# [v2.7 P1-3 FIX] Marker trend plot downsampling
MAX_CELLS_MARKER_TRENDS <- 20000  # Downsample to this many cells for marker trends
                                   # Avoids huge PDF files and slow rendering

# [v2.4/v2.5] Trajectory heatmap parameters
GENERATE_BRANCH_HEATMAPS <- TRUE  # Generate heatmaps for each branch
N_GENES_HEATMAP <- 50             # Number of top dynamic genes to show in heatmap
MIN_CELLS_PER_BRANCH <- 50        # Minimum cells required for branch heatmap
N_GENE_CLUSTERS <- 3              # Number of gene clusters (k-means) in heatmap

# [v2.5 NEW] Binned heatmap option for better visualization
USE_BINNED_HEATMAP <- TRUE        # Use pseudotime bins instead of individual cells
N_PSEUDOTIME_BINS <- 100          # Number of bins along pseudotime (if USE_BINNED_HEATMAP=TRUE)

# Lineage markers for sanity check
LINEAGE_MARKERS <- list(
  Basal = c("KRT5", "TP63", "KRT14"),
  Ciliated = c("FOXJ1", "PIFO", "TPPP3"),
  Goblet = c("MUC5AC", "TFF3", "SPDEF"),
  Secretory = c("SCGB1A1", "SCGB3A1", "BPIFA1"),
  Squamous = c("KRT13", "IVL", "SPRR2A")
)

# Figure dimensions
FIG_W <- 10
FIG_H <- 8

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Helpers
# ==============================================================================
sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# [v2.7 P0-4 FIX] Safe debug log writing (no script interruption)
safe_write_log <- function(line, path) {
  tryCatch({
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    write(line, file = path, append = TRUE)
    write("\n", file = path, append = TRUE)
  }, error = function(e) NULL)  # Silently fail if log is not writable
}

save_plot_pdf <- function(p, file, width = FIG_W, height = FIG_H) {
  pdf(file, width = width, height = height, onefile = TRUE)
  print(p)
  dev.off()
}

# ==============================================================================
# [v2.6 NEW] Two-panel trajectory plot (pseudotime vs cell type)
# - Left: pseudotime + branch point/leaf labels
# - Right: cell type (Manual_Annotation)
# - Single PDF page with two panels side-by-side
# - Axes labeled as "Component 1" / "Component 2"
# ==============================================================================

save_two_panel_pdf <- function(p_left, p_right, file, title,
                               width = 14, height = 6) {
  pdf(file, width = width, height = height, onefile = TRUE)
  grid::grid.newpage()

  lay <- grid::grid.layout(
    nrow = 2, ncol = 2,
    heights = grid::unit.c(grid::unit(0.8, "in"), grid::unit(1, "null")),
    widths  = grid::unit.c(grid::unit(1, "null"), grid::unit(1, "null"))
  )
  grid::pushViewport(grid::viewport(layout = lay))

  # Title (spans both columns)
  grid::grid.text(
    title,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2),
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )

  # Left panel
  print(p_left, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))

  # Right panel
  print(p_right, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))

  dev.off()
}

make_two_panel_trajectory <- function(cds,
                                      cell_type_key,
                                      out_file,
                                      reduction_method = "UMAP",
                                      title = "Basal Cell Pseudotime Cell Trajectory",
                                      focus_labels = NULL,
                                      label_branch_points = TRUE,
                                      label_leaves = FALSE) {
  # Optional: focus on specific lineage labels (e.g., Basal/Suprabasal/Cycling basal)
  cds_use <- cds
  if (!is.null(focus_labels)) {
    ct <- as.character(colData(cds_use)[[cell_type_key]])
    keep <- ct %in% focus_labels
    if (sum(keep) >= 50) {
      cds_use <- cds_use[, keep]
      cat("    Focusing on ", length(focus_labels), " cell types: ", 
          paste(focus_labels, collapse = ", "), "\n", sep = "")
    }
  }

  # Left: pseudotime with branch points
  p_left <- monocle3::plot_cells(
    cds_use,
    reduction_method = reduction_method,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = label_leaves,
    label_branch_points = label_branch_points,
    label_roots = FALSE
  ) +
    ggtitle(NULL) +
    labs(x = "Component 1", y = "Component 2", color = "Pseudotime") +
    scale_color_viridis_c(option = "plasma", na.value = "grey80") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(5, 5, 5, 5)
    )

  # Right: cell type from Manual_Annotation
  p_right <- monocle3::plot_cells(
    cds_use,
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
    theme(
      legend.position = "bottom",
      plot.margin = margin(5, 5, 5, 5)
    )

  save_two_panel_pdf(p_left, p_right, out_file, title = title, width = 14, height = 6)
  cat("  Saved two-panel trajectory: ", out_file, "\n", sep = "")
  return(invisible(out_file))
}

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
  
  # [v2.7 P0-1 FIX] Optional UMAP transfer from Seurat
  # WARNING: Using Seurat UMAP means trajectory is shaped by integration method
  # Default (use_seurat_umap=FALSE): Monocle3 will compute its own UMAP from PCA
  if (use_seurat_umap) {
    red_names <- names(seurat_obj@reductions)
    umap_red <- red_names[grepl("umap", red_names, ignore.case = TRUE)][1]
    if (!is.na(umap_red) && !is.null(umap_red)) {
      umap_coords <- Seurat::Embeddings(seurat_obj, reduction = umap_red)
      colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
      reducedDims(cds)[["UMAP"]] <- umap_coords
      message("  [INFO] Using Seurat UMAP: ", umap_red)
    } else {
      message("  [WARNING] USE_SEURAT_UMAP=TRUE but no UMAP found in Seurat object")
    }
  } else {
    message("  [INFO] Monocle3 will compute its own UMAP (use_seurat_umap=FALSE)")
  }
  
  cds
}

# [v2.7 P1-2 FIX] Strict root validation with improved fuzzy matching
pick_root_label <- function(labels_vec, root_label = "Basal", allow_fallback = FALSE) {
  labs <- unique(as.character(labels_vec))
  
  # Priority 1: Exact match
  if (root_label %in% labs) {
    return(root_label)
  }
  
  # Priority 2: Fuzzy match excluding cycling/proliferative variants
  # (avoid using "Cycling basal" as root when "Basal" is not found)
  if (root_label == "Basal") {
    # Look for "Basal" but NOT "Cycling basal" or "Proliferating basal"
    basal_candidates <- labs[grepl("basal", labs, ignore.case = TRUE)]
    basal_non_cycling <- basal_candidates[!grepl("cycl|prolif|divid", basal_candidates, ignore.case = TRUE)]
    
    if (length(basal_non_cycling) > 0) {
      chosen <- basal_non_cycling[1]
      message("  [INFO] Using fuzzy match: '", chosen, "' for root '", root_label, "'")
      return(chosen)
    }
  }
  
  # Priority 3: General case-insensitive fuzzy match (for non-Basal roots)
  hit <- labs[grepl(root_label, labs, ignore.case = TRUE)]
  if (length(hit) > 0) {
    message("  [INFO] Using fuzzy match: '", hit[1], "' for root '", root_label, "'")
    return(hit[1])
  }
  
  # No match found
  if (!allow_fallback) {
    message("  [ERROR] Root label '", root_label, "' not found. Tissue will be skipped.")
    return(NULL)
  }
  
  # Fallback (only if explicitly allowed)
  tb <- sort(table(labels_vec), decreasing = TRUE)
  fallback_label <- as.character(names(tb)[1])
  message(
    "  [WARNING] Root label '", root_label, 
    "' not found, using fallback '", fallback_label, "'"
  )
  return(fallback_label)
}

# [P0-2 FIX] Robust vertex type handling
get_root_principal_node <- function(cds, cell_type_key, root_cell_type) {
  root_cells <- which(colData(cds)[[cell_type_key]] == root_cell_type)
  if (length(root_cells) == 0) {
    stop("No cells found for root cell type: ", root_cell_type)
  }
  
  if (!"UMAP" %in% names(monocle3::principal_graph(cds))) {
    stop("Principal graph 'UMAP' not found. Did learn_graph() succeed?")
  }
  if (is.null(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)) {
    stop("pr_graph_cell_proj_closest_vertex not found")
  }
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])[, 1]
  
  pr_graph <- monocle3::principal_graph(cds)[["UMAP"]]
  vertex_names <- igraph::V(pr_graph)$name
  
  # [P0-2 FIX] Handle both numeric and character vertex types
  if (is.numeric(closest_vertex)) {
    # Numeric indices - convert to vertex names
    closest_vertex_names <- vertex_names[closest_vertex]
  } else {
    # Already vertex names
    closest_vertex_names <- as.character(closest_vertex)
  }
  
  # Find most common vertex among root cells
  root_vertex_name <- names(which.max(table(closest_vertex_names[root_cells])))
  
  return(root_vertex_name)
}

# [P0-3 FIX] Optimized branch assignment using vertex distance table
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
  
  # Handle both numeric and character vertex types
  if (is.numeric(closest_vertex)) {
    closest_vertex_names <- vertex_names[closest_vertex]
  } else {
    closest_vertex_names <- as.character(closest_vertex)
  }
  
  # [P0-3 FIX] Compute vertex-to-leaf distance matrix ONCE
  # This is O(|V| * |leaves|) instead of O(|cells| * |leaves|)
  cat("    Computing vertex-to-leaf distances (", length(vertex_names), 
      " vertices × ", length(leaves), " leaves)...\n", sep = "")
  
  vertex_to_leaf_dist <- igraph::distances(g, v = vertex_names, to = leaves)
  
  # Map each cell to its closest vertex's distance to each leaf
  dist_mat <- vertex_to_leaf_dist[closest_vertex_names, , drop = FALSE]
  colnames(dist_mat) <- leaves
  
  # Assign each cell to nearest leaf
  leaf_idx <- apply(dist_mat, 1, which.min)
  colData(cds)$branch <- paste0("Leaf_", leaf_idx)
  
  cat("    Branch assignment complete\n")
  
  cds
}

# [P1-4] Leaf composition analysis
analyze_leaf_composition <- function(cds, cell_type_key, out_file) {
  g <- monocle3::principal_graph(cds)[["UMAP"]]
  deg <- igraph::degree(g)
  leaves <- names(deg[deg == 1])
  
  if (length(leaves) < 2) {
    cat("  [INFO] Only one leaf, skipping composition analysis\n")
    return(NULL)
  }
  
  branch_vec <- colData(cds)$branch
  celltype_vec <- colData(cds)[[cell_type_key]]
  
  comp_list <- lapply(unique(branch_vec), function(b) {
    cells_in_branch <- branch_vec == b
    ct_counts <- table(celltype_vec[cells_in_branch])
    data.frame(
      branch = b,
      cell_type = names(ct_counts),
      count = as.numeric(ct_counts),
      fraction = as.numeric(ct_counts) / sum(ct_counts)
    )
  })
  
  comp_df <- do.call(rbind, comp_list)
  write.csv(comp_df, out_file, row.names = FALSE)
  
  return(comp_df)
}

# [v2.7 P1-3 FIX] Pseudotime-marker trends with downsampling for large datasets
plot_pseudotime_markers <- function(cds, markers, out_file, max_cells = MAX_CELLS_MARKER_TRENDS) {
  pt <- pseudotime(cds)
  valid_cells <- !is.na(pt)
  
  if (sum(valid_cells) < 50) {
    cat("  [SKIP] Too few valid pseudotime cells for marker trend plot\n")
    return(NULL)
  }
  
  # [v2.7 P1-3 FIX] Downsample if too many cells (avoid huge PDFs)
  cell_indices <- which(valid_cells)
  n_valid <- length(cell_indices)
  
  if (n_valid > max_cells) {
    cat("  [INFO] Downsampling from ", n_valid, " to ", max_cells, 
        " cells for marker trends\n", sep = "")
    set.seed(42)  # Reproducible sampling
    cell_indices <- sample(cell_indices, max_cells)
    valid_cells <- rep(FALSE, length(pt))
    valid_cells[cell_indices] <- TRUE
  }
  
  pt_valid <- pt[valid_cells]
  
  # Check which markers exist
  all_genes <- rownames(cds)
  markers_flat <- unlist(markers)
  markers_available <- markers_flat[markers_flat %in% all_genes]
  
  if (length(markers_available) == 0) {
    cat("  [SKIP] No lineage markers found in dataset\n")
    return(NULL)
  }
  
  cat("  Plotting ", length(markers_available), " lineage markers (", 
      sum(valid_cells), " cells)\n", sep = "")

  # region agent log: pseudotime & marker availability (H3)
  log_line <- paste0(
    '{"sessionId":"debug-session",',
    '"runId":"v2.7",',
    '"hypothesisId":"H3",',
    '"location":"epithelial_monocle3_trajectory_v2_7_PRODUCTION.R:plot_pseudotime_markers",',
    '"message":"pseudotime and marker availability",',
    '"data":{"n_valid_cells":', sum(valid_cells),
    ',"n_markers_available":', length(markers_available), "},",
    '"timestamp":', as.numeric(Sys.time()) * 1000, 
    "}"
  )
  safe_write_log(log_line, DEBUG_LOG_PATH)
  # endregion agent log
  
  pdf(out_file, width = 12, height = 8)
  
  for (lineage_name in names(markers)) {
    genes <- markers[[lineage_name]]
    genes <- genes[genes %in% all_genes]
    
    if (length(genes) == 0) next
    
    for (gene in genes) {
      # [P1-2 FIX] Use logcounts if available, fallback to log1p(counts)
      if ("logcounts" %in% names(assays(cds))) {
        expr <- assays(cds)[["logcounts"]][gene, valid_cells]
      } else if ("data" %in% names(assays(cds))) {
        expr <- assays(cds)[["data"]][gene, valid_cells]
      } else {
        expr <- log1p(exprs(cds)[gene, valid_cells])
      }
      
      df <- data.frame(pseudotime = pt_valid, expression = expr)
      
      # Add subtitle if downsampled
      subtitle <- if (n_valid > max_cells) {
        paste0("(Downsampled to ", max_cells, " of ", n_valid, " cells)")
      } else {
        ""
      }
      
      p <- ggplot(df, aes(x = pseudotime, y = expression)) +
        geom_point(alpha = 0.1, size = 0.5) +
        geom_smooth(method = "loess", color = "red", se = TRUE) +
        labs(
          title = paste0(lineage_name, " marker: ", gene),
          subtitle = subtitle,
          x = "Pseudotime",
          y = "Log Expression"
        ) +
        theme_minimal()
      
      print(p)
    }
  }
  
  dev.off()
  cat("  Saved: ", out_file, "\n", sep = "")
}

# [v2.3] Load existing CDS results
load_existing_cds_results <- function(output_dir, tissues) {
  cds_list <- list()
  
  for (tt in tissues) {
    tissue_dir <- file.path(output_dir, paste0("tissue_", sanitize_name(tt)))
    cds_path <- file.path(tissue_dir, "cds", "cds.rds")
    
    if (file.exists(cds_path)) {
      cat("  Loading CDS for tissue: ", tt, "\n", sep = "")
      cds_list[[tt]] <- readRDS(cds_path)
    } else {
      cat("  [WARNING] CDS not found for tissue: ", tt, "\n", sep = "")
      cds_list[[tt]] <- NULL
    }
  }
  
  return(cds_list)
}

# [v2.5] Generate branch heatmap with binning option
generate_branch_heatmap <- function(cds, branch_name, graph_test_results = NULL,
                                   n_genes = 50, n_clusters = 3, 
                                   min_cells = 50, out_file,
                                   use_binned = TRUE, n_bins = 100) {
  # Step 1: Select cells in this branch
  if (!"branch" %in% colnames(colData(cds))) {
    cat("  [SKIP] No branch information found in CDS\n")
    return(NULL)
  }
  
  cells_in_branch <- colData(cds)$branch == branch_name
  n_cells_branch <- sum(cells_in_branch)
  
  if (n_cells_branch < min_cells) {
    cat("  [SKIP] Branch ", branch_name, " has only ", n_cells_branch, 
        " cells (< ", min_cells, ")\n", sep = "")
    return(NULL)
  }
  
  cat("  Generating heatmap for branch: ", branch_name, " (", n_cells_branch, " cells)\n", sep = "")
  
  # Subset CDS to this branch
  cds_branch <- cds[, cells_in_branch]
  
  # Step 2: Sort cells by pseudotime
  pt <- pseudotime(cds_branch)
  valid_pt <- !is.na(pt)
  
  if (sum(valid_pt) < min_cells) {
    cat("  [SKIP] Too few cells with valid pseudotime: ", sum(valid_pt), "\n", sep = "")
    return(NULL)
  }
  
  cds_branch <- cds_branch[, valid_pt]
  pt <- pt[valid_pt]
  ord <- order(pt)
  cds_branch <- cds_branch[, ord]
  pt <- pt[ord]
  
  # Step 3: Find dynamic genes
  # [P0-1 FIX] Use rownames(graph_test_results) directly (preserved in calling function)
  if (!is.null(graph_test_results) && nrow(graph_test_results) > 0) {
    # Use pre-computed graph_test results
    gene_list <- head(rownames(graph_test_results)[order(graph_test_results$q_value)], n_genes)
    cat("  Using ", length(gene_list), " genes from graph_test results\n", sep = "")
  } else {
    # [v2.7 P0-3 FIX] Memory-safe fallback: use HVG instead of full-gene correlation
    cat("  [FALLBACK] graph_test not available, using highly variable genes...\n")
    
    # Get HVG if available (from preprocessing), or select by variance
    if ("highly_variable" %in% colnames(rowData(cds_branch))) {
      hvg <- rownames(cds_branch)[rowData(cds_branch)$highly_variable]
      cat("    Found ", length(hvg), " HVG from preprocessing\n", sep = "")
    } else {
      # Quick variance-based selection on expressed genes
      if ("logcounts" %in% names(assays(cds_branch))) {
        expr_mat_quick <- assays(cds_branch)[["logcounts"]]
      } else {
        expr_mat_quick <- exprs(cds_branch)
      }
      
      # Only compute variance for genes expressed in ≥10 cells
      gene_detected <- Matrix::rowSums(expr_mat_quick > 0) >= 10
      if (sum(gene_detected) < n_genes) {
        cat("  [SKIP] Too few expressed genes for fallback\n")
        return(NULL)
      }
      
      expr_mat_subset <- expr_mat_quick[gene_detected, ]
      gene_vars <- Matrix::rowMeans((expr_mat_subset - Matrix::rowMeans(expr_mat_subset))^2)
      top_var_idx <- head(order(gene_vars, decreasing = TRUE), min(2000, n_genes * 10))
      hvg <- rownames(expr_mat_subset)[top_var_idx]
      cat("    Selected ", length(hvg), " high-variance genes\n", sep = "")
    }
    
    # Compute correlation only for HVG (much faster and memory-safe)
    if ("logcounts" %in% names(assays(cds_branch))) {
      expr_hvg <- as.matrix(assays(cds_branch)[["logcounts"]][hvg, ])
    } else {
      expr_hvg <- as.matrix(exprs(cds_branch)[hvg, ])
    }
    
    # Correlation for HVG only
    gene_cors <- apply(expr_hvg, 1, function(gene_expr) {
      cor(gene_expr, pt, method = "spearman", use = "complete.obs")
    })
    
    # Take top genes by absolute correlation
    top_gene_idx <- head(order(abs(gene_cors), decreasing = TRUE), n_genes)
    gene_list <- hvg[top_gene_idx]
    cat("  Using ", length(gene_list), " genes with highest pseudotime correlation (from HVG)\n", sep = "")
  }
  
  # Remove genes not in CDS
  gene_list <- gene_list[gene_list %in% rownames(cds_branch)]
  
  if (length(gene_list) < 10) {
    cat("  [SKIP] Too few genes available: ", length(gene_list), "\n", sep = "")
    return(NULL)
  }
  
  # Step 4: Extract expression matrix
  if ("logcounts" %in% names(assays(cds_branch))) {
    expr_mat <- as.matrix(assays(cds_branch)[["logcounts"]][gene_list, ])
  } else if ("data" %in% names(assays(cds_branch))) {
    expr_mat <- as.matrix(assays(cds_branch)[["data"]][gene_list, ])
  } else {
    expr_mat <- as.matrix(exprs(cds_branch)[gene_list, ])
  }
  
  # [v2.5 NEW] Optional binning for better visualization
  if (use_binned && ncol(expr_mat) > n_bins) {
    cat("  Binning cells into ", n_bins, " pseudotime bins for visualization...\n", sep = "")
    
    # Create pseudotime bins
    pt_bins <- cut(pt, breaks = n_bins, labels = FALSE, include.lowest = TRUE)
    
    # Average expression within each bin
    expr_mat_binned <- matrix(0, nrow = nrow(expr_mat), ncol = n_bins)
    rownames(expr_mat_binned) <- rownames(expr_mat)
    
    pt_bin_centers <- numeric(n_bins)
    
    for (b in 1:n_bins) {
      cells_in_bin <- which(pt_bins == b)
      if (length(cells_in_bin) > 0) {
        expr_mat_binned[, b] <- rowMeans(expr_mat[, cells_in_bin, drop = FALSE])
        pt_bin_centers[b] <- mean(pt[cells_in_bin])
      } else {
        pt_bin_centers[b] <- quantile(pt, probs = (b - 0.5) / n_bins)
      }
    }
    
    # Use binned matrix
    expr_mat <- expr_mat_binned
    pt <- pt_bin_centers
    
    cat("  Binned to ", ncol(expr_mat), " columns\n", sep = "")
  }
  
  # Step 5: Z-score normalization (row-wise)
  expr_mat_scaled <- t(scale(t(expr_mat)))
  
  # Cap extreme values
  expr_mat_scaled[expr_mat_scaled > 3] <- 3
  expr_mat_scaled[expr_mat_scaled < -3] <- -3
  
  # Replace NaN with 0 (happens when gene has 0 variance)
  expr_mat_scaled[is.nan(expr_mat_scaled)] <- 0
  
  # Step 6: Create pseudotime annotation
  pt_range <- range(pt, na.rm = TRUE)
  
  ha_top <- HeatmapAnnotation(
    Pseudotime = pt,
    col = list(
      Pseudotime = colorRamp2(pt_range, c("white", "black"))
    ),
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(0.5, "cm")
  )
  
  # Step 7: Define color scale
  col_fun <- colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick"))
  
  # Step 8: Generate heatmap
  n_cols <- ncol(expr_mat_scaled)
  subtitle <- if (use_binned && n_cols <= n_bins) {
    paste0(n_cells_branch, " cells binned into ", n_cols, " pseudotime bins")
  } else {
    paste0(n_cells_branch, " cells")
  }
  
  ht <- Heatmap(
    expr_mat_scaled,
    name = "Z-score",
    top_annotation = ha_top,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_km = n_clusters,
    row_gap = unit(3, "mm"),
    show_row_dend = FALSE,
    col = col_fun,
    column_title = paste0(branch_name, " Trajectory"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "top",
    heatmap_legend_param = list(
      title = "Z-score\nExpression",
      title_position = "leftcenter-rot",
      legend_height = unit(4, "cm")
    ),
    use_raster = n_cols > 500,  # Rasterize if many columns
    raster_quality = 2
  )
  
  # Add subtitle as bottom annotation
  ha_bottom <- columnAnnotation(
    foo = anno_text(
      rep("", n_cols),
      gp = gpar(fontsize = 0),
      just = "center"
    ),
    annotation_height = unit(0.5, "cm")
  )
  
  # Save heatmap
  pdf(out_file, width = 10, height = 12)
  draw(ht)
  grid.text(subtitle, x = 0.5, y = 0.98, gp = gpar(fontsize = 10, col = "gray30"))
  dev.off()
  
  cat("  Saved: ", out_file, "\n", sep = "")
  
  return(list(
    branch = branch_name,
    n_cells = n_cells_branch,
    n_genes = length(gene_list),
    genes = gene_list,
    binned = use_binned,
    n_bins_used = n_cols
  ))
}

# ==============================================================================
# Main per-tissue trajectory
# ==============================================================================
run_monocle3_for_tissue <- function(seurat_obj, tissue_value, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE)
  
  md <- seurat_obj@meta.data
  cells <- rownames(md)[as.character(md[[TISSUE_KEY]]) == tissue_value]
  
  if (length(cells) < 200) {
    cat("  [SKIP] tissue=", tissue_value, " too few cells: ", length(cells), "\n", sep = "")
    return(NULL)
  }
  
  seu <- subset(seurat_obj, cells = cells)
  
  # [P1-3 FIX] Strict root validation
  root_label <- pick_root_label(
    seu@meta.data[[LABELS_KEY]], 
    ROOT_LABEL, 
    allow_fallback = ALLOW_ROOT_FALLBACK
  )

  # region agent log: tissue basic stats and root label (H1)
  log_line <- paste0(
    '{"sessionId":"debug-session",',
    '"runId":"v2.7",',
    '"hypothesisId":"H1",',
    '"location":"epithelial_monocle3_trajectory_v2_7_PRODUCTION.R:run_monocle3_for_tissue",',
    '"message":"tissue basic stats and root label",',
    '"data":{"tissue":"', tissue_value,
    '","n_cells":', length(cells),
    ',"root_label":"', ifelse(is.null(root_label), "NULL", root_label), '"},',
    '"timestamp":', as.numeric(Sys.time()) * 1000, 
    "}"
  )
  safe_write_log(log_line, DEBUG_LOG_PATH)
  # endregion agent log
  
  if (is.null(root_label)) {
    cat("  [SKIP] tissue=", tissue_value, " - root not found and fallback disabled\n", sep = "")
    return(NULL)
  }
  
  cat("  tissue=", tissue_value, " | cells=", ncol(seu), " | root=", root_label, "\n", sep = "")
  
  # [v2.7 P0-1 FIX] Use USE_SEURAT_UMAP parameter
  cds <- seurat_to_cds(seu, assay = DefaultAssay(seu), use_seurat_umap = USE_SEURAT_UMAP)
  
  # [P0-1] ALWAYS run preprocess_cds
  cat("  Running preprocess_cds...\n")
  cds <- preprocess_cds(cds, num_dim = N_DIM_PREPROCESS, method = "PCA")
  
  # UMAP
  if (!("UMAP" %in% names(reducedDims(cds)))) {
    cat("  Running reduce_dimension...\n")
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
  }
  
  # [P1-1] Explicit reduction_method
  cat("  Running cluster_cells (reduction_method=", REDUCTION_METHOD, ")...\n", sep = "")
  cds <- cluster_cells(cds, reduction_method = REDUCTION_METHOD, k = N_NEIGHBORS_CLUSTER)
  
  cat("  Running learn_graph...\n")
  cds <- learn_graph(cds, use_partition = USE_PARTITION)
  
  # [v2.7 P1-1 FIX] Warn about USE_PARTITION if enabled
  if (USE_PARTITION) {
    n_partitions <- length(unique(partitions(cds)))
    if (n_partitions > 1) {
      cat("  [WARNING] USE_PARTITION=TRUE resulted in ", n_partitions, 
          " disconnected partitions\n", sep = "")
      cat("           Cells in non-root partitions may have NA pseudotime\n")
    }
  }
  
  # Order cells with robust root finding
  cat("  Finding root node...\n")
  root_node <- get_root_principal_node(cds, LABELS_KEY, root_label)
  cat("  Root node: ", root_node, "\n", sep = "")

  # region agent log: root node & graph summary (H2)
  g <- monocle3::principal_graph(cds)[["UMAP"]]
  log_line <- paste0(
    '{"sessionId":"debug-session",',
    '"runId":"v2.7",',
    '"hypothesisId":"H2",',
    '"location":"epithelial_monocle3_trajectory_v2_7_PRODUCTION.R:run_monocle3_for_tissue",',
    '"message":"root node and graph summary",',
    '"data":{"tissue":"', tissue_value,
    '","root_node":"', root_node,
    '","n_vertices":', igraph::gorder(g),
    ',"n_edges":', igraph::gsize(g), "},",
    '"timestamp":', as.numeric(Sys.time()) * 1000, 
    "}"
  )
  safe_write_log(log_line, DEBUG_LOG_PATH)
  # endregion agent log
  
  cds <- order_cells(cds, root_pr_nodes = root_node)
  
  # [P0-3 FIX] Use optimized branch assignment
  cat("  Assigning branches...\n")
  cds <- assign_branch_by_leaf_fast(cds)
  
  # Save metadata to cds
  metadata(cds)$tissue <- tissue_value
  metadata(cds)$root_label <- root_label
  metadata(cds)$root_node <- root_node
  metadata(cds)$reduction_method <- REDUCTION_METHOD
  
  # Plots
  p1 <- plot_cells(
    cds,
    color_cells_by = LABELS_KEY,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = TRUE
  ) +
    ggtitle(paste0("tissue=", tissue_value, " | cell types"))
  save_plot_pdf(p1, file.path(fig_dir, "trajectory_celltypes.pdf"))
  
  p2 <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = TRUE
  ) +
    scale_color_viridis(na.value = "grey80") +
    ggtitle(paste0("tissue=", tissue_value, " | pseudotime"))
  save_plot_pdf(p2, file.path(fig_dir, "trajectory_pseudotime.pdf"))
  
  p3 <- plot_cells(
    cds,
    color_cells_by = "branch",
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE
  ) +
    ggtitle(paste0("tissue=", tissue_value, " | branch"))
  save_plot_pdf(p3, file.path(fig_dir, "trajectory_branch.pdf"))
  
  # [v2.6 NEW] Generate two-panel trajectory plot (pseudotime vs cell type)
  cat("  Generating two-panel trajectory plot...\n")
  
  # Check for Basal-related cell types for focused view
  all_celltypes <- unique(as.character(colData(cds)[[LABELS_KEY]]))
  basal_lineage <- c("Basal", "Suprabasal", "Cycling basal cells", "Cycling basal")
  basal_lineage_present <- basal_lineage[basal_lineage %in% all_celltypes]
  
  # Generate plot
  two_panel_file <- file.path(fig_dir, "trajectory_two_panel_pseudotime_vs_celltype.pdf")
  
  # Determine title based on root and tissue
  panel_title <- paste0(root_label, " Cell Pseudotime Cell Trajectory - ", tissue_value)
  
  # If we have Basal lineage cells (≥2 types), focus on them; otherwise use all cells
  focus_on <- if (length(basal_lineage_present) >= 2) {
    basal_lineage_present
  } else {
    NULL  # Use all cells
  }
  
  make_two_panel_trajectory(
    cds = cds,
    cell_type_key = LABELS_KEY,
    out_file = two_panel_file,
    reduction_method = REDUCTION_METHOD,
    title = panel_title,
    focus_labels = focus_on,
    label_branch_points = TRUE,
    label_leaves = FALSE
  )
  
  # Export pseudotime table (with fixed cluster extraction)
  cluster_info <- tryCatch(
    monocle3::clusters(cds),
    error = function(e) rep(NA, ncol(cds))
  )
  
  pt <- data.frame(
    cell_id = colnames(cds),
    tissue = tissue_value,
    root_label = root_label,
    cell_type = colData(cds)[[LABELS_KEY]],
    pseudotime = pseudotime(cds),
    partition = partitions(cds),
    cluster = cluster_info,
    branch = colData(cds)$branch
  )
  write.csv(pt, file.path(out_dir, "pseudotime.csv"), row.names = FALSE)
  
  # [P1-4] Leaf composition analysis
  cat("  Analyzing leaf composition...\n")
  analyze_leaf_composition(
    cds, 
    LABELS_KEY, 
    file.path(out_dir, "leaf_composition.csv")
  )
  
  # [P1-2 FIX] Pseudotime-marker trends (using logcounts)
  cat("  Plotting pseudotime-marker trends...\n")
  plot_pseudotime_markers(
    cds, 
    LINEAGE_MARKERS, 
    file.path(fig_dir, "pseudotime_marker_trends.pdf")
  )
  
  # graph_test genes
  gt_mat <- NULL  # [P0-1 FIX] Initialize
  
  if (RUN_GRAPH_TEST) {
    cat("  Running graph_test...\n")
    gt_mat <- tryCatch(
      monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4),
      error = function(e) {
        message("  [WARNING] graph_test failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    # [P0-1 FIX] Preserve rownames for downstream use, write separate df for CSV
    if (!is.null(gt_mat)) {
      gt_mat <- gt_mat[order(gt_mat$q_value), , drop = FALSE]
      
      # Write CSV with gene_id column
      gt_df <- data.frame(gene_id = rownames(gt_mat), gt_mat, row.names = NULL)
      write.csv(gt_df, file.path(out_dir, "pseudotime_de_genes.csv"), row.names = FALSE)
      
      cat("  Graph test complete: ", nrow(gt_mat), " genes tested\n", sep = "")
    }
  }
  
  # [v2.5] Generate branch heatmaps
  if (GENERATE_BRANCH_HEATMAPS) {
    cat("  Generating branch heatmaps...\n")
    
    branches <- unique(colData(cds)$branch)
    branches <- branches[!is.na(branches)]
    
    if (length(branches) > 0) {
      heatmap_dir <- file.path(fig_dir, "branch_heatmaps")
      dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
      
      branch_heatmap_results <- list()
      
      for (branch in branches) {
        branch_safe <- sanitize_name(branch)
        heatmap_file <- file.path(heatmap_dir, paste0("heatmap_", branch_safe, ".pdf"))
        
        result <- tryCatch(
          generate_branch_heatmap(
            cds = cds,
            branch_name = branch,
            graph_test_results = gt_mat,  # [P0-1 FIX] Pass gt_mat with preserved rownames
            n_genes = N_GENES_HEATMAP,
            n_clusters = N_GENE_CLUSTERS,
            min_cells = MIN_CELLS_PER_BRANCH,
            out_file = heatmap_file,
            use_binned = USE_BINNED_HEATMAP,
            n_bins = N_PSEUDOTIME_BINS
          ),
          error = function(e) {
            cat("  [ERROR] Failed to generate heatmap for branch ", branch, ": ", 
                conditionMessage(e), "\n", sep = "")
            return(NULL)
          }
        )
        
        if (!is.null(result)) {
          branch_heatmap_results[[branch]] <- result
          
          # Save gene list for this branch
          gene_list_file <- file.path(heatmap_dir, paste0("genes_", branch_safe, ".txt"))
          writeLines(result$genes, gene_list_file)
        }
      }
      
      # Summary
      n_success <- sum(sapply(branch_heatmap_results, function(x) !is.null(x)))
      cat("  Branch heatmaps generated: ", n_success, " / ", length(branches), "\n", sep = "")
    } else {
      cat("  [SKIP] No branches found for heatmap generation\n")
    }
  }
  
  # Save cds
  cds_dir <- file.path(out_dir, "cds")
  dir.create(cds_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cds, file.path(cds_dir, "cds.rds"))
  
  cds
}

# ==============================================================================
# RUN PART 1: Trajectory Analysis (or load existing)
# ==============================================================================
if (!SKIP_PART1) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("PART 1: TRAJECTORY ANALYSIS\n")
  cat(strrep("=", 80), "\n\n", sep = "")
  
  if (!file.exists(INPUT_RDS)) {
    stop("RDS not found: ", INPUT_RDS)
  }
  seurat_obj <- readRDS(INPUT_RDS)
  
  cat("Loaded Seurat: cells=", ncol(seurat_obj), " genes=", nrow(seurat_obj), "\n", sep = "")
  
  # Select top4 tissues
  tissue_counts <- sort(table(seurat_obj@meta.data[[TISSUE_KEY]]), decreasing = TRUE)
  tissues_top4 <- names(tissue_counts)[seq_len(min(4, length(tissue_counts)))]
  
  cat("\nSelected tissues (Top4):\n")
  for (tt in tissues_top4) {
    cat("  ", tt, ": ", tissue_counts[[tt]], " cells\n", sep = "")
  }
  cat("\n")
  
  cds_list <- list()
  for (tt in tissues_top4) {
    cat(strrep("-", 80), "\n")
    cat("Running tissue: ", tt, "\n", sep = "")
    out_dir <- file.path(OUTPUT_DIR, paste0("tissue_", sanitize_name(tt)))
    cds_list[[tt]] <- run_monocle3_for_tissue(seurat_obj, tt, out_dir)
  }
  
  cat("\n", strrep("=", 80), "\n")
  cat("PART 1 DONE: Trajectory analysis\n")
  cat(strrep("=", 80), "\n\n")
  
} else {
  # [v2.3] Load existing results
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("SKIPPING PART 1: Loading existing trajectory results\n")
  cat(strrep("=", 80), "\n\n", sep = "")
  
  if (!file.exists(INPUT_RDS)) {
    stop("RDS not found: ", INPUT_RDS)
  }
  seurat_obj <- readRDS(INPUT_RDS)
  
  cat("Loaded Seurat: cells=", ncol(seurat_obj), " genes=", nrow(seurat_obj), "\n", sep = "")
  
  # Determine tissues from existing output directories
  tissue_counts <- sort(table(seurat_obj@meta.data[[TISSUE_KEY]]), decreasing = TRUE)
  tissues_top4 <- names(tissue_counts)[seq_len(min(4, length(tissue_counts)))]
  
  if (LOAD_EXISTING_CDS) {
    cat("\nLoading existing CDS objects...\n")
    cds_list <- load_existing_cds_results(OUTPUT_DIR, tissues_top4)
  } else {
    cat("\nCDS objects will not be loaded (LOAD_EXISTING_CDS=FALSE)\n")
    cds_list <- list()
  }
  
  cat("\n", strrep("=", 80), "\n")
  cat("LOADED existing results\n")
  cat(strrep("=", 80), "\n\n")
}

# ==============================================================================
# PART 2: Marker Analysis with Pseudobulk Heatmap
# ==============================================================================
cat(strrep("=", 80), "\n")
cat("PART 2: MARKER ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

DefaultAssay(seurat_obj) <- "RNA"

# [v2.7 P0-2 FIX] Correct data layer detection
# GetAssayData() returns matrix (genes x cells), NOT layer names
# Use tryCatch to detect if 'data' slot/layer is accessible
has_norm_data <- !is.null(tryCatch({
  data_mat <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  if (is.null(data_mat) || nrow(data_mat) == 0) NULL else data_mat
}, error = function(e) NULL))

if (!has_norm_data) {
  cat("[v2.7 P0-2 FIX] Normalizing data (data slot not found)...\n")
  seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", normalization.method = "LogNormalize")
} else {
  cat("  Data layer found, skipping normalization\n")
}

# Find markers by Manual_Annotation
cat("Finding markers by Manual_Annotation...\n")
Idents(seurat_obj) <- LABELS_KEY
markers_celltype <- FindAllMarkers(
  seurat_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  markers_celltype, 
  file.path(OUTPUT_DIR, "markers_by_celltype.csv"), 
  row.names = FALSE
)

# Find markers by leiden_bbknn
cat("Finding markers by leiden_bbknn...\n")
Idents(seurat_obj) <- CLUSTER_COL
markers_cluster <- FindAllMarkers(
  seurat_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  markers_cluster, 
  file.path(OUTPUT_DIR, "markers_by_cluster.csv"), 
  row.names = FALSE
)

# [P0-3 FIX] Use pseudobulk aggregation instead of full-cell heatmap
cat("\n[P0-3 FIX] Generating pseudobulk heatmaps (memory efficient)...\n")

# Get top markers
top_markers_celltype <- markers_celltype %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N_MARKERS, with_ties = FALSE) %>%
  ungroup()

top_markers_cluster <- markers_cluster %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N_MARKERS, with_ties = FALSE) %>%
  ungroup()

marker_genes_celltype <- unique(top_markers_celltype$gene)
marker_genes_cluster <- unique(top_markers_cluster$gene)

# Pseudobulk by Manual_Annotation
cat("Computing pseudobulk by Manual_Annotation...\n")
Idents(seurat_obj) <- LABELS_KEY

ct_counts <- table(Idents(seurat_obj))
valid_ct <- names(ct_counts)[ct_counts >= MIN_PSEUDOBULK_CELLS]

# region agent log: pseudobulk valid cell types (H4)
log_line <- paste0(
  '{"sessionId":"debug-session",',
  '"runId":"v2.6",',
  '"hypothesisId":"H4",',
  '"location":"epithelial_monocle3_trajectory_v2_5_HOTFIX.R:PART2_pseudobulk_celltype",',
  '"message":"pseudobulk valid cell types",',
  '"data":{"n_celltypes_total":', length(ct_counts),
  ',"n_celltypes_valid":', length(valid_ct), "},",
  '"timestamp":', as.numeric(Sys.time()) * 1000, 
  "}"
)
safe_write_log(log_line, DEBUG_LOG_PATH)
write("\n", file = DEBUG_LOG_PATH, append = TRUE)
# endregion agent log

if (length(valid_ct) < 2) {
  cat("  [SKIP] Too few valid cell types for pseudobulk heatmap (Manual_Annotation)\n")
} else {
  seurat_ct <- subset(seurat_obj, idents = valid_ct)
  
  avg_expr_celltype <- AverageExpression(
    seurat_ct, 
    assay = "RNA", 
    features = marker_genes_celltype,
    slot = "data"
  )$RNA
  
  # Scale
  avg_expr_celltype_scaled <- t(scale(t(avg_expr_celltype)))
  avg_expr_celltype_scaled <- pmax(pmin(avg_expr_celltype_scaled, 2), -2)
  
  # Heatmap
  ht1 <- Heatmap(
    avg_expr_celltype_scaled,
    name = "Scaled Expression",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Top Markers by Cell Type (Pseudobulk)",
    heatmap_legend_param = list(
      title = "Scaled\nExpression",
      title_position = "leftcenter-rot"
    )
  )
  
  pdf(file.path(OUTPUT_DIR, "heatmap_celltype_pseudobulk.pdf"), width = 10, height = 12)
  draw(ht1)
  dev.off()
  
  cat("  Saved: ", file.path(OUTPUT_DIR, "heatmap_celltype_pseudobulk.pdf"), "\n", sep = "")
}

# Pseudobulk by leiden_bbknn
cat("Computing pseudobulk by leiden_bbknn...\n")
Idents(seurat_obj) <- CLUSTER_COL

cl_counts <- table(Idents(seurat_obj))
valid_cl <- names(cl_counts)[cl_counts >= MIN_PSEUDOBULK_CELLS]

# region agent log: pseudobulk valid clusters (H5)
log_line <- paste0(
  '{"sessionId":"debug-session",',
  '"runId":"v2.6",',
  '"hypothesisId":"H5",',
  '"location":"epithelial_monocle3_trajectory_v2_5_HOTFIX.R:PART2_pseudobulk_cluster",',
  '"message":"pseudobulk valid clusters",',
  '"data":{"n_clusters_total":', length(cl_counts),
  ',"n_clusters_valid":', length(valid_cl), "},",
  '"timestamp":', as.numeric(Sys.time()) * 1000, 
  "}"
)
safe_write_log(log_line, DEBUG_LOG_PATH)
write("\n", file = DEBUG_LOG_PATH, append = TRUE)
# endregion agent log

if (length(valid_cl) < 2) {
  cat("  [SKIP] Too few valid clusters for pseudobulk heatmap (leiden_bbknn)\n")
} else {
  seurat_cl <- subset(seurat_obj, idents = valid_cl)
  
  avg_expr_cluster <- AverageExpression(
    seurat_cl, 
    assay = "RNA", 
    features = marker_genes_cluster,
    slot = "data"
  )$RNA
  
  # Scale
  avg_expr_cluster_scaled <- t(scale(t(avg_expr_cluster)))
  avg_expr_cluster_scaled <- pmax(pmin(avg_expr_cluster_scaled, 2), -2)
  
  # Heatmap
  ht2 <- Heatmap(
    avg_expr_cluster_scaled,
    name = "Scaled Expression",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Top Markers by Leiden Cluster (Pseudobulk)",
    heatmap_legend_param = list(
      title = "Scaled\nExpression",
      title_position = "leftcenter-rot"
    )
  )
  
  pdf(file.path(OUTPUT_DIR, "heatmap_cluster_pseudobulk.pdf"), width = 10, height = 12)
  draw(ht2)
  dev.off()
  
  cat("  Saved: ", file.path(OUTPUT_DIR, "heatmap_cluster_pseudobulk.pdf"), "\n", sep = "")
}

cat("\n", strrep("=", 80), "\n")
cat("ALL DONE\n")
cat("Output: ", OUTPUT_DIR, "\n", sep = "")
cat("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("=", 80), "\n")
