#!/usr/bin/env Rscript
################################################################################
# Epithelial Trajectory Analysis - PRODUCTION VERSION v2.0
# 
# KEY IMPROVEMENTS:
# - [P0-1] FIXED: Always run preprocess_cds() before graph_test
# - [P0-2] FIXED: Robust vertex type handling (numeric vs character)
# - [P0-3] FIXED: Use pseudobulk aggregation for heatmap (no full-cell clustering)
# - [P1-1] FIXED: Explicit reduction_method specification
# - [P1-3] FIXED: Strict root validation (skip tissue if no Basal found)
# - [P1-4] ADDED: Leaf composition analysis
# - ADDED: Biological sanity checks (pseudotime-marker trends)
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
  library(scales)
  library(grid)
})

cat("\n", strrep("=", 80), "\n", sep = "")
cat("EPITHELIAL TRAJECTORY ANALYSIS - PRODUCTION v2.0\n")
cat(strrep("=", 80), "\n\n", sep = "")
cat("monocle3:", as.character(packageVersion("monocle3")), "\n")
cat("Seurat  :", as.character(packageVersion("Seurat")), "\n")
cat("igraph  :", as.character(packageVersion("igraph")), "\n")
cat("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# region agent log: script entry (H0)
DEBUG_LOG_PATH <- "/home/h2048/.cursor/debug.log"
log_line <- paste0(
  '{"sessionId":"debug-session",',
  '"runId":"pre-fix",',
  '"hypothesisId":"H0",',
  '"location":"epithelial_monocle3_trajectory_20251218_v2_2.R:global",',
  '"message":"script entry and package versions",',
  '"data":{"monocle3":"', as.character(packageVersion("monocle3")),
  '","Seurat":"', as.character(packageVersion("Seurat")),
  '","igraph":"', as.character(packageVersion("igraph")), '"},',
  '"timestamp":', as.numeric(Sys.time()) * 1000,
  "}"
)
write(log_line, file = DEBUG_LOG_PATH, append = TRUE)
write("\n", file = DEBUG_LOG_PATH, append = TRUE)
# endregion agent log

# ==============================================================================
# CONFIG
# ==============================================================================
INPUT_RDS <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/py/1218/Epi_monocle3_PRODUCTION_v2_20251218"

LABELS_KEY <- "Manual_Annotation"
TISSUE_KEY <- "tissue"
ROOT_LABEL <- "Basal"
CLUSTER_COL <- "leiden_bbknn"

# [P1-3 FIX] Strict root validation - skip tissue if no Basal found
ALLOW_ROOT_FALLBACK <- FALSE  # Set TRUE to use old fallback behavior

# Monocle3 parameters
N_DIM_PREPROCESS <- 50
N_NEIGHBORS_CLUSTER <- 15
USE_PARTITION <- TRUE
RUN_GRAPH_TEST <- TRUE
REDUCTION_METHOD <- "UMAP"  # [P1-1] Explicit specification; required for learn_graph in monocle3>=1.4

# Marker analysis
TOP_N_MARKERS <- 10

# Lineage markers for sanity check
LINEAGE_MARKERS <- list(
  Basal = c("KRT5", "TP63", "KRT14"),
  Ciliated = c("FOXJ1", "PIFO", "TPPP3"),
  Goblet = c("MUC5AC", "TFF3", "SPDEF"),
  Secretory = c("SCGB1A1", "SCGB3A1", "BPIFA1"),
  Squamous = c("KRT13", "IVL", "SPRR2A")
)

# Debug-mode constant (for log ingestion)
DEBUG_LOG_PATH <- "/home/h2048/.cursor/debug.log"

FIG_W <- 10
FIG_H <- 8

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Helpers
# ==============================================================================
sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

save_plot_pdf <- function(p, file, width = FIG_W, height = FIG_H) {
  pdf(file, width = width, height = height, onefile = TRUE)
  print(p)
  dev.off()
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

seurat_to_cds <- function(seurat_obj, assay = "RNA") {
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
  
  # Transfer UMAP if present
  red_names <- names(seurat_obj@reductions)
  umap_red <- red_names[grepl("umap", red_names, ignore.case = TRUE)][1]
  if (!is.na(umap_red) && !is.null(umap_red)) {
    umap_coords <- Seurat::Embeddings(seurat_obj, reduction = umap_red)
    colnames(umap_coords) <- paste0("UMAP_", seq_len(ncol(umap_coords)))
    reducedDims(cds)[["UMAP"]] <- umap_coords
  }
  
  cds
}

# [P1-3 FIX] Strict root validation
pick_root_label <- function(labels_vec, root_label = "Basal", allow_fallback = FALSE) {
  labs <- unique(as.character(labels_vec))
  
  # Exact match
  if (root_label %in% labs) {
    return(root_label)
  }
  
  # Case-insensitive fuzzy match
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

# [P0-2 FIX] Robust branch assignment
assign_branch_by_leaf <- function(cds) {
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
  
  # [P0-2 FIX] Handle both numeric and character vertex types
  if (is.numeric(closest_vertex)) {
    closest_vertex_names <- vertex_names[closest_vertex]
  } else {
    closest_vertex_names <- as.character(closest_vertex)
  }
  
  # Compute distance from each cell to each leaf
  dist_list <- lapply(leaves, function(lf) {
    d <- igraph::distances(g, v = closest_vertex_names, to = lf)
    as.numeric(d[, 1])
  })
  dist_mat <- do.call(cbind, dist_list)
  colnames(dist_mat) <- leaves
  
  leaf_idx <- apply(dist_mat, 1, which.min)
  colData(cds)$branch <- paste0("Leaf_", leaf_idx)
  
  cds
}

# [P1-4 NEW] Leaf composition analysis
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

# [SANITY CHECK] Pseudotime-marker trends
plot_pseudotime_markers <- function(cds, markers, out_file) {
  pt <- pseudotime(cds)
  valid_cells <- !is.na(pt)
  
  if (sum(valid_cells) < 50) {
    cat("  [SKIP] Too few valid pseudotime cells for marker trend plot\n")
    return(NULL)
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
  
  cat("  Plotting ", length(markers_available), " lineage markers\n", sep = "")

  # region agent log: pseudotime & marker availability (H3)
  log_line <- paste0(
    '{"sessionId":"debug-session",',
    '"runId":"pre-fix",',
    '"hypothesisId":"H3",',
    '"location":"epithelial_monocle3_trajectory_20251218_v2_2.R:plot_pseudotime_markers",',
    '"message":"pseudotime and marker availability",',
    '"data":{"n_valid_cells":', sum(valid_cells),
    ',"n_markers_available":', length(markers_available), "},",
    '"timestamp":', as.numeric(Sys.time()) * 1000, 
    "}"
  )
  write(log_line, file = DEBUG_LOG_PATH, append = TRUE)
  write("\n", file = DEBUG_LOG_PATH, append = TRUE)
  # endregion agent log
  
  pdf(out_file, width = 12, height = 8)
  
  for (lineage_name in names(markers)) {
    genes <- markers[[lineage_name]]
    genes <- genes[genes %in% all_genes]
    
    if (length(genes) == 0) next
    
    for (gene in genes) {
      expr <- exprs(cds)[gene, valid_cells]
      
      df <- data.frame(pseudotime = pt_valid, expression = expr)
      
      p <- ggplot(df, aes(x = pseudotime, y = expression)) +
        geom_point(alpha = 0.1, size = 0.5) +
        geom_smooth(method = "loess", color = "red", se = TRUE) +
        labs(
          title = paste0(lineage_name, " marker: ", gene),
          x = "Pseudotime",
          y = "Expression"
        ) +
        theme_minimal()
      
      print(p)
    }
  }
  
  dev.off()
  cat("  Saved: ", out_file, "\n", sep = "")
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
    '"runId":"pre-fix",',
    '"hypothesisId":"H1",',
    '"location":"epithelial_monocle3_trajectory_20251218_v2_2.R:run_monocle3_for_tissue",',
    '"message":"tissue basic stats and root label",',
    '"data":{"tissue":"', tissue_value,
    '","n_cells":', length(cells),
    ',"root_label":"', ifelse(is.null(root_label), "NULL", root_label), '"},',
    '"timestamp":', as.numeric(Sys.time()) * 1000, 
    "}"
  )
  write(log_line, file = DEBUG_LOG_PATH, append = TRUE)
  write("\n", file = DEBUG_LOG_PATH, append = TRUE)
  # endregion agent log
  
  if (is.null(root_label)) {
    cat("  [SKIP] tissue=", tissue_value, " - root not found and fallback disabled\n", sep = "")
    return(NULL)
  }
  
  cat("  tissue=", tissue_value, " | cells=", ncol(seu), " | root=", root_label, "\n", sep = "")
  
  cds <- seurat_to_cds(seu, assay = DefaultAssay(seu))
  
  # [P0-1 FIX] ALWAYS run preprocess_cds
  cat("  Running preprocess_cds...\n")
  cds <- preprocess_cds(cds, num_dim = N_DIM_PREPROCESS, method = "PCA")
  
  # Optional: use scANVI latent as PCA if available
  red_names <- names(reducedDims(cds))
  if (any(grepl("scanvi|scvi", red_names, ignore.case = TRUE))) {
    latent_red <- red_names[grepl("scanvi|scvi", red_names, ignore.case = TRUE)][1]
    cat("  Found ", latent_red, ", using as PCA replacement\n", sep = "")
    reducedDims(cds)[["PCA"]] <- reducedDims(cds)[[latent_red]]
  }
  
  # UMAP
  if (!("UMAP" %in% names(reducedDims(cds)))) {
    cat("  Running reduce_dimension...\n")
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
  }
  
  # [P1-1 FIX] Explicit reduction_method
  cat("  Running cluster_cells (reduction_method=", REDUCTION_METHOD, ")...\n", sep = "")
  cds <- cluster_cells(cds, reduction_method = REDUCTION_METHOD, k = N_NEIGHBORS_CLUSTER)
  
  cat("  Running learn_graph...\n")
  cds <- learn_graph(cds, use_partition = USE_PARTITION)
  
  # Order cells with robust root finding
  cat("  Finding root node...\n")
  root_node <- get_root_principal_node(cds, LABELS_KEY, root_label)
  cat("  Root node: ", root_node, "\n", sep = "")

  # region agent log: root node & graph summary (H2)
  g <- monocle3::principal_graph(cds)[["UMAP"]]
  log_line <- paste0(
    '{"sessionId":"debug-session",',
    '"runId":"pre-fix",',
    '"hypothesisId":"H2",',
    '"location":"epithelial_monocle3_trajectory_20251218_v2_2.R:run_monocle3_for_tissue",',
    '"message":"root node and graph summary",',
    '"data":{"tissue":"', tissue_value,
    '","root_node":"', root_node,
    '","n_vertices":', igraph::gorder(g),
    ',"n_edges":', igraph::gsize(g), "},",
    '"timestamp":', as.numeric(Sys.time()) * 1000, 
    "}"
  )
  write(log_line, file = DEBUG_LOG_PATH, append = TRUE)
  write("\n", file = DEBUG_LOG_PATH, append = TRUE)
  # endregion agent log
  
  cds <- order_cells(cds, root_pr_nodes = root_node)
  
  # Branch assignment
  cds <- assign_branch_by_leaf(cds)
  
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
  
  # [SANITY CHECK] Pseudotime-marker trends
  cat("  Plotting pseudotime-marker trends...\n")
  plot_pseudotime_markers(
    cds, 
    LINEAGE_MARKERS, 
    file.path(fig_dir, "pseudotime_marker_trends.pdf")
  )
  
  # graph_test genes
  if (RUN_GRAPH_TEST) {
    cat("  Running graph_test...\n")
    gt <- tryCatch(
      monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4),
      error = function(e) {
        message("  [WARNING] graph_test failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    if (!is.null(gt)) {
      gt <- gt[order(gt$q_value), , drop = FALSE]
      gt$gene_id <- rownames(gt)
      rownames(gt) <- NULL
      write.csv(gt, file.path(out_dir, "pseudotime_de_genes.csv"), row.names = FALSE)
    }
  }
  
  # Save cds
  cds_dir <- file.path(out_dir, "cds")
  dir.create(cds_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cds, file.path(cds_dir, "cds.rds"))
  
  cds
}

# ==============================================================================
# RUN PART 1: Trajectory Analysis
# ==============================================================================
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

# ==============================================================================
# PART 2: Marker Analysis with Pseudobulk Heatmap
# ==============================================================================
cat(strrep("=", 80), "\n")
cat("PART 2: MARKER ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

DefaultAssay(seurat_obj) <- "RNA"

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
  '"runId":"pre-fix",',
  '"hypothesisId":"H4",',
  '"location":"epithelial_monocle3_trajectory_20251218_v2_2.R:PART2_pseudobulk_celltype",',
  '"message":"pseudobulk valid cell types",',
  '"data":{"n_celltypes_total":', length(ct_counts),
  ',"n_celltypes_valid":', length(valid_ct), "},",
  '"timestamp":', as.numeric(Sys.time()) * 1000, 
  "}"
)
write(log_line, file = DEBUG_LOG_PATH, append = TRUE)
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
    name = "Scaled",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Markers by Cell Type (Pseudobulk)"
  )
  
  pdf(file.path(OUTPUT_DIR, "heatmap_celltype_pseudobulk.pdf"), width = 10, height = 12)
  draw(ht1)
  dev.off()
}

# Pseudobulk by leiden_bbknn
cat("Computing pseudobulk by leiden_bbknn...\n")
Idents(seurat_obj) <- CLUSTER_COL

cl_counts <- table(Idents(seurat_obj))
valid_cl <- names(cl_counts)[cl_counts >= MIN_PSEUDOBULK_CELLS]

# region agent log: pseudobulk valid clusters (H5)
log_line <- paste0(
  '{"sessionId":"debug-session",',
  '"runId":"pre-fix",',
  '"hypothesisId":"H5",',
  '"location":"epithelial_monocle3_trajectory_20251218_v2_2.R:PART2_pseudobulk_cluster",',
  '"message":"pseudobulk valid clusters",',
  '"data":{"n_clusters_total":', length(cl_counts),
  ',"n_clusters_valid":', length(valid_cl), "},",
  '"timestamp":', as.numeric(Sys.time()) * 1000, 
  "}"
)
write(log_line, file = DEBUG_LOG_PATH, append = TRUE)
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
    name = "Scaled",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Markers by Leiden Cluster (Pseudobulk)"
  )
  
  pdf(file.path(OUTPUT_DIR, "heatmap_cluster_pseudobulk.pdf"), width = 10, height = 12)
  draw(ht2)
  dev.off()
}

cat("\n", strrep("=", 80), "\n")
cat("ALL DONE\n")
cat("Output: ", OUTPUT_DIR, "\n", sep = "")
cat("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("=", 80), "\n")
