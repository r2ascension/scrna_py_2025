#!/usr/bin/env Rscript
################################################################################
# Epithelial Trajectory by Tissue with Monocle3
# - Load Seurat from RDS
# - Split into Top4 tissues (by cell count)
# - Root = Basal (auto match, case-insensitive)
# - Save per-tissue outputs (pdf + csv + cds)
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
cat("EPITHELIAL TRAJECTORY BY TISSUE (ROOT = BASAL)\n")
cat(strrep("=", 80), "\n\n", sep = "")
cat("monocle3:", as.character(packageVersion("monocle3")), "\n")
cat("Seurat  :", as.character(packageVersion("Seurat")), "\n")
cat("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# CONFIG
# ==============================================================================
INPUT_RDS <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/py/1218/Epi_monocle3_by_tissue_rootBasal_20251218"

# 你需要确认这两个 key 在 meta.data 里真实存在
LABELS_KEY <- "Manual_Annotation" # 必须含 Basal（或类似名称）
TISSUE_KEY <- "tissue" # tissue 分组列
ROOT_LABEL <- "Basal"

# Marker 分析相关配置（用于后半部分热图）
ID_COL <- "Manual_Annotation" # 细胞种类（与 LABELS_KEY 相同）
CLUSTER_COL <- "leiden_bbknn" # 聚类列名
TOP_N <- 10 # 每个 cluster 取前 N 个 marker genes

N_DIM_PREPROCESS <- 50
N_NEIGHBORS_CLUSTER <- 15
USE_PARTITION <- TRUE
RUN_GRAPH_TEST <- TRUE # 是否运行 graph_test（可能较慢）

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
  # Seurat v5 layer vs v4 slot 兼容
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    x <- tryCatch(
      Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts"),
      error = function(e) NULL
    )
    if (!is.null(x)) return(x)
    # 有些对象 counts layer 不叫 counts：尽量回退到 slot=counts
  }
  x <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts"),
    error = function(e) NULL
  )
  if (is.null(x)) {
    stop(
      "Cannot find raw counts in assay=",
      assay,
      " (layer/slot 'counts' missing)."
    )
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

  # Transfer scVI/scANVI latent if present (optional)
  latent_red <- red_names[grepl("scanvi|scvi", red_names, ignore.case = TRUE)][
    1
  ]
  if (!is.na(latent_red) && !is.null(latent_red)) {
    latent <- Seurat::Embeddings(seurat_obj, reduction = latent_red)
    reducedDims(cds)[["scANVI"]] <- latent
  }

  cds
}

pick_root_label <- function(labels_vec, root_label = "Basal") {
  labs <- unique(as.character(labels_vec))
  if (root_label %in% labs) {
    return(root_label)
  }

  hit <- labs[grepl(root_label, labs, ignore.case = TRUE)]
  if (length(hit) > 0) {
    message(
      "  [WARNING] Root label '",
      root_label,
      "' not found exactly, using '",
      hit[1],
      "' instead"
    )
    return(hit[1])
  }

  tb <- sort(table(labels_vec), decreasing = TRUE)
  fallback_label <- as.character(names(tb)[1])
  message(
    "  [WARNING] Root label '",
    root_label,
    "' not found, using most frequent label '",
    fallback_label,
    "' as fallback"
  )
  return(fallback_label)
}

get_root_principal_node <- function(cds, cell_type_key, root_cell_type) {
  root_cells <- which(colData(cds)[[cell_type_key]] == root_cell_type)
  if (length(root_cells) == 0) {
    stop("No cells found for root cell type: ", root_cell_type)
  }

  # 检查 principal graph 是否存在
  if (!"UMAP" %in% names(monocle3::principal_graph(cds))) {
    stop("Principal graph 'UMAP' not found. Did learn_graph() succeed?")
  }
  if (!"UMAP" %in% names(cds@principal_graph_aux)) {
    stop("Principal graph aux 'UMAP' not found. Did learn_graph() succeed?")
  }
  if (
    is.null(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
  ) {
    stop("pr_graph_cell_proj_closest_vertex not found in principal_graph_aux")
  }

  closest_vertex <- cds@principal_graph_aux[[
    "UMAP"
  ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])[,
    1
  ]

  pr_graph <- monocle3::principal_graph(cds)[["UMAP"]]
  root_vertex_numeric <- as.numeric(names(which.max(table(closest_vertex[
    root_cells
  ]))))
  igraph::V(pr_graph)$name[root_vertex_numeric]
}

assign_branch_by_leaf <- function(cds) {
  g <- monocle3::principal_graph(cds)[["UMAP"]]
  deg <- igraph::degree(g)
  leaves <- names(deg[deg == 1])

  if (length(leaves) < 2) {
    colData(cds)$branch <- "Main"
    return(cds)
  }

  closest_vertex <- cds@principal_graph_aux[[
    "UMAP"
  ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])[,
    1
  ]
  closest_vertex <- igraph::V(g)$name[as.numeric(closest_vertex)]

  # For each leaf, compute shortest distance from each cell's closest vertex
  dist_list <- lapply(leaves, function(lf) {
    d <- igraph::distances(g, v = closest_vertex, to = lf)
    as.numeric(d[, 1])
  })
  dist_mat <- do.call(cbind, dist_list)
  colnames(dist_mat) <- leaves

  leaf_idx <- apply(dist_mat, 1, which.min)
  colData(cds)$branch <- paste0("Leaf_", leaf_idx)
  cds
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
    cat(
      "  [SKIP] tissue=",
      tissue_value,
      " too few cells: ",
      length(cells),
      "\n",
      sep = ""
    )
    return(NULL)
  }

  seu <- subset(seurat_obj, cells = cells)

  root_label <- pick_root_label(seu@meta.data[[LABELS_KEY]], ROOT_LABEL)
  cat(
    "  tissue=",
    tissue_value,
    " | cells=",
    ncol(seu),
    " | root_label=",
    root_label,
    "\n",
    sep = ""
  )

  cds <- seurat_to_cds(seu, assay = DefaultAssay(seu))

  # preprocess
  if ("scANVI" %in% names(reducedDims(cds))) {
    reducedDims(cds)[["PCA"]] <- reducedDims(cds)[["scANVI"]]
  } else {
    cds <- preprocess_cds(cds, num_dim = N_DIM_PREPROCESS)
  }

  # UMAP
  if (!("UMAP" %in% names(reducedDims(cds)))) {
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
  }

  # cluster + graph
  cds <- cluster_cells(cds, k = N_NEIGHBORS_CLUSTER)
  cds <- learn_graph(cds, use_partition = USE_PARTITION)

  # order cells
  root_node <- get_root_principal_node(cds, LABELS_KEY, root_label)
  cds <- order_cells(cds, root_pr_nodes = root_node)

  # branch label
  cds <- assign_branch_by_leaf(cds)

  # plots
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
    ggtitle(paste0("tissue=", tissue_value, " | branch (leaf assignment)"))
  save_plot_pdf(p3, file.path(fig_dir, "trajectory_branch.pdf"))

  # export pseudotime table
  pt <- data.frame(
    cell_id = colnames(cds),
    tissue = tissue_value,
    cell_type = colData(cds)[[LABELS_KEY]],
    pseudotime = pseudotime(cds),
    partition = partitions(cds),
    cluster = clusters(cds),
    branch = colData(cds)$branch
  )
  write.csv(pt, file.path(out_dir, "pseudotime.csv"), row.names = FALSE)

  # graph_test genes (optional)
  if (RUN_GRAPH_TEST) {
    cat("  graph_test ...\n")
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
      write.csv(
        gt,
        file.path(out_dir, "pseudotime_de_genes.csv"),
        row.names = FALSE
      )
    }
  } else {
    cat("  [SKIP] graph_test (RUN_GRAPH_TEST=FALSE)\n")
  }

  # save cds
  cds_dir <- file.path(out_dir, "cds")
  dir.create(cds_dir, recursive = TRUE, showWarnings = FALSE)
  if ("save_monocle_objects" %in% getNamespaceExports("monocle3")) {
    monocle3::save_monocle_objects(cds, directory = cds_dir)
  } else {
    saveRDS(cds, file.path(cds_dir, "cds.rds"))
  }

  cds
}

# ==============================================================================
# RUN
# ==============================================================================
if (!file.exists(INPUT_RDS)) {
  stop("RDS not found: ", INPUT_RDS)
}
seurat_obj <- readRDS(INPUT_RDS)
if (!inherits(seurat_obj, "Seurat")) {
  stop("Loaded object is not a Seurat object.")
}

cat(
  "Loaded Seurat: cells=",
  ncol(seurat_obj),
  " genes=",
  nrow(seurat_obj),
  "\n",
  sep = ""
)

if (!TISSUE_KEY %in% colnames(seurat_obj@meta.data)) {
  stop("TISSUE_KEY not in meta.data: ", TISSUE_KEY)
}
if (!LABELS_KEY %in% colnames(seurat_obj@meta.data)) {
  stop("LABELS_KEY not in meta.data: ", LABELS_KEY)
}

# choose top4 tissues by cell count
tissue_counts <- sort(
  table(seurat_obj@meta.data[[TISSUE_KEY]]),
  decreasing = TRUE
)
tissues_top4 <- names(tissue_counts)[seq_len(min(4, length(tissue_counts)))]

cat("Selected tissues (Top4 by cell count):\n")
for (tt in tissues_top4) {
  cat("  ", tt, ": ", tissue_counts[[tt]], "\n", sep = "")
}
cat("\n")

cds_list <- list()
for (tt in tissues_top4) {
  cat(strrep("-", 80), "\n")
  cat("Running tissue: ", tt, "\n", sep = "")
  out_dir <- file.path(OUTPUT_DIR, paste0("tissue_", sanitize_name(tt)))
  cds_list[[tt]] <- run_monocle3_for_tissue(seurat_obj, tt, out_dir)
}

cat("\n", strrep("=", 80), "\n", sep = "")
cat("DONE\n")
cat("Output: ", OUTPUT_DIR, "\n", sep = "")
cat("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("=", 80), "\n")

# -----------------------------
# Marker analysis and heatmap
# -----------------------------
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- ID_COL

# -----------------------------
# Find markers
# -----------------------------
markers <- FindAllMarkers(
  seurat_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, file.path(OUTPUT_DIR, "all_markers.csv"), row.names = FALSE)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

marker_genes <- top_markers$gene
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]
marker_genes <- marker_genes[!duplicated(marker_genes)]

# 检查 marker_genes 是否为空或过少
if (length(marker_genes) == 0) {
  stop("No marker genes found after filtering. Check FindAllMarkers results and gene names.")
}
if (length(marker_genes) < 5) {
  warning("Only ", length(marker_genes), " marker genes found. Heatmap may not be informative.")
}

cat("Using ", length(marker_genes), " marker genes for heatmap\n", sep = "")

# 只对 marker genes 做 scale（省时省内存）
seurat_obj <- ScaleData(
  seurat_obj,
  assay = "RNA",
  features = marker_genes,
  verbose = FALSE
)

# -----------------------------
# Get matrix (Seurat v5 优先 layer；失败则 fallback slot)
# -----------------------------
.get_scaled <- function(obj, assay = "RNA") {
  out <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "scale.data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "scale.data")
  )
  out
}

mat <- .get_scaled(seurat_obj, assay = "RNA")[marker_genes, , drop = FALSE]
mat <- as.matrix(mat)

# 可选：截断极端值，让热图更可读
mat <- pmax(pmin(mat, 2), -2)

cluster_vec <- seurat_obj[[CLUSTER_COL]][, 1]
celltype_vec <- seurat_obj[[ID_COL]][, 1]

cluster_fac <- factor(cluster_vec)
celltype_fac <- factor(celltype_vec)

# 颜色（自动）
cluster_cols <- setNames(hue_pal()(nlevels(cluster_fac)), levels(cluster_fac))
celltype_cols <- setNames(
  hue_pal()(nlevels(celltype_fac)),
  levels(celltype_fac)
)

ha_top <- HeatmapAnnotation(
  Cluster = cluster_fac,
  col = list(Cluster = cluster_cols),
  show_annotation_name = TRUE
)

ha_bottom <- HeatmapAnnotation(
  CellType = celltype_fac,
  col = list(CellType = celltype_cols),
  show_annotation_name = TRUE
)

# -----------------------------
# Heatmap: 全局列聚类（不分块）
# -----------------------------
ht <- Heatmap(
  mat,
  name = "Scaled",
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  cluster_rows = TRUE,
  cluster_columns = TRUE, # 全局对所有细胞做聚类
  clustering_method_columns = "ward.D2",
  clustering_distance_columns = "pearson",
  top_annotation = ha_top,
  bottom_annotation = ha_bottom,
  use_raster = TRUE
)

# 保存热图到输出目录
heatmap_file <- file.path(OUTPUT_DIR, "epi_subtype_marker.pdf")
pdf(heatmap_file, width = 14, height = 10, onefile = TRUE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Heatmap saved: ", heatmap_file, "\n", sep = "")
