#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

scalar_chr <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("%s must be a non-empty character scalar", arg), call. = FALSE)
  }
  x
}

dir_create_safe <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_tsv_safe <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

split_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  vals[nzchar(vals)]
}

wrap_value <- function(x, width = NULL) {
  if (is.null(width) || is.na(width) || width <= 0L || !nzchar(x)) return(x)
  paste(strwrap(x, width = width), collapse = "\n")
}

display_width_chars <- function(values) {
  if (!length(values)) return(0)
  max(vapply(
    strsplit(as.character(values), "\n", fixed = TRUE),
    function(parts) {
      if (!length(parts)) return(0L)
      max(nchar(parts, type = "width"), na.rm = TRUE)
    },
    integer(1L)
  ), na.rm = TRUE)
}

get_assay_matrix_safe <- function(seurat_obj, assay = NULL, slot = "data", layer = NULL) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object", call. = FALSE)
  assay <- assay %||% SeuratObject::DefaultAssay(seurat_obj)
  get_args <- names(formals(SeuratObject::GetAssayData))

  if (!is.null(layer) && "layer" %in% get_args) {
    return(SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer))
  }
  if (!is.null(slot) && "slot" %in% get_args) {
    return(SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot))
  }
  if (!is.null(layer) && "slot" %in% get_args) {
    return(SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = layer))
  }
  SeuratObject::GetAssayData(seurat_obj, assay = assay)
}

extract_group_labels <- function(seurat_obj, group_col, sort_labels = TRUE, drop_tip_labels = NULL) {
  meta <- seurat_obj@meta.data
  if (!group_col %in% colnames(meta)) {
    stop(sprintf("group_col '%s' not found in meta.data", group_col), call. = FALSE)
  }
  labels <- as.character(meta[[group_col]])
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) stop(sprintf("group_col '%s' has no non-empty labels", group_col), call. = FALSE)
  groups <- unique(labels)
  if (sort_labels) groups <- sort(groups)
  if (!is.null(drop_tip_labels) && length(drop_tip_labels)) {
    groups <- setdiff(groups, as.character(drop_tip_labels))
  }
  if (!length(groups)) stop("No groups remain after filtering drop_tip_labels", call. = FALSE)
  groups
}

write_hierarchy_annotation_template <- function(seurat_obj,
                                                group_col = "short.label",
                                                out_tsv,
                                                drop_tip_labels = NULL,
                                                sort_labels = TRUE,
                                                extra_cols = c("display_label", "level_1", "level_2", "level_3", "meaning")) {
  groups <- extract_group_labels(
    seurat_obj = seurat_obj,
    group_col = group_col,
    sort_labels = sort_labels,
    drop_tip_labels = drop_tip_labels
  )
  template <- data.frame(label = groups, stringsAsFactors = FALSE)
  for (col_name in extra_cols) template[[col_name]] <- ""
  template$display_label <- groups
  write_tsv_safe(template, out_tsv)
}

compute_cluster_means <- function(seurat_obj,
                                  group_col = "short.label",
                                  assay = NULL,
                                  slot = "data",
                                  layer = NULL,
                                  drop_tip_labels = NULL,
                                  sort_labels = TRUE) {
  meta <- seurat_obj@meta.data
  groups <- extract_group_labels(
    seurat_obj = seurat_obj,
    group_col = group_col,
    sort_labels = sort_labels,
    drop_tip_labels = drop_tip_labels
  )
  labels <- as.character(meta[[group_col]])
  expr_mat <- get_assay_matrix_safe(seurat_obj, assay = assay, slot = slot, layer = layer)
  row_mean_fn <- if (inherits(expr_mat, "Matrix")) Matrix::rowMeans else base::rowMeans

  mean_list <- lapply(groups, function(group_id) {
    idx <- which(labels == group_id)
    if (!length(idx)) stop(sprintf("No cells found for label '%s'", group_id), call. = FALSE)
    row_mean_fn(expr_mat[, idx, drop = FALSE])
  })

  clus_means <- do.call(rbind, mean_list)
  rownames(clus_means) <- groups
  as.matrix(clus_means)
}

build_cluster_phylo <- function(cluster_means,
                                tree_method = c("hclust", "nj"),
                                dist_method = "euclidean",
                                outgroup = NULL) {
  tree_method <- match.arg(tree_method)
  dist_obj <- stats::dist(cluster_means, method = dist_method)
  phy <- if (identical(tree_method, "nj")) {
    ape::nj(dist_obj)
  } else {
    ape::as.phylo(stats::hclust(dist_obj))
  }

  if (!is.null(outgroup)) {
    outgroup <- as.character(outgroup)
    if (!outgroup %in% phy$tip.label) {
      stop(sprintf("outgroup '%s' not present in tree tip labels", outgroup), call. = FALSE)
    }
    phy <- ape::root(phy, outgroup = outgroup)
  }

  list(phy = phy, dist = dist_obj)
}

prepare_annotation_table <- function(annotation_df,
                                     tip_order,
                                     tip_col = "label",
                                     tip_label_col = "display_label",
                                     annotation_cols = c("level_1", "level_2", "level_3", "meaning"),
                                     wrap_cols = "meaning",
                                     wrap_width = 28L) {
  if (!tip_col %in% colnames(annotation_df)) {
    stop(sprintf("annotation_df must contain tip_col '%s'", tip_col), call. = FALSE)
  }

  annotation_df[[tip_col]] <- as.character(annotation_df[[tip_col]])
  if (anyDuplicated(annotation_df[[tip_col]])) {
    dup <- unique(annotation_df[[tip_col]][duplicated(annotation_df[[tip_col]])])
    stop(sprintf("annotation_df has duplicated labels: %s", paste(dup, collapse = ", ")), call. = FALSE)
  }

  missing_labels <- setdiff(tip_order, annotation_df[[tip_col]])
  if (length(missing_labels)) {
    stop(sprintf(
      "annotation_df is missing labels for %d tree tips: %s",
      length(missing_labels),
      paste(utils::head(missing_labels, 10L), collapse = ", ")
    ), call. = FALSE)
  }

  annotation_df <- annotation_df[match(tip_order, annotation_df[[tip_col]]), , drop = FALSE]
  rownames(annotation_df) <- annotation_df[[tip_col]]

  if (!tip_label_col %in% colnames(annotation_df)) {
    annotation_df[[tip_label_col]] <- annotation_df[[tip_col]]
  }

  keep_cols <- unique(c(tip_col, tip_label_col, annotation_cols))
  keep_cols <- keep_cols[keep_cols %in% colnames(annotation_df)]
  annotation_df <- annotation_df[, keep_cols, drop = FALSE]
  for (col_name in colnames(annotation_df)) annotation_df[[col_name]] <- as.character(annotation_df[[col_name]])
  annotation_df[[tip_label_col]][!nzchar(annotation_df[[tip_label_col]])] <- annotation_df[[tip_col]][!nzchar(annotation_df[[tip_label_col]])]

  wrap_cols <- intersect(wrap_cols, colnames(annotation_df))
  if (length(wrap_cols)) {
    for (col_name in wrap_cols) {
      annotation_df[[col_name]] <- vapply(
        annotation_df[[col_name]],
        wrap_value,
        character(1L),
        width = as.integer(wrap_width)
      )
    }
  }

  annotation_df
}

build_annotation_text_data <- function(tree_plot,
                                       annotation_df,
                                       tip_label_col = "display_label",
                                       annotation_cols = c("level_1", "level_2", "level_3", "meaning"),
                                       tip_offset = 0.35,
                                       column_gap = 0.8,
                                       char_width = 0.16,
                                       header_offset = 1) {
  tip_df <- subset(tree_plot$data, isTip)
  tip_df <- tip_df[match(rownames(annotation_df), tip_df$label), c("label", "y"), drop = FALSE]

  text_columns <- annotation_cols[annotation_cols %in% colnames(annotation_df)]
  tip_width <- display_width_chars(annotation_df[[tip_label_col]]) * char_width
  current_x <- max(tree_plot$data$x, na.rm = TRUE) + tip_offset + tip_width + column_gap
  column_pos <- numeric(length(text_columns))
  names(column_pos) <- text_columns
  if (length(text_columns)) {
    for (idx in seq_along(text_columns)) {
      col_name <- text_columns[[idx]]
      column_pos[[col_name]] <- current_x
      current_x <- current_x + display_width_chars(annotation_df[[col_name]]) * char_width + column_gap
    }
  }

  ann_long <- do.call(rbind, lapply(text_columns, function(col_name) {
    data.frame(
      x = rep(column_pos[[col_name]], nrow(tip_df)),
      y = tip_df$y,
      column = col_name,
      value = annotation_df[[col_name]],
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(ann_long)) {
    ann_long <- data.frame(x = numeric(0), y = numeric(0), column = character(0), value = character(0), stringsAsFactors = FALSE)
  }
  ann_long$value[is.na(ann_long$value)] <- ""

  headers <- data.frame(
    x = unname(column_pos),
    y = max(tip_df$y, na.rm = TRUE) + header_offset,
    label = names(column_pos),
    stringsAsFactors = FALSE
  )

  list(annotation_text = ann_long, headers = headers, x_max = current_x)
}

plot_cluster_hierarchy_tree <- function(phy,
                                        annotation_df,
                                        tip_label_col = "display_label",
                                        annotation_cols = c("level_1", "level_2", "level_3", "meaning"),
                                        branch.length = "none",
                                        layout = "rectangular",
                                        tip_offset = 0.35,
                                        column_gap = 0.8,
                                        char_width = 0.16,
                                        header_offset = 1,
                                        tip_label_size = 3,
                                        annotation_text_size = 3,
                                        header_text_size = 3.2,
                                        point_size = 2.2,
                                        right_margin_lines = 16) {
  plot_df <- annotation_df
  plot_df$label <- rownames(annotation_df)

  tree_plot <- ggtree::ggtree(phy, layout = layout, branch.length = branch.length) %<+% plot_df
  tree_plot <- tree_plot +
    ggtree::geom_tiplab(ggplot2::aes_string(label = tip_label_col), align = TRUE, hjust = 0, offset = tip_offset, size = tip_label_size) +
    ggtree::geom_tippoint(size = point_size)

  ann_pack <- build_annotation_text_data(
    tree_plot = tree_plot,
    annotation_df = annotation_df,
    tip_label_col = tip_label_col,
    annotation_cols = annotation_cols,
    tip_offset = tip_offset,
    column_gap = column_gap,
    char_width = char_width,
    header_offset = header_offset
  )

  tree_plot <- tree_plot +
    ggplot2::geom_text(
      data = ann_pack$annotation_text,
      mapping = ggplot2::aes(x = x, y = y, label = value),
      inherit.aes = FALSE,
      hjust = 0,
      lineheight = 0.95,
      size = annotation_text_size
    ) +
    ggplot2::geom_text(
      data = ann_pack$headers,
      mapping = ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      fontface = "bold",
      size = header_text_size
    ) +
    ggplot2::coord_cartesian(clip = "off", xlim = c(0, ann_pack$x_max)) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = grid::unit(c(0.5, right_margin_lines, 0.5, 0.5), "lines")
    )

  tree_plot
}

run_hierarchy_tree_pipeline <- function(seurat_obj,
                                        group_col = "short.label",
                                        annotation_df = NULL,
                                        annotation_tsv = NULL,
                                        assay = NULL,
                                        slot = "data",
                                        layer = NULL,
                                        out_prefix,
                                        tree_method = "hclust",
                                        dist_method = "euclidean",
                                        outgroup = NULL,
                                        drop_tip_labels = NULL,
                                        sort_labels = TRUE,
                                        tip_col = "label",
                                        tip_label_col = "display_label",
                                        annotation_cols = c("level_1", "level_2", "level_3", "meaning"),
                                        wrap_cols = "meaning",
                                        wrap_width = 28L,
                                        branch.length = "none",
                                        layout = "rectangular",
                                        width = 12,
                                        height = NULL,
                                        dpi = 300) {
  out_prefix <- scalar_chr(out_prefix, "out_prefix")
  out_dir <- dir_create_safe(dirname(out_prefix))
  out_prefix <- file.path(out_dir, basename(out_prefix))

  groups <- extract_group_labels(
    seurat_obj = seurat_obj,
    group_col = group_col,
    sort_labels = sort_labels,
    drop_tip_labels = drop_tip_labels
  )

  if (is.null(annotation_df) && is.null(annotation_tsv)) {
    template_path <- paste0(out_prefix, "_annotation_template.tsv")
    write_hierarchy_annotation_template(
      seurat_obj = seurat_obj,
      group_col = group_col,
      out_tsv = template_path,
      drop_tip_labels = drop_tip_labels,
      sort_labels = sort_labels,
      extra_cols = unique(c(tip_label_col, annotation_cols))
    )
    return(list(status = "template_written", annotation_template = template_path, labels = groups))
  }

  if (!is.null(annotation_tsv)) {
    annotation_df <- utils::read.delim(annotation_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }

  cluster_means <- compute_cluster_means(
    seurat_obj = seurat_obj,
    group_col = group_col,
    assay = assay,
    slot = slot,
    layer = layer,
    drop_tip_labels = drop_tip_labels,
    sort_labels = sort_labels
  )

  tree_pack <- build_cluster_phylo(
    cluster_means = cluster_means,
    tree_method = tree_method,
    dist_method = dist_method,
    outgroup = outgroup
  )

  matched_annotation <- prepare_annotation_table(
    annotation_df = annotation_df,
    tip_order = tree_pack$phy$tip.label,
    tip_col = tip_col,
    tip_label_col = tip_label_col,
    annotation_cols = annotation_cols,
    wrap_cols = wrap_cols,
    wrap_width = wrap_width
  )

  tree_plot <- plot_cluster_hierarchy_tree(
    phy = tree_pack$phy,
    annotation_df = matched_annotation,
    tip_label_col = tip_label_col,
    annotation_cols = annotation_cols,
    branch.length = branch.length,
    layout = layout
  )

  height <- height %||% max(7, 1.5 + nrow(matched_annotation) * 0.35)
  pdf_path <- paste0(out_prefix, ".pdf")
  png_path <- paste0(out_prefix, ".png")
  newick_path <- paste0(out_prefix, ".newick")
  means_path <- paste0(out_prefix, "_cluster_means.rds")
  ann_path <- paste0(out_prefix, "_matched_annotations.tsv")

  ggplot2::ggsave(filename = pdf_path, plot = tree_plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(filename = png_path, plot = tree_plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  ape::write.tree(phy = tree_pack$phy, file = newick_path)
  saveRDS(cluster_means, file = means_path)
  write_tsv_safe(data.frame(label = rownames(matched_annotation), matched_annotation, row.names = NULL, check.names = FALSE), ann_path)

  list(
    status = "ok",
    plot = tree_plot,
    phy = tree_pack$phy,
    dist = tree_pack$dist,
    cluster_means = cluster_means,
    matched_annotation = matched_annotation,
    files = list(
      pdf = pdf_path,
      png = png_path,
      newick = newick_path,
      cluster_means = means_path,
      annotations = ann_path
    )
  )
}

parse_cli_args <- function(args) {
  parsed <- list()
  idx <- 1L
  while (idx <= length(args)) {
    token <- args[[idx]]
    if (!startsWith(token, "--")) stop(sprintf("Unexpected argument: %s", token), call. = FALSE)
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      key <- strsplit(token, "=", fixed = TRUE)[[1L]][1L]
      value <- sub(sprintf("^%s=", key), "", token)
      parsed[[key]] <- value
      idx <- idx + 1L
      next
    }
    if (idx == length(args) || startsWith(args[[idx + 1L]], "--")) {
      parsed[[token]] <- TRUE
      idx <- idx + 1L
    } else {
      parsed[[token]] <- args[[idx + 1L]]
      idx <- idx + 2L
    }
  }
  parsed
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help) || is.null(args$seurat) || is.null(args$`group-col`) || is.null(args$`out-prefix`)) {
    cat(
      paste(
        "Usage:",
        "Rscript script/R/cluster_hierarchy_tree_helper_20260520.R \\",
        "  --seurat path/to/object.rds \\",
        "  --group-col short.label \\",
        "  --annotation-tsv path/to/cluster_annotation.tsv \\",
        "  --assay SCT \\",
        "  --slot data \\",
        "  --out-prefix output/cluster_hierarchy_tree",
        "",
        "If --annotation-tsv is omitted, the script writes an annotation template TSV and exits.",
        "This helper is lineage-agnostic and can be reused for macrophage, B/T/NK, epithelial, stromal, or any other Seurat subset.",
        sep = "\n"
      )
    )
    quit(save = "no", status = 0L)
  }

  if (!file.exists(args$seurat)) stop(sprintf("Seurat RDS not found: %s", args$seurat), call. = FALSE)
  seurat_obj <- readRDS(args$seurat)
  res <- run_hierarchy_tree_pipeline(
    seurat_obj = seurat_obj,
    group_col = args$`group-col`,
    annotation_tsv = args$`annotation-tsv`,
    assay = args$assay,
    slot = args$slot %||% "data",
    layer = args$layer,
    out_prefix = args$`out-prefix`,
    tree_method = args$`tree-method` %||% "hclust",
    dist_method = args$`dist-method` %||% "euclidean",
    outgroup = args$outgroup,
    drop_tip_labels = split_csv_arg(args$`drop-tip-labels`),
    tip_col = args$`tip-col` %||% "label",
    tip_label_col = args$`tip-label-col` %||% "display_label",
    annotation_cols = split_csv_arg(args$`annotation-cols` %||% "level_1,level_2,level_3,meaning"),
    wrap_cols = split_csv_arg(args$`wrap-cols` %||% "meaning"),
    wrap_width = as.integer(args$`wrap-width` %||% 28L),
    branch.length = args$`branch-length` %||% "none",
    layout = args$layout %||% "rectangular",
    width = as.numeric(args$width %||% 12),
    height = if (is.null(args$height)) NULL else as.numeric(args$height),
    dpi = as.integer(args$dpi %||% 300L)
  )

  if (identical(res$status, "template_written")) {
    cat(sprintf("Template written: %s\n", res$annotation_template))
  } else {
    cat(sprintf("Tree written: %s\n", res$files$pdf))
    cat(sprintf("Tree written: %s\n", res$files$png))
    cat(sprintf("Newick written: %s\n", res$files$newick))
    cat(sprintf("Cluster means written: %s\n", res$files$cluster_means))
    cat(sprintf("Matched annotations written: %s\n", res$files$annotations))
  }
}

if (sys.nframe() == 0L) main()