#!/usr/bin/env Rscript

source('/home/h2048/script/R/cnmf_l3_finestk_gep_hierarchy_tree_20260524.R')

RUN_ROOT_DEFAULT <- '/home/h2048/output/program_full_parallel_methods_20260507'
OUTPUT_SUBDIR_DEFAULT <- 'cnmf_by_celltype'
OUT_DIR_DEFAULT <- file.path('/home/h2048/data/R', format(Sys.Date(), '%m%d'), 'cnmf_l3_finestk_celltype_hierarchy_tree_20260531')
OUT_PREFIX_DEFAULT <- file.path(OUT_DIR_DEFAULT, 'cnmf_l3_finestk_celltype_tree')

aggregate_gep_to_celltype_profile <- function(gep_mat,
                                              aggregation_method = c('mean', 'max')) {
  aggregation_method <- match.arg(aggregation_method)
  if (is.null(gep_mat) || !is.matrix(gep_mat) || nrow(gep_mat) == 0L) {
    stop('gep_mat must be a non-empty matrix', call. = FALSE)
  }
  if (aggregation_method == 'mean') {
    profile <- colMeans(gep_mat)
  } else {
    profile <- apply(gep_mat, 2L, max)
  }
  profile <- as.numeric(profile)
  names(profile) <- colnames(gep_mat)
  profile[!is.finite(profile)] <- 0
  profile
}

collect_l3_finestk_celltype_inputs <- function(run_root = RUN_ROOT_DEFAULT,
                                               output_subdir = OUTPUT_SUBDIR_DEFAULT,
                                               include_lineages = character(),
                                               include_celltypes = character(),
                                               aggregation_method = c('mean', 'max')) {
  aggregation_method <- match.arg(aggregation_method)
  status_files <- discover_l3_status_files(run_root = run_root, output_subdir = output_subdir)
  manifest_rows <- list()
  profile_vectors <- list()
  skipped_rows <- list()

  for (idx in seq_along(status_files)) {
    status_path <- status_files[[idx]]
    if (idx == 1L || idx %% 10L == 0L || idx == length(status_files)) {
      message(sprintf('[celltype-tree] scanning source %d/%d: %s', idx, length(status_files), dirname(status_path)))
    }

    status <- read_json_safe_tree(status_path, default = list())
    source_dir <- dirname(status_path)
    lineage <- scalar_or_na_chr(status$lineage)
    if (is.na(lineage)) {
      lineage <- basename(dirname(dirname(dirname(source_dir))))
    }
    if (length(include_lineages) > 0L && !lineage %in% include_lineages) next

    hierarchy <- infer_l3_hierarchy_tree(source_dir, status = status)
    celltype_l3 <- scalar_or_na_chr(hierarchy$celltype_l3 %||% status$celltype_l3)
    safe_celltype <- scalar_or_na_chr(status$safe_celltype)
    if (is.na(safe_celltype) && !is.na(celltype_l3)) safe_celltype <- safe_unit_id(celltype_l3)
    if (length(include_celltypes) > 0L && !any(c(celltype_l3, safe_celltype) %in% include_celltypes)) next

    if (!identical(scalar_or_na_chr(status$status), 'ok')) {
      skipped_rows[[length(skipped_rows) + 1L]] <- data.frame(
        lineage = lineage,
        celltype_l3 = celltype_l3,
        safe_celltype = safe_celltype,
        reason = paste0('status=', scalar_or_na_chr(status$status)),
        status_json = normalizePath(status_path, winslash = '/', mustWork = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      next
    }

    finest <- select_finest_k(source_dir)
    if (is.na(finest$k)) {
      skipped_rows[[length(skipped_rows) + 1L]] <- data.frame(
        lineage = lineage,
        celltype_l3 = celltype_l3,
        safe_celltype = safe_celltype,
        reason = 'missing_finest_k_score_file',
        status_json = normalizePath(status_path, winslash = '/', mustWork = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      next
    }

    loaded <- tryCatch(load_full_gep_score_matrix(source_dir, finest$k), error = identity)
    if (inherits(loaded, 'error')) {
      skipped_rows[[length(skipped_rows) + 1L]] <- data.frame(
        lineage = lineage,
        celltype_l3 = celltype_l3,
        safe_celltype = safe_celltype,
        reason = paste0('load_error: ', conditionMessage(loaded)),
        status_json = normalizePath(status_path, winslash = '/', mustWork = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      next
    }

    profile <- aggregate_gep_to_celltype_profile(
      gep_mat = loaded$gep_mat,
      aggregation_method = aggregation_method
    )
    keep <- is.finite(profile) & profile != 0
    profile_vectors[[length(profile_vectors) + 1L]] <- stats::setNames(as.numeric(profile[keep]), names(profile)[keep])

    label <- paste(lineage, safe_celltype, paste0('k', finest$k), sep = '__')
    display_label <- celltype_l3 %||% safe_celltype %||% 'Unknown_L3'
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      label = label,
      display_label = display_label,
      lineage = lineage,
      celltype_l1 = scalar_or_na_chr(hierarchy$celltype_l1),
      celltype_l2 = scalar_or_na_chr(hierarchy$celltype_l2),
      celltype_l3 = celltype_l3,
      k = as.integer(finest$k),
      n_geps_collapsed = nrow(loaded$gep_mat),
      aggregation_method = aggregation_method,
      all_k = paste(finest$all_k, collapse = ','),
      source_dir = normalizePath(source_dir, winslash = '/', mustWork = FALSE),
      spectra_file = loaded$spectra_path,
      status_json = normalizePath(status_path, winslash = '/', mustWork = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  manifest_df <- if (length(manifest_rows) > 0L) do.call(rbind, manifest_rows) else data.frame()
  skipped_df <- if (length(skipped_rows) > 0L) do.call(rbind, skipped_rows) else data.frame()
  list(manifest_df = manifest_df, profile_vectors = profile_vectors, skipped_df = skipped_df)
}

run_l3_finestk_celltype_tree <- function(run_root = RUN_ROOT_DEFAULT,
                                         output_subdir = OUTPUT_SUBDIR_DEFAULT,
                                         out_prefix = OUT_PREFIX_DEFAULT,
                                         include_lineages = character(),
                                         include_celltypes = character(),
                                         aggregation_method = c('mean', 'max'),
                                         distance_method = c('cosine', 'euclidean'),
                                         linkage_method = c('complete', 'average', 'single', 'ward.D2'),
                                         width = 18,
                                         height = NULL,
                                         dpi = 180) {
  aggregation_method <- match.arg(aggregation_method)
  distance_method <- match.arg(distance_method)
  linkage_method <- match.arg(linkage_method)
  out_prefix <- normalizePath(out_prefix, winslash = '/', mustWork = FALSE)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  message(sprintf('[celltype-tree] collecting L3 finest-k celltype profiles from %s ...', run_root))
  collected <- collect_l3_finestk_celltype_inputs(
    run_root = run_root,
    output_subdir = output_subdir,
    include_lineages = include_lineages,
    include_celltypes = include_celltypes,
    aggregation_method = aggregation_method
  )
  manifest_df <- collected$manifest_df
  skipped_df <- collected$skipped_df
  if (!is.data.frame(manifest_df) || nrow(manifest_df) < 2L) {
    stop('Need at least two cell types to build a hierarchy tree', call. = FALSE)
  }

  message(sprintf('[celltype-tree] included %d celltype leaves across %d lineages',
                  nrow(manifest_df), length(unique(na.omit(manifest_df$lineage)))))
  if (is.data.frame(skipped_df) && nrow(skipped_df) > 0L) {
    message(sprintf('[celltype-tree] skipped %d inputs that were not ready/valid', nrow(skipped_df)))
  }

  message('[celltype-tree] building sparse union-gene celltype matrix ...')
  profile_matrix <- build_sparse_matrix_from_leaf_vectors(
    leaf_vectors = collected$profile_vectors,
    row_labels = manifest_df$label
  )
  message(sprintf('[celltype-tree] sparse matrix dims: %d x %d', nrow(profile_matrix), ncol(profile_matrix)))

  message(sprintf('[celltype-tree] computing %s distance matrix ...', distance_method))
  dist_obj <- compute_program_distance(profile_matrix, distance_method = distance_method)
  message(sprintf('[celltype-tree] clustering with %s linkage ...', linkage_method))
  phy <- ape::as.phylo(stats::hclust(dist_obj, method = linkage_method))

  annotation_df <- manifest_df[, c('label', 'display_label', 'lineage', 'celltype_l2', 'k', 'n_geps_collapsed', 'aggregation_method'), drop = FALSE]
  matched_annotation <- prepare_annotation_table(
    annotation_df = annotation_df,
    tip_order = phy$tip.label,
    tip_col = 'label',
    tip_label_col = 'display_label',
    annotation_cols = c('lineage', 'celltype_l2', 'k', 'n_geps_collapsed', 'aggregation_method'),
    wrap_cols = c('display_label', 'celltype_l2'),
    wrap_width = 28L
  )

  tree_plot <- plot_cluster_hierarchy_tree(
    phy = phy,
    annotation_df = matched_annotation,
    tip_label_col = 'display_label',
    annotation_cols = c('lineage', 'celltype_l2', 'k', 'n_geps_collapsed', 'aggregation_method'),
    branch.length = 'none',
    layout = 'rectangular',
    tip_offset = 0.18,
    column_gap = 0.45,
    char_width = 0.11,
    header_offset = 1,
    tip_label_size = 2.0,
    annotation_text_size = 1.9,
    header_text_size = 2.2,
    point_size = 1.1,
    right_margin_lines = 22
  ) +
    ggplot2::labs(
      title = 'L3 cNMF finest-k celltype hierarchy tree',
      subtitle = sprintf('Leaves=%d L3 celltypes | aggregation=%s | distance=%s | linkage=%s',
                         nrow(matched_annotation), aggregation_method, distance_method, linkage_method)
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'bold', size = 11),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  height <- height %||% max(10, 3 + nrow(matched_annotation) * 0.24)

  pdf_path <- paste0(out_prefix, '.pdf')
  png_path <- paste0(out_prefix, '.png')
  newick_path <- paste0(out_prefix, '.newick')
  manifest_path <- paste0(out_prefix, '_input_manifest.tsv')
  matched_ann_path <- paste0(out_prefix, '_matched_annotations.tsv')
  matrix_path <- paste0(out_prefix, '_profile_matrix.rds')
  dist_path <- paste0(out_prefix, '_distance.rds')
  skipped_path <- paste0(out_prefix, '_skipped.tsv')
  summary_path <- paste0(out_prefix, '_summary.json')

  message('[celltype-tree] writing manifest / matrix / distance intermediates ...')
  utils::write.table(manifest_df, file = manifest_path, sep = '\t', quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(label = rownames(matched_annotation), matched_annotation, row.names = NULL, check.names = FALSE),
                     file = matched_ann_path, sep = '\t', quote = FALSE, row.names = FALSE)
  saveRDS(profile_matrix, file = matrix_path)
  saveRDS(dist_obj, file = dist_path)
  if (is.data.frame(skipped_df) && nrow(skipped_df) > 0L) {
    utils::write.table(skipped_df, file = skipped_path, sep = '\t', quote = FALSE, row.names = FALSE)
  }

  message('[celltype-tree] rendering PDF tree ...')
  ggplot2::ggsave(filename = pdf_path, plot = tree_plot, width = width, height = height, limitsize = FALSE)
  message('[celltype-tree] rendering PNG preview ...')
  ggplot2::ggsave(filename = png_path, plot = tree_plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  message('[celltype-tree] writing Newick tree ...')
  ape::write.tree(phy = phy, file = newick_path)

  summary_obj <- list(
    run_root = normalizePath(run_root, winslash = '/', mustWork = FALSE),
    output_subdir = output_subdir,
    out_prefix = out_prefix,
    aggregation_method = aggregation_method,
    distance_method = distance_method,
    linkage_method = linkage_method,
    n_celltype_leaves = nrow(manifest_df),
    n_lineages = length(unique(na.omit(manifest_df$lineage))),
    n_union_genes = ncol(profile_matrix),
    finest_k_distribution = as.list(stats::setNames(as.integer(table(manifest_df$k)), names(table(manifest_df$k)))),
    lineage_distribution = as.list(stats::setNames(as.integer(table(manifest_df$lineage)), names(table(manifest_df$lineage)))),
    n_geps_collapsed_summary = list(
      min = min(manifest_df$n_geps_collapsed),
      median = as.numeric(stats::median(manifest_df$n_geps_collapsed)),
      max = max(manifest_df$n_geps_collapsed)
    ),
    skipped_inputs = if (is.data.frame(skipped_df)) nrow(skipped_df) else 0L,
    files = list(
      pdf = pdf_path,
      png = png_path,
      newick = newick_path,
      input_manifest_tsv = manifest_path,
      matched_annotations_tsv = matched_ann_path,
      profile_matrix_rds = matrix_path,
      distance_rds = dist_path,
      skipped_tsv = if (file.exists(skipped_path)) skipped_path else NULL
    )
  )
  jsonlite::write_json(summary_obj, path = summary_path, pretty = TRUE, auto_unbox = TRUE, null = 'null')
  message(sprintf('[celltype-tree] complete: %s', pdf_path))

  list(
    plot = tree_plot,
    phy = phy,
    manifest_df = manifest_df,
    matched_annotation = matched_annotation,
    skipped_df = skipped_df,
    profile_matrix = profile_matrix,
    dist = dist_obj,
    files = c(summary_obj$files, list(summary_json = summary_path))
  )
}

main <- function() {
  args <- parse_cli_args_tree(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    cat(
      paste(
        'Usage:',
        'Rscript script/R/modules/cnmf/programs/cnmf_l3_finestk_celltype_hierarchy_tree_20260531.R',
        '  --run-root /home/h2048/output/program_full_parallel_methods_20260507',
        '  --out-prefix /home/h2048/data/R/0531/cnmf_l3_finestk_celltype_hierarchy_tree_20260531/cnmf_l3_finestk_celltype_tree',
        '  --aggregation-method mean',
        '  --distance-method cosine',
        '  --linkage-method complete',
        '',
        'Optional filters:',
        '  --include-lineages bcell,tnk',
        '  --include-celltypes Memory_B,AT2',
        sep = '\n'
      )
    )
    quit(save = 'no', status = 0L)
  }

  run_root <- args$`run-root` %||% RUN_ROOT_DEFAULT
  output_subdir <- args$`output-subdir` %||% OUTPUT_SUBDIR_DEFAULT
  out_prefix <- args$`out-prefix` %||% OUT_PREFIX_DEFAULT
  include_lineages <- split_csv_arg(args$`include-lineages` %||% '')
  include_celltypes <- split_csv_arg(args$`include-celltypes` %||% '')
  aggregation_method <- args$`aggregation-method` %||% 'mean'
  distance_method <- args$`distance-method` %||% 'cosine'
  linkage_method <- args$`linkage-method` %||% 'complete'
  width <- as.numeric(args$width %||% 18)
  height <- if (is.null(args$height)) NULL else as.numeric(args$height)
  dpi <- as.integer(args$dpi %||% 180L)

  res <- run_l3_finestk_celltype_tree(
    run_root = run_root,
    output_subdir = output_subdir,
    out_prefix = out_prefix,
    include_lineages = include_lineages,
    include_celltypes = include_celltypes,
    aggregation_method = aggregation_method,
    distance_method = distance_method,
    linkage_method = linkage_method,
    width = width,
    height = height,
    dpi = dpi
  )

  cat(sprintf('Tree written: %s\n', res$files$pdf))
  cat(sprintf('Preview written: %s\n', res$files$png))
  cat(sprintf('Newick written: %s\n', res$files$newick))
  cat(sprintf('Manifest written: %s\n', res$files$input_manifest_tsv))
  cat(sprintf('Summary written: %s\n', res$files$summary_json))
  cat(sprintf('Included celltype leaves: %d\n', nrow(res$manifest_df)))
  if (is.data.frame(res$skipped_df) && nrow(res$skipped_df) > 0L) {
    cat(sprintf('Skipped inputs: %d\n', nrow(res$skipped_df)))
  }
}

if (sys.nframe() == 0L) main()
