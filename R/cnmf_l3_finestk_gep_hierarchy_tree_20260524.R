#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(Matrix)
  library(jsonlite)
  library(ggplot2)
  library(ggtree)
})

source('/home/h2048/script/R/cluster_hierarchy_tree_helper_20260520.R')

RUN_ROOT_DEFAULT <- '/home/h2048/output/program_full_parallel_methods_20260507'
OUTPUT_SUBDIR_DEFAULT <- 'cnmf_by_celltype'
OUT_DIR_DEFAULT <- file.path('/home/h2048/data/R', format(Sys.Date(), '%m%d'), 'cnmf_l3_finestk_gep_hierarchy_tree_20260524')
OUT_PREFIX_DEFAULT <- file.path(OUT_DIR_DEFAULT, 'cnmf_l3_finestk_gep_tree')

scalar_or_na_chr <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  out <- as.character(x)[1]
  if (is.na(out) || !nzchar(trimws(out))) return(NA_character_)
  trimws(out)
}

read_json_safe_tree <- function(path, default = list()) {
  if (!file.exists(path)) return(default)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) default)
}

pick_primary_value_tree <- function(values) {
  values <- as.character(values)
  values <- trimws(values)
  values <- values[!is.na(values) & nzchar(values) & !values %in% c('NA', 'nan', 'None', 'NULL')]
  if (length(values) == 0L) return(NA_character_)
  tab <- sort(table(values), decreasing = TRUE)
  names(tab)[1]
}

parse_cli_args_tree <- function(args) {
  parsed <- list()
  idx <- 1L
  while (idx <= length(args)) {
    token <- args[[idx]]
    if (!startsWith(token, '--')) stop(sprintf('Unexpected argument: %s', token), call. = FALSE)
    token <- substring(token, 3L)
    if (grepl('=', token, fixed = TRUE)) {
      parts <- strsplit(token, '=', fixed = TRUE)[[1L]]
      key <- parts[[1L]]
      value <- substring(token, nchar(key) + 2L)
      parsed[[key]] <- value
      idx <- idx + 1L
      next
    }
    if (idx == length(args) || startsWith(args[[idx + 1L]], '--')) {
      parsed[[token]] <- TRUE
      idx <- idx + 1L
    } else {
      parsed[[token]] <- args[[idx + 1L]]
      idx <- idx + 2L
    }
  }
  parsed
}

discover_l3_status_files <- function(run_root = RUN_ROOT_DEFAULT,
                                     output_subdir = OUTPUT_SUBDIR_DEFAULT) {
  sort(unique(Sys.glob(file.path(run_root, '*', output_subdir, 'cnmf_*', 'cnmf_full', 'cnmf_l3_status.json'))))
}

infer_l3_hierarchy_tree <- function(source_dir, status = list()) {
  l1 <- scalar_or_na_chr(status$celltype_l1 %||% status$cell_type_L1)
  l2 <- scalar_or_na_chr(status$celltype_l2 %||% status$cell_type_L2)
  l3 <- scalar_or_na_chr(status$celltype_l3 %||% status$cell_type_L3)

  meta_path <- file.path(source_dir, 'cell_metadata_for_l3_cnmf.csv')
  if ((is.na(l1) || is.na(l2) || is.na(l3)) && file.exists(meta_path)) {
    meta <- tryCatch(utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
    if (is.data.frame(meta) && nrow(meta) > 0L) {
      if (is.na(l1) && 'cell_type_L1' %in% colnames(meta)) l1 <- pick_primary_value_tree(meta$cell_type_L1)
      if (is.na(l2) && 'cell_type_L2' %in% colnames(meta)) l2 <- pick_primary_value_tree(meta$cell_type_L2)
      if (is.na(l3) && 'cell_type_L3' %in% colnames(meta)) l3 <- pick_primary_value_tree(meta$cell_type_L3)
    }
  }

  list(
    celltype_l1 = l1,
    celltype_l2 = l2,
    celltype_l3 = l3,
    metadata_path = if (file.exists(meta_path)) normalizePath(meta_path, winslash = '/', mustWork = FALSE) else NA_character_
  )
}

select_finest_k <- function(source_dir) {
  gep_dir <- file.path(source_dir, 'gep_gene_tables')
  score_files <- sort(list.files(gep_dir, pattern = '^gep_gene_scores_k[0-9]+\\.tsv$', full.names = TRUE))
  if (length(score_files) == 0L) {
    return(list(k = NA_integer_, score_file = NA_character_, all_k = integer()))
  }
  k_vals <- suppressWarnings(as.integer(sub('^.*_k([0-9]+)\\.tsv$', '\\1', score_files)))
  keep <- !is.na(k_vals)
  score_files <- score_files[keep]
  k_vals <- k_vals[keep]
  if (length(k_vals) == 0L) {
    return(list(k = NA_integer_, score_file = NA_character_, all_k = integer()))
  }
  best_idx <- which.max(k_vals)
  list(k = as.integer(k_vals[[best_idx]]), score_file = score_files[[best_idx]], all_k = sort(unique(as.integer(k_vals))))
}

find_cnmf_run_dir <- function(source_dir, run_summary = list()) {
  cnmf_output_dir <- file.path(source_dir, 'cnmf_output')
  if (!dir.exists(cnmf_output_dir)) {
    stop(sprintf('cnmf_output directory not found: %s', cnmf_output_dir), call. = FALSE)
  }
  run_name <- scalar_or_na_chr(run_summary$run_name)
  if (!is.na(run_name)) {
    candidate <- file.path(cnmf_output_dir, run_name)
    if (dir.exists(candidate)) return(normalizePath(candidate, winslash = '/', mustWork = FALSE))
  }
  subdirs <- list.dirs(cnmf_output_dir, recursive = FALSE, full.names = TRUE)
  if (length(subdirs) == 1L) return(normalizePath(subdirs[[1L]], winslash = '/', mustWork = FALSE))
  if (length(subdirs) > 1L && !is.na(run_name)) {
    hit <- subdirs[basename(subdirs) == run_name]
    if (length(hit) >= 1L) return(normalizePath(hit[[1L]], winslash = '/', mustWork = FALSE))
  }
  stop(sprintf('Unable to resolve unique cNMF run directory under: %s', cnmf_output_dir), call. = FALSE)
}

find_gep_spectra_file_tree <- function(cnmf_run_dir, k) {
  patterns <- c(
    sprintf('*.gene_spectra_score.k_%d.dt_*.txt', k),
    sprintf('*.gene_spectra_score.k_%d.*.txt', k),
    sprintf('*.spectra.k_%d.dt_*.consensus.txt', k),
    sprintf('*.spectra.k_%d.*consensus.txt', k)
  )
  for (pat in patterns) {
    hits <- sort(unique(Sys.glob(file.path(cnmf_run_dir, pat))))
    if (length(hits) > 0L) return(normalizePath(hits[[1L]], winslash = '/', mustWork = FALSE))
  }
  NA_character_
}

load_full_gep_score_matrix <- function(source_dir, k) {
  run_summary <- read_json_safe_tree(file.path(source_dir, 'run_summary.json'), default = list())
  cnmf_run_dir <- find_cnmf_run_dir(source_dir, run_summary = run_summary)
  spectra_path <- find_gep_spectra_file_tree(cnmf_run_dir, k)
  if (is.na(spectra_path) || !file.exists(spectra_path)) {
    stop(sprintf('Full spectra file not found for k=%s in %s', k, cnmf_run_dir), call. = FALSE)
  }

  if (requireNamespace('data.table', quietly = TRUE)) {
    dt <- data.table::fread(spectra_path, sep = '\t', data.table = FALSE, check.names = FALSE)
    if (ncol(dt) < 2L) {
      stop(sprintf('Spectra table has too few columns: %s', spectra_path), call. = FALSE)
    }
    row_ids <- as.character(dt[[1L]])
    mat <- as.matrix(dt[, -1L, drop = FALSE])
    rownames(mat) <- row_ids
  } else {
    df <- utils::read.delim(spectra_path, sep = '\t', row.names = 1, check.names = FALSE)
    mat <- as.matrix(df)
  }
  storage.mode(mat) <- 'double'
  mat[!is.finite(mat)] <- 0

  if (nrow(mat) == k) {
    gep_mat <- mat
  } else if (ncol(mat) == k) {
    gep_mat <- t(mat)
  } else {
    stop(sprintf('Spectra matrix shape %s is inconsistent with k=%s: %s', paste(dim(mat), collapse = 'x'), k, spectra_path), call. = FALSE)
  }

  rownames(gep_mat) <- paste0('GEP_', seq_len(k))
  colnames(gep_mat) <- as.character(colnames(gep_mat))

  list(
    gep_mat = gep_mat,
    spectra_path = normalizePath(spectra_path, winslash = '/', mustWork = FALSE),
    cnmf_run_dir = cnmf_run_dir,
    run_name = scalar_or_na_chr(run_summary$run_name)
  )
}

build_sparse_matrix_from_leaf_vectors <- function(leaf_vectors, row_labels) {
  if (length(leaf_vectors) == 0L) stop('No leaf vectors available', call. = FALSE)
  all_genes <- sort(unique(unlist(lapply(leaf_vectors, names), use.names = FALSE)))
  gene_index <- setNames(seq_along(all_genes), all_genes)

  i_list <- vector('list', length(leaf_vectors))
  j_list <- vector('list', length(leaf_vectors))
  x_list <- vector('list', length(leaf_vectors))

  for (idx in seq_along(leaf_vectors)) {
    vec <- leaf_vectors[[idx]]
    vals <- as.numeric(vec)
    genes <- names(vec)
    keep <- is.finite(vals) & vals != 0 & !is.na(genes) & nzchar(genes)
    if (!any(keep)) {
      i_list[[idx]] <- integer()
      j_list[[idx]] <- integer()
      x_list[[idx]] <- numeric()
      next
    }
    i_list[[idx]] <- rep.int(idx, sum(keep))
    j_list[[idx]] <- unname(gene_index[genes[keep]])
    x_list[[idx]] <- vals[keep]
  }

  Matrix::sparseMatrix(
    i = unlist(i_list, use.names = FALSE),
    j = unlist(j_list, use.names = FALSE),
    x = unlist(x_list, use.names = FALSE),
    dims = c(length(leaf_vectors), length(all_genes)),
    dimnames = list(row_labels, all_genes)
  )
}

compute_program_distance <- function(program_matrix, distance_method = c('cosine', 'euclidean')) {
  distance_method <- match.arg(distance_method)
  if (distance_method == 'euclidean') {
    return(stats::dist(as.matrix(program_matrix), method = 'euclidean'))
  }
  row_norms <- sqrt(as.numeric(Matrix::rowSums(program_matrix * program_matrix)))
  row_norms[!is.finite(row_norms) | row_norms <= 0] <- 1
  normed <- Matrix::Diagonal(x = 1 / row_norms) %*% program_matrix
  sim <- as.matrix(Matrix::tcrossprod(normed))
  dimnames(sim) <- list(rownames(program_matrix), rownames(program_matrix))
  sim[!is.finite(sim)] <- 0
  sim[sim > 1] <- 1
  sim[sim < -1] <- -1
  dist_mat <- 1 - sim
  dimnames(dist_mat) <- dimnames(sim)
  diag(dist_mat) <- 0
  stats::as.dist(dist_mat)
}

collect_l3_finestk_gep_inputs <- function(run_root = RUN_ROOT_DEFAULT,
                                          output_subdir = OUTPUT_SUBDIR_DEFAULT,
                                          include_lineages = character(),
                                          include_celltypes = character()) {
  status_files <- discover_l3_status_files(run_root = run_root, output_subdir = output_subdir)
  manifest_rows <- list()
  leaf_vectors <- list()
  skipped_rows <- list()

  for (idx in seq_along(status_files)) {
    status_path <- status_files[[idx]]
    if (idx == 1L || idx %% 10L == 0L || idx == length(status_files)) {
      message(sprintf('[tree] scanning source %d/%d: %s', idx, length(status_files), dirname(status_path)))
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

    gep_mat <- loaded$gep_mat
    for (row_idx in seq_len(nrow(gep_mat))) {
      gep_name <- rownames(gep_mat)[[row_idx]]
      leaf_label <- paste(lineage, safe_celltype, paste0('k', finest$k), gep_name, sep = '__')
      display_label <- paste(celltype_l3 %||% safe_celltype %||% 'Unknown_L3', gep_name, sep = '\n')
      vec <- gep_mat[row_idx, ]
      names(vec) <- colnames(gep_mat)
      keep <- is.finite(vec) & vec != 0
      leaf_vectors[[length(leaf_vectors) + 1L]] <- stats::setNames(as.numeric(vec[keep]), colnames(gep_mat)[keep])

      manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
        label = leaf_label,
        display_label = display_label,
        lineage = lineage,
        celltype_l1 = scalar_or_na_chr(hierarchy$celltype_l1),
        celltype_l2 = scalar_or_na_chr(hierarchy$celltype_l2),
        celltype_l3 = celltype_l3,
        gep = gep_name,
        k = as.integer(finest$k),
        all_k = paste(finest$all_k, collapse = ','),
        source_dir = normalizePath(source_dir, winslash = '/', mustWork = FALSE),
        spectra_file = loaded$spectra_path,
        status_json = normalizePath(status_path, winslash = '/', mustWork = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  manifest_df <- if (length(manifest_rows) > 0L) do.call(rbind, manifest_rows) else data.frame()
  skipped_df <- if (length(skipped_rows) > 0L) do.call(rbind, skipped_rows) else data.frame()
  list(manifest_df = manifest_df, leaf_vectors = leaf_vectors, skipped_df = skipped_df)
}

run_l3_finestk_gep_tree <- function(run_root = RUN_ROOT_DEFAULT,
                                    output_subdir = OUTPUT_SUBDIR_DEFAULT,
                                    out_prefix = OUT_PREFIX_DEFAULT,
                                    include_lineages = character(),
                                    include_celltypes = character(),
                                    distance_method = c('cosine', 'euclidean'),
                                    linkage_method = c('complete', 'average', 'single', 'ward.D2'),
                                    width = 22,
                                    height = NULL,
                                    dpi = 180) {
  distance_method <- match.arg(distance_method)
  linkage_method <- match.arg(linkage_method)
  out_prefix <- normalizePath(out_prefix, winslash = '/', mustWork = FALSE)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  message(sprintf('[tree] collecting L3 finest-k inputs from %s ...', run_root))

  collected <- collect_l3_finestk_gep_inputs(
    run_root = run_root,
    output_subdir = output_subdir,
    include_lineages = include_lineages,
    include_celltypes = include_celltypes
  )
  manifest_df <- collected$manifest_df
  skipped_df <- collected$skipped_df
  if (!is.data.frame(manifest_df) || nrow(manifest_df) < 2L) {
    stop('Need at least two GEP leaves to build a hierarchy tree', call. = FALSE)
  }

  message(sprintf('[tree] included %d GEP leaves across %d L3 cell types',
                  nrow(manifest_df), length(unique(na.omit(manifest_df$celltype_l3)))))
  if (is.data.frame(skipped_df) && nrow(skipped_df) > 0L) {
    message(sprintf('[tree] skipped %d inputs that were not ready/valid', nrow(skipped_df)))
  }

  message('[tree] building sparse union-gene program matrix ...')
  program_matrix <- build_sparse_matrix_from_leaf_vectors(collected$leaf_vectors, row_labels = manifest_df$label)
  message(sprintf('[tree] sparse matrix dims: %d x %d', nrow(program_matrix), ncol(program_matrix)))

  message(sprintf('[tree] computing %s distance matrix ...', distance_method))
  dist_obj <- compute_program_distance(program_matrix, distance_method = distance_method)
  message(sprintf('[tree] clustering with %s linkage ...', linkage_method))
  phy <- ape::as.phylo(stats::hclust(dist_obj, method = linkage_method))

  annotation_df <- manifest_df[, c('label', 'display_label', 'lineage', 'celltype_l2', 'celltype_l3', 'gep', 'k'), drop = FALSE]
  matched_annotation <- prepare_annotation_table(
    annotation_df = annotation_df,
    tip_order = phy$tip.label,
    tip_col = 'label',
    tip_label_col = 'display_label',
    annotation_cols = c('lineage', 'celltype_l2', 'celltype_l3', 'gep', 'k'),
    wrap_cols = c('display_label', 'celltype_l2', 'celltype_l3'),
    wrap_width = 28L
  )

  tree_plot <- plot_cluster_hierarchy_tree(
    phy = phy,
    annotation_df = matched_annotation,
    tip_label_col = 'display_label',
    annotation_cols = c('lineage', 'celltype_l2', 'celltype_l3', 'gep', 'k'),
    branch.length = 'none',
    layout = 'rectangular',
    tip_offset = 0.18,
    column_gap = 0.45,
    char_width = 0.11,
    header_offset = 1,
    tip_label_size = 1.7,
    annotation_text_size = 1.7,
    header_text_size = 2.0,
    point_size = 0.9,
    right_margin_lines = 24
  ) +
    ggplot2::labs(
      title = 'L3 cNMF finest-k GEP hierarchy tree',
      subtitle = sprintf('Leaves=%d GEPs across %d L3 cell types | distance=%s | linkage=%s',
                         nrow(matched_annotation), length(unique(na.omit(matched_annotation$celltype_l3))), distance_method, linkage_method)
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 'bold', size = 11),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  height <- height %||% max(24, 3 + nrow(matched_annotation) * 0.18)

  pdf_path <- paste0(out_prefix, '.pdf')
  png_path <- paste0(out_prefix, '.png')
  newick_path <- paste0(out_prefix, '.newick')
  manifest_path <- paste0(out_prefix, '_input_manifest.tsv')
  matched_ann_path <- paste0(out_prefix, '_matched_annotations.tsv')
  matrix_path <- paste0(out_prefix, '_sparse_matrix.rds')
  dist_path <- paste0(out_prefix, '_distance.rds')
  skipped_path <- paste0(out_prefix, '_skipped.tsv')
  summary_path <- paste0(out_prefix, '_summary.json')

  message('[tree] writing manifest / matrix / distance intermediates ...')
  utils::write.table(manifest_df, file = manifest_path, sep = '\t', quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(label = rownames(matched_annotation), matched_annotation, row.names = NULL, check.names = FALSE),
                     file = matched_ann_path, sep = '\t', quote = FALSE, row.names = FALSE)
  saveRDS(program_matrix, file = matrix_path)
  saveRDS(dist_obj, file = dist_path)
  if (is.data.frame(skipped_df) && nrow(skipped_df) > 0L) {
    utils::write.table(skipped_df, file = skipped_path, sep = '\t', quote = FALSE, row.names = FALSE)
  }

  message('[tree] rendering PDF tree ...')
  ggplot2::ggsave(filename = pdf_path, plot = tree_plot, width = width, height = height, limitsize = FALSE)
  message('[tree] rendering PNG preview ...')
  ggplot2::ggsave(filename = png_path, plot = tree_plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  message('[tree] writing Newick tree ...')
  ape::write.tree(phy = phy, file = newick_path)

  summary_obj <- list(
    run_root = normalizePath(run_root, winslash = '/', mustWork = FALSE),
    output_subdir = output_subdir,
    out_prefix = out_prefix,
    distance_method = distance_method,
    linkage_method = linkage_method,
    n_gep_leaves = nrow(manifest_df),
    n_l3_celltypes = length(unique(na.omit(manifest_df$celltype_l3))),
    n_lineages = length(unique(na.omit(manifest_df$lineage))),
    n_union_genes = ncol(program_matrix),
    finest_k_distribution = as.list(stats::setNames(as.integer(table(manifest_df$k)), names(table(manifest_df$k)))),
    lineage_distribution = as.list(stats::setNames(as.integer(table(manifest_df$lineage)), names(table(manifest_df$lineage)))),
    skipped_inputs = if (is.data.frame(skipped_df)) nrow(skipped_df) else 0L,
    files = list(
      pdf = pdf_path,
      png = png_path,
      newick = newick_path,
      input_manifest_tsv = manifest_path,
      matched_annotations_tsv = matched_ann_path,
      sparse_matrix_rds = matrix_path,
      distance_rds = dist_path,
      skipped_tsv = if (file.exists(skipped_path)) skipped_path else NULL
    )
  )
  jsonlite::write_json(summary_obj, path = summary_path, pretty = TRUE, auto_unbox = TRUE, null = 'null')
  message(sprintf('[tree] complete: %s', pdf_path))

  list(
    plot = tree_plot,
    phy = phy,
    manifest_df = manifest_df,
    matched_annotation = matched_annotation,
    skipped_df = skipped_df,
    program_matrix = program_matrix,
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
        'Rscript script/R/cnmf_l3_finestk_gep_hierarchy_tree_20260524.R',
        '  --run-root /home/h2048/output/program_full_parallel_methods_20260507',
        '  --out-prefix /home/h2048/data/R/0524/cnmf_l3_finestk_gep_hierarchy_tree_20260524/cnmf_l3_finestk_gep_tree',
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
  distance_method <- args$`distance-method` %||% 'cosine'
  linkage_method <- args$`linkage-method` %||% 'complete'
  width <- as.numeric(args$width %||% 22)
  height <- if (is.null(args$height)) NULL else as.numeric(args$height)
  dpi <- as.integer(args$dpi %||% 180L)

  res <- run_l3_finestk_gep_tree(
    run_root = run_root,
    output_subdir = output_subdir,
    out_prefix = out_prefix,
    include_lineages = include_lineages,
    include_celltypes = include_celltypes,
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
  cat(sprintf('Included finest-k GEP leaves: %d\n', nrow(res$manifest_df)))
  cat(sprintf('Included L3 cell types: %d\n', length(unique(na.omit(res$manifest_df$celltype_l3)))))
  if (is.data.frame(res$skipped_df) && nrow(res$skipped_df) > 0L) {
    cat(sprintf('Skipped inputs: %d\n', nrow(res$skipped_df)))
  }
}

if (sys.nframe() == 0L) main()
