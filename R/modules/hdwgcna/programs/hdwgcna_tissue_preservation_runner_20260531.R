#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

OUT_DIR_DEFAULT <- file.path(
  '/home/h2048/data/R',
  format(Sys.Date(), '%Y%m%d'),
  'hdwgcna_tissue_preservation_20260531'
)

scalar_chr <- function(x, default = NA_character_) {
  if (length(x) == 0L || is.null(x) || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) return(default)
  as.character(x[[1]])
}

scalar_int <- function(x, default = NA_integer_) {
  out <- suppressWarnings(as.integer(x[[1]]))
  if (length(out) == 0L || is.na(out)) default else out
}

scalar_num <- function(x, default = NA_real_) {
  out <- suppressWarnings(as.numeric(x[[1]]))
  if (length(out) == 0L || is.na(out)) default else out
}

split_csv <- function(x) {
  x <- scalar_chr(x, '')
  if (!nzchar(x)) return(character())
  vals <- trimws(strsplit(x, ',', fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

split_int_csv <- function(x) {
  vals <- split_csv(x)
  if (length(vals) == 0L) return(integer())
  out <- suppressWarnings(as.integer(vals))
  out[is.finite(out)]
}

parse_args <- function(args) {
  out <- list(
    rds_path = NA_character_,
    lineage_name = 'lineage',
    out_dir = OUT_DIR_DEFAULT,
    celltype_col = 'cell_type_final_l3',
    tissue_col = 'tissue',
    sample_col = 'sample',
    include_celltypes = character(),
    include_tissues = character(),
    min_cells_per_tissue = 80L,
    min_samples_per_tissue = 2L,
    gene_select_mode = GENE_SELECT_MODE,
    gene_fraction = GENE_FRACTION,
    gene_n_top = GENE_N_TOP,
    metacell_k = N_METACELL_K,
    max_shared = NA_integer_,
    target_metacells = 250L,
    metacell_target_use = METACELL_TARGET,
    metacell_reduction = NA_character_,
    metacell_dims = integer(),
    metacell_min_cells = NA_integer_,
    soft_power = NA_integer_,
    soft_power_r2_cutoff = 0.85,
    soft_power_fallback = 12L,
    network_type = 'signed hybrid',
    tom_type = 'signed',
    cor_type = 'pearson',
    deep_split = 4L,
    min_module_size = 20L,
    merge_cut_height = 0.25,
    detect_cut_height = NA_real_,
    preservation_permutations = 25L,
    preservation_random_seed = 1L,
    wgcna_name = 'hdWGCNA_tissue'
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    next_arg <- function() {
      i <<- i + 1L
      if (i > length(args)) stop(sprintf('Missing value for argument: %s', arg), call. = FALSE)
      args[[i]]
    }
    if (identical(arg, '--rds-path')) out$rds_path <- normalizePath(next_arg(), winslash = '/', mustWork = FALSE)
    else if (identical(arg, '--lineage-name')) out$lineage_name <- scalar_chr(next_arg(), out$lineage_name)
    else if (identical(arg, '--out-dir')) out$out_dir <- normalizePath(next_arg(), winslash = '/', mustWork = FALSE)
    else if (identical(arg, '--celltype-col')) out$celltype_col <- scalar_chr(next_arg(), out$celltype_col)
    else if (identical(arg, '--tissue-col')) out$tissue_col <- scalar_chr(next_arg(), out$tissue_col)
    else if (identical(arg, '--sample-col')) out$sample_col <- scalar_chr(next_arg(), out$sample_col)
    else if (identical(arg, '--include-celltypes')) out$include_celltypes <- split_csv(next_arg())
    else if (identical(arg, '--include-tissues')) out$include_tissues <- split_csv(next_arg())
    else if (identical(arg, '--min-cells-per-tissue')) out$min_cells_per_tissue <- scalar_int(next_arg(), out$min_cells_per_tissue)
    else if (identical(arg, '--min-samples-per-tissue')) out$min_samples_per_tissue <- scalar_int(next_arg(), out$min_samples_per_tissue)
    else if (identical(arg, '--gene-select-mode')) out$gene_select_mode <- scalar_chr(next_arg(), out$gene_select_mode)
    else if (identical(arg, '--gene-fraction')) out$gene_fraction <- scalar_num(next_arg(), out$gene_fraction)
    else if (identical(arg, '--gene-n-top')) out$gene_n_top <- scalar_int(next_arg(), out$gene_n_top)
    else if (identical(arg, '--metacell-k')) out$metacell_k <- scalar_int(next_arg(), out$metacell_k)
    else if (identical(arg, '--max-shared')) out$max_shared <- scalar_int(next_arg(), out$max_shared)
    else if (identical(arg, '--target-metacells')) out$target_metacells <- scalar_int(next_arg(), out$target_metacells)
    else if (identical(arg, '--metacell-target-use')) out$metacell_target_use <- scalar_num(next_arg(), out$metacell_target_use)
    else if (identical(arg, '--metacell-reduction')) out$metacell_reduction <- scalar_chr(next_arg(), out$metacell_reduction)
    else if (identical(arg, '--metacell-dims')) out$metacell_dims <- split_int_csv(next_arg())
    else if (identical(arg, '--metacell-min-cells')) out$metacell_min_cells <- scalar_int(next_arg(), out$metacell_min_cells)
    else if (identical(arg, '--soft-power')) out$soft_power <- scalar_int(next_arg(), out$soft_power)
    else if (identical(arg, '--soft-power-r2-cutoff')) out$soft_power_r2_cutoff <- scalar_num(next_arg(), out$soft_power_r2_cutoff)
    else if (identical(arg, '--soft-power-fallback')) out$soft_power_fallback <- scalar_int(next_arg(), out$soft_power_fallback)
    else if (identical(arg, '--network-type')) out$network_type <- scalar_chr(next_arg(), out$network_type)
    else if (identical(arg, '--tom-type')) out$tom_type <- scalar_chr(next_arg(), out$tom_type)
    else if (identical(arg, '--cor-type')) out$cor_type <- scalar_chr(next_arg(), out$cor_type)
    else if (identical(arg, '--deep-split')) out$deep_split <- scalar_int(next_arg(), out$deep_split)
    else if (identical(arg, '--min-module-size')) out$min_module_size <- scalar_int(next_arg(), out$min_module_size)
    else if (identical(arg, '--merge-cut-height')) out$merge_cut_height <- scalar_num(next_arg(), out$merge_cut_height)
    else if (identical(arg, '--detect-cut-height')) out$detect_cut_height <- scalar_num(next_arg(), out$detect_cut_height)
    else if (identical(arg, '--preservation-permutations')) out$preservation_permutations <- scalar_int(next_arg(), out$preservation_permutations)
    else if (identical(arg, '--preservation-random-seed')) out$preservation_random_seed <- scalar_int(next_arg(), out$preservation_random_seed)
    else if (identical(arg, '--wgcna-name')) out$wgcna_name <- scalar_chr(next_arg(), out$wgcna_name)
    else if (identical(arg, '--help')) {
      cat(paste(
        'Usage:',
        '  /usr/bin/Rscript script/R/modules/hdwgcna/programs/hdwgcna_tissue_preservation_runner_20260531.R \\\n',
        '    --rds-path <lineage_seurat_rds> [--lineage-name <lineage>] [--out-dir <dir>]\\\n',
        '    [--celltype-col <col>] [--tissue-col <col>] [--sample-col <col>]\\\n',
        '    [--include-celltypes A,B] [--include-tissues nose,sinus,bronchial]\\\n',
        '    [--min-cells-per-tissue 80] [--min-samples-per-tissue 2]\\\n',
        '    [--metacell-k 25] [--max-shared 12] [--target-metacells 250]\\\n',
        '    [--soft-power-r2-cutoff 0.85] [--network-type "signed hybrid"] [--tom-type signed]\\\n',
        '    [--deep-split 4] [--min-module-size 20] [--merge-cut-height 0.25]\\\n',
        '    [--preservation-permutations 25]\\n',
        sep = ''
      ))
      quit(save = 'no', status = 0L)
    } else {
      stop(sprintf('Unknown argument: %s', arg), call. = FALSE)
    }
    i <- i + 1L
  }

  if (!is.finite(out$max_shared)) out$max_shared <- max(1L, floor(out$metacell_k / 2L))
  if (!is.finite(out$soft_power)) out$soft_power <- NA_integer_
  if (!is.finite(out$detect_cut_height)) out$detect_cut_height <- NA_real_
  if (!is.finite(out$metacell_min_cells)) out$metacell_min_cells <- out$metacell_k
  out
}

write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  utils::write.table(df, file = path, sep = '\t', quote = FALSE, row.names = FALSE, na = '')
  invisible(path)
}

prepare_tissue_grid <- function(seurat_obj,
                                celltype_col,
                                tissue_col,
                                sample_col,
                                include_celltypes = character(),
                                include_tissues = character(),
                                min_cells_per_tissue = 80L,
                                min_samples_per_tissue = 2L) {
  meta <- seurat_obj@meta.data
  required_cols <- c(celltype_col, tissue_col, sample_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0L) {
    stop(sprintf('Missing metadata columns: %s', paste(missing_cols, collapse = ', ')), call. = FALSE)
  }

  grid_df <- data.frame(
    cell_id = rownames(meta),
    celltype = as.character(meta[[celltype_col]]),
    tissue = as.character(meta[[tissue_col]]),
    sample = as.character(meta[[sample_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  grid_df <- grid_df[!is.na(grid_df$celltype) & nzchar(grid_df$celltype) &
                       !is.na(grid_df$tissue) & nzchar(grid_df$tissue) &
                       !is.na(grid_df$sample) & nzchar(grid_df$sample), , drop = FALSE]
  if (length(include_celltypes) > 0L) {
    grid_df <- grid_df[grid_df$celltype %in% include_celltypes, , drop = FALSE]
  }
  if (length(include_tissues) > 0L) {
    grid_df <- grid_df[grid_df$tissue %in% include_tissues, , drop = FALSE]
  }
  if (nrow(grid_df) == 0L) return(data.frame())

  split_key <- interaction(grid_df$celltype, grid_df$tissue, drop = TRUE, lex.order = TRUE)
  split_rows <- split(grid_df, split_key)
  rows <- lapply(split_rows, function(df) {
    data.frame(
      celltype = df$celltype[[1]],
      tissue = df$tissue[[1]],
      n_cells = nrow(df),
      n_samples = length(unique(df$sample)),
      cells = I(list(df$cell_id)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[out$n_cells >= min_cells_per_tissue & out$n_samples >= min_samples_per_tissue, , drop = FALSE]
  out[order(out$celltype, out$tissue), , drop = FALSE]
}

run_single_tissue_network <- function(seurat_obj,
                                      celltype_label,
                                      tissue_label,
                                      cells,
                                      args) {
  tissue_dir <- file.path(
    args$out_dir,
    'celltypes',
    hdwgcna_safe_id(celltype_label),
    'tissue_networks',
    hdwgcna_safe_id(tissue_label)
  )
  ensure_dir(tissue_dir)
  status_path <- file.path(tissue_dir, 'hdwgcna_tissue_network_status.json')

  sub_obj <- subset(seurat_obj, cells = cells)
  sub_obj@meta.data <- droplevels(sub_obj@meta.data)
  n_cells <- ncol(sub_obj)
  n_samples <- length(unique(as.character(sub_obj@meta.data[[args$sample_col]])))
  configured_min_cells <- if (is.na(args$metacell_min_cells)) args$metacell_k else as.integer(args$metacell_min_cells)
  max_shared_effective <- if (is.finite(args$max_shared)) as.integer(args$max_shared) else max(0L, floor(as.integer(args$metacell_k) / 2L))
  minimum_cells_for_two_metacells <- max(
    as.integer(args$metacell_k) + 1L,
    (2L * as.integer(args$metacell_k)) - max_shared_effective
  )
  effective_min_cells <- max(configured_min_cells, minimum_cells_for_two_metacells)
  sample_counts <- table(as.character(sub_obj@meta.data[[args$sample_col]]))
  valid_samples <- names(sample_counts)[sample_counts >= effective_min_cells]
  n_valid_metacell_groups <- sum(sample_counts >= effective_min_cells)
  status_rec <- list(
    lineage = args$lineage_name,
    celltype = celltype_label,
    tissue = tissue_label,
    n_cells = n_cells,
    n_samples = n_samples,
    configured_metacell_min_cells = configured_min_cells,
    n_valid_metacell_groups = n_valid_metacell_groups,
    metacell_min_cells = effective_min_cells,
    max_shared = max_shared_effective,
    output_dir = tissue_dir
  )
  hdwgcna_write_json(c(status_rec, list(status = 'started')), status_path)

  if (n_valid_metacell_groups < 2L) {
    rec <- c(
      status_rec,
      list(
        status = 'skipped',
        reason = sprintf('Need >= 2 sample groups with >= %d cells to support at least two metacells (k=%d, max_shared=%d); found %d', effective_min_cells, as.integer(args$metacell_k), max_shared_effective, n_valid_metacell_groups)
      )
    )
    hdwgcna_write_json(rec, status_path)
    rm(sub_obj)
    return(invisible(list(
      status = 'skipped',
      record = rec,
      modules = NULL,
      module_ids = character(),
      hub_genes = data.frame(),
      dat_expr = NULL,
      output_dir = tissue_dir,
      soft_power = NA_integer_
    )))
  }

  if (length(valid_samples) < length(sample_counts)) {
    keep_cells <- rownames(sub_obj@meta.data)[as.character(sub_obj@meta.data[[args$sample_col]]) %in% valid_samples]
    sub_obj <- subset(sub_obj, cells = keep_cells)
    sub_obj@meta.data <- droplevels(sub_obj@meta.data)
  }

  result <- tryCatch({
    setup_obj <- hdwgcna_global_setup(
      seurat_obj = sub_obj,
      out_dir = tissue_dir,
      group_by = unique(c(args$sample_col, args$tissue_col)),
      gene_select_mode = args$gene_select_mode,
      gene_fraction = args$gene_fraction,
      gene_n_top = args$gene_n_top,
      metacell_k = args$metacell_k,
      metacell_target = args$metacell_target_use,
      max_shared = args$max_shared,
      target_metacells = args$target_metacells,
      metacell_reduction = if (is.na(args$metacell_reduction)) NULL else args$metacell_reduction,
      metacell_dims = if (length(args$metacell_dims) == 0L) NULL else args$metacell_dims,
      metacell_min_cells = effective_min_cells,
      wgcna_name = args$wgcna_name
    )

    setdatexpr_res <- tryCatch(
      list(
        skipped = FALSE,
        value = hdwgcna_set_datexpr(
          seurat_obj = setup_obj,
          group_name = tissue_label,
          group_by = args$tissue_col,
          wgcna_name = args$wgcna_name
        )
      ),
      error = function(e) {
        err_msg <- conditionMessage(e)
        if (grepl('Too few genes with valid expression levels in the required number of samples[.]?', err_msg, ignore.case = TRUE)) {
          rec <- c(
            status_rec,
            list(
              status = 'skipped',
              reason = sprintf('SetDatExpr skipped due to insufficient valid metacell gene coverage: %s', err_msg)
            )
          )
          hdwgcna_write_json(rec, status_path)
          return(list(
            skipped = TRUE,
            value = list(
              status = 'skipped',
              record = rec,
              modules = NULL,
              module_ids = character(),
              hub_genes = data.frame(),
              dat_expr = NULL,
              output_dir = tissue_dir,
              soft_power = NA_integer_
            )
          ))
        }
        stop(e)
      }
    )
    if (isTRUE(setdatexpr_res$skipped)) {
      return(setdatexpr_res$value)
    }
    network_obj <- setdatexpr_res$value

    soft_power_res <- hdwgcna_test_soft_power(
      seurat_obj = network_obj,
      out_dir = tissue_dir,
      wgcna_name = args$wgcna_name,
      r2_cutoff = args$soft_power_r2_cutoff,
      fallback_power = args$soft_power_fallback
    )
    network_obj <- soft_power_res$obj
    soft_power_used <- if (is.finite(args$soft_power)) args$soft_power else soft_power_res$recommended_power

    network_obj <- hdwgcna_construct_network(
      seurat_obj = network_obj,
      soft_power = soft_power_used,
      net_type = args$network_type,
      tom_type = args$tom_type,
      cor_type = args$cor_type,
      deep_split = args$deep_split,
      min_module_size = args$min_module_size,
      merge_cut_height = args$merge_cut_height,
      detect_cut_height = if (is.finite(args$detect_cut_height)) args$detect_cut_height else NULL,
      wgcna_name = args$wgcna_name
    )

    modules_df <- GetModules(network_obj, wgcna_name = args$wgcna_name)
    module_ids <- setdiff(unique(as.character(modules_df$module)), 'grey')
    hdwgcna_export_modules(network_obj, out_dir = tissue_dir, wgcna_name = args$wgcna_name)

    hub_df <- data.frame()
    if (length(module_ids) > 0L) {
      network_obj <- hdwgcna_module_eigengenes(
        seurat_obj = network_obj,
        group_by = args$tissue_col,
        harmonize_by = NULL,
        wgcna_name = args$wgcna_name
      )
      hub_res <- hdwgcna_hub_genes(network_obj, wgcna_name = args$wgcna_name)
      network_obj <- hub_res$obj
      hub_df <- hub_res$hub_genes
      if (nrow(hub_df) > 0L) {
        utils::write.csv(hub_df, file.path(tissue_dir, 'hdwgcna_hub_genes.csv'), row.names = FALSE)
      }
    }

    dat_expr <- tryCatch(GetDatExpr(network_obj, wgcna_name = args$wgcna_name), error = function(e) NULL)

    rec <- c(
      status_rec,
      list(
        status = if (length(module_ids) > 0L) 'ok' else 'no_modules',
        soft_power = soft_power_used,
        n_modules = length(module_ids),
        module_ids = module_ids,
        module_membership_csv = file.path(tissue_dir, 'hdwgcna_module_membership.csv'),
        hub_genes_csv = if (nrow(hub_df) > 0L) file.path(tissue_dir, 'hdwgcna_hub_genes.csv') else NA_character_,
        soft_power_pdf = file.path(tissue_dir, 'hdwgcna_soft_power.pdf')
      )
    )
    hdwgcna_write_json(rec, status_path)

    list(
      status = rec$status,
      record = rec,
      modules = modules_df,
      module_ids = module_ids,
      hub_genes = hub_df,
      dat_expr = dat_expr,
      output_dir = tissue_dir,
      soft_power = soft_power_used
    )
  }, error = function(e) {
    rec <- c(status_rec, list(status = 'error', error = conditionMessage(e)))
    hdwgcna_write_json(rec, status_path)
    list(
      status = 'error',
      record = rec,
      modules = NULL,
      module_ids = character(),
      hub_genes = data.frame(),
      dat_expr = NULL,
      output_dir = tissue_dir,
      soft_power = NA_integer_
    )
  })

  rm(sub_obj)
  invisible(result)
}

summarize_preservation_pairs <- function(preservation_df) {
  if (!is.data.frame(preservation_df) || nrow(preservation_df) == 0L) return(data.frame())
  keep <- !preservation_df$module %in% c('grey', 'gold')
  if (all(!keep)) return(data.frame())
  df <- preservation_df[keep, , drop = FALSE]
  split_key <- interaction(df$lineage, df$celltype, df$reference_network, df$query_network, drop = TRUE)
  pieces <- split(df, split_key)
  rows <- lapply(pieces, function(x) {
    data.frame(
      lineage = x$lineage[[1]],
      celltype = x$celltype[[1]],
      reference_network = x$reference_network[[1]],
      query_network = x$query_network[[1]],
      n_modules_considered = nrow(x),
      modules_z_gt_2 = sum(x$Zsummary.pres > 2, na.rm = TRUE),
      modules_z_gt_10 = sum(x$Zsummary.pres > 10, na.rm = TRUE),
      median_Zsummary = stats::median(x$Zsummary.pres, na.rm = TRUE),
      mean_Zsummary = mean(x$Zsummary.pres, na.rm = TRUE),
      median_medianRank = if ('medianRank.pres' %in% colnames(x)) stats::median(x$medianRank.pres, na.rm = TRUE) else NA_real_,
      n_common_genes = x$n_common_genes[[1]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_preservation_heatmap <- function(pair_summary, out_dir) {
  if (!is.data.frame(pair_summary) || nrow(pair_summary) == 0L) return(invisible(NULL))
  p <- ggplot2::ggplot(pair_summary, ggplot2::aes(x = query_network, y = reference_network, fill = median_Zsummary)) +
    ggplot2::geom_tile(color = 'white') +
    ggplot2::geom_text(ggplot2::aes(label = sprintf('%d/%d', modules_z_gt_2, n_modules_considered)), size = 2.6) +
    ggplot2::facet_wrap(~celltype) +
    ggplot2::scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 2, name = 'Median Z') +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = 'bold')
    ) +
    ggplot2::labs(
      title = 'hdWGCNA tissue module preservation',
      subtitle = 'Tile label = modules with Zsummary > 2 / total reference modules',
      x = 'Query tissue',
      y = 'Reference tissue'
    )

  pdf_path <- file.path(out_dir, 'hdwgcna_tissue_preservation_heatmap.pdf')
  png_path <- file.path(out_dir, 'hdwgcna_tissue_preservation_heatmap.png')
  ggplot2::ggsave(filename = pdf_path, plot = p, width = 12, height = max(6, 2 + 2.2 * length(unique(pair_summary$celltype))), limitsize = FALSE)
  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = max(6, 2 + 2.2 * length(unique(pair_summary$celltype))), dpi = 180, limitsize = FALSE)
  invisible(list(pdf = pdf_path, png = png_path))
}

run_celltype_preservation <- function(celltype_label, tissue_networks, args) {
  ok_tissues <- names(tissue_networks)[vapply(tissue_networks, function(x) {
    identical(x$status, 'ok') && length(x$module_ids) > 0L && !is.null(x$dat_expr)
  }, logical(1))]
  if (length(ok_tissues) < 2L) {
    return(list(preservation = data.frame(), overlap = data.frame(), hub_overlap = data.frame()))
  }

  preservation_rows <- list()
  overlap_rows <- list()
  hub_overlap_rows <- list()
  p_idx <- 1L
  o_idx <- 1L
  h_idx <- 1L

  unordered_pairs <- combn(ok_tissues, 2L, simplify = FALSE)
  for (pair in unordered_pairs) {
    left <- pair[[1]]
    right <- pair[[2]]

    ov <- hdwgcna_pairwise_module_overlap(
      reference_modules = tissue_networks[[left]]$modules,
      query_modules = tissue_networks[[right]]$modules,
      reference_label = left,
      query_label = right
    )
    if (nrow(ov) > 0L) {
      ov$lineage <- args$lineage_name
      ov$celltype <- celltype_label
      overlap_rows[[o_idx]] <- ov
      o_idx <- o_idx + 1L
    }

    hov <- hdwgcna_pairwise_hub_overlap(
      reference_hubs = tissue_networks[[left]]$hub_genes,
      query_hubs = tissue_networks[[right]]$hub_genes,
      reference_label = left,
      query_label = right,
      top_n = N_HUB_GENES
    )
    if (nrow(hov) > 0L) {
      hov$lineage <- args$lineage_name
      hov$celltype <- celltype_label
      hub_overlap_rows[[h_idx]] <- hov
      h_idx <- h_idx + 1L
    }

    for (direction in list(c(left, right), c(right, left))) {
      ref_label <- direction[[1]]
      query_label <- direction[[2]]
      pres <- tryCatch(
        hdwgcna_run_module_preservation(
          reference_expr = tissue_networks[[ref_label]]$dat_expr,
          reference_modules = tissue_networks[[ref_label]]$modules,
          query_expr = tissue_networks[[query_label]]$dat_expr,
          reference_label = ref_label,
          query_label = query_label,
          network_type = args$network_type,
          n_permutations = args$preservation_permutations,
          random_seed = args$preservation_random_seed,
          verbose = 0
        ),
        error = function(e) list(summary = data.frame(), common_genes = character(), error = conditionMessage(e))
      )
      if (is.data.frame(pres$summary) && nrow(pres$summary) > 0L) {
        pres_df <- pres$summary
        pres_df$lineage <- args$lineage_name
        pres_df$celltype <- celltype_label
        pres_df$n_common_genes <- length(unique(pres$common_genes))
        preservation_rows[[p_idx]] <- pres_df
        p_idx <- p_idx + 1L
      }
    }
  }

  list(
    preservation = if (length(preservation_rows) > 0L) do.call(rbind, preservation_rows) else data.frame(),
    overlap = if (length(overlap_rows) > 0L) do.call(rbind, overlap_rows) else data.frame(),
    hub_overlap = if (length(hub_overlap_rows) > 0L) do.call(rbind, hub_overlap_rows) else data.frame()
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!file.exists(args$rds_path)) {
    stop(sprintf('Input RDS not found: %s', args$rds_path), call. = FALSE)
  }

  ensure_dir(args$out_dir)
  load_coexpr_libs()
  seurat_obj <- readRDS(args$rds_path)

  grid <- prepare_tissue_grid(
    seurat_obj = seurat_obj,
    celltype_col = args$celltype_col,
    tissue_col = args$tissue_col,
    sample_col = args$sample_col,
    include_celltypes = args$include_celltypes,
    include_tissues = args$include_tissues,
    min_cells_per_tissue = args$min_cells_per_tissue,
    min_samples_per_tissue = args$min_samples_per_tissue
  )
  if (!is.data.frame(grid) || nrow(grid) == 0L) {
    stop('No celltype × tissue groups passed the minimum filters', call. = FALSE)
  }

  manifest_rows <- list()
  all_hub_rows <- list()
  all_preservation_rows <- list()
  all_overlap_rows <- list()
  all_hub_overlap_rows <- list()
  m_idx <- 1L
  hub_idx <- 1L
  pres_idx <- 1L
  ov_idx <- 1L
  hov_idx <- 1L

  for (celltype_label in unique(grid$celltype)) {
    celltype_grid <- grid[grid$celltype == celltype_label, , drop = FALSE]
    message(sprintf('[hdWGCNA tissue] celltype=%s | tissues=%d', celltype_label, nrow(celltype_grid)))
    tissue_networks <- list()

    for (i in seq_len(nrow(celltype_grid))) {
      tissue_label <- celltype_grid$tissue[[i]]
      cells <- celltype_grid$cells[[i]]
      network_res <- run_single_tissue_network(
        seurat_obj = seurat_obj,
        celltype_label = celltype_label,
        tissue_label = tissue_label,
        cells = cells,
        args = args
      )
      tissue_networks[[tissue_label]] <- network_res
      manifest_rows[[m_idx]] <- data.frame(
        lineage = args$lineage_name,
        celltype = celltype_label,
        tissue = tissue_label,
        n_cells = celltype_grid$n_cells[[i]],
        n_samples = celltype_grid$n_samples[[i]],
        status = network_res$status,
        soft_power = network_res$soft_power %||% NA_integer_,
        n_modules = length(network_res$module_ids %||% character()),
        output_dir = network_res$output_dir,
        status_json = file.path(network_res$output_dir, 'hdwgcna_tissue_network_status.json'),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      m_idx <- m_idx + 1L

      if (is.data.frame(network_res$hub_genes) && nrow(network_res$hub_genes) > 0L) {
        hub_df <- network_res$hub_genes
        hub_df$lineage <- args$lineage_name
        hub_df$celltype <- celltype_label
        hub_df$tissue <- tissue_label
        all_hub_rows[[hub_idx]] <- hub_df
        hub_idx <- hub_idx + 1L
      }
    }

    preservation_res <- run_celltype_preservation(celltype_label, tissue_networks, args)
    if (is.data.frame(preservation_res$preservation) && nrow(preservation_res$preservation) > 0L) {
      all_preservation_rows[[pres_idx]] <- preservation_res$preservation
      pres_idx <- pres_idx + 1L
    }
    if (is.data.frame(preservation_res$overlap) && nrow(preservation_res$overlap) > 0L) {
      all_overlap_rows[[ov_idx]] <- preservation_res$overlap
      ov_idx <- ov_idx + 1L
    }
    if (is.data.frame(preservation_res$hub_overlap) && nrow(preservation_res$hub_overlap) > 0L) {
      all_hub_overlap_rows[[hov_idx]] <- preservation_res$hub_overlap
      hov_idx <- hov_idx + 1L
    }

    rm(tissue_networks)
    gc(verbose = FALSE)
  }

  manifest_df <- do.call(rbind, manifest_rows)
  preservation_df <- if (length(all_preservation_rows) > 0L) do.call(rbind, all_preservation_rows) else data.frame()
  overlap_df <- if (length(all_overlap_rows) > 0L) do.call(rbind, all_overlap_rows) else data.frame()
  hub_overlap_df <- if (length(all_hub_overlap_rows) > 0L) do.call(rbind, all_hub_overlap_rows) else data.frame()
  hub_df <- if (length(all_hub_rows) > 0L) do.call(rbind, all_hub_rows) else data.frame()
  pair_summary_df <- summarize_preservation_pairs(preservation_df)

  write_tsv(manifest_df, file.path(args$out_dir, 'hdwgcna_tissue_network_manifest.tsv'))
  if (nrow(preservation_df) > 0L) write_tsv(preservation_df, file.path(args$out_dir, 'hdwgcna_tissue_preservation_stats.tsv'))
  if (nrow(pair_summary_df) > 0L) write_tsv(pair_summary_df, file.path(args$out_dir, 'hdwgcna_tissue_preservation_pair_summary.tsv'))
  if (nrow(overlap_df) > 0L) write_tsv(overlap_df, file.path(args$out_dir, 'hdwgcna_tissue_module_overlap.tsv'))
  if (nrow(hub_overlap_df) > 0L) write_tsv(hub_overlap_df, file.path(args$out_dir, 'hdwgcna_tissue_hub_overlap.tsv'))
  if (nrow(hub_df) > 0L) write_tsv(hub_df, file.path(args$out_dir, 'hdwgcna_tissue_hub_genes.tsv'))
  plot_preservation_heatmap(pair_summary_df, args$out_dir)

  summary_obj <- list(
    lineage_name = args$lineage_name,
    input_rds = normalizePath(args$rds_path, winslash = '/', mustWork = FALSE),
    out_dir = normalizePath(args$out_dir, winslash = '/', mustWork = FALSE),
    celltype_col = args$celltype_col,
    tissue_col = args$tissue_col,
    sample_col = args$sample_col,
    n_celltype_tissue_networks = nrow(manifest_df),
    n_ok_networks = sum(manifest_df$status == 'ok'),
    n_no_module_networks = sum(manifest_df$status == 'no_modules'),
    n_skipped_networks = sum(manifest_df$status == 'skipped'),
    n_error_networks = sum(manifest_df$status == 'error'),
    n_preservation_rows = if (is.data.frame(preservation_df)) nrow(preservation_df) else 0L,
    n_pair_summary_rows = if (is.data.frame(pair_summary_df)) nrow(pair_summary_df) else 0L,
    n_module_overlap_rows = if (is.data.frame(overlap_df)) nrow(overlap_df) else 0L,
    n_hub_overlap_rows = if (is.data.frame(hub_overlap_df)) nrow(hub_overlap_df) else 0L,
    parameters = list(
      metacell_k = args$metacell_k,
      max_shared = args$max_shared,
      target_metacells = args$target_metacells,
      soft_power = args$soft_power,
      soft_power_r2_cutoff = args$soft_power_r2_cutoff,
      network_type = args$network_type,
      tom_type = args$tom_type,
      cor_type = args$cor_type,
      deep_split = args$deep_split,
      min_module_size = args$min_module_size,
      merge_cut_height = args$merge_cut_height,
      preservation_permutations = args$preservation_permutations
    )
  )
  hdwgcna_write_json(summary_obj, file.path(args$out_dir, 'hdwgcna_tissue_preservation_summary.json'))

  message(sprintf('[hdWGCNA tissue] done: %s', args$out_dir))
}

if (sys.nframe() == 0L) {
  main()
}
