#!/usr/bin/env Rscript

HELPER_PATH <- "/home/h2048/script/R/cell_communication_helper_20260506_v1_0.R"
source(HELPER_PATH)

DEFAULT_CELLCHAT_DIR <- "/home/h2048/data/R/0525/non_unified_airway/communication/cellchat_lineage_native_export_v2"
DEFAULT_LIANA_DIR <- "/home/h2048/data/py/20260527/non_unified_airway/communication_lineage_native_export_v3/liana_20260527_v3"
DEFAULT_CELLPHONEDB_DIR <- "/home/h2048/data/py/20260527/non_unified_airway/communication_lineage_native_export_v3/cellphonedb_20260527_v3"
DEFAULT_CONSENSUS_DIR <- "/home/h2048/data/py/20260527/non_unified_airway/communication_lineage_native_export_v3/consensus_20260527_v3"
DEFAULT_OUTPUT_DIR <- "/home/h2048/data/R/20260530/non_unified_airway/communication/chord_panels_lineage_native_export_v3_20260530_v1"

RIGHT_HIGHER_COLOR <- grDevices::adjustcolor("#D55E00", alpha.f = 0.72)
LEFT_HIGHER_COLOR <- grDevices::adjustcolor("#2C7BB6", alpha.f = 0.72)
NEUTRAL_COLOR <- grDevices::adjustcolor("#7F7F7F", alpha.f = 0.55)

parse_bool_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || !nzchar(trimws(x[[1]]))) return(default)
  value <- tolower(trimws(as.character(x[[1]])))
  if (value %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  default
}

parse_args_simple <- function(args) {
  out <- list(
    cellchat_dir = DEFAULT_CELLCHAT_DIR,
    liana_dir = DEFAULT_LIANA_DIR,
    cellphonedb_dir = DEFAULT_CELLPHONEDB_DIR,
    consensus_dir = DEFAULT_CONSENSUS_DIR,
    output_dir = DEFAULT_OUTPUT_DIR,
    top_n = 24L,
    overwrite = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--cellchat-dir" && i < length(args)) {
      out$cellchat_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--liana-dir" && i < length(args)) {
      out$liana_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--cellphonedb-dir" && i < length(args)) {
      out$cellphonedb_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--consensus-dir" && i < length(args)) {
      out$consensus_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--output-dir" && i < length(args)) {
      out$output_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--top-n" && i < length(args)) {
      out$top_n <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--overwrite" && i < length(args)) {
      out$overwrite <- parse_bool_flag(args[[i + 1L]], default = TRUE)
      i <- i + 2L
    } else if (key == "--overwrite") {
      out$overwrite <- TRUE
      i <- i + 1L
    } else {
      stop(sprintf("Unknown or incomplete argument: %s", key), call. = FALSE)
    }
  }
  out
}

read_tsv_auto <- function(path) {
  if (!file.exists(path)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE))
  }
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file = gzfile(path), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

normalize_path_safe <- function(path) {
  if (exists("cci_path_normalize", mode = "function")) {
    return(cci_path_normalize(path))
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_head_df <- function(df, n) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  utils::head(df, max(0L, as.integer(n)))
}

list_run_keys_from_base <- function(base_dir) {
  if (!dir.exists(base_dir)) return(data.frame())
  contrasts <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  rows <- list()
  for (contrast_id in contrasts) {
    contrast_path <- file.path(base_dir, contrast_id)
    levels <- list.dirs(contrast_path, recursive = FALSE, full.names = FALSE)
    if (length(levels) == 0L) next
    for (level in levels) {
      rows[[length(rows) + 1L]] <- data.frame(
        contrast_id = contrast_id,
        level = level,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  unique(out)
}

collect_run_keys <- function(cellchat_dir, liana_dir, cellphonedb_dir) {
  parts <- list(
    list_run_keys_from_base(cellchat_dir),
    list_run_keys_from_base(liana_dir),
    list_run_keys_from_base(cellphonedb_dir)
  )
  parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, parts)
  if (length(parts) == 0L) return(data.frame())
  out <- do.call(rbind, parts)
  out <- unique(out[, c("contrast_id", "level"), drop = FALSE])
  out <- out[order(out$contrast_id, out$level), , drop = FALSE]
  rownames(out) <- NULL
  out
}

resolve_liana_weight_col <- function(df) {
  candidates <- c("scaled_weight", "lr_means", "expr_prod", "lr_logfc")
  found <- candidates[candidates %in% colnames(df)]
  if (length(found) == 0L) return(NULL)
  found[[1L]]
}

summarise_sender_receiver_edges <- function(df,
                                            sender_col,
                                            receiver_col,
                                            weight_col,
                                            top_n = 24L,
                                            valid_labels = NULL,
                                            aggregate_abs = FALSE,
                                            dominant_side_col = NULL,
                                            use_delta_colors = FALSE) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(data.frame())
  required_cols <- c(sender_col, receiver_col, weight_col)
  if (!all(required_cols %in% colnames(df))) return(data.frame())

  out <- data.frame(
    sender = trimws(as.character(df[[sender_col]])),
    receiver = trimws(as.character(df[[receiver_col]])),
    raw_weight = safe_numeric(df[[weight_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(dominant_side_col) && dominant_side_col %in% colnames(df)) {
    out$dominant_side_raw <- as.character(df[[dominant_side_col]])
  }
  if (!is.null(valid_labels)) {
    valid_labels <- unique(trimws(as.character(valid_labels)))
    valid_labels <- valid_labels[nzchar(valid_labels)]
    out <- out[out$sender %in% valid_labels & out$receiver %in% valid_labels, , drop = FALSE]
  }
  out <- out[nzchar(out$sender) & nzchar(out$receiver) & is.finite(out$raw_weight), , drop = FALSE]
  if (nrow(out) == 0L) return(data.frame())

  out$edge_weight <- if (isTRUE(aggregate_abs)) abs(out$raw_weight) else out$raw_weight
  out <- out[is.finite(out$edge_weight) & out$edge_weight > 0, , drop = FALSE]
  if (nrow(out) == 0L) return(data.frame())

  agg_weight <- aggregate(edge_weight ~ sender + receiver, data = out, FUN = sum, na.rm = TRUE)
  agg_signed <- aggregate(raw_weight ~ sender + receiver, data = out, FUN = sum, na.rm = TRUE)
  agg_counts <- aggregate(edge_weight ~ sender + receiver, data = out, FUN = length)
  colnames(agg_counts)[colnames(agg_counts) == "edge_weight"] <- "n_source_rows"

  merged <- merge(agg_weight, agg_signed, by = c("sender", "receiver"), all = TRUE)
  merged <- merge(merged, agg_counts, by = c("sender", "receiver"), all = TRUE)
  colnames(merged)[colnames(merged) == "raw_weight"] <- "signed_weight"
  if (isTRUE(use_delta_colors)) {
    merged$dominant_side <- ifelse(
      merged$signed_weight > 0,
      "right_higher",
      ifelse(merged$signed_weight < 0, "left_higher", "neutral")
    )
    merged$link_color <- ifelse(
      merged$dominant_side == "right_higher",
      RIGHT_HIGHER_COLOR,
      ifelse(merged$dominant_side == "left_higher", LEFT_HIGHER_COLOR, NEUTRAL_COLOR)
    )
  }
  merged <- merged[order(merged$edge_weight, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  merged <- safe_head_df(merged, top_n)
  merged$edge_rank <- seq_len(nrow(merged))
  merged$edge_label <- paste(merged$sender, "->", merged$receiver)
  rownames(merged) <- NULL
  merged
}

plot_sender_receiver_chord_from_edges <- function(edge_df,
                                                  output_prefix,
                                                  title,
                                                  subtitle = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize is not available; chord plot skipped.")
    return(NULL)
  }
  if (is.null(edge_df) || !is.data.frame(edge_df) || nrow(edge_df) == 0L) return(NULL)

  chord_df <- edge_df[, c("sender", "receiver", "edge_weight"), drop = FALSE]
  sectors <- sort(unique(c(chord_df$sender, chord_df$receiver)))
  if (length(sectors) < 2L) return(NULL)
  grid_col <- grDevices::hcl(
    h = seq(15, 375, length.out = length(sectors) + 1L)[seq_along(sectors)],
    c = 70,
    l = 65
  )
  names(grid_col) <- sectors
  link_cols <- if ("link_color" %in% colnames(edge_df)) {
    edge_df$link_color
  } else {
    grDevices::adjustcolor(grid_col[chord_df$sender], alpha.f = 0.55)
  }

  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)
  draw_once <- function(device_fun) {
    device_opened <- FALSE
    tryCatch({
      device_fun()
      device_opened <- TRUE
      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op), add = TRUE)
      circlize::circos.clear()
      circlize::circos.par(
        start.degree = 90,
        gap.after = rep(if (length(sectors) > 18L) 1.1 else 2, length(sectors)),
        track.margin = c(0.01, 0.01),
        canvas.xlim = c(-1.4, 1.4),
        canvas.ylim = c(-1.4, 1.4),
        points.overflow.warning = FALSE
      )
      circlize::chordDiagram(
        x = chord_df,
        grid.col = grid_col,
        col = link_cols,
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        diffHeight = -0.04,
        link.arr.type = "big.arrow",
        link.arr.length = 0.12,
        link.sort = TRUE,
        link.largest.ontop = TRUE,
        transparency = 0,
        reduce = 0,
        annotationTrack = "grid",
        preAllocateTracks = list(track.height = 0.18)
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          sector_name <- circlize::get.cell.meta.data("sector.index")
          xlim <- circlize::get.cell.meta.data("xlim")
          ylim <- circlize::get.cell.meta.data("ylim")
          circlize::circos.text(
            x = mean(xlim),
            y = ylim[[1]],
            labels = sector_name,
            facing = "clockwise",
            niceFacing = TRUE,
            adj = c(0, 0.5),
            cex = if (length(sectors) > 18L) 0.36 else if (length(sectors) > 12L) 0.44 else 0.52
          )
        }
      )
      graphics::title(main = title, cex.main = 1.02)
      if (!is.null(subtitle) && nzchar(subtitle)) {
        graphics::mtext(subtitle, side = 1, line = -1, cex = 0.72)
      }
      if ("dominant_side" %in% colnames(edge_df) && any(edge_df$dominant_side %in% c("left_higher", "right_higher"))) {
        graphics::legend(
          "topleft",
          legend = c("right higher", "left higher"),
          fill = c(RIGHT_HIGHER_COLOR, LEFT_HIGHER_COLOR),
          border = NA,
          bty = "n",
          cex = 0.8
        )
      }
      grDevices::dev.off()
      device_opened <- FALSE
    }, finally = {
      try(circlize::circos.clear(), silent = TRUE)
      if (isTRUE(device_opened)) try(grDevices::dev.off(), silent = TRUE)
    })
  }

  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  draw_once(function() grDevices::pdf(pdf_path, width = 12, height = 12, onefile = FALSE))
  draw_once(function() grDevices::png(png_path, width = 3200, height = 3200, res = 240))
  list(
    pdf = normalize_path_safe(pdf_path),
    png = normalize_path_safe(png_path)
  )
}

load_cellchat_side_edges <- function(path, top_n) {
  df <- read_tsv_auto(path)
  sender_col <- if ("sender" %in% colnames(df)) "sender" else if ("source" %in% colnames(df)) "source" else NULL
  receiver_col <- if ("receiver" %in% colnames(df)) "receiver" else if ("target" %in% colnames(df)) "target" else NULL
  if (is.null(sender_col) || is.null(receiver_col)) return(data.frame())
  summarise_sender_receiver_edges(df, sender_col = sender_col, receiver_col = receiver_col, weight_col = "prob", top_n = top_n, aggregate_abs = FALSE)
}

load_liana_side_edges <- function(path, top_n) {
  df <- read_tsv_auto(path)
  weight_col <- resolve_liana_weight_col(df)
  if (is.null(weight_col)) return(data.frame())
  edge_df <- summarise_sender_receiver_edges(df, sender_col = "sender", receiver_col = "receiver", weight_col = weight_col, top_n = top_n, aggregate_abs = FALSE)
  if (nrow(edge_df) > 0L) edge_df$weight_metric <- weight_col
  edge_df
}

load_cellphonedb_side_edges <- function(side_dir, top_n) {
  means_path <- file.path(side_dir, "means_long.tsv.gz")
  meta_path <- file.path(side_dir, "meta.tsv")
  means_df <- read_tsv_auto(means_path)
  meta_df <- read_tsv_auto(meta_path)
  valid_labels <- if (is.data.frame(meta_df) && "cell_type" %in% colnames(meta_df)) unique(as.character(meta_df$cell_type)) else NULL
  summarise_sender_receiver_edges(means_df, sender_col = "sender", receiver_col = "receiver", weight_col = "mean", top_n = top_n, valid_labels = valid_labels, aggregate_abs = FALSE)
}

load_delta_edges <- function(path, top_n, weight_col = "delta_score_right_minus_left") {
  df <- read_tsv_auto(path)
  summarise_sender_receiver_edges(df, sender_col = "sender", receiver_col = "receiver", weight_col = weight_col, top_n = top_n, aggregate_abs = TRUE, use_delta_colors = TRUE)
}

load_consensus_edges <- function(consensus_sender_receiver_df, contrast_id, level, top_n) {
  if (is.null(consensus_sender_receiver_df) || !is.data.frame(consensus_sender_receiver_df) || nrow(consensus_sender_receiver_df) == 0L) return(data.frame())
  sub <- consensus_sender_receiver_df[
    as.character(consensus_sender_receiver_df$contrast_id) %in% contrast_id &
      as.character(consensus_sender_receiver_df$level) %in% level,
    ,
    drop = FALSE
  ]
  if (nrow(sub) == 0L) return(data.frame())
  summarise_sender_receiver_edges(sub, sender_col = "sender", receiver_col = "receiver", weight_col = "sum_abs_mean_delta", top_n = top_n, aggregate_abs = FALSE)
}

render_chord_bundle <- function(edge_df,
                                output_dir,
                                contrast_id,
                                level,
                                method,
                                scope,
                                title,
                                source_path,
                                weight_metric,
                                overwrite = FALSE) {
  fig_dir <- file.path(output_dir, contrast_id, level, "figures")
  table_dir <- file.path(output_dir, contrast_id, level, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- file.path(fig_dir, sprintf("%s_%s_chord", method, scope))
  data_tsv <- file.path(table_dir, sprintf("%s_%s_chord_edges.tsv.gz", method, scope))
  data_csv <- file.path(table_dir, sprintf("%s_%s_chord_edges.csv", method, scope))
  pdf_path <- sprintf("%s.pdf", prefix)
  png_path <- sprintf("%s.png", prefix)

  if (!isTRUE(overwrite) && file.exists(pdf_path) && file.exists(png_path) && file.exists(data_tsv)) {
    return(data.frame(
      contrast_id = contrast_id,
      level = level,
      method = method,
      scope = scope,
      status = "skipped_exists",
      n_edges = if (is.data.frame(edge_df)) nrow(edge_df) else 0L,
      weight_metric = weight_metric,
      source_path = normalize_path_safe(source_path),
      data_tsv = normalize_path_safe(data_tsv),
      data_csv = normalize_path_safe(data_csv),
      figure_pdf = normalize_path_safe(pdf_path),
      figure_png = normalize_path_safe(png_path),
      title = title,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  if (is.null(edge_df) || !is.data.frame(edge_df) || nrow(edge_df) == 0L) {
    return(data.frame(
      contrast_id = contrast_id,
      level = level,
      method = method,
      scope = scope,
      status = "skipped_empty",
      n_edges = 0L,
      weight_metric = weight_metric,
      source_path = normalize_path_safe(source_path),
      data_tsv = NA_character_,
      data_csv = NA_character_,
      figure_pdf = NA_character_,
      figure_png = NA_character_,
      title = title,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  write_tsv_gz(edge_df, data_tsv)
  utils::write.csv(edge_df, data_csv, row.names = FALSE)
  subtitle <- sprintf("Top %d sender→receiver edges by %s", nrow(edge_df), weight_metric)
  paths <- plot_sender_receiver_chord_from_edges(edge_df, output_prefix = prefix, title = title, subtitle = subtitle)
  status <- if (is.null(paths)) "blocked_missing_circlize" else "ok"
  data.frame(
    contrast_id = contrast_id,
    level = level,
    method = method,
    scope = scope,
    status = status,
    n_edges = nrow(edge_df),
    weight_metric = weight_metric,
    source_path = normalize_path_safe(source_path),
    data_tsv = normalize_path_safe(data_tsv),
    data_csv = normalize_path_safe(data_csv),
    figure_pdf = if (is.null(paths)) NA_character_ else paths$pdf,
    figure_png = if (is.null(paths)) NA_character_ else paths$png,
    title = title,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_task_specs <- function(contrast_id, level, cellchat_dir, liana_dir, cellphonedb_dir, consensus_dir, consensus_sender_receiver_df, top_n) {
  specs <- list(
    list(
      method = "cellchat",
      scope = "left",
      title = sprintf("CellChat chord | %s | %s | left", contrast_id, level),
      source_path = file.path(cellchat_dir, contrast_id, level, "cellchat_left_candidates.tsv.gz"),
      weight_metric = "prob",
      edges = load_cellchat_side_edges(file.path(cellchat_dir, contrast_id, level, "cellchat_left_candidates.tsv.gz"), top_n)
    ),
    list(
      method = "cellchat",
      scope = "right",
      title = sprintf("CellChat chord | %s | %s | right", contrast_id, level),
      source_path = file.path(cellchat_dir, contrast_id, level, "cellchat_right_candidates.tsv.gz"),
      weight_metric = "prob",
      edges = load_cellchat_side_edges(file.path(cellchat_dir, contrast_id, level, "cellchat_right_candidates.tsv.gz"), top_n)
    ),
    list(
      method = "cellchat",
      scope = "delta",
      title = sprintf("CellChat delta chord | %s | %s", contrast_id, level),
      source_path = file.path(cellchat_dir, contrast_id, level, "cellchat_delta.tsv.gz"),
      weight_metric = "abs(delta_score_right_minus_left)",
      edges = load_delta_edges(file.path(cellchat_dir, contrast_id, level, "cellchat_delta.tsv.gz"), top_n)
    ),
    list(
      method = "liana",
      scope = "left",
      title = sprintf("LIANA chord | %s | %s | left", contrast_id, level),
      source_path = file.path(liana_dir, contrast_id, level, "liana_left.tsv.gz"),
      weight_metric = "scaled_weight",
      edges = load_liana_side_edges(file.path(liana_dir, contrast_id, level, "liana_left.tsv.gz"), top_n)
    ),
    list(
      method = "liana",
      scope = "right",
      title = sprintf("LIANA chord | %s | %s | right", contrast_id, level),
      source_path = file.path(liana_dir, contrast_id, level, "liana_right.tsv.gz"),
      weight_metric = "scaled_weight",
      edges = load_liana_side_edges(file.path(liana_dir, contrast_id, level, "liana_right.tsv.gz"), top_n)
    ),
    list(
      method = "liana",
      scope = "delta",
      title = sprintf("LIANA delta chord | %s | %s", contrast_id, level),
      source_path = file.path(liana_dir, contrast_id, level, "liana_delta.tsv.gz"),
      weight_metric = "abs(delta_score_right_minus_left)",
      edges = load_delta_edges(file.path(liana_dir, contrast_id, level, "liana_delta.tsv.gz"), top_n)
    ),
    list(
      method = "cellphonedb",
      scope = "left",
      title = sprintf("CellPhoneDB chord | %s | %s | left", contrast_id, level),
      source_path = file.path(cellphonedb_dir, contrast_id, level, "left", "means_long.tsv.gz"),
      weight_metric = "mean",
      edges = load_cellphonedb_side_edges(file.path(cellphonedb_dir, contrast_id, level, "left"), top_n)
    ),
    list(
      method = "cellphonedb",
      scope = "right",
      title = sprintf("CellPhoneDB chord | %s | %s | right", contrast_id, level),
      source_path = file.path(cellphonedb_dir, contrast_id, level, "right", "means_long.tsv.gz"),
      weight_metric = "mean",
      edges = load_cellphonedb_side_edges(file.path(cellphonedb_dir, contrast_id, level, "right"), top_n)
    ),
    list(
      method = "cellphonedb",
      scope = "delta",
      title = sprintf("CellPhoneDB delta chord | %s | %s", contrast_id, level),
      source_path = file.path(cellphonedb_dir, contrast_id, level, "cellphonedb_delta.tsv.gz"),
      weight_metric = "abs(delta_score_right_minus_left)",
      edges = load_delta_edges(file.path(cellphonedb_dir, contrast_id, level, "cellphonedb_delta.tsv.gz"), top_n)
    ),
    list(
      method = "multimethod",
      scope = "consensus",
      title = sprintf("Consensus chord | %s | %s", contrast_id, level),
      source_path = file.path(consensus_dir, "sender_receiver_consensus.tsv.gz"),
      weight_metric = "sum_abs_mean_delta",
      edges = load_consensus_edges(consensus_sender_receiver_df, contrast_id, level, top_n)
    )
  )
  specs
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  run_keys <- collect_run_keys(args$cellchat_dir, args$liana_dir, args$cellphonedb_dir)
  if (!is.data.frame(run_keys) || nrow(run_keys) == 0L) {
    stop("No contrast × level run directories were found across CellChat/LIANA/CellPhoneDB roots.", call. = FALSE)
  }

  consensus_sender_receiver_df <- read_tsv_auto(file.path(args$consensus_dir, "sender_receiver_consensus.tsv.gz"))
  manifest_rows <- list()
  for (idx in seq_len(nrow(run_keys))) {
    contrast_id <- as.character(run_keys$contrast_id[[idx]])
    level <- as.character(run_keys$level[[idx]])
    specs <- build_task_specs(
      contrast_id = contrast_id,
      level = level,
      cellchat_dir = args$cellchat_dir,
      liana_dir = args$liana_dir,
      cellphonedb_dir = args$cellphonedb_dir,
      consensus_dir = args$consensus_dir,
      consensus_sender_receiver_df = consensus_sender_receiver_df,
      top_n = args$top_n
    )
    for (spec in specs) {
      manifest_rows[[length(manifest_rows) + 1L]] <- render_chord_bundle(
        edge_df = spec$edges,
        output_dir = args$output_dir,
        contrast_id = contrast_id,
        level = level,
        method = spec$method,
        scope = spec$scope,
        title = spec$title,
        source_path = spec$source_path,
        weight_metric = spec$weight_metric,
        overwrite = isTRUE(args$overwrite)
      )
    }
  }

  manifest <- do.call(rbind, manifest_rows)
  rownames(manifest) <- NULL
  manifest_path <- file.path(args$output_dir, "per_run_chord_manifest.tsv.gz")
  write_tsv_gz(manifest, manifest_path)
  utils::write.csv(manifest, file.path(args$output_dir, "per_run_chord_manifest.csv"), row.names = FALSE)

  summary <- list(
    status = "ok",
    output_dir = normalize_path_safe(args$output_dir),
    manifest_tsv = normalize_path_safe(manifest_path),
    n_runs = nrow(run_keys),
    n_total_specs = nrow(manifest),
    n_ok = sum(manifest$status == "ok"),
    n_skipped_exists = sum(manifest$status == "skipped_exists"),
    n_skipped_empty = sum(manifest$status == "skipped_empty"),
    n_blocked = sum(grepl("^blocked", manifest$status)),
    methods = sort(unique(as.character(manifest$method))),
    scopes = sort(unique(as.character(manifest$scope))),
    top_n = args$top_n,
    overwrite = isTRUE(args$overwrite)
  )
  cci_write_json(summary, file.path(args$output_dir, "summary.json"))
  message(sprintf(
    "[06c_v1] runs=%d specs=%d ok=%d skipped_exists=%d skipped_empty=%d",
    summary$n_runs,
    summary$n_total_specs,
    summary$n_ok,
    summary$n_skipped_exists,
    summary$n_skipped_empty
  ))
}

main()
