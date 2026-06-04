#!/usr/bin/env Rscript

HELPER_PATH <- "/home/h2048/script/R/cell_communication_helper_20260506_v1_0.R"
source(HELPER_PATH)

DEFAULT_PY_COMM_DIR <- "/home/h2048/data/py/0525/non_unified_airway/communication_lineage_native_export_v2"
DEFAULT_CELLCHAT_DIR <- "/home/h2048/data/R/0525/non_unified_airway/communication/cellchat_lineage_native_export_v2"
DEFAULT_LIANA_DIR <- file.path(DEFAULT_PY_COMM_DIR, "liana_20260525_v2")
DEFAULT_CELLPHONEDB_DIR <- file.path(DEFAULT_PY_COMM_DIR, "cellphonedb_20260525_v2")
DEFAULT_CONSENSUS_DIR <- file.path(DEFAULT_PY_COMM_DIR, "consensus_20260525_v2")
DEFAULT_OUTPUT_DIR <- "/home/h2048/data/R/0525/non_unified_airway/communication/summary_lineage_native_export_v2_20260525_v2"

parse_bool_flag <- function(x, default = TRUE) {
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
    top_n = 25L,
    chord_only = FALSE,
    llm_enabled = TRUE,
    llm_model = "deepseek-reasoner",
    llm_timeout_sec = 180L
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
    } else if (key == "--chord-only" && i < length(args)) {
      out$chord_only <- parse_bool_flag(args[[i + 1L]], default = FALSE)
      i <- i + 2L
    } else if (key == "--llm-enabled" && i < length(args)) {
      out$llm_enabled <- parse_bool_flag(args[[i + 1L]], default = TRUE)
      i <- i + 2L
    } else if (key == "--llm-model" && i < length(args)) {
      out$llm_model <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--llm-timeout-sec" && i < length(args)) {
      out$llm_timeout_sec <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else {
      stop(sprintf("Unknown or incomplete argument: %s", key), call. = FALSE)
    }
  }
  out
}

read_tsv_auto <- function(path) {
  if (!file.exists(path)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, sep = "\t", data.table = FALSE))
  }
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file = gzfile(path), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

list_files_recursive <- function(base_dir, pattern) {
  if (!dir.exists(base_dir)) return(character())
  list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

build_interaction_label <- function(df) {
  paste(
    df$contrast_id,
    df$level,
    paste(df$sender, "->", df$receiver),
    sprintf("[%s-%s]", df$ligand, df$receptor),
    sep = " | "
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_head_df <- function(df, n) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  utils::head(df, max(0L, as.integer(n)))
}

collect_cellchat_pathway_strength <- function(cellchat_dir) {
  files <- list_files_recursive(cellchat_dir, "cellchat_pathway_strength\\.tsv\\.gz$")
  if (length(files) == 0L) return(data.frame())
  dfs <- lapply(files, read_tsv_auto)
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, dfs)
  if (length(dfs) == 0L) return(data.frame())
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out$total_prob <- safe_numeric(out$total_prob)
  out$n_edges <- safe_numeric(out$n_edges)
  out
}

prepare_method_delta_top <- function(all_methods_df, method, top_n = 25L) {
  if (is.null(all_methods_df) || !is.data.frame(all_methods_df) || nrow(all_methods_df) == 0L) return(data.frame())
  df <- all_methods_df[all_methods_df$method %in% method, , drop = FALSE]
  if (nrow(df) == 0L) return(data.frame())
  df$delta_score_right_minus_left <- safe_numeric(df$delta_score_right_minus_left)
  df$abs_delta <- abs(df$delta_score_right_minus_left)
  df$direction_label <- ifelse(
    is.na(df$delta_score_right_minus_left),
    "missing",
    ifelse(df$delta_score_right_minus_left > 0, "right_higher", ifelse(df$delta_score_right_minus_left < 0, "left_higher", "neutral"))
  )
  df$plot_label <- build_interaction_label(df)
  df <- df[order(df$abs_delta, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  safe_head_df(df, top_n)
}

summarise_method_sender_receiver <- function(method_df, top_n = 25L) {
  if (is.null(method_df) || !is.data.frame(method_df) || nrow(method_df) == 0L) return(data.frame())
  method_df$delta_score_right_minus_left <- safe_numeric(method_df$delta_score_right_minus_left)
  method_df$abs_delta <- abs(method_df$delta_score_right_minus_left)
  sr <- aggregate(
    abs_delta ~ contrast_id + level + sender + receiver,
    data = method_df,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  counts <- aggregate(
    interaction_key ~ contrast_id + level + sender + receiver,
    data = method_df,
    FUN = function(x) length(unique(as.character(x)))
  )
  colnames(counts)[colnames(counts) == "interaction_key"] <- "n_interactions"
  sr <- merge(sr, counts, by = c("contrast_id", "level", "sender", "receiver"), all = TRUE)
  sr$plot_label <- paste(sr$contrast_id, sr$level, paste(sr$sender, "->", sr$receiver), sep = " | ")
  sr <- sr[order(sr$abs_delta, sr$n_interactions, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  safe_head_df(sr, top_n)
}

summarise_cellchat_pathway_delta <- function(pathway_df, top_n = 25L) {
  if (is.null(pathway_df) || !is.data.frame(pathway_df) || nrow(pathway_df) == 0L) return(data.frame())
  pathway_df$total_prob <- safe_numeric(pathway_df$total_prob)
  agg <- aggregate(total_prob ~ contrast_id + level + pathway + side, data = pathway_df, FUN = sum, na.rm = TRUE)
  wide <- reshape(agg, idvar = c("contrast_id", "level", "pathway"), timevar = "side", direction = "wide")
  if (!"total_prob.left" %in% colnames(wide)) wide$total_prob.left <- 0
  if (!"total_prob.right" %in% colnames(wide)) wide$total_prob.right <- 0
  wide$left_total_prob <- safe_numeric(wide$total_prob.left)
  wide$right_total_prob <- safe_numeric(wide$total_prob.right)
  wide$delta_total_prob <- wide$right_total_prob - wide$left_total_prob
  wide$abs_delta_total_prob <- abs(wide$delta_total_prob)
  wide$plot_label <- paste(wide$contrast_id, wide$level, wide$pathway, sep = " | ")
  wide <- wide[order(wide$abs_delta_total_prob, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  safe_head_df(wide, top_n)
}

summarise_overlap_pairs <- function(overlap_df) {
  if (is.null(overlap_df) || !is.data.frame(overlap_df) || nrow(overlap_df) == 0L) return(data.frame())
  overlap_df$jaccard <- safe_numeric(overlap_df$jaccard)
  overlap_df$direction_agreement_rate <- safe_numeric(overlap_df$direction_agreement_rate)
  overlap_df <- overlap_df[is.finite(overlap_df$jaccard) | is.finite(overlap_df$direction_agreement_rate), , drop = FALSE]
  if (nrow(overlap_df) == 0L) return(data.frame())
  pair_key <- paste(overlap_df$method_a, overlap_df$method_b, sep = "||")
  split_df <- split(overlap_df, pair_key, drop = TRUE)
  rows <- lapply(split_df, function(df) {
    data.frame(
      method_a = as.character(df$method_a[[1]]),
      method_b = as.character(df$method_b[[1]]),
      jaccard = mean(df$jaccard, na.rm = TRUE),
      direction_agreement_rate = if (all(!is.finite(df$direction_agreement_rate))) NA_real_ else mean(df$direction_agreement_rate, na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

build_overlap_matrix <- function(overlap_pair_df, value_col) {
  methods <- sort(unique(c(as.character(overlap_pair_df$method_a), as.character(overlap_pair_df$method_b))))
  if (length(methods) == 0L) return(matrix(numeric(), nrow = 0L, ncol = 0L))
  mat <- matrix(NA_real_, nrow = length(methods), ncol = length(methods), dimnames = list(methods, methods))
  diag(mat) <- 1
  for (i in seq_len(nrow(overlap_pair_df))) {
    a <- as.character(overlap_pair_df$method_a[[i]])
    b <- as.character(overlap_pair_df$method_b[[i]])
    value <- safe_numeric(overlap_pair_df[[value_col]][[i]])
    mat[a, b] <- value
    mat[b, a] <- value
  }
  mat
}

matrix_to_long <- function(mat, value_name = "value") {
  if (length(mat) == 0L) return(data.frame())
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("row_label", "col_label", value_name)
  df
}

plot_heatmap_matrix <- function(mat, output_prefix, title, zlab) {
  if (is.null(mat) || length(mat) == 0L || nrow(mat) == 0L || ncol(mat) == 0L) return(NULL)
  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)
  plot_df <- matrix_to_long(mat, value_name = zlab)
  value_range <- range(mat, finite = TRUE, na.rm = TRUE)
  if (!all(is.finite(value_range))) value_range <- c(0, 1)
  if (diff(value_range) == 0) value_range <- c(value_range[[1]], value_range[[2]] + 1e-6)
  palette_vals <- grDevices::colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100)

  draw_once <- function(device_fun) {
    device_fun()
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mar = c(6, 6, 4, 5) + 0.1)
    nr <- nrow(mat)
    nc <- ncol(mat)
    plot.new()
    plot.window(xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5), xaxs = "i", yaxs = "i")
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        value <- mat[i, j]
        color_idx <- if (is.na(value)) {
          NA_integer_
        } else {
          floor((value - value_range[[1]]) / diff(value_range) * 99) + 1L
        }
        rect(
          xleft = j - 0.5,
          ybottom = nr - i + 0.5,
          xright = j + 0.5,
          ytop = nr - i + 1.5,
          col = if (is.na(color_idx)) "grey95" else palette_vals[pmax(1L, pmin(100L, color_idx))],
          border = "white"
        )
        text(j, nr - i + 1, labels = ifelse(is.na(value), "NA", format(round(value, 2), nsmall = 2)), cex = 0.9)
      }
    }
    axis(1, at = seq_len(nc), labels = colnames(mat), las = 2)
    axis(2, at = seq_len(nr), labels = rev(rownames(mat)), las = 2)
    box()
    title(main = title)
    par(xpd = TRUE)
    legend_x <- nc + 0.95
    legend_y <- seq(0.75, nr + 0.25, length.out = 100)
    rect(legend_x, head(legend_y, -1), legend_x + 0.25, tail(legend_y, -1), col = palette_vals, border = NA)
    text(legend_x + 0.45, 0.75, labels = sprintf("%.2f", value_range[[1]]), adj = c(0, 0.5), cex = 0.8)
    text(legend_x + 0.45, nr + 0.25, labels = sprintf("%.2f", value_range[[2]]), adj = c(0, 0.5), cex = 0.8)
    text(legend_x + 0.12, nr + 0.8, labels = zlab, srt = 90, cex = 0.8)
    grDevices::dev.off()
  }

  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  draw_once(function() grDevices::pdf(pdf_path, width = 8, height = 7))
  draw_once(function() grDevices::png(png_path, width = 1800, height = 1500, res = 180))
  list(pdf = cci_path_normalize(pdf_path), png = cci_path_normalize(png_path), data = plot_df)
}

read_csv_auto <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_sender_receiver_chord_edges <- function(sr_df,
                                                weight_col,
                                                method,
                                                source_table,
                                                top_n = 25L) {
  if (is.null(sr_df) || !is.data.frame(sr_df) || nrow(sr_df) == 0L) return(data.frame())
  required_cols <- c("sender", "receiver", weight_col)
  if (!all(required_cols %in% colnames(sr_df))) return(data.frame())
  df <- as.data.frame(sr_df, stringsAsFactors = FALSE, check.names = FALSE)
  df$sender <- trimws(as.character(df$sender))
  df$receiver <- trimws(as.character(df$receiver))
  df$edge_weight <- safe_numeric(df[[weight_col]])
  keep <- nzchar(df$sender) & nzchar(df$receiver) & is.finite(df$edge_weight) & df$edge_weight > 0
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0L) return(data.frame())
  df <- df[order(df$edge_weight, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  df <- safe_head_df(df, top_n)
  edge_weight <- aggregate(edge_weight ~ sender + receiver, data = df, FUN = sum, na.rm = TRUE)
  n_source_rows <- aggregate(edge_weight ~ sender + receiver, data = df, FUN = length)
  colnames(n_source_rows)[colnames(n_source_rows) == "edge_weight"] <- "n_source_rows"
  out <- merge(edge_weight, n_source_rows, by = c("sender", "receiver"), all = TRUE)
  if ("n_interactions" %in% colnames(df)) {
    interaction_counts <- aggregate(safe_numeric(df$n_interactions) ~ sender + receiver, data = df, FUN = sum, na.rm = TRUE)
    colnames(interaction_counts)[colnames(interaction_counts) == "safe_numeric(df$n_interactions)"] <- "n_interactions"
    out <- merge(out, interaction_counts, by = c("sender", "receiver"), all = TRUE)
  }
  if ("n_consensus_interactions" %in% colnames(df)) {
    consensus_counts <- aggregate(safe_numeric(df$n_consensus_interactions) ~ sender + receiver, data = df, FUN = sum, na.rm = TRUE)
    colnames(consensus_counts)[colnames(consensus_counts) == "safe_numeric(df$n_consensus_interactions)"] <- "n_consensus_interactions"
    out <- merge(out, consensus_counts, by = c("sender", "receiver"), all = TRUE)
  }
  out <- out[order(out$edge_weight, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  out$method <- method
  out$source_table <- source_table
  out$weight_col <- weight_col
  out$edge_rank <- seq_len(nrow(out))
  out$edge_label <- paste(out$sender, "->", out$receiver)
  max_weight <- max(out$edge_weight, na.rm = TRUE)
  out$edge_weight_scaled <- if (is.finite(max_weight) && max_weight > 0) out$edge_weight / max_weight else NA_real_
  rownames(out) <- NULL
  out
}

plot_sender_receiver_chord <- function(sr_df,
                                       weight_col,
                                       output_prefix,
                                       title,
                                       method,
                                       source_table,
                                       top_n = 25L) {
  if (!requireNamespace("circlize", quietly = TRUE)) return(NULL)
  edge_df <- prepare_sender_receiver_chord_edges(
    sr_df = sr_df,
    weight_col = weight_col,
    method = method,
    source_table = source_table,
    top_n = top_n
  )
  if (nrow(edge_df) == 0L) return(NULL)
  chord_df <- edge_df[, c("sender", "receiver", "edge_weight"), drop = FALSE]
  sectors <- sort(unique(c(chord_df$sender, chord_df$receiver)))
  grid_col <- grDevices::hcl(
    h = seq(15, 375, length.out = length(sectors) + 1L)[seq_along(sectors)],
    c = 70,
    l = 65
  )
  names(grid_col) <- sectors
  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)
  draw_once <- function(device_fun) {
    device_opened <- FALSE
    tryCatch({
      device_fun()
      device_opened <- TRUE
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      circlize::circos.clear()
      circlize::circos.par(
        start.degree = 90,
        gap.after = rep(if (length(sectors) > 18L) 1.2 else 2, length(sectors)),
        track.margin = c(0.01, 0.01),
        canvas.xlim = c(-1.35, 1.35),
        canvas.ylim = c(-1.35, 1.35),
        points.overflow.warning = FALSE
      )
      circlize::chordDiagram(
        x = chord_df,
        grid.col = grid_col,
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        diffHeight = -0.04,
        link.arr.type = "big.arrow",
        link.arr.length = 0.18,
        link.sort = TRUE,
        link.largest.ontop = TRUE,
        transparency = 0.35,
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
            cex = if (length(sectors) > 18L) 0.38 else if (length(sectors) > 12L) 0.46 else 0.52
          )
        }
      )
      title(main = title, cex.main = 1.05)
      mtext(
        sprintf("Directed chord diagram; edge width is proportional to %s; top %d rows aggregated by sender-receiver pair.", weight_col, top_n),
        side = 1,
        line = -1,
        cex = 0.72
      )
      grDevices::dev.off()
      device_opened <- FALSE
    }, finally = {
      try(circlize::circos.clear(), silent = TRUE)
      if (isTRUE(device_opened)) try(grDevices::dev.off(), silent = TRUE)
    })
  }
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  draw_once(function() grDevices::pdf(pdf_path, width = 12, height = 12, onefile = FALSE))
  draw_once(function() grDevices::png(png_path, width = 3000, height = 3000, res = 240))
  list(pdf = cci_path_normalize(pdf_path), png = cci_path_normalize(png_path), data = edge_df)
}

read_figure_manifest <- function(manifest_path) {
  if (!file.exists(manifest_path)) return(data.frame())
  utils::read.delim(
    manifest_path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_figure_manifest <- function(manifest, manifest_path) {
  if (is.null(manifest) || !is.data.frame(manifest) || nrow(manifest) == 0L) {
    unlink(manifest_path)
  } else {
    utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  }
  invisible(manifest_path)
}

deduplicate_figure_manifest <- function(manifest_path) {
  if (!file.exists(manifest_path)) return(invisible(manifest_path))
  manifest <- read_figure_manifest(manifest_path)
  if (!is.data.frame(manifest) || nrow(manifest) == 0L || !all(c("figure_type", "figure_path") %in% colnames(manifest))) return(invisible(manifest_path))
  key <- paste(as.character(manifest$figure_type), as.character(manifest$figure_path), sep = "||")
  manifest <- manifest[!duplicated(key, fromLast = TRUE), , drop = FALSE]
  write_figure_manifest(manifest, manifest_path)
}

remove_existing_chord_manifest_rows <- function(manifest_path) {
  if (!file.exists(manifest_path)) return(invisible(manifest_path))
  manifest <- read_figure_manifest(manifest_path)
  if (!is.data.frame(manifest) || nrow(manifest) == 0L || !"figure_type" %in% colnames(manifest)) return(invisible(manifest_path))
  manifest <- manifest[!grepl("chord", as.character(manifest$figure_type), ignore.case = TRUE), , drop = FALSE]
  write_figure_manifest(manifest, manifest_path)
  invisible(manifest_path)
}

count_figure_manifest_rows <- function(manifest_path) {
  if (!file.exists(manifest_path)) return(0L)
  manifest <- read_figure_manifest(manifest_path)
  if (!is.data.frame(manifest) || nrow(manifest) == 0L) return(0L)
  nrow(manifest)
}

update_summary_registered_figures <- function(output_dir, figure_manifest_path, chord_only = FALSE) {
  if (file.exists(figure_manifest_path)) {
    deduplicate_figure_manifest(figure_manifest_path)
  }
  n_registered_figures <- count_figure_manifest_rows(figure_manifest_path)
  summary_path <- file.path(output_dir, "summary.json")
  summary <- list(status = "ok")
  if (file.exists(summary_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    summary <- tryCatch(jsonlite::fromJSON(summary_path, simplifyVector = FALSE), error = function(e) summary)
  }
  summary$n_registered_figures <- n_registered_figures
  summary$output_dir <- cci_path_normalize(output_dir)
  summary$figure_manifest <- cci_path_normalize(figure_manifest_path)
  summary$chord_diagrams_added <- TRUE
  summary$chord_only_last_update <- if (isTRUE(chord_only)) format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z") else summary$chord_only_last_update
  cci_write_json(summary, summary_path)
  n_registered_figures
}

register_sender_receiver_chord <- function(sr_df,
                                           weight_col,
                                           output_prefix,
                                           title,
                                           method,
                                           figure_type,
                                           source_table,
                                           top_n,
                                           llm_config,
                                           manifest_path) {
  chord_paths <- plot_sender_receiver_chord(
    sr_df = sr_df,
    weight_col = weight_col,
    output_prefix = output_prefix,
    title = title,
    method = method,
    source_table = source_table,
    top_n = top_n
  )
  register_plot_result(
    chord_paths,
    method = method,
    figure_type = figure_type,
    title = title,
    extra_context = c(
      "Chord links are directed from sender to receiver.",
      sprintf("Edge width uses `%s`; top rows are aggregated by unique sender-receiver pair.", weight_col),
      "This is a visualization of computational CCI candidate burden, not protein-level proof of direct signaling."
    ),
    llm_config = llm_config,
    manifest_path = manifest_path
  )
}

run_chord_only_outputs <- function(tables_dir, figures_dir, figure_manifest_path, top_n, llm_config, output_dir) {
  remove_existing_chord_manifest_rows(figure_manifest_path)
  statuses <- list()
  for (method in c("cellchat", "liana", "cellphonedb")) {
    sr_path <- file.path(tables_dir, sprintf("%s_sender_receiver_top.csv", method))
    sr_top <- read_csv_auto(sr_path)
    if (nrow(sr_top) > 0L) {
      statuses[[length(statuses) + 1L]] <- register_sender_receiver_chord(
        sr_df = sr_top,
        weight_col = "abs_delta",
        output_prefix = file.path(figures_dir, sprintf("%s_sender_receiver_top_chord", method)),
        title = sprintf("Chord diagram: top %s sender-receiver burdens", method),
        method = method,
        figure_type = "sender_receiver_burden_chord",
        source_table = basename(sr_path),
        top_n = top_n,
        llm_config = llm_config,
        manifest_path = figure_manifest_path
      )
    }
  }
  consensus_sr_path <- file.path(tables_dir, "sender_receiver_consensus_top.csv")
  consensus_sr_top <- read_csv_auto(consensus_sr_path)
  if (nrow(consensus_sr_top) > 0L) {
    statuses[[length(statuses) + 1L]] <- register_sender_receiver_chord(
      sr_df = consensus_sr_top,
      weight_col = "sum_abs_mean_delta",
      output_prefix = file.path(figures_dir, "sender_receiver_consensus_top_chord"),
      title = "Chord diagram: top consensus sender-receiver programs",
      method = "multimethod",
      figure_type = "sender_receiver_consensus_chord",
      source_table = basename(consensus_sr_path),
      top_n = top_n,
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )
  }
  n_registered_figures <- update_summary_registered_figures(output_dir, figure_manifest_path, chord_only = TRUE)
  message(sprintf("[06b_v2] chord_only=true chord_figures=%d manifest_figures=%d", length(Filter(Negate(is.null), statuses)), n_registered_figures))
  invisible(statuses)
}

export_table_with_llm <- function(df, tsv_path, csv_path, title, scenario, extra_context = NULL, llm_config = list()) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(NULL)
  write_tsv_gz(df, tsv_path)
  cci_write_csv(df, csv_path)
  prefix <- sub("\\.csv$", "", basename(csv_path))
  cci_run_table_llm_interpretation(
    df = df,
    output_dir = dirname(csv_path),
    prefix = prefix,
    title = title,
    scenario = scenario,
    extra_context = extra_context,
    llm_config = llm_config
  )
}

register_plot_result <- function(paths, method, figure_type, title, extra_context = NULL, llm_config = list(), manifest_path) {
  if (is.null(paths) || is.null(paths$png) || is.null(paths$data)) return(NULL)
  cci_register_figure(
    figure_path = paths$png,
    data = paths$data,
    method = method,
    figure_type = figure_type,
    title = title,
    extra_context = extra_context,
    llm_config = llm_config,
    manifest_path = manifest_path,
    scenario = "A_scrna | non_unified_airway multimethod"
  )
}

method_fill_color <- function(method) {
  switch(
    tolower(method),
    "cellchat" = "#4E79A7",
    "liana" = "#F28E2B",
    "cellphonedb" = "#E15759",
    "#76B7B2"
  )
}

prepare_llm_table_spec <- function(df, method_label, evidence_scope, score_col, confidence_default = "putative") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(NULL)
  out <- df
  out$score <- safe_numeric(out[[score_col]])
  out$interaction <- if ("interaction" %in% colnames(out)) out$interaction else build_interaction_label(out)
  out$confidence_label <- if (!"confidence_label" %in% colnames(out)) confidence_default else out$confidence_label
  list(
    data = out,
    method = method_label,
    evidence_scope = evidence_scope,
    statistical_unit = "cells for screening only",
    confidence_default = confidence_default
  )
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  tables_dir <- file.path(args$output_dir, "tables")
  figures_dir <- file.path(args$output_dir, "figures")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  figure_manifest_path <- file.path(figures_dir, "figure_companion_manifest.tsv")
  if (!isTRUE(args$chord_only) && file.exists(figure_manifest_path)) unlink(figure_manifest_path)

  llm_config <- list(
    enabled = isTRUE(args$llm_enabled),
    model = args$llm_model,
    timeout_sec = args$llm_timeout_sec
  )
  if (isTRUE(args$llm_enabled) && requireNamespace("fanyi", quietly = TRUE)) {
    key_info <- cci_resolve_deepseek_key(write_env_placeholder = FALSE)
    if (isTRUE(key_info$has_live_key)) {
      try(
        fanyi::set_translate_option(key = key_info$api_key, source = "deepseek"),
        silent = TRUE
      )
    }
  }

  if (isTRUE(args$chord_only)) {
    run_chord_only_outputs(
      tables_dir = tables_dir,
      figures_dir = figures_dir,
      figure_manifest_path = figure_manifest_path,
      top_n = args$top_n,
      llm_config = llm_config,
      output_dir = args$output_dir
    )
    return(invisible(NULL))
  }

  all_methods_df <- read_tsv_auto(file.path(args$consensus_dir, "all_method_delta_rows.tsv.gz"))
  method_summary_df <- read_tsv_auto(file.path(args$consensus_dir, "method_summary.tsv.gz"))
  consensus_df <- read_tsv_auto(file.path(args$consensus_dir, "consensus_interactions.tsv.gz"))
  sender_receiver_df <- read_tsv_auto(file.path(args$consensus_dir, "sender_receiver_consensus.tsv.gz"))
  overlap_df <- read_tsv_auto(file.path(args$consensus_dir, "method_overlap_long.tsv.gz"))
  support_distribution_df <- read_tsv_auto(file.path(args$consensus_dir, "support_distribution.tsv.gz"))
  cellchat_pathway_df <- collect_cellchat_pathway_strength(args$cellchat_dir)

  figure_statuses <- list()
  table_statuses <- list()

  if (nrow(method_summary_df) > 0L) {
    table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
      method_summary_df,
      file.path(tables_dir, "method_summary.tsv.gz"),
      file.path(tables_dir, "method_summary.csv"),
      title = "Method-level communication delta summary",
      scenario = "A_scrna | non_unified_airway multimethod | method summary",
      extra_context = c("Each row summarizes one method within one contrast and one annotation level.", "This table is primarily for execution and coverage auditing."),
      llm_config = llm_config
    )
  }

  llm_specs <- list()
  for (method in c("cellchat", "liana", "cellphonedb")) {
    method_all <- if (nrow(all_methods_df) > 0L) all_methods_df[all_methods_df$method %in% method, , drop = FALSE] else data.frame()
    method_top <- prepare_method_delta_top(all_methods_df, method, top_n = args$top_n)
    if (nrow(method_top) > 0L) {
      table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
        method_top,
        file.path(tables_dir, sprintf("%s_delta_top.tsv.gz", method)),
        file.path(tables_dir, sprintf("%s_delta_top.csv", method)),
        title = sprintf("%s top delta interactions", toupper(substring(method, 1L, 1L)) %+% substring(method, 2L)),
        scenario = "A_scrna | non_unified_airway multimethod | method-specific delta top table",
        extra_context = c(
          sprintf("Method = %s", method),
          "The ranking metric is absolute right-minus-left delta score.",
          "Interpret these as computationally inferred candidate differences, not direct validated signaling."
        ),
        llm_config = llm_config
      )
      plot_paths <- cci_plot_bar(
        method_top,
        label_col = "plot_label",
        value_col = "abs_delta",
        output_prefix = file.path(figures_dir, sprintf("%s_delta_top", method)),
        title = sprintf("Top %s delta interactions", method),
        xlab = "Absolute right-minus-left delta score",
        top_n = args$top_n,
        fill = method_fill_color(method)
      )
      figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
        plot_paths,
        method = method,
        figure_type = "delta_top_barplot",
        title = sprintf("Top %s delta interactions", method),
        extra_context = c("Bar height encodes absolute delta magnitude.", "Rows come from method-specific left vs right contrast runs."),
        llm_config = llm_config,
        manifest_path = figure_manifest_path
      )

      sr_top <- summarise_method_sender_receiver(method_all, top_n = args$top_n)
      if (nrow(sr_top) > 0L) {
        table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
          sr_top,
          file.path(tables_dir, sprintf("%s_sender_receiver_top.tsv.gz", method)),
          file.path(tables_dir, sprintf("%s_sender_receiver_top.csv", method)),
          title = sprintf("%s sender-receiver burden summary", method),
          scenario = "A_scrna | non_unified_airway multimethod | sender receiver burden",
          extra_context = c("This is an aggregate of interaction-level delta magnitudes per sender-receiver pair."),
          llm_config = llm_config
        )
        sr_paths <- cci_plot_bar(
          sr_top,
          label_col = "plot_label",
          value_col = "abs_delta",
          output_prefix = file.path(figures_dir, sprintf("%s_sender_receiver_top", method)),
          title = sprintf("Top %s sender-receiver burdens", method),
          xlab = "Summed absolute delta score",
          top_n = args$top_n,
          fill = method_fill_color(method)
        )
        figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
          sr_paths,
          method = method,
          figure_type = "sender_receiver_burden_barplot",
          title = sprintf("Top %s sender-receiver burdens", method),
          extra_context = c("This summarizes multiple ligand-receptor records into sender-receiver burden."),
          llm_config = llm_config,
          manifest_path = figure_manifest_path
        )
        figure_statuses[[length(figure_statuses) + 1L]] <- register_sender_receiver_chord(
          sr_df = sr_top,
          weight_col = "abs_delta",
          output_prefix = file.path(figures_dir, sprintf("%s_sender_receiver_top_chord", method)),
          title = sprintf("Chord diagram: top %s sender-receiver burdens", method),
          method = method,
          figure_type = "sender_receiver_burden_chord",
          source_table = sprintf("%s_sender_receiver_top.csv", method),
          top_n = args$top_n,
          llm_config = llm_config,
          manifest_path = figure_manifest_path
        )
      }

      llm_specs[[method]] <- prepare_llm_table_spec(
        method_top,
        method_label = sprintf("%s delta", method),
        evidence_scope = "site-contrast delta screening",
        score_col = "abs_delta"
      )
    }
  }

  pathway_top <- summarise_cellchat_pathway_delta(cellchat_pathway_df, top_n = args$top_n)
  if (nrow(pathway_top) > 0L) {
    table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
      pathway_top,
      file.path(tables_dir, "cellchat_pathway_delta_top.tsv.gz"),
      file.path(tables_dir, "cellchat_pathway_delta_top.csv"),
      title = "CellChat pathway-level delta summary",
      scenario = "A_scrna | non_unified_airway multimethod | CellChat pathway delta",
      extra_context = c("Pathway delta is computed from total pathway probability: right minus left."),
      llm_config = llm_config
    )
    pathway_paths <- cci_plot_bar(
      pathway_top,
      label_col = "plot_label",
      value_col = "abs_delta_total_prob",
      output_prefix = file.path(figures_dir, "cellchat_pathway_delta_top"),
      title = "Top CellChat pathway deltas",
      xlab = "Absolute right-minus-left total pathway probability",
      top_n = args$top_n,
      fill = "#76B7B2"
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
      pathway_paths,
      method = "cellchat",
      figure_type = "pathway_delta_barplot",
      title = "Top CellChat pathway deltas",
      extra_context = c("Pathway-level totals summarize CellChat pathway probability burden by contrast side."),
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )
  }

  if (nrow(consensus_df) > 0L) {
    consensus_df$n_methods <- safe_numeric(consensus_df$n_methods)
    consensus_df$mean_delta_score <- safe_numeric(consensus_df$mean_delta_score)
    consensus_df$max_abs_delta <- safe_numeric(consensus_df$max_abs_delta)
    consensus_df$plot_label <- build_interaction_label(consensus_df)
    consensus_top <- consensus_df[order(consensus_df$n_methods, consensus_df$max_abs_delta, consensus_df$mean_delta_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    consensus_top <- safe_head_df(consensus_top, args$top_n)
    table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
      consensus_top,
      file.path(tables_dir, "consensus_interactions_top.tsv.gz"),
      file.path(tables_dir, "consensus_interactions_top.csv"),
      title = "Top multi-method consensus interactions",
      scenario = "A_scrna | non_unified_airway multimethod | interaction consensus",
      extra_context = c("Primary ranking: number of supporting methods, then effect size.", "Direction consistency is reported explicitly and should be used as a guardrail."),
      llm_config = llm_config
    )
    consensus_paths <- cci_plot_bar(
      consensus_top,
      label_col = "plot_label",
      value_col = "n_methods",
      output_prefix = file.path(figures_dir, "consensus_support_top"),
      title = "Top multi-method consensus interactions",
      xlab = "Number of supporting methods",
      top_n = args$top_n,
      fill = "#9C755F"
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
      consensus_paths,
      method = "multimethod",
      figure_type = "consensus_support_barplot",
      title = "Top multi-method consensus interactions",
      extra_context = c("Higher bar means more methods support the same interaction key.", "Use direction consistency before claiming robust left/right shift."),
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )

    llm_specs[["consensus"]] <- prepare_llm_table_spec(
      consensus_top,
      method_label = "multi-method consensus",
      evidence_scope = "cross-method consensus support",
      score_col = "n_methods",
      confidence_default = "candidate"
    )
  }

  if (nrow(sender_receiver_df) > 0L) {
    sender_receiver_df$max_support <- safe_numeric(sender_receiver_df$max_support)
    sender_receiver_df$sum_abs_mean_delta <- safe_numeric(sender_receiver_df$sum_abs_mean_delta)
    sender_receiver_df$sender_receiver_label <- if ("sender_receiver_label" %in% colnames(sender_receiver_df)) {
      sender_receiver_df$sender_receiver_label
    } else {
      paste(sender_receiver_df$contrast_id, sender_receiver_df$level, paste(sender_receiver_df$sender, "->", sender_receiver_df$receiver), sep = " | ")
    }
    sender_receiver_top <- sender_receiver_df[order(sender_receiver_df$max_support, sender_receiver_df$sum_abs_mean_delta, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    sender_receiver_top <- safe_head_df(sender_receiver_top, args$top_n)
    table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
      sender_receiver_top,
      file.path(tables_dir, "sender_receiver_consensus_top.tsv.gz"),
      file.path(tables_dir, "sender_receiver_consensus_top.csv"),
      title = "Top consensus sender-receiver programs",
      scenario = "A_scrna | non_unified_airway multimethod | sender receiver consensus",
      extra_context = c("This table summarizes interaction-level consensus into sender-receiver programs."),
      llm_config = llm_config
    )
    sender_receiver_paths <- cci_plot_bar(
      sender_receiver_top,
      label_col = "sender_receiver_label",
      value_col = "sum_abs_mean_delta",
      output_prefix = file.path(figures_dir, "sender_receiver_consensus_top"),
      title = "Top consensus sender-receiver programs",
      xlab = "Summed consensus absolute delta",
      top_n = args$top_n,
      fill = "#59A14F"
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
      sender_receiver_paths,
      method = "multimethod",
      figure_type = "sender_receiver_consensus_barplot",
      title = "Top consensus sender-receiver programs",
      extra_context = c("Sender-receiver burden aggregates multiple consensus-supported interactions."),
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_sender_receiver_chord(
      sr_df = sender_receiver_top,
      weight_col = "sum_abs_mean_delta",
      output_prefix = file.path(figures_dir, "sender_receiver_consensus_top_chord"),
      title = "Chord diagram: top consensus sender-receiver programs",
      method = "multimethod",
      figure_type = "sender_receiver_consensus_chord",
      source_table = "sender_receiver_consensus_top.csv",
      top_n = args$top_n,
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )
  }

  if (nrow(support_distribution_df) > 0L) {
    support_distribution_df$n_interactions <- safe_numeric(support_distribution_df$n_interactions)
    support_plot_df <- aggregate(n_interactions ~ support_tier + direction_consistency, data = support_distribution_df, FUN = sum, na.rm = TRUE)
    support_plot_df$plot_label <- paste(support_plot_df$support_tier, support_plot_df$direction_consistency, sep = " | ")
    support_plot_df <- support_plot_df[order(support_plot_df$n_interactions, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
      support_plot_df,
      file.path(tables_dir, "support_distribution_summary.tsv.gz"),
      file.path(tables_dir, "support_distribution_summary.csv"),
      title = "Support-tier and direction-consistency distribution",
      scenario = "A_scrna | non_unified_airway multimethod | support distribution",
      extra_context = c("This summarizes how many interactions fall into single/two/three-method support tiers and whether directions agree across methods."),
      llm_config = llm_config
    )
    support_paths <- cci_plot_bar(
      support_plot_df,
      label_col = "plot_label",
      value_col = "n_interactions",
      output_prefix = file.path(figures_dir, "support_distribution_summary"),
      title = "Support-tier and direction-consistency distribution",
      xlab = "Number of interactions",
      top_n = args$top_n,
      fill = "#EDC948"
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
      support_paths,
      method = "multimethod",
      figure_type = "support_distribution_barplot",
      title = "Support-tier and direction-consistency distribution",
      extra_context = c("A healthy multi-method result should not be dominated by discordant directions when support tier is high."),
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )
  }

  overlap_pair_df <- summarise_overlap_pairs(overlap_df)
  if (nrow(overlap_pair_df) > 0L) {
    table_statuses[[length(table_statuses) + 1L]] <- export_table_with_llm(
      overlap_pair_df,
      file.path(tables_dir, "method_overlap_pair_summary.tsv.gz"),
      file.path(tables_dir, "method_overlap_pair_summary.csv"),
      title = "Pairwise method overlap summary",
      scenario = "A_scrna | non_unified_airway multimethod | method overlap",
      extra_context = c("Mean Jaccard is averaged across contrast-level combinations.", "Direction agreement summarizes whether overlapping interactions point in the same left/right direction."),
      llm_config = llm_config
    )

    jaccard_mat <- build_overlap_matrix(overlap_pair_df, "jaccard")
    jaccard_paths <- plot_heatmap_matrix(
      jaccard_mat,
      output_prefix = file.path(figures_dir, "method_overlap_jaccard_heatmap"),
      title = "Method-overlap mean Jaccard heatmap",
      zlab = "mean_jaccard"
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
      jaccard_paths,
      method = "multimethod",
      figure_type = "method_overlap_jaccard_heatmap",
      title = "Method-overlap mean Jaccard heatmap",
      extra_context = c("Diagonal equals 1 by construction.", "Off-diagonal values reflect average overlap between methods across contrasts and levels."),
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )

    agreement_mat <- build_overlap_matrix(overlap_pair_df, "direction_agreement_rate")
    agreement_paths <- plot_heatmap_matrix(
      agreement_mat,
      output_prefix = file.path(figures_dir, "method_overlap_direction_agreement_heatmap"),
      title = "Method-overlap direction-agreement heatmap",
      zlab = "direction_agreement"
    )
    figure_statuses[[length(figure_statuses) + 1L]] <- register_plot_result(
      agreement_paths,
      method = "multimethod",
      figure_type = "method_overlap_direction_agreement_heatmap",
      title = "Method-overlap direction-agreement heatmap",
      extra_context = c("This heatmap asks: when two methods overlap on the same interaction, do they agree on which side is higher?"),
      llm_config = llm_config,
      manifest_path = figure_manifest_path
    )
  }

  llm_specs <- Filter(Negate(is.null), llm_specs)
  llm_interaction_df <- if (length(llm_specs) > 0L) {
    cci_build_llm_interaction_table(
      llm_specs,
      scenario = "non_unified_airway_lineage_native_export_v2_multimethod",
      default_statistical_unit = "cells for screening only",
      top_n_per_method = args$top_n
    )
  } else {
    data.frame()
  }

  if (nrow(llm_interaction_df) > 0L) {
    write_tsv_gz(llm_interaction_df, file.path(tables_dir, "llm_interaction_review_table.tsv.gz"))
    cci_write_csv(llm_interaction_df, file.path(tables_dir, "llm_interaction_review_table.csv"))
    table_statuses[[length(table_statuses) + 1L]] <- cci_run_table_llm_interpretation(
      df = llm_interaction_df,
      output_dir = tables_dir,
      prefix = "llm_interaction_review_table",
      title = "Integrated multi-method communication review table",
      scenario = "A_scrna | non_unified_airway multimethod | integrated LLM review",
      extra_context = c(
        "This table harmonizes method-specific top interactions and multi-method consensus rows.",
        "Use conservative language: candidate / putative / computationally inferred."
      ),
      llm_config = llm_config
    )
  }

  n_registered_figures <- if (file.exists(figure_manifest_path)) nrow(read_tsv_auto(figure_manifest_path)) else 0L
  summary <- list(
    status = "ok",
    methods_available = sort(unique(as.character(all_methods_df$method))),
    n_all_method_rows = if (nrow(all_methods_df) > 0L) nrow(all_methods_df) else 0L,
    n_consensus_rows = if (nrow(consensus_df) > 0L) nrow(consensus_df) else 0L,
    n_sender_receiver_rows = if (nrow(sender_receiver_df) > 0L) nrow(sender_receiver_df) else 0L,
    n_overlap_rows = if (nrow(overlap_df) > 0L) nrow(overlap_df) else 0L,
    n_cellchat_pathway_rows = if (nrow(cellchat_pathway_df) > 0L) nrow(cellchat_pathway_df) else 0L,
    n_registered_figures = n_registered_figures,
    llm_enabled = isTRUE(args$llm_enabled),
    output_dir = cci_path_normalize(args$output_dir),
    figure_manifest = cci_path_normalize(figure_manifest_path)
  )
  cci_write_json(summary, file.path(args$output_dir, "summary.json"))
  message(sprintf("[06b_v2] all_method_rows=%d consensus_rows=%d figures=%d", summary$n_all_method_rows, summary$n_consensus_rows, summary$n_registered_figures))
}

`%+%` <- function(a, b) paste0(a, b)

main()
