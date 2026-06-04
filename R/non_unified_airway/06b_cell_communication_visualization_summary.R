#!/usr/bin/env Rscript

HELPER_PATH <- "/home/h2048/script/R/cell_communication_helper_20260506_v1_0.R"
source(HELPER_PATH)

DEFAULT_CELLCHAT_DIR <- "/home/h2048/data/R/0525/non_unified_airway/communication/cellchat"
DEFAULT_CONSENSUS_DIR <- "/home/h2048/data/py/0525/non_unified_airway/communication/consensus"
DEFAULT_OUTPUT_DIR <- "/home/h2048/data/R/0525/non_unified_airway/communication/summary"

parse_args_simple <- function(args) {
  out <- list(
    cellchat_dir = DEFAULT_CELLCHAT_DIR,
    consensus_dir = DEFAULT_CONSENSUS_DIR,
    output_dir = DEFAULT_OUTPUT_DIR,
    top_n = 25L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--cellchat-dir" && i < length(args)) {
      out$cellchat_dir <- args[[i + 1L]]
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
    } else {
      stop(sprintf("Unknown or incomplete argument: %s", key), call. = FALSE)
    }
  }
  out
}

read_tsv_auto <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing TSV input: %s", path), call. = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, sep = "\t", data.table = FALSE))
  }
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_local <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file = gzfile(path), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_json_local <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    writeLines(utils::capture.output(str(x)), con = path)
  }
}

list_files_recursive <- function(base_dir, pattern) {
  if (!dir.exists(base_dir)) return(character())
  list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

collect_cellchat_delta <- function(cellchat_dir) {
  files <- list_files_recursive(cellchat_dir, "cellchat_delta\\.tsv\\.gz$")
  if (length(files) == 0L) return(data.frame())
  dfs <- lapply(files, read_tsv_auto)
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

plot_top_delta <- function(df, output_prefix, title, top_n = 25L, value_col = "abs_delta") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(NULL)
  cci_plot_bar(
    df,
    label_col = "plot_label",
    value_col = value_col,
    output_prefix = output_prefix,
    title = title,
    xlab = value_col,
    top_n = top_n,
    fill = "#4E79A7"
  )
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  cellchat_df <- collect_cellchat_delta(args$cellchat_dir)
  consensus_path <- file.path(args$consensus_dir, "consensus_interactions.tsv.gz")
  consensus_df <- if (file.exists(consensus_path)) read_tsv_auto(consensus_path) else data.frame()

  if (nrow(cellchat_df) > 0L) {
    cellchat_df$left_score <- suppressWarnings(as.numeric(cellchat_df$left_score))
    cellchat_df$right_score <- suppressWarnings(as.numeric(cellchat_df$right_score))
    cellchat_df$delta_score_right_minus_left <- suppressWarnings(as.numeric(cellchat_df$delta_score_right_minus_left))
    cellchat_df$abs_delta <- abs(cellchat_df$delta_score_right_minus_left)
    cellchat_df$plot_label <- paste(cellchat_df$contrast_id, cellchat_df$level, cellchat_df$sender, "->", cellchat_df$receiver, paste0("[", cellchat_df$ligand, "-", cellchat_df$receptor, "]"))
    top_all <- cellchat_df[order(cellchat_df$abs_delta, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    top_all <- head(top_all, args$top_n)
    write_tsv_local(cellchat_df, file.path(args$output_dir, "cellchat_delta_all.tsv.gz"))
    write_tsv_local(top_all, file.path(args$output_dir, "cellchat_delta_top.tsv.gz"))
    plot_top_delta(top_all, file.path(args$output_dir, "cellchat_delta_top"), "Top CellChat delta interactions", top_n = args$top_n)

    sr <- aggregate(
      abs_delta ~ contrast_id + level + sender + receiver,
      data = cellchat_df,
      FUN = function(x) sum(x, na.rm = TRUE)
    )
    sr$plot_label <- paste(sr$contrast_id, sr$level, sr$sender, "->", sr$receiver)
    sr <- sr[order(sr$abs_delta, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    sr <- head(sr, args$top_n)
    write_tsv_local(sr, file.path(args$output_dir, "cellchat_sender_receiver_top.tsv.gz"))
    cci_plot_bar(
      sr,
      label_col = "plot_label",
      value_col = "abs_delta",
      output_prefix = file.path(args$output_dir, "cellchat_sender_receiver_top"),
      title = "Top sender-receiver CellChat delta burdens",
      xlab = "Summed absolute delta probability",
      top_n = args$top_n,
      fill = "#59A14F"
    )
  }

  if (nrow(consensus_df) > 0L) {
    consensus_df$mean_delta_score <- suppressWarnings(as.numeric(consensus_df$mean_delta_score))
    consensus_df$plot_label <- paste(consensus_df$contrast_id, consensus_df$level, consensus_df$sender, "->", consensus_df$receiver, paste0("[", consensus_df$ligand, "-", consensus_df$receptor, "]"))
    top_consensus <- consensus_df[order(consensus_df$n_methods, consensus_df$mean_delta_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    top_consensus <- head(top_consensus, args$top_n)
    write_tsv_local(consensus_df, file.path(args$output_dir, "consensus_interactions_all.tsv.gz"))
    write_tsv_local(top_consensus, file.path(args$output_dir, "consensus_interactions_top.tsv.gz"))
    cci_plot_bar(
      top_consensus,
      label_col = "plot_label",
      value_col = "n_methods",
      output_prefix = file.path(args$output_dir, "consensus_support_top"),
      title = "Top multi-method communication consensus interactions",
      xlab = "Number of supporting methods",
      top_n = args$top_n,
      fill = "#F28E2B"
    )
  }

  summary <- list(
    status = "ok",
    n_cellchat_rows = if (exists("cellchat_df")) nrow(cellchat_df) else 0L,
    n_consensus_rows = if (exists("consensus_df")) nrow(consensus_df) else 0L
  )
  write_json_local(summary, file.path(args$output_dir, "summary.json"))
  message(sprintf("[06b] cellchat_rows=%d consensus_rows=%d", summary$n_cellchat_rows, summary$n_consensus_rows))
}

main()
