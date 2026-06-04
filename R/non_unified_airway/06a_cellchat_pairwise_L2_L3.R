#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

HELPER_PATH <- "/home/h2048/script/R/cell_communication_helper_20260506_v1_0.R"
source(HELPER_PATH)

DEFAULT_PREPARED_MANIFEST <- "/home/h2048/data/py/0525/non_unified_airway/communication/prepared_inputs_manifest.tsv"
DEFAULT_OUTPUT_DIR <- "/home/h2048/data/R/0525/non_unified_airway/communication/cellchat"
DEFAULT_PREP_SCRIPT <- "/home/h2048/script/py/non_unified_airway/05a_prepare_cell_communication_inputs.py"

parse_args_simple <- function(args) {
  out <- list(
    prepared_manifest = DEFAULT_PREPARED_MANIFEST,
    output_dir = DEFAULT_OUTPUT_DIR,
    prep_script = DEFAULT_PREP_SCRIPT,
    prepare_if_missing = TRUE,
    force = FALSE,
    nboot = 100L,
    workers = 4L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--prepared-manifest" && i < length(args)) {
      out$prepared_manifest <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--output-dir" && i < length(args)) {
      out$output_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--prep-script" && i < length(args)) {
      out$prep_script <- args[[i + 1L]]
      i <- i + 2L
    } else if (key == "--nboot" && i < length(args)) {
      out$nboot <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--workers" && i < length(args)) {
      out$workers <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (key == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else if (key == "--no-prepare-if-missing") {
      out$prepare_if_missing <- FALSE
      i <- i + 1L
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

ensure_prepared_manifest <- function(prepared_manifest, prep_script, prepare_if_missing = TRUE) {
  if (file.exists(prepared_manifest)) return(invisible(prepared_manifest))
  if (!isTRUE(prepare_if_missing)) {
    stop(sprintf("Prepared manifest does not exist: %s", prepared_manifest), call. = FALSE)
  }
  message("[06a] Prepared manifest missing; invoking Python 05a export helper...")
  cmd <- c(
    prep_script,
    "--export-h5ad",
    "--export-mtx",
    "--out-dir",
    dirname(prepared_manifest)
  )
  status <- system2("python3", cmd)
  if (!identical(status, 0L) || !file.exists(prepared_manifest)) {
    stop(sprintf("Failed to generate prepared manifest via %s", prep_script), call. = FALSE)
  }
  invisible(prepared_manifest)
}

read_mtx_bundle <- function(counts_mtx, genes_tsv, barcodes_tsv, metadata_tsv) {
  counts <- Matrix::readMM(counts_mtx)
  genes <- read_tsv_auto(genes_tsv)
  barcodes <- read_tsv_auto(barcodes_tsv)
  meta <- read_tsv_auto(metadata_tsv)

  gene_label <- if ("gene_symbol" %in% colnames(genes)) genes$gene_symbol else genes$gene_id
  gene_label <- as.character(gene_label)
  gene_label[!nzchar(gene_label)] <- as.character(genes$gene_id[!nzchar(gene_label)])
  rownames(counts) <- make.unique(gene_label)
  colnames(counts) <- as.character(barcodes$cell_id)

  meta$cell_id <- as.character(meta$cell_id)
  meta <- meta[match(colnames(counts), meta$cell_id), , drop = FALSE]
  rownames(meta) <- meta$cell_id
  list(counts = counts, metadata = meta)
}

make_seurat_bundle <- function(bundle) {
  seu <- Seurat::CreateSeuratObject(counts = bundle$counts, meta.data = bundle$metadata, assay = "RNA")
  seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  seu
}

run_cellchat_side <- function(seu, side, nboot = 100L, workers = 4L) {
  sub_obj <- subset(seu, subset = contrast_side == side)
  meta <- sub_obj@meta.data
  if (nrow(meta) == 0L) return(NULL)
  if (length(unique(meta$comm_celltype)) < 2L) return(NULL)
  cc <- cci_prepare_cellchat(
    sub_obj,
    celltype_col = "comm_celltype",
    species = "human",
    assay = "RNA",
    slot_use = "data",
    min_cells = 10L
  )
  cci_run_cellchat_screening(cc, nboot = nboot, workers = workers)
}

standardize_candidates <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(data.frame(
      sender = character(), receiver = character(), ligand = character(), receptor = character(),
      interaction_key = character(), prob = numeric(), stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  out <- df
  out$sender <- as.character(out$source)
  out$receiver <- as.character(out$target)
  if (!"ligand" %in% colnames(out)) out$ligand <- ""
  if (!"receptor" %in% colnames(out)) out$receptor <- ""
  out$ligand <- as.character(out$ligand)
  out$receptor <- as.character(out$receptor)
  out$interaction_key <- paste(out$sender, out$receiver, out$ligand, out$receptor, sep = "|")
  out
}

build_candidate_delta <- function(left_df, right_df, contrast_id, level, left_label, right_label) {
  left_std <- standardize_candidates(left_df)
  right_std <- standardize_candidates(right_df)
  left_keep <- intersect(c("interaction_key", "sender", "receiver", "ligand", "receptor", "prob", "pval", "interaction_name"), colnames(left_std))
  right_keep <- intersect(c("interaction_key", "sender", "receiver", "ligand", "receptor", "prob", "pval", "interaction_name"), colnames(right_std))
  merged <- merge(
    left_std[, left_keep, drop = FALSE],
    right_std[, right_keep, drop = FALSE],
    by = c("interaction_key", "sender", "receiver", "ligand", "receptor"),
    all = TRUE,
    suffixes = c("_left", "_right")
  )
  merged$left_score <- suppressWarnings(as.numeric(merged$prob_left))
  merged$right_score <- suppressWarnings(as.numeric(merged$prob_right))
  merged$delta_score_right_minus_left <- ifelse(is.na(merged$right_score), 0, merged$right_score) - ifelse(is.na(merged$left_score), 0, merged$left_score)
  merged$method <- "cellchat"
  merged$contrast_id <- contrast_id
  merged$level <- level
  merged$left_label <- left_label
  merged$right_label <- right_label
  merged <- merged[order(merged$delta_score_right_minus_left, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  rownames(merged) <- NULL
  merged
}

extract_pathway_strength <- function(cc_obj, side, contrast_id, level) {
  if (is.null(cc_obj)) return(data.frame())
  prob <- tryCatch(cc_obj@netP$prob, error = function(e) NULL)
  if (is.null(prob) || length(dim(prob)) != 3L) return(data.frame())
  pathways <- dimnames(prob)[[3]]
  if (is.null(pathways)) pathways <- paste0("pathway_", seq_len(dim(prob)[3]))
  rows <- lapply(seq_along(pathways), function(idx) {
    slice <- prob[, , idx, drop = TRUE]
    data.frame(
      contrast_id = contrast_id,
      level = level,
      side = side,
      pathway = pathways[[idx]],
      total_prob = sum(slice, na.rm = TRUE),
      n_edges = sum(slice > 0, na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

run_one_prepared_row <- function(row, output_dir, force = FALSE, nboot = 100L, workers = 4L) {
  contrast_id <- as.character(row$contrast_id)
  level <- as.character(row$level)
  run_dir <- file.path(output_dir, contrast_id, level)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  status_path <- file.path(run_dir, "cellchat_status.json")
  result_path <- file.path(run_dir, "cellchat_delta.tsv.gz")
  if (file.exists(result_path) && !isTRUE(force)) {
    status <- list(contrast_id = contrast_id, level = level, status = "skipped_exists", result = result_path)
    write_json_local(status, status_path)
    return(status)
  }

  required_paths <- c(row$counts_mtx, row$genes_tsv, row$barcodes_tsv, row$cell_metadata_tsv)
  if (any(!file.exists(required_paths))) {
    status <- list(
      contrast_id = contrast_id,
      level = level,
      status = "blocked_missing_mtx_bundle",
      missing = required_paths[!file.exists(required_paths)]
    )
    write_json_local(status, status_path)
    return(status)
  }

  tryCatch({
    bundle <- read_mtx_bundle(row$counts_mtx, row$genes_tsv, row$barcodes_tsv, row$cell_metadata_tsv)
    seu <- make_seurat_bundle(bundle)
    left_cc <- run_cellchat_side(seu, "left", nboot = nboot, workers = workers)
    right_cc <- run_cellchat_side(seu, "right", nboot = nboot, workers = workers)

    left_candidates <- if (!is.null(left_cc)) cci_extract_cellchat_candidates(left_cc) else data.frame()
    right_candidates <- if (!is.null(right_cc)) cci_extract_cellchat_candidates(right_cc) else data.frame()
    delta_df <- build_candidate_delta(left_candidates, right_candidates, contrast_id, level, as.character(row$left_label), as.character(row$right_label))

    left_pathway <- extract_pathway_strength(left_cc, "left", contrast_id, level)
    right_pathway <- extract_pathway_strength(right_cc, "right", contrast_id, level)
    pathway_df <- rbind(left_pathway, right_pathway)

    if (!is.null(left_cc)) saveRDS(left_cc, file.path(run_dir, "cellchat_left.rds"))
    if (!is.null(right_cc)) saveRDS(right_cc, file.path(run_dir, "cellchat_right.rds"))
    write_tsv_local(left_candidates, file.path(run_dir, "cellchat_left_candidates.tsv.gz"))
    write_tsv_local(right_candidates, file.path(run_dir, "cellchat_right_candidates.tsv.gz"))
    write_tsv_local(delta_df, result_path)
    write_tsv_local(pathway_df, file.path(run_dir, "cellchat_pathway_strength.tsv.gz"))

    status <- list(
      contrast_id = contrast_id,
      level = level,
      status = "ok",
      n_left = nrow(left_candidates),
      n_right = nrow(right_candidates),
      n_delta = nrow(delta_df),
      result = result_path
    )
    write_json_local(status, status_path)
    status
  }, error = function(e) {
    status <- list(
      contrast_id = contrast_id,
      level = level,
      status = "error",
      error = conditionMessage(e)
    )
    write_json_local(status, status_path)
    status
  })
}

main <- function() {
  args <- parse_args_simple(commandArgs(trailingOnly = TRUE))
  ensure_prepared_manifest(args$prepared_manifest, args$prep_script, args$prepare_if_missing)
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- read_tsv_auto(args$prepared_manifest)
  manifest <- manifest[manifest$status == "ready", , drop = FALSE]
  if (nrow(manifest) == 0L) {
    write_json_local(list(status = "blocked_no_ready_inputs", n_inputs = 0L), file.path(args$output_dir, "run_summary.json"))
    message("[06a] No ready inputs found in prepared manifest.")
    return(invisible(NULL))
  }

  statuses <- lapply(seq_len(nrow(manifest)), function(idx) {
    run_one_prepared_row(manifest[idx, , drop = FALSE], args$output_dir, force = args$force, nboot = args$nboot, workers = args$workers)
  })
  status_df <- do.call(rbind, lapply(statuses, as.data.frame, stringsAsFactors = FALSE))
  write_tsv_local(status_df, file.path(args$output_dir, "cellchat_status.tsv.gz"))
  summary <- list(
    status = "ok",
    n_inputs = nrow(manifest),
    n_ok = sum(status_df$status == "ok"),
    n_error = sum(status_df$status == "error"),
    n_blocked = sum(grepl("^blocked", status_df$status))
  )
  write_json_local(summary, file.path(args$output_dir, "run_summary.json"))
  message(sprintf("[06a] inputs=%d ok=%d error=%d", summary$n_inputs, summary$n_ok, summary$n_error))
}

main()
