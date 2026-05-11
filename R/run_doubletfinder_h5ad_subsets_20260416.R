#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DoubletFinder)
  library(reticulate)
})

options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v3")
Sys.setenv(PYTHONNOUSERSITE = "1")

DEFAULT_H5AD_PATHS <- c(
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/b_cells.h5ad",
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/endothelial_cells.h5ad",
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/epithelial_cells.h5ad",
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/fibroblast_cells.h5ad",
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/myeloid_cells.h5ad",
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/smc_cells.h5ad",
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/t_cells.h5ad"
)
DEFAULT_OUTPUT_DIR <- file.path(
  "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets",
  paste0("doubletfinder_", format(Sys.Date(), "%Y%m%d"))
)
LOCAL_GETSEURAT <- "/home/h2048/script/R/GetSeurat.R"
PREFERRED_PYTHON <- "/home/h2048/miniconda3/envs/bbknn_env/bin/python"

DOUBLETFINDER_RATE <- 0.06
MIN_CELLS_PER_SAMPLE <- 50L
MAX_PCS <- 30L
DEFAULT_PK <- 0.09
PK_SWEEP_MAX_CELLS <- 20000L
VERIFY_WRITTEN_H5AD <- TRUE
PREFERRED_ASSAY <- "RNA"
ENABLE_PK_SWEEP <- identical(Sys.getenv("DOUBLETFINDER_ENABLE_PK_SWEEP", "0"), "1")

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_DIR <- if (length(args) >= 1 && nzchar(args[[1]])) {
  normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
} else {
  DEFAULT_OUTPUT_DIR
}
H5AD_PATHS <- if (length(args) >= 2) args[-1] else DEFAULT_H5AD_PATHS
FULL_DIR <- file.path(OUTPUT_DIR, "full")
SINGLET_DIR <- file.path(OUTPUT_DIR, "singlets")
RDS_DIR <- file.path(OUTPUT_DIR, "rds")
SUMMARY_DIR <- file.path(OUTPUT_DIR, "summary")

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(fmt, ...)))
}

safe_chr <- function(x, fallback = "") {
  x <- as.character(x)
  x[is.na(x)] <- fallback
  x
}

normalize_meta_values <- function(x, fallback = "unknown") {
  x <- safe_chr(x, fallback = fallback)
  x <- trimws(x)
  x[!nzchar(x)] <- fallback
  x
}

get_assay_names_safe <- function(obj) {
  nm <- tryCatch(names(obj@assays), error = function(e) NULL)
  if (is.null(nm) || length(nm) == 0) {
    nm <- tryCatch(names(SeuratObject::Assays(obj)), error = function(e) character(0))
  }
  nm <- as.character(nm)
  nm[!is.na(nm) & nzchar(nm)]
}

pick_assay <- function(obj, preferred = "RNA") {
  assay_names <- get_assay_names_safe(obj)
  if (preferred %in% assay_names) return(preferred)
  active <- tryCatch(DefaultAssay(obj), error = function(e) NA_character_)
  if (!is.na(active) && nzchar(active) && active %in% assay_names) return(active)
  if (length(assay_names) == 0) stop("No assay found in Seurat object")
  assay_names[[1]]
}

get_assay_matrix <- function(obj, assay = "RNA", layer = c("counts", "data")) {
  layer <- match.arg(layer)
  mat <- tryCatch(
    LayerData(obj, assay = assay, layer = layer),
    error = function(e1) {
      tryCatch(
        GetAssayData(obj, assay = assay, slot = layer),
        error = function(e2) NULL
      )
    }
  )
  mat
}

fix_dimnames_counts <- function(counts_mat, obj, assay = "RNA") {
  if (is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
    feats <- tryCatch(SeuratObject::Features(obj, assay = assay), error = function(e) rownames(obj))
    cells <- colnames(obj)
    dimnames(counts_mat) <- list(feats, cells)
  }
  counts_mat
}

ensure_joined_layers <- function(obj, assay = "RNA") {
  assay_obj <- obj[[assay]]
  if (!inherits(assay_obj, "Assay5")) return(obj)
  layer_names <- tryCatch(Layers(assay_obj), error = function(e) character(0))
  if (length(layer_names) <= 1) return(obj)
  obj[[assay]] <- JoinLayers(obj[[assay]])
  obj
}

rebuild_assay_from_scratch <- function(obj, assay = "RNA") {
  if (!inherits(obj[[assay]], "Assay5")) return(obj)
  obj <- ensure_joined_layers(obj, assay = assay)
  counts_mat <- get_assay_matrix(obj, assay = assay, layer = "counts")
  if (is.null(counts_mat)) {
    stop(sprintf("Cannot extract counts matrix for assay '%s'", assay))
  }
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
  }
  counts_mat <- fix_dimnames_counts(counts_mat, obj, assay = assay)
  obj[[assay]] <- CreateAssayObject(counts = counts_mat)
  obj
}

ensure_normalize_command <- function(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
) {
  DefaultAssay(obj) <- assay
  cmd_key <- paste0("NormalizeData.", assay)
  if (is.null(obj@commands) || is.null(obj@commands[[cmd_key]])) {
    obj <- NormalizeData(
      obj,
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      verbose = FALSE
    )
  }
  obj
}

resolve_doubletfinder_fn <- function() {
  if (exists("doubletFinder", mode = "function")) return(get("doubletFinder", mode = "function"))
  if (exists("doubletFinder_v3", mode = "function")) return(get("doubletFinder_v3", mode = "function"))
  stop("Could not find doubletFinder() or doubletFinder_v3() in current DoubletFinder package")
}

resolve_param_sweep_fn <- function() {
  if (exists("paramSweep", mode = "function")) return(get("paramSweep", mode = "function"))
  if (exists("paramSweep_v3", mode = "function")) return(get("paramSweep_v3", mode = "function"))
  stop("Could not find paramSweep() or paramSweep_v3() in current DoubletFinder package")
}

DOUBLETFINDER_FN <- resolve_doubletfinder_fn()
PARAM_SWEEP_FN <- resolve_param_sweep_fn()

find_optimal_pk <- function(obj, pcs_use) {
  if (!ENABLE_PK_SWEEP) {
    return(DEFAULT_PK)
  }
  if (ncol(obj) > PK_SWEEP_MAX_CELLS) {
    return(DEFAULT_PK)
  }
  tryCatch({
    sweep_res <- PARAM_SWEEP_FN(obj, PCs = seq_len(pcs_use), sct = FALSE)
    sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    pk <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
    rm(sweep_res, sweep_stats, bcmvn)
    gc(verbose = FALSE)
    if (!is.finite(pk) || is.na(pk) || pk < 0.01 || pk > 0.3) DEFAULT_PK else pk
  }, error = function(e) {
    DEFAULT_PK
  })
}

detect_sample_col <- function(meta_df) {
  candidates <- c(
    "sample", "sample_id", "orig.ident", "donor_id", "donor",
    "patient_id", "patient", "library_id", "batch", "study"
  )
  for (col in candidates) {
    if (!col %in% colnames(meta_df)) next
    vals <- normalize_meta_values(meta_df[[col]], fallback = "")
    vals <- vals[nzchar(vals)]
    if (length(vals) >= 1) return(col)
  }
  NULL
}

detect_homotypic_col <- function(meta_df) {
  candidates <- c(
    "ann_finest_level", "ann_level_4", "ann_level_3", "ann_level_2",
    "cell_type_L3", "cell_type_L2", "cell_type", "Annotation_2",
    "Annotation", "seurat_clusters"
  )
  for (col in candidates) {
    if (!col %in% colnames(meta_df)) next
    vals <- normalize_meta_values(meta_df[[col]], fallback = "Unknown")
    if (length(unique(vals)) >= 2) return(col)
  }
  NULL
}

sanitize_obs_for_h5ad <- function(df) {
  out <- df
  for (col in colnames(out)) {
    vec <- out[[col]]
    if (is.factor(vec)) {
      out[[col]] <- as.character(vec)
    } else if (is.logical(vec)) {
      out[[col]] <- ifelse(is.na(vec), NA_integer_, as.integer(vec))
    } else if (inherits(vec, c("POSIXct", "POSIXt"))) {
      out[[col]] <- format(vec, tz = "UTC", usetz = TRUE)
    } else if (inherits(vec, "Date")) {
      out[[col]] <- as.character(vec)
    } else if (is.list(vec)) {
      out[[col]] <- vapply(vec, function(x) {
        if (length(x) == 0 || all(is.na(x))) return("")
        paste(as.character(x), collapse = "; ")
      }, character(1))
    } else if (is.character(vec)) {
      out[[col]] <- trimws(vec)
    } else if (!(is.integer(vec) || is.numeric(vec))) {
      out[[col]] <- as.character(vec)
    }

    if (is.character(out[[col]])) {
      out[[col]][is.na(out[[col]])] <- ""
      out[[col]] <- trimws(out[[col]])
    }
  }
  out
}

matrix_to_scipy_csr <- function(mat, scipy_sparse, np) {
  mat_csc <- as(mat, "dgCMatrix")
  scipy_sparse$csc_matrix(
    reticulate::tuple(
      np$array(as.numeric(mat_csc@x), dtype = np$float32),
      np$array(as.integer(mat_csc@i), dtype = np$int32),
      np$array(as.integer(mat_csc@p), dtype = np$int32)
    ),
    shape = reticulate::tuple(as.integer(nrow(mat_csc)), as.integer(ncol(mat_csc)))
  )$transpose()$tocsr()
}

write_seurat_to_h5ad <- function(obj, out_path, assay = "RNA") {
  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  try({
    anndata$settings$allow_write_nullable_strings <- TRUE
  }, silent = TRUE)

  counts_mat <- get_assay_matrix(obj, assay = assay, layer = "counts")
  if (is.null(counts_mat)) {
    stop(sprintf("Failed to extract counts matrix for export from assay '%s'", assay))
  }
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
  }
  counts_mat <- fix_dimnames_counts(counts_mat, obj, assay = assay)

  meta_export <- sanitize_obs_for_h5ad(obj@meta.data)
  meta_export <- meta_export[colnames(obj), , drop = FALSE]
  var_export <- data.frame(gene_symbol = rownames(obj), row.names = rownames(obj), stringsAsFactors = FALSE)

  expr_scipy <- matrix_to_scipy_csr(counts_mat, scipy_sparse, np)
  obs_py <- reticulate::r_to_py(meta_export)
  var_py <- reticulate::r_to_py(var_export)

  adata <- anndata$AnnData(X = expr_scipy, obs = obs_py, var = var_py)
  adata$layers$`__setitem__`("counts", expr_scipy$copy())

  for (red in Reductions(obj)) {
    emb <- tryCatch(Embeddings(obj, reduction = red), error = function(e) NULL)
    if (is.null(emb) || nrow(emb) != ncol(obj)) next
    emb <- emb[colnames(obj), , drop = FALSE]
    adata$obsm$`__setitem__`(paste0("X_", red), np$array(unname(emb), dtype = np$float32))
  }

  adata$write_h5ad(out_path, compression = "gzip")

  if (VERIFY_WRITTEN_H5AD) {
    adata_check <- anndata$read_h5ad(out_path, backed = "r")
    on.exit(try(adata_check$file$close(), silent = TRUE), add = TRUE)
    if (reticulate::py_to_r(adata_check$n_obs) != ncol(obj)) {
      stop(sprintf("Verification failed for %s: n_obs mismatch", basename(out_path)))
    }
    if (reticulate::py_to_r(adata_check$n_vars) != nrow(obj)) {
      stop(sprintf("Verification failed for %s: n_vars mismatch", basename(out_path)))
    }
  }

  invisible(out_path)
}

run_doubletfinder_single_sample <- function(obj, assay = "RNA", sample_name = "sample", homotypic_col = NULL) {
  cell_ids <- colnames(obj)
  if (length(cell_ids) < MIN_CELLS_PER_SAMPLE) {
    return(list(
      class = stats::setNames(rep("Skipped_too_few_cells", length(cell_ids)), cell_ids),
      pann = stats::setNames(rep(NA_real_, length(cell_ids)), cell_ids),
      status = "Skipped_too_few_cells",
      pK = NA_real_,
      nExp_poi = NA_integer_,
      nExp_adj = NA_integer_,
      homotypic_col = ifelse(is.null(homotypic_col), "", homotypic_col),
      homotypic_prop = NA_real_
    ))
  }

  DefaultAssay(obj) <- assay
  obj <- rebuild_assay_from_scratch(obj, assay = assay)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- ensure_normalize_command(obj, assay = assay)

  nfeatures_use <- min(2000L, max(50L, nrow(obj) - 1L))
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures_use, verbose = FALSE)
  if (length(VariableFeatures(obj)) < 50) {
    return(list(
      class = stats::setNames(rep("Skipped_low_feature_variance", length(cell_ids)), cell_ids),
      pann = stats::setNames(rep(NA_real_, length(cell_ids)), cell_ids),
      status = "Skipped_low_feature_variance",
      pK = NA_real_,
      nExp_poi = NA_integer_,
      nExp_adj = NA_integer_,
      homotypic_col = ifelse(is.null(homotypic_col), "", homotypic_col),
      homotypic_prop = NA_real_
    ))
  }

  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

  pcs_use <- min(
    as.integer(MAX_PCS),
    as.integer(max(5L, ncol(obj) - 1L)),
    as.integer(max(5L, length(VariableFeatures(obj)) - 1L))
  )
  if (!is.finite(pcs_use) || is.na(pcs_use) || pcs_use < 5) {
    return(list(
      class = stats::setNames(rep("Skipped_insufficient_pcs", length(cell_ids)), cell_ids),
      pann = stats::setNames(rep(NA_real_, length(cell_ids)), cell_ids),
      status = "Skipped_insufficient_pcs",
      pK = NA_real_,
      nExp_poi = NA_integer_,
      nExp_adj = NA_integer_,
      homotypic_col = ifelse(is.null(homotypic_col), "", homotypic_col),
      homotypic_prop = NA_real_
    ))
  }

  obj <- RunPCA(obj, npcs = pcs_use, features = VariableFeatures(obj), verbose = FALSE)
  pK <- find_optimal_pk(obj, pcs_use)
  nExp_poi <- max(1L, round(DOUBLETFINDER_RATE * ncol(obj)))

  homotypic_prop <- NA_real_
  if (!is.null(homotypic_col) && homotypic_col %in% colnames(obj@meta.data)) {
    annotations <- normalize_meta_values(obj@meta.data[[homotypic_col]], fallback = "Unknown")
    if (length(unique(annotations)) >= 2) {
      homotypic_prop <- tryCatch(modelHomotypic(annotations), error = function(e) NA_real_)
    }
  }
  nExp_adj <- if (is.finite(homotypic_prop) && !is.na(homotypic_prop)) {
    max(1L, round(nExp_poi * (1 - homotypic_prop)))
  } else {
    nExp_poi
  }

  nExp_final <- if (is.finite(nExp_adj) && !is.na(nExp_adj)) {
    nExp_adj
  } else {
    nExp_poi
  }

  obj <- DOUBLETFINDER_FN(
    obj,
    PCs = seq_len(pcs_use),
    pN = 0.25,
    pK = pK,
    nExp = nExp_final,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  final_pann_col <- tail(grep("^pANN_", colnames(obj@meta.data), value = TRUE), 1)
  final_class_col <- tail(grep("^DF.classifications_", colnames(obj@meta.data), value = TRUE), 1)
  if (!length(final_class_col)) {
    stop(sprintf("DoubletFinder finished for sample '%s' but no DF.classifications_* column was found", sample_name))
  }

  final_class <- as.character(obj@meta.data[[final_class_col]])
  final_class <- ifelse(final_class == "Singlet", "Singlet", "Doublet")
  names(final_class) <- colnames(obj)

  final_pann <- if (length(final_pann_col)) as.numeric(obj@meta.data[[final_pann_col]]) else rep(NA_real_, ncol(obj))
  names(final_pann) <- colnames(obj)

  list(
    class = final_class,
    pann = final_pann,
    status = "Success",
    pK = pK,
    nExp_poi = nExp_poi,
    nExp_adj = nExp_adj,
    homotypic_col = ifelse(is.null(homotypic_col), "", homotypic_col),
    homotypic_prop = homotypic_prop
  )
}

process_one_h5ad <- function(h5ad_path) {
  h5ad_path <- normalizePath(h5ad_path, winslash = "/", mustWork = TRUE)
  file_stem <- tools::file_path_sans_ext(basename(h5ad_path))
  log_msg("Loading %s", basename(h5ad_path))

  obj <- GetSeurat(
    h5ad_path = h5ad_path,
    assay = PREFERRED_ASSAY,
    prefer_raw = TRUE,
    prefer_layer_counts = TRUE,
    validate_counts = TRUE,
    debug = FALSE
  )
  assay_use <- pick_assay(obj, preferred = PREFERRED_ASSAY)
  DefaultAssay(obj) <- assay_use

  sample_col <- detect_sample_col(obj@meta.data)
  if (is.null(sample_col)) {
    obj$doubletfinder_sample_proxy <- "all_cells"
    sample_col <- "doubletfinder_sample_proxy"
  }

  sample_values <- normalize_meta_values(obj@meta.data[[sample_col]], fallback = "missing_sample")
  n_samples <- length(unique(sample_values))
  homotypic_col <- detect_homotypic_col(obj@meta.data)

  log_msg("  [%s] cells=%d | sample_col=%s | n_samples=%d | homotypic_col=%s",
          file_stem,
          ncol(obj),
          sample_col,
          n_samples,
          ifelse(is.null(homotypic_col), "<none>", homotypic_col))

  sample_levels <- unique(sample_values)
  sample_groups <- lapply(sample_levels, function(sample_id) {
    colnames(obj)[sample_values == sample_id]
  })
  names(sample_groups) <- sample_levels

  final_class <- stats::setNames(rep(NA_character_, ncol(obj)), colnames(obj))
  final_pann <- stats::setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  final_status <- stats::setNames(rep(NA_character_, ncol(obj)), colnames(obj))
  final_pK <- stats::setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  final_nExp_poi <- stats::setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  final_nExp_adj <- stats::setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  final_homotypic_col <- stats::setNames(rep(ifelse(is.null(homotypic_col), "", homotypic_col), ncol(obj)), colnames(obj))
  final_homotypic_prop <- stats::setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  summary_rows <- vector("list", length(sample_groups))

  for (idx in seq_along(sample_groups)) {
    sample_name <- names(sample_groups)[[idx]]
    sample_cells <- sample_groups[[idx]]

    log_msg("  [%s] sample %d/%d | sample=%s | cells=%d",
            file_stem,
            idx,
            length(sample_groups),
            sample_name,
            length(sample_cells))

    sample_obj <- subset(obj, cells = sample_cells)
    df_res <- tryCatch(
      run_doubletfinder_single_sample(
        sample_obj,
        assay = assay_use,
        sample_name = sample_name,
        homotypic_col = homotypic_col
      ),
      error = function(e) {
        list(
          class = stats::setNames(rep("DF_error", length(sample_cells)), sample_cells),
          pann = stats::setNames(rep(NA_real_, length(sample_cells)), sample_cells),
          status = paste0("Error: ", conditionMessage(e)),
          pK = NA_real_,
          nExp_poi = NA_integer_,
          nExp_adj = NA_integer_,
          homotypic_col = ifelse(is.null(homotypic_col), "", homotypic_col),
          homotypic_prop = NA_real_
        )
      }
    )

    final_class[sample_cells] <- df_res$class[sample_cells]
    final_pann[sample_cells] <- df_res$pann[sample_cells]
    final_status[sample_cells] <- df_res$status
    final_pK[sample_cells] <- df_res$pK
    final_nExp_poi[sample_cells] <- df_res$nExp_poi
    final_nExp_adj[sample_cells] <- df_res$nExp_adj
    final_homotypic_col[sample_cells] <- df_res$homotypic_col
    final_homotypic_prop[sample_cells] <- df_res$homotypic_prop

    summary_rows[[idx]] <- data.frame(
      file_stem = file_stem,
      h5ad_path = h5ad_path,
      sample = sample_name,
      sample_col = sample_col,
      n_samples = n_samples,
      homotypic_col = ifelse(is.null(df_res$homotypic_col), "", df_res$homotypic_col),
      n_cells = length(sample_cells),
      n_singlets = sum(df_res$class == "Singlet", na.rm = TRUE),
      n_doublets = sum(df_res$class == "Doublet", na.rm = TRUE),
      status = df_res$status,
      pK = df_res$pK,
      nExp_poi = df_res$nExp_poi,
      nExp_adj = df_res$nExp_adj,
      homotypic_prop = df_res$homotypic_prop,
      stringsAsFactors = FALSE
    )

    rm(sample_obj, df_res)
    gc(verbose = FALSE)
  }

  obj$doubletfinder_class <- final_class[colnames(obj)]
  obj$doubletfinder_pANN <- final_pann[colnames(obj)]
  obj$doubletfinder_status <- final_status[colnames(obj)]
  obj$doubletfinder_pK <- final_pK[colnames(obj)]
  obj$doubletfinder_nExp_poi <- final_nExp_poi[colnames(obj)]
  obj$doubletfinder_nExp_adj <- final_nExp_adj[colnames(obj)]
  obj$doubletfinder_homotypic_col <- final_homotypic_col[colnames(obj)]
  obj$doubletfinder_homotypic_prop <- final_homotypic_prop[colnames(obj)]
  obj$doubletfinder_is_singlet <- obj$doubletfinder_class == "Singlet"

  summary_df <- do.call(rbind, summary_rows)
  summary_path <- file.path(SUMMARY_DIR, paste0(file_stem, "_doubletfinder_summary.csv"))
  write.csv(summary_df, summary_path, row.names = FALSE)

  full_rds_path <- file.path(RDS_DIR, paste0(file_stem, "_doubletfinder_annotated.rds"))
  full_h5ad_path <- file.path(FULL_DIR, paste0(file_stem, "_doubletfinder_annotated.h5ad"))
  saveRDS(obj, full_rds_path)
  write_seurat_to_h5ad(obj, full_h5ad_path, assay = assay_use)

  singlet_cells <- colnames(obj)[obj$doubletfinder_class == "Singlet"]
  singlet_rds_path <- NA_character_
  singlet_h5ad_path <- NA_character_
  if (length(singlet_cells) > 0) {
    singlet_obj <- subset(obj, cells = singlet_cells)
    singlet_rds_path <- file.path(RDS_DIR, paste0(file_stem, "_doubletfinder_singlets.rds"))
    singlet_h5ad_path <- file.path(SINGLET_DIR, paste0(file_stem, "_doubletfinder_singlets.h5ad"))
    saveRDS(singlet_obj, singlet_rds_path)
    write_seurat_to_h5ad(singlet_obj, singlet_h5ad_path, assay = assay_use)
  }

  file_manifest <- data.frame(
    file_stem = file_stem,
    input_h5ad = h5ad_path,
    output_summary_csv = summary_path,
    output_full_rds = full_rds_path,
    output_full_h5ad = full_h5ad_path,
    output_singlet_rds = singlet_rds_path,
    output_singlet_h5ad = singlet_h5ad_path,
    total_cells = ncol(obj),
    total_singlets = sum(obj$doubletfinder_class == "Singlet", na.rm = TRUE),
    total_doublets = sum(obj$doubletfinder_class == "Doublet", na.rm = TRUE),
    sample_col = sample_col,
    stringsAsFactors = FALSE
  )

  manifest_path <- file.path(SUMMARY_DIR, paste0(file_stem, "_doubletfinder_manifest.csv"))
  write.csv(file_manifest, manifest_path, row.names = FALSE)

  log_msg("Finished %s | singlets=%d | doublets=%d",
          file_stem,
          sum(obj$doubletfinder_class == "Singlet", na.rm = TRUE),
          sum(obj$doubletfinder_class == "Doublet", na.rm = TRUE))

  summary_df
}

configure_python <- function() {
  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
  } else if (file.exists(PREFERRED_PYTHON)) {
    Sys.setenv(RETICULATE_PYTHON = PREFERRED_PYTHON)
    use_python(PREFERRED_PYTHON, required = TRUE)
  } else {
    use_condaenv("bbknn_env", required = TRUE)
  }
  invisible(py_config())
}

main <- function() {
  if (!file.exists(LOCAL_GETSEURAT)) {
    stop(sprintf("GetSeurat helper not found: %s", LOCAL_GETSEURAT))
  }
  source(LOCAL_GETSEURAT)
  configure_python()

  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  for (subdir in c(FULL_DIR, SINGLET_DIR, RDS_DIR, SUMMARY_DIR)) {
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  }

  log_msg("Python configured: %s", py_config()$python)
  log_msg("Output directory: %s", OUTPUT_DIR)
  log_msg("Processing %d h5ad files", length(H5AD_PATHS))

  combined_rows <- list()
  for (idx in seq_along(H5AD_PATHS)) {
    h5ad_path <- H5AD_PATHS[[idx]]
    if (!file.exists(h5ad_path)) {
      log_msg("[SKIP] Missing file: %s", h5ad_path)
      combined_rows[[length(combined_rows) + 1L]] <- data.frame(
        file_stem = tools::file_path_sans_ext(basename(h5ad_path)),
        h5ad_path = h5ad_path,
        sample = NA_character_,
        sample_col = NA_character_,
        homotypic_col = NA_character_,
        n_cells = NA_integer_,
        n_singlets = NA_integer_,
        n_doublets = NA_integer_,
        status = "Missing file",
        pK = NA_real_,
        nExp_poi = NA_integer_,
        nExp_adj = NA_integer_,
        homotypic_prop = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    combined_rows[[length(combined_rows) + 1L]] <- process_one_h5ad(h5ad_path)
  }

  combined_summary <- do.call(rbind, combined_rows)
  combined_path <- file.path(OUTPUT_DIR, "summary", "combined_doubletfinder_summary.csv")
  write.csv(combined_summary, combined_path, row.names = FALSE)
  log_msg("Combined summary saved: %s", combined_path)
}

main()
