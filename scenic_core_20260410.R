#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# SCENIC Module for Seurat Objects
# ==============================================================================
# Purpose:
#   1. Read final Seurat object
#   2. Run SCENIC regulon inference on stratified downsampled cells
#   3. Score regulon activity (AUCell) on full object
#   4. Write continuous regulon AUC back to Seurat as SCENIC assay
#   5. Export heatmaps / RSS / network / tables
#
# Best practice in this version:
#   - Network inference: downsampled subset only
#   - AUCell scoring   : full object
#   - Seurat v4/v5 compatible counts extraction
#
# Author : OpenAI
# Date   : 2026-04-08
# ==============================================================================

# ==============================================================================
# 0. Thread Control
# ==============================================================================
Sys.setenv(
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS  = "1"
)

# ==============================================================================
# 1. Packages
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SCopeLoomR)
  library(BiocParallel)
  library(visNetwork)
})

# ==============================================================================
# 2. Configuration
# ==============================================================================

# ----- Input / Output -----
if (!exists("INPUT_RDS")) {
  INPUT_RDS <- "Basal_analyzed.rds"
}
if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- "./SCENIC_Basal"
}
if (!exists("DATABASE_DIR")) {
  compatible_db_dir <- "/home/h2048/data/index_genome/cisTarget_databases_rscenic"
  raw_db_dir <- "/home/h2048/data/index_genome/cisTarget_databases"
  DATABASE_DIR <- if (dir.exists(compatible_db_dir)) compatible_db_dir else raw_db_dir
}
if (!exists("SCENIC_DB_10KB")) {
  SCENIC_DB_10KB <- NULL
}
if (!exists("SCENIC_DB_500BP")) {
  SCENIC_DB_500BP <- NULL
}
if (!exists("SCENIC_DB_INDEX_COL")) {
  SCENIC_DB_INDEX_COL <- NULL
}
if (!exists("ORGANISM")) {
  ORGANISM <- "hgnc"                        # human = hgnc, mouse = mgi
}

# ----- Seurat columns -----
if (!exists("ASSAY_USE")) {
  ASSAY_USE <- "RNA"
}
if (!exists("CELL_TYPE_COL")) {
  CELL_TYPE_COL <- "Annotation"              # <-- modify
}
if (!exists("SAMPLE_COL")) {
  SAMPLE_COL <- "sample"                     # <-- modify if needed
}

# ----- Basic filters -----
if (!exists("MIN_GENES_PER_CELL")) {
  MIN_GENES_PER_CELL <- 200
}
if (!exists("MIN_CELLS_PER_GENE")) {
  MIN_CELLS_PER_GENE <- 10
}
if (!exists("MIN_GENE_PCT")) {
  MIN_GENE_PCT <- 0.01
}

# ----- SCENIC inference subset -----
if (!exists("INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE")) {
  INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE <- 200
}
if (!exists("INFERENCE_MAX_CELLS_PER_CELLTYPE")) {
  INFERENCE_MAX_CELLS_PER_CELLTYPE <- 3000
}
if (!exists("INFERENCE_GLOBAL_MAX_CELLS")) {
  INFERENCE_GLOBAL_MAX_CELLS <- 30000
}

# ----- SCENIC parameters -----
if (!exists("N_CORES")) {
  N_CORES <- 8
}
if (!exists("MIN_GENES_PER_REGULON")) {
  MIN_GENES_PER_REGULON <- 20
}
if (!exists("SCENIC_DATASET_TITLE")) {
  SCENIC_DATASET_TITLE <- "Seurat_SCENIC"
}
if (!exists("TOP_RSS_PER_CELLTYPE")) {
  TOP_RSS_PER_CELLTYPE <- 8
}
if (!exists("TOP_VAR_REGULONS_FOR_PLOTS")) {
  TOP_VAR_REGULONS_FOR_PLOTS <- 12
}
if (!exists("NETWORK_TOP_REGULONS")) {
  NETWORK_TOP_REGULONS <- 30
}
if (!exists("NETWORK_TOP_TARGETS_PER_REGULON")) {
  NETWORK_TOP_TARGETS_PER_REGULON <- 30
}

# ----- Plot options -----
if (!exists("HEATMAP_WIDTH")) {
  HEATMAP_WIDTH <- 12
}
if (!exists("HEATMAP_HEIGHT")) {
  HEATMAP_HEIGHT <- 10
}
if (!exists("FEATUREPLOT_PDF_WIDTH")) {
  FEATUREPLOT_PDF_WIDTH <- 12
}
if (!exists("FEATUREPLOT_PDF_HEIGHT")) {
  FEATUREPLOT_PDF_HEIGHT <- 10
}
if (!exists("REDUCTION_CANDIDATES")) {
  REDUCTION_CANDIDATES <- c("umap", "umap_harmony", "umap_scanvi", "harmony", "pca")
}
if (!exists("SCENIC_SOURCE_ONLY")) {
  SCENIC_SOURCE_ONLY <- FALSE
}
if (!exists("MAX_INFERENCE_DENSE_GB")) {
  MAX_INFERENCE_DENSE_GB <- 8
}
if (!exists("EXPORT_FULL_AUC_MATRIX_CSV")) {
  EXPORT_FULL_AUC_MATRIX_CSV <- FALSE
}
if (!exists("RANDOM_SEED")) {
  RANDOM_SEED <- 42
}

# ----- Reproducibility -----
set.seed(RANDOM_SEED)
options(stringsAsFactors = FALSE)
options(scipen = 999)

# ==============================================================================
# 3. Helpers
# ==============================================================================

DIR_INT   <- NULL
DIR_TABLE <- NULL
DIR_FIG   <- NULL
DIR_RDS   <- NULL

msg <- function(...) {
  txt <- sprintf(...)
  cat(txt)
  try(flush(stdout()), silent = TRUE)
  try(flush.console(), silent = TRUE)
  invisible(txt)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) {
    return(y)
  }
  x
}

initialize_output_dirs <- function(output_dir = OUTPUT_DIR) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  DIR_INT <<- file.path(output_dir, "int")
  DIR_TABLE <<- file.path(output_dir, "tables")
  DIR_FIG <<- file.path(output_dir, "figures")
  DIR_RDS <<- file.path(output_dir, "rds")

  for (d in c(DIR_INT, DIR_TABLE, DIR_FIG, DIR_RDS)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  invisible(list(
    int = DIR_INT,
    tables = DIR_TABLE,
    figures = DIR_FIG,
    rds = DIR_RDS
  ))
}

normalize_path_safe <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

hash_serialized_object <- function(x) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf)
  unname(tools::md5sum(tf))
}

vector_signature <- function(x) {
  x <- as.character(x)
  list(
    length = length(x),
    hash = hash_serialized_object(x)
  )
}

matrix_signature <- function(mat) {
  list(
    nrow = nrow(mat),
    ncol = ncol(mat),
    row_hash = hash_serialized_object(rownames(mat)),
    col_hash = hash_serialized_object(colnames(mat)),
    row_sum_hash = hash_serialized_object(as.numeric(Matrix::rowSums(mat))),
    col_sum_hash = hash_serialized_object(as.numeric(Matrix::colSums(mat)))
  )
}

regulon_signature <- function(regulons) {
  gene_sets <- lapply(regulons, extract_regulon_gene_vector)
  list(
    n_regulons = length(gene_sets),
    name_hash = hash_serialized_object(names(gene_sets)),
    size_hash = hash_serialized_object(as.integer(lengths(gene_sets))),
    gene_hash = hash_serialized_object(gene_sets)
  )
}

file_signature <- function(paths) {
  lapply(paths, function(path) {
    path_norm <- normalize_path_safe(path)
    info <- file.info(path_norm)
    list(
      path = path_norm,
      size = unname(info$size),
      mtime = as.character(info$mtime)
    )
  })
}

manifest_matches <- function(manifest_file, current_manifest) {
  if (!file.exists(manifest_file)) {
    return(FALSE)
  }
  previous_manifest <- tryCatch(readRDS(manifest_file), error = function(e) NULL)
  identical(previous_manifest, current_manifest)
}

write_manifest <- function(manifest_file, manifest) {
  saveRDS(manifest, manifest_file)
  invisible(manifest_file)
}

database_pair_key <- function(path) {
  x <- tolower(basename(path))
  x <- gsub("10kb_up_and_down_tss|500bp_up_and_100bp_down_tss", "{region}", x)
  x <- gsub("10kb|500bp", "{region}", x)
  x
}

resolve_scenic_database_files <- function(
  organism,
  db_dir,
  db_10kb = SCENIC_DB_10KB,
  db_500bp = SCENIC_DB_500BP
) {
  explicit_10kb <- !is.null(db_10kb) && nzchar(db_10kb)
  explicit_500bp <- !is.null(db_500bp) && nzchar(db_500bp)

  if (xor(explicit_10kb, explicit_500bp)) {
    stop("SCENIC_DB_10KB and SCENIC_DB_500BP must be set together.")
  }

  if (explicit_10kb && explicit_500bp) {
    selected <- c(
      "10kb" = normalize_path_safe(db_10kb),
      "500bp" = normalize_path_safe(db_500bp)
    )
  } else {
    db_files <- list.files(db_dir, pattern = "\\.feather$", full.names = TRUE)
    if (length(db_files) == 0) {
      stop(sprintf("No .feather SCENIC databases found in: %s", db_dir))
    }

    organism_pattern <- switch(
      organism,
      hgnc = "^hg[0-9]+__",
      mgi = "^mm[0-9]+__",
      ""
    )
    if (nzchar(organism_pattern)) {
      db_files <- db_files[grepl(organism_pattern, basename(db_files), ignore.case = TRUE)]
    }

    db_10kb_candidates <- db_files[grepl("10kb", basename(db_files), ignore.case = TRUE)]
    db_500bp_candidates <- db_files[grepl("500bp", basename(db_files), ignore.case = TRUE)]

    if (length(db_10kb_candidates) != 1 || length(db_500bp_candidates) != 1) {
      stop(sprintf(
        paste(
          "Database selection is ambiguous in '%s'.",
          "Found %d 10kb candidate(s) and %d 500bp candidate(s).",
          "Please set SCENIC_DB_10KB and SCENIC_DB_500BP explicitly."
        ),
        db_dir,
        length(db_10kb_candidates),
        length(db_500bp_candidates)
      ))
    }

    selected <- c(
      "10kb" = normalize_path_safe(db_10kb_candidates[[1]]),
      "500bp" = normalize_path_safe(db_500bp_candidates[[1]])
    )
  }

  if (any(!file.exists(selected))) {
    missing_files <- selected[!file.exists(selected)]
    stop(sprintf(
      "The following SCENIC database files do not exist: %s",
      paste(missing_files, collapse = ", ")
    ))
  }

  pair_keys <- vapply(selected, database_pair_key, character(1))
  if (length(unique(pair_keys)) != 1) {
    stop(sprintf(
      paste(
        "Selected 10kb/500bp databases do not appear to belong to the same database set:",
        "%s"
      ),
      paste(basename(selected), collapse = " vs ")
    ))
  }

  selected
}

resolve_db_index_col <- function(db_files, preferred = SCENIC_DB_INDEX_COL) {
  db_file <- unname(db_files[[1]])

  if (!is.null(preferred) && nzchar(preferred)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      return(preferred)
    }

    rf <- arrow::ReadableFile$create(db_file)
    fr <- arrow::FeatherReader$create(rf)
    on.exit(try(rf$close(), silent = TRUE), add = TRUE)
    col_names <- names(fr)

    if (!preferred %in% col_names) {
      stop(sprintf(
        "Configured SCENIC_DB_INDEX_COL '%s' not found in database '%s'. Available columns include: %s",
        preferred,
        basename(db_file),
        paste(utils::head(col_names, 10), collapse = ", ")
      ))
    }
    return(preferred)
  }

  if (!requireNamespace("arrow", quietly = TRUE)) {
    msg("[WARN] 'arrow' not available; defaulting SCENIC dbIndexCol to 'features'.\n")
    return("features")
  }

  rf <- arrow::ReadableFile$create(db_file)
  fr <- arrow::FeatherReader$create(rf)
  on.exit(try(rf$close(), silent = TRUE), add = TRUE)
  col_names <- names(fr)

  if ("features" %in% col_names) {
    return("features")
  }
  if ("motifs" %in% col_names) {
    msg("[WARN] Database '%s' does not expose 'features'; falling back to dbIndexCol='motifs'.\n", basename(db_file))
    return("motifs")
  }

  stop(sprintf(
    "Could not infer SCENIC database index column for '%s'. Please set SCENIC_DB_INDEX_COL explicitly.",
    basename(db_file)
  ))
}

build_input_fingerprint <- function() {
  if (file.exists(INPUT_RDS)) {
    return(file_signature(INPUT_RDS))
  }
  list(list(path = normalize_path_safe(INPUT_RDS), size = NA_real_, mtime = NA_character_))
}

build_inference_manifest <- function(expr_mat_dense, scenic_options) {
  list(
    stage = "inference",
    input = build_input_fingerprint(),
    assay = ASSAY_USE,
    organism = ORGANISM,
    cell_type_col = CELL_TYPE_COL,
    sample_col = SAMPLE_COL,
    db_files = file_signature(unname(scenic_options@settings$dbFilesSelected)),
    db_index_col = scenic_options@settings$dbIndexCol %||% SCENIC_DB_INDEX_COL %||% "features",
    counts = matrix_signature(expr_mat_dense),
    params = list(
      min_genes_per_cell = MIN_GENES_PER_CELL,
      min_cells_per_gene = MIN_CELLS_PER_GENE,
      min_gene_pct = MIN_GENE_PCT,
      inference_max_cells_per_sample_celltype = INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE,
      inference_max_cells_per_celltype = INFERENCE_MAX_CELLS_PER_CELLTYPE,
      inference_global_max_cells = INFERENCE_GLOBAL_MAX_CELLS
    )
  )
}

build_regulon_manifest <- function(inference_manifest, min_genes_per_regulon) {
  utils::modifyList(
    inference_manifest,
    list(
      stage = "regulon",
      params = utils::modifyList(
        inference_manifest$params,
        list(min_genes_per_regulon = min_genes_per_regulon)
      )
    )
  )
}

build_scoring_manifest <- function(full_counts_sparse, regulons) {
  list(
    stage = "scoring",
    input = build_input_fingerprint(),
    assay = ASSAY_USE,
    counts = matrix_signature(full_counts_sparse),
    regulons = regulon_signature(regulons)
  )
}

build_rss_manifest <- function(auc_mat, cell_types) {
  list(
    stage = "rss",
    auc = matrix_signature(auc_mat),
    cell_types = vector_signature(cell_types)
  )
}

get_scenic_int_dir <- function(scenic_options) {
  file.path(scenic_options@settings$outDir, "int")
}

run_scenic_coexpression_modules <- function(scenic_options) {
  if (exists("runSCENIC_1_coexpressionModules", where = asNamespace("SCENIC"), inherits = FALSE)) {
    runSCENIC_1_coexpressionModules(scenic_options)
  } else if (exists("runSCENIC_1_coexNetwork2modules", where = asNamespace("SCENIC"), inherits = FALSE)) {
    runSCENIC_1_coexNetwork2modules(scenic_options)
  } else {
    stop("Neither runSCENIC_1_coexpressionModules nor runSCENIC_1_coexNetwork2modules is available in the installed SCENIC package.")
  }
}

stop_if_missing <- function(x, name) {
  if (is.null(x)) stop(sprintf("'%s' is NULL.", name))
}

get_counts_matrix <- function(obj, assay = "RNA") {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(assay %in% names(obj@assays))

  assay_obj <- obj[[assay]]

  if (inherits(assay_obj, "Assay5")) {
    lyr <- Layers(assay_obj)
    if (!"counts" %in% lyr) {
      stop(sprintf("Assay5 '%s' has no counts layer.", assay))
    }
    mat <- LayerData(obj, assay = assay, layer = "counts")
  } else {
    mat <- GetAssayData(obj, assay = assay, slot = "counts")
  }

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

choose_reduction <- function(obj, candidates = REDUCTION_CANDIDATES) {
  reds <- Reductions(obj)
  hit  <- candidates[candidates %in% reds]
  if (length(hit) > 0) return(hit[1])
  if (length(reds) > 0) return(reds[1])
  NULL
}

extract_tf_name <- function(regulon_name) {
  sub("\\s*\\(.*$", "", regulon_name)
}

extract_regulon_gene_vector <- function(x) {
  if (is.character(x)) return(unique(x))
  if (methods::is(x, "GeneSet")) {
    return(unique(GSEABase::geneIds(x)))
  }
  if (is.list(x)) {
    return(unique(as.character(unlist(x, use.names = FALSE))))
  }
  unique(as.character(x))
}

make_parallel_param <- function(n_cores) {
  if (.Platform$OS.type == "windows") {
    SnowParam(workers = n_cores, type = "SOCK", progressbar = TRUE)
  } else {
    MulticoreParam(workers = n_cores, progressbar = TRUE)
  }
}

estimate_dense_matrix_gb <- function(nrow_mat, ncol_mat, bytes_per_entry = 8) {
  as.numeric(nrow_mat) * as.numeric(ncol_mat) * bytes_per_entry / (1024^3)
}

assert_dense_conversion_safe <- function(mat_sparse, max_dense_gb = MAX_INFERENCE_DENSE_GB) {
  est_gb <- estimate_dense_matrix_gb(nrow(mat_sparse), ncol(mat_sparse))
  msg(
    "[INFO] Dense inference matrix estimate: %.2f GB for %d genes x %d cells (limit=%.2f GB)\n",
    est_gb,
    nrow(mat_sparse),
    ncol(mat_sparse),
    max_dense_gb
  )
  if (is.finite(max_dense_gb) && est_gb > max_dense_gb) {
    stop(sprintf(
      paste(
        "Refusing to densify inference matrix: estimated dense size %.2f GB exceeds MAX_INFERENCE_DENSE_GB=%.2f GB.",
        "Reduce inference cells/genes or raise MAX_INFERENCE_DENSE_GB if you really want to spend the RAM."
      ),
      est_gb,
      max_dense_gb
    ))
  }
  invisible(est_gb)
}

stratified_downsample_cells <- function(
  obj,
  cell_type_col,
  sample_col,
  max_cells_per_sample_celltype = 200,
  max_cells_per_celltype = 3000,
  global_max_cells = 30000,
  seed = RANDOM_SEED
) {
  meta <- obj@meta.data
  stopifnot(cell_type_col %in% colnames(meta))
  stopifnot(sample_col %in% colnames(meta))

  if (!is.null(seed) && length(seed) == 1 && !is.na(seed)) {
    set.seed(seed)
  }

  ct <- as.character(meta[[cell_type_col]])
  sp <- as.character(meta[[sample_col]])

  keep <- !is.na(ct) & nzchar(trimws(ct)) & !is.na(sp) & nzchar(trimws(sp))
  meta <- meta[keep, , drop = FALSE]

  meta$cell_id    <- rownames(meta)
  meta$cell_type  <- as.character(meta[[cell_type_col]])
  meta$sample_id  <- as.character(meta[[sample_col]])
  meta$stratum_id <- paste(meta$cell_type, meta$sample_id, sep = "___")

  picked <- unlist(lapply(split(meta$cell_id, meta$stratum_id), function(v) {
    sample(v, min(length(v), max_cells_per_sample_celltype))
  }), use.names = FALSE)

  meta2 <- meta[picked, , drop = FALSE]
  picked2 <- unlist(lapply(split(meta2$cell_id, meta2$cell_type), function(cells) {
    sample(cells, min(length(cells), max_cells_per_celltype))
  }), use.names = FALSE)

  picked2 <- unique(picked2)

  if (length(picked2) > global_max_cells) {
    picked2 <- sample(picked2, global_max_cells)
  }

  if (length(picked2) < 2) {
    stop(sprintf(
      "Too few cells selected for SCENIC inference after downsampling: %d. Check cell_type/sample columns or downsampling thresholds.",
      length(picked2)
    ))
  }

  picked2
}

aggregate_matrix_by_group <- function(mat, groups) {
  groups <- as.character(groups)
  keep   <- !is.na(groups) & nzchar(trimws(groups))

  mat    <- mat[, keep, drop = FALSE]
  groups <- groups[keep]

  levs <- unique(groups)
  out  <- vapply(levs, function(g) {
    rowMeans(mat[, groups == g, drop = FALSE])
  }, numeric(nrow(mat)))

  if (is.vector(out)) {
    out <- matrix(out, nrow = nrow(mat), ncol = 1)
    colnames(out) <- levs[1]
  }

  rownames(out) <- rownames(mat)
  colnames(out) <- levs
  out
}

pick_top_variable_regulons <- function(auc_mat, n_top = 12) {
  vars <- apply(auc_mat, 1, var, na.rm = TRUE)
  vars <- sort(vars, decreasing = TRUE)
  names(vars)[seq_len(min(n_top, length(vars)))]
}

compute_rss_safe <- function(auc_mat, cell_types) {
  if (!"calcRSS" %in% getNamespaceExports("SCENIC")) {
    msg("[WARN] calcRSS not exported by SCENIC; skip RSS.\n")
    return(NULL)
  }

  keep <- !is.na(cell_types) & nzchar(trimws(as.character(cell_types)))
  auc_use <- auc_mat[, keep, drop = FALSE]
  ct_use  <- as.character(cell_types[keep])

  tryCatch({
    rss <- calcRSS(AUC = auc_use, cellAnnotation = ct_use)
    as.matrix(rss)
  }, error = function(e) {
    msg("[WARN] RSS failed: %s\n", e$message)
    NULL
  })
}

pick_top_rss_regulons <- function(rss_mat, top_n_per_type = 8) {
  if (is.null(rss_mat) || nrow(rss_mat) == 0 || ncol(rss_mat) == 0) {
    return(character(0))
  }

  regs <- unique(unlist(lapply(seq_len(ncol(rss_mat)), function(i) {
    ord <- order(rss_mat[, i], decreasing = TRUE)
    rownames(rss_mat)[ord][seq_len(min(top_n_per_type, nrow(rss_mat)))]
  }), use.names = FALSE))

  regs
}

save_matrix_csv <- function(mat, file) {
  fwrite(
    data.table(feature = rownames(mat), as.data.frame(mat), check.names = FALSE),
    file
  )
}

validate_scenic_runtime_inputs <- function(
  seurat_obj,
  assay_use,
  cell_type_col,
  sample_col,
  database_dir,
  output_dir = OUTPUT_DIR
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("run_scenic_module() requires a Seurat object.")
  }
  if (!assay_use %in% names(seurat_obj@assays)) {
    stop(sprintf("ASSAY_USE '%s' not found in Seurat object.", assay_use))
  }
  missing_cols <- setdiff(c(cell_type_col, sample_col), colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required metadata columns: %s", paste(missing_cols, collapse = ", ")))
  }
  if (!dir.exists(database_dir)) {
    stop(sprintf("DATABASE_DIR not found: %s", database_dir))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  invisible(TRUE)
}

build_scenic_feature_map <- function(original_features, assay_features) {
  if (length(original_features) != length(assay_features)) {
    stop("Feature map construction failed: original and assay feature counts differ.")
  }
  data.frame(
    original_regulon = as.character(original_features),
    assay_feature = as.character(assay_features),
    stringsAsFactors = FALSE
  )
}

map_regulons_to_assay_features <- function(seurat_obj, regulon_names, assay = "SCENIC") {
  regulon_names <- as.character(regulon_names)
  if (length(regulon_names) == 0 || !assay %in% names(seurat_obj@assays)) {
    return(regulon_names)
  }

  assay_features <- rownames(seurat_obj[[assay]])
  mapped <- regulon_names

  scenic_map <- tryCatch(seurat_obj@misc$scenic_feature_map, error = function(e) NULL)
  if (is.data.frame(scenic_map) && all(c("original_regulon", "assay_feature") %in% colnames(scenic_map))) {
    lookup <- setNames(as.character(scenic_map$assay_feature), as.character(scenic_map$original_regulon))
    mapped_lookup <- unname(lookup[regulon_names])
    valid_lookup <- !is.na(mapped_lookup) & nzchar(mapped_lookup)
    mapped[valid_lookup] <- mapped_lookup[valid_lookup]
  }

  missing <- !(mapped %in% assay_features)
  if (any(missing)) {
    fallback <- gsub("_", "-", regulon_names[missing], fixed = TRUE)
    fallback_ok <- fallback %in% assay_features
    mapped[which(missing)[fallback_ok]] <- fallback[fallback_ok]
  }

  mapped
}

# ==============================================================================
# 4. SCENIC Core
# ==============================================================================

check_scenic_requirements <- function(
  db_dir,
  organism = ORGANISM,
  db_10kb = SCENIC_DB_10KB,
  db_500bp = SCENIC_DB_500BP
) {
  db_files <- resolve_scenic_database_files(
    organism = organism,
    db_dir = db_dir,
    db_10kb = db_10kb,
    db_500bp = db_500bp
  )
  invisible(db_files)
}

ensure_rcistarget_annotations <- function(organism) {
  obj_name <- sprintf("motifAnnotations_%s", organism)
  if (exists(obj_name, inherits = TRUE)) {
    return(invisible(TRUE))
  }

  pkg_dir <- find.package("RcisTarget")
  data_dir <- file.path(pkg_dir, "data")
  primary_rdata <- file.path(data_dir, sprintf("%s.RData", obj_name))
  versioned_obj <- sprintf("%s_v9", obj_name)
  versioned_rdata <- file.path(data_dir, sprintf("%s.RData", versioned_obj))

  target_env <- .GlobalEnv

  if (file.exists(primary_rdata)) {
    load(primary_rdata, envir = target_env)
  }
  if (!exists(obj_name, envir = target_env, inherits = FALSE) && file.exists(versioned_rdata)) {
    load(versioned_rdata, envir = target_env)
    if (exists(versioned_obj, envir = target_env, inherits = FALSE)) {
      assign(obj_name, get(versioned_obj, envir = target_env), envir = target_env)
    }
  }

  if (!exists(obj_name, envir = target_env, inherits = FALSE)) {
    stop(sprintf("Failed to load RcisTarget motif annotation object: %s", obj_name))
  }

  invisible(TRUE)
}

initialize_scenic_options_local <- function(
  organism,
  db_dir,
  out_dir,
  dataset_title,
  db_10kb = SCENIC_DB_10KB,
  db_500bp = SCENIC_DB_500BP,
  db_index_col = SCENIC_DB_INDEX_COL
) {
  db_files <- resolve_scenic_database_files(
    organism = organism,
    db_dir = db_dir,
    db_10kb = db_10kb,
    db_500bp = db_500bp
  )
  db_10kb  <- db_files[["10kb"]]
  db_500bp <- db_files[["500bp"]]
  resolved_db_index_col <- resolve_db_index_col(db_files, preferred = db_index_col)

  ensure_rcistarget_annotations(organism)

  scenic_options <- initializeScenic(
    org = organism,
    dbDir = db_dir,
    dbs = c("10kb" = basename(db_10kb), "500bp" = basename(db_500bp)),
    datasetTitle = dataset_title
  )

  scenic_options@settings$outDir   <- out_dir
  scenic_options@settings$verbose  <- TRUE
  scenic_options@settings$nCores   <- N_CORES
  scenic_options@settings$dbFilesSelected <- unname(db_files)
  scenic_options@settings$dbIndexCol <- resolved_db_index_col

  scenic_options
}

prepare_scenic_input <- function(
  seurat_obj,
  assay_use,
  cell_type_col,
  sample_col,
  min_genes_per_cell,
  min_cells_per_gene,
  min_gene_pct
) {
  stopifnot(cell_type_col %in% colnames(seurat_obj@meta.data))
  stopifnot(sample_col %in% colnames(seurat_obj@meta.data))

  counts_all <- get_counts_matrix(seurat_obj, assay = assay_use)
  meta_all   <- seurat_obj@meta.data

  ct <- as.character(meta_all[[cell_type_col]])
  sp <- as.character(meta_all[[sample_col]])

  keep_cells <- !is.na(ct) & nzchar(trimws(ct)) &
    !is.na(sp) & nzchar(trimws(sp))

  counts_metadata_valid <- counts_all[, keep_cells, drop = FALSE]
  meta_metadata_valid   <- meta_all[keep_cells, , drop = FALSE]

  msg("[INFO] After metadata filter: cells=%d genes=%d\n", ncol(counts_metadata_valid), nrow(counts_metadata_valid))

  keep_cells2 <- Matrix::colSums(counts_metadata_valid > 0) >= min_genes_per_cell
  counts_metadata_valid <- counts_metadata_valid[, keep_cells2, drop = FALSE]
  meta_metadata_valid   <- meta_metadata_valid[keep_cells2, , drop = FALSE]
  meta_inference <- meta_metadata_valid

  min_cells_by_pct <- ceiling(ncol(counts_metadata_valid) * min_gene_pct)
  min_cells_final  <- max(min_cells_per_gene, min_cells_by_pct)

  keep_genes <- Matrix::rowSums(counts_metadata_valid > 0) >= min_cells_final
  counts_inference <- counts_metadata_valid[keep_genes, , drop = FALSE]

  msg("[INFO] After expression filter: cells=%d genes=%d\n", ncol(counts_inference), nrow(counts_inference))
  msg("[INFO] Gene detection threshold: %d cells\n", min_cells_final)

  list(
    counts_scoring = counts_all,
    meta_all = meta_all,
    counts_metadata_valid = counts_metadata_valid,
    meta_metadata_valid = meta_metadata_valid,
    counts_inference = counts_inference,
    meta_inference = meta_inference
  )
}

run_scenic_regulon_inference <- function(
  expr_mat_dense,
  scenic_options,
  n_cores,
  min_genes_per_regulon
) {
  bp <- make_parallel_param(n_cores)
  register(bp, default = TRUE)

  int_dir <- get_scenic_int_dir(scenic_options)
  dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)

  genes_kept_file <- file.path(int_dir, "1.1_genesKept.Rds")
  genie3_link_file <- file.path(int_dir, "1.4_GENIE3_linkList.Rds")
  inference_manifest_file <- file.path(int_dir, "manifest_inference.rds")
  regulon_manifest_file <- file.path(int_dir, "manifest_regulon.rds")
  db_index_col <- scenic_options@settings$dbIndexCol %||% SCENIC_DB_INDEX_COL %||% "features"
  regulon_candidates <- c(
    file.path(int_dir, "3.4_regulons_forAUCell.Rds"),
    file.path(int_dir, "2.6_regulons_asGeneSet.Rds")
  )
  existing_regulon_file <- regulon_candidates[file.exists(regulon_candidates)][1]
  inference_manifest <- build_inference_manifest(expr_mat_dense, scenic_options)
  regulon_manifest <- build_regulon_manifest(inference_manifest, min_genes_per_regulon)
  can_reuse_inference <- manifest_matches(inference_manifest_file, inference_manifest) &&
    file.exists(genes_kept_file) && file.exists(genie3_link_file)
  can_reuse_regulons <- manifest_matches(regulon_manifest_file, regulon_manifest) &&
    !is.na(existing_regulon_file) && nzchar(existing_regulon_file)

  if (file.exists(genes_kept_file) && file.exists(genie3_link_file) && !can_reuse_inference) {
    msg("[SCENIC] Existing geneFiltering/GENIE3 artifacts found but manifest mismatch; recomputing.\n")
  }

  if (can_reuse_inference) {
    msg("[SCENIC] Reusing existing geneFiltering + GENIE3 results from %s\n", int_dir)
    genes_kept <- readRDS(genes_kept_file)
    expr_filt <- expr_mat_dense[genes_kept, , drop = FALSE]
  } else {
    msg("[SCENIC] Step 1/5 geneFiltering...\n")
    genes_kept <- geneFiltering(expr_mat_dense, scenic_options)
    expr_filt  <- expr_mat_dense[genes_kept, , drop = FALSE]
    saveRDS(genes_kept, file.path(DIR_RDS, "genes_kept_after_geneFiltering.rds"))

    msg("[SCENIC] Step 2/5 runCorrelation...\n")
    runCorrelation(expr_filt, scenic_options)

    msg("[SCENIC] Step 3/5 runGenie3...\n")
    runGenie3(expr_filt, scenic_options, nParts = n_cores)
    write_manifest(inference_manifest_file, inference_manifest)
  }

  if (!is.na(existing_regulon_file) && nzchar(existing_regulon_file) && !can_reuse_regulons) {
    msg("[SCENIC] Existing regulon artifact found but manifest mismatch; rebuilding regulons.\n")
  }

  if (can_reuse_regulons) {
    msg("[SCENIC] Reusing existing regulon file: %s\n", existing_regulon_file)
    regulons <- readRDS(existing_regulon_file)
  } else {
    msg("[SCENIC] Step 4/5 build coexpression modules...\n")
    run_scenic_coexpression_modules(scenic_options)

    msg("[SCENIC] Step 5/5 runSCENIC_2_createRegulons...\n")
    runSCENIC_2_createRegulons(
      scenic_options,
      minGenes = min_genes_per_regulon,
      dbIndexCol = db_index_col
    )

    existing_regulon_file <- regulon_candidates[file.exists(regulon_candidates)][1]
    if (is.na(existing_regulon_file) || !nzchar(existing_regulon_file)) {
      stop(sprintf(
        "SCENIC regulon file not found after runSCENIC_2_createRegulons(). Looked for: %s",
        paste(basename(regulon_candidates), collapse = ", ")
      ))
    }
    regulons <- readRDS(existing_regulon_file)
    write_manifest(regulon_manifest_file, regulon_manifest)
  }

  list(
    scenic_options = scenic_options,
    genes_kept = genes_kept,
    regulons = regulons
  )
}

score_regulons_full_object <- function(full_counts_sparse, regulons, n_cores) {
  rankings_file <- file.path(DIR_RDS, "aucell_rankings_full_object.rds")
  regulon_auc_file <- file.path(DIR_RDS, "regulon_auc_full_object.rds")
  auc_mat_file <- file.path(DIR_RDS, "regulon_auc_matrix_full_object.rds")
  scoring_manifest_file <- file.path(DIR_RDS, "manifest_scoring.rds")
  scoring_manifest <- build_scoring_manifest(full_counts_sparse, regulons)

  if (
    manifest_matches(scoring_manifest_file, scoring_manifest) &&
      file.exists(rankings_file) && file.exists(regulon_auc_file) && file.exists(auc_mat_file)
  ) {
    msg("[AUCell] Reusing existing rankings/AUC results from %s\n", DIR_RDS)
    return(list(
      rankings = readRDS(rankings_file),
      regulon_auc = readRDS(regulon_auc_file),
      auc_mat = readRDS(auc_mat_file)
    ))
  }

  if (file.exists(rankings_file) || file.exists(regulon_auc_file) || file.exists(auc_mat_file)) {
    msg("[AUCell] Existing scoring artifacts found but manifest mismatch; recalculating.\n")
  }

  msg("[AUCell] Building rankings on full object...\n")

  rankings <- AUCell_buildRankings(
    full_counts_sparse,
    nCores = n_cores,
    plotStats = FALSE,
    splitByBlocks = TRUE,
    verbose = TRUE
  )

  saveRDS(rankings, rankings_file)

  regulon_gene_sets <- lapply(regulons, extract_regulon_gene_vector)
  regulon_gene_sets <- regulon_gene_sets[lengths(regulon_gene_sets) > 0]

  msg("[AUCell] Calculating AUC...\n")
  regulon_auc <- AUCell_calcAUC(
    geneSets = regulon_gene_sets,
    rankings = rankings,
    nCores = n_cores,
    normAUC = TRUE
  )

  auc_mat <- getAUC(regulon_auc)
  saveRDS(regulon_auc, regulon_auc_file)
  saveRDS(auc_mat, auc_mat_file)
  write_manifest(scoring_manifest_file, scoring_manifest)

  list(
    rankings = rankings,
    regulon_auc = regulon_auc,
    auc_mat = auc_mat
  )
}

integrate_scenic_to_seurat <- function(seurat_obj, auc_mat) {
  common_cells <- intersect(colnames(seurat_obj), colnames(auc_mat))
  if (length(common_cells) == 0) {
    stop("No common cells between Seurat object and AUC matrix.")
  }

  auc_use <- auc_mat[, common_cells, drop = FALSE]
  original_regulons <- rownames(auc_use)
  scenic_assay <- withCallingHandlers(
    CreateAssayObject(data = auc_use),
    warning = function(w) {
      if (grepl("Feature names cannot have underscores", conditionMessage(w), fixed = TRUE)) {
        msg("[INFO] Seurat normalized SCENIC regulon feature names from '_' to '-'; feature map recorded in obj@misc$scenic_feature_map.\n")
        invokeRestart("muffleWarning")
      }
    }
  )
  scenic_feature_map <- build_scenic_feature_map(original_regulons, rownames(scenic_assay))
  assay_feature_lookup <- setNames(
    scenic_feature_map$assay_feature,
    scenic_feature_map$original_regulon
  )
  seurat_obj[["SCENIC"]] <- scenic_assay
  seurat_obj@misc$scenic_feature_map <- scenic_feature_map

  cell_ids <- colnames(auc_use)
  top_reg_idx <- apply(auc_use, 2, which.max)
  top_reg_name <- rownames(auc_use)[top_reg_idx]
  names(top_reg_name) <- cell_ids
  top_reg_name_assay <- stats::setNames(
    unname(assay_feature_lookup[top_reg_name]),
    cell_ids
  )
  top_reg_score <- apply(auc_use, 2, max)
  names(top_reg_score) <- cell_ids

  seurat_obj$scenic_top_regulon <- NA_character_
  seurat_obj$scenic_top_regulon[names(top_reg_name)] <- top_reg_name

  seurat_obj$scenic_top_regulon_assay_feature <- NA_character_
  seurat_obj$scenic_top_regulon_assay_feature[names(top_reg_name_assay)] <- top_reg_name_assay

  seurat_obj$scenic_top_regulon_auc <- NA_real_
  seurat_obj$scenic_top_regulon_auc[names(top_reg_score)] <- top_reg_score

  seurat_obj
}

# ==============================================================================
# 5. Post-processing
# ==============================================================================

plot_regulon_activity_heatmap <- function(auc_mat, cell_types, rss_mat = NULL) {
  avg_auc <- aggregate_matrix_by_group(auc_mat, cell_types)
  save_matrix_csv(avg_auc, file.path(DIR_TABLE, "regulon_mean_auc_by_celltype.csv"))

  regs_rss <- pick_top_rss_regulons(rss_mat, top_n_per_type = TOP_RSS_PER_CELLTYPE)
  regs_var <- pick_top_variable_regulons(avg_auc, n_top = TOP_VAR_REGULONS_FOR_PLOTS)
  regs_use <- unique(c(regs_rss, regs_var))
  regs_use <- regs_use[regs_use %in% rownames(avg_auc)]

  if (length(regs_use) < 2) {
    regs_use <- rownames(avg_auc)[seq_len(min(20, nrow(avg_auc)))]
  }

  hm <- avg_auc[regs_use, , drop = FALSE]
  hm <- hm[order(rownames(hm)), , drop = FALSE]

  pheatmap(
    hm,
    scale = "row",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 10,
    main = "Regulon activity by cell type",
    filename = file.path(DIR_FIG, "regulon_activity_heatmap.pdf"),
    width = HEATMAP_WIDTH,
    height = HEATMAP_HEIGHT
  )

  fwrite(
    data.table(regulon = regs_use),
    file.path(DIR_TABLE, "regulons_used_in_heatmap.csv")
  )
}

plot_rss_heatmap <- function(rss_mat) {
  if (is.null(rss_mat)) return(invisible(NULL))

  save_matrix_csv(rss_mat, file.path(DIR_TABLE, "regulon_rss_by_celltype.csv"))

  regs <- pick_top_rss_regulons(rss_mat, top_n_per_type = TOP_RSS_PER_CELLTYPE)
  regs <- regs[regs %in% rownames(rss_mat)]
  if (length(regs) == 0) return(invisible(NULL))

  hm <- rss_mat[regs, , drop = FALSE]

  pheatmap(
    hm,
    scale = "none",
    color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 10,
    main = "Regulon Specificity Score (RSS)",
    filename = file.path(DIR_FIG, "regulon_rss_heatmap.pdf"),
    width = HEATMAP_WIDTH,
    height = HEATMAP_HEIGHT
  )
}

plot_top_regulon_featureplots <- function(seurat_obj, auc_mat, rss_mat = NULL) {
  red <- choose_reduction(seurat_obj)
  if (is.null(red)) {
    msg("[WARN] No reduction found; skip FeaturePlot.\n")
    return(invisible(NULL))
  }

  regs_rss <- pick_top_rss_regulons(rss_mat, top_n_per_type = 2)
  regs_var <- pick_top_variable_regulons(auc_mat, n_top = TOP_VAR_REGULONS_FOR_PLOTS)
  regs_use <- unique(c(regs_rss, regs_var))
  regs_use <- regs_use[regs_use %in% rownames(auc_mat)]
  regs_use <- regs_use[seq_len(min(length(regs_use), TOP_VAR_REGULONS_FOR_PLOTS))]

  if (length(regs_use) == 0) return(invisible(NULL))

  old_assay <- DefaultAssay(seurat_obj)
  DefaultAssay(seurat_obj) <- "SCENIC"
  on.exit(DefaultAssay(seurat_obj) <- old_assay, add = TRUE)

  assay_features <- rownames(seurat_obj[["SCENIC"]])
  plot_features <- map_regulons_to_assay_features(seurat_obj, regs_use, assay = "SCENIC")
  keep <- !is.na(plot_features) & nzchar(plot_features) & plot_features %in% assay_features
  regs_use <- regs_use[keep]
  plot_features <- plot_features[keep]

  if (length(plot_features) == 0) {
    msg("[WARN] No regulon features remained after Seurat feature-name normalization; skip FeaturePlot.\n")
    return(invisible(NULL))
  }

  pdf(file.path(DIR_FIG, "top_regulon_featureplots.pdf"),
      width = FEATUREPLOT_PDF_WIDTH,
      height = FEATUREPLOT_PDF_HEIGHT)

  for (i in seq_along(regs_use)) {
    reg <- regs_use[i]
    plot_reg <- plot_features[i]
    p <- FeaturePlot(
      seurat_obj,
      features = plot_reg,
      reduction = red,
      raster = TRUE
    ) + ggtitle(reg)
    print(p)
  }

  dev.off()

  fwrite(
    data.table(regulon = regs_use, plotted_feature = plot_features),
    file.path(DIR_TABLE, "top_regulons_featureplot.csv")
  )
}

save_regulon_gene_tables <- function(regulons) {
  dt_list <- lapply(names(regulons), function(reg) {
    genes <- extract_regulon_gene_vector(regulons[[reg]])
    data.table(
      regulon = reg,
      tf = extract_tf_name(reg),
      target_gene = genes
    )
  })

  edges_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  fwrite(edges_dt, file.path(DIR_TABLE, "regulon_targets_full.csv"))
}

plot_regulon_network <- function(regulons, rss_mat = NULL) {
  regs_rss <- pick_top_rss_regulons(rss_mat, top_n_per_type = 3)
  regs_all <- names(regulons)

  if (length(regs_rss) == 0) {
    regs_use <- regs_all[seq_len(min(NETWORK_TOP_REGULONS, length(regs_all)))]
  } else {
    regs_use <- unique(regs_rss)
    if (length(regs_use) > NETWORK_TOP_REGULONS) {
      regs_use <- regs_use[seq_len(NETWORK_TOP_REGULONS)]
    }
  }

  edge_list <- rbindlist(lapply(regs_use, function(reg) {
    genes <- extract_regulon_gene_vector(regulons[[reg]])
    genes <- genes[seq_len(min(length(genes), NETWORK_TOP_TARGETS_PER_REGULON))]
    data.table(
      from = extract_tf_name(reg),
      to   = genes,
      regulon = reg
    )
  }), use.names = TRUE, fill = TRUE)

  if (nrow(edge_list) == 0) {
    msg("[WARN] No regulon-target edges available; skip regulon network export.\n")
    return(invisible(NULL))
  }

  fwrite(edge_list, file.path(DIR_TABLE, "regulon_network_edges.csv"))

  node_names <- unique(c(edge_list$from, edge_list$to))
  nodes <- data.table(
    id = node_names,
    label = node_names,
    group = ifelse(node_names %in% unique(edge_list$from), "TF", "Target")
  )

  vis <- visNetwork(
    nodes = nodes,
    edges = edge_list[, .(from, to)]
  ) |>
    visGroups(groupname = "TF", color = list(background = "#D55E00")) |>
    visGroups(groupname = "Target", color = list(background = "#0072B2")) |>
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
    visPhysics(stabilization = TRUE)

  htmlwidgets::saveWidget(
    vis,
    file = file.path(DIR_FIG, "regulon_network.html"),
    selfcontained = FALSE
  )
}

# ==============================================================================
# 6. Main Runner
# ==============================================================================

run_scenic_module <- function(seurat_obj) {
  msg("============================================================\n")
  msg("SCENIC module started\n")
  msg("============================================================\n")

  initialize_output_dirs(OUTPUT_DIR)
  validate_scenic_runtime_inputs(
    seurat_obj = seurat_obj,
    assay_use = ASSAY_USE,
    cell_type_col = CELL_TYPE_COL,
    sample_col = SAMPLE_COL,
    database_dir = DATABASE_DIR,
    output_dir = OUTPUT_DIR
  )
  old_wd <- getwd()
  setwd(OUTPUT_DIR)
  on.exit(setwd(old_wd), add = TRUE)

  # ----- checks -----
  if (!CELL_TYPE_COL %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("CELL_TYPE_COL '%s' not found.", CELL_TYPE_COL))
  }
  if (!SAMPLE_COL %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("SAMPLE_COL '%s' not found.", SAMPLE_COL))
  }

  selected_db_files <- check_scenic_requirements(
    db_dir = DATABASE_DIR,
    organism = ORGANISM,
    db_10kb = SCENIC_DB_10KB,
    db_500bp = SCENIC_DB_500BP
  )
  msg("[INFO] Selected SCENIC DBs: %s | %s\n", basename(selected_db_files[["10kb"]]), basename(selected_db_files[["500bp"]]))

  # ----- prepare full input -----
  prep <- prepare_scenic_input(
    seurat_obj = seurat_obj,
    assay_use = ASSAY_USE,
    cell_type_col = CELL_TYPE_COL,
    sample_col = SAMPLE_COL,
    min_genes_per_cell = MIN_GENES_PER_CELL,
    min_cells_per_gene = MIN_CELLS_PER_GENE,
    min_gene_pct = MIN_GENE_PCT
  )

  counts_scoring <- prep$counts_scoring
  counts_inference <- prep$counts_inference
  meta_inference <- prep$meta_inference

  msg("[INFO] Full scoring cells: %d\n", ncol(counts_scoring))
  msg("[INFO] Metadata/QC eligible inference cells: %d\n", ncol(counts_inference))

  fwrite(
    data.table(
      cell = rownames(meta_inference),
      cell_type = as.character(meta_inference[[CELL_TYPE_COL]]),
      sample = as.character(meta_inference[[SAMPLE_COL]])
    ),
    file.path(DIR_TABLE, "cell_metadata_used_for_scenic.csv")
  )

  # ----- inference subset -----
  obj_use <- subset(seurat_obj, cells = colnames(counts_inference))

  infer_cells <- stratified_downsample_cells(
    obj = obj_use,
    cell_type_col = CELL_TYPE_COL,
    sample_col = SAMPLE_COL,
    max_cells_per_sample_celltype = INFERENCE_MAX_CELLS_PER_SAMPLE_CELLTYPE,
    max_cells_per_celltype = INFERENCE_MAX_CELLS_PER_CELLTYPE,
    global_max_cells = INFERENCE_GLOBAL_MAX_CELLS,
    seed = RANDOM_SEED
  )

  msg("[INFO] Inference subset cells: %d\n", length(infer_cells))

  counts_infer_sparse <- counts_inference[, infer_cells, drop = FALSE]
  assert_dense_conversion_safe(counts_infer_sparse, max_dense_gb = MAX_INFERENCE_DENSE_GB)
  counts_infer_dense  <- as.matrix(counts_infer_sparse)

  saveRDS(infer_cells, file.path(DIR_RDS, "scenic_inference_cells.rds"))

  # ----- scenic options -----
  scenic_options <- initialize_scenic_options_local(
    organism = ORGANISM,
    db_dir = DATABASE_DIR,
    out_dir = OUTPUT_DIR,
    dataset_title = SCENIC_DATASET_TITLE,
    db_10kb = selected_db_files[["10kb"]],
    db_500bp = selected_db_files[["500bp"]],
    db_index_col = SCENIC_DB_INDEX_COL
  )

  scenic_options@inputDatasetInfo$cellInfo <- data.frame(
    cell_id   = colnames(counts_infer_dense),
    cell_type = as.character(obj_use@meta.data[colnames(counts_infer_dense), CELL_TYPE_COL]),
    sample    = as.character(obj_use@meta.data[colnames(counts_infer_dense), SAMPLE_COL]),
    row.names = colnames(counts_infer_dense),
    stringsAsFactors = FALSE
  )

  saveRDS(scenic_options, file.path(DIR_RDS, "scenic_options_initialized.rds"))

  # ----- network inference -----
  scenic_core <- run_scenic_regulon_inference(
    expr_mat_dense = counts_infer_dense,
    scenic_options = scenic_options,
    n_cores = N_CORES,
    min_genes_per_regulon = MIN_GENES_PER_REGULON
  )

  regulons <- scenic_core$regulons
  saveRDS(regulons, file.path(DIR_RDS, "regulons_final.rds"))
  save_regulon_gene_tables(regulons)

  # ----- full object scoring -----
  scoring <- score_regulons_full_object(
    full_counts_sparse = counts_scoring,
    regulons = regulons,
    n_cores = N_CORES
  )

  auc_mat <- scoring$auc_mat

  # ----- RSS -----
  cell_types_full <- as.character(seurat_obj@meta.data[colnames(auc_mat), CELL_TYPE_COL])
  rss_keep <- !is.na(cell_types_full) & nzchar(trimws(cell_types_full))
  msg("[INFO] RSS eligible cells: %d / %d\n", sum(rss_keep), length(cell_types_full))
  rss_file <- file.path(DIR_RDS, "regulon_rss.rds")
  rss_manifest_file <- file.path(DIR_RDS, "manifest_rss.rds")
  rss_manifest <- build_rss_manifest(auc_mat, cell_types_full)
  if (manifest_matches(rss_manifest_file, rss_manifest) && file.exists(rss_file)) {
    msg("[SCENIC] Reusing existing RSS matrix from %s\n", rss_file)
    rss_mat <- readRDS(rss_file)
  } else {
    if (file.exists(rss_file)) {
      msg("[SCENIC] Existing RSS artifact found but manifest mismatch; recalculating.\n")
    }
    rss_mat <- compute_rss_safe(auc_mat, cell_types_full)
    if (!is.null(rss_mat)) {
      saveRDS(rss_mat, rss_file)
      write_manifest(rss_manifest_file, rss_manifest)
    }
  }

  # ----- integrate back -----
  seurat_obj_out <- seurat_obj
  seurat_obj_out <- integrate_scenic_to_seurat(seurat_obj_out, auc_mat)

  # ----- plots -----
  plot_regulon_activity_heatmap(
    auc_mat = auc_mat,
    cell_types = cell_types_full,
    rss_mat = rss_mat
  )
  plot_rss_heatmap(rss_mat)
  plot_top_regulon_featureplots(
    seurat_obj = seurat_obj_out,
    auc_mat = auc_mat,
    rss_mat = rss_mat
  )
  plot_regulon_network(
    regulons = regulons,
    rss_mat = rss_mat
  )

  # ----- save outputs -----
  saveRDS(seurat_obj_out, file.path(DIR_RDS, "seurat_with_scenic.rds"))
  if (isTRUE(EXPORT_FULL_AUC_MATRIX_CSV)) {
    save_matrix_csv(auc_mat, file.path(DIR_TABLE, "regulon_auc_matrix.csv"))
  } else {
    msg("[INFO] Skipping full regulon_auc_matrix.csv export (EXPORT_FULL_AUC_MATRIX_CSV=FALSE).\n")
  }

  top_reg_table <- data.table(
    cell = colnames(auc_mat),
    top_regulon = seurat_obj_out$scenic_top_regulon[colnames(auc_mat)],
    top_regulon_assay_feature = seurat_obj_out$scenic_top_regulon_assay_feature[colnames(auc_mat)],
    top_regulon_auc = seurat_obj_out$scenic_top_regulon_auc[colnames(auc_mat)],
    cell_type = as.character(seurat_obj_out@meta.data[colnames(auc_mat), CELL_TYPE_COL]),
    sample = as.character(seurat_obj_out@meta.data[colnames(auc_mat), SAMPLE_COL])
  )
  fwrite(top_reg_table, file.path(DIR_TABLE, "cellwise_top_regulon.csv"))

  result <- list(
    seurat_obj = seurat_obj_out,
    scenic_options = scenic_options,
    regulons = regulons,
    auc_mat = auc_mat,
    rss_mat = rss_mat,
    inference_cells = infer_cells
  )

  saveRDS(result, file.path(DIR_RDS, "scenic_full_result.rds"))

  msg("============================================================\n")
  msg("SCENIC module completed\n")
  msg("Final cells : %d\n", ncol(seurat_obj_out))
  msg("Regulons    : %d\n", nrow(auc_mat))
  msg("Output dir  : %s\n", OUTPUT_DIR)
  msg("============================================================\n")

  invisible(result)
}

# ==============================================================================
# 7. Execute
# ==============================================================================
if (!isTRUE(SCENIC_SOURCE_ONLY)) {
  if (!file.exists(INPUT_RDS)) {
    stop(sprintf("INPUT_RDS not found: %s", INPUT_RDS))
  }
  if (!dir.exists(DATABASE_DIR)) {
    stop(sprintf("DATABASE_DIR not found: %s", DATABASE_DIR))
  }
  msg("[LOAD] %s\n", INPUT_RDS)
  seurat_obj_input <- readRDS(INPUT_RDS)

  result <- run_scenic_module(seurat_obj_input)
}