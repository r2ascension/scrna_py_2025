#!/usr/bin/env Rscript
# ==============================================================================
# Stromal scHPL Downstream Postprocess Pipeline v1.0
# ==============================================================================
#
# Purpose:
#   Split the downstream post-processing out of the original notebook and keep
#   the notebook focused on scHPL inference only.
#
# Methodology note:
#   This R script belongs to the downstream-summary layer of the stromal schpl
#   workflow. It does not define the main scHPL inference method; instead it
#   audits outputs, compares scANVI vs scHPL labels, and prepares report-ready
#   artifacts. For the newer branch-wise stromal workflow, the primary reject
#   interpretation step is `stromal_schpl_reject_followup_20260408_v1_0.py`.
#
# Main tasks:
#   1. Read stromal scHPL h5ad -> Seurat via local GetSeurat()
#   2. Audit AnnData structure (obs / layers / raw / obsm)
#   3. Compare scANVI vs scHPL predictions and export summary CSVs
#   4. Optionally run CHOIR on the latent space (scanvi/scvi)
#   5. Generate REPORT.md + LLM-ready JSON/prompt (+ optional DeepSeek call)
#   6. Save h5ad safely by preserving original layers/raw/obsm and only writing
#      back cleaned metadata / uns fields
#
# Smoke test:
#   PIPELINE_TEST_MODE=true Rscript stromal_schpl_postprocess_v1_0_20260401.R
#   -> exits after load + audit + comparison CSV generation
#
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

PIPELINE_VERSION <- "v1.0"
PIPELINE_DATE    <- "2026-04-01"
PIPELINE_TAG     <- "stromal_schpl_postprocess_v1_0_20260401"

INPUT_H5AD <- "/home/h2048/data/py/0330/stromal_schpl_v1_0/adata_stromal_query_schpl_v1_0.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0401/stromal_schpl_postprocess_v1_0"
LOCAL_GETSEURAT <- "/home/h2048/script/R/GetSeurat.R"
PYTHON_BIN <- "/home/h2048/miniconda3/envs/scarches_stable/bin/python"
PYTHON_CONDA_ENV <- "scarches_stable"
CONDA_BIN <- "/home/h2048/miniconda3/bin/conda"
DOTENV_PATH <- "/home/h2048/.env"

PIPELINE_TEST_MODE <- identical(tolower(Sys.getenv("PIPELINE_TEST_MODE", "false")), "true")
RUN_CHOIR <- identical(tolower(Sys.getenv("RUN_CHOIR", "true")), "true")
RUN_OPTIONAL_LLM_CALL <- identical(tolower(Sys.getenv("RUN_OPTIONAL_LLM_CALL", "false")), "true")
CHOIR_ALPHA <- 0.05
CHOIR_N_CORES <- 4
DEEPSEEK_MODEL <- "deepseek-chat"
LLM_TIMEOUT_SEC <- 120

KEYS <- list(
  sample = "sample",
  query_subset = "query_subset",
  scanvi_pred = "cell_type_scarches_pred",
  scanvi_final = "cell_type_scarches_final",
  scanvi_conf = "scarches_confidence",
  schpl_raw = "schpl_pred_raw",
  schpl_pred = "schpl_pred",
  schpl_prob = "schpl_prob",
  schpl_rejected = "schpl_rejected",
  schpl_reject_type = "schpl_reject_type",
  leiden = "leiden_schpl_qc",
  novel = "schpl_novel_candidate"
)

UMAP_REDUCTION_PREFERRED <- c("umap", "umapscanvi", "umapscvi")
CHOIR_REDUCTION_CANDIDATES <- c("scanvi", "scvi", "pca")

FIG_DIR <- file.path(OUTPUT_DIR, "figures")
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")
CHOIR_DIR <- file.path(REPORT_DIR, "choir")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHOIR_DIR, recursive = TRUE, showWarnings = FALSE)

ARTIFACTS <- list(
  structure_json = file.path(REPORT_DIR, "anndata_structure_v1_0.json"),
  comparison_csv = file.path(REPORT_DIR, "schpl_label_comparison_v1_0.csv"),
  agreement_csv = file.path(REPORT_DIR, "schpl_agreement_summary_v1_0.csv"),
  disagreement_csv = file.path(REPORT_DIR, "schpl_disagreement_pairs_v1_0.csv"),
  rejected_scanvi_csv = file.path(REPORT_DIR, "schpl_rejected_scanvi_distribution_v1_0.csv"),
  confidence_csv = file.path(REPORT_DIR, "schpl_confidence_summary_v1_0.csv"),
  cluster_summary_csv = file.path(REPORT_DIR, "cluster_rejection_summary_v1_0.csv"),
  choir_clusters_csv = file.path(CHOIR_DIR, "choir_clusters.csv"),
  choir_context_json = file.path(CHOIR_DIR, "choir_ready_context_v1_0.json"),
  llm_input_json = file.path(REPORT_DIR, "llm_input_v1_0.json"),
  llm_prompt_md = file.path(REPORT_DIR, "LLM_PROMPT_v1_0.md"),
  llm_response_md = file.path(REPORT_DIR, "LLM_RESPONSE_v1_0.md"),
  report_md = file.path(OUTPUT_DIR, "REPORT.md"),
  output_h5ad = file.path(OUTPUT_DIR, "adata_stromal_query_schpl_postprocess_v1_0.h5ad"),
  summary_json = file.path(OUTPUT_DIR, "stromal_schpl_postprocess_summary_v1_0.json")
)

load_dotenv_if_present <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!startsWith(lines, "#")]
  for (ln in lines) {
    kv <- strsplit(ln, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) next
    key <- trimws(kv[1])
    val <- trimws(paste(kv[-1], collapse = "="))
    if ((startsWith(val, '"') && endsWith(val, '"')) ||
        (startsWith(val, "'") && endsWith(val, "'"))) {
      val <- substring(val, 2, nchar(val) - 1)
    }
    if (!nzchar(Sys.getenv(key, ""))) Sys.setenv(structure(val, names = key))
  }
  invisible(TRUE)
}

load_dotenv_if_present(DOTENV_PATH)
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY", "")
ENABLE_LLM <- nchar(DEEPSEEK_API_KEY) >= 10

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
})

CHOIR_AVAILABLE <- requireNamespace("CHOIR", quietly = TRUE)
HTTR_AVAILABLE  <- requireNamespace("httr", quietly = TRUE)

configure_python_bridge <- function(py_bin, conda_env, conda_bin) {
  if (nzchar(Sys.getenv("RETICULATE_PYTHON", "")) && file.exists(Sys.getenv("RETICULATE_PYTHON"))) {
    return(Sys.getenv("RETICULATE_PYTHON"))
  }
  if (file.exists(py_bin)) {
    Sys.setenv(RETICULATE_PYTHON = py_bin)
    return(py_bin)
  }

  if (file.exists(conda_bin)) {
    fallback_python <- file.path(dirname(conda_bin), "..", "envs", conda_env, "bin", "python")
    fallback_python <- normalizePath(fallback_python, winslash = "/", mustWork = FALSE)
    if (file.exists(fallback_python)) {
      Sys.setenv(RETICULATE_PYTHON = fallback_python)
      return(fallback_python)
    }
  }

  stop(sprintf(
    "Failed to configure reticulate Python. Tried PYTHON_BIN=%s and conda env=%s",
    py_bin, conda_env
  ))
}

python_configured <- configure_python_bridge(PYTHON_BIN, PYTHON_CONDA_ENV, CONDA_BIN)

anndata  <- reticulate::import("anndata", convert = FALSE)
pd       <- reticulate::import("pandas", convert = FALSE)
builtins <- reticulate::import("builtins", convert = FALSE)

if (!file.exists(LOCAL_GETSEURAT)) stop(sprintf("Missing helper: %s", LOCAL_GETSEURAT))
source(LOCAL_GETSEURAT)

safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", as.character(x))

save_plot <- function(plot_obj, path_no_ext, width = 10, height = 8) {
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot_obj, width = width, height = height, dpi = 300)
}

pick_reduction <- function(obj, preferred) {
  reds <- Seurat::Reductions(obj)
  hit  <- preferred[preferred %in% reds]
  if (length(hit) > 0) hit[[1]] else NULL
}

coerce_to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  y <- trimws(tolower(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

sanitize_obs_for_h5ad <- function(df) {
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  if ("_index" %in% colnames(out)) {
    colnames(out)[colnames(out) == "_index"] <- "orig_index"
  }
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
        if (length(x) == 0 || all(is.na(x))) return(NA_character_)
        paste(as.character(x), collapse = "; ")
      }, character(1))
    } else if (is.character(vec)) {
      out[[col]] <- trimws(vec)
    } else if (!(is.integer(vec) || is.numeric(vec))) {
      out[[col]] <- as.character(vec)
    }
  }
  out
}

write_json_pretty <- function(x, path) {
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

audit_anndata_structure <- function(adata) {
  layers <- reticulate::py_to_r(builtins$list(adata$layers$keys()))
  obsm   <- reticulate::py_to_r(builtins$list(adata$obsm$keys()))
  uns    <- reticulate::py_to_r(builtins$list(adata$uns$keys()))
  obs_columns <- reticulate::py_to_r(builtins$list(adata$obs$columns))
  raw_info <- NULL
  if (!reticulate::py_is_null_xptr(adata$raw)) {
    raw_info <- list(
      n_obs = reticulate::py_to_r(adata$raw$n_obs),
      n_vars = reticulate::py_to_r(adata$raw$n_vars)
    )
  }
  list(
    shape = list(
      n_obs = reticulate::py_to_r(adata$n_obs),
      n_vars = reticulate::py_to_r(adata$n_vars)
    ),
    layers = as.list(layers),
    raw = raw_info,
    obsm = as.list(obsm),
    obs_columns = as.list(obs_columns),
    uns_keys = as.list(uns)
  )
}

validate_required_columns <- function(cols, required) {
  missing <- setdiff(required, cols)
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

call_deepseek_chat <- function(prompt, api_key, model = DEEPSEEK_MODEL, timeout_sec = LLM_TIMEOUT_SEC) {
  if (!HTTR_AVAILABLE) {
    return(list(status = "skipped", reason = "httr not installed"))
  }
  if (!nzchar(api_key)) {
    return(list(status = "skipped", reason = "DEEPSEEK_API_KEY not available"))
  }
  resp <- tryCatch(
    httr::POST(
      url = "https://api.deepseek.com/v1/chat/completions",
      httr::add_headers(
        Authorization = paste("Bearer", api_key),
        `Content-Type` = "application/json"
      ),
      body = list(
        model = model,
        messages = list(
          list(role = "system", content = "你是一个谨慎的单细胞分析助手。请用简体中文 Markdown 输出。"),
          list(role = "user", content = prompt)
        ),
        temperature = 0.2
      ),
      encode = "json",
      httr::timeout(timeout_sec)
    ),
    error = function(e) e
  )
  if (inherits(resp, "error")) {
    return(list(status = "error", reason = conditionMessage(resp)))
  }
  if (httr::http_error(resp)) {
    return(list(
      status = "error",
      reason = sprintf("HTTP %s: %s", httr::status_code(resp), httr::content(resp, "text", encoding = "UTF-8"))
    ))
  }
  payload <- httr::content(resp, as = "parsed", simplifyVector = FALSE)
  content <- payload$choices[[1]]$message$content
  list(status = "ok", content = content, raw = payload)
}

update_obs_in_anndata <- function(adata, meta_export) {
  obs_names <- as.character(reticulate::py_to_r(builtins$list(adata$obs_names)))
  if (!all(obs_names %in% rownames(meta_export))) {
    stop("Some AnnData obs_names are missing from metadata export")
  }
  meta_export <- meta_export[obs_names, , drop = FALSE]
  for (col in colnames(meta_export)) {
    adata$obs$`__setitem__`(col, reticulate::r_to_py(unname(meta_export[[col]])))
  }
  invisible(meta_export)
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("Stromal scHPL downstream postprocess %s\n", PIPELINE_VERSION))
cat(sprintf("Date: %s\n", PIPELINE_DATE))
cat(strrep("=", 70), "\n\n", sep = "")
cat(sprintf("[INFO] Input h5ad: %s\n", INPUT_H5AD))
cat(sprintf("[INFO] Output dir : %s\n", OUTPUT_DIR))
cat(sprintf("[INFO] LLM enabled: %s (run=%s)\n", ENABLE_LLM, RUN_OPTIONAL_LLM_CALL))
cat(sprintf("[INFO] CHOIR available: %s (run=%s)\n\n", CHOIR_AVAILABLE, RUN_CHOIR))

if (!file.exists(INPUT_H5AD)) stop(sprintf("Input h5ad not found: %s", INPUT_H5AD))

cat("=== Step 1: Audit AnnData structure ===\n")
adata_audit <- anndata$read_h5ad(INPUT_H5AD, backed = "r")
structure_info <- audit_anndata_structure(adata_audit)
required_obs <- unname(unlist(KEYS))
validate_required_columns(unlist(structure_info$obs_columns), required_obs)
if (!("counts" %in% unlist(structure_info$layers))) {
  stop("Input h5ad is missing layers['counts']")
}
if (is.null(structure_info$raw)) {
  stop("Input h5ad is missing .raw; refusing to continue because layer alignment matters")
}
write_json_pretty(structure_info, ARTIFACTS$structure_json)
cat(sprintf("[OK] Structure JSON saved: %s\n\n", ARTIFACTS$structure_json))
rm(adata_audit)
gc()

cat("=== Step 2: Load h5ad -> Seurat ===\n")
obj <- GetSeurat(
  h5ad_path = INPUT_H5AD,
  prefer_raw = FALSE,
  prefer_layer_counts = TRUE,
  validate_counts = TRUE,
  debug = TRUE
)
cat(sprintf("[OK] Seurat object: %d genes x %d cells\n", nrow(obj), ncol(obj)))
counts_ok <- tryCatch({
  tmp <- Seurat::LayerData(obj, assay = "RNA", layer = "counts")
  !is.null(tmp)
}, error = function(e) FALSE)
if (!counts_ok) stop("RNA assay missing counts layer after GetSeurat()")
cat("[OK] counts layer verified in Seurat\n\n")

meta <- obj@meta.data
validate_required_columns(colnames(meta), required_obs)

meta[[KEYS$schpl_rejected]] <- coerce_to_logical(meta[[KEYS$schpl_rejected]])
meta[[KEYS$novel]] <- coerce_to_logical(meta[[KEYS$novel]])
meta[[KEYS$scanvi_conf]] <- suppressWarnings(as.numeric(meta[[KEYS$scanvi_conf]]))
meta[[KEYS$schpl_prob]] <- suppressWarnings(as.numeric(meta[[KEYS$schpl_prob]]))
obj@meta.data <- meta

cat("=== Step 3: Compare scANVI vs scHPL ===\n")
comparison_df <- data.frame(
  cell = rownames(meta),
  sample = meta[[KEYS$sample]],
  query_subset = meta[[KEYS$query_subset]],
  cell_type_scarches_pred = meta[[KEYS$scanvi_pred]],
  cell_type_scarches_final = meta[[KEYS$scanvi_final]],
  scarches_confidence = meta[[KEYS$scanvi_conf]],
  schpl_pred_raw = meta[[KEYS$schpl_raw]],
  schpl_pred = meta[[KEYS$schpl_pred]],
  schpl_prob = meta[[KEYS$schpl_prob]],
  schpl_rejected = meta[[KEYS$schpl_rejected]],
  schpl_reject_type = meta[[KEYS$schpl_reject_type]],
  leiden_schpl_qc = meta[[KEYS$leiden]],
  schpl_novel_candidate = meta[[KEYS$novel]],
  stringsAsFactors = FALSE,
  row.names = rownames(meta)
)

accepted_df <- comparison_df[!comparison_df$schpl_rejected, , drop = FALSE]
rejected_df <- comparison_df[comparison_df$schpl_rejected, , drop = FALSE]
agreed <- accepted_df$cell_type_scarches_pred == accepted_df$schpl_pred
agree_rate <- if (nrow(accepted_df) == 0) NA_real_ else mean(agreed) * 100

disagreement_pairs <- if (nrow(accepted_df) > 0) {
  as.data.frame(
    comparison_df[!comparison_df$schpl_rejected, c("cell_type_scarches_pred", "schpl_pred")] |>
      table(useNA = "ifany"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::rename(n = Freq) |>
    dplyr::filter(n > 0) |>
    dplyr::arrange(dplyr::desc(n))
} else {
  data.frame(cell_type_scarches_pred = character(), schpl_pred = character(), n = integer())
}

rejected_scanvi_dist <- if (nrow(rejected_df) > 0) {
  rejected_df |>
    dplyr::count(cell_type_scarches_final, name = "n") |>
    dplyr::mutate(pct = round(n / sum(n) * 100, 2)) |>
    dplyr::arrange(dplyr::desc(n))
} else {
  data.frame(cell_type_scarches_final = character(), n = integer(), pct = numeric())
}

confidence_summary <- comparison_df |>
  dplyr::group_by(schpl_pred) |>
  dplyr::summarise(
    count = dplyr::n(),
    scarches_conf_mean = mean(scarches_confidence, na.rm = TRUE),
    scarches_conf_median = median(scarches_confidence, na.rm = TRUE),
    schpl_prob_mean = mean(schpl_prob, na.rm = TRUE),
    schpl_prob_median = median(schpl_prob, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(count))

cluster_rej <- comparison_df |>
  dplyr::group_by(leiden_schpl_qc) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_rejected = sum(schpl_rejected, na.rm = TRUE),
    pct_rejected = round(100 * mean(schpl_rejected, na.rm = TRUE), 2),
    n_novel_candidate = sum(schpl_novel_candidate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(pct_rejected), dplyr::desc(n_rejected))

agreement_summary <- data.frame(
  total_cells = nrow(comparison_df),
  accepted_cells = nrow(accepted_df),
  rejected_cells = nrow(rejected_df),
  agreement_pct = round(agree_rate, 4),
  novel_candidate_cells = sum(comparison_df$schpl_novel_candidate, na.rm = TRUE),
  novel_candidate_clusters = sum(cluster_rej$n_novel_candidate > 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)

fwrite(comparison_df, ARTIFACTS$comparison_csv)
fwrite(agreement_summary, ARTIFACTS$agreement_csv)
fwrite(disagreement_pairs, ARTIFACTS$disagreement_csv)
fwrite(rejected_scanvi_dist, ARTIFACTS$rejected_scanvi_csv)
fwrite(confidence_summary, ARTIFACTS$confidence_csv)
fwrite(cluster_rej, ARTIFACTS$cluster_summary_csv)
cat(sprintf("[OK] Agreement summary saved: %s\n", ARTIFACTS$agreement_csv))
cat(sprintf("[OK] Disagreement pairs saved: %s\n", ARTIFACTS$disagreement_csv))
cat(sprintf("[OK] Cluster rejection summary saved: %s\n\n", ARTIFACTS$cluster_summary_csv))

cat("=== Step 4: Visualizations ===\n")
umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
generated_figures <- character()
if (!is.null(umap_reduction)) {
  p_scanvi <- Seurat::DimPlot(
    obj, reduction = umap_reduction, group.by = KEYS$scanvi_final,
    pt.size = 0.3, shuffle = TRUE, label = FALSE
  ) + ggplot2::ggtitle("scANVI final label") + ggplot2::theme_classic(base_size = 12)
  scanvi_path <- file.path(FIG_DIR, "umap_scanvi_final")
  save_plot(p_scanvi, scanvi_path, width = 10, height = 8)
  generated_figures <- c(generated_figures, paste0(scanvi_path, c(".pdf", ".png")))

  p_schpl <- Seurat::DimPlot(
    obj, reduction = umap_reduction, group.by = KEYS$schpl_pred,
    pt.size = 0.3, shuffle = TRUE, label = FALSE
  ) + ggplot2::ggtitle("scHPL prediction") + ggplot2::theme_classic(base_size = 12)
  schpl_path <- file.path(FIG_DIR, "umap_schpl_pred")
  save_plot(p_schpl, schpl_path, width = 10, height = 8)
  generated_figures <- c(generated_figures, paste0(schpl_path, c(".pdf", ".png")))

  p_rej <- Seurat::DimPlot(
    obj, reduction = umap_reduction, group.by = KEYS$schpl_reject_type,
    pt.size = 0.3, shuffle = TRUE, label = FALSE
  ) + ggplot2::ggtitle("scHPL rejection type") + ggplot2::theme_classic(base_size = 12)
  rej_path <- file.path(FIG_DIR, "umap_schpl_reject_type")
  save_plot(p_rej, rej_path, width = 10, height = 8)
  generated_figures <- c(generated_figures, paste0(rej_path, c(".pdf", ".png")))
}

p_cluster <- ggplot(cluster_rej[seq_len(min(20, nrow(cluster_rej))), , drop = FALSE],
                    aes(x = reorder(as.character(leiden_schpl_qc), pct_rejected), y = pct_rejected)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = "Leiden cluster", y = "% rejected", title = "Top clusters by scHPL rejection rate")
cluster_path <- file.path(FIG_DIR, "cluster_rejection_rate_top20")
save_plot(p_cluster, cluster_path, width = 9, height = 7)
generated_figures <- c(generated_figures, paste0(cluster_path, c(".pdf", ".png")))
cat(sprintf("[OK] Figures saved under: %s\n\n", FIG_DIR))

choir_col <- NULL
choir_status <- "skipped"
choir_note <- "CHOIR not run"
if (RUN_CHOIR && CHOIR_AVAILABLE) {
  cat("=== Step 5: CHOIR clustering ===\n")
  choir_reduction <- pick_reduction(obj, CHOIR_REDUCTION_CANDIDATES)
  if (!is.null(choir_reduction)) {
    choir_embedding <- tryCatch(as.matrix(Seurat::Embeddings(obj, reduction = choir_reduction)), error = function(e) NULL)
    if (!is.null(choir_embedding)) {
      choir_var_features <- tryCatch(Seurat::VariableFeatures(obj), error = function(e) character())
      if (length(choir_var_features) == 0) choir_var_features <- rownames(obj)
      obj_choir <- tryCatch(
        CHOIR::CHOIR(
          obj,
          use_assay = "RNA",
          reduction = choir_embedding,
          var_features = choir_var_features,
          n_cores = CHOIR_N_CORES,
          alpha = CHOIR_ALPHA,
          random_seed = 42
        ),
        error = function(e) {
          message(sprintf("[WARN] CHOIR failed: %s", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(obj_choir)) obj <- obj_choir
      choir_col <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
      if (!choir_col %in% colnames(obj@meta.data)) {
        fallback <- grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
        choir_col <- if (length(fallback) > 0) fallback[[1]] else NULL
      }
      if (!is.null(choir_col)) {
        choir_status <- "completed"
        choir_note <- sprintf("CHOIR completed using reduction '%s'", choir_reduction)
        choir_export <- obj@meta.data |>
          tibble::rownames_to_column("cell") |>
          dplyr::select(
            cell,
            dplyr::all_of(c(KEYS$sample, KEYS$query_subset, KEYS$scanvi_final,
                            KEYS$schpl_pred, KEYS$schpl_rejected, KEYS$novel, choir_col))
          )
        data.table::fwrite(choir_export, ARTIFACTS$choir_clusters_csv)
        if (!is.null(umap_reduction)) {
          p_choir <- Seurat::DimPlot(
            obj, reduction = umap_reduction, group.by = choir_col,
            pt.size = 0.3, shuffle = TRUE, label = TRUE
          ) + ggplot2::ggtitle(sprintf("CHOIR clusters (alpha=%.2f)", CHOIR_ALPHA)) +
            ggplot2::theme_classic(base_size = 12)
          choir_path <- file.path(FIG_DIR, "choir_umap")
          save_plot(p_choir, choir_path, width = 12, height = 9)
          generated_figures <- c(generated_figures, paste0(choir_path, c(".pdf", ".png")))
        }
      } else {
        choir_status <- "failed"
        choir_note <- "CHOIR ran but no CHOIR_clusters_* column was found"
      }
    } else {
      choir_status <- "failed"
      choir_note <- sprintf("Reduction '%s' could not be extracted for CHOIR", choir_reduction)
    }
  } else {
    choir_status <- "skipped"
    choir_note <- "No suitable scanvi/scvi/pca reduction found for CHOIR"
  }
} else if (!CHOIR_AVAILABLE) {
  choir_note <- "CHOIR package not installed in current R environment"
}

choir_context <- list(
  status = choir_status,
  note = choir_note,
  choir_column = choir_col,
  preferred_reductions = CHOIR_REDUCTION_CANDIDATES,
  available_reductions = as.list(Seurat::Reductions(obj)),
  output_csv = if (file.exists(ARTIFACTS$choir_clusters_csv)) ARTIFACTS$choir_clusters_csv else NULL
)
write_json_pretty(choir_context, ARTIFACTS$choir_context_json)
cat(sprintf("[OK] CHOIR context saved: %s\n\n", ARTIFACTS$choir_context_json))

summary_payload <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_h5ad = INPUT_H5AD,
  output_dir = OUTPUT_DIR,
  structure = structure_info,
  metrics = list(
    total_cells = nrow(comparison_df),
    accepted_cells = nrow(accepted_df),
    rejected_cells = nrow(rejected_df),
    agreement_pct = round(agree_rate, 4),
    novel_candidate_cells = sum(comparison_df$schpl_novel_candidate, na.rm = TRUE),
    novel_candidate_clusters = sum(cluster_rej$n_novel_candidate > 0, na.rm = TRUE)
  ),
  top_disagreement_pairs = utils::head(disagreement_pairs, 15),
  top_rejection_clusters = utils::head(cluster_rej, 10),
  choir = choir_context,
  figures = generated_figures
)
write_json_pretty(summary_payload, ARTIFACTS$llm_input_json)

llm_prompt <- paste(
  "# Stromal scHPL downstream interpretation prompt",
  "",
  "请基于以下结构化结果，用简体中文总结 stromal query 的映射质量、潜在新状态证据、以及后续验证建议。",
  "请重点关注：",
  "1. scANVI 与 scHPL 在接受细胞上的一致性；",
  "2. 哪些 cluster 的 rejection rate 高，是否具有一致的 scANVI 标签背景；",
  "3. 是否值得继续做 CHOIR / marker review / DE / pathway 分析；",
  "4. 明确指出技术噪音与真实 novel state 的边界。",
  "",
  "## 输入文件",
  sprintf("- 结构化 JSON: `%s`", ARTIFACTS$llm_input_json),
  sprintf("- 比较总表: `%s`", ARTIFACTS$comparison_csv),
  sprintf("- disagreement pairs: `%s`", ARTIFACTS$disagreement_csv),
  sprintf("- cluster rejection summary: `%s`", ARTIFACTS$cluster_summary_csv),
  sprintf("- CHOIR context: `%s`", ARTIFACTS$choir_context_json),
  "",
  "## 期望输出",
  "1. 一段执行摘要；",
  "2. 关键证据点（novel candidate / disagreement / rejection hotspots）；",
  "3. 主要 caveats；",
  "4. 推荐的下一步分析。",
  sep = "\n"
)
writeLines(llm_prompt, ARTIFACTS$llm_prompt_md, useBytes = TRUE)

llm_result <- list(status = "skipped", reason = "RUN_OPTIONAL_LLM_CALL=false or API unavailable")
if (ENABLE_LLM && RUN_OPTIONAL_LLM_CALL) {
  llm_result <- call_deepseek_chat(llm_prompt, DEEPSEEK_API_KEY)
}
if (identical(llm_result$status, "ok")) {
  writeLines(as.character(llm_result$content), ARTIFACTS$llm_response_md, useBytes = TRUE)
} else {
  writeLines(
    paste0("LLM call skipped or failed.\n\n", jsonlite::toJSON(llm_result, pretty = TRUE, auto_unbox = TRUE)),
    ARTIFACTS$llm_response_md,
    useBytes = TRUE
  )
}
cat(sprintf("[OK] LLM input JSON: %s\n", ARTIFACTS$llm_input_json))
cat(sprintf("[OK] LLM prompt MD : %s\n", ARTIFACTS$llm_prompt_md))
cat(sprintf("[OK] LLM response  : %s\n\n", ARTIFACTS$llm_response_md))

report_lines <- c(
  "# Stromal scHPL Downstream Report",
  "",
  sprintf("- Generated at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Input h5ad: `%s`", INPUT_H5AD),
  sprintf("- Output dir: `%s`", OUTPUT_DIR),
  sprintf("- UMAP reduction: `%s`", ifelse(is.null(umap_reduction), "(not available)", umap_reduction)),
  sprintf("- CHOIR status: `%s`", choir_status),
  sprintf("- CHOIR note: %s", choir_note),
  "",
  "## Dataset overview",
  "",
  sprintf("- Cells / HVGs: %d / %d", nrow(comparison_df), nrow(obj)),
  sprintf("- layers: %s", paste(unlist(structure_info$layers), collapse = ", ")),
  sprintf("- raw present: %s", !is.null(structure_info$raw)),
  sprintf("- obsm: %s", paste(unlist(structure_info$obsm), collapse = ", ")),
  "",
  "## Key metrics",
  "",
  sprintf("- Accepted cells: %d", nrow(accepted_df)),
  sprintf("- Rejected cells: %d", nrow(rejected_df)),
  sprintf("- Agreement (accepted only): %.2f%%", agree_rate),
  sprintf("- Novel candidate cells: %d", sum(comparison_df$schpl_novel_candidate, na.rm = TRUE)),
  sprintf("- Novel candidate clusters: %d", sum(cluster_rej$n_novel_candidate > 0, na.rm = TRUE)),
  "",
  "## Top rejection-enriched clusters",
  ""
)
if (nrow(cluster_rej) > 0) {
  for (i in seq_len(min(10, nrow(cluster_rej)))) {
    row <- cluster_rej[i, ]
    report_lines <- c(
      report_lines,
      sprintf(
        "- Cluster %s: %d/%d rejected (%.2f%%), novel candidates=%d",
        as.character(row$leiden_schpl_qc), row$n_rejected, row$n_cells, row$pct_rejected, row$n_novel_candidate
      )
    )
  }
}
report_lines <- c(
  report_lines,
  "",
  "## Artifacts",
  "",
  sprintf("- comparison CSV: `%s`", ARTIFACTS$comparison_csv),
  sprintf("- agreement CSV: `%s`", ARTIFACTS$agreement_csv),
  sprintf("- disagreement CSV: `%s`", ARTIFACTS$disagreement_csv),
  sprintf("- confidence CSV: `%s`", ARTIFACTS$confidence_csv),
  sprintf("- cluster summary CSV: `%s`", ARTIFACTS$cluster_summary_csv),
  sprintf("- CHOIR context JSON: `%s`", ARTIFACTS$choir_context_json),
  sprintf("- structure JSON: `%s`", ARTIFACTS$structure_json),
  sprintf("- LLM input JSON: `%s`", ARTIFACTS$llm_input_json),
  sprintf("- output h5ad: `%s`", ARTIFACTS$output_h5ad),
  "",
  "## Figures",
  ""
)
if (length(generated_figures) > 0) {
  report_lines <- c(report_lines, sprintf("- `%s`", generated_figures))
} else {
  report_lines <- c(report_lines, "- No figures generated.")
}
writeLines(report_lines, ARTIFACTS$report_md, useBytes = TRUE)
cat(sprintf("[OK] Report saved: %s\n\n", ARTIFACTS$report_md))

if (PIPELINE_TEST_MODE) {
  cat("[TEST MODE] Stopping after load + audit + comparison/report generation.\n")
  quit(save = "no", status = 0)
}

cat("=== Step 6: Safe h5ad export (preserve layers/raw/obsm) ===\n")
meta_export <- sanitize_obs_for_h5ad(obj@meta.data)
adata_out <- anndata$read_h5ad(INPUT_H5AD)
updated_meta <- update_obs_in_anndata(adata_out, meta_export)
adata_out$uns$`__setitem__`(
  "stromal_postprocess_r",
  reticulate::r_to_py(list(
    version = PIPELINE_VERSION,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    report_md = ARTIFACTS$report_md,
    llm_input_json = ARTIFACTS$llm_input_json,
    choir_context_json = ARTIFACTS$choir_context_json,
    agreement_pct = round(agree_rate, 4)
  ))
)
adata_out$write_h5ad(ARTIFACTS$output_h5ad, compression = "gzip")

adata_check <- anndata$read_h5ad(ARTIFACTS$output_h5ad, backed = "r")
check_layers <- as.character(reticulate::py_to_r(builtins$list(adata_check$layers$keys())))
check_obsm <- as.character(reticulate::py_to_r(builtins$list(adata_check$obsm$keys())))
check_obs_cols <- as.character(reticulate::py_to_r(builtins$list(adata_check$obs$columns)))
check_raw <- !reticulate::py_is_null_xptr(adata_check$raw)
validate_required_columns(check_obs_cols, required_obs)
if (!("counts" %in% check_layers)) stop("Reloaded h5ad is missing layers['counts']")
if (!check_raw) stop("Reloaded h5ad lost .raw")
if (!("X_scanvi" %in% check_obsm)) stop("Reloaded h5ad is missing obsm['X_scanvi']")
if (!("X_umap" %in% check_obsm)) warning("Reloaded h5ad does not contain obsm['X_umap']")
cat(sprintf("[OK] h5ad verified: %d cells x %d genes, layers=%s\n\n",
            reticulate::py_to_r(adata_check$n_obs),
            reticulate::py_to_r(adata_check$n_vars),
            paste(check_layers, collapse = ", ")))

final_summary <- list(
  version = PIPELINE_VERSION,
  input_h5ad = INPUT_H5AD,
  output_h5ad = ARTIFACTS$output_h5ad,
  report_md = ARTIFACTS$report_md,
  llm_input_json = ARTIFACTS$llm_input_json,
  llm_response_md = ARTIFACTS$llm_response_md,
  choir_context_json = ARTIFACTS$choir_context_json,
  generated_figures = generated_figures,
  metrics = agreement_summary,
  layers_verified = check_layers,
  obsm_verified = check_obsm,
  raw_preserved = check_raw
)
write_json_pretty(final_summary, ARTIFACTS$summary_json)
cat(sprintf("[OK] Summary JSON saved: %s\n\n", ARTIFACTS$summary_json))

cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("%s COMPLETE\n", toupper(PIPELINE_TAG)))
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("- output h5ad : %s\n", ARTIFACTS$output_h5ad))
cat(sprintf("- report md   : %s\n", ARTIFACTS$report_md))
cat(sprintf("- choir state : %s\n", choir_status))
cat(sprintf("- agreement   : %.2f%%\n", agree_rate))
cat(sprintf("- rejected    : %d\n", nrow(rejected_df)))
cat(sprintf("- novel cells : %d\n", sum(comparison_df$schpl_novel_candidate, na.rm = TRUE)))
