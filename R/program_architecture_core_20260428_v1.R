#!/usr/bin/env Rscript
# ==============================================================================
# Program Architecture Core Helper (2026-04-28 v1)
# ==============================================================================

PA_HELPER_VERSION_20260428_V1 <- "20260428_v1"

pa_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

pa_safe_trim <- function(x) {
  out <- trimws(as.character(x))
  out[is.na(out)] <- ""
  out
}

pa_unique_chr <- function(x) {
  x <- pa_safe_trim(x)
  x[nzchar(x)] |> unique()
}

pa_scalar_chr <- function(x, arg_name) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(as.character(x)))) {
    stop(sprintf("%s must be one non-empty string", arg_name), call. = FALSE)
  }
  as.character(x)
}

pa_validate_required_columns <- function(df, required_cols, object_name = "data.frame") {
  if (!is.data.frame(df)) {
    stop(sprintf("%s must be a data.frame", object_name), call. = FALSE)
  }
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "%s is missing required columns: %s",
      object_name,
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

pa_prepare_output_dir <- function(path) {
  path <- pa_scalar_chr(path, "path")
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

pa_stringify_value <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.list(x)) {
    if (length(x) == 0L) return(NA_character_)
    if (is.data.frame(x)) {
      return(sprintf("<data.frame:%d x %d>", nrow(x), ncol(x)))
    }
    if (all(vapply(x, function(item) length(item) <= 1L && !is.list(item), logical(1)))) {
      vals <- vapply(x, function(item) if (length(item) == 0L || is.null(item)) NA_character_ else as.character(item[[1]]), character(1))
      vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
      return(if (length(vals) == 0L) NA_character_ else paste(vals, collapse = "|"))
    }
    return("<list>")
  }
  if (length(x) == 0L) return(NA_character_)
  vals <- as.character(x)
  vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
  if (length(vals) == 0L) NA_character_ else paste(vals, collapse = "|")
}

pa_tsv_safe_df <- function(df) {
  if (!is.data.frame(df)) {
    stop("pa_tsv_safe_df expects a data.frame", call. = FALSE)
  }
  out <- df
  for (nm in colnames(out)) {
    if (is.list(out[[nm]])) {
      out[[nm]] <- vapply(out[[nm]], pa_stringify_value, character(1))
    }
  }
  out
}

pa_write_tsv <- function(df, path) {
  pa_validate_required_columns(data.frame(dummy = 1), character(), "dummy")
  path <- pa_scalar_chr(path, "path")
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    pa_tsv_safe_df(df),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

pa_write_markdown <- function(lines, path) {
  path <- pa_scalar_chr(path, "path")
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(lines), con = path, useBytes = TRUE)
  invisible(path)
}

pa_write_json <- function(x, path, pretty = TRUE, auto_unbox = TRUE) {
  path <- pa_scalar_chr(path, "path")
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not installed; JSON file was not written", call. = FALSE)
    return(invisible(NULL))
  }
  jsonlite::write_json(x, path = path, pretty = pretty, auto_unbox = auto_unbox, null = "null")
  invisible(path)
}

pa_write_rds <- function(x, path) {
  path <- pa_scalar_chr(path, "path")
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path)
  invisible(path)
}

pa_json_read <- function(path, default = NULL) {
  path <- pa_scalar_chr(path, "path")
  if (!file.exists(path)) return(default)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not installed; returning default from pa_json_read()", call. = FALSE)
    return(default)
  }
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) default)
}

pa_read_table_auto <- function(path) {
  path <- pa_scalar_chr(path, "path")
  if (!file.exists(path)) stop(sprintf("Table file does not exist: %s", path), call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "tsv") || identical(ext, "txt")) {
    return(utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (identical(ext, "csv")) {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  utils::read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
}

pa_resolve_python_executable <- function(python_cmd = NULL) {
  candidates <- unique(pa_safe_trim(c(
    python_cmd,
    Sys.getenv("PA_PYTHON_CMD", unset = ""),
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    "/home/h2048/miniconda3/bin/python",
    Sys.which("python3"),
    Sys.which("python")
  )))
  candidates <- candidates[nzchar(candidates)]

  for (cand in candidates) {
    expanded <- path.expand(cand)
    if (file.exists(expanded)) {
      return(normalizePath(expanded, winslash = "/", mustWork = FALSE))
    }
    resolved <- Sys.which(cand)
    if (nzchar(resolved)) {
      return(normalizePath(resolved, winslash = "/", mustWork = FALSE))
    }
  }

  stop(
    "Could not resolve a Python executable. Set PA_PYTHON_CMD or RETICULATE_PYTHON.",
    call. = FALSE
  )
}

pa_run_system_command <- function(command,
                                  args = character(),
                                  wd = NULL,
                                  env = character(),
                                  fail_on_error = TRUE) {
  command <- pa_scalar_chr(command, "command")
  resolved_cmd <- if (file.exists(path.expand(command))) {
    normalizePath(path.expand(command), winslash = "/", mustWork = FALSE)
  } else {
    sys_cmd <- Sys.which(command)
    if (!nzchar(sys_cmd)) command else normalizePath(sys_cmd, winslash = "/", mustWork = FALSE)
  }

  old_wd <- NULL
  if (!is.null(wd)) {
    wd <- pa_scalar_chr(wd, "wd")
    if (!dir.exists(wd)) stop(sprintf("Working directory does not exist: %s", wd), call. = FALSE)
    old_wd <- getwd()
    setwd(wd)
    on.exit(setwd(old_wd), add = TRUE)
  }

  output <- tryCatch(
    system2(resolved_cmd, args = args, stdout = TRUE, stderr = TRUE, env = env),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  result <- list(
    command = resolved_cmd,
    args = args,
    wd = pa_null_coalesce(wd, getwd()),
    status = as.integer(status),
    output = as.character(output)
  )

  if (isTRUE(fail_on_error) && !identical(result$status, 0L)) {
    stop(
      sprintf(
        "Command failed [%s]: %s\n%s",
        result$status,
        paste(c(resolved_cmd, args), collapse = " "),
        paste(result$output, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  result
}

pa_python_module_available <- function(module_name, python_cmd = NULL) {
  module_name <- pa_scalar_chr(module_name, "module_name")
  py_exec <- pa_resolve_python_executable(python_cmd)
  check_code <- sprintf(
    "import importlib.util, sys; sys.stdout.write('1' if importlib.util.find_spec(%s) else '0')",
    dQuote(module_name)
  )
  res <- pa_run_system_command(py_exec, args = c("-c", check_code), fail_on_error = FALSE)
  identical(res$status, 0L) && any(trimws(res$output) == "1")
}

pa_make_unit_id <- function(lineage, state_level, contrast_id, aggregation_type = "cell") {
  paste(
    pa_scalar_chr(lineage, "lineage"),
    pa_scalar_chr(state_level, "state_level"),
    pa_scalar_chr(contrast_id, "contrast_id"),
    pa_scalar_chr(aggregation_type, "aggregation_type"),
    sep = "__"
  )
}

pa_normalize_gene_vector <- function(genes) {
  pa_unique_chr(toupper(genes))
}

pa_new_analysis_unit <- function(
  lineage,
  state_level,
  contrast_id,
  condition_a,
  condition_b,
  sample_scope = "within-lineage",
  matrix_type = "counts",
  aggregation_type = "cell",
  program_sources = c("cNMF", "hdWGCNA", "CoVarNet", "ssGSEA"),
  rewiring_enabled = FALSE,
  trajectory_enabled = FALSE,
  output_dir
) {
  lineage <- pa_scalar_chr(lineage, "lineage")
  state_level <- pa_scalar_chr(state_level, "state_level")
  contrast_id <- pa_scalar_chr(contrast_id, "contrast_id")
  condition_a <- pa_scalar_chr(condition_a, "condition_a")
  condition_b <- pa_scalar_chr(condition_b, "condition_b")
  sample_scope <- pa_scalar_chr(sample_scope, "sample_scope")
  matrix_type <- pa_scalar_chr(matrix_type, "matrix_type")
  aggregation_type <- pa_scalar_chr(aggregation_type, "aggregation_type")
  output_dir <- pa_prepare_output_dir(output_dir)

  data.frame(
    unit_id = pa_make_unit_id(lineage, state_level, contrast_id, aggregation_type),
    lineage = lineage,
    state_level = state_level,
    contrast_id = contrast_id,
    condition_a = condition_a,
    condition_b = condition_b,
    sample_scope = sample_scope,
    matrix_type = matrix_type,
    aggregation_type = aggregation_type,
    program_sources = I(list(pa_unique_chr(program_sources))),
    rewiring_enabled = isTRUE(rewiring_enabled),
    trajectory_enabled = isTRUE(trajectory_enabled),
    output_dir = output_dir,
    helper_version = PA_HELPER_VERSION_20260428_V1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

if (sys.nframe() == 0) {
  cat("Program Architecture Core Helper (2026-04-28 v1) loaded.\n")
}