#!/usr/bin/env Rscript
# ==============================================================================
# Tissue Comparison Advanced Helper
# ==============================================================================

tc_safe_trim <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x[[1]])
}

tc_to_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return("")
  if (length(x) == 1) return(x[[1]])
  paste(x, collapse = "; ")
}

tc_regex_escape <- function(x) {
  gsub("([][{}()+*^$.|\\?])", "\\\\\\1", x)
}

tc_strip_path_prefix <- function(path, prefix) {
  path <- as.character(path)
  prefix <- as.character(prefix)
  sub(sprintf("^%s/?", tc_regex_escape(prefix)), "", path)
}

tc_path_exists <- function(path) {
  !is.null(path) && length(path) == 1 && nzchar(path) && file.exists(path)
}

tc_list_files <- function(dir_path, pattern = NULL, recursive = FALSE, full.names = TRUE) {
  if (!tc_path_exists(dir_path) || !dir.exists(dir_path)) return(character())
  list.files(
    dir_path,
    pattern = pattern,
    recursive = recursive,
    full.names = full.names,
    ignore.case = TRUE
  )
}

tc_safe_read_text <- function(path, max_lines = 200L) {
  if (!tc_path_exists(path)) return(character())
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  if (is.finite(max_lines)) {
    lines <- readLines(con, n = max_lines + 1L, warn = FALSE)
    if (length(lines) > max_lines) {
      return(c(lines[seq_len(max_lines)], sprintf("... [truncated more lines after %d shown]", max_lines)))
    }
    return(lines)
  }
  readLines(con, warn = FALSE)
}

tc_safe_read_table <- function(path) {
  if (!tc_path_exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "tsv") {
    dt_base <- tryCatch(
      utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "\""),
      error = function(e) NULL
    )
    if (!is.null(dt_base)) return(dt_base)
  }
  if (ext == "csv") {
    dt_base <- tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "\""),
      error = function(e) NULL
    )
    if (!is.null(dt_base)) return(dt_base)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    sep <- if (ext == "tsv") "\t" else if (ext == "csv") "," else NULL
    if (!is.null(sep)) {
      dt_sep <- tryCatch(
        data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE, fill = TRUE),
        error = function(e) NULL
      )
      if (!is.null(dt_sep)) return(dt_sep)
    }
    dt_auto <- tryCatch(data.table::fread(path, data.table = FALSE, fill = TRUE), error = function(e) NULL)
    if (!is.null(dt_auto)) return(dt_auto)
    return(NULL)
  }
  tryCatch({
    if (ext == "tsv") {
      utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "\"")
    } else if (ext == "csv") {
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "\"")
    } else {
      utils::read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
    }
  }, error = function(e) NULL)
}

tc_safe_read_rds <- function(path) {
  if (!tc_path_exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

tc_first_existing_path <- function(paths) {
  paths <- as.character(paths)
  paths <- paths[nzchar(paths)]
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

tc_record_to_row <- function(rec) {
  if (is.null(rec) || !is.list(rec)) return(NULL)
  fields <- c(
    "celltype_level", "celltype_label", "celltype_l2", "comparison", "direction",
    "source_db", "status", "warnings", "error", "overview", "key_mechanisms",
    "hypothesis", "narrative", "key_drivers", "evidence", "limitations", "raw_text"
  )
  row <- as.list(setNames(rep("", length(fields)), fields))
  for (nm in intersect(fields, names(rec))) {
    if (identical(nm, "warnings")) {
      row[[nm]] <- paste(as.character(unlist(rec[[nm]], use.names = FALSE)), collapse = " | ")
    } else if (identical(nm, "key_drivers")) {
      row[[nm]] <- paste(as.character(unique(unlist(rec[[nm]], use.names = FALSE))), collapse = ", ")
    } else {
      row[[nm]] <- tc_to_scalar(rec[[nm]])
    }
  }
  as.data.frame(row, stringsAsFactors = FALSE)
}

tc_flatten_structured_records <- function(x) {
  if (is.null(x)) return(data.frame())
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
  rows <- list()
  walk <- function(node) {
    if (is.null(node)) return(invisible(NULL))
    if (is.data.frame(node)) {
      rows[[length(rows) + 1L]] <<- as.data.frame(node, stringsAsFactors = FALSE)
      return(invisible(NULL))
    }
    if (is.list(node)) {
      node_names <- names(node)
      if (!is.null(node_names) && any(node_names %in% c("comparison", "direction", "overview", "status", "celltype_level", "celltype_label"))) {
        row <- tc_record_to_row(node)
        if (!is.null(row)) rows[[length(rows) + 1L]] <<- row
        return(invisible(NULL))
      }
      invisible(lapply(node, walk))
    }
  }
  walk(x)
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

tc_read_previous_artifacts <- function(paths, artifact_type = c("table", "rds")) {
  artifact_type <- match.arg(artifact_type)
  paths <- unique(as.character(paths))
  paths <- paths[nzchar(paths)]
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(if (identical(artifact_type, "table")) data.frame() else NULL)
  }
  if (identical(artifact_type, "table")) {
    flat_files <- existing[tolower(tools::file_ext(existing)) %in% c("tsv", "csv", "txt")]
    if (length(flat_files) > 0) existing <- flat_files
  }
  loaded <- lapply(existing, function(path) {
    ext <- tolower(tools::file_ext(path))
    if (identical(ext, "rds")) {
      obj <- tc_safe_read_rds(path)
      if (identical(artifact_type, "table")) return(tc_flatten_structured_records(obj))
      return(obj)
    }
    if (identical(artifact_type, "table")) return(tc_safe_read_table(path))
    NULL
  })
  loaded <- Filter(function(x) !is.null(x) && (!is.data.frame(x) || nrow(x) > 0), loaded)
  if (length(loaded) == 0) {
    return(if (identical(artifact_type, "table")) data.frame() else NULL)
  }
  if (identical(artifact_type, "rds")) return(loaded[[1]])
  out <- do.call(rbind, lapply(loaded, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

tc_is_placeholder_secret <- function(x) {
  x <- trimws(as.character(x))
  if (!length(x) || !nzchar(x)) return(TRUE)
  lowered <- tolower(x[[1]])
  lowered %in% c(
    "your-deepseek-api-key",
    "your_deepseek_api_key_here",
    "your_deepseek_api_key",
    "your-key",
    "replace_me",
    "changeme"
  ) || grepl("^your[-_a-z]*api[-_a-z]*key", lowered)
}

tc_load_env_file <- function(path, overwrite_placeholder = TRUE) {
  if (!tc_path_exists(path)) return(invisible(FALSE))
  lines <- readLines(path, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- trimws(gsub("^export\\s+", "", sub("=.*$", "", line)))
    val <- trimws(gsub("^['\"]|['\"]$", "", sub("^[^=]*=", "", line)))
    cur <- Sys.getenv(key, unset = "")
    should_set <- !nzchar(cur)
    if (!should_set && isTRUE(overwrite_placeholder)) {
      should_set <- tc_is_placeholder_secret(cur)
    }
    if (nzchar(key) && should_set) {
      Sys.setenv(structure(val, names = key))
    }
  }
  invisible(TRUE)
}

tc_ensure_env_placeholder <- function(path, key) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    if (any(grepl(sprintf("^%s\\s*=", tc_regex_escape(key)), trimws(lines)))) return(invisible(FALSE))
    write(sprintf("%s=YOUR_%s_HERE", key, key), file = path, append = TRUE)
    return(invisible(TRUE))
  }
  writeLines(sprintf("%s=YOUR_%s_HERE", key, key), path)
  invisible(TRUE)
}

tc_path_equal <- function(path_a, path_b) {
  if (is.null(path_a) || is.null(path_b)) return(FALSE)
  if (length(path_a) != 1 || length(path_b) != 1) return(FALSE)
  path_a <- tc_safe_trim(path_a)
  path_b <- tc_safe_trim(path_b)
  if (!nzchar(path_a) || !nzchar(path_b)) return(FALSE)
  norm_a <- normalizePath(path.expand(path_a), winslash = "/", mustWork = FALSE)
  norm_b <- normalizePath(path.expand(path_b), winslash = "/", mustWork = FALSE)
  identical(norm_a, norm_b)
}

tc_normalize_meta_values <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "<NA>"
  x
}

tc_find_ambiguous_sample_tissue_keys <- function(meta_df,
                                                 sample_col,
                                                 tissue_col,
                                                 resolved_cols = character(),
                                                 candidate_cols = character()) {
  compare_cols <- unique(c(resolved_cols, candidate_cols))
  if (is.null(meta_df) || !is.data.frame(meta_df) || nrow(meta_df) == 0 || length(compare_cols) == 0) {
    return(character())
  }
  if (!all(c(sample_col, tissue_col) %in% colnames(meta_df))) return(character())
  compare_cols <- compare_cols[compare_cols %in% colnames(meta_df)]
  if (length(compare_cols) == 0) return(character())
  key_df <- unique(meta_df[, unique(c(sample_col, tissue_col, compare_cols)), drop = FALSE])
  key_df[] <- lapply(key_df, tc_normalize_meta_values)
  resolved_key <- do.call(
    paste,
    c(key_df[, unique(c(sample_col, tissue_col, resolved_cols)), drop = FALSE], sep = "__")
  )
  base_key <- paste(key_df[[sample_col]], key_df[[tissue_col]], sep = "__")
  unique(base_key[duplicated(resolved_key) | duplicated(resolved_key, fromLast = TRUE)])
}

tc_resolve_pseudobulk_group_columns <- function(meta_df,
                                                sample_col,
                                                tissue_col,
                                                candidate_cols = character()) {
  if (is.null(meta_df) || !is.data.frame(meta_df) || nrow(meta_df) == 0) {
    return(list(group_cols = c(sample_col, tissue_col), disambiguation_cols = character(), ambiguous_keys_remaining = character()))
  }
  group_cols <- c(sample_col, tissue_col)
  candidate_cols <- unique(candidate_cols[candidate_cols %in% colnames(meta_df)])
  ambiguous_keys <- tc_find_ambiguous_sample_tissue_keys(
    meta_df,
    sample_col = sample_col,
    tissue_col = tissue_col,
    candidate_cols = candidate_cols
  )
  disambiguation_cols <- character()

  while (length(ambiguous_keys) > 0 && length(candidate_cols) > 0) {
    informative_cols <- candidate_cols[vapply(candidate_cols, function(col) {
      probe_df <- unique(meta_df[, c(sample_col, tissue_col, col), drop = FALSE])
      probe_df[] <- lapply(probe_df, tc_normalize_meta_values)
      probe_key <- paste(probe_df[[sample_col]], probe_df[[tissue_col]], sep = "__")
      probe_df <- probe_df[probe_key %in% ambiguous_keys, , drop = FALSE]
      probe_key <- probe_key[probe_key %in% ambiguous_keys]
      if (nrow(probe_df) == 0) return(FALSE)
      any(vapply(split(probe_df[[col]], probe_key), function(v) length(unique(v)) > 1, logical(1)))
    }, logical(1))]
    if (length(informative_cols) == 0) break
    chosen_col <- informative_cols[[1]]
    group_cols <- c(group_cols, chosen_col)
    disambiguation_cols <- c(disambiguation_cols, chosen_col)
    candidate_cols <- setdiff(candidate_cols, chosen_col)
    ambiguous_keys <- tc_find_ambiguous_sample_tissue_keys(
      meta_df,
      sample_col = sample_col,
      tissue_col = tissue_col,
      resolved_cols = disambiguation_cols,
      candidate_cols = setdiff(candidate_cols, disambiguation_cols)
    )
  }

  list(
    group_cols = unique(group_cols),
    disambiguation_cols = unique(disambiguation_cols),
    ambiguous_keys_remaining = ambiguous_keys
  )
}

tc_write_table_file <- function(df, path, sep = "\t") {
  if (is.null(path) || !length(path) || !nzchar(path)) return(invisible(FALSE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(df, path, sep = sep)
  } else {
    utils::write.table(df, path, sep = sep, row.names = FALSE, quote = TRUE)
  }
  invisible(TRUE)
}

tc_placeholder_noncoding_gene_regex <- function() {
  paste(
    c(
      "^(AC|AL|AP|BX|Z)[0-9]+\\.",
      "^CT[BCD]-",
      "-OT[0-9]+$",
      "^LOC[0-9]+"
    ),
    collapse = "|"
  )
}

tc_default_technical_gene_regex <- function() {
  paste(
    c(
      "^ENSG",
      "^LINC",
      "^MT-",
      "^MT\\.",
      "^MTRNR",
      "^RPS",
      "^RPL",
      "^MRPS",
      "^MRPL",
      "^RP[0-9]+-",
      "^RP[0-9]+$"
    ),
    collapse = "|"
  )
}

tc_noninformative_gene_rule_text <- function() {
  paste(
    "For non-B-cell analyses, do not use IG genes as primary evidence;",
    "for all analyses, treat MT-, ribosomal (RPS/RPL/MRPS/MRPL), ENSG/LINC,",
    "and placeholder non-coding loci (for example AC123456.1, AL123456.1, AP000000.1,",
    "CTB-/CTC-/CTD- clone-style loci, LOC genes, and *-OT transcripts) as non-informative",
    "unless no stronger lineage-relevant evidence exists."
  )
}

tc_is_placeholder_noncoding_gene <- function(genes,
                                             regex = tc_placeholder_noncoding_gene_regex()) {
  genes <- toupper(trimws(as.character(genes)))
  genes[is.na(genes)] <- ""
  nzchar(genes) & grepl(regex, genes, perl = TRUE)
}

tc_is_technical_gene_generic <- function(genes,
                                         technical_regex = tc_default_technical_gene_regex()) {
  genes <- toupper(trimws(as.character(genes)))
  genes[is.na(genes)] <- ""
  nzchar(genes) & grepl(technical_regex, genes, perl = TRUE)
}

tc_is_noninformative_gene <- function(genes,
                                      technical_regex = tc_default_technical_gene_regex(),
                                      placeholder_regex = tc_placeholder_noncoding_gene_regex(),
                                      extra_regex = NULL) {
  flag <- tc_is_technical_gene_generic(genes, technical_regex = technical_regex) |
    tc_is_placeholder_noncoding_gene(genes, regex = placeholder_regex)
  if (!is.null(extra_regex) && length(extra_regex) == 1 && nzchar(tc_safe_trim(extra_regex))) {
    genes_norm <- toupper(trimws(as.character(genes)))
    genes_norm[is.na(genes_norm)] <- ""
    flag <- flag | (nzchar(genes_norm) & grepl(extra_regex, genes_norm, perl = TRUE))
  }
  flag
}

tc_filter_noninformative_gene_symbols <- function(genes,
                                                  technical_regex = tc_default_technical_gene_regex(),
                                                  placeholder_regex = tc_placeholder_noncoding_gene_regex(),
                                                  extra_regex = NULL,
                                                  unique_only = TRUE) {
  genes <- toupper(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes <- genes[!tc_is_noninformative_gene(
    genes,
    technical_regex = technical_regex,
    placeholder_regex = placeholder_regex,
    extra_regex = extra_regex
  )]
  if (isTRUE(unique_only)) unique(genes) else genes
}

tc_build_annotation_count_tables <- function(meta_df,
                                             group_col,
                                             annotation_cols,
                                             group_label = "group") {
  if (is.null(meta_df) || !is.data.frame(meta_df) || nrow(meta_df) == 0) {
    return(list(long = data.frame(), wide = list()))
  }
  if (is.null(group_col) || length(group_col) != 1 || !group_col %in% colnames(meta_df)) {
    return(list(long = data.frame(), wide = list()))
  }
  if (is.null(annotation_cols) || length(annotation_cols) == 0) {
    return(list(long = data.frame(), wide = list()))
  }
  if (is.null(names(annotation_cols)) || any(!nzchar(names(annotation_cols)))) {
    stop("annotation_cols must be a named character vector or list.")
  }

  group_label <- tc_safe_trim(group_label)
  if (!nzchar(group_label)) group_label <- "group"

  meta_use <- meta_df
  meta_use[[group_label]] <- as.character(meta_use[[group_col]])
  keep_group <- !is.na(meta_use[[group_label]]) & nzchar(trimws(meta_use[[group_label]]))
  meta_use <- meta_use[keep_group, , drop = FALSE]
  if (nrow(meta_use) == 0) return(list(long = data.frame(), wide = list()))

  group_sizes <- dplyr::count(meta_use, .data[[group_label]], name = "group_size")
  long_rows <- list()
  wide_tables <- list()

  for (level_name in names(annotation_cols)) {
    annotation_col <- annotation_cols[[level_name]]
    if (is.null(annotation_col) || !annotation_col %in% colnames(meta_use)) next

    count_df <- meta_use |>
      dplyr::mutate(annotation = as.character(.data[[annotation_col]])) |>
      dplyr::filter(!is.na(annotation), nzchar(trimws(annotation))) |>
      dplyr::count(.data[[group_label]], annotation, name = "n_cells") |>
      dplyr::left_join(group_sizes, by = group_label) |>
      dplyr::mutate(
        annotation_level = as.character(level_name),
        pct_in_group = ifelse(group_size > 0, 100 * n_cells / group_size, NA_real_)
      ) |>
      dplyr::arrange(suppressWarnings(as.numeric(.data[[group_label]])), .data[[group_label]], dplyr::desc(n_cells), annotation)
    if (nrow(count_df) == 0) next

    count_df[[group_label]] <- as.character(count_df[[group_label]])
    names(count_df)[names(count_df) == "group_size"] <- paste0(group_label, "_size")
    names(count_df)[names(count_df) == "pct_in_group"] <- paste0("pct_in_", group_label)

    long_rows[[length(long_rows) + 1L]] <- count_df[, c(
      group_label,
      "annotation_level",
      "annotation",
      "n_cells",
      paste0(group_label, "_size"),
      paste0("pct_in_", group_label)
    ), drop = FALSE]

    wide_tables[[level_name]] <- count_df |>
      dplyr::select(dplyr::all_of(group_label), annotation, n_cells) |>
      tidyr::pivot_wider(names_from = annotation, values_from = n_cells, values_fill = 0) |>
      dplyr::left_join(
        group_sizes |> dplyr::rename(!!paste0(group_label, "_size") := group_size),
        by = group_label
      ) |>
      dplyr::relocate(dplyr::all_of(paste0(group_label, "_size")), .after = dplyr::all_of(group_label))
  }

  list(
    long = if (length(long_rows) > 0) dplyr::bind_rows(long_rows) else data.frame(),
    wide = wide_tables
  )
}

tc_build_choir_annotation_count_tables <- function(meta_df,
                                                   choir_col,
                                                   annotation_cols = c(L2 = "cell_type_L2", L3 = "cell_type_L3")) {
  tc_build_annotation_count_tables(
    meta_df = meta_df,
    group_col = choir_col,
    annotation_cols = annotation_cols,
    group_label = "choir_cluster"
  )
}

tc_write_choir_annotation_count_tables <- function(meta_df,
                                                   choir_col,
                                                   choir_dir,
                                                   annotation_cols = c(L2 = "cell_type_L2", L3 = "cell_type_L3"),
                                                   verbose = TRUE) {
  choir_dir <- path.expand(tc_safe_trim(choir_dir))
  if (!nzchar(choir_dir)) {
    return(list(files = list(), tables = list(long = data.frame(), wide = list())))
  }
  dir.create(choir_dir, recursive = TRUE, showWarnings = FALSE)

  tables <- tc_build_choir_annotation_count_tables(
    meta_df = meta_df,
    choir_col = choir_col,
    annotation_cols = annotation_cols
  )

  files <- list(
    l2 = file.path(choir_dir, "choir_cluster_L2_counts.csv"),
    l3 = file.path(choir_dir, "choir_cluster_L3_counts.csv"),
    long = file.path(choir_dir, "choir_cluster_L2_L3_counts_long.csv")
  )

  if (is.data.frame(tables$long) && nrow(tables$long) > 0) {
    tc_write_table_file(tables$long, files$long, sep = ",")
  }
  if (is.list(tables$wide) && !is.null(tables$wide[["L2"]]) && nrow(tables$wide[["L2"]]) > 0) {
    tc_write_table_file(tables$wide[["L2"]], files$l2, sep = ",")
  }
  if (is.list(tables$wide) && !is.null(tables$wide[["L3"]]) && nrow(tables$wide[["L3"]]) > 0) {
    tc_write_table_file(tables$wide[["L3"]], files$l3, sep = ",")
  }

  if (isTRUE(verbose) && is.data.frame(tables$long) && nrow(tables$long) > 0) {
    cat("[OK] CHOIR cluster L2/L3 count tables saved\n")
  }

  list(files = files, tables = tables)
}

tc_pct_label <- function(x, digits = 1L) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x[[1]])) return("NA")
  sprintf(paste0("%.", as.integer(digits), "f%%"), x[[1]])
}

tc_normalize_pct_values <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  finite_vals <- x_num[is.finite(x_num)]
  if (length(finite_vals) > 0 && all(finite_vals >= 0 & finite_vals <= 1.0001)) {
    x_num <- x_num * 100
  }
  x_num
}

tc_pick_deg_pct_cols <- function(df,
                                 pct1_candidates = c(
                                   "pct.1", "pct_1", "pct1", "pct_expr_1", "pct_expr1",
                                   "pct_case", "pct_target", "pct_cluster", "pct_focal", "pct_group_1"
                                 ),
                                 pct2_candidates = c(
                                   "pct.2", "pct_2", "pct2", "pct_expr_2", "pct_expr2",
                                   "pct_reference", "pct_ref", "pct_control", "pct_rest", "pct_group_2"
                                 )) {
  list(
    pct1_col = tc_llm_pick_existing_col(df, pct1_candidates),
    pct2_col = tc_llm_pick_existing_col(df, pct2_candidates)
  )
}

tc_get_seurat_layer_matrix <- function(obj, assay = "RNA", layer = "counts") {
  if (is.null(obj)) return(NULL)
  mat <- NULL
  if (requireNamespace("Seurat", quietly = TRUE)) {
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (is.null(mat) && !is.null(assay)) {
      assay_obj <- tryCatch(obj[[assay]], error = function(e) NULL)
      if (!is.null(assay_obj)) {
        mat <- tryCatch(Seurat::GetAssayData(assay_obj, layer = layer), error = function(e) NULL)
      }
    }
    if (is.null(mat)) {
      slot_name <- if (identical(layer, "counts")) "counts" else "data"
      mat <- tryCatch(
        Seurat::GetAssayData(obj, assay = assay, slot = slot_name),
        error = function(e) NULL
      )
    }
    if (is.null(mat) && !identical(layer, "data")) {
      mat <- tryCatch(
        Seurat::GetAssayData(obj, assay = assay, layer = "data"),
        error = function(e) NULL
      )
    }
  }
  mat
}

tc_compute_deg_expression_pct_lookup <- function(obj,
                                                 genes,
                                                 cells_1,
                                                 cells_2,
                                                 assay = "RNA",
                                                 layer = "counts") {
  if (is.null(obj) || is.null(genes) || length(genes) == 0) return(data.frame())
  mat <- tc_get_seurat_layer_matrix(obj, assay = assay, layer = layer)
  if (is.null(mat)) return(data.frame())

  genes_use <- tc_match_features_to_available(genes, rownames(mat))
  if (length(genes_use) == 0) return(data.frame())

  all_cells <- colnames(mat)
  cells_1 <- intersect(unique(as.character(cells_1)), all_cells)
  cells_2 <- intersect(unique(as.character(cells_2)), all_cells)
  if (length(cells_1) == 0 && length(cells_2) == 0) return(data.frame())

  compute_pct <- function(cell_ids) {
    if (length(cell_ids) == 0) return(rep(NA_real_, length(genes_use)))
    det_mat <- mat[genes_use, cell_ids, drop = FALSE] > 0
    pct <- if (requireNamespace("Matrix", quietly = TRUE) && methods::is(det_mat, "Matrix")) {
      Matrix::rowMeans(det_mat) * 100
    } else {
      rowMeans(as.matrix(det_mat)) * 100
    }
    suppressWarnings(as.numeric(pct))
  }

  data.frame(
    gene = genes_use,
    pct_expr_1 = compute_pct(cells_1),
    pct_expr_2 = compute_pct(cells_2),
    stringsAsFactors = FALSE
  )
}

tc_prepare_deg_pct_columns <- function(df,
                                       gene_col,
                                       pct_lookup = NULL,
                                       pct1_col = NULL,
                                       pct2_col = NULL,
                                       output_pct1_col = "pct_expr_1",
                                       output_pct2_col = "pct_expr_2") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col) || !gene_col %in% colnames(df)) {
    return(df)
  }
  pct_cols <- tc_pick_deg_pct_cols(df)
  if (is.null(pct1_col)) pct1_col <- pct_cols$pct1_col
  if (is.null(pct2_col)) pct2_col <- pct_cols$pct2_col

  out <- df
  out[[output_pct1_col]] <- if (!is.null(pct1_col) && pct1_col %in% colnames(out)) {
    tc_normalize_pct_values(out[[pct1_col]])
  } else {
    rep(NA_real_, nrow(out))
  }
  out[[output_pct2_col]] <- if (!is.null(pct2_col) && pct2_col %in% colnames(out)) {
    tc_normalize_pct_values(out[[pct2_col]])
  } else {
    rep(NA_real_, nrow(out))
  }

  if (is.data.frame(pct_lookup) && nrow(pct_lookup) > 0 && "gene" %in% colnames(pct_lookup)) {
    row_key <- toupper(trimws(as.character(out[[gene_col]])))
    lookup_key <- toupper(trimws(as.character(pct_lookup$gene)))
    idx <- match(row_key, lookup_key)
    if ("pct_expr_1" %in% colnames(pct_lookup)) {
      lookup_pct1 <- tc_normalize_pct_values(pct_lookup$pct_expr_1[idx])
      fill1 <- is.na(out[[output_pct1_col]])
      out[[output_pct1_col]][fill1] <- lookup_pct1[fill1]
    }
    if ("pct_expr_2" %in% colnames(pct_lookup)) {
      lookup_pct2 <- tc_normalize_pct_values(pct_lookup$pct_expr_2[idx])
      fill2 <- is.na(out[[output_pct2_col]])
      out[[output_pct2_col]][fill2] <- lookup_pct2[fill2]
    }
  }
  out
}

tc_format_deg_pct_pair <- function(pct1,
                                   pct2,
                                   pct1_label = "group_1",
                                   pct2_label = "group_2") {
  parts <- character()
  pct1 <- suppressWarnings(as.numeric(pct1))
  pct2 <- suppressWarnings(as.numeric(pct2))
  pct1_label <- tc_safe_trim(pct1_label)
  pct2_label <- tc_safe_trim(pct2_label)
  if (!nzchar(pct1_label)) pct1_label <- "group_1"
  if (!nzchar(pct2_label)) pct2_label <- "group_2"
  if (!is.na(pct1)) parts <- c(parts, sprintf("pct[%s]=%s", pct1_label, tc_pct_label(pct1)))
  if (!is.na(pct2)) parts <- c(parts, sprintf("pct[%s]=%s", pct2_label, tc_pct_label(pct2)))
  if (!is.na(pct1) && !is.na(pct2)) {
    parts <- c(parts, sprintf("Δpct=%s", tc_pct_label(pct1 - pct2)))
  }
  if (length(parts) == 0) return("")
  paste0(", ", paste(parts, collapse = ", "))
}

tc_choose_alignment_degree <- function(primary_pct, secondary_pct = NA_real_) {
  primary_pct <- suppressWarnings(as.numeric(primary_pct))
  secondary_pct <- suppressWarnings(as.numeric(secondary_pct))
  if (length(primary_pct) == 0 || is.na(primary_pct[[1]])) return("unknown")
  primary_pct <- primary_pct[[1]]
  gap <- primary_pct - ifelse(is.na(secondary_pct[[1]]), 0, secondary_pct[[1]])
  if (primary_pct >= 95) return("very_high")
  if (primary_pct >= 85 && gap >= 50) return("high")
  if (primary_pct >= 70 && gap >= 25) return("moderate")
  if (primary_pct >= 50) return("mixed")
  "low"
}

tc_extract_cluster_annotation_profile <- function(annotation_long_df,
                                                  cluster_id,
                                                  annotation_level = "L3",
                                                  top_n = 5L) {
  if (is.null(annotation_long_df) || !is.data.frame(annotation_long_df) || nrow(annotation_long_df) == 0) {
    return(data.frame())
  }
  required_cols <- c("choir_cluster", "annotation_level", "annotation", "n_cells")
  if (!all(required_cols %in% colnames(annotation_long_df))) return(data.frame())
  df <- annotation_long_df
  df$choir_cluster <- as.character(df$choir_cluster)
  df$annotation_level <- as.character(df$annotation_level)
  df <- df[
    !is.na(df$choir_cluster) & df$choir_cluster == as.character(cluster_id) &
      !is.na(df$annotation_level) & toupper(df$annotation_level) == toupper(annotation_level),
    , drop = FALSE
  ]
  if (nrow(df) == 0) return(data.frame())
  if (!"pct_in_choir_cluster" %in% colnames(df) && all(c("n_cells", "choir_cluster_size") %in% colnames(df))) {
    df$pct_in_choir_cluster <- ifelse(df$choir_cluster_size > 0, 100 * df$n_cells / df$choir_cluster_size, NA_real_)
  }
  df <- df[order(-suppressWarnings(as.numeric(df$n_cells)), df$annotation), , drop = FALSE]
  utils::head(df, max(1L, as.integer(top_n)))
}

tc_summarize_cluster_annotation_profile <- function(annotation_long_df,
                                                    cluster_id,
                                                    annotation_level = "L3",
                                                    top_n = 5L) {
  df <- tc_extract_cluster_annotation_profile(
    annotation_long_df = annotation_long_df,
    cluster_id = cluster_id,
    annotation_level = annotation_level,
    top_n = top_n
  )
  if (nrow(df) == 0) {
    return(list(
      cluster_id = as.character(cluster_id),
      annotation_level = annotation_level,
      cluster_size = NA_real_,
      dominant_label = NA_character_,
      dominant_n = NA_real_,
      dominant_pct = NA_real_,
      secondary_label = NA_character_,
      secondary_n = NA_real_,
      secondary_pct = NA_real_,
      purity_class = "unknown",
      top_table = data.frame(),
      top_labels_text = "NA",
      summary_text = sprintf("No %s annotation-count evidence available for CHOIR cluster %s.", annotation_level, cluster_id)
    ))
  }
  dominant <- df[1, , drop = FALSE]
  secondary <- if (nrow(df) >= 2) df[2, , drop = FALSE] else NULL
  cluster_size <- if ("choir_cluster_size" %in% colnames(df)) suppressWarnings(as.numeric(df$choir_cluster_size[1])) else NA_real_
  if (is.na(cluster_size) && "cluster_size" %in% colnames(df)) cluster_size <- suppressWarnings(as.numeric(df$cluster_size[1]))
  top_labels_text <- paste(vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      "%s=%s/%s (%s)",
      as.character(df$annotation[i]),
      as.character(df$n_cells[i]),
      ifelse(is.na(cluster_size), "NA", as.character(as.integer(cluster_size))),
      tc_pct_label(df$pct_in_choir_cluster[i])
    )
  }, character(1)), collapse = "; ")
  purity_class <- tc_choose_alignment_degree(
    primary_pct = dominant$pct_in_choir_cluster[1],
    secondary_pct = if (is.null(secondary)) NA_real_ else secondary$pct_in_choir_cluster[1]
  )
  summary_text <- sprintf(
    paste(
      "Annotated %s composition for CHOIR cluster %s (n=%s):",
      "dominant=%s (%s), secondary=%s (%s), purity_class=%s.",
      "Top labels: %s"
    ),
    annotation_level,
    as.character(cluster_id),
    ifelse(is.na(cluster_size), "NA", as.character(as.integer(cluster_size))),
    as.character(dominant$annotation[1]),
    tc_pct_label(dominant$pct_in_choir_cluster[1]),
    if (is.null(secondary)) "NA" else as.character(secondary$annotation[1]),
    if (is.null(secondary)) "NA" else tc_pct_label(secondary$pct_in_choir_cluster[1]),
    purity_class,
    top_labels_text
  )
  list(
    cluster_id = as.character(cluster_id),
    annotation_level = annotation_level,
    cluster_size = cluster_size,
    dominant_label = as.character(dominant$annotation[1]),
    dominant_n = suppressWarnings(as.numeric(dominant$n_cells[1])),
    dominant_pct = suppressWarnings(as.numeric(dominant$pct_in_choir_cluster[1])),
    secondary_label = if (is.null(secondary)) NA_character_ else as.character(secondary$annotation[1]),
    secondary_n = if (is.null(secondary)) NA_real_ else suppressWarnings(as.numeric(secondary$n_cells[1])),
    secondary_pct = if (is.null(secondary)) NA_real_ else suppressWarnings(as.numeric(secondary$pct_in_choir_cluster[1])),
    purity_class = purity_class,
    top_table = df,
    top_labels_text = top_labels_text,
    summary_text = summary_text
  )
}

tc_build_choir_annotation_evidence <- function(annotation_tables,
                                               cluster_id,
                                               top_n_l2 = 3L,
                                               top_n_l3 = 5L) {
  annotation_long_df <- NULL
  if (is.list(annotation_tables) && !is.null(annotation_tables$long) && is.data.frame(annotation_tables$long)) {
    annotation_long_df <- annotation_tables$long
  } else if (is.data.frame(annotation_tables)) {
    annotation_long_df <- annotation_tables
  }
  if (is.null(annotation_long_df) || nrow(annotation_long_df) == 0) {
    empty <- list(
      l2 = tc_summarize_cluster_annotation_profile(data.frame(), cluster_id, annotation_level = "L2", top_n = top_n_l2),
      l3 = tc_summarize_cluster_annotation_profile(data.frame(), cluster_id, annotation_level = "L3", top_n = top_n_l3),
      text = sprintf("No CHOIR annotation-count table available for cluster %s.", cluster_id)
    )
    return(empty)
  }
  l2_summary <- tc_summarize_cluster_annotation_profile(annotation_long_df, cluster_id, annotation_level = "L2", top_n = top_n_l2)
  l3_summary <- tc_summarize_cluster_annotation_profile(annotation_long_df, cluster_id, annotation_level = "L3", top_n = top_n_l3)
  list(
    l2 = l2_summary,
    l3 = l3_summary,
    text = paste(
      "CHOIR annotation-count evidence (derived from user L2/L3 labels):",
      l2_summary$summary_text,
      l3_summary$summary_text,
      sep = "\n"
    )
  )
}

tc_lookup_in_caller <- function(name, default = NULL, caller_env = parent.frame()) {
  if (!is.character(name) || length(name) != 1 || !nzchar(tc_safe_trim(name))) return(default)
  if (exists(name, envir = caller_env, inherits = TRUE)) {
    return(get(name, envir = caller_env, inherits = TRUE))
  }
  default
}

tc_llm_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

tc_placeholder_text <- function(x, default = "Not available from current evidence.") {
  x <- tc_safe_trim(x)
  if (nzchar(x)) x else default
}

tc_truncate_text <- function(x, max_chars = 12000L) {
  x <- tc_safe_trim(x)
  if (!nzchar(x)) return("")
  if (nchar(x, type = "chars") <= max_chars) return(x)
  paste0(substr(x, 1, max_chars), "\n...[truncated]")
}

tc_append_nonempty_blocks <- function(..., sep = "\n\n") {
  blocks <- unlist(list(...), use.names = FALSE)
  blocks <- blocks[!is.na(blocks) & nzchar(trimws(blocks))]
  paste(blocks, collapse = sep)
}

tc_collapse_driver_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  paste(unique(x), collapse = ", ")
}

tc_collapse_source_db_label <- function(db_names, representative_db = NULL) {
  db_names <- unique(as.character(db_names))
  db_names <- db_names[!is.na(db_names) & nzchar(trimws(db_names))]
  if (length(db_names) == 0) {
    if (!is.null(representative_db) && nzchar(tc_safe_trim(representative_db))) return(representative_db)
    return(NA_character_)
  }
  label <- paste(db_names, collapse = " + ")
  if (!is.null(representative_db) && nzchar(tc_safe_trim(representative_db))) {
    return(sprintf("multi_db[%s] | seed=%s", label, representative_db))
  }
  sprintf("multi_db[%s]", label)
}

tc_llm_pick_existing_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

tc_format_top_deg_entries <- function(df,
                                      gene_col,
                                      logfc_col = NULL,
                                      padj_col = NULL,
                                      top_n = 10L,
                                      decreasing = TRUE,
                                      pct_lookup = NULL,
                                      pct1_col = NULL,
                                      pct2_col = NULL,
                                      pct1_label = "group_1",
                                      pct2_label = "group_2") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col)) return("None")
  df <- tc_prepare_deg_pct_columns(
    df,
    gene_col = gene_col,
    pct_lookup = pct_lookup,
    pct1_col = pct1_col,
    pct2_col = pct2_col
  )
  lfc <- if (!is.null(logfc_col)) tc_llm_num(df[[logfc_col]]) else rep(NA_real_, nrow(df))
  padj <- if (!is.null(padj_col)) tc_llm_num(df[[padj_col]]) else rep(Inf, nrow(df))
  ord <- if (!all(is.na(lfc))) {
    order(if (decreasing) -lfc else lfc, padj, na.last = TRUE)
  } else {
    order(padj, na.last = TRUE)
  }
  df <- df[ord, , drop = FALSE]
  df <- utils::head(df, max(1L, as.integer(top_n)))
  genes <- as.character(df[[gene_col]])
  genes[is.na(genes) | !nzchar(trimws(genes))] <- "NA"
  if (!is.null(logfc_col) && !is.null(padj_col)) {
    return(paste0(
      genes, "(", round(tc_llm_num(df[[logfc_col]]), 2), ", padj=",
      signif(tc_llm_num(df[[padj_col]]), 2), ")",
      collapse = ", "
    ))
  }
  if (!is.null(logfc_col)) {
    return(paste0(genes, "(", round(tc_llm_num(df[[logfc_col]]), 2), ")", collapse = ", "))
  }
  paste(genes, collapse = ", ")
}

tc_format_integrated_deg_entries <- function(df,
                                             gene_col,
                                             logfc_col,
                                             padj_col = NULL,
                                             top_n = 10L) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col) || is.null(logfc_col)) return("None")
  lfc <- tc_llm_num(df[[logfc_col]])
  padj <- if (!is.null(padj_col)) tc_llm_num(df[[padj_col]]) else rep(Inf, nrow(df))
  ord <- order(-abs(lfc), padj, na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  df <- utils::head(df, max(1L, as.integer(top_n)))
  genes <- as.character(df[[gene_col]])
  genes[is.na(genes) | !nzchar(trimws(genes))] <- "NA"
  paste0(
    genes,
    "(log2FC=", round(tc_llm_num(df[[logfc_col]]), 2),
    if (!is.null(padj_col)) paste0(", padj=", signif(tc_llm_num(df[[padj_col]]), 2)) else "",
    ")",
    collapse = ", "
  )
}

tc_build_integrated_directional_enrichment_text <- function(evidence_df,
                                                            direction_map = NULL,
                                                            caller_env = parent.frame()) {
  if (is.null(evidence_df) || !is.data.frame(evidence_df) || nrow(evidence_df) == 0) return("")
  if (!"direction" %in% colnames(evidence_df)) return("")
  sanitize_enrichment_df <- tc_lookup_in_caller(
    "tc_sanitize_enrichment_df",
    function(df, pathway_col = "Description", geneid_col = "geneID") {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)
      if (!is.null(pathway_col) && pathway_col %in% colnames(df)) {
        df[[pathway_col]] <- trimws(as.character(df[[pathway_col]]))
      }
      if (!is.null(geneid_col) && geneid_col %in% colnames(df)) {
        df[[geneid_col]] <- as.character(df[[geneid_col]])
        df[[geneid_col]][is.na(df[[geneid_col]])] <- ""
      }
      df
    },
    caller_env
  )
  df <- evidence_df
  if (!"db" %in% colnames(df)) df$db <- "NA"
  if (!"Description" %in% colnames(df)) df$Description <- rownames(df)
  if (!"geneID" %in% colnames(df)) df$geneID <- ""
  if (!"Count" %in% colnames(df)) df$Count <- NA_integer_
  pct_text <- vapply(seq_len(nrow(df)), function(i) {
    tc_format_deg_pct_pair(
      df$pct_expr_1[i],
      df$pct_expr_2[i],
      pct1_label = pct1_label,
      pct2_label = pct2_label
    )
  }, character(1))
  if (!"p.adjust" %in% colnames(df)) df$p.adjust <- NA_real_
  df <- sanitize_enrichment_df(df)
      genes, "(log2FC=", round(tc_llm_num(df[[logfc_col]]), 2), ", padj=",
      signif(tc_llm_num(df[[padj_col]]), 2), pct_text, ")",
  lines <- vapply(seq_len(nrow(df)), function(i) {
    dir_key <- df$direction[i]
    dir_label <- if (!is.null(direction_map) && dir_key %in% names(direction_map)) direction_map[[dir_key]] else toupper(dir_key)
    sprintf(
    return(paste0(
      genes,
      "(log2FC=",
      round(tc_llm_num(df[[logfc_col]]), 2),
      pct_text,
      ")",
      collapse = ", "
    ))
      dir_label,
      tc_safe_trim(df$db[i]),
      tc_safe_trim(gsub("^\\[[^]]+\\]\\s*", "", tc_safe_trim(df$Description[i]))),
      format(df$p.adjust[i], scientific = TRUE, digits = 3),
      tc_safe_trim(df$Count[i]),
      tc_truncate_text(gsub("/", ", ", tc_safe_trim(df$geneID[i])), 220)
    )
  }, character(1))
                                             top_n = 10L,
                                             pct_lookup = NULL,
                                             pct1_col = NULL,
                                             pct2_col = NULL,
                                             pct1_label = "group_1",
                                             pct2_label = "group_2") {
    "Integrated cross-database enrichment summary (all comparison directions are provided together; synthesize them in one judgment):",
  df <- tc_prepare_deg_pct_columns(
    df,
    gene_col = gene_col,
    pct_lookup = pct_lookup,
    pct1_col = pct1_col,
    pct2_col = pct2_col
  )
    lines
  ), collapse = "\n")
}

tc_build_choir_deg_evidence_text <- function(de_df,
                                             top_n_abs = NULL,
                                             top_n_up = 12L,
  pct_text <- vapply(seq_len(nrow(df)), function(i) {
    tc_format_deg_pct_pair(
      df$pct_expr_1[i],
      df$pct_expr_2[i],
      pct1_label = pct1_label,
      pct2_label = pct2_label
    )
  }, character(1))
                                             top_n_down = 12L,
                                             padj_thr = NULL,
                                             lfc_thr = NULL,
                                             caller_env = parent.frame()) {
    pct_text,
  if (is.null(top_n_abs)) top_n_abs <- tc_lookup_in_caller("CHOIR_LLM_TOP_DEG_N", 20L, caller_env)
  if (is.null(padj_thr)) padj_thr <- tc_lookup_in_caller("OFA_PADJ_THR", 0.05, caller_env)
  if (is.null(lfc_thr)) lfc_thr <- tc_lookup_in_caller("OFA_LFC_THR", 0.25, caller_env)
  if (is.null(de_df) || !is.data.frame(de_df) || nrow(de_df) == 0) return("")

  pick_col <- tc_lookup_in_caller("llm_pick_existing_col", tc_llm_pick_existing_col, caller_env)
                                             cluster_id = NULL,
  filter_df_by_gene_col <- tc_lookup_in_caller(
    "tc_filter_df_by_gene_col",
    function(df, gene_col) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || is.null(gene_col) || !gene_col %in% colnames(df)) return(df)
                                             lfc_thr = NULL,
                                             pct1_label = "cluster",
                                             pct2_label = "rest",
      df[keep, , drop = FALSE]
    },
    caller_env
  )

  gene_col <- pick_col(de_df, c("gene", "symbol", "feature", "features", "genes"))
  logfc_col <- pick_col(de_df, c("avg_log2FC", "log2FoldChange", "log2FC", "avg_logFC"))
  padj_col <- pick_col(de_df, c("padj", "p_val_adj", "FDR", "qvalue", "p_adj"))
  if (is.null(gene_col) || is.null(logfc_col)) return("")
  df <- filter_df_by_gene_col(de_df, gene_col)
  if (nrow(df) == 0) return("")

  lfc <- tc_llm_num(df[[logfc_col]])
  padj <- if (!is.null(padj_col)) tc_llm_num(df[[padj_col]]) else rep(NA_real_, nrow(df))
  sig_flag <- !is.na(lfc) & abs(lfc) >= lfc_thr
  if (!is.null(padj_col)) sig_flag <- sig_flag & !is.na(padj) & padj <= padj_thr
  sig_df <- df[sig_flag, , drop = FALSE]
  up_df <- sig_df[tc_llm_num(sig_df[[logfc_col]]) > 0, , drop = FALSE]
  down_df <- sig_df[tc_llm_num(sig_df[[logfc_col]]) < 0, , drop = FALSE]
  abs_df <- if (nrow(sig_df) > 0) sig_df else df[!is.na(lfc), , drop = FALSE]
  pct_cols <- tc_pick_deg_pct_cols(df)
  pct_lookup <- NULL
  if ((is.null(pct_cols$pct1_col) || is.null(pct_cols$pct2_col)) && !is.null(cluster_id) && nzchar(tc_safe_trim(cluster_id))) {
    obj <- tc_lookup_in_caller("obj", NULL, caller_env)
    choir_col <- tc_lookup_in_caller("choir_col", NULL, caller_env)
    if (!is.null(obj) && !is.null(choir_col) && choir_col %in% colnames(obj@meta.data)) {
      meta <- obj@meta.data
      focal_cells <- rownames(meta)[as.character(meta[[choir_col]]) == as.character(cluster_id)]
      rest_cells <- rownames(meta)[as.character(meta[[choir_col]]) != as.character(cluster_id)]
      genes_for_pct <- unique(c(
        as.character(abs_df[[gene_col]]),
        as.character(up_df[[gene_col]]),
        as.character(down_df[[gene_col]])
      ))
      pct_lookup <- tc_compute_deg_expression_pct_lookup(
        obj = obj,
        genes = genes_for_pct,
        cells_1 = focal_cells,
        cells_2 = rest_cells
      )
    }
  }

  paste(
    sprintf(
      "Complete DEG evidence for CHOIR cluster vs rest (sig genes: up=%d, down=%d, total=%d; thresholds: padj<=%s and |log2FC|>=%s; expression proportions shown as pct[%s] and pct[%s]).",
      nrow(up_df), nrow(down_df), nrow(sig_df), padj_thr, lfc_thr, pct1_label, pct2_label
    ),
    paste0(
      "- Top absolute-effect DEG: ",
      tc_format_integrated_deg_entries(
        abs_df,
        gene_col,
        logfc_col,
        padj_col,
        top_n_abs,
        pct_lookup = pct_lookup,
        pct1_col = pct_cols$pct1_col,
        pct2_col = pct_cols$pct2_col,
        pct1_label = pct1_label,
        pct2_label = pct2_label
      )
    ),
    paste0(
      "- Top cluster-enriched DEG (up): ",
      tc_format_top_deg_entries(
        up_df,
        gene_col,
        logfc_col,
        padj_col,
        top_n_up,
        decreasing = TRUE,
        pct_lookup = pct_lookup,
        pct1_col = pct_cols$pct1_col,
        pct2_col = pct_cols$pct2_col,
        pct1_label = pct1_label,
        pct2_label = pct2_label
      )
    ),
    paste0(
      "- Top rest-enriched DEG (down): ",
      tc_format_top_deg_entries(
        down_df,
        gene_col,
        logfc_col,
        padj_col,
        top_n_down,
        decreasing = FALSE,
        pct_lookup = pct_lookup,
        pct1_col = pct_cols$pct1_col,
        pct2_col = pct_cols$pct2_col,
        pct1_label = pct1_label,
        pct2_label = pct2_label
      )
    ),
    sep = "\n"
  )
}

tc_build_choir_cluster_ssgsea_evidence <- function(ssgsea_choir_all,
                                                   cluster_id,
                                                   n_top_per_method = NULL,
                                                   caller_env = parent.frame()) {
  if (is.null(n_top_per_method)) {
    n_top_per_method <- tc_lookup_in_caller("CHOIR_LLM_TOP_SSGSEA_TERMS_PER_DIRECTION", 8L, caller_env)
  }
  if (is.null(ssgsea_choir_all) || length(ssgsea_choir_all) == 0) return("")
  top_df_all <- dplyr::bind_rows(lapply(ssgsea_choir_all, function(x) x[["top_df"]]))
  if (nrow(top_df_all) == 0 || !"cluster" %in% colnames(top_df_all)) return("")
  top_df_cluster <- top_df_all[as.character(top_df_all$cluster) == as.character(cluster_id), , drop = FALSE]
  if (nrow(top_df_cluster) == 0) return("")

  build_ssgsea_evidence <- tc_lookup_in_caller("build_ssgsea_evidence_text", NULL, caller_env)
  evidence_body <- if (is.function(build_ssgsea_evidence)) {
    build_ssgsea_evidence(top_df_cluster, n_top_per_method = n_top_per_method)
  } else {
    if (!"direction" %in% colnames(top_df_cluster)) {
      top_df_cluster$direction <- ifelse(top_df_cluster$z_score < 0, "negative", "positive")
    }
    top_df_cluster$direction <- ifelse(
      tolower(as.character(top_df_cluster$direction)) %in% c("negative", "down", "neg", "suppressed"),
      "negative",
      "positive"
    )
    format_lines <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      label <- if ("method" %in% colnames(df) && length(unique(as.character(df$method))) > 1) {
        sprintf("[%s] %s", as.character(df$method), tc_safe_trim(df$pathway))
      } else {
        tc_safe_trim(df$pathway)
      }
      vapply(seq_len(nrow(df)), function(i) {
        sprintf("- %s | z_score=%.3f | score=%.3f | rank=%d",
                label[i], df$z_score[i], df$score[i], df$rank[i])
      }, character(1))
    }
    pos_df <- top_df_cluster[top_df_cluster$direction == "positive", , drop = FALSE]
    neg_df <- top_df_cluster[top_df_cluster$direction == "negative", , drop = FALSE]
    pos_df <- pos_df[order(pos_df$z_score, pos_df$score, decreasing = TRUE), , drop = FALSE]
    neg_df <- neg_df[order(neg_df$z_score, neg_df$score, decreasing = FALSE), , drop = FALSE]
    pos_df <- utils::head(pos_df, max(1L, as.integer(n_top_per_method)))
    neg_df <- utils::head(neg_df, max(1L, as.integer(n_top_per_method)))
    paste(c(
      "Integrated ssGSEA pathway summary across all selected gene set databases (do not interpret database-by-database; integrate positive and negative signals together):",
      if (nrow(pos_df) > 0) c("Integrated positive / identity-supporting pathways:", format_lines(pos_df)) else NULL,
      if (nrow(neg_df) > 0) c("Integrated negative / suppressed pathways:", format_lines(neg_df)) else NULL
    ), collapse = "\n")
  }

  paste(
    "Cluster-level ssGSEA evidence (use mainly for cell-type/state judgment):",
    evidence_body,
    sep = "\n"
  )
}

tc_merge_directional_gene_pathway_maps <- function(bundle_list) {
  rows <- lapply(names(bundle_list), function(dir_name) {
    bundle <- bundle_list[[dir_name]]
    if (is.null(bundle) || is.null(bundle$gene_pathway_map) || nrow(bundle$gene_pathway_map) == 0) return(NULL)
    df <- bundle$gene_pathway_map
    df$pathway <- sprintf("[%s] %s", toupper(dir_name), df$pathway)
    df
  })
  dplyr::bind_rows(rows)
}

tc_run_choir_cluster_review_llm <- function(cluster_id,
                                            choir_ctx,
                                            evidence_text,
                                            source_db_label,
                                            warnings = character(),
                                            error_message = NULL,
                                            caller_env = parent.frame()) {
  enable_llm <- isTRUE(tc_lookup_in_caller("ENABLE_LLM", FALSE, caller_env))
  if (!enable_llm) {
    return(list(result = NULL, warnings = character(), error = "LLM disabled"))
  }

  lineage_context_lower <- tc_lookup_in_caller("LINEAGE_CONTEXT_LOWER", "tissue", caller_env)
  standardize_model <- tc_lookup_in_caller("STANDARDIZE_LLM_MODEL", "deepseek-chat", caller_env)
  deepseek_api_key <- tc_lookup_in_caller("DEEPSEEK_API_KEY", Sys.getenv("DEEPSEEK_API_KEY", unset = ""), caller_env)
  standardize_retries <- as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_RETRIES", 3L, caller_env))
  standardize_retry_sleep <- tc_lookup_in_caller("STANDARDIZE_LLM_RETRY_SLEEP_SEC", 2, caller_env)
  max_input_chars <- max(as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_INPUT_CHARS", 12000L, caller_env)), 24000L)

  prompt <- paste(
    sprintf("You are reviewing a CHOIR one-vs-rest cluster in %s tissue comparison.", lineage_context_lower),
    "Return valid JSON only. No markdown, no code fences, no commentary.",
    "Use exactly these keys:",
    "cell_type_judgment, confidence, annotation_match_degree, annotated_l3_correspondence, outlier_assessment, discovery_assessment, integrated_diagnostic_comment, overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "",
    "Rules:",
    "1. All narrative fields must be Simplified Chinese strings; key_drivers must be an array of English gene symbols.",
    "2. confidence must be one of: high, medium, low.",
    "3. annotation_match_degree must be one of: high, moderate, mixed, low.",
    "4. annotated_l3_correspondence must explicitly compare the inferred cell type/state with the dominant user L3 annotation and mention the quantitative annotation composition when relevant.",
    "5. outlier_assessment must state whether this cluster is more consistent with annotation error / mislabel / strong ambient contamination / stress artifact; if not, say that current evidence does not support annotation error.",
    "6. discovery_assessment must state whether this cluster is more consistent with a plausible biological finding or substate; if not, say so explicitly.",
    "7. integrated_diagnostic_comment must synthesize annotation composition, complete DEG, ssGSEA, and OFA enrichment together in 2-4 Chinese sentences, and explain whether the main interpretation is aligned annotation, mixed/transition state, outlier/misannotation, tissue contamination, stress program, or plausible biological discovery.",
    "8. overview must start with a concise cell-type/state judgment, then summarize the cluster in 1-2 sentences.",
    "9. key_mechanisms should prioritize lineage-informative and state-informative genes/pathways, not housekeeping signals.",
    "10. hypothesis should be brief; use 'Not applicable.' if no strong comparative or discovery-oriented hypothesis is justified.",
    "11. evidence must explicitly cite DEG, ssGSEA, OFA, and annotation-count clues when they are present.",
    sprintf("12. %s", tc_noninformative_gene_rule_text()),
    "13. If annotation purity is low or the molecular evidence strongly conflicts with the dominant user L3 label, downgrade annotation_match_degree and explain why.",
    "14. If annotation purity is high but the cluster shows a coherent unexpected program that is not well explained by contamination or technical noise, consider it as a possible biological discovery rather than an annotation error.",
    "15. Be conservative: do not call a discovery if the pattern is better explained by stress, plasma spillover, epithelial contamination, or other ambient RNA.",
    "",
    sprintf("Cluster: CHOIR_%s", cluster_id),
    sprintf("Source DB summary: %s", source_db_label),
    "Context:",
    tc_truncate_text(choir_ctx, max_chars = max_input_chars),
    "Evidence:",
    tc_truncate_text(evidence_text, max_chars = max_input_chars),
    "Warnings and errors:",
    tc_truncate_text(paste(c(warnings, error_message), collapse = "\n"), max_chars = 3000L),
    sep = "\n"
  )

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(max(1L, standardize_retries))) {
    if (attempt > 1) cat(sprintf("    [INFO] CHOIR review LLM retry %d/%d\n", attempt, standardize_retries))
    response_text <- tryCatch(
      tc_safe_trim(fanyi::chat_request(prompt, model = standardize_model, api_key = deepseek_api_key)),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )
    parsed <- tc_parse_json(response_text)
    if (!is.null(parsed) && is.list(parsed)) {
      return(list(result = parsed, warnings = unique(collected_warnings), error = NULL))
    }
    if (nzchar(response_text)) {
      collected_warnings <- c(collected_warnings, sprintf("choir review attempt %d returned non-JSON text", attempt))
    }
    if (attempt < standardize_retries) Sys.sleep(standardize_retry_sleep)
  }

  list(
    result = NULL,
    warnings = unique(collected_warnings),
    error = ifelse(is.null(last_error) || !nzchar(last_error), "CHOIR review LLM failed", last_error)
  )
}

tc_run_ssgsea_group_review_llm <- function(group_id,
                                           ss_ctx,
                                           evidence_text,
                                           source_db_label,
                                           annotated_label,
                                           annotated_level,
                                           warnings = character(),
                                           error_message = NULL,
                                           caller_env = parent.frame()) {
  enable_llm <- isTRUE(tc_lookup_in_caller("ENABLE_LLM", FALSE, caller_env))
  if (!enable_llm) {
    return(list(result = NULL, warnings = character(), error = "LLM disabled"))
  }

  lineage_context_lower <- tc_lookup_in_caller("LINEAGE_CONTEXT_LOWER", "tissue", caller_env)
  standardize_model <- tc_lookup_in_caller("STANDARDIZE_LLM_MODEL", "deepseek-chat", caller_env)
  deepseek_api_key <- tc_lookup_in_caller("DEEPSEEK_API_KEY", Sys.getenv("DEEPSEEK_API_KEY", unset = ""), caller_env)
  standardize_retries <- as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_RETRIES", 3L, caller_env))
  standardize_retry_sleep <- tc_lookup_in_caller("STANDARDIZE_LLM_RETRY_SLEEP_SEC", 2, caller_env)
  max_input_chars <- max(as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_INPUT_CHARS", 12000L, caller_env)), 18000L)
  parse_json <- tc_lookup_in_caller("parse_standardized_json", tc_parse_json, caller_env)

  prompt <- paste(
    sprintf("You are reviewing a grouped ssGSEA profile in %s tissue comparison.", lineage_context_lower),
    "Return valid JSON only. No markdown, no code fences, no commentary.",
    "Use exactly these keys:",
    "cell_type_judgment, confidence, annotation_match_degree, annotated_l3_correspondence, outlier_assessment, discovery_assessment, integrated_diagnostic_comment, overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "",
    "Rules:",
    "1. All narrative fields must be Simplified Chinese strings; key_drivers must be an array of English gene symbols or marker names if needed.",
    "2. confidence must be one of: high, medium, low.",
    "3. annotation_match_degree must be one of: high, moderate, mixed, low.",
    "4. annotated_l3_correspondence must explicitly compare the ssGSEA-inferred state with the provided user annotation label.",
    "5. outlier_assessment must state whether the ssGSEA pattern suggests likely misannotation, severe contamination, or technical artifact; if not, say current evidence does not support that.",
    "6. discovery_assessment must state whether the ssGSEA pattern suggests a plausible biological substate / tissue-adapted program / activation state worth follow-up; if not, say so explicitly.",
    "7. integrated_diagnostic_comment must integrate positive and negative ssGSEA signals together and explain whether this group is annotation-aligned, mixed/transition-like, likely outlier, or plausible biological discovery.",
    "8. Use ssGSEA as the main evidence; do not overstate mechanism beyond the pathways provided.",
    "9. overview must start with a concise cell-type/state judgment.",
    "10. hypothesis should be brief; use 'Not applicable.' if no strong hypothesis is justified.",
    sprintf("11. %s", tc_noninformative_gene_rule_text()),
    "12. Be conservative: broad stress/metabolic signatures alone do not automatically imply discovery.",
    "",
    sprintf("Group: %s", group_id),
    sprintf("Annotated label: %s (%s)", annotated_label, annotated_level),
    sprintf("Source DB summary: %s", source_db_label),
    "Context:",
    tc_truncate_text(ss_ctx, max_chars = max_input_chars),
    "Evidence:",
    tc_truncate_text(evidence_text, max_chars = max_input_chars),
    "Warnings and errors:",
    tc_truncate_text(paste(c(warnings, error_message), collapse = "\n"), max_chars = 3000L),
    sep = "\n"
  )

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(max(1L, standardize_retries))) {
    if (attempt > 1) cat(sprintf("    [INFO] ssGSEA review LLM retry %d/%d\n", attempt, standardize_retries))
    response_text <- tryCatch(
      tc_safe_trim(fanyi::chat_request(prompt, model = standardize_model, api_key = deepseek_api_key)),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )
    parsed <- parse_json(response_text)
    if (!is.null(parsed) && is.list(parsed)) {
      return(list(result = parsed, warnings = unique(collected_warnings), error = NULL))
    }
    if (nzchar(response_text)) {
      collected_warnings <- c(collected_warnings, sprintf("ssgsea review attempt %d returned non-JSON text", attempt))
    }
    if (attempt < standardize_retries) Sys.sleep(standardize_retry_sleep)
  }

  list(
    result = NULL,
    warnings = unique(collected_warnings),
    error = ifelse(is.null(last_error) || !nzchar(last_error), "ssGSEA review LLM failed", last_error)
  )
}

tc_build_record_screening_prompt <- function(batch_df, family_label) {
  required_cols <- c("record_id", "record_label", "annotation_label", "primary_text", "supporting_text")
  for (nm in setdiff(required_cols, names(batch_df))) batch_df[[nm]] <- ""
  payload <- jsonlite::toJSON(
    batch_df[, required_cols, drop = FALSE],
    dataframe = "rows",
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )
  paste(
    sprintf("You are screening %s records to identify likely biological discoveries and likely outliers.", family_label),
    "Return valid JSON only. No markdown, no commentary.",
    "Return a JSON array with the same number of items and the same record_id values as the input.",
    "Each item must contain exactly these keys:",
    "record_id, biological_signal_class, confidence, short_call, discovery_flag, outlier_flag, evidence_summary, followup.",
    "Rules:",
    "1. biological_signal_class must be one of: aligned_state, potential_discovery, likely_outlier, mixed_or_uncertain.",
    "2. confidence must be one of: high, medium, low.",
    "3. discovery_flag and outlier_flag must be yes or no.",
    "4. short_call, evidence_summary, and followup must be concise Simplified Chinese.",
    "5. potential_discovery means a coherent state/substate or tissue-adapted program that is not better explained by contamination or obvious misannotation.",
    "6. likely_outlier means evidence favors misannotation, contamination, doublet, or technical/stress artifact over a coherent lineage-consistent state.",
    "7. aligned_state means mostly matches annotation with no strong discovery or outlier signal.",
    "8. mixed_or_uncertain means evidence is insufficient or conflicting.",
    sprintf("9. %s", tc_noninformative_gene_rule_text()),
    "Input records:",
    payload,
    sep = "\n"
  )
}

tc_run_record_screening_llm_batch <- function(batch_df,
                                              family_label,
                                              caller_env = parent.frame()) {
  enable_llm <- isTRUE(tc_lookup_in_caller("ENABLE_LLM", FALSE, caller_env))
  if (!enable_llm) stop("LLM disabled")
  prompt <- tc_build_record_screening_prompt(batch_df, family_label)
  standardize_model <- tc_lookup_in_caller("STANDARDIZE_LLM_MODEL", "deepseek-chat", caller_env)
  deepseek_api_key <- tc_lookup_in_caller("DEEPSEEK_API_KEY", Sys.getenv("DEEPSEEK_API_KEY", unset = ""), caller_env)
  standardize_retries <- as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_RETRIES", 3L, caller_env))
  standardize_retry_sleep <- tc_lookup_in_caller("STANDARDIZE_LLM_RETRY_SLEEP_SEC", 2, caller_env)

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(max(1L, standardize_retries))) {
    response_text <- tryCatch(
      tc_safe_trim(fanyi::chat_request(prompt, model = standardize_model, api_key = deepseek_api_key)),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )
    parsed <- tc_parse_json(response_text)
    if (is.list(parsed) && length(parsed) > 0) {
      if (!is.null(names(parsed)) && all(c("record_id", "biological_signal_class") %in% names(parsed))) parsed <- list(parsed)
      out <- lapply(parsed, function(item) {
        list(
          record_id = tc_safe_trim(item$record_id),
          biological_signal_class = tc_safe_trim(item$biological_signal_class),
          confidence = tc_safe_trim(item$confidence),
          short_call = tc_safe_trim(item$short_call),
          discovery_flag = tc_safe_trim(item$discovery_flag),
          outlier_flag = tc_safe_trim(item$outlier_flag),
          evidence_summary = tc_safe_trim(item$evidence_summary),
          followup = tc_safe_trim(item$followup)
        )
      })
      ids <- vapply(out, function(x) x$record_id, character(1))
      if (length(out) == nrow(batch_df) && setequal(ids, batch_df$record_id)) {
        return(out[match(batch_df$record_id, ids)])
      }
      collected_warnings <- c(collected_warnings, sprintf("screening attempt %d returned mismatched ids", attempt))
    } else if (nzchar(response_text)) {
      collected_warnings <- c(collected_warnings, sprintf("screening attempt %d returned non-JSON text", attempt))
    }
    if (attempt < standardize_retries) Sys.sleep(standardize_retry_sleep)
  }
  stop(ifelse(is.null(last_error) || !nzchar(last_error), paste(collected_warnings, collapse = " | "), last_error))
}

tc_run_record_screening_llm <- function(records_df,
                                        family_label,
                                        batch_size = 10L,
                                        caller_env = parent.frame()) {
  if (is.null(records_df) || !is.data.frame(records_df) || nrow(records_df) == 0) return(data.frame())
  batch_idx <- split(seq_len(nrow(records_df)), ceiling(seq_len(nrow(records_df)) / max(1L, as.integer(batch_size))))
  retry_sleep <- tc_lookup_in_caller("STANDARDIZE_LLM_RETRY_SLEEP_SEC", 2, caller_env)
  rows <- list()
  for (i in seq_along(batch_idx)) {
    idx <- batch_idx[[i]]
    cat(sprintf("[LLM screen] %s batch %d/%d (%d records)\n", family_label, i, length(batch_idx), length(idx)))
    batch_out <- tc_run_record_screening_llm_batch(
      records_df[idx, , drop = FALSE],
      family_label = family_label,
      caller_env = caller_env
    )
    for (j in seq_along(idx)) {
      rows[[length(rows) + 1L]] <- data.frame(
        record_id = batch_out[[j]]$record_id,
        biological_signal_class = batch_out[[j]]$biological_signal_class,
        confidence = batch_out[[j]]$confidence,
        short_call = batch_out[[j]]$short_call,
        discovery_flag = batch_out[[j]]$discovery_flag,
        outlier_flag = batch_out[[j]]$outlier_flag,
        evidence_summary = batch_out[[j]]$evidence_summary,
        followup = batch_out[[j]]$followup,
        stringsAsFactors = FALSE
      )
    }
    Sys.sleep(retry_sleep)
  }
  out <- dplyr::bind_rows(rows)
  dplyr::left_join(records_df, out, by = "record_id")
}

tc_write_discovery_screen_markdown <- function(screen_df, path, title, family_label) {
  md <- c(title, "", sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
  if (is.null(screen_df) || !is.data.frame(screen_df) || nrow(screen_df) == 0) {
    md <- c(md, "No screening records available.")
    writeLines(md, path)
    return(invisible(path))
  }
  class_counts <- table(screen_df$biological_signal_class, useNA = "ifany")
  md <- c(md, sprintf("**Family:** %s", family_label), "", "## Class Summary", "")
  for (nm in names(class_counts)) {
    md <- c(md, sprintf("- %s: %d", nm, class_counts[[nm]]))
  }
  md <- c(md, "")
  top_section <- function(df, heading, class_value) {
    sub <- df[df$biological_signal_class == class_value, , drop = FALSE]
    if (nrow(sub) == 0) return(character())
    conf_rank <- match(sub$confidence, c("high", "medium", "low"), nomatch = 99L)
    sub <- sub[order(conf_rank, sub$record_id), , drop = FALSE]
    sub <- utils::head(sub, 15)
    out <- c(sprintf("## %s", heading), "")
    for (i in seq_len(nrow(sub))) {
      out <- c(
        out,
        sprintf("### %s", sub$record_label[i]),
        "",
        sprintf("- **Class:** %s | **Confidence:** %s", sub$biological_signal_class[i], sub$confidence[i]),
        if ("annotation_label" %in% colnames(sub) && nzchar(tc_safe_trim(sub$annotation_label[i]))) sprintf("- **Annotated label:** %s", sub$annotation_label[i]) else NULL,
        sprintf("- **Short call:** %s", sub$short_call[i]),
        sprintf("- **Evidence:** %s", sub$evidence_summary[i]),
        sprintf("- **Follow-up:** %s", sub$followup[i]),
        ""
      )
    }
    out
  }
  md <- c(
    md,
    top_section(screen_df, "Top potential discoveries", "potential_discovery"),
    top_section(screen_df, "Top likely outliers", "likely_outlier"),
    top_section(screen_df, "Mixed or uncertain records", "mixed_or_uncertain")
  )
  writeLines(md, path)
  invisible(path)
}

tc_empty_interpretation_records_df <- function() {
  df <- data.frame(
    celltype_level = character(),
    celltype_label = character(),
    celltype_l2 = character(),
    comparison = character(),
    direction = character(),
    source_db = character(),
    status = character(),
    warnings = character(),
    error = character(),
    overview = character(),
    key_mechanisms = character(),
    hypothesis = character(),
    narrative = character(),
    key_drivers = character(),
    evidence = character(),
    limitations = character(),
    raw_text = character(),
    stringsAsFactors = FALSE
  )
  optional_char_fields <- c(
    "confidence", "cell_type_judgment", "annotation_match_degree", "annotated_l3_correspondence",
    "outlier_assessment", "discovery_assessment", "integrated_diagnostic_comment",
    "annotated_l2_dominant", "annotated_l2_purity_class",
    "annotated_l3_dominant", "annotated_l3_secondary", "annotated_l3_purity_class",
    "annotation_count_summary"
  )
  optional_num_fields <- c("annotated_l2_dominant_pct", "annotated_l3_dominant_pct", "annotated_l3_secondary_pct")
  for (nm in optional_char_fields) df[[nm]] <- character()
  for (nm in optional_num_fields) df[[nm]] <- numeric()
  df
}

tc_empty_discovery_input_df <- function() {
  data.frame(
    record_id = character(),
    record_label = character(),
    annotation_label = character(),
    primary_text = character(),
    supporting_text = character(),
    stringsAsFactors = FALSE
  )
}

tc_empty_discovery_screen_df <- function() {
  df <- tc_empty_discovery_input_df()
  df$biological_signal_class <- character()
  df$confidence <- character()
  df$short_call <- character()
  df$discovery_flag <- character()
  df$outlier_flag <- character()
  df$evidence_summary <- character()
  df$followup <- character()
  df
}

tc_build_choir_ofa_screening_table <- function(choir_ofa_all,
                                               caller_env = parent.frame()) {
  if (is.null(choir_ofa_all) || length(choir_ofa_all) == 0) return(tc_empty_discovery_input_df())
  prepare_bundle <- tc_lookup_in_caller("prepare_llm_enrichment_bundle", NULL, caller_env)
  build_context_text <- tc_lookup_in_caller("build_multi_db_context_text", NULL, caller_env)
  build_deg_text <- tc_lookup_in_caller("build_choir_deg_evidence_text", tc_build_choir_deg_evidence_text, caller_env)
  if (!is.function(prepare_bundle) || !is.function(build_context_text) || !is.function(build_deg_text)) {
    stop("prepare_llm_enrichment_bundle(), build_multi_db_context_text(), and build_choir_deg_evidence_text() must be available in caller environment.")
  }

  rows <- list()
  for (cl_name in names(choir_ofa_all)) {
    rec <- choir_ofa_all[[cl_name]]
    de_df <- rec$de_table
    enrich <- rec$enrich
    l2_label <- if (!is.null(rec$llm$integrated$celltype_l2)) tc_safe_trim(rec$llm$integrated$celltype_l2) else ""
    annotation_label <- if (!is.null(rec$llm$integrated$annotated_l3_dominant) && nzchar(tc_safe_trim(rec$llm$integrated$annotated_l3_dominant))) {
      tc_safe_trim(rec$llm$integrated$annotated_l3_dominant)
    } else if (nzchar(l2_label)) l2_label else paste0("CHOIR_", cl_name)
    enrich_parts <- lapply(names(enrich), function(dir_name) {
      bundle <- prepare_bundle(enrich[[dir_name]], max_dbs = 4L, n_terms_per_db = 3L)
      text <- build_context_text(bundle$evidence_df)
      if (!nzchar(tc_safe_trim(text))) return(NULL)
      paste0("[", dir_name, "]\n", text)
    })
    rows[[length(rows) + 1L]] <- data.frame(
      record_id = paste0("ofa_", cl_name),
      record_label = paste0("cluster_", cl_name, "_vs_rest"),
      annotation_label = annotation_label,
      primary_text = tc_append_nonempty_blocks(
        if (!is.null(rec$llm$integrated$annotation_count_summary)) tc_safe_trim(rec$llm$integrated$annotation_count_summary) else "",
        build_deg_text(
          de_df,
          cluster_id = cl_name,
          top_n_abs = 12L,
          top_n_up = 8L,
          top_n_down = 8L,
          pct1_label = paste0("CHOIR_", cl_name),
          pct2_label = "rest"
        )
      ),
      supporting_text = tc_append_nonempty_blocks(unlist(enrich_parts, use.names = FALSE)),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(tc_empty_discovery_input_df())
  dplyr::bind_rows(rows)
}

tc_run_and_write_discovery_screen <- function(records_df,
                                              family_label,
                                              output_dir,
                                              prefix,
                                              title,
                                              batch_size = 10L,
                                              caller_env = parent.frame()) {
  if (is.null(records_df) || !is.data.frame(records_df)) records_df <- tc_empty_discovery_input_df()
  required_input_cols <- c("record_id", "record_label", "annotation_label", "primary_text", "supporting_text")
  for (nm in setdiff(required_input_cols, colnames(records_df))) records_df[[nm]] <- character(nrow(records_df))
  records_df <- records_df[, required_input_cols, drop = FALSE]

  tsv_path <- file.path(output_dir, "reports", sprintf("%s_discovery_screen.tsv", prefix))
  md_path <- file.path(output_dir, sprintf("%s_DISCOVERY_REVIEW.md", toupper(prefix)))

  if (nrow(records_df) == 0) {
    screened <- tc_empty_discovery_screen_df()
    tc_write_table_file(screened, tsv_path, sep = "\t")
    tc_write_discovery_screen_markdown(
      screened,
      path = md_path,
      title = title,
      family_label = family_label
    )
    return(screened)
  }

  screened <- tc_run_record_screening_llm(
    records_df,
    family_label = family_label,
    batch_size = batch_size,
    caller_env = caller_env
  )
  tc_write_table_file(screened, tsv_path, sep = "\t")
  tc_write_discovery_screen_markdown(
    screened,
    path = md_path,
    title = title,
    family_label = family_label
  )
  screened
}

tc_run_choir_cluster_llm <- function(cluster_id,
                                     de_df,
                                     ofa_enrich,
                                     ssgsea_choir_all = NULL,
                                     choir_cluster_l2 = NA_character_,
                                     choir_annotation_tables = NULL,
                                     caller_env = parent.frame()) {
  enable_llm <- isTRUE(tc_lookup_in_caller("ENABLE_LLM", FALSE, caller_env))
  if (!enable_llm || is.null(de_df) || !is.data.frame(de_df) || nrow(de_df) == 0) return(NULL)

  annotation_bundle <- tc_build_choir_annotation_evidence(choir_annotation_tables, cluster_id = cluster_id)
  if (is.na(choir_cluster_l2) || !nzchar(trimws(choir_cluster_l2))) {
    choir_cluster_l2 <- annotation_bundle$l2$dominant_label
  }

  pick_col <- tc_lookup_in_caller("llm_pick_existing_col", tc_llm_pick_existing_col, caller_env)
  filter_named_gene_fc <- tc_lookup_in_caller(
    "tc_filter_named_gene_fc",
    function(gene_fc) {
      if (is.null(gene_fc) || length(gene_fc) == 0) return(numeric())
      keep <- !is.na(names(gene_fc)) & nzchar(trimws(names(gene_fc)))
      stats::setNames(as.numeric(gene_fc[keep]), toupper(names(gene_fc[keep])))
    },
    caller_env
  )
  prepare_bundle <- tc_lookup_in_caller("prepare_llm_enrichment_bundle", NULL, caller_env)
  if (!is.function(prepare_bundle)) {
    stop("prepare_llm_enrichment_bundle() is required in caller environment for tc_run_choir_cluster_llm().")
  }
  format_gene_pathway_map <- tc_lookup_in_caller(
    "format_gene_pathway_map_text",
    function(gene_pathway_map) {
      if (is.null(gene_pathway_map) || nrow(gene_pathway_map) == 0) return(character())
      apply(gene_pathway_map, 1, function(row) {
        sprintf(
          "- %s(log2FC=%.2f) -> %s | pathway_padj=%s",
          row[["gene"]], as.numeric(row[["log2FC"]]), row[["pathway"]],
          format(as.numeric(row[["pathway_padj"]]), scientific = TRUE, digits = 3)
        )
      })
    },
    caller_env
  )
  allow_comparative_hypothesis <- isTRUE(tc_lookup_in_caller("LLM_ALLOW_COMPARATIVE_HYPOTHESIS", TRUE, caller_env))
  choi_top_deg_n <- tc_lookup_in_caller("CHOIR_LLM_TOP_DEG_N", 20L, caller_env)
  choir_top_enrich_terms <- tc_lookup_in_caller("CHOIR_LLM_TOP_ENRICH_TERMS_PER_DB", 6L, caller_env)
  choir_top_ssgsea_terms <- tc_lookup_in_caller("CHOIR_LLM_TOP_SSGSEA_TERMS_PER_DIRECTION", 8L, caller_env)
  choir_max_dbs <- tc_lookup_in_caller("CHOIR_LLM_MAX_DBS", 8L, caller_env)
  line_base_context <- tc_lookup_in_caller("LINEAGE_BASE_CONTEXT", "", caller_env)

  deg_text <- tc_build_choir_deg_evidence_text(
    de_df = de_df,
    cluster_id = cluster_id,
    top_n_abs = choi_top_deg_n,
    pct1_label = paste0("CHOIR_", cluster_id),
    pct2_label = "rest",
    caller_env = caller_env
  )
  gene_col <- pick_col(de_df, c("gene", "symbol", "feature", "features", "genes"))
  logfc_col <- pick_col(de_df, c("avg_log2FC", "log2FoldChange", "log2FC", "avg_logFC"))
  gene_fc_all <- if (!is.null(gene_col) && !is.null(logfc_col)) {
    filter_named_gene_fc(stats::setNames(tc_llm_num(de_df[[logfc_col]]), toupper(as.character(de_df[[gene_col]]))))
  } else numeric()

  bundle_list <- list()
  for (dir_name in c("up", "down")) {
    if (is.null(ofa_enrich[[dir_name]])) next
    gene_fc_dir <- if (dir_name == "up") gene_fc_all[gene_fc_all > 0] else gene_fc_all[gene_fc_all < 0]
    bundle_list[[dir_name]] <- prepare_bundle(
      ofa_enrich[[dir_name]],
      gene_fc = gene_fc_dir,
      max_dbs = choir_max_dbs,
      n_terms_per_db = choir_top_enrich_terms
    )
  }

  evidence_rows <- lapply(names(bundle_list), function(dir_name) {
    df <- bundle_list[[dir_name]]$evidence_df
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$direction <- dir_name
    df
  })
  evidence_df <- dplyr::bind_rows(evidence_rows)
  evidence_text <- tc_append_nonempty_blocks(
    annotation_bundle$text,
    deg_text,
    tc_build_choir_cluster_ssgsea_evidence(
      ssgsea_choir_all,
      cluster_id,
      n_top_per_method = choir_top_ssgsea_terms,
      caller_env = caller_env
    ),
    tc_build_integrated_directional_enrichment_text(
      evidence_df,
      direction_map = c(up = "cluster_gt_rest", down = "rest_gt_cluster"),
      caller_env = caller_env
    )
  )
  if (!nzchar(tc_safe_trim(evidence_text))) return(NULL)

  source_dbs <- unique(unlist(lapply(bundle_list, function(x) x$selected_dbs), use.names = FALSE))
  source_db_label <- if (length(source_dbs) > 0) {
    paste0(tc_collapse_source_db_label(source_dbs), " + ssGSEA + annotation_counts")
  } else {
    "ssGSEA + annotation_counts"
  }
  gene_pathway_map <- tc_merge_directional_gene_pathway_maps(bundle_list)
  gene_pathway_map_text <- paste(format_gene_pathway_map(gene_pathway_map), collapse = "\n")

  choir_ctx <- sprintf(
    paste(
      "%s",
      "\nCHOIR cluster %s vs rest (integrated one-vs-rest interpretation).",
      "\nDominant coarse label: %s.",
      "\nDominant user-annotated L3 label: %s (%s); secondary L3 label: %s (%s); purity class=%s.",
      "\nBiological question: What cell type/state is CHOIR cluster %s most consistent with,",
      "how well does that judgment correspond to the user-annotated L3 label, and is this cluster better explained as aligned annotation, outlier/misannotation, mixed transition state, stress/contamination artifact, or a plausible biological discovery when integrating annotation composition, complete DEG, ssGSEA, and OFA evidence?",
      "\n\nOutput requirements:",
      "1. Write ALL interpretive content in Chinese (Simplified Chinese characters).",
      "2. key_drivers must remain as English gene symbols.",
      "3. Integrate up and down evidence into one cell-type/state judgment.",
      "4. Use cluster-level ssGSEA mainly for cell-type/state judgment.",
      "5. Explicitly compare the inferred identity with the user-annotated L3 composition and state the correspondence degree.",
      "6. Explicitly distinguish likely annotation outliers from plausible biological discoveries.",
      "7. Keep the interpretation concise but evidence-dense; only propose a brief scientific hypothesis when the comparison supports it.",
      "8. Do not invent gene->pathway links beyond the validated map."
    ),
    line_base_context,
    cluster_id,
    ifelse(is.na(choir_cluster_l2) || !nzchar(trimws(choir_cluster_l2)), "NA", choir_cluster_l2),
    ifelse(is.null(annotation_bundle$l3$dominant_label) || is.na(annotation_bundle$l3$dominant_label), "NA", annotation_bundle$l3$dominant_label),
    tc_pct_label(annotation_bundle$l3$dominant_pct),
    ifelse(is.null(annotation_bundle$l3$secondary_label) || is.na(annotation_bundle$l3$secondary_label), "NA", annotation_bundle$l3$secondary_label),
    tc_pct_label(annotation_bundle$l3$secondary_pct),
    ifelse(is.null(annotation_bundle$l3$purity_class) || !nzchar(annotation_bundle$l3$purity_class), "unknown", annotation_bundle$l3$purity_class),
    cluster_id
  )

  std_payload <- tc_run_choir_cluster_review_llm(
    cluster_id = cluster_id,
    choir_ctx = tc_append_nonempty_blocks(
      choir_ctx,
      if (nzchar(gene_pathway_map_text)) paste("Validated gene-pathway map:\n", gene_pathway_map_text, sep = "") else ""
    ),
    evidence_text = evidence_text,
    source_db_label = source_db_label,
    warnings = character(),
    error_message = NULL,
    caller_env = caller_env
  )

  if (!is.null(std_payload$result)) {
    std <- std_payload$result
    return(list(
      celltype_level = "CHOIR",
      celltype_label = paste0("CHOIR_", cluster_id),
      celltype_l2 = choir_cluster_l2,
      comparison = paste0("cluster_", cluster_id, "_vs_rest"),
      direction = "integrated_up_down",
      source_db = source_db_label,
      status = "structured",
      warnings = unique(std_payload$warnings),
      error = "",
      confidence = tc_placeholder_text(std$confidence, default = "medium"),
      cell_type_judgment = tc_placeholder_text(std$cell_type_judgment),
      annotation_match_degree = tc_placeholder_text(std$annotation_match_degree),
      annotated_l3_correspondence = tc_placeholder_text(std$annotated_l3_correspondence),
      outlier_assessment = tc_placeholder_text(std$outlier_assessment),
      discovery_assessment = tc_placeholder_text(std$discovery_assessment),
      integrated_diagnostic_comment = tc_placeholder_text(std$integrated_diagnostic_comment),
      annotated_l2_dominant = ifelse(is.null(annotation_bundle$l2$dominant_label), "", annotation_bundle$l2$dominant_label),
      annotated_l2_dominant_pct = ifelse(is.null(annotation_bundle$l2$dominant_pct), NA_real_, annotation_bundle$l2$dominant_pct),
      annotated_l2_purity_class = ifelse(is.null(annotation_bundle$l2$purity_class), "", annotation_bundle$l2$purity_class),
      annotated_l3_dominant = ifelse(is.null(annotation_bundle$l3$dominant_label), "", annotation_bundle$l3$dominant_label),
      annotated_l3_dominant_pct = ifelse(is.null(annotation_bundle$l3$dominant_pct), NA_real_, annotation_bundle$l3$dominant_pct),
      annotated_l3_secondary = ifelse(is.null(annotation_bundle$l3$secondary_label), "", annotation_bundle$l3$secondary_label),
      annotated_l3_secondary_pct = ifelse(is.null(annotation_bundle$l3$secondary_pct), NA_real_, annotation_bundle$l3$secondary_pct),
      annotated_l3_purity_class = ifelse(is.null(annotation_bundle$l3$purity_class), "", annotation_bundle$l3$purity_class),
      annotation_count_summary = ifelse(is.null(annotation_bundle$text), "", annotation_bundle$text),
      overview = tc_placeholder_text(std$overview),
      key_mechanisms = tc_placeholder_text(std$key_mechanisms),
      hypothesis = tc_placeholder_text(std$hypothesis),
      narrative = tc_placeholder_text(std$narrative),
      key_drivers = tc_placeholder_text(tc_collapse_driver_field(std$key_drivers)),
      evidence = tc_placeholder_text(std$evidence),
      limitations = tc_placeholder_text(std$limitations),
      raw_text = evidence_text,
      raw_result = std_payload$result
    ))
  }

  list(
    celltype_level = "CHOIR",
    celltype_label = paste0("CHOIR_", cluster_id),
    celltype_l2 = choir_cluster_l2,
    comparison = paste0("cluster_", cluster_id, "_vs_rest"),
    direction = "integrated_up_down",
    source_db = source_db_label,
    status = "error",
    warnings = unique(std_payload$warnings),
    error = ifelse(is.null(std_payload$error), "CHOIR LLM standardization failed", std_payload$error),
    confidence = "low",
    cell_type_judgment = "Not available from current evidence.",
    annotation_match_degree = "low",
    annotated_l3_correspondence = tc_placeholder_text(annotation_bundle$l3$summary_text),
    outlier_assessment = "当前因 LLM 结构化失败，无法可靠判断是否为标注异常。",
    discovery_assessment = "当前因 LLM 结构化失败，无法可靠判断是否为潜在生物学发现。",
    integrated_diagnostic_comment = "当前仅保留原始证据汇总；需要重新运行 LLM 以完成异常值/发现判别。",
    annotated_l2_dominant = ifelse(is.null(annotation_bundle$l2$dominant_label), "", annotation_bundle$l2$dominant_label),
    annotated_l2_dominant_pct = ifelse(is.null(annotation_bundle$l2$dominant_pct), NA_real_, annotation_bundle$l2$dominant_pct),
    annotated_l2_purity_class = ifelse(is.null(annotation_bundle$l2$purity_class), "", annotation_bundle$l2$purity_class),
    annotated_l3_dominant = ifelse(is.null(annotation_bundle$l3$dominant_label), "", annotation_bundle$l3$dominant_label),
    annotated_l3_dominant_pct = ifelse(is.null(annotation_bundle$l3$dominant_pct), NA_real_, annotation_bundle$l3$dominant_pct),
    annotated_l3_secondary = ifelse(is.null(annotation_bundle$l3$secondary_label), "", annotation_bundle$l3$secondary_label),
    annotated_l3_secondary_pct = ifelse(is.null(annotation_bundle$l3$secondary_pct), NA_real_, annotation_bundle$l3$secondary_pct),
    annotated_l3_purity_class = ifelse(is.null(annotation_bundle$l3$purity_class), "", annotation_bundle$l3$purity_class),
    annotation_count_summary = ifelse(is.null(annotation_bundle$text), "", annotation_bundle$text),
    overview = tc_placeholder_text(""),
    key_mechanisms = tc_placeholder_text(""),
    hypothesis = tc_placeholder_text(if (isTRUE(allow_comparative_hypothesis)) "" else "Not applicable."),
    narrative = tc_placeholder_text(""),
    key_drivers = tc_placeholder_text(""),
    evidence = tc_placeholder_text(evidence_text),
    limitations = tc_placeholder_text("Evidence integration failed during LLM standardization."),
    raw_text = evidence_text,
    raw_result = NULL
  )
}

tc_collect_choir_llm_records <- function(choir_ofa) {
  rows <- list()
  for (cl_name in names(choir_ofa)) {
    llm_list <- choir_ofa[[cl_name]][["llm"]]
    if (is.null(llm_list) || length(llm_list) == 0) next
    for (dir_name in names(llm_list)) {
      rec <- llm_list[[dir_name]]
      if (is.null(rec)) next
      row_df <- data.frame(
        celltype_level = tc_safe_trim(rec$celltype_level),
        celltype_label = tc_safe_trim(rec$celltype_label),
        celltype_l2 = tc_safe_trim(rec$celltype_l2),
        comparison = tc_safe_trim(rec$comparison),
        direction = tc_safe_trim(rec$direction),
        source_db = tc_safe_trim(rec$source_db),
        status = tc_safe_trim(rec$status),
        warnings = paste(rec$warnings, collapse = " | "),
        error = tc_safe_trim(rec$error),
        overview = tc_safe_trim(rec$overview),
        key_mechanisms = tc_safe_trim(rec$key_mechanisms),
        hypothesis = tc_safe_trim(rec$hypothesis),
        narrative = tc_safe_trim(rec$narrative),
        key_drivers = tc_safe_trim(rec$key_drivers),
        evidence = tc_safe_trim(rec$evidence),
        limitations = tc_safe_trim(rec$limitations),
        raw_text = tc_safe_trim(rec$raw_text),
        stringsAsFactors = FALSE
      )
      optional_char_fields <- c(
        "confidence", "cell_type_judgment", "annotation_match_degree", "annotated_l3_correspondence",
        "outlier_assessment", "discovery_assessment", "integrated_diagnostic_comment",
        "annotated_l2_dominant", "annotated_l2_purity_class",
        "annotated_l3_dominant", "annotated_l3_secondary", "annotated_l3_purity_class",
        "annotation_count_summary"
      )
      optional_num_fields <- c("annotated_l2_dominant_pct", "annotated_l3_dominant_pct", "annotated_l3_secondary_pct")
      for (nm in optional_char_fields) {
        row_df[[nm]] <- if (!is.null(rec[[nm]])) tc_safe_trim(rec[[nm]]) else ""
      }
      for (nm in optional_num_fields) {
        row_df[[nm]] <- if (!is.null(rec[[nm]])) suppressWarnings(as.numeric(rec[[nm]])) else NA_real_
      }
      rows[[length(rows) + 1L]] <- row_df
    }
  }
  if (length(rows) == 0) return(tc_empty_interpretation_records_df())
  dplyr::bind_rows(rows)
}

tc_safe_slug <- function(x) {
  x <- tolower(tc_safe_trim(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "panel" else x
}

tc_match_features_to_available <- function(features, available_features) {
  features <- as.character(features)
  features <- features[!is.na(features) & nzchar(trimws(features))]
  available_features <- as.character(available_features)
  available_features <- available_features[!is.na(available_features) & nzchar(trimws(available_features))]
  if (length(features) == 0 || length(available_features) == 0) return(character())
  feature_map <- setNames(available_features, toupper(available_features))
  matched <- unname(feature_map[toupper(features)])
  unique(matched[!is.na(matched) & nzchar(trimws(matched))])
}

tc_extract_custom_panel_markers <- function(custom_markers_db, subtype_pattern = NULL) {
  if (is.null(custom_markers_db) || !is.data.frame(custom_markers_db) || nrow(custom_markers_db) == 0) return(character())
  if (!all(c("subtype", "markers") %in% colnames(custom_markers_db))) return(character())
  db_use <- custom_markers_db
  if (!is.null(subtype_pattern) && length(subtype_pattern) == 1 && nzchar(tc_safe_trim(subtype_pattern))) {
    keep <- grepl(subtype_pattern, as.character(db_use$subtype), ignore.case = TRUE, perl = TRUE)
    db_use <- db_use[keep, , drop = FALSE]
  }
  if (nrow(db_use) == 0) return(character())
  markers <- unlist(strsplit(as.character(db_use$markers), ",", fixed = TRUE), use.names = FALSE)
  markers <- trimws(toupper(markers))
  unique(markers[!is.na(markers) & nzchar(markers)])
}

tc_default_bcell_marker_panels <- function(known_markers = NULL, custom_markers_db = NULL) {
  panel_genes <- function(base_genes, subtype_pattern = NULL) {
    unique(c(
      toupper(trimws(as.character(base_genes))),
      tc_extract_custom_panel_markers(custom_markers_db, subtype_pattern = subtype_pattern)
    ))
  }

  list(
    pan_b = list(
      label = "Pan-B",
      genes = unique(toupper(c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "CD74", "LAPTM5", "HLA-DRA")))
    ),
    naive_b = list(
      label = "Naive / Transitional",
      genes = panel_genes(
        c("IGHD", "IGHM", "TCL1A", "FCER2", "IL4R", "CD24", "CD38", "MME"),
        subtype_pattern = "NAIVE|TRANSITIONAL"
      )
    ),
    memory_b = list(
      label = "Memory / Atypical",
      genes = panel_genes(
        c("CD27", "TNFRSF13B", "AIM2", "ITGAX", "TBX21", "FCRL4", "FCRL5", "CCR6"),
        subtype_pattern = "MEMORY|ATYPICAL"
      )
    ),
    gc_b = list(
      label = "GC B",
      genes = panel_genes(
        c("BCL6", "AICDA", "RGS13", "MEF2B", "MME", "MKI67", "CXCR4", "CD83"),
        subtype_pattern = "^GC_|GERMINAL"
      )
    ),
    plasma = list(
      label = "Plasma",
      genes = panel_genes(
        c("PRDM1", "XBP1", "MZB1", "JCHAIN", "SDC1", "DERL3", "SSR4", "IGHA1", "IGHG1", "IGHE"),
        subtype_pattern = "PLASMA|PLASMABLAST"
      )
    )
  )
}

tc_default_tnk_marker_panels <- function(known_markers = NULL, custom_markers_db = NULL) {
  panel_genes <- function(base_genes, subtype_pattern = NULL) {
    unique(c(
      toupper(trimws(as.character(base_genes))),
      tc_extract_custom_panel_markers(custom_markers_db, subtype_pattern = subtype_pattern)
    ))
  }

  list(
    pan_t = list(
      label = "Pan-T / Lymphocyte",
      genes = unique(toupper(c("CD3D", "CD3E", "TRAC", "CD247", "IL7R", "LTB")))
    ),
    cd4_helper = list(
      label = "CD4 helper / regulatory",
      genes = panel_genes(
        c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "CXCR5", "PDCD1", "ICOS", "BCL6", "FOXP3", "IL2RA", "CTLA4"),
        subtype_pattern = "CD4|TREG|TFH|TFR|TH1|TH17|TRM"
      )
    ),
    cd8_cytotoxic = list(
      label = "CD8 / cytotoxic",
      genes = panel_genes(
        c("CD8A", "CD8B", "NKG7", "CCL5", "PRF1", "GZMB", "GZMK", "FGFBP2", "KLRG1", "CX3CR1"),
        subtype_pattern = "CD8|CYTOTOXIC|TEM|TEMRA|TEFF|TRM|GDT|MAIT"
      )
    ),
    nk_ilc = list(
      label = "NK / innate-like",
      genes = panel_genes(
        c("FCGR3A", "KLRD1", "NCR1", "GNLY", "XCL1", "XCL2", "KLRB1", "TRDC", "TRGC1", "KIT", "IL23R"),
        subtype_pattern = "NK|ILC|MAIT|GDT"
      )
    ),
    dysfunction_residency = list(
      label = "Residency / dysfunction",
      genes = panel_genes(
        c("CD69", "ITGAE", "CXCR6", "ZNF683", "TIGIT", "LAG3", "HAVCR2", "TOX", "LAYN"),
        subtype_pattern = "TRM|EXHAUST|TREG|TFR"
      )
    )
  )
}

tc_default_epithelial_marker_panels <- function(known_markers = NULL, custom_markers_db = NULL) {
  panel_genes <- function(base_genes, subtype_pattern = NULL) {
    unique(c(
      toupper(trimws(as.character(base_genes))),
      tc_extract_custom_panel_markers(custom_markers_db, subtype_pattern = subtype_pattern)
    ))
  }

  list(
    pan_epithelial = list(
      label = "Pan-epithelial",
      genes = unique(toupper(c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19")))
    ),
    alveolar = list(
      label = "Alveolar",
      genes = panel_genes(
        c("AGER", "HOPX", "CAV1", "SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NAPSA", "SLC34A2"),
        subtype_pattern = "AT1|AT2|ALVEOLAR"
      )
    ),
    basal_suprabasal = list(
      label = "Basal / suprabasal",
      genes = panel_genes(
        c("TP63", "KRT5", "KRT14", "KRT15", "ITGA6", "NGFR", "KRT4", "KRT13", "KRT17"),
        subtype_pattern = "BASAL|SUPRABASAL"
      )
    ),
    ciliated = list(
      label = "Ciliated",
      genes = panel_genes(
        c("FOXJ1", "TPPP3", "DNAH5", "DNAH9", "RSPH1", "DEUP1", "CCNO", "MCIDAS"),
        subtype_pattern = "CILIATED|CILIOGENESIS|DEUTEROSOMAL"
      )
    ),
    secretory_smg = list(
      label = "Secretory / SMG",
      genes = panel_genes(
        c("SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5AC", "MUC5B", "SPDEF", "FOXA3", "LTF", "LYZ", "SLPI", "DMBT1", "BPIFA1", "PIGR"),
        subtype_pattern = "GOBLET|CLUB|SMG|SECRETORY|DUOX2"
      )
    ),
    rare_specialized = list(
      label = "Rare specialized / squamous",
      genes = panel_genes(
        c("FOXI1", "ASCL3", "CFTR", "ATP6V0D2", "CLCNKA", "CLCNKB", "SPRR1A", "IVL", "KRT6A", "S100A7"),
        subtype_pattern = "IONOCYTE|BRUSH|SQUAMOUS"
      )
    )
  )
}

tc_normalize_marker_panels <- function(marker_panels, available_features = NULL) {
  if (is.null(marker_panels) || length(marker_panels) == 0) return(list())
  panel_ids <- names(marker_panels)
  if (is.null(panel_ids) || any(!nzchar(panel_ids))) {
    panel_ids <- sprintf("panel_%02d", seq_along(marker_panels))
  }
  normalized <- list()
  for (i in seq_along(marker_panels)) {
    panel <- marker_panels[[i]]
    panel_id <- panel_ids[[i]]
    if (is.list(panel) && !is.null(panel$genes)) {
      panel_label <- if (!is.null(panel$label) && nzchar(tc_safe_trim(panel$label))) tc_safe_trim(panel$label) else tc_safe_trim(panel_id)
      genes <- panel$genes
    } else {
      panel_label <- tc_safe_trim(panel_id)
      genes <- panel
    }
    genes <- unique(toupper(trimws(as.character(genes))))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!is.null(available_features)) {
      genes <- tc_match_features_to_available(genes, available_features)
    }
    if (length(genes) == 0) next
    normalized[[length(normalized) + 1L]] <- list(
      panel_id = panel_id,
      panel_slug = tc_safe_slug(panel_id),
      panel_label = if (nzchar(panel_label)) panel_label else panel_id,
      genes = genes
    )
  }
  names(normalized) <- vapply(normalized, function(x) x$panel_id, character(1))
  normalized
}

tc_save_plot <- function(plot_obj, path_no_ext, width = 10, height = 8, dpi = 300) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot_obj, width = width, height = height, dpi = dpi)
  invisible(list(
    pdf = paste0(path_no_ext, ".pdf"),
    png = paste0(path_no_ext, ".png")
  ))
}

tc_downsample_object_by_group <- function(obj, group_col, n_per_group = 100L, seed = 42L) {
  if (is.null(obj) || is.null(obj@meta.data) || !group_col %in% colnames(obj@meta.data)) return(obj)
  set.seed(seed)
  meta <- obj@meta.data
  groups <- unique(as.character(meta[[group_col]]))
  groups <- groups[!is.na(groups) & nzchar(trimws(groups))]
  keep_cells <- unlist(lapply(groups, function(g) {
    hits <- rownames(meta)[as.character(meta[[group_col]]) == g]
    if (length(hits) <= n_per_group) return(hits)
    sample(hits, n_per_group)
  }), use.names = FALSE)
  keep_cells <- unique(keep_cells)
  if (length(keep_cells) == 0) return(obj)
  obj[, keep_cells]
}

tc_generate_marker_panel_visualizations <- function(obj,
                                                    marker_panels,
                                                    fig_dir,
                                                    group_col,
                                                    lineage_display = "Lineage",
                                                    reduction_name = NULL,
                                                    report_dir = NULL,
                                                    panel_prefix = "marker_panel",
                                                    downsample_n = 100L,
                                                    verbose = TRUE) {
  if (is.null(obj) || is.null(obj@meta.data)) stop("obj must be a Seurat-like object with meta.data.")
  if (missing(group_col) || !nzchar(tc_safe_trim(group_col)) || !group_col %in% colnames(obj@meta.data)) {
    stop("group_col must exist in obj@meta.data.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(report_dir) && nzchar(tc_safe_trim(report_dir))) {
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  }

  panels <- tc_normalize_marker_panels(marker_panels, available_features = rownames(obj))
  if (length(panels) == 0) {
    return(list(summary = data.frame(), panels = list(), figure_dir = fig_dir))
  }

  groups_present <- unique(as.character(obj@meta.data[[group_col]]))
  groups_present <- groups_present[!is.na(groups_present) & nzchar(trimws(groups_present))]
  n_groups <- length(groups_present)
  obj_heatmap <- tc_downsample_object_by_group(obj, group_col = group_col, n_per_group = downsample_n)

  summary_rows <- list()
  panel_results <- list()

  for (panel in panels) {
    genes <- panel$genes
    if (length(genes) < 2) next
    panel_slug <- panel$panel_slug
    panel_label <- panel$panel_label
    title_stub <- sprintf("%s - %s markers", lineage_display, panel_label)

    dotplot_stem <- file.path(fig_dir, sprintf("%s_dotplot_%s", panel_prefix, panel_slug))
    featureplot_stem <- file.path(fig_dir, sprintf("%s_featureplot_%s", panel_prefix, panel_slug))
    heatmap_stem <- file.path(fig_dir, sprintf("%s_heatmap_%s", panel_prefix, panel_slug))

    dotplot_paths <- tryCatch({
      p_dot <- Seurat::DotPlot(obj, features = genes, group.by = group_col) +
        Seurat::RotatedAxis() +
        ggplot2::ggtitle(sprintf("%s (%s)", title_stub, group_col)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7))
      tc_save_plot(
        p_dot,
        dotplot_stem,
        width = min(20, max(10, length(genes) * 0.6)),
        height = min(18, max(5, n_groups * 0.4 + 2))
      )
    }, error = function(e) {
      if (isTRUE(verbose)) cat(sprintf("[WARN] Marker-panel dotplot failed for %s: %s\n", panel_label, e$message))
      list(pdf = NA_character_, png = NA_character_)
    })

    featureplot_paths <- list(pdf = NA_character_, png = NA_character_)
    reduction_available <- !is.null(reduction_name) && nzchar(tc_safe_trim(reduction_name)) &&
      reduction_name %in% names(obj@reductions)
    if (isTRUE(reduction_available)) {
      featureplot_paths <- tryCatch({
        ncol_plot <- min(4L, length(genes))
        nrow_plot <- ceiling(length(genes) / ncol_plot)
        p_feature <- Seurat::FeaturePlot(
          obj,
          features = genes,
          reduction = reduction_name,
          combine = TRUE,
          ncol = ncol_plot,
          order = TRUE
        ) + patchwork::plot_annotation(title = sprintf("%s (%s)", title_stub, reduction_name))
        tc_save_plot(
          p_feature,
          featureplot_stem,
          width = min(18, max(10, ncol_plot * 3.8)),
          height = min(18, max(6, nrow_plot * 3.8 + 1))
        )
      }, error = function(e) {
        if (isTRUE(verbose)) cat(sprintf("[WARN] Marker-panel featureplot failed for %s: %s\n", panel_label, e$message))
        list(pdf = NA_character_, png = NA_character_)
      })
    }

    heatmap_paths <- tryCatch({
      hm_obj <- obj_heatmap
      hm_obj <- Seurat::ScaleData(hm_obj, features = genes, verbose = FALSE)
      p_heat <- Seurat::DoHeatmap(hm_obj, features = genes, group.by = group_col, size = 3) +
        ggplot2::ggtitle(sprintf("%s (downsampled, visualization only)", title_stub))
      tc_save_plot(
        p_heat,
        heatmap_stem,
        width = min(20, max(12, length(genes) * 0.55)),
        height = min(18, max(6, n_groups * 0.45 + 2))
      )
    }, error = function(e) {
      if (isTRUE(verbose)) cat(sprintf("[WARN] Marker-panel heatmap failed for %s: %s\n", panel_label, e$message))
      list(pdf = NA_character_, png = NA_character_)
    })

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      panel_id = panel$panel_id,
      panel_slug = panel_slug,
      panel_label = panel_label,
      n_genes = length(genes),
      genes = paste(genes, collapse = ", "),
      dotplot_png = tc_to_scalar(dotplot_paths$png),
      featureplot_png = tc_to_scalar(featureplot_paths$png),
      heatmap_png = tc_to_scalar(heatmap_paths$png),
      stringsAsFactors = FALSE
    )
    panel_results[[panel$panel_id]] <- list(
      genes = genes,
      dotplot = dotplot_paths,
      featureplot = featureplot_paths,
      heatmap = heatmap_paths
    )
    if (isTRUE(verbose)) cat(sprintf("[OK] Marker panel saved: %s (%d genes)\n", panel_label, length(genes)))
  }

  summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else data.frame()
  if (!is.null(report_dir) && nzchar(tc_safe_trim(report_dir)) && nrow(summary_df) > 0) {
    tc_write_table_file(summary_df, file.path(report_dir, "marker_panel_visualization_summary.csv"), sep = ",")
  }
  list(summary = summary_df, panels = panel_results, figure_dir = fig_dir)
}

tc_marker_panel_report_lines <- function(summary_df, figure_dir_rel = "figures/marker_panels") {
  if (is.null(summary_df) || !is.data.frame(summary_df) || nrow(summary_df) == 0) return(character())
  lines <- c("### 2.1 Class-based marker panels", "")
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, , drop = FALSE]
    lines <- c(
      lines,
      sprintf("#### %s (%d markers)", row$panel_label, row$n_genes),
      ""
    )
    if (nzchar(tc_safe_trim(row$dotplot_png))) {
      lines <- c(lines, sprintf("![%s dotplot](%s/%s)", row$panel_label, figure_dir_rel, basename(row$dotplot_png)), "")
    }
    if (nzchar(tc_safe_trim(row$featureplot_png))) {
      lines <- c(lines, sprintf("![%s featureplot](%s/%s)", row$panel_label, figure_dir_rel, basename(row$featureplot_png)), "")
    }
    if (nzchar(tc_safe_trim(row$heatmap_png))) {
      lines <- c(lines, sprintf("![%s heatmap](%s/%s)", row$panel_label, figure_dir_rel, basename(row$heatmap_png)), "")
    }
  }
  lines
}

tc_load_env_candidates <- function(env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env")) {
  env_candidates <- unique(path.expand(as.character(env_candidates)))
  env_candidates <- env_candidates[nzchar(env_candidates)]
  loaded <- env_candidates[file.exists(env_candidates)]
  invisible(lapply(loaded, tc_load_env_file))
  loaded
}

tc_assert_named_list_keys <- function(config, required_keys, object_name = "config") {
  if (!is.list(config) || is.null(names(config)) || any(!nzchar(names(config)))) {
    stop(sprintf("%s must be a named list.", object_name))
  }
  missing_keys <- setdiff(required_keys, names(config))
  if (length(missing_keys) > 0) {
    stop(sprintf("%s is missing required keys: %s", object_name, paste(missing_keys, collapse = ", ")))
  }
  invisible(config)
}

tc_apply_named_list <- function(config, envir = parent.frame(), drop_null = TRUE) {
  tc_assert_named_list_keys(config, required_keys = character(), object_name = deparse(substitute(config)))
  for (nm in names(config)) {
    if (isTRUE(drop_null) && is.null(config[[nm]])) next
    assign(nm, config[[nm]], envir = envir)
  }
  invisible(config)
}

tc_preflight_previous_run <- function(previous_output_dir,
                                      current_output_dir = NULL,
                                      allow_in_place_resume = FALSE,
                                      require_previous_dir = FALSE,
                                      reuse_previous_output_summary = TRUE,
                                      summary_output_path = NULL,
                                      include_rds = FALSE,
                                      include_markdown = FALSE,
                                      max_markdown_lines = 200L,
                                      verbose = TRUE) {
  previous_output_dir <- path.expand(tc_safe_trim(previous_output_dir))
  current_output_dir <- if (is.null(current_output_dir)) "" else path.expand(tc_safe_trim(current_output_dir))
  same_dir <- nzchar(current_output_dir) && tc_path_equal(previous_output_dir, current_output_dir)

  if (same_dir && !isTRUE(allow_in_place_resume)) {
    stop(sprintf(
      "PREVIOUS_OUTPUT_DIR and OUTPUT_DIR resolve to the same path (%s). Set ALLOW_IN_PLACE_RESUME=TRUE only if you intentionally want overwrite-prone in-place resume.",
      previous_output_dir
    ))
  }

  if (same_dir && isTRUE(allow_in_place_resume) && isTRUE(verbose)) {
    cat(sprintf("[WARN] In-place resume enabled: previous and current output directories are identical (%s)\n", previous_output_dir))
  }

  if (!dir.exists(previous_output_dir)) {
    if (isTRUE(require_previous_dir)) {
      stop(sprintf("Required PREVIOUS_OUTPUT_DIR does not exist: %s", previous_output_dir))
    }
    if (isTRUE(verbose)) {
      cat(sprintf("[INFO] Previous output directory not found; continue without historical preflight: %s\n", previous_output_dir))
    }
    return(list(
      previous_output_dir = previous_output_dir,
      current_output_dir = if (nzchar(current_output_dir)) current_output_dir else NA_character_,
      same_dir = same_dir,
      exists = FALSE,
      previous_run = NULL,
      summary = data.frame(),
      summary_output_path = summary_output_path
    ))
  }

  previous_run <- if (isTRUE(reuse_previous_output_summary) || !is.null(summary_output_path)) {
    tc_read_previous_run(
      output_dir = previous_output_dir,
      include_rds = include_rds,
      include_markdown = include_markdown,
      max_markdown_lines = max_markdown_lines
    )
  } else {
    NULL
  }
  summary_df <- if (!is.null(previous_run)) previous_run$summary else data.frame()

  if (isTRUE(verbose)) {
    cat(sprintf("[INFO] Previous output preflight: %s\n", previous_output_dir))
    if (!is.null(previous_run)) tc_print_previous_run_summary(previous_run)
  }

  if (!is.null(summary_output_path) && nzchar(tc_safe_trim(summary_output_path)) && !is.null(summary_df) && nrow(summary_df) > 0) {
    tc_write_table_file(summary_df, summary_output_path, sep = "\t")
    if (isTRUE(verbose)) cat(sprintf("[OK] Wrote previous-run summary: %s\n", summary_output_path))
  }

  list(
    previous_output_dir = previous_output_dir,
    current_output_dir = if (nzchar(current_output_dir)) current_output_dir else NA_character_,
    same_dir = same_dir,
    exists = TRUE,
    previous_run = previous_run,
    summary = summary_df,
    summary_output_path = summary_output_path
  )
}

tc_extract_balanced_json_block <- function(text, start_idx) {
  if (start_idx < 1 || start_idx > nchar(text)) return("")
  opener <- substr(text, start_idx, start_idx)
  closer <- if (identical(opener, "[")) "]" else if (identical(opener, "{")) "}" else return("")
  depth <- 0L
  in_string <- FALSE
  escaped <- FALSE
  for (idx in seq.int(start_idx, nchar(text))) {
    ch <- substr(text, idx, idx)
    if (in_string) {
      if (escaped) {
        escaped <- FALSE
      } else if (identical(ch, "\\")) {
        escaped <- TRUE
      } else if (identical(ch, '"')) {
        in_string <- FALSE
      }
      next
    }
    if (identical(ch, '"')) {
      in_string <- TRUE
      next
    }
    if (identical(ch, opener)) depth <- depth + 1L
    if (identical(ch, closer)) {
      depth <- depth - 1L
      if (depth == 0L) return(substr(text, start_idx, idx))
    }
  }
  ""
}

tc_extract_json_string <- function(text) {
  text <- tc_safe_trim(text)
  if (!nzchar(text)) return("")
  text <- gsub("^```(?:json)?\\s*", "", text, perl = TRUE)
  text <- gsub("\\s*```$", "", text, perl = TRUE)
  start_idx <- regexpr("\\[|\\{", text, perl = TRUE)[1]
  if (start_idx < 1) return("")
  tc_extract_balanced_json_block(text, start_idx)
}

tc_parse_json <- function(text) {
  candidates <- unique(c(tc_safe_trim(text), tc_extract_json_string(text)))
  candidates <- candidates[nzchar(candidates)]
  if (!requireNamespace("jsonlite", quietly = TRUE) || length(candidates) == 0) return(NULL)
  for (candidate in candidates) {
    parsed <- tryCatch(jsonlite::fromJSON(candidate, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(parsed)) return(parsed)
  }
  NULL
}

tc_prepare_reasoner_api <- function(
    env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env"),
    model = "deepseek-reasoner",
    write_env_placeholder = FALSE) {
  env_candidates <- env_candidates[file.exists(env_candidates)]
  invisible(lapply(env_candidates, tc_load_env_file))
  key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  if (nchar(key) < 10) {
    if (isTRUE(write_env_placeholder)) {
      tc_ensure_env_placeholder("/home/h2048/.env", "DEEPSEEK_API_KEY")
    }
    stop("DEEPSEEK_API_KEY not set; cannot run reasoner reformat.")
  }
  if (!requireNamespace("fanyi", quietly = TRUE)) {
    stop("Package 'fanyi' is required for reasoner reformat.")
  }
  fanyi::set_translate_option(key = key, source = "deepseek")
  list(api_key = key, model = model, env_files_loaded = env_candidates)
}

tc_pick_table_by_family <- function(previous_run, family) {
  stopifnot(is.list(previous_run), !is.null(previous_run$tables))
  if (identical(family, "pairwise")) {
    pieces <- Filter(function(x) !is.null(x) && nrow(x) > 0, list(
      previous_run$tables$interpret_agent_structured,
      previous_run$tables$interpret_agent_structured_l3
    ))
    if (length(pieces) == 0) return(data.frame())
    return(do.call(rbind, pieces))
  }
  if (identical(family, "ssgsea")) {
    x <- previous_run$tables$ssgsea_llm_structured
    return(if (is.null(x)) data.frame() else x)
  }
  if (identical(family, "choir")) {
    x <- previous_run$tables$choir_llm_structured
    return(if (is.null(x)) data.frame() else x)
  }
  data.frame()
}

tc_reasoner_prompt_for_batch <- function(batch_df, family) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")
  required_cols <- c(
    "record_id", "celltype_level", "celltype_label", "celltype_l2",
    "comparison", "direction", "source_db", "overview", "key_mechanisms",
    "hypothesis", "narrative", "key_drivers", "evidence", "limitations", "raw_text"
  )
  missing_cols <- setdiff(required_cols, names(batch_df))
  if (length(missing_cols) > 0) {
    for (nm in missing_cols) batch_df[[nm]] <- ""
  }
  payload_df <- batch_df[, required_cols, drop = FALSE]
  payload_json <- jsonlite::toJSON(payload_df, dataframe = "rows", auto_unbox = TRUE, pretty = TRUE, na = "null")
  paste(
    "You are reformatting historical single-cell interpretation records into a stricter 2026 schema.",
    "Return valid JSON only. No markdown, no code fences, no commentary.",
    "Return a JSON array with the same number of items and the same record_id values as the input array.",
    "Each output item must contain exactly these keys:",
    "record_id, cell_type_judgment, confidence, overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "Rules:",
    "1. All narrative fields must be concise Simplified Chinese strings.",
    "2. key_drivers must be an array of English gene symbols or regulator names.",
    "3. cell_type_judgment must directly state the most likely cell type/state; if ambiguous, say so.",
    "4. confidence must be one of: high, medium, low.",
    "5. Integrate up and down evidence together; do not discuss them as disconnected conclusions.",
    "6. Use ssGSEA mainly for cell-type/state judgment.",
    "7. hypothesis must be brief and comparative only when comparison context supports it; otherwise use 'Not applicable.'.",
    "8. evidence must be concise but specific, and should cite top DEG / ssGSEA / enrichment clues when present in the input.",
    sprintf("9. %s", tc_noninformative_gene_rule_text()),
    "10. For non-epithelial analyses, do not use cilia/ciliogenesis-related genes (for example FOXJ1, TPPP3, DNAH*, RSPH*, DEUP1, CCNO, MCIDAS, CDC20B) as primary evidence unless the input explicitly indicates an epithelial lineage.",
    "11. Energy metabolism, mitochondrial, ribosomal/translation, ENSG-like placeholder, clone-style non-coding locus, and non-epithelial cilia-related pathways must not be promoted as top pathways when more lineage-informative pathways are available.",
    sprintf("12. These records come from family: %s.", family),
    "13. overview should open with the cell-type/state judgment, then summarize the biological state in 1-2 sentences.",
    "14. Keep every field compact: no separate paragraphs for up/down or per-database evidence unless the field would otherwise become ambiguous.",
    "Input records:",
    payload_json,
    sep = "\n"
  )
}

tc_validate_reasoner_output <- function(out, expected_ids) {
  if (length(out) != length(expected_ids)) {
    stop(sprintf("reasoner returned %d items for %d expected ids", length(out), length(expected_ids)))
  }
  returned_ids <- vapply(out, function(x) tc_to_scalar(x$record_id), character(1))
  if (anyDuplicated(returned_ids)) {
    stop("reasoner returned duplicated record_id values")
  }
  if (!setequal(returned_ids, expected_ids)) {
    stop("reasoner returned mismatched record_id set")
  }
  out[match(expected_ids, returned_ids)]
}

tc_reasoner_reformat_batch <- function(batch_df,
                                       family,
                                       model = "deepseek-reasoner",
                                       max_retries = 3L,
                                       retry_sleep_sec = 2) {
  prompt <- tc_reasoner_prompt_for_batch(batch_df, family)
  last_error <- NULL
  warn_msgs <- character()
  for (attempt in seq_len(max(1L, as.integer(max_retries)))) {
    response_text <- tryCatch(
      tc_safe_trim(fanyi::chat_request(prompt, model = model, api_key = Sys.getenv("DEEPSEEK_API_KEY", unset = ""))),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )
    parsed <- tc_parse_json(response_text)
    if (is.list(parsed) && !is.null(parsed[[1]])) {
      if (!is.null(names(parsed)) && all(c("record_id", "cell_type_judgment") %in% names(parsed))) {
        parsed <- list(parsed)
      }
      out <- lapply(parsed, function(item) {
        if (is.null(item$record_id)) return(NULL)
        list(
          record_id = tc_to_scalar(item$record_id),
          cell_type_judgment = tc_to_scalar(item$cell_type_judgment),
          confidence = tc_to_scalar(item$confidence),
          overview = tc_to_scalar(item$overview),
          key_mechanisms = tc_to_scalar(item$key_mechanisms),
          hypothesis = tc_to_scalar(item$hypothesis),
          narrative = tc_to_scalar(item$narrative),
          key_drivers = if (is.null(item$key_drivers)) character() else as.character(unlist(item$key_drivers, use.names = FALSE)),
          evidence = tc_to_scalar(item$evidence),
          limitations = tc_to_scalar(item$limitations)
        )
      })
      out <- Filter(Negate(is.null), out)
      out <- tryCatch(
        tc_validate_reasoner_output(out, batch_df$record_id),
        error = function(e) {
          warn_msgs <<- c(warn_msgs, sprintf("attempt %d failed validation: %s", attempt, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(out)) return(out)
      warn_msgs <- c(warn_msgs, sprintf("attempt %d returned %d items for %d records", attempt, length(out), nrow(batch_df)))
    } else if (nzchar(response_text)) {
      warn_msgs <- c(warn_msgs, sprintf("attempt %d returned non-JSON text", attempt))
    }
    if (attempt < max_retries) Sys.sleep(retry_sleep_sec)
  }
  stop(sprintf("reasoner reformat failed for family '%s': %s", family, ifelse(is.null(last_error), paste(warn_msgs, collapse = " | "), last_error)))
}

tc_split_batches <- function(n, batch_size = 4L) {
  split(seq_len(n), ceiling(seq_len(n) / max(1L, as.integer(batch_size))))
}

tc_write_reasoner_markdown <- function(records, path, title, include_raw = FALSE) {
  md <- c(title, "", sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
  if (is.null(records) || nrow(records) == 0) {
    md <- c(md, "No reasoner-formatted records available.")
    writeLines(md, path)
    return(invisible(path))
  }
  for (i in seq_len(nrow(records))) {
    rec <- records[i, , drop = FALSE]
    md <- c(
      md,
      sprintf("## %s | %s | %s | %s", rec$celltype_level, rec$celltype_label, rec$comparison, rec$direction),
      "",
      sprintf("- **Status:** %s", rec$status),
      sprintf("- **Source DB:** %s", ifelse(nzchar(rec$source_db), rec$source_db, "NA")),
      sprintf("- **Reasoner Model:** %s", rec$reasoner_model),
      sprintf("- **Confidence:** %s", rec$confidence),
      "",
      "### Cell Type Judgment",
      "",
      tc_to_scalar(rec$cell_type_judgment),
      "",
      "### Overview",
      "",
      tc_to_scalar(rec$overview),
      "",
      "### Key Mechanisms",
      "",
      tc_to_scalar(rec$key_mechanisms),
      "",
      "### Comparative Hypothesis",
      "",
      tc_to_scalar(rec$hypothesis),
      "",
      "### Narrative",
      "",
      tc_to_scalar(rec$narrative),
      "",
      "### Key Drivers",
      "",
      tc_to_scalar(rec$key_drivers),
      "",
      "### Evidence",
      "",
      tc_to_scalar(rec$evidence),
      "",
      "### Limitations",
      "",
      tc_to_scalar(rec$limitations),
      ""
    )
    if (isTRUE(include_raw) && nzchar(tc_to_scalar(rec$raw_text))) {
      md <- c(md, "### Historical Raw Text", "", tc_to_scalar(rec$raw_text), "")
    }
  }
  writeLines(md, path)
  invisible(path)
}

tc_reformat_structured_table_with_reasoner <- function(df,
                                                       family,
                                                       output_dir,
                                                       prefix,
                                                       title,
                                                       model = "deepseek-reasoner",
                                                       batch_size = 4L,
                                                       max_retries = 3L,
                                                       retry_sleep_sec = 2,
                                                       include_raw_markdown = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  tc_prepare_reasoner_api(model = model)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$record_id <- sprintf("%s_%04d", family, seq_len(nrow(df)))

  batches <- tc_split_batches(nrow(df), batch_size = batch_size)
  reformatted <- vector("list", length = nrow(df))
  for (batch_idx in seq_along(batches)) {
    idx <- batches[[batch_idx]]
    cat(sprintf("[Reasoner] %s batch %d/%d (%d records)\n", family, batch_idx, length(batches), length(idx)))
    batch_out <- tryCatch(
      tc_reasoner_reformat_batch(
        batch_df = df[idx, , drop = FALSE],
        family = family,
        model = model,
        max_retries = max_retries,
        retry_sleep_sec = retry_sleep_sec
      ),
      error = function(e) {
        if (length(idx) == 1L) stop(e)
        cat(sprintf("[Reasoner] %s batch %d fallback to single-record retries: %s\n",
                    family, batch_idx, conditionMessage(e)))
        unlist(lapply(idx, function(single_idx) {
          tc_reasoner_reformat_batch(
            batch_df = df[single_idx, , drop = FALSE],
            family = family,
            model = model,
            max_retries = max_retries,
            retry_sleep_sec = retry_sleep_sec
          )
        }), recursive = FALSE)
      }
    )
    batch_map <- setNames(batch_out, vapply(batch_out, function(x) x$record_id, character(1)))
    for (row_idx in idx) {
      rec <- batch_map[[df$record_id[[row_idx]]]]
      reformatted[[row_idx]] <- data.frame(
        celltype_level = tc_to_scalar(df$celltype_level[[row_idx]]),
        celltype_label = tc_to_scalar(df$celltype_label[[row_idx]]),
        celltype_l2 = if ("celltype_l2" %in% names(df)) tc_to_scalar(df$celltype_l2[[row_idx]]) else "",
        comparison = tc_to_scalar(df$comparison[[row_idx]]),
        direction = tc_to_scalar(df$direction[[row_idx]]),
        source_db = tc_to_scalar(df$source_db[[row_idx]]),
        status = "reasoner_reformatted",
        confidence = tc_to_scalar(rec$confidence),
        cell_type_judgment = tc_to_scalar(rec$cell_type_judgment),
        overview = tc_to_scalar(rec$overview),
        key_mechanisms = tc_to_scalar(rec$key_mechanisms),
        hypothesis = tc_to_scalar(rec$hypothesis),
        narrative = tc_to_scalar(rec$narrative),
        key_drivers = paste(unique(rec$key_drivers[nzchar(trimws(rec$key_drivers))]), collapse = ", "),
        evidence = tc_to_scalar(rec$evidence),
        limitations = tc_to_scalar(rec$limitations),
        raw_text = tc_to_scalar(df$raw_text[[row_idx]]),
        reasoner_model = model,
        reformatted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
    }
  }
  out_df <- do.call(rbind, reformatted)
  reports_dir <- file.path(output_dir, "reports")
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  tsv_path <- file.path(reports_dir, sprintf("%s_reasoner_structured.tsv", prefix))
  rds_path <- file.path(reports_dir, sprintf("%s_reasoner_structured.rds", prefix))
  md_path <- file.path(output_dir, sprintf("%s_REASONER.md", toupper(prefix)))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(out_df, tsv_path, sep = "\t")
  } else {
    utils::write.table(out_df, tsv_path, sep = "\t", row.names = FALSE, quote = TRUE)
  }
  saveRDS(out_df, rds_path)
  tc_write_reasoner_markdown(out_df, md_path, title = title, include_raw = include_raw_markdown)
  out_df
}

tc_reformat_previous_run_with_reasoner <- function(output_dir,
                                                   families = c("ssgsea", "choir"),
                                                   model = "deepseek-reasoner",
                                                   batch_size = 4L,
                                                   max_retries = 3L,
                                                   retry_sleep_sec = 2,
                                                   include_raw_markdown = FALSE) {
  previous_run <- tc_read_previous_run(output_dir, include_rds = FALSE, include_markdown = FALSE)
  if (!isTRUE(previous_run$exists)) stop(sprintf("Output dir does not exist: %s", output_dir))

  family_spec <- list(
    pairwise = list(prefix = "llm_interpretation", title = "# LLM Interpretation Reformatted by DeepSeek Reasoner"),
    ssgsea = list(prefix = "llm_ssgsea_interpretation", title = "# LLM ssGSEA Interpretation Reformatted by DeepSeek Reasoner"),
    choir = list(prefix = "llm_choir_interpretation", title = "# LLM CHOIR Interpretation Reformatted by DeepSeek Reasoner")
  )

  outputs <- list()
  for (family in families) {
    if (!family %in% names(family_spec)) next
    df <- tc_pick_table_by_family(previous_run, family)
    if (is.null(df) || nrow(df) == 0) {
      cat(sprintf("[Reasoner] Skip family '%s': no rows found\n", family))
      next
    }
    outputs[[family]] <- tc_reformat_structured_table_with_reasoner(
      df = df,
      family = family,
      output_dir = output_dir,
      prefix = family_spec[[family]]$prefix,
      title = family_spec[[family]]$title,
      model = model,
      batch_size = batch_size,
      max_retries = max_retries,
      retry_sleep_sec = retry_sleep_sec,
      include_raw_markdown = include_raw_markdown
    )
  }
  outputs
}

tc_apply_advanced_shared_overrides <- function(envir = parent.frame(), overrides = list()) {
  defaults <- list(
    INTERPRET_MULTI_DB_MIN_TERMS = 1L,
    INTERPRET_MULTI_DB_MAX_DBS = 5L,
    INTERPRET_MULTI_DB_TERMS_PER_DB = 4L,
    LLM_INCLUDE_TOP_DEG = TRUE,
    LLM_TOP_DEG_N = 10L,
    LLM_DEG_PADJ_THR = 0.10,
    LLM_DEG_LFC_THR = 0.15,
    LLM_SSGSEA_TERMS_PER_DIRECTION = 6L,
    LLM_REQUIRE_INTEGRATED_UP_DOWN = TRUE,
    LLM_SSGSEA_PRIMARY_USE = "cell_type_judgment",
    LLM_ALLOW_COMPARATIVE_HYPOTHESIS = TRUE,
    LLM_EXTRA_RULES = c(
      "Start overview with a concise cell-type/state judgment.",
      "Use ssGSEA mainly to support cell-type judgment.",
      "Integrate up/down evidence into one compact interpretation.",
      "Treat clone-style non-coding loci (such as AC/AL/AP/BX/CTB/CTC/CTD/LOC/*-OT) as non-informative unless no stronger evidence exists."
    )
  )
  merged <- utils::modifyList(defaults, overrides)
  for (nm in names(merged)) assign(nm, merged[[nm]], envir = envir)
  invisible(merged)
}

tc_detect_output_paths <- function(output_dir) {
  output_dir <- path.expand(output_dir)
  reports_dir <- file.path(output_dir, "reports")
  choir_dir <- file.path(reports_dir, "choir")
  list(
    output_dir = output_dir,
    reports_dir = reports_dir,
    figures_dir = file.path(output_dir, "figures"),
    choir_dir = choir_dir,
    markdown = list(
      report = file.path(output_dir, "REPORT.md"),
      llm_main = file.path(output_dir, "LLM_INTERPRETATION.md"),
      llm_ssgsea = file.path(output_dir, "LLM_SSGSEA_INTERPRETATION.md"),
      llm_choir = file.path(output_dir, "LLM_CHOIR_INTERPRETATION.md")
    ),
    tables = list(
      interpret_agent_structured = file.path(reports_dir, "interpret_agent_structured.tsv"),
      interpret_agent_structured_l3 = file.path(reports_dir, "interpret_agent_structured_L3.tsv"),
      ssgsea_llm_structured = file.path(reports_dir, "ssgsea_llm_structured.tsv"),
      choir_llm_structured = file.path(reports_dir, "choir_llm_structured.tsv"),
      l3_l2_mapping_summary = file.path(reports_dir, "L3_to_L2_mapping_summary.csv"),
      choir_cluster_l2_counts = file.path(choir_dir, "choir_cluster_L2_counts.csv"),
      choir_cluster_l3_counts = file.path(choir_dir, "choir_cluster_L3_counts.csv"),
      choir_cluster_l2_l3_counts_long = file.path(choir_dir, "choir_cluster_L2_L3_counts_long.csv")
    ),
    rds = list(
      interpret_agent_all = file.path(reports_dir, "interpret_agent_all.rds"),
      interpret_agent_structured_all = file.path(reports_dir, "interpret_agent_structured_all.rds"),
      ssgsea_results_all = file.path(reports_dir, "ssgsea_results_all.rds"),
      ssgsea_results_l3_all = file.path(reports_dir, "ssgsea_results_L3_all.rds"),
      choir_results_all = file.path(choir_dir, "choir_results_all.rds"),
      choir_ofa_all = file.path(choir_dir, "choir_ofa_all.rds")
    )
  )
}

tc_collect_standard_tables <- function(output_dir) {
  paths <- tc_detect_output_paths(output_dir)
  reports_dir <- paths$reports_dir
  choir_dir <- paths$choir_dir
  table_candidates <- list(
    interpret_agent_structured = c(
      paths$tables$interpret_agent_structured,
      file.path(reports_dir, "interpret_agent_structured_all.rds"),
      file.path(reports_dir, "interpret_agent_all.rds")
    ),
    interpret_agent_structured_l3 = c(
      paths$tables$interpret_agent_structured_l3,
      file.path(reports_dir, "interpret_agent_structured_L3_all.rds"),
      file.path(reports_dir, "interpret_agent_L3_all.rds")
    ),
    ssgsea_llm_structured = c(
      paths$tables$ssgsea_llm_structured,
      file.path(reports_dir, "ssgsea_llm_structured_all.rds"),
      file.path(reports_dir, "ssgsea_llm_structured_L2_all.rds"),
      file.path(reports_dir, "ssgsea_llm_structured_L3_all.rds")
    ),
    choir_llm_structured = c(
      paths$tables$choir_llm_structured,
      file.path(reports_dir, "choir_llm_structured_all.rds"),
      file.path(choir_dir, "choir_llm_structured_all.rds")
    ),
    l3_l2_mapping_summary = c(paths$tables$l3_l2_mapping_summary),
    choir_cluster_l2_counts = c(paths$tables$choir_cluster_l2_counts),
    choir_cluster_l3_counts = c(paths$tables$choir_cluster_l3_counts),
    choir_cluster_l2_l3_counts_long = c(paths$tables$choir_cluster_l2_l3_counts_long)
  )
  out <- lapply(table_candidates, tc_read_previous_artifacts, artifact_type = "table")
  out[!vapply(out, function(x) is.null(x) || (is.data.frame(x) && nrow(x) == 0), logical(1))]
}

tc_collect_markdown_outputs <- function(output_dir, max_lines = 200L) {
  md_paths <- tc_detect_output_paths(output_dir)$markdown
  out <- lapply(md_paths, function(path) {
    resolved <- tc_first_existing_path(path)
    if (is.na(resolved) || !nzchar(resolved)) return(character())
    tc_safe_read_text(resolved, max_lines = max_lines)
  })
  out[vapply(out, length, integer(1)) > 0]
}

tc_collect_ssgsea_top_tables <- function(output_dir) {
  paths <- tc_detect_output_paths(output_dir)
  report_files <- tc_list_files(paths$reports_dir, pattern = "^ssgsea(_l3)?_top_pathways_.*\\.(csv|tsv)$", recursive = FALSE)
  choir_files <- tc_list_files(paths$choir_dir, pattern = "^ssgsea_top_pathways_.*\\.(csv|tsv)$", recursive = FALSE)
  all_files <- c(report_files, choir_files)
  if (length(all_files) == 0) return(list())
  names(all_files) <- basename(all_files)
  out <- lapply(all_files, tc_safe_read_table)
  out[!vapply(out, is.null, logical(1))]
}

tc_collect_choir_marker_tables <- function(output_dir) {
  choir_dir <- tc_detect_output_paths(output_dir)$choir_dir
  marker_files <- tc_list_files(choir_dir, pattern = "markers\\.(csv|tsv)$", recursive = TRUE)
  if (length(marker_files) == 0) return(list())
  names(marker_files) <- tc_strip_path_prefix(marker_files, choir_dir)
  out <- lapply(marker_files, tc_safe_read_table)
  out[!vapply(out, is.null, logical(1))]
}

tc_collect_rds_outputs <- function(output_dir) {
  rds_paths <- tc_detect_output_paths(output_dir)$rds
  rds_candidates <- list(
    interpret_agent_all = c(rds_paths$interpret_agent_all, file.path(output_dir, "reports", "interpret_agent_structured_all.rds")),
    interpret_agent_structured_all = c(rds_paths$interpret_agent_structured_all, file.path(output_dir, "reports", "interpret_agent_all.rds")),
    ssgsea_results_all = c(rds_paths$ssgsea_results_all),
    ssgsea_results_l3_all = c(rds_paths$ssgsea_results_l3_all, file.path(output_dir, "reports", "ssgsea_results_l3_all.rds")),
    choir_results_all = c(rds_paths$choir_results_all, file.path(output_dir, "reports", "choir_results_all.rds")),
    choir_ofa_all = c(rds_paths$choir_ofa_all, file.path(output_dir, "reports", "choir_ofa_all.rds"))
  )
  out <- lapply(rds_candidates, tc_read_previous_artifacts, artifact_type = "rds")
  out[!vapply(out, is.null, logical(1))]
}

tc_file_manifest <- function(output_dir) {
  if (!dir.exists(output_dir)) return(data.frame())
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (length(files) == 0) return(data.frame())
  info <- file.info(files)
  data.frame(
    file = tc_strip_path_prefix(files, output_dir),
    size_bytes = unname(info$size),
    modified = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
}

tc_summarize_previous_run <- function(x) {
  if (is.character(x) && length(x) == 1) x <- tc_read_previous_run(x)
  if (!is.list(x) || is.null(x$output_dir)) stop("Input must be an output_dir path or tc_read_previous_run() result.")

  tables <- x$tables
  summary_list <- list(
    output_dir = x$output_dir,
    exists = isTRUE(x$exists),
    n_files = if (is.null(x$manifest) || nrow(x$manifest) == 0) 0L else nrow(x$manifest),
    has_report_md = "report" %in% names(x$markdown),
    has_llm_main = "llm_main" %in% names(x$markdown),
    has_llm_ssgsea = "llm_ssgsea" %in% names(x$markdown),
    has_llm_choir = "llm_choir" %in% names(x$markdown),
    n_pairwise_llm_rows = if (!is.null(tables$interpret_agent_structured)) nrow(tables$interpret_agent_structured) else 0L,
    n_pairwise_l3_llm_rows = if (!is.null(tables$interpret_agent_structured_l3)) nrow(tables$interpret_agent_structured_l3) else 0L,
    n_ssgsea_llm_rows = if (!is.null(tables$ssgsea_llm_structured)) nrow(tables$ssgsea_llm_structured) else 0L,
    n_choir_llm_rows = if (!is.null(tables$choir_llm_structured)) nrow(tables$choir_llm_structured) else 0L,
    n_choir_annotation_tables = sum(vapply(c(
      "choir_cluster_l2_counts",
      "choir_cluster_l3_counts",
      "choir_cluster_l2_l3_counts_long"
    ), function(nm) !is.null(tables[[nm]]) && nrow(tables[[nm]]) > 0, logical(1))),
    n_ssgsea_top_tables = length(x$ssgsea_top_tables),
    n_choir_marker_tables = length(x$choir_marker_tables),
    n_rds_loaded = length(x$rds)
  )
  as.data.frame(summary_list, stringsAsFactors = FALSE)
}

tc_print_previous_run_summary <- function(x) {
  summary_df <- tc_summarize_previous_run(x)
  cat(sprintf("Output dir: %s\n", summary_df$output_dir[[1]]))
  cat(sprintf("Exists: %s\n", summary_df$exists[[1]]))
  cat(sprintf("Files indexed: %d\n", summary_df$n_files[[1]]))
  cat(sprintf("Pairwise LLM rows: %d\n", summary_df$n_pairwise_llm_rows[[1]]))
  cat(sprintf("Pairwise L3 LLM rows: %d\n", summary_df$n_pairwise_l3_llm_rows[[1]]))
  cat(sprintf("ssGSEA LLM rows: %d\n", summary_df$n_ssgsea_llm_rows[[1]]))
  cat(sprintf("CHOIR LLM rows: %d\n", summary_df$n_choir_llm_rows[[1]]))
  cat(sprintf("CHOIR annotation tables: %d\n", summary_df$n_choir_annotation_tables[[1]]))
  cat(sprintf("ssGSEA top tables: %d\n", summary_df$n_ssgsea_top_tables[[1]]))
  cat(sprintf("CHOIR marker tables: %d\n", summary_df$n_choir_marker_tables[[1]]))
  invisible(summary_df)
}

tc_read_previous_run <- function(output_dir,
                                 include_rds = FALSE,
                                 include_markdown = TRUE,
                                 max_markdown_lines = 200L) {
  output_dir <- path.expand(output_dir)
  exists_flag <- dir.exists(output_dir)
  result <- list(
    output_dir = output_dir,
    exists = exists_flag,
    manifest = if (exists_flag) tc_file_manifest(output_dir) else data.frame(),
    tables = if (exists_flag) tc_collect_standard_tables(output_dir) else list(),
    ssgsea_top_tables = if (exists_flag) tc_collect_ssgsea_top_tables(output_dir) else list(),
    choir_marker_tables = if (exists_flag) tc_collect_choir_marker_tables(output_dir) else list(),
    markdown = if (exists_flag && isTRUE(include_markdown)) tc_collect_markdown_outputs(output_dir, max_lines = max_markdown_lines) else list(),
    rds = if (exists_flag && isTRUE(include_rds)) tc_collect_rds_outputs(output_dir) else list()
  )
  result$summary <- tc_summarize_previous_run(result)
  class(result) <- c("tc_previous_run", class(result))
  result
}

tc_find_previous_runs <- function(parent_dir,
                                  pattern = "tissue_comparison",
                                  recursive = TRUE) {
  parent_dir <- path.expand(parent_dir)
  if (!dir.exists(parent_dir)) return(data.frame())
  all_dirs <- list.dirs(parent_dir, recursive = recursive, full.names = TRUE)
  all_dirs <- unique(c(parent_dir, all_dirs))
  keep <- vapply(all_dirs, function(dir_path) {
    basename_ok <- grepl(pattern, basename(dir_path), ignore.case = TRUE)
    report_ok <- file.exists(file.path(dir_path, "reports")) || file.exists(file.path(dir_path, "REPORT.md"))
    basename_ok && report_ok
  }, logical(1))
  hit_dirs <- all_dirs[keep]
  if (length(hit_dirs) == 0) return(data.frame())
  info <- file.info(hit_dirs)
  out <- data.frame(
    output_dir = hit_dirs,
    modified = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  out[order(info$mtime, decreasing = TRUE), , drop = FALSE]
}

tc_read_latest_previous_run <- function(parent_dir,
                                        pattern = "tissue_comparison",
                                        recursive = TRUE,
                                        include_rds = FALSE,
                                        include_markdown = TRUE,
                                        max_markdown_lines = 200L) {
  runs <- tc_find_previous_runs(parent_dir = parent_dir, pattern = pattern, recursive = recursive)
  if (nrow(runs) == 0) return(NULL)
  tc_read_previous_run(
    runs$output_dir[[1]],
    include_rds = include_rds,
    include_markdown = include_markdown,
    max_markdown_lines = max_markdown_lines
  )
}

if (sys.nframe() == 0) {
  cat("Tissue Comparison Advanced Helper loaded.\n")
  cat("Use tc_apply_advanced_shared_overrides() before sourcing shared engines.\n")
  cat("Use tc_read_previous_run(output_dir) to inspect historical results.\n")
}