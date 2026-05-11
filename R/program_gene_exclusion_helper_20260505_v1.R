#!/usr/bin/env Rscript
# ==============================================================================
# Program Gene Exclusion Helper (2026-05-05 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)

PA_GENE_EXCLUSION_HELPER_VERSION_20260505_V1 <- "20260505_v1"
PA_GENE_EXCLUSION_SPEC_PATH_20260505_V1 <- "/home/h2048/script/config/program_gene_exclusion_spec_20260505_v1.json"

pa_as_logical <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) return(isTRUE(default))
  isTRUE(x[[1]])
}

pa_upper_chr <- function(x) {
  out <- pa_safe_trim(x)
  toupper(out)
}

pa_match_prefix_any <- function(values, prefixes) {
  values <- pa_upper_chr(values)
  prefixes <- pa_upper_chr(prefixes)
  if (length(values) == 0L || length(prefixes) == 0L) return(rep(FALSE, length(values)))
  Reduce(`|`, lapply(prefixes, function(prefix) startsWith(values, prefix)), init = rep(FALSE, length(values)))
}

pa_match_regex_any <- function(values, patterns) {
  values <- pa_safe_trim(values)
  patterns <- pa_safe_trim(patterns)
  if (length(values) == 0L || length(patterns) == 0L) return(rep(FALSE, length(values)))
  Reduce(`|`, lapply(patterns, function(pattern) grepl(pattern, values, perl = TRUE, ignore.case = TRUE)), init = rep(FALSE, length(values)))
}

pa_load_gene_exclusion_spec <- function(spec_path = PA_GENE_EXCLUSION_SPEC_PATH_20260505_V1,
                                        overrides = list()) {
  spec_path <- pa_scalar_chr(spec_path, "spec_path")
  spec <- pa_json_read(spec_path, default = NULL)
  if (is.null(spec)) {
    stop(sprintf("Gene exclusion spec does not exist or could not be parsed: %s", spec_path), call. = FALSE)
  }
  if (!is.list(overrides)) overrides <- list()
  utils::modifyList(spec, overrides)
}

pa_normalize_annotation_df <- function(annotation_df,
                                       spec) {
  if (is.null(annotation_df) || !is.data.frame(annotation_df) || nrow(annotation_df) == 0L) return(NULL)
  gene_candidates <- pa_safe_trim(pa_null_coalesce(spec$annotation_gene_column_candidates, character()))
  biotype_candidates <- pa_safe_trim(pa_null_coalesce(spec$annotation_column_candidates, character()))

  gene_col <- gene_candidates[gene_candidates %in% colnames(annotation_df)][1]
  biotype_col <- biotype_candidates[biotype_candidates %in% colnames(annotation_df)][1]
  if (is.na(gene_col) || !nzchar(gene_col) || is.na(biotype_col) || !nzchar(biotype_col)) return(NULL)

  out <- data.frame(
    gene_symbol = pa_safe_trim(annotation_df[[gene_col]]),
    annotation_biotype = pa_safe_trim(annotation_df[[biotype_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[nzchar(out$gene_symbol), , drop = FALSE]
  if (nrow(out) == 0L) return(NULL)
  out$gene_symbol_upper <- pa_upper_chr(out$gene_symbol)
  out <- out[!duplicated(out$gene_symbol_upper), , drop = FALSE]
  rownames(out) <- out$gene_symbol_upper
  out
}

pa_read_annotation_df <- function(annotation_path, spec) {
  if (is.null(annotation_path) || length(annotation_path) == 0L) return(NULL)
  annotation_path <- pa_safe_trim(annotation_path)[1]
  if (!nzchar(annotation_path) || !file.exists(annotation_path)) return(NULL)
  df <- tryCatch(pa_read_table_auto(annotation_path), error = function(e) NULL)
  pa_normalize_annotation_df(df, spec)
}

pa_infer_lineage_context_from_seurat <- function(seurat_obj,
                                                 preferred_meta_cols = c("lineage", "lineage_tag", "cell_type_L2", "cell_type_L3", "cell_type_final_l2", "cell_type_final_l3")) {
  ctx <- character()
  meta <- tryCatch(seurat_obj@meta.data, error = function(e) NULL)
  if (is.null(meta) || !is.data.frame(meta) || nrow(meta) == 0L) return(NA_character_)
  for (col in preferred_meta_cols) {
    if (!col %in% colnames(meta)) next
    vals <- pa_unique_chr(meta[[col]])
    if (length(vals) == 0L) next
    ctx <- c(ctx, vals)
  }
  ctx <- unique(ctx)
  if (length(ctx) == 0L) return(NA_character_)
  paste(head(ctx, 25L), collapse = " | ")
}

pa_lineage_context_keeps_ig <- function(lineage_context,
                                        spec) {
  policy <- tolower(pa_safe_trim(pa_null_coalesce(spec$ig_policy, "auto"))[1])
  if (identical(policy, "keep")) return(TRUE)
  if (identical(policy, "exclude")) return(FALSE)
  keep_patterns <- pa_null_coalesce(spec$ig_keep_lineage_patterns, character())
  if (is.null(lineage_context) || length(lineage_context) == 0L || is.na(lineage_context[[1]]) || !nzchar(pa_safe_trim(lineage_context[[1]]))) {
    return(FALSE)
  }
  any(grepl(paste(keep_patterns, collapse = "|"), pa_upper_chr(lineage_context[[1]]), perl = TRUE, ignore.case = TRUE))
}

pa_build_gene_exclusion_packet <- function(feature_names,
                                           spec = pa_load_gene_exclusion_spec(),
                                           lineage_context = NA_character_,
                                           annotation_df = NULL,
                                           annotation_path = NULL) {
  feature_names <- pa_safe_trim(feature_names)
  feature_names <- feature_names[nzchar(feature_names)]
  if (length(feature_names) == 0L) {
    stop("feature_names must contain at least one non-empty gene symbol", call. = FALSE)
  }

  annotation_tbl <- pa_normalize_annotation_df(annotation_df, spec)
  if (is.null(annotation_tbl)) {
    annotation_tbl <- pa_read_annotation_df(annotation_path, spec)
  }

  genes_upper <- pa_upper_chr(feature_names)
  annotation_biotype <- rep(NA_character_, length(feature_names))
  if (!is.null(annotation_tbl)) {
    idx <- match(genes_upper, rownames(annotation_tbl))
    hit <- !is.na(idx)
    annotation_biotype[hit] <- annotation_tbl$annotation_biotype[idx[hit]]
  }

  mt_match <- pa_match_prefix_any(feature_names, pa_null_coalesce(spec$mitochondrial_prefixes, character()))
  ribo_match <- pa_match_prefix_any(feature_names, pa_null_coalesce(spec$ribosomal_prefixes, character()))
  ig_match <- pa_match_regex_any(feature_names, pa_null_coalesce(spec$ig_patterns, character()))

  lnc_biotypes <- tolower(pa_safe_trim(pa_null_coalesce(spec$lncrna_biotypes, character())))
  lnc_annotation_match <- !is.na(annotation_biotype) & tolower(annotation_biotype) %in% lnc_biotypes
  lnc_heuristic_match <- pa_match_regex_any(feature_names, pa_null_coalesce(spec$lncrna_symbol_patterns, character()))
  lnc_match <- lnc_annotation_match | lnc_heuristic_match

  keep_ig <- pa_lineage_context_keeps_ig(lineage_context, spec)

  exclude_mt <- pa_as_logical(spec$exclude_mt, TRUE)
  exclude_ribo <- pa_as_logical(spec$exclude_ribo, TRUE)
  exclude_lnc <- pa_as_logical(spec$exclude_lncRNA, TRUE)

  should_exclude <- rep(FALSE, length(feature_names))
  exclude_reason <- rep("kept", length(feature_names))
  rule_source <- rep("none", length(feature_names))

  if (exclude_mt) {
    should_exclude[mt_match] <- TRUE
    exclude_reason[mt_match] <- "mitochondrial"
    rule_source[mt_match] <- "heuristic"
  }
  if (exclude_ribo) {
    update_idx <- ribo_match & !should_exclude
    should_exclude[update_idx] <- TRUE
    exclude_reason[update_idx] <- "ribosomal"
    rule_source[update_idx] <- "heuristic"
  }
  ig_exclude <- ig_match & !keep_ig
  update_idx <- ig_exclude & !should_exclude
  should_exclude[update_idx] <- TRUE
  exclude_reason[update_idx] <- "immunoglobulin"
  rule_source[update_idx] <- "heuristic"
  kept_ig_idx <- ig_match & keep_ig & !should_exclude
  exclude_reason[kept_ig_idx] <- "kept_immunoglobulin_due_to_lineage"
  rule_source[kept_ig_idx] <- "heuristic"

  lnc_exclude <- lnc_match & exclude_lnc
  update_idx <- lnc_exclude & !should_exclude
  should_exclude[update_idx] <- TRUE
  exclude_reason[update_idx] <- "lncrna"
  rule_source[update_idx] <- ifelse(lnc_annotation_match[update_idx], "annotation", "heuristic")

  matched_categories <- vapply(seq_along(feature_names), function(i) {
    cats <- character()
    if (mt_match[[i]]) cats <- c(cats, "mitochondrial")
    if (ribo_match[[i]]) cats <- c(cats, "ribosomal")
    if (ig_match[[i]]) cats <- c(cats, "immunoglobulin")
    if (lnc_match[[i]]) cats <- c(cats, "lncrna")
    if (length(cats) == 0L) cats <- "other"
    paste(cats, collapse = "|")
  }, character(1))

  primary_category <- ifelse(mt_match, "mitochondrial",
    ifelse(ribo_match, "ribosomal",
      ifelse(ig_match, "immunoglobulin",
        ifelse(lnc_match, "lncrna", "other")
      )
    )
  )

  audit_table <- data.frame(
    gene_symbol = feature_names,
    gene_symbol_upper = genes_upper,
    primary_category = primary_category,
    matched_categories = matched_categories,
    should_exclude = should_exclude,
    exclude_reason = exclude_reason,
    rule_source = rule_source,
    annotation_biotype = annotation_biotype,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  summary_categories <- c("mitochondrial", "ribosomal", "immunoglobulin", "lncrna", "other")
  max_examples <- suppressWarnings(as.integer(pa_null_coalesce(spec$max_example_genes, 12L)))
  if (is.na(max_examples) || max_examples < 1L) max_examples <- 12L
  summary_table <- do.call(rbind, lapply(summary_categories, function(cat_name) {
    idx <- audit_table$primary_category == cat_name
    exclude_idx <- idx & audit_table$should_exclude
    data.frame(
      category = cat_name,
      detected_n = sum(idx),
      excluded_n = sum(exclude_idx),
      kept_n = sum(idx & !audit_table$should_exclude),
      pct_of_all_features = round(100 * sum(idx) / nrow(audit_table), 3),
      pct_of_all_excluded = round(100 * sum(exclude_idx) / max(1, sum(audit_table$should_exclude)), 3),
      example_genes = paste(head(audit_table$gene_symbol[idx], max_examples), collapse = ";"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  keep_genes <- audit_table$gene_symbol[!audit_table$should_exclude]
  exclude_genes <- audit_table$gene_symbol[audit_table$should_exclude]

  manifest <- list(
    spec_version = pa_null_coalesce(spec$spec_version, PA_GENE_EXCLUSION_HELPER_VERSION_20260505_V1),
    lineage_context = pa_null_coalesce(lineage_context, NA_character_),
    ig_policy = pa_null_coalesce(spec$ig_policy, "auto"),
    keep_ig = keep_ig,
    n_total_features = nrow(audit_table),
    n_excluded_features = length(exclude_genes),
    n_kept_features = length(keep_genes),
    annotation_used = !is.null(annotation_tbl)
  )

  list(
    audit_table = audit_table,
    summary_table = summary_table,
    keep_genes = keep_genes,
    exclude_genes = exclude_genes,
    manifest = manifest
  )
}

pa_plot_gene_exclusion_summary <- function(summary_table,
                                           output_prefix,
                                           title = "Gene exclusion summary") {
  summary_table <- as.data.frame(summary_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_table) == 0L) return(invisible(NULL))
  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)
  mat <- t(as.matrix(summary_table[, c("detected_n", "excluded_n", "kept_n"), drop = FALSE]))
  colnames(mat) <- summary_table$category
  rownames(mat) <- c("Detected", "Excluded", "Kept")

  draw_once <- function(device_fun) {
    device_fun()
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mar = c(7, 6, 4, 2) + 0.1)
    cols <- c("#9ECAE1", "#FC9272", "#A1D99B")
    bp <- barplot(mat,
                  beside = TRUE,
                  col = cols,
                  border = NA,
                  las = 2,
                  cex.names = 0.9,
                  main = title,
                  ylab = "Gene count",
                  ylim = c(0, max(1, mat, na.rm = TRUE) * 1.18))
    legend("topright", legend = rownames(mat), fill = cols, bty = "n", cex = 0.9)
    text(x = as.vector(bp), y = as.vector(mat), labels = as.vector(mat), pos = 3, cex = 0.75)
    grDevices::dev.off()
  }

  draw_once(function() grDevices::pdf(pdf_path, width = 10, height = 7))
  draw_once(function() grDevices::png(png_path, width = 1800, height = 1200, res = 180))
  invisible(list(pdf = pdf_path, png = png_path))
}

pa_plot_gene_exclusion_dashboard <- function(packet,
                                             output_prefix,
                                             title = "Gene exclusion QC dashboard") {
  if (is.null(packet) || !is.list(packet) || is.null(packet$summary_table)) return(invisible(NULL))
  summary_table <- as.data.frame(packet$summary_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(summary_table) == 0L) return(invisible(NULL))
  pdf_path <- sprintf("%s.pdf", output_prefix)
  png_path <- sprintf("%s.png", output_prefix)

  detected <- suppressWarnings(as.numeric(summary_table$detected_n))
  excluded <- suppressWarnings(as.numeric(summary_table$excluded_n))
  kept <- suppressWarnings(as.numeric(summary_table$kept_n))
  detected[is.na(detected)] <- 0
  excluded[is.na(excluded)] <- 0
  kept[is.na(kept)] <- 0
  pct_excluded <- ifelse(detected > 0, 100 * excluded / detected, 0)
  names(detected) <- names(excluded) <- names(kept) <- names(pct_excluded) <- summary_table$category

  manifest <- packet$manifest
  audit <- as.data.frame(packet$audit_table, stringsAsFactors = FALSE, check.names = FALSE)
  top_excluded <- audit[audit$should_exclude, , drop = FALSE]
  if (nrow(top_excluded) > 0L) {
    top_excluded <- head(top_excluded[, c("gene_symbol", "exclude_reason", "rule_source"), drop = FALSE], 12L)
    example_lines <- apply(top_excluded, 1, function(row) paste(row, collapse = " / "))
  } else {
    example_lines <- "No excluded genes in this feature set."
  }

  draw_once <- function(device_fun) {
    device_fun()
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), widths = c(1.1, 1), heights = c(1, 1))

    par(mar = c(7, 5, 4, 1) + 0.1)
    mat_keep_exclude <- rbind(Excluded = excluded, Kept = kept)
    barplot(mat_keep_exclude,
            beside = FALSE,
            col = c("#FC9272", "#A1D99B"),
            border = NA,
            las = 2,
            main = "Excluded vs kept by category",
            ylab = "Gene count",
            ylim = c(0, max(1, detected, na.rm = TRUE) * 1.12))
    legend("topright", legend = rownames(mat_keep_exclude), fill = c("#FC9272", "#A1D99B"), bty = "n", cex = 0.8)

    par(mar = c(7, 5, 4, 1) + 0.1)
    bp <- barplot(pct_excluded,
                  col = "#9ECAE1",
                  border = NA,
                  las = 2,
                  main = "Excluded fraction within category",
                  ylab = "% excluded",
                  ylim = c(0, max(100, pct_excluded, na.rm = TRUE) * 1.05))
    text(x = bp, y = pct_excluded, labels = sprintf("%.1f%%", pct_excluded), pos = 3, cex = 0.75)

    par(mar = c(7, 5, 4, 1) + 0.1)
    bp2 <- barplot(detected,
                   col = "#BCBDDC",
                   border = NA,
                   las = 2,
                   main = "Detected technical categories",
                   ylab = "Detected genes",
                   ylim = c(0, max(1, detected, na.rm = TRUE) * 1.18))
    text(x = bp2, y = detected, labels = detected, pos = 3, cex = 0.75)

    par(mar = c(1, 1, 4, 1) + 0.1)
    plot.new()
    text_lines <- c(
      title,
      "",
      sprintf("Total features: %s", pa_null_coalesce(manifest$n_total_features, NA_integer_)),
      sprintf("Excluded: %s", pa_null_coalesce(manifest$n_excluded_features, NA_integer_)),
      sprintf("Kept: %s", pa_null_coalesce(manifest$n_kept_features, NA_integer_)),
      sprintf("Lineage context: %s", pa_null_coalesce(manifest$lineage_context, "NA")),
      sprintf("IG policy: %s; keep IG: %s", pa_null_coalesce(manifest$ig_policy, "NA"), isTRUE(manifest$keep_ig)),
      sprintf("Annotation used: %s", isTRUE(manifest$annotation_used)),
      "",
      "Example excluded genes:",
      paste0("- ", example_lines)
    )
    text(0, 1, paste(text_lines, collapse = "\n"), adj = c(0, 1), cex = 0.78, family = "mono")
    title(main = "Audit notes", line = 1)
    grDevices::dev.off()
  }

  draw_once(function() grDevices::pdf(pdf_path, width = 13, height = 9))
  draw_once(function() grDevices::png(png_path, width = 2340, height = 1620, res = 180))
  invisible(list(pdf = pdf_path, png = png_path))
}

pa_plot_gene_exclusion_visualizations <- function(packet,
                                                  output_dir,
                                                  prefix = "gene_exclusion") {
  output_dir <- pa_prepare_output_dir(output_dir)
  prefix <- pa_scalar_chr(prefix, "prefix")
  title_context <- pa_null_coalesce(packet$manifest$lineage_context, "lineage:auto")
  summary_paths <- pa_plot_gene_exclusion_summary(
    summary_table = packet$summary_table,
    output_prefix = file.path(output_dir, sprintf("%s_summary", prefix)),
    title = sprintf("Gene exclusion summary (%s)", title_context)
  )
  dashboard_paths <- pa_plot_gene_exclusion_dashboard(
    packet = packet,
    output_prefix = file.path(output_dir, sprintf("%s_dashboard", prefix)),
    title = sprintf("Gene exclusion dashboard (%s)", title_context)
  )
  manifest <- list(
    summary_plot = summary_paths,
    dashboard_plot = dashboard_paths
  )
  pa_write_json(manifest, file.path(output_dir, sprintf("%s_visualization_manifest.json", prefix)))
  manifest$visualization_manifest_json <- file.path(output_dir, sprintf("%s_visualization_manifest.json", prefix))
  manifest
}

pa_ge_is_placeholder_secret <- function(x) {
  x <- pa_safe_trim(x)
  if (length(x) == 0L || !nzchar(x[[1]])) return(TRUE)
  lowered <- tolower(x[[1]])
  lowered %in% c(
    "your-deepseek-api-key",
    "your_deepseek_api_key_here",
    "your_deepseek_api_key",
    "your-key",
    "replace_me",
    "changeme"
  ) || grepl("^your[-_a-z]*api[-_a-z]*key", lowered) || grepl("placeholder|replace|changeme", lowered)
}

pa_ge_load_env_file <- function(path, overwrite_placeholder = TRUE) {
  if (is.null(path) || length(path) == 0L || !file.exists(path[[1]])) return(invisible(FALSE))
  lines <- readLines(path[[1]], warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- trimws(gsub("^export\\s+", "", sub("=.*$", "", line)))
    val <- trimws(gsub("^['\"]|['\"]$", "", sub("^[^=]*=", "", line)))
    if (!nzchar(key)) next
    cur <- Sys.getenv(key, unset = "")
    should_set <- !nzchar(cur)
    if (!should_set && isTRUE(overwrite_placeholder)) should_set <- pa_ge_is_placeholder_secret(cur)
    if (should_set) do.call(Sys.setenv, stats::setNames(list(val), key))
  }
  invisible(TRUE)
}

pa_ge_ensure_env_placeholder <- function(path = "/home/h2048/.env", key = "DEEPSEEK_API_KEY") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  line <- sprintf("%s=your_deepseek_api_key_here", key)
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    if (any(grepl(sprintf("^%s\\s*=", key), trimws(lines)))) return(invisible(FALSE))
    write(line, file = path, append = TRUE)
    return(invisible(TRUE))
  }
  writeLines(line, path)
  invisible(TRUE)
}

pa_ge_resolve_deepseek_key <- function(api_key = NULL,
                                       env_candidates = c("/home/h2048/.env", "/home/h2048/script/.env"),
                                       write_env_placeholder = FALSE) {
  env_candidates <- as.character(env_candidates)
  env_candidates <- env_candidates[file.exists(env_candidates)]
  invisible(lapply(env_candidates, pa_ge_load_env_file))
  key <- if (is.null(api_key) || length(api_key) == 0L || !nzchar(pa_safe_trim(api_key)[[1]])) {
    Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  } else {
    api_key[[1]]
  }
  key <- pa_safe_trim(key)[[1]]
  has_live_key <- nchar(key) >= 10L && !pa_ge_is_placeholder_secret(key)
  if (!has_live_key && isTRUE(write_env_placeholder)) pa_ge_ensure_env_placeholder("/home/h2048/.env", "DEEPSEEK_API_KEY")
  list(
    api_key = if (has_live_key) key else "",
    has_live_key = has_live_key,
    env_files_loaded = env_candidates,
    status = if (has_live_key) "ready" else "missing_or_placeholder_key"
  )
}

pa_ge_markdown_table <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) return("No rows available.")
  df[] <- lapply(df, function(col) gsub("\\|", "/", pa_safe_trim(col)))
  header <- paste(c("", colnames(df), ""), collapse = "|")
  sep <- paste(c("", rep("---", ncol(df)), ""), collapse = "|")
  rows <- apply(df, 1, function(row) paste(c("", row, ""), collapse = "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

pa_build_gene_exclusion_llm_prompt <- function(packet,
                                               extra_context = NULL,
                                               max_examples = 10L) {
  if (is.null(packet) || !is.list(packet) || is.null(packet$summary_table) || is.null(packet$audit_table)) {
    stop("packet must be a gene exclusion packet", call. = FALSE)
  }
  max_examples <- suppressWarnings(as.integer(max_examples))
  if (is.na(max_examples) || max_examples < 1L) max_examples <- 10L
  manifest <- packet$manifest
  summary_table <- as.data.frame(packet$summary_table, stringsAsFactors = FALSE, check.names = FALSE)
  audit <- as.data.frame(packet$audit_table, stringsAsFactors = FALSE, check.names = FALSE)
  excluded_examples <- do.call(rbind, lapply(unique(summary_table$category), function(cat_name) {
    idx <- audit$primary_category == cat_name & audit$should_exclude
    genes <- paste(head(audit$gene_symbol[idx], max_examples), collapse = ", ")
    data.frame(category = cat_name, excluded_examples = genes, stringsAsFactors = FALSE)
  }))

  system_prompt <- paste(
    "You are a single-cell RNA-seq gene-program QC expert.",
    "Interpret gene-exclusion audit results conservatively.",
    "Respond in concise Chinese Markdown.",
    "Do not treat mitochondrial, ribosomal, lncRNA, or immunoglobulin genes as biological program drivers without lineage context."
  )
  user_prompt <- paste(
    "请解读下面的单细胞基因过滤/排除审计结果，并给出后续 program discovery / trajectory / enrichment 分析建议。",
    "",
    "## Manifest",
    sprintf("- spec_version: %s", pa_null_coalesce(manifest$spec_version, NA_character_)),
    sprintf("- lineage_context: %s", pa_null_coalesce(manifest$lineage_context, NA_character_)),
    sprintf("- ig_policy: %s", pa_null_coalesce(manifest$ig_policy, NA_character_)),
    sprintf("- keep_ig: %s", isTRUE(manifest$keep_ig)),
    sprintf("- annotation_used: %s", isTRUE(manifest$annotation_used)),
    sprintf("- n_total_features: %s", pa_null_coalesce(manifest$n_total_features, NA_integer_)),
    sprintf("- n_excluded_features: %s", pa_null_coalesce(manifest$n_excluded_features, NA_integer_)),
    sprintf("- n_kept_features: %s", pa_null_coalesce(manifest$n_kept_features, NA_integer_)),
    "",
    "## Category summary",
    pa_ge_markdown_table(summary_table[, c("category", "detected_n", "excluded_n", "kept_n", "pct_of_all_features", "pct_of_all_excluded", "example_genes"), drop = FALSE]),
    "",
    "## Excluded examples by category",
    pa_ge_markdown_table(excluded_examples),
    "",
    if (!is.null(extra_context) && length(extra_context) > 0L) paste("## Extra context\n", paste(pa_safe_trim(extra_context), collapse = "\n"), sep = "") else "",
    "",
    "请输出：",
    "1. 过滤是否符合该 lineage context；",
    "2. 可能误删/应保留的风险点，特别是 IG 与 lncRNA；",
    "3. 对下游 cNMF/hdWGCNA/CoVarNet/ssGSEA/trajectory 的影响；",
    "4. 建议的人工复核清单。",
    sep = "\n"
  )
  list(system_prompt = system_prompt, user_prompt = user_prompt)
}

pa_build_gene_exclusion_rule_based_interpretation <- function(packet,
                                                              status = "skipped") {
  manifest <- packet$manifest
  summary_table <- as.data.frame(packet$summary_table, stringsAsFactors = FALSE, check.names = FALSE)
  total <- suppressWarnings(as.numeric(pa_null_coalesce(manifest$n_total_features, NA_real_)))
  excluded <- suppressWarnings(as.numeric(pa_null_coalesce(manifest$n_excluded_features, NA_real_)))
  excluded_pct <- if (is.finite(total) && total > 0) 100 * excluded / total else NA_real_
  dominant <- summary_table[order(-summary_table$excluded_n), , drop = FALSE]
  dominant <- dominant[dominant$excluded_n > 0, , drop = FALSE]
  dominant_txt <- if (nrow(dominant) > 0L) {
    paste(sprintf("- `%s`: %s excluded (%s)", dominant$category, dominant$excluded_n, dominant$example_genes), collapse = "\n")
  } else {
    "- No excluded categories."
  }
  c(
    "# Gene exclusion interpretation",
    "",
    sprintf("**LLM status:** `%s`", status),
    "",
    "## Rule-based summary",
    sprintf("- Total features: `%s`", pa_null_coalesce(manifest$n_total_features, NA_integer_)),
    sprintf("- Excluded features: `%s`%s", pa_null_coalesce(manifest$n_excluded_features, NA_integer_), if (is.na(excluded_pct)) "" else sprintf(" (`%.2f%%`)", excluded_pct)),
    sprintf("- Kept features: `%s`", pa_null_coalesce(manifest$n_kept_features, NA_integer_)),
    sprintf("- Lineage context: `%s`", pa_null_coalesce(manifest$lineage_context, "NA")),
    sprintf("- IG policy: `%s`; keep IG: `%s`", pa_null_coalesce(manifest$ig_policy, "NA"), isTRUE(manifest$keep_ig)),
    sprintf("- Annotation used: `%s`", isTRUE(manifest$annotation_used)),
    "",
    "## Main excluded categories",
    dominant_txt,
    "",
    "## QC guidance",
    "- For non-B/plasma contexts, excluding immunoglobulin genes is usually appropriate to avoid ambient B-cell/plasma-cell signal dominating program discovery.",
    "- For B/plasma contexts, IG genes are retained by the auto policy and should be interpreted as lineage biology rather than technical contamination.",
    "- Mitochondrial and ribosomal genes should generally not anchor cNMF/hdWGCNA/CoVarNet programs unless a stress/translation-focused analysis is explicitly intended.",
    "- lncRNA exclusions should be reviewed if the project has a specific non-coding RNA hypothesis or high-quality gene-biotype annotations."
  )
}

pa_run_gene_exclusion_llm_interpretation <- function(packet,
                                                     output_dir,
                                                     prefix = "gene_exclusion",
                                                     llm_config = list()) {
  output_dir <- pa_prepare_output_dir(output_dir)
  prefix <- pa_scalar_chr(prefix, "prefix")
  if (is.null(llm_config) || !is.list(llm_config)) llm_config <- list()
  enabled <- pa_as_logical(pa_null_coalesce(llm_config$enabled, TRUE), default = TRUE)
  model <- pa_safe_trim(pa_null_coalesce(llm_config$model, "deepseek-reasoner"))[[1]]
  if (!nzchar(model)) model <- "deepseek-reasoner"
  env_candidates <- pa_null_coalesce(llm_config$env_candidates, c("/home/h2048/.env", "/home/h2048/script/.env"))
  write_env_placeholder <- pa_as_logical(pa_null_coalesce(llm_config$write_env_placeholder, FALSE), default = FALSE)
  timeout_sec <- suppressWarnings(as.numeric(pa_null_coalesce(llm_config$timeout_sec, 180)))
  if (is.na(timeout_sec) || timeout_sec <= 0) timeout_sec <- 180

  prompt <- pa_build_gene_exclusion_llm_prompt(
    packet = packet,
    extra_context = llm_config$extra_context,
    max_examples = pa_null_coalesce(llm_config$max_examples, 10L)
  )
  prompt_path <- file.path(output_dir, sprintf("%s_LLM_prompt.md", prefix))
  interpretation_path <- file.path(output_dir, sprintf("%s_LLM_interpretation.md", prefix))
  raw_path <- file.path(output_dir, sprintf("%s_LLM_raw_response.txt", prefix))
  status_path <- file.path(output_dir, sprintf("%s_LLM_status.json", prefix))
  pa_write_markdown(c("# Gene exclusion LLM prompt", "", "## System", prompt$system_prompt, "", "## User", prompt$user_prompt), prompt_path)

  status <- list(
    enabled = enabled,
    status = "skipped_disabled",
    model = model,
    prompt_md = prompt_path,
    interpretation_md = interpretation_path,
    raw_response_txt = raw_path,
    error = NULL
  )

  if (!enabled) {
    lines <- pa_build_gene_exclusion_rule_based_interpretation(packet, status = "skipped_disabled")
    pa_write_markdown(lines, interpretation_path)
    pa_write_json(status, status_path)
    status$status_json <- status_path
    return(status)
  }

  key_info <- pa_ge_resolve_deepseek_key(
    api_key = llm_config$api_key,
    env_candidates = env_candidates,
    write_env_placeholder = write_env_placeholder
  )
  status$env_files_loaded <- key_info$env_files_loaded
  if (!isTRUE(key_info$has_live_key)) {
    status$status <- "skipped_no_live_key"
    lines <- pa_build_gene_exclusion_rule_based_interpretation(packet, status = status$status)
    lines <- c(lines, "", "## How to enable live LLM", "Set a real `DEEPSEEK_API_KEY` in `/home/h2048/.env` or pass `llm_config = list(api_key = ...)`.")
    pa_write_markdown(lines, interpretation_path)
    pa_write_json(status, status_path)
    status$status_json <- status_path
    return(status)
  }

  response <- tryCatch({
    if (!exists("tc_deepseek_chat_request", mode = "function")) {
      source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R")
    }
    tc_deepseek_chat_request(
      prompt = prompt$user_prompt,
      model = model,
      api_key = key_info$api_key,
      base_url = llm_config$base_url,
      system_prompt = prompt$system_prompt,
      thinking = llm_config$thinking,
      reasoning_effort = llm_config$reasoning_effort,
      timeout_sec = timeout_sec
    )
  }, error = function(e) {
    structure(conditionMessage(e), class = "pa_gene_exclusion_llm_error")
  })

  if (inherits(response, "pa_gene_exclusion_llm_error")) {
    status$status <- "error"
    status$error <- as.character(response[[1]])
    lines <- c(
      pa_build_gene_exclusion_rule_based_interpretation(packet, status = "error"),
      "",
      "## LLM error",
      sprintf("`%s`", status$error)
    )
    pa_write_markdown(lines, interpretation_path)
  } else {
    status$status <- "ok"
    pa_write_markdown(c("# Gene exclusion LLM interpretation", "", response), interpretation_path)
    writeLines(response, raw_path, useBytes = TRUE)
  }
  pa_write_json(status, status_path)
  status$status_json <- status_path
  status
}

pa_export_gene_exclusion_packet <- function(packet,
                                            output_dir,
                                            prefix = "gene_exclusion",
                                            llm_config = list()) {
  output_dir <- pa_prepare_output_dir(output_dir)
  prefix <- pa_scalar_chr(prefix, "prefix")
  summary_path <- file.path(output_dir, sprintf("%s_summary.csv", prefix))
  audit_path <- file.path(output_dir, sprintf("%s_audit.csv", prefix))
  manifest_path <- file.path(output_dir, sprintf("%s_manifest.json", prefix))

  utils::write.csv(packet$summary_table, summary_path, row.names = FALSE, quote = TRUE)
  utils::write.csv(packet$audit_table, audit_path, row.names = FALSE, quote = TRUE)
  visualization_manifest <- pa_plot_gene_exclusion_visualizations(packet, output_dir = output_dir, prefix = prefix)
  llm_manifest <- pa_run_gene_exclusion_llm_interpretation(packet, output_dir = output_dir, prefix = prefix, llm_config = llm_config)

  manifest <- utils::modifyList(packet$manifest, list(
    output_dir = output_dir,
    summary_csv = summary_path,
    audit_csv = audit_path,
    visualizations = visualization_manifest,
    llm = llm_manifest,
    helper_version = PA_GENE_EXCLUSION_HELPER_VERSION_20260505_V1
  ))
  pa_write_json(manifest, manifest_path)
  manifest$manifest_json <- manifest_path
  manifest
}

pa_apply_gene_exclusion_to_seurat <- function(seurat_obj,
                                              output_dir,
                                              gene_exclusion_config = list(),
                                              prefix = "gene_exclusion") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object", call. = FALSE)
  }
  if (is.null(gene_exclusion_config) || !is.list(gene_exclusion_config)) gene_exclusion_config <- list()
  spec <- pa_load_gene_exclusion_spec(overrides = gene_exclusion_config)
  lineage_context <- pa_null_coalesce(gene_exclusion_config$lineage_context, pa_infer_lineage_context_from_seurat(seurat_obj))
  packet <- pa_build_gene_exclusion_packet(
    feature_names = rownames(seurat_obj),
    spec = spec,
    lineage_context = lineage_context,
    annotation_df = gene_exclusion_config$annotation_df,
    annotation_path = gene_exclusion_config$annotation_path
  )
  manifest <- pa_export_gene_exclusion_packet(
    packet,
    output_dir = output_dir,
    prefix = prefix,
    llm_config = pa_null_coalesce(gene_exclusion_config$llm_config, gene_exclusion_config$llm)
  )

  keep_genes <- intersect(packet$keep_genes, rownames(seurat_obj))
  filtered_obj <- seurat_obj[keep_genes, , drop = FALSE]
  filtered_obj@misc$program_gene_exclusion <- list(
    manifest = manifest,
    summary_table = packet$summary_table
  )
  list(seurat_obj = filtered_obj, packet = packet, manifest = manifest)
}

if (sys.nframe() == 0) {
  cat("Program Gene Exclusion Helper (2026-05-05 v1) loaded.\n")
}
