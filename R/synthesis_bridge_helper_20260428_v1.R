#!/usr/bin/env Rscript
# ==============================================================================
# Synthesis Bridge Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)

PA_TISSUE_COMPARISON_ADVANCED_HELPER_PATH_20260408 <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R"

pa_resolve_orgdb <- function(species = "human") {
  species_key <- tolower(pa_scalar_chr(species, "species"))
  pkg <- switch(
    species_key,
    human = "org.Hs.eg.db",
    homo_sapiens = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    mus_musculus = "org.Mm.eg.db",
    stop(sprintf("Unsupported species for enrichment runner: %s", species), call. = FALSE)
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required OrgDb package is not installed: %s", pkg), call. = FALSE)
  }
  list(
    package = pkg,
    object = get(pkg, envir = asNamespace(pkg)),
    msigdbr_species = if (identical(pkg, "org.Hs.eg.db")) "Homo sapiens" else "Mus musculus"
  )
}

pa_prepare_gene_stats <- function(gene_stats) {
  if (is.null(gene_stats)) return(NULL)
  if (is.numeric(gene_stats) && !is.null(names(gene_stats))) {
    out <- gene_stats
    out <- out[is.finite(out)]
    return(sort(out, decreasing = TRUE))
  }
  if (is.data.frame(gene_stats)) {
    gene_col <- c("gene", "symbol", "feature")[c("gene", "symbol", "feature") %in% colnames(gene_stats)][1]
    stat_col <- c("log2FC", "avg_log2FC", "avg_logFC", "stat")[c("log2FC", "avg_log2FC", "avg_logFC", "stat") %in% colnames(gene_stats)][1]
    if (is.na(gene_col) || is.na(stat_col)) return(NULL)
    vals <- suppressWarnings(as.numeric(gene_stats[[stat_col]]))
    names(vals) <- as.character(gene_stats[[gene_col]])
    vals <- vals[is.finite(vals) & nzchar(names(vals))]
    return(sort(vals, decreasing = TRUE))
  }
  NULL
}

pa_top_enrichment_terms_text <- function(enrichment_packet, n = 12L) {
  ora_tbl <- if (is.list(enrichment_packet)) enrichment_packet$ora_table else NULL
  if (is.null(ora_tbl) || !is.data.frame(ora_tbl) || nrow(ora_tbl) == 0L) return("No enrichment terms available.")
  desc_col <- c("Description", "term", "ID")[c("Description", "term", "ID") %in% colnames(ora_tbl)][1]
  db_col <- if ("db" %in% colnames(ora_tbl)) "db" else NA_character_
  padj_col <- c("p.adjust", "qvalue", "pvalue")[c("p.adjust", "qvalue", "pvalue") %in% colnames(ora_tbl)][1]
  rows <- head(ora_tbl, max(1L, as.integer(n)))
  vapply(seq_len(nrow(rows)), function(i) {
    bits <- c()
    if (!is.na(db_col)) bits <- c(bits, sprintf("[%s]", pa_safe_trim(rows[[db_col]][i])))
    bits <- c(bits, pa_safe_trim(rows[[desc_col]][i]))
    if (!is.na(padj_col)) bits <- c(bits, sprintf("padj=%s", signif(as.numeric(rows[[padj_col]][i]), 3)))
    paste(bits, collapse = " ")
  }, character(1)) |> paste(collapse = "\n")
}

pa_top_gene_fc_text <- function(gene_fc, n = 12L) {
  if (is.null(gene_fc)) return("No gene-level fold-change evidence provided.")
  stats <- pa_prepare_gene_stats(gene_fc)
  if (is.null(stats) || length(stats) == 0L) return("No gene-level fold-change evidence provided.")
  top_stats <- head(stats, max(1L, as.integer(n)))
  paste(sprintf("%s(%.2f)", names(top_stats), as.numeric(top_stats)), collapse = ", ")
}

pa_strip_json_fence <- function(text) {
  text <- paste(as.character(text), collapse = "\n")
  text <- trimws(text)
  text <- gsub("^```[A-Za-z0-9_+-]*\\s*", "", text, perl = TRUE)
  text <- gsub("\\s*```$", "", text, perl = TRUE)
  trimws(text)
}

pa_parse_interpret_json <- function(text) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  clean <- pa_strip_json_fence(text)
  tryCatch(jsonlite::fromJSON(clean, simplifyVector = TRUE), error = function(e) NULL)
}

pa_interpret_payload_to_table <- function(payload, raw_text, status) {
  if (is.null(payload) || !is.list(payload)) {
    return(data.frame(
      overview = if (identical(status, "raw_text")) pa_safe_trim(raw_text) else NA_character_,
      key_mechanisms = NA_character_,
      hypothesis = NA_character_,
      narrative = if (identical(status, "raw_text")) pa_safe_trim(raw_text) else NA_character_,
      key_drivers = NA_character_,
      evidence = NA_character_,
      limitations = NA_character_,
      status = status,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  key_drivers <- payload$key_drivers
  if (is.list(key_drivers)) key_drivers <- unlist(key_drivers, use.names = FALSE)
  data.frame(
    overview = pa_null_coalesce(pa_safe_trim(payload$overview), NA_character_),
    key_mechanisms = pa_null_coalesce(pa_safe_trim(payload$key_mechanisms), NA_character_),
    hypothesis = pa_null_coalesce(pa_safe_trim(payload$hypothesis), NA_character_),
    narrative = pa_null_coalesce(pa_safe_trim(payload$narrative), NA_character_),
    key_drivers = if (length(key_drivers) == 0L) NA_character_ else paste(as.character(key_drivers), collapse = ", "),
    evidence = pa_null_coalesce(pa_safe_trim(payload$evidence), NA_character_),
    limitations = pa_null_coalesce(pa_safe_trim(payload$limitations), NA_character_),
    status = status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

pa_build_interpret_agent_prompt <- function(enrichment_packet,
                                            context_str = "",
                                            gene_fc = NULL) {
  paste(
    "You are a computational biologist summarizing pathway enrichment results from a single-cell analysis.",
    "Return strict JSON with keys: overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "Write concise Simplified Chinese.",
    "Do not invent pathways or genes that are not present in the evidence.",
    "",
    sprintf("Context: %s", pa_null_coalesce(context_str, "")),
    "",
    "Top enrichment terms:",
    pa_top_enrichment_terms_text(enrichment_packet),
    "",
    "Top genes / fold changes:",
    pa_top_gene_fc_text(gene_fc),
    sep = "\n"
  )
}

pa_source_deepseek_helper <- function(helper_path = PA_TISSUE_COMPARISON_ADVANCED_HELPER_PATH_20260408) {
  helper_path <- pa_scalar_chr(helper_path, "helper_path")
  if (!file.exists(helper_path)) stop(sprintf("DeepSeek helper file does not exist: %s", helper_path), call. = FALSE)
  source(helper_path, local = .GlobalEnv)
  if (exists("tc_load_env_file", mode = "function", inherits = TRUE)) {
    tc_load_env_file("/home/h2048/.env")
  }
  invisible(helper_path)
}

pa_build_enrichment_packet <- function(
  ora_table = NULL,
  gsea_table = NULL,
  source_dbs = NULL,
  source_paths = NULL,
  summary_note = NULL
) {
  list(
    packet_type = "enrichment",
    ora_table = ora_table,
    gsea_table = gsea_table,
    source_dbs = pa_null_coalesce(pa_unique_chr(source_dbs), character()),
    source_paths = pa_null_coalesce(pa_unique_chr(source_paths), character()),
    summary_note = pa_null_coalesce(summary_note, NA_character_),
    support_status = if (!is.null(ora_table) || !is.null(gsea_table)) "available" else "empty",
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_build_interpret_agent_packet <- function(
  structured_table = NULL,
  raw_payload = NULL,
  model = NA_character_,
  source_paths = NULL,
  status = NULL
) {
  list(
    packet_type = "interpret_agent",
    structured_table = structured_table,
    raw_payload = raw_payload,
    model = pa_null_coalesce(model, NA_character_),
    source_paths = pa_null_coalesce(pa_unique_chr(source_paths), character()),
    support_status = pa_null_coalesce(status, if (is.data.frame(structured_table) && nrow(structured_table) > 0L) "available" else "empty"),
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_build_report_packet <- function(
  report_md_path = NA_character_,
  figure_paths = NULL,
  table_paths = NULL,
  sections = NULL,
  report_title = NA_character_
) {
  report_exists <- !is.na(report_md_path) && nzchar(report_md_path) && file.exists(report_md_path)
  list(
    packet_type = "report",
    report_md_path = pa_null_coalesce(report_md_path, NA_character_),
    figure_paths = pa_null_coalesce(pa_unique_chr(figure_paths), character()),
    table_paths = pa_null_coalesce(pa_unique_chr(table_paths), character()),
    sections = pa_null_coalesce(pa_unique_chr(sections), character()),
    report_title = pa_null_coalesce(report_title, NA_character_),
    support_status = if (report_exists) "available" else "missing",
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_build_synthesis_packet <- function(
  enrichment_packet = NULL,
  interpret_agent_packet = NULL,
  report_packet = NULL
) {
  support_status <- if (!is.null(enrichment_packet) || !is.null(interpret_agent_packet) || !is.null(report_packet)) "available" else "empty"
  list(
    enrichment_packet = enrichment_packet,
    interpret_agent_packet = interpret_agent_packet,
    report_packet = report_packet,
    support_status = support_status,
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_synthesis_packet_summary_lines <- function(synthesis_packet) {
  if (is.null(synthesis_packet) || !is.list(synthesis_packet)) return("- 未提供 synthesis packet。")
  enrich_status <- if (is.list(synthesis_packet$enrichment_packet)) synthesis_packet$enrichment_packet$support_status else "无"
  interpret_status <- if (is.list(synthesis_packet$interpret_agent_packet)) synthesis_packet$interpret_agent_packet$support_status else "无"
  report_status <- if (is.list(synthesis_packet$report_packet)) synthesis_packet$report_packet$support_status else "无"
  c(
    sprintf("- enrichment：%s", pa_null_coalesce(enrich_status, "无")),
    sprintf("- interpret_agent：%s", pa_null_coalesce(interpret_status, "无")),
    sprintf("- REPORT：%s", pa_null_coalesce(report_status, "无")),
    sprintf("- synthesis status：%s", pa_null_coalesce(synthesis_packet$support_status, "无"))
  )
}

pa_run_enrichment_runner <- function(gene_vector,
                                     output_dir,
                                     species = "human",
                                     sources = c("GO_BP", "Hallmark"),
                                     universe = NULL,
                                     gene_stats = NULL,
                                     max_terms = 50L) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Package 'clusterProfiler' is required for enrichment runner.", call. = FALSE)
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package 'msigdbr' is required for enrichment runner.", call. = FALSE)

  output_dir <- pa_prepare_output_dir(output_dir)
  genes <- unique(pa_safe_trim(gene_vector))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0L) stop("gene_vector must contain at least one gene symbol", call. = FALSE)

  org_info <- pa_resolve_orgdb(species)
  sources <- unique(pa_safe_trim(sources))
  stats_vec <- pa_prepare_gene_stats(gene_stats)
  ora_rows <- list()
  gsea_rows <- list()
  idx <- 1L

  for (src in sources) {
    if (src %in% c("GO_BP", "GO_MF", "GO_CC")) {
      ont <- sub("GO_", "", src)
      res <- tryCatch(
        clusterProfiler::enrichGO(
          gene = genes,
          OrgDb = org_info$object,
          keyType = "SYMBOL",
          universe = universe,
          ont = ont,
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        df <- as.data.frame(res)
        if (nrow(df) > 0L) {
          df$db <- src
          ora_rows[[idx]] <- head(df, max(1L, as.integer(max_terms)))
          idx <- idx + 1L
        }
      }
      next
    }

    if (identical(src, "Hallmark")) {
      hallmark_df <- msigdbr::msigdbr(species = org_info$msigdbr_species, category = "H")
      term2gene <- unique(hallmark_df[, c("gs_name", "gene_symbol")])
      term2name <- unique(hallmark_df[, c("gs_name", "gs_name")])
      colnames(term2name) <- c("gs_name", "gs_description")
      ora_res <- tryCatch(
        clusterProfiler::enricher(gene = genes, universe = universe, TERM2GENE = term2gene, TERM2NAME = term2name),
        error = function(e) NULL
      )
      if (!is.null(ora_res)) {
        df <- as.data.frame(ora_res)
        if (nrow(df) > 0L) {
          df$db <- src
          ora_rows[[idx]] <- head(df, max(1L, as.integer(max_terms)))
          idx <- idx + 1L
        }
      }
      if (!is.null(stats_vec) && length(stats_vec) >= 10L) {
        gsea_res <- tryCatch(
          clusterProfiler::GSEA(geneList = stats_vec, TERM2GENE = term2gene, TERM2NAME = term2name, verbose = FALSE),
          error = function(e) NULL
        )
        if (!is.null(gsea_res)) {
          df <- as.data.frame(gsea_res)
          if (nrow(df) > 0L) {
            df$db <- src
            gsea_rows[[length(gsea_rows) + 1L]] <- head(df, max(1L, as.integer(max_terms)))
          }
        }
      }
    }
  }

  ora_table <- if (length(ora_rows) > 0L) do.call(rbind, ora_rows) else data.frame()
  gsea_table <- if (length(gsea_rows) > 0L) do.call(rbind, gsea_rows) else data.frame()
  rownames(ora_table) <- NULL
  rownames(gsea_table) <- NULL

  ora_path <- file.path(output_dir, "enrichment_ora.tsv")
  gsea_path <- file.path(output_dir, "enrichment_gsea.tsv")
  if (nrow(ora_table) > 0L) pa_write_tsv(ora_table, ora_path)
  if (nrow(gsea_table) > 0L) pa_write_tsv(gsea_table, gsea_path)

  enrichment_packet <- pa_build_enrichment_packet(
    ora_table = if (nrow(ora_table) > 0L) ora_table else NULL,
    gsea_table = if (nrow(gsea_table) > 0L) gsea_table else NULL,
    source_dbs = sources,
    source_paths = c(if (nrow(ora_table) > 0L) ora_path else NULL, if (nrow(gsea_table) > 0L) gsea_path else NULL),
    summary_note = sprintf("genes=%d | sources=%s", length(genes), paste(sources, collapse = ", "))
  )

  list(
    status = "ok",
    enrichment_packet = enrichment_packet,
    manifest = list(output_dir = output_dir, ora_path = if (nrow(ora_table) > 0L) ora_path else NA_character_, gsea_path = if (nrow(gsea_table) > 0L) gsea_path else NA_character_)
  )
}

pa_run_interpret_agent_runner <- function(enrichment_packet,
                                          output_dir,
                                          context_str = "",
                                          gene_fc = NULL,
                                          model = NULL,
                                          api_key = NULL,
                                          dry_run = FALSE) {
  output_dir <- pa_prepare_output_dir(output_dir)
  prompt <- pa_build_interpret_agent_prompt(enrichment_packet = enrichment_packet, context_str = context_str, gene_fc = gene_fc)
  prompt_path <- file.path(output_dir, "interpret_agent_prompt.txt")
  pa_write_markdown(prompt, prompt_path)

  if (isTRUE(dry_run)) {
    packet <- pa_build_interpret_agent_packet(
      structured_table = data.frame(
        overview = "Dry-run only; no LLM request executed.",
        key_mechanisms = NA_character_,
        hypothesis = NA_character_,
        narrative = "Prompt generated successfully.",
        key_drivers = NA_character_,
        evidence = NA_character_,
        limitations = "Dry-run mode.",
        status = "dry_run",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      raw_payload = list(prompt = prompt),
      model = pa_null_coalesce(model, "deepseek-v4-flash"),
      source_paths = prompt_path,
      status = "dry_run"
    )
    return(list(status = "dry_run", prompt_path = prompt_path, interpret_packet = packet))
  }

  pa_source_deepseek_helper()
  resolved_model <- if (!is.null(model) && nzchar(pa_safe_trim(model))) pa_safe_trim(model) else tc_deepseek_default_chat_model()
  raw_text <- tryCatch(
    tc_deepseek_chat_request(prompt, model = resolved_model, api_key = api_key),
    error = function(e) structure(conditionMessage(e), class = "pa_interpret_error")
  )

  raw_path <- file.path(output_dir, "interpret_agent_raw.txt")
  pa_write_markdown(as.character(raw_text), raw_path)

  parsed <- if (!inherits(raw_text, "pa_interpret_error")) pa_parse_interpret_json(raw_text) else NULL
  status <- if (inherits(raw_text, "pa_interpret_error")) "error" else if (!is.null(parsed)) "available" else "raw_text"
  structured_table <- pa_interpret_payload_to_table(parsed, raw_text = as.character(raw_text), status = status)
  structured_path <- file.path(output_dir, "interpret_agent_structured.tsv")
  pa_write_tsv(structured_table, structured_path)

  packet <- pa_build_interpret_agent_packet(
    structured_table = structured_table,
    raw_payload = list(prompt = prompt, response = as.character(raw_text)),
    model = resolved_model,
    source_paths = c(prompt_path, raw_path, structured_path),
    status = status
  )

  list(
    status = status,
    prompt_path = prompt_path,
    raw_path = raw_path,
    structured_path = structured_path,
    interpret_packet = packet
  )
}

pa_run_report_runner <- function(unit,
                                 output_dir,
                                 program_registry = NULL,
                                 trajectory_branch_packet = NULL,
                                 synthesis_packet = NULL,
                                 report_title = NULL) {
  output_dir <- pa_prepare_output_dir(output_dir)
  if (!is.data.frame(unit) || nrow(unit) == 0L) stop("unit must be a non-empty analysis unit data.frame", call. = FALSE)

  if (is.null(report_title) || !nzchar(pa_safe_trim(report_title))) {
    report_title <- sprintf("Program Architecture Report — %s", unit$unit_id[[1]])
  }

  lines <- c(
    sprintf("# %s", report_title),
    "",
    sprintf("- **unit_id**: `%s`", unit$unit_id[[1]]),
    sprintf("- **lineage**: `%s`", unit$lineage[[1]]),
    sprintf("- **state_level**: `%s`", unit$state_level[[1]]),
    sprintf("- **contrast**: `%s` (`%s` vs `%s`)", unit$contrast_id[[1]], unit$condition_a[[1]], unit$condition_b[[1]]),
    "",
    "## Program sources",
    ""
  )

  if (is.data.frame(program_registry) && nrow(program_registry) > 0L && "source_type" %in% colnames(program_registry)) {
    src_counts <- sort(table(program_registry$source_type), decreasing = TRUE)
    lines <- c(lines, sprintf("- %s: %d", names(src_counts), as.integer(src_counts)), "")
  } else {
    lines <- c(lines, "- No program registry provided.", "")
  }

  lines <- c(lines, "## Trajectory", "", pa_trajectory_packet_summary_lines(trajectory_branch_packet), "")
  lines <- c(lines, "## Synthesis", "", pa_synthesis_packet_summary_lines(synthesis_packet), "")

  ora_tbl <- NULL
  interpret_tbl <- NULL
  if (is.list(synthesis_packet) && is.list(synthesis_packet$enrichment_packet)) ora_tbl <- synthesis_packet$enrichment_packet$ora_table
  if (is.list(synthesis_packet) && is.list(synthesis_packet$interpret_agent_packet)) interpret_tbl <- synthesis_packet$interpret_agent_packet$structured_table

  lines <- c(lines, "## Top enrichment terms", "")
  if (is.data.frame(ora_tbl) && nrow(ora_tbl) > 0L) {
    desc_col <- c("Description", "term", "ID")[c("Description", "term", "ID") %in% colnames(ora_tbl)][1]
    db_col <- if ("db" %in% colnames(ora_tbl)) "db" else NA_character_
    top_tbl <- head(ora_tbl, 10)
    lines <- c(lines, vapply(seq_len(nrow(top_tbl)), function(i) {
      prefix <- if (!is.na(db_col)) sprintf("[%s] ", pa_safe_trim(top_tbl[[db_col]][i])) else ""
      sprintf("- %s%s", prefix, pa_safe_trim(top_tbl[[desc_col]][i]))
    }, character(1)), "")
  } else {
    lines <- c(lines, "- No enrichment rows available.", "")
  }

  lines <- c(lines, "## Interpret agent", "")
  if (is.data.frame(interpret_tbl) && nrow(interpret_tbl) > 0L) {
    rec <- interpret_tbl[1, , drop = FALSE]
    for (field in c("overview", "key_mechanisms", "hypothesis", "narrative", "key_drivers", "evidence", "limitations")) {
      if (field %in% colnames(rec)) {
        lines <- c(lines, sprintf("### %s", gsub("_", " ", field)), "", pa_null_coalesce(pa_safe_trim(rec[[field]][1]), "NA"), "")
      }
    }
  } else {
    lines <- c(lines, "- No interpret_agent output available.", "")
  }

  report_path <- file.path(output_dir, "REPORT.md")
  pa_write_markdown(lines, report_path)
  report_packet <- pa_build_report_packet(
    report_md_path = report_path,
    sections = c("program_sources", "trajectory", "synthesis", "top_enrichment_terms", "interpret_agent"),
    report_title = report_title
  )

  list(
    status = "ok",
    report_packet = report_packet,
    report_md_path = report_path
  )
}

pa_run_synthesis_runner <- function(unit,
                                    gene_vector,
                                    output_dir,
                                    context_str = "",
                                    gene_fc = NULL,
                                    enrichment_sources = c("GO_BP", "Hallmark"),
                                    interpret_dry_run = TRUE) {
  output_dir <- pa_prepare_output_dir(output_dir)
  enrich_res <- pa_run_enrichment_runner(
    gene_vector = gene_vector,
    output_dir = file.path(output_dir, "enrichment"),
    sources = enrichment_sources,
    gene_stats = gene_fc
  )
  interpret_res <- pa_run_interpret_agent_runner(
    enrichment_packet = enrich_res$enrichment_packet,
    output_dir = file.path(output_dir, "interpret_agent"),
    context_str = context_str,
    gene_fc = gene_fc,
    dry_run = interpret_dry_run
  )
  synthesis_packet <- pa_build_synthesis_packet(
    enrichment_packet = enrich_res$enrichment_packet,
    interpret_agent_packet = interpret_res$interpret_packet
  )
  report_res <- pa_run_report_runner(
    unit = unit,
    output_dir = file.path(output_dir, "report"),
    synthesis_packet = synthesis_packet
  )
  synthesis_packet$report_packet <- report_res$report_packet
  list(
    status = "ok",
    enrichment = enrich_res,
    interpret = interpret_res,
    report = report_res,
    synthesis_packet = synthesis_packet
  )
}

if (sys.nframe() == 0) {
  cat("Synthesis Bridge Helper (2026-04-28 v1) loaded.\n")
}
