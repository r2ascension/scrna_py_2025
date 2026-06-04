#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

DATE_TAG <- "20260525"
DEFAULT_OUTPUT_DIR <- file.path("/home/h2048/output", paste0("tissue_deg_llm_crosslineage_", DATE_TAG))
PADJ_THR <- 0.05
LFC_THR <- 1
TOP_N <- 10

LINEAGE_DIRS <- list(
  bcell = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
  epithelial = "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
  endothelial = "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
  tnk = "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508",
  myeloid = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
  fibroblast = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414",
  smc = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

safe_trim <- function(x) {
  x <- as.character(x %||% "")
  x[is.na(x)] <- ""
  trimws(x[[1]])
}

collapse_field <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- x[nzchar(x)]
  paste(unique(x), collapse = ", ")
}

pick_col <- function(dt, candidates) {
  hit <- intersect(candidates, colnames(dt))
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_fread <- function(path) {
  if (!file.exists(path)) return(data.table())
  tryCatch(fread(path, data.table = TRUE, showProgress = FALSE), error = function(e) data.table())
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

safe_read_text <- function(path, n = 60L) {
  if (!file.exists(path)) return(character())
  out <- tryCatch(readLines(path, warn = FALSE, n = n), error = function(e) character())
  Encoding(out) <- "UTF-8"
  out
}

fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x[[1]]) || !is.finite(x[[1]])) return("NA")
  format(round(x[[1]], digits), nsmall = digits, trim = TRUE)
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x[[1]]) || !is.finite(x[[1]])) return("NA")
  format(signif(x[[1]], 3), scientific = TRUE, trim = TRUE)
}

md_escape <- function(x) {
  x <- gsub("\\|", "\\\\|", as.character(x))
  x
}

comparison_display <- function(x) {
  x <- safe_trim(x)
  if (!nzchar(x)) return("")
  gsub("_", " ", x, fixed = TRUE)
}

summarize_deg <- function(dt, padj_thr = PADJ_THR, lfc_thr = LFC_THR, top_n = TOP_N) {
  if (!nrow(dt)) {
    return(list(
      n_genes = 0L, n_sig = 0L, n_up = 0L, n_down = 0L,
      top_up = "", top_down = "", top_abs = ""
    ))
  }
  gene_col <- pick_col(dt, c("gene", "symbol", "feature", "features"))
  lfc_col <- pick_col(dt, c("log2FoldChange", "avg_log2FC", "log2FC", "avg_logFC"))
  padj_col <- pick_col(dt, c("padj", "p_val_adj", "FDR", "qvalue"))
  if (is.na(gene_col) || is.na(lfc_col)) {
    return(list(
      n_genes = nrow(dt), n_sig = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
      top_up = "", top_down = "", top_abs = ""
    ))
  }
  dt <- copy(dt)
  dt[, gene_symbol := as.character(get(gene_col))]
  dt[, log2fc_num := suppressWarnings(as.numeric(get(lfc_col)))]
  dt[, padj_num := if (!is.na(padj_col)) suppressWarnings(as.numeric(get(padj_col))) else NA_real_]
  dt <- dt[!is.na(gene_symbol) & nzchar(trimws(gene_symbol)) & is.finite(log2fc_num)]
  sig_flag <- abs(dt$log2fc_num) >= lfc_thr
  if (!is.na(padj_col)) {
    sig_flag <- sig_flag & is.finite(dt$padj_num) & dt$padj_num < padj_thr
  }
  sig_dt <- dt[sig_flag]
  up_dt <- sig_dt[log2fc_num > 0][order(padj_num, -log2fc_num)]
  down_dt <- sig_dt[log2fc_num < 0][order(padj_num, log2fc_num)]
  abs_dt <- sig_dt[order(-abs(log2fc_num), padj_num)]
  fmt_gene <- function(subdt) {
    if (!nrow(subdt)) return("")
    subdt <- subdt[seq_len(min(nrow(subdt), top_n))]
    paste(sprintf("%s(log2FC=%s,padj=%s)", subdt$gene_symbol, fmt_num(subdt$log2fc_num), fmt_p(subdt$padj_num)), collapse = "; ")
  }
  list(
    n_genes = nrow(dt),
    n_sig = nrow(sig_dt),
    n_up = nrow(up_dt),
    n_down = nrow(down_dt),
    top_up = fmt_gene(up_dt),
    top_down = fmt_gene(down_dt),
    top_abs = fmt_gene(abs_dt)
  )
}

extract_llm_fields <- function(obj) {
  if (is.null(obj) || !is.list(obj)) {
    return(list(
      status = "missing", overview = "", key_mechanisms = "", hypothesis = "",
      narrative = "", key_drivers = "", evidence = "", limitations = "", error = ""
    ))
  }
  list(
    status = safe_trim(obj$status %||% "structured"),
    overview = safe_trim(obj$overview),
    key_mechanisms = safe_trim(obj$key_mechanisms),
    hypothesis = safe_trim(obj$hypothesis),
    narrative = safe_trim(obj$narrative),
    key_drivers = collapse_field(obj$key_drivers),
    evidence = safe_trim(obj$evidence),
    limitations = safe_trim(obj$limitations),
    error = safe_trim(obj$error)
  )
}

parse_contrast_record <- function(lineage, root_dir, level_label, csv_path) {
  rel <- sub(paste0("^", normalizePath(root_dir, winslash = "/", mustWork = FALSE), "/"), "", normalizePath(csv_path, winslash = "/", mustWork = FALSE))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(NULL)
  celltype <- parts[[2]]
  comparison <- parts[[3]]
  contrast_dir <- dirname(csv_path)
  llm_rds <- file.path(contrast_dir, "interpret_agent_integrated_up_down_structured.rds")
  llm_raw <- file.path(contrast_dir, "interpret_agent_integrated_up_down_raw.txt")
  evidence_csv <- file.path(contrast_dir, "llm_multi_db_evidence_integrated_up_down.csv")
  pathway_map_csv <- file.path(contrast_dir, "validated_gene_pathway_map_integrated_up_down.csv")

  de_dt <- safe_fread(csv_path)
  deg <- summarize_deg(de_dt)
  llm_obj <- safe_read_rds(llm_rds)
  llm <- extract_llm_fields(llm_obj)
  raw_preview <- paste(safe_read_text(llm_raw, n = 20L), collapse = "\n")
  evidence_dt <- safe_fread(evidence_csv)
  pathway_dt <- safe_fread(pathway_map_csv)

  data.table(
    lineage = lineage,
    root_dir = normalizePath(root_dir, winslash = "/", mustWork = FALSE),
    level = level_label,
    celltype = celltype,
    comparison = comparison,
    comparison_display = comparison_display(comparison),
    contrast_dir = normalizePath(contrast_dir, winslash = "/", mustWork = FALSE),
    deg_csv = normalizePath(csv_path, winslash = "/", mustWork = FALSE),
    llm_rds = normalizePath(llm_rds, winslash = "/", mustWork = FALSE),
    llm_raw = normalizePath(llm_raw, winslash = "/", mustWork = FALSE),
    evidence_csv = normalizePath(evidence_csv, winslash = "/", mustWork = FALSE),
    pathway_map_csv = normalizePath(pathway_map_csv, winslash = "/", mustWork = FALSE),
    n_genes = deg$n_genes,
    n_sig = deg$n_sig,
    n_up = deg$n_up,
    n_down = deg$n_down,
    top_abs = deg$top_abs,
    top_up = deg$top_up,
    top_down = deg$top_down,
    llm_status = llm$status,
    llm_overview = llm$overview,
    llm_key_mechanisms = llm$key_mechanisms,
    llm_hypothesis = llm$hypothesis,
    llm_narrative = llm$narrative,
    llm_key_drivers = llm$key_drivers,
    llm_evidence = llm$evidence,
    llm_limitations = llm$limitations,
    llm_error = llm$error,
    llm_raw_preview = raw_preview,
    n_evidence_rows = nrow(evidence_dt),
    n_pathway_map_rows = nrow(pathway_dt)
  )
}

collect_lineage_records <- function(lineage, root_dir) {
  if (!dir.exists(root_dir)) return(data.table())
  level_specs <- list(L2 = "pseudobulk_de", L3 = "pseudobulk_de_L3")
  rows <- list()
  for (level_label in names(level_specs)) {
    level_dir <- file.path(root_dir, level_specs[[level_label]])
    if (!dir.exists(level_dir)) next
    csv_files <- list.files(level_dir, pattern = "DESeq2_results\\.csv$", recursive = TRUE, full.names = TRUE)
    if (!length(csv_files)) next
    for (csv_path in csv_files) {
      rec <- parse_contrast_record(lineage, root_dir, level_label, csv_path)
      if (!is.null(rec)) rows[[length(rows) + 1L]] <- rec
    }
  }
  if (!length(rows)) return(data.table())
  rbindlist(rows, fill = TRUE)
}

build_overview_tables <- function(master_dt) {
  if (!nrow(master_dt)) {
    return(list(
      lineage_level = data.table(),
      lineage_celltype = data.table(),
      top_contrasts = data.table()
    ))
  }
  safe_median_num <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(stats::median(x, na.rm = TRUE))
  }
  safe_max_num <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(max(x, na.rm = TRUE))
  }
  safe_sum_num <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) return(0)
    as.numeric(sum(x, na.rm = TRUE))
  }
  safe_mean_num <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(mean(x, na.rm = TRUE))
  }
  lineage_level <- master_dt[, .(
    n_contrasts = as.integer(.N),
    n_celltypes = as.integer(uniqueN(celltype)),
    contrasts_with_llm = as.integer(sum(llm_status == "structured", na.rm = TRUE)),
    median_sig = safe_median_num(n_sig),
    max_sig = safe_max_num(n_sig)
  ), by = .(lineage, level)][order(lineage, level)]

  lineage_celltype <- master_dt[, .(
    n_contrasts = as.integer(.N),
    total_sig = safe_sum_num(n_sig),
    max_sig = safe_max_num(n_sig),
    mean_sig = round(safe_mean_num(n_sig), 1)
  ), by = .(lineage, level, celltype)][order(lineage, level, -total_sig)]

  top_contrasts <- master_dt[order(-n_sig, lineage, level, celltype)][, .(
    lineage, level, celltype, comparison, n_sig, n_up, n_down, llm_overview
  )]
  top_contrasts <- top_contrasts[seq_len(min(nrow(top_contrasts), 80L))]

  list(
    lineage_level = lineage_level,
    lineage_celltype = lineage_celltype,
    top_contrasts = top_contrasts
  )
}

write_markdown_report <- function(master_dt, overview, out_path) {
  lines <- c(
    "# Cross-lineage tissue DEG + LLM detailed summary",
    "",
    sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("**Threshold for DEG counts:** padj < %.2f and |log2FC| >= %d", PADJ_THR, LFC_THR),
    "",
    "## Scope",
    "",
    paste0("- Included lineages: ", paste(unique(master_dt$lineage), collapse = ", ")),
    "- Source evidence comes from existing tissue-comparison outputs (`pseudobulk_de/`, `pseudobulk_de_L3/`, and structured LLM `.rds` files).",
    "- This report does **not** rerun Seurat/DESeq2/LLM; it consolidates current finished outputs into a delivery-oriented summary.",
    "",
    "## Lineage-level overview",
    "",
    "| Lineage | Level | #Contrasts | #Celltypes | Structured LLM | Median sig DEG | Max sig DEG |",
    "|---|---:|---:|---:|---:|---:|---:|"
  )

  if (nrow(overview$lineage_level)) {
    for (i in seq_len(nrow(overview$lineage_level))) {
      row <- overview$lineage_level[i]
      lines <- c(lines, sprintf(
        "| %s | %s | %s | %s | %s | %s | %s |",
        md_escape(row$lineage), md_escape(row$level), row$n_contrasts, row$n_celltypes,
        row$contrasts_with_llm, fmt_num(row$median_sig, 1), row$max_sig
      ))
    }
  }

  lines <- c(lines, "", "## Highest-burden contrasts (top by significant DEG count)", "", "| Lineage | Level | Cell type | Comparison | Significant | Up | Down | LLM overview |", "|---|---|---|---|---:|---:|---:|---|")
  if (nrow(overview$top_contrasts)) {
    for (i in seq_len(nrow(overview$top_contrasts))) {
      row <- overview$top_contrasts[i]
      lines <- c(lines, sprintf(
        "| %s | %s | %s | %s | %s | %s | %s | %s |",
        md_escape(row$lineage), md_escape(row$level), md_escape(row$celltype), md_escape(comparison_display(row$comparison)),
        row$n_sig, row$n_up, row$n_down, md_escape(row$llm_overview)
      ))
    }
  }

  for (lineage_name in unique(master_dt$lineage)) {
    lineage_dt <- master_dt[lineage == lineage_name][order(level, celltype, comparison)]
    lines <- c(lines, "", sprintf("## %s", lineage_name), "")
    lineage_root <- unique(lineage_dt$root_dir)
    lines <- c(lines, sprintf("**Source dir:** `%s`", lineage_root[[1]]), "")

    celltype_summary <- overview$lineage_celltype[lineage == lineage_name]
    if (nrow(celltype_summary)) {
      lines <- c(lines, "### Cell-type burden summary", "", "| Level | Cell type | #Contrasts | Total sig DEG | Max sig DEG | Mean sig DEG |", "|---|---|---:|---:|---:|---:|")
      for (i in seq_len(nrow(celltype_summary))) {
        row <- celltype_summary[i]
        lines <- c(lines, sprintf(
          "| %s | %s | %s | %s | %s | %s |",
          md_escape(row$level), md_escape(row$celltype), row$n_contrasts, row$total_sig, row$max_sig, fmt_num(row$mean_sig, 1)
        ))
      }
      lines <- c(lines, "")
    }

    for (level_label in unique(lineage_dt$level)) {
      level_dt <- lineage_dt[level == level_label]
      lines <- c(lines, sprintf("### %s detailed contrasts", level_label), "")
      for (celltype_name in unique(level_dt$celltype)) {
        ct_dt <- level_dt[celltype == celltype_name][order(-n_sig, comparison)]
        lines <- c(lines, sprintf("#### %s", celltype_name), "")
        for (i in seq_len(nrow(ct_dt))) {
          row <- ct_dt[i]
          lines <- c(
            lines,
            sprintf("##### %s", comparison_display(row$comparison)),
            "",
            sprintf("- **Significant DEG:** %s (up=%s, down=%s)", row$n_sig, row$n_up, row$n_down),
            sprintf("- **Top absolute-effect genes:** %s", ifelse(nzchar(row$top_abs), row$top_abs, "NA")),
            sprintf("- **Top up genes:** %s", ifelse(nzchar(row$top_up), row$top_up, "NA")),
            sprintf("- **Top down genes:** %s", ifelse(nzchar(row$top_down), row$top_down, "NA")),
            sprintf("- **LLM overview:** %s", ifelse(nzchar(row$llm_overview), row$llm_overview, "NA")),
            sprintf("- **LLM key mechanisms:** %s", ifelse(nzchar(row$llm_key_mechanisms), row$llm_key_mechanisms, "NA")),
            sprintf("- **LLM narrative:** %s", ifelse(nzchar(row$llm_narrative), row$llm_narrative, "NA")),
            sprintf("- **LLM evidence:** %s", ifelse(nzchar(row$llm_evidence), row$llm_evidence, "NA")),
            sprintf("- **LLM key drivers:** %s", ifelse(nzchar(row$llm_key_drivers), row$llm_key_drivers, "NA")),
            sprintf("- **LLM limitations:** %s", ifelse(nzchar(row$llm_limitations), row$llm_limitations, "NA")),
            sprintf("- **Files:** DEG `%s`; LLM `%s`", row$deg_csv, row$llm_rds),
            ""
          )
        }
      }
    }
  }

  writeLines(enc2utf8(lines), out_path, useBytes = TRUE)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args) >= 1L && nzchar(trimws(args[[1]]))) args[[1]] else DEFAULT_OUTPUT_DIR
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "lineage_reports"), recursive = TRUE, showWarnings = FALSE)

  master_rows <- lapply(names(LINEAGE_DIRS), function(lineage) collect_lineage_records(lineage, LINEAGE_DIRS[[lineage]]))
  master_dt <- rbindlist(master_rows, fill = TRUE)
  if (!nrow(master_dt)) stop("No DEG/LLM records found in configured lineage directories.")
  setorder(master_dt, lineage, level, celltype, comparison)

  overview <- build_overview_tables(master_dt)

  fwrite(master_dt, file.path(output_dir, "all_lineages_tissue_deg_llm_master.tsv"), sep = "\t")
  fwrite(overview$lineage_level, file.path(output_dir, "all_lineages_tissue_deg_llm_lineage_level_summary.tsv"), sep = "\t")
  fwrite(overview$lineage_celltype, file.path(output_dir, "all_lineages_tissue_deg_llm_celltype_summary.tsv"), sep = "\t")
  fwrite(overview$top_contrasts, file.path(output_dir, "all_lineages_tissue_deg_llm_top_contrasts.tsv"), sep = "\t")

  for (lineage_name in unique(master_dt$lineage)) {
    lineage_dt <- master_dt[lineage == lineage_name]
    lineage_overview <- list(
      lineage_level = overview$lineage_level[lineage == lineage_name],
      lineage_celltype = overview$lineage_celltype[lineage == lineage_name],
      top_contrasts = overview$top_contrasts[lineage == lineage_name]
    )
    write_markdown_report(
      lineage_dt,
      lineage_overview,
      file.path(output_dir, "lineage_reports", sprintf("%s_tissue_deg_llm_detailed.md", lineage_name))
    )
  }

  write_markdown_report(master_dt, overview, file.path(output_dir, "all_lineages_tissue_deg_llm_detailed.md"))

  cat("[DONE] tissue DEG + LLM cross-lineage summary generated\n")
  cat(sprintf("Output dir: %s\n", normalizePath(output_dir, winslash = "/", mustWork = FALSE)))
  cat(sprintf("Master TSV : %s\n", file.path(output_dir, "all_lineages_tissue_deg_llm_master.tsv")))
  cat(sprintf("Markdown   : %s\n", file.path(output_dir, "all_lineages_tissue_deg_llm_detailed.md")))
}

main()
