#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

save_plot <- function(plot_obj, path_no_ext, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), plot_obj, width = width, height = height, dpi = dpi)
  invisible(list(pdf = paste0(path_no_ext, ".pdf"), png = paste0(path_no_ext, ".png")))
}

is_meaningful_gene <- function(x) {
  x <- as.character(x)
  !grepl("^(ENSG|LOC[0-9]|CIM|SPMIP|HOATZ)", x, perl = TRUE)
}

build_labels <- function(df, forced_genes, top_n = 12) {
  df <- df %>%
    mutate(
      gene = as.character(gene),
      padj_plot = ifelse(is.na(padj) | padj <= 0, NA_real_, -log10(padj)),
      sig_flag = ifelse("sig" %in% names(.), sig == "sig", !is.na(padj) & padj < 0.05),
      forced_flag = gene %in% forced_genes,
      meaningful_flag = is_meaningful_gene(gene)
    )

  top_genes <- df %>%
    filter(sig_flag, meaningful_flag) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n) %>%
    pull(gene) %>%
    unique()

  df %>%
    mutate(
      label = ifelse(forced_flag | gene %in% top_genes, gene, ""),
      point_class = case_when(
        forced_flag ~ "forced",
        sig_flag ~ "sig",
        TRUE ~ "ns"
      )
    )
}

plot_one <- function(csv_path, out_stem, title, forced_genes, padj_thr = 0.05, lfc_thr = 0.15) {
  df <- fread(csv_path)
  stopifnot(all(c("gene", "log2FoldChange", "padj") %in% names(df)))

  df2 <- build_labels(df, forced_genes = forced_genes, top_n = 12)

  p <- ggplot(df2, aes(x = log2FoldChange, y = padj_plot)) +
    geom_point(aes(color = point_class), alpha = 0.55, size = 0.9, na.rm = TRUE) +
    scale_color_manual(
      values = c(ns = "grey75", sig = "#e34a33", forced = "#2166ac"),
      breaks = c("forced", "sig", "ns"),
      labels = c("Forced labels", "Significant", "Non-significant"),
      name = NULL
    ) +
    geom_text_repel(
      data = df2 %>% filter(nzchar(label), !is.na(padj_plot)),
      aes(label = label),
      size = 3,
      max.overlaps = 40,
      min.segment.length = 0,
      seed = 42,
      box.padding = 0.3,
      point.padding = 0.15
    ) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "blue") +
    labs(
      title = title,
      subtitle = paste("Forced genes:", paste(forced_genes, collapse = ", ")),
      x = "log2 Fold Change",
      y = "-log10(padj)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  save_plot(p, out_stem, width = 8.2, height = 6.4, dpi = 320)
}

main <- function() {
  input_dir <- Sys.getenv("VOLCANO_INPUT_DIR", unset = ".")
  output_dir <- Sys.getenv("VOLCANO_OUTPUT_DIR", unset = file.path(input_dir, "annotated"))

  jobs <- list(
    list(
      csv = file.path(input_dir, "basal_lineage_respiratory_airway_vs_nose_DESeq2_results.csv"),
      stem = file.path(output_dir, "basal_lineage_respiratory_airway_vs_nose_annotated_volcano"),
      title = "Basal lineage: bronchus vs nose",
      forced = c("IRX2", "HOTAIRM1", "AQP4-AS1", "NKX2-1", "LYPD2", "C9orf24")
    ),
    list(
      csv = file.path(input_dir, "ciliated_lineage_sinus_vs_respiratory_airway_DESeq2_results.csv"),
      stem = file.path(output_dir, "ciliated_lineage_sinus_vs_respiratory_airway_annotated_volcano"),
      title = "Ciliated lineage: sinus vs bronchus",
      forced = c("C20orf85", "C9orf24", "PIFO", "IRX2", "HOTAIRM1", "AQP4-AS1", "CFAP144", "CFAP276")
    ),
    list(
      csv = file.path(input_dir, "secretory_lineage_sinus_vs_nose_DESeq2_results.csv"),
      stem = file.path(output_dir, "secretory_lineage_sinus_vs_nose_annotated_volcano"),
      title = "Secretory lineage: sinus vs nose",
      forced = c("AQP4-AS1", "C20orf85", "HOTAIRM1", "C9orf24", "PIFO", "MEG3", "IRX2")
    )
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (job in jobs) {
    message("Rendering: ", basename(job$csv))
    plot_one(
      csv_path = job$csv,
      out_stem = job$stem,
      title = job$title,
      forced_genes = job$forced
    )
  }
}

if (sys.nframe() == 0) {
  main()
}
