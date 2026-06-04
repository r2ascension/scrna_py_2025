#!/usr/bin/env Rscript
# fgsea_tissue_deg_runner_20260602.R
# FGSEA batch runner for all tissue-comparison DEG tables across 7 lineages.
#
# For each DESeq2_results.csv (L2 + L3, all celltypes × contrasts):
#   1. Rank ALL genes by DESeq2 Wald statistic ("stat" column)
#   2. Run fgsea::fgsea() on MSigDB gene sets (Hallmark, KEGG, GO:BP/CC/MF)
#   3. Split results by NES sign → fgsea_up/ and fgsea_down/
#   4. Write TSV tables + dotplot PDFs
#
# Output location: <contrast_dir>/fgsea_up/   <contrast_dir>/fgsea_down/
#   alongside existing enrichment_up/ enrichment_down/.
#
# Resume-safe: skips contrasts that already have fgsea output.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(parallel)
})

# ── Configuration ───────────────────────────────────────────────────────────

# Gene set collections to run (msigdbr categories / subcategories)
# msigdbr >= 10.0.0 uses collection/subcollection (not category/subcategory)
GS_CONFIG <- list(
  Hallmark = list(collection = "H",  subcollection = NULL,            label = "MSigDB Hallmark"),
  KEGG     = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS", label = "KEGG"),
  GO_BP    = list(collection = "C5", subcollection = "GO:BP",          label = "GO Biological Process"),
  GO_CC    = list(collection = "C5", subcollection = "GO:CC",          label = "GO Cellular Component"),
  GO_MF    = list(collection = "C5", subcollection = "GO:MF",          label = "GO Molecular Function")
)

# Lineage root directories (latest versions from inventory 20260530)
LINEAGE_ROOTS <- list(
  epithelial  = "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508",
  bcell       = "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508",
  tnk         = "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508",
  myeloid     = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
  endothelial = "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508",
  fibroblast  = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414",
  smc         = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
)

FGSEA_N_PERM    <- 10000L
FGSEA_MIN_SIZE  <- 10L
FGSEA_MAX_SIZE  <- 500L
FGSEA_SEED      <- 20260602L
TOP_N_DOTPLOT   <- 20L
MIN_N_GENES     <- 50L   # skip contrasts with too few expressed genes

# ── Helpers ─────────────────────────────────────────────────────────────────

#' Build a named ranked vector from DESeq2 results using ALL genes
build_ranked_list <- function(de) {
  # de is a data.table from fread(DESeq2_results.csv)
  if (!"stat" %in% names(de)) {
    if (all(c("log2FoldChange", "lfcSE") %in% names(de))) {
      de[, stat := log2FoldChange / lfcSE]
    } else {
      stop("DESeq2_results.csv missing 'stat' column and cannot compute it")
    }
  }
  de <- de[!is.na(stat) & !is.infinite(stat)]
  # Sort descending by stat (positive stat = up in first condition)
  setorder(de, -stat)
  rnk <- de[["stat"]]
  names(rnk) <- de[["gene"]]
  rnk
}

#' Load and cache MSigDB gene set lists
load_gene_sets <- function(gs_config, species = "Homo sapiens") {
  gs_list <- list()
  for (nm in names(gs_config)) {
    cfg <- gs_config[[nm]]
    args <- list(species = species, collection = cfg$collection)
    if (!is.null(cfg$subcollection)) {
      args$subcollection <- cfg$subcollection
    }
    msig <- do.call(msigdbr, args)
    gs <- split(msig$gene_symbol, msig$gs_name)
    gs_list[[nm]] <- gs
    cat(sprintf("  %-8s: %d gene sets, %d unique genes\n",
                nm, length(gs), length(unique(msig$gene_symbol))))
  }
  gs_list
}

#' Run fgsea on ranked list, split results by NES sign
run_fgsea_directional <- function(rnk, gs_cache, nperm = FGSEA_N_PERM,
                                   min_size = FGSEA_MIN_SIZE,
                                   max_size = FGSEA_MAX_SIZE,
                                   seed = FGSEA_SEED) {
  # Returns list(up = list(...), down = list(...))
  results_up   <- list()
  results_down <- list()

  for (nm in names(gs_cache)) {
    set.seed(seed)
    res <- fgsea(
      pathways     = gs_cache[[nm]],
      stats        = rnk,
      minSize      = min_size,
      maxSize      = max_size,
      nPermSimple  = nperm
    )
    if (nrow(res) == 0L) next

    res <- as.data.table(res)
    # Flatten leadingEdge for TSV output
    res[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), "")]
    res[, collection := nm]

    # Split by NES direction
    res_up <- res[NES > 0]
    res_down <- res[NES < 0]

    if (nrow(res_up) > 0L) {
      setorder(res_up, pval)
      results_up[[nm]] <- res_up
    }
    if (nrow(res_down) > 0L) {
      setorder(res_down, pval)
      results_down[[nm]] <- res_down
    }
  }
  list(up = results_up, down = results_down)
}

#' Write FGSEA results to output directory with dotplots
write_fgsea_outputs <- function(results, out_dir, top_n = TOP_N_DOTPLOT) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  n_sets <- 0L

  for (nm in names(results)) {
    res <- results[[nm]]
    if (is.null(res) || nrow(res) == 0L) next

    # TSV
    fwrite(res, file.path(out_dir, paste0(nm, ".tsv")), sep = "\t")
    n_sets <- n_sets + 1L

    # Dotplot: top N by pval, show both up+down enriched together
    top_terms <- res[order(pval)][seq_len(min(top_n, nrow(res)))]
    top_terms <- top_terms[!is.na(pathway)]
    if (nrow(top_terms) == 0L) next

    top_terms[, `:=`(
      direction  = ifelse(NES > 0, "up", "down"),
      path_label = pathway
    )]

    p <- ggplot(top_terms, aes(x = NES, y = reorder(path_label, NES))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
      geom_point(aes(size = size, fill = direction), shape = 21, alpha = 0.85) +
      scale_fill_manual(values = c(up = "#E74C3C", down = "#2980B9")) +
      labs(
        title    = paste0(nm, " — Top ", nrow(top_terms)),
        subtitle = paste0("n = ", nrow(res), " enriched sets"),
        x        = "Normalized Enrichment Score (NES)",
        y        = "",
        size     = "Set Size",
        fill     = "Direction"
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 7, color = "grey50"),
        legend.position  = "bottom",
        panel.grid.major.y = element_line(linewidth = 0.2)
      )

    pdf_path <- file.path(out_dir, paste0(nm, "_dotplot.pdf"))
    ggsave(pdf_path, p, width = 10, height = max(4, nrow(top_terms) * 0.28),
           limitsize = FALSE)
  }
  n_sets
}

#' Process one DESeq2 table: one FGSEA run → split into fgsea_up / fgsea_down
process_one_contrast <- function(de_path, gs_cache) {
  contrast_dir <- dirname(de_path)
  contrast_name <- basename(contrast_dir)
  celltype_dir  <- dirname(contrast_dir)
  celltype_name <- basename(celltype_dir)

  result <- data.table(
    lineage       = NA_character_,
    level         = NA_character_,
    celltype      = celltype_name,
    contrast      = contrast_name,
    direction     = c("up", "down"),
    n_total_genes = 0L,
    n_sig_deg     = 0L,
    n_fgsea_sets  = 0L,
    status        = "pending",
    output_dir    = NA_character_
  )

  # Check skip conditions
  up_dir   <- file.path(contrast_dir, "fgsea_up")
  down_dir <- file.path(contrast_dir, "fgsea_down")
  up_done  <- length(list.files(up_dir,   pattern = "\\.tsv$")) >= length(gs_cache)
  down_done <- length(list.files(down_dir, pattern = "\\.tsv$")) >= length(gs_cache)

  if (up_done && down_done) {
    result[, status := "already_done"]
    result[, n_fgsea_sets := length(gs_cache)]
    result[ direction == "up",   output_dir := up_dir]
    result[ direction == "down", output_dir := down_dir]
    return(result[])
  }

  # Read DESeq2 results
  de <- fread(de_path)
  if (!"stat" %in% names(de)) {
    if (all(c("log2FoldChange", "lfcSE") %in% names(de))) {
      de[, stat := log2FoldChange / lfcSE]
    } else {
      result[, status := "error:no_stat_column"]
      return(result[])
    }
  }

  # Count sig DEGs per direction
  if ("sig" %in% names(de) && "direction" %in% names(de)) {
    n_up   <- de[sig == "sig" & direction == "up",   .N]
    n_down <- de[sig == "sig" & direction == "down", .N]
    result[ direction == "up",   n_sig_deg := n_up]
    result[ direction == "down", n_sig_deg := n_down]
  }

  # Build ranked list from ALL genes
  de_clean <- de[!is.na(stat) & !is.infinite(stat)]
  result[, n_total_genes := nrow(de_clean)]

  if (nrow(de_clean) < MIN_N_GENES) {
    result[, status := "skipped_too_few_genes"]
    return(result[])
  }

  setorder(de_clean, -stat)
  rnk <- de_clean[["stat"]]
  names(rnk) <- de_clean[["gene"]]

  # Run FGSEA once with full ranked list
  fgsea_res <- tryCatch({
    run_fgsea_directional(rnk, gs_cache)
  }, error = function(e) {
    result[, status := paste0("error:", conditionMessage(e))]
    return(NULL)
  })

  if (is.null(fgsea_res)) return(result[])

  # Write up direction
  if (!up_done) {
    result[ direction == "up", output_dir := up_dir]
    n_up_sets <- write_fgsea_outputs(fgsea_res$up, up_dir)
    result[direction == "up", `:=`(n_fgsea_sets = n_up_sets, status = "ok")]
  } else {
    result[direction == "up", `:=`(output_dir = up_dir, n_fgsea_sets = length(gs_cache), status = "already_done")]
  }

  # Write down direction
  if (!down_done) {
    result[ direction == "down", output_dir := down_dir]
    n_down_sets <- write_fgsea_outputs(fgsea_res$down, down_dir)
    result[direction == "down", `:=`(n_fgsea_sets = n_down_sets, status = "ok")]
  } else {
    result[direction == "down", `:=`(output_dir = down_dir, n_fgsea_sets = length(gs_cache), status = "already_done")]
  }

  result[]
}

# ── Main ────────────────────────────────────────────────────────────────────

main <- function() {
  cat("═══════════════════════════════════════════════\n")
  cat("  FGSEA Tissue DEG Batch Runner\n")
  cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("═══════════════════════════════════════════════\n\n")

  # Load gene sets once
  cat("─ Loading MSigDB gene sets ─\n")
  gs_cache <- load_gene_sets(GS_CONFIG)
  cat("\n")

  all_results <- list()
  total_tables <- 0L
  total_up_ok  <- 0L
  total_down_ok <- 0L
  start_time   <- Sys.time()

  for (lineage_name in names(LINEAGE_ROOTS)) {
    lineage_root <- LINEAGE_ROOTS[[lineage_name]]
    cat(sprintf("━━━ %s ━━━\n", toupper(lineage_name)))
    cat(sprintf("  %s\n", lineage_root))

    for (level_dir_name in c("pseudobulk_de", "pseudobulk_de_L3")) {
      level <- if (level_dir_name == "pseudobulk_de") "L2" else "L3"
      de_root <- file.path(lineage_root, level_dir_name)

      if (!dir.exists(de_root)) {
        cat(sprintf("  %s: directory not found\n", level))
        next
      }

      de_files <- list.files(
        de_root,
        pattern = "DESeq2_results\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
      de_files <- grep("/fgsea_", de_files, invert = TRUE, value = TRUE)

      cat(sprintf("  %s: %d contrast(s)\n", level, length(de_files)))
      total_tables <- total_tables + length(de_files)

      for (j in seq_along(de_files)) {
        de_path <- de_files[j]
        rel <- sub(paste0(lineage_root, "/"), "", de_path)
        rel <- sub("/DESeq2_results\\.csv$", "", rel)

        res <- process_one_contrast(de_path, gs_cache)
        res[, `:=`(lineage = lineage_name, level = level)]
        all_results[[length(all_results) + 1L]] <- res

        n_up_ok   <- sum(res[direction == "up",   status %in% c("ok", "already_done")])
        n_down_ok <- sum(res[direction == "down", status %in% c("ok", "already_done")])
        n_skip    <- sum(grepl("skipped", res$status))

        cat(sprintf("    [%3d/%3d] %s  (up=%d down=%d skip=%d)\n",
                    j, length(de_files), rel, n_up_ok, n_down_ok, n_skip))
        total_up_ok   <- total_up_ok + n_up_ok
        total_down_ok <- total_down_ok + n_down_ok
      }
    }
    cat("\n")
  }

  # ── Write manifest ──
  manifest <- rbindlist(all_results, fill = TRUE)
  out_base <- "/home/h2048/data/R/20260602/fgsea_tissue_deg_20260602"
  dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

  manifest_file <- file.path(out_base, "fgsea_manifest.tsv")
  fwrite(manifest, manifest_file, sep = "\t")

  # ── Summary ──
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat("═══════════════════════════════════════════════\n")
  cat("  SUMMARY\n")
  cat("═══════════════════════════════════════════════\n")
  cat(sprintf("  Total contrasts       : %d\n", total_tables))
  cat(sprintf("  Manifest rows         : %d\n", nrow(manifest)))
  cat(sprintf("  FGSEA up   completed  : %d\n", total_up_ok))
  cat(sprintf("  FGSEA down completed  : %d\n", total_down_ok))
  cat(sprintf("  Elapsed               : %.1f min\n", elapsed))
  cat(sprintf("  Manifest              : %s\n", manifest_file))
  cat("\n  Status breakdown:\n")
  tab <- manifest[, .N, by = .(status)]
  print(tab)
  cat(sprintf("\n  Done: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

main()
