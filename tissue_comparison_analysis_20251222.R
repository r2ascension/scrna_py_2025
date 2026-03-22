# =============================================================================
# Production-Ready Tissue Comparison Analysis
# =============================================================================
# Version: 3.0 Final (Offline GMT + All Stability Fixes)
#
# Key Features:
# - Completely offline (uses local msigdb GMT file)
# - Fixed DESeq2 factor levels issue
# - GSVA on pseudobulk groups (not celltype means)
# - Enhanced error handling and logging
# - All English labels and comments
#
# Author: Clinical Bioinformatics Team
# Date: 2024-12-23
# =============================================================================

#' Load MSigDB Gene Sets from Local GMT File
#'
#' @param gmt_file Path to local msigdb GMT file
#' @param gene_universe Gene names to filter against (e.g., rownames(expr))
#' @param min_size Minimum gene set size (default: 10)
#' @param max_size Maximum gene set size (default: 5000)
#'
#' @return List with hallmark, kegg_medicus, and kegg_legacy gene sets
#'
#' @details
#' Parses local GMT file and extracts:
#' - HALLMARK_* gene sets
#' - KEGG_MEDICUS_* gene sets (new naming in v2025.1)
#' - Legacy KEGG_* gene sets (fallback)
#'
#' @export
load_msigdb_subsets_from_gmt <- function(
  gmt_file,
  gene_universe = NULL,
  min_size = 10,
  max_size = 5000
) {
  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Please install fgsea: BiocManager::install('fgsea')")
  }

  message(sprintf("Parsing GMT file: %s", basename(gmt_file)))
  gs_all <- fgsea::gmtPathways(gmt_file)
  message(sprintf("Total gene sets in GMT: %d", length(gs_all)))

  # Extract subsets by name pattern
  hallmark <- gs_all[grepl("^HALLMARK_", names(gs_all))]
  kegg_medicus <- gs_all[grepl("^KEGG_MEDICUS_", names(gs_all))]
  kegg_legacy <- gs_all[
    grepl("^KEGG_", names(gs_all)) & !grepl("^KEGG_MEDICUS_", names(gs_all))
  ]

  message(sprintf("  Hallmark: %d", length(hallmark)))
  message(sprintf("  KEGG_MEDICUS: %d", length(kegg_medicus)))
  message(sprintf("  KEGG_LEGACY: %d", length(kegg_legacy)))

  # Filter by gene universe and size
  filter_sets <- function(sets, label) {
    if (!is.null(gene_universe)) {
      sets <- lapply(sets, function(x) intersect(unique(x), gene_universe))
    }

    sz <- lengths(sets)
    sets_filtered <- sets[sz >= min_size & sz <= max_size]

    if (length(sets) > 0) {
      message(sprintf(
        "  %s: %d → %d after filtering (size %d-%d)",
        label,
        length(sets),
        length(sets_filtered),
        min_size,
        max_size
      ))
    }

    sets_filtered
  }

  list(
    hallmark = filter_sets(hallmark, "Hallmark"),
    kegg_medicus = filter_sets(kegg_medicus, "KEGG_MEDICUS"),
    kegg_legacy = filter_sets(kegg_legacy, "KEGG_LEGACY")
  )
}


#' Tissue-Specific Differential Analysis and Pathway Activity Assessment
#'
#' @description
#' Production-ready tissue comparison analysis (OFFLINE VERSION):
#' - Uses local msigdb GMT file (no internet required)
#' - Robust pseudobulk aggregation (no parsing bugs)
#' - True parallel processing with BiocParallel
#' - GSVA on pseudobulk groups (handles single cell type)
#' - Fixed DESeq2 factor levels issue
#' - Enhanced QC and error handling
#'
#' @param seurat_obj Seurat object containing single-cell data
#' @param cell_anno_col Cell type annotation column (default: "Annotation_2")
#' @param tissue_col Tissue type column (default: "tissue")
#' @param sample_col Sample ID column (default: "sample")
#' @param min_cell_per_sample Minimum cells per pseudobulk (default: 3)
#' @param min_sample_per_tissue Minimum samples per tissue (default: 3)
#' @param min_pb_libsize Minimum pseudobulk library size (default: 1000)
#' @param min_pb_detected_genes Minimum detected genes per pseudobulk (default: 500)
#' @param run_gsva Perform GSVA analysis (default: TRUE)
#' @param run_go Perform GO enrichment (default: TRUE)
#' @param output_dir Output directory (default: "./analysis_results")
#' @param species Species for GO enrichment (default: "Homo sapiens")
#' @param msigdb_gmt Path to local msigdb GMT file (REQUIRED if run_gsva=TRUE)
#' @param padj_thr Adjusted p-value threshold (default: 0.05)
#' @param lfc_thr Log2 fold change threshold (default: 1.0)
#' @param min_gene_total_counts Minimum gene total counts (default: 10)
#' @param top_deg_heatmap Top N DEGs for heatmap (default: 50)
#' @param gsva_top_var_h Top variable Hallmark pathways (default: 20)
#' @param gsva_top_var_k Top variable KEGG pathways (default: 30)
#' @param n_cores CPU cores for parallel processing (default: 1)
#' @param use_lfc_shrink Use apeglm LFC shrinkage (default: TRUE)
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Download GMT file first (do once):
#' # wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/msigdb.v2025.1.Hs.symbols.gmt
#'
#' run_tissue_comparison_analysis(
#'   seurat_obj = T_object,
#'   msigdb_gmt = "/path/to/msigdb.v2025.1.Hs.symbols.gmt",
#'   n_cores = 48,
#'   output_dir = "./results"
#' )
#' }
#'
#' @export
run_tissue_comparison_analysis <- function(
  seurat_obj,
  cell_anno_col = "Annotation_2",
  tissue_col = "tissue",
  sample_col = "sample",
  min_cell_per_sample = 3,
  min_sample_per_tissue = 3,
  min_pb_libsize = 1000,
  min_pb_detected_genes = 500,
  run_gsva = TRUE,
  run_go = TRUE,
  output_dir = "./analysis_results",
  species = "Homo sapiens",
  msigdb_gmt = NULL,
  padj_thr = 0.05,
  lfc_thr = 1.0,
  min_gene_total_counts = 10,
  top_deg_heatmap = 50,
  gsva_top_var_h = 20,
  gsva_top_var_k = 30,
  n_cores = 1,
  use_lfc_shrink = TRUE
) {
  # ===== 1. Dependency Check =====
  required_pkgs <- c(
    "Seurat",
    "DESeq2",
    "ggplot2",
    "pheatmap",
    "dplyr",
    "fgsea",
    "GSVA",
    "tidyr",
    "BiocParallel"
  )

  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Please install missing packages: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }

  # Check LFC shrinkage
  if (use_lfc_shrink && !requireNamespace("apeglm", quietly = TRUE)) {
    message("⚠️ apeglm not installed, disabling LFC shrinkage")
    use_lfc_shrink <- FALSE
  }

  # Check GO enrichment
  if (run_go) {
    go_pkgs <- c("clusterProfiler", "org.Hs.eg.db")
    if (!all(vapply(go_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
      message(
        "⚠️ clusterProfiler or org.Hs.eg.db not installed, skipping GO enrichment"
      )
      run_go <- FALSE
    }
  }

  # Check GSVA GMT file
  if (run_gsva) {
    if (is.null(msigdb_gmt) || !nzchar(msigdb_gmt)) {
      stop(
        "run_gsva=TRUE but msigdb_gmt is NULL/empty. Please provide local GMT file path."
      )
    }
    if (!file.exists(msigdb_gmt)) {
      stop("msigdb_gmt file not found: ", msigdb_gmt)
    }
    message(sprintf("✓ Using local GMT file: %s", basename(msigdb_gmt)))
  }

  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  # Check metadata columns
  metadata <- seurat_obj@meta.data
  required_cols <- c(cell_anno_col, tissue_col, sample_col)
  missing_cols <- required_cols[!required_cols %in% colnames(metadata)]

  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  # ===== 2. Setup Parallel Processing =====
  if (n_cores > 1) {
    message(sprintf("Setting up parallel processing with %d cores", n_cores))

    BPPARAM <- if (.Platform$OS.type == "unix") {
      BiocParallel::MulticoreParam(workers = n_cores)
    } else {
      BiocParallel::SnowParam(workers = n_cores)
    }

    BiocParallel::register(BPPARAM)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }

  # ===== 3. Directory Setup =====
  out_root <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  pb_dir <- file.path(out_root, "pseudobulk_analysis")
  de_dir <- file.path(pb_dir, "DESeq2_pairwise")
  sum_dir <- file.path(pb_dir, "summary")
  qc_dir <- file.path(pb_dir, "qc_metrics")
  gsva_dir <- file.path(out_root, "gsva_analysis")

  dir.create(pb_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  if (run_gsva) {
    dir.create(gsva_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # ===== 4. Dataset Overview =====
  message("==== Dataset Overview ====")
  message(sprintf("Total cells: %d", ncol(seurat_obj)))
  message(sprintf("Total genes: %d", nrow(seurat_obj)))
  message("\nTissue distribution:")
  print(table(metadata[[tissue_col]]))
  message("\nCell type distribution:")
  print(table(metadata[[cell_anno_col]]))

  # ===== 5. Robust Pseudobulk ID Generation =====
  # Clean and validate metadata columns
  metadata[[tissue_col]] <- trimws(as.character(metadata[[tissue_col]]))
  metadata[[sample_col]] <- trimws(as.character(metadata[[sample_col]]))
  metadata[[cell_anno_col]] <- trimws(as.character(metadata[[cell_anno_col]]))

  # Remove cells with NA or empty values (critical for tissue comparisons)
  bad_meta <- is.na(metadata[[tissue_col]]) |
    metadata[[tissue_col]] == "" |
    is.na(metadata[[sample_col]]) |
    metadata[[sample_col]] == "" |
    is.na(metadata[[cell_anno_col]]) |
    metadata[[cell_anno_col]] == ""

  if (any(bad_meta)) {
    message(sprintf(
      "⚠️ Dropping %d cells with NA/empty tissue/sample/celltype",
      sum(bad_meta)
    ))
    seurat_obj <- Seurat::subset(
      seurat_obj,
      cells = rownames(metadata)[!bad_meta]
    )
    metadata <- seurat_obj@meta.data

    # Re-clean after subset
    metadata[[tissue_col]] <- trimws(as.character(metadata[[tissue_col]]))
    metadata[[sample_col]] <- trimws(as.character(metadata[[sample_col]]))
    metadata[[cell_anno_col]] <- trimws(as.character(metadata[[cell_anno_col]]))
  }

  sep_token <- "|||"
  metadata$.tissue <- metadata[[tissue_col]]
  metadata$.sample <- metadata[[sample_col]]
  metadata$.celltype <- metadata[[cell_anno_col]]
  metadata$.pb_id <- paste(
    metadata$.tissue,
    metadata$.sample,
    metadata$.celltype,
    sep = sep_token
  )

  seurat_obj$.pb_id <- metadata$.pb_id

  # ===== 6. Enhanced QC =====
  message("\n==== Pseudobulk Quality Control ====")

  pb_cell_counts <- table(metadata$.pb_id)
  message(sprintf(
    "Total pseudobulk groups before QC: %d",
    length(pb_cell_counts)
  ))

  valid_by_cells <- names(pb_cell_counts[pb_cell_counts >= min_cell_per_sample])
  message(sprintf(
    "Groups with >= %d cells: %d",
    min_cell_per_sample,
    length(valid_by_cells)
  ))

  if (length(valid_by_cells) == 0) {
    stop("All pseudobulk groups have < min_cell_per_sample cells")
  }

  keep_cells <- rownames(metadata)[metadata$.pb_id %in% valid_by_cells]
  seurat_obj_filtered <- subset(seurat_obj, cells = keep_cells)
  metadata_filtered <- seurat_obj_filtered@meta.data

  dist_table <- table(
    metadata_filtered[[tissue_col]],
    metadata_filtered[[cell_anno_col]]
  )
  utils::write.csv(
    dist_table,
    file.path(qc_dir, "tissue_celltype_distribution_cells.csv")
  )

  # ===== 7. Pseudobulk Aggregation =====
  message("\n==== Pseudobulk Aggregation ====")

  aggregate_counts_safe <- function(obj, group.by, assay = "RNA") {
    agg_formals <- names(formals(Seurat::AggregateExpression))

    args <- list(
      object = obj,
      group.by = group.by,
      assays = assay,
      return.seurat = FALSE,
      verbose = FALSE
    )

    if ("layer" %in% agg_formals) {
      args$layer <- "counts"
    }
    if ("slot" %in% agg_formals) {
      args$slot <- "counts"
    }

    result <- do.call(Seurat::AggregateExpression, args)
    result[[assay]]
  }

  pb_counts <- aggregate_counts_safe(seurat_obj_filtered, group.by = ".pb_id")
  message(sprintf(
    "Pseudobulk matrix: %d genes × %d groups",
    nrow(pb_counts),
    ncol(pb_counts)
  ))

  # ===== 8. Additional QC: Library Size and Gene Detection =====
  pb_metadata <- unique(metadata_filtered[, c(
    ".pb_id",
    tissue_col,
    sample_col,
    cell_anno_col
  )])
  colnames(pb_metadata) <- c("pb_id", "tissue", "sample", "celltype")
  pb_metadata$pb_id <- as.character(pb_metadata$pb_id)
  rownames(pb_metadata) <- pb_metadata$pb_id
  pb_metadata <- pb_metadata[colnames(pb_counts), , drop = FALSE]

  # Convert to character to avoid factor issues in DESeq2
  pb_metadata$tissue <- as.character(pb_metadata$tissue)
  pb_metadata$sample <- as.character(pb_metadata$sample)
  pb_metadata$celltype <- as.character(pb_metadata$celltype)

  # Calculate QC metrics
  pb_metadata$library_size <- colSums(pb_counts)
  pb_metadata$detected_genes <- colSums(pb_counts > 0)
  pb_metadata$n_cells <- as.numeric(pb_cell_counts[pb_metadata$pb_id])

  # Filter by library size and gene detection
  qc_pass <- (pb_metadata$library_size >= min_pb_libsize) &
    (pb_metadata$detected_genes >= min_pb_detected_genes)

  message(sprintf(
    "QC filtering: library size >= %d, detected genes >= %d",
    min_pb_libsize,
    min_pb_detected_genes
  ))
  message(sprintf(
    "Groups passing all QC: %d / %d (%.1f%%)",
    sum(qc_pass),
    length(qc_pass),
    100 * mean(qc_pass)
  ))

  # Apply QC filter
  pb_counts <- pb_counts[, qc_pass, drop = FALSE]
  pb_metadata <- pb_metadata[qc_pass, , drop = FALSE]

  # Save QC metrics
  utils::write.csv(
    pb_metadata,
    file.path(qc_dir, "pseudobulk_qc_metrics.csv"),
    row.names = FALSE
  )

  # QC distribution plots
  pdf(file.path(qc_dir, "qc_distributions.pdf"), width = 12, height = 8)
  par(mfrow = c(2, 2))
  hist(
    pb_metadata$n_cells,
    breaks = 30,
    main = "Cells per Pseudobulk",
    xlab = "Number of cells",
    col = "steelblue"
  )
  abline(v = min_cell_per_sample, col = "red", lty = 2, lwd = 2)

  hist(
    log10(pb_metadata$library_size),
    breaks = 30,
    main = "Library Size (log10)",
    xlab = "log10(Total counts)",
    col = "steelblue"
  )
  abline(v = log10(min_pb_libsize), col = "red", lty = 2, lwd = 2)

  hist(
    pb_metadata$detected_genes,
    breaks = 30,
    main = "Detected Genes",
    xlab = "Number of genes",
    col = "steelblue"
  )
  abline(v = min_pb_detected_genes, col = "red", lty = 2, lwd = 2)

  boxplot(
    library_size ~ celltype,
    data = pb_metadata,
    las = 2,
    main = "Library Size by Cell Type",
    ylab = "Total counts",
    col = "lightblue"
  )
  dev.off()

  # Save pseudobulk data
  saveRDS(pb_counts, file.path(pb_dir, "pseudobulk_counts_matrix.rds"))
  saveRDS(pb_metadata, file.path(pb_dir, "pseudobulk_metadata.rds"))
  utils::write.csv(
    pb_metadata,
    file.path(pb_dir, "pseudobulk_metadata.csv"),
    row.names = FALSE
  )

  message(sprintf("Final pseudobulk groups: %d", nrow(pb_metadata)))
  message(sprintf("  Unique tissues: %d", length(unique(pb_metadata$tissue))))
  message(sprintf("  Unique samples: %d", length(unique(pb_metadata$sample))))
  message(sprintf(
    "  Unique cell types: %d",
    length(unique(pb_metadata$celltype))
  ))

  # ===== 9. Helper Functions =====

  sanitize_name <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
  }

  make_volcano_plot <- function(res_df, title, output_pdf, padj_thr, lfc_thr) {
    res_df$padj_plot <- pmax(res_df$padj, .Machine$double.xmin)

    p <- ggplot2::ggplot(
      res_df,
      ggplot2::aes(x = log2FoldChange, y = -log10(padj_plot))
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = regulation),
        size = 1.2,
        alpha = 0.7
      ) +
      ggplot2::scale_color_manual(
        values = c(down = "#00468B", stable = "gray70", up = "#E64B35"),
        name = "Regulation"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(padj_thr),
        linetype = "dashed",
        linewidth = 0.4,
        color = "black"
      ) +
      ggplot2::geom_vline(
        xintercept = c(-lfc_thr, lfc_thr),
        linetype = "dashed",
        linewidth = 0.4,
        color = "black"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        )
      ) +
      ggplot2::labs(
        x = "Log2(Fold Change)",
        y = "-log10(Adjusted p-value)",
        title = title
      )

    ggplot2::ggsave(
      output_pdf,
      plot = p,
      width = 10,
      height = 8,
      device = "pdf"
    )
    invisible(p)
  }

  make_deg_heatmap <- function(counts_mat, col_anno, genes, title, output_pdf) {
    if (length(genes) == 0) {
      return(invisible(NULL))
    }

    genes <- intersect(genes, rownames(counts_mat))
    if (length(genes) == 0) {
      return(invisible(NULL))
    }

    expr <- counts_mat[genes, , drop = FALSE]
    expr <- log2(expr + 1)
    expr <- t(scale(t(expr)))
    expr[!is.finite(expr)] <- 0

    pheatmap::pheatmap(
      expr,
      annotation_col = col_anno,
      show_colnames = FALSE,
      fontsize_row = 7,
      main = title,
      filename = output_pdf,
      width = 12,
      height = max(6, 0.12 * nrow(expr) + 4),
      color = colorRampPalette(c("#00468B", "white", "#E64B35"))(100)
    )

    invisible(NULL)
  }

  run_go_enrichment <- function(
    gene_vec,
    tested_genes,
    out_csv,
    out_pdf,
    title
  ) {
    if (!run_go) {
      return(invisible(NULL))
    }
    if (length(gene_vec) < 10) {
      message(sprintf("Too few genes (%d) for GO enrichment", length(gene_vec)))
      return(invisible(NULL))
    }

    gene_vec <- unique(gene_vec)
    tested_genes <- unique(tested_genes)

    tryCatch(
      {
        is_ensembl <- mean(grepl("^ENSG", gene_vec)) > 0.5
        key_type <- if (is_ensembl) "ENSEMBL" else "SYMBOL"

        if (is_ensembl) {
          gene_vec <- sub("\\..*", "", gene_vec)
          tested_genes <- sub("\\..*", "", tested_genes)
        }

        message(sprintf("GO enrichment using %s IDs", key_type))

        ego <- clusterProfiler::enrichGO(
          gene = gene_vec,
          universe = tested_genes,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
          keyType = key_type,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )

        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
          message("No significant GO terms found")
          return(invisible(NULL))
        }

        utils::write.csv(as.data.frame(ego), out_csv, row.names = FALSE)

        p <- clusterProfiler::dotplot(ego, showCategory = 20, title = title)
        ggplot2::ggsave(
          out_pdf,
          plot = p,
          width = 10,
          height = 8,
          device = "pdf"
        )

        message(sprintf(
          "GO enrichment: %d significant terms",
          nrow(as.data.frame(ego))
        ))
        invisible(ego)
      },
      error = function(e) {
        message(sprintf("GO enrichment error: %s", e$message))
        invisible(NULL)
      }
    )
  }

  # ===== 10. DESeq2 Pairwise Comparisons =====
  message("\n==== DESeq2 Pairwise Tissue Comparisons ====")

  all_celltypes <- sort(unique(pb_metadata$celltype))
  all_celltypes <- all_celltypes[!is.na(all_celltypes) & all_celltypes != ""]

  message(sprintf("Analyzing %d cell types", length(all_celltypes)))

  summary_rows <- list()

  for (ct in all_celltypes) {
    message(sprintf("\n--- Processing cell type: %s ---", ct))

    ct_indices <- which(pb_metadata$celltype == ct)

    if (length(ct_indices) < (min_sample_per_tissue * 2)) {
      message(sprintf(
        "[SKIP] %s: Insufficient samples (%d)",
        ct,
        length(ct_indices)
      ))
      next
    }

    ct_counts <- pb_counts[, ct_indices, drop = FALSE]
    ct_meta <- pb_metadata[ct_indices, , drop = FALSE]

    tissue_counts <- table(ct_meta$tissue)
    message(sprintf("Tissue distribution for %s:", ct))
    print(tissue_counts)

    valid_tissues <- names(tissue_counts[
      tissue_counts >= min_sample_per_tissue
    ])

    if (length(valid_tissues) < 2) {
      message(sprintf(
        "[SKIP] %s: < 2 tissues with >= %d samples",
        ct,
        min_sample_per_tissue
      ))
      next
    }

    keep_samples <- ct_meta$tissue %in% valid_tissues
    ct_counts <- ct_counts[, keep_samples, drop = FALSE]
    ct_meta <- ct_meta[keep_samples, , drop = FALSE]

    tissue_pairs <- combn(valid_tissues, 2, simplify = FALSE)
    message(sprintf("Performing %d pairwise comparisons", length(tissue_pairs)))

    for (pair in tissue_pairs) {
      tissue1 <- pair[1]
      tissue2 <- pair[2]

      message(sprintf("  Comparing: %s vs %s", tissue1, tissue2))

      pair_samples <- ct_meta$tissue %in% c(tissue1, tissue2)
      pair_counts <- ct_counts[, pair_samples, drop = FALSE]
      pair_meta <- ct_meta[pair_samples, , drop = FALSE]

      # FIX: Force character to avoid factor level issues
      pair_tissue_counts <- table(as.character(pair_meta$tissue))

      if (any(pair_tissue_counts < min_sample_per_tissue)) {
        message(sprintf(
          "  [SKIP] %s vs %s: insufficient samples: %s",
          tissue1,
          tissue2,
          paste(
            names(pair_tissue_counts),
            pair_tissue_counts,
            sep = "=",
            collapse = ", "
          )
        ))
        next
      }

      pair_counts <- round(as.matrix(pair_counts))
      mode(pair_counts) <- "integer"

      tryCatch(
        {
          dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = pair_counts,
            colData = data.frame(pair_meta),
            design = ~tissue
          )

          dds$tissue <- stats::relevel(factor(dds$tissue), ref = tissue1)

          keep_genes <- rowSums(DESeq2::counts(dds)) >= min_gene_total_counts
          dds <- dds[keep_genes, ]

          message(sprintf("  Testing %d genes", nrow(dds)))

          # Run DESeq2 with explicit fitType to avoid warnings
          dds <- DESeq2::DESeq(
            dds,
            parallel = (n_cores > 1),
            BPPARAM = BPPARAM,
            fitType = "local",
            quiet = TRUE
          )

          res <- DESeq2::results(dds, contrast = c("tissue", tissue2, tissue1))

          # CRITICAL: Create res_df FIRST (ensures it always exists)
          res_df <- as.data.frame(res)
          res_df$gene <- rownames(res_df)

          # OPTIONAL: Apply LFC shrinkage (modifies res_df in place if successful)
          if (use_lfc_shrink) {
            # Auto-detect coefficient name (handles spaces/special chars in tissue names)
            rn <- DESeq2::resultsNames(dds)

            # Build pattern using make.names() to match DESeq2's internal naming
            tissue2_safe <- make.names(tissue2)
            tissue1_safe <- make.names(tissue1)
            coef_pattern <- paste0(
              "^tissue_",
              tissue2_safe,
              "_vs_",
              tissue1_safe,
              "$"
            )

            coef_matches <- grep(coef_pattern, rn, value = TRUE)

            if (length(coef_matches) == 1) {
              shrink_result <- tryCatch(
                {
                  DESeq2::lfcShrink(
                    dds,
                    coef = coef_matches[1],
                    type = "apeglm",
                    quiet = TRUE
                  )
                },
                error = function(e) NULL
              )

              if (!is.null(shrink_result)) {
                res_df$log2FoldChange_original <- res_df$log2FoldChange
                res_df$log2FoldChange <- shrink_result$log2FoldChange
                message("  LFC shrinkage applied successfully")
              } else {
                message(
                  "  LFC shrinkage failed, using unshrunken log2FoldChange"
                )
              }
            } else {
              message(sprintf(
                "  LFC shrinkage skipped: found %d matching coefficients (expected 1)",
                length(coef_matches)
              ))
              if (length(coef_matches) > 0) {
                message(sprintf(
                  "  Candidates: %s",
                  paste(coef_matches, collapse = ", ")
                ))
              }
            }
          }

          res_df$gene <- rownames(res_df)
          res_df <- res_df[!is.na(res_df$padj), ]
          res_df <- res_df[order(res_df$padj), ]

          # Classify regulation
          res_df$regulation <- "stable"
          sig_up <- (res_df$padj <= padj_thr) &
            (res_df$log2FoldChange >= lfc_thr)
          sig_down <- (res_df$padj <= padj_thr) &
            (res_df$log2FoldChange <= -lfc_thr)
          res_df$regulation[sig_up] <- "up"
          res_df$regulation[sig_down] <- "down"

          n_up <- sum(res_df$regulation == "up")
          n_down <- sum(res_df$regulation == "down")
          n_total_deg <- n_up + n_down

          message(sprintf(
            "  DEGs: %d total (%d up, %d down)",
            n_total_deg,
            n_up,
            n_down
          ))

          comp_name <- paste0(tissue2, "_vs_", tissue1)
          out_subdir <- file.path(
            de_dir,
            sanitize_name(ct),
            sanitize_name(comp_name)
          )
          dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)

          utils::write.csv(
            res_df,
            file.path(out_subdir, "DEGs.csv"),
            row.names = FALSE
          )

          stat_row <- data.frame(
            CellType = ct,
            Comparison = comp_name,
            Tissue_1 = tissue1,
            Tissue_2 = tissue2,
            Total_DEGs = n_total_deg,
            Up_regulated = n_up,
            Down_regulated = n_down,
            Total_genes = nrow(res_df),
            stringsAsFactors = FALSE
          )

          summary_rows[[length(summary_rows) + 1]] <- stat_row
          utils::write.csv(
            stat_row,
            file.path(out_subdir, "summary_stats.csv"),
            row.names = FALSE
          )

          make_volcano_plot(
            res_df,
            sprintf("%s: %s vs %s", ct, tissue2, tissue1),
            file.path(out_subdir, "volcano_plot.pdf"),
            padj_thr,
            lfc_thr
          )

          if (n_total_deg > 0) {
            top_deg_genes <- res_df %>%
              dplyr::filter(regulation != "stable") %>%
              dplyr::arrange(padj) %>%
              dplyr::pull(gene) %>%
              head(top_deg_heatmap)

            anno_col <- data.frame(Tissue = pair_meta$tissue)
            rownames(anno_col) <- rownames(pair_meta)

            make_deg_heatmap(
              pair_counts,
              anno_col,
              top_deg_genes,
              sprintf("%s DEGs: %s vs %s", ct, tissue2, tissue1),
              file.path(out_subdir, "DEG_heatmap.pdf")
            )
          }

          # GO enrichment
          if (run_go && n_total_deg > 0) {
            tested_genes <- rownames(dds)

            up_genes <- res_df$gene[res_df$regulation == "up"]
            if (length(up_genes) >= 10) {
              run_go_enrichment(
                up_genes,
                tested_genes,
                file.path(out_subdir, "GO_upregulated.csv"),
                file.path(out_subdir, "GO_upregulated_dotplot.pdf"),
                "GO Enrichment: Up-regulated Genes"
              )
            }

            down_genes <- res_df$gene[res_df$regulation == "down"]
            if (length(down_genes) >= 10) {
              run_go_enrichment(
                down_genes,
                tested_genes,
                file.path(out_subdir, "GO_downregulated.csv"),
                file.path(out_subdir, "GO_downregulated_dotplot.pdf"),
                "GO Enrichment: Down-regulated Genes"
              )
            }
          }
        },
        error = function(e) {
          message(sprintf("  [ERROR] DESeq2 failed: %s", e$message))
        }
      )
    }
  }

  # ===== 11. Summary Visualizations =====
  message("\n==== Creating Summary Visualizations ====")

  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)

    saveRDS(summary_df, file.path(sum_dir, "all_comparisons_summary.rds"))
    utils::write.csv(
      summary_df,
      file.path(sum_dir, "all_comparisons_summary.csv"),
      row.names = FALSE
    )

    message(sprintf(
      "Summary: %d comparisons across %d cell types",
      nrow(summary_df),
      length(unique(summary_df$CellType))
    ))

    # Bar plot
    p1 <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = CellType, y = Total_DEGs, fill = Comparison)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      ) +
      ggplot2::labs(
        title = "DEG Counts by Cell Type and Tissue Comparison",
        x = "Cell Type",
        y = "Number of DEGs",
        fill = "Comparison"
      )

    ggplot2::ggsave(
      file.path(sum_dir, "DEG_counts_barplot.pdf"),
      p1,
      width = 15,
      height = 10,
      device = "pdf"
    )

    # Stacked bar plot
    summary_long <- summary_df %>%
      tidyr::pivot_longer(
        cols = c("Up_regulated", "Down_regulated"),
        names_to = "Direction",
        values_to = "Count"
      )

    p2 <- ggplot2::ggplot(
      summary_long,
      ggplot2::aes(x = Comparison, y = Count, fill = Direction)
    ) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::facet_wrap(~CellType, scales = "free_y", ncol = 3) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      ) +
      ggplot2::scale_fill_manual(
        values = c(Up_regulated = "#E64B35", Down_regulated = "#00468B"),
        labels = c(
          Up_regulated = "Up-regulated",
          Down_regulated = "Down-regulated"
        )
      ) +
      ggplot2::labs(
        title = "Up- and Down-regulated Genes by Comparison",
        x = "Comparison",
        y = "Number of Genes",
        fill = "Regulation"
      )

    ggplot2::ggsave(
      file.path(sum_dir, "up_down_genes_by_comparison.pdf"),
      p2,
      width = 15,
      height = 10,
      device = "pdf"
    )

    # Heatmap
    heatmap_data <- summary_df %>%
      tidyr::pivot_wider(
        id_cols = "CellType",
        names_from = "Comparison",
        values_from = "Total_DEGs",
        values_fill = 0
      )

    hm_matrix <- as.matrix(heatmap_data[, -1])
    rownames(hm_matrix) <- heatmap_data$CellType

    if (nrow(hm_matrix) > 1 && ncol(hm_matrix) > 1) {
      pheatmap::pheatmap(
        hm_matrix,
        display_numbers = TRUE,
        fontsize_number = 8,
        main = "DEG Counts Heatmap: Cell Type × Comparison",
        filename = file.path(sum_dir, "comparison_heatmap.pdf"),
        width = 12,
        height = 10,
        color = colorRampPalette(c("white", "#E64B35"))(100)
      )
    }
  } else {
    message(
      "⚠️ No comparisons completed - check QC thresholds and tissue/sample distribution"
    )
  }

  # ===== 12. GSVA Analysis (OFFLINE MODE) =====
  if (run_gsva) {
    message("\n==== GSVA Pathway Activity Analysis (Offline Mode) ====")
    message("NOTE: Using local GMT file (no internet required)")

    # GSVA version compatibility wrapper
    run_gsva_safe <- function(expr_mat, gene_sets) {
      if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
        message("Using GSVA new API (gsvaParam)")
        param <- GSVA::gsvaParam(expr_mat, gene_sets, kcdf = "Gaussian")
        GSVA::gsva(param)
      } else {
        message("Using GSVA legacy API")
        GSVA::gsva(
          expr_mat,
          gene_sets,
          method = "gsva",
          kcdf = "Gaussian",
          mx.diff = TRUE,
          verbose = FALSE,
          parallel.sz = if (n_cores > 1) n_cores else 1
        )
      }
    }

    # Convert to log-CPM
    counts_to_logcpm <- function(cnt_mat) {
      lib_sizes <- colSums(cnt_mat)
      lib_sizes[lib_sizes == 0] <- 1
      cpm_mat <- t(t(cnt_mat) / lib_sizes * 1e6)
      log2(cpm_mat + 1)
    }

    pb_logcpm <- counts_to_logcpm(as.matrix(pb_counts))

    # Load gene sets from local GMT
    message("Loading gene sets from local GMT file...")

    geneset_cache_dir <- file.path(gsva_dir, "genesets_cache")
    dir.create(geneset_cache_dir, showWarnings = FALSE, recursive = TRUE)

    hallmark_rds <- file.path(geneset_cache_dir, "msigdb_HALLMARK.rds")
    keggmed_rds <- file.path(geneset_cache_dir, "msigdb_KEGG_MEDICUS.rds")
    keggleg_rds <- file.path(geneset_cache_dir, "msigdb_KEGG_LEGACY.rds")

    gene_universe <- rownames(pb_logcpm)

    if (
      file.exists(hallmark_rds) &&
        file.exists(keggmed_rds) &&
        file.exists(keggleg_rds)
    ) {
      message("Loading cached gene sets from RDS...")
      hallmark_sets <- readRDS(hallmark_rds)
      keggmed_sets <- readRDS(keggmed_rds)
      keggleg_sets <- readRDS(keggleg_rds)
    } else {
      message("Parsing GMT file (first time - will be cached)...")
      gs_sub <- load_msigdb_subsets_from_gmt(
        gmt_file = msigdb_gmt,
        gene_universe = gene_universe,
        min_size = 10,
        max_size = 5000
      )

      hallmark_sets <- gs_sub$hallmark
      keggmed_sets <- gs_sub$kegg_medicus
      keggleg_sets <- gs_sub$kegg_legacy

      saveRDS(hallmark_sets, hallmark_rds)
      saveRDS(keggmed_sets, keggmed_rds)
      saveRDS(keggleg_sets, keggleg_rds)
      message("Gene sets cached for future runs")
    }

    # Select KEGG sets (prefer KEGG_MEDICUS)
    if (length(keggmed_sets) >= 5) {
      kegg_sets <- keggmed_sets
      message(sprintf("Using KEGG_MEDICUS gene sets: %d", length(kegg_sets)))
    } else if (length(keggleg_sets) >= 5) {
      kegg_sets <- keggleg_sets
      message(sprintf("Using legacy KEGG gene sets: %d", length(kegg_sets)))
    } else {
      stop("No usable KEGG gene sets found after filtering")
    }

    message(sprintf(
      "Final gene sets: Hallmark=%d, KEGG=%d",
      length(hallmark_sets),
      length(kegg_sets)
    ))

    # FIX: Use pseudobulk groups as input (not celltype means)
    # This handles single cell type datasets robustly
    message("\n==== Running GSVA on Pseudobulk Groups ====")
    message(sprintf(
      "Input matrix: %d genes × %d pseudobulk groups",
      nrow(pb_logcpm),
      ncol(pb_logcpm)
    ))

    # Column annotation for heatmaps
    gsva_anno_col <- data.frame(
      Tissue = pb_metadata$tissue,
      CellType = pb_metadata$celltype,
      stringsAsFactors = FALSE
    )
    rownames(gsva_anno_col) <- colnames(pb_logcpm)

    # GSVA helper function
    run_gsva_and_plot <- function(
      expr_mat,
      gene_sets,
      prefix,
      top_var,
      anno_col = NULL
    ) {
      common_genes <- intersect(rownames(expr_mat), unique(unlist(gene_sets)))

      if (length(common_genes) < 50) {
        message(sprintf(
          "[SKIP] %s: Too few common genes (%d)",
          prefix,
          length(common_genes)
        ))
        return(NULL)
      }

      expr_subset <- expr_mat[common_genes, , drop = FALSE]
      sets_subset <- lapply(gene_sets, function(g) intersect(g, common_genes))
      sets_subset <- sets_subset[vapply(sets_subset, length, integer(1)) >= 10]

      if (length(sets_subset) < 5) {
        message(sprintf(
          "[SKIP] %s: Too few gene sets (%d)",
          prefix,
          length(sets_subset)
        ))
        return(NULL)
      }

      message(sprintf(
        "Running GSVA for %s: %d genes, %d pathways, %d samples",
        prefix,
        nrow(expr_subset),
        length(sets_subset),
        ncol(expr_subset)
      ))

      gsva_result <- tryCatch(
        {
          run_gsva_safe(expr_subset, sets_subset)
        },
        error = function(e) {
          message(sprintf("[ERROR] GSVA failed for %s: %s", prefix, e$message))
          return(NULL)
        }
      )

      if (is.null(gsva_result) || nrow(gsva_result) < 2) {
        message(sprintf(
          "[SKIP] %s: GSVA returned insufficient pathways",
          prefix
        ))
        return(NULL)
      }

      message(sprintf(
        "GSVA complete: %d pathways × %d samples",
        nrow(gsva_result),
        ncol(gsva_result)
      ))

      # Save results
      saveRDS(
        gsva_result,
        file.path(gsva_dir, paste0("gsva_", prefix, "_scores.rds"))
      )
      utils::write.csv(
        gsva_result,
        file.path(gsva_dir, paste0("gsva_", prefix, "_scores.csv"))
      )

      # Select top variable pathways
      row_vars <- apply(gsva_result, 1, var)
      row_vars <- row_vars[is.finite(row_vars)]

      if (length(row_vars) == 0) {
        message(sprintf("[SKIP] %s: No finite variances", prefix))
        return(gsva_result)
      }

      n_select <- min(top_var, length(row_vars))
      top_pathways <- names(tail(sort(row_vars), n_select))

      message(sprintf(
        "Creating heatmap with top %d variable pathways",
        n_select
      ))

      # Create heatmap
      hm_args <- list(
        mat = gsva_result[top_pathways, , drop = FALSE],
        scale = "row",
        show_colnames = FALSE,
        main = sprintf(
          "%s GSVA Scores (Top %d Variable Pathways)",
          gsub("_", " ", prefix),
          n_select
        ),
        filename = file.path(gsva_dir, paste0(prefix, "_gsva_heatmap.pdf")),
        width = 12,
        height = max(6, 0.12 * n_select + 4),
        color = colorRampPalette(c("#00468B", "white", "#E64B35"))(100)
      )

      if (!is.null(anno_col)) {
        hm_args$annotation_col <- anno_col
      }

      tryCatch(
        {
          do.call(pheatmap::pheatmap, hm_args)
          message(sprintf("Heatmap saved: %s_gsva_heatmap.pdf", prefix))
        },
        error = function(e) {
          message(sprintf("Heatmap creation failed: %s", e$message))
        }
      )

      return(gsva_result)
    }

    # Run global GSVA
    gsva_hallmark_global <- run_gsva_and_plot(
      pb_logcpm,
      hallmark_sets,
      "hallmark_global",
      gsva_top_var_h,
      gsva_anno_col
    )

    gsva_kegg_global <- run_gsva_and_plot(
      pb_logcpm,
      kegg_sets,
      "kegg_global",
      gsva_top_var_k,
      gsva_anno_col
    )

    # Per-tissue GSVA (tissue-specific pseudobulk groups)
    message("\n==== Per-Tissue GSVA Analysis ====")
    tissues_for_gsva <- sort(unique(pb_metadata$tissue))
    tissue_gsva_results <- list()

    for (tissue in tissues_for_gsva) {
      tissue_indices <- which(pb_metadata$tissue == tissue)

      if (length(tissue_indices) < 3) {
        message(sprintf(
          "[SKIP] %s: Too few pseudobulk groups (%d)",
          tissue,
          length(tissue_indices)
        ))
        next
      }

      tissue_expr <- pb_logcpm[, tissue_indices, drop = FALSE]
      tissue_anno <- gsva_anno_col[tissue_indices, , drop = FALSE]

      message(sprintf(
        "Processing %s: %d pseudobulk groups",
        tissue,
        ncol(tissue_expr)
      ))

      tissue_prefix_h <- paste0("hallmark_", sanitize_name(tissue))
      tissue_prefix_k <- paste0("kegg_", sanitize_name(tissue))

      tissue_gsva_h <- run_gsva_and_plot(
        tissue_expr,
        hallmark_sets,
        tissue_prefix_h,
        gsva_top_var_h,
        tissue_anno
      )

      tissue_gsva_k <- run_gsva_and_plot(
        tissue_expr,
        kegg_sets,
        tissue_prefix_k,
        gsva_top_var_k,
        tissue_anno
      )

      tissue_gsva_results[[tissue]] <- list(
        hallmark = tissue_gsva_h,
        kegg = tissue_gsva_k
      )
    }

    # Tissue-wise pathway differences using limma (statistical testing)
    if (length(tissue_gsva_results) >= 2) {
      message("\n==== Tissue-wise Pathway Activity Differences (limma) ====")
      message("Statistical testing of pathway differences between tissues")

      # Check if limma is available
      if (!requireNamespace("limma", quietly = TRUE)) {
        message(
          "⚠️ limma not installed, skipping statistical pathway comparison"
        )
        message("   Install with: BiocManager::install('limma')")
      } else {
        # Run limma on global GSVA scores (pseudobulk-level)
        run_limma_pathway_comparison <- function(
          gsva_mat,
          pb_meta,
          pathway_type
        ) {
          # Ensure no NA tissues
          valid_idx <- !is.na(pb_meta$tissue) & pb_meta$tissue != ""
          gsva_mat <- gsva_mat[, valid_idx, drop = FALSE]
          pb_meta <- pb_meta[valid_idx, , drop = FALSE]

          # Sanitize tissue names for formula
          pb_meta$tissue_safe <- factor(make.names(pb_meta$tissue))

          # Design matrix (no intercept for easier contrasts)
          design <- stats::model.matrix(~ 0 + tissue_safe, data = pb_meta)
          colnames(design) <- sub("^tissue_safe", "", colnames(design))

          # Fit linear model
          fit <- limma::lmFit(gsva_mat, design)

          # All pairwise contrasts
          tissues_safe <- levels(pb_meta$tissue_safe)
          if (length(tissues_safe) < 2) {
            return(invisible(NULL))
          }

          tissue_pairs <- combn(tissues_safe, 2, simplify = FALSE)

          for (pair in tissue_pairs) {
            t1 <- pair[1]
            t2 <- pair[2]

            # Original tissue names (for output)
            t1_orig <- unique(pb_meta$tissue[make.names(pb_meta$tissue) == t1])[
              1
            ]
            t2_orig <- unique(pb_meta$tissue[make.names(pb_meta$tissue) == t2])[
              1
            ]

            contrast_str <- paste0(t2, " - ", t1)

            tryCatch(
              {
                cont_mat <- limma::makeContrasts(
                  contrasts = contrast_str,
                  levels = design
                )
                fit2 <- limma::contrasts.fit(fit, cont_mat)
                fit2 <- limma::eBayes(fit2)

                # Get results
                top_table <- limma::topTable(fit2, number = Inf, sort.by = "P")

                # Save results
                out_prefix <- paste0(
                  "limma_",
                  pathway_type,
                  "_",
                  sanitize_name(t2_orig),
                  "_vs_",
                  sanitize_name(t1_orig)
                )
                utils::write.csv(
                  top_table,
                  file.path(gsva_dir, paste0(out_prefix, ".csv")),
                  row.names = TRUE
                )

                # Count significant pathways
                n_sig <- sum(top_table$adj.P.Val < 0.05, na.rm = TRUE)
                message(sprintf(
                  "  %s vs %s (%s): %d/%d pathways FDR<0.05",
                  t2_orig,
                  t1_orig,
                  pathway_type,
                  n_sig,
                  nrow(top_table)
                ))

                # Create volcano plot for pathways
                if (nrow(top_table) > 0) {
                  top_table$pathway <- rownames(top_table)
                  top_table$neg_log10_padj <- pmax(
                    -log10(top_table$adj.P.Val),
                    0
                  )
                  top_table$neg_log10_padj[
                    !is.finite(top_table$neg_log10_padj)
                  ] <-
                    max(
                      top_table$neg_log10_padj[is.finite(
                        top_table$neg_log10_padj
                      )],
                      na.rm = TRUE
                    )

                  p_volcano <- ggplot2::ggplot(
                    top_table,
                    ggplot2::aes(x = logFC, y = neg_log10_padj)
                  ) +
                    ggplot2::geom_point(
                      ggplot2::aes(color = adj.P.Val < 0.05),
                      size = 2,
                      alpha = 0.6
                    ) +
                    ggplot2::scale_color_manual(
                      values = c("TRUE" = "#E64B35", "FALSE" = "gray70"),
                      labels = c("TRUE" = "FDR<0.05", "FALSE" = "NS"),
                      name = ""
                    ) +
                    ggplot2::geom_hline(
                      yintercept = -log10(0.05),
                      linetype = "dashed",
                      color = "black",
                      linewidth = 0.4
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                      panel.grid = ggplot2::element_blank(),
                      plot.title = ggplot2::element_text(
                        hjust = 0.5,
                        face = "bold"
                      )
                    ) +
                    ggplot2::labs(
                      title = sprintf(
                        "Pathway Activity: %s vs %s (%s)",
                        t2_orig,
                        t1_orig,
                        toupper(pathway_type)
                      ),
                      x = "Log2(Fold Change)",
                      y = "-log10(Adjusted p-value)"
                    )

                  ggplot2::ggsave(
                    file.path(gsva_dir, paste0(out_prefix, "_volcano.pdf")),
                    plot = p_volcano,
                    width = 10,
                    height = 8,
                    device = "pdf"
                  )
                }

                # Heatmap of top significant pathways
                if (n_sig > 0) {
                  top_sig <- rownames(top_table)[top_table$adj.P.Val < 0.05]
                  n_show <- min(25, length(top_sig))

                  if (n_show > 0) {
                    # Get data for this tissue pair only
                    pair_idx <- pb_meta$tissue_safe %in% c(t1, t2)
                    pair_gsva <- gsva_mat[
                      top_sig[1:n_show],
                      pair_idx,
                      drop = FALSE
                    ]
                    pair_anno <- data.frame(
                      Tissue = pb_meta$tissue[pair_idx],
                      stringsAsFactors = FALSE
                    )
                    rownames(pair_anno) <- colnames(pair_gsva)

                    pheatmap::pheatmap(
                      pair_gsva,
                      annotation_col = pair_anno,
                      scale = "row",
                      show_colnames = FALSE,
                      main = sprintf(
                        "Top %d Differential Pathways (%s): %s vs %s\n(FDR<0.05)",
                        n_show,
                        toupper(pathway_type),
                        t2_orig,
                        t1_orig
                      ),
                      filename = file.path(
                        gsva_dir,
                        paste0(out_prefix, "_heatmap.pdf")
                      ),
                      width = 12,
                      height = max(6, 0.15 * n_show + 4),
                      color = colorRampPalette(c(
                        "#00468B",
                        "white",
                        "#E64B35"
                      ))(100)
                    )
                  }
                }
              },
              error = function(e) {
                message(sprintf(
                  "  [ERROR] limma failed for %s vs %s: %s",
                  t2_orig,
                  t1_orig,
                  e$message
                ))
              }
            )
          }

          invisible(NULL)
        }

        # Run limma on Hallmark
        if (!is.null(gsva_hallmark_global) && nrow(gsva_hallmark_global) > 0) {
          message("\nHallmark pathway comparisons:")
          run_limma_pathway_comparison(
            gsva_hallmark_global,
            pb_metadata,
            "hallmark"
          )
        }

        # Run limma on KEGG
        if (!is.null(gsva_kegg_global) && nrow(gsva_kegg_global) > 0) {
          message("\nKEGG pathway comparisons:")
          run_limma_pathway_comparison(gsva_kegg_global, pb_metadata, "kegg")
        }
      }
    }

    message("GSVA analysis complete")
  }

  # ===== 13. Save Session Info =====
  message("\n==== Saving Session Info ====")
  session_info <- capture.output(sessionInfo())
  writeLines(session_info, file.path(out_root, "session_info.txt"))

  params <- list(
    cell_anno_col = cell_anno_col,
    tissue_col = tissue_col,
    sample_col = sample_col,
    min_cell_per_sample = min_cell_per_sample,
    min_sample_per_tissue = min_sample_per_tissue,
    min_pb_libsize = min_pb_libsize,
    min_pb_detected_genes = min_pb_detected_genes,
    padj_thr = padj_thr,
    lfc_thr = lfc_thr,
    min_gene_total_counts = min_gene_total_counts,
    run_gsva = run_gsva,
    run_go = run_go,
    msigdb_gmt = if (run_gsva) msigdb_gmt else NA,
    use_lfc_shrink = use_lfc_shrink,
    n_cores = n_cores,
    analysis_date = Sys.Date(),
    R_version = R.version.string,
    Seurat_version = as.character(packageVersion("Seurat")),
    DESeq2_version = as.character(packageVersion("DESeq2")),
    GSVA_version = if (run_gsva) as.character(packageVersion("GSVA")) else NA
  )

  saveRDS(params, file.path(out_root, "analysis_parameters.rds"))

  # ===== 14. Final Summary =====
  message("\n==== Analysis Complete ====")
  message(sprintf("All results saved to: %s", out_root))
  message("\nOutput structure:")
  message(sprintf("  - Pseudobulk data: %s", pb_dir))
  message(sprintf("  - QC metrics: %s", qc_dir))
  message(sprintf("  - DESeq2 results: %s", de_dir))
  message(sprintf("  - Summary visualizations: %s", sum_dir))
  if (run_gsva) {
    message(sprintf("  - GSVA results: %s", gsva_dir))
  }
  message(sprintf("  - Session info: %s/session_info.txt", out_root))

  invisible(TRUE)
}

# =============================================================================
# One-vs-Rest MAST Comparison v1.0
# =============================================================================
# Version: 1.0 (2024-12-23)
#
# Purpose:
# Compare each group vs all other groups combined
# (e.g., cluster_0 vs rest, cluster_1 vs rest, ...)
#
# Use Cases:
# - Find cluster-specific markers
# - Identify cell type signature genes
# - Discover treatment-specific responses
# - Detect condition-specific patterns
#
# Method:
# For each group:
#   group_X vs (all other groups combined as "rest")
#
# Based on: v4.2 MAST optimization (sparse matrices, ref/case encoding)
#
# Author: Clinical Bioinformatics Team
# =============================================================================

#' One-vs-Rest MAST Comparison
#'
#' @description
#' For each group in group_col, compare it against all other groups combined.
#' This is useful for finding group-specific marker genes.
#'
#' Example: For clusters 0, 1, 2, 3:
#'   - cluster_0 vs rest (rest = clusters 1+2+3 combined)
#'   - cluster_1 vs rest (rest = clusters 0+2+3 combined)
#'   - cluster_2 vs rest (rest = clusters 0+1+3 combined)
#'   - cluster_3 vs rest (rest = clusters 0+1+2 combined)
#'
#' @param seurat_obj Seurat object with normalized data
#' @param group_col Column defining groups (e.g., "seurat_clusters", "celltype")
#' @param batch_col Optional batch column (default: NULL)
#' @param groups_to_test Optional vector of specific groups to test.
#'                       If NULL (default), tests all groups.
#' @param assay Assay to use (default: "RNA")
#' @param slot_use Slot to use (default: "data")
#' @param min_cells_in_group Min cells in focal group (default: 50)
#' @param min_cells_in_rest Min total cells in "rest" (default: 100)
#' @param max_cells_downsample Max cells to use (downsamples if more, NULL=no limit)
#' @param min_expr_cells_per_gene Min cells expressing gene (default: 10)
#' @param padj_thr Adjusted p-value threshold (default: 0.05)
#' @param lfc_thr Log fold change threshold in log2 scale (default: 0.25)
#' @param scale_covariates Scale covariates (default: TRUE)
#' @param convert_to_log2 Convert to log2 (default: TRUE)
#' @param control_mt Include percent.mt (default: FALSE)
#' @param output_dir Output directory (default: "./mast_one_vs_rest")
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' \dontrun{
#' # Example 1: Find markers for all clusters
#' run_one_vs_rest_mast(
#'   seurat_obj = pbmc,
#'   group_col = "seurat_clusters"
#' )
#' # Output: cluster_0_vs_rest, cluster_1_vs_rest, ...
#'
#' # Example 2: Find markers for specific clusters only
#' run_one_vs_rest_mast(
#'   seurat_obj = pbmc,
#'   group_col = "seurat_clusters",
#'   groups_to_test = c("0", "1", "3")  # Only these clusters
#' )
#'
#' # Example 3: With batch correction
#' run_one_vs_rest_mast(
#'   seurat_obj = lung,
#'   group_col = "celltype",
#'   batch_col = "sample"
#' )
#'
#' # Example 4: Downsample large datasets
#' run_one_vs_rest_mast(
#'   seurat_obj = large_obj,
#'   group_col = "seurat_clusters",
#'   max_cells_downsample = 50000  # Use max 50k cells total
#' )
#' }
#'
#' @export
run_one_vs_rest_mast <- function(
  seurat_obj,
  group_col,
  batch_col = NULL,
  groups_to_test = NULL,
  assay = "RNA",
  slot_use = "data",
  min_cells_in_group = 50,
  min_cells_in_rest = 100,
  max_cells_downsample = NULL,
  min_expr_cells_per_gene = 10,
  padj_thr = 0.05,
  lfc_thr = 0.25,
  scale_covariates = TRUE,
  convert_to_log2 = TRUE,
  control_mt = FALSE,
  output_dir = "./mast_one_vs_rest"
) {
  # ===== 1. Dependency Check =====
  required_pkgs <- c("MAST", "Seurat", "data.table", "Matrix")
  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {
    stop("Please install: ", paste(missing_pkgs, collapse = ", "))
  }

  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  # ===== 2. Setup =====
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  sum_dir <- file.path(output_dir, "summary")
  qc_dir <- file.path(output_dir, "qc_metrics")
  dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  message("==== One-vs-Rest MAST Comparison v1.0 ====")
  message(sprintf("Date: %s", Sys.Date()))
  message(sprintf("MAST version: %s", packageVersion("MAST")))

  # ===== 3. Metadata Validation =====
  meta <- seurat_obj@meta.data

  if (!group_col %in% colnames(meta)) {
    stop("group_col '", group_col, "' not found")
  }

  use_batch_correction <- FALSE
  if (!is.null(batch_col) && nzchar(batch_col)) {
    if (!batch_col %in% colnames(meta)) {
      message("⚠️ batch_col '", batch_col, "' not found, skipping")
      batch_col <- NULL
    } else {
      use_batch_correction <- TRUE
      message(sprintf("✓ Batch correction: '%s'", batch_col))
    }
  }

  # Clean metadata
  message("\n==== Cleaning Metadata ====")
  meta[[group_col]] <- trimws(as.character(meta[[group_col]]))

  bad_idx <- is.na(meta[[group_col]]) | meta[[group_col]] == ""

  if (use_batch_correction) {
    meta[[batch_col]] <- trimws(as.character(meta[[batch_col]]))
    bad_idx <- bad_idx | is.na(meta[[batch_col]]) | meta[[batch_col]] == ""
  }

  if (sum(bad_idx) > 0) {
    message(sprintf("Dropping %d cells with NA/empty", sum(bad_idx)))
    seurat_obj <- Seurat::subset(seurat_obj, cells = rownames(meta)[!bad_idx])
    meta <- seurat_obj@meta.data
  }

  message(sprintf("Total cells: %d", nrow(meta)))
  message(sprintf("Total genes: %d", nrow(seurat_obj)))

  if (control_mt && !"percent.mt" %in% colnames(meta)) {
    message("⚠️ percent.mt not found, skipping")
    control_mt <- FALSE
  }

  # ===== 4. Downsampling (if needed) =====
  if (!is.null(max_cells_downsample) && nrow(meta) > max_cells_downsample) {
    message(sprintf(
      "\n==== Downsampling: %d → %d cells ====",
      nrow(meta),
      max_cells_downsample
    ))

    # Stratified downsampling (proportional to group sizes)
    set.seed(42)
    group_props <- table(meta[[group_col]]) / nrow(meta)
    target_per_group <- round(group_props * max_cells_downsample)
    target_per_group <- pmax(target_per_group, 10) # At least 10 per group

    sampled_cells <- unlist(lapply(names(target_per_group), function(g) {
      cells_in_g <- rownames(meta)[meta[[group_col]] == g]
      if (length(cells_in_g) <= target_per_group[g]) {
        return(cells_in_g)
      } else {
        return(sample(cells_in_g, target_per_group[g]))
      }
    }))

    seurat_obj <- Seurat::subset(seurat_obj, cells = sampled_cells)
    meta <- seurat_obj@meta.data

    message(sprintf("After downsampling: %d cells", nrow(meta)))
  }

  # ===== 5. Get Expression Data (SPARSE) =====
  message("\n==== Loading Expression Data ====")

  expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot_use)
  if (!inherits(expr_data, "dgCMatrix")) {
    expr_data <- as(expr_data, "dgCMatrix")
  }

  message(sprintf(
    "Expression: %d genes × %d cells (sparse)",
    nrow(expr_data),
    ncol(expr_data)
  ))

  expr_counts <- Seurat::GetAssayData(
    seurat_obj,
    assay = assay,
    slot = "counts"
  )

  # ===== 6. Calculate Covariates =====
  message("\n==== Calculating Covariates ====")

  if ("nFeature_RNA" %in% colnames(meta)) {
    meta$cngeneson <- meta$nFeature_RNA
  } else {
    meta$cngeneson <- Matrix::colSums(expr_counts > 0)
  }

  if ("nCount_RNA" %in% colnames(meta)) {
    meta$log10_umi <- log10(meta$nCount_RNA + 1)
  } else {
    meta$log10_umi <- log10(Matrix::colSums(expr_counts) + 1)
  }

  # Save QC
  qc_df <- data.frame(
    cell_id = rownames(meta),
    group = meta[[group_col]],
    cngeneson = meta$cngeneson,
    log10_umi = meta$log10_umi,
    stringsAsFactors = FALSE
  )

  if (use_batch_correction) {
    qc_df$batch <- meta[[batch_col]]
  }
  if (control_mt) {
    qc_df$percent.mt <- meta$percent.mt
  }

  utils::write.csv(
    qc_df,
    file.path(qc_dir, "cell_qc_metrics.csv"),
    row.names = FALSE
  )

  # QC plots
  pdf(file.path(qc_dir, "covariate_distributions.pdf"), width = 12, height = 8)
  par(mfrow = c(2, 2))

  hist(
    meta$cngeneson,
    breaks = 50,
    main = "Genes Detected",
    xlab = "cngeneson",
    col = "steelblue",
    border = "white"
  )

  hist(
    meta$log10_umi,
    breaks = 50,
    main = "Library Size",
    xlab = "log10(UMI + 1)",
    col = "steelblue",
    border = "white"
  )

  boxplot(
    meta$cngeneson ~ meta[[group_col]],
    las = 2,
    main = paste("Detection by", group_col),
    ylab = "cngeneson",
    col = "lightblue",
    cex.axis = 0.7
  )

  boxplot(
    meta$log10_umi ~ meta[[group_col]],
    las = 2,
    main = paste("Library by", group_col),
    ylab = "log10_umi",
    col = "lightblue",
    cex.axis = 0.7
  )

  dev.off()

  # ===== 7. Determine Groups =====
  message("\n==== Determining Groups ====")

  all_groups <- sort(unique(meta[[group_col]]))
  all_groups <- all_groups[!is.na(all_groups) & all_groups != ""]

  if (!is.null(groups_to_test)) {
    missing <- setdiff(groups_to_test, all_groups)
    if (length(missing) > 0) {
      stop(
        "groups_to_test contains non-existent: ",
        paste(missing, collapse = ", ")
      )
    }
    groups <- groups_to_test
    message(sprintf("User-specified groups: %d", length(groups)))
  } else {
    groups <- all_groups
    message(sprintf("All groups: %d", length(groups)))
  }

  # Check cell counts
  group_counts <- table(meta[[group_col]])
  valid_groups <- names(group_counts[group_counts >= min_cells_in_group])
  valid_groups <- intersect(groups, valid_groups)

  if (length(valid_groups) == 0) {
    stop("No groups with >= ", min_cells_in_group, " cells")
  }

  message(sprintf(
    "Valid groups (>= %d cells): %d",
    min_cells_in_group,
    length(valid_groups)
  ))
  message("Group cell counts:")
  print(group_counts[valid_groups])

  # ===== 8. Helper Functions =====
  sanitize_name <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
  }

  # ===== 9. One-vs-Rest Comparisons =====
  message("\n==== Performing One-vs-Rest Comparisons ====")
  message(sprintf("Total comparisons: %d", length(valid_groups)))

  all_results <- list()

  for (focal_group in valid_groups) {
    message(sprintf("\n--- Testing: %s vs rest ---", focal_group))

    # Get cells
    focal_cells <- rownames(meta)[meta[[group_col]] == focal_group]
    rest_cells <- rownames(meta)[meta[[group_col]] != focal_group]

    n_focal <- length(focal_cells)
    n_rest <- length(rest_cells)

    message(sprintf("Cells: %s=%d, rest=%d", focal_group, n_focal, n_rest))

    if (n_focal < min_cells_in_group) {
      message("[SKIP] Focal group too small")
      next
    }

    if (n_rest < min_cells_in_rest) {
      message("[SKIP] Rest group too small")
      next
    }

    tryCatch(
      {
        # Combine focal + rest
        all_cells <- c(focal_cells, rest_cells)
        subset_meta <- meta[all_cells, , drop = FALSE]

        # Create binary group variable
        subset_meta$group_binary <- ifelse(
          subset_meta[[group_col]] == focal_group,
          focal_group,
          "rest"
        )

        # ===== 10. Extract & Filter Expression =====
        subset_expr <- expr_data[, all_cells, drop = FALSE]

        # Gene filtering
        gene_expr_count <- Matrix::rowSums(subset_expr > 0)
        keep_genes <- gene_expr_count >= min_expr_cells_per_gene

        message(sprintf("Genes: %d → %d", nrow(subset_expr), sum(keep_genes)))

        subset_expr <- subset_expr[keep_genes, , drop = FALSE]

        if (nrow(subset_expr) < 100) {
          message("[SKIP] Too few genes")
          next
        }

        # Scale conversion
        if (convert_to_log2) {
          subset_expr <- subset_expr / log(2)
          lfc_thr_use <- lfc_thr
        } else {
          lfc_thr_use <- lfc_thr * log(2)
        }

        # Densify
        subset_expr <- as.matrix(subset_expr)

        # ===== 11. Prepare Metadata =====
        subset_meta$cngeneson <- meta[all_cells, "cngeneson"]
        subset_meta$log10_umi <- meta[all_cells, "log10_umi"]

        # Scale
        if (scale_covariates) {
          subset_meta$cngeneson_s <- as.numeric(scale(subset_meta$cngeneson))
          subset_meta$log10_umi_s <- as.numeric(scale(subset_meta$log10_umi))
          cov_suffix <- "_s"
        } else {
          subset_meta$cngeneson_s <- subset_meta$cngeneson
          subset_meta$log10_umi_s <- subset_meta$log10_umi
          cov_suffix <- "_s"
        }

        # ref/case encoding (rest=ref, focal=case)
        subset_meta$group_safe <- factor(
          ifelse(subset_meta$group_binary == "rest", "ref", "case"),
          levels = c("ref", "case")
        )

        if (use_batch_correction) {
          subset_meta$batch_safe <- factor(make.names(subset_meta[[batch_col]]))
        }

        if (control_mt) {
          subset_meta$percent.mt <- meta[all_cells, "percent.mt"]
          if (scale_covariates) {
            subset_meta$percent.mt_s <- as.numeric(scale(
              subset_meta$percent.mt
            ))
          } else {
            subset_meta$percent.mt_s <- subset_meta$percent.mt
          }
        }

        # ===== 12. Build SCA =====
        sca <- MAST::FromMatrix(
          exprsArray = subset_expr,
          cData = subset_meta,
          fData = data.frame(
            primerid = rownames(subset_expr),
            stringsAsFactors = FALSE
          )
        )

        # ===== 13. Build Formula =====
        formula_parts <- c(
          "~ group_safe",
          paste0("cngeneson", cov_suffix),
          paste0("log10_umi", cov_suffix)
        )

        if (control_mt) {
          formula_parts <- c(formula_parts, paste0("percent.mt", cov_suffix))
        }

        if (use_batch_correction) {
          formula_parts <- c(formula_parts, "batch_safe")
        }

        formula_str <- paste(formula_parts, collapse = " + ")
        form <- stats::as.formula(formula_str)

        message(sprintf("Formula: %s", formula_str))

        fit <- MAST::zlm(form, sca, method = "glm", ebayes = FALSE)

        # ===== 14. Extract Results =====
        coef_name <- "group_safecase"

        summ <- summary(fit, doLRT = coef_name)
        dt <- data.table::as.data.table(summ$datatable)

        if (!coef_name %in% unique(dt$contrast)) {
          message("[ERROR] Coefficient not found")
          next
        }

        p_dt <- dt[
          contrast == coef_name & component == "H",
          .(primerid, pvalue = `Pr(>Chisq)`)
        ]

        fc_dt <- dt[
          contrast == coef_name & component == "logFC",
          .(primerid, logFC = coef)
        ]

        res <- merge(fc_dt, p_dt, by = "primerid", all = TRUE)

        res[, padj := p.adjust(pvalue, method = "BH")]

        res[,
          regulation := data.table::fcase(
            !is.na(padj) & padj <= padj_thr & logFC >= lfc_thr_use  , "up"   ,
            !is.na(padj) & padj <= padj_thr & logFC <= -lfc_thr_use , "down" ,
            default = "stable"
          )
        ]

        # Sort by logFC (descending) for markers
        data.table::setorder(res, -logFC)

        n_up <- sum(res$regulation == "up", na.rm = TRUE)
        n_down <- sum(res$regulation == "down", na.rm = TRUE)

        message(sprintf(
          "Markers: %d up (enriched in %s), %d down",
          n_up,
          focal_group,
          n_down
        ))

        # ===== 15. Save Results =====
        comp_name <- paste0(focal_group, "_vs_rest")
        out_subdir <- file.path(output_dir, sanitize_name(comp_name))
        dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)

        data.table::fwrite(res, file.path(out_subdir, "MAST_markers.csv"))

        # Save top markers separately (easier to view)
        top_markers <- res[regulation == "up"][1:min(100, .N)]
        data.table::fwrite(
          top_markers,
          file.path(out_subdir, "top_markers.csv")
        )

        stat_row <- data.table::data.table(
          Focal_group = focal_group,
          Comparison = comp_name,
          N_markers_up = n_up,
          N_markers_down = n_down,
          N_genes_tested = nrow(res),
          N_cells_focal = n_focal,
          N_cells_rest = n_rest,
          Formula = formula_str,
          LogFC_scale = ifelse(convert_to_log2, "log2", "ln"),
          stringsAsFactors = FALSE
        )

        data.table::fwrite(stat_row, file.path(out_subdir, "summary_stats.csv"))

        all_results[[length(all_results) + 1]] <- stat_row

        message(sprintf("✓ Saved: %s/", basename(out_subdir)))
      },
      error = function(e) {
        message(sprintf("[ERROR] %s", e$message))
      }
    )
  }

  # ===== 16. Global Summary =====
  message("\n==== Creating Summary ====")

  if (length(all_results) > 0) {
    summary_df <- data.table::rbindlist(all_results)

    data.table::fwrite(summary_df, file.path(sum_dir, "all_groups_summary.csv"))
    saveRDS(summary_df, file.path(sum_dir, "all_groups_summary.rds"))

    message(sprintf("Completed: %d groups", nrow(summary_df)))

    print(summary_df[, .(Focal_group, N_markers_up, N_cells_focal)])
  } else {
    message("⚠️ No comparisons completed")
  }

  # ===== 17. Save Parameters =====
  session_info <- capture.output(sessionInfo())
  writeLines(session_info, file.path(output_dir, "session_info.txt"))

  params <- list(
    method = "One-vs-Rest MAST",
    version = "1.0",
    group_col = group_col,
    batch_col = batch_col,
    groups_tested = valid_groups,
    min_cells_in_group = min_cells_in_group,
    min_cells_in_rest = min_cells_in_rest,
    max_cells_downsample = max_cells_downsample,
    analysis_date = as.character(Sys.Date())
  )

  saveRDS(params, file.path(output_dir, "analysis_parameters.rds"))

  # ===== 18. Final Message =====
  message("\n==== Analysis Complete ====")
  message(sprintf("Results: %s", output_dir))
  message(sprintf("Groups tested: %d", length(all_results)))
  message(
    "\nFor each group, top markers are in: [group]_vs_rest/top_markers.csv"
  )

  invisible(TRUE)
}


# =============================================================================
# One-vs-Rest MAST with GO/GSEA Enrichment v1.1.1 - HOTFIX
# =============================================================================
# Version: 1.1.1 (2024-12-23)
#
# CRITICAL FIXES from v1.1:
# P0 (will crash):
#   - Complete dependency check (data.table, tidyr, enrichplot)
#   - orgdb <- get() with asNamespace()
#   - enrichplot::barplot() instead of graphics::barplot()
#   - pivot_wider values_fill as list
# P1 (statistical rigor):
#   - GO enrichment with universe (background gene set)
#   - Gene ID type handling (SYMBOL vs ENSEMBL detection)
#   - Local GMT support for offline use
# P2 (performance):
#   - TERM2GENE built once outside loop
#   - sink() with on.exit()
#   - GO clustering uses ID instead of Description
#
# Based on: One-vs-Rest v1.0 + expert review feedback
#
# Author: Clinical Bioinformatics Team
# =============================================================================

#' One-vs-Rest MAST with Enrichment Analysis (PATCHED)
#'
#' @description
#' Extended version with GO/GSEA enrichment analysis and clustering.
#' This version includes critical fixes for production stability.
#'
#' @param seurat_obj Seurat object
#' @param group_col Group column
#' @param batch_col Batch column (optional)
#' @param groups_to_test Groups to test (NULL = all)
#' @param run_go Run GO enrichment (default: TRUE)
#' @param run_gsea Run GSEA (default: TRUE)
#' @param run_clustering Cluster pathways across groups (default: TRUE)
#' @param organism Organism OrgDb (default: "org.Hs.eg.db")
#' @param msigdb_species Species for MSigDB (default: "Homo sapiens")
#' @param msigdb_gmt_file Local GMT file path (PRIORITY, for offline use)
#' @param go_ont GO ontology (BP/MF/CC/ALL, default: "BP")
#' @param min_pathway_size Min genes per pathway (default: 10)
#' @param max_pathway_size Max genes per pathway (default: 500)
#' @param pathway_pval_cutoff Pathway significance (default: 0.05)
#' @param min_genes_for_enrichment Min DE genes for enrichment (default: 10)
#' @param gene_id_type Gene ID type: "auto"/"SYMBOL"/"ENSEMBL" (default: "auto")
#' @param aggregate_duplicates Aggregate duplicate symbols (default: TRUE)
#' @param ... Additional parameters for run_one_vs_rest_mast
#'
#' @return Invisible TRUE
#'
#' @export
run_one_vs_rest_enrichment <- function(
  seurat_obj,
  group_col,
  batch_col = NULL,
  groups_to_test = NULL,
  run_go = TRUE,
  run_gsea = TRUE,
  run_clustering = TRUE,
  organism = "org.Hs.eg.db",
  msigdb_species = "Homo sapiens",
  msigdb_gmt_file = NULL,
  go_ont = "BP",
  min_pathway_size = 10,
  max_pathway_size = 500,
  pathway_pval_cutoff = 0.05,
  min_genes_for_enrichment = 10,
  gene_id_type = "auto",
  aggregate_duplicates = TRUE,
  output_dir = "./mast_one_vs_rest_enriched",
  ...
) {
  # ===== 0. Check Enrichment Dependencies (FIXED: Complete list) =====
  message("\n==== Checking Enrichment Dependencies ====")
  safe_enrich_barplot <- function(x, showCategory = 20, title = NULL) {
    p <- NULL

    # 1) enrichplot::barplot (if exported)
    if (
      exists("barplot", where = asNamespace("enrichplot"), inherits = FALSE)
    ) {
      p <- enrichplot::barplot(x, showCategory = showCategory)
      # 2) clusterProfiler::barplot (often available even when enrichplot doesn't export)
    } else if (
      exists(
        "barplot",
        where = asNamespace("clusterProfiler"),
        inherits = FALSE
      )
    ) {
      p <- clusterProfiler::barplot(x, showCategory = showCategory)
      # 3) fallback: manual ggplot
    } else {
      df <- tryCatch(as.data.frame(x@result), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      df <- df[order(df$p.adjust), , drop = FALSE]
      df <- df[seq_len(min(showCategory, nrow(df))), , drop = FALSE]
      df$Description <- factor(df$Description, levels = rev(df$Description))

      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = Description, y = -log10(p.adjust))
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = NULL, y = "-log10(FDR)")
    }

    if (!is.null(p) && !is.null(title)) {
      p <- p + ggplot2::ggtitle(title)
    }
    p
  }

  base_pkgs <- c("data.table", "tidyr") # FIXED: Added base packages
  enrichment_pkgs <- base_pkgs

  if (run_go) {
    enrichment_pkgs <- c(
      enrichment_pkgs,
      "clusterProfiler",
      "enrichplot",
      organism
    )
  }

  if (run_gsea) {
    enrichment_pkgs <- c(enrichment_pkgs, "clusterProfiler", "enrichplot")
    if (is.null(msigdb_gmt_file)) {
      enrichment_pkgs <- c(enrichment_pkgs, "msigdbr")
    }
  }

  if (run_clustering) {
    enrichment_pkgs <- c(enrichment_pkgs, "pheatmap", "RColorBrewer")
  }

  missing_pkgs <- enrichment_pkgs[
    !vapply(
      enrichment_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Missing packages for enrichment: ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall with: BiocManager::install(c('",
      paste(missing_pkgs, collapse = "', '"),
      "'))"
    )
  }

  # ===== 1. Run Base One-vs-Rest Analysis =====
  message("\n==== Running Base MAST Analysis ====")

  if (!exists("run_one_vs_rest_mast")) {
    stop(
      "run_one_vs_rest_mast() not found. Please source one_vs_rest_mast_v1.0.R first"
    )
  }

  run_one_vs_rest_mast(
    seurat_obj = seurat_obj,
    group_col = group_col,
    batch_col = batch_col,
    groups_to_test = groups_to_test,
    output_dir = output_dir,
    ...
  )

  # ===== 2. Setup Enrichment Directories =====
  enrich_dir <- file.path(output_dir, "enrichment_analysis")
  go_dir <- file.path(enrich_dir, "go_enrichment")
  gsea_dir <- file.path(enrich_dir, "gsea")
  cluster_dir <- file.path(enrich_dir, "pathway_clustering")

  dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)
  if (run_go) {
    dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
  }
  if (run_gsea) {
    dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
  }
  if (run_clustering) {
    dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # ===== 3. Gene ID Type Detection and Mapping =====
  message("\n==== Gene ID Mapping ====")

  all_features <- rownames(seurat_obj)

  # Auto-detect gene ID type (FIXED: Handle mixed SYMBOL/ENSEMBL)
  if (gene_id_type == "auto") {
    n_ensg <- sum(grepl("^ENSG[0-9]+", all_features))
    n_ensmusg <- sum(grepl("^ENSMUSG[0-9]+", all_features))

    if (n_ensg > length(all_features) * 0.5) {
      gene_id_type <- "ENSEMBL"
      message("Auto-detected: ENSEMBL IDs (Human)")
    } else if (n_ensmusg > length(all_features) * 0.5) {
      gene_id_type <- "ENSEMBL"
      message("Auto-detected: ENSEMBL IDs (Mouse)")
    } else {
      gene_id_type <- "SYMBOL"
      message("Auto-detected: SYMBOL IDs (or mixed)")
    }
  }

  message(sprintf("Gene ID type: %s", gene_id_type))

  # Build gene mapping table (FIXED: Use asNamespace)
  tryCatch(
    {
      orgdb <- get(organism, envir = asNamespace(organism)) # FIXED: asNamespace()

      # Clean ENSEMBL IDs (remove version numbers)
      features_clean <- gsub("\\.[0-9]+$", "", all_features)

      # Map to SYMBOL and ENTREZ
      if (gene_id_type == "ENSEMBL") {
        gene_mapping <- suppressMessages(
          clusterProfiler::bitr(
            features_clean,
            fromType = "ENSEMBL",
            toType = c("SYMBOL", "ENTREZID"),
            OrgDb = orgdb
          )
        )
        gene_mapping$feature_raw <- all_features[match(
          gene_mapping$ENSEMBL,
          features_clean
        )]
      } else {
        # SYMBOL or mixed
        gene_mapping <- suppressMessages(
          clusterProfiler::bitr(
            features_clean,
            fromType = "SYMBOL",
            toType = c("ENTREZID", "ENSEMBL"),
            OrgDb = orgdb
          )
        )
        gene_mapping$feature_raw <- all_features[match(
          gene_mapping$SYMBOL,
          features_clean
        )]

        # Try to map remaining as ENSEMBL
        unmapped <- setdiff(features_clean, gene_mapping$SYMBOL)
        if (length(unmapped) > 0) {
          ensg_map <- suppressMessages(
            clusterProfiler::bitr(
              unmapped,
              fromType = "ENSEMBL",
              toType = c("SYMBOL", "ENTREZID"),
              OrgDb = orgdb
            )
          )
          if (nrow(ensg_map) > 0) {
            ensg_map$feature_raw <- all_features[match(
              ensg_map$ENSEMBL,
              features_clean
            )]
            gene_mapping <- rbind(gene_mapping, ensg_map)
          }
        }
      }

      # Standardize column names
      if (
        !"SYMBOL" %in% colnames(gene_mapping) &&
          "gene" %in% colnames(gene_mapping)
      ) {
        gene_mapping$SYMBOL <- gene_mapping$gene
      }

      message(sprintf(
        "Mapped %d / %d genes to SYMBOL/ENTREZ",
        nrow(gene_mapping),
        length(all_features)
      ))
    },
    error = function(e) {
      message("⚠️ Gene ID mapping failed: ", e$message)
      gene_mapping <- NULL
    }
  )

  # ===== 4. Load MSigDB Gene Sets =====
  if (run_gsea) {
    message("\n==== Loading MSigDB Gene Sets ====")

    hallmark_list <- NULL
    kegg_list <- NULL

    # FIXED: Priority to local GMT file
    if (!is.null(msigdb_gmt_file) && file.exists(msigdb_gmt_file)) {
      message(sprintf("Loading local GMT: %s", basename(msigdb_gmt_file)))

      tryCatch(
        {
          gmt_lines <- readLines(msigdb_gmt_file)

          hallmark_lines <- gmt_lines[grepl("^HALLMARK_", gmt_lines)]
          kegg_lines <- gmt_lines[grepl("^KEGG_", gmt_lines)]

          # Parse GMT format
          parse_gmt_line <- function(line) {
            parts <- strsplit(line, "\t")[[1]]
            pathway_name <- parts[1]
            genes <- parts[-(1:2)] # Skip name and description
            genes <- genes[genes != ""]
            list(name = pathway_name, genes = genes)
          }

          hallmark_parsed <- lapply(hallmark_lines, parse_gmt_line)
          hallmark_list <- setNames(
            lapply(hallmark_parsed, function(x) x$genes),
            sapply(hallmark_parsed, function(x) x$name)
          )

          kegg_parsed <- lapply(kegg_lines, parse_gmt_line)
          kegg_list <- setNames(
            lapply(kegg_parsed, function(x) x$genes),
            sapply(kegg_parsed, function(x) x$name)
          )

          message(sprintf(
            "Loaded %d Hallmark + %d KEGG pathways (local GMT)",
            length(hallmark_list),
            length(kegg_list)
          ))
        },
        error = function(e) {
          message("⚠️ Local GMT loading failed: ", e$message)
        }
      )
    }

    # Fallback to msigdbr if no local GMT
    if (is.null(hallmark_list) || is.null(kegg_list)) {
      message("Loading from msigdbr...")

      tryCatch(
        {
          msigdb_all <- msigdbr::msigdbr(species = msigdb_species)

          hallmark_sets <- msigdb_all[msigdb_all$gs_cat == "H", ]
          hallmark_list <- split(
            hallmark_sets$gene_symbol,
            hallmark_sets$gs_name
          )

          kegg_sets <- msigdb_all[
            msigdb_all$gs_cat == "C2" &
              grepl("^KEGG_", msigdb_all$gs_name),
          ]
          kegg_list <- split(kegg_sets$gene_symbol, kegg_sets$gs_name)

          message(sprintf(
            "Loaded %d Hallmark + %d KEGG pathways (msigdbr)",
            length(hallmark_list),
            length(kegg_list)
          ))
        },
        error = function(e) {
          message("⚠️ msigdbr failed: ", e$message)
          run_gsea <- FALSE
        }
      )
    }

    # Filter by size
    if (!is.null(hallmark_list)) {
      hallmark_list <- hallmark_list[
        sapply(hallmark_list, length) >= min_pathway_size &
          sapply(hallmark_list, length) <= max_pathway_size
      ]
    }
    if (!is.null(kegg_list)) {
      kegg_list <- kegg_list[
        sapply(kegg_list, length) >= min_pathway_size &
          sapply(kegg_list, length) <= max_pathway_size
      ]
    }

    # FIXED: Pre-build TERM2GENE outside loop
    if (!is.null(hallmark_list) && length(hallmark_list) > 0) {
      term2gene_hallmark <- data.frame(
        term = rep(names(hallmark_list), sapply(hallmark_list, length)),
        gene = unlist(hallmark_list, use.names = FALSE),
        stringsAsFactors = FALSE
      )
    } else {
      term2gene_hallmark <- NULL
    }

    if (!is.null(kegg_list) && length(kegg_list) > 0) {
      term2gene_kegg <- data.frame(
        term = rep(names(kegg_list), sapply(kegg_list, length)),
        gene = unlist(kegg_list, use.names = FALSE),
        stringsAsFactors = FALSE
      )
    } else {
      term2gene_kegg <- NULL
    }
  }

  # ===== 5. Find Group Result Directories =====
  group_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  group_dirs <- group_dirs[
    !grepl("summary|qc_metrics|enrichment", basename(group_dirs))
  ]

  if (length(group_dirs) == 0) {
    message("⚠️ No group results found for enrichment")
    return(invisible(TRUE))
  }

  message(sprintf(
    "\n==== Performing Enrichment for %d Groups ====",
    length(group_dirs)
  ))

  # Storage for clustering
  all_go_results <- list()
  all_gsea_hallmark <- list()
  all_gsea_kegg <- list()

  # ===== 6. Loop Through Groups =====
  for (group_dir in group_dirs) {
    focal_group <- basename(group_dir)
    focal_group <- gsub("_vs_rest$", "", focal_group)

    message(sprintf("\n--- Enrichment for: %s ---", focal_group))

    # Read MAST results (FIXED: Flexible column names)
    markers_file <- file.path(group_dir, "MAST_markers.csv")
    if (!file.exists(markers_file)) {
      message("  [SKIP] MAST_markers.csv not found")
      next
    }

    markers <- data.table::fread(markers_file)

    # Standardize column names
    if ("gene" %in% colnames(markers)) {
      markers$primerid <- markers$gene
    }
    if ("feature" %in% colnames(markers)) {
      markers$primerid <- markers$feature
    }
    if ("avg_log2FC" %in% colnames(markers)) {
      markers$logFC <- markers$avg_log2FC
    }
    if ("coef" %in% colnames(markers)) {
      markers$logFC <- markers$coef
    }
    if ("p_val_adj" %in% colnames(markers)) {
      markers$padj <- markers$p_val_adj
    }
    if ("FDR" %in% colnames(markers)) {
      markers$padj <- markers$FDR
    }

    # Get significant up-regulated genes
    sig_up <- markers[regulation == "up" & !is.na(padj)]

    if (nrow(sig_up) < min_genes_for_enrichment) {
      message(sprintf(
        "  [SKIP] Only %d up-regulated genes (need >= %d)",
        nrow(sig_up),
        min_genes_for_enrichment
      ))
      next
    }

    message(sprintf("  Markers: %d up-regulated", nrow(sig_up)))

    # ===== 6a. GO Enrichment (FIXED: universe + mapping) =====
    if (run_go && !is.null(gene_mapping)) {
      message("  Running GO enrichment...")

      tryCatch(
        {
          # Map significant genes to ENTREZ
          sig_genes <- sig_up$primerid
          sig_entrez <- unique(gene_mapping$ENTREZID[
            gene_mapping$feature_raw %in%
              sig_genes |
              gene_mapping$SYMBOL %in% sig_genes
          ])
          sig_entrez <- sig_entrez[!is.na(sig_entrez)]

          # FIXED: Build universe (all tested genes)
          all_tested <- markers$primerid
          universe_entrez <- unique(gene_mapping$ENTREZID[
            gene_mapping$feature_raw %in%
              all_tested |
              gene_mapping$SYMBOL %in% all_tested
          ])
          universe_entrez <- universe_entrez[!is.na(universe_entrez)]

          message(sprintf(
            "    Genes: %d sig / %d universe",
            length(sig_entrez),
            length(universe_entrez)
          ))

          if (length(sig_entrez) < min_genes_for_enrichment) {
            message(sprintf(
              "    [SKIP] Only %d genes with ENTREZ IDs",
              length(sig_entrez)
            ))
          } else {
            # Run enrichGO with universe (FIXED)
            ego <- clusterProfiler::enrichGO(
              gene = sig_entrez,
              universe = universe_entrez, # FIXED: Added universe
              OrgDb = orgdb,
              ont = go_ont,
              pAdjustMethod = "BH",
              pvalueCutoff = pathway_pval_cutoff,
              qvalueCutoff = pathway_pval_cutoff,
              readable = TRUE,
              minGSSize = min_pathway_size,
              maxGSSize = max_pathway_size
            )

            if (!is.null(ego) && nrow(ego@result) > 0) {
              n_pathways <- sum(ego@result$p.adjust < pathway_pval_cutoff)
              message(sprintf("    Found %d enriched GO terms", n_pathways))

              go_result <- as.data.frame(ego@result)
              go_result$focal_group <- focal_group

              group_go_dir <- file.path(go_dir, focal_group)
              dir.create(group_go_dir, showWarnings = FALSE)

              utils::write.csv(
                go_result,
                file.path(group_go_dir, "go_enrichment.csv"),
                row.names = FALSE
              )

              all_go_results[[focal_group]] <- go_result

              # Visualization (FIXED: Use enrichplot)
              if (n_pathways > 0) {
                pdf(
                  file.path(group_go_dir, "go_dotplot.pdf"),
                  width = 10,
                  height = 8
                )
                print(enrichplot::dotplot(
                  ego,
                  showCategory = 20,
                  font.size = 10
                ))
                dev.off()

                pdf(
                  file.path(group_go_dir, "go_barplot.pdf"),
                  width = 10,
                  height = 8
                )
                pdf(
                  file.path(group_go_dir, "go_barplot.pdf"),
                  width = 10,
                  height = 8
                )
                p_bar <- safe_enrich_barplot(
                  ego,
                  showCategory = 20,
                  title = "GO Enrichment"
                )
                if (!is.null(p_bar)) {
                  print(p_bar)
                }
                dev.off()
              }
            } else {
              message("    No significant GO terms")
            }
          }
        },
        error = function(e) {
          message(sprintf("    [ERROR] GO failed: %s", e$message))
        }
      )
    }

    # ===== 6b. GSEA (FIXED: Pre-built TERM2GENE + duplicate aggregation) =====
    if (run_gsea) {
      message("  Running GSEA...")

      tryCatch(
        {
          # Prepare ranked gene list
          gene_list_raw <- markers$logFC
          names(gene_list_raw) <- markers$primerid
          gene_list_raw <- gene_list_raw[!is.na(gene_list_raw)]

          # FIXED: Aggregate duplicates to SYMBOL level
          if (aggregate_duplicates && !is.null(gene_mapping)) {
            # Map to SYMBOL
            gene_to_symbol <- setNames(
              gene_mapping$SYMBOL[match(
                names(gene_list_raw),
                gene_mapping$feature_raw
              )],
              names(gene_list_raw)
            )

            # Aggregate: take max(abs(logFC)) with original sign
            symbol_logfc <- tapply(gene_list_raw, gene_to_symbol, function(x) {
              x[which.max(abs(x))]
            })

            # drop NA symbols
            symbol_logfc <- symbol_logfc[!is.na(names(symbol_logfc))]

            # CRITICAL: force to plain numeric vector (not "array")
            gene_list <- setNames(as.numeric(symbol_logfc), names(symbol_logfc))
            storage.mode(gene_list) <- "numeric"

            # keep only finite
            gene_list <- gene_list[is.finite(gene_list)]

            # sort decreasing as required by GSEA
            gene_list <- sort(gene_list, decreasing = TRUE)

            # ensure unique names
            gene_list <- gene_list[!duplicated(names(gene_list))]

            message(sprintf(
              "    Aggregated: %d features → %d symbols",
              length(gene_list_raw),
              length(gene_list)
            ))
          } else {
            gene_list <- sort(gene_list_raw, decreasing = TRUE)
          }

          # Remove duplicates
          gene_list <- gene_list[!duplicated(names(gene_list))]

          message(sprintf("    Ranked gene list: %d genes", length(gene_list)))

          # GSEA Hallmark (FIXED: Use pre-built TERM2GENE)
          if (!is.null(term2gene_hallmark)) {
            gsea_h <- clusterProfiler::GSEA(
              geneList = gene_list,
              TERM2GENE = term2gene_hallmark, # FIXED: Pre-built
              pvalueCutoff = 1,
              pAdjustMethod = "BH",
              minGSSize = min_pathway_size,
              maxGSSize = max_pathway_size
            )

            if (!is.null(gsea_h) && nrow(gsea_h@result) > 0) {
              n_sig <- sum(gsea_h@result$p.adjust < pathway_pval_cutoff)
              message(sprintf("    Hallmark: %d significant pathways", n_sig))

              gsea_h_result <- as.data.frame(gsea_h@result)
              gsea_h_result$focal_group <- focal_group

              group_gsea_dir <- file.path(gsea_dir, focal_group)
              dir.create(group_gsea_dir, showWarnings = FALSE)

              utils::write.csv(
                gsea_h_result,
                file.path(group_gsea_dir, "gsea_hallmark.csv"),
                row.names = FALSE
              )

              all_gsea_hallmark[[focal_group]] <- gsea_h_result

              if (n_sig > 0) {
                pdf(
                  file.path(group_gsea_dir, "gsea_hallmark_dotplot.pdf"),
                  width = 10,
                  height = 8
                )
                print(enrichplot::dotplot(
                  gsea_h,
                  showCategory = 20,
                  font.size = 10
                ))
                dev.off()
              }
            }
          }

          # GSEA KEGG (FIXED: Use pre-built TERM2GENE)
          if (!is.null(term2gene_kegg)) {
            gsea_k <- clusterProfiler::GSEA(
              geneList = gene_list,
              TERM2GENE = term2gene_kegg, # FIXED: Pre-built
              pvalueCutoff = 1,
              pAdjustMethod = "BH",
              minGSSize = min_pathway_size,
              maxGSSize = max_pathway_size
            )

            if (!is.null(gsea_k) && nrow(gsea_k@result) > 0) {
              n_sig <- sum(gsea_k@result$p.adjust < pathway_pval_cutoff)
              message(sprintf("    KEGG: %d significant pathways", n_sig))

              gsea_k_result <- as.data.frame(gsea_k@result)
              gsea_k_result$focal_group <- focal_group

              group_gsea_dir <- file.path(gsea_dir, focal_group)
              dir.create(group_gsea_dir, showWarnings = FALSE)

              utils::write.csv(
                gsea_k_result,
                file.path(group_gsea_dir, "gsea_kegg.csv"),
                row.names = FALSE
              )

              all_gsea_kegg[[focal_group]] <- gsea_k_result

              if (n_sig > 0) {
                pdf(
                  file.path(group_gsea_dir, "gsea_kegg_dotplot.pdf"),
                  width = 10,
                  height = 8
                )
                print(enrichplot::dotplot(
                  gsea_k,
                  showCategory = 20,
                  font.size = 10
                ))
                dev.off()
              }
            }
          }
        },
        error = function(e) {
          message(sprintf("    [ERROR] GSEA failed: %s", e$message))
        }
      )
    }
  }

  # ===== 7. Cross-Group Pathway Clustering (FIXED: ID not Description, list values_fill) =====
  if (run_clustering) {
    message("\n==== Pathway Clustering Across Groups ====")

    # ===== 7a. GO Clustering =====
    if (length(all_go_results) >= 2) {
      message("Clustering GO pathways...")

      tryCatch(
        {
          go_combined <- data.table::rbindlist(all_go_results, fill = TRUE)
          go_combined <- go_combined[p.adjust < pathway_pval_cutoff, ]

          if (nrow(go_combined) > 0) {
            # FIXED: Use ID instead of Description to avoid duplicates
            pathway_group_mat <- tidyr::pivot_wider(
              go_combined[, c("ID", "focal_group", "p.adjust")],
              names_from = "focal_group",
              values_from = "p.adjust",
              values_fill = list(p.adjust = 1) # FIXED: list form
            )

            pathway_ids <- pathway_group_mat$ID
            pathway_group_mat <- as.matrix(pathway_group_mat[, -1])
            rownames(pathway_group_mat) <- pathway_ids

            # -log10(padj)
            pathway_group_mat <- -log10(pathway_group_mat)
            pathway_group_mat[pathway_group_mat > 10] <- 10

            # Filter: at least 2 groups
            keep_pathways <- rowSums(
              pathway_group_mat > -log10(pathway_pval_cutoff)
            ) >=
              2
            pathway_group_mat <- pathway_group_mat[
              keep_pathways,
              ,
              drop = FALSE
            ]

            if (nrow(pathway_group_mat) >= 3) {
              message(sprintf(
                "  Heatmap: %d pathways × %d groups",
                nrow(pathway_group_mat),
                ncol(pathway_group_mat)
              ))

              # FIXED: Get descriptions for display
              id_to_desc <- setNames(go_combined$Description, go_combined$ID)
              row_labels <- id_to_desc[rownames(pathway_group_mat)]
              row_labels <- ifelse(
                is.na(row_labels),
                rownames(pathway_group_mat),
                row_labels
              )

              pdf(
                file.path(cluster_dir, "go_pathways_heatmap.pdf"),
                width = 12,
                height = max(8, nrow(pathway_group_mat) * 0.2)
              )

              pheatmap::pheatmap(
                pathway_group_mat,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                labels_row = row_labels,
                color = grDevices::colorRampPalette(
                  RColorBrewer::brewer.pal(9, "YlOrRd")
                )(100),
                main = "GO Pathway Enrichment Across Groups",
                fontsize_row = 8,
                fontsize_col = 10,
                cellwidth = 20,
                cellheight = 10
              )

              dev.off()

              utils::write.csv(
                pathway_group_mat,
                file.path(cluster_dir, "go_pathways_matrix.csv")
              )
            }
          }
        },
        error = function(e) {
          message(sprintf("  [ERROR] GO clustering failed: %s", e$message))
        }
      )
    }

    # ===== 7b. GSEA Hallmark Clustering (FIXED: list values_fill) =====
    if (length(all_gsea_hallmark) >= 2) {
      message("Clustering Hallmark pathways...")

      tryCatch(
        {
          gsea_h_combined <- data.table::rbindlist(
            all_gsea_hallmark,
            fill = TRUE
          )
          gsea_h_combined <- gsea_h_combined[p.adjust < pathway_pval_cutoff, ]

          if (nrow(gsea_h_combined) > 0) {
            nes_mat <- tidyr::pivot_wider(
              gsea_h_combined[, c("ID", "focal_group", "NES")],
              names_from = "focal_group",
              values_from = "NES",
              values_fill = list(NES = 0) # FIXED: list form
            )

            pathway_names <- nes_mat$ID
            nes_mat <- as.matrix(nes_mat[, -1])
            rownames(nes_mat) <- pathway_names

            keep_pathways <- rowSums(abs(nes_mat) > 0) >= 2
            nes_mat <- nes_mat[keep_pathways, , drop = FALSE]

            if (nrow(nes_mat) >= 3) {
              message(sprintf(
                "  Heatmap: %d Hallmark pathways × %d groups",
                nrow(nes_mat),
                ncol(nes_mat)
              ))

              pdf(
                file.path(cluster_dir, "gsea_hallmark_heatmap.pdf"),
                width = 12,
                height = max(8, nrow(nes_mat) * 0.2)
              )

              pheatmap::pheatmap(
                nes_mat,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                color = grDevices::colorRampPalette(
                  rev(RColorBrewer::brewer.pal(11, "RdBu"))
                )(100),
                breaks = seq(-3, 3, length.out = 101),
                main = "GSEA Hallmark NES Across Groups",
                fontsize_row = 8,
                fontsize_col = 10,
                cellwidth = 20,
                cellheight = 10
              )

              dev.off()

              utils::write.csv(
                nes_mat,
                file.path(cluster_dir, "gsea_hallmark_nes_matrix.csv")
              )
            }
          }
        },
        error = function(e) {
          message(sprintf(
            "  [ERROR] Hallmark clustering failed: %s",
            e$message
          ))
        }
      )
    }

    # ===== 7c. GSEA KEGG Clustering (FIXED: list values_fill) =====
    if (length(all_gsea_kegg) >= 2) {
      message("Clustering KEGG pathways...")

      tryCatch(
        {
          gsea_k_combined <- data.table::rbindlist(all_gsea_kegg, fill = TRUE)
          gsea_k_combined <- gsea_k_combined[p.adjust < pathway_pval_cutoff, ]

          if (nrow(gsea_k_combined) > 0) {
            nes_mat <- tidyr::pivot_wider(
              gsea_k_combined[, c("ID", "focal_group", "NES")],
              names_from = "focal_group",
              values_from = "NES",
              values_fill = list(NES = 0) # FIXED: list form
            )

            pathway_names <- nes_mat$ID
            nes_mat <- as.matrix(nes_mat[, -1])
            rownames(nes_mat) <- pathway_names

            keep_pathways <- rowSums(abs(nes_mat) > 0) >= 2
            nes_mat <- nes_mat[keep_pathways, , drop = FALSE]

            if (nrow(nes_mat) >= 3) {
              message(sprintf(
                "  Heatmap: %d KEGG pathways × %d groups",
                nrow(nes_mat),
                ncol(nes_mat)
              ))

              pdf(
                file.path(cluster_dir, "gsea_kegg_heatmap.pdf"),
                width = 14,
                height = max(10, nrow(nes_mat) * 0.15)
              )

              pheatmap::pheatmap(
                nes_mat,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                color = grDevices::colorRampPalette(
                  rev(RColorBrewer::brewer.pal(11, "RdBu"))
                )(100),
                breaks = seq(-3, 3, length.out = 101),
                main = "GSEA KEGG NES Across Groups",
                fontsize_row = 6,
                fontsize_col = 10,
                cellwidth = 20,
                cellheight = 8
              )

              dev.off()

              utils::write.csv(
                nes_mat,
                file.path(cluster_dir, "gsea_kegg_nes_matrix.csv")
              )
            }
          }
        },
        error = function(e) {
          message(sprintf("  [ERROR] KEGG clustering failed: %s", e$message))
        }
      )
    }
  }

  # ===== 8. Summary Report (FIXED: sink with on.exit) =====
  message("\n==== Creating Enrichment Summary ====")

  summary_file <- file.path(enrich_dir, "enrichment_summary.txt")

  sink(summary_file)
  on.exit(sink(), add = TRUE) # FIXED: on.exit()

  cat("=== One-vs-Rest Enrichment Analysis Summary ===\n\n")
  cat(sprintf("Analysis date: %s\n", Sys.Date()))
  cat(sprintf("Groups analyzed: %d\n", length(group_dirs)))
  cat(sprintf("Organism: %s\n", organism))
  cat(sprintf("Gene ID type: %s\n", gene_id_type))
  cat("\n")

  if (run_go) {
    cat("GO Enrichment:\n")
    cat(sprintf("  - Groups with GO results: %d\n", length(all_go_results)))
    if (length(all_go_results) > 0) {
      for (g in names(all_go_results)) {
        n_terms <- sum(all_go_results[[g]]$p.adjust < pathway_pval_cutoff)
        cat(sprintf("  - %s: %d enriched GO terms\n", g, n_terms))
      }
    }
    cat("\n")
  }

  if (run_gsea) {
    cat("GSEA Hallmark:\n")
    cat(sprintf("  - Groups with results: %d\n", length(all_gsea_hallmark)))
    if (length(all_gsea_hallmark) > 0) {
      for (g in names(all_gsea_hallmark)) {
        n_sig <- sum(all_gsea_hallmark[[g]]$p.adjust < pathway_pval_cutoff)
        cat(sprintf("  - %s: %d significant pathways\n", g, n_sig))
      }
    }
    cat("\n")

    cat("GSEA KEGG:\n")
    cat(sprintf("  - Groups with results: %d\n", length(all_gsea_kegg)))
    if (length(all_gsea_kegg) > 0) {
      for (g in names(all_gsea_kegg)) {
        n_sig <- sum(all_gsea_kegg[[g]]$p.adjust < pathway_pval_cutoff)
        cat(sprintf("  - %s: %d significant pathways\n", g, n_sig))
      }
    }
    cat("\n")
  }

  if (run_clustering) {
    cat("Pathway Clustering:\n")
    cat("  - Heatmaps generated for cross-group comparison\n")
    cat("  - See pathway_clustering/ directory\n")
  }

  sink() # Explicit sink() before on.exit()

  message(sprintf("Summary saved: %s", summary_file))

  # ===== 9. Final Message =====
  message("\n==== Enrichment Analysis Complete ====")
  message(sprintf("Results: %s", enrich_dir))

  if (run_go) {
    message("  ✓ GO enrichment (with universe correction)")
  }
  if (run_gsea) {
    message("  ✓ GSEA Hallmark + KEGG (aggregated symbols)")
  }
  if (run_clustering) {
    message("  ✓ Pathway clustering heatmaps")
  }

  message("\nKey outputs:")
  message("  - enrichment_analysis/go_enrichment/[group]/")
  message("  - enrichment_analysis/gsea/[group]/")
  message("  - enrichment_analysis/pathway_clustering/*.pdf")

  invisible(TRUE)
}
