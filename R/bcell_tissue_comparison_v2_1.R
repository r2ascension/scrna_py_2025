#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.1 - STANDALONE PRODUCTION
# ==============================================================================
#
# Purpose:
#   B cell lineage tissue comparison:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. Visualization (UMAP, marker dotplot/heatmap, composition)
#     3. Pseudobulk DESeq2 per L2 subtype x tissue pair
#     4. Exploratory cell-level wilcox (marker discovery only)
#     5. GMT-based enrichment (GO BP/MF/CC, KEGG, Hallmark) - symbol direct
#     6. interpret_agent (DeepSeek) on original enrichResult (optional)
#     7. Structured REPORT.md (png embeds)
#
# Standalone Notes:
#   - This script processes B cells ONLY; copy and adapt for other lineages
#   - All GMT enrichment is symbol-based (no bitr / no gene ID loss)
#   - Primary DE = pseudobulk DESeq2 (Squair et al. 2021)
#   - Cell-level wilcox is supplementary (NOT for inference)
#   - LLM (interpret_agent) is optional; runs if DEEPSEEK_API_KEY is set
#
# v2.1 Fixes (vs v2.0):
#   [P0]   sample col added to required column check + NA filter
#   [P1-2] enrichment universe = DESeq2 tested genes (not GMT full set)
#   [P1-3] DESeq2 design formula configurable (default ~ tissue)
#   [P1-4] WILCOX_LFC_THR now actually applied to filtering
#   [P1-5] LLM optional — no API key = skip interpret_agent, not crash
#   [P2-7] set.seed(42) for reproducible heatmap/downsampling
#   [P2-8] Added sample-level composition boxplot
#   [+]    counts layer check after GetSeurat
#
# Author: r2end
# Date:   2026-03-23
# Memory: < 30GB (B cell subset)
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
# 1. Configuration
# ==============================================================================

# ----- Input / Output -----
# UPDATE this to your actual B cell final h5ad
H5AD_PATH  <- "/home/h2048/data/py/0208/merged_scanvi_L2_prod_v1/merged_scanvi_L2_prod.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0323/bcell_tissue_comparison_v2"

# ----- GMT Databases (symbol-based, zero gene loss) -----
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL      <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# ----- DeepSeek API (OPTIONAL - pipeline runs without it) -----
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
ENABLE_LLM <- nchar(DEEPSEEK_API_KEY) >= 10
if (!ENABLE_LLM) {
  cat("[WARN] DEEPSEEK_API_KEY not set. interpret_agent will be SKIPPED.\n")
  cat("       Set via: export DEEPSEEK_API_KEY='your-key' to enable LLM.\n")
}
INTERPRET_AGENT_MODEL      <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS <- 50
INTERPRET_AGENT_ADD_PPI    <- FALSE

# ----- Python / reticulate -----
PYTHON_CONDA_ENV <- "bbknn_env"

# ----- Metadata Column Names -----
TISSUE_COL      <- "tissue"
SAMPLE_COL      <- "sample"
CELLTYPE_L2_COL <- "cell_type_L2"      # broad B cell subtype
CELLTYPE_L3_COL <- "cell_type_L3"      # fine subtype (display label)
LABEL_COL       <- "cell_type_L3"      # annotation column for visualization

CELLTYPE_L2_COL_CANDIDATES <- c(
  "Cell_Type_L2",
  "Cell_Type_L2_final",
  "cell_type_level_2",
  "Cell_type_annotation_level2",
  "Cell_Type_L2_pred",
  "cell_type_scanvi_corrected",
  "cell_type"
)
CELLTYPE_L3_COL_CANDIDATES <- c(
  "cell_type_L3",
  "cell_type_level_3",
  "Cell_type_annotation_level3",
  "cell_type_expert",
  "cell_type_scanvi_pred",
  "predicted_labels"
)
LABEL_COL_CANDIDATES <- c(
  "cell_type_L3",
  "cell_type_level_3",
  "Cell_type_annotation_level3",
  "cell_type_expert",
  "Cell_Type_L2",
  "cell_type_level_2"
)

UMAP_REDUCTION_PREFERRED <- c("umap_scanvi", "umap")
UMAP_PT_SIZE <- 0.35
UMAP_ALPHA   <- 0.9
UMAP_TISSUE_COLORS <- c(
  "lung parenchyma"   = "#D55E00",
  "nose"              = "#009E73",
  "respiratory airway" = "#0072B2",
  "sinus"             = "#CC79A7"
)

# ----- Pseudobulk DE Parameters -----
MIN_CELLS_PER_PSEUDOBULK <- 10
MIN_SAMPLES_PER_TISSUE   <- 3
PADJ_THR                 <- 0.05
LFC_THR                  <- 1.0    # log2FC for DESeq2 (stringent)

# DESeq2 design formula (configurable for multi-study data)
# Options:
#   ~ tissue                      (single study, or tissue not confounded with batch)
#   ~ dataset + tissue            (multi-study: block on dataset/center)
#   ~ donor + tissue              (paired design: same donor, multiple tissues)
# Set DESIGN_COVARIATE to NULL for simple ~ tissue, or to a column name to add it.
DESIGN_COVARIATE <- NULL  # e.g., "dataset" for ~ dataset + tissue

# ----- Reproducibility -----
set.seed(42)

# ----- Exploratory Wilcox Parameters -----
WILCOX_MIN_CELLS <- 50
WILCOX_LFC_THR   <- 0.25
WILCOX_TOP_N     <- 200

# ----- Enrichment -----
TOP_N_DEG_ENRICHMENT <- 200

# ----- Visualization -----
HEATMAP_CELLS_PER_TYPE <- 100

# ----- B Cell Known Markers -----
# From bcell_markers_comprehensive.csv + standard panels
KNOWN_MARKERS <- c(
  # Pan-B
  "CD79A", "CD79B", "MS4A1", "CD19", "PAX5",
  # Naive
  "IGHD", "IGHM", "TCL1A", "FCER2", "IL4R",
  # Memory
  "CD27", "TNFRSF13B", "AIM2",
  # GC
  "BCL6", "AICDA", "RGS13", "MEF2B", "MME", "MKI67",
  # Plasma / Plasmablast
  "PRDM1", "XBP1", "MZB1", "JCHAIN", "SDC1", "DERL3", "SSR4",
  # Isotypes
  "IGHA1", "IGHG1", "IGHE",
  # Atypical / ABC
  "ITGAX", "TBX21", "FCRL4", "FCRL5",
  # Activation / IFN
  "CD69", "CD86", "ISG15",
  # CRSwNP-specific
  "LAPTM5", "CD74"
)

# ----- Disease Context for LLM -----
BCELL_CONTEXT <- paste(
  "B cells and plasma cells from CRSwNP nasal polyps and healthy nasal tissue.",
  "Key subtypes: Naive B, Memory B (IgM/IgG/IgA/IgE-switched), GC B (DZ/LZ),",
  "Plasmablast, Plasma cell (IgA/IgG/IgE), Atypical/ABC, Breg.",
  "CRSwNP-specific: local IgE production (1/5 plasma cells in NP),",
  "germinal center reactions in polyp tissue, MZB1+ population,",
  "FCRL4+ tissue-resident atypical B cells, extrafollicular IgE switching.",
  "Focus: B cell differentiation, local antibody production,",
  "class switch recombination, mucosal immunity in Type 2 inflammation."
)

# ==============================================================================
# 2. Load Libraries
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("B Cell Tissue Comparison Pipeline v2.0 (STANDALONE)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

suppressPackageStartupMessages({
  library(reticulate)
  library(SCNT)
  library(Seurat)
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(fanyi)
})

# Configure Python bridge for h5ad import
if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
  py_bin <- Sys.getenv("RETICULATE_PYTHON")
} else {
  use_condaenv(PYTHON_CONDA_ENV, required = TRUE)
  py_bin <- tryCatch(py_config()$python, error = function(e) paste0("conda:", PYTHON_CONDA_ENV))
}

# Prefer local GetSeurat implementation (handles h5ad obsm -> matrix safely)
LOCAL_GETSEURAT <- "/home/h2048/script/R/GetSeurat.R"
if (file.exists(LOCAL_GETSEURAT)) {
  source(LOCAL_GETSEURAT)
}

# Configure DeepSeek (only if enabled)
if (ENABLE_LLM) {
  fanyi::set_translate_option(key = DEEPSEEK_API_KEY, source = "deepseek")
  cat("[OK] LLM enabled (DeepSeek)\n")
} else {
  cat("[INFO] LLM disabled (no API key)\n")
}
cat(sprintf("[OK] Python configured: %s\n", py_bin))
if (file.exists(LOCAL_GETSEURAT)) cat(sprintf("[OK] Local GetSeurat loaded: %s\n", LOCAL_GETSEURAT))
cat("[OK] Libraries loaded\n\n")

# ==============================================================================
# 3. Load GMT Databases
# ==============================================================================

cat("=== Loading GMT Databases ===\n")

gmt_all <- tryCatch({
  g <- read.gmt(MSIGDB_GMT_PATH); g$gene <- toupper(g$gene)
  cat(sprintf("[OK] MSigDB: %d gene sets\n", length(unique(g$term)))); g
}, error = function(e) { cat("[WARN] MSigDB failed\n"); NULL })

hallmark_t2g <- if (!is.null(gmt_all)) {
  h <- gmt_all %>% filter(grepl("^HALLMARK_", term))
  if (nrow(h) > 0) { cat(sprintf("[OK] Hallmark: %d\n", length(unique(h$term)))); h } else NULL
} else NULL

kegg_t2g <- if (!is.null(gmt_all)) {
  k <- gmt_all %>% filter(grepl("^KEGG_", term))
  if (nrow(k) > 0) { cat(sprintf("[OK] KEGG: %d\n", length(unique(k$term)))); k } else NULL
} else NULL

go_bp_t2g <- NULL; go_mf_t2g <- NULL; go_cc_t2g <- NULL
tryCatch({
  gmt_go <- read.gmt(GMT_GO_ALL); gmt_go$gene <- toupper(gmt_go$gene)
  go_bp_t2g <- gmt_go %>% filter(grepl("^GOBP_", term))
  go_mf_t2g <- gmt_go %>% filter(grepl("^GOMF_", term))
  go_cc_t2g <- gmt_go %>% filter(grepl("^GOCC_", term))
  cat(sprintf("[OK] GO: BP=%d, MF=%d, CC=%d\n",
              length(unique(go_bp_t2g$term)),
              length(unique(go_mf_t2g$term)),
              length(unique(go_cc_t2g$term))))
}, error = function(e) cat("[WARN] GO GMT failed\n"))

cat("\n")

# ==============================================================================
# 4. Helper Functions
# ==============================================================================

safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

first_nonempty_column <- function(meta, candidates, min_unique = 2) {
  for (nm in candidates) {
    if (!nm %in% colnames(meta)) next
    vals <- meta[[nm]]
    vals <- vals[!is.na(vals) & trimws(as.character(vals)) != ""]
    if (length(unique(as.character(vals))) >= min_unique) return(nm)
  }
  NULL
}

pick_reduction <- function(obj, preferred) {
  red <- Reductions(obj)
  hit <- preferred[preferred %in% red]
  if (length(hit) > 0) hit[[1]] else NULL
}

build_umap_plot <- function(obj, reduction_name, group_col, title,
                            split_col = NULL, label = FALSE, width = 10, height = 8,
                            cols = NULL) {
  p <- DimPlot(
    obj,
    reduction = reduction_name,
    group.by = group_col,
    split.by = split_col,
    pt.size = UMAP_PT_SIZE,
    raster = FALSE,
    shuffle = TRUE,
    alpha = UMAP_ALPHA,
    label = label,
    repel = label,
    cols = cols
  ) +
    ggtitle(title) +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )

  list(plot = p, width = width, height = height)
}

# ----- 4.1 GMT enrichment (symbol-direct, universe = tested genes) -----
run_gmt_enrichment <- function(gene_list, t2g, db_name, tested_genes = NULL) {
  if (is.null(t2g) || nrow(t2g) == 0 || length(gene_list) < 5) return(NULL)

  # P1-2 fix: universe = genes actually tested in DESeq2, intersected with GMT
  # This avoids inflated significance from using the entire GMT gene space
  if (!is.null(tested_genes)) {
    universe_use <- intersect(toupper(tested_genes), unique(t2g$gene))
  } else {
    universe_use <- unique(t2g$gene)  # fallback if no tested genes provided
  }

  tryCatch(
    enricher(
      gene = toupper(unique(gene_list)), TERM2GENE = t2g,
      universe = universe_use,
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH",
      minGSSize = 10, maxGSSize = 500
    ),
    error = function(e) { cat(sprintf("    [WARN] %s: %s\n", db_name, e$message)); NULL }
  )
}

# ----- 4.2 Pseudobulk aggregation -----
aggregate_pseudobulk <- function(obj, celltype_l2) {
  cells <- colnames(obj)[obj@meta.data[[CELLTYPE_L2_COL]] == celltype_l2]
  if (length(cells) < 30) return(NULL)

  sub <- subset(obj, cells = cells)
  meta_sub <- sub@meta.data
  counts_mat <- GetAssayData(sub, layer = "counts")

  meta_sub$pb_group <- paste0(meta_sub[[SAMPLE_COL]], "__", meta_sub[[TISSUE_COL]])
  groups <- unique(meta_sub$pb_group)

  pb_counts <- matrix(0, nrow = nrow(counts_mat), ncol = length(groups))
  rownames(pb_counts) <- rownames(counts_mat)
  colnames(pb_counts) <- groups
  pb_meta <- data.frame(row.names = groups, stringsAsFactors = FALSE)

  for (g in groups) {
    cells_g <- rownames(meta_sub)[meta_sub$pb_group == g]
    if (length(cells_g) < MIN_CELLS_PER_PSEUDOBULK) next

    pb_counts[, g] <- if (length(cells_g) == 1) {
      as.numeric(counts_mat[, cells_g])
    } else {
      Matrix::rowSums(counts_mat[, cells_g])
    }
    pb_meta[g, "tissue"]  <- meta_sub[[TISSUE_COL]][meta_sub$pb_group == g][1]
    pb_meta[g, "sample"]  <- meta_sub[[SAMPLE_COL]][meta_sub$pb_group == g][1]
    pb_meta[g, "n_cells"] <- length(cells_g)
  }

  valid <- colnames(pb_counts)[colSums(pb_counts) > 0 & !is.na(pb_meta$tissue)]
  if (length(valid) < 4) return(NULL)
  list(counts = pb_counts[, valid, drop = FALSE], meta = pb_meta[valid, , drop = FALSE])
}

# ----- 4.3 DESeq2 pairwise -----
run_deseq2_pairwise <- function(pb, t1, t2) {
  keep <- pb$meta$tissue %in% c(t1, t2)
  if (sum(keep) < 4) return(NULL)

  counts_sub <- pb$counts[, keep, drop = FALSE]
  meta_sub   <- pb$meta[keep, , drop = FALSE]
  tissue_levels_safe <- make.names(c(t1, t2), unique = TRUE)
  meta_sub$tissue <- droplevels(factor(meta_sub$tissue, levels = c(t1, t2)))
  meta_sub$tissue_safe <- factor(
    ifelse(meta_sub$tissue == t1, tissue_levels_safe[1], tissue_levels_safe[2]),
    levels = tissue_levels_safe
  )

  n1 <- sum(meta_sub$tissue == t1); n2 <- sum(meta_sub$tissue == t2)
  if (n1 < MIN_SAMPLES_PER_TISSUE || n2 < MIN_SAMPLES_PER_TISSUE) return(NULL)

  keep_genes <- rowSums(counts_sub >= 1) >= max(3, ncol(counts_sub) * 0.2)
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  if (nrow(counts_sub) < 100) return(NULL)

  # P1-2: record tested genes for enrichment universe
  tested_genes <- rownames(counts_sub)

  counts_int <- round(counts_sub); storage.mode(counts_int) <- "integer"

  # P1-3: build design formula (configurable covariate blocking)
  design_formula <- ~ tissue_safe
  if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
    n_levels <- length(unique(meta_sub[[DESIGN_COVARIATE]]))
    if (n_levels >= 2 && n_levels < nrow(meta_sub)) {
      meta_sub[[DESIGN_COVARIATE]] <- factor(meta_sub[[DESIGN_COVARIATE]])
      design_formula <- as.formula(paste("~", DESIGN_COVARIATE, "+ tissue_safe"))
      cat(sprintf("    Design: %s (covariate: %d levels)\n",
                  deparse(design_formula), n_levels))
    } else {
      cat(sprintf("    Design: ~ tissue_safe (covariate '%s' has %d level(s), skipped)\n",
                  DESIGN_COVARIATE, n_levels))
    }
  }

  dds <- tryCatch({
    dds <- DESeqDataSetFromMatrix(counts_int, meta_sub, design_formula)
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    # Fallback to simple ~ tissue_safe if design is rank-deficient
    if (!is.null(DESIGN_COVARIATE)) {
      cat(sprintf("    [WARN] Design failed (%s), falling back to ~ tissue_safe\n", e$message))
      tryCatch({
        dds2 <- DESeqDataSetFromMatrix(counts_int, meta_sub, ~ tissue_safe)
        DESeq(dds2, quiet = TRUE)
      }, error = function(e2) { cat(sprintf("    [ERROR] DESeq2: %s\n", e2$message)); NULL })
    } else {
      cat(sprintf("    [ERROR] DESeq2: %s\n", e$message)); NULL
    }
  })
  if (is.null(dds)) return(NULL)

  res <- results(dds, contrast = c("tissue_safe", tissue_levels_safe[2], tissue_levels_safe[1]), alpha = PADJ_THR)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(
      sig = ifelse(padj < PADJ_THR & abs(log2FoldChange) > LFC_THR, "sig", "ns"),
      direction = ifelse(log2FoldChange > 0, "up", "down")
    )

  list(
    de_table = res_df,
    tested_genes = tested_genes,
    n_up = sum(res_df$sig == "sig" & res_df$direction == "up"),
    n_down = sum(res_df$sig == "sig" & res_df$direction == "down"),
    n_samples_1 = n1, n_samples_2 = n2, tissue_1 = t1, tissue_2 = t2,
    design = deparse(design_formula)
  )
}

# ----- 4.4 Exploratory wilcox -----
run_wilcox_exploratory <- function(obj, celltype_l2) {
  cells <- colnames(obj)[obj@meta.data[[CELLTYPE_L2_COL]] == celltype_l2]
  if (length(cells) < WILCOX_MIN_CELLS) return(NULL)

  sub <- subset(obj, cells = cells)
  tissues <- sort(unique(na.omit(sub@meta.data[[TISSUE_COL]])))
  if (length(tissues) < 2) return(NULL)

  Idents(sub) <- TISSUE_COL
  results <- list()
  for (pair in combn(tissues, 2, simplify = FALSE)) {
    t1 <- pair[1]; t2 <- pair[2]
    n1 <- sum(sub@meta.data[[TISSUE_COL]] == t1, na.rm = TRUE)
    n2 <- sum(sub@meta.data[[TISSUE_COL]] == t2, na.rm = TRUE)
    if (n1 < WILCOX_MIN_CELLS || n2 < WILCOX_MIN_CELLS) next

    de <- tryCatch(
      FindMarkers(sub, ident.1 = t2, ident.2 = t1, test.use = "wilcox",
                  min.pct = 0.1, logfc.threshold = 0, max.cells.per.ident = 5000),
      error = function(e) NULL
    )
    if (!is.null(de) && nrow(de) > 0) {
      de$gene <- rownames(de)
      de$padj <- p.adjust(de$p_val, method = "BH")
      # P1-4 fix: actually apply WILCOX_LFC_THR filter
      results[[paste0(t2, "_vs_", t1)]] <- de %>%
        filter(abs(avg_log2FC) >= WILCOX_LFC_THR) %>%
        arrange(padj) %>%
        head(WILCOX_TOP_N)
    }
  }
  results
}

# ----- 4.5 interpret_agent (original enrichResult, no re-enrichment) -----
run_interpret_agent_safe <- function(enrich_obj, context_str, gene_fc = NULL) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(NULL)
  tryCatch(
    clusterProfiler::interpret_agent(
      x = enrich_obj, context = context_str,
      n_pathways = INTERPRET_AGENT_N_PATHWAYS,
      model = INTERPRET_AGENT_MODEL, api_key = DEEPSEEK_API_KEY,
      add_ppi = INTERPRET_AGENT_ADD_PPI, gene_fold_change = gene_fc
    ),
    error = function(e) { cat(sprintf("    [ERROR] agent: %s\n", e$message)); NULL }
  )
}

to_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.character(x) && length(x) == 1) return(x)
  if (is.character(x)) return(paste(x, collapse = "; "))
  if (is.list(x)) return(paste(capture.output(str(x, max.level = 2)), collapse = "\n"))
  as.character(x)
}

# ----- 4.6 Stratified downsampling -----
stratified_downsample <- function(obj, group_col, n_per = HEATMAP_CELLS_PER_TYPE) {
  meta <- obj@meta.data
  cells <- unlist(lapply(unique(meta[[group_col]]), function(g) {
    gc <- rownames(meta)[meta[[group_col]] == g]
    sample(gc, min(length(gc), n_per))
  }))
  subset(obj, cells = cells)
}

# ----- 4.7 Save pdf + png -----
save_plot <- function(p, path_no_ext, width = 10, height = 8) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 150)
}

# ==============================================================================
# 5. Load Data
# ==============================================================================

cat("=== Loading B Cell Data ===\n")

if (!file.exists(H5AD_PATH)) stop(sprintf("File not found: %s", H5AD_PATH))
obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
cat(sprintf("[OK] %d cells x %d genes\n", ncol(obj), nrow(obj)))

# Verify counts layer exists
if (!"counts" %in% Layers(obj[["RNA"]])) {
  stop("RNA assay missing 'counts' layer. Check GetSeurat output.")
}
cat("[OK] counts layer verified\n\n")

# Create output dirs
fig_dir <- file.path(OUTPUT_DIR, "figures")
rpt_dir <- file.path(OUTPUT_DIR, "reports")
de_dir  <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir  <- file.path(OUTPUT_DIR, "wilcox_exploratory")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rpt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(wx_dir,  recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. Validate & Clean Metadata
# ==============================================================================

cat("=== Validating Metadata ===\n")

meta <- obj@meta.data

resolved_l2_col <- first_nonempty_column(meta, c(CELLTYPE_L2_COL, CELLTYPE_L2_COL_CANDIDATES))
resolved_l3_col <- first_nonempty_column(meta, c(CELLTYPE_L3_COL, CELLTYPE_L3_COL_CANDIDATES))
resolved_label_col <- first_nonempty_column(meta, c(LABEL_COL, LABEL_COL_CANDIDATES, resolved_l3_col, resolved_l2_col), min_unique = 1)

if (!is.null(resolved_l2_col) && resolved_l2_col != CELLTYPE_L2_COL) {
  cat(sprintf("[INFO] Using L2 column: %s (instead of %s)\n", resolved_l2_col, CELLTYPE_L2_COL))
  CELLTYPE_L2_COL <- resolved_l2_col
}
if (!is.null(resolved_l3_col) && resolved_l3_col != CELLTYPE_L3_COL) {
  cat(sprintf("[INFO] Using L3 column: %s (instead of %s)\n", resolved_l3_col, CELLTYPE_L3_COL))
  CELLTYPE_L3_COL <- resolved_l3_col
}
if (!is.null(resolved_label_col) && resolved_label_col != LABEL_COL) {
  cat(sprintf("[INFO] Using label column: %s (instead of %s)\n", resolved_label_col, LABEL_COL))
  LABEL_COL <- resolved_label_col
}

# P0 fix: check ALL required columns including SAMPLE_COL
for (col in c(TISSUE_COL, SAMPLE_COL, CELLTYPE_L2_COL, LABEL_COL)) {
  if (!col %in% colnames(meta)) stop(sprintf("Missing required column: %s", col))
}

# P0 fix: drop NA/empty for tissue + sample + L2 + label
bad_idx <- is.na(meta[[TISSUE_COL]]) | trimws(as.character(meta[[TISSUE_COL]])) == "" |
           is.na(meta[[SAMPLE_COL]]) | trimws(as.character(meta[[SAMPLE_COL]])) == "" |
           is.na(meta[[CELLTYPE_L2_COL]]) | trimws(as.character(meta[[CELLTYPE_L2_COL]])) == "" |
           is.na(meta[[LABEL_COL]])  | trimws(as.character(meta[[LABEL_COL]])) == ""
if (sum(bad_idx) > 0) {
  cat(sprintf("[INFO] Dropping %d cells with NA tissue/sample/L2/label\n", sum(bad_idx)))
  obj <- subset(obj, cells = rownames(meta)[!bad_idx])
}

# Verify sufficient samples
n_samples <- length(unique(obj@meta.data[[SAMPLE_COL]]))
if (n_samples < 2) stop(sprintf("Only %d unique sample(s). Need >= 2 for pseudobulk DE.", n_samples))
cat(sprintf("[OK] %d unique samples\n", n_samples))

# Ensure L2 column
l2_is_fallback <- FALSE
if (!CELLTYPE_L2_COL %in% colnames(obj@meta.data)) {
  cat(sprintf("[WARN] %s not found; using %s as L2 proxy\n", CELLTYPE_L2_COL, LABEL_COL))
  obj@meta.data[[CELLTYPE_L2_COL]] <- obj@meta.data[[LABEL_COL]]
  l2_is_fallback <- TRUE
}

# Re-fetch meta after modifications
meta <- obj@meta.data

tissues  <- sort(unique(na.omit(meta[[TISSUE_COL]])))
l2_types <- sort(unique(na.omit(meta[[CELLTYPE_L2_COL]])))

cat(sprintf("[OK] Cells: %d\n", ncol(obj)))
cat(sprintf("[OK] Tissues: %s\n", paste(tissues, collapse = ", ")))
cat(sprintf("[OK] L2 types (%d): %s\n", length(l2_types), paste(l2_types, collapse = ", ")))
cat(sprintf("[OK] Active columns: L2=%s | L3=%s | Label=%s\n", CELLTYPE_L2_COL, CELLTYPE_L3_COL, LABEL_COL))
if (l2_is_fallback) cat("[WARN] L2 = L3 fallback active\n")
cat("\n")

# Ensure normalized data layer exists for visualization / marker finding
cat("[INFO] Running NormalizeData() for visualization and marker analyses\n")
obj <- NormalizeData(obj, verbose = FALSE)
cat("[OK] data layer ready\n\n")

# ==============================================================================
# 7. Visualization
# ==============================================================================

cat("=== Visualization ===\n")

# 7.1 UMAP
umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
if (!is.null(umap_reduction)) {
  tissue_cols_use <- UMAP_TISSUE_COLORS[names(UMAP_TISSUE_COLORS) %in% tissues]
  p1 <- build_umap_plot(
    obj, umap_reduction, TISSUE_COL,
    title = sprintf("B Cell - Tissue (%s)", umap_reduction),
    cols = tissue_cols_use,
    width = 11, height = 8
  )
  save_plot(p1$plot, file.path(fig_dir, "umap_tissue"), width = p1$width, height = p1$height)

  p2 <- build_umap_plot(
    obj, umap_reduction, CELLTYPE_L2_COL,
    title = sprintf("B Cell - Cell Type L2 (%s)", umap_reduction),
    label = TRUE,
    width = 16, height = 10
  )
  save_plot(p2$plot, file.path(fig_dir, "umap_celltype"), width = p2$width, height = p2$height)

  p2s <- build_umap_plot(
    obj, umap_reduction, CELLTYPE_L2_COL,
    title = sprintf("B Cell - Cell Type L2 by Tissue (%s)", umap_reduction),
    split_col = TISSUE_COL,
    label = TRUE,
    width = max(12, 4.5 * length(tissues)), height = 8
  )
  save_plot(p2s$plot, file.path(fig_dir, "umap_celltype_split_tissue"),
            width = p2s$width, height = p2s$height)

  if (LABEL_COL != CELLTYPE_L2_COL) {
    p2_l3 <- build_umap_plot(
      obj, umap_reduction, LABEL_COL,
      title = sprintf("B Cell - Cell Type L3 (%s)", umap_reduction),
      label = TRUE,
      width = 18, height = 12
    )
    save_plot(p2_l3$plot, file.path(fig_dir, "umap_subtype"),
              width = p2_l3$width, height = p2_l3$height)
  }

  cat(sprintf("[OK] UMAP saved using reduction: %s\n", umap_reduction))
} else {
  cat("[WARN] No UMAP reduction\n")
}

# 7.2 Dotplot
markers_present <- intersect(KNOWN_MARKERS, rownames(obj))
if (length(markers_present) >= 3) {
  Idents(obj) <- LABEL_COL
  p3 <- DotPlot(obj, features = markers_present) + RotatedAxis() +
    ggtitle("B Cell - Known Markers") +
    theme(axis.text.x = element_text(size = 7))
  save_plot(p3, file.path(fig_dir, "dotplot_markers"),
            width = max(12, length(markers_present) * 0.45),
            height = max(6, length(l2_types) * 0.4))
  cat("[OK] Dotplot saved\n")
}

# 7.3 Composition barplot
comp_df <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL))) %>%
  count(!!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(TISSUE_COL)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p4 <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = pct,
                           fill = !!sym(CELLTYPE_L2_COL))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Tissue", y = "Percentage (%)", title = "B Cell - Composition by Tissue") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4, file.path(fig_dir, "composition_tissue_celltype"))

# Absolute count version
p4b <- ggplot(comp_df, aes(x = !!sym(TISSUE_COL), y = n,
                            fill = !!sym(CELLTYPE_L2_COL))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Tissue", y = "Cell Count", title = "B Cell - Absolute Count by Tissue") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4b, file.path(fig_dir, "count_tissue_celltype"))
cat("[OK] Composition plots saved\n")

# 7.3b Sample-level composition (P2-8: avoids cell-count dominance)
sample_comp <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL)), !is.na(!!sym(SAMPLE_COL))) %>%
  count(!!sym(SAMPLE_COL), !!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(SAMPLE_COL)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

if (nrow(sample_comp) > 0) {
  facet_ncol <- min(4, length(l2_types))
  facet_nrow <- ceiling(length(l2_types) / facet_ncol)
  p4c <- ggplot(sample_comp, aes(x = !!sym(TISSUE_COL), y = pct,
                                  fill = !!sym(TISSUE_COL))) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(as.formula(paste("~", CELLTYPE_L2_COL)), scales = "free_y", ncol = facet_ncol) +
    labs(x = "Tissue", y = "Proportion per Sample (%)",
         title = "B Cell - Sample-level Composition by Tissue") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                            legend.position = "none")
  save_plot(p4c, file.path(fig_dir, "composition_sample_level"),
            width = max(10, facet_ncol * 4), height = max(6, facet_nrow * 3.5))
  cat("[OK] Sample-level composition saved\n")
}

# 7.4 Heatmap (top markers per L2, stratified downsampling)
Idents(obj) <- CELLTYPE_L2_COL
top_mk <- tryCatch(
  FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                 max.cells.per.ident = 500, test.use = "wilcox"),
  error = function(e) { cat("[WARN] FindAllMarkers failed\n"); NULL }
)

if (!is.null(top_mk) && nrow(top_mk) > 0) {
  fwrite(top_mk, file.path(rpt_dir, "all_markers_per_L2.csv"))
  top10 <- top_mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()

  obj_ds <- stratified_downsample(obj, CELLTYPE_L2_COL, HEATMAP_CELLS_PER_TYPE)
  obj_ds <- ScaleData(obj_ds, features = unique(top10$gene), verbose = FALSE)
  p5 <- DoHeatmap(obj_ds, features = unique(top10$gene), size = 3) +
    ggtitle("B Cell - Top Markers per L2")
  save_plot(p5, file.path(fig_dir, "heatmap_top_markers"), width = 14, height = 10)
  rm(obj_ds); gc()
  cat("[OK] Heatmap saved\n")
}

cat("\n")

# ==============================================================================
# 8. Pseudobulk DESeq2 Tissue Comparison
# ==============================================================================

cat("=== Pseudobulk DESeq2 Tissue Comparison ===\n")
cat(sprintf("Thresholds: padj < %s, |log2FC| > %s\n\n", PADJ_THR, LFC_THR))

pb_de_all   <- list()
enrich_all  <- list()
agent_all   <- list()

for (ct_l2 in l2_types) {
  cat(sprintf("\n>> L2: %s\n", ct_l2))

  pb <- aggregate_pseudobulk(obj, ct_l2)
  if (is.null(pb)) { cat("  [SKIP] Insufficient pseudobulk\n"); next }

  pb_tissues <- unique(na.omit(pb$meta$tissue))
  if (length(pb_tissues) < 2) { cat("  [SKIP] < 2 tissues\n"); next }

  cat(sprintf("  Pseudobulk samples: %d (%s)\n",
              nrow(pb$meta), paste(pb_tissues, collapse = ", ")))

  pb_de_all[[ct_l2]]  <- list()
  enrich_all[[ct_l2]] <- list()
  agent_all[[ct_l2]]  <- list()

  for (pair in combn(as.character(pb_tissues), 2, simplify = FALSE)) {
    t1 <- pair[1]; t2 <- pair[2]
    comp_name <- paste0(t2, "_vs_", t1)
    cat(sprintf("  DESeq2: %s\n", comp_name))

    res <- run_deseq2_pairwise(pb, t1, t2)
    if (is.null(res)) { cat("    [SKIP] Too few samples per tissue\n"); next }

    pb_de_all[[ct_l2]][[comp_name]] <- res
    cat(sprintf("    DEGs: %d up, %d down (samples: %d vs %d)\n",
                res$n_up, res$n_down, res$n_samples_1, res$n_samples_2))

    # Save
    comp_dir <- file.path(de_dir, safe_name(ct_l2), safe_name(comp_name))
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(res$de_table, file.path(comp_dir, "DESeq2_results.csv"))

    # Volcano
    vp <- res$de_table %>%
      mutate(label = ifelse(sig == "sig" & rank(padj) <= 20, gene, ""))
    pv <- ggplot(vp, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("sig" = "red", "ns" = "grey70")) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
      geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", color = "blue") +
      labs(title = sprintf("DESeq2: %s (%s)", comp_name, ct_l2),
           x = "log2 Fold Change", y = "-log10(padj)") +
      theme_minimal()
    save_plot(pv, file.path(comp_dir, "volcano"), width = 8, height = 6)

    # --- Enrichment ---
    enrich_all[[ct_l2]][[comp_name]] <- list()
    agent_all[[ct_l2]][[comp_name]]  <- list()

    for (dir_name in c("up", "down")) {
      genes <- if (dir_name == "up") {
        res$de_table %>% filter(sig == "sig", log2FoldChange > 0) %>%
          arrange(desc(log2FoldChange)) %>% head(TOP_N_DEG_ENRICHMENT) %>% pull(gene)
      } else {
        res$de_table %>% filter(sig == "sig", log2FoldChange < 0) %>%
          arrange(log2FoldChange) %>% head(TOP_N_DEG_ENRICHMENT) %>% pull(gene)
      }

      if (length(genes) < 5) { enrich_all[[ct_l2]][[comp_name]][[dir_name]] <- NULL; next }
      cat(sprintf("    Enrichment [%s]: %d genes\n", dir_name, length(genes)))

      # P1-2 fix: pass DESeq2 tested genes as universe
      tested <- res$tested_genes
      enr_list <- list(
        GO_BP = run_gmt_enrichment(genes, go_bp_t2g, "GO_BP", tested),
        GO_MF = run_gmt_enrichment(genes, go_mf_t2g, "GO_MF", tested),
        GO_CC = run_gmt_enrichment(genes, go_cc_t2g, "GO_CC", tested),
        KEGG  = run_gmt_enrichment(genes, kegg_t2g,  "KEGG",  tested),
        Hallmark = run_gmt_enrichment(genes, hallmark_t2g, "Hallmark", tested)
      )
      enr_list <- enr_list[!sapply(enr_list, is.null)]
      enrich_all[[ct_l2]][[comp_name]][[dir_name]] <- enr_list

      # Save enrichment
      enr_out <- file.path(comp_dir, paste0("enrichment_", dir_name))
      dir.create(enr_out, showWarnings = FALSE)
      for (db in names(enr_list)) {
        er <- enr_list[[db]]
        if (nrow(as.data.frame(er)) > 0) {
          fwrite(as.data.frame(er), file.path(enr_out, paste0(db, ".csv")))
          tryCatch({
            pdf(file.path(enr_out, paste0(db, "_dotplot.pdf")), width = 10, height = 8)
            print(dotplot(er, showCategory = 15,
                          title = sprintf("%s %s (%s %s)", db, dir_name, ct_l2, comp_name)))
            dev.off()
          }, error = function(e) NULL)
        }
      }
      saveRDS(enr_list, file.path(enr_out, "all_enrichment.rds"))

      # --- interpret_agent (P1-5: only if LLM enabled) ---
      if (ENABLE_LLM) {
        best_er <- NULL
        for (db in c("GO_BP", "Hallmark", "KEGG", "GO_MF")) {
          if (!is.null(enr_list[[db]]) && nrow(as.data.frame(enr_list[[db]])) >= 3) {
            best_er <- enr_list[[db]]; break
          }
        }

        if (!is.null(best_er)) {
          tissue_ctx <- sprintf(
            "%s\nComparing %s vs %s for %s cells. Direction: %s-regulated in %s.",
            BCELL_CONTEXT, t2, t1, ct_l2, dir_name, t2
          )
          sig_de <- res$de_table %>% filter(sig == "sig")
          gene_fc <- setNames(sig_de$log2FoldChange, toupper(sig_de$gene))

          cat(sprintf("    interpret_agent [%s]\n", dir_name))
          ia <- run_interpret_agent_safe(best_er, tissue_ctx, gene_fc)
          agent_all[[ct_l2]][[comp_name]][[dir_name]] <- ia

          if (!is.null(ia)) {
            saveRDS(ia, file.path(comp_dir, paste0("interpret_agent_", dir_name, ".rds")))
            writeLines(capture.output(print(ia)),
                       file.path(comp_dir, paste0("interpret_agent_", dir_name, ".txt")))
          }
        }
      }  # end ENABLE_LLM
    }
  }
}

# ==============================================================================
# 9. Exploratory Wilcoxon (marker discovery only)
# ==============================================================================

cat("\n=== Exploratory Cell-level Wilcoxon ===\n")
cat("NOTE: P-values are inflated (pseudoreplication). For marker discovery ONLY.\n\n")

wilcox_all <- list()
for (ct_l2 in l2_types) {
  wx <- run_wilcox_exploratory(obj, ct_l2)
  if (!is.null(wx) && length(wx) > 0) {
    wilcox_all[[ct_l2]] <- wx
    for (comp_name in names(wx)) {
      wx_out <- file.path(wx_dir, safe_name(ct_l2))
      dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(wx[[comp_name]], file.path(wx_out, paste0(safe_name(comp_name), "_wilcox.csv")))
    }
    cat(sprintf("  [OK] %s: %d comparisons\n", ct_l2, length(wx)))
  }
}

# ==============================================================================
# 10. Save Summary RDS
# ==============================================================================

cat("\n=== Saving Summary ===\n")
saveRDS(pb_de_all,   file.path(rpt_dir, "pseudobulk_de_all.rds"))
saveRDS(wilcox_all,  file.path(rpt_dir, "wilcox_exploratory_all.rds"))
saveRDS(enrich_all,  file.path(rpt_dir, "enrichment_all.rds"))
saveRDS(agent_all,   file.path(rpt_dir, "interpret_agent_all.rds"))

# ==============================================================================
# 11. Generate REPORT.md
# ==============================================================================

cat("\n=== Generating REPORT.md ===\n")

md <- character()
add <- function(...) md <<- c(md, paste0(...))

add("# B Cell Tissue Comparison Report")
add("")
add("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M"))
add("")
add("**Pipeline:** B Cell Tissue Comparison v2.0 (pseudobulk DESeq2)")
add("")
add("---")
add("")

# Overview
add("## 1. Data Overview")
add("")
add(sprintf("- **Input:** `%s`", basename(H5AD_PATH)))
add(sprintf("- **Total cells:** %s", format(ncol(obj), big.mark = ",")))
add(sprintf("- **Tissues:** %s", paste(tissues, collapse = ", ")))
add(sprintf("- **L2 subtypes (%d):** %s", length(l2_types), paste(l2_types, collapse = ", ")))
if (l2_is_fallback) add("- **WARNING:** `cell_type_L2` missing; L3 used as proxy")
add("")

# Viz
add("## 2. Visualization")
add("")
add("### 2.1 UMAP")
add("![UMAP tissue](figures/umap_tissue.png)")
add("")
add("![UMAP celltype](figures/umap_celltype.png)")
add("")
add("![UMAP split](figures/umap_celltype_split_tissue.png)")
add("")
add("### 2.2 Marker Dotplot")
add("![Dotplot](figures/dotplot_markers.png)")
add("")
add("### 2.3 Cell Composition")
add("![Composition](figures/composition_tissue_celltype.png)")
add("")
add("![Counts](figures/count_tissue_celltype.png)")
add("")
add("### 2.4 Sample-level Composition")
add("![Sample composition](figures/composition_sample_level.png)")
add("")
add("### 2.5 Top Marker Heatmap")
add("![Heatmap](figures/heatmap_top_markers.png)")
add("")

# DESeq2 summary
add("## 3. Pseudobulk DESeq2 (Primary Inference)")
add("")
add(sprintf("Statistical unit: pseudobulk (sum counts per sample x tissue x L2). Thresholds: padj < %s, |log2FC| > %s.", PADJ_THR, LFC_THR))
add("")
add("| L2 Subtype | Comparison | Samples (ref/case) | Up | Down | Total |")
add("|---|---|---|---|---|---|")
for (ct_l2 in names(pb_de_all)) {
  for (comp_name in names(pb_de_all[[ct_l2]])) {
    r <- pb_de_all[[ct_l2]][[comp_name]]
    if (is.null(r)) next
    add(sprintf("| %s | %s | %d / %d | %d | %d | %d |",
                ct_l2, comp_name, r$n_samples_1, r$n_samples_2,
                r$n_up, r$n_down, r$n_up + r$n_down))
  }
}
add("")

# Enrichment
add("## 4. GO / KEGG / Hallmark Enrichment")
add("")
for (ct_l2 in names(enrich_all)) {
  for (comp_name in names(enrich_all[[ct_l2]])) {
    for (dir_name in names(enrich_all[[ct_l2]][[comp_name]])) {
      enr_l <- enrich_all[[ct_l2]][[comp_name]][[dir_name]]
      if (is.null(enr_l) || length(enr_l) == 0) next
      add(sprintf("### %s | %s | %s", ct_l2, comp_name, dir_name))
      add("")
      for (db in names(enr_l)) {
        er <- enr_l[[db]]
        if (is.null(er) || nrow(as.data.frame(er)) == 0) next
        top5 <- head(as.data.frame(er), 5)
        add(sprintf("**%s (top 5):**", db))
        add("")
        add("| Term | p.adjust | Count |")
        add("|---|---|---|")
        for (j in seq_len(nrow(top5))) {
          add(sprintf("| %s | %.2e | %s |", top5$Description[j], top5$p.adjust[j], top5$Count[j]))
        }
        add("")
      }
    }
  }
}

# LLM
add("## 5. LLM Interpretation (interpret_agent)")
add("")
if (!ENABLE_LLM) {
  add("**SKIPPED:** DEEPSEEK_API_KEY not set. Re-run with API key to enable.")
  add("")
} else {
  has_any_agent <- FALSE
  for (ct_l2 in names(agent_all)) {
    for (comp_name in names(agent_all[[ct_l2]])) {
      for (dir_name in names(agent_all[[ct_l2]][[comp_name]])) {
        ia <- agent_all[[ct_l2]][[comp_name]][[dir_name]]
        if (is.null(ia)) next
        has_any_agent <- TRUE
        add(sprintf("### %s | %s | %s", ct_l2, comp_name, dir_name))
        add("")
        res_one <- if (inherits(ia, "interpretation_list") && length(ia) >= 1) ia[[1]] else ia
        if (is.list(res_one)) {
          add(sprintf("- **Overview:** %s", to_scalar(res_one$overview)))
          add(sprintf("- **Key Mechanisms:** %s", to_scalar(res_one$key_mechanisms)))
          add(sprintf("- **Hypothesis:** %s", to_scalar(res_one$hypothesis)))
        }
        add("")
      }
    }
  }
  if (!has_any_agent) add("No interpret_agent results generated (insufficient enrichment).")
  add("")
}

# Methods
add("## 6. Methods")
add("")
add("- **DE:** Pseudobulk DESeq2 (sample-level aggregation; Squair et al. 2021 Nat Commun)")
add("- **Exploratory:** Cell-level Wilcoxon rank-sum (marker discovery only, NOT inference)")
add("- **Enrichment:** clusterProfiler::enricher() + MSigDB GMT (symbol-based, zero ID loss)")
add("- **Enrichment universe:** DESeq2-tested genes intersected with GMT (not full GMT)")
if (ENABLE_LLM) {
  add("- **LLM:** interpret_agent with DeepSeek deepseek-reasoner")
} else {
  add("- **LLM:** SKIPPED (no API key)")
}
if (!is.null(DESIGN_COVARIATE)) {
  add(sprintf("- **DESeq2 design:** ~ %s + tissue", DESIGN_COVARIATE))
} else {
  add("- **DESeq2 design:** ~ tissue")
}
add(sprintf("- **Context:** %s", BCELL_CONTEXT))
add("")
add("---")
add("*Generated by bcell_tissue_comparison_v2_1.R*")

writeLines(md, file.path(OUTPUT_DIR, "REPORT.md"))
cat("[OK] REPORT.md written\n")

# ==============================================================================
# 12. Final Summary
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("B CELL TISSUE COMPARISON COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Directory structure:\n")
cat("  REPORT.md                <- structured report (png embeds)\n")
cat("  figures/                 <- UMAP, dotplot, heatmap, composition (pdf+png)\n")
cat("  pseudobulk_de/           <- DESeq2 per L2 x tissue pair + volcano + enrichment\n")
cat("  wilcox_exploratory/      <- cell-level wilcox (marker discovery only)\n")
cat("  reports/                 <- RDS summary objects\n")

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
cat("\n[OK] Done\n")
