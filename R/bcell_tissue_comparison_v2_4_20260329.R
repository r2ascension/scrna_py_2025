#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison Pipeline v2.2.2 - PRODUCTION
# ==============================================================================
#
# Purpose:
#   B cell lineage tissue comparison across normal respiratory tract sites:
#     1. Load h5ad -> Seurat via SCNT::GetSeurat()
#     2. L3 -> L2 remapping (merge Memory/GC/Plasma subtypes)
#     3. Visualization (UMAP, marker dotplot/heatmap, composition)
#     4. Pseudobulk DESeq2 per L2 subtype x tissue pair
#     5. Exploratory cell-level wilcox (marker discovery only)
#     6. Multi-database enrichment (GO BP/MF/CC, KEGG, Hallmark,
#        CellMarker, PanglaoDB, custom B cell markers) - symbol direct
#     7. interpret_agent (DeepSeek) with tissue-pair-specific context
#     8. Structured REPORT.md (png embeds)
#
# v2.2.2 Changes (vs v2.2.1):
#   [LLM-1] build_tissue_pair_context(): Chinese output instruction + gene-pathway
#            emphasis appended to every LLM context string
#   [LLM-2] standardize_result_with_llm(): two new rules added to prompt:
#            (a) all narrative fields must be written in Chinese
#            (b) emphasize genes with largest absolute log2FC and their pathways
#
# v2.2.1 Fixes (vs v2.2):
#   [P1-1] UMAP output without raster / alpha; png dpi 300
#   [P1-2] aggregate_pseudobulk: sparse-friendly list-cbind (no dense pre-alloc)
#   [P1-3] Unmapped L3 values -> stop() fail-fast (not silent drop)
#   [P1-4] CellMarker/PanglaoDB: strict B cell / plasma / GC term filter
#   [P1-5] interpret_agent retry + atomic-result-tolerant normalization
#   [P1-6] Secondary DeepSeek standardization to fixed JSON schema
#   [P2-6] Heatmap title annotated "visualization only"
#   [P2-7] set.seed(42) inside stratified_downsample for cross-session reproducibility
#
# v2.2 Changes (vs v2.1):
#   [NEW] L3 -> L2 remapping via L3_TO_L2_REMAP (no more column candidates)
#   [NEW] Tissue-pair-specific LLM context (normal tissue comparison)
#   [NEW] Multi-database enrichment: CellMarker, PanglaoDB, custom B cell markers
#   [NEW] Custom B cell marker database (12 subtypes, 80+ markers)
#   [DEL] Removed all *_CANDIDATES lists and first_nonempty_column()
#
# Author: r2end
# Date:   2026-03-29
# ==============================================================================

# ==============================================================================
# 0. Thread Control
# ==============================================================================
library(CHOIR)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ==============================================================================
# 1. Configuration
# ==============================================================================

# ----- Input / Output -----
H5AD_PATH <- "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0329/bcell_tissue_comparison_v2_2_2"

# ----- L3 -> L2 Remapping -----
# Source column: cell_type_scanvi_pred (L3, 9 subtypes)
# Target column: cell_type_L2 (created by this pipeline)
L3_SOURCE_COL <- "cell_type_scanvi_pred"

L3_TO_L2_REMAP <- c(
  "Atypical_Memory_B" = "Memory_B",
  "IGHEplus_Atypical_Memory_B" = "Memory_B",
  "Memory_B" = "Memory_B",
  "GC_B_Dark_Zone_Centroblast_Cycling" = "GC_B",
  "GC_B_Light_Zone_Centrocyte" = "GC_B",
  "GC_B_Transitional" = "GC_B",
  "Plasma_IgA" = "Plasma",
  "Plasma_IgG" = "Plasma",
  "Naive_B" = "Naive_B"
)

# ----- Reference Database Paths -----
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"

# ----- DeepSeek API (OPTIONAL - pipeline runs without it) -----
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
ENABLE_LLM <- nchar(DEEPSEEK_API_KEY) >= 10
if (!ENABLE_LLM) {
  cat("[WARN] DEEPSEEK_API_KEY not set. interpret_agent will be SKIPPED.\n")
  cat("       Set via: export DEEPSEEK_API_KEY='your-key' to enable LLM.\n")
}
INTERPRET_AGENT_MODEL <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS <- 50
INTERPRET_AGENT_ADD_PPI <- TRUE
INTERPRET_AGENT_MAX_RETRIES <- 3L
INTERPRET_AGENT_RETRY_SLEEP_SEC <- 2
STANDARDIZE_LLM_OUTPUT <- TRUE
STANDARDIZE_LLM_MODEL <- "deepseek-chat"
STANDARDIZE_LLM_MAX_RETRIES <- 3L
STANDARDIZE_LLM_RETRY_SLEEP_SEC <- 2
STANDARDIZE_LLM_MAX_INPUT_CHARS <- 12000L

# ----- Python / reticulate -----
PYTHON_CONDA_ENV <- "bbknn_env"

# ----- Metadata Column Names (FIXED, no candidates) -----
TISSUE_COL <- "tissue"
SAMPLE_COL <- "sample"
CELLTYPE_L2_COL <- "cell_type_L2" # created by L3_TO_L2_REMAP
CELLTYPE_L3_COL <- "cell_type_scanvi_pred"
LABEL_COL <- "cell_type_scanvi_pred" # fine label for visualization

UMAP_REDUCTION_PREFERRED <- c(
  "umap_scanvi",
  "umap_scanvi_corrected",
  "umap_scvi",
  "umap"
)
UMAP_PT_SIZE <- 0.35
UMAP_TISSUE_COLORS <- c(
  "lung parenchyma" = "#D55E00",
  "nose" = "#009E73",
  "respiratory airway" = "#0072B2",
  "sinus" = "#CC79A7"
)

# ----- Pseudobulk DE Parameters -----
MIN_CELLS_PER_PSEUDOBULK <- 10
MIN_SAMPLES_PER_TISSUE <- 3
PADJ_THR <- 0.05
LFC_THR <- 1.0

DESIGN_COVARIATE <- NULL # e.g., "dataset" for ~ dataset + tissue

# ----- Reproducibility -----
set.seed(42)

# ----- Exploratory Wilcox Parameters -----
WILCOX_MIN_CELLS <- 50
WILCOX_LFC_THR <- 0.25
WILCOX_TOP_N <- 200

# ----- Enrichment -----
TOP_N_DEG_ENRICHMENT <- 200

# ----- Visualization -----
HEATMAP_CELLS_PER_TYPE <- 100

# ----- B Cell Known Markers -----
KNOWN_MARKERS <- c(
  # Pan-B
  "CD79A",
  "CD79B",
  "MS4A1",
  "CD19",
  "PAX5",
  # Naive
  "IGHD",
  "IGHM",
  "TCL1A",
  "FCER2",
  "IL4R",
  # Memory
  "CD27",
  "TNFRSF13B",
  "AIM2",
  # GC
  "BCL6",
  "AICDA",
  "RGS13",
  "MEF2B",
  "MME",
  "MKI67",
  # Plasma / Plasmablast
  "PRDM1",
  "XBP1",
  "MZB1",
  "JCHAIN",
  "SDC1",
  "DERL3",
  "SSR4",
  # Isotypes
  "IGHA1",
  "IGHG1",
  "IGHE",
  # Atypical / ABC
  "ITGAX",
  "TBX21",
  "FCRL4",
  "FCRL5",
  # Activation / IFN
  "CD69",
  "CD86",
  "ISG15",
  # CRSwNP-specific
  "LAPTM5",
  "CD74"
)

# ----- Custom B Cell Marker Database (for enrichment) -----
BCELL_MARKERS_DB <- data.frame(
  subtype = c(
    "Naive_B",
    "Transitional_B",
    "Memory_B_Unswitched",
    "Memory_B_Switched",
    "Atypical_Memory_B",
    "GC_B_Dark_Zone",
    "GC_B_Light_Zone",
    "Plasmablast",
    "Plasma_IgA",
    "Plasma_IgG",
    "Plasma_IgE",
    "Breg"
  ),
  markers = c(
    # Naive
    "IGHD,IGHM,TCL1A,FCER2,IL4R,CD38,CD24",
    # Transitional
    "CD24,CD38,IGHM,IGHD,MME,SOX4",
    # Memory unswitched
    "CD27,IGHM,IGHD,TNFRSF13B",
    # Memory switched
    "CD27,IGHA1,IGHG1,IGHG2,AIM2,TNFRSF13B",
    # Atypical / ABC
    "ITGAX,TBX21,FCRL4,FCRL5,ZEB2,CXCR3,FGR",
    # GC DZ
    "BCL6,AICDA,MKI67,TOP2A,CXCR4,FOXO1",
    # GC LZ
    "BCL6,LMO2,RGS13,MEF2B,CD83,CXCR5",
    # Plasmablast
    "PRDM1,XBP1,MZB1,JCHAIN,MKI67,IRF4",
    # Plasma IgA
    "PRDM1,XBP1,MZB1,JCHAIN,SDC1,IGHA1,IGHA2,DERL3,SSR4",
    # Plasma IgG
    "PRDM1,XBP1,MZB1,JCHAIN,SDC1,IGHG1,IGHG2,IGHG3,DERL3,SSR4",
    # Plasma IgE
    "PRDM1,XBP1,MZB1,JCHAIN,SDC1,IGHE,DERL3,SSR4",
    # Breg
    "IL10,CD24,CD27,GZMB,TGFB1"
  ),
  stringsAsFactors = FALSE
)

# ----- Tissue-Pair-Specific Context for LLM -----
# These are NORMAL tissues; no disease comparison
BCELL_BASE_CONTEXT <- paste(
  "B cells and plasma cells from NORMAL (non-diseased) human respiratory tract tissues.",
  "This is a cross-site anatomical comparison, NOT a disease vs healthy comparison.",
  "Key B cell subtypes: Naive B, Memory B (unswitched/switched/atypical),",
  "GC B (dark zone centroblast / light zone centrocyte / transitional),",
  "Plasma cell (IgA/IgG).",
  "Focus: regional variation in B cell composition, differentiation state,",
  "and mucosal humoral immunity along the respiratory tract."
)

# Tissue-specific biology for context enrichment
TISSUE_CONTEXT <- list(
  "nose" = paste(
    "Nasal cavity: first-line mucosal barrier, high antigen exposure,",
    "NALT (nasal-associated lymphoid tissue) supports local B cell responses,",
    "IgA-dominant secretory immunity, resident memory B cells."
  ),
  "sinus" = paste(
    "Paranasal sinus: semi-enclosed mucosal cavity, lower antigen load than nose,",
    "dependent on drainage for antigen clearance,",
    "potential site for ectopic lymphoid structures in inflammation,",
    "IgA secretion for mucosal defense."
  ),
  "respiratory airway" = paste(
    "Conducting airways (trachea/bronchi): ciliated epithelium with mucociliary clearance,",
    "BALT (bronchus-associated lymphoid tissue) in some individuals,",
    "IgA/IgG transport across epithelium, tissue-resident memory B cells,",
    "interface between upper and lower respiratory immunity."
  ),
  "lung parenchyma" = paste(
    "Lung parenchyma (alveolar region): gas exchange surface,",
    "thin epithelial barrier, alveolar macrophage-dominated immunity,",
    "lower B cell density than conducting airways,",
    "IgG-dominant (vs IgA in upper airways), systemic-like immune features,",
    "interstitial B cells and plasma cells near bronchovascular bundles."
  )
)

# ==============================================================================
# [LLM-1] build_tissue_pair_context
# Changes vs v2.2.1:
#   - Added Chinese output instruction (请用中文输出所有解释内容)
#   - Added gene-pathway emphasis instruction at the end of the context string
# ==============================================================================
build_tissue_pair_context <- function(t1, t2, ct_l2, dir_name) {
  ctx1 <- TISSUE_CONTEXT[[t1]]
  ctx2 <- TISSUE_CONTEXT[[t2]]
  if (is.null(ctx1)) {
    ctx1 <- paste(t1, ": no specific context available.")
  }
  if (is.null(ctx2)) {
    ctx2 <- paste(t2, ": no specific context available.")
  }

  sprintf(
    paste(
      "%s",
      "\nComparing %s vs %s for %s cells.",
      "Direction: genes %s-regulated in %s relative to %s.",
      "\n--- %s ---\n%s",
      "\n--- %s ---\n%s",
      "\nBiological question: What regional differences in B cell biology,",
      "mucosal immunity, or tissue microenvironment explain the observed",
      "transcriptional divergence between these two anatomical sites?",
      # [LLM-1] Chinese output instruction
      "\n\nOutput requirements:",
      "1. Write ALL interpretive content in Chinese (Chinese characters).",
      "2. For key_drivers and key_mechanisms, explicitly identify the genes with the",
      "   largest absolute log2FoldChange and explain WHICH enriched pathway(s) they",
      "   belong to and WHY those pathways are biologically meaningful in this tissue",
      "   comparison context.",
      "3. For each highlighted gene, state its log2FC value and the top 1-2 pathways",
      "   it contributes to, using the format: gene(log2FC=X.X) -> pathway name."
    ),
    BCELL_BASE_CONTEXT,
    t2,
    t1,
    ct_l2,
    dir_name,
    t2,
    t1,
    t2,
    ctx2,
    t1,
    ctx1
  )
}

# ==============================================================================
# 2. Load Libraries
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("B Cell Tissue Comparison Pipeline v2.2.2 (PRODUCTION)\n")
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
  py_bin <- tryCatch(py_config()$python, error = function(e) {
    paste0("conda:", PYTHON_CONDA_ENV)
  })
}

# Prefer local GetSeurat implementation
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
if (file.exists(LOCAL_GETSEURAT)) {
  cat(sprintf("[OK] Local GetSeurat loaded: %s\n", LOCAL_GETSEURAT))
}
cat("[OK] Libraries loaded\n\n")

# ==============================================================================
# 3. Load GMT + Reference Databases
# ==============================================================================

cat("=== Loading Reference Databases ===\n")

# ----- 3.1 MSigDB (Hallmark + KEGG) -----
gmt_all <- tryCatch(
  {
    g <- read.gmt(MSIGDB_GMT_PATH)
    g$gene <- toupper(g$gene)
    cat(sprintf("[OK] MSigDB: %d gene sets\n", length(unique(g$term))))
    g
  },
  error = function(e) {
    cat("[WARN] MSigDB failed\n")
    NULL
  }
)

hallmark_t2g <- if (!is.null(gmt_all)) {
  h <- gmt_all %>% filter(grepl("^HALLMARK_", term))
  if (nrow(h) > 0) {
    cat(sprintf("[OK] Hallmark: %d\n", length(unique(h$term))))
    h
  } else {
    NULL
  }
} else {
  NULL
}

kegg_t2g <- if (!is.null(gmt_all)) {
  k <- gmt_all %>% filter(grepl("^KEGG_", term))
  if (nrow(k) > 0) {
    cat(sprintf("[OK] KEGG: %d\n", length(unique(k$term))))
    k
  } else {
    NULL
  }
} else {
  NULL
}

# ----- 3.2 GO (BP/MF/CC) -----
go_bp_t2g <- NULL
go_mf_t2g <- NULL
go_cc_t2g <- NULL
tryCatch(
  {
    gmt_go <- read.gmt(GMT_GO_ALL)
    gmt_go$gene <- toupper(gmt_go$gene)
    go_bp_t2g <- gmt_go %>% filter(grepl("^GOBP_", term))
    go_mf_t2g <- gmt_go %>% filter(grepl("^GOMF_", term))
    go_cc_t2g <- gmt_go %>% filter(grepl("^GOCC_", term))
    cat(sprintf(
      "[OK] GO: BP=%d, MF=%d, CC=%d\n",
      length(unique(go_bp_t2g$term)),
      length(unique(go_mf_t2g$term)),
      length(unique(go_cc_t2g$term))
    ))
  },
  error = function(e) cat("[WARN] GO GMT failed\n")
)

# ----- 3.3 CellMarker (B cell / plasma / GC specific only) -----
cellmarker_t2g <- NULL
tryCatch(
  {
    cm_db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
    cm_db <- cm_db[grepl("Human", species, ignore.case = TRUE)]
    cat(sprintf("[OK] CellMarker raw: %d human entries\n", nrow(cm_db)))

    cm_bcell <- cm_db[
      grepl(
        "B cell|B-cell|B lymphocyte|Plasma cell|Plasmablast|Germinal center|Memory B|Naive B|Follicular B|Marginal zone B|Breg|Age-associated B|Atypical B|Pre-B|Pro-B|Transitional B",
        cell_name,
        ignore.case = TRUE
      )
    ]
    cat(sprintf(
      "[OK] CellMarker B-cell-specific: %d entries\n",
      nrow(cm_bcell)
    ))

    if (nrow(cm_bcell) < 10) {
      cat(
        "[WARN] Too few B cell entries in CellMarker; skipping this database\n"
      )
    } else {
      t2g_list <- list()
      for (i in seq_len(nrow(cm_bcell))) {
        markers_raw <- cm_bcell$marker[i]
        cell_type <- cm_bcell$cell_name[i]
        if (is.na(markers_raw) || markers_raw == "") {
          next
        }
        markers <- unlist(strsplit(markers_raw, "[,;\\s]+"))
        markers <- gsub('["\r\n\\[\\]]', '', markers)
        markers <- trimws(toupper(markers))
        markers <- unique(markers[markers != "" & !is.na(markers)])
        if (length(markers) > 0) {
          t2g_list[[length(t2g_list) + 1]] <- data.frame(
            term = cell_type,
            gene = markers,
            stringsAsFactors = FALSE
          )
        }
      }
      cellmarker_t2g <- dplyr::bind_rows(t2g_list) %>% distinct()
      cat(sprintf(
        "[OK] CellMarker TERM2GENE: %d pairs (%d B cell types)\n",
        nrow(cellmarker_t2g),
        length(unique(cellmarker_t2g$term))
      ))
    }
  },
  error = function(e) cat(sprintf("[WARN] CellMarker failed: %s\n", e$message))
)

# ----- 3.4 PanglaoDB (B cell / plasma specific only) -----
panglaodb_t2g <- NULL
tryCatch(
  {
    pdb <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)
    setnames(
      pdb,
      old = c("official gene symbol", "cell type"),
      new = c("gene_symbol", "cell_type"),
      skip_absent = TRUE
    )
    pdb <- pdb[grepl("Hs", species, fixed = TRUE)]
    cat(sprintf("[OK] PanglaoDB raw: %d human markers\n", nrow(pdb)))

    pdb_bcell <- pdb[
      grepl(
        "B cell|Plasma cell|Plasmablast|Germinal center|Memory B|Naive B|Follicular",
        cell_type,
        ignore.case = TRUE
      )
    ]
    cat(sprintf(
      "[OK] PanglaoDB B-cell-specific: %d markers (%d cell types)\n",
      nrow(pdb_bcell),
      length(unique(pdb_bcell$cell_type))
    ))

    if (nrow(pdb_bcell) < 5) {
      cat(
        "[WARN] Too few B cell markers in PanglaoDB; skipping this database\n"
      )
    } else {
      panglaodb_t2g <- pdb_bcell %>%
        dplyr::select(cell_type, gene_symbol) %>%
        mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
        filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
        distinct() %>%
        dplyr::rename(term = cell_type, gene = gene_symbol)
      cat(sprintf(
        "[OK] PanglaoDB TERM2GENE: %d pairs (%d B cell types)\n",
        nrow(panglaodb_t2g),
        length(unique(panglaodb_t2g$term))
      ))
    }
  },
  error = function(e) cat(sprintf("[WARN] PanglaoDB failed: %s\n", e$message))
)

# ----- 3.5 Custom B Cell Markers -----
bcell_custom_t2g <- BCELL_MARKERS_DB %>%
  tidyr::separate_rows(markers, sep = ",") %>%
  dplyr::mutate(markers = trimws(toupper(markers))) %>%
  dplyr::filter(markers != "") %>%
  dplyr::select(term = subtype, gene = markers) %>%
  distinct()
cat(sprintf(
  "[OK] Custom B cell markers: %d subtypes, %d pairs\n",
  length(unique(bcell_custom_t2g$term)),
  nrow(bcell_custom_t2g)
))

cat("\n")

# ==============================================================================
# 4. Helper Functions
# ==============================================================================

safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

pick_reduction <- function(obj, preferred) {
  red <- Reductions(obj)
  hit <- preferred[preferred %in% red]
  if (length(hit) > 0) hit[[1]] else NULL
}

build_umap_plot <- function(
  obj,
  reduction_name,
  group_col,
  title,
  split_col = NULL,
  label = FALSE,
  width = 10,
  height = 8,
  cols = NULL
) {
  p <- DimPlot(
    obj,
    reduction = reduction_name,
    group.by = group_col,
    split.by = split_col,
    pt.size = UMAP_PT_SIZE,
    shuffle = TRUE,
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
  if (is.null(t2g) || nrow(t2g) == 0 || length(gene_list) < 5) {
    return(NULL)
  }
  if (!is.null(tested_genes)) {
    universe_use <- intersect(toupper(tested_genes), unique(t2g$gene))
  } else {
    universe_use <- unique(t2g$gene)
  }
  tryCatch(
    enricher(
      gene = toupper(unique(gene_list)),
      TERM2GENE = t2g,
      universe = universe_use,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    ),
    error = function(e) {
      cat(sprintf("    [WARN] %s: %s\n", db_name, e$message))
      NULL
    }
  )
}

# ----- 4.2 Pseudobulk aggregation (sparse-friendly) -----
aggregate_pseudobulk <- function(obj, celltype_l2) {
  cells <- colnames(obj)[obj@meta.data[[CELLTYPE_L2_COL]] == celltype_l2]
  if (length(cells) < 30) {
    return(NULL)
  }

  sub <- subset(obj, cells = cells)
  meta_sub <- sub@meta.data
  counts_mat <- GetAssayData(sub, layer = "counts")

  meta_sub$pb_group <- paste0(
    meta_sub[[SAMPLE_COL]],
    "__",
    meta_sub[[TISSUE_COL]]
  )
  groups <- unique(meta_sub$pb_group)

  pb_cols <- list()
  pb_meta <- data.frame(row.names = character(), stringsAsFactors = FALSE)

  for (g in groups) {
    cells_g <- rownames(meta_sub)[meta_sub$pb_group == g]
    if (length(cells_g) < MIN_CELLS_PER_PSEUDOBULK) {
      next
    }

    col_sum <- if (length(cells_g) == 1) {
      counts_mat[, cells_g, drop = FALSE]
    } else {
      Matrix::rowSums(counts_mat[, cells_g, drop = FALSE])
    }

    if (is.numeric(col_sum) && is.null(dim(col_sum))) {
      col_sum <- Matrix::Matrix(col_sum, ncol = 1, sparse = TRUE)
      rownames(col_sum) <- rownames(counts_mat)
    }
    colnames(col_sum) <- g
    pb_cols[[g]] <- col_sum

    pb_meta[g, "tissue"] <- meta_sub[[TISSUE_COL]][meta_sub$pb_group == g][1]
    pb_meta[g, "sample"] <- meta_sub[[SAMPLE_COL]][meta_sub$pb_group == g][1]
    pb_meta[g, "n_cells"] <- length(cells_g)
  }

  if (length(pb_cols) < 4) {
    return(NULL)
  }

  pb_counts <- as.matrix(do.call(cbind, pb_cols))
  valid <- colnames(pb_counts)[
    colSums(pb_counts) > 0 & !is.na(pb_meta[colnames(pb_counts), "tissue"])
  ]
  if (length(valid) < 4) {
    return(NULL)
  }
  list(
    counts = pb_counts[, valid, drop = FALSE],
    meta = pb_meta[valid, , drop = FALSE]
  )
}

# ----- 4.3 DESeq2 pairwise -----
run_deseq2_pairwise <- function(pb, t1, t2) {
  keep <- pb$meta$tissue %in% c(t1, t2)
  if (sum(keep) < 4) {
    return(NULL)
  }

  counts_sub <- pb$counts[, keep, drop = FALSE]
  meta_sub <- pb$meta[keep, , drop = FALSE]
  tissue_levels_safe <- make.names(c(t1, t2), unique = TRUE)
  meta_sub$tissue <- droplevels(factor(meta_sub$tissue, levels = c(t1, t2)))
  meta_sub$tissue_safe <- factor(
    ifelse(meta_sub$tissue == t1, tissue_levels_safe[1], tissue_levels_safe[2]),
    levels = tissue_levels_safe
  )

  n1 <- sum(meta_sub$tissue == t1)
  n2 <- sum(meta_sub$tissue == t2)
  if (n1 < MIN_SAMPLES_PER_TISSUE || n2 < MIN_SAMPLES_PER_TISSUE) {
    return(NULL)
  }

  keep_genes <- rowSums(counts_sub >= 1) >= max(3, ncol(counts_sub) * 0.2)
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  if (nrow(counts_sub) < 100) {
    return(NULL)
  }

  tested_genes <- rownames(counts_sub)
  counts_int <- round(counts_sub)
  storage.mode(counts_int) <- "integer"

  design_formula <- ~tissue_safe
  if (!is.null(DESIGN_COVARIATE) && DESIGN_COVARIATE %in% colnames(meta_sub)) {
    n_levels <- length(unique(meta_sub[[DESIGN_COVARIATE]]))
    if (n_levels >= 2 && n_levels < nrow(meta_sub)) {
      meta_sub[[DESIGN_COVARIATE]] <- factor(meta_sub[[DESIGN_COVARIATE]])
      design_formula <- as.formula(paste(
        "~",
        DESIGN_COVARIATE,
        "+ tissue_safe"
      ))
      cat(sprintf(
        "    Design: %s (covariate: %d levels)\n",
        deparse(design_formula),
        n_levels
      ))
    }
  }

  dds <- tryCatch(
    {
      dds <- DESeqDataSetFromMatrix(counts_int, meta_sub, design_formula)
      DESeq(dds, quiet = TRUE)
    },
    error = function(e) {
      if (!is.null(DESIGN_COVARIATE)) {
        cat(sprintf(
          "    [WARN] Design failed (%s), falling back to ~ tissue_safe\n",
          e$message
        ))
        tryCatch(
          {
            dds2 <- DESeqDataSetFromMatrix(counts_int, meta_sub, ~tissue_safe)
            DESeq(dds2, quiet = TRUE)
          },
          error = function(e2) {
            cat(sprintf("    [ERROR] DESeq2: %s\n", e2$message))
            NULL
          }
        )
      } else {
        cat(sprintf("    [ERROR] DESeq2: %s\n", e$message))
        NULL
      }
    }
  )
  if (is.null(dds)) {
    return(NULL)
  }

  res <- results(
    dds,
    contrast = c("tissue_safe", tissue_levels_safe[2], tissue_levels_safe[1]),
    alpha = PADJ_THR
  )
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(
      sig = ifelse(
        padj < PADJ_THR & abs(log2FoldChange) > LFC_THR,
        "sig",
        "ns"
      ),
      direction = ifelse(log2FoldChange > 0, "up", "down")
    )

  list(
    de_table = res_df,
    tested_genes = tested_genes,
    n_up = sum(res_df$sig == "sig" & res_df$direction == "up"),
    n_down = sum(res_df$sig == "sig" & res_df$direction == "down"),
    n_samples_1 = n1,
    n_samples_2 = n2,
    tissue_1 = t1,
    tissue_2 = t2,
    design = deparse(design_formula)
  )
}

# ----- 4.4 Exploratory wilcox -----
run_wilcox_exploratory <- function(obj, celltype_l2) {
  cells <- colnames(obj)[obj@meta.data[[CELLTYPE_L2_COL]] == celltype_l2]
  if (length(cells) < WILCOX_MIN_CELLS) {
    return(NULL)
  }

  sub <- subset(obj, cells = cells)
  tissues <- sort(unique(na.omit(sub@meta.data[[TISSUE_COL]])))
  if (length(tissues) < 2) {
    return(NULL)
  }

  Idents(sub) <- TISSUE_COL
  results <- list()
  for (pair in combn(tissues, 2, simplify = FALSE)) {
    t1 <- pair[1]
    t2 <- pair[2]
    n1 <- sum(sub@meta.data[[TISSUE_COL]] == t1, na.rm = TRUE)
    n2 <- sum(sub@meta.data[[TISSUE_COL]] == t2, na.rm = TRUE)
    if (n1 < WILCOX_MIN_CELLS || n2 < WILCOX_MIN_CELLS) {
      next
    }

    de <- tryCatch(
      FindMarkers(
        sub,
        ident.1 = t2,
        ident.2 = t1,
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0,
        max.cells.per.ident = 5000
      ),
      error = function(e) NULL
    )
    if (!is.null(de) && nrow(de) > 0) {
      de$gene <- rownames(de)
      de$padj <- p.adjust(de$p_val, method = "BH")
      results[[paste0(t2, "_vs_", t1)]] <- de %>%
        filter(abs(avg_log2FC) >= WILCOX_LFC_THR) %>%
        arrange(padj) %>%
        head(WILCOX_TOP_N)
    }
  }
  results
}

# ----- 4.5 interpret_agent + standardize -----

safe_trim <- function(x) {
  x <- to_scalar(x)
  trimws(x)
}

capture_object_text <- function(x) {
  if (is.null(x)) {
    return("")
  }
  txt <- tryCatch(
    paste(capture.output(print(x)), collapse = "\n"),
    error = function(e) ""
  )
  if (!nzchar(trimws(txt))) {
    txt <- tryCatch(
      paste(capture.output(str(x, max.level = 3)), collapse = "\n"),
      error = function(e) ""
    )
  }
  trimws(txt)
}

unwrap_interpret_agent_result <- function(x) {
  if (inherits(x, "interpretation_list") && length(x) >= 1) {
    return(x[[1]])
  }
  if (
    is.list(x) && length(x) == 1 && (is.list(x[[1]]) || !is.null(names(x[[1]])))
  ) {
    return(x[[1]])
  }
  x
}

has_named_entries <- function(x) !is.null(names(x)) && any(nzchar(names(x)))

extract_named_text <- function(x, candidates) {
  for (nm in candidates) {
    val <- NULL
    if (is.list(x) && !is.null(x[[nm]])) {
      val <- x[[nm]]
    } else if (!is.list(x) && has_named_entries(x) && nm %in% names(x)) {
      val <- x[[nm]]
    }
    if (!is.null(val)) {
      txt <- safe_trim(val)
      if (nzchar(txt)) return(txt)
    }
  }
  ""
}

is_retryable_interpret_agent_issue <- function(
  warnings = character(),
  error_message = NULL
) {
  signals <- c(warnings, error_message)
  if (length(signals) == 0) {
    return(FALSE)
  }
  any(grepl(
    "Failed to parse JSON response|premature EOF|invalid for atomic vectors|429|5[0-9]{2}|timeout|temporar",
    signals,
    ignore.case = TRUE,
    perl = TRUE
  ))
}

looks_structured_interpret_agent_result <- function(x) {
  core <- unwrap_interpret_agent_result(x)
  any(nzchar(c(
    extract_named_text(
      core,
      c("overview", "summary", "interpretation", "narrative")
    ),
    extract_named_text(
      core,
      c("key_mechanisms", "mechanisms", "keyMechanisms")
    ),
    extract_named_text(core, c("hypothesis", "model", "working_hypothesis")),
    extract_named_text(
      core,
      c("key_drivers", "drivers", "genes", "gene_drivers")
    )
  )))
}

placeholder_text <- function(
  x,
  default = "Not available from current evidence."
) {
  x <- safe_trim(x)
  if (nzchar(x)) x else default
}

collapse_driver_field <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("")
  }
  if (is.list(x)) {
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  paste(unique(x), collapse = ", ")
}

truncate_text <- function(x, max_chars = STANDARDIZE_LLM_MAX_INPUT_CHARS) {
  x <- safe_trim(x)
  if (!nzchar(x)) {
    return("")
  }
  if (nchar(x, type = "chars") <= max_chars) {
    return(x)
  }
  paste0(substr(x, 1, max_chars), "\n...[truncated]")
}

build_enrichment_evidence_text <- function(
  enrich_obj,
  gene_fc = NULL,
  n_terms = 8L,
  n_genes = 15L
) {
  if (is.null(enrich_obj)) {
    return("")
  }
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return("")
  }

  top_terms <- utils::head(
    enrich_df[order(enrich_df$p.adjust), , drop = FALSE],
    n_terms
  )
  term_lines <- vapply(
    seq_len(nrow(top_terms)),
    function(i) {
      sprintf(
        "- %s | padj=%s | Count=%s | Genes=%s",
        safe_trim(top_terms$Description[i]),
        format(top_terms$p.adjust[i], scientific = TRUE, digits = 3),
        safe_trim(top_terms$Count[i]),
        truncate_text(gsub("/", ", ", safe_trim(top_terms$geneID[i])), 200)
      )
    },
    character(1)
  )

  gene_lines <- character()
  if (!is.null(gene_fc) && length(gene_fc) > 0) {
    gene_fc <- sort(gene_fc[!is.na(gene_fc)], decreasing = TRUE)
    top_up <- utils::head(gene_fc, n_genes)
    top_down <- utils::head(sort(gene_fc, decreasing = FALSE), n_genes)
    gene_lines <- c(
      sprintf(
        "Top up genes: %s",
        paste(sprintf("%s(%.2f)", names(top_up), top_up), collapse = ", ")
      ),
      sprintf(
        "Top down genes: %s",
        paste(sprintf("%s(%.2f)", names(top_down), top_down), collapse = ", ")
      )
    )
  }

  paste(c("Top enrichment evidence:", term_lines, gene_lines), collapse = "\n")
}

build_multi_enrichment_evidence_text <- function(
  enrich_list,
  gene_fc = NULL,
  db_order = c(
    "GO_BP",
    "GO_MF",
    "GO_CC",
    "KEGG",
    "Hallmark",
    "CellMarker",
    "PanglaoDB",
    "B_cell_custom"
  ),
  n_terms_per_db = 5L,
  n_genes = 15L
) {
  if (is.null(enrich_list) || length(enrich_list) == 0) {
    return("")
  }

  dbs <- intersect(db_order, names(enrich_list))
  if (length(dbs) == 0) {
    dbs <- names(enrich_list)
  }

  blocks <- list()
  for (db in dbs) {
    enrich_obj <- enrich_list[[db]]
    if (is.null(enrich_obj)) {
      next
    }
    enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
    if (is.null(enrich_df) || nrow(enrich_df) == 0) {
      next
    }
    blocks[[length(blocks) + 1L]] <- paste(
      sprintf("[%s]", db),
      build_enrichment_evidence_text(
        enrich_obj,
        gene_fc = gene_fc,
        n_terms = n_terms_per_db,
        n_genes = n_genes
      ),
      sep = "\n"
    )
  }

  if (length(blocks) == 0) {
    return("")
  }

  paste(blocks, collapse = "\n\n")
}

extract_json_string <- function(text) {
  text <- safe_trim(text)
  if (!nzchar(text)) {
    return("")
  }
  text <- gsub("^```(?:json)?\\s*", "", text, perl = TRUE)
  text <- gsub("\\s*```$", "", text, perl = TRUE)
  start_idx <- regexpr("\\{", text, perl = TRUE)[1]
  end_positions <- gregexpr("\\}", text, perl = TRUE)[[1]]
  if (start_idx < 1 || length(end_positions) == 0 || end_positions[1] < 1) {
    return("")
  }
  end_idx <- end_positions[length(end_positions)]
  substr(text, start_idx, end_idx)
}

parse_standardized_json <- function(text) {
  candidates <- unique(c(safe_trim(text), extract_json_string(text)))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) == 0) {
    return(NULL)
  }
  for (candidate in candidates) {
    parsed <- tryCatch(
      jsonlite::fromJSON(candidate, simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (!is.null(parsed) && is.list(parsed)) return(parsed)
  }
  NULL
}

# ==============================================================================
# [LLM-2] standardize_result_with_llm
# Changes vs v2.2.1:
#   Rule 6 (NEW): All narrative fields must be written in Chinese
#   Rule 7 (NEW): Explicitly highlight top-|log2FC| genes and their pathways
# ==============================================================================
standardize_result_with_llm <- function(
  raw_text,
  context_str,
  evidence_text,
  ct_l2,
  comp_name,
  dir_name,
  source_db,
  warnings = character(),
  error_message = NULL
) {
  if (!ENABLE_LLM || !STANDARDIZE_LLM_OUTPUT) {
    return(list(
      result = NULL,
      warnings = character(),
      error = "standardization disabled"
    ))
  }

  prompt <- paste(
    "You are standardizing a biological interpretation for single-cell B-cell tissue comparison.",
    "Return valid JSON only. No markdown, no code fences, no extra commentary.",
    "Use exactly these keys:",
    "overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "",
    "Rules:",
    "1. overview/key_mechanisms/hypothesis/narrative/evidence/limitations must be strings.",
    "2. key_drivers must be an array of short gene or regulator names.",
    "3. Every field must be present.",
    "4. If evidence is weak, keep the content conservative and say so in limitations.",
    "5. If a field cannot be recovered, use 'Not available from current evidence.'.",
    # [LLM-2a] Chinese output rule
    "6. Write ALL narrative content (overview, key_mechanisms, hypothesis, narrative,",
    "   evidence, limitations) in Chinese (Simplified Chinese characters).",
    "   key_drivers should remain as gene symbol strings (e.g. [\"BCL6\", \"AICDA\"]).",
    # [LLM-2b] Gene-pathway emphasis rule
    "7. In key_mechanisms and key_drivers, prioritize genes with the largest absolute",
    "   log2FoldChange values. For each of the top 5 such genes, explicitly state:",
    "   (a) the gene name and its log2FC value,",
    "   (b) which enriched pathway(s) it contributes to (from the enrichment evidence),",
    "   (c) why this gene-pathway relationship is biologically meaningful for this",
    "       specific tissue comparison.",
    "   Use the format in key_mechanisms: gene(log2FC=X.X) -> pathway: explanation.",
    "",
    sprintf("Cell type L2: %s", ct_l2),
    sprintf("Comparison: %s", comp_name),
    sprintf("Direction: %s", dir_name),
    sprintf("Source DB: %s", source_db),
    "Context:",
    truncate_text(context_str),
    "Original interpret_agent output:",
    truncate_text(raw_text),
    "Enrichment evidence:",
    truncate_text(evidence_text),
    "Warnings and errors:",
    truncate_text(paste(c(warnings, error_message), collapse = "\n")),
    sep = "\n"
  )

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(STANDARDIZE_LLM_MAX_RETRIES)) {
    if (attempt > 1) {
      cat(sprintf(
        "    [INFO] standardize_llm retry %d/%d\n",
        attempt,
        STANDARDIZE_LLM_MAX_RETRIES
      ))
    } else {
      cat(sprintf("    [INFO] standardize_llm (%s)\n", STANDARDIZE_LLM_MODEL))
    }

    response_text <- tryCatch(
      safe_trim(fanyi::chat_request(
        prompt,
        model = STANDARDIZE_LLM_MODEL,
        api_key = DEEPSEEK_API_KEY
      )),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )

    parsed <- parse_standardized_json(response_text)
    if (!is.null(parsed)) {
      return(list(
        result = parsed,
        warnings = unique(collected_warnings),
        error = NULL
      ))
    }

    if (nzchar(response_text)) {
      collected_warnings <- c(
        collected_warnings,
        sprintf("standardize attempt %d returned non-JSON text", attempt)
      )
    }

    if (attempt < STANDARDIZE_LLM_MAX_RETRIES) {
      Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
    }
  }

  list(
    result = NULL,
    warnings = unique(collected_warnings),
    error = ifelse(
      is.null(last_error) || !nzchar(last_error),
      "failed to standardize interpret_agent output",
      last_error
    )
  )
}

normalize_interpret_agent_result <- function(
  raw_result,
  ct_l2,
  comp_name,
  dir_name,
  source_db = NA_character_,
  warnings = character(),
  error_message = NULL,
  enrich_obj = NULL,
  enrich_list = NULL,
  context_str = NULL,
  gene_fc = NULL,
  force_standardization = FALSE
) {
  core <- unwrap_interpret_agent_result(raw_result)
  raw_text <- capture_object_text(raw_result)

  overview <- extract_named_text(
    core,
    c("overview", "summary", "interpretation", "narrative")
  )
  key_mechanisms <- extract_named_text(
    core,
    c("key_mechanisms", "mechanisms", "keyMechanisms")
  )
  hypothesis <- extract_named_text(
    core,
    c("hypothesis", "model", "working_hypothesis")
  )
  narrative <- extract_named_text(core, c("narrative", "story", "details"))
  key_drivers <- extract_named_text(
    core,
    c("key_drivers", "drivers", "genes", "gene_drivers")
  )
  evidence <- extract_named_text(
    core,
    c("evidence", "supporting_evidence", "rationale")
  )
  limitations <- extract_named_text(core, c("limitations", "caveats", "notes"))

  if (!nzchar(overview) && nzchar(raw_text)) {
    overview <- raw_text
  }

  major_fields_present <- c(
    overview,
    key_mechanisms,
    hypothesis,
    narrative,
    evidence
  )
  needs_standardization <- isTRUE(STANDARDIZE_LLM_OUTPUT) &&
    (isTRUE(force_standardization) ||
      !all(nzchar(major_fields_present)) ||
      !nzchar(key_drivers) ||
      !is.null(error_message) ||
      !is.list(core))

  if (needs_standardization) {
    evidence_text <- if (!is.null(enrich_list) && length(enrich_list) > 0) {
      build_multi_enrichment_evidence_text(
        enrich_list,
        gene_fc = gene_fc
      )
    } else {
      build_enrichment_evidence_text(
        enrich_obj,
        gene_fc = gene_fc
      )
    }
    standardized_payload <- standardize_result_with_llm(
      raw_text = raw_text,
      context_str = context_str,
      evidence_text = evidence_text,
      ct_l2 = ct_l2,
      comp_name = comp_name,
      dir_name = dir_name,
      source_db = source_db,
      warnings = warnings,
      error_message = error_message
    )

    if (!is.null(standardized_payload$result)) {
      std <- standardized_payload$result
      overview <- placeholder_text(std$overview)
      key_mechanisms <- placeholder_text(std$key_mechanisms)
      hypothesis <- placeholder_text(std$hypothesis)
      narrative <- placeholder_text(std$narrative)
      key_drivers <- placeholder_text(collapse_driver_field(std$key_drivers))
      evidence <- placeholder_text(std$evidence)
      limitations <- placeholder_text(std$limitations)
      warnings <- unique(c(
        warnings,
        standardized_payload$warnings,
        if (!is.null(error_message) && nzchar(error_message)) {
          sprintf("initial interpret_agent error: %s", error_message)
        }
      ))
      error_message <- NULL
    } else {
      warnings <- unique(c(warnings, standardized_payload$warnings))
      limitations <- placeholder_text(limitations)
      if (
        !is.null(standardized_payload$error) &&
          nzchar(standardized_payload$error)
      ) {
        error_message <- paste(
          unique(c(error_message, standardized_payload$error)),
          collapse = " | "
        )
      }
    }
  }

  overview <- placeholder_text(overview)
  key_mechanisms <- placeholder_text(key_mechanisms)
  hypothesis <- placeholder_text(hypothesis)
  narrative <- placeholder_text(narrative)
  key_drivers <- placeholder_text(key_drivers)
  evidence <- placeholder_text(evidence)
  limitations <- placeholder_text(limitations)

  status <- dplyr::case_when(
    !is.null(error_message) ~ "error",
    is.null(raw_result) ~ "missing",
    any(nzchar(c(
      overview,
      key_mechanisms,
      hypothesis,
      narrative,
      key_drivers,
      evidence,
      limitations
    ))) ~ "structured",
    nzchar(raw_text) ~ "raw_text",
    TRUE ~ "empty"
  )

  list(
    celltype_l2 = ct_l2,
    comparison = comp_name,
    direction = dir_name,
    source_db = source_db,
    status = status,
    warnings = unique(warnings),
    error = if (is.null(error_message)) "" else error_message,
    overview = overview,
    key_mechanisms = key_mechanisms,
    hypothesis = hypothesis,
    narrative = narrative,
    key_drivers = key_drivers,
    evidence = evidence,
    limitations = limitations,
    raw_text = raw_text,
    raw_result = raw_result
  )
}

flatten_interpretation_records <- function(x) {
  rows <- list()
  idx <- 1L
  for (ct_l2 in names(x)) {
    for (comp_name in names(x[[ct_l2]])) {
      for (dir_name in names(x[[ct_l2]][[comp_name]])) {
        rec <- x[[ct_l2]][[comp_name]][[dir_name]]
        if (is.null(rec)) {
          next
        }
        rows[[idx]] <- data.frame(
          celltype_l2 = safe_trim(rec$celltype_l2),
          comparison = safe_trim(rec$comparison),
          direction = safe_trim(rec$direction),
          source_db = safe_trim(rec$source_db),
          status = safe_trim(rec$status),
          warnings = paste(rec$warnings, collapse = " | "),
          error = safe_trim(rec$error),
          overview = safe_trim(rec$overview),
          key_mechanisms = safe_trim(rec$key_mechanisms),
          hypothesis = safe_trim(rec$hypothesis),
          narrative = safe_trim(rec$narrative),
          key_drivers = safe_trim(rec$key_drivers),
          evidence = safe_trim(rec$evidence),
          limitations = safe_trim(rec$limitations),
          raw_text = safe_trim(rec$raw_text),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(rows)
}

write_interpretation_markdown <- function(
  records,
  path,
  title,
  include_raw = FALSE
) {
  md <- c(title, "")
  md <- c(
    md,
    sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
    ""
  )
  if (nrow(records) == 0) {
    md <- c(md, "No interpret_agent records available.")
    writeLines(md, path)
    return(invisible(path))
  }
  for (i in seq_len(nrow(records))) {
    rec <- records[i, , drop = FALSE]
    md <- c(
      md,
      sprintf(
        "## %s | %s | %s",
        rec$celltype_l2,
        rec$comparison,
        rec$direction
      ),
      ""
    )
    md <- c(md, sprintf("- **Status:** %s", rec$status))
    md <- c(
      md,
      sprintf(
        "- **Source DB:** %s",
        ifelse(nzchar(rec$source_db), rec$source_db, "NA")
      )
    )
    if (nzchar(rec$warnings)) {
      md <- c(md, sprintf("- **Warnings:** %s", rec$warnings))
    }
    if (nzchar(rec$error)) {
      md <- c(md, sprintf("- **Error:** %s", rec$error))
    }
    md <- c(md, "")
    field_block <- list(
      Overview = rec$overview,
      `Key Mechanisms` = rec$key_mechanisms,
      Hypothesis = rec$hypothesis,
      Narrative = rec$narrative,
      `Key Drivers` = rec$key_drivers,
      Evidence = rec$evidence,
      Limitations = rec$limitations
    )
    for (nm in names(field_block)) {
      val <- placeholder_text(field_block[[nm]])
      md <- c(md, sprintf("### %s", nm), "", val, "")
    }
    if (include_raw && nzchar(rec$raw_text)) {
      md <- c(md, "### Raw Output", "", rec$raw_text, "")
    }
  }
  writeLines(md, path)
  invisible(path)
}

run_interpret_agent_once <- function(enrich_obj, context_str, gene_fc = NULL) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(list(
      result = NULL,
      warnings = character(),
      error = "empty enrichment"
    ))
  }
  warn_msgs <- character()
  res <- tryCatch(
    withCallingHandlers(
      clusterProfiler::interpret_agent(
        x = enrich_obj,
        context = context_str,
        n_pathways = INTERPRET_AGENT_N_PATHWAYS,
        model = INTERPRET_AGENT_MODEL,
        api_key = DEEPSEEK_API_KEY,
        add_ppi = INTERPRET_AGENT_ADD_PPI,
        gene_fold_change = gene_fc
      ),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat(sprintf("    [ERROR] agent: %s\n", e$message))
      structure(list(message = e$message), class = "interpret_agent_error")
    }
  )
  if (inherits(res, "interpret_agent_error")) {
    return(list(
      result = NULL,
      warnings = unique(warn_msgs),
      error = res$message
    ))
  }
  list(result = res, warnings = unique(warn_msgs), error = NULL)
}

run_interpret_agent_safe <- function(enrich_obj, context_str, gene_fc = NULL) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(list(
      result = NULL,
      warnings = character(),
      error = "empty enrichment"
    ))
  }

  combined_warnings <- character()
  last_payload <- list(
    result = NULL,
    warnings = character(),
    error = "interpret_agent did not run"
  )

  for (attempt in seq_len(INTERPRET_AGENT_MAX_RETRIES)) {
    if (attempt > 1) {
      cat(sprintf(
        "    [INFO] interpret_agent retry %d/%d\n",
        attempt,
        INTERPRET_AGENT_MAX_RETRIES
      ))
    }
    payload <- run_interpret_agent_once(enrich_obj, context_str, gene_fc)
    last_payload <- payload

    if (length(payload$warnings) > 0) {
      combined_warnings <- c(
        combined_warnings,
        sprintf("attempt %d: %s", attempt, payload$warnings)
      )
    }

    retryable <- is_retryable_interpret_agent_issue(
      payload$warnings,
      payload$error
    )
    structured <- !is.null(payload$result) &&
      looks_structured_interpret_agent_result(payload$result)
    raw_text_valid <- !is.null(payload$result) &&
      nzchar(capture_object_text(payload$result))

    if (
      is.null(payload$error) && (structured || (raw_text_valid && !retryable))
    ) {
      payload$warnings <- unique(c(combined_warnings, payload$warnings))
      return(payload)
    }
    if (
      attempt < INTERPRET_AGENT_MAX_RETRIES &&
        (retryable || is.null(payload$result))
    ) {
      Sys.sleep(INTERPRET_AGENT_RETRY_SLEEP_SEC)
      next
    }
    break
  }

  last_payload$warnings <- unique(c(combined_warnings, last_payload$warnings))
  if (is.null(last_payload$error) && !is.null(last_payload$result)) {
    return(last_payload)
  }
  if (is.null(last_payload$error) || !nzchar(last_payload$error)) {
    last_payload$error <- sprintf(
      "interpret_agent failed after %d attempts",
      INTERPRET_AGENT_MAX_RETRIES
    )
  }
  last_payload
}

to_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("")
  }
  if (is.character(x) && length(x) == 1) {
    return(x)
  }
  if (is.character(x)) {
    return(paste(x, collapse = "; "))
  }
  if (is.list(x)) {
    return(paste(capture.output(str(x, max.level = 2)), collapse = "\n"))
  }
  as.character(x)
}

# ----- 4.6 Stratified downsampling (visualization only) -----
stratified_downsample <- function(
  obj,
  group_col,
  n_per = HEATMAP_CELLS_PER_TYPE
) {
  set.seed(42)
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
  ggsave(
    paste0(path_no_ext, ".png"),
    p,
    width = width,
    height = height,
    dpi = 300
  )
}

# ==============================================================================
# 5. Load Data
# ==============================================================================

cat("=== Loading B Cell Data ===\n")

if (!file.exists(H5AD_PATH)) {
  stop(sprintf("File not found: %s", H5AD_PATH))
}
obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
cat(sprintf("[OK] %d cells x %d genes\n", ncol(obj), nrow(obj)))

if (!"counts" %in% Layers(obj[["RNA"]])) {
  stop("RNA assay missing 'counts' layer. Check GetSeurat output.")
}
cat("[OK] counts layer verified\n\n")

fig_dir <- file.path(OUTPUT_DIR, "figures")
rpt_dir <- file.path(OUTPUT_DIR, "reports")
de_dir <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir <- file.path(OUTPUT_DIR, "wilcox_exploratory")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rpt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wx_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. Validate Metadata & Apply L2 Remapping
# ==============================================================================

cat("=== Validating Metadata & L2 Remapping ===\n")

meta <- obj@meta.data

for (col in c(TISSUE_COL, SAMPLE_COL, L3_SOURCE_COL)) {
  if (!col %in% colnames(meta)) {
    stop(sprintf("Missing required column: %s", col))
  }
}

l3_vals <- as.character(meta[[L3_SOURCE_COL]])
l2_vals <- L3_TO_L2_REMAP[l3_vals]

unmapped <- unique(l3_vals[is.na(l2_vals)])
if (length(unmapped) > 0) {
  stop(sprintf(
    "Unmapped L3 values found in '%s': %s\nUpdate L3_TO_L2_REMAP before running.",
    L3_SOURCE_COL,
    paste(unmapped, collapse = ", ")
  ))
}

obj@meta.data[[CELLTYPE_L2_COL]] <- l2_vals

cat("\nL3 -> L2 Remapping Summary:\n")
print(table(L3 = l3_vals, L2 = l2_vals, useNA = "ifany"))
cat("\nL2 Distribution:\n")
print(table(obj@meta.data[[CELLTYPE_L2_COL]], useNA = "ifany"))
cat("\n")

bad_idx <- is.na(meta[[TISSUE_COL]]) |
  trimws(as.character(meta[[TISSUE_COL]])) == "" |
  is.na(meta[[SAMPLE_COL]]) |
  trimws(as.character(meta[[SAMPLE_COL]])) == "" |
  is.na(obj@meta.data[[CELLTYPE_L2_COL]]) |
  is.na(meta[[LABEL_COL]]) |
  trimws(as.character(meta[[LABEL_COL]])) == ""
if (sum(bad_idx) > 0) {
  cat(sprintf(
    "[INFO] Dropping %d cells with NA tissue/sample/L2/label\n",
    sum(bad_idx)
  ))
  obj <- subset(obj, cells = colnames(obj)[!bad_idx])
}

n_samples <- length(unique(obj@meta.data[[SAMPLE_COL]]))
if (n_samples < 2) {
  stop(sprintf(
    "Only %d unique sample(s). Need >= 2 for pseudobulk DE.",
    n_samples
  ))
}
cat(sprintf("[OK] %d unique samples\n", n_samples))

meta <- obj@meta.data
tissues <- sort(unique(na.omit(meta[[TISSUE_COL]])))
l2_types <- sort(unique(na.omit(meta[[CELLTYPE_L2_COL]])))

cat(sprintf("[OK] Cells: %d\n", ncol(obj)))
cat(sprintf("[OK] Tissues: %s\n", paste(tissues, collapse = ", ")))
cat(sprintf(
  "[OK] L2 types (%d): %s\n",
  length(l2_types),
  paste(l2_types, collapse = ", ")
))
cat(sprintf(
  "[OK] L3 column (label): %s (%d unique)\n",
  LABEL_COL,
  length(unique(meta[[LABEL_COL]]))
))
cat("\n")

cat("[INFO] Running NormalizeData() for visualization and marker analyses\n")
obj <- NormalizeData(obj, verbose = FALSE)
cat("[OK] data layer ready\n\n")

# ==============================================================================
# 7. Visualization
# ==============================================================================

cat("=== Visualization ===\n")

umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
if (!is.null(umap_reduction)) {
  tissue_cols_use <- UMAP_TISSUE_COLORS[names(UMAP_TISSUE_COLORS) %in% tissues]

  p1 <- build_umap_plot(
    obj,
    umap_reduction,
    TISSUE_COL,
    title = sprintf("B Cell - Tissue (%s)", umap_reduction),
    cols = tissue_cols_use,
    width = 11,
    height = 8
  )
  save_plot(
    p1$plot,
    file.path(fig_dir, "umap_tissue"),
    width = p1$width,
    height = p1$height
  )

  p2 <- build_umap_plot(
    obj,
    umap_reduction,
    CELLTYPE_L2_COL,
    title = sprintf("B Cell - Cell Type L2 (%s)", umap_reduction),
    label = TRUE,
    width = 14,
    height = 10
  )
  save_plot(
    p2$plot,
    file.path(fig_dir, "umap_celltype_L2"),
    width = p2$width,
    height = p2$height
  )

  p2s <- build_umap_plot(
    obj,
    umap_reduction,
    CELLTYPE_L2_COL,
    title = sprintf("B Cell - L2 by Tissue (%s)", umap_reduction),
    split_col = TISSUE_COL,
    label = TRUE,
    width = max(12, 4.5 * length(tissues)),
    height = 8
  )
  save_plot(
    p2s$plot,
    file.path(fig_dir, "umap_L2_split_tissue"),
    width = p2s$width,
    height = p2s$height
  )

  p2_l3 <- build_umap_plot(
    obj,
    umap_reduction,
    LABEL_COL,
    title = sprintf("B Cell - Cell Type L3 (%s)", umap_reduction),
    label = TRUE,
    width = 18,
    height = 12
  )
  save_plot(
    p2_l3$plot,
    file.path(fig_dir, "umap_celltype_L3"),
    width = p2_l3$width,
    height = p2_l3$height
  )

  cat(sprintf("[OK] UMAP saved using reduction: %s\n", umap_reduction))
} else {
  cat("[WARN] No UMAP reduction\n")
}

markers_present <- intersect(KNOWN_MARKERS, rownames(obj))
if (length(markers_present) >= 3) {
  Idents(obj) <- LABEL_COL
  p3 <- DotPlot(obj, features = markers_present) +
    RotatedAxis() +
    ggtitle("B Cell - Known Markers (L3)") +
    theme(axis.text.x = element_text(size = 7))
  save_plot(
    p3,
    file.path(fig_dir, "dotplot_markers"),
    width = max(12, length(markers_present) * 0.45),
    height = max(6, length(unique(meta[[LABEL_COL]])) * 0.4)
  )
  cat("[OK] Dotplot saved\n")
}

comp_df <- meta %>%
  filter(!is.na(!!sym(TISSUE_COL)), !is.na(!!sym(CELLTYPE_L2_COL))) %>%
  count(!!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(TISSUE_COL)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p4 <- ggplot(
  comp_df,
  aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(CELLTYPE_L2_COL))
) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Tissue",
    y = "Percentage (%)",
    title = "B Cell - L2 Composition by Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4, file.path(fig_dir, "composition_tissue_L2"))

p4b <- ggplot(
  comp_df,
  aes(x = !!sym(TISSUE_COL), y = n, fill = !!sym(CELLTYPE_L2_COL))
) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Tissue",
    y = "Cell Count",
    title = "B Cell - Absolute Count by Tissue"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4b, file.path(fig_dir, "count_tissue_L2"))
cat("[OK] Composition plots saved\n")

sample_comp <- meta %>%
  filter(
    !is.na(!!sym(TISSUE_COL)),
    !is.na(!!sym(CELLTYPE_L2_COL)),
    !is.na(!!sym(SAMPLE_COL))
  ) %>%
  count(!!sym(SAMPLE_COL), !!sym(TISSUE_COL), !!sym(CELLTYPE_L2_COL)) %>%
  group_by(!!sym(SAMPLE_COL)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

if (nrow(sample_comp) > 0) {
  facet_ncol <- min(4, length(l2_types))
  facet_nrow <- ceiling(length(l2_types) / facet_ncol)
  p4c <- ggplot(
    sample_comp,
    aes(x = !!sym(TISSUE_COL), y = pct, fill = !!sym(TISSUE_COL))
  ) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    facet_wrap(
      as.formula(paste("~", CELLTYPE_L2_COL)),
      scales = "free_y",
      ncol = facet_ncol
    ) +
    labs(
      x = "Tissue",
      y = "Proportion per Sample (%)",
      title = "B Cell - Sample-level Composition by Tissue"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  save_plot(
    p4c,
    file.path(fig_dir, "composition_sample_level"),
    width = max(10, facet_ncol * 4),
    height = max(6, facet_nrow * 3.5)
  )
  cat("[OK] Sample-level composition saved\n")
}

Idents(obj) <- CELLTYPE_L2_COL
top_mk <- tryCatch(
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    test.use = "wilcox"
  ),
  error = function(e) {
    cat("[WARN] FindAllMarkers failed\n")
    NULL
  }
)

if (!is.null(top_mk) && nrow(top_mk) > 0) {
  fwrite(top_mk, file.path(rpt_dir, "all_markers_per_L2.csv"))
  top10 <- top_mk %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10) %>%
    ungroup()
  obj_ds <- stratified_downsample(obj, CELLTYPE_L2_COL, HEATMAP_CELLS_PER_TYPE)
  obj_ds <- ScaleData(obj_ds, features = unique(top10$gene), verbose = FALSE)
  p5 <- DoHeatmap(obj_ds, features = unique(top10$gene), size = 3) +
    ggtitle("B Cell - Top Markers per L2 (downsampled, visualization only)")
  save_plot(
    p5,
    file.path(fig_dir, "heatmap_top_markers"),
    width = 14,
    height = 10
  )
  rm(obj_ds)
  gc()
  cat("[OK] Heatmap saved\n")
}
cat("\n")

# ==============================================================================
# 8. Pseudobulk DESeq2 + Multi-Database Enrichment + interpret_agent
# ==============================================================================

cat("=== Pseudobulk DESeq2 Tissue Comparison ===\n")
cat(sprintf("Thresholds: padj < %s, |log2FC| > %s\n", PADJ_THR, LFC_THR))
cat(
  "Enrichment databases: GO_BP, GO_MF, GO_CC, KEGG, Hallmark, CellMarker, PanglaoDB, B_cell_custom\n\n"
)

pb_de_all <- list()
enrich_all <- list()
agent_all <- list()
agent_structured_all <- list()

for (ct_l2 in l2_types) {
  cat(sprintf("\n>> L2: %s\n", ct_l2))

  pb <- aggregate_pseudobulk(obj, ct_l2)
  if (is.null(pb)) {
    cat("  [SKIP] Insufficient pseudobulk\n")
    next
  }

  pb_tissues <- unique(na.omit(pb$meta$tissue))
  if (length(pb_tissues) < 2) {
    cat("  [SKIP] < 2 tissues\n")
    next
  }

  cat(sprintf(
    "  Pseudobulk samples: %d (%s)\n",
    nrow(pb$meta),
    paste(pb_tissues, collapse = ", ")
  ))

  pb_de_all[[ct_l2]] <- list()
  enrich_all[[ct_l2]] <- list()
  agent_all[[ct_l2]] <- list()
  agent_structured_all[[ct_l2]] <- list()

  for (pair in combn(as.character(pb_tissues), 2, simplify = FALSE)) {
    t1 <- pair[1]
    t2 <- pair[2]
    comp_name <- paste0(t2, "_vs_", t1)
    cat(sprintf("  DESeq2: %s\n", comp_name))

    res <- run_deseq2_pairwise(pb, t1, t2)
    if (is.null(res)) {
      cat("    [SKIP] Too few samples per tissue\n")
      next
    }

    pb_de_all[[ct_l2]][[comp_name]] <- res
    cat(sprintf(
      "    DEGs: %d up, %d down (samples: %d vs %d)\n",
      res$n_up,
      res$n_down,
      res$n_samples_1,
      res$n_samples_2
    ))

    comp_dir <- file.path(de_dir, safe_name(ct_l2), safe_name(comp_name))
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(res$de_table, file.path(comp_dir, "DESeq2_results.csv"))

    vp <- res$de_table %>%
      mutate(label = ifelse(sig == "sig" & rank(padj) <= 20, gene, ""))
    pv <- ggplot(vp, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("sig" = "red", "ns" = "grey70")) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
      geom_hline(
        yintercept = -log10(PADJ_THR),
        linetype = "dashed",
        color = "blue"
      ) +
      geom_vline(
        xintercept = c(-LFC_THR, LFC_THR),
        linetype = "dashed",
        color = "blue"
      ) +
      labs(
        title = sprintf("DESeq2: %s (%s)", comp_name, ct_l2),
        x = "log2 Fold Change",
        y = "-log10(padj)"
      ) +
      theme_minimal()
    save_plot(pv, file.path(comp_dir, "volcano"), width = 8, height = 6)

    enrich_all[[ct_l2]][[comp_name]] <- list()
    agent_all[[ct_l2]][[comp_name]] <- list()
    agent_structured_all[[ct_l2]][[comp_name]] <- list()

    for (dir_name in c("up", "down")) {
      genes <- if (dir_name == "up") {
        res$de_table %>%
          filter(sig == "sig", log2FoldChange > 0) %>%
          arrange(desc(log2FoldChange)) %>%
          head(TOP_N_DEG_ENRICHMENT) %>%
          pull(gene)
      } else {
        res$de_table %>%
          filter(sig == "sig", log2FoldChange < 0) %>%
          arrange(log2FoldChange) %>%
          head(TOP_N_DEG_ENRICHMENT) %>%
          pull(gene)
      }

      if (length(genes) < 5) {
        enrich_all[[ct_l2]][[comp_name]][[dir_name]] <- NULL
        next
      }
      cat(sprintf("    Enrichment [%s]: %d genes\n", dir_name, length(genes)))

      tested <- res$tested_genes
      enr_list <- list(
        GO_BP = run_gmt_enrichment(genes, go_bp_t2g, "GO_BP", tested),
        GO_MF = run_gmt_enrichment(genes, go_mf_t2g, "GO_MF", tested),
        GO_CC = run_gmt_enrichment(genes, go_cc_t2g, "GO_CC", tested),
        KEGG = run_gmt_enrichment(genes, kegg_t2g, "KEGG", tested),
        Hallmark = run_gmt_enrichment(genes, hallmark_t2g, "Hallmark", tested),
        CellMarker = run_gmt_enrichment(
          genes,
          cellmarker_t2g,
          "CellMarker",
          tested
        ),
        PanglaoDB = run_gmt_enrichment(
          genes,
          panglaodb_t2g,
          "PanglaoDB",
          tested
        ),
        B_cell_custom = run_gmt_enrichment(
          genes,
          bcell_custom_t2g,
          "B_cell_custom",
          tested
        )
      )
      enr_list <- enr_list[!sapply(enr_list, is.null)]
      enrich_all[[ct_l2]][[comp_name]][[dir_name]] <- enr_list

      enr_out <- file.path(comp_dir, paste0("enrichment_", dir_name))
      dir.create(enr_out, showWarnings = FALSE)
      for (db in names(enr_list)) {
        er <- enr_list[[db]]
        if (nrow(as.data.frame(er)) > 0) {
          fwrite(as.data.frame(er), file.path(enr_out, paste0(db, ".csv")))
          tryCatch(
            {
              pdf(
                file.path(enr_out, paste0(db, "_dotplot.pdf")),
                width = 10,
                height = 8
              )
              print(dotplot(
                er,
                showCategory = 15,
                title = sprintf("%s %s (%s %s)", db, dir_name, ct_l2, comp_name)
              ))
              dev.off()
            },
            error = function(e) NULL
          )
        }
      }
      saveRDS(enr_list, file.path(enr_out, "all_enrichment.rds"))

      if (ENABLE_LLM) {
        best_er <- NULL
        best_db <- NA_character_
        for (db in c(
          "GO_BP",
          "GO_MF",
          "GO_CC",
          "Hallmark",
          "KEGG",
          "CellMarker",
          "PanglaoDB",
          "B_cell_custom"
        )) {
          if (
            !is.null(enr_list[[db]]) && nrow(as.data.frame(enr_list[[db]])) >= 3
          ) {
            best_er <- enr_list[[db]]
            best_db <- db
            break
          }
        }

        if (!is.null(best_er)) {
          tissue_ctx <- build_tissue_pair_context(t1, t2, ct_l2, dir_name)
          sig_de <- res$de_table %>% filter(sig == "sig")
          gene_fc <- setNames(sig_de$log2FoldChange, toupper(sig_de$gene))

          cat(sprintf("    interpret_agent [%s] (db: %s)\n", dir_name, best_db))
          ia_payload <- run_interpret_agent_safe(best_er, tissue_ctx, gene_fc)
          ia <- ia_payload$result
          agent_all[[ct_l2]][[comp_name]][[dir_name]] <- ia

          ia_structured <- normalize_interpret_agent_result(
            raw_result = ia,
            ct_l2 = ct_l2,
            comp_name = comp_name,
            dir_name = dir_name,
            source_db = sprintf(
              "primary=%s; evidence=%s",
              best_db,
              paste(names(enr_list), collapse = ", ")
            ),
            warnings = ia_payload$warnings,
            error_message = ia_payload$error,
            enrich_obj = best_er,
            enrich_list = enr_list,
            context_str = tissue_ctx,
            gene_fc = gene_fc,
            force_standardization = TRUE
          )
          agent_structured_all[[ct_l2]][[comp_name]][[
            dir_name
          ]] <- ia_structured

          saveRDS(
            ia_payload,
            file.path(
              comp_dir,
              paste0("interpret_agent_", dir_name, "_payload.rds")
            )
          )
          saveRDS(
            ia_structured,
            file.path(
              comp_dir,
              paste0("interpret_agent_", dir_name, "_structured.rds")
            )
          )
          if (nzchar(ia_structured$raw_text)) {
            writeLines(
              ia_structured$raw_text,
              file.path(
                comp_dir,
                paste0("interpret_agent_", dir_name, "_raw.txt")
              )
            )
          }
        }
      }
    }
  }
}

# ==============================================================================
# 9. Exploratory Wilcoxon (marker discovery only)
# ==============================================================================

cat("\n=== Exploratory Cell-level Wilcoxon ===\n")
cat(
  "NOTE: P-values are inflated (pseudoreplication). For marker discovery ONLY.\n\n"
)

wilcox_all <- list()
for (ct_l2 in l2_types) {
  wx <- run_wilcox_exploratory(obj, ct_l2)
  if (!is.null(wx) && length(wx) > 0) {
    wilcox_all[[ct_l2]] <- wx
    for (comp_name in names(wx)) {
      wx_out <- file.path(wx_dir, safe_name(ct_l2))
      dir.create(wx_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(
        wx[[comp_name]],
        file.path(wx_out, paste0(safe_name(comp_name), "_wilcox.csv"))
      )
    }
    cat(sprintf("  [OK] %s: %d comparisons\n", ct_l2, length(wx)))
  }
}

# ==============================================================================
# 10. Save Summary RDS
# ==============================================================================

cat("\n=== Saving Summary ===\n")
saveRDS(pb_de_all, file.path(rpt_dir, "pseudobulk_de_all.rds"))
saveRDS(wilcox_all, file.path(rpt_dir, "wilcox_exploratory_all.rds"))
saveRDS(enrich_all, file.path(rpt_dir, "enrichment_all.rds"))
saveRDS(agent_all, file.path(rpt_dir, "interpret_agent_all.rds"))
saveRDS(
  agent_structured_all,
  file.path(rpt_dir, "interpret_agent_structured_all.rds")
)
saveRDS(L3_TO_L2_REMAP, file.path(rpt_dir, "L3_to_L2_remap.rds"))

agent_structured_df <- flatten_interpretation_records(agent_structured_all)
if (nrow(agent_structured_df) > 0) {
  fwrite(
    agent_structured_df,
    file.path(rpt_dir, "interpret_agent_structured.tsv"),
    sep = "\t"
  )
}

write_interpretation_markdown(
  agent_structured_df,
  file.path(OUTPUT_DIR, "LLM_INTERPRETATION.md"),
  "# LLM Interpretation Summary (Normal Tissue Comparison)",
  include_raw = FALSE
)
write_interpretation_markdown(
  agent_structured_df,
  file.path(rpt_dir, "LLM_INTERPRETATION_FOR_LLM.md"),
  "# LLM Interpretation Structured Input",
  include_raw = TRUE
)

# ==============================================================================
# 11. Generate REPORT.md
# ==============================================================================

cat("\n=== Generating REPORT.md ===\n")

md <- character()
add <- function(...) md <<- c(md, paste0(...))

add("# B Cell Tissue Comparison Report (Normal Respiratory Tract)")
add("")
add("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M"))
add("")
add(
  "**Pipeline:** B Cell Tissue Comparison v2.2.2 (pseudobulk DESeq2, multi-database enrichment)"
)
add("")
add(
  "**Note:** This is a cross-site anatomical comparison of NORMAL tissues, NOT disease vs healthy."
)
add("")
add("---")
add("")

add("## 1. Data Overview")
add("")
add(sprintf("- **Input:** `%s`", basename(H5AD_PATH)))
add(sprintf("- **Total cells:** %s", format(ncol(obj), big.mark = ",")))
add(sprintf("- **Tissues:** %s", paste(tissues, collapse = ", ")))
add(sprintf(
  "- **L2 subtypes (%d):** %s",
  length(l2_types),
  paste(l2_types, collapse = ", ")
))
add(sprintf(
  "- **L3 source column:** `%s` (%d unique)",
  L3_SOURCE_COL,
  length(unique(meta[[LABEL_COL]]))
))
add("")
add("### L3 -> L2 Remapping")
add("")
add("| L3 (scanvi_pred) | L2 (merged) |")
add("|---|---|")
for (i in seq_along(L3_TO_L2_REMAP)) {
  add(sprintf("| %s | %s |", names(L3_TO_L2_REMAP)[i], L3_TO_L2_REMAP[i]))
}
add("")

add("## 2. Visualization")
add("")
add("### 2.1 UMAP")
add("![UMAP tissue](figures/umap_tissue.png)")
add("")
add("![UMAP L2](figures/umap_celltype_L2.png)")
add("")
add("![UMAP L3](figures/umap_celltype_L3.png)")
add("")
add("![UMAP split](figures/umap_L2_split_tissue.png)")
add("")
add("### 2.2 Marker Dotplot")
add("![Dotplot](figures/dotplot_markers.png)")
add("")
add("### 2.3 Cell Composition")
add("![Composition](figures/composition_tissue_L2.png)")
add("")
add("![Counts](figures/count_tissue_L2.png)")
add("")
add("### 2.4 Sample-level Composition")
add("![Sample composition](figures/composition_sample_level.png)")
add("")
add("### 2.5 Top Marker Heatmap")
add("![Heatmap](figures/heatmap_top_markers.png)")
add("")

add("## 3. Pseudobulk DESeq2 (Primary Inference)")
add("")
add(sprintf(
  "Statistical unit: pseudobulk (sum counts per sample x tissue x L2). Thresholds: padj < %s, |log2FC| > %s.",
  PADJ_THR,
  LFC_THR
))
add("")
add("| L2 Subtype | Comparison | Samples (ref/case) | Up | Down | Total |")
add("|---|---|---|---|---|---|")
for (ct_l2 in names(pb_de_all)) {
  for (comp_name in names(pb_de_all[[ct_l2]])) {
    r <- pb_de_all[[ct_l2]][[comp_name]]
    if (is.null(r)) {
      next
    }
    add(sprintf(
      "| %s | %s | %d / %d | %d | %d | %d |",
      ct_l2,
      comp_name,
      r$n_samples_1,
      r$n_samples_2,
      r$n_up,
      r$n_down,
      r$n_up + r$n_down
    ))
  }
}
add("")

add("## 4. Multi-Database Enrichment")
add("")
add(
  "Databases: GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, B cell custom markers"
)
add("")
for (ct_l2 in names(enrich_all)) {
  for (comp_name in names(enrich_all[[ct_l2]])) {
    for (dir_name in names(enrich_all[[ct_l2]][[comp_name]])) {
      enr_l <- enrich_all[[ct_l2]][[comp_name]][[dir_name]]
      if (is.null(enr_l) || length(enr_l) == 0) {
        next
      }
      add(sprintf("### %s | %s | %s", ct_l2, comp_name, dir_name))
      add("")
      for (db in names(enr_l)) {
        er <- enr_l[[db]]
        if (is.null(er) || nrow(as.data.frame(er)) == 0) {
          next
        }
        top5 <- head(as.data.frame(er), 5)
        add(sprintf("**%s (top 5):**", db))
        add("")
        add("| Term | p.adjust | Count |")
        add("|---|---|---|")
        for (j in seq_len(nrow(top5))) {
          add(sprintf(
            "| %s | %.2e | %s |",
            top5$Description[j],
            top5$p.adjust[j],
            top5$Count[j]
          ))
        }
        add("")
      }
    }
  }
}

add("## 5. LLM Interpretation (interpret_agent)")
add("")
if (!ENABLE_LLM) {
  add("**SKIPPED:** DEEPSEEK_API_KEY not set. Re-run with API key to enable.")
  add("")
} else {
  if (nrow(agent_structured_df) == 0) {
    add(
      "No interpret_agent results generated (insufficient enrichment or all calls failed)."
    )
    add("")
  } else {
    add("Structured outputs:")
    add(sprintf(
      "- `%s`",
      file.path("reports", "interpret_agent_structured.tsv")
    ))
    add(sprintf("- `%s`", "LLM_INTERPRETATION.md"))
    add(
      "- Output normalization: interpret_agent result retried first, then standardized by"
    )
    add("  a secondary DeepSeek pass to a fixed JSON schema when needed.")
    add("- LLM output language: Chinese (Simplified).")
    add(
      "- LLM emphasis: top-|log2FC| genes linked to specific enriched pathways."
    )
    add("")
    for (i in seq_len(nrow(agent_structured_df))) {
      rec <- agent_structured_df[i, , drop = FALSE]
      add(sprintf(
        "### %s | %s | %s",
        rec$celltype_l2,
        rec$comparison,
        rec$direction
      ))
      add("")
      add(sprintf("- **Status:** %s", rec$status))
      add(sprintf(
        "- **Source DB:** %s",
        ifelse(nzchar(rec$source_db), rec$source_db, "NA")
      ))
      if (nzchar(rec$overview)) {
        add(sprintf("- **Overview:** %s", rec$overview))
      }
      if (nzchar(rec$key_mechanisms)) {
        add(sprintf("- **Key Mechanisms:** %s", rec$key_mechanisms))
      }
      if (nzchar(rec$hypothesis)) {
        add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
      }
      add("")
    }
  }
}

add("## 6. Methods")
add("")
add(
  "- **DE:** Pseudobulk DESeq2 (sample-level aggregation; Squair et al. 2021 Nat Commun)"
)
add(
  "- **Exploratory:** Cell-level Wilcoxon rank-sum (marker discovery only, NOT inference)"
)
add(
  "- **Enrichment:** clusterProfiler::enricher() + MSigDB GMT + CellMarker + PanglaoDB + custom B cell markers (symbol-based, zero ID loss)"
)
add(
  "- **Enrichment universe:** DESeq2-tested genes intersected with GMT (not full GMT)"
)
add(
  "- **Enrichment databases (8):** GO BP, GO MF, GO CC, KEGG, Hallmark, CellMarker, PanglaoDB, B cell custom"
)
if (ENABLE_LLM) {
  add(
    "- **LLM:** interpret_agent with DeepSeek deepseek-reasoner (tissue-pair-specific context)"
  )
  add(
    "- **LLM post-processing:** secondary deepseek-chat pass standardizes outputs to fixed JSON schema when interpret_agent output is partial, raw, or malformed"
  )
  add(
    "- **LLM output language:** Chinese (Simplified); gene symbols retained in English"
  )
  add(
    "- **LLM gene emphasis:** top-|log2FC| genes explicitly linked to enriched pathways in key_mechanisms"
  )
} else {
  add("- **LLM:** SKIPPED (no API key)")
}
if (!is.null(DESIGN_COVARIATE)) {
  add(sprintf("- **DESeq2 design:** ~ %s + tissue", DESIGN_COVARIATE))
} else {
  add("- **DESeq2 design:** ~ tissue")
}
add(
  "- **L2 remapping:** L3 (cell_type_scanvi_pred) merged into 4 L2 categories"
)
add(sprintf(
  "- **Context:** Normal tissue comparison across %s",
  paste(tissues, collapse = ", ")
))
add("")
add("### Output Objects")
add("")
add(
  "- `bcell_tissue_comparison_final.rds` -- Seurat object with `cell_type_L2` + `cell_type_L3`"
)
add(
  "- `bcell_tissue_comparison_final.h5ad` -- AnnData object with `cell_type_L2` + `cell_type_L3`"
)
add("")
add("---")
add("*Generated by bcell_tissue_comparison_v2_2_2.R*")

writeLines(md, file.path(OUTPUT_DIR, "REPORT.md"))
cat("[OK] REPORT.md written\n")

# ==============================================================================
# 12. Save Final Object (RDS + h5ad) with Standardized Columns
# ==============================================================================

cat("\n=== Saving Final Object (RDS + h5ad) ===\n")

if (!"cell_type_L3" %in% colnames(obj@meta.data)) {
  obj@meta.data[["cell_type_L3"]] <- as.character(obj@meta.data[[
    L3_SOURCE_COL
  ]])
  cat(sprintf("[OK] Created cell_type_L3 from '%s'\n", L3_SOURCE_COL))
} else {
  cat("[OK] cell_type_L3 already exists\n")
}

stopifnot("cell_type_L2" %in% colnames(obj@meta.data))
stopifnot("cell_type_L3" %in% colnames(obj@meta.data))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L2"]])))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L3"]])))

cat(sprintf(
  "[OK] cell_type_L2: %d unique (%s)\n",
  length(unique(obj@meta.data[["cell_type_L2"]])),
  paste(sort(unique(obj@meta.data[["cell_type_L2"]])), collapse = ", ")
))
cat(sprintf(
  "[OK] cell_type_L3: %d unique (%s)\n",
  length(unique(obj@meta.data[["cell_type_L3"]])),
  paste(sort(unique(obj@meta.data[["cell_type_L3"]])), collapse = ", ")
))

rds_path <- file.path(OUTPUT_DIR, "bcell_tissue_comparison_final.rds")
saveRDS(obj, rds_path)
cat(sprintf(
  "[OK] RDS saved: %s (%.1f MB)\n",
  basename(rds_path),
  file.size(rds_path) / 1e6
))

h5ad_out_path <- file.path(OUTPUT_DIR, "bcell_tissue_comparison_final.h5ad")

tryCatch(
  {
    anndata <- reticulate::import("anndata")
    scipy_sparse <- reticulate::import("scipy.sparse")
    np <- reticulate::import("numpy")

    counts_mat <- GetAssayData(obj, layer = "counts")

    if (inherits(counts_mat, "dgCMatrix")) {
      counts_t <- Matrix::t(counts_mat)
    } else {
      counts_t <- Matrix::t(as(counts_mat, "CsparseMatrix"))
    }

    counts_csc <- as(counts_t, "CsparseMatrix")
    X_scipy <- scipy_sparse$csr_matrix(
      reticulate::tuple(
        np$array(counts_csc@x, dtype = np$float32),
        np$array(counts_csc@i, dtype = np$int32),
        np$array(counts_csc@p, dtype = np$int32)
      ),
      shape = reticulate::tuple(
        as.integer(nrow(counts_csc)),
        as.integer(ncol(counts_csc))
      )
    )

    meta_export <- obj@meta.data
    for (col in colnames(meta_export)) {
      if (is.factor(meta_export[[col]])) {
        meta_export[[col]] <- as.character(meta_export[[col]])
      }
    }
    obs_df <- reticulate::r_to_py(meta_export)
    var_py <- reticulate::r_to_py(data.frame(row.names = rownames(obj)))

    adata <- anndata$AnnData(X = X_scipy, obs = obs_df, var = var_py)

    for (red_name in Reductions(obj)) {
      emb <- Embeddings(obj, reduction = red_name)
      adata$obsm[[paste0("X_", red_name)]] <- np$array(emb, dtype = np$float32)
    }

    adata$write_h5ad(h5ad_out_path, compression = "gzip")
    cat(sprintf(
      "[OK] h5ad saved: %s (%.1f MB)\n",
      basename(h5ad_out_path),
      file.size(h5ad_out_path) / 1e6
    ))

    adata_check <- anndata$read_h5ad(h5ad_out_path)
    stopifnot(
      "cell_type_L2" %in% reticulate::py_to_r(adata_check$obs$columns$tolist())
    )
    stopifnot(
      "cell_type_L3" %in% reticulate::py_to_r(adata_check$obs$columns$tolist())
    )
    cat(sprintf(
      "[OK] h5ad verified: %d cells x %d genes, cell_type_L2 + cell_type_L3 present\n",
      reticulate::py_to_r(adata_check$n_obs),
      reticulate::py_to_r(adata_check$n_vars)
    ))
    rm(
      adata,
      adata_check,
      counts_mat,
      counts_t,
      counts_csc,
      X_scipy,
      meta_export
    )
    gc()
  },
  error = function(e) {
    cat(sprintf("[ERROR] h5ad export failed: %s\n", e$message))
    cat("[INFO] RDS saved successfully; convert manually if h5ad needed.\n")
  }
)

cat("\n")

# ==============================================================================
# 13. Final Summary
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("B CELL TISSUE COMPARISON COMPLETE (v2.2.2)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Directory structure:\n")
cat(
  "  REPORT.md                              <- structured report (png embeds)\n"
)
cat(
  "  LLM_INTERPRETATION.md                  <- structured LLM interpretation summary (Chinese)\n"
)
cat(
  "  bcell_tissue_comparison_final.rds       <- Seurat object (cell_type_L2 + cell_type_L3)\n"
)
cat(
  "  bcell_tissue_comparison_final.h5ad      <- AnnData object (cell_type_L2 + cell_type_L3)\n"
)
cat(
  "  figures/                               <- UMAP, dotplot, heatmap, composition (pdf+png)\n"
)
cat(
  "  pseudobulk_de/                         <- DESeq2 per L2 x tissue pair + volcano + enrichment (8 DBs)\n"
)
cat(
  "  wilcox_exploratory/                    <- cell-level wilcox (marker discovery only)\n"
)
cat(
  "  reports/                               <- RDS summary objects + structured LLM files\n"
)
cat("\n")

cat("v2.2.2 Changes:\n")
cat(
  "  [LLM-1] build_tissue_pair_context(): Chinese output instruction + gene-pathway emphasis\n"
)
cat(
  "  [LLM-2] standardize_result_with_llm(): Rule 6 (Chinese output) + Rule 7 (gene->pathway format)\n"
)
cat("\n")
cat("v2.2.1 Fixes (inherited):\n")
cat("  [P1-1] UMAP output without raster / alpha; png dpi 300\n")
cat("  [P1-2] aggregate_pseudobulk sparse-friendly (no dense pre-alloc)\n")
cat("  [P1-3] Unmapped L3 -> stop() fail-fast\n")
cat("  [P1-4] CellMarker/PanglaoDB strict B cell term filter\n")
cat("  [P1-5] interpret_agent retry + atomic-result-tolerant normalization\n")
cat("  [P1-6] Secondary DeepSeek standardization to fixed JSON schema\n")
cat("  [P2-6] Heatmap title: 'visualization only'\n")
cat("  [P2-7] set.seed(42) in stratified_downsample\n")
cat("\n")

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)
cat("[OK] Done\n")
