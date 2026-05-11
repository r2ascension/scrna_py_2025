# ==============================================================================
# Epithelial Subcluster Enrichment Interpretation Analysis
# Using clusterProfiler::interpret() for Deep Biological Insight
# ==============================================================================
#
# Purpose:
#   - Load epithelial subcluster results from Python analysis
#   - Perform CellMarker + GO enrichment analysis
#   - Use LLM-powered interpret() for:
#     * Cell subtype annotation
#     * Functional state characterization
#     * Mechanistic narrative generation
#
# Author: r2end
# Date: 2025-01-26
# Version: v1.0
#
# ==============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input files
H5AD_PATH <- "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0126/epithelial_interpret"

# Cell type/subcluster columns
CELLTYPE_COL <- "celltypist_pred"
SUBCLUSTER_COL <- "subcluster" # Format: "celltype_clusterID"
LEIDEN_COL <- "subcluster_leiden" # Pure cluster ID

# Analysis parameters
MIN_CELLS_FOR_ANALYSIS <- 50 # Minimum cells for reliable enrichment
TOP_N_MARKERS <- 50 # Top N markers per cluster for enrichment
MARKER_LOGFC_THRESHOLD <- 0.5
MARKER_PADJ_THRESHOLD <- 0.05

# LLM API settings (DeepSeek recommended)
USE_DEEPSEEK <- TRUE
API_KEY_ENV <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05" # Or set directly: API_KEY <- "your_key"

# CellMarker database
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.xlsx"

# Output options
SAVE_INDIVIDUAL_REPORTS <- TRUE # Save per-celltype HTML reports
GENERATE_COMBINED_REPORT <- TRUE # Combined markdown report
SAVE_ENRICHMENT_OBJECTS <- TRUE # Save .rds for re-analysis

# =============================================================================
# LOAD LIBRARIES
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Epithelial Subcluster Interpretation Analysis\n")
cat(
  "================================================================================\n\n"
)

cat("[INFO] Loading required libraries...\n")

suppressPackageStartupMessages({
  library(reticulate) # Python interface
  library(SCNT) # GetSeurat function
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  # library(rio) # For reading Excel
  library(fanyi) # For LLM API setup
})

cat("\n[OK] Libraries loaded\n\n")

# =============================================================================
# PYTHON / CONDA ENVIRONMENT (for GetSeurat)
# =============================================================================
library(reticulate)
cat("[INFO] Configuring Python environment...\n")
use_condaenv("bbknn_env", required = TRUE)
py_config()
cat("[OK] Python environment configured\n\n")

# =============================================================================
# SETUP OUTPUT DIRECTORY
# =============================================================================

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("[INFO] Created output directory: %s\n", OUTPUT_DIR))
}

# Create subdirectories
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "enrichment_objects"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "marker_genes"), showWarnings = FALSE)

# =============================================================================
# SETUP LLM API
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("LLM API Configuration\n")
cat(
  "================================================================================\n\n"
)

if (USE_DEEPSEEK) {
  # Try to get API key from environment or set directly
  api_key <- Sys.getenv(API_KEY_ENV)

  if (api_key == "") {
    cat("[WARN] No API key found in environment variable:", API_KEY_ENV, "\n")
    cat(
      "[WARN] You can set it via: Sys.setenv(",
      API_KEY_ENV,
      " = 'your_key')\n"
    )
    cat("[WARN] Or pass it directly in interpret() calls\n")
    USE_LLM <- FALSE
  } else {
    cat("[OK] Found API key in environment\n")

    # Configure fanyi for DeepSeek
    tryCatch(
      {
        set_translate_api("deepseek", api_key = api_key)
        cat("[OK] DeepSeek API configured via fanyi\n")
        USE_LLM <- TRUE
      },
      error = function(e) {
        cat("[WARN] Failed to configure fanyi:", conditionMessage(e), "\n")
        USE_LLM <- FALSE
      }
    )
  }
} else {
  cat("[INFO] LLM interpretation disabled\n")
  USE_LLM <- FALSE
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Loading Epithelial Subcluster Data\n")
cat(
  "================================================================================\n\n"
)

cat("[INFO] Loading h5ad using GetSeurat()...\n")
cat("  Input:", H5AD_PATH, "\n")

# Load Seurat object directly from h5ad
seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat("\n[OK] Data loaded successfully\n")
cat(sprintf("  Cells: %d\n", ncol(seurat_obj)))
cat(sprintf("  Genes: %d\n", nrow(seurat_obj)))
cat(sprintf("  Assays: %s\n", paste(names(seurat_obj@assays), collapse = ", ")))

# Check required columns
required_cols <- c(CELLTYPE_COL, SUBCLUSTER_COL, LEIDEN_COL)
missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))

if (length(missing_cols) > 0) {
  stop(
    "[ERROR] Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("\n[INFO] Cell type distribution:\n")
celltype_counts <- table(seurat_obj@meta.data[[CELLTYPE_COL]])
print(celltype_counts)

cat("\n[INFO] Subcluster distribution:\n")
subcluster_counts <- table(seurat_obj@meta.data[[SUBCLUSTER_COL]])
cat(sprintf("  Total subclusters: %f\n", length(subcluster_counts)))
cat(sprintf("  Median size: %f cells\n", median(subcluster_counts)))
cat(sprintf(
  "  Range: %f - %f cells\n",
  min(subcluster_counts),
  max(subcluster_counts)
))

# =============================================================================
# LOAD CELLMARKER DATABASE
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Loading CellMarker Database\n")
cat(
  "================================================================================\n\n"
)
library(rio)
cat("[INFO] Reading CellMarker database from:", CELLMARKER_PATH, "\n")

cellmarker_db <- import(CELLMARKER_PATH)

cat(sprintf("[OK] Loaded %d records\n", nrow(cellmarker_db)))

# Check required columns
if (!all(c("cell_name", "marker") %in% colnames(cellmarker_db))) {
  stop("[ERROR] CellMarker database must have 'cell_name' and 'marker' columns")
}

# Filter for relevant tissue types (optional)
cat("\n[INFO] Tissue types in database:\n")
if ("tissue_type" %in% colnames(cellmarker_db)) {
  tissue_counts <- table(cellmarker_db$tissue_type)
  print(head(sort(tissue_counts, decreasing = TRUE), 20))

  # Optional: filter for respiratory tissues
  respiratory_keywords <- c(
    "Airway",
    "Lung",
    "Nasal",
    "Respiratory",
    "Bronch",
    "Alveol"
  )
  cellmarker_filtered <- cellmarker_db %>%
    filter(grepl(
      paste(respiratory_keywords, collapse = "|"),
      tissue_type,
      ignore.case = TRUE
    ))

  if (nrow(cellmarker_filtered) > 0) {
    cat(sprintf(
      "\n[INFO] Filtered to %d respiratory-relevant records\n",
      nrow(cellmarker_filtered)
    ))
    cellmarker_db <- cellmarker_filtered
  }
}

# Prepare TERM2GENE format
cellmarker_term2gene <- cellmarker_db %>%
  select(cell_name, marker) %>%
  distinct()

cat(sprintf(
  "[OK] Prepared TERM2GENE with %d cell type - gene pairs\n",
  nrow(cellmarker_term2gene)
))
library(reticulate)
cat("[INFO] Configuring Python environment...\n")
use_condaenv("bbknn_env", required = TRUE)
py_config()
cat("[OK] Python environment configured\n\n")

# =============================================================================
# FIND MARKER GENES FOR ALL SUBCLUSTERS
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Computing Marker Genes\n")
cat(
  "================================================================================\n\n"
)

cat("[INFO] Setting default assay and running FindAllMarkers...\n")

# Set RNA assay as default (contains raw counts)
DefaultAssay(seurat_obj) <- "RNA"

# Ensure data is normalized
if (!"data" %in% names(seurat_obj@assays$RNA@layers)) {
  cat("[INFO] Normalizing data...\n")
  seurat_obj <- NormalizeData(seurat_obj)
}

cat("[INFO] Running FindAllMarkers...\n")
cat(sprintf("  Grouping by: %s\n", SUBCLUSTER_COL))
cat(sprintf(
  "  Thresholds: logFC >= %.2f, padj <= %.3f\n",
  MARKER_LOGFC_THRESHOLD,
  MARKER_PADJ_THRESHOLD
))

all_markers <- FindAllMarkers(
  seurat_obj,
  group.by = SUBCLUSTER_COL,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = MARKER_LOGFC_THRESHOLD,
  test.use = "wilcoxon",
  verbose = FALSE
)

cat(sprintf(
  "[OK] Found %d marker genes across all subclusters\n",
  nrow(all_markers)
))

# Filter by adjusted p-value
all_markers <- all_markers %>%
  filter(p_val_adj <= MARKER_PADJ_THRESHOLD)

cat(sprintf("[OK] After p-value filter: %d marker genes\n", nrow(all_markers)))

# Save marker genes
marker_output <- file.path(
  OUTPUT_DIR,
  "marker_genes",
  "all_subcluster_markers.csv"
)
write.csv(all_markers, marker_output, row.names = FALSE)
cat(sprintf("[OK] Saved markers to: %s\n", marker_output))

# Prepare top markers per cluster for enrichment
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup()

cat(sprintf(
  "[INFO] Selected top %d markers per cluster for enrichment\n",
  TOP_N_MARKERS
))

# =============================================================================
# ENRICHMENT ANALYSIS: CELLMARKER + GO
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Enrichment Analysis (CellMarker + GO)\n")
cat(
  "================================================================================\n\n"
)

# Get unique subclusters with sufficient cells
subcluster_sizes <- table(seurat_obj@meta.data[[SUBCLUSTER_COL]])
valid_subclusters <- names(subcluster_sizes[
  subcluster_sizes >= MIN_CELLS_FOR_ANALYSIS
])

cat(sprintf(
  "[INFO] Analyzing %d subclusters with >= %d cells\n",
  length(valid_subclusters),
  MIN_CELLS_FOR_ANALYSIS
))

# Prepare data for compareCluster
top_markers_filtered <- top_markers %>%
  filter(cluster %in% valid_subclusters) %>%
  select(gene, cluster)

# CellMarker enrichment
cat("\n[INFO] Running CellMarker enrichment...\n")

cellmarker_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers_filtered,
  fun = enricher,
  TERM2GENE = cellmarker_term2gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

cat(sprintf("[OK] CellMarker enrichment complete\n"))
cat(sprintf(
  "  Significant terms: %d\n",
  sum(cellmarker_enrich@compareClusterResult$p.adjust < 0.05)
))

# Save enrichment object
if (SAVE_ENRICHMENT_OBJECTS) {
  saveRDS(
    cellmarker_enrich,
    file.path(OUTPUT_DIR, "enrichment_objects", "cellmarker_enrichment.rds")
  )
}

# GO Biological Process enrichment
cat("\n[INFO] Running GO Biological Process enrichment...\n")

go_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers_filtered,
  fun = enrichGO,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

cat(sprintf("[OK] GO enrichment complete\n"))
cat(sprintf(
  "  Significant terms: %d\n",
  sum(go_enrich@compareClusterResult$p.adjust < 0.05)
))

# Save enrichment object
if (SAVE_ENRICHMENT_OBJECTS) {
  saveRDS(
    go_enrich,
    file.path(OUTPUT_DIR, "enrichment_objects", "go_enrichment.rds")
  )
}

# =============================================================================
# VISUALIZATION: ENRICHMENT RESULTS
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Generating Enrichment Visualizations\n")
cat(
  "================================================================================\n\n"
)

# CellMarker dotplot
cat("[INFO] Creating CellMarker dotplot...\n")

pdf(
  file.path(OUTPUT_DIR, "figures", "cellmarker_enrichment_dotplot.pdf"),
  width = 14,
  height = 10
)
print(
  dotplot(cellmarker_enrich, showCategory = 10, font.size = 8) +
    ggtitle("CellMarker Enrichment by Subcluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

cat("[OK] Saved: cellmarker_enrichment_dotplot.pdf\n")

# GO dotplot
cat("[INFO] Creating GO dotplot...\n")

pdf(
  file.path(OUTPUT_DIR, "figures", "go_enrichment_dotplot.pdf"),
  width = 14,
  height = 12
)
print(
  dotplot(go_enrich, showCategory = 15, font.size = 7) +
    ggtitle("GO Biological Process Enrichment by Subcluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

cat("[OK] Saved: go_enrichment_dotplot.pdf\n")

# =============================================================================
# LLM INTERPRETATION: CELL SUBTYPE ANNOTATION
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("LLM-Powered Interpretation\n")
cat(
  "================================================================================\n\n"
)

if (!USE_LLM) {
  cat("[WARN] LLM interpretation skipped (API not configured)\n")
  cat("[INFO] You can still manually inspect enrichment results\n")
} else {
  # Task 1: Cell subtype annotation
  cat("[INFO] Task 1: Cell subtype annotation (using CellMarker + GO)\n")
  cat("  This will take several minutes...\n\n")

  annotation_results <- tryCatch(
    {
      interpret(
        list(cellmarker_enrich, go_enrich),
        context = "Epithelial cells from chronic rhinosinusitis with nasal polyps (CRSwNP) study. These are subclusters of major epithelial cell types including basal, secretory, ciliated, and goblet cells.",
        task = "annotation"
      )
    },
    error = function(e) {
      cat("[ERROR] Annotation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(annotation_results)) {
    cat("\n[OK] Annotation complete!\n\n")

    # Print summary
    cat(
      "================================================================================\n"
    )
    cat("ANNOTATION SUMMARY\n")
    cat(
      "================================================================================\n\n"
    )

    for (i in seq_along(annotation_results)) {
      result <- annotation_results[[i]]
      cat(sprintf("Cluster: %s\n", names(annotation_results)[i]))
      cat(sprintf("  Cell Type: %s\n", result$cell_type))
      cat(sprintf("  Confidence: %s\n", result$confidence))
      cat("\n")
    }

    # Save annotation results
    annotation_df <- data.frame(
      subcluster = names(annotation_results),
      cell_type = sapply(annotation_results, function(x) x$cell_type),
      confidence = sapply(annotation_results, function(x) x$confidence),
      stringsAsFactors = FALSE
    )

    write.csv(
      annotation_df,
      file.path(OUTPUT_DIR, "subcluster_annotations_llm.csv"),
      row.names = FALSE
    )

    cat("[OK] Saved: subcluster_annotations_llm.csv\n")

    # Save full report
    saveRDS(
      annotation_results,
      file.path(OUTPUT_DIR, "reports", "annotation_results_full.rds")
    )

    cat("[OK] Saved: annotation_results_full.rds\n")
  }

  # Task 2: Functional phenotyping
  cat("\n[INFO] Task 2: Functional phenotyping\n")
  cat("  This will take several minutes...\n\n")

  phenotype_results <- tryCatch(
    {
      interpret(
        list(cellmarker_enrich, go_enrich),
        context = "Epithelial cells from CRSwNP patients showing Type 2 inflammation. Looking for functional states like inflammatory, activated, remodeling, or damaged states.",
        task = "phenotyping"
      )
    },
    error = function(e) {
      cat("[ERROR] Phenotyping failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(phenotype_results)) {
    cat("\n[OK] Phenotyping complete!\n\n")

    # Print summary
    cat(
      "================================================================================\n"
    )
    cat("PHENOTYPE SUMMARY\n")
    cat(
      "================================================================================\n\n"
    )

    for (i in seq_along(phenotype_results)) {
      result <- phenotype_results[[i]]
      cat(sprintf("Cluster: %s\n", names(phenotype_results)[i]))
      cat(sprintf("  Phenotype: %s\n", result$phenotype))
      cat(sprintf("  Confidence: %s\n", result$confidence))
      cat("\n")
    }

    # Save phenotype results
    phenotype_df <- data.frame(
      subcluster = names(phenotype_results),
      phenotype = sapply(phenotype_results, function(x) x$phenotype),
      confidence = sapply(phenotype_results, function(x) x$confidence),
      stringsAsFactors = FALSE
    )

    write.csv(
      phenotype_df,
      file.path(OUTPUT_DIR, "subcluster_phenotypes_llm.csv"),
      row.names = FALSE
    )

    cat("[OK] Saved: subcluster_phenotypes_llm.csv\n")

    # Save full report
    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results_full.rds")
    )

    cat("[OK] Saved: phenotype_results_full.rds\n")
  }

  # Task 3: Mechanistic interpretation (example for one major cell type)
  cat("\n[INFO] Task 3: Mechanistic interpretation (example for Basal cells)\n")

  # Filter for basal cell subclusters
  basal_subclusters <- grep("^Basal", names(annotation_results), value = TRUE)

  if (length(basal_subclusters) > 0) {
    cat(sprintf(
      "  Analyzing %d Basal cell subclusters...\n",
      length(basal_subclusters)
    ))

    # Extract results for basal cells only
    basal_cellmarker <- cellmarker_enrich@compareClusterResult %>%
      filter(Cluster %in% basal_subclusters)

    basal_go <- go_enrich@compareClusterResult %>%
      filter(Cluster %in% basal_subclusters)

    # Create subset enrichment objects
    basal_cm_obj <- cellmarker_enrich
    basal_cm_obj@compareClusterResult <- basal_cellmarker

    basal_go_obj <- go_enrich
    basal_go_obj@compareClusterResult <- basal_go

    mechanism_results <- tryCatch(
      {
        interpret(
          list(basal_cm_obj, basal_go_obj),
          context = "Basal epithelial cells from CRSwNP. Basal cells are progenitor cells that can differentiate into secretory and ciliated cells. Interested in understanding differentiation trajectories and inflammatory responses.",
          task = "interpretation"
        )
      },
      error = function(e) {
        cat(
          "[ERROR] Mechanism interpretation failed:",
          conditionMessage(e),
          "\n"
        )
        return(NULL)
      }
    )

    if (!is.null(mechanism_results)) {
      cat("\n[OK] Mechanistic interpretation complete!\n")

      # Save
      saveRDS(
        mechanism_results,
        file.path(OUTPUT_DIR, "reports", "basal_mechanism_interpretation.rds")
      )

      cat("[OK] Saved: basal_mechanism_interpretation.rds\n")

      # Print narrative (if available)
      if (!is.null(mechanism_results$narrative)) {
        cat("\n")
        cat(
          "================================================================================\n"
        )
        cat("NARRATIVE FOR PAPER (Basal Cells)\n")
        cat(
          "================================================================================\n\n"
        )
        cat(mechanism_results$narrative)
        cat("\n\n")
      }
    }
  } else {
    cat(
      "  [INFO] No Basal cell subclusters found, skipping mechanistic interpretation\n"
    )
  }
}

# =============================================================================
# GENERATE COMBINED REPORT
# =============================================================================

if (GENERATE_COMBINED_REPORT) {
  cat("\n")
  cat(
    "================================================================================\n"
  )
  cat("Generating Combined Report\n")
  cat(
    "================================================================================\n\n"
  )

  report_file <- file.path(OUTPUT_DIR, "INTERPRETATION_REPORT.md")

  cat("[INFO] Writing report to:", report_file, "\n")

  sink(report_file)

  cat("# Epithelial Subcluster Interpretation Report\n\n")
  cat("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("---\n\n")

  cat("## Overview\n\n")
  cat(
    "This report contains LLM-powered interpretations of epithelial cell subclusters\n"
  )
  cat("from chronic rhinosinusitis with nasal polyps (CRSwNP) analysis.\n\n")

  cat("### Dataset Summary\n\n")
  cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("- Subclusters analyzed: %d\n", length(valid_subclusters)))
  cat(sprintf("- Minimum cells per cluster: %d\n", MIN_CELLS_FOR_ANALYSIS))
  cat("\n")

  if (USE_LLM && !is.null(annotation_results)) {
    cat("---\n\n")
    cat("## Cell Subtype Annotations\n\n")

    for (i in seq_along(annotation_results)) {
      result <- annotation_results[[i]]
      cat(sprintf("### %s\n\n", names(annotation_results)[i]))
      cat(sprintf("**Cell Type:** %s  \n", result$cell_type))
      cat(sprintf("**Confidence:** %s  \n\n", result$confidence))

      if (!is.null(result$reasoning)) {
        cat("**Reasoning:**\n\n")
        cat(result$reasoning, "\n\n")
      }

      if (!is.null(result$supporting_markers)) {
        cat("**Supporting Markers/Pathways:**\n\n")
        for (marker in result$supporting_markers) {
          cat(sprintf("- %s\n", marker))
        }
        cat("\n")
      }

      cat("---\n\n")
    }
  }

  if (USE_LLM && !is.null(phenotype_results)) {
    cat("## Functional Phenotypes\n\n")

    for (i in seq_along(phenotype_results)) {
      result <- phenotype_results[[i]]
      cat(sprintf("### %s\n\n", names(phenotype_results)[i]))
      cat(sprintf("**Phenotype:** %s  \n", result$phenotype))
      cat(sprintf("**Confidence:** %s  \n\n", result$confidence))

      if (!is.null(result$reasoning)) {
        cat("**Reasoning:**\n\n")
        cat(result$reasoning, "\n\n")
      }

      if (!is.null(result$key_processes)) {
        cat("**Key Processes:**\n\n")
        for (process in result$key_processes) {
          cat(sprintf("- %s\n", process))
        }
        cat("\n")
      }

      cat("---\n\n")
    }
  }

  cat("## Files Generated\n\n")
  cat("- `subcluster_annotations_llm.csv` - LLM cell type annotations\n")
  cat("- `subcluster_phenotypes_llm.csv` - LLM functional phenotypes\n")
  cat("- `all_subcluster_markers.csv` - Marker genes for all clusters\n")
  cat("- `cellmarker_enrichment.rds` - CellMarker enrichment object\n")
  cat("- `go_enrichment.rds` - GO enrichment object\n")
  cat("- `annotation_results_full.rds` - Complete annotation results\n")
  cat("- `phenotype_results_full.rds` - Complete phenotype results\n")
  cat("\n")

  sink()

  cat("[OK] Report saved\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("ANALYSIS COMPLETE\n")
cat(
  "================================================================================\n\n"
)

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("Key files:\n")
cat("  - INTERPRETATION_REPORT.md\n")
cat("  - subcluster_annotations_llm.csv\n")
cat("  - subcluster_phenotypes_llm.csv\n")
cat("  - figures/cellmarker_enrichment_dotplot.pdf\n")
cat("  - figures/go_enrichment_dotplot.pdf\n")
cat("  - marker_genes/all_subcluster_markers.csv\n\n")

if (USE_LLM) {
  cat("LLM interpretation: SUCCESS\n")
  cat("  You can load .rds files to explore detailed results:\n")
  cat("    results <- readRDS('reports/annotation_results_full.rds')\n")
  cat("    print(results[[1]])  # View first cluster's interpretation\n\n")
} else {
  cat("LLM interpretation: SKIPPED (configure API key to enable)\n\n")
}

cat("To re-run with custom parameters, edit the CONFIGURATION section\n")
cat("and source this script again.\n\n")

cat(
  "================================================================================\n"
)
cat("END OF SCRIPT\n")
cat(
  "================================================================================\n"
)
