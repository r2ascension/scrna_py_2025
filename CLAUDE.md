# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a single-cell RNA sequencing (scRNA-seq) analysis repository focused on T cell and immune cell analysis in nasal/respiratory tissues. The codebase uses **Seurat v5** as the primary framework, with integration of Python-based tools (scanpy/scANVI) via `reticulate`.

## Key Analysis Frameworks

### Core Stack
- **Seurat v5**: Primary single-cell analysis framework
- **Python integration**: Uses `reticulate` with conda environment `bbknn_env` for scanpy/scANVI workflows
- **Monocle3**: Trajectory analysis (T cell differentiation)
- **DoubletFinder + celda/DecontX**: Quality control
- **DESeq2**: Pseudobulk differential expression
- **MASC (Mixed models)**: Cell composition analysis

### File Format Conversions
- Primary interchange format: **h5ad** (AnnData) ↔ **RDS** (Seurat)
- Use `GetSeurat()` function in `GetSeurat.R` for h5ad → Seurat conversion with raw counts validation
- Use `SeuratDisk::Convert()` or custom functions for Seurat → h5ad

## Critical Analysis Functions

### 1. Data Import/Export

**GetSeurat() - h5ad to Seurat conversion** (`GetSeurat.R`)
```r
# Reads h5ad with automatic raw counts detection
seurat_obj <- GetSeurat(
  h5ad_path = "path/to/file.h5ad",
  assay = "RNA",
  prefer_raw = TRUE,           # Try adata.raw.X first
  prefer_layer_counts = TRUE,   # Then adata.layers['counts']
  validate_counts = TRUE,       # Validate raw counts (not log-transformed)
  debug = TRUE
)
```
- Automatically searches for raw counts in: `adata.raw.X` → `adata.layers['counts']` → `adata.X`
- Validates counts are truly raw (integer-like, reasonable range)
- Imports `obsm` as Seurat reductions (UMAP, PCA)
- Imports `obs` as metadata

### 2. Quality Control Pipeline

**Complete QC + Smart Merge** (`rds_folder_qc_smartmerge_20251213_v1.R`)

Command-line usage:
```bash
Rscript rds_folder_qc_smartmerge_20251213_v1.R /path/to/input_rds_dir /path/to/output_dir
```

Pipeline:
1. Gene name standardization (SYMBOL → ALIAS mapping via `org.Hs.eg.db`)
2. QC filters: `nFeature_RNA` (200-6000), `percent.mt` (<20%), `percent.ribo` (<40%)
3. Regression-based outlier detection
4. DoubletFinder (default 6% doublet rate)
5. DecontX (ambient RNA removal, contamination <25%)
6. Smart merge across samples with batch-level gene availability matrix

Key parameters:
```r
NFEATURE_MIN <- 200
NFEATURE_MAX <- 6000
MT_MAX <- 20
RB_MAX <- 40
DOUBLET_RATE <- 0.06
DECONTX_CONTAM_MAX <- 0.25
```

### 3. Trajectory Analysis

**T Cell Trajectory with Monocle3** (`tcell_monocle3_trajectory_analysis_v1.1.R`)

Features:
- CD4 trajectory: `CD4_Naive → CD4_CM → CD4_EM → Treg`
- CD8 trajectory: `CD8_Naive → CD8_CM → CD8_EM → CD8_TEMRA`
- T cell purity QC (validates expression of CD3D/E/G)
- Branch-specific gene analysis
- Pseudotime validation plots

Key configuration:
```r
CD4_CELLTYPES <- c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg")
CD8_CELLTYPES <- c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA")
TCELL_MARKERS <- c("CD3D", "CD3E", "CD3G")
NAIVE_MARKERS <- c("CCR7", "SELL", "LEF1", "TCF7", "IL7R")
EFFECTOR_MARKERS <- c("GZMK", "GZMB", "PRF1", "GNLY", "NKG7")
```

### 4. Differential Analysis

**Pseudobulk DESeq2 Analysis** (`pseudobulk.R`)
```r
run_tissue_comparison_analysis(
  seurat_obj = T_object,
  cell_anno_col = "Annotation_2",
  tissue_col = "tissue",
  sample_col = "sample",
  min_cell_per_sample = 3,
  min_sample_per_tissue = 3,
  run_gsva = TRUE,
  run_go = TRUE,
  output_dir = "./results"
)
```
- Aggregates counts by tissue × sample × cell type
- Runs DESeq2 for tissue comparisons
- GSVA pathway analysis (MSigDB)
- GO enrichment

**MASC Cell Composition Analysis** (`scMASC.R`)
```r
scPairwiseMASCAnalysis(
  seurat_obj,
  cell_type_col = "Annotation",
  sample_col = "sample",       # Random effect
  contrast_col = "tissue",      # Fixed effect
  min_cells = 10,
  min_samples = 2,
  p_threshold = 0.05,
  output_dir = "pairwise_MASC_analysis"
)
```
- Mixed-effects models for cell composition shifts
- Pairwise tissue comparisons
- Handles batch effects via random effects

### 5. Enrichment Analysis

**Multi-format GO/KEGG enrichment** (`enrichment_functions.R`)
```r
# Accepts CSV files or data.frames with DESeq2/Seurat output
# Auto-maps column names: p_val, logFC, p_val_adj, gene_names
# Outputs GO/KEGG plots + interactive HTML visualizations
```

## Environment Setup

### Required Conda Environment
```r
library(reticulate)
use_condaenv("bbknn_env", required = TRUE)
py_config()
```

Python packages in `bbknn_env`:
- scanpy
- scvi-tools (scANVI)
- bbknn
- anndata

### Key R Packages (`Installation.R`)

Bioconductor:
- `org.Hs.eg.db`, `clusterProfiler`, `KEGGREST`
- `DropletUtils`, `SingleCellExperiment`
- `SCENIC`, `infercnv`

GitHub:
- `mojaveazure/seurat-disk` (SeuratDisk)
- `cellgeni/sceasy`
- `PaulingLiu/ROGUE`
- `Danko-Lab/BayesPrism/BayesPrism`
- `campbio/celda` (DecontX)

CRAN:
- `Seurat`, `monocle3`, `DoubletFinder`
- `DESeq2`, `GSVA`, `msigdbr`
- `lme4` (MASC models)

## Common Workflows

### Workflow 1: Import Multi-dataset from Python
```r
source("GetSeurat.R")

# Read h5ad files
seurat_obj1 <- GetSeurat("dataset1.h5ad", debug = TRUE)
seurat_obj2 <- GetSeurat("dataset2.h5ad", debug = TRUE)

# Add metadata
seurat_obj1$dataset <- "Study1"
seurat_obj1$tissue <- "nose"

# Merge
merged <- merge(seurat_obj1, seurat_obj2,
                add.cell.ids = c("Study1", "Study2"))
saveRDS(merged, "merged.rds")
```

### Workflow 2: QC Pipeline for Multiple Samples
```bash
# Input: folder with sample1.rds, sample2.rds, ...
# Output: cleaned samples + merged object
Rscript rds_folder_qc_smartmerge_20251213_v1.R \
  /data/raw_samples \
  /data/cleaned_output
```

### Workflow 3: T Cell Trajectory Analysis
```r
# Input: h5ad with T cell subset + cell type annotations
# Edit configuration section in tcell_monocle3_trajectory_analysis_v1.1.R
# Set INPUT_H5AD, OUTPUT_DIR, LABELS_KEY, CD4_CELLTYPES, CD8_CELLTYPES
Rscript tcell_monocle3_trajectory_analysis_v1.1.R
```

## Code Architecture Patterns

### Seurat v5 Compatibility
The codebase handles Seurat v4/v5 differences with fallback patterns:

```r
# Get counts matrix (works across Seurat versions)
get_counts_matrix <- function(seurat_obj, assay = "RNA") {
  counts <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),  # v5
    error = function(e1) tryCatch(
      GetAssayData(seurat_obj, assay = assay, layer = "counts"),  # v5 alt
      error = function(e2) tryCatch(
        GetAssayData(seurat_obj, assay = assay, slot = "counts"),  # v4
        error = function(e3) seurat_obj[[assay]]@counts  # v3
      )
    )
  )
  return(counts)
}

# Create assay (v5 vs v4)
create_assay_object <- function(counts_matrix) {
  assay <- tryCatch(
    CreateAssay5Object(counts = counts_matrix),  # v5
    error = function(e) CreateAssayObject(counts = counts_matrix)  # v4
  )
  return(assay)
}
```

### Gene Name Standardization
Always map to official HGNC symbols before merging datasets:

```r
# Uses org.Hs.eg.db
# Layer 1: Direct SYMBOL match
# Layer 2: ALIAS match
# Layer 3: Keep original if no match
# Handles synonyms: merge duplicate genes by summing counts
```

## Data Locations and File Naming

Typical analysis structure (based on existing scripts):
```
/home/h2048/data/
├── py/                          # Python/scanpy outputs
│   └── YYYYMMDD/
│       └── analysis_name/
│           └── adata_*.h5ad
├── R/                           # R/Seurat analyses
│   └── YYYYMMDD/
│       └── sample_name.rds
└── script/R/                    # This repository
```

File naming conventions:
- h5ad files: `adata_<dataset>_<method>_<annotation>.h5ad`
- RDS files: `<dataset>_<date>.rds` or `<study_author>_<year>.rds`
- Output folders: Use ISO date prefix `YYYYMMDD` or versioned suffix `_v1`, `_v2`

## Key Metadata Columns

Standard columns in analysis scripts:
- **Cell identity**: `Annotation`, `Annotation_2`, `Multinomial_Label`, `celltype`
- **Sample/Batch**: `sample`, `sample_id`, `dataset`, `batch`
- **Tissue**: `tissue`, `Location`
- **Disease**: `disease`, `COVID_status`
- **QC metrics**: `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`
- **Doublet filtering**: `DF.classifications`, `decontX_contamination`

## Disease Focus

Primary research focus: **Chronic Rhinosinusitis with Nasal Polyps (CRSwNP)**

Common tissue types:
- Nasal brush/biopsy samples
- Healthy vs disease comparisons
- Multi-dataset integration across different studies

Cell type hierarchies emphasize:
- T cell subsets (CD4/CD8 differentiation states)
- Myeloid cells (monocytes, macrophages)
- Epithelial cells
- Fibroblasts
- Endothelial cells

## Best Practices

1. **Always validate raw counts** when importing h5ad files
   - Use `validate_counts = TRUE` in `GetSeurat()`
   - Check max/mean values are in expected UMI range

2. **Gene standardization before merging**
   - Different datasets may use different gene IDs (Ensembl vs symbols)
   - Use `map_gene_names()` function before integration

3. **Quality control order**
   - EmptyDrops (if not already done in Python)
   - Basic QC (nFeature, MT%, Ribo%)
   - Doublet removal
   - Ambient RNA removal (DecontX)

4. **Batch effects**
   - For visualization: Harmony, scVI, BBKNN (Python)
   - For DE analysis: Include batch as random effect (pseudobulk/MASC)

5. **Pseudobulk requirements**
   - Minimum 3 cells per sample per cell type
   - Minimum 3 samples per condition
   - Use sample-level aggregation (not cell-level)

6. **Trajectory analysis prerequisites**
   - Subset to relevant cell types only (e.g., T cells)
   - Ensure raw counts are available
   - Validate marker gene expression for QC
