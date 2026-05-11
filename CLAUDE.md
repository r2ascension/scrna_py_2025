# CLAUDE.md - Single-Cell RNA-seq Analysis Repository

This file provides guidance to Claude Code when working with this repository.

## Repository Overview

This is a **single-cell RNA sequencing (scRNA-seq) analysis repository** focused on immune cell analysis in nasal/respiratory tissues, particularly **Chronic Rhinosinusitis with Nasal Polyps (CRSwNP)**.

The repository implements a **dual-language workflow**:
- **Python** (`py/`): Deep learning-based integration and annotation using scanpy/scvi-tools
- **R** (`R/`): Statistical analysis, QC, and interpretation using Seurat v5

---

## Quick Navigation

| Folder | Purpose | Key Technologies |
|--------|---------|------------------|
| `py/` | Batch correction, cell annotation, visualization | scVI, scANVI, CellTypist, BBKNN, scanpy |
| `R/` | QC pipeline, statistical analysis, DE, trajectory | Seurat v5, Monocle3, DESeq2, MASC |

---

## Python Pipeline (`py/`)

### Core Three-Stage Architecture

All major analysis scripts follow this pattern:

1. **scVI** - Variational autoencoder for batch effect removal
   - Uses raw UMI counts (`layers['counts']`)
   - Generates latent space (`X_scvi`)

2. **CellTypist** - Automated cell type classification
   - Requires log-normalized counts
   - Uses Human_Lung_Atlas model (gene symbols)

3. **scANVI** - Semi-supervised label refinement
   - Initialized from trained scVI model
   - Low-confidence cells (<0.5) marked as "Unknown"

### Key Data Structure

```python
# AnnData object conventions
adata.X                    # Processed (log-normalized)
adata.layers['counts']     # Raw UMI counts (CRITICAL for scVI)
adata.raw.X                # Full gene set raw counts

# Common obs columns
adata.obs['celltype']              # Cell type annotations
adata.obs['Multinomial_Label']     # scANVI predictions
adata.obs['majority_voting']       # CellTypist predictions
```

### File Naming Conventions

- `*_scvi_celltypist_scanvi_pipeline_YYYYMMDD_v*.py` - Full pipelines
- `*_bbknn_*.py` - BBKNN integration scripts
- `*_subcluster_*.py` - Cell type-specific subclustering
- `*_visualization_*.ipynb` - Jupyter notebooks for figures

---

## R Pipeline (`R/`)

### Core Analysis Framework

- **Seurat v5**: Primary single-cell analysis framework
- **Python integration**: Uses `reticulate` with conda env `bbknn_env`
- **Monocle3**: Trajectory analysis (T cell differentiation)
- **DoubletFinder + DecontX**: Quality control
- **DESeq2 + MASC**: Differential expression and composition analysis

### Critical Functions

#### 1. Data Import: `GetSeurat()`
```r
seurat_obj <- GetSeurat(
  h5ad_path = "path/to/file.h5ad",
  prefer_raw = TRUE,           # Try adata.raw.X first
  prefer_layer_counts = TRUE,   # Then adata.layers['counts']
  validate_counts = TRUE,       # Validate raw counts
  debug = TRUE
)
```

#### 2. QC Pipeline: `rds_folder_qc_smartmerge_*.R`
```bash
Rscript rds_folder_qc_smartmerge_20251213_v1.R \
  /data/raw_samples /data/cleaned_output
```

Pipeline steps:
1. Gene name standardization (SYMBOL → ALIAS via `org.Hs.eg.db`)
2. QC filters: nFeature_RNA (200-6000), percent.mt (<20%), percent.ribo (<40%)
3. DoubletFinder (6% doublet rate)
4. DecontX (ambient RNA removal, contamination <25%)
5. Smart merge across samples

#### 3. Trajectory Analysis: `tcell_monocle3_trajectory_*.R`
- CD4 trajectory: `CD4_Naive → CD4_CM → CD4_EM → Treg`
- CD8 trajectory: `CD8_Naive → CD8_CM → CD8_EM → CD8_TEMRA`

#### 4. Differential Analysis
```r
# Pseudobulk DESeq2
run_tissue_comparison_analysis(
  seurat_obj = T_object,
  cell_anno_col = "Annotation_2",
  tissue_col = "tissue",
  sample_col = "sample",
  min_cell_per_sample = 3
)

# MASC Cell Composition
scPairwiseMASCAnalysis(
  seurat_obj,
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue"
)
```

### Seurat v5 Compatibility Pattern

```r
# Get counts matrix (works across versions)
get_counts_matrix <- function(seurat_obj, assay = "RNA") {
  counts <- tryCatch(
    LayerData(seurat_obj, assay = assay, layer = "counts"),  # v5
    error = function(e1) tryCatch(
      GetAssayData(seurat_obj, assay = assay, slot = "counts"),  # v4
      error = function(e2) seurat_obj[[assay]]@counts  # v3
    )
  )
  return(counts)
}
```

---

## Data Locations

### Directory Structure (Latest by Category)

```
/home/h2048/data/
├── py/                    # Python scanpy/scVI outputs (by date)
│   ├── 0214/              # Latest: Stromal/T cell visualization
│   │   └── stromal_marker_visualization/
│   ├── 0212/              # T cell annotation viz
│   │   └── tcell_annotation_viz/
│   ├── 0209/              # Myeloid validation optimized
│   │   ├── myeloid_validation_optimized/
│   │   └── epithelial_viz_L3_v1_0/
│   ├── 0208/              # Merged scANVI L2 (v2.5.5 HOTFIX)
│   │   └── merged_scanvi_L2_prod_v2_5_5_HOTFIX/
│   ├── 0204/              # scArches mapping L2
│   │   └── scarches_mapping_L2_v2_5_3/
│   ├── 0129/              # T/NK unified analysis
│   │   └── tnk_analysis_unified/
│   └── (earlier dates: 0128, 0127, 0121, 0119, 0118...)
│
├── R/                     # R Seurat outputs (by date)
│   ├── 0210/              # Latest: Epithelial/Stromal interpretation
│   │   ├── epithelial_interpret_v3_2_1/
│   │   └── stromal_interpret_v3_2_1/
│   ├── 0209/              # Myeloid/B cell analysis
│   ├── 0205/              # Stromal analysis
│   ├── 0204/              # Integrated analysis
│   ├── 0131/              # B cell merge viz
│   ├── 0130/              # T cell analysis
│   └── (earlier dates: 0129, 0128, 0127, 0126...)
│
├── bulk/                  # Bulk RNA-seq data
├── core_data/             # Core reference datasets
│   └── cellranger/        # Cell Ranger reference genomes
├── index_genome/          # Reference genomes and indexes
│   └── cisTarget_databases/
└── source/                # Raw data sources
```

### Key Analysis Outputs by Cell Type

| Cell Type | Latest Python Output | Latest R Output |
|-----------|---------------------|-----------------|
| **T cells** | `py/0212/tcell_annotation_viz/` | `R/0210/` (interpretation) |
| **B cells** | `py/0203/bcell_scarches_v4_1/` | `R/0131/` |
| **Myeloid** | `py/0209/myeloid_validation_optimized/` | `R/0209/` |
| **Epithelial** | `py/0209/epithelial_viz_L3_v1_0/` | `R/0210/epithelial_interpret_v3_2_1/` |
| **Stromal** | `py/0214/stromal_marker_visualization/` | `R/0210/stromal_interpret_v3_2_1/` |
| **All cells** | `py/0208/merged_scanvi_L2_prod_v2_5_5_HOTFIX/` | - |

### Output File Patterns

**Python outputs (`py/YYYYMMDD/`):**
- `adata_*_results.h5ad` - Main AnnData with scVI/scANVI latents
- `*_scanvi_model/` - Saved scANVI models
- `*_umap_operator.joblib` - UMAP transformers for projection
- `*_markers.csv` - Differential expression results
- `figures/*.pdf` - Visualization outputs

**R outputs (`R/YYYYMMDD/`):**
- `*_seurat.rds` - Seurat objects
- `*_markers.csv` - FindAllMarkers results
- `*_degs.csv` - DESeq2 differential expression
- `figures/*.pdf` - Plots and visualizations

---

## Environment Setup

### Python Environment (`bbknn_env`)
```python
# Required packages
scanpy, scvi-tools, celltypist, bbknn, anndata
```

### R Environment
```r
library(reticulate)
use_condaenv("bbknn_env", required = TRUE)

# Key packages
# Bioconductor: org.Hs.eg.db, clusterProfiler, DropletUtils
# GitHub: mojaveazure/seurat-disk, cellgeni/sceasy, PaulingLiu/ROGUE
# CRAN: Seurat, monocle3, DoubletFinder, DESeq2, GSVA
```

---

## Common Workflows

### 1. Python → R Handoff
```python
# In Python
adata.write_h5ad("output.h5ad")
```
```r
# In R
source("R/GetSeurat.R")
seurat_obj <- GetSeurat("output.h5ad", validate_counts = TRUE)
```

### 2. Complete Analysis Pipeline
1. **Python**: QC + scVI + CellTypist + scANVI → h5ad
2. **R**: Import with GetSeurat → Subset → DE/Trajectory

### 3. Cell Type Subclustering
1. **Python**: Subset cell type → BBKNN + scANVI
2. **R**: Import → Annotation refinement → Interpretation

---

## Key Metadata Columns

| Column | Description |
|--------|-------------|
| `celltype` / `Annotation` / `Annotation_2` | Cell identity |
| `sample` / `sample_id` | Sample identifier |
| `tissue` / `Location` | Tissue type |
| `disease` / `COVID_status` | Disease status |
| `nFeature_RNA`, `nCount_RNA` | QC metrics |
| `percent.mt`, `percent.ribo` | QC metrics |
| `DF.classifications` | Doublet status |
| `decontX_contamination` | Ambient RNA level |

---

## Best Practices

1. **Always validate raw counts** when importing h5ad files
2. **Use gene standardization** before merging datasets
3. **QC order**: Basic QC → Doublet removal → Ambient RNA removal
4. **Batch effects**: Use Harmony/scVI for visualization; include batch as random effect for DE
5. **Pseudobulk**: Minimum 3 cells per sample per cell type; minimum 3 samples per condition
6. **File naming**: Use ISO date prefix `YYYYMMDD` and version suffix `_v1`, `_v2`

---

## Cell Type Focus

Primary cell types analyzed:
- **T cells**: CD4/CD8 subsets, differentiation states
- **B cells**: Subclustering and activation states
- **Myeloid cells**: Monocytes, macrophages
- **Epithelial cells**: AT1, AT2, ciliated, goblet, club, basal
- **Stromal cells**: Fibroblasts, endothelial

---

## See Also

- `py/CLAUDE.md` - Detailed Python pipeline documentation
- `R/CLAUDE.md` - Detailed R pipeline documentation
