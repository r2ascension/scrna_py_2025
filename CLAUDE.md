# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a **single-cell RNA-seq analysis pipeline** repository for lung tissue data, implementing state-of-the-art computational biology workflows using Python/scanpy/scvi-tools ecosystem. The primary focus is on **cell type annotation and batch correction** for large-scale single-cell datasets.

## Core Pipeline Architecture

### Three-Stage Deep Learning Pipeline Pattern

All major analysis scripts follow this standardized architecture:

1. **scVI (Batch Correction)** - Variational autoencoder for batch effect removal
   - Trains on raw UMI counts (integer values, NOT normalized)
   - Generates latent space (`X_scvi`) for batch-corrected representation
   - CRITICAL: Always use `layer='counts'` or `layers['counts']`

2. **CellTypist (Automated Annotation)** - Machine learning-based cell type classification
   - Requires log-normalized counts (NOT raw counts)
   - Expects gene symbols (NOT ENSEMBL IDs) for Human_Lung_Atlas model
   - Produces `predicted_labels`, `majority_voting`, and confidence scores

3. **scANVI (Semi-supervised Refinement)** - VAE with classifier for label refinement
   - Initialized from trained scVI model
   - Uses CellTypist predictions as reference labels
   - Low-confidence cells (<0.5) marked as "Unknown" for refinement

### Data Structure Conventions

**AnnData object layers:**
- `.X` - Processed expression (typically log-normalized from BBKNN preprocessing)
- `.layers['counts']` - **CRITICAL**: Raw UMI counts (required for scVI/scANVI)
- `.raw.X` - Full gene set raw counts (often shared memory with `layers['counts']`)

**Important metadata columns:**
- `dataset` or `sample` - Batch identifier for correction
- `tissue` or `tissue_sampling_method` - Tissue type annotation
- `cell_type` - Major cell type assignments (from BBKNN preprocessing)

## Common Development Commands

### Running Analysis Pipelines

```bash
# Cell-type-specific pipelines (recommended for >100k cells per type)
python epithelial_scvi_celltypist_scanvi_20251212_v3_4.py
python tcell_scvi_celltypist_scanvi_pipeline_20251206_v2.py
python stromal_vascular_scvi_scanvi_pipeline_20251212_v2.py
python myeloid_scvi_celltypist_scanvi_pipeline_v1.py
python b_scvi_celltypist_scanvi_pipeline_v1.py

# All-cells pipeline (uses HVG optimization for memory efficiency)
python allcells_scvi_celltypist_scanvi_pipeline_20251208_v2.1.py

# Legacy BBKNN-based analysis (older approach, not recommended)
python bbknn_annotation_universal.py
```

### GPU Verification

```bash
# Check GPU availability (REQUIRED before running pipelines)
python -c "import torch; print(f'GPU: {torch.cuda.is_available()}')"
python gputest.py
```

### Jupyter Notebook Development

```bash
# Launch notebooks for interactive analysis
jupyter notebook complete_epithelial_pipeline_v2.ipynb
jupyter notebook tcell_analysis_with_starcat_v1.4.ipynb
jupyter notebook Epithelial_scANVI_Training.ipynb
```

## Critical Implementation Patterns

### 1. Gene Name Handling

**Gene conversion priority (ENFORCED across all v2+ pipelines):**

```python
# CRITICAL: Gene conversion BEFORE creating .raw
normalize_gene_names(adata)  # Adds 'symbol_base' column
preserve_full_raw_if_missing(adata)  # Then create .raw

# Priority order:
# 1. Local metadata columns: 'gene_symbol', 'gene_symbols', 'symbol'
# 2. mygene.org conversion (ENSEMBL → HGNC symbols)
# 3. Fallback: use var_names as-is

# For CellTypist compatibility:
# - Human_Lung_Atlas.pkl expects HGNC gene symbols
# - Use base symbols WITHOUT -1/-2 suffixes (deduplication artifacts)
# - Build overlap using adata.raw.var['symbol_base']
```

### 2. HVG (Highly Variable Genes) Selection

**Two approaches based on dataset size:**

```python
# For large datasets (>400k cells) - MEMORY OPTIMIZED
N_HVG = 4000  # 2k-5k recommended
sc.pp.highly_variable_genes(adata, layer='counts', n_top_genes=N_HVG,
                             batch_key='dataset', flavor='seurat_v3', subset=False)

# Critical: Create SEPARATE adata_model for training
adata_model = sc.AnnData(
    X=adata.layers['counts'][:, hvg_mask].copy(),
    obs=adata.obs[['dataset']].copy(),
    var=adata.var.loc[hvg_mask].copy()
)
# Benefits: ~90% memory reduction, 10-20x training speedup

# For smaller datasets (<100k cells) - ALL GENES
# Use full gene set directly (no subsetting)
```

### 3. scVI Model Training/Loading Pattern

```python
# ALWAYS preserve HVG list for model reuse
hvg_file = output_dir / "hvg_genes.txt"
pd.Series(hvg_genes).to_csv(hvg_file, index=False, header=False)

# Reuse existing model (production pattern):
if PRETRAINED_SCVI_MODEL and Path(PRETRAINED_SCVI_MODEL).exists():
    genes_for_training = load_gene_list_for_pretrained(PRETRAINED_SCVI_MODEL)
    # Must subset new data to EXACT same genes
    hvg_mask = adata.var['symbol_base'].isin(genes_for_training)
    scvi_model = scvi.model.SCVI.load(PRETRAINED_SCVI_MODEL, adata=adata_model)
else:
    # Train new model
    scvi_model = scvi.model.SCVI(adata_model, n_latent=75, n_layers=3, ...)
    scvi_model.train(max_epochs=800, batch_size=256, ...)
    scvi_model.save(model_path)
```

### 4. Dual UMAP Management (CRITICAL FIX in v2.x)

**Problem:** scANVI UMAP was overwriting scVI UMAP in older versions.

**Solution:** Use `neighbors_key` to maintain separate embeddings:

```python
# scVI UMAP
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30, key_added='neighbors_scvi')
sc.tl.umap(adata, neighbors_key='neighbors_scvi')
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

# scANVI UMAP (separate)
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=30, key_added='neighbors_scanvi')
sc.tl.umap(adata, neighbors_key='neighbors_scanvi')
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()

# Result: Both UMAPs preserved in separate .obsm slots
```

### 5. Rare Cell Type Filtering

```python
# Stability improvement: Merge rare types to prevent scANVI training issues
MIN_CELLS_PER_TYPE = 10  # Threshold for rare types

def merge_rare_types(labels, min_cells=10, unknown_label="Unknown"):
    vc = labels.value_counts()
    rare_types = vc[vc < min_cells].index
    labels[labels.isin(rare_types)] = unknown_label
    return labels

# Apply BEFORE scANVI training
adata.obs['labels_for_scanvi'] = merge_rare_types(
    adata.obs['cell_type_celltypist_filt'],
    min_cells=MIN_CELLS_PER_TYPE
)
```

### 6. Index-Aligned Result Writing (v2.1+ Critical Fix)

```python
# WRONG (v1.x - position-based, prone to reordering bugs):
adata.obs['cell_type'] = predictions.predicted_labels['predicted_labels'].values

# CORRECT (v2.x - index-aligned):
pred_df = predictions.predicted_labels
pred_df = pred_df.reindex(adata.obs_names)  # Force alignment
adata.obs['cell_type'] = pred_df['predicted_labels'].astype(str).values
```

## File Naming Conventions

**Versioning pattern:** `{celltype}_scvi_celltypist_scanvi_{date}_v{version}.py`

Examples:
- `epithelial_scvi_celltypist_scanvi_20251212_v3_4.py` - Production (v3.4)
- `tcell_scvi_celltypist_scanvi_pipeline_20251206_v2.py` - Stable (v2.0)
- `stromal_vascular_scvi_scanvi_pipeline_20251212_v2.py` - Latest

**Checkpoint naming:** `adata_{celltype}_FINAL.h5ad` for production outputs

## Data Paths (HPC Environment)

**Input data:**
- `/home/h2048/data/py/1128/bbknn_celltype_analysis/{CellType}/adata_{CellType}_bbknn.h5ad`
- `/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad` (all cells)

**Output structure:**
```
/home/h2048/data/py/{YYMMDD}/{analysis_name}/
├── adata_{celltype}_FINAL.h5ad
├── models/
│   ├── scvi_model/
│   └── scanvi_model/
├── hvg_genes.txt (CRITICAL for model reuse)
├── figures/
└── README.md (auto-generated summary)
```

**Reference models:**
- `/home/h2048/data/source/reference/celltypist_models/Human_Lung_Atlas.pkl`
- `/home/h2048/data/source/reference/celltypist_models/Immune_All_Low.pkl`

## GPU Configuration

**Single-GPU enforcement (v2.1+ fix):**
```python
# CRITICAL: Explicit single GPU control (prevents multi-GPU conflicts)
if torch.cuda.is_available():
    train_kwargs.update({'accelerator': 'gpu', 'devices': 1})
else:
    train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})
```

**Reproducibility setup:**
```python
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():  # CPU-safe (v2.1 fix)
    torch.cuda.manual_seed_all(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED
scvi.settings.dl_num_workers = 0
```

## Major Cell Type Pipelines

### Epithelial Cells (v3.4 - Production)
- **File:** `epithelial_scvi_celltypist_scanvi_20251212_v3_4.py`
- **Features:** Robust gene conversion, HVG coordination with pretrained models
- **Expected subtypes:** Basal, AT1, AT2, Ciliated, Secretory, Ionocytes

### T/NK Cells (v2.0)
- **File:** `tcell_scvi_celltypist_scanvi_pipeline_20251206_v2.py`
- **Model:** Immune_All_Low.pkl
- **Expected subtypes:** CD4+ T, CD8+ T, Tregs, NK cells, NKT, MAIT

### Stromal/Vascular Cells (v2.3)
- **File:** `stromal_vascular_scvi_scanvi_pipeline_20251212_v2.py`
- **Target types:** Endothelial, Fibroblast, SMC (smooth muscle cells)
- **Marker genes:** PECAM1/CDH5 (endothelial), COL1A1/PDGFRA (fibroblasts), ACTA2/MYH11 (SMC)

### Myeloid Cells
- **File:** `myeloid_scvi_celltypist_scanvi_pipeline_v1.py`
- **Expected subtypes:** Macrophages (M1/M2), Monocytes, DCs, Mast cells

### B Cells
- **File:** `b_scvi_celltypist_scanvi_pipeline_v1.py`
- **Expected subtypes:** Naive B, Memory B, Plasma cells

## Common Issues and Solutions

### Issue: "layers['counts'] not found"
**Solution:** Check if raw counts are in `.raw.X` or `.X`. Extract to `layers['counts']`:
```python
if adata.raw is not None:
    adata.layers['counts'] = adata.raw.X.copy()
```

### Issue: "Model gene list mismatch"
**Solution:** Always load `hvg_genes.txt` when reusing pretrained models:
```python
hvg_file = model_dir.parent / "hvg_genes.txt"
hvg_genes = pd.read_csv(hvg_file, header=None)[0].tolist()
hvg_mask = adata.var['symbol_base'].isin(hvg_genes)
```

### Issue: "Low CellTypist gene overlap (<50%)"
**Solution:** Check gene name format. CellTypist expects symbols, not ENSEMBL IDs:
```python
sample_gene = str(adata.var_names[0])
if sample_gene.startswith('ENSG'):
    # Need conversion via mygene or local metadata
    normalize_gene_names(adata)
```

### Issue: "scANVI training unstable"
**Solution:** Filter rare cell types before training:
```python
merge_rare_types(adata.obs['labels_for_scanvi'], min_cells=10, other='Unknown')
```

### Issue: "Out of memory during training"
**Solution:** Use HVG optimization for large datasets:
```python
USE_HVG_FOR_SCVI = True
N_HVG = 4000  # Reduce from 58k genes
# Expected memory reduction: ~60-70%
```

## Critical Fixes in v2+ Pipelines

All modern pipelines (v2.0+) include these MANDATORY fixes:

1. **Gene alignment** - `normalize_gene_names()` BEFORE `.raw` creation
2. **UMAP separation** - Use `neighbors_key` to prevent overwriting
3. **CPU-safe CUDA** - Conditional `torch.cuda.manual_seed_all()`
4. **Filter on counts** - Gene filtering uses `layers['counts']`, not `.X`
5. **Index-aligned writes** - Use `.reindex()` when transferring predictions
6. **Single GPU control** - Explicit `devices=1` for scVI/scANVI
7. **Rare type handling** - Merge types with <10 cells to "Unknown"

## Dependencies

**Core packages:**
```
scanpy>=1.9
scvi-tools>=1.0
celltypist>=1.6
torch>=2.0 (with CUDA support)
mygene>=3.2
```

**Optional but recommended:**
```
bbknn>=1.5 (for BBKNN preprocessing)
harmony-pytorch (alternative batch correction)
```

## Performance Guidelines

**Dataset size vs. approach:**
- **<50k cells:** All genes, standard parameters
- **50k-200k cells:** HVG (4k genes), standard epochs
- **>200k cells:** HVG (4k genes), reduce epochs (scVI: 400→200, scANVI: 600→200)
- **>500k cells:** Consider cell-type-specific split pipelines

**Expected runtimes (on V100 GPU):**
- scVI training (100k cells, 4k HVG): ~30-60 min
- scANVI training (100k cells): ~20-40 min
- CellTypist annotation: ~5-10 min
- Total pipeline (100k cells): ~2-3 hours

## Code Quality Standards

When modifying pipelines:
1. Preserve version history in docstrings (list critical fixes)
2. Use explicit layer specifications (`layer='counts'` not `layer=None` with assumptions)
3. Always save HVG lists when training new models (`hvg_genes.txt`)
4. Add sanity checks (gene count validation, confidence score ranges)
5. Export statistics to CSV for audit trails
6. Use `adata.uns['pipeline_info']` to document parameters
