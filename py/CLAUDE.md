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
# Latest v3.5.x series with P0 fixes and production optimizations
python bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py
python myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py
python t_scvi_celltypist_scanvi_20260108_v3_5_1.py
python stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py

# All-cells pipeline (latest with covariates system)
python allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py

# Epithelial cells (dual scANVI architecture)
python allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py  # v2.7-PRODUCTION

# Subcluster analysis (post-integration, uses BBKNN graph)
python epithelial_subcluster_bbknn_pipeline_20260120_v3_4.py
python bcell_subcluster_analysis_v2_20260119.py
python myeloid_subcluster_analysis_v2_20260121.py
python stromal_subcluster_analysis_v2_20260121.py

# Quality control utilities
python celltypist_doublet_detection_20260114_v1_1.py
python remove_doublets_and_gse299751_20260114.py

# Complete integrated pipelines (BBKNN + Harmony + cNMF)
python complete_epithelial_analysis_pipeline_20260121_v2_1.py

# cNMF analysis
python label_guided_cnmf_pipeline_20260114_v1_1.py
python cnmf_results_analysis_20260114_v1_1.py

# Legacy BBKNN-based analysis (older approach, not recommended)
python bbknn_annotation_universal.py

# scArches Reference Mapping (2-step workflow)
# Step 1: Train reference model with UMAP operator
python step1_train_bcell_L2_20260204_v2_5_3.py
# Step 2: Map query and merge
python step2_map_query_20260204_v2_5_4.py

# Alternative: Generic scArches mapping
python scarches_mapping_20260127_v1_2_1.py

# Merged training (reference + query after mapping)
python merged_scanvi_training_20260208_v2_5_5.py

# Visualization and analysis (2026-02 series)
python visualize_scanvi_results_20260208.py
python visualize_query_results_20260204.py
python stromal_marker_visualization_20260214_v1_1.py
python epithelial_subcluster_visualization_L3_20260205_v1_0.py
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

### 7. Covariates System (v2.3.1+ NEW)

**Problem:** Technical variation (MT%, stress, cell cycle) confounds biological signals.

**Solution:** Pass as continuous covariates to scVI/scANVI models:

```python
# Calculate covariates
adata.obs['pct_counts_mt'] = (adata[:, adata.var_names.str.startswith('MT-')].X.sum(1).A1 /
                               adata.X.sum(1).A1 * 100)

# Stress signature score
stress_genes = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
sc.tl.score_genes(adata, stress_genes, score_name='stress_score')

# Cell cycle scoring
s_genes = [g for g in S_GENES if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Setup model with covariates
scvi.model.SCVI.setup_anndata(
    adata_model,
    layer='counts',
    batch_key='dataset',
    continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
)
```

**Stress signature genes** (from allcells v2.3.1):
```python
STRESS_SIGNATURE_GENES = [
    "ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3",
    "BANF1","CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2",
    # ... (full list in allcells_20260115_v2_3_1.py lines 148-162)
]
```

### 8. Unknown Cleaning (v2.3.1+ NEW)

**Problem:** Low-confidence scANVI predictions ("Unknown") may include low-quality cells.

**Solution:** Multi-criteria gating using neighborhood purity and quality metrics:

```python
# Thresholds
EXISTING_LABEL_PURITY_THRESHOLD = 0.5  # Neighborhood purity
EXISTING_LABEL_MT_THRESHOLD = 20       # mt% cutoff
EXISTING_LABEL_STRESS_PERCENTILE = 95  # stress percentile

# Calculate neighborhood purity (fraction of neighbors with same label)
from sklearn.neighbors import NearestNeighbors
nn = NearestNeighbors(n_neighbors=30, metric='euclidean')
nn.fit(adata.obsm['X_scanvi'])
_, indices = nn.kneighbors(adata.obsm['X_scanvi'])

purity = []
for i, neighbors in enumerate(indices):
    my_label = adata.obs['cell_type_scanvi'].iloc[i]
    neighbor_labels = adata.obs['cell_type_scanvi'].iloc[neighbors]
    purity.append((neighbor_labels == my_label).mean())

adata.obs['neighborhood_purity'] = purity

# Filter low-quality "Unknown" cells
unknown_mask = adata.obs['cell_type_scanvi'] == 'Unknown'
low_quality = (
    (adata.obs['neighborhood_purity'] < EXISTING_LABEL_PURITY_THRESHOLD) |
    (adata.obs['pct_counts_mt'] > EXISTING_LABEL_MT_THRESHOLD) |
    (adata.obs['stress_score'] > adata.obs['stress_score'].quantile(EXISTING_LABEL_STRESS_PERCENTILE/100))
)

# Remove low-quality unknowns
adata = adata[~(unknown_mask & low_quality)].copy()
```

### 9. Hierarchical Labeling (v2.3.1+ NEW)

**Problem:** Single-level scANVI training may over-split rare types or under-resolve major lineages.

**Solution:** Dual scANVI training at two granularity levels:

```python
# Step 1: Define major lineage mapping
MAJOR_LINEAGE_MAP = {
    'Basal': 'Basal_Lineage',
    'Suprabasal': 'Basal_Lineage',
    'Dividing_Basal': 'Basal_Lineage',

    'Ciliated': 'Ciliated_Lineage',
    'Deuterosome': 'Ciliated_Lineage',

    'Secretory_Goblet': 'Secretory_Lineage',
    'Secretory_Club': 'Secretory_Lineage',
    # ...
}

adata.obs['major_lineage'] = adata.obs['cell_type_celltypist'].map(MAJOR_LINEAGE_MAP)

# Step 2: Train scANVI at major lineage level
scanvi_major = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key='major_lineage',
    unlabeled_category='Unknown'
)
scanvi_major.train(max_epochs=200, ...)
scanvi_major.save(output_dir / "scanvi_major_model")

# Step 3: Train scANVI at fine-grained level
scanvi_fine = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key='cell_type_celltypist_filt',
    unlabeled_category='Unknown'
)
scanvi_fine.train(max_epochs=200, ...)
scanvi_fine.save(output_dir / "scanvi_fine_model")

# Result: Two annotation levels available
adata.obsm['X_scanvi_major'] = scanvi_major.get_latent_representation()
adata.obsm['X_scanvi_fine'] = scanvi_fine.get_latent_representation()
```

### 10. Class Imbalance Handling (v2.3.1+ NEW)

**Problem:** Rare cell types (<1% of dataset) cause scANVI training instability.

**Solution:** Use `n_samples_per_label` to balance training batches:

```python
# Calculate cell type frequencies
celltype_counts = adata.obs['labels_for_scanvi'].value_counts()
min_cells = celltype_counts.min()

# Set balanced sampling (prevents rare types from dominating gradients)
N_SAMPLES_PER_LABEL = min(min_cells, 50)  # Cap at 50 to avoid huge batch sizes

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key='labels_for_scanvi',
    unlabeled_category='Unknown'
)

scanvi_model.train(
    max_epochs=200,
    n_samples_per_label=N_SAMPLES_PER_LABEL,  # KEY PARAMETER
    batch_size=512,
    ...
)
```

### 11. Counts Recovery from External Files (v3.5.2+ NEW)

**Problem:** Preprocessing pipelines may lose `layers['counts']` after subsetting/filtering.

**Solution:** Recover counts by matching cell IDs from external h5ad:

```python
def recover_counts_from_external(adata, external_file_path, counts_layer_name='counts'):
    """
    Recover counts layer by matching cell IDs from external h5ad file
    """
    adata_ext = sc.read_h5ad(external_file_path)

    # Match cells
    matching_cells = set(adata.obs.index) & set(adata_ext.obs.index)

    # Match genes (handle ENSEMBL vs symbols)
    gene_map = {}  # current_index -> external_index
    for i, curr_gene in enumerate(adata.var_names):
        for j, ext_gene in enumerate(adata_ext.var_names):
            if curr_gene == ext_gene:
                gene_map[i] = j
                break

    # Build recovered matrix
    from scipy import sparse
    recovered = sparse.lil_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)

    for cid in matching_cells:
        ci = list(adata.obs.index).index(cid)
        ei = list(adata_ext.obs.index).index(cid)

        counts_vec = adata_ext.layers['counts'][ei, :].toarray().flatten()
        for curr_g, ext_g in gene_map.items():
            recovered[ci, curr_g] = counts_vec[ext_g]

    adata.layers['counts'] = recovered.tocsr()
    return adata

# Usage
adata = recover_counts_from_external(
    adata,
    "/path/to/original/data_with_counts.h5ad"
)
```

### 12. HPC Stability and Thread Control (v2.5+ NEW)

**Problem:** Multi-threading conflicts on HPC clusters cause deadlocks and memory issues.

**Solution:** Explicit single-thread configuration:

```python
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# scVI-specific stability
scvi.settings.dl_num_workers = 0  # Avoid DataLoader multiprocessing issues
```

### 13. L2/L3 Hierarchical Labeling (v2.5.3+ NEW)

**Problem:** Fine-grained L3 labels are too specific for stable scANVI training, but coarse L2 labels lose biological detail.

**Solution:** Two-level hierarchy with explicit mapping:

```python
# Step 1: Define L3 → L2 mapping
L3_TO_L2_MAP = {
    # Fine-grained L3 labels → Coarse L2 labels
    "GC_B_Light_Zone_Centrocyte": "GC_B",
    "GC_B_Transitional": "GC_B",
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC_B",
    "GC_B": "GC_B",
    "Plasma_IgA": "Plasma",
    "Plasma_IgG": "Plasma",
    # ... etc
}

# Step 2: Create L2 labels from L3
adata.obs['Cell_Type_L2'] = adata.obs['cell_type_expert'].map(L3_TO_L2_MAP)

# Step 3: Train scANVI at L2 level for stability
scanvi_L2 = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key='Cell_Type_L2',
    unlabeled_category='Unknown'
)
scanvi_L2.train(max_epochs=200, ...)

# Step 4: Train separate scANVI at L3 level for granularity
scanvi_L3 = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key='cell_type_expert',
    unlabeled_category='Unknown'
)
scanvi_L3.train(max_epochs=200, ...)
```

**Key benefits:**
- L2 provides stable training with sufficient cells per type
- L3 captures biological subtypes for final annotation
- Explicit mapping ensures consistency

### 14. UMAP Operator Saving for Query Projection (v2.5.3+ NEW)

**Problem:** Query data needs to be projected into the same UMAP space as reference for visualization.

**Solution:** Save UMAP operator during reference training:

```python
import umap
import joblib

# After training, fit UMAP on reference latent
umap_operator = umap.UMAP(
    n_neighbors=30,
    min_dist=0.5,
    metric='euclidean',
    random_state=42
)
umap_operator.fit(adata.obsm['X_scanvi'])

# Save operator
joblib.dump(umap_operator, output_dir / "umap_operator.joblib")

# Save reference with UMAP for merging
adata.obsm['X_umap'] = umap_operator.transform(adata.obsm['X_scanvi'])
adata.write(output_dir / "reference_with_umap.h5ad")

# In query mapping step:
umap_op = joblib.load(model_dir / "umap_operator.joblib")
adata_query.obsm['X_umap'] = umap_op.transform(adata_query.obsm['X_scanvi'])
```

```python
def recover_counts_from_external(adata, external_file_path, counts_layer_name='counts'):
    """
    Recover counts layer by matching cell IDs from external h5ad file
    """
    adata_ext = sc.read_h5ad(external_file_path)

    # Match cells
    matching_cells = set(adata.obs.index) & set(adata_ext.obs.index)

    # Match genes (handle ENSEMBL vs symbols)
    gene_map = {}  # current_index -> external_index
    for i, curr_gene in enumerate(adata.var_names):
        for j, ext_gene in enumerate(adata_ext.var_names):
            if curr_gene == ext_gene:
                gene_map[i] = j
                break

    # Build recovered matrix
    from scipy import sparse
    recovered = sparse.lil_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)

    for cid in matching_cells:
        ci = list(adata.obs.index).index(cid)
        ei = list(adata_ext.obs.index).index(cid)

        counts_vec = adata_ext.layers['counts'][ei, :].toarray().flatten()
        for curr_g, ext_g in gene_map.items():
            recovered[ci, curr_g] = counts_vec[ext_g]

    adata.layers['counts'] = recovered.tocsr()
    return adata

# Usage
adata = recover_counts_from_external(
    adata,
    "/path/to/original/data_with_counts.h5ad"
)
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

### All Cells (v2.5 - PRODUCTION)
- **File:** `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5.py`
- **Features:** Complete covariates system, unknown cleaning, scArches-ready
- **New in v2.5:**
  - Fixed HVG workflow (correct execution order)
  - Simplified model control (RETRAIN_SCVI flags)
  - Safety mechanisms (HVG check, error handling)
  - Thread control for HPC stability
- **Model:** Human_Lung_Atlas.pkl
- **Parameters:** n_latent=100, n_layers=2, HVG=4000

### All Cells (v2.3.1 - STABLE)
- **File:** `allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py`
- **Features:** Complete covariates system, unknown cleaning, scArches-ready
- **New in v2.3.1:**
  - Continuous covariates (MT%, stress, cell cycle)
  - Neighborhood purity-based unknown filtering
  - Class imbalance handling (`n_samples_per_label`)
  - Expanded marker gene forcing

### Epithelial Cells (v2.7-PRODUCTION)
- **File:** `allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py`
- **Architecture:** Dual scANVI (major lineage + fine-grained)
- **Features:**
  - Hierarchical labeling (Basal_Lineage → Basal/Suprabasal/Dividing)
  - Whitelist filtering from downstream h5ad
  - Force-included epithelial markers (45+ genes)
- **Expected subtypes:** Basal, Suprabasal, Ciliated, Secretory (Goblet/Club), SMG types, Ionocytes, Brush

### T/NK Cells (v3.5.1 - PRODUCTION)
- **File:** `t_scvi_celltypist_scanvi_20260108_v3_5_1.py`
- **Model:** Immune_All_Low.pkl
- **P0 Fixes:** scANVI dedicated adata_model, predict() fix, majority voting integration
- **Expected subtypes:** CD4+ (naive/CM/EM/TRM), CD8+ (EM/EMRA/TRM), Tregs, NK, NKT, MAIT, ILC

### B Cells (v3.5.2 - PRODUCTION)
- **File:** `bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py`
- **Features:**
  - Counts recovery function (NEW in v3.5.2)
  - Immunoglobulin isotype markers (IgM/IgG/IgA/IgE)
  - Plasma cell maturation tracking
- **Expected subtypes:** Naive B, Memory B (IgM/IgG/IgA), Plasma cells, Plasmablasts, GC B cells, Bregs

### B Cells - L2 scArches Workflow (v2.5.3 - PRODUCTION)
- **Step 1:** `step1_train_bcell_L2_20260204_v2_5_3.py` - Train L2 reference with UMAP operator
- **Step 2:** `step2_map_query_20260204_v2_5_4.py` - Map query and merge with reference
- **Merged Training:** `merged_scanvi_training_20260208_v2_5_5.py` - Train unified model on merged data
- **Features:**
  - Two-level hierarchy (L2 coarse: GC_B/Plasma/Naive_B/Memory_B, L3 fine: subtypes)
  - UMAP operator saving for consistent query projection
  - Covariate matching between reference and query
  - Post-merge training for unified latent space

### Myeloid Cells (v3.5.1 - PRODUCTION)
- **File:** `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py`
- **Features:** All v3.5.x P0 fixes
- **Expected subtypes:**
  - Macrophages (alveolar, interstitial, CCL+, CHIT1+, CX3CR1+)
  - Monocytes (CD14+, CD16+)
  - DCs (DC1, DC2, pDC, activated)
  - Mast cells

### Stromal/Vascular Cells (v3.5.1 - PRODUCTION)
- **File:** `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py`
- **Target types:**
  - Endothelial (arterial/venous pulmonary/systemic, lymphatic, capillary)
  - Fibroblasts (adventitial, alveolar, peribronchial, myofibroblast)
  - Smooth muscle (airway, arterial, perivascular)
  - Pericytes, Mesothelial, Chondrocytes
- **Marker genes:** PECAM1/CDH5 (endothelial), COL1A1/PDGFRA (fibroblasts), ACTA2/MYH11 (SMC)

## Advanced Analysis Workflows

### Subcluster Analysis (Post-Integration)

**Purpose:** Fine-grained clustering within cell types while preserving global BBKNN batch correction.

**Key principle:** Use `restrict_to` parameter to leverage existing BBKNN graph without rebuilding neighbors.

**Files:**
- `epithelial_subcluster_bbknn_pipeline_20260120_v3_4.py` (v4.1-HOTFIX)
- `bcell_subcluster_analysis_v2_20260119.py` (v2.0-PRODUCTION)
- `myeloid_subcluster_analysis_v2_20260121.py`
- `stromal_subcluster_analysis_v2_20260121.py`

**Workflow:**
```python
# Load BBKNN-integrated data (BBKNN already run)
adata = sc.read_h5ad("epithelial_celltypist_filtered_final.h5ad")

# For each cell type, perform subclustering using global graph
for celltype in celltypes:
    # DO NOT rebuild neighbors - use restrict_to
    sc.tl.leiden(
        adata,
        restrict_to=('celltypist_pred', [celltype]),  # KEY PARAMETER
        resolution=0.2,  # Conservative to avoid over-splitting
        key_added=f'leiden_{celltype}_subcluster',
        neighbors_key='neighbors'  # Use existing BBKNN graph
    )

# Marker analysis with filtering
sc.tl.rank_genes_groups(
    adata,
    groupby=f'leiden_{celltype}_subcluster',
    method='wilcoxon',
    use_raw=True,  # Or False depending on data structure
    pts=True
)

# Filter technical artifacts
# - MT-* genes
# - Ribosomal (RPS/RPL)
# - Immediate early genes (FOS/JUN/EGR)
# - Dissociation stress (HSPA1A/HSPA1B)
# - Unannotated (LINC/AC/AL/RP)
```

**Critical design choices:**
- **Fixed low resolution** (0.2-0.3) - intentionally conservative to avoid batch-driven splits
- **No neighbor rebuild** - `restrict_to` preserves batch correction
- **Original UMAP** - visualize on existing UMAP coordinates
- **Marker filtering** - remove technical artifacts before interpretation

### Doublet Detection (CellTypist-based)

**Purpose:** Identify doublets based on multi-lineage signatures without removing all transitional states.

**File:** `celltypist_doublet_detection_20260114_v1_1.py` (v1.1-PRODUCTION)

**Method:** Top2 lineage + margin logic (more robust than simple thresholding)

**Workflow:**
```python
# Run CellTypist annotation
predictions = celltypist.annotate(
    adata,
    model=CELLTYPIST_MODEL,
    majority_voting=True
)

# Calculate top2 lineage scores
lineage_probs = predictions.probability_matrix.groupby(LINEAGE_GROUPS, axis=1).sum()
top1_prob = lineage_probs.max(axis=1)
top2_prob = lineage_probs.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
margin = top1_prob - top2_prob

# Identify doublets (high top2 + low margin)
is_doublet = (top2_prob > 0.5) & (margin < 0.15)

# Multi-evidence gating (optional)
# - Neighborhood consistency (cells surrounded by same lineage less likely doublets)
# - Mapping confidence (high-confidence assignments more reliable)
# - Marker gene expression (validate with known lineage markers)

adata.obs['doublet_score'] = top2_prob
adata.obs['lineage_margin'] = margin
adata.obs['is_doublet_candidate'] = is_doublet
```

**Lineage groups** (example for lung data):
- Epithelial (Basal, Ciliated, Secretory, AT1/AT2)
- T_NK (CD4+, CD8+, NK, NKT, MAIT)
- B_Plasma (B naive/memory, Plasma)
- Myeloid (Macrophages, DCs, Monocytes, Mast)
- Stromal (Fibroblasts, SMC, Pericytes)
- Endothelial (Arterial, Venous, Lymphatic, Capillary)

**When to use:**
- High cell density regions with mixed signals
- Epithelial cells with immune markers (e.g., Tier C in nasal data)
- Post-dissociation datasets with high doublet rates

### cNMF Integration (Gene Expression Programs)

**Purpose:** Identify gene expression programs (GEPs) and map to cell states/clusters.

**Files:**
- `label_guided_cnmf_pipeline_20260114_v1_1.py` (v1.1-HOTFIX)
- `cnmf_results_analysis_20260114_v1_1.py`
- `complete_epithelial_analysis_pipeline_20260121_v2_1.py` (BBKNN + Harmony + cNMF)

**Workflow:**
```python
# Step 1: Gene filtering (CRITICAL for biological GEPs)
FILTER_GENES = {
    'remove_mt': True,           # MT-* genes
    'remove_ribo': True,         # RPS/RPL/MRPS/MRPL
    'remove_histone': True,      # H1/H2A/H2B/H3/H4
    'remove_pseudogenes': True,  # *P* genes
    'remove_ensg': True,         # ENSG* IDs
    'remove_unannotated': True   # AC/AL/RP/CTD/LINC
}

# Step 2: Run cNMF per cell type
from cnmf import cNMF

for celltype in celltypes:
    adata_subset = adata[adata.obs['cell_type'] == celltype].copy()

    # Confidence filtering (optional but recommended)
    if USE_CONFIDENCE_FILTER:
        adata_subset = adata_subset[
            adata_subset.obs['scanvi_confidence'] > 0.5
        ].copy()

    # Initialize cNMF
    cnmf_obj = cNMF(
        output_dir=f"cnmf_{celltype}",
        name=f"{celltype}_cnmf"
    )

    # Prepare counts (use raw integer counts)
    cnmf_obj.prepare(
        counts_fn=adata_subset,
        components=np.arange(10, 55, 5),  # K range: 10-50
        n_iter=100,
        seed=14,
        num_hvg=2000
    )

    # Factorize (can use multi-processing)
    cnmf_obj.factorize(
        worker_i=0,
        total_workers=16  # Multi-processing speedup
    )

    # Combine and select K
    cnmf_obj.combine()
    cnmf_obj.k_selection_plot()  # Choose optimal K

    # Get GEP loadings
    usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=25)

# Step 3: Map GEPs to clusters
# - Correlate GEP usage with cluster centroids
# - Identify GEP-dominant clusters
# - Validate with known marker genes
```

**Complete pipeline integration:**
```bash
# Full workflow: BBKNN + Harmony + cNMF + mapping
python complete_epithelial_analysis_pipeline_20260121_v2_1.py

# Pipeline steps:
# 1. BBKNN batch integration (5-10 min)
# 2. Harmony batch correction (10-20 min)
# 3. cNMF on uncorrected data (2-6 hours with multi-processing)
# 4. cNMF on Harmony-corrected data (2-6 hours)
# 5. GEP-cluster mapping and comparison (10-20 min)
```

**Key parameters:**
- **K range:** 10-50 for large cell types (>10k cells), 10-25 for small (<5k)
- **n_iter:** 100 (balance speed vs. stability)
- **num_hvg:** 2000 (after technical gene filtering)
- **Multi-processing:** 4-16 workers for 4-8x speedup

### scArches Reference Mapping (Query Projection)

**Purpose:** Map new query data onto a pre-trained reference model for annotation transfer.

**Files:**
- `step1_train_bcell_L2_20260204_v2_5_3.py` (v2.5.3-PRODUCTION) - Reference training
- `step2_map_query_20260204_v2_5_4.py` (v2.5.3-PRODUCTION) - Query mapping + merge
- `scarches_mapping_20260127_v1_2_1.py` (v1.2.1-PRODUCTION) - Generic scArches mapping
- `merged_scanvi_training_20260208_v2_5_5.py` (v2.5.5-HOTFIX) - Post-merge training

**Workflow (2-Step Pattern):**

```python
# Step 1: Train Reference Model with UMAP operator
# ================================================
# - Train scVI on reference data
# - Train scANVI on reference labels
# - Fit and save UMAP operator for query projection
# - Save reference with UMAP coordinates

import umap
import joblib

# After scANVI training
umap_op = umap.UMAP(n_neighbors=30, min_dist=0.5, random_state=42)
umap_op.fit(adata_ref.obsm['X_scanvi'])
joblib.dump(umap_op, "umap_operator.joblib")
adata_ref.obsm['X_umap'] = umap_op.transform(adata_ref.obsm['X_scanvi'])
adata_ref.write("reference_with_umap.h5ad")

# Step 2: Map Query Data
# ======================
# - Load pretrained scVI/scANVI models
# - Calculate covariates (MUST match reference training)
# - Use scVI.prepare_query_anndata() for gene alignment
# - Load query into reference model with load_query_data()
# - Fine-tune model on query batches
# - Predict labels with confidence scores
# - Project UMAP using saved operator

from scvi.model import SCANVI

# Prepare query (gene alignment)
scvi.model.SCVI.prepare_query_anndata(adata_query, ref_scvi_model_dir)

# Load into reference
scanvi_model = SCANVI.load_query_data(
    adata_query,
    ref_scanvi_model_dir,
    freeze_classifier=True  # Keep reference classification
)

# Fine-tune
scanvi_model.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0},
    check_val_every_n_epoch=10
)

# Predict
predictions = scanvi_model.predict(adata_query)
proba = scanvi_model.predict(adata_query, soft=True)
confidence = np.asarray(proba).max(axis=1)

# UMAP projection using saved operator
umap_op = joblib.load("umap_operator.joblib")
adata_query.obsm['X_umap'] = umap_op.transform(adata_query.obsm['X_scanvi'])

# Step 3: Merge Reference + Query
# ===============================
# - Concatenate reference and query AnnData objects
# - Align on HVG intersection in reference order
# - Preserve reference UMAP and add query cells
# - Retrain scANVI on merged data for unified latent space

# Merge on HVG intersection
hvg_ref = pd.read_csv(hvg_file, header=None)[0].tolist()
common_genes = list(set(adata_ref.var_names) & set(adata_query.var_names))
hvg_order = [g for g in hvg_ref if g in common_genes]  # Reference order

adata_ref_subset = adata_ref[:, hvg_order].copy()
adata_query_subset = adata_query[:, hvg_order].copy()
adata_merged = ad.concat([adata_ref_subset, adata_query_subset])

# Retrain on merged data for unified space
scanvi_merged = SCANVI.from_scvi_model(
    scvi_merged,
    labels_key='Cell_Type_L2_train',
    unlabeled_category='Unknown'
)
scanvi_merged.train(max_epochs=200, ...)
```

**Key Configuration Parameters:**
```python
# Covariates (MUST match between reference and query)
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"
CONTINUOUS_COVARIATES = ['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']

# Required for prepare_query_anndata()
scvi.model.SCVI.setup_anndata(
    adata_ref,
    layer='counts',
    batch_key=BATCH_KEY,
    categorical_covariate_keys=[TISSUE_KEY],  # If used
    continuous_covariate_keys=CONTINUOUS_COVARIATES
)
```

**Critical Success Factors:**
1. **Covariate Consistency** - Query MUST have same covariates as reference
2. **Gene Name Compatibility** - Use `prepare_query_anndata()` for alignment
3. **UMAP Operator** - Save during reference training for consistent projection
4. **Fine-tuning** - Use `freeze_classifier=True` to preserve reference annotations
5. **Confidence Threshold** - Filter low-confidence predictions (<0.5) post-mapping

**Common Errors and Solutions:**
- **"Category not in categories"** → Pre-add 'unknown_tissue' to tissue categories
- **"NaN in predictions"** → Force-overwrite result columns in metadata restore step
- **"Gene mismatch"** → Ensure HVG intersection preserves reference gene order

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

### Issue: "Covariates not reducing technical variation"
**Solution:** Check covariate calculation and verify model setup:
```python
# Ensure covariates are numeric and not NaN
print(adata.obs[['pct_counts_mt', 'stress_score']].describe())

# Verify model setup
scvi.model.SCVI.setup_anndata(
    adata_model,
    layer='counts',
    batch_key='dataset',
    continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
)

# Check that covariates were registered
print(adata_model.uns['_scvi']['extra_continuous_keys'])
```

### Issue: "scANVI predictions inconsistent with CellTypist"
**Solution:** Use hierarchical labeling or adjust confidence threshold:
```python
# Option 1: Start with major lineages (more stable)
adata.obs['major_lineage'] = adata.obs['cell_type_celltypist'].map(MAJOR_LINEAGE_MAP)
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key='major_lineage',  # Not fine-grained labels
    unlabeled_category='Unknown'
)

# Option 2: Filter low-confidence CellTypist predictions before scANVI
confident_mask = adata.obs['celltypist_conf_max'] > 0.6
adata.obs['labels_for_scanvi'] = adata.obs['cell_type_celltypist'].copy()
adata.obs.loc[~confident_mask, 'labels_for_scanvi'] = 'Unknown'
```

### Issue: "Subcluster analysis reintroduces batch effects"
**Solution:** Use `restrict_to` parameter, DO NOT rebuild neighbors:
```python
# WRONG - rebuilds neighbors, loses batch correction
adata_subset = adata[adata.obs['cell_type'] == 'Basal'].copy()
sc.pp.neighbors(adata_subset)  # BAD
sc.tl.leiden(adata_subset)

# CORRECT - uses existing BBKNN graph
sc.tl.leiden(
    adata,  # Full dataset, not subset
    restrict_to=('cell_type', ['Basal']),  # Restrict clustering only
    resolution=0.2,
    key_added='leiden_basal_subcluster',
    neighbors_key='neighbors'  # Existing BBKNN graph
)
```

### Issue: "Doublet detection removes transitional states"
**Solution:** Use multi-evidence gating instead of hard threshold:
```python
# Calculate multiple lines of evidence
lineage_probs = predictions.probability_matrix.groupby(LINEAGE_GROUPS, axis=1).sum()
top2_prob = lineage_probs.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
margin = lineage_probs.max(axis=1) - top2_prob

# Calculate neighborhood purity
from sklearn.neighbors import NearestNeighbors
nn = NearestNeighbors(n_neighbors=30)
nn.fit(adata.obsm['X_scanvi'])
_, indices = nn.kneighbors(adata.obsm['X_scanvi'])
purity = [
    (adata.obs['predicted_lineage'].iloc[neighbors] ==
     adata.obs['predicted_lineage'].iloc[i]).mean()
    for i, neighbors in enumerate(indices)
]

# Multi-evidence scoring (requires multiple red flags)
adata.obs['doublet_score'] = (
    (top2_prob > 0.5).astype(int) +          # High second lineage
    (margin < 0.15).astype(int) +            # Small margin
    (np.array(purity) < 0.5).astype(int)     # Low neighborhood purity
)

# Only remove cells with 2+ red flags
adata = adata[adata.obs['doublet_score'] < 2].copy()
```

### Issue: "cNMF factorization fails or produces poor GEPs"
**Solution:** Check gene filtering and ensure counts are integers:
```python
# Verify counts are integers
print(f"Counts dtype: {adata.layers['counts'].dtype}")
print(f"Counts range: {adata.layers['counts'].min():.2f} - {adata.layers['counts'].max():.2f}")

# Apply comprehensive gene filtering
filter_genes = (
    ~adata.var_names.str.startswith('MT-') &
    ~adata.var_names.str.match(r'^RP[SL]') &
    ~adata.var_names.str.match(r'^MRP[SL]') &
    ~adata.var_names.str.match(r'^H[1234]') &
    ~adata.var_names.str.match(r'^HIST') &
    ~adata.var_names.str.match(r'ENSG\d+') &
    ~adata.var_names.str.match(r'^(LINC|AC\d+|AL\d+|RP11-|CTD-|CTB-)')
)
adata_for_cnmf = adata[:, filter_genes].copy()

# Use appropriate K range based on cell count
n_cells = adata_for_cnmf.n_obs
if n_cells > 10000:
    k_range = [25, 30, 35, 40, 45, 50]
elif n_cells > 5000:
    k_range = [15, 20, 25, 30, 35]
else:
    k_range = [10, 15, 20, 25]
```

### Issue: "Counts recovery fails with gene name mismatch"
**Solution:** Use both ENSEMBL IDs and symbols for matching:
```python
# Build flexible gene matching
if 'symbol_base' in adata.var.columns and 'symbol_base' in adata_ext.var.columns:
    # Match by symbol_base first
    curr_symbols = adata.var['symbol_base'].to_dict()
    ext_symbols = adata_ext.var['symbol_base'].to_dict()
    gene_map = {}
    for i, sym in curr_symbols.items():
        for j, ext_sym in ext_symbols.items():
            if sym == ext_sym:
                gene_map[i] = j
                break

# Fallback: match by var_names (may be ENSEMBL IDs)
for i in range(adata.n_vars):
    if i not in gene_map:
        curr_id = adata.var_names[i]
        if curr_id in adata_ext.var_names:
            gene_map[i] = list(adata_ext.var_names).index(curr_id)

print(f"Matched {len(gene_map)}/{adata.n_vars} genes ({len(gene_map)/adata.n_vars*100:.1f}%)")
```

## Critical Fixes in Modern Pipelines

### v2.x Series (All-cells baseline)
All pipelines (v2.0+) include these MANDATORY fixes:

1. **Gene alignment** - `normalize_gene_names()` BEFORE `.raw` creation
2. **UMAP separation** - Use `neighbors_key` to prevent overwriting
3. **CPU-safe CUDA** - Conditional `torch.cuda.manual_seed_all()`
4. **Filter on counts** - Gene filtering uses `layers['counts']`, not `.X`
5. **Index-aligned writes** - Use `.reindex()` when transferring predictions
6. **Single GPU control** - Explicit `devices=1` for scVI/scANVI
7. **Rare type handling** - Merge types with <10 cells to "Unknown"

### v2.3.1 Series (Covariates + Unknown Cleaning)
Additional features in all-cells v2.3.1:

8. **Covariates system** - MT%, stress, cell cycle as continuous covariates
9. **Unknown cleaning** - Neighborhood purity + quality gating
10. **Class imbalance** - `n_samples_per_label` for stable training
11. **Unified n_latent** - scVI and scANVI use same latent dimension (100)
12. **scArches parameters** - `encode_covariates=True`, `use_layer_norm="both"`, `n_layers=2`

### v2.5.x Series (scArches + HPC Stability)
Production enhancements for reference mapping and cluster environments (v2.5+):

20. **HPC Thread Control** - Explicit `OMP_NUM_THREADS=1` prevents deadlocks
21. **scVI dl_num_workers=0** - Avoids DataLoader multiprocessing issues
22. **UMAP Operator Saving** - `joblib.dump()` for query projection consistency
23. **L2/L3 Hierarchy** - Two-level labeling with explicit mapping dictionaries
24. **Reference+Query Merge** - HVG intersection in reference gene order
25. **Model Control Flags** - `RETRAIN_SCVI`, `RETRAIN_SCANVI` for workflow management
26. **Robust Integer Validation** - Float32 tolerance for count verification

### v3.5.x Series (P0 Fixes for Cell-type Pipelines)
Critical bug fixes in cell-type-specific pipelines (v3.5.1+):

13. **P0-1: scANVI dedicated adata_model** - Labels alignment fix (prevents index mismatch)
14. **P0-2: predict() method fix** - `predict(soft=True).max()` → `np.asarray()` (dtype safety)
15. **P0-3: Majority voting integration** - CellTypist `majority_voting` correctly merged
16. **P0-4: Feature validation** - Hard check for CellTypist gene overlap (prevent silent failures)
17. **P0-5: Pretrained HVG uses var_names** - Not `symbol_base` (fixes gene list loading)
18. **P1-6: Memory optimization** - No full X copy, use layers directly
19. **P1-7: Minimal adata_model** - Only essential obs/var for training (60-70% memory reduction)

## Dependencies

**Core packages:**
```
scanpy>=1.9
scvi-tools>=1.0
celltypist>=1.6
torch>=2.0 (with CUDA support)
mygene>=3.2
numpy>=1.21
pandas>=1.3
scipy>=1.7
scikit-learn>=1.0
```

**Optional but recommended:**
```
bbknn>=1.5 (for BBKNN preprocessing)
harmony-pytorch (alternative batch correction)
cnmf>=1.4 (for gene expression program analysis)
scArches (for model transfer learning)
```

**System requirements:**
```
GPU: NVIDIA GPU with CUDA support (recommended: V100/A100)
RAM: 64-256GB depending on dataset size
CPU: 16-48 cores for multi-processing
Storage: 500GB+ for large datasets and intermediate files
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
- scArches query mapping (50k cells): ~15-30 min
- Merged training (ref+query, 150k cells): ~45-90 min

## Code Quality Standards

When modifying pipelines:

### Version Control
1. Preserve version history in docstrings (list critical fixes)
2. Use date-stamped filenames: `{celltype}_{analysis}_{YYYYMMDD}_v{X}_{Y}.py`
3. Document major changes in P0/P1/P2 priority levels
4. Keep old versions for reproducibility (don't overwrite)

### Layer and Data Specifications
5. Use explicit layer specifications (`layer='counts'` not `layer=None` with assumptions)
6. Always verify data types before training (counts must be integers)
7. Check for `layers['counts']` existence before scVI training
8. Document which representation each analysis step uses (`.X`, `.raw.X`, `layers['counts']`)

### Model Management
9. Always save HVG lists when training new models (`hvg_genes.txt`)
10. Save model architecture parameters in `model_config.json`
11. Use descriptive model paths: `{date}/{celltype}/scvi_model/`
12. Document scArches-compatible parameters in model metadata

### Quality Control
13. Add sanity checks (gene count validation, confidence score ranges)
14. Export statistics to CSV for audit trails
15. Use `adata.uns['pipeline_info']` to document parameters
16. Save intermediate checkpoints for large datasets (>200k cells)

### Code Organization
17. Separate configuration section at top (UPPERCASE variables)
18. Use helper functions for repeated operations (gene filtering, covariate calculation)
19. Add progress bars for long-running operations
20. Clear memory after large operations (`gc.collect()`)
21. Thread control at module import (HPC stability)

### Reproducibility
22. Set random seeds at start (numpy, torch, scanpy, scvi)
23. Document exact package versions used (`pip freeze > requirements.txt`)
24. Save command-line arguments if using argparse
25. Log all major parameters to output file
26. Save UMAP operator with reference model (scArches)

### Performance Optimization
27. Use HVG optimization for datasets >200k cells
28. Leverage multi-processing where available (cNMF, BBKNN)
29. Profile memory usage for large operations
30. Use sparse matrices throughout (avoid `.toarray()` unless necessary)
31. Use `dl_num_workers=0` for scVI stability on HPC

### New Best Practices (v2.3.1+)
29. **Always calculate covariates** (MT%, stress, cell cycle) for scVI/scANVI
30. **Filter Unknown cells** using neighborhood purity before downstream analysis
31. **Use hierarchical labeling** for complex tissues (major lineage → fine-grained)
32. **Balance class imbalance** with `n_samples_per_label` in scANVI training
33. **Validate gene overlap** before CellTypist annotation (>70% recommended)
34. **Preserve BBKNN graphs** in subcluster analysis (use `restrict_to`)
35. **Filter technical genes** before cNMF (MT/Ribo/Histone/IEG/Unannotated)
36. **Multi-evidence doublet detection** (don't rely on single metric)
37. **HPC thread control** - Set `OMP_NUM_THREADS=1` and `dl_num_workers=0`
38. **Save UMAP operator** during reference training for query projection
39. **Use L2/L3 hierarchy** - Coarse for training stability, fine for final annotation
40. **scArches preparation** - Match covariates exactly between reference and query
41. **Model control flags** - Use RETRAIN_SCVI/RETRAIN_SCANVI for workflow management
42. **Merged training** - Retrain on reference+query for unified latent space

### Anti-patterns to Avoid
37. ❌ Don't subset data before normalizing (breaks cross-sample comparisons)
38. ❌ Don't rebuild neighbors graph after BBKNN (loses batch correction)
39. ❌ Don't use `.X` for scVI when `layers['counts']` exists
40. ❌ Don't ignore low CellTypist overlap (<50% genes matched)
41. ❌ Don't skip gene name normalization before CellTypist
42. ❌ Don't use position-based indexing for predictions (`.values` without `.reindex()`)
43. ❌ Don't train scANVI on rare types (<10 cells) without merging
44. ❌ Don't use high resolution (>0.5) for subcluster analysis (causes over-splitting)
45. ❌ Don't forget to set `dl_num_workers=0` on HPC clusters
46. ❌ Don't skip covariate calculation in scArches query mapping
47. ❌ Don't project query UMAP without using saved reference operator
48. ❌ Don't use default threading on HPC (causes deadlocks)
49. ❌ Don't merge reference+query without HVG intersection in reference order
