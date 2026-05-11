#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Epithelial Cell scVI-scANVI Training Pipeline (PRODUCTION v2.3)
With Hierarchical Label Consolidation

New in v2.3:
- Hierarchical label consolidation (13 fine types → 6 major lineages)
- Coarse-to-fine strategy for robust scANVI training
- Enhanced visualization comparing major lineages vs fine types
- Complete annotation workflow for downstream analysis

Based on CellTypist Annotations

Author: r2end
Date: 2025-01-10
Version: 2.3 - Production with hierarchical labels
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from scipy import sparse
from sklearn.metrics import confusion_matrix

warnings.filterwarnings('ignore')

# Set plotting style
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

print("="*80)
print("Epithelial scVI-scANVI Pipeline (PRODUCTION v2.3)")
print("With Hierarchical Label Consolidation")
print("="*80)
print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"scanpy version: {sc.__version__}")
print(f"scvi-tools version: {scvi.__version__}")
print(f"PyTorch version: {torch.__version__}")
print(f"Device: {'GPU' if torch.cuda.is_available() else 'CPU'}")
print("="*80)

# ===== CONFIGURATION =====

# Input/Output paths
INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_celltypist_filtered_final.h5ad"
OUTPUT_DIR = (f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/celltypist_epithelial")
CHECKPOINT_DIR = f"{OUTPUT_DIR}/checkpoints"

# Create directories
for dir_path in [OUTPUT_DIR, f"{OUTPUT_DIR}/figures", f"{OUTPUT_DIR}/scvi_model", 
                 f"{OUTPUT_DIR}/scanvi_model", CHECKPOINT_DIR]:
    os.makedirs(dir_path, exist_ok=True)

# Key column names
BATCH_KEY = 'sample'  # Options: 'dataset', 'Sample', 'batch', etc.
                       # IMPORTANT: Use 'Sample' for per-sample batch correction
                       # Avoid using 'study' if study = biological difference
CELLTYPE_KEY = 'celltypist_pred'
UNLABELED_CATEGORY = 'Unknown'

# ===== Hierarchical Label Configuration =====
USE_HIERARCHICAL_LABELS = True  # ⭐ Enable hierarchical labels
VISUALIZE_HIERARCHY = True      # Generate hierarchy visualizations

# Hierarchical label mapping (complete epithelial cell types)
MAJOR_LINEAGE_MAP = {
    # Alveolar lineage (lung-specific)
    'AT1': 'Alveolar',
    'AT1 ': 'Alveolar',  # Handle trailing space
    'AT2': 'Alveolar',
    
    # Basal lineage (stem/progenitor cells)
    'Basal': 'Basal_Lineage',
    'Suprabasal': 'Basal_Lineage',
    'SMG_Basal': 'Basal_Lineage',
    'Dividing_Basal': 'Basal_Lineage',  # ⭐ Proliferating basal cells
    
    # Ciliated lineage (including precursors)
    'Ciliated': 'Ciliated_Lineage',
    'Deuterosome': 'Ciliated_Lineage',
    'Deuterosomal': 'Ciliated_Lineage',  # ⭐ Alternative spelling
    
    # Secretory lineage (mucus-producing cells)
    'Secretory_Goblet': 'Secretory_Lineage',
    'Secretory_Club': 'Secretory_Lineage',  # ⭐ Club cells
    'SMG_Mucous': 'Secretory_Lineage',
    'SMG_Serous': 'Secretory_Lineage',
    'SCGB1A1+': 'Secretory_Lineage',  # Alternative club cell annotation
    
    # Duct cells
    'SMG_Duct': 'Duct',
    
    # Rare specialized cells
    'Ionocyte_n_Brush': 'Rare_Specialized',
    'Ionocyte': 'Rare_Specialized',
    'Brush': 'Rare_Specialized'
}

# ===== Force-Include Marker Genes (P1 Strategy) =====
# These critical marker genes will be FORCED into the model training
# even if they're not in top HVG, preventing AT2/rare cell type "collapse"
FORCE_INCLUDE_MARKERS = [
    # Alveolar markers (AT1/AT2 separation)
    'SFTPC', 'SFTPB', 'ABCA3', 'SLC34A2', 'ETV5',  # AT2
    'AGER', 'PDPN', 'CLDN18', 'CAV1',  # AT1
    
    # Basal markers
    'TP63', 'KRT5', 'KRT14', 'KRT15',
    
    # Suprabasal markers
    'KRT13', 'KRT4',
    
    # Ciliated markers
    'FOXJ1', 'RSPH1', 'PIFO', 'DNAH5',
    
    # Deuterosome markers
    'DEUP1', 'CCNO',
    
    # Secretory markers
    'MUC5AC', 'MUC5B', 'SCGB1A1', 'SPDEF',
    'BPIFA1', 'WFDC2',  # Club cells
    
    # Rare cells
    'CFTR', 'FOXI1',  # Ionocyte
    'DCLK1', 'TRPM5',  # Brush
    
    # SMG markers
    'LYZ', 'LTF',  # Serous
    'KRT7', 'KRT19'  # Duct
]

# ===== Stratified HVG Selection (Advanced Option) =====
USE_STRATIFIED_HVG = False  # ⭐ Enable lineage-aware HVG selection
# When True: Select HVG within each major lineage, then union
# Prevents rare lineages from being "drowned out" by abundant ones
# Recommended if you have very imbalanced cell types (e.g., AT2 << Basal)

# Color schemes
MAJOR_LINEAGE_COLORS = {
    'Alveolar': '#E74C3C',
    'Basal_Lineage': '#3498DB',
    'Ciliated_Lineage': '#2ECC71',
    'Secretory_Lineage': '#F39C12',
    'Duct': '#9B59B6',
    'Rare_Specialized': '#95A5A6'
}

FINE_TYPE_COLORS = {
    'AT1': '#E74C3C', 'AT2': '#C0392B',
    'Basal': '#3498DB', 'Suprabasal': '#5DADE2', 'SMG_Basal': '#2874A6',
    'Ciliated': '#2ECC71', 'Deuterosome': '#58D68D',
    'Secretory_Goblet': '#F39C12', 'SMG_Mucous': '#F8C471', 'SMG_Serous': '#D68910',
    'SMG_Duct': '#9B59B6',
    'Ionocyte_n_Brush': '#95A5A6'
}

# scVI parameters
N_HVG = 4000  # More genes for complex epithelial heterogeneity
N_LATENT_SCVI = 100  # Higher for 270k cells
N_LAYERS = 3
N_HIDDEN = 256
DROPOUT_RATE = 0.1
GENE_LIKELIHOOD = "nb"
DISPERSION = "gene-batch"  # Better for multi-dataset integration

# scANVI parameters
N_LATENT_SCANVI = 75

# Training parameters
SCVI_MAX_EPOCHS = 400
SCANVI_MAX_EPOCHS = 200
EARLY_STOPPING = True
BATCH_SIZE = 1024  # Larger batch for 270k cells

# QC thresholds
MIN_CELLS_PER_BATCH = 10
LOW_CONFIDENCE_THRESHOLD = 0.3
MAX_LOW_CONFIDENCE_PCT = 80

# Options
SAVE_CHECKPOINTS = True
MARK_LOW_CONFIDENCE_AS_UNKNOWN = True

# Random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

# Device
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

print("\n" + "="*80)
print("CONFIGURATION")
print("="*80)
print(f"Input file: {INPUT_H5AD}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Device: {DEVICE}")
print(f"\n⭐ Hierarchical Labels: {USE_HIERARCHICAL_LABELS}")
if USE_HIERARCHICAL_LABELS:
    print(f"  Fine types: {len(MAJOR_LINEAGE_MAP)}")
    print(f"  Major lineages: {len(set(MAJOR_LINEAGE_MAP.values()))}")
print(f"\nBatch key: {BATCH_KEY}")
print(f"Min cells per batch: {MIN_CELLS_PER_BATCH}")
print(f"Low confidence threshold: {LOW_CONFIDENCE_THRESHOLD}")
print(f"\nscVI parameters:")
print(f"  HVG: {N_HVG}, Latent: {N_LATENT_SCVI}, Batch size: {BATCH_SIZE}")
print(f"  Dispersion: {DISPERSION}")
print(f"\nscANVI parameters:")
print(f"  Latent: {N_LATENT_SCANVI}, Max epochs: {SCANVI_MAX_EPOCHS}")
print("="*80)


# ===== STEP 1: Load and Validate Data =====

print("\n" + "="*80)
print("STEP 1: Loading and Validating Data")
print("="*80)

adata = sc.read_h5ad(INPUT_H5AD)

print(f"\n📊 Data Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

# Ensure unique cell IDs
if not adata.obs_names.is_unique:
    print("\n⚠️  Duplicate cell IDs found. Making unique...")
    adata.obs_names_make_unique()
    print("✓ Cell IDs are now unique")

# Check required columns
required_columns = [CELLTYPE_KEY, 'celltypist_conf']
missing_columns = [col for col in required_columns if col not in adata.obs.columns]
if missing_columns:
    raise ValueError(f"Missing required columns: {missing_columns}")

# Check batch key
if BATCH_KEY not in adata.obs.columns:
    print(f"\n⚠️  WARNING: '{BATCH_KEY}' not found")
    alternative_keys = ['batch', 'Sample', 'sample']
    for alt_key in alternative_keys:
        if alt_key in adata.obs.columns:
            BATCH_KEY = alt_key
            print(f"   Using '{alt_key}' as batch key")
            break
    else:
        print("   Creating dummy batch")
        adata.obs[BATCH_KEY] = 'batch1'

# Filter tiny batches
print("\n🔍 Checking batch sizes...")
batch_counts = adata.obs[BATCH_KEY].value_counts()
tiny_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index

if len(tiny_batches) > 0:
    print(f"\n⚠️  Found {len(tiny_batches)} batch(es) with < {MIN_CELLS_PER_BATCH} cells")
    for batch in tiny_batches:
        print(f"   - {batch}: {batch_counts[batch]} cells")
    
    n_before = adata.n_obs
    adata = adata[~adata.obs[BATCH_KEY].isin(tiny_batches)].copy()
    n_after = adata.n_obs
    
    print(f"\n✓ Removed {n_before - n_after:,} cells")
    print(f"✓ Remaining cells: {n_after:,}")
    
    if pd.api.types.is_categorical_dtype(adata.obs[BATCH_KEY]):
        adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].cat.remove_unused_categories()
else:
    print(f"✓ All batches have ≥ {MIN_CELLS_PER_BATCH} cells")

n_batches = adata.obs[BATCH_KEY].nunique()
n_celltypes = adata.obs[CELLTYPE_KEY].nunique()

print(f"\n  Batches: {n_batches}")
print(f"  Cell types: {n_celltypes}")

print(f"\n📋 Cell type distribution (top 10):")
for celltype, count in adata.obs[CELLTYPE_KEY].value_counts().head(10).items():
    pct = count / adata.n_obs * 100
    print(f"  {str(celltype):35s}: {count:7,} ({pct:5.2f}%)")

# Check if Unknown exists
if UNLABELED_CATEGORY in adata.obs[CELLTYPE_KEY].unique():
    print(f"\n⚠️  '{UNLABELED_CATEGORY}' exists in CellTypist labels")
    UNLABELED_CATEGORY = 'Unlabeled_scANVI'
    print(f"   Using '{UNLABELED_CATEGORY}' instead")

# Confidence statistics
avg_conf = adata.obs['celltypist_conf'].mean()
median_conf = adata.obs['celltypist_conf'].median()
low_conf_count = (adata.obs['celltypist_conf'] < 0.5).sum()
low_conf_pct = low_conf_count / adata.n_obs * 100

print(f"\n📊 Confidence statistics:")
print(f"  Mean: {avg_conf:.3f}")
print(f"  Median: {median_conf:.3f}")
print(f"  Low confidence (<0.5): {low_conf_count:,} ({low_conf_pct:.1f}%)")

very_low_conf_count = (adata.obs['celltypist_conf'] < LOW_CONFIDENCE_THRESHOLD).sum()
very_low_conf_pct = very_low_conf_count / adata.n_obs * 100

if very_low_conf_pct > MAX_LOW_CONFIDENCE_PCT:
    print(f"\n⚠️  WARNING: {very_low_conf_pct:.1f}% have confidence < {LOW_CONFIDENCE_THRESHOLD}")
    print(f"   scANVI will be mostly unsupervised")

print("\n✓ Data loaded and validated")


# ===== STEP 2: Data Preprocessing with Strict Validation =====

print("\n" + "="*80)
print("STEP 2: Data Preprocessing with Strict Validation")
print("="*80)

print("\n2.1 Validating data structure...")

# Strict counts layer validation
if 'counts' not in adata.layers:
    print("\n❌ ERROR: 'counts' layer not found!")
    raise ValueError("Missing 'counts' layer. Cannot proceed without raw count data.")
else:
    print("✓ 'counts' layer found")
    
    counts_min = adata.layers['counts'].min()
    counts_max = adata.layers['counts'].max()
    
    print(f"   Range: [{counts_min:.2f}, {counts_max:.2f}]")
    
    if counts_min < 0:
        raise ValueError(f"counts layer contains negative values (min={counts_min})")
    
    if counts_max < 50:
        print(f"   ⚠️  WARNING: Maximum count unusually low ({counts_max:.1f})")
    
    if not sparse.issparse(adata.layers['counts']):
        print(f"   Converting counts to sparse matrix...")
        adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])
        print(f"   ✓ Converted to sparse")

# Ensure .X is log-normalized
if 'log1p' in adata.layers:
    print("\n✓ Using existing 'log1p' layer for .X")
    adata.X = adata.layers['log1p'].copy()
else:
    print("\nComputing log-normalization...")
    adata.X = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()
    print("✓ Log-normalization complete")

print(f"\n✓ Data structure validated")

# Memory check
import psutil
mem_gb = psutil.Process(os.getpid()).memory_info().rss / 1e9
print(f"\n💾 Current memory usage: {mem_gb:.1f} GB")

# Save checkpoint
if SAVE_CHECKPOINTS:
    checkpoint_path = f"{CHECKPOINT_DIR}/01_preprocessed.h5ad"
    adata.write_h5ad(checkpoint_path, compression='gzip')
    print(f"💾 Checkpoint saved: {checkpoint_path}")


# ===== STEP 3: Highly Variable Genes Selection =====

print("\n" + "="*80)
print("STEP 3: Highly Variable Genes Selection (P1 Strategy)")
print("="*80)

if USE_STRATIFIED_HVG:
    print(f"\n⭐ Using Stratified HVG Selection (lineage-aware)")
    print(f"   Prevents rare lineages from being drowned out")
    
    # Temporarily assign major lineages if not done yet
    if 'major_lineage' not in adata.obs.columns:
        adata.obs['temp_major_lineage'] = adata.obs[CELLTYPE_KEY].map(MAJOR_LINEAGE_MAP)
        lineage_key = 'temp_major_lineage'
    else:
        lineage_key = 'major_lineage'
    
    # Select HVG within each major lineage, then union
    all_hvg_genes = set()
    n_genes_per_lineage = max(500, N_HVG // adata.obs[lineage_key].nunique())
    
    for lineage in adata.obs[lineage_key].unique():
        if pd.isna(lineage):
            continue
        
        print(f"\n   Processing {lineage}...")
        adata_lineage = adata[adata.obs[lineage_key] == lineage].copy()
        
        try:
            sc.pp.highly_variable_genes(
                adata_lineage,
                layer='counts',
                n_top_genes=n_genes_per_lineage,
                batch_key=BATCH_KEY if BATCH_KEY in adata_lineage.obs else None,
                flavor='seurat_v3',
                subset=False
            )
            lineage_hvg = set(adata_lineage.var_names[adata_lineage.var['highly_variable']])
            all_hvg_genes.update(lineage_hvg)
            print(f"     Added {len(lineage_hvg):,} HVG from {lineage}")
        except Exception as e:
            print(f"     ⚠️ Failed: {str(e)[:80]}")
            continue
    
    # Mark HVG in original adata
    adata.var['highly_variable'] = adata.var_names.isin(all_hvg_genes)
    hvg_method = "stratified-by-lineage"
    print(f"\n✓ Stratified HVG selection: {len(all_hvg_genes):,} genes")
    
    # Clean up temporary column
    if 'temp_major_lineage' in adata.obs.columns:
        adata.obs.drop('temp_major_lineage', axis=1, inplace=True)

else:
    print(f"\n3.1 Selecting top {N_HVG} HVGs (standard method)...")
    
    # Robust HVG selection
    try:
        sc.pp.highly_variable_genes(
            adata,
            layer='counts',
            n_top_genes=N_HVG,
            batch_key=BATCH_KEY,
            flavor='seurat_v3',
            subset=False
        )
        hvg_method = "batch-aware"
        print(f"✓ Batch-aware HVG selection successful")
    except Exception as e:
        print(f"⚠️  Batch-aware HVG failed: {str(e)[:100]}")
        print("   Falling back to non-batch-aware method...")
        sc.pp.highly_variable_genes(
            adata,
            layer='counts',
            n_top_genes=N_HVG,
            flavor='seurat_v3',
            subset=False
        )
        hvg_method = "non-batch-aware"
        print(f"✓ Non-batch-aware HVG selection successful")

# ===== P1 Strategy: Force-Include Critical Marker Genes =====
print(f"\n3.2 Force-including critical marker genes...")
print(f"   Strategy: Prevent AT2/rare cell collapse by including key markers")

# Find markers that exist in data but weren't selected as HVG
available_markers = [m for m in FORCE_INCLUDE_MARKERS if m in adata.var_names]
missing_from_hvg = [m for m in available_markers if not adata.var.loc[m, 'highly_variable']]

if len(missing_from_hvg) > 0:
    print(f"   Found {len(missing_from_hvg)} critical markers missing from HVG:")
    
    # Group by lineage for display
    marker_groups = {
        'Alveolar (AT1/AT2)': ['SFTPC', 'SFTPB', 'ABCA3', 'AGER', 'PDPN'],
        'Basal': ['TP63', 'KRT5', 'KRT14'],
        'Ciliated': ['FOXJ1', 'RSPH1', 'PIFO'],
        'Secretory': ['MUC5AC', 'MUC5B', 'SCGB1A1'],
        'Rare': ['CFTR', 'FOXI1', 'DCLK1']
    }
    
    for group_name, group_markers in marker_groups.items():
        group_missing = [m for m in group_markers if m in missing_from_hvg]
        if group_missing:
            print(f"     {group_name}: {', '.join(group_missing)}")
    
    # Force-include these markers
    adata.var.loc[missing_from_hvg, 'highly_variable'] = True
    print(f"\n✓ Force-included {len(missing_from_hvg)} critical markers")
    hvg_method += " + force-markers"
else:
    print(f"✓ All {len(available_markers)} critical markers already in HVG")

# Note missing markers
missing_markers = [m for m in FORCE_INCLUDE_MARKERS if m not in adata.var_names]
if len(missing_markers) > 0:
    print(f"\n⚠️  {len(missing_markers)} markers not found in data:")
    print(f"   {', '.join(missing_markers[:10])}...")

adata.uns['hvg_method'] = hvg_method
n_hvg = adata.var['highly_variable'].sum()
print(f"\n✓ Final gene set: {n_hvg:,} genes")
print(f"  Method: {hvg_method}")
print(f"  HVG: {n_hvg - len(missing_from_hvg):,}")
print(f"  Force-included: {len(missing_from_hvg):,}")

# Preserve full genes to .raw
print("\n3.3 Preserving full gene data in .raw (shared memory)...")
adata.raw = sc.AnnData(
    X=adata.layers['counts'],
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ adata.raw saved: {adata.raw.n_vars:,} genes (0 GB extra)")

# Subset to HVG + force-included markers
print("\n3.4 Subsetting to selected genes for model training...")
adata = adata[:, adata.var['highly_variable']].copy()
print(f"✓ Subset to {adata.n_vars:,} genes (HVG + critical markers)")
print(f"✓ Full gene access preserved in adata.raw: {adata.raw.n_vars:,} genes")

# Memory check
mem_gb = psutil.Process(os.getpid()).memory_info().rss / 1e9
print(f"\n💾 Memory after HVG subset: {mem_gb:.1f} GB")

# Save checkpoint
if SAVE_CHECKPOINTS:
    checkpoint_path = f"{CHECKPOINT_DIR}/02_hvg_selected.h5ad"
    adata.write_h5ad(checkpoint_path, compression='gzip')
    print(f"💾 Checkpoint saved: {checkpoint_path}")


# ===== STEP 4: scVI Model Training =====

print("\n" + "="*80)
print("STEP 4: scVI Model Training")
print("="*80)

print("\n4.1 Setting up scVI model...")

scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key=BATCH_KEY
)

print(f"✓ AnnData registered with scVI")

# Create scVI model
print("\n4.2 Creating scVI model...")
scvi_model = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT_SCVI,
    n_layers=N_LAYERS,
    n_hidden=N_HIDDEN,
    dropout_rate=DROPOUT_RATE,
    gene_likelihood=GENE_LIKELIHOOD,
    dispersion=DISPERSION
)

# Robust parameter counting
try:
    n_params = scvi_model.module.n_params
except AttributeError:
    n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)

print(f"✓ scVI model created")
print(f"  Parameters: {n_params:,}")
print(f"  Latent: {N_LATENT_SCVI}, Layers: {N_LAYERS}, Hidden: {N_HIDDEN}")

# Train scVI
print("\n4.3 Training scVI model...")
print(f"  Max epochs: {SCVI_MAX_EPOCHS}")
print(f"  Batch size: {BATCH_SIZE}")
print(f"  Device: {DEVICE}")

train_start = datetime.now()

scvi_model.train(
    max_epochs=SCVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=EARLY_STOPPING,
    train_size=0.9,
    plan_kwargs={'lr': 1e-3}
)

train_end = datetime.now()
train_duration = (train_end - train_start).total_seconds() / 60

print(f"\n✓ scVI training completed in {train_duration:.1f} minutes")

# Save scVI model
scvi_model_path = f"{OUTPUT_DIR}/scvi_model"
scvi_model.save(scvi_model_path, overwrite=True)
print(f"✓ scVI model saved to: {scvi_model_path}")

# Extract latent
print("\n4.4 Extracting scVI latent representation...")
adata.obsm['X_scvi'] = scvi_model.get_latent_representation()
print(f"✓ Shape: {adata.obsm['X_scvi'].shape}")

# Compute neighbors and UMAP
print("\n4.5 Computing scVI-based neighbors and UMAP...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
sc.tl.umap(adata)

# Leiden clustering
print("\n4.6 Computing Leiden clustering...")
for res in [0.5, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_scvi_r{res}')
    n_clusters = adata.obs[f'leiden_scvi_r{res}'].nunique()
    print(f"  Resolution {res}: {n_clusters} clusters")

print("✓ scVI analysis complete")


# ===== STEP 5: scVI Visualization =====

print("\n" + "="*80)
print("STEP 5: scVI Results Visualization")
print("="*80)

print("\n5.1 Creating scVI UMAP overview...")

fig, axes = plt.subplots(2, 3, figsize=(30, 20))

sc.pl.umap(adata, color=CELLTYPE_KEY, title='CellTypist Predictions',
           legend_loc='right margin', ax=axes[0, 0], show=False)

sc.pl.umap(adata, color='leiden_scvi_r1.0', title='Leiden (scVI, r=1.0)',
           legend_loc='right margin', ax=axes[0, 1], show=False)

sc.pl.umap(adata, color=BATCH_KEY, title=f'Batch Distribution ({n_batches} batches)',
           legend_loc='right margin', ax=axes[0, 2], show=False)

sc.pl.umap(adata, color='celltypist_conf', title='CellTypist Confidence',
           cmap='viridis', vmin=0, vmax=1, ax=axes[1, 0], show=False)

# Marker genes
key_markers = ['TP63', 'KRT5', 'FOXJ1', 'MUC5AC', 'SCGB1A1', 'SFTPC']
available_markers = [m for m in key_markers if m in adata.raw.var_names]

if len(available_markers) >= 2:
    sc.pl.umap(adata, color=available_markers[0], use_raw=True,
               title=f'{available_markers[0]} Expression',
               cmap='viridis', ax=axes[1, 1], show=False)
    sc.pl.umap(adata, color=available_markers[1], use_raw=True,
               title=f'{available_markers[1]} Expression',
               cmap='viridis', ax=axes[1, 2], show=False)
else:
    axes[1, 1].axis('off')
    axes[1, 2].axis('off')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/scvi_umap_overview.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: scvi_umap_overview.png")

# Training history
print("\n5.2 Plotting scVI training history...")
fig, ax = plt.subplots(figsize=(10, 6))
train_elbo = scvi_model.history['elbo_train'][1:]
val_elbo = scvi_model.history['elbo_validation'][1:]
epochs = range(1, len(train_elbo) + 1)

ax.plot(epochs, train_elbo, label='Training ELBO', linewidth=2)
ax.plot(epochs, val_elbo, label='Validation ELBO', linewidth=2)
ax.set_xlabel('Epoch')
ax.set_ylabel('ELBO')
ax.set_title('scVI Training History')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/scvi_training_history.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: scvi_training_history.png")

# Save checkpoint
if SAVE_CHECKPOINTS:
    checkpoint_path = f"{CHECKPOINT_DIR}/03_scvi_complete.h5ad"
    adata.write_h5ad(checkpoint_path, compression='gzip')
    print(f"\n💾 Checkpoint saved: {checkpoint_path}")

# Clean up
del scvi_model
gc.collect()
print("\n💾 Cleaned up scVI model from memory")


# ===== STEP 6: Prepare Hierarchical Labels for scANVI =====

print("\n" + "="*80)
print("STEP 6: Preparing Hierarchical Labels for scANVI")
print("="*80)

if USE_HIERARCHICAL_LABELS:
    print("\n⭐ Applying Hierarchical Label Consolidation")
    print("   Strategy: Fine types → 6 major lineages")
    
    # Store original predictions as fine_type (strip whitespace)
    adata.obs['fine_type'] = adata.obs[CELLTYPE_KEY].astype(str).str.strip()
    
    # Apply major lineage mapping
    adata.obs['major_lineage'] = adata.obs['fine_type'].map(MAJOR_LINEAGE_MAP)
    
    # Check for unmapped labels
    unmapped = adata.obs[adata.obs['major_lineage'].isna()]['fine_type'].unique()
    if len(unmapped) > 0:
        print(f"\n⚠️  WARNING: {len(unmapped)} labels not in mapping:")
        for label in unmapped:
            print(f"  - '{label}' (length: {len(label)})")
        
        print(f"\n💡 Suggestion: Add these to MAJOR_LINEAGE_MAP:")
        for label in unmapped:
            # Try to suggest appropriate lineage based on name
            suggested_lineage = "Unknown_Lineage"
            label_lower = label.lower()
            
            if 'basal' in label_lower or 'dividing' in label_lower:
                suggested_lineage = "Basal_Lineage"
            elif 'ciliat' in label_lower or 'deutero' in label_lower:
                suggested_lineage = "Ciliated_Lineage"
            elif 'secret' in label_lower or 'goblet' in label_lower or 'club' in label_lower or 'mucous' in label_lower or 'serous' in label_lower:
                suggested_lineage = "Secretory_Lineage"
            elif 'at1' in label_lower or 'at2' in label_lower or 'alveolar' in label_lower:
                suggested_lineage = "Alveolar"
            elif 'duct' in label_lower:
                suggested_lineage = "Duct"
            elif 'ionocyte' in label_lower or 'brush' in label_lower:
                suggested_lineage = "Rare_Specialized"
            
            print(f"  '{label}': '{suggested_lineage}',  # ⭐ Suggested")
        
        raise ValueError("Please update MAJOR_LINEAGE_MAP dictionary with above suggestions")
    
    # Convert to categorical
    adata.obs['major_lineage'] = pd.Categorical(adata.obs['major_lineage'])
    adata.obs['fine_type'] = pd.Categorical(adata.obs['fine_type'])
    
    # Statistics
    print("\n📊 Hierarchical Label Statistics:")
    print("\n=== Level 1: Major Lineages (for scANVI) ===")
    major_counts = adata.obs['major_lineage'].value_counts()
    for lineage, count in major_counts.items():
        pct = count / adata.n_obs * 100
        print(f"  {lineage:25s}: {count:7,} ({pct:5.2f}%)")
    
    print("\n=== Level 2: Fine Types (for detailed analysis) ===")
    fine_counts = adata.obs['fine_type'].value_counts()
    for ftype, count in fine_counts.items():
        pct = count / adata.n_obs * 100
        lineage = MAJOR_LINEAGE_MAP.get(ftype, 'Unknown')
        print(f"  {ftype:25s} → {lineage:25s}: {count:7,} ({pct:5.2f}%)")
    
    # Use major lineages for scANVI
    source_key = 'major_lineage'
    print(f"\n✓ Using major lineages for scANVI training")
    
else:
    # No hierarchical labels, use CellTypist predictions directly
    adata.obs['fine_type'] = adata.obs[CELLTYPE_KEY].astype(str)
    source_key = 'fine_type'
    print(f"\n✓ Using fine types directly (no hierarchical consolidation)")

# Create scanvi_labels
print(f"\n6.2 Preparing scANVI labels...")
adata.obs['scanvi_labels'] = adata.obs[source_key].astype(str)

# Mark low confidence as Unknown
if MARK_LOW_CONFIDENCE_AS_UNKNOWN:
    low_conf_mask = adata.obs['celltypist_conf'] < LOW_CONFIDENCE_THRESHOLD
    n_low_conf = low_conf_mask.sum()
    
    if n_low_conf > 0:
        adata.obs.loc[low_conf_mask, 'scanvi_labels'] = UNLABELED_CATEGORY
        pct = n_low_conf / adata.n_obs * 100
        print(f"✓ Marked {n_low_conf:,} low-confidence cells ({pct:.2f}%) as '{UNLABELED_CATEGORY}'")
        
        if pct > MAX_LOW_CONFIDENCE_PCT:
            print(f"\n⚠️  WARNING: {pct:.1f}% marked as Unknown")
            print(f"   scANVI will be mostly unsupervised")
    else:
        print(f"✓ No cells below confidence threshold {LOW_CONFIDENCE_THRESHOLD}")
else:
    # Random 5% as Unknown
    n_unknown = int(adata.n_obs * 0.05)
    unknown_idx = np.random.choice(adata.n_obs, n_unknown, replace=False)
    adata.obs['scanvi_labels'] = adata.obs['scanvi_labels'].astype('object')
    adata.obs.loc[adata.obs.index[unknown_idx], 'scanvi_labels'] = UNLABELED_CATEGORY
    print(f"✓ Randomly marked {n_unknown:,} cells (5%) as '{UNLABELED_CATEGORY}'")

# Convert to categorical
adata.obs['scanvi_labels'] = pd.Categorical(adata.obs['scanvi_labels'])

# Summary
n_labeled = (adata.obs['scanvi_labels'] != UNLABELED_CATEGORY).sum()
n_unlabeled = (adata.obs['scanvi_labels'] == UNLABELED_CATEGORY).sum()
n_categories = adata.obs['scanvi_labels'].nunique()

print(f"\n📊 scANVI Label Summary:")
print(f"  Total: {adata.n_obs:,}")
print(f"  Labeled: {n_labeled:,} ({n_labeled/adata.n_obs*100:.2f}%)")
print(f"  Unknown: {n_unlabeled:,} ({n_unlabeled/adata.n_obs*100:.2f}%)")
print(f"  Categories: {n_categories}")

print(f"\n📋 Label distribution:")
for label, count in adata.obs['scanvi_labels'].value_counts().items():
    pct = count / adata.n_obs * 100
    print(f"  {str(label):25s}: {count:7,} ({pct:5.2f}%)")


# ===== STEP 6.5: Visualize Hierarchical Structure (Optional) =====

if USE_HIERARCHICAL_LABELS and VISUALIZE_HIERARCHY:
    print("\n" + "="*80)
    print("STEP 6.5: Visualizing Hierarchical Label Structure")
    print("="*80)
    
    hierarchy_dir = f"{OUTPUT_DIR}/figures/hierarchical_structure"
    os.makedirs(hierarchy_dir, exist_ok=True)
    
    # Figure 1: UMAP comparison
    print("\n6.5.1 Creating hierarchical UMAP comparison...")
    fig, axes = plt.subplots(2, 2, figsize=(24, 20))
    
    sc.pl.umap(adata, color='major_lineage', 
               title='Level 1: Major Lineages (for scANVI)',
               palette=MAJOR_LINEAGE_COLORS,
               legend_loc='right margin', ax=axes[0, 0], show=False)
    
    sc.pl.umap(adata, color='fine_type', 
               title='Level 2: Fine Types (CellTypist)',
               palette=FINE_TYPE_COLORS,
               legend_loc='right margin', ax=axes[0, 1], show=False)
    
    sc.pl.umap(adata, color='celltypist_conf',
               title='CellTypist Confidence',
               cmap='viridis', vmin=0, vmax=1, ax=axes[1, 0], show=False)
    
    sc.pl.umap(adata, color=BATCH_KEY,
               title='Dataset Distribution', ax=axes[1, 1], show=False)
    
    plt.tight_layout()
    plt.savefig(f'{hierarchy_dir}/hierarchical_umap_comparison.png', dpi=300)
    plt.close()
    print("  ✓ Saved: hierarchical_umap_comparison.png")
    
    # Figure 2: Hierarchical structure
    print("\n6.5.2 Creating hierarchical structure plot...")
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))
    
    # Major lineages bar plot
    major_counts = adata.obs['major_lineage'].value_counts()
    colors_l1 = [MAJOR_LINEAGE_COLORS[x] for x in major_counts.index]
    axes[0].barh(range(len(major_counts)), major_counts.values, color=colors_l1)
    axes[0].set_yticks(range(len(major_counts)))
    axes[0].set_yticklabels(major_counts.index)
    axes[0].set_xlabel('Number of Cells')
    axes[0].set_title('Level 1: Major Lineages (Use for scANVI)', 
                      fontsize=14, fontweight='bold')
    axes[0].invert_yaxis()
    
    for i, (lineage, count) in enumerate(major_counts.items()):
        pct = count / adata.n_obs * 100
        axes[0].text(count, i, f'  {count:,} ({pct:.1f}%)', 
                    va='center', fontsize=10)
    
    # Fine types grouped bar plot
    fine_by_major = adata.obs.groupby(['major_lineage', 'fine_type']).size().reset_index(name='count')
    
    y_pos = 0
    yticks = []
    yticklabels = []
    
    for major in major_counts.index:
        subset = fine_by_major[fine_by_major['major_lineage'] == major]
        n_subtypes = len(subset)
        
        for i, (_, row) in enumerate(subset.iterrows()):
            color = FINE_TYPE_COLORS.get(row['fine_type'], '#95A5A6')
            axes[1].barh(y_pos, row['count'], color=color, alpha=0.8)
            yticks.append(y_pos)
            yticklabels.append(f"  {row['fine_type']}")
            
            pct = row['count'] / adata.n_obs * 100
            axes[1].text(row['count'], y_pos, f'  {row["count"]:,} ({pct:.1f}%)', 
                        va='center', fontsize=9)
            
            y_pos += 1
        
        if n_subtypes > 0:
            mid_pos = y_pos - n_subtypes/2 - 0.5
            axes[1].text(-adata.n_obs*0.05, mid_pos, major, 
                        va='center', ha='right', fontsize=11, fontweight='bold',
                        color=MAJOR_LINEAGE_COLORS[major])
        
        y_pos += 0.5
    
    axes[1].set_yticks(yticks)
    axes[1].set_yticklabels(yticklabels, fontsize=9)
    axes[1].set_xlabel('Number of Cells')
    axes[1].set_title('Level 2: Fine Types (For Detailed Analysis)', 
                      fontsize=14, fontweight='bold')
    axes[1].invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(f'{hierarchy_dir}/hierarchical_structure.png', dpi=300)
    plt.close()
    print("  ✓ Saved: hierarchical_structure.png")
    
    # Figure 3: Transition matrix
    print("\n6.5.3 Creating lineage-type matrix...")
    transition_df = pd.crosstab(
        adata.obs['major_lineage'], 
        adata.obs['fine_type'],
        normalize='index'
    ) * 100
    
    fig, ax = plt.subplots(figsize=(14, 8))
    sns.heatmap(transition_df, annot=True, fmt='.1f', 
                cmap='YlOrRd', cbar_kws={'label': 'Percentage (%)'},
                ax=ax)
    ax.set_title('Fine Type Distribution within Major Lineages', 
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Fine Type')
    ax.set_ylabel('Major Lineage')
    plt.tight_layout()
    plt.savefig(f'{hierarchy_dir}/lineage_type_matrix.png', dpi=300)
    plt.close()
    print("  ✓ Saved: lineage_type_matrix.png")
    
    # Save mapping table
    mapping_df = pd.DataFrame([
        {'Fine_Type': k, 'Major_Lineage': v, 'Color': FINE_TYPE_COLORS.get(k, '#95A5A6')}
        for k, v in MAJOR_LINEAGE_MAP.items()
    ])
    mapping_df = mapping_df.sort_values(['Major_Lineage', 'Fine_Type'])
    mapping_df.to_csv(f'{OUTPUT_DIR}/hierarchical_label_mapping.csv', index=False)
    print("\n  ✓ Saved: hierarchical_label_mapping.csv")
    
    print("\n✓ Hierarchical structure visualization complete")


# ===== STEP 7: scANVI Model Training (Version-Safe) =====

print("\n" + "="*80)
print("STEP 7: scANVI Model Training (Version-Safe)")
print("="*80)

print("\n7.1 Loading scVI model for scANVI initialization...")

try:
    scvi_model = scvi.model.SCVI.load(
        f"{OUTPUT_DIR}/scvi_model", 
        adata=adata,
        map_location=DEVICE
    )
    print(f"✓ scVI model loaded (device: {DEVICE})")
except Exception as e:
    raise RuntimeError(f"Failed to load scVI model. Error: {e}")

print("\n7.2 Creating scANVI model from scVI...")

try:
    print("   Attempting standard initialization (scvi-tools v1.0+)...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        unlabeled_category=UNLABELED_CATEGORY,
        labels_key='scanvi_labels'
    )
    print("✓ scANVI initialized successfully")
    
except (ValueError, RuntimeError) as e:
    print(f"⚠️  Standard init failed: {str(e)[:100]}")
    print("   Attempting explicit setup fallback...")
    
    scvi.model.SCANVI.setup_anndata(
        adata,
        layer='counts',
        batch_key=BATCH_KEY,
        labels_key='scanvi_labels',
        unlabeled_category=UNLABELED_CATEGORY
    )
    
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        unlabeled_category=UNLABELED_CATEGORY,
        labels_key='scanvi_labels'
    )
    print("✓ scANVI initialized with fallback method")

# Robust parameter counting
try:
    n_params = scanvi_model.module.n_params
except AttributeError:
    n_params = sum(p.numel() for p in scanvi_model.module.parameters() if p.requires_grad)

print(f"\n✓ scANVI model created")
print(f"  Parameters: {n_params:,}")
print(f"  Unlabeled category: {UNLABELED_CATEGORY}")
print(f"  Labeled cells: {n_labeled:,}")
print(f"  Unlabeled cells: {n_unlabeled:,}")

# Clean up scVI model
del scvi_model
gc.collect()

# Train scANVI
print("\n7.3 Training scANVI model...")
print(f"  Max epochs: {SCANVI_MAX_EPOCHS}")
print(f"  Batch size: {BATCH_SIZE}")
print(f"  Device: {DEVICE}")

train_start = datetime.now()

scanvi_model.train(
    max_epochs=SCANVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=EARLY_STOPPING,
    train_size=0.9,
    plan_kwargs={'lr': 5e-4}
)

train_end = datetime.now()
scanvi_train_duration = (train_end - train_start).total_seconds() / 60

print(f"\n✓ scANVI training completed in {scanvi_train_duration:.1f} minutes")

# Save scANVI model
scanvi_model_path = f"{OUTPUT_DIR}/scanvi_model"
scanvi_model.save(scanvi_model_path, overwrite=True)
print(f"✓ scANVI model saved to: {scanvi_model_path}")


# ===== STEP 8: scANVI Predictions =====

print("\n" + "="*80)
print("STEP 8: scANVI Predictions")
print("="*80)

print("\n8.1 Predicting cell types...")
adata.obs['scanvi_pred'] = scanvi_model.predict()
print("✓ Predictions complete")

print("\n8.2 Computing prediction probabilities...")
scanvi_predictions = scanvi_model.predict(soft=True)
adata.obsm['scanvi_prob'] = scanvi_predictions
adata.obs['scanvi_conf'] = scanvi_predictions.max(axis=1)
print("✓ Probabilities computed")

print("\n8.3 Extracting scANVI latent representation...")
adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
print(f"✓ Shape: {adata.obsm['X_scanvi'].shape}")

print("\n8.4 Computing scANVI neighbors and UMAP...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15, key_added='scanvi')
sc.tl.umap(adata, neighbors_key='scanvi')
print("✓ UMAP computed")

print("\n8.5 Computing Leiden clustering...")
for res in [0.5, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_scanvi_r{res}', 
                 neighbors_key='scanvi')
    n_clusters = adata.obs[f'leiden_scanvi_r{res}'].nunique()
    print(f"  Resolution {res}: {n_clusters} clusters")

# Summary
print("\n📊 scANVI Prediction Summary:")
print(f"\nPredicted categories: {adata.obs['scanvi_pred'].nunique()}")

if USE_HIERARCHICAL_LABELS:
    print(f"\n⭐ scANVI Predictions (Major Lineages):")
else:
    print(f"\nscANVI Predictions:")

for celltype, count in adata.obs['scanvi_pred'].value_counts().items():
    pct = count / adata.n_obs * 100
    avg_conf = adata.obs[adata.obs['scanvi_pred'] == celltype]['scanvi_conf'].mean()
    print(f"  {str(celltype):35s}: {count:7,} ({pct:5.2f}%) - Conf: {avg_conf:.3f}")

print(f"\nConfidence statistics:")
print(f"  Mean: {adata.obs['scanvi_conf'].mean():.3f}")
print(f"  Median: {adata.obs['scanvi_conf'].median():.3f}")

# Compare with input labels
changed_cells = (adata.obs[source_key].astype(str) != 
                 adata.obs['scanvi_pred'].astype(str)).sum()
print(f"\n8.6 Comparison with input labels:")
print(f"✓ scANVI changed {changed_cells:,} cells ({changed_cells/adata.n_obs*100:.1f}%)")

# Clean up
del scanvi_model
gc.collect()
print("\n💾 Cleaned up scANVI model from memory")


# ===== STEP 9: Comprehensive Visualization =====

print("\n" + "="*80)
print("STEP 9: Comprehensive Visualization")
print("="*80)

print("\n9.1 Creating scANVI UMAP overview...")

fig, axes = plt.subplots(2, 3, figsize=(30, 20))

if USE_HIERARCHICAL_LABELS:
    # Show both major lineages (input) and scANVI predictions
    sc.pl.umap(adata, color='major_lineage', 
               title='Input: Major Lineages',
               palette=MAJOR_LINEAGE_COLORS,
               legend_loc='right margin', ax=axes[0, 0], show=False)
    
    sc.pl.umap(adata, color='scanvi_pred', 
               title='scANVI: Refined Major Lineages',
               palette=MAJOR_LINEAGE_COLORS,
               legend_loc='right margin', ax=axes[0, 1], show=False)
else:
    sc.pl.umap(adata, color=CELLTYPE_KEY, 
               title='Input: CellTypist Predictions',
               legend_loc='right margin', ax=axes[0, 0], show=False)
    
    sc.pl.umap(adata, color='scanvi_pred', 
               title='scANVI: Refined Predictions',
               legend_loc='right margin', ax=axes[0, 1], show=False)

sc.pl.umap(adata, color='leiden_scanvi_r1.0', 
           title='Leiden (scANVI, r=1.0)',
           legend_loc='right margin', ax=axes[0, 2], show=False)

sc.pl.umap(adata, color='scanvi_conf', 
           title='scANVI Confidence',
           cmap='viridis', vmin=0, vmax=1, ax=axes[1, 0], show=False)

sc.pl.umap(adata, color=BATCH_KEY, 
           title='Batch Distribution',
           legend_loc='right margin', ax=axes[1, 1], show=False)

if len(available_markers) >= 1:
    sc.pl.umap(adata, color=available_markers[0], use_raw=True,
               title=f'{available_markers[0]} Expression',
               cmap='viridis', ax=axes[1, 2], show=False)
else:
    axes[1, 2].axis('off')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/scanvi_umap_overview.png', dpi=300)
plt.close()
print("✓ Saved: scanvi_umap_overview.png")

# Comparison plots
if USE_HIERARCHICAL_LABELS:
    print("\n9.2 Creating hierarchical comparison plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(24, 20))
    
    # Major lineages (input)
    sc.pl.umap(adata, color='major_lineage', 
               title='Input: Major Lineages (Level 1)',
               palette=MAJOR_LINEAGE_COLORS,
               legend_loc='right margin', ax=axes[0, 0], show=False)
    
    # Fine types (original CellTypist)
    sc.pl.umap(adata, color='fine_type', 
               title='Original: Fine Types (Level 2)',
               palette=FINE_TYPE_COLORS,
               legend_loc='right margin', ax=axes[0, 1], show=False)
    
    # scANVI predictions (major lineages)
    sc.pl.umap(adata, color='scanvi_pred', 
               title='scANVI: Refined Major Lineages',
               palette=MAJOR_LINEAGE_COLORS,
               legend_loc='right margin', ax=axes[1, 0], show=False)
    
    # Confidence comparison
    conf_df = pd.DataFrame({
        'CellTypist': adata.obs['celltypist_conf'],
        'scANVI': adata.obs['scanvi_conf']
    })
    conf_df.plot.hist(bins=50, alpha=0.7, ax=axes[1, 1])
    axes[1, 1].set_xlabel('Confidence Score')
    axes[1, 1].set_ylabel('Number of Cells')
    axes[1, 1].set_title('Confidence Comparison')
    axes[1, 1].axvline(0.5, color='red', linestyle='--', alpha=0.5)
    axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/figures/hierarchical_comparison.png', dpi=300)
    plt.close()
    print("✓ Saved: hierarchical_comparison.png")

else:
    print("\n9.2 Creating standard comparison plots...")
    
    fig, axes = plt.subplots(1, 3, figsize=(30, 10))
    
    sc.pl.umap(adata, color=CELLTYPE_KEY, 
               title='CellTypist Annotations',
               legend_loc='right margin', ax=axes[0], show=False)
    
    sc.pl.umap(adata, color='scanvi_pred', 
               title='scANVI Refined Annotations',
               legend_loc='right margin', ax=axes[1], show=False)
    
    conf_df = pd.DataFrame({
        'CellTypist': adata.obs['celltypist_conf'],
        'scANVI': adata.obs['scanvi_conf']
    })
    conf_df.plot.hist(bins=50, alpha=0.7, ax=axes[2])
    axes[2].set_xlabel('Confidence Score')
    axes[2].set_ylabel('Number of Cells')
    axes[2].set_title('Confidence Comparison')
    axes[2].axvline(0.5, color='red', linestyle='--', alpha=0.5)
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/figures/comparison.png', dpi=300)
    plt.close()
    print("✓ Saved: comparison.png")

# Confusion matrix
print("\n9.3 Creating confusion matrix...")

ct_labels = adata.obs[source_key].astype(str).unique()
scanvi_labels = adata.obs['scanvi_pred'].astype(str).unique()
all_labels = sorted(set(list(ct_labels) + list(scanvi_labels)))

cm = confusion_matrix(
    adata.obs[source_key].astype(str),
    adata.obs['scanvi_pred'].astype(str),
    labels=all_labels
)

fig, ax = plt.subplots(figsize=(max(12, len(all_labels)*0.5), 
                                 max(10, len(all_labels)*0.4)))
sns.heatmap(cm, annot=False, fmt='d', cmap='YlOrRd',
            xticklabels=all_labels, yticklabels=all_labels,
            cbar_kws={'label': 'Number of Cells'}, ax=ax)
ax.set_xlabel('scANVI Predictions')
ax.set_ylabel('Input Labels')
ax.set_title('Confusion Matrix')
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/confusion_matrix.png', dpi=300)
plt.close()
print("✓ Saved: confusion_matrix.png")


# ===== STEP 10: Save Final Results =====

print("\n" + "="*80)
print("STEP 10: Saving Final Results")
print("="*80)

# Save final data
output_h5ad = f"{OUTPUT_DIR}/epithelial_scvi_scanvi_hierarchical_final.h5ad"
adata.write_h5ad(output_h5ad, compression='gzip')
print(f"✓ Saved: {output_h5ad}")
print(f"  Size: {os.path.getsize(output_h5ad) / 1024**3:.2f} GB")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes (HVG): {adata.n_vars:,}")
print(f"  Full genes (.raw): {adata.raw.n_vars:,}")

# Save statistics
stats_file = f"{OUTPUT_DIR}/annotation_statistics.csv"

if USE_HIERARCHICAL_LABELS:
    stats_df = pd.DataFrame({
        'Fine_Type': adata.obs['fine_type'].value_counts(),
        'Major_Lineage': adata.obs['major_lineage'].value_counts(),
        'scANVI_Prediction': adata.obs['scanvi_pred'].value_counts()
    })
else:
    stats_df = pd.DataFrame({
        'CellTypist': adata.obs[CELLTYPE_KEY].value_counts(),
        'scANVI': adata.obs['scanvi_pred'].value_counts()
    })

for col in stats_df.columns:
    stats_df[f'{col}_pct'] = (stats_df[col] / adata.n_obs * 100).round(2)

stats_df = stats_df.sort_values(stats_df.columns[0], ascending=False)
stats_df.to_csv(stats_file)
print(f"✓ Saved: {stats_file}")

# Save report
report_file = f"{OUTPUT_DIR}/analysis_report.txt"
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("Epithelial scVI-scANVI Analysis Report (v2.3)\n")
    if USE_HIERARCHICAL_LABELS:
        f.write("With Hierarchical Label Consolidation\n")
    f.write("="*80 + "\n\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"scvi-tools version: {scvi.__version__}\n")
    f.write(f"PyTorch version: {torch.__version__}\n")
    f.write(f"Device: {DEVICE}\n\n")
    
    f.write("Configuration:\n")
    f.write(f"  Input: {INPUT_H5AD}\n")
    f.write(f"  Batch key: {BATCH_KEY}\n")
    f.write(f"  Min cells per batch: {MIN_CELLS_PER_BATCH}\n")
    f.write(f"  HVG: {N_HVG} ({hvg_method})\n")
    f.write(f"  scVI latent: {N_LATENT_SCVI}, Batch size: {BATCH_SIZE}\n")
    f.write(f"  Dispersion: {DISPERSION}\n\n")
    
    if USE_HIERARCHICAL_LABELS:
        f.write("Hierarchical Labels:\n")
        f.write(f"  Fine types: {len(MAJOR_LINEAGE_MAP)}\n")
        f.write(f"  Major lineages: {len(set(MAJOR_LINEAGE_MAP.values()))}\n")
        f.write(f"  Strategy: Coarse-to-fine (major lineages for scANVI)\n\n")
    
    f.write("Training Time:\n")
    f.write(f"  scVI: {train_duration:.1f} min\n")
    f.write(f"  scANVI: {scanvi_train_duration:.1f} min\n")
    f.write(f"  Total: {train_duration + scanvi_train_duration:.1f} min\n\n")
    
    f.write("Dataset:\n")
    f.write(f"  Cells: {adata.n_obs:,}\n")
    f.write(f"  Genes (HVG): {adata.n_vars:,}\n")
    f.write(f"  Full genes: {adata.raw.n_vars:,}\n")
    f.write(f"  Batches: {n_batches}\n\n")
    
    if USE_HIERARCHICAL_LABELS:
        f.write("Major Lineages (Level 1 - for scANVI):\n")
        f.write(f"  Categories: {adata.obs['major_lineage'].nunique()}\n\n")
        for lineage, count in adata.obs['major_lineage'].value_counts().items():
            pct = count / adata.n_obs * 100
            f.write(f"  {str(lineage):30s}: {count:7,} ({pct:6.2f}%)\n")
        
        f.write("\nFine Types (Level 2 - for detailed analysis):\n")
        f.write(f"  Categories: {adata.obs['fine_type'].nunique()}\n\n")
        for ftype, count in adata.obs['fine_type'].value_counts().items():
            pct = count / adata.n_obs * 100
            lineage = MAJOR_LINEAGE_MAP.get(ftype, 'Unknown')
            f.write(f"  {str(ftype):30s} → {lineage:20s}: {count:7,} ({pct:6.2f}%)\n")
    else:
        f.write("CellTypist Annotations:\n")
        f.write(f"  Cell types: {adata.obs[CELLTYPE_KEY].nunique()}\n")
        f.write(f"  Avg confidence: {adata.obs['celltypist_conf'].mean():.3f}\n\n")
        for ct, count in adata.obs[CELLTYPE_KEY].value_counts().items():
            pct = count / adata.n_obs * 100
            avg_conf = adata.obs[adata.obs[CELLTYPE_KEY] == ct]['celltypist_conf'].mean()
            f.write(f"  {str(ct):40s}: {count:7,} ({pct:6.2f}%) - {avg_conf:.3f}\n")
    
    f.write("\nscANVI Predictions:\n")
    f.write(f"  Cell types: {adata.obs['scanvi_pred'].nunique()}\n")
    f.write(f"  Avg confidence: {adata.obs['scanvi_conf'].mean():.3f}\n")
    f.write(f"  Changed: {changed_cells:,} ({changed_cells/adata.n_obs*100:.2f}%)\n\n")
    for ct, count in adata.obs['scanvi_pred'].value_counts().items():
        pct = count / adata.n_obs * 100
        avg_conf = adata.obs[adata.obs['scanvi_pred'] == ct]['scanvi_conf'].mean()
        f.write(f"  {str(ct):40s}: {count:7,} ({pct:6.2f}%) - {avg_conf:.3f}\n")

print(f"✓ Saved: {report_file}")

# Final memory
mem_gb = psutil.Process(os.getpid()).memory_info().rss / 1e9
print(f"\n💾 Final memory: {mem_gb:.1f} GB")

print("\n" + "="*80)
print("✅ PIPELINE COMPLETED SUCCESSFULLY!")
print("="*80)
print(f"\nCompleted at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time: {train_duration + scanvi_train_duration:.1f} min")
print(f"\nKey results:")
print(f"  - {adata.n_obs:,} cells analyzed")
if USE_HIERARCHICAL_LABELS:
    print(f"  - {adata.obs['major_lineage'].nunique()} major lineages identified")
    print(f"  - {adata.obs['fine_type'].nunique()} fine types preserved")
else:
    print(f"  - {adata.obs['scanvi_pred'].nunique()} cell types identified")
print(f"  - {n_batches} batches integrated")
print(f"  - {adata.raw.n_vars:,} genes available for analysis")
print(f"\nAnnotations:")
if USE_HIERARCHICAL_LABELS:
    print(f"  - adata.obs['major_lineage']: Level 1 (6 major lineages)")
    print(f"  - adata.obs['fine_type']: Level 2 (13 fine types)")
    print(f"  - adata.obs['scanvi_pred']: scANVI refined (major lineages)")
else:
    print(f"  - adata.obs['{CELLTYPE_KEY}']: CellTypist")
    print(f"  - adata.obs['scanvi_pred']: scANVI refined")
print(f"  - adata.obs['scanvi_conf']: Confidence scores")
print(f"  - adata.obsm['X_scvi']: scVI latent")
print(f"  - adata.obsm['X_scanvi']: scANVI latent")
print(f"  - adata.raw: Full gene access")
print(f"\nAll results in: {OUTPUT_DIR}/")
if USE_HIERARCHICAL_LABELS and VISUALIZE_HIERARCHY:
    print(f"Hierarchical visualizations in: {OUTPUT_DIR}/figures/hierarchical_structure/")
if SAVE_CHECKPOINTS:
    print(f"Checkpoints in: {CHECKPOINT_DIR}/")
print("="*80)

if USE_HIERARCHICAL_LABELS:
    print("\n" + "="*80)
    print("⭐ NEXT STEPS: Sub-Analysis within Major Lineages")
    print("="*80)
    print("\nNow that scANVI has learned robust major lineage structure,")
    print("you can perform detailed sub-analysis within each lineage:")
    print("\nFor each major lineage:")
    print("  1. Subset cells: adata[adata.obs['scanvi_pred'] == lineage]")
    print("  2. Re-cluster: sc.tl.leiden(adata_sub, resolution=0.5)")
    print("  3. Validate fine types: Compare with adata.obs['fine_type']")
    print("  4. Differential expression: sc.tl.rank_genes_groups(use_raw=True)")
    print("  5. Marker validation: sc.pl.umap(color=marker_genes, use_raw=True)")
    print("\nThis coarse-to-fine strategy ensures:")
    print("  ✓ Robust batch correction at major lineage level")
    print("  ✓ Preserved fine-grained biological variation")
    print("  ✓ Clean latent space for downstream analysis")
    print("="*80)