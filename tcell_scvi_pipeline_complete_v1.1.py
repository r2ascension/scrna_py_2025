# %% [markdown]
# # T/NK Cell scVI Integration - Complete Pipeline
# 
# **Objective**: Batch integration using scVI for T/NK cell analysis
# 
# **Pipeline Overview**:
# 1. Import libraries and configuration
# 2. Data loading and raw counts preparation
# 3. Preprocessing and HVG selection
# 4. scVI model training
# 5. Clustering and UMAP
# 6. Quality assessment
# 7. Visualization
# 8. Results saving

# %% [markdown]
# ---
# ## Step 1: Import Libraries

# %%
# Import libraries
import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns
import time
from datetime import datetime
from scipy.stats import entropy
from sklearn.neighbors import NearestNeighbors

# Single-cell analysis
import scanpy as sc
import scvi

warnings.filterwarnings('ignore')

# Check versions
print(f"scanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"Python: {sys.version}")

# Check GPU
import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")

# Plot settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
plt.rcParams['figure.figsize'] = (8, 6)

print("\n✓ All libraries loaded successfully")

# %% [markdown]
# ---
# ## Step 2: Configuration

# %%
# ==================== Configuration ====================

# Input/Output
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1203/scvi_integration"

# Batch configuration
BATCH_KEY = "dataset"
CONFOUNDER_KEY = ["tissue"]  # Biological variables to preserve

# Gene filtering
MIN_CELLS_PER_GENE = 10
N_TOP_GENES = 4000

# scVI model parameters
SCVI_N_LATENT = 30       # Latent dimensions
SCVI_N_LAYERS = 2        # Neural network depth
SCVI_DROPOUT = 0.1       # Dropout rate
MAX_EPOCHS = 400         # Training epochs
EARLY_STOPPING = True    # Auto-stop if converged
LEARNING_RATE = 1e-3     # Learning rate

# UMAP parameters
UMAP_MIN_DIST = 0.15     # Smaller = tighter clusters
UMAP_SPREAD = 1.3        # Larger = more dispersed
UMAP_N_NEIGHBORS = 30

# Clustering
LEIDEN_RESOLUTIONS = [2.0, 2.5, 3.0, 3.5, 4.0]
DEFAULT_RESOLUTION = 3.0

# Marker genes
MARKER_GENES = {
    'Pan_T': ['CD3D', 'CD3E', 'CD3G'],
    'NK': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],
    'CD4': ['CD4', 'IL7R'],
    'CD8': ['CD8A', 'CD8B'],
    'Naive': ['CCR7', 'SELL', 'LEF1'],
    'Effector': ['GZMA', 'GZMB', 'PRF1'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4'],
    'Exhaustion': ['PDCD1', 'LAG3', 'TIGIT'],
    'Proliferation': ['MKI67', 'TOP2A']
}

# Visualization
DPI = 300
FIGURE_FORMAT = "pdf"

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)
model_dir = output_dir / "scvi_model"

print("✓ Configuration loaded")
print(f"  Input: {INPUT_H5AD}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Batch key: {BATCH_KEY}")

# %% [markdown]
# ---
# ## Step 3: Load Data and Prepare Raw Counts

# %%
print("\n" + "="*70)
print("Step 3: Loading Data")
print("="*70)

# Load data
adata = sc.read_h5ad(INPUT_H5AD)
print(f"\nLoaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Check batch info
n_batches = adata.obs[BATCH_KEY].nunique()
print(f"\nBatch information:")
print(f"  Number of batches: {n_batches}")
batch_counts = adata.obs[BATCH_KEY].value_counts()
for batch, count in list(batch_counts.items())[:10]:
    print(f"    {batch}: {count:,} cells")

# Check for confounders
has_confounder = all(c in adata.obs.columns for c in CONFOUNDER_KEY)
if has_confounder:
    print(f"\nConfounder keys found: {CONFOUNDER_KEY}")
else:
    print(f"\nWarning: Confounder keys {CONFOUNDER_KEY} not all found")
    print("Will proceed without confounders")

# %%
# Prepare raw counts for scVI
print("\nPreparing raw counts...")

if adata.raw is not None:
    print("✓ Using .raw for counts")
    adata.layers['counts'] = adata.raw.X.copy()
elif 'counts' in adata.layers:
    print("✓ Using layers['counts']")
else:
    print("⚠️  Using .X as counts (verify this is raw counts!)")
    adata.layers['counts'] = adata.X.copy()

# Store original for later use
adata_full = adata.copy()
print("✓ Raw counts prepared")

# %% [markdown]
# ---
# ## Step 4: Preprocessing

# %%
print("\n" + "="*70)
print("Step 4: Preprocessing")
print("="*70)

# Filter genes
print(f"\nFiltering genes (min_cells={MIN_CELLS_PER_GENE})...")
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
print(f"  Genes: {n_genes_before:,} → {adata.n_vars:,}")

# Normalize and select HVGs
print(f"\nSelecting HVGs...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP_GENES,
    flavor='seurat_v3',
    batch_key=BATCH_KEY,
    subset=True
)

print(f"  Selected: {adata.n_vars:,} HVGs")
print(f"\n✓ Preprocessing complete")
print(f"  Final: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# %% [markdown]
# ---
# ## Step 5: Setup and Train scVI

# %%
print("\n" + "="*70)
print("Step 5: Setting up scVI")
print("="*70)

# Setup AnnData for scVI
setup_kwargs = {
    'layer': 'counts',
    'batch_key': BATCH_KEY
}

if has_confounder:
    setup_kwargs['categorical_covariate_keys'] = CONFOUNDER_KEY
    print(f"  Including confounders: {CONFOUNDER_KEY}")

scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)
print("✓ AnnData setup complete")

# Create model
print("\nCreating scVI model...")
vae = scvi.model.SCVI(
    adata,
    n_latent=SCVI_N_LATENT,
    n_layers=SCVI_N_LAYERS,
    dropout_rate=SCVI_DROPOUT,
    gene_likelihood='nb'
)
print(f"✓ Model created (n_latent={SCVI_N_LATENT}, n_layers={SCVI_N_LAYERS})")

# %%
# Train model
print("\n" + "="*70)
print("Training scVI Model")
print("="*70)
print(f"Max epochs: {MAX_EPOCHS}")
print(f"Early stopping: {EARLY_STOPPING}")
print(f"GPU available: {gpu_available}")
print("\nTraining started (this may take 10-30 minutes)...\n")

start_time = time.time()

vae.train(
    max_epochs=MAX_EPOCHS,
    early_stopping=EARLY_STOPPING,
    early_stopping_patience=15,
    plan_kwargs={'lr': LEARNING_RATE}
)

elapsed = time.time() - start_time
train_elbo = vae.history['elbo_train']

print(f"\n{'='*70}")
print(f"✓ Training completed")
print(f"  Time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
print(f"  Epochs: {len(train_elbo)}")
print(f"  Final ELBO: {train_elbo[-1]:.2f}")
print("="*70)

# %%
# Plot training history
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(train_elbo, linewidth=2)
axes[0].set_xlabel('Epoch')
axes[0].set_ylabel('ELBO')
axes[0].set_title('Training History')
axes[0].grid(alpha=0.3)

axes[1].plot(train_elbo[-50:], linewidth=2)
axes[1].set_xlabel('Epoch (last 50)')
axes[1].set_ylabel('ELBO')
axes[1].set_title('Convergence Detail')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(fig_dir / f'training_history.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Training history saved")

# %% [markdown]
# ---
# ## Step 6: Extract Latent Representation

# %%
print("\n" + "="*70)
print("Step 6: Extracting Latent Representation")
print("="*70)

# Get latent representation
latent = vae.get_latent_representation()
adata.obsm['X_scvi'] = latent

print(f"✓ Latent representation extracted")
print(f"  Shape: {latent.shape}")

# Save model
vae.save(model_dir, overwrite=True)
print(f"✓ Model saved to: {model_dir}")

# %% [markdown]
# ---
# ## Step 7: Clustering and UMAP

# %%
print("\n" + "="*70)
print("Step 7: Clustering and UMAP")
print("="*70)

# Compute neighbors
print("\nComputing neighbors...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS)
print(f"✓ Neighbors computed (n={UMAP_N_NEIGHBORS})")

# UMAP
print("\nComputing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=42)
print(f"✓ UMAP computed (min_dist={UMAP_MIN_DIST}, spread={UMAP_SPREAD})")

# Clustering
print(f"\nClustering at {len(LEIDEN_RESOLUTIONS)} resolutions...")
for res in LEIDEN_RESOLUTIONS:
    key = f'leiden_res{res}'
    sc.tl.leiden(adata, resolution=res, key_added=key)
    n_clust = adata.obs[key].nunique()
    print(f"  Resolution {res}: {n_clust} clusters")

# Set default
adata.obs['leiden'] = adata.obs[f'leiden_res{DEFAULT_RESOLUTION}']
n_clusters = adata.obs['leiden'].nunique()
print(f"\n✓ Default clustering: resolution={DEFAULT_RESOLUTION}, n={n_clusters}")

# %% [markdown]
# ---
# ## Step 8: Batch Mixing Assessment

# %%
print("\n" + "="*70)
print("Step 8: Batch Mixing Assessment")
print("="*70)

def calculate_batch_mixing(adata, use_rep='X_scvi', batch_key='dataset', n_neighbors=50):
    knn = NearestNeighbors(n_neighbors=n_neighbors)
    knn.fit(adata.obsm[use_rep])
    _, indices = knn.kneighbors(adata.obsm[use_rep])
    
    entropies = []
    for idx in indices:
        batch_dist = adata.obs[batch_key].iloc[idx].value_counts(normalize=True)
        entropies.append(entropy(batch_dist))
    
    return np.mean(entropies), np.std(entropies)

mean_entropy, std_entropy = calculate_batch_mixing(adata, batch_key=BATCH_KEY)
max_entropy = np.log(n_batches)
mixing_ratio = mean_entropy / max_entropy * 100

print(f"\nBatch Mixing Metrics:")
print(f"  Entropy: {mean_entropy:.3f} ± {std_entropy:.3f}")
print(f"  Maximum: {max_entropy:.3f}")
print(f"  Mixing ratio: {mixing_ratio:.1f}%")

if mean_entropy > 2.0:
    print("  ✅ Excellent batch mixing")
elif mean_entropy > 1.5:
    print("  ✅ Good batch mixing")
elif mean_entropy > 1.0:
    print("  ⚠️  Moderate batch mixing")
else:
    print("  ❌ Limited batch mixing")

adata.uns['batch_mixing_entropy'] = mean_entropy

# %% [markdown]
# ---
# ## Step 9: Visualization

# %%
print("\n" + "="*70)
print("Step 9: Visualization")
print("="*70)

# UMAP: batch and clustering
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0],
           title='Batch Mixing', show=False, legend_loc='right margin', frameon=False)

sc.pl.umap(adata, color='leiden', ax=axes[1],
           title=f'Clustering (n={n_clusters})',
           show=False, legend_loc='on data', legend_fontsize=7, frameon=False)

plt.tight_layout()
plt.savefig(fig_dir / f'umap_batch_clustering.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ UMAP plots saved")

# %%
# Multi-resolution comparison
cluster_keys = [f'leiden_res{res}' for res in LEIDEN_RESOLUTIONS]

fig = sc.pl.umap(adata, color=cluster_keys, ncols=3,
                 return_fig=True, show=False, legend_loc='on data',
                 legend_fontsize=6, frameon=False)

plt.tight_layout()
plt.savefig(fig_dir / f'umap_multiresolution.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Multi-resolution plots saved")

# %%
# Marker genes
print("\nPlotting marker genes...")

# Set raw data for markers
adata.raw = adata_full[adata.obs_names, :].copy()

for category, genes in MARKER_GENES.items():
    available = [g for g in genes if g in adata.raw.var_names]
    if len(available) == 0:
        continue
    
    print(f"  {category}: {len(available)} genes")
    
    fig = sc.pl.umap(adata, color=available, use_raw=True,
                     ncols=3, vmax='p99', cmap='Reds',
                     return_fig=True, show=False, frameon=False)
    
    fig.suptitle(f'{category} Markers', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(fig_dir / f'markers_{category}.{FIGURE_FORMAT}',
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()

print("\n✓ Marker plots saved")

# %% [markdown]
# ---
# ## Step 10: Save Results

# %%
print("\n" + "="*70)
print("Step 10: Saving Results")
print("="*70)

# Transfer results to full dataset
adata_full.obsm['X_scvi'] = adata.obsm['X_scvi']
adata_full.obsm['X_umap'] = adata.obsm['X_umap']
adata_full.obsp['connectivities'] = adata.obsp['connectivities']
adata_full.obsp['distances'] = adata.obsp['distances']
adata_full.uns = adata.uns.copy()

# Transfer clustering
for col in adata.obs.columns:
    if 'leiden' in col:
        adata_full.obs[col] = adata.obs[col]

# Save
output_hvg = output_dir / "adata_tcell_scvi_hvg.h5ad"
output_full = output_dir / "adata_tcell_scvi_full.h5ad"

print(f"\nSaving datasets...")
adata.write_h5ad(output_hvg)
print(f"  HVG: {output_hvg}")

adata_full.write_h5ad(output_full)
print(f"  Full: {output_full}")

print("\n✓ All results saved")

# %% [markdown]
# ---
# ## Step 11: Summary

# %%
print("\n" + "="*70)
print("ANALYSIS SUMMARY")
print("="*70)

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\nDataset: {adata_full.shape[0]:,} cells × {adata_full.shape[1]:,} genes")
print(f"HVGs: {adata.shape[1]:,}")
print(f"Batches: {n_batches}")

print(f"\nscVI Configuration:")
print(f"  Latent dimensions: {SCVI_N_LATENT}")
print(f"  Training epochs: {len(train_elbo)}")
print(f"  Training time: {elapsed/60:.1f} min")

print(f"\nBatch Mixing:")
print(f"  Entropy: {mean_entropy:.3f}")
print(f"  Mixing ratio: {mixing_ratio:.1f}%")

print(f"\nClustering:")
print(f"  Resolution: {DEFAULT_RESOLUTION}")
print(f"  Clusters: {n_clusters}")

print(f"\nOutput Files:")
print(f"  {output_hvg}")
print(f"  {output_full}")
print(f"  {model_dir}")
print(f"  {fig_dir}/*.{FIGURE_FORMAT}")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)


