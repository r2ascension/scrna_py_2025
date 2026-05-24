#!/usr/bin/env python3
"""
SCVI-based integration pipeline (linear version, enriched).

All parameters are explicit below for easy editing.
Run directly: `python tcell_scvi_integration.py`.

Adds:
- Optional preprocessing + HVG selection
- Optional categorical covariate (e.g., tissue)
- Training history plot & model save
- Multi-resolution clustering, marker plots, batch heatmap
- Summary CSV with key metrics

Notes:
- Prefer integer raw counts in `layers['counts']` for SCVI. If absent, falls back to .X.
- Requires scvi-tools (`pip install scvi-tools`). Uses GPU automatically if available.
"""

import sys
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import scanpy as sc
except Exception:
    print("scanpy is required. Install with: pip install scanpy", file=sys.stderr)
    raise

try:
    import scvi
except Exception:
    print("scvi-tools is required. Install with: pip install scvi-tools", file=sys.stderr)
    raise


# =====================
# Configuration (edit!)
# =====================
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1203/scvi_integration"

BATCH_KEY = "dataset"            # obs column for batches
COUNTS_LAYER = "counts"           # set to None to force using .X
CATEGORICAL_COVARIATE = "tissue"  # set to None to disable

N_LATENT = 30                     # SCVI latent dim
N_LAYERS = 2                      # NN layers
DROPOUT_RATE = 0.1                # Dropout rate
EPOCHS = 300                      # training epochs (with early stopping plot)
LEARNING_RATE = 1e-3              # LR for training (used in plan_kwargs)

N_NEIGHBORS = 30                  # graph neighbors
UMAP_MIN_DIST = 0.3               # UMAP min_dist
UMAP_SPREAD = 1.3                 # UMAP spread
RESOLUTIONS = [2.0, 2.5, 3.0, 3.5]  # Leiden resolutions to run
DEFAULT_RESOLUTION = 3.0

SAVE_H5AD = True
QUIET = False


# =============
# Print config
# =============
inp = Path(INPUT_H5AD).expanduser().resolve()
outdir = Path(OUTPUT_DIR).expanduser().resolve()
outdir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("SCVI INTEGRATION PIPELINE (linear)")
print("=" * 70)
print(f"Input:        {inp}")
print(f"Output dir:   {outdir}")
print(f"Batch key:    {BATCH_KEY}")
print(f"Counts layer: {COUNTS_LAYER}")
print(f"Latent dim:   {N_LATENT}; Epochs: {EPOCHS}; LR: {LEARNING_RATE}")
print(f"Neighbors:    {N_NEIGHBORS}; UMAP min_dist: {UMAP_MIN_DIST}")
print(f"UMAP spread:  {UMAP_SPREAD}")
print(f"Resolutions:  {RESOLUTIONS}; default={DEFAULT_RESOLUTION}")
print("=" * 70)

if QUIET:
    sc.settings.verbosity = 1
    warnings.filterwarnings("ignore")
else:
    sc.settings.verbosity = 2

if not inp.exists():
    print(f"Input not found: {inp}", file=sys.stderr)
    sys.exit(2)


# ===============
# 1) Load AnnData
# ===============
print("\n[1/9] Loading data...")
adata = sc.read_h5ad(str(inp))
print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

if BATCH_KEY not in adata.obs.columns:
    print(f"Error: obs missing batch key '{BATCH_KEY}'", file=sys.stderr)
    sys.exit(2)
print(f"Batches in '{BATCH_KEY}': {adata.obs[BATCH_KEY].nunique()}")
has_covariate = CATEGORICAL_COVARIATE is not None and CATEGORICAL_COVARIATE in adata.obs.columns
if has_covariate:
    print(f"Categorical covariate available: '{CATEGORICAL_COVARIATE}'")
else:
    if CATEGORICAL_COVARIATE:
        print(f"Warning: '{CATEGORICAL_COVARIATE}' not found in obs; proceeding without categorical covariate")


# ======================
# 2) Setup for SCVI
# ======================
print("\n[2/9] Preprocessing for scVI (optional HVG)...")
# Prepare counts layer if missing
if COUNTS_LAYER and COUNTS_LAYER not in adata.layers:
    print("Counts layer missing; creating from .X or .raw...")
    if adata.raw is not None and hasattr(adata.raw, "X"):
        adata.layers[COUNTS_LAYER] = adata.raw.X.copy()
    else:
        adata.layers[COUNTS_LAYER] = adata.X.copy()

# Keep a copy for marker visualization
adata_original = adata.copy()

# Basic gene filtering + HVG selection (batch-aware)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
try:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=4000,
        flavor="seurat_v3",
        batch_key=BATCH_KEY,
        subset=True,
    )
    print(f"Selected HVGs: {adata.n_vars}")
except Exception:
    print("Skipping HVG selection due to error; proceeding with current genes")

print("\n[3/9] Setting up SCVI...")
use_layer = None
if COUNTS_LAYER and COUNTS_LAYER in adata.layers:
    use_layer = COUNTS_LAYER
    print(f"Using counts from layer '{use_layer}' for SCVI.")
else:
    print("Warning: counts layer not found; SCVI will use .X (may be log-normalized).")

if has_covariate:
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=BATCH_KEY,
        layer=use_layer,
        categorical_covariate_keys=[CATEGORICAL_COVARIATE],
    )
else:
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=BATCH_KEY,
        layer=use_layer,
    )


# ======================
# 3) Train SCVI model
# ======================
print("\n[4/9] Creating & training SCVI model...")
model = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    dropout_rate=DROPOUT_RATE,
)
import time
start_time = time.time()
try:
    model.train(max_epochs=EPOCHS, plan_kwargs={"lr": LEARNING_RATE}, early_stopping=True, early_stopping_patience=15)
except TypeError:
    model.train(max_epochs=EPOCHS)
elapsed = time.time() - start_time
print(f"✓ Training finished in {elapsed:.1f}s ({elapsed/60:.1f} min)")

# Try plotting training history if available
try:
    hist_train = model.history_["elbo_train"]
    hist_val = model.history_.get("elbo_validation", [])
    plt.figure(figsize=(8, 4))
    plt.plot(hist_train, label="Training ELBO", alpha=0.7)
    if len(hist_val) > 0:
        plt.plot(hist_val, label="Validation ELBO", alpha=0.7)
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("scVI Training History")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "scvi_training_history.pdf", dpi=150, bbox_inches="tight")
    plt.close()
    print("✓ Training history plot saved")
except Exception:
    print("Training history not available for plotting")


# ==================================
# 4) Latent representation (X_scVI)
# ==================================
print("\n[5/9] Extracting latent representation...")
adata.obsm["X_scVI"] = model.get_latent_representation()
print(f"Latent shape: {adata.obsm['X_scVI'].shape}")


# ================================
# 5) Neighbors graph and UMAP
# ================================
print("\n[6/9] Building neighbors / UMAP...")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBORS, random_state=0)
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=0)
print("✓ Neighbors/UMAP done")


# ======================================
# 6) Leiden clustering for RESOLUTIONS
# ======================================
print(f"\n[7/9] Leiden clustering at resolutions: {RESOLUTIONS}")
for res in RESOLUTIONS:
    key = f"leiden_scvi_res{res}"
    sc.tl.leiden(adata, resolution=float(res), key_added=key)
    ncl = adata.obs[key].nunique()
    print(f"  res={res}: {ncl} clusters")

# Set default clustering
adata.obs["leiden_scvi"] = adata.obs[f"leiden_scvi_res{DEFAULT_RESOLUTION}"]
n_clusters_default = adata.obs["leiden_scvi"].nunique()
print(f"✓ Default clustering: resolution={DEFAULT_RESOLUTION}, n_clusters={n_clusters_default}")


# ===================
# Save UMAP figures
# ===================
print("\n[8/9] Saving figures...")
sc.set_figure_params(dpi=100, facecolor="white")

# UMAP by batch
fig1 = outdir / "umap_by_batch.pdf"
sc.pl.umap(adata, color=BATCH_KEY, legend_loc="none", show=False)
plt.tight_layout()
plt.savefig(fig1, dpi=150, bbox_inches="tight")
plt.close()

# UMAP by first resolution
first_res = float(RESOLUTIONS[0])
key0 = f"leiden_scvi_res{first_res}"
fig2 = outdir / f"umap_leiden_res{first_res}.pdf"
sc.pl.umap(adata, color=key0, legend_loc="on data", legend_fontsize=6, show=False)
plt.tight_layout()
plt.savefig(fig2, dpi=150, bbox_inches="tight")
plt.close()

# Marker genes (optional quick check)
KEY_MARKERS = ["NCAM1", "CD4", "CD8A", "CD3D", "GZMB", "FOXP3"]
available_markers = [m for m in KEY_MARKERS if m in adata_original.var_names]
if len(available_markers) > 0:
    print(f"Plotting marker genes: {available_markers}")
    adata.raw = adata_original[adata.obs_names, :].copy()
    fig = sc.pl.umap(
        adata,
        color=available_markers[:6],
        use_raw=True,
        ncols=3,
        vmax="p99",
        cmap="Reds",
        return_fig=True,
        show=False,
    )
    fig.savefig(outdir / "umap_markers.pdf", dpi=150, bbox_inches="tight")
    plt.close()
    print("✓ Marker gene plots saved")
else:
    print("No key markers found for marker plots")


# ==========================
# Batch mixing entropy (kNN)
# ==========================
print("[9/9] Computing batch mixing entropy (kNN on X_scVI)...")
from sklearn.neighbors import NearestNeighbors
from scipy.stats import entropy

X = adata.obsm["X_scVI"]
knn = NearestNeighbors(n_neighbors=min(50, X.shape[0]))
knn.fit(X)
_, indices = knn.kneighbors(X)

batches = adata.obs[BATCH_KEY].astype("category")
entropies = []
for idx in indices:
    dist = batches.iloc[idx].value_counts(normalize=True)
    entropies.append(entropy(dist))
entropies = np.asarray(entropies, dtype=float)
mean_ent = float(np.mean(entropies))
std_ent = float(np.std(entropies))

metrics = []
for res in RESOLUTIONS:
    key = f"leiden_scvi_res{float(res)}"
    ncl = int(adata.obs[key].nunique())
    metrics.append({
        "resolution": float(res),
        "n_clusters": ncl,
        "batch_mixing_entropy": mean_ent,
        "entropy_std": std_ent,
    })

metrics_df = pd.DataFrame(metrics)
metrics_path = outdir / "scvi_metrics.csv"
metrics_df.to_csv(metrics_path, index=False)
print(f"✓ Metrics saved to: {metrics_path}")

# Batch distribution per cluster heatmap
try:
    import seaborn as sns
    batch_cluster = pd.crosstab(
        adata.obs["leiden_scvi"],
        adata.obs[BATCH_KEY],
        normalize="index",
    )
    plt.figure(figsize=(12, max(6, n_clusters_default * 0.3)))
    sns.heatmap(batch_cluster, cmap="viridis", cbar_kws={"label": "Proportion"})
    plt.xlabel("Dataset (Batch)")
    plt.ylabel("Cluster")
    plt.title("Batch Distribution per Cluster")
    plt.tight_layout()
    plt.savefig(outdir / "batch_per_cluster.pdf", dpi=150, bbox_inches="tight")
    plt.close()
    print("✓ Batch distribution heatmap saved")
    batch_dominant = (batch_cluster > 0.5).sum(axis=1)
    problematic_clusters = batch_dominant[batch_dominant > 0].index.tolist()
except Exception:
    problematic_clusters = []


# ========================
# Save integrated AnnData
# ========================
if SAVE_H5AD:
    out_h5ad = outdir / "adata_scvi_integrated.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"✓ Integrated AnnData saved to: {out_h5ad}")

# Summary CSV
summary = {
    "method": "scVI",
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_batches": int(adata.obs[BATCH_KEY].nunique()),
    "n_latent": int(N_LATENT),
    "default_resolution": float(DEFAULT_RESOLUTION),
    "n_clusters_default": int(n_clusters_default),
    "batch_mixing_entropy": float(mean_ent),
    "batch_mixing_std": float(std_ent),
    "umap_min_dist": float(UMAP_MIN_DIST),
    "umap_spread": float(UMAP_SPREAD),
    "has_covariate": bool(has_covariate),
}
pd.DataFrame([summary]).to_csv(outdir / "scvi_summary.csv", index=False)
print(f"✓ Summary saved to: {outdir / 'scvi_summary.csv'}")


print("\n" + "=" * 70)
print("SCVI integration completed")
print("Figures:")
print(f"  - {fig1}")
print(f"  - {fig2}")
print(f"Metrics: {metrics_path}")
if SAVE_H5AD:
    print(f"AnnData: {out_h5ad}")
if len(available_markers) > 0:
    print(f"Markers: {outdir / 'umap_markers.pdf'}")
print(f"Batch heatmap: {outdir / 'batch_per_cluster.pdf'}")
print("=" * 70)

