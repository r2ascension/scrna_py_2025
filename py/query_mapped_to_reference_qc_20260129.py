#!/usr/bin/env python3
"""
Robust marker diagnostic for scArches query mapping
Fixed version with extensive error handling
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
from pathlib import Path

# ============== CONFIG ==============
H5AD = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/query_mapped_to_reference.h5ad"
OUTDIR = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/markers_diagnostic"

GROUP_KEY = "cell_type_final"
MIN_CELLS = 100

# Key markers
MARKERS = [
    "EPCAM",   # Epithelial
    "KRT5",    # Basal
    "FOXJ1",   # Ciliated
    "MUC5AC",  # Goblet
    "CD3D",    # T
    "CD8A",    # CD8 T
    "NKG7",    # NK
    "MS4A1",   # B
    "MZB1",    # Plasma
    "LYZ",     # Myeloid
    "APOE",    # Macrophage
    "PECAM1",  # Endothelial
    "COL1A1",  # Fibroblast
    "ACTA2",   # SMC
]

# ============== MAIN ==============
print("="*60)
print("Marker Diagnostic for scArches Query Mapping")
print("="*60)

# Create output dir
os.makedirs(OUTDIR, exist_ok=True)

# Load
print(f"\n[1] Loading data...")
print(f"    {H5AD}")
adata = sc.read_h5ad(H5AD)
print(f"    Cells: {adata.n_obs:,}")
print(f"    Genes: {adata.n_vars:,}")

# Check group key
if GROUP_KEY not in adata.obs:
    print(f"ERROR: {GROUP_KEY} not found in obs!")
    print(f"Available keys: {adata.obs.columns.tolist()}")
    exit(1)

print(f"\n[2] Cell type distribution:")
print(adata.obs[GROUP_KEY].value_counts())

# Check confidence
if "mapping_confidence" in adata.obs:
    conf = adata.obs["mapping_confidence"]
    print(f"\n[3] Mapping quality:")
    print(f"    Median confidence: {conf.median():.4f}")
    print(f"    Mean confidence: {conf.mean():.4f}")
    low = (conf < 0.5).sum()
    print(f"    Low confidence (<0.5): {low} ({100*low/len(adata):.2f}%)")

# Check UMAP
print(f"\n[4] Checking UMAP coordinates...")
umap_key = None
for key in ["X_umap", "X_umap_mapped"]:
    if key in adata.obsm:
        umap_key = key
        print(f"    Found: {key}")
        break

if umap_key is None:
    print("    WARNING: No UMAP found! Will skip UMAP plots.")
else:
    # Ensure it's named X_umap for plotting
    if umap_key != "X_umap":
        adata.obsm["X_umap"] = adata.obsm[umap_key]
        print(f"    Copied {umap_key} -> X_umap")

# Check markers
print(f"\n[5] Checking marker genes...")
present = [g for g in MARKERS if g in adata.var_names]
missing = [g for g in MARKERS if g not in adata.var_names]
print(f"    Present: {len(present)}/{len(MARKERS)}")
if missing:
    print(f"    Missing: {missing}")

if len(present) == 0:
    print("ERROR: No marker genes found!")
    exit(1)

# ============== VISUALIZATION ==============
sc.settings.figdir = OUTDIR
sc.settings.set_figure_params(dpi=100, dpi_save=300)

# CRITICAL: Enable rasterization for large scatter plots
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts
matplotlib.rcParams['ps.fonttype'] = 42
sc.settings.set_figure_params(vector_friendly=False)  # Enable rasterization

# 1. UMAP by cell type
if umap_key:
    print(f"\n[6] UMAP by cell type...")
    try:
        fig = sc.pl.umap(
            adata,
            color=GROUP_KEY,
            legend_loc="on data",
            legend_fontsize=6,
            show=False,
            return_fig=True,
            frameon=False
        )
        # Rasterize scatter points
        for ax in fig.axes:
            for coll in ax.collections:
                coll.set_rasterized(True)
        fig.savefig(f"{OUTDIR}/umap_celltype.pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)
        print("    ✓ Saved: umap_celltype.pdf")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

# 2. UMAP by confidence
if umap_key and "mapping_confidence" in adata.obs:
    print(f"\n[7] UMAP by confidence...")
    try:
        fig = sc.pl.umap(
            adata,
            color="mapping_confidence",
            show=False,
            return_fig=True,
            frameon=False
        )
        # Rasterize scatter points
        for ax in fig.axes:
            for coll in ax.collections:
                coll.set_rasterized(True)
        fig.savefig(f"{OUTDIR}/umap_confidence.pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)
        print("    ✓ Saved: umap_confidence.pdf")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

# 3. Individual marker UMAPs
if umap_key:
    print(f"\n[8] Individual marker UMAPs...")
    
    # CRITICAL: Ensure we have log1p data for visualization
    # Check if .X is already log1p or if we need to use a layer
    is_log1p = adata.X.max() < 20  # heuristic: log1p data should be <20
    
    if is_log1p:
        print("    Using .X (already log1p)")
        layer_to_use = None
        use_raw_flag = False
    elif "log1p" in adata.layers:
        print("    Using layer='log1p'")
        layer_to_use = "log1p"
        use_raw_flag = False
    else:
        print("    WARNING: No log1p data found, will normalize on-the-fly")
        # Create temporary log1p layer for visualization
        if "counts" in adata.layers:
            adata.layers["log1p_temp"] = adata.layers["counts"].copy()
            sc.pp.normalize_total(adata, layer="log1p_temp", target_sum=1e4)
            sc.pp.log1p(adata, layer="log1p_temp")
            layer_to_use = "log1p_temp"
        else:
            adata.layers["log1p_temp"] = adata.X.copy()
            sc.pp.normalize_total(adata, layer="log1p_temp", target_sum=1e4)
            sc.pp.log1p(adata, layer="log1p_temp")
            layer_to_use = "log1p_temp"
        use_raw_flag = False
    
    for i, gene in enumerate(present[:8], 1):  # First 8 markers
        try:
            fig = sc.pl.umap(
                adata,
                color=gene,
                layer=layer_to_use,
                use_raw=use_raw_flag,
                show=False,
                return_fig=True,
                title=f"{gene} (log1p)",
                frameon=False,
                vmax="p99"  # Cap at 99th percentile for better visualization
            )
            # Rasterize scatter points
            for ax in fig.axes:
                for coll in ax.collections:
                    coll.set_rasterized(True)
            fig.savefig(f"{OUTDIR}/umap_{gene}.pdf", bbox_inches='tight', dpi=300)
            plt.close(fig)
            print(f"    ✓ [{i}/8] {gene}")
        except Exception as e:
            print(f"    ✗ [{i}/8] {gene}: {e}")

# 4. Dotplot - need to prepare data carefully
print(f"\n[9] Dotplot preparation...")
try:
    # Subset to markers
    adata_m = adata[:, present].copy()
    
    # Check data type
    print(f"    Data type: {adata_m.X.dtype}")
    print(f"    Data range: {adata_m.X.min():.2f} to {adata_m.X.max():.2f}")
    
    # Ensure numeric data
    if adata_m.X.dtype == object:
        print("    Converting to float32...")
        adata_m.X = adata_m.X.astype(np.float32)
    
    # If data looks like counts, normalize
    if adata_m.X.max() > 20:
        print("    Normalizing (looks like counts)...")
        sc.pp.normalize_total(adata_m, target_sum=1e4)
        sc.pp.log1p(adata_m)
    
    # Filter small groups
    vc = adata_m.obs[GROUP_KEY].value_counts()
    keep_groups = vc[vc >= MIN_CELLS].index.tolist()
    print(f"    Filtering groups (>={MIN_CELLS} cells): {len(keep_groups)}/{len(vc)}")
    
    adata_m = adata_m[adata_m.obs[GROUP_KEY].isin(keep_groups)].copy()
    adata_m.obs[GROUP_KEY] = adata_m.obs[GROUP_KEY].astype('category')
    
    # Generate dotplot
    print("    Generating dotplot...")
    
    # Method 1: Try with dendrogram
    try:
        sc.tl.dendrogram(adata_m, groupby=GROUP_KEY)
        fig = sc.pl.dotplot(
            adata_m,
            var_names=present,
            groupby=GROUP_KEY,
            dendrogram=True,
            standard_scale="var",
            show=False,
            return_fig=True
        )
        fig.savefig(f"{OUTDIR}/dotplot_markers.pdf", bbox_inches='tight')
        plt.close(fig)
        print("    ✓ Saved: dotplot_markers.pdf")
    except Exception as e1:
        print(f"    ✗ Dendrogram method failed: {e1}")
        
        # Method 2: Try without dendrogram
        try:
            fig = sc.pl.dotplot(
                adata_m,
                var_names=present,
                groupby=GROUP_KEY,
                standard_scale="var",
                show=False,
                return_fig=True
            )
            fig.savefig(f"{OUTDIR}/dotplot_markers_nodendrogram.pdf", bbox_inches='tight')
            plt.close(fig)
            print("    ✓ Saved: dotplot_markers_nodendrogram.pdf")
        except Exception as e2:
            print(f"    ✗ Simple dotplot failed: {e2}")

except Exception as e:
    print(f"    ✗ Dotplot preparation failed: {e}")

# 5. Violin plot as alternative
print(f"\n[10] Violin plot (alternative to dotplot)...")
try:
    # Use subset for speed
    adata_v = adata[:, present[:6]].copy()
    
    # Filter groups
    vc = adata_v.obs[GROUP_KEY].value_counts()
    keep = vc[vc >= MIN_CELLS].index.tolist()[:10]  # Top 10 groups
    adata_v = adata_v[adata_v.obs[GROUP_KEY].isin(keep)].copy()
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, gene in enumerate(present[:6]):
        ax = axes[idx]
        sc.pl.violin(
            adata_v,
            keys=gene,
            groupby=GROUP_KEY,
            rotation=90,
            ax=ax,
            show=False
        )
        ax.set_title(gene)
    
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/violin_markers.pdf", bbox_inches='tight')
    plt.close()
    print("    ✓ Saved: violin_markers.pdf")
    
except Exception as e:
    print(f"    ✗ Violin plot failed: {e}")

# ============== SUMMARY ==============
print("\n" + "="*60)
print("✅ Diagnostic complete!")
print("="*60)
print(f"\n📁 Output directory: {OUTDIR}")
print("\nGenerated files:")
for f in sorted(Path(OUTDIR).glob("*.pdf")):
    print(f"  - {f.name}")

print("\n📊 Next steps:")
print("  1. Check umap_celltype.pdf for cell type distribution")
print("  2. Check umap_confidence.pdf for mapping quality")
print("  3. Check individual marker UMAPs")
print("  4. Check dotplot or violin for marker specificity")