# -*- coding: utf-8 -*-
# %% [markdown]
# # Epithelial Subcluster Visualization with L3 Annotations v1.0
#
# **Purpose:**
# - Apply refined L3 cell type annotations to epithelial subclusters
# - Generate comprehensive publication-quality visualizations
# - Create marker dotplots, heatmaps, and UMAPs
#
# **Author:** r2end
# **Date:** 2025-02-05
# **Version:** 1.0

# %% [markdown]
# ## Configuration

# %%
import os
import sys
import gc
import warnings
from pathlib import Path
from datetime import datetime
import time

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
from matplotlib.patches import Rectangle

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# =============================================================================
# PATHS
# =============================================================================
INPUT_H5AD = "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad"
MARKER_DIR = Path("/home/h2048/data/R/0129/epithelial_interpret_v2_7_FIXED")
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/epithelial_viz_L3_v1_0")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# L3 annotation mapping CSV (from your document 4)
L3_MAPPING_CSV = OUTPUT_DIR / "cluster_to_L3_mapping.csv"

# =============================================================================
# MARKER DEFINITIONS
# =============================================================================
# Core markers for each L3 category (from your document 4)
L3_MARKER_PANELS = {
    # Alveolar
    "AT1_Canonical": ["AGER", "HOPX", "CAV1", "AQP4", "RTKN2", "CLDN18", "EMP2"],
    "AT1_MatrixRemodeling": ["AGER", "CAV1", "SPARC", "COL4A1", "COL4A2", "SPOCK2"],
    "AT2_Canonical": ["SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NAPSA", "SLC34A2"],
    "AT2_Inflammatory_Repair": ["SFTPC", "ABCA3", "CHI3L1", "CXCL8", "SAA1", "LCN2"],
    "Epithelial_Cycling": ["MKI67", "TOP2A", "UBE2C", "BIRC5", "AURKB", "CENPA", "CCNB1"],
    
    # Basal
    "Basal_Progenitor": ["KRT5", "KRT14", "TP63", "KRT15", "KRT19", "ITGA6", "NGFR"],
    "Basal_Cycling": ["KRT14", "KRT5", "TP63", "MKI67", "TOP2A", "BIRC5"],
    "Basal_Inflammatory": ["KRT17", "CXCL8", "CXCL1", "CXCL2", "TNFAIP3", "FOS", "JUN"],
    "Basal_EMT_ECM": ["KRT14", "TP63", "NGFR", "FN1", "COL17A1", "MMP2", "VIM"],
    
    # Ciliated
    "Ciliated_Mature": ["FOXJ1", "TPPP3", "DNAH5", "DNAH9", "RSPH1", "RFX2", "RFX3"],
    "Ciliogenesis_Deuterosomal": ["DEUP1", "CCNO", "FOXN4", "MCIDAS", "CDC20B", "E2F7", "PLK4"],
    "Ciliated_Cycling_Immature": ["TPPP3", "RSPH1", "MKI67", "TOP2A", "FOXN4"],
    
    # Secretory
    "Secretory_Club": ["SCGB1A1", "SCGB3A1", "SCGB3A2", "AGR2", "AGR3", "CYP2F1"],
    "Secretory_Club_AT2_Transitional": ["SCGB1A1", "SCGB3A1", "SFTPB", "NAPSA", "GPR116", "CLDN18"],
    "Goblet_Mucin": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
    "Goblet_Defense_DUOX2": ["DUOX2", "DUOXA2", "LCN2", "BPIFA2", "CEACAM5"],
    
    # SMG
    "SMG_Serous": ["LTF", "LYZ", "SLPI", "DMBT1", "BPIFA1", "AZGP1", "WFDC2"],
    "SMG_Duct_Secretory_Defense": ["PIGR", "SCGB3A1", "TCN1", "WFDC2", "DMBT1", "SLPI"],
    
    # Special
    "Squamous_Metaplasia": ["SPRR1A", "SPRR2A", "SPRR2E", "IVL", "KRT6A", "KLK7", "S100A7"],
    "Ionocyte_Brush": ["FOXI1", "ASCL3", "CFTR", "ATP6V0D2", "CLCNKA", "CLCNKB", "BSND"],
    "Mesenchymal_Contaminant": ["COL4A1", "COL4A2", "LAMA1", "DCN", "COL1A1", "COL1A2"]
}

# =============================================================================
# VISUALIZATION PARAMETERS
# =============================================================================
FIGURE_DPI = 300
FIGURE_FORMAT = "pdf"
UMAP_SIZE = 3
UMAP_ALPHA = 0.6

# Color schemes
L3_COLORS = {
    # Alveolar - Blues/Cyans
    "AT1_Canonical": "#1f77b4",
    "AT1_MatrixRemodeling": "#aec7e8",
    "AT2_Canonical": "#17becf",
    "AT2_Inflammatory_Repair": "#9edae5",
    "Epithelial_Cycling": "#c49c94",
    
    # Basal - Reds/Oranges
    "Basal_Progenitor": "#d62728",
    "Basal_Cycling": "#ff7f0e",
    "Basal_Inflammatory": "#ff9896",
    "Basal_EMT_ECM": "#ffbb78",
    
    # Ciliated - Greens
    "Ciliated_Mature": "#2ca02c",
    "Ciliogenesis_Deuterosomal": "#98df8a",
    "Ciliated_Cycling_Immature": "#8c564b",
    
    # Secretory - Purples/Pinks
    "Secretory_Club": "#9467bd",
    "Secretory_Club_AT2_Transitional": "#c5b0d5",
    "Goblet_Mucin": "#e377c2",
    "Goblet_Defense_DUOX2": "#f7b6d2",
    
    # SMG - Yellows/Browns
    "SMG_Serous": "#bcbd22",
    "SMG_Duct_Secretory_Defense": "#dbdb8d",
    
    # Special - Grays
    "Squamous_Metaplasia": "#7f7f7f",
    "Ionocyte_Brush": "#c7c7c7",
    "Mesenchymal_Contaminant": "#3f3f3f"
}

# =============================================================================
# REPRODUCIBILITY
# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

sc.settings.figdir = OUTPUT_DIR / "figures"
sc.settings.figdir.mkdir(exist_ok=True)

PIPELINE_START = time.time()

print("=" * 80)
print("Epithelial Subcluster Visualization Pipeline v1.0")
print("=" * 80)
print(f"Input:  {INPUT_H5AD}")
print(f"Markers: {MARKER_DIR}")
print(f"Output: {OUTPUT_DIR}")
print(f"Date:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# %% [markdown]
# ## Step 1: Create L3 Mapping File

# %%
print("\n" + "=" * 80)
print("STEP 1: Creating L3 Annotation Mapping")
print("=" * 80)

# Mapping from document 4 (corrected based on actual analysis)
cluster_to_l3_mapping = {
    # AT lineage
    "AT1 _0": "AT1_Canonical",
    "AT1 _1": "AT1_MatrixRemodeling",
    "AT1 _2": "AT2_Canonical",  # Actually AT2-like based on markers
    
    "AT2_0": "AT2_Canonical",
    "AT2_1": "AT2_Canonical",
    "AT2_2": "AT2_Inflammatory_Repair",
    "AT2_3": "Epithelial_Cycling",
    
    # Basal lineage
    "Basal_0": "Basal_EMT_ECM",
    "Basal_1": "Basal_Cycling",
    "Basal_2": "Basal_Inflammatory",
    "Basal_3": "Goblet_Mucin",
    "Basal_4": "Goblet_Mucin",
    "Basal_5": "Basal_Progenitor",
    
    # Ciliated lineage
    "Ciliated_0": "Ciliated_Mature",
    "Ciliated_1": "Ciliated_Mature",
    "Ciliated_2": "Ciliated_Cycling_Immature",
    "Ciliated_3": "Ciliated_Mature",
    "Ciliated_4": "Ciliated_Mature",
    "Ciliated_5": "Goblet_Mucin",
    
    "Deuterosomal_0": "Ciliogenesis_Deuterosomal",
    "Deuterosomal_1": "Ciliogenesis_Deuterosomal",
    
    # Dividing basal
    "Dividing_Basal_0": "Basal_Cycling",
    "Dividing_Basal_1": "Basal_Cycling",
    
    # Ionocyte
    "Ionocyte_n_Brush_0": "Ionocyte_Brush",
    "Ionocyte_n_Brush_1": "Ionocyte_Brush",
    
    # SMG lineage
    "SMG_Basal_0": "Basal_EMT_ECM",
    "SMG_Basal_1": "Mesenchymal_Contaminant",
    "SMG_Basal_2": "Basal_EMT_ECM",
    
    "SMG_Duct_0": "SMG_Duct_Secretory_Defense",
    "SMG_Duct_1": "Squamous_Metaplasia",
    "SMG_Duct_2": "Squamous_Metaplasia",
    "SMG_Duct_3": "Squamous_Metaplasia",
    "SMG_Duct_4": "Squamous_Metaplasia",
    
    "SMG_Mucous_0": "Goblet_Mucin",
    "SMG_Mucous_1": "Ionocyte_Brush",  # Verify this with POU2F3/TRPM5
    
    "SMG_Serous_0": "SMG_Serous",
    "SMG_Serous_1": "SMG_Serous",
    "SMG_Serous_2": "SMG_Serous",
    
    # Secretory
    "Secretory_Club_0": "Secretory_Club_AT2_Transitional",
    
    "Secretory_Goblet_0": "Goblet_Defense_DUOX2",
    "Secretory_Goblet_1": "Secretory_Club",
    "Secretory_Goblet_2": "Basal_Progenitor",
    "Secretory_Goblet_3": "Squamous_Metaplasia",
    "Secretory_Goblet_4": "Ciliated_Cycling_Immature",
    
    # Suprabasal
    "Suprabasal_0": "Basal_Progenitor",
    "Suprabasal_1": "Basal_Cycling",
    "Suprabasal_2": "Squamous_Metaplasia",
    "Suprabasal_3": "Squamous_Metaplasia"
}

# Save mapping to CSV
mapping_df = pd.DataFrame(list(cluster_to_l3_mapping.items()), 
                         columns=["Cluster", "cell_type_L3_refined"])
mapping_df.to_csv(L3_MAPPING_CSV, index=False)

print(f"[OK] Created L3 mapping: {L3_MAPPING_CSV.name}")
print(f"  Total clusters: {len(cluster_to_l3_mapping)}")
print(f"  Unique L3 labels: {len(set(cluster_to_l3_mapping.values()))}")
print("\nL3 label distribution:")
l3_counts = pd.Series(cluster_to_l3_mapping.values()).value_counts()
for label, count in l3_counts.head(15).items():
    print(f"  {label}: {count} clusters")

# %% [markdown]
# ## Step 2: Load Data and Apply L3 Annotations

# %%
print("\n" + "=" * 80)
print("STEP 2: Loading Data and Applying L3 Annotations")
print("=" * 80)

t0 = time.time()
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[OK] Loaded in {time.time() - t0:.1f}s")

print(f"\nDataset summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Has .raw: {adata.raw is not None}")

# Check subcluster column
if 'subcluster' not in adata.obs.columns:
    raise ValueError("'subcluster' column not found in adata.obs!")

print(f"\nOriginal subcluster distribution:")
subcluster_counts = adata.obs['subcluster'].value_counts()
print(f"  Total subclusters: {len(subcluster_counts)}")
print(f"  Range: {subcluster_counts.min()} - {subcluster_counts.max()} cells")

# Apply L3 annotations
print("\n[INFO] Applying L3 annotations...")

# Normalize subcluster labels to avoid hidden whitespace issues
subcluster_clean = adata.obs['subcluster'].astype(str).str.strip()
n_changed = (subcluster_clean != adata.obs['subcluster'].astype(str)).sum()
if n_changed > 0:
    print(f"[WARN] Found {n_changed} subcluster labels with leading/trailing spaces. Using stripped values for mapping.")

# Normalize mapping keys as well
cluster_to_l3_mapping_clean = {str(k).strip(): v for k, v in cluster_to_l3_mapping.items()}
adata.obs['cell_type_L3'] = subcluster_clean.map(cluster_to_l3_mapping_clean)

# Check for unmapped clusters
unmapped = adata.obs['cell_type_L3'].isna().sum()
if unmapped > 0:
    print(f"[WARN] {unmapped} cells with unmapped subclusters:")
    unmapped_clusters = adata.obs.loc[
        adata.obs['cell_type_L3'].isna(),
        'subcluster'
    ].astype(str).str.strip().unique()
    for cluster in unmapped_clusters:
        n_cells = (subcluster_clean == cluster).sum()
        print(f"  {cluster}: {n_cells} cells")
    
    # Fill unmapped with cleaned subcluster name
    adata.obs.loc[adata.obs['cell_type_L3'].isna(), 'cell_type_L3'] = \
        subcluster_clean[adata.obs['cell_type_L3'].isna()]
    print("[INFO] Filled unmapped cells with original subcluster names")

# Convert to categorical
adata.obs['cell_type_L3'] = adata.obs['cell_type_L3'].astype('category')

print(f"\n[OK] L3 annotations applied")
l3_counts = adata.obs['cell_type_L3'].value_counts()
print(f"  Total L3 labels: {len(l3_counts)}")
print("\nL3 cell type distribution:")
for label, count in l3_counts.head(20).items():
    pct = count / adata.n_obs * 100
    print(f"  {label}: {count:,} cells ({pct:.1f}%)")
if len(l3_counts) > 20:
    print(f"  ... and {len(l3_counts) - 20} more")

# %% [markdown]
# ## Step 3: Load Marker Data

# %%
print("\n" + "=" * 80)
print("STEP 3: Loading Marker Gene Data")
print("=" * 80)

# Load marker files
markers_all = MARKER_DIR / "all_markers.csv"
markers_top30 = MARKER_DIR / "top30_per_cluster.csv"
markers_filtered = MARKER_DIR / "top_markers_filtered.csv"

if not markers_all.exists():
    raise FileNotFoundError(f"Marker file not found: {markers_all}")

df_markers = pd.read_csv(markers_all)
print(f"[OK] Loaded markers: {markers_all.name}")
print(f"  Total genes: {len(df_markers):,}")
print(f"  Clusters: {df_markers['cluster'].nunique()}")

if 'avg_log2FC' in df_markers.columns:
    print(f"  Log2FC range: [{df_markers['avg_log2FC'].min():.2f}, {df_markers['avg_log2FC'].max():.2f}]")
if 'p_val_adj' in df_markers.columns:
    sig_markers = (df_markers['p_val_adj'] < 0.05).sum()
    print(f"  Significant markers (padj < 0.05): {sig_markers:,}")

# Load filtered markers if available
if markers_filtered.exists():
    df_markers_filt = pd.read_csv(markers_filtered)
    print(f"\n[OK] Loaded filtered markers: {markers_filtered.name}")
    print(f"  Total genes: {len(df_markers_filt):,}")
else:
    print(f"\n[WARN] Filtered markers not found: {markers_filtered.name}")
    df_markers_filt = df_markers.copy()

# %% [markdown]
# ## Step 4: Generate L3 Overview UMAP

# %%
print("\n" + "=" * 80)
print("STEP 4: Generating L3 Overview UMAP")
print("=" * 80)

if 'X_umap' not in adata.obsm:
    print("[WARN] No UMAP found in adata.obsm, computing...")
    if 'X_pca' not in adata.obsm:
        print("[INFO] Computing PCA first...")
        sc.pp.pca(adata, n_comps=50, random_state=RANDOM_SEED)
    sc.pp.neighbors(adata, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    print("[OK] UMAP computed")

# Create color mapping
l3_labels = sorted(adata.obs['cell_type_L3'].unique())
n_labels = len(l3_labels)

# Use predefined colors where available, generate for others
color_map = {}
for label in l3_labels:
    if label in L3_COLORS:
        color_map[label] = L3_COLORS[label]
    else:
        # Generate color for unmapped labels
        idx = l3_labels.index(label)
        color_map[label] = plt.cm.tab20(idx % 20)

adata.uns['cell_type_L3_colors'] = [color_map[label] for label in adata.obs['cell_type_L3'].cat.categories]

# Plot L3 UMAP
fig, ax = plt.subplots(figsize=(14, 10))
sc.pl.umap(
    adata,
    color='cell_type_L3',
    ax=ax,
    show=False,
    legend_loc='right margin',
    legend_fontsize=8,
    frameon=False,
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    title='Epithelial Cells - L3 Annotations'
)
plt.tight_layout()
output_file = sc.settings.figdir / f'01_umap_L3_overview.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] Saved: {output_file.name}")

# Plot by major lineage (color-coded subsets)
print("\n[INFO] Generating lineage-specific UMAPs...")

lineage_groups = {
    'Alveolar': ['AT1_Canonical', 'AT1_MatrixRemodeling', 'AT2_Canonical', 'AT2_Inflammatory_Repair'],
    'Basal': ['Basal_Progenitor', 'Basal_Cycling', 'Basal_Inflammatory', 'Basal_EMT_ECM'],
    'Ciliated': ['Ciliated_Mature', 'Ciliogenesis_Deuterosomal', 'Ciliated_Cycling_Immature'],
    'Secretory': ['Secretory_Club', 'Secretory_Club_AT2_Transitional', 'Goblet_Mucin', 'Goblet_Defense_DUOX2'],
    'SMG': ['SMG_Serous', 'SMG_Duct_Secretory_Defense'],
    'Other': ['Squamous_Metaplasia', 'Ionocyte_Brush', 'Epithelial_Cycling', 'Mesenchymal_Contaminant']
}

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, (lineage, labels) in enumerate(lineage_groups.items()):
    # Create mask for cells in this lineage
    mask = adata.obs['cell_type_L3'].isin(labels)
    
    # Highlight cells in lineage
    adata.obs[f'is_{lineage}'] = 'Other'
    adata.obs.loc[mask, f'is_{lineage}'] = adata.obs.loc[mask, 'cell_type_L3']
    
    # Plot
    sc.pl.umap(
        adata,
        color=f'is_{lineage}',
        ax=axes[idx],
        show=False,
        frameon=False,
        size=UMAP_SIZE * 0.7,
        alpha=UMAP_ALPHA * 0.8,
        title=f'{lineage} Lineage',
        legend_loc='none' if idx < 5 else 'right margin',
        legend_fontsize=6
    )

plt.tight_layout()
output_file = sc.settings.figdir / f'02_umap_L3_by_lineage.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] Saved: {output_file.name}")

# Clean up temporary columns
for lineage in lineage_groups.keys():
    if f'is_{lineage}' in adata.obs.columns:
        adata.obs.drop(columns=f'is_{lineage}', inplace=True)

# %% [markdown]
# ## Step 5: Generate Marker Dotplots

# %%
print("\n" + "=" * 80)
print("STEP 5: Generating Marker Dotplots")
print("=" * 80)

# Check which markers are available in dataset
def get_available_markers(marker_dict, var_names):
    """Get markers that are actually in the dataset"""
    available = {}
    missing_all = []
    
    for cell_type, markers in marker_dict.items():
        available_for_type = [m for m in markers if m in var_names]
        missing = [m for m in markers if m not in var_names]
        
        if len(available_for_type) > 0:
            available[cell_type] = available_for_type
        if len(missing) > 0:
            missing_all.extend(missing)
    
    return available, list(set(missing_all))

# Determine gene universe (prefer .raw)
if adata.raw is not None:
    gene_universe = adata.raw.var_names
    use_raw_for_plot = True
    print(f"[INFO] Using .raw for markers ({len(gene_universe):,} genes)")
else:
    gene_universe = adata.var_names
    use_raw_for_plot = False
    print(f"[INFO] Using .X for markers ({len(gene_universe):,} genes)")

available_markers, missing_markers = get_available_markers(L3_MARKER_PANELS, gene_universe)

print(f"\n[INFO] Marker availability:")
print(f"  L3 categories with markers: {len(available_markers)}/{len(L3_MARKER_PANELS)}")
print(f"  Total available markers: {sum(len(v) for v in available_markers.values())}")
print(f"  Missing markers: {len(missing_markers)}")

if len(missing_markers) > 0 and len(missing_markers) < 20:
    print(f"  Missing: {', '.join(sorted(missing_markers)[:20])}")
    if len(missing_markers) > 20:
        print(f"  ... and {len(missing_markers) - 20} more")

# Generate dotplot for each major lineage
print("\n[INFO] Generating lineage-specific dotplots...")

lineage_marker_groups = {
    'Alveolar': ['AT1_Canonical', 'AT1_MatrixRemodeling', 'AT2_Canonical', 'AT2_Inflammatory_Repair', 'Epithelial_Cycling'],
    'Basal': ['Basal_Progenitor', 'Basal_Cycling', 'Basal_Inflammatory', 'Basal_EMT_ECM'],
    'Ciliated': ['Ciliated_Mature', 'Ciliogenesis_Deuterosomal', 'Ciliated_Cycling_Immature'],
    'Secretory_SMG': ['Secretory_Club', 'Secretory_Club_AT2_Transitional', 'Goblet_Mucin', 
                      'Goblet_Defense_DUOX2', 'SMG_Serous', 'SMG_Duct_Secretory_Defense'],
    'Special': ['Squamous_Metaplasia', 'Ionocyte_Brush', 'Mesenchymal_Contaminant']
}

for lineage, cell_types in lineage_marker_groups.items():
    print(f"\n  Processing {lineage}...")
    
    # Collect markers for this lineage
    markers_for_lineage = []
    for ct in cell_types:
        if ct in available_markers:
            markers_for_lineage.extend(available_markers[ct])
    
    # Remove duplicates while preserving order
    markers_for_lineage = list(dict.fromkeys(markers_for_lineage))
    
    if len(markers_for_lineage) == 0:
        print(f"    [WARN] No markers available for {lineage}")
        continue
    
    print(f"    Markers: {len(markers_for_lineage)}")
    
    # Filter to L3 labels present in this lineage
    adata_subset = adata[adata.obs['cell_type_L3'].isin(cell_types)].copy()
    
    if adata_subset.n_obs == 0:
        print(f"    [WARN] No cells found for {lineage}")
        continue
    
    print(f"    Cells: {adata_subset.n_obs:,}")
    
    # Create dotplot
    try:
        fig_width = max(12, len(markers_for_lineage) * 0.3)
        fig_height = max(6, len(cell_types) * 0.4)
        
        sc.pl.dotplot(
            adata_subset,
            var_names=markers_for_lineage,
            groupby='cell_type_L3',
            use_raw=use_raw_for_plot,
            standard_scale='var',
            show=False,
            figsize=(fig_width, fig_height),
            dendrogram=True,
            cmap='Reds',
            vmin=-2,
            vmax=2
        )
        
        plt.tight_layout()
        output_file = sc.settings.figdir / f'03_dotplot_{lineage}.{FIGURE_FORMAT}'
        plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
        plt.close()
        print(f"    [OK] Saved: {output_file.name}")
        
    except Exception as e:
        print(f"    [ERROR] Dotplot failed: {e}")
        import traceback
        traceback.print_exc()
    
    del adata_subset
    gc.collect()

# %% [markdown]
# ## Step 6: Generate Comprehensive Marker Heatmap

# %%
print("\n" + "=" * 80)
print("STEP 6: Generating Marker Expression Heatmap")
print("=" * 80)

# Collect all core markers
all_core_markers = []
for markers in available_markers.values():
    all_core_markers.extend(markers)
all_core_markers = list(dict.fromkeys(all_core_markers))  # Remove duplicates, preserve order

print(f"[INFO] Total core markers: {len(all_core_markers)}")

if len(all_core_markers) > 0:
    print("[INFO] Computing mean expression per L3 type...")
    
    # Compute mean expression per L3 type
    l3_types = sorted(adata.obs['cell_type_L3'].unique())
    
    expression_matrix = []
    available_genes = []
    
    for gene in all_core_markers:
        if use_raw_for_plot and gene in adata.raw.var_names:
            gene_expr = adata.raw[:, gene].X
        elif gene in adata.var_names:
            gene_expr = adata[:, gene].X
        else:
            continue
        
        # Convert to dense if sparse
        if sparse.issparse(gene_expr):
            gene_expr = gene_expr.toarray().ravel()
        else:
            gene_expr = np.asarray(gene_expr).ravel()
        
        # Compute mean per L3 type
        mean_expr = []
        for l3_type in l3_types:
            mask = adata.obs['cell_type_L3'] == l3_type
            mean_expr.append(gene_expr[mask].mean())
        
        expression_matrix.append(mean_expr)
        available_genes.append(gene)
    
    if len(expression_matrix) == 0:
        print("[WARN] No expression data available for heatmap")
    else:
        expression_df = pd.DataFrame(
            expression_matrix,
            index=available_genes,
            columns=l3_types
        )
        
        print(f"[OK] Expression matrix: {expression_df.shape[0]} genes × {expression_df.shape[1]} cell types")
        
        # Plot heatmap
        fig_width = max(12, len(l3_types) * 0.4)
        fig_height = max(10, len(available_genes) * 0.15)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        sns.heatmap(
            expression_df,
            cmap='RdYlBu_r',
            center=0,
            robust=True,
            yticklabels=True,
            xticklabels=True,
            cbar_kws={'label': 'Mean Expression (scaled)'},
            linewidths=0.1,
            linecolor='lightgray',
            ax=ax
        )
        
        plt.title('Core Marker Expression Across L3 Cell Types', fontsize=14, pad=20)
        plt.xlabel('L3 Cell Type', fontsize=12)
        plt.ylabel('Marker Gene', fontsize=12)
        plt.xticks(rotation=45, ha='right', fontsize=9)
        plt.yticks(fontsize=8)
        plt.tight_layout()
        
        output_file = sc.settings.figdir / f'04_heatmap_core_markers.{FIGURE_FORMAT}'
        plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
        plt.close()
        print(f"[OK] Saved: {output_file.name}")
        
        del expression_df, expression_matrix
        gc.collect()

# %% [markdown]
# ## Step 7: Generate Marker Expression UMAPs

# %%
print("\n" + "=" * 80)
print("STEP 7: Generating Marker Expression UMAPs")
print("=" * 80)

# Select representative markers for each lineage
representative_markers = {
    'Alveolar': ['AGER', 'HOPX', 'SFTPC', 'NAPSA', 'CHI3L1'],
    'Basal': ['KRT5', 'TP63', 'KRT14', 'KRT17', 'NGFR'],
    'Ciliated': ['FOXJ1', 'DNAH5', 'DEUP1', 'CCNO'],
    'Secretory': ['SCGB1A1', 'MUC5AC', 'DUOX2', 'SPDEF'],
    'SMG': ['LTF', 'PIGR', 'WFDC2'],
    'Special': ['SPRR2A', 'FOXI1', 'CFTR', 'MKI67', 'TOP2A'],
    'ECM': ['COL4A1', 'SPARC', 'FN1', 'VIM']
}

for category, markers in representative_markers.items():
    print(f"\n  Processing {category}...")
    
    available = [m for m in markers if m in gene_universe]
    
    if len(available) == 0:
        print(f"    [WARN] No markers available")
        continue
    
    print(f"    Markers: {', '.join(available)}")
    
    n_markers = len(available)
    n_cols = min(3, n_markers)
    n_rows = int(np.ceil(n_markers / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_rows > 1 else axes
    
    for idx, marker in enumerate(available):
        sc.pl.umap(
            adata,
            color=marker,
            ax=axes[idx],
            show=False,
            use_raw=use_raw_for_plot,
            vmax='p99',
            frameon=False,
            size=UMAP_SIZE * 0.8,
            alpha=UMAP_ALPHA,
            cmap='Reds',
            title=marker
        )
    
    # Hide extra subplots
    for idx in range(len(available), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle(f'{category} Markers', fontsize=14, y=1.02)
    plt.tight_layout()
    
    output_file = sc.settings.figdir / f'05_umap_markers_{category}.{FIGURE_FORMAT}'
    plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"    [OK] Saved: {output_file.name}")

# %% [markdown]
# ## Step 8: Generate Cell Type Proportion Analysis

# %%
print("\n" + "=" * 80)
print("STEP 8: Cell Type Proportion Analysis")
print("=" * 80)

# Overall proportions
print("[INFO] Computing L3 proportions...")
l3_counts = adata.obs['cell_type_L3'].value_counts()
l3_props = (l3_counts / adata.n_obs * 100).sort_values(ascending=True)

# Bar plot
fig, ax = plt.subplots(figsize=(10, 12))
colors = [L3_COLORS.get(label, '#cccccc') for label in l3_props.index]
l3_props.plot(kind='barh', ax=ax, color=colors)
ax.set_xlabel('Percentage of Cells', fontsize=12)
ax.set_ylabel('L3 Cell Type', fontsize=12)
ax.set_title('L3 Cell Type Distribution', fontsize=14)
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
output_file = sc.settings.figdir / f'06_barplot_L3_proportions.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] Saved: {output_file.name}")

# If batch/sample information exists, create stacked bar plot
if 'dataset' in adata.obs.columns:
    print("\n[INFO] Computing per-dataset proportions...")
    
    prop_df = pd.crosstab(
        adata.obs['dataset'],
        adata.obs['cell_type_L3'],
        normalize='index'
    ) * 100
    
    # Transpose for better visualization
    prop_df_t = prop_df.T
    
    # Plot top 15 L3 types
    top_l3 = l3_counts.head(15).index
    prop_df_plot = prop_df_t.loc[top_l3]
    
    fig, ax = plt.subplots(figsize=(14, 8))
    prop_df_plot.plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_xlabel('L3 Cell Type', fontsize=12)
    ax.set_ylabel('Percentage', fontsize=12)
    ax.set_title('L3 Cell Type Distribution Across Datasets (Top 15)', fontsize=14)
    ax.legend(title='Dataset', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    output_file = sc.settings.figdir / f'07_stacked_bar_L3_by_dataset.{FIGURE_FORMAT}'
    plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"[OK] Saved: {output_file.name}")
    
    # Save proportion table
    prop_csv = OUTPUT_DIR / "L3_proportions_by_dataset.csv"
    prop_df.to_csv(prop_csv)
    print(f"[OK] Saved: {prop_csv.name}")

# %% [markdown]
# ## Step 9: Save Annotated Data

# %%
print("\n" + "=" * 80)
print("STEP 9: Saving Annotated Data")
print("=" * 80)

# Save updated h5ad
output_h5ad = OUTPUT_DIR / "epithelial_with_L3_annotations.h5ad"
print(f"[INFO] Writing: {output_h5ad.name}")
adata.write_h5ad(output_h5ad, compression="gzip", compression_opts=9)
file_size_gb = output_h5ad.stat().st_size / 1e9
print(f"[OK] Saved: {output_h5ad.name} ({file_size_gb:.2f} GB)")

# Save L3 annotations as separate CSV
l3_annotations = adata.obs[['subcluster', 'cell_type_L3']].copy()
l3_csv = OUTPUT_DIR / "L3_annotations.csv"
l3_annotations.to_csv(l3_csv)
print(f"[OK] Saved: {l3_csv.name}")

# Save summary statistics
summary_file = OUTPUT_DIR / "L3_annotation_summary.txt"
with open(summary_file, "w") as f:
    f.write("=" * 80 + "\n")
    f.write("Epithelial L3 Annotation Summary\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input: {INPUT_H5AD}\n")
    f.write(f"Output: {OUTPUT_DIR}\n\n")
    
    f.write(f"Dataset Statistics:\n")
    f.write(f"  Total cells: {adata.n_obs:,}\n")
    f.write(f"  Total genes: {adata.n_vars:,}\n")
    f.write(f"  Original subclusters: {adata.obs['subcluster'].nunique()}\n")
    f.write(f"  L3 cell types: {adata.obs['cell_type_L3'].nunique()}\n\n")
    
    f.write("L3 Cell Type Distribution:\n")
    for label, count in l3_counts.items():
        pct = count / adata.n_obs * 100
        f.write(f"  {label}: {count:,} cells ({pct:.2f}%)\n")
    
    f.write("\n" + "=" * 80 + "\n")
    f.write("Key Outputs:\n")
    f.write("=" * 80 + "\n")
    f.write(f"  - Annotated data: epithelial_with_L3_annotations.h5ad\n")
    f.write(f"  - L3 mapping: cluster_to_L3_mapping.csv\n")
    f.write(f"  - L3 annotations: L3_annotations.csv\n")
    f.write(f"  - Figures: figures/ directory\n")
    f.write("\nFigures generated:\n")
    f.write(f"  - 01_umap_L3_overview.{FIGURE_FORMAT}\n")
    f.write(f"  - 02_umap_L3_by_lineage.{FIGURE_FORMAT}\n")
    f.write(f"  - 03_dotplot_[lineage].{FIGURE_FORMAT}\n")
    f.write(f"  - 04_heatmap_core_markers.{FIGURE_FORMAT}\n")
    f.write(f"  - 05_umap_markers_[category].{FIGURE_FORMAT}\n")
    f.write(f"  - 06_barplot_L3_proportions.{FIGURE_FORMAT}\n")
    if 'dataset' in adata.obs.columns:
        f.write(f"  - 07_stacked_bar_L3_by_dataset.{FIGURE_FORMAT}\n")

print(f"[OK] Saved: {summary_file.name}")

# %% [markdown]
# ## Summary

# %%
elapsed_time = time.time() - PIPELINE_START
print("\n" + "=" * 80)
print("VISUALIZATION PIPELINE COMPLETE")
print("=" * 80)
print(f"Output directory: {OUTPUT_DIR}")
print(f"Annotated data: epithelial_with_L3_annotations.h5ad")
print(f"Total L3 types: {adata.obs['cell_type_L3'].nunique()}")
print(f"Elapsed time: {elapsed_time / 60:.1f} min")
print(f"\nKey outputs:")
print(f"  - L3 overview UMAP")
print(f"  - Lineage-specific UMAPs")
print(f"  - Marker dotplots per lineage")
print(f"  - Core marker heatmap")
print(f"  - Marker expression UMAPs")
print(f"  - Cell type proportion plots")
print("\nAll figures saved in: figures/")
print("=" * 80)

# %%
