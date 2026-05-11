#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Query Results Visualization and Analysis
=========================================
Analyzes query_mapped_L2.h5ad results without needing merge

Usage:
    python visualize_query_results.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc

print("=" * 70)
print("Query Results Visualization")
print("=" * 70)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

QUERY_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/query_mapped_L2.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/analysis"

CONFIDENCE_THRESHOLD = 0.5
DPI = 300

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

sc.settings.set_figure_params(dpi=DPI, facecolor='white', format='pdf')

print(f"\nConfiguration:")
print(f"  Query: {QUERY_H5AD}")
print(f"  Output: {OUTPUT_DIR}")

# ==============================================================================
# Load Data
# ==============================================================================

print("\n" + "=" * 70)
print("Load and Inspect Query Results")
print("=" * 70)

print(f"\nLoading query...")
adata = sc.read_h5ad(QUERY_H5AD)

print(f"\nDataset info:")
print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  obs columns: {len(adata.obs.columns)}")
print(f"  obsm keys: {list(adata.obsm.keys())}")
print(f"  layers: {list(adata.layers.keys())}")

# Check critical columns
critical_cols = ['Cell_Type_L2_pred', 'Cell_Type_L2_final', 'mapping_confidence']
missing = [c for c in critical_cols if c not in adata.obs.columns]
if missing:
    print(f"\nWARNING: Missing columns: {missing}")
else:
    print(f"\nOK All critical columns present")

# ==============================================================================
# Summary Statistics
# ==============================================================================

print("\n" + "=" * 70)
print("Summary Statistics")
print("=" * 70)

# L2 distribution
print(f"\nL2 Cell Type Distribution:")
l2_counts = adata.obs['Cell_Type_L2_final'].value_counts()
for celltype, count in l2_counts.items():
    pct = count / adata.n_obs * 100
    print(f"  {celltype}: {count:,} ({pct:.1f}%)")

# Confidence
conf = adata.obs['mapping_confidence']
print(f"\nMapping Confidence:")
print(f"  Mean: {conf.mean():.3f}")
print(f"  Median: {conf.median():.3f}")
print(f"  Min: {conf.min():.3f}")
print(f"  Max: {conf.max():.3f}")

high_conf = (conf >= CONFIDENCE_THRESHOLD).sum()
print(f"\nHigh confidence (>={CONFIDENCE_THRESHOLD}): {high_conf:,} ({high_conf/adata.n_obs*100:.1f}%)")

# ==============================================================================
# Visualization 1: UMAP Overview
# ==============================================================================

print("\n" + "=" * 70)
print("Generating Visualizations")
print("=" * 70)

print("\n1. UMAP overview (4 panels)...")

fig, axes = plt.subplots(2, 2, figsize=(16, 16))

# Panel 1: L2 predictions (all)
sc.pl.umap(adata, color='Cell_Type_L2_pred', ax=axes[0, 0], show=False,
           title='L2 Predictions (All Cells)', legend_loc='right margin', s=10)

# Panel 2: L2 final (high conf only)
sc.pl.umap(adata, color='Cell_Type_L2_final', ax=axes[0, 1], show=False,
           title=f'L2 Final (Confidence >= {CONFIDENCE_THRESHOLD})', 
           legend_loc='right margin', s=10)

# Panel 3: Confidence
sc.pl.umap(adata, color='mapping_confidence', ax=axes[1, 0], show=False,
           title='Mapping Confidence', cmap='viridis', vmin=0, vmax=1, s=10)

# Panel 4: Confidence histogram
conf_vals = adata.obs['mapping_confidence'].values
axes[1, 1].hist(conf_vals, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
axes[1, 1].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--', 
                   linewidth=2, label=f'Threshold={CONFIDENCE_THRESHOLD}')
axes[1, 1].set_xlabel('Mapping Confidence', fontsize=12)
axes[1, 1].set_ylabel('Number of Cells', fontsize=12)
axes[1, 1].set_title('Confidence Distribution', fontsize=14, fontweight='bold')
axes[1, 1].legend()
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
fig_path = output_dir / "query_umap_overview.pdf"
plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"   Saved: {fig_path.name}")

# ==============================================================================
# Visualization 2: Cell Type Composition
# ==============================================================================

print("\n2. Cell type composition...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Bar plot
l2_counts_sorted = l2_counts.sort_values(ascending=True)
colors = plt.cm.Set3(np.linspace(0, 1, len(l2_counts_sorted)))

axes[0].barh(range(len(l2_counts_sorted)), l2_counts_sorted.values, color=colors)
axes[0].set_yticks(range(len(l2_counts_sorted)))
axes[0].set_yticklabels(l2_counts_sorted.index)
axes[0].set_xlabel('Number of Cells', fontsize=12)
axes[0].set_title('L2 Cell Type Counts', fontsize=14, fontweight='bold')
axes[0].grid(axis='x', alpha=0.3)

# Add count labels
for i, (celltype, count) in enumerate(l2_counts_sorted.items()):
    axes[0].text(count + max(l2_counts_sorted)*0.01, i, f'{count:,}', 
                va='center', fontsize=10)

# Pie chart
axes[1].pie(l2_counts.values, labels=l2_counts.index, autopct='%1.1f%%',
            colors=colors, startangle=90)
axes[1].set_title('L2 Cell Type Proportions', fontsize=14, fontweight='bold')

plt.tight_layout()
fig_path = output_dir / "celltype_composition.pdf"
plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"   Saved: {fig_path.name}")

# ==============================================================================
# Visualization 3: Confidence by Cell Type
# ==============================================================================

print("\n3. Confidence by cell type...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Violin plot
df = adata.obs[['Cell_Type_L2_final', 'mapping_confidence']].copy()
df = df[df['Cell_Type_L2_final'] != 'Unknown']  # Exclude Unknown

celltypes_sorted = df['Cell_Type_L2_final'].value_counts().index.tolist()

sns.violinplot(data=df, y='Cell_Type_L2_final', x='mapping_confidence', 
               order=celltypes_sorted, ax=axes[0], orient='h')
axes[0].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--', linewidth=2)
axes[0].set_xlabel('Mapping Confidence', fontsize=12)
axes[0].set_ylabel('L2 Cell Type', fontsize=12)
axes[0].set_title('Confidence Distribution by Cell Type', fontsize=14, fontweight='bold')
axes[0].grid(axis='x', alpha=0.3)

# Box plot with stats
conf_by_type = df.groupby('Cell_Type_L2_final')['mapping_confidence'].describe()
conf_by_type = conf_by_type.sort_values('mean', ascending=True)

axes[1].barh(range(len(conf_by_type)), conf_by_type['mean'], 
             xerr=conf_by_type['std'], capsize=5, color='skyblue', 
             edgecolor='black', alpha=0.7)
axes[1].set_yticks(range(len(conf_by_type)))
axes[1].set_yticklabels(conf_by_type.index)
axes[1].set_xlabel('Mean Confidence ± Std', fontsize=12)
axes[1].set_title('Mean Confidence by Cell Type', fontsize=14, fontweight='bold')
axes[1].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--', linewidth=2)
axes[1].grid(axis='x', alpha=0.3)

plt.tight_layout()
fig_path = output_dir / "confidence_by_celltype.pdf"
plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"   Saved: {fig_path.name}")

# ==============================================================================
# Visualization 4: Marker Gene Expression
# ==============================================================================

print("\n4. B cell marker expression...")

# Define key B cell markers
bcell_markers = {
    'Pan-B': ['CD19', 'MS4A1', 'CD79A'],
    'Naive': ['TCL1A', 'FCER2', 'IGHD'],
    'Memory': ['CD27', 'TNFRSF13B'],
    'GC': ['AICDA', 'BCL6', 'MME'],
    'Plasma': ['SDC1', 'XBP1', 'JCHAIN', 'MZB1'],
    'Atypical': ['ITGAX', 'TBX21', 'FCRL5']
}

# Check which markers are in data
available_markers = []
for category, markers in bcell_markers.items():
    for marker in markers:
        if marker in adata.var_names:
            available_markers.append(marker)

if len(available_markers) > 0:
    print(f"   Found {len(available_markers)} markers in data")
    
    # UMAP with markers (max 12)
    n_markers = min(12, len(available_markers))
    markers_to_plot = available_markers[:n_markers]
    
    n_cols = 4
    n_rows = (n_markers + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes
    
    for i, marker in enumerate(markers_to_plot):
        sc.pl.umap(adata, color=marker, ax=axes[i], show=False, 
                   cmap='Reds', s=10, title=marker)
    
    # Hide extra subplots
    for i in range(n_markers, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    fig_path = output_dir / "marker_expression_umap.pdf"
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"   Saved: {fig_path.name}")
else:
    print("   WARNING: No standard B cell markers found in data")

# ==============================================================================
# Summary Report
# ==============================================================================

print("\n" + "=" * 70)
print("Generate Summary Report")
print("=" * 70)

report_lines = []
report_lines.append("=" * 70)
report_lines.append("Query Mapping Results - Summary Report")
report_lines.append("=" * 70)

report_lines.append(f"\nDataset: {QUERY_H5AD}")
report_lines.append(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")

report_lines.append(f"\n[ Dataset Size ]")
report_lines.append(f"  Cells: {adata.n_obs:,}")
report_lines.append(f"  Genes: {adata.n_vars:,}")

report_lines.append(f"\n[ L2 Cell Type Distribution ]")
for celltype, count in l2_counts.items():
    pct = count / adata.n_obs * 100
    report_lines.append(f"  {celltype}: {count:,} ({pct:.1f}%)")

report_lines.append(f"\n[ Mapping Quality ]")
report_lines.append(f"  Confidence mean: {conf.mean():.3f}")
report_lines.append(f"  Confidence median: {conf.median():.3f}")
report_lines.append(f"  High confidence (>={CONFIDENCE_THRESHOLD}): {high_conf:,} ({high_conf/adata.n_obs*100:.1f}%)")

report_lines.append(f"\n[ Confidence by Cell Type ]")
for celltype in conf_by_type.index:
    mean_conf = conf_by_type.loc[celltype, 'mean']
    std_conf = conf_by_type.loc[celltype, 'std']
    report_lines.append(f"  {celltype}: {mean_conf:.3f} ± {std_conf:.3f}")

report_lines.append(f"\n[ Generated Files ]")
report_lines.append(f"  - query_umap_overview.pdf")
report_lines.append(f"  - celltype_composition.pdf")
report_lines.append(f"  - confidence_by_celltype.pdf")
if len(available_markers) > 0:
    report_lines.append(f"  - marker_expression_umap.pdf")

report_lines.append(f"\n" + "=" * 70)
report_lines.append("Analysis Complete")
report_lines.append("=" * 70)

report_text = '\n'.join(report_lines)
print(report_text)

report_file = output_dir / "analysis_summary.txt"
with open(report_file, 'w') as f:
    f.write(report_text)

print(f"\nReport saved: {report_file}")

print("\n" + "=" * 70)
print("SUCCESS - All Visualizations Generated")
print("=" * 70)
print(f"\nOutput directory: {output_dir}")
print(f"Generated files:")
print(f"  - query_umap_overview.pdf")
print(f"  - celltype_composition.pdf")
print(f"  - confidence_by_celltype.pdf")
if len(available_markers) > 0:
    print(f"  - marker_expression_umap.pdf")
print(f"  - analysis_summary.txt")
