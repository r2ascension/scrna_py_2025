#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Subset Cells by scArches Mapped Cell Types (Production Version)

Purpose:
- Read scArches mapped query data
- Subset to major cell types
- Save each cell type as separate h5ad file
- Generate cell type distribution report

Input:
- query_mapped_to_reference.h5ad (from step2_scarches_mapping.py)

Output:
- Individual h5ad files for each major cell type
- Cell type distribution summary

Author: r2end
Date: 2025-01-27
Version: 1.0
"""

# ==============================================================================
# Step 0: Imports and Configuration
# ==============================================================================

import sys
import os
from pathlib import Path
import warnings
import json
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc

warnings.filterwarnings('ignore')

print("=" * 70)
print("Cell Type Subsetting Pipeline")
print("=" * 70)
print(f"\nscanpy: {sc.__version__}")
print(f"Python: {sys.version}")

# ==============================================================================
# CONFIGURATION SECTION (⚠️ MODIFY THIS)
# ==============================================================================

# ===== Input/Output Paths =====
INPUT_H5AD = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/query_mapped_to_reference.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets"

# ===== Cell Type Key =====
CELLTYPE_KEY = "cell_type_final"  # 验证一下实际列名

# ===== Subsetting Strategy =====
SUBSET_ALL = True       # 每个类型单独一个文件（推荐）
SUBSET_MAJOR = False

# ===== Quality Control =====
REMOVE_UNKNOWN = True              # 移除86个Unknown细胞
CONFIDENCE_THRESHOLD = 0.5         # 已经很高了，保持默认

# ===== Minimum Cells =====
MIN_CELLS_PER_TYPE = 50

# ===== Quality Control Filters =====
REMOVE_UNKNOWN = True              # Remove "Unknown" cells before subsetting
CONFIDENCE_THRESHOLD = 0.5         # Remove cells below this confidence (if available)

# ===== Compression =====
COMPRESSION = 'gzip'
COMPRESSION_OPTS = 9

# ===== Visualization =====
DPI = 300
FIGURE_FORMAT = 'pdf'

# ===== Reproducibility =====
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)

# Validate input path
if not Path(INPUT_H5AD).exists():
    print(f"\n❌ ERROR: Input file not found at {INPUT_H5AD}")
    print("   Please update INPUT_H5AD in the configuration section")
    sys.exit(1)

# Create output directory
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)

print(f"\nOutput directory: {output_dir}")

# Save configuration
config = {
    'pipeline': 'celltype_subsetting',
    'version': '1.0',
    'input_h5ad': str(INPUT_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'celltype_key': CELLTYPE_KEY,
    'subset_strategy': 'major' if SUBSET_MAJOR else 'all',
    'min_cells_per_type': MIN_CELLS_PER_TYPE,
    'remove_unknown': REMOVE_UNKNOWN,
    'confidence_threshold': CONFIDENCE_THRESHOLD,
    'random_seed': RANDOM_SEED,
    'timestamp': datetime.now().isoformat()
}

with open(output_dir / 'subset_config.json', 'w') as f:
    json.dump(config, f, indent=2)

print("✓ Configuration saved")

# ==============================================================================
# Step 1: Load Mapped Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1: Loading Mapped Data")
print("=" * 70)

print(f"\nLoading from: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

print(f"\nData loaded:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Obs columns: {list(adata.obs.columns)}")

# Check if cell type key exists
if CELLTYPE_KEY not in adata.obs.columns:
    # Try alternative keys
    alternative_keys = [
        'cell_type_mapped',
        'cell_type',
        'celltype_mapped',
        'celltype',
        'predicted_labels',
        'cell_type_pred'
    ]
    
    found = False
    for alt_key in alternative_keys:
        if alt_key in adata.obs.columns:
            print(f"\n⚠️ '{CELLTYPE_KEY}' not found, using '{alt_key}' instead")
            CELLTYPE_KEY = alt_key
            found = True
            break
    
    if not found:
        print(f"\n❌ ERROR: No cell type annotation found in obs columns")
        print(f"   Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
else:
    print(f"\n✓ Found cell type key: {CELLTYPE_KEY}")

# Display cell type distribution
print(f"\nCell Type Distribution:")
celltype_counts = adata.obs[CELLTYPE_KEY].value_counts()
print(celltype_counts.head(30))

# ==============================================================================
# Step 2: Quality Control Filtering (Optional)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Quality Control Filtering")
print("=" * 70)

n_cells_before = adata.n_obs

# 1. Remove Unknown cells
if REMOVE_UNKNOWN:
    unknown_labels = ['Unknown', 'unknown', 'Uncertain', 'Unassigned', 'NA']
    unknown_mask = adata.obs[CELLTYPE_KEY].isin(unknown_labels)
    n_unknown = unknown_mask.sum()
    
    if n_unknown > 0:
        print(f"\n1. Removing Unknown cells:")
        print(f"   Unknown cells: {n_unknown:,} ({n_unknown/n_cells_before*100:.1f}%)")
        adata = adata[~unknown_mask].copy()
        print(f"   Remaining: {adata.n_obs:,} cells")
    else:
        print(f"\n1. No Unknown cells found")

# 2. Remove low-confidence predictions (if confidence score exists)
confidence_keys = ['mapping_confidence', 'scanvi_prediction_confidence', 
                   'prediction_confidence', 'confidence']
confidence_key = None

for ck in confidence_keys:
    if ck in adata.obs.columns:
        confidence_key = ck
        break

if confidence_key is not None and CONFIDENCE_THRESHOLD > 0:
    low_conf_mask = adata.obs[confidence_key] < CONFIDENCE_THRESHOLD
    n_low_conf = low_conf_mask.sum()
    
    if n_low_conf > 0:
        print(f"\n2. Removing low-confidence cells (< {CONFIDENCE_THRESHOLD}):")
        print(f"   Low confidence: {n_low_conf:,} ({n_low_conf/adata.n_obs*100:.1f}%)")
        adata = adata[~low_conf_mask].copy()
        print(f"   Remaining: {adata.n_obs:,} cells")
    else:
        print(f"\n2. All cells above confidence threshold")
else:
    print(f"\n2. No confidence filtering (key not found or threshold=0)")

n_cells_after = adata.n_obs
print(f"\nTotal filtered: {n_cells_before:,} → {n_cells_after:,} cells "
      f"({(1 - n_cells_after/n_cells_before)*100:.1f}% removed)")

# Update cell type distribution after filtering
celltype_counts_filtered = adata.obs[CELLTYPE_KEY].value_counts()
print(f"\nCell Type Distribution After Filtering:")
print(celltype_counts_filtered.head(30))

print("\n✓ Quality control completed")

# ==============================================================================
# Step 3: Define Cell Type Groups
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Defining Cell Type Groups")
print("=" * 70)

# Get all unique cell types
all_celltypes = adata.obs[CELLTYPE_KEY].unique().tolist()
print(f"\nTotal unique cell types: {len(all_celltypes)}")

if SUBSET_MAJOR:
    print("\nUsing MAJOR cell type grouping strategy")
    print("\nDefined major groups:")
    
    # Build reverse mapping: celltype -> group
    celltype_to_group = {}
    for group_name, keywords in MAJOR_CELLTYPE_GROUPS.items():
        print(f"\n  {group_name}:")
        matched_types = []
        
        for celltype in all_celltypes:
            # Case-insensitive matching
            celltype_lower = celltype.lower()
            for keyword in keywords:
                keyword_lower = keyword.lower()
                if keyword_lower in celltype_lower:
                    celltype_to_group[celltype] = group_name
                    matched_types.append(celltype)
                    break
        
        print(f"    Matched: {len(matched_types)} types")
        if matched_types:
            print(f"    Examples: {matched_types[:5]}")
    
    # Identify ungrouped cell types
    ungrouped = [ct for ct in all_celltypes if ct not in celltype_to_group]
    if ungrouped:
        print(f"\n⚠️ Ungrouped cell types ({len(ungrouped)}):")
        for ct in ungrouped[:10]:
            count = (adata.obs[CELLTYPE_KEY] == ct).sum()
            print(f"    - {ct}: {count} cells")
        
        if len(ungrouped) > 10:
            print(f"    ... and {len(ungrouped) - 10} more")
        
        # Ask user what to do with ungrouped
        print(f"\n  These will be saved as 'Other' group")
        celltype_to_group.update({ct: 'Other' for ct in ungrouped})
    
    # Add group annotation to adata
    adata.obs['celltype_group'] = adata.obs[CELLTYPE_KEY].map(celltype_to_group)
    
    # Get list of groups to subset
    groups_to_subset = list(MAJOR_CELLTYPE_GROUPS.keys())
    if ungrouped:
        groups_to_subset.append('Other')
    
    print(f"\nWill create {len(groups_to_subset)} subset files:")
    for group in groups_to_subset:
        n_cells = (adata.obs['celltype_group'] == group).sum()
        print(f"  - {group}: {n_cells:,} cells")

else:
    print("\nUsing ALL cell types strategy (one file per type)")
    
    # Filter by minimum cells
    celltypes_to_subset = []
    for celltype in all_celltypes:
        n_cells = (adata.obs[CELLTYPE_KEY] == celltype).sum()
        if n_cells >= MIN_CELLS_PER_TYPE:
            celltypes_to_subset.append(celltype)
        else:
            print(f"  Skipping {celltype}: only {n_cells} cells (< {MIN_CELLS_PER_TYPE})")
    
    print(f"\nWill create {len(celltypes_to_subset)} subset files:")
    for celltype in celltypes_to_subset[:20]:
        n_cells = (adata.obs[CELLTYPE_KEY] == celltype).sum()
        print(f"  - {celltype}: {n_cells:,} cells")
    
    if len(celltypes_to_subset) > 20:
        print(f"  ... and {len(celltypes_to_subset) - 20} more")

print("\n✓ Cell type groups defined")

# ==============================================================================
# Step 4: Subset and Save
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Subsetting and Saving")
print("=" * 70)

subset_summary = []

if SUBSET_MAJOR:
    # Subset by major groups
    for group_name in groups_to_subset:
        print(f"\n-" * 70)
        print(f"Processing: {group_name}")
        print(f"-" * 70)
        
        # Subset data
        mask = adata.obs['celltype_group'] == group_name
        adata_subset = adata[mask].copy()
        
        n_cells = adata_subset.n_obs
        print(f"  Cells: {n_cells:,}")
        
        if n_cells < MIN_CELLS_PER_TYPE:
            print(f"  ⚠️ Skipping (< {MIN_CELLS_PER_TYPE} cells)")
            continue
        
        # Cell type composition within group
        celltype_dist = adata_subset.obs[CELLTYPE_KEY].value_counts()
        print(f"  Unique cell types: {len(celltype_dist)}")
        print(f"  Top types:")
        for ct, count in celltype_dist.head(5).items():
            pct = count / n_cells * 100
            print(f"    - {ct}: {count:,} ({pct:.1f}%)")
        
        # Save file
        output_file = output_dir / f"{group_name.lower()}_cells.h5ad"
        print(f"\n  Saving to: {output_file.name}")
        
        try:
            adata_subset.write_h5ad(
                output_file,
                compression=COMPRESSION,
                compression_opts=COMPRESSION_OPTS
            )
            file_size = output_file.stat().st_size / (1024**2)  # MB
            print(f"  ✓ Saved ({file_size:.1f} MB)")
            
            subset_summary.append({
                'group': group_name,
                'cells': n_cells,
                'unique_types': len(celltype_dist),
                'file': str(output_file.name),
                'size_mb': round(file_size, 2)
            })
            
        except Exception as e:
            print(f"  ❌ ERROR saving {group_name}: {e}")

else:
    # Subset by individual cell types
    for i, celltype in enumerate(celltypes_to_subset, 1):
        print(f"\n-" * 70)
        print(f"Processing [{i}/{len(celltypes_to_subset)}]: {celltype}")
        print(f"-" * 70)
        
        # Subset data
        mask = adata.obs[CELLTYPE_KEY] == celltype
        adata_subset = adata[mask].copy()
        
        n_cells = adata_subset.n_obs
        print(f"  Cells: {n_cells:,}")
        
        # Save file (sanitize filename)
        safe_name = celltype.replace(' ', '_').replace('/', '_').replace('+', 'pos')
        output_file = output_dir / f"{safe_name.lower()}_cells.h5ad"
        print(f"  Saving to: {output_file.name}")
        
        try:
            adata_subset.write_h5ad(
                output_file,
                compression=COMPRESSION,
                compression_opts=COMPRESSION_OPTS
            )
            file_size = output_file.stat().st_size / (1024**2)  # MB
            print(f"  ✓ Saved ({file_size:.1f} MB)")
            
            subset_summary.append({
                'celltype': celltype,
                'cells': n_cells,
                'file': str(output_file.name),
                'size_mb': round(file_size, 2)
            })
            
        except Exception as e:
            print(f"  ❌ ERROR saving {celltype}: {e}")

print("\n✓ Subsetting completed")

# ==============================================================================
# Step 5: Generate Visualization
# ==============================================================================

print("\n" + "=" * 70)
print("Step 5: Generating Visualization")
print("=" * 70)

# Create summary dataframe
summary_df = pd.DataFrame(subset_summary)

if len(summary_df) > 0:
    # Sort by cell count
    summary_df = summary_df.sort_values('cells', ascending=False)
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Cell counts
    if SUBSET_MAJOR:
        x_col = 'group'
        title1 = 'Cell Counts by Major Group'
    else:
        x_col = 'celltype'
        title1 = 'Cell Counts by Cell Type'
    
    axes[0].bar(range(len(summary_df)), summary_df['cells'])
    axes[0].set_xticks(range(len(summary_df)))
    axes[0].set_xticklabels(summary_df[x_col], rotation=45, ha='right')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title(title1)
    axes[0].grid(axis='y', alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(summary_df['cells']):
        axes[0].text(i, v, f'{v:,}', ha='center', va='bottom', fontsize=8)
    
    # Plot 2: File sizes
    axes[1].bar(range(len(summary_df)), summary_df['size_mb'], color='orange')
    axes[1].set_xticks(range(len(summary_df)))
    axes[1].set_xticklabels(summary_df[x_col], rotation=45, ha='right')
    axes[1].set_ylabel('File Size (MB)')
    axes[1].set_title('Output File Sizes')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    plot_file = output_dir / "figures" / f"subset_summary.{FIGURE_FORMAT}"
    plt.savefig(plot_file, dpi=DPI, bbox_inches='tight')
    plt.savefig(output_dir / "figures" / "subset_summary.png", 
                dpi=150, bbox_inches='tight')
    print(f"✓ Saved: subset_summary.{FIGURE_FORMAT}")
    plt.close()

print("\n✓ Visualization completed")

# ==============================================================================
# Step 6: Export Summary Report
# ==============================================================================

print("\n" + "=" * 70)
print("Step 6: Generating Summary Report")
print("=" * 70)

# Save summary table
summary_csv = output_dir / "subset_summary.csv"
summary_df.to_csv(summary_csv, index=False)
print(f"✓ Saved CSV: {summary_csv}")

# Generate text report
report = []
report.append("=" * 70)
report.append("Cell Type Subsetting Summary")
report.append("=" * 70)
report.append(f"\nTimestamp: {datetime.now().isoformat()}")
report.append(f"\nInput: {INPUT_H5AD}")
report.append(f"Output Directory: {output_dir}")

report.append(f"\n\n[ Input Data ]")
report.append(f"  Total cells before filtering: {n_cells_before:,}")
report.append(f"  Total cells after filtering: {n_cells_after:,}")
report.append(f"  Filtered out: {n_cells_before - n_cells_after:,} "
              f"({(1 - n_cells_after/n_cells_before)*100:.1f}%)")

report.append(f"\n[ Subsetting Strategy ]")
if SUBSET_MAJOR:
    report.append(f"  Strategy: Major cell type groups")
    report.append(f"  Number of groups: {len(groups_to_subset)}")
else:
    report.append(f"  Strategy: Individual cell types")
    report.append(f"  Minimum cells per type: {MIN_CELLS_PER_TYPE}")

report.append(f"\n[ Output Files ]")
report.append(f"  Number of files created: {len(summary_df)}")
report.append(f"  Total size: {summary_df['size_mb'].sum():.1f} MB")

report.append(f"\n[ Cell Distribution ]")
for _, row in summary_df.iterrows():
    if SUBSET_MAJOR:
        name = row['group']
        extra = f" ({row['unique_types']} types)"
    else:
        name = row['celltype']
        extra = ""
    
    report.append(f"  {name}{extra}:")
    report.append(f"    Cells: {row['cells']:,}")
    report.append(f"    File: {row['file']} ({row['size_mb']} MB)")

report.append("\n" + "=" * 70)
report.append("Subsetting completed successfully")
report.append("=" * 70)

report_text = '\n'.join(report)
print("\n" + report_text)

report_file = output_dir / "subset_summary.txt"
with open(report_file, 'w') as f:
    f.write(report_text)
print(f"\n✓ Saved text report: {report_file}")

# ==============================================================================
# PIPELINE COMPLETE
# ==============================================================================

print("\n" + "=" * 70)
print("✅ CELL TYPE SUBSETTING COMPLETE")
print("=" * 70)

print(f"\n📊 Summary:")
print(f"  Input cells: {n_cells_after:,}")
print(f"  Output files: {len(summary_df)}")
print(f"  Output directory: {output_dir}")

print(f"\n📈 Next Steps:")
print(f"  - Review cell type distribution in subset_summary.csv")
print(f"  - Perform cell-type-specific analyses on each subset")
print(f"  - Consider further QC for individual cell types")
print(f"  - Run downstream analyses (DEG, trajectory, etc.)")

print("\n" + "=" * 70)
print("Thank you for using the cell type subsetting pipeline!")
print("=" * 70)
