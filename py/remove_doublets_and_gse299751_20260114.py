"""
Data Cleaning: Remove Doublets and Problematic Dataset
======================================================

Remove:
1. High-confidence doublets (28,406 cells, 4.38%)
2. GSE299751 dataset (problematic geometric compatibility)

Author: r2end
Date: 2025-01-14
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# ==============================================================================
# Configuration
# ==============================================================================

INPUT_H5AD = "/home/h2048/data/py/0112/celltypist_stromal/celltypist_doublet_analysis/adata_with_doublet_predictions.h5ad"
OUTPUT_DIR = Path("/home/h2048/data/py/0112/celltypist_stromal/cleaned_data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Datasets to remove (in addition to doublets)
DATASETS_TO_REMOVE = ['Matthew_C_Altman_2025']

# ==============================================================================
# Load Data
# ==============================================================================

print("\n" + "=" * 70)
print("Loading Data with Doublet Predictions")
print("=" * 70)

adata = sc.read_h5ad(INPUT_H5AD)
print(f"\nOriginal data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# ==============================================================================
# Generate Filtering Report
# ==============================================================================

print("\n" + "=" * 70)
print("Pre-Filtering Analysis")
print("=" * 70)

# Count doublets
n_doublets = adata.obs['doublet_high_confidence'].sum()
print(f"\n1. Doublets to remove:")
print(f"   High-confidence doublets: {n_doublets:,} cells ({n_doublets/adata.n_obs*100:.2f}%)")

# Show doublet type distribution
if n_doublets > 0:
    doublet_types = adata[adata.obs['doublet_high_confidence']].obs['doublet_type'].value_counts()
    print(f"\n   Top doublet types:")
    for dtype, count in doublet_types.head(5).items():
        if dtype != 'Singlet':
            print(f"     - {dtype}: {count:,} cells")

# Count cells in problematic datasets
print(f"\n2. Problematic datasets to remove:")
for ds in DATASETS_TO_REMOVE:
    if ds in adata.obs['dataset'].values:
        n_cells = (adata.obs['dataset'] == ds).sum()
        print(f"   - {ds}: {n_cells:,} cells ({n_cells/adata.n_obs*100:.2f}%)")
    else:
        print(f"   - {ds}: Not found in data")

# Check overlap (cells that are both doublets AND in problematic datasets)
mask_doublet = adata.obs['doublet_high_confidence']
mask_dataset = adata.obs['dataset'].isin(DATASETS_TO_REMOVE)
overlap = (mask_doublet & mask_dataset).sum()

print(f"\n3. Overlap:")
print(f"   Cells that are BOTH doublets AND in problematic datasets: {overlap:,}")
print(f"   (These will only be counted once in removal)")

# Total cells to remove
total_to_remove = (mask_doublet | mask_dataset).sum()
print(f"\n4. Total cells to remove: {total_to_remove:,} ({total_to_remove/adata.n_obs*100:.2f}%)")

# ==============================================================================
# Apply Filters
# ==============================================================================

print("\n" + "=" * 70)
print("Applying Filters")
print("=" * 70)

# Create combined filter mask
filter_mask = (
    (adata.obs['doublet_high_confidence'] == False) &  # Not a doublet
    (~adata.obs['dataset'].isin(DATASETS_TO_REMOVE))   # Not in problematic datasets
)

# Apply filter
adata_clean = adata[filter_mask].copy()

print(f"\n✓ Filtering complete:")
print(f"  Original: {adata.n_obs:,} cells")
print(f"  Cleaned: {adata_clean.n_obs:,} cells")
print(f"  Removed: {adata.n_obs - adata_clean.n_obs:,} cells ({(1-adata_clean.n_obs/adata.n_obs)*100:.2f}%)")

# ==============================================================================
# Post-Filtering Analysis
# ==============================================================================

print("\n" + "=" * 70)
print("Post-Filtering Analysis")
print("=" * 70)

# Dataset distribution
print("\n1. Remaining datasets:")
remaining_datasets = adata_clean.obs['dataset'].value_counts()
print(f"   Total datasets: {len(remaining_datasets)}")
print(f"\n   Top 10 datasets by cell count:")
for ds, count in remaining_datasets.head(10).items():
    pct = count / adata_clean.n_obs * 100
    print(f"     {ds}: {count:,} cells ({pct:.2f}%)")

# Cell type distribution (if available)
if 'cell_type_mapped' in adata_clean.obs.columns:
    print("\n2. Cell type distribution:")
    cell_types = adata_clean.obs['cell_type_mapped'].value_counts()
    for ct, count in cell_types.head(10).items():
        pct = count / adata_clean.n_obs * 100
        print(f"     {ct}: {count:,} cells ({pct:.2f}%)")

# Doublet status (should all be False now)
remaining_doublets = adata_clean.obs['doublet_high_confidence'].sum()
if remaining_doublets > 0:
    print(f"\n⚠️  Warning: {remaining_doublets} doublets remain (unexpected!)")
else:
    print(f"\n✓ All doublets removed")

# ==============================================================================
# Save Cleaned Data
# ==============================================================================

print("\n" + "=" * 70)
print("Saving Cleaned Data")
print("=" * 70)

output_h5ad = OUTPUT_DIR / "adata_cleaned_no_doublets_no_GSE299751.h5ad"
print(f"\nSaving to: {output_h5ad}")

adata_clean.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)

file_size = output_h5ad.stat().st_size / (1024**3)
print(f"✓ Saved ({file_size:.2f} GB)")

# Save filtering summary
summary = {
    'original_cells': int(adata.n_obs),
    'cleaned_cells': int(adata_clean.n_obs),
    'removed_cells': int(adata.n_obs - adata_clean.n_obs),
    'removal_rate': float((adata.n_obs - adata_clean.n_obs) / adata.n_obs),
    'removal_breakdown': {
        'doublets': int(n_doublets),
        'problematic_datasets': {ds: int((adata.obs['dataset'] == ds).sum()) for ds in DATASETS_TO_REMOVE if ds in adata.obs['dataset'].values},
        'overlap': int(overlap)
    },
    'remaining_datasets': int(len(remaining_datasets)),
    'datasets_removed': DATASETS_TO_REMOVE
}

import json
summary_file = OUTPUT_DIR / "cleaning_summary.json"
with open(summary_file, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"✓ Saved summary: {summary_file}")

# Save dataset statistics
dataset_stats = pd.DataFrame({
    'n_cells': remaining_datasets,
    'percentage': (remaining_datasets / adata_clean.n_obs * 100).round(2)
})
dataset_stats = dataset_stats.sort_values('n_cells', ascending=False)
dataset_stats.to_csv(OUTPUT_DIR / "remaining_datasets.csv")

print(f"✓ Saved dataset stats: {OUTPUT_DIR / 'remaining_datasets.csv'}")

# ==============================================================================
# Generate Comparison Plots (Optional)
# ==============================================================================

print("\n" + "=" * 70)
print("Generating Comparison Plots (Optional)")
print("=" * 70)

try:
    import matplotlib
    matplotlib.use('Agg')  # For headless environments
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Plot 1: Removal breakdown
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    removal_counts = {
        'Doublets': n_doublets,
        'GSE299751': (adata.obs['dataset'] == 'GSE299751').sum() if 'GSE299751' in adata.obs['dataset'].values else 0,
        'Overlap': overlap,
        'Retained': adata_clean.n_obs
    }
    
    colors = ['#e74c3c', '#f39c12', '#95a5a6', '#27ae60']
    ax.bar(removal_counts.keys(), removal_counts.values(), color=colors, edgecolor='black')
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('Cell Filtering Summary', fontsize=14, fontweight='bold')
    
    # Add count labels on bars
    for i, (key, val) in enumerate(removal_counts.items()):
        ax.text(i, val + adata.n_obs*0.01, f'{val:,}', 
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'filtering_summary.pdf', dpi=300)
    plt.savefig(OUTPUT_DIR / 'filtering_summary.png', dpi=150)
    plt.close()
    
    print("✓ Saved plots: filtering_summary.pdf/png")
    
except Exception as e:
    print(f"⚠️  Could not generate plots: {e}")
    print("   (Not critical, data cleaning was successful)")

# ==============================================================================
# Final Summary
# ==============================================================================

print("\n" + "=" * 70)
print("CLEANING COMPLETE")
print("=" * 70)

print(f"""
Summary:
  Original cells:     {adata.n_obs:>10,}
  Cleaned cells:      {adata_clean.n_obs:>10,}
  Removed cells:      {adata.n_obs - adata_clean.n_obs:>10,} ({(1-adata_clean.n_obs/adata.n_obs)*100:.2f}%)

Breakdown:
  - Doublets removed: {n_doublets:>10,}
  - GSE299751 cells:  {(adata.obs['dataset'] == 'GSE299751').sum() if 'GSE299751' in adata.obs['dataset'].values else 0:>10,}
  - Overlap:          {overlap:>10,}

Output files:
  - {output_h5ad.name}
  - cleaning_summary.json
  - remaining_datasets.csv
  - filtering_summary.pdf (if generated)

Next steps:
  1. Load cleaned data for re-integration
  2. Verify cell type distributions are reasonable
  3. Re-run BBKNN/scVI integration if needed
  4. Proceed with downstream analysis
""")

print("✓ All done! Your data is ready for re-integration.")
