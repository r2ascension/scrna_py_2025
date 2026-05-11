#!/usr/bin/env python3
"""
Transfer raw counts from BBKNN h5ad to scVI h5ad

This script extracts true raw counts from the BBKNN integrated file
and adds them to the scVI integrated file for starCAT analysis.
"""

import scanpy as sc
import numpy as np
from scipy.sparse import issparse
from pathlib import Path
import time

print("="*70)
print("TRANSFER RAW COUNTS: BBKNN → scVI")
print("="*70)

# File paths
BBKNN_H5AD = "/home/h2048/data/py/1128/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"
SCVI_H5AD = "/home/h2048/data/py/1203/T_scvi_integration/adata_tcell_scvi_full.h5ad"
OUTPUT_H5AD = "/home/h2048/data/py/1203/T_scvi_integration/adata_tcell_scvi_full_with_raw.h5ad"

print(f"\nInput files:")
print(f"  BBKNN (source): {BBKNN_H5AD}")
print(f"  scVI (target):  {SCVI_H5AD}")
print(f"Output:")
print(f"  New file:       {OUTPUT_H5AD}")

# ==================== Step 1: Load BBKNN data ====================
print("\n" + "="*70)
print("Step 1: Loading BBKNN data (source of raw counts)")
print("="*70)

print(f"\nLoading BBKNN h5ad...")
adata_bbknn = sc.read_h5ad(BBKNN_H5AD)
print(f"✓ Loaded: {adata_bbknn.shape[0]:,} cells × {adata_bbknn.shape[1]:,} genes")

# Check for raw counts in BBKNN
print(f"\nChecking BBKNN data structure:")
print(f"  .raw: {'Yes' if adata_bbknn.raw is not None else 'No'}")
print(f"  .layers: {list(adata_bbknn.layers.keys()) if adata_bbknn.layers else 'None'}")

# Function to check if data is raw counts
def check_raw_counts(data, name):
    x_sample = data[:100, :100]
    if issparse(x_sample):
        x_sample = x_sample.toarray()
    
    x_max = np.max(x_sample)
    x_mean = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0
    
    print(f"\n{name} statistics:")
    print(f"  Max: {x_max:.2f}")
    print(f"  Mean (non-zero): {x_mean:.2f}")
    
    is_raw = x_max > 20 or (x_max > 10 and x_mean > 2)
    if is_raw:
        print(f"  → Appears to be RAW COUNTS ✓")
    else:
        print(f"  → May be log-transformed ⚠️")
    
    return is_raw

# Check different sources in BBKNN
raw_counts_source = None
raw_counts_data = None

# Priority 1: .raw
if adata_bbknn.raw is not None:
    is_raw = check_raw_counts(adata_bbknn.raw.X, "BBKNN .raw")
    if is_raw:
        raw_counts_source = ".raw"
        raw_counts_data = adata_bbknn.raw.X.copy()
        raw_var_names = adata_bbknn.raw.var_names.copy()

# Priority 2: .X
if raw_counts_source is None:
    is_raw = check_raw_counts(adata_bbknn.X, "BBKNN .X")
    if is_raw:
        raw_counts_source = ".X"
        raw_counts_data = adata_bbknn.X.copy()
        raw_var_names = adata_bbknn.var_names.copy()

if raw_counts_source is None:
    print("\n✗ ERROR: No raw counts found in BBKNN file!")
    print("BBKNN file also appears to be log-transformed.")
    print("\nPlease use the h5ad file BEFORE any normalization/log transformation.")
    import sys
    sys.exit(1)

print(f"\n✓ Found raw counts in BBKNN: {raw_counts_source}")
print(f"  Shape: {raw_counts_data.shape}")

# ==================== Step 2: Load scVI data ====================
print("\n" + "="*70)
print("Step 2: Loading scVI data (target)")
print("="*70)

print(f"\nLoading scVI h5ad...")
adata_scvi = sc.read_h5ad(SCVI_H5AD)
print(f"✓ Loaded: {adata_scvi.shape[0]:,} cells × {adata_scvi.shape[1]:,} genes")

print(f"\nscVI data structure:")
print(f"  obs columns: {len(adata_scvi.obs.columns)}")
print(f"  obsm keys: {list(adata_scvi.obsm.keys())}")
print(f"  layers: {list(adata_scvi.layers.keys()) if adata_scvi.layers else 'None'}")

# ==================== Step 3: Verify cell matching ====================
print("\n" + "="*70)
print("Step 3: Verifying cell matching")
print("="*70)

# Check if cell barcodes match
bbknn_barcodes = set(adata_bbknn.obs_names)
scvi_barcodes = set(adata_scvi.obs_names)

n_bbknn = len(bbknn_barcodes)
n_scvi = len(scvi_barcodes)
n_common = len(bbknn_barcodes & scvi_barcodes)

print(f"\nCell barcode matching:")
print(f"  BBKNN cells: {n_bbknn:,}")
print(f"  scVI cells:  {n_scvi:,}")
print(f"  Common:      {n_common:,}")
print(f"  Match rate:  {n_common/n_scvi*100:.1f}%")

if n_common < n_scvi * 0.95:
    print(f"\n⚠️  Warning: Only {n_common/n_scvi*100:.1f}% of cells match")
    print(f"   BBKNN and scVI may be from different subsets")
    
    response = input("\nContinue anyway? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("Aborted.")
        import sys
        sys.exit(0)

if n_common != n_scvi:
    print(f"\n⚠️  Cell counts differ: BBKNN has {n_bbknn:,}, scVI has {n_scvi:,}")
    print(f"   Will only transfer raw counts for {n_common:,} matching cells")

# ==================== Step 4: Transfer raw counts ====================
print("\n" + "="*70)
print("Step 4: Transferring raw counts")
print("="*70)

print(f"\nPreparing to transfer raw counts...")

# Find matching cells (preserve order from scVI)
scvi_indices = []
bbknn_indices = []

bbknn_barcode_to_idx = {barcode: idx for idx, barcode in enumerate(adata_bbknn.obs_names)}

for idx, barcode in enumerate(adata_scvi.obs_names):
    if barcode in bbknn_barcode_to_idx:
        scvi_indices.append(idx)
        bbknn_indices.append(bbknn_barcode_to_idx[barcode])

print(f"✓ Found {len(scvi_indices):,} matching cells")

# Extract raw counts for matching cells
print(f"\nExtracting raw counts from BBKNN...")
if issparse(raw_counts_data):
    raw_counts_subset = raw_counts_data[bbknn_indices, :]
else:
    raw_counts_subset = raw_counts_data[bbknn_indices, :]

print(f"✓ Extracted raw counts: {raw_counts_subset.shape}")

# ==================== Step 5: Match genes ====================
print("\n" + "="*70)
print("Step 5: Matching genes between BBKNN and scVI")
print("="*70)

print(f"\nGene matching:")
print(f"  BBKNN genes: {len(raw_var_names):,}")
print(f"  scVI genes:  {len(adata_scvi.var_names):,}")

# Find common genes
bbknn_genes = set(raw_var_names)
scvi_genes = set(adata_scvi.var_names)
common_genes = bbknn_genes & scvi_genes

print(f"  Common:      {len(common_genes):,}")
print(f"  Match rate:  {len(common_genes)/len(scvi_genes)*100:.1f}%")

if len(common_genes) < len(scvi_genes) * 0.8:
    print(f"\n⚠️  Warning: Only {len(common_genes)/len(scvi_genes)*100:.1f}% of genes match")

# Create gene mapping
bbknn_gene_to_idx = {gene: idx for idx, gene in enumerate(raw_var_names)}
scvi_gene_order = []
bbknn_gene_order = []

for scvi_idx, gene in enumerate(adata_scvi.var_names):
    if gene in bbknn_gene_to_idx:
        scvi_gene_order.append(scvi_idx)
        bbknn_gene_order.append(bbknn_gene_to_idx[gene])

print(f"\n✓ Will transfer {len(scvi_gene_order):,} genes")

# Reorder raw counts to match scVI gene order
print(f"\nReordering genes to match scVI...")
if issparse(raw_counts_subset):
    raw_counts_ordered = raw_counts_subset[:, bbknn_gene_order]
else:
    raw_counts_ordered = raw_counts_subset[:, bbknn_gene_order]

print(f"✓ Reordered: {raw_counts_ordered.shape}")

# Create full raw counts matrix (for all scVI cells)
print(f"\nCreating full raw counts matrix...")
from scipy.sparse import csr_matrix, lil_matrix

if issparse(raw_counts_ordered):
    # Use sparse matrix
    full_raw_counts = lil_matrix((adata_scvi.n_obs, len(scvi_gene_order)), dtype=raw_counts_ordered.dtype)
    full_raw_counts[scvi_indices, :] = raw_counts_ordered
    full_raw_counts = full_raw_counts.tocsr()
else:
    # Use dense matrix
    full_raw_counts = np.zeros((adata_scvi.n_obs, len(scvi_gene_order)), dtype=raw_counts_ordered.dtype)
    full_raw_counts[scvi_indices, :] = raw_counts_ordered

print(f"✓ Created: {full_raw_counts.shape}")

# Verify
x_sample = full_raw_counts[:100, :100]
if issparse(x_sample):
    x_sample = x_sample.toarray()
x_max = np.max(x_sample)
x_mean = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0

print(f"\nVerification (first 100×100):")
print(f"  Max: {x_max:.2f}")
print(f"  Mean (non-zero): {x_mean:.2f}")

if x_max > 20:
    print(f"  ✓ Looks like raw counts!")
else:
    print(f"  ⚠️  Still looks log-transformed (max={x_max:.2f})")
    print(f"  This might not work for starCAT...")

# ==================== Step 6: Add to scVI h5ad ====================
print("\n" + "="*70)
print("Step 6: Adding raw counts to scVI h5ad")
print("="*70)

# Create subset of scVI with matched genes
adata_scvi_subset = adata_scvi[:, [adata_scvi.var_names[i] for i in scvi_gene_order]].copy()

print(f"\nCreating new h5ad with raw counts...")
print(f"  Original scVI: {adata_scvi.shape}")
print(f"  Subset (matched genes): {adata_scvi_subset.shape}")

# Add raw counts as a layer
adata_scvi_subset.layers['raw_counts_from_bbknn'] = full_raw_counts

print(f"✓ Added layer: 'raw_counts_from_bbknn'")

# Also store full raw data in .raw
print(f"\nStoring raw counts in .raw...")

# Create AnnData for raw
from anndata import AnnData
adata_raw = AnnData(
    X=full_raw_counts,
    obs=adata_scvi_subset.obs.copy(),
    var=adata_scvi_subset.var.copy()
)
adata_scvi_subset.raw = adata_raw

print(f"✓ Stored in .raw: {adata_scvi_subset.raw.shape}")

# ==================== Step 7: Save ====================
print("\n" + "="*70)
print("Step 7: Saving new h5ad")
print("="*70)

print(f"\nSaving to: {OUTPUT_H5AD}")
adata_scvi_subset.write_h5ad(OUTPUT_H5AD)

# Get file size
file_size = Path(OUTPUT_H5AD).stat().st_size / (1024**2)
print(f"✓ Saved: {file_size:.1f} MB")

# ==================== Summary ====================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"\n✅ Raw counts successfully transferred!")

print(f"\nSource:")
print(f"  BBKNN file: {BBKNN_H5AD}")
print(f"  Raw counts from: {raw_counts_source}")

print(f"\nTarget:")
print(f"  scVI file: {SCVI_H5AD}")
print(f"  New file: {OUTPUT_H5AD}")

print(f"\nTransferred:")
print(f"  Cells: {len(scvi_indices):,} / {adata_scvi.n_obs:,} ({len(scvi_indices)/adata_scvi.n_obs*100:.1f}%)")
print(f"  Genes: {len(scvi_gene_order):,} / {len(adata_scvi.var_names):,} ({len(scvi_gene_order)/len(adata_scvi.var_names)*100:.1f}%)")

print(f"\nNew h5ad structure:")
print(f"  .layers['raw_counts_from_bbknn']: Raw counts for starCAT")
print(f"  .raw: Full raw data")

print(f"\nVerification:")
print(f"  Max value: {x_max:.2f}")
print(f"  Mean (non-zero): {x_mean:.2f}")
if x_max > 20:
    print(f"  Status: ✓ Ready for starCAT")
else:
    print(f"  Status: ⚠️  May still be log-transformed")

print(f"\n" + "="*70)
print("NEXT STEPS")
print("="*70)

print(f"\n1. Update starCAT script INPUT_H5AD:")
print(f"   INPUT_H5AD = \"{OUTPUT_H5AD}\"")

print(f"\n2. Run starCAT analysis:")
print(f"   python tcell_scvi_starcat_analysis.py")

print(f"\n3. starCAT will use .raw for raw counts")

print(f"\n" + "="*70)
print("✅ DONE!")
print("="*70)
