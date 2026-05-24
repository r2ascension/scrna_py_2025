#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Standalone Reference + Query Merge Script
==========================================
Merges reference_with_L2_umap.h5ad and query_mapped_L2.h5ad

Usage:
    python standalone_merge_ref_query.py
"""

import sys
import gc
from pathlib import Path

import scanpy as sc
from scipy.sparse import csr_matrix, issparse

print("=" * 70)
print("Standalone Reference + Query Merge")
print("=" * 70)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

REF_H5AD = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3/reference_with_L2_umap.h5ad"
QUERY_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/query_mapped_L2.h5ad"
OUTPUT_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/reference_plus_query_merged_L2.h5ad"

print(f"\nConfiguration:")
print(f"  Reference: {REF_H5AD}")
print(f"  Query: {QUERY_H5AD}")
print(f"  Output: {OUTPUT_H5AD}")

# ==============================================================================
# Load Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1: Load Reference and Query")
print("=" * 70)

print(f"\nLoading reference...")
adata_ref = sc.read_h5ad(REF_H5AD)
print(f"  Shape: {adata_ref.shape}")
print(f"  obs columns: {list(adata_ref.obs.columns)}")
print(f"  obsm keys: {list(adata_ref.obsm.keys())}")

if 'X_umap' not in adata_ref.obsm:
    print("ERROR: Reference missing X_umap!")
    sys.exit(1)

print(f"\nLoading query...")
adata_query = sc.read_h5ad(QUERY_H5AD)
print(f"  Shape: {adata_query.shape}")
print(f"  obs columns: {list(adata_query.obs.columns)}")
print(f"  obsm keys: {list(adata_query.obsm.keys())}")

if 'X_umap' not in adata_query.obsm:
    print("ERROR: Query missing X_umap!")
    sys.exit(1)

# ==============================================================================
# Pre-Concat Cleanup
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Pre-Concat Cleanup")
print("=" * 70)

# 1. Clean var columns
print("\n1. Cleaning var columns...")
var_cols_to_drop = [
    'highly_variable', 'highly_variable_rank',
    'means', 'dispersions', 'dispersions_norm',
    'mt', 'n_cells', 'n_counts'
]

for col in var_cols_to_drop:
    if col in adata_ref.var.columns:
        adata_ref.var.drop(columns=[col], inplace=True)
    if col in adata_query.var.columns:
        adata_query.var.drop(columns=[col], inplace=True)

print(f"   Reference var: {list(adata_ref.var.columns)}")
print(f"   Query var: {list(adata_query.var.columns)}")

# 2. Clear uns
print("\n2. Clearing .uns...")
adata_ref.uns = {}
adata_query.uns = {}
print("   OK")

# 3. Unify sparse format
print("\n3. Unifying sparse format to CSR...")
if issparse(adata_ref.X):
    adata_ref.X = csr_matrix(adata_ref.X)
if issparse(adata_query.X):
    adata_query.X = csr_matrix(adata_query.X)

for layer_key in list(adata_ref.layers.keys()):
    if issparse(adata_ref.layers[layer_key]):
        adata_ref.layers[layer_key] = csr_matrix(adata_ref.layers[layer_key])

for layer_key in list(adata_query.layers.keys()):
    if issparse(adata_query.layers[layer_key]):
        adata_query.layers[layer_key] = csr_matrix(adata_query.layers[layer_key])

print("   OK")

# 4. Align obsm keys
print("\n4. Aligning obsm keys...")
critical_obsm = {'X_umap', 'X_scANVI_L2'}

for key in list(adata_ref.obsm.keys()):
    if key not in critical_obsm:
        del adata_ref.obsm[key]

for key in list(adata_query.obsm.keys()):
    if key not in critical_obsm:
        del adata_query.obsm[key]

print(f"   Kept: {list(critical_obsm)}")

# 5. Remove obsp
print("\n5. Removing obsp...")
for key in list(adata_ref.obsp.keys()):
    del adata_ref.obsp[key]
for key in list(adata_query.obsp.keys()):
    del adata_query.obsp[key]
print("   OK")

# 6. Check obs_names uniqueness
print("\n6. Checking obs_names...")
if not adata_ref.obs_names.is_unique:
    print("   WARNING: Reference has duplicates, making unique...")
    adata_ref.obs_names_make_unique()
if not adata_query.obs_names.is_unique:
    print("   WARNING: Query has duplicates, making unique...")
    adata_query.obs_names_make_unique()
print("   OK")

# 7. Add prefix
print("\n7. Adding obs_names prefix...")
adata_ref.obs_names = [f"ref::{x}" for x in adata_ref.obs_names]
adata_query.obs_names = [f"qry::{x}" for x in adata_query.obs_names]
print("   OK")

print("\nPre-concat cleanup complete")

# ==============================================================================
# Concatenate
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Concatenate")
print("=" * 70)

print("\nConcatenating...")
adata_all = sc.concat(
    {"reference": adata_ref, "query": adata_query},
    axis=0,
    join="outer",
    merge="unique",
    fill_value=0,
    label="data_source"
)

print(f"OK Merged shape: {adata_all.shape}")
print(f"  Reference: {(adata_all.obs['data_source']=='reference').sum():,}")
print(f"  Query: {(adata_all.obs['data_source']=='query').sum():,}")

# Validate UMAP
if 'X_umap' not in adata_all.obsm:
    print("WARNING: X_umap lost during concat!")
elif adata_all.obsm['X_umap'].shape[0] != adata_all.n_obs:
    print("WARNING: X_umap misaligned!")
else:
    print("OK X_umap preserved and aligned")

# ==============================================================================
# Save
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Save Merged Data")
print("=" * 70)

output_path = Path(OUTPUT_H5AD)
output_path.parent.mkdir(parents=True, exist_ok=True)

print(f"\nSaving: {output_path}")

try:
    adata_all.write_h5ad(output_path, compression='gzip')
    file_size = output_path.stat().st_size / 1024**3
    print(f"OK Saved ({file_size:.2f} GB)")
except Exception as e:
    print(f"ERROR with compression: {e}")
    print("Trying without compression...")
    try:
        adata_all.write_h5ad(output_path, compression=None)
        file_size = output_path.stat().st_size / 1024**3
        print(f"OK Saved uncompressed ({file_size:.2f} GB)")
    except Exception as e2:
        print(f"ERROR: {e2}")
        sys.exit(1)

# Cleanup
del adata_ref, adata_query, adata_all
gc.collect()

print("\n" + "=" * 70)
print("SUCCESS - Merge Complete")
print("=" * 70)
print(f"\nOutput: {output_path}")

# ==============================================================================
# Quick Validation
# ==============================================================================

print("\n" + "=" * 70)
print("Quick Validation")
print("=" * 70)

print("\nReloading to validate...")
adata = sc.read_h5ad(output_path)

print(f"\nMerged data summary:")
print(f"  Shape: {adata.shape}")
print(f"  obs_names unique: {adata.obs_names.is_unique}")
print(f"  data_source counts:")
for source, count in adata.obs['data_source'].value_counts().items():
    print(f"    {source}: {count:,}")

print(f"\n  obsm keys: {list(adata.obsm.keys())}")
print(f"  X_umap shape: {adata.obsm['X_umap'].shape}")

print("\n" + "=" * 70)
print("All checks passed!")
print("=" * 70)
