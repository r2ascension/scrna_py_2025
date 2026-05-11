#!/usr/bin/env python3
"""
Quick verification of raw counts in h5ad files
"""

import scanpy as sc
import numpy as np
from scipy.sparse import issparse

def check_h5ad_raw_counts(h5ad_path):
    """Check if h5ad file has raw counts"""
    print(f"\n{'='*70}")
    print(f"File: {h5ad_path.split('/')[-1]}")
    print(f"{'='*70}")
    
    # Load
    adata = sc.read_h5ad(h5ad_path)
    print(f"\nShape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    # Check structure
    print(f"\nData structure:")
    print(f"  .raw: {'Yes' if adata.raw is not None else 'No'}")
    print(f"  .layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
    
    # Function to check if data is raw counts
    def check_data(data, name):
        if data.shape[0] == 0:
            return False, 0, 0
        
        sample_size = 2000
        x = data[:sample_size, :sample_size]
        
        if issparse(x):
            x = x.toarray()
        
        x_max = np.max(x)
        x_mean = np.mean(x[x > 0]) if np.any(x > 0) else 0
        
        is_raw = x_max > 20 or (x_max > 10 and x_mean > 2)
        
        return is_raw, x_max, x_mean
    
    # Check all sources
    print(f"\n{'Source':<25} {'Max':>10} {'Mean':>10} {'Status'}")
    print(f"{'-'*60}")
    
    # Check .raw
    if adata.raw is not None:
        is_raw, x_max, x_mean = check_data(adata.raw.X, ".raw")
        status = "✓ RAW COUNTS" if is_raw else "✗ log-transformed"
        print(f"{'.raw':<25} {x_max:>10.2f} {x_mean:>10.2f} {status}")
    
    # Check .X
    is_raw, x_max, x_mean = check_data(adata.X, ".X")
    status = "✓ RAW COUNTS" if is_raw else "✗ log-transformed"
    print(f"{'.X':<25} {x_max:>10.2f} {x_mean:>10.2f} {status}")
    
    # Check layers
    if adata.layers:
        for layer_name in ['raw_counts', 'counts']:
            if layer_name in adata.layers:
                is_raw, x_max, x_mean = check_data(adata.layers[layer_name], layer_name)
                status = "✓ RAW COUNTS" if is_raw else "✗ log-transformed"
                print(f"layers['{layer_name}']"[:25].ljust(25) + f" {x_max:>10.2f} {x_mean:>10.2f} {status}")
    
    # Summary
    has_raw_in_raw = False
    if adata.raw is not None:
        is_raw, _, _ = check_data(adata.raw.X, ".raw")
        has_raw_in_raw = is_raw
    
    print(f"\n{'='*70}")
    if has_raw_in_raw:
        print(f"✅ This file HAS raw counts in .raw")
        print(f"   → Can be used for starCAT")
    else:
        print(f"❌ This file DOES NOT have raw counts")
        print(f"   → Cannot be used for starCAT")
    
    return has_raw_in_raw


# File paths
file1 = "/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/annotation_results/adata_with_manual_annotations.h5ad"
file2 = "/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/adata_Epithelial_bbknn.h5ad"

print("\n" + "="*70)
print("CHECKING TWO H5AD FILES FOR RAW COUNTS")
print("="*70)

# Check both files
result1 = check_h5ad_raw_counts(file1)
# result2 = check_h5ad_raw_counts(file2)

# Final summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"\nFile 1 (annotated_corrected): {'✅ Has raw counts' if result1 else '❌ No raw counts'}")
# print(f"File 2 (tcell_bbknn_starcat): {'✅ Has raw counts' if result2 else '❌ No raw counts'}")

if result1 or result2:
    print(f"\n✅ At least one file has raw counts!")
    print(f"\nRecommended for starCAT:")
    if result1:
        print(f"  Use: {file1}")
    if result2:
        print(f"  Use: {file2}")
else:
    print(f"\n❌ Neither file has raw counts")
    print(f"   Need to find the original h5ad file before any processing")

print(f"\n" + "="*70)
