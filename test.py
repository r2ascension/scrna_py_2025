#!/usr/bin/env python3
"""
Quick Data Status Check for starCAT Analysis
检查你的h5ad文件是否包含raw counts
"""

import scanpy as sc
import numpy as np
from scipy.sparse import issparse

# =============================================================================
# 配置：修改为你的文件路径
# =============================================================================
INPUT_FILE = "/home/h2048/data/py/1128/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"

# =============================================================================
# 加载数据
# =============================================================================
print("=" * 70)
print("Loading data...")
print("=" * 70)
adata = sc.read_h5ad(INPUT_FILE)
print(f"✓ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# =============================================================================
# 1. 检查数据结构
# =============================================================================
print("\n" + "=" * 70)
print("1. Data Structure")
print("=" * 70)

print(f"\n.raw exists: {adata.raw is not None}")
if adata.raw is not None:
    print(f"  → .raw shape: {adata.raw.shape}")

print(f"\n.layers exists: {bool(adata.layers)}")
if adata.layers:
    print(f"  → Available layers: {list(adata.layers.keys())}")
    if 'counts' in adata.layers:
        print(f"  → ✓ 'counts' layer found!")

print(f"\nTransformation markers in .uns:")
log_markers = [k for k in adata.uns.keys() if 'log' in k.lower()]
if log_markers:
    print(f"  → ⚠️  Found: {log_markers}")
else:
    print(f"  → None found")

# =============================================================================
# 2. 检查 .X 的统计特征
# =============================================================================
print("\n" + "=" * 70)
print("2. .X Statistics (Sample: first 1000 cells × 1000 genes)")
print("=" * 70)

# Sample data
n_sample_cells = min(1000, adata.n_obs)
n_sample_genes = min(1000, adata.n_vars)
x_sample = adata.X[:n_sample_cells, :n_sample_genes]

if issparse(x_sample):
    x_sample = x_sample.toarray()

# Calculate statistics
x_max = np.max(x_sample)
x_min = np.min(x_sample)
x_mean = np.mean(x_sample)
x_mean_nonzero = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0
x_std = np.std(x_sample)
x_nonzero_ratio = np.count_nonzero(x_sample) / x_sample.size

print(f"\nValue range:")
print(f"  Min:  {x_min:.4f}")
print(f"  Max:  {x_max:.2f}")
print(f"  Mean: {x_mean:.4f}")
print(f"  Mean (non-zero): {x_mean_nonzero:.2f}")
print(f"  Std:  {x_std:.4f}")

print(f"\nSparsity:")
print(f"  Non-zero ratio: {x_nonzero_ratio:.2%}")
print(f"  Zero ratio:     {(1-x_nonzero_ratio):.2%}")

# =============================================================================
# 3. 判断数据类型
# =============================================================================
print("\n" + "=" * 70)
print("3. Data Type Assessment")
print("=" * 70)

is_raw = False
is_log = False
is_normalized = False

# Check for raw counts characteristics
if x_max > 100:
    print("\n✓ HIGH MAX VALUE (>100)")
    print(f"  → Max = {x_max:.1f}")
    print(f"  → This is VERY LIKELY raw counts")
    is_raw = True
elif x_max > 20 and x_mean_nonzero > 1:
    print("\n✓ MODERATE VALUES")
    print(f"  → Max = {x_max:.1f}, Mean (non-zero) = {x_mean_nonzero:.2f}")
    print(f"  → This MIGHT BE raw counts (borderline)")
    is_raw = True
elif x_max < 20 and x_mean < 2:
    print("\n✗ LOW VALUES")
    print(f"  → Max = {x_max:.1f}, Mean = {x_mean:.4f}")
    print(f"  → This is LIKELY log-transformed")
    is_log = True
else:
    print("\n? UNCLEAR")
    print(f"  → Max = {x_max:.1f}, Mean = {x_mean:.4f}")
    print(f"  → Manual verification needed")

# Check for log transformation markers
if 'log1p' in adata.uns:
    print("\n⚠️  LOG TRANSFORMATION MARKER")
    print(f"  → 'log1p' found in .uns")
    print(f"  → Data was likely log-transformed")
    is_log = True

# =============================================================================
# 4. 给出建议
# =============================================================================
print("\n" + "=" * 70)
print("4. Recommendation for starCAT")
print("=" * 70)

if adata.raw is not None:
    print("\n✓✓✓ EXCELLENT")
    print("  → .raw exists")
    print("  → starCAT can use .raw.X")
    print("  → No action needed!")
    
elif 'counts' in adata.layers:
    print("\n✓✓ GOOD")
    print("  → .layers['counts'] exists")
    print("  → starCAT can use this layer")
    print("  → No action needed!")
    
elif is_raw and not is_log:
    print("\n✓ ACCEPTABLE")
    print("  → .X appears to contain raw counts")
    print("  → starCAT can try using .X")
    print("  → RECOMMENDATION: Verify manually if possible")
    
elif is_log:
    print("\n✗✗ PROBLEMATIC")
    print("  → .X appears to be log-transformed")
    print("  → starCAT requires raw counts (NOT log-transformed)")
    print("\n  OPTIONS:")
    print("  1. Provide original h5ad with .raw preserved")
    print("  2. Re-generate h5ad before log transformation")
    print("  3. Skip starCAT (set RUN_STARCAT = False)")
    print("  4. Continue with other analyses (BBKNN, clustering, markers)")
    
else:
    print("\n? UNCLEAR")
    print("  → Cannot determine data type with certainty")
    print("  → Suggest manual verification")

# =============================================================================
# 5. 检查具体的几个基因
# =============================================================================
print("\n" + "=" * 70)
print("5. Sample Gene Expression Values")
print("=" * 70)

# Look for common T cell markers
marker_genes = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B']
found_markers = [g for g in marker_genes if g in adata.var_names]

if found_markers:
    print(f"\nFound T cell markers: {found_markers}")
    print("\nExpression values (first 10 cells):")
    print("-" * 50)
    
    for gene in found_markers[:3]:  # Show first 3
        gene_idx = adata.var_names.get_loc(gene)
        expr = adata.X[:10, gene_idx]
        if issparse(expr):
            expr = expr.toarray().flatten()
        
        print(f"\n{gene}:")
        print(f"  Values: {expr}")
        print(f"  Max: {np.max(expr):.2f}, Mean: {np.mean(expr[expr>0]):.2f}")
        
        if np.max(expr) > 20:
            print(f"  → Looks like RAW counts")
        elif np.max(expr) < 10:
            print(f"  → Looks like LOG-transformed")
else:
    print("\nNo common T cell markers found in first check")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("\nData sources for starCAT (in priority order):")
print("1. .raw:              ", "✓ EXISTS" if adata.raw is not None else "✗ NOT FOUND")
print("2. .layers['counts']: ", "✓ EXISTS" if 'counts' in adata.layers else "✗ NOT FOUND")
print("3. .X (verified):     ", "✓ LIKELY RAW" if is_raw else "✗ LIKELY LOG-TRANSFORMED" if is_log else "? UNCLEAR")

print("\n" + "=" * 70)
if adata.raw is not None or 'counts' in adata.layers:
    print("✓✓✓ starCAT CAN RUN")
    print("=" * 70)
    print("You're all set! Run the v1.3 notebook.")
elif is_raw:
    print("✓ starCAT CAN PROBABLY RUN")
    print("=" * 70)
    print("The v1.3 notebook will attempt to use .X")
    print("Verify results carefully!")
else:
    print("✗ starCAT WILL BE SKIPPED")
    print("=" * 70)
    print("The v1.3 notebook will skip starCAT analysis")
    print("All other analyses will complete normally!")

print("\n" + "=" * 70)
print("CHECK COMPLETED")
print("=" * 70)