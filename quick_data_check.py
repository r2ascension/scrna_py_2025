#!/usr/bin/env python3
"""
Quick AnnData Structure Checker
快速检查h5ad文件的数据结构和raw counts状态

Usage: 
  方式1（命令行）: python quick_data_check.py <path_to_h5ad>
  方式2（直接运行）: 修改下面的H5AD_PATH后直接运行
"""

# ==================== 配置区域 ====================
# 直接修改这里的路径，然后运行脚本
H5AD_PATH = '/home/h2048/data/R/1215/merge/cleaned_samples_COMPLETE_v4.1/merged_seurat_standardized_filterd.h5ad'

# 采样参数
N_SAMPLE_CELLS = 10000  # 采样细胞数
N_SAMPLE_GENES = 10000    # 采样基因数

# ==================== 程序主体 ====================

import sys
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse

def check_raw_counts(X, name):
    """检查是否为raw counts"""
    if X is None:
        return None, None, None, "N/A"
    
    # 采样
    n_cells = X.shape[0]
    n_genes = X.shape[1]
    
    sample_cells = min(N_SAMPLE_CELLS, n_cells)
    sample_genes = min(N_SAMPLE_GENES, n_genes)
    
    sample = X[:sample_cells, :sample_genes]
    if issparse(sample):
        sample = sample.toarray()
    
    # 统计
    max_val = np.max(sample)
    mean_val = np.mean(sample[sample > 0]) if np.any(sample > 0) else 0
    
    # 判断（更严格的标准）
    is_raw = False
    confidence = "Unknown"
    
    if max_val > 50:  # 明确的raw counts
        is_raw = True
        confidence = "High"
    elif max_val > 20 and mean_val > 2:  # 可能是raw counts
        is_raw = True
        confidence = "Medium"
    elif max_val < 15:  # 明确的log-transformed
        is_raw = False
        confidence = "High"
    else:
        confidence = "Low"
    
    status = "✓ RAW COUNTS" if is_raw else "✗ log-transformed"
    
    return max_val, mean_val, status, confidence

def check_metadata(adata):
    """检查metadata信息"""
    print("\n" + "="*70)
    print("Metadata Check")
    print("="*70)
    
    # 检查obs (细胞metadata)
    print(f"\nCell metadata (.obs): {adata.n_obs:,} cells")
    print(f"  Columns: {len(adata.obs.columns)}")
    
    if len(adata.obs.columns) > 0:
        print(f"\n  Available columns:")
        for col in adata.obs.columns:
            dtype = adata.obs[col].dtype
            n_unique = adata.obs[col].nunique()
            
            # 分类还是连续变量
            if dtype == 'object' or dtype.name == 'category':
                var_type = "categorical"
                # 如果类别少，显示分布
                if n_unique <= 10:
                    counts = adata.obs[col].value_counts().head(5)
                    dist = ", ".join([f"{k}:{v}" for k, v in counts.items()])
                    print(f"    • {col:<30} {var_type:<12} (n={n_unique}) - {dist}")
                else:
                    print(f"    • {col:<30} {var_type:<12} (n={n_unique})")
            else:
                var_type = "continuous"
                mean_val = adata.obs[col].mean() if pd.api.types.is_numeric_dtype(adata.obs[col]) else "N/A"
                if isinstance(mean_val, (int, float)):
                    print(f"    • {col:<30} {var_type:<12} (mean={mean_val:.2f})")
                else:
                    print(f"    • {col:<30} {var_type:<12}")
    
    # 检查var (基因metadata)
    print(f"\n\nGene metadata (.var): {adata.n_vars:,} genes")
    print(f"  Columns: {len(adata.var.columns)}")
    
    if len(adata.var.columns) > 0:
        print(f"\n  Available columns:")
        for col in adata.var.columns:
            dtype = adata.var[col].dtype
            n_unique = adata.var[col].nunique()
            
            if dtype == 'object' or dtype.name == 'category':
                var_type = "categorical"
                print(f"    • {col:<30} {var_type:<12} (n={n_unique})")
            elif dtype == 'bool':
                n_true = adata.var[col].sum()
                pct = n_true / adata.n_vars * 100
                print(f"    • {col:<30} {'boolean':<12} ({n_true}/{adata.n_vars} = {pct:.1f}%)")
            else:
                var_type = "continuous"
                mean_val = adata.var[col].mean() if pd.api.types.is_numeric_dtype(adata.var[col]) else "N/A"
                if isinstance(mean_val, (int, float)):
                    print(f"    • {col:<30} {var_type:<12} (mean={mean_val:.2f})")
                else:
                    print(f"    • {col:<30} {var_type:<12}")
    
    # 检查obsm (embeddings等)
    if len(adata.obsm.keys()) > 0:
        print(f"\n\nEmbeddings (.obsm): {len(adata.obsm.keys())} items")
        for key in adata.obsm.keys():
            shape = adata.obsm[key].shape
            print(f"    • {key:<30} shape={shape}")
    
    # 检查uns (其他信息)
    if len(adata.uns.keys()) > 0:
        print(f"\n\nOther info (.uns): {len(adata.uns.keys())} items")
        for key in list(adata.uns.keys())[:20]:  # 只显示前20个
            print(f"    • {key}")
        if len(adata.uns.keys()) > 20:
            print(f"    ... and {len(adata.uns.keys()) - 20} more")

def main(h5ad_path):
    print("="*70)
    print("AnnData Structure Check")
    print("="*70)
    print(f"\nFile: {h5ad_path}")
    
    # 读取数据
    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception as e:
        print(f"\n❌ Error reading file: {e}")
        return
    
    print(f"\nShape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    # Metadata检查
    check_metadata(adata)
    
    # 检查各个来源
    print("\n" + "="*70)
    print("Data Sources Check")
    print("="*70)
    
    results = []
    
    # 检查 .raw
    if adata.raw is not None:
        max_v, mean_v, status, conf = check_raw_counts(adata.raw.X, ".raw")
        results.append({
            'source': '.raw',
            'exists': True,
            'max': max_v,
            'mean': mean_v,
            'status': status,
            'confidence': conf
        })
    else:
        results.append({
            'source': '.raw',
            'exists': False,
            'max': None,
            'mean': None,
            'status': 'N/A',
            'confidence': 'N/A'
        })
    
    # 检查 .X
    max_v, mean_v, status, conf = check_raw_counts(adata.X, ".X")
    results.append({
        'source': '.X',
        'exists': True,
        'max': max_v,
        'mean': mean_v,
        'status': status,
        'confidence': conf
    })
    
    # 检查 layers
    for layer_name in ['counts', 'raw_counts', 'log1p']:
        if layer_name in adata.layers:
            max_v, mean_v, status, conf = check_raw_counts(
                adata.layers[layer_name], f"layers['{layer_name}']"
            )
            results.append({
                'source': f"layers['{layer_name}']",
                'exists': True,
                'max': max_v,
                'mean': mean_v,
                'status': status,
                'confidence': conf
            })
    
    # 打印结果表格
    print(f"\n{'Source':<25} {'Exists':<8} {'Max':<10} {'Mean':<10} {'Status':<20} {'Confidence':<12}")
    print("-"*90)
    
    for r in results:
        exists = "Yes" if r['exists'] else "No"
        max_str = f"{r['max']:.2f}" if r['max'] is not None else "N/A"
        mean_str = f"{r['mean']:.2f}" if r['mean'] is not None else "N/A"
        
        print(f"{r['source']:<25} {exists:<8} {max_str:<10} {mean_str:<10} {r['status']:<20} {r['confidence']:<12}")
    
    # 综合判断
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    
    has_raw = False
    raw_sources = []
    
    for r in results:
        if r['exists'] and r['status'] == "✓ RAW COUNTS" and r['confidence'] in ['High', 'Medium']:
            has_raw = True
            raw_sources.append(r['source'])
    
    if has_raw:
        print(f"\n✅ Raw counts AVAILABLE in: {', '.join(raw_sources)}")
        print("\nYour data is READY for:")
        print("  - Differential expression analysis")
        print("  - scVI/scANVI training")
        print("  - starCAT annotation")
        print("  - Any analysis requiring raw counts")
    else:
        print("\n⚠️  No clear raw counts detected")
        print("\nPossible reasons:")
        print("  1. Data only contains normalized/log-transformed values")
        print("  2. Raw counts stored under different layer name")
        print("  3. Need to load data with .raw preserved")
    
    # 数据流建议
    print("\n" + "="*70)
    print("Recommended Data Flow")
    print("="*70)
    
    if has_raw:
        print("\nFor downstream analysis, use:")
        if '.raw' in raw_sources:
            print("  - Raw counts: adata.raw.X")
        if any('counts' in s for s in raw_sources):
            counts_layer = [s for s in raw_sources if 'counts' in s][0]
            print(f"  - Raw counts: adata.{counts_layer}")
        print("  - Normalized: adata.X (if max < 15)")
    else:
        print("\n⚠️  Consider re-loading data to preserve raw counts")

if __name__ == "__main__":
    # 优先使用命令行参数，否则使用配置的路径
    if len(sys.argv) > 1:
        h5ad_path = sys.argv[1]
        print("Using path from command line argument")
    else:
        h5ad_path = H5AD_PATH
        print("Using path from configuration")
    
    print(f"Checking: {h5ad_path}\n")
    main(h5ad_path)