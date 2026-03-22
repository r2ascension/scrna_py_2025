#!/usr/bin/env python3
"""
Quick AnnData Structure Checker (Detailed Version)
快速检查h5ad文件的数据结构和raw counts状态
v2.0 - 显示metadata每列的所有可能值

Usage: 
  方式1（命令行）: python quick_data_check_detailed.py <path_to_h5ad>
  方式2（直接运行）: 修改下面的H5AD_PATH后直接运行
"""

# ==================== 配置区域 ====================
# 直接修改这里的路径，然后运行脚本
H5AD_PATH = '/home/h2048/data/py/0209/myeloid_validation_optimized/adata_myeloid_refined_optimized.h5ad'

# 采样参数
N_SAMPLE_CELLS = 10000  # 采样细胞数
N_SAMPLE_GENES = 10000  # 采样基因数

# 显示控制
MAX_CATEGORICAL_VALUES = 100  # 分类变量最多显示多少个唯一值（防止输出过长）
SHOW_CONTINUOUS_VALUES_THRESHOLD = 50  # 连续变量唯一值少于此数时也显示所有值

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
    """检查metadata信息，显示每列的所有可能值"""
    print("\n" + "="*90)
    print("Metadata Check (Detailed - All Possible Values)")
    print("="*90)
    
    # 检查obs (细胞metadata)
    print(f"\n{'='*90}")
    print(f"Cell metadata (.obs): {adata.n_obs:,} cells")
    print(f"  Total columns: {len(adata.obs.columns)}")
    print("="*90)
    
    if len(adata.obs.columns) > 0:
        for col in adata.obs.columns:
            dtype = adata.obs[col].dtype
            n_unique = adata.obs[col].nunique()
            n_missing = adata.obs[col].isna().sum()
            
            print(f"\n📊 Column: {col}")
            print(f"   Type: {dtype}")
            print(f"   Unique values: {n_unique:,}")
            if n_missing > 0:
                print(f"   Missing values: {n_missing:,} ({n_missing/adata.n_obs*100:.2f}%)")
            
            # 分类变量或布尔变量：显示所有唯一值及计数
            if dtype == 'object' or dtype.name == 'category' or dtype == 'bool':
                counts = adata.obs[col].value_counts().sort_index()
                
                if n_unique <= MAX_CATEGORICAL_VALUES:
                    print(f"   All values and counts:")
                    for val, count in counts.items():
                        pct = count / adata.n_obs * 100
                        print(f"      • {val}: {count:,} ({pct:.2f}%)")
                else:
                    print(f"   ⚠️ Too many unique values ({n_unique}), showing top {MAX_CATEGORICAL_VALUES}:")
                    for val, count in counts.head(MAX_CATEGORICAL_VALUES).items():
                        pct = count / adata.n_obs * 100
                        print(f"      • {val}: {count:,} ({pct:.2f}%)")
                    print(f"      ... and {n_unique - MAX_CATEGORICAL_VALUES} more values")
            
            # 连续变量：显示统计摘要
            else:
                if pd.api.types.is_numeric_dtype(adata.obs[col]):
                    series = adata.obs[col].dropna()
                    if len(series) > 0:
                        print(f"   Statistics:")
                        print(f"      • Min: {series.min():.4f}")
                        print(f"      • 25th percentile: {series.quantile(0.25):.4f}")
                        print(f"      • Median: {series.median():.4f}")
                        print(f"      • Mean: {series.mean():.4f}")
                        print(f"      • 75th percentile: {series.quantile(0.75):.4f}")
                        print(f"      • Max: {series.max():.4f}")
                        print(f"      • Std: {series.std():.4f}")
                        
                        # 如果唯一值不多（<threshold），也显示所有值
                        if n_unique <= SHOW_CONTINUOUS_VALUES_THRESHOLD:
                            print(f"   All values and counts:")
                            counts = adata.obs[col].value_counts().sort_index()
                            for val, count in counts.items():
                                pct = count / adata.n_obs * 100
                                print(f"      • {val:.4f}: {count:,} ({pct:.2f}%)")
                else:
                    print(f"   ⚠️ Non-numeric continuous variable")
    
    # 检查var (基因metadata)
    print(f"\n\n{'='*90}")
    print(f"Gene metadata (.var): {adata.n_vars:,} genes")
    print(f"  Total columns: {len(adata.var.columns)}")
    print("="*90)
    
    if len(adata.var.columns) > 0:
        for col in adata.var.columns:
            dtype = adata.var[col].dtype
            n_unique = adata.var[col].nunique()
            n_missing = adata.var[col].isna().sum()
            
            print(f"\n📊 Column: {col}")
            print(f"   Type: {dtype}")
            print(f"   Unique values: {n_unique:,}")
            if n_missing > 0:
                print(f"   Missing values: {n_missing:,} ({n_missing/adata.n_vars*100:.2f}%)")
            
            # 分类变量或布尔变量：显示所有唯一值及计数
            if dtype == 'object' or dtype.name == 'category' or dtype == 'bool':
                counts = adata.var[col].value_counts().sort_index()
                
                if n_unique <= MAX_CATEGORICAL_VALUES:
                    print(f"   All values and counts:")
                    for val, count in counts.items():
                        pct = count / adata.n_vars * 100
                        print(f"      • {val}: {count:,} ({pct:.2f}%)")
                else:
                    print(f"   ⚠️ Too many unique values ({n_unique}), showing top {MAX_CATEGORICAL_VALUES}:")
                    for val, count in counts.head(MAX_CATEGORICAL_VALUES).items():
                        pct = count / adata.n_vars * 100
                        print(f"      • {val}: {count:,} ({pct:.2f}%)")
                    print(f"      ... and {n_unique - MAX_CATEGORICAL_VALUES} more values")
            
            # 连续变量：显示统计摘要
            else:
                if pd.api.types.is_numeric_dtype(adata.var[col]):
                    series = adata.var[col].dropna()
                    if len(series) > 0:
                        print(f"   Statistics:")
                        print(f"      • Min: {series.min():.4f}")
                        print(f"      • 25th percentile: {series.quantile(0.25):.4f}")
                        print(f"      • Median: {series.median():.4f}")
                        print(f"      • Mean: {series.mean():.4f}")
                        print(f"      • 75th percentile: {series.quantile(0.75):.4f}")
                        print(f"      • Max: {series.max():.4f}")
                        print(f"      • Std: {series.std():.4f}")
                        
                        # 如果唯一值不多（<threshold），也显示所有值
                        if n_unique <= SHOW_CONTINUOUS_VALUES_THRESHOLD:
                            print(f"   All values and counts:")
                            counts = adata.var[col].value_counts().sort_index()
                            for val, count in counts.items():
                                pct = count / adata.n_vars * 100
                                print(f"      • {val:.4f}: {count:,} ({pct:.2f}%)")
                else:
                    print(f"   ⚠️ Non-numeric continuous variable")
    
    # 检查obsm (embeddings等)
    if len(adata.obsm.keys()) > 0:
        print(f"\n\n{'='*90}")
        print(f"Embeddings (.obsm): {len(adata.obsm.keys())} items")
        print("="*90)
        for key in adata.obsm.keys():
            shape = adata.obsm[key].shape
            print(f"   • {key:<30} shape={shape}")
    
    # 检查obsp (邻居图等)
    if len(adata.obsp.keys()) > 0:
        print(f"\n\n{'='*90}")
        print(f"Pairwise data (.obsp): {len(adata.obsp.keys())} items")
        print("="*90)
        for key in adata.obsp.keys():
            shape = adata.obsp[key].shape
            print(f"   • {key:<30} shape={shape}")
    
    # 检查uns (其他信息)
    if len(adata.uns.keys()) > 0:
        print(f"\n\n{'='*90}")
        print(f"Other info (.uns): {len(adata.uns.keys())} items")
        print("="*90)
        for key in list(adata.uns.keys())[:30]:  # 显示前30个
            print(f"   • {key}")
        if len(adata.uns.keys()) > 30:
            print(f"   ... and {len(adata.uns.keys()) - 30} more items")

def main(h5ad_path):
    print("="*90)
    print("AnnData Structure Check (Detailed Version)")
    print("="*90)
    print(f"\nFile: {h5ad_path}")
    
    # 读取数据
    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception as e:
        print(f"\n❌ Error reading file: {e}")
        return
    
    print(f"\nShape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    # Metadata检查（详细版本，显示所有可能值）
    check_metadata(adata)
    
    # 检查各个来源
    print("\n" + "="*90)
    print("Data Sources Check (Raw Counts Detection)")
    print("="*90)
    
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
    for layer_name in ['counts', 'raw_counts', 'log1p', 'data']:
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
    print("-"*100)
    
    for r in results:
        exists = "Yes" if r['exists'] else "No"
        max_str = f"{r['max']:.2f}" if r['max'] is not None else "N/A"
        mean_str = f"{r['mean']:.2f}" if r['mean'] is not None else "N/A"
        
        print(f"{r['source']:<25} {exists:<8} {max_str:<10} {mean_str:<10} {r['status']:<20} {r['confidence']:<12}")
    
    # 综合判断
    print("\n" + "="*90)
    print("Summary")
    print("="*90)
    
    has_raw = False
    raw_sources = []
    
    for r in results:
        if r['exists'] and r['status'] == "✓ RAW COUNTS" and r['confidence'] in ['High', 'Medium']:
            has_raw = True
            raw_sources.append(r['source'])
    
    if has_raw:
        print(f"\n✅ Raw counts AVAILABLE in: {', '.join(raw_sources)}")
        print("\n📋 Your data is READY for:")
        print("   • Differential expression analysis")
        print("   • scVI/scANVI training")
        print("   • CellTypist annotation")
        print("   • Any analysis requiring raw counts")
    else:
        print("\n⚠️  No clear raw counts detected")
        print("\n📋 Possible reasons:")
        print("   1. Data only contains normalized/log-transformed values")
        print("   2. Raw counts stored under different layer name")
        print("   3. Need to load data with .raw preserved")
    
    # 数据流建议
    print("\n" + "="*90)
    print("Recommended Data Flow")
    print("="*90)
    
    if has_raw:
        print("\n📂 For downstream analysis, use:")
        if '.raw' in raw_sources:
            print("   • Raw counts: adata.raw.X")
        if any('counts' in s for s in raw_sources):
            counts_layer = [s for s in raw_sources if 'counts' in s][0]
            print(f"   • Raw counts: adata.{counts_layer}")
        print("   • Normalized: adata.X (if max < 15)")
    else:
        print("\n⚠️  Consider re-loading data to preserve raw counts")

if __name__ == "__main__":
    # 优先使用命令行参数，否则使用配置的路径
    if len(sys.argv) > 1:
        h5ad_path = sys.argv[1]
        print("📁 Using path from command line argument")
    else:
        h5ad_path = H5AD_PATH
        print("📁 Using path from configuration")
    
    print(f"🔍 Checking: {h5ad_path}\n")
    main(h5ad_path)
