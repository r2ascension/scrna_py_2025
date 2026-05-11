#!/usr/bin/env python3
"""
bbknn_integration_optimized.py

BBKNN批次整合脚本 - 逐行执行版本
- 基于图的快速批次整合
- 高变基因筛选
- 批次整合效果评估
- 多分辨率聚类
- 增强可视化

作者:临床-生信团队
日期:2025-11-01
版本:v1.0
"""

import sys
import os
from pathlib import Path
import warnings
import numpy as np
import time
warnings.filterwarnings('ignore')

# ==================== 配置基础 ====================

# ========== 输入输出配置 ==========
INPUT_H5AD_PATH = "/home/h2048/data/R/1029/seurat_without_partial_dual_contamination_doublets_2_1029.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1029/bbknn_output_optimized"
OVERWRITE_EXISTING = True

# ========== BBKNN整合配置 ==========
BATCH_KEY = "dataset"  # 批次变量名称

# BBKNN核心参数
BBKNN_NEIGHBORS_WITHIN_BATCH = 5  # 每个批次内选择的邻居数量(建议3-5)
BBKNN_N_PCS = 50  # 使用的主成分数量

# ========== 高变基因筛选 ==========
USE_HVG = True  # 是否筛选高变基因
N_TOP_GENES = 3000  # 保留的高变基因数量
HVG_FLAVOR = "seurat_v3"  # 高变基因筛选方法

# ========== 基因过滤配置 ==========
MIN_CELLS_PER_GENE = 3

# ========== 标准化配置 ==========
NORMALIZE_TOTAL = True  # 是否进行总数标准化
TARGET_SUM = 1e4  # 标准化目标总数
LOG_TRANSFORM = True  # 是否进行log转换
SCALE_DATA = True  # 是否scale数据
MAX_VALUE = 10  # scale的最大值

# ========== PCA配置 ==========
N_PCS = 50  # 计算的主成分数量

# ========== 降维和可视化配置 ==========
RUN_UMAP = True
UMAP_MIN_DIST = 0.5
UMAP_N_NEIGHBORS = 100  # UMAP使用的邻居数(会被BBKNN的邻居图覆盖)

# ========== 聚类配置(多分辨率) ==========
RUN_CLUSTERING = True
LEIDEN_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8]  # 多个分辨率
DEFAULT_RESOLUTION = 0.4  # 默认使用的分辨率

# ========== 批次整合评估 ==========
EVALUATE_INTEGRATION = False  # 是否评估批次整合效果

# ========== 可视化变量 ==========
VISUALIZATION_VARS = [
    "dataset",
    "Annotation",
    "tissue_sampling_method",
]

# ========== 增强可视化 ==========
GENERATE_FACET_PLOTS = True  # 是否生成分批次显示图

VERBOSE = True

# ==================== 辅助函数 ====================

def log_msg(msg):
    """打印日志"""
    if VERBOSE:
        print(msg)


def log_step(step_num, step_name):
    """打印步骤标题"""
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)


# ==================== 主要功能函数 ====================

def load_anndata(h5ad_path):
    """读取h5ad数据"""
    import scanpy as sc
    
    log_step(1, "Loading Data")
    log_msg(f"\nReading file: {h5ad_path}")
    
    if not Path(h5ad_path).exists():
        raise FileNotFoundError(f"File not found: {h5ad_path}")
    
    # 读取数据
    adata = sc.read_h5ad(h5ad_path)
    
    log_msg(f"Data loaded successfully")
    log_msg(f"   Cells: {adata.n_obs:,}")
    log_msg(f"   Genes: {adata.n_vars:,}")
    
    # 显示obs列
    log_msg(f"\nAvailable metadata columns:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        log_msg(f"   - {col}: {n_unique} unique values")
    
    # 检查batch key
    if BATCH_KEY not in adata.obs.columns:
        raise ValueError(f"Batch key '{BATCH_KEY}' not found in adata.obs")
    
    # 显示batch分布
    log_msg(f"\nBatch distribution (key: {BATCH_KEY}):")
    batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")
    
    return adata


def preprocess_for_bbknn(adata):
    """预处理数据用于BBKNN"""
    import scanpy as sc
    
    log_step(2, "Preprocessing for BBKNN")
    
    # 复制数据避免修改原始数据
    adata = adata.copy()
    
    # 保存原始counts
    if 'counts' not in adata.layers:
        log_msg("\nSaving raw counts to layers['counts']...")
        adata.layers['counts'] = adata.X.copy()
    
    # 基因过滤
    log_msg(f"\nGene filtering (min_cells={MIN_CELLS_PER_GENE})...")
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    n_genes_after = adata.n_vars
    log_msg(f"   Genes before: {n_genes_before:,}")
    log_msg(f"   Genes after: {n_genes_after:,}")
    log_msg(f"   Removed: {n_genes_before - n_genes_after:,}")
    
    # 标准化
    if NORMALIZE_TOTAL:
        log_msg(f"\nNormalizing total counts (target_sum={TARGET_SUM})...")
        sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    
    # Log转换
    if LOG_TRANSFORM:
        log_msg("Applying log1p transformation...")
        sc.pp.log1p(adata)
    
    # 高变基因筛选
    if USE_HVG:
        log_msg(f"\nSelecting highly variable genes (n={N_TOP_GENES})...")
        log_msg(f"   Method: {HVG_FLAVOR}")
        
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=N_TOP_GENES,
            flavor=HVG_FLAVOR,
            batch_key=BATCH_KEY,
            subset=False  # 先不subset,保留所有基因信息
        )
        
        n_hvg = adata.var['highly_variable'].sum()
        log_msg(f"   Selected {n_hvg} highly variable genes")
        
        # 子集到高变基因
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
        log_msg(f"   Final genes for downstream analysis: {adata_hvg.n_vars:,}")
    else:
        adata_hvg = adata.copy()
    
    # Scale数据(仅对HVG)
    if SCALE_DATA:
        log_msg(f"\nScaling data (max_value={MAX_VALUE})...")
        sc.pp.scale(adata_hvg, max_value=MAX_VALUE)
    
    log_msg("\nPreprocessing completed")
    log_msg(f"   Final dimensions: {adata_hvg.n_obs:,} cells x {adata_hvg.n_vars:,} genes")
    
    return adata, adata_hvg


def run_pca(adata):
    """运行PCA降维"""
    import scanpy as sc
    
    log_msg(f"\nRunning PCA (n_comps={N_PCS})...")
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver='arpack')
    
    # 计算解释方差比例
    var_ratio = adata.uns['pca']['variance_ratio']
    cumsum_var = np.cumsum(var_ratio)
    
    log_msg(f"   PC1-10 explained variance: {cumsum_var[9]:.2%}")
    log_msg(f"   PC1-20 explained variance: {cumsum_var[19]:.2%}")
    log_msg(f"   PC1-50 explained variance: {cumsum_var[49]:.2%}")
    
    return adata


def run_bbknn_integration(adata):
    """运行BBKNN批次整合"""
    import bbknn
    
    log_step(3, "BBKNN Batch Integration")
    
    log_msg("\nBBKNN parameters:")
    log_msg(f"   batch_key: {BATCH_KEY}")
    log_msg(f"   neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    log_msg(f"   n_pcs: {BBKNN_N_PCS}")
    
    log_msg("\nRunning BBKNN...")
    log_msg("   BBKNN constructs batch-balanced k-nearest neighbor graph")
    log_msg("   This may take a few minutes depending on data size...")
    
    start_time = time.time()
    
    # 运行BBKNN
    # BBKNN会直接修改adata的neighbors信息
    bbknn.bbknn(
        adata,
        batch_key=BATCH_KEY,
        neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
        n_pcs=BBKNN_N_PCS,
        copy=False  # 直接修改adata
    )
    
    elapsed_time = time.time() - start_time
    log_msg(f"   BBKNN completed in {elapsed_time:.1f} seconds")
    
    log_msg("\n   Batch-balanced neighbor graph constructed")
    log_msg("   Ready for downstream analysis (UMAP, clustering)")
    
    return adata


def run_umap_and_clustering(adata):
    """运行UMAP降维和聚类"""
    import scanpy as sc
    
    log_step(4, "UMAP and Clustering")
    
    # UMAP降维
    if RUN_UMAP:
        log_msg("\nComputing UMAP...")
        log_msg(f"   min_dist: {UMAP_MIN_DIST}")
        log_msg("   Using BBKNN neighbor graph")
        
        # BBKNN已经计算了neighbors,直接用UMAP
        sc.tl.umap(adata, min_dist=UMAP_MIN_DIST)
        log_msg("   UMAP completed")
    
    # 多分辨率聚类
    if RUN_CLUSTERING:
        log_msg("\nRunning multi-resolution Leiden clustering...")
        log_msg(f"   Resolutions: {LEIDEN_RESOLUTIONS}")
        
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_bbknn_res{res}'
            sc.tl.leiden(adata, resolution=res, key_added=key)
            n_clusters = adata.obs[key].nunique()
            log_msg(f"   Resolution {res}: {n_clusters} clusters")
        
        # 设置默认聚类
        default_key = f'leiden_bbknn_res{DEFAULT_RESOLUTION}'
        if default_key in adata.obs.columns:
            adata.obs['leiden_bbknn'] = adata.obs[default_key]
            log_msg(f"\n   Default clustering: {default_key}")
    
    return adata


def generate_visualizations(adata, output_dir):
    """生成可视化图表"""
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    log_step(5, "Generating Visualizations")
    
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    if not RUN_UMAP:
        log_msg("Skipping visualizations (UMAP not computed)")
        return
    
    # 设置scanpy图表保存路径
    sc.settings.figdir = fig_dir
    
    # 基础UMAP图
    log_msg("\nGenerating UMAP plots...")
    for var in VISUALIZATION_VARS:
        if var in adata.obs.columns:
            log_msg(f"   - UMAP colored by {var}")
            sc.pl.umap(
                adata,
                color=var,
                show=False,
                title=f'UMAP - {var}',
                save=f'_{var}.png'
            )
    
    # 聚类结果可视化
    if RUN_CLUSTERING:
        log_msg("\nGenerating clustering plots...")
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_bbknn_res{res}'
            if key in adata.obs.columns:
                log_msg(f"   - UMAP colored by {key}")
                sc.pl.umap(
                    adata,
                    color=key,
                    show=False,
                    title=f'UMAP - Leiden (res={res})',
                    save=f'_{key}.png'
                )
    
    # 分批次显示(facet plots)
    if GENERATE_FACET_PLOTS and BATCH_KEY in adata.obs.columns:
        log_msg("\nGenerating facet plots...")
        batches = adata.obs[BATCH_KEY].unique()
        
        if len(batches) <= 10:
            try:
                log_msg(f"   - Facet plot by {BATCH_KEY}")
                
                n_batches = len(batches)
                n_cols = min(3, n_batches)
                n_rows = (n_batches + n_cols - 1) // n_cols
                
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
                if n_batches == 1:
                    axes = [axes]
                else:
                    axes = axes.flatten()
                
                color_var = 'leiden_bbknn' if 'leiden_bbknn' in adata.obs.columns else None
                
                for i, batch in enumerate(sorted(batches)):
                    adata_batch = adata[adata.obs[BATCH_KEY] == batch]
                    sc.pl.umap(
                        adata_batch,
                        color=color_var,
                        ax=axes[i],
                        show=False,
                        title=f'{batch}'
                    )
                
                # 隐藏多余的子图
                for i in range(n_batches, len(axes)):
                    axes[i].axis('off')
                
                plt.tight_layout()
                plt.savefig(fig_dir / f'umap_facet_by_{BATCH_KEY}.png', dpi=150, bbox_inches='tight')
                plt.close()
                
            except Exception as e:
                log_msg(f"   Warning: Facet plot failed - {e}")
        else:
            log_msg(f"   Skipping facet plots (too many batches: {len(batches)})")
    
    log_msg(f"\nFigures saved to: {fig_dir}")


def evaluate_integration(adata, output_dir):
    """评估批次整合效果"""
    import json
    from sklearn.metrics import silhouette_score
    
    log_step(6, "Evaluating Integration Quality")
    
    results = {}
    
    # Silhouette score (batch) - 使用PCA空间
    log_msg("\nComputing silhouette scores...")
    
    if BATCH_KEY in adata.obs.columns and 'X_pca' in adata.obsm:
        batch_labels = adata.obs[BATCH_KEY].astype('category').cat.codes
        
        # 使用PCA空间评估
        sil_batch = silhouette_score(adata.obsm['X_pca'][:, :BBKNN_N_PCS], batch_labels)
        results['silhouette_batch'] = float(sil_batch)
        log_msg(f"   Silhouette (batch): {sil_batch:.4f}")
        log_msg(f"   Interpretation: closer to 0 = better batch mixing")
    
    # 如果有细胞类型信息,计算biological conservation
    if 'cell_type' in adata.obs.columns and 'X_pca' in adata.obsm:
        log_msg("\nComputing biological conservation...")
        cell_type_labels = adata.obs['cell_type'].astype('category').cat.codes
        
        sil_bio = silhouette_score(adata.obsm['X_pca'][:, :BBKNN_N_PCS], cell_type_labels)
        results['silhouette_biology'] = float(sil_bio)
        log_msg(f"   Silhouette (cell type): {sil_bio:.4f}")
        log_msg(f"   Interpretation: closer to 1 = better biological separation")
    
    # 计算batch mixing entropy
    if 'leiden_bbknn' in adata.obs.columns:
        log_msg("\nComputing batch mixing metrics...")
        
        from scipy.stats import entropy
        
        # 对每个cluster计算batch分布的熵
        clusters = adata.obs['leiden_bbknn'].unique()
        entropies = []
        
        for cluster in clusters:
            cluster_mask = adata.obs['leiden_bbknn'] == cluster
            batch_dist = adata.obs.loc[cluster_mask, BATCH_KEY].value_counts(normalize=True)
            ent = entropy(batch_dist)
            entropies.append(ent)
        
        mean_entropy = np.mean(entropies)
        results['mean_cluster_batch_entropy'] = float(mean_entropy)
        log_msg(f"   Mean cluster batch entropy: {mean_entropy:.4f}")
        log_msg(f"   Interpretation: higher = better batch mixing within clusters")
    
    # 保存结果
    eval_path = output_dir / "integration_evaluation.json"
    with open(eval_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    log_msg(f"\nEvaluation results saved to: {eval_path}")
    
    return results


def analyze_unintegrated_data(adata_original, output_dir):
    """分析未整合数据(评估批次效应)"""
    import scanpy as sc
    
    log_step("3.5", "Analyzing Unintegrated Data (Baseline)")
    
    # 创建临时副本
    adata_temp = adata_original.copy()
    
    # 如果已经预处理过,使用现有数据
    if 'counts' in adata_temp.layers:
        log_msg("\nUsing preprocessed data...")
        adata_temp.X = adata_temp.layers['counts'].copy()
        
        # 标准化
        if NORMALIZE_TOTAL:
            sc.pp.normalize_total(adata_temp, target_sum=TARGET_SUM)
        if LOG_TRANSFORM:
            sc.pp.log1p(adata_temp)
    
    # 如果使用HVG,只用高变基因
    if USE_HVG and 'highly_variable' in adata_temp.var:
        log_msg("Subsetting to highly variable genes...")
        adata_temp = adata_temp[:, adata_temp.var['highly_variable']].copy()
    
    # Scale和PCA
    log_msg("\nRunning PCA on unintegrated data...")
    if SCALE_DATA:
        sc.pp.scale(adata_temp, max_value=MAX_VALUE)
    sc.tl.pca(adata_temp, n_comps=50)
    
    # 标准neighbors和UMAP(不使用BBKNN)
    log_msg("Computing standard neighbors and UMAP...")
    sc.pp.neighbors(adata_temp, n_neighbors=UMAP_N_NEIGHBORS, use_rep="X_pca")
    sc.tl.umap(adata_temp, min_dist=UMAP_MIN_DIST)
    
    # 绘图
    log_msg("\nGenerating unintegrated visualizations...")
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    sc.settings.figdir = fig_dir
    
    # UMAP by batch
    if BATCH_KEY in adata_temp.obs.columns:
        sc.pl.umap(
            adata_temp,
            color=BATCH_KEY,
            show=False,
            title='UMAP - Unintegrated (by batch)',
            save='_unintegrated_batch.png'
        )
    
    # UMAP by cell type (if available)
    if 'cell_type' in adata_temp.obs.columns:
        sc.pl.umap(
            adata_temp,
            color='cell_type',
            show=False,
            title='UMAP - Unintegrated (by cell type)',
            save='_unintegrated_celltype.png'
        )
    
    log_msg("   Unintegrated analysis completed")
    
    del adata_temp


def generate_summary_report(adata, output_dir, processing_time=None, eval_results=None):
    """生成分析总结报告"""
    from datetime import datetime
    
    log_step(7, "Generating Summary Report")
    
    report = []
    report.append("="*70)
    report.append("BBKNN Batch Integration - Analysis Summary")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if processing_time:
        report.append(f"Total processing time: {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    report.append("")
    
    # 数据信息
    report.append("[ Data Overview ]")
    report.append(f"  Cells: {adata.n_obs:,}")
    report.append(f"  Genes: {adata.n_vars:,}")
    report.append("")
    
    # 批次信息
    if BATCH_KEY in adata.obs.columns:
        report.append(f"[ Batch Distribution (key: {BATCH_KEY}) ]")
        batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
        for batch, count in batch_counts.items():
            pct = count / adata.n_obs * 100
            report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
        report.append(f"  Total batches: {len(batch_counts)}")
        report.append("")
    
    # 批次整合评估
    if eval_results:
        report.append("[ Integration Quality Metrics ]")
        if 'silhouette_batch' in eval_results:
            report.append(f"  Silhouette Score (batch): {eval_results['silhouette_batch']:.4f}")
            report.append("  (Interpretation: closer to 0 = better batch mixing)")
        if 'silhouette_biology' in eval_results:
            report.append(f"  Silhouette Score (cell type): {eval_results['silhouette_biology']:.4f}")
            report.append("  (Interpretation: closer to 1 = better biological separation)")
        if 'mean_cluster_batch_entropy' in eval_results:
            report.append(f"  Mean cluster batch entropy: {eval_results['mean_cluster_batch_entropy']:.4f}")
            report.append("  (Interpretation: higher = better batch mixing)")
        report.append("")
    
    # 多分辨率聚类
    if RUN_CLUSTERING:
        report.append("[ Multi-resolution Clustering Results ]")
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_bbknn_res{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                default_marker = " (default)" if res == DEFAULT_RESOLUTION else ""
                report.append(f"  Resolution {res}: {n_clusters} clusters{default_marker}")
        report.append("")
    
    # Embeddings
    report.append("[ Available Embeddings ]")
    for key in adata.obsm.keys():
        report.append(f"  - {key}: {adata.obsm[key].shape}")
    report.append("")
    
    # BBKNN配置
    report.append("[ BBKNN Configuration ]")
    report.append(f"  neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    report.append(f"  n_pcs: {BBKNN_N_PCS}")
    report.append(f"  batch_key: {BATCH_KEY}")
    report.append(f"  High variable genes: {USE_HVG} (n={N_TOP_GENES if USE_HVG else 'N/A'})")
    report.append("")
    
    # 输出文件
    report.append("[ Output Files ]")
    report.append(f"  - Data: {output_dir / 'adata_bbknn_integrated.h5ad'}")
    report.append(f"  - Figures: {output_dir / 'figures/'}*.png")
    if EVALUATE_INTEGRATION:
        report.append(f"  - Evaluation: {output_dir / 'integration_evaluation.json'}")
    report.append("")
    
    # 方法学说明
    report.append("[ Method Notes ]")
    report.append("  BBKNN (Batch Balanced k-Nearest Neighbors):")
    report.append("  - Graph-based batch correction method")
    report.append("  - Constructs batch-balanced neighbor graph")
    report.append("  - Fast, no GPU required")
    report.append("  - Suitable for well-separated batches")
    report.append("")
    
    report.append("="*70)
    report.append("Analysis completed successfully")
    report.append("="*70)
    
    report_text = '\n'.join(report)
    
    # 保存
    report_path = output_dir / "analysis_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    log_msg(f"\nReport saved to: {report_path}")
    log_msg("\n" + report_text)


# ==================== 主程序 ====================

def main():
    """主函数"""
    
    print("\n" + "="*70)
    print("BBKNN Batch Integration Pipeline (Optimized v1.0)")
    print("="*70)
    
    start_time = time.time()
    
    # 检查环境
    log_msg("\nChecking Python environment...")
    try:
        import scanpy as sc
        import bbknn
        log_msg(f"   scanpy: {sc.__version__}")
        #log_msg(f"   bbknn: {bbknn.__version__}")
    except ImportError as e:
        print(f"\nError: {e}", file=sys.stderr)
        print("\nInstallation command: pip install scanpy bbknn", file=sys.stderr)
        sys.exit(1)
    
    # 创建输出目录
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_msg(f"\nOutput directory: {output_dir}")
    
    # Step 1: 读取数据
    adata = load_anndata(INPUT_H5AD_PATH)
    
    # Step 2: 预处理
    adata_full, adata_hvg = preprocess_for_bbknn(adata)
    
    # Step 3.5: 分析未整合数据
    analyze_unintegrated_data(adata_full, output_dir)
    
    # Step 3: PCA
    adata_hvg = run_pca(adata_hvg)
    
    # Step 4: BBKNN整合
    adata_hvg = run_bbknn_integration(adata_hvg)
    
    # Step 5: UMAP和聚类
    if RUN_UMAP or RUN_CLUSTERING:
        adata_hvg = run_umap_and_clustering(adata_hvg)
    
    # Step 6: 可视化
    if RUN_UMAP:
        generate_visualizations(adata_hvg, output_dir)
    
    # Step 7: 评估
    eval_results = None
    if EVALUATE_INTEGRATION:
        eval_results = evaluate_integration(adata_hvg, output_dir)
    
    # 将结果传回完整数据集
    log_msg("\nTransferring results to full dataset...")
    for key in ['X_pca', 'X_umap']:
        if key in adata_hvg.obsm:
            adata_full.obsm[key] = adata_hvg.obsm[key]
    
    for col in adata_hvg.obs.columns:
        if col.startswith('leiden_bbknn'):
            adata_full.obs[col] = adata_hvg.obs[col]
    
    # 复制neighbors信息
    if 'neighbors' in adata_hvg.uns:
        adata_full.uns['neighbors'] = adata_hvg.uns['neighbors']
    if 'connectivities' in adata_hvg.obsp:
        adata_full.obsp['connectivities'] = adata_hvg.obsp['connectivities']
    if 'distances' in adata_hvg.obsp:
        adata_full.obsp['distances'] = adata_hvg.obsp['distances']
    
    # Step 8: 保存
    log_msg("\nSaving integrated data...")
    final_path = output_dir / "adata_bbknn_integrated.h5ad"
    
    log_msg("   Using gzip compression (level 9)...")
    adata_full.write_h5ad(final_path, compression='gzip', compression_opts=9)
    
    file_size = final_path.stat().st_size / (1024**3)
    log_msg(f"   Saved: {final_path} ({file_size:.2f} GB)")
    
    # Step 9: 生成报告
    total_time = time.time() - start_time
    generate_summary_report(adata_full, output_dir, total_time, eval_results)
    
    # 最终总结
    print("\n" + "="*70)
    print("All analyses completed successfully")
    print("="*70)
    print(f"\nOutput directory: {output_dir}")
    print(f"Data: {final_path} ({file_size:.2f} GB)")
    if RUN_UMAP:
        print(f"Figures: {output_dir / 'figures/'}*.png")
    if EVALUATE_INTEGRATION:
        print(f"Evaluation: {output_dir / 'integration_evaluation.json'}")
    print(f"\nTotal time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    
    print(f"\nKey parameters:")
    print(f"   BBKNN neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    print(f"   High variable genes: {N_TOP_GENES if USE_HVG else 'Not used'}")
    print(f"   Multi-resolution clustering: {LEIDEN_RESOLUTIONS}")
    print()
    
    return adata_full, output_dir


if __name__ == "__main__":
    try:
        adata, output_dir = main()
    except KeyboardInterrupt:
        print("\n\nUser interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)