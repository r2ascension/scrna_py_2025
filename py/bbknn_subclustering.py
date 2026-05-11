#!/usr/bin/env python3
"""
bbknn_celltype_subanalysis_v1.1.py

BBKNN 注释后的细胞类型分层分析(整合修复版)
- 修复 HVG(seurat_v3) 在小样本 batch 下 LOESS 奇异报错的鲁棒回退
- 修复 BBKNN 小批次剔除后未赋值 neighbors_within_batch 的报错
- 统一文件名安全(避免 / 空格 导致保存失败)
- 稀疏原始 counts 保留，避免 densify 爆内存
- PCA 解释方差安全打印；BBKNN/Neighbors 使用可用 n_pcs
- UMAP/Leiden 增加 random_state 以可复现
- 图像输出统一 dpi；分面图颜色列回退

作者：临床-生信团队(整合修复 by ChatGPT)
日期：2025-11-07
版本：v1.1
"""

import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import re
from tqdm import tqdm
from scipy.sparse import issparse, csr_matrix
warnings.filterwarnings('ignore')

# ==================== 配置部分 ====================

# ========== 输入输出配置 ==========
input_h5ad_path = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
output_base_dir = "/home/h2048/data/py/1029/bbknn_celltype_analysis"

# ========== 原始counts来源配置 ==========
raw_counts_source = "auto"  # "auto" | "raw" | "layers" | "X"

# ========== 细胞类型筛选配置 ==========
target_cell_types = []  # 例如: ["Epithelial", "T", "Myeloid"]
min_cells_per_celltype = 100

# ========== 批次配置 ==========
batch_key = "dataset"
cell_type_key = "cell_type"

# ========== BBKNN整合配置 ==========
bbknn_neighbors_within_batch = 3
bbknn_n_pcs = 50

# ========== 小batch处理策略 ==========
remove_small_batches = False  # True: 移除细胞数不足的batch (推荐)
# 每个 batch 最少细胞数：至少 neighbors_within_batch + 1；同时不低于 3
min_cells_for_bbknn = max(bbknn_neighbors_within_batch + 1, 3)

# ========== 高变基因筛选 ==========
use_hvg = True
n_top_genes = 4000
hvg_flavor = "seurat_v3"  # 默认使用 seurat_v3；内部做鲁棒回退
# 若使用 batch_key 进行 HVG 计算，要求每个 batch 至少多少细胞(过小会导致 LOESS 奇异)
min_cells_hvg_batch = 20

# ========== 基因过滤 ==========
min_cells_per_gene = 3

# ========== 标准化/变换 ==========
normalize_total = True
target_sum = 1e4
log_transform = True
scale_data = True
max_value = 10

# ========== PCA/UMAP/聚类 ==========
n_pcs = 50
run_umap = True
umap_min_dist = 0.3
run_clustering = True
leiden_resolutions = [0.2, 0.4, 0.6, 0.8]
default_resolution = 0.4
random_state = 0  # 复现性

# ========== 可视化 ==========
visualization_vars = ["dataset"]
generate_facet_plots = True
dpi = 300
figure_format = "pdf"

verbose = True

# ==================== 工具函数 ====================

def _safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")

def log_msg(msg):
    if verbose:
        print(msg)

def log_step(step_num, step_name):
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)

def create_output_dir(cell_type):
    clean_name = _safe_name(cell_type)
    output_dir = Path(output_base_dir) / clean_name
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    return output_dir, fig_dir

# ==================== 预处理 ====================

def preprocess_celltype_data(adata_subset, cell_type):
    import scanpy as sc
    import numpy as np

    log_step(2, f"Preprocessing {cell_type}")
    adata = adata_subset.copy()

    # ===== 选择并保存原始 counts 到 layers['raw_counts'] =====
    log_msg("\n检查数据状态并选择counts来源...")
    log_msg(f"   配置: raw_counts_source = '{raw_counts_source}'")

    has_raw_counts = False
    raw_layer_key = 'raw_counts'

    def _store_raw_counts(matrix):
        if raw_layer_key not in adata.layers:
            if issparse(matrix):
                adata.layers[raw_layer_key] = matrix.copy()
            else:
                adata.layers[raw_layer_key] = csr_matrix(matrix)
        return adata.layers[raw_layer_key].copy()

    if raw_counts_source == "raw":
        if adata.raw is not None:
            log_msg("   ✓ 使用 adata.raw (用户指定)")
            adata.X = _store_raw_counts(adata.raw.X)
            has_raw_counts = True
        else:
            raise ValueError("配置要求使用 adata.raw，但数据中不存在。")

    elif raw_counts_source == "layers":
        if 'counts' in adata.layers:
            log_msg("   ✓ 使用 layers['counts'] (用户指定)")
            adata.X = _store_raw_counts(adata.layers['counts'])
            has_raw_counts = True
        else:
            raise ValueError("配置要求使用 layers['counts']，但数据中不存在。")

    elif raw_counts_source == "X":
        log_msg("   ✓ 直接使用 X (用户指定为原始counts)")
        adata.X = _store_raw_counts(adata.X)
        has_raw_counts = True

    else:  # auto
        log_msg("   使用自动检测模式...")
        if adata.raw is not None:
            log_msg("   ✓ 检测到 adata.raw，将使用raw counts")
            adata.X = _store_raw_counts(adata.raw.X)
            has_raw_counts = True
        elif 'counts' in adata.layers:
            log_msg("   ✓ 检测到 layers['counts']，将使用原始counts")
            adata.X = _store_raw_counts(adata.layers['counts'])
            has_raw_counts = True
        else:
            # 粗略判断 X 是否像 log 数据
            try:
                max_val = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
            except Exception:
                max_val = np.max(adata.X.A) if issparse(adata.X) else np.max(adata.X)
            if max_val < 20:
                log_msg(f"   ⚠️  似为已预处理/log数据 (max_value={max_val:.2f})；无 raw counts 可用")
                has_raw_counts = False
            else:
                log_msg(f"   似为原始counts (max_value={max_val:.2f})")
                adata.X = _store_raw_counts(adata.X)
                has_raw_counts = True

    if has_raw_counts and 'log1p' in adata.uns:
        log_msg("   已移除 adata.uns['log1p'] 以反映当前为原始counts")
        adata.uns.pop('log1p', None)

    # 基因过滤
    log_msg(f"\n基因过滤 (min_cells={min_cells_per_gene})...")
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    n_genes_after = adata.n_vars
    log_msg(f"   过滤前: {n_genes_before:,}")
    log_msg(f"   过滤后: {n_genes_after:,}")
    log_msg(f"   移除: {n_genes_before - n_genes_after:,}")

    # 标准化 & Log1p(仅当我们有 raw counts 且需要)
    if has_raw_counts:
        if normalize_total:
            log_msg(f"\n标准化总counts (target_sum={target_sum})...")
            sc.pp.normalize_total(adata, target_sum=target_sum)
        if log_transform:
            log_msg("Log1p转换...")
            sc.pp.log1p(adata)
    else:
        log_msg("\n⚠️  跳过标准化和log转换(无原始counts，假定数据已预处理)")

    # ===== HVG(鲁棒回退策略) =====
    if use_hvg:
        log_msg(f"\n筛选高变基因 (n={n_top_genes})...")
        log_msg(f"   方法: {hvg_flavor}")

        n_batches = adata.obs[batch_key].nunique() if batch_key in adata.obs.columns else 1
        batch_counts = adata.obs[batch_key].value_counts() if batch_key in adata.obs.columns else pd.Series({"all": adata.n_obs})
        min_cells_per_batch = batch_counts.min()

        # seurat_v3 期望使用 count 数据，优先用我们保留的 raw_counts 层
        layer_for_hvg = 'raw_counts' if (hvg_flavor == 'seurat_v3' and has_raw_counts) else None

        # 判定是否用 batch_key
        use_batch_key = (n_batches > 1) and (min_cells_per_batch >= min_cells_hvg_batch)
        if n_batches == 1:
            log_msg("   注意: 仅有1个batch, 不使用 batch_key")
        elif not use_batch_key:
            log_msg(f"   提示: 存在较小 batch (最小={min_cells_per_batch}), 为避免 LOESS 奇异，先不使用 batch_key")

        def _try_hvg(flavor, use_batch):
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor=flavor,
                subset=False,
                batch_key=(batch_key if use_batch else None),
                layer=(layer_for_hvg if flavor == 'seurat_v3' else None)
            )

        try:
            if hvg_flavor == 'seurat_v3':
                _try_hvg('seurat_v3', use_batch_key)
            else:
                _try_hvg(hvg_flavor, use_batch_key)
        except Exception as e:
            log_msg(f"   ⚠️  HVG(seurat_v3) 失败: {e}\n       → 回退至不带 batch_key 的 seurat_v3")
            try:
                _try_hvg('seurat_v3', False)
            except Exception as e2:
                log_msg(f"   ⚠️  seurat_v3 仍失败: {e2}\n       → 最终回退到 flavor='seurat' (基于 log 数据)")
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top_genes,
                    flavor='seurat',
                    subset=False
                )

        n_hvg = int(adata.var.get('highly_variable', pd.Series([], dtype=bool)).sum())
        log_msg(f"   选中 {n_hvg} 个高变基因")
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
        log_msg(f"   最终基因数: {adata_hvg.n_vars:,}")
    else:
        adata_hvg = adata.copy()

    # Scale
    if scale_data:
        log_msg(f"\nScale数据 (max_value={max_value})...")
        import scanpy as sc
        sc.pp.scale(adata_hvg, max_value=max_value)

    log_msg("\n预处理完成")
    log_msg(f"   最终维度: {adata_hvg.n_obs:,} cells x {adata_hvg.n_vars:,} genes")
    return adata, adata_hvg

# ==================== PCA ====================

def run_pca_for_celltype(adata):
    import scanpy as sc
    log_msg(f"\n运行PCA (n_comps={n_pcs})...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    var_ratio = adata.uns['pca']['variance_ratio']
    cumsum_var = np.cumsum(var_ratio)
    k10 = min(10, len(cumsum_var))
    k30 = min(30, len(cumsum_var))
    log_msg(f"   PC1-{k10}累计解释方差: {cumsum_var[k10-1]:.2%}")
    log_msg(f"   PC1-{k30}累计解释方差: {cumsum_var[k30-1]:.2%}")
    return adata

# ==================== BBKNN ====================

def run_bbknn_for_celltype(adata, cell_type):
    from bbknn import bbknn as bbknn_func
    import scanpy as sc

    log_step(3, f"BBKNN Integration for {cell_type}")

    # 默认值，防止未赋值
    bbknn_neighbors_within_batch_to_use = bbknn_neighbors_within_batch

    batch_counts = adata.obs[batch_key].value_counts()
    n_batches = len(batch_counts)

    log_msg(f"\nBatch信息:")
    log_msg(f"   Batch数量: {n_batches}")
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")

    # 识别小 batch
    small_batches = batch_counts[batch_counts < min_cells_for_bbknn]
    if len(small_batches) > 0:
        log_msg(f"\n⚠️  检测到 {len(small_batches)} 个batch细胞数 < {min_cells_for_bbknn}:")
        for batch, count in small_batches.items():
            log_msg(f"   {batch}: {count} cells")

        if remove_small_batches:
            log_msg(f"\n根据配置 (remove_small_batches=True)，移除这些batch的细胞...")
            cells_before = adata.n_obs
            adata = adata[~adata.obs[batch_key].isin(small_batches.index)].copy()
            cells_after = adata.n_obs
            log_msg(f"   移除前: {cells_before:,} cells")
            log_msg(f"   移除后: {cells_after:,} cells")
            log_msg(f"   移除: {cells_before - cells_after} cells ({(cells_before - cells_after)/cells_before*100:.2f}%)")

            # 更新 batch 统计
            batch_counts = adata.obs[batch_key].value_counts()
            n_batches = len(batch_counts)
            log_msg(f"\n更新后的Batch信息:")
            log_msg(f"   剩余batch数: {n_batches}")
            for batch, count in batch_counts.items():
                pct = count / adata.n_obs * 100
                log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")

            # 保持 neighbors_within_batch_to_use 为默认设定
            bbknn_neighbors_within_batch_to_use = bbknn_neighbors_within_batch
        else:
            min_cells = batch_counts.min()
            bbknn_neighbors_within_batch_to_use = min(bbknn_neighbors_within_batch, int(min_cells))
            log_msg(f"\n根据配置 (remove_small_batches=False)，动态调整 neighbors_within_batch -> {bbknn_neighbors_within_batch_to_use}")

    # 计算可用 PCs
    n_pcs_available = adata.obsm['X_pca'].shape[1]
    n_pcs_to_use = min(bbknn_n_pcs, n_pcs_available)

    if n_batches == 1:
        log_msg(f"\n⚠️  注意: 仅有1个batch，跳过BBKNN整合，使用标准 neighbors 计算")
        sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca', n_pcs=n_pcs_to_use)
    else:
        log_msg(f"\nBBKNN参数:")
        log_msg(f"   batch_key: {batch_key}")
        log_msg(f"   neighbors_within_batch: {bbknn_neighbors_within_batch_to_use}")
        log_msg(f"   n_pcs: {n_pcs_to_use}")
        start_time = time.time()
        bbknn_func(
            adata,
            batch_key=batch_key,
            neighbors_within_batch=bbknn_neighbors_within_batch_to_use,
            n_pcs=n_pcs_to_use,
            copy=False
        )
        log_msg(f"   BBKNN完成，耗时 {time.time() - start_time:.1f} 秒")

    return adata

# ==================== UMAP & 聚类 ====================

def run_umap_and_clustering_for_celltype(adata, cell_type):
    import scanpy as sc
    log_step(4, f"UMAP and Clustering for {cell_type}")

    if run_umap:
        log_msg("\n计算UMAP...")
        log_msg(f"   min_dist: {umap_min_dist}")
        sc.tl.umap(adata, min_dist=umap_min_dist, random_state=random_state)
        log_msg("   UMAP完成")

    if run_clustering:
        log_msg("\n运行多分辨率Leiden聚类...")
        log_msg(f"   分辨率: {leiden_resolutions}")
        for res in leiden_resolutions:
            key = f'leiden_{cell_type}_res{res}'
            sc.tl.leiden(adata, resolution=res, key_added=key, random_state=random_state)
            n_clusters = adata.obs[key].nunique()
            log_msg(f"   分辨率 {res}: {n_clusters} 个clusters")
        default_key = f'leiden_{cell_type}_res{default_resolution}'
        if default_key in adata.obs.columns:
            adata.obs[f'leiden_{cell_type}'] = adata.obs[default_key]
            log_msg(f"\n   默认聚类: {default_key}")

    return adata

# ==================== 可视化 ====================

def generate_visualizations_for_celltype(adata, cell_type, output_dir, fig_dir):
    import scanpy as sc
    log_step(5, f"Generating Visualizations for {cell_type}")

    if not run_umap:
        log_msg("跳过可视化(UMAP未计算)")
        return

    sc.settings.figdir = fig_dir
    sc.settings.dpi = dpi

    safe_cell = _safe_name(cell_type)

    # 1) 基础 UMAP
    log_msg("\n生成UMAP图...")
    for var in visualization_vars:
        if var in adata.obs.columns:
            log_msg(f"   - UMAP按 {var} 着色")
            sc.pl.umap(
                adata,
                color=var,
                show=False,
                title=f'{cell_type} - {var}',
                save=f'_{safe_cell}_{var}.{figure_format}'
            )

    # 2) 聚类图
    if run_clustering:
        log_msg("\n生成聚类图...")
        for res in leiden_resolutions:
            key = f'leiden_{cell_type}_res{res}'
            if key in adata.obs.columns:
                sc.pl.umap(
                    adata,
                    color=key,
                    show=False,
                    title=f'{cell_type} - Leiden (res={res})',
                    save=f'_{safe_cell}_leiden_res{res}.{figure_format}'
                )

    # 3) 分批次分面图
    if generate_facet_plots and batch_key in adata.obs.columns:
        batches = adata.obs[batch_key].unique()
        if 1 < len(batches) <= 10:
            log_msg("\n生成分批次UMAP图...")
            try:
                n_batches = len(batches)
                n_cols = min(3, n_batches)
                n_rows = (n_batches + n_cols - 1) // n_cols
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
                axes = np.array(axes).reshape(-1) if n_batches > 1 else np.array([axes])

                color_var = f'leiden_{cell_type}' if f'leiden_{cell_type}' in adata.obs.columns else (batch_key if batch_key in adata.obs.columns else None)

                for i, batch in enumerate(sorted(batches)):
                    adata_batch = adata[adata.obs[batch_key] == batch]
                    sc.pl.umap(
                        adata_batch,
                        color=color_var,
                        ax=axes[i],
                        show=False,
                        title=f'{batch}'
                    )
                for i in range(n_batches, len(axes)):
                    axes[i].axis('off')
                plt.tight_layout()
                plt.savefig(fig_dir / f'umap_facet_{safe_cell}_by_batch.{figure_format}', dpi=dpi, bbox_inches='tight')
                plt.close()
                log_msg("   分批次UMAP图已保存")
            except Exception as e:
                log_msg(f"   ⚠️  分批次图生成失败: {e}")

    log_msg(f"\n可视化已保存至: {fig_dir}")

# ==================== 摘要 ====================

def generate_celltype_summary(adata, cell_type, output_dir, processing_time):
    from datetime import datetime
    log_step(6, f"Generating Summary for {cell_type}")

    report = []
    report.append("="*70)
    report.append(f"{cell_type} - BBKNN Subanalysis Summary")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Processing time: {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    report.append("")

    report.append("[ Data Overview ]")
    report.append(f"  Cell type: {cell_type}")
    report.append(f"  Cells: {adata.n_obs:,}")
    report.append(f"  Genes: {adata.n_vars:,}")
    report.append("")

    if batch_key in adata.obs.columns:
        report.append(f"[ Batch Distribution ]")
        batch_counts = adata.obs[batch_key].value_counts().sort_index()
        for batch, count in batch_counts.items():
            pct = count / adata.n_obs * 100
            report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
        report.append(f"  Total batches: {len(batch_counts)}")
        report.append("")

    if run_clustering:
        report.append("[ Clustering Results ]")
        for res in leiden_resolutions:
            key = f'leiden_{cell_type}_res{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                default_marker = " (default)" if res == default_resolution else ""
                report.append(f"  Resolution {res}: {n_clusters} clusters{default_marker}")
        report.append("")

    report.append("[ Available Embeddings ]")
    for key in adata.obsm.keys():
        report.append(f"  - {key}: {adata.obsm[key].shape}")
    report.append("")

    report.append("[ Analysis Configuration ]")
    report.append(f"  BBKNN neighbors_within_batch: {bbknn_neighbors_within_batch}")
    report.append(f"  BBKNN n_pcs: {bbknn_n_pcs}")
    report.append(f"  High variable genes: {use_hvg} (n={n_top_genes if use_hvg else 'N/A'})")
    report.append(f"  Clustering resolutions: {leiden_resolutions}")
    report.append("")

    safe_cell = _safe_name(cell_type)
    report.append("[ Output Files ]")
    report.append(f"  - Data: {output_dir / f'adata_{safe_cell}_bbknn.h5ad'}")
    report.append(f"  - Figures: {output_dir / 'figures'}/*.{figure_format}")
    report.append("")

    report.append("="*70)
    report.append("Analysis completed successfully")
    report.append("="*70)

    report_text = '\n'.join(report)
    report_path = output_dir / f"{safe_cell}_analysis_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)

    log_msg(f"\n摘要已保存: {report_path}")
    log_msg("\n" + report_text)

# ==================== 单细胞类型处理 ====================

def process_single_celltype(adata_full, cell_type):
    print("\n" + "="*70)
    print(f"🔬 Processing Cell Type: {cell_type}")
    print("="*70)

    start_time = time.time()
    output_dir, fig_dir = create_output_dir(cell_type)
    log_msg(f"\n输出目录: {output_dir}")

    # Step 1: 子集
    log_step(1, f"Extracting {cell_type} Cells")
    adata_subset = adata_full[adata_full.obs[cell_type_key] == cell_type].copy()

    log_msg(f"\n提取的细胞数: {adata_subset.n_obs:,}")
    log_msg(f"基因数: {adata_subset.n_vars:,}")

    if batch_key in adata_subset.obs.columns:
        log_msg(f"\nBatch分布:")
        batch_counts = adata_subset.obs[batch_key].value_counts().sort_index()
        for batch, count in batch_counts.items():
            pct = count / adata_subset.n_obs * 100
            log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")

    # Step 2: 预处理
    adata_full_processed, adata_hvg = preprocess_celltype_data(adata_subset, cell_type)

    # Step 3: PCA
    adata_hvg = run_pca_for_celltype(adata_hvg)

    # Step 4: BBKNN/Neighbors
    adata_hvg = run_bbknn_for_celltype(adata_hvg, cell_type)

    # Step 5: UMAP/聚类
    if run_umap or run_clustering:
        adata_hvg = run_umap_and_clustering_for_celltype(adata_hvg, cell_type)

    # Step 6: 可视化
    if run_umap:
        generate_visualizations_for_celltype(adata_hvg, cell_type, output_dir, fig_dir)

    # 回写嵌入与聚类到子集(保持与原逻辑一致)
    log_msg("\n将结果传回完整数据集...")
    for key in ['X_pca', 'X_umap']:
        if key in adata_hvg.obsm:
            adata_full_processed.obsm[key] = adata_hvg.obsm[key]
    for col in adata_hvg.obs.columns:
        if col.startswith(f'leiden_{cell_type}'):
            adata_full_processed.obs[col] = adata_hvg.obs[col]
    if 'neighbors' in adata_hvg.uns:
        adata_full_processed.uns['neighbors'] = adata_hvg.uns['neighbors']
    if 'connectivities' in adata_hvg.obsp:
        adata_full_processed.obsp['connectivities'] = adata_hvg.obsp['connectivities']
    if 'distances' in adata_hvg.obsp:
        adata_full_processed.obsp['distances'] = adata_hvg.obsp['distances']

    # Step 7: 保存数据
    log_msg("\n保存数据...")
    safe_cell = _safe_name(cell_type)
    output_file = output_dir / f"adata_{safe_cell}_bbknn.h5ad"
    log_msg("   使用gzip压缩 (level 9)...")
    adata_full_processed.write_h5ad(output_file, compression='gzip', compression_opts=9)
    file_size = output_file.stat().st_size / (1024**3)
    log_msg(f"   已保存: {output_file} ({file_size:.2f} GB)")

    # Step 8: 摘要
    processing_time = time.time() - start_time
    generate_celltype_summary(adata_full_processed, cell_type, output_dir, processing_time)

    print("\n" + "="*70)
    print(f"✅ {cell_type} Analysis Completed")
    print("="*70)
    print(f"Output: {output_dir}")
    print(f"Data: {output_file} ({file_size:.2f} GB)")
    if run_umap:
        print(f"Figures: {fig_dir}/*.{figure_format}")
    print(f"Time: {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    print()

    return adata_full_processed

# ==================== 主程序 ====================

def main():
    print("\n" + "="*70)
    print("🔬 BBKNN Cell Type Subanalysis Pipeline (v1.1)")
    print("="*70)

    overall_start_time = time.time()

    # 环境检查
    print("\n检查环境...")
    try:
        import scanpy as sc
        from bbknn import bbknn as bbknn_func
        print(f"   scanpy: {sc.__version__}")
    except ImportError as e:
        print(f"\n错误: {e}", file=sys.stderr)
        sys.exit(1)

    # 输出目录
    output_base = Path(output_base_dir)
    output_base.mkdir(parents=True, exist_ok=True)
    print(f"\n输出基础目录: {output_base}")

    # Step 1: 读取数据
    print("\n" + "="*70)
    print("Step 1: Loading Annotated Data")
    print("="*70)

    print(f"\n读取文件: {input_h5ad_path}")
    if not Path(input_h5ad_path).exists():
        raise FileNotFoundError(f"文件未找到: {input_h5ad_path}")

    import scanpy as sc
    adata = sc.read_h5ad(input_h5ad_path)

    print(f"数据加载成功")
    print(f"   Cells: {adata.n_obs:,}")
    print(f"   Genes: {adata.n_vars:,}")

    if cell_type_key not in adata.obs.columns:
        raise ValueError(f"细胞类型列 '{cell_type_key}' 未找到")
    if batch_key not in adata.obs.columns:
        raise ValueError(f"Batch列 '{batch_key}' 未找到")

    # Step 2: 筛选细胞类型
    print("\n" + "="*70)
    print("Step 2: Determining Cell Types to Analyze")
    print("="*70)

    celltype_counts = adata.obs[cell_type_key].value_counts()
    print(f"\n所有细胞类型及细胞数:")
    for ct, count in celltype_counts.items():
        pct = count / adata.n_obs * 100
        print(f"   {ct}: {count:,} cells ({pct:.1f}%)")

    if target_cell_types:
        celltypes_to_analyze = [ct for ct in target_cell_types if ct in celltype_counts.index]
        if len(celltypes_to_analyze) < len(target_cell_types):
            missing = set(target_cell_types) - set(celltypes_to_analyze)
            print(f"\n⚠️  指定的细胞类型未找到: {', '.join(missing)}")
    else:
        celltypes_to_analyze = [ct for ct, count in celltype_counts.items() if count >= min_cells_per_celltype]

    celltypes_to_analyze = sorted(celltypes_to_analyze, key=lambda x: celltype_counts[x], reverse=True)

    print(f"\n将分析的细胞类型 ({len(celltypes_to_analyze)}):")
    for ct in celltypes_to_analyze:
        count = celltype_counts[ct]
        pct = count / adata.n_obs * 100
        print(f"   ✓ {ct}: {count:,} cells ({pct:.1f}%)")

    if len(celltypes_to_analyze) == 0:
        print("\n⚠️  没有符合条件的细胞类型")
        sys.exit(0)

    # Step 3: 迭代处理
    print("\n" + "="*70)
    print("Step 3: Processing Each Cell Type")
    print("="*70)

    results = {}
    for i, cell_type in enumerate(celltypes_to_analyze, 1):
        print(f"\n{'='*70}")
        print(f"Processing {i}/{len(celltypes_to_analyze)}: {cell_type}")
        print(f"{'='*70}")
        try:
            adata_celltype = process_single_celltype(adata, cell_type)
            results[cell_type] = {
                'status': 'success',
                'n_cells': adata_celltype.n_obs,
                'output_dir': Path(output_base_dir) / _safe_name(cell_type)
            }
        except Exception as e:
            print(f"\n❌ 错误: {cell_type} 处理失败")
            print(f"   {e}")
            import traceback
            traceback.print_exc()
            results[cell_type] = {'status': 'failed', 'error': str(e)}

    # Step 4: 总结报告
    print("\n" + "="*70)
    print("Step 4: Generating Overall Summary")
    print("="*70)

    total_time = time.time() - overall_start_time
    from datetime import datetime

    report = []
    report.append("="*70)
    report.append("BBKNN Cell Type Subanalysis - Overall Summary")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Total processing time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    report.append("")

    report.append("[ Analysis Overview ]")
    report.append(f"  Total cell types analyzed: {len(celltypes_to_analyze)}")
    report.append(f"  Successful: {sum(1 for r in results.values() if r['status'] == 'success')}")
    report.append(f"  Failed: {sum(1 for r in results.values() if r['status'] == 'failed')}")
    report.append("")

    report.append("[ Results by Cell Type ]")
    for cell_type in celltypes_to_analyze:
        result = results[cell_type]
        if result['status'] == 'success':
            report.append(f"  ✓ {cell_type}:")
            report.append(f"    - Cells: {result['n_cells']:,}")
            report.append(f"    - Output: {result['output_dir']}")
        else:
            report.append(f"  ✗ {cell_type}: FAILED")
            report.append(f"    - Error: {result['error']}")
        report.append("")

    report.append("[ Configuration ]")
    report.append(f"  BBKNN neighbors_within_batch: {bbknn_neighbors_within_batch}")
    report.append(f"  BBKNN n_pcs: {bbknn_n_pcs}")
    report.append(f"  High variable genes: {n_top_genes if use_hvg else 'Not used'}")
    report.append(f"  Clustering resolutions: {leiden_resolutions}")
    report.append(f"  Min cells per cell type: {min_cells_per_celltype}")
    report.append(f"  Min cells per batch for HVG (batch_key): {min_cells_hvg_batch}")
    report.append("")

    report.append("="*70)
    report.append("Overall analysis completed")
    report.append("="*70)

    report_text = '\n'.join(report)
    report_path = output_base / "overall_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)

    print(f"\n总结报告已保存: {report_path}")
    print("\n" + report_text)

    print("\n" + "="*70)
    print("✅ All Cell Types Processed")
    print("="*70)
    print(f"\nOutput directory: {output_base}")
    print(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    print(f"\nSuccessful: {sum(1 for r in results.values() if r['status'] == 'success')}/{len(celltypes_to_analyze)}")
    if any(r['status'] == 'failed' for r in results.values()):
        print(f"\n⚠️  Some cell types failed. Check the summary for details.")
    print()

    return results

if __name__ == "__main__":
    try:
        results = main()
    except KeyboardInterrupt:
        print("\n\n⚠️  用户中断", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n❌ 错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
