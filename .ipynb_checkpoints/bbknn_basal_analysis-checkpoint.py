#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bbknn_epithelial_subtype_analysis.py

针对特定上皮细胞亚型的深入分析
目标细胞类型: Goblet, Basal, Secretory, AT0

该脚本从注释后的上皮细胞数据中提取特定亚型,进行重新聚类和深入分析
重点关注各亚型的异质性和功能状态

作者: 临床-生信团队
日期: 2025-11-12
版本: v1.1 (添加BBKNN整合)
"""

import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import time
from tqdm import tqdm
warnings.filterwarnings('ignore')

# 尝试加载 BBKNN
try:
    import bbknn
    HAS_BBKNN = True
except Exception:
    HAS_BBKNN = False

# ==================== 配置部分 ====================

# ========== 输入输出配置 ==========
# 注释后的上皮细胞h5ad文件路径
INPUT_H5AD_PATH = "/home/h2048/data/py/1029/bbknn_annotation/Epithelial/res_leiden_bbknn_res2.4/epithelial_annotated_resleiden_bbknn_res2.4.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1029/bbknn_celltype_analysis/Epithelial/subtype_analysis"
OVERWRITE_EXISTING = True

# ========== 目标细胞类型配置 ==========
# 需要提取和分析的细胞类型
TARGET_CELL_TYPES = ["goblet", "basal", "secretory"]
# 注释列名 (根据实际数据调整)
ANNOTATION_KEY = "cell_type"  # 可能是 "cell_type", "Annotation" 等

# ========== 批次和组织信息 ==========
BATCH_KEY = "dataset"
TISSUE_KEY = "tissue"  # 组织来源
TISSUE_SAMPLING_KEY = "tissue_sampling_method"  # 组织采样方法

# 需要可视化的元数据列
METADATA_COLS_TO_PLOT = [BATCH_KEY, TISSUE_KEY, TISSUE_SAMPLING_KEY]

# ========== 数据预处理配置 ==========
# 基因过滤
MIN_CELLS = 10  # 基因至少在多少个细胞中表达
MIN_GENES = 200  # 细胞至少表达多少个基因

# 标准化
NORMALIZE_TOTAL = True
TARGET_SUM = 1e4
LOG_TRANSFORM = True

# 高变基因
USE_HVG = True
N_TOP_GENES = 4000
HVG_FLAVOR = "seurat_v3"

# Scaling
SCALE_DATA = True
MAX_VALUE = 10

# ========== 降维 / 邻居 / UMAP 配置 ==========
N_PCS = 50  # PCA主成分数
N_NEIGHBORS = 15  # 若不使用BBKNN时的邻居数 (用于sc.pp.neighbors)
UMAP_MIN_DIST = 0.3

# ========== BBKNN 配置（新增） ==========
USE_BBKNN = True                    # 是否启用 BBKNN
BBKNN_NEIGHBORS_WITHIN_BATCH = 5    # 每个batch内的邻居数
BBKNN_TRIM = 15                   # 可设为整数进行裁剪；None 表示不裁剪
BBKNN_METRIC = "cosine"          # 距离度量: 'euclidean'/'angular' 等

# 多分辨率聚类
LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
DEFAULT_RESOLUTION = 2.4

# ========== Marker基因分析 ==========
RUN_FIND_MARKERS = True
MARKER_MIN_PCT = 0.25
MARKER_LOGFC_THRESHOLD = 0.25
TOP_N_MARKERS = 10

# ========== 特定亚型的marker基因 ==========
SUBTYPE_MARKERS = {
    'Basal': {
        'core': ['TP63', 'KRT5', 'KRT14', 'KRT15'],
        'proliferation': ['MKI67', 'TOP2A', 'PCNA'],
        'differentiation': ['KRT8', 'KRT18', 'KRT19'],
        'stem': ['ITGA6', 'NGFR', 'SOX2']
    },
    'Secretory': {
        'club': ['SCGB1A1', 'SCGB3A1', 'SCGB3A2'],
        'mucous': ['MUC5B', 'MUC5AC'],
        'serous': ['LTF', 'LYZ', 'DMBT1'],
        'antimicrobial': ['BPIFA1', 'BPIFB1', 'PIGR']
    },
    'Goblet': {
        'mucin': ['MUC5AC', 'MUC5B', 'MUC2'],
        'secretory': ['TFF1', 'TFF3', 'SPDEF'],
        'processing': ['AGR2', 'LYPD2'],
        'surface': ['ITLN1', 'FCGBP']
    },
    'AT0': {
        # AT0是肺泡上皮祖细胞
        'progenitor': ['AGER', 'SFTPC', 'HOPX'],
        'at1_signature': ['AGER', 'PDPN', 'RTKN2'],
        'at2_signature': ['SFTPC', 'SFTPB', 'LAMP3'],
        'proliferation': ['MKI67', 'TOP2A']
    }
}

# ========== 功能评分基因集 ==========
FUNCTIONAL_GENE_SETS = {
    'Proliferation': ['MKI67', 'TOP2A', 'TK1', 'CENPW', 'PCNA'],
    'Stress_Response': ['HSPA1A', 'HSPA1B', 'HSPB1', 'HSP90AA1'],
    'Inflammatory': ['IL6', 'IL8', 'CXCL1', 'CXCL2', 'CXCL10'],
    'Ciliated_Transition': ['FOXJ1', 'RSPH1', 'DEUP1', 'FOXN4'],
    'EMT': ['VIM', 'CDH2', 'TWIST1', 'SNAI1', 'ZEB1'],
    'Metaplasia': ['KRT7', 'KRT13', 'KRT4']
}

# ========== 可视化配置 ==========
DPI = 300
FIGURE_FORMAT = "png"
COLOR_PALETTE = "tab10"

VERBOSE = True

# ==================== 辅助函数 ====================

def log_msg(msg):
    """打印日志信息"""
    if VERBOSE:
        print(msg)


def log_step(step_num, step_name):
    """打印步骤标题"""
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)


def setup_output_dir(output_dir, overwrite=False):
    """创建输出目录"""
    output_path = Path(output_dir)
    
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output directory exists: {output_dir}\n"
            f"Set OVERWRITE_EXISTING=True to overwrite"
        )
    
    output_path.mkdir(parents=True, exist_ok=True)
    log_msg(f"✓ Output directory ready: {output_dir}")
    return output_path


def check_data_layers(adata):
    """检查数据结构"""
    log_msg("\n数据结构检查:")
    log_msg(f"  .X: {adata.X.shape}, dtype={adata.X.dtype}")
    
    if hasattr(adata, 'layers') and len(adata.layers) > 0:
        log_msg(f"  .layers available:")
        for key in adata.layers.keys():
            log_msg(f"    - {key}: {adata.layers[key].shape}")
    else:
        log_msg("  .layers: None")
    
    if adata.raw is not None:
        log_msg(f"  .raw: {adata.raw.X.shape}")
    else:
        log_msg("  .raw: None")


def plot_initial_umap(adata, output_dir):
    """
    绘制数据加载后的初始UMAP图
    使用原始的UMAP坐标和注释
    """
    log_step("1b", "Plotting Initial UMAP (Before Re-analysis)")
    
    # 检查是否有UMAP坐标
    if 'X_umap' not in adata.obsm:
        log_msg("⚠️  警告: 数据中没有原始UMAP坐标,跳过初始可视化")
        return
    
    output_path = Path(output_dir)
    figures_dir = output_path / "figures_initial"
    figures_dir.mkdir(exist_ok=True)
    
    log_msg(f"\n初始UMAP图保存目录: {figures_dir}")
    
    # 1. 原始细胞类型
    log_msg("\n1. 原始细胞类型分布")
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=ANNOTATION_KEY,
        ax=ax,
        show=False,
        legend_loc='right margin',
        title='Initial Cell Types (Before Re-analysis)'
    )
    plt.tight_layout()
    plt.savefig(figures_dir / f"initial_umap_cell_types.{FIGURE_FORMAT}", dpi=DPI)
    plt.close()
    
    # 2-4. 元数据列 (dataset, tissue, tissue_sampling_method)
    for idx, col_name in enumerate(METADATA_COLS_TO_PLOT, start=2):
        if col_name not in adata.obs.columns:
            log_msg(f"⚠️  列 '{col_name}' 不存在,跳过")
            continue
        
        log_msg(f"{idx}. {col_name}")
        
        # 统计该列的值
        unique_vals = adata.obs[col_name].unique()
        n_unique = len(unique_vals)
        log_msg(f"   唯一值数量: {n_unique}")
        
        if n_unique > 50:
            log_msg(f"   ⚠️  类别过多({n_unique}),跳过可视化")
            continue
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata,
            color=col_name,
            ax=ax,
            show=False,
            legend_loc='right margin',
            title=f'Initial {col_name} (Before Re-analysis)'
        )
        plt.tight_layout()
        
        # 文件名处理
        safe_col_name = col_name.replace(' ', '_').replace('/', '_')
        plt.savefig(figures_dir / f"initial_umap_{safe_col_name}.{FIGURE_FORMAT}", dpi=DPI)
        plt.close()
    
    # 5. 组合图
    log_msg(f"{len(METADATA_COLS_TO_PLOT)+1}. 组合图")
    available_cols = [ANNOTATION_KEY] + [col for col in METADATA_COLS_TO_PLOT if col in adata.obs.columns]
    n_plots = len(available_cols)
    
    if n_plots > 0:
        ncols = 2
        nrows = (n_plots + 1) // 2
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*6, nrows*5))
        if n_plots == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for idx, col_name in enumerate(available_cols):
            n_unique = adata.obs[col_name].nunique()
            if n_unique > 50:
                axes[idx].text(0.5, 0.5, f'{col_name}\nToo many categories ({n_unique})',
                             ha='center', va='center', transform=axes[idx].transAxes)
                axes[idx].axis('off')
                continue
            
            sc.pl.umap(
                adata,
                color=col_name,
                ax=axes[idx],
                show=False,
                title=col_name,
                legend_loc='right margin',
                frameon=False
            )
        
        for idx in range(n_plots, len(axes)):
            fig.delaxes(axes[idx])
        
        plt.tight_layout()
        plt.savefig(figures_dir / f"initial_umap_combined.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
        plt.close()
    
    log_msg(f"\n✓ 初始UMAP图已保存到: {figures_dir}")
    log_msg(f"   共生成 {len(list(figures_dir.glob('*.png')))} 个图表")


# ==================== 主要分析函数 ====================

def load_and_subset_data(input_path, target_types, annotation_key):
    """
    读取数据并提取目标细胞类型
    """
    log_step(1, "Loading Data and Extracting Target Cell Types")
    
    log_msg(f"\n读取文件: {input_path}")
    if not Path(input_path).exists():
        raise FileNotFoundError(f"文件未找到: {input_path}")
    
    adata = sc.read_h5ad(input_path)
    log_msg(f"✓ 数据加载成功")
    log_msg(f"   总细胞数: {adata.n_obs:,}")
    log_msg(f"   总基因数: {adata.n_vars:,}")
    
    if annotation_key not in adata.obs.columns:
        available_keys = [col for col in adata.obs.columns if 'anno' in col.lower() or 'type' in col.lower()]
        raise KeyError(
            f"注释列 '{annotation_key}' 不存在\n"
            f"可用的注释相关列: {available_keys}\n"
            f"请修改ANNOTATION_KEY参数"
        )
    
    available_types = adata.obs[annotation_key].unique().tolist()
    log_msg(f"\n可用的细胞类型: {available_types}")
    
    missing_types = [t for t in target_types if t not in available_types]
    if missing_types:
        log_msg(f"\n⚠️  警告: 以下目标类型未找到: {missing_types}")
        log_msg("\n尝试模糊匹配...")
        matched_types = []
        for target in target_types:
            target_lower = target.lower()
            for avail in available_types:
                if target_lower in str(avail).lower() or str(avail).lower() in target_lower:
                    matched_types.append(avail)
                    log_msg(f"  '{target}' 匹配到 '{avail}'")
        if matched_types:
            target_types = list(dict.fromkeys(matched_types))  # 去重并保持顺序
            log_msg(f"\n✓ 使用匹配后的类型: {target_types}")
        else:
            raise ValueError(
                f"未找到目标细胞类型\n"
                f"目标: {target_types}\n"
                f"可用: {available_types}"
            )
    
    # 提取子集
    log_msg(f"\n提取目标细胞类型...")
    mask = adata.obs[annotation_key].isin(target_types)
    adata_subset = adata[mask].copy()
    
    log_msg(f"✓ 子集提取完成")
    log_msg(f"   提取细胞数: {adata_subset.n_obs:,}")
    
    type_counts = adata_subset.obs[annotation_key].value_counts()
    log_msg(f"\n各类型细胞数分布:")
    for cell_type, count in type_counts.items():
        pct = 100 * count / adata_subset.n_obs
        log_msg(f"   {cell_type}: {count:,} ({pct:.1f}%)")
    
    check_data_layers(adata_subset)
    return adata_subset


def preprocess_subset(adata, output_dir):
    """
    对子集数据进行预处理
    - 不将数据写入 .raw；以 layers['counts'] 和 layers['log1p'] 管理
    """
    log_step(2, "Preprocessing Subset Data")
    
    # 保存原始数据到layers (如果还没有)
    if 'counts' not in adata.layers:
        if adata.raw is not None:
            log_msg("从.raw复制原始counts...")
            adata.layers['counts'] = adata.raw.X.copy()
        else:
            log_msg("⚠️  警告: 没有原始counts,使用当前.X作为counts")
            adata.layers['counts'] = adata.X.copy()
    
    # 基因过滤
    log_msg(f"\n基因过滤 (min_cells={MIN_CELLS})...")
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    log_msg(f"✓ 保留基因数: {adata.n_vars:,}")
    
    # 细胞质控
    log_msg(f"\n细胞质控 (min_genes={MIN_GENES})...")
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    n_after = adata.n_obs
    log_msg(f"✓ 过滤前: {n_before:,} 细胞")
    log_msg(f"  过滤后: {n_after:,} 细胞")
    log_msg(f"  移除: {n_before - n_after:,} 细胞")
    
    # 标准化
    if NORMALIZE_TOTAL:
        log_msg(f"\n标准化 (target_sum={TARGET_SUM:.0e})...")
        adata.X = adata.layers['counts'].copy()
        sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
        if LOG_TRANSFORM:
            log_msg("Log1p转换...")
            sc.pp.log1p(adata)
        adata.layers['log1p'] = adata.X.copy()
        log_msg("✓ 标准化完成")
    
    # 高变基因筛选
    if USE_HVG:
        log_msg(f"\n高变基因筛选 (n_top_genes={N_TOP_GENES})...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=N_TOP_GENES,
            flavor=HVG_FLAVOR,
            batch_key=BATCH_KEY if BATCH_KEY in adata.obs.columns else None
        )
        n_hvg = int(adata.var['highly_variable'].sum())
        log_msg(f"✓ 高变基因数: {n_hvg}")
        
        # 强制保留关键marker基因
        all_markers = []
        for markers_dict in SUBTYPE_MARKERS.values():
            for marker_list in markers_dict.values():
                all_markers.extend(marker_list)
        all_markers = list(set(all_markers))
        markers_in_data = [m for m in all_markers if m in adata.var_names]
        log_msg(f"\n强制保留marker基因: {len(markers_in_data)}/{len(all_markers)}")
        if markers_in_data:
            adata.var.loc[markers_in_data, 'highly_variable'] = True
            n_hvg_final = int(adata.var['highly_variable'].sum())
            log_msg(f"✓ 最终高变基因数: {n_hvg_final}")
    
    # Scaling
    if SCALE_DATA:
        log_msg(f"\nScaling (max_value={MAX_VALUE})...")
        sc.pp.scale(adata, max_value=MAX_VALUE)
        log_msg("✓ Scaling完成")
    
    log_msg(f"\n预处理完成")
    log_msg(f"   最终细胞数: {adata.n_obs:,}")
    log_msg(f"   最终基因数: {adata.n_vars:,}")
    return adata


def _ensure_categorical_batch(adata, key):
    """保证 batch 列为分类类型"""
    if key in adata.obs.columns:
        if not pd.api.types.is_categorical_dtype(adata.obs[key]):
            adata.obs[key] = adata.obs[key].astype('category')


def run_dimensionality_reduction(adata):
    """
    降维分析: PCA + (BBKNN 或标准邻居) + UMAP
    """
    log_step(3, "Dimensionality Reduction")
    
    # PCA
    log_msg(f"\nPCA分析 (n_pcs={N_PCS})...")
    if USE_HVG:
        sc.tl.pca(adata, n_comps=N_PCS, use_highly_variable=True)
    else:
        sc.tl.pca(adata, n_comps=N_PCS)
    log_msg("✓ PCA完成")
    
    # 构建邻居图（BBKNN 优先）
    use_bbknn_now = (
        USE_BBKNN and
        HAS_BBKNN and
        (BATCH_KEY in adata.obs.columns) and
        (adata.obs[BATCH_KEY].nunique() > 1)
    )
    _ensure_categorical_batch(adata, BATCH_KEY)
    
    if use_bbknn_now:
        log_msg(f"\n使用 BBKNN 构建邻居图 (batch_key='{BATCH_KEY}', "
                f"neighbors_within_batch={BBKNN_NEIGHBORS_WITHIN_BATCH}, n_pcs={N_PCS})...")
        try:
            bbknn.bbknn(
                adata,
                batch_key=BATCH_KEY,
                neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
                n_pcs=N_PCS,
                metric=BBKNN_METRIC,
                trim=BBKNN_TRIM
            )
            adata.uns['integration_method'] = 'bbknn'
            log_msg("✓ BBKNN 邻居图构建完成")
        except Exception as e:
            log_msg(f"⚠️  BBKNN 失败，改用标准 neighbors。原因: {e}")
            sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
            adata.uns['integration_method'] = 'neighbors'
            log_msg("✓ 标准邻居图构建完成")
    else:
        if not HAS_BBKNN and USE_BBKNN:
            log_msg("⚠️  未安装 bbknn 包，将改用标准 neighbors。")
        elif BATCH_KEY not in adata.obs.columns or adata.obs[BATCH_KEY].nunique() <= 1:
            log_msg("⚠️  缺少有效的 batch 列或batch数量≤1，跳过BBKNN，使用标准 neighbors。")
        log_msg(f"\n构建标准邻居图 (n_neighbors={N_NEIGHBORS})...")
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        adata.uns['integration_method'] = 'neighbors'
        log_msg("✓ 标准邻居图构建完成")
    
    # UMAP
    log_msg(f"\nUMAP降维 (min_dist={UMAP_MIN_DIST})...")
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST)
    log_msg("✓ UMAP完成")
    return adata


def run_clustering(adata):
    """
    多分辨率Leiden聚类
    """
    log_step(4, "Multi-Resolution Clustering")
    
    log_msg(f"\nLeiden聚类 (resolutions: {LEIDEN_RESOLUTIONS})...")
    for res in tqdm(LEIDEN_RESOLUTIONS, desc="Clustering"):
        key = f'leiden_r{res}'
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        log_msg(f"  Resolution {res}: {n_clusters} clusters")
    
    # 设置默认聚类结果
    default_key = f'leiden_r{DEFAULT_RESOLUTION}'
    adata.obs['leiden'] = adata.obs[default_key]
    
    log_msg(f"\n✓ 聚类完成")
    log_msg(f"   默认分辨率: {DEFAULT_RESOLUTION}")
    log_msg(f"   默认cluster数: {adata.obs['leiden'].nunique()}")
    return adata


def calculate_functional_scores(adata):
    """
    计算功能评分
    """
    log_step(5, "Functional Scoring")
    log_msg("\n计算功能基因集评分...")
    
    for gene_set_name, genes in FUNCTIONAL_GENE_SETS.items():
        genes_in_data = [g for g in genes if g in adata.var_names]
        if len(genes_in_data) == 0:
            log_msg(f"⚠️  {gene_set_name}: 没有基因在数据中")
            continue
        if len(genes_in_data) < len(genes):
            log_msg(f"  {gene_set_name}: {len(genes_in_data)}/{len(genes)} 基因可用")
        score_key = f'score_{gene_set_name}'
        sc.tl.score_genes(adata, genes_in_data, score_name=score_key)
        log_msg(f"✓ {gene_set_name}: {score_key}")
    
    log_msg("\n✓ 功能评分计算完成")
    return adata


def find_cluster_markers(adata, output_dir):
    """
    查找cluster marker基因
    """
    if not RUN_FIND_MARKERS:
        log_msg("\n跳过marker基因分析 (RUN_FIND_MARKERS=False)")
        return adata
    
    log_step(6, "Finding Cluster Markers")
    
    # 确保使用log1p数据
    if 'log1p' in adata.layers:
        log_msg("使用log1p层进行差异分析...")
        adata.X = adata.layers['log1p'].copy()
    
    log_msg(f"\n查找差异基因 (method=wilcoxon)...")
    log_msg(f"   min_pct: {MARKER_MIN_PCT}")
    log_msg(f"   logfc_threshold: {MARKER_LOGFC_THRESHOLD}")
    
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby='leiden',
            method='wilcoxon',
            key_added='rank_genes_leiden',
            pts=True
        )
        log_msg("✓ 差异基因分析完成")
        
        # 保存结果
        marker_file = Path(output_dir) / "cluster_markers.xlsx"
        log_msg(f"\n保存marker基因到: {marker_file}")
        
        result = adata.uns['rank_genes_leiden']
        groups = result['names'].dtype.names
        
        try:
            with pd.ExcelWriter(marker_file, engine='openpyxl') as writer:
                for group in groups:
                    df = pd.DataFrame({
                        'gene': result['names'][group][:50],
                        'scores': result['scores'][group][:50],
                        'logfoldchanges': result['logfoldchanges'][group][:50],
                        'pvals': result['pvals'][group][:50],
                        'pvals_adj': result['pvals_adj'][group][:50],
                        'pct_in_group': result['pts'][group][:50]
                    })
                    df.to_excel(writer, sheet_name=f'Cluster_{group}', index=False)
            log_msg("✓ Marker基因已保存 (xlsx)")
        except Exception as e:
            # 如果没有 openpyxl，则退回 CSV
            csv_dir = Path(output_dir) / "markers_csv"
            csv_dir.mkdir(exist_ok=True)
            for group in groups:
                df = pd.DataFrame({
                    'gene': result['names'][group][:50],
                    'scores': result['scores'][group][:50],
                    'logfoldchanges': result['logfoldchanges'][group][:50],
                    'pvals': result['pvals'][group][:50],
                    'pvals_adj': result['pvals_adj'][group][:50],
                    'pct_in_group': result['pts'][group][:50]
                })
                df.to_csv(csv_dir / f"Cluster_{group}.csv", index=False)
            log_msg(f"⚠️  保存xlsx失败({e})，已改为CSV输出到: {csv_dir}")
        
    except Exception as e:
        log_msg(f"⚠️  差异分析失败: {e}")
    
    return adata


def generate_visualizations(adata, output_dir):
    """
    生成各种可视化图表
    """
    log_step(7, "Generating Visualizations (After Re-analysis)")
    
    output_path = Path(output_dir)
    figures_dir = output_path / "figures"
    figures_dir.mkdir(exist_ok=True)
    
    log_msg(f"\n图表保存目录: {figures_dir}")
    
    # 1. UMAP - 细胞类型
    log_msg("\n1. UMAP - 原始细胞类型")
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=ANNOTATION_KEY,
        ax=ax,
        show=False,
        legend_loc='right margin',
        title='Cell Types (After Re-analysis)'
    )
    plt.tight_layout()
    plt.savefig(figures_dir / f"umap_cell_types.{FIGURE_FORMAT}", dpi=DPI)
    plt.close()
    
    # 2. UMAP - Leiden聚类
    log_msg("2. UMAP - Leiden聚类")
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color='leiden',
        ax=ax,
        show=False,
        legend_loc='right margin',
        title='Leiden Clusters (After Re-analysis)'
    )
    plt.tight_layout()
    plt.savefig(figures_dir / f"umap_leiden.{FIGURE_FORMAT}", dpi=DPI)
    plt.close()
    
    # 3-5. UMAP - 元数据列
    plot_counter = 3
    for col_name in METADATA_COLS_TO_PLOT:
        if col_name not in adata.obs.columns:
            log_msg(f"⚠️  列 '{col_name}' 不存在,跳过")
            continue
        
        log_msg(f"{plot_counter}. UMAP - {col_name}")
        unique_vals = adata.obs[col_name].unique()
        n_unique = len(unique_vals)
        log_msg(f"   唯一值数量: {n_unique}")
        if n_unique > 50:
            log_msg(f"   ⚠️  类别过多({n_unique}),跳过可视化")
            plot_counter += 1
            continue
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata,
            color=col_name,
            ax=ax,
            show=False,
            legend_loc='right margin',
            title=f'{col_name} (After Re-analysis)'
        )
        plt.tight_layout()
        safe_col_name = col_name.replace(' ', '_').replace('/', '_')
        plt.savefig(figures_dir / f"umap_{safe_col_name}.{FIGURE_FORMAT}", dpi=DPI)
        plt.close()
        plot_counter += 1
    
    # 6. UMAP - 功能评分
    log_msg(f"{plot_counter}. UMAP - 功能评分")
    score_cols = [col for col in adata.obs.columns if col.startswith('score_')]
    if score_cols:
        n_scores = len(score_cols)
        ncols = 3
        nrows = (n_scores + ncols - 1) // ncols
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*5, nrows*4))
        axes = axes.flatten() if n_scores > 1 else [axes]
        
        for idx, score_col in enumerate(score_cols):
            sc.pl.umap(
                adata,
                color=score_col,
                ax=axes[idx],
                show=False,
                title=score_col.replace('score_', ''),
                cmap='viridis'
            )
        for idx in range(n_scores, len(axes)):
            fig.delaxes(axes[idx])
        
        plt.tight_layout()
        plt.savefig(figures_dir / f"umap_functional_scores.{FIGURE_FORMAT}", dpi=DPI)
        plt.close()
    plot_counter += 1
    
    # 7. 小提琴图 - 关键marker基因
    log_msg(f"{plot_counter}. 小提琴图 - 关键marker基因")
    all_markers = []
    for cell_type, markers_dict in SUBTYPE_MARKERS.items():
        for category, genes in markers_dict.items():
            all_markers.extend(genes)
    all_markers = list(set(all_markers))
    markers_in_data = [m for m in all_markers if m in adata.var_names]
    if markers_in_data:
        markers_to_plot = markers_in_data[:min(20, len(markers_in_data))]
        fig = plt.figure(figsize=(4*len(markers_to_plot), 6))
        sc.pl.violin(
            adata,
            keys=markers_to_plot,
            groupby='leiden',
            rotation=90,
            show=False
        )
        plt.tight_layout()
        plt.savefig(figures_dir / f"violin_markers.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
        plt.close()
    plot_counter += 1
    
    # 8. Dotplot - marker基因
    log_msg(f"{plot_counter}. Dotplot - marker基因表达")
    if markers_in_data:
        fig = plt.figure(figsize=(max(12, len(markers_in_data)*0.3), 8))
        sc.pl.dotplot(
            adata,
            var_names=markers_in_data,
            groupby='leiden',
            show=False,
            standard_scale='var'
        )
        plt.tight_layout()
        plt.savefig(figures_dir / f"dotplot_markers.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
        plt.close()
    plot_counter += 1
    
    # 9. 细胞类型组成热图
    log_msg(f"{plot_counter}. 细胞类型组成热图")
    ct_cluster = pd.crosstab(
        adata.obs[ANNOTATION_KEY],
        adata.obs['leiden'],
        normalize='columns'
    ) * 100
    fig, ax = plt.subplots(figsize=(max(10, ct_cluster.shape[1]*0.6), 6))
    sns.heatmap(
        ct_cluster,
        annot=True,
        fmt='.1f',
        cmap='YlOrRd',
        ax=ax,
        cbar_kws={'label': 'Percentage (%)'}
    )
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Cell Type')
    ax.set_title('Cell Type Composition by Cluster')
    plt.tight_layout()
    plt.savefig(figures_dir / f"heatmap_composition.{FIGURE_FORMAT}", dpi=DPI)
    plt.close()
    plot_counter += 1
    
    # 10. 细胞数统计柱状图
    log_msg(f"{plot_counter}. 细胞数统计")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    type_counts = adata.obs[ANNOTATION_KEY].value_counts()
    type_counts.plot(kind='bar', ax=axes[0], color='steelblue')
    axes[0].set_xlabel('Cell Type')
    axes[0].set_ylabel('Cell Count')
    axes[0].set_title('Cell Count by Type')
    axes[0].tick_params(axis='x', rotation=45)
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    cluster_counts.plot(kind='bar', ax=axes[1], color='coral')
    axes[1].set_xlabel('Leiden Cluster')
    axes[1].set_ylabel('Cell Count')
    axes[1].set_title('Cell Count by Cluster')
    plt.tight_layout()
    plt.savefig(figures_dir / f"barplot_cell_counts.{FIGURE_FORMAT}", dpi=DPI)
    plt.close()
    plot_counter += 1
    
    # 11. 多分辨率聚类比较
    log_msg(f"{plot_counter}. 多分辨率聚类比较")
    resolution_keys = [f'leiden_r{res}' for res in LEIDEN_RESOLUTIONS]
    existing_keys = [k for k in resolution_keys if k in adata.obs.columns]
    if len(existing_keys) > 1:
        n_res = len(existing_keys)
        ncols = min(3, n_res)
        nrows = (n_res + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*5, nrows*4))
        axes = axes.flatten() if n_res > 1 else [axes]
        for idx, res_key in enumerate(existing_keys):
            res_value = res_key.replace('leiden_r', '')
            sc.pl.umap(
                adata,
                color=res_key,
                ax=axes[idx],
                show=False,
                title=f'Resolution {res_value}',
                legend_loc='right margin'
            )
        for idx in range(n_res, len(axes)):
            fig.delaxes(axes[idx])
        plt.tight_layout()
        plt.savefig(figures_dir / f"umap_multi_resolution.{FIGURE_FORMAT}", dpi=DPI)
        plt.close()
    plot_counter += 1
    
    # 12. 分析前后对比组合图
    log_msg(f"{plot_counter}. 分析前后对比组合图")
    cols_to_show = [ANNOTATION_KEY, 'leiden']
    for col in METADATA_COLS_TO_PLOT:
        if col in adata.obs.columns and adata.obs[col].nunique() <= 50:
            cols_to_show.append(col)
    n_plots = len(cols_to_show)
    if n_plots > 0:
        ncols = 3
        nrows = (n_plots + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*5, nrows*4))
        if n_plots == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        for idx, col_name in enumerate(cols_to_show):
            sc.pl.umap(
                adata,
                color=col_name,
                ax=axes[idx],
                show=False,
                title=col_name,
                legend_loc='right margin',
                frameon=False
            )
        for idx in range(n_plots, len(axes)):
            fig.delaxes(axes[idx])
        plt.tight_layout()
        plt.savefig(figures_dir / f"umap_combined_after_analysis.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
        plt.close()
    
    log_msg(f"\n✓ 所有图表已保存到: {figures_dir}")
    log_msg(f"   共生成 {len(list(figures_dir.glob('*.png')))} 个图表")


def generate_summary_report(adata, output_dir, target_types):
    """
    生成分析总结报告
    """
    log_step(8, "Generating Summary Report")
    
    report_file = Path(output_dir) / "analysis_summary.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("上皮细胞亚型深入分析 - 总结报告\n")
        f.write("="*70 + "\n\n")
        f.write("1. 基本信息\n")
        f.write("-"*70 + "\n")
        f.write(f"分析日期: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"目标细胞类型: {', '.join(target_types)}\n")
        f.write(f"总细胞数: {adata.n_obs:,}\n")
        f.write(f"总基因数: {adata.n_vars:,}\n\n")
        f.write("2. 细胞类型分布\n")
        f.write("-"*70 + "\n")
        type_counts = adata.obs[ANNOTATION_KEY].value_counts()
        for cell_type, count in type_counts.items():
            pct = 100 * count / adata.n_obs
            f.write(f"   {cell_type}: {count:,} ({pct:.2f}%)\n")
        f.write("\n")
        f.write("3. Leiden聚类结果\n")
        f.write("-"*70 + "\n")
        f.write(f"默认分辨率: {DEFAULT_RESOLUTION}\n")
        f.write(f"Cluster数量: {adata.obs['leiden'].nunique()}\n\n")
        cluster_counts = adata.obs['leiden'].value_counts().sort_index()
        f.write("各cluster细胞数:\n")
        for cluster, count in cluster_counts.items():
            pct = 100 * count / adata.n_obs
            f.write(f"   Cluster {cluster}: {count:,} ({pct:.2f}%)\n")
        f.write("\n")
        f.write("4. Cluster细胞类型组成\n")
        f.write("-"*70 + "\n")
        ct_cluster = pd.crosstab(
            adata.obs[ANNOTATION_KEY],
            adata.obs['leiden']
        )
        for cluster in sorted(adata.obs['leiden'].unique()):
            f.write(f"\nCluster {cluster}:\n")
            cluster_data = ct_cluster[cluster]
            cluster_total = cluster_data.sum()
            for cell_type, count in cluster_data.items():
                if count > 0:
                    pct = 100 * count / cluster_total
                    f.write(f"   {cell_type}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
        if BATCH_KEY in adata.obs.columns:
            f.write("5. 批次分布\n")
            f.write("-"*70 + "\n")
            batch_counts = adata.obs[BATCH_KEY].value_counts()
            for batch, count in batch_counts.items():
                pct = 100 * count / adata.n_obs
                f.write(f"   {batch}: {count:,} ({pct:.2f}%)\n")
            f.write("\n")
        if TISSUE_KEY in adata.obs.columns:
            f.write("6. 组织来源分布\n")
            f.write("-"*70 + "\n")
            tissue_counts = adata.obs[TISSUE_KEY].value_counts()
            for tissue, count in tissue_counts.items():
                pct = 100 * count / adata.n_obs
                f.write(f"   {tissue}: {count:,} ({pct:.2f}%)\n")
            f.write("\n")
        score_cols = [col for col in adata.obs.columns if col.startswith('score_')]
        if score_cols:
            f.write("7. 功能评分统计\n")
            f.write("-"*70 + "\n")
            for score_col in score_cols:
                score_name = score_col.replace('score_', '')
                mean_score = adata.obs[score_col].mean()
                std_score = adata.obs[score_col].std()
                f.write(f"   {score_name}:\n")
                f.write(f"      Mean: {mean_score:.3f}\n")
                f.write(f"      Std: {std_score:.3f}\n")
            f.write("\n")
        f.write("8. 分析参数\n")
        f.write("-"*70 + "\n")
        f.write(f"高变基因数: {N_TOP_GENES}\n")
        f.write(f"PCA成分数: {N_PCS}\n")
        f.write(f"邻居方法: {adata.uns.get('integration_method','unknown')}\n")
        if adata.uns.get('integration_method','') == 'bbknn':
            f.write(f"BBKNN neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}\n")
            f.write(f"BBKNN metric: {BBKNN_METRIC}\n")
            f.write(f"BBKNN trim: {BBKNN_TRIM}\n")
        else:
            f.write(f"标准邻居 n_neighbors: {N_NEIGHBORS}\n")
        f.write(f"聚类分辨率: {LEIDEN_RESOLUTIONS}\n")
        f.write(f"默认分辨率: {DEFAULT_RESOLUTION}\n")
        f.write("\n")
        f.write("="*70 + "\n")
        f.write("报告结束\n")
        f.write("="*70 + "\n")
    log_msg(f"✓ 总结报告已保存: {report_file}")


def save_results(adata, output_dir):
    """
    保存分析结果
    """
    log_step(9, "Saving Results")
    
    output_path = Path(output_dir)
    
    # 1. 保存完整的AnnData对象
    h5ad_file = output_path / "adata_subtype_analyzed.h5ad"
    log_msg(f"\n保存AnnData对象: {h5ad_file}")
    adata.write_h5ad(h5ad_file, compression='gzip', compression_opts=9)
    log_msg("✓ AnnData已保存")
    
    # 2. 保存metadata
    metadata_file = output_path / "metadata.csv"
    log_msg(f"\n保存metadata: {metadata_file}")
    adata.obs.to_csv(metadata_file)
    log_msg("✓ Metadata已保存")
    
    # 3. 保存聚类结果
    clustering_file = output_path / "clustering_results.csv"
    log_msg(f"\n保存聚类结果: {clustering_file}")
    cluster_cols = [ANNOTATION_KEY] + [col for col in adata.obs.columns if 'leiden' in col.lower()]
    if BATCH_KEY in adata.obs.columns:
        cluster_cols.append(BATCH_KEY)
    if TISSUE_KEY in adata.obs.columns:
        cluster_cols.append(TISSUE_KEY)
    cluster_df = adata.obs[cluster_cols].copy()
    cluster_df.to_csv(clustering_file)
    log_msg("✓ 聚类结果已保存")
    
    # 4. 保存功能评分
    score_cols = [col for col in adata.obs.columns if col.startswith('score_')]
    if score_cols:
        scores_file = output_path / "functional_scores.csv"
        log_msg(f"\n保存功能评分: {scores_file}")
        score_df = adata.obs[[ANNOTATION_KEY, 'leiden'] + score_cols].copy()
        score_df.to_csv(scores_file)
        log_msg("✓ 功能评分已保存")
    
    log_msg(f"\n✓ 所有结果已保存到: {output_dir}")


# ==================== 主程序 ====================

def main():
    """
    主函数
    """
    log_msg("\n" + "="*70)
    log_msg("上皮细胞亚型深入分析流程")
    log_msg("="*70)
    log_msg(f"\n目标细胞类型: {TARGET_CELL_TYPES}")
    log_msg(f"输入文件: {INPUT_H5AD_PATH}")
    log_msg(f"输出目录: {OUTPUT_DIR}")
    
    # 设置输出目录
    output_path = setup_output_dir(OUTPUT_DIR, OVERWRITE_EXISTING)
    
    # 记录开始时间
    start_time = time.time()
    
    try:
        # Step 1: 加载数据并提取子集
        adata = load_and_subset_data(
            INPUT_H5AD_PATH,
            TARGET_CELL_TYPES,
            ANNOTATION_KEY
        )
        
        # Step 1b: 绘制初始UMAP (使用原始坐标)
        plot_initial_umap(adata, OUTPUT_DIR)
        
        # Step 2: 预处理
        adata = preprocess_subset(adata, OUTPUT_DIR)
        
        # Step 3: 降维 + BBKNN/Neighbors + UMAP
        adata = run_dimensionality_reduction(adata)
        
        # Step 4: 聚类
        adata = run_clustering(adata)
        
        # Step 5: 功能评分
        adata = calculate_functional_scores(adata)
        
        # Step 6: 查找marker基因
        adata = find_cluster_markers(adata, OUTPUT_DIR)
        
        # Step 7: 生成可视化
        generate_visualizations(adata, OUTPUT_DIR)
        
        # Step 8: 生成总结报告
        generate_summary_report(adata, OUTPUT_DIR, TARGET_CELL_TYPES)
        
        # Step 9: 保存结果
        save_results(adata, OUTPUT_DIR)
        
        # 计算总耗时
        elapsed = time.time() - start_time
        log_msg("\n" + "="*70)
        log_msg(f"✓ 分析完成!")
        log_msg(f"   总耗时: {elapsed/60:.2f} 分钟")
        log_msg(f"   结果保存在: {OUTPUT_DIR}")
        log_msg("="*70)
        
    except Exception as e:
        log_msg(f"\n❌ 分析失败: {e}")
        import traceback
        log_msg(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
