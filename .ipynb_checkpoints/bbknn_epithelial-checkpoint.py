#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
上皮细胞 BBKNN 批次整合与分析流程 - 统一优化版

特点：
1. 完整的 BBKNN 批次校正流程，并提供事后参数校验信息
2. 强制保留指定 marker（避免在 HVG/过滤中被丢弃）
3. 差异分析/绘图默认从 .raw 的全基因矩阵读取（use_raw=True）
4. 多分辨率 Leiden 聚类；UMAP 使用 BBKNN 构建的邻居图
5. 产出详细且可复现的可视化与报告

作者：临床-生信团队
版本：v1.3
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
import time
from tqdm import tqdm
import pickle
warnings.filterwarnings('ignore')

# ==================== 配置区 ====================

# ---------- 输入/输出 ----------
INPUT_H5AD_PATH = "/home/h2048/data/py/1029/bbknn_celltype_analysis/Epithelial/adata_Epithelial_bbknn.h5ad"  # 输入 h5ad 文件路径
OUTPUT_DIR = "/home/h2048/data/py/1029/bbknn_celltype_analysis/Epithelial/output_optimized"  # 输出目录
OVERWRITE_EXISTING = True  # 若已存在产物是否覆盖

# ---------- BBKNN 整合参数 ----------
BATCH_KEY = "dataset"                 # 批次信息列名（位于 adata.obs）
BBKNN_NEIGHBORS_WITHIN_BATCH = 3      # 每个批次内选取的邻居数
BBKNN_N_PCS = 50                      # 参与 BBKNN 的 PC 数

# ---------- 高变基因（HVG） ----------
USE_HVG = True
N_TOP_GENES = 3000
HVG_FLAVOR = "seurat_v3"

# ---------- 基因过滤 ----------
MIN_CELLS_PER_GENE = 3  # 至少在 N 个细胞中出现

# ---------- 原始 counts 来源 ----------
# 可选："auto" | "raw" | "layers" | "X"
RAW_COUNTS_SOURCE = "auto"

# ---------- 标准化 ----------
NORMALIZE_TOTAL = True
TARGET_SUM = 1e4
LOG_TRANSFORM = True
SCALE_DATA = True
MAX_VALUE = 10

# ---------- PCA ----------
N_PCS = 50

# ---------- 降维/可视化 ----------
RUN_UMAP = True
UMAP_MIN_DIST = 0.3
# UMAP 在 BBKNN 模式下不需要显式设置 n_neighbors

# ---------- 聚类 ----------
RUN_CLUSTERING = True
LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
DEFAULT_RESOLUTION = 2.0

# ---------- 上皮细胞 Marker 基因 ----------
MARKER_EPITHELIAL = [
    # Basal
    'TP63','KRT5','KRT14',
    # Secretory
    'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1',
    # Ciliated
    'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24',
    # Differentiation
    'KRT19','NOTCH3','KRT16','KRT23',
    # Goblet
    'MUC5AC','SPDEF','LYPD2','ITLN1',
    # Neuroendocrine
    'ASCL1','GRP','KRT8',
    # Tuft / Ionocytes
    'POU2F3','ASCL2','CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB','PDE1C',
    # AT1
    'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
    # AT2
    'SFTPC','LAMP3','MFSD2A','C8orf4','C11orf96','SFTPB','SFTA2',
    # Mesenchymal
    'VIM','SOX9','MYH11','ACTA2','MYLK',
    # Secretory/Immune
    'DMBT1','RNASE1','MUC5B','LYZ','LTF','PIP','CCL28',
    # Ciliogenesis
    'DEUP1','FOXN4','CDC20B','CCNO',
    # Proliferation
    'MKI67','TOP2A','TK1','CENPW'
]

# 关键 marker（用于分型/展示）
KEY_MARKERS = {
    'Basal': ['TP63','KRT5'],
    'Secretory': ['SCGB1A1','SCGB3A2'],
    'Ciliated': ['FOXJ1','RSPH1'],
    'Goblet': ['MUC5AC','SPDEF'],
    'Tuft': ['POU2F3','ASCL2'],
    'Ionocyte': ['FOXI1','CFTR'],
    'AT1': ['AGER','RTKN2'],
    'AT2': ['SFTPC','SFTPB'],
    'Proliferating': ['MKI67','TOP2A']
}

# ---------- Marker 基因分析 ----------
RUN_FIND_MARKERS = True
MARKER_MIN_PCT = 0.25
MARKER_LOGFC_THRESHOLD = 0.25
TOP_N_MARKERS = 10

# ---------- 性能/缓存 ----------
USE_CACHE = True
CHUNK_SIZE = 100

# ---------- 可视化 ----------
GENERATE_DOTPLOT = True
GENERATE_HEATMAP = True
GENERATE_FACET_PLOTS = True
DPI = 300
FIGURE_FORMAT = "png"

VERBOSE = True

# ---------- 强制保留（即便不在 HVG/或不满足 min_cells 也不丢弃） ----------
FORCE_INCLUDE_MARKERS = sorted({
    *MARKER_EPITHELIAL,
    *[g for genes in KEY_MARKERS.values() for g in genes],
})

# ==================== 工具函数 ====================

def log_msg(msg):
    """轻量日志打印（受 VERBOSE 控制）。"""
    if VERBOSE:
        print(msg)

def log_step(step_num, step_name):
    """统一格式打印阶段标题。"""
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)

def save_checkpoint(data, filename, output_dir):
    """保存中间结果到缓存文件（启用 USE_CACHE 时生效）。"""
    if USE_CACHE:
        checkpoint_path = Path(output_dir) / f".cache_{filename}"
        with open(checkpoint_path, 'wb') as f:
            pickle.dump(data, f)
        log_msg(f"   Checkpoint saved: {checkpoint_path}")

def load_checkpoint(filename, output_dir):
    """读取缓存文件（若存在且启用 USE_CACHE）。"""
    if USE_CACHE:
        checkpoint_path = Path(output_dir) / f".cache_{filename}"
        if checkpoint_path.exists():
            with open(checkpoint_path, 'rb') as f:
                return pickle.load(f)
    return None

# ==================== 主功能函数 ====================


def load_anndata(h5ad_path):
    """读取 AnnData，并检查批次列与分布。"""
    import scanpy as sc
    log_step(1, "加载数据")
    log_msg(f"\n读取文件: {h5ad_path}")
    if not Path(h5ad_path).exists():
        raise FileNotFoundError(f"File not found: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    log_msg(f"数据读取完成")
    log_msg(f"   细胞数: {adata.n_obs:,}")
    log_msg(f"   基因数: {adata.n_vars:,}")

    log_msg(f"\n可用的元数据列:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        log_msg(f"   - {col}: {n_unique} unique values")

    if BATCH_KEY not in adata.obs.columns:
        raise ValueError(f"批次列 '{BATCH_KEY}' 未在 adata.obs 中找到")

    log_msg(f"\n批次分布 (key: {BATCH_KEY}):")
    batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")
    return adata


def check_marker_genes(adata):
    """检查 marker 可用性（同时检查 var 与 raw）。"""
    log_step(2, "检查 Marker 基因")
    raw_names = set(adata.raw.var_names) if adata.raw is not None else set()
    var_names = set(adata.var_names)

    def present(g):
        return (g in var_names) or (g in raw_names)

    available_markers = [g for g in MARKER_EPITHELIAL if present(g)]
    missing_markers   = [g for g in MARKER_EPITHELIAL if not present(g)]

    log_msg(f"\nMarker 基因可用性:")
    log_msg(f"   总数: {len(MARKER_EPITHELIAL)}")
    log_msg(f"   可用: {len(available_markers)} ({len(available_markers)/len(MARKER_EPITHELIAL)*100:.1f}%)")
    log_msg(f"   缺失: {len(missing_markers)} ({len(missing_markers)/len(MARKER_EPITHELIAL)*100:.1f}%)")

    if missing_markers:
        log_msg(f"\n   缺失基因示例: {', '.join(missing_markers[:10])}")
        if len(missing_markers) > 10:
            log_msg(f"   ... 以及另外 {len(missing_markers)-10} 个")

    log_msg(f"\n关键细胞类型 marker:")
    for celltype, genes in KEY_MARKERS.items():
        available = [g for g in genes if present(g)]
        log_msg(f"   {celltype}: {len(available)}/{len(genes)} 可用")
    return available_markers

import numpy as np
from scipy.sparse import issparse, csr_matrix
import scanpy as sc


def preprocess_for_bbknn(adata, n_hvg=3000, batch_key=None, FORCE_INCLUDE_MARKERS=None):
    """
    目的：
      - 将原始 counts 保存到 layers['counts']
      - 归一化 + log1p 后的表达保存到 layers['log1p']
      - 返回 (adata_full, adata_hvg)

    说明：
      - adata_full：完整基因矩阵，X=log1p-normalized，layers['counts']=原始计数，layers['log1p']=与 X 一致
      - adata_hvg：HVG（+强制保留）子集，X=log1p-normalized，layers['counts']=子集对应原始计数
      - adata_hvg.raw = adata_full：在 .raw 中保留完整基因，用于后续差异分析/可视化（use_raw=True）
    """
    # —— 轻量日志 ——
    def _log(msg):
        try:
            log_msg(msg)
        except NameError:
            print(msg)

    _log("="*70)
    _log(f"筛选高变基因 (n={n_hvg}) 并设置 layers，不直接写入 .raw")
    _log("-"*70)

    # 1) 确保原始 counts 在 layers['counts']
    # 优先使用已有 counts；否则从 .raw 或 X 提取
    if 'counts' in adata.layers:
        counts = adata.layers['counts']
        _log("   使用已有 layers['counts'] 作为原始计数")
    elif getattr(adata, "raw", None) is not None:
        counts = adata.raw.X
        _log("   从 .raw 提取原始计数 -> layers['counts']")
    else:
        counts = adata.X
        _log("   未检测到 .raw，使用当前 X 作为原始计数 -> layers['counts']")

    # 稀疏化以节省内存
    if not issparse(counts):
        counts = csr_matrix(counts)
    adata.layers['counts'] = counts

    # 2) HVG 选择（带稳健回退）
    # 为保证 HVG 计算稳定，scanpy 默认用 X；此处先令 X 指向原始 counts 的副本
    adata.X = adata.layers['counts'].copy()

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvg,
            flavor="seurat_v3",
            batch_key=batch_key
        )
        _log(f"   HVG(seurat_v3, batch_key={batch_key}) 成功")
    except Exception as e:
        _log(f"   ⚠️  HVG(seurat_v3, batch_key={batch_key}) 失败: {repr(e)}")
        _log("      → 回退至不带 batch_key 的 seurat_v3")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvg,
            flavor="seurat_v3",
            batch_key=None
        )

    hv_mask = adata.var['highly_variable'].to_numpy(dtype=bool)

    # 强制保留 marker
    if FORCE_INCLUDE_MARKERS is None:
        FORCE_INCLUDE_MARKERS = []
    force_mask = adata.var_names.isin(FORCE_INCLUDE_MARKERS)

    keep_mask = hv_mask | force_mask
    n_force = int(force_mask.sum())
    _log(f"   选中 HVG: {int(hv_mask.sum())}; 强制保留: {n_force}; 合计: {int(keep_mask.sum())}")

    # 3) 在“完整矩阵”上做一次标准化，得到 adata_full（X=log1p）
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

    # 4) 构建 HVG 子集对象（不直接覆写 .raw）
    adata_hvg = adata[:, keep_mask].copy()
    # 子集的 counts 层
    adata_hvg.layers['counts'] = adata.layers['counts'][:, keep_mask].copy()

    # 令 HVG 子集的 X = 子集 counts 的 log1p 结果
    adata_hvg.X = adata_hvg.layers['counts'].copy()
    sc.pp.normalize_total(adata_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hvg)
    adata_hvg.layers['log1p'] = adata_hvg.X.copy()

    # 5) 基本元信息记录
    adata_hvg.uns['hvg_n_top_genes'] = int(n_hvg)
    adata_hvg.uns['force_include_markers'] = np.array(FORCE_INCLUDE_MARKERS, dtype=str)
    adata_hvg.uns['keep_mask_sum'] = int(keep_mask.sum())

    # 6) 将完整基因数据设为 HVG 子集的 .raw（后续 use_raw=True 可访问）
    adata_hvg.raw = adata
    _log(f"   设置 adata_hvg.raw = adata_full（包含 {adata.n_vars} 个基因）")

    _log("   ✓ 完成：原始计数存于 layers['counts']，log1p 存于 layers['log1p']")
    _log("="*70)

    # 返回：完整矩阵 + HVG 子集
    return adata, adata_hvg


def run_pca(adata, n_pcs=50):
    """
    执行 PCA 降维。

    参数：
    adata : AnnData
        数据对象
    n_pcs : int
        计算的主成分数

    返回：
    adata : AnnData
        在 .obsm['X_pca'] 中附加 PCA 结果
    """
    import scanpy as sc
    log_msg(f"\n运行 PCA (n_comps={n_pcs})...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    var_ratio = adata.uns['pca']['variance_ratio']
    cumsum_var = np.cumsum(var_ratio)
    if len(cumsum_var) >= 50:
        log_msg(f"   PC1-10 解释率: {cumsum_var[9]:.2%}")
        log_msg(f"   PC1-20 解释率: {cumsum_var[19]:.2%}")
        log_msg(f"   PC1-50 解释率: {cumsum_var[49]:.2%}")
    return adata


def _bbknn_sanity_check(adata):
    """打印 BBKNN 邻居数的事后校验信息。"""
    try:
        n_batches = adata.obs[BATCH_KEY].nunique()
        expected = BBKNN_NEIGHBORS_WITHIN_BATCH * n_batches
        neigh = adata.uns.get('neighbors', {}).get('params', {})
        actual = neigh.get('n_neighbors', None)
        log_msg(
            f"   BBKNN 检查 → n_batches={n_batches}, "
            f"neighbors_within_batch={BBKNN_NEIGHBORS_WITHIN_BATCH} → "
            f"期望总邻居≈{expected}, 实际 {actual}"
        )
    except Exception as e:
        log_msg(f"   ⚠️ BBKNN 事后校验失败: {e}")


def run_bbknn_integration(adata, batch_key='dataset', neighbors_within_batch=3, n_pcs=50):
    """
    执行 BBKNN 批次整合。

    参数：
    adata : AnnData
        已完成 PCA 的数据对象
    batch_key : str
        批次信息列名
    neighbors_within_batch : int
        每个批次取的邻居数
    n_pcs : int
        BBKNN 使用的 PC 数

    返回：
    adata : AnnData
        在 .uns/.obsp 中写入邻居图
    """
    from bbknn import bbknn as bbknn_func
    log_step(4, "BBKNN 批次整合")
    log_msg("\nBBKNN 参数:")
    log_msg(f"   batch_key: {batch_key}")
    log_msg(f"   neighbors_within_batch: {neighbors_within_batch}")
    log_msg(f"   n_pcs: {n_pcs}")

    log_msg("\n运行 BBKNN...")
    start_time = time.time()
    bbknn_func(
        adata,
        batch_key=batch_key,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        copy=False
    )
    elapsed_time = time.time() - start_time
    log_msg(f"\n   BBKNN 完成，用时 {elapsed_time:.1f} 秒")

    _bbknn_sanity_check(adata)
    return adata


def run_umap_and_clustering(adata, run_umap=True, min_dist=0.3,
                            run_clustering=True, resolutions=[1.2, 1.6, 2.0, 2.4, 2.8],
                            default_resolution=2.0):
    """
    使用 BBKNN 邻居图运行 UMAP，并进行多分辨率 Leiden 聚类。

    参数：
    adata : AnnData
        已计算邻居图的数据对象
    run_umap : bool
        是否运行 UMAP
    min_dist : float
        UMAP 的 min_dist 参数
    run_clustering : bool
        是否运行 Leiden 聚类
    resolutions : list[float]
        需要测试的分辨率列表
    default_resolution : float
        默认分辨率，结果写入 adata.obs['leiden_bbknn']

    返回：
    adata : AnnData
        写入 UMAP 坐标与多分辨率聚类结果
    """
    import scanpy as sc
    log_step(5, "UMAP 与多分辨率聚类")

    if run_umap:
        log_msg(f"\n运行 UMAP (min_dist={min_dist}，使用 BBKNN 邻居图)...")
        sc.tl.umap(adata, min_dist=min_dist)
        log_msg("   UMAP 完成")

    if run_clustering:
        log_msg(f"\n运行多分辨率 Leiden 聚类...")
        log_msg(f"   分辨率列表: {resolutions}")
        for res in resolutions:
            cluster_key = f'leiden_bbknn_res{res}'
            try:
                sc.tl.leiden(adata, resolution=res, key_added=cluster_key, flavor='igraph',
                             n_iterations=2, directed=False)
            except TypeError:
                sc.tl.leiden(adata, resolution=res, key_added=cluster_key, n_iterations=2)
            n_clusters = adata.obs[cluster_key].nunique()
            log_msg(f"   分辨率 {res}: {n_clusters} 个簇")
        default_key = f'leiden_bbknn_res{default_resolution}'
        adata.obs['leiden_bbknn'] = adata.obs[default_key]
        log_msg(f"\n   默认聚类: {default_key}")
    return adata


def compute_pct_expressed_vectorized(adata, cluster_key, cluster, genes):
    """向量化计算基因在簇内/外的表达比例；优先使用 var，不在则回退至 raw。"""
    cluster_mask = (adata.obs[cluster_key] == cluster).values
    n_in = int(cluster_mask.sum()); n_out = int((~cluster_mask).sum())
    if n_in == 0 or n_out == 0 or len(genes) == 0:
        return [], []

    genes_in_var = [g for g in genes if g in adata.var_names]
    rem = [g for g in genes if g not in adata.var_names]
    X_in = X_out = None

    if genes_in_var:
        Xin = adata[cluster_mask, genes_in_var].X
        Xout = adata[~cluster_mask, genes_in_var].X
        Xin = Xin.toarray() if hasattr(Xin, 'toarray') else np.asarray(Xin)
        Xout = Xout.toarray() if hasattr(Xout, 'toarray') else np.asarray(Xout)
        X_in, X_out = Xin, Xout

    if rem and (adata.raw is not None):
        # 构建 raw 名称→索引映射，避免多次 list.index
        raw_names = np.array(adata.raw.var_names)
        raw_map = {g:i for i,g in enumerate(raw_names)}
        rem_exist = [g for g in rem if g in raw_map]
        if rem_exist:
            ridx = [raw_map[g] for g in rem_exist]
            R = adata.raw.X
            Rin = R[cluster_mask][:, ridx]
            Rout = R[~cluster_mask][:, ridx]
            Rin = Rin.toarray() if hasattr(Rin, 'toarray') else np.asarray(Rin)
            Rout = Rout.toarray() if hasattr(Rout, 'toarray') else np.asarray(Rout)
            if X_in is None:
                X_in, X_out = Rin, Rout
                genes_in_var = rem_exist
            else:
                X_in = np.concatenate([X_in, Rin], axis=1)
                X_out = np.concatenate([X_out, Rout], axis=1)
                genes_in_var = genes_in_var + rem_exist

    if X_in is None:
        return [], []
    pct_in  = (X_in  > 0).sum(axis=0) / n_in
    pct_out = (X_out > 0).sum(axis=0) / n_out
    return pct_in.tolist(), pct_out.tolist()


def find_all_markers_optimized(adata, cluster_key, output_dir,
                               min_pct=0.25,
                               logfc_threshold=0.25):
    """
    优化版 FindAllMarkers：
      - 结果为空时安全返回
      - 控制每簇导出数量上限（TOP_N_MARKERS）
      - 差异分析使用 .raw 全基因矩阵（use_raw=True）
    """
    import scanpy as sc
    log_step(6, "计算聚类差异基因（优化版）")

    cache_file = "epithelial_markers_cache.pkl"
    cached_markers = load_checkpoint(cache_file, output_dir)
    if cached_markers is not None:
        log_msg("\n   从缓存加载...")
        log_msg(f"   Markers 总数: {len(cached_markers)}")
        return cached_markers

    markers_file = Path(output_dir) / "cluster_markers.csv"

    log_msg(f"\n运行差异分析（Wilcoxon）...")
    log_msg(f"   min_pct: {min_pct}")
    log_msg(f"   logfc_threshold: {logfc_threshold}")

    # 差异分析（显式使用 .raw）
    log_msg("\n   计算差异...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method='wilcoxon',
        key_added='rank_genes_groups',
        use_raw=True,
        layer=None
    )
    log_msg("   ✓ 差异分析完成")

    log_msg("\n   提取并过滤差异结果...")
    clusters = adata.obs[cluster_key].unique()
    markers_list = []

    for cluster in tqdm(clusters, desc="Processing clusters"):
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges","pvals_adj"])
        df = df[(df["logfoldchanges"] > logfc_threshold) & (df["pvals_adj"] < 0.05)].copy()
        if df.empty:
            continue

        genes = df["names"].tolist()
        pct_in, pct_out = compute_pct_expressed_vectorized(adata, cluster_key, cluster, genes)
        if len(pct_in) != len(genes):
            n = len(genes)
            pct_in  = (pct_in  + [np.nan]*(n-len(pct_in)))[:n]
            pct_out = (pct_out + [np.nan]*(n-len(pct_out)))[:n]

        df["cluster"] = cluster
        df["pct_in_cluster"] = pct_in
        df["pct_out_cluster"] = pct_out
        df = df[df["pct_in_cluster"] > min_pct]
        if df.empty:
            continue

        if TOP_N_MARKERS is not None and TOP_N_MARKERS > 0:
            df = df.sort_values("pvals_adj", ascending=True).head(TOP_N_MARKERS)
        markers_list.append(df)

    if len(markers_list) == 0:
        log_msg("\n   ⚠️ 无符合当前阈值的差异基因，返回空结果。")
        empty_cols = ["names","scores","logfoldchanges","pvals","pvals_adj","cluster","pct_in_cluster","pct_out_cluster"]
        empty_df = pd.DataFrame(columns=empty_cols)
        empty_df.to_csv(markers_file, index=False)
        save_checkpoint(empty_df, cache_file, output_dir)
        return empty_df

    all_markers = pd.concat(markers_list, ignore_index=True)
    log_msg(f"\n   ✓ 共获得 {len(all_markers)} 个 marker，覆盖 {len(clusters)} 个簇")
    save_checkpoint(all_markers, cache_file, output_dir)
    all_markers.to_csv(markers_file, index=False)
    log_msg(f"   ✓ 已保存: {markers_file}")
    return all_markers


def generate_visualizations(adata, available_markers, output_dir):
    """生成 UMAP 概览、多分辨率聚类、marker UMAP、Dotplot、Heatmap、按批次分面等图。"""
    import scanpy as sc
    log_step(7, "生成可视化图件")

    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(exist_ok=True)
    sc.settings.figdir = fig_dir

    # 1) 概览图
    log_msg("\n   生成整合质量概览图...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0], show=False, title='Batch')
    sc.pl.umap(adata, color='leiden_bbknn', ax=axes[1], show=False,
               title='Clusters', legend_loc='on data', legend_fontsize=8)
    if 'cell_type' in adata.obs.columns:
        sc.pl.umap(adata, color='cell_type', ax=axes[2], show=False, title='Cell Type')
    else:
        axes[2].axis('off')
    plt.tight_layout()
    plt.savefig(fig_dir / f'umap_overview.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    log_msg(f"   ✓ 概览 UMAP 已保存")

    # 2) 多分辨率聚类
    if RUN_CLUSTERING and len(LEIDEN_RESOLUTIONS) > 1:
        log_msg("\n   生成多分辨率聚类图...")
        n_res = len(LEIDEN_RESOLUTIONS)
        n_cols = 3
        n_rows = (n_res + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 6*n_rows))
        axes = axes.flatten() if n_res > 1 else [axes]
        for i, res in enumerate(LEIDEN_RESOLUTIONS):
            cluster_key = f'leiden_bbknn_res{res}'
            sc.pl.umap(
                adata, color=cluster_key, ax=axes[i], show=False,
                title=f'Resolution {res}', legend_loc='on data', legend_fontsize=8
            )
        for i in range(n_res, len(axes)):
            axes[i].axis('off')
        plt.tight_layout()
        plt.savefig(fig_dir / f'umap_resolutions.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
        plt.close()
        log_msg(f"   ✓ 多分辨率 UMAP 已保存")

    # 3) 关键 marker 表达 UMAP（use_raw=True）
    log_msg("\n   生成 marker 表达 UMAP...")
    for celltype, genes in KEY_MARKERS.items():
        available = [g for g in genes if (g in (adata.raw.var_names if adata.raw is not None else adata.var_names))]
        if available:
            sc.pl.umap(
                adata,
                color=available,
                use_raw=True,
                cmap='Reds',
                ncols=len(available),
                vmax='p99',
                save=f'_{celltype}_markers.{FIGURE_FORMAT}'
            )
            log_msg(f"   ✓ {celltype} markers UMAP 已保存")

    # 4) Dotplot（use_raw=True）
    if GENERATE_DOTPLOT and len(available_markers) > 0:
        log_msg("\n   生成 marker Dotplot...")
        try:
            markers_to_plot = available_markers[:50] if len(available_markers) > 50 else available_markers
            sc.pl.dotplot(
                adata,
                var_names=markers_to_plot,
                groupby='leiden_bbknn',
                standard_scale='var',
                use_raw=True,
                save=f'_epithelial_markers.{FIGURE_FORMAT}',
                figsize=(max(20, len(markers_to_plot)*0.4), 10)
            )
            log_msg(f"   ✓ Dotplot 已保存（展示 {len(markers_to_plot)} 个基因）")
        except Exception as e:
            log_msg(f"   ⚠️ Dotplot 生成失败: {e}")

    # 5) Heatmap（按簇均值表达；优先 raw）
    if GENERATE_HEATMAP and len(available_markers) > 0:
        log_msg("\n   生成表达热图...")
        try:
            cluster_expr = pd.DataFrame()
            for gene in available_markers:
                if (adata.raw is not None) and (gene in adata.raw.var_names):
                    Xg = adata.raw[:, gene].X
                elif gene in adata.var_names:
                    Xg = adata[:, gene].X
                else:
                    continue
                Xg = Xg.toarray().ravel() if hasattr(Xg, 'toarray') else np.asarray(Xg).ravel()
                cluster_expr[gene] = Xg
            cluster_expr['leiden_bbknn'] = adata.obs['leiden_bbknn'].values
            cluster_mean_expr = cluster_expr.groupby('leiden_bbknn').mean()
            plt.figure(figsize=(max(20, len(available_markers)*0.3), 10))
            sns.heatmap(
                cluster_mean_expr.T,
                cmap='RdYlBu_r',
                center=0,
                robust=True,
                yticklabels=True,
                xticklabels=True,
                cbar_kws={'label': 'Mean Expression'}
            )
            plt.title('Epithelial Marker Expression Across Clusters')
            plt.xlabel('Cluster')
            plt.ylabel('Marker Genes')
            plt.tight_layout()
            plt.savefig(fig_dir / f'heatmap_marker_expression.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
            plt.close()
            log_msg(f"   ✓ Heatmap 已保存")
        except Exception as e:
            log_msg(f"   ⚠️ Heatmap 生成失败: {e}")

    # 6) 按批次分面
    if GENERATE_FACET_PLOTS and BATCH_KEY in adata.obs.columns:
        log_msg("\n   生成分批次 UMAP...")
        try:
            batches = sorted(adata.obs[BATCH_KEY].unique())
            n_batches = len(batches)
            n_cols = min(3, n_batches)
            n_rows = (n_batches + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
            axes = [axes] if n_batches == 1 else axes.flatten()
            
            for i, batch in enumerate(batches):
                # 只选择该 batch 的细胞
                adata_batch = adata[adata.obs[BATCH_KEY] == batch].copy()
                
                sc.pl.umap(
                    adata_batch,
                    color='leiden_bbknn',
                    ax=axes[i],
                    show=False,
                    title=f'Batch: {batch}',
                    legend_loc='right margin',
                    legend_fontsize=6
                )
            
            for i in range(n_batches, len(axes)):
                axes[i].axis('off')
            
            plt.tight_layout()
            plt.savefig(fig_dir / f'umap_by_batch.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
            plt.close()
            log_msg(f"   ✓ 分批次 UMAP 已保存")
        except Exception as e:
            log_msg(f"   ⚠️ 分批次 UMAP 生成失败: {e}")

    log_msg(f"\n   所有图件已保存至: {fig_dir}")


def generate_summary_report(adata, available_markers, output_dir, processing_time=None):
    """生成纯文本总结报告。"""
    from datetime import datetime
    log_step(8, "生成总结报告")
    report = []
    report.append("="*70)
    report.append("Epithelial Cell Analysis - BBKNN Integration Summary")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if processing_time:
        report.append(f"Processing time: {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    report.append("")
    report.append("[ Data Overview ]")
    report.append(f"  Cells: {adata.n_obs:,}")
    report.append(f"  Genes: {adata.n_vars:,}")
    report.append("")
    if BATCH_KEY in adata.obs.columns:
        report.append(f"[ Batch Distribution (key: {BATCH_KEY}) ]")
        batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
        for batch, count in batch_counts.items():
            pct = count / adata.n_obs * 100
            report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
        report.append("")
    report.append("[ Epithelial Marker Genes ]")
    report.append(f"  Total markers: {len(MARKER_EPITHELIAL)}")
    report.append(f"  Available: {len(available_markers)} ({len(available_markers)/len(MARKER_EPITHELIAL)*100:.1f}%)")
    report.append("")
    if RUN_CLUSTERING:
        report.append("[ Multi-resolution Clustering Results ]")
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_bbknn_res{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                default_marker = " (default)" if res == DEFAULT_RESOLUTION else ""
                report.append(f"  Resolution {res}: {n_clusters} clusters{default_marker}")
        report.append("")
    report.append("[ BBKNN Configuration ]")
    report.append(f"  neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    report.append(f"  n_pcs: {BBKNN_N_PCS}")
    report.append(f"  batch_key: {BATCH_KEY}")
    report.append(f"  High variable genes: {USE_HVG} (n={N_TOP_GENES if USE_HVG else 'N/A'})")
    report.append("")
    report.append("[ Output Files ]")
    report.append(f"  - Data: {Path(output_dir) / 'epithelial_bbknn_integrated.h5ad'}")
    report.append(f"  - Figures: {Path(output_dir) / 'figures/'}*.{FIGURE_FORMAT}")
    if RUN_FIND_MARKERS:
        report.append(f"  - Markers: {Path(output_dir) / 'cluster_markers.csv'}")
    report.append("")
    report.append("="*70)
    report.append("Analysis completed successfully")
    report.append("="*70)
    report_text = '\n'.join(report)
    report_path = Path(output_dir) / "analysis_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    log_msg(f"\nReport saved to: {report_path}")
    log_msg("\n" + report_text)

# ==================== 主程序 ====================

def main():
    """主入口：按顺序执行加载→检查→预处理→PCA→BBKNN→UMAP/聚类→差异→可视化→报告。"""
    print("\n" + "="*70)
    print("Epithelial Cell BBKNN Integration Pipeline (Optimized v1.3)")
    print("="*70)

    start_time = time.time()

    log_msg("\n检查 Python 环境...")
    try:
        import scanpy as sc   # noqa
        from bbknn import bbknn as bbknn_func  # noqa
        log_msg(f"   scanpy: {sc.__version__}")
        log_msg(f"   bbknn: installed")
    except ImportError as e:
        print(f"\nError: {e}", file=sys.stderr)
        print("\nInstallation: pip install scanpy bbknn", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_msg(f"\n输出目录: {output_dir}")

    # Step 1 加载
    adata = load_anndata(INPUT_H5AD_PATH)

    # Step 2 Marker 检查
    available_markers = check_marker_genes(adata)

    # Step 3 预处理（不直接覆写 .raw，完整矩阵保存在返回的 adata_full）
    adata_full, adata_hvg = preprocess_for_bbknn(
        adata,
        n_hvg=N_TOP_GENES,
        batch_key=BATCH_KEY,
        FORCE_INCLUDE_MARKERS=FORCE_INCLUDE_MARKERS
    )

    # Step 3.5 PCA
    adata_hvg = run_pca(adata_hvg, n_pcs=N_PCS)

    # Step 4 BBKNN 整合
    adata_hvg = run_bbknn_integration(
        adata_hvg,
        batch_key=BATCH_KEY,
        neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
        n_pcs=BBKNN_N_PCS
    )

    # Step 5 UMAP 与聚类
    if RUN_UMAP or RUN_CLUSTERING:
        adata_hvg = run_umap_and_clustering(
            adata_hvg,
            run_umap=RUN_UMAP,
            min_dist=UMAP_MIN_DIST,
            run_clustering=RUN_CLUSTERING,
            resolutions=LEIDEN_RESOLUTIONS,
            default_resolution=DEFAULT_RESOLUTION
        )

    # Step 6 差异基因
    all_markers = None
    if RUN_FIND_MARKERS:
        all_markers = find_all_markers_optimized(
            adata_hvg,
            'leiden_bbknn',
            output_dir,
            min_pct=MARKER_MIN_PCT,
            logfc_threshold=MARKER_LOGFC_THRESHOLD
        )

    # Step 7 可视化
    if RUN_UMAP:
        generate_visualizations(adata_hvg, available_markers, output_dir)

    # 结果回填至完整对象（以便统一保存/下游使用）
    log_msg("\n回填结果到完整数据集...")
    for key in ['X_pca','X_umap']:
        if key in adata_hvg.obsm:
            adata_full.obsm[key] = adata_hvg.obsm[key]
    for col in adata_hvg.obs.columns:
        if col.startswith('leiden_bbknn'):
            adata_full.obs[col] = adata_hvg.obs[col]
    if 'neighbors' in adata_hvg.uns:
        adata_full.uns['neighbors'] = adata_hvg.uns['neighbors']
    if 'connectivities' in adata_hvg.obsp:
        adata_full.obsp['connectivities'] = adata_hvg.obsp['connectivities']
    if 'distances' in adata_hvg.obsp:
        adata_full.obsp['distances'] = adata_hvg.obsp['distances']

    # 保存整合结果
    log_msg("\n保存整合后的 AnnData...")
    final_path = output_dir / "epithelial_bbknn_integrated.h5ad"
    log_msg("   使用 gzip 压缩...")
    adata_full.write_h5ad(final_path, compression='gzip', compression_opts=9)
    file_size = final_path.stat().st_size / (1024**3)
    log_msg(f"   Saved: {final_path} ({file_size:.2f} GB)")

    # 报告
    total_time = time.time() - start_time
    generate_summary_report(adata_full, available_markers, output_dir, total_time)

    # 终端摘要
    print("\n" + "="*70)
    print("✓ All analyses completed successfully")
    print("="*70)
    print(f"\nOutput directory: {output_dir}")
    print(f"Data: {final_path} ({file_size:.2f} GB)")
    print(f"Figures: {output_dir / 'figures/'}*.{FIGURE_FORMAT}")
    if RUN_FIND_MARKERS:
        print(f"Markers: {output_dir / 'cluster_markers.csv'}")
    print(f"\nTotal time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")

    print(f"\n📊 关键信息:")
    print(f"   Cells: {adata_full.n_obs:,}")
    print(f"   Available markers: {len(available_markers)}/{len(MARKER_EPITHELIAL)}")
    print(f"   Default clusters: {adata_full.obs['leiden_bbknn'].nunique()}")

    print(f"\n🔧 BBKNN 参数:")
    print(f"   neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    print(f"   n_pcs: {BBKNN_N_PCS}")
    print(f"   HVG: {N_TOP_GENES}")
    print()

    return adata_full, output_dir


if __name__ == "__main__":
    try:
        adata, output_dir = main()
    except KeyboardInterrupt:
        print("\n\n⚠️  用户终止", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n✖ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
