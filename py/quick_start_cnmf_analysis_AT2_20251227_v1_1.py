#!/usr/bin/env python3
"""
cNMF Downstream Analysis - Quick Start Template (IMPROVED v1.1)
================================================================

改进点：
1. 使用 safe_name() 确保与生产管线一致
2. 增强文件模式匹配（spectra.*consensus）
3. GEP矩阵方向自检和自动修正
4. 更详细的日志和错误提示

Author: r2end
Date: 2024-12-24
Version: 1.1 (IMPROVED)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import entropy, mannwhitneyu, kruskal
from scipy.cluster.hierarchy import dendrogram, linkage
import gc
import json
import warnings
import logging
warnings.filterwarnings('ignore')

# ============================================================================
# SETUP LOGGING
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# UTILITIES (从生产管线复用)
# ============================================================================

def safe_name(x: str, max_len: int = 180) -> str:
    """
    Sanitize names for filesystem compatibility
    与生产管线保持一致的命名规则
    """
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t']:
        s = s.replace(ch, '_')
    return s[:max_len]


def find_file(base_dir: Path, patterns: list) -> Path:
    """
    Find file matching any of the patterns
    按优先级尝试多种文件名模式
    """
    for pattern in patterns:
        files = list(base_dir.glob(pattern))
        if files:
            logger.info(f"  Found file: {files[0].name}")
            return files[0]
    return None


def validate_gep_matrix(gep_scores: pd.DataFrame, adata_genes: pd.Index) -> pd.DataFrame:
    """
    Validate and auto-correct GEP matrix orientation
    
    Expected:
        - Rows: GEPs (K, small number like 25-50)
        - Columns: Genes (large number like 20000)
    
    Returns:
        Corrected DataFrame with proper orientation
    """
    n_rows, n_cols = gep_scores.shape
    logger.info(f"  GEP matrix shape: {n_rows} × {n_cols}")
    
    # 方向检测逻辑
    needs_transpose = False
    
    # 策略1：如果columns数量远小于index，很可能是转置了
    if n_cols < n_rows and n_rows > 100:
        logger.warning(f"  ⚠️  Detected potential transpose (cols={n_cols} < rows={n_rows})")
        needs_transpose = True
    
    # 策略2：检查columns是否像基因名
    if not needs_transpose:
        # 检查columns是否与adata的基因名有重叠
        overlap_cols = len(set(gep_scores.columns) & set(adata_genes))
        overlap_idx = len(set(gep_scores.index) & set(adata_genes))
        
        logger.info(f"  Gene overlap: columns={overlap_cols}, index={overlap_idx}")
        
        if overlap_idx > overlap_cols * 2:  # index的重叠远大于columns
            logger.warning(f"  ⚠️  Detected transpose (index has more gene overlap)")
            needs_transpose = True
    
    # 执行转置（如果需要）
    if needs_transpose:
        logger.info(f"  → Transposing GEP matrix...")
        gep_scores = gep_scores.T
        logger.info(f"  → New shape: {gep_scores.shape}")
    
    # 最终验证：确保columns包含大量基因
    final_overlap = len(set(gep_scores.columns) & set(adata_genes))
    if final_overlap < 1000:
        logger.warning(f"  ⚠️  Low gene overlap ({final_overlap}), may need manual check")
    else:
        logger.info(f"  ✓ Validated: {final_overlap} genes overlap")
    
    return gep_scores
# ============================================================================
# GEP VISUALIZATION FUNCTIONS
# ============================================================================

def plot_gep_visualizations(
    adata: sc.AnnData,
    gep_name: str,
    gep_idx: int,
    celltype_col: str,
    cluster_col: str,
    output_dir: Path,
    figsize: tuple = (18, 5)
):
    """Create three visualization panels for a single GEP"""
    
    gep_usage = adata.obsm['X_cnmf_usages'][:, gep_idx]
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Panel A: UMAP
    ax = axes[0]
    if 'X_umap' not in adata.obsm:
        ax.text(0.5, 0.5, 'UMAP not available', ha='center', va='center', fontsize=12)
        ax.axis('off')
    else:
        umap_coords = adata.obsm['X_umap']
        scatter = ax.scatter(
            umap_coords[:, 0], umap_coords[:, 1],
            c=gep_usage, cmap='viridis', s=1, alpha=0.6, rasterized=True
        )
        ax.set_xlabel('UMAP 1', fontsize=10)
        ax.set_ylabel('UMAP 2', fontsize=10)
        ax.set_title(f'{gep_name} Usage on UMAP', fontsize=11, weight='bold')
        ax.axis('off')
        
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Usage', fontsize=9)
        cbar.ax.tick_params(labelsize=8)
    
    # Panel B: Cluster violin plot
    ax = axes[1]
    if cluster_col not in adata.obs.columns:
        ax.text(0.5, 0.5, f'Column "{cluster_col}" not found', ha='center', va='center')
        ax.axis('off')
    else:
        df_cluster = pd.DataFrame({
            'usage': gep_usage,
            'cluster': adata.obs[cluster_col].values
        }).dropna()
        
        cluster_medians = df_cluster.groupby('cluster')['usage'].median().sort_values(ascending=False)
        cluster_order = cluster_medians.index.tolist()[:20]  # Top 20
        
        df_cluster = df_cluster[df_cluster['cluster'].isin(cluster_order)]
        
        parts = ax.violinplot(
            [df_cluster[df_cluster['cluster'] == c]['usage'].values for c in cluster_order],
            positions=range(len(cluster_order)),
            widths=0.7, showmeans=True, showmedians=True
        )
        
        norm = plt.Normalize(vmin=cluster_medians.min(), vmax=cluster_medians.max())
        cmap = plt.cm.RdYlBu_r
        
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(cmap(norm(cluster_medians[cluster_order[i]])))
            pc.set_alpha(0.7)
        
        ax.set_xticks(range(len(cluster_order)))
        ax.set_xticklabels(cluster_order, rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('Cluster', fontsize=10)
        ax.set_ylabel('GEP Usage', fontsize=10)
        ax.set_title(f'{gep_name} by Cluster', fontsize=11, weight='bold')
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Panel C: Celltype violin plot
    ax = axes[2]
    if celltype_col not in adata.obs.columns:
        ax.text(0.5, 0.5, f'Column "{celltype_col}" not found', ha='center', va='center')
        ax.axis('off')
    else:
        df_celltype = pd.DataFrame({
            'usage': gep_usage,
            'celltype': adata.obs[celltype_col].values
        }).dropna()
        
        celltype_medians = df_celltype.groupby('celltype')['usage'].median().sort_values(ascending=False)
        celltype_order = celltype_medians.index.tolist()[:20]
        
        df_celltype = df_celltype[df_celltype['celltype'].isin(celltype_order)]
        
        parts = ax.violinplot(
            [df_celltype[df_celltype['celltype'] == c]['usage'].values for c in celltype_order],
            positions=range(len(celltype_order)),
            widths=0.7, showmeans=True, showmedians=True
        )
        
        norm = plt.Normalize(vmin=celltype_medians.min(), vmax=celltype_medians.max())
        cmap = plt.cm.RdYlBu_r
        
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(cmap(norm(celltype_medians[celltype_order[i]])))
            pc.set_alpha(0.7)
        
        ax.set_xticks(range(len(celltype_order)))
        ax.set_xticklabels(celltype_order, rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('Cell Type', fontsize=10)
        ax.set_ylabel('GEP Usage', fontsize=10)
        ax.set_title(f'{gep_name} by Cell Type', fontsize=11, weight='bold')
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    
    plt.tight_layout()
    
    output_file = output_dir / f"{gep_name}_detailed_visualization.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file


# ============================================================================
# CONFIGURATION - 请根据你的数据修改这部分
# ============================================================================

# 基本设置
DATASET_ID = "AT2"              # 你的数据集ID
TRACK = "batch_aware"                # 推荐使用batch_aware
K_VALUE = 35                        # 你选定的K值

# 路径设置
BASE_DIR = Path("/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2")
CNMF_OUTPUT_DIR = BASE_DIR / DATASET_ID / TRACK / "cnmf_output"

# ⭐ 关键改进：使用 safe_name() 确保与生产管线一致
CNMF_NAME = safe_name(f"{DATASET_ID}_{TRACK}")

# 原始数据路径（包含细胞类型注释）
ORIGINAL_H5AD = Path('/home/h2048/data/R/1223/per_celltype_harmony_rogue/h5ad_objects/Club_harmony.h5ad')

# 输出目录
OUTPUT_BASE = BASE_DIR / DATASET_ID / TRACK / f"cnmf_analysis_k{K_VALUE}"

# 列名设置（根据你的数据调整）
CELLTYPE_COL = 'RNA_snn_res.1.5'           # 细胞类型列名
DISEASE_COL = 'disease_status'       # 疾病状态列名（如果有）
BATCH_COL = 'sample'                 # 批次列名
CLUSTER_COL = 'RNA_snn_res.1.5'  # cluster列名，会自动检测

# 分析参数
TOP_N_GENES = 50                     # 每个GEP提取的top基因数
MIN_CELLS_PER_TYPE = 10              # 最小细胞数阈值
P_VALUE_THRESHOLD = 0.05             # 显著性阈值
VISUALIZE_TOP_N_GEPS = None  # None=全部，或设置数字如10, 20
VISUALIZATION_BATCH_SIZE = 5
# ============================================================================
# 不需要修改以下代码
# ============================================================================

logger.info("="*70)
logger.info("cNMF DOWNSTREAM ANALYSIS - IMPROVED VERSION")
logger.info("="*70)
logger.info(f"\nDataset: {DATASET_ID}")
logger.info(f"Track: {TRACK}")
logger.info(f"K value: {K_VALUE}")
logger.info(f"Safe name: {CNMF_NAME}")
logger.info(f"\n" + "="*70)

# 创建输出目录
OUTPUT_BASE.mkdir(exist_ok=True, parents=True)
(OUTPUT_BASE / "figures").mkdir(exist_ok=True)
(OUTPUT_BASE / "tables").mkdir(exist_ok=True)

# ============================================================================
# PART 1: 加载和整合cNMF结果
# ============================================================================

logger.info(f"\n{'='*70}")
logger.info("PART 1: Loading cNMF Results")
logger.info("="*70)

# cNMF输出子目录（使用safe_name后的名称）
cnmf_name_dir = CNMF_OUTPUT_DIR / CNMF_NAME

if not cnmf_name_dir.exists():
    logger.error(f"\n❌ ERROR: cNMF output directory not found")
    logger.error(f"  Expected: {cnmf_name_dir}")
    logger.error(f"\nPossible reasons:")
    logger.error(f"  1. cNMF pipeline not yet run")
    logger.error(f"  2. DATASET_ID or TRACK incorrect")
    logger.error(f"  3. Name mismatch (check if safe_name() was used in production)")
    
    # 尝试列出可用目录
    if CNMF_OUTPUT_DIR.exists():
        available = list(CNMF_OUTPUT_DIR.glob("*"))
        if available:
            logger.info(f"\n  Available directories in {CNMF_OUTPUT_DIR}:")
            for d in available[:10]:
                logger.info(f"    - {d.name}")
    
    exit(1)

logger.info(f"\n✓ cNMF output directory found: {cnmf_name_dir.name}")

# ⭐ 改进：增强文件模式（包含 spectra.*consensus）
logger.info(f"\nLoading GEP-gene scores...")
gep_file = find_file(
    cnmf_name_dir,
    [
        # 优先级1：标准gene_spectra_score文件
        f"{CNMF_NAME}.gene_spectra_score.k_{K_VALUE}.dt_0.1.txt",
        f"{CNMF_NAME}.gene_spectra_score.k_{K_VALUE}.dt_0_1.txt",
        
        # 优先级2：spectra consensus文件（某些cNMF配置）
        f"{CNMF_NAME}.spectra.k_{K_VALUE}.dt_0.1.consensus.txt",
        f"{CNMF_NAME}.spectra.k_{K_VALUE}.dt_0_1.consensus.txt",
        
        # 优先级3：通配符兜底
        f"*gene_spectra_score.k_{K_VALUE}*.txt",
        f"*spectra.k_{K_VALUE}*consensus*.txt",
    ]
)

if not gep_file:
    logger.error(f"\n❌ ERROR: GEP file not found in {cnmf_name_dir}")
    logger.error(f"\nSearched for patterns:")
    logger.error(f"  - {CNMF_NAME}.gene_spectra_score.k_{K_VALUE}.dt_*.txt")
    logger.error(f"  - {CNMF_NAME}.spectra.k_{K_VALUE}.dt_*.consensus.txt")
    
    # 列出实际存在的文件
    actual_files = list(cnmf_name_dir.glob(f"*k_{K_VALUE}*.txt"))
    if actual_files:
        logger.info(f"\nActual files found for K={K_VALUE}:")
        for f in actual_files[:10]:
            logger.info(f"  - {f.name}")
    
    exit(1)

gep_scores = pd.read_csv(gep_file, sep='\t', index_col=0)
logger.info(f"  ✓ Loaded GEP matrix (raw): {gep_scores.shape}")

# 加载细胞-GEP使用
logger.info(f"\nLoading cell-GEP usage...")
usage_file = find_file(
    cnmf_name_dir,
    [
        f"{CNMF_NAME}.usages.k_{K_VALUE}.dt_0.1.consensus.txt",
        f"{CNMF_NAME}.usages.k_{K_VALUE}.dt_0_1.consensus.txt",
        f"*usages.k_{K_VALUE}*consensus.txt"
    ]
)

if not usage_file:
    logger.error(f"❌ ERROR: Usage file not found in {cnmf_name_dir}")
    exit(1)

usages = pd.read_csv(usage_file, sep='\t', index_col=0)
logger.info(f"  ✓ Cell-GEP usage matrix: {usages.shape}")

# 加载原始数据
logger.info(f"\nLoading original data...")

if not ORIGINAL_H5AD.exists():
    logger.error(f"❌ ERROR: Original h5ad file not found")
    logger.error(f"  Expected: {ORIGINAL_H5AD}")
    exit(1)

adata = sc.read_h5ad(ORIGINAL_H5AD)
logger.info(f"  ✓ Original data: {adata.n_obs} cells × {adata.n_vars} genes")

# ⭐ 关键改进：GEP矩阵方向验证和自动修正
logger.info(f"\nValidating GEP matrix orientation...")
gep_scores = validate_gep_matrix(gep_scores, adata.var_names)

# 对齐和整合
logger.info(f"\nIntegrating results...")
common_cells = adata.obs_names.intersection(usages.index)
overlap_pct = len(common_cells) / adata.n_obs * 100

logger.info(f"  Overlapping cells: {len(common_cells)} ({overlap_pct:.1f}%)")

if overlap_pct < 90:
    logger.warning(f"  ⚠️  Warning: <90% cells overlap, check cell IDs")

usages_aligned = usages.loc[common_cells]
adata = adata[common_cells].copy()

# 添加到AnnData
adata.obsm['X_cnmf_usages'] = usages_aligned.values
adata.uns['cnmf_k'] = K_VALUE
adata.uns['cnmf_gep_names'] = usages_aligned.columns.tolist()

# 对齐基因
common_genes = adata.var_names.intersection(gep_scores.columns)
gene_overlap_pct = len(common_genes) / adata.n_vars * 100

logger.info(f"  Overlapping genes: {len(common_genes)} ({gene_overlap_pct:.1f}%)")

if gene_overlap_pct < 50:
    logger.warning(f"  ⚠️  Warning: <50% genes overlap, check gene IDs")

gep_scores_aligned = gep_scores[common_genes].T  # 转置：(n_genes, K)
adata = adata[:, common_genes].copy()

adata.varm['cnmf_spectra'] = gep_scores_aligned.values

# 计算衍生指标
adata.obs['dominant_gep'] = usages_aligned.idxmax(axis=1).values
adata.obs['max_gep_usage'] = usages_aligned.max(axis=1).values
adata.obs['gep_entropy'] = usages_aligned.apply(
    lambda row: entropy(row / row.sum()), axis=1
).values

# 保存
output_h5ad = OUTPUT_BASE / f"{DATASET_ID}_with_cnmf_k{K_VALUE}.h5ad"
adata.write_h5ad(output_h5ad, compression='gzip')
logger.info(f"\n✓ Saved integrated data: {output_h5ad.name}")

# 保存矩阵（方便后续分析）
gep_scores_aligned.to_csv(OUTPUT_BASE / "tables" / f"gep_gene_scores_k{K_VALUE}.csv")
usages_aligned.to_csv(OUTPUT_BASE / "tables" / f"cell_gep_usages_k{K_VALUE}.csv")

# ============================================================================
# PART 2: GEP功能注释
# ============================================================================

logger.info(f"\n{'='*70}")
logger.info("PART 2: GEP Functional Annotation")
logger.info("="*70)

# 提取top基因
logger.info(f"\nExtracting top {TOP_N_GENES} genes per GEP...")
gep_scores_df = pd.DataFrame(
    adata.varm['cnmf_spectra'],
    index=adata.var_names,
    columns=adata.uns['cnmf_gep_names']
)

top_genes_dict = {}
for gep in gep_scores_df.columns:
    top_genes = gep_scores_df[gep].sort_values(ascending=False).head(TOP_N_GENES)
    top_genes_dict[gep] = {
        'genes': top_genes.index.tolist(),
        'scores': top_genes.values.tolist()
    }

# 保存top基因
with open(OUTPUT_BASE / "tables" / "top_genes_per_gep.json", 'w') as f:
    json.dump(top_genes_dict, f, indent=2)

# 打印前3个GEP的top 10基因
logger.info(f"\nTop 10 genes for first 3 GEPs:")
for i, gep in enumerate(list(gep_scores_df.columns)[:3]):
    logger.info(f"\n{gep}:")
    logger.info(f"  {', '.join(top_genes_dict[gep]['genes'][:10])}")

# 手动功能评分
logger.info(f"\nCalculating functional signature scores...")

FUNCTIONAL_SIGNATURES = {
    'Cell Cycle': ['MKI67', 'TOP2A', 'PCNA', 'CCNA2', 'CCNB1'],
    'Interferon': ['IFIT1', 'IFIT2', 'IFIT3', 'ISG15', 'MX1'],
    'Inflammatory': ['IL1B', 'IL6', 'TNF', 'CXCL8', 'CCL2'],
    'Cytotoxicity': ['GZMB', 'GZMA', 'PRF1', 'GNLY', 'NKG7'],
}

gep_functional_scores = pd.DataFrame(index=gep_scores_df.columns)

for func_name, genes in FUNCTIONAL_SIGNATURES.items():
    available_genes = [g for g in genes if g in gep_scores_df.index]
    if len(available_genes) > 0:
        scores = gep_scores_df.loc[available_genes].mean(axis=0)
        gep_functional_scores[func_name] = scores

# 归一化
if len(gep_functional_scores.columns) > 0:
    gep_functional_scores = (gep_functional_scores - gep_functional_scores.min()) / \
                            (gep_functional_scores.max() - gep_functional_scores.min())
    
    # 可视化
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(gep_functional_scores.T, cmap='RdYlBu_r', 
                cbar_kws={'label': 'Signature Score'},
                linewidths=0.5, ax=ax)
    ax.set_xlabel('GEP', fontsize=11)
    ax.set_ylabel('Functional Signature', fontsize=11)
    ax.set_title(f'GEP Functional Annotation (K={K_VALUE})', fontsize=13, weight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_BASE / "figures" / "gep_functional_heatmap.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved functional heatmap")
    
    # 保存
    gep_functional_scores.to_csv(OUTPUT_BASE / "tables" / "gep_functional_scores.csv")

# ============================================================================
# PART 3: GEP-细胞类型关联
# ============================================================================

logger.info(f"\n{'='*70}")
logger.info("PART 3: GEP-Cell Type Association")
logger.info("="*70)

if CELLTYPE_COL not in adata.obs.columns:
    logger.warning(f"  ⚠️  Column '{CELLTYPE_COL}' not found in data")
    logger.info(f"  Available columns: {list(adata.obs.columns[:20])}")
    logger.info(f"  Skipping cell type analysis...")
else:
    logger.info(f"\nCalculating mean GEP usage by cell type...")
    
    gep_usage = pd.DataFrame(
        adata.obsm['X_cnmf_usages'],
        index=adata.obs_names,
        columns=adata.uns['cnmf_gep_names']
    )
    gep_usage['cell_type'] = adata.obs[CELLTYPE_COL].values
    
    mean_usage = gep_usage.groupby('cell_type').mean()
    
    # 可视化
    fig, ax = plt.subplots(figsize=(14, max(8, len(mean_usage) * 0.4)))
    sns.heatmap(mean_usage, cmap='RdYlBu_r', 
                cbar_kws={'label': 'Mean Usage'},
                linewidths=0.5, ax=ax)
    ax.set_xlabel('GEP', fontsize=11)
    ax.set_ylabel('Cell Type', fontsize=11)
    ax.set_title(f'Mean GEP Usage by Cell Type (K={K_VALUE})', fontsize=13, weight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_BASE / "figures" / "gep_celltype_heatmap.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved cell type heatmap")
    
    # 保存
    mean_usage.to_csv(OUTPUT_BASE / "tables" / "mean_gep_usage_by_celltype.csv")
    
    # 统计检验
    logger.info(f"\nTesting GEP-cell type specificity...")
    
    specificity_results = []
    
    for gep in adata.uns['cnmf_gep_names']:
        groups = [gep_usage[gep_usage['cell_type'] == ct][gep].values 
                  for ct in mean_usage.index]
        
        stat, pval = kruskal(*groups)
        
        max_celltype = mean_usage[gep].idxmax()
        max_usage = mean_usage[gep].max()
        specificity_score = max_usage / mean_usage[gep].mean()
        
        specificity_results.append({
            'GEP': gep,
            'max_cell_type': max_celltype,
            'max_usage': max_usage,
            'specificity_score': specificity_score,
            'p_value': pval
        })
    
    specificity_df = pd.DataFrame(specificity_results)
    
    # FDR校正
    from statsmodels.stats.multitest import multipletests
    specificity_df['p_adj'] = multipletests(specificity_df['p_value'], method='fdr_bh')[1]
    specificity_df = specificity_df.sort_values('specificity_score', ascending=False)
    
    logger.info(f"\nTop 10 cell-type-specific GEPs:")
    for idx, row in specificity_df.head(10).iterrows():
        logger.info(f"  {row['GEP']}: {row['max_cell_type']} "
                   f"(score={row['specificity_score']:.2f}, p_adj={row['p_adj']:.2e})")
    
    # 保存
    specificity_df.to_csv(OUTPUT_BASE / "tables" / "gep_celltype_specificity.csv", index=False)

# ============================================================================
# PART 4: GEP-疾病状态关联（如果有）
# ============================================================================

logger.info(f"\n{'='*70}")
logger.info("PART 4: GEP-Disease Association (Optional)")
logger.info("="*70)

if DISEASE_COL not in adata.obs.columns:
    logger.info(f"  Column '{DISEASE_COL}' not found, skipping disease analysis")
elif pd.isna(adata.obs[DISEASE_COL]).all():
    logger.info(f"  Column '{DISEASE_COL}' is all NaN, skipping disease analysis")
else:
    logger.info(f"\nComparing GEP usage across disease groups...")
    
    gep_usage = pd.DataFrame(
        adata.obsm['X_cnmf_usages'],
        index=adata.obs_names,
        columns=adata.uns['cnmf_gep_names']
    )
    gep_usage['disease'] = adata.obs[DISEASE_COL].values
    
    disease_groups = gep_usage['disease'].unique()
    disease_groups = disease_groups[~pd.isna(disease_groups)]
    
    if len(disease_groups) == 2:
        group1, group2 = disease_groups
        
        diff_results = []
        
        for gep in adata.uns['cnmf_gep_names']:
            g1_values = gep_usage[gep_usage['disease'] == group1][gep]
            g2_values = gep_usage[gep_usage['disease'] == group2][gep]
            
            stat, pval = mannwhitneyu(g1_values, g2_values, alternative='two-sided')
            
            mean_diff = g1_values.mean() - g2_values.mean()
            pooled_std = np.sqrt((g1_values.std()**2 + g2_values.std()**2) / 2)
            cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
            
            diff_results.append({
                'GEP': gep,
                f'{group1}_mean': g1_values.mean(),
                f'{group2}_mean': g2_values.mean(),
                'mean_diff': mean_diff,
                'cohens_d': cohens_d,
                'p_value': pval
            })
        
        diff_df = pd.DataFrame(diff_results)
        diff_df['p_adj'] = multipletests(diff_df['p_value'], method='fdr_bh')[1]
        diff_df = diff_df.sort_values('mean_diff', key=abs, ascending=False)
        
        logger.info(f"\nTop 10 differential GEPs ({group1} vs {group2}):")
        for idx, row in diff_df.head(10).iterrows():
            direction = "↑" if row['mean_diff'] > 0 else "↓"
            logger.info(f"  {row['GEP']}: {direction} "
                       f"(d={row['cohens_d']:.2f}, p_adj={row['p_adj']:.2e})")
        
        # 保存
        diff_df.to_csv(OUTPUT_BASE / "tables" / "gep_disease_differential.csv", index=False)
        
        # Volcano plot
        plt.figure(figsize=(10, 8))
        
        significant = diff_df['p_adj'] < P_VALUE_THRESHOLD
        large_effect = abs(diff_df['cohens_d']) > 0.5
        
        colors = []
        for sig, large in zip(significant, large_effect):
            if sig and large:
                colors.append('red')
            elif sig:
                colors.append('orange')
            else:
                colors.append('gray')
        
        plt.scatter(diff_df['cohens_d'], -np.log10(diff_df['p_adj']), 
                    c=colors, alpha=0.6, s=50)
        
        # 标注
        for idx, row in diff_df.head(5).iterrows():
            plt.annotate(row['GEP'], 
                        (row['cohens_d'], -np.log10(row['p_adj'])),
                        fontsize=8)
        
        plt.axhline(-np.log10(P_VALUE_THRESHOLD), linestyle='--', 
                    color='black', alpha=0.3)
        plt.axvline(-0.5, linestyle='--', color='blue', alpha=0.3)
        plt.axvline(0.5, linestyle='--', color='blue', alpha=0.3)
        
        plt.xlabel("Cohen's d", fontsize=11)
        plt.ylabel('-log10(FDR-adjusted p)', fontsize=11)
        plt.title(f'Differential GEP Usage: {group1} vs {group2}', 
                  fontsize=13, weight='bold')
        plt.tight_layout()
        plt.savefig(OUTPUT_BASE / "figures" / "gep_disease_volcano.png", 
                    dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  ✓ Saved volcano plot")


# ============================================================================
# PART 4: Per-GEP Detailed Visualization ⭐ NEW
# ============================================================================

logger.info(f"\n{'='*70}")
logger.info("PART 4: Per-GEP Detailed Visualization")
logger.info("="*70)

# Create output directory
gep_viz_dir = OUTPUT_BASE / "gep_visualizations"
gep_viz_dir.mkdir(exist_ok=True, parents=True)

# Get GEP names
gep_names = adata.uns['cnmf_gep_names']
n_geps = len(gep_names)

# Select GEPs to visualize
if VISUALIZE_TOP_N_GEPS and VISUALIZE_TOP_N_GEPS < n_geps:
    logger.info(f"\nSelecting top {VISUALIZE_TOP_N_GEPS} most variable GEPs...")
    gep_usage = adata.obsm['X_cnmf_usages']
    gep_vars = np.var(gep_usage, axis=0)
    top_indices = np.argsort(gep_vars)[::-1][:VISUALIZE_TOP_N_GEPS]
    selected_geps = [(i, gep_names[i]) for i in top_indices]
    logger.info(f"  Selected: {', '.join([g for _, g in selected_geps])}")
else:
    selected_geps = list(enumerate(gep_names))
    logger.info(f"\nVisualizing all {n_geps} GEPs...")

# Ensure UMAP exists
if 'X_umap' not in adata.obsm:
    logger.info(f"\n  Computing UMAP...")
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_cnmf_usages')
    sc.tl.umap(adata)
    logger.info(f"  ✓ UMAP computed")

# Find or compute cluster column
if CLUSTER_COL not in adata.obs.columns:
    cluster_candidates = ['leiden', 'louvain', 'seurat_clusters', 'cluster']
    found = False
    for candidate in cluster_candidates:
        if candidate in adata.obs.columns:
            CLUSTER_COL = candidate
            found = True
            logger.info(f"  Using '{CLUSTER_COL}' as cluster column")
            break
    
    if not found:
        logger.info(f"  Computing Leiden clustering...")
        sc.tl.leiden(adata, resolution=1.0, key_added='leiden')
        CLUSTER_COL = 'leiden'

# Generate visualizations
logger.info(f"\nGenerating visualizations...")

successful = 0
failed = []

for batch_start in range(0, len(selected_geps), VISUALIZATION_BATCH_SIZE):
    batch_end = min(batch_start + VISUALIZATION_BATCH_SIZE, len(selected_geps))
    batch_geps = selected_geps[batch_start:batch_end]
    
    logger.info(f"  Batch {batch_start//VISUALIZATION_BATCH_SIZE + 1}: "
               f"Processing {batch_start+1}-{batch_end}/{len(selected_geps)}...")
    
    for gep_idx, gep_name in batch_geps:
        try:
            plot_gep_visualizations(
                adata=adata,
                gep_name=gep_name,
                gep_idx=gep_idx,
                celltype_col=CELLTYPE_COL,
                cluster_col=CLUSTER_COL,
                output_dir=gep_viz_dir
            )
            successful += 1
        except Exception as e:
            logger.warning(f"    ⚠️  Failed {gep_name}: {str(e)}")
            failed.append(gep_name)

# Summary
logger.info(f"\n{'='*70}")
logger.info(f"Visualization Summary:")
logger.info(f"  ✓ Successful: {successful}/{len(selected_geps)}")

if failed:
    logger.warning(f"  ⚠️  Failed: {len(failed)} - {', '.join(failed)}")

logger.info(f"\n  Output: {gep_viz_dir}")

# Create HTML index
logger.info(f"\nCreating HTML index...")

html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>GEP Visualizations - {DATASET_ID}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }}
        h1 {{ color: #333; border-bottom: 3px solid #4CAF50; padding-bottom: 10px; }}
        .summary {{ background: white; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .gallery {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(350px, 1fr)); gap: 20px; }}
        .gep-card {{ background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .gep-card img {{ width: 100%; cursor: pointer; transition: transform 0.2s; }}
        .gep-card img:hover {{ transform: scale(1.02); }}
        .gep-card h3 {{ margin: 10px; color: #4CAF50; }}
        .modal {{ display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; 
                 background-color: rgba(0,0,0,0.9); }}
        .modal-content {{ margin: auto; display: block; max-width: 95%; max-height: 95%; padding-top: 50px; }}
        .close {{ position: absolute; top: 15px; right: 35px; color: #f1f1f1; font-size: 40px; 
                 font-weight: bold; cursor: pointer; }}
    </style>
</head>
<body>
    <h1>🧬 {DATASET_ID} - GEP Detailed Visualizations (K={K_VALUE})</h1>
    
    <div class="summary">
        <h2>Summary</h2>
        <p><strong>Total visualized:</strong> {successful} GEPs</p>
        <p><strong>Each figure contains 3 panels:</strong></p>
        <ul>
            <li><strong>Panel A:</strong> UMAP colored by GEP usage (viridis colormap)</li>
            <li><strong>Panel B:</strong> Violin plot distribution across clusters (sorted by median)</li>
            <li><strong>Panel C:</strong> Violin plot distribution across cell types (sorted by median)</li>
        </ul>
        <p><em>💡 Tip: Click any image to view full size. Use browser back button to return.</em></p>
    </div>
    
    <div class="gallery">
"""

for gep_idx, gep_name in selected_geps:
    img_file = f"{gep_name}_detailed_visualization.png"
    if (gep_viz_dir / img_file).exists():
        html_content += f"""
        <div class="gep-card">
            <img src="{img_file}" alt="{gep_name}" onclick="openModal('{img_file}')">
            <h3>{gep_name}</h3>
        </div>
"""

html_content += """
    </div>
    
    <div id="myModal" class="modal" onclick="closeModal()">
        <span class="close">&times;</span>
        <img class="modal-content" id="modalImg">
    </div>
    
    <script>
        function openModal(src) {
            document.getElementById('myModal').style.display = 'block';
            document.getElementById('modalImg').src = src;
        }
        function closeModal() {
            document.getElementById('myModal').style.display = 'none';
        }
    </script>
</body>
</html>
"""

index_file = gep_viz_dir / "index.html"
with open(index_file, 'w', encoding='utf-8') as f:
    f.write(html_content)

logger.info(f"  ✓ Created: {index_file}")
logger.info(f"    → Open in browser to browse all visualizations")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

logger.info(f"\n{'='*70}")
logger.info("ANALYSIS COMPLETE!")
logger.info("="*70)

logger.info(f"\nOutput directory: {OUTPUT_BASE}")
logger.info(f"\nGenerated files:")
logger.info(f"  ✓ {DATASET_ID}_with_cnmf_k{K_VALUE}.h5ad")
logger.info(f"  ✓ figures/*.png")
logger.info(f"  ✓ tables/*.csv")
logger.info(f"  ✓ gep_visualizations/*.png ({successful} GEPs)")
logger.info(f"  ✓ gep_visualizations/index.html")

logger.info(f"\nQuick access:")
logger.info(f"  HTML Index: file://{index_file}")

logger.info(f"\n{'='*70}")






