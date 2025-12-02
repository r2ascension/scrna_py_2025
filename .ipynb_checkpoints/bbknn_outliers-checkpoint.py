#!/usr/bin/env python3
"""
rare_cells_subset_analysis.py

提取在dataset中细胞数<10的所有细胞进行子集分析
- 识别所有稀有细胞(dataset中<10个)
- 提取这些细胞进行重新聚类和注释
- 生成marker基因分析和可视化
- 提供重新注释建议

作者: 临床-生信团队  
日期: 2025-11-07
版本: v1.0
"""

import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats
warnings.filterwarnings('ignore')

# ==================== 配置部分 ====================

# ========== 运行模式配置 ==========
# Mode 1: "pre_annotation"  - 第一次运行，生成marker分析和annotation模板
# Mode 2: "post_annotation" - 第二次运行，应用手动annotation
# Mode 3: "full"            - 完整运行（包含自动注释建议，但不推荐）
RUN_MODE = "post_annotation"  # 选项: "pre_annotation", "post_annotation", "full"

# ========== 输入输出配置 ==========
input_h5ad_path = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated.h5ad"
output_dir = "/home/h2048/data/py/1029/rare_cells_analysis"
annotation_file = "/home/h2048/data/py/1029/rare_cells_analysis/Annotation.csv"  # 第一列:Cluster, 第二列:CellType (无表头)

# ========== 关键列名 ==========
cell_type_key = "cell_type"
cluster_key = "leiden_bbknn"
batch_key = "dataset"
tissue_key = "tissue_sampling_method"

# ========== 稀有细胞阈值 ==========
max_cells_threshold = 10  # dataset中少于此数量的细胞类型视为稀有

# ========== 预处理配置 ==========
use_hvg = True
n_top_genes = 2000
min_cells_per_gene = 3
scale_max_value = 10

# ========== PCA配置 ==========
n_pcs = 50

# ========== 聚类配置 ==========
leiden_resolutions = [1.3, 1.7, 2.0, 2.5]
default_resolution = 2.0

# ========== Marker基因配置 ==========
marker_min_pct = 0.25
marker_logfc_threshold = 0.25
top_n_markers = 20

# ========== Marker基因列表 ==========
celltype_markers = {
    "Epithelial": ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1', 'FXYD3', 'ELF3'],
    "T": ['CD3D', 'CD3E', 'CD3G', 'CD2', 'TRBC2', 'CD4', 'CD8A', 'CD8B'],
    "B": ['CD19', 'CD79A', 'MS4A1', 'CD79B', 'TNFRSF17'],
    "Myeloid": ['CD14', 'CD68', 'LYZ', 'CD163', 'FCGR3A', 'FCER1G', 'CD86'],
    "Endothelial": ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'ECSCR'],
    "Fibroblast": ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'COL6A3'],
    "SMC": ['ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'DES'],
    "Proliferation": ['MKI67', 'TOP2A', 'PCNA', 'TK1']
}

# ========== 可视化配置 ==========
dpi = 300
figure_format = "pdf"

# 细胞类型颜色
celltype_colors = {
    "Epithelial": "#E41A1C",
    "Fibroblast": "#377EB8",
    "T": "#4DAF4A",
    "Myeloid": "#984EA3",
    "Endothelial": "#FF7F00",
    "SMC": "#FFFF33",
    "B": "#A65628",
    "Proliferation": "#F781BF"
}

verbose = True

# ==================== 辅助函数 ====================

def log_msg(msg):
    if verbose:
        print(msg)

def log_step(step_num, step_name):
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)

# ==================== Step 1: 加载数据 ====================

log_step(1, "Loading Data")

output_path = Path(output_dir)
output_path.mkdir(parents=True, exist_ok=True)

fig_dir = output_path / "figures"
fig_dir.mkdir(exist_ok=True)

log_msg(f"\n读取文件: {input_h5ad_path}")
adata_full = sc.read_h5ad(input_h5ad_path)

log_msg(f"✓ 数据加载成功")
log_msg(f"   细胞总数: {adata_full.n_obs:,}")
log_msg(f"   基因总数: {adata_full.n_vars:,}")

# ==================== Step 2: 识别稀有细胞 ====================

log_step(2, "Identifying Rare Cells")

# 统计每个dataset-celltype组合的细胞数
log_msg(f"\n统计每个dataset中各细胞类型的细胞数...")
dataset_celltype_counts = pd.crosstab(
    adata_full.obs[batch_key],
    adata_full.obs[cell_type_key]
)

log_msg(f"  Dataset × CellType矩阵: {dataset_celltype_counts.shape}")

# 保存完整统计
stats_file = output_path / "dataset_celltype_counts.csv"
dataset_celltype_counts.to_csv(stats_file)
log_msg(f"  ✓ 保存: {stats_file}")

# 识别稀有细胞
log_msg(f"\n识别细胞数 < {max_cells_threshold} 的dataset-celltype组合...")

rare_combinations = []
for dataset in dataset_celltype_counts.index:
    for celltype in dataset_celltype_counts.columns:
        count = dataset_celltype_counts.loc[dataset, celltype]
        if 0 < count < max_cells_threshold:
            rare_combinations.append({
                'Dataset': dataset,
                'CellType': celltype,
                'Count': count
            })

rare_df = pd.DataFrame(rare_combinations)

log_msg(f"\n发现 {len(rare_df)} 个稀有组合:")
log_msg(f"  涉及数据集: {rare_df['Dataset'].nunique()}")
log_msg(f"  涉及细胞类型: {rare_df['CellType'].nunique()}")

# 按细胞类型分组显示
log_msg(f"\n稀有组合按细胞类型统计:")
for celltype in sorted(rare_df['CellType'].unique()):
    ct_data = rare_df[rare_df['CellType'] == celltype]
    n_datasets = len(ct_data)
    total_cells = ct_data['Count'].sum()
    log_msg(f"  {celltype}: {n_datasets} datasets, 共 {total_cells} 个稀有细胞")

rare_file = output_path / "rare_combinations.csv"
rare_df.to_csv(rare_file, index=False)
log_msg(f"\n✓ 保存稀有组合: {rare_file}")

# ==================== Step 3: 提取稀有细胞 ====================

log_step(3, "Extracting Rare Cells")

# 获取所有稀有细胞的cell_id
rare_cell_ids = []
cell_info_list = []

log_msg(f"\n提取稀有细胞...")

for _, row in rare_df.iterrows():
    dataset = row['Dataset']
    celltype = row['CellType']
    
    # 找到符合条件的细胞
    mask = (adata_full.obs[batch_key] == dataset) & \
           (adata_full.obs[cell_type_key] == celltype)
    
    cells = adata_full.obs[mask].index.tolist()
    rare_cell_ids.extend(cells)
    
    # 记录细胞信息
    for cell_id in cells:
        cell_obs = adata_full.obs.loc[cell_id]
        cell_info_list.append({
            'Cell_ID': cell_id,
            'Original_CellType': celltype,
            'Original_Cluster': cell_obs[cluster_key],
            'Dataset': dataset,
            'Tissue': cell_obs[tissue_key],
            'Is_Rare': True
        })

log_msg(f"  提取到 {len(rare_cell_ids)} 个稀有细胞")

# 保存稀有细胞详细信息
cell_info_df = pd.DataFrame(cell_info_list)
cell_info_file = output_path / "rare_cells_info.csv"
cell_info_df.to_csv(cell_info_file, index=False)
log_msg(f"  ✓ 保存细胞信息: {cell_info_file}")

# 创建子集
log_msg(f"\n创建稀有细胞子集...")
adata_rare = adata_full[rare_cell_ids, :].copy()

log_msg(f"  子集维度: {adata_rare.n_obs:,} cells × {adata_rare.n_vars:,} genes")

# 统计子集中的细胞类型分布
log_msg(f"\n子集中的细胞类型分布:")
subset_celltype_counts = adata_rare.obs[cell_type_key].value_counts()
for celltype, count in subset_celltype_counts.items():
    pct = count / adata_rare.n_obs * 100
    log_msg(f"  {celltype}: {count} cells ({pct:.1f}%)")

# ==================== Step 4: 预处理 ====================

log_step(4, "Preprocessing Rare Cells Subset")

# 4.1 基因过滤
log_msg(f"\n过滤低表达基因 (min_cells={min_cells_per_gene})...")
sc.pp.filter_genes(adata_rare, min_cells=min_cells_per_gene)
log_msg(f"  保留基因数: {adata_rare.n_vars:,}")

# 4.2 标准化
log_msg(f"\n标准化数据...")
sc.pp.normalize_total(adata_rare, target_sum=1e4)
sc.pp.log1p(adata_rare)

# 4.3 高变基因选择
if use_hvg:
    log_msg(f"\n选择高变基因 (n_top_genes={n_top_genes})...")
    
    # 检查batch大小
    batch_counts = adata_rare.obs[batch_key].value_counts()
    min_batch_size = batch_counts.min()
    
    if min_batch_size < 3:
        log_msg(f"  警告: 存在batch细胞数<3, 不使用batch_key")
        sc.pp.highly_variable_genes(
            adata_rare,
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
            subset=False
        )
    else:
        sc.pp.highly_variable_genes(
            adata_rare,
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
            batch_key=batch_key,
            subset=False
        )
    
    n_hvg = adata_rare.var['highly_variable'].sum()
    log_msg(f"  选中 {n_hvg} 个高变基因")
    
    # 子集到高变基因用于降维
    adata_hvg = adata_rare[:, adata_rare.var['highly_variable']].copy()
else:
    adata_hvg = adata_rare.copy()

# 4.4 Scale
log_msg(f"\nScale数据...")
sc.pp.scale(adata_hvg, max_value=scale_max_value)

log_msg(f"\n预处理完成")
log_msg(f"  最终维度: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")

# ==================== Step 5: 降维和聚类 ====================

log_step(5, "Dimensionality Reduction and Clustering")

# 5.1 PCA
log_msg(f"\n运行PCA (n_comps={n_pcs})...")
sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver='arpack')

var_ratio = adata_hvg.uns['pca']['variance_ratio']
cumsum_var = np.cumsum(var_ratio)
log_msg(f"  PC1-10 累计解释方差: {cumsum_var[9]:.2%}")
log_msg(f"  PC1-30 累计解释方差: {cumsum_var[29]:.2%}")

# 5.2 Neighbors
log_msg(f"\n计算邻居图...")
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=n_pcs)

# 5.3 UMAP
log_msg(f"\n运行UMAP...")
sc.tl.umap(adata_hvg, min_dist=0.3)

# 5.4 Leiden聚类(多个分辨率)
log_msg(f"\n运行Leiden聚类...")
for res in leiden_resolutions:
    key_added = f"leiden_rare_res{res}"
    sc.tl.leiden(adata_hvg, resolution=res, key_added=key_added)
    n_clusters = adata_hvg.obs[key_added].nunique()
    log_msg(f"  Resolution {res}: {n_clusters} clusters")

# 将结果复制回完整数据
for key in adata_hvg.obs.columns:
    if key.startswith('leiden_rare'):
        adata_rare.obs[key] = adata_hvg.obs[key]

adata_rare.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata_rare.obsm['X_umap'] = adata_hvg.obsm['X_umap']

log_msg(f"\n✓ 降维和聚类完成")

# ==================== Step 6: 可视化 ====================

log_step(6, "Visualization")

# 6.1 原始细胞类型UMAP
log_msg(f"\n生成原始细胞类型UMAP...")

fig, ax = plt.subplots(figsize=(10, 8))

# 获取颜色
celltypes_in_subset = adata_rare.obs[cell_type_key].unique()
colors = [celltype_colors.get(ct, '#808080') for ct in celltypes_in_subset]

sc.pl.umap(
    adata_rare,
    color=cell_type_key,
    palette=colors,
    legend_loc='right margin',
    title='Original Cell Type Annotation',
    ax=ax,
    show=False
)

plt.tight_layout()
umap_file = fig_dir / f"rare_cells_original_celltype.{figure_format}"
plt.savefig(umap_file, dpi=dpi, bbox_inches='tight')
plt.close()
log_msg(f"  ✓ 保存: {umap_file}")

# 6.2 新聚类UMAP
log_msg(f"\n生成新聚类UMAP...")

n_res = len(leiden_resolutions)
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
axes = axes.flatten()

for i, res in enumerate(leiden_resolutions):
    key = f"leiden_rare_res{res}"
    
    sc.pl.umap(
        adata_rare,
        color=key,
        title=f'Leiden Clustering (res={res})',
        ax=axes[i],
        show=False,
        legend_loc='right margin'
    )

plt.tight_layout()
cluster_umap_file = fig_dir / f"rare_cells_new_clusters.{figure_format}"
plt.savefig(cluster_umap_file, dpi=dpi, bbox_inches='tight')
plt.close()
log_msg(f"  ✓ 保存: {cluster_umap_file}")

# 6.3 按Dataset和Tissue的UMAP
log_msg(f"\n生成Dataset和Tissue UMAP...")

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

sc.pl.umap(
    adata_rare,
    color=batch_key,
    title='Dataset',
    ax=axes[0],
    show=False,
    legend_loc='right margin'
)

sc.pl.umap(
    adata_rare,
    color=tissue_key,
    title='Tissue',
    ax=axes[1],
    show=False,
    legend_loc='right margin'
)

plt.tight_layout()
meta_umap_file = fig_dir / f"rare_cells_metadata.{figure_format}"
plt.savefig(meta_umap_file, dpi=dpi, bbox_inches='tight')
plt.close()
log_msg(f"  ✓ 保存: {meta_umap_file}")

# ==================== Step 7: Marker基因分析 ====================

log_step(7, "Marker Gene Analysis")

# 使用默认分辨率进行marker分析
cluster_key_rare = f"leiden_rare_res{default_resolution}"

log_msg(f"\n使用聚类: {cluster_key_rare}")
log_msg(f"  聚类数: {adata_rare.obs[cluster_key_rare].nunique()}")

# 7.1 FindAllMarkers
log_msg(f"\n运行FindAllMarkers...")

sc.tl.rank_genes_groups(
    adata_rare,
    groupby=cluster_key_rare,
    method='wilcoxon',
    key_added='rank_genes_groups'
)

# 提取marker基因
log_msg(f"\n提取top marker基因...")

all_markers = []
for cluster in adata_rare.obs[cluster_key_rare].unique():
    markers_cluster = sc.get.rank_genes_groups_df(
        adata_rare,
        group=cluster,
        key='rank_genes_groups'
    )
    
    # 添加cluster信息
    markers_cluster['cluster'] = cluster
    
    # 过滤
    markers_cluster = markers_cluster[
        (markers_cluster['pvals_adj'] < 0.05) &
        (markers_cluster['logfoldchanges'] > marker_logfc_threshold)
    ]
    
    all_markers.append(markers_cluster)

markers_df = pd.concat(all_markers, ignore_index=True)

log_msg(f"  共识别 {len(markers_df)} 个marker基因")

# 保存marker结果
markers_file = output_path / "rare_cells_markers.csv"
markers_df.to_csv(markers_file, index=False)
log_msg(f"  ✓ 保存: {markers_file}")

# 显示每个cluster的top marker
log_msg(f"\n每个cluster的top 5 marker基因:")
for cluster in sorted(adata_rare.obs[cluster_key_rare].unique()):
    cluster_markers = markers_df[markers_df['cluster'] == cluster].head(5)
    genes_str = ', '.join(cluster_markers['names'].tolist())
    log_msg(f"  Cluster {cluster}: {genes_str}")

# 7.2 核心Marker基因DotPlot
log_msg(f"\n生成核心marker基因DotPlot...")

# 收集所有marker基因
all_celltype_markers = []
for markers in celltype_markers.values():
    all_celltype_markers.extend(markers)
all_celltype_markers = list(set(all_celltype_markers))

# 过滤存在的基因
available_markers = [g for g in all_celltype_markers if g in adata_rare.var_names]
log_msg(f"  可用marker基因: {len(available_markers)}/{len(all_celltype_markers)}")

if len(available_markers) > 0:
    # 按细胞类型分组marker
    marker_dict_filtered = {}
    for ct, markers in celltype_markers.items():
        filtered = [g for g in markers if g in adata_rare.var_names]
        if filtered:
            marker_dict_filtered[ct] = filtered
    
    # 生成DotPlot
    fig, ax = plt.subplots(figsize=(max(12, len(available_markers)*0.4), 8))
    
    sc.pl.dotplot(
        adata_rare,
        var_names=marker_dict_filtered,
        groupby=cluster_key_rare,
        dendrogram=True,
        standard_scale='var',
        ax=ax,
        show=False
    )
    
    plt.tight_layout()
    dotplot_file = fig_dir / f"rare_cells_markers_dotplot.{figure_format}"
    plt.savefig(dotplot_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    log_msg(f"  ✓ 保存: {dotplot_file}")

# ==================== Step 8: 聚类-细胞类型交叉分析 ====================

log_step(8, "Cluster vs Original CellType Analysis")

# 创建交叉表
log_msg(f"\n创建聚类×原始细胞类型交叉表...")

cross_tab = pd.crosstab(
    adata_rare.obs[cluster_key_rare],
    adata_rare.obs[cell_type_key],
    margins=True
)

log_msg(f"  维度: {cross_tab.shape}")

# 保存交叉表
cross_file = output_path / "cluster_vs_celltype_crosstab.csv"
cross_tab.to_csv(cross_file)
log_msg(f"  ✓ 保存: {cross_file}")

# 显示交叉表
log_msg(f"\n聚类×原始细胞类型交叉表:")
log_msg(cross_tab.to_string())

# 可视化交叉表
log_msg(f"\n生成交叉表热图...")

fig, ax = plt.subplots(figsize=(12, 8))

# 去掉Total行列
plot_data = cross_tab.iloc[:-1, :-1]

sns.heatmap(
    plot_data,
    annot=True,
    fmt='g',
    cmap='YlOrRd',
    cbar_kws={'label': 'Cell Count'},
    ax=ax
)

plt.title('Cluster vs Original Cell Type', fontsize=14, fontweight='bold')
plt.xlabel('Original Cell Type', fontsize=12)
plt.ylabel(f'New Cluster ({cluster_key_rare})', fontsize=12)
plt.tight_layout()

heatmap_file = fig_dir / f"cluster_celltype_crosstab_heatmap.{figure_format}"
plt.savefig(heatmap_file, dpi=dpi, bbox_inches='tight')
plt.close()
log_msg(f"  ✓ 保存: {heatmap_file}")

# ==================== Step 9: 应用手动注释或生成模板 ====================

log_step(9, "Applying Manual Annotation")

annotation_path = output_path / annotation_file
has_annotation = False

# 为每个聚类生成参考信息(用于辅助手动注释)
cluster_reference = []

for cluster in sorted(adata_rare.obs[cluster_key_rare].unique()):
    cluster_cells = adata_rare.obs[adata_rare.obs[cluster_key_rare] == cluster]
    celltype_dist = cluster_cells[cell_type_key].value_counts()
    total_cells = len(cluster_cells)
    
    # 获取top marker
    cluster_markers = markers_df[markers_df['cluster'] == cluster].head(10)
    top_markers = ', '.join(cluster_markers['names'].tolist())
    
    cluster_reference.append({
        'Cluster': cluster,
        'N_Cells': total_cells,
        'Original_Types': ', '.join([f"{ct}:{cnt}" for ct, cnt in celltype_dist.items()]),
        'Top_Markers': top_markers
    })

reference_df = pd.DataFrame(cluster_reference)

# 保存参考信息
reference_file = output_path / "cluster_reference_for_annotation.csv"
reference_df.to_csv(reference_file, index=False)
log_msg(f"\n✓ 聚类参考信息已保存: {reference_file}")

# 模式检查
if RUN_MODE == "pre_annotation":
    log_msg(f"\n📋 Pre-annotation模式: 仅生成annotation模板")
    log_msg(f"   Annotation文件检查: 跳过")
    
    # 生成模板
    template_path = output_path / "Annotation_template.csv"
    clusters = sorted(adata_rare.obs[cluster_key_rare].unique())
    template_df = pd.DataFrame({
        'Cluster': clusters,
        'CellType': ['Unknown'] * len(clusters)
    })
    template_df.to_csv(template_path, index=False, header=False)
    
    log_msg(f"\n✓ Annotation模板已创建: {template_path}")
    log_msg(f"✓ 聚类总数: {len(clusters)}")
    log_msg(f"\n📝 下一步操作:")
    log_msg(f"   1. 查看DotPlot: {fig_dir / f'rare_cells_markers_dotplot.{figure_format}'}")
    log_msg(f"   2. 查看参考信息: {reference_file}")
    log_msg(f"   3. 编辑 {template_path} 添加细胞类型注释")
    log_msg(f"   4. 重命名为 {annotation_file}")
    log_msg(f"   5. 运行脚本,设置 RUN_MODE = 'post_annotation'")
    
    # 显示聚类参考信息帮助注释
    log_msg(f"\n聚类参考信息(辅助注释):")
    for _, row in reference_df.iterrows():
        log_msg(f"\n  Cluster {row['Cluster']} ({row['N_Cells']} cells):")
        log_msg(f"    原始类型: {row['Original_Types']}")
        log_msg(f"    Top markers: {row['Top_Markers']}")

elif not annotation_path.exists():
    log_msg(f"\n❌ Annotation文件未找到: {annotation_file}")
    
    if RUN_MODE == "post_annotation":
        log_msg(f"\n❌ 错误: post_annotation模式需要Annotation.csv")
        log_msg(f"   请先创建 {annotation_file} 或切换到 'pre_annotation' 模式")
        sys.exit(1)
    
    log_msg(f"   跳过annotation步骤")
    log_msg(f"\n请创建annotation文件,格式:")
    log_msg(f"   第1列: Cluster编号")
    log_msg(f"   第2列: 细胞类型名称")
    log_msg(f"   不要表头")

else:
    # 读取annotation文件
    log_msg(f"\n读取annotation文件: {annotation_file}")
    annotations = pd.read_csv(annotation_path, header=None, names=['Cluster', 'CellType'])
    
    log_msg(f"   加载了 {len(annotations)} 个注释")
    
    # 显示annotation内容
    log_msg(f"\n   Annotation映射:")
    for _, row in annotations.iterrows():
        log_msg(f"   Cluster {row['Cluster']} -> {row['CellType']}")
    
    # 创建映射字典
    annotation_dict = dict(zip(annotations['Cluster'].astype(str), annotations['CellType']))
    
    # 应用annotation
    adata_rare.obs['manual_celltype'] = adata_rare.obs[cluster_key_rare].astype(str).map(annotation_dict)
    
    # 检查未匹配的cluster
    unmatched = adata_rare.obs['manual_celltype'].isna().sum()
    if unmatched > 0:
        log_msg(f"\n   ⚠️  警告: {unmatched} 个细胞没有匹配的annotation")
        unmatched_clusters = adata_rare.obs[adata_rare.obs['manual_celltype'].isna()][cluster_key_rare].unique()
        log_msg(f"   未匹配的cluster: {', '.join(map(str, unmatched_clusters))}")
    
    # 显示细胞类型分布
    log_msg(f"\n   手动注释后的细胞类型分布:")
    manual_celltype_counts = adata_rare.obs['manual_celltype'].value_counts()
    for celltype, count in manual_celltype_counts.items():
        pct = count / adata_rare.n_obs * 100
        log_msg(f"   {celltype}: {count} cells ({pct:.1f}%)")
    
    has_annotation = True
    log_msg(f"\n✓ 手动annotation应用成功")

# ==================== Step 10: 生成细胞级注释表 ====================

log_step(10, "Generating Cell-Level Annotation Table")

# 模式检查
if RUN_MODE == "pre_annotation":
    log_msg(f"\n跳过细胞级注释表(pre_annotation模式)")
    log_msg(f"   请在应用手动annotation后运行此步骤")
    
elif has_annotation:
    log_msg(f"\n为每个细胞生成注释信息...")
    
    # 为每个细胞添加注释信息
    cell_annotation_table = []
    
    for cell_id in adata_rare.obs.index:
        cell_obs = adata_rare.obs.loc[cell_id]
        
        cell_annotation_table.append({
            'Cell_ID': cell_id,
            'Original_CellType': cell_obs[cell_type_key],
            'Original_Cluster': cell_obs[cluster_key],
            'New_Cluster': cell_obs[cluster_key_rare],
            'Manual_CellType': cell_obs['manual_celltype'],
            'Dataset': cell_obs[batch_key],
            'Tissue': cell_obs[tissue_key],
            'Annotation_Changed': cell_obs[cell_type_key] != cell_obs['manual_celltype']
        })
    
    cell_annotation_df = pd.DataFrame(cell_annotation_table)
    
    # 保存细胞级注释
    cell_annotation_file = output_path / "cell_level_annotation.csv"
    cell_annotation_df.to_csv(cell_annotation_file, index=False)
    log_msg(f"  ✓ 保存: {cell_annotation_file}")
    
    # 统计修改情况
    n_changed = cell_annotation_df['Annotation_Changed'].sum()
    pct_changed = n_changed / len(cell_annotation_df) * 100
    
    log_msg(f"\n注释修改统计:")
    log_msg(f"  修改的细胞: {n_changed}/{len(cell_annotation_df)} ({pct_changed:.1f}%)")
    
    if n_changed > 0:
        log_msg(f"\n按 原始类型 → 新类型 统计:")
        change_summary = cell_annotation_df[cell_annotation_df['Annotation_Changed']].groupby(
            ['Original_CellType', 'Manual_CellType']
        ).size().reset_index(name='N_Cells')
        
        for _, row in change_summary.iterrows():
            log_msg(f"  {row['Original_CellType']} → {row['Manual_CellType']}: {row['N_Cells']} cells")
    
    # 生成新旧注释对比UMAP
    log_msg(f"\n生成注释对比UMAP...")
    
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    
    # 原始注释
    celltypes_orig = adata_rare.obs[cell_type_key].unique()
    colors_orig = [celltype_colors.get(ct, '#808080') for ct in celltypes_orig]
    
    sc.pl.umap(
        adata_rare,
        color=cell_type_key,
        palette=colors_orig,
        title='Original Annotation',
        ax=axes[0],
        show=False,
        legend_loc='right margin'
    )
    
    # 手动注释
    celltypes_manual = adata_rare.obs['manual_celltype'].unique()
    colors_manual = [celltype_colors.get(ct, '#808080') for ct in celltypes_manual]
    
    sc.pl.umap(
        adata_rare,
        color='manual_celltype',
        palette=colors_manual,
        title='Manual Annotation',
        ax=axes[1],
        show=False,
        legend_loc='right margin'
    )
    
    plt.tight_layout()
    compare_umap_file = fig_dir / f"annotation_comparison.{figure_format}"
    plt.savefig(compare_umap_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    log_msg(f"  ✓ 保存: {compare_umap_file}")
else:
    log_msg(f"\n跳过细胞级注释表(无手动annotation)")
    log_msg(f"   请创建 {annotation_file} 并重新运行")

# ==================== Step 11: 保存数据 ====================

log_step(11, "Saving Processed Data")

if RUN_MODE == "pre_annotation":
    # pre_annotation模式: 保存带有新聚类的数据
    log_msg(f"\n保存带有新聚类的数据...")
    output_h5ad = output_path / "rare_cells_with_new_clusters.h5ad"
    adata_rare.write_h5ad(output_h5ad, compression='gzip')
    
    file_size = output_h5ad.stat().st_size / (1024**2)
    log_msg(f"✓ 保存: {output_h5ad} ({file_size:.2f} MB)")
    log_msg(f"   包含: 新聚类结果, 用于手动注释")

elif has_annotation:
    # post_annotation模式: 保存带有手动注释的数据
    log_msg(f"\n保存带有手动注释的数据...")
    output_h5ad = output_path / "rare_cells_manually_annotated.h5ad"
    adata_rare.write_h5ad(output_h5ad, compression='gzip')
    
    file_size = output_h5ad.stat().st_size / (1024**2)
    log_msg(f"✓ 保存: {output_h5ad} ({file_size:.2f} MB)")
    log_msg(f"   包含: 手动注释 (manual_celltype列)")

else:
    # full模式但没有annotation
    log_msg(f"\n保存分析结果...")
    output_h5ad = output_path / "rare_cells_analyzed.h5ad"
    adata_rare.write_h5ad(output_h5ad, compression='gzip')
    
    file_size = output_h5ad.stat().st_size / (1024**2)
    log_msg(f"✓ 保存: {output_h5ad} ({file_size:.2f} MB)")

# ==================== Step 12: 生成总结报告 ====================

log_step(12, "Generating Summary Report")

from datetime import datetime

report = []
report.append("="*70)
report.append("Rare Cells Subset Analysis - Summary Report")
report.append("="*70)
report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report.append(f"Run Mode: {RUN_MODE.upper()}")
report.append("")

report.append("[ Analysis Parameters ]")
report.append(f"  Rare cell threshold: < {max_cells_threshold} cells per dataset")
report.append(f"  Clustering resolution: {default_resolution}")
report.append(f"  N_top_genes: {n_top_genes}")
report.append("")

report.append("[ Rare Cell Statistics ]")
report.append(f"  Total rare cells: {adata_rare.n_obs:,}")
report.append(f"  Rare combinations: {len(rare_df)}")
report.append(f"  Cell types involved: {adata_rare.obs[cell_type_key].nunique()}")
report.append(f"  Datasets involved: {adata_rare.obs[batch_key].nunique()}")
report.append(f"  Tissues involved: {adata_rare.obs[tissue_key].nunique()}")
report.append("")

report.append("[ Cell Type Distribution in Subset ]")
for celltype, count in subset_celltype_counts.items():
    pct = count / adata_rare.n_obs * 100
    report.append(f"  {celltype}: {count} cells ({pct:.1f}%)")
report.append("")

report.append("[ New Clustering Results ]")
report.append(f"  Clustering method: Leiden")
report.append(f"  Default resolution: {default_resolution}")
report.append(f"  Number of clusters: {adata_rare.obs[cluster_key_rare].nunique()}")
report.append("")

if RUN_MODE == "pre_annotation":
    report.append("[ Pre-Annotation Mode ]")
    report.append(f"  Annotation template created: Yes")
    report.append(f"  DotPlot generated: Yes")
    report.append(f"  Cluster reference saved: Yes")
    report.append("")
    report.append("[ Next Steps ]")
    report.append(f"  1. Review DotPlot: figures/rare_cells_markers_dotplot.{figure_format}")
    report.append(f"  2. Review reference: cluster_reference_for_annotation.csv")
    report.append(f"  3. Edit Annotation_template.csv")
    report.append(f"  4. Rename to {annotation_file}")
    report.append(f"  5. Run with RUN_MODE='post_annotation'")
    report.append("")

elif has_annotation:
    report.append("[ Manual Annotation Applied ]")
    manual_celltype_counts = adata_rare.obs['manual_celltype'].value_counts()
    for celltype, count in manual_celltype_counts.items():
        pct = count / adata_rare.n_obs * 100
        report.append(f"  {celltype}: {count} cells ({pct:.1f}%)")
    report.append("")
    
    if 'cell_annotation_df' in locals():
        n_changed = cell_annotation_df['Annotation_Changed'].sum()
        pct_changed = n_changed / len(cell_annotation_df) * 100
        report.append("[ Annotation Changes ]")
        report.append(f"  Cells with changed annotation: {n_changed}/{len(cell_annotation_df)} ({pct_changed:.1f}%)")
        report.append("")
        
        if n_changed > 0:
            report.append("[ Change Summary ]")
            change_summary = cell_annotation_df[cell_annotation_df['Annotation_Changed']].groupby(
                ['Original_CellType', 'Manual_CellType']
            ).size().reset_index(name='N_Cells')
            for _, row in change_summary.iterrows():
                report.append(f"  {row['Original_CellType']} → {row['Manual_CellType']}: {row['N_Cells']} cells")
            report.append("")

report.append("[ Output Files ]")
report.append(f"  - Dataset counts: dataset_celltype_counts.csv")
report.append(f"  - Rare combinations: rare_combinations.csv")
report.append(f"  - Rare cells info: rare_cells_info.csv")
report.append(f"  - Marker genes: rare_cells_markers.csv")
report.append(f"  - Cluster crosstab: cluster_vs_celltype_crosstab.csv")
report.append(f"  - Cluster reference: cluster_reference_for_annotation.csv")

if RUN_MODE == "pre_annotation":
    report.append(f"  - Annotation template: Annotation_template.csv")
    report.append(f"  - Data with clusters: rare_cells_with_new_clusters.h5ad")
elif has_annotation:
    report.append(f"  - Cell annotation table: cell_level_annotation.csv")
    report.append(f"  - Annotated data: rare_cells_manually_annotated.h5ad")

report.append(f"  - Figures: figures/")
report.append("")

report.append("="*70)

# 打印和保存报告
report_text = "\n".join(report)
print("\n" + report_text)

report_file = output_path / "summary_report.txt"
with open(report_file, 'w') as f:
    f.write(report_text)

log_msg(f"\n✓ 总结报告已保存: {report_file}")

if RUN_MODE == "pre_annotation":
    log_msg(f"\n{'='*70}")
    log_msg("Pre-annotation阶段完成!")
    log_msg(f"{'='*70}")
    log_msg(f"\n📋 下一步: 手动注释聚类")
    log_msg(f"   1. 打开DotPlot查看marker表达")
    log_msg(f"   2. 参考 cluster_reference_for_annotation.csv")
    log_msg(f"   3. 编辑 Annotation_template.csv")
    log_msg(f"   4. 重命名为 {annotation_file}")
    log_msg(f"   5. 设置 RUN_MODE = 'post_annotation' 并重新运行")
elif has_annotation:
    log_msg(f"\n{'='*70}")
    log_msg("手动注释已应用!")
    log_msg(f"{'='*70}")
    log_msg(f"\n✓ 可以使用 cell_level_annotation.csv 更新主数据集")
else:
    log_msg(f"\n{'='*70}")
    log_msg("分析完成!")
    log_msg(f"{'='*70}")

log_msg(f"\n所有结果保存在: {output_dir}")