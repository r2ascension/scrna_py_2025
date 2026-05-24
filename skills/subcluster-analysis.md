---
name: subcluster-analysis
description: Fine-grained Leiden clustering within cell types to identify functional subpopulations while preserving BBKNN batch correction.
license: Complete terms in LICENSE.txt
metadata:
  title: Subcluster Analysis with BBKNN Integration
  version: "4.1"
  author: Single-cell Analysis Pipeline Team
  category: analysis
  difficulty: intermediate
  estimated_time: "30-60 minutes"
  tags:
    - leiden-clustering
    - bbknn
    - marker-genes
    - subclustering
    - batch-correction
    - wilcoxon-test
  prerequisites:
    - Completed cell type annotation with BBKNN integration
    - "BBKNN neighbor graph exists (uns['neighbors'])"
    - UMAP coordinates available
    - "Cell type labels in obs['cell_type']"
  outputs:
    - adata_{celltype}_subclustered.h5ad with subcluster labels
    - Filtered marker gene lists per subcluster
    - UMAP visualizations of subclusters
    - Batch effect quality control reports
    - Marker gene heatmaps and dotplots
  related_skills:
    - celltype-specific-analysis
    - cnmf-analysis
    - allcells-integration
---

# Skill: Subcluster Analysis

## 功能描述

在细胞类型内部进行精细聚类，识别功能亚群和细胞状态，同时保留全局BBKNN批次校正效果。

## 核心技术

- **Leiden聚类**: 使用 `restrict_to` 参数保留BBKNN图
- **标记基因分析**: Wilcoxon秩和检验
- **技术伪影过滤**: 去除MT/Ribo/IEG/应激基因
- **批次效应监控**: 确保亚群不由批次驱动

## 适用场景

- [OK] 已完成细胞类型注释，需要识别功能亚群
- [OK] 细胞类型内部异质性高
- [OK] 需要识别疾病相关细胞状态
- [OK] 需要保留批次校正效果
- [NO] 不适用于初始聚类 (应先用全细胞整合)

## 核心设计原则

### 关键原则1: 保留BBKNN图
**不要重建neighbors图!** 使用 `restrict_to` 参数:

```python
# [NO] WRONG - 重建图会丢失批次校正
adata_subset = adata[adata.obs['cell_type'] == 'Basal'].copy()
sc.pp.neighbors(adata_subset)  # 重建图 - 批次效应回来了!
sc.tl.leiden(adata_subset)

# [OK] CORRECT - 使用restrict_to保留BBKNN图
sc.tl.leiden(
    adata,  # 使用完整数据集
    restrict_to=('cell_type', ['Basal']),  # 只在Basal细胞内聚类
    resolution=0.2,
    key_added='leiden_basal_subcluster',
    neighbors_key='neighbors'  # 使用现有BBKNN图
)
```

### 关键原则2: 固定低分辨率
**避免过度分裂**: 使用保守的分辨率参数

```python
# 推荐分辨率范围
RESOLUTION_RANGE = [0.1, 0.2, 0.3]  # 不要超过0.5

# 原因:
# - 高分辨率容易产生批次驱动的伪亚群
# - 低分辨率更稳定，识别真实生物学差异
# - 可以通过多个分辨率测试找到最优值
```

### 关键原则3: 技术基因过滤
**标记基因分析前必须过滤**:

```python
# 技术伪影基因列表
TECHNICAL_GENES = {
    'mt': r'^MT-',                    # 线粒体基因
    'ribo': r'^RP[SL]|^MRP[SL]',     # 核糖体基因
    'histone': r'^H[1234]|^HIST',    # 组蛋白基因
    'ieg': ['FOS', 'JUN', 'EGR1', 'EGR2'],  # 即早基因
    'stress': ['HSPA1A', 'HSPA1B', 'DNAJB1'],  # 应激基因
    'unannotated': r'^(LINC|AC\d+|AL\d+|RP11-|CTD-)'  # 未注释基因
}

# 过滤函数
def filter_technical_genes(marker_df):
    mask = ~(
        marker_df['names'].str.match(r'^MT-') |
        marker_df['names'].str.match(r'^RP[SL]') |
        marker_df['names'].isin(['FOS', 'JUN', 'EGR1'])
    )
    return marker_df[mask]
```

## 推荐脚本

### 通用流程 (推荐)
**脚本**: `universal_celltype_subcluster_pipeline_20260121_v4_1.py` (v4.1 [OK])

**特性**:
- 参数化细胞类型输入
- 自动化标记基因分析
- 批次效应检测
- 多分辨率测试

### 细胞类型特异性脚本

| 细胞类型 | 脚本 | 版本 |
|---------|------|------|
| Epithelial | `epithelial_subcluster_bbknn_pipeline_20260122_v4_5_2.py` | v4.5.2 [OK] |
| B cells | `bcell_subcluster_analysis_v2_20260119.py` | v2.0 |
| Myeloid | `myeloid_subcluster_analysis_v2_20260121.py` | v2.0 |
| Stromal | `stromal_subcluster_analysis_v2_20260121.py` | v2.0 |
| T/NK | `tcell_subcluster_analysis_20260126_v1_0.py` | v1.0 |
| T/NK (统一) | `tnk_subcluster_unified_20260129_v2_0.py` | v2.0 |

## 输入要求

### 数据来源
通常从细胞类型特异性分析结果:
```python
# 从细胞类型分析结果加载
adata = sc.read_h5ad("adata_epithelial_FINAL.h5ad")

# 必需字段
adata.obs['cell_type']         # 细胞类型标签
adata.uns['neighbors']         # BBKNN邻居图 (关键!)
adata.obsm['X_umap']           # UMAP坐标
```

### 必需的图结构
```python
# 检查BBKNN图是否存在
assert 'neighbors' in adata.uns, "BBKNN graph not found!"
assert 'connectivities' in adata.obsp, "Connectivity matrix missing!"
```

## 输出结果

### 文件结构
```
/home/h2048/data/py/{YYMMDD}/{celltype}_subcluster/
├── adata_{celltype}_subclustered.h5ad
├── figures/
│   ├── umap_subclusters_by_celltype.pdf
│   ├── umap_batch_check.pdf
│   ├── marker_heatmap_filtered.pdf
│   ├── marker_dotplot_top10.pdf
│   └── resolution_comparison.pdf
├── markers/
│   ├── {celltype}_subcluster_markers.csv
│   └── {celltype}_subcluster_markers_filtered.csv
└── README.md
```

### 新增注释
```python
# 亚群标签 (多分辨率)
adata.obs['leiden_basal_subcluster_0.1']
adata.obs['leiden_basal_subcluster_0.2']
adata.obs['leiden_basal_subcluster_0.3']

# 标记基因结果
adata.uns['rank_genes_groups_basal']
```

## 工作流程

### Step 1: 多分辨率聚类
```python
# 测试多个分辨率
for celltype in celltypes:
    for res in [0.1, 0.2, 0.3]:
        sc.tl.leiden(
            adata,
            restrict_to=('cell_type', [celltype]),
            resolution=res,
            key_added=f'leiden_{celltype}_res{res}',
            neighbors_key='neighbors'  # 使用BBKNN图
        )
```

### Step 2: 标记基因分析
```python
# 对每个细胞类型的亚群
sc.tl.rank_genes_groups(
    adata,
    groupby=f'leiden_{celltype}_subcluster',
    method='wilcoxon',
    use_raw=True,  # 或False，取决于数据结构
    pts=True,      # 计算表达百分比
    key_added=f'rank_genes_{celltype}'
)

# 提取结果
marker_df = sc.get.rank_genes_groups_df(
    adata,
    group=None,
    key=f'rank_genes_{celltype}'
)
```

### Step 3: 技术基因过滤
```python
# 过滤技术伪影
marker_df_filtered = marker_df[
    ~marker_df['names'].str.match(r'^MT-') &
    ~marker_df['names'].str.match(r'^RP[SL]') &
    ~marker_df['names'].str.match(r'^MRP[SL]') &
    ~marker_df['names'].str.match(r'^H[1234]') &
    ~marker_df['names'].isin(['FOS', 'JUN', 'EGR1', 'EGR2', 'EGR3']) &
    ~marker_df['names'].isin(['HSPA1A', 'HSPA1B', 'DNAJB1']) &
    ~marker_df['names'].str.match(r'^(LINC|AC\d+|AL\d+|RP11-)')
]

# 保存过滤后的标记基因
marker_df_filtered.to_csv(f"markers/{celltype}_markers_filtered.csv")
```

### Step 4: 批次效应检查
```python
# 可视化批次分布
sc.pl.umap(
    adata,
    color=[f'leiden_{celltype}_subcluster', 'batch', 'dataset'],
    ncols=3,
    save=f'_{celltype}_batch_check.pdf'
)

# 统计检验
from scipy.stats import chi2_contingency
contingency = pd.crosstab(
    adata.obs[f'leiden_{celltype}_subcluster'],
    adata.obs['batch']
)
chi2, p_value, dof, expected = chi2_contingency(contingency)

if p_value < 0.01:
    print(f"WARNING: Subclusters may be batch-driven (p={p_value:.2e})")
```

### Step 5: 可视化
```python
# 在原始UMAP上可视化亚群
sc.pl.umap(
    adata,
    color=f'leiden_{celltype}_subcluster',
    legend_loc='on data',
    save=f'_{celltype}_subclusters.pdf'
)

# 标记基因热图 (过滤后)
sc.pl.rank_genes_groups_heatmap(
    adata,
    n_genes=10,
    groupby=f'leiden_{celltype}_subcluster',
    key=f'rank_genes_{celltype}',
    show_gene_labels=True,
    save=f'_{celltype}_markers_heatmap.pdf'
)

# 标记基因点图
top_markers = marker_df_filtered.groupby('group').head(5)['names'].unique()
sc.pl.dotplot(
    adata,
    var_names=top_markers,
    groupby=f'leiden_{celltype}_subcluster',
    save=f'_{celltype}_markers_dotplot.pdf'
)
```

## 使用示例

### 基本用法 (通用流程)
```bash
# 运行通用亚群分析
python universal_celltype_subcluster_pipeline_20260121_v4_1.py \
    --input adata_epithelial_FINAL.h5ad \
    --celltype Basal \
    --resolutions 0.1,0.2,0.3 \
    --output /home/h2048/data/py/20260203/epithelial_subcluster
```

### 上皮细胞亚群分析
```bash
# 使用上皮专用脚本
python epithelial_subcluster_bbknn_pipeline_20260122_v4_5_2.py

# 预期亚群 (Basal细胞示例):
# - Basal_Proliferating (MKI67+, TOP2A+)
# - Basal_Suprabasal_transition (KRT13+, KRT4+)
# - Basal_Stress (HSPA1A+, FOS+) - 可能需要过滤
```

### B细胞亚群分析
```bash
python bcell_subcluster_analysis_v2_20260119.py

# 预期亚群:
# - Naive_B (IGHD+, IGHM+, TCL1A+)
# - Memory_IgM (CD27+, IGHM+)
# - Memory_IgG (CD27+, IGHG1+)
# - Memory_IgA (CD27+, IGHA1+)
# - Plasma_cells (MZB1+, SDC1+, JCHAIN+)
# - Plasmablasts (MZB1+, XBP1+, MKI67+)
```

### 髓系细胞亚群分析
```bash
python myeloid_subcluster_analysis_v2_20260121.py

# 预期亚群 (Macrophage示例):
# - AM_FABP4 (FABP4+, MARCO+)
# - AM_PPARG (PPARG+, APOE+)
# - IM_CCL (CCL2+, CCL3+, CCL4+)
# - IM_CHIT1 (CHIT1+, CHI3L1+)
# - IM_CX3CR1 (CX3CR1+, FCGR3A+)
```

## 关键参数

### Leiden聚类参数
```python
RESOLUTION_RANGE = [0.1, 0.2, 0.3]  # 保守范围
N_NEIGHBORS = 30                     # 使用BBKNN的默认值
MIN_CLUSTER_SIZE = 50                # 最小亚群大小
```

### 标记基因参数
```python
METHOD = 'wilcoxon'                  # 推荐方法
N_GENES = 100                        # 每个亚群的标记基因数
MIN_LOGFOLDCHANGE = 0.5              # 最小log fold change
MIN_PCT_DIFF = 0.1                   # 最小表达百分比差异
```

### 过滤阈值
```python
# 标记基因过滤
MIN_PVAL_ADJ = 0.05                  # 校正后p值
MIN_LOGFC = 0.5                      # 最小log fold change
MIN_PCT = 0.2                        # 最小表达百分比

# 批次效应检测
BATCH_PVAL_THRESHOLD = 0.01          # 批次关联p值阈值
```

## 质量检查清单

### 运行前
- [ ] 确认BBKNN图存在 (`'neighbors' in adata.uns`)
- [ ] 检查细胞类型标签无缺失
- [ ] 验证UMAP坐标存在
- [ ] 确认批次标识可用

### 运行后
- [ ] 亚群数量合理 (3-10个，不超过15个)
- [ ] 批次效应检验p值 >0.01
- [ ] 标记基因生物学意义明确
- [ ] 无技术伪影基因在top标记中
- [ ] UMAP可视化亚群边界清晰
- [ ] 亚群在批次间分布均匀

## 常见问题

### Q1: 亚群过多 (>15个)
**原因**: 分辨率过高或批次效应
**解决方案**:
```python
# 降低分辨率
resolution = 0.1  # 从0.3降到0.1

# 检查批次效应
sc.pl.umap(adata, color=['leiden_subcluster', 'batch'])
```

### Q2: 亚群由批次驱动
**检测**:
```python
# 卡方检验
from scipy.stats import chi2_contingency
contingency = pd.crosstab(adata.obs['leiden_subcluster'], adata.obs['batch'])
chi2, p_value, _, _ = chi2_contingency(contingency)
print(f"Batch association p-value: {p_value:.2e}")
```

**解决方案**:
- 确认使用了 `neighbors_key='neighbors'` (BBKNN图)
- 降低分辨率
- 检查BBKNN参数是否合适

### Q3: 标记基因都是技术伪影
**原因**: 未过滤技术基因
**解决方案**: 应用技术基因过滤 (见Step 3)

### Q4: 某些亚群细胞数极少 (<50 cells)
**解决方案**:
```python
# 合并小亚群
cluster_sizes = adata.obs['leiden_subcluster'].value_counts()
small_clusters = cluster_sizes[cluster_sizes < 50].index

# 重新标记为"Other"或合并到最近的大亚群
adata.obs['leiden_subcluster_filtered'] = adata.obs['leiden_subcluster'].copy()
adata.obs.loc[
    adata.obs['leiden_subcluster'].isin(small_clusters),
    'leiden_subcluster_filtered'
] = 'Other'
```

### Q5: 不同分辨率结果差异大
**分析**:
```python
# 比较不同分辨率
for res in [0.1, 0.2, 0.3]:
    n_clusters = adata.obs[f'leiden_res{res}'].nunique()
    print(f"Resolution {res}: {n_clusters} clusters")

# 可视化比较
sc.pl.umap(
    adata,
    color=['leiden_res0.1', 'leiden_res0.2', 'leiden_res0.3'],
    ncols=3
)
```

**选择标准**:
- 生物学可解释性 (标记基因是否有意义)
- 亚群数量合理 (不过多也不过少)
- 批次效应最小
- 与已知生物学一致

## 高级用法

### 层级聚类
```python
# 先粗聚类，再在每个粗聚类内细聚类
# Step 1: 粗聚类
sc.tl.leiden(
    adata,
    restrict_to=('cell_type', ['Macrophage']),
    resolution=0.1,
    key_added='leiden_macro_coarse'
)

# Step 2: 在每个粗聚类内细聚类
for coarse_cluster in adata.obs['leiden_macro_coarse'].unique():
    sc.tl.leiden(
        adata,
        restrict_to=('leiden_macro_coarse', [coarse_cluster]),
        resolution=0.2,
        key_added=f'leiden_macro_{coarse_cluster}_fine'
    )
```

### 标记基因富集分析
```python
# 使用gseapy进行通路富集
import gseapy as gp

for cluster in adata.obs['leiden_subcluster'].unique():
    marker_genes = marker_df_filtered[
        marker_df_filtered['group'] == cluster
    ]['names'].head(100).tolist()

    enr = gp.enrichr(
        gene_list=marker_genes,
        gene_sets=['GO_Biological_Process_2021', 'KEGG_2021_Human'],
        organism='human'
    )

    enr.results.to_csv(f"enrichment/cluster_{cluster}_enrichment.csv")
```

### 亚群稳定性测试
```python
# Bootstrap重采样测试亚群稳定性
from sklearn.metrics import adjusted_rand_score

ari_scores = []
for i in range(10):
    # 随机采样80%细胞
    sample_idx = np.random.choice(
        adata.n_obs,
        size=int(adata.n_obs * 0.8),
        replace=False
    )
    adata_sample = adata[sample_idx].copy()

    # 重新聚类
    sc.tl.leiden(
        adata_sample,
        resolution=0.2,
        key_added='leiden_bootstrap'
    )

    # 计算ARI
    ari = adjusted_rand_score(
        adata.obs.loc[sample_idx, 'leiden_subcluster'],
        adata_sample.obs['leiden_bootstrap']
    )
    ari_scores.append(ari)

print(f"Mean ARI: {np.mean(ari_scores):.3f} ± {np.std(ari_scores):.3f}")
# ARI >0.8 表示稳定
```

## 最佳实践总结

1. [OK] **始终使用 `restrict_to`** - 保留BBKNN批次校正
2. [OK] **使用低分辨率** - 0.1-0.3范围，避免过度分裂
3. [OK] **过滤技术基因** - MT/Ribo/IEG/Stress基因
4. [OK] **检查批次效应** - 卡方检验 + UMAP可视化
5. [OK] **在原始UMAP可视化** - 不要重新计算UMAP
6. [OK] **验证标记基因** - 与已知生物学一致
7. [OK] **测试多个分辨率** - 选择最优解
8. [OK] **记录参数** - 便于重现和比较

## 相关Skills

- `celltype-specific-analysis` - 细胞类型分析 (上游)
- `cnmf-analysis` - 基因表达程序 (互补)
- `allcells-integration` - 全细胞整合 (上游)

## 参考文献

- Leiden algorithm: Traag et al., Scientific Reports 2019
- BBKNN: Polański et al., Bioinformatics 2020
- Marker gene analysis: Soneson & Robinson, Nature Methods 2018
