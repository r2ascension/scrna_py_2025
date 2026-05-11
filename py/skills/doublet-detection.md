---
name: doublet-detection
description: Identify doublets using CellTypist multi-lineage predictions and perform comprehensive quality control for technical artifacts.
license: Complete terms in LICENSE.txt
metadata:
  title: Doublet Detection and Quality Control
  version: "1.1"
  author: Single-cell Analysis Pipeline Team
  category: quality-control
  difficulty: beginner
  estimated_time: "30-60 minutes"
  tags:
    - doublet-detection
    - quality-control
    - celltypist
    - multi-evidence-gating
    - lineage-analysis
    - data-cleaning
  prerequisites:
    - Log-normalized expression matrix (adata.X)
    - "Batch identifiers in obs['dataset']"
    - CellTypist annotations or ability to run CellTypist
    - UMAP coordinates for visualization (optional)
  outputs:
    - adata_doublet_filtered.h5ad with doublets removed
    - Doublet scores and multi-evidence metrics
    - Quality control metric distributions
    - UMAP visualizations of doublet scores
    - Statistical summaries by batch
  related_skills:
    - allcells-integration
    - celltype-specific-analysis
    - subcluster-analysis
---

# Skill: Doublet Detection and Quality Control

## 功能描述

基于CellTypist多谱系预测识别双细胞(doublets)，并进行全面的质量控制，去除低质量细胞和技术伪影。

## 核心技术

- **Top2谱系逻辑**: 基于第二高谱系概率识别双细胞
- **Margin分析**: 评估谱系预测的确定性
- **多证据门控**: 结合邻域纯度、映射置信度等多个指标
- **质量指标**: MT%, stress score, cell cycle, doublet score

## 适用场景

- [OK] 高细胞密度区域 (易产生双细胞)
- [OK] 解离后数据集 (双细胞率通常5-10%)
- [OK] 混合谱系信号的细胞 (如上皮+免疫标记)
- [OK] 质量控制流程的一部分
- [NO] 不适用于真实的过渡状态细胞 (需谨慎区分)

## 核心概念

### 什么是双细胞(Doublet)?
双细胞是两个细胞被错误地识别为一个细胞，导致:
- 混合的基因表达谱
- 多个谱系的标记基因共表达
- 异常高的基因/UMI计数

### 双细胞 vs 过渡状态
**关键区别**:
- **双细胞**: 两个不相关谱系的混合 (如T细胞+成纤维细胞)
- **过渡状态**: 同一谱系内的中间状态 (如Naive B → Memory B)

**检测策略**: 使用谱系级别而非细胞类型级别

## 推荐脚本

### 双细胞检测
**脚本**: `celltypist_doublet_detection_20260114_v1_1.py` (v1.1 [OK])

**特性**:
- Top2谱系 + margin逻辑
- 多证据门控
- 保留过渡状态
- 批次效应考虑

### 数据清洗
**脚本**: `remove_doublets_and_gse299751_20260114.py` (v1.0)

**功能**:
- 移除检测到的双细胞
- 移除特定数据集 (如GSE299751)
- 质量指标过滤

### 质量检查
**脚本**: `check_data_quality.py`

**功能**:
- 全面质量指标计算
- 批次效应检测
- 异常值识别

**脚本**: `quick_data_check_detailed_20260127.py` (v1.0)

**功能**:
- 快速质量概览
- 统计摘要
- 可视化报告

## 输入要求

### 数据格式
```python
# AnnData对象
adata.X                        # 表达矩阵 (log-normalized)
adata.layers['counts']         # 原始counts (可选)
adata.obs['dataset']           # 批次标识
adata.obsm['X_umap']           # UMAP坐标 (可选，用于可视化)
```

### CellTypist注释
如果已有CellTypist结果:
```python
adata.obs['celltypist_pred']           # 预测标签
adata.obs['celltypist_conf_max']       # 最大置信度
adata.uns['celltypist_prob_matrix']    # 概率矩阵 (关键!)
```

如果没有，脚本会自动运行CellTypist。

## 输出结果

### 文件结构
```
/home/h2048/data/py/{YYMMDD}/doublet_detection/
├── adata_doublet_filtered.h5ad        # 过滤后数据
├── doublet_scores.csv                 # 双细胞评分
├── figures/
│   ├── doublet_score_distribution.pdf
│   ├── doublet_umap.pdf
│   ├── top2_lineage_scatter.pdf
│   ├── margin_distribution.pdf
│   └── qc_metrics_violin.pdf
├── statistics/
│   ├── doublet_summary.txt
│   ├── removed_cells_by_batch.csv
│   └── lineage_confusion_matrix.csv
└── README.md
```

### 新增注释
```python
# 双细胞评分
adata.obs['doublet_score']             # 综合双细胞评分
adata.obs['top1_lineage']              # 第一高谱系
adata.obs['top2_lineage']              # 第二高谱系
adata.obs['top1_prob']                 # 第一高概率
adata.obs['top2_prob']                 # 第二高概率
adata.obs['lineage_margin']            # top1 - top2
adata.obs['is_doublet_candidate']      # 双细胞候选 (bool)

# 质量指标
adata.obs['pct_counts_mt']             # 线粒体百分比
adata.obs['stress_score']              # 应激评分
adata.obs['neighborhood_purity']       # 邻域纯度
```

## 工作流程

### Step 1: 运行CellTypist (如果需要)
```python
import celltypist
import scanpy as sc

# 加载数据
adata = sc.read_h5ad("adata_input.h5ad")

# 运行CellTypist
model = celltypist.models.Model.load(
    model='Human_Lung_Atlas.pkl'
)

predictions = celltypist.annotate(
    adata,
    model=model,
    majority_voting=True
)

# 获取概率矩阵 (关键!)
prob_matrix = predictions.probability_matrix
```

### Step 2: 定义谱系分组
```python
# 将细胞类型分组到主要谱系
LINEAGE_GROUPS = {
    'Epithelial': [
        'Basal', 'Suprabasal', 'Ciliated', 'Secretory_Goblet',
        'Secretory_Club', 'AT1', 'AT2', 'Ionocyte', 'Brush'
    ],
    'T_NK': [
        'CD4_T', 'CD8_T', 'Treg', 'NK', 'NKT', 'MAIT', 'ILC'
    ],
    'B_Plasma': [
        'B_naive', 'B_memory', 'Plasma', 'Plasmablast', 'GC_B'
    ],
    'Myeloid': [
        'Macrophage_alveolar', 'Macrophage_interstitial',
        'Monocyte_classical', 'Monocyte_nonclassical',
        'DC1', 'DC2', 'pDC', 'Mast'
    ],
    'Stromal': [
        'Fibroblast_adventitial', 'Fibroblast_alveolar',
        'Myofibroblast', 'SMC_airway', 'SMC_vascular',
        'Pericyte', 'Mesothelial'
    ],
    'Endothelial': [
        'Endothelial_arterial', 'Endothelial_venous',
        'Endothelial_capillary', 'Endothelial_lymphatic'
    ]
}

# 创建反向映射
celltype_to_lineage = {}
for lineage, celltypes in LINEAGE_GROUPS.items():
    for ct in celltypes:
        celltype_to_lineage[ct] = lineage
```

### Step 3: 计算谱系概率
```python
# 将细胞类型概率聚合到谱系级别
lineage_probs = prob_matrix.copy()
lineage_probs.columns = lineage_probs.columns.map(
    lambda x: celltype_to_lineage.get(x, 'Other')
)
lineage_probs = lineage_probs.groupby(level=0, axis=1).sum()

# 计算top1和top2谱系
top1_prob = lineage_probs.max(axis=1)
top2_prob = lineage_probs.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
margin = top1_prob - top2_prob

top1_lineage = lineage_probs.idxmax(axis=1)
top2_lineage = lineage_probs.apply(
    lambda x: x.nlargest(2).index[1], axis=1
)

# 添加到adata
adata.obs['top1_lineage'] = top1_lineage.values
adata.obs['top2_lineage'] = top2_lineage.values
adata.obs['top1_prob'] = top1_prob.values
adata.obs['top2_prob'] = top2_prob.values
adata.obs['lineage_margin'] = margin.values
```

### Step 4: 双细胞检测 (基础版)
```python
# 基础双细胞检测逻辑
TOP2_THRESHOLD = 0.5        # 第二高谱系阈值
MARGIN_THRESHOLD = 0.15     # margin阈值

is_doublet_basic = (
    (adata.obs['top2_prob'] > TOP2_THRESHOLD) &
    (adata.obs['lineage_margin'] < MARGIN_THRESHOLD)
)

adata.obs['is_doublet_basic'] = is_doublet_basic

print(f"Basic doublet detection: {is_doublet_basic.sum()} cells "
      f"({is_doublet_basic.sum()/adata.n_obs*100:.2f}%)")
```

### Step 5: 多证据门控 (推荐)
```python
# 计算邻域纯度
from sklearn.neighbors import NearestNeighbors

nn = NearestNeighbors(n_neighbors=30, metric='euclidean')
nn.fit(adata.obsm['X_scanvi'])  # 或 X_scvi
_, indices = nn.kneighbors(adata.obsm['X_scanvi'])

purity = []
for i, neighbors in enumerate(indices):
    my_lineage = adata.obs['top1_lineage'].iloc[i]
    neighbor_lineages = adata.obs['top1_lineage'].iloc[neighbors]
    purity.append((neighbor_lineages == my_lineage).mean())

adata.obs['neighborhood_purity'] = purity

# 多证据评分
adata.obs['doublet_score'] = (
    (adata.obs['top2_prob'] > 0.5).astype(int) +           # 高top2
    (adata.obs['lineage_margin'] < 0.15).astype(int) +     # 低margin
    (adata.obs['neighborhood_purity'] < 0.5).astype(int)   # 低纯度
)

# 只移除2+个红旗的细胞
is_doublet_multi = adata.obs['doublet_score'] >= 2

print(f"Multi-evidence doublet detection: {is_doublet_multi.sum()} cells "
      f"({is_doublet_multi.sum()/adata.n_obs*100:.2f}%)")
```

### Step 6: 质量指标过滤
```python
# 计算质量指标
# MT%
mt_genes = adata.var_names.str.startswith('MT-')
adata.obs['pct_counts_mt'] = (
    adata[:, mt_genes].X.sum(1).A1 / adata.X.sum(1).A1 * 100
)

# Stress score
STRESS_GENES = ['HSPA1A', 'HSPA1B', 'DNAJB1', 'FOS', 'JUN', 'EGR1']
stress_genes = [g for g in STRESS_GENES if g in adata.var_names]
sc.tl.score_genes(adata, stress_genes, score_name='stress_score')

# 质量过滤阈值
MT_THRESHOLD = 20
STRESS_PERCENTILE = 95

low_quality = (
    (adata.obs['pct_counts_mt'] > MT_THRESHOLD) |
    (adata.obs['stress_score'] > adata.obs['stress_score'].quantile(0.95))
)

print(f"Low quality cells: {low_quality.sum()} cells "
      f"({low_quality.sum()/adata.n_obs*100:.2f}%)")
```

### Step 7: 综合过滤
```python
# 综合过滤: 双细胞 OR 低质量
cells_to_remove = is_doublet_multi | low_quality

print(f"\nTotal cells to remove: {cells_to_remove.sum()} "
      f"({cells_to_remove.sum()/adata.n_obs*100:.2f}%)")
print(f"  - Doublets: {is_doublet_multi.sum()}")
print(f"  - Low quality: {low_quality.sum()}")
print(f"  - Overlap: {(is_doublet_multi & low_quality).sum()}")

# 过滤
adata_filtered = adata[~cells_to_remove].copy()

print(f"\nCells remaining: {adata_filtered.n_obs} "
      f"({adata_filtered.n_obs/adata.n_obs*100:.2f}%)")
```

### Step 8: 可视化
```python
# Top2 vs Margin散点图
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 散点图
ax = axes[0]
scatter = ax.scatter(
    adata.obs['lineage_margin'],
    adata.obs['top2_prob'],
    c=adata.obs['doublet_score'],
    cmap='RdYlBu_r',
    s=1,
    alpha=0.5
)
ax.axhline(y=0.5, color='red', linestyle='--', label='Top2 threshold')
ax.axvline(x=0.15, color='red', linestyle='--', label='Margin threshold')
ax.set_xlabel('Lineage Margin (Top1 - Top2)')
ax.set_ylabel('Top2 Lineage Probability')
ax.set_title('Doublet Detection')
ax.legend()
plt.colorbar(scatter, ax=ax, label='Doublet Score')

# UMAP
ax = axes[1]
sc.pl.umap(
    adata,
    color='doublet_score',
    cmap='RdYlBu_r',
    ax=ax,
    show=False
)
ax.set_title('Doublet Score on UMAP')

plt.tight_layout()
plt.savefig('figures/doublet_detection_overview.pdf')
```

## 使用示例

### 基本用法
```bash
# 运行双细胞检测
python celltypist_doublet_detection_20260114_v1_1.py \
    --input adata_allcells_FINAL.h5ad \
    --output /home/h2048/data/py/20260203/doublet_detection \
    --top2-threshold 0.5 \
    --margin-threshold 0.15

# 预期输出:
# - 双细胞比例: 5-10% (正常范围)
# - 过滤后数据: adata_doublet_filtered.h5ad
```

### 移除双细胞和特定数据集
```bash
# 同时移除双细胞和GSE299751数据集
python remove_doublets_and_gse299751_20260114.py \
    --input adata_allcells_FINAL.h5ad \
    --output adata_cleaned.h5ad \
    --remove-dataset GSE299751
```

### 快速质量检查
```bash
# 生成质量报告
python quick_data_check_detailed_20260127.py \
    --input adata_allcells_FINAL.h5ad \
    --output qc_report.html
```

## 关键参数

### 双细胞检测阈值
```python
TOP2_THRESHOLD = 0.5        # 第二高谱系概率阈值
MARGIN_THRESHOLD = 0.15     # Top1-Top2 margin阈值
PURITY_THRESHOLD = 0.5      # 邻域纯度阈值
MIN_EVIDENCE = 2            # 最少证据数 (0-3)
```

### 质量过滤阈值
```python
MT_THRESHOLD = 20           # 线粒体百分比上限
STRESS_PERCENTILE = 95      # 应激评分百分位数
MIN_GENES = 200             # 最少基因数
MAX_GENES = 6000            # 最多基因数 (可能是双细胞)
MIN_COUNTS = 500            # 最少UMI数
```

### 邻域分析参数
```python
N_NEIGHBORS = 30            # 邻居数量
METRIC = 'euclidean'        # 距离度量
```

## 质量检查

### 运行前
- [ ] 确认CellTypist已运行或概率矩阵可用
- [ ] 检查谱系分组定义是否完整
- [ ] 验证UMAP坐标存在 (用于可视化)

### 运行后
- [ ] 双细胞比例合理 (5-15%)
- [ ] 双细胞在UMAP上位于谱系边界
- [ ] 移除的细胞MT%或stress score高
- [ ] 保留的细胞质量指标正常
- [ ] 批次间双细胞比例相似

## 常见问题

### Q1: 双细胞比例过高 (>20%)
**可能原因**:
- 阈值过于宽松
- 数据质量问题
- 真实的过渡状态被误判

**解决方案**:
```python
# 提高阈值
TOP2_THRESHOLD = 0.6  # 从0.5提高到0.6
MARGIN_THRESHOLD = 0.1  # 从0.15降低到0.1

# 或增加最少证据数
MIN_EVIDENCE = 3  # 从2提高到3
```

### Q2: 双细胞比例过低 (<2%)
**可能原因**:
- 阈值过于严格
- 数据质量很好
- 谱系分组过于粗糙

**解决方案**:
```python
# 降低阈值
TOP2_THRESHOLD = 0.4
MARGIN_THRESHOLD = 0.2

# 或使用更细的谱系分组
```

### Q3: 过渡状态被误判为双细胞
**识别**:
- 查看被标记为双细胞的细胞的top1和top2谱系
- 如果是相关谱系 (如Naive B和Memory B)，可能是过渡状态

**解决方案**:
```python
# 定义相关谱系对 (不应被标记为双细胞)
RELATED_LINEAGES = [
    ('B_Plasma', 'B_Plasma'),  # B细胞内部过渡
    ('T_NK', 'T_NK'),          # T/NK细胞内部过渡
]

# 在过滤时排除相关谱系对
for i, row in adata.obs.iterrows():
    if (row['top1_lineage'], row['top2_lineage']) in RELATED_LINEAGES:
        adata.obs.loc[i, 'is_doublet_candidate'] = False
```

### Q4: 某些批次双细胞比例异常高
**分析**:
```python
# 按批次统计双细胞比例
doublet_by_batch = adata.obs.groupby('dataset')['is_doublet_candidate'].mean()
print(doublet_by_batch.sort_values(ascending=False))
```

**可能原因**:
- 该批次细胞密度高
- 解离条件不同
- 批次效应导致误判

**解决方案**:
- 检查该批次的质量指标
- 考虑批次特异性阈值
- 检查批次校正是否充分

### Q5: 邻域纯度计算失败
**原因**: 缺少潜在空间表示
**解决方案**:
```python
# 确保有scVI或scANVI潜在空间
if 'X_scanvi' not in adata.obsm:
    if 'X_scvi' not in adata.obsm:
        # 使用PCA作为替代
        sc.pp.pca(adata, n_comps=50)
        representation = adata.obsm['X_pca']
    else:
        representation = adata.obsm['X_scvi']
else:
    representation = adata.obsm['X_scanvi']

# 使用该表示计算邻域纯度
nn.fit(representation)
```

## 高级用法

### 批次特异性阈值
```python
# 为每个批次计算自适应阈值
for batch in adata.obs['dataset'].unique():
    batch_mask = adata.obs['dataset'] == batch
    batch_data = adata[batch_mask]

    # 使用该批次的95th percentile作为阈值
    top2_thresh = batch_data.obs['top2_prob'].quantile(0.95)
    margin_thresh = batch_data.obs['lineage_margin'].quantile(0.05)

    # 应用批次特异性阈值
    is_doublet_batch = (
        (batch_data.obs['top2_prob'] > top2_thresh) &
        (batch_data.obs['lineage_margin'] < margin_thresh)
    )

    adata.obs.loc[batch_mask, 'is_doublet_adaptive'] = is_doublet_batch
```

### 与计算方法比较
```python
# 使用Scrublet作为对照
import scrublet as scr

scrub = scr.Scrublet(adata.X)
doublet_scores_scrublet, predicted_doublets_scrublet = scrub.scrub_doublets()

adata.obs['doublet_score_scrublet'] = doublet_scores_scrublet
adata.obs['is_doublet_scrublet'] = predicted_doublets_scrublet

# 比较两种方法
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(
    adata.obs['is_doublet_candidate'],
    adata.obs['is_doublet_scrublet']
)
print("Confusion matrix (CellTypist vs Scrublet):")
print(cm)
```

### 双细胞模拟验证
```python
# 人工创建双细胞验证检测性能
import numpy as np

# 随机选择两个不同谱系的细胞
lineages = adata.obs['top1_lineage'].unique()
n_synthetic = 1000

synthetic_doublets = []
for _ in range(n_synthetic):
    # 随机选择两个不同谱系
    lin1, lin2 = np.random.choice(lineages, 2, replace=False)

    # 随机选择每个谱系的一个细胞
    cell1 = adata[adata.obs['top1_lineage'] == lin1].obs_names[0]
    cell2 = adata[adata.obs['top1_lineage'] == lin2].obs_names[0]

    # 混合表达谱 (50:50)
    mixed_expr = (adata[cell1].X + adata[cell2].X) / 2

    synthetic_doublets.append(mixed_expr)

# 在合成双细胞上测试检测性能
# ...
```

## 最佳实践

1. [OK] **使用谱系级别** - 不要用细胞类型级别 (避免误判过渡状态)
2. [OK] **多证据门控** - 不要只依赖单一指标
3. [OK] **可视化验证** - 检查双细胞在UMAP上的位置
4. [OK] **批次效应检查** - 确保双细胞比例在批次间一致
5. [OK] **保守过滤** - 宁可保留可疑细胞，也不要过度过滤
6. [OK] **记录统计** - 保存移除细胞的详细信息
7. [OK] **与其他方法比较** - 如Scrublet, DoubletFinder

## 相关Skills

- `allcells-integration` - 全细胞整合 (下游)
- `celltype-specific-analysis` - 细胞类型分析 (下游)
- `subcluster-analysis` - 亚群聚类 (下游)

## 参考文献

- CellTypist: Domínguez Conde et al., Science 2022
- Scrublet: Wolock et al., Cell Systems 2019
- DoubletFinder: McGinnis et al., Cell Systems 2019
- Doublet detection review: Xi & Li, Genome Biology 2021
