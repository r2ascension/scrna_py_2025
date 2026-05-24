---
name: celltype-specific-analysis
description: Identify cell subtypes and states within pre-classified populations using scVI-CellTypist-scANVI deep learning pipeline.
license: Complete terms in LICENSE.txt
metadata:
  title: Cell Type-Specific Deep Analysis
  version: "3.5.1"
  author: Single-cell Analysis Pipeline Team
  category: analysis
  difficulty: intermediate
  estimated_time: "2-4 hours"
  tags:
    - scvi
    - celltypist
    - scanvi
    - cell-type-annotation
    - subtype-identification
    - immune-cells
    - deep-learning
  prerequisites:
    - Completed all-cells integration with cell type labels
    - "Raw counts layer available (layers['counts'])"
    - Minimum 5k cells per cell type (10k+ recommended)
    - "Batch identifiers in obs['dataset']"
    - GPU access for model training
  outputs:
    - adata_{celltype}_FINAL.h5ad with subtype annotations
    - Trained scVI and scANVI models
    - HVG gene list for model reuse
    - UMAP visualizations of subtypes
    - Marker gene analysis results
  related_skills:
    - allcells-integration
    - subcluster-analysis
    - cnmf-analysis
    - doublet-detection
---
# Skill: Cell-Type Specific Analysis

## 功能描述

针对特定细胞类型进行深度分析，使用scVI-CellTypist-scANVI流程识别细胞亚型和状态。适用于已完成初步分类的细胞类型子集。

## 核心技术栈

- **scVI**: 细胞类型内批次校正
- **CellTypist**: 细胞亚型自动注释 (专用模型)
- **scANVI**: 亚型标签精炼
- **P0修复系列**: v3.5.x生产级bug修复

## 适用场景

- [OK] 已完成全细胞整合，需要细胞类型深度分析
- [OK] 细胞类型数量 >10k cells (推荐)
- [OK] 需要识别细胞亚型和功能状态
- [OK] 免疫细胞精细分型 (T/B/Myeloid)
- [OK] 基质/血管细胞异质性分析

## 支持的细胞类型

### 1. T/NK细胞
**脚本**: `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` (v3.5.1 [OK])

**CellTypist模型**: `Immune_All_Low.pkl`

**预期亚型**:
- CD4+ T cells: Naive, Central Memory (CM), Effector Memory (EM), Tissue-Resident Memory (TRM)
- CD8+ T cells: EM, EMRA, TRM, Exhausted
- Regulatory T cells (Tregs)
- NK cells, NKT cells, MAIT cells
- ILC (Innate Lymphoid Cells)

**特异性标记**:
```python
CD4_MARKERS = ['CD4', 'IL7R', 'TCF7', 'LEF1']
CD8_MARKERS = ['CD8A', 'CD8B', 'GZMK', 'GZMB']
TREG_MARKERS = ['FOXP3', 'IL2RA', 'CTLA4']
NK_MARKERS = ['NCAM1', 'KLRD1', 'NKG7', 'GNLY']
```

### 2. B细胞
**脚本**: `bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` (v3.5.1 [OK])

**CellTypist模型**: `Immune_All_Low.pkl`

**预期亚型**:
- Naive B cells
- Memory B cells: IgM+, IgG+, IgA+
- Plasma cells, Plasmablasts
- Germinal Center (GC) B cells
- Regulatory B cells (Bregs)

**特异性标记**:
```python
B_MARKERS = ['CD19', 'MS4A1', 'CD79A', 'CD79B']
NAIVE_MARKERS = ['IGHD', 'IGHM', 'TCL1A']
MEMORY_MARKERS = ['CD27', 'IGHG1', 'IGHG3', 'IGHA1']
PLASMA_MARKERS = ['MZB1', 'SDC1', 'JCHAIN', 'XBP1']
```

**v3.5.2新增**: Counts恢复功能 (从外部文件恢复丢失的counts层)

### 3. 髓系细胞
**脚本**: `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` (v3.5.1 [OK])

**CellTypist模型**: `Immune_All_Low.pkl`

**预期亚型**:
- Alveolar Macrophages (AM)
- Interstitial Macrophages (IM)
- Monocyte-derived Macrophages
- Classical Monocytes (CD14+)
- Non-classical Monocytes (CD16+)
- Dendritic Cells: DC1, DC2, pDC, Activated DC
- Mast Cells

**特异性标记**:
```python
MACRO_MARKERS = ['CD68', 'CD163', 'MSR1', 'MRC1']
AM_MARKERS = ['FABP4', 'MARCO', 'PPARG']
MONO_MARKERS = ['CD14', 'FCGR3A', 'S100A8', 'S100A9']
DC_MARKERS = ['CLEC9A', 'CD1C', 'CLEC10A', 'LILRA4']
MAST_MARKERS = ['TPSAB1', 'CPA3', 'KIT']
```

### 4. 基质/血管细胞
**脚本**: `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py` (v3.5.1 [OK])

**CellTypist模型**: `Human_Lung_Atlas.pkl`

**预期亚型**:
- **内皮细胞**: Arterial, Venous (pulmonary/systemic), Lymphatic, Capillary
- **成纤维细胞**: Adventitial, Alveolar, Peribronchial, Myofibroblast
- **平滑肌细胞**: Airway SMC, Arterial SMC, Perivascular
- **周细胞** (Pericytes)
- **间皮细胞** (Mesothelial)
- **软骨细胞** (Chondrocytes)

**特异性标记**:
```python
ENDO_MARKERS = ['PECAM1', 'CDH5', 'VWF', 'KDR']
FIBRO_MARKERS = ['COL1A1', 'COL1A2', 'PDGFRA', 'DCN']
SMC_MARKERS = ['ACTA2', 'MYH11', 'TAGLN', 'CNN1']
PERICYTE_MARKERS = ['RGS5', 'PDGFRB', 'CSPG4']
```

### 5. 上皮细胞
**脚本**: `epithelial_scvi_scanvi_20260118_v2_6.py` (v2.6)

**CellTypist模型**: `Human_Lung_Atlas.pkl`

**预期亚型**:
- Basal cells, Suprabasal, Dividing Basal
- Ciliated cells, Deuterosome cells
- Secretory cells: Goblet, Club
- AT1, AT2 (alveolar epithelial)
- Submucosal Gland (SMG) cells
- Ionocytes, Brush cells

**特殊架构**: 双scANVI训练 (major lineage + fine-grained)

## v3.5.x系列P0修复 (Critical Bug Fixes)

### P0-1: scANVI专用adata_model
**问题**: 标签索引不对齐导致训练错误
**修复**:
```python
# WRONG (v3.4及之前)
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,  # 使用完整adata
    labels_key='cell_type'
)

# CORRECT (v3.5.1+)
adata_model_scanvi = sc.AnnData(
    X=adata_model.X,
    obs=adata.obs[['cell_type']].copy(),  # 确保索引对齐
    var=adata_model.var.copy()
)
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata_model_scanvi,
    labels_key='cell_type'
)
```

### P0-2: predict()方法dtype修复
**问题**: `predict(soft=True).max()` 返回非标量导致崩溃
**修复**:
```python
# WRONG
conf = predictions.predict(soft=True).max(axis=1)

# CORRECT
conf = np.asarray(predictions.predict(soft=True).max(axis=1)).flatten()
```

### P0-3: 多数投票整合
**问题**: CellTypist的 `majority_voting` 结果未正确合并
**修复**:
```python
# 确保使用majority_voting结果
predictions = celltypist.annotate(
    adata,
    model=model,
    majority_voting=True  # 必须启用
)
adata.obs['celltypist_majority'] = predictions.predicted_labels['majority_voting']
```

### P0-4: 特征验证
**问题**: CellTypist基因重叠率低时静默失败
**修复**:
```python
# 硬性检查基因重叠
overlap_genes = set(adata.var_names) & set(model.features)
overlap_ratio = len(overlap_genes) / len(model.features)

if overlap_ratio < 0.5:
    raise ValueError(f"Gene overlap too low: {overlap_ratio:.2%}")
elif overlap_ratio < 0.7:
    print(f"WARNING: Gene overlap suboptimal: {overlap_ratio:.2%}")
```

### P0-5: 预训练HVG加载修复
**问题**: 使用 `symbol_base` 导致基因列表不匹配
**修复**:
```python
# WRONG
hvg_mask = adata.var['symbol_base'].isin(hvg_genes)

# CORRECT
hvg_mask = adata.var_names.isin(hvg_genes)
```

### P1-6: 内存优化
**改进**: 避免完整X矩阵复制
```python
# 直接使用layers，不复制到X
adata_model = sc.AnnData(
    X=adata.layers['counts'][:, hvg_mask],  # 直接引用
    obs=adata.obs[['dataset']].copy(),
    var=adata.var.loc[hvg_mask].copy()
)
# 内存减少: 60-70%
```

## 输入要求

### 数据来源
通常从全细胞整合结果中子集化:
```python
# 从全细胞数据中提取特定类型
adata_allcells = sc.read_h5ad("adata_allcells_FINAL.h5ad")
adata_tcell = adata_allcells[adata_allcells.obs['cell_type'] == 'T'].copy()
```

### 必需字段
```python
adata.layers['counts']     # 原始counts (整数)
adata.obs['dataset']       # 批次标识
adata.obs['cell_type']     # 粗分类标签 (可选)
```

### 推荐细胞数
- **最小**: 5k cells
- **推荐**: 10k-100k cells
- **大规模**: >100k cells (启用HVG优化)

## 输出结果

### 文件结构
```
/home/h2048/data/py/{YYMMDD}/{celltype}_analysis/
├── adata_{celltype}_FINAL.h5ad
├── models/
│   ├── scvi_model/
│   └── scanvi_model/
├── hvg_genes.txt
├── figures/
│   ├── umap_scvi_batch.pdf
│   ├── umap_scanvi_subtype.pdf
│   ├── marker_dotplot.pdf
│   └── confidence_distribution.pdf
└── README.md
```

### 新增注释
```python
adata.obs['subtype_celltypist']    # CellTypist亚型预测
adata.obs['subtype_scanvi']        # scANVI精炼亚型
adata.obs['subtype_final']         # 最终推荐亚型
adata.obs['subtype_confidence']    # 预测置信度
```

## 关键参数

### scVI参数 (细胞类型特异性)
```python
N_LATENT = 75               # 潜在维度 (比全细胞小)
N_LAYERS = 3                # 网络层数
N_HVG = 3000                # HVG数量 (细胞类型特异性)
MAX_EPOCHS = 800            # 训练轮数 (更多轮数)
```

### scANVI参数
```python
MAX_EPOCHS = 600            # 训练轮数
MIN_CELLS_PER_TYPE = 10     # 稀有亚型阈值
N_SAMPLES_PER_LABEL = 50    # 类别平衡
```

## 使用示例

### T细胞分析
```bash
# 1. 准备输入数据 (从全细胞结果子集化)
python -c "
import scanpy as sc
adata = sc.read_h5ad('adata_allcells_FINAL.h5ad')
adata_t = adata[adata.obs['cell_type'].isin(['T', 'NK'])].copy()
adata_t.write_h5ad('adata_tcell_input.h5ad')
"

# 2. 运行T细胞分析
python t_scvi_celltypist_scanvi_20260108_v3_5_1.py

# 3. 检查结果
python -c "
import scanpy as sc
adata = sc.read_h5ad('adata_tcell_FINAL.h5ad')
print(adata.obs['subtype_final'].value_counts())
"
```

### B细胞分析 (含counts恢复)
```bash
# 如果counts层丢失，使用v3.5.2的恢复功能
python bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py \
    --recover-counts \
    --external-file /path/to/original_with_counts.h5ad
```

## 质量检查

### 运行前
- [ ] 细胞数 >5k
- [ ] `layers['counts']` 存在
- [ ] 批次标识无缺失
- [ ] 基因名为symbols格式

### 运行后
- [ ] CellTypist基因重叠 >70%
- [ ] scVI/scANVI训练收敛
- [ ] 亚型数量合理 (5-20个)
- [ ] 标记基因表达与亚型一致
- [ ] 批次混合良好 (UMAP可视化)

## 常见问题

### Q1: 亚型过度分裂 (>30个亚型)
**原因**: 批次效应未完全校正
**解决方案**:
- 增加scVI训练轮数
- 检查批次标识是否正确
- 考虑使用更保守的scANVI参数

### Q2: 某些亚型细胞数极少 (<10 cells)
**解决方案**: 脚本会自动合并到"Unknown"

### Q3: counts层丢失
**解决方案** (v3.5.2):
```python
# 使用counts恢复功能
adata = recover_counts_from_external(
    adata,
    external_file="/path/to/original.h5ad"
)
```

### Q4: 免疫细胞注释不准确
**检查**:
- 确认使用 `Immune_All_Low.pkl` 模型
- 检查基因重叠率
- 验证标记基因表达

## 下游分析

完成细胞类型特异性分析后:

1. **亚群聚类**: `universal_celltype_subcluster_pipeline_20260121_v4_1.py`
2. **cNMF分析**: `quick_start_cnmf_analysis_20251225_{CellType}_v1_1.py`
3. **差异表达**: 使用scanpy的 `sc.tl.rank_genes_groups()`
4. **轨迹推断**: 使用scVelo或CellRank

## 批量运行

```bash
# 并行运行多个细胞类型
parallel -j 4 python {} ::: \
    t_scvi_celltypist_scanvi_20260108_v3_5_1.py \
    bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py \
    myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py \
    stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py
```

## 版本对比

| 版本 | 关键特性 | 推荐使用 |
|------|---------|---------|
| v3.5.1 | P0修复完整版 | [OK] 生产环境 |
| v3.4 | 稳定版 (有P0 bugs) | [NO] 已弃用 |
| v3.2 | 早期版本 | [NO] 已弃用 |

## 相关Skills

- `allcells-integration` - 全细胞整合 (上游)
- `subcluster-analysis` - 亚群聚类 (下游)
- `cnmf-analysis` - 基因表达程序 (下游)
- `doublet-detection` - 质量控制

## 参考文献

- scVI: Lopez et al., Nature Methods 2018
- CellTypist: Domínguez Conde et al., Science 2022
- scANVI: Xu et al., Molecular Systems Biology 2021
- Human Lung Atlas: Sikkema et al., Nature Medicine 2023
