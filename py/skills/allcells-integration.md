---
name: allcells-integration
description: Integrate multi-batch single-cell data using scVI-CellTypist-scANVI pipeline for batch correction and automated cell type annotation.
license: Complete terms in LICENSE.txt
metadata:
  title: All-Cells Integration Analysis
  version: "2.3.1"
  author: Single-cell Analysis Pipeline Team
  category: integration
  difficulty: intermediate
  estimated_time: "2-3 hours"
  tags:
    - single-cell
    - integration
    - batch-correction
    - scvi
    - celltypist
    - scanvi
  prerequisites:
    - Raw count matrix in AnnData format
    - Batch identifiers
    - GPU availability recommended
  outputs:
    - Batch-corrected embeddings (X_scvi, X_scanvi)
    - Cell type annotations
    - UMAP visualizations
    - Trained models
  related_skills:
    - celltype-specific-analysis
    - doublet-detection
    - subcluster-analysis
---

# Skill: All-Cells Integration Analysis

## 功能描述

对所有细胞类型进行批次校正和自动注释的完整流程，使用scVI-CellTypist-scANVI三阶段深度学习架构。

## 核心技术栈

- **scVI**: 批次效应去除 (Variational Autoencoder)
- **CellTypist**: 机器学习自动注释 (Human_Lung_Atlas模型)
- **scANVI**: 半监督标签精炼 (VAE with classifier)

## 适用场景

- [OK] 新数据集的初始整合分析
- [OK] 多批次/多样本数据的批次校正
- [OK] 需要自动化细胞类型注释
- [OK] 大规模数据集 (>100k cells)
- [OK] 作为下游分析的基础流程

## 推荐脚本

**生产版本**: `allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py` (v2.3.1)

### 版本特性

**v2.3.1 新增功能**:
- [OK] 协变量系统 (MT%, stress score, cell cycle)
- [OK] Unknown细胞质量过滤 (邻域纯度 + 质量门控)
- [OK] 类别不平衡处理 (`n_samples_per_label`)
- [OK] scArches兼容参数
- [OK] 统一n_latent维度 (scVI和scANVI均为100)

## 输入要求

### 数据格式
```python
# AnnData对象要求
adata.X                    # 处理后的表达矩阵 (通常为log-normalized)
adata.layers['counts']     # 原始UMI counts (整数, 必需!)
adata.obs['dataset']       # 批次标识符
adata.obs['sample']        # 样本标识符 (可选)
adata.var_names            # 基因名 (推荐HGNC symbols)
```

### 基因名要求
- **必须**: HGNC gene symbols (如 "CD3D", "CD8A")
- **不推荐**: ENSEMBL IDs (如 "ENSG00000167286")
- **处理**: 脚本会自动调用 `normalize_gene_names()` 进行转换

### 数据规模建议
- **小数据集** (<100k cells): 使用全基因集
- **中等数据集** (100k-400k cells): HVG=4000
- **大数据集** (>400k cells): HVG=4000 + 内存优化模式

## 输出结果

### 主要输出文件

```
/home/h2048/data/py/{YYMMDD}/allcells_integration/
├── adata_allcells_FINAL.h5ad          # 最终整合数据
├── models/
│   ├── scvi_model/                     # scVI模型 (可复用)
│   └── scanvi_model/                   # scANVI模型
├── hvg_genes.txt                       # HVG列表 (模型复用必需!)
├── figures/
│   ├── umap_scvi_batch.pdf            # scVI UMAP (按批次)
│   ├── umap_scanvi_celltype.pdf       # scANVI UMAP (按细胞类型)
│   ├── celltypist_confidence.pdf      # CellTypist置信度分布
│   └── qc_metrics.pdf                 # 质量控制指标
└── README.md                           # 自动生成的分析摘要
```

### AnnData对象新增内容

```python
# 嵌入表示
adata.obsm['X_scvi']          # scVI潜在空间 (批次校正)
adata.obsm['X_scanvi']        # scANVI潜在空间 (标签引导)
adata.obsm['X_umap_scvi']     # scVI UMAP坐标
adata.obsm['X_umap_scanvi']   # scANVI UMAP坐标

# 细胞类型注释
adata.obs['celltypist_pred']           # CellTypist预测
adata.obs['celltypist_conf_max']       # CellTypist置信度
adata.obs['celltypist_majority']       # 多数投票结果
adata.obs['scanvi_pred']               # scANVI预测
adata.obs['scanvi_confidence']         # scANVI置信度
adata.obs['cell_type_final']           # 最终推荐标签

# 质量控制
adata.obs['pct_counts_mt']             # 线粒体基因百分比
adata.obs['stress_score']              # 应激基因评分
adata.obs['S_score']                   # S期评分
adata.obs['G2M_score']                 # G2M期评分
adata.obs['neighborhood_purity']       # 邻域纯度 (Unknown过滤用)

# 图结构
adata.uns['neighbors_scvi']            # scVI邻居图
adata.uns['neighbors_scanvi']          # scANVI邻居图
```

## 关键参数配置

### scVI训练参数
```python
N_LATENT = 100              # 潜在维度 (推荐100)
N_LAYERS = 2                # 网络层数 (推荐2)
N_HVG = 4000                # 高变基因数 (大数据集推荐)
MAX_EPOCHS_SCVI = 400       # 最大训练轮数

# 协变量 (v2.3.1新增)
CONTINUOUS_COVARIATES = [
    'pct_counts_mt',        # 线粒体百分比
    'stress_score',         # 应激评分
    'S_score',              # 细胞周期S期
    'G2M_score'             # 细胞周期G2M期
]
```

### scANVI训练参数
```python
MAX_EPOCHS_SCANVI = 200     # 最大训练轮数
MIN_CELLS_PER_TYPE = 10     # 稀有类型阈值
N_SAMPLES_PER_LABEL = 50    # 类别平衡采样 (v2.3.1)

# Unknown过滤阈值 (v2.3.1)
NEIGHBORHOOD_PURITY_THRESHOLD = 0.5
MT_THRESHOLD = 20
STRESS_PERCENTILE = 95
```

### CellTypist参数
```python
MODEL = 'Human_Lung_Atlas.pkl'
MAJORITY_VOTING = True
MIN_GENE_OVERLAP = 0.7      # 最小基因重叠率
```

## 使用示例

### 基本用法

```bash
# 直接运行脚本
python allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py

# 预期运行时间 (V100 GPU, 200k cells):
# - scVI训练: ~40-60分钟
# - CellTypist注释: ~5-10分钟
# - scANVI训练: ~20-30分钟
# - 总计: ~2-3小时
```

### 修改配置

```python
# 在脚本顶部修改配置变量
INPUT_FILE = "/path/to/your/adata.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/20260203/allcells_integration"
BATCH_KEY = "dataset"       # 批次列名
SAMPLE_KEY = "sample"       # 样本列名

# 内存优化 (大数据集)
USE_HVG_FOR_SCVI = True
N_HVG = 4000

# GPU设置
USE_GPU = True
GPU_DEVICE = 0              # 单GPU强制
```

### 复用已训练模型

```python
# 如果已有训练好的scVI模型
PRETRAINED_SCVI_MODEL = "/path/to/scvi_model"

# 脚本会自动:
# 1. 加载 hvg_genes.txt
# 2. 子集化新数据到相同基因
# 3. 加载模型并跳过训练
```

## 质量检查清单

### 运行前检查
- [ ] 确认 `layers['counts']` 存在且为整数
- [ ] 检查基因名格式 (symbols vs ENSEMBL)
- [ ] 验证批次标识符无缺失值
- [ ] 确认GPU可用 (`torch.cuda.is_available()`)

### 运行后检查
- [ ] CellTypist基因重叠率 >70%
- [ ] scVI训练loss收敛
- [ ] scANVI训练loss收敛
- [ ] UMAP可视化批次混合良好
- [ ] 细胞类型注释合理 (与已知标记基因一致)
- [ ] Unknown细胞比例 <10%

## 常见问题

### Q1: "layers['counts'] not found"
**解决方案**:
```python
# 从raw恢复counts
if adata.raw is not None:
    adata.layers['counts'] = adata.raw.X.copy()
```

### Q2: CellTypist基因重叠率低 (<50%)
**原因**: 基因名格式不匹配 (ENSEMBL vs symbols)
**解决方案**: 脚本会自动调用 `normalize_gene_names()` 转换

### Q3: scANVI训练不稳定
**原因**: 稀有细胞类型 (<10 cells)
**解决方案**: 脚本会自动合并稀有类型到 "Unknown"

### Q4: 内存不足 (OOM)
**解决方案**:
```python
# 启用HVG优化
USE_HVG_FOR_SCVI = True
N_HVG = 4000  # 从58k基因减少到4k
# 预期内存减少: 60-70%
```

### Q5: Unknown细胞过多 (>20%)
**原因**:
- 数据质量问题 (高MT%, 应激)
- 新细胞类型未在参考中
- 批次效应未完全校正

**解决方案**: v2.3.1会自动过滤低质量Unknown细胞

## 下游分析建议

完成全细胞整合后，推荐的下游分析流程:

1. **质量控制**: `celltypist_doublet_detection_20260114_v1_1.py`
2. **细胞类型特异性分析**:
   - T/NK: `t_scvi_celltypist_scanvi_20260108_v3_5_1.py`
   - B cells: `bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py`
   - Myeloid: `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py`
   - Stromal: `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py`
3. **亚群聚类**: `universal_celltype_subcluster_pipeline_20260121_v4_1.py`
4. **基因表达程序**: `label_guided_cnmf_pipeline_20260114_v1_1.py`

## 技术细节

### 三阶段架构原理

```
原始数据 (Raw counts)
    ↓
[Stage 1: scVI]
    - 输入: layers['counts'] (整数)
    - 输出: X_scvi (批次校正潜在空间)
    - 目的: 去除批次效应
    ↓
[Stage 2: CellTypist]
    - 输入: log-normalized counts
    - 输出: predicted_labels + confidence
    - 目的: 自动化细胞类型注释
    ↓
[Stage 3: scANVI]
    - 输入: X_scvi + CellTypist labels
    - 输出: X_scanvi (标签引导潜在空间)
    - 目的: 精炼低置信度预测
    ↓
最终注释 (cell_type_final)
```

## 参考文献

- scVI: Lopez et al., Nature Methods 2018
- CellTypist: Domínguez Conde et al., Science 2022
- scANVI: Xu et al., Molecular Systems Biology 2021
