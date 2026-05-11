---
name: scarches-mapping
description: Map new query data to pre-trained reference models using scArches transfer learning for rapid annotation without retraining.
license: Complete terms in LICENSE.txt
metadata:
  title: scArches Reference-to-Query Mapping
  version: "1.2.1"
  author: Single-cell Analysis Pipeline Team
  category: transfer-learning
  difficulty: advanced
  estimated_time: "30-60 minutes"
  tags:
    - scarches
    - transfer-learning
    - model-mapping
    - label-transfer
    - architecture-surgery
    - reference-mapping
  prerequisites:
    - Pre-trained reference scVI/scANVI model
    - "Query data with raw counts (layers['counts'])"
    - Gene overlap >80% between query and reference
    - Query data similar to reference (same tissue/species)
    - Reference gene list (var_names.csv)
  outputs:
    - adata_query_mapped.h5ad with transferred labels
    - Query-specific fine-tuned model
    - Label transfer confidence scores
    - Joint reference-query UMAP visualizations
    - Mapping quality metrics and outlier detection
  related_skills:
    - allcells-integration
    - celltype-specific-analysis
    - doublet-detection
---

# Skill: scArches Model Transfer Learning

## 功能描述

使用scArches框架将新查询数据映射到已训练的参考模型，实现快速注释和批次校正，无需重新训练完整模型。

## 核心技术

- **scArches**: 架构手术(Architecture Surgery)框架
- **Transfer Learning**: 迁移学习
- **Surgical Fine-tuning**: 手术式微调
- **Label Transfer**: 标签迁移

## 适用场景

- [OK] 有高质量参考数据集和训练好的模型
- [OK] 新数据集需要快速注释
- [OK] 新数据集与参考数据集相似 (同组织/物种)
- [OK] 计算资源有限 (比从头训练快10-100x)
- [OK] 需要保持与参考数据的一致性
- [NO] 新数据集与参考差异很大 (不同组织/物种)

## 核心概念

### 什么是scArches?
scArches允许将新数据"映射"到已有的参考模型:

```
参考数据 (Reference)
    ↓ 训练
参考模型 (scVI/scANVI)
    ↓ 架构手术
查询模型 (Query-specific)
    ↓ 微调 (少量epochs)
映射结果 (Query in reference space)
```

### 优势
- **速度快**: 微调只需10-50 epochs (vs 400-800 epochs从头训练)
- **一致性**: 查询数据在参考空间中，便于比较
- **标签迁移**: 自动继承参考数据的注释
- **批次校正**: 同时校正查询数据的批次效应

## 推荐脚本

### scArches映射
**脚本**: `scarches_mapping_20260127_v1_2_1.py` (v1.2.1 [OK])

**特性**:
- 自动化映射流程
- 标签迁移
- 质量控制
- 映射置信度评估

### 映射质量控制
**脚本**: `query_mapped_to_reference_qc_20260129.py` (v1.0)

**功能**:
- 映射质量评估
- 批次混合检查
- 标签迁移验证
- 异常值检测

### 基础映射
**脚本**: `step2_scarches_mapping_20260112.py` (v1.0)

**功能**:
- 基础scArches映射
- 简化流程

## 输入要求

### 参考模型
```python
# 必需文件
reference_model/
├── model.pt                    # PyTorch模型权重
├── adata.h5ad                  # 参考AnnData (用于setup)
├── var_names.csv               # 基因列表 (关键!)
└── model_config.json           # 模型配置 (可选)
```

### 查询数据
```python
# AnnData对象要求
adata_query.layers['counts']    # 原始counts (整数, 必需!)
adata_query.obs['batch']        # 批次标识 (可选)
adata_query.var_names           # 基因名 (必须与参考匹配)
```

### 基因匹配
**关键**: 查询数据必须包含参考模型的基因

```python
# 加载参考基因列表
ref_genes = pd.read_csv("reference_model/var_names.csv", header=None)[0].tolist()

# 检查重叠
overlap = set(adata_query.var_names) & set(ref_genes)
overlap_ratio = len(overlap) / len(ref_genes)

print(f"Gene overlap: {overlap_ratio:.2%}")

if overlap_ratio < 0.8:
    raise ValueError(f"Insufficient gene overlap: {overlap_ratio:.2%}")

# 子集化到参考基因
adata_query = adata_query[:, ref_genes].copy()
```

## 输出结果

### 文件结构
```
/home/h2048/data/py/{YYMMDD}/scarches_mapping/
├── adata_query_mapped.h5ad             # 映射后的查询数据
├── query_model/                        # 查询特异性模型
│   ├── model.pt
│   └── adata.h5ad
├── figures/
│   ├── umap_reference_and_query.pdf    # 参考+查询联合UMAP
│   ├── umap_batch_integration.pdf      # 批次混合检查
│   ├── label_transfer_confidence.pdf   # 标签迁移置信度
│   └── mapping_quality_metrics.pdf     # 映射质量指标
├── statistics/
│   ├── mapping_summary.txt
│   ├── label_transfer_stats.csv
│   └── batch_mixing_scores.csv
└── README.md
```

### 新增注释
```python
# 映射结果
adata_query.obsm['X_scvi_ref']          # 参考空间中的scVI表示
adata_query.obsm['X_scanvi_ref']        # 参考空间中的scANVI表示
adata_query.obsm['X_umap_ref']          # 参考空间中的UMAP

# 标签迁移
adata_query.obs['predicted_celltype']   # 迁移的细胞类型标签
adata_query.obs['prediction_confidence'] # 预测置信度
adata_query.obs['knn_majority_vote']    # KNN多数投票结果

# 质量指标
adata_query.obs['mapping_score']        # 映射质量评分
adata_query.obs['is_outlier']           # 异常值标记
```

## 工作流程

### Step 1: 准备参考模型
```python
import scanpy as sc
import scvi
import scarches as sca

# 加载参考数据和模型
adata_ref = sc.read_h5ad("reference_data.h5ad")
scvi_model_ref = scvi.model.SCVI.load("reference_model/scvi_model", adata=adata_ref)

# 如果有scANVI模型
scanvi_model_ref = scvi.model.SCANVI.load("reference_model/scanvi_model", adata=adata_ref)

# 提取参考基因列表
ref_genes = adata_ref.var_names.tolist()
pd.Series(ref_genes).to_csv("reference_model/var_names.csv", index=False, header=False)
```

### Step 2: 准备查询数据
```python
# 加载查询数据
adata_query = sc.read_h5ad("query_data.h5ad")

# 基因名标准化 (如果需要)
from utils import normalize_gene_names
normalize_gene_names(adata_query)

# 子集化到参考基因
missing_genes = set(ref_genes) - set(adata_query.var_names)
if missing_genes:
    print(f"WARNING: {len(missing_genes)} genes missing in query")
    # 添加缺失基因 (全0)
    for gene in missing_genes:
        adata_query.var[gene] = 0

# 重排序到参考基因顺序
adata_query = adata_query[:, ref_genes].copy()

# 确认counts层存在
if 'counts' not in adata_query.layers:
    if adata_query.raw is not None:
        adata_query.layers['counts'] = adata_query.raw.X.copy()
    else:
        raise ValueError("counts layer not found!")
```

### Step 3: scVI映射 (Architecture Surgery)
```python
# 使用scArches进行架构手术
scvi_model_query = sca.models.SCVI.load_query_data(
    adata_query,
    reference_model=scvi_model_ref,
    freeze_dropout=True,        # 冻结dropout层
    freeze_expression=True,     # 冻结表达解码器
    freeze_batchnorm=True       # 冻结batch normalization
)

# 手术式微调 (surgical fine-tuning)
scvi_model_query.train(
    max_epochs=50,              # 少量epochs (vs 400+ for training)
    plan_kwargs={'weight_decay': 0.0},
    check_val_every_n_epoch=10
)

# 保存查询模型
scvi_model_query.save("query_model/scvi_model")

# 获取查询数据在参考空间的表示
adata_query.obsm['X_scvi_ref'] = scvi_model_query.get_latent_representation()
```

### Step 4: scANVI映射 (如果有参考scANVI模型)
```python
# 从查询scVI模型初始化scANVI
scanvi_model_query = sca.models.SCANVI.load_query_data(
    adata_query,
    reference_model=scanvi_model_ref,
    freeze_dropout=True,
    freeze_expression=True,
    freeze_batchnorm=True
)

# 微调
scanvi_model_query.train(
    max_epochs=50,
    plan_kwargs={'weight_decay': 0.0}
)

# 保存
scanvi_model_query.save("query_model/scanvi_model")

# 获取表示
adata_query.obsm['X_scanvi_ref'] = scanvi_model_query.get_latent_representation()
```

### Step 5: 标签迁移
```python
# 方法1: scANVI预测 (如果有scANVI模型)
adata_query.obs['predicted_celltype'] = scanvi_model_query.predict()

# 方法2: KNN标签迁移 (基于scVI潜在空间)
from sklearn.neighbors import KNeighborsClassifier

# 训练KNN分类器
knn = KNeighborsClassifier(n_neighbors=15)
knn.fit(
    adata_ref.obsm['X_scvi'],
    adata_ref.obs['cell_type']
)

# 预测查询数据
adata_query.obs['knn_predicted_celltype'] = knn.predict(
    adata_query.obsm['X_scvi_ref']
)

# 预测概率 (置信度)
pred_proba = knn.predict_proba(adata_query.obsm['X_scvi_ref'])
adata_query.obs['prediction_confidence'] = pred_proba.max(axis=1)
```

### Step 6: 联合UMAP可视化
```python
# 合并参考和查询数据
adata_ref.obs['dataset_type'] = 'Reference'
adata_query.obs['dataset_type'] = 'Query'

adata_combined = sc.concat([adata_ref, adata_query], label='dataset_type')

# 使用参考空间的scVI表示计算UMAP
sc.pp.neighbors(adata_combined, use_rep='X_scvi_ref', n_neighbors=30)
sc.tl.umap(adata_combined)

# 可视化
sc.pl.umap(
    adata_combined,
    color=['dataset_type', 'cell_type', 'batch'],
    ncols=3,
    save='_reference_and_query.pdf'
)
```

### Step 7: 映射质量评估
```python
# 批次混合评分
from scib.metrics import silhouette_batch

batch_score = silhouette_batch(
    adata_combined,
    batch_key='batch',
    group_key='cell_type',
    embed='X_scvi_ref'
)
print(f"Batch mixing score: {batch_score:.3f}")

# 标签迁移准确性 (如果有真实标签)
if 'true_celltype' in adata_query.obs.columns:
    from sklearn.metrics import accuracy_score, confusion_matrix

    acc = accuracy_score(
        adata_query.obs['true_celltype'],
        adata_query.obs['predicted_celltype']
    )
    print(f"Label transfer accuracy: {acc:.3f}")

    cm = confusion_matrix(
        adata_query.obs['true_celltype'],
        adata_query.obs['predicted_celltype']
    )
    # 可视化混淆矩阵

# 异常值检测
from scipy.spatial.distance import cdist

# 计算查询细胞到参考细胞的最近距离
distances = cdist(
    adata_query.obsm['X_scvi_ref'],
    adata_ref.obsm['X_scvi'],
    metric='euclidean'
)
min_distances = distances.min(axis=1)

# 标记异常值 (距离过大)
threshold = np.percentile(min_distances, 95)
adata_query.obs['is_outlier'] = min_distances > threshold

print(f"Outliers: {adata_query.obs['is_outlier'].sum()} cells "
      f"({adata_query.obs['is_outlier'].sum()/adata_query.n_obs*100:.2f}%)")
```

## 使用示例

### 基本用法
```bash
# 运行scArches映射
python scarches_mapping_20260127_v1_2_1.py \
    --reference-model /path/to/reference_model \
    --query-data query_data.h5ad \
    --output /home/h2048/data/py/20260203/scarches_mapping

# 预期运行时间 (50k query cells):
# - 数据准备: ~5分钟
# - scVI微调: ~10-20分钟
# - scANVI微调: ~10-20分钟
# - 标签迁移: ~2分钟
# - 总计: ~30-50分钟 (vs 2-3小时从头训练)
```

### 映射质量控制
```bash
# 评估映射质量
python query_mapped_to_reference_qc_20260129.py \
    --mapped-data adata_query_mapped.h5ad \
    --reference-data reference_data.h5ad \
    --output qc_report.html
```

## 关键参数

### Architecture Surgery参数
```python
# 冻结参数 (推荐全部冻结)
FREEZE_DROPOUT = True           # 冻结dropout层
FREEZE_EXPRESSION = True        # 冻结表达解码器
FREEZE_BATCHNORM = True         # 冻结batch normalization
FREEZE_DECODER = False          # 通常不冻结 (需要适应新数据)
```

### 微调参数
```python
MAX_EPOCHS = 50                 # 微调轮数 (vs 400+ for training)
LEARNING_RATE = 1e-3            # 学习率 (可以比训练时高)
WEIGHT_DECAY = 0.0              # 权重衰减
BATCH_SIZE = 128                # 批次大小
```

### 标签迁移参数
```python
KNN_NEIGHBORS = 15              # KNN邻居数
CONFIDENCE_THRESHOLD = 0.5      # 置信度阈值
```

### 质量控制阈值
```python
MIN_GENE_OVERLAP = 0.8          # 最小基因重叠率
MAX_OUTLIER_RATIO = 0.1         # 最大异常值比例
MIN_BATCH_MIXING = 0.7          # 最小批次混合评分
```

## 质量检查

### 运行前
- [ ] 参考模型可用且完整
- [ ] 基因重叠率 >80%
- [ ] 查询数据有counts层
- [ ] 查询数据与参考数据相似 (同组织)

### 运行后
- [ ] 微调loss收敛
- [ ] 查询数据在UMAP上与参考数据混合良好
- [ ] 标签迁移置信度 >0.5 (大部分细胞)
- [ ] 异常值比例 <10%
- [ ] 批次混合评分 >0.7

## 常见问题

### Q1: 基因重叠率低 (<80%)
**原因**: 基因名格式不匹配或测序平台不同
**解决方案**:
```python
# 基因名转换
from utils import normalize_gene_names
normalize_gene_names(adata_query)

# 如果仍然不够，考虑:
# 1. 使用基因ID而非symbol
# 2. 添加缺失基因 (全0)
# 3. 重新训练参考模型 (使用共同基因)
```

### Q2: 微调loss不收敛
**原因**: 学习率过高或查询数据与参考差异大
**解决方案**:
```python
# 降低学习率
scvi_model_query.train(
    max_epochs=50,
    plan_kwargs={'lr': 1e-4}  # 从1e-3降到1e-4
)

# 或增加微调轮数
MAX_EPOCHS = 100  # 从50增加到100
```

### Q3: 查询数据在UMAP上与参考分离
**原因**: 批次效应强或数据差异大
**解决方案**:
```python
# 检查批次标识是否正确
print(adata_query.obs['batch'].value_counts())

# 增加微调轮数
MAX_EPOCHS = 100

# 或考虑不冻结某些层
FREEZE_EXPRESSION = False  # 允许表达解码器适应
```

### Q4: 标签迁移置信度低
**分析**:
```python
# 查看低置信度细胞的分布
low_conf = adata_query[adata_query.obs['prediction_confidence'] < 0.5]
print(low_conf.obs['predicted_celltype'].value_counts())

# 可视化
sc.pl.umap(
    adata_query,
    color='prediction_confidence',
    cmap='viridis'
)
```

**可能原因**:
- 新细胞类型 (不在参考中)
- 过渡状态细胞
- 低质量细胞

**解决方案**:
- 标记低置信度细胞为"Unknown"
- 使用更细粒度的参考模型
- 结合CellTypist等其他方法

### Q5: 异常值比例高 (>20%)
**原因**: 查询数据与参考数据差异大
**解决方案**:
```python
# 检查异常值的特征
outliers = adata_query[adata_query.obs['is_outlier']]

# 质量指标
print(outliers.obs[['n_genes', 'n_counts', 'pct_counts_mt']].describe())

# 可能需要:
# 1. 更严格的质量控制
# 2. 使用更相似的参考数据
# 3. 从头训练而非迁移学习
```

## 高级用法

### 多参考模型集成
```python
# 使用多个参考模型进行映射
references = [
    'reference_model_1',
    'reference_model_2',
    'reference_model_3'
]

predictions = []
for ref_path in references:
    model_ref = scvi.model.SCVI.load(ref_path, adata=adata_ref)
    model_query = sca.models.SCVI.load_query_data(
        adata_query,
        reference_model=model_ref
    )
    model_query.train(max_epochs=50)

    pred = model_query.predict()
    predictions.append(pred)

# 多数投票
from scipy.stats import mode
adata_query.obs['ensemble_prediction'] = mode(predictions, axis=0)[0]
```

### 增量学习
```python
# 将查询数据添加到参考中，用于未来映射
adata_ref_updated = sc.concat([adata_ref, adata_query])

# 重新训练 (或微调)
scvi_model_updated = scvi.model.SCVI(adata_ref_updated, ...)
scvi_model_updated.train(max_epochs=200)

# 保存更新的参考
scvi_model_updated.save("reference_model_updated")
```

### 跨物种映射
```python
# 使用同源基因进行跨物种映射
import mygene

mg = mygene.MyGeneInfo()

# 人类基因 → 小鼠同源基因
human_genes = adata_ref.var_names.tolist()
query = mg.querymany(
    human_genes,
    scopes='symbol',
    fields='homologene.genes',
    species='human'
)

# 构建人类-小鼠基因映射
human_to_mouse = {}
for result in query:
    if 'homologene' in result:
        for gene in result['homologene']['genes']:
            if gene['taxid'] == 10090:  # 小鼠
                human_to_mouse[result['query']] = gene['symbol']

# 映射小鼠数据到人类参考
# ...
```

## 最佳实践

1. [OK] **确保基因重叠 >80%** - 否则考虑重新训练
2. [OK] **冻结大部分层** - 只微调必要的部分
3. [OK] **使用少量epochs** - 避免过拟合
4. [OK] **验证批次混合** - 确保查询数据与参考整合良好
5. [OK] **检查异常值** - 标记并分析异常细胞
6. [OK] **多方法验证** - 结合KNN、scANVI等多种标签迁移方法
7. [OK] **保存查询模型** - 便于未来使用
8. [OK] **记录参数** - 便于重现

## 相关Skills

- `allcells-integration` - 全细胞整合 (替代方案)
- `celltype-specific-analysis` - 细胞类型分析 (下游)
- `doublet-detection` - 质量控制 (上游)

## 参考文献

- scArches: Lotfollahi et al., Nature Biotechnology 2022
- scVI: Lopez et al., Nature Methods 2018
- scANVI: Xu et al., Molecular Systems Biology 2021
- Transfer learning review: Zhuang et al., Proceedings of the IEEE 2021
