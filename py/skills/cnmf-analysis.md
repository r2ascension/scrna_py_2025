---
name: cnmf-analysis
description: Identify gene expression programs using consensus NMF to reveal cell states, functional modules, and transcriptional patterns.
license: Complete terms in LICENSE.txt
metadata:
  title: cNMF Gene Expression Program Analysis
  version: "1.1"
  author: Single-cell Analysis Pipeline Team
  category: analysis
  difficulty: advanced
  estimated_time: "2-6 hours"
  tags:
    - cnmf
    - gene-expression-programs
    - nmf
    - transcriptional-modules
    - cell-states
    - functional-analysis
  prerequisites:
    - "Raw integer counts in layers['counts']"
    - Cell type labels for label-guided analysis
    - Minimum 1000 cells (5k+ recommended)
    - Technical genes filtered (MT, Ribo, Histone)
    - Multi-processing capability for speedup
  outputs:
    - "GEP usage matrix (obsm['X_cnmf_usages'])"
    - GEP gene weights and top genes
    - K selection plots for optimal GEP number
    - GEP-cluster association heatmaps
    - Pathway enrichment results per GEP
  related_skills:
    - celltype-specific-analysis
    - subcluster-analysis
    - allcells-integration
---
# Skill: cNMF Gene Expression Program Analysis

## 功能描述

使用共识非负矩阵分解(cNMF)识别基因表达程序(GEPs)，揭示细胞状态、功能模块和转录调控模式。

## 核心技术

- **cNMF**: Consensus Non-negative Matrix Factorization
- **基因过滤**: 去除技术伪影基因
- **K值选择**: 稳定性分析确定最优GEP数量
- **GEP映射**: 将GEP关联到细胞聚类和状态

## 适用场景

- [OK] 识别细胞状态和功能模块
- [OK] 发现转录调控程序
- [OK] 比较不同批次校正方法 (BBKNN vs Harmony)
- [OK] 疾病相关基因表达模式
- [OK] 发育轨迹和细胞分化
- [NO] 不适用于极小数据集 (<1000 cells)

## 核心概念

### 什么是GEP (Gene Expression Program)?
GEP是一组共表达的基因集合，代表特定的生物学功能或细胞状态:

- **细胞周期GEP**: MKI67, TOP2A, PCNA等
- **应激响应GEP**: HSPA1A, HSPA1B, FOS, JUN等
- **免疫激活GEP**: IFNG, TNF, IL2等
- **代谢GEP**: 糖酵解、氧化磷酸化相关基因

### cNMF原理
```
表达矩阵 (Cells × Genes)
    ↓ 分解
Usage矩阵 (Cells × K) × GEP矩阵 (K × Genes)

其中:
- K = GEP数量 (需要选择)
- Usage = 每个细胞对每个GEP的使用程度
- GEP = 每个GEP中基因的权重
```

## 推荐脚本

### 标签引导cNMF (推荐)
**脚本**: `label_guided_cnmf_pipeline_20260114_v1_1.py` (v1.1 [OK])

**特性**:
- 按细胞类型分别运行cNMF
- 置信度过滤 (>0.5)
- 自动化K值选择
- GEP-聚类映射

### 结果分析
**脚本**: `cnmf_results_analysis_20260114_v1_1.py` (v1.1)

**功能**:
- GEP可视化
- Top基因提取
- 通路富集分析
- GEP-聚类相关性

### 完整流程 (BBKNN + Harmony + cNMF)
**脚本**: `complete_epithelial_analysis_pipeline_20260122_v2_1.py` (v2.1)

**特性**:
- BBKNN批次校正
- Harmony批次校正
- 在两种校正数据上运行cNMF
- 比较批次校正方法对GEP的影响

### 快速启动脚本 (按细胞类型)

| 细胞类型 | 脚本 | K范围 |
|---------|------|-------|
| B cells | `quick_start_cnmf_analysis_20251225_B_v1_1.py` | 15-35 |
| T cells | `quick_start_cnmf_analysis_20251225_T_v1_1.py` | 20-40 |
| Myeloid | `quick_start_cnmf_analysis_20251225_Myeloid_v1_1.py` | 20-40 |
| Stromal | `quick_start_cnmf_analysis_20251227_Stromal_v1_1.py` | 15-35 |
| Ciliated | `quick_start_cnmf_analysis_ciliated_20251227_v1_1.py` | 10-25 |
| Club | `quick_start_cnmf_analysis_club_20251227_v1_1.py` | 10-25 |
| AT2 | `quick_start_cnmf_analysis_AT2_20251227_v1_1.py` | 10-25 |
| Goblet | `quick_start_cnmf_analysis_goblet_20260105_v1_1.py` | 10-25 |

## 输入要求

### 数据格式
```python
# AnnData对象
adata.layers['counts']         # 原始整数counts (必需!)
adata.obs['cell_type']         # 细胞类型标签
adata.obs['leiden']            # 聚类标签 (用于GEP映射)
adata.obs['scanvi_confidence'] # 置信度 (可选，用于过滤)
```

### 数据规模建议
- **最小**: 1000 cells
- **推荐**: 5k-50k cells per cell type
- **大规模**: >50k cells (使用多进程加速)

### 基因过滤 (关键!)
**必须过滤技术伪影基因**，否则GEP会被技术变异主导:

```python
FILTER_GENES = {
    'remove_mt': True,           # MT-* 线粒体基因
    'remove_ribo': True,         # RPS/RPL/MRPS/MRPL 核糖体
    'remove_histone': True,      # H1/H2A/H2B/H3/H4 组蛋白
    'remove_pseudogenes': True,  # *P* 假基因
    'remove_ensg': True,         # ENSG* ENSEMBL IDs
    'remove_unannotated': True   # AC/AL/RP/CTD/LINC 未注释
}

# 过滤函数
def filter_technical_genes(adata):
    keep_genes = (
        ~adata.var_names.str.match(r'^MT-') &
        ~adata.var_names.str.match(r'^RP[SL]') &
        ~adata.var_names.str.match(r'^MRP[SL]') &
        ~adata.var_names.str.match(r'^H[1234]') &
        ~adata.var_names.str.match(r'^HIST') &
        ~adata.var_names.str.match(r'ENSG\d+') &
        ~adata.var_names.str.match(r'^(LINC|AC\d+|AL\d+|RP11-|CTD-|CTB-)')
    )
    return adata[:, keep_genes].copy()
```

## 输出结果

### 文件结构
```
/home/h2048/data/py/{YYMMDD}/cnmf_{celltype}/
├── cnmf_output/
│   ├── {celltype}_cnmf.spectra.k_{K}.dt_0_01.consensus.txt    # GEP矩阵
│   ├── {celltype}_cnmf.usages.k_{K}.dt_0_01.consensus.txt     # Usage矩阵
│   ├── {celltype}_cnmf.gene_spectra_score.k_{K}.dt_0_01.txt   # 基因权重
│   ├── {celltype}_cnmf.gene_spectra_tpm.k_{K}.dt_0_01.txt     # TPM标准化
│   └── {celltype}_cnmf.k_selection.png                         # K选择图
├── figures/
│   ├── gep_umap_overlay.pdf                # GEP在UMAP上的分布
│   ├── gep_cluster_heatmap.pdf             # GEP-聚类热图
│   ├── top_genes_per_gep.pdf               # 每个GEP的top基因
│   └── gep_correlation_matrix.pdf          # GEP间相关性
├── results/
│   ├── gep_top_genes.csv                   # 每个GEP的top基因列表
│   ├── gep_cluster_association.csv         # GEP与聚类的关联
│   └── gep_enrichment/                     # 通路富集结果
│       ├── GEP_1_enrichment.csv
│       ├── GEP_2_enrichment.csv
│       └── ...
└── README.md
```

### AnnData新增内容
```python
# GEP usage scores
adata.obsm['X_cnmf_usages']        # (n_cells, K) usage矩阵

# GEP信息
adata.uns['cnmf_geps']             # GEP基因权重
adata.uns['cnmf_top_genes']        # 每个GEP的top基因

# GEP-聚类关联
adata.obs['dominant_gep']          # 每个细胞的主导GEP
adata.obs['gep_entropy']           # GEP使用的熵 (多样性)
```

## 工作流程

### Step 1: 数据准备和基因过滤
```python
import scanpy as sc
from cnmf import cNMF

# 加载数据
adata = sc.read_h5ad("adata_epithelial_FINAL.h5ad")

# 按细胞类型子集化
celltype = "Basal"
adata_subset = adata[adata.obs['cell_type'] == celltype].copy()

# 置信度过滤 (可选但推荐)
if 'scanvi_confidence' in adata_subset.obs.columns:
    adata_subset = adata_subset[
        adata_subset.obs['scanvi_confidence'] > 0.5
    ].copy()

# 基因过滤 (关键!)
adata_filtered = filter_technical_genes(adata_subset)

print(f"Cells: {adata_filtered.n_obs}")
print(f"Genes after filtering: {adata_filtered.n_vars}")
```

### Step 2: cNMF初始化和准备
```python
# 设置输出目录
output_dir = f"/home/h2048/data/py/20260203/cnmf_{celltype}"
os.makedirs(output_dir, exist_ok=True)

# 初始化cNMF对象
cnmf_obj = cNMF(
    output_dir=output_dir,
    name=f"{celltype}_cnmf"
)

# K值范围 (根据细胞数量调整)
n_cells = adata_filtered.n_obs
if n_cells > 10000:
    k_range = np.arange(25, 55, 5)      # [25, 30, 35, 40, 45, 50]
elif n_cells > 5000:
    k_range = np.arange(15, 40, 5)      # [15, 20, 25, 30, 35]
else:
    k_range = np.arange(10, 30, 5)      # [10, 15, 20, 25]

# 准备数据
cnmf_obj.prepare(
    counts_fn=adata_filtered,           # 可以直接传AnnData
    components=k_range,                 # K值范围
    n_iter=100,                         # 迭代次数 (100推荐)
    seed=14,                            # 随机种子
    num_hvg=2000                        # HVG数量 (过滤后)
)
```

### Step 3: 矩阵分解 (可多进程)
```python
# 单进程运行
cnmf_obj.factorize(worker_i=0, total_workers=1)

# 或多进程加速 (推荐)
total_workers = 16
for worker_i in range(total_workers):
    cnmf_obj.factorize(
        worker_i=worker_i,
        total_workers=total_workers
    )

# 预期加速: 4-16x (取决于worker数量)
```

### Step 4: 合并结果和K值选择
```python
# 合并所有worker的结果
cnmf_obj.combine()

# 生成K选择图
cnmf_obj.k_selection_plot()

# 查看K选择图，选择最优K
# 标准:
# 1. 稳定性高 (stability score高)
# 2. 误差低 (reconstruction error低)
# 3. 生物学可解释性

# 手动选择K (示例)
selected_k = 25
```

### Step 5: 加载结果
```python
# 加载选定K的结果
usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(
    K=selected_k,
    density_threshold=0.01  # 稀疏性阈值
)

# usage_norm: (n_cells, K) - 每个细胞的GEP使用程度
# gep_scores: (K, n_genes) - 每个GEP的基因权重
# gep_tpm: (K, n_genes) - TPM标准化的基因权重
# topgenes: dict - 每个GEP的top基因
```

### Step 6: 将结果添加到AnnData
```python
# 添加usage矩阵
adata_subset.obsm['X_cnmf_usages'] = usage_norm

# 添加主导GEP
adata_subset.obs['dominant_gep'] = usage_norm.argmax(axis=1).astype(str)

# 计算GEP熵 (多样性)
from scipy.stats import entropy
adata_subset.obs['gep_entropy'] = entropy(usage_norm, axis=1)

# 保存GEP信息
adata_subset.uns['cnmf_geps'] = gep_scores
adata_subset.uns['cnmf_top_genes'] = topgenes
```

### Step 7: GEP可视化
```python
# 在UMAP上可视化GEP usage
for k in range(selected_k):
    adata_subset.obs[f'GEP_{k+1}'] = usage_norm[:, k]

sc.pl.umap(
    adata_subset,
    color=[f'GEP_{k+1}' for k in range(selected_k)],
    cmap='viridis',
    ncols=5,
    save='_gep_umap_overlay.pdf'
)

# 可视化主导GEP
sc.pl.umap(
    adata_subset,
    color='dominant_gep',
    legend_loc='on data',
    save='_dominant_gep.pdf'
)
```

### Step 8: GEP-聚类关联分析
```python
# 计算每个聚类的平均GEP usage
cluster_gep_mean = adata_subset.obs.groupby('leiden')[
    [f'GEP_{k+1}' for k in range(selected_k)]
].mean()

# 热图可视化
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 8))
sns.heatmap(
    cluster_gep_mean.T,
    cmap='RdBu_r',
    center=0,
    cbar_kws={'label': 'Mean GEP usage'}
)
plt.xlabel('Cluster')
plt.ylabel('GEP')
plt.title(f'{celltype} - GEP-Cluster Association')
plt.tight_layout()
plt.savefig(f'{output_dir}/figures/gep_cluster_heatmap.pdf')
```

### Step 9: Top基因提取和注释
```python
# 提取每个GEP的top基因
n_top_genes = 50
gep_top_genes_df = []

for k in range(selected_k):
    top_gene_indices = gep_scores[k, :].argsort()[-n_top_genes:][::-1]
    top_genes = adata_filtered.var_names[top_gene_indices]
    top_scores = gep_scores[k, top_gene_indices]

    for gene, score in zip(top_genes, top_scores):
        gep_top_genes_df.append({
            'GEP': f'GEP_{k+1}',
            'Gene': gene,
            'Score': score
        })

gep_top_genes_df = pd.DataFrame(gep_top_genes_df)
gep_top_genes_df.to_csv(f'{output_dir}/results/gep_top_genes.csv', index=False)
```

### Step 10: 通路富集分析 (可选)
```python
import gseapy as gp

for k in range(selected_k):
    # 获取GEP的top基因
    top_genes = gep_top_genes_df[
        gep_top_genes_df['GEP'] == f'GEP_{k+1}'
    ]['Gene'].head(100).tolist()

    # 运行富集分析
    enr = gp.enrichr(
        gene_list=top_genes,
        gene_sets=[
            'GO_Biological_Process_2021',
            'KEGG_2021_Human',
            'Reactome_2022'
        ],
        organism='human',
        outdir=f'{output_dir}/results/gep_enrichment/GEP_{k+1}'
    )

    # 保存结果
    enr.results.to_csv(
        f'{output_dir}/results/gep_enrichment/GEP_{k+1}_enrichment.csv'
    )
```

## 使用示例

### 基本用法 (标签引导)
```bash
# 运行标签引导cNMF
python label_guided_cnmf_pipeline_20260114_v1_1.py \
    --input adata_epithelial_FINAL.h5ad \
    --celltype Basal \
    --k-range 15,20,25,30,35 \
    --output /home/h2048/data/py/20260203/cnmf_basal

# 预期运行时间 (10k cells, 16 workers):
# - Prepare: ~5分钟
# - Factorize: ~30-60分钟
# - Combine: ~5分钟
# - 总计: ~1-2小时
```

### 快速启动 (T细胞)
```bash
python quick_start_cnmf_analysis_20251225_T_v1_1.py

# 自动化流程:
# 1. 加载T细胞数据
# 2. 基因过滤
# 3. cNMF运行 (K=20-40)
# 4. K选择
# 5. 结果可视化
```

### 完整流程 (BBKNN + Harmony比较)
```bash
# 运行完整上皮细胞分析 (包含cNMF)
python complete_epithelial_analysis_pipeline_20260122_v2_1.py

# 流程:
# 1. BBKNN批次校正
# 2. Harmony批次校正
# 3. 在BBKNN数据上运行cNMF
# 4. 在Harmony数据上运行cNMF
# 5. 比较两种方法的GEP差异
```

## 关键参数

### cNMF核心参数
```python
# K值范围 (根据细胞数量)
K_RANGE_SMALL = [10, 15, 20, 25]           # <5k cells
K_RANGE_MEDIUM = [15, 20, 25, 30, 35]      # 5k-10k cells
K_RANGE_LARGE = [25, 30, 35, 40, 45, 50]   # >10k cells

# 迭代次数
N_ITER = 100                                # 推荐值 (平衡速度和稳定性)

# HVG数量
NUM_HVG = 2000                              # 过滤技术基因后

# 稀疏性阈值
DENSITY_THRESHOLD = 0.01                    # 默认值
```

### 多进程参数
```python
TOTAL_WORKERS = 16                          # 推荐4-16
# 加速比: ~0.7 * TOTAL_WORKERS (有overhead)
```

### K值选择标准
```python
# 查看K选择图，综合考虑:
# 1. Stability score (越高越好)
# 2. Reconstruction error (越低越好)
# 3. 生物学可解释性 (top基因是否有意义)

# 经验法则:
# - 小细胞类型 (<5k): K=10-20
# - 中等细胞类型 (5k-20k): K=15-30
# - 大细胞类型 (>20k): K=25-50
```

## 质量检查

### 运行前
- [ ] 确认 `layers['counts']` 为整数
- [ ] 技术基因已过滤
- [ ] 细胞数 >1000
- [ ] 基因数 >1000 (过滤后)

### 运行后
- [ ] K选择图显示明确的最优K
- [ ] Top基因生物学意义明确
- [ ] GEP之间相关性不过高 (<0.8)
- [ ] GEP-聚类关联合理
- [ ] 无技术伪影GEP (MT/Ribo主导)

## 常见问题

### Q1: cNMF运行失败 "counts not integers"
**原因**: counts层不是整数
**解决方案**:
```python
# 确保counts为整数
adata.layers['counts'] = adata.layers['counts'].astype(int)

# 或从raw恢复
if adata.raw is not None:
    adata.layers['counts'] = adata.raw.X.astype(int)
```

### Q2: GEP主要是MT/Ribo基因
**原因**: 未过滤技术基因
**解决方案**: 应用技术基因过滤 (见输入要求部分)

### Q3: K选择图没有明确最优值
**分析**:
```python
# 查看多个K的top基因
for k in [15, 20, 25, 30]:
    _, _, _, topgenes = cnmf_obj.load_results(K=k)
    print(f"\nK={k}:")
    for gep_id, genes in topgenes.items():
        print(f"  GEP {gep_id}: {', '.join(genes[:10])}")
```

**选择策略**:
- 选择生物学可解释性最好的K
- 偏向较小的K (避免过拟合)
- 可以运行多个K并比较

### Q4: GEP之间高度相关
**原因**: K值过大或数据异质性不足
**解决方案**:
```python
# 计算GEP相关性
from scipy.stats import spearmanr
gep_corr = spearmanr(usage_norm)[0]

# 如果相关性>0.8，减小K值
if (gep_corr > 0.8).sum() > selected_k:
    print("GEPs highly correlated, consider reducing K")
```

### Q5: 多进程运行内存不足
**解决方案**:
```python
# 减少worker数量
TOTAL_WORKERS = 4  # 从16减到4

# 或减少HVG数量
NUM_HVG = 1500  # 从2000减到1500
```

## 高级用法

### 比较不同批次校正方法
```python
# 在BBKNN数据上运行cNMF
adata_bbknn = sc.read_h5ad("adata_bbknn.h5ad")
cnmf_bbknn = run_cnmf(adata_bbknn, name="bbknn")

# 在Harmony数据上运行cNMF
adata_harmony = sc.read_h5ad("adata_harmony.h5ad")
cnmf_harmony = run_cnmf(adata_harmony, name="harmony")

# 比较GEP相似性
from scipy.spatial.distance import cosine
for k in range(selected_k):
    similarity = 1 - cosine(
        cnmf_bbknn.gep_scores[k, :],
        cnmf_harmony.gep_scores[k, :]
    )
    print(f"GEP {k+1} similarity: {similarity:.3f}")
```

### GEP轨迹分析
```python
# 使用GEP usage进行轨迹推断
import scvelo as scv

# 将GEP usage作为"表达"数据
adata_gep = sc.AnnData(
    X=usage_norm,
    obs=adata_subset.obs,
    var=pd.DataFrame(index=[f'GEP_{k+1}' for k in range(selected_k)])
)

# 运行RNA velocity (如果有spliced/unspliced数据)
scv.tl.velocity(adata_gep)
scv.tl.velocity_graph(adata_gep)
scv.pl.velocity_embedding_stream(adata_gep, basis='umap')
```

### GEP调控网络推断
```python
# 使用SCENIC推断GEP的转录因子调控
import pyscenic

# 对每个GEP的top基因运行SCENIC
for k in range(selected_k):
    top_genes = topgenes[k][:100]

    # 运行SCENIC (需要预先准备的数据库)
    regulons = pyscenic.run_scenic(
        adata_subset[:, top_genes],
        grn_output='regulons.csv',
        ctx_output='ctx.csv'
    )
```

## 解释GEP的策略

### 1. Top基因分析
```python
# 查看每个GEP的top 20基因
for k in range(selected_k):
    print(f"\nGEP {k+1}:")
    print(", ".join(topgenes[k][:20]))

# 手动注释GEP功能
gep_annotations = {
    'GEP_1': 'Cell cycle / Proliferation',
    'GEP_2': 'Interferon response',
    'GEP_3': 'Stress response',
    # ...
}
```

### 2. 通路富集
使用GO/KEGG/Reactome富集分析识别GEP功能

### 3. 与已知标记基因比较
```python
# 定义已知功能模块
KNOWN_MODULES = {
    'cell_cycle': ['MKI67', 'TOP2A', 'PCNA', 'CDK1'],
    'interferon': ['ISG15', 'MX1', 'IFIT1', 'IFIT3'],
    'stress': ['HSPA1A', 'HSPA1B', 'FOS', 'JUN']
}

# 计算GEP与已知模块的重叠
for gep_id, genes in topgenes.items():
    for module_name, module_genes in KNOWN_MODULES.items():
        overlap = len(set(genes[:50]) & set(module_genes))
        if overlap > 2:
            print(f"GEP {gep_id} overlaps with {module_name}: {overlap} genes")
```

### 4. GEP-聚类关联
查看哪些聚类高表达特定GEP，结合聚类的生物学特征推断GEP功能

## 最佳实践

1. [OK] **必须过滤技术基因** - MT/Ribo/Histone/IEG
2. [OK] **按细胞类型分别运行** - 避免细胞类型差异主导GEP
3. [OK] **使用置信度过滤** - 只保留高质量细胞
4. [OK] **测试多个K值** - 不要只运行单个K
5. [OK] **验证生物学意义** - 通路富集 + 已知标记基因
6. [OK] **使用多进程** - 4-16 workers加速
7. [OK] **保存中间结果** - 便于重新分析
8. [OK] **比较批次校正方法** - BBKNN vs Harmony

## 相关Skills

- `celltype-specific-analysis` - 细胞类型分析 (上游)
- `subcluster-analysis` - 亚群聚类 (互补)
- `allcells-integration` - 全细胞整合 (上游)

## 参考文献

- cNMF: Kotliar et al., eLife 2019
- NMF: Lee & Seung, Nature 1999
- Gene expression programs: Jerby-Arnon et al., Cell 2018
