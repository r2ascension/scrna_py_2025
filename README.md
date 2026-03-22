# 🫁 scRNA-seq 肺组织单细胞分析流程 (scrna_py_2025)

> 一套面向人类肺组织的**生产级** GPU 加速单细胞 RNA 测序 (scRNA-seq) 分析管道，集批次校正、自动细胞类型注释、亚群精细分析与迁移学习于一体。

---

## 📖 项目简介

本仓库包含 43 个 Python 脚本与 11 个 Jupyter Notebook（共 54 个文件），构成一套完整的 scRNA-seq 分析工作流，主要用于：

- 多样本/批次人类肺组织单细胞转录组数据整合
- 基于 scVI/scANVI 的深度学习批次校正
- 使用 CellTypist 进行自动化细胞类型注释（Human Lung Atlas 模型）
- 各细胞谱系的细粒度亚群鉴定
- 基于 scArches 的迁移学习，将新数据映射到已有参考图谱

**数据来源**: 人类肺组织 scRNA-seq（多批次多样本，h5ad 格式）  
**主要细胞谱系**: 上皮细胞、T/NK 细胞、B 细胞、髓系细胞、基质/血管细胞

---

## ✨ 核心特性

| 特性 | 说明 |
|------|------|
| **批次校正** | scVI (n_latent=100, n_layers=2)、BBKNN、Harmony |
| **自动注释** | CellTypist 多数投票（Human_Lung_Atlas + Immune_All_Low） |
| **半监督精炼** | scANVI 标签引导细化 |
| **亚群聚类** | Leiden + restrict_to 保留 BBKNN 图 |
| **基因程序分析** | cNMF（K 值选择 + GEP 映射） |
| **双细胞检测** | Top2 谱系 + margin 逻辑 |
| **迁移学习** | scArches surgical fine-tuning |
| **内存优化** | HVG 训练减少 60–70% 内存占用 |
| **GPU 加速** | 单卡 CUDA 训练，CPU 安全回退 |
| **完全可重复** | 固定随机种子 + HVG 列表保存 |

---

## 🛠️ 依赖环境

```bash
# 核心单细胞分析
scvi-tools>=1.1.4
scanpy
anndata
celltypist

# 批次校正
bbknn

# 深度学习
torch
lightning>=2.2,<3

# 数据处理 & 可视化
numpy
scipy
pandas
matplotlib
seaborn
umap-learn
h5py
mygene
statsmodels
scikit-learn
```

> **推荐安装方式**（conda 环境）：
> ```bash
> conda create -n scrna python=3.10
> conda activate scrna
> pip install scvi-tools celltypist bbknn scanpy anndata lightning mygene
> ```

---

## 📂 脚本目录

> 详细分类与版本说明请参阅 **[SCRIPT_CATEGORIZATION.md](./SCRIPT_CATEGORIZATION.md)**

### 1️⃣ 全细胞整合分析
| 脚本 | 版本 | 状态 |
|------|------|------|
| `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py` | v2.5.2 | ✅ 生产推荐 |

### 2️⃣ 细胞类型特异性分析
| 细胞谱系 | 推荐脚本 | 版本 |
|---------|---------|------|
| 上皮细胞 | `complete_epithelial_analysis_pipeline_20260122_v2_1.py` | v2.1 |
| T/NK 细胞 | `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` | v3.5.1 |
| B 细胞 | `b_scvi_celltypist_scanvi_pipeline_v1.py` | v1.0 |
| 髓系细胞 | `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` | v3.5.1 |
| 基质/血管 | `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py` | v3.5.1 |

### 3️⃣ 亚群聚类 & 通用流程
| 脚本 | 版本 | 说明 |
|------|------|------|
| `universal_celltype_subcluster_pipeline_20260121_v4_1.py` | v4.1 | ✅ 推荐通用亚群流程 |

### 4️⃣ 质量控制
| 脚本 | 说明 |
|------|------|
| `remove_doublets_and_gse299751_20260114.py` | 双细胞检测与移除 |
| `quick_data_check_detailed_20260127.py` | h5ad 结构检查 |

### 5️⃣ 迁移学习 (scArches)
| 脚本 | 版本 | 说明 |
|------|------|------|
| `scarches_mapping_20260127_v1_2_1.py` | v1.2.1 | ✅ 推荐参考模型映射 |

---

## 🚀 快速开始

### 场景 1：新数据全流程分析

```bash
# Step 1: 全细胞整合 + 自动注释
python allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py

# Step 2: 双细胞检测
python remove_doublets_and_gse299751_20260114.py

# Step 3: 各谱系细化分析（按需选择）
python t_scvi_celltypist_scanvi_20260108_v3_5_1.py        # T/NK 细胞
python b_scvi_celltypist_scanvi_pipeline_v1.py             # B 细胞
python myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py   # 髓系
python stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py   # 基质

# Step 4: 亚群精细分析
python universal_celltype_subcluster_pipeline_20260121_v4_1.py
```

### 场景 2：上皮细胞一站式分析（BBKNN + Harmony + cNMF）

```bash
python complete_epithelial_analysis_pipeline_20260122_v2_1.py
```

### 场景 3：将新数据映射到已有参考图谱

```bash
python scarches_mapping_20260127_v1_2_1.py
```

---

## 📥 数据格式

**输入**: h5ad（AnnData）格式，需包含：
- `layers['counts']` — 原始基因表达计数
- `layers['log1p']` — Log 标准化表达
- `.raw` — 全基因备份
- `obs` 元数据列 — `batch`、`sample`、`cell_type` 等

**输出**: 带注释、嵌入与亚群标签的 h5ad 文件

---

## 📈 版本演进

| 系列 | 关键改进 |
|------|---------|
| v2.x (All-cells) | 基因对齐、UMAP 分离、CPU 安全 CUDA、索引对齐写入 |
| v2.3.x | 协变量系统 (MT%, stress, cell cycle)、Unknown 清洗、类别不平衡处理 |
| v3.5.x | P0 修复 (scANVI 标签对齐)、dtype 安全、内存优化 60–70% |
| v4.x | `restrict_to` 保留 BBKNN 图、技术基因过滤、批次效应监控 |

---

## 📝 文件命名规范

```
{celltype}_{analysis}_{YYYYMMDD}_v{X}_{Y}.py

示例：
  allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py
  stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py
```

---

## 📊 仓库统计

| 指标 | 数量 |
|------|------|
| Python 脚本 (.py) | 43 |
| Jupyter Notebook (.ipynb) | 11 |
| 总代码行数 | ~56,000 |
| 覆盖细胞谱系 | 5 (上皮/T-NK/B/髓系/基质) |

---

## 📄 文档

- **[SCRIPT_CATEGORIZATION.md](./SCRIPT_CATEGORIZATION.md)** — 全部脚本详细分类、推荐版本与快速查找指南

---

## 🔗 相关资源

- [scVI-tools 文档](https://scvi-tools.org/)
- [CellTypist 文档](https://celltypist.readthedocs.io/)
- [Scanpy 文档](https://scanpy.readthedocs.io/)
- [scArches 文档](https://scarches.readthedocs.io/)
- [Human Lung Cell Atlas](https://www.humancellatlas.org/)
