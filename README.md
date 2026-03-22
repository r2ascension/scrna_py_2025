# 🫁 scRNA-seq 肺组织单细胞分析流程 (scrna_py_2025)

> 一套面向人类肺组织的**生产级** GPU 加速单细胞 RNA 测序 (scRNA-seq) 分析管道，集批次校正、自动细胞类型注释、亚群精细分析与迁移学习于一体。

---

## 📖 项目简介

本仓库包含 **42 个 Python 脚本**与 **11 个 Jupyter Notebook**（共 53 个文件），构成一套完整的 scRNA-seq 分析工作流，主要用于：

- 多样本/批次人类肺组织单细胞转录组数据整合
- 基于 scVI/scANVI 的深度学习批次校正
- 使用 CellTypist 进行自动化细胞类型注释（Human Lung Atlas 模型）
- 各细胞谱系的细粒度亚群鉴定
- 基于 scArches 的迁移学习，将新数据映射到已有参考图谱

**数据来源**: 人类肺组织 scRNA-seq（多批次多样本，h5ad 格式）  
**主要细胞谱系**: 上皮/基底细胞、T/NK 细胞、B 细胞、髓系细胞、基质/血管细胞

---

## ✨ 核心特性

| 特性 | 说明 |
|------|------|
| **批次校正** | scVI (n_latent=100, n_layers=2)、BBKNN、Harmony |
| **自动注释** | CellTypist 多数投票（Human_Lung_Atlas + Immune_All_Low） |
| **半监督精炼** | scANVI 标签引导细化 |
| **亚群聚类** | Leiden + restrict_to 保留 BBKNN 图 |
| **双细胞检测** | Top2 谱系 + margin 逻辑 |
| **迁移学习** | scArches surgical fine-tuning |
| **内存优化** | HVG 训练减少 60–70% 内存占用 |
| **GPU 加速** | 单卡 CUDA 训练，CPU 安全回退 |
| **完全可重复** | 固定随机种子 + HVG 列表保存 + scArches dry-run 验证 |

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

## 📂 完整脚本目录

> 推荐版本以 ✅ 标注；详细分类说明另见 **[SCRIPT_CATEGORIZATION.md](./SCRIPT_CATEGORIZATION.md)**

### 1️⃣ 全细胞整合分析

| 脚本 | 版本 | 说明 |
|------|------|------|
| `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py` | v2.5.2 | ✅ 生产推荐：scVI → CellTypist → scANVI 完整流程 |
| `allcells_integrate_scanvi_annotations_pipeline_20251216_v1_1.py` | v1.1 | 将各谱系注释回写到全细胞数据集 |
| `allcells_integrate_scanvi_annotations_pipeline_v1.py` | v1.0 | 注释整合初始版本 |
| `allcells_combined_scanvi_20260319_v1_5.ipynb` | v1.5 | 全细胞联合分析与可视化（Notebook） |

### 2️⃣ 细胞类型特异性分析

#### T/NK 细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` | v3.5.1 | ✅ 生产推荐：scVI → CellTypist → scANVI（P0 修复） |
| `tcell_only_merged_pipeline_20260318_v2_4.py` | v2.4 | ✅ 最新：T/NK 合并管道（CellTypist 过滤后重训练） |
| `tnk_scvi_scanvi_scarches_ref_20260315_v1_2.py` | v1.2 | T/NK scArches 参考模型训练（encode_covariates=True） |
| `tcell_scvi_celltypist_scanvi_pipeline_20251206_v2.py` | v2.0 | T 细胞早期版本管道 |
| `tcell_scvi_pipeline_v1.4.py` | v1.4 | T 细胞 scVI 流程（附 UMAP 管理） |
| `tcell_scvi_integration.py` | — | T 细胞基础 scVI 训练脚本 |
| `tcell_subcluster_analysis_20260126_v1_0.py` | v1.0 | T 细胞亚群分析（CD4/CD8 精细分型） |
| `tcell_bbknn_marker_analysis_20260106.py` | v1.0 | T 细胞 BBKNN 批次校正 + 标记基因分析 |
| `tcell_annotation_visualization_20260214.py` | v1.1 | T 细胞注释可视化与验证 |
| `tcell_myeloid_ref_query_merge_pipeline_20260225_v1_0.py` | v1.0 | T 细胞 + 髓系联合参考/查询合并管道 |
| `tcell_myeloid_subcluster_scanvi_based_20260225_v1_0.py` | v1.0 | T 细胞 + 髓系联合亚群分析（scANVI 基础） |
| `tcell_scvi_pipeline_v1.4.ipynb` | v1.4 | T 细胞 scVI 流程（Notebook） |
| `tcell_bbknn_marker_analysis_filtered_20260107.ipynb` | — | T 细胞 BBKNN 标记分析（过滤版，Notebook） |
| `tcell_cd4cd8_scanvi_20260109_v1_3_1.ipynb` | v1.3.1 | CD4/CD8 scANVI 分析（Notebook） |
| `tcell_cd4cd8_scanvi_retraining_20260109_v1_0.ipynb` | v1.0 | CD4/CD8 重训练（Notebook） |
| `tcell_only_ref_query_merge_pipeline_20260225_v1_3.ipynb` | v1.3 | T 细胞参考/查询合并（Notebook） |
| `tnk_subcluster_annotate_scanvi_retrain_20260318_v1_2_1.ipynb` | v1.2.1 | T/NK 亚群标注与 scANVI 重训练（Notebook） |
| `tnk_scanvi_ref_debug_write_20260315.ipynb` | — | T/NK 参考模型调试（Notebook） |

#### B 细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `b_scvi_celltypist_scanvi_pipeline_v1.py` | v1.0 | ✅ B 细胞主流程：scVI → CellTypist → scANVI |
| `b_bbknn_20260104.py` | v1.0 | B 细胞 BBKNN 批次校正 |
| `step1_train_bcell_L2_20260204_v2_5_3.py` | v2.5.3 | B 细胞 L2 参考模型训练（scArches 用） |
| `step2_map_query_20260204_v2_5_4.py` | v2.5.4 | 查询数据映射到 B 细胞 L2 参考 |

#### 髓系细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `myeloid_scvi_scanvi_v2_4_20260322.py` | **v2.4** 🆕 | ✅ 最新生产版：scArches-ready，scArches dry-run 验证，UMAP 栅格化 |
| `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` | v3.5.1 | 髓系主流程（P0 修复，含 CellTypist） |
| `myeloid_scvi_celltypist_scanvi_pipeline_v1.py` | v1.0 | 髓系初始版本管道 |
| `myeloid_subcluster_analysis_20260121_v2_1.py` | v2.1 | 髓系亚群聚类分析 |
| `stromal_celltype_subcluster_pipeline_20260121_v4_1.py` | v4.1 | 细胞类型亚群详细流程（含基质） |

#### 基质/血管细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py` | v3.5.1 | ✅ 生产推荐：基质/血管主流程（P0 修复） |
| `stromal_vascular_scvi_scanvi_pipeline_20251212_v2.py` | v2.0 | 基质 + 血管联合分析 |
| `stromal_subcluster_analysis_20260121_v2_1.py` | v2.1 | 基质亚群聚类 |
| `stromal_subcluster_analysis_v2_20260121.py` | v2.0 | 基质亚群分析（另一版本） |
| `stromal_pure_scanvi_reanalysis_20260114.py` | v1.0 | 纯 scANVI 基质重分析 |
| `stromal_scarches_query_20260315_v2_1.py` | v2.1 | 基质 scArches 查询映射 |
| `stromal_marker_visualization_20260214_v1_1.py` | v1.1 | 基质标记基因可视化 |
| `stromal_reintegration_20260303_v1_1.ipynb` | v1.1 | 基质重整合（Notebook） |
| `stromal_reintegration_scvi_scanvi_20260312_v1_4.ipynb` | v1.4 | 基质重整合 scVI/scANVI（Notebook） |

#### 上皮 / 基底细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `basal_bbknn_20260104.py` | v1.0 | 基底细胞 BBKNN 批次校正与聚类 |
| `secretory_lineage_bbknn_marker_analysis_20260118_v3.ipynb` | v3.0 | 分泌谱系（Secretory）BBKNN + 标记基因分析（Notebook） |

### 3️⃣ 通用亚群流程

| 脚本 | 版本 | 说明 |
|------|------|------|
| `universal_celltype_subcluster_pipeline_20260121_v4_1.py` | v4.1 | ✅ 推荐：可参数化应用于任意细胞类型的通用亚群流程 |

### 4️⃣ 下游分析 & 可视化

| 脚本 | 版本 | 说明 |
|------|------|------|
| `scanvi_downstream_pipeline_20251218_v1_6_2.py` | v1.6.2 | scANVI 下游全面分析（QC + 可视化） |
| `scanvi_bbknn_analysis_20251217_v1_6.py` | v1.6 | BBKNN 批次校正与下游分析 |
| `visualize_scanvi_results_20260208.py` | v1.0 | scANVI 结果可视化与导出 |

### 5️⃣ 质量控制

| 脚本 | 版本 | 说明 |
|------|------|------|
| `remove_doublets_and_gse299751_20260114.py` | v1.0 | 双细胞检测与移除、特定数据集过滤 |
| `quick_data_check_detailed_20260127.py` | v1.0 | h5ad 结构与元数据详细检查 |

### 6️⃣ 迁移学习 (scArches)

| 脚本 | 版本 | 说明 |
|------|------|------|
| `scarches_mapping_20260127_v1_2_1.py` | v1.2.1 | ✅ 推荐：新数据映射到参考模型（surgical fine-tuning） |
| `step2_scarches_mapping_20260112.py` | v1.0 | scArches 基础映射脚本 |

### 7️⃣ 工具脚本

| 脚本 | 说明 |
|------|------|
| `subset_by_celltype_20260127_v1_0.py` | 从全细胞数据集中按类型提取子集 |
| `step1_add_umap_operator_20260112.py` | 添加或更新 UMAP 嵌入 |
| `sync_counts_from_starcat_to_scanvi.py` | 从 STARcat 输出同步原始计数到 scANVI 数据 |

---

## 🚀 快速开始

### 场景 1：新数据全流程分析

```bash
# Step 1: 全细胞整合 + 自动注释
python allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py

# Step 2: 双细胞检测
python remove_doublets_and_gse299751_20260114.py

# Step 3: 各谱系细化分析（按需选择）
python t_scvi_celltypist_scanvi_20260108_v3_5_1.py                    # T/NK 细胞
python b_scvi_celltypist_scanvi_pipeline_v1.py                         # B 细胞
python myeloid_scvi_scanvi_v2_4_20260322.py                            # 髓系（最新）
python stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py      # 基质/血管

# Step 4: 亚群精细分析
python universal_celltype_subcluster_pipeline_20260121_v4_1.py
```

### 场景 2：将新数据映射到已有参考图谱

```bash
python scarches_mapping_20260127_v1_2_1.py
```

### 场景 3：T/NK 细胞合并管道（含污染细胞过滤）

```bash
python tcell_only_merged_pipeline_20260318_v2_4.py
```

### 场景 4：构建 T/NK scArches 参考模型

```bash
python tnk_scvi_scanvi_scarches_ref_20260315_v1_2.py
```

---

## 📥 数据格式

**输入**: h5ad（AnnData）格式，需包含：
- `layers['counts']` — 原始基因表达计数
- `layers['log1p']` — Log 标准化表达
- `.raw` — 全基因备份
- `obs` 元数据列 — `batch`、`sample`、`cell_type` 等

**输出**: 带注释、嵌入与亚群标签的 h5ad 文件，以及 scVI/scANVI 模型目录（含 `var_names.csv`，供 scArches 使用）

---

## 📈 版本演进

| 系列 | 关键改进 |
|------|---------|
| v2.x (All-cells) | 基因对齐、UMAP 分离、CPU 安全 CUDA、索引对齐写入 |
| v2.3.x | 协变量系统 (MT%, stress, cell cycle)、Unknown 清洗、类别不平衡处理 |
| v3.5.x | P0 修复 (scANVI 标签对齐)、dtype 安全、内存优化 60–70% |
| v4.x | `restrict_to` 保留 BBKNN 图、技术基因过滤、批次效应监控 |
| **v2.4 PRODUCTION** 🆕 | Categorical crash 修复、stress_score 归一化统一、scArches dry-run 验证、var_names.csv 双目录保存、UMAP 栅格化 |

---

## 📝 文件命名规范

```
{celltype}_{analysis}_{YYYYMMDD}_v{X}_{Y}.py

示例：
  allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py
  myeloid_scvi_scanvi_v2_4_20260322.py
  stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py
```

---

## 📊 仓库统计

| 指标 | 数量 |
|------|------|
| Python 脚本 (.py) | 42 |
| Jupyter Notebook (.ipynb) | 11 |
| 总文件数 | 53 |
| 总代码行数 | ~56,000 |
| 覆盖细胞谱系 | 6（上皮/基底、T-NK、B、髓系、基质/血管、联合分析） |
| 最新更新 | 2026-03-22（`myeloid_scvi_scanvi_v2_4_20260322.py`） |

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
