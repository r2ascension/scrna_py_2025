# 🫁 scRNA-seq 肺组织单细胞分析流程 (scrna_py_2025)

> 一套面向人类肺组织的**生产级**双语言（R + Python）单细胞 RNA 测序 (scRNA-seq) 分析管道，覆盖从原始数据预处理到 GPU 加速深度学习整合、细胞类型注释、轨迹分析、解卷积等完整分析链路。

---

## 📖 项目简介

本仓库包含 **139 个 R 脚本**、**58 个 Python 脚本**与 **16 个 Jupyter Notebook**（共 **218 个文件**），构成一套完整的双语言 scRNA-seq 分析工作流，主要用于：

- 多样本/批次人类肺组织单细胞转录组数据预处理与整合（R: Seurat/Harmony；Python: scVI/scANVI）
- 基于 scVI/scANVI 的深度学习批次校正（Python）
- Seurat + Harmony + ROGUE 质量驱动整合（R）
- 使用 CellTypist 进行自动化细胞类型注释（Python）及基于 marker 基因的人工注释（R）
- 各细胞谱系的细粒度亚群鉴定（R: Seurat 亚群 + LLM 解读；Python: scANVI 精炼）
- 基于 scHPL/treeArches 的层级细胞类型校验与新类型发现（Python）
- 差异表达与富集分析（R: MAST + clusterProfiler + GO/KEGG）
- 多谱系组织间比较分析（R: DESeq2 pseudobulk + GSVA + MASC）
- 轨迹与拟时序分析（R: Monocle3）；差异丰度分析（Python: Milo/PAGA）
- SCENIC 基因调控网络推断（R + Python）
- 表型–基因关联多模块管道（R + Python: pseudobulk HPO 富集、监督学习排序、图链接预测、模块评分）
- 基于 BayesPrism 的 bulk RNA 解卷积（R）
- 基于 scArches 的迁移学习，将新数据映射到已有参考图谱（Python）

**数据来源**: 人类肺组织 scRNA-seq（多批次多样本，h5ad / RDS 格式）  
**主要细胞谱系**: 上皮/基底细胞、纤毛/杯状/黏液分泌谱系、T/NK 细胞、B 细胞、髓系细胞、基质/血管/成纤维细胞

---

## ✨ 核心特性

| 特性 | 语言 | 说明 |
|------|------|------|
| **批次校正** | R / Python | Seurat/Harmony；scVI (n_latent=100, n_layers=2)；BBKNN |
| **自动注释** | Python | CellTypist 多数投票（Human_Lung_Atlas + Immune_All_Low） |
| **半监督精炼** | Python | scANVI 标签引导细化 |
| **亚群聚类** | R / Python | Seurat Leiden；restrict_to 保留 BBKNN 图 |
| **LLM 亚群解读** | R | `interpret()` 函数驱动 AI 自动注释亚群 |
| **差异表达** | R | MAST 回归框架 + DESeq2 pseudobulk |
| **富集分析** | R | clusterProfiler（GO/KEGG/GSVA）、msigdbdf |
| **轨迹分析** | R | Monocle3 multi-tissue trajectory |
| **质量评估** | R | ROGUE 同质性评分、DoubletFinder、DecontX |
| **解卷积** | R | BayesPrism bulk RNA 解卷积 |
| **双细胞检测** | Python | Top2 谱系 + margin 逻辑 |
| **迁移学习** | Python | scArches surgical fine-tuning |
| **层级校验** | Python | scHPL/treeArches 层级细胞类型校验（global lineage-wise + branch-wise） |
| **组织比较** | R | DESeq2 pseudobulk + GSVA + MASC 多谱系组织间差异分析框架 |
| **SCENIC/GRN** | R / Python | SCENIC 基因调控网络推断（cisTarget 数据库 + R 核心引擎） |
| **差异丰度** | Python | Milo/PAGA 邻域差异丰度分析 |
| **表型–基因关联** | R / Python | HPO 富集 + 监督学习排序 + 图链接预测 + 模块评分 |
| **内存优化** | Python | HVG 训练减少 60–70% 内存占用 |
| **GPU 加速** | Python | 单卡 CUDA 训练，CPU 安全回退 |
| **完全可重复** | R / Python | 固定随机种子 + HVG 列表保存 + scArches dry-run 验证 |

---

## 🛠️ 依赖环境

### R 环境

```r
# 核心单细胞分析
Seurat (≥5.0)
SeuratDisk
SingleCellExperiment

# 批次校正 & 质量评估
harmony
ROGUE
bbknnR
DoubletFinder
DecontX (celda)

# 差异表达
MAST
DESeq2
limma
speckle

# 富集分析
clusterProfiler
org.Hs.eg.db
KEGGREST
msigdbdf
enrichit

# 轨迹分析
monocle3

# 解卷积
BayesPrism

# 格式转换
SCNT          # Seurat → H5AD
DropletUtils  # 10x H5 处理
SeuratData
reticulate

# LLM 辅助注释
# interpret() 函数（自定义封装）

# 可视化
ggplot2
ggalluvial
ggrastr
ComplexHeatmap
pheatmap
cowplot
```

> **推荐安装**：
> ```r
> install.packages(c("Seurat","harmony","ggplot2","cowplot","ggrastr"))
> BiocManager::install(c("MAST","DESeq2","clusterProfiler","org.Hs.eg.db","monocle3"))
> remotes::install_github("PaulingLiu/ROGUE")
> install_github("Danko-Lab/BayesPrism/BayesPrism")
> ```

### Python 环境

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
| `tcell_merge_schpl_20260402_v1.py` | v1.0 | 🆕 T/NK scHPL 层级校验与合并 |
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
| `bcell_merge_schpl_20260412_v1.py` | v1.0 | 🆕 B 细胞 scHPL 层级校验与合并 |
| `bcell_paga_milo_20260412_v1_1.py` | v1.1 | 🆕 B 细胞 PAGA 轨迹 + Milo 差异丰度分析 |
| `bcell_scanvi_validation_20260412.py` | v1.0 | 🆕 B 细胞 scANVI 注释验证 |
| `bcell_schpl_reject_followup_20260412_v1.py` | v1.0 | 🆕 B 细胞 scHPL Rejected 细胞跟进分析 |
| `prepare_bcell_scanvi_umap_input_20260413.py` | v1.0 | 🆕 B 细胞 scANVI UMAP 输入准备 |

#### 髓系细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `myeloid_scvi_scanvi_v2_4_20260322.py` | **v2.4** 🆕 | ✅ 最新生产版：scArches-ready，scArches dry-run 验证，UMAP 栅格化 |
| `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` | v3.5.1 | 髓系主流程（P0 修复，含 CellTypist） |
| `myeloid_scvi_celltypist_scanvi_pipeline_v1.py` | v1.0 | 髓系初始版本管道 |
| `myeloid_subcluster_analysis_20260121_v2_1.py` | v2.1 | 髓系亚群聚类分析 |
| `stromal_celltype_subcluster_pipeline_20260121_v4_1.py` | v4.1 | 细胞类型亚群详细流程（含基质） |
| `myeloid_merge_schpl_20260402_v1.py` | v1.0 | 🆕 髓系 scHPL 层级校验与合并 |
| `myeloid_milopy_tissue_celltype_20260410_v1.py` | v1.0 | 🆕 髓系 Milo 邻域差异丰度（按组织 × 细胞类型） |

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
| `stromal_reintegration_branchwise_scvi_scanvi_20260407_v1_5.py` | v1.5 | 🆕 基质分支式（branchwise）scVI/scANVI 重整合 |
| `stromal_scarches_query_branchwise_20260407_v1_0.py` | v1.0 | 🆕 基质分支式 scArches 查询映射 |
| `stromal_schpl_reject_followup_20260408_v1_0.py` | v1.0 | 🆕 基质 scHPL Rejected 细胞跟进分析 |
| `stromal_reintegration_20260303_v1_1.ipynb` | v1.1 | 基质重整合（Notebook） |
| `stromal_reintegration_scvi_scanvi_20260312_v1_4.ipynb` | v1.4 | 基质重整合 scVI/scANVI（Notebook） |

#### 上皮 / 基底细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `basal_bbknn_20260104.py` | v1.0 | 基底细胞 BBKNN 批次校正与聚类 |
| `epithelial_scarches_query_mapping_20260315_v1_1.py` | v1.1 | 🆕 上皮细胞 scArches 查询映射 |
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

### 8️⃣ scHPL / treeArches 层级校验 🆕

**功能**: 在 base model（scANVI/scArches）之后作为独立层级校验器，判断 query 细胞是否被 reference hierarchy 稳定吸收

| 脚本 | 版本 | 说明 |
|------|------|------|
| `lineage_merge_schpl_core_20260402.py` | v1.0 | ✅ global lineage-wise 复用核心（kNN classifier, dimred=True） |
| `bcell_merge_schpl_20260412_v1.py` | v1.0 | B 细胞 scHPL 谱系合并 |
| `tcell_merge_schpl_20260402_v1.py` | v1.0 | T/NK scHPL 谱系合并 |
| `myeloid_merge_schpl_20260402_v1.py` | v1.0 | 髓系 scHPL 谱系合并 |
| `bcell_schpl_reject_followup_20260412_v1.py` | v1.0 | B 细胞 Rejected 细胞跟进分析 |
| `stromal_schpl_reject_followup_20260408_v1_0.py` | v1.0 | 基质 Rejected 细胞跟进分析 |
| `bcell_merge_schpl_20260330_v3_0.ipynb` | v3.0 | B 细胞 scHPL 合并（Notebook） |
| `stromal_schpl_treeArches_20260330_v1_0.ipynb` | v1.0 | 基质 treeArches branch-wise 主流程（Notebook） |
| `myeloid_L3refined_patch_20260329_v1_0.ipynb` | v1.0 | 髓系 L3 精修补丁（Notebook） |

> 方法学详见 **[SCHPL_METHODOLOGY_UNIFIED_20260410.md](./SCHPL_METHODOLOGY_UNIFIED_20260410.md)**

### 9️⃣ 表型–基因关联管道 🆕

**功能**: 多模块管道，从 pseudobulk DE 到图神经网络链接预测，系统性关联表型与基因

| 脚本 | 模块 | 说明 |
|------|------|------|
| `phenotype_gene_pipeline_A_pseudobulk_hpo_enrichment.R` | A | Pseudobulk + HPO 表型富集 |
| `phenotype_gene_pipeline_B_supervised_ranking.py` | B | 监督学习基因排序 |
| `phenotype_gene_pipeline_C_graph_link_prediction.py` | C | 图链接预测（基因–表型关联） |
| `phenotype_gene_pipeline_D_module_score_hpo.R` | D | 模块评分 + HPO 映射 |

### 🔟 SCENIC 基因调控网络 🆕

**功能**: 基于 SCENIC 框架推断转录因子调控网络

| 脚本 | 说明 |
|------|------|
| `prepare_rscenic_cistarget_db_20260410.py` | cisTarget 数据库准备（Python 下载 + 格式化） |
| `configure_scenic_r_env_20260410.R` | SCENIC R 环境配置 |
| `scenic_core_20260410.R` | ✅ SCENIC 核心分析引擎（GRNBoost2 + cisTarget + AUCell） |
| `tnk_tcell_scenic_L3_20260410.R` | T/NK L3 精细级别 SCENIC 分析 |
| `bcell_scenic_L3_20260412.R` | B 细胞 L3 精细级别 SCENIC 分析 |

---

## 📂 R 脚本目录（107 个）

> R 脚本覆盖从 RDS 读取、QC、整合到差异表达、轨迹、解卷积的完整分析链路。

### R-1️⃣ 环境配置

| 脚本 | 说明 |
|------|------|
| `Installation.R` | 一键安装所有 R 依赖（BiocManager / devtools / remotes） |

### R-2️⃣ 数据读取与格式转换

| 脚本 | 说明 |
|------|------|
| `GetSeurat.R` | 创建/加载 Seurat 对象 |
| `readRDS.R` | 批量读取 RDS 文件 |
| `H5To10x.R` | H5 格式 → 10x Cell Ranger 格式转换 |
| `ProcessH5AD.R` | H5AD 文件处理与 AnnData ↔ Seurat 互转 |
| `ProcessBatch.r` | 批量数据处理入口 |
| `rds_to_h5ad.R` | Seurat RDS → H5AD 批量转换（SCNT::GetH5ad） |
| `create_seurat_from_csv_20251214_v1.R` | 从 CSV/矩阵文件构建 Seurat 对象 |
| `gsm_one_by_one_debug.R` | 逐 GSM 样本调试加载 |

### R-3️⃣ QC、双细胞检测与智能合并

| 脚本 | 版本 | 说明 |
|------|------|------|
| `rds_folder_qc_smartmerge_20251216_v4_memory_optimized.R` | v4 | ✅ 内存优化：惰性加载 + 增量合并 + DoubletFinder + DecontX |
| `rds_folder_qc_smartmerge_20251216_v5.R` | v5 | 最新智能合并管道 |
| `filter_dataset_rds.R` | — | RDS 数据集过滤 |
| `filter_dataset_rds_20251227.R` | — | 过滤脚本更新版 |

### R-4️⃣ 基因标准化与整合

| 脚本 | 版本 | 说明 |
|------|------|------|
| `gene_standardization_merge_20251221_v4_2.R` | v4.2 | 多批次基因名标准化合并 |
| `merge_stage_gene_standardization_20251218.R` | — | 合并阶段基因标准化 |
| `merge_stage_gene_standardization_single_20251220_v2.R` | v2 | 单样本基因标准化 |
| `gene_mapping_batch_fix.R` | — | 基因 ID 映射批次修复 |
| `RenameENSG.R` | v1 | ENSEMBL ID → Gene Symbol 重命名 |
| `RenameENSG_v2.R` | v2 | ENSG 重命名更新版 |
| `Integrating-1011.R` | — | 早期多批次整合（2024-10） |
| `Integrating_1124.R` | — | 整合管道（2024-11） |
| `Intergrated-1003.R` | — | 整合初始版本（2024-10） |

### R-5️⃣ 细胞类型注释

| 脚本 | 说明 |
|------|------|
| `CellMarker.R` | 基于 CellMarker 数据库的 marker 基因匹配 |
| `Marker.R` | 自定义 marker 基因分析与打分 |
| `SubsetAnalysis.R` | 基于注释结果的子集提取与分析 |
| `intersection.R` | 跨数据集基因集交集分析 |

### R-6️⃣ 上皮细胞系分析

#### 综合分析

| 脚本 | 版本 | 说明 |
|------|------|------|
| `Epithelial_0316.r` | — | 早期上皮分析（2024-03） |
| `all_epithelial_20251222.R` | — | 全上皮细胞整合分析 |
| `epithelial_20251217.R` | — | 上皮细胞分析（2025-12） |
| `epithelial_celltype_harmony_rogue_20251220_v3_1.R` | v3.1 | ✅ Harmony + ROGUE 同质性评估 |
| `epithelial_harmony_rogue_analysis_20251220_v2_6.R` | v2.6 | Harmony + ROGUE 分析 |
| `epithelial_rogue_20251218.R` | — | ROGUE 质量评估 |
| `epithelial_subset_rogue_20251218.R` | — | 亚集 ROGUE 分析 |
| `harmony_subset_rogue_v1.2.R` | v1.2 | Harmony 亚集 ROGUE |
| `harmony_subset_rogue_20251219_v1_4.R` | v1.4 | Harmony 亚集 ROGUE 更新版 |
| `tissue_comparison_analysis_20251222.R` | v3.0 | 多组织比较（DESeq2 + GSVA + pseudobulk） |

#### 亚群解读（LLM 辅助）

| 脚本 | 版本 | 说明 |
|------|------|------|
| `epithelial_subcluster_interpret_analysis_20260128_v2_7.R` | v2.7 | 亚群 LLM 解读 |
| `epithelial_subcluster_interpret_analysis_20260205_v2_8.R` | v2.8 | 亚群解读更新版 |
| `epithelial_subcluster_interpret_20260209_v3_1_EPITHELIAL copy.R` | v3.1 | v3.1 生产版 |
| `epithelial_subcluster_interpret_20260210_v4.R` | v4 | ✅ 最新：上皮亚群解读 v4 |

#### 组织间比较 🆕

| 脚本 | 版本 | 说明 |
|------|------|------|
| `epithelial_tissue_comparison_v1_3_2_20260413.R` | v1.3.2 | ✅ 上皮组织间比较最新版 |
| `epithelial_tissue_comparison_20260407_v1_3_1.R` | v1.3.1 | 上皮组织比较 v1.3.1 |
| `epithelial_tissue_comparison_20260407_v1_3.R` | v1.3 | 上皮组织比较 v1.3 |
| `epithelial_tissue_comparison_20260402_v1_2.R` | v1.2 | 上皮组织比较 v1.2 |
| `epithelial_tissue_comparison_20260331_v1_1.R` | v1.1 | 上皮组织比较 v1.1 |

#### 轨迹分析

| 脚本 | 版本 | 说明 |
|------|------|------|
| `epithelial_monocle3_by_tissue.R` | — | Monocle3 按组织分类轨迹 |
| `epithelial_monocle3_trajectory_20251220_v2_7.R` | v2.7 | ✅ Monocle3 多组织拟时序轨迹 |

### R-7️⃣ 基底细胞（Basal）

| 脚本 | 版本 | 说明 |
|------|------|------|
| `basal_20251221.R` | — | 基底细胞分析（2025-12） |
| `basal_20251228.R` | — | 基底细胞分析更新 |
| `basal_20260105.R` | — | 基底细胞分析（2026-01） |
| `basal_multitissue_trajectory_20260103_v2_8.R` | v2.8 | ✅ 多组织基底细胞轨迹分析 |

### R-8️⃣ 分泌上皮亚型（Ciliated / Club / Goblet / Serous）

| 脚本 | 版本 | 说明 |
|------|------|------|
| `ciliated_20251222.R` | — | 纤毛细胞分析 |
| `ciliated_20260103.R` | — | 纤毛细胞（2026-01） |
| `ciliated_20260103_2.R` | — | 纤毛细胞分析（第二版） |
| `ciliated_20260104_3.R` | — | 纤毛细胞分析（第三版） |
| `ciliated_bbknnR_native_analysis_20251230.R` | — | bbknnR 原生 BBKNN 分析 |
| `ciliated_bbknn_wilcox_analysis_20251230.R` | — | BBKNN + Wilcox 差异分析 |
| `club_20251230.R` | — | Club 细胞分析 |
| `goblet_20251222.R` | — | 杯状细胞分析 |
| `serous_20251221.R` | — | 浆液腺细胞分析 |

### R-9️⃣ T/NK 细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `T_0317.r` | — | 早期 T 细胞分析（2024-03） |
| `tcell_subcluster_interpret_analysis_20260204_v1_2.R` | v1.2 | T 细胞亚群 LLM 解读 v1.2 |
| `tcell_subcluster_interpret_20260209_v3_0.R` | v3.0 | ✅ T 细胞亚群解读 v3.0 |
| `tnk_tissue_comparison_v2_6_1_20260413.R` | v2.6.1 | 🆕 ✅ T/NK 组织间比较最新版 |
| `tnk_tissue_comparison_v2_6_20260407.R` | v2.6 | T/NK 组织比较 v2.6 |
| `tnk_tissue_comparison_v2_4_20260401.R` | v2.4 | T/NK 组织比较 v2.4 |
| `tnk_tcell_scenic_L3_20260410.R` | v1.0 | 🆕 T/NK L3 SCENIC 基因调控网络分析 |
| `tnk_llm_focus_dotplot_20260413.R` | v1.0 | 🆕 T/NK LLM 聚焦 DotPlot 可视化 |

### R-🔟 B 细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `B_0626.r` | — | 早期 B 细胞分析（2024-06） |
| `B_20260104.R` | — | B 细胞分析（2026-01） |
| `bcell_subcluster_interpret_analysis_20260126_v1_0.R` | v1.0 | B 细胞亚群解读初始版 |
| `bcell_subcluster_interpret_analysis_20260126_v1_1.R` | v1.1 | B 细胞亚群解读 v1.1 |
| `bcell_subcluster_interpret_analysis_20260126_v1_2.R` | v1.2 | B 细胞亚群解读 v1.2 |
| `bcell_subcluster_interpret_analysis_20260127_v2.R` | v2.0 | B 细胞亚群解读 v2 |
| `bcell_subcluster_interpret_analysis_20260127_v2_1.R` | v2.1 | B 细胞亚群解读 v2.1 |
| `bcell_interpret_PRODUCTION_v2_0_20260127.R` | v2.0 | ✅ 生产推荐：B 细胞 LLM 解读（QUICK_REFERENCE_MEMORY v2.13 合规） |
| `bcell_tissue_comparison_v2_6_4_20260412.R` | v2.6.4 | 🆕 ✅ B 细胞组织间比较最新版（DESeq2 + GSVA + MASC） |
| `bcell_tissue_comparison_v2_6_3_20260411.R` | v2.6.3 | B 细胞组织比较 v2.6.3 |
| `bcell_tissue_comparison_v2_6_2_20260410.R` | v2.6.2 | B 细胞组织比较 v2.6.2 |
| `bcell_tissue_comparison_v2_6_1_20260410.R` | v2.6.1 | B 细胞组织比较 v2.6.1 |
| `bcell_tissue_comparison_v2_6_20260410.R` | v2.6 | B 细胞组织比较 v2.6 |
| `bcell_tissue_comparison_v2_6_20260410_llm_resume.R` | v2.6 | B 细胞组织比较（LLM 续跑版） |
| `bcell_tissue_comparison_v2_6_20260406.R` | v2.6 | B 细胞组织比较 v2.6（早期） |
| `bcell_tissue_comparison_v2_5_20260402.R` | v2.5 | B 细胞组织比较 v2.5 |
| `bcell_tissue_comparison_v2_4_20260390.R` | v2.4 | B 细胞组织比较 v2.4 |
| `bcell_scenic_L3_20260412.R` | v1.0 | 🆕 B 细胞 L3 SCENIC 基因调控网络分析 |

### R-1️⃣1️⃣ 髓系细胞（Myeloid）

| 脚本 | 版本 | 说明 |
|------|------|------|
| `Myeloid_0317.r` | — | 早期髓系分析（2024-03） |
| `Myeloid_20251228.R` | — | 髓系分析（2025-12） |
| `myeloid_subcluster_interpret_analysis_20260131_v2_6.R` | v2.6 | 髓系亚群解读 v2.6 |
| `myeloid_subcluster_interpret_DUAL_20260205_v2_7.R` | v2.7 | 髓系双模式（DUAL）解读 v2.7 |
| `myeloid_subcluster_interpret_20260209_v3_0.R` | v3.0 | ✅ 髓系亚群解读 v3.0 |
| `myeloid_tissue_comparison_20260407_v1_1.R` | v1.1 | 🆕 ✅ 髓系组织间比较最新版 |
| `myeloid_tissue_comparison_20260401_v1_0.R` | v1.0 | 髓系组织比较 v1.0 |

### R-1️⃣2️⃣ 基质 / 血管 / 成纤维细胞

| 脚本 | 版本 | 说明 |
|------|------|------|
| `Fibroblast_0318.r` | — | 早期成纤维细胞分析（2024-03） |
| `SMC_0627.r` | — | 平滑肌细胞分析（2024-06） |
| `Endothelial_0627.r` | — | 内皮细胞分析（2024-06） |
| `stromal_20251216.R` | — | 基质细胞整合分析 |
| `stromal_subcluster_interpret_20260209_v3_1_STROMAL.R` | v3.1 | 基质亚群解读 v3.1 |
| `stromal_subcluster_interpret_20260210_v4.R` | v4 | ✅ 基质亚群解读 v4（PRODUCTION） |
| `stromal_endothelial_tissue_comparison_20260408_v1_0.R` | v1.0 | 🆕 内皮细胞亚型组织间比较 |
| `stromal_fibroblast_tissue_comparison_20260408_v1_0.R` | v1.0 | 🆕 成纤维细胞亚型组织间比较 |
| `stromal_smc_tissue_comparison_20260408_v1_0.R` | v1.0 | 🆕 平滑肌细胞亚型组织间比较 |
| `stromal_schpl_postprocess_v1_0_20260401.R` | v1.0 | 🆕 基质 scHPL 输出后处理（R 侧下游整理） |

### R-1️⃣3️⃣ 差异表达（MAST）

| 脚本 | 说明 |
|------|------|
| `MAST_pipeline_20251223.R` | ✅ 通用 MAST DEG 管道（批次 NMF + 稀疏优化） |
| `MAST_pipeline_single_h5ad_AT1_20251224.R` | AT1 细胞 MAST 管道 |
| `MAST_pipeline_single_h5ad_AT2_20251227.R` | AT2 细胞 MAST 管道 |
| `MAST_pipeline_single_h5ad_B_20251225.R` | B 细胞 MAST 管道 |
| `MAST_pipeline_single_h5ad_Myeloid_20251225.R` | 髓系 MAST 管道 |
| `MAST_pipeline_single_h5ad_Stromal_20251227_0.R` | 基质 MAST 管道 |
| `MAST_pipeline_single_h5ad_T_20251225.R` | T 细胞 MAST 管道 |
| `MAST_pipeline_single_h5ad_basal_20251225.R` | 基底细胞 MAST 管道 |
| `MAST_pipeline_single_h5ad_ciliated_20251227.R` | 纤毛细胞 MAST 管道 |
| `MAST_pipeline_single_h5ad_club_20251227.R` | Club 细胞 MAST 管道 |

### R-1️⃣4️⃣ GO/KEGG 富集分析

| 脚本 | 版本 | 说明 |
|------|------|------|
| `h5ad_mast_go_kegg_single_B_20251228_v1_2.R` | v1.2 | B 细胞 MAST + GO + KEGG 联合分析 |
| `h5ad_mast_go_kegg_single_basal_20251228_v1_1.R` | v1.1 | 基底细胞 MAST + GO + KEGG |
| `h5ad_mast_go_kegg_single_ciliated_20251228_v1_1.R` | v1.1 | 纤毛细胞 MAST + GO + KEGG |
| `h5ad_mast_go_kegg_single_club_20251230_v1_1.R` | v1.1 | Club 细胞 MAST + GO + KEGG |
| `enrichment_functions.R` | — | 富集分析公共函数库 |
| `pseudobulk.R` | — | Pseudobulk 差异表达框架 |
| `RdsLimma.R` | — | 基于 limma 的 RDS 差异分析 |
| `scMASC.R` | — | scMASC 注释辅助分析 |

### R-1️⃣4️⃣-b 组织间比较框架 🆕

| 脚本 | 版本 | 说明 |
|------|------|------|
| `tissue_comparison_generic_wrapper_20260412.R` | v1.0 | ✅ 通用组织比较封装（一键调用任意谱系比较） |
| `tissue_comparison_advanced_helper_20260408.R` | v1.0 | 组织比较高级辅助函数（热图、富集联合可视化） |
| `tissue_comparison_analysis_20251222.R` | v3.0 | 多组织比较基础版（DESeq2 + GSVA + pseudobulk） |

### R-1️⃣4️⃣-c SCENIC 基因调控网络 🆕

| 脚本 | 说明 |
|------|------|
| `configure_scenic_r_env_20260410.R` | SCENIC R 环境与依赖配置 |
| `scenic_core_20260410.R` | ✅ SCENIC 核心分析引擎（GRNBoost2 + cisTarget + AUCell） |
| `tnk_tcell_scenic_L3_20260410.R` | T/NK L3 SCENIC 分析 |
| `bcell_scenic_L3_20260412.R` | B 细胞 L3 SCENIC 分析 |

### R-1️⃣4️⃣-d 表型–基因关联管道 🆕

| 脚本 | 模块 | 说明 |
|------|------|------|
| `phenotype_gene_pipeline_A_pseudobulk_hpo_enrichment.R` | A | Pseudobulk + HPO 表型富集 |
| `phenotype_gene_pipeline_D_module_score_hpo.R` | D | 模块评分 + HPO 映射 |

### R-1️⃣5️⃣ 解卷积（BayesPrism）

| 脚本 | 版本 | 说明 |
|------|------|------|
| `deconvolution.R` | — | 解卷积基础脚本 |
| `bayesprism_deconvolution_production_20260124.R` | v2.1 | ✅ BayesPrism bulk RNA 解卷积生产版（单/双层注释） |
| `bayesprism_deconvolution_production_20260303.R` | — | BayesPrism 解卷积（2026-03） |
| `bayesprism_deconvolution_production_20260308.R` | — | BayesPrism 解卷积（2026-03 更新） |
| `bayesprism_quick_recovery_20260126.R` | — | BayesPrism 快速恢复/续跑脚本 |

### R-1️⃣6️⃣ LLM 亚群解读工具

| 脚本 | 说明 |
|------|------|
| `interpret.R` | 核心 LLM interpret() 函数封装 |
| `interpret_agent_hotfix.R` | LLM 解读代理热修复 |

### R-1️⃣7️⃣ 疾病特异性分析（鼻息肉 / SARS）

| 脚本 | 版本 | 说明 |
|------|------|------|
| `PolypDiagnosis20251212.R` | — | 鼻息肉诊断特征分析 |
| `PolypRead20251213.R` | — | 鼻息肉数据读取 |
| `PolypRead_20251218_v2_3.R` | v2.3 | 鼻息肉读取更新版 |
| `PolypRead_20251218_v2_3_SARS.R` | v2.3 | SARS 相关数据版本 |

### R-1️⃣8️⃣ 其他工具

| 脚本 | 说明 |
|------|------|
| `test_quick_20260127_v_2.R` | 快速测试 / 调试脚本 |

---

## 🚀 快速开始

### Python 工作流

#### 场景 1：新数据全流程分析（Python）

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

#### 场景 2：将新数据映射到已有参考图谱

```bash
python scarches_mapping_20260127_v1_2_1.py
```

#### 场景 3：T/NK 细胞合并管道（含污染细胞过滤）

```bash
python tcell_only_merged_pipeline_20260318_v2_4.py
```

#### 场景 4：构建 T/NK scArches 参考模型

```bash
python tnk_scvi_scanvi_scarches_ref_20260315_v1_2.py
```

#### 场景 4b：scHPL 层级校验（新数据 vs 参考）🆕

```bash
# 通用 lineage-wise 核心
python lineage_merge_schpl_core_20260402.py

# 谱系特异性
python bcell_merge_schpl_20260412_v1.py     # B 细胞
python tcell_merge_schpl_20260402_v1.py     # T/NK
python myeloid_merge_schpl_20260402_v1.py   # 髓系

# Rejected 细胞跟进
python bcell_schpl_reject_followup_20260412_v1.py
```

#### 场景 4c：B 细胞 PAGA + Milo 差异丰度分析 🆕

```bash
python bcell_paga_milo_20260412_v1_1.py
```

### R 工作流

#### 场景 5：RDS 批量 QC + 智能合并

```r
source("rds_folder_qc_smartmerge_20251216_v5.R")
```

#### 场景 6：上皮细胞 Harmony + ROGUE 整合与轨迹分析

```r
source("epithelial_celltype_harmony_rogue_20251220_v3_1.R")   # 整合 + 质量
source("epithelial_monocle3_trajectory_20251220_v2_7.R")       # Monocle3 轨迹
```

#### 场景 7：各谱系亚群 LLM 解读

```r
source("interpret.R")                                          # 加载解读函数
source("bcell_interpret_PRODUCTION_v2_0_20260127.R")          # B 细胞
source("myeloid_subcluster_interpret_20260209_v3_0.R")        # 髓系
source("stromal_subcluster_interpret_20260210_v4.R")          # 基质/血管
source("tcell_subcluster_interpret_20260209_v3_0.R")          # T 细胞
source("epithelial_subcluster_interpret_20260210_v4.R")       # 上皮
```

#### 场景 8：差异表达 + GO/KEGG 富集

```r
source("MAST_pipeline_20251223.R")                            # 批量 MAST DEG
source("h5ad_mast_go_kegg_single_B_20251228_v1_2.R")         # B 细胞 MAST+GO+KEGG
```

#### 场景 9：Bulk RNA 解卷积（BayesPrism）

```r
source("bayesprism_deconvolution_production_20260124.R")
```

#### 场景 10：RDS → H5AD 格式转换（R → Python 桥接）

```r
source("rds_to_h5ad.R")   # 批量 Seurat RDS → H5AD，供 Python 管道使用
```

#### 场景 11：全谱系组织间比较分析 🆕

```r
# 通用封装器（可参数化调用任意谱系）
source("tissue_comparison_generic_wrapper_20260412.R")

# 或按谱系单独调用
source("bcell_tissue_comparison_v2_6_4_20260412.R")              # B 细胞
source("tnk_tissue_comparison_v2_6_1_20260413.R")                # T/NK
source("myeloid_tissue_comparison_20260407_v1_1.R")              # 髓系
source("epithelial_tissue_comparison_v1_3_2_20260413.R")         # 上皮
source("stromal_endothelial_tissue_comparison_20260408_v1_0.R")  # 内皮
source("stromal_fibroblast_tissue_comparison_20260408_v1_0.R")   # 成纤维
source("stromal_smc_tissue_comparison_20260408_v1_0.R")          # 平滑肌
```

#### 场景 12：SCENIC 基因调控网络推断 🆕

```r
# Step 1: 配置环境
source("configure_scenic_r_env_20260410.R")

# Step 2: 核心 SCENIC 分析
source("scenic_core_20260410.R")

# Step 3: 谱系特异性 GRN（按需选择）
source("tnk_tcell_scenic_L3_20260410.R")   # T/NK L3
source("bcell_scenic_L3_20260412.R")       # B 细胞 L3
```

#### 场景 13：表型–基因关联多模块管道 🆕

```r
# Module A: Pseudobulk HPO 富集
source("phenotype_gene_pipeline_A_pseudobulk_hpo_enrichment.R")

# Module D: 模块评分 + HPO 映射
source("phenotype_gene_pipeline_D_module_score_hpo.R")
```

```bash
# Module B: 监督学习基因排序
python phenotype_gene_pipeline_B_supervised_ranking.py

# Module C: 图链接预测
python phenotype_gene_pipeline_C_graph_link_prediction.py
```

---

## 📥 数据格式

**Python 管道输入** (h5ad / AnnData)：
- `layers['counts']` — 原始基因表达计数
- `layers['log1p']` — Log 标准化表达
- `.raw` — 全基因备份
- `obs` 元数据列 — `batch`、`sample`、`cell_type` 等

**R 管道输入** (RDS / Seurat)：
- Seurat 对象（v4/v5），含 `RNA` assay 的 `counts` 和 `data` 层
- 可通过 `rds_to_h5ad.R` 转换为 h5ad 格式供 Python 管道使用

**输出**：
- 带注释、嵌入与亚群标签的 h5ad 文件
- scVI/scANVI 模型目录（含 `var_names.csv`，供 scArches 使用）
- 带 Harmony/ROGUE 嵌入的 Seurat RDS 文件（R 管道）

---

## 📈 版本演进

| 系列 | 关键改进 |
|------|---------|
| v2.x (All-cells) | 基因对齐、UMAP 分离、CPU 安全 CUDA、索引对齐写入 |
| v2.3.x | 协变量系统 (MT%, stress, cell cycle)、Unknown 清洗、类别不平衡处理 |
| v3.5.x | P0 修复 (scANVI 标签对齐)、dtype 安全、内存优化 60–70% |
| v4.x | `restrict_to` 保留 BBKNN 图、技术基因过滤、批次效应监控 |
| **v2.4 PRODUCTION** | Categorical crash 修复、stress_score 归一化统一、scArches dry-run 验证、var_names.csv 双目录保存、UMAP 栅格化 |
| **scHPL/treeArches** 🆕 | 层级细胞类型校验、global lineage-wise + branch-wise 双模式、标准 Rejected follow-up 框架 |
| **Tissue Comparison v2.6** 🆕 | 全谱系组织间比较框架（B/T-NK/Myeloid/Epithelial/Stromal 各亚型）、通用封装器、LLM 续跑支持 |
| **SCENIC/GRN** 🆕 | SCENIC 基因调控网络推断管道（cisTarget DB 准备 + GRNBoost2 + AUCell）、L3 精细级别分析 |
| **Phenotype–Gene** 🆕 | 四模块管道（pseudobulk HPO 富集 → 监督学习排序 → 图链接预测 → 模块评分） |
| **Milo/PAGA** 🆕 | 邻域差异丰度分析（myeloid_milopy）、PAGA 轨迹 + Milo 联合（bcell_paga_milo） |

---

## 📝 文件命名规范

```
# Python
{celltype}_{analysis}_{YYYYMMDD}_v{X}_{Y}.py

# R
{celltype}_{analysis}_{YYYYMMDD}_v{X}_{Y}.R   # 命名规则一致

示例：
  allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py
  myeloid_scvi_scanvi_v2_4_20260322.py
  bcell_interpret_PRODUCTION_v2_0_20260127.R
  epithelial_subcluster_interpret_20260210_v4.R
```

---

## 📊 仓库统计

| 指标 | 数量 |
|------|------|
| R 脚本 (.R / .r) | 139 |
| Python 脚本 (.py) | 58 |
| Jupyter Notebook (.ipynb) | 16 |
| Markdown 文档 (.md) | 5 |
| 总文件数 | 218 |
| 总代码行数 | ~150,000+ |
| 覆盖细胞谱系 | 8（上皮/基底/纤毛/分泌、T-NK、B、髓系、基质/血管、联合分析） |
| 主要分析框架 | R: Seurat + Harmony + MAST + Monocle3 + BayesPrism + SCENIC；Python: scVI + scANVI + CellTypist + scArches + scHPL/treeArches + Milo |
| 最新更新 | 2026-04-13（`tnk_tissue_comparison_v2_6_1_20260413.R`、`epithelial_tissue_comparison_v1_3_2_20260413.R`、`prepare_bcell_scanvi_umap_input_20260413.py`） |

---

## 📄 文档

- **[SCRIPT_CATEGORIZATION.md](./SCRIPT_CATEGORIZATION.md)** — 全部脚本详细分类、推荐版本与快速查找指南
- **[SCHPL_METHODOLOGY_UNIFIED_20260410.md](./SCHPL_METHODOLOGY_UNIFIED_20260410.md)** — scHPL / treeArches 统一方法学说明（global lineage-wise vs branch-wise）
- **[DENORM_INTEGRATION_SUMMARY.md](./DENORM_INTEGRATION_SUMMARY.md)** — 去对数归一化集成说明（DecontX 兼容性）
- **[CLAUDE.md](./CLAUDE.md)** — AI 辅助编码指引（分析框架、关键函数、工作流参考）

---

## 🔗 相关资源

**Python**
- [scVI-tools 文档](https://scvi-tools.org/)
- [CellTypist 文档](https://celltypist.readthedocs.io/)
- [Scanpy 文档](https://scanpy.readthedocs.io/)
- [scArches 文档](https://scarches.readthedocs.io/)

**R**
- [Seurat 文档](https://satijalab.org/seurat/)
- [Harmony 文档](https://portals.broadinstitute.org/harmony/)
- [ROGUE 文档](https://github.com/PaulingLiu/ROGUE)
- [Monocle3 文档](https://cole-trapnell-lab.github.io/monocle3/)
- [BayesPrism 文档](https://github.com/Danko-Lab/BayesPrism)
- [clusterProfiler 文档](https://bioconductor.org/packages/clusterProfiler/)

**数据集**
- [Human Lung Cell Atlas](https://www.humancellatlas.org/)
