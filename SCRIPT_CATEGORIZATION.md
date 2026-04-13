# 单细胞分析脚本分类总结

## 📊 分类概览

本文档对 `/home/h2048/script/py` 目录下的较新版本脚本进行分类整理，便于快速定位和使用。

---

## 1️⃣ 全细胞整合分析 (All-cells Integration)

**功能**: 对所有细胞类型进行批次校正和注释的完整流程

| 脚本名称 | 版本 | 日期 | 特性 |
|---------|------|------|------|
| `allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py` | **v2.3.1** | 2026-01-15 | ✅ **推荐** - 协变量系统、Unknown清洗、类别不平衡处理 |
| `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5.py` | v2.5 | 2026-01-21 | 最新实验版本 |
| `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py` | v2.5.2 | 2026-01-21 | 最新实验版本（修订） |

**核心技能**:
- scVI批次校正 (n_latent=100, n_layers=2)
- CellTypist自动注释 (Human_Lung_Atlas模型)
- scANVI半监督精炼
- 协变量整合 (MT%, stress, cell cycle)
- Unknown细胞质量过滤

---

## 2️⃣ 细胞类型特异性分析 (Cell-type Specific Pipelines)

### 2.1 上皮细胞 (Epithelial)

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `epithelial_scvi_scanvi_20260118_v2_6.py` | v2.6 | 最新版本 |
| `epithelial_scvi_scanvi_NO_AT_20260118_v2_5.py` | v2.5 | 排除AT细胞版本 |
| `epithelial_celltypist_cnmf_pipeline_20260115_v1_4.py` | v1.4 | CellTypist + cNMF整合 |
| `complete_epithelial_analysis_pipeline_20260122_v2_1.py` | v2.1 | ✅ **完整流程** - BBKNN + Harmony + cNMF |

**预期亚型**: Basal, Suprabasal, Ciliated, Secretory (Goblet/Club), SMG, Ionocytes, Brush

### 2.2 T/NK细胞 (T/NK Cells)

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` | **v3.5.1** | ✅ **生产版本** - P0修复 |
| `tcell_cd4cd8_RETRAIN_by_sample_20260111_v1_0.py` | v1.0 | CD4/CD8按样本重训练 |
| `tcell_subcluster_analysis_20260126_v1_0.py` | v1.0 | 亚群分析 |
| `tnk_subcluster_unified_20260129_v2_0.py` | v2.0 | 统一亚群分析 |

**预期亚型**: CD4+ (naive/CM/EM/TRM), CD8+ (EM/EMRA/TRM), Tregs, NK, NKT, MAIT, ILC

### 2.3 B细胞 (B Cells)

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` | **v3.5.1** | ✅ **生产版本** - P0修复 |
| `bcell_subcluster_analysis_v2_20260119.py` | v2.0 | 亚群分析 |
| `bcell_analysis_notebook_20260114.py` | v1.0 | 交互式分析 |

**预期亚型**: Naive B, Memory B (IgM/IgG/IgA), Plasma cells, Plasmablasts, GC B cells, Bregs

### 2.4 髓系细胞 (Myeloid)

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` | **v3.5.1** | ✅ **生产版本** - P0修复 |
| `myeloid_subcluster_analysis_v2_20260121.py` | v2.0 | 亚群分析 |
| `myeloid_celltype_subcluster_pipeline_20260121_v4_1.py` | v4.1 | 细胞类型亚群流程 |
| `myeloid_subcluster_converter_20260128_v2.py` | v2.0 | 亚群标签转换 |

**预期亚型**: Macrophages (alveolar/interstitial/CCL+/CHIT1+), Monocytes (CD14+/CD16+), DCs (DC1/DC2/pDC), Mast cells

### 2.5 基质/血管细胞 (Stromal/Vascular)

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py` | **v3.5.1** | ✅ **生产版本** - P0修复 |
| `stromal_subcluster_analysis_v2_20260121.py` | v2.0 | 亚群分析 |
| `stromal_celltype_subcluster_pipeline_20260121_v4_1.py` | v4.1 | 细胞类型亚群流程 |
| `stromal_pure_scanvi_reanalysis_20260114.py` | v1.0 | 纯scANVI重分析 |
| `stromal_subcluster_converter_20260128_v2.py` | v2.0 | 亚群标签转换 |

**预期亚型**: Endothelial (arterial/venous/lymphatic/capillary), Fibroblasts (adventitial/alveolar/peribronchial), SMC, Pericytes, Mesothelial

---

## 3️⃣ 亚群聚类分析 (Subcluster Analysis)

**功能**: 在细胞类型内部进行精细聚类，保留BBKNN批次校正

| 脚本名称 | 细胞类型 | 版本 | 特性 |
|---------|---------|------|------|
| `epithelial_subcluster_bbknn_pipeline_20260122_v4_5_2.py` | Epithelial | **v4.5.2** | ✅ **最新** - HOTFIX |
| `bcell_subcluster_analysis_v2_20260119.py` | B cells | v2.0 | 生产版本 |
| `myeloid_subcluster_analysis_v2_20260121.py` | Myeloid | v2.0 | 生产版本 |
| `stromal_subcluster_analysis_v2_20260121.py` | Stromal | v2.0 | 生产版本 |
| `tcell_subcluster_analysis_20260126_v1_0.py` | T cells | v1.0 | 初始版本 |
| `tnk_subcluster_unified_20260129_v2_0.py` | T/NK | v2.0 | 统一流程 |

**核心技能**:
- 使用 `restrict_to` 参数保留BBKNN图
- 固定低分辨率 (0.2-0.3) 避免过度分裂
- 标记基因过滤 (去除MT/Ribo/IEG/技术伪影)
- 在原始UMAP上可视化

---

## 4️⃣ 通用亚群流程 (Universal Subcluster Pipelines)

**功能**: 可应用于任意细胞类型的通用亚群分析框架

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `universal_celltype_subcluster_pipeline_20260121_v4_1.py` | **v4.1** | ✅ **推荐** - 最新通用流程 |
| `universal_celltype_subcluster_pipeline_20260121_v3_2.py` | v3.2 | 稳定版本 |
| `universal_celltype_subcluster_pipeline_v3_0.py` | v3.0 | 基础版本 |

**核心技能**:
- 参数化细胞类型输入
- 自动化标记基因分析
- 批次效应检测
- 多分辨率聚类测试

---

## 5️⃣ cNMF基因表达程序分析 (cNMF - Gene Expression Programs)

**功能**: 识别基因表达程序(GEPs)并映射到细胞状态

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `label_guided_cnmf_pipeline_20260114_v1_1.py` | **v1.1** | ✅ **推荐** - 标签引导cNMF |
| `cnmf_results_analysis_20260114_v1_1.py` | v1.1 | cNMF结果分析 |
| `epithelial_celltypist_cnmf_pipeline_20260115_v1_4.py` | v1.4 | 上皮细胞专用 |
| `complete_epithelial_analysis_pipeline_20260122_v2_1.py` | v2.1 | BBKNN + Harmony + cNMF完整流程 |

**快速启动脚本** (按细胞类型):
- `quick_start_cnmf_analysis_20251225_B_v1_1.py` - B细胞
- `quick_start_cnmf_analysis_20251225_T_v1_1.py` - T细胞
- `quick_start_cnmf_analysis_20251225_Myeloid_v1_1.py` - 髓系
- `quick_start_cnmf_analysis_20251227_Stromal_v1_1.py` - 基质
- `quick_start_cnmf_analysis_ciliated_20251227_v1_1.py` - 纤毛细胞
- `quick_start_cnmf_analysis_club_20251227_v1_1.py` - Club细胞
- `quick_start_cnmf_analysis_AT2_20251227_v1_1.py` - AT2细胞
- `quick_start_cnmf_analysis_goblet_20260105_v1_1.py` - 杯状细胞

**核心技能**:
- 技术基因过滤 (MT/Ribo/Histone/IEG)
- K值选择 (10-50)
- 多进程加速 (4-16 workers)
- GEP-聚类映射

---

## 6️⃣ 质量控制与数据清洗 (Quality Control & Data Cleaning)

| 脚本名称 | 功能 | 版本 |
|---------|------|------|
| `celltypist_doublet_detection_20260114_v1_1.py` | **双细胞检测** | v1.1 ✅ |
| `remove_doublets_and_gse299751_20260114.py` | 移除双细胞和特定数据集 | v1.0 |
| `check_data_quality.py` | 数据质量检查 | - |
| `quick_data_check_detailed_20260127.py` | 详细数据检查 | v1.0 |
| `query_mapped_to_reference_qc_20260129.py` | 映射质量控制 | v1.0 |

**核心技能**:
- Top2谱系 + margin逻辑双细胞检测
- 多证据门控 (邻域纯度、映射置信度)
- 技术变异评估 (MT%, stress, cell cycle)
- 批次效应检测

---

## 7️⃣ scArches模型迁移 (scArches - Transfer Learning)

**功能**: 将新数据映射到已训练的参考模型

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `scarches_mapping_20260127_v1_2_1.py` | **v1.2.1** | ✅ **推荐** - 最新映射流程 |
| `step2_scarches_mapping_20260112.py` | v1.0 | 基础映射 |

**核心技能**:
- 参考模型加载
- 查询数据预处理
- 手术式微调 (surgical fine-tuning)
- 标签迁移

---

## 8️⃣ 数据转换与工具 (Data Conversion & Utilities)

| 脚本名称 | 功能 |
|---------|------|
| `subset_by_celltype_20260127_v1_0.py` | 按细胞类型子集化数据 |
| `myeloid_subcluster_converter_20260128_v2.py` | 髓系亚群标签转换 |
| `stromal_subcluster_converter_20260128_v2.py` | 基质亚群标签转换 |

---

## 9️⃣ scHPL / treeArches 层级校验 🆕

**功能**: 在 base model（scANVI/scArches）之后作为独立层级校验器，判断 query 是否被 reference hierarchy 稳定吸收

| 脚本名称 | 模式 | 版本 | 特性 |
|---------|------|------|------|
| `lineage_merge_schpl_core_20260402.py` | global | **v1.0** | ✅ **推荐** - global lineage-wise 复用核心 |
| `bcell_merge_schpl_20260412_v1.py` | global | v1.0 | B 细胞 scHPL 谱系合并 |
| `tcell_merge_schpl_20260402_v1.py` | global | v1.0 | T/NK scHPL 谱系合并 |
| `myeloid_merge_schpl_20260402_v1.py` | global | v1.0 | 髓系 scHPL 谱系合并 |
| `bcell_schpl_reject_followup_20260412_v1.py` | follow-up | v1.0 | B 细胞 Rejected 跟进 |
| `stromal_schpl_reject_followup_20260408_v1_0.py` | follow-up | v1.0 | 基质 Rejected 跟进 |
| `stromal_schpl_postprocess_v1_0_20260401.R` | downstream | v1.0 | 基质 scHPL R 侧后处理 |

**核心技能**:
- kNN classifier (dimred=True, useRE=True, n_neighbors=50)
- 标准输出列：`schpl_pred`、`schpl_prob`、`schpl_rejected`、`schpl_reject_type`
- Rejected 细胞跟进：marker DE、邻域分析、新类型候选识别
- 方法学详见 `SCHPL_METHODOLOGY_UNIFIED_20260410.md`

---

## 🔟 组织间比较分析 (Tissue Comparison) 🆕

**功能**: 各谱系亚型的组织间差异分析，包含 DESeq2 pseudobulk + GSVA + MASC 框架

### 通用框架

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `tissue_comparison_generic_wrapper_20260412.R` | **v1.0** | ✅ **推荐** - 通用封装器 |
| `tissue_comparison_advanced_helper_20260408.R` | v1.0 | 高级辅助函数 |
| `tissue_comparison_analysis_20251222.R` | v3.0 | 基础多组织比较 |

### B 细胞组织比较

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `bcell_tissue_comparison_v2_6_4_20260412.R` | **v2.6.4** | ✅ **最新** |
| `bcell_tissue_comparison_v2_6_3_20260411.R` | v2.6.3 | |
| `bcell_tissue_comparison_v2_6_2_20260410.R` | v2.6.2 | |
| `bcell_tissue_comparison_v2_6_1_20260410.R` | v2.6.1 | |
| `bcell_tissue_comparison_v2_6_20260410.R` | v2.6 | |
| `bcell_tissue_comparison_v2_6_20260410_llm_resume.R` | v2.6 | LLM 续跑版 |
| `bcell_tissue_comparison_v2_6_20260406.R` | v2.6 | 早期版 |
| `bcell_tissue_comparison_v2_5_20260402.R` | v2.5 | |
| `bcell_tissue_comparison_v2_4_20260390.R` | v2.4 | |

### T/NK 组织比较

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `tnk_tissue_comparison_v2_6_1_20260413.R` | **v2.6.1** | ✅ **最新** |
| `tnk_tissue_comparison_v2_6_20260407.R` | v2.6 | |
| `tnk_tissue_comparison_v2_4_20260401.R` | v2.4 | |

### 髓系组织比较

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `myeloid_tissue_comparison_20260407_v1_1.R` | **v1.1** | ✅ **最新** |
| `myeloid_tissue_comparison_20260401_v1_0.R` | v1.0 | |

### 上皮组织比较

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `epithelial_tissue_comparison_v1_3_2_20260413.R` | **v1.3.2** | ✅ **最新** |
| `epithelial_tissue_comparison_20260407_v1_3_1.R` | v1.3.1 | |
| `epithelial_tissue_comparison_20260407_v1_3.R` | v1.3 | |
| `epithelial_tissue_comparison_20260402_v1_2.R` | v1.2 | |
| `epithelial_tissue_comparison_20260331_v1_1.R` | v1.1 | |

### 基质/血管亚型组织比较

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `stromal_endothelial_tissue_comparison_20260408_v1_0.R` | v1.0 | 内皮细胞 |
| `stromal_fibroblast_tissue_comparison_20260408_v1_0.R` | v1.0 | 成纤维细胞 |
| `stromal_smc_tissue_comparison_20260408_v1_0.R` | v1.0 | 平滑肌细胞 |

**核心技能**:
- DESeq2 pseudobulk 差异表达
- GSVA 通路活性评分
- MASC 混合效应细胞组成分析
- LLM 辅助结果解读与续跑

---

## 1️⃣1️⃣ SCENIC 基因调控网络 🆕

**功能**: 基于 SCENIC 框架推断转录因子调控网络

| 脚本名称 | 语言 | 特性 |
|---------|------|------|
| `prepare_rscenic_cistarget_db_20260410.py` | Python | cisTarget 数据库下载与格式化 |
| `configure_scenic_r_env_20260410.R` | R | SCENIC R 环境配置 |
| `scenic_core_20260410.R` | R | ✅ **核心引擎** - GRNBoost2 + cisTarget + AUCell |
| `tnk_tcell_scenic_L3_20260410.R` | R | T/NK L3 精细级别 SCENIC |
| `bcell_scenic_L3_20260412.R` | R | B 细胞 L3 精细级别 SCENIC |

---

## 1️⃣2️⃣ 表型–基因关联管道 🆕

**功能**: 多模块管道，从 pseudobulk DE 到图神经网络链接预测

| 脚本名称 | 语言 | 模块 | 特性 |
|---------|------|------|------|
| `phenotype_gene_pipeline_A_pseudobulk_hpo_enrichment.R` | R | A | Pseudobulk + HPO 表型富集 |
| `phenotype_gene_pipeline_B_supervised_ranking.py` | Python | B | 监督学习基因排序 |
| `phenotype_gene_pipeline_C_graph_link_prediction.py` | Python | C | 图链接预测 |
| `phenotype_gene_pipeline_D_module_score_hpo.R` | R | D | 模块评分 + HPO |

---

## 1️⃣3️⃣ 差异丰度分析 (Milo/PAGA) 🆕

**功能**: 基于邻域的差异丰度检验

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `bcell_paga_milo_20260412_v1_1.py` | v1.1 | B 细胞 PAGA 轨迹 + Milo 差异丰度 |
| `myeloid_milopy_tissue_celltype_20260410_v1.py` | v1.0 | 髓系 Milo（按组织 × 细胞类型） |

---

## 1️⃣4️⃣ scANVI 验证与准备 🆕

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `bcell_scanvi_validation_20260412.py` | v1.0 | B 细胞 scANVI 注释验证 |
| `prepare_bcell_scanvi_umap_input_20260413.py` | v1.0 | B 细胞 scANVI UMAP 输入准备 |

---

## 1️⃣5️⃣ 基质重整合 (Stromal Reintegration) 🆕

| 脚本名称 | 版本 | 特性 |
|---------|------|------|
| `stromal_reintegration_branchwise_scvi_scanvi_20260407_v1_5.py` | **v1.5** | ✅ 分支式 scVI/scANVI 重整合 |
| `stromal_scarches_query_branchwise_20260407_v1_0.py` | v1.0 | 分支式 scArches 查询映射 |
| `epithelial_scarches_query_mapping_20260315_v1_1.py` | v1.1 | 上皮 scArches 查询映射 |

### 场景1: 新数据全流程分析
```bash
# Step 1: 全细胞整合
python allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py

# Step 2: 双细胞检测
python celltypist_doublet_detection_20260114_v1_1.py

# Step 3: 细胞类型特异性分析 (选择对应类型)
python t_scvi_celltypist_scanvi_20260108_v3_5_1.py
python bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py
python myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py
python stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py

# Step 4: 亚群分析
python universal_celltype_subcluster_pipeline_20260121_v4_1.py

# Step 5: cNMF分析
python label_guided_cnmf_pipeline_20260114_v1_1.py
```

### 场景2: 上皮细胞完整分析
```bash
# 一站式流程 (BBKNN + Harmony + cNMF)
python complete_epithelial_analysis_pipeline_20260122_v2_1.py
```

### 场景3: 模型迁移学习
```bash
# 将新数据映射到已有参考
python scarches_mapping_20260127_v1_2_1.py
```

---

## 📈 版本演进关键特性

### v2.x系列 (All-cells基线)
- ✅ 基因对齐 (`normalize_gene_names()`)
- ✅ UMAP分离 (`neighbors_key`)
- ✅ CPU安全CUDA
- ✅ 索引对齐写入
- ✅ 单GPU控制
- ✅ 稀有类型处理

### v2.3.1系列 (协变量 + Unknown清洗)
- ✅ 协变量系统 (MT%, stress, cell cycle)
- ✅ Unknown清洗 (邻域纯度 + 质量门控)
- ✅ 类别不平衡处理 (`n_samples_per_label`)
- ✅ scArches参数优化

### v3.5.x系列 (P0修复 - 细胞类型流程)
- ✅ P0-1: scANVI专用adata_model (标签对齐修复)
- ✅ P0-2: predict()方法修复 (dtype安全)
- ✅ P0-3: 多数投票整合
- ✅ P0-4: 特征验证 (CellTypist基因重叠检查)
- ✅ P0-5: 预训练HVG使用var_names
- ✅ P1-6: 内存优化 (60-70%减少)

### v4.x系列 (亚群分析)
- ✅ `restrict_to` 参数保留BBKNN图
- ✅ 固定低分辨率 (0.2-0.3)
- ✅ 技术基因过滤
- ✅ 批次效应监控

### scHPL/treeArches系列 (层级校验) 🆕
- ✅ scHPL 作为 base model 后的 secondary hierarchical validator
- ✅ global lineage-wise 与 branch-wise 双模式
- ✅ 标准 Rejected follow-up 框架
- ✅ 标准输出列：schpl_pred, schpl_prob, schpl_rejected, schpl_reject_type

### Tissue Comparison v2.6系列 (组织比较) 🆕
- ✅ 全谱系覆盖 (B/TNK/Myeloid/Epithelial/Stromal各亚型)
- ✅ 通用封装器 (tissue_comparison_generic_wrapper)
- ✅ LLM 续跑支持
- ✅ DESeq2 + GSVA + MASC 三维分析

### SCENIC系列 (基因调控) 🆕
- ✅ cisTarget 数据库自动准备
- ✅ GRNBoost2 + cisTarget + AUCell 完整管道
- ✅ L3 精细级别谱系特异性分析

---

## 🔧 技能总结 (Skills Summary)

### 核心分析技能
1. **批次校正** - scVI, BBKNN, Harmony
2. **自动注释** - CellTypist (Human_Lung_Atlas, Immune_All_Low)
3. **半监督精炼** - scANVI (单层/双层架构)
4. **亚群聚类** - Leiden (restrict_to保留批次校正)
5. **基因表达程序** - cNMF (K选择, GEP映射)
6. **双细胞检测** - Top2谱系 + margin逻辑
7. **模型迁移** - scArches (surgical fine-tuning)
8. **质量控制** - 协变量整合, Unknown过滤
9. **层级校验** 🆕 - scHPL/treeArches (global lineage-wise + branch-wise)
10. **组织间比较** 🆕 - DESeq2 pseudobulk + GSVA + MASC (全谱系覆盖)
11. **基因调控网络** 🆕 - SCENIC (GRNBoost2 + cisTarget + AUCell)
12. **差异丰度** 🆕 - Milo 邻域分析 + PAGA 轨迹
13. **表型–基因关联** 🆕 - HPO 富集 + 监督学习 + 图链接预测 + 模块评分

### 数据处理技能
9. **基因名标准化** - ENSEMBL → HGNC symbols
10. **HVG优化** - 内存减少60-70%
11. **稀有类型合并** - 防止训练不稳定
12. **索引对齐** - 防止重排序bug
13. **Counts恢复** - 从外部文件匹配恢复

### 可视化技能
14. **双UMAP管理** - scVI vs scANVI分离
15. **标记基因热图** - 过滤技术伪影
16. **批次效应检测** - UMAP叠加可视化
17. **GEP-聚类映射** - 相关性分析

---

## 📝 文件命名规范

**格式**: `{celltype}_{analysis}_{YYYYMMDD}_v{X}_{Y}.py`

- `celltype`: allcells, epithelial, tcell, bcell, myeloid, stromal
- `analysis`: scvi_celltypist_scanvi, subcluster, cnmf, bbknn
- `YYYYMMDD`: 日期戳
- `v{X}_{Y}`: 主版本.次版本

**示例**:
- `allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py` - 全细胞v2.3.1
- `epithelial_subcluster_bbknn_pipeline_20260122_v4_5_2.py` - 上皮亚群v4.5.2

---

## 🚀 快速查找指南

| 需求 | 推荐脚本 |
|------|---------|
| 全细胞整合 | `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py` |
| T细胞分析 | `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` |
| B细胞分析 | `b_scvi_celltypist_scanvi_pipeline_v1.py` |
| 髓系分析 | `myeloid_scvi_scanvi_v2_4_20260322.py` |
| 基质分析 | `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py` |
| 上皮完整流程 | `complete_epithelial_analysis_pipeline_20260122_v2_1.py` |
| 亚群聚类 | `universal_celltype_subcluster_pipeline_20260121_v4_1.py` |
| cNMF分析 | `label_guided_cnmf_pipeline_20260114_v1_1.py` |
| 双细胞检测 | `celltypist_doublet_detection_20260114_v1_1.py` |
| 模型迁移 | `scarches_mapping_20260127_v1_2_1.py` |
| scHPL 层级校验 🆕 | `lineage_merge_schpl_core_20260402.py` |
| B 细胞组织比较 🆕 | `bcell_tissue_comparison_v2_6_4_20260412.R` |
| T/NK 组织比较 🆕 | `tnk_tissue_comparison_v2_6_1_20260413.R` |
| 髓系组织比较 🆕 | `myeloid_tissue_comparison_20260407_v1_1.R` |
| 上皮组织比较 🆕 | `epithelial_tissue_comparison_v1_3_2_20260413.R` |
| 通用组织比较 🆕 | `tissue_comparison_generic_wrapper_20260412.R` |
| SCENIC GRN 🆕 | `scenic_core_20260410.R` |
| Milo 差异丰度 🆕 | `bcell_paga_milo_20260412_v1_1.py` |
| 表型–基因关联 🆕 | `phenotype_gene_pipeline_A_pseudobulk_hpo_enrichment.R` |

---

**生成时间**: 2026-04-13
**脚本总数**: 58 个 Python 脚本 + 139 个 R 脚本 + 16 个 Notebook = 218 个文件
**推荐版本**: 标记为 ✅ 的脚本为生产环境推荐版本
