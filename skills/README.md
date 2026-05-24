---
name: single-cell-skills-index
description: Comprehensive index of single-cell RNA-seq analysis skills covering integration, cell type annotation, clustering, quality control, and advanced analysis methods
version: 1.0
author: Single-cell Analysis Pipeline Team
tags:
  - index
  - documentation
  - single-cell
  - workflow
  - pipeline
category: documentation
---

# Single-Cell RNA-seq Analysis Skills Index

本目录包含单细胞RNA-seq分析的核心技能文档，每个skill对应一个主要分析流程。

## 📚 Skills列表

### 1. [All-Cells Integration](./allcells-integration.md)
**功能**: 全细胞整合分析 - scVI + CellTypist + scANVI三阶段流程

**适用场景**:
- 新数据集的初始整合
- 多批次数据批次校正
- 自动化细胞类型注释
- 大规模数据集 (>100k cells)

**推荐脚本**: `allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py` (v2.3.1)

**元数据**:
- Category: `integration`
- Difficulty: `intermediate`
- Estimated Time: `2-3 hours`

---

### 2. [Cell-Type Specific Analysis](./celltype-specific-analysis.md)
**功能**: 细胞类型特异性深度分析 - 识别细胞亚型和功能状态

**适用场景**:
- 已完成全细胞整合，需要细胞类型深度分析
- 细胞类型数量 >10k cells
- 免疫细胞精细分型 (T/B/Myeloid)
- 基质/血管细胞异质性分析

**支持的细胞类型**:
- T/NK cells: `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` (v3.5.1 ✅)
- B cells: `bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` (v3.5.1 ✅)
- Myeloid: `myeloid_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py` (v3.5.1 ✅)
- Stromal: `stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py` (v3.5.1 ✅)
- Epithelial: `epithelial_scvi_scanvi_20260118_v2_6.py` (v2.6)

**元数据**:
- Category: `analysis`
- Difficulty: `intermediate`
- Estimated Time: `2-4 hours`

---

### 3. [Subcluster Analysis](./subcluster-analysis.md)
**功能**: 细胞类型内精细聚类 - 保留BBKNN批次校正

**适用场景**:
- 已完成细胞类型注释，需要识别功能亚群
- 细胞类型内部异质性高
- 需要识别疾病相关细胞状态
- 需要保留批次校正效果

**推荐脚本**:
- 通用: `universal_celltype_subcluster_pipeline_20260121_v4_1.py` (v4.1 ✅)
- 上皮: `epithelial_subcluster_bbknn_pipeline_20260122_v4_5_2.py` (v4.5.2 ✅)
- B细胞: `bcell_subcluster_analysis_v2_20260119.py` (v2.0)
- 髓系: `myeloid_subcluster_analysis_v2_20260121.py` (v2.0)
- 基质: `stromal_subcluster_analysis_v2_20260121.py` (v2.0)

**元数据**:
- Category: `analysis`
- Difficulty: `intermediate`
- Estimated Time: `30-60 minutes`

---

### 4. [cNMF Gene Expression Program Analysis](./cnmf-analysis.md)
**功能**: 识别基因表达程序(GEPs) - 揭示细胞状态和功能模块

**适用场景**:
- 识别细胞状态和功能模块
- 发现转录调控程序
- 比较不同批次校正方法
- 疾病相关基因表达模式
- 发育轨迹和细胞分化

**推荐脚本**:
- 标签引导: `label_guided_cnmf_pipeline_20260114_v1_1.py` (v1.1 ✅)
- 结果分析: `cnmf_results_analysis_20260114_v1_1.py` (v1.1)
- 完整流程: `complete_epithelial_analysis_pipeline_20260122_v2_1.py` (v2.1)

**快速启动** (按细胞类型):
- B cells: `quick_start_cnmf_analysis_20251225_B_v1_1.py`
- T cells: `quick_start_cnmf_analysis_20251225_T_v1_1.py`
- Myeloid: `quick_start_cnmf_analysis_20251225_Myeloid_v1_1.py`
- Stromal: `quick_start_cnmf_analysis_20251227_Stromal_v1_1.py`

**元数据**:
- Category: `analysis`
- Difficulty: `advanced`
- Estimated Time: `2-6 hours`

---

### 5. [Doublet Detection and Quality Control](./doublet-detection.md)
**功能**: 双细胞检测和质量控制 - 基于CellTypist多谱系预测

**适用场景**:
- 高细胞密度区域
- 解离后数据集 (双细胞率5-10%)
- 混合谱系信号的细胞
- 质量控制流程的一部分

**推荐脚本**:
- 双细胞检测: `celltypist_doublet_detection_20260114_v1_1.py` (v1.1 ✅)
- 数据清洗: `remove_doublets_and_gse299751_20260114.py` (v1.0)
- 质量检查: `check_data_quality.py`
- 详细检查: `quick_data_check_detailed_20260127.py` (v1.0)

**元数据**:
- Category: `quality-control`
- Difficulty: `beginner`
- Estimated Time: `30-60 minutes`

---

### 6. [scArches Model Transfer Learning](./scarches-mapping.md)
**功能**: 模型迁移学习 - 将新数据映射到已训练参考模型

**适用场景**:
- 有高质量参考数据集和模型
- 新数据集需要快速注释
- 新数据与参考相似 (同组织/物种)
- 计算资源有限 (比从头训练快10-100x)
- 需要保持与参考的一致性

**推荐脚本**:
- scArches映射: `scarches_mapping_20260127_v1_2_1.py` (v1.2.1 ✅)
- 映射质控: `query_mapped_to_reference_qc_20260129.py` (v1.0)
- 基础映射: `step2_scarches_mapping_20260112.py` (v1.0)

**元数据**:
- Category: `transfer-learning`
- Difficulty: `advanced`
- Estimated Time: `30-60 minutes`

---

## 🔄 推荐工作流程

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

## 📊 技能依赖关系

```
[原始数据]
    ↓
[Doublet Detection] ← 质量控制
    ↓
[All-Cells Integration] ← 全细胞整合
    ↓
[Cell-Type Specific] ← 细胞类型深度分析
    ↓
    ├─→ [Subcluster Analysis] ← 亚群聚类
    └─→ [cNMF Analysis] ← 基因表达程序

[scArches Mapping] ← 替代路径 (如果有参考模型)
```

---

## 🎯 快速查找指南

| 需求 | 推荐Skill | Category | Difficulty |
|------|----------|----------|-----------|
| 新数据集初始分析 | All-Cells Integration | integration | intermediate |
| T细胞亚型鉴定 | Cell-Type Specific (T/NK) | analysis | intermediate |
| B细胞亚型鉴定 | Cell-Type Specific (B cells) | analysis | intermediate |
| 髓系细胞分型 | Cell-Type Specific (Myeloid) | analysis | intermediate |
| 基质/血管细胞分析 | Cell-Type Specific (Stromal) | analysis | intermediate |
| 上皮细胞完整流程 | Cell-Type Specific (Epithelial) + cNMF | analysis | intermediate |
| 细胞类型内亚群 | Subcluster Analysis | analysis | intermediate |
| 基因表达程序 | cNMF Analysis | analysis | advanced |
| 双细胞检测 | Doublet Detection | quality-control | beginner |
| 质量控制 | Doublet Detection | quality-control | beginner |
| 快速注释新数据 | scArches Mapping | transfer-learning | advanced |
| 模型迁移学习 | scArches Mapping | transfer-learning | advanced |

---

## 📋 Skills元数据总览

| Skill Name | Version | Category | Difficulty | Time |
|-----------|---------|----------|-----------|------|
| allcells-integration | v2.3.1 | integration | intermediate | 2-3h |
| celltype-specific-analysis | v3.5.1 | analysis | intermediate | 2-4h |
| subcluster-analysis | v4.1 | analysis | intermediate | 30-60min |
| cnmf-analysis | v1.1 | analysis | advanced | 2-6h |
| doublet-detection | v1.1 | quality-control | beginner | 30-60min |
| scarches-mapping | v1.2.1 | transfer-learning | advanced | 30-60min |

---

## 🔧 核心技能总结

### 分析技能 (8项)
1. **批次校正** - scVI, BBKNN, Harmony
2. **自动注释** - CellTypist (Human_Lung_Atlas, Immune_All_Low)
3. **半监督精炼** - scANVI (单层/双层架构)
4. **亚群聚类** - Leiden (restrict_to保留批次校正)
5. **基因表达程序** - cNMF (K选择, GEP映射)
6. **双细胞检测** - Top2谱系 + margin逻辑
7. **模型迁移** - scArches (surgical fine-tuning)
8. **质量控制** - 协变量整合, Unknown过滤

### 数据处理技能 (5项)
9. **基因名标准化** - ENSEMBL → HGNC symbols
10. **HVG优化** - 内存减少60-70%
11. **稀有类型合并** - 防止训练不稳定
12. **索引对齐** - 防止重排序bug
13. **Counts恢复** - 从外部文件匹配恢复

### 可视化技能 (4项)
14. **双UMAP管理** - scVI vs scANVI分离
15. **标记基因热图** - 过滤技术伪影
16. **批次效应检测** - UMAP叠加可视化
17. **GEP-聚类映射** - 相关性分析

---

## 📝 版本演进

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

---

## 🚀 使用建议

### 新手入门
1. 从 **Doublet Detection** 开始进行质量控制
2. 学习 **All-Cells Integration** 进行全细胞整合
3. 根据研究兴趣选择 **Cell-Type Specific** 分析

### 进阶分析
4. 使用 **Subcluster Analysis** 识别功能亚群
5. 使用 **cNMF Analysis** 发现基因表达程序
6. 学习 **scArches Mapping** 进行模型迁移

### 最佳实践
- ✅ 始终从质量控制开始
- ✅ 保存中间结果和模型
- ✅ 记录所有参数和版本
- ✅ 验证批次校正效果
- ✅ 使用多种方法交叉验证

---

## 📚 参考文献

### 核心方法
- **scVI**: Lopez et al., Nature Methods 2018
- **scANVI**: Xu et al., Molecular Systems Biology 2021
- **CellTypist**: Domínguez Conde et al., Science 2022
- **cNMF**: Kotliar et al., eLife 2019
- **scArches**: Lotfollahi et al., Nature Biotechnology 2022
- **BBKNN**: Polański et al., Bioinformatics 2020

### 数据集
- **Human Lung Atlas**: Sikkema et al., Nature Medicine 2023

---

**生成时间**: 2026-02-03
**文档版本**: v1.0
**维护者**: Single-cell Analysis Pipeline Team
