# scrna_py_2025

面向**鼻腔 / 鼻窦 / 呼吸道组织单细胞 RNA-seq**的分析仓库，覆盖从 QC、整合、注释、亚群分析，到参考映射、组织比较、cNMF 程序分析与结果可视化的完整工作流。

仓库以 **Python + R 双栈协作**为核心：

- `py/`：scVI / scANVI / CellTypist / scArches / cNMF / 可视化
- `R/`：Seurat v5 / QC / pseudobulk / MASC / Monocle3 / 富集分析
- 其余顶层目录：按专题拆分的模块化流程、配置与测试

> 当前主入口脚本统一保留在子目录中，根目录以说明文档和模块目录为主。

---

## 这是什么仓库

这个仓库主要服务于以下几类任务：

1. **全细胞整合与注释**：scVI → CellTypist → scANVI
2. **谱系专项分析**：T / B / Myeloid / Stromal / Epithelial
3. **Python ↔ R 桥接**：AnnData (`.h5ad`) 与 Seurat (`.rds`) 互转
4. **下游统计**：pseudobulk、MASC、差异表达、富集分析
5. **高级专题**：scArches 映射、non-unified airway 比较、normal airway ML、细胞通讯、SCENIC / cNMF

---

## 仓库结构

```text
scrna_py_2025/
├── README.md
├── R/                      # R/Seurat 主流程与统计分析
├── py/                     # Python/scanpy/scvi-tools 主流程与可视化
├── tests/                  # 主要 Python 单元测试
├── bcell/                  # B 细胞专题模块
├── epithelial/             # 上皮细胞专题模块与测试
├── non_unified_airway/     # 上下气道/非统一气道比较流程
├── normal_airway_ml/       # normal airway 机器学习与可视化模块
├── config/                 # YAML / JSON 运行配置
├── core/                   # 通用核心逻辑
└── 其他专题目录/文档
```

### 目录职责

#### `py/`

偏向深度学习整合与 AnnData 工作流，常见内容包括：

- scVI / scANVI 训练与重训练
- CellTypist 自动注释与标签精炼
- scArches 参考构建、query 映射、合并可视化
- cNMF、MILO、marker 可视化、结果审查
- 针对主要谱系的专项流程脚本

#### `R/`

偏向统计分析与生物学解释，常见内容包括：

- RDS / h5ad 导入与格式桥接
- QC、DoubletFinder、DecontX、智能合并
- Seurat v5 亚群分析与 marker 解释
- pseudobulk、MAST、MASC、轨迹分析
- BayesPrism、SCENIC、富集分析、报告型可视化

#### 模块化目录

- `bcell/`：B 细胞整合、内部 signature、可视化子模块
- `epithelial/`：上皮细胞细分流程与局部测试
- `non_unified_airway/`：气道组织分层、site group、通讯分析
- `normal_airway_ml/`：normal airway 机器学习训练/验证/可视化
- `config/`：为较新的流程提供参数文件，避免把配置硬编码进脚本

---

## 推荐入口

如果你第一次接手这个仓库，建议按下面顺序理解：

### 1. R 侧数据导入与 QC

- `R/GetSeurat.R`
- `R/rds_folder_qc_smartmerge_20251216_v5.R`
- `R/gene_standardization_merge_20251221_v4_2.R`

### 2. Python 侧全细胞整合

- `py/allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py`
- `py/universal_celltype_subcluster_pipeline_20260121_v4_1.py`

### 3. 谱系专项流程

- T / NK：`py/t_scvi_celltypist_scanvi_20260108_v3_5_1.py`
- B：`py/bcell_scvi_celltypist_scanvi_pipeline_20260111_v3_5_1.py`
- Myeloid：`py/myeloid_scvi_scanvi_v2_4_20260322.py`
- Stromal：`py/stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py`
- Epithelial：`py/epithelial_scvi_scanvi_20260317_v2_7.py`

### 4. 参考映射 / 合并 / 可视化

- `py/scarches_mapping_20260127_v1_2_1.py`
- `py/20260204_standalone_merge_ref_query.py`
- `py/20260204_complete_merge_and_visualize_v1_2.py`
- `py/visualize_query_results_20260204.py`

### 5. 下游统计与解释

- `R/pseudobulk.R`
- `R/MAST_pipeline_20251223.R`
- `R/tissue_comparison_analysis_20251222.R`
- `R/interpret.R`
- `R/bayesprism_deconvolution_production_20260124.R`

---

## 数据对象约定

### Python / AnnData

- `adata.X`：当前用于分析的表达矩阵
- `adata.layers["counts"]`：原始 UMI counts（scVI / scANVI 必需）
- `adata.raw`：保留全基因表达或原始矩阵

常见 `obs` 字段：

- `sample` / `sample_id`
- `dataset`
- `tissue`
- `celltype` 或相关标签列
- `majority_voting`
- `Multinomial_Label`

### R / Seurat

推荐通过 `GetSeurat()` 导入 h5ad：

```r
source("R/GetSeurat.R")

seurat_obj <- GetSeurat(
  h5ad_path = "path/to/file.h5ad",
  prefer_raw = TRUE,
  prefer_layer_counts = TRUE,
  validate_counts = TRUE,
  debug = TRUE
)
```

---

## 常见工作流

### 工作流 1：从 Python 整合到 R 统计

1. 在 `py/` 运行整合 / 注释主流程  
2. 导出 `.h5ad`  
3. 在 `R/` 使用 `GetSeurat()` 读入  
4. 继续做 pseudobulk、富集、轨迹或 tissue comparison

Python 导出：

```python
adata.write_h5ad("output.h5ad")
```

R 导入：

```r
source("R/GetSeurat.R")
seurat_obj <- GetSeurat("output.h5ad", validate_counts = TRUE)
```

### 工作流 2：新数据直接做全细胞注释

推荐从以下脚本起步：

- `py/allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py`
- 或按谱系拆分到对应专项流程

### 工作流 3：已有参考模型做 query 映射

推荐从以下脚本起步：

- `py/scarches_mapping_20260127_v1_2_1.py`
- `py/20260204_standalone_merge_ref_query.py`
- `py/20260205_visualize_merged_highquality.py`

### 工作流 4：R 侧组织比较与解释

推荐入口：

- `R/pseudobulk.R`
- `R/scMASC.R`
- `R/tissue_comparison_analysis_20251222.R`
- `R/enrichment_functions.R`

---

## 环境建议

### Python

仓库内脚本通常假设使用单独的 conda 环境（文档中常见名称为 `bbknn_env`）：

```bash
conda create -n bbknn_env python=3.10
conda activate bbknn_env
pip install scanpy scvi-tools celltypist bbknn anndata pandas numpy scipy h5py
```

常见依赖：

- `scanpy`
- `scvi-tools`
- `celltypist`
- `bbknn`
- `anndata`
- `pandas`
- `h5py`
- `torch`
- `lightning`

### R

常见包包括：

- `Seurat`
- `reticulate`
- `DoubletFinder`
- `celda` / `DecontX`
- `DESeq2`
- `MAST`
- `clusterProfiler`
- `monocle3`
- `BayesPrism`

---

## 测试

仓库当前包含以 Python 为主的回归/单元测试，主要位于：

- `tests/`
- `epithelial/tests/`

在当前沙箱环境中，直接运行：

```bash
python3 -m unittest discover -s tests -q
```

会因为缺少 `anndata`、`pandas`、`h5py` 等科学计算依赖而失败；在完整分析环境安装后再执行更合适。

---

## 文件组织约定

- 优先使用子目录中的脚本，不在根目录堆放重复分析脚本
- `R/` 与 `py/` 仍然是最主要的历史脚本入口
- 较新的专题流程逐步拆到独立模块目录中
- 文件名常包含日期与版本号，便于追踪迭代
- 若同类脚本存在多个版本，应优先查看日期更近、带明确版本号或配置文件配套的实现

---

## 相关文档

- `py/CLAUDE.md`
- `R/CLAUDE.md`
- `py/SCRIPT_CATEGORIZATION.md`
- `py/skills/README.md`
- `PROJECT_EXPERIENCE_CONSOLIDATED_20260415.md`

---

## 一句话建议

如果你只是想尽快找到主流程：

- **全细胞整合**：`py/allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py`
- **Python → R 桥接**：`R/GetSeurat.R`
- **QC 合并**：`R/rds_folder_qc_smartmerge_20251216_v5.R`
- **query 映射**：`py/scarches_mapping_20260127_v1_2_1.py`
- **组织比较**：`R/pseudobulk.R`
