# scrna_py_2025

面向 **鼻腔 / 呼吸道组织单细胞 RNA-seq** 的双语言分析仓库，核心场景包括免疫细胞分析、CRSwNP（Chronic Rhinosinusitis with Nasal Polyps）相关研究，以及参考图谱映射与下游统计解释。

仓库采用 **Python + R** 双轨工作流：

- `py/`：负责 scVI / scANVI / CellTypist / scArches 等深度学习整合、注释与可视化
- `R/`：负责 Seurat v5 质控、统计分析、差异表达、轨迹分析与结果解释

> 当前仓库以 `R/` 和 `py/` 目录中的脚本为主；根目录仅保留共享入口、桥接脚本和少量历史文件。

---

## 仓库结构

```text
scrna_py_2025/
├── R/                     # Seurat、QC、DE、trajectory、interpretation
├── py/                    # scVI/scANVI、CellTypist、BBKNN、scArches
├── GetSeurat.R            # h5ad → Seurat 读取桥接
├── README.md
└── CLAUDE.md
```

### `py/` 主要内容

Python 侧主线流程通常遵循三阶段：

1. **scVI**：基于原始 counts 进行批次校正与潜空间学习
2. **CellTypist**：自动细胞类型注释
3. **scANVI**：结合标签进行半监督精炼

常见任务包括：

- 全细胞整合与注释
- T / B / Myeloid / Stromal / Epithelial 等谱系专项分析
- BBKNN / marker 分析
- scArches 参考构建与 query 映射
- Notebook 可视化与结果核查

### `R/` 主要内容

R 侧重点在统计与解释：

- RDS / h5ad 导入与格式桥接
- QC、DoubletFinder、DecontX、智能合并
- Seurat v5 整合、亚群分析与 marker 解读
- MAST / DESeq2 / pseudobulk / MASC
- Monocle3 轨迹分析
- BayesPrism bulk RNA 解卷积
- LLM 辅助亚群解释脚本

---

## 数据对象约定

### Python / AnnData

- `adata.X`：处理后的表达矩阵
- `adata.layers['counts']`：原始 UMI counts（scVI 必需）
- `adata.raw`：保留全基因表达信息

常见 `obs` 字段：

- `celltype`
- `majority_voting`
- `Multinomial_Label`
- `sample` / `sample_id`
- `tissue`

### R / Seurat

推荐通过 `GetSeurat()` 导入 h5ad：

```r
source("GetSeurat.R")
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

### 1. Python：全细胞整合与注释

推荐从 `py/` 目录中的主流程脚本开始，例如：

- `py/allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py`
- `py/t_scvi_celltypist_scanvi_20260108_v3_5_1.py`
- `py/myeloid_scvi_scanvi_v2_4_20260322.py`
- `py/stromal_scvi_celltypist_scanvi_pipeline_20260110_v3_5_1.py`

### 2. R：QC + 智能合并

- `R/rds_folder_qc_smartmerge_20251216_v5.R`
- `R/gene_standardization_merge_20251221_v4_2.R`

### 3. Python → R 桥接

Python 侧导出：

```python
adata.write_h5ad("output.h5ad")
```

R 侧导入：

```r
source("GetSeurat.R")
seurat_obj <- GetSeurat("output.h5ad", validate_counts = TRUE)
```

### 4. R：下游统计与解释

常见入口：

- `R/MAST_pipeline_20251223.R`
- `R/pseudobulk.R`
- `R/tissue_comparison_analysis_20251222.R`
- `R/interpret.R`
- `R/bayesprism_deconvolution_production_20260124.R`

---

## 环境建议

### Python

建议使用独立 conda 环境（仓库内文档常以 `bbknn_env` 为例）：

```bash
conda create -n bbknn_env python=3.10
conda activate bbknn_env
pip install scanpy scvi-tools celltypist bbknn anndata
```

常见依赖：

- `scanpy`
- `scvi-tools`
- `celltypist`
- `bbknn`
- `anndata`
- `torch`
- `lightning`

### R

常用包包括：

- `Seurat`
- `harmony`
- `DoubletFinder`
- `celda` / `DecontX`
- `MAST`
- `DESeq2`
- `clusterProfiler`
- `monocle3`
- `reticulate`

---

## 文件组织说明

- **优先使用 `R/` 与 `py/` 下的脚本**，它们是当前的主目录
- 根目录保留少量共享脚本与兼容入口
- 脚本通常带日期与版本号，便于追踪迭代
- 若根目录与子目录存在同名且内容完全一致的文件，应以子目录版本为准

命名示例：

- `allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py`
- `myeloid_scvi_scanvi_v2_4_20260322.py`
- `bcell_interpret_PRODUCTION_v2_0_20260127.R`
- `epithelial_subcluster_interpret_20260210_v4.R`

---

## 相关文档

- `CLAUDE.md`
- `py/CLAUDE.md`
- `R/CLAUDE.md`
- `py/SCRIPT_CATEGORIZATION.md`
- `py/skills/README.md`

---

## 备注

当前仓库包含较多历史脚本和实验版本。若你是首次接手，建议优先从以下入口理解整体流程：

1. `GetSeurat.R`
2. `R/rds_folder_qc_smartmerge_20251216_v5.R`
3. `py/allcells_scvi_celltypist_scanvi_pipeline_20260121_v2_5_2.py`
4. `py/universal_celltype_subcluster_pipeline_20260121_v4_1.py`
5. `py/scarches_mapping_20260127_v1_2_1.py`

