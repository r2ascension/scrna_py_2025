---
name: latest-processed-files
description: "Use when you need to find the latest processed versions in /home/h2048/script/py, especially 202602-202603 merged training, ref-query merge, scArches mapping, reintegration, and visualization entry files. Trigger phrases: 最新处理版本, 文件索引, 哪个脚本最新, merge pipeline, scArches mapping, ref query."
license: Complete terms in LICENSE.txt
metadata:
  title: Latest Processed Files Index
  version: "2026-03-17"
  author: GitHub Copilot
  category: documentation
  difficulty: beginner
  tags:
    - index
    - latest-files
    - scvi
    - scanvi
    - scarches
    - merge
    - visualization
  outputs:
    - Recommended latest entry files by workflow
    - Mapping between older and newer script families
    - Notes about filename-version mismatches
  related_skills:
    - allcells-integration
    - celltype-specific-analysis
    - scarches-mapping
    - subcluster-analysis
---

# Skill: Latest Processed Files Index

## 功能描述

这个索引专门补充 [SCRIPT_CATEGORIZATION.md](/home/h2048/script/py/SCRIPT_CATEGORIZATION.md) 和 [skills/README.md](/home/h2048/script/py/skills/README.md) 尚未覆盖的 2026-02 到 2026-03 新流程。重点不是替代 2026-01 的生产基线，而是帮助快速定位后续新增的 merge、reference/query、filtered retrain、reintegration 和 visualization 文件。

## 适用场景

- 想知道某条分析链目前应从哪个最新脚本开始。
- 需要区分 reference model、query mapping、merge、merged retrain 这几个阶段。
- 需要快速找到 202602-202603 新增而旧分类文档未收录的处理版本。
- 需要先做文件检索，再决定是否深入阅读 notebook 或 pipeline 脚本。

## 使用原则

1. 2026-01 的基础生产流程仍以 [CLAUDE.md](/home/h2048/script/py/CLAUDE.md) 和 [SCRIPT_CATEGORIZATION.md](/home/h2048/script/py/SCRIPT_CATEGORIZATION.md) 为主。
2. 2026-02 到 2026-03 的新增流程优先看本索引。
3. 文件名中的版本号和文件头注释不总是一致，优先同时核对日期、用途、输入输出路径。

## 最新处理版本索引

### 1. B 细胞 L2 reference/query 链路

| 阶段 | 推荐文件 | 旧文件或相邻文件 | 用途 |
|---|---|---|---|
| Reference training | `step1_train_bcell_L2_20260204_v2_5_3.py` | 2026-01 的 B cell 主流程 | 训练 B 细胞 L2 reference，保存 UMAP operator、reference_with_umap 和 HVG 列表。 |
| Query mapping | `step2_map_query_20260204_v2_5_4.py` | `step2_map_query_20260204_v2_5_3.py` | 将 query 映射到 L2 reference，并修正 metadata restore 与 merge 稳定性问题。 |
| Full-gene merge + visualization | `20260204_complete_merge_and_visualize_v1_2.py` | `20260204_standalone_merge_ref_query.py`, `20260204_complete_merge_and_visualize_v1_1.py` | 在 full-gene 保留前提下完成 ref+query 合并，适合后续 merged training 之前的结果检查。 |
| Merged retrain | `merged_scanvi_training_20260208_v2_5_5.py` | `bcell_merged_scanvi_training_20260207_v2_5_3.py` | 2 月这条线的生产热修版本，统一 label 构建、强化 counts 校验，并补充诊断输出。 |

### 2. 上皮细胞 3 月组合流程

| 阶段 | 推荐文件 | 旧文件或相邻文件 | 用途 |
|---|---|---|---|
| L3 可视化入口 | `epithelial_subcluster_visualization_L3_20260205_v1_0.py` | 2026-01 的 subcluster pipeline | 把 L3 映射、marker 面板和出版级图整合到一个可视化入口。 |
| SELF + REF 联合训练 | `epithelial_scvi_scanvi_20260308_v2_7.py` | `epithelial_scvi_scanvi_20260114_v2_4.py`, `epithelial_scvi_scanvi_NO_AT_20260118_v2_6.py` | 3 月上皮主线，单次 scVI 训练后分 SELF 和 REF 两个分支，给后续 scArches query mapping 提供 reference 基础。 |
| Query mapping | `epithelial_scarches_query_mapping_20260315_v1_1.py` | `scarches_mapping_20260127_v1_2_1.py` | 基于 v2.6 combined reference 的 query 映射脚本，分别处理 scVI、major scANVI、fine scANVI，并保存完整 posterior 概率矩阵。 |

补充说明：
- [epithelial_subcluster_visualization_L3_20260205_v1_0.py](/home/h2048/script/py/epithelial_subcluster_visualization_L3_20260205_v1_0.py) 需要注意 L3 映射时先清理 subcluster 和映射键中的隐藏空白。
- [epithelial_scvi_scanvi_20260308_v2_7.py](/home/h2048/script/py/epithelial_scvi_scanvi_20260308_v2_7.py) 的文件名写 v2_7，但文件头标的是 v2.6-COMBINED，应以用途和日期为主。

### 3. Myeloid merged 和 filtered retrain 链路

| 阶段 | 推荐文件 | 旧文件或相邻文件 | 用途 |
|---|---|---|---|
| Ref-query merge | `myeloid_only_ref_query_merge_pipeline_20260225_v1_1.py` | `myeloid_only_ref_query_merge_pipeline_20260225_v1_0.py` | 先把 myeloid reference 和 query 合并到同一训练对象，修复了 v1.0 误指向 T 细胞路径的问题。 |
| Filtered retrain | `myeloid_only_merged_pipeline_20260315_v2_2.py` | `myeloid_only_merged_pipeline_20260315_v2_1.py` | 基于 CellTypist 过滤污染细胞后重训 scVI/scANVI，是 3 月 myeloid 主线的后处理版本。 |

补充说明：
- [myeloid_only_merged_pipeline_20260315_v2_2.py](/home/h2048/script/py/myeloid_only_merged_pipeline_20260315_v2_2.py) 的文件头仍写 v2.1.1 HOTFIX，路径名则已经到 v2_2，检索时要按文件路径和更新日期识别。

### 4. T/NK reference、merge 和 filtered retrain 链路

| 阶段 | 推荐文件 | 旧文件或相邻文件 | 用途 |
|---|---|---|---|
| scArches-ready reference | `tnk_scvi_scanvi_scarches_ref_20260315_v1_1.py` | 2026-01 的 `t_scvi_celltypist_scanvi_20260108_v3_5_1.py` | 构建 T/NK reference，补充全局 finite 检查、sanitize_sparse_counts 和 4-tier HVG fallback。 |
| Ref-query merge | `tcell_only_ref_query_merge_pipeline_20260225_v1_1.py` | `tcell_only_ref_query_merge_pipeline_20260225_v1_0.py` | T cell only 的 ref-query merge 主脚本；文件头已升级到 v1.2，但路径名仍停留在 v1_1。 |
| Filtered retrain | `tcell_only_merged_pipeline_20260315_v2_2.py` | 2026-02 的 merged 预处理链路 | 基于 CellTypist 去污染后重训 T/NK-only 模型，并修复 cycling NK、raw 附着时机和标签审计逻辑。 |

### 5. Stromal reintegration 和 query mapping

| 阶段 | 推荐文件 | 旧文件或相邻文件 | 用途 |
|---|---|---|---|
| Reintegration notebook 参考 | `stromal_reintegration_20260303_v1_1.ipynb` | 2026-01 stromal subcluster 系列 | 3 月初基质 reintegration 的关键 notebook，后续 query mapping 依赖它产出的 reference 标签体系。 |
| Query mapping | `stromal_scarches_query_20260315_v2_1.py` | 2026-01 的 `scarches_mapping_20260127_v1_2_1.py` | 基于 reintegration v1.3 reference 执行 stromal query mapping，修复了 SCANVI label space 和 labels_key 不匹配问题。 |

补充说明：
- `stromal_reintegration_20260303_v1_1.ipynb` 中真实的 L3 细分需要由 `cell_type_L2 + '_c' + subcluster_id` 重建，不能直接依赖 `cell_type_L3`。
- tissue 列要做规范化和关键词兜底，不能只靠精确匹配肺组织值。

### 6. 通用 merged visualization 和结果检查

| 阶段 | 推荐文件 | 旧文件或相邻文件 | 用途 |
|---|---|---|---|
| Merge quality check | `20260204_complete_merge_and_visualize_v1_2.py` | `20260204_complete_merge_and_visualize.py` | 适合检查 ref/query 合并是否保持 full-gene universe、HVG 交集和 query-only 颜色列。 |
| Query result review | `visualize_query_results_20260204.py` | 同目录 2026-01 可视化脚本 | 用于 query 映射结果的快速复核。 |
| Merged result review | `visualize_scanvi_results_20260208.py` | `20260205_visualize_merged_highquality.py` | 用于 merged scanVI 结果检查和图形输出。 |

## 当前索引缺口已经补到哪里

本 skill 已补入以下旧文档未系统覆盖的主题：

- 2 月 B cell L2 reference/query/merge/merged retrain 全链路。
- 3 月 epithelial combined SELF+REF 和独立 query mapping。
- 3 月 myeloid 与 T/NK 的 filtered retrain 主线。
- 3 月 stromal reintegration 后的 scArches query mapping。

## 检索建议

如果只是想尽快找到入口文件，优先按下面顺序搜：

1. 文件名里包含 `merged_pipeline` 或 `merged_scanvi_training`，通常表示 ref+query 合并后的重训练阶段。
2. 文件名里包含 `scarches_query`、`map_query` 或 `ref_query_merge`，通常表示 reference/query 映射或并置阶段。
3. 文件名里包含 `visualize`，通常表示结果复核或出图入口，不是训练主流程。
4. 文件名里包含 `ref` 或 `reference`，通常表示可复用 reference model 构建入口。

## 命名陷阱

- 文件路径版本号可能落后于文件头版本号，例如 `tcell_only_ref_query_merge_pipeline_20260225_v1_1.py` 的文件头实际是 v1.2。
- 文件路径版本号也可能领先于文件头热修描述，例如 `myeloid_only_merged_pipeline_20260315_v2_2.py` 的说明仍沿用 v2.1.1 HOTFIX 文案。
- 遇到这种情况时，先看输入输出路径是否指向同一条数据链，再决定使用哪个文件。