# 单细胞项目经验总表

更新时间：2026-04-15  
适用范围：`/home/h2048` 工作区当前已验证的 scRNA-seq 分析经验、方法学约定、工程坑点与可复用结论。  
本文件目的：把分散在脚本头注释、repo memory、近期调试记录和产出目录中的经验，收敛成一份可直接接手使用的总说明。

> 建议把本文件视为“项目经验总表”；后续新增经验优先继续追加到这里，而不是再散落到多个零散笔记里。

## 0. 仓库整理记录（2026-06-04）

- 根目录历史重复脚本已清理，重复内容仅在 `R/` 与 `py/` 子目录保留。
- 后续引用脚本路径时，默认使用子目录完整路径（如 `R/GetSeurat.R`、`py/...`）。

## 1. 仓库整体判断

这是一个以呼吸道/鼻部组织单细胞分析为主的仓库，核心是双语言工作流：

- Python 侧负责整合、映射、latent 建模、自动注释、schpl/treeArches 层级校验
- R 侧负责 QC 后下游统计、tissue comparison、CHOIR/ssGSEA/DE 解释与 LLM 汇总

当前稳定思路不是“谁都重写一遍”，而是：

- **Python 统一潜空间与层级校验**
- **R 统一共享 tissue-comparison engine + 薄 wrapper**
- **所有 lineage 尽量共享方法学与输出命名**

## 2. 总体方法学主线

### 2.1 Python 主线

推荐理解为：

1. `scVI`：在原始 counts 上做整合与 batch effect removal
2. `CellTypist` / `scArches`：给出迁移或初始标签
3. `scANVI`：半监督 refine label
4. `scHPL / treeArches`：在共享 latent 上做 **secondary hierarchical validator**

关键原则：

- `scHPL` 不是替代 base model，而是用来回答“这个 query label 能否被 hierarchy 稳定吸收”。
- `scHPL` 一律运行在 **整合后的 latent** 上，不在 counts 或 UMAP 上训练。
- UMAP 只负责可视化，**不参与** `scHPL` 训练。

### 2.2 R 主线

当前更推荐：

- 用共享 engine 统一处理 pseudobulk DE、ssGSEA、CHOIR、LLM bundle
- 每个 lineage 只做薄 wrapper，覆写最少的参数
- 旧 wrapper/helper 尽量保持不可变；新增逻辑通过新 wrapper 或小版本 engine 扩展

## 3. schpl / treeArches 统一经验

### 3.1 统一口径

当前 schpl 相关实现分两类：

1. **global lineage-wise**：一条谱系一棵树，适用于 B/TNK/Myeloid
2. **branch-wise**：大分支拆开后分别训练/预测，适用于 Stromal

统一原则：

- `Rejected` 只是 **follow-up flag**，不是自动“新细胞类型”结论
- 标准输出列尽量统一为：
  - `schpl_pred_raw`
  - `schpl_pred`
  - `schpl_prob`
  - `schpl_rejected`
  - `schpl_reject_type`
  - `leiden_schpl_qc`
  - `schpl_novel_candidate`
- branch-wise 流程额外保留：
  - `schpl_branch`
  - `schpl_branch_source`

解释 reject 的标准顺序应是：

1. 看 reject 是否集中在某个 query label / cluster
2. 对比 base label 与 schpl 输出
3. 做 rejected vs accepted marker / DE
4. 再判断是 novel state、污染、技术漂移还是 hierarchy 覆盖不足

### 3.2 B 细胞 schpl 经验

B 细胞这条线近期已经形成比较明确的结论：

- downstream basis 应继续信任 **`X_scANVI_L2`**，不要因为局部标签不稳就把 latent 一起推翻。
- B 细胞 Naive / Memory / Plasma 的混淆，不一定是 schpl 新引入的问题；更大概率是上游 label 边界本身就存在不稳定区。
- 历史 `0208` 合并产物中的 `L2_scanvi_pred` 对 query 侧是**退化的**（几乎全 Plasma），**不能**当作正常 gold standard。
- 如果要做“旧 scanvi 结果验证”，应优先使用旧映射输出：
  - `Cell_Type_L2_pred`
  - `mapping_confidence`
- 当前最终结论是：schpl 重新分配的细胞多数是旧映射里置信度较低的边界细胞，而不是大面积新错误。

B 细胞 scanvi-UMAP 实战经验：

- 如果新 merged 输入里只有 `X_umap`，但历史参考里有 `X_umap_scanvi`，想强制“所有 UMAP 都用 scanvi 的”，**最稳妥做法不是改所有 plotting 脚本**，而是先把输入 h5ad 中的 `X_umap` 直接替换为 scanvi UMAP。
- 已验证脚本：
  - `script/py/prepare_bcell_scanvi_umap_input_20260413.py`
- 该脚本做法：
  - 通过 `barcode` 对齐当前 merged input 和历史 scanvi artifact
  - 写入 `obsm['X_umap_scanvi']`
  - 同时把 `obsm['X_umap']` 设成 scanvi 坐标
  - 旧 UMAP 另存为 `X_umap_legacy_pre_scanvi_swap`
- 这样后续所有默认吃 `X_umap` 的主流程和 follow-up 脚本都会自动走 scanvi UMAP，而无需逐个脚本打补丁。

B 细胞当前权威产物：

- scanvi UMAP 注入后的输入：
  - `/home/h2048/data/py/0413/bcell_merge_schpl_scanvi_umap_input_v1_0/reference_plus_query_merged_scanvi_umap_v1_0.h5ad`
- 主 schpl 结果：
  - `/home/h2048/data/py/0413/bcell_merge_schpl_v1_2_scanvi_umap_20260413`
- reject follow-up：
  - `/home/h2048/data/py/0413/bcell_schpl_reject_followup_v1_2_scanvi_umap_20260413`

该次 rerun 已验证：

- final output 的 `X_umap` 与注入的 scanvi UMAP **完全一致**（`maxabs = 0`, `meanabs = 0`）
- 主流程结果：
  - reference 11,570
  - query 41,217
  - accepted 40,181（97.5%）
  - rejected 1,036（2.5%）
  - 无 novel candidate cluster
- reject follow-up 主要候选热点仍是：
  - `Naive_B`
  - `Atypical_Memory_B`

### 3.3 Stromal branch-wise 经验

Stromal 是当前 branch-wise schpl 的标准范例。

关键经验：

- query 端必须按分支处理：`fibroblast`、`endothelial`、`smc`
- `smc` 分支实际对应 `Smooth_Muscle + Pericyte`
- `Schwann` 不在当前 modeled branches 里，应 **ignore**，不要硬塞进 branch-wise stromal schpl
- 训练 reference 时必须保留原始细 cluster ID（如 `source_cluster_id`），否则 DROP / REASSIGN 逻辑会悄悄错误引用 remapped 后的 `cell_type_L3`
- branch-wise stromal **不能**继续吃旧全局 query latent；必须使用与 branch reference 匹配的 branch-specific scArches/scANVI latent，否则 reject 会被 latent mismatch 人为放大
- `scvi-tools 1.3.3` 下，stromal branchwise：
  - `SCVI.train()` 可在 GPU 跑
  - `SCANVI.train()` 可能在 GPU 上直接因 `torchmetrics` 崩掉（`multiclass_accuracy` / `torch.AcceleratorError`）
  - 因此流程必须允许：先保存 SCVI 权重，再回退到 CPU 重跑 SCANVI

已验证结论：

- 2026-04-07 完整 rerun 后，branch-wise stromal 总 reject 降到 **23.8%**（82,296 modeled cells）
- endothelial 仍然偏高（42.6%），但 fibroblast / smc 明显低很多（16.1% / 7.9%）
- 说明在修掉“全局 latent 错配”后，剩余主要问题是 endothelial 分支内部问题，而不是跨分支路由错误

当前权威 stromal 结果：

- 主结果：`/home/h2048/data/py/0406/stromal_schpl_v1_1_branchwise`
- branch export：`branch_h5ad/` 和 `figures/branches/`
- reject follow-up：`/home/h2048/data/py/0408/stromal_schpl_reject_followup_v1_0`
- 旧 `v1.0` 输出：`/home/h2048/data/py/0330/stromal_schpl_v1_0` 已移除，不建议回头继续引用

## 4. Tissue comparison / LLM / CHOIR 经验

### 4.1 共享 engine 是当前主线

当前 tissue-comparison wrapper 的共享 engine 思路已经比较成熟：

- 共享 engine 基线：`script/R/bcell_tissue_comparison_v2_6_20260406.R`
- 对 myeloid 的 CHOIR tuning 扩展版：`script/R/bcell_tissue_comparison_v2_6_2_20260414.R`
- wrapper 的角色应尽量保持为：
  - 指定输入 / 输出
  - 指定 lineage preset
  - 做少量 override

### 4.2 LLM bundle 经验

当前更稳定的 LLM 模式：

- DESeq2 tissue-pair enrichment 和 CHOIR OFA enrichment 都已经通过 `prepare_llm_enrichment_bundle()` 统一成多数据库证据包
- 之前唯一仍保留“每数据库各跑一次 LLM”的地方是 grouped ssGSEA；现已修正为：
  - 先把 `top_df` 按 `SSGSEA_METHODS` 聚合
  - 再对每个 tissue × cell-type group 调一次 LLM
- 经验上，这样比逐数据库零碎调用稳定得多，也更容易比较证据一致性

### 4.3 CHOIR 参数经验

CHOIR 在当前仓库里很容易出现“过细切分”或重计算过重，因此需要保守调参。

共享层经验：

- 为减少当前 wrapper 的 CHOIR over-splitting，默认建议：`CHOIR_ALPHA <- 0.20`

Myeloid 特化经验：

- 当前薄入口脚本：
  - `script/R/myeloid_tissue_comparison_v1_2_2_20260414.R`
- 对应通用 wrapper 扩展：
  - `script/R/tissue_comparison_generic_wrapper_20260414_v3.R`
- 该 wrapper 通过覆写 shared engine，恢复并显式传递 CHOIR tuning，避免 `v1.2.1` 的 pathological heavy-default branch
- 当前 myeloid 保守 CHOIR 控制为：
  - `CHOIR_N_CORES = 1`
  - `CHOIR_SAMPLE_MAX = 5000`
  - `CHOIR_DOWNSAMPLING_RATE = 0.05`
  - `CHOIR_SUBTREE_REDUCTIONS = FALSE`
- 经验解释：在这里应优先追求 **稳定、可复现** 的 cluster state call，而不是超细粒度拆分

当前 myeloid 输出目录：

- `/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414`

### 4.4 TNK wrapper 经验

TNK 这条线有一个容易踩的语义坑：

- `cell_type_L2` 不能直接复用输入原始列
- 应从 `scanvi_label_refined` 重新推回 L2
- 否则像 `CD8 Teff`、`gdT` 之类容易在原 metadata 里同时落到 `CD8 T cells` 和 `NK cells`

当前已验证的 TNK L3 → L2 粗映射：

- 所有 `CD4*` → `CD4 T cells`
- `CD8*` + `gdT` + `MAIT` → `CD8 T cells`
- `NK` + `NK Exhausted` + `ILC3` → `NK cells`

TNK 当前验证通过的输出：

- `/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0`

成功导出包括：

- `tnk_tissue_comparison_final.rds`
- `tnk_tissue_comparison_final.h5ad`
- `interpret_agent_structured.tsv`
- `ssgsea_llm_structured.tsv`
- `choir_llm_structured.tsv`

## 5. h5ad / AnnData / reticulate / 稀疏矩阵经验

这是当前仓库最容易反复踩坑的一大块。

### 5.1 h5ad 写盘最常见问题

如果写 `.h5ad` 时在 `/obs` 或 `/raw/var` 报：

- `Can't implicitly convert non-string objects to strings`

优先怀疑：

- 混合 `object` 列
- `bool + NA`
- `str + int` 混合 ID 列
- factor / list 派生列含有全 NA

已验证的处理经验：

- `bool + NaN` 优先转成 `int8`
- `str + int` 这类 ID 列先统一转字符串
- 对 character / factor / list-derived 列中的 NA，先替换为空字符串再写盘
- tissue-comparison shared engine 中需要有 `sanitize_obs_for_h5ad()` 这一层，特别是 branch CHOIR 可能留下全 NA 的 `subcluster_*` 列

### 5.2 R → Python 稀疏导出

R 的 `Matrix::CsparseMatrix` 直接喂给 `scipy.sparse.csr_matrix` 容易报：

- `index pointer size ... should be ...`

原因是 R 的 `dgCMatrix` 本质是 CSC，不是 CSR。

已验证做法：

- 直接用 `scipy.sparse.csc_matrix`
- 或者先转 `RsparseMatrix` 再导成 `csr_matrix`
- triplet 导出时可用：
  - `as.data.frame(Matrix::summary(Matrix::t(as(mat, "dgCMatrix"))))`
- 并把 1-based 的 `i/j` 改成 0-based

### 5.3 reticulate 经验

- `import("numpy"|"scipy.sparse"|"anndata", convert = FALSE)` 更稳；否则对象可能被自动转回 R 的 S4，导致 `$tocsr()`、`layers`、`obsm` 写入异常
- 给 `AnnData$layers / uns / obsm` 赋值时，优先用 ``$`__setitem__`()`` 显式赋值
- 校验 `obs` 列名时，`py_to_r(adata$obs)` 后再 `colnames()` 更稳；不要过分依赖 `obs$columns$tolist()`

### 5.4 可视化小坑

在当前环境下：

- `scanpy 1.11.5` 的 `sc.pl.umap()` / `sc.pl.embedding()` 会在内部处理 `rasterized`
- 如果再显式传 `rasterized = TRUE`，可能报：
  - `Axes.scatter() got multiple values for keyword argument 'rasterized'`

因此：

- 当前环境里不要重复显式传 `rasterized = TRUE`

## 6. scvi-tools / scanvi / 环境经验

### 6.1 scvi-tools API 经验

在 `scvi-tools 1.x`：

- `model.train()` 不再推荐 `use_gpu`
- 应写成：
  - `accelerator = 'gpu' | 'cpu'`
  - `devices = ...`

如果在 notebook 里有多个 Step 分别训练 `scVI` / `scANVI`，要一起统一改，不要只改其中一个。

### 6.2 `SCANVI.from_scvi_model()` 经验

`scvi-tools 1.3.3` 下，如果在子集数据上做：

- `SCANVI.from_scvi_model()`

容易因为 batch category 数量变少而触发 size mismatch。

稳妥做法：

- 先把 `adata.obs[batch_key]` 的 `categories` 对齐到原 `SCVI` 模型 registry
- 再迁移权重

另外，即使做 fully supervised，并且已经过滤掉 `Unknown` 细胞：

- `unlabeled_category` 这个类别名也可能仍然需要保留在 `adata.obs[labels_key].cat.categories` 里
- 即使没有对应细胞，也建议保留一个未使用 category，避免 API 兼容问题

## 7. R / Python 跨环境调试经验

### 7.1 conda / ABI / 动态库污染

以下问题在当前机器上都是真实踩过的坑：

- `conda install scvi-tools` 在混合环境里可能通过 PyPI bridge 塞入大量 `pypi/pypi` 包
- 如果同时保留 conda 的 `numpy`，又装了旧版 pip `pandas`，容易报：
  - `numpy.dtype size changed`
- 更稳策略：
  - conda 固定底层二进制栈
  - pip 只用 `--no-deps` 装 `scvi-tools` / `scArches`

如果 `conda create --clone` 复制混合环境：

- 容易因大量 `pypi_0` 记录报 `PackagesNotFoundError`
- 更稳的备份方法是：**直接复制整个 env 目录** 做沙箱，再用新前缀下的 `python` 自检关键导入

如果运行时错误地加载了 base 环境里的 `libtorch_python.so`，报：

- `undefined symbol: PyObject_ClearManagedDict`

优先用：

- `env -u LD_LIBRARY_PATH -u PYTHONPATH conda run -p <env> ...`

清掉外层动态库污染。

### 7.2 base 环境不可信时的经验

如果 `/home/h2048/miniconda3` base 环境里 `scanpy` 导入失败，并且 `PIL/` 目录残缺：

- 不要优先在 base 上硬修
- 当前经验是：直接切到完整的 `scarches_stable` / `scvi_env` 更稳

### 7.3 R 安装源码包时的 curl 污染

Linux 上如果 R 装源码包时 `configure` 调 `curl` 报：

- `undefined symbol: curl_easy_header`
- 或 `curl/libcurl versions do not match`

大概率是 conda 的 `libcurl` 污染了系统 `curl`。

可临时用：

- `PATH=/usr/bin:/bin:$PATH R CMD INSTALL ...`

绕过 conda 的 curl。

### 7.4 非致命告警不要过度惊慌

例如在 `bbknn_env` 里跑 `GSVA` + `BiocParallel::SnowParam(type='SOCK')` 时，worker 启动可能出现：

- `/bin/bash: ... libtinfo.so.6: no version information available`

当前环境下这通常只是**非致命告警**，只要 worker 真能起来，ssGSEA / CHOIR 仍可跑完。

## 8. 其他已验证的 R 侧经验

### 8.1 CHOIR 调用方式

- `CHOIR::CHOIR` 如果手动提供 `reduction`，传入的必须是：
  - `Embeddings(obj, reduction = ...)` 得到的**矩阵**
- 不能直接传 reduction 名字符串
- 同时必须显式提供 `var_features`
- 否则会分别报：
  - “reduction 不是 matrix”
  - 或 `var_features cannot be NULL`

### 8.2 clusterProfiler / STRING 超时

`clusterProfiler::interpret_agent()` 如果开 `add_ppi = TRUE`，当前网络环境里容易卡在 STRINGdb：

- `stringdb-static.org/species.v12.txt`

经验上：

- 第一次超时后就应自动关闭后续 `add_ppi`
- 否则整条分析线每次会白白多拖约 60 秒

### 8.3 `sitecustomize.py` 的隐蔽坑

如果使用 `sitecustomize.py` 通过：

- `os.execve(sys.executable, [sys.executable, *sys.argv], env)`

修补环境，务必特判：

- `sys.argv[:1] == ['-c']`

否则会把 `python -c ...` 错改成裸 `python -c`，从而误伤 reticulate / 探针脚本。

## 9. 当前推荐的“接手顺序”

如果后面要继续接手这个仓库，推荐按下面顺序判断：

1. **先看本文件**，确定当前主线和已踩过的坑
2. 再看 `script/CLAUDE.md`，确认仓库全貌与目录
3. 如果是 schpl 相关，再看：
   - `script/py/SCHPL_METHODOLOGY_UNIFIED_20260410.md`
4. 如果是 tissue comparison 相关，优先确认：
   - 当前是共享 engine 哪个版本
   - wrapper 只覆写了哪些参数
5. 如果要做 B 细胞 UMAP 相关判断，默认使用 `/0413/` 的 scanvi-UMAP rerun 结果
6. 如果要做 stromal branch-wise 判断，默认使用 `/0406/` + `/0408/` 的保留结果
7. 如果遇到 h5ad/reticulate/export 问题，先检查本文件第 5 节，而不是先怀疑数据本身坏掉
8. 如果遇到 scvi / torch / scanpy 导入问题，先检查本文件第 6-7 节，再决定是否重建环境

## 10. 当前最值得保留的权威路径

### Python

- B 细胞 schpl（scanvi UMAP 版）：
  - `/home/h2048/data/py/0413/bcell_merge_schpl_v1_2_scanvi_umap_20260413`
- B 细胞 reject follow-up：
  - `/home/h2048/data/py/0413/bcell_schpl_reject_followup_v1_2_scanvi_umap_20260413`
- Stromal branch-wise schpl：
  - `/home/h2048/data/py/0406/stromal_schpl_v1_1_branchwise`
- Stromal reject follow-up：
  - `/home/h2048/data/py/0408/stromal_schpl_reject_followup_v1_0`

### R

- Myeloid tissue comparison（CHOIR tuned）：
  - `/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414`
- TNK tissue comparison：
  - `/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0`

## 11. 一句话结论

当前仓库最重要的经验不是“再造一套新流程”，而是：

- **沿用统一方法学**
- **保住共享 engine / 共享 latent / 共享命名**
- **把 branch 特例和环境坑点显式写出来**
- **优先复用已验证输出，而不是反复从旧版本半成品重启**

如果后续继续补经验，建议直接在本文件追加新小节，保持它作为唯一总表。