# SCHPL / treeArches 统一方法学说明

更新时间：2026-04-10

本目录下的 `schpl` 相关文件目前分成两种实现形态，但应遵循同一套解释框架：

1. **global lineage-wise**：在单一谱系内训练一棵 scHPL 树，再对 query 做层级预测。代表文件：
   - `bcell_merge_schpl_20260330_v3_0.ipynb`
   - `lineage_merge_schpl_core_20260402.py`
   - `tcell_merge_schpl_20260402_v1.py`
   - `myeloid_merge_schpl_20260402_v1.py`
2. **branch-wise**：当 query 已在上游流程中按大类分支拆开，且不同大类共享一棵树不合适时，对每个 branch 单独训练/预测。代表文件：
   - `stromal_schpl_treeArches_20260330_v1_0.ipynb`
   - `stromal_schpl_treeArches_20260330_v1_0 copy.ipynb`
   - `stromal_schpl_reject_followup_20260408_v1_0.py`

## 统一原则

### 1. scHPL 的角色

- scHPL / treeArches 是 **base model（通常为 scANVI / scArches）之后的独立层级校验器**。
- 它的目标不是替代 base model，而是回答两个问题：
  - query 是否能被 reference hierarchy 稳定吸收；
  - 哪些细胞虽然有 transferred label，但在层级结构里仍应被保留为可疑/待跟进对象。

### 2. 输入空间

- scHPL 一律运行在 **上游整合后的 latent embedding** 上，而不是原始 counts。
- 具体键名可以因流水线而异（例如 `X_scANVI`、`X_scanvi`、`X_scANVI_L2`），但语义应保持一致：
  - 必须是 reference 与 query 共享的整合潜空间；
  - UMAP 只用于可视化，不参与 scHPL 训练。

### 3. 推荐参数

- classifier：`knn`
- `dimred=True`
- `useRE=True`
- `n_neighbors=50`
- `FN=0.5`
- `rej_threshold=0.5`

说明：kNN 更适合低维或中维 integrated latent；`Rejected` 应由距离、重建误差、概率三类标准共同定义。

### 4. 标准输出列

无论 global 还是 branch-wise，最终都应优先写出以下标准列：

- `schpl_pred_raw`：原始 scHPL 输出
- `schpl_pred`：清洗后的最终输出（接受标签或 `Rejected` / `Ignored`）
- `schpl_prob`：posterior probability
- `schpl_rejected`：是否被 reject
- `schpl_reject_type`：`dist` / `RE` / `prob` / `accepted` / `branch_skip`
- `leiden_schpl_qc`：基于 latent 的 QC cluster
- `schpl_novel_candidate`：高 reject cluster 中的候选新群体

仅 branch-wise 流程额外需要：

- `schpl_branch`
- `schpl_branch_source`

### 5. 对 Rejected 的统一解释

- `Rejected` **不是** “已经证明是新细胞类型”。
- `Rejected` 表示该细胞当前 **不能被 reference hierarchy 稳定吸收**，后续需要：
  1. 看是否集中于某个 query label / cluster；
  2. 比较 base model 标签与 scHPL 结果；
  3. 在同一 branch / label 内对 `Rejected` vs `Accepted` 做 marker / DE；
  4. 再判断是 novel state、污染、技术漂移还是层级未覆盖。

### 6. 何时用 global，何时用 branch-wise

- **优先 global lineage-wise**：当 reference / query 本身就在同一条清晰谱系内，标签空间可由一棵树表达时。
- **只在必要时 branch-wise**：当上游已经按大类拆分，或不同宏观谱系混在同一树中会导致训练语义不稳定时。
- 当前目录中：
  - B / TNK / Myeloid 属于 global lineage-wise；
  - Stromal 属于 branch-wise，因为其标签空间横跨 Endothelial / Fibroblast / Smooth Muscle / Pericyte 等大分支。

### 7. 文件职责分层

- `lineage_merge_schpl_core_20260402.py`：global lineage-wise 的复用核心。
- `*_merge_schpl_*.py`：不同谱系的配置入口。
- `stromal_schpl_treeArches_*.ipynb`：stromal 的 branch-wise 主流程。
- `stromal_schpl_reject_followup_20260408_v1_0.py`：reject 之后的标准 follow-up。
- `stromal_schpl_postprocess_v1_0_20260401.R`：旧版 stromal 输出的 R 侧下游整理；方法学上属于 downstream summary，不是主推断阶段。

## 实际落地建议

后续若继续新增 `schpl` 文件，建议默认沿用以下口径：

1. 开头先说明自己属于 **global lineage-wise** 还是 **branch-wise**。
2. 明确写出：scHPL 是 base model 之后的 **secondary hierarchical validator**。
3. 明确写出：`Rejected` 只是 follow-up flag，不是自动新类型结论。
4. 标准输出列尽量复用现有 `schpl_*` 命名，避免另起别名。
5. 若是 branch-wise 文件，必须说明 branch assignment 的来源和 `Ignored` 的定义。