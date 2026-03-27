#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
表型→基因关联 路线B：基因多视角特征 + 监督学习 → 基因排序
Phenotype→Gene Pipeline B: Gene Multi-view Features + Supervised Learning

工作流摘要 (Workflow Summary):
  1. 从 AnnData 对象按 (sample × cell_type) 聚合 → pseudobulk 计数矩阵
  2. 按细胞类型计算每个基因的多视角特征：
       avg_expr     : 各细胞类型平均表达（log1p 归一化后）
       pct_expr     : 各细胞类型中表达该基因的细胞比例
       logFC        : case vs control 的 log2 折叠变化
       module_corr  : 与 HPO 基因集的共表达相关性
  3. 用 HPO 基因集构建正/负例标签（"表型相关基因" vs 背景）
  4. 训练随机森林（或梯度提升）分类器，输出每个基因的"表型关联概率"
  5. 按概率降序输出 Top-K 候选基因排序

输入 (Inputs):
  adata_path      : AnnData h5ad 文件路径（含 raw counts、细胞类型注释、
                    样本 ID 与分组信息）
  hpo_ids         : 目标 HPO term ID 列表（用于构建正例标签）
  hpo_gene_path   : HPO phenotype_to_genes.txt 本地路径

输出 (Outputs):
  <output_dir>/gene_features.csv    : 基因特征矩阵
  <output_dir>/gene_labels.csv      : 基因标签（正/负例）
  <output_dir>/gene_ranking_B.csv   : 监督学习基因排序结果

依赖包 (Required Packages):
  anndata, scanpy, pandas, numpy, scipy, scikit-learn, joblib

参考 (References):
  - scanpy.get.aggregate():
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.aggregate.html
  - muscat pseudobulk best practice:
    https://doi.org/10.1038/s41467-020-19894-4
  - HPO genes_to_phenotype.txt:
    https://hpo.jax.org/app/download/annotation

版本 (Version): v1.0  2026-03-27
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ==============================================================================
# 0. 配置参数 Configuration
# ==============================================================================

# --- 真实路径替换区 (Replace paths for real analysis) ---
ADATA_PATH      = None   # 替换为真实 h5ad 路径，如 "/data/T_object.h5ad"
HPO_GENE_PATH   = None   # HPO phenotype_to_genes.txt 本地路径（None=下载）
OUTPUT_DIR      = "./results_pipeline_B"

# --- 分析参数 ---
CELL_TYPE_KEY   = "Annotation_2"   # adata.obs 中的细胞类型列
SAMPLE_KEY      = "sample"         # adata.obs 中的样本 ID 列
CONDITION_KEY   = "disease"        # adata.obs 中的条件列
CONDITION_CASE  = "CRSwNP"         # 实验组值
CONDITION_CTRL  = "Control"        # 对照组值

# HPO terms（定义目标表型；None=使用 Toy 演示）
TARGET_HPO_IDS  = None  # 如 ["HP:0000118", "HP:0001250", "HP:0002110"]

# 模型参数
N_ESTIMATORS    = 200   # 随机森林树数量
N_TOP_GENES     = 200   # 返回 top-k 基因数
RANDOM_STATE    = 42
TEST_SIZE       = 0.2   # 训练/测试集拆分比例
MIN_CELLS_PER_PSEUDOBULK = 3  # pseudobulk 最小细胞数

# ==============================================================================
# 1. 工具函数 Utilities
# ==============================================================================

def _check_import(module: str, install_hint: str) -> bool:
    """检查模块是否可用"""
    import importlib
    if importlib.util.find_spec(module) is None:
        log.warning("模块 %s 未安装: %s", module, install_hint)
        return False
    return True


def _safe_import_scanpy():
    import importlib
    spec = importlib.util.find_spec("scanpy")
    if spec is None:
        raise ImportError("请安装 scanpy: pip install scanpy")
    import scanpy as sc
    return sc


# ==============================================================================
# 2. HPO 基因集加载 Load HPO Gene Sets
# ==============================================================================

def load_hpo_gene_sets(local_path: Optional[str] = None) -> pd.DataFrame:
    """
    加载 HPO phenotype_to_genes.txt

    Returns
    -------
    pd.DataFrame with columns: hpo_id, gene_symbol
    """
    if local_path and Path(local_path).exists():
        log.info("从本地加载 HPO 基因集: %s", local_path)
        df = pd.read_csv(local_path, sep="\t", comment="#",
                         header=None, low_memory=False)
        # 尝试自动检测列格式
        # 新格式: hpo_id, hpo_name, entrez_id, gene_symbol, ...
        if df.shape[1] >= 4:
            df = df.iloc[:, [0, 3]]
            df.columns = ["hpo_id", "gene_symbol"]
        else:
            df = df.iloc[:, :2]
            df.columns = ["hpo_id", "gene_symbol"]
    else:
        log.info("尝试下载 HPO phenotype_to_genes.txt...")
        url = ("https://purl.obolibrary.org/obo/hp/hpoa/"
               "phenotype_to_genes.txt")
        try:
            df = pd.read_csv(url, sep="\t", comment="#",
                             header=None, low_memory=False)
            df = df.iloc[:, [0, 3]]
            df.columns = ["hpo_id", "gene_symbol"]
        except Exception as e:
            log.warning("HPO 下载失败: %s\n  使用 Toy HPO", e)
            df = _make_toy_hpo()

    df["gene_symbol"] = df["gene_symbol"].str.upper()
    df = df.dropna(subset=["hpo_id", "gene_symbol"])
    log.info("HPO 基因集: %d 行 (terms × genes)", len(df))
    return df


def get_hpo_gene_sets(hpo_df: pd.DataFrame) -> dict[str, set[str]]:
    """将 HPO df 转为 {hpo_id: {gene_symbol, ...}} 字典"""
    return {
        hpo_id: set(group["gene_symbol"].tolist())
        for hpo_id, group in hpo_df.groupby("hpo_id")
    }


# ==============================================================================
# 3. AnnData 加载与 Pseudobulk 聚合
# ==============================================================================

def load_adata(path: str):
    """加载 AnnData h5ad 文件"""
    sc = _safe_import_scanpy()
    log.info("加载 AnnData: %s", path)
    adata = sc.read_h5ad(path)
    log.info("  细胞: %d  基因: %d", adata.n_obs, adata.n_vars)
    return adata


def aggregate_pseudobulk(
    adata,
    cell_type_key: str = "Annotation_2",
    sample_key: str = "sample",
    condition_key: str = "disease",
    min_cells: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    按 (cell_type, sample) 聚合 raw counts → pseudobulk 矩阵

    Returns
    -------
    pb_counts : pd.DataFrame, shape (n_genes, n_pseudobulk_samples)
    pb_meta   : pd.DataFrame, columns [sample, cell_type, condition, n_cells]
    """
    sc = _safe_import_scanpy()
    import scipy.sparse as sp

    log.info("[Step 3] 聚合 pseudobulk 计数矩阵...")

    # 确保使用 raw counts
    if "counts" in adata.layers:
        X = adata.layers["counts"]
    elif adata.raw is not None:
        X = adata.raw.X
        # raw 可能是子集，需对齐到 adata.var_names
        X = X[:, [list(adata.raw.var_names).index(g)
                  for g in adata.var_names
                  if g in adata.raw.var_names]]
    else:
        X = adata.X
        log.warning("  未找到 counts layer / raw，使用 adata.X（可能已归一化）")

    if sp.issparse(X):
        X = X.toarray()

    obs = adata.obs[[cell_type_key, sample_key, condition_key]].copy()
    obs["_idx"] = np.arange(len(obs))

    pb_list = []
    pb_meta_list = []

    for ct, ct_df in obs.groupby(cell_type_key):
        for samp, samp_df in ct_df.groupby(sample_key):
            if len(samp_df) < min_cells:
                continue
            cell_idx = samp_df["_idx"].values
            pb_sum = X[cell_idx, :].sum(axis=0)
            cond = samp_df[condition_key].iloc[0]
            pb_list.append(pb_sum)
            pb_meta_list.append({
                sample_key:      samp,
                "cell_type":     ct,
                condition_key:   cond,
                "n_cells":       len(samp_df),
            })

    if not pb_list:
        raise ValueError("无任何 pseudobulk 样本通过最小细胞数过滤")

    pb_counts = pd.DataFrame(
        np.vstack(pb_list),
        columns=adata.var_names,
        index=[f"{r['cell_type']}__{r[sample_key]}" for r in pb_meta_list],
    ).T  # shape: (n_genes, n_pseudobulk_samples)

    pb_meta = pd.DataFrame(pb_meta_list,
                           index=pb_counts.columns)
    log.info("  Pseudobulk 矩阵: %d 基因 × %d 样本",
             pb_counts.shape[0], pb_counts.shape[1])
    return pb_counts, pb_meta


# ==============================================================================
# 4. 基因多视角特征工程 Gene Feature Engineering
# ==============================================================================

def compute_gene_features(
    adata,
    pb_counts: pd.DataFrame,
    pb_meta: pd.DataFrame,
    cell_type_key: str = "Annotation_2",
    condition_key: str = "disease",
    condition_case: str = "CRSwNP",
    condition_ctrl: str = "Control",
) -> pd.DataFrame:
    """
    按细胞类型计算每个基因的多视角特征

    特征（Feature）说明：
      avg_expr_{celltype}   : 该细胞类型 log1p 归一化后的平均表达
      pct_expr_{celltype}   : 该细胞类型中基因表达 > 0 的细胞比例
      logFC_{celltype}      : case vs control pseudobulk log2 FC
                              （log2(mean_case + 1) - log2(mean_ctrl + 1)）

    Returns
    -------
    pd.DataFrame, shape (n_genes, n_features)，index = gene_symbol
    """
    import scipy.sparse as sp
    log.info("[Step 4] 计算基因多视角特征...")

    cell_types = adata.obs[cell_type_key].unique().tolist()
    genes = list(adata.var_names)
    feat_dict: dict[str, np.ndarray] = {}

    # 获取 log1p 归一化表达（如已归一化则直接用 adata.X）
    if "counts" in adata.layers:
        import scanpy as sc
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        X_log = adata_norm.X
    else:
        X_log = adata.X

    if sp.issparse(X_log):
        X_log = X_log.toarray()

    for ct in cell_types:
        ct_mask = (adata.obs[cell_type_key] == ct).values
        if ct_mask.sum() < 5:
            continue

        ct_X = X_log[ct_mask, :]  # shape: (n_cells_ct, n_genes)
        feat_dict[f"avg_expr_{ct}"]  = ct_X.mean(axis=0)
        feat_dict[f"pct_expr_{ct}"]  = (ct_X > 0).mean(axis=0)

        # logFC（pseudobulk 层面）
        ct_pb = pb_counts.loc[:, pb_meta["cell_type"] == ct]
        ct_meta_ct = pb_meta[pb_meta["cell_type"] == ct]

        case_cols = ct_meta_ct.index[ct_meta_ct[condition_key] == condition_case]
        ctrl_cols = ct_meta_ct.index[ct_meta_ct[condition_key] == condition_ctrl]

        if len(case_cols) > 0 and len(ctrl_cols) > 0:
            mean_case = ct_pb[case_cols].mean(axis=1).values
            mean_ctrl = ct_pb[ctrl_cols].mean(axis=1).values
            logfc = np.log2(mean_case + 1) - np.log2(mean_ctrl + 1)
        else:
            logfc = np.zeros(len(genes))
        feat_dict[f"logFC_{ct}"] = logfc

    feature_df = pd.DataFrame(feat_dict, index=genes)
    feature_df = feature_df.fillna(0.0)
    log.info("  基因特征矩阵: %d 基因 × %d 特征",
             feature_df.shape[0], feature_df.shape[1])
    return feature_df


# ==============================================================================
# 5. 构建监督学习标签 Build Supervised Labels
# ==============================================================================

def build_gene_labels(
    genes: list[str],
    hpo_gene_sets: dict[str, set[str]],
    target_hpo_ids: list[str],
    negative_ratio: float = 5.0,
    seed: int = RANDOM_STATE,
) -> pd.Series:
    """
    基于目标 HPO terms 的基因集构建正/负例标签

    策略：
      正例 (+1) = 在任意一个目标 HPO term 的基因集中出现的基因
      负例 (0)  = 不在任意目标 HPO term 基因集中，且不与正例邻近

    Parameters
    ----------
    genes         : 所有基因列表
    hpo_gene_sets : get_hpo_gene_sets() 的返回值
    target_hpo_ids: 目标 HPO term ID 列表
    negative_ratio: 负例与正例的比例（默认 5:1）

    Returns
    -------
    pd.Series, index=gene, values={0, 1}（仅包含选中的正/负例）
    """
    log.info("[Step 5] 构建监督学习标签...")

    # 正例基因集（目标 HPO terms 的基因并集）
    positive_genes: set[str] = set()
    for hpo_id in target_hpo_ids:
        positive_genes.update(hpo_gene_sets.get(hpo_id, set()))

    # 只保留在特征矩阵中存在的基因
    genes_upper = [g.upper() for g in genes]
    gene_map = {g.upper(): g for g in genes}

    pos_in_data = [gene_map[g] for g in positive_genes
                   if g in gene_map]
    neg_candidates = [g for g in genes if g.upper() not in positive_genes]

    n_pos = len(pos_in_data)
    n_neg = min(len(neg_candidates), int(n_pos * negative_ratio))

    rng = np.random.default_rng(seed)
    neg_selected = rng.choice(neg_candidates, size=n_neg, replace=False).tolist()

    labels = pd.Series(
        [1] * n_pos + [0] * n_neg,
        index=pos_in_data + neg_selected,
        name="label",
    )
    log.info("  正例: %d  负例: %d  (比例 1:%.1f)",
             n_pos, n_neg, n_neg / max(n_pos, 1))
    return labels


# ==============================================================================
# 6. 训练分类器 Train Classifier
# ==============================================================================

def train_classifier(
    feature_df: pd.DataFrame,
    labels: pd.Series,
    n_estimators: int = N_ESTIMATORS,
    test_size: float = TEST_SIZE,
    seed: int = RANDOM_STATE,
) -> tuple:
    """
    训练随机森林分类器，返回 (model, eval_metrics, feature_importance)

    Parameters
    ----------
    feature_df : pd.DataFrame, shape (n_genes, n_features)
    labels     : pd.Series，只包含标注基因（正/负例）
    test_size  : 验证集比例

    Returns
    -------
    (model, metrics_dict, feature_importance_df)
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score, average_precision_score

    log.info("[Step 6] 训练随机森林分类器...")

    # 对齐特征矩阵与标签
    common_genes = feature_df.index.intersection(labels.index)
    X = feature_df.loc[common_genes].values
    y = labels.loc[common_genes].values

    if len(np.unique(y)) < 2:
        raise ValueError("标签只有一类，无法训练（检查 HPO 基因集是否与数据基因重叠）")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=seed
    )

    clf = RandomForestClassifier(
        n_estimators=n_estimators,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )
    clf.fit(X_train, y_train)

    # 评估
    y_score = clf.predict_proba(X_test)[:, 1]
    auc  = roc_auc_score(y_test, y_score)
    aupr = average_precision_score(y_test, y_score)
    log.info("  验证集 AUROC: %.4f  AUPR: %.4f", auc, aupr)

    # 特征重要性
    feat_imp = pd.DataFrame({
        "feature":   feature_df.columns.tolist(),
        "importance": clf.feature_importances_,
    }).sort_values("importance", ascending=False)

    metrics = {"AUROC": auc, "AUPR": aupr,
               "n_train": len(y_train), "n_test": len(y_test)}
    return clf, metrics, feat_imp


# ==============================================================================
# 7. 全基因组预测与排序 Genome-wide Prediction & Ranking
# ==============================================================================

def predict_and_rank(
    clf,
    feature_df: pd.DataFrame,
    top_k: int = N_TOP_GENES,
) -> pd.DataFrame:
    """
    对所有基因预测"表型关联概率"并排序

    Returns
    -------
    pd.DataFrame with columns: gene, phenotype_prob, rank
    """
    log.info("[Step 7] 全基因组预测与排序...")

    probs = clf.predict_proba(feature_df.values)[:, 1]
    result = pd.DataFrame({
        "gene":           feature_df.index,
        "phenotype_prob": probs,
    })
    result = result.sort_values("phenotype_prob", ascending=False).reset_index(drop=True)
    result["rank"] = result.index + 1
    result = result.head(top_k)
    log.info("  Top-%d 候选基因排序完成", len(result))
    return result


# ==============================================================================
# 8. Toy 数据生成 Toy Data Helpers
# ==============================================================================

def _make_toy_adata():
    """生成模拟 AnnData 用于演示"""
    import anndata as ad
    import scipy.sparse as sp

    np.random.seed(42)
    n_cells, n_genes = 600, 500
    X = np.random.negative_binomial(2, 0.5, size=(n_cells, n_genes)).astype(np.float32)

    # 注入差异表达：前 30 个基因在 case 细胞中高表达
    case_mask = np.array([i < n_cells // 2 for i in range(n_cells)])
    X[case_mask, :30] += 8

    cell_types = (["CD4_T"] * 200 + ["CD8_T"] * 200 + ["Macrophage"] * 200)
    samples    = ["S1"] * 100 + ["S2"] * 100 + ["S3"] * 100 + ["S4"] * 100 + \
                 ["S5"] * 100 + ["S6"] * 100
    conditions = (["CRSwNP"] * 300 + ["Control"] * 300)

    obs = pd.DataFrame({
        "Annotation_2": cell_types,
        "sample":       samples,
        "disease":      conditions,
    })
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)])

    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    log.info("Toy AnnData: %d 细胞 × %d 基因", adata.n_obs, adata.n_vars)
    return adata


def _make_toy_hpo() -> pd.DataFrame:
    """生成模拟 HPO 基因集"""
    np.random.seed(1)
    n_terms, n_genes = 50, 500
    rows = []
    for i in range(n_terms):
        hpo_id = f"HP:{i:07d}"
        gene_idx = np.random.choice(n_genes, 20, replace=False)
        for gi in gene_idx:
            rows.append({"hpo_id": hpo_id, "gene_symbol": f"GENE{gi}"})
    return pd.DataFrame(rows)


# ==============================================================================
# 9. 主流程 Main Pipeline
# ==============================================================================

def run_pipeline_B(
    adata_path: Optional[str] = None,
    hpo_gene_path: Optional[str] = None,
    target_hpo_ids: Optional[list[str]] = None,
    cell_type_key: str = CELL_TYPE_KEY,
    sample_key: str = SAMPLE_KEY,
    condition_key: str = CONDITION_KEY,
    condition_case: str = CONDITION_CASE,
    condition_ctrl: str = CONDITION_CTRL,
    output_dir: str = OUTPUT_DIR,
    top_k: int = N_TOP_GENES,
) -> dict:
    """
    运行路线 B 完整流程

    Returns
    -------
    dict with keys: feature_df, labels, model, metrics, gene_ranking
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- 真实数据 or Toy Demo ---
    if adata_path is None:
        log.info("=== [Toy Demo] 未提供 AnnData 路径，使用模拟数据 ===")
        adata = _make_toy_adata()
        cell_type_key   = "Annotation_2"
        sample_key      = "sample"
        condition_key   = "disease"
        condition_case  = "CRSwNP"
        condition_ctrl  = "Control"
    else:
        adata = load_adata(adata_path)

    # Step 3: Pseudobulk 聚合
    pb_counts, pb_meta = aggregate_pseudobulk(
        adata,
        cell_type_key = cell_type_key,
        sample_key    = sample_key,
        condition_key = condition_key,
        min_cells     = MIN_CELLS_PER_PSEUDOBULK,
    )

    # Step 4: 基因特征工程
    feature_df = compute_gene_features(
        adata, pb_counts, pb_meta,
        cell_type_key   = cell_type_key,
        condition_key   = condition_key,
        condition_case  = condition_case,
        condition_ctrl  = condition_ctrl,
    )
    feature_df.to_csv(os.path.join(output_dir, "gene_features.csv"))
    log.info("  特征矩阵保存至: %s/gene_features.csv", output_dir)

    # Step 5: HPO 基因集 & 标签
    hpo_df = load_hpo_gene_sets(hpo_gene_path)
    hpo_gene_sets = get_hpo_gene_sets(hpo_df)

    # 如果没有指定目标 HPO，使用 Toy HPO（取前 5 个 term）
    if target_hpo_ids is None:
        target_hpo_ids = list(hpo_gene_sets.keys())[:5]
        log.info("  未指定目标 HPO，使用演示 HPO IDs: %s", target_hpo_ids)

    labels = build_gene_labels(
        genes          = feature_df.index.tolist(),
        hpo_gene_sets  = hpo_gene_sets,
        target_hpo_ids = target_hpo_ids,
    )
    labels.to_csv(os.path.join(output_dir, "gene_labels.csv"))

    # Step 6: 训练分类器
    try:
        clf, metrics, feat_imp = train_classifier(feature_df, labels)
        feat_imp.to_csv(os.path.join(output_dir, "feature_importance.csv"),
                        index=False)
        log.info("  特征重要性保存至: %s/feature_importance.csv", output_dir)
    except Exception as e:
        log.error("分类器训练失败: %s", e)
        return {}

    # Step 7: 排序
    gene_ranking = predict_and_rank(clf, feature_df, top_k=top_k)
    out_file = os.path.join(output_dir, "gene_ranking_B.csv")
    gene_ranking.to_csv(out_file, index=False)
    log.info("[Done] 基因排序结果保存至: %s", out_file)

    # 打印 top-20
    log.info("\n=== Top-20 候选基因（监督学习表型关联概率） ===\n%s",
             gene_ranking.head(20).to_string(index=False))

    return {
        "feature_df":   feature_df,
        "labels":       labels,
        "model":        clf,
        "metrics":      metrics,
        "gene_ranking": gene_ranking,
    }


# ==============================================================================
# 脚本入口 Script Entry Point
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline B: Supervised gene ranking from scRNA-seq features"
    )
    parser.add_argument("--adata",      default=None, help="AnnData h5ad 路径")
    parser.add_argument("--hpo_genes",  default=None, help="HPO phenotype_to_genes.txt 路径")
    parser.add_argument("--hpo_ids",    default=None, nargs="+",
                        help="目标 HPO term IDs，如 HP:0000118 HP:0001250")
    parser.add_argument("--cell_type_key",   default=CELL_TYPE_KEY)
    parser.add_argument("--sample_key",      default=SAMPLE_KEY)
    parser.add_argument("--condition_key",   default=CONDITION_KEY)
    parser.add_argument("--condition_case",  default=CONDITION_CASE)
    parser.add_argument("--condition_ctrl",  default=CONDITION_CTRL)
    parser.add_argument("--output_dir",      default=OUTPUT_DIR)
    parser.add_argument("--top_k",     type=int, default=N_TOP_GENES)
    args = parser.parse_args()

    run_pipeline_B(
        adata_path      = args.adata,
        hpo_gene_path   = args.hpo_genes,
        target_hpo_ids  = args.hpo_ids,
        cell_type_key   = args.cell_type_key,
        sample_key      = args.sample_key,
        condition_key   = args.condition_key,
        condition_case  = args.condition_case,
        condition_ctrl  = args.condition_ctrl,
        output_dir      = args.output_dir,
        top_k           = args.top_k,
    )


if __name__ == "__main__":
    main()
