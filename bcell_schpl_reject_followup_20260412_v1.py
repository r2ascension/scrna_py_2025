#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B-cell scHPL rejected-cell follow-up analysis.

This is the B-cell counterpart of the stromal reject-followup workflow:
1. summarize where rejected cells accumulate (label / cluster),
2. compare transferred labels to scHPL outcomes,
3. visualize rejected vs accepted cells on UMAP,
4. run rejected-vs-accepted DE within candidate labels,
5. write a concise report for downstream interpretation.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

for _k in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(_k, "4")

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy import sparse

PIPELINE_VERSION = "v1.0-bcell-reject-followup"
PIPELINE_DATE = "2026-04-12"
RANDOM_SEED = 42

DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/py/0412/bcell_merge_schpl_unified_v1_0/"
    "bcell_reference_plus_query_schpl_v1_0_unified_20260412.h5ad"
)
DEFAULT_CONFIG_JSON = Path(
    "/home/h2048/data/py/0412/bcell_merge_schpl_unified_v1_0/"
    "config_v1_0_unified_20260412.json"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/h2048/data/py/0412/bcell_schpl_reject_followup_v1_0"
)

COL_SOURCE = "data_source"
COL_BASE_PRED = "Cell_Type_L2_pred"
COL_BASE_FINAL = "Cell_Type_L2_final"
COL_BASE_CONF = "mapping_confidence"
COL_SCHPL_RAW = "schpl_pred_raw"
COL_SCHPL_PRED = "schpl_pred"
COL_SCHPL_PROB = "schpl_prob"
COL_REJECTED = "schpl_rejected"
COL_REJECT_TYPE = "schpl_reject_type"
COL_CLUSTER = "leiden_schpl_qc"
COL_NOVEL = "schpl_novel_candidate"
UMAP_KEY = "X_umap"

BCELL_MARKERS = [
    "CD19", "MS4A1", "CD79A", "CD79B", "CD27", "IGHD", "IGHM",
    "MZB1", "JCHAIN", "SDC1", "TCL1A", "HLA-DRA", "CD74",
    "BANK1", "AICDA", "XBP1", "PRDM1",
]

WEB_REFERENCES = [
    {
        "title": "scHPL tips",
        "url": "https://schpl.readthedocs.io/en/latest/scHPL_tips.html",
        "note": "kNN is recommended for lower-dimensional integrated embeddings such as scVI/scArches latents.",
    },
    {
        "title": "treeArches basic tutorial",
        "url": "https://schpl.readthedocs.io/en/latest/treeArches_pbmc.html",
        "note": "After prediction, compare transferred labels with scHPL outcomes using heatmap / visualization.",
    },
    {
        "title": "treeArches identifying new cell types",
        "url": "https://schpl.readthedocs.io/en/latest/treeArches_identifying_new_ct.html",
        "note": "When one query label is enriched for rejected cells, compare rejected vs accepted cells with markers / DE.",
    },
]


@dataclass
class AnalysisConfig:
    input_h5ad: Path
    config_json: Path
    output_dir: Path
    source_key: str = COL_SOURCE
    query_source_value: str = "query"
    base_pred_key: str = COL_BASE_PRED
    base_final_key: str = COL_BASE_FINAL
    base_conf_key: str = COL_BASE_CONF
    min_cells_per_label: int = 100
    min_rejected_per_label: int = 20
    min_reject_rate_pct: float = 5.0
    min_cells_per_cluster: int = 100
    min_rejected_per_cluster: int = 20
    min_cluster_reject_rate_pct: float = 10.0
    max_de_labels: int = 5
    top_n_de_genes: int = 50


def set_seed(seed: int = RANDOM_SEED) -> None:
    np.random.seed(seed)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=150, dpi_save=300, frameon=False, facecolor="white")
    sc.settings.n_jobs = 4
    sns.set_context("talk")


def slugify(text: str) -> str:
    safe = [c if c.isalnum() else "_" for c in str(text)]
    out = "".join(safe)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")[:120] or "unknown"


def ensure_dirs(cfg: AnalysisConfig) -> Dict[str, Path]:
    out = cfg.output_dir
    dirs = {
        "root": out,
        "fig": out / "figures",
        "fig_de": out / "figures" / "de",
        "tables": out / "tables",
        "de": out / "de",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def load_inputs(cfg: AnalysisConfig) -> tuple[AnnData, dict]:
    if not cfg.input_h5ad.exists():
        raise FileNotFoundError(f"[ERROR] Input h5ad not found: {cfg.input_h5ad}")
    if not cfg.config_json.exists():
        raise FileNotFoundError(f"[ERROR] Config JSON not found: {cfg.config_json}")

    adata = sc.read_h5ad(cfg.input_h5ad)
    with open(cfg.config_json) as f:
        meta = json.load(f)

    required_obs = [
        cfg.source_key,
        cfg.base_pred_key,
        cfg.base_final_key,
        cfg.base_conf_key,
        COL_SCHPL_PRED,
        COL_SCHPL_RAW,
        COL_SCHPL_PROB,
        COL_REJECTED,
        COL_REJECT_TYPE,
        COL_CLUSTER,
        COL_NOVEL,
    ]
    missing = [c for c in required_obs if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"[ERROR] Missing required obs columns: {missing}")
    if UMAP_KEY not in adata.obsm:
        raise KeyError(f"[ERROR] Missing obsm['{UMAP_KEY}']")

    return adata, meta


def add_helper_columns(adata: AnnData, cfg: AnalysisConfig) -> AnnData:
    qry_mask = adata.obs[cfg.source_key].astype(str) == cfg.query_source_value
    adata = adata[qry_mask].copy()
    adata.obs[cfg.base_pred_key] = adata.obs[cfg.base_pred_key].astype(str)
    adata.obs[cfg.base_final_key] = adata.obs[cfg.base_final_key].astype(str)
    adata.obs[COL_SCHPL_PRED] = adata.obs[COL_SCHPL_PRED].astype(str)
    adata.obs[COL_SCHPL_RAW] = adata.obs[COL_SCHPL_RAW].astype(str)
    adata.obs[COL_REJECT_TYPE] = adata.obs[COL_REJECT_TYPE].astype(str)
    adata.obs["reject_group"] = np.where(adata.obs[COL_REJECTED].to_numpy(), "Rejected", "Accepted")
    adata.obs["schpl_pred_collapsed"] = np.where(
        adata.obs[COL_REJECTED].to_numpy(),
        "Rejected",
        adata.obs[COL_SCHPL_PRED].astype(str).to_numpy(),
    )
    return adata


def summarize_labels(adata: AnnData, cfg: AnalysisConfig) -> pd.DataFrame:
    df = (
        adata.obs.groupby(cfg.base_final_key, observed=False)
        .agg(
            n_total=(COL_SCHPL_PRED, "size"),
            n_rejected=(COL_REJECTED, "sum"),
            n_accepted=("reject_group", lambda s: int((s == "Accepted").sum())),
        )
        .reset_index()
    )
    df["rej_rate_pct"] = np.where(df["n_total"] > 0, df["n_rejected"] / df["n_total"] * 100.0, np.nan)
    return df.sort_values(["n_rejected", "rej_rate_pct"], ascending=[False, False])


def summarize_clusters(adata: AnnData) -> pd.DataFrame:
    df = (
        adata.obs.groupby(COL_CLUSTER, observed=False)
        .agg(
            n_total=(COL_SCHPL_PRED, "size"),
            n_rejected=(COL_REJECTED, "sum"),
            n_novel=(COL_NOVEL, "sum"),
        )
        .reset_index()
    )
    df["n_accepted"] = df["n_total"] - df["n_rejected"]
    df["rej_rate_pct"] = np.where(df["n_total"] > 0, df["n_rejected"] / df["n_total"] * 100.0, np.nan)
    return df.sort_values(["rej_rate_pct", "n_rejected"], ascending=[False, False])


def select_label_candidates(df: pd.DataFrame, cfg: AnalysisConfig) -> pd.DataFrame:
    keep = df[
        (df[cfg.base_final_key] != "Unknown")
        & (df["n_total"] >= cfg.min_cells_per_label)
        & (df["n_rejected"] >= cfg.min_rejected_per_label)
        & (df["rej_rate_pct"] >= cfg.min_reject_rate_pct)
    ].copy()
    return keep.sort_values(["n_rejected", "rej_rate_pct"], ascending=[False, False])


def select_cluster_candidates(df: pd.DataFrame, cfg: AnalysisConfig) -> pd.DataFrame:
    keep = df[
        (df["n_total"] >= cfg.min_cells_per_cluster)
        & (df["n_rejected"] >= cfg.min_rejected_per_cluster)
        & (df["rej_rate_pct"] >= cfg.min_cluster_reject_rate_pct)
    ].copy()
    return keep.sort_values(["rej_rate_pct", "n_rejected"], ascending=[False, False])


def plot_overall_umap(adata: AnnData, cfg: AnalysisConfig, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(18, 15))
    axes = axes.ravel()
    sc.pl.umap(adata, color=cfg.base_final_key, ax=axes[0], show=False, frameon=False, size=6,
               legend_loc="right margin", legend_fontsize=6, title="Base model final label")
    sc.pl.umap(adata, color=COL_SCHPL_PRED, ax=axes[1], show=False, frameon=False, size=6,
               legend_loc="right margin", legend_fontsize=6, title="scHPL prediction")
    sc.pl.umap(adata, color=COL_REJECT_TYPE, ax=axes[2], show=False, frameon=False, size=6,
               legend_loc="right margin", legend_fontsize=7, title="Reject type")
    sc.pl.umap(adata, color=COL_SCHPL_PROB, ax=axes[3], show=False, frameon=False, size=6,
               cmap="viridis", vmin=0, vmax=1, title="scHPL probability")
    plt.suptitle("B-cell scHPL rejected-cell follow-up: overall view", fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")


def plot_transition_heatmap(adata: AnnData, cfg: AnalysisConfig, out_path: Path, max_rows: int = 20) -> None:
    ct = pd.crosstab(
        adata.obs[cfg.base_final_key].astype(str),
        adata.obs["schpl_pred_collapsed"].astype(str),
        normalize="index",
    )
    row_order = adata.obs[cfg.base_final_key].astype(str).value_counts().head(max_rows).index
    ct = ct.loc[[r for r in row_order if r in ct.index]]
    fig, ax = plt.subplots(figsize=(max(8, ct.shape[1] * 0.8), max(6, ct.shape[0] * 0.5)))
    sns.heatmap(ct, cmap="Reds", linewidths=0.3, cbar_kws={"label": "row fraction"}, ax=ax)
    ax.set_title("B cell: base final vs scHPL outcome")
    ax.set_xlabel("scHPL outcome")
    ax.set_ylabel("Base final label")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")


def make_expr_subset(adata_sub: AnnData) -> AnnData:
    if adata_sub.raw is not None:
        X = adata_sub.raw.X
        var = adata_sub.raw.var.copy()
    elif "counts" in adata_sub.layers:
        X = adata_sub.layers["counts"]
        var = adata_sub.var.copy()
    else:
        X = adata_sub.X
        var = adata_sub.var.copy()

    X = sparse.csr_matrix(X) if sparse.issparse(X) else np.asarray(X, dtype=np.float32)
    expr = AnnData(X=X.copy(), obs=adata_sub.obs.copy(), var=var)
    expr.var_names_make_unique()

    min_cells = max(3, min(10, int(np.ceil(expr.n_obs * 0.01))))
    sc.pp.filter_genes(expr, min_cells=min_cells)
    sc.pp.normalize_total(expr, target_sum=1e4)
    sc.pp.log1p(expr)
    return expr


def save_rank_genes(expr: AnnData, out_csv: Path, top_n: int) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        expr,
        groupby="reject_group",
        groups=["Rejected"],
        reference="Accepted",
        method="wilcoxon",
        pts=True,
    )
    df = sc.get.rank_genes_groups_df(expr, group="Rejected")
    df.to_csv(out_csv, index=False)
    return df.head(top_n).copy()


def plot_reject_vs_accept_dotplot(expr: AnnData, label: str, top_de: List[str], out_path: Path) -> None:
    genes = []
    for g in BCELL_MARKERS + top_de:
        if g in expr.var_names and g not in genes:
            genes.append(g)
    genes = genes[:20]
    if not genes:
        return
    fig = sc.pl.dotplot(
        expr,
        var_names=genes,
        groupby="reject_group",
        standard_scale="var",
        colorbar_title="Scaled\nexpression",
        show=False,
        return_fig=True,
        title=f"B cell | {label}\nRejected vs Accepted",
    )
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")


def plot_candidate_umap(adata_sub: AnnData, label: str, out_path: Path) -> None:
    ad = adata_sub.copy()
    ad.obs["candidate_group"] = np.where(ad.obs[COL_REJECTED].to_numpy(), f"{label} | Rejected", f"{label} | Accepted")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    sc.pl.umap(ad, color="candidate_group", ax=axes[0], show=False, frameon=False, size=12,
               title=f"{label}: rejected vs accepted")
    sc.pl.umap(ad, color=COL_SCHPL_PROB, ax=axes[1], show=False, frameon=False, size=12,
               cmap="viridis", vmin=0, vmax=1, title=f"{label}: scHPL probability")
    plt.suptitle("B-cell candidate label follow-up", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")


def run_de_followup(adata: AnnData, label_candidates: pd.DataFrame, dirs: Dict[str, Path], cfg: AnalysisConfig) -> pd.DataFrame:
    rows = []
    for _, row in label_candidates.head(cfg.max_de_labels).iterrows():
        label = row[cfg.base_final_key]
        mask = adata.obs[cfg.base_final_key].astype(str) == str(label)
        ad_sub = adata[mask].copy()
        n_rej = int(ad_sub.obs[COL_REJECTED].sum())
        n_acc = int((~ad_sub.obs[COL_REJECTED].to_numpy()).sum())
        if n_rej < 20 or n_acc < 20:
            rows.append({
                "label": label,
                "n_total": int(ad_sub.n_obs),
                "n_rejected": n_rej,
                "n_accepted": n_acc,
                "status": "skip_too_few_cells",
                "de_csv": None,
                "dotplot": None,
                "umap": None,
            })
            continue

        expr = make_expr_subset(ad_sub)
        expr.obs["reject_group"] = ad_sub.obs["reject_group"].values

        stem = slugify(label)
        de_csv = dirs["de"] / f"de_rejected_vs_accepted__{stem}.csv"
        top_df = save_rank_genes(expr, de_csv, cfg.top_n_de_genes)
        top_genes = top_df["names"].dropna().astype(str).head(10).tolist()

        dotplot_path = dirs["fig_de"] / f"dotplot_rejected_vs_accepted__{stem}.pdf"
        plot_reject_vs_accept_dotplot(expr, label, top_genes, dotplot_path)

        umap_path = dirs["fig_de"] / f"umap_rejected_vs_accepted__{stem}.pdf"
        plot_candidate_umap(ad_sub, label, umap_path)

        rows.append({
            "label": label,
            "n_total": int(ad_sub.n_obs),
            "n_rejected": n_rej,
            "n_accepted": n_acc,
            "status": "ok",
            "de_csv": str(de_csv),
            "dotplot": str(dotplot_path),
            "umap": str(umap_path),
            "top_de_genes": ", ".join(top_genes[:8]),
        })
    return pd.DataFrame(rows)


def write_report(
    cfg: AnalysisConfig,
    meta: dict,
    label_df: pd.DataFrame,
    label_candidates: pd.DataFrame,
    cluster_candidates: pd.DataFrame,
    de_df: pd.DataFrame,
    out_path: Path,
) -> None:
    lines: List[str] = []
    lines.append(f"# B-cell scHPL reject follow-up report ({PIPELINE_VERSION})")
    lines.append("")
    lines.append(f"- Date: {PIPELINE_DATE}")
    lines.append(f"- Input h5ad: `{cfg.input_h5ad}`")
    lines.append(f"- Source config: `{cfg.config_json}`")
    lines.append("")
    lines.append("## 官方建议（基于 scHPL / treeArches 文档）")
    lines.append("")
    lines.append("scHPL/treeArches 的官方教程不建议在看到 reject 之后就停下；推荐的下一步是：")
    lines.append("")
    lines.append("1. 先看 rejected 细胞是不是集中在某些 query 标签或 cluster。")
    lines.append("2. 比较 transferred labels 与 scHPL 预测结果（heatmap、可视化）。")
    lines.append("3. 如果某个 query 标签整体匹配到 B-cell hierarchy，但其中一部分持续被 reject，优先怀疑新状态/亚群，而不是直接判定技术噪音。")
    lines.append("4. 在该标签内部，把 `Rejected` 与 `Accepted` 拆开做 DE 和 marker 验证。")
    lines.append("")
    lines.append("参考链接：")
    for ref in WEB_REFERENCES:
        lines.append(f"- [{ref['title']}]({ref['url']}): {ref['note']}")
    lines.append("")
    lines.append("## 本次 B-cell 数据的 reject 热点")
    lines.append("")
    for _, row in label_df.head(15).iterrows():
        lines.append(
            f"- `{row[cfg.base_final_key]}`: total={int(row['n_total']):,}, rejected={int(row['n_rejected']):,}, rej_rate={row['rej_rate_pct']:.1f}%"
        )
    lines.append("")
    if not label_candidates.empty:
        lines.append("## 高 reject 候选标签")
        lines.append("")
        for _, row in label_candidates.iterrows():
            lines.append(
                f"- `{row[cfg.base_final_key]}`: n={int(row['n_total']):,}, rejected={int(row['n_rejected']):,}, rej_rate={row['rej_rate_pct']:.1f}%"
            )
        lines.append("")
    if not cluster_candidates.empty:
        lines.append("## 高 reject 候选 cluster")
        lines.append("")
        for _, row in cluster_candidates.head(15).iterrows():
            lines.append(
                f"- cluster `{row[COL_CLUSTER]}`: n={int(row['n_total']):,}, rejected={int(row['n_rejected']):,}, rej_rate={row['rej_rate_pct']:.1f}%"
            )
        lines.append("")
    lines.append("## 已生成的 follow-up 分析")
    lines.append("")
    lines.append(f"- label / cluster 汇总表：`{cfg.output_dir / 'tables'}`")
    lines.append(f"- UMAP 和 heatmap：`{cfg.output_dir / 'figures'}`")
    lines.append(f"- rejected vs accepted DE：`{cfg.output_dir / 'de'}`")
    lines.append("")
    if not de_df.empty:
        lines.append("## rejected-vs-accepted DE 已完成的候选标签")
        lines.append("")
        for _, row in de_df.iterrows():
            if row.get("status") != "ok":
                continue
            lines.append(f"- `{row['label']}`: top DE genes = {row.get('top_de_genes', '')}")
        lines.append("")
    lines.append("## 推荐下一步")
    lines.append("")
    lines.append("1. 先看 `candidate_labels.csv` 锁定 reject 富集的 B-cell 标签。")
    lines.append("2. 看 `transition_heatmap.pdf` 判断是整体不匹配还是局部 reject。")
    lines.append("3. 看 `de_rejected_vs_accepted__*.csv` 和 dotplot，判断 rejected 是否表现为浆细胞化、GC 激活、干扰素/stress、或非 B 污染程序。")
    lines.append("4. 只有在 rejected 形成稳定表达程序并有局部空间位置时，再考虑定义为新状态/亚群。")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def save_json(obj: dict, path: Path) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def parse_args() -> AnalysisConfig:
    parser = argparse.ArgumentParser(description="Follow-up analysis for B-cell scHPL rejected cells")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--config-json", type=Path, default=DEFAULT_CONFIG_JSON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-cells-per-label", type=int, default=100)
    parser.add_argument("--min-rejected-per-label", type=int, default=20)
    parser.add_argument("--min-reject-rate-pct", type=float, default=5.0)
    parser.add_argument("--max-de-labels", type=int, default=5)
    parser.add_argument("--top-n-de-genes", type=int, default=50)
    args = parser.parse_args()
    return AnalysisConfig(
        input_h5ad=args.input_h5ad,
        config_json=args.config_json,
        output_dir=args.output_dir,
        min_cells_per_label=args.min_cells_per_label,
        min_rejected_per_label=args.min_rejected_per_label,
        min_reject_rate_pct=args.min_reject_rate_pct,
        max_de_labels=args.max_de_labels,
        top_n_de_genes=args.top_n_de_genes,
    )


def main() -> None:
    cfg = parse_args()
    set_seed()
    dirs = ensure_dirs(cfg)
    adata, meta = load_inputs(cfg)

    cfg.source_key = meta.get("source_key", cfg.source_key)
    cfg.query_source_value = meta.get("query_source_value", cfg.query_source_value)
    cfg.base_pred_key = meta.get("query_pred_key", cfg.base_pred_key)
    cfg.base_final_key = meta.get("query_final_key", cfg.base_final_key)
    cfg.base_conf_key = meta.get("query_conf_key", cfg.base_conf_key)

    adata = add_helper_columns(adata, cfg)

    print("=" * 80)
    print("B-cell scHPL rejected-cell follow-up")
    print("=" * 80)
    print(f"Input : {cfg.input_h5ad}")
    print(f"Output: {cfg.output_dir}")
    print(f"Query-only shape : {adata.shape}")
    print(f"Raw genes        : {None if adata.raw is None else adata.raw.n_vars}")

    label_df = summarize_labels(adata, cfg)
    cluster_df = summarize_clusters(adata)
    label_candidates = select_label_candidates(label_df, cfg)
    cluster_candidates = select_cluster_candidates(cluster_df, cfg)

    label_df.to_csv(dirs["tables"] / "label_rejection_summary.csv", index=False)
    cluster_df.to_csv(dirs["tables"] / "cluster_rejection_summary.csv", index=False)
    label_candidates.to_csv(dirs["tables"] / "candidate_labels.csv", index=False)
    cluster_candidates.to_csv(dirs["tables"] / "candidate_clusters.csv", index=False)

    plot_overall_umap(adata, cfg, dirs["fig"] / "overall_reject_followup_umap.pdf")
    plot_transition_heatmap(adata, cfg, dirs["fig"] / "transition_heatmap.pdf")

    de_df = run_de_followup(adata, label_candidates, dirs, cfg)
    de_df.to_csv(dirs["tables"] / "de_followup_summary.csv", index=False)

    report_path = cfg.output_dir / "report.md"
    write_report(cfg, meta, label_df, label_candidates, cluster_candidates, de_df, report_path)

    save_json(
        {
            "version": PIPELINE_VERSION,
            "date": PIPELINE_DATE,
            "input_h5ad": str(cfg.input_h5ad),
            "config_json": str(cfg.config_json),
            "output_dir": str(cfg.output_dir),
            "candidate_label_thresholds": {
                "min_cells_per_label": cfg.min_cells_per_label,
                "min_rejected_per_label": cfg.min_rejected_per_label,
                "min_reject_rate_pct": cfg.min_reject_rate_pct,
            },
            "candidate_cluster_thresholds": {
                "min_cells_per_cluster": cfg.min_cells_per_cluster,
                "min_rejected_per_cluster": cfg.min_rejected_per_cluster,
                "min_cluster_reject_rate_pct": cfg.min_cluster_reject_rate_pct,
            },
            "label_summary": label_df.to_dict(orient="records"),
        },
        cfg.output_dir / "analysis_config.json",
    )

    print("\nLabel rejection summary:")
    print(label_df.head(12).to_string(index=False))
    print("\nTop candidate labels:")
    if label_candidates.empty:
        print("  [INFO] No label candidates passed thresholds")
    else:
        print(label_candidates.head(12).to_string(index=False))

    print("\nTop candidate clusters:")
    if cluster_candidates.empty:
        print("  [INFO] No cluster candidates passed thresholds")
    else:
        print(cluster_candidates.head(12).to_string(index=False))

    print("\n[OK] report:", report_path)
    print("[OK] tables:", dirs["tables"])
    print("[OK] figures:", dirs["fig"])
    print("[OK] de:", dirs["de"])


if __name__ == "__main__":
    main()
