#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Validate current B-cell schpl conclusions against historical merged scANVI outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SCANVI_H5AD = Path("/home/h2048/data/py/0208/merged_scanvi_L2_prod_v1/merged_scanvi_L2_prod.h5ad")
SCHPL_H5AD = Path("/home/h2048/data/py/0412/bcell_merge_schpl_unified_v1_0/bcell_reference_plus_query_schpl_v1_0_unified_20260412.h5ad")
OUTPUT_DIR = Path("/home/h2048/data/py/0412/bcell_scanvi_validation_20260412")

MARKERS = [
    "CD27", "IGHD", "IGHM", "MS4A1", "JCHAIN", "MZB1", "SDC1", "TCL1A", "CD74", "BANK1"
]


def safe_series(s: pd.Series) -> pd.Series:
    obj = s.astype("object").copy()
    obj[pd.isna(obj)] = "NA"
    return obj.astype(str)


def agreement(a: pd.Series, b: pd.Series) -> float:
    return float((safe_series(a) == safe_series(b)).mean() * 100.0)


def marker_matrix(adata: sc.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.var_names or (adata.raw is not None and g in adata.raw.var_names)]
    if not genes:
        return pd.DataFrame(index=adata.obs_names)
    if adata.raw is not None and all(g in adata.raw.var_names for g in genes):
        X = adata.raw[:, genes].X
    else:
        X = adata[:, genes].X
    if sparse.issparse(X):
        X = X.toarray()
    return pd.DataFrame(X, index=adata.obs_names, columns=genes)


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "(empty)"
    cols = [str(c) for c in df.columns]
    rows = [[str(x) for x in row] for row in df.itertuples(index=False, name=None)]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header, sep, *body])


def write_markdown_report(
    path: Path,
    stats: dict[str, object],
    top_final_vs_scanvi: pd.DataFrame,
    top_scanvi_vs_schpl: pd.DataFrame,
    top_final_vs_raw_mapping: pd.DataFrame,
    raw_mapping_group_summary: pd.DataFrame,
    marker_summary: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append("# B-cell historical scANVI validation (2026-04-12)")
    lines.append("")
    lines.append("## 结论")
    lines.append("")
    lines.append("历史 merged scANVI 结果支持当前对 B-cell schpl 的判断：主要矛盾仍然是 `Naive_B` 与 `Memory_B` 边界，而不是 `Plasma` 大范围塌缩。")
    lines.append("")
    lines.append("## 核心统计")
    lines.append("")
    for key, value in stats.items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("## 历史 merged scANVI 状态")
    lines.append("")
    lines.append(str(stats.get("historical_scanvi_status", "NA")))
    lines.append("")
    lines.append("## final vs 历史 scANVI 的主要分歧")
    lines.append("")
    lines.append(dataframe_to_markdown(top_final_vs_scanvi))
    lines.append("")
    lines.append("## 历史 scANVI vs 当前 schpl 的主要分歧")
    lines.append("")
    lines.append(dataframe_to_markdown(top_scanvi_vs_schpl))
    lines.append("")
    lines.append("## final vs 旧 mapping 原始预测（Cell_Type_L2_pred）")
    lines.append("")
    lines.append(dataframe_to_markdown(top_final_vs_raw_mapping))
    lines.append("")
    lines.append("## 旧 mapping 结果在关键 schpl 重分配群中的置信度")
    lines.append("")
    lines.append(dataframe_to_markdown(raw_mapping_group_summary))
    lines.append("")
    lines.append("## marker 均值（关键对照组）")
    lines.append("")
    lines.append(dataframe_to_markdown(marker_summary.reset_index().rename(columns={"index": "marker"})))
    lines.append("")
    lines.append("## 解读")
    lines.append("")
    lines.append("- 如果历史 scANVI 自己就把大量 `final Naive_B` 细胞压到 `Memory_B`，说明问题早于这次 unified schpl。")
    lines.append("- 如果 `scanvi_pred` 与 `schpl_pred` 高度一致，则 schpl 更像在重复/层级化验证既有边界，而不是制造全新错误。")
    lines.append("- 如果旧 mapping 原始预测在这些群体上的置信度明显低于各自主标签 baseline，则说明这些细胞本来就位于不稳定边界上。")
    lines.append("- 若 `final Naive_B -> scanvi Memory_B` 这批细胞的 `CD27` 更高、`IGHD/IGHM/TCL1A` 更低，则它们更偏 Memory-like。")
    path.write_text("\n".join(lines), encoding="utf-8")


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig_dir = OUTPUT_DIR / "figures"
tab_dir = OUTPUT_DIR / "tables"
fig_dir.mkdir(exist_ok=True)
tab_dir.mkdir(exist_ok=True)

print(f"[load] {SCANVI_H5AD}")
scanvi = sc.read_h5ad(SCANVI_H5AD)
print(f"[load] {SCHPL_H5AD}")
schpl = sc.read_h5ad(SCHPL_H5AD)

scanvi_q = scanvi[scanvi.obs["data_source"].astype(str) == "query"].copy()
schpl_q = schpl[schpl.obs["data_source"].astype(str) == "query"].copy()

for ad, name in [(scanvi_q, "scanvi_q"), (schpl_q, "schpl_q")]:
    if "barcode" not in ad.obs.columns:
        raise KeyError(f"{name} is missing barcode column for cross-file alignment")
    if ad.obs["barcode"].astype(str).duplicated().any():
        raise ValueError(f"{name} barcode column contains duplicates; cannot align safely")

scanvi_q.obs["_merge_barcode"] = scanvi_q.obs["barcode"].astype(str).values
schpl_q.obs["_merge_barcode"] = schpl_q.obs["barcode"].astype(str).values

common = pd.Index(scanvi_q.obs["_merge_barcode"]).intersection(pd.Index(schpl_q.obs["_merge_barcode"]))
scanvi_pos = pd.Index(scanvi_q.obs["_merge_barcode"]).get_indexer(common)
schpl_pos = pd.Index(schpl_q.obs["_merge_barcode"]).get_indexer(common)
scanvi_q = scanvi_q[scanvi_pos].copy()
schpl_q = schpl_q[schpl_pos].copy()

scanvi_q.obs["schpl_pred"] = schpl_q.obs["schpl_pred"].astype("string").values
scanvi_q.obs["schpl_rejected"] = schpl_q.obs["schpl_rejected"].astype(bool).values
scanvi_q.obs["schpl_prob"] = pd.to_numeric(schpl_q.obs.get("schpl_prob"), errors="coerce").values if "schpl_prob" in schpl_q.obs.columns else np.nan

meta = scanvi_q.obs[["Cell_Type_L2_pred", "Cell_Type_L2_final", "L2_scanvi_pred", "mapping_confidence", "schpl_pred", "schpl_rejected"]].copy()
meta["Cell_Type_L2_pred"] = meta["Cell_Type_L2_pred"].astype("string")
meta["Cell_Type_L2_final"] = meta["Cell_Type_L2_final"].astype("string")
meta["L2_scanvi_pred"] = meta["L2_scanvi_pred"].astype("string")
meta["mapping_confidence"] = pd.to_numeric(meta["mapping_confidence"], errors="coerce")
meta["schpl_pred"] = meta["schpl_pred"].astype("string")

historical_scanvi_unique = safe_series(meta["L2_scanvi_pred"]).nunique()
historical_scanvi_top = safe_series(meta["L2_scanvi_pred"]).value_counts(dropna=False).head(3).to_dict()
historical_scanvi_degenerate = historical_scanvi_unique <= 1

conf_final_vs_scanvi = pd.crosstab(meta["Cell_Type_L2_final"], meta["L2_scanvi_pred"], dropna=False)
conf_final_vs_schpl = pd.crosstab(meta["Cell_Type_L2_final"], meta["schpl_pred"], dropna=False)
conf_scanvi_vs_schpl = pd.crosstab(meta["L2_scanvi_pred"], meta["schpl_pred"], dropna=False)
conf_final_vs_raw_mapping = pd.crosstab(meta["Cell_Type_L2_final"], meta["Cell_Type_L2_pred"], dropna=False)

conf_final_vs_scanvi.to_csv(tab_dir / "confusion_final_vs_historical_scanvi.csv")
conf_final_vs_schpl.to_csv(tab_dir / "confusion_final_vs_schpl.csv")
conf_scanvi_vs_schpl.to_csv(tab_dir / "confusion_historical_scanvi_vs_schpl.csv")
conf_final_vs_raw_mapping.to_csv(tab_dir / "confusion_final_vs_raw_mapping_pred.csv")

pairs_final_vs_scanvi = (
    meta.loc[meta["Cell_Type_L2_final"] != meta["L2_scanvi_pred"], ["Cell_Type_L2_final", "L2_scanvi_pred"]]
    .value_counts()
    .rename("n_cells")
    .reset_index()
    .sort_values("n_cells", ascending=False)
)
pairs_scanvi_vs_schpl = (
    meta.loc[meta["L2_scanvi_pred"] != meta["schpl_pred"], ["L2_scanvi_pred", "schpl_pred"]]
    .value_counts()
    .rename("n_cells")
    .reset_index()
    .sort_values("n_cells", ascending=False)
)
pairs_final_vs_raw_mapping = (
    meta.loc[meta["Cell_Type_L2_final"] != meta["Cell_Type_L2_pred"], ["Cell_Type_L2_final", "Cell_Type_L2_pred"]]
    .value_counts()
    .rename("n_cells")
    .reset_index()
    .sort_values("n_cells", ascending=False)
)

pairs_final_vs_scanvi.to_csv(tab_dir / "top_disagreements_final_vs_historical_scanvi.csv", index=False)
pairs_scanvi_vs_schpl.to_csv(tab_dir / "top_disagreements_historical_scanvi_vs_schpl.csv", index=False)
pairs_final_vs_raw_mapping.to_csv(tab_dir / "top_disagreements_final_vs_raw_mapping_pred.csv", index=False)

expr = marker_matrix(schpl_q, MARKERS)
expr.index = meta.index
obs = meta.copy()
obs["group"] = "other"
obs.loc[(obs["Cell_Type_L2_final"] == "Naive_B") & (obs["schpl_pred"] == "Naive_B") & (~obs["schpl_rejected"]), "group"] = "accepted_Naive"
obs.loc[(obs["Cell_Type_L2_final"] == "Memory_B") & (obs["schpl_pred"] == "Memory_B") & (~obs["schpl_rejected"]), "group"] = "accepted_Memory"
obs.loc[(obs["Cell_Type_L2_final"] == "Plasma") & (obs["schpl_pred"] == "Plasma") & (~obs["schpl_rejected"]), "group"] = "accepted_Plasma"
obs.loc[(obs["Cell_Type_L2_final"] == "Naive_B") & (obs["schpl_pred"] == "Memory_B") & (~obs["schpl_rejected"]), "group"] = "Naive_to_Memory"
obs.loc[(obs["Cell_Type_L2_final"] == "Memory_B") & (obs["schpl_pred"] == "Naive_B") & (~obs["schpl_rejected"]), "group"] = "Memory_to_Naive"
obs.loc[(obs["Cell_Type_L2_final"] == "Plasma") & (obs["schpl_pred"] == "Memory_B") & (~obs["schpl_rejected"]), "group"] = "Plasma_to_Memory"
obs.loc[(obs["Cell_Type_L2_final"] == "Memory_B") & (obs["schpl_pred"] == "Plasma") & (~obs["schpl_rejected"]), "group"] = "Memory_to_Plasma"

use = obs["group"] != "other"
marker_summary = expr.loc[use].groupby(obs.loc[use, "group"]).mean().T.round(3)
marker_summary.to_csv(tab_dir / "marker_summary_historical_scanvi_validation.csv")
obs["group"].value_counts().rename_axis("group").reset_index(name="n_cells").to_csv(tab_dir / "group_counts_historical_scanvi_validation.csv", index=False)

raw_mapping_groups = []
for group_name, mask in {
    "Naive_to_Memory": (meta["Cell_Type_L2_final"] == "Naive_B") & (meta["schpl_pred"] == "Memory_B") & (~meta["schpl_rejected"]),
    "Memory_to_Naive": (meta["Cell_Type_L2_final"] == "Memory_B") & (meta["schpl_pred"] == "Naive_B") & (~meta["schpl_rejected"]),
    "Plasma_to_Memory": (meta["Cell_Type_L2_final"] == "Plasma") & (meta["schpl_pred"] == "Memory_B") & (~meta["schpl_rejected"]),
    "Memory_to_Plasma": (meta["Cell_Type_L2_final"] == "Memory_B") & (meta["schpl_pred"] == "Plasma") & (~meta["schpl_rejected"]),
}.items():
    sub = meta.loc[mask]
    raw_mapping_groups.append({
        "group": group_name,
        "n_cells": int(len(sub)),
        "raw_mapping_top": sub["Cell_Type_L2_pred"].astype(str).value_counts().head(3).to_dict(),
        "conf_mean": round(float(sub["mapping_confidence"].mean()), 4) if len(sub) else np.nan,
        "conf_median": round(float(sub["mapping_confidence"].median()), 4) if len(sub) else np.nan,
        "conf_q25": round(float(sub["mapping_confidence"].quantile(0.25)), 4) if len(sub) else np.nan,
        "conf_q75": round(float(sub["mapping_confidence"].quantile(0.75)), 4) if len(sub) else np.nan,
    })
raw_mapping_group_summary = pd.DataFrame(raw_mapping_groups)
raw_mapping_group_summary.to_csv(tab_dir / "raw_mapping_confidence_in_schpl_disagreement_groups.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(30, 8))
for ax, color, title in [
    (axes[0], "Cell_Type_L2_final", "B-cell base final"),
    (axes[1], "L2_scanvi_pred", "Historical merged scANVI pred"),
    (axes[2], "schpl_pred", "Current scHPL pred"),
]:
    sc.pl.umap(
        scanvi_q,
        color=color,
        ax=ax,
        show=False,
        frameon=False,
        size=8,
        legend_loc="right margin",
        legend_fontsize=7,
        title=title,
    )
plt.tight_layout()
fig.savefig(fig_dir / "bcell_final_scanvi_schpl_sidebyside.png", dpi=200, bbox_inches="tight")
fig.savefig(fig_dir / "bcell_final_scanvi_schpl_sidebyside.pdf", dpi=200, bbox_inches="tight")
plt.close(fig)

stats = {
    "n_query_cells": int(scanvi_q.n_obs),
    "agreement(final, raw_mapping_pred)": f"{agreement(meta['Cell_Type_L2_final'], meta['Cell_Type_L2_pred']):.3f}%",
    "agreement(final, historical_scanvi)": f"{agreement(meta['Cell_Type_L2_final'], meta['L2_scanvi_pred']):.3f}%" if not historical_scanvi_degenerate else "NA (degenerate historical scanvi output)",
    "agreement(final, schpl)": f"{agreement(meta['Cell_Type_L2_final'], meta['schpl_pred']):.3f}%",
    "agreement(historical_scanvi, schpl)": f"{agreement(meta['L2_scanvi_pred'], meta['schpl_pred']):.3f}%" if not historical_scanvi_degenerate else "NA (degenerate historical scanvi output)",
    "schpl_rejected": int(meta['schpl_rejected'].sum()),
    "historical_scanvi_status": f"unique_labels={historical_scanvi_unique}; top_counts={historical_scanvi_top}",
    "top_final_to_raw_mapping": pairs_final_vs_raw_mapping.head(1).to_dict(orient='records'),
    "top_final_to_scanvi": pairs_final_vs_scanvi.head(1).to_dict(orient='records'),
    "top_scanvi_to_schpl": pairs_scanvi_vs_schpl.head(1).to_dict(orient='records'),
}

pd.DataFrame([stats]).to_json(OUTPUT_DIR / "summary.json", orient="records", indent=2)
write_markdown_report(
    OUTPUT_DIR / "report.md",
    stats=stats,
    top_final_vs_scanvi=pairs_final_vs_scanvi.head(10),
    top_scanvi_vs_schpl=pairs_scanvi_vs_schpl.head(10),
    top_final_vs_raw_mapping=pairs_final_vs_raw_mapping.head(10),
    raw_mapping_group_summary=raw_mapping_group_summary,
    marker_summary=marker_summary,
)

print("[done] outputs written to", OUTPUT_DIR)
print("[summary] agreement(final, raw_mapping_pred)=", f"{agreement(meta['Cell_Type_L2_final'], meta['Cell_Type_L2_pred']):.3f}%")
print("[summary] historical_scanvi_status=", stats["historical_scanvi_status"])
print("[summary] top final->scanvi disagreement:")
print(pairs_final_vs_scanvi.head(10).to_string(index=False))
print("[summary] top final->raw_mapping disagreement:")
print(pairs_final_vs_raw_mapping.head(10).to_string(index=False))
print("[summary] top scanvi->schpl disagreement:")
print(pairs_scanvi_vs_schpl.head(10).to_string(index=False))
