#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PY_ROOT = Path(__file__).resolve().parents[2]
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from normal_airway_ml_common_20260527 import deep_get, ensure_dir, load_config, write_json  # noqa: E402

CONFIG_KEY = "bcell_internal_site_core_validation"
DEFAULT_CONFIG_PATH = "/home/h2048/script/config/normal_airway_ml_bcell_internal_site_core_validation_20260531.yaml"
DPI = 300


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create PNG/PDF visualizations for B-cell internal-only site-core validation.")
    parser.add_argument("--config", type=Path, default=Path(DEFAULT_CONFIG_PATH), help="YAML config path.")
    parser.add_argument("--run-root", type=Path, default=None, help="Optional explicit run root override.")
    return parser.parse_args()


def resolve_run_root(cfg: dict, run_root: Path | None) -> Path:
    if run_root is not None:
        return run_root
    output_root = Path(str(deep_get(cfg, "run", "output_root", default="/home/h2048/data/py/20260531")))
    prefix = str(deep_get(cfg, "run", "run_name_prefix", default="bcell_internal_site_core_validation_20260531"))
    return output_root / prefix


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def save_figure(fig: plt.Figure, stem: Path) -> list[str]:
    outputs: list[str] = []
    for suffix in (".pdf", ".png"):
        out = stem.with_suffix(suffix)
        fig.savefig(out, dpi=DPI, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def load_tables(run_root: Path) -> dict[str, pd.DataFrame]:
    return {
        "labels": pd.read_csv(run_root / "labels" / "query_label_transfer.tsv", sep="\t"),
        "qc": pd.read_csv(run_root / "labels" / "query_label_transfer_qc.tsv", sep="\t"),
        "sample_scores": pd.read_csv(run_root / "validation" / "sample_scores.tsv", sep="\t"),
        "retention": pd.read_csv(run_root / "validation" / "retention_summary.tsv", sep="\t"),
        "disease_shift": pd.read_csv(run_root / "validation" / "disease_shift_summary.tsv", sep="\t"),
        "panel_summary": pd.read_csv(run_root / "panels" / "panel_summary.tsv", sep="\t"),
    }


def plot_label_qc(labels: pd.DataFrame, figure_dir: Path, subtype_order: list[str]) -> list[str]:
    counts = (
        labels.groupby(["transferred_subtype_label", "keep_for_validation"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    order = [x for x in subtype_order if x in set(counts["transferred_subtype_label"].astype(str))]
    tail = [x for x in counts["transferred_subtype_label"].astype(str).unique().tolist() if x not in order]
    order.extend(sorted(tail))

    pivot = counts.pivot(index="transferred_subtype_label", columns="keep_for_validation", values="n_cells").fillna(0)
    for col in [True, False]:
        if col not in pivot.columns:
            pivot[col] = 0
    pivot = pivot.loc[order]

    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.6 * len(pivot))))
    y = np.arange(len(pivot))
    ax.barh(y, pivot[False], color="#cfcfcf", edgecolor="black", label="excluded")
    ax.barh(y, pivot[True], left=pivot[False], color="#4daf4a", edgecolor="black", label="kept_for_validation")
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index.astype(str))
    ax.set_xlabel("Cells")
    ax.set_title("Query label-transfer QC: kept vs excluded cells")
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.legend(loc="lower right")
    for yi, idx in enumerate(pivot.index):
        total = int(pivot.loc[idx, False] + pivot.loc[idx, True])
        ax.text(total + max(50, total * 0.01), yi, str(total), va="center", fontsize=8)
    fig.tight_layout()
    return save_figure(fig, figure_dir / "bcell_internal_label_transfer_qc")


def plot_retention_heatmap(retention: pd.DataFrame, figure_dir: Path, subtype_order: list[str], site_order: list[str]) -> list[str]:
    if retention.empty:
        return []
    plot_df = retention.copy()
    plot_df["subtype"] = pd.Categorical(plot_df["subtype"], categories=subtype_order, ordered=True)
    plot_df["site_anchor"] = pd.Categorical(plot_df["site_anchor"], categories=site_order, ordered=True)
    plot_df = plot_df.sort_values(["subtype", "site_anchor"]).reset_index(drop=True)

    matrix = plot_df.pivot(index="subtype", columns="site_anchor", values="delta_match_minus_other")
    matrix = matrix.reindex(index=[x for x in subtype_order if x in matrix.index], columns=[x for x in site_order if x in matrix.columns])
    values = matrix.to_numpy(dtype=float)

    vmax = np.nanmax(np.abs(values)) if np.isfinite(values).any() else 1.0
    vmax = max(vmax, 0.5)

    fig, ax = plt.subplots(figsize=(1.8 * matrix.shape[1] + 3, 0.8 * matrix.shape[0] + 3))
    im = ax.imshow(values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.astype(str))
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.astype(str))
    ax.set_title("Retention of healthy site-core signatures\n(sample mean score: matched tissue minus other tissues)")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Delta (matched - other)")
    fig.tight_layout()
    return save_figure(fig, figure_dir / "bcell_internal_retention_heatmap")


def plot_disease_shift_heatmap(disease_shift: pd.DataFrame, figure_dir: Path, subtype_order: list[str]) -> list[str]:
    if disease_shift.empty:
        return []
    plot_df = disease_shift.copy()
    plot_df["row_label"] = plot_df["subtype"].astype(str) + " | " + plot_df["site_anchor"].astype(str)

    disease_order = ["healthy", "CRSsNP", "CRSwNP", "asthma", "COVID-19"]
    existing = [x for x in disease_order if x in set(plot_df["disease_level_1"].astype(str))]
    existing.extend(sorted([x for x in plot_df["disease_level_1"].astype(str).unique().tolist() if x not in existing]))

    row_order = []
    for subtype in subtype_order:
        for _, row in plot_df[plot_df["subtype"].astype(str) == subtype].drop_duplicates(["subtype", "site_anchor"]).sort_values("site_anchor").iterrows():
            row_order.append(row["row_label"])
    row_order.extend([x for x in plot_df["row_label"].unique().tolist() if x not in row_order])

    matrix = plot_df.pivot(index="row_label", columns="disease_level_1", values="score_mean")
    matrix = matrix.reindex(index=row_order, columns=existing)
    values = matrix.to_numpy(dtype=float)
    vmin = np.nanmin(values) if np.isfinite(values).any() else -1.0
    vmax = np.nanmax(values) if np.isfinite(values).any() else 1.0
    if not np.isfinite(vmin):
        vmin = -1.0
    if not np.isfinite(vmax):
        vmax = 1.0
    if math.isclose(vmin, vmax):
        vmin -= 0.5
        vmax += 0.5

    fig, ax = plt.subplots(figsize=(1.5 * matrix.shape[1] + 3.5, 0.75 * matrix.shape[0] + 3))
    im = ax.imshow(values, cmap="viridis", aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.astype(str), rotation=20, ha="right")
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.astype(str))
    ax.set_title("Matched-tissue disease shift of panel scores\n(sample-level mean score)")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="white" if val < (vmin + vmax) / 2 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Sample mean panel score")
    fig.tight_layout()
    return save_figure(fig, figure_dir / "bcell_internal_disease_shift_heatmap")


def plot_sample_boxplots(sample_scores: pd.DataFrame, figure_dir: Path, subtype_order: list[str], max_panels_per_subtype: int) -> list[str]:
    if sample_scores.empty:
        return []
    disease_palette = {
        "healthy": "#4daf4a",
        "CRSsNP": "#377eb8",
        "CRSwNP": "#e41a1c",
        "asthma": "#984ea3",
        "COVID-19": "#ff7f00",
    }
    panel_order = (
        sample_scores[["panel_key", "subtype", "site_anchor"]]
        .drop_duplicates()
        .sort_values(["subtype", "site_anchor"])
    )
    selected_panels = []
    for subtype in subtype_order:
        sub = panel_order[panel_order["subtype"].astype(str) == subtype].head(max_panels_per_subtype)
        selected_panels.extend(sub["panel_key"].astype(str).tolist())
    selected_panels.extend([x for x in panel_order["panel_key"].astype(str).tolist() if x not in selected_panels])
    panel_order = panel_order[panel_order["panel_key"].astype(str).isin(selected_panels)].reset_index(drop=True)

    n_panels = panel_order.shape[0]
    n_cols = 2 if n_panels > 1 else 1
    n_rows = math.ceil(n_panels / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.2 * n_cols, 4.0 * n_rows), squeeze=False)
    axes_flat = axes.flatten()

    for ax, row in zip(axes_flat, panel_order.itertuples(index=False)):
        panel_df = sample_scores[sample_scores["panel_key"].astype(str) == str(row.panel_key)].copy()
        tissues = panel_df["tissue"].astype(str).unique().tolist()
        tissues = sorted(tissues)
        tissue_to_x = {t: i for i, t in enumerate(tissues)}
        box_data = [panel_df.loc[panel_df["tissue"].astype(str) == tissue, "score_mean"].astype(float).to_numpy() for tissue in tissues]
        ax.boxplot(box_data, positions=list(tissue_to_x.values()), widths=0.5, patch_artist=True, boxprops={"facecolor": "#d9d9d9", "edgecolor": "black"}, medianprops={"color": "black"})
        for _, r in panel_df.iterrows():
            x = tissue_to_x[str(r["tissue"])] + np.random.uniform(-0.12, 0.12)
            color = disease_palette.get(str(r["disease_level_1"]), "#555555")
            ax.scatter(x, float(r["score_mean"]), color=color, s=28, alpha=0.85, edgecolor="black", linewidth=0.3)
        ax.set_xticks(list(tissue_to_x.values()))
        ax.set_xticklabels(tissues)
        ax.set_ylabel("Sample mean panel score")
        ax.set_title(f"{row.subtype} | anchor={row.site_anchor}")
        ax.grid(axis="y", linestyle="--", alpha=0.25)

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    handles = [
        plt.Line2D([0], [0], marker="o", color="white", markerfacecolor=color, markeredgecolor="black", markersize=7, linestyle="", label=label)
        for label, color in disease_palette.items()
        if label in set(sample_scores["disease_level_1"].astype(str))
    ]
    if handles:
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.01),
            ncol=min(5, len(handles)),
            frameon=False,
            title="disease_level_1",
        )
    fig.suptitle("Sample-level panel scores by tissue and disease group", y=1.06, fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return save_figure(fig, figure_dir / "bcell_internal_sample_score_boxplots")


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    run_root = resolve_run_root(cfg, args.run_root)
    figure_dir = run_root / "validation" / "figures"
    ensure_dirs([figure_dir])

    viz_cfg = deep_get(cfg, CONFIG_KEY, "viz", default={}) or {}
    subtype_order = [str(x) for x in deep_get(viz_cfg, "subtype_order", default=["Naive_B", "Memory_B", "Plasma"])]
    site_order = [str(x) for x in deep_get(viz_cfg, "site_order", default=["nose", "sinus"])]
    max_panels_per_subtype = int(deep_get(viz_cfg, "max_panels_per_subtype", default=8))

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "figure.facecolor": "white",
    })

    tables = load_tables(run_root)
    manifest: list[dict[str, object]] = []

    manifest.append({
        "name": "label_transfer_qc",
        "outputs": plot_label_qc(tables["labels"], figure_dir, subtype_order=subtype_order),
        "sources": [str(run_root / "labels" / "query_label_transfer.tsv")],
    })
    manifest.append({
        "name": "retention_heatmap",
        "outputs": plot_retention_heatmap(tables["retention"], figure_dir, subtype_order=subtype_order, site_order=site_order),
        "sources": [str(run_root / "validation" / "retention_summary.tsv")],
    })
    manifest.append({
        "name": "disease_shift_heatmap",
        "outputs": plot_disease_shift_heatmap(tables["disease_shift"], figure_dir, subtype_order=subtype_order),
        "sources": [str(run_root / "validation" / "disease_shift_summary.tsv")],
    })
    manifest.append({
        "name": "sample_score_boxplots",
        "outputs": plot_sample_boxplots(tables["sample_scores"], figure_dir, subtype_order=subtype_order, max_panels_per_subtype=max_panels_per_subtype),
        "sources": [str(run_root / "validation" / "sample_scores.tsv")],
    })

    manifest = [item for item in manifest if item["outputs"]]
    write_json({"run_root": str(run_root), "figures": manifest}, figure_dir / "figure_manifest.json")
    print(json.dumps({"n_figures": len(manifest)}, ensure_ascii=False, indent=2))
    for item in manifest:
        print(f"[ok] {item['name']} -> {', '.join(item['outputs'])}")


if __name__ == "__main__":
    main()
