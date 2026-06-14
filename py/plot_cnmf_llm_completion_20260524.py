#!/usr/bin/env python3
"""Plot completion overview for cNMF LLM manifests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "lineage",
    "celltype_value",
    "total_units",
    "ok_units",
    "error_units",
}
LEVEL_ORDER = ["L2", "L3"]
PANEL_LABELS = ["A", "B", "C", "D"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate status plots for L2/L3 cNMF LLM completion manifests."
    )
    parser.add_argument(
        "--l2-manifest",
        default="/home/h2048/output/program_full_parallel_methods_20260507/cnmf_l2_celltype_LLM_combined_manifest.tsv",
        help="Path to the L2 combined manifest TSV.",
    )
    parser.add_argument(
        "--l3-manifest",
        default="/home/h2048/output/program_full_parallel_methods_20260507/cnmf_l3_celltype_LLM_combined_manifest.tsv",
        help="Path to the L3 combined manifest TSV.",
    )
    parser.add_argument(
        "--output-prefix",
        default="/home/h2048/output/program_full_parallel_methods_20260507/figures/cnmf_llm_completion_20260524",
        help="Output prefix for PNG/SVG and summary tables.",
    )
    parser.add_argument(
        "--title",
        default="cNMF LLM completion overview",
        help="Main title for the figure.",
    )
    return parser.parse_args()


def load_manifest(path: str | Path, level: str) -> pd.DataFrame:
    manifest_path = Path(path)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    df = pd.read_csv(manifest_path, sep="\t")
    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(
            f"Manifest {manifest_path} missing required columns: {sorted(missing)}"
        )

    out = df.copy()
    out["level"] = level
    out["total_units"] = pd.to_numeric(out["total_units"], errors="raise")
    out["ok_units"] = pd.to_numeric(out["ok_units"], errors="raise")
    out["error_units"] = pd.to_numeric(out["error_units"], errors="raise")
    out["completion_pct"] = np.where(
        out["total_units"] > 0,
        out["ok_units"] / out["total_units"] * 100.0,
        np.nan,
    )
    return out


def summarise(combined: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_level = (
        combined.groupby("level", as_index=False)
        .agg(
            celltype_groups=("celltype_value", "size"),
            total_units=("total_units", "sum"),
            ok_units=("ok_units", "sum"),
            error_units=("error_units", "sum"),
        )
        .assign(
            completion_pct=lambda d: np.where(
                d["total_units"] > 0,
                d["ok_units"] / d["total_units"] * 100.0,
                np.nan,
            )
        )
    )

    by_lineage_level = (
        combined.groupby(["lineage", "level"], as_index=False)
        .agg(
            celltype_groups=("celltype_value", "size"),
            total_units=("total_units", "sum"),
            ok_units=("ok_units", "sum"),
            error_units=("error_units", "sum"),
        )
        .assign(
            completion_pct=lambda d: np.where(
                d["total_units"] > 0,
                d["ok_units"] / d["total_units"] * 100.0,
                np.nan,
            )
        )
    )

    level_rank = {level: idx for idx, level in enumerate(LEVEL_ORDER)}
    by_level = by_level.sort_values("level", key=lambda s: s.map(level_rank)).reset_index(drop=True)
    by_lineage_level = by_lineage_level.sort_values(
        ["lineage", "level"], key=lambda s: s.map(level_rank) if s.name == "level" else s
    ).reset_index(drop=True)
    return by_level, by_lineage_level


def draw_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def plot_completion_overview(
    by_level: pd.DataFrame,
    by_lineage_level: pd.DataFrame,
    output_prefix: Path,
    title: str,
) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)
    fig.suptitle(title, fontsize=18, fontweight="bold")

    ok_color = "#2E8B57"
    err_color = "#D55E00"
    level_colors = {"L2": "#4C78A8", "L3": "#8E63CE"}

    # Panel A: overall completion by level.
    ax = axes[0, 0]
    y = np.arange(len(by_level))
    ax.barh(y, by_level["ok_units"], color=ok_color, label="OK units")
    ax.barh(
        y,
        by_level["error_units"],
        left=by_level["ok_units"],
        color=err_color,
        label="Error units",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(by_level["level"])
    ax.set_xlabel("Unit count")
    ax.set_title("Overall completion by level")
    max_units = max(float(by_level["total_units"].max()), 1.0)
    ax.set_xlim(0, max_units * 1.18)
    for idx, row in by_level.iterrows():
        label = (
            f"{int(row['ok_units'])}/{int(row['total_units'])} OK  "
            f"({row['celltype_groups']} celltype groups, {row['completion_pct']:.1f}%)"
        )
        ax.text(
            row["total_units"] + max_units * 0.02,
            idx,
            label,
            va="center",
            ha="left",
            fontsize=10,
        )
    ax.legend(frameon=False, loc="lower right")
    draw_panel_label(ax, PANEL_LABELS[0])

    # Panel B: total units by lineage and level.
    ax = axes[0, 1]
    pivot_units = (
        by_lineage_level.pivot(index="lineage", columns="level", values="total_units")
        .fillna(0)
        .reindex(columns=LEVEL_ORDER)
    )
    x = np.arange(len(pivot_units.index))
    width = 0.36
    for offset, level in enumerate(LEVEL_ORDER):
        values = pivot_units[level].to_numpy()
        ax.bar(
            x + (offset - 0.5) * width,
            values,
            width=width,
            label=level,
            color=level_colors[level],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(pivot_units.index, rotation=35, ha="right")
    ax.set_ylabel("Unit count")
    ax.set_title("Units covered per lineage")
    ax.legend(frameon=False)
    draw_panel_label(ax, PANEL_LABELS[1])

    # Panel C: celltype-group counts by lineage and level.
    ax = axes[1, 0]
    pivot_groups = (
        by_lineage_level.pivot(index="lineage", columns="level", values="celltype_groups")
        .fillna(0)
        .reindex(columns=LEVEL_ORDER)
    )
    for offset, level in enumerate(LEVEL_ORDER):
        values = pivot_groups[level].to_numpy()
        ax.bar(
            x + (offset - 0.5) * width,
            values,
            width=width,
            label=level,
            color=level_colors[level],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(pivot_groups.index, rotation=35, ha="right")
    ax.set_ylabel("Celltype-group count")
    ax.set_title("Celltype groups with combined markdown")
    ax.legend(frameon=False)
    draw_panel_label(ax, PANEL_LABELS[2])

    # Panel D: completion heatmap.
    ax = axes[1, 1]
    pivot_pct = (
        by_lineage_level.pivot(index="lineage", columns="level", values="completion_pct")
        .reindex(columns=LEVEL_ORDER)
    )
    im = ax.imshow(pivot_pct.to_numpy(), aspect="auto", cmap="Greens", vmin=0, vmax=100)
    ax.set_xticks(np.arange(len(pivot_pct.columns)))
    ax.set_xticklabels(pivot_pct.columns)
    ax.set_yticks(np.arange(len(pivot_pct.index)))
    ax.set_yticklabels(pivot_pct.index)
    ax.set_title("Completion rate by lineage × level")
    for i in range(pivot_pct.shape[0]):
        for j in range(pivot_pct.shape[1]):
            value = pivot_pct.iloc[i, j]
            text = "NA" if pd.isna(value) else f"{value:.1f}%"
            ax.text(j, i, text, ha="center", va="center", color="black", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Completion (%)")
    draw_panel_label(ax, PANEL_LABELS[3])

    fig.savefig(output_prefix.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_prefix = Path(args.output_prefix)

    l2 = load_manifest(args.l2_manifest, level="L2")
    l3 = load_manifest(args.l3_manifest, level="L3")
    combined = pd.concat([l2, l3], ignore_index=True)

    by_level, by_lineage_level = summarise(combined)
    plot_completion_overview(by_level, by_lineage_level, output_prefix, args.title)

    by_level.to_csv(output_prefix.with_suffix(".by_level.tsv"), sep="\t", index=False)
    by_lineage_level.to_csv(
        output_prefix.with_suffix(".by_lineage_level.tsv"), sep="\t", index=False
    )

    metadata = {
        "title": args.title,
        "l2_manifest": str(Path(args.l2_manifest).resolve()),
        "l3_manifest": str(Path(args.l3_manifest).resolve()),
        "png": str(output_prefix.with_suffix(".png")),
        "svg": str(output_prefix.with_suffix(".svg")),
        "by_level_tsv": str(output_prefix.with_suffix(".by_level.tsv")),
        "by_lineage_level_tsv": str(output_prefix.with_suffix(".by_lineage_level.tsv")),
        "levels": by_level.to_dict(orient="records"),
    }
    output_prefix.with_suffix(".meta.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"[OK] Wrote {output_prefix.with_suffix('.png')}")
    print(f"[OK] Wrote {output_prefix.with_suffix('.svg')}")
    print(f"[OK] Wrote {output_prefix.with_suffix('.by_level.tsv')}")
    print(f"[OK] Wrote {output_prefix.with_suffix('.by_lineage_level.tsv')}")


if __name__ == "__main__":
    main()
