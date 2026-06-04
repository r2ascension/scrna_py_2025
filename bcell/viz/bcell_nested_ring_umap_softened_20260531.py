#!/usr/bin/env python3
"""Softened B-cell nested-ring UMAP preview.

This version keeps the refined layout logic but pushes the rings much farther
outward, softens the outer ring edge with feathered radial fades, and switches
all palettes to gentler pastel tones.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd

BASE_SCRIPT_PATH = Path(__file__).with_name("bcell_nested_ring_umap_refined_20260531.py")
DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/"
    "bcell_tissue_comparison_final_fullgene.h5ad"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/h2048/data/py/20260531/bcell_nested_ring_umap_softened_20260531"
)

SOFT_L2_COLOR_MAP: Dict[str, str] = {
    "GC_B": "#F2A8A0",
    "Memory_B": "#B7CE84",
    "Naive_B": "#7FD2D9",
    "Plasma": "#D9B8F6",
}

SOFT_L3_COLOR_MAP: Dict[str, str] = {
    "Atypical_Memory_B": "#F5B4AE",
    "GC_B_Dark_Zone_Centroblast_Cycling": "#E9BF8A",
    "GC_B_Light_Zone_Centrocyte": "#B9D68A",
    "GC_B_Transitional": "#98D4A5",
    "Memory_B": "#8BD3DD",
    "Naive_B": "#94C2F2",
    "Plasma_IgA": "#D7B8F5",
    "Plasma_IgG": "#E8B8DD",
}

SOFT_TISSUE_COLOR_MAP: Dict[str, str] = {
    "lung parenchyma": "#E4A07C",
    "nose": "#87C7AA",
    "respiratory airway": "#8FAFE6",
    "sinus": "#D7A6D8",
}


def load_base_module():
    spec = importlib.util.spec_from_file_location("bcell_nested_ring_umap_refined_20260531", BASE_SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load base script from {BASE_SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create softened B-cell nested-ring UMAP preview.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--prefix", default="bcell_nested_ring_umap_softened_20260531")
    parser.add_argument("--l2-col", default="cell_type_L2")
    parser.add_argument("--l3-col", default="cell_type_L3")
    parser.add_argument("--tissue-col", default="tissue")
    parser.add_argument("--umap-key", default=None)
    parser.add_argument("--title", default="B-cell nested ring UMAP softened preview")
    return parser.parse_args()


def patch_soft_palettes(base) -> None:
    base.L2_COLOR_MAP = SOFT_L2_COLOR_MAP.copy()
    base.L3_COLOR_MAP = SOFT_L3_COLOR_MAP.copy()
    base.TISSUE_COLOR_MAP = SOFT_TISSUE_COLOR_MAP.copy()


def softened_style_main_axes(ax, df, title: str, umap_key: str) -> None:
    x = df["umap_1"].to_numpy()
    y = df["umap_2"].to_numpy()
    xr = float(x.max() - x.min())
    yr = float(y.max() - y.min())
    pad = 0.58 * max(xr, yr)
    ax.set_xlim(float(x.min() - pad), float(x.max() + pad))
    ax.set_ylim(float(y.min() - pad), float(y.max() + pad))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=16.0, fontweight="bold", pad=16)
    ax.text(
        0.5,
        1.01,
        f"center = {umap_key} with L2-colored outer hull | enlarged soft rings = tissue composition within L2 / L3",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.0,
        color="#666666",
    )


def draw_feathered_outer_segment(base, polar_ax, theta: float, width: float, bottom: float, core_height: float, feather_height: float, color: str, zorder: float) -> None:
    polar_ax.bar(
        x=theta,
        height=core_height,
        width=width,
        bottom=bottom,
        color=base.adjust_lightness(color, amount=1.01, alpha=0.90),
        edgecolor="none",
        linewidth=0.0,
        align="center",
        zorder=zorder,
    )
    if feather_height <= 0:
        return
    n_steps = 4
    step = feather_height / n_steps
    alphas = [0.22, 0.14, 0.08, 0.04]
    for i, alpha in enumerate(alphas):
        polar_ax.bar(
            x=theta,
            height=step,
            width=width * (1.002 + 0.002 * i),
            bottom=bottom + core_height + i * step,
            color=base.adjust_lightness(color, amount=1.06 + 0.03 * i, alpha=alpha),
            edgecolor="none",
            linewidth=0.0,
            align="center",
            zorder=zorder + 0.15 + i * 0.01,
        )


def softened_draw_nested_rings(polar_ax, df, l2_col: str, l3_col: str, tissue_col: str, l2_order, l3_order_by_l2):
    tissue_order = [t for t in SOFT_TISSUE_COLOR_MAP if t in set(df[tissue_col])]

    inner_r0, inner_r1 = 0.76, 0.94
    outer_core_r0, outer_core_r1 = 1.01, 1.16
    outer_feather = 0.12
    l2_gap_fraction, l3_gap_fraction = 0.030, 0.045

    width_l2 = 2 * base_module.np.pi / max(len(l2_order), 1)
    cursor = 0.0
    sector_rows = []

    for l2_label in l2_order:
        start = cursor
        end = cursor + width_l2
        l2_draw_start, l2_draw_end, l2_draw_width, _ = base_module.sector_span_with_gap(
            start,
            end,
            gap_fraction=l2_gap_fraction,
            max_gap_deg=3.2,
        )
        sector_df = df.loc[df[l2_col].eq(l2_label)].copy()
        tissue_counts = (
            sector_df.groupby(tissue_col, observed=False)
            .size()
            .reindex(tissue_order, fill_value=0)
        )
        total = int(tissue_counts.sum())
        tissue_cursor = l2_draw_start
        for tissue in tissue_order:
            frac = 0.0 if total == 0 else float(tissue_counts[tissue]) / float(total)
            seg_width = l2_draw_width * frac
            if seg_width <= 0:
                continue
            polar_ax.bar(
                x=tissue_cursor + seg_width / 2,
                height=inner_r1 - inner_r0,
                width=seg_width,
                bottom=inner_r0,
                color=base_module.adjust_lightness(SOFT_TISSUE_COLOR_MAP[tissue], amount=1.01, alpha=0.88),
                edgecolor="none",
                linewidth=0.0,
                align="center",
                zorder=6,
            )
            sector_rows.append(
                {
                    "level": "L2",
                    "label": l2_label,
                    "parent": "",
                    "tissue": tissue,
                    "start_rad": tissue_cursor,
                    "end_rad": tissue_cursor + seg_width,
                    "fraction_within_group": frac,
                    "n_cells": int(tissue_counts[tissue]),
                }
            )
            tissue_cursor += seg_width

        if l2_draw_width > 0:
            polar_ax.bar(
                x=(l2_draw_start + l2_draw_end) / 2,
                height=inner_r1 - inner_r0,
                width=l2_draw_width,
                bottom=inner_r0,
                color="none",
                edgecolor=base_module.adjust_lightness(SOFT_L2_COLOR_MAP.get(l2_label, "#777777"), amount=0.82, alpha=0.92),
                linewidth=1.45,
                align="center",
                zorder=9,
            )
            polar_ax.bar(
                x=(l2_draw_start + l2_draw_end) / 2,
                height=0.013,
                width=l2_draw_width,
                bottom=inner_r1 - 0.013,
                color=base_module.adjust_lightness(SOFT_L2_COLOR_MAP.get(l2_label, "#777777"), amount=1.03, alpha=0.90),
                edgecolor="none",
                align="center",
                zorder=10,
            )

        base_module.add_sector_label(
            polar_ax,
            theta=(start + end) / 2,
            radius=inner_r0 - 0.036,
            text=base_module.L2_SHORT_LABELS.get(l2_label, l2_label),
            color=SOFT_L2_COLOR_MAP.get(l2_label, "#555555"),
            fontsize=9.9,
            weight="bold",
            alpha=0.90,
        )

        l3_labels = l3_order_by_l2.get(l2_label, [])
        width_l3 = width_l2 / max(len(l3_labels), 1)
        l3_cursor = start
        for l3_label in l3_labels:
            l3_start = l3_cursor
            l3_end = l3_cursor + width_l3
            l3_draw_start, l3_draw_end, l3_draw_width, _ = base_module.sector_span_with_gap(
                l3_start,
                l3_end,
                gap_fraction=l3_gap_fraction,
                max_gap_deg=2.4,
            )
            l3_df = sector_df.loc[sector_df[l3_col].eq(l3_label)].copy()
            l3_tissue_counts = (
                l3_df.groupby(tissue_col, observed=False)
                .size()
                .reindex(tissue_order, fill_value=0)
            )
            l3_total = int(l3_tissue_counts.sum())
            sub_cursor = l3_draw_start
            for tissue in tissue_order:
                frac = 0.0 if l3_total == 0 else float(l3_tissue_counts[tissue]) / float(l3_total)
                seg_width = l3_draw_width * frac
                if seg_width <= 0:
                    continue
                draw_feathered_outer_segment(
                    base_module,
                    polar_ax,
                    theta=sub_cursor + seg_width / 2,
                    width=seg_width,
                    bottom=outer_core_r0,
                    core_height=outer_core_r1 - outer_core_r0,
                    feather_height=outer_feather,
                    color=SOFT_TISSUE_COLOR_MAP[tissue],
                    zorder=7,
                )
                sector_rows.append(
                    {
                        "level": "L3",
                        "label": l3_label,
                        "parent": l2_label,
                        "tissue": tissue,
                        "start_rad": sub_cursor,
                        "end_rad": sub_cursor + seg_width,
                        "fraction_within_group": frac,
                        "n_cells": int(l3_tissue_counts[tissue]),
                    }
                )
                sub_cursor += seg_width

            if l3_draw_width > 0:
                polar_ax.bar(
                    x=(l3_draw_start + l3_draw_end) / 2,
                    height=(outer_core_r1 - outer_core_r0) + outer_feather * 0.55,
                    width=l3_draw_width,
                    bottom=outer_core_r0,
                    color="none",
                    edgecolor=base_module.adjust_lightness(SOFT_L3_COLOR_MAP.get(l3_label, "#777777"), amount=0.84, alpha=0.85),
                    linewidth=1.05,
                    align="center",
                    zorder=12,
                )

            base_module.add_sector_label(
                polar_ax,
                theta=(l3_start + l3_end) / 2,
                radius=outer_core_r1 + outer_feather + 0.035,
                text=base_module.L3_SHORT_LABELS.get(l3_label, l3_label),
                color=SOFT_L3_COLOR_MAP.get(l3_label, "#555555"),
                fontsize=7.8,
                weight="normal",
                alpha=0.86,
            )
            l3_cursor = l3_end

        cursor = end

    polar_ax.set_theta_offset(base_module.np.pi / 2.0)
    polar_ax.set_theta_direction(-1)
    polar_ax.set_ylim(0, 1.34)
    polar_ax.set_axis_off()

    sector_df = pd.DataFrame(sector_rows)
    if not sector_df.empty:
        sector_df["start_deg_clockwise_from_top"] = base_module.np.degrees(sector_df["start_rad"])
        sector_df["end_deg_clockwise_from_top"] = base_module.np.degrees(sector_df["end_rad"])
    return sector_df


def softened_annotate_l2_centroids(ax, df, l2_col: str) -> None:
    centroids = base_module.compute_group_centroids(df, l2_col)
    for row in centroids.itertuples(index=False):
        label = getattr(row, l2_col)
        text = ax.text(
            row.cx,
            row.cy,
            label,
            fontsize=11.2,
            fontweight="bold",
            color="#2C2C2C",
            ha="center",
            va="center",
            zorder=6,
        )
        text.set_path_effects([
            base_module.pe.Stroke(linewidth=3.1, foreground="white"),
            base_module.pe.Normal(),
        ])


def main() -> None:
    global base_module
    base_module = load_base_module()
    patch_soft_palettes(base_module)
    base_module.style_main_axes = softened_style_main_axes
    base_module.draw_nested_rings = softened_draw_nested_rings
    base_module.annotate_l2_centroids = softened_annotate_l2_centroids

    args = parse_args()
    plot_df, umap_key = base_module.load_plot_df(
        input_h5ad=args.input_h5ad,
        l2_col=args.l2_col,
        l3_col=args.l3_col,
        tissue_col=args.tissue_col,
        umap_key=args.umap_key,
    )
    summary = base_module.make_figure(
        df=plot_df,
        input_h5ad=args.input_h5ad,
        umap_key=umap_key,
        l2_col=args.l2_col,
        l3_col=args.l3_col,
        tissue_col=args.tissue_col,
        output_dir=args.output_dir,
        prefix=args.prefix,
        title=args.title,
    )
    summary["notes"] = [
        "Outer ring radii are pushed to >=1.5x the previous refined outer radius and the UMAP center is compressed inward with extra axis padding.",
        "Outer ring segments use feathered radial fades instead of hard outer edges.",
        "All palettes are softened to pastel tones while preserving category identity.",
    ]
    summary_path = args.output_dir / f"{args.prefix}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
