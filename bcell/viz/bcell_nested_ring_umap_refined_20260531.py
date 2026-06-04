#!/usr/bin/env python3
"""Refined B-cell nested-ring UMAP preview.

Refinements over the initial prototype:
  - larger / thicker nested rings
  - visible white gaps between L2 and L3 cell-type sectors
  - outer UMAP boundary colored by L2 identity while still avoiding
    internal borders between adjacent cell groups
"""

from __future__ import annotations

import argparse
import colorsys
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import anndata as ad
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter

DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/"
    "bcell_tissue_comparison_final_fullgene.h5ad"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/h2048/data/py/20260531/bcell_nested_ring_umap_refined_20260531"
)

L2_COLOR_MAP: Dict[str, str] = {
    "GC_B": "#F8766D",
    "Memory_B": "#7CAE00",
    "Naive_B": "#00BFC4",
    "Plasma": "#C77CFF",
}

L3_COLOR_MAP: Dict[str, str] = {
    "Atypical_Memory_B": "#F8766D",
    "GC_B_Dark_Zone_Centroblast_Cycling": "#D89000",
    "GC_B_Light_Zone_Centrocyte": "#7CAE00",
    "GC_B_Transitional": "#00BA38",
    "Memory_B": "#00BFC4",
    "Naive_B": "#00A6FF",
    "Plasma_IgA": "#C77CFF",
    "Plasma_IgG": "#FF61CC",
}

TISSUE_COLOR_MAP: Dict[str, str] = {
    "lung parenchyma": "#D55E00",
    "nose": "#009E73",
    "respiratory airway": "#0072B2",
    "sinus": "#CC79A7",
}

L3_SHORT_LABELS: Dict[str, str] = {
    "Atypical_Memory_B": "AtypMem",
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC-DZ",
    "GC_B_Light_Zone_Centrocyte": "GC-LZ",
    "GC_B_Transitional": "GC-Trans",
    "Memory_B": "Memory",
    "Naive_B": "Naive",
    "Plasma_IgA": "Pl-IgA",
    "Plasma_IgG": "Pl-IgG",
}

L2_SHORT_LABELS: Dict[str, str] = {
    "GC_B": "GC",
    "Memory_B": "Memory",
    "Naive_B": "Naive",
    "Plasma": "Plasma",
}

UMAP_FALLBACKS: Tuple[str, ...] = (
    "X_umap_scanvi",
    "X_umap_scANVI",
    "X_umap",
    "X_umap_scvi",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create refined nested-ring B-cell UMAP preview.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--prefix", default="bcell_nested_ring_umap_refined_20260531")
    parser.add_argument("--l2-col", default="cell_type_L2")
    parser.add_argument("--l3-col", default="cell_type_L3")
    parser.add_argument("--tissue-col", default="tissue")
    parser.add_argument("--umap-key", default=None)
    parser.add_argument("--title", default="B-cell nested ring UMAP refined preview")
    return parser.parse_args()


def choose_umap_key(adata: ad.AnnData, requested: str | None) -> str:
    if requested:
        if requested not in adata.obsm_keys():
            raise KeyError(f"Requested UMAP key not found: {requested}")
        return requested
    for key in UMAP_FALLBACKS:
        if key in adata.obsm_keys():
            return key
    raise KeyError(f"No supported UMAP key found. Tried: {', '.join(UMAP_FALLBACKS)}")


def load_plot_df(
    input_h5ad: Path,
    l2_col: str,
    l3_col: str,
    tissue_col: str,
    umap_key: str | None,
) -> tuple[pd.DataFrame, str]:
    adata = ad.read_h5ad(input_h5ad)
    chosen_umap_key = choose_umap_key(adata, umap_key)
    coords = np.asarray(adata.obsm[chosen_umap_key])
    if coords.shape[1] < 2:
        raise ValueError(f"UMAP key {chosen_umap_key} does not have at least 2 dimensions")

    required_cols = [l2_col, l3_col, tissue_col]
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    obs = adata.obs[required_cols].copy()
    plot_df = pd.DataFrame(
        {
            "umap_1": coords[:, 0],
            "umap_2": coords[:, 1],
        },
        index=adata.obs_names.copy(),
    ).join(obs)

    plot_df = plot_df.dropna(subset=[l2_col, l3_col, tissue_col, "umap_1", "umap_2"]).copy()
    for col in [l2_col, l3_col, tissue_col]:
        plot_df[col] = plot_df[col].astype(str).str.strip()
        plot_df = plot_df.loc[plot_df[col].ne("")].copy()

    return plot_df, chosen_umap_key


def compute_group_centroids(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    centroids = (
        df.groupby(group_col, observed=False)[["umap_1", "umap_2"]]
        .mean()
        .rename(columns={"umap_1": "cx", "umap_2": "cy"})
        .reset_index()
    )
    return centroids


def compute_plot_center(df: pd.DataFrame) -> tuple[float, float]:
    return (
        0.5 * (df["umap_1"].min() + df["umap_1"].max()),
        0.5 * (df["umap_2"].min() + df["umap_2"].max()),
    )


def add_clockwise_order(
    centroids: pd.DataFrame,
    center_xy: tuple[float, float],
    group_col: str,
) -> pd.DataFrame:
    center_x, center_y = center_xy
    angle_deg = np.degrees(np.arctan2(centroids["cy"] - center_y, centroids["cx"] - center_x))
    angle_deg = np.mod(angle_deg, 360.0)
    clockwise_from_top = np.mod(90.0 - angle_deg, 360.0)
    centroids = centroids.copy()
    centroids["angle_deg"] = angle_deg
    centroids["clockwise_from_top"] = clockwise_from_top
    centroids = centroids.sort_values(["clockwise_from_top", group_col]).reset_index(drop=True)
    return centroids


def build_orderings(
    df: pd.DataFrame,
    l2_col: str,
    l3_col: str,
) -> tuple[list[str], dict[str, list[str]], pd.DataFrame, pd.DataFrame]:
    center_xy = compute_plot_center(df)
    l2_centroids = add_clockwise_order(compute_group_centroids(df, l2_col), center_xy, l2_col)
    l2_order = l2_centroids[l2_col].tolist()

    l3_order_by_l2: dict[str, list[str]] = {}
    l3_centroid_frames: list[pd.DataFrame] = []
    l3_parent_rows: list[pd.DataFrame] = []
    for l2_label in l2_order:
        sub_df = df.loc[df[l2_col].eq(l2_label)].copy()
        l3_centroids = add_clockwise_order(compute_group_centroids(sub_df, l3_col), center_xy, l3_col)
        l3_centroids["parent_l2"] = l2_label
        l3_order_by_l2[l2_label] = l3_centroids[l3_col].tolist()
        l3_centroid_frames.append(l3_centroids)
        l3_parent_rows.append(
            sub_df.groupby(l3_col, observed=False)[l2_col].first().rename("parent_l2").reset_index()
        )

    l3_centroids_all = pd.concat(l3_centroid_frames, ignore_index=True)
    l3_parent_df = pd.concat(l3_parent_rows, ignore_index=True).drop_duplicates(subset=[l3_col])
    return (
        l2_order,
        l3_order_by_l2,
        l2_centroids,
        l3_centroids_all.merge(l3_parent_df, on=l3_col, how="left", suffixes=("", "_resolved")),
    )


def composition_table(
    df: pd.DataFrame,
    group_col: str,
    tissue_col: str,
    group_order: Iterable[str],
    tissue_order: Iterable[str],
) -> pd.DataFrame:
    table = (
        df.groupby([group_col, tissue_col], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    table[group_col] = pd.Categorical(table[group_col], categories=list(group_order), ordered=True)
    table[tissue_col] = pd.Categorical(table[tissue_col], categories=list(tissue_order), ordered=True)
    table = table.sort_values([group_col, tissue_col]).reset_index(drop=True)
    totals = table.groupby(group_col, observed=False)["n_cells"].transform("sum")
    table["fraction_within_group"] = np.where(totals > 0, table["n_cells"] / totals, 0.0)
    table["percent_within_group"] = 100.0 * table["fraction_within_group"]
    return table


def build_sector_manifest(
    df: pd.DataFrame,
    l2_col: str,
    l3_col: str,
    tissue_col: str,
    l2_order: list[str],
    l3_order_by_l2: dict[str, list[str]],
    l2_centroids: pd.DataFrame,
    l3_centroids: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    width_l2 = 2 * np.pi / max(len(l2_order), 1)
    cursor = 0.0

    l2_centroids_lookup = l2_centroids.set_index(l2_col).to_dict(orient="index")
    l3_centroids_lookup = l3_centroids.set_index(l3_col).to_dict(orient="index")

    for l2_label in l2_order:
        start = cursor
        end = cursor + width_l2
        sub_df = df.loc[df[l2_col].eq(l2_label)].copy()
        rows.append(
            {
                "level": "L2",
                "label": l2_label,
                "parent": "",
                "start_rad": start,
                "end_rad": end,
                "mid_rad": 0.5 * (start + end),
                "n_cells": int(sub_df.shape[0]),
                "centroid_x": float(l2_centroids_lookup[l2_label]["cx"]),
                "centroid_y": float(l2_centroids_lookup[l2_label]["cy"]),
                "angle_deg": float(l2_centroids_lookup[l2_label]["angle_deg"]),
            }
        )

        l3_labels = l3_order_by_l2.get(l2_label, [])
        width_l3 = width_l2 / max(len(l3_labels), 1)
        l3_cursor = start
        for l3_label in l3_labels:
            l3_start = l3_cursor
            l3_end = l3_cursor + width_l3
            l3_df = sub_df.loc[sub_df[l3_col].eq(l3_label)].copy()
            rows.append(
                {
                    "level": "L3",
                    "label": l3_label,
                    "parent": l2_label,
                    "start_rad": l3_start,
                    "end_rad": l3_end,
                    "mid_rad": 0.5 * (l3_start + l3_end),
                    "n_cells": int(l3_df.shape[0]),
                    "centroid_x": float(l3_centroids_lookup[l3_label]["cx"]),
                    "centroid_y": float(l3_centroids_lookup[l3_label]["cy"]),
                    "angle_deg": float(l3_centroids_lookup[l3_label]["angle_deg"]),
                }
            )
            l3_cursor = l3_end

        cursor = end

    manifest = pd.DataFrame(rows)
    if not manifest.empty:
        manifest["start_deg_clockwise_from_top"] = np.degrees(manifest["start_rad"])
        manifest["end_deg_clockwise_from_top"] = np.degrees(manifest["end_rad"])
        manifest["mid_deg_clockwise_from_top"] = np.degrees(manifest["mid_rad"])
    return manifest


def adjust_lightness(color: str, amount: float = 1.0, alpha: float | None = None) -> tuple[float, float, float, float]:
    r, g, b = mcolors.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = np.clip(l * amount, 0.0, 1.0)
    rr, gg, bb = colorsys.hls_to_rgb(h, l, s)
    if alpha is None:
        alpha = 1.0
    return (rr, gg, bb, alpha)


def sector_span_with_gap(
    start: float,
    end: float,
    gap_fraction: float,
    max_gap_deg: float,
) -> tuple[float, float, float, float]:
    width = max(end - start, 0.0)
    if width <= 0:
        return start, end, 0.0, 0.0
    gap = min(width * gap_fraction, np.deg2rad(max_gap_deg))
    gap = min(gap, width * 0.34)
    draw_start = start + 0.5 * gap
    draw_end = end - 0.5 * gap
    return draw_start, draw_end, max(draw_end - draw_start, 0.0), gap


def add_group_colored_outer_boundary(
    ax: plt.Axes,
    df: pd.DataFrame,
    group_col: str,
    group_order: list[str],
    color_map: Dict[str, str],
    bins: int = 380,
    sigma: float = 3.0,
    level_fraction: float = 0.085,
) -> None:
    x = df["umap_1"].to_numpy()
    y = df["umap_2"].to_numpy()
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))
    hist_union, xedges, yedges = np.histogram2d(
        x,
        y,
        bins=bins,
        range=[[xmin, xmax], [ymin, ymax]],
    )
    smooth_union = gaussian_filter(hist_union, sigma=sigma)
    union_max = float(smooth_union.max())
    if union_max <= 0:
        return

    level = max(union_max * level_fraction, 0.015)
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])

    group_smooth: dict[str, np.ndarray] = {}
    for label in group_order:
        sub = df.loc[df[group_col].eq(label)]
        if sub.empty:
            continue
        hist, _, _ = np.histogram2d(
            sub["umap_1"].to_numpy(),
            sub["umap_2"].to_numpy(),
            bins=[xedges, yedges],
        )
        group_smooth[label] = gaussian_filter(hist, sigma=sigma)

    contour = ax.contour(
        xcenters,
        ycenters,
        smooth_union.T,
        levels=[level],
        colors="none",
        linewidths=0.0,
        alpha=0.0,
        zorder=4,
    )

    segments: list[np.ndarray] = []
    segment_colors: list[tuple[float, float, float, float]] = []
    group_labels = [label for label in group_order if label in group_smooth]
    if not group_labels:
        return

    for polyline in contour.allsegs[0]:
        if polyline is None or len(polyline) < 2:
            continue
        p0 = polyline[:-1]
        p1 = polyline[1:]
        mids = 0.5 * (p0 + p1)
        ix = np.clip(np.searchsorted(xedges, mids[:, 0], side="right") - 1, 0, bins - 1)
        iy = np.clip(np.searchsorted(yedges, mids[:, 1], side="right") - 1, 0, bins - 1)
        score_matrix = np.column_stack([group_smooth[label][ix, iy] for label in group_labels])
        best_idx = np.argmax(score_matrix, axis=1)
        segments.extend(np.stack([p0, p1], axis=1))
        segment_colors.extend([
            adjust_lightness(color_map.get(group_labels[idx], "#444444"), amount=0.92, alpha=1.0)
            for idx in best_idx
        ])

    if hasattr(contour, "collections"):
        for coll in contour.collections:
            coll.remove()
    elif hasattr(contour, "remove"):
        try:
            contour.remove()
        except Exception:
            pass

    if not segments:
        return

    underlay = LineCollection(segments, colors=[(1, 1, 1, 0.98)] * len(segments), linewidths=4.6, zorder=4.35)
    underlay.set_capstyle("round")
    underlay.set_joinstyle("round")
    ax.add_collection(underlay)

    overlay = LineCollection(segments, colors=segment_colors, linewidths=2.8, zorder=4.55)
    overlay.set_capstyle("round")
    overlay.set_joinstyle("round")
    ax.add_collection(overlay)


def add_sector_label(
    polar_ax: plt.Axes,
    theta: float,
    radius: float,
    text: str,
    color: str,
    fontsize: float,
    weight: str = "normal",
    alpha: float = 0.9,
) -> None:
    polar_ax.text(
        theta,
        radius,
        text,
        color=color,
        fontsize=fontsize,
        fontweight=weight,
        ha="center",
        va="center",
        clip_on=False,
        bbox=dict(
            boxstyle="round,pad=0.17",
            facecolor=(1, 1, 1, alpha),
            edgecolor="none",
        ),
        zorder=14,
    )


def draw_nested_rings(
    polar_ax: plt.Axes,
    df: pd.DataFrame,
    l2_col: str,
    l3_col: str,
    tissue_col: str,
    l2_order: list[str],
    l3_order_by_l2: dict[str, list[str]],
) -> pd.DataFrame:
    tissue_order = [t for t in TISSUE_COLOR_MAP if t in set(df[tissue_col])]

    inner_r0, inner_r1 = 0.405, 0.515
    outer_r0, outer_r1 = 0.535, 0.655
    l2_gap_fraction, l3_gap_fraction = 0.055, 0.075

    width_l2 = 2 * np.pi / max(len(l2_order), 1)
    cursor = 0.0
    sector_rows: list[dict[str, object]] = []

    for l2_label in l2_order:
        start = cursor
        end = cursor + width_l2
        l2_draw_start, l2_draw_end, l2_draw_width, _ = sector_span_with_gap(
            start,
            end,
            gap_fraction=l2_gap_fraction,
            max_gap_deg=4.4,
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
                color=TISSUE_COLOR_MAP[tissue],
                edgecolor="white",
                linewidth=0.65,
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
                edgecolor=adjust_lightness(L2_COLOR_MAP.get(l2_label, "#444444"), amount=0.72, alpha=0.95),
                linewidth=1.65,
                align="center",
                zorder=9,
            )
            polar_ax.bar(
                x=(l2_draw_start + l2_draw_end) / 2,
                height=0.012,
                width=l2_draw_width,
                bottom=inner_r1 - 0.012,
                color=adjust_lightness(L2_COLOR_MAP.get(l2_label, "#777777"), amount=0.96, alpha=0.95),
                edgecolor="none",
                align="center",
                zorder=10,
            )

        add_sector_label(
            polar_ax,
            theta=(start + end) / 2,
            radius=inner_r0 - 0.033,
            text=L2_SHORT_LABELS.get(l2_label, l2_label),
            color=L2_COLOR_MAP.get(l2_label, "#333333"),
            fontsize=10.0,
            weight="bold",
            alpha=0.93,
        )

        l3_labels = l3_order_by_l2.get(l2_label, [])
        width_l3 = width_l2 / max(len(l3_labels), 1)
        l3_cursor = start
        for l3_label in l3_labels:
            l3_start = l3_cursor
            l3_end = l3_cursor + width_l3
            l3_draw_start, l3_draw_end, l3_draw_width, _ = sector_span_with_gap(
                l3_start,
                l3_end,
                gap_fraction=l3_gap_fraction,
                max_gap_deg=3.2,
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
                polar_ax.bar(
                    x=sub_cursor + seg_width / 2,
                    height=outer_r1 - outer_r0,
                    width=seg_width,
                    bottom=outer_r0,
                    color=TISSUE_COLOR_MAP[tissue],
                    edgecolor="white",
                    linewidth=0.55,
                    align="center",
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
                    height=outer_r1 - outer_r0,
                    width=l3_draw_width,
                    bottom=outer_r0,
                    color="none",
                    edgecolor=adjust_lightness(L3_COLOR_MAP.get(l3_label, "#444444"), amount=0.76, alpha=0.95),
                    linewidth=1.15,
                    align="center",
                    zorder=11,
                )
                polar_ax.bar(
                    x=(l3_draw_start + l3_draw_end) / 2,
                    height=0.010,
                    width=l3_draw_width,
                    bottom=outer_r1 - 0.010,
                    color=adjust_lightness(L3_COLOR_MAP.get(l3_label, "#777777"), amount=0.96, alpha=0.92),
                    edgecolor="none",
                    align="center",
                    zorder=12,
                )

            add_sector_label(
                polar_ax,
                theta=(l3_start + l3_end) / 2,
                radius=outer_r1 + 0.035,
                text=L3_SHORT_LABELS.get(l3_label, l3_label),
                color=L3_COLOR_MAP.get(l3_label, "#333333"),
                fontsize=7.9,
                weight="normal",
                alpha=0.88,
            )
            l3_cursor = l3_end

        cursor = end

    polar_ax.set_theta_offset(np.pi / 2.0)
    polar_ax.set_theta_direction(-1)
    polar_ax.set_ylim(0, 0.735)
    polar_ax.set_axis_off()

    sector_df = pd.DataFrame(sector_rows)
    if not sector_df.empty:
        sector_df["start_deg_clockwise_from_top"] = np.degrees(sector_df["start_rad"])
        sector_df["end_deg_clockwise_from_top"] = np.degrees(sector_df["end_rad"])
    return sector_df


def annotate_l2_centroids(ax: plt.Axes, df: pd.DataFrame, l2_col: str) -> None:
    centroids = compute_group_centroids(df, l2_col)
    for row in centroids.itertuples(index=False):
        label = getattr(row, l2_col)
        text = ax.text(
            row.cx,
            row.cy,
            label,
            fontsize=12.3,
            fontweight="bold",
            color="#222222",
            ha="center",
            va="center",
            zorder=6,
        )
        text.set_path_effects([pe.Stroke(linewidth=3.4, foreground="white"), pe.Normal()])


def style_main_axes(ax: plt.Axes, df: pd.DataFrame, title: str, umap_key: str) -> None:
    x = df["umap_1"].to_numpy()
    y = df["umap_2"].to_numpy()
    xr = float(np.max(x) - np.min(x))
    yr = float(np.max(y) - np.min(y))
    pad = 0.20 * max(xr, yr)
    ax.set_xlim(float(np.min(x) - pad), float(np.max(x) + pad))
    ax.set_ylim(float(np.min(y) - pad), float(np.max(y) + pad))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=16.5, fontweight="bold", pad=16)
    ax.text(
        0.5,
        1.01,
        f"center = {umap_key} with L2-colored outer hull | inner / outer rings = tissue composition within L2 / L3",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.2,
        color="#555555",
    )


def make_figure(
    df: pd.DataFrame,
    input_h5ad: Path,
    umap_key: str,
    l2_col: str,
    l3_col: str,
    tissue_col: str,
    output_dir: Path,
    prefix: str,
    title: str,
) -> dict[str, object]:
    l2_order, l3_order_by_l2, l2_centroids, l3_centroids = build_orderings(df, l2_col, l3_col)
    tissue_order = [t for t in TISSUE_COLOR_MAP if t in set(df[tissue_col])]

    l2_table = composition_table(df, l2_col, tissue_col, l2_order, tissue_order)
    l3_order_flat = [label for parent in l2_order for label in l3_order_by_l2.get(parent, [])]
    l3_table = composition_table(df, l3_col, tissue_col, l3_order_flat, tissue_order)
    sector_manifest = build_sector_manifest(df, l2_col, l3_col, tissue_col, l2_order, l3_order_by_l2, l2_centroids, l3_centroids)

    output_dir.mkdir(parents=True, exist_ok=True)
    l2_table.to_csv(output_dir / f"{prefix}_l2_tissue.tsv", sep="\t", index=False)
    l3_table.to_csv(output_dir / f"{prefix}_l3_tissue.tsv", sep="\t", index=False)
    sector_manifest.to_csv(output_dir / f"{prefix}_group_sector_manifest.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.0, 11.5), facecolor="white")
    ax = fig.add_axes([0.035, 0.065, 0.795, 0.865])

    for label in l2_order:
        sub = df.loc[df[l2_col].eq(label)]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=8.2,
            c=L2_COLOR_MAP.get(label, "#999999"),
            linewidths=0.0,
            alpha=0.94,
            rasterized=True,
            zorder=3,
        )

    add_group_colored_outer_boundary(ax, df, group_col=l2_col, group_order=l2_order, color_map=L2_COLOR_MAP)
    annotate_l2_centroids(ax, df, l2_col)
    style_main_axes(ax, df, title, umap_key)

    ring_ax = fig.add_axes(ax.get_position(), projection="polar", facecolor="none")
    ring_sector_df = draw_nested_rings(ring_ax, df, l2_col, l3_col, tissue_col, l2_order, l3_order_by_l2)
    if not ring_sector_df.empty:
        ring_sector_df.to_csv(output_dir / f"{prefix}_ring_segments.tsv", sep="\t", index=False)

    tissue_handles = [
        Patch(facecolor=color, edgecolor="none", label=label)
        for label, color in TISSUE_COLOR_MAP.items()
        if label in tissue_order
    ]
    fig.legend(
        handles=tissue_handles,
        loc="center right",
        bbox_to_anchor=(0.985, 0.50),
        frameon=False,
        title="Tissue",
        title_fontsize=12,
        fontsize=10.5,
    )

    fig.text(
        0.83,
        0.102,
        "Outer-ring abbreviations:\nAtypMem, GC-DZ, GC-LZ, GC-Trans,\nMemory, Naive, Pl-IgA, Pl-IgG",
        ha="left",
        va="bottom",
        fontsize=8.3,
        color="#555555",
    )

    png_path = output_dir / f"{prefix}.png"
    pdf_path = output_dir / f"{prefix}.pdf"
    fig.savefig(png_path, dpi=340, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, dpi=340, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    summary = {
        "input_h5ad": str(input_h5ad),
        "output_dir": str(output_dir),
        "png": str(png_path),
        "pdf": str(pdf_path),
        "n_obs": int(df.shape[0]),
        "umap_key": umap_key,
        "l2_order": l2_order,
        "l3_order_by_l2": l3_order_by_l2,
        "tissue_order": tissue_order,
        "l2_levels": sorted(df[l2_col].astype(str).unique().tolist()),
        "l3_levels": sorted(df[l3_col].astype(str).unique().tolist()),
        "notes": [
            "UMAP outer hull is colored by L2 identity while internal inter-group borders remain unoutlined.",
            "Both L2 and L3 sectors use fixed-width cell-type slots with intentional white gaps between sectors.",
            "Ring bar thickness is enlarged relative to the initial prototype.",
        ],
    }
    (output_dir / f"{prefix}_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    args = parse_args()
    plot_df, umap_key = load_plot_df(
        input_h5ad=args.input_h5ad,
        l2_col=args.l2_col,
        l3_col=args.l3_col,
        tissue_col=args.tissue_col,
        umap_key=args.umap_key,
    )
    summary = make_figure(
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
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
