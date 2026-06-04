#!/usr/bin/env python3
"""B-cell nested-ring UMAP preview.

Center:
  - UMAP colored by `cell_type_L2`
  - only global outer contour lines are drawn, so internal borders between
    adjacent cell groups are not outlined
Inner ring:
  - tissue composition within each L2 group
Outer ring:
  - tissue composition within each L3 group, nested inside the parent L2 sector

Outputs:
  - PNG / PDF figure
  - TSV tables for L2 and L3 tissue composition
  - sector manifest for reproducibility / future refinements
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import anndata as ad
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter

DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/"
    "bcell_tissue_comparison_final_fullgene.h5ad"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/h2048/data/py/20260531/bcell_nested_ring_umap_preview_20260531"
)

# ggplot2 default discrete palette values used by current lineage plots
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
    parser = argparse.ArgumentParser(description="Create nested-ring B-cell UMAP preview.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--prefix", default="bcell_nested_ring_umap_20260531")
    parser.add_argument("--l2-col", default="cell_type_L2")
    parser.add_argument("--l3-col", default="cell_type_L3")
    parser.add_argument("--tissue-col", default="tissue")
    parser.add_argument("--umap-key", default=None)
    parser.add_argument("--title", default="B-cell nested ring UMAP preview")
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

    l3_to_l2 = (
        plot_df.groupby(l3_col)[l2_col]
        .agg(lambda x: x.mode().iat[0])
        .to_dict()
    )
    unexpected_l3 = [label for label in plot_df[l3_col].unique() if label not in l3_to_l2]
    if unexpected_l3:
        raise ValueError(f"Could not resolve L3->L2 mapping for: {unexpected_l3}")

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
    return l2_order, l3_order_by_l2, l2_centroids, l3_centroids_all.merge(l3_parent_df, on=l3_col, how="left", suffixes=("", "_resolved"))


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


def add_global_contour(ax: plt.Axes, x: np.ndarray, y: np.ndarray) -> None:
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    bins = 360
    hist, xedges, yedges = np.histogram2d(
        x,
        y,
        bins=bins,
        range=[[xmin, xmax], [ymin, ymax]],
    )
    smooth = gaussian_filter(hist, sigma=3.2)
    if float(smooth.max()) <= 0:
        return
    level = max(float(smooth.max()) * 0.085, 0.015)
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    ax.contour(
        xcenters,
        ycenters,
        smooth.T,
        levels=[level],
        colors=["#2F2F2F"],
        linewidths=1.3,
        alpha=0.95,
        zorder=4,
    )


def add_sector_label(
    polar_ax: plt.Axes,
    theta: float,
    radius: float,
    text: str,
    color: str,
    fontsize: float,
    weight: str = "normal",
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
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.82),
        zorder=10,
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
    inner_r0, inner_r1 = 0.38, 0.445
    outer_r0, outer_r1 = 0.445, 0.50

    width_l2 = 2 * np.pi / max(len(l2_order), 1)
    cursor = 0.0
    sector_rows: list[dict[str, object]] = []

    for l2_label in l2_order:
        start = cursor
        end = cursor + width_l2
        sector_df = df.loc[df[l2_col].eq(l2_label)].copy()
        tissue_counts = (
            sector_df.groupby(tissue_col, observed=False)
            .size()
            .reindex(tissue_order, fill_value=0)
        )
        total = tissue_counts.sum()
        tissue_cursor = start
        for tissue in tissue_order:
            frac = 0.0 if total == 0 else float(tissue_counts[tissue]) / float(total)
            seg_width = width_l2 * frac
            if seg_width <= 0:
                continue
            polar_ax.bar(
                x=tissue_cursor + seg_width / 2,
                height=inner_r1 - inner_r0,
                width=seg_width,
                bottom=inner_r0,
                color=TISSUE_COLOR_MAP[tissue],
                edgecolor="white",
                linewidth=0.55,
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

        polar_ax.bar(
            x=(start + end) / 2,
            height=inner_r1 - inner_r0,
            width=width_l2,
            bottom=inner_r0,
            color="none",
            edgecolor="#444444",
            linewidth=0.9,
            align="center",
            zorder=7,
        )
        add_sector_label(
            polar_ax,
            theta=(start + end) / 2,
            radius=inner_r0 - 0.028,
            text=L2_SHORT_LABELS.get(l2_label, l2_label),
            color=L2_COLOR_MAP.get(l2_label, "#333333"),
            fontsize=9.2,
            weight="bold",
        )

        l3_labels = l3_order_by_l2.get(l2_label, [])
        width_l3 = width_l2 / max(len(l3_labels), 1)
        l3_cursor = start
        for l3_label in l3_labels:
            l3_start = l3_cursor
            l3_end = l3_cursor + width_l3
            l3_df = sector_df.loc[sector_df[l3_col].eq(l3_label)].copy()
            l3_tissue_counts = (
                l3_df.groupby(tissue_col, observed=False)
                .size()
                .reindex(tissue_order, fill_value=0)
            )
            l3_total = l3_tissue_counts.sum()
            sub_cursor = l3_start
            for tissue in tissue_order:
                frac = 0.0 if l3_total == 0 else float(l3_tissue_counts[tissue]) / float(l3_total)
                seg_width = width_l3 * frac
                if seg_width <= 0:
                    continue
                polar_ax.bar(
                    x=sub_cursor + seg_width / 2,
                    height=outer_r1 - outer_r0,
                    width=seg_width,
                    bottom=outer_r0,
                    color=TISSUE_COLOR_MAP[tissue],
                    edgecolor="white",
                    linewidth=0.45,
                    align="center",
                    zorder=8,
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

            polar_ax.bar(
                x=(l3_start + l3_end) / 2,
                height=outer_r1 - outer_r0,
                width=width_l3,
                bottom=outer_r0,
                color="none",
                edgecolor="#444444",
                linewidth=0.78,
                align="center",
                zorder=9,
            )
            add_sector_label(
                polar_ax,
                theta=(l3_start + l3_end) / 2,
                radius=outer_r1 + 0.032,
                text=L3_SHORT_LABELS.get(l3_label, l3_label),
                color=L3_COLOR_MAP.get(l3_label, "#333333"),
                fontsize=7.3,
            )
            l3_cursor = l3_end

        cursor = end

    polar_ax.set_theta_offset(np.pi / 2.0)
    polar_ax.set_theta_direction(-1)
    polar_ax.set_ylim(0, 0.57)
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
            fontsize=11.5,
            fontweight="bold",
            color="#1F1F1F",
            ha="center",
            va="center",
            zorder=6,
        )
        text.set_path_effects([pe.Stroke(linewidth=3.0, foreground="white"), pe.Normal()])


def style_main_axes(ax: plt.Axes, df: pd.DataFrame, title: str, umap_key: str) -> None:
    x = df["umap_1"].to_numpy()
    y = df["umap_2"].to_numpy()
    xr = float(np.max(x) - np.min(x))
    yr = float(np.max(y) - np.min(y))
    pad = 0.28 * max(xr, yr)
    ax.set_xlim(float(np.min(x) - pad), float(np.max(x) + pad))
    ax.set_ylim(float(np.min(y) - pad), float(np.max(y) + pad))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=18, fontweight="bold", pad=18)
    ax.text(
        0.5,
        1.01,
        f"center = {umap_key} colored by L2 | inner ring = tissue fraction within L2 | outer ring = tissue fraction within L3",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.8,
        color="#4A4A4A",
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

    fig = plt.figure(figsize=(12, 11), facecolor="white")
    ax = fig.add_axes([0.04, 0.08, 0.78, 0.84])

    for label in l2_order:
        sub = df.loc[df[l2_col].eq(label)]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=7,
            c=L2_COLOR_MAP.get(label, "#999999"),
            linewidths=0.0,
            alpha=0.96,
            rasterized=True,
            zorder=3,
        )

    add_global_contour(ax, df["umap_1"].to_numpy(), df["umap_2"].to_numpy())
    annotate_l2_centroids(ax, df, l2_col)
    style_main_axes(ax, df, title, umap_key)

    ring_ax = fig.add_axes(ax.get_position(), projection="polar", facecolor="none")
    ring_sector_df = draw_nested_rings(ring_ax, df, l2_col, l3_col, tissue_col, l2_order, l3_order_by_l2)
    if not ring_sector_df.empty:
        ring_sector_df.to_csv(output_dir / f"{prefix}_ring_segments.tsv", sep="\t", index=False)

    tissue_handles = [Patch(facecolor=color, edgecolor="none", label=label) for label, color in TISSUE_COLOR_MAP.items() if label in tissue_order]
    fig.legend(
        handles=tissue_handles,
        loc="center right",
        bbox_to_anchor=(0.98, 0.53),
        frameon=False,
        title="Tissue",
        title_fontsize=12,
        fontsize=10,
    )

    fig.text(
        0.83,
        0.13,
        "Outer labels:\nAtypMem, GC-DZ, GC-LZ, GC-Trans,\nMemory, Naive, Pl-IgA, Pl-IgG",
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#4A4A4A",
    )

    png_path = output_dir / f"{prefix}.png"
    pdf_path = output_dir / f"{prefix}.pdf"
    fig.savefig(png_path, dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, dpi=320, bbox_inches="tight", facecolor="white")
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
            "UMAP center uses L2 colors and only global occupancy contours.",
            "Ring sectors use fixed-width label slots; within-sector angle encodes tissue proportion.",
        ],
    }
    (output_dir / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
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
