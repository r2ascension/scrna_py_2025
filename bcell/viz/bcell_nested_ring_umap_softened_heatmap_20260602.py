#!/usr/bin/env python3
"""B-cell nested-ring UMAP with enlarged rotated rings and sampled-cell heatmap.

This version extends the softened preview by:
  - pushing the tissue rings farther outward again
  - narrowing the L2/L3 ring bands slightly
  - rotating the rendered ring system 90 degrees clockwise
  - adding an outer sampled-cell heatmap ring based on per-cell expression

The center UMAP keeps the existing design principle:
  - L2-colored outer hull
  - no explicit internal borders between neighboring cell groups
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Dict, Sequence

import anndata as ad
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb, to_hex
from matplotlib.patches import Patch
from scipy import sparse
from scipy.ndimage import (
    binary_closing,
    binary_fill_holes,
    binary_opening,
    gaussian_filter,
    label as ndi_label,
)

BASE_SCRIPT_PATH = Path(__file__).with_name("bcell_nested_ring_umap_refined_20260531.py")
DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/"
    "bcell_tissue_comparison_final_fullgene.h5ad"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/h2048/data/py/20260602/bcell_nested_ring_umap_softened_heatmap_20260602"
)
DEFAULT_PREFIX = "bcell_nested_ring_umap_softened_heatmap_20260602"
DEFAULT_ROTATION_DEG = 90.0
DEFAULT_HEATMAP_GENES: tuple[str, ...] = (
    "TCL1A",
    "CD27",
    "FCRL5",
    "RGS13",
    "MKI67",
    "AICDA",
    "IGHA1",
    "IGHG1",
)
HEATMAP_FALLBACKS: Dict[str, tuple[str, ...]] = {
    "CD27": ("BANK1",),
    "RGS13": ("CD83", "BCL6"),
    "AICDA": ("BCL6",),
    "IGHA1": ("IGHA2",),
    "IGHG1": ("IGHG3", "IGHG4"),
}
DEFAULT_HEATMAP_LAYER = "log1p"
DEFAULT_HEATMAP_TOTAL_BINS = 420
DEFAULT_MIN_BINS_PER_SEGMENT = 3

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
HEATMAP_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "bcell_soft_expression",
    ["#FCFBFD", "#F1E5F5", "#DABEE7", "#AF80CA", "#7A49A6"],
)

INNER_R0 = 1.482
INNER_R1 = 1.625
OUTER_CORE_R0 = 1.781
OUTER_CORE_R1 = 1.911
OUTER_FEATHER = 0.08
HEATMAP_GAP = 0.08
HEATMAP_R0 = OUTER_CORE_R1 + OUTER_FEATHER + HEATMAP_GAP
HEATMAP_R1 = 2.652
L2_GAP_FRACTION = 0.030
L3_GAP_FRACTION = 0.045
POLAR_RMAX = HEATMAP_R1 + 0.23

base_module = None


def _darken(color: str, factor: float = 0.65) -> str:
    rgb = np.array(to_rgb(color))
    return to_hex(np.clip(rgb * factor, 0, 1))


def _sample_idx(n_total: int, max_points: int, rng: np.random.Generator) -> np.ndarray:
    if max_points <= 0 or n_total <= max_points:
        return np.arange(n_total)
    return np.sort(rng.choice(n_total, size=max_points, replace=False))


def _choose_density_level(hist: np.ndarray, n_points: int) -> float | None:
    nonzero = hist[hist > 0]
    if nonzero.size == 0:
        return None
    hist_max = float(nonzero.max())
    quantile = 0.35 if n_points >= 1000 else 0.25
    level = max(float(np.quantile(nonzero, quantile)), hist_max * 0.12)
    if level >= hist_max:
        level = hist_max * 0.50
    if level <= 0:
        return None
    return level


def _build_component_masks(hist: np.ndarray, level: float) -> list[np.ndarray]:
    dense_mask = hist >= level
    if not np.any(dense_mask):
        return []
    structure = np.ones((3, 3), dtype=bool)
    dense_mask = binary_opening(dense_mask, structure=structure, iterations=1)
    dense_mask = binary_closing(dense_mask, structure=structure, iterations=1)
    dense_mask = binary_fill_holes(dense_mask)
    if not np.any(dense_mask):
        return []
    labeled, n_components = ndi_label(dense_mask)
    masks: list[np.ndarray] = []
    for cid in range(1, n_components + 1):
        cm = labeled == cid
        if int(cm.sum()) >= 6:
            masks.append(cm)
    return masks


def _nearest_winner(
    midpoint: np.ndarray,
    winner: np.ndarray,
    occupied: np.ndarray,
    xcenters: np.ndarray,
    ycenters: np.ndarray,
    categories: list[str],
) -> str | None:
    ix = int(np.clip(np.searchsorted(xcenters, midpoint[0]), 0, len(xcenters) - 1))
    iy = int(np.clip(np.searchsorted(ycenters, midpoint[1]), 0, len(ycenters) - 1))
    for radius in range(0, 5):
        x0, x1 = max(ix - radius, 0), min(ix + radius + 1, len(xcenters))
        y0, y1 = max(iy - radius, 0), min(iy + radius + 1, len(ycenters))
        local = occupied[x0:x1, y0:y1]
        if not np.any(local):
            continue
        ly, lx = np.where(local)
        gx, gy = xcenters[x0:x1][lx], ycenters[y0:y1][ly]
        dists = (gx - midpoint[0]) ** 2 + (gy - midpoint[1]) ** 2
        wi = int(winner[x0 + ly[int(np.argmin(dists))], y0 + lx[int(np.argmin(dists))]])
        if 0 <= wi < len(categories):
            return categories[wi]
    return None


def add_per_group_outer_boundary(
    ax: plt.Axes,
    df: pd.DataFrame,
    group_col: str,
    group_order: list[str],
    color_map: dict[str, str],
    max_points: int = 20000,
    grid_size: int = 260,
    seed: int = 42,
) -> None:
    """Boundary hull based on `prepare_outer_edge_segments` from the epithelial viz
    pipeline: per-group quantile-based density levels, morphological cleanup, connected-
    component detection, union occupied mask contour, and winner-based segment coloring."""
    rng = np.random.default_rng(seed)
    x = df["umap_1"].to_numpy(dtype=float)
    y = df["umap_2"].to_numpy(dtype=float)
    labels_arr = df[group_col].to_numpy()
    categories = [c for c in group_order if c in set(labels_arr)]
    if len(categories) < 2:
        return

    xr = float(x.max() - x.min())
    yr = float(y.max() - y.min())
    if xr <= 0 or yr <= 0:
        return

    x_range = (float(x.min()) - xr * 0.03, float(x.max()) + xr * 0.03)
    y_range = (float(y.min()) - yr * 0.03, float(y.max()) + yr * 0.03)
    gs = min(260, grid_size)

    density_stack: list[np.ndarray] = []
    active_masks: list[np.ndarray] = []
    for label in categories:
        pts = np.column_stack([x[labels_arr == label], y[labels_arr == label]])
        if pts.shape[0] > max_points:
            pts = pts[_sample_idx(pts.shape[0], max_points, rng)]
        hist, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], bins=gs, range=[list(x_range), list(y_range)])
        hist = gaussian_filter(hist, sigma=1.2)
        level = _choose_density_level(hist, int(pts.shape[0]))
        masks = _build_component_masks(hist, level) if level is not None else []
        active = np.zeros_like(hist, dtype=bool)
        for m in masks:
            active |= m
        if not masks and level is not None:
            active = hist >= level
        density_stack.append(hist)
        active_masks.append(active)

    occupied = np.zeros_like(density_stack[0], dtype=bool)
    for am in active_masks:
        occupied |= am
    occupied = binary_closing(occupied, structure=np.ones((3, 3), dtype=bool), iterations=1)
    occupied = binary_fill_holes(occupied)

    score_stack = np.stack(
        [np.where(am, h, -np.inf) for h, am in zip(density_stack, active_masks)], axis=0,
    )
    winner_grid = np.argmax(score_stack, axis=0)

    xedges = np.linspace(x_range[0], x_range[1], gs + 1)
    yedges = np.linspace(y_range[0], y_range[1], gs + 1)
    xcenters = (xedges[:-1] + xedges[1:]) / 2.0
    ycenters = (yedges[:-1] + yedges[1:]) / 2.0

    contour_field = occupied.astype(float)
    contour_field = gaussian_filter(contour_field, sigma=1.0)

    # extract contour via invisible axes
    tmp_fig, tmp_ax = plt.subplots(figsize=(2, 2))
    try:
        contour = tmp_ax.contour(xcenters, ycenters, contour_field.T, levels=[0.5], colors=["black"], linewidths=1.0)
        paths = [seg for seg in contour.allsegs[0] if len(seg) >= 2]
    finally:
        plt.close(tmp_fig)

    segments: list[np.ndarray] = []
    seg_colors: list[str] = []
    for path in paths:
        for p0, p1 in zip(path[:-1], path[1:]):
            mid = (p0 + p1) / 2.0
            winner_label = _nearest_winner(mid, winner_grid, occupied, xcenters, ycenters, categories)
            if winner_label is None:
                continue
            segments.append(np.vstack([p0, p1]))
            seg_colors.append(_darken(color_map.get(winner_label, "#444444"), 0.82))

    if not segments:
        return

    underlay = LineCollection(
        segments, colors=[(1, 1, 1, 0.96)] * len(segments),
        linewidths=3.6, capstyle="round", joinstyle="round", zorder=4.35,
    )
    ax.add_collection(underlay)
    overlay = LineCollection(
        segments, colors=seg_colors,
        linewidths=2.2, capstyle="round", joinstyle="round", zorder=4.55, alpha=0.94,
    )
    ax.add_collection(overlay)


def load_base_module():
    spec = importlib.util.spec_from_file_location(
        "bcell_nested_ring_umap_refined_20260531",
        BASE_SCRIPT_PATH,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load base script from {BASE_SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create enlarged rotated B-cell nested-ring UMAP with sampled-cell heatmap ring."
    )
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--l2-col", default="cell_type_L2")
    parser.add_argument("--l3-col", default="cell_type_L3")
    parser.add_argument("--tissue-col", default="tissue")
    parser.add_argument("--umap-key", default=None)
    parser.add_argument("--title", default="B-cell nested ring UMAP sampled heatmap preview")
    parser.add_argument(
        "--heatmap-genes",
        default=",".join(DEFAULT_HEATMAP_GENES),
        help="Comma-separated heatmap genes in inner-to-outer order.",
    )
    parser.add_argument(
        "--heatmap-layer",
        default=DEFAULT_HEATMAP_LAYER,
        help="Expression layer used for heatmap extraction. Use X to read adata.X.",
    )
    parser.add_argument(
        "--heatmap-total-bins",
        type=int,
        default=DEFAULT_HEATMAP_TOTAL_BINS,
        help="Approximate total number of sampled-cell angular bins across the full heatmap ring.",
    )
    parser.add_argument(
        "--heatmap-min-bins-per-segment",
        type=int,
        default=DEFAULT_MIN_BINS_PER_SEGMENT,
        help="Minimum sampled cells per non-empty L3×tissue segment.",
    )
    parser.add_argument(
        "--rotation-deg",
        type=float,
        default=DEFAULT_ROTATION_DEG,
        help="Clockwise rotation applied to the rendered ring system.",
    )
    return parser.parse_args()


def parse_gene_list(raw_value: str) -> list[str]:
    genes = [gene.strip() for gene in raw_value.split(",") if gene.strip()]
    if not genes:
        raise ValueError("At least one heatmap gene is required")
    return genes


def patch_soft_palettes() -> None:
    base_module.L2_COLOR_MAP = SOFT_L2_COLOR_MAP.copy()
    base_module.L3_COLOR_MAP = SOFT_L3_COLOR_MAP.copy()
    base_module.TISSUE_COLOR_MAP = SOFT_TISSUE_COLOR_MAP.copy()


def load_adata_and_plot_df(
    input_h5ad: Path,
    l2_col: str,
    l3_col: str,
    tissue_col: str,
    umap_key: str | None,
) -> tuple[ad.AnnData, pd.DataFrame, str]:
    adata = ad.read_h5ad(input_h5ad)
    chosen_umap_key = base_module.choose_umap_key(adata, umap_key)
    coords = np.asarray(adata.obsm[chosen_umap_key])
    if coords.shape[1] < 2:
        raise ValueError(f"UMAP key {chosen_umap_key} does not have at least 2 dimensions")

    required_cols = [l2_col, l3_col, tissue_col]
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    plot_df = pd.DataFrame(
        {"umap_1": coords[:, 0], "umap_2": coords[:, 1]},
        index=adata.obs_names.copy(),
    ).join(adata.obs[required_cols].copy())

    plot_df = plot_df.dropna(subset=[l2_col, l3_col, tissue_col, "umap_1", "umap_2"]).copy()
    for col in [l2_col, l3_col, tissue_col]:
        plot_df[col] = plot_df[col].astype(str).str.strip()
        plot_df = plot_df.loc[plot_df[col].ne("")].copy()

    adata = adata[plot_df.index].copy()
    return adata, plot_df, chosen_umap_key


def resolve_heatmap_genes(
    requested_genes: Sequence[str],
    var_names: Sequence[str],
) -> tuple[list[str], dict[str, str]]:
    available = set(map(str, var_names))
    resolved: list[str] = []
    mapping: dict[str, str] = {}
    used: set[str] = set()

    for gene in requested_genes:
        candidates = (gene, *HEATMAP_FALLBACKS.get(gene, ()))
        chosen = next(
            (
                candidate
                for candidate in candidates
                if candidate in available and candidate not in used
            ),
            None,
        )
        if chosen is None:
            continue
        resolved.append(chosen)
        mapping[gene] = chosen
        used.add(chosen)

    if not resolved:
        raise ValueError("None of the requested heatmap genes are available in the AnnData object")
    return resolved, mapping


def get_expression_matrix(adata: ad.AnnData, gene_names: Sequence[str], layer: str) -> pd.DataFrame:
    var_index = pd.Index(map(str, adata.var_names))
    gene_positions = var_index.get_indexer(gene_names)
    if (gene_positions < 0).any():
        missing = [gene for gene, idx in zip(gene_names, gene_positions) if idx < 0]
        raise KeyError(f"Missing resolved heatmap genes in AnnData: {missing}")

    if layer.upper() == "X":
        source = adata.X[:, gene_positions]
    else:
        if layer not in adata.layers:
            raise KeyError(f"Requested heatmap layer not found: {layer}")
        source = adata.layers[layer][:, gene_positions]

    if sparse.issparse(source):
        matrix = source.toarray()
    else:
        matrix = np.asarray(source)

    return pd.DataFrame(matrix, index=adata.obs_names.copy(), columns=list(gene_names))


def scale_expression_matrix(expr_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scaled = pd.DataFrame(index=expr_df.index)
    rows: list[dict[str, float | str]] = []
    for gene in expr_df.columns:
        values = expr_df[gene].to_numpy(dtype=float)
        lower = float(np.nanpercentile(values, 1.0))
        upper = float(np.nanpercentile(values, 99.0))
        if upper <= lower:
            upper = lower + 1e-6
        scaled[gene] = np.clip((values - lower) / (upper - lower), 0.0, 1.0)
        rows.append(
            {
                "gene": gene,
                "q01": lower,
                "q99": upper,
                "min": float(np.nanmin(values)),
                "max": float(np.nanmax(values)),
            }
        )
    return scaled, pd.DataFrame(rows)


def softened_style_main_axes(ax, df: pd.DataFrame, title: str, umap_key: str) -> None:
    x = df["umap_1"].to_numpy()
    y = df["umap_2"].to_numpy()
    xr = float(x.max() - x.min())
    yr = float(y.max() - y.min())
    pad = 0.93 * max(xr, yr)
    ax.set_xlim(float(x.min() - pad), float(x.max() + pad))
    ax.set_ylim(float(y.min() - pad), float(y.max() + pad))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=15.8, fontweight="bold", pad=16)
    ax.text(
        0.5,
        1.01,
        f"center = {umap_key} with L2-colored outer hull | enlarged rotated soft rings + sampled-cell outer heatmap",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.9,
        color="#666666",
    )


def apply_rotation(theta: float | np.ndarray, rotation_deg: float) -> float | np.ndarray:
    return np.asarray(theta) + np.deg2rad(rotation_deg)


def draw_feathered_outer_segment(
    polar_ax: plt.Axes,
    theta: float,
    width: float,
    bottom: float,
    core_height: float,
    feather_height: float,
    color: str,
    zorder: float,
) -> None:
    polar_ax.bar(
        x=theta,
        height=core_height,
        width=width,
        bottom=bottom,
        color=base_module.adjust_lightness(color, amount=1.01, alpha=0.88),
        edgecolor="none",
        linewidth=0.0,
        align="center",
        zorder=zorder,
    )
    if feather_height <= 0:
        return
    n_steps = 4
    step = feather_height / n_steps
    alphas = [0.19, 0.12, 0.07, 0.035]
    for idx, alpha in enumerate(alphas):
        polar_ax.bar(
            x=theta,
            height=step,
            width=width * (1.002 + 0.002 * idx),
            bottom=bottom + core_height + idx * step,
            color=base_module.adjust_lightness(color, amount=1.05 + 0.03 * idx, alpha=alpha),
            edgecolor="none",
            linewidth=0.0,
            align="center",
            zorder=zorder + 0.12 + idx * 0.01,
        )


def draw_nested_rings(
    polar_ax: plt.Axes,
    df: pd.DataFrame,
    l2_col: str,
    l3_col: str,
    tissue_col: str,
    l2_order: Sequence[str],
    l3_order_by_l2: dict[str, list[str]],
    rotation_deg: float,
) -> pd.DataFrame:
    tissue_order = [tissue for tissue in SOFT_TISSUE_COLOR_MAP if tissue in set(df[tissue_col])]
    width_l2 = 2 * np.pi / max(len(l2_order), 1)
    cursor = 0.0
    sector_rows: list[dict[str, object]] = []

    for l2_label in l2_order:
        start = cursor
        end = cursor + width_l2
        l2_draw_start, l2_draw_end, l2_draw_width, _ = base_module.sector_span_with_gap(
            start,
            end,
            gap_fraction=L2_GAP_FRACTION,
            max_gap_deg=3.0,
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
            seg_start = tissue_cursor
            seg_end = tissue_cursor + seg_width
            polar_ax.bar(
                x=float(apply_rotation(seg_start + seg_width / 2, rotation_deg)),
                height=INNER_R1 - INNER_R0,
                width=seg_width,
                bottom=INNER_R0,
                color=base_module.adjust_lightness(
                    SOFT_TISSUE_COLOR_MAP[tissue],
                    amount=1.01,
                    alpha=0.88,
                ),
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
                    "start_rad_raw": seg_start,
                    "end_rad_raw": seg_end,
                    "start_rad_render": float(apply_rotation(seg_start, rotation_deg)),
                    "end_rad_render": float(apply_rotation(seg_end, rotation_deg)),
                    "width_rad": seg_width,
                    "fraction_within_group": frac,
                    "n_cells": int(tissue_counts[tissue]),
                }
            )
            tissue_cursor = seg_end

        if l2_draw_width > 0:
            polar_ax.bar(
                x=float(apply_rotation((l2_draw_start + l2_draw_end) / 2, rotation_deg)),
                height=INNER_R1 - INNER_R0,
                width=l2_draw_width,
                bottom=INNER_R0,
                color="none",
                edgecolor=base_module.adjust_lightness(
                    SOFT_L2_COLOR_MAP.get(l2_label, "#777777"),
                    amount=0.82,
                    alpha=0.92,
                ),
                linewidth=1.34,
                align="center",
                zorder=9,
            )
            polar_ax.bar(
                x=float(apply_rotation((l2_draw_start + l2_draw_end) / 2, rotation_deg)),
                height=0.010,
                width=l2_draw_width,
                bottom=INNER_R1 - 0.010,
                color=base_module.adjust_lightness(
                    SOFT_L2_COLOR_MAP.get(l2_label, "#777777"),
                    amount=1.03,
                    alpha=0.90,
                ),
                edgecolor="none",
                align="center",
                zorder=10,
            )

        base_module.add_sector_label(
            polar_ax,
            theta=float(apply_rotation((start + end) / 2, rotation_deg)),
            radius=INNER_R0 - 0.040,
            text=base_module.L2_SHORT_LABELS.get(l2_label, l2_label),
            color=SOFT_L2_COLOR_MAP.get(l2_label, "#555555"),
            fontsize=9.6,
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
                gap_fraction=L3_GAP_FRACTION,
                max_gap_deg=2.3,
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
                seg_start = sub_cursor
                seg_end = sub_cursor + seg_width
                draw_feathered_outer_segment(
                    polar_ax=polar_ax,
                    theta=float(apply_rotation(seg_start + seg_width / 2, rotation_deg)),
                    width=seg_width,
                    bottom=OUTER_CORE_R0,
                    core_height=OUTER_CORE_R1 - OUTER_CORE_R0,
                    feather_height=OUTER_FEATHER,
                    color=SOFT_TISSUE_COLOR_MAP[tissue],
                    zorder=7,
                )
                sector_rows.append(
                    {
                        "level": "L3",
                        "label": l3_label,
                        "parent": l2_label,
                        "tissue": tissue,
                        "start_rad_raw": seg_start,
                        "end_rad_raw": seg_end,
                        "start_rad_render": float(apply_rotation(seg_start, rotation_deg)),
                        "end_rad_render": float(apply_rotation(seg_end, rotation_deg)),
                        "width_rad": seg_width,
                        "fraction_within_group": frac,
                        "n_cells": int(l3_tissue_counts[tissue]),
                    }
                )
                sub_cursor = seg_end

            if l3_draw_width > 0:
                polar_ax.bar(
                    x=float(apply_rotation((l3_draw_start + l3_draw_end) / 2, rotation_deg)),
                    height=(OUTER_CORE_R1 - OUTER_CORE_R0) + OUTER_FEATHER * 0.55,
                    width=l3_draw_width,
                    bottom=OUTER_CORE_R0,
                    color="none",
                    edgecolor=base_module.adjust_lightness(
                        SOFT_L3_COLOR_MAP.get(l3_label, "#777777"),
                        amount=0.84,
                        alpha=0.84,
                    ),
                    linewidth=0.98,
                    align="center",
                    zorder=12,
                )

            base_module.add_sector_label(
                polar_ax,
                theta=float(apply_rotation((l3_start + l3_end) / 2, rotation_deg)),
                radius=HEATMAP_R1 + 0.070,
                text=base_module.L3_SHORT_LABELS.get(l3_label, l3_label),
                color=SOFT_L3_COLOR_MAP.get(l3_label, "#555555"),
                fontsize=7.5,
                weight="normal",
                alpha=0.86,
            )
            l3_cursor = l3_end

        cursor = end

    polar_ax.set_theta_offset(np.pi / 2.0)
    polar_ax.set_theta_direction(-1)
    polar_ax.set_ylim(0, POLAR_RMAX)
    polar_ax.set_axis_off()

    sector_df = pd.DataFrame(sector_rows)
    if not sector_df.empty:
        sector_df["start_deg_clockwise_from_top"] = np.degrees(sector_df["start_rad_raw"])
        sector_df["end_deg_clockwise_from_top"] = np.degrees(sector_df["end_rad_raw"])
        sector_df["start_deg_rendered"] = np.mod(
            np.degrees(sector_df["start_rad_render"]),
            360.0,
        )
        sector_df["end_deg_rendered"] = np.mod(
            np.degrees(sector_df["end_rad_render"]),
            360.0,
        )
    return sector_df


def annotate_l2_centroids(ax: plt.Axes, df: pd.DataFrame, l2_col: str) -> None:
    centroids = base_module.compute_group_centroids(df, l2_col)
    for row in centroids.itertuples(index=False):
        label = getattr(row, l2_col)
        text = ax.text(
            row.cx,
            row.cy,
            label,
            fontsize=11.1,
            fontweight="bold",
            color="#2C2C2C",
            ha="center",
            va="center",
            zorder=6,
        )
        text.set_path_effects(
            [
                base_module.pe.Stroke(linewidth=3.1, foreground="white"),
                base_module.pe.Normal(),
            ]
        )


def clockwise_angle_from_top(x: np.ndarray, y: np.ndarray, cx: float, cy: float) -> np.ndarray:
    angle = np.degrees(np.arctan2(y - cy, x - cx))
    angle = np.mod(angle, 360.0)
    return np.mod(90.0 - angle, 360.0)


def sample_cells_by_ring_segment(
    plot_df: pd.DataFrame,
    expr_log1p_df: pd.DataFrame,
    expr_scaled_df: pd.DataFrame,
    l3_centroids: pd.DataFrame,
    l3_segments: pd.DataFrame,
    l3_col: str,
    tissue_col: str,
    heatmap_genes: Sequence[str],
    target_total_bins: int,
    min_bins_per_segment: int,
    rotation_deg: float,
) -> pd.DataFrame:
    if l3_segments.empty:
        return pd.DataFrame()

    centroid_lookup = l3_centroids.set_index(l3_col)[["cx", "cy"]].to_dict(orient="index")
    total_width = float(l3_segments["width_rad"].sum())
    if total_width <= 0:
        total_width = 1.0
    sampled_rows: list[dict[str, object]] = []

    ordered_segments = l3_segments.sort_values(["start_rad_raw", "label", "tissue"]).reset_index(drop=True)
    for segment_order, segment in enumerate(ordered_segments.itertuples(index=False), start=1):
        segment_cells = plot_df.loc[
            plot_df[l3_col].eq(segment.label) & plot_df[tissue_col].eq(segment.tissue)
        ].copy()
        if segment_cells.empty:
            continue

        width_fraction = float(segment.width_rad) / total_width
        target_bins = int(round(target_total_bins * width_fraction))
        target_bins = max(target_bins, min(min_bins_per_segment, int(segment.n_cells)))
        n_bins = min(int(segment.n_cells), max(target_bins, 1))

        centroid = centroid_lookup[str(segment.label)]
        segment_cells["local_angle_clockwise_from_top"] = clockwise_angle_from_top(
            segment_cells["umap_1"].to_numpy(),
            segment_cells["umap_2"].to_numpy(),
            float(centroid["cx"]),
            float(centroid["cy"]),
        )
        segment_cells = segment_cells.sort_values(
            ["local_angle_clockwise_from_top", "umap_1", "umap_2"]
        ).copy()

        sample_positions = np.floor(
            (np.arange(n_bins) + 0.5) * float(segment_cells.shape[0]) / float(n_bins)
        ).astype(int)
        sample_positions = np.clip(sample_positions, 0, segment_cells.shape[0] - 1)
        sampled_cells = segment_cells.iloc[sample_positions].copy()

        theta_edges_raw = np.linspace(float(segment.start_rad_raw), float(segment.end_rad_raw), n_bins + 1)
        theta_edges_render = apply_rotation(theta_edges_raw, rotation_deg)

        for bin_index, (cell_id, row) in enumerate(sampled_cells.iterrows(), start=1):
            sampled_row: dict[str, object] = {
                "sample_id": f"seg{segment_order:02d}_bin{bin_index:03d}",
                "cell_id": cell_id,
                "parent_l2": segment.parent,
                "l3_label": segment.label,
                "tissue": segment.tissue,
                "segment_order": segment_order,
                "segment_n_cells": int(segment.n_cells),
                "segment_width_rad": float(segment.width_rad),
                "segment_width_deg": float(np.degrees(segment.width_rad)),
                "n_bins_segment": n_bins,
                "bin_index_within_segment": bin_index,
                "sample_rank_from_sorted_segment": int(sample_positions[bin_index - 1]) + 1,
                "theta_start_raw": float(theta_edges_raw[bin_index - 1]),
                "theta_end_raw": float(theta_edges_raw[bin_index]),
                "theta_mid_raw": float(0.5 * (theta_edges_raw[bin_index - 1] + theta_edges_raw[bin_index])),
                "theta_start_render": float(theta_edges_render[bin_index - 1]),
                "theta_end_render": float(theta_edges_render[bin_index]),
                "theta_mid_render": float(0.5 * (theta_edges_render[bin_index - 1] + theta_edges_render[bin_index])),
                "theta_start_deg_rendered": float(np.mod(np.degrees(theta_edges_render[bin_index - 1]), 360.0)),
                "theta_end_deg_rendered": float(np.mod(np.degrees(theta_edges_render[bin_index]), 360.0)),
                "theta_mid_deg_rendered": float(
                    np.mod(
                        np.degrees(0.5 * (theta_edges_render[bin_index - 1] + theta_edges_render[bin_index])),
                        360.0,
                    )
                ),
                "local_angle_clockwise_from_top": float(row["local_angle_clockwise_from_top"]),
            }
            for gene in heatmap_genes:
                sampled_row[f"expr_{gene}"] = float(expr_log1p_df.at[cell_id, gene])
                sampled_row[f"scaled_{gene}"] = float(expr_scaled_df.at[cell_id, gene])
            sampled_rows.append(sampled_row)

    sampled_df = pd.DataFrame(sampled_rows)
    if sampled_df.empty:
        return sampled_df
    return sampled_df.sort_values(["segment_order", "bin_index_within_segment", "sample_id"]).reset_index(drop=True)


def build_heatmap_column_payload(
    sampled_cells_df: pd.DataFrame,
    heatmap_genes: Sequence[str],
) -> tuple[list[tuple[float, float]], np.ndarray]:
    if sampled_cells_df.empty:
        return [], np.empty((len(heatmap_genes), 0), dtype=float)

    column_edges: list[tuple[float, float]] = []
    column_values: list[np.ndarray] = []
    prev_end: float | None = None

    scaled_cols = [f"scaled_{gene}" for gene in heatmap_genes]
    payload_df = sampled_cells_df.loc[:, ["theta_start_render", "theta_end_render", *scaled_cols]]

    for row in payload_df.itertuples(index=False, name=None):
        start = float(row[0])
        end = float(row[1])
        if prev_end is not None and start > prev_end + 1e-9:
            column_edges.append((prev_end, start))
            column_values.append(np.full(len(heatmap_genes), np.nan, dtype=float))
        column_edges.append((start, end))
        column_values.append(np.asarray(row[2:], dtype=float))
        prev_end = end

    if not column_edges:
        return [], np.empty((len(heatmap_genes), 0), dtype=float)
    matrix = np.column_stack(column_values)
    return column_edges, matrix


def draw_outer_heatmap_ring(
    polar_ax: plt.Axes,
    sampled_cells_df: pd.DataFrame,
    heatmap_genes: Sequence[str],
):
    column_edges, heatmap_matrix = build_heatmap_column_payload(sampled_cells_df, heatmap_genes)
    if not column_edges:
        return None

    theta_edges = np.asarray([column_edges[0][0], *[end for _, end in column_edges]], dtype=float)
    radial_edges = np.linspace(HEATMAP_R0, HEATMAP_R1, len(heatmap_genes) + 1)
    cmap = HEATMAP_CMAP.copy()
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    mesh = polar_ax.pcolormesh(
        theta_edges,
        radial_edges,
        heatmap_matrix,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        shading="flat",
        antialiased=False,
        zorder=13,
    )
    mesh.set_rasterized(True)
    return mesh


def format_gene_order_text(genes: Sequence[str]) -> str:
    chunks = [", ".join(genes[idx : idx + 4]) for idx in range(0, len(genes), 4)]
    return "\n".join(chunks)


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
    expr_log1p_df: pd.DataFrame,
    expr_scaled_df: pd.DataFrame,
    heatmap_genes_requested: Sequence[str],
    heatmap_genes_resolved: Sequence[str],
    gene_resolution: dict[str, str],
    scaling_df: pd.DataFrame,
    heatmap_layer: str,
    heatmap_total_bins: int,
    heatmap_min_bins_per_segment: int,
    rotation_deg: float,
) -> dict[str, object]:
    l2_order, l3_order_by_l2, l2_centroids, l3_centroids = base_module.build_orderings(df, l2_col, l3_col)
    tissue_order = [t for t in SOFT_TISSUE_COLOR_MAP if t in set(df[tissue_col])]

    l2_table = base_module.composition_table(df, l2_col, tissue_col, l2_order, tissue_order)
    l3_order_flat = [label for parent in l2_order for label in l3_order_by_l2.get(parent, [])]
    l3_table = base_module.composition_table(df, l3_col, tissue_col, l3_order_flat, tissue_order)
    sector_manifest = base_module.build_sector_manifest(
        df,
        l2_col,
        l3_col,
        tissue_col,
        l2_order,
        l3_order_by_l2,
        l2_centroids,
        l3_centroids,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    l2_table.to_csv(output_dir / f"{prefix}_l2_tissue.tsv", sep="\t", index=False)
    l3_table.to_csv(output_dir / f"{prefix}_l3_tissue.tsv", sep="\t", index=False)
    sector_manifest.to_csv(output_dir / f"{prefix}_group_sector_manifest.tsv", sep="\t", index=False)
    scaling_df.to_csv(output_dir / f"{prefix}_heatmap_scaling.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(18, 16), facecolor="white")
    ax = fig.add_axes([0.03, 0.05, 0.77, 0.88])

    for label in l2_order:
        sub = df.loc[df[l2_col].eq(label)]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=8.0,
            c=SOFT_L2_COLOR_MAP.get(label, "#999999"),
            linewidths=0.0,
            alpha=0.93,
            rasterized=True,
            zorder=3,
        )

    add_per_group_outer_boundary(
        ax,
        df,
        group_col=l2_col,
        group_order=l2_order,
        color_map=SOFT_L2_COLOR_MAP,
        max_points=20000,
    )
    annotate_l2_centroids(ax, df, l2_col)
    softened_style_main_axes(ax, df, title, umap_key)

    ring_ax = fig.add_axes(ax.get_position(), projection="polar", facecolor="none", frameon=False)
    ring_ax.patch.set_visible(False)
    ring_ax.set_axis_off()
    ring_ax.set_zorder(0)
    ax.set_facecolor("none")
    ax.set_zorder(1)
    ring_sector_df = draw_nested_rings(
        polar_ax=ring_ax,
        df=df,
        l2_col=l2_col,
        l3_col=l3_col,
        tissue_col=tissue_col,
        l2_order=l2_order,
        l3_order_by_l2=l3_order_by_l2,
        rotation_deg=rotation_deg,
    )
    ring_sector_df.to_csv(output_dir / f"{prefix}_ring_segments.tsv", sep="\t", index=False)

    sampled_cells_df = sample_cells_by_ring_segment(
        plot_df=df,
        expr_log1p_df=expr_log1p_df,
        expr_scaled_df=expr_scaled_df,
        l3_centroids=l3_centroids,
        l3_segments=ring_sector_df.loc[ring_sector_df["level"].eq("L3")].copy(),
        l3_col=l3_col,
        tissue_col=tissue_col,
        heatmap_genes=heatmap_genes_resolved,
        target_total_bins=heatmap_total_bins,
        min_bins_per_segment=heatmap_min_bins_per_segment,
        rotation_deg=rotation_deg,
    )
    if sampled_cells_df.empty:
        raise ValueError("No sampled cells were generated for the outer heatmap ring")

    mesh = draw_outer_heatmap_ring(ring_ax, sampled_cells_df, heatmap_genes_resolved)
    if mesh is None:
        raise ValueError("Failed to draw outer heatmap ring")

    heatmap_matrix_log1p_df = pd.DataFrame(
        sampled_cells_df[[f"expr_{gene}" for gene in heatmap_genes_resolved]].to_numpy(dtype=float).T,
        index=list(heatmap_genes_resolved),
        columns=sampled_cells_df["sample_id"].tolist(),
    )
    heatmap_matrix_scaled_df = pd.DataFrame(
        sampled_cells_df[[f"scaled_{gene}" for gene in heatmap_genes_resolved]].to_numpy(dtype=float).T,
        index=list(heatmap_genes_resolved),
        columns=sampled_cells_df["sample_id"].tolist(),
    )
    sampled_cells_df.to_csv(output_dir / f"{prefix}_sampled_cells.tsv", sep="\t", index=False)
    heatmap_matrix_log1p_df.to_csv(output_dir / f"{prefix}_heatmap_matrix.tsv", sep="\t")
    heatmap_matrix_scaled_df.to_csv(output_dir / f"{prefix}_heatmap_matrix_scaled.tsv", sep="\t")

    tissue_handles = [
        Patch(facecolor=color, edgecolor="none", label=label)
        for label, color in SOFT_TISSUE_COLOR_MAP.items()
        if label in tissue_order
    ]
    fig.legend(
        handles=tissue_handles,
        loc="center right",
        bbox_to_anchor=(0.985, 0.54),
        frameon=False,
        title="Tissue",
        title_fontsize=12,
        fontsize=10.4,
    )

    cax = fig.add_axes([0.855, 0.28, 0.018, 0.18])
    colorbar = fig.colorbar(
        ScalarMappable(norm=mcolors.Normalize(vmin=0.0, vmax=1.0), cmap=HEATMAP_CMAP),
        cax=cax,
    )
    colorbar.set_label("Scaled\nlog1p", fontsize=9.8)
    colorbar.ax.tick_params(labelsize=8.4)
    colorbar.outline.set_visible(False)

    fig.text(
        0.826,
        0.20,
        f"Heatmap genes (inner→outer):\n{format_gene_order_text(heatmap_genes_resolved)}",
        ha="left",
        va="top",
        fontsize=8.7,
        color="#555555",
    )
    fig.text(
        0.826,
        0.11,
        "Outer-ring abbreviations:\nAtypMem, GC-DZ, GC-LZ, GC-Trans,\nMemory, Naive, Pl-IgA, Pl-IgG",
        ha="left",
        va="top",
        fontsize=8.3,
        color="#555555",
    )

    png_path = output_dir / f"{prefix}.png"
    pdf_path = output_dir / f"{prefix}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.5, facecolor="white")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.5, facecolor="white")
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
        "rotation_deg": rotation_deg,
        "heatmap_layer": heatmap_layer,
        "heatmap_genes_requested": list(heatmap_genes_requested),
        "heatmap_genes_resolved": list(heatmap_genes_resolved),
        "heatmap_gene_resolution": gene_resolution,
        "heatmap_total_bins_target": heatmap_total_bins,
        "heatmap_min_bins_per_segment": heatmap_min_bins_per_segment,
        "n_sampled_cells": int(sampled_cells_df.shape[0]),
        "heatmap_radius": {"r0": HEATMAP_R0, "r1": HEATMAP_R1},
        "ring_radius": {
            "inner_r0": INNER_R0,
            "inner_r1": INNER_R1,
            "outer_core_r0": OUTER_CORE_R0,
            "outer_core_r1": OUTER_CORE_R1,
            "outer_feather": OUTER_FEATHER,
        },
        "companion_files": {
            "l2_tissue": f"{prefix}_l2_tissue.tsv",
            "l3_tissue": f"{prefix}_l3_tissue.tsv",
            "group_sector_manifest": f"{prefix}_group_sector_manifest.tsv",
            "ring_segments": f"{prefix}_ring_segments.tsv",
            "sampled_cells": f"{prefix}_sampled_cells.tsv",
            "heatmap_matrix_log1p": f"{prefix}_heatmap_matrix.tsv",
            "heatmap_matrix_scaled": f"{prefix}_heatmap_matrix_scaled.tsv",
            "heatmap_scaling": f"{prefix}_heatmap_scaling.tsv",
        },
        "notes": [
            "Ring bands are pushed farther outward again and rendered clockwise by 90 degrees relative to the softened preview.",
            "The outermost heatmap ring is built from sampled single cells aligned to the displayed L3×tissue angular segments.",
            "The center UMAP still only keeps an L2-colored outer hull and avoids explicit internal borders between adjacent groups.",
        ],
    }
    summary_path = output_dir / f"{prefix}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def main() -> None:
    global base_module
    base_module = load_base_module()
    patch_soft_palettes()

    args = parse_args()
    requested_genes = parse_gene_list(args.heatmap_genes)
    adata, plot_df, umap_key = load_adata_and_plot_df(
        input_h5ad=args.input_h5ad,
        l2_col=args.l2_col,
        l3_col=args.l3_col,
        tissue_col=args.tissue_col,
        umap_key=args.umap_key,
    )

    resolved_heatmap_genes, gene_resolution = resolve_heatmap_genes(requested_genes, adata.var_names)
    expr_log1p_df = get_expression_matrix(adata, resolved_heatmap_genes, args.heatmap_layer)
    expr_scaled_df, scaling_df = scale_expression_matrix(expr_log1p_df)

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
        expr_log1p_df=expr_log1p_df,
        expr_scaled_df=expr_scaled_df,
        heatmap_genes_requested=requested_genes,
        heatmap_genes_resolved=resolved_heatmap_genes,
        gene_resolution=gene_resolution,
        scaling_df=scaling_df,
        heatmap_layer=args.heatmap_layer,
        heatmap_total_bins=args.heatmap_total_bins,
        heatmap_min_bins_per_segment=args.heatmap_min_bins_per_segment,
        rotation_deg=args.rotation_deg,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
