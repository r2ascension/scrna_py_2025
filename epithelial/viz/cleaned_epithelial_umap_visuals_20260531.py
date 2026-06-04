#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import to_hex, to_rgb
from matplotlib.lines import Line2D

try:
    from scipy.ndimage import (
        binary_closing,
        binary_fill_holes,
        binary_opening,
        gaussian_filter,
        label as ndi_label,
    )
except Exception:  # pragma: no cover - optional dependency fallback
    gaussian_filter = None
    binary_opening = None
    binary_closing = None
    binary_fill_holes = None
    ndi_label = None

DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/py/20260531/main_lineage_h5ad_cleanup_20260531/epithelial/"
    "epithelial_tissue_comparison_final_fullgene_cleaned_l3refined_20260531.h5ad"
)
DEFAULT_REFERENCE_H5AD = Path(
    "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/"
    "epithelial_tissue_comparison_final_fullgene.h5ad"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/h2048/data/py/20260531/main_lineage_h5ad_cleanup_20260531/epithelial/figures"
)
DEFAULT_FORMATS = ("png", "pdf")
LABEL_COLUMN_HINTS = (
    "ann_",
    "cell_type",
    "label",
    "scanvi",
    "annotation",
    "pred",
    "class",
    "cluster",
    "ident",
)


@dataclass(frozen=True)
class PlotSpec:
    obs_column: str
    label: str
    filename_stub: str
    basis_candidates: tuple[str, ...]


PLOT_SPECS: tuple[PlotSpec, ...] = (
    PlotSpec(
        obs_column="tissue",
        label="Tissue",
        filename_stub="tissue",
        basis_candidates=("X_umap_scanvi", "X_umap"),
    ),
    PlotSpec(
        obs_column="cell_type_L2",
        label="Cell type L2",
        filename_stub="cell_type_l2",
        basis_candidates=("X_umap_major", "X_umap_scanvi", "X_umap"),
    ),
    PlotSpec(
        obs_column="cell_type_scanvi_pred",
        label="Cell type L3",
        filename_stub="cell_type_l3",
        basis_candidates=("X_umap_fine", "X_umap_scanvi", "X_umap"),
    ),
)


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize cleaned epithelial h5ad with boundary UMAPs and count panels.")
    parser.add_argument("--input-h5ad", default=str(DEFAULT_INPUT_H5AD))
    parser.add_argument("--reference-h5ad", default=str(DEFAULT_REFERENCE_H5AD))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--formats", default=",".join(DEFAULT_FORMATS), help="Comma-separated output formats, e.g. png,pdf")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--background-point-size", type=float, default=0.35)
    parser.add_argument("--highlight-point-size", type=float, default=1.0)
    parser.add_argument("--background-alpha", type=float, default=0.70)
    parser.add_argument("--highlight-alpha", type=float, default=0.95)
    parser.add_argument("--legend-ncol", type=int, default=1)
    parser.add_argument("--max-panel-background-points", type=int, default=80000)
    parser.add_argument("--max-highlight-points", type=int, default=25000)
    parser.add_argument("--max-hull-points", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def first_existing(candidates: Iterable[str], available: Iterable[str]) -> str:
    available_set = set(available)
    for key in candidates:
        if key in available_set:
            return key
    raise KeyError(f"None of the candidates exist: {list(candidates)}")


def normalize_labels(series: pd.Series) -> pd.Series:
    labels = series.astype("string").fillna("<NA>").str.strip()
    labels = labels.replace("", "<NA>")
    return labels


def sample_indices(n_total: int, max_points: int, rng: np.random.Generator) -> np.ndarray:
    if max_points <= 0 or n_total <= max_points:
        return np.arange(n_total)
    return np.sort(rng.choice(n_total, size=max_points, replace=False))


def darken(color: str, factor: float = 0.65) -> str:
    rgb = np.array(to_rgb(color))
    return to_hex(np.clip(rgb * factor, 0, 1))


def create_palette(categories: list[str]) -> dict[str, str]:
    n = len(categories)
    if n <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [to_hex(cmap(i)) for i in range(n)]
    elif n <= 40:
        cmap1 = plt.get_cmap("tab20")
        cmap2 = plt.get_cmap("tab20b")
        colors = [to_hex(cmap1(i)) for i in range(20)] + [to_hex(cmap2(i)) for i in range(n - 20)]
    else:
        colors = [to_hex(plt.get_cmap("hsv")(i / max(n, 1))) for i in range(n)]
    return dict(zip(categories, colors))


def select_label_columns(columns: Iterable[str]) -> list[str]:
    selected = [col for col in columns if any(hint in col.lower() for hint in LABEL_COLUMN_HINTS)]
    return selected


def load_selected_obs(input_h5ad: Path, columns: list[str]) -> pd.DataFrame:
    backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        if not columns:
            return pd.DataFrame(index=backed.obs_names.copy())
        return backed.obs[columns].copy()
    finally:
        file_manager = getattr(backed, "file", None)
        if file_manager is not None:
            file_manager.close()


def get_selected_label_columns(input_h5ad: Path) -> list[str]:
    backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        return select_label_columns(backed.obs.columns)
    finally:
        file_manager = getattr(backed, "file", None)
        if file_manager is not None:
            file_manager.close()


def build_special_label_audit(dataset_obs_map: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_name, obs in dataset_obs_map.items():
        for column in obs.columns:
            series = obs[column].astype("string").fillna("").str.strip()
            match_specs = {
                "contains_pnec": series.str.contains("PNEC", case=False, regex=False),
                "contains_neuro": series.str.contains("neuro", case=False, regex=False),
                "exact_l3": series.str.fullmatch("L3", case=False),
            }
            for match_type, mask in match_specs.items():
                n = int(mask.sum())
                top_values: dict[str, int] = {}
                if n:
                    counts = series[mask].value_counts(dropna=False).head(10)
                    top_values = {str(key): int(value) for key, value in counts.items()}
                rows.append(
                    {
                        "dataset_name": dataset_name,
                        "column": column,
                        "match_type": match_type,
                        "match_count": n,
                        "top_values_json": json.dumps(top_values, ensure_ascii=False, sort_keys=True),
                    }
                )
    return pd.DataFrame(rows)


def summarize_special_label_audit(audit_df: pd.DataFrame, dataset_columns: dict[str, list[str]]) -> dict[str, object]:
    summary_rows: list[dict[str, object]] = []
    for dataset_name, columns in dataset_columns.items():
        dataset_df = audit_df[audit_df["dataset_name"] == dataset_name].copy()
        counts = {
            match_type: int(dataset_df.loc[dataset_df["match_type"] == match_type, "match_count"].sum())
            for match_type in ("contains_pnec", "contains_neuro", "exact_l3")
        }
        pnec_hit_columns = (
            dataset_df[(dataset_df["match_type"] == "contains_pnec") & (dataset_df["match_count"] > 0)]["column"]
            .drop_duplicates()
            .tolist()
        )
        neuro_hit_columns = (
            dataset_df[(dataset_df["match_type"] == "contains_neuro") & (dataset_df["match_count"] > 0)]["column"]
            .drop_duplicates()
            .tolist()
        )
        l3_hit_columns = (
            dataset_df[(dataset_df["match_type"] == "exact_l3") & (dataset_df["match_count"] > 0)]["column"]
            .drop_duplicates()
            .tolist()
        )
        summary_rows.append(
            {
                "dataset_name": dataset_name,
                "n_label_columns_scanned": int(len(columns)),
                "contains_pnec_total": counts["contains_pnec"],
                "contains_neuro_total": counts["contains_neuro"],
                "exact_l3_total": counts["exact_l3"],
                "pnec_hit_columns": pnec_hit_columns,
                "neuro_hit_columns": neuro_hit_columns,
                "l3_hit_columns": l3_hit_columns,
            }
        )
    return {"datasets": summary_rows}


def build_special_label_note(audit_summary: dict[str, object]) -> str:
    rows = audit_summary.get("datasets", [])
    if not rows:
        return "PNEC audit: no label columns scanned"
    total_pnec = sum(int(row["contains_pnec_total"]) for row in rows)
    total_neuro = sum(int(row["contains_neuro_total"]) for row in rows)
    parts = []
    for row in rows:
        dataset_name = str(row["dataset_name"])
        parts.append(
            f"{dataset_name}: PNEC={int(row['contains_pnec_total'])}, "
            f"neuro={int(row['contains_neuro_total'])}, exact L3={int(row['exact_l3_total'])}"
        )
    prefix = "Special-label audit"
    if total_pnec == 0:
        prefix += " | no literal PNEC string"
    if total_neuro > 0:
        prefix += " | neuroendocrine labels present"
    return prefix + " | " + " | ".join(parts)


def add_corner_note(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.01,
        0.01,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="#303030",
        bbox={"facecolor": "white", "edgecolor": "#bdbdbd", "alpha": 0.92, "boxstyle": "round,pad=0.35"},
    )


def collect_special_label_groups(obs: pd.DataFrame) -> list[dict[str, object]]:
    groups: dict[bytes, dict[str, object]] = {}
    for column in obs.columns:
        series = obs[column].astype("string").fillna("").str.strip()
        mask = (
            series.str.contains("PNEC", case=False, regex=False)
            | series.str.contains("neuro", case=False, regex=False)
        ).to_numpy()
        n = int(mask.sum())
        if n == 0:
            continue
        key = np.packbits(mask.astype(np.uint8)).tobytes()
        top_values = series[mask].value_counts(dropna=False).head(5)
        group = groups.setdefault(
            key,
            {
                "mask": mask,
                "n_cells": n,
                "columns": [],
                "top_values": {},
            },
        )
        group["columns"].append(column)
        for label, value in top_values.items():
            label_str = str(label)
            group["top_values"][label_str] = max(group["top_values"].get(label_str, 0), int(value))

    ordered_groups = sorted(groups.values(), key=lambda item: (-int(item["n_cells"]), ",".join(item["columns"])))
    for idx, group in enumerate(ordered_groups, start=1):
        top_value_items = sorted(group["top_values"].items(), key=lambda kv: (-kv[1], kv[0]))
        label_display = "; ".join(f"{label}={count}" for label, count in top_value_items[:2])
        group["group_id"] = f"G{idx}"
        group["label_display"] = label_display
    return ordered_groups


def make_special_label_group_df(groups: list[dict[str, object]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group in groups:
        rows.append(
            {
                "group_id": str(group["group_id"]),
                "n_cells": int(group["n_cells"]),
                "label_display": str(group["label_display"]),
                "source_columns": json.dumps(group["columns"], ensure_ascii=False),
                "top_values_json": json.dumps(group["top_values"], ensure_ascii=False, sort_keys=True),
            }
        )
    return pd.DataFrame(rows)


def plot_special_label_panels(
    coords: np.ndarray,
    groups: list[dict[str, object]],
    title: str,
    basis_name: str,
    output_stem: Path,
    formats: list[str],
    dpi: int,
    background_size: float,
    highlight_size: float,
    background_alpha: float,
    highlight_alpha: float,
    max_background_points: int,
    max_highlight_points: int,
    max_hull_points: int,
    rng: np.random.Generator,
    figure_note: str | None = None,
) -> list[str]:
    if not groups:
        return []
    n_groups = len(groups)
    ncols = min(3, n_groups)
    nrows = int(math.ceil(n_groups / max(1, ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.5 * nrows))
    axes_array = np.atleast_1d(axes).reshape(nrows, ncols)

    background_idx = sample_indices(coords.shape[0], max_background_points, rng)
    background_coords = coords[background_idx]
    highlight_color = "#c51b7d"

    for ax in axes_array.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_aspect("equal")

    for ax, group in zip(axes_array.flat, groups):
        highlight_coords = coords[np.asarray(group["mask"], dtype=bool)]
        if highlight_coords.shape[0] > max_highlight_points:
            idx = sample_indices(highlight_coords.shape[0], max_highlight_points, rng)
            highlight_coords = highlight_coords[idx]

        background = ax.scatter(
            background_coords[:, 0],
            background_coords[:, 1],
            c="#d9d9d9",
            s=background_size,
            alpha=background_alpha,
            linewidths=0,
        )
        background.set_rasterized(True)
        if highlight_coords.size > 0:
            highlight = ax.scatter(
                highlight_coords[:, 0],
                highlight_coords[:, 1],
                c=highlight_color,
                s=highlight_size,
                alpha=highlight_alpha,
                linewidths=0,
            )
            highlight.set_rasterized(True)
            draw_density_boundary(ax=ax, points=highlight_coords, color=highlight_color, max_points=max_hull_points, rng=rng, linewidth=1.4)

        source_columns = list(group["columns"])
        source_preview = "; ".join(source_columns[:2])
        if len(source_columns) > 2:
            source_preview += f"; +{len(source_columns) - 2} more"
        ax.set_title(
            f"{group['label_display']}\nn={int(group['n_cells']):,} | {source_preview}",
            fontsize=10,
            fontweight="bold",
        )

    for ax in axes_array.flat[n_groups:]:
        ax.axis("off")

    fig.suptitle(f"{title}\nBasis: {basis_name}", fontsize=17, fontweight="bold", y=0.995)
    if figure_note:
        fig.text(
            0.012,
            0.012,
            figure_note,
            ha="left",
            va="bottom",
            fontsize=9,
            color="#303030",
            bbox={"facecolor": "white", "edgecolor": "#bdbdbd", "alpha": 0.92, "boxstyle": "round,pad=0.35"},
        )
    fig.tight_layout(rect=(0, 0.03 if figure_note else 0, 1, 0.98))
    return save_figure(fig, output_stem, formats, dpi)


def build_component_masks(hist: np.ndarray, level: float) -> list[np.ndarray]:
    dense_mask = hist >= level
    if not np.any(dense_mask):
        return []
    if binary_opening is not None and binary_closing is not None:
        structure = np.ones((3, 3), dtype=bool)
        dense_mask = binary_opening(dense_mask, structure=structure, iterations=1)
        dense_mask = binary_closing(dense_mask, structure=structure, iterations=1)
        if binary_fill_holes is not None:
            dense_mask = binary_fill_holes(dense_mask)
    if not np.any(dense_mask):
        return []

    if ndi_label is None:
        return [dense_mask]

    labeled, n_components = ndi_label(dense_mask)
    component_masks: list[np.ndarray] = []
    for component_id in range(1, n_components + 1):
        component_mask = labeled == component_id
        if int(component_mask.sum()) < 6:
            continue
        component_masks.append(component_mask)
    return component_masks


def compute_density_hist(
    points: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    grid_size: int,
    max_points: int,
    rng: np.random.Generator,
) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.shape[0] == 0:
        return np.zeros((grid_size, grid_size), dtype=float)
    if points.shape[0] > max_points:
        idx = sample_indices(points.shape[0], max_points, rng)
        points = points[idx]

    hist, _, _ = np.histogram2d(
        points[:, 0],
        points[:, 1],
        bins=grid_size,
        range=[list(x_range), list(y_range)],
    )
    if gaussian_filter is not None:
        hist = gaussian_filter(hist, sigma=1.2)
    return hist


def choose_density_level(hist: np.ndarray, n_points: int) -> float | None:
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


def extract_contour_paths(
    xcenters: np.ndarray,
    ycenters: np.ndarray,
    contour_field: np.ndarray,
    level: float,
) -> list[np.ndarray]:
    tmp_fig, tmp_ax = plt.subplots(figsize=(2, 2))
    try:
        contour = tmp_ax.contour(
            xcenters,
            ycenters,
            contour_field.T,
            levels=[level],
            colors=["black"],
            linewidths=1.0,
        )
        paths = [segment for segment in contour.allsegs[0] if len(segment) >= 2]
    finally:
        plt.close(tmp_fig)
    return paths


def nearest_winner_label(
    midpoint: np.ndarray,
    winner: np.ndarray,
    occupied_mask: np.ndarray,
    xcenters: np.ndarray,
    ycenters: np.ndarray,
    categories: list[str],
) -> str | None:
    ix = int(np.clip(np.searchsorted(xcenters, midpoint[0]), 0, len(xcenters) - 1))
    iy = int(np.clip(np.searchsorted(ycenters, midpoint[1]), 0, len(ycenters) - 1))

    for radius in range(0, 5):
        x0 = max(ix - radius, 0)
        x1 = min(ix + radius + 1, len(xcenters))
        y0 = max(iy - radius, 0)
        y1 = min(iy + radius + 1, len(ycenters))
        local_mask = occupied_mask[x0:x1, y0:y1]
        if not np.any(local_mask):
            continue

        local_x, local_y = np.where(local_mask)
        grid_x = xcenters[x0:x1][local_x]
        grid_y = ycenters[y0:y1][local_y]
        distances = (grid_x - midpoint[0]) ** 2 + (grid_y - midpoint[1]) ** 2
        best = int(np.argmin(distances))
        winner_idx = int(winner[x0 + local_x[best], y0 + local_y[best]])
        if 0 <= winner_idx < len(categories):
            return categories[winner_idx]
    return None


def prepare_outer_edge_segments(
    coords: np.ndarray,
    labels: pd.Series,
    categories: list[str],
    max_points: int,
    rng: np.random.Generator,
) -> tuple[list[np.ndarray], list[str]]:
    coords = np.asarray(coords, dtype=float)
    labels_array = labels.to_numpy()
    if coords.shape[0] == 0:
        return [], []

    x = coords[:, 0]
    y = coords[:, 1]
    xr = float(x.max() - x.min())
    yr = float(y.max() - y.min())
    if xr <= 0 or yr <= 0:
        return [], []

    x_pad = xr * 0.03
    y_pad = yr * 0.03
    x_range = (float(x.min() - x_pad), float(x.max() + x_pad))
    y_range = (float(y.min() - y_pad), float(y.max() + y_pad))
    grid_size = 260 if coords.shape[0] >= 100000 else 220 if coords.shape[0] >= 30000 else 180

    density_stack: list[np.ndarray] = []
    active_masks: list[np.ndarray] = []
    for label in categories:
        points = coords[labels_array == label]
        hist = compute_density_hist(
            points=points,
            x_range=x_range,
            y_range=y_range,
            grid_size=grid_size,
            max_points=max_points,
            rng=rng,
        )
        level = choose_density_level(hist, int(points.shape[0]))
        if level is None:
            active_mask = np.zeros_like(hist, dtype=bool)
        else:
            component_masks = build_component_masks(hist, level)
            if component_masks:
                active_mask = np.zeros_like(hist, dtype=bool)
                for component_mask in component_masks:
                    active_mask |= component_mask
            else:
                active_mask = hist >= level
        density_stack.append(hist)
        active_masks.append(active_mask)

    occupied_mask = np.zeros_like(density_stack[0], dtype=bool)
    for active_mask in active_masks:
        occupied_mask |= active_mask

    if binary_closing is not None:
        occupied_mask = binary_closing(occupied_mask, structure=np.ones((3, 3), dtype=bool), iterations=1)
    if binary_fill_holes is not None:
        occupied_mask = binary_fill_holes(occupied_mask)

    score_stack = np.stack(
        [np.where(active_mask, hist, -np.inf) for hist, active_mask in zip(density_stack, active_masks)],
        axis=0,
    )
    winner = np.argmax(score_stack, axis=0)
    xedges = np.linspace(x_range[0], x_range[1], grid_size + 1)
    yedges = np.linspace(y_range[0], y_range[1], grid_size + 1)
    xcenters = (xedges[:-1] + xedges[1:]) / 2.0
    ycenters = (yedges[:-1] + yedges[1:]) / 2.0

    contour_field = occupied_mask.astype(float)
    if gaussian_filter is not None:
        contour_field = gaussian_filter(contour_field, sigma=1.0)

    paths = extract_contour_paths(xcenters=xcenters, ycenters=ycenters, contour_field=contour_field, level=0.5)
    segments: list[np.ndarray] = []
    segment_labels: list[str] = []
    for path in paths:
        for start, end in zip(path[:-1], path[1:]):
            midpoint = (start + end) / 2.0
            label = nearest_winner_label(
                midpoint=midpoint,
                winner=winner,
                occupied_mask=occupied_mask,
                xcenters=xcenters,
                ycenters=ycenters,
                categories=categories,
            )
            if label is None:
                continue
            segments.append(np.vstack([start, end]))
            segment_labels.append(label)
    return segments, segment_labels


def draw_density_boundary(
    ax: plt.Axes,
    points: np.ndarray,
    color: str,
    max_points: int,
    rng: np.random.Generator,
    linewidth: float = 1.4,
) -> None:
    points = np.asarray(points, dtype=float)
    if points.shape[0] < 12:
        return
    if points.shape[0] > max_points:
        idx = sample_indices(points.shape[0], max_points, rng)
        points = points[idx]

    x = points[:, 0]
    y = points[:, 1]
    xr = float(x.max() - x.min())
    yr = float(y.max() - y.min())
    if xr <= 0 or yr <= 0:
        return

    x_pad = xr * 0.03
    y_pad = yr * 0.03
    grid_size = 200 if points.shape[0] >= 5000 else 140
    hist, xedges, yedges = np.histogram2d(
        x,
        y,
        bins=grid_size,
        range=[[x.min() - x_pad, x.max() + x_pad], [y.min() - y_pad, y.max() + y_pad]],
    )
    if gaussian_filter is not None:
        hist = gaussian_filter(hist, sigma=1.2)

    nonzero = hist[hist > 0]
    if nonzero.size == 0:
        return
    hist_max = float(nonzero.max())
    quantile = 0.35 if points.shape[0] >= 1000 else 0.25
    level = max(float(np.quantile(nonzero, quantile)), hist_max * 0.12)
    if level >= hist_max:
        level = hist_max * 0.50
    if level <= 0:
        return

    xcenters = (xedges[:-1] + xedges[1:]) / 2.0
    ycenters = (yedges[:-1] + yedges[1:]) / 2.0
    component_masks = build_component_masks(hist, level)
    if not component_masks:
        component_masks = [hist >= level]

    for component_mask in component_masks:
        contour_field = component_mask.astype(float)
        if gaussian_filter is not None:
            contour_field = gaussian_filter(contour_field, sigma=0.8)
        ax.contour(
            xcenters,
            ycenters,
            contour_field.T,
            levels=[0.5],
            colors=[darken(color, 0.72)],
            linewidths=linewidth,
            alpha=0.98,
        )


def save_figure(fig: plt.Figure, stem: Path, formats: list[str], dpi: int) -> list[str]:
    outputs: list[str] = []
    for fmt in formats:
        path = stem.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def make_counts_df(labels: pd.Series) -> pd.DataFrame:
    counts = labels.value_counts(dropna=False)
    total = int(counts.sum())
    df = counts.rename_axis("label").reset_index(name="count")
    df["percentage"] = df["count"] / total * 100.0
    return df


def plot_boundary_umap(
    coords: np.ndarray,
    labels: pd.Series,
    counts_df: pd.DataFrame,
    palette: dict[str, str],
    title: str,
    basis_name: str,
    output_stem: Path,
    formats: list[str],
    dpi: int,
    point_size: float,
    alpha: float,
    max_hull_points: int,
    legend_ncol: int,
    rng: np.random.Generator,
    corner_note: str | None = None,
) -> list[str]:
    categories = counts_df["label"].tolist()
    cat = pd.Categorical(labels, categories=categories, ordered=True)
    color_values = np.array([palette.get(label, "#bdbdbd") for label in cat], dtype=object)

    fig, ax = plt.subplots(figsize=(14, 12))
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=color_values,
        s=point_size,
        alpha=alpha,
        linewidths=0,
    )
    scatter.set_rasterized(True)
    edge_segments, edge_segment_labels = prepare_outer_edge_segments(
        coords=coords,
        labels=labels,
        categories=categories,
        max_points=max_hull_points,
        rng=rng,
    )

    if edge_segments:
        edge_colors = [darken(palette[label], 0.72) for label in edge_segment_labels]
        edge_collection = LineCollection(
            edge_segments,
            colors=edge_colors,
            linewidths=max(1.2, point_size * 1.8),
            alpha=0.98,
            capstyle="round",
            joinstyle="round",
            zorder=4,
        )
        ax.add_collection(edge_collection)

    handles: list[Line2D] = []
    for row in counts_df.itertuples(index=False):
        label = str(row.label)
        color = palette[label]
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color=darken(color, 0.72),
                markerfacecolor=color,
                markeredgewidth=0,
                markersize=7,
                label=f"{label} (n={int(row.count):,}, {row.percentage:.1f}%)",
            )
        )

    ax.set_title(f"{title} boundary UMAP\nBasis: {basis_name}", fontsize=16, fontweight="bold")
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=9,
        ncol=max(1, legend_ncol),
    )
    if corner_note:
        add_corner_note(ax, corner_note)
    return save_figure(fig, output_stem, formats, dpi)


def plot_highlight_panels(
    coords: np.ndarray,
    labels: pd.Series,
    counts_df: pd.DataFrame,
    palette: dict[str, str],
    title: str,
    basis_name: str,
    output_stem: Path,
    formats: list[str],
    dpi: int,
    background_size: float,
    highlight_size: float,
    background_alpha: float,
    highlight_alpha: float,
    max_background_points: int,
    max_highlight_points: int,
    max_hull_points: int,
    rng: np.random.Generator,
    figure_note: str | None = None,
) -> list[str]:
    n_categories = len(counts_df)
    ncols = 4 if n_categories > 4 else n_categories
    nrows = int(math.ceil(n_categories / max(1, ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 4.3 * nrows))
    axes_array = np.atleast_1d(axes).reshape(nrows, ncols)

    background_idx = sample_indices(coords.shape[0], max_background_points, rng)
    background_coords = coords[background_idx]
    labels_array = labels.to_numpy()
    total = int(counts_df["count"].sum())

    for ax in axes_array.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_aspect("equal")

    for ax, row in zip(axes_array.flat, counts_df.itertuples(index=False)):
        label = str(row.label)
        color = palette[label]
        mask = labels_array == label
        highlight_coords = coords[mask]
        if highlight_coords.shape[0] > max_highlight_points:
            idx = sample_indices(highlight_coords.shape[0], max_highlight_points, rng)
            highlight_coords = highlight_coords[idx]

        background = ax.scatter(
            background_coords[:, 0],
            background_coords[:, 1],
            c="#d9d9d9",
            s=background_size,
            alpha=background_alpha,
            linewidths=0,
        )
        background.set_rasterized(True)
        if highlight_coords.size > 0:
            highlight = ax.scatter(
                highlight_coords[:, 0],
                highlight_coords[:, 1],
                c=color,
                s=highlight_size,
                alpha=highlight_alpha,
                linewidths=0,
            )
            highlight.set_rasterized(True)
            draw_density_boundary(ax=ax, points=highlight_coords, color=color, max_points=max_hull_points, rng=rng, linewidth=1.3)
        ax.set_title(f"{label}\nn={int(row.count):,} ({float(row.count) / total * 100:.1f}%)", fontsize=11, fontweight="bold")

    for ax in axes_array.flat[n_categories:]:
        ax.axis("off")

    fig.suptitle(f"{title} highlight UMAP panels with counts\nBasis: {basis_name}", fontsize=18, fontweight="bold", y=0.995)
    if figure_note:
        fig.text(
            0.012,
            0.012,
            figure_note,
            ha="left",
            va="bottom",
            fontsize=9,
            color="#303030",
            bbox={"facecolor": "white", "edgecolor": "#bdbdbd", "alpha": 0.92, "boxstyle": "round,pad=0.35"},
        )
    fig.tight_layout(rect=(0, 0.025 if figure_note else 0, 1, 0.98))
    return save_figure(fig, output_stem, formats, dpi)


def plot_count_barplot(
    counts_df: pd.DataFrame,
    title: str,
    output_stem: Path,
    formats: list[str],
    dpi: int,
    palette: dict[str, str],
) -> list[str]:
    df = counts_df.sort_values("count", ascending=True).reset_index(drop=True)
    fig_height = max(4.5, 0.42 * len(df) + 2.0)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    colors = [palette[str(label)] for label in df["label"]]
    bars = ax.barh(df["label"], df["count"], color=colors, edgecolor="none", alpha=0.9)
    ax.set_title(f"{title} counts", fontsize=16, fontweight="bold")
    ax.set_xlabel("Number of cells")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.22)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    max_count = float(df["count"].max()) if len(df) else 0.0
    x_offset = max_count * 0.01 if max_count > 0 else 1.0
    for bar, row in zip(bars, df.itertuples(index=False)):
        ax.text(
            bar.get_width() + x_offset,
            bar.get_y() + bar.get_height() / 2.0,
            f"{int(row.count):,} ({row.percentage:.1f}%)",
            va="center",
            ha="left",
            fontsize=9,
        )

    fig.tight_layout()
    return save_figure(fig, output_stem, formats, dpi)


def load_minimal_frame(input_h5ad: Path, specs: tuple[PlotSpec, ...]) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, str]]:
    backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        available_bases = list(backed.obsm.keys())
        chosen_bases = {spec.obs_column: first_existing(spec.basis_candidates, available_bases) for spec in specs}
        obs_cols = [spec.obs_column for spec in specs]
        obs = backed.obs[obs_cols].copy()
        coords_map: dict[str, np.ndarray] = {}
        for basis in sorted(set(chosen_bases.values())):
            coords_map[basis] = np.asarray(backed.obsm[basis], dtype=np.float32)
    finally:
        file_manager = getattr(backed, "file", None)
        if file_manager is not None:
            file_manager.close()
    return obs, coords_map, chosen_bases


def main() -> None:
    args = parse_args()
    input_h5ad = Path(args.input_h5ad)
    reference_h5ad = Path(args.reference_h5ad) if str(args.reference_h5ad).strip() else None
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    formats = split_csv(args.formats)
    rng = np.random.default_rng(args.seed)

    obs, coords_map, chosen_bases = load_minimal_frame(input_h5ad, PLOT_SPECS)

    clean_label_columns = get_selected_label_columns(input_h5ad)
    cleaned_label_obs = load_selected_obs(input_h5ad, clean_label_columns)
    label_dataset_obs: dict[str, pd.DataFrame] = {"cleaned": cleaned_label_obs}
    dataset_columns: dict[str, list[str]] = {"cleaned": clean_label_columns}
    if reference_h5ad is not None and reference_h5ad.exists():
        reference_label_columns = get_selected_label_columns(reference_h5ad)
        label_dataset_obs["reference"] = load_selected_obs(reference_h5ad, reference_label_columns)
        dataset_columns["reference"] = reference_label_columns

    audit_df = build_special_label_audit(label_dataset_obs)
    audit_summary = summarize_special_label_audit(audit_df, dataset_columns)
    special_label_note = build_special_label_note(audit_summary)
    special_label_audit_tsv = output_dir / "epithelial_cleaned_special_label_audit_20260531.tsv"
    special_label_audit_json = output_dir / "epithelial_cleaned_special_label_audit_20260531.json"
    audit_df.to_csv(special_label_audit_tsv, sep="\t", index=False)
    special_label_audit_json.write_text(json.dumps(audit_summary, indent=2, ensure_ascii=False), encoding="utf-8")
    special_label_groups = collect_special_label_groups(cleaned_label_obs)
    special_label_group_tsv = output_dir / "epithelial_cleaned_special_label_groups_20260531.tsv"
    special_label_group_df = make_special_label_group_df(special_label_groups)
    special_label_group_df.to_csv(special_label_group_tsv, sep="\t", index=False)
    special_label_basis = first_existing(("X_umap_fine", "X_umap_scanvi", "X_umap"), coords_map.keys())
    special_label_panel_outputs = plot_special_label_panels(
        coords=coords_map[special_label_basis],
        groups=special_label_groups,
        title="PNEC / neuroendocrine-related label panels",
        basis_name=special_label_basis,
        output_stem=output_dir / "epithelial_cleaned_special_label_panels",
        formats=formats,
        dpi=args.dpi,
        background_size=args.background_point_size,
        highlight_size=args.highlight_point_size,
        background_alpha=args.background_alpha,
        highlight_alpha=args.highlight_alpha,
        max_background_points=args.max_panel_background_points,
        max_highlight_points=args.max_highlight_points,
        max_hull_points=args.max_hull_points,
        rng=rng,
        figure_note=special_label_note,
    )

    manifest_rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "input_h5ad": str(input_h5ad),
        "reference_h5ad": str(reference_h5ad) if reference_h5ad is not None and reference_h5ad.exists() else None,
        "output_dir": str(output_dir),
        "n_obs": int(obs.shape[0]),
        "plot_specs": {},
        "special_label_audit": {
            "audit_tsv": str(special_label_audit_tsv),
            "audit_json": str(special_label_audit_json),
            "group_tsv": str(special_label_group_tsv),
            "panel_basis": special_label_basis,
            "panel_outputs": special_label_panel_outputs,
            "summary": audit_summary,
        },
    }

    manifest_rows.extend(
        [
            {
                "obs_column": "special_label_audit",
                "label": "Special label audit",
                "basis": "NA",
                "artifact_type": "special_label_audit_tsv",
                "path": str(special_label_audit_tsv),
            },
            {
                "obs_column": "special_label_audit",
                "label": "Special label audit",
                "basis": "NA",
                "artifact_type": "special_label_audit_json",
                "path": str(special_label_audit_json),
            },
            {
                "obs_column": "special_label_audit",
                "label": "Special label audit",
                "basis": special_label_basis,
                "artifact_type": "special_label_group_tsv",
                "path": str(special_label_group_tsv),
            },
            *[
                {
                    "obs_column": "special_label_audit",
                    "label": "Special label audit",
                    "basis": special_label_basis,
                    "artifact_type": "special_label_panels",
                    "path": path,
                }
                for path in special_label_panel_outputs
            ],
        ]
    )

    for spec in PLOT_SPECS:
        labels = normalize_labels(obs[spec.obs_column])
        counts_df = make_counts_df(labels)
        counts_path = output_dir / f"epithelial_cleaned_{spec.filename_stub}_counts.tsv"
        counts_df.to_csv(counts_path, sep="\t", index=False)

        categories = counts_df["label"].tolist()
        palette = create_palette(categories)
        basis_name = chosen_bases[spec.obs_column]
        coords = coords_map[basis_name]

        boundary_outputs = plot_boundary_umap(
            coords=coords,
            labels=labels,
            counts_df=counts_df,
            palette=palette,
            title=spec.label,
            basis_name=basis_name,
            output_stem=output_dir / f"epithelial_cleaned_{spec.filename_stub}_boundary_umap",
            formats=formats,
            dpi=args.dpi,
            point_size=args.point_size,
            alpha=args.highlight_alpha,
            max_hull_points=args.max_hull_points,
            legend_ncol=args.legend_ncol,
            rng=rng,
            corner_note=special_label_note if spec.obs_column == "cell_type_scanvi_pred" else None,
        )
        panel_outputs = plot_highlight_panels(
            coords=coords,
            labels=labels,
            counts_df=counts_df,
            palette=palette,
            title=spec.label,
            basis_name=basis_name,
            output_stem=output_dir / f"epithelial_cleaned_{spec.filename_stub}_highlight_panels",
            formats=formats,
            dpi=args.dpi,
            background_size=args.background_point_size,
            highlight_size=args.highlight_point_size,
            background_alpha=args.background_alpha,
            highlight_alpha=args.highlight_alpha,
            max_background_points=args.max_panel_background_points,
            max_highlight_points=args.max_highlight_points,
            max_hull_points=args.max_hull_points,
            rng=rng,
            figure_note=special_label_note if spec.obs_column == "cell_type_scanvi_pred" else None,
        )
        barplot_outputs = plot_count_barplot(
            counts_df=counts_df,
            title=spec.label,
            output_stem=output_dir / f"epithelial_cleaned_{spec.filename_stub}_counts_barplot",
            formats=formats,
            dpi=args.dpi,
            palette=palette,
        )

        summary["plot_specs"][spec.obs_column] = {
            "label": spec.label,
            "basis": basis_name,
            "n_categories": int(len(counts_df)),
            "counts_tsv": str(counts_path),
            "boundary_outputs": boundary_outputs,
            "panel_outputs": panel_outputs,
            "barplot_outputs": barplot_outputs,
        }

        manifest_rows.extend(
            [
                {
                    "obs_column": spec.obs_column,
                    "label": spec.label,
                    "basis": basis_name,
                    "artifact_type": "counts_tsv",
                    "path": str(counts_path),
                },
                *[
                    {
                        "obs_column": spec.obs_column,
                        "label": spec.label,
                        "basis": basis_name,
                        "artifact_type": "boundary_umap",
                        "path": path,
                    }
                    for path in boundary_outputs
                ],
                *[
                    {
                        "obs_column": spec.obs_column,
                        "label": spec.label,
                        "basis": basis_name,
                        "artifact_type": "highlight_panels",
                        "path": path,
                    }
                    for path in panel_outputs
                ],
                *[
                    {
                        "obs_column": spec.obs_column,
                        "label": spec.label,
                        "basis": basis_name,
                        "artifact_type": "counts_barplot",
                        "path": path,
                    }
                    for path in barplot_outputs
                ],
            ]
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "epithelial_cleaned_visual_manifest_20260531.tsv"
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    summary["manifest_tsv"] = str(manifest_path)

    summary_path = output_dir / "epithelial_cleaned_visual_summary_20260531.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[viz] wrote manifest: {manifest_path}")
    print(f"[viz] wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
