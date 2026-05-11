#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Epithelial-specific analysis helper.

This module extracts reusable epithelial L3 annotation + visualization logic from
`epithelial_scvi_scanvi_20260308_v2_6.ipynb` while delegating generic expression
access and plotting utilities to `anndata_expression_viz_helper_20260423_v1.py`.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd


def _is_categorical_series(series: pd.Series) -> bool:
    return isinstance(series.dtype, pd.CategoricalDtype)

try:
    from anndata_expression_viz_helper_20260423_v1 import (
        build_mean_expression_matrix,
        build_review_dotplot_table,
        compute_dynamic_figsize,
        deduplicate_preserve_order,
        ensure_categorical_with_order,
        extract_gene_vector,
        gene_available,
        get_available_markers,
        normalize_string_series,
        prepare_output_dirs,
        render_scanpy_dotplot,
        render_seaborn_heatmap,
        resolve_expression_source,
        safe_obs_column,
        save_current_figure,
        write_tsv,
    )
except ImportError:
    _HELPER_PATH = Path(__file__).resolve().with_name("anndata_expression_viz_helper_20260423_v1.py")
    _HELPER_SPEC = importlib.util.spec_from_file_location(
        "anndata_expression_viz_helper_20260423_v1",
        _HELPER_PATH,
    )
    if _HELPER_SPEC is None or _HELPER_SPEC.loader is None:
        raise
    _helper = importlib.util.module_from_spec(_HELPER_SPEC)
    _HELPER_SPEC.loader.exec_module(_helper)
    build_mean_expression_matrix = _helper.build_mean_expression_matrix
    build_review_dotplot_table = _helper.build_review_dotplot_table
    compute_dynamic_figsize = _helper.compute_dynamic_figsize
    deduplicate_preserve_order = _helper.deduplicate_preserve_order
    ensure_categorical_with_order = _helper.ensure_categorical_with_order
    extract_gene_vector = _helper.extract_gene_vector
    gene_available = _helper.gene_available
    get_available_markers = _helper.get_available_markers
    normalize_string_series = _helper.normalize_string_series
    prepare_output_dirs = _helper.prepare_output_dirs
    render_scanpy_dotplot = _helper.render_scanpy_dotplot
    render_seaborn_heatmap = _helper.render_seaborn_heatmap
    resolve_expression_source = _helper.resolve_expression_source
    safe_obs_column = _helper.safe_obs_column
    save_current_figure = _helper.save_current_figure
    write_tsv = _helper.write_tsv

__all__ = [
    "build_default_l3_mapping",
    "build_default_marker_panels",
    "build_default_l3_colors",
    "build_default_lineage_groups",
    "build_misannotation_review_specs",
    "apply_l3_annotations",
    "make_overview_umap",
    "make_lineage_umaps",
    "make_lineage_dotplots",
    "make_core_marker_heatmap",
    "build_misannotation_review_table",
    "run_misannotation_review",
]


def build_default_marker_panels() -> dict[str, list[str]]:
    return {
        "AT1_Canonical": ["AGER", "HOPX", "CAV1", "AQP4", "RTKN2", "CLDN18", "EMP2"],
        "AT1_MatrixRemodeling": ["AGER", "CAV1", "SPARC", "COL4A1", "COL4A2", "SPOCK2"],
        "AT2_Canonical": ["SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NAPSA", "SLC34A2"],
        "AT2_Inflammatory_Repair": ["SFTPC", "ABCA3", "CHI3L1", "CXCL8", "SAA1", "LCN2"],
        "Epithelial_Cycling": ["MKI67", "TOP2A", "UBE2C", "BIRC5", "AURKB", "CENPA", "CCNB1"],
        "Basal_Progenitor": ["KRT5", "KRT14", "TP63", "KRT15", "KRT19", "ITGA6", "NGFR"],
        "Basal_Cycling": ["KRT14", "KRT5", "TP63", "MKI67", "TOP2A", "BIRC5"],
        "Basal_Inflammatory": ["KRT17", "CXCL8", "CXCL1", "CXCL2", "TNFAIP3", "FOS", "JUN"],
        "Basal_EMT_ECM": ["KRT14", "TP63", "NGFR", "FN1", "COL17A1", "MMP2", "VIM"],
        "Ciliated_Mature": ["FOXJ1", "TPPP3", "DNAH5", "DNAH9", "RSPH1", "RFX2", "RFX3"],
        "Ciliogenesis_Deuterosomal": ["DEUP1", "CCNO", "FOXN4", "MCIDAS", "CDC20B", "E2F7", "PLK4"],
        "Ciliated_Cycling_Immature": ["TPPP3", "RSPH1", "MKI67", "TOP2A", "FOXN4"],
        "Secretory_Club": ["SCGB1A1", "SCGB3A1", "SCGB3A2", "AGR2", "AGR3", "CYP2F1"],
        "Secretory_Club_AT2_Transitional": ["SCGB1A1", "SCGB3A1", "SFTPB", "NAPSA", "GPR116", "CLDN18"],
        "Goblet_Mucin": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
        "Goblet_Defense_DUOX2": ["DUOX2", "DUOXA2", "LCN2", "BPIFA2", "CEACAM5"],
        "SMG_Serous": ["LTF", "LYZ", "SLPI", "DMBT1", "BPIFA1", "AZGP1", "WFDC2"],
        "SMG_Duct_Secretory_Defense": ["PIGR", "SCGB3A1", "TCN1", "WFDC2", "DMBT1", "SLPI"],
        "Squamous_Metaplasia": ["SPRR1A", "SPRR2A", "SPRR2E", "IVL", "KRT6A", "KLK7", "S100A7"],
        "Ionocyte_Brush": ["FOXI1", "ASCL3", "CFTR", "ATP6V0D2", "CLCNKA", "CLCNKB", "BSND"],
        "Mesenchymal_Contaminant": ["COL4A1", "COL4A2", "LAMA1", "DCN", "COL1A1", "COL1A2"],
    }


def build_default_l3_mapping() -> dict[str, str]:
    return {
        "AT1 _0": "AT1_Canonical",
        "AT1 _1": "AT1_MatrixRemodeling",
        "AT1 _2": "AT2_Canonical",
        "AT2_0": "AT2_Canonical",
        "AT2_1": "AT2_Canonical",
        "AT2_2": "AT2_Inflammatory_Repair",
        "AT2_3": "Epithelial_Cycling",
        "Basal_0": "Basal_EMT_ECM",
        "Basal_1": "Basal_Cycling",
        "Basal_2": "Basal_Inflammatory",
        "Basal_3": "Goblet_Mucin",
        "Basal_4": "Goblet_Mucin",
        "Basal_5": "Basal_Progenitor",
        "Ciliated_0": "Ciliated_Mature",
        "Ciliated_1": "Ciliated_Mature",
        "Ciliated_2": "Ciliated_Cycling_Immature",
        "Ciliated_3": "Ciliated_Mature",
        "Ciliated_4": "Ciliated_Mature",
        "Ciliated_5": "Goblet_Mucin",
        "Deuterosomal_0": "Ciliogenesis_Deuterosomal",
        "Deuterosomal_1": "Ciliogenesis_Deuterosomal",
        "Dividing_Basal_0": "Basal_Cycling",
        "Dividing_Basal_1": "Basal_Cycling",
        "Ionocyte_n_Brush_0": "Ionocyte_Brush",
        "Ionocyte_n_Brush_1": "Ionocyte_Brush",
        "SMG_Basal_0": "Basal_EMT_ECM",
        "SMG_Basal_1": "Mesenchymal_Contaminant",
        "SMG_Basal_2": "Basal_EMT_ECM",
        "SMG_Duct_0": "SMG_Duct_Secretory_Defense",
        "SMG_Duct_1": "Squamous_Metaplasia",
        "SMG_Duct_2": "Squamous_Metaplasia",
        "SMG_Duct_3": "Squamous_Metaplasia",
        "SMG_Duct_4": "Squamous_Metaplasia",
        "SMG_Mucous_0": "Goblet_Mucin",
        "SMG_Mucous_1": "Ionocyte_Brush",
        "SMG_Serous_0": "SMG_Serous",
        "SMG_Serous_1": "SMG_Serous",
        "SMG_Serous_2": "SMG_Serous",
        "Secretory_Club_0": "Secretory_Club_AT2_Transitional",
        "Secretory_Goblet_0": "Goblet_Defense_DUOX2",
        "Secretory_Goblet_1": "Secretory_Club",
        "Secretory_Goblet_2": "Basal_Progenitor",
        "Secretory_Goblet_3": "Squamous_Metaplasia",
        "Secretory_Goblet_4": "Ciliated_Cycling_Immature",
        "Suprabasal_0": "Basal_Progenitor",
        "Suprabasal_1": "Basal_Cycling",
        "Suprabasal_2": "Squamous_Metaplasia",
        "Suprabasal_3": "Squamous_Metaplasia",
    }


def build_default_l3_colors() -> dict[str, str]:
    return {
        "AT1_Canonical": "#1f77b4",
        "AT1_MatrixRemodeling": "#aec7e8",
        "AT2_Canonical": "#17becf",
        "AT2_Inflammatory_Repair": "#9edae5",
        "Epithelial_Cycling": "#c49c94",
        "Basal_Progenitor": "#d62728",
        "Basal_Cycling": "#ff7f0e",
        "Basal_Inflammatory": "#ff9896",
        "Basal_EMT_ECM": "#ffbb78",
        "Ciliated_Mature": "#2ca02c",
        "Ciliogenesis_Deuterosomal": "#98df8a",
        "Ciliated_Cycling_Immature": "#8c564b",
        "Secretory_Club": "#9467bd",
        "Secretory_Club_AT2_Transitional": "#c5b0d5",
        "Goblet_Mucin": "#e377c2",
        "Goblet_Defense_DUOX2": "#f7b6d2",
        "SMG_Serous": "#bcbd22",
        "SMG_Duct_Secretory_Defense": "#dbdb8d",
        "Squamous_Metaplasia": "#7f7f7f",
        "Ionocyte_Brush": "#c7c7c7",
        "Mesenchymal_Contaminant": "#3f3f3f",
    }


def build_default_lineage_groups(kind: str = "umap") -> dict[str, list[str]]:
    if kind == "dotplot":
        return {
            "Alveolar": ["AT1_Canonical", "AT1_MatrixRemodeling", "AT2_Canonical", "AT2_Inflammatory_Repair", "Epithelial_Cycling"],
            "Basal": ["Basal_Progenitor", "Basal_Cycling", "Basal_Inflammatory", "Basal_EMT_ECM"],
            "Ciliated": ["Ciliated_Mature", "Ciliogenesis_Deuterosomal", "Ciliated_Cycling_Immature"],
            "Secretory_SMG": ["Secretory_Club", "Secretory_Club_AT2_Transitional", "Goblet_Mucin", "Goblet_Defense_DUOX2", "SMG_Serous", "SMG_Duct_Secretory_Defense"],
            "Special": ["Squamous_Metaplasia", "Ionocyte_Brush", "Mesenchymal_Contaminant"],
        }
    return {
        "Alveolar": ["AT1_Canonical", "AT1_MatrixRemodeling", "AT2_Canonical", "AT2_Inflammatory_Repair"],
        "Basal": ["Basal_Progenitor", "Basal_Cycling", "Basal_Inflammatory", "Basal_EMT_ECM"],
        "Ciliated": ["Ciliated_Mature", "Ciliogenesis_Deuterosomal", "Ciliated_Cycling_Immature"],
        "Secretory": ["Secretory_Club", "Secretory_Club_AT2_Transitional", "Goblet_Mucin", "Goblet_Defense_DUOX2"],
        "SMG": ["SMG_Serous", "SMG_Duct_Secretory_Defense"],
        "Other": ["Squamous_Metaplasia", "Ionocyte_Brush", "Epithelial_Cycling", "Mesenchymal_Contaminant"],
    }


def build_misannotation_review_specs() -> list[dict[str, Any]]:
    return [
        {
            "cluster": "AT1 _2",
            "wrong_label": "AT1_Canonical",
            "correct_label": "AT2_Canonical",
            "wrong_genes": ["AGER", "HOPX", "CAV1", "AQP4", "EMP2"],
            "correct_genes": ["SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NAPSA"],
            "note": "AT1-like cluster reassigned to AT2-like surfactant program",
        },
        {
            "cluster": "Ciliated_5",
            "wrong_label": "Ciliated_Mature",
            "correct_label": "Goblet_Mucin",
            "wrong_genes": ["FOXJ1", "TPPP3", "DNAH5", "DNAH9", "RSPH1"],
            "correct_genes": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
            "note": "Ciliated cluster reassigned to goblet/mucin program",
        },
        {
            "cluster": "SMG_Mucous_1",
            "wrong_label": "Goblet_Mucin",
            "correct_label": "Ionocyte_Brush",
            "wrong_genes": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
            "correct_genes": ["FOXI1", "ASCL3", "CFTR", "ATP6V0D2", "POU2F3", "TRPM5"],
            "note": "Verify ionocyte/brush identity against mucous program",
        },
        {
            "cluster": "SMG_Basal_1",
            "wrong_label": "Basal_Progenitor",
            "correct_label": "Mesenchymal_Contaminant",
            "wrong_genes": ["KRT5", "KRT14", "TP63", "KRT15", "ITGA6", "NGFR"],
            "correct_genes": ["COL4A1", "COL4A2", "LAMA1", "DCN", "COL1A1", "COL1A2"],
            "note": "Basal-looking cluster suspected to be mesenchymal contamination",
        },
        {
            "cluster": "Secretory_Goblet_2",
            "wrong_label": "Goblet_Mucin",
            "correct_label": "Basal_Progenitor",
            "wrong_genes": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
            "correct_genes": ["KRT5", "KRT14", "TP63", "KRT15", "ITGA6", "NGFR"],
            "note": "Secretory/goblet cluster reassigned to basal progenitor",
        },
        {
            "cluster": "Secretory_Goblet_3",
            "wrong_label": "Goblet_Mucin",
            "correct_label": "Squamous_Metaplasia",
            "wrong_genes": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
            "correct_genes": ["SPRR1A", "SPRR2A", "SPRR2E", "IVL", "KRT6A", "S100A7"],
            "note": "Secretory/goblet cluster reassigned to squamous metaplasia",
        },
        {
            "cluster": "Secretory_Goblet_4",
            "wrong_label": "Goblet_Mucin",
            "correct_label": "Ciliated_Cycling_Immature",
            "wrong_genes": ["SPDEF", "FOXA3", "MUC5AC", "MUC5B", "FCGBP", "TFF3"],
            "correct_genes": ["TPPP3", "RSPH1", "MKI67", "TOP2A", "FOXN4", "CDC20B"],
            "note": "Secretory/goblet cluster reassigned to immature cycling ciliated state",
        },
    ]


def apply_l3_annotations(
    adata: Any,
    mapping: Mapping[str, str] | None = None,
    subcluster_col: str = "subcluster",
    output_col: str = "cell_type_L3",
    inplace: bool = False,
) -> Any:
    work = adata if inplace else adata.copy()
    mapping = dict(mapping or build_default_l3_mapping())
    subcluster = normalize_string_series(safe_obs_column(work, subcluster_col))
    mapping_clean = {str(key).strip(): value for key, value in mapping.items()}
    mapped = subcluster.map(mapping_clean)

    unmapped_mask = mapped.isna()
    if unmapped_mask.any():
        mapped.loc[unmapped_mask] = subcluster.loc[unmapped_mask]

    categories = deduplicate_preserve_order(mapped.astype(str).tolist())
    work.obs[output_col] = ensure_categorical_with_order(mapped.astype(str), categories)
    work.uns[f"{output_col}_summary"] = {
        "input_column": subcluster_col,
        "output_column": output_col,
        "n_total": int(work.n_obs),
        "n_unmapped": int(unmapped_mask.sum()),
        "unmapped_clusters": deduplicate_preserve_order(subcluster.loc[unmapped_mask].astype(str).tolist()),
    }
    return work


def _ensure_umap(adata: Any, random_seed: int = 42) -> None:
    import scanpy as sc

    if "X_umap" in adata.obsm:
        return
    if "X_pca" not in adata.obsm:
        n_comps = max(2, min(50, adata.n_vars, adata.n_obs - 1 if adata.n_obs > 1 else 2))
        sc.pp.pca(adata, n_comps=n_comps, random_state=random_seed)
    sc.pp.neighbors(adata, random_state=random_seed)
    sc.tl.umap(adata, random_state=random_seed)


def _assign_l3_colors(adata: Any, color_map: Mapping[str, str], color_col: str = "cell_type_L3") -> None:
    labels = adata.obs[color_col]
    categories = list(labels.cat.categories) if _is_categorical_series(labels) else sorted(labels.astype(str).unique())
    fallback = list(np.linspace(0, 1, max(len(categories), 1)))
    colors = []
    for idx, label in enumerate(categories):
        if label in color_map:
            colors.append(color_map[label])
        else:
            import matplotlib.pyplot as plt

            colors.append(plt.cm.tab20(fallback[idx % len(fallback)]))
    adata.uns[f"{color_col}_colors"] = colors


def make_overview_umap(
    adata: Any,
    output_path: str | Path,
    color_col: str = "cell_type_L3",
    l3_colors: Mapping[str, str] | None = None,
    random_seed: int = 42,
    figure_dpi: int = 300,
    umap_size: float = 3.0,
    umap_alpha: float = 0.6,
    title: str = "Epithelial Cells - L3 Annotations",
) -> Path:
    import matplotlib.pyplot as plt
    import scanpy as sc

    _ensure_umap(adata, random_seed=random_seed)
    _assign_l3_colors(adata, l3_colors or build_default_l3_colors(), color_col=color_col)

    fig, ax = plt.subplots(figsize=(14, 10))
    sc.pl.umap(
        adata,
        color=color_col,
        ax=ax,
        show=False,
        legend_loc="right margin",
        legend_fontsize=8,
        frameon=False,
        size=umap_size,
        alpha=umap_alpha,
        title=title,
    )
    plt.tight_layout()
    return save_current_figure(output_path, dpi=figure_dpi, close=True, figure=fig)


def make_lineage_umaps(
    adata: Any,
    output_path: str | Path,
    lineage_groups: Mapping[str, Sequence[str]] | None = None,
    color_col: str = "cell_type_L3",
    random_seed: int = 42,
    figure_dpi: int = 300,
    umap_size: float = 3.0,
    umap_alpha: float = 0.6,
) -> Path:
    import matplotlib.pyplot as plt
    import scanpy as sc

    _ensure_umap(adata, random_seed=random_seed)
    groups = dict(lineage_groups or build_default_lineage_groups(kind="umap"))
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    temp_cols: list[str] = []
    for idx, (lineage, labels) in enumerate(groups.items()):
        temp_col = f"is_{lineage}"
        temp_cols.append(temp_col)
        mask = adata.obs[color_col].astype(str).isin(list(labels))
        adata.obs[temp_col] = "Other"
        adata.obs.loc[mask, temp_col] = adata.obs.loc[mask, color_col].astype(str)
        sc.pl.umap(
            adata,
            color=temp_col,
            ax=axes[idx],
            show=False,
            frameon=False,
            size=umap_size * 0.7,
            alpha=umap_alpha * 0.8,
            title=f"{lineage} Lineage",
            legend_loc="none" if idx < len(groups) - 1 else "right margin",
            legend_fontsize=6,
        )

    for idx in range(len(groups), len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()
    out = save_current_figure(output_path, dpi=figure_dpi, close=True, figure=fig)
    adata.obs.drop(columns=temp_cols, inplace=True, errors="ignore")
    return out


def make_lineage_dotplots(
    adata: Any,
    output_dir: str | Path,
    marker_panels: Mapping[str, Sequence[str]] | None = None,
    lineage_groups: Mapping[str, Sequence[str]] | None = None,
    color_col: str = "cell_type_L3",
    figure_format: str = "pdf",
    figure_dpi: int = 300,
) -> dict[str, Path]:
    dirs = prepare_output_dirs(output_dir)
    source = resolve_expression_source(adata)
    marker_panels = dict(marker_panels or build_default_marker_panels())
    lineage_groups = dict(lineage_groups or build_default_lineage_groups(kind="dotplot"))
    available_markers, _ = get_available_markers(marker_panels, source["gene_universe"])
    results: dict[str, Path] = {}

    for lineage, cell_types in lineage_groups.items():
        markers_for_lineage: list[str] = []
        for cell_type in cell_types:
            markers_for_lineage.extend(available_markers.get(cell_type, []))
        markers_for_lineage = deduplicate_preserve_order(markers_for_lineage)
        if not markers_for_lineage:
            continue
        adata_subset = adata[adata.obs[color_col].astype(str).isin(list(cell_types))].copy()
        if adata_subset.n_obs == 0:
            continue
        fig_w, fig_h = compute_dynamic_figsize(len(cell_types), len(markers_for_lineage), kind="dotplot")
        output_path = dirs["figures"] / f"03_dotplot_{lineage}.{figure_format}"
        try:
            render_scanpy_dotplot(
                adata_subset,
                var_names=markers_for_lineage,
                groupby=color_col,
                output_path=output_path,
                use_raw=source["use_raw"],
                figsize=(fig_w, fig_h),
                dpi=figure_dpi,
                standard_scale="var",
                dendrogram=True,
                cmap="Reds",
                vmin=-2,
                vmax=2,
            )
        except Exception:
            render_scanpy_dotplot(
                adata_subset,
                var_names=markers_for_lineage,
                groupby=color_col,
                output_path=output_path,
                use_raw=source["use_raw"],
                figsize=(fig_w, fig_h),
                dpi=figure_dpi,
                standard_scale="var",
                dendrogram=False,
                cmap="Reds",
                vmin=-2,
                vmax=2,
            )
        results[lineage] = output_path
    return results


def make_core_marker_heatmap(
    adata: Any,
    output_path: str | Path,
    marker_panels: Mapping[str, Sequence[str]] | None = None,
    color_col: str = "cell_type_L3",
    figure_dpi: int = 300,
) -> Path | None:
    marker_panels = dict(marker_panels or build_default_marker_panels())
    source = resolve_expression_source(adata)
    available_markers, _ = get_available_markers(marker_panels, source["gene_universe"])
    all_core_markers = deduplicate_preserve_order(
        gene for markers in available_markers.values() for gene in markers
    )
    if not all_core_markers:
        return None

    expression_df = build_mean_expression_matrix(
        adata,
        all_core_markers,
        groupby=color_col,
        use_raw=source["use_raw"],
    )
    expression_df.index.name = "L3 Cell Type"
    expression_df.columns.name = "Marker Gene"
    figsize = compute_dynamic_figsize(expression_df.shape[0], expression_df.shape[1], kind="heatmap")
    render_seaborn_heatmap(
        expression_df,
        output_path=output_path,
        title="Core Marker Expression Across L3 Cell Types",
        cbar_label="Mean Expression (scaled)",
        figsize=figsize,
        dpi=figure_dpi,
    )
    return Path(output_path)


def _build_review_panel_context(
    adata: Any,
    specs: Sequence[Mapping[str, Any]],
) -> tuple[pd.DataFrame, list[dict[str, Any]], list[dict[str, Any]]]:
    section_rows: list[dict[str, Any]] = []
    panel_specs: list[dict[str, Any]] = []
    section_layout: list[dict[str, Any]] = []
    current_x = 0

    for spec in specs:
        for side in ("wrong", "correct"):
            requested_genes = [str(gene) for gene in spec[f"{side}_genes"]]
            available_genes = [gene for gene in requested_genes if gene_available(adata, gene)]
            missing_genes = [gene for gene in requested_genes if gene not in available_genes]
            label = str(spec[f"{side}_label"])
            cluster = str(spec["cluster"])
            note = str(spec.get("note", ""))
            section_rows.append(
                {
                    "cluster": cluster,
                    "side": side,
                    "label": label,
                    "note": note,
                    "requested_genes": ";".join(requested_genes),
                    "available_genes": ";".join(available_genes),
                    "missing_genes": ";".join(missing_genes),
                }
            )
            if not available_genes:
                continue
            section_start = current_x
            for _gene in available_genes:
                current_x += 1
            section_layout.append(
                {
                    "cluster": cluster,
                    "side": side,
                    "label": label,
                    "start": section_start,
                    "end": current_x - 1,
                    "n_genes": len(available_genes),
                }
            )
            panel_specs.append(
                {
                    "cluster": cluster,
                    "side": side,
                    "label": label,
                    "note": note,
                    "genes": available_genes,
                }
            )
    return pd.DataFrame(section_rows), panel_specs, section_layout


def build_misannotation_review_table(
    adata: Any,
    specs: Sequence[Mapping[str, Any]] | None = None,
    cluster_col: str = "subcluster",
) -> pd.DataFrame:
    specs = list(specs or build_misannotation_review_specs())
    review_clusters = [str(spec["cluster"]) for spec in specs]
    _section_df, panel_specs, _layout = _build_review_panel_context(adata, specs)
    return build_review_dotplot_table(
        adata,
        review_clusters=review_clusters,
        panel_specs=panel_specs,
        groupby_col=cluster_col,
    )


def run_misannotation_review(
    adata: Any,
    output_dir: str | Path,
    specs: Sequence[Mapping[str, Any]] | None = None,
    cluster_col: str = "subcluster",
    figure_format: str = "pdf",
    figure_dpi: int = 300,
) -> dict[str, Any]:
    import matplotlib.pyplot as plt

    specs = list(specs or build_misannotation_review_specs())
    dirs = prepare_output_dirs(output_dir)
    cluster_series = normalize_string_series(safe_obs_column(adata, cluster_col))
    review_clusters = [str(spec["cluster"]) for spec in specs]
    review_counts = cluster_series.value_counts(dropna=False).reindex(review_clusters, fill_value=0)
    section_df, panel_specs, section_layout = _build_review_panel_context(adata, specs)
    plot_df = build_review_dotplot_table(
        adata,
        review_clusters=review_clusters,
        panel_specs=panel_specs,
        groupby_col=cluster_col,
    )
    if plot_df.empty:
        raise ValueError("Misannotation review dotplot has no data to plot")

    counts_path = write_tsv(review_counts.rename("n_cells").to_frame(), dirs["tables"] / "misannotation_review_cluster_counts.tsv")
    plan_path = write_tsv(section_df, dirs["tables"] / "misannotation_review_plan.tsv", index=False)
    data_path = write_tsv(plot_df, dirs["tables"] / "misannotation_review_dotplot_data.tsv", index=False)

    plot_label_order = deduplicate_preserve_order(plot_df["plot_label"].tolist())
    fig_width, fig_height = compute_dynamic_figsize(len(review_clusters), len(plot_label_order), kind="review")
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=np.clip(plot_df["pct_expr"] * 500, 10, 500),
        c=plot_df["mean_expr_scaled"],
        cmap="Reds",
        vmin=0,
        vmax=1,
        edgecolors="black",
        linewidths=0.15,
    )
    ax.set_xticks(range(len(plot_label_order)))
    ax.set_xticklabels([label.split("|")[-1] for label in plot_label_order], rotation=90, fontsize=8)
    ax.set_yticks(range(len(review_clusters)))
    ax.set_yticklabels(review_clusters, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Marker genes (wrong vs correct panels)", fontsize=12)
    ax.set_ylabel("Suspected misannotation subclusters", fontsize=12)
    ax.set_title("Misannotation review dotplot: wrong vs correct marker programs", fontsize=14, pad=18)

    for section in section_layout:
        if section["start"] > 0:
            ax.axvline(section["start"] - 0.5, color="lightgray", linewidth=0.8)
        center = (section["start"] + section["end"]) / 2
        ax.text(
            center,
            1.02,
            f"{section['cluster']}\n{section['side']} → {section['label']}",
            rotation=60,
            ha="left",
            va="bottom",
            fontsize=8,
            transform=ax.get_xaxis_transform(),
        )

    cbar = plt.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("Scaled mean log1p expression", rotation=90)
    size_breaks = [0.1, 0.25, 0.5, 0.75, 1.0]
    legend_handles = [
        ax.scatter([], [], s=np.clip(value * 500, 10, 500), c="lightgray", edgecolors="black", linewidths=0.15)
        for value in size_breaks
    ]
    ax.legend(
        legend_handles,
        [f"{int(value * 100)}%" for value in size_breaks],
        title="% expressing",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
    )
    plt.subplots_adjust(top=0.72, bottom=0.28, right=0.88)
    plot_path = save_current_figure(
        dirs["figures"] / f"03b_dotplot_misannotation_review.{figure_format}",
        dpi=figure_dpi,
        close=True,
        figure=fig,
    )

    return {
        "counts_path": counts_path,
        "plan_path": plan_path,
        "data_path": data_path,
        "plot_path": plot_path,
        "plot_df": plot_df,
        "section_df": section_df,
        "review_counts": review_counts,
    }


if __name__ == "__main__":
    print("[INFO] epithelial_analysis_helper_20260423_v1.py loaded successfully")
