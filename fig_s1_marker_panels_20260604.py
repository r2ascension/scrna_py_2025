#!/usr/bin/env python3
"""
Fig S1: Comprehensive marker gene distribution panels.
Creates pan-lineage + per-lineage dotplots, heatmaps, and UMAP overviews.

Output: /home/h2048/data/py/0604/fig_s1_marker_panels/
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import warnings
import os
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# ── Output ──────────────────────────────────────────────────
OUT = Path("/home/h2048/data/py/0604/fig_s1_marker_panels")
OUT.mkdir(parents=True, exist_ok=True)

# ── Colours ─────────────────────────────────────────────────
L1_COLORS = {
    "Epithelial": "#1f77b4",
    "Endothelial": "#ff7f0e",
    "Fibroblast": "#2ca02c",
    "SMC": "#d62728",
    "Myeloid": "#9467bd",
    "T": "#8c564b",
    "B": "#e377c2",
}

TISSUE_COLORS = {
    "nose": "#fdae61",
    "sinus": "#abd9e9",
    "respiratory airway": "#2c7bb6",
    "lung parenchyma": "#d7191c",
}

# ── Canonical marker gene sets ──────────────────────────────
# Pan-L1 canonical markers
L1_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "SMC": ["ACTA2", "MYH11", "TAGLN", "CNN1"],
    "Myeloid": ["CD68", "CD14", "AIF1", "FCGR3A"],
    "T": ["CD3D", "CD3E", "CD2", "CD7"],
    "B": ["CD19", "CD79A", "MS4A1", "PAX5"],
}

# L2 canonical markers (broad)
L2_CANONICAL = {
    # B cells
    "Naive_B": ["MS4A1", "CD19", "IGHD", "CD79A"],
    "Memory_B": ["MS4A1", "CD27", "CD79A"],
    "GC_B": ["MS4A1", "BCL6", "AICDA", "MKI67"],
    "Plasma": ["SDC1", "MZB1", "XBP1", "JCHAIN"],
    "Atypical_Memory_B": ["MS4A1", "FCRL4", "CD27"],
    # T/NK
    "CD4 T cells": ["CD4", "CD3D", "IL7R"],
    "CD8 T cells": ["CD8A", "CD8B", "CD3D", "GZMK"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "PRF1"],
    # Myeloid
    "Alveolar macrophages": ["MARCO", "FABP4", "CD68", "PPARG"],
    "Interstitial macrophages": ["CD68", "CD163", "MRC1", "MAF"],
    "Macrophages": ["CD68", "CD163", "C1QA", "C1QB"],
    "Classical monocytes": ["S100A8", "S100A9", "CD14", "VCAN"],
    "Non-classical monocytes": ["FCGR3A", "CDKN1C", "LILRB2"],
    "DC2": ["CD1C", "FCER1A", "CLEC10A"],
    "DC": ["CLEC9A", "XCR1", "BATF3"],
    "Mast cells": ["TPSAB1", "CPA3", "KIT", "HPGDS"],
    "pDC": ["CLEC4C", "LILRA4", "IRF7", "TCF4"],
    # Epithelial
    "Alveolar": ["SFTPC", "SFTPB", "SFTPD", "AGER", "AQP4"],
    "Basal_Lineage": ["KRT5", "KRT14", "TP63", "KRT15"],
    "Ciliated_Lineage": ["FOXJ1", "SNTN", "TPPP3", "DNAH5"],
    "Secretory_Lineage": ["MUC5AC", "MUC5B", "SCGB1A1", "SCGB3A2"],
    "Rare_Specialized": ["FOXI1", "CFTR", "ASCL3", "POU2F3"],
    "SMG": ["LTF", "LYZ", "PIP", "STATH"],
    # Endothelial
    "Endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5"],
    # Fibroblast
    "Fibroblast": ["COL1A1", "DCN", "LUM", "PDGFRA"],
    # SMC
    "Smooth_Muscle": ["ACTA2", "MYH11", "TAGLN"],
    "Pericyte": ["RGS5", "CSPG4", "PDGFRB"],
}

# ── Helper functions ────────────────────────────────────────
def safe_load(path, backed=False):
    """Load h5ad, handling large files."""
    try:
        if backed:
            return sc.read_h5ad(path, backed="r")
        return sc.read_h5ad(path)
    except Exception as e:
        print(f"  WARNING: {e}")
        return None


def get_top_markers(adata, groupby, n_top=5, min_fc=1.5):
    """Compute top marker genes per group using Wilcoxon rank-sum."""
    sc.tl.rank_genes_groups(
        adata, groupby, method="wilcoxon", n_genes=n_top * 2,
        pts=True, tie_correct=True
    )
    result = adata.uns["rank_genes_groups"]
    markers = {}
    for group in result["names"].dtype.names:
        genes = result["names"][group][:n_top * 2]
        scores = result["scores"][group][:n_top * 2]
        pvals = result["pvals"][group][:n_top * 2]
        logfcs = result.get("logfoldchanges", None)
        if logfcs is not None:
            logfcs = logfcs[group][:n_top * 2]
        filtered = []
        for i, g in enumerate(genes):
            if logfcs is not None and logfcs[i] < np.log2(min_fc):
                continue
            if pvals[i] > 0.01:
                continue
            filtered.append(g)
            if len(filtered) >= n_top:
                break
        markers[group] = filtered
    return markers


def save_figure(fig, name, dpi=300):
    """Save figure as both PDF and PNG."""
    for fmt, d in [("pdf", 150), ("png", dpi)]:
        path = OUT / f"{name}.{fmt}"
        fig.savefig(path, dpi=d, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {path}")


# ═══════════════════════════════════════════════════════════
# Panel A: UMAP overview — L1 cell types
# ═══════════════════════════════════════════════════════════
print("=" * 60)
print("PANEL A: UMAP L1 cell types (total h5ad)")
print("=" * 60)

adata_total = sc.read_h5ad(
    "/home/h2048/data/coredata0524/coredata0524_allcells_scanvi_l1_20260524.h5ad"
)
print(f"  Loaded: {adata_total.shape}")

# Subsample for plotting speed (use all cells but rasterize)
fig, ax = plt.subplots(figsize=(10, 8))
for ct in sorted(adata_total.obs["cell_type"].unique()):
    mask = adata_total.obs["cell_type"] == ct
    # plot background first
    bg = ~mask
    ax.scatter(
        adata_total.obsm["X_umap"][bg, 0][::5],
        adata_total.obsm["X_umap"][bg, 1][::5],
        s=0.3, c="lightgray", rasterized=True, alpha=0.3,
    )
    ax.scatter(
        adata_total.obsm["X_umap"][mask, 0],
        adata_total.obsm["X_umap"][mask, 1],
        s=1, c=L1_COLORS.get(ct, "gray"),
        label=ct, rasterized=True, alpha=0.7,
    )
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title("Cell Type Overview (L1)")
ax.legend(
    markerscale=8, loc="lower left",
    bbox_to_anchor=(1.01, 0), frameon=True,
    title="Cell Type",
)
fig.tight_layout()
save_figure(fig, "figS1_panelA_umap_L1_overview")
plt.close()

del adata_total

# ═══════════════════════════════════════════════════════════
# Panel B: UMAP split by tissue
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("PANEL B: UMAP split by tissue")
print("=" * 60)

adata_total = sc.read_h5ad(
    "/home/h2048/data/coredata0524/coredata0524_allcells_scanvi_l1_20260524.h5ad"
)

tissues = sorted(adata_total.obs["tissue"].dropna().unique())
n_tissues = len(tissues)
fig, axes = plt.subplots(1, n_tissues, figsize=(5 * n_tissues, 5))

for ax, tissue in zip(axes, tissues):
    mask = adata_total.obs["tissue"] == tissue
    # background
    bg_mask = ~mask
    ax.scatter(
        adata_total.obsm["X_umap"][:, 0][bg_mask][::10],
        adata_total.obsm["X_umap"][:, 1][bg_mask][::10],
        s=0.3, c="lightgray", rasterized=True, alpha=0.2,
    )
    ax.scatter(
        adata_total.obsm["X_umap"][:, 0][mask],
        adata_total.obsm["X_umap"][:, 1][mask],
        s=1, c=TISSUE_COLORS.get(tissue, "gray"),
        rasterized=True, alpha=0.7,
    )
    ax.set_title(f"{tissue}\n({mask.sum():,} cells)")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

fig.suptitle("Tissue Distribution", fontsize=14, y=1.02)
fig.tight_layout()
save_figure(fig, "figS1_panelB_umap_tissue_split")
plt.close()

del adata_total

# ═══════════════════════════════════════════════════════════
# Panel C: Pan-L1 canonical marker dotplot
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("PANEL C: Pan-L1 canonical marker dotplot")
print("=" * 60)

adata_total = sc.read_h5ad(
    "/home/h2048/data/coredata0524/coredata0524_allcells_scanvi_l1_20260524.h5ad"
)

# Build gene list: flatten L1_MARKERS
all_markers = []
gene_to_ct = {}
for ct, genes in L1_MARKERS.items():
    for g in genes:
        if g not in gene_to_ct:
            all_markers.append(g)
            gene_to_ct[g] = ct

# Filter to genes present in adata
available_markers = [g for g in all_markers if g in adata_total.var_names]
print(f"  {len(available_markers)}/{len(all_markers)} markers available")

# Create dotplot
sc.pl.dotplot(
    adata_total, available_markers, groupby="cell_type",
    dendrogram=False, standard_scale="var",
    title="Canonical Marker Genes Across Cell Types (L1)",
    show=False,
)
fig = plt.gcf()
fig.set_size_inches(18, 5)
fig.tight_layout()
save_figure(fig, "figS1_panelC_dotplot_L1_canonical")
plt.close()

# Also do L2 canonical dotplot with subsampling for performance
print("\n  Creating L2 canonical dotplot...")
# Subsample to 30k cells for performance
sc.settings.n_jobs = 4
adata_L2 = adata_total[
    adata_total.obs["cell_type_L2"].isin(
        adata_total.obs["cell_type_L2"].value_counts().head(25).index
    )
].copy()

l2_flat = []
for ct, genes in L2_CANONICAL.items():
    for g in genes:
        if g not in l2_flat:
            l2_flat.append(g)
l2_available = [g for g in l2_flat if g in adata_L2.var_names]

sc.pl.dotplot(
    adata_L2, l2_available, groupby="cell_type_L2",
    dendrogram=True, standard_scale="var",
    title="Canonical Marker Genes Across Cell Types (L2)",
    show=False,
)
fig = plt.gcf()
fig.set_size_inches(22, 8)
fig.tight_layout()
save_figure(fig, "figS1_panelC2_dotplot_L2_canonical")
plt.close()

del adata_total, adata_L2

# ═══════════════════════════════════════════════════════════
# Panel D: Heatmap of top markers per L2 cell type
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("PANEL D: Heatmap of top differentially expressed genes per L2")
print("=" * 60)

adata_total = sc.read_h5ad(
    "/home/h2048/data/coredata0524/coredata0524_allcells_scanvi_l1_20260524.h5ad"
)

# Subsample to max 500 cells per L2 type for manageable heatmap
target_per_group = 300
sc.pp.subsample(adata_total, groupby="cell_type_L2", n_obs=target_per_group, copy=False)
print(f"  Subsampled to {adata_total.shape[0]} cells for heatmap")

# Compute top markers
print("  Computing differential expression...")
sc.tl.rank_genes_groups(
    adata_total, groupby="cell_type_L2", method="wilcoxon",
    n_genes=3, pts=True, tie_correct=True,
)
sc.pl.rank_genes_groups_heatmap(
    adata_total, n_genes=5, groupby="cell_type_L2",
    standard_scale="var", cmap="RdBu_r",
    show_gene_labels=True, show=False,
)
fig = plt.gcf()
fig.set_size_inches(20, 14)
save_figure(fig, "figS1_panelD_heatmap_L2_top_markers")
plt.close()

del adata_total

# ═══════════════════════════════════════════════════════════
# Panels E-K: Per-lineage refined dotplots (full-gene h5ads)
# ═══════════════════════════════════════════════════════════

LINEAGE_CONFIGS = [
    {
        "name": "epithelial",
        "path": "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/epithelial_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L2",
        "umap_key": "X_umap",
        "title": "Epithelial Lineage",
        "panel": "E",
    },
    {
        "name": "bcell",
        "path": "/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L2",
        "umap_key": "X_umap",
        "title": "B Cell Lineage",
        "panel": "F",
    },
    {
        "name": "tnk",
        "path": "/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508/tnk_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L2",
        "umap_key": "X_umap",
        "title": "T/NK Cell Lineage",
        "panel": "G",
    },
    {
        "name": "myeloid",
        "path": "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L2",
        "umap_key": "X_umap",
        "title": "Myeloid Lineage",
        "panel": "H",
    },
    {
        "name": "endothelial",
        "path": "/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L3",
        "umap_key": "X_umap",
        "title": "Endothelial Lineage",
        "panel": "I",
    },
    {
        "name": "fibroblast",
        "path": "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L3",
        "umap_key": "X_umap",
        "title": "Fibroblast Lineage",
        "panel": "J",
    },
    {
        "name": "smc",
        "path": "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final_fullgene.h5ad",
        "groupby": "cell_type_L3",
        "umap_key": "X_umap",
        "title": "SMC / Pericyte / Neuron-like Lineage",
        "panel": "K",
    },
]

for cfg in LINEAGE_CONFIGS:
    print(f"\n{'=' * 60}")
    print(f"PANEL {cfg['panel']}: {cfg['title']}")
    print(f"{'=' * 60}")

    if not os.path.exists(cfg["path"]):
        print(f"  SKIP: file not found: {cfg['path']}")
        continue

    try:
        adata = sc.read_h5ad(cfg["path"])
        print(f"  Loaded: {adata.shape}")
    except Exception as e:
        print(f"  SKIP: {e}")
        continue

    # Find best groupby column
    gb = cfg["groupby"]
    if gb not in adata.obs.columns:
        # Try alternatives
        for alt in ["cell_type_L2", "cell_type_L3", "CHOIR_clusters_0.2", "leiden"]:
            if alt in adata.obs.columns:
                gb = alt
                break

    n_groups = adata.obs[gb].nunique()
    print(f"  Groupby: {gb} ({n_groups} groups)")
    print(f"  Groups: {sorted(adata.obs[gb].unique().tolist())}")

    # --- UMAP panel ---
    if cfg["umap_key"] in adata.obsm:
        fig, ax = plt.subplots(figsize=(8, 7))
        groups = sorted(adata.obs[gb].unique())
        colors = plt.cm.tab20(np.linspace(0, 1, len(groups)))
        for i, group in enumerate(groups):
            mask = adata.obs[gb] == group
            ax.scatter(
                adata.obsm[cfg["umap_key"]][:, 0][mask],
                adata.obsm[cfg["umap_key"]][:, 1][mask],
                s=2, c=[colors[i]], label=group,
                rasterized=True, alpha=0.8,
            )
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"{cfg['title']} — {gb}")
        ax.legend(markerscale=5, bbox_to_anchor=(1.01, 1),
                  loc="upper left", fontsize=7, frameon=True)
        fig.tight_layout()
        save_figure(fig, f"figS1_panel{cfg['panel']}_umap_{cfg['name']}")
        plt.close()

    # --- Dotplot ---
    n_top = min(4, max(2, 20 // max(n_groups, 1)))
    try:
        sc.pl.dotplot(
            adata, var_names=None, groupby=gb,
            n_genes=n_top, dendrogram=(n_groups > 3),
            standard_scale="var",
            title=f"{cfg['title']} — Top {n_top} Markers per Cluster",
            show=False,
        )
        fig = plt.gcf()
        fig.set_size_inches(max(14, n_groups * 1.2), 6)
        fig.tight_layout()
        save_figure(fig, f"figS1_panel{cfg['panel']}_dotplot_{cfg['name']}")
        plt.close()
    except Exception as e:
        print(f"  Dotplot failed: {e}")

    # --- Heatmap ---
    if adata.shape[0] > 10000:
        sc.pp.subsample(adata, n_obs=min(adata.shape[0], 200 * n_groups), copy=False)

    try:
        sc.tl.rank_genes_groups(
            adata, groupby=gb, method="wilcoxon",
            n_genes=n_top, pts=True, tie_correct=True,
        )
        sc.pl.rank_genes_groups_heatmap(
            adata, n_genes=n_top, groupby=gb,
            standard_scale="var", cmap="RdBu_r",
            show_gene_labels=True, show=False,
        )
        fig = plt.gcf()
        fig.set_size_inches(max(12, n_groups * 0.8), max(8, n_groups * 0.5))
        fig.tight_layout()
        save_figure(fig, f"figS1_panel{cfg['panel']}_heatmap_{cfg['name']}")
        plt.close()
    except Exception as e:
        print(f"  Heatmap failed: {e}")

    # --- UMAP split by tissue (if tissue column exists) ---
    tissue_col = None
    for col in ["tissue", "Tissue", "condition"]:
        if col in adata.obs.columns:
            tissue_col = col
            break

    if tissue_col and cfg["umap_key"] in adata.obsm:
        tissue_vals = sorted(adata.obs[tissue_col].dropna().unique())
        n_tis = len(tissue_vals)
        if 2 <= n_tis <= 6:
            fig, axes = plt.subplots(
                1, n_tis, figsize=(5 * n_tis, 5),
                squeeze=False,
            )
            for ax, tv in zip(axes[0], tissue_vals):
                mask = adata.obs[tissue_col] == tv
                bg = ~mask
                ax.scatter(
                    adata.obsm[cfg["umap_key"]][:, 0][bg][::5],
                    adata.obsm[cfg["umap_key"]][:, 1][bg][::5],
                    s=0.3, c="lightgray", rasterized=True, alpha=0.2,
                )
                ax.scatter(
                    adata.obsm[cfg["umap_key"]][:, 0][mask],
                    adata.obsm[cfg["umap_key"]][:, 1][mask],
                    s=2, rasterized=True, alpha=0.7,
                )
                ax.set_title(f"{tv}\n({mask.sum():,} cells)")
                ax.set_xlabel("UMAP 1")
                ax.set_ylabel("UMAP 2")
            fig.suptitle(f"{cfg['title']} — Split by {tissue_col}",
                        fontsize=12, y=1.02)
            fig.tight_layout()
            save_figure(
                fig,
                f"figS1_panel{cfg['panel']}_umap_{cfg['name']}_split_{tissue_col}",
            )
            plt.close()

    del adata

# ═══════════════════════════════════════════════════════════
# Summary manifest
# ═══════════════════════════════════════════════════════════
print(f"\n{'=' * 60}")
print("FIG S1 GENERATION COMPLETE")
print(f"{'=' * 60}")
print(f"\nOutput directory: {OUT}")
print("\nFiles generated:")
for f in sorted(OUT.iterdir()):
    size_mb = f.stat().st_size / 1e6
    print(f"  {f.name} ({size_mb:.1f} MB)")
