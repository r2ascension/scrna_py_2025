#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stromal/Vascular Cell Marker Visualization for Cluster Confirmation
====================================================================
Author: r2end
Date: 2026-02-14
Version: v1.0

Purpose:
    Comprehensive marker-based visualization to confirm subcluster identities
    in the merged stromal/vascular dataset. Covers all 31 cell type categories
    from the curated marker panel, organized into hierarchical dotplots and
    UMAP feature plots.

Input:
    adata_stromal_subclustered_FINAL_v2_20260128.h5ad

Output:
    figures/
        01_umap_L2_overview.pdf
        02_umap_L3_subclusters.pdf
        02b_umap_L2_L3_combined.pdf
        03_umap_batch.pdf
        04_dotplot_endothelium_{scaled,raw}.pdf  -- 11 endothelial subtypes
        05_dotplot_fibroblast_{scaled,raw}.pdf    -- 8 fibroblast subtypes
        06_dotplot_pericyte_smc_{scaled,raw}.pdf  -- 7 pericyte/SMC/Schwann subtypes
        07_dotplot_contaminants_{scaled,raw}.pdf  -- 5 contaminant types (B/Plasma split)
        08_dotplot_global_overview_{scaled,raw}.pdf
        09_feature_umap_endothelium.pdf
        10_feature_umap_fibroblast.pdf
        11_feature_umap_pericyte_smc.pdf
        12_feature_umap_discriminatory.pdf        -- colorbars enabled

Key Markers Source: Curated discriminatory panel (r2end, 2026-02-14)
"""

# ============================================================================
# BACKEND SETUP (must be before pyplot import)
# ============================================================================
import os
# Prevent BLAS/numexpr thread over-subscription; visualization tasks are
# single-threaded by nature and 48 threads causes resource contention.
for _k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[_k] = "1"

import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import gc
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, dpi_save=300, frameon=False,
                               facecolor='white', fontsize=11)
sc.settings.n_jobs = 1    # visualization tasks; CPU-bound parallelism unhelpful
np.random.seed(42)

# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_H5AD = Path("/home/h2048/data/py/0120/stromal_analysis_unified/results/"
                  "subcluster_unified_v2_20260128/"
                  "adata_stromal_subclustered_FINAL_v2_20260128.h5ad")

OUTPUT_DIR = Path("/home/h2048/data/py/0214/stromal_marker_visualization")
FIG_DIR    = OUTPUT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Annotation keys
CELLTYPE_L2   = 'cell_type_L2'    # Major cell type (Fibroblast, Endothelial, etc.)
CELLTYPE_L3   = 'cell_type_L3'    # Subcluster label  (Fibroblast_c0, etc.)
BATCH_KEY     = 'sample'

# Visualization settings
DPI         = 300
FIG_FORMAT  = 'pdf'
UMAP_SIZE   = 3        # pt, adjust if needed

# ============================================================================
# MARKER PANEL DEFINITION
# ============================================================================
# Each entry: (group_label, [gene list])
# Genes selected for maximum discriminatory power; 3-5 per type

ENDOTHELIAL_MARKERS = {
    # ---- Lymphatic ----
    'Quiescent Lymphatic\n(PROX1+LYVE1+)':
        ['PROX1', 'FLT4', 'PDPN', 'LYVE1', 'CCL21', 'FOXC2'],

    'Venous-like Lymphatic\n(FOXC2+MME+NR2F2+)':
        ['FOXC2', 'MME', 'FLT4', 'PDPN', 'LYVE1', 'NR2F2'],

    # PROX1 included as explicit negative contrast (should be low/absent here)
    'LYVE1+PROX1- Immune\nInteracting':
        ['LYVE1', 'CCL14', 'KDR', 'ICAM1', 'SELE', 'PROX1'],

    # ---- Capillary ----
    'Aerocyte-like Capillary\n(aCap: CA4+AGER+)':
        ['CA4', 'AGER', 'EDNRB', 'ICAM2', 'KDR', 'RAMP2'],

    'General Capillary\n(gCap: GPIHBP1+PLVAP+)':
        ['GPIHBP1', 'RGCC', 'PLVAP', 'EMCN', 'TMEM100'],

    'IFN-primed Capillary\n(ISG15+MX1+)':
        ['ISG15', 'IFIT1', 'MX1', 'OAS1', 'STAT1', 'B2M'],

    # ---- Arterial ----
    'Quiescent Arterial\n(GJA5+HEY1+SEMA3G+)':
        ['GJA5', 'HEY1', 'EFNB2', 'SOX17', 'DLL4', 'SEMA3G'],

    # GJA5/HEY1 included to anchor EC identity; avoids mis-reading as APC/immune
    'Immunomodulatory Arterial\n(IDO1+IRF1+GJA5+)':
        ['IDO1', 'IRF1', 'STAT1', 'B2M', 'IFITM1', 'GJA5', 'HEY1'],

    # ---- Venous ----
    # VWF removed: not specific to venous (broad EC); PLVAP+NR2F2 more discriminatory
    'Homeostatic Venous\n(ACKR1+NR2F2+PLVAP+)':
        ['ACKR1', 'NR2F2', 'PLVAP', 'SEPP1', 'EMCN'],

    'Immune-recruiting Venous\n(ACKR1+SELE+VCAM1+)':
        ['ACKR1', 'SELE', 'VCAM1', 'ICAM1', 'CCL2'],

    # ---- Angiogenic ----
    'Angiogenic Tip\n(CXCR4+APLN+ESM1+)':
        ['CXCR4', 'APLN', 'ANGPT2', 'ESM1', 'PGF', 'KDR'],
}

FIBROBLAST_MARKERS = {
    # ---- Adventitial ----
    'Perivascular Adventitial\n(PI16+MFAP5+CD34+)':
        ['PI16', 'MFAP5', 'CD34', 'COL15A1', 'DPT', 'C7'],

    'Niche-supporting Adventitial\n(PI16+SFRP4+CXCL14+)':
        ['PI16', 'SFRP4', 'MFAP5', 'C7', 'CXCL14'],

    # ---- Alveolar/Tissue ----
    # DCN/LUM retained as one anchor but COL13A1/CFD are more alveolar-specific;
    # dropping LUM prevents over-overlap with adventitial fibro
    'Homeostatic Alveolar Fibro\n(PDGFRA+TCF21+COL13A1+)':
        ['PDGFRA', 'TCF21', 'COL13A1', 'CFD', 'DCN'],

    'Activated ECM-high Fibro\n(POSTN+CTHRC1+FN1+)':
        ['FN1', 'POSTN', 'CTHRC1', 'COL3A1', 'COL1A1'],

    # DCN as fibro-lineage anchor; prevents mis-assignment of FOS/JUN+ EC or SMC
    'Stress-activated Fibro\n(FOS+JUN+MRTFB+)':
        ['MRTFB', 'NIBAN2', 'FOS', 'JUN', 'HSPA1A', 'DCN'],

    # ---- Specialized ----
    # PDGFRB (pericyte) + MYH11 (vSMC) as negative contrasts; distinguishes
    # myofibro contractile program from true SMC/pericyte contractility
    'Myofibroblasts\n(ACTA2+TAGLN+POSTN+)':
        ['ACTA2', 'TAGLN', 'MYL9', 'CNN1', 'TPM2', 'POSTN', 'PDGFRB', 'MYH11'],

    'Lipofibroblasts\n(APOE+FABP4+PPARG+)':
        ['APOE', 'APOD', 'FABP4', 'PLIN2', 'PPARG', 'IGFBP5'],

    'Developmental-signaling Fibro\n(FGF10+TWIST1+WNT2+)':
        ['FGF10', 'TWIST1', 'WNT2', 'BMP4', 'TBX4'],
}

PERICYTE_SMC_MARKERS = {
    # ---- Pericyte ----
    'Pericyte Type 1\n(RGS5+ABCC9+KCNJ8+)':
        ['RGS5', 'PDGFRB', 'ABCC9', 'KCNJ8', 'NOTCH3', 'MCAM'],

    'Pericyte Type 2\n(CSPG4+DES+RERGL+)':
        ['CSPG4', 'PDGFRB', 'MCAM', 'DES', 'RERGL'],

    'Contractile Pericyte\n(ACTA2+RGS5+ABCC9+)':
        ['ACTA2', 'TAGLN', 'MYL9', 'RGS5', 'PDGFRB', 'ABCC9'],

    # CD34 removed: also high in PI16+ adventitial fibro; RERGL/NOTCH3 more specific
    'Precursor-like Pericyte\n(CSPG4+RERGL+NOTCH3+)':
        ['CSPG4', 'PDGFRB', 'MCAM', 'RERGL', 'NOTCH3'],

    # ---- Smooth Muscle ----
    'Airway Smooth Muscle\n(PRKG1+CASQ2+RGS4+)':
        ['PRKG1', 'CASQ2', 'ACTA2', 'MYLK', 'RGS4'],

    'Vascular Smooth Muscle\n(MYH11+SMTN+CNN1+)':
        ['MYH11', 'SMTN', 'CNN1', 'TAGLN', 'ACTA2'],

    # ---- Neural ----
    'Schwann / Neural-like\n(SOX10+S100B+ERBB3+)':
        ['SOX10', 'S100B', 'ERBB3', 'GFRA3', 'PLP1', 'MPZ'],
}

CONTAMINANT_MARKERS = {
    'IFN-activated Immune\n(PTPRC+ISG15+LST1+)':
        ['PTPRC', 'LST1', 'TYROBP', 'ISG15', 'IFIT1', 'FCER1A'],

    'MHC-II-high APC-like\n(HLA-DRA+CIITA+CD74+)':
        ['HLA-DRA', 'CD74', 'CIITA', 'HLA-DRB1', 'CLEC10A'],

    # Split B vs Plasma: merged panel conflates two distinct populations
    # whose dotplot peaks are often on different subclusters
    'B-cell contamination\n(MS4A1+CD79A+CD37+)':
        ['MS4A1', 'CD79A', 'CD74', 'HLA-DRA', 'CD37'],

    'Plasma-cell contamination\n(MZB1+JCHAIN+XBP1+)':
        ['MZB1', 'JCHAIN', 'XBP1', 'SDC1', 'IGKC'],

    'Osteogenic\n(ALPL+SP7+RUNX2+)':
        ['ALPL', 'SP7', 'BGLAP', 'IBSP', 'RUNX2'],
}

# Global minimal panel: 3 genes per category (for the combined overview)
GLOBAL_MINIMAL = {
    # Endothelium
    'Lymphatic EC':         ['PROX1', 'FLT4', 'LYVE1'],
    'Capillary EC (gCap)':  ['GPIHBP1', 'PLVAP', 'RGCC'],
    'Capillary EC (aCap)':  ['CA4', 'AGER', 'EDNRB'],
    'IFN-EC':               ['ISG15', 'IFIT1', 'MX1'],
    'Arterial EC':          ['GJA5', 'HEY1', 'EFNB2'],
    'Venous EC':            ['ACKR1', 'NR2F2', 'PLVAP'],
    'Angiogenic tip EC':    ['CXCR4', 'APLN', 'ESM1'],
    # Fibroblast
    'Adventitial Fibro':    ['PI16', 'MFAP5', 'COL15A1'],
    'Alveolar Fibro':       ['PDGFRA', 'TCF21', 'COL13A1'],
    'Activated Fibro':      ['POSTN', 'CTHRC1', 'FN1'],
    'Lipofibro':            ['APOE', 'FABP4', 'PPARG'],
    'Myofibro':             ['ACTA2', 'TAGLN', 'CNN1'],
    # Pericyte / SMC
    'Pericyte':             ['RGS5', 'PDGFRB', 'ABCC9'],
    'Airway SMC':           ['PRKG1', 'CASQ2', 'MYLK'],
    'Vascular SMC':         ['MYH11', 'SMTN', 'NOTCH3'],
    'Schwann':              ['SOX10', 'S100B', 'ERBB3'],
    # Contamination markers
    'Immune contamination': ['PTPRC', 'LST1', 'TYROBP'],
    'B cells':              ['MS4A1', 'CD79A', 'CD37'],
    'Plasma cells':         ['MZB1', 'JCHAIN', 'IGKC'],
}

# Key discriminatory genes for feature UMAP (organized by lineage)
FEATURE_UMAP_ENDO = [
    'PECAM1',   # Pan-endothelial
    'PROX1',    # Lymphatic
    'LYVE1',    # Lymphatic / LYVE1+
    'CA4',      # aCap
    'GPIHBP1',  # gCap
    'GJA5',     # Arterial
    'ACKR1',    # Venous
    'ISG15',    # IFN-EC
    'CXCR4',    # Angiogenic tip
    'SELE',     # Immune-recruiting venous
    'KDR',      # Pan-VEGFR / tip
    'IDO1',     # Immunomodulatory arterial
]

FEATURE_UMAP_FIBRO = [
    'DCN',      # Pan-fibroblast
    'PI16',     # Adventitial
    'MFAP5',    # Adventitial
    'CD34',     # Perivascular / adventitial
    'PDGFRA',   # Alveolar
    'TCF21',    # Alveolar
    'POSTN',    # Activated
    'CTHRC1',   # Activated / ECM
    'FN1',      # ECM-high
    'FABP4',    # Lipofibro
    'APOE',     # Lipofibro
    'FOS',      # Stress-activated
]

FEATURE_UMAP_PERI_SMC = [
    'RGS5',     # Pericyte type 1
    'PDGFRB',   # Pan-pericyte
    'CSPG4',    # Pericyte type 2
    'ABCC9',    # Pericyte
    'ACTA2',    # Contractile / myofibro / SMC
    'MYH11',    # Vascular SMC
    'PRKG1',    # Airway SMC
    'CASQ2',    # Airway SMC
    'SOX10',    # Schwann
    'S100B',    # Schwann
    'ERBB3',    # Schwann
    'MKI67',    # Proliferation
]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def filter_available(gene_list, var_names):
    """Return genes that exist in the data (use .raw if available)."""
    available = [g for g in gene_list if g in var_names]
    missing   = [g for g in gene_list if g not in var_names]
    if missing:
        print(f"  [warn] not in data: {missing}")
    return available


def make_dotplot(adata, marker_dict, groupby, title, figsize,
                 use_raw=True, vmin=None, vmax=3, dendrogram=False,
                 standard_scale='var', save_path=None):
    """
    Create a grouped dotplot.

    marker_dict: OrderedDict-like  {group_label: [gene1, gene2, ...]}
    standard_scale: 'var' (normalize per gene, shows relative specificity) or
                    None (raw expression scale, shows absolute strength).
    When save_path is given, outputs both a *_scaled.pdf and a *_raw.pdf so
    that relative specificity and absolute expression strength can both be read.
    """
    has_raw = adata.raw is not None
    var_names_ref = adata.raw.var_names if (use_raw and has_raw) else adata.var_names

    # Build filtered var_names dict (drop empty groups silently)
    var_names_filtered = {
        k: filter_available(v, var_names_ref)
        for k, v in marker_dict.items()
    }
    var_names_filtered = {k: v for k, v in var_names_filtered.items() if len(v) > 0}

    if not var_names_filtered:
        print("[error] No genes available for dotplot!")
        return

    n_genes = sum(len(v) for v in var_names_filtered.values())
    print(f"  Dotplot: {n_genes} genes x {adata.obs[groupby].nunique()} groups")

    def _draw(ss, path):
        try:
            fig = sc.pl.dotplot(
                adata,
                var_names=var_names_filtered,
                groupby=groupby,
                use_raw=use_raw and has_raw,
                vmin=vmin,
                vmax=vmax,
                standard_scale=ss,
                colorbar_title='Mean\nexpression' if ss is None else 'Scaled\nexpression',
                figsize=figsize,
                title=title,
                dendrogram=dendrogram,
                show=False,
                return_fig=True
            )
            if path:
                fig.savefig(path, dpi=DPI, bbox_inches='tight')
                print(f"  Saved: {path.name}")
        except Exception as e:
            print(f"  [error] Dotplot failed: {e}")
        finally:
            plt.close('all')

    if save_path:
        p = Path(save_path)
        stem, suffix = p.stem, p.suffix
        # scaled version (relative specificity)
        _draw('var', p.parent / f'{stem}_scaled{suffix}')
        # raw version (absolute expression strength)
        _draw(None,  p.parent / f'{stem}_raw{suffix}')
    else:
        _draw(standard_scale, None)


def make_feature_umap(adata, genes, title, ncols=4, use_raw=True,
                       basis='X_umap', size=UMAP_SIZE, show_colorbar=False,
                       save_path=None):
    """Feature UMAPs for a list of genes in a grid."""
    has_raw = adata.raw is not None
    var_names_ref = adata.raw.var_names if (use_raw and has_raw) else adata.var_names
    avail = filter_available(genes, var_names_ref)
    if not avail:
        print("[warn] No genes available for feature UMAP")
        return

    nrows = int(np.ceil(len(avail) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 3.5, nrows * 3.2))
    axes = np.array(axes).flatten()

    colorbar_loc = 'right' if show_colorbar else None

    for i, gene in enumerate(avail):
        sc.pl.embedding(
            adata, basis=basis.replace('X_', ''),
            color=gene,
            ax=axes[i],
            show=False,
            use_raw=use_raw and has_raw,
            vmax='p99',
            frameon=False,
            title=gene,
            size=size,
            colorbar_loc=colorbar_loc,
        )

    # Turn off extra axes
    for j in range(len(avail), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(title, fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=DPI, bbox_inches='tight')
        plt.close('all')
        print(f"  Saved: {save_path.name}")
    else:
        plt.close('all')


# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("=" * 80)
print("STEP 1: LOADING STROMAL/VASCULAR DATA")
print("=" * 80)

adata = sc.read_h5ad(INPUT_H5AD)

print(f"  Cells:  {adata.n_obs:,}")
print(f"  Genes:  {adata.n_vars:,}")
if adata.raw is not None:
    print(f"  Raw genes: {adata.raw.n_vars:,}")
print(f"  Obs keys:  {list(adata.obs.columns)}")
print(f"  Obsm keys: {list(adata.obsm.keys())}")

# Print cell type counts
print(f"\n{CELLTYPE_L2} distribution:")
for ct, n in adata.obs[CELLTYPE_L2].value_counts().items():
    print(f"    {ct}: {n:,}")

print(f"\n{CELLTYPE_L3} distribution ({adata.obs[CELLTYPE_L3].nunique()} subclusters):")
for ct, n in adata.obs[CELLTYPE_L3].value_counts().items():
    pct = n / adata.n_obs * 100
    print(f"    {ct}: {n:,} ({pct:.1f}%)")

# Verify / resolve UMAP basis — always normalize to X_umap so downstream is uniform
umap_source = 'X_umap'
if 'X_umap' not in adata.obsm:
    _fallback_order = ['X_umap_scanvi', 'X_umap_scvi', 'umap_scanvi', 'umap_scvi']
    for key in _fallback_order:
        if key in adata.obsm:
            adata.obsm['X_umap'] = adata.obsm[key]
            umap_source = key
            print(f"  [info] X_umap copied from: {key}")
            break
    else:
        print("  [warn] No UMAP found - computing from HVG/PCA")
        sc.pp.highly_variable_genes(
            adata,
            layer='log1p' if 'log1p' in adata.layers else None,
            n_top_genes=3000, subset=False
        )
        sc.pp.pca(adata, use_highly_variable=True, n_comps=50)
        sc.pp.neighbors(adata, n_pcs=50)
        sc.tl.umap(adata, random_state=42)
        umap_source = 'computed_from_HVG_PCA'
        print("  [info] UMAP computed from scratch")
else:
    print("  [info] Using existing X_umap")

# umap_key: scanpy embedding() requires 'umap' not 'X_umap'
umap_key = 'umap'

print(f"  UMAP source: {umap_source}")

gc.collect()

# ============================================================================
# STEP 1b: STANDARDIZE L3 LABELS
# L3 is currently numeric (0, 1, 2 ...) — convert to hierarchical format:
#   {cell_type_L2}_c{subcluster_id}
# e.g. "Endothelia_vascular_venous_systemic_c0"
# Mirrors the R standardization logic in the downstream interpret pipeline.
# ============================================================================

print("\n--- Standardizing cell_type_L3 labels ---")
print(f"  L3 before: {sorted(adata.obs[CELLTYPE_L3].unique().tolist())}")

# Preserve the original numeric subcluster id for sorting / reference
adata.obs['subcluster_id'] = (
    adata.obs[CELLTYPE_L3]
    .astype(str)
    .str.extract(r'^(\d+)$', expand=False)   # only numeric IDs need standardization
)

needs_standardize = adata.obs['subcluster_id'].notna().all()

if needs_standardize:
    # Build hierarchical label
    adata.obs[CELLTYPE_L3] = (
        adata.obs[CELLTYPE_L2].astype(str)
        + '_c'
        + adata.obs['subcluster_id'].astype(str)
    )

    # Order factor: alphabetical L2, then numeric subcluster_id
    tmp = (
        adata.obs[[CELLTYPE_L2, 'subcluster_id', CELLTYPE_L3]]
        .drop_duplicates()
        .assign(subcluster_id=lambda d: d['subcluster_id'].astype(int))
        .sort_values([CELLTYPE_L2, 'subcluster_id'])
    )
    ordered_levels = tmp[CELLTYPE_L3].tolist()
    adata.obs[CELLTYPE_L3] = pd.Categorical(
        adata.obs[CELLTYPE_L3], categories=ordered_levels
    )
    print(f"  L3 after ({adata.obs[CELLTYPE_L3].nunique()} subclusters):")
    for label, n in adata.obs[CELLTYPE_L3].value_counts().sort_index().items():
        pct = n / adata.n_obs * 100
        print(f"    {label}: {n:,} ({pct:.1f}%)")
else:
    print("  L3 labels already hierarchical — skipping standardization")
    print(f"  L3 unique: {sorted(adata.obs[CELLTYPE_L3].unique().tolist())[:10]} ...")

print("[ok] Data loaded")

# ============================================================================
# STEP 2: OVERVIEW UMAP PLOTS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 2: OVERVIEW UMAP PLOTS")
print("=" * 80)

# umap_key is 'umap'; scanpy embedding() prepends 'X_' internally to find obsm['X_umap']

# 2-a  L2 overview
fig, ax = plt.subplots(1, 1, figsize=(8, 7))
sc.pl.embedding(adata, basis=umap_key,
                color=CELLTYPE_L2, ax=ax, show=False,
                legend_loc='right margin', frameon=False,
                size=UMAP_SIZE, title='Stromal/Vascular - Cell Type L2')
plt.tight_layout()
fig.savefig(FIG_DIR / f'01_umap_L2_overview.{FIG_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close('all')
print("  Saved: 01_umap_L2_overview")

# 2-b  L3 subclusters
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
sc.pl.embedding(adata, basis=umap_key,
                color=CELLTYPE_L3, ax=ax, show=False,
                legend_loc='right margin', frameon=False,
                size=UMAP_SIZE, title='Stromal/Vascular - Subclusters L3',
                legend_fontsize=6)
plt.tight_layout()
fig.savefig(FIG_DIR / f'02_umap_L3_subclusters.{FIG_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close('all')
print("  Saved: 02_umap_L3_subclusters")

# 2-c  batch
if BATCH_KEY in adata.obs.columns:
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    sc.pl.embedding(adata, basis=umap_key,
                    color=BATCH_KEY, ax=ax, show=False,
                    legend_loc='right margin', frameon=False,
                    size=UMAP_SIZE, title='Stromal/Vascular - Batch',
                    legend_fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f'03_umap_batch.{FIG_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close('all')
    print("  Saved: 03_umap_batch")

# 2-d  side-by-side L2 + L3  (one combined figure)
fig, axes = plt.subplots(1, 2, figsize=(20, 8))
sc.pl.embedding(adata, basis=umap_key, color=CELLTYPE_L2,
                ax=axes[0], show=False, legend_loc='right margin',
                frameon=False, size=UMAP_SIZE, title='Level 2 (Major Type)')
sc.pl.embedding(adata, basis=umap_key, color=CELLTYPE_L3,
                ax=axes[1], show=False, legend_loc='right margin',
                frameon=False, size=UMAP_SIZE, title='Level 3 (Subcluster)',
                legend_fontsize=6)
plt.suptitle('Stromal/Vascular Cells — Cell Type Overview', fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(FIG_DIR / f'02b_umap_L2_L3_combined.{FIG_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close('all')
print("  Saved: 02b_umap_L2_L3_combined")

# ============================================================================
# STEP 3: DOTPLOTS — ENDOTHELIUM
# ============================================================================

print("\n" + "=" * 80)
print("STEP 3: DOTPLOT - ENDOTHELIUM MARKERS (11 subtypes)")
print("=" * 80)

make_dotplot(
    adata,
    marker_dict=ENDOTHELIAL_MARKERS,
    groupby=CELLTYPE_L3,
    title='Endothelial Subtypes — Discriminatory Marker Panel',
    figsize=(22, 9),
    use_raw=True,
    vmax=3,
    save_path=FIG_DIR / f'04_dotplot_endothelium.{FIG_FORMAT}'
)

# ============================================================================
# STEP 4: DOTPLOTS — FIBROBLAST
# ============================================================================

print("\n" + "=" * 80)
print("STEP 4: DOTPLOT - FIBROBLAST MARKERS (8 subtypes)")
print("=" * 80)

make_dotplot(
    adata,
    marker_dict=FIBROBLAST_MARKERS,
    groupby=CELLTYPE_L3,
    title='Fibroblast Subtypes — Discriminatory Marker Panel',
    figsize=(22, 9),
    use_raw=True,
    vmax=3,
    save_path=FIG_DIR / f'05_dotplot_fibroblast.{FIG_FORMAT}'
)

# ============================================================================
# STEP 5: DOTPLOTS — PERICYTE / SMC / SCHWANN
# ============================================================================

print("\n" + "=" * 80)
print("STEP 5: DOTPLOT - PERICYTE/SMC/SCHWANN MARKERS (7 subtypes)")
print("=" * 80)

make_dotplot(
    adata,
    marker_dict=PERICYTE_SMC_MARKERS,
    groupby=CELLTYPE_L3,
    title='Pericyte / Smooth Muscle / Schwann — Discriminatory Marker Panel',
    figsize=(22, 9),
    use_raw=True,
    vmax=3,
    save_path=FIG_DIR / f'06_dotplot_pericyte_smc.{FIG_FORMAT}'
)

# ============================================================================
# STEP 6: DOTPLOTS — CONTAMINANTS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 6: DOTPLOT - CONTAMINANT MARKERS")
print("=" * 80)

make_dotplot(
    adata,
    marker_dict=CONTAMINANT_MARKERS,
    groupby=CELLTYPE_L3,
    title='Non-stromal Contaminant Markers',
    figsize=(18, 8),
    use_raw=True,
    vmax=3,
    save_path=FIG_DIR / f'07_dotplot_contaminants.{FIG_FORMAT}'
)

# ============================================================================
# STEP 7: GLOBAL OVERVIEW DOTPLOT (minimal panel, all types)
# ============================================================================

print("\n" + "=" * 80)
print("STEP 7: GLOBAL OVERVIEW DOTPLOT")
print("=" * 80)

make_dotplot(
    adata,
    marker_dict=GLOBAL_MINIMAL,
    groupby=CELLTYPE_L3,
    title='Stromal/Vascular — Global Marker Overview (3 genes/type)',
    figsize=(32, 10),
    use_raw=True,
    vmax=3,
    save_path=FIG_DIR / f'08_dotplot_global_overview.{FIG_FORMAT}'
)

# Also make the same dotplot grouped by L2 (coarser but cleaner)
make_dotplot(
    adata,
    marker_dict=GLOBAL_MINIMAL,
    groupby=CELLTYPE_L2,
    title='Stromal/Vascular — Global Marker Overview by L2 Major Type',
    figsize=(30, 7),
    use_raw=True,
    vmax=3,
    save_path=FIG_DIR / f'08b_dotplot_global_L2.{FIG_FORMAT}'
)

# ============================================================================
# STEP 8: FEATURE UMAP — ENDOTHELIAL MARKERS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 8: FEATURE UMAP - ENDOTHELIAL MARKERS")
print("=" * 80)

make_feature_umap(
    adata,
    genes=FEATURE_UMAP_ENDO,
    title='Endothelial Marker Expression — UMAP',
    ncols=4,
    use_raw=True,
    basis='X_umap',
    size=UMAP_SIZE,
    save_path=FIG_DIR / f'09_feature_umap_endothelium.{FIG_FORMAT}'
)

# ============================================================================
# STEP 9: FEATURE UMAP — FIBROBLAST MARKERS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 9: FEATURE UMAP - FIBROBLAST MARKERS")
print("=" * 80)

make_feature_umap(
    adata,
    genes=FEATURE_UMAP_FIBRO,
    title='Fibroblast Marker Expression — UMAP',
    ncols=4,
    use_raw=True,
    basis='X_umap',
    size=UMAP_SIZE,
    save_path=FIG_DIR / f'10_feature_umap_fibroblast.{FIG_FORMAT}'
)

# ============================================================================
# STEP 10: FEATURE UMAP — PERICYTE / SMC / SCHWANN MARKERS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 10: FEATURE UMAP - PERICYTE/SMC/SCHWANN MARKERS")
print("=" * 80)

make_feature_umap(
    adata,
    genes=FEATURE_UMAP_PERI_SMC,
    title='Pericyte / SMC / Schwann Marker Expression — UMAP',
    ncols=4,
    use_raw=True,
    basis='X_umap',
    size=UMAP_SIZE,
    save_path=FIG_DIR / f'11_feature_umap_pericyte_smc.{FIG_FORMAT}'
)

# ============================================================================
# STEP 11: KEY DISCRIMINATORY GENES — ONE-LINER FEATURE UMAP
# ============================================================================

print("\n" + "=" * 80)
print("STEP 11: KEY DISCRIMINATORY FEATURE UMAP (Top differentiators)")
print("=" * 80)

# These 16 genes maximally distinguish all major lineages in one figure
KEY_DISCRIMINATORY = [
    # EC lineage discriminators
    'PROX1',     # Lymphatic EC
    'GJA5',      # Arterial EC
    'GPIHBP1',   # gCap
    'CA4',       # aCap
    'ACKR1',     # Venous EC
    # Fibroblast lineage
    'PI16',      # Adventitial
    'PDGFRA',    # Alveolar
    'POSTN',     # Activated / ECM
    'FABP4',     # Lipofibro
    # Pericyte / SMC
    'RGS5',      # Pericyte
    'CSPG4',     # Pericyte type 2
    'MYH11',     # Vascular SMC
    'PRKG1',     # Airway SMC
    'SOX10',     # Schwann
    # Contamination
    'PTPRC',     # Immune
    'MS4A1',     # B/Plasma
]

make_feature_umap(
    adata,
    genes=KEY_DISCRIMINATORY,
    title='Top Discriminatory Markers — All Stromal/Vascular Lineages',
    ncols=4,
    use_raw=True,
    basis='X_umap',
    size=UMAP_SIZE,
    show_colorbar=True,    # retain colorbar here: dynamic range varies hugely across lineages
    save_path=FIG_DIR / f'12_feature_umap_discriminatory.{FIG_FORMAT}'
)

# ============================================================================
# STEP 12: SUMMARY TABLE — gene availability check
# ============================================================================

print("\n" + "=" * 80)
print("STEP 12: GENE AVAILABILITY SUMMARY")
print("=" * 80)

var_ref = adata.raw.var_names if adata.raw is not None else adata.var_names

all_marker_genes = []
for d in [ENDOTHELIAL_MARKERS, FIBROBLAST_MARKERS, PERICYTE_SMC_MARKERS,
          CONTAMINANT_MARKERS, GLOBAL_MINIMAL]:
    for genes in d.values():
        all_marker_genes.extend(genes)

all_marker_genes = sorted(set(all_marker_genes))
found   = [g for g in all_marker_genes if g in var_ref]
missing = [g for g in all_marker_genes if g not in var_ref]

summary_df = pd.DataFrame({
    'gene': all_marker_genes,
    'in_data': [g in var_ref for g in all_marker_genes]
})
summary_df.to_csv(OUTPUT_DIR / 'gene_availability.csv', index=False)

print(f"  Total marker genes queried: {len(all_marker_genes)}")
print(f"  Found in data:              {len(found)} ({len(found)/len(all_marker_genes)*100:.1f}%)")
print(f"  Not found:                  {len(missing)}")
if missing:
    print(f"  Missing: {missing}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("VISUALIZATION COMPLETE")
print("=" * 80)

figures = sorted(FIG_DIR.glob(f'*.{FIG_FORMAT}'))
print(f"\nOutput figures ({len(figures)} files):")
for f in figures:
    size_kb = f.stat().st_size / 1024
    print(f"  {f.name}  ({size_kb:.0f} KB)")

print(f"\nOutput directory: {OUTPUT_DIR}")
print("""
Next Steps:
  1. Review 08_dotplot_global_overview.pdf  -- high-level confirmation
  2. Check 04-07 dotplots per lineage       -- confirm specific subtypes
  3. Use 12_feature_umap_discriminatory.pdf -- spatial context on UMAP
  4. If a subcluster shows unexpected marker pattern -> re-annotate L3 label
""")
