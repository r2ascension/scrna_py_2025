"""
Title   : B Cell Reference-Only PAGA + Optional Milo Analysis (R-metadata aligned)
Version : 1.3
Date    : 2026-04-13
Author  : r2end

Purpose :
    - PAGA trajectory analysis on the B-cell reference h5ad only
    - Harmonize metadata to match `bcell_tissue_comparison_v2_6_1_20260410.R`
    - Optional Milo (milopy) neighborhood-level differential abundance testing
    - Automatically select a compatible latent/UMAP embedding from the reference object

Input   :
    - bcell_reference_20260203.h5ad  (reference only)

Output  :
    - adata_paga_milo_reference_only_rmeta_v1_3.h5ad
    - PAGA graph figures
    - metadata harmonization summary tables
    - Optional Milo figures/results if a valid condition column is configured

Memory  : ~8-15 GB peak
Runtime : ~10-30 min (CPU only, depending on Milo)
"""

# %% =========================================================================
# Section 1: Imports
# =============================================================================

import sys
import os
import gc
import time
import warnings
import json
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import scanpy as sc
import anndata as ad

# milopy for Milo differential abundance
try:
    import milopy
    import milopy.core as milo
    import milopy.utils as milo_utils
    MILOPY_AVAILABLE = True
    print(f"milopy : {milopy.__version__}")
except ImportError:
    MILOPY_AVAILABLE = False
    print("[WARN] milopy not found. Install: pip install milopy")
    print("       Milo sections will be skipped.")

warnings.filterwarnings("ignore")
PIPELINE_START = time.time()

print("=" * 80)
print("B Cell PAGA + Optional Milo Analysis v1.3")
print("=" * 80)
print(f"scanpy : {sc.__version__}")
print(f"anndata: {ad.__version__}")
print(f"Python : {sys.version}")


# %% =========================================================================
# Section 2: Configuration  --  EDIT THIS SECTION
# =============================================================================

# ----- Input paths -----
PATH_REFERENCE = (
    "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/"
    "bcell_reference_20260203.h5ad"
)

PATH_OUT = "/home/h2048/data/py/0413/bcell_paga_reference_rmeta_v1_3"

# ----- Key column names -----
BATCH_KEY      = "sample"             # sample/replicate column (for Milo)
SUBCLUSTER_KEY = "cell_type_L3"       # standardized fine-grained annotation (R metadata)
L2_KEY         = "cell_type_L2"       # standardized coarse annotation (R metadata)
SOURCE_KEY     = "data_source"        # kept for compatibility if present
TISSUE_KEY     = "tissue"             # tissue column

# ----- R script metadata harmonization -----
L3_SOURCE_COL = "cell_type_scanvi_pred"
USE_EXISTING_L2 = False
L2_SOURCE_COL = "cell_type_L2"
L3_TO_L2_REMAP = {
    "Atypical_Memory_B": "Memory_B",
    "IGHEplus_Atypical_Memory_B": "Memory_B",
    "Memory_B": "Memory_B",
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC_B",
    "GC_B_Light_Zone_Centrocyte": "GC_B",
    "GC_B_Transitional": "GC_B",
    "Plasma_IgA": "Plasma",
    "Plasma_IgG": "Plasma",
    "Naive_B": "Naive_B",
}

# Condition column for Milo DA test.
# Must be a column in obs with >=2 biologically meaningful levels.
# Reference h5ad's 'status' is mostly 'Unknown', so Milo is disabled by default.
# Set CONDITION_KEY manually (e.g. a curated disease/group column) to enable Milo.
CONDITION_KEY  = None
CONDITION_REF  = None

# Preferred embeddings to use for PAGA + Milo neighborhoods / plotting.
LATENT_KEY_CANDIDATES = ["X_scanvi_corrected", "X_scanvi", "X_scANVI_L2", "X_scvi", "X_harmony"]
UMAP_KEY_CANDIDATES   = ["X_umap_scanvi_corrected", "X_umap_scanvi", "X_umap", "X_umap_scvi"]

# ----- Which cells to include -----
# reference-only script: source subsetting is only used if such a column exists.
SOURCE_SUBSET  = "all"

# ----- PAGA parameters -----
PAGA_N_NEIGHBORS = 30                 # neighbors computed on LATENT_KEY
PAGA_RESOLUTION  = 0.5               # Leiden used internally by PAGA

# ----- Milo parameters -----
MILO_PROP        = 0.1               # proportion of cells used as index cells
MILO_K           = 30                # k for KNN nhood graph
MILO_D           = 30                # n_pcs (latent dims) for nhoods
MILO_SPATIAL_FDR = 0.1              # significance threshold

# ----- Visualization -----
DPI           = 300
FIG_FORMAT    = "pdf"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
sc.settings.set_figure_params(dpi=DPI, facecolor="white",
                               format=FIG_FORMAT, vector_friendly=False)

# ----- Reproducibility -----
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# ----- Output dirs -----
out_dir  = Path(PATH_OUT)
fig_dir  = out_dir / "figures"
for _d in [out_dir, fig_dir]:
    _d.mkdir(parents=True, exist_ok=True)

print(f"\nReference input : {PATH_REFERENCE}")
print(f"Output       : {PATH_OUT}")
print(f"Latent key candidates : {LATENT_KEY_CANDIDATES}")
print(f"UMAP key candidates   : {UMAP_KEY_CANDIDATES}")
print(f"PAGA groupby : {SUBCLUSTER_KEY}")
print(f"R metadata L3 source : {L3_SOURCE_COL}")
print(f"Condition    : {CONDITION_KEY}  (ref={CONDITION_REF})")
print(f"Source subset: {SOURCE_SUBSET}")


def resolve_first_available_key(available_keys, candidates, label):
    """Return the first candidate that exists; raise a helpful error otherwise."""
    for key in candidates:
        if key in available_keys:
            print(f"[OK] Using {label}: {key}")
            return key
    raise KeyError(
        f"No compatible {label} found. Candidates={candidates}; available={list(available_keys)}"
    )


def harmonize_metadata_from_r_script(adata_obj):
    """Apply the same metadata standardization used by the R tissue comparison pipeline."""
    required_cols = [TISSUE_KEY, BATCH_KEY, L3_SOURCE_COL]
    if USE_EXISTING_L2:
        required_cols.append(L2_SOURCE_COL)
    missing = [col for col in required_cols if col not in adata_obj.obs.columns]
    if missing:
        raise KeyError(f"Missing required metadata columns for R harmonization: {missing}")

    l3_vals = adata_obj.obs[L3_SOURCE_COL].astype("string").fillna("").str.strip()
    if USE_EXISTING_L2:
        l2_vals = adata_obj.obs[L2_SOURCE_COL].astype("string").fillna("").str.strip()
        if (l2_vals == "").all():
            raise ValueError(f"Existing L2 column '{L2_SOURCE_COL}' is empty.")
        print(f"[OK] Using existing L2 from '{L2_SOURCE_COL}'")
    else:
        l2_vals = l3_vals.map(L3_TO_L2_REMAP)
        unmapped = sorted(pd.unique(l3_vals[l2_vals.isna() & (l3_vals != "")]))
        if unmapped:
            raise ValueError(
                f"Unmapped L3 values in '{L3_SOURCE_COL}': {unmapped}. Update L3_TO_L2_REMAP."
            )
        backup_l2_col = f"{L2_KEY}_input"
        if L2_KEY in adata_obj.obs.columns and backup_l2_col not in adata_obj.obs.columns:
            adata_obj.obs[backup_l2_col] = adata_obj.obs[L2_KEY].astype("string")
            print(f"[OK] Preserved original '{L2_KEY}' as '{backup_l2_col}'")

    backup_l3_col = f"{SUBCLUSTER_KEY}_input"
    if SUBCLUSTER_KEY in adata_obj.obs.columns and backup_l3_col not in adata_obj.obs.columns:
        adata_obj.obs[backup_l3_col] = adata_obj.obs[SUBCLUSTER_KEY].astype("string")
        print(f"[OK] Preserved original '{SUBCLUSTER_KEY}' as '{backup_l3_col}'")

    adata_obj.obs[L2_KEY] = pd.Categorical(l2_vals)
    adata_obj.obs[SUBCLUSTER_KEY] = pd.Categorical(l3_vals)

    summary_df = (
        adata_obj.obs[[SUBCLUSTER_KEY, L2_KEY]]
        .astype("string")
        .value_counts(dropna=True)
        .rename("Freq")
        .reset_index()
        .sort_values([SUBCLUSTER_KEY, "Freq", L2_KEY], ascending=[True, False, True])
    )

    bad_idx = (
        adata_obj.obs[TISSUE_KEY].astype("string").fillna("").str.strip().eq("")
        | adata_obj.obs[BATCH_KEY].astype("string").fillna("").str.strip().eq("")
        | adata_obj.obs[L2_KEY].astype("string").fillna("").str.strip().eq("")
        | adata_obj.obs[SUBCLUSTER_KEY].astype("string").fillna("").str.strip().eq("")
    )
    n_bad = int(bad_idx.sum())
    if n_bad > 0:
        print(f"[INFO] Dropping {n_bad} cells with NA/blank tissue/sample/L2/L3")
        adata_obj._inplace_subset_obs((~bad_idx).to_numpy())

    print("\nR-metadata L3 / L2 summary:")
    print(summary_df.to_string(index=False))
    print("\nR-metadata L2 distribution:")
    print(adata_obj.obs[L2_KEY].astype("string").value_counts().to_string())
    print("\nR-metadata L3 distribution:")
    print(adata_obj.obs[SUBCLUSTER_KEY].astype("string").value_counts().to_string())
    return summary_df


def align_umap_for_plotting(adata_obj, umap_key):
    """Make the chosen UMAP visible to sc.pl.umap via adata.obsm['X_umap']."""
    if umap_key != "X_umap":
        adata_obj.obsm["X_umap"] = adata_obj.obsm[umap_key].copy()
        print(f"[INFO] Copied {umap_key} -> X_umap for plotting")


# %% =========================================================================
# Section 3: Load Reference h5ad
# =============================================================================

print("\n" + "=" * 80)
print("STEP 1: LOAD REFERENCE H5AD")
print("=" * 80)

p = Path(PATH_REFERENCE)
if not p.exists():
    raise FileNotFoundError(f"[ERROR] Not found: {p}")

adata = sc.read_h5ad(p)
adata.var_names_make_unique()
adata.obs_names_make_unique()

print(f"Loaded shape : {adata.shape}")
print(f"obs cols     : {list(adata.obs.columns)}")
print(f"obsm keys    : {list(adata.obsm.keys())}")

LATENT_KEY = resolve_first_available_key(adata.obsm.keys(), LATENT_KEY_CANDIDATES, "latent embedding")
UMAP_KEY = resolve_first_available_key(adata.obsm.keys(), UMAP_KEY_CANDIDATES, "UMAP embedding")
align_umap_for_plotting(adata, UMAP_KEY)

# Source distribution
if SOURCE_KEY in adata.obs.columns:
    print(f"\nSource dist  : {adata.obs[SOURCE_KEY].value_counts().to_dict()}")
else:
    print("\n[INFO] No data_source column detected; treating input as pure reference object")


# %% =========================================================================
# Section 4: Harmonize Metadata from R Tissue Comparison Script
# =============================================================================

print("\n" + "=" * 80)
print("STEP 2: HARMONIZE R-SCRIPT METADATA")
print("=" * 80)

mapping_summary_df = harmonize_metadata_from_r_script(adata)
mapping_summary_df.to_csv(out_dir / "l3_l2_mapping_summary.csv", index=False)
GROUPBY_KEY = SUBCLUSTER_KEY
print(f"[OK] Using harmonized '{GROUPBY_KEY}' for PAGA groupby")

print(f"\n{GROUPBY_KEY} distribution:")
for lbl, n in adata.obs[GROUPBY_KEY].value_counts().head(12).items():
    print(f"  {lbl}: {n:,}")


# %% =========================================================================
# Section 5: Subset by Source (if requested)
# =============================================================================

print("\n" + "=" * 80)
print("STEP 3: SUBSET BY SOURCE")
print("=" * 80)

if SOURCE_SUBSET != "all" and SOURCE_KEY in adata.obs.columns:
    n_before = adata.n_obs
    adata = adata[adata.obs[SOURCE_KEY] == SOURCE_SUBSET].copy()
    gc.collect()
    print(f"Kept '{SOURCE_SUBSET}': {n_before:,} -> {adata.n_obs:,} cells")
elif SOURCE_KEY not in adata.obs.columns:
    print("Reference-only input has no data_source column; skipping source subset")
else:
    print(f"Using all cells: {adata.n_obs:,}")

# Drop cells with missing GROUPBY_KEY
n_before = adata.n_obs
adata = adata[adata.obs[GROUPBY_KEY].notna()].copy()
if adata.n_obs < n_before:
    print(f"Dropped {n_before - adata.n_obs:,} cells with NA in '{GROUPBY_KEY}'")

# Validate condition column for Milo
MILO_READY = False
if CONDITION_KEY is None:
    print("[INFO] CONDITION_KEY=None. Milo disabled by default for reference-only mode.")
elif MILOPY_AVAILABLE and CONDITION_KEY is not None:
    if CONDITION_KEY not in adata.obs.columns:
        print(f"[WARN] CONDITION_KEY='{CONDITION_KEY}' not in obs. Milo will be skipped.")
        print(f"       Available: {list(adata.obs.columns)}")
    else:
        cond_vals = adata.obs[CONDITION_KEY].dropna().unique()
        print(f"\nCondition '{CONDITION_KEY}' levels: {list(cond_vals)}")
        if len(cond_vals) < 2:
            print("[WARN] <2 condition levels. Milo needs >=2 groups.")
        elif CONDITION_REF is None:
            print("[WARN] CONDITION_REF is None. Set a reference level to enable Milo.")
        elif CONDITION_REF not in cond_vals:
            print(f"[WARN] CONDITION_REF='{CONDITION_REF}' not found in condition levels.")
        else:
            MILO_READY = True
            print(f"[OK] Milo ready. Contrast: ~{CONDITION_KEY}  ref={CONDITION_REF}")
elif not MILOPY_AVAILABLE:
    print("[INFO] milopy not available; PAGA will still run.")

# Validate batch column
if BATCH_KEY not in adata.obs.columns:
    print(f"[WARN] BATCH_KEY='{BATCH_KEY}' not in obs.")
    n_samples = 0
else:
    n_samples = adata.obs[BATCH_KEY].nunique()
    print(f"\nBatch '{BATCH_KEY}': {n_samples} unique samples")

print(f"\nFinal adata for analysis: {adata.shape}")


# %% =========================================================================
# Section 6: Compute Neighbors on Latent Space (shared for PAGA + Milo)
# =============================================================================

print("\n" + "=" * 80)
print("STEP 4: COMPUTE NEIGHBORS ON LATENT SPACE")
print("=" * 80)

# Use the selected latent embedding directly as pre-computed latent; skip PCA
NEIGHBORS_KEY = "neighbors_milo"

print(f"  use_rep='{LATENT_KEY}'  n_neighbors={PAGA_N_NEIGHBORS}")
sc.pp.neighbors(
    adata,
    use_rep=LATENT_KEY,
    n_neighbors=PAGA_N_NEIGHBORS,
    random_state=RANDOM_SEED,
    key_added=NEIGHBORS_KEY,
)
print(f"  [OK] Neighbors computed (key='{NEIGHBORS_KEY}')")

# Copy to default neighbors slot for PAGA (sc.tl.paga uses uns['neighbors'])
adata.uns["neighbors"]  = adata.uns[NEIGHBORS_KEY].copy()
adata.obsp["distances"] = adata.obsp[f"{NEIGHBORS_KEY}_distances"]
adata.obsp["connectivities"] = adata.obsp[f"{NEIGHBORS_KEY}_connectivities"]
print("  [OK] Default neighbors slot updated for PAGA")

# Compute UMAP on current neighbors if X_umap already exists (skip recompute)
if "X_umap" not in adata.obsm:
    print("  Computing UMAP (X_umap not found)...")
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    print("  [OK] UMAP computed")
else:
    print(f"  [INFO] Using plotting UMAP '{UMAP_KEY}': {adata.obsm['X_umap'].shape}")


# %% =========================================================================
# Section 7: PAGA — Partition-based Graph Abstraction
# =============================================================================

print("\n" + "=" * 80)
print("STEP 5: PAGA TRAJECTORY ANALYSIS")
print("=" * 80)

# Make sure GROUPBY_KEY is categorical
adata.obs[GROUPBY_KEY] = pd.Categorical(adata.obs[GROUPBY_KEY])

print(f"  Groups for PAGA: {adata.obs[GROUPBY_KEY].nunique()} categories")

# Leiden clustering (PAGA needs uns['leiden'] for internal use; we compute on
# the groupby column directly instead)
# sc.tl.paga accepts any categorical obs column directly.
sc.tl.paga(adata, groups=GROUPBY_KEY)
print(f"  [OK] PAGA computed on '{GROUPBY_KEY}'")

# Save PAGA connectivity matrix
paga_conn = adata.uns["paga"]["connectivities"].toarray()
paga_cats = list(adata.obs[GROUPBY_KEY].cat.categories)
df_paga = pd.DataFrame(paga_conn, index=paga_cats, columns=paga_cats)
df_paga.to_csv(out_dir / "paga_connectivities.csv")
print(f"  [OK] paga_connectivities.csv  shape={df_paga.shape}")

# ----- Figure 1: PAGA graph layout -----
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.paga(
    adata,
    ax=ax,
    show=False,
    threshold=0.05,              # only show edges above this connectivity
    node_size_scale=1.2,
    edge_width_scale=1.5,
    fontsize=8,
    title=f"PAGA — B Cell ({GROUPBY_KEY})",
    frameon=False,
)
plt.tight_layout()
fig.savefig(fig_dir / f"01_paga_graph.{FIG_FORMAT}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] 01_paga_graph.{FIG_FORMAT}")

# ----- Figure 2: UMAP with PAGA overlay -----
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.paga(
    adata,
    ax=ax,
    show=False,
    threshold=0.05,
    title=f"PAGA + UMAP — {GROUPBY_KEY}",
    frameon=False,
    plot=False,         # compute positions only
)
# PAGA positions stored in uns["paga"]["pos"]
paga_pos = adata.uns["paga"]["pos"]   # shape (n_cats, 2)

# Scatter cells colored by GROUPBY_KEY
sc.pl.umap(
    adata,
    color=GROUPBY_KEY,
    ax=ax,
    show=False,
    legend_loc="right margin",
    legend_fontsize=7,
    frameon=False,
    s=10,
    alpha=0.6,
    title=f"UMAP with PAGA overlay — {GROUPBY_KEY}",
)
# Overlay PAGA edges on top
for i, lbl_i in enumerate(paga_cats):
    for j, lbl_j in enumerate(paga_cats):
        if i >= j:
            continue
        w = paga_conn[i, j]
        if w > 0.05:
            x = [paga_pos[i, 0], paga_pos[j, 0]]
            y = [paga_pos[i, 1], paga_pos[j, 1]]
            ax.plot(x, y, "k-", linewidth=w * 4, alpha=0.5, zorder=3)
for i, lbl in enumerate(paga_cats):
    ax.scatter(paga_pos[i, 0], paga_pos[i, 1],
               s=80, zorder=4, color="white", edgecolors="black", linewidth=0.8)

plt.tight_layout()
fig.savefig(fig_dir / f"02_paga_umap_overlay.{FIG_FORMAT}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] 02_paga_umap_overlay.{FIG_FORMAT}")

# ----- Figure 3: PAGA heatmap (connectivity matrix) -----
fig, ax = plt.subplots(figsize=(max(8, len(paga_cats) * 0.6),
                                max(6, len(paga_cats) * 0.5)))
mask = np.eye(len(paga_cats), dtype=bool)
sns.heatmap(
    df_paga,
    ax=ax,
    cmap="Blues",
    annot=(len(paga_cats) <= 20),
    fmt=".2f",
    linewidths=0.4,
    mask=mask,
    cbar_kws={"label": "PAGA Connectivity"},
)
ax.set_title(f"PAGA Connectivity Matrix — {GROUPBY_KEY}", fontsize=12, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=7)
plt.yticks(rotation=0, fontsize=7)
plt.tight_layout()
fig.savefig(fig_dir / f"03_paga_connectivity_heatmap.{FIG_FORMAT}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] 03_paga_connectivity_heatmap.{FIG_FORMAT}")

# ----- Figure 4: UMAP coloured by key obs columns -----
color_cols = [GROUPBY_KEY, L2_KEY, TISSUE_KEY, BATCH_KEY, SOURCE_KEY]
color_cols = [c for c in color_cols if c in adata.obs.columns]
if CONDITION_KEY and CONDITION_KEY in adata.obs.columns:
    color_cols.append(CONDITION_KEY)

n_panels = len(color_cols)
ncols = min(n_panels, 3)
nrows = int(np.ceil(n_panels / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
axes_flat = np.array(axes).flatten() if n_panels > 1 else [axes]

for i, col in enumerate(color_cols):
    sc.pl.umap(
        adata,
        color=col,
        ax=axes_flat[i],
        show=False,
        legend_loc="right margin",
        legend_fontsize=7,
        frameon=False,
        s=10,
        title=col,
    )
for j in range(i + 1, len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.suptitle("B Cell UMAP — Key Metadata", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(fig_dir / f"04_umap_metadata.{FIG_FORMAT}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] 04_umap_metadata.{FIG_FORMAT}")


# %% =========================================================================
# Section 8: Milo — Differential Abundance Testing
# =============================================================================

print("\n" + "=" * 80)
print("STEP 6: MILO DIFFERENTIAL ABUNDANCE")
print("=" * 80)

if not MILO_READY:
    print("[SKIP] Milo conditions not met (see warnings above).")
    print("       Set CONDITION_KEY and CONDITION_REF correctly in Section 2.")
else:
    # ---- 8a. Make neighborhoods ----
    print(f"  make_nhoods: prop={MILO_PROP}  k={MILO_K}  d={MILO_D}")
    milo.make_nhoods(adata, prop=MILO_PROP)
    print(f"  [OK] Neighborhoods: {adata.obsm['nhoods'].shape[1]:,} nhoods")

    # ---- 8b. Count cells per sample per nhood ----
    print(f"  count_nhoods: sample_col='{BATCH_KEY}'")
    milo.count_nhoods(adata, sample_col=BATCH_KEY)
    print(f"  [OK] nhoodCounts shape: {adata.uns['nhood_adata'].X.shape}")

    # ---- 8c. Build sample-level metadata for design ----
    # One row per sample; must have CONDITION_KEY and BATCH_KEY
    sample_meta = (
        adata.obs[[BATCH_KEY, CONDITION_KEY]]
        .drop_duplicates(subset=BATCH_KEY)
        .set_index(BATCH_KEY)
        .sort_index()
    )
    # Add TISSUE_KEY if available
    if TISSUE_KEY in adata.obs.columns:
        tissue_meta = (
            adata.obs[[BATCH_KEY, TISSUE_KEY]]
            .drop_duplicates(subset=BATCH_KEY)
            .set_index(BATCH_KEY)
        )
        sample_meta = sample_meta.join(tissue_meta, how="left")

    adata.uns["nhood_adata"].obs = (
        adata.uns["nhood_adata"].obs
        .join(sample_meta, on=BATCH_KEY, how="left")
    )

    print(f"\n  Sample metadata ({sample_meta.shape[0]} samples):")
    print(f"  Condition distribution: {sample_meta[CONDITION_KEY].value_counts().to_dict()}")
    sample_meta.to_csv(out_dir / "milo_sample_metadata.csv")

    # ---- 8d. DA testing ----
    # Design formula: ~CONDITION_KEY
    # If confounders are needed (e.g., tissue), extend design below.
    design_formula = f"~{CONDITION_KEY}"
    print(f"\n  DA_nhoods: design='{design_formula}'")

    milo.DA_nhoods(
        adata,
        design=design_formula,
    )

    nhood_adata = adata.uns["nhood_adata"]
    print(f"\n  DA results (nhood_adata.obs columns): {list(nhood_adata.obs.columns)}")

    # Number of significant nhoods
    sig_col = "SpatialFDR"
    if sig_col in nhood_adata.obs.columns:
        n_sig  = (nhood_adata.obs[sig_col] < MILO_SPATIAL_FDR).sum()
        n_tot  = nhood_adata.obs.shape[0]
        print(f"  Significant nhoods (SpatialFDR<{MILO_SPATIAL_FDR}): {n_sig}/{n_tot}")
        print(f"  logFC range: "
              f"{nhood_adata.obs['logFC'].min():.2f} to "
              f"{nhood_adata.obs['logFC'].max():.2f}")

    # ---- 8e. Annotate nhoods with cell type label ----
    print(f"\n  Annotating nhoods with '{GROUPBY_KEY}'...")
    milo_utils.annotate_nhoods(adata, anno_col=GROUPBY_KEY)
    print(f"  [OK] nhood_annotation column added")

    # ---- 8f. Save DA results ----
    da_results = nhood_adata.obs.copy()
    da_results.to_csv(out_dir / "milo_DA_results.csv")
    print(f"  [OK] milo_DA_results.csv  ({da_results.shape[0]} nhoods)")

    # ---- 8g. Visualize: Beeswarm ----
    print(f"\n  Plotting Milo beeswarm...")
    try:
        fig, ax = plt.subplots(figsize=(14, 8))
        milo_utils.plot_nhood_graph(
            adata,
            alpha=MILO_SPATIAL_FDR,
            min_size=2,
        )
        plt.suptitle(
            f"Milo DA — B Cells ({CONDITION_KEY}: {CONDITION_REF} vs others)",
            fontsize=12, fontweight="bold"
        )
        plt.tight_layout()
        fig.savefig(fig_dir / f"05_milo_nhood_graph.{FIG_FORMAT}",
                    dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  [OK] 05_milo_nhood_graph.{FIG_FORMAT}")
    except Exception as e:
        print(f"  [WARN] plot_nhood_graph failed: {e}")

    # ---- 8h. Beeswarm per cell type ----
    try:
        from milopy.utils import plot_DA_beeswarm

        fig, ax = plt.subplots(figsize=(10, 8))
        plot_DA_beeswarm(
            adata,
            alpha=MILO_SPATIAL_FDR,
            ax=ax,
        )
        ax.set_title(
            f"Milo DA Beeswarm — {GROUPBY_KEY}  (SpatialFDR<{MILO_SPATIAL_FDR})",
            fontsize=12, fontweight="bold"
        )
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        plt.tight_layout()
        fig.savefig(fig_dir / f"06_milo_beeswarm.{FIG_FORMAT}",
                    dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  [OK] 06_milo_beeswarm.{FIG_FORMAT}")
    except Exception as e:
        print(f"  [WARN] plot_DA_beeswarm failed: {e}")
        # Manual beeswarm fallback
        if "nhood_annotation" in nhood_adata.obs.columns and sig_col in nhood_adata.obs.columns:
            df_bee = nhood_adata.obs[["logFC", sig_col, "nhood_annotation"]].copy()
            df_bee["significant"] = df_bee[sig_col] < MILO_SPATIAL_FDR
            ct_order = (
                df_bee.groupby("nhood_annotation")["logFC"]
                .median().sort_values().index
            )
            fig, ax = plt.subplots(figsize=(8, max(6, len(ct_order) * 0.5)))
            palette = {True: "#E64B35", False: "lightgrey"}
            sns.stripplot(
                data=df_bee,
                y="nhood_annotation",
                x="logFC",
                hue="significant",
                palette=palette,
                order=ct_order,
                dodge=False,
                size=4,
                alpha=0.7,
                ax=ax,
            )
            ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
            ax.set_xlabel("log2 Fold Change")
            ax.set_ylabel("Cell Type")
            ax.set_title(
                f"Milo DA — {CONDITION_KEY}  (SpatialFDR<{MILO_SPATIAL_FDR} in red)",
                fontsize=11, fontweight="bold"
            )
            ax.legend(title="Significant", loc="lower right")
            ax.grid(axis="x", alpha=0.3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            plt.tight_layout()
            fig.savefig(fig_dir / f"06_milo_beeswarm_manual.{FIG_FORMAT}",
                        dpi=DPI, bbox_inches="tight")
            plt.close(fig)
            print(f"  [OK] 06_milo_beeswarm_manual.{FIG_FORMAT} (fallback)")

    # ---- 8i. Volcano plot: logFC vs -log10(SpatialFDR) ----
    if sig_col in nhood_adata.obs.columns:
        df_vol = nhood_adata.obs[["logFC", sig_col]].copy()
        df_vol["-log10_fdr"] = -np.log10(df_vol[sig_col].clip(lower=1e-300))
        df_vol["significant"] = df_vol[sig_col] < MILO_SPATIAL_FDR

        fig, ax = plt.subplots(figsize=(10, 7))
        colors = df_vol["significant"].map({True: "#E64B35", False: "lightgrey"})
        ax.scatter(df_vol["logFC"], df_vol["-log10_fdr"],
                   c=colors, s=15, alpha=0.6, rasterized=True)
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.axhline(-np.log10(MILO_SPATIAL_FDR), color="red",
                   linewidth=0.8, linestyle="--",
                   label=f"SpatialFDR={MILO_SPATIAL_FDR}")
        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10(SpatialFDR)")
        ax.set_title(f"Milo Volcano — {CONDITION_KEY}", fontsize=12, fontweight="bold")
        ax.legend()
        ax.grid(alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(fig_dir / f"07_milo_volcano.{FIG_FORMAT}",
                    dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  [OK] 07_milo_volcano.{FIG_FORMAT}")

    # ---- 8j. UMAP colored by Milo logFC ----
    # Transfer per-cell logFC (mean logFC of nhoods the cell belongs to)
    nhoods_matrix = adata.obsm["nhoods"]   # sparse (cells x nhoods)
    logfc_arr     = nhood_adata.obs["logFC"].values.astype(np.float32)
    fdr_arr       = nhood_adata.obs[sig_col].values.astype(np.float32) if sig_col in nhood_adata.obs.columns else np.ones(len(logfc_arr))

    # Mean logFC across all nhoods each cell belongs to
    with np.errstate(invalid="ignore"):
        nhoods_dense = nhoods_matrix.toarray() if issparse(nhoods_matrix) else nhoods_matrix
        cell_membership_count = nhoods_dense.sum(axis=1)
        cell_logfc = nhoods_dense.dot(logfc_arr) / np.maximum(cell_membership_count, 1)
        cell_logfc[cell_membership_count == 0] = np.nan

    adata.obs["milo_logFC"] = cell_logfc

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    sc.pl.umap(adata, color=GROUPBY_KEY,
               ax=axes[0], show=False, legend_loc="right margin",
               legend_fontsize=7, frameon=False, s=8,
               title=f"Cell Type ({GROUPBY_KEY})")
    sc.pl.umap(adata, color="milo_logFC",
               ax=axes[1], show=False, frameon=False, s=8,
               cmap="RdBu_r", vcenter=0,
               title=f"Milo logFC (mean over nhoods)")
    plt.suptitle(f"Milo DA — {CONDITION_KEY}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(fig_dir / f"08_milo_logfc_umap.{FIG_FORMAT}",
                dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] 08_milo_logfc_umap.{FIG_FORMAT}")


# %% =========================================================================
# Section 9: Save h5ad
# =============================================================================

print("\n" + "=" * 80)
print("STEP 7: SAVE H5AD")
print("=" * 80)

# Sanitize boolean/object columns before write
for col in adata.obs.columns:
    s = adata.obs[col]
    if s.dtype == object:
        non_null = s.dropna()
        if len(non_null) == 0:
            continue
        py_types = set(non_null.map(lambda x: type(x).__name__))
        if py_types <= {"bool"}:
            adata.obs[col] = s.fillna(False).astype(bool)
        elif py_types <= {"str"}:
            adata.obs[col] = s.where(s.notna(), "").astype(str)

import anndata
try:
    anndata.settings.allow_write_nullable_strings = True
except AttributeError:
    pass

out_h5ad = out_dir / "adata_paga_milo_reference_only_v1_2.h5ad"
out_h5ad = out_dir / "adata_paga_milo_reference_only_rmeta_v1_3.h5ad"
adata.write_h5ad(out_h5ad, compression="gzip")
print(f"  [OK] {out_h5ad}  ({out_h5ad.stat().st_size / 1024**3:.2f} GB)")


# %% =========================================================================
# Section 10: Final Summary
# =============================================================================

elapsed = (time.time() - PIPELINE_START) / 60

config = {
    "pipeline":        "bcell_paga_milo_reference_only_rmeta",
    "version":         "1.3",
    "timestamp":       datetime.now().isoformat(),
    "input_reference": PATH_REFERENCE,
    "output_dir":      PATH_OUT,
    "latent_key":      LATENT_KEY,
    "umap_key":        UMAP_KEY,
    "groupby_key":     GROUPBY_KEY,
    "batch_key":       BATCH_KEY,
    "tissue_key":      TISSUE_KEY,
    "l3_source_col":   L3_SOURCE_COL,
    "l2_key":          L2_KEY,
    "metadata_aligned_to_r_script": True,
    "l3_to_l2_remap":  L3_TO_L2_REMAP,
    "condition_key":   CONDITION_KEY,
    "condition_ref":   CONDITION_REF,
    "source_subset":   SOURCE_SUBSET,
    "milo_ready":      MILO_READY,
    "paga": {
        "n_neighbors":  PAGA_N_NEIGHBORS,
        "n_groups":     adata.obs[GROUPBY_KEY].nunique(),
        "threshold":    0.05,
    },
    "milo": {
        "prop":         MILO_PROP,
        "k":            MILO_K,
        "d":            MILO_D,
        "spatial_fdr":  MILO_SPATIAL_FDR,
        "design":       f"~{CONDITION_KEY}" if MILO_READY else None,
    } if MILO_READY else None,
}
cfg_path = out_dir / "config_v1_3.json"
with open(cfg_path, "w") as f:
    json.dump(config, f, indent=2, default=str)

print("\n" + "=" * 80)
print("FINAL SUMMARY — B Cell PAGA + Optional Milo v1.3 (R-metadata aligned)")
print("=" * 80)
print(f"Runtime       : {elapsed:.1f} min")
print(f"Final shape   : {adata.shape}")
print(f"PAGA groups   : {adata.obs[GROUPBY_KEY].nunique()}  ({GROUPBY_KEY})")
print(f"PAGA L2 groups: {adata.obs[L2_KEY].nunique()}  ({L2_KEY})")
if MILO_READY:
    n_sig_final = (
        (adata.uns["nhood_adata"].obs["SpatialFDR"] < MILO_SPATIAL_FDR).sum()
        if "SpatialFDR" in adata.uns["nhood_adata"].obs.columns else "N/A"
    )
    print(f"Milo nhoods   : {adata.uns['nhood_adata'].X.shape[1]:,}")
    print(f"Significant   : {n_sig_final}  (SpatialFDR<{MILO_SPATIAL_FDR})")
print(f"\nOutputs:")
print(f"  {out_h5ad}")
print(f"  {out_dir/'paga_connectivities.csv'}")
print(f"  {out_dir/'l3_l2_mapping_summary.csv'}")
print(f"  {fig_dir}/  ({len(list(fig_dir.glob('*')))} figures)")
if MILO_READY:
    print(f"  {out_dir/'milo_DA_results.csv'}")
    print(f"  {out_dir/'milo_sample_metadata.csv'}")
print(f"  {cfg_path}")
print("\n" + "=" * 80)
print("SUCCESS")
print("=" * 80)
