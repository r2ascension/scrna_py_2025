# %% [markdown]
# # Epithelial Cell Subcluster Analysis (Post-BBKNN Integration)
#
# **Purpose:**
# - Perform subcluster analysis on already-integrated epithelial cells
# - Visualize subclusters on original UMAP (preserves global structure)
# - Identify marker genes with filtering of technical artifacts
#
# **Author:** r2end
# **Date:** 2025-01-20
# **Version:** 4.1 - RESTRICT_TO Method (TRUE Post-BBKNN) + HOTFIXES
#
# **Key Features:**
# - ⭐ Uses global BBKNN graph via `restrict_to` parameter
# - Does NOT rebuild neighbors graph (preserves batch correction)
# - True post-integration subclustering (no batch effect re-introduction)
# - Original UMAP preservation
# - Marker gene filtering (MT/Ribo/Histone/Unannotated + Dissociation/IEG-only)
# - Reproducible (random seeds)
#
# **HOTFIXES v4.1:**
# - P0: Lock gene-universe sources (raw/var) and write flags safely (no mismatch)
# - P0: Stable `rank_genes_groups` result retrieval (no repeated slicing)
# - P1: Apply `MARKER_MIN_PCT` + delta pct filter (reduces low-expression noise)
# - P1: DE input sanity check (raw log1p vs counts warning + auto-selection)
# - P2: Version consistency + numeric sorting for subcluster ids + remove emoji encoding noise

# %% [markdown]
# ## Configuration Parameters

# %%
import os
import sys
import gc
import warnings
from pathlib import Path
from datetime import datetime
import time
import random
import re

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# =============================================================================
# PIPELINE METADATA
# =============================================================================
PIPELINE_VERSION = "v4.2-RESTRICT_TO-HOTFIX"

# =============================================================================
# GLOBAL NEIGHBORS GRAPH SELECTION (BBKNN)
# =============================================================================
# If your BBKNN used key_added (e.g. "neighbors_bbknn"), set it here.
NEIGHBORS_KEY = "neighbors"

# =============================================================================
# REPRODUCIBILITY
# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
try:
    from anndata import settings as anndata_settings
    anndata_settings.seed = RANDOM_SEED
except Exception:
    pass

print(f"[INFO] Random seed set to {RANDOM_SEED}")

# =============================================================================
# INPUT/OUTPUT
# =============================================================================
INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_celltypist_filtered_final.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/epithelial_subcluster_bbknn_restrict_to")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# ANNOTATION KEYS
# =============================================================================
BATCH_KEY = "dataset"             # visualization only; set None if not present
CELLTYPE_COLUMN = "celltypist_pred"

# =============================================================================
# SUBCLUSTERING (INTENTIONALLY CONSERVATIVE)
# =============================================================================
MIN_CELLS_FOR_SUBCLUSTER = 100

RESOLUTION_FIXED = 0.2
def calculate_resolution(n_cells, fixed_res=RESOLUTION_FIXED):
    """
    Conservative design: fixed low resolution to avoid over-splitting.
    This pipeline intentionally uses a single fixed resolution across cell types.
    """
    return float(fixed_res)

# =============================================================================
# MARKER GENE SETTINGS
# =============================================================================
MARKER_METHOD = "wilcoxon"
MARKER_LOGFC_THRESHOLD = 0.5
MARKER_PADJ_THRESHOLD = 0.05

# pct filters (recommended for high drop-out epithelial data)
MARKER_MIN_PCT = 0.25
MARKER_MIN_DELTA_PCT = 0.10

TOP_N_MARKERS = 20

# =============================================================================
# MARKER GENE FILTERING
# =============================================================================
FILTER_MT_GENES = True
FILTER_RIBO_GENES = True
FILTER_UNANNOTATED = True
FILTER_HISTONE_GENES = True

# Reduced stress filtering: dissociation/IEG artifacts only
FILTER_STRESS_GENES = True
STRESS_MATCH_THRESHOLD = 0.5

STRESS_SIGNATURE_GENES = [
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "EGR1", "EGR2", "EGR3", "EGR4",
    "ATF3",
    "DUSP1", "DUSP2", "DUSP4", "DUSP5",
    "HSPA1A", "HSPA1B",
    "IER2", "IER3", "IER5",
    "NR4A1", "NR4A2", "NR4A3",
    "ZFP36", "ZFP36L1", "ZFP36L2"
]

# =============================================================================
# GENE UNIVERSE / DE INPUT SELECTION
# =============================================================================
# Lock gene-universe sources to avoid raw/var mismatch.
USE_RAW_FOR_GENE_FILTERING = True   # creates raw-gene flags if raw exists
ALSO_BUILD_VAR_FILTERS = True       # always build var-name filter dict too (recommended)

# DE input selection: "auto" (recommended), or force:
#   DE_MODE = "raw"    -> use_raw=True (assumes raw is log1p)
#   DE_MODE = "X"      -> use_raw=False (uses adata.X)
#   DE_MODE = "layer"  -> use_raw=False + layer=DE_LAYER
DE_MODE = "auto"
DE_LAYER = None  # if DE_MODE == "layer", set e.g. "log1p" / "lognorm"

# =============================================================================
# VISUALIZATION
# =============================================================================
FIGURE_DPI = 300
FIGURE_FORMAT = "pdf"
UMAP_SIZE = 5
UMAP_ALPHA = 0.6

# =============================================================================
# PERFORMANCE
# =============================================================================
N_JOBS = 48
MEMORY_LIMIT_GB = 256

# =============================================================================
# MARKER DE SETTINGS (performance/robustness)
# =============================================================================
TOP_DE_GENES = 4000  # limit DE retrieval to reduce runtime/uns size (safe for top markers)

# =============================================================================
# REGEX / SETS (avoid repeated compilation)
# =============================================================================
STRESS_SET = set(STRESS_SIGNATURE_GENES)
UNANNOTATED_REGEX = r"(LINC|RP11-|^RP\d+-|^AC\d+|^AL\d+|CTD-|CTB-|pseudogene)"

sc.settings.n_jobs = N_JOBS
sc.settings.figdir = OUTPUT_DIR / "figures"
sc.settings.figdir.mkdir(exist_ok=True)

print("=" * 80)
print(f"Epithelial Subcluster Analysis Pipeline {PIPELINE_VERSION}")
print("=" * 80)
print(f"Input:      {INPUT_H5AD}")
print(f"Output:     {OUTPUT_DIR}")
print(f"Cell type:  {CELLTYPE_COLUMN}")
print(f"Batch key:  {BATCH_KEY if BATCH_KEY else 'None'}")
print(f"Resolution: FIXED {RESOLUTION_FIXED}")
print(f"DE mode:    {DE_MODE} (layer={DE_LAYER})")
print(f"Date:       {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# %% [markdown]
# ## Step 1: Load Data

# %%
print("\n" + "=" * 80)
print("STEP 1: Loading Data")
print("=" * 80)

t0 = time.time()
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[OK] Loaded in {time.time() - t0:.1f}s")

print("\nDataset Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print("\nStructure:")
print(f"  .X shape: {adata.X.shape} ({type(adata.X)})")
print(f"  .layers: {list(adata.layers.keys())}")
print(f"  .raw: {'present' if adata.raw is not None else 'missing'}")
if adata.raw is not None:
    print(f"  .raw.X shape: {adata.raw.X.shape} ({type(adata.raw.X)})")

print("\nEmbeddings:")
print(f"  obsm keys: {list(adata.obsm.keys())}")
if "X_umap" not in adata.obsm:
    if "X_umap_scvi" in adata.obsm:
        print("[WARN] X_umap not found; using X_umap_scvi as X_umap")
        adata.obsm["X_umap"] = adata.obsm["X_umap_scvi"]
    else:
        print("[WARN] X_umap not found; will compute UMAP later if needed")

print("\nAnnotation Check:")
if CELLTYPE_COLUMN not in adata.obs.columns:
    raise ValueError(f"Cell type column '{CELLTYPE_COLUMN}' not found in adata.obs")

celltype_counts = adata.obs[CELLTYPE_COLUMN].value_counts()
print(f"  Cell types: {len(celltype_counts)}")
print("  Distribution:")
for ct, n in celltype_counts.items():
    print(f"    {ct}: {n:,}")

if BATCH_KEY and BATCH_KEY in adata.obs.columns:
    print(f"\n  Batch key '{BATCH_KEY}': {adata.obs[BATCH_KEY].nunique()} batches")
else:
    print(f"\n  Batch key '{BATCH_KEY}' not found; disabling batch visualization")
    BATCH_KEY = None

# %% [markdown]
# ## Step 2: Build Gene Filters (raw + var, mismatch-safe)

# %% [markdown]
# ## Step 2: Build Gene Filters (raw + var, mismatch-safe)

# %%
print("\n" + "=" * 80)
print("STEP 2: Building Gene Filters (Mismatch-safe)")
print("=" * 80)

# -------------------------------------------------------------------------
# Local safe constants (in case you didn't apply Patch 0)
# -------------------------------------------------------------------------
try:
    STRESS_SET
except NameError:
    STRESS_SET = set(STRESS_SIGNATURE_GENES)

try:
    UNANNOTATED_REGEX
except NameError:
    UNANNOTATED_REGEX = r"(LINC|RP11-|^RP\d+-|^AC\d+|^AL\d+|CTD-|CTB-|pseudogene)"

def _ensure_index_str(idx):
    # idx can be pandas Index; ensure safe string operations
    return pd.Index(idx.astype(str))

def _build_filter_flags(gene_names: pd.Index):
    """
    Build mismatch-safe filter flags AND write-consistent category masks.

    Returns:
      to_filter (np.ndarray bool),
      summary (dict),
      flags (dict[str, np.ndarray bool])  # consistent with actual filtering rules
    """
    gene_names = _ensure_index_str(gene_names)

    # default masks
    is_mt = np.zeros(len(gene_names), dtype=bool)
    is_ribo = np.zeros(len(gene_names), dtype=bool)
    is_ensg = np.zeros(len(gene_names), dtype=bool)
    is_unanno = np.zeros(len(gene_names), dtype=bool)
    is_histone = np.zeros(len(gene_names), dtype=bool)
    is_stress = np.zeros(len(gene_names), dtype=bool)

    if FILTER_MT_GENES:
        is_mt = gene_names.str.startswith(("MT-", "mt-"))

    if FILTER_RIBO_GENES:
        is_ribo = gene_names.str.startswith(("RPL", "RPS", "MRPL", "MRPS"))

    if FILTER_UNANNOTATED:
        is_ensg = gene_names.str.startswith("ENSG")
        is_unanno = (
            gene_names.str.contains("LINC", case=False, na=False)
            | gene_names.str.contains(r"^RP\d+-", regex=True, na=False)
            | gene_names.str.contains("RP11-", na=False)
            | gene_names.str.contains(r"^AC\d+", regex=True, na=False)
            | gene_names.str.contains(r"^AL\d+", regex=True, na=False)
            | gene_names.str.contains("CTD-", na=False)
            | gene_names.str.contains("CTB-", na=False)
            | gene_names.str.contains("pseudogene", case=False, na=False)
        )

    if FILTER_HISTONE_GENES:
        is_histone = (
            gene_names.str.startswith(("HIST1", "HIST2", "HIST3", "HIST4"))
            | gene_names.str.startswith("H1-")
            | gene_names.str.startswith("H2A")
            | gene_names.str.startswith("H2B")
            | gene_names.str.startswith(("H3-", "H3C"))
            | gene_names.str.startswith(("H4-", "H4C"))
        )

    if FILTER_STRESS_GENES:
        is_stress = gene_names.isin(STRESS_SET)

    to_filter = is_mt | is_ribo | is_ensg | is_unanno | is_histone | is_stress

    summary = {
        "n_genes": len(gene_names),
        "n_filter": int(to_filter.sum()),
        "mt": int(is_mt.sum()),
        "ribo": int(is_ribo.sum()),
        "ensg": int(is_ensg.sum()),
        "unannotated": int(is_unanno.sum()),
        "histone": int(is_histone.sum()),
        "stress_ieg": int(is_stress.sum()),
    }

    flags = {
        "is_mt": is_mt,
        "is_ribo": is_ribo,
        "is_ensg": is_ensg,
        "is_unannotated": is_unanno,
        "is_histone": is_histone,
        "is_stress_ieg": is_stress,
        "filter_from_markers": to_filter,
    }
    return to_filter, summary, flags

# -------------------------------------------------------------------------
# Build RAW filters (if raw exists and enabled)
# -------------------------------------------------------------------------
gene_filter_dict_raw = None
raw_filter_flags = None
raw_summary = None
raw_gene_names = None

if USE_RAW_FOR_GENE_FILTERING and adata.raw is not None:
    raw_gene_names = _ensure_index_str(adata.raw.var_names)
    raw_filter_flags, raw_summary, raw_flags = _build_filter_flags(raw_gene_names)

    # write flags to raw.var safely (CONSISTENT with actual filter logic)
    rv = adata.raw.var
    for k, v in raw_flags.items():
        rv[k] = v

    gene_filter_dict_raw = dict(zip(raw_gene_names, raw_filter_flags))

    print("\n[RAW] Gene filter summary:")
    for k, v in raw_summary.items():
        print(f"  {k}: {v}")
    if raw_summary["stress_ieg"] < len(STRESS_SIGNATURE_GENES) * STRESS_MATCH_THRESHOLD:
        missing_pct = (1 - raw_summary["stress_ieg"] / len(STRESS_SIGNATURE_GENES)) * 100
        print(f"  [WARN] {missing_pct:.1f}% of IEG list not found in .raw.var_names")

# -------------------------------------------------------------------------
# Build VAR filters (recommended always; used for filtered_genes_list.csv export)
# -------------------------------------------------------------------------
gene_filter_dict_var = None
var_filter_flags = None
var_summary = None

var_gene_names = _ensure_index_str(adata.var_names)
var_filter_flags, var_summary, var_flags = _build_filter_flags(var_gene_names)

# write flags to adata.var (CONSISTENT with actual filter logic)
for k, v in var_flags.items():
    adata.var[k] = v

gene_filter_dict_var = dict(zip(var_gene_names, var_filter_flags))

print("\n[VAR] Gene filter summary:")
for k, v in var_summary.items():
    print(f"  {k}: {v}")
if var_summary["stress_ieg"] < len(STRESS_SIGNATURE_GENES) * STRESS_MATCH_THRESHOLD:
    missing_pct = (1 - var_summary["stress_ieg"] / len(STRESS_SIGNATURE_GENES)) * 100
    print(f"  [WARN] {missing_pct:.1f}% of IEG list not found in adata.var_names")

# -------------------------------------------------------------------------
# Export filtered gene list (based on VAR by default, because it always exists)
# -------------------------------------------------------------------------
def _gene_category(g: str) -> str:
    g = str(g)
    if g.startswith(("MT-", "mt-")):
        return "MT"
    if g.startswith(("RPL", "RPS", "MRPL", "MRPS")):
        return "Ribo"
    if g.startswith("ENSG"):
        return "ENSG"
    if g.startswith(("HIST1", "HIST2", "HIST3", "HIST4", "H1-", "H2A", "H2B", "H3-", "H3C", "H4-", "H4C")):
        return "Histone"
    if g in STRESS_SET:
        return "IEG_Dissociation"
    if re.search(UNANNOTATED_REGEX, g, flags=re.IGNORECASE):
        return "Unannotated"
    return "Other"

filtered_genes = var_gene_names[var_filter_flags]
filtered_genes_df = pd.DataFrame({
    "gene": filtered_genes.values,
    "category": [_gene_category(g) for g in filtered_genes.values],
    "reason": "Technical/state artifact excluded from identity marker reporting",
    "source": "adata.var_names",  # explicit to avoid raw/var confusion
})

filtered_genes_csv = OUTPUT_DIR / "filtered_genes_list.csv"
filtered_genes_df.to_csv(filtered_genes_csv, index=False)
print(f"\n[OK] Saved filtered gene list: {filtered_genes_csv.name} ({len(filtered_genes_df):,} genes)")


# %% [markdown]
# ## Step 3: Preserve Original UMAP

# %%
print("\n" + "=" * 80)
print("STEP 3: Preserving Original UMAP")
print("=" * 80)

if "X_umap" in adata.obsm:
    adata.obsm["X_umap_original"] = adata.obsm["X_umap"].copy()
    print("[OK] Saved X_umap_original")
else:
    print("[WARN] No X_umap found; will compute later if needed")

for key in ["X_pca", "X_scvi", "X_scanvi"]:
    if key in adata.obsm:
        adata.obsm[f"{key}_original"] = adata.obsm[key].copy()
        print(f"[OK] Preserved {key} -> {key}_original")

# %% [markdown]
# ## Step 4: Subcluster Each Cell Type (restrict_to on global BBKNN graph)

# %%
print("\n" + "=" * 80)
print("STEP 4: Subclustering (restrict_to on global graph)")
print("=" * 80)
def resolve_neighbors_graph(adata, neighbors_key="neighbors"):
    if neighbors_key not in adata.uns:
        raise ValueError(
            f"neighbors_key='{neighbors_key}' not found in adata.uns. "
            f"Available keys: {list(adata.uns.keys())}"
        )
    nuns = adata.uns[neighbors_key]
    conn_key = nuns.get("connectivities_key", "connectivities")
    dist_key = nuns.get("distances_key", "distances")

    if conn_key not in adata.obsp:
        raise ValueError(
            f"connectivities_key='{conn_key}' not found in adata.obsp. "
            f"Available keys: {list(adata.obsp.keys())}"
        )
    if dist_key not in adata.obsp:
        print(f"[WARN] distances_key='{dist_key}' not found in adata.obsp (OK if not needed).")

    params = nuns.get("params", {})
    return conn_key, dist_key, params

CONN_KEY, DIST_KEY, params = resolve_neighbors_graph(adata, NEIGHBORS_KEY)

print("\nGlobal neighbors graph (resolved):")
print(f"s  neighbors_key:       {NEIGHBORS_KEY}")
print(f"  connectivities_key:   {CONN_KEY}")
print(f"  distances_key:        {DIST_KEY}")
print(f"  method:               {params.get('method', 'unknown')}")
print(f"  n_neighbors:          {params.get('n_neighbors', 'unknown')}")
print(f"  n_pcs:                {params.get('n_pcs', 'unknown')}")

celltype_counts = adata.obs[CELLTYPE_COLUMN].astype(str).value_counts()
celltypes = celltype_counts.index.tolist()  # stable order (desc by size)
print(f"\nCell types to process: {len(celltypes)}")

# Initialize outputs
adata.obs["subcluster"] = "unassigned"
adata.obs["subcluster_leiden"] = -1  # int id for stable sorting

processed_celltypes = []

for celltype in celltypes:
    print("\n" + "-" * 80)
    print(f"Processing celltype: {celltype}")

    mask = (adata.obs[CELLTYPE_COLUMN].astype(str) == str(celltype))
    n_cells = int(mask.sum())
    print(f"Cells: {n_cells:,}")

    if n_cells < MIN_CELLS_FOR_SUBCLUSTER:
        print(f"[SKIP] < {MIN_CELLS_FOR_SUBCLUSTER} cells; assign single cluster C0")
        adata.obs.loc[mask, "subcluster"] = f"{celltype}_C0"
        adata.obs.loc[mask, "subcluster_leiden"] = 0
        continue

    res = calculate_resolution(n_cells)
    print(f"Resolution: {res:.3f} (fixed)")

    key_temp = f"leiden_temp_{re.sub(r'[^A-Za-z0-9_]+', '_', str(celltype))}"
    if key_temp in adata.obs.columns:
        adata.obs.drop(columns=[key_temp], inplace=True)

    try:
                # run leiden on GLOBAL graph restricted to this celltype
        leiden_kwargs = dict(
            resolution=res,
            restrict_to=(CELLTYPE_COLUMN, [celltype]),
            key_added=key_temp,
            random_state=RANDOM_SEED,
        )

        # prefer explicit neighbors_key if supported by scanpy version
        try:
            sc.tl.leiden(adata, neighbors_key=NEIGHBORS_KEY, **leiden_kwargs)
        except TypeError:
            # older scanpy: no neighbors_key param; relies on default adata.uns['neighbors']
            if NEIGHBORS_KEY != "neighbors":
                print(f"[WARN] scanpy.leiden has no neighbors_key; requires NEIGHBORS_KEY='neighbors' to be default.")
            sc.tl.leiden(adata, **leiden_kwargs)

        # Extract this celltype's restricted Leiden labels
        leiden_vals = adata.obs.loc[mask, key_temp].astype(str).values

        # parse "celltype,cluster_id" -> cluster_id (use rsplit for robustness)
        sub_ids = []
        for v in leiden_vals:
            if "," in v:
                try:
                    _, cid = v.rsplit(",", 1)
                    sub_ids.append(int(cid))
                except Exception:
                    sub_ids.append(0)
            elif v.lower() == "nan":
                sub_ids.append(0)
            else:
                try:
                    sub_ids.append(int(v))
                except Exception:
                    sub_ids.append(0)

        sub_ids = np.array(sub_ids, dtype=int)
        adata.obs.loc[mask, "subcluster_leiden"] = sub_ids
        adata.obs.loc[mask, "subcluster"] = [f"{celltype}_C{c}" for c in sub_ids]

        n_clusters = int(pd.Series(sub_ids).nunique())
        print(f"[OK] Found {n_clusters} subclusters")

        # small cluster warning
        vc = pd.Series(sub_ids).value_counts()
        tiny = vc[vc < 20]
        if len(tiny) > 0:
            print(f"[WARN] {len(tiny)} subclusters have <20 cells:")
            for cid, cnt in tiny.items():
                print(f"  C{cid}: {cnt} cells")

        processed_celltypes.append(celltype)
    except Exception as e:
        print(f"[FAIL] Leiden failed for {celltype}: {e}")
        print("[FALLBACK] Assign single cluster C0")
        adata.obs.loc[mask, "subcluster"] = f"{celltype}_C0"
        adata.obs.loc[mask, "subcluster_leiden"] = 0

    # cleanup temp column
    # cleanup temp column + uns to keep object clean
    if key_temp in adata.obs.columns:
        adata.obs.drop(columns=[key_temp], inplace=True)
    if key_temp in adata.uns:
        del adata.uns[key_temp]


    gc.collect()

print("\n" + "=" * 80)
print("Subclustering Summary")
print("=" * 80)
print(f"Processed cell types: {len(processed_celltypes)}")
print(f"Total cells: {adata.n_obs:,}")
print(f"Total subclusters: {adata.obs['subcluster'].nunique()}")

# %% [markdown]
# ## Step 5: Visualization

# %%
print("\n" + "=" * 80)
print("STEP 5: Visualization")
print("=" * 80)

# restore original UMAP if present
if "X_umap_original" in adata.obsm:
    adata.obsm["X_umap"] = adata.obsm["X_umap_original"].copy()
    print("[OK] Using original UMAP")
elif "X_umap" not in adata.obsm:
    print("[INFO] Computing UMAP (no original available)")
    sc.tl.umap(adata)

# 5.1 All subclusters
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(
    adata,
    color="subcluster",
    title="Epithelial Subclusters (Post-BBKNN, restrict_to)",
    frameon=False,
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    legend_loc="right margin",
    show=False,
    ax=ax,
)
plt.tight_layout()
plt.savefig(sc.settings.figdir / f"umap_subclusters_all.{FIGURE_FORMAT}", dpi=FIGURE_DPI, bbox_inches="tight")
plt.close()
print(f"[OK] Saved umap_subclusters_all.{FIGURE_FORMAT}")

# 5.2 Comparison
fig, axes = plt.subplots(1, 2, figsize=(24, 10))
sc.pl.umap(
    adata,
    color=CELLTYPE_COLUMN,
    title="Original Cell Types",
    frameon=False,
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    show=False,
    ax=axes[0],
)
sc.pl.umap(
    adata,
    color="subcluster",
    title="Subclusters",
    frameon=False,
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    legend_loc="right margin",
    show=False,
    ax=axes[1],
)
plt.tight_layout()
plt.savefig(sc.settings.figdir / f"umap_comparison_celltype_vs_subcluster.{FIGURE_FORMAT}", dpi=FIGURE_DPI, bbox_inches="tight")
plt.close()
print(f"[OK] Saved umap_comparison_celltype_vs_subcluster.{FIGURE_FORMAT}")

# 5.3 Batch distribution (optional)
if BATCH_KEY and BATCH_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(
        adata,
        color=BATCH_KEY,
        title="Batch Distribution (BBKNN already applied)",
        frameon=False,
        size=UMAP_SIZE,
        alpha=UMAP_ALPHA,
        show=False,
        ax=ax,
    )
    plt.tight_layout()
    plt.savefig(sc.settings.figdir / f"umap_batch_distribution.{FIGURE_FORMAT}", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved umap_batch_distribution.{FIGURE_FORMAT}")

# 5.4 Per-celltype highlight plots (in-place, memory-safe)
for celltype in processed_celltypes:
    mask = (adata.obs[CELLTYPE_COLUMN].astype(str) == str(celltype))
    adata.obs["_temp_highlight"] = "Other"
    adata.obs.loc[mask, "_temp_highlight"] = adata.obs.loc[mask, "subcluster"].astype(str)

    highlighted_groups = adata.obs.loc[mask, "subcluster"].astype(str).unique().tolist()
    all_groups = highlighted_groups + ["Other"]

    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(
        adata,
        color="_temp_highlight",
        title=f"{celltype} Subclusters",
        frameon=False,
        size=UMAP_SIZE,
        alpha=UMAP_ALPHA,
        groups=all_groups,
        show=False,
        ax=ax,
    )
    plt.tight_layout()
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(celltype))
    plt.savefig(sc.settings.figdir / f"umap_subcluster_{safe}.{FIGURE_FORMAT}", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close()

    adata.obs.drop(columns=["_temp_highlight"], inplace=True)

print("[OK] Per-celltype plots completed")

# %% [markdown]
# ## Step 6: Marker Gene Analysis (stable result retrieval + pct filters)

# %%
print("\n" + "=" * 80)
print("STEP 6: Marker Gene Analysis")
print("=" * 80)

def _sample_matrix_values(X, n=20000, seed=42):
    rng = np.random.default_rng(seed)
    if sparse.issparse(X):
        data = X.data
        if data.size == 0:
            return np.array([], dtype=float)
        if data.size <= n:
            return np.array(data, dtype=float)
        idx = rng.choice(data.size, size=n, replace=False)
        return np.array(data[idx], dtype=float)
    else:
        arr = np.asarray(X)
        flat = arr.ravel()
        if flat.size == 0:
            return np.array([], dtype=float)
        if flat.size <= n:
            return flat.astype(float)
        idx = rng.choice(flat.size, size=n, replace=False)
        return flat[idx].astype(float)

def _looks_like_counts(X):
    vals = _sample_matrix_values(X)
    if vals.size == 0:
        return False
    vmax = float(np.nanmax(vals))
    # many exact integers + large max suggests counts
    frac_int = float(np.mean(np.isclose(vals, np.round(vals))))
    return (vmax > 30) or (frac_int > 0.90 and vmax > 15)

def _choose_de_settings(adata):
    """
    Return (use_raw_bool, layer_or_None, filter_dict_to_use, de_note)
    """
    # forced modes
    if DE_MODE == "raw":
        if adata.raw is None:
            return (False, None, gene_filter_dict_var, "DE_MODE=raw but .raw missing -> fallback to X")
        return (True, None, gene_filter_dict_raw if gene_filter_dict_raw is not None else gene_filter_dict_var, "DE_MODE=raw")
    if DE_MODE == "X":
        return (False, None, gene_filter_dict_var, "DE_MODE=X")
    if DE_MODE == "layer":
        if DE_LAYER is None:
            return (False, None, gene_filter_dict_var, "DE_MODE=layer but DE_LAYER=None -> fallback to X")
        if DE_LAYER not in adata.layers:
            return (False, None, gene_filter_dict_var, f"DE_MODE=layer but layer '{DE_LAYER}' missing -> fallback to X")
        return (False, DE_LAYER, gene_filter_dict_var, f"DE_MODE=layer ({DE_LAYER})")

    # auto mode
    raw_counts_like = False
    x_counts_like = False
    if adata.raw is not None:
        raw_counts_like = _looks_like_counts(adata.raw.X)
    x_counts_like = _looks_like_counts(adata.X)

    # If raw looks like counts, avoid use_raw unless user forced
    if adata.raw is not None and (not raw_counts_like):
        # raw seems log1p-ish; prefer raw (consistent with many pipelines)
        return (True, None, gene_filter_dict_raw if gene_filter_dict_raw is not None else gene_filter_dict_var, "AUTO: use_raw=True (raw appears log1p)")
    else:
        # prefer a log1p layer if available
        for cand in ["log1p", "lognorm", "X_log1p", "normalized", "lognormalized"]:
            if cand in adata.layers:
                return (False, cand, gene_filter_dict_var, f"AUTO: use layer='{cand}' (raw counts-like or missing)")
        # fallback to X
        note = "AUTO: use_raw=False (raw counts-like or missing); using X"
        if x_counts_like:
            note += " [WARN: X also looks counts-like; DE may be distorted unless X is log1p]"
        return (False, None, gene_filter_dict_var, note)

DE_USE_RAW, DE_USE_LAYER, marker_filter_dict, de_note = _choose_de_settings(adata)
print(f"[INFO] DE selection: use_raw={DE_USE_RAW}, layer={DE_USE_LAYER} | {de_note}")

def filter_marker_results(marker_df, filter_dict):
    marker_df = marker_df.copy()
    marker_df["to_filter"] = marker_df["names"].map(filter_dict).fillna(False)
    n_total = len(marker_df)
    n_filtered = int(marker_df["to_filter"].sum())
    kept = marker_df.loc[~marker_df["to_filter"]].drop(columns=["to_filter"]).copy()
    if n_filtered > 0:
        print(f"    Filtered {n_filtered}/{n_total} markers ({n_filtered/n_total*100:.1f}%)")
    return kept

all_markers = {}
all_markers_filtered = {}

for celltype in processed_celltypes:
    print("\n" + "-" * 80)
    print(f"Marker analysis: {celltype}")

    mask = (adata.obs[CELLTYPE_COLUMN].astype(str) == str(celltype))
    n_cells = int(mask.sum())
    n_sub = int(adata.obs.loc[mask, "subcluster_leiden"].nunique())
    print(f"  Cells: {n_cells:,} | Subclusters: {n_sub}")

    if n_sub <= 1:
        print("  [SKIP] Only 1 subcluster")
        continue

    # Important: keep a single slice object reference for result retrieval
    adata_ct = adata[mask].copy()

    key_de = "rank_genes_groups_temp"

    try:
        # run DE
        rg_kwargs = dict(
            groupby="subcluster_leiden",
            method=MARKER_METHOD,
            pts=True,
            key_added=key_de,
            n_genes=TOP_DE_GENES,
        )
        if DE_USE_RAW:
            rg_kwargs["use_raw"] = True
        else:
            rg_kwargs["use_raw"] = False
            if DE_USE_LAYER is not None:
                rg_kwargs["layer"] = DE_USE_LAYER

        sc.tl.rank_genes_groups(adata_ct, **rg_kwargs)

        result = adata_ct.uns[key_de]
        groups = result["names"].dtype.names

        markers_list = []
        for g in groups:
            df = pd.DataFrame({
                "celltype": celltype,
                "subcluster": [f"{celltype}_C{g}"] * len(result["names"][g]),
                "leiden_cluster": [str(g)] * len(result["names"][g]),
                "names": result["names"][g],
                "scores": result["scores"][g],
                "pvals": result["pvals"][g],
                "pvals_adj": result["pvals_adj"][g],
                "logfoldchanges": result["logfoldchanges"][g],
            })
            if "pts" in result:
                df["pct_in_group"] = result["pts"][g]
            if "pts_rest" in result:
                df["pct_in_others"] = result["pts_rest"][g]
            markers_list.append(df)

        markers_df = pd.concat(markers_list, ignore_index=True)

        # quality filters (statistics + pct)
        cond = (
            (markers_df["pvals_adj"] <= MARKER_PADJ_THRESHOLD)
            & (markers_df["logfoldchanges"] >= MARKER_LOGFC_THRESHOLD)
        )
        if "pct_in_group" in markers_df.columns:
            cond = cond & (markers_df["pct_in_group"] >= MARKER_MIN_PCT)
        if ("pct_in_group" in markers_df.columns) and ("pct_in_others" in markers_df.columns):
            cond = cond & ((markers_df["pct_in_group"] - markers_df["pct_in_others"]) >= MARKER_MIN_DELTA_PCT)

        markers_df = markers_df.loc[cond].copy()
        print(f"  [OK] Significant markers (pre-filter): {len(markers_df):,}")

        all_markers[celltype] = markers_df.copy()

        print("  Filtering technical/IEG markers...")
        markers_df_f = filter_marker_results(markers_df, marker_filter_dict)
        print(f"  [OK] Clean markers (post-filter): {len(markers_df_f):,}")

        all_markers_filtered[celltype] = markers_df_f

    except Exception as e:
        print(f"  [FAIL] Marker analysis failed: {e}")

    # cleanup: remove temp key from the slice uns to keep objects small
    try:
        if key_de in adata_ct.uns:
            del adata_ct.uns[key_de]
    except Exception:
        pass

    gc.collect()

print("\n" + "=" * 80)
print(f"Marker analysis completed for {len(all_markers)} cell types")
print("=" * 80)

# %% [markdown]
# ## Step 7: Marker Visualization (Dotplot/Heatmap)

# %%
print("\n" + "=" * 80)
print("STEP 7: Marker Visualization")
print("=" * 80)

markers_for_viz = all_markers_filtered

for celltype, markers_df in markers_for_viz.items():
    print(f"\nVisualizing: {celltype}")

    mask = (adata.obs[CELLTYPE_COLUMN].astype(str) == str(celltype))
    adata_temp = adata[mask].copy()

    # Ensure numeric ordering for groupby
    adata_temp.obs["subcluster_leiden"] = pd.to_numeric(adata_temp.obs["subcluster_leiden"], errors="coerce").fillna(0).astype(int)
    adata_temp.obs["subcluster_leiden"] = pd.Categorical(
        adata_temp.obs["subcluster_leiden"].astype(str),
        categories=[str(i) for i in sorted(adata_temp.obs["subcluster_leiden"].unique())],
        ordered=True,
    )

    # select top markers per subcluster (by smallest padj)
    top_markers = []
    for sub in markers_df["subcluster"].unique():
        sub_df = markers_df.loc[markers_df["subcluster"] == sub]
        top_n = sub_df.nsmallest(10, "pvals_adj")["names"].astype(str).tolist()
        top_markers.extend(top_n)

    top_markers = list(dict.fromkeys(top_markers))[:50]

    # keep only genes available in the same namespace used for plots
    # align plotting use_raw with DE_USE_RAW, to avoid gene mismatch
    if DE_USE_RAW and adata_temp.raw is not None:
        top_markers = [g for g in top_markers if g in adata_temp.raw.var_names]
        plot_use_raw = True
    else:
        top_markers = [g for g in top_markers if g in adata_temp.var_names]
        plot_use_raw = False

    if len(top_markers) == 0:
        print("  [SKIP] No valid markers after gene availability check")
        del adata_temp
        gc.collect()
        continue

    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(celltype))

    # dotplot
    try:
        dp = sc.pl.dotplot(
            adata_temp,
            var_names=top_markers,
            groupby="subcluster_leiden",
            standard_scale="var",
            use_raw=plot_use_raw,
            show=False,
            return_fig=True,
        )
        dp.savefig(sc.settings.figdir / f"dotplot_markers_{safe}_filtered.{FIGURE_FORMAT}")
        plt.close("all")
        print(f"  [OK] Saved dotplot_markers_{safe}_filtered.{FIGURE_FORMAT}")
    except Exception as e:
        print(f"  [WARN] Dotplot failed: {e}")

    # heatmap (top 20)
    try:
        hm_genes = top_markers[:20]
        if len(hm_genes) > 0:
            hm = sc.pl.heatmap(
                adata_temp,
                var_names=hm_genes,
                groupby="subcluster_leiden",
                standard_scale="var",
                use_raw=plot_use_raw,
                show=False,
                return_fig=True,
            )
            hm.savefig(sc.settings.figdir / f"heatmap_markers_{safe}_filtered.{FIGURE_FORMAT}")
            plt.close("all")
            print(f"  [OK] Saved heatmap_markers_{safe}_filtered.{FIGURE_FORMAT}")
    except Exception as e:
        print(f"  [WARN] Heatmap failed: {e}")

    del adata_temp
    gc.collect()

# Example markers on UMAP (first celltype)
if len(markers_for_viz) > 0:
    first_ct = list(markers_for_viz.keys())[0]
    ex = markers_for_viz[first_ct].nsmallest(6, "pvals_adj")["names"].astype(str).tolist()
    if (DE_USE_RAW and adata.raw is not None):
        ex = [g for g in ex if g in adata.raw.var_names][:6]
        plot_use_raw = True
    else:
        ex = [g for g in ex if g in adata.var_names][:6]
        plot_use_raw = False

    if len(ex) > 0:
        sc.pl.umap(
            adata,
            color=ex,
            use_raw=plot_use_raw,
            ncols=3,
            frameon=False,
            size=UMAP_SIZE * 0.8,
            alpha=UMAP_ALPHA,
            show=False,
        )
        plt.suptitle("Example Marker Genes (Filtered)", y=1.02)
        plt.tight_layout()
        plt.savefig(sc.settings.figdir / f"umap_example_markers_filtered.{FIGURE_FORMAT}", dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close()
        print(f"\n[OK] Saved umap_example_markers_filtered.{FIGURE_FORMAT}")

# %% [markdown]
# ## Step 8: Save Results

# %%
print("\n" + "=" * 80)
print("STEP 8: Saving Results")
print("=" * 80)

output_h5ad = OUTPUT_DIR / "epithelial_with_subclusters.h5ad"
print(f"[INFO] Writing h5ad: {output_h5ad}")
adata.write_h5ad(output_h5ad, compression="gzip", compression_opts=9)
file_size_gb = output_h5ad.stat().st_size / 1e9
print(f"[OK] Saved: {output_h5ad.name} ({file_size_gb:.2f} GB)")

# Subcluster annotations
subcluster_df = adata.obs[[CELLTYPE_COLUMN, "subcluster", "subcluster_leiden"]].copy()
if BATCH_KEY and BATCH_KEY in adata.obs.columns:
    subcluster_df[BATCH_KEY] = adata.obs[BATCH_KEY]
subcluster_csv = OUTPUT_DIR / "subcluster_annotations.csv"
subcluster_df.to_csv(subcluster_csv, index=False)
print(f"[OK] Saved: {subcluster_csv.name}")

# Markers
all_markers_df = None
all_markers_filtered_df = None

if len(all_markers) > 0:
    all_markers_df = pd.concat(all_markers.values(), ignore_index=True)
    markers_all_csv = OUTPUT_DIR / "markers_all_subclusters_unfiltered.csv"
    all_markers_df.to_csv(markers_all_csv, index=False)
    print(f"[OK] Saved: {markers_all_csv.name} ({len(all_markers_df):,} rows)")

if len(all_markers_filtered) > 0:
    all_markers_filtered_df = pd.concat(all_markers_filtered.values(), ignore_index=True)
    markers_filtered_csv = OUTPUT_DIR / "markers_all_subclusters_filtered.csv"
    all_markers_filtered_df.to_csv(markers_filtered_csv, index=False)
    print(f"[OK] Saved: {markers_filtered_csv.name} ({len(all_markers_filtered_df):,} rows)")

    # Top N per subcluster
    top_list = []
    for sub in all_markers_filtered_df["subcluster"].unique():
        sub_df = all_markers_filtered_df.loc[all_markers_filtered_df["subcluster"] == sub]
        top_list.append(sub_df.nsmallest(TOP_N_MARKERS, "pvals_adj"))
    top_df = pd.concat(top_list, ignore_index=True) if len(top_list) > 0 else None

    if top_df is not None:
        top_csv = OUTPUT_DIR / f"markers_top{TOP_N_MARKERS}_per_subcluster_filtered.csv"
        top_df.to_csv(top_csv, index=False)
        print(f"[OK] Saved: {top_csv.name} ({len(top_df):,} rows)")

    # per-celltype
    for ct, df in all_markers_filtered.items():
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(ct))
        df.to_csv(OUTPUT_DIR / f"markers_{safe}_filtered.csv", index=False)
    print("[OK] Saved per-celltype marker files")

# Summary report
summary_file = OUTPUT_DIR / "subcluster_summary.txt"
with open(summary_file, "w") as f:
    f.write("=" * 80 + "\n")
    f.write(f"Epithelial Subcluster Analysis Summary ({PIPELINE_VERSION})\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input: {INPUT_H5AD}\n")
    f.write(f"Cell type column: {CELLTYPE_COLUMN}\n")
    f.write(f"Batch key: {BATCH_KEY if BATCH_KEY else 'None'}\n")
    f.write(f"Random seed: {RANDOM_SEED}\n\n")

    f.write("Subclustering:\n")
    f.write(f"  neighbors_key: {NEIGHBORS_KEY}\n")
    f.write(f"  connectivities_key: {CONN_KEY}\n")
    f.write(f"  marker_filter_dict source: {'raw' if (DE_USE_RAW and gene_filter_dict_raw is not None) else 'var'}\n")

    f.write(f"  Resolution: FIXED {RESOLUTION_FIXED}\n")
    f.write(f"  Total subclusters: {adata.obs['subcluster'].nunique()}\n")
    f.write(f"  Processed cell types: {len(processed_celltypes)}\n\n")

    f.write("Marker DE:\n")
    f.write(f"  Method: {MARKER_METHOD}\n")
    f.write(f"  p_adj <= {MARKER_PADJ_THRESHOLD}\n")
    f.write(f"  logFC >= {MARKER_LOGFC_THRESHOLD}\n")
    f.write(f"  pct_in_group >= {MARKER_MIN_PCT}\n")
    f.write(f"  delta_pct >= {MARKER_MIN_DELTA_PCT}\n")
    f.write(f"  DE selection: use_raw={DE_USE_RAW}, layer={DE_USE_LAYER}\n")
    f.write(f"  Note: {de_note}\n\n")

    f.write("Gene filtering (identity marker exclusion):\n")
    f.write(f"  MT: {FILTER_MT_GENES}\n")
    f.write(f"  Ribo: {FILTER_RIBO_GENES}\n")
    f.write(f"  Unannotated: {FILTER_UNANNOTATED}\n")
    f.write(f"  Histone: {FILTER_HISTONE_GENES}\n")
    f.write(f"  IEG/Dissociation: {FILTER_STRESS_GENES}\n\n")

    f.write("Subcluster distribution:\n")
    for sub, cnt in adata.obs["subcluster"].value_counts().items():
        f.write(f"  {sub}: {cnt}\n")

    if all_markers_filtered_df is not None:
        f.write("\nMarkers (filtered):\n")
        f.write(f"  Total rows: {len(all_markers_filtered_df):,}\n")

print(f"[OK] Saved: {summary_file.name}")

# Quick reference table
ref_table = (
    adata.obs.groupby("subcluster")
    .agg({CELLTYPE_COLUMN: "first", "subcluster_leiden": "first"})
    .reset_index()
)
ref_table["n_cells"] = adata.obs["subcluster"].value_counts().loc[ref_table["subcluster"]].values
ref_table["pct_total"] = ref_table["n_cells"] / adata.n_obs * 100

if all_markers_filtered_df is not None and len(all_markers_filtered_df) > 0:
    ref_table["top_markers_filtered"] = ""
    for i, row in ref_table.iterrows():
        sub = row["subcluster"]
        sub_df = all_markers_filtered_df.loc[all_markers_filtered_df["subcluster"] == sub]
        if len(sub_df) > 0:
            top3 = sub_df.nsmallest(3, "pvals_adj")["names"].astype(str).tolist()
            ref_table.at[i, "top_markers_filtered"] = ", ".join(top3)

ref_table = ref_table.sort_values([CELLTYPE_COLUMN, "subcluster_leiden"])
ref_csv = OUTPUT_DIR / "subcluster_reference_table.csv"
ref_table.to_csv(ref_csv, index=False)
print(f"[OK] Saved: {ref_csv.name}")

print("\nReference table preview:")
print(ref_table.head(15).to_string(index=False))

# %% [markdown]
# ## Step 9: Output Validation

# %%
print("\n" + "=" * 80)
print("STEP 9: Output Validation")
print("=" * 80)

assert len(subcluster_df) == adata.n_obs, "Annotation row count mismatch"
print("[OK] Annotation rows match")

sizes = adata.obs["subcluster"].value_counts()
print(f"Subclusters: {adata.obs['subcluster'].nunique()}")
print(f"Size range: {sizes.min()} - {sizes.max()} | median={int(sizes.median())}")

tiny = sizes[sizes < 20]
if len(tiny) > 0:
    print(f"[WARN] {len(tiny)} subclusters have <20 cells (showing up to 10):")
    print(tiny.head(10).to_string())

if all_markers_filtered_df is not None:
    mt_left = all_markers_filtered_df["names"].astype(str).str.startswith(("MT-", "mt-")).sum()
    ribo_left = all_markers_filtered_df["names"].astype(str).str.startswith(("RPL", "RPS", "MRPL", "MRPS")).sum()
    if mt_left > 0 or ribo_left > 0:
        print(f"[WARN] Filtered markers still contain MT/Ribo: MT={mt_left}, Ribo={ribo_left}")
    else:
        print("[OK] Filtered markers contain no MT/Ribo (by prefix check)")

print("\n" + "=" * 80)
print(f"PIPELINE COMPLETE ({PIPELINE_VERSION})")
print("=" * 80)
print(f"Output directory: {OUTPUT_DIR}")
print("Main outputs:")
print("  - epithelial_with_subclusters.h5ad")
print("  - subcluster_annotations.csv")
print("  - subcluster_reference_table.csv")
print("  - subcluster_summary.txt")
print("  - filtered_genes_list.csv")
if len(all_markers) > 0:
    print("  - markers_all_subclusters_unfiltered.csv")
if len(all_markers_filtered) > 0:
    print("  - markers_all_subclusters_filtered.csv")
    print(f"  - markers_top{TOP_N_MARKERS}_per_subcluster_filtered.csv")
    print("  - markers_[celltype]_filtered.csv")
print("Figures in figures/:")
print("  - umap_subclusters_all")
print("  - umap_comparison_celltype_vs_subcluster")
if BATCH_KEY:
    print("  - umap_batch_distribution")
print("  - umap_subcluster_[celltype]")
print("  - dotplot/heatmap markers (if available)")
