# %% [markdown]
# # T/NK Cell Subcluster Annotation + scANVI Retrain
# Version: v1.2 (2026-03-18)
# 
# **Changes vs v1.1**:
#   [v1.2-1] SPLIT_KEY changed to "scanvi_label" (L3-level subclustering).
#            Each L3 label is subclustered independently so users can
#            split/merge/relabel at fine-grained level within each known type.
#   [v1.2-2] TARGET_CELL_TYPES=None auto-detects all L3 categories.
#            Guard added: auto-expand before verification check (fixes NoneType error).
#   [v1.2-3] BACKUP_LABELS_KEY: original scanvi_label backed up before overwrite.
#   [v1.2-4] SUBCLUSTER_REP = "X_scanvi" (supervision-shaped latent, better
#            within-lineage separation for state discovery).
#   [v1.2-5] matrixplot: color_map= -> cmap= (matplotlib compatibility fix).
# 
# **Workflow**:
#   1. Load h5ad (post-scVI/scANVI reference)
#   2. For each L3 label in scanvi_label:
#      - Extract cells, compute neighbors on X_scanvi, run Leiden r=1.0
#      - Compute markers (use_raw=True), save CSV
#      - Save UMAP + known-marker dotplot + featureplot + data CSVs
#   3. USER reviews outputs and fills annotation maps (Part 5)
#   4. Original scanvi_label backed up; new labels applied from map
#   5. scVI reloaded (not retrained); scANVI retrained with new labels
#   6. Output: refined h5ad + new scANVI model
# 
# **Input**:
#   - adata_tnk_scanvi_ref_20260315_v1_2.h5ad
#   - tnk_scvi_ref_model/  (original scVI; not retrained)
# 
# **Output**:
#   - subcluster_markers/<L3_type>/markers_leiden_r1.0.csv
#   - subcluster_markers/<L3_type>/dotplot_*.pdf + matrixplot_*.pdf
#   - subcluster_markers/<L3_type>/featureplot_known_markers.pdf
#   - subcluster_markers/<L3_type>/dotplot_data_r1.0_*.csv
#   - adata_tnk_scanvi_ref_retrain_v1_2.h5ad
#   - tnk_scanvi_ref_model_retrain/
#   NOTE: scVI is NOT retrained. Reuse original SCVI_MODEL_DIR for scArches queries.
# 
# **QRM rules applied**:
#   - basis='umap' (NOT 'X_umap') in sc.pl.embedding
#   - sc.settings.vector_friendly=True before plot blocks; no rasterized= kwarg
#   - use_raw=True for rank_genes_groups
#   - QRM 15.1 : from_scvi_model for reference training
#   - QRM 13   : category dtype before write_h5ad

# %% [markdown]
# ## 0. Configuration

# %%
# ============================================================================
# CONFIGURATION  (edit here only)
# ============================================================================

INPUT_H5AD     = "/home/h2048/data/py/0315/tnk_scarches_ref/adata_tnk_scanvi_ref_20260315_v1_2.h5ad"
SCVI_MODEL_DIR = "/home/h2048/data/py/0315/tnk_scarches_ref/tnk_scvi_ref_model"
OUTPUT_DIR     = "/home/h2048/data/py/0318/tnk_subcluster_retrain"

# Split by L3 label so each known cell type is subclustered independently
SPLIT_KEY      = "scanvi_label"

# None = auto-detect all categories in SPLIT_KEY
TARGET_CELL_TYPES = None

# Latent space for within-type subclustering
# X_scanvi: supervision-shaped latent, better within-lineage separation
SUBCLUSTER_REP = "X_scanvi"

# Single resolution
LEIDEN_RESOLUTION = 0.3

# Existing scANVI labels (backed up, then overwritten with refined labels)
OLD_LABELS_KEY    = "scanvi_label"
BACKUP_LABELS_KEY = "scanvi_label_original"   # backup column name
NEW_LABELS_KEY    = "scanvi_label_refined"
UNLABELED         = "Unknown"

# Batch key (for scVI reload)
BATCH_KEY      = "sample"

# scANVI retrain
SCANVI_EPOCHS  = 200
BATCH_SIZE     = 256
N_SAMPLES_PER_LABEL = None   # None = auto

# Marker detection
MARKER_METHOD         = "wilcoxon"
N_TOP_MARKERS         = 20
MIN_IN_GROUP_FRACTION = 0.1

# Visualization
DPI  = 300
SEED = 42

print("Configuration loaded")
print(f"  Input h5ad    : {INPUT_H5AD}")
print(f"  scVI model    : {SCVI_MODEL_DIR}")
print(f"  Output        : {OUTPUT_DIR}")
print(f"  Split key     : {SPLIT_KEY}  (L3-level subclustering)")
print(f"  Subcluster rep: {SUBCLUSTER_REP}")

# %% [markdown]
# ## 1. Imports & Setup

# %%
import os, gc, warnings
warnings.filterwarnings("ignore")

os.environ["OMP_NUM_THREADS"]      = "4"
os.environ["MKL_NUM_THREADS"]      = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import scanpy as sc
import scvi
from pathlib import Path

sc.settings.verbosity = 2
sc.settings.n_jobs    = 8
sc.settings.set_figure_params(dpi=DPI, facecolor="white", frameon=False)

np.random.seed(SEED)
scvi.settings.seed = SEED

output_dir = Path(OUTPUT_DIR)
marker_dir = output_dir / "subcluster_markers"
fig_dir    = output_dir / "figures"
for d in [output_dir, marker_dir, fig_dir]:
    d.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(fig_dir)

print(f"scanpy : {sc.__version__}")
print(f"scvi   : {scvi.__version__}")

# %% [markdown]
# ## 2. Load Data

# %%
print("Loading h5ad...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  Shape  : {adata.shape}")
print(f"  obsm   : {list(adata.obsm.keys())}")
print(f"  layers : {list(adata.layers.keys()) if adata.layers else 'None'}")
if adata.raw is not None:
    print(f"  .raw   : {adata.raw.n_obs} x {adata.raw.n_vars} genes")

assert SUBCLUSTER_REP in adata.obsm, \
    f"'{SUBCLUSTER_REP}' not in obsm. Available: {list(adata.obsm.keys())}"
assert SPLIT_KEY in adata.obs.columns, \
    f"'{SPLIT_KEY}' not in obs.columns"
assert OLD_LABELS_KEY in adata.obs.columns, \
    f"'{OLD_LABELS_KEY}' not in obs.columns"

# Auto-detect TARGET_CELL_TYPES before any verification
if TARGET_CELL_TYPES is None:
    TARGET_CELL_TYPES = sorted(
        adata.obs[SPLIT_KEY].dropna().astype(str).unique().tolist()
    )
    # Exclude UNLABELED from subclustering targets
    TARGET_CELL_TYPES = [ct for ct in TARGET_CELL_TYPES if ct != UNLABELED]
    print(f"\nAuto-detected {len(TARGET_CELL_TYPES)} target types from '{SPLIT_KEY}':")
    for ct in TARGET_CELL_TYPES:
        n = (adata.obs[SPLIT_KEY].astype(str) == ct).sum()
        print(f"  {ct}: {n:,} cells")
else:
    missing = [ct for ct in TARGET_CELL_TYPES
               if ct not in adata.obs[SPLIT_KEY].values]
    if missing:
        raise ValueError(
            f"TARGET_CELL_TYPES not found in obs['{SPLIT_KEY}']: {missing}\n"
            f"Available: {sorted(adata.obs[SPLIT_KEY].unique().tolist())}"
        )
    print(f"\nTarget cell types ({len(TARGET_CELL_TYPES)}):")
    for ct in TARGET_CELL_TYPES:
        n = (adata.obs[SPLIT_KEY].astype(str) == ct).sum()
        print(f"  {ct}: {n:,} cells")

# %% [markdown]
# ## 3. Per-L3-Type Subclustering on X_scanvi Space
# 
# For each L3 label in scanvi_label:
#   1. Extract cells, compute within-type neighbors on X_scanvi
#   2. Run Leiden r=1.0
#   3. Compute rank_genes_groups markers (use_raw=True)
#   4. Save marker CSV, UMAP, known-marker dotplot + matrixplot + featureplot
#   5. Store cluster IDs back into main adata.obs

# %%
def sanitize_colname(s):
    return (s.lower()
             .replace(" ", "_")
             .replace("/", "_")
             .replace("+", "plus")
             .replace("-", "minus")
             .replace("(", "")
             .replace(")", ""))

LEIDEN_KEY = f"leiden_r{LEIDEN_RESOLUTION}"

subcluster_summary = {}   # {cell_type: {annotation_key, n_subclusters}}
res_results_global = {}   # {cell_type: n_clusters}

KNOWN_MARKERS_ORDERED = [
    # Lineage
    "CD3D", "CD3E",
    "CD4",
    "CD8A", "CD8B",
    "GNLY", "NKG7", "NCAM1", "FCGR3A",
    "TRDC", "TRGC1",
    "KLRB1", "ZBTB16",
    "KIT", "IL7R",
    # Naive / central memory
    "CCR7", "SELL", "TCF7", "LEF1", "MAL", "IL6R",
    # Effector / TEM
    "GZMK", "GZMA", "GZMH", "GZMB", "PRF1",
    "IFNG", "TNF",
    # TRM
    "ITGAE", "ITGA1", "CD69", "CXCR6", "RGS1", "ZNF683",
    # TEMRA / terminal
    "KLRG1", "CX3CR1", "S1PR5", "ZEB2", "FGFBP2",
    # Exhaustion / inhibitory
    "PDCD1", "HAVCR2", "TIGIT", "LAG3", "TOX", "ENTPD1",
    # Tfh / Treg
    "CXCR5", "CXCL13", "TOX2",
    "FOXP3", "IL2RA", "CTLA4",
    # Th1 / Th17
    "TBX21", "RORC", "IL23R", "CCR6",
    # NK functional
    "XCL1", "XCL2", "AREG",
    "ISG15", "IFIT1", "MX1",
    # Proliferation
    "MKI67", "TOP2A",
    # MT-high flag
    "MT1E", "MT1X",
]
_seen = set()
KNOWN_MARKERS_ORDERED = [g for g in KNOWN_MARKERS_ORDERED
                          if not (g in _seen or _seen.add(g))]

FEATURE_PRIORITY = [
    "CD3D", "CD4", "CD8A", "GNLY", "NKG7", "FCGR3A",
    "CCR7", "TCF7", "GZMB", "GZMK", "PRF1",
    "ITGAE", "CD69", "CXCR6", "RGS1",
    "KLRG1", "CX3CR1",
    "PDCD1", "HAVCR2", "FOXP3",
    "CXCL13", "TOX",
    "TRDC", "KLRB1", "ZBTB16", "KIT", "MKI67",
]

for cell_type in TARGET_CELL_TYPES:
    print(f"\n{'='*60}")
    print(f"Subclustering: {cell_type}")
    print(f"{'='*60}")

    mask    = adata.obs[SPLIT_KEY].astype(str) == cell_type
    adata_s = adata[mask].copy()
    n_cells = adata_s.n_obs
    print(f"  Cells: {n_cells:,}")

    if n_cells < 30:
        print(f"  [SKIP] < 30 cells -- too few for meaningful subclustering")
        continue

    ct_safe = sanitize_colname(cell_type)
    ct_mdir = marker_dir / ct_safe
    ct_mdir.mkdir(exist_ok=True)

    # Within-subset neighbors on X_scanvi
    n_neighbors = min(30, max(5, n_cells // 5))
    print(f"  Neighbors on {SUBCLUSTER_REP} (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata_s, use_rep=SUBCLUSTER_REP, n_neighbors=n_neighbors)
    sc.tl.umap(adata_s, min_dist=0.3, spread=1.0)

    # Leiden at fixed resolution
    try:
        sc.tl.leiden(adata_s, resolution=LEIDEN_RESOLUTION, key_added=LEIDEN_KEY)
        n_clust = adata_s.obs[LEIDEN_KEY].nunique()
        res_results_global[cell_type] = n_clust
        print(f"  Leiden r={LEIDEN_RESOLUTION}: {n_clust} clusters")
    except Exception as e:
        print(f"  Leiden failed: {e} -- skipping {cell_type}")
        del adata_s; gc.collect()
        continue

    if n_clust <= 1:
        print(f"  [SKIP] Only {n_clust} cluster -- no subclustering possible")
        del adata_s; gc.collect()
        continue

    # Marker genes
    print(f"  Computing markers ({n_clust} clusters)...")
    try:
        sc.tl.rank_genes_groups(
            adata_s,
            groupby   = LEIDEN_KEY,
            method    = MARKER_METHOD,
            use_raw   = True,
            pts       = True,
            key_added = f"rgg_{LEIDEN_KEY}",
        )
        all_markers = []
        for clust in adata_s.obs[LEIDEN_KEY].cat.categories:
            df = sc.get.rank_genes_groups_df(
                adata_s, group=clust, key=f"rgg_{LEIDEN_KEY}",
                pval_cutoff=0.05, log2fc_min=0.25,
            )
            if df.empty:
                continue
            if "pct_nz_group" in df.columns:
                df = df[df["pct_nz_group"] >= MIN_IN_GROUP_FRACTION]
            df.insert(0, "cluster",   clust)
            df.insert(1, "cell_type", cell_type)
            all_markers.append(df.head(N_TOP_MARKERS))

        if all_markers:
            markers_df = pd.concat(all_markers, ignore_index=True)
            markers_df.to_csv(ct_mdir / f"markers_leiden_r{LEIDEN_RESOLUTION}.csv",
                              index=False)
            print(f"    Saved: markers_leiden_r{LEIDEN_RESOLUTION}.csv")

    except Exception as e:
        print(f"    Markers failed: {e}")

    # UMAP overview
    sc.settings.vector_friendly = True
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sc.pl.embedding(adata_s, basis="umap", color=OLD_LABELS_KEY,
                    title=f"Original label: {cell_type}",
                    ax=axes[0], show=False,
                    legend_loc="right margin", legend_fontsize=7)

    sc.pl.embedding(adata_s, basis="umap", color=LEIDEN_KEY,
                    title=f"Leiden r={LEIDEN_RESOLUTION} ({n_clust} clusters)",
                    ax=axes[1], show=False,
                    legend_loc="on data", legend_fontsize=8)

    plt.suptitle(f"{cell_type} -- Subclusters", fontsize=13, y=1.01)
    plt.tight_layout()
    fig.savefig(ct_mdir / "umap_subclusters.pdf", bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print(f"  Saved: umap_subclusters.pdf")

    # Data-driven dotplot
    rgg_key = f"rgg_{LEIDEN_KEY}"
    if rgg_key in adata_s.uns:
        try:
            sc.settings.vector_friendly = True
            dp = sc.pl.rank_genes_groups_dotplot(
                adata_s, n_genes=5, key=rgg_key,
                groupby=LEIDEN_KEY, use_raw=True,
                show=False, return_fig=True,
            )
            dp.savefig(ct_mdir / f"dotplot_datadriven_r{LEIDEN_RESOLUTION}.pdf",
                       bbox_inches="tight", dpi=DPI)
            plt.close()
            print(f"  Saved: dotplot_datadriven_r{LEIDEN_RESOLUTION}.pdf")
        except Exception as e:
            print(f"  Data-driven dotplot failed: {e}")

    # Known-marker dotplot + matrixplot
    avail       = (set(adata_s.raw.var_names) if adata_s.raw is not None
                   else set(adata_s.var_names))
    known_valid = [g for g in KNOWN_MARKERS_ORDERED if g in avail]
    known_miss  = [g for g in KNOWN_MARKERS_ORDERED if g not in avail]
    if known_miss:
        print(f"  Markers absent from .raw ({len(known_miss)}): "
              f"{known_miss[:6]}{'...' if len(known_miss) > 6 else ''}")
    print(f"  Known markers: {len(known_valid)}/{len(KNOWN_MARKERS_ORDERED)}")

    if known_valid:
        fig_w = max(14, len(known_valid) * 0.38)
        fig_h = max(4,  n_clust * 0.6)
        try:
            sc.settings.vector_friendly = True
            fig_dp = sc.pl.dotplot(
                adata_s, var_names=known_valid, groupby=LEIDEN_KEY,
                use_raw=True, standard_scale="var",
                dot_max=1.0, color_map="Blues",
                show=False, return_fig=True,
                figsize=(fig_w, fig_h),
                title=f"{cell_type}  |  Known Markers  |  Leiden r={LEIDEN_RESOLUTION}",
            )
            fig_dp.savefig(ct_mdir / f"dotplot_known_markers_r{LEIDEN_RESOLUTION}.pdf",
                           bbox_inches="tight", dpi=DPI)
            plt.close()
            print(f"  Saved: dotplot_known_markers_r{LEIDEN_RESOLUTION}.pdf")
        except Exception as e:
            print(f"  Known-marker dotplot failed: {e}")

        try:
            sc.settings.vector_friendly = True
            fig_mp = sc.pl.matrixplot(
                adata_s, var_names=known_valid, groupby=LEIDEN_KEY,
                use_raw=True, standard_scale="var",
                cmap="RdBu_r",
                show=False, return_fig=True,
                figsize=(fig_w, fig_h),
                title=f"{cell_type}  |  Known Markers (matrix)  |  Leiden r={LEIDEN_RESOLUTION}",
            )
            fig_mp.savefig(ct_mdir / f"matrixplot_known_markers_r{LEIDEN_RESOLUTION}.pdf",
                           bbox_inches="tight", dpi=DPI)
            plt.close()
            print(f"  Saved: matrixplot_known_markers_r{LEIDEN_RESOLUTION}.pdf")
        except Exception as e:
            print(f"  Matrixplot failed: {e}")

    # Feature plots
    feat_valid = [g for g in FEATURE_PRIORITY if g in avail]
    if feat_valid:
        try:
            sc.settings.vector_friendly = True
            ncols_f = 6
            nrows_f = int(np.ceil(len(feat_valid) / ncols_f))
            fig_f, axes_f = plt.subplots(
                nrows_f, ncols_f, figsize=(ncols_f * 3.5, nrows_f * 3.2))
            axes_ff = np.array(axes_f).ravel()
            for idx, gene in enumerate(feat_valid):
                sc.pl.embedding(
                    adata_s, basis="umap", color=gene,
                    ax=axes_ff[idx], show=False,
                    use_raw=(adata_s.raw is not None),
                    color_map="viridis", vmin=0,
                )
                axes_ff[idx].set_title(gene, fontsize=9, fontweight="bold")
            for idx in range(len(feat_valid), len(axes_ff)):
                axes_ff[idx].set_visible(False)
            plt.suptitle(f"{cell_type} -- Key Markers (use_raw=True)",
                         fontsize=12, y=1.01)
            plt.tight_layout()
            fig_f.savefig(ct_mdir / "featureplot_known_markers.pdf",
                          bbox_inches="tight", dpi=DPI)
            plt.close(fig_f)
            print(f"  Saved: featureplot_known_markers.pdf")
        except Exception as e:
            print(f"  Feature plots failed: {e}")

    # Store cluster IDs into main adata.obs
    main_col = f"subcluster_{ct_safe}"
    adata.obs.loc[mask, main_col] = (
        adata_s.obs[LEIDEN_KEY].astype(str)
        .apply(lambda x: f"{ct_safe}_{x}").values
    )
    subcluster_summary[cell_type] = {
        "n_cells"        : n_cells,
        "annotation_key" : main_col,
        "n_subclusters"  : n_clust,
    }
    print(f"  Stored in adata.obs['{main_col}']: {n_clust} subclusters")

    del adata_s
    gc.collect()

print("\n" + "="*60)
print("Subclustering complete. Summary:")
for ct, info in subcluster_summary.items():
    print(f"  {ct}: {info['n_subclusters']} subclusters -> obs['{info['annotation_key']}']")

# %% [markdown]
# ## 3b. Export Dotplot Underlying Data as CSV
# 
# For each L3 type, saves three tables (rows=clusters, cols=marker genes):
#   - dotplot_data_r1.0_mean_expr_raw.csv     : log1p mean (unscaled)
#   - dotplot_data_r1.0_mean_expr_scaled.csv  : 0-1 per gene (= dotplot color)
#   - dotplot_data_r1.0_pct_expressing.csv    : fraction cells > 0 (= dot size)

# %%
def _get_expr_matrix(adata_s, genes):
    if adata_s.raw is not None:
        src  = adata_s.raw
        idx  = src.var_names.get_indexer(genes)
        ok   = idx >= 0
        X    = src.X[:, idx[ok]]
        genes_used = [genes[i] for i, v in enumerate(ok) if v]
    else:
        idx  = adata_s.var_names.get_indexer(genes)
        ok   = idx >= 0
        X    = adata_s.X[:, idx[ok]]
        genes_used = [genes[i] for i, v in enumerate(ok) if v]
    if sparse.issparse(X):
        X = np.asarray(X.todense())
    return X.astype(np.float32), genes_used

print("Exporting dotplot data as CSV...")

for cell_type, info in subcluster_summary.items():
    ct_safe = sanitize_colname(cell_type)
    ct_mdir = marker_dir / ct_safe
    mask    = adata.obs[SPLIT_KEY].astype(str) == cell_type
    adata_s = adata[mask].copy()

    avail_s     = (set(adata_s.raw.var_names) if adata_s.raw is not None
                   else set(adata_s.var_names))
    known_valid = [g for g in KNOWN_MARKERS_ORDERED if g in avail_s]

    if not known_valid:
        print(f"  {cell_type}: no known markers -- skipping CSV")
        del adata_s; gc.collect()
        continue

    main_col = info["annotation_key"]
    if main_col not in adata_s.obs.columns:
        # Reconstruct from main adata
        adata_s.obs[main_col] = adata.obs.loc[adata_s.obs_names, main_col]

    if main_col not in adata_s.obs.columns or adata_s.obs[main_col].isna().all():
        print(f"  {cell_type}: cluster column missing -- skipping CSV")
        del adata_s; gc.collect()
        continue

    clusters = sorted(adata_s.obs[main_col].dropna().unique().tolist())
    if len(clusters) <= 1:
        del adata_s; gc.collect()
        continue

    X_full, genes_used = _get_expr_matrix(adata_s, known_valid)

    rows_mean_raw = {}
    rows_pct      = {}
    for clust in clusters:
        cmask = (adata_s.obs[main_col] == clust).values
        X_c   = X_full[cmask, :]
        rows_mean_raw[clust] = X_c.mean(axis=0)
        rows_pct[clust]      = (X_c > 0).mean(axis=0)

    df_mean_raw = pd.DataFrame(rows_mean_raw, index=genes_used).T
    df_pct      = pd.DataFrame(rows_pct,      index=genes_used).T
    gene_min    = df_mean_raw.min(axis=0)
    gene_range  = (df_mean_raw.max(axis=0) - gene_min).replace(0, 1)
    df_scaled   = (df_mean_raw - gene_min) / gene_range

    prefix = ct_mdir / f"dotplot_data_r{LEIDEN_RESOLUTION}"
    for df, suffix in [
        (df_mean_raw.round(4), "mean_expr_raw"),
        (df_scaled.round(4),   "mean_expr_scaled"),
        (df_pct.round(4),      "pct_expressing"),
    ]:
        df.index.name = "cluster"
        df.to_csv(f"{prefix}_{suffix}.csv")

    print(f"  {cell_type}: {len(clusters)} clusters x {len(genes_used)} genes -- saved")

    del adata_s, X_full
    gc.collect()

print(f"\nCSV export complete -> {marker_dir}/<cell_type>/dotplot_data_*.csv")

# %%
# ============================================================
# Merge all dotplot CSVs into a single long-format table
# Columns: cell_type | cluster | gene | mean_expr_scaled | pct_expressing
# One row per (cluster x gene) -- suitable for manual dotplot in R/Python
# ============================================================
import pandas as pd
from pathlib import Path

merged_rows = []

for cell_type, info in subcluster_summary.items():
    ct_safe = sanitize_colname(cell_type)
    ct_mdir = marker_dir / ct_safe
    res     = LEIDEN_RESOLUTION

    f_scaled = ct_mdir / f"dotplot_data_r{res}_mean_expr_scaled.csv"
    f_pct    = ct_mdir / f"dotplot_data_r{res}_pct_expressing.csv"

    if not f_scaled.exists() or not f_pct.exists():
        print(f"  [SKIP] {cell_type}: CSV files not found")
        continue

    df_scaled = pd.read_csv(f_scaled, index_col="cluster")
    df_pct    = pd.read_csv(f_pct,    index_col="cluster")

    # Melt to long format
    df_s_long = df_scaled.reset_index().melt(
        id_vars="cluster", var_name="gene", value_name="mean_expr_scaled"
    )
    df_p_long = df_pct.reset_index().melt(
        id_vars="cluster", var_name="gene", value_name="pct_expressing"
    )

    df_merged = df_s_long.merge(df_p_long, on=["cluster", "gene"])
    df_merged.insert(0, "cell_type", cell_type)
    merged_rows.append(df_merged)

    print(f"  {cell_type}: {df_merged.shape[0]} rows")

df_all = pd.concat(merged_rows, ignore_index=True)

# Sort for readability
df_all = df_all.sort_values(["cell_type", "cluster", "gene"]).reset_index(drop=True)

out_csv = output_dir / f"dotplot_data_all_r{LEIDEN_RESOLUTION}_combined.csv"
df_all.to_csv(out_csv, index=False)

print(f"\nCombined CSV saved: {out_csv}")
print(f"  Shape : {df_all.shape}  ({df_all['cell_type'].nunique()} cell types, "
      f"{df_all['cluster'].nunique()} clusters, {df_all['gene'].nunique()} genes)")
print(f"\nPreview:")
print(df_all.head(10).to_string(index=False))

# %% [markdown]
# ## 4. Print Annotation Templates
# 
# Run this cell after Part 3.
# Each subcluster ID follows the pattern: {l3_type_safe}_{cluster_number}
# Review marker CSVs and dotplots, then fill in labels in Part 5.

# %%
print("="*70)
print("ANNOTATION TEMPLATES")
print("Review subcluster_markers/<type>/ then fill labels below.")
print("="*70)
print()
print("# Paste this into Part 5 and fill in the label strings.")
print()

for cell_type, info in subcluster_summary.items():
    ct_safe  = sanitize_colname(cell_type)
    anno_key = info["annotation_key"]
    n_sub    = info["n_subclusters"]

    print(f"# --- {cell_type} ({n_sub} subclusters) ---")
    print(f"ANNOTATION_MAP_{ct_safe} = {{")
    if anno_key in adata.obs.columns:
        for c in sorted(adata.obs[anno_key].dropna().unique().tolist()):
            print(f'    "{c}": "",   # fill in label')
    print("}")
    print()

# %% [markdown]
# ## 5. USER FILLS ANNOTATION MAPS
# 
# Instructions:
#   1. Review PDFs and CSVs in subcluster_markers/<l3_type>/
#   2. For each cluster, fill in the new label string
#   3. Leave "" or set to UNLABELED = "Unknown" to keep as unlabeled
#   4. If splitting one L3 into two new types, use different label strings
#   5. If merging clusters, use the same label string for multiple entries
#   6. Run this cell, then continue to Part 6

# %%
# ============================================================
# USER FILLS: annotation map (revised from combined table)
# ============================================================

ANNOTATION_MAP_cd4_tem = {
    "cd4_tem_0": "CD4 Trm",
    "cd4_tem_1": "CD4 Tcm",
}

ANNOTATION_MAP_cd4_tfh = {
    "cd4_tfh_0": "CD4 Tfh",
    "cd4_tfh_1": "CD4 Tfh",
}

ANNOTATION_MAP_cd4_treg = {
    "cd4_treg_0": "CD4 Treg",
    "cd4_treg_1": "CD4 Treg",
    "cd4_treg_2": "CD4 Tfr",
}

ANNOTATION_MAP_cd8_gammaminusdelta = {
    "cd8_gammaminusdelta_0": "gdT",
    "cd8_gammaminusdelta_1": "gdT",
}

ANNOTATION_MAP_cd8_mait = {
    "cd8_mait_0": "MAIT",
    "cd8_mait_1": "MAIT",
}

ANNOTATION_MAP_cd8_naive_tcm = {
    "cd8_naive_tcm_0": "CD8 Trm",
    "cd8_naive_tcm_1": "CD8 Naive",
}

ANNOTATION_MAP_cd8_temra = {
    "cd8_temra_0": "CD8 Temra",
    "cd8_temra_1": "CD8 Tem",
    "cd8_temra_2": "CD8 Teff",
}

ANNOTATION_MAP_cd8_trm = {
    "cd8_trm_0": "CD8 Trm",
    "cd8_trm_1": "CD8 Trm ",
    "cd8_trm_2": "CD8 Trm",
    "cd8_trm_3": "MAIT",
}

ANNOTATION_MAP_nk_cytotoxic = {
    "nk_cytotoxic_0": "NK",
    "nk_cytotoxic_1": "NK",
    "nk_cytotoxic_2": "gdT",
    "nk_cytotoxic_3": "CD8 Teff",
}

ANNOTATION_MAP_nk_inflammatory = {
    "nk_inflammatory_0": "CD8 Trm",
    "nk_inflammatory_1": "CD8 Trm",
    "nk_inflammatory_2": "gdT",
    "nk_inflammatory_3": "CD8 Teff",
    "nk_inflammatory_4": "CD8 Teff",
}

ANNOTATION_MAP_nk_resting = {
    "nk_resting_0": "NK",
    "nk_resting_1": "NK",
    "nk_resting_2": "NK",
}

COMBINED_ANNOTATION_MAP = {}
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd4_tem)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd4_tfh)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd4_treg)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd8_gammaminusdelta)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd8_mait)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd8_naive_tcm)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd8_temra)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_cd8_trm)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_nk_cytotoxic)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_nk_inflammatory)
COMBINED_ANNOTATION_MAP.update(ANNOTATION_MAP_nk_resting)

n_filled = sum(1 for v in COMBINED_ANNOTATION_MAP.values() if v.strip() != "")
n_total  = len(COMBINED_ANNOTATION_MAP)
print(f"Annotation map: {n_filled}/{n_total} entries filled")

# %% [markdown]
# ## 6. Apply Annotations -> Build NEW_LABELS_KEY
# 
# Strategy:
#   1. Backup OLD_LABELS_KEY -> BACKUP_LABELS_KEY
#   2. NEW_LABELS_KEY inherits OLD_LABELS_KEY (all cells start with original label)
#   3. Only cells in subcluster_summary types are cleared then re-filled from map
#   4. Non-subclustered types keep their original labels unchanged

# %%
print(f"Backing up '{OLD_LABELS_KEY}' -> '{BACKUP_LABELS_KEY}'...")
adata.obs[BACKUP_LABELS_KEY] = adata.obs[OLD_LABELS_KEY].astype(str).copy()
print(f"  Backup done: {adata.obs[BACKUP_LABELS_KEY].nunique()} original categories")

print(f"\nBuilding '{NEW_LABELS_KEY}'...")

# Step 1: Inherit original labels as starting point
adata.obs[NEW_LABELS_KEY] = adata.obs[OLD_LABELS_KEY].astype(str)

# Step 2: Clear only the subclustered types (will be re-filled from map)
for cell_type in subcluster_summary:
    mask_ct = adata.obs[SPLIT_KEY].astype(str) == cell_type
    adata.obs.loc[mask_ct, NEW_LABELS_KEY] = UNLABELED

n_cleared = (adata.obs[NEW_LABELS_KEY] == UNLABELED).sum()
print(f"  Cleared to '{UNLABELED}' (subclustered types only): {n_cleared:,} cells")

# Step 3: Fill from COMBINED_ANNOTATION_MAP
n_mapped_labeled   = 0
n_mapped_unlabeled = 0
n_missing_keys     = 0

for cell_type, info in subcluster_summary.items():
    anno_key = info["annotation_key"]
    ct_safe  = sanitize_colname(cell_type)

    if anno_key not in adata.obs.columns:
        print(f"  [WARN] '{anno_key}' not in obs -- {cell_type} stays Unknown")
        continue

    for cluster_id, label in COMBINED_ANNOTATION_MAP.items():
        if not cluster_id.startswith(ct_safe + "_"):
            continue
        label = label.strip() if label.strip() != "" else UNLABELED

        mask = adata.obs[anno_key] == cluster_id
        n    = int(mask.sum())
        if n == 0:
            print(f"  [WARN] Cluster '{cluster_id}' not found in obs['{anno_key}']")
            n_missing_keys += 1
            continue

        adata.obs.loc[mask, NEW_LABELS_KEY] = label
        if label != UNLABELED:
            n_mapped_labeled   += n
        else:
            n_mapped_unlabeled += n

n_total_labeled   = (adata.obs[NEW_LABELS_KEY] != UNLABELED).sum()
n_total_unlabeled = (adata.obs[NEW_LABELS_KEY] == UNLABELED).sum()

print(f"\nAssignment summary:")
print(f"  Total cells           : {adata.n_obs:,}")
print(f"  Labeled (non-Unknown) : {n_total_labeled:,}")
print(f"  Unknown               : {n_total_unlabeled:,}")
print(f"  Mapped labeled        : {n_mapped_labeled:,}")
print(f"  Mapped Unknown        : {n_mapped_unlabeled:,}")
print(f"  Unresolved map keys   : {n_missing_keys}")
print(f"\nFull label distribution:")
print(adata.obs[NEW_LABELS_KEY].value_counts())

adata.obs[NEW_LABELS_KEY] = adata.obs[NEW_LABELS_KEY].astype("category")
if UNLABELED not in adata.obs[NEW_LABELS_KEY].cat.categories:
    adata.obs[NEW_LABELS_KEY] = adata.obs[NEW_LABELS_KEY].cat.add_categories([UNLABELED])

# %% [markdown]
# ## 7. Sanity Check UMAP Before Retrain

# %%
if "X_umap" not in adata.obsm:
    print("X_umap not found -- computing from X_scanvi...")
    sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.3)

sc.settings.vector_friendly = True
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

sc.pl.embedding(adata, basis="umap", color=BACKUP_LABELS_KEY,
                title=f"Original labels ({BACKUP_LABELS_KEY})",
                ax=axes[0], show=False,
                legend_loc="right margin", legend_fontsize=6)

sc.pl.embedding(adata, basis="umap", color=NEW_LABELS_KEY,
                title=f"New labels ({NEW_LABELS_KEY})",
                ax=axes[1], show=False,
                legend_loc="right margin", legend_fontsize=6)

# Unknown cells highlighted
unknown_col = (adata.obs[NEW_LABELS_KEY].astype(str) == UNLABELED).map(
    {True: "Unknown", False: "Labeled"}).astype("category")
adata.obs["_unknown_check"] = unknown_col
sc.pl.embedding(adata, basis="umap", color="_unknown_check",
                title="Unknown vs Labeled",
                ax=axes[2], show=False,
                palette={"Unknown": "#e74c3c", "Labeled": "#95a5a6"})
adata.obs.drop(columns=["_unknown_check"], inplace=True)

plt.suptitle("Annotation Sanity Check Before scANVI Retrain", fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig(fig_dir / "annotation_sanity_check.pdf", bbox_inches="tight", dpi=DPI)
plt.close(fig)
print("Saved: annotation_sanity_check.pdf")
print("\nVerify labels look correct before continuing to Part 8.")

# %% [markdown]
# ## 8. Guard: Validate Label Distribution

# %%
label_counts = adata.obs[NEW_LABELS_KEY].value_counts()
labeled_only = label_counts[label_counts.index != UNLABELED]

print(f"Label validation:")
print(f"  Total classes (incl. Unknown) : {len(label_counts)}")
print(f"  Labeled classes               : {len(labeled_only)}")

if len(labeled_only) == 0:
    raise ValueError(
        f"No labeled classes in '{NEW_LABELS_KEY}'. "
        "Fill annotation maps in Part 5 before retraining."
    )
if len(labeled_only) < 2:
    raise ValueError(
        f"Only {len(labeled_only)} labeled class in '{NEW_LABELS_KEY}'. "
        "scANVI requires at least 2 supervised classes."
    )

smallest_class = labeled_only.idxmin()
min_class_size = int(labeled_only.min())
print(f"  Smallest class: '{smallest_class}' ({min_class_size} cells)")
print(f"\nPer-class distribution:")
print(labeled_only)

if N_SAMPLES_PER_LABEL is None:
    n_samples_per_label = min(100, max(2, int(min_class_size * 0.8)))
    n_samples_per_label = min(n_samples_per_label, min_class_size)
else:
    n_samples_per_label = N_SAMPLES_PER_LABEL
print(f"\nn_samples_per_label (adaptive): {n_samples_per_label}")
print("Validation passed.")

# %% [markdown]
# ## 9. Reload scVI + Train scANVI with New Labels

# %%
import torch

print(f"Loading scVI model from: {SCVI_MODEL_DIR}")
assert "counts" in adata.layers, \
    "'counts' layer not found. Expected HVG-subset h5ad."

if UNLABELED not in adata.obs[NEW_LABELS_KEY].cat.categories:
    adata.obs[NEW_LABELS_KEY] = adata.obs[NEW_LABELS_KEY].cat.add_categories([UNLABELED])

vae = scvi.model.SCVI.load(SCVI_MODEL_DIR, adata=adata)
print(f"  scVI loaded: n_latent={vae.module.n_latent}")

# Latent consistency check vs X_scvi (not X_scanvi -- different spaces)
if "X_scvi" not in adata.obsm:
    print("  [INFO] 'X_scvi' not found -- skipping consistency check")
else:
    n_check     = min(500, adata.n_obs)
    test_latent = vae.get_latent_representation(indices=list(range(n_check)))
    existing    = adata.obsm["X_scvi"][:n_check]
    corr        = float(np.corrcoef(test_latent.ravel(), existing.ravel())[0, 1])
    print(f"  X_scvi correlation (reloaded vs stored): {corr:.4f}")
    if corr < 0.95:
        print(f"  [WARN] Low correlation ({corr:.3f}). Verify SCVI_MODEL_DIR is correct.")

# %%
print("\nBuilding scANVI from reloaded scVI...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category = UNLABELED,
    labels_key         = NEW_LABELS_KEY,
)

label_cats  = sorted([c for c in adata.obs[NEW_LABELS_KEY].cat.categories
                       if c != UNLABELED])
n_labeled   = int((adata.obs[NEW_LABELS_KEY] != UNLABELED).sum())
n_unlabeled = int((adata.obs[NEW_LABELS_KEY] == UNLABELED).sum())
print(f"  labels_key       : {NEW_LABELS_KEY}")
print(f"  labeled          : {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"  unknown          : {n_unlabeled:,}")
print(f"  supervision classes ({len(label_cats)}): {label_cats}")

# %%
print(f"\nTraining scANVI "
      f"(max_epochs={SCANVI_EPOCHS}, n_samples_per_label={n_samples_per_label})...")
lvae.train(
    max_epochs              = SCANVI_EPOCHS,
    batch_size              = BATCH_SIZE,
    train_size              = 0.9,
    early_stopping          = True,
    early_stopping_patience = 20,
    n_samples_per_label     = n_samples_per_label,
)
print("scANVI training complete")

# %%
adata.obsm["X_scanvi_refined"]       = lvae.get_latent_representation()
adata.obs["scanvi_pred_refined"]      = lvae.predict()
scanvi_proba                          = lvae.predict(soft=True)
adata.obs["scanvi_pred_prob_refined"] = (
    scanvi_proba.max(axis=1).values.astype(np.float32)
)

concordance = (
    adata.obs["scanvi_pred_refined"].astype(str) ==
    adata.obs[NEW_LABELS_KEY].astype(str)
).mean()
print(f"Overall concordance: {concordance:.3f}")
print(f"Confidence:\n{adata.obs['scanvi_pred_prob_refined'].describe()}")

low_conf = (adata.obs["scanvi_pred_prob_refined"] < 0.5).sum()
print(f"Low-confidence (<0.5): {low_conf:,} ({low_conf/adata.n_obs*100:.1f}%)")

print("\nPer-class concordance:")
for cls in label_cats:
    mask = adata.obs[NEW_LABELS_KEY].astype(str) == cls
    if mask.sum() == 0:
        continue
    acc = (adata.obs.loc[mask, "scanvi_pred_refined"].astype(str) == cls).mean()
    print(f"  {cls:<40}: {acc:.3f}  (n={mask.sum():,})")

# %%
try:
    fig_loss, ax_loss = plt.subplots(figsize=(8, 4))
    train_hist = lvae.history.get("train_loss_epoch", None)
    if train_hist is not None:
        ax_loss.plot(train_hist.values.flatten(), label="train ELBO")
    ax_loss.set_xlabel("Epoch"); ax_loss.set_ylabel("ELBO"); ax_loss.legend()
    ax_loss.set_title("scANVI Retrain -- Training Loss")
    fig_loss.savefig(fig_dir / "scanvi_retrain_training_loss.pdf",
                     bbox_inches="tight", dpi=DPI)
    plt.close(fig_loss)
    print("Saved: scanvi_retrain_training_loss.pdf")
except Exception as e:
    print(f"  Training loss plot failed: {e}")

scanvi_model_dir_new = str(output_dir / "tnk_scanvi_ref_model_retrain")
lvae.save(scanvi_model_dir_new, overwrite=True)
pd.Series(adata.var_names.tolist()).to_csv(
    Path(scanvi_model_dir_new) / "var_names.csv", index=False, header=False)
print(f"scANVI model saved: {scanvi_model_dir_new}/")

del vae, lvae
gc.collect()
if torch.cuda.is_available():
    torch.cuda.empty_cache()
    print("GPU memory released")

# %% [markdown]
# ## 10. Post-Retrain Visualization

# %%
print("Computing neighbors + UMAP on X_scanvi_refined...")
sc.pp.neighbors(adata, use_rep="X_scanvi_refined", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3, spread=1.0, key_added="X_umap_refined")

# %%
sc.settings.vector_friendly = True
fig, axes = plt.subplots(2, 3, figsize=(24, 14))

sc.pl.embedding(adata, basis="umap_refined", color=BACKUP_LABELS_KEY,
                title=f"Original labels ({BACKUP_LABELS_KEY})",
                ax=axes[0, 0], show=False,
                legend_loc="right margin", legend_fontsize=6)

sc.pl.embedding(adata, basis="umap_refined", color=NEW_LABELS_KEY,
                title=f"New labels ({NEW_LABELS_KEY})",
                ax=axes[0, 1], show=False,
                legend_loc="right margin", legend_fontsize=6)

sc.pl.embedding(adata, basis="umap_refined", color="scanvi_pred_refined",
                title="scANVI Predictions",
                ax=axes[0, 2], show=False,
                legend_loc="right margin", legend_fontsize=6)

sc.pl.embedding(adata, basis="umap_refined", color="scanvi_pred_prob_refined",
                title="Prediction Confidence",
                ax=axes[1, 0], show=False, color_map="RdYlGn", vmin=0, vmax=1)

sc.pl.embedding(adata, basis="umap_refined", color=BATCH_KEY,
                title=f"Batch ({BATCH_KEY})",
                ax=axes[1, 1], show=False,
                legend_loc="right margin", legend_fontsize=5)

concordance_col = (
    adata.obs["scanvi_pred_refined"].astype(str) ==
    adata.obs[NEW_LABELS_KEY].astype(str)
).map({True: "correct", False: "mismatch"}).astype("category")
adata.obs["_pred_concordance"] = concordance_col
sc.pl.embedding(adata, basis="umap_refined", color="_pred_concordance",
                title="Prediction Concordance",
                ax=axes[1, 2], show=False,
                palette={"correct": "#2ecc71", "mismatch": "#e74c3c"})
adata.obs.drop(columns=["_pred_concordance"], inplace=True)

plt.suptitle("Post-Retrain Overview -- scANVI Refined", fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(fig_dir / "umap_post_retrain_overview.pdf", bbox_inches="tight", dpi=DPI)
plt.close(fig)
print("Saved: umap_post_retrain_overview.pdf")

# %%
MARKER_GENES = [
    "CD3D", "CD4", "CD8A", "GNLY", "NKG7", "FCGR3A",
    "CCR7", "SELL", "TCF7",
    "GZMB", "PRF1", "GZMK",
    "ITGAE", "CD69", "CXCR6", "RGS1",
    "KLRG1", "CX3CR1",
    "PDCD1", "HAVCR2", "FOXP3",
    "CXCL13", "TOX",
    "TRDC", "KLRB1", "ZBTB16", "KIT", "MKI67",
]
avail_raw    = set(adata.raw.var_names) if adata.raw is not None else set(adata.var_names)
marker_valid = [g for g in MARKER_GENES if g in avail_raw]

sc.settings.vector_friendly = True
ncols = 5
nrows = int(np.ceil(len(marker_valid) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
axes_flat  = np.array(axes).ravel()
for i, gene in enumerate(marker_valid):
    sc.pl.embedding(adata, basis="umap_refined", color=gene,
                    ax=axes_flat[i], show=False,
                    use_raw=(adata.raw is not None),
                    color_map="viridis", vmin=0)
    axes_flat[i].set_title(gene, fontsize=10, fontweight="bold")
for i in range(len(marker_valid), len(axes_flat)):
    axes_flat[i].set_visible(False)
plt.suptitle("T/NK Key Markers -- Refined UMAP (use_raw=True)", fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig(fig_dir / "featureplot_refined_umap.pdf", bbox_inches="tight", dpi=DPI)
plt.close(fig)
print("Saved: featureplot_refined_umap.pdf")

# %% [markdown]
# ## 11. uns Metadata + Save

# %%
import anndata as _anndata
_anndata.settings.allow_write_nullable_strings = True

adata.uns["retrain_params"] = {
    "script_version"          : "v1.2",
    "base_input_h5ad"         : INPUT_H5AD,
    "scvi_model_dir"          : SCVI_MODEL_DIR,
    "scvi_retrained"          : False,
    "scanvi_model_dir_new"    : scanvi_model_dir_new,
    "split_key"               : SPLIT_KEY,
    "leiden_resolution"       : LEIDEN_RESOLUTION,
    "subcluster_rep"          : SUBCLUSTER_REP,
    "old_labels_key"          : OLD_LABELS_KEY,
    "backup_labels_key"       : BACKUP_LABELS_KEY,
    "new_labels_key"          : NEW_LABELS_KEY,
    "unlabeled_category"      : UNLABELED,
    "scanvi_epochs"           : SCANVI_EPOCHS,
    "n_samples_per_label"     : n_samples_per_label,
    "batch_key"               : BATCH_KEY,
    "scanvi_label_categories" : sorted(
        [c for c in adata.obs[NEW_LABELS_KEY].cat.categories if c != UNLABELED]
    ),
    "subcluster_summary"      : dict(subcluster_summary),
    "combined_annotation_map" : COMBINED_ANNOTATION_MAP,
}

cat_cols = [NEW_LABELS_KEY, BACKUP_LABELS_KEY, "scanvi_pred_refined", BATCH_KEY]
if OLD_LABELS_KEY in adata.obs.columns:
    cat_cols.append(OLD_LABELS_KEY)
for col in cat_cols:
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

for attr in ("obs", "var"):
    df = getattr(adata, attr)
    if "_index" in df.columns:
        setattr(adata, attr, df.rename(columns={"_index": "orig_index"}))

# %%
out_h5ad = output_dir / "adata_tnk_scanvi_ref_retrain_v1_2.h5ad"
print(f"Writing {out_h5ad.name}...")
adata.write_h5ad(out_h5ad, compression="gzip", compression_opts=9)
print(f"Saved: {out_h5ad}")

# %% [markdown]
# ## 12. Final Summary

# %%
print("="*70)
print("PIPELINE COMPLETE (v1.2)")
print("="*70)
print(f"\nOutput h5ad      : {out_h5ad}")
print(f"scVI model       : {SCVI_MODEL_DIR}/  (UNCHANGED)")
print(f"scANVI model     : {scanvi_model_dir_new}/")
print(f"Subcluster files : {marker_dir}/")
print(f"Figures          : {fig_dir}/")
print(f"\nKey obs columns:")
print(f"  '{BACKUP_LABELS_KEY}' : original labels (backup)")
print(f"  '{NEW_LABELS_KEY}'    : new refined labels")
print(f"  'scanvi_pred_refined' : scANVI predictions")
print(f"  'X_scanvi_refined'    : refined latent (obsm)")
print(f"\nLabel distribution:")
print(adata.obs[NEW_LABELS_KEY].value_counts())
print(f"\nscArches query checklist:")
print(f"  labels_key     = '{NEW_LABELS_KEY}'")
print(f"  scvi_dir       = '{SCVI_MODEL_DIR}'")
print(f"  scanvi_dir     = '{scanvi_model_dir_new}'")
print(f"  label_space    = {sorted([c for c in adata.obs[NEW_LABELS_KEY].cat.categories if c != UNLABELED])}")


