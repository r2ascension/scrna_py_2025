# %% [markdown]
# # Epithelial Cell BBKNN Integration Pipeline - Linear Execution Version
# 
# Complete BBKNN batch correction workflow with:
# 1. Complete BBKNN batch correction workflow with post-hoc parameter validation
# 2. Force retention of specified markers (avoid removal in HVG/filtering)
# 3. Differential analysis/plotting defaults to reading from .raw full gene matrix (use_raw=True)
# 4. Multi-resolution Leiden clustering; UMAP uses BBKNN-constructed neighbor graph
# 5. Detailed and reproducible visualization and reporting
# 
# Author: Clinical-Bioinformatics Team  
# Version: v1.3 - Linear Execution

# %% [markdown]
# ## 🚀 Quick Start Guide
# 
# ### Option 1: Run Complete Pipeline (from beginning)
# Run all cells sequentially from top to bottom
# 
# ### Option 2: Start from Harmony Section
# 1. Find the section: **"Alternative: Harmony Batch Correction"**
# 2. Uncomment and run the checkpoint cell to load BBKNN results
# 3. Continue from there
# 
# ### Option 3: Start from cNMF Section ⭐
# 1. Find the section: **"cNMF (consensus Non-negative Matrix Factorization) Analysis"**
# 2. Uncomment ALL lines in the checkpoint cell (remove all `#` at the start)
# 3. Run that cell to load previous results
# 4. Continue with cNMF analysis
# 
# ---
# 
# ### 📍 Checkpoint Locations:
# - **Section 1-8**: BBKNN Integration Pipeline (baseline)
# - **🔄 Checkpoint 1**: Before Harmony section
# - **Section 9**: Harmony Batch Correction
# - **🔄 Checkpoint 2**: Before cNMF section ⭐ **← Most commonly used**
# - **Section 10**: cNMF Analysis (with/without batch correction)
# 
# ---

# %% [markdown]
# ## Import Libraries

# %%
import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import time
from tqdm import tqdm
import pickle
import scanpy as sc
from bbknn import bbknn as bbknn_func
from datetime import datetime

warnings.filterwarnings('ignore')

print(f"scanpy version: {sc.__version__}")
print("bbknn: installed")

# %% [markdown]
# ## Configuration Parameters

# %%
# ---------- Input/Output ----------
INPUT_H5AD_PATH = "/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/adata_Epithelial_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/output_optimized"
OVERWRITE_EXISTING = True

# ---------- BBKNN Integration Parameters ----------
BATCH_KEY = "dataset"
BBKNN_NEIGHBORS_WITHIN_BATCH = 3
BBKNN_N_PCS = 50

# ---------- High Variable Genes (HVG) ----------
USE_HVG = True
N_TOP_GENES = 3000
HVG_FLAVOR = "seurat_v3"

# ---------- Gene Filtering ----------
MIN_CELLS_PER_GENE = 3

# ---------- Raw Counts Source ----------
RAW_COUNTS_SOURCE = "auto"

# ---------- Normalization ----------
NORMALIZE_TOTAL = True
TARGET_SUM = 1e4
LOG_TRANSFORM = True
SCALE_DATA = True
MAX_VALUE = 10

# ---------- PCA ----------
N_PCS = 50

# ---------- Dimensionality Reduction/Visualization ----------
RUN_UMAP = True
UMAP_MIN_DIST = 0.5

# ---------- Clustering ----------
RUN_CLUSTERING = True
LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
DEFAULT_RESOLUTION = 2.0

# ---------- Marker Gene Analysis ----------
RUN_FIND_MARKERS = True
MARKER_MIN_PCT = 0.25
MARKER_LOGFC_THRESHOLD = 0.25
TOP_N_MARKERS = 10

# ---------- Performance/Caching ----------
USE_CACHE = True
CHUNK_SIZE = 100

# ---------- Visualization ----------
GENERATE_DOTPLOT = True
GENERATE_HEATMAP = True
GENERATE_FACET_PLOTS = True
DPI = 300
FIGURE_FORMAT = "png"

VERBOSE = True

print("Configuration loaded successfully")

# %%
# ---------- Epithelial Cell Marker Genes ----------
MARKER_EPITHELIAL = [
    # Basal
    'TP63','KRT5','KRT14',
    # Secretory
    'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1',
    # Ciliated
    'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24',
    # Differentiation
    'KRT19','NOTCH3','KRT16','KRT23',
    # Goblet
    'MUC5AC','SPDEF','LYPD2','ITLN1',
    # Neuroendocrine
    'ASCL1','GRP','KRT8',
    # Tuft / Ionocytes
    'POU2F3','ASCL2','CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB','PDE1C',
    # AT1
    'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
    # AT2
    'SFTPC','LAMP3','MFSD2A','C8orf4','C11orf96','SFTPB','SFTA2',
    # Mesenchymal
    'VIM','SOX9','MYH11','ACTA2','MYLK',
    # Secretory/Immune
    'DMBT1','RNASE1','MUC5B','LYZ','LTF','PIP','CCL28',
    # Ciliogenesis
    'DEUP1','FOXN4','CDC20B','CCNO',
    # Proliferation
    'MKI67','TOP2A','TK1','CENPW'
]

# Key markers for cell typing/visualization
KEY_MARKERS = {
    'Basal': ['TP63','KRT5'],
    'Secretory': ['SCGB1A1','SCGB3A2'],
    'Ciliated': ['FOXJ1','RSPH1'],
    'Goblet': ['MUC5AC','SPDEF'],
    'Tuft': ['POU2F3','ASCL2'],
    'Ionocyte': ['FOXI1','CFTR'],
    'AT1': ['AGER','RTKN2'],
    'AT2': ['SFTPC','SFTPB'],
    'Proliferating': ['MKI67','TOP2A']
}

# Force include markers
FORCE_INCLUDE_MARKERS = sorted({
    *MARKER_EPITHELIAL,
    *[g for genes in KEY_MARKERS.values() for g in genes],
})

print(f"Total markers defined: {len(MARKER_EPITHELIAL)}")
print(f"Force include markers: {len(FORCE_INCLUDE_MARKERS)}")

# %% [markdown]
# ## Initialize Output Directory and Timing

# %%
print("\n" + "="*70)
print("Epithelial Cell BBKNN Integration Pipeline (Optimized v1.3)")
print("="*70)

start_time = time.time()

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"\nOutput directory: {output_dir}")

# %% [markdown]
# ## Step 1: Load Data

# %%
print("\n" + "="*70)
print("Step 1: Load Data")
print("="*70)

print(f"\nReading file: {INPUT_H5AD_PATH}")
if not Path(INPUT_H5AD_PATH).exists():
    raise FileNotFoundError(f"File not found: {INPUT_H5AD_PATH}")

adata = sc.read_h5ad(INPUT_H5AD_PATH)
print(f"Data loaded successfully")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")

# %%
print(f"\nAvailable metadata columns:")
for col in adata.obs.columns:
    n_unique = adata.obs[col].nunique()
    print(f"   - {col}: {n_unique} unique values")

if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch column '{BATCH_KEY}' not found in adata.obs")

# %%
print(f"\nBatch distribution (key: {BATCH_KEY}):")
batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
for batch, count in batch_counts.items():
    pct = count / adata.n_obs * 100
    print(f"   {batch}: {count:,} cells ({pct:.1f}%)")

# %% [markdown]
# ## Step 2: Check Marker Genes

# %%
print("\n" + "="*70)
print("Step 2: Check Marker Genes")
print("="*70)

raw_names = set(adata.raw.var_names) if adata.raw is not None else set()
var_names = set(adata.var_names)

def present(g):
    return (g in var_names) or (g in raw_names)

available_markers = [g for g in MARKER_EPITHELIAL if present(g)]
missing_markers   = [g for g in MARKER_EPITHELIAL if not present(g)]

print(f"\nMarker gene availability:")
print(f"   Total: {len(MARKER_EPITHELIAL)}")
print(f"   Available: {len(available_markers)} ({len(available_markers)/len(MARKER_EPITHELIAL)*100:.1f}%)")
print(f"   Missing: {len(missing_markers)} ({len(missing_markers)/len(MARKER_EPITHELIAL)*100:.1f}%)")

if missing_markers:
    print(f"\n   Missing genes (sample): {', '.join(missing_markers[:10])}")
    if len(missing_markers) > 10:
        print(f"   ... and {len(missing_markers)-10} more")

print(f"\nKey cell type markers:")
for celltype, genes in KEY_MARKERS.items():
    available = [g for g in genes if present(g)]
    print(f"   {celltype}: {len(available)}/{len(genes)} available")

# %% [markdown]
# ## Step 3: Preprocessing - Save Raw Counts to Layers

# %%
print("\n" + "="*70)
print(f"Step 3: Select HVGs (n={N_TOP_GENES}) and set layers")
print("-"*70)

# 1) Ensure raw counts in layers['counts']
if 'counts' in adata.layers:
    counts = adata.layers['counts']
    print("   Using existing layers['counts'] as raw counts")
elif getattr(adata, "raw", None) is not None:
    counts = adata.raw.X
    print("   Extracting raw counts from .raw -> layers['counts']")
else:
    counts = adata.X
    print("   No .raw detected, using current X as raw counts -> layers['counts']")

if not issparse(counts):
    counts = csr_matrix(counts)
adata.layers['counts'] = counts

print(f"   Raw counts saved to layers['counts']")

# %% [markdown]
# ## Step 3.1: HVG Selection

# %%
# 2) HVG selection (with robust fallback)
adata.X = adata.layers['counts'].copy()

try:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        flavor="seurat_v3",
        batch_key=BATCH_KEY
    )
    print(f"   HVG(seurat_v3, batch_key={BATCH_KEY}) successful")
except Exception as e:
    print(f"   ⚠️  HVG(seurat_v3, batch_key={BATCH_KEY}) failed: {repr(e)}")
    print("      → Fallback to seurat_v3 without batch_key")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        flavor="seurat_v3",
        batch_key=None
    )

hv_mask = adata.var['highly_variable'].to_numpy(dtype=bool)

# Force retain markers
force_mask = adata.var_names.isin(FORCE_INCLUDE_MARKERS)

keep_mask = hv_mask | force_mask
n_force = int(force_mask.sum())
print(f"   Selected HVGs: {int(hv_mask.sum())}; Force retained: {n_force}; Total: {int(keep_mask.sum())}")

# %% [markdown]
# ## Step 3.2: Normalize Full Matrix

# %%
# 3) Normalize on full matrix to get adata_full (X=log1p)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()

print(f"   Full matrix normalized and log-transformed")
print(f"   layers['log1p'] created")

# %% [markdown]
# ## Step 3.3: Create HVG Subset

# %%
# 4) Create HVG subset object
adata_full = adata.copy()
adata_hvg = adata[:, keep_mask].copy()
adata_hvg.layers['counts'] = adata.layers['counts'][:, keep_mask].copy()

print(f"   Created HVG subset with {adata_hvg.n_vars} genes")

# %%
# Make HVG subset X = log1p of subset counts
adata_hvg.X = adata_hvg.layers['counts'].copy()
sc.pp.normalize_total(adata_hvg, target_sum=1e4)
sc.pp.log1p(adata_hvg)
adata_hvg.layers['log1p'] = adata_hvg.X.copy()

print(f"   HVG subset normalized and log-transformed")

# %%
# 5) Record basic metadata
adata_hvg.uns['hvg_n_top_genes'] = int(N_TOP_GENES)
adata_hvg.uns['force_include_markers'] = np.array(FORCE_INCLUDE_MARKERS, dtype=str)
adata_hvg.uns['keep_mask_sum'] = int(keep_mask.sum())

# 6) Set full gene data as HVG subset's .raw
adata_hvg.raw = adata_full
print(f"   Set adata_hvg.raw = adata_full (contains {adata_full.n_vars} genes)")
print("   ✓ Complete: raw counts in layers['counts'], log1p in layers['log1p']")
print("="*70)

# %% [markdown]
# ## Step 3.5: PCA

# %%
print(f"\nRunning PCA (n_comps={N_PCS})...")
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')

var_ratio = adata_hvg.uns['pca']['variance_ratio']
cumsum_var = np.cumsum(var_ratio)
if len(cumsum_var) >= 50:
    print(f"   PC1-10 explained variance: {cumsum_var[9]:.2%}")
    print(f"   PC1-20 explained variance: {cumsum_var[19]:.2%}")
    print(f"   PC1-50 explained variance: {cumsum_var[49]:.2%}")

# %% [markdown]
# ## Step 4: BBKNN Batch Integration

# %%
print("\n" + "="*70)
print("Step 4: BBKNN Batch Integration")
print("="*70)

print("\nBBKNN parameters:")
print(f"   batch_key: {BATCH_KEY}")
print(f"   neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
print(f"   n_pcs: {BBKNN_N_PCS}")

print("\nRunning BBKNN...")
bbknn_start = time.time()
bbknn_func(
    adata_hvg,
    batch_key=BATCH_KEY,
    neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
    n_pcs=BBKNN_N_PCS,
    copy=False
)
bbknn_time = time.time() - bbknn_start
print(f"\n   BBKNN complete, time elapsed: {bbknn_time:.1f} seconds")

# %%
# BBKNN sanity check
try:
    n_batches = adata_hvg.obs[BATCH_KEY].nunique()
    expected = BBKNN_NEIGHBORS_WITHIN_BATCH * n_batches
    neigh = adata_hvg.uns.get('neighbors', {}).get('params', {})
    actual = neigh.get('n_neighbors', None)
    print(
        f"   BBKNN check → n_batches={n_batches}, "
        f"neighbors_within_batch={BBKNN_NEIGHBORS_WITHIN_BATCH} → "
        f"expected total neighbors≈{expected}, actual {actual}"
    )
except Exception as e:
    print(f"   ⚠️ BBKNN post-hoc validation failed: {e}")

# %% [markdown]
# ## Step 5: UMAP and Multi-resolution Clustering

# %%
print("\n" + "="*70)
print("Step 5: UMAP and Multi-resolution Clustering")
print("="*70)

# %%
if RUN_UMAP:
    print(f"\nRunning UMAP (min_dist={UMAP_MIN_DIST}, using BBKNN neighbor graph)...")
    sc.tl.umap(adata_hvg, min_dist=UMAP_MIN_DIST)
    print("   UMAP complete")

# %%
if RUN_CLUSTERING:
    print(f"\nRunning multi-resolution Leiden clustering...")
    print(f"   Resolution list: {LEIDEN_RESOLUTIONS}")
    for res in LEIDEN_RESOLUTIONS:
        cluster_key = f'leiden_bbknn_res{res}'
        try:
            sc.tl.leiden(adata_hvg, resolution=res, key_added=cluster_key, flavor='igraph',
                         n_iterations=2, directed=False)
        except TypeError:
            sc.tl.leiden(adata_hvg, resolution=res, key_added=cluster_key, n_iterations=2)
        n_clusters = adata_hvg.obs[cluster_key].nunique()
        print(f"   Resolution {res}: {n_clusters} clusters")
    
    default_key = f'leiden_bbknn_res{DEFAULT_RESOLUTION}'
    adata_hvg.obs['leiden_bbknn'] = adata_hvg.obs[default_key]
    print(f"\n   Default clustering: {default_key}")

# %% [markdown]
# ## Step 6: Compute Cluster Differential Genes

# %%
if RUN_FIND_MARKERS:
    print("\n" + "="*70)
    print("Step 6: Compute Cluster Differential Genes (Optimized)")
    print("="*70)
    
    cache_file = "epithelial_markers_cache.pkl"
    cache_path = output_dir / f".cache_{cache_file}"
    
    if USE_CACHE and cache_path.exists():
        print("\n   Loading from cache...")
        with open(cache_path, 'rb') as f:
            all_markers = pickle.load(f)
        print(f"   Total markers: {len(all_markers)}")
    else:
        markers_file = output_dir / "cluster_markers.csv"
        
        print(f"\nRunning differential analysis (Wilcoxon)...")
        print(f"   min_pct: {MARKER_MIN_PCT}")
        print(f"   logfc_threshold: {MARKER_LOGFC_THRESHOLD}")
        
        print("\n   Computing differentials...")
        sc.tl.rank_genes_groups(
            adata_hvg,
            groupby='leiden_bbknn',
            method='wilcoxon',
            key_added='rank_genes_groups',
            use_raw=True,
            layer=None
        )
        print("   ✓ Differential analysis complete")

# %%
if RUN_FIND_MARKERS and not (USE_CACHE and cache_path.exists()):
    print("\n   Extracting and filtering differential results...")
    
    # Helper function for computing pct expressed
    def compute_pct_expressed_vectorized(adata, cluster_key, cluster, genes):
        cluster_mask = (adata.obs[cluster_key] == cluster).values
        n_in = int(cluster_mask.sum()); n_out = int((~cluster_mask).sum())
        if n_in == 0 or n_out == 0 or len(genes) == 0:
            return [], []

        genes_in_var = [g for g in genes if g in adata.var_names]
        rem = [g for g in genes if g not in adata.var_names]
        X_in = X_out = None

        if genes_in_var:
            Xin = adata[cluster_mask, genes_in_var].X
            Xout = adata[~cluster_mask, genes_in_var].X
            Xin = Xin.toarray() if hasattr(Xin, 'toarray') else np.asarray(Xin)
            Xout = Xout.toarray() if hasattr(Xout, 'toarray') else np.asarray(Xout)
            X_in, X_out = Xin, Xout

        if rem and (adata.raw is not None):
            raw_names = np.array(adata.raw.var_names)
            raw_map = {g:i for i,g in enumerate(raw_names)}
            rem_exist = [g for g in rem if g in raw_map]
            if rem_exist:
                ridx = [raw_map[g] for g in rem_exist]
                R = adata.raw.X
                Rin = R[cluster_mask][:, ridx]
                Rout = R[~cluster_mask][:, ridx]
                Rin = Rin.toarray() if hasattr(Rin, 'toarray') else np.asarray(Rin)
                Rout = Rout.toarray() if hasattr(Rout, 'toarray') else np.asarray(Rout)
                if X_in is None:
                    X_in, X_out = Rin, Rout
                    genes_in_var = rem_exist
                else:
                    X_in = np.concatenate([X_in, Rin], axis=1)
                    X_out = np.concatenate([X_out, Rout], axis=1)
                    genes_in_var = genes_in_var + rem_exist

        if X_in is None:
            return [], []
        pct_in  = (X_in  > 0).sum(axis=0) / n_in
        pct_out = (X_out > 0).sum(axis=0) / n_out
        return pct_in.tolist(), pct_out.tolist()
    
    clusters = adata_hvg.obs['leiden_bbknn'].unique()
    markers_list = []

    for cluster in tqdm(clusters, desc="Processing clusters"):
        df = sc.get.rank_genes_groups_df(adata_hvg, group=cluster)
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges","pvals_adj"])
        df = df[(df["logfoldchanges"] > MARKER_LOGFC_THRESHOLD) & (df["pvals_adj"] < 0.05)].copy()
        if df.empty:
            continue

        genes = df["names"].tolist()
        pct_in, pct_out = compute_pct_expressed_vectorized(adata_hvg, 'leiden_bbknn', cluster, genes)
        if len(pct_in) != len(genes):
            n = len(genes)
            pct_in  = (pct_in  + [np.nan]*(n-len(pct_in)))[:n]
            pct_out = (pct_out + [np.nan]*(n-len(pct_out)))[:n]

        df["cluster"] = cluster
        df["pct_in_cluster"] = pct_in
        df["pct_out_cluster"] = pct_out
        df = df[df["pct_in_cluster"] > MARKER_MIN_PCT]
        if df.empty:
            continue

        if TOP_N_MARKERS is not None and TOP_N_MARKERS > 0:
            df = df.sort_values("pvals_adj", ascending=True).head(TOP_N_MARKERS)
        markers_list.append(df)

    if len(markers_list) == 0:
        print("\n   ⚠️ No differential genes meet current thresholds, returning empty results.")
        empty_cols = ["names","scores","logfoldchanges","pvals","pvals_adj","cluster","pct_in_cluster","pct_out_cluster"]
        all_markers = pd.DataFrame(columns=empty_cols)
    else:
        all_markers = pd.concat(markers_list, ignore_index=True)
        print(f"\n   ✓ Total {len(all_markers)} markers obtained, covering {len(clusters)} clusters")
    
    # Save to cache and CSV
    if USE_CACHE:
        with open(cache_path, 'wb') as f:
            pickle.dump(all_markers, f)
        print(f"   Checkpoint saved: {cache_path}")
    
    all_markers.to_csv(markers_file, index=False)
    print(f"   ✓ Saved: {markers_file}")

# %% [markdown]
# ## Step 7: Generate Visualizations

# %%
if RUN_UMAP:
    print("\n" + "="*70)
    print("Step 7: Generate Visualizations")
    print("="*70)

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    sc.settings.figdir = fig_dir
    
    print(f"\nFigures will be saved to: {fig_dir}")

# %%
if RUN_UMAP:
    # 1) Overview plot
    print("\n   Generating integration quality overview...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata_hvg, color=BATCH_KEY, ax=axes[0], show=False, title='Batch')
    sc.pl.umap(adata_hvg, color='leiden_bbknn', ax=axes[1], show=False,
               title='Clusters', legend_loc='on data', legend_fontsize=8)
    if 'cell_type' in adata_hvg.obs.columns:
        sc.pl.umap(adata_hvg, color='cell_type', ax=axes[2], show=False, title='Cell Type')
    else:
        axes[2].axis('off')
    plt.tight_layout()
    plt.savefig(fig_dir / f'umap_overview.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Overview UMAP saved")

# %%
if RUN_UMAP and RUN_CLUSTERING and len(LEIDEN_RESOLUTIONS) > 1:
    # 2) Multi-resolution clustering
    print("\n   Generating multi-resolution clustering plots...")
    n_res = len(LEIDEN_RESOLUTIONS)
    n_cols = 3
    n_rows = (n_res + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 6*n_rows))
    axes = axes.flatten() if n_res > 1 else [axes]
    for i, res in enumerate(LEIDEN_RESOLUTIONS):
        cluster_key = f'leiden_bbknn_res{res}'
        sc.pl.umap(
            adata_hvg, color=cluster_key, ax=axes[i], show=False,
            title=f'Resolution {res}', legend_loc='on data', legend_fontsize=8
        )
    for i in range(n_res, len(axes)):
        axes[i].axis('off')
    plt.tight_layout()
    plt.savefig(fig_dir / f'umap_resolutions.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Multi-resolution UMAP saved")

# %%
if RUN_UMAP:
    # 3) Key marker expression UMAPs
    print("\n   Generating marker expression UMAPs...")
    for celltype, genes in KEY_MARKERS.items():
        available = [g for g in genes if (g in (adata_hvg.raw.var_names if adata_hvg.raw is not None else adata_hvg.var_names))]
        if available:
            sc.pl.umap(
                adata_hvg,
                color=available,
                use_raw=True,
                cmap='Reds',
                ncols=len(available),
                vmax='p99',
                save=f'_{celltype}_markers.{FIGURE_FORMAT}'
            )
            print(f"   ✓ {celltype} markers UMAP saved")

# %%
if RUN_UMAP and GENERATE_DOTPLOT and len(available_markers) > 0:
    # 4) Dotplot
    print("\n   Generating marker Dotplot...")
    try:
        markers_to_plot = available_markers[:50] if len(available_markers) > 50 else available_markers
        sc.pl.dotplot(
            adata_hvg,
            var_names=markers_to_plot,
            groupby='leiden_bbknn',
            standard_scale='var',
            use_raw=True,
            save=f'_epithelial_markers.{FIGURE_FORMAT}',
            figsize=(max(20, len(markers_to_plot)*0.4), 10)
        )
        print(f"   ✓ Dotplot saved (showing {len(markers_to_plot)} genes)")
    except Exception as e:
        print(f"   ⚠️ Dotplot generation failed: {e}")

# %%
if RUN_UMAP and GENERATE_HEATMAP and len(available_markers) > 0:
    # 5) Heatmap
    print("\n   Generating expression heatmap...")
    try:
        cluster_expr = pd.DataFrame()
        for gene in available_markers:
            if (adata_hvg.raw is not None) and (gene in adata_hvg.raw.var_names):
                Xg = adata_hvg.raw[:, gene].X
            elif gene in adata_hvg.var_names:
                Xg = adata_hvg[:, gene].X
            else:
                continue
            Xg = Xg.toarray().ravel() if hasattr(Xg, 'toarray') else np.asarray(Xg).ravel()
            cluster_expr[gene] = Xg
        cluster_expr['leiden_bbknn'] = adata_hvg.obs['leiden_bbknn'].values
        cluster_mean_expr = cluster_expr.groupby('leiden_bbknn').mean()
        plt.figure(figsize=(max(20, len(available_markers)*0.3), 10))
        sns.heatmap(
            cluster_mean_expr.T,
            cmap='RdYlBu_r',
            center=0,
            robust=True,
            yticklabels=True,
            xticklabels=True,
            cbar_kws={'label': 'Mean Expression'}
        )
        plt.title('Epithelial Marker Expression Across Clusters')
        plt.xlabel('Cluster')
        plt.ylabel('Marker Genes')
        plt.tight_layout()
        plt.savefig(fig_dir / f'heatmap_marker_expression.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Heatmap saved")
    except Exception as e:
        print(f"   ⚠️ Heatmap generation failed: {e}")

# %%
if RUN_UMAP and GENERATE_FACET_PLOTS and BATCH_KEY in adata_hvg.obs.columns:
    # 6) Batch facets
    print("\n   Generating batch-wise UMAP...")
    try:
        batches = sorted(adata_hvg.obs[BATCH_KEY].unique())
        n_batches = len(batches)
        n_cols = min(3, n_batches)
        n_rows = (n_batches + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
        axes = [axes] if n_batches == 1 else axes.flatten()
        
        for i, batch in enumerate(batches):
            adata_batch = adata_hvg[adata_hvg.obs[BATCH_KEY] == batch].copy()
            
            sc.pl.umap(
                adata_batch,
                color='leiden_bbknn',
                ax=axes[i],
                show=False,
                title=f'Batch: {batch}',
                legend_loc='right margin',
                legend_fontsize=6
            )
        
        for i in range(n_batches, len(axes)):
            axes[i].axis('off')
        
        plt.tight_layout()
        plt.savefig(fig_dir / f'umap_by_batch.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Batch-wise UMAP saved")
    except Exception as e:
        print(f"   ⚠️ Batch-wise UMAP generation failed: {e}")

# %%
if RUN_UMAP:
    print(f"\n   All figures saved to: {fig_dir}")

# %% [markdown]
# ## Backfill Results to Full Dataset

# %%
print("\nBackfilling results to full dataset...")

# Transfer results from HVG subset to full dataset
for key in ['X_pca','X_umap']:
    if key in adata_hvg.obsm:
        adata_full.obsm[key] = adata_hvg.obsm[key]

for col in adata_hvg.obs.columns:
    if col.startswith('leiden_bbknn'):
        adata_full.obs[col] = adata_hvg.obs[col]

if 'neighbors' in adata_hvg.uns:
    adata_full.uns['neighbors'] = adata_hvg.uns['neighbors']

if 'connectivities' in adata_hvg.obsp:
    adata_full.obsp['connectivities'] = adata_hvg.obsp['connectivities']

if 'distances' in adata_hvg.obsp:
    adata_full.obsp['distances'] = adata_hvg.obsp['distances']

print("   ✓ Results backfilled")

# %% [markdown]
# ## Save Integrated Data

# %%
print("\nSaving integrated AnnData...")
final_path = output_dir / "epithelial_bbknn_integrated.h5ad"
print("   Using gzip compression...")
adata_full.write_h5ad(final_path, compression='gzip', compression_opts=9)
file_size = final_path.stat().st_size / (1024**3)
print(f"   Saved: {final_path} ({file_size:.2f} GB)")

# %% [markdown]
# ## Step 8: Generate Summary Report

# %%
print("\n" + "="*70)
print("Step 8: Generate Summary Report")
print("="*70)

total_time = time.time() - start_time

report = []
report.append("="*70)
report.append("Epithelial Cell Analysis - BBKNN Integration Summary")
report.append("="*70)
report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report.append(f"Processing time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
report.append("")
report.append("[ Data Overview ]")
report.append(f"  Cells: {adata_full.n_obs:,}")
report.append(f"  Genes: {adata_full.n_vars:,}")
report.append("")

if BATCH_KEY in adata_full.obs.columns:
    report.append(f"[ Batch Distribution (key: {BATCH_KEY}) ]")
    batch_counts = adata_full.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata_full.n_obs * 100
        report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
    report.append("")

report.append("[ Epithelial Marker Genes ]")
report.append(f"  Total markers: {len(MARKER_EPITHELIAL)}")
report.append(f"  Available: {len(available_markers)} ({len(available_markers)/len(MARKER_EPITHELIAL)*100:.1f}%)")
report.append("")

if RUN_CLUSTERING:
    report.append("[ Multi-resolution Clustering Results ]")
    for res in LEIDEN_RESOLUTIONS:
        key = f'leiden_bbknn_res{res}'
        if key in adata_full.obs.columns:
            n_clusters = adata_full.obs[key].nunique()
            default_marker = " (default)" if res == DEFAULT_RESOLUTION else ""
            report.append(f"  Resolution {res}: {n_clusters} clusters{default_marker}")
    report.append("")

report.append("[ BBKNN Configuration ]")
report.append(f"  neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
report.append(f"  n_pcs: {BBKNN_N_PCS}")
report.append(f"  batch_key: {BATCH_KEY}")
report.append(f"  High variable genes: {USE_HVG} (n={N_TOP_GENES if USE_HVG else 'N/A'})")
report.append("")

report.append("[ Output Files ]")
report.append(f"  - Data: {final_path}")
if RUN_UMAP:
    report.append(f"  - Figures: {output_dir / 'figures/'}*.{FIGURE_FORMAT}")
if RUN_FIND_MARKERS:
    report.append(f"  - Markers: {output_dir / 'cluster_markers.csv'}")
report.append("")

report.append("="*70)
report.append("Analysis completed successfully")
report.append("="*70)

report_text = '\n'.join(report)
report_path = output_dir / "analysis_summary.txt"
with open(report_path, 'w') as f:
    f.write(report_text)

print(f"\nReport saved to: {report_path}")
print("\n" + report_text)

# %% [markdown]
# ## Final Summary

# %%
print("\n" + "="*70)
print("✓ All analyses completed successfully")
print("="*70)
print(f"\nOutput directory: {output_dir}")
print(f"Data: {final_path} ({file_size:.2f} GB)")
if RUN_UMAP:
    print(f"Figures: {output_dir / 'figures/'}*.{FIGURE_FORMAT}")
if RUN_FIND_MARKERS:
    print(f"Markers: {output_dir / 'cluster_markers.csv'}")
print(f"\nTotal time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")

print(f"\n📊 Key Information:")
print(f"   Cells: {adata_full.n_obs:,}")
print(f"   Available markers: {len(available_markers)}/{len(MARKER_EPITHELIAL)}")
if RUN_CLUSTERING:
    print(f"   Default clusters: {adata_full.obs['leiden_bbknn'].nunique()}")

print(f"\n🔧 BBKNN Parameters:")
print(f"   neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
print(f"   n_pcs: {BBKNN_N_PCS}")
print(f"   HVG: {N_TOP_GENES}")
print()

# %% [markdown]
# ---
# 
# # Alternative: Harmony Batch Correction
# 
# This section demonstrates an alternative batch correction approach using Harmony.  
# Unlike BBKNN which builds batch-balanced neighbor graphs, Harmony corrects the principal components directly.

# %% [markdown]
# ## 🔄 Checkpoint: Load Previous Results (Optional)
# 
# Run this cell if you want to start from here without running previous steps

# %%
# Uncomment and run this cell to start from Harmony section

# import scanpy as sc
# import pandas as pd
# import numpy as np
# from pathlib import Path
# import matplotlib.pyplot as plt
# import seaborn as sns
# import time

# # Load saved results from BBKNN analysis
# output_dir = Path("/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/output_optimized")
# final_path = output_dir / "epithelial_bbknn_integrated.h5ad"

# print("Loading BBKNN results...")
# adata_full = sc.read_h5ad(final_path)
# adata_hvg = adata_full[:, adata_full.var_names.isin(adata_full.uns.get('force_include_markers', []))].copy()
# adata_hvg.raw = adata_full

# # Restore configuration
# BATCH_KEY = "dataset"
# BBKNN_N_PCS = 50
# UMAP_MIN_DIST = 0.5
# DEFAULT_RESOLUTION = 2.0
# LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
# DPI = 300
# FIGURE_FORMAT = "png"
# RUN_CLUSTERING = True
# GENERATE_FACET_PLOTS = True

# print(f"✓ Loaded: {adata_full.n_obs} cells, {adata_full.n_vars} genes")
# print(f"   BBKNN clusters: {adata_full.obs['leiden_bbknn'].nunique()}")
# start_time = time.time()

# %% [markdown]
# ## Harmony Configuration

# %%
# Harmony-specific parameters
HARMONY_BATCH_KEY = BATCH_KEY  # Use same batch key as BBKNN
HARMONY_MAX_ITER = 20          # Maximum number of iterations

# Theta: Diversity clustering penalty parameter (ridge regression penalty)
# - Range: 0 to infinity
# - Default: 2.0
# - Higher values (e.g., 3-5): More aggressive batch correction, risk of over-correction
# - Lower values (e.g., 0.5-1): Gentler correction, may retain batch effects
# - Can be a list for multiple batch keys: e.g., [2.0, 1.5] for ['dataset', 'sex']
HARMONY_THETA = 1.0

# Lambda: Ridge regression penalty parameter
# - Range: 0 to infinity  
# - Default: 1.0
# - Higher values (e.g., 2-5): Stronger regularization, more conservative correction
# - Lower values (e.g., 0.1-0.5): Less regularization, more flexible correction
# - Can be a list for multiple batch keys
HARMONY_LAMBDA = 5.0

# Sigma: Width of soft kmeans clusters
# - Range: 0 to 1
# - Default: 0.1
# - Typically doesn't need adjustment
HARMONY_SIGMA = 0.1

# Number of clusters for Harmony
# - Default: None (auto-determined)
# - Can set manually: e.g., 100
HARMONY_NCLUST = None

# Tau: Protection against overclustering
# - Range: 0 to infinity
# - Default: 0
HARMONY_TAU = 0

# Whether to run Harmony analysis
RUN_HARMONY = True

print(f"Harmony configuration:")
print(f"   batch_key: {HARMONY_BATCH_KEY}")
print(f"   max_iter: {HARMONY_MAX_ITER}")
print(f"\n   Key parameters:")
print(f"   theta (diversity penalty): {HARMONY_THETA}")
print(f"   lambda (ridge penalty): {HARMONY_LAMBDA}")
print(f"   sigma (cluster width): {HARMONY_SIGMA}")
print(f"   nclust (num clusters): {HARMONY_NCLUST}")
print(f"   tau (overclustering protection): {HARMONY_TAU}")

# %% [markdown]
# ### 🧪 Parameter Testing Examples
# 
# **Example 1: Gentle correction (preserve more biological variation)**
# ```python
# HARMONY_THETA = 1.0
# HARMONY_LAMBDA = 2.0
# ```
# 
# **Example 2: Aggressive correction (remove strong batch effects)**
# ```python
# HARMONY_THETA = 4.0
# HARMONY_LAMBDA = 0.5
# ```
# 
# **Example 3: Balanced correction (recommended starting point)**
# ```python
# HARMONY_THETA = 2.0
# HARMONY_LAMBDA = 1.0
# ```
# 
# **Example 4: Multiple batch variables with different strengths**
# ```python
# HARMONY_BATCH_KEY = ['dataset', 'sex']  # Correct for both
# HARMONY_THETA = [2.0, 1.0]               # Stronger for dataset, gentler for sex
# HARMONY_LAMBDA = [1.0, 2.0]              # More conservative for sex
# ```
# 
# **Parameter Selection Strategy:**
# 1. Start with default (theta=2, lambda=1)
# 2. Check batch mixing plots
# 3. If under-corrected: increase theta OR decrease lambda
# 4. If over-corrected: decrease theta OR increase lambda
# 5. Compare clustering consistency across batches

# %% [markdown]
# ## Install and Import Harmony (if needed)

# %%
if RUN_HARMONY:
    try:
        import harmonypy as hm
        print(f"harmonypy version: {hm.__version__}")
    except ImportError:
        print("harmonypy not found. Installing...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy"])
        import harmonypy as hm
        print(f"harmonypy installed successfully: {hm.__version__}")

# %% [markdown]
# ## Prepare Data for Harmony
# 
# Harmony works on the PCA space, so we'll use the same preprocessed data

# %%
if RUN_HARMONY:
    print("\n" + "="*70)
    print("Harmony Batch Correction Pipeline")
    print("="*70)
    
    # Create a copy for Harmony analysis
    adata_harmony = adata_hvg.copy()
    print(f"\nCreated Harmony working copy:")
    print(f"   Cells: {adata_harmony.n_obs:,}")
    print(f"   Genes: {adata_harmony.n_vars:,}")
    print(f"   PCA components: {adata_harmony.obsm['X_pca'].shape[1]}")

# %% [markdown]
# ## Run Harmony on PCA

# %%
if RUN_HARMONY:
    print("\nRunning Harmony batch correction...")
    print(f"   Correcting for: {HARMONY_BATCH_KEY}")
    print(f"   theta={HARMONY_THETA}, lambda={HARMONY_LAMBDA}")
    
    harmony_start = time.time()
    
    # Prepare parameters
    harmony_params = {
        'key': HARMONY_BATCH_KEY,
        'basis': 'X_pca',
        'adjusted_basis': 'X_pca_harmony',
        'max_iter_harmony': HARMONY_MAX_ITER,
    }
    
    # Add theta (can be single value or list)
    if isinstance(HARMONY_THETA, (list, tuple)):
        harmony_params['theta'] = HARMONY_THETA
    else:
        harmony_params['theta'] = [HARMONY_THETA] if isinstance(HARMONY_BATCH_KEY, list) else HARMONY_THETA
    
    # Add lambda (can be single value or list)
    if isinstance(HARMONY_LAMBDA, (list, tuple)):
        harmony_params['lamb'] = HARMONY_LAMBDA  # Note: parameter name is 'lamb' not 'lambda'
    else:
        harmony_params['lamb'] = [HARMONY_LAMBDA] if isinstance(HARMONY_BATCH_KEY, list) else HARMONY_LAMBDA
    
    # Add sigma
    harmony_params['sigma'] = HARMONY_SIGMA
    
    # Add nclust if specified
    if HARMONY_NCLUST is not None:
        harmony_params['nclust'] = HARMONY_NCLUST
    
    # Add tau
    harmony_params['tau'] = HARMONY_TAU
    
    # Run Harmony with all parameters
    sc.external.pp.harmony_integrate(
        adata_harmony,
        **harmony_params
    )
    
    harmony_time = time.time() - harmony_start
    print(f"\n   Harmony complete, time elapsed: {harmony_time:.1f} seconds")
    print(f"   Corrected PCs stored in: .obsm['X_pca_harmony']")
    print(f"\n   Final parameters used:")
    print(f"   - theta: {harmony_params.get('theta')}")
    print(f"   - lambda: {harmony_params.get('lamb')}")
    print(f"   - sigma: {harmony_params.get('sigma')}")
    print(f"   - tau: {harmony_params.get('tau')}")

# %% [markdown]
# ## Compute UMAP on Harmony-corrected PCs

# %%
if RUN_HARMONY:
    print("\nComputing neighbor graph on Harmony-corrected PCs...")
    
    # Compute neighbors using Harmony-corrected PCs
    sc.pp.neighbors(
        adata_harmony,
        n_neighbors=15,
        n_pcs=BBKNN_N_PCS,
        use_rep='X_pca_harmony'
    )
    
    print("   Neighbor graph computed")
    
    # Compute UMAP
    print(f"\nRunning UMAP on Harmony PCs (min_dist={UMAP_MIN_DIST})...")
    sc.tl.umap(adata_harmony, min_dist=UMAP_MIN_DIST)
    print("   UMAP complete")

# %% [markdown]
# ## Clustering on Harmony-corrected Data

# %%
if RUN_HARMONY and RUN_CLUSTERING:
    print("\nRunning multi-resolution Leiden clustering on Harmony data...")
    print(f"   Resolution list: {LEIDEN_RESOLUTIONS}")
    
    for res in LEIDEN_RESOLUTIONS:
        cluster_key = f'leiden_harmony_res{res}'
        try:
            sc.tl.leiden(
                adata_harmony,
                resolution=res,
                key_added=cluster_key,
                flavor='igraph',
                n_iterations=2,
                directed=False
            )
        except TypeError:
            sc.tl.leiden(
                adata_harmony,
                resolution=res,
                key_added=cluster_key,
                n_iterations=2
            )
        n_clusters = adata_harmony.obs[cluster_key].nunique()
        print(f"   Resolution {res}: {n_clusters} clusters")
    
    # Set default clustering
    default_key = f'leiden_harmony_res{DEFAULT_RESOLUTION}'
    adata_harmony.obs['leiden_harmony'] = adata_harmony.obs[default_key]
    print(f"\n   Default clustering: {default_key}")

# %% [markdown]
# ## Compare BBKNN vs Harmony Integration Quality

# %%
if RUN_HARMONY:
    print("\n" + "="*70)
    print("Comparison: BBKNN vs Harmony")
    print("="*70)
    
    # Compare cluster numbers
    bbknn_clusters = adata_hvg.obs['leiden_bbknn'].nunique()
    harmony_clusters = adata_harmony.obs['leiden_harmony'].nunique()
    
    print(f"\nClustering results (resolution={DEFAULT_RESOLUTION}):")
    print(f"   BBKNN:   {bbknn_clusters} clusters")
    print(f"   Harmony: {harmony_clusters} clusters")
    
    # Compare batch distribution
    print(f"\nBatch distribution across clusters:")
    print(f"\nBBKNN:")
    bbknn_batch_dist = pd.crosstab(
        adata_hvg.obs['leiden_bbknn'],
        adata_hvg.obs[BATCH_KEY],
        normalize='index'
    )
    print(bbknn_batch_dist.round(3))
    
    print(f"\nHarmony:")
    harmony_batch_dist = pd.crosstab(
        adata_harmony.obs['leiden_harmony'],
        adata_harmony.obs[BATCH_KEY],
        normalize='index'
    )
    print(harmony_batch_dist.round(3))

# %% [markdown]
# ## Visualize Harmony Results

# %%
if RUN_HARMONY:
    print("\nGenerating Harmony visualization...")
    
    harmony_fig_dir = output_dir / "figures_harmony"
    harmony_fig_dir.mkdir(exist_ok=True)
    
    # Overview comparison plot
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Row 1: BBKNN results
    sc.pl.umap(adata_hvg, color=BATCH_KEY, ax=axes[0,0], show=False, title='BBKNN - Batch')
    sc.pl.umap(adata_hvg, color='leiden_bbknn', ax=axes[0,1], show=False,
               title='BBKNN - Clusters', legend_loc='on data', legend_fontsize=6)
    if 'cell_type' in adata_hvg.obs.columns:
        sc.pl.umap(adata_hvg, color='cell_type', ax=axes[0,2], show=False, title='BBKNN - Cell Type')
    else:
        axes[0,2].axis('off')
    
    # Row 2: Harmony results
    sc.pl.umap(adata_harmony, color=BATCH_KEY, ax=axes[1,0], show=False, title='Harmony - Batch')
    sc.pl.umap(adata_harmony, color='leiden_harmony', ax=axes[1,1], show=False,
               title='Harmony - Clusters', legend_loc='on data', legend_fontsize=6)
    if 'cell_type' in adata_harmony.obs.columns:
        sc.pl.umap(adata_harmony, color='cell_type', ax=axes[1,2], show=False, title='Harmony - Cell Type')
    else:
        axes[1,2].axis('off')
    
    plt.tight_layout()
    plt.savefig(harmony_fig_dir / f'comparison_bbknn_vs_harmony.{FIGURE_FORMAT}',
                dpi=DPI, bbox_inches='tight')
    plt.close()
    
    print(f"   ✓ Comparison plot saved to {harmony_fig_dir}")

# %% [markdown]
# ## Harmony Multi-resolution Clustering Visualization

# %%
if RUN_HARMONY and RUN_CLUSTERING and len(LEIDEN_RESOLUTIONS) > 1:
    print("\nGenerating Harmony multi-resolution clustering plots...")
    
    n_res = len(LEIDEN_RESOLUTIONS)
    n_cols = 3
    n_rows = (n_res + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 6*n_rows))
    axes = axes.flatten() if n_res > 1 else [axes]
    
    for i, res in enumerate(LEIDEN_RESOLUTIONS):
        cluster_key = f'leiden_harmony_res{res}'
        sc.pl.umap(
            adata_harmony,
            color=cluster_key,
            ax=axes[i],
            show=False,
            title=f'Harmony - Resolution {res}',
            legend_loc='on data',
            legend_fontsize=8
        )
    
    for i in range(n_res, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    plt.savefig(harmony_fig_dir / f'harmony_resolutions.{FIGURE_FORMAT}',
                dpi=DPI, bbox_inches='tight')
    plt.close()
    
    print(f"   ✓ Multi-resolution plot saved")

# %% [markdown]
# ## Harmony Marker Expression Visualization

# %%
if RUN_HARMONY:
    print("\nGenerating Harmony marker expression UMAPs...")
    
    # Set figure directory for scanpy
    sc.settings.figdir = harmony_fig_dir
    
    for celltype, genes in KEY_MARKERS.items():
        available = [g for g in genes if (g in (adata_harmony.raw.var_names if adata_harmony.raw is not None else adata_harmony.var_names))]
        if available:
            sc.pl.umap(
                adata_harmony,
                color=available,
                use_raw=True,
                cmap='Reds',
                ncols=len(available),
                vmax='p99',
                save=f'_harmony_{celltype}_markers.{FIGURE_FORMAT}'
            )
            print(f"   ✓ {celltype} markers UMAP saved")
    
    # Reset figdir
    sc.settings.figdir = output_dir / "figures"

# %% [markdown]
# ## Harmony Batch-wise Facet Plot

# %%
if RUN_HARMONY and GENERATE_FACET_PLOTS and BATCH_KEY in adata_harmony.obs.columns:
    print("\nGenerating Harmony batch-wise UMAP...")
    
    try:
        batches = sorted(adata_harmony.obs[BATCH_KEY].unique())
        n_batches = len(batches)
        n_cols = min(3, n_batches)
        n_rows = (n_batches + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
        axes = [axes] if n_batches == 1 else axes.flatten()
        
        for i, batch in enumerate(batches):
            adata_batch = adata_harmony[adata_harmony.obs[BATCH_KEY] == batch].copy()
            
            sc.pl.umap(
                adata_batch,
                color='leiden_harmony',
                ax=axes[i],
                show=False,
                title=f'Harmony - Batch: {batch}',
                legend_loc='right margin',
                legend_fontsize=6
            )
        
        for i in range(n_batches, len(axes)):
            axes[i].axis('off')
        
        plt.tight_layout()
        plt.savefig(harmony_fig_dir / f'harmony_umap_by_batch.{FIGURE_FORMAT}',
                    dpi=DPI, bbox_inches='tight')
        plt.close()
        
        print(f"   ✓ Batch-wise UMAP saved")
    except Exception as e:
        print(f"   ⚠️ Batch-wise UMAP generation failed: {e}")

# %% [markdown]
# ## Save Harmony Results

# %%
if RUN_HARMONY:
    print("\nBackfilling Harmony results to full dataset...")
    
    # Transfer Harmony results to full dataset
    adata_full.obsm['X_pca_harmony'] = adata_harmony.obsm['X_pca_harmony']
    adata_full.obsm['X_umap_harmony'] = adata_harmony.obsm['X_umap']
    
    for col in adata_harmony.obs.columns:
        if col.startswith('leiden_harmony'):
            adata_full.obs[col] = adata_harmony.obs[col]
    
    print("   ✓ Harmony results backfilled")
    
    # Save Harmony-integrated data
    print("\nSaving Harmony-integrated AnnData...")
    harmony_path = output_dir / "epithelial_harmony_integrated.h5ad"
    print("   Using gzip compression...")
    adata_full.write_h5ad(harmony_path, compression='gzip', compression_opts=9)
    harmony_file_size = harmony_path.stat().st_size / (1024**3)
    print(f"   Saved: {harmony_path} ({harmony_file_size:.2f} GB)")

# %% [markdown]
# ## Harmony Summary Report

# %%
if RUN_HARMONY:
    print("\n" + "="*70)
    print("Harmony Integration Summary")
    print("="*70)
    
    harmony_report = []
    harmony_report.append("="*70)
    harmony_report.append("Epithelial Cell Analysis - Harmony Integration Summary")
    harmony_report.append("="*70)
    harmony_report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    harmony_report.append(f"Harmony processing time: {harmony_time:.1f} seconds")
    harmony_report.append("")
    
    harmony_report.append("[ Harmony Configuration ]")
    harmony_report.append(f"  batch_key: {HARMONY_BATCH_KEY}")
    harmony_report.append(f"  max_iter: {HARMONY_MAX_ITER}")
    harmony_report.append(f"  n_pcs: {BBKNN_N_PCS}")
    harmony_report.append("")
    harmony_report.append("[ Harmony Parameters ]")
    harmony_report.append(f"  theta (diversity penalty): {HARMONY_THETA}")
    harmony_report.append(f"  lambda (ridge penalty): {HARMONY_LAMBDA}")
    harmony_report.append(f"  sigma (cluster width): {HARMONY_SIGMA}")
    harmony_report.append(f"  tau (overclustering protection): {HARMONY_TAU}")
    if HARMONY_NCLUST is not None:
        harmony_report.append(f"  nclust (number of clusters): {HARMONY_NCLUST}")
    harmony_report.append("")
    
    harmony_report.append("[ Clustering Results ]")
    harmony_report.append(f"  Default resolution: {DEFAULT_RESOLUTION}")
    harmony_report.append(f"  Harmony clusters: {adata_full.obs['leiden_harmony'].nunique()}")
    harmony_report.append(f"  BBKNN clusters: {adata_full.obs['leiden_bbknn'].nunique()}")
    harmony_report.append("")
    
    harmony_report.append("[ Parameter Guidelines ]")
    harmony_report.append("  theta: Higher = more aggressive correction (try 2-5)")
    harmony_report.append("  lambda: Higher = more conservative (try 0.5-5)")
    harmony_report.append("  sigma: Typically 0.1, adjust if convergence issues")
    harmony_report.append("")
    
    harmony_report.append("[ Output Files ]")
    harmony_report.append(f"  - Data: {harmony_path}")
    harmony_report.append(f"  - Figures: {harmony_fig_dir}/*.{FIGURE_FORMAT}")
    harmony_report.append("")
    
    harmony_report.append("="*70)
    harmony_report.append("Harmony analysis completed successfully")
    harmony_report.append("="*70)
    
    harmony_report_text = '\n'.join(harmony_report)
    harmony_report_path = output_dir / "harmony_summary.txt"
    with open(harmony_report_path, 'w') as f:
        f.write(harmony_report_text)
    
    print(f"\nHarmony report saved to: {harmony_report_path}")
    print("\n" + harmony_report_text)

# %% [markdown]
# ## Final Summary: BBKNN + Harmony

# %%
if RUN_HARMONY:
    print("\n" + "="*70)
    print("✓ Complete Pipeline (BBKNN + Harmony) Finished")
    print("="*70)
    
    total_pipeline_time = time.time() - start_time
    
    print(f"\n📊 Final Results Summary:")
    print(f"   Total cells: {adata_full.n_obs:,}")
    print(f"   Total genes: {adata_full.n_vars:,}")
    print(f"   Available markers: {len(available_markers)}/{len(MARKER_EPITHELIAL)}")
    print(f"\n   BBKNN clusters (res={DEFAULT_RESOLUTION}): {adata_full.obs['leiden_bbknn'].nunique()}")
    print(f"   Harmony clusters (res={DEFAULT_RESOLUTION}): {adata_full.obs['leiden_harmony'].nunique()}")
    
    print(f"\n⏱️  Processing Times:")
    print(f"   BBKNN integration: {bbknn_time:.1f} seconds")
    print(f"   Harmony integration: {harmony_time:.1f} seconds")
    print(f"   Total pipeline: {total_pipeline_time:.1f} seconds ({total_pipeline_time/60:.1f} minutes)")
    
    print(f"\n💾 Output Files:")
    print(f"   BBKNN results: {final_path}")
    print(f"   Harmony results: {harmony_path}")
    print(f"   BBKNN figures: {output_dir / 'figures/'}")
    print(f"   Harmony figures: {harmony_fig_dir}")
    
    print(f"\n📝 Recommendation:")
    print(f"   Compare the batch mixing and clustering quality between BBKNN and Harmony")
    print(f"   Check batch distribution plots and UMAP visualizations")
    print(f"   Choose the method that best preserves biological signal while removing batch effects")
    print()

# %% [markdown]
# ---
# 
# # cNMF (consensus Non-negative Matrix Factorization) Analysis
# 
# This section runs cNMF to identify gene expression programs (GEPs) in epithelial cells.  
# We'll run cNMF on both:
# 1. **Uncorrected data** - to see batch effects
# 2. **Harmony-corrected data** - using the two-stage approach (corrected for usage, uncorrected for spectra)

# %% [markdown]
# ## 🔄 Checkpoint: Load Previous Results (Start from cNMF)
# 
# **Run this cell if you want to start directly from cNMF analysis without running previous sections**

# %%
# ============================================================================
# CHECKPOINT: Start from cNMF Analysis
# Uncomment ALL lines below to load previous results and skip BBKNN/Harmony
# ============================================================================

# import scanpy as sc
# import pandas as pd
# import numpy as np
# from pathlib import Path
# import matplotlib.pyplot as plt
# import seaborn as sns
# import time
# from datetime import datetime
# import warnings
# warnings.filterwarnings('ignore')

# print("="*70)
# print("Loading Previous Results for cNMF Analysis")
# print("="*70)

# # Configuration
# output_dir = Path("/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/output_optimized")
# BATCH_KEY = "dataset"
# BBKNN_N_PCS = 50
# N_TOP_GENES = 3000
# UMAP_MIN_DIST = 0.5
# DEFAULT_RESOLUTION = 2.0
# LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
# DPI = 300
# FIGURE_FORMAT = "png"
# RUN_CLUSTERING = True
# RUN_UMAP = True
# GENERATE_FACET_PLOTS = True

# # Marker genes
# KEY_MARKERS = {
#     'Basal': ['TP63','KRT5'],
#     'Secretory': ['SCGB1A1','SCGB3A2'],
#     'Ciliated': ['FOXJ1','RSPH1'],
#     'Goblet': ['MUC5AC','SPDEF'],
#     'Tuft': ['POU2F3','ASCL2'],
#     'Ionocyte': ['FOXI1','CFTR'],
#     'AT1': ['AGER','RTKN2'],
#     'AT2': ['SFTPC','SFTPB'],
#     'Proliferating': ['MKI67','TOP2A']
# }

# # Load integrated data
# print("\nLoading integrated data...")
# integrated_file = output_dir / "epithelial_bbknn_integrated.h5ad"
# if not integrated_file.exists():
#     # Try harmony version
#     integrated_file = output_dir / "epithelial_harmony_integrated.h5ad"

# if integrated_file.exists():
#     adata_full = sc.read_h5ad(integrated_file)
#     print(f"   ✓ Loaded: {adata_full.n_obs} cells × {adata_full.n_vars} genes")
#     print(f"   Batches: {adata_full.obs[BATCH_KEY].nunique()}")
#     
#     # Check what's available
#     if 'leiden_bbknn' in adata_full.obs.columns:
#         print(f"   BBKNN clusters: {adata_full.obs['leiden_bbknn'].nunique()}")
#     if 'leiden_harmony' in adata_full.obs.columns:
#         print(f"   Harmony clusters: {adata_full.obs['leiden_harmony'].nunique()}")
#         RUN_HARMONY = True
#     else:
#         RUN_HARMONY = False
#     
#     # Create HVG subset if needed
#     if 'highly_variable' in adata_full.var.columns:
#         hvg_mask = adata_full.var['highly_variable'].values
#         adata_hvg = adata_full[:, hvg_mask].copy()
#         adata_hvg.raw = adata_full
#         print(f"   HVG subset: {adata_hvg.n_vars} genes")
#     else:
#         adata_hvg = adata_full.copy()
#     
#     # Create harmony copy if exists
#     if 'X_pca_harmony' in adata_full.obsm.keys():
#         adata_harmony = adata_hvg.copy()
#         print(f"   ✓ Harmony results available")
#     
#     start_time = time.time()
#     print(f"\n✓ Ready to run cNMF analysis")
#     print("="*70)
# else:
#     print(f"\n✗ Error: Could not find integrated data file")
#     print(f"   Expected: {integrated_file}")
#     print(f"   Please run BBKNN or Harmony integration first")

# %% [markdown]
# ## Install and Import cNMF

# %%
# Install cnmf if needed
try:
    from cnmf import cNMF
    print("cNMF already installed")
except ImportError:
    print("Installing cNMF...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "cnmf"])
    from cnmf import cNMF
    print("cNMF installed successfully")

# %% [markdown]
# ## cNMF Configuration

# %%
# cNMF parameters
RUN_CNMF = True
CNMF_OUTPUT_DIR = output_dir / "cnmf_results"
CNMF_OUTPUT_DIR.mkdir(exist_ok=True)

# cNMF settings
CNMF_COMPONENTS = [10, 15, 20, 25, 30]  # Test multiple K values
CNMF_N_ITER = 100  # Number of NMF replicates
CNMF_SEED = 14
CNMF_NUM_HVG = 2000  # Number of high variance genes
CNMF_DENSITY_THRESHOLD = 0.1  # For consensus clustering

print(f"cNMF Configuration:")
print(f"   Output directory: {CNMF_OUTPUT_DIR}")
print(f"   Components to test: {CNMF_COMPONENTS}")
print(f"   Iterations per K: {CNMF_N_ITER}")
print(f"   HVG count: {CNMF_NUM_HVG}")
print(f"   Density threshold: {CNMF_DENSITY_THRESHOLD}")

# %% [markdown]
# ## Prepare Data for cNMF
# 
# We need to save the data in formats cNMF can read

# %%
if RUN_CNMF:
    print("\n" + "="*70)
    print("Preparing data for cNMF")
    print("="*70)
    
    # Use the original full data (before any batch correction)
    # We need raw counts for cNMF
    cnmf_data_path = CNMF_OUTPUT_DIR / "epithelial_raw_counts.h5ad"
    
    # Create a clean AnnData with raw counts
    adata_cnmf = adata_full.copy()
    
    # Make sure we have raw counts in X
    if 'counts' in adata_cnmf.layers:
        adata_cnmf.X = adata_cnmf.layers['counts'].copy()
        print(f"   Using counts from layers['counts']")
    else:
        print(f"   Warning: No counts layer found, using current X")
    
    # Save for cNMF
    adata_cnmf.write_h5ad(cnmf_data_path)
    print(f"   Saved raw counts to: {cnmf_data_path}")
    print(f"   Shape: {adata_cnmf.n_obs} cells × {adata_cnmf.n_vars} genes")

# %% [markdown]
# ## Run cNMF without Batch Correction (Baseline)
# 
# First, let's run cNMF on uncorrected data to see the batch effects

# %%
if RUN_CNMF:
    print("\n" + "="*70)
    print("cNMF Analysis - WITHOUT Batch Correction")
    print("="*70)
    
    cnmf_nocorr = cNMF(output_dir=str(CNMF_OUTPUT_DIR), name='Epithelial_NoBatchCorrection')
    
    print("\nStep 1: Prepare...")
    cnmf_nocorr.prepare(
        counts_fn=str(cnmf_data_path),
        components=CNMF_COMPONENTS,
        n_iter=CNMF_N_ITER,
        seed=CNMF_SEED,
        num_highvar_genes=CNMF_NUM_HVG
    )
    print("   ✓ Preparation complete")

# %%
if RUN_CNMF:
    print("\nStep 2: Factorize (this may take a while)...")
    print(f"   Running {CNMF_N_ITER} iterations for each K in {CNMF_COMPONENTS}")
    
    factorize_start = time.time()
    cnmf_nocorr.factorize(worker_i=0, total_workers=1)
    factorize_time = time.time() - factorize_start
    
    print(f"   ✓ Factorization complete ({factorize_time:.1f} seconds)")

# %%
if RUN_CNMF:
    print("\nStep 3: Combine results...")
    cnmf_nocorr.combine()
    print("   ✓ Results combined")

# %%
if RUN_CNMF:
    print("\nStep 4: Consensus clustering for each K...")
    
    for k in CNMF_COMPONENTS:
        print(f"\n   Processing K={k}...")
        cnmf_nocorr.consensus(
            k=k,
            density_threshold=CNMF_DENSITY_THRESHOLD,
            show_clustering=True,
            close_clustergram_fig=True
        )
        print(f"   ✓ K={k} consensus complete")

# %% [markdown]
# ## Select Optimal K and Load Results (No Correction)

# %%
if RUN_CNMF:
    # Choose K based on stability metrics
    # For this example, we'll use K=15 (you should examine the plots to choose)
    SELECTED_K = 15
    
    print(f"\nLoading cNMF results for K={SELECTED_K} (no batch correction)...")
    
    (usage_nocorr, spectra_scores_nocorr, 
     spectra_tpm_nocorr, top_genes_nocorr) = cnmf_nocorr.load_results(
        K=SELECTED_K,
        density_threshold=CNMF_DENSITY_THRESHOLD
    )
    
    print(f"   Usage matrix shape: {usage_nocorr.shape}")
    print(f"   Spectra shape: {spectra_tpm_nocorr.shape}")
    print(f"\n   First few cells:")
    print(usage_nocorr.head())

# %% [markdown]
# ## Visualize GEP Usage by Batch (No Correction)
# 
# This will show batch effects in the uncorrected data

# %%
if RUN_CNMF:
    print("\nVisualizing GEP usage by batch (uncorrected data)...")
    
    # Prepare data for plotting
    usage_plot_nocorr = usage_nocorr.unstack().reset_index()
    usage_plot_nocorr = pd.merge(
        left=usage_plot_nocorr,
        right=adata_full.obs[[BATCH_KEY, 'leiden_bbknn']],
        left_on='level_1',
        right_index=True
    )
    usage_plot_nocorr.columns = ['GEP', 'cell', 'Usage', 'Batch', 'Cluster']
    
    # Get unique clusters
    clusters = sorted(adata_full.obs['leiden_bbknn'].unique())
    
    # Plot
    fig, axes = plt.subplots(
        len(clusters), 1,
        figsize=(12, len(clusters)*1.5),
        dpi=150,
        gridspec_kw={'hspace': 0.8}
    )
    
    if len(clusters) == 1:
        axes = [axes]
    
    for i, cluster in enumerate(clusters):
        cluster_data = usage_plot_nocorr[usage_plot_nocorr['Cluster'] == cluster]
        sns.boxplot(
            x='GEP', y='Usage', hue='Batch',
            data=cluster_data,
            ax=axes[i],
            fliersize=1
        )
        axes[i].set_title(f'Cluster {cluster}')
        axes[i].set_xlabel('')
        axes[i].set_ylabel('Usage')
        
        if i != (len(clusters)-1):
            axes[i].set_xticklabels([])
        
        if i != 0:
            axes[i].legend().remove()
        else:
            axes[i].legend(bbox_to_anchor=(1, 1), title='Batch')
    
    plt.tight_layout()
    plt.savefig(
        CNMF_OUTPUT_DIR / f'gep_usage_by_batch_nocorr_K{SELECTED_K}.{FIGURE_FORMAT}',
        dpi=DPI, bbox_inches='tight'
    )
    plt.close()
    
    print(f"   ✓ Plot saved")
    print(f"\n   Look for GEPs that are batch-specific (different distributions across batches)")

# %% [markdown]
# ## Prepare Harmony-Corrected Data for cNMF
# 
# Now we'll use the Harmony approach: corrected counts for learning usage, uncorrected TP10K for fitting spectra

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\n" + "="*70)
    print("Preparing Harmony-corrected data for cNMF")
    print("="*70)
    
    # We need to create:
    # 1. Corrected HVG counts (for learning usage)
    # 2. Uncorrected TP10K matrix (for fitting spectra)
    # 3. HVG gene list
    
    print("\nNote: Creating corrected counts from Harmony PCs...")
    print("This is a simplified approach - ideally use harmonypy's full preprocessing pipeline")

# %%
if RUN_CNMF and RUN_HARMONY:
    # Step 1: Create variance-normalized corrected data
    # Use the harmony-corrected PCA to reconstruct expression
    
    # Get HVG subset
    adata_hvg_harmony = adata_harmony.copy()
    
    # For simplicity, we'll use variance-normalized log-transformed data
    # This mimics what Harmony preprocessing does
    adata_corrected_hvg = adata_hvg_harmony.copy()
    
    # Scale the data
    sc.pp.scale(adata_corrected_hvg, max_value=10)
    
    # Save corrected HVG data
    corrected_hvg_path = CNMF_OUTPUT_DIR / "epithelial_harmony_corrected_hvg.h5ad"
    adata_corrected_hvg.write_h5ad(corrected_hvg_path)
    print(f"   Corrected HVG data saved: {corrected_hvg_path}")
    print(f"   Shape: {adata_corrected_hvg.n_obs} cells × {adata_corrected_hvg.n_vars} genes")

# %%
if RUN_CNMF and RUN_HARMONY:
    # Step 2: Create uncorrected TP10K matrix (all genes)
    adata_tp10k = adata_full.copy()
    
    # Get raw counts
    if 'counts' in adata_tp10k.layers:
        adata_tp10k.X = adata_tp10k.layers['counts'].copy()
    
    # Normalize to TP10K (TPM × 10,000)
    sc.pp.normalize_total(adata_tp10k, target_sum=1e4)
    
    # Save TP10K data
    tp10k_path = CNMF_OUTPUT_DIR / "epithelial_tp10k.h5ad"
    adata_tp10k.write_h5ad(tp10k_path)
    print(f"   TP10K data saved: {tp10k_path}")
    print(f"   Shape: {adata_tp10k.n_obs} cells × {adata_tp10k.n_vars} genes")

# %%
if RUN_CNMF and RUN_HARMONY:
    # Step 3: Save HVG list
    hvg_list = adata_corrected_hvg.var_names.tolist()
    hvg_file = CNMF_OUTPUT_DIR / "harmony_hvgs.txt"
    
    with open(hvg_file, 'w') as f:
        for gene in hvg_list:
            f.write(f"{gene}\n")
    
    print(f"   HVG list saved: {hvg_file}")
    print(f"   Number of HVGs: {len(hvg_list)}")

# %% [markdown]
# ## Run cNMF with Harmony Batch Correction
# 
# Using the two-stage approach: corrected data for usage, uncorrected for spectra

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\n" + "="*70)
    print("cNMF Analysis - WITH Harmony Batch Correction")
    print("="*70)
    
    cnmf_harmony = cNMF(output_dir=str(CNMF_OUTPUT_DIR), name='Epithelial_HarmonyBatchCorrected')
    
    print("\nStep 1: Prepare...")
    cnmf_harmony.prepare(
        counts_fn=str(corrected_hvg_path),  # Corrected counts for usage
        tpm_fn=str(tp10k_path),             # Uncorrected TP10K for spectra
        genes_file=str(hvg_file),           # HVG list
        components=CNMF_COMPONENTS,
        n_iter=CNMF_N_ITER,
        seed=CNMF_SEED,
        num_highvar_genes=CNMF_NUM_HVG
    )
    print("   ✓ Preparation complete")

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\nStep 2: Factorize (this may take a while)...")
    
    factorize_harmony_start = time.time()
    cnmf_harmony.factorize(worker_i=0, total_workers=1)
    factorize_harmony_time = time.time() - factorize_harmony_start
    
    print(f"   ✓ Factorization complete ({factorize_harmony_time:.1f} seconds)")

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\nStep 3: Combine results...")
    cnmf_harmony.combine()
    print("   ✓ Results combined")

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\nStep 4: Consensus clustering for each K...")
    
    for k in CNMF_COMPONENTS:
        print(f"\n   Processing K={k}...")
        cnmf_harmony.consensus(
            k=k,
            density_threshold=CNMF_DENSITY_THRESHOLD,
            show_clustering=True,
            close_clustergram_fig=True
        )
        print(f"   ✓ K={k} consensus complete")

# %% [markdown]
# ## Load Harmony-Corrected cNMF Results

# %%
if RUN_CNMF and RUN_HARMONY:
    print(f"\nLoading cNMF results for K={SELECTED_K} (Harmony-corrected)...")
    
    (usage_harmony, spectra_scores_harmony,
     spectra_tpm_harmony, top_genes_harmony) = cnmf_harmony.load_results(
        K=SELECTED_K,
        density_threshold=CNMF_DENSITY_THRESHOLD
    )
    
    print(f"   Usage matrix shape: {usage_harmony.shape}")
    print(f"   Spectra shape: {spectra_tpm_harmony.shape}")
    print(f"\n   First few cells:")
    print(usage_harmony.head())

# %% [markdown]
# ## Visualize GEP Usage by Batch (Harmony-Corrected)
# 
# Compare with uncorrected to see the improvement

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\nVisualizing GEP usage by batch (Harmony-corrected)...")
    
    # Prepare data for plotting
    usage_plot_harmony = usage_harmony.unstack().reset_index()
    usage_plot_harmony = pd.merge(
        left=usage_plot_harmony,
        right=adata_full.obs[[BATCH_KEY, 'leiden_bbknn']],
        left_on='level_1',
        right_index=True
    )
    usage_plot_harmony.columns = ['GEP', 'cell', 'Usage', 'Batch', 'Cluster']
    
    # Plot
    fig, axes = plt.subplots(
        len(clusters), 1,
        figsize=(12, len(clusters)*1.5),
        dpi=150,
        gridspec_kw={'hspace': 0.8}
    )
    
    if len(clusters) == 1:
        axes = [axes]
    
    for i, cluster in enumerate(clusters):
        cluster_data = usage_plot_harmony[usage_plot_harmony['Cluster'] == cluster]
        sns.boxplot(
            x='GEP', y='Usage', hue='Batch',
            data=cluster_data,
            ax=axes[i],
            fliersize=1
        )
        axes[i].set_title(f'Cluster {cluster} (Harmony-corrected)')
        axes[i].set_xlabel('')
        axes[i].set_ylabel('Usage')
        
        if i != (len(clusters)-1):
            axes[i].set_xticklabels([])
        
        if i != 0:
            axes[i].legend().remove()
        else:
            axes[i].legend(bbox_to_anchor=(1, 1), title='Batch')
    
    plt.tight_layout()
    plt.savefig(
        CNMF_OUTPUT_DIR / f'gep_usage_by_batch_harmony_K{SELECTED_K}.{FIGURE_FORMAT}',
        dpi=DPI, bbox_inches='tight'
    )
    plt.close()
    
    print(f"   ✓ Plot saved")
    print(f"\n   GEPs should be more consistent across batches now")

# %% [markdown]
# ## Compare Top Genes for Each GEP

# %%
if RUN_CNMF:
    print("\n" + "="*70)
    print("Top Genes for Each GEP (No Correction)")
    print("="*70)
    
    for gep in range(1, SELECTED_K + 1):
        print(f"\nGEP {gep}:")
        top_genes_gep = top_genes_nocorr[str(gep)][:10]
        print(f"   {', '.join(top_genes_gep)}")

# %%
if RUN_CNMF and RUN_HARMONY:
    print("\n" + "="*70)
    print("Top Genes for Each GEP (Harmony-Corrected)")
    print("="*70)
    
    for gep in range(1, SELECTED_K + 1):
        print(f"\nGEP {gep}:")
        top_genes_gep = top_genes_harmony[str(gep)][:10]
        print(f"   {', '.join(top_genes_gep)}")

# %% [markdown]
# ## Add GEP Usage to AnnData Objects

# %%
if RUN_CNMF:
    print("\nAdding GEP usage to AnnData objects...")
    
    # Add uncorrected GEP usage
    for gep in range(1, SELECTED_K + 1):
        gep_col = f'GEP_{gep}_nocorr'
        adata_full.obs[gep_col] = usage_nocorr[str(gep)].reindex(adata_full.obs_names).values
    
    print(f"   ✓ Added {SELECTED_K} GEP columns (uncorrected) to adata_full")
    
    if RUN_HARMONY:
        # Add Harmony-corrected GEP usage
        for gep in range(1, SELECTED_K + 1):
            gep_col = f'GEP_{gep}_harmony'
            adata_full.obs[gep_col] = usage_harmony[str(gep)].reindex(adata_full.obs_names).values
        
        print(f"   ✓ Added {SELECTED_K} GEP columns (Harmony-corrected) to adata_full")

# %% [markdown]
# ## Visualize GEP Usage on UMAP

# %%
if RUN_CNMF and RUN_UMAP:
    print("\nGenerating GEP usage UMAPs (uncorrected)...")
    
    cnmf_fig_dir = output_dir / "figures_cnmf"
    cnmf_fig_dir.mkdir(exist_ok=True)
    
    # Select a few interesting GEPs to plot
    geps_to_plot = [f'GEP_{i}_nocorr' for i in range(1, min(7, SELECTED_K + 1))]
    
    sc.pl.umap(
        adata_full,
        color=geps_to_plot,
        cmap='Reds',
        ncols=3,
        vmax='p99',
        save=f'_gep_usage_nocorr.{FIGURE_FORMAT}'
    )
    
    print(f"   ✓ UMAP plots saved")

# %%
if RUN_CNMF and RUN_HARMONY and RUN_UMAP:
    print("\nGenerating GEP usage UMAPs (Harmony-corrected)...")
    
    # Select a few interesting GEPs to plot
    geps_to_plot = [f'GEP_{i}_harmony' for i in range(1, min(7, SELECTED_K + 1))]
    
    sc.pl.umap(
        adata_full,
        color=geps_to_plot,
        cmap='Reds',
        ncols=3,
        vmax='p99',
        save=f'_gep_usage_harmony.{FIGURE_FORMAT}'
    )
    
    print(f"   ✓ UMAP plots saved")

# %% [markdown]
# ## Save Final Integrated Data with cNMF Results

# %%
if RUN_CNMF:
    print("\nSaving final integrated data with cNMF results...")
    
    final_integrated_path = output_dir / "epithelial_complete_integrated.h5ad"
    adata_full.write_h5ad(final_integrated_path, compression='gzip', compression_opts=9)
    
    final_size = final_integrated_path.stat().st_size / (1024**3)
    print(f"   Saved: {final_integrated_path} ({final_size:.2f} GB)")
    print(f"\n   This file contains:")
    print(f"   - BBKNN integration results")
    if RUN_HARMONY:
        print(f"   - Harmony integration results")
    print(f"   - cNMF GEP usage scores")
    print(f"   - All clustering results")
    print(f"   - All visualizations")

# %% [markdown]
# ## cNMF Summary Report

# %%
if RUN_CNMF:
    print("\n" + "="*70)
    print("cNMF Analysis Summary")
    print("="*70)
    
    cnmf_report = []
    cnmf_report.append("="*70)
    cnmf_report.append("Epithelial Cell cNMF Analysis Summary")
    cnmf_report.append("="*70)
    cnmf_report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    cnmf_report.append("")
    
    cnmf_report.append("[ cNMF Configuration ]")
    cnmf_report.append(f"  Selected K: {SELECTED_K}")
    cnmf_report.append(f"  Tested K values: {CNMF_COMPONENTS}")
    cnmf_report.append(f"  Iterations per K: {CNMF_N_ITER}")
    cnmf_report.append(f"  High variance genes: {CNMF_NUM_HVG}")
    cnmf_report.append(f"  Density threshold: {CNMF_DENSITY_THRESHOLD}")
    cnmf_report.append("")
    
    cnmf_report.append("[ Analysis Performed ]")
    cnmf_report.append("  1. cNMF without batch correction (baseline)")
    if RUN_HARMONY:
        cnmf_report.append("  2. cNMF with Harmony batch correction")
    cnmf_report.append("")
    
    cnmf_report.append("[ Key Findings ]")
    cnmf_report.append("  - Check boxplots to compare batch effects before/after correction")
    cnmf_report.append("  - GEPs should be more consistent across batches after correction")
    cnmf_report.append("  - Top genes reveal biological programs in epithelial cells")
    cnmf_report.append("")
    
    cnmf_report.append("[ Output Files ]")
    cnmf_report.append(f"  - cNMF results directory: {CNMF_OUTPUT_DIR}")
    cnmf_report.append(f"  - GEP usage plots: {CNMF_OUTPUT_DIR}/gep_usage_*.{FIGURE_FORMAT}")
    cnmf_report.append(f"  - Final integrated data: {final_integrated_path}")
    cnmf_report.append("")
    
    cnmf_report.append("="*70)
    cnmf_report.append("cNMF analysis completed successfully")
    cnmf_report.append("="*70)
    
    cnmf_report_text = '\n'.join(cnmf_report)
    cnmf_report_path = output_dir / "cnmf_summary.txt"
    with open(cnmf_report_path, 'w') as f:
        f.write(cnmf_report_text)
    
    print(f"\ncNMF report saved to: {cnmf_report_path}")
    print("\n" + cnmf_report_text)

# %% [markdown]
# ## Complete Pipeline Summary

# %%
print("\n" + "="*70)
print("🎉 COMPLETE PIPELINE FINISHED 🎉")
print("="*70)

final_pipeline_time = time.time() - start_time

print(f"\n📊 Analysis Summary:")
print(f"   Cells analyzed: {adata_full.n_obs:,}")
print(f"   Genes: {adata_full.n_vars:,}")
print(f"   Batches: {adata_full.obs[BATCH_KEY].nunique()}")

print(f"\n🔬 Methods Applied:")
print(f"   ✓ BBKNN batch integration")
if RUN_HARMONY:
    print(f"   ✓ Harmony batch correction")
if RUN_CNMF:
    print(f"   ✓ cNMF gene expression programs (K={SELECTED_K})")
print(f"   ✓ Multi-resolution clustering")
print(f"   ✓ Differential gene analysis")

print(f"\n⏱️  Total Processing Time: {final_pipeline_time:.1f} seconds ({final_pipeline_time/60:.1f} minutes)")

print(f"\n💾 Key Output Files:")
print(f"   Main results: {final_integrated_path}")
print(f"   BBKNN figures: {output_dir / 'figures/'}")
if RUN_HARMONY:
    print(f"   Harmony figures: {harmony_fig_dir}")
if RUN_CNMF:
    print(f"   cNMF results: {CNMF_OUTPUT_DIR}")

print(f"\n📈 Next Steps:")
print(f"   1. Compare BBKNN vs Harmony batch mixing quality")
if RUN_CNMF:
    print(f"   2. Examine GEP boxplots for batch effect reduction")
    print(f"   3. Interpret biological meaning of top genes in each GEP")
print(f"   4. Validate findings with marker genes")
print(f"   5. Consider downstream analyses (trajectory, cell-cell communication, etc.)")

print(f"\n✅ All analyses completed successfully!")
print()


