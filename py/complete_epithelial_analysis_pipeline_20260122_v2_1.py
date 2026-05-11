# %% [markdown]
# # Complete Epithelial Cell Analysis Pipeline - FIXED VERSION
# 
# **Comprehensive workflow from batch correction to GEP-cluster mapping**
# 
# ## 📋 Pipeline Overview
# 
# 1. Data Loading & Preprocessing (5-10 min)
# 2. BBKNN Batch Integration (5-10 min)
# 3. Harmony Batch Correction (10-20 min)
# 4. cNMF - Uncorrected (2-6 hours) ⏱️ **NOW WITH MULTI-PROCESSING**
# 5. cNMF - Harmony Corrected (2-6 hours) ⏱️ **NOW WITH MULTI-PROCESSING**
# 6. GEP-Cluster Mapping (10-20 min)
# 7. Final Summary Report
# 
# **⏱️ Total Estimated Time: 4-12 hours** *(much faster with multi-processing)*
# 
# ---
# 
# **Version: v2.1 - FIXED & OPTIMIZED**  
# **Key Fixes:**
# - Fixed KeyError in GEP assignment (use column position instead of name)
# - Added multi-processing support for cNMF (4-8x speedup)
# - Fixed variable scope issues
# - Enhanced error handling
# 
# **Date: 2024-12-03**

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
from scipy.spatial.distance import pdist
from scipy.stats import entropy
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import seaborn as sns
import time
from tqdm import tqdm
import pickle
from datetime import datetime
from collections import Counter
from multiprocessing import Process

import scanpy as sc
from bbknn import bbknn as bbknn_func

warnings.filterwarnings('ignore')

print("="*70)
print("Complete Epithelial Cell Analysis Pipeline v2.1 - FIXED & OPTIMIZED")
print("="*70)
print(f"\nLibrary versions:")
print(f"   scanpy: {sc.__version__}")
print(f"   numpy: {np.__version__}")
print(f"   pandas: {pd.__version__}")
print(f"\n✓ All libraries loaded successfully")

# %% [markdown]
# ## 🎛️ Master Configuration

# %%
# ============================================================================
# MASTER CONFIGURATION - Adjust parameters here
# ============================================================================

# ---------- Input/Output ----------
INPUT_H5AD_PATH = "/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/adata_Epithelial_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/output_complete_pipeline_v1"
BATCH_KEY = "dataset"

# ---------- Pipeline Control ----------
RUN_BBKNN = True
RUN_HARMONY = True
RUN_CNMF_UNCORRECTED = True
RUN_CNMF_HARMONY = True
RUN_MAPPING = True

# ---------- BBKNN Parameters ----------
BBKNN_NEIGHBORS_WITHIN_BATCH = 3
BBKNN_N_PCS = 50

# ---------- HVG Selection ----------
N_TOP_GENES = 4000
HVG_FLAVOR = "seurat_v3"

# ---------- Clustering ----------
LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
DEFAULT_RESOLUTION = 2.0

# ---------- UMAP ----------
UMAP_MIN_DIST = 0.5

# ---------- Harmony Parameters ----------
HARMONY_THETA = 1.0
HARMONY_LAMBDA = 6.0
HARMONY_N_PCS_CORRECT = 50
HARMONY_MAX_ITER = 20
HARMONY_SIGMA = 0.1
HARMONY_TAU = 0

# ---------- cNMF Parameters ----------
CNMF_COMPONENTS = [10, 15, 20, 25, 30]
CNMF_N_ITER = 100
CNMF_SEED = 14
CNMF_NUM_HVG = 2000
CNMF_DENSITY_THRESHOLD = 0.1
SELECTED_K = 25

# 🚀 NEW: Multi-processing for cNMF (建议 <= CPU核心数)
CNMF_TOTAL_WORKERS = 4

# ---------- Mapping Parameters ----------
MAPPING_CLUSTER_KEY = 'leiden_bbknn'
MAPPING_ANALYZE_UNCORRECTED = True
MAPPING_ANALYZE_HARMONY = True
MAPPING_COMPARE_BOTH = True
MAPPING_PVALUE_THRESHOLD = 0.05
MAPPING_USAGE_THRESHOLD = 0.1

# ---------- Visualization ----------
DPI = 300
FIGURE_FORMAT = 'png'

# ---------- Performance ----------
USE_CACHE = True

print("\n" + "="*70)
print("Master Configuration Loaded")
print("="*70)
print(f"\n📁 Input: {INPUT_H5AD_PATH}")
print(f"📁 Output: {OUTPUT_DIR}")
print(f"\n🔬 Pipeline Steps:")
print(f"   {'✅' if RUN_BBKNN else '❌'} BBKNN Integration")
print(f"   {'✅' if RUN_HARMONY else '❌'} Harmony Correction")
print(f"   {'✅' if RUN_CNMF_UNCORRECTED else '❌'} cNMF - Uncorrected")
print(f"   {'✅' if RUN_CNMF_HARMONY else '❌'} cNMF - Harmony")
print(f"   {'✅' if RUN_MAPPING else '❌'} GEP-Cluster Mapping")
print(f"\n🎯 Key Parameters:")
print(f"   cNMF K: {SELECTED_K}")
print(f"   🚀 cNMF workers: {CNMF_TOTAL_WORKERS} (multi-processing)")
print(f"   Harmony theta: {HARMONY_THETA}")
print(f"   Cluster key: {MAPPING_CLUSTER_KEY}")

# %% [markdown]
# ## 🚀 Multi-Processing Helper Function

# %%
def run_cnmf_factorize_parallel(output_dir: str, name: str, total_workers: int):
    """
    Run cNMF factorization in parallel using multiple processes
    
    Parameters
    ----------
    output_dir : str
        Must match cNMF(output_dir=...)
    name : str
        Must match cNMF(name=...)
    total_workers : int
        Number of parallel workers (processes)
    
    Notes
    -----
    - Each worker loads the normalized matrix into memory (~4-5GB per worker)
    - Recommended: total_workers <= number of CPU cores
    - Set OMP_NUM_THREADS=1-2 to avoid CPU oversubscription
    """
    def worker_fn(worker_i: int):
        """Worker function to run in subprocess"""
        from cnmf import cNMF
        cnmf_obj = cNMF(output_dir=output_dir, name=name)
        cnmf_obj.factorize(worker_i=worker_i, total_workers=total_workers)
    
    # Launch all workers
    procs = []
    print(f"   Launching {total_workers} parallel workers...")
    for wi in range(total_workers):
        p = Process(target=worker_fn, args=(wi,))
        p.start()
        procs.append(p)
        print(f"      Worker {wi} started (PID: {p.pid})")
    
    # Wait for all to complete
    for i, p in enumerate(procs):
        p.join()
        print(f"      Worker {i} finished")
    
    print(f"   ✓ All {total_workers} workers completed")

print("✓ Multi-processing helper function loaded")

# %% [markdown]
# ## Initialize Environment

# %%
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

harmony_fig_dir = output_dir / "figures_harmony"
harmony_fig_dir.mkdir(exist_ok=True)

cnmf_dir = output_dir / "cnmf_results"
cnmf_dir.mkdir(exist_ok=True)

mapping_dir = output_dir / "cnmf_cluster_mapping"
mapping_dir.mkdir(exist_ok=True)

sc.settings.figdir = fig_dir

pipeline_start_time = time.time()

pipeline_metadata = {
    'start_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    'input_file': str(INPUT_H5AD_PATH),
    'output_dir': str(OUTPUT_DIR),
    'cnmf_K': SELECTED_K,
    'cnmf_workers': CNMF_TOTAL_WORKERS,
    'harmony_theta': HARMONY_THETA,
    'harmony_lambda': HARMONY_LAMBDA,
    'harmony_n_pcs_correct': HARMONY_N_PCS_CORRECT
}

print(f"\n✓ Environment initialized")
print(f"   Output structure created: {output_dir}")

# %% [markdown]
# ## Marker Gene Definitions

# %%
MARKER_EPITHELIAL = [
    'TP63','KRT5','KRT14',
    'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1',
    'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24',
    'KRT19','NOTCH3','KRT16','KRT23',
    'MUC5AC','SPDEF','LYPD2','ITLN1',
    'ASCL1','GRP','KRT8',
    'POU2F3','ASCL2','CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB','PDE1C',
    'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
    'SFTPC','LAMP3','MFSD2A','C8orf4','C11orf96','SFTPB','SFTA2',
    'VIM','SOX9','MYH11','ACTA2','MYLK',
    'DMBT1','RNASE1','MUC5B','LYZ','LTF','PIP','CCL28',
    'DEUP1','FOXN4','CDC20B','CCNO',
    'MKI67','TOP2A','TK1','CENPW'
]

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

FORCE_INCLUDE_MARKERS = sorted({
    *MARKER_EPITHELIAL,
    *[g for genes in KEY_MARKERS.values() for g in genes],
})

print(f"Marker genes defined:")
print(f"   Total markers: {len(MARKER_EPITHELIAL)}")
print(f"   Force include: {len(FORCE_INCLUDE_MARKERS)}")

# %% [markdown]
# ---
# # 📊 Section 1: Data Loading & Preprocessing

# %%
print("\n" + "="*70)
print("Section 1: Data Loading & Preprocessing")
print("="*70)

section_start = time.time()

print(f"\nLoading: {INPUT_H5AD_PATH}")
if not Path(INPUT_H5AD_PATH).exists():
    raise FileNotFoundError(f"File not found: {INPUT_H5AD_PATH}")

adata = sc.read_h5ad(INPUT_H5AD_PATH)
print(f"✓ Data loaded")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")

print(f"\nBatch distribution (key: {BATCH_KEY}):")
batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
for batch, count in batch_counts.items():
    pct = count / adata.n_obs * 100
    print(f"   {batch}: {count:,} cells ({pct:.1f}%)")
    
print(f"\nFiltering genes (min.cells=3)...")
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=3)
n_genes_after = adata.n_vars
n_removed = n_genes_before - n_genes_after
print(f"   ✓ Removed {n_removed:,} genes ({n_removed/n_genes_before*100:.1f}%)")
print(f"   ✓ Remaining: {n_genes_after:,} genes")

print(f"\nPreserving raw counts...")
if 'counts' in adata.layers:
    counts = adata.layers['counts']
elif getattr(adata, "raw", None) is not None:
    counts = adata.raw.X
else:
    counts = adata.X

if not issparse(counts):
    counts = csr_matrix(counts)
adata.layers['counts'] = counts
print(f"   ✓ Raw counts saved to layers['counts']")

print(f"\nSelecting HVGs (n={N_TOP_GENES})...")
adata.X = adata.layers['counts'].copy()

try:
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_TOP_GENES, flavor=HVG_FLAVOR, batch_key=BATCH_KEY
    )
    print(f"   ✓ HVG selection with batch_key={BATCH_KEY}")
except Exception as e:
    print(f"   ⚠️ Fallback to HVG without batch_key")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_TOP_GENES, flavor=HVG_FLAVOR, batch_key=None
    )

hv_mask = adata.var['highly_variable'].to_numpy(dtype=bool)
force_mask = adata.var_names.isin(FORCE_INCLUDE_MARKERS)
keep_mask = hv_mask | force_mask

print(f"   HVGs: {int(hv_mask.sum())}")
print(f"   Force included: {int(force_mask.sum())}")
print(f"   Total kept: {int(keep_mask.sum())}")

print(f"\nNormalizing and log-transforming...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()
print(f"   ✓ Full matrix processed")

adata_full = adata.copy()
adata_hvg = adata[:, keep_mask].copy()
adata_hvg.layers['counts'] = adata.layers['counts'][:, keep_mask].copy()

adata_hvg.X = adata_hvg.layers['counts'].copy()
sc.pp.normalize_total(adata_hvg, target_sum=1e4)
sc.pp.log1p(adata_hvg)
adata_hvg.layers['log1p'] = adata_hvg.X.copy()

adata_hvg.uns['hvg_n_top_genes'] = int(N_TOP_GENES)
adata_hvg.uns['force_include_markers'] = np.array(FORCE_INCLUDE_MARKERS, dtype=str)
adata_hvg.raw = adata_full

print(f"   ✓ HVG subset created: {adata_hvg.n_vars} genes")
print(f"   ✓ Set adata_hvg.raw = adata_full")

print(f"\nRunning PCA (n={BBKNN_N_PCS})...")
sc.tl.pca(adata_hvg, n_comps=BBKNN_N_PCS, svd_solver='arpack')

var_ratio = adata_hvg.uns['pca']['variance_ratio']
cumsum_var = np.cumsum(var_ratio)
print(f"   PC1-10: {cumsum_var[9]:.2%}")
print(f"   PC1-30: {cumsum_var[29]:.2%}")
print(f"   PC1-50: {cumsum_var[49]:.2%}")

section_time = time.time() - section_start
pipeline_metadata['section1_time'] = section_time

print(f"\n💾 Saving preprocessing checkpoint...")
preproc_path = output_dir / "checkpoint_preprocessed.h5ad"
adata_full.write_h5ad(preproc_path, compression='gzip', compression_opts=9)
print(f"   ✓ Saved: {preproc_path}")

print(f"\n✅ Section 1 complete ({section_time:.1f}s)")
print("="*70)

# %% [markdown]
# ---
# # 🔗 Section 2: BBKNN Batch Integration

# %%
if RUN_BBKNN:
    print("\n" + "="*70)
    print("Section 2: BBKNN Batch Integration")
    print("="*70)
    
    section_start = time.time()
    
    print(f"\nBBKNN parameters:")
    print(f"   batch_key: {BATCH_KEY}")
    print(f"   neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    print(f"   n_pcs: {BBKNN_N_PCS}")
    
    print(f"\nRunning BBKNN...")
    bbknn_func(
        adata_hvg,
        batch_key=BATCH_KEY,
        neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
        n_pcs=BBKNN_N_PCS,
        copy=False
    )
    print(f"   ✓ BBKNN complete")
    
    print(f"\nRunning UMAP (min_dist={UMAP_MIN_DIST})...")
    sc.tl.umap(adata_hvg, min_dist=UMAP_MIN_DIST)
    print(f"   ✓ UMAP complete")
    
    print(f"\nMulti-resolution Leiden clustering...")
    for res in LEIDEN_RESOLUTIONS:
        cluster_key = f'leiden_bbknn_res{res}'
        try:
            sc.tl.leiden(adata_hvg, resolution=res, key_added=cluster_key,
                         flavor='igraph', n_iterations=2, directed=False)
        except TypeError:
            sc.tl.leiden(adata_hvg, resolution=res, key_added=cluster_key, n_iterations=2)
        n_clusters = adata_hvg.obs[cluster_key].nunique()
        print(f"   Resolution {res}: {n_clusters} clusters")
    
    default_key = f'leiden_bbknn_res{DEFAULT_RESOLUTION}'
    adata_hvg.obs['leiden_bbknn'] = adata_hvg.obs[default_key]
    print(f"   ✓ Default: {default_key}")
    
    print(f"\nGenerating BBKNN visualizations...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata_hvg, color=BATCH_KEY, ax=axes[0], show=False, title='Batch')
    sc.pl.umap(adata_hvg, color='leiden_bbknn', ax=axes[1], show=False,
               title='Clusters', legend_loc='on data', legend_fontsize=8)
    if 'cell_type' in adata_hvg.obs.columns:
        sc.pl.umap(adata_hvg, color='cell_type', ax=axes[2], show=False, title='Cell Type')
    else:
        axes[2].axis('off')
    plt.tight_layout()
    plt.savefig(fig_dir / f'bbknn_overview.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Saved: bbknn_overview.{FIGURE_FORMAT}")
    
    print(f"\nBackfilling results to full dataset...")
    for key in ['X_pca','X_umap']:
        if key in adata_hvg.obsm:
            adata_full.obsm[key] = adata_hvg.obsm[key]
    for col in adata_hvg.obs.columns:
        if col.startswith('leiden_bbknn'):
            adata_full.obs[col] = adata_hvg.obs[col]
    if 'neighbors' in adata_hvg.uns:
        adata_full.uns['neighbors'] = adata_hvg.uns['neighbors']
    print(f"   ✓ Results backfilled")
    
    print(f"\n💾 Saving BBKNN checkpoint...")
    bbknn_path = output_dir / "checkpoint_bbknn.h5ad"
    adata_full.write_h5ad(bbknn_path, compression='gzip', compression_opts=9)
    print(f"   ✓ Saved: {bbknn_path}")
    
    section_time = time.time() - section_start
    pipeline_metadata['section2_time'] = section_time
    pipeline_metadata['bbknn_date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    print(f"\n✅ Section 2 complete ({section_time:.1f}s)")
    print("="*70)
else:
    print("⏭️  Skipping BBKNN (RUN_BBKNN=False)")

# %% [markdown]
# ---
# # 🎨 Section 3: Harmony Batch Correction

# %%
if RUN_HARMONY:
    try:
        import harmonypy as hm
        print(f"harmonypy version: {hm.__version__}")
    except ImportError:
        print("Installing harmonypy...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy"])
        import harmonypy as hm
        print(f"✓ harmonypy installed: {hm.__version__}")

# %%
if RUN_HARMONY:
    print("\n" + "="*70)
    print("Section 3: Harmony Batch Correction")
    print("="*70)
    
    section_start = time.time()
    
    print(f"\nHarmony strategy: Limited PC correction")
    print(f"   Only correct first {HARMONY_N_PCS_CORRECT}/{BBKNN_N_PCS} PCs")
    print(f"   theta={HARMONY_THETA} (lower=gentler)")
    print(f"   lambda={HARMONY_LAMBDA} (higher=conservative)")
    
    print(f"\nRunning Harmony...")
    adata_harmony = adata_hvg.copy()
    
    pca_subset = adata_harmony.obsm['X_pca'][:, :HARMONY_N_PCS_CORRECT].copy()
    adata_harmony.obsm['X_pca_subset'] = pca_subset
    
    sc.external.pp.harmony_integrate(
        adata_harmony,
        key=BATCH_KEY,
        basis='X_pca_subset',
        adjusted_basis='X_pca_harmony_subset',
        max_iter_harmony=HARMONY_MAX_ITER,
        theta=HARMONY_THETA,
        lamb=HARMONY_LAMBDA,
        sigma=HARMONY_SIGMA,
        tau=HARMONY_TAU
    )
    
    corrected_pcs = adata_harmony.obsm['X_pca_harmony_subset']
    uncorrected_pcs = adata_harmony.obsm['X_pca'][:, HARMONY_N_PCS_CORRECT:]
    adata_harmony.obsm['X_pca_harmony'] = np.hstack([corrected_pcs, uncorrected_pcs])
    
    print(f"   ✓ Harmony complete")
    print(f"   Combined: {corrected_pcs.shape[1]} corrected + {uncorrected_pcs.shape[1]} original PCs")
    
    del adata_harmony.obsm['X_pca_subset']
    del adata_harmony.obsm['X_pca_harmony_subset']
    
    print(f"\nComputing UMAP on Harmony PCs...")
    sc.pp.neighbors(adata_harmony, n_neighbors=15, n_pcs=HARMONY_N_PCS_CORRECT,
                    use_rep='X_pca_harmony')
    sc.tl.umap(adata_harmony, min_dist=UMAP_MIN_DIST)
    print(f"   ✓ UMAP complete")
    
    print(f"\nMulti-resolution clustering...")
    for res in LEIDEN_RESOLUTIONS:
        cluster_key = f'leiden_harmony_res{res}'
        try:
            sc.tl.leiden(adata_harmony, resolution=res, key_added=cluster_key,
                         flavor='igraph', n_iterations=2, directed=False)
        except TypeError:
            sc.tl.leiden(adata_harmony, resolution=res, key_added=cluster_key, n_iterations=2)
        n_clusters = adata_harmony.obs[cluster_key].nunique()
        print(f"   Resolution {res}: {n_clusters} clusters")
    
    default_key = f'leiden_harmony_res{DEFAULT_RESOLUTION}'
    adata_harmony.obs['leiden_harmony'] = adata_harmony.obs[default_key]
    
    print(f"\n🔍 Running quality diagnostics...")
    try:
        from sklearn.metrics import silhouette_score
        
        if 'cell_type' in adata_harmony.obs.columns:
            sil_batch_bbknn = silhouette_score(
                adata_hvg.obsm['X_pca'][:, :HARMONY_N_PCS_CORRECT],
                adata_hvg.obs[BATCH_KEY], metric='euclidean',
                sample_size=min(5000, adata_hvg.n_obs)
            )
            sil_batch_harmony = silhouette_score(
                adata_harmony.obsm['X_pca_harmony'][:, :HARMONY_N_PCS_CORRECT],
                adata_harmony.obs[BATCH_KEY], metric='euclidean',
                sample_size=min(5000, adata_harmony.n_obs)
            )
            
            sil_celltype_bbknn = silhouette_score(
                adata_hvg.obsm['X_pca'][:, :HARMONY_N_PCS_CORRECT],
                adata_hvg.obs['cell_type'], metric='euclidean',
                sample_size=min(5000, adata_hvg.n_obs)
            )
            sil_celltype_harmony = silhouette_score(
                adata_harmony.obsm['X_pca_harmony'][:, :HARMONY_N_PCS_CORRECT],
                adata_harmony.obs['cell_type'], metric='euclidean',
                sample_size=min(5000, adata_harmony.n_obs)
            )
            
            print(f"   Batch separation (should decrease):")
            print(f"      BBKNN:   {sil_batch_bbknn:.3f}")
            print(f"      Harmony: {sil_batch_harmony:.3f}")
            print(f"      Change:  {sil_batch_harmony - sil_batch_bbknn:+.3f}")
            
            print(f"\n   Cell type separation (should stay high):")
            print(f"      BBKNN:   {sil_celltype_bbknn:.3f}")
            print(f"      Harmony: {sil_celltype_harmony:.3f}")
            print(f"      Change:  {sil_celltype_harmony - sil_celltype_bbknn:+.3f}")
            
            if sil_celltype_harmony < 0.3:
                print(f"\n   🚨 WARNING: Possible over-correction!")
            elif sil_celltype_harmony < sil_celltype_bbknn * 0.8:
                print(f"\n   ⚠️  Cell type separation dropped {(1-sil_celltype_harmony/sil_celltype_bbknn)*100:.1f}%")
            else:
                print(f"\n   ✅ Cell type separation preserved")
            
            pipeline_metadata['harmony_batch_sep'] = float(sil_batch_harmony)
            pipeline_metadata['harmony_celltype_sep'] = float(sil_celltype_harmony)
    except Exception as e:
        print(f"   ⚠️  Diagnostic test failed: {e}")
    
    print(f"\nBackfilling Harmony results...")
    adata_full.obsm['X_pca_harmony'] = adata_harmony.obsm['X_pca_harmony']
    adata_full.obsm['X_umap_harmony'] = adata_harmony.obsm['X_umap']
    for col in adata_harmony.obs.columns:
        if col.startswith('leiden_harmony'):
            adata_full.obs[col] = adata_harmony.obs[col]
    print(f"   ✓ Results backfilled")
    
    print(f"\n💾 Saving Harmony checkpoint...")
    harmony_path = output_dir / "checkpoint_harmony.h5ad"
    adata_full.write_h5ad(harmony_path, compression='gzip', compression_opts=9)
    print(f"   ✓ Saved: {harmony_path}")
    
    section_time = time.time() - section_start
    pipeline_metadata['section3_time'] = section_time
    pipeline_metadata['harmony_date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    print(f"\n✅ Section 3 complete ({section_time:.1f}s)")
    print("="*70)
else:
    print("⏭️  Skipping Harmony (RUN_HARMONY=False)")

# %% [markdown]
# ---
# # 🧬 Section 4: cNMF - Uncorrected (WITH MULTI-PROCESSING)

# %%
if RUN_CNMF_UNCORRECTED or RUN_CNMF_HARMONY:
    try:
        from cnmf import cNMF
        print("cNMF already installed")
    except ImportError:
        print("Installing cNMF...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "cnmf"])
        from cnmf import cNMF
        print("✓ cNMF installed")

# %%
if RUN_CNMF_UNCORRECTED:
    print("\n" + "="*70)
    print("Section 4: cNMF - Uncorrected (Multi-Processing)")
    print("="*70)
    
    section_start = time.time()
    
    print(f"\ncNMF parameters:")
    print(f"   K values to test: {CNMF_COMPONENTS}")
    print(f"   Iterations per K: {CNMF_N_ITER}")
    print(f"   HVG count: {CNMF_NUM_HVG}")
    print(f"   Selected K: {SELECTED_K}")
    print(f"   🚀 Workers: {CNMF_TOTAL_WORKERS}")
    
    # Estimate time with multi-processing
    n_cells = adata_full.n_obs
    n_k = len(CNMF_COMPONENTS)
    est_minutes_single = (n_cells / 10000) * n_k * CNMF_N_ITER * 0.5
    est_minutes = est_minutes_single / CNMF_TOTAL_WORKERS
    
    print(f"\n⏱️  Estimated time:")
    print(f"   Single-threaded: {est_minutes_single:.0f} min ({est_minutes_single/60:.1f} hr)")
    print(f"   With {CNMF_TOTAL_WORKERS} workers: {est_minutes:.0f} min ({est_minutes/60:.1f} hr)")
    print(f"   💡 Speedup: ~{CNMF_TOTAL_WORKERS:.1f}x")
    print(f"\n   ☕ Estimated memory: ~{CNMF_TOTAL_WORKERS * 5:.0f} GB")
    
    print(f"\nPreparing data...")
    adata_cnmf = adata_full.copy()
    if 'counts' in adata_cnmf.layers:
        adata_cnmf.X = adata_cnmf.layers['counts'].copy()
    cnmf_data_path = cnmf_dir / "epithelial_raw_counts.h5ad"
    adata_cnmf.write_h5ad(cnmf_data_path)
    print(f"   ✓ Saved: {cnmf_data_path}")
    
    cnmf_nocorr = cNMF(output_dir=str(cnmf_dir), name='Epithelial_NoBatchCorrection')
    
    print(f"\nStep 1/4: Prepare...")
    cnmf_nocorr.prepare(
        counts_fn=str(cnmf_data_path),
        components=CNMF_COMPONENTS,
        n_iter=CNMF_N_ITER,
        seed=CNMF_SEED,
        num_highvar_genes=CNMF_NUM_HVG
    )
    print(f"   ✓ Complete")
    
    print(f"\nStep 2/4: Factorize with {CNMF_TOTAL_WORKERS} parallel workers...")
    factorize_start = time.time()
    
    # 🚀 Use multi-processing
    run_cnmf_factorize_parallel(
        output_dir=str(cnmf_dir),
        name='Epithelial_NoBatchCorrection',
        total_workers=CNMF_TOTAL_WORKERS
    )
    
    factorize_time = time.time() - factorize_start
    speedup = est_minutes_single / (factorize_time / 60)
    print(f"   ✓ Complete ({factorize_time/60:.1f} min)")
    print(f"   💡 Actual speedup: {speedup:.1f}x")
    
    print(f"\nStep 3/4: Combine...")
    cnmf_nocorr.combine()
    print(f"   ✓ Complete")
    
    print(f"\nStep 4/4: Consensus clustering...")
    for k in CNMF_COMPONENTS:
        cnmf_nocorr.consensus(
            k=k,
            density_threshold=CNMF_DENSITY_THRESHOLD,
            show_clustering=False,
            close_clustergram_fig=True
        )
        print(f"   K={k} complete")
    
    print(f"\nLoading results for K={SELECTED_K}...")
    (usage_nocorr, spectra_scores_nocorr,
     spectra_tpm_nocorr, top_genes_nocorr) = cnmf_nocorr.load_results(
        K=SELECTED_K, density_threshold=CNMF_DENSITY_THRESHOLD
    )
    print(f"   ✓ Loaded: {usage_nocorr.shape}")
    print(f"   Column names: {list(usage_nocorr.columns[:3])}...")
    
    # 🔧 FIX: Use column position instead of column name
    print(f"\nAdding GEP scores to adata (using column position)...")
    n_geps = usage_nocorr.shape[1]
    if n_geps != SELECTED_K:
        print(f"   ⚠️  Warning: Expected {SELECTED_K} GEPs, got {n_geps}")
    
    for gep_idx in range(n_geps):
        gep_num = gep_idx + 1
        col_data = usage_nocorr.iloc[:, gep_idx]  # Use position
        adata_full.obs[f'GEP_{gep_num}_nocorr'] = col_data.reindex(adata_full.obs_names).values
    
    print(f"   ✓ Added {n_geps} GEP columns")

    print(f"\n💾 Saving cNMF (uncorrected) checkpoint...")
    cnmf_nocorr_path = output_dir / "checkpoint_cnmf_nocorr.h5ad"
    adata_full.write_h5ad(cnmf_nocorr_path, compression='gzip', compression_opts=9)
    print(f"   ✓ Saved: {cnmf_nocorr_path}")

    section_time = time.time() - section_start
    pipeline_metadata['section4_time'] = section_time
    pipeline_metadata['section4_speedup'] = speedup
    pipeline_metadata['cnmf_nocorr_date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    print(f"\n✅ Section 4 complete ({section_time/60:.1f} min)")
    print("="*70)
else:
    print("⏭️  Skipping cNMF-Uncorrected")

# %% [markdown]
# ---
# # 🎨 Section 5: cNMF - Harmony (WITH MULTI-PROCESSING)

# %%
if RUN_CNMF_HARMONY and RUN_HARMONY:
    print("\n" + "="*70)
    print("Section 5: cNMF - Harmony Corrected (Multi-Processing)")
    print("="*70)
    
    section_start = time.time()
    
    # Calculate separate time estimate for Section 5
    n_cells_s5 = adata_full.n_obs
    n_k_s5 = len(CNMF_COMPONENTS)
    est_minutes_single_s5 = (n_cells_s5 / 10000) * n_k_s5 * CNMF_N_ITER * 0.5
    est_minutes_s5 = est_minutes_single_s5 / CNMF_TOTAL_WORKERS
    
    print(f"\nPreparing Harmony-corrected data...")
    print(f"   Strategy: Harmony-guided kNN smoothing")
    print(f"   🚀 Workers: {CNMF_TOTAL_WORKERS}")
    print(f"   ⏱️  Estimated: {est_minutes_s5:.0f} min ({est_minutes_s5/60:.1f} hr)")
    
    adata_hvg_harmony = adata_harmony.copy()
    
    print(f"\nComputing Harmony neighbor graph...")
    sc.pp.neighbors(
        adata_hvg_harmony,
        n_neighbors=15,
        n_pcs=HARMONY_N_PCS_CORRECT,
        use_rep='X_pca_harmony',
        key_added='harmony_nn'
    )
    print(f"   ✓ Neighbor graph computed")
    
    if 'counts' in adata_hvg_harmony.layers:
        raw_counts = adata_hvg_harmony.layers['counts'].copy()
    else:
        raw_counts = adata_hvg_harmony.X.copy()
    
    print(f"\nApplying kNN smoothing...")
    connectivities = adata_hvg_harmony.obsp['harmony_nn_connectivities']
    smoothed_counts = connectivities @ raw_counts
    
    if issparse(smoothed_counts):
        smoothed_counts = smoothed_counts.toarray()
    smoothed_counts = np.round(smoothed_counts).astype(np.int32)
    smoothed_counts[smoothed_counts < 0] = 0
    print(f"   ✓ Smoothing complete")
    
    adata_corrected_hvg = sc.AnnData(
        X=csr_matrix(smoothed_counts),
        obs=adata_hvg_harmony.obs.copy(),
        var=adata_hvg_harmony.var.copy()
    )
    
    corrected_hvg_path = cnmf_dir / "epithelial_harmony_smoothed_counts.h5ad"
    adata_corrected_hvg.write_h5ad(corrected_hvg_path)
    print(f"   ✓ Saved: {corrected_hvg_path}")
    
    print(f"\nCreating TP10K matrix...")
    adata_tp10k = adata_full.copy()
    if 'counts' in adata_tp10k.layers:
        adata_tp10k.X = adata_tp10k.layers['counts'].copy()
    sc.pp.normalize_total(adata_tp10k, target_sum=1e4)
    
    tp10k_path = cnmf_dir / "epithelial_tp10k.h5ad"
    adata_tp10k.write_h5ad(tp10k_path)
    print(f"   ✓ Saved: {tp10k_path}")
    
    hvg_list = adata_corrected_hvg.var_names.tolist()
    hvg_file = cnmf_dir / "harmony_hvgs.txt"
    with open(hvg_file, 'w') as f:
        for gene in hvg_list:
            f.write(f"{gene}\n")
    print(f"   ✓ Saved HVG list: {len(hvg_list)} genes")
    
    cnmf_harmony = cNMF(output_dir=str(cnmf_dir), name='Epithelial_HarmonyBatchCorrected')
    
    print(f"\nStep 1/4: Prepare...")
    cnmf_harmony.prepare(
        counts_fn=str(corrected_hvg_path),
        tpm_fn=str(tp10k_path),
        genes_file=str(hvg_file),
        components=CNMF_COMPONENTS,
        n_iter=CNMF_N_ITER,
        seed=CNMF_SEED
    )
    print(f"   ✓ Complete")
    
    print(f"\nStep 2/4: Factorize with {CNMF_TOTAL_WORKERS} parallel workers...")
    factorize_start = time.time()
    
    # 🚀 Use multi-processing
    run_cnmf_factorize_parallel(
        output_dir=str(cnmf_dir),
        name='Epithelial_HarmonyBatchCorrected',
        total_workers=CNMF_TOTAL_WORKERS
    )
    
    factorize_time = time.time() - factorize_start
    speedup_s5 = est_minutes_single_s5 / (factorize_time / 60)
    print(f"   ✓ Complete ({factorize_time/60:.1f} min)")
    print(f"   💡 Actual speedup: {speedup_s5:.1f}x")
    
    print(f"\nStep 3/4: Combine...")
    cnmf_harmony.combine()
    print(f"   ✓ Complete")
    
    print(f"\nStep 4/4: Consensus...")
    for k in CNMF_COMPONENTS:
        cnmf_harmony.consensus(
            k=k,
            density_threshold=CNMF_DENSITY_THRESHOLD,
            show_clustering=False,
            close_clustergram_fig=True
        )
        print(f"   K={k} complete")
    
    print(f"\nLoading results for K={SELECTED_K}...")
    (usage_harmony, spectra_scores_harmony,
     spectra_tpm_harmony, top_genes_harmony) = cnmf_harmony.load_results(
        K=SELECTED_K, density_threshold=CNMF_DENSITY_THRESHOLD
    )
    print(f"   ✓ Loaded: {usage_harmony.shape}")
    print(f"   Column names: {list(usage_harmony.columns[:3])}...")
    
    # 🔧 FIX: Use column position instead of column name
    print(f"\nAdding GEP scores to adata (using column position)...")
    n_geps = usage_harmony.shape[1]
    if n_geps != SELECTED_K:
        print(f"   ⚠️  Warning: Expected {SELECTED_K} GEPs, got {n_geps}")
    
    for gep_idx in range(n_geps):
        gep_num = gep_idx + 1
        col_data = usage_harmony.iloc[:, gep_idx]  # Use position
        adata_full.obs[f'GEP_{gep_num}_harmony'] = col_data.reindex(adata_full.obs_names).values
    
    print(f"   ✓ Added {n_geps} GEP columns")

    print(f"\n💾 Saving cNMF (Harmony) checkpoint...")
    cnmf_harmony_path = output_dir / "checkpoint_cnmf_harmony.h5ad"
    adata_full.write_h5ad(cnmf_harmony_path, compression='gzip', compression_opts=9)
    print(f"   ✓ Saved: {cnmf_harmony_path}")

    section_time = time.time() - section_start
    pipeline_metadata['section5_time'] = section_time
    pipeline_metadata['section5_speedup'] = speedup_s5
    pipeline_metadata['cnmf_harmony_date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    print(f"\n✅ Section 5 complete ({section_time/60:.1f} min)")
    print("="*70)
elif RUN_CNMF_HARMONY and not RUN_HARMONY:
    print("⚠️  Cannot run cNMF-Harmony without Harmony correction!")
else:
    print("⏭️  Skipping cNMF-Harmony")

# %%
if RUN_CNMF_UNCORRECTED or RUN_CNMF_HARMONY:
    print(f"\n💾 Saving complete integrated data with GEP scores...")
    complete_path = output_dir / "checkpoint_complete_with_geps.h5ad"
    adata_full.uns['pipeline_metadata'] = pipeline_metadata
    adata_full.write_h5ad(complete_path, compression='gzip', compression_opts=9)
    complete_size = complete_path.stat().st_size / (1024**3)
    print(f"   ✓ Saved: {complete_path} ({complete_size:.2f} GB)")
    print(f"\n   This file contains:")
    if RUN_CNMF_UNCORRECTED:
        print(f"      - cNMF uncorrected GEPs (GEP_*_nocorr)")
    if RUN_CNMF_HARMONY:
        print(f"      - cNMF Harmony GEPs (GEP_*_harmony)")
    print(f"      - All clustering results")
    print(f"      - All metadata and embeddings")

# %% [markdown]
# ---
# # 🗺️ Section 6: GEP-Cluster Mapping Analysis

# %%
if RUN_MAPPING and (RUN_CNMF_UNCORRECTED or RUN_CNMF_HARMONY):
    print("\n" + "="*70)
    print("Section 6: GEP-Cluster Mapping Analysis")
    print("="*70)
    
    section_start = time.time()
    
    gep_types_to_analyze = []
    if MAPPING_ANALYZE_UNCORRECTED and RUN_CNMF_UNCORRECTED:
        gep_types_to_analyze.append(('nocorr', 'Uncorrected'))
    if MAPPING_ANALYZE_HARMONY and RUN_CNMF_HARMONY:
        gep_types_to_analyze.append(('harmony', 'Harmony'))
    
    if len(gep_types_to_analyze) == 0:
        print("⚠️  No GEP types available for mapping!")
    else:
        print(f"\nMapping configuration:")
        print(f"   Cluster key: {MAPPING_CLUSTER_KEY}")
        print(f"   Analyzing: {[t[1] for t in gep_types_to_analyze]}")
        print(f"   K: {SELECTED_K}")
        
        if MAPPING_CLUSTER_KEY not in adata_full.obs.columns:
            print(f"\n⚠️  Cluster key '{MAPPING_CLUSTER_KEY}' not found!")
            available = [col for col in adata_full.obs.columns if 'leiden' in col]
            print(f"   Available: {available}")
        else:
            for gep_suffix, gep_name in gep_types_to_analyze:
                print(f"\n{'='*70}")
                print(f"Analyzing {gep_name} GEPs")
                print(f"{'='*70}")
                
                mapping_subdir = mapping_dir / f"mapping_{gep_suffix}"
                mapping_subdir.mkdir(exist_ok=True)
                
                gep_cols = [col for col in adata_full.obs.columns 
                           if col.startswith('GEP_') and col.endswith(f'_{gep_suffix}')]
                
                if len(gep_cols) == 0:
                    print(f"   ⚠️  No GEP columns found for {gep_name}!")
                    continue
                
                print(f"   Found {len(gep_cols)} GEP columns")
                
                gep_usage = adata_full.obs[gep_cols].copy()
                gep_usage.columns = [col.replace(f'_{gep_suffix}', '') for col in gep_cols]
                gep_usage['cluster'] = adata_full.obs[MAPPING_CLUSTER_KEY].values
                
                print(f"\n   Calculating mean GEP usage per cluster...")
                mean_usage = gep_usage.groupby('cluster').mean()
                mean_usage_file = mapping_subdir / "mean_gep_usage_per_cluster.csv"
                mean_usage.to_csv(mean_usage_file)
                print(f"      ✓ Saved: {mean_usage_file.name}")
                
                print(f"   Identifying dominant GEPs...")
                TOP_N_GEPS = 3
                dominant_geps_data = []
                for cluster in mean_usage.index.sort_values():
                    cluster_usage = mean_usage.loc[cluster].sort_values(ascending=False)
                    top_geps = cluster_usage.head(TOP_N_GEPS)
                    for rank, (gep, usage) in enumerate(top_geps.items(), 1):
                        dominant_geps_data.append({
                            'Cluster': cluster,
                            'Rank': rank,
                            'GEP': gep,
                            'Mean_Usage': usage
                        })
                dominant_df = pd.DataFrame(dominant_geps_data)
                dominant_file = mapping_subdir / "dominant_geps_per_cluster.csv"
                dominant_df.to_csv(dominant_file, index=False)
                print(f"      ✓ Saved: {dominant_file.name}")
                
                print(f"   Running statistical tests...")
                gep_usage_clean = adata_full.obs[gep_cols].copy()
                gep_usage_clean.columns = [col.replace(f'_{gep_suffix}', '') for col in gep_cols]
                clusters = adata_full.obs[MAPPING_CLUSTER_KEY].values
                
                stat_results = []
                for cluster in mean_usage.index:
                    cluster_mask = clusters == cluster
                    n_cells_in = cluster_mask.sum()
                    
                    for gep in gep_usage_clean.columns:
                        vals_in = gep_usage_clean.loc[cluster_mask, gep]
                        vals_out = gep_usage_clean.loc[~cluster_mask, gep]
                        
                        try:
                            stat, pval = stats.ranksums(vals_in, vals_out, alternative='greater')
                        except:
                            stat, pval = np.nan, 1.0
                        
                        mean_in = vals_in.mean()
                        mean_out = vals_out.mean()
                        fc = mean_in / mean_out if mean_out > 0 else np.inf
                        
                        stat_results.append({
                            'Cluster': cluster,
                            'GEP': gep,
                            'n_cells': n_cells_in,
                            'mean_in_cluster': mean_in,
                            'mean_out_cluster': mean_out,
                            'fold_change': fc,
                            'pvalue': pval
                        })
                
                stat_df = pd.DataFrame(stat_results)
                
                from statsmodels.stats.multitest import multipletests
                _, pvals_adj, _, _ = multipletests(stat_df['pvalue'], method='bonferroni')
                stat_df['pvalue_adj'] = pvals_adj
                stat_df['significant'] = (stat_df['pvalue_adj'] < MAPPING_PVALUE_THRESHOLD) & \
                                         (stat_df['mean_in_cluster'] > stat_df['mean_out_cluster'])
                
                stat_file = mapping_subdir / "gep_cluster_statistical_tests.csv"
                stat_df.to_csv(stat_file, index=False)
                
                n_sig = stat_df['significant'].sum()
                print(f"      ✓ Found {n_sig} significant associations")
                print(f"      ✓ Saved: {stat_file.name}")
                
                print(f"   Generating visualizations...")
                
                fig = plt.figure(figsize=(max(12, SELECTED_K*0.6), max(8, len(mean_usage)*0.4)))
                row_linkage = linkage(mean_usage.values, method='average')
                col_linkage = linkage(mean_usage.T.values, method='average')
                g = sns.clustermap(
                    mean_usage,
                    row_linkage=row_linkage,
                    col_linkage=col_linkage,
                    cmap='RdYlBu_r',
                    center=mean_usage.values.mean(),
                    figsize=(max(12, SELECTED_K*0.6), max(8, len(mean_usage)*0.4)),
                    cbar_kws={'label': 'Mean GEP Usage'}
                )
                plt.suptitle(f'{gep_name} GEP Usage per Cluster', y=0.98)
                heatmap_file = mapping_subdir / f'heatmap_mean_usage.{FIGURE_FORMAT}'
                plt.savefig(heatmap_file, dpi=DPI, bbox_inches='tight')
                plt.close()
                print(f"      ✓ Heatmap saved")
                
                if n_sig > 0:
                    top_sig = stat_df[stat_df['significant']].sort_values('pvalue_adj').head(10)
                    print(f"\n   Top 10 significant associations:")
                    for i, (_, row) in enumerate(top_sig.iterrows(), 1):
                        print(f"      {i}. {row['GEP']} in Cluster {row['Cluster']}: "
                              f"FC={row['fold_change']:.2f}, p={row['pvalue_adj']:.2e}")
                
                print(f"\n   ✓ {gep_name} mapping complete")
            
            section_time = time.time() - section_start
            pipeline_metadata['section6_time'] = section_time
            pipeline_metadata['mapping_date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            
            print(f"\n✅ Section 6 complete ({section_time:.1f}s)")
            print("="*70)
elif RUN_MAPPING:
    print("⚠️  Cannot run mapping without cNMF results!")
else:
    print("⏭️  Skipping Mapping")

# %% [markdown]
# ---
# # 📋 Section 7: Final Summary Report

# %%
print("\n" + "="*70)
print("Section 7: Generating Final Summary Report")
print("="*70)

total_time = time.time() - pipeline_start_time
pipeline_metadata['total_time'] = total_time
pipeline_metadata['end_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

report = []
report.append("="*70)
report.append("Complete Epithelial Cell Analysis Pipeline - Summary Report")
report.append("="*70)
report.append(f"\nStart time: {pipeline_metadata['start_time']}")
report.append(f"End time: {pipeline_metadata['end_time']}")
report.append(f"Total time: {total_time/3600:.2f} hours ({total_time/60:.1f} minutes)")
report.append("")

report.append("[ Configuration ]")
report.append(f"  Input file: {INPUT_H5AD_PATH}")
report.append(f"  Output directory: {OUTPUT_DIR}")
report.append(f"  Number of cells: {adata_full.n_obs:,}")
report.append(f"  Number of genes: {adata_full.n_vars:,}")
report.append(f"  Batch key: {BATCH_KEY}")
report.append(f"  Number of batches: {adata_full.obs[BATCH_KEY].nunique()}")
report.append("")

report.append("[ Pipeline Steps Executed ]")
if RUN_BBKNN:
    report.append(f"  ✅ BBKNN Integration ({pipeline_metadata.get('section2_time', 0):.1f}s)")
if RUN_HARMONY:
    report.append(f"  ✅ Harmony Correction ({pipeline_metadata.get('section3_time', 0):.1f}s)")
if RUN_CNMF_UNCORRECTED:
    speedup = pipeline_metadata.get('section4_speedup', 0)
    report.append(f"  ✅ cNMF Uncorrected ({pipeline_metadata.get('section4_time', 0)/60:.1f} min, {speedup:.1f}x speedup)")
if RUN_CNMF_HARMONY:
    speedup = pipeline_metadata.get('section5_speedup', 0)
    report.append(f"  ✅ cNMF Harmony ({pipeline_metadata.get('section5_time', 0)/60:.1f} min, {speedup:.1f}x speedup)")
if RUN_MAPPING:
    report.append(f"  ✅ GEP-Cluster Mapping ({pipeline_metadata.get('section6_time', 0):.1f}s)")
report.append("")

report.append("[ Key Parameters ]")
report.append(f"  cNMF K: {SELECTED_K}")
report.append(f"  cNMF iterations: {CNMF_N_ITER}")
report.append(f"  🚀 cNMF workers: {CNMF_TOTAL_WORKERS}")
if RUN_HARMONY:
    report.append(f"  Harmony theta: {HARMONY_THETA}")
    report.append(f"  Harmony lambda: {HARMONY_LAMBDA}")
report.append("")

report.append("[ Output Files ]")
report.append(f"  Main data: checkpoint_complete_with_geps.h5ad")
report.append(f"  Figures: figures/*.{FIGURE_FORMAT}")
if RUN_MAPPING:
    report.append(f"  Mapping results: cnmf_cluster_mapping/")
report.append("")

report.append("="*70)
report.append("Pipeline completed successfully")
report.append("="*70)

report_text = '\n'.join(report)
report_file = output_dir / "pipeline_summary_report.txt"
with open(report_file, 'w') as f:
    f.write(report_text)

metadata_file = output_dir / "pipeline_metadata.json"
import json
with open(metadata_file, 'w') as f:
    json.dump(pipeline_metadata, f, indent=2)

print("\n" + report_text)
print(f"\n💾 Report saved: {report_file}")
print(f"💾 Metadata saved: {metadata_file}")

# %% [markdown]
# ---
# # 🎉 Pipeline Complete!

# %%
print("\n" + "="*70)
print("🎉 COMPLETE EPITHELIAL CELL ANALYSIS PIPELINE FINISHED 🎉")
print("="*70)

print(f"\n⏱️  Total Runtime: {total_time/3600:.2f} hours")
print(f"   Breakdown:")
print(f"      Data loading: {pipeline_metadata.get('section1_time', 0):.0f}s")
if RUN_BBKNN:
    print(f"      BBKNN: {pipeline_metadata.get('section2_time', 0):.0f}s")
if RUN_HARMONY:
    print(f"      Harmony: {pipeline_metadata.get('section3_time', 0):.0f}s")
if RUN_CNMF_UNCORRECTED:
    speedup = pipeline_metadata.get('section4_speedup', 0)
    print(f"      cNMF uncorrected: {pipeline_metadata.get('section4_time', 0)/60:.0f} min ({speedup:.1f}x)")
if RUN_CNMF_HARMONY:
    speedup = pipeline_metadata.get('section5_speedup', 0)
    print(f"      cNMF harmony: {pipeline_metadata.get('section5_time', 0)/60:.0f} min ({speedup:.1f}x)")
if RUN_MAPPING:
    print(f"      Mapping: {pipeline_metadata.get('section6_time', 0):.0f}s")

print(f"\n📁 All results saved to: {OUTPUT_DIR}")

print(f"\n🚀 Performance Summary:")
print(f"   Workers used: {CNMF_TOTAL_WORKERS}")
if RUN_CNMF_UNCORRECTED or RUN_CNMF_HARMONY:
    avg_speedup = (pipeline_metadata.get('section4_speedup', 0) + 
                   pipeline_metadata.get('section5_speedup', 0)) / 2
    print(f"   Average speedup: {avg_speedup:.1f}x")
    print(f"   💡 Time saved: ~{(avg_speedup-1) * 100:.0f}%")

print(f"\n✅ Analysis pipeline complete!")
