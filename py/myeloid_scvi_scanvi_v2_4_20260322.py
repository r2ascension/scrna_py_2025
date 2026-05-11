# %% [markdown]
# # Myeloid Cell Refinement Pipeline v2.4 PRODUCTION
# 
# **Author:** r2end | **Date:** 2026-03-22 | **Version:** 2.4 PRODUCTION
# 
# ## Changes from v2.3
# - P0-4: `.cat.categories` crash — restore Categorical after map/fillna
# - P1-1: stress_score NO [0,1] normalization (QRM rule; match query pipelines)
# - P1-2: scArches dry-run validation after scVI/scANVI save (QRM v2.8 Sec 9.2)
# - P1-3: var_names.csv saved in BOTH model directories (QRM v2.8 Sec 9.3)
# - P1-4: categorical_covariate_levels saved AFTER batch filtering (stale levels fix)
# - P2-1: adata_scvi construction simplified (redundant HVG indexing removed)
# - P2-2: UMAP rasterization via sc.settings.vector_friendly (QRM v3.6)
# 
# **Specs:** 50-70GB peak | 2-3h GPU | scArches-ready output
# 

# %% [markdown]
# ---
# ## 0. Setup
# 

# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
for var in ["OMP", "OPENBLAS", "MKL", "NUMEXPR"]:
    os.environ[f"{var}_NUM_THREADS"] = "1"

import sys, gc
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import scanpy as sc
import scvi
import torch

print(f"scanpy {sc.__version__} | scvi {scvi.__version__} | GPU: {torch.cuda.is_available()}")


# %%
# CONFIGURATION
INPUT_H5AD = "/home/h2048/data/py/0128/myeloid_analysis_unified/results/subcluster_unified_v2_20260128/adata_myeloid_subclustered_FINAL_v2_20260128.h5ad"
OUTPUT_DIR = Path("/home/h2048/data/py/0322/myeloid_validation_optimized")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CELLTYPE_L2_COL = 'cell_type_L2'
CELLTYPE_L3_COL = 'cell_type_L3'
CELLTYPE_L3_REFINED_COL = 'cell_type_L3_refined'
BATCH_KEY = 'sample'
TISSUE_KEY = 'tissue'

# scArches covariates (CRITICAL)
REQUIRED_COVARIATES = [BATCH_KEY, TISSUE_KEY, 'pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']

# scVI/scANVI (scArches-optimized)
RANDOM_SEED = 42
N_HVG = 4000
N_LATENT = 150
N_LAYERS = 2
DROPOUT_RATE = 0.2
MAX_EPOCHS_SCVI = 400
MAX_EPOCHS_SCANVI = 200
BATCH_SIZE = 256
LEARNING_RATE = 1e-3
ENCODE_COVARIATES = True
USE_LAYER_NORM = "both"
USE_BATCH_NORM = "none"
N_JOBS = 48

print(f"Config: n_latent={N_LATENT} (scArches), encode_covariates={ENCODE_COVARIATES}")


# %%
sc.settings.verbosity = 3
sc.settings.n_jobs = N_JOBS
sc.settings.set_figure_params(dpi=100, facecolor='white')
np.random.seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED
torch.manual_seed(RANDOM_SEED)
print("Settings configured")


# %%
# Stress genes (147)
STRESS_SIGNATURE_GENES = [
    "ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3","BANF1",
    "CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2","CNOT4","CNOT6",
    "CXXC1","DCP1A","DCP2","DCTN1","DDIT4","DDX10","DKC1","DNAJA4","DNAJB9",
    "DNAJC3","EDC4","EDEM1","EEF2","EIF2AK3","EIF2S1","EIF4A1","EIF4A2","EIF4A3",
    "EIF4E","EIF4EBP1","EIF4G1","ERN1","ERO1A","EXOC2","EXOSC1","EXOSC10",
    "EXOSC2","EXOSC4","EXOSC5","EXOSC9","FKBP14","FUS","GEMIN4","GOSR2","H2AX",
    "HERPUD1","HSP90B1","HSPA5","HSPA9","HYOU1","IARS1","IFIT1","IGFBP1","IMP3",
    "KDELR3","KHSRP","KIF5B","LSM1","LSM4","MTHFD2","NFYA","NFYB","NHP2","NOLC1",
    "NOP14","NOP56","NPM1","NABP1","PAIP1","PARN","PDIA5","PDIA6","POP4","PREB",
    "PSAT1","RPS14","RRP9","SDAD1","SEC11A","SEC31A","SERP1","SHC1","MTREX",
    "SLC1A4","SLC30A5","SLC7A5","SPCS1","SPCS3","SRPRA","SRPRB","SSR1","STC2",
    "TARS1","TATDN2","TSPYL2","SKIC3","TUBB2A","VEGFA","WFS1","WIPI1","XBP1",
    "XPOT","YIF1A","YWHAZ","ZBTB17"
]

# Cell cycle genes
S_GENES = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP',
    'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2',
    'RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
    'USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8'
]

G2M_GENES = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
    'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2',
    'CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B',
    'HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2',
    'DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN',
    'LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA'
]

print(f"Gene lists: Stress={len(STRESS_SIGNATURE_GENES)}, S={len(S_GENES)}, G2M={len(G2M_GENES)}")


# %% [markdown]
# ---
# ## 1. Load Data
# 

# %%
print("="*80)
print("1. LOADING DATA")
print("="*80)

adata = sc.read_h5ad(INPUT_H5AD)

print(f"Cells: {adata.n_obs:,} | Genes: {adata.n_vars:,}")
print(f"Layers: {list(adata.layers.keys())}")
print(f"Batches: {adata.obs[BATCH_KEY].nunique()}")

if 'counts' not in adata.layers:
    raise ValueError("Missing 'counts' layer!")
print("Data loaded")


# %% [markdown]
# ---
# ## 1.5 Calculate Covariates (on full gene space, BEFORE any subsetting)
# 

# %%
print("\n" + "-"*70)
print("Calculating MT%...")
print("-"*70)

mt_genes = (adata.var_names.str.startswith('MT-') |
            adata.var_names.str.startswith('Mt-') |
            adata.var_names.str.startswith('mt-'))
n_mt = mt_genes.sum()
print(f"  MT genes: {n_mt}")

if n_mt > 0:
    counts = adata.layers['counts']
    if issparse(counts):
        total = np.array(counts.sum(axis=1)).flatten()
        mt = np.array(counts[:, mt_genes].sum(axis=1)).flatten()
    else:
        total = counts.sum(axis=1)
        mt = counts[:, mt_genes].sum(axis=1)

    mt_pct = np.zeros_like(total, dtype=np.float32)
    mask = total > 0
    mt_pct[mask] = (mt[mask] / total[mask]) * 100
    adata.obs['pct_counts_mt'] = mt_pct
    print(f"  Mean MT%: {mt_pct.mean():.2f}%")
else:
    adata.obs['pct_counts_mt'] = np.float32(0.0)


# %%
# ========================================================================
# P1-1 FIX (v2.4): stress_score = mean(log1p(CPM)), NO [0,1] normalization
# QRM rule: "Stress score = mean of log1p(CPM) across available stress genes;
# no [0,1] normalization". This ensures consistency with query pipelines.
# ========================================================================
print("\n" + "-"*70)
print("Calculating stress score (NO [0,1] normalization)...")
print("-"*70)

stress_genes = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
n_stress = len(stress_genes)
print(f"  Stress genes: {n_stress}/{len(STRESS_SIGNATURE_GENES)}")

if n_stress >= 10:
    counts = adata.layers['counts']
    stress_idx = np.array([i for i, g in enumerate(adata.var_names) if g in stress_genes])

    if issparse(counts):
        stress_counts = counts[:, stress_idx]
        total = np.array(counts.sum(axis=1)).flatten().astype(np.float64)
        from scipy.sparse import diags
        scale = diags(1e4 / np.where(total > 0, total, 1))
        stress_cpm = (scale @ stress_counts).tocsr()
        stress_cpm.data = np.log1p(stress_cpm.data).astype(np.float32)
        raw_score = np.array(stress_cpm.mean(axis=1)).flatten()
    else:
        stress_counts = counts[:, stress_idx].astype(np.float64)
        total = counts.sum(axis=1, keepdims=True)
        cpm = stress_counts * (1e4 / np.where(total > 0, total, 1))
        raw_score = np.log1p(cpm).mean(axis=1)

    # [v2.4] NO [0,1] normalization — raw mean(log1p(CPM)) as-is
    adata.obs['stress_score'] = raw_score.astype(np.float32)

    del stress_counts
    if issparse(counts):
        del stress_cpm
    del raw_score
    gc.collect()
    print(f"  Mean: {adata.obs['stress_score'].mean():.4f} (raw, no [0,1])")
else:
    adata.obs['stress_score'] = np.float32(0.0)


# %%
print("\n" + "-"*70)
print("Calculating cell cycle (MEMORY OPTIMIZED)...")
print("-"*70)

s_genes = [g for g in S_GENES if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]
print(f"  S: {len(s_genes)}/{len(S_GENES)} | G2M: {len(g2m_genes)}/{len(G2M_GENES)}")

if len(s_genes) >= 10 and len(g2m_genes) >= 10:
    all_cc = list(set(s_genes + g2m_genes))
    cc_mask = adata.var_names.isin(all_cc)

    adata_temp = sc.AnnData(
        X=adata.layers['counts'][:, cc_mask].copy(),
        var=adata.var.iloc[cc_mask].copy()
    )
    adata_temp.var_names = adata.var_names[cc_mask]

    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)

    s_filt = [g for g in s_genes if g in adata_temp.var_names]
    g2m_filt = [g for g in g2m_genes if g in adata_temp.var_names]

    sc.tl.score_genes_cell_cycle(adata_temp, s_genes=s_filt, g2m_genes=g2m_filt)

    adata.obs['S_score'] = adata_temp.obs['S_score'].values.astype(np.float32)
    adata.obs['G2M_score'] = adata_temp.obs['G2M_score'].values.astype(np.float32)
    adata.obs['phase'] = adata_temp.obs['phase'].values

    del adata_temp
    gc.collect()

    print(f"  S: {adata.obs['S_score'].mean():.3f} | G2M: {adata.obs['G2M_score'].mean():.3f}")
    print(f"  MEMORY: only {len(all_cc)} genes used")
else:
    adata.obs['S_score'] = np.float32(0.0)
    adata.obs['G2M_score'] = np.float32(0.0)
    adata.obs['phase'] = 'G1'


# %%
# Create tissue if missing
if TISSUE_KEY not in adata.obs.columns:
    adata.obs[TISSUE_KEY] = 'myeloid'
    print(f"  Created '{TISSUE_KEY}' column")

# Force dtypes for categorical covariates
print("\n" + "-"*70)
print("Enforcing covariate dtypes...")
print("-"*70)

for cat in [BATCH_KEY, TISSUE_KEY]:
    adata.obs[cat] = adata.obs[cat].astype('category')
    print(f"  {cat}: category ({adata.obs[cat].nunique()} levels)")

for cont in ['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']:
    adata.obs[cont] = adata.obs[cont].astype(np.float32)
    print(f"  {cont}: float32")

# NOTE: categorical_covariate_levels saved AFTER batch filtering (Step 4)
# to avoid stale levels from removed small batches [P1-4 FIX v2.4]

# Validation
missing = [c for c in REQUIRED_COVARIATES if c not in adata.obs.columns]
if missing:
    raise ValueError(f"Missing: {missing}")

print("\nAll covariates ready")


# %% [markdown]
# ---
# ## 2. Standardize L3 Labels
# 

# %%
print("="*80)
print("2. STANDARDIZING L3 LABELS")
print("="*80)

def standardize_l3(adata, l2=CELLTYPE_L2_COL, l3=CELLTYPE_L3_COL):
    sample = str(adata.obs[l3].iloc[0])
    if '_c' in sample:
        print("  Already hierarchical")
        return adata

    print("  Standardizing...")
    adata.obs['cell_type_L3_original'] = adata.obs[l3].copy()
    adata.obs['subcluster_id'] = pd.to_numeric(adata.obs[l3].astype(str), errors='coerce').astype('Int64')

    if isinstance(adata.obs[l3].dtype, pd.CategoricalDtype):
        adata.obs[l3] = adata.obs[l3].astype('object')

    l2_vals = adata.obs[l2].astype(str)
    sid_vals = adata.obs['subcluster_id']
    mask = l2_vals.notna() & sid_vals.notna()
    adata.obs.loc[mask, l3] = l2_vals[mask] + '_c' + sid_vals[mask].astype('Int64').astype(str)

    all_labels = pd.unique(adata.obs[l3].astype(str)).tolist()
    hierarchical = [l for l in all_labels if '_c' in l]
    other = [l for l in all_labels if '_c' not in l]

    if hierarchical:
        df = pd.DataFrame({
            'label': hierarchical,
            'l2': [l.rsplit('_c', 1)[0] for l in hierarchical],
            'cid': [int(l.rsplit('_c', 1)[1]) if l.rsplit('_c', 1)[1].isdigit() else 999 for l in hierarchical]
        })
        df = df.sort_values(['l2', 'cid'])
        ordered = df['label'].tolist() + other
    else:
        ordered = all_labels

    adata.obs[l3] = pd.Categorical(adata.obs[l3].astype(str), categories=ordered, ordered=True)
    print(f"  Standardized: {len(hierarchical)} hierarchical + {len(other)} other")
    return adata

adata = standardize_l3(adata)


# %% [markdown]
# ---
# ## 3. Apply Refined Annotation
# 

# %%
print("="*80)
print("3. APPLY REFINED ANNOTATION")
print("="*80)

CELLTYPE_REANNOTATION = {
    # Alveolar (3+1 -> 2)
    'Alveolar macrophages_c0': 'Resident Alveolar macrophages',
    'Alveolar macrophages_c1': 'Resident Alveolar macrophages',
    'Alveolar macrophages_c2': 'Resident Alveolar macrophages',
    'Alveolar macrophages_c3': 'Resting Alveolar macrophages',

    # Classical mono (3 -> 3, detect neutrophils)
    'Classical monocytes_c0': 'Neutrophils',
    'Classical monocytes_c1': 'Typical Classical monocytes',
    'Classical monocytes_c2': 'Inflammatory Classical monocytes',

    # DC2 (2+2 -> 2)
    'DC2_c0': 'Conventional cDC2',
    'DC2_c1': 'Langerhans-like cDC2',
    'DC_c0': 'Conventional cDC2',
    'DC_c1': 'Conventional cDC2',

    # pDC (2 -> 1)
    'pDC_c0': 'pDC',
    'pDC_c1': 'pDC',

    # Mast (2 -> 1)
    'Mast cells_c0': 'Mast cells',
    'Mast cells_c1': 'Mast cells',

    # Macrophages (3 -> 3)
    'Macrophages_c0': 'Inflammatory Interstitial macrophages',
    'Macrophages_c1': 'M2-like Interstitial macrophages',
    'Macrophages_c2': 'Atypically activated Interstitial macrophages',

    # Intestinal (5 -> 2+1 preserved)
    'Intestinal macrophages_c0': 'CD163L1+ Interstitial macrophages',
    'Intestinal macrophages_c1': 'Low-quality Interstitial macrophages',
    'Intestinal macrophages_c2': 'CD163L1+ Interstitial macrophages',
    'Intestinal macrophages_c3': 'Immunoregulatory Interstitial macrophages',
    'Intestinal macrophages_c4': 'CD163L1+ Interstitial macrophages',
}

adata.obs[CELLTYPE_L3_REFINED_COL] = adata.obs[CELLTYPE_L3_COL].map(CELLTYPE_REANNOTATION).fillna(adata.obs[CELLTYPE_L3_COL])

# ========================================================================
# P0-4 FIX (v2.4): .map().fillna() strips categorical dtype.
# Restore it immediately so downstream .cat.categories works.
# ========================================================================
adata.obs[CELLTYPE_L3_REFINED_COL] = pd.Categorical(adata.obs[CELLTYPE_L3_REFINED_COL])

print(f"Labels: {adata.obs[CELLTYPE_L3_COL].nunique()} -> {adata.obs[CELLTYPE_L3_REFINED_COL].nunique()}")
print("Refined annotation applied")


# %% [markdown]
# ---
# ## 4. Preprocessing
# 

# %%
print("="*80)
print("4. PREPROCESSING")
print("="*80)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Filter small batches
batch_counts = adata.obs[BATCH_KEY].value_counts()
small = batch_counts[batch_counts < 3].index
if len(small) > 0:
    print(f"  Removing {len(small)} small batches")
    adata = adata[~adata.obs[BATCH_KEY].isin(small)].copy()

# ========================================================================
# P1-4 FIX (v2.4): Save categorical_covariate_levels AFTER batch filtering
# so levels don't include removed small batches (stale levels).
# ========================================================================
adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].cat.remove_unused_categories()
adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].cat.remove_unused_categories()

adata.uns['categorical_covariate_levels'] = {
    BATCH_KEY: adata.obs[BATCH_KEY].cat.categories.tolist(),
    TISSUE_KEY: adata.obs[TISSUE_KEY].cat.categories.tolist(),
}
print(f"  Levels saved to .uns (post-filter): {BATCH_KEY}={adata.obs[BATCH_KEY].nunique()}, "
      f"{TISSUE_KEY}={adata.obs[TISSUE_KEY].nunique()}")

print(f"After filtering: {adata.n_obs:,} cells | {adata.n_vars:,} genes")


# %%
# HVG selection
print(f"Selecting {N_HVG} HVGs...")
try:
    sc.pp.highly_variable_genes(
        adata, layer='counts', n_top_genes=N_HVG,
        batch_key=BATCH_KEY, flavor='seurat_v3', subset=False
    )
    method = "batch-aware (seurat_v3)"
except Exception:
    sc.pp.highly_variable_genes(
        adata, layer='counts', n_top_genes=N_HVG,
        flavor='seurat_v3', subset=False
    )
    method = "non-batch-aware (seurat_v3)"

print(f"  {adata.var['highly_variable'].sum()} HVGs ({method})")


# %%
# Create log1p layer for .raw
print("Creating log1p layer...")

if 'log1p' not in adata.layers:
    counts = adata.layers['counts']

    if issparse(counts):
        total = np.array(counts.sum(axis=1)).flatten().astype(np.float64)
        from scipy.sparse import diags as sp_diags
        scale = sp_diags(1e4 / np.where(total > 0, total, 1))
        log1p_data = (scale @ counts).tocsr()
        log1p_data.data = np.log1p(log1p_data.data).astype(np.float32)
        adata.layers['log1p'] = log1p_data
    else:
        total = counts.sum(axis=1, keepdims=True)
        cpm = counts * (1e4 / np.where(total > 0, total, 1))
        adata.layers['log1p'] = np.log1p(cpm).astype(np.float32)

    del counts
    gc.collect()
    print("  log1p created")

# Save full gene space to .raw (shared memory)
print("Preserving full matrix to .raw...")
adata.raw = sc.AnnData(
    X=adata.layers['log1p'],    # shared memory, no .copy()
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"  {adata.raw.n_vars} genes (log1p, optimal for dotplot)")


# %%
# Subset to HVG
adata = adata[:, adata.var['highly_variable']].copy()
print(f"  Subset to {adata.n_vars} HVGs")


# %% [markdown]
# ---
# ## 5. Create adata_scvi
# 

# %%
# ========================================================================
# P2-1 FIX (v2.4): adata is already HVG-subset, so build adata_scvi
# directly from its counts layer (no redundant hvg_idx indexing).
# ========================================================================
print("="*80)
print("5. CREATE adata_scvi")
print("="*80)

adata_scvi = sc.AnnData(
    X=adata.layers['counts'].copy(),
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"  adata_scvi: {adata_scvi.shape}")


# %% [markdown]
# ---
# ## 6. Train scVI
# 

# %%
print("="*80)
print("6. TRAIN scVI")
print("="*80)

print("Setting up scVI with covariates...")
scvi.model.SCVI.setup_anndata(
    adata_scvi,
    layer=None,
    batch_key=BATCH_KEY,
    categorical_covariate_keys=[TISSUE_KEY],
    continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
)
print("  Setup complete")


# %%
print(f"Creating scVI (n_latent={N_LATENT}, encode_covariates={ENCODE_COVARIATES})...")

try:
    vae = scvi.model.SCVI(
        adata_scvi, n_latent=N_LATENT, n_layers=N_LAYERS, dropout_rate=DROPOUT_RATE,
        gene_likelihood="nb", encode_covariates=ENCODE_COVARIATES,
        use_layer_norm=USE_LAYER_NORM, use_batch_norm=USE_BATCH_NORM
    )
    print("  Created with scVI 1.1+ parameters")
except TypeError:
    print("  New parameters not supported, using fallback...")
    vae = scvi.model.SCVI(
        adata_scvi, n_latent=N_LATENT, n_layers=N_LAYERS,
        dropout_rate=DROPOUT_RATE, gene_likelihood="nb"
    )
    print("  Created with fallback parameters")


# %%
print("Training scVI...")
vae.train(max_epochs=MAX_EPOCHS_SCVI, batch_size=BATCH_SIZE, early_stopping=True)
print("  Training complete")

vae.save(OUTPUT_DIR / "scvi_model", overwrite=True)
print("  Model saved")


# %%
# ========================================================================
# P1-2 FIX (v2.4): scArches dry-run validation — scVI
# QRM v2.8 Section 9.2: verify model artifact is query-compatible.
# ========================================================================
print("[scArches validation] Dry-run prepare_query_anndata on scVI...")
try:
    adata_check = adata_scvi[:100].copy()
    scvi.model.SCVI.prepare_query_anndata(adata_check, str(OUTPUT_DIR / "scvi_model"))
    del adata_check; gc.collect()
    print("  scVI passes scArches dry-run")
except Exception as e:
    print(f"  [WARN] scVI dry-run FAILED: {e}")


# %%
adata.obsm['X_scvi'] = vae.get_latent_representation()
print("  Latent saved to .obsm['X_scvi']")


# %% [markdown]
# ---
# ## 7. Train scANVI
# 

# %%
print("="*80)
print("7. TRAIN scANVI")
print("="*80)

print("Preparing labels for scANVI...")
adata_scvi.obs[CELLTYPE_L3_REFINED_COL] = adata_scvi.obs[CELLTYPE_L3_REFINED_COL].astype("category")

if "Unknown" not in adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories:
    adata_scvi.obs[CELLTYPE_L3_REFINED_COL] = \
        adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.add_categories(["Unknown"])
    print("  Added 'Unknown' to categories")
else:
    print("  'Unknown' already in categories")

print(f"  Total: {len(adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories)} categories")


# %%
print("Training scANVI...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae, adata=adata_scvi, labels_key=CELLTYPE_L3_REFINED_COL, unlabeled_category="Unknown"
)
lvae.train(max_epochs=MAX_EPOCHS_SCANVI, batch_size=BATCH_SIZE, n_samples_per_label=2000)
print("  Training complete")

lvae.save(OUTPUT_DIR / "scanvi_model", overwrite=True)
print("  Model saved")


# %%
# ========================================================================
# P1-2 FIX (v2.4): scArches dry-run validation — scANVI
# ========================================================================
print("[scArches validation] Dry-run prepare_query_anndata on scANVI...")
try:
    adata_check = adata_scvi[:100].copy()
    scvi.model.SCANVI.prepare_query_anndata(adata_check, str(OUTPUT_DIR / "scanvi_model"))
    del adata_check; gc.collect()
    print("  scANVI passes scArches dry-run")
except Exception as e:
    print(f"  [WARN] scANVI dry-run FAILED: {e}")


# %%
# Get predictions + probabilities
print("Generating predictions...")

pred_probs = lvae.predict(soft=True)
adata.obsm['scanvi_proba'] = pred_probs.astype(np.float32)

try:
    label_map = lvae.adata_manager.get_state_registry('labels')['categorical_mapping']
    label_order = [label_map[i] for i in range(len(label_map))]
except Exception:
    label_order = list(adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories)

adata.uns['scanvi_label_order'] = label_order

scanvi_pred = lvae.predict()
adata.obs['scanvi_predictions'] = pd.Categorical(scanvi_pred, categories=label_order)
adata.obs['scanvi_confidence'] = pred_probs.max(axis=1).astype(np.float32)

print(f"  Probs: {pred_probs.shape}")
print(f"  Mean confidence: {adata.obs['scanvi_confidence'].mean():.3f}")

# scANVI latent
adata.obsm['X_scanvi'] = lvae.get_latent_representation()


# %% [markdown]
# ---
# ## 8. Dimensionality Reduction
# 

# %%
print("="*80)
print("8. DIMENSIONALITY REDUCTION")
print("="*80)

sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30)
sc.tl.umap(adata, min_dist=0.5)
print("  UMAP on scVI latent")


# %% [markdown]
# ---
# ## 9. Validation Markers
# 

# %%
VALIDATION_MARKERS = {
    'Resident_Alveolar_Mac': ['MARCO','FABP4','PPARG','SIGLEC1','CHIT1','MSR1'],
    'Resting_Alveolar_Mac': ['IL10','TGFB1','MERTK','MRC1','APOE','TREM2'],
    'Neutrophils': ['FCGR3B','CSF3R','MPO','ELANE','S100A8','S100A9'],
    'Typical_Classical_Mono': ['FCN1','S100A8','S100A9','LILRB1','CTSS'],
    'Inflammatory_Classical_Mono': ['IL1B','CXCL8','PTX3','NFKBIA','TNF'],
    'Conventional_cDC2': ['CD1C','FCER1A','CLEC10A','IRF4','HLA-DRA'],
    'Langerhans_like_cDC2': ['CD207','CD1A','CCR6','EPCAM','CXCL14'],
    'pDC': ['CLEC4C','IL3RA','GZMB','TCF4','IRF7'],
    'Mast_cells': ['TPSB2','TPSAB1','CPA3','MS4A2','KIT'],
    'Inflammatory_Interstitial_Mac': ['CCL20','PTX3','TIMP1','IL1B','SPP1'],
    'M2_like_Interstitial_Mac': ['C1QA','C1QB','C1QC','FOLR2','MRC1'],
    'CD163L1_Interstitial_Mac': ['CD163L1','SELENOP','LYVE1','F13A1','CXCL12'],
    'Immunoregulatory_Interstitial_Mac': ['IL10','TGFB1','MERTK','GAS6','STAB1'],
}

all_markers = [g for genes in VALIDATION_MARKERS.values() for g in genes]
available = {k: [g for g in v if g in adata.raw.var_names] for k, v in VALIDATION_MARKERS.items()}
available_flat = [g for genes in available.values() for g in genes]

print(f"  Markers: {len(available_flat)}/{len(all_markers)} available")


# %% [markdown]
# ---
# ## 10. Comprehensive Dotplot
# 

# %%
celltype_order = list(adata.obs[CELLTYPE_L3_REFINED_COL].cat.categories)

n_genes = len(available_flat)
n_types = len(celltype_order)
fig_w = max(30, n_genes * 0.3)
fig_h = max(12, n_types * 0.7)

print(f"Generating dotplot ({n_genes} markers x {n_types} types)...")

dp = sc.pl.dotplot(
    adata, var_names=available_flat, groupby=CELLTYPE_L3_REFINED_COL,
    standard_scale='var', use_raw=True, show=False, figsize=(fig_w, fig_h)
)

dp.fig.suptitle('Myeloid Validation - Refined Annotation', fontsize=18, y=0.998, weight='bold')
dp.fig.tight_layout()

dp.savefig(OUTPUT_DIR / 'dotplot_VALIDATION.pdf', bbox_inches='tight')
dp.savefig(OUTPUT_DIR / 'dotplot_VALIDATION.png', dpi=300, bbox_inches='tight')

print("  Dotplot saved")


# %% [markdown]
# ---
# ## 11. UMAP Visualizations
# 

# %%
# ========================================================================
# P2-2 FIX (v2.4): Enable rasterization for large UMAP plots.
# QRM v3.6: use sc.settings.vector_friendly = True; do NOT pass
# rasterized=True to sc.pl.embedding (causes TypeError).
# ========================================================================
sc.settings.vector_friendly = True

fig, axes = plt.subplots(1, 3, figsize=(24, 7))

sc.pl.umap(adata, color=CELLTYPE_L3_REFINED_COL, ax=axes[0], show=False, title='Refined Annotation')
sc.pl.umap(adata, color='scanvi_predictions', ax=axes[1], show=False, title='scANVI Predictions')
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[2], show=False, title='Batch')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'umap_overview.pdf')
plt.savefig(OUTPUT_DIR / 'umap_overview.png', dpi=300)
plt.close()

print("  UMAPs saved")


# %% [markdown]
# ---
# ## 12. Save scArches Files
# 

# %%
print("="*80)
print("12. SAVING scArches-REQUIRED FILES")
print("="*80)

# HVG list
hvg_genes = adata_scvi.var_names.tolist()
pd.Series(hvg_genes).to_csv(OUTPUT_DIR / 'hvg_genes.txt', index=False, header=False)
print(f"  HVG list: {len(hvg_genes)} genes")

# ========================================================================
# P1-3 FIX (v2.4): Save var_names.csv in BOTH model directories
# QRM v2.8 Section 9.3: SCANVI dir also needs var_names.csv for
# direct load_query_data without going through scVI path.
# ========================================================================
for model_dir in ['scvi_model', 'scanvi_model']:
    pd.Series(hvg_genes).to_csv(OUTPUT_DIR / model_dir / 'var_names.csv', index=False, header=False)
    print(f"  var_names.csv -> {model_dir}/")

# Covariate config
import json
config = {
    'required_covariates': REQUIRED_COVARIATES,
    'categorical_covariates': {
        BATCH_KEY: {
            'levels': adata.obs[BATCH_KEY].cat.categories.tolist(),
            'n_levels': adata.obs[BATCH_KEY].nunique()
        },
        TISSUE_KEY: {
            'levels': adata.obs[TISSUE_KEY].cat.categories.tolist(),
            'n_levels': adata.obs[TISSUE_KEY].nunique()
        }
    },
    'continuous_covariates': {
        cov: {
            'min': float(adata.obs[cov].min()),
            'max': float(adata.obs[cov].max()),
            'mean': float(adata.obs[cov].mean()),
            'std': float(adata.obs[cov].std())
        }
        for cov in ['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
    },
    'model_params': {'n_hvg': len(hvg_genes), 'n_latent': N_LATENT, 'n_layers': N_LAYERS},
    'stress_score_normalization': 'raw_mean_log1p_cpm',  # [v2.4] document no [0,1]
}

with open(OUTPUT_DIR / 'covariate_config.json', 'w') as f:
    json.dump(config, f, indent=2)

print("  Covariate config saved")
print("  CRITICAL for scArches query mapping")


# %% [markdown]
# ---
# ## 13. Save Final H5AD
# 

# %%
output_h5ad = OUTPUT_DIR / 'adata_myeloid_refined_FINAL.h5ad'
adata.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)

print("="*80)
print("PIPELINE COMPLETE")
print("="*80)
print(f"\nMain output: {output_h5ad.name}")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,} (HVG) | {adata.raw.n_vars:,} (total in .raw)")
print(f"  Layers: {list(adata.layers.keys())}")
print(f"  Embeddings: {list(adata.obsm.keys())}")
print(f"\nModels:")
print(f"  scvi_model/ (for scArches)")
print(f"  scanvi_model/ (for scArches)")
print(f"\nscArches files:")
print(f"  hvg_genes.txt")
print(f"  var_names.csv (in BOTH model dirs)")
print(f"  covariate_config.json")
print(f"\nVisualizations:")
print(f"  dotplot_VALIDATION.pdf/png")
print(f"  umap_overview.pdf/png")
print("="*80)
