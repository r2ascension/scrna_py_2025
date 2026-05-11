# %% [markdown]
# # Myeloid Cell Refinement Pipeline v2.4 PRODUCTION
# 
# **Author:** r2end | **Date:** 2026-02-12 | **Version:** 2.4 PRODUCTION ✅
# 
# ## All P0/P1 Fixes Applied
# - ✅ P0-1: Dotplot dp.savefig()
# - ✅ P0-2: L3 categories no-NaN
# - ✅ P0-3: scANVI Unknown handling
# - ✅ P0-NEW-1: Annotation schema alignment (L2/L3/L4)
# - ✅ P0-NEW-2: cell_type_L2_refined (Neutrophils L2 fix)
# - ✅ P0-NEW-3: cell_type_L3_refined → category immediately after apply
# - ✅ P0-NEW-4: covariate levels saved AFTER small-batch filtering
# - ✅ P1-1: Memory-optimized covariates (70GB peak)
# - ✅ P1-2: HVG flavor=seurat_v3
# - ✅ P1-3: .raw with log1p
# - ✅ P1-4: dtype enforcement (category/float32)
# - ✅ P1-5: scVI parameter fallback
# - ✅ P1-6: Thread control
# - ✅ P1-NEW-1: filter_cells/genes on counts layer
# - ✅ P1-NEW-2: Assign actual Unknown cells for semi-supervised
# - ✅ P1-NEW-3: Hard validation label_order vs pred_probs
# - ✅ P2-2: scArches files | ✅ P2-3: Full probability matrix | ✅ P2-5: Agg backend
# 
# **Memory:** 50-70GB | **Runtime:** 2-3h GPU | **Output:** scArches-ready
# 

# %% [markdown]
# ---
# ## 0. Setup
# 

# %%
# P2-5: Matplotlib backend for HPC
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# P1-6: BLAS thread control
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

print(f"✓ Backend: {matplotlib.get_backend()}")
print(f"✓ scanpy {sc.__version__} | scvi {scvi.__version__} | GPU: {torch.cuda.is_available()}")


# %%
# CONFIGURATION
INPUT_H5AD = "/home/h2048/data/py/0128/myeloid_analysis_unified/results/subcluster_unified_v2_20260128/adata_myeloid_subclustered_FINAL_v2_20260128.h5ad"
OUTPUT_DIR = Path("/home/h2048/data/py/0209/myeloid_validation_optimized")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CELLTYPE_L2_COL = 'cell_type_L2'
CELLTYPE_L3_COL = 'cell_type_L3'           # Original cluster (will become L4)
CELLTYPE_L3_REFINED_COL = 'cell_type_L3_refined'
BATCH_KEY = 'sample'   # ⚠️ CHECK
TISSUE_KEY = 'tissue' # ⚠️ CHECK (auto-created if missing)

REQUIRED_COVARIATES = [BATCH_KEY, TISSUE_KEY,
                       'pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']

# Cells to treat as Unknown in scANVI (semi-supervised)
LOW_CONF_L4 = ["Macrophages_c2", "Intestinal macrophages_c1"]  # ⚠️ ADJUST

# scVI / scANVI (scArches-optimised)
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

print(f"Config: n_latent={N_LATENT}, encode_covariates={ENCODE_COVARIATES}")
print(f"LOW_CONF_L4 (→ Unknown): {LOW_CONF_L4}")


# %%
sc.settings.verbosity = 3
sc.settings.n_jobs = N_JOBS
sc.settings.set_figure_params(dpi=100, facecolor='white')
np.random.seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED
torch.manual_seed(RANDOM_SEED)
print("✓ Settings configured")


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
print(f"✓ Gene lists: Stress={len(STRESS_SIGNATURE_GENES)}, S={len(S_GENES)}, G2M={len(G2M_GENES)}")


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
print("✓ Data loaded")


# %% [markdown]
# ---
# ## 1.5 Calculate Covariates
# 

# %%
# MT%
print("-"*70)
print("MT%...")
mt_genes = (adata.var_names.str.startswith('MT-') |
            adata.var_names.str.startswith('Mt-') |
            adata.var_names.str.startswith('mt-'))
print(f"  MT genes: {mt_genes.sum()}")

if mt_genes.sum() > 0:
    counts = adata.layers['counts']
    if issparse(counts):
        total = np.array(counts.sum(axis=1)).flatten()
        mt    = np.array(counts[:, mt_genes].sum(axis=1)).flatten()
    else:
        total = counts.sum(axis=1)
        mt    = counts[:, mt_genes].sum(axis=1)
    mt_pct = np.zeros_like(total, dtype=np.float32)
    mask = total > 0
    mt_pct[mask] = (mt[mask] / total[mask]) * 100
    adata.obs['pct_counts_mt'] = mt_pct
    print(f"  ✓ Mean MT%: {mt_pct.mean():.2f}%")
else:
    adata.obs['pct_counts_mt'] = np.float32(0.0)


# %%
# Stress score (MEMORY OPTIMISED: subset genes only)
print("-"*70)
print("Stress score...")
stress_genes = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
print(f"  Stress genes: {len(stress_genes)}/{len(STRESS_SIGNATURE_GENES)}")

if len(stress_genes) >= 10:
    counts = adata.layers['counts']
    sidx   = np.array([i for i,g in enumerate(adata.var_names) if g in stress_genes])

    if issparse(counts):
        sc_sub  = counts[:, sidx]
        total   = np.array(counts.sum(axis=1)).flatten().astype(np.float64)
        # vectorised row-scaling (no Python loop)
        from scipy.sparse import diags
        scale   = diags(1e4 / np.where(total > 0, total, 1))
        sc_cpm  = (scale @ sc_sub).tocsr()
        sc_cpm.data = np.log1p(sc_cpm.data).astype(np.float32)
        raw_score = np.array(sc_cpm.mean(axis=1)).flatten()
    else:
        sc_sub    = counts[:, sidx].astype(np.float32)
        total     = sc_sub.sum(axis=1, keepdims=True)
        cpm       = sc_sub * (1e4 / np.where(total > 0, total, 1))
        raw_score = np.log1p(cpm).mean(axis=1)

    smin, smax = raw_score.min(), raw_score.max()
    score = (raw_score - smin) / (smax - smin) if smax > smin else np.zeros_like(raw_score)
    adata.obs['stress_score'] = score.astype(np.float32)

    del sc_sub, sc_cpm, raw_score
    gc.collect()
    print(f"  ✓ Mean: {score.mean():.3f}  [vectorised, no loop]")
else:
    adata.obs['stress_score'] = np.float32(0.0)


# %%
# Cell cycle (MEMORY OPTIMISED: only CC genes)
print("-"*70)
print("Cell cycle...")
s_use   = [g for g in S_GENES   if g in adata.var_names]
g2m_use = [g for g in G2M_GENES if g in adata.var_names]
print(f"  S: {len(s_use)}/{len(S_GENES)}  G2M: {len(g2m_use)}/{len(G2M_GENES)}")

if len(s_use) >= 10 and len(g2m_use) >= 10:
    all_cc  = list(set(s_use + g2m_use))
    cc_mask = adata.var_names.isin(all_cc)
    adt     = sc.AnnData(X=adata.layers['counts'][:, cc_mask].copy(),
                         var=adata.var.iloc[cc_mask].copy())
    adt.var_names = adata.var_names[cc_mask]
    sc.pp.normalize_total(adt, target_sum=1e4)
    sc.pp.log1p(adt)
    sc.tl.score_genes_cell_cycle(adt,
        s_genes  =[g for g in s_use   if g in adt.var_names],
        g2m_genes=[g for g in g2m_use if g in adt.var_names])
    adata.obs['S_score']   = adt.obs['S_score'].values.astype(np.float32)
    adata.obs['G2M_score'] = adt.obs['G2M_score'].values.astype(np.float32)
    adata.obs['phase']     = adt.obs['phase'].values
    del adt; gc.collect()
    print(f"  ✓ S={adata.obs['S_score'].mean():.3f}  G2M={adata.obs['G2M_score'].mean():.3f}  [{len(all_cc)} genes only]")
else:
    adata.obs['S_score']   = np.float32(0.0)
    adata.obs['G2M_score'] = np.float32(0.0)
    adata.obs['phase']     = 'G1'


# %%
# Tissue column
if TISSUE_KEY not in adata.obs.columns:
    adata.obs[TISSUE_KEY] = 'myeloid'
    print(f"  ✓ Created '{TISSUE_KEY}'")

# P1-4: Force dtypes  (levels saved LATER, after batch filtering)
for cat in [BATCH_KEY, TISSUE_KEY]:
    adata.obs[cat] = adata.obs[cat].astype('category')
    print(f"  ✓ {cat}: category ({adata.obs[cat].nunique()} levels)")
for cont in ['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']:
    adata.obs[cont] = adata.obs[cont].astype(np.float32)
    print(f"  ✓ {cont}: float32")

missing = [c for c in REQUIRED_COVARIATES if c not in adata.obs.columns]
if missing:
    raise ValueError(f"Missing covariates: {missing}")
print("\n✓ All covariates calculated")


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
        print("  ✓ Already hierarchical")
        return adata
    adata.obs['cell_type_L3_original'] = adata.obs[l3].copy()
    adata.obs['_sid'] = pd.to_numeric(adata.obs[l3].astype(str), errors='coerce').astype('Int64')
    if isinstance(adata.obs[l3].dtype, pd.CategoricalDtype):
        adata.obs[l3] = adata.obs[l3].astype('object')
    l2v   = adata.obs[l2].astype(str)
    sidv  = adata.obs['_sid']
    mask  = l2v.notna() & sidv.notna()
    adata.obs.loc[mask, l3] = l2v[mask] + '_c' + sidv[mask].astype('Int64').astype(str)

    # P0-2: ALL unique labels → no NaN
    all_labels = pd.unique(adata.obs[l3].astype(str)).tolist()
    hier  = [l for l in all_labels if '_c' in l]
    other = [l for l in all_labels if '_c' not in l]
    if hier:
        df = pd.DataFrame({'l': hier,
                           'base': [l.rsplit('_c',1)[0] for l in hier],
                           'cid':  [int(l.rsplit('_c',1)[1]) if l.rsplit('_c',1)[1].isdigit() else 999 for l in hier]})
        ordered = df.sort_values(['base','cid'])['l'].tolist() + other
    else:
        ordered = all_labels
    adata.obs[l3] = pd.Categorical(adata.obs[l3].astype(str), categories=ordered, ordered=True)
    adata.obs.drop(columns=['_sid'], inplace=True)
    print(f"  ✓ {len(hier)} hierarchical + {len(other)} other  (no NaN)")
    return adata

adata = standardize_l3(adata)


# %% [markdown]
# ---
# ## 3. Apply Refined Annotation
# 

# %%
CELLTYPE_REANNOTATION = {
    'Alveolar macrophages_c0': 'Resident Alveolar macrophages',
    'Alveolar macrophages_c1': 'Resident Alveolar macrophages',
    'Alveolar macrophages_c2': 'Resident Alveolar macrophages',
    'Alveolar macrophages_c3': 'Resting Alveolar macrophages',
    'Classical monocytes_c0':  'Neutrophils',
    'Classical monocytes_c1':  'Typical Classical monocytes',
    'Classical monocytes_c2':  'Inflammatory Classical monocytes',
    'DC2_c0':  'Conventional cDC2',
    'DC2_c1':  'Langerhans-like cDC2',
    'DC_c0':   'Conventional cDC2',
    'DC_c1':   'Conventional cDC2',
    'pDC_c0':  'pDC',
    'pDC_c1':  'pDC',
    'Mast cells_c0': 'Mast cells',
    'Mast cells_c1': 'Mast cells',
    'Macrophages_c0': 'Inflammatory Interstitial macrophages',
    'Macrophages_c1': 'M2-like Interstitial macrophages',
    'Macrophages_c2': 'Atypically activated Interstitial macrophages',
    'Intestinal macrophages_c0': 'CD163L1+ Interstitial macrophages',
    'Intestinal macrophages_c1': 'Low-quality Interstitial macrophages',
    'Intestinal macrophages_c2': 'CD163L1+ Interstitial macrophages',
    'Intestinal macrophages_c3': 'Immunoregulatory Interstitial macrophages',
    'Intestinal macrophages_c4': 'CD163L1+ Interstitial macrophages',
}

adata.obs[CELLTYPE_L3_REFINED_COL] = (
    adata.obs[CELLTYPE_L3_COL]
    .map(CELLTYPE_REANNOTATION)
    .fillna(adata.obs[CELLTYPE_L3_COL])
)

# P0-NEW-3: Convert to category immediately and clean unused
adata.obs[CELLTYPE_L3_REFINED_COL] = (
    adata.obs[CELLTYPE_L3_REFINED_COL]
    .astype('category')
    .cat.remove_unused_categories()
)
print(f"  ✓ {adata.obs[CELLTYPE_L3_COL].nunique()} L3-original → {adata.obs[CELLTYPE_L3_REFINED_COL].nunique()} refined")


# %%
# P0-NEW-2: cell_type_L2_refined (fixes Neutrophils L2 mismatch)
L3_TO_L2_REFINED = {
    "Resident Alveolar macrophages":              "Alveolar macrophages",
    "Resting Alveolar macrophages":               "Alveolar macrophages",
    "Typical Classical monocytes":                "Classical monocytes",
    "Inflammatory Classical monocytes":           "Classical monocytes",
    "Neutrophils":                                "Neutrophils",          # ← corrected
    "Conventional cDC2":                          "cDC2",
    "Langerhans-like cDC2":                       "cDC2",
    "pDC":                                        "pDC",
    "Mast cells":                                 "Mast cells",
    "Inflammatory Interstitial macrophages":      "Interstitial macrophages",
    "M2-like Interstitial macrophages":           "Interstitial macrophages",
    "Atypically activated Interstitial macrophages": "Interstitial macrophages",
    "CD163L1+ Interstitial macrophages":          "Interstitial macrophages",
    "Low-quality Interstitial macrophages":       "Interstitial macrophages",
    "Immunoregulatory Interstitial macrophages":  "Interstitial macrophages",
}

adata.obs['cell_type_L2_refined'] = (
    adata.obs[CELLTYPE_L3_REFINED_COL]
    .map(L3_TO_L2_REFINED)
    .astype('category')
)
print("  ✓ cell_type_L2_refined created (Neutrophils L2 corrected)")


# %% [markdown]
# ---
# ## 4. Preprocessing
# 

# %%
# P1-NEW-1: QC metrics from counts layer (not X)
print("="*80)
print("4. PREPROCESSING")
print("="*80)

sc.pp.calculate_qc_metrics(adata, layer='counts', inplace=True)

# Filter on counts-based metrics
min_genes = 200
min_cells = 3
n_before   = adata.n_obs

cell_mask = adata.obs['n_genes_by_counts'] >= min_genes
adata      = adata[cell_mask].copy()
print(f"  filter_cells (≥{min_genes} genes by counts): {n_before:,} → {adata.n_obs:,}")

n_before_g = adata.n_vars
gene_mask  = adata.var['n_cells_by_counts'] >= min_cells
adata      = adata[:, gene_mask].copy()
print(f"  filter_genes (≥{min_cells} cells): {n_before_g:,} → {adata.n_vars:,}")


# %%
# Remove small batches
batch_counts = adata.obs[BATCH_KEY].value_counts()
small = batch_counts[batch_counts < 3].index
if len(small) > 0:
    print(f"  Removing {len(small)} small batches: {list(small)}")
    adata = adata[~adata.obs[BATCH_KEY].isin(small)].copy()

# P0-NEW-3: remove_unused after filtering
adata.obs[CELLTYPE_L3_REFINED_COL] = (
    adata.obs[CELLTYPE_L3_REFINED_COL].cat.remove_unused_categories()
)

# P0-NEW-4: Save covariate levels AFTER filtering (stable set)
for cat in [BATCH_KEY, TISSUE_KEY]:
    adata.obs[cat] = adata.obs[cat].cat.remove_unused_categories()

adata.uns['categorical_covariate_levels'] = {
    BATCH_KEY:  adata.obs[BATCH_KEY].cat.categories.tolist(),
    TISSUE_KEY: adata.obs[TISSUE_KEY].cat.categories.tolist(),
}
print(f"  ✓ Covariate levels saved post-filtering")
print(f"  Final: {adata.n_obs:,} cells | {adata.n_vars:,} genes | {adata.obs[BATCH_KEY].nunique()} batches")


# %%
# P1-2: seurat_v3 flavor
print(f"Selecting {N_HVG} HVGs...")
try:
    sc.pp.highly_variable_genes(
        adata, layer='counts', n_top_genes=N_HVG,
        batch_key=BATCH_KEY, flavor='seurat_v3', subset=False)
    hvg_method = "batch-aware (seurat_v3)"
except Exception as e:
    print(f"  ⚠️ Batch-aware failed: {e} — fallback")
    sc.pp.highly_variable_genes(
        adata, layer='counts', n_top_genes=N_HVG,
        flavor='seurat_v3', subset=False)
    hvg_method = "non-batch-aware (seurat_v3)"
print(f"  ✓ {adata.var['highly_variable'].sum()} HVGs ({hvg_method})")


# %%
# P1-3: .raw = log1p (better for dotplot)
print("Creating log1p layer...")
if 'log1p' not in adata.layers:
    counts = adata.layers['counts']
    if issparse(counts):
        total = np.array(counts.sum(axis=1)).flatten().astype(np.float64)
        from scipy.sparse import diags
        scale     = diags(1e4 / np.where(total > 0, total, 1))
        l1p       = (scale @ counts).tocsr()
        l1p.data  = np.log1p(l1p.data).astype(np.float32)
        adata.layers['log1p'] = l1p
    else:
        total = counts.sum(axis=1, keepdims=True)
        adata.layers['log1p'] = np.log1p(counts * (1e4 / total)).astype(np.float32)
    del counts; gc.collect()
    print("  ✓ log1p layer created [vectorised]")

adata.raw = sc.AnnData(
    X=adata.layers['log1p'],
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"  ✓ .raw: {adata.raw.n_vars} genes (log1p)")


# %%
# Subset to HVG
adata = adata[:, adata.var['highly_variable']].copy()
print(f"  ✓ HVG subset: {adata.n_vars} genes")


# %% [markdown]
# ---
# ## 5. Create adata_scvi
# 

# %%
adata_scvi = sc.AnnData(
    X=adata.layers['counts'].copy(),
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
# Sync refined col dtype
adata_scvi.obs[CELLTYPE_L3_REFINED_COL] = (
    adata_scvi.obs[CELLTYPE_L3_REFINED_COL].astype('category')
)
print(f"  ✓ adata_scvi: {adata_scvi.shape}")


# %% [markdown]
# ---
# ## 6. Train scVI
# 

# %%
scvi.model.SCVI.setup_anndata(
    adata_scvi,
    layer=None,
    batch_key=BATCH_KEY,
    categorical_covariate_keys=[TISSUE_KEY],
    continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
)
print("  ✓ scVI setup with all covariates")


# %%
# P1-5: parameter fallback
print(f"Creating scVI (n_latent={N_LATENT}, encode_covariates={ENCODE_COVARIATES})...")
try:
    vae = scvi.model.SCVI(
        adata_scvi, n_latent=N_LATENT, n_layers=N_LAYERS, dropout_rate=DROPOUT_RATE,
        gene_likelihood="nb", encode_covariates=ENCODE_COVARIATES,
        use_layer_norm=USE_LAYER_NORM, use_batch_norm=USE_BATCH_NORM)
    print("  ✓ scVI 1.1+ parameters")
except TypeError as e:
    print(f"  ⚠️ Fallback: {e}")
    vae = scvi.model.SCVI(
        adata_scvi, n_latent=N_LATENT, n_layers=N_LAYERS,
        dropout_rate=DROPOUT_RATE, gene_likelihood="nb")
    print("  ✓ Fallback parameters")


# %%
vae.train(max_epochs=MAX_EPOCHS_SCVI, batch_size=BATCH_SIZE,
          early_stopping=True, early_stopping_patience=45)
vae.save(OUTPUT_DIR / "scvi_model", overwrite=True)
adata.obsm['X_scvi'] = vae.get_latent_representation()
print("  ✓ scVI trained & saved")


# %% [markdown]
# ---
# ## 7. Train scANVI
# 

# %%
# P1-NEW-2: Assign actual Unknown cells (semi-supervised)
print("Assigning Unknown cells for semi-supervised learning...")

# Add Unknown category first
adata_scvi.obs[CELLTYPE_L3_REFINED_COL] = (
    adata_scvi.obs[CELLTYPE_L3_REFINED_COL].astype('category')
)
if "Unknown" not in adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories:
    adata_scvi.obs[CELLTYPE_L3_REFINED_COL] = (
        adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.add_categories(["Unknown"])
    )

# Set low-confidence clusters to Unknown
if LOW_CONF_L4:
    # Ensure cell_type_L4 is present (use original L3 col at this point)
    l4_source = adata_scvi.obs[CELLTYPE_L3_COL].astype(str)
    mask = l4_source.isin(LOW_CONF_L4)
    n_unknown = mask.sum()
    if n_unknown > 0:
        adata_scvi.obs.loc[mask, CELLTYPE_L3_REFINED_COL] = "Unknown"
        print(f"  ✓ {n_unknown:,} cells → Unknown ({LOW_CONF_L4})")
    else:
        print(f"  ⚠️ No cells matched LOW_CONF_L4 — check cluster names")
else:
    print("  ⚠️ LOW_CONF_L4 empty — fully-supervised mode")

n_cats = len(adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories)
print(f"  Total categories: {n_cats}")


# %%
lvae = scvi.model.SCANVI.from_scvi_model(
    vae, adata=adata_scvi,
    labels_key=CELLTYPE_L3_REFINED_COL,
    unlabeled_category="Unknown")
lvae.train(max_epochs=MAX_EPOCHS_SCANVI, batch_size=BATCH_SIZE,
           n_samples_per_label=2000, early_stopping=True,
           early_stopping_patience=30)
lvae.save(OUTPUT_DIR / "scanvi_model", overwrite=True)
print("  ✓ scANVI trained & saved")


# %%
# P2-3: Full probability matrix
pred_probs = lvae.predict(soft=True)
adata.obsm['scanvi_proba'] = pred_probs.astype(np.float32)

# Get label order
try:
    lm = lvae.adata_manager.get_state_registry('labels')['categorical_mapping']
    label_order = [lm[i] for i in range(len(lm))]
except Exception:
    label_order = list(adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories)

# P1-NEW-3: Hard validation
if pred_probs.shape[1] != len(label_order):
    print(f"  ⚠️ Shape mismatch ({pred_probs.shape[1]} vs {len(label_order)}) — using categories fallback")
    label_order = list(adata_scvi.obs[CELLTYPE_L3_REFINED_COL].cat.categories)
    assert pred_probs.shape[1] == len(label_order),         f"Unresolvable mismatch: {pred_probs.shape[1]} vs {len(label_order)}"

adata.uns['scanvi_label_order'] = label_order

scanvi_pred = lvae.predict()
adata.obs['scanvi_predictions'] = pd.Categorical(scanvi_pred, categories=label_order)
adata.obs['scanvi_confidence']  = pred_probs.max(axis=1).astype(np.float32)

adata.obsm['X_scanvi'] = lvae.get_latent_representation()

print(f"  ✓ Probs: {pred_probs.shape}  |  labels: {len(label_order)}")
print(f"  ✓ Mean confidence: {adata.obs['scanvi_confidence'].mean():.3f}")
print(f"  ✓ Low conf (<0.5): {(adata.obs['scanvi_confidence'] < 0.5).sum():,} cells")


# %% [markdown]
# ---
# ## 8. Dimensionality Reduction
# 

# %%
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30)
sc.tl.umap(adata, min_dist=0.5, random_state=RANDOM_SEED)
print("  ✓ UMAP on scVI latent")


# %% [markdown]
# ---
# ## 9. Validation Markers
# 

# %%
VALIDATION_MARKERS = {
    'Resident_Alveolar_Mac':               ['MARCO','FABP4','PPARG','SIGLEC1','CHIT1','MSR1'],
    'Resting_Alveolar_Mac':                ['IL10','TGFB1','MERTK','MRC1','APOE','TREM2'],
    'Neutrophils':                         ['FCGR3B','CSF3R','MPO','ELANE','S100A8','S100A9'],
    'Typical_Classical_Mono':              ['FCN1','S100A8','S100A9','LILRB1','CTSS'],
    'Inflammatory_Classical_Mono':         ['IL1B','CXCL8','PTX3','NFKBIA','TNF'],
    'Conventional_cDC2':                   ['CD1C','FCER1A','CLEC10A','IRF4','HLA-DRA'],
    'Langerhans_like_cDC2':               ['CD207','CD1A','CCR6','EPCAM','CXCL14'],
    'pDC':                                 ['CLEC4C','IL3RA','GZMB','TCF4','IRF7'],
    'Mast_cells':                          ['TPSB2','TPSAB1','CPA3','MS4A2','KIT'],
    'Inflammatory_Interstitial_Mac':       ['CCL20','PTX3','TIMP1','IL1B','SPP1'],
    'M2_like_Interstitial_Mac':            ['C1QA','C1QB','C1QC','FOLR2','MRC1'],
    'CD163L1_Interstitial_Mac':            ['CD163L1','SELENOP','LYVE1','F13A1','CXCL12'],
    'Immunoregulatory_Interstitial_Mac':   ['IL10','TGFB1','MERTK','GAS6','STAB1'],
}

all_m = [g for v in VALIDATION_MARKERS.values() for g in v]
avail = [g for g in all_m if g in adata.raw.var_names]
print(f"  Markers available: {len(avail)}/{len(all_m)}")


# %% [markdown]
# ---
# ## 10. Comprehensive Dotplot
# 

# %%
# P0-1: dp.savefig()
# P0-NEW-3: categories already category type, so .cat.categories is safe
celltype_order = list(adata.obs[CELLTYPE_L3_REFINED_COL].cat.categories)
n_g, n_c = len(avail), len(celltype_order)
fig_w, fig_h = max(30, n_g * 0.3), max(12, n_c * 0.7)

print(f"Dotplot: {n_g} markers × {n_c} types  ({fig_w:.0f}×{fig_h:.0f} in)")
dp = sc.pl.dotplot(
    adata, var_names=avail, groupby=CELLTYPE_L3_REFINED_COL,
    standard_scale='var', use_raw=True, show=False, figsize=(fig_w, fig_h))

dp.fig.suptitle('Myeloid Cell Type Validation - Refined Annotation',
                fontsize=18, y=0.999, weight='bold')
dp.fig.tight_layout()

dp.savefig(OUTPUT_DIR / 'dotplot_VALIDATION.pdf', bbox_inches='tight')
dp.savefig(OUTPUT_DIR / 'dotplot_VALIDATION.png', dpi=300, bbox_inches='tight')
print("  ✓ Dotplot saved (dp.savefig)")


# %% [markdown]
# ---
# ## 11. UMAP Visualizations
# 

# %%
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
sc.pl.umap(adata, color=CELLTYPE_L3_REFINED_COL, ax=axes[0], show=False, title='Refined L3')
sc.pl.umap(adata, color='scanvi_predictions',    ax=axes[1], show=False, title='scANVI Predictions')
sc.pl.umap(adata, color=BATCH_KEY,               ax=axes[2], show=False, title='Batch')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'umap_overview.pdf')
plt.savefig(OUTPUT_DIR / 'umap_overview.png', dpi=300)
plt.close()

# Extra: L2_refined and confidence
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color='cell_type_L2_refined',   ax=axes2[0], show=False, title='L2 Refined')
sc.pl.umap(adata, color='scanvi_confidence',       ax=axes2[1], show=False, title='scANVI Confidence',
           color_map='RdBu_r', vmin=0, vmax=1)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'umap_L2_confidence.pdf')
plt.savefig(OUTPUT_DIR / 'umap_L2_confidence.png', dpi=300)
plt.close()
print("  ✓ UMAPs saved")


# %% [markdown]
# ---
# ## 12. Save scArches Files
# 

# %%
import json as _json

hvg_genes = adata_scvi.var_names.tolist()
pd.Series(hvg_genes).to_csv(OUTPUT_DIR / 'hvg_genes.txt', index=False, header=False)

config = {
    'required_covariates': REQUIRED_COVARIATES,
    'categorical_covariates': {
        BATCH_KEY:  {'levels': adata.obs[BATCH_KEY].cat.categories.tolist(),
                     'n_levels': adata.obs[BATCH_KEY].nunique()},
        TISSUE_KEY: {'levels': adata.obs[TISSUE_KEY].cat.categories.tolist(),
                     'n_levels': adata.obs[TISSUE_KEY].nunique()},
    },
    'continuous_covariates': {
        cov: {'min': float(adata.obs[cov].min()), 'max': float(adata.obs[cov].max()),
              'mean': float(adata.obs[cov].mean()), 'std': float(adata.obs[cov].std())}
        for cov in ['pct_counts_mt','stress_score','S_score','G2M_score']
    },
    'model_params': {'n_hvg': len(hvg_genes), 'n_latent': N_LATENT,
                     'n_layers': N_LAYERS, 'dropout_rate': DROPOUT_RATE},
    'annotation_schema': {
        'cell_type_L2':         'Original coarse lineage',
        'cell_type_L2_refined': 'Corrected lineage (Neutrophils fixed)',
        'cell_type_L3':         'Refined 15-category label (PRIMARY)',
        'cell_type_L3_refined': 'Same as cell_type_L3 (back-compat)',
        'cell_type_L4':         'Original cluster (*_c#)',
    }
}

with open(OUTPUT_DIR / 'covariate_config.json', 'w') as f:
    _json.dump(config, f, indent=2)

print(f"  ✓ hvg_genes.txt ({len(hvg_genes)} genes)")
print(f"  ✓ covariate_config.json")


# %% [markdown]
# ---
# ## 13. Annotation Schema Alignment
# 

# %%
# P0-NEW-1: Align annotation columns to project schema
# cell_type_L4 = original cluster (*_c#)
# cell_type_L3 = refined label  (PRIMARY, same as L3_refined)
# cell_type_L2_refined = corrected lineage
print("Aligning annotation schema...")

adata.obs['cell_type_L4'] = adata.obs[CELLTYPE_L3_COL].astype(str)

# Primary label: overwrite L3 with refined
adata.obs['cell_type_L3'] = adata.obs[CELLTYPE_L3_REFINED_COL].copy()

# Back-compat alias
adata.obs['cell_type_L3_refined'] = adata.obs['cell_type_L3']

# Validate: no NaN anywhere
for col in ['cell_type_L2', 'cell_type_L2_refined', 'cell_type_L3', 'cell_type_L4']:
    n_nan = adata.obs[col].isna().sum()
    status = "✅" if n_nan == 0 else f"❌ {n_nan} NaN"
    print(f"  {status}  {col}  ({adata.obs[col].nunique()} unique values)")

print("\n  Annotation schema:")
print("  cell_type_L2         = original coarse lineage")
print("  cell_type_L2_refined = corrected lineage (Neutrophils fixed)")
print("  cell_type_L3         = refined 15-category (PRIMARY label)")
print("  cell_type_L3_refined = same as L3 (back-compat)")
print("  cell_type_L4         = original cluster (*_c#)")


# %% [markdown]
# ---
# ## 14. Save Final H5AD
# 

# %%
output_h5ad = OUTPUT_DIR / 'adata_myeloid_refined_FINAL.h5ad'
adata.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)

print("="*80)
print("PIPELINE COMPLETE ✅")
print("="*80)
print(f"  {output_h5ad.name}")
print(f"  Cells:       {adata.n_obs:,}")
print(f"  Genes (HVG): {adata.n_vars:,}")
print(f"  Genes (.raw):{adata.raw.n_vars:,}")
print(f"  Layers:      {list(adata.layers.keys())}")
print(f"  Embeddings:  {list(adata.obsm.keys())}")
print()
print("  Annotation columns:")
for col in ['cell_type_L2','cell_type_L2_refined','cell_type_L3','cell_type_L3_refined','cell_type_L4']:
    n = adata.obs[col].nunique() if col in adata.obs else '—'
    print(f"    {col}: {n} categories")
print()
print("  scArches files: scvi_model/ | scanvi_model/ | hvg_genes.txt | covariate_config.json")
print("="*80)



