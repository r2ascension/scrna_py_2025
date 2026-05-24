# %% [markdown]
# # Epithelial scArches Query Mapping Pipeline v1.2 — PRODUCTION
# #
# # **Changelog from v1.1:**
# #   [NEW] Default reference updated to v2.7-HOTFIX
# #   [NEW] Reference model paths loaded from REF/manifest.json when available
# #
# # **Preserved from v1.0/v1.1:**
# #   [P0] Actually call prepare_query_anndata per-model; three separate adata objects
# #   [P1] No manual SCANVI.setup_anndata() — let registry drive everything
# #   [P1] query-only UMAP explicitly labelled; reference-space mapping attempted
# #   [P1] adata_main preserved (full gene) for marker visualization
# #   [P2] stress_score matches v2.6 exactly (log1p mean, no [0,1] normalization)
# #   [P2] major / fine adata objects are independent (no registry cross-contamination)
# #   [P2] Full posterior probability matrices saved (obsm + uns label order)
# #
# # **Reference models (v2.7-HOTFIX):**
# #   scVI  : SELF/scvi_model        (trained on ALL epithelial cells)
# #   major : REF/scanvi_major_model  (AT/Alveolar removed)
# #   fine  : REF/scanvi_fine_model   (AT/Alveolar removed)
# #
# # **Author:** r2end  **Date:** 2026-03-31  **Version:** 1.2

# %% [markdown]
# ## Environment Setup

# %%
import os
os.environ["OMP_NUM_THREADS"]        = "8"
os.environ["OPENBLAS_NUM_THREADS"]   = "8"
os.environ["MKL_NUM_THREADS"]        = "8"
os.environ["VECLIB_MAXIMUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"]    = "8"

import gc
import json
import time
import warnings
import traceback
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
from scipy import sparse

from query_covid_filter_helper_20260523 import (
    DEFAULT_COVID_PATTERNS,
    build_covid_filter_debug_frame,
    build_removed_kept_value_counts,
    summarize_covid_filter_debug_frame,
)

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def _env_flag(name: str, default: bool = False) -> bool:
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def _env_csv(name: str, default: tuple[str, ...]) -> tuple[str, ...]:
    value = os.getenv(name)
    if value is None:
        return default
    parsed = tuple(item.strip() for item in value.split(',') if item.strip())
    return parsed or default

# =============================================================================
# PATHS
# =============================================================================
QUERY_H5AD = os.getenv(
    "EPITHELIAL_QUERY_H5AD",
    "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/epithelial_cells.h5ad",
)

V27_BASE         = Path("/home/h2048/data/py/0317/epithelial_v2_7_HOTFIX")
REF_MANIFEST     = V27_BASE / "REF" / "manifest.json"
SCVI_REF_DIR     = str(V27_BASE / "SELF" / "scvi_model")
HVG_FILE         = str(V27_BASE / "SELF" / "scvi_model" / "hvg_genes.txt")
SCANVI_MAJOR_REF = str(V27_BASE / "REF"  / "scanvi_major_model")
SCANVI_FINE_REF  = str(V27_BASE / "REF"  / "scanvi_fine_model")

TODAY    = datetime.now().strftime('%m%d')
OUT_BASE = Path(
    os.getenv(
        "EPITHELIAL_SCARCHES_OUT_DIR",
        f"/home/h2048/data/py/{TODAY}/epithelial_scarches_query_v1_2",
    )
)
FIG_DIR  = OUT_BASE / "figures"
CKPT_DIR = OUT_BASE / "checkpoints"
for d in [OUT_BASE, FIG_DIR, CKPT_DIR,
          OUT_BASE / "query_scvi_model",
          OUT_BASE / "query_scanvi_major_model",
          OUT_BASE / "query_scanvi_fine_model"]:
    d.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = FIG_DIR

# =============================================================================
# KEYS — must match v2.6 training exactly
# =============================================================================
BATCH_KEY          = 'sample'
UNLABELED_CATEGORY = 'Unknown'
REQUIRED_COVS      = ['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
 
# scArches fine-tuning
SCVI_MAX_EPOCHS   = 40
SCANVI_MAX_EPOCHS = 20
BATCH_SIZE        = 512
WEIGHT_DECAY      = 0.0   # REQUIRED for scArches surgery

# Checkpoints / IO
CHECKPOINT_COMPRESSION = None  # faster than gzip for large intermediate states
FINAL_H5AD_COMPRESSION = None  # practical for very large mapped query objects
 
# Visualization
FIGURE_DPI        = 300
FIGURE_FORMAT     = "pdf"
UMAP_SIZE         = 3
UMAP_ALPHA        = 0.6
LOW_CONF_THRESHOLD = 0.5
MIN_GENE_OVERLAP_PCT = 70.0

# Optional metadata-driven COVID filtering for extended/query atlas reruns.
FILTER_COVID_RELATED = _env_flag("EPITHELIAL_FILTER_COVID_RELATED", False)
COVID_FILTER_PRIMARY_COLS = _env_csv("EPITHELIAL_FILTER_PRIMARY_COLS", ("disease",))
COVID_FILTER_SECONDARY_COLS = _env_csv("EPITHELIAL_FILTER_SECONDARY_COLS", ("COVID_status",))
COVID_FILTER_PATTERNS = _env_csv("EPITHELIAL_FILTER_PATTERNS", DEFAULT_COVID_PATTERNS)
COVID_FILTER_INCLUDE_DATASET_FALLBACK = _env_flag(
    "EPITHELIAL_FILTER_INCLUDE_DATASET_FALLBACK",
    False,
)
 
# =============================================================================
# COVARIATE SIGNATURES — byte-for-byte identical to v2.6 training code
# =============================================================================
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
    "XPOT","YIF1A","YWHAZ","ZBTB17",
]
S_GENES = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP',
    'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2',
    'RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
    'USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8',
]
G2M_GENES = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
    'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2',
    'CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B',
    'HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2',
    'DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN',
    'LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA',
]
 
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
 
PIPELINE_START = time.time()
print("=" * 80)
print("Epithelial scArches Query Mapping Pipeline v1.2 — PRODUCTION")
print("=" * 80)
print(f"Query:         {QUERY_H5AD}")
print(f"Manifest:      {REF_MANIFEST}")
print(f"scVI ref:      {SCVI_REF_DIR}")
print(f"scANVI major:  {SCANVI_MAJOR_REF}")
print(f"scANVI fine:   {SCANVI_FINE_REF}")
print(f"Output:        {OUT_BASE}")
print(f"Device:        {DEVICE}")
print(f"Date:          {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"COVID filter:  {FILTER_COVID_RELATED}")
print("=" * 80)
 
 
# %% [markdown]
# ## Step 1: Validate Reference Model Paths + Manifest + Load HVG List
 
# %%
print("\n" + "=" * 80)
print("STEP 1: Validating Paths")
print("=" * 80)

if REF_MANIFEST.exists():
    ref_manifest = json.loads(REF_MANIFEST.read_text(encoding="utf-8"))
    SCVI_REF_DIR     = ref_manifest.get("scvi_model_path", SCVI_REF_DIR)
    HVG_FILE         = ref_manifest.get("hvg_file", HVG_FILE)
    SCANVI_MAJOR_REF = ref_manifest.get("scanvi_major_model", SCANVI_MAJOR_REF)
    SCANVI_FINE_REF  = ref_manifest.get("scanvi_fine_model", SCANVI_FINE_REF)
    print(f"  [OK] REF manifest loaded: {REF_MANIFEST}")
else:
    print(f"  [WARN] REF manifest not found, using hardcoded v2.7 paths: {REF_MANIFEST}")
 
for label, path in [
    ("REF manifest",      str(REF_MANIFEST)),
    ("scVI model dir",    SCVI_REF_DIR),
    ("HVG file",          HVG_FILE),
    ("scANVI major dir",  SCANVI_MAJOR_REF),
    ("scANVI fine dir",   SCANVI_FINE_REF),
    ("Query h5ad",        QUERY_H5AD),
]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"[ERROR] Required path not found: {path}")
    print(f"  [OK] {label}: {path}")
 
hvg_genes_ref = pd.read_csv(HVG_FILE, header=None)[0].astype(str).tolist()
print(f"\n[OK] Reference HVG: {len(hvg_genes_ref)} genes")
print(f"     Sample: {hvg_genes_ref[:5]}")
 
 
# %% [markdown]
# ## Step 2: Load Query → adata_main  (full gene, kept for visualization)
 
# %%
print("\n" + "=" * 80)
print("STEP 2: Loading Query Data into adata_main")
print("=" * 80)
 
t0 = time.time()
adata_main = sc.read_h5ad(QUERY_H5AD)
print(f"[OK] Loaded in {time.time()-t0:.1f}s  —  shape: {adata_main.shape}")
print(f"  Layers:  {list(adata_main.layers.keys())}")
print(f"  obs cols: {list(adata_main.obs.columns)}")
 
# Ensure counts layer on adata_main
if 'counts' not in adata_main.layers:
    x_max = float(adata_main.X.max() if not sparse.issparse(adata_main.X)
                  else adata_main.X.max())
    if x_max > 20:
        print("[INFO] 'counts' layer absent; treating .X as raw counts")
        adata_main.layers['counts'] = (
            sparse.csr_matrix(adata_main.X) if not sparse.issparse(adata_main.X)
            else adata_main.X.copy()
        )
    else:
        raise ValueError(
            f"No 'counts' layer found and .X.max()={x_max:.2f} — not raw counts."
        )
else:
    if not sparse.issparse(adata_main.layers['counts']):
        adata_main.layers['counts'] = sparse.csr_matrix(adata_main.layers['counts'])
 
# Ensure BATCH_KEY
if BATCH_KEY not in adata_main.obs.columns:
    for alt in ['batch', 'Sample', 'Batch']:
        if alt in adata_main.obs.columns:
            adata_main.obs[BATCH_KEY] = adata_main.obs[alt].astype(str)
            print(f"[INFO] BATCH_KEY remapped from '{alt}'")
            break
    else:
        raise ValueError(
            f"BATCH_KEY='{BATCH_KEY}' not found. "
            f"Available: {list(adata_main.obs.columns)}"
        )
adata_main.obs[BATCH_KEY] = adata_main.obs[BATCH_KEY].astype(str)
 
# Gene overlap check (on full gene space, before any alignment)
ref_genes_set   = set(hvg_genes_ref)
query_genes_set = set(adata_main.var_names.astype(str))
overlap_pct = len(ref_genes_set & query_genes_set) / len(ref_genes_set) * 100
missing_in_query = ref_genes_set - query_genes_set
print(f"\n  Reference HVGs: {len(ref_genes_set):,}")
print(f"  Query genes:    {len(query_genes_set):,}")
print(f"  Overlap:        {overlap_pct:.1f}%  "
      f"(missing from query: {len(missing_in_query):,})")
if overlap_pct < MIN_GENE_OVERLAP_PCT:
    raise ValueError(
        f"[ERROR] Gene overlap {overlap_pct:.1f}% < {MIN_GENE_OVERLAP_PCT}%.\n"
        f"  Check gene naming. Ref sample: {hvg_genes_ref[:5]}\n"
        f"  Query sample: {list(query_genes_set)[:5]}"
    )
print(f"[OK] Gene overlap sufficient")
print(f"\n  Batches: {adata_main.obs[BATCH_KEY].nunique()}")
print(adata_main.obs[BATCH_KEY].value_counts().to_string())

covid_filter_summary = None
covid_removed_n = 0
covid_removed_fraction = 0.0
if FILTER_COVID_RELATED:
    print("\n" + "=" * 80)
    print("STEP 2A: Filtering COVID-related Cells from Query")
    print("=" * 80)
    filter_debug = build_covid_filter_debug_frame(
        adata_main.obs.copy(),
        primary_cols=COVID_FILTER_PRIMARY_COLS,
        secondary_cols=COVID_FILTER_SECONDARY_COLS,
        dataset_cols=("dataset",),
        patterns=COVID_FILTER_PATTERNS,
        allow_dataset_fallback=COVID_FILTER_INCLUDE_DATASET_FALLBACK,
    )
    filter_mask = filter_debug["remove_covid_related"].to_numpy(dtype=bool)
    covid_removed_n = int(filter_mask.sum())
    covid_removed_fraction = float(covid_removed_n / max(adata_main.n_obs, 1))
    covid_filter_summary = summarize_covid_filter_debug_frame(filter_debug)

    audit_path = OUT_BASE / "covid_filter_audit.tsv.gz"
    filter_debug.to_csv(audit_path, sep='\t', compression='gzip')
    summary_path = OUT_BASE / "covid_filter_summary.tsv"
    covid_filter_summary.to_csv(summary_path, sep='\t', index=False)
    print(f"[OK] COVID filter audit: {audit_path}")
    print(f"[OK] COVID filter summary: {summary_path}")

    value_count_tables = build_removed_kept_value_counts(
        adata_main.obs,
        filter_mask,
        columns=("disease", "COVID_status", "dataset", "condition", "tissue"),
    )
    for col, table in value_count_tables.items():
        table_path = OUT_BASE / f"covid_filter_{col}_counts.tsv"
        table.to_csv(table_path, sep='\t')
        print(f"[OK] {table_path.name}")

    print(f"[INFO] Removing {covid_removed_n:,} / {adata_main.n_obs:,} cells ({covid_removed_fraction:.1%})")
    print(covid_filter_summary.to_string(index=False))

    if covid_removed_n > 0:
        adata_main = adata_main[~filter_mask].copy()
        print(f"[OK] Query shape after COVID filter: {adata_main.shape}")
    else:
        print("[INFO] No COVID-related cells matched the configured metadata filter")
 
 
# %% [markdown]
# ## Step 3: Calculate Covariates on FULL Gene Space  [P0: before any gene reordering]
# #
# # Must match epithelial reference training exactly:
# #   pct_counts_mt : MT-gene percentage from raw counts
# #   stress_score  : mean of log1p(CPM) across stress genes (no [0,1] rescaling)
# #   S_score / G2M_score : scanpy cell cycle scoring on log1p(CPM)
 
# %%
print("\n" + "=" * 80)
print("STEP 3: Calculating Covariates (full gene space)")
print("=" * 80)
 
# Temporary log1p-normalized object on full genes (matches v2.6 _tmp pattern exactly)
_tmp = sc.AnnData(
    X=adata_main.layers['counts'].copy(),
    var=adata_main.var.copy(),
    obs=adata_main.obs[[BATCH_KEY]].copy(),
)
sc.pp.normalize_total(_tmp, target_sum=1e4)
sc.pp.log1p(_tmp)
 
# --- MT% from raw counts (handle zero-total cells) ---
mt_mask = adata_main.var_names.str.startswith('MT-')
counts_mat = adata_main.layers['counts']
total_counts = np.asarray(counts_mat.sum(axis=1)).ravel().astype(float)
mt_counts = (
    np.asarray(counts_mat[:, mt_mask].sum(axis=1)).ravel().astype(float)
    if mt_mask.sum() > 0 else np.zeros(adata_main.n_obs, dtype=float)
)
mt_pct = np.zeros_like(total_counts)
nz = total_counts > 0
mt_pct[nz] = mt_counts[nz] / total_counts[nz] * 100
adata_main.obs['pct_counts_mt'] = mt_pct.astype(np.float32)
print(f"[OK] pct_counts_mt  mean={mt_pct.mean():.2f}%  MT genes={mt_mask.sum()}")
 
# --- Stress score: mean of log1p(CPM) across available stress genes ---
# Identical to v2.6 training (no [0,1] normalization)
sg_avail = [g for g in STRESS_SIGNATURE_GENES if g in _tmp.var_names]
if len(sg_avail) > 10:
    sg_idx = [list(_tmp.var_names).index(g) for g in sg_avail]
    sg_x   = _tmp.X[:, sg_idx]
    stress = (np.asarray(sg_x.mean(axis=1)).ravel()
              if sparse.issparse(sg_x) else sg_x.mean(axis=1))
else:
    stress = np.zeros(adata_main.n_obs, dtype=float)
    print(f"[WARN] Only {len(sg_avail)} stress genes found; stress_score=0")
adata_main.obs['stress_score'] = stress.astype(np.float32)
print(f"[OK] stress_score  mean={stress.mean():.4f}  n_genes={len(sg_avail)}")
 
# --- Cell cycle scoring: identical to v2.6 ---
s_avail   = [g for g in S_GENES   if g in _tmp.var_names]
g2m_avail = [g for g in G2M_GENES if g in _tmp.var_names]
if len(s_avail) > 10 and len(g2m_avail) > 10:
    sc.tl.score_genes_cell_cycle(_tmp, s_genes=s_avail, g2m_genes=g2m_avail)
    adata_main.obs['S_score']   = _tmp.obs['S_score'].values.astype(np.float32)
    adata_main.obs['G2M_score'] = _tmp.obs['G2M_score'].values.astype(np.float32)
    adata_main.obs['phase']     = _tmp.obs['phase'].values
    print(f"[OK] Cell cycle: {adata_main.obs['phase'].value_counts().to_dict()}")
else:
    adata_main.obs['S_score']   = np.float32(0.0)
    adata_main.obs['G2M_score'] = np.float32(0.0)
    adata_main.obs['phase']     = 'G1'
    print(f"[WARN] CC genes insufficient (S={len(s_avail)}, G2M={len(g2m_avail)}); set to 0")
 
del _tmp; gc.collect()
 
for cov in REQUIRED_COVS:
    assert cov in adata_main.obs.columns, f"[P0 ASSERT] Covariate '{cov}' missing!"
print(f"[OK] All covariates present: {REQUIRED_COVS}")
 
 
# %% [markdown]
# ## Step 4: Set .raw on adata_main  (full-gene counts, for marker visualization)
# #
# # MUST happen after covariates are calculated but before any gene subsetting.
# # adata_main itself is never subsetted — it stays at full gene width.
 
# %%
print("\n" + "=" * 80)
print("STEP 4: Preserving Full-Gene .raw on adata_main")
print("=" * 80)
 
# Shared memory: no .copy() on X — zero extra cost
adata_main.raw = sc.AnnData(
    X=adata_main.layers['counts'],    # shared memory, 0 extra cost
    obs=adata_main.obs.copy(),
    var=adata_main.var.copy(),
)
print(f"[OK] adata_main.raw set: {adata_main.raw.shape}")
print(f"[OK] adata_main kept at full gene width: {adata_main.shape}")
 
# Checkpoint
adata_main.write_h5ad(
    str(CKPT_DIR / "adata_main_with_covariates.h5ad"),
    compression=CHECKPOINT_COMPRESSION,
)
print(f"[OK] Checkpoint: adata_main_with_covariates.h5ad")
 
 
# %% [markdown]
# ## Step 5: Branch A — scVI prepare_query_anndata + load_query_data
# #
# # Each branch gets its own independent adata object.
# # prepare_query_anndata handles: zero-padding + gene reordering + registry alignment.
# # DO NOT call setup_anndata() manually — let the saved model registry drive everything.
 
# %%
print("\n" + "=" * 80)
print("STEP 5: Branch A — scVI Architecture Surgery")
print("=" * 80)
 
# Independent copy for scVI branch (starts at full gene width)
adata_q_scvi = sc.AnnData(
    X=adata_main.layers['counts'].copy(),
    obs=adata_main.obs.copy(),
    var=adata_main.var.copy(),
    layers={'counts': adata_main.layers['counts'].copy()},
)
print(f"[INFO] adata_q_scvi before prepare: {adata_q_scvi.shape}")
 
# prepare_query_anndata: registry-level alignment
# Handles zero-padding for missing HVG genes and reorders to reference var space
scvi.model.SCVI.prepare_query_anndata(adata_q_scvi, SCVI_REF_DIR)
print(f"[OK]  adata_q_scvi after prepare:  {adata_q_scvi.shape}")
 
# Critical post-prepare covariate check
for cov in REQUIRED_COVS:
    assert cov in adata_q_scvi.obs.columns, (
        f"[P0 ASSERT] Covariate '{cov}' lost during prepare_query_anndata (scVI)! "
        "Check that adata_main.obs had correct dtypes."
    )
print(f"[OK] Covariates intact after scVI prepare: {REQUIRED_COVS}")
 
# Architecture surgery + fine-tuning (weight_decay=0 required for scArches)
# CRITICAL: load reference WITHOUT adata — passing query adata to .load() triggers
# category mismatch for new batches not seen during training.
# load_query_data accepts a directory path directly and handles new categories internally.
query_scvi = scvi.model.SCVI.load_query_data(
    adata_q_scvi, SCVI_REF_DIR, freeze_dropout=True
)
print(f"[INFO] scVI n_latent={query_scvi.module.n_latent}")
 
query_scvi.train(
    max_epochs=SCVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    train_size=0.9,
    plan_kwargs={'weight_decay': WEIGHT_DECAY, 'lr': 1e-3},
)
query_scvi.save(str(OUT_BASE / "query_scvi_model"), overwrite=True)
 
# Write latent back to adata_main (index-aligned)
latent_scvi = query_scvi.get_latent_representation()
adata_main.obsm['X_scvi'] = pd.DataFrame(
    latent_scvi, index=adata_q_scvi.obs_names
).reindex(adata_main.obs_names).values
print(f"[OK] scVI latent saved: {adata_main.obsm['X_scvi'].shape}")
 
del query_scvi; gc.collect()
torch.cuda.empty_cache() if DEVICE == 'cuda' else None
 
 
# %% [markdown]
# ## Step 6: Branch B — scANVI Major prepare_query_anndata + load_query_data
 
# %%
print("\n" + "=" * 80)
print("STEP 6: Branch B — scANVI Major Architecture Surgery")
print("=" * 80)
 
# Independent copy for major branch — never reuses adata_q_scvi
adata_q_major = sc.AnnData(
    X=adata_main.layers['counts'].copy(),
    obs=adata_main.obs.copy(),
    var=adata_main.var.copy(),
    layers={'counts': adata_main.layers['counts'].copy()},
)
print(f"[INFO] adata_q_major before prepare: {adata_q_major.shape}")
 
# prepare_query_anndata for scANVI major — uses saved model registry
# No manual setup_anndata: the saved model defines labels_key / unlabeled_category
scvi.model.SCANVI.prepare_query_anndata(adata_q_major, SCANVI_MAJOR_REF)
print(f"[OK]  adata_q_major after prepare:  {adata_q_major.shape}")
 
for cov in REQUIRED_COVS:
    assert cov in adata_q_major.obs.columns, (
        f"[P0 ASSERT] Covariate '{cov}' lost during prepare_query_anndata (major)!"
    )
print(f"[OK] Covariates intact after scANVI_major prepare: {REQUIRED_COVS}")
 
# CRITICAL: same reason as scVI — load reference WITHOUT query adata
query_scanvi_major = scvi.model.SCANVI.load_query_data(
    adata_q_major, SCANVI_MAJOR_REF, freeze_dropout=True
)
query_scanvi_major.train(
    max_epochs=SCANVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    train_size=0.9,
    plan_kwargs={'weight_decay': WEIGHT_DECAY, 'lr': 1e-3},
)
query_scanvi_major.save(str(OUT_BASE / "query_scanvi_major_model"), overwrite=True)
 
# Predictions
latent_major = query_scanvi_major.get_latent_representation()
pred_major   = query_scanvi_major.predict()
soft_major   = query_scanvi_major.predict(soft=True)
 
# Write back to adata_main (index-aligned)
idx = adata_q_major.obs_names
adata_main.obsm['X_scanvi_major'] = pd.DataFrame(latent_major, index=idx).reindex(
    adata_main.obs_names).values
adata_main.obs['scanvi_major_pred'] = pd.Series(pred_major, index=idx).reindex(
    adata_main.obs_names)
conf_major = (soft_major.values.max(axis=1) if isinstance(soft_major, pd.DataFrame)
              else np.asarray(soft_major).max(axis=1))
adata_main.obs['scanvi_major_conf'] = pd.Series(conf_major, index=idx).reindex(
    adata_main.obs_names)
 
# Save full posterior matrix (P2: needed for ambiguous state analysis)
if isinstance(soft_major, pd.DataFrame):
    adata_main.obsm['scanvi_major_probabilities'] = soft_major.reindex(
        adata_main.obs_names).values.astype(np.float32)
    adata_main.uns['scanvi_major_celltype_order'] = list(soft_major.columns)
    print(f"[OK] Full posterior saved: shape={soft_major.shape}")
 
print(f"\n[OK] scANVI_major predictions:")
print(adata_main.obs['scanvi_major_pred'].value_counts().to_string())
print(f"  Mean conf: {adata_main.obs['scanvi_major_conf'].mean():.3f}")
 
del query_scanvi_major; gc.collect()
torch.cuda.empty_cache() if DEVICE == 'cuda' else None
 
 
# %% [markdown]
# ## Step 7: Branch C — scANVI Fine prepare_query_anndata + load_query_data
 
# %%
print("\n" + "=" * 80)
print("STEP 7: Branch C — scANVI Fine Architecture Surgery")
print("=" * 80)
 
# Independent copy for fine branch — does not share registry with major branch
adata_q_fine = sc.AnnData(
    X=adata_main.layers['counts'].copy(),
    obs=adata_main.obs.copy(),
    var=adata_main.var.copy(),
    layers={'counts': adata_main.layers['counts'].copy()},
)
print(f"[INFO] adata_q_fine before prepare: {adata_q_fine.shape}")
 
scvi.model.SCANVI.prepare_query_anndata(adata_q_fine, SCANVI_FINE_REF)
print(f"[OK]  adata_q_fine after prepare:  {adata_q_fine.shape}")
 
for cov in REQUIRED_COVS:
    assert cov in adata_q_fine.obs.columns, (
        f"[P0 ASSERT] Covariate '{cov}' lost during prepare_query_anndata (fine)!"
    )
print(f"[OK] Covariates intact after scANVI_fine prepare: {REQUIRED_COVS}")
 
# CRITICAL: same reason as scVI — load reference WITHOUT query adata
query_scanvi_fine = scvi.model.SCANVI.load_query_data(
    adata_q_fine, SCANVI_FINE_REF, freeze_dropout=True
)
query_scanvi_fine.train(
    max_epochs=SCANVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    train_size=0.9,
    plan_kwargs={'weight_decay': WEIGHT_DECAY, 'lr': 1e-3},
)
query_scanvi_fine.save(str(OUT_BASE / "query_scanvi_fine_model"), overwrite=True)
 
# Predictions
latent_fine = query_scanvi_fine.get_latent_representation()
pred_fine   = query_scanvi_fine.predict()
soft_fine   = query_scanvi_fine.predict(soft=True)
 
idx = adata_q_fine.obs_names
adata_main.obsm['X_scanvi_fine'] = pd.DataFrame(latent_fine, index=idx).reindex(
    adata_main.obs_names).values
adata_main.obs['scanvi_fine_pred'] = pd.Series(pred_fine, index=idx).reindex(
    adata_main.obs_names)
conf_fine = (soft_fine.values.max(axis=1) if isinstance(soft_fine, pd.DataFrame)
             else np.asarray(soft_fine).max(axis=1))
adata_main.obs['scanvi_fine_conf'] = pd.Series(conf_fine, index=idx).reindex(
    adata_main.obs_names)
 
# Save full posterior matrix
if isinstance(soft_fine, pd.DataFrame):
    adata_main.obsm['scanvi_fine_probabilities'] = soft_fine.reindex(
        adata_main.obs_names).values.astype(np.float32)
    adata_main.uns['scanvi_fine_celltype_order'] = list(soft_fine.columns)
    print(f"[OK] Full posterior saved: shape={soft_fine.shape}")
 
# Confidence filter: low-conf cells -> Unknown
low_conf_mask = adata_main.obs['scanvi_fine_conf'] < LOW_CONF_THRESHOLD
adata_main.obs['scanvi_fine_pred_filtered'] = adata_main.obs['scanvi_fine_pred'].copy()
adata_main.obs.loc[low_conf_mask, 'scanvi_fine_pred_filtered'] = UNLABELED_CATEGORY
 
print(f"\n[OK] scANVI_fine predictions:")
print(adata_main.obs['scanvi_fine_pred'].value_counts().to_string())
print(f"\n  Mean conf: {adata_main.obs['scanvi_fine_conf'].mean():.3f}")
print(f"  Low conf (<{LOW_CONF_THRESHOLD}): {low_conf_mask.sum():,} ({low_conf_mask.mean():.1%})")
print(f"\nFiltered predictions:")
print(adata_main.obs['scanvi_fine_pred_filtered'].value_counts().to_string())
 
del query_scanvi_fine
del adata_q_scvi, adata_q_major, adata_q_fine
gc.collect()
torch.cuda.empty_cache() if DEVICE == 'cuda' else None
 
 
# %% [markdown]
# ## Step 8: UMAP Computation
# #
# # v2.7-HOTFIX did not explicitly save a UMAP operator, so we compute query-only UMAPs.
# # These are NOT reference-space projections — they reflect the query latent geometry.
# # UMAPs are clearly labelled as "query-only" in all figure titles.
 
# %%
print("\n" + "=" * 80)
print("STEP 8: UMAP Computation (query-only, from scANVI latents)")
print("=" * 80)
 
for rep_key, nbrs_key, umap_key, label in [
    ('X_scanvi_major', 'nbrs_major', 'X_umap_major', 'major'),
    ('X_scanvi_fine',  'nbrs_fine',  'X_umap_fine',  'fine'),
]:
    sc.pp.neighbors(adata_main, use_rep=rep_key, n_neighbors=15, key_added=nbrs_key)
    sc.tl.umap(adata_main, neighbors_key=nbrs_key, min_dist=0.5)
    adata_main.obsm[umap_key] = adata_main.obsm['X_umap'].copy()
    print(f"[OK] UMAP_{label} computed")
 
sc.tl.leiden(adata_main, neighbors_key='nbrs_fine', resolution=0.5,
             key_added='leiden_fine_0.5')
print(f"[OK] Leiden (res=0.5): {adata_main.obs['leiden_fine_0.5'].nunique()} clusters")
 
# Checkpoint post-UMAP
adata_main.write_h5ad(
    str(CKPT_DIR / "adata_main_post_umap.h5ad"),
    compression=CHECKPOINT_COMPRESSION,
)
print(f"[OK] Checkpoint: adata_main_post_umap.h5ad")
 
 
# %% [markdown]
# ## Step 9: Visualizations
# #
# # All marker plots use adata_main.raw (full gene, ~53k genes) — not HVG space.
 
# %%
print("\n" + "=" * 80)
print("STEP 9: Generating Visualizations")
print("=" * 80)
 
# Marker universe: full gene set via .raw
use_raw = adata_main.raw is not None
gene_universe = adata_main.raw.var_names if use_raw else adata_main.var_names
print(f"[INFO] Gene universe for markers: {len(gene_universe):,} "
      f"({'raw' if use_raw else 'HVG only'})")
 
# ------------------------------------------------------------------
# 9-A: Overview UMAPs (query-only coordinate space)
# ------------------------------------------------------------------
print("[INFO] 9-A: Overview UMAPs...")
 
fig, axes = plt.subplots(2, 3, figsize=(24, 14))
axes_flat = axes.flatten()
 
plot_specs = [
    ('X_umap_fine', 'scanvi_fine_pred_filtered', 'scANVI Fine Pred (filtered) [query-only UMAP]',  'tab20', None),
    ('X_umap_fine', 'scanvi_major_pred',          'scANVI Major Pred [query-only UMAP]',             'Set2',  None),
    ('X_umap_fine', 'leiden_fine_0.5',            'Leiden res=0.5 [query-only UMAP]',               'tab20', None),
    ('X_umap_fine', 'scanvi_fine_conf',           'Fine Confidence',                                 'RdYlGn', (0, 1)),
    ('X_umap_fine', 'scanvi_major_conf',          'Major Confidence',                                'RdYlGn', (0, 1)),
    ('X_umap_fine', BATCH_KEY,                    'Sample (Batch)',                                  'tab20', None),
]
 
for i, (basis, color, title, cmap, vrange) in enumerate(plot_specs):
    kw = dict(
        basis=basis, color=color, ax=axes_flat[i],
        show=False, frameon=False,
        size=UMAP_SIZE, alpha=UMAP_ALPHA,
        legend_loc='right margin', legend_fontsize=7,
        title=title,
    )
    if vrange is not None:
        kw['vmin'], kw['vmax'] = vrange
        kw['cmap'] = cmap
    sc.pl.embedding(adata_main, **kw)
    # Rasterize scatter layer post-hoc (sc.pl.embedding does not accept rasterized= in all scanpy versions)
    for coll in axes_flat[i].collections:
        coll.set_rasterized(True)
 
plt.tight_layout()
out = FIG_DIR / f"01_umap_overview.{FIGURE_FORMAT}"
plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] {out.name}")
 
# ------------------------------------------------------------------
# 9-B: Confidence Distribution
# ------------------------------------------------------------------
print("[INFO] 9-B: Confidence distributions...")
 
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, col, title in [
    (axes[0], 'scanvi_major_conf', 'scANVI Major Confidence'),
    (axes[1], 'scanvi_fine_conf',  'scANVI Fine Confidence'),
]:
    vals = adata_main.obs[col]
    ax.hist(vals, bins=50, color='steelblue', edgecolor='white', linewidth=0.3)
    ax.axvline(LOW_CONF_THRESHOLD, color='red', linestyle='--', lw=1.5,
               label=f'Threshold = {LOW_CONF_THRESHOLD}')
    ax.set_xlabel('Confidence', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    pct_low = (vals < LOW_CONF_THRESHOLD).mean() * 100
    ax.text(0.05, 0.92, f'Low conf: {pct_low:.1f}%', transform=ax.transAxes,
            color='red', fontsize=11)
plt.tight_layout()
out = FIG_DIR / f"02_confidence_distribution.{FIGURE_FORMAT}"
plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] {out.name}")
 
# ------------------------------------------------------------------
# 9-C: Per-sample stacked bar
# ------------------------------------------------------------------
print("[INFO] 9-C: Per-sample proportion bar...")
 
prop_df  = pd.crosstab(adata_main.obs[BATCH_KEY],
                        adata_main.obs['scanvi_fine_pred_filtered'],
                        normalize='index') * 100
top_types = adata_main.obs['scanvi_fine_pred_filtered'].value_counts().head(15).index
plot_cols  = [c for c in top_types if c in prop_df.columns]
 
fig, ax = plt.subplots(figsize=(max(14, len(prop_df) * 0.45), 7))
prop_df[plot_cols].plot(kind='bar', stacked=True, ax=ax, width=0.8)
ax.set_ylabel('Percentage of Cells (%)', fontsize=12)
ax.set_xlabel('Sample', fontsize=12)
ax.set_title('Predicted L3 Cell Type Distribution per Sample (Top 15, low-conf filtered)',
             fontsize=13)
ax.legend(title='L3 Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.tight_layout()
out = FIG_DIR / f"03_stacked_bar_per_sample.{FIGURE_FORMAT}"
plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] {out.name}")
 
# ------------------------------------------------------------------
# 9-D: Marker Expression UMAPs — uses adata_main.raw (full gene)
# ------------------------------------------------------------------
print("[INFO] 9-D: Marker expression UMAPs (full-gene .raw)...")
 
representative_markers = {
    'Basal':     ['KRT5', 'TP63', 'KRT14', 'KRT17', 'NGFR'],
    'Ciliated':  ['FOXJ1', 'DNAH5', 'DEUP1', 'CCNO'],
    'Secretory': ['SCGB1A1', 'MUC5AC', 'MUC5B', 'DUOX2', 'SPDEF'],
    'SMG':       ['LTF', 'LYZ', 'PIGR', 'WFDC2'],
    'Special':   ['SPRR2A', 'FOXI1', 'CFTR', 'MKI67', 'VIM'],
}
 
for category, markers in representative_markers.items():
    avail = [m for m in markers if m in gene_universe]
    if not avail:
        print(f"    [WARN] No markers in gene universe for {category}")
        continue
    n_cols = min(3, len(avail))
    n_rows = int(np.ceil(len(avail) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes_flat = np.array(axes).flatten()
    for i, gene in enumerate(avail):
        sc.pl.embedding(adata_main, basis='X_umap_fine', color=gene,
                        ax=axes_flat[i], show=False, use_raw=use_raw,
                        vmax='p99', frameon=False,
                        size=UMAP_SIZE * 0.8, alpha=UMAP_ALPHA,
                        cmap='Reds', title=gene)
        for coll in axes_flat[i].collections:
            coll.set_rasterized(True)
    for i in range(len(avail), len(axes_flat)):
        axes_flat[i].axis('off')
    plt.suptitle(f'{category} Markers (full-gene .raw)', fontsize=14, y=1.02)
    plt.tight_layout()
    out = FIG_DIR / f"04_markers_{category}.{FIGURE_FORMAT}"
    plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"    [OK] {out.name}")
 
# ------------------------------------------------------------------
# 9-E: Top-2 margin plot (ambiguous state detection)
# ------------------------------------------------------------------
print("[INFO] 9-E: Top-2 margin plot...")
 
if 'scanvi_fine_probabilities' in adata_main.obsm:
    probs = adata_main.obsm['scanvi_fine_probabilities']
    sorted_probs = np.sort(probs, axis=1)[:, ::-1]
    top2_margin = sorted_probs[:, 0] - sorted_probs[:, 1]
    adata_main.obs['scanvi_fine_top2_margin'] = top2_margin.astype(np.float32)
 
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(top2_margin, bins=60, color='steelblue', edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Top-2 Probability Margin', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('scANVI Fine — Top-2 Margin (Ambiguity Score)', fontsize=13)
    ambig = (top2_margin < 0.2).mean() * 100
    ax.text(0.05, 0.92, f'Ambiguous (<0.2 margin): {ambig:.1f}%',
            transform=ax.transAxes, color='red', fontsize=11)
    plt.tight_layout()
    out = FIG_DIR / f"05_top2_margin.{FIGURE_FORMAT}"
    plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"[OK] {out.name}")
else:
    print("[WARN] posterior matrix not available; skipping top-2 margin plot")
 
 
# %% [markdown]
# ## Step 10: Save Final Output
 
# %%
print("\n" + "=" * 80)
print("STEP 10: Saving Final Output")
print("=" * 80)
 
# Ensure category dtype on all label columns before write (prevents h5ad warnings)
for col in ['scanvi_major_pred', 'scanvi_fine_pred',
            'scanvi_fine_pred_filtered']:
    if col in adata_main.obs.columns:
        adata_main.obs[col] = adata_main.obs[col].astype('category')
 
adata_main.uns['scarches_mapping'] = {
    'pipeline_version':       'epithelial_scarches_query_v1.2',
    'reference_scvi':         SCVI_REF_DIR,
    'reference_scanvi_major': SCANVI_MAJOR_REF,
    'reference_scanvi_fine':  SCANVI_FINE_REF,
    'hvg_file':               HVG_FILE,
    'n_hvg_genes':            len(hvg_genes_ref),
    'gene_overlap_pct':       round(overlap_pct, 2),
    'n_missing_genes_padded': len(missing_in_query),
    'low_conf_threshold':     LOW_CONF_THRESHOLD,
    'scvi_max_epochs':        SCVI_MAX_EPOCHS,
    'scanvi_max_epochs':      SCANVI_MAX_EPOCHS,
    'weight_decay':           WEIGHT_DECAY,
    'covid_filter_enabled':   FILTER_COVID_RELATED,
    'covid_filter_primary_cols': list(COVID_FILTER_PRIMARY_COLS),
    'covid_filter_secondary_cols': list(COVID_FILTER_SECONDARY_COLS),
    'covid_filter_patterns': list(COVID_FILTER_PATTERNS),
    'covid_filter_include_dataset_fallback': COVID_FILTER_INCLUDE_DATASET_FALLBACK,
    'covid_removed_cells':    covid_removed_n,
    'covid_removed_fraction': round(covid_removed_fraction, 6),
    'date':                   datetime.now().strftime('%Y-%m-%d'),
    'umap_note':              (
        'UMAPs are query-only (computed from query latent); '
        'not projected to reference UMAP space (v2.7-HOTFIX ref does not save UMAP operator)'
    ),
}
 
h5ad_out = OUT_BASE / "epithelial_query_mapped_v1_2.h5ad"
adata_main.write_h5ad(str(h5ad_out), compression=FINAL_H5AD_COMPRESSION)
size_gb = h5ad_out.stat().st_size / 1e9
print(f"[OK] Saved: {h5ad_out.name} ({size_gb:.2f} GB)  shape: {adata_main.shape}")
 
# Annotation CSV
adata_main.obs[[
    BATCH_KEY,
    'scanvi_major_pred',            'scanvi_major_conf',
    'scanvi_fine_pred',             'scanvi_fine_conf',
    'scanvi_fine_pred_filtered',
    'pct_counts_mt', 'stress_score',
]].to_csv(OUT_BASE / "query_cell_annotations.csv")
print("[OK] query_cell_annotations.csv")
 
# Fine statistics
fine_stats = (
    adata_main.obs['scanvi_fine_pred_filtered']
    .value_counts().rename_axis('Cell_Type').reset_index(name='Count')
)
fine_stats['Percentage'] = 100 * fine_stats['Count'] / fine_stats['Count'].sum()
conf_agg = (adata_main.obs
            .groupby('scanvi_fine_pred_filtered')['scanvi_fine_conf']
            .agg(['mean', 'std']))
fine_stats = fine_stats.merge(conf_agg, left_on='Cell_Type', right_index=True, how='left')
fine_stats.to_csv(OUT_BASE / "scanvi_fine_statistics.csv", index=False)
print("[OK] scanvi_fine_statistics.csv")
 
# Per-sample proportion table
pd.crosstab(adata_main.obs[BATCH_KEY],
            adata_main.obs['scanvi_fine_pred_filtered'],
            normalize='index').round(4).to_csv(
    OUT_BASE / "cell_type_proportions_per_sample.csv"
)
print("[OK] cell_type_proportions_per_sample.csv")
 
 
# %% [markdown]
# ## Step 11: Run Log + AnnData Structure Summary
 
# %%
print("\n" + "=" * 80)
print("STEP 11: Writing Run Log and AnnData Structure Summary")
print("=" * 80)
 
elapsed_min = (time.time() - PIPELINE_START) / 60
 
def _mat_line(name, obj):
    return (f"  {name}: type={type(obj).__name__}  "
            f"shape={getattr(obj,'shape','NA')}  "
            f"dtype={getattr(obj,'dtype','NA')}  "
            f"sparse={sparse.issparse(obj)}")
 
run_lines = [
    "=" * 80,
    "Epithelial scArches Query Mapping Run Log — v1.2",
    "=" * 80,
    f"written_at:              {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    f"query_h5ad:              {QUERY_H5AD}",
    f"reference_scvi:          {SCVI_REF_DIR}",
    f"reference_scanvi_major:  {SCANVI_MAJOR_REF}",
    f"reference_scanvi_fine:   {SCANVI_FINE_REF}",
    f"output_dir:              {OUT_BASE}",
    f"saved_h5ad:              {h5ad_out}",
    f"saved_h5ad_size_gb:      {size_gb:.3f}",
    f"cells:                   {adata_main.n_obs:,}",
    f"full_genes_main:         {adata_main.n_vars:,}",
    f"hvg_genes_ref:           {len(hvg_genes_ref):,}",
    f"batches:                 {adata_main.obs[BATCH_KEY].nunique()}",
    f"gene_overlap_pct:        {overlap_pct:.1f}%",
    f"n_missing_genes_padded:  {len(missing_in_query)}",
    f"low_conf_threshold:      {LOW_CONF_THRESHOLD}",
    f"covid_filter_enabled:    {FILTER_COVID_RELATED}",
    f"covid_removed_cells:     {covid_removed_n:,}",
    f"covid_removed_fraction:  {covid_removed_fraction:.1%}",
    f"elapsed_minutes:         {elapsed_min:.2f}",
    f"device:                  {DEVICE}",
    "",
    "Predicted fine labels (filtered):",
    adata_main.obs['scanvi_fine_pred_filtered'].value_counts().to_string(),
    "",
    "Predicted major labels:",
    adata_main.obs['scanvi_major_pred'].value_counts().to_string(),
    "",
    f"Low-confidence cells (<{LOW_CONF_THRESHOLD}): "
    f"{low_conf_mask.sum():,} ({low_conf_mask.mean():.1%})",
    "",
    "UMAP note: query-only coordinate space; reference UMAP operator not saved in v2.7-HOTFIX.",
]
(OUT_BASE / "pipeline_run_log.txt").write_text('\n'.join(run_lines) + '\n', encoding='utf-8')
 
struct_lines = [
    "=" * 80,
    "AnnData Structure Summary — epithelial_query_mapped_v1_2",
    "=" * 80,
    f"shape: {adata_main.shape}",
    _mat_line('X', adata_main.X),
    "",
    f"obs ({len(adata_main.obs.columns)} columns):",
]
for col in adata_main.obs.columns:
    struct_lines.append(
        f"  {col}: dtype={adata_main.obs[col].dtype}  "
        f"unique={adata_main.obs[col].nunique(dropna=False)}"
    )
struct_lines += [
    "", f"layers ({len(adata_main.layers)}):",
    *[_mat_line(k, v) for k, v in adata_main.layers.items()],
    "", f"obsm ({len(adata_main.obsm)}):",
    *[f"  {k}: shape={getattr(v,'shape','NA')}" for k, v in adata_main.obsm.items()],
    "", f"uns scarches_mapping keys: {list(adata_main.uns.get('scarches_mapping', {}).keys())}",
    "",
    f"raw: {adata_main.raw is not None}",
]
if adata_main.raw is not None:
    struct_lines.append(f"  raw.shape: {adata_main.raw.shape}")
(OUT_BASE / "anndata_structure.txt").write_text('\n'.join(struct_lines) + '\n', encoding='utf-8')
 
print(f"[OK] pipeline_run_log.txt")
print(f"[OK] anndata_structure.txt")
 
elapsed_min = (time.time() - PIPELINE_START) / 60
print(f"\n{'='*80}")
print(f"PIPELINE COMPLETE  |  Elapsed: {elapsed_min:.1f} min")
print(f"{'='*80}")
print(f"Output:  {h5ad_out}")
print(f"Figures: {FIG_DIR}")
print("=" * 80)