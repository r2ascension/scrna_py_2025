#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stromal scArches Query Mapping v2.1
=====================================
Author  : r2end
Date    : 2026-03-15
Version : v2.1

Fixes over v2.0 (self-review):

  [P0-A] SCANVI label space error (critical correctness bug in v2.0):
    v2.0 used SCANVI.from_scvi_model(qry_scvi) with all cells = Unknown.
    This creates a SCANVI that has NO knowledge of reference classifier head.
    predict() would return an empty/random label space, NOT reference L3 labels.

    Fix: Use SCANVI.load_query_data(adata, ref_scanvi) which inherits the
    reference classifier head (correct L3 label space) while still performing
    encoder architecture surgery.

    Correct two-step workflow (Route B):
      Step 6 : SCVI.load_query_data  -> fine-tune -> X_scvi  (encoder surgery)
      Step 7a: SCANVI.prepare_query_anndata (labels setup on top of SCVI setup)
      Step 7b: SCANVI.load_query_data(adata, ref_scanvi) -> fine-tune -> predict
               ref_scanvi loaded with adata that has labels_key='scanvi_label'
               (reference v1.3 reintegration used exactly this key name)

  [P0-B] Reference SCANVI labels_key mismatch:
    v2.0 defined LABELS_KEY = 'cell_type_scarches_label'.
    Reference v1.3 reintegration used labels_key='scanvi_label'.
    SCANVI.load(..., adata) checks for the registered labels_key column in adata.

    Fix: Use REF_SCANVI_LABELS_KEY = 'scanvi_label' to match reference exactly.
    Create this column with all UNLABELED values before loading ref_scanvi.

  [P1]  Gene alignment - unreliable lil_matrix bulk column assignment:
    scipy lil_matrix[:, list_of_ints] = sparse_submatrix is not guaranteed
    to produce correct results across scipy versions.

    Fix: COO scatter approach - extract present columns as COO, remap column
    indices in-place, build final CSR via COO constructor. O(nnz), correct.

  [P2-A] ref_scvi released too late:
    v2.0 kept ref_scvi alive through step 6 fine-tuning despite being unneeded
    after load_query_data returns.
    Fix: del ref_scvi immediately after load_query_data(adata, ref_scvi).

  [P2-B] Unused variable absent_in_ref in v2.0.
    Fix: removed.

Reference (stromal_reintegration v1.3):
  scVI   : <REF_BASE_DIR>/models/scvi_stromal_v1_3
  scANVI : <REF_BASE_DIR>/models/scanvi_stromal_v1_3
  HVG    : scvi_stromal_v1_3/hvg_genes.csv  (no header, one gene per line)
  Reference SCANVI labels_key = 'scanvi_label'  (L3 cluster prefix labels)
  Reference SCANVI unlabeled_category = 'Unknown'

Query inputs (merged into single object before mapping):
  fibroblast_cells.h5ad   expected_L2 = Fibroblast
  endothelial_cells.h5ad  expected_L2 = Endothelial
  smc_cells.h5ad          expected_L2 = Smooth_Muscle (incl. Pericyte; L3 resolves)

Outputs:
  <OUTPUT_DIR>/
    adata_stromal_query_mapped_v2_1.h5ad  <- main: HVG .X + full-gene .raw
    models/
      query_scvi_stromal_v2_1/
      query_scanvi_stromal_v2_1/
    figures/
      query_overview.pdf
      query_per_subset_labels.pdf
      query_confidence_histogram.pdf
      query_expected_L2_audit.pdf
    expected_L2_audit.csv
    label_summary.csv
    mapping_report.txt
"""

# ============================================================
# 0. Environment  (before torch / scvi imports)
# ============================================================
import os
for _k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[_k] = "8"

import matplotlib; matplotlib.use('Agg')
import warnings; warnings.filterwarnings('ignore')
import gc, time, inspect
import numpy as np
import pandas as pd
from scipy import sparse
from pathlib import Path

import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt

PIPELINE_START = time.time()

np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)
scvi.settings.seed = 42
scvi.settings.dl_num_workers = 0
sc.settings.verbosity = 2
sc.settings.n_jobs = 16
USE_GPU = torch.cuda.is_available()
ACCEL  = 'gpu' if USE_GPU else 'cpu'
print(f"GPU: {USE_GPU}" + (f"  ({torch.cuda.get_device_name(0)})" if USE_GPU else ""))

# %%
# ============================================================
# CONFIGURATION
# ============================================================

# Reference model directory (stromal_reintegration v1.3 output)
REF_BASE_DIR   = Path("/home/h2048/data/py/0308/stromal_reintegration_v1_3")
SCVI_REF_DIR   = REF_BASE_DIR / "models" / "scvi_stromal_v1_3"
SCANVI_REF_DIR = REF_BASE_DIR / "models" / "scanvi_stromal_v1_3"
HVG_FILE       = SCVI_REF_DIR / "hvg_genes.csv"    # no header, one gene per line

# Query data (three stromal subsets; merged before mapping)
QUERY_DIR = Path("/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets")
QUERY_CONFIGS = [
    {"name": "fibroblast",  "h5ad": QUERY_DIR / "fibroblast_cells.h5ad",
     "expected_L2": "Fibroblast"},
    {"name": "endothelial", "h5ad": QUERY_DIR / "endothelial_cells.h5ad",
     "expected_L2": "Endothelial"},
    {"name": "smc",         "h5ad": QUERY_DIR / "smc_cells.h5ad",
     "expected_L2": "Smooth_Muscle"},   # includes Pericyte; L3 resolves
]

OUTPUT_DIR = Path("/home/h2048/data/py/0315/stromal_scarches_query_v2_1")
for _d in [OUTPUT_DIR, OUTPUT_DIR / "models", OUTPUT_DIR / "figures"]:
    _d.mkdir(parents=True, exist_ok=True)

# Keys (must match reference training)
BATCH_KEY              = "sample"
UNLABELED              = "Unknown"
TISSUE_KEY             = "tissue"     # metadata only; not used in model training

# CRITICAL: must match the labels_key used in reference v1.3 reintegration SCANVI.
# Reference code: scvi.model.SCANVI.from_scvi_model(..., labels_key='scanvi_label')
REF_SCANVI_LABELS_KEY  = "scanvi_label"

# Fine-tuning parameters
MAX_EPOCHS_SCVI_QRY   = 200    # scVI query surgery fine-tune
MAX_EPOCHS_SCANVI_QRY = 100    # scANVI query surgery fine-tune
BATCH_SIZE            = 256
WEIGHT_DECAY          = 0.0    # CRITICAL: must be 0 for scArches surgery
LR_QRY                = 5e-4

CONFIDENCE_THRESHOLD  = 0.5    # label = UNLABELED if max soft prob < threshold
MIN_GENE_OVERLAP      = 0.70   # fail-fast if gene overlap below this fraction

# L3 -> L2 mapping for QC audit (from reference ANNOTATION_TABLE, non-contaminant)
L3_TO_L2 = {
    "Endothelia_Lymphatic":                   "Endothelial",
    "Endothelia_vascular_Cap_a":              "Endothelial",
    "Endothelia_vascular_Cap_g":              "Endothelial",
    "Endothelia_vascular_arterial_pulmonary": "Endothelial",
    "Endothelia_vascular_arterial_systemic":  "Endothelial",
    "Endothelia_vascular_venous_pulmonary":   "Endothelial",
    "Endothelia_vascular_venous_systemic":    "Endothelial",
    "Fibro_adventitial":                      "Fibroblast",
    "Fibro_alveolar":                         "Fibroblast",
    "Fibro_myofibroblast":                    "Fibroblast",
    "Fibro_peribronchial":                    "Fibroblast",
    "Fibro_stress_activated":                 "Fibroblast",
    "Muscle_pericyte_pulmonary":              "Pericyte",
    "Muscle_pericyte_systemic":               "Pericyte",
    "Muscle_perivascular_immune_recruiting":  "Smooth_Muscle",
    "Muscle_smooth_arterial_systemic":        "Smooth_Muscle",
    "Muscle_smooth_pulmonary":                "Smooth_Muscle",
    "Schwann_nonmyelinating":                 "Schwann",
}

# ============================================================
# COVARIATE GENE SETS (for obs-level QC; NOT registered in model)
# ============================================================

STRESS_GENES = [
    "ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3","BANF1",
    "CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2","CNOT4","CNOT6",
    "DDIT4","DNAJA4","DNAJB9","DNAJC3","EDC4","EDEM1","EEF2","EIF2AK3","EIF2S1",
    "EIF4A1","EIF4A2","EIF4A3","EIF4E","EIF4EBP1","EIF4G1","ERN1","ERO1A",
    "FKBP14","GEMIN4","HERPUD1","HSP90B1","HSPA5","HSPA9","HYOU1","IFIT1",
    "KDELR3","MTHFD2","NFYA","NFYB","NPM1","PDIA5","PDIA6","PSAT1",
    "SLC1A4","SLC7A5","STC2","VEGFA","WFS1","XBP1","YWHAZ","ZBTB17",
]

S_GENES = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','HELLS','RFC2','RPA2','NASP','RAD51AP1',
    'GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51',
    'RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1',
    'CLSPN','POLA1','CHAF1B','BRIP1','E2F8',
]

G2M_GENES = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
    'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','SMC4','CCNB2','CKAP2L',
    'CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP',
    'CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5',
    'CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN','LBR',
    'CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA',
]


# %%
# ============================================================
# STEP 0: Validate reference models + load HVG gene list
# ============================================================
print("\n" + "="*70)
print("STEP 0: VALIDATE REFERENCE MODELS + LOAD HVG LIST")
print("="*70)

for _p, _lbl in [(SCVI_REF_DIR,   "scVI ref dir"),
                 (SCANVI_REF_DIR,  "scANVI ref dir"),
                 (HVG_FILE,        "HVG file")]:
    if not Path(_p).exists():
        raise FileNotFoundError(f"[ERROR] Missing {_lbl}: {_p}")
    print(f"[OK] {_lbl}: {_p}")

# Load reference HVG gene list (no header, one gene per line)
ref_hvg_genes = pd.read_csv(HVG_FILE, header=None)[0].tolist()
ref_hvg_set   = set(ref_hvg_genes)
assert len(ref_hvg_genes) == len(ref_hvg_set), (
    f"[ERROR] Duplicate genes in HVG file: "
    f"{len(ref_hvg_genes)} rows vs {len(ref_hvg_set)} unique"
)
print(f"\nReference HVG : {len(ref_hvg_genes)} genes (verified unique)")
print(f"Sample        : {ref_hvg_genes[:5]}")

# Verify SCANVI var_names consistency if present
_vn_path = SCANVI_REF_DIR / "SCANVI_var_names.csv"
if _vn_path.exists():
    _vn = pd.read_csv(_vn_path, header=None)[0].tolist()
    assert _vn == ref_hvg_genes, (
        f"[ERROR] SCANVI_var_names.csv ({len(_vn)}) != hvg_genes.csv "
        f"({len(ref_hvg_genes)}) or gene order differs. "
        "Inconsistent reference; stop and verify."
    )
    print("[OK] SCANVI var_names consistent with HVG list (exact match)")


# %%
# ============================================================
# STEP 1: Load all query subsets + merge into ONE object
# ============================================================
print("\n" + "="*70)
print("STEP 1: LOAD + MERGE QUERY SUBSETS")
print("="*70)

subset_adatas = []
for cfg in QUERY_CONFIGS:
    name = cfg["name"]
    path = cfg["h5ad"]
    print(f"\n  [{name}]  {path}")
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] Not found: {path}")

    ad = sc.read_h5ad(path)
    print(f"    shape : {ad.n_obs:,} x {ad.n_vars:,}")
    print(f"    obs   : {list(ad.obs.columns)}")

    # Guard: var_names uniqueness
    if not ad.var_names.is_unique:
        print("    [WARN] Duplicate var_names; making unique")
        ad.var_names_make_unique()
    assert ad.var_names.is_unique, f"[ERROR] var_names still not unique in {name}"

    # Resolve counts layer
    if 'counts' not in ad.layers:
        if ad.raw is not None:
            print("    [INFO] No counts layer; recovering from .raw.X")
            ad.layers['counts'] = sparse.csr_matrix(
                ad.raw[:, ad.var_names].X
            ).astype(np.float32)
        elif float(np.asarray(ad.X.max()).ravel()[0]) > 25:
            print(f"    [INFO] No counts layer; treating .X as raw counts")
            ad.layers['counts'] = sparse.csr_matrix(ad.X).astype(np.float32)
        else:
            raise ValueError(
                f"[ERROR] {name}: no counts layer and .X appears log-normalized "
                f"(max={float(np.asarray(ad.X.max()).ravel()[0]):.3f}). "
                "Raw counts required for scArches mapping."
            )

    # Ensure BATCH_KEY exists
    if BATCH_KEY not in ad.obs.columns:
        _cands = [c for c in ad.obs.columns
                  if any(k in c.lower() for k in ('sample', 'batch', 'donor'))]
        if _cands:
            print(f"    [INFO] '{BATCH_KEY}' absent; aliasing from '{_cands[0]}'")
            ad.obs[BATCH_KEY] = ad.obs[_cands[0]].astype(str)
        else:
            print(f"    [WARN] No batch column found; assigning 'query_{name}'")
            ad.obs[BATCH_KEY] = f"query_{name}"

    # Apply subset prefix to prevent cross-subset sample ID collision
    ad.obs[BATCH_KEY] = "qry_" + name + "_" + ad.obs[BATCH_KEY].astype(str)

    # Ensure TISSUE_KEY exists (metadata only)
    if TISSUE_KEY not in ad.obs.columns:
        ad.obs[TISSUE_KEY] = f"query_{name}"
        print(f"    [INFO] '{TISSUE_KEY}' absent; set to 'query_{name}'")

    ad.obs['query_subset'] = name
    subset_adatas.append(ad)
    print(f"    [OK] {ad.n_obs:,} cells | batches={ad.obs[BATCH_KEY].nunique()}")

# Concatenate with outer join (pads missing genes with 0)
print(f"\n  Concatenating {len(subset_adatas)} subsets (join='outer') ...")
adata_full = sc.concat(
    subset_adatas,
    join='outer',
    merge='same',
    label='query_subset',
    keys=[c['name'] for c in QUERY_CONFIGS],
)

# Sanitize NaN introduced by outer join
if sparse.issparse(adata_full.X):
    _d = adata_full.X.data
    _d[np.isnan(_d)] = 0.0
    adata_full.X.eliminate_zeros()

# Rebuild counts layer from merged .X (outer-join layer propagation is unreliable)
# At this point .X holds concatenated counts from the original subset counts layers
# (each subset had .X = log1p or counts depending on source; we set counts into layers
#  before concat, but .X was whatever the original file had). Safest: use the
# layer if it propagated correctly, otherwise rebuild from .X.
if ('counts' in adata_full.layers and
        adata_full.layers['counts'].shape == adata_full.shape):
    _layer_ok = True
    _lyr_max  = float(np.asarray(adata_full.layers['counts'].max()).ravel()[0])
    if _lyr_max < 1.0:
        print("  [WARN] counts layer looks normalized after concat; rebuilding from .X")
        _layer_ok = False
else:
    _layer_ok = False

if not _layer_ok:
    print("  [INFO] Rebuilding counts layer from merged .X")
    adata_full.layers['counts'] = sparse.csr_matrix(adata_full.X).astype(np.float32)
else:
    adata_full.layers['counts'] = sparse.csr_matrix(
        adata_full.layers['counts']
    ).astype(np.float32)

# Category dtype
for col in [BATCH_KEY, TISSUE_KEY, 'query_subset']:
    adata_full.obs[col] = adata_full.obs[col].astype(str).astype('category')

del subset_adatas; gc.collect()
print(f"\n  Merged: {adata_full.n_obs:,} cells x {adata_full.n_vars:,} genes")
print(f"  Batches: {adata_full.obs[BATCH_KEY].nunique()}")
print(adata_full.obs['query_subset'].value_counts().to_string())


# %%
# ============================================================
# STEP 2: Covariate calculation on FULL gene space
# CRITICAL ORDER: must run BEFORE gene alignment (Step 4) and
# prepare_query_anndata (Step 5a). Gene index operations below
# depend on adata_full.var_names being the full gene set.
# Note: ref v1.3 setup_anndata used only layer='counts' + batch_key.
# Covariates below are stored in .obs for downstream analysis only,
# NOT registered in the model.
# ============================================================
print("\n" + "="*70)
print("STEP 2: COVARIATE CALCULATION  (full-gene space)")
print("="*70)

counts_mat = adata_full.layers['counts']    # (n_cells x n_full_genes), sparse
gene_names  = adata_full.var_names

# --- 2a. MT percentage (safe: handles zero-count cells) ---
mt_mask = (gene_names.str.upper().str.startswith('MT-') |
           gene_names.str.upper().str.startswith('MT.'))
mt_genes = gene_names[mt_mask]
print(f"  MT genes found : {len(mt_genes)}")

if len(mt_genes) > 0:
    mt_idx   = np.where(mt_mask)[0]
    total    = np.asarray(counts_mat.sum(axis=1)).ravel().astype(float)
    mt_sum   = np.asarray(counts_mat[:, mt_idx].sum(axis=1)).ravel().astype(float)
    pct_mt   = np.zeros_like(total)
    nz       = total > 0
    pct_mt[nz] = mt_sum[nz] / total[nz] * 100.0
else:
    print("  [WARN] No MT genes detected; pct_counts_mt set to 0")
    pct_mt = np.zeros(adata_full.n_obs, dtype=float)
adata_full.obs['pct_counts_mt'] = pct_mt.astype(np.float32)
print(f"  pct_counts_mt  : mean={pct_mt.mean():.2f}%  max={pct_mt.max():.2f}%")

# --- 2b. Stress score (mean of stress signature, normalized to [0,1]) ---
stress_present = [g for g in STRESS_GENES if g in gene_names]
print(f"  Stress genes   : {len(stress_present)} / {len(STRESS_GENES)} present")

if len(stress_present) >= 5:
    s_idx       = np.where(gene_names.isin(stress_present))[0]
    stress_raw  = np.asarray(
        counts_mat[:, s_idx].mean(axis=1)
    ).ravel().astype(float)
    s_min, s_max = stress_raw.min(), stress_raw.max()
    stress_norm  = (stress_raw - s_min) / (s_max - s_min + 1e-12)
else:
    print("  [WARN] <5 stress genes; stress_score set to 0")
    stress_norm = np.zeros(adata_full.n_obs, dtype=float)
adata_full.obs['stress_score'] = stress_norm.astype(np.float32)
print(f"  stress_score   : mean={stress_norm.mean():.4f}")

# --- 2c. Cell cycle scores via scanpy (requires log-normalized .X) ---
# Create a lightweight temporary object; share obs only for batch key
print("  Computing cell cycle scores ...")
adata_tmp = sc.AnnData(
    X   = sparse.csr_matrix(counts_mat).astype(np.float32),
    obs = adata_full.obs[[BATCH_KEY]].copy(),
    var = adata_full.var.copy(),
)
sc.pp.normalize_total(adata_tmp, target_sum=1e4)
sc.pp.log1p(adata_tmp)

s_present   = [g for g in S_GENES   if g in gene_names]
g2m_present = [g for g in G2M_GENES if g in gene_names]
print(f"  S genes        : {len(s_present)} / {len(S_GENES)} present")
print(f"  G2M genes      : {len(g2m_present)} / {len(G2M_GENES)} present")

if len(s_present) >= 5 and len(g2m_present) >= 5:
    try:
        sc.tl.score_genes_cell_cycle(adata_tmp,
                                      s_genes=s_present,
                                      g2m_genes=g2m_present)
        adata_full.obs['S_score']   = adata_tmp.obs['S_score'].values.astype(np.float32)
        adata_full.obs['G2M_score'] = adata_tmp.obs['G2M_score'].values.astype(np.float32)
        adata_full.obs['phase']     = adata_tmp.obs['phase'].astype('category')
        print(f"  S_score        : mean={adata_full.obs['S_score'].mean():.4f}")
        print(f"  G2M_score      : mean={adata_full.obs['G2M_score'].mean():.4f}")
    except Exception as _e:
        print(f"  [WARN] Cell cycle scoring failed: {_e}; scores set to 0")
        adata_full.obs['S_score']   = np.float32(0.0)
        adata_full.obs['G2M_score'] = np.float32(0.0)
        adata_full.obs['phase']     = pd.Categorical(['G1'] * adata_full.n_obs)
else:
    print("  [WARN] Insufficient cell cycle genes; S_score/G2M_score set to 0")
    adata_full.obs['S_score']   = np.float32(0.0)
    adata_full.obs['G2M_score'] = np.float32(0.0)
    adata_full.obs['phase']     = pd.Categorical(['G1'] * adata_full.n_obs)

del adata_tmp; gc.collect()
print("  [OK] All covariates stored in .obs")


# Step 3 removed: adata_full.raw is not used as an intermediate.
# Full-gene .raw is constructed directly on adata_hvg in Step 4,
# using adata_full.layers['counts'] (shared memory, no copy).
# Background: AnnData .raw setter only accepts AnnData objects;
# assigning a Raw object from another AnnData raises ValueError.


# %%
# ============================================================
# STEP 4: Gene alignment to reference HVG space (COO scatter)
# COO approach: extract present columns as COO matrix, remap
# column indices to their reference positions, build CSR in one shot.
# O(nnz) time; no per-column loop; no lil_matrix assignment ambiguity.
# ============================================================
print("\n" + "="*70)
print("STEP 4: GENE ALIGNMENT TO REFERENCE HVG SPACE  (COO scatter)")
print("="*70)

assert adata_full.var_names.is_unique, "[ERROR] var_names not unique after concat"

qry_genes   = set(adata_full.var_names)
overlap     = qry_genes & ref_hvg_set
missing     = ref_hvg_set - qry_genes
overlap_pct = len(overlap) / len(ref_hvg_genes) * 100

print(f"  Query genes    : {len(qry_genes):,}")
print(f"  Ref HVG genes  : {len(ref_hvg_genes):,}")
print(f"  Overlap        : {len(overlap):,} ({overlap_pct:.1f}%)")
print(f"  Missing in qry : {len(missing):,}  (padded with zeros)")

if overlap_pct < MIN_GENE_OVERLAP * 100:
    raise ValueError(
        f"[ERROR] Gene overlap {overlap_pct:.1f}% < minimum "
        f"{MIN_GENE_OVERLAP*100:.0f}%. "
        "Check gene naming (symbols vs ENSEMBL) or source data."
    )

# Build bidirectional index maps
gene_to_qry_col = {g: i for i, g in enumerate(adata_full.var_names)}

present_ref_idx = []   # column positions in ref HVG space
present_qry_idx = []   # corresponding column positions in query space
for j, g in enumerate(ref_hvg_genes):
    if g in gene_to_qry_col:
        present_ref_idx.append(j)
        present_qry_idx.append(gene_to_qry_col[g])

present_ref_idx = np.asarray(present_ref_idx, dtype=np.int32)
present_qry_idx = np.asarray(present_qry_idx, dtype=np.int32)

# COO scatter: slice present query columns, remap to ref positions
src_mat = adata_full.layers['counts']    # (n_cells, n_full_genes) CSR/CSC
src_sub_coo = sparse.coo_matrix(
    src_mat[:, present_qry_idx]          # (n_cells, n_present)
)
# src_sub_coo.col is in [0, n_present); remap to ref column positions
new_col_indices = present_ref_idx[src_sub_coo.col]

n_cells = adata_full.n_obs
n_ref   = len(ref_hvg_genes)
aligned_counts = sparse.csr_matrix(
    (src_sub_coo.data.astype(np.float32),
     (src_sub_coo.row, new_col_indices)),
    shape=(n_cells, n_ref),
)

var_df = pd.DataFrame(index=ref_hvg_genes)
var_df['in_query'] = var_df.index.isin(qry_genes)

# Build HVG-space AnnData; .X = counts (required by prepare_query_anndata)
adata_hvg = sc.AnnData(
    X   = aligned_counts,
    obs = adata_full.obs.copy(),
    var = var_df,
)
adata_hvg.layers['counts'] = aligned_counts.copy()

# Preserve full-gene .raw on adata_hvg directly.
# IMPORTANT: .raw setter only accepts AnnData objects; a Raw object from another
# AnnData cannot be assigned (raises ValueError). Must construct a new AnnData.
# adata_full.layers['counts'] is shared by reference (no .copy()) -> zero extra memory.
adata_hvg.raw = sc.AnnData(
    X   = adata_full.layers['counts'],   # shared memory; NO .copy()
    obs = adata_full.obs.copy(),
    var = adata_full.var.copy(),
)

del src_sub_coo, aligned_counts; gc.collect()
print(f"  Aligned shape  : {adata_hvg.shape}")
print(f"  .raw preserved : {adata_hvg.raw.n_vars:,} full genes (shared memory)")


# %%
# ============================================================
# STEP 5a: SCVI prepare_query_anndata
# Reorders HVG space to match reference gene order exactly.
# Injects _scvi_batch and other model-manager columns into obs.
# MUST run before obs cleanup / labels column creation to avoid
# overwriting scArches-injected fields.
# ============================================================
print("\n" + "="*70)
print("STEP 5a: SCVI prepare_query_anndata")
print("="*70)

# Pass the model directory path directly -- no explicit SCVI.load() needed.
# Both prepare_query_anndata() and load_query_data() accept a path string.
# This avoids the transfer_field() category-validation path that fires when
# SCVI.load(..., adata=query_adata) is called and query batch categories
# are absent from the reference registry.
scvi.model.SCVI.prepare_query_anndata(adata_hvg, str(SCVI_REF_DIR))
print(f"  [OK] SCVI prepare_query_anndata | shape: {adata_hvg.shape}")

# Restore counts as .X (prepare may have modified .X)
if 'counts' in adata_hvg.layers:
    adata_hvg.X = adata_hvg.layers['counts']

# Verify batch key survived prepare
assert BATCH_KEY in adata_hvg.obs.columns, \
    f"[ERROR] '{BATCH_KEY}' lost during SCVI prepare_query_anndata"
adata_hvg.obs[BATCH_KEY] = adata_hvg.obs[BATCH_KEY].astype(str).astype('category')


# %%
# ============================================================
# STEP 5b: Add REF_SCANVI_LABELS_KEY column (all UNLABELED)
# CRITICAL: SCANVI.load() validates that the registered labels_key
# column exists in adata. Reference v1.3 used labels_key='scanvi_label'.
# All query cells are Unknown (pure label-transfer mode).
# This must happen AFTER SCVI prepare_query_anndata (step 5a) to avoid
# the column being overwritten, and BEFORE SCANVI.load (step 7a).
# ============================================================
print("\n" + "="*70)
print("STEP 5b: ADD REFERENCE SCANVI LABELS_KEY COLUMN")
print("="*70)

adata_hvg.obs[REF_SCANVI_LABELS_KEY] = UNLABELED
adata_hvg.obs[REF_SCANVI_LABELS_KEY] = adata_hvg.obs[REF_SCANVI_LABELS_KEY].astype('category')
print(f"  Column '{REF_SCANVI_LABELS_KEY}': {adata_hvg.n_obs:,} cells = '{UNLABELED}'")
print(f"  (This matches the labels_key used in reference v1.3 SCANVI training)")


# %%
# ============================================================
# STEP 6: scVI query fine-tune (architecture surgery on encoder)
# Purpose: obtain a query-adapted scVI latent (X_scvi) for visualization.
# The fine-tuned encoder weights feed into the SCANVI surgery below via
# the shared latent architecture. weight_decay=0 is mandatory.
# ============================================================
print("\n" + "="*70)
print("STEP 6: scVI QUERY FINE-TUNE  (weight_decay=0)")
print("="*70)

SCVI_QRY_DIR = str(OUTPUT_DIR / "models" / "query_scvi_stromal_v2_1")

# Pass path string directly; ref model object is not needed
qry_scvi = scvi.model.SCVI.load_query_data(adata_hvg, str(SCVI_REF_DIR))
gc.collect()
if torch.cuda.is_available():
    torch.cuda.empty_cache()

print(f"  n_latent     = {qry_scvi.module.n_latent}")
print(f"  cells        = {adata_hvg.n_obs:,}")
print(f"  HVGs         = {adata_hvg.n_vars:,}")
print(f"  max_epochs   = {MAX_EPOCHS_SCVI_QRY}")
print(f"  weight_decay = {WEIGHT_DECAY}  (CRITICAL for scArches)")

t0 = time.time()
qry_scvi.train(
    max_epochs             = MAX_EPOCHS_SCVI_QRY,
    batch_size             = BATCH_SIZE,
    early_stopping         = True,
    early_stopping_patience= 20,
    train_size             = 0.9,
    accelerator            = ACCEL,
    devices                = 1,
    plan_kwargs            = {'weight_decay': WEIGHT_DECAY, 'lr': LR_QRY},
)
print(f"  [OK] scVI fine-tune: {(time.time()-t0)/60:.1f} min")
qry_scvi.save(SCVI_QRY_DIR, overwrite=True)

adata_hvg.obsm['X_scvi'] = qry_scvi.get_latent_representation()
print(f"  [OK] X_scvi: {adata_hvg.obsm['X_scvi'].shape}")

del qry_scvi; gc.collect()
if torch.cuda.is_available():
    torch.cuda.empty_cache()


# %%
# ============================================================
# STEP 7a: SCANVI prepare_query_anndata (labels setup)
# Extends SCVI's gene/batch setup to also register the labels field.
# Must run before SCANVI.load_query_data.
# ============================================================
print("\n" + "="*70)
print("STEP 7a: SCANVI prepare_query_anndata  (labels field setup)")
print("="*70)

# Pass the model directory path directly -- no explicit SCANVI.load() needed.
# Same rationale as Step 5a: avoids transfer_field() category validation.
# Label class info (cell_type_mapping) will be read from the query model
# after load_query_data() in Step 7b.
scvi.model.SCANVI.prepare_query_anndata(adata_hvg, str(SCANVI_REF_DIR))
print(f"  [OK] SCANVI prepare_query_anndata done")


# %%
# ============================================================
# STEP 7b: scANVI query fine-tune (architecture surgery preserving classifier)
# SCANVI.load_query_data inherits reference classifier head (correct L3 label space).
# All query cells are Unknown; scANVI runs in label-transfer mode using
# reference-learned class structure to project query into label space.
# This is the fix for v2.0 from_scvi_model which had no reference label classes.
# ============================================================
print("\n" + "="*70)
print("STEP 7b: scANVI QUERY FINE-TUNE  (load_query_data, weight_decay=0)")
print("="*70)

SCANVI_QRY_DIR = str(OUTPUT_DIR / "models" / "query_scanvi_stromal_v2_1")

# load_query_data performs encoder surgery from ref_scanvi while keeping
# the reference classifier head -> predict() returns reference L3 labels
# Pass path string directly
qry_scanvi = scvi.model.SCANVI.load_query_data(adata_hvg, str(SCANVI_REF_DIR))
print(f"  [OK] SCANVI surgery from reference (classifier head inherited)")
if hasattr(qry_scanvi, 'cell_type_mapping'):
    _classes = qry_scanvi.cell_type_mapping
    print(f"  Label classes (ref L3): {len(_classes)} classes")
    for _lbl in sorted(_classes)[:8]:
        print(f"    {_lbl}")
    if len(_classes) > 8:
        print(f"    ... ({len(_classes)} total)")

gc.collect()
if torch.cuda.is_available():
    torch.cuda.empty_cache()

_train_kw = dict(
    max_epochs             = MAX_EPOCHS_SCANVI_QRY,
    batch_size             = BATCH_SIZE,
    early_stopping         = True,
    early_stopping_patience= 15,
    train_size             = 0.9,
    accelerator            = ACCEL,
    devices                = 1,
    plan_kwargs            = {'weight_decay': WEIGHT_DECAY, 'lr': LR_QRY},
)
if 'n_samples_per_label' in inspect.signature(qry_scanvi.train).parameters:
    _train_kw['n_samples_per_label'] = 100

t0 = time.time()
qry_scanvi.train(**_train_kw)
print(f"  [OK] scANVI fine-tune: {(time.time()-t0)/60:.1f} min")
qry_scanvi.save(SCANVI_QRY_DIR, overwrite=True)

# Extract latent representation + predictions
adata_hvg.obsm['X_scanvi']           = qry_scanvi.get_latent_representation()
adata_hvg.obs['cell_type_scarches_pred'] = qry_scanvi.predict()

# Soft predictions -> confidence and margin
soft_df = qry_scanvi.predict(soft=True)
soft_arr = soft_df.to_numpy() if hasattr(soft_df, 'to_numpy') else np.asarray(soft_df)

confidence   = soft_arr.max(axis=1).astype(np.float32)
sorted_soft  = np.sort(soft_arr, axis=1)[:, ::-1]
margin       = (sorted_soft[:, 0] - sorted_soft[:, 1]).astype(np.float32)

adata_hvg.obs['scarches_confidence']     = confidence
adata_hvg.obs['scarches_margin']         = margin

# Final label: Unknown if below confidence threshold
adata_hvg.obs['cell_type_scarches_final'] = np.where(
    confidence >= CONFIDENCE_THRESHOLD,
    adata_hvg.obs['cell_type_scarches_pred'].astype(str),
    UNLABELED
)

# Predicted L2 (for QC audit)
adata_hvg.obs['predicted_L2'] = (
    adata_hvg.obs['cell_type_scarches_pred']
    .astype(str)
    .map(L3_TO_L2)
    .fillna('Unknown')
)

# Category dtype
for _col in ['cell_type_scarches_pred', 'cell_type_scarches_final', 'predicted_L2']:
    adata_hvg.obs[_col] = adata_hvg.obs[_col].astype('category')

if hasattr(soft_df, 'columns'):
    adata_hvg.uns['scarches_label_classes'] = soft_df.columns.tolist()

low_conf_pct = (confidence < CONFIDENCE_THRESHOLD).sum() / adata_hvg.n_obs * 100
print(f"\n  X_scanvi shape      : {adata_hvg.obsm['X_scanvi'].shape}")
print(f"  Confidence median   : {np.median(confidence):.3f}")
print(f"  Margin median       : {np.median(margin):.3f}")
print(f"  Low confidence (<{CONFIDENCE_THRESHOLD}): {low_conf_pct:.1f}%")
if low_conf_pct > 30:
    print("  [WARN] >30% low-confidence predictions. "
          "Check reference coverage and gene overlap.")

del qry_scanvi; gc.collect()


# %%
# ============================================================
# STEP 8: UMAP on unified X_scanvi (single comparable latent space)
# All three query subsets were mapped by ONE model; X_scanvi is
# in a single coordinate frame -> combined UMAP is geometrically valid.
# ============================================================
print("\n" + "="*70)
print("STEP 8: UMAP  (unified X_scanvi)")
print("="*70)

sc.pp.neighbors(adata_hvg, use_rep='X_scanvi', n_neighbors=30,
                random_state=42, key_added='neighbors_scanvi')
sc.tl.umap(adata_hvg, neighbors_key='neighbors_scanvi',
           random_state=42, min_dist=0.5)
print(f"  [OK] UMAP on {adata_hvg.n_obs:,} cells (unified latent)")


# %%
# ============================================================
# STEP 9: Figures
# ============================================================
print("\n" + "="*70)
print("STEP 9: FIGURES")
print("="*70)

FIG_DIR = OUTPUT_DIR / "figures"
sc.settings.vector_friendly = True   # rasterizes scatter points

# --- 9a. 2x2 overview: subset | final label | confidence | margin ---
fig, axes = plt.subplots(2, 2, figsize=(22, 18))
sc.pl.embedding(adata_hvg, basis='umap', color='query_subset',
                ax=axes[0, 0], show=False, frameon=False, size=2,
                legend_loc='right margin', legend_fontsize=8,
                title='Query Subset')
sc.pl.embedding(adata_hvg, basis='umap', color='cell_type_scarches_final',
                ax=axes[0, 1], show=False, frameon=False, size=2,
                legend_loc='right margin', legend_fontsize=6,
                title=f'scArches Final Label (L3, conf>={CONFIDENCE_THRESHOLD})')
sc.pl.embedding(adata_hvg, basis='umap', color='scarches_confidence',
                ax=axes[1, 0], show=False, frameon=False, size=2,
                cmap='RdYlGn', vmin=0, vmax=1,
                title='Prediction Confidence (max soft prob)')
sc.pl.embedding(adata_hvg, basis='umap', color='scarches_margin',
                ax=axes[1, 1], show=False, frameon=False, size=2,
                cmap='RdYlGn', vmin=0, vmax=0.8,
                title='Prediction Margin (top1 - top2)')
plt.suptitle('Stromal scArches Query Mapping v2.1  (unified latent)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
fig.savefig(FIG_DIR / 'query_overview.pdf', dpi=300, bbox_inches='tight')
plt.close('all')
print("  [OK] query_overview.pdf")

# --- 9b. Per-subset label panels ---
_n = len(QUERY_CONFIGS)
fig, axes = plt.subplots(1, _n, figsize=(11 * _n, 9))
if _n == 1: axes = [axes]
for ax, cfg in zip(axes, QUERY_CONFIGS):
    _sub = adata_hvg[adata_hvg.obs['query_subset'] == cfg['name']]
    sc.pl.embedding(_sub, basis='umap', color='cell_type_scarches_final',
                    ax=ax, show=False, frameon=False, size=3,
                    legend_loc='right margin', legend_fontsize=6,
                    title=f"{cfg['name']} (n={_sub.n_obs:,})")
plt.suptitle('scArches Final Labels per Subset', fontsize=13)
plt.tight_layout()
fig.savefig(FIG_DIR / 'query_per_subset_labels.pdf', dpi=300, bbox_inches='tight')
plt.close('all')
print("  [OK] query_per_subset_labels.pdf")

# --- 9c. Confidence histogram ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.hist(confidence, bins=50, color='steelblue', edgecolor='white', alpha=0.85)
ax.axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--',
           label=f'Threshold = {CONFIDENCE_THRESHOLD}')
ax.set_xlabel('Prediction Confidence (max soft probability)', fontsize=12)
ax.set_ylabel('Cell Count', fontsize=12)
ax.set_title('scArches Confidence Distribution', fontsize=13)
ax.legend()
fig.tight_layout()
fig.savefig(FIG_DIR / 'query_confidence_histogram.pdf', dpi=300, bbox_inches='tight')
plt.close('all')
print("  [OK] query_confidence_histogram.pdf")

gc.collect()


# %%
# ============================================================
# STEP 10: Expected L2 consistency audit
# Predicted L3 -> L2 (via L3_TO_L2) compared against expected_L2.
# Flags unexpected cross-lineage assignments (>5% of cells).
# ============================================================
print("\n" + "="*70)
print("STEP 10: EXPECTED L2 CONSISTENCY AUDIT")
print("="*70)

audit_rows = []
for cfg in QUERY_CONFIGS:
    sub_name    = cfg['name']
    expected_L2 = cfg['expected_L2']
    mask        = adata_hvg.obs['query_subset'] == sub_name
    sub_obs     = adata_hvg.obs[mask]
    n_total     = len(sub_obs)

    pred_L2   = sub_obs['predicted_L2']
    n_match   = (pred_L2 == expected_L2).sum()
    n_unknown = (sub_obs['cell_type_scarches_final'] == UNLABELED).sum()
    pct_match = n_match / n_total * 100

    print(f"\n  [{sub_name}]  expected_L2={expected_L2}")
    print(f"    n_cells             : {n_total:,}")
    print(f"    Predicted L2 match  : {n_match:,} ({pct_match:.1f}%)")
    print(f"    Low-conf Unknown    : {n_unknown:,} ({n_unknown/n_total*100:.1f}%)")
    print(f"    Predicted L2 distribution:")
    for lbl, n in pred_L2.value_counts().items():
        flag = ""
        if lbl == expected_L2:
            flag = "  <- EXPECTED"
        elif lbl != "Unknown" and n / n_total > 0.05:
            flag = "  [WARN] unexpected >5%"
        print(f"      {lbl}: {n:,} ({n/n_total*100:.1f}%){flag}")

    if pct_match < 50:
        print(f"    [WARN] <50% L2 match for {sub_name}. "
              "Check reference coverage or query data purity.")

    audit_rows.append({
        'subset':       sub_name,
        'expected_L2':  expected_L2,
        'n_cells':      n_total,
        'n_match_L2':   int(n_match),
        'pct_match_L2': round(pct_match, 2),
        'n_unknown':    int(n_unknown),
        'pct_unknown':  round(n_unknown / n_total * 100, 2),
    })

audit_df = pd.DataFrame(audit_rows)
audit_df.to_csv(OUTPUT_DIR / 'expected_L2_audit.csv', index=False)
print(f"\n  [OK] expected_L2_audit.csv")

# Audit barplot
fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(audit_df['subset'], audit_df['pct_match_L2'], color='steelblue', alpha=0.85)
ax.axhline(50, color='orange', linestyle='--', label='50% threshold')
ax.set_ylabel('% Cells Matching Expected L2', fontsize=12)
ax.set_title('Expected L2 Consistency per Query Subset', fontsize=13)
ax.set_ylim(0, 110)
for i, row in audit_df.iterrows():
    ax.text(i, row['pct_match_L2'] + 1.5, f"{row['pct_match_L2']:.0f}%",
            ha='center', fontsize=11)
ax.legend()
fig.tight_layout()
fig.savefig(FIG_DIR / 'query_expected_L2_audit.pdf', dpi=300, bbox_inches='tight')
plt.close('all')
print("  [OK] query_expected_L2_audit.pdf")


# %%
# ============================================================
# STEP 11: Save outputs
# ============================================================
print("\n" + "="*70)
print("STEP 11: SAVE OUTPUTS")
print("="*70)

adata_hvg.uns['stromal_scarches_mapping'] = {
    'version':                'v2.1',
    'date':                   '2026-03-15',
    'architecture':           'Route-B: SCVI.load_query_data + SCANVI.load_query_data',
    'reference_scvi':         str(SCVI_REF_DIR),
    'reference_scanvi':       str(SCANVI_REF_DIR),
    'ref_scanvi_labels_key':  REF_SCANVI_LABELS_KEY,
    'n_hvg':                  int(len(ref_hvg_genes)),
    'gene_overlap_pct':       round(overlap_pct, 2),
    'n_missing_genes':        len(missing),
    'n_cells_query':          int(adata_hvg.n_obs),
    'n_batches':              int(adata_hvg.obs[BATCH_KEY].nunique()),
    'max_epochs_scvi_qry':    MAX_EPOCHS_SCVI_QRY,
    'max_epochs_scanvi_qry':  MAX_EPOCHS_SCANVI_QRY,
    'weight_decay':           WEIGHT_DECAY,
    'confidence_threshold':   CONFIDENCE_THRESHOLD,
    'subsets':                [c['name'] for c in QUERY_CONFIGS],
    'fixes_over_v2_0': [
        'P0-A: SCANVI.load_query_data preserves ref classifier head (from_scvi_model had none)',
        'P0-B: REF_SCANVI_LABELS_KEY=scanvi_label matches reference v1.3',
        'P1: COO scatter gene alignment (no lil_matrix column assignment)',
        'P2-A: ref_scvi released immediately after load_query_data',
        'P2-B: removed unused variable absent_in_ref',
    ],
    'note_covariates': (
        'MT%, stress, S/G2M computed on full-gene space before alignment. '
        'Stored in .obs only; NOT registered in model '
        '(ref v1.3 used layer=counts + batch_key only).'
    ),
}

# Label summary CSV
adata_hvg.obs[[
    BATCH_KEY, TISSUE_KEY, 'query_subset',
    'pct_counts_mt', 'stress_score', 'S_score', 'G2M_score', 'phase',
    'cell_type_scarches_pred', 'predicted_L2',
    'scarches_confidence', 'scarches_margin',
    'cell_type_scarches_final',
]].to_csv(OUTPUT_DIR / 'label_summary.csv')
print("  [OK] label_summary.csv")

OUT_H5AD = OUTPUT_DIR / "adata_stromal_query_mapped_v2_1.h5ad"
adata_hvg.write_h5ad(OUT_H5AD, compression='gzip', compression_opts=9)
print(f"  [OK] {OUT_H5AD}")

# Human-readable report
elapsed = (time.time() - PIPELINE_START) / 60
report = [
    "Stromal scArches Query Mapping v2.1",
    f"Date        : 2026-03-15",
    f"Runtime     : {elapsed:.1f} min",
    f"Ref scVI    : {SCVI_REF_DIR}",
    f"Ref scANVI  : {SCANVI_REF_DIR}",
    f"n_cells     : {adata_hvg.n_obs:,}",
    f"n_hvg       : {adata_hvg.n_vars:,}",
    f"n_raw_genes : {adata_hvg.raw.n_vars:,}",
    f"Gene overlap: {overlap_pct:.1f}%",
    "",
    "Per-subset summary:",
]
for cfg in QUERY_CONFIGS:
    nm   = cfg['name']
    sub  = adata_hvg.obs[adata_hvg.obs['query_subset'] == nm]
    arow = audit_df[audit_df['subset'] == nm].iloc[0]
    report += [
        f"  [{nm}]  n={len(sub):,}",
        f"    pct_mt mean       : {sub['pct_counts_mt'].mean():.2f}%",
        f"    confidence median : {sub['scarches_confidence'].median():.3f}",
        f"    low-conf Unknown  : {(sub['cell_type_scarches_final']==UNLABELED).sum():,}",
        f"    expected_L2 match : {arow['pct_match_L2']:.1f}%",
        f"    top L3 predictions:",
    ]
    for lbl, n in sub['cell_type_scarches_final'].value_counts().head(6).items():
        report.append(f"      {lbl}: {n:,}")
    report.append("")

report += [
    "Output .obs columns:",
    "  cell_type_scarches_pred  : raw L3-level prediction (all cells)",
    "  predicted_L2             : L3->L2 mapped (for QC)",
    "  scarches_confidence      : max soft probability [0,1]",
    "  scarches_margin          : top1 - top2 soft probability",
    f"  cell_type_scarches_final : prediction if conf>={CONFIDENCE_THRESHOLD}, else {UNLABELED}",
    "  pct_counts_mt            : MT% (full gene space)",
    "  stress_score             : stress signature score [0,1]",
    "  S_score / G2M_score      : cell cycle phase scores",
    "  phase                    : G1 / S / G2M",
    "",
    "AnnData structure:",
    f"  .X                : counts (HVG space, {adata_hvg.n_vars:,} genes)",
    f"  .layers['counts'] : same counts (HVG space)",
    f"  .raw.X            : full-gene counts ({adata_hvg.raw.n_vars:,} genes)",
    f"  .obsm['X_scvi']   : query-adapted scVI latent (encoder surgery)",
    f"  .obsm['X_scanvi'] : query-adapted scANVI latent (unified; use for downstream)",
    f"  .obsm['X_umap']   : UMAP from unified X_scanvi",
]

with open(OUTPUT_DIR / 'mapping_report.txt', 'w') as f:
    f.write('\n'.join(report))
print("  [OK] mapping_report.txt")

print(f"\n{'='*70}")
print(f"PIPELINE COMPLETE  ({elapsed:.1f} min)")
print(f"  Output : {OUT_H5AD}")
print(f"  Cells  : {adata_hvg.n_obs:,}  |  HVG: {adata_hvg.n_vars:,}  |  .raw: {adata_hvg.raw.n_vars:,}")
print(f"{'='*70}")