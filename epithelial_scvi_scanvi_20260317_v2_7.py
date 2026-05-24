# %% [markdown]
# # Epithelial Subcluster scVI-scANVI Pipeline v2.7-HOTFIX — COMBINED (SELF + REF)
# #
# # **Purpose:**
# # - Apply refined L3 annotations to epithelial subclusters (with label renames)
# # - Train scVI **once** on all post-remap cells
# # - Branch A (SELF): AT2 retained — downstream analysis
# # - Branch B (REF):  AT/Alveolar removed — scArches reference building
# #
# # **HOTFIX changes vs v2.7:**
# # - P0-1: .raw now stores full-gene log1p (not counts) — fixes all viz
# # - P0-2: _init_scanvi no longer silently falls back to fresh SCANVI; raises instead
# # - P0-3: GLOBAL_STRESS_THRESHOLD computed once after Step 4, reused in both branches
# #         REF branch re-computes kNN purity after AT removal
# # - P1-2: AT removal now uses explicit label set (not regex)
# # - P1-3: pre-training mapping validation raises on unmapped labels
# # - P1-4: removed stale Basal_EMT_ECC typo from MAJOR_MAP_SELF
# # - P1-5: REF outputs manifest.json listing all model dependencies
# # - P2-4: all CSV exports use index=False
# # - QRM v3.4-30: sanitize_sparse_counts + assert finite + assert >= 0
# # - QRM v3.4-31: filter_genes(min_cells=3) + filter_cells(min_genes=200) before .raw
# # - QRM v3.4-32: 4-tier HVG fallback with explicit flavor per tier
# #
# # **Author:** r2end  **Date:** 2026-03-17  **Version:** 2.7-HOTFIX

# %% [markdown]
# ## Configuration

# %%
import os
os.environ["OMP_NUM_THREADS"]        = "8"
os.environ["OPENBLAS_NUM_THREADS"]   = "8"
os.environ["MKL_NUM_THREADS"]        = "8"
os.environ["VECLIB_MAXIMUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"]    = "8"

import gc
import json
import warnings
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from pathlib import Path
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import traceback

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# =============================================================================
# PATHS
# =============================================================================
INPUT_H5AD  = "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad"
MARKER_DIR  = Path("/home/h2048/data/R/0129/epithelial_interpret_v2_7_FIXED")
BASE_DIR    = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/epithelial_v2_7_HOTFIX")
SELF_DIR    = BASE_DIR / "SELF"
REF_DIR     = BASE_DIR / "REF"
CKPT_DIR    = BASE_DIR / "checkpoints"

for d in [BASE_DIR, CKPT_DIR,
          SELF_DIR/"figures", SELF_DIR/"scvi_model",
          SELF_DIR/"scanvi_major_model", SELF_DIR/"scanvi_fine_model",
          REF_DIR/"figures",
          REF_DIR/"scanvi_major_model", REF_DIR/"scanvi_fine_model"]:
    d.mkdir(parents=True, exist_ok=True)

SCVI_MODEL_PATH = str(SELF_DIR / "scvi_model")
HVG_FILE_PATH   = str(SELF_DIR / "scvi_model" / "hvg_genes.txt")

sc.settings.figdir = SELF_DIR / "figures"

# =============================================================================
# KEYS
# =============================================================================
BATCH_KEY            = 'sample'
UNLABELED_CATEGORY   = 'Unknown'
MIN_CELLS_PER_BATCH  = 30

# =============================================================================
# ANNOTATION RENAMES
# =============================================================================
LABEL_REMAP = {
    'AT2_Canonical':                   'AT2',
    'AT2_Inflammatory_Repair':         'AT2',
    'Epithelial_Cycling':              'AT2_Cycling',
    'Goblet_Mucin':                    'SMG_Mucous',
    'Secretory_Club':                  'Goblet',
    'Secretory_Club_AT2_Transitional': 'Club',
}
EXCLUDE_LABELS = {'Mesenchymal_Contaminant'}

# P1-2 fix: explicit label set instead of regex (robust against future renames)
AT_LABELS_TO_REMOVE = {
    'AT1_Canonical',
    'AT1_MatrixRemodeling',
    'AT2',
    'AT2_Cycling',
}

# =============================================================================
# MAJOR LINEAGE MAPS
# P1-4 fix: removed stale 'Basal_EMT_ECC' typo
# =============================================================================
MAJOR_MAP_SELF = {
    # Alveolar (SELF only)
    'AT2':                       'Alveolar',
    'AT2_Cycling':               'Alveolar',
    'AT1_Canonical':             'Alveolar',
    'AT1_MatrixRemodeling':      'Alveolar',
    # Basal — includes suprabasal as transitional basal-lineage cells
    'Basal_Progenitor':          'Basal_Lineage',
    'Basal_Cycling':             'Basal_Lineage',
    'Basal_Inflammatory':        'Basal_Lineage',
    'Basal_EMT_ECM':             'Basal_Lineage',
    'Suprabasal_Progenitor':     'Basal_Lineage',
    'Suprabasal_Cycling':        'Basal_Lineage',
    # Ciliated
    'Ciliated_Mature':           'Ciliated_Lineage',
    'Ciliogenesis_Deuterosomal': 'Ciliated_Lineage',
    'Ciliated_Cycling_Immature': 'Ciliated_Lineage',
    # Secretory
    'Goblet':                    'Secretory_Lineage',
    'Club':                      'Secretory_Lineage',
    'SMG_Mucous':                'Secretory_Lineage',
    'Goblet_Defense_DUOX2':      'Secretory_Lineage',
    # SMG
    'SMG_Serous':                'SMG',
    'SMG_Duct_Secretory_Defense':'SMG',
    # Rare / specialized
    'Squamous_Metaplasia':       'Rare_Specialized',
    'Ionocyte_Brush':            'Rare_Specialized',
}

MAJOR_MAP_REF = {k: v for k, v in MAJOR_MAP_SELF.items() if v != 'Alveolar'}

# =============================================================================
# MARKER PANELS
# =============================================================================
L3_MARKER_PANELS = {
    # Alveolar
    "AT1_Canonical":             ["AGER","HOPX","CAV1","AQP4","RTKN2","CLDN18","EMP2"],
    "AT1_MatrixRemodeling":      ["AGER","CAV1","SPARC","COL4A1","COL4A2","SPOCK2"],
    "AT2":                       ["SFTPC","SFTPA1","SFTPA2","SFTPB","ABCA3","NAPSA",
                                   "SLC34A2","CHI3L1","CXCL8","SAA1"],
    "AT2_Cycling":               ["MKI67","TOP2A","UBE2C","BIRC5","AURKB","CENPA","CCNB1"],
    # Basal
    "Basal_Progenitor":          ["KRT5","KRT14","TP63","KRT15","KRT19","ITGA6","NGFR"],
    "Basal_Cycling":             ["KRT14","KRT5","TP63","MKI67","TOP2A","BIRC5"],
    "Basal_Inflammatory":        ["KRT17","CXCL8","CXCL1","CXCL2","TNFAIP3","FOS","JUN"],
    "Basal_EMT_ECM":             ["KRT14","TP63","NGFR","FN1","COL17A1","MMP2","VIM"],
    # Suprabasal — non-keratinizing transitional layer (KRT4/KRT13+, SPRR-)
    "Suprabasal_Progenitor":     ["KRT4","KRT13","KRT19","NOTCH1","NOTCH3","TP63","KRT5"],
    "Suprabasal_Cycling":        ["KRT4","KRT13","MKI67","TOP2A","BIRC5","TP63","KRT5"],
    # Ciliated
    "Ciliated_Mature":           ["FOXJ1","TPPP3","DNAH5","DNAH9","RSPH1","RFX2","RFX3"],
    "Ciliogenesis_Deuterosomal": ["DEUP1","CCNO","FOXN4","MCIDAS","CDC20B","E2F7","PLK4"],
    "Ciliated_Cycling_Immature": ["TPPP3","RSPH1","MKI67","TOP2A","FOXN4"],
    # Secretory (new names)
    "Goblet":                    ["SCGB1A1","SCGB3A1","SCGB3A2","AGR2","AGR3","CYP2F1"],
    "Club":                      ["SCGB1A1","SCGB3A1","SFTPB","NAPSA","GPR116","CLDN18"],
    "SMG_Mucous":                ["SPDEF","FOXA3","MUC5AC","MUC5B","FCGBP","TFF3",
                                   "BPIFB2","AZGP1"],
    "Goblet_Defense_DUOX2":      ["DUOX2","DUOXA2","LCN2","BPIFA2","CEACAM5"],
    # SMG
    "SMG_Serous":                ["LTF","LYZ","SLPI","DMBT1","BPIFA1","AZGP1","WFDC2"],
    "SMG_Duct_Secretory_Defense":["PIGR","SCGB3A1","TCN1","WFDC2","DMBT1","SLPI"],
    # Special
    "Squamous_Metaplasia":       ["SPRR1A","SPRR2A","SPRR2E","IVL","KRT6A","KLK7","S100A7"],
    "Ionocyte_Brush":            ["FOXI1","ASCL3","CFTR","ATP6V0D2","CLCNKA","CLCNKB","BSND"],
}

# =============================================================================
# VISUALIZATION
# =============================================================================
FIGURE_DPI    = 300
FIGURE_FORMAT = "pdf"
UMAP_SIZE     = 3
UMAP_ALPHA    = 0.6

L3_COLORS = {
    "AT1_Canonical":             "#1f77b4",
    "AT1_MatrixRemodeling":      "#aec7e8",
    "AT2":                       "#17becf",
    "AT2_Cycling":               "#c49c94",
    "Basal_Progenitor":          "#d62728",
    "Basal_Cycling":             "#ff7f0e",
    "Basal_Inflammatory":        "#ff9896",
    "Basal_EMT_ECM":             "#ffbb78",
    "Suprabasal_Progenitor":     "#e6550d",
    "Suprabasal_Cycling":        "#fdae6b",
    "Ciliated_Mature":           "#2ca02c",
    "Ciliogenesis_Deuterosomal": "#98df8a",
    "Ciliated_Cycling_Immature": "#8c564b",
    "Goblet":                    "#9467bd",
    "Club":                      "#c5b0d5",
    "SMG_Mucous":                "#e377c2",
    "Goblet_Defense_DUOX2":      "#f7b6d2",
    "SMG_Serous":                "#bcbd22",
    "SMG_Duct_Secretory_Defense":"#dbdb8d",
    "Squamous_Metaplasia":       "#7f7f7f",
    "Ionocyte_Brush":            "#c7c7c7",
}

# =============================================================================
# STRESS / CELL-CYCLE SIGNATURES
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
FORCE_INCLUDE_MARKERS = [
    'EPCAM','CDH1','ELF3','FXYD3','CLDN4','CLDN7',
    'TP63','KRT5','KRT14','KRT15','KRT17','ITGA6','ITGB4','COL17A1','KRT13','KRT4',
    'FOXJ1','RSPH1','PIFO','DNAH5','DNAH11','TUBA1A','TUBB4B','DEUP1','CCNO','CEP78',
    'MUC5AC','MUC5B','MUC4','TFF3','SPDEF','AGR2',
    'SCGB1A1','SCGB3A1','SCGB3A2','BPIFA1','WFDC2',
    'LYZ','LTF','DMBT1','MUC7',
    'KRT7','KRT19','KRT8','KRT18','AQP3',
    'CFTR','FOXI1','ATP6V1B1','DCLK1','TRPM5','POU2F3','GFI1B',
    'TP73','NOTCH1','NOTCH2','NOTCH3','DLL1','HES1','WNT5A','WNT7B','SOX2','SOX9',
    'SFTPC','SFTPB','SFTPA1','SFTPA2','ABCA3','SLC34A2','ETV5','LAMP3','LRRK2','MFSD2A',
    'MKI67','TOP2A','STMN1','CENPW','TK1',
]

# =============================================================================
# MODEL HYPERPARAMETERS
# =============================================================================
N_HVG                      = 4000
N_LATENT_SCVI              = 100
N_LATENT_SCANVI            = 100
N_LAYERS                   = 2
N_HIDDEN                   = 256
DROPOUT_RATE               = 0.2
GENE_LIKELIHOOD            = "nb"
DISPERSION                 = "gene-batch"
SCVI_ENCODE_COVARIATES     = True
SCVI_USE_LAYER_NORM        = "both"
SCVI_USE_BATCH_NORM        = "none"
SCVI_MAX_EPOCHS            = 400
SCANVI_MAX_EPOCHS          = 200
EARLY_STOPPING             = True
BATCH_SIZE                 = 1024
SCANVI_N_SAMPLES_PER_LABEL = 2000

PURITY_THRESHOLD  = 0.5
MT_THRESHOLD      = 20
STRESS_PERCENTILE = 95

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

PIPELINE_START = time.time()
print("=" * 80)
print("Epithelial Subcluster scVI-scANVI Pipeline v2.7-HOTFIX — COMBINED")
print("=" * 80)
print(f"Input:   {INPUT_H5AD}")
print(f"SELF:    {SELF_DIR}")
print(f"REF:     {REF_DIR}")
print(f"Device:  {DEVICE}")
print(f"Date:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# %% [markdown]
# ## Step 1: Create L3 Annotation Mapping

# %%
print("\n" + "=" * 80)
print("STEP 1: Creating L3 Annotation Mapping")
print("=" * 80)

cluster_to_l3_mapping = {
    # AT lineage
    "AT1 _0": "AT1_Canonical",
    "AT1 _1": "AT1_MatrixRemodeling",
    "AT1 _2": "AT2_Canonical",
    "AT2_0":  "AT2_Canonical",
    "AT2_1":  "AT2_Canonical",
    "AT2_2":  "AT2_Inflammatory_Repair",
    "AT2_3":  "Epithelial_Cycling",
    # Basal lineage
    "Basal_0": "Basal_EMT_ECM",
    "Basal_1": "Basal_Cycling",
    "Basal_2": "Basal_Inflammatory",
    "Basal_3": "Goblet_Mucin",
    "Basal_4": "Goblet_Mucin",
    "Basal_5": "Basal_Progenitor",
    # Ciliated lineage
    "Ciliated_0": "Ciliated_Mature",
    "Ciliated_1": "Ciliated_Mature",
    "Ciliated_2": "Ciliated_Cycling_Immature",
    "Ciliated_3": "Ciliated_Mature",
    "Ciliated_4": "Ciliated_Mature",
    "Ciliated_5": "Goblet_Mucin",
    "Deuterosomal_0": "Ciliogenesis_Deuterosomal",
    "Deuterosomal_1": "Ciliogenesis_Deuterosomal",
    # Dividing basal
    "Dividing_Basal_0": "Basal_Cycling",
    "Dividing_Basal_1": "Basal_Cycling",
    # Ionocyte
    "Ionocyte_n_Brush_0": "Ionocyte_Brush",
    "Ionocyte_n_Brush_1": "Ionocyte_Brush",
    # SMG lineage
    "SMG_Basal_0": "Basal_EMT_ECM",
    "SMG_Basal_1": "Mesenchymal_Contaminant",
    "SMG_Basal_2": "Basal_EMT_ECM",
    "SMG_Duct_0":  "SMG_Duct_Secretory_Defense",
    "SMG_Duct_1":  "Squamous_Metaplasia",
    "SMG_Duct_2":  "Squamous_Metaplasia",
    "SMG_Duct_3":  "Squamous_Metaplasia",
    "SMG_Duct_4":  "Squamous_Metaplasia",
    "SMG_Mucous_0": "Goblet_Mucin",
    "SMG_Mucous_1": "Ionocyte_Brush",
    "SMG_Serous_0": "SMG_Serous",
    "SMG_Serous_1": "SMG_Serous",
    "SMG_Serous_2": "SMG_Serous",
    # Secretory
    "Secretory_Club_0":    "Secretory_Club_AT2_Transitional",
    "Secretory_Goblet_0":  "Goblet_Defense_DUOX2",
    "Secretory_Goblet_1":  "Secretory_Club",
    "Secretory_Goblet_2":  "Basal_Progenitor",
    "Secretory_Goblet_3":  "Squamous_Metaplasia",
    "Secretory_Goblet_4":  "Ciliated_Cycling_Immature",
    # Suprabasal — v2.7: non-keratinizing vs pathological split
    "Suprabasal_0": "Suprabasal_Progenitor",   # KRT4+/KRT13+, non-keratinizing
    "Suprabasal_1": "Suprabasal_Cycling",       # KRT4+/MKI67+, proliferating transitional
    "Suprabasal_2": "Squamous_Metaplasia",      # SPRR+/IVL+, pathological
    "Suprabasal_3": "Squamous_Metaplasia",      # SPRR+/IVL+, pathological
}

mapping_df = pd.DataFrame(list(cluster_to_l3_mapping.items()),
                          columns=["Cluster", "cell_type_L3_original"])
mapping_df["cell_type_L3_refined"] = mapping_df["cell_type_L3_original"].map(
    lambda x: LABEL_REMAP.get(x, x) if x not in EXCLUDE_LABELS else "EXCLUDED"
)
mapping_df.to_csv(BASE_DIR / "cluster_to_L3_mapping.csv", index=False)

print(f"[OK] Total clusters: {len(cluster_to_l3_mapping)}")
print(f"     Unique original labels: {len(set(cluster_to_l3_mapping.values()))}")
print(f"     Unique final labels (post-remap): {mapping_df['cell_type_L3_refined'].nunique()}")

# %% [markdown]
# ## Step 2: Load Data and Apply L3 Annotations

# %%
print("\n" + "=" * 80)
print("STEP 2: Loading Data and Applying L3 Annotations")
print("=" * 80)

t0 = time.time()
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[OK] Loaded in {time.time() - t0:.1f}s  |  {adata.shape}")
print(f"  Layers: {list(adata.layers.keys())}  |  Has .raw: {adata.raw is not None}")

if 'subcluster' not in adata.obs.columns:
    raise ValueError("'subcluster' column not found in adata.obs!")
if 'counts' not in adata.layers:
    raise ValueError("'counts' layer is required but not found.")

# Map subcluster --> original L3 (strip hidden whitespace)
subcluster_clean = adata.obs['subcluster'].astype(str).str.strip()
n_ws = (subcluster_clean != adata.obs['subcluster'].astype(str)).sum()
if n_ws > 0:
    print(f"[WARN] {n_ws} subcluster labels had whitespace — stripped.")

mapping_clean = {str(k).strip(): v for k, v in cluster_to_l3_mapping.items()}
adata.obs['cell_type_L3'] = subcluster_clean.map(mapping_clean)

unmapped = adata.obs['cell_type_L3'].isna().sum()
if unmapped > 0:
    print(f"[WARN] {unmapped} cells with unmapped subclusters:")
    for cl in subcluster_clean[adata.obs['cell_type_L3'].isna()].unique():
        print(f"  '{cl}': {(subcluster_clean==cl).sum()} cells")
    adata.obs.loc[adata.obs['cell_type_L3'].isna(), 'cell_type_L3'] = \
        subcluster_clean[adata.obs['cell_type_L3'].isna()]

# Apply LABEL_REMAP
n_remapped = 0
for old, new in LABEL_REMAP.items():
    mask = adata.obs['cell_type_L3'] == old
    n = mask.sum()
    if n > 0:
        adata.obs.loc[mask, 'cell_type_L3'] = new
        n_remapped += n
        print(f"  Renamed: {old:<35} --> {new}  ({n:,} cells)")
print(f"  Total renamed: {n_remapped:,} cells")

# Hard-exclude Mesenchymal_Contaminant
excl_mask = adata.obs['cell_type_L3'].isin(EXCLUDE_LABELS)
n_excl = excl_mask.sum()
if n_excl > 0:
    print(f"[INFO] Removing {n_excl:,} Mesenchymal_Contaminant cells")
    adata = adata[~excl_mask].copy()

adata.obs['cell_type_L3'] = adata.obs['cell_type_L3'].astype('category')

# Remove small batches
bc = adata.obs[BATCH_KEY].value_counts()
sm = bc[bc < MIN_CELLS_PER_BATCH].index
if len(sm) > 0:
    print(f"[INFO] Removing {len(sm)} small batches (<{MIN_CELLS_PER_BATCH} cells)")
    adata = adata[~adata.obs[BATCH_KEY].isin(sm)].copy()

print(f"\n[OK] Final: {adata.n_obs:,} cells")
print("\nL3 cell type distribution:")
l3_counts = adata.obs['cell_type_L3'].value_counts()
for label, count in l3_counts.items():
    print(f"  {label:<35} {count:>8,}  ({count/adata.n_obs*100:.1f}%)")

# %% [markdown]
# ## Step 3: Load Marker Gene Data

# %%
print("\n" + "=" * 80)
print("STEP 3: Loading Marker Gene Data")
print("=" * 80)

markers_all      = MARKER_DIR / "all_markers.csv"
markers_filtered = MARKER_DIR / "top_markers_filtered.csv"

if not markers_all.exists():
    raise FileNotFoundError(f"Marker file not found: {markers_all}")

df_markers = pd.read_csv(markers_all)
print(f"[OK] Loaded: {markers_all.name}  ({len(df_markers):,} rows, {df_markers['cluster'].nunique()} clusters)")
if 'p_val_adj' in df_markers.columns:
    print(f"     Sig. markers (padj<0.05): {(df_markers['p_val_adj']<0.05).sum():,}")

if markers_filtered.exists():
    df_markers_filt = pd.read_csv(markers_filtered)
    print(f"[OK] Loaded filtered markers: {markers_filtered.name}  ({len(df_markers_filt):,} rows)")
else:
    print(f"[WARN] Filtered markers not found — using all_markers as fallback")
    df_markers_filt = df_markers.copy()

# %% [markdown]
# ## Step 4: Preprocessing & Covariates  *(shared)*

# %%
print("\n" + "=" * 80)
print("STEP 4: Preprocessing & Covariates")
print("=" * 80)

# --- QRM v3.4-30: sanitize_sparse_counts with assert ---
def _layer_flat_data(X):
    return X.data if sparse.issparse(X) else np.asarray(X).ravel()

def sanitize_sparse_counts(X, label="counts"):
    """Replace non-finite entries with 0 in-place; return (X_clean, n_bad)."""
    X = sparse.csr_matrix(X, dtype=np.float32, copy=True)
    bad = ~np.isfinite(X.data)
    n_bad = int(bad.sum())
    if n_bad > 0:
        X.data[bad] = 0.0
        X.eliminate_zeros()
        print(f"  [WARN] {label}: replaced {n_bad} non-finite entries with 0")
    return X, n_bad

if not sparse.issparse(adata.layers['counts']):
    adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])

# Global finite check (not sample check — QRM v3.4-28)
arr = _layer_flat_data(adata.layers['counts'])
print(f"[INFO] counts global check: finite={np.isfinite(arr).all()}, "
      f"inf={np.isinf(arr).sum()}, nan={np.isnan(arr).sum()}")

adata.layers['counts'], n_bad = sanitize_sparse_counts(adata.layers['counts'])

# Hard asserts after sanitize (QRM v3.4-30)
arr = _layer_flat_data(adata.layers['counts'])
assert np.isfinite(arr).all(), "counts still contains non-finite values after sanitize"
assert np.all(arr >= 0),       "counts contains negative values"
print("[OK] counts sanitize passed")

# Build log1p layer (needed for .raw and covariates)
if 'log1p' not in adata.layers:
    _t = sc.AnnData(X=adata.layers['counts'].copy(), var=adata.var.copy())
    sc.pp.normalize_total(_t, target_sum=1e4); sc.pp.log1p(_t)
    adata.layers['log1p'] = _t.X.copy()
    del _t; gc.collect()
adata.X = adata.layers['log1p'].copy()

# MT% — reuse existing column if valid
if 'pct_counts_mt' not in adata.obs.columns or (adata.obs['pct_counts_mt']==0).mean()>0.1:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], layer='counts',
                               percent_top=None, inplace=True)
print(f"[OK] pct_counts_mt mean: {adata.obs['pct_counts_mt'].mean():.2f}%")

# Stress score and cell-cycle — build ONE shared temp object to avoid double construction
print("[INFO] Computing stress score and cell-cycle on shared temp object...")
_t = sc.AnnData(X=adata.layers['counts'].copy(), var=adata.var.copy())
sc.pp.normalize_total(_t, target_sum=1e4); sc.pp.log1p(_t)

sg = [g for g in STRESS_SIGNATURE_GENES if g in _t.var_names]
if len(sg) > 10:
    idx_ = [list(_t.var_names).index(g) for g in sg]
    x_   = _t.X[:, idx_]
    adata.obs['stress_score'] = (np.asarray(x_.mean(axis=1)).flatten()
                                  if sparse.issparse(x_) else x_.mean(axis=1))
else:
    adata.obs['stress_score'] = 0.0
print(f"  stress_score mean: {adata.obs['stress_score'].mean():.4f}  (genes used: {len(sg)})")

s_p   = [g for g in S_GENES   if g in _t.var_names]
g2m_p = [g for g in G2M_GENES if g in _t.var_names]
if len(s_p) > 10 and len(g2m_p) > 10:
    sc.tl.score_genes_cell_cycle(_t, s_genes=s_p, g2m_genes=g2m_p)
    adata.obs['S_score']   = _t.obs['S_score'].values
    adata.obs['G2M_score'] = _t.obs['G2M_score'].values
    adata.obs['phase']     = _t.obs['phase'].values
    print(f"  Phase: {adata.obs['phase'].value_counts().to_dict()}")
else:
    adata.obs['S_score'] = adata.obs['G2M_score'] = 0.0
    adata.obs['phase'] = 'G1'
    print("  [WARN] Insufficient cell-cycle genes — set to 0/G1")

del _t; gc.collect()
print("[OK] Covariates ready")

# P0-3 fix: compute GLOBAL stress threshold here, once.
# Both branches will reuse this single value — never re-quantile on a subset.
GLOBAL_STRESS_THRESHOLD = float(adata.obs['stress_score'].quantile(STRESS_PERCENTILE / 100))
print(f"[OK] GLOBAL_STRESS_THRESHOLD (p{STRESS_PERCENTILE}): {GLOBAL_STRESS_THRESHOLD:.4f}")

# %% [markdown]
# ## Step 5: Gene Filtering + HVG Selection  *(shared)*

# %%
print("\n" + "=" * 80)
print("STEP 5: Gene Filtering + HVG Selection")
print("=" * 80)

# QRM v3.4-31: filter_genes + filter_cells BEFORE .raw
# (removes near-zero genes that destabilize seurat_v3 loess; .raw captures post-filter space)
n_cells_pre, n_genes_pre = adata.shape
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
print(f"[OK] filter_genes + filter_cells: {n_cells_pre:,}x{n_genes_pre:,} "
      f"-> {adata.n_obs:,}x{adata.n_vars:,}")

# QRM v3.4-32: 4-tier HVG fallback with explicit flavor per tier
# seurat_v3 <-> counts (raw UMI); seurat <-> log1p (normalized). Do NOT cross.
hvg_method = None
for tier, layer_, bkey, flav in [
    (1, 'counts', BATCH_KEY, 'seurat_v3'),
    (2, 'counts', None,      'seurat_v3'),
    (3, 'log1p',  BATCH_KEY, 'seurat'),
    (4, 'log1p',  None,      'seurat'),
]:
    try:
        kw = dict(layer=layer_, n_top_genes=N_HVG, flavor=flav, subset=False)
        if bkey:
            kw['batch_key'] = bkey
        sc.pp.highly_variable_genes(adata, **kw)
        hvg_method = f"tier{tier}: {flav} ({layer_}, batch={'yes' if bkey else 'no'})"
        break
    except Exception as e:
        print(f"  [WARN] HVG tier {tier} failed: {e}")

n_hvg = int(adata.var['highly_variable'].sum())
if n_hvg == 0:
    raise ValueError("No HVGs selected after all 4 fallback tiers — check input data")
adata.uns['hvg_method'] = hvg_method
print(f"[OK] {hvg_method}  |  HVGs: {n_hvg}")

# Force-include marker genes
n_added = 0
for g in FORCE_INCLUDE_MARKERS:
    if g in adata.var_names and not adata.var.loc[g, 'highly_variable']:
        adata.var.loc[g, 'highly_variable'] = True
        n_added += 1
print(f"[OK] Force-included {n_added} additional marker genes. Final HVG: {adata.var['highly_variable'].sum()}")

# P0-1 fix: .raw stores full-gene log1p (not counts).
# All downstream use_raw=True plots will read normalized expression.
adata.raw = sc.AnnData(
    X   = adata.layers['log1p'],   # normalized, NOT counts
    obs = adata.obs.copy(),
    var = adata.var.copy()
)
print(f"[OK] .raw set to full-gene log1p  |  shape: {adata.raw.shape}")

# Subset to HVG for model training
adata = adata[:, adata.var['highly_variable']].copy()
print(f"[OK] HVG subset for training: {adata.shape}")

os.makedirs(SCVI_MODEL_PATH, exist_ok=True)
pd.Series(adata.var_names.astype(str)).to_csv(HVG_FILE_PATH, index=False, header=False)

# Re-check small batches after HVG subset
bs = adata.obs[BATCH_KEY].value_counts()
sm = bs[bs < MIN_CELLS_PER_BATCH].index
if len(sm):
    adata = adata[~adata.obs[BATCH_KEY].isin(sm)].copy()
    adata.raw = adata.raw[adata.obs_names]
    print(f"[INFO] Removed {len(sm)} newly-small batches after HVG subset")

adata.write_h5ad(str(CKPT_DIR/"adata_preprocessed.h5ad"), compression='gzip')
print(f"[OK] Checkpoint saved.  Final training shape: {adata.shape}")

# %% [markdown]
# ## Step 6: scVI Training  *(shared — trained ONCE)*

# %%
print("\n" + "=" * 80)
print("STEP 6: scVI Training (shared)")
print("=" * 80)

scvi.model.SCVI.setup_anndata(
    adata, layer='counts', batch_key=BATCH_KEY,
    continuous_covariate_keys=['pct_counts_mt','stress_score','S_score','G2M_score']
)
scvi_model = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT_SCVI, n_layers=N_LAYERS, n_hidden=N_HIDDEN,
    dropout_rate=DROPOUT_RATE, gene_likelihood=GENE_LIKELIHOOD,
    dispersion=DISPERSION, encode_covariates=SCVI_ENCODE_COVARIATES,
    use_layer_norm=SCVI_USE_LAYER_NORM, use_batch_norm=SCVI_USE_BATCH_NORM
)
try:
    n_p = scvi_model.module.n_params
except AttributeError:
    n_p = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)
print(f"[INFO] scVI params: {n_p:,}")

scvi_model.train(max_epochs=SCVI_MAX_EPOCHS, batch_size=BATCH_SIZE,
                  early_stopping=EARLY_STOPPING, train_size=0.9,
                  plan_kwargs={'lr':1e-3})

scvi_model.save(SCVI_MODEL_PATH, overwrite=True)
adata.obsm['X_scvi'] = scvi_model.get_latent_representation()
print(f"[OK] scVI saved: {SCVI_MODEL_PATH}  |  Latent: {adata.obsm['X_scvi'].shape}")

# %% [markdown]
# ## Step 7: Unknown Cleaning  *(shared — using global thresholds)*

# %%
print("\n" + "=" * 80)
print("STEP 7: Unknown Cleaning (kNN purity + global quality gates)")
print("=" * 80)

# Build initial training labels for SELF (used to compute shared purity scores)
adata.obs['major_lineage']       = adata.obs['cell_type_L3'].map(MAJOR_MAP_SELF)
adata.obs.loc[adata.obs['major_lineage'].isna(), 'major_lineage'] = 'Other'
adata.obs['scanvi_labels_major'] = adata.obs['major_lineage'].astype('category')
adata.obs['scanvi_labels_fine']  = adata.obs['cell_type_L3'].copy()

for col in ['scanvi_labels_major', 'scanvi_labels_fine']:
    if UNLABELED_CATEGORY not in adata.obs[col].cat.categories:
        adata.obs[col] = adata.obs[col].cat.add_categories([UNLABELED_CATEGORY])

# kNN purity (k=30) in shared scVI latent space
print("[INFO] Computing kNN purity (k=30)...")
knn = NearestNeighbors(n_neighbors=30)
knn.fit(adata.obsm['X_scvi'])
indices = knn.kneighbors(return_distance=False)

lm = adata.obs['scanvi_labels_major'].values
lf = adata.obs['scanvi_labels_fine'].values
pm, pf = [], []
for i, nb in enumerate(indices):
    om  = lm[i]
    of_ = lf[i]
    pm.append(1.0 if om  == UNLABELED_CATEGORY else (lm[nb[1:]]==om ).mean())
    pf.append(1.0 if of_ == UNLABELED_CATEGORY else (lf[nb[1:]]==of_).mean())

adata.obs['label_purity_major'] = pm
adata.obs['label_purity_fine']  = pf

low_pm = adata.obs['label_purity_major'] < PURITY_THRESHOLD
low_pf = adata.obs['label_purity_fine']  < PURITY_THRESHOLD
print(f"  Major low purity (<{PURITY_THRESHOLD}): {low_pm.sum():,} ({low_pm.mean():.1%})")
print(f"  Fine  low purity (<{PURITY_THRESHOLD}): {low_pf.sum():,} ({low_pf.mean():.1%})")
adata.obs.loc[low_pm, 'scanvi_labels_major'] = UNLABELED_CATEGORY
adata.obs.loc[low_pf, 'scanvi_labels_fine']  = UNLABELED_CATEGORY

# Quality gates — use GLOBAL thresholds (not subset-derived)
high_mt = adata.obs['pct_counts_mt'] > MT_THRESHOLD
hs      = adata.obs['stress_score']  > GLOBAL_STRESS_THRESHOLD
print(f"  High mt% (>{MT_THRESHOLD}%): {high_mt.sum():,}  "
      f"High stress (>p{STRESS_PERCENTILE}={GLOBAL_STRESS_THRESHOLD:.4f}): {hs.sum():,}")
for m in [high_mt, hs]:
    adata.obs.loc[m, 'scanvi_labels_major'] = UNLABELED_CATEGORY
    adata.obs.loc[m, 'scanvi_labels_fine']  = UNLABELED_CATEGORY

unk_m = (adata.obs['scanvi_labels_major']==UNLABELED_CATEGORY).sum()
unk_f = (adata.obs['scanvi_labels_fine'] ==UNLABELED_CATEGORY).sum()
print(f"[OK] Shared cleaning done. Unknown major: {unk_m:,} ({unk_m/adata.n_obs:.1%})  "
      f"Unknown fine: {unk_f:,} ({unk_f/adata.n_obs:.1%})")

# %% [markdown]
# ## Step 8: Helper Functions

# %%
def _validate_l3_to_major_mapping(adata_b, major_map, branch_name):
    """P1-3 fix: raise on any unmapped L3 label before training begins."""
    present_labels = set(adata_b.obs['cell_type_L3'].cat.categories)
    missing = sorted(present_labels - set(major_map.keys()))
    if missing:
        raise ValueError(
            f"[{branch_name}] Unmapped L3 labels in MAJOR_MAP — "
            f"fix MAJOR_MAP before training: {missing}"
        )


def _build_labels_for_branch(adata_b, major_map, branch_name):
    """Rebuild scanvi label columns using global thresholds (not subset quantiles)."""
    # P1-3 validation first
    _validate_l3_to_major_mapping(adata_b, major_map, branch_name)

    adata_b.obs['scanvi_labels_major'] = (
        adata_b.obs['cell_type_L3'].map(major_map).astype('category')
    )
    adata_b.obs['scanvi_labels_fine'] = adata_b.obs['cell_type_L3'].copy()

    for col in ['scanvi_labels_major', 'scanvi_labels_fine']:
        if UNLABELED_CATEGORY not in adata_b.obs[col].cat.categories:
            adata_b.obs[col] = adata_b.obs[col].cat.add_categories([UNLABELED_CATEGORY])

    # Re-apply cleaning: purity from stored column + GLOBAL thresholds
    for mask in [
        adata_b.obs['label_purity_major'] < PURITY_THRESHOLD,
        adata_b.obs['label_purity_fine']  < PURITY_THRESHOLD,
        adata_b.obs['pct_counts_mt']       > MT_THRESHOLD,
        adata_b.obs['stress_score']        > GLOBAL_STRESS_THRESHOLD,  # P0-3 fix: global
    ]:
        adata_b.obs.loc[mask, 'scanvi_labels_major'] = UNLABELED_CATEGORY
        adata_b.obs.loc[mask, 'scanvi_labels_fine']  = UNLABELED_CATEGORY

    return adata_b


def _recompute_purity_for_branch(adata_b, branch_name):
    """P0-3 fix: Re-run kNN purity in branch-specific latent space.
    Must be called AFTER at_mask removal so the neighborhood reflects
    the actual cell composition the scANVI will see.
    """
    print(f"  [{branch_name}] Re-computing kNN purity (k=30) in branch latent space...")
    knn_b = NearestNeighbors(n_neighbors=30)
    knn_b.fit(adata_b.obsm['X_scvi'])
    indices_b = knn_b.kneighbors(return_distance=False)

    lm_b = adata_b.obs['scanvi_labels_major'].values
    lf_b = adata_b.obs['scanvi_labels_fine'].values
    pm_b, pf_b = [], []
    for i, nb in enumerate(indices_b):
        om  = lm_b[i]
        of_ = lf_b[i]
        pm_b.append(1.0 if om  == UNLABELED_CATEGORY else (lm_b[nb[1:]]==om ).mean())
        pf_b.append(1.0 if of_ == UNLABELED_CATEGORY else (lf_b[nb[1:]]==of_).mean())
    adata_b.obs['label_purity_major'] = pm_b
    adata_b.obs['label_purity_fine']  = pf_b

    # Re-apply purity gate with updated scores
    low_pm_b = adata_b.obs['label_purity_major'] < PURITY_THRESHOLD
    low_pf_b = adata_b.obs['label_purity_fine']  < PURITY_THRESHOLD
    print(f"    Major low purity: {low_pm_b.sum():,}  Fine low purity: {low_pf_b.sum():,}")
    adata_b.obs.loc[low_pm_b, 'scanvi_labels_major'] = UNLABELED_CATEGORY
    adata_b.obs.loc[low_pf_b, 'scanvi_labels_fine']  = UNLABELED_CATEGORY
    return adata_b


def _align_batch_categories_for_transfer(base_model, adata_b):
    """Restore full batch category list from trained scVI registry before transfer."""
    try:
        batch_registry   = (base_model.adata_manager.registry
                            .get('field_registries', {}).get('batch', {}))
        batch_state      = batch_registry.get('state_registry', {})
        batch_categories = batch_state.get('categorical_mapping')
        batch_obs_key    = batch_state.get('original_key', BATCH_KEY)
    except Exception:
        return adata_b

    if batch_categories is None or batch_obs_key not in adata_b.obs:
        return adata_b

    batch_categories = list(batch_categories)
    batch_series     = adata_b.obs[batch_obs_key]
    if not isinstance(batch_series.dtype, pd.CategoricalDtype):
        batch_series = batch_series.astype('category')

    if list(batch_series.cat.categories) != batch_categories:
        adata_b.obs[batch_obs_key] = pd.Categorical(
            adata_b.obs[batch_obs_key], categories=batch_categories
        )
        print(f"  [INFO] Restored {len(batch_categories)} batch categories for transfer "
              f"({adata_b.obs[batch_obs_key].nunique()} observed in subset)")
    return adata_b


def _init_scanvi(base_model, adata_b, labels_col):
    """P0-2 fix: no silent fallback to fresh SCANVI. Fail hard on error."""
    adata_b = _align_batch_categories_for_transfer(base_model, adata_b)
    try:
        return scvi.model.SCANVI.from_scvi_model(
            base_model,
            unlabeled_category=UNLABELED_CATEGORY,
            labels_key=labels_col,
            adata=adata_b
        )
    except Exception as e:
        raise RuntimeError(
            f"SCANVI.from_scvi_model failed for labels_key='{labels_col}'. "
            f"This pipeline relies on shared scVI — do NOT fall back to fresh SCANVI. "
            f"Fix the registry/category mismatch instead.\n"
            f"Original error: {e}"
        ) from e


def train_scanvi_pair(base_model, adata_b, major_path, fine_path, tag):
    for col, path, rep_key, pred_key, conf_key in [
        ('scanvi_labels_major', major_path, 'X_scanvi_major', 'scanvi_major_pred', 'scanvi_major_conf'),
        ('scanvi_labels_fine',  fine_path,  'X_scanvi_fine',  'scanvi_fine_pred',  'scanvi_fine_conf'),
    ]:
        label = 'major' if 'major' in col else 'fine'
        print(f"\n[{tag}] scANVI_{label}  "
              f"({adata_b.obs[col].nunique()} types, {adata_b.n_obs:,} cells)")
        m = _init_scanvi(base_model, adata_b, col)
        m.train(max_epochs=SCANVI_MAX_EPOCHS, batch_size=BATCH_SIZE,
                early_stopping=EARLY_STOPPING, train_size=0.9,
                n_samples_per_label=SCANVI_N_SAMPLES_PER_LABEL,
                plan_kwargs={'lr':1e-3})
        m.save(path, overwrite=True)
        adata_b.obsm[rep_key] = m.get_latent_representation()
        adata_b.obs[pred_key] = m.predict()
        p = m.predict(soft=True)
        adata_b.obs[conf_key] = (p.values.max(axis=1) if isinstance(p, pd.DataFrame)
                                  else np.asarray(p).max(axis=1))
        print(f"  Saved: {path}")
        print(f"  Predictions:\n{adata_b.obs[pred_key].value_counts().to_string()}")
        del m; gc.collect()
    return adata_b


def compute_umaps_leiden(adata_b, tag):
    sc.pp.neighbors(adata_b, use_rep='X_scanvi_major', n_neighbors=15, key_added='nbrs_major')
    sc.tl.umap(adata_b, neighbors_key='nbrs_major')
    adata_b.obsm['X_umap_major'] = adata_b.obsm['X_umap'].copy()

    sc.pp.neighbors(adata_b, use_rep='X_scanvi_fine', n_neighbors=15, key_added='nbrs_fine')
    sc.tl.umap(adata_b, neighbors_key='nbrs_fine')
    adata_b.obsm['X_umap_fine'] = adata_b.obsm['X_umap'].copy()

    sc.tl.leiden(adata_b, neighbors_key='nbrs_fine', resolution=0.5, key_added='leiden_fine')
    print(f"  [{tag}] Leiden clusters (res=0.5): {adata_b.obs['leiden_fine'].nunique()}")
    return adata_b


def get_available_markers(marker_dict, var_names):
    available, missing_all = {}, []
    for ct, markers in marker_dict.items():
        avail = [m for m in markers if m in var_names]
        miss  = [m for m in markers if m not in var_names]
        if avail: available[ct] = avail
        missing_all.extend(miss)
    return available, list(set(missing_all))


print("[OK] Helper functions defined")

# %% [markdown]
# ## Step 9: Branch A — SELF  *(AT2 retained)*

# %%
print("\n" + "=" * 80)
print("STEP 9: Branch A — SELF (AT2 retained)")
print("=" * 80)

adata_self = adata.copy()
# SELF uses shared purity (whole latent space including AT2) — no purity recompute needed
adata_self = _build_labels_for_branch(adata_self, MAJOR_MAP_SELF, 'SELF')

print("[SELF] Label distribution (major):")
print(adata_self.obs['scanvi_labels_major'].value_counts().to_string())
print("\n[SELF] Label distribution (fine):")
print(adata_self.obs['scanvi_labels_fine'].value_counts().to_string())

adata_self = train_scanvi_pair(
    scvi_model, adata_self,
    major_path=str(SELF_DIR / "scanvi_major_model"),
    fine_path=str(SELF_DIR  / "scanvi_fine_model"),
    tag='SELF'
)
adata_self = compute_umaps_leiden(adata_self, 'SELF')

l3_cats = list(adata_self.obs['cell_type_L3'].cat.categories)
adata_self.uns['cell_type_L3_colors'] = [
    L3_COLORS.get(c, plt.cm.tab20(i % 20)) for i, c in enumerate(l3_cats)
]
print("[OK] SELF scANVI complete")

# %% [markdown]
# ## Step 10: Branch B — REF  *(AT/Alveolar removed)*

# %%
print("\n" + "=" * 80)
print("STEP 10: Branch B — REF (AT/Alveolar removed)")
print("=" * 80)

# P1-2 fix: explicit label set (not regex) for AT removal
at_mask = adata.obs['cell_type_L3'].isin(AT_LABELS_TO_REMOVE)
n_at = at_mask.sum()
print(f"[INFO] AT_LABELS_TO_REMOVE: {sorted(AT_LABELS_TO_REMOVE)}")
print(f"[INFO] Removing {n_at:,} AT/Alveolar cells ({n_at/adata.n_obs:.2%})")
print(adata.obs.loc[at_mask, 'cell_type_L3'].value_counts().to_string())

adata_ref = adata[~at_mask].copy()

# Remove batches newly too small after AT removal
bs_r = adata_ref.obs[BATCH_KEY].value_counts()
sm_r = bs_r[bs_r < MIN_CELLS_PER_BATCH].index
if len(sm_r):
    print(f"[INFO] Removing {len(sm_r)} newly-small batches post AT-filter")
    adata_ref = adata_ref[~adata_ref.obs[BATCH_KEY].isin(sm_r)].copy()
print(f"[OK] REF cells: {adata_ref.n_obs:,}")

# Build REF labels first (needed for purity computation)
adata_ref = _build_labels_for_branch(adata_ref, MAJOR_MAP_REF, 'REF')

# P0-3 fix: re-compute kNN purity after AT removal.
# The AT-removed REF has a different neighborhood structure;
# using shared purity would penalize cells that were only "near AT" cells.
adata_ref = _recompute_purity_for_branch(adata_ref, 'REF')

print("\n[REF] Label distribution (major):")
print(adata_ref.obs['scanvi_labels_major'].value_counts().to_string())
print("\n[REF] Label distribution (fine):")
print(adata_ref.obs['scanvi_labels_fine'].value_counts().to_string())

adata_ref = train_scanvi_pair(
    scvi_model, adata_ref,
    major_path=str(REF_DIR / "scanvi_major_model"),
    fine_path=str(REF_DIR  / "scanvi_fine_model"),
    tag='REF'
)
adata_ref = compute_umaps_leiden(adata_ref, 'REF')

l3_cats_r = list(adata_ref.obs['cell_type_L3'].cat.categories)
adata_ref.uns['cell_type_L3_colors'] = [
    L3_COLORS.get(c, plt.cm.tab20(i % 20)) for i, c in enumerate(l3_cats_r)
]

# Free scVI model — no longer needed
del scvi_model; gc.collect()
torch.cuda.empty_cache() if torch.cuda.is_available() else None
print("[OK] REF scANVI complete  |  scVI model freed from memory")

# %% [markdown]
# ## Step 11: Visualizations
# 
# All figures generated for **SELF** branch (AT2 retained).
# To generate REF figures: set `VIZ_ADATA = adata_ref` and `VIZ_DIR = REF_DIR`.

# %%
VIZ_ADATA = adata_self
VIZ_DIR   = SELF_DIR
sc.settings.figdir = VIZ_DIR / "figures"
sc.settings.figdir.mkdir(exist_ok=True)

if VIZ_ADATA.raw is not None:
    gene_universe    = VIZ_ADATA.raw.var_names
    use_raw_for_plot = True
    print(f"[INFO] Using .raw for markers ({len(gene_universe):,} genes, log1p normalized)")
else:
    gene_universe    = VIZ_ADATA.var_names
    use_raw_for_plot = False
    print(f"[INFO] Using .X for markers ({len(gene_universe):,} genes)")

available_markers, missing_markers = get_available_markers(L3_MARKER_PANELS, gene_universe)
print(f"[INFO] Marker availability: {len(available_markers)}/{len(L3_MARKER_PANELS)} L3 types")
if 0 < len(missing_markers) < 30:
    print(f"  Missing: {', '.join(sorted(missing_markers))}")

# %% [markdown]
# ### Step 11-A: L3 Overview UMAP

# %%
print("\n" + "=" * 80)
print("STEP 11-A: L3 Overview UMAP")
print("=" * 80)

fig, ax = plt.subplots(figsize=(14, 10))
sc.pl.embedding(VIZ_ADATA, basis='X_umap_fine', color='cell_type_L3',
                ax=ax, show=False, legend_loc='right margin', legend_fontsize=8,
                frameon=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                title='Epithelial Cells - L3 Annotations (SELF)')
plt.tight_layout()
out = sc.settings.figdir / f'01_umap_L3_overview.{FIGURE_FORMAT}'
plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] Saved: {out.name}")

# Lineage-specific UMAPs — Basal group includes suprabasal subtypes
lineage_groups = {
    'Alveolar':  ['AT1_Canonical','AT1_MatrixRemodeling','AT2','AT2_Cycling'],
    'Basal':     ['Basal_Progenitor','Basal_Cycling','Basal_Inflammatory','Basal_EMT_ECM',
                  'Suprabasal_Progenitor','Suprabasal_Cycling'],
    'Ciliated':  ['Ciliated_Mature','Ciliogenesis_Deuterosomal','Ciliated_Cycling_Immature'],
    'Secretory': ['Goblet','Club','SMG_Mucous','Goblet_Defense_DUOX2'],
    'SMG':       ['SMG_Serous','SMG_Duct_Secretory_Defense'],
    'Other':     ['Squamous_Metaplasia','Ionocyte_Brush'],
}

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
for idx, (lin, labels) in enumerate(lineage_groups.items()):
    mask = VIZ_ADATA.obs['cell_type_L3'].isin(labels)
    col_ = f'_lin_{lin}'
    VIZ_ADATA.obs[col_] = 'Other'
    VIZ_ADATA.obs.loc[mask, col_] = VIZ_ADATA.obs.loc[mask, 'cell_type_L3']
    sc.pl.embedding(VIZ_ADATA, basis='X_umap_fine', color=col_,
                    ax=axes.flatten()[idx], show=False, frameon=False,
                    size=UMAP_SIZE*0.7, alpha=UMAP_ALPHA*0.8,
                    title=f'{lin} Lineage', legend_loc='right margin', legend_fontsize=6)
    VIZ_ADATA.obs.drop(columns=col_, inplace=True)
plt.tight_layout()
out = sc.settings.figdir / f'02_umap_L3_by_lineage.{FIGURE_FORMAT}'
plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] Saved: {out.name}")

# %% [markdown]
# ### Step 11-B: Marker Dotplots

# %%
print("\n" + "=" * 80)
print("STEP 11-B: Generating Marker Dotplots")
print("=" * 80)

lineage_marker_groups = {
    'Alveolar':     ['AT1_Canonical','AT1_MatrixRemodeling','AT2','AT2_Cycling'],
    'Basal':        ['Basal_Progenitor','Basal_Cycling','Basal_Inflammatory','Basal_EMT_ECM',
                     'Suprabasal_Progenitor','Suprabasal_Cycling'],
    'Ciliated':     ['Ciliated_Mature','Ciliogenesis_Deuterosomal','Ciliated_Cycling_Immature'],
    'Secretory_SMG':['Goblet','Club','SMG_Mucous','Goblet_Defense_DUOX2',
                     'SMG_Serous','SMG_Duct_Secretory_Defense'],
    'Special':      ['Squamous_Metaplasia','Ionocyte_Brush'],
}

for lineage, cell_types in lineage_marker_groups.items():
    print(f"\n  Processing {lineage}...")
    markers_for_lineage = list(dict.fromkeys(
        m for ct in cell_types if ct in available_markers
        for m in available_markers[ct]
    ))
    if not markers_for_lineage:
        print(f"    [WARN] No markers available for {lineage}"); continue

    adata_sub = VIZ_ADATA[VIZ_ADATA.obs['cell_type_L3'].isin(cell_types)].copy()
    if adata_sub.n_obs == 0:
        print(f"    [WARN] No cells found for {lineage}"); continue

    print(f"    Cells: {adata_sub.n_obs:,}  |  Markers: {len(markers_for_lineage)}")
    try:
        fig_w = max(12, len(markers_for_lineage) * 0.3)
        fig_h = max(6,  len(cell_types) * 0.4)
        dp = sc.pl.dotplot(adata_sub, var_names=markers_for_lineage, groupby='cell_type_L3',
                           use_raw=use_raw_for_plot, standard_scale='var', show=False,
                           return_fig=True, figsize=(fig_w, fig_h),
                           dendrogram=True, cmap='Reds', vmin=-2, vmax=2)
        out = sc.settings.figdir / f'03_dotplot_{lineage}.{FIGURE_FORMAT}'
        dp.savefig(str(out), dpi=FIGURE_DPI, bbox_inches='tight')
        plt.close('all')
        print(f"    [OK] Saved: {out.name}")
    except Exception as e:
        print(f"    [ERROR] Dotplot failed: {e}")
        traceback.print_exc()
    del adata_sub; gc.collect()

# %% [markdown]
# ### Step 11-C: Core Marker Heatmap

# %%
print("\n" + "=" * 80)
print("STEP 11-C: Generating Marker Expression Heatmap (log1p, from .raw)")
print("=" * 80)

all_core_markers = list(dict.fromkeys(m for ms in available_markers.values() for m in ms))
print(f"[INFO] Total core markers: {len(all_core_markers)}")

l3_types = sorted(VIZ_ADATA.obs['cell_type_L3'].unique())
expr_rows, row_names = [], []

for gene in all_core_markers:
    if use_raw_for_plot and gene in VIZ_ADATA.raw.var_names:
        e = VIZ_ADATA.raw[:, gene].X
    elif gene in VIZ_ADATA.var_names:
        e = VIZ_ADATA[:, gene].X
    else:
        continue
    if sparse.issparse(e):
        e = e.toarray().ravel()
    else:
        e = np.asarray(e).ravel()
    expr_rows.append([e[VIZ_ADATA.obs['cell_type_L3']==t].mean() for t in l3_types])
    row_names.append(gene)

if expr_rows:
    expression_df = pd.DataFrame(expr_rows, index=row_names, columns=l3_types)
    print(f"[OK] Expression matrix: {expression_df.shape}")

    fig, ax = plt.subplots(figsize=(max(12, len(l3_types)*0.4),
                                     max(10, len(row_names)*0.15)))
    sns.heatmap(expression_df, cmap='RdYlBu_r', center=0, robust=True,
                yticklabels=True, xticklabels=True,
                cbar_kws={'label': 'Mean Expression (log1p normalized)'},
                linewidths=0.1, linecolor='lightgray', ax=ax)
    plt.title('Core Marker Expression Across L3 Cell Types', fontsize=14, pad=20)
    plt.xlabel('L3 Cell Type', fontsize=12); plt.ylabel('Marker Gene', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=9); plt.yticks(fontsize=8)
    plt.tight_layout()
    out = sc.settings.figdir / f'04_heatmap_core_markers.{FIGURE_FORMAT}'
    plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"[OK] Saved: {out.name}")
    del expression_df, expr_rows; gc.collect()

# %% [markdown]
# ### Step 11-D: Marker Expression UMAPs

# %%
print("\n" + "=" * 80)
print("STEP 11-D: Marker Expression UMAPs")
print("=" * 80)

representative_markers = {
    'Alveolar':   ['AGER','HOPX','SFTPC','NAPSA','CHI3L1','MKI67'],
    'Basal':      ['KRT5','TP63','KRT14','KRT17','NGFR'],
    'Suprabasal': ['KRT4','KRT13','NOTCH1','NOTCH3','MKI67'],
    'Ciliated':   ['FOXJ1','DNAH5','DEUP1','CCNO'],
    'Secretory':  ['SCGB1A1','MUC5AC','MUC5B','DUOX2','SPDEF'],
    'SMG':        ['LTF','LYZ','PIGR','WFDC2'],
    'Special':    ['SPRR2A','FOXI1','CFTR','TOP2A','VIM'],
}

for category, markers in representative_markers.items():
    available = [m for m in markers if m in gene_universe]
    if not available: continue
    n_cols  = min(3, len(available))
    n_rows  = int(np.ceil(len(available) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes_flat = np.array(axes).flatten()
    for i, marker in enumerate(available):
        sc.pl.embedding(VIZ_ADATA, basis='X_umap_fine', color=marker,
                        ax=axes_flat[i], show=False, use_raw=use_raw_for_plot,
                        vmax='p99', frameon=False, size=UMAP_SIZE*0.8,
                        alpha=UMAP_ALPHA, cmap='Reds', title=marker)
    for i in range(len(available), len(axes_flat)):
        axes_flat[i].axis('off')
    plt.suptitle(f'{category} Markers', fontsize=14, y=1.02)
    plt.tight_layout()
    out = sc.settings.figdir / f'05_umap_markers_{category}.{FIGURE_FORMAT}'
    plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"[OK] Saved: {out.name}  ({', '.join(available)})")

# %% [markdown]
# ### Step 11-E: Cell Type Proportion Analysis

# %%
print("\n" + "=" * 80)
print("STEP 11-E: Cell Type Proportion Analysis")
print("=" * 80)

l3_counts_viz = VIZ_ADATA.obs['cell_type_L3'].value_counts()
l3_props = (l3_counts_viz / VIZ_ADATA.n_obs * 100).sort_values(ascending=True)

fig, ax = plt.subplots(figsize=(10, max(8, len(l3_props)*0.35)))
colors = [L3_COLORS.get(lb, '#cccccc') for lb in l3_props.index]
l3_props.plot(kind='barh', ax=ax, color=colors)
ax.set_xlabel('Percentage of Cells', fontsize=12)
ax.set_ylabel('L3 Cell Type', fontsize=12)
ax.set_title('L3 Cell Type Distribution (SELF)', fontsize=14)
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
out = sc.settings.figdir / f'06_barplot_L3_proportions.{FIGURE_FORMAT}'
plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()
print(f"[OK] Saved: {out.name}")

if 'dataset' in VIZ_ADATA.obs.columns:
    prop_df = pd.crosstab(VIZ_ADATA.obs['dataset'], VIZ_ADATA.obs['cell_type_L3'],
                           normalize='index') * 100
    top_l3  = l3_counts_viz.head(15).index
    fig, ax = plt.subplots(figsize=(14, 8))
    prop_df[top_l3].plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_ylabel('Percentage'); ax.set_xlabel('Dataset')
    ax.set_title('L3 Cell Type Distribution Across Datasets (Top 15)', fontsize=14)
    ax.legend(title='L3 Type', bbox_to_anchor=(1.05,1), loc='upper left', fontsize=8)
    plt.xticks(rotation=45, ha='right'); plt.tight_layout()
    out = sc.settings.figdir / f'07_stacked_bar_L3_by_dataset.{FIGURE_FORMAT}'
    plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"[OK] Saved: {out.name}")
    prop_df.to_csv(VIZ_DIR / "L3_proportions_by_dataset.csv", index=False)

# %% [markdown]
# ## Step 12: Save Outputs

# %%
print("\n" + "=" * 80)
print("STEP 12: Saving Outputs")
print("=" * 80)

for branch, adata_b, out_dir in [('SELF', adata_self, SELF_DIR),
                                   ('REF',  adata_ref,  REF_DIR)]:
    print(f"\n[{branch}] Saving...")
    h5ad_path = out_dir / f"epithelial_scanvi_v2_7_HOTFIX_{branch}_final.h5ad"
    adata_b.write_h5ad(str(h5ad_path), compression='gzip', compression_opts=9)
    size_gb = h5ad_path.stat().st_size / 1e9
    print(f"  h5ad: {h5ad_path.name} ({size_gb:.2f} GB)  shape: {adata_b.shape}")

    # P2-4 fix: index=False on all CSV exports
    adata_b.obs[['subcluster','cell_type_L3',
                  'scanvi_fine_pred','scanvi_major_pred',
                  'scanvi_fine_conf','scanvi_major_conf']].to_csv(
        out_dir / "L3_scanvi_annotations.csv", index=False
    )

    fs = (adata_b.obs['scanvi_fine_pred'].value_counts()
          .rename_axis('Cell_Type').reset_index(name='Count'))
    fs['Percentage'] = 100 * fs['Count'] / fs['Count'].sum()
    conf = adata_b.obs.groupby('scanvi_fine_pred')['scanvi_fine_conf'].agg(['mean','std'])
    fs = fs.merge(conf, left_on='Cell_Type', right_index=True, how='left')
    fs.to_csv(out_dir / "scanvi_fine_statistics.csv", index=False)

    with open(out_dir / "annotation_summary.txt", "w") as f:
        f.write("=" * 80 + "\n")
        f.write(f"Epithelial L3 scANVI Annotation Summary — {branch}\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input: {INPUT_H5AD}\n")
        f.write(f"scVI model: {SCVI_MODEL_PATH}\n\n")
        f.write(f"Cells: {adata_b.n_obs:,}  |  Genes (HVG): {adata_b.n_vars:,}\n")
        f.write(f"L3 types: {adata_b.obs['cell_type_L3'].nunique()}\n\n")
        f.write("L3 Cell Type Distribution:\n")
        for lbl, cnt in adata_b.obs['cell_type_L3'].value_counts().items():
            f.write(f"  {lbl}: {cnt:,} ({cnt/adata_b.n_obs*100:.2f}%)\n")

    print(f"  [OK] Outputs saved for {branch}")

# P1-5 fix: REF manifest.json — documents all model dependencies
ref_manifest = {
    "branch":                "REF",
    "version":               "2.7-HOTFIX",
    "created_at":            datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    "input_h5ad":            INPUT_H5AD,
    "scvi_model_path":       SCVI_MODEL_PATH,
    "hvg_file":              HVG_FILE_PATH,
    "scanvi_major_model":    str(REF_DIR / "scanvi_major_model"),
    "scanvi_fine_model":     str(REF_DIR / "scanvi_fine_model"),
    "scanvi_labels_key":     "scanvi_labels_fine",
    "at_labels_removed":     sorted(AT_LABELS_TO_REMOVE),
    "n_cells":               int(adata_ref.n_obs),
    "n_hvg":                 int(adata_ref.n_vars),
    "n_l3_types":            int(adata_ref.obs['cell_type_L3'].nunique()),
    "note": ("scVI model is shared with SELF branch. "
             "REF scANVI was trained on AT-removed subset of the same scVI latent space. "
             "Use SCVI.prepare_query_anndata + SCANVI.load_query_data (pass path strings) "
             "for scArches query mapping — do NOT call .load() explicitly.")
}
with open(REF_DIR / "manifest.json", "w") as f:
    json.dump(ref_manifest, f, indent=2)
print(f"[OK] REF manifest.json saved: {REF_DIR}/manifest.json")

elapsed = (time.time() - PIPELINE_START) / 60
print(f"\n{'='*80}")
print(f"PIPELINE COMPLETE  |  Elapsed: {elapsed:.1f} min")
print(f"{'='*80}")
print(f"Shared scVI:     {SCVI_MODEL_PATH}")
print(f"SELF final h5ad: {SELF_DIR}/epithelial_scanvi_v2_7_HOTFIX_SELF_final.h5ad")
print(f"REF  final h5ad: {REF_DIR}/epithelial_scanvi_v2_7_HOTFIX_REF_final.h5ad")
print(f"REF  manifest:   {REF_DIR}/manifest.json")
print("=" * 80)

# %% [markdown]
# ## Step 13: Run Logs and AnnData Structure Summaries

# %%
print("\n" + "=" * 80)
print("STEP 13: Writing Run Logs and AnnData Structure Summaries")
print("=" * 80)


def _fmt_shape(obj):
    s = getattr(obj, 'shape', None)
    return 'NA' if s is None else ' x '.join(str(x) for x in s)

def _matrix_line(name, obj):
    return (f"- {name}: type={type(obj).__name__}, shape={_fmt_shape(obj)}, "
            f"dtype={getattr(obj,'dtype','NA')}, sparse={sparse.issparse(obj) if obj is not None else False}")

def _mapping_lines(title, mapping):
    lines = [f"{title} ({len(mapping)}):"]
    if not mapping:
        lines.append("  (none)"); return lines
    for key, value in mapping.items():
        lines.append(f"  - {key}: type={type(value).__name__}, shape={getattr(value,'shape','NA')}, "
                     f"dtype={getattr(value,'dtype','NA')}")
    return lines

def _df_column_lines(df, title):
    lines = [f"{title} ({len(df.columns)} columns):"]
    if not len(df.columns):
        lines.append("  (none)"); return lines
    for col in df.columns:
        lines.append(f"  - {col}: dtype={df[col].dtype}, non_null={df[col].notna().sum():,}, "
                     f"unique={df[col].nunique(dropna=False):,}")
    return lines

def _write_branch_logs(branch, adata_b, out_dir):
    h5ad_path      = out_dir / f"epithelial_scanvi_v2_7_HOTFIX_{branch}_final.h5ad"
    run_log_path   = out_dir / "pipeline_run_log.txt"
    structure_path = out_dir / "anndata_structure.txt"
    elapsed_min    = (time.time() - PIPELINE_START) / 60

    run_lines = [
        "=" * 80,
        f"Epithelial scVI/scANVI Pipeline Run Log — {branch} (v2.7-HOTFIX)",
        "=" * 80,
        f"written_at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"input_h5ad: {INPUT_H5AD}",
        f"saved_h5ad: {h5ad_path}",
        f"saved_h5ad_exists: {h5ad_path.exists()}",
        f"cells: {adata_b.n_obs:,}",
        f"genes_hvg: {adata_b.n_vars:,}",
        f"unique_L3_labels: {adata_b.obs['cell_type_L3'].nunique():,}",
        f"batches: {adata_b.obs[BATCH_KEY].nunique():,}",
        f"unknown_major: {(adata_b.obs['scanvi_labels_major']==UNLABELED_CATEGORY).sum():,}",
        f"unknown_fine: {(adata_b.obs['scanvi_labels_fine']==UNLABELED_CATEGORY).sum():,}",
        f"global_stress_threshold_p{STRESS_PERCENTILE}: {GLOBAL_STRESS_THRESHOLD:.4f}",
        f"elapsed_minutes_at_write: {elapsed_min:.2f}",
        f"figure_dir: {out_dir/'figures'}",
        f"major_model_dir: {out_dir/'scanvi_major_model'}",
        f"fine_model_dir: {out_dir/'scanvi_fine_model'}",
    ]
    if h5ad_path.exists():
        run_lines.append(f"saved_h5ad_size_gb: {h5ad_path.stat().st_size / 1e9:.3f}")
    run_lines += [
        "", "Predicted fine labels:",
        adata_b.obs['scanvi_fine_pred'].value_counts().to_string(),
        "", "Predicted major labels:",
        adata_b.obs['scanvi_major_pred'].value_counts().to_string(),
    ]
    run_log_path.write_text('\n'.join(run_lines) + '\n', encoding='utf-8')

    structure_lines = [
        "=" * 80, f"AnnData Structure Summary — {branch}", "=" * 80,
        f"written_at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"shape: {adata_b.shape}",
        _matrix_line('X', adata_b.X), "",
    ]
    structure_lines.extend(_df_column_lines(adata_b.obs, 'obs')); structure_lines.append("")
    structure_lines.extend(_df_column_lines(adata_b.var, 'var')); structure_lines.append("")
    structure_lines.extend(_mapping_lines('layers', adata_b.layers)); structure_lines.append("")
    structure_lines.extend(_mapping_lines('obsm', adata_b.obsm)); structure_lines.append("")
    structure_lines.extend(_mapping_lines('obsp', adata_b.obsp)); structure_lines.append("")
    structure_lines.extend(_mapping_lines('uns', adata_b.uns)); structure_lines.append("")
    if adata_b.raw is not None:
        structure_lines.append(f"raw: shape={adata_b.raw.shape} (log1p normalized)")
        raw_cols = ', '.join(adata_b.raw.var.columns.astype(str)) if len(adata_b.raw.var.columns) else '(none)'
        structure_lines.append(f"raw.var columns: {raw_cols}")
    else:
        structure_lines.append("raw: None")

    structure_path.write_text('\n'.join(structure_lines) + '\n', encoding='utf-8')
    print(f"[{branch}] Wrote {run_log_path.name} and {structure_path.name}")


for branch, adata_b, out_dir in [
    ('SELF', adata_self, SELF_DIR),
    ('REF',  adata_ref,  REF_DIR),
]:
    _write_branch_logs(branch, adata_b, out_dir)

print("[OK] All run logs and structure summaries saved")

# %%
import inspect
print('scvi-tools version:', scvi.__version__)
print('SCANVI.__init__:', inspect.signature(scvi.model.SCANVI.__init__))
print('SCANVI.from_scvi_model:', inspect.signature(scvi.model.SCANVI.from_scvi_model))


