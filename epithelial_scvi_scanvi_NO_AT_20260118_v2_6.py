#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Epithelial scVI-scANVI DUAL Pipeline v2.6.2 - WITH MODEL LOADING
==================================================================

Complete Implementation of Allcells v2.3.1 Best Practices + Model Loading:
- ✅ P0 Fixes (stats alignment, n_latent, Unknown, HVG save, expanded markers)
- ✅ P0-5 Fix (Categorical Unknown assignment after .copy())
- ✅ Covariates System (mt%, stress, cell cycle as continuous covariates)
- ✅ Unknown Cleaning (neighborhood purity + quality gating)
- ✅ Whitelist Filtering (from downstream h5ad)
- ✅ scArches-ready parameters
- ✅ Class imbalance handling (n_samples_per_label)
- ⭐ NEW: Model Loading (skip retraining if models exist)
- ✅ HOTFIX: Path object type consistency
- ⭐ NEW: HVG loading from saved model (ensures compatibility)

Reference: allcells_scvi_celltypist_scanvi_pipeline_20260115_v2_3_1.py
Based on: epithelial_scvi_scanvi_NO_AT_20260115_v2_5_FIXED.py

Author: r2end
Date: 2026-01-18
Version: 2.6.2-PRODUCTION - Model loading feature + path hotfix
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from scipy import sparse
from pathlib import Path
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings('ignore')

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

print("="*80)
print("Epithelial scVI-scANVI Dual Training Pipeline (v2.6.2-PRODUCTION)")
print("✅ Complete Allcells v2.3.1 Architecture + Model Loading")
print("="*80)
print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"scanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"PyTorch: {torch.__version__}")
print(f"Device: {'GPU' if torch.cuda.is_available() else 'CPU'}")
print("="*80)

# ===== CONFIGURATION =====

# Input/Output
INPUT_H5AD = "/home/h2048/data/py/0115/epithelial_DUAL_SCANVI_v2_7_PRODUCTION/checkpoints/adata_preprocessed.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0115/epithelial_DUAL_SCANVI_v2_7_PRODUCTION"
CHECKPOINT_DIR = f"{OUTPUT_DIR}/checkpoints"

# Whitelist filtering (from downstream h5ad)
DOWNSTREAM_H5AD_DIR = "/home/h2048/data/py/1217/downstream_analysis_v1_6_2"
DOWNSTREAM_H5AD_GLOB = "**/*epithelial*.h5ad"  # Only epithelial downstream files
CELLID_SET_MODE = "union"  # "union" or "intersection"
ENABLE_WHITELIST_FILTERING = False

# ===== Model Loading Configuration (NEW in v2.6) =====
# Set to True to skip retraining if valid model exists
# Useful for debugging downstream steps without full retraining
# Time savings: ~70-105 min training → ~70 sec loading

LOAD_SCVI_IF_EXISTS = True              # Load scVI if model.pt exists
LOAD_SCANVI_MAJOR_IF_EXISTS = False      # Load scANVI-major if model.pt exists  
LOAD_SCANVI_FINE_IF_EXISTS = False       # Load scANVI-fine if model.pt exists

# Note: HVG genes must match for model loading to work
# If HVG selection changes, existing models become incompatible

# Create directories
for dir_path in [OUTPUT_DIR, f"{OUTPUT_DIR}/figures", f"{OUTPUT_DIR}/scvi_model", 
                 f"{OUTPUT_DIR}/scanvi_major_model", f"{OUTPUT_DIR}/scanvi_fine_model", 
                 CHECKPOINT_DIR]:
    os.makedirs(dir_path, exist_ok=True)

# Keys
BATCH_KEY = 'sample'
CELLTYPE_KEY = 'celltypist_pred'
UNLABELED_CATEGORY = 'Unknown'

USE_HIERARCHICAL_LABELS = True

# Hierarchical labels
MAJOR_LINEAGE_MAP = {
    'Basal': 'Basal_Lineage',
    'Suprabasal': 'Basal_Lineage',
    'SMG_Basal': 'Basal_Lineage',
    'Dividing_Basal': 'Basal_Lineage',
    
    'Ciliated': 'Ciliated_Lineage',
    'Deuterosome': 'Ciliated_Lineage',
    'Deuterosomal': 'Ciliated_Lineage',
    
    'Secretory_Goblet': 'Secretory_Lineage',
    'Secretory_Club': 'Secretory_Lineage',
    'SMG_Mucous': 'Secretory_Lineage',
    'SMG_Serous': 'Secretory_Lineage',
    'SCGB1A1+': 'Secretory_Lineage',
    
    'SMG_Duct': 'Duct',
    
    'Ionocyte_n_Brush': 'Rare_Specialized',
    'Ionocyte': 'Rare_Specialized',
    'Brush': 'Rare_Specialized'
}

# Epithelial-specific forced markers
FORCE_INCLUDE_MARKERS = [
    # Core epithelial identity
    'EPCAM', 'CDH1', 'ELF3', 'FXYD3', 'CLDN4', 'CLDN7',
    
    # Basal lineage
    'TP63', 'KRT5', 'KRT14', 'KRT15', 'KRT17',
    'ITGA6', 'ITGB4', 'COL17A1',
    'KRT13', 'KRT4',
    
    # Ciliated lineage
    'FOXJ1', 'RSPH1', 'PIFO', 'DNAH5', 'DNAH11',
    'TUBA1A', 'TUBB4B',
    'DEUP1', 'CCNO', 'CEP78',
    
    # Secretory lineage
    'MUC5AC', 'MUC5B', 'MUC4', 'TFF3', 'SPDEF', 'AGR2',
    'SCGB1A1', 'SCGB3A1', 'SCGB3A2', 'BPIFA1', 'WFDC2',
    'LYZ', 'LTF', 'DMBT1',
    'MUC7', 'PRH1',
    
    # Duct
    'KRT7', 'KRT19', 'KRT8', 'KRT18', 'AQP3',
    
    # Rare/Specialized
    'CFTR', 'FOXI1', 'ATP6V1B1',
    'DCLK1', 'TRPM5', 'POU2F3', 'GFI1B',
    
    # Differentiation/Transition
    'TP73', 'NOTCH1', 'NOTCH2', 'DLL1', 'HES1',
    'WNT5A', 'WNT7B',
    'SOX2', 'SOX9',
]

# â­ NEW: Stress signature genes (from allcells line 165-179)
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

# â­ NEW: Cell cycle genes (from allcells line 182-197)
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

# â­ NEW: Unknown cleaning thresholds (from allcells line 203-206)
EXISTING_LABEL_PURITY_THRESHOLD = 0.5  # Neighborhood purity
EXISTING_LABEL_MT_THRESHOLD = 20  # mt% cutoff
EXISTING_LABEL_STRESS_PERCENTILE = 95  # stress percentile cutoff

# Model parameters (unified n_latent)
N_HVG = 4000
N_LATENT_SCVI = 100
N_LATENT_SCANVI = 100  # Must match scVI
N_LAYERS = 2  # scArches-ready
N_HIDDEN = 256
DROPOUT_RATE = 0.2  # scArches-ready
GENE_LIKELIHOOD = "nb"
DISPERSION = "gene-batch"

# â­ NEW: scArches-ready parameters (from allcells line 218-221)
SCVI_ENCODE_COVARIATES = True
SCVI_USE_LAYER_NORM = "both"
SCVI_USE_BATCH_NORM = "none"

# Training
SCVI_MAX_EPOCHS = 400
SCANVI_MAX_EPOCHS = 200
EARLY_STOPPING = True
BATCH_SIZE = 1024

# â­ NEW: Class imbalance handling (from allcells line 228)
SCANVI_N_SAMPLES_PER_LABEL = 2000

# QC thresholds
MIN_CELLS_PER_BATCH = 30
LOW_CONFIDENCE_THRESHOLD = 0.3

SAVE_CHECKPOINTS = True
MARK_LOW_CONFIDENCE_AS_UNKNOWN = True

# Random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# ===== Model Path Configuration (v2.6 NEW) =====
# Convert to Path objects for easier existence checking

scvi_model_path = Path(f"{OUTPUT_DIR}/scvi_model")
scanvi_major_model_path = Path(f"{OUTPUT_DIR}/scanvi_major_model")
scanvi_fine_model_path = Path(f"{OUTPUT_DIR}/scanvi_fine_model")
hvg_file_path = scvi_model_path / "hvg_genes.txt"

def check_model_exists(model_dir):
    """Check if a valid scVI/scANVI model exists"""
    if not isinstance(model_dir, Path):
        model_dir = Path(model_dir)
    return model_dir.exists() and (model_dir / "model.pt").exists()

# Model existence flags (will be set after HVG selection)
scvi_model_exists = False
scanvi_major_exists = False
scanvi_fine_exists = False

print("\n" + "="*80)
print("CONFIGURATION")
print("="*80)
print(f"Input: {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")

print(f"\nâœ… PRODUCTION FEATURES:")
print(f"  1. P0 fixes (stats, n_latent, Unknown, HVG save, markers)")
print(f"  2. Whitelist filtering: {ENABLE_WHITELIST_FILTERING}")
print(f"  3. Covariates: mt%, stress, S_score, G2M_score")
print(f"  4. Unknown cleaning: purity + quality gating")
print(f"  5. scArches-ready: encode_cov=True, n_layers=2")
print(f"  6. Class balance: n_samples_per_label={SCANVI_N_SAMPLES_PER_LABEL}")

print(f"\nâ­ Dual scANVI Training:")
print(f"  scANVI_major: {len(set(MAJOR_LINEAGE_MAP.values()))} major lineages")
print(f"  scANVI_fine: {len(MAJOR_LINEAGE_MAP)} fine types")
print("="*80)


# ===== UTILITY FUNCTIONS =====

def _read_obs_index_from_h5ad(h5ad_path: str):
    """Read obs_names from h5ad without loading the whole AnnData."""
    try:
        import h5py
        with h5py.File(h5ad_path, "r") as f:
            if "obs" in f and "_index" in f["obs"]:
                arr = f["obs"]["_index"][()]
                if hasattr(arr, "dtype") and arr.dtype.kind in ("S", "O"):
                    return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in arr]
                return [str(x) for x in arr]
    except Exception as e:
        print(f"  Warning: Could not read {h5ad_path}: {e}")
    return []


# ===== STEP 0: Whitelist Filtering =====

print("\n" + "="*80)
print("STEP 0: Whitelist Filtering (from Downstream H5AD)")
print("="*80)

whitelist_cells = None

if ENABLE_WHITELIST_FILTERING:
    print(f"\nCollecting cell IDs from: {DOWNSTREAM_H5AD_DIR}")
    print(f"Pattern: {DOWNSTREAM_H5AD_GLOB}")
    
    import glob
    h5ad_files = glob.glob(f"{DOWNSTREAM_H5AD_DIR}/{DOWNSTREAM_H5AD_GLOB}", recursive=True)
    
    print(f"Found {len(h5ad_files)} h5ad files")
    
    cell_id_sets = []
    for h5ad_file in h5ad_files:
        cell_ids = _read_obs_index_from_h5ad(h5ad_file)
        if cell_ids:
            cell_id_sets.append(set(cell_ids))
            print(f"  {Path(h5ad_file).name}: {len(cell_ids)} cells")
    
    if cell_id_sets:
        if CELLID_SET_MODE == "union":
            whitelist_cells = set.union(*cell_id_sets)
            print(f"\nâœ“ Whitelist (union): {len(whitelist_cells)} cells")
        elif CELLID_SET_MODE == "intersection":
            whitelist_cells = set.intersection(*cell_id_sets)
            print(f"\nâœ“ Whitelist (intersection): {len(whitelist_cells)} cells")
        else:
            raise ValueError(f"Invalid CELLID_SET_MODE: {CELLID_SET_MODE}")
    else:
        print("\nâš ï¸  No cell IDs collected, skipping whitelist filtering")
        ENABLE_WHITELIST_FILTERING = False
else:
    print("\nâš ï¸  Whitelist filtering disabled")


# ===== STEP 1: Load Data =====

print("\n" + "="*80)
print("STEP 1: Loading Data")
print("="*80)

adata = sc.read_h5ad(INPUT_H5AD)
print(f"Initial: {adata.shape}")
print(f"Layers: {list(adata.layers.keys())}")

# Apply whitelist
if ENABLE_WHITELIST_FILTERING and whitelist_cells:
    before_n = adata.n_obs
    adata = adata[adata.obs_names.isin(whitelist_cells)].copy()
    after_n = adata.n_obs
    removed_n = before_n - after_n
    
    print(f"\nâœ“ Whitelist filtered:")
    print(f"  Before: {before_n}")
    print(f"  After: {after_n}")
    print(f"  Removed: {removed_n} ({100*removed_n/before_n:.1f}%)")

if 'counts' not in adata.layers:
    raise ValueError("âŒ 'counts' layer required")


# ===== STEP 1.5: Remove AT1/AT2 =====

print("\n" + "="*80)
print("STEP 1.5: Filtering AT1/AT2")
print("="*80)

at_cell_patterns = ['AT1', 'AT2', 'Alveolar']
at_cells_mask = adata.obs[CELLTYPE_KEY].str.contains('|'.join(at_cell_patterns), case=False, na=False)
n_at_cells = at_cells_mask.sum()

print(f"AT1/AT2 cells: {n_at_cells} ({100*n_at_cells/adata.n_obs:.2f}%)")

if n_at_cells > 0:
    adata = adata[~at_cells_mask].copy()
    print(f"âœ“ Removed {n_at_cells} cells, remaining: {adata.n_obs}")
    
    batch_counts = adata.obs[BATCH_KEY].value_counts()
    small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index
    
    if len(small_batches) > 0:
        print(f"Removing {len(small_batches)} small batches (<{MIN_CELLS_PER_BATCH} cells)")
        adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()


# ===== STEP 2: Data Preprocessing & Covariates =====

print("\n" + "="*80)
print("STEP 2: Preprocessing & Covariates Calculation")
print("="*80)

# 2.1 Validate counts
print("\n2.1 Validating counts layer...")
if not sparse.issparse(adata.layers['counts']):
    adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])

if np.any(np.isinf(adata.layers['counts'].data)) or np.any(np.isnan(adata.layers['counts'].data)):
    adata.layers['counts'].data[np.isinf(adata.layers['counts'].data)] = 0
    adata.layers['counts'].data[np.isnan(adata.layers['counts'].data)] = 0
    adata.layers['counts'].eliminate_zeros()

# 2.2 Create log1p layer
if 'log1p' not in adata.layers:
    adata_temp = adata.copy()
    adata_temp.X = adata_temp.layers['counts'].copy()
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    adata.layers['log1p'] = adata_temp.X.copy()
    del adata_temp
    gc.collect()

adata.X = adata.layers['log1p'].copy()

# â­ 2.3 MT% recalculation (from allcells line 664-696)
print("\nâ­ 2.3 MT% calculation...")

if 'pct_counts_mt' in adata.obs.columns:
    zero_mt_pct = (adata.obs['pct_counts_mt'] == 0).sum() / adata.n_obs
    print(f"  Existing pct_counts_mt: {zero_mt_pct*100:.1f}% are zero")
    
    if zero_mt_pct > 0.1:
        print("  âš ï¸  High zero rate, recalculating...")
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], layer='counts', 
                                   percent_top=None, inplace=True)
        print(f"  âœ“ Recalculated, mean: {adata.obs['pct_counts_mt'].mean():.2f}%")
else:
    print("  Computing pct_counts_mt...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], layer='counts',
                               percent_top=None, inplace=True)
    print(f"  âœ“ Computed, mean: {adata.obs['pct_counts_mt'].mean():.2f}%")

# â­ 2.4 Stress score (from allcells line 708-734)
print("\nâ­ 2.4 Calculating stress_score...")

adata_temp = sc.AnnData(X=adata.layers['counts'].copy(), var=adata.var.copy())
sc.pp.normalize_total(adata_temp, target_sum=1e4)
sc.pp.log1p(adata_temp)

stress_genes_in_data = [g for g in STRESS_SIGNATURE_GENES if g in adata_temp.var_names]
print(f"  Found {len(stress_genes_in_data)}/{len(STRESS_SIGNATURE_GENES)} stress genes")

if len(stress_genes_in_data) > 10:
    stress_gene_idx = [list(adata_temp.var_names).index(g) for g in stress_genes_in_data]
    if sparse.issparse(adata_temp.X):
        stress_expr = np.asarray(adata_temp.X[:, stress_gene_idx].mean(axis=1)).flatten()
    else:
        stress_expr = adata_temp.X[:, stress_gene_idx].mean(axis=1)
    adata.obs['stress_score'] = stress_expr
    print(f"  âœ“ Mean: {adata.obs['stress_score'].mean():.3f}, Median: {adata.obs['stress_score'].median():.3f}")
else:
    adata.obs['stress_score'] = 0.0
    print("  âš ï¸  Too few stress genes, set to 0")

del adata_temp
gc.collect()

# â­ 2.5 Cell cycle scores (from allcells line 737-770)
print("\nâ­ 2.5 Calculating cell cycle scores...")

s_genes = [x for x in S_GENES if x in adata.var_names]
g2m_genes = [x for x in G2M_GENES if x in adata.var_names]

print(f"  S genes: {len(s_genes)}/{len(S_GENES)}")
print(f"  G2M genes: {len(g2m_genes)}/{len(G2M_GENES)}")

if len(s_genes) > 10 and len(g2m_genes) > 10:
    adata_temp = sc.AnnData(X=adata.layers['counts'].copy(), var=adata.var.copy())
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    
    sc.tl.score_genes_cell_cycle(adata_temp, s_genes=s_genes, g2m_genes=g2m_genes)
    adata.obs['S_score'] = adata_temp.obs['S_score'].values
    adata.obs['G2M_score'] = adata_temp.obs['G2M_score'].values
    adata.obs['phase'] = adata_temp.obs['phase'].values
    
    del adata_temp
    gc.collect()
    
    print(f"  âœ“ S_score mean: {adata.obs['S_score'].mean():.3f}")
    print(f"  âœ“ G2M_score mean: {adata.obs['G2M_score'].mean():.3f}")
    print(f"  Phase distribution:")
    print(adata.obs['phase'].value_counts())
else:
    adata.obs['S_score'] = 0.0
    adata.obs['G2M_score'] = 0.0
    adata.obs['phase'] = 'G1'
    print("  âš ï¸  Too few cell cycle genes, set to 0")

print("\nâœ“ All covariates calculated")


# ===== STEP 3: Hierarchical Labels =====

print("\n" + "="*80)
print("STEP 3: Hierarchical Labels")
print("="*80)

adata.obs['fine_type'] = adata.obs[CELLTYPE_KEY].copy()
adata.obs['major_lineage'] = adata.obs['fine_type'].map(MAJOR_LINEAGE_MAP)

unmapped = adata.obs['major_lineage'].isna()
if unmapped.any():
    unmapped_types = adata.obs.loc[unmapped, 'fine_type'].unique()
    print(f"âš ï¸  Unmapped types: {len(unmapped_types)}")
    adata.obs.loc[unmapped, 'major_lineage'] = 'Other'

print("\nMajor lineages:")
print(adata.obs['major_lineage'].value_counts())

# Dual label columns
adata.obs['scanvi_labels_major'] = adata.obs['major_lineage'].copy()
adata.obs['scanvi_labels_fine'] = adata.obs['fine_type'].copy()

# Force Unknown into categories
print("\n3.1 Adding Unknown to label categories...")
for k in ['scanvi_labels_major', 'scanvi_labels_fine']:
    adata.obs[k] = adata.obs[k].astype('category')
    if UNLABELED_CATEGORY not in adata.obs[k].cat.categories:
        adata.obs[k] = adata.obs[k].cat.add_categories([UNLABELED_CATEGORY])

# Mark low confidence as Unknown (initial)
if MARK_LOW_CONFIDENCE_AS_UNKNOWN and 'celltypist_conf_score' in adata.obs.columns:
    low_conf_mask = adata.obs['celltypist_conf_score'] < LOW_CONFIDENCE_THRESHOLD
    n_low_conf = low_conf_mask.sum()
    print(f"\n3.2 Low confidence (<{LOW_CONFIDENCE_THRESHOLD}): {n_low_conf} ({100*n_low_conf/adata.n_obs:.1f}%)")
    
    adata.obs.loc[low_conf_mask, 'scanvi_labels_major'] = UNLABELED_CATEGORY
    adata.obs.loc[low_conf_mask, 'scanvi_labels_fine'] = UNLABELED_CATEGORY


# ===== STEP 4: HVG Selection =====

print("\n" + "="*80)
print("STEP 4: HVG Selection")
print("="*80)

# Check if we should load existing HVG list (for model compatibility)
hvg_loaded_from_file = False

if LOAD_SCVI_IF_EXISTS and hvg_file_path.exists():
    print("\n⭐ Loading HVG list from existing model...")
    try:
        saved_hvg = pd.read_csv(hvg_file_path, header=None)[0].tolist()
        
        # Verify all genes exist in current data
        available_saved_hvg = [g for g in saved_hvg if g in adata.var_names]
        missing_genes = [g for g in saved_hvg if g not in adata.var_names]
        
        if len(missing_genes) > 0:
            print(f"⚠️  {len(missing_genes)} genes from saved HVG not in current data")
            print(f"   Will select new HVG to ensure consistency")
        else:
            # Mark saved genes as highly variable
            adata.var['highly_variable'] = False
            adata.var.loc[available_saved_hvg, 'highly_variable'] = True
            
            hvg_method = f"loaded from file ({len(available_saved_hvg)} genes)"
            hvg_loaded_from_file = True
            
            print(f"✓ Loaded {len(available_saved_hvg)} HVG from: {hvg_file_path}")
            print(f"  Matches saved model genes: {len(saved_hvg)}")
            
    except Exception as e:
        print(f"⚠️  Failed to load HVG list: {e}")
        print(f"   Will select new HVG")

# If not loaded from file, select new HVG
if not hvg_loaded_from_file:
    print("\n⭐ Selecting new HVG...")
    try:
        sc.pp.highly_variable_genes(
            adata,
            layer='counts',
            n_top_genes=N_HVG,
            batch_key=BATCH_KEY,
            flavor='seurat_v3',
            subset=False
        )
        hvg_method = "batch-aware (seurat_v3)"
    except:
        sc.pp.highly_variable_genes(
            adata,
            layer='counts',
            n_top_genes=N_HVG,
            flavor='seurat_v3',
            subset=False
        )
        hvg_method = "non-batch-aware (seurat_v3)"
    
    print(f"✓ Selected {adata.var['highly_variable'].sum()} HVG")

adata.uns['hvg_method'] = hvg_method
print(f"✓ HVG method: {hvg_method}")
# Force-include markers
print("\n4.1 Force-including epithelial markers...")
available_markers = [g for g in FORCE_INCLUDE_MARKERS if g in adata.var_names]
missing_markers = [g for g in FORCE_INCLUDE_MARKERS if g not in adata.var_names]

if missing_markers:
    print(f"  Missing: {len(missing_markers)}")

n_force_added = 0
for gene in available_markers:
    if not adata.var.loc[gene, 'highly_variable']:
        adata.var.loc[gene, 'highly_variable'] = True
        n_force_added += 1

print(f"âœ“ Force-included {len(available_markers)} markers ({n_force_added} new)")
print(f"  Final HVG: {adata.var['highly_variable'].sum()}")

# Preserve full genes in .raw
print("\n4.2 Preserving full genes in .raw...")
adata.raw = sc.AnnData(
    X=adata.layers['counts'],
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"âœ“ Saved to .raw: {adata.raw.shape}")

# Subset to HVG
adata = adata[:, adata.var['highly_variable']].copy()
print(f"âœ“ Training shape: {adata.shape}")


# ===== STEP 5: QC After Subsetting =====

print("\n" + "="*80)
print("STEP 5: QC After Subsetting")
print("="*80)

batch_sizes = adata.obs[BATCH_KEY].value_counts()
print(f"Batch sizes: min={batch_sizes.min()}, max={batch_sizes.max()}")

small_batches = batch_sizes[batch_sizes < MIN_CELLS_PER_BATCH].index
if len(small_batches) > 0:
    print(f"Removing {len(small_batches)} small batches")
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    adata.raw = adata.raw[adata.obs_names]

if SAVE_CHECKPOINTS:
    checkpoint_file = f"{CHECKPOINT_DIR}/adata_preprocessed.h5ad"
    adata.write_h5ad(checkpoint_file, compression='gzip')
    print(f"âœ“ Checkpoint: {checkpoint_file}")


# ===== Model Existence Check (v2.6 NEW) =====

print("\n" + "="*80)
print("⭐ Model Status Check (v2.6)")
print("="*80)

scvi_model_exists = check_model_exists(scvi_model_path)
scanvi_major_exists = check_model_exists(scanvi_major_model_path)
scanvi_fine_exists = check_model_exists(scanvi_fine_model_path)
hvg_file_exists = hvg_file_path.exists()

print(f"\nExisting Models:")
print(f"  scVI model:         {scvi_model_exists} ({scvi_model_path})")
print(f"  scANVI-major model: {scanvi_major_exists} ({scanvi_major_model_path})")
print(f"  scANVI-fine model:  {scanvi_fine_exists} ({scanvi_fine_model_path})")
print(f"  HVG file:           {hvg_file_exists} ({hvg_file_path})")

print(f"\nLoad Settings:")
print(f"  LOAD_SCVI_IF_EXISTS:              {LOAD_SCVI_IF_EXISTS}")
print(f"  LOAD_SCANVI_MAJOR_IF_EXISTS:      {LOAD_SCANVI_MAJOR_IF_EXISTS}")
print(f"  LOAD_SCANVI_FINE_IF_EXISTS:       {LOAD_SCANVI_FINE_IF_EXISTS}")

# Warn if trying to load but HVG file missing
if (LOAD_SCVI_IF_EXISTS and scvi_model_exists) and not hvg_file_exists:
    print("\n⚠️  WARNING: scVI model exists but HVG file missing!")
    print("   Model may be incompatible if current HVG selection differs")

print("\nPlan:")
if LOAD_SCVI_IF_EXISTS and scvi_model_exists:
    print("  → Will attempt to LOAD scVI model")
else:
    print("  → Will TRAIN new scVI model")
    
if LOAD_SCANVI_MAJOR_IF_EXISTS and scanvi_major_exists:
    print("  → Will attempt to LOAD scANVI-major model")
else:
    print("  → Will TRAIN new scANVI-major model")
    
if LOAD_SCANVI_FINE_IF_EXISTS and scanvi_fine_exists:
    print("  → Will attempt to LOAD scANVI-fine model")
else:
    print("  → Will TRAIN new scANVI-fine model")


# ===== STEP 6: scVI Model (Load or Train) =====

print("\n" + "="*80)
print("STEP 6: scVI Model (Load or Train)")
print("="*80)

# Load existing model
if scvi_model_exists and LOAD_SCVI_IF_EXISTS:
    print("\n⭐ Loading existing scVI model...")
    try:
        # Setup anndata (required for loading)
        scvi.model.SCVI.setup_anndata(
            adata,
            layer='counts',
            batch_key=BATCH_KEY,
            continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
        )
        
        scvi_model = scvi.model.SCVI.load(str(scvi_model_path), adata=adata)
        print("✓ scVI model loaded successfully")
        
        # Verify HVG match
        if hvg_file_exists:
            saved_hvg = pd.read_csv(hvg_file_path, header=None)[0].tolist()
            current_hvg = adata.var_names.tolist()
            if saved_hvg != current_hvg:
                print("\n⚠️  WARNING: HVG mismatch detected!")
                print(f"   Saved HVG: {len(saved_hvg)} genes")
                print(f"   Current HVG: {len(current_hvg)} genes")
                print("   Retraining model with current HVG...")
                scvi_model_exists = False
    
    except Exception as e:
        print(f"\n⚠️  Model loading failed: {e}")
        print("   Falling back to training new model...")
        scvi_model_exists = False

# Train new model
if not scvi_model_exists or not LOAD_SCVI_IF_EXISTS:
    print("\n⭐ Training new scVI model...")
    
    # Setup with covariates (from allcells line 902-908)
    print("\n6.1 Setting up scVI with covariates...")
    
    scvi.model.SCVI.setup_anndata(
        adata,
        layer='counts',
        batch_key=BATCH_KEY,
        continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
    )
    
    print("✓ Registered covariates:")
    print("  - pct_counts_mt (continuous)")
    print("  - stress_score (continuous)")
    print("  - S_score (continuous)")
    print("  - G2M_score (continuous)")
    
    # scArches-ready parameters (from allcells line 911-920)
    print("\n6.2 Creating scVI model (scArches-ready)...")
    
    scvi_model = scvi.model.SCVI(
        adata,
        n_latent=N_LATENT_SCVI,
        n_layers=N_LAYERS,
        n_hidden=N_HIDDEN,
        dropout_rate=DROPOUT_RATE,
        gene_likelihood=GENE_LIKELIHOOD,
        dispersion=DISPERSION,
        encode_covariates=SCVI_ENCODE_COVARIATES,
        use_layer_norm=SCVI_USE_LAYER_NORM,
        use_batch_norm=SCVI_USE_BATCH_NORM
    )
    
    print(f"✓ Model created:")
    print(f"  n_latent: {N_LATENT_SCVI}")
    print(f"  n_layers: {N_LAYERS}")
    print(f"  encode_covariates: {SCVI_ENCODE_COVARIATES}")
    print(f"  use_layer_norm: {SCVI_USE_LAYER_NORM}")
    
    print("\n6.3 Training scVI...")
    
    scvi_model.train(
        max_epochs=SCVI_MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOPPING,
        train_size=0.9,
        plan_kwargs={'lr': 1e-3}
    )
    
    print("✓ Training complete")
    
    # Save model
    scvi_model.save(str(scvi_model_path), overwrite=True)
    print(f"✓ Model saved: {scvi_model_path}")
    
    # Save HVG list
    scvi_model_path.mkdir(parents=True, exist_ok=True)
    pd.Series(adata.var_names.astype(str)).to_csv(hvg_file_path, index=False, header=False)
    print(f"✓ HVG list saved: {hvg_file_path}")

# Extract latent representation (whether loaded or trained)
adata.obsm['X_scvi'] = scvi_model.get_latent_representation()
print(f"\n✓ Latent representation: {adata.obsm['X_scvi'].shape}")

# â­ STEP 6.5: Unknown Cleaning (Strategy B+C) =====

print("\n" + "="*80)
print("â­ STEP 6.5: Unknown Cleaning (Before scANVI Training)")
print("="*80)


# Store original labels for comparison
adata.obs['scanvi_labels_major_original'] = adata.obs['scanvi_labels_major'].copy()
adata.obs['scanvi_labels_fine_original'] = adata.obs['scanvi_labels_fine'].copy()

# ⭐ CRITICAL FIX: Re-ensure Unknown is in categories (may have been removed during copy)
print("\nRe-adding Unknown to label categories after copy...")
for k in ['scanvi_labels_major', 'scanvi_labels_fine']:
    if UNLABELED_CATEGORY not in adata.obs[k].cat.categories:
        adata.obs[k] = adata.obs[k].cat.add_categories([UNLABELED_CATEGORY])
        print(f"  Added {UNLABELED_CATEGORY} to {k}")
    else:
        print(f"  {UNLABELED_CATEGORY} already in {k}")

# â­ Strategy B: Neighborhood purity (from allcells line 1100-1150)
print("\n6.5.1 Strategy B: Neighborhood label agreement purity...")

knn = NearestNeighbors(n_neighbors=30)
knn.fit(adata.obsm['X_scvi'])
indices = knn.kneighbors(return_distance=False)

# Compute purity for major labels
labels_major = adata.obs['scanvi_labels_major'].values
purity_major = []

for i, neighs in enumerate(indices):
    own_label = labels_major[i]
    if own_label == UNLABELED_CATEGORY:
        purity_major.append(1.0)
    else:
        neighbor_labels = labels_major[neighs[1:]]  # Exclude self
        purity_major.append((neighbor_labels == own_label).mean())

adata.obs['label_purity_major'] = purity_major

# Compute purity for fine labels
labels_fine = adata.obs['scanvi_labels_fine'].values
purity_fine = []

for i, neighs in enumerate(indices):
    own_label = labels_fine[i]
    if own_label == UNLABELED_CATEGORY:
        purity_fine.append(1.0)
    else:
        neighbor_labels = labels_fine[neighs[1:]]
        purity_fine.append((neighbor_labels == own_label).mean())

adata.obs['label_purity_fine'] = purity_fine

# Mark low purity as Unknown
low_purity_major_mask = adata.obs['label_purity_major'] < EXISTING_LABEL_PURITY_THRESHOLD
low_purity_fine_mask = adata.obs['label_purity_fine'] < EXISTING_LABEL_PURITY_THRESHOLD

n_low_purity_major = low_purity_major_mask.sum()
n_low_purity_fine = low_purity_fine_mask.sum()

print(f"  Major low purity (<{EXISTING_LABEL_PURITY_THRESHOLD}): {n_low_purity_major} ({100*n_low_purity_major/adata.n_obs:.1f}%)")
print(f"  Fine low purity (<{EXISTING_LABEL_PURITY_THRESHOLD}): {n_low_purity_fine} ({100*n_low_purity_fine/adata.n_obs:.1f}%)")

adata.obs.loc[low_purity_major_mask, 'scanvi_labels_major'] = UNLABELED_CATEGORY
adata.obs.loc[low_purity_fine_mask, 'scanvi_labels_fine'] = UNLABELED_CATEGORY

# â­ Strategy C: Quality gating (from allcells line 1200-1250)
print("\n6.5.2 Strategy C: Quality gating (mt% + stress)...")

# High mt%
high_mt_mask = adata.obs['pct_counts_mt'] > EXISTING_LABEL_MT_THRESHOLD
n_high_mt = high_mt_mask.sum()
print(f"  High mt% (>{EXISTING_LABEL_MT_THRESHOLD}): {n_high_mt} ({100*n_high_mt/adata.n_obs:.1f}%)")

adata.obs.loc[high_mt_mask, 'scanvi_labels_major'] = UNLABELED_CATEGORY
adata.obs.loc[high_mt_mask, 'scanvi_labels_fine'] = UNLABELED_CATEGORY

# High stress
stress_threshold = adata.obs['stress_score'].quantile(EXISTING_LABEL_STRESS_PERCENTILE/100)
high_stress_mask = adata.obs['stress_score'] > stress_threshold
n_high_stress = high_stress_mask.sum()
print(f"  High stress (>p{EXISTING_LABEL_STRESS_PERCENTILE}): {n_high_stress} ({100*n_high_stress/adata.n_obs:.1f}%)")

adata.obs.loc[high_stress_mask, 'scanvi_labels_major'] = UNLABELED_CATEGORY
adata.obs.loc[high_stress_mask, 'scanvi_labels_fine'] = UNLABELED_CATEGORY

# Summary
total_unknown_major = (adata.obs['scanvi_labels_major'] == UNLABELED_CATEGORY).sum()
total_unknown_fine = (adata.obs['scanvi_labels_fine'] == UNLABELED_CATEGORY).sum()

print(f"\nâœ“ Unknown cleaning complete:")
print(f"  Major: {total_unknown_major} Unknown ({100*total_unknown_major/adata.n_obs:.1f}%)")
print(f"  Fine: {total_unknown_fine} Unknown ({100*total_unknown_fine/adata.n_obs:.1f}%)")

print("\nFinal label distribution:")

# ===== STEP 7: scANVI-Major Model (Load or Train) =====

print("\n" + "="*80)
print("STEP 7: scANVI-Major Model (Load or Train)")
print("="*80)

print(f"Label column: scanvi_labels_major ({adata.obs['scanvi_labels_major'].nunique()} types)")

# Load existing model
if scanvi_major_exists and LOAD_SCANVI_MAJOR_IF_EXISTS:
    print("\n⭐ Loading existing scANVI-major model...")
    try:
        # Setup anndata (required)
        scvi.model.SCANVI.setup_anndata(
            adata,
            layer='counts',
            batch_key=BATCH_KEY,
            labels_key='scanvi_labels_major',
            unlabeled_category=UNLABELED_CATEGORY,
            continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
        )
        
        scanvi_major_model = scvi.model.SCANVI.load(str(scanvi_major_model_path), adata=adata)
        print("✓ scANVI-major model loaded successfully")
        
    except Exception as e:
        print(f"\n⚠️  Model loading failed: {e}")
        print("   Falling back to training new model...")
        scanvi_major_exists = False

# Train new model
if not scanvi_major_exists or not LOAD_SCANVI_MAJOR_IF_EXISTS:
    print("\n⭐ Training new scANVI-major model...")
    
    try:
        scanvi_major_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category=UNLABELED_CATEGORY,
            labels_key='scanvi_labels_major'
        )
        print("✓ Initialized from scVI")
        
    except RuntimeError as e:
        if "size mismatch" in str(e):
            print("⚠️ Using fallback initialization")
            
            scvi.model.SCANVI.setup_anndata(
                adata,
                layer='counts',
                batch_key=BATCH_KEY,
                labels_key='scanvi_labels_major',
                unlabeled_category=UNLABELED_CATEGORY,
                continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
            )
            
            scanvi_major_model = scvi.model.SCANVI(
                adata,
                n_latent=N_LATENT_SCANVI,
                unlabeled_category=UNLABELED_CATEGORY,
                labels_key='scanvi_labels_major'
            )
            print("✓ Fallback initialization")
        else:
            raise
    
    # Train with class balancing (from allcells line 228)
    print("\n7.1 Training scANVI_major (with class balancing)...")
    print(f"  n_samples_per_label: {SCANVI_N_SAMPLES_PER_LABEL}")
    
    scanvi_major_model.train(
        max_epochs=SCANVI_MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOPPING,
        train_size=0.9,
        n_samples_per_label=SCANVI_N_SAMPLES_PER_LABEL,  # ⭐ Class balancing
        plan_kwargs={'lr': 1e-3}
    )
    
    print("✓ Training complete")
    
    # Save model
    scanvi_major_model.save(str(scanvi_major_model_path), overwrite=True)
    print(f"✓ Model saved: {scanvi_major_model_path}")

# Extract predictions (whether loaded or trained)
adata.obsm['X_scanvi_major'] = scanvi_major_model.get_latent_representation()
adata.obs['scanvi_major_pred'] = scanvi_major_model.predict()

pred_probs_major = scanvi_major_model.predict(soft=True)
if isinstance(pred_probs_major, np.ndarray):
    adata.obs['scanvi_major_conf'] = pred_probs_major.max(axis=1)
else:
    adata.obs['scanvi_major_conf'] = np.asarray(pred_probs_major).max(axis=1)

print("\n✓ Predictions extracted")


# ===== STEP 8: scANVI-Fine Model (Load or Train) =====

print("\n" + "="*80)
print("STEP 8: scANVI-Fine Model (Load or Train)")
print("="*80)

print(f"Label column: scanvi_labels_fine ({adata.obs['scanvi_labels_fine'].nunique()} types)")

# Load existing model
if scanvi_fine_exists and LOAD_SCANVI_FINE_IF_EXISTS:
    print("\n⭐ Loading existing scANVI-fine model...")
    try:
        # Setup anndata (required)
        scvi.model.SCANVI.setup_anndata(
            adata,
            layer='counts',
            batch_key=BATCH_KEY,
            labels_key='scanvi_labels_fine',
            unlabeled_category=UNLABELED_CATEGORY,
            continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
        )
        
        scanvi_fine_model = scvi.model.SCANVI.load(str(scanvi_fine_model_path), adata=adata)
        print("✓ scANVI-fine model loaded successfully")
        
    except Exception as e:
        print(f"\n⚠️  Model loading failed: {e}")
        print("   Falling back to training new model...")
        scanvi_fine_exists = False

# Train new model
if not scanvi_fine_exists or not LOAD_SCANVI_FINE_IF_EXISTS:
    print("\n⭐ Training new scANVI-fine model...")
    
    try:
        scanvi_fine_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category=UNLABELED_CATEGORY,
            labels_key='scanvi_labels_fine'
        )
        print("✓ Initialized from scVI")
        
    except RuntimeError as e:
        if "size mismatch" in str(e):
            print("⚠️ Using fallback initialization")
            
            scvi.model.SCANVI.setup_anndata(
                adata,
                layer='counts',
                batch_key=BATCH_KEY,
                labels_key='scanvi_labels_fine',
                unlabeled_category=UNLABELED_CATEGORY,
                continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']
            )
            
            scanvi_fine_model = scvi.model.SCANVI(
                adata,
                n_latent=N_LATENT_SCANVI,
                unlabeled_category=UNLABELED_CATEGORY,
                labels_key='scanvi_labels_fine'
            )
            print("✓ Fallback initialization")
        else:
            raise
    
    print("\n8.1 Training scANVI_fine (with class balancing)...")
    print(f"  n_samples_per_label: {SCANVI_N_SAMPLES_PER_LABEL}")
    
    scanvi_fine_model.train(
        max_epochs=SCANVI_MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOPPING,
        train_size=0.9,
        n_samples_per_label=SCANVI_N_SAMPLES_PER_LABEL,
        plan_kwargs={'lr': 1e-3}
    )
    
    print("✓ Training complete")
    
    # Save model
    scanvi_fine_model.save(str(scanvi_fine_model_path), overwrite=True)
    print(f"✓ Model saved: {scanvi_fine_model_path}")

# Extract predictions (whether loaded or trained)
adata.obsm['X_scanvi_fine'] = scanvi_fine_model.get_latent_representation()
adata.obs['scanvi_fine_pred'] = scanvi_fine_model.predict()

pred_probs_fine = scanvi_fine_model.predict(soft=True)
if isinstance(pred_probs_fine, np.ndarray):
    adata.obs['scanvi_fine_conf'] = pred_probs_fine.max(axis=1)
else:
    adata.obs['scanvi_fine_conf'] = np.asarray(pred_probs_fine).max(axis=1)

print("\n✓ Predictions extracted")

# Cleanup
del scvi_model, scanvi_major_model, scanvi_fine_model
gc.collect()

# ===== STEP 9: UMAPs =====

print("\n" + "="*80)
print("STEP 9: Computing UMAPs")
print("="*80)

sc.pp.neighbors(adata, use_rep='X_scanvi_major', n_neighbors=15, key_added='scanvi_major')
sc.tl.umap(adata, neighbors_key='scanvi_major')
adata.obsm['X_umap_major'] = adata.obsm['X_umap'].copy()
print("âœ“ scANVI_major UMAP")

sc.pp.neighbors(adata, use_rep='X_scanvi_fine', n_neighbors=15, key_added='scanvi_fine')
sc.tl.umap(adata, neighbors_key='scanvi_fine')
adata.obsm['X_umap_fine'] = adata.obsm['X_umap'].copy()
print("âœ“ scANVI_fine UMAP")


# ===== STEP 10: Statistics & Visualization =====

print("\n" + "="*80)
print("STEP 10: Statistics & Visualization")
print("="*80)

# Fixed stats alignment
print("\n10.1 Generating statistics (with correct alignment)...")

vc_major = adata.obs['scanvi_major_pred'].value_counts()
mean_conf_major = adata.obs.groupby('scanvi_major_pred')['scanvi_major_conf'].mean()

stats_major = pd.DataFrame({
    'Cell_Type': vc_major.index,
    'Count': vc_major.values,
    'Percentage': 100 * (vc_major / vc_major.sum()).values,
    'Mean_Confidence': mean_conf_major.reindex(vc_major.index).values
})
stats_major = stats_major.sort_values('Count', ascending=False)
stats_major.to_csv(f"{OUTPUT_DIR}/scanvi_major_statistics.csv", index=False)
print("âœ“ Major stats")

vc_fine = adata.obs['scanvi_fine_pred'].value_counts()
mean_conf_fine = adata.obs.groupby('scanvi_fine_pred')['scanvi_fine_conf'].mean()

stats_fine = pd.DataFrame({
    'Cell_Type': vc_fine.index,
    'Count': vc_fine.values,
    'Percentage': 100 * (vc_fine / vc_fine.sum()).values,
    'Mean_Confidence': mean_conf_fine.reindex(vc_fine.index).values
})
stats_fine = stats_fine.sort_values('Count', ascending=False)
stats_fine.to_csv(f"{OUTPUT_DIR}/scanvi_fine_statistics.csv", index=False)
print("âœ“ Fine stats")

# Unknown cleaning summary
unknown_summary = pd.DataFrame({
    'Metric': [
        'Initial Unknown (low conf)',
        'Low purity (major)',
        'Low purity (fine)',
        'High mt%',
        'High stress',
        'Total Unknown (major)',
        'Total Unknown (fine)'
    ],
    'Count': [
        (adata.obs['scanvi_labels_major_original'] == UNLABELED_CATEGORY).sum() if 'celltypist_conf_score' in adata.obs else 0,
        n_low_purity_major,
        n_low_purity_fine,
        n_high_mt,
        n_high_stress,
        total_unknown_major,
        total_unknown_fine
    ]
})
unknown_summary.to_csv(f"{OUTPUT_DIR}/unknown_cleaning_summary.csv", index=False)
print("âœ“ Unknown summary")

# Comparison plot
print("\n10.2 Creating comparison visualizations...")

fig, axes = plt.subplots(3, 2, figsize=(16, 22))

adata.obsm['X_umap'] = adata.obsm['X_umap_major']
sc.pl.umap(adata, color='major_lineage', ax=axes[0, 0], title='Original Major Lineages', s=5, show=False)
adata.obsm['X_umap'] = adata.obsm['X_umap_fine']
sc.pl.umap(adata, color='fine_type', ax=axes[0, 1], title='Original Fine Types', legend_fontsize=7, s=5, show=False)

adata.obsm['X_umap'] = adata.obsm['X_umap_major']
sc.pl.umap(adata, color='scanvi_major_pred', ax=axes[1, 0], title='scANVI_major Predictions', s=5, show=False)
adata.obsm['X_umap'] = adata.obsm['X_umap_fine']
sc.pl.umap(adata, color='scanvi_fine_pred', ax=axes[1, 1], title='scANVI_fine Predictions', legend_fontsize=7, s=5, show=False)

adata.obsm['X_umap'] = adata.obsm['X_umap_major']
sc.pl.umap(adata, color='scanvi_major_conf', ax=axes[2, 0], title='scANVI_major Confidence', cmap='viridis', vmin=0, vmax=1, s=5, show=False)
adata.obsm['X_umap'] = adata.obsm['X_umap_fine']
sc.pl.umap(adata, color='scanvi_fine_conf', ax=axes[2, 1], title='scANVI_fine Confidence', cmap='viridis', vmin=0, vmax=1, s=5, show=False)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/dual_scanvi_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("âœ“ Comparison plot")

# Label purity distribution
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].hist(adata.obs['label_purity_major'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
axes[0].axvline(EXISTING_LABEL_PURITY_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Threshold: {EXISTING_LABEL_PURITY_THRESHOLD}')
axes[0].set_xlabel('Label Purity (Major)')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Major Label Purity Distribution')
axes[0].legend()

axes[1].hist(adata.obs['label_purity_fine'], bins=50, color='coral', alpha=0.7, edgecolor='black')
axes[1].axvline(EXISTING_LABEL_PURITY_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Threshold: {EXISTING_LABEL_PURITY_THRESHOLD}')
axes[1].set_xlabel('Label Purity (Fine)')
axes[1].set_ylabel('Number of Cells')
axes[1].set_title('Fine Label Purity Distribution')
axes[1].legend()

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/label_purity_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print("âœ“ Purity distribution plot")


# ===== STEP 11: Save Final Dataset =====

print("\n" + "="*80)
print("STEP 11: Saving Final Dataset")
print("="*80)

final_file = f"{OUTPUT_DIR}/epithelial_DUAL_SCANVI_v2_7_PRODUCTION_final.h5ad"

print("\nðŸ“Š Final structure:")
print(f"  adata.X: {adata.X.shape}")
print(f"  adata.raw.X: {adata.raw.X.shape}")
print(f"  Layers: {list(adata.layers.keys())}")
print(f"  Obsm: {list(adata.obsm.keys())}")
print(f"  Key obs columns:")
print(f"    - Covariates: pct_counts_mt, stress_score, S_score, G2M_score")
print(f"    - Purity: label_purity_major, label_purity_fine")
print(f"    - Predictions: scanvi_major_pred, scanvi_fine_pred")
print(f"    - Confidence: scanvi_major_conf, scanvi_fine_conf")

adata.obsm['X_umap'] = adata.obsm['X_umap_major']

adata.write_h5ad(final_file, compression='gzip')
print(f"\nâœ“ Saved: {final_file}")

print("\n" + "="*80)
print("âœ… PIPELINE COMPLETE")
print("="*80)
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

print(f"\nðŸŽ¯ Key outputs:")
print(f"  1. scVI model: {scvi_model_path}")
print(f"  2. scANVI_major: {scanvi_major_model_path}")
print(f"  3. scANVI_fine: {scanvi_fine_model_path}")
print(f"  4. HVG list: {hvg_file_path}")
print(f"  5. Final h5ad: {final_file}")
print(f"  6. Statistics: {OUTPUT_DIR}/*.csv")
print(f"  7. Figures: {OUTPUT_DIR}/figures/")

print(f"\nâœ… PRODUCTION FEATURES APPLIED:")
print(f"  1. P0 fixes (5 critical fixes)")
print(f"  2. Whitelist filtering: {ENABLE_WHITELIST_FILTERING}")
print(f"  3. Covariates: 4 continuous (mt%, stress, S, G2M)")
print(f"  4. Unknown cleaning:")
print(f"     - Purity-based: {n_low_purity_major} + {n_low_purity_fine}")
print(f"     - Quality-based: {n_high_mt} (mt) + {n_high_stress} (stress)")
print(f"     - Total Unknown: {total_unknown_major} (major), {total_unknown_fine} (fine)")
print(f"  5. scArches-ready: encode_cov=True, n_layers=2")
print(f"  6. Class balancing: n_samples_per_label={SCANVI_N_SAMPLES_PER_LABEL}")

print(f"\nðŸ“Š Model Quality Metrics:")
print(f"  Major predictions:")
print(f"    Mean confidence: {adata.obs['scanvi_major_conf'].mean():.3f}")
print(f"    Mean purity: {adata.obs['label_purity_major'].mean():.3f}")
print(f"  Fine predictions:")
print(f"    Mean confidence: {adata.obs['scanvi_fine_conf'].mean():.3f}")
print(f"    Mean purity: {adata.obs['label_purity_fine'].mean():.3f}")

print("\nâš ï¸  NEXT STEPS:")
print("  1. Run add_umap_operator_epithelial.py to save UMAP operators")
print("  2. Test with query mapping using scArches")
print("  3. Check epithelial mixing using guide in EPITHELIAL_SCARCHES_READY_GUIDE.md")

print("="*80)