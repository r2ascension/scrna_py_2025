#!/usr/bin/env python3
"""
stromal_vascular_scvi_celltypist_scanvi_pipeline.py

Stromal/Vascular Cells Complete Pipeline - VERSION 2.3 PRODUCTION
scVI → CellTypist → scANVI Workflow

Critical Fixes in v2.3 (Must-fix):
- ✅ FIXED: MARKER_GENES syntax error (quote mismatch)
- ✅ FIXED: pd.Categorical usage in Step 7 (use Series.astype('category'))
- ✅ FIXED: preserve_full_raw_if_missing early-return (patch symbol_base)
- ✅ FIXED: Model reuse HVG consistency (load hvg_genes.txt first)
- ✅ FIXED: Marker validation gene name mapping (symbol_base → var_names)
- ✅ IMPROVED: Rare-type statistics export to CSV

All v2.2 features retained:
- ✅ Dual UMAP visualization (scVI + scANVI)
- ✅ Rare type filtering (<10 cells → Unknown)
- ✅ neighbors_key to avoid UMAP collision
- ✅ Separate adata_model from main adata

All v2.1 fixes retained:
- ✅ Gene alignment (normalize_gene_names BEFORE raw creation)
- ✅ CPU-safe CUDA seed
- ✅ filter_genes on counts layer
- ✅ Index-aligned result writing
- ✅ Single GPU control

Author: Clinical-Bioinformatics Team
Date: 2024-12-13
Version: v2.3-PRODUCTION-STABLE
"""

import os
import sys
import warnings
from pathlib import Path
from datetime import datetime
import time
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
from scipy import sparse
from typing import Optional, Dict, List

# CellTypist
import celltypist
from celltypist import models

# Optional: mygene for gene conversion
try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    MYGENE_AVAILABLE = False
    print("⚠️  mygene not available, will use local gene names only")

# Suppress warnings
warnings.filterwarnings('ignore')

# ============================================================================
# TIME TRACKING - START
# ============================================================================
PIPELINE_START = time.time()

# ============================================================================
# REPRODUCIBILITY & GPU AUTO-DETECTION
# ============================================================================

# --- Unified Random Seeds ---
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)

# ⚠️ CRITICAL FIX: CPU-safe CUDA seed
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

scvi.settings.seed = RANDOM_SEED
scvi.settings.dl_num_workers = 0

# --- GPU Auto-Detection ---
GPU_AVAILABLE = torch.cuda.is_available()
print("\n" + "="*80)
print("GPU & REPRODUCIBILITY CHECK")
print("="*80)
print(f"Random seed: {RANDOM_SEED}")
print(f"GPU available: {GPU_AVAILABLE}")
if GPU_AVAILABLE:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")
    print(f"GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
else:
    print("⚠️  Running on CPU (training will be slow)")

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# ========== Input/Output Paths ==========
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1214/stromal_vascular_celltypist_v2_3_production"

# ========== Pretrained Model (Optional) ==========
PRETRAINED_SCVI_MODEL = '/home/h2048/data/py/1212/stromal_vascular_celltypist_v2_1_production/scvi_models/scvi_model/' # Set path to reuse model, or None to train
# Example: "/path/to/previous_run/scvi_models/scvi_model"

# ========== Cell Type Selection ==========
MAJOR_CELLTYPE_KEY = "cell_type"
TARGET_CELLTYPES = ["Endothelial", "Fibroblast", "SMC"]

# ========== CellTypist Settings ==========
CELLTYPIST_MODEL = '/home/h2048/data/source/reference/celltypist_models/Human_Lung_Atlas.pkl'
CELLTYPIST_MAJORITY_VOTING = True
CELLTYPIST_MIN_CONFIDENCE = 0.5  # Cells below this → Unknown for scANVI

# ========== Rare Type Filtering ==========
MIN_CELLS_PER_TYPE = 200  # Types with <10 cells → Unknown (for stability)

# ========== Batch Correction Settings ==========
PREFERRED_BATCH_KEYS = ["sample", "dataset", "batch"]
N_HVG = 4000
MIN_CELLS_PER_GENE = 3

# ========== scVI Model Parameters ==========
SCVI_N_LATENT = 75
SCVI_N_LAYERS = 3
SCVI_N_HIDDEN = 128
SCVI_DROPOUT_RATE = 0.1
SCVI_GENE_LIKELIHOOD = "nb"

# ========== scVI Training Parameters ==========
SCVI_MAX_EPOCHS = 800
SCVI_BATCH_SIZE = 256
SCVI_TRAIN_SIZE = 0.9
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 45

# ========== scANVI Parameters ==========
SCANVI_N_LAYERS = 3
SCANVI_N_LATENT = 75
SCANVI_DROPOUT_RATE = 0.1
SCANVI_MAX_EPOCHS = 600
SCANVI_BATCH_SIZE = 256
SCANVI_TRAIN_SIZE = 0.9
SCANVI_EARLY_STOPPING = True
SCANVI_UNLABELED_CATEGORY = "Unknown"

# ========== Model Reuse Flags ==========
USE_EXISTING_SCVI_MODEL = True
USE_EXISTING_SCANVI_MODEL = True
FORCE_RETRAIN = False

# ========== Other Settings ==========
N_NEIGHBORS = 30
DPI = 300
FIGURE_FORMAT = "pdf"
UMAP_SIZE = 2
N_JOBS = 48

print(f"\n⚙️  Configuration:")
print(f"   Pretrained scVI: {PRETRAINED_SCVI_MODEL or 'None (train fresh)'}")
print(f"   CellTypist model: {Path(CELLTYPIST_MODEL).name}")
print(f"   Rare type threshold: {MIN_CELLS_PER_TYPE} cells")
print(f"   Model reuse: scVI={USE_EXISTING_SCVI_MODEL}, scANVI={USE_EXISTING_SCANVI_MODEL}")

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

print("\n" + "="*80)
print("STROMAL/VASCULAR scVI-CellTypist-scANVI PIPELINE - v2.3 PRODUCTION-STABLE")
print("="*80)
print(f"\nPipeline Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Target Cell Types: {', '.join(TARGET_CELLTYPES)}")

# Create output directory structure
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

model_dir = output_dir / "scvi_models"
model_dir.mkdir(exist_ok=True)

scanvi_model_dir = output_dir / "scanvi_models"
scanvi_model_dir.mkdir(exist_ok=True)

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.n_jobs = N_JOBS
sc.settings.figdir = fig_dir
sc.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)

print(f"\n✓ Output directories created")

# ============================================================================
# MARKER GENE DEFINITIONS
# ============================================================================

print("\n" + "="*80)
print("MARKER GENE DEFINITIONS")
print("="*80)

MARKER_GENES = {
    "Endothelial": {
        "PanEndothelial":  ['PECAM1', 'CDH5', 'VWF', 'KDR', 'CLDN5', 'ENG'],
        "ArterialEC":      ['GJA5', 'EFNB2', 'BMX', 'HEY1', 'DLL4'],
        "VenousEC":        ['NR2F2', 'ACKR1', 'SELP', 'VCAM1'],
        "CapillaryEC":     ['CA4', 'RGCC', 'APLN', 'PLVAP'],
        "LymphaticEC":     ['PROX1', 'PDPN', 'LYVE1', 'FLT4'],
        "Proliferating":   ['MKI67', 'TOP2A']
    },
    "Fibroblast": {
        "PanFibroblast":   ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRA'],
        "PI16_Adventitial": ['PI16', 'DPT', 'MFAP5', 'CXCL14'],
        "COL15A1_Parenchymal": ['COL15A1', 'LAMA2', 'THY1'],
        "FRClike":         ['CCL19', 'CCL21', 'CXCL13', 'CD74'],
        "Myofibroblast":   ['ACTA2', 'TAGLN', 'MYL9', 'POSTN'],
        "IFNResponse":     ['ISG15', 'IFI6', 'IFIT1', 'CXCL10'],  # ⭐ FIXED: quote mismatch
        "Proliferating":   ['MKI67', 'TOP2A']
    },
    "SMC": {
        "VascularSMC":     ['ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'MYLK'],
        "ContractileSMC":  ['MYH11', 'ACTA2', 'TAGLN', 'CARMN'],
        "SyntheticSMC":    ['FN1', 'COL1A1', 'COL3A1', 'MMP2'],
        "Pericyte":        ['RGS5', 'PDGFRB', 'NOTCH3', 'MCAM'],
        "Proliferating":   ['MKI67', 'TOP2A']
    }
}

ALL_MARKERS = sorted({
    gene 
    for cell_type_markers in MARKER_GENES.values() 
    for marker_list in cell_type_markers.values() 
    for gene in marker_list
})

print(f"✓ Total unique markers: {len(ALL_MARKERS)}")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def normalize_gene_names(adata: sc.AnnData) -> None:
    """
    Add symbol_base column for CellTypist alignment.
    Priority: local columns > mygene > var_names
    
    CRITICAL: Must be called BEFORE creating .raw to ensure raw.var has symbol_base
    """
    # Make var_names unique (AnnData requirement)
    if adata.var_names.duplicated().sum() > 0:
        adata.var_names_make_unique()
    
    # Try to find symbol column
    symbol_cols = ['gene_symbol', 'gene_symbols', 'symbol', 'feature_name']
    for col in symbol_cols:
        if col in adata.var.columns:
            vals = adata.var[col].astype(str)
            non_empty = (vals.str.len() > 0).sum() / len(vals)
            if non_empty > 0.5:
                adata.var['symbol_base'] = vals.values
                print(f"✓ Using local column: '{col}'")
                return
    
    # Try mygene if ENSEMBL detected
    first_gene = str(adata.var_names[0])
    if MYGENE_AVAILABLE and first_gene.startswith('ENSG'):
        print("Attempting mygene conversion (ENSEMBL detected)...")
        try:
            mg = mygene.MyGeneInfo()
            results = mg.querymany(
                adata.var_names.tolist(),
                scopes='ensembl.gene',
                fields='symbol',
                species='human',
                as_dataframe=False
            )
            
            symbol_map = {}
            for res in results:
                if 'symbol' in res:
                    symbol_map[res['query']] = res['symbol']
            
            converted = [symbol_map.get(g, g) for g in adata.var_names]
            frac_converted = np.mean([not str(x).startswith('ENSG') for x in converted])
            
            if frac_converted > 0.3:
                adata.var['symbol_base'] = np.array(converted, dtype=object)
                print(f"✓ mygene: {int(frac_converted*len(converted))} genes converted")
                return
        except Exception as e:
            print(f"⚠️  mygene failed: {str(e)[:80]}")
    
    # Fallback: use var_names
    adata.var['symbol_base'] = adata.var_names.astype(str).values
    print("✓ Using var_names as symbol_base")


def preserve_full_raw_if_missing(adata: sc.AnnData, counts_layer: str = 'counts') -> None:
    """
    Preserve full genes to .raw if not already present.
    
    IMPORTANT: Uses shared memory reference (X=adata.layers[counts_layer])
    DO NOT modify adata.layers[counts_layer] in-place after this call.
    
    - If raw exists and has >= current genes → keep it (but patch symbol_base if missing)
    - Otherwise create raw from counts layer (shared memory, no copy)
    """
    # ⭐ CRITICAL FIX: Patch symbol_base into raw.var if missing
    if adata.raw is not None and adata.raw.n_vars >= adata.n_vars:
        if 'symbol_base' in adata.var.columns and 'symbol_base' not in adata.raw.var.columns:
            # Only patch genes that exist in both
            common_genes = adata.raw.var_names.intersection(adata.var_names)
            if len(common_genes) > 0:
                adata.raw.var.loc[common_genes, 'symbol_base'] = adata.var.loc[common_genes, 'symbol_base'].values
                print(f"✓ Patched adata.raw.var['symbol_base'] from adata.var ({len(common_genes):,} genes)")
        print(f"✓ Keep existing adata.raw: {adata.raw.n_vars:,} genes (full)")
        return
    
    # Create raw from counts (shared memory reference)
    adata.raw = sc.AnnData(
        X=adata.layers[counts_layer],  # ⚠️ Shared reference, no .copy()
        obs=adata.obs.copy(),
        var=adata.var.copy(),  # ✓ This includes symbol_base if already added
    )
    print(f"✓ Created adata.raw: {adata.raw.n_vars:,} genes (shared memory)")
    print(f"   ⚠️  Do NOT modify layers['{counts_layer}'] in-place after this")


def load_gene_list_for_pretrained(pretrained_dir: str) -> List[str]:
    """Load gene list from pretrained scVI model directory"""
    p = Path(pretrained_dir).parent  # Go up from scvi_model to models/
    
    candidates = [
        p / "hvg_genes.txt",
        p / "genes_used.txt",
        Path(pretrained_dir) / "genes.txt",
    ]
    
    for fp in candidates:
        if fp.exists():
            genes = pd.read_csv(fp, header=None)[0].astype(str).tolist()
            print(f"✓ Loaded pretrained gene list: {fp}")
            print(f"  Genes: {len(genes):,}")
            return genes
    
    raise FileNotFoundError(
        f"Cannot find gene list for pretrained model: {pretrained_dir}\n"
        f"Expected: hvg_genes.txt or genes_used.txt in {p}\n"
        f"Please save gene list when training."
    )


def build_celltypist_input(
    adata: sc.AnnData,
    model_features: Optional[pd.Index] = None,
    counts_layer: str = 'counts'
) -> sc.AnnData:
    """
    Build minimal AnnData for CellTypist from .raw or counts layer.
    Restrict to model features if provided.
    
    CRITICAL: Tries raw.var['symbol_base'] first, fallback to adata.var['symbol_base']
    """
    # Use .raw if available (full genes), else counts layer
    if adata.raw is not None:
        X = adata.raw.X
        var = adata.raw.var.copy()
    else:
        X = adata.layers[counts_layer]
        var = adata.var.copy()
    
    # Get base symbols - try raw.var first, fallback to main adata.var
    if 'symbol_base' in var.columns:
        base_symbols = var['symbol_base'].astype(str).values
    elif 'symbol_base' in adata.var.columns:
        # Fallback: use main adata.var mapping
        print("⚠️  raw.var missing symbol_base, using main adata.var")
        if adata.raw is not None:
            # Map raw.var_names to adata.var symbol_base
            mapping = adata.var['symbol_base'].to_dict()
            base_symbols = np.array([mapping.get(g, str(g)) for g in var.index])
        else:
            base_symbols = adata.var['symbol_base'].astype(str).values
    else:
        base_symbols = var.index.astype(str).values
    
    base_index = pd.Index(base_symbols)
    
    # Restrict to model features if provided
    if model_features is not None:
        in_model = base_index.isin(model_features)
        not_dup = ~base_index.duplicated(keep='first')
        keep = in_model & not_dup
    else:
        keep = ~base_index.duplicated(keep='first')
    
    # Build minimal AnnData
    X_ct = X[:, keep].copy()
    var_ct = pd.DataFrame(index=base_index[keep])
    obs_ct = pd.DataFrame(index=adata.obs_names)
    
    ad_ct = sc.AnnData(X=X_ct, obs=obs_ct, var=var_ct)
    
    # Normalize
    sc.pp.normalize_total(ad_ct, target_sum=1e4)
    sc.pp.log1p(ad_ct)
    
    return ad_ct


def extract_celltypist_confidence(pred) -> np.ndarray:
    """Extract confidence scores from CellTypist predictions"""
    if hasattr(pred, 'probability_matrix') and pred.probability_matrix is not None:
        try:
            return np.asarray(pred.probability_matrix.max(axis=1)).ravel()
        except Exception:
            pass
    
    df = getattr(pred, 'predicted_labels', None)
    if isinstance(df, pd.DataFrame):
        for col in ('conf_score', 'confidence', 'confidence_score', 'prob'):
            if col in df.columns:
                return df[col].to_numpy()
    
    return np.ones(shape=(len(pred.predicted_labels),), dtype=float)


def merge_rare_types(s: pd.Series, min_cells: int = 10, other: str = "Unknown") -> pd.Series:
    """
    Merge rare cell types (< min_cells) into 'other' category.
    
    This improves model stability and visualization clarity.
    """
    s = s.astype(str).copy()
    vc = s.value_counts()
    rare = vc[vc < min_cells].index
    
    if len(rare) > 0:
        print(f"   Merging {len(rare)} rare types (<{min_cells} cells) → {other}")
        for rt in rare[:5]:  # Show first 5
            print(f"      {rt}: {vc[rt]} cells")
        if len(rare) > 5:
            print(f"      ... and {len(rare)-5} more")
        s.loc[s.isin(rare)] = other
    
    return s


def save_celltype_counts(counts_dict: Dict[str, pd.Series], output_dir: Path) -> None:
    """Save cell type count statistics to CSV files"""
    for name, series in counts_dict.items():
        df = pd.DataFrame({
            'cell_type': series.index,
            'count': series.values,
            'percentage': (series.values / series.sum() * 100).round(2)
        })
        df = df.sort_values('count', ascending=False)
        filepath = output_dir / f"{name}.csv"
        df.to_csv(filepath, index=False)
        print(f"   ✓ Saved: {filepath.name}")


# ============================================================================
# STEP 1: DATA LOADING & SUBSET
# ============================================================================

print("\n" + "="*80)
print("STEP 1: DATA LOADING & SUBSETTING")
print("="*80)

print(f"\nLoading main dataset: {INPUT_H5AD}")
if not Path(INPUT_H5AD).exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_H5AD}")

adata_full = sc.read_h5ad(INPUT_H5AD)
print(f"✓ Full dataset loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars:,} genes")

# Auto-detect batch key
print(f"\n⭐ Auto-detecting batch key...")
BATCH_KEY = None
for key in PREFERRED_BATCH_KEYS:
    if key in adata_full.obs.columns:
        BATCH_KEY = key
        print(f"✓ Found batch key: '{BATCH_KEY}'")
        break

if BATCH_KEY is None:
    raise ValueError(f"No batch key found. Tried: {PREFERRED_BATCH_KEYS}")

# Check cell type distribution
if MAJOR_CELLTYPE_KEY not in adata_full.obs.columns:
    raise ValueError(f"Cell type key '{MAJOR_CELLTYPE_KEY}' not found")

print(f"\nCell type distribution:")
print(adata_full.obs[MAJOR_CELLTYPE_KEY].value_counts())

# Subset to target cell types
print(f"\nSubsetting to: {TARGET_CELLTYPES}")
mask = adata_full.obs[MAJOR_CELLTYPE_KEY].isin(TARGET_CELLTYPES)
adata = adata_full[mask].copy()

print(f"✓ Subset: {adata.n_obs:,} cells ({adata.n_obs/adata_full.n_obs*100:.1f}%)")

# Clean up
del adata_full
gc.collect()

# ============================================================================
# STEP 2: DATA PREPROCESSING - CRITICAL ORDER
# ============================================================================

print("\n" + "="*80)
print("STEP 2: DATA PREPROCESSING (CORRECTED ORDER)")
print("="*80)

# ⭐ CRITICAL FIX 1: Ensure counts layer exists FIRST
if 'counts' not in adata.layers:
    if adata.raw is not None:
        adata.layers['counts'] = adata.raw.X
        print("✓ layers['counts'] <- .raw.X")
    else:
        if sparse.issparse(adata.X):
            adata.layers['counts'] = adata.X.copy()
        else:
            adata.layers['counts'] = sparse.csr_matrix(adata.X)
        print("✓ layers['counts'] <- .X (assumed counts)")
else:
    print("✓ layers['counts'] already exists")

# ⭐ CRITICAL FIX 2: Filter genes on counts layer
print(f"\nFiltering genes (min_cells={MIN_CELLS_PER_GENE})...")
n_genes_before = adata.n_vars

# Ensure we filter based on counts
counts_for_filter = adata.layers['counts']
gene_counts = np.asarray((counts_for_filter > 0).sum(axis=0)).ravel()
keep_genes = gene_counts >= MIN_CELLS_PER_GENE

adata = adata[:, keep_genes].copy()
print(f"✓ Genes: {n_genes_before:,} → {adata.n_vars:,}")

# ⭐⭐⭐ CRITICAL FIX 3: Normalize gene names BEFORE creating .raw ⭐⭐⭐
print(f"\n⭐⭐⭐ CRITICAL: Normalizing gene names BEFORE .raw creation...")
normalize_gene_names(adata)
print(f"✓ Gene names normalized (symbol_base added to adata.var)")

# ⭐ CRITICAL FIX 4: Preserve full genes to .raw (now includes symbol_base)
print(f"\n⭐ Preserving full genes to .raw (shared memory)...")
preserve_full_raw_if_missing(adata, counts_layer='counts')
print(f"   Full genes available: {adata.raw.n_vars:,}")
print(f"   ✓ raw.var includes symbol_base: {'symbol_base' in adata.raw.var.columns}")

# ============================================================================
# STEP 3: HVG SELECTION WITH MODEL REUSE LOGIC
# ============================================================================

print("\n" + "="*80)
print("STEP 3: HVG SELECTION (WITH MODEL REUSE CONSISTENCY)")
print("="*80)

hvg_file = model_dir / "hvg_genes.txt"
scvi_model_path = model_dir / "scvi_model"

# ⭐ CRITICAL FIX: Model reuse HVG consistency
if (USE_EXISTING_SCVI_MODEL and 
    scvi_model_path.exists() and 
    hvg_file.exists() and 
    not FORCE_RETRAIN):
    
    print(f"\n✅ REUSE MODE: Loading HVG list from hvg_genes.txt for model consistency")
    genes_for_training = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
    
    # Create HVG mask
    if 'symbol_base' in adata.var.columns:
        gene_col = adata.var['symbol_base'].astype(str)
    else:
        gene_col = adata.var_names.astype(str)
    
    genes_in_data = [g for g in genes_for_training if g in gene_col.values]
    missing = len(genes_for_training) - len(genes_in_data)
    
    if missing > 0:
        print(f"⚠️  {missing}/{len(genes_for_training)} genes not in current data")
        if missing / len(genes_for_training) > 0.2:
            print(f"⚠️  Too many missing genes (>{20}%), will re-select HVG")
            use_reuse_hvg = False
        else:
            use_reuse_hvg = True
    else:
        use_reuse_hvg = True
    
    if use_reuse_hvg:
        hvg_mask = gene_col.isin(genes_in_data).to_numpy()
        hvg_method = "reuse_hvg_file"
        print(f"✓ Using HVG from file: {len(genes_in_data):,} genes")
    else:
        use_reuse_hvg = False
        
else:
    use_reuse_hvg = False

# If not reusing, select HVG
if not use_reuse_hvg:
    if PRETRAINED_SCVI_MODEL and Path(PRETRAINED_SCVI_MODEL).exists():
        # ⭐ PRETRAINED MODE: Load gene list from model
        print(f"\n⭐ PRETRAINED MODE")
        genes_for_training = load_gene_list_for_pretrained(PRETRAINED_SCVI_MODEL)
        
        # Create HVG mask
        if 'symbol_base' in adata.var.columns:
            gene_col = adata.var['symbol_base'].astype(str)
        else:
            gene_col = adata.var_names.astype(str)
        
        genes_in_data = [g for g in genes_for_training if g in gene_col.values]
        missing = len(genes_for_training) - len(genes_in_data)
        
        if missing > 0:
            print(f"⚠️  {missing}/{len(genes_for_training)} genes not in data")
            if missing / len(genes_for_training) > 0.2:
                raise ValueError("Too many missing genes for pretrained model")
        
        # Create mask
        hvg_mask = gene_col.isin(genes_in_data).to_numpy()
        hvg_method = "pretrained_gene_list"
        
    else:
        # ⭐ TRAINING MODE: Select HVG (but don't subset yet)
        print(f"\n⭐ TRAINING MODE: Selecting {N_HVG} HVGs")
        
        try:
            print(f"   Attempting batch-aware HVG...")
            sc.pp.highly_variable_genes(
                adata,
                layer='counts',
                n_top_genes=N_HVG,
                batch_key=BATCH_KEY,
                subset=False,  # ⚠️ NEVER subset main adata
                flavor='seurat_v3'
            )
            hvg_method = "batch-aware"
            print(f"   ✓ Success")
        except Exception as e:
            print(f"   ⚠️  Failed: {str(e)[:80]}")
            print(f"   Falling back to non-batch-aware...")
            sc.pp.highly_variable_genes(
                adata,
                layer='counts',
                n_top_genes=N_HVG,
                subset=False,  # ⚠️ NEVER subset
                flavor='seurat_v3'
            )
            hvg_method = "non-batch-aware"
        
        hvg_mask = adata.var['highly_variable'].to_numpy()
        n_hvg = int(hvg_mask.sum())
        print(f"✓ Selected {n_hvg:,} HVGs ({hvg_method})")
        
        # Save gene list for future reuse
        if 'symbol_base' in adata.var.columns:
            genes_to_save = adata.var.loc[hvg_mask, 'symbol_base'].astype(str)
        else:
            genes_to_save = adata.var_names[hvg_mask].astype(str)
        
        pd.Series(genes_to_save.values).to_csv(hvg_file, index=False, header=False)
        print(f"✓ Gene list saved: {hvg_file}")

print(f"\n✓ Main adata preserved: {adata.n_vars:,} genes (full)")
print(f"✓ HVG mask created: {int(hvg_mask.sum()):,} genes for training")

# ============================================================================
# STEP 4: CREATE SEPARATE adata_model FOR scVI/scANVI
# ============================================================================

print("\n" + "="*80)
print("STEP 4: CREATE SEPARATE adata_model (HVG ONLY)")
print("="*80)

print(f"\n⭐⭐⭐ Creating dedicated adata_model for scVI/scANVI...")
print(f"   This preserves main adata (full genes) for CellTypist & visualization")

# Extract HVG counts
X_hvg = adata.layers['counts'][:, hvg_mask].copy()
obs_hvg = adata.obs[[BATCH_KEY]].copy()
var_hvg = adata.var.loc[hvg_mask].copy()

adata_model = sc.AnnData(X=X_hvg, obs=obs_hvg, var=var_hvg)
adata_model.layers['counts'] = adata_model.X.copy()

# Normalize for training (log1p only, no scaling)
print(f"   Normalizing adata_model...")
sc.pp.normalize_total(adata_model, target_sum=1e4, inplace=True)
sc.pp.log1p(adata_model)

print(f"\n✓ adata_model created:")
print(f"   Shape: {adata_model.n_obs:,} cells × {adata_model.n_vars:,} genes")
print(f"   .X: log1p normalized")
print(f"   layers['counts']: raw counts (for scVI)")

# Save metadata
adata.uns['hvg_method'] = hvg_method
adata.uns['n_hvg'] = int(hvg_mask.sum())

# ============================================================================
# STEP 5: scVI MODEL TRAINING/LOADING ON adata_model
# ============================================================================

print("\n" + "="*80)
print("STEP 5: scVI BATCH CORRECTION (on adata_model)")
print("="*80)

trained_scvi_this_run = False

# Setup AnnData for scVI
print(f"\nSetting up adata_model for scVI...")
print(f"   Using layer='counts' (HVG only)")
print(f"   Batch key: {BATCH_KEY}")

scvi.model.SCVI.setup_anndata(
    adata_model,
    layer='counts',
    batch_key=BATCH_KEY
)
print(f"✓ scVI data setup complete")

# Load or train scVI
if PRETRAINED_SCVI_MODEL and Path(PRETRAINED_SCVI_MODEL).exists():
    print(f"\n⭐ Loading pretrained scVI: {PRETRAINED_SCVI_MODEL}")
    scvi_model = scvi.model.SCVI.load(str(PRETRAINED_SCVI_MODEL), adata=adata_model)
    print(f"✓ Model loaded")
    
elif USE_EXISTING_SCVI_MODEL and not FORCE_RETRAIN and scvi_model_path.exists():
    model_files = list(scvi_model_path.glob("*"))
    if len(model_files) > 0:
        print(f"\n✅ Loading existing scVI: {scvi_model_path}")
        try:
            scvi_model = scvi.model.SCVI.load(str(scvi_model_path), adata=adata_model)
            print(f"✓ Model loaded")
        except Exception as e:
            print(f"⚠️  Load failed: {e}")
            scvi_model = None
    else:
        scvi_model = None
else:
    scvi_model = None

# Train if needed
if scvi_model is None:
    print(f"\n⭐ Training new scVI model...")
    scvi_model = scvi.model.SCVI(
        adata_model,
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        n_hidden=SCVI_N_HIDDEN,
        dropout_rate=SCVI_DROPOUT_RATE,
        gene_likelihood=SCVI_GENE_LIKELIHOOD
    )
    
    # Robust parameter counting
    try:
        n_params = scvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)
    print(f"   Parameters: {n_params:,}")
    
    train_kwargs = {
        'max_epochs': SCVI_MAX_EPOCHS,
        'batch_size': SCVI_BATCH_SIZE,
        'train_size': SCVI_TRAIN_SIZE,
        'early_stopping': SCVI_EARLY_STOPPING,
        'early_stopping_patience': SCVI_EARLY_STOPPING_PATIENCE,
    }
    
    # ⚠️ FIXED: Explicit single GPU control
    if GPU_AVAILABLE:
        train_kwargs.update({'accelerator': 'gpu', 'devices': 1})  # Single GPU
    else:
        train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})
    
    train_start = time.time()
    scvi_model.train(**train_kwargs)
    trained_scvi_this_run = True
    train_elapsed = time.time() - train_start
    
    print(f"\n✓ Training complete: {train_elapsed/60:.1f} min")
    
    # Save model
    scvi_model.save(str(scvi_model_path), overwrite=True)
    print(f"✓ Model saved: {scvi_model_path}")

# ⭐ CRITICAL: Write latent back to MAIN adata (index-aligned)
print(f"\n⭐ Writing X_scvi to main adata (index-aligned)...")
latent_scvi = scvi_model.get_latent_representation(adata_model)
adata.obsm['X_scvi'] = latent_scvi  # Index-aligned by construction
print(f"✓ X_scvi: {adata.obsm['X_scvi'].shape}")

# ============================================================================
# STEP 6: CELLTYPIST ANNOTATION ON MAIN ADATA
# ============================================================================

print("\n" + "="*80)
print("STEP 6: CELLTYPIST ANNOTATION (ON FULL GENES)")
print("="*80)

if not Path(CELLTYPIST_MODEL).exists():
    raise FileNotFoundError(f"CellTypist model not found: {CELLTYPIST_MODEL}")

print(f"\nLoading CellTypist model: {Path(CELLTYPIST_MODEL).name}")
ct_model = models.Model.load(CELLTYPIST_MODEL)

# Get model features if available
model_features = None
for attr in ('features', 'genes', 'feature_names'):
    if hasattr(ct_model, attr):
        feats = getattr(ct_model, attr)
        try:
            model_features = pd.Index([str(x) for x in feats])
            print(f"✓ Model features: {len(model_features):,}")
            break
        except Exception:
            pass

if model_features is None:
    print("⚠️  Model features not accessible")

# Build CellTypist input from MAIN adata (full genes via .raw)
print(f"\nBuilding CellTypist input from main adata.raw...")
ad_ct = build_celltypist_input(adata, model_features=model_features, counts_layer='counts')
print(f"✓ CellTypist input: {ad_ct.n_obs:,} cells × {ad_ct.n_vars:,} genes")

# Verify gene alignment
if model_features is not None:
    overlap = len(set(ad_ct.var_names) & set(model_features))
    print(f"   Gene overlap with model: {overlap}/{len(model_features)} ({overlap/len(model_features)*100:.1f}%)")

# Run CellTypist
print(f"\nRunning CellTypist (majority_voting={CELLTYPIST_MAJORITY_VOTING})...")
pred = celltypist.annotate(
    ad_ct,
    model=ct_model,
    majority_voting=CELLTYPIST_MAJORITY_VOTING
)

# ⭐ CRITICAL FIX: Index-aligned result writing
print(f"\n⭐ Writing CellTypist results (index-aligned)...")
pred_df = pred.predicted_labels.copy()
pred_df.index = ad_ct.obs_names  # Should match adata.obs_names

# Use reindex for safety
adata.obs['cell_type_celltypist_raw'] = (
    pred_df['predicted_labels'].astype(str).reindex(adata.obs_names).values
)

if CELLTYPIST_MAJORITY_VOTING and 'majority_voting' in pred_df.columns:
    adata.obs['cell_type_celltypist_majority_raw'] = (
        pred_df['majority_voting'].astype(str).reindex(adata.obs_names).values
    )

# Confidence scores
confidence = extract_celltypist_confidence(pred)
adata.obs['celltypist_confidence'] = confidence

# Cleanup
del ad_ct, pred, pred_df
gc.collect()

# Summary
ct_types = adata.obs['cell_type_celltypist_raw'].nunique()
ct_mean_conf = float(np.mean(adata.obs['celltypist_confidence']))
print(f"\n✓ CellTypist complete")
print(f"   Unique types (raw): {ct_types}")
print(f"   Mean confidence: {ct_mean_conf:.3f}")

# ⭐ NEW: Save raw counts for audit
print(f"\nSaving CellTypist raw counts...")
save_celltype_counts(
    {
        'celltypist_counts_raw': adata.obs['cell_type_celltypist_raw'].value_counts()
    },
    output_dir
)

# ============================================================================
# STEP 7: PREPARE scANVI LABELS (WITH RARE TYPE FILTERING - FIXED)
# ============================================================================

print("\n" + "="*80)
print("STEP 7: PREPARE scANVI LABELS (WITH RARE TYPE FILTERING)")
print("="*80)

print(f"\n⭐ Creating scANVI training labels from CellTypist...")
print(f"   Low confidence threshold: {CELLTYPIST_MIN_CONFIDENCE}")
print(f"   Rare type threshold: {MIN_CELLS_PER_TYPE} cells")

# Start with raw CellTypist predictions
ct_raw = adata.obs['cell_type_celltypist_raw'].astype(str)

# Step 1: Low confidence → Unknown
print(f"\nStep 1: Handling low confidence cells...")
ct_lc = ct_raw.where(
    adata.obs['celltypist_confidence'] >= CELLTYPIST_MIN_CONFIDENCE,
    SCANVI_UNLABELED_CATEGORY
)
n_low_conf = (ct_lc == SCANVI_UNLABELED_CATEGORY).sum()
print(f"   Low confidence → Unknown: {n_low_conf:,} cells")

# Step 2: Rare types → Unknown (KEY FOR STABILITY)
print(f"\nStep 2: Merging rare types...")
ct_filt = merge_rare_types(
    ct_lc,
    min_cells=MIN_CELLS_PER_TYPE,
    other=SCANVI_UNLABELED_CATEGORY
)

# Save both versions
adata.obs['cell_type_celltypist_filt'] = ct_filt

# ⭐ CRITICAL FIX: Use Series.astype('category'), not pd.Categorical
labels = ct_filt.astype('category')
if SCANVI_UNLABELED_CATEGORY not in labels.cat.categories:
    labels = labels.cat.add_categories([SCANVI_UNLABELED_CATEGORY])

adata.obs['labels_for_scanvi'] = labels
adata_model.obs['labels_for_scanvi'] = labels.values

# Statistics
n_unknown = int((labels == SCANVI_UNLABELED_CATEGORY).sum())
n_labeled = int(adata.n_obs - n_unknown)
n_unique_types = int(labels.nunique()) - (1 if SCANVI_UNLABELED_CATEGORY in labels.cat.categories else 0)

print(f"\n✓ scANVI labels prepared:")
print(f"   Labeled cells: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"   Unknown cells: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")
print(f"   Unique types (excl. Unknown): {n_unique_types}")

# ⭐ NEW: Save filtered counts for audit
print(f"\nSaving CellTypist filtered counts...")
save_celltype_counts(
    {
        'celltypist_counts_filt': adata.obs['cell_type_celltypist_filt'].value_counts()
    },
    output_dir
)

# ============================================================================
# STEP 8: scANVI MODEL TRAINING (WITH FILTERED LABELS)
# ============================================================================

print("\n" + "="*80)
print("STEP 8: scANVI TRAINING (on adata_model with filtered labels)")
print("="*80)

scanvi_model_path = scanvi_model_dir / "scanvi_model"
trained_scanvi_this_run = False

print(f"\nscANVI configuration:")
print(f"   Base: scVI (pre-trained)")
print(f"   Labels: labels_for_scanvi (filtered CellTypist)")
print(f"   Unknown: {SCANVI_UNLABELED_CATEGORY}")

# Load or train scANVI
if USE_EXISTING_SCANVI_MODEL and not FORCE_RETRAIN and scanvi_model_path.exists():
    model_files = list(scanvi_model_path.glob("*"))
    if len(model_files) > 0:
        print(f"\n✅ Loading existing scANVI: {scanvi_model_path}")
        try:
            scanvi_model = scvi.model.SCANVI.load(str(scanvi_model_path), adata=adata_model)
            print(f"✓ Model loaded")
        except Exception as e:
            print(f"⚠️  Load failed: {e}")
            scanvi_model = None
    else:
        scanvi_model = None
else:
    scanvi_model = None

if scanvi_model is None:
    print(f"\n⭐ Initializing scANVI from scVI...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key='labels_for_scanvi',
        unlabeled_category=SCANVI_UNLABELED_CATEGORY
    )
    
    # Robust parameter counting
    try:
        n_params = scanvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scanvi_model.module.parameters() if p.requires_grad)
    print(f"   Parameters: {n_params:,}")
    
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'batch_size': SCANVI_BATCH_SIZE,
        'train_size': SCANVI_TRAIN_SIZE,
    }
    
    # ⚠️ FIXED: Explicit single GPU control
    if GPU_AVAILABLE:
        train_kwargs.update({'accelerator': 'gpu', 'devices': 1})  # Single GPU
    else:
        train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})
    
    train_start = time.time()
    scanvi_model.train(**train_kwargs)
    trained_scanvi_this_run = True
    train_elapsed = time.time() - train_start
    
    print(f"\n✓ Training complete: {train_elapsed/60:.1f} min")
    
    # Save model
    scanvi_model.save(str(scanvi_model_path), overwrite=True)
    print(f"✓ Model saved: {scanvi_model_path}")

# ⭐ CRITICAL: Write predictions back to MAIN adata (index-aligned)
print(f"\n⭐ Writing scANVI results to main adata (index-aligned)...")
adata.obs['cell_type_scanvi_raw'] = scanvi_model.predict(adata_model).astype(str)
predictions_probs = scanvi_model.predict(adata_model, soft=True)
adata.obs['scanvi_confidence'] = predictions_probs.max(axis=1)
adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation(adata_model)

# Filter rare types in scANVI predictions
print(f"\n⭐ Filtering rare types in scANVI predictions...")
adata.obs['cell_type_scanvi_filt'] = merge_rare_types(
    adata.obs['cell_type_scanvi_raw'],
    min_cells=MIN_CELLS_PER_TYPE,
    other=SCANVI_UNLABELED_CATEGORY
).astype('category')

print(f"✓ scANVI predictions complete")
print(f"   X_scanvi: {adata.obsm['X_scanvi'].shape}")
print(f"   Unique types (raw): {adata.obs['cell_type_scanvi_raw'].nunique()}")
print(f"   Unique types (filtered): {adata.obs['cell_type_scanvi_filt'].nunique()}")

# ⭐ NEW: Save scANVI counts for audit
print(f"\nSaving scANVI counts...")
save_celltype_counts(
    {
        'scanvi_counts_raw': adata.obs['cell_type_scanvi_raw'].value_counts(),
        'scanvi_counts_filt': adata.obs['cell_type_scanvi_filt'].value_counts()
    },
    output_dir
)

# Cleanup
del adata_model
gc.collect()

# ============================================================================
# STEP 9: VISUALIZATION WITH BOTH scVI AND scANVI UMAPs
# ============================================================================

print("\n" + "="*80)
print("STEP 9: VISUALIZATION (BOTH scVI AND scANVI UMAPs)")
print("="*80)

# Normalize main adata for visualization
print(f"\nNormalizing main adata for visualization...")
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)
print(f"✓ adata.X now contains log1p normalized data")

# ⭐⭐⭐ CRITICAL: Compute TWO separate UMAPs with neighbors_key ⭐⭐⭐
print(f"\n⭐⭐⭐ Computing scVI UMAP (with neighbors_key)...")
sc.pp.neighbors(
    adata,
    use_rep='X_scvi',
    n_neighbors=N_NEIGHBORS,
    key_added='neighbors_scvi'
)
sc.tl.umap(adata, neighbors_key='neighbors_scvi')
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
print(f"✓ UMAP(scVI) stored in X_umap_scvi: {adata.obsm['X_umap_scvi'].shape}")

print(f"\n⭐⭐⭐ Computing scANVI UMAP (with neighbors_key)...")
sc.pp.neighbors(
    adata,
    use_rep='X_scanvi',
    n_neighbors=N_NEIGHBORS,
    key_added='neighbors_scanvi'
)
sc.tl.umap(adata, neighbors_key='neighbors_scanvi')
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
print(f"✓ UMAP(scANVI) stored in X_umap_scanvi: {adata.obsm['X_umap_scanvi'].shape}")

# ========== Figure A: scVI space (includes CellTypist) ==========
print(f"\nGenerating scVI comparison plot...")
fig, axes = plt.subplots(1, 4, figsize=(22, 6))

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=MAJOR_CELLTYPE_KEY,
    ax=axes[0],
    show=False,
    title='Original Labels (scVI)',
    size=UMAP_SIZE,
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color='cell_type_celltypist_filt',
    ax=axes[1],
    show=False,
    title='CellTypist (filtered, scVI)',
    size=UMAP_SIZE,
    legend_loc='right margin',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color='celltypist_confidence',
    ax=axes[2],
    show=False,
    title='CellTypist Confidence (scVI)',
    size=UMAP_SIZE,
    cmap='viridis',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=BATCH_KEY,
    ax=axes[3],
    show=False,
    title='Batch (scVI)',
    size=1,
    frameon=False
)

plt.tight_layout()
plt.savefig(fig_dir / f'comparison_scvi.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close()
print(f"✓ Saved: comparison_scvi.{FIGURE_FORMAT}")

# ========== Figure B: scANVI space (final) ==========
print(f"\nGenerating scANVI comparison plot...")
fig, axes = plt.subplots(1, 4, figsize=(22, 6))

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='cell_type_scanvi_filt',
    ax=axes[0],
    show=False,
    title='scANVI Final (filtered)',
    size=UMAP_SIZE,
    legend_loc='right margin',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='scanvi_confidence',
    ax=axes[1],
    show=False,
    title='scANVI Confidence',
    size=UMAP_SIZE,
    cmap='viridis',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color=BATCH_KEY,
    ax=axes[2],
    show=False,
    title='Batch (scANVI)',
    size=1,
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='cell_type_celltypist_filt',
    ax=axes[3],
    show=False,
    title='CellTypist (filtered, scANVI space)',
    size=UMAP_SIZE,
    legend_loc='right margin',
    frameon=False
)

plt.tight_layout()
plt.savefig(fig_dir / f'comparison_scanvi.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close()
print(f"✓ Saved: comparison_scanvi.{FIGURE_FORMAT}")

# ⭐ IMPROVED: Marker validation with proper gene name mapping
print(f"\nGenerating marker validation plot...")

# Build symbol → var_names mapping
if adata.raw is not None and 'symbol_base' in adata.raw.var.columns:
    sym = adata.raw.var['symbol_base'].astype(str)
    sym_to_var = pd.Series(adata.raw.var_names.values, index=sym).drop_duplicates()
    markers = [m for m in ALL_MARKERS[:20] if m in sym_to_var.index]
    markers_var = [sym_to_var[m] for m in markers]
    use_gene_symbols = 'symbol_base'
else:
    markers_var = [m for m in ALL_MARKERS[:20] if m in adata.raw.var_names]
    use_gene_symbols = None

if len(markers_var) > 0:
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.dotplot(
        adata,
        markers_var,
        groupby='cell_type_scanvi_filt',
        use_raw=True,
        gene_symbols=use_gene_symbols,
        ax=ax,
        show=False,
        title='Top Markers (from .raw)'
    )
    plt.tight_layout()
    plt.savefig(fig_dir / f'marker_validation.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"✓ Marker validation plot saved ({len(markers_var)} markers)")
else:
    print(f"⚠️  No markers available in .raw for validation plot")

# ============================================================================
# STEP 10: SAVE FINAL RESULTS
# ============================================================================

print("\n" + "="*80)
print("STEP 10: SAVE FINAL RESULTS")
print("="*80)

checkpoint_file = output_dir / "adata_stromal_vascular_FINAL.h5ad"

# Add metadata
adata.uns['pipeline_info'] = {
    'version': '2.3-PRODUCTION-STABLE',
    'timestamp': datetime.now().isoformat(),
    'input_file': str(INPUT_H5AD),
    'target_celltypes': TARGET_CELLTYPES,
    'batch_key': BATCH_KEY,
    'n_hvg': adata.uns['n_hvg'],
    'hvg_method': hvg_method,
    'celltypist_model': str(CELLTYPIST_MODEL),
    'celltypist_confidence_threshold': CELLTYPIST_MIN_CONFIDENCE,
    'min_cells_per_type': MIN_CELLS_PER_TYPE,
    'scanvi_label_source': 'celltypist_filtered',
    'scvi_trained_this_run': trained_scvi_this_run,
    'scanvi_trained_this_run': trained_scanvi_this_run,
    'critical_fixes_v2_3': [
        'MARKER_GENES syntax error fixed (quote mismatch)',
        'pd.Categorical usage fixed (use Series.astype)',
        'preserve_full_raw_if_missing early-return fixed (patch symbol_base)',
        'Model reuse HVG consistency (load hvg_genes.txt first)',
        'Marker validation gene mapping fixed (symbol_base → var_names)',
        'Rare-type statistics export to CSV'
    ],
    'all_v2_2_features': [
        'Dual UMAP visualization (scVI + scANVI)',
        'Rare type filtering (<10 cells → Unknown)',
        'neighbors_key to avoid UMAP collision',
        'Separate adata_model from main adata'
    ],
    'all_v2_1_fixes': [
        'Gene name normalization BEFORE .raw creation',
        'CPU-safe CUDA seed setup',
        'filter_genes on counts layer',
        'Index-aligned result writing',
        'Single GPU device control'
    ]
}

adata.write_h5ad(checkpoint_file, compression='gzip')
size_gb = checkpoint_file.stat().st_size / 1e9
print(f"\n✓ Saved: {checkpoint_file}")
print(f"   Size: {size_gb:.2f} GB")

# Export annotations
annotations = adata.obs[[
    MAJOR_CELLTYPE_KEY,
    'cell_type_celltypist_raw',
    'cell_type_celltypist_filt',
    'celltypist_confidence',
    'cell_type_scanvi_raw',
    'cell_type_scanvi_filt',
    'scanvi_confidence',
    BATCH_KEY
]].copy()
annotations.to_csv(output_dir / "annotations_complete.csv")
print(f"✓ Annotations exported")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

total_time = time.time() - PIPELINE_START

print(f"\n{'='*80}")
print("🎉 PIPELINE COMPLETE - v2.3 PRODUCTION-STABLE")
print("="*80)

print(f"\n📊 Analysis Summary:")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (full): {adata.n_vars:,}")
print(f"   Genes (raw): {adata.raw.n_vars:,}")
print(f"   Genes (HVG for training): {adata.uns['n_hvg']}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

print(f"\n⭐ CellTypist Results:")
print(f"   Unique types (raw): {adata.obs['cell_type_celltypist_raw'].nunique()}")
print(f"   Unique types (filtered): {adata.obs['cell_type_celltypist_filt'].nunique()}")
print(f"   Mean confidence: {float(np.mean(adata.obs['celltypist_confidence'])):.3f}")

print(f"\n⭐ scANVI Results:")
print(f"   Unique types (raw): {adata.obs['cell_type_scanvi_raw'].nunique()}")
print(f"   Unique types (filtered): {adata.obs['cell_type_scanvi_filt'].nunique()}")
print(f"   Mean confidence: {float(np.mean(adata.obs['scanvi_confidence'])):.3f}")

# Agreement stats (using filtered versions)
ct_series = adata.obs['cell_type_celltypist_filt'].astype('object')
scanvi_series = adata.obs['cell_type_scanvi_filt'].astype('object')
agreement_mask = (
    (ct_series == scanvi_series)
    | (pd.isna(ct_series) & pd.isna(scanvi_series))
)
ct_to_scanvi = int(agreement_mask.sum())
agreement_pct = ct_to_scanvi / adata.n_obs * 100
print(f"   CellTypist-scANVI agreement (filtered): {agreement_pct:.1f}%")

print(f"\n🔧 Critical Fixes Applied (v2.3 - Must-fix):")
print(f"   ✅ MARKER_GENES syntax error (quote mismatch)")
print(f"   ✅ pd.Categorical usage (use Series.astype)")
print(f"   ✅ preserve_full_raw early-return (patch symbol_base)")
print(f"   ✅ Model reuse HVG consistency")
print(f"   ✅ Marker validation gene mapping")
print(f"   ✅ Rare-type statistics export")

print(f"\n📁 Output Files:")
print(f"   h5ad: {checkpoint_file.name}")
print(f"   Annotations: annotations_complete.csv")
print(f"   Cell type counts: celltypist_counts_*.csv, scanvi_counts_*.csv")
print(f"   scVI model: {scvi_model_path}")
print(f"   scANVI model: {scanvi_model_path}")
print(f"   Gene list: hvg_genes.txt")

print(f"\n📈 Figures:")
print(f"   comparison_scvi.{FIGURE_FORMAT} - scVI UMAP with CellTypist")
print(f"   comparison_scanvi.{FIGURE_FORMAT} - scANVI final results")
print(f"   marker_validation.{FIGURE_FORMAT} - Marker gene expression")

print(f"\n⏱️  Total time: {total_time/60:.1f} min")
print("="*80)