#!/usr/bin/env python3
"""
scANVI Downstream Analysis - Production Pipeline v1.3
======================================================

Critical fixes (P0/P1) + optional tissue-stratified runs:
1. ✅ Robust model loading with proper gene set matching
2. ✅ Explicit DE loop without comparison string dependency
3. ✅ Schema validation for DE results (cross-version tolerant)
4. ✅ Smart expression source resolution
5. ✅ Memory-safe denoised expression (no layers write)
6. ✅ Pseudobulk export for donor-level statistics (counts-safe)
7. ✅ Tissue comparison support
8. ✅ Optional tissue-stratified downstream runs (toggle)

Author: r2end
Date: 2025-12-17
Version: v1.3 (Production)
"""

import os
import sys
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scvi
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import scipy.sparse as sp
import gc
import json

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("/home/h2048/data/core_data")
OUTPUT_BASE = Path("/home/h2048/data/downstream_analysis_v1_3")

# Cell types configuration
CELL_TYPES = {
    'T_cells': {
        'h5ad': 'adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_model': 'models/T_models/scanvi_model',
        'naming_mode': 'standard',
    },
    'B_cells': {
        'h5ad': 'adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_model': 'models/B_models/scanvi_model',
        'naming_mode': 'standard',
    },
    'Myeloid': {
        'h5ad': 'adata_myeloid_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_model': 'models/Myeloid_models/scanvi_model',
        'naming_mode': 'standard',
    },
    'Stromal_Vascular': {
        'h5ad': 'adata_stromal_vascular_FINAL.h5ad',
        'scanvi_model': 'models/Stromal_models/scanvi_model',
        'naming_mode': 'semantic',
    }
}

# Analysis parameters
DE_PARAMS = {
    'mode': 'change',
    'delta': 0.5,
    'batch_correction': True,
    'n_samples': 5000,
    # Filtering thresholds
    'bayes_factor_threshold': 3.0,
    'lfc_threshold': 0.5,
    'min_proportion': 0.1,
    # when lfc_mean not available
    'proba_de_threshold': 0.95,
}

DENOISED_PARAMS = {
    'n_samples': 25,
    'library_size': 1e4,
    'max_elements_in_memory': 2e8,   # n_obs * n_vars upper bound
    'save_separate': True,
    'save_in_layers': False  # ⭐ Never write to layers by default
}

SUBCLUSTER_PARAMS = {
    'resolutions': [1.2, 1.4, 1.6, 1.8, 2.0],
    'default_resolution': 1.6,
    'n_neighbors': 30,
    'min_cells_for_subcluster': 500
}

PSEUDOBULK_PARAMS = {
    'min_cells_per_sample': 10,
    'aggregation': 'sum'
}

# Tissue comparison (if applicable)
TISSUE_COMPARISON_PAIRS = [
    # ('ethmoid_sinus', 'inferior_turbinate'),
    # ('nasal_cavity', 'sinus'),
    # Add your actual tissue pairs here
]

# Tissue-stratified downstream (recommended for multi-tissue atlases)
TISSUE_STRATIFY_PARAMS = {
    "enabled": True,              # toggle tissue stratification on/off
    "run_global_all": True,       # also run on full dataset
    "min_cells_per_tissue": 1000, # skip tiny tissues
    "min_samples_per_tissue": 3,  # skip tissues with too few donors/samples
    "tissue_whitelist": None,     # e.g. ["ethmoid sinus", "inferior turbinate"]
}

# Visualization
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
np.random.seed(42)

# ============================================================================
# SMALL HELPERS
# ============================================================================

def safe_name(x, max_len=180):
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t']:
        s = s.replace(ch, '_')
    return s[:max_len]


def _get_proba_de_col(df: pd.DataFrame):
    # scvi-tools versions may expose proba_de or proba_m2 (forum reports in 1.4.x)
    if 'proba_de' in df.columns:
        return 'proba_de'
    if 'proba_m2' in df.columns:
        return 'proba_m2'
    return None


def _compute_log2fc_from_means(df: pd.DataFrame, eps=1e-8):
    """
    Best-effort log2FC approximation if lfc_mean is absent.
    Preference: raw_normalized_mean1/2, else raw_mean1/2, else scale1/2.
    """
    pairs = [
        ('raw_normalized_mean1', 'raw_normalized_mean2'),
        ('raw_mean1', 'raw_mean2'),
        ('scale1', 'scale2'),
    ]
    for a, b in pairs:
        if a in df.columns and b in df.columns:
            x1 = df[a].astype(float).to_numpy()
            x2 = df[b].astype(float).to_numpy()
            return np.log2((x1 + eps) / (x2 + eps))
    return None


# ============================================================================
# UTILITY: HVG EXTRACTION FROM MODEL
# ============================================================================

def extract_var_names_from_model(model_path, model_type='SCANVI'):
    """
    Extract gene names (training var_names) from saved model.

    This extracts the gene set used during training.
    If you trained on HVG-only, this IS your HVG list.
    If you trained on full genes, this is the full gene set.
    """
    try:
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI

        # Method 1: Try registry first (lightest)
        print(f"  Attempting registry extraction...")
        reg = Model.load_registry(model_path)

        def find_in_dict(obj, target_key="var_names", path=""):
            results = []
            if isinstance(obj, dict):
                for k, v in obj.items():
                    p = f"{path}.{k}" if path else k
                    if target_key in str(k).lower():
                        results.append((p, v))
                    if isinstance(v, (dict, list)):
                        results.extend(find_in_dict(v, target_key, p))
            elif isinstance(obj, list):
                for i, item in enumerate(obj):
                    results.extend(find_in_dict(item, target_key, f"{path}[{i}]"))
            return results

        var_name_candidates = find_in_dict(reg, "var_names")

        if var_name_candidates:
            for path, var_names in var_name_candidates:
                if isinstance(var_names, (list, np.ndarray)) and len(var_names) > 0:
                    print(f"  ✓ Found var_names at: {path}")
                    print(f"  ✓ Extracted {len(var_names)} genes from registry")
                    return list(var_names)

        # Method 2: Load model without adata
        print(f"  Registry method incomplete, loading model...")
        model = Model.load(model_path, adata=False)
        vn_dict = model.get_var_names()

        var_names = list(vn_dict.values())[0]
        print(f"  ✓ Extracted {len(var_names)} genes from model.get_var_names()")

        return list(var_names)

    except Exception as e:
        print(f"  ⚠️ Failed to extract var_names: {e}")
        return None


def save_hvg_list(genes, output_path):
    """Save gene list to text file for reproducibility"""
    with open(output_path, 'w') as f:
        for gene in genes:
            f.write(f"{gene}\n")
    print(f"  ✓ Saved gene list: {output_path}")


# ============================================================================
# UTILITY: EXPRESSION SOURCE RESOLUTION
# ============================================================================

def resolve_expr_source(adata):
    """
    Smart resolution of expression matrix source.

    Priority:
    1. adata.raw (if exists and has data)
    2. adata.layers['log1p']
    3. adata.layers['counts'] (with warning)
    4. adata.X (fallback)

    Returns dict for scanpy: {use_raw, layer}
    """
    if adata.raw is not None and getattr(adata.raw, "X", None) is not None and adata.raw.n_vars > 0:
        return {'use_raw': True, 'layer': None}

    if 'log1p' in adata.layers:
        return {'use_raw': False, 'layer': 'log1p'}

    if 'counts' in adata.layers:
        print("  ⚠️ Using 'counts' layer - ensure normalization is appropriate for the chosen test")
        return {'use_raw': False, 'layer': 'counts'}

    return {'use_raw': False, 'layer': None}


def _looks_integer_like(x, n_check=2000, tol=1e-6):
    """Heuristic: check if values are close to integers (sparse/dense)."""
    try:
        if sp.issparse(x):
            vals = x.data
        else:
            vals = np.asarray(x).ravel()
        if vals.size == 0:
            return False
        vals = vals[: min(vals.size, n_check)]
        return np.all(np.abs(vals - np.round(vals)) < tol)
    except Exception:
        return False


def get_counts_matrix(adata):
    """
    Get raw counts matrix for pseudobulk/DESeq2/edgeR.

    Priority:
    1) adata.layers['counts']
    2) adata.raw.X only if it looks integer-like (avoid log1p in .raw)
    3) fallback: adata.X (warn)
    """
    if 'counts' in adata.layers:
        return adata.layers['counts'], adata.var_names

    if adata.raw is not None and getattr(adata.raw, "X", None) is not None:
        if _looks_integer_like(adata.raw.X):
            return adata.raw.X, adata.raw.var_names
        else:
            print("  ⚠️ adata.raw exists but does not look like integer counts; not using for pseudobulk")

    print("  ⚠️ No counts layer found; falling back to adata.X (may not be counts)")
    return adata.X, adata.var_names


# ============================================================================
# UTILITY: COLUMN NAME RESOLUTION WITH VALIDATION
# ============================================================================

def get_column_names(adata, naming_mode):
    """
    Get standardized column names with auto-detection and fallbacks.
    """
    cols = {}

    # Cell type columns (mode-specific)
    if naming_mode == 'standard':
        cols['celltype_col'] = 'scanvi_predictions'
        cols['confidence_col'] = 'scanvi_confidence'
    elif naming_mode == 'semantic':
        cols['celltype_col'] = 'cell_type_scanvi_filt'
        cols['confidence_col'] = 'scanvi_confidence'
    else:
        raise ValueError(f"Unknown naming_mode: {naming_mode}")

    # Batch column (required)
    batch_candidates = ['batch', 'dataset', 'study', 'platform']
    cols['batch_col'] = None
    for candidate in batch_candidates:
        if candidate in adata.obs.columns:
            cols['batch_col'] = candidate
            break

    # Sample column (required)
    sample_candidates = ['sample_id', 'Sample', 'donor_id', 'patient_id', 'subject']
    cols['sample_col'] = None
    for candidate in sample_candidates:
        if candidate in adata.obs.columns:
            cols['sample_col'] = candidate
            break

    # Disease status (optional)
    disease_candidates = ['disease_status', 'Disease', 'condition', 'group', 'phenotype', 'disease_level_1']
    cols['disease_col'] = None
    for candidate in disease_candidates:
        if candidate in adata.obs.columns:
            cols['disease_col'] = candidate
            break

    # Tissue location (optional)
    tissue_candidates = ['tissue', 'organ', 'anatomical_site', 'location', 'site', 'organ__ontology_label']
    cols['tissue_col'] = None
    for candidate in tissue_candidates:
        if candidate in adata.obs.columns:
            cols['tissue_col'] = candidate
            break

    return cols


def validate_columns_strict(adata, cols):
    """
    Strict validation with detailed error messages.
    Raises ValueError if critical columns missing.
    """
    errors = []
    warnings_list = []

    required = {
        'celltype_col': 'Cell type labels',
        'confidence_col': 'scANVI confidence scores',
        'batch_col': 'Batch/dataset identifier',
        'sample_col': 'Sample/donor identifier'
    }

    for key, desc in required.items():
        col = cols.get(key)
        if col is None:
            errors.append(f"Required column '{key}' ({desc}) not configured")
        elif col not in adata.obs.columns:
            errors.append(f"Configured column '{col}' for {desc} not found in adata.obs")
        else:
            if key == 'confidence_col':
                if not np.issubdtype(adata.obs[col].dtype, np.number):
                    errors.append(f"Column '{col}' should be numeric (confidence scores)")

    if 'X_scanvi' not in adata.obsm:
        warnings_list.append("X_scanvi not in adata.obsm - subclustering will be skipped")

    counts_available = ('counts' in adata.layers) or (adata.raw is not None)
    if not counts_available:
        warnings_list.append("No counts matrix (layers['counts'] or adata.raw) - pseudobulk/DE may be limited")

    if errors:
        print("\n" + "="*70)
        print("❌ COLUMN VALIDATION FAILED")
        print("="*70)
        for err in errors:
            print(f"  - {err}")
        print(f"\nAvailable obs columns (head): {list(adata.obs.columns)[:30]}")
        raise ValueError("Critical columns missing or invalid")

    if warnings_list:
        print("\n⚠️ Column Validation Warnings:")
        for warn in warnings_list:
            print(f"  - {warn}")

    print("✓ Column validation passed")


# ============================================================================
# UTILITY: DE SCHEMA VALIDATION (ROBUST)
# ============================================================================

def validate_and_fix_de_schema(de_df, source="scANVI", var_names=None):
    """
    Validate and standardize DE results schema (robust to scvi-tools version differences).

    - Avoids aggressive reset_index that can turn RangeIndex into 'gene' column.
    - For scANVI/scVI DE, requires bayes_factor OR proba_de/proba_m2 (depending on version).
    """
    if not isinstance(de_df, pd.DataFrame):
        raise ValueError(f"DE result is not a DataFrame (got {type(de_df)})")

    # only reset index if it looks like gene names (overlap with var_names)
    if var_names is not None and de_df.index is not None:
        if (de_df.index.dtype == object) and (de_df.index.name in [None, 'gene', 'genes']):
            idx_vals = de_df.index.astype(str)
            overlap = np.intersect1d(idx_vals[: min(len(idx_vals), 2000)], np.asarray(var_names).astype(str)).size
            if overlap >= 50:
                de_df = de_df.reset_index().rename(columns={'index': 'gene'})

    if source == "scANVI":
        if 'gene' not in de_df.columns:
            if de_df.index.dtype == object:
                de_df = de_df.reset_index().rename(columns={'index': 'gene'})
            else:
                raise ValueError("Cannot identify gene column in DE results")

        if ('bayes_factor' not in de_df.columns) and (_get_proba_de_col(de_df) is None):
            raise ValueError(
                "DE schema validation failed (scANVI). "
                "Expected 'bayes_factor' or 'proba_de/proba_m2'. "
                f"Available: {list(de_df.columns)[:30]}"
            )

        return de_df

    else:
        required = {'names', 'scores'}
        missing = required - set(de_df.columns)
        if missing:
            raise ValueError(
                f"DE schema validation failed ({source}). Missing columns: {missing}. "
                f"Available: {list(de_df.columns)[:30]}"
            )
        return de_df


# ============================================================================
# STEP 1: ROBUST MODEL LOADING
# ============================================================================

def log_step(step_num, step_name, celltype=""):
    suffix = f" - {celltype}" if celltype else ""
    print("\n" + "="*70)
    print(f"STEP {step_num}: {step_name}{suffix}")
    print("="*70)


def create_model_compatible_adata(adata_full, model_genes):
    """Create subset of adata matching model's gene set (order = model order)."""
    common_genes = [g for g in model_genes if g in adata_full.var_names]
    missing_genes = set(model_genes) - set(adata_full.var_names)

    print(f"  Model genes: {len(model_genes)}")
    print(f"  Full data genes: {adata_full.n_vars}")
    print(f"  Common genes: {len(common_genes)} ({len(common_genes)/len(model_genes)*100:.1f}%)")

    if missing_genes:
        print(f"  ⚠️ {len(missing_genes)} genes from model not in data")
        if len(missing_genes) <= 10:
            print(f"     Missing: {list(missing_genes)}")

    if len(common_genes) < len(model_genes) * 0.9:
        print(f"  ⚠️ WARNING: Only {len(common_genes)/len(model_genes)*100:.1f}% match")

    adata_model = adata_full[:, common_genes].copy()
    return adata_model


def load_data_and_model_production(celltype_name, config, output_dir):
    """
    Production-grade model loading with HVG extraction.

    Returns:
    adata_full, adata_model, lvae, cols, model_genes
    """
    log_step(1, "Loading Data and Model", celltype_name)

    h5ad_path = BASE_DIR / config['h5ad']
    print(f"\nLoading full data: {h5ad_path}")
    adata_full = sc.read_h5ad(h5ad_path)

    print(f"  Cells: {adata_full.n_obs:,}")
    print(f"  Genes: {adata_full.n_vars:,}")

    cols = get_column_names(adata_full, config['naming_mode'])
    print(f"\nColumn mapping:")
    for key, val in cols.items():
        print(f"  {key}: {val}")

    try:
        validate_columns_strict(adata_full, cols)
    except ValueError as e:
        print(f"\n❌ Column validation failed: {e}")
        return adata_full, None, None, cols, None

    scanvi_path = BASE_DIR / config['scanvi_model']
    print(f"\nAttempting scANVI model load: {scanvi_path}")

    lvae = None
    adata_model = None
    model_genes = None

    try:
        print("\n--- Extracting training genes from model ---")
        model_genes = extract_var_names_from_model(scanvi_path, model_type='SCANVI')
        if model_genes is None:
            raise ValueError("Could not extract gene list from model")

        gene_list_path = output_dir / f'{celltype_name}_model_genes.txt'
        save_hvg_list(model_genes, gene_list_path)

        print("\n--- Creating model-compatible subset ---")
        adata_model = create_model_compatible_adata(adata_full, model_genes)

        print("\n--- Loading scANVI model ---")
        lvae = scvi.model.SCANVI.load(scanvi_path, adata=adata_model)

        print("  ✓ Model loaded successfully")
        print(f"  Model adata: {adata_model.n_obs:,} cells × {adata_model.n_vars} genes")

    except Exception as e:
        print(f"\n⚠️ Model loading failed: {e}")
        print("  Proceeding with scanpy-only analysis")
        lvae = None
        adata_model = None

    return adata_full, adata_model, lvae, cols, model_genes


# ============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION (ROBUST)
# ============================================================================

def _filter_de_significant(de: pd.DataFrame):
    """Cross-version robust significance filter for scvi-tools DE outputs."""
    bf_ok = None
    if 'bayes_factor' in de.columns:
        bf_ok = de['bayes_factor'].astype(float) > float(DE_PARAMS['bayes_factor_threshold'])

    proba_col = _get_proba_de_col(de)
    proba_ok = None
    if proba_col is not None:
        proba_ok = de[proba_col].astype(float) > float(DE_PARAMS['proba_de_threshold'])

    # LFC
    if 'lfc_mean' in de.columns:
        lfc = de['lfc_mean'].astype(float)
        lfc_ok = lfc > float(DE_PARAMS['lfc_threshold'])
    else:
        lfc_array = _compute_log2fc_from_means(de)
        if lfc_array is not None:
            # ⭐ FIX: Preserve index alignment by creating Series
            lfc = pd.Series(lfc_array, index=de.index, name='lfc_approx')
            lfc_ok = lfc > float(DE_PARAMS['lfc_threshold'])
        else:
            lfc = None
            lfc_ok = None

    # Combine
    mask = None
    for m in [bf_ok, proba_ok, lfc_ok]:
        if m is None:
            continue
        mask = m if mask is None else (mask & m)

    if mask is None:
        # nothing to filter on; return empty
        return de.iloc[0:0].copy()

    de_sig = de.loc[mask].copy()

    if 'non_zeros_proportion1' in de_sig.columns:
        de_sig = de_sig[de_sig['non_zeros_proportion1'].astype(float) > float(DE_PARAMS['min_proportion'])]

    # Attach computed lfc if it was derived
    if 'lfc_mean' not in de_sig.columns and lfc is not None:
        # ⭐ FIX: lfc is already a Series with proper index, just join it
        if isinstance(lfc, pd.Series):
            de_sig = de_sig.join(lfc.rename('lfc_approx_log2'), how='left')
        else:
            # Fallback for ndarray (should not happen with above fix)
            tmp = pd.Series(lfc, index=de.index, name='lfc_approx_log2')
            de_sig = de_sig.join(tmp.loc[de_sig.index], how='left')

    return de_sig


def perform_de_one_vs_rest_explicit(adata_model, lvae, celltype_col, output_dir):
    """
    Explicit loop DE without comparison string dependency.

    Note: per scvi-tools API, group2=None is one-vs-rest behaviour. citeturn0search0
    """
    print("\n--- One-vs-Rest DE (Explicit Loop) ---")

    celltypes = (adata_model.obs[celltype_col].cat.categories
                 if hasattr(adata_model.obs[celltype_col], 'cat')
                 else pd.unique(adata_model.obs[celltype_col]))

    all_results = []

    for ct in celltypes:
        print(f"\n  Processing: {ct}")
        try:
            de = lvae.differential_expression(
                groupby=celltype_col,
                group1=ct,
                group2=None,
                mode=DE_PARAMS['mode'],
                delta=DE_PARAMS['delta'],
                batch_correction=DE_PARAMS['batch_correction'],
                n_samples=DE_PARAMS['n_samples']
            )

            de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_model.var_names)

            de['cell_type'] = ct
            de['comparison'] = f'{ct}_vs_Rest'

            all_results.append(de)

            de_sig = _filter_de_significant(de)
            print(f"    Significant markers: {len(de_sig)}")

            safe_ct = safe_name(ct)
            de_sig.to_csv(output_dir / f'{safe_ct}_markers.csv', index=False)

        except Exception as e:
            print(f"    ⚠️ DE failed: {e}")
            continue

    if all_results:
        de_all = pd.concat(all_results, ignore_index=True)
        de_all.to_csv(output_dir / 'celltype_markers_all_scANVI.csv', index=False)
        print(f"\n  ✓ Saved: celltype_markers_all_scANVI.csv")
        return de_all

    return None


def perform_de_scanpy_fallback(adata, celltype_col, output_dir):
    """Scanpy fallback with proper expression source"""
    print("\n--- Scanpy Fallback DE ---")

    expr_source = resolve_expr_source(adata)
    print(f"  Expression source: {expr_source}")

    sc.tl.rank_genes_groups(
        adata,
        groupby=celltype_col,
        method='wilcoxon',
        n_genes=100,
        **expr_source
    )

    de_results = sc.get.rank_genes_groups_df(adata, group=None)
    de_results.to_csv(output_dir / 'celltype_markers_all_scanpy.csv', index=False)
    print("  ✓ Saved: celltype_markers_all_scanpy.csv")

    return de_results


def perform_de_disease_within_celltype(adata_for_de, lvae, cols, output_dir):
    """
    Disease vs Healthy DE within each cell type.

    P0 fix: enforce explicit two-group comparison (not one-vs-rest).
    In scvi-tools, group2=None performs one-vs-rest behaviour. citeturn0search0
    """
    if not cols['disease_col'] or cols['disease_col'] not in adata_for_de.obs.columns:
        print("\n⚠️ No disease column - skipping")
        return

    print("\n--- Disease DE (Within Cell Type; explicit 2-group) ---")
    print("⚠️ NOTE: Cell-level statistics - pseudobulk recommended for final inference")

    celltype_col = cols['celltype_col']
    disease_col = cols['disease_col']

    for ct in pd.unique(adata_for_de.obs[celltype_col]):
        print(f"\n  {ct}...")
        adata_ct = adata_for_de[adata_for_de.obs[celltype_col] == ct].copy()

        if adata_ct.n_obs < 100:
            print(f"    Skipping (n={adata_ct.n_obs})")
            continue

        states = list(pd.unique(adata_ct.obs[disease_col].astype(str)))
        if len(states) != 2:
            print(f"    Skipping (need exactly 2 disease states, got {states})")
            continue

        g1, g2 = states[0], states[1]
        safe_ct = safe_name(ct)

        if lvae is not None:
            try:
                de = lvae.differential_expression(
                    adata=adata_ct,
                    groupby=disease_col,
                    group1=g1,
                    group2=g2,
                    mode=DE_PARAMS['mode'],
                    delta=0.25,
                    batch_correction=True,
                    n_samples=DE_PARAMS['n_samples']
                )

                de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_ct.var_names)
                de.to_csv(output_dir / f'{safe_ct}_{safe_name(g1)}_vs_{safe_name(g2)}_scANVI.csv', index=False)

                de_sig = _filter_de_significant(de)
                print(f"    ✓ Significant: {len(de_sig)} (robust filter)")

            except Exception as e:
                print(f"    ⚠️ Failed: {e}")
        else:
            expr_source = resolve_expr_source(adata_ct)
            sc.tl.rank_genes_groups(adata_ct, groupby=disease_col, method='wilcoxon', **expr_source)
            de = sc.get.rank_genes_groups_df(adata_ct, group=None)
            de.to_csv(output_dir / f'{safe_ct}_{safe_name(g1)}_vs_{safe_name(g2)}_scanpy.csv', index=False)
            print(f"    ✓ Scanpy DE done")


def perform_de_tissue_within_celltype(adata_for_de, lvae, cols, tissue_pairs, output_dir):
    """
    Tissue comparison within each cell type.
    """
    if not cols['tissue_col'] or cols['tissue_col'] not in adata_for_de.obs.columns:
        print("\n⚠️ No tissue column - skipping")
        return

    if not tissue_pairs:
        print("\n⚠️ No tissue pairs configured - skipping")
        return

    print("\n--- Tissue Comparison DE (Within Cell Type) ---")
    print("⚠️ NOTE: Tissue often confounded with batch - pseudobulk recommended")

    celltype_col = cols['celltype_col']
    tissue_col = cols['tissue_col']

    tissue_dir = output_dir / 'tissue_comparison'
    tissue_dir.mkdir(exist_ok=True, parents=True)

    for tissue1, tissue2 in tissue_pairs:
        print(f"\n  Comparing: {tissue1} vs {tissue2}")

        for ct in pd.unique(adata_for_de.obs[celltype_col]):
            adata_ct = adata_for_de[adata_for_de.obs[celltype_col] == ct].copy()

            tissues_in_ct = pd.unique(adata_ct.obs[tissue_col].astype(str))
            if str(tissue1) not in set(tissues_in_ct) or str(tissue2) not in set(tissues_in_ct):
                continue

            if adata_ct.n_obs < 100:
                continue

            print(f"    {ct}...")

            if lvae is not None:
                try:
                    de = lvae.differential_expression(
                        adata=adata_ct,
                        groupby=tissue_col,
                        group1=tissue1,
                        group2=tissue2,
                        mode=DE_PARAMS['mode'],
                        delta=0.25,
                        batch_correction=True
                    )

                    de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_ct.var_names)

                    safe_ct = safe_name(ct)
                    filename = f'{safe_ct}_{safe_name(tissue1)}_vs_{safe_name(tissue2)}.csv'
                    de.to_csv(tissue_dir / filename, index=False)

                    de_sig = _filter_de_significant(de)
                    print(f"      ✓ Significant: {len(de_sig)}")

                except Exception as e:
                    print(f"      ⚠️ Failed: {e}")


def perform_differential_expression(adata_full, adata_model, lvae, cols, output_dir):
    """Complete DE workflow"""
    log_step(2, "Differential Expression Analysis", "")

    de_dir = output_dir / 'differential_expression'
    de_dir.mkdir(exist_ok=True, parents=True)

    celltype_col = cols['celltype_col']
    adata_for_de = adata_model if adata_model is not None else adata_full

    if lvae is not None and adata_model is not None:
        perform_de_one_vs_rest_explicit(adata_model, lvae, celltype_col, de_dir)
    else:
        perform_de_scanpy_fallback(adata_full, celltype_col, de_dir)

    perform_de_disease_within_celltype(adata_for_de, lvae, cols, de_dir)

    if TISSUE_COMPARISON_PAIRS:
        perform_de_tissue_within_celltype(adata_for_de, lvae, cols, TISSUE_COMPARISON_PAIRS, de_dir)

    print("\n✓ Differential expression complete")
    return de_dir


# ============================================================================
# STEP 3: DENOISED EXPRESSION (MEMORY-SAFE)
# ============================================================================

def generate_denoised_expression_safe(adata_model, lvae, output_dir):
    """Memory-safe denoised expression - NEVER writes to layers"""
    log_step(3, "Generating Denoised Expression", "")

    if lvae is None or adata_model is None:
        print("⚠️ Model not available - skipping")
        return None

    n_elements = adata_model.n_obs * adata_model.n_vars
    if n_elements > DENOISED_PARAMS['max_elements_in_memory']:
        print(f"⚠️ Dataset too large ({n_elements:.2e} > {DENOISED_PARAMS['max_elements_in_memory']:.2e})")
        print("  Skipping to prevent memory issues")
        return None

    print(f"\nGenerating denoised expression...")
    print(f"  n_samples: {DENOISED_PARAMS['n_samples']}")
    print(f"  library_size: {DENOISED_PARAMS['library_size']}")

    try:
        denoised = lvae.get_normalized_expression(
            adata=adata_model,
            n_samples=DENOISED_PARAMS['n_samples'],
            return_mean=True,
            library_size=DENOISED_PARAMS['library_size']
        )

        # ⭐ FIX: Handle both DataFrame and ndarray return types
        if isinstance(denoised, pd.DataFrame):
            denoised_array = denoised.values.astype(np.float32)
            # DataFrame columns are gene names, use them to subset var
            denoised_var_names = denoised.columns
            var_subset = adata_model.var.loc[denoised_var_names].copy()
        else:
            # Already ndarray - assume same order as adata_model.var_names
            denoised_array = denoised.astype(np.float32)
            var_subset = adata_model.var.copy()

        print(f"  ✓ Generated: {denoised_array.shape}")
        print(f"  Memory: {denoised_array.nbytes / 1e9:.2f} GB")

        print("\n  Saving as separate h5ad...")
        adata_denoised = sc.AnnData(
            X=denoised_array,
            obs=adata_model.obs.copy(),
            var=var_subset
        )

        denoised_path = output_dir / 'denoised_expression.h5ad'
        adata_denoised.write_h5ad(denoised_path, compression='gzip')
        print(f"  ✓ Saved: {denoised_path}")

        del adata_denoised
        gc.collect()

        return denoised_path

    except Exception as e:
        print(f"\n  ⚠️ Failed: {e}")
        return None


# ============================================================================
# STEP 4: CELL COMPOSITION ANALYSIS
# ============================================================================

def perform_cell_composition_analysis(adata, cols, output_dir):
    """Cell composition with tissue support"""
    log_step(4, "Cell Composition Analysis", "")

    comp_dir = output_dir / 'cell_composition'
    comp_dir.mkdir(exist_ok=True, parents=True)

    celltype_col = cols['celltype_col']
    sample_col = cols['sample_col']

    print("\n--- Overall Composition ---")
    composition = pd.crosstab(
        adata.obs[sample_col],
        adata.obs[celltype_col],
        normalize='index'
    ) * 100

    composition.to_csv(comp_dir / 'composition_by_sample.csv')

    fig, ax = plt.subplots(figsize=(12, 6))
    composition.mean().sort_values(ascending=False).plot(kind='bar', ax=ax)
    ax.set_ylabel('Cell Proportion (%)')
    ax.set_xlabel('Cell Type')
    ax.set_title('Average Cell Type Composition')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(comp_dir / 'composition_overall.png', dpi=300, bbox_inches='tight')
    plt.close()

    if cols['disease_col'] and cols['disease_col'] in adata.obs.columns:
        print("\n--- Disease Composition ---")

        sample_disease = adata.obs[[sample_col, cols['disease_col']]].drop_duplicates().set_index(sample_col)
        composition_with_disease = composition.join(sample_disease)

        comp_by_disease = composition_with_disease.groupby(cols['disease_col']).mean()
        comp_by_disease.to_csv(comp_dir / 'composition_by_disease.csv')

        comp_stats = []
        disease_states = list(pd.unique(adata.obs[cols['disease_col']].astype(str)))

        if len(disease_states) == 2:
            state1, state2 = disease_states

            for ct in composition.columns:
                group1 = composition_with_disease[composition_with_disease[cols['disease_col']] == state1][ct]
                group2 = composition_with_disease[composition_with_disease[cols['disease_col']] == state2][ct]

                if len(group1) > 0 and len(group2) > 0:
                    u_stat, p_val = mannwhitneyu(group1, group2)

                    comp_stats.append({
                        'cell_type': ct,
                        f'{safe_name(state1)}_mean': group1.mean(),
                        f'{safe_name(state2)}_mean': group2.mean(),
                        'fold_change': group1.mean() / (group2.mean() + 1e-10),
                        'p_value': p_val
                    })

            if comp_stats:
                comp_stats_df = pd.DataFrame(comp_stats)
                comp_stats_df['p_adj'] = multipletests(comp_stats_df['p_value'], method='fdr_bh')[1]
                comp_stats_df = comp_stats_df.sort_values('p_adj')
                comp_stats_df.to_csv(comp_dir / 'composition_disease_statistics.csv', index=False)

                sig = comp_stats_df[comp_stats_df['p_adj'] < 0.05]
                if len(sig) > 0:
                    print(f"\n  Significant changes (FDR < 0.05): {len(sig)}")

    if cols['tissue_col'] and cols['tissue_col'] in adata.obs.columns:
        print("\n--- Tissue Composition ---")

        sample_tissue = adata.obs[[sample_col, cols['tissue_col']]].drop_duplicates().set_index(sample_col)
        composition_with_tissue = composition.join(sample_tissue)

        comp_by_tissue = composition_with_tissue.groupby(cols['tissue_col']).mean()
        comp_by_tissue.to_csv(comp_dir / 'composition_by_tissue.csv')
        print(f"  ✓ Saved: composition_by_tissue.csv")

    print("\n✓ Cell composition complete")
    return comp_dir


# ============================================================================
# STEP 5: PSEUDOBULK EXPORT (COUNTS-SAFE + NO SIDE EFFECTS)
# ============================================================================

def export_pseudobulk_enhanced(adata, cols, output_dir):
    """
    Enhanced pseudobulk export with metadata preservation.
    P0 fix: counts-safe matrix selection (prefer layers['counts']).
    P2 fix: no writes to adata.obs (no pseudobulk_key side effect).
    """
    log_step(5, "Pseudobulk Export (Donor-Level)", "")

    pseudobulk_dir = output_dir / 'pseudobulk'
    pseudobulk_dir.mkdir(exist_ok=True, parents=True)

    celltype_col = cols['celltype_col']
    sample_col = cols['sample_col']

    counts, genes = get_counts_matrix(adata)
    if counts is None:
        print("⚠️ No counts matrix - skipping")
        return None

    print(f"\nAggregating by {sample_col} × {celltype_col}...")

    key_series = (
        adata.obs[sample_col].astype(str) + '|' +
        adata.obs[celltype_col].astype(str)
    )

    pseudobulk_list = []
    metadata_list = []

    for key in pd.unique(key_series):
        mask = (key_series == key).values
        n_cells = int(mask.sum())

        if n_cells < PSEUDOBULK_PARAMS['min_cells_per_sample']:
            continue

        # ⭐ FIX: Handle sparse matrix properly
        counts_subset = counts[mask, :]
        if sp.issparse(counts_subset):
            if PSEUDOBULK_PARAMS['aggregation'] == 'sum':
                pseudo_counts = np.array(counts_subset.sum(axis=0)).ravel()
            else:
                pseudo_counts = np.array(counts_subset.mean(axis=0)).ravel()
        else:
            if PSEUDOBULK_PARAMS['aggregation'] == 'sum':
                pseudo_counts = np.array(counts_subset.sum(axis=0)).ravel()
            else:
                pseudo_counts = np.array(counts_subset.mean(axis=0)).ravel()

        pseudobulk_list.append(pseudo_counts)

        meta = adata.obs.loc[mask].iloc[0]
        metadata_list.append({
            'pseudobulk_id': key,
            'sample_id': meta[sample_col],
            'cell_type': meta[celltype_col],
            'n_cells': n_cells,
            'batch': meta[cols['batch_col']] if cols['batch_col'] in meta.index else None,
            'disease_status': meta[cols['disease_col']] if cols['disease_col'] and cols['disease_col'] in meta.index else None,
            'tissue': meta[cols['tissue_col']] if cols['tissue_col'] and cols['tissue_col'] in meta.index else None
        })

    if len(pseudobulk_list) == 0:
        print("⚠️ No pseudobulk samples passed filters")
        return None

    pseudobulk_matrix = np.vstack(pseudobulk_list)
    pseudobulk_metadata = pd.DataFrame(metadata_list)

    adata_pseudo = sc.AnnData(
        X=pseudobulk_matrix.astype(np.int32),
        obs=pseudobulk_metadata,
        var=pd.DataFrame(index=genes)
    )

    print(f"\n  Generated {adata_pseudo.n_obs} pseudobulk samples")
    print(f"  Genes: {adata_pseudo.n_vars}")

    adata_pseudo.write_h5ad(pseudobulk_dir / 'pseudobulk_counts.h5ad', compression='gzip')
    pseudobulk_metadata.to_csv(pseudobulk_dir / 'pseudobulk_metadata.csv', index=False)

    summary = []
    summary.append("Pseudobulk Export Summary")
    summary.append("="*50)
    summary.append(f"Aggregation: {PSEUDOBULK_PARAMS['aggregation']}")
    summary.append(f"Min cells/sample: {PSEUDOBULK_PARAMS['min_cells_per_sample']}")
    summary.append(f"Total samples: {adata_pseudo.n_obs}")
    summary.append(f"\nSamples per cell type:")
    for ct, count in pseudobulk_metadata['cell_type'].value_counts().items():
        summary.append(f"  {ct}: {count}")

    summary_text = '\n'.join(summary)
    with open(pseudobulk_dir / 'pseudobulk_summary.txt', 'w') as f:
        f.write(summary_text)

    print(f"\n{summary_text}")
    print("\n⚠️ CRITICAL NOTE:")
    print("  For final statistical inference, use donor-level analysis:")
    print("  - edgeR/DESeq2 (R) or pydeseq2 (Python) on pseudobulk_counts.h5ad")
    print("  Cell-level DE p-values are inflated and should be considered exploratory")

    print("\n✓ Pseudobulk export complete")
    return pseudobulk_dir


# ============================================================================
# STEP 6: SUBCLUSTERING (MEMORY-SAFE)
# ============================================================================

def perform_subclustering_safe(adata, cols, output_dir):
    """Memory-safe subclustering - stores paths not objects"""
    log_step(6, "Subclustering on scANVI Latent Space", "")

    subcluster_dir = output_dir / 'subclustering'
    subcluster_dir.mkdir(exist_ok=True, parents=True)

    celltype_col = cols['celltype_col']

    if 'X_scanvi' not in adata.obsm:
        print("⚠️ X_scanvi not found - skipping")
        return None

    celltype_counts = adata.obs[celltype_col].value_counts()
    celltypes_to_process = celltype_counts[celltype_counts >= SUBCLUSTER_PARAMS['min_cells_for_subcluster']].index.tolist()

    print(f"\nCell types (≥{SUBCLUSTER_PARAMS['min_cells_for_subcluster']} cells):")
    for ct in celltypes_to_process:
        print(f"  {ct}: {celltype_counts[ct]:,}")

    if len(celltypes_to_process) == 0:
        print("\n⚠️ No cell types meet threshold")
        return None

    results_summary = {}

    for ct in celltypes_to_process:
        print(f"\n{'='*70}")
        print(f"Subclustering: {ct}")
        print(f"{'='*70}")

        ct_safe = safe_name(ct)
        ct_dir = subcluster_dir / ct_safe
        ct_dir.mkdir(exist_ok=True, parents=True)

        adata_sub = adata[adata.obs[celltype_col] == ct].copy()
        print(f"  Cells: {adata_sub.n_obs:,}")

        n_neighbors = min(SUBCLUSTER_PARAMS['n_neighbors'], max(10, adata_sub.n_obs // 10))
        sc.pp.neighbors(adata_sub, use_rep='X_scanvi', n_neighbors=n_neighbors)

        print("  Clustering...")
        for res in SUBCLUSTER_PARAMS['resolutions']:
            key = f'leiden_sub_res{res}'
            sc.tl.leiden(adata_sub, resolution=res, key_added=key)
            n_clusters = adata_sub.obs[key].nunique()
            print(f"    Res {res}: {n_clusters} clusters")

        sc.tl.umap(adata_sub, min_dist=0.3)

        default_key = f'leiden_sub_res{SUBCLUSTER_PARAMS["default_resolution"]}'
        expr_source = resolve_expr_source(adata_sub)

        sc.tl.rank_genes_groups(adata_sub, groupby=default_key, method='wilcoxon', n_genes=100, **expr_source)
        markers_df = sc.get.rank_genes_groups_df(adata_sub, group=None)
        markers_df.to_csv(ct_dir / 'subcluster_markers.csv', index=False)

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        sc.pl.umap(adata_sub, color=celltype_col, ax=axes[0, 0], show=False, title='Original', frameon=False, size=50)
        sc.pl.umap(adata_sub, color=cols['batch_col'], ax=axes[0, 1], show=False, title='Batch', frameon=False, size=50)

        if cols['disease_col'] and cols['disease_col'] in adata_sub.obs.columns:
            sc.pl.umap(adata_sub, color=cols['disease_col'], ax=axes[0, 2], show=False, title='Disease', frameon=False, size=50)
        else:
            axes[0, 2].axis('off')

        for i, res in enumerate([0.2, 0.4, 0.6]):
            key = f'leiden_sub_res{res}'
            sc.pl.umap(adata_sub, color=key, ax=axes[1, i], show=False,
                       title=f'Res {res}', legend_loc='right margin', frameon=False, size=50)

        plt.tight_layout()
        plt.savefig(ct_dir / 'subclustering_overview.png', dpi=300, bbox_inches='tight')
        plt.close()

        output_path = ct_dir / f'{ct_safe}_subclustered.h5ad'
        adata_sub.write_h5ad(output_path, compression='gzip')

        results_summary[str(ct)] = {
            'path': str(output_path),
            'n_cells': int(adata_sub.n_obs),
            'n_subclusters': int(adata_sub.obs[default_key].nunique())
        }

        print(f"  ✓ Complete")

        del adata_sub
        gc.collect()

    with open(subcluster_dir / 'subclustering_summary.json', 'w') as f:
        json.dump(results_summary, f, indent=2)

    print(f"\n{'='*70}")
    print(f"✓ Subclustering complete: {len(results_summary)} cell types")
    print(f"{'='*70}")

    return results_summary


# ============================================================================
# SUMMARY GENERATION
# ============================================================================

def save_summary_stats(adata, cols, model_genes, output_dir):
    """Enhanced summary with model info"""
    summary = []
    summary.append("="*70)
    summary.append("Dataset Summary Statistics")
    summary.append("="*70)
    summary.append(f"\nTotal cells: {adata.n_obs:,}")
    summary.append(f"Total genes: {adata.n_vars:,}")

    if model_genes:
        summary.append(f"Model training genes: {len(model_genes)}")
        overlap = len(set(model_genes) & set(adata.var_names))
        summary.append(f"  Overlap with data: {overlap} ({overlap/len(model_genes)*100:.1f}%)")

    summary.append(f"\nCell Type Distribution ({cols['celltype_col']}):")
    for ct, count in adata.obs[cols['celltype_col']].value_counts().items():
        pct = count / adata.n_obs * 100
        summary.append(f"  {ct}: {count:,} ({pct:.1f}%)")

    summary.append(f"\nscANVI Confidence:")
    conf = adata.obs[cols['confidence_col']].describe()
    summary.append(f"  Mean: {conf['mean']:.3f}")
    summary.append(f"  Median: {conf['50%']:.3f}")
    low_conf = (adata.obs[cols['confidence_col']] < 0.5).sum()
    summary.append(f"  Low (<0.5): {low_conf:,} ({low_conf/adata.n_obs*100:.1f}%)")

    if cols['batch_col'] in adata.obs.columns:
        summary.append(f"\nBatches: {adata.obs[cols['batch_col']].nunique()}")

    if cols['tissue_col'] and cols['tissue_col'] in adata.obs.columns:
        summary.append(f"\nTissue Distribution:")
        for tissue, count in adata.obs[cols['tissue_col']].value_counts().items():
            pct = count / adata.n_obs * 100
            summary.append(f"  {tissue}: {count:,} ({pct:.1f}%)")

    if cols['disease_col'] and cols['disease_col'] in adata.obs.columns:
        summary.append(f"\nDisease Status:")
        for status, count in adata.obs[cols['disease_col']].value_counts().items():
            pct = count / adata.n_obs * 100
            summary.append(f"  {status}: {count:,} ({pct:.1f}%)")

    summary.append("\n" + "="*70)

    summary_text = '\n'.join(summary)
    with open(output_dir / 'summary_statistics.txt', 'w') as f:
        f.write(summary_text)

    print(summary_text)


# ============================================================================
# TISSUE-STRATIFIED RUNNER
# ============================================================================

def iter_tissue_subsets(adata_full, adata_model, cols):
    """
    Yield (subset_name, adata_full_subset, adata_model_subset).
    Uses thresholds to avoid tiny/non-informative tissues.
    """
    tissue_col = cols.get('tissue_col')
    sample_col = cols.get('sample_col')

    if (not TISSUE_STRATIFY_PARAMS["enabled"]) or (not tissue_col) or (tissue_col not in adata_full.obs.columns):
        return [("ALL", adata_full, adata_model)]

    tissues = pd.unique(adata_full.obs[tissue_col].astype(str))
    if TISSUE_STRATIFY_PARAMS["tissue_whitelist"] is not None:
        allowed = set(map(str, TISSUE_STRATIFY_PARAMS["tissue_whitelist"]))
        tissues = [t for t in tissues if str(t) in allowed]

    subsets = []
    for t in tissues:
        mask = (adata_full.obs[tissue_col].astype(str) == str(t)).values
        n_cells = int(mask.sum())
        if n_cells < int(TISSUE_STRATIFY_PARAMS["min_cells_per_tissue"]):
            continue

        if sample_col and sample_col in adata_full.obs.columns:
            n_samples = int(adata_full.obs.loc[mask, sample_col].nunique())
            if n_samples < int(TISSUE_STRATIFY_PARAMS["min_samples_per_tissue"]):
                continue

        # ⭐ FIX: Use copy instead of view to avoid issues with subsequent operations
        ad_full_t = adata_full[mask].copy()
        ad_model_t = adata_model[mask].copy() if adata_model is not None else None
        subsets.append((t, ad_full_t, ad_model_t))

    if TISSUE_STRATIFY_PARAMS["run_global_all"]:
        subsets = [("ALL", adata_full, adata_model)] + subsets

    return subsets


def run_downstream_for_subset(subset_name, adata_full_sub, adata_model_sub, lvae, cols, model_genes, out_dir):
    """Run steps 2-6 for a given subset (e.g., tissue)."""
    out_dir.mkdir(exist_ok=True, parents=True)

    save_summary_stats(adata_full_sub, cols, model_genes, out_dir)
    perform_differential_expression(adata_full_sub, adata_model_sub, lvae, cols, out_dir)
    generate_denoised_expression_safe(adata_model_sub, lvae, out_dir)
    perform_cell_composition_analysis(adata_full_sub, cols, out_dir)
    export_pseudobulk_enhanced(adata_full_sub, cols, out_dir)
    perform_subclustering_safe(adata_full_sub, cols, out_dir)

    final_path = out_dir / f'{safe_name(subset_name)}_processed_final.h5ad'
    adata_full_sub.copy().write_h5ad(final_path, compression='gzip')
    print(f"  ✓ Saved subset final: {final_path}")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def process_celltype(celltype_name, config):
    """Main processing pipeline"""
    print("\n" + "="*70)
    print(f"PROCESSING: {celltype_name}")
    print("="*70)

    output_dir = OUTPUT_BASE / celltype_name
    output_dir.mkdir(exist_ok=True, parents=True)

    try:
        # Step 1: Load
        adata_full, adata_model, lvae, cols, model_genes = load_data_and_model_production(
            celltype_name, config, output_dir
        )

        # Tissue-stratified (and/or global) downstream
        subsets = iter_tissue_subsets(adata_full, adata_model, cols)

        print("\n" + "="*70)
        print(f"SUBSETS TO PROCESS ({celltype_name}): {len(subsets)}")
        print("="*70)

        for tname, ad_full_sub, ad_model_sub in subsets:
            sub_dir = output_dir / ("ALL" if tname == "ALL" else f"tissue_{safe_name(tname)}")
            print("\n" + "-"*70)
            print(f"Running subset: {tname}  ->  {sub_dir}")
            print("-"*70)

            run_downstream_for_subset(
                subset_name=tname,
                adata_full_sub=ad_full_sub,
                adata_model_sub=ad_model_sub,
                lvae=lvae,
                cols=cols,
                model_genes=model_genes,
                out_dir=sub_dir
            )

        print("\n" + "="*70)
        print(f"✓ {celltype_name} COMPLETE")
        print("="*70)

        # Cleanup
        del adata_full
        if adata_model is not None:
            del adata_model
        if lvae is not None:
            del lvae
        gc.collect()

        return True

    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main entry point"""
    print("="*70)
    print("scANVI DOWNSTREAM ANALYSIS - v1.3 PRODUCTION")
    print("="*70)
    print(f"\nBase: {BASE_DIR}")
    print(f"Output: {OUTPUT_BASE}")

    OUTPUT_BASE.mkdir(exist_ok=True, parents=True)

    results = {}
    for celltype_name, config in CELL_TYPES.items():
        success = process_celltype(celltype_name, config)
        results[celltype_name] = success

    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)

    for celltype, success in results.items():
        status = "✓ SUCCESS" if success else "❌ FAILED"
        print(f"  {celltype}: {status}")

    successful = sum(results.values())
    print(f"\nTotal: {successful}/{len(results)} successful")
    print(f"\nOutputs: {OUTPUT_BASE}")
    print("="*70)


if __name__ == "__main__":
    main()
