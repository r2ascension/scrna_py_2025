#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Step 1: Add UMAP Operator to scanvi_existing_model (ONE-TIME SETUP)

Purpose:
- Load reference scANVI model and data
- Fit UMAP on reference latent space
- Save UMAP operator to model for future query projections

⚠️  Run this script ONCE before any query mapping
✓  After running, all queries can use .transform() for aligned UMAP

Author: r2end
Date: 2025-01-13
Version: 1.0
"""

import scanpy as sc
import scvi
import umap
from pathlib import Path
import numpy as np
import sys

# ==============================================================================
# Configuration Section (MODIFY THESE PATHS)
# ==============================================================================

# ⚠️  CRITICAL: Update these paths to match your setup
REF_H5AD = "/home/h2048/data/py/1207/allcells_scvi_analysis/adata_allcells_dual_scanvi_final_v2.2.0.h5ad"
SCANVI_MODEL_DIR = "/home/h2048/data/py/1207/allcells_scvi_analysis/models/scanvi_existing_model"

# Keys (must match your training script)
LATENT_KEY = "X_scanvi_existing"
UMAP_KEY = "X_umap_scanvi_existing"

# UMAP parameters (must match your training script)
UMAP_N_NEIGHBORS = 50
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0

# ==============================================================================
# Main Script
# ==============================================================================

def main():
    print("=" * 70)
    print("Adding UMAP Operator to scanvi_existing_model")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Reference data: {REF_H5AD}")
    print(f"  Model directory: {SCANVI_MODEL_DIR}")
    print(f"  Latent key: {LATENT_KEY}")
    print(f"  UMAP key: {UMAP_KEY}")
    
    # Validate paths
    if not Path(REF_H5AD).exists():
        print(f"\n❌ ERROR: Reference h5ad not found at {REF_H5AD}")
        print("   Please update REF_H5AD in the configuration section")
        sys.exit(1)
    
    if not Path(SCANVI_MODEL_DIR).exists():
        print(f"\n❌ ERROR: Model directory not found at {SCANVI_MODEL_DIR}")
        print("   Please update SCANVI_MODEL_DIR in the configuration section")
        sys.exit(1)
    
    # ==============================================================================
    # Step 1: Load Reference Data
    # ==============================================================================
    
    print(f"\n[1/6] Loading reference data...")
    print(f"      Path: {REF_H5AD}")
    
    try:
        adata_ref = sc.read_h5ad(REF_H5AD)
        print(f"      ✓ Loaded: {adata_ref.n_obs:,} cells × {adata_ref.n_vars:,} genes")
    except Exception as e:
        print(f"\n❌ ERROR loading reference data: {e}")
        sys.exit(1)
    
    # Verify required keys
    if LATENT_KEY not in adata_ref.obsm:
        print(f"\n❌ ERROR: Latent key '{LATENT_KEY}' not found in reference data")
        print(f"   Available obsm keys: {list(adata_ref.obsm.keys())}")
        sys.exit(1)
    
    if UMAP_KEY not in adata_ref.obsm:
        print(f"\n❌ ERROR: UMAP key '{UMAP_KEY}' not found in reference data")
        print(f"   Available obsm keys: {list(adata_ref.obsm.keys())}")
        sys.exit(1)
    
    print(f"      ✓ Found latent: {LATENT_KEY} {adata_ref.obsm[LATENT_KEY].shape}")
    print(f"      ✓ Found UMAP: {UMAP_KEY} {adata_ref.obsm[UMAP_KEY].shape}")
    
    # ==============================================================================
    # Step 2: Load scANVI Model
    # ==============================================================================
    
    print(f"\n[2/6] Loading scANVI model...")
    print(f"      Path: {SCANVI_MODEL_DIR}")
    
    try:
        # Important: pass adata to load() to ensure registry alignment
        scanvi_model = scvi.model.SCANVI.load(SCANVI_MODEL_DIR, adata=adata_ref)
        print(f"      ✓ Model loaded successfully")
    except Exception as e:
        print(f"\n❌ ERROR loading model: {e}")
        sys.exit(1)
    
    # ==============================================================================
    # Step 3: Extract Reference Latent
    # ==============================================================================
    
    print(f"\n[3/6] Extracting reference latent representation...")
    
    # Use existing latent from h5ad (already computed during training)
    Z_ref = adata_ref.obsm[LATENT_KEY]
    print(f"      ✓ Latent shape: {Z_ref.shape}")
    
    # Sanity check: no NaNs or Infs
    if np.isnan(Z_ref).any():
        print(f"\n❌ ERROR: Reference latent contains NaN values")
        sys.exit(1)
    if np.isinf(Z_ref).any():
        print(f"\n❌ ERROR: Reference latent contains Inf values")
        sys.exit(1)
    
    print(f"      ✓ Latent validation passed")
    
    # ==============================================================================
    # Step 4: Fit UMAP Operator
    # ==============================================================================
    
    print(f"\n[4/6] Fitting UMAP operator...")
    print(f"      Parameters:")
    print(f"        - n_neighbors: {UMAP_N_NEIGHBORS}")
    print(f"        - min_dist: {UMAP_MIN_DIST}")
    print(f"        - spread: {UMAP_SPREAD}")
    
    try:
        umap_op = umap.UMAP(
            n_neighbors=UMAP_N_NEIGHBORS,
            min_dist=UMAP_MIN_DIST,
            spread=UMAP_SPREAD,
            random_state=42,
            verbose=False
        )
        
        # Fit on reference latent
        print(f"      Fitting (this may take a few minutes)...")
        umap_embedding = umap_op.fit_transform(Z_ref)
        print(f"      ✓ UMAP fitted: {umap_embedding.shape}")
    except Exception as e:
        print(f"\n❌ ERROR fitting UMAP: {e}")
        sys.exit(1)
    
    # Verify fit matches existing UMAP (should be very close)
    existing_umap = adata_ref.obsm[UMAP_KEY]
    correlation = np.corrcoef(umap_embedding.flatten(), existing_umap.flatten())[0, 1]
    print(f"      ✓ Correlation with existing UMAP: {correlation:.4f}")
    
    if correlation < 0.90:
        print(f"\n⚠️  WARNING: Low correlation with existing UMAP ({correlation:.4f})")
        print(f"   This might indicate parameter mismatch with training")
        print(f"   Proceeding anyway, but verify UMAP parameters match training")
    else:
        print(f"      ✓ High correlation - parameters match training")
    
    # ==============================================================================
    # Step 5: Save UMAP Operator to Model
    # ==============================================================================
    
    print(f"\n[5/6] Saving UMAP operator to model...")
    
    try:
        # Attach UMAP operator to model (with trailing underscore to persist)
        scanvi_model.umap_op_ = umap_op
        
        # Save model (overwrite with new attribute)
        scanvi_model.save(SCANVI_MODEL_DIR, overwrite=True)
        print(f"      ✓ Model saved with UMAP operator")
    except Exception as e:
        print(f"\n❌ ERROR saving model: {e}")
        sys.exit(1)
    
    # ==============================================================================
    # Step 6: Verification
    # ==============================================================================
    
    print(f"\n[6/6] Verification...")
    
    try:
        # Reload model to verify UMAP operator persisted
        scanvi_model_reloaded = scvi.model.SCANVI.load(SCANVI_MODEL_DIR)
        has_umap_op = hasattr(scanvi_model_reloaded, 'umap_op_')
        
        if has_umap_op:
            print(f"      ✓ UMAP operator successfully persisted")
            
            # Test transform on a small subset
            test_latent = Z_ref[:100]
            test_umap = scanvi_model_reloaded.umap_op_.transform(test_latent)
            print(f"      ✓ Test transform successful: {test_umap.shape}")
            
            # Verify transformation is deterministic
            test_umap2 = scanvi_model_reloaded.umap_op_.transform(test_latent)
            is_deterministic = np.allclose(test_umap, test_umap2)
            print(f"      ✓ Transform is deterministic: {is_deterministic}")
        else:
            print(f"\n❌ ERROR: UMAP operator not found after reload")
            sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERROR during verification: {e}")
        sys.exit(1)
    
    # ==============================================================================
    # Summary
    # ==============================================================================
    
    print("\n" + "=" * 70)
    print("✅ UMAP OPERATOR SUCCESSFULLY ADDED")
    print("=" * 70)
    print(f"\nModel location: {SCANVI_MODEL_DIR}")
    print(f"\nWhat's new:")
    print(f"  - Model now has 'umap_op_' attribute")
    print(f"  - All future queries can use .transform() for aligned UMAP")
    print(f"  - No need to recompute UMAP from scratch")
    print(f"\nNext step:")
    print(f"  - Run step2_scarches_mapping.py to map your query data")
    print("=" * 70)


if __name__ == "__main__":
    main()
