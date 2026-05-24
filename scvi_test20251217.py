#!/usr/bin/env python3
"""
Diagnostic Script: Extract Genes from scVI/scANVI Models
=========================================================

Test different methods to extract gene lists from saved models
"""

import sys
from pathlib import Path
import numpy as np
import scanpy as sc
import scvi

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("/home/h2048/data/core_data")

# Test both model types
MODELS_TO_TEST = {
    'T_scVI': BASE_DIR / 'models/T_models/scvi_model',
    'T_scANVI': BASE_DIR / 'models/T_models/scanvi_model',
}

# Also load the h5ad to compare
H5AD_PATH = BASE_DIR / 'adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad'

# ============================================================================
# METHOD 1: Load model without adata and use get_var_names()
# ============================================================================

def method1_get_var_names(model_path, model_type='SCANVI'):
    """Load model and use get_var_names() API"""
    print(f"\n{'='*70}")
    print(f"METHOD 1: model.get_var_names()")
    print(f"{'='*70}")
    
    try:
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
        
        print(f"  Loading {model_type} model...")
        model = Model.load(model_path, adata=None)
        
        print(f"  Calling get_var_names()...")
        var_names_dict = model.get_var_names()
        
        print(f"  Result type: {type(var_names_dict)}")
        print(f"  Keys: {list(var_names_dict.keys()) if isinstance(var_names_dict, dict) else 'N/A'}")
        
        if isinstance(var_names_dict, dict):
            for modality, genes in var_names_dict.items():
                print(f"\n  Modality: {modality}")
                print(f"    Type: {type(genes)}")
                print(f"    Length: {len(genes) if hasattr(genes, '__len__') else 'N/A'}")
                if hasattr(genes, '__len__') and len(genes) > 0:
                    print(f"    First 5: {list(genes)[:5]}")
                    print(f"    Last 5: {list(genes)[-5:]}")
                    return list(genes)
        
        print(f"  ❌ Cannot extract genes from get_var_names() output")
        return None
        
    except Exception as e:
        print(f"  ❌ Failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# METHOD 2: Load registry and search for var_names
# ============================================================================

def method2_registry_search(model_path, model_type='SCANVI'):
    """Load registry and recursively search for var_names"""
    print(f"\n{'='*70}")
    print(f"METHOD 2: Load registry and search")
    print(f"{'='*70}")
    
    try:
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
        
        print(f"  Loading registry...")
        reg = Model.load_registry(model_path)
        
        print(f"  Registry type: {type(reg)}")
        print(f"  Registry keys (if dict): {list(reg.keys()) if isinstance(reg, dict) else 'N/A'}")
        
        # Recursive search
        def find_var_names(obj, path="root", depth=0, max_depth=10):
            if depth > max_depth:
                return []
            
            results = []
            
            if isinstance(obj, dict):
                for k, v in obj.items():
                    current_path = f"{path}.{k}"
                    
                    # Check if this looks like var_names
                    if 'var' in str(k).lower() and isinstance(v, (list, np.ndarray)):
                        if len(v) > 100:  # Likely gene list
                            results.append((current_path, v, len(v)))
                    
                    # Recurse
                    results.extend(find_var_names(v, current_path, depth + 1))
                    
            elif isinstance(obj, list) and len(obj) < 100:  # Don't recurse into large lists
                for i, item in enumerate(obj):
                    results.extend(find_var_names(item, f"{path}[{i}]", depth + 1))
            
            return results
        
        print(f"\n  Searching for var_names-like structures...")
        candidates = find_var_names(reg)
        
        if candidates:
            print(f"\n  Found {len(candidates)} candidate(s):")
            for path, genes, length in candidates:
                print(f"\n    Path: {path}")
                print(f"    Length: {length}")
                if hasattr(genes, '__iter__'):
                    genes_list = list(genes)
                    print(f"    First 5: {genes_list[:5]}")
                    print(f"    Last 5: {genes_list[-5:]}")
                    
                    # Return the first valid one
                    if length > 100:
                        return genes_list
        else:
            print(f"  ❌ No var_names candidates found")
        
        return None
        
    except Exception as e:
        print(f"  ❌ Failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# METHOD 3: Load model with adata and check adata_manager
# ============================================================================

def method3_adata_manager(model_path, adata, model_type='SCANVI'):
    """Load model with adata and inspect adata_manager"""
    print(f"\n{'='*70}")
    print(f"METHOD 3: Load with adata and check adata_manager")
    print(f"{'='*70}")
    
    try:
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
        
        print(f"  Loading model with full adata...")
        model = Model.load(model_path, adata=adata)
        
        print(f"  Model loaded successfully")
        
        # Check adata_manager
        if hasattr(model, 'adata_manager'):
            print(f"\n  Checking adata_manager...")
            manager = model.adata_manager
            print(f"    Type: {type(manager)}")
            
            # Check data_registry
            if hasattr(manager, 'data_registry'):
                print(f"\n    data_registry found")
                reg = manager.data_registry
                print(f"      Type: {type(reg)}")
                
                # Check for var_names
                if hasattr(reg, 'var_names'):
                    var_names = reg.var_names
                    print(f"      var_names type: {type(var_names)}")
                    print(f"      Length: {len(var_names)}")
                    print(f"      First 5: {list(var_names)[:5]}")
                    return list(var_names)
        
        # Alternative: check adata used by model
        if hasattr(model, 'adata'):
            print(f"\n  Checking model.adata...")
            model_adata = model.adata
            print(f"    Genes in model.adata: {model_adata.n_vars}")
            print(f"    First 5: {list(model_adata.var_names)[:5]}")
            return list(model_adata.var_names)
        
        print(f"  ❌ Cannot find var_names in model")
        return None
        
    except Exception as e:
        print(f"  ❌ Failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# METHOD 4: Direct file inspection (model.pt)
# ============================================================================

def method4_direct_file(model_path):
    """Try to load model.pt file directly"""
    print(f"\n{'='*70}")
    print(f"METHOD 4: Direct model.pt inspection")
    print(f"{'='*70}")
    
    try:
        import torch
        
        model_pt = Path(model_path) / 'model.pt'
        print(f"  Loading: {model_pt}")
        
        checkpoint = torch.load(model_pt, map_location='cpu')
        
        print(f"  Checkpoint type: {type(checkpoint)}")
        print(f"  Keys: {list(checkpoint.keys())[:20]}")
        
        # Search for var_names
        def find_in_checkpoint(obj, path="root", depth=0, max_depth=5):
            if depth > max_depth:
                return []
            
            results = []
            
            if isinstance(obj, dict):
                for k, v in list(obj.items())[:50]:  # Limit iteration
                    if 'var' in str(k).lower():
                        if isinstance(v, (list, np.ndarray, torch.Tensor)):
                            if len(v) > 100:
                                results.append((f"{path}.{k}", v))
                    
                    if depth < 3:  # Only recurse shallow
                        results.extend(find_in_checkpoint(v, f"{path}.{k}", depth + 1))
            
            return results
        
        candidates = find_in_checkpoint(checkpoint)
        
        if candidates:
            print(f"\n  Found {len(candidates)} candidate(s):")
            for path, genes in candidates:
                print(f"\n    Path: {path}")
                print(f"    Type: {type(genes)}")
                print(f"    Length: {len(genes)}")
                
                if hasattr(genes, 'numpy'):
                    genes = genes.numpy()
                genes_list = list(genes)[:5]
                print(f"    First 5: {genes_list}")
        else:
            print(f"  ❌ No var_names found in checkpoint")
        
        return None
        
    except Exception as e:
        print(f"  ❌ Failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# METHOD 5: Compare with h5ad genes
# ============================================================================

def method5_compare_h5ad(model_path, adata, model_type='SCANVI'):
    """
    Load model with different gene subsets and see which works
    This helps identify which genes the model was trained on
    """
    print(f"\n{'='*70}")
    print(f"METHOD 5: Trial loading with gene subsets")
    print(f"{'='*70}")
    
    try:
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
        
        # Try loading with full adata
        print(f"\n  Test 1: Load with full adata ({adata.n_vars} genes)...")
        try:
            model = Model.load(model_path, adata=adata)
            print(f"    ✓ Success!")
            
            # If successful, the model likely trained on all genes
            print(f"\n    Model likely trained on all {adata.n_vars} genes")
            print(f"    Genes: {list(adata.var_names)[:10]}")
            return list(adata.var_names)
            
        except Exception as e:
            print(f"    ❌ Failed: {e}")
        
        # Try with HVG if available
        if 'highly_variable' in adata.var.columns:
            hvg_genes = adata.var_names[adata.var['highly_variable']].tolist()
            print(f"\n  Test 2: Load with HVG subset ({len(hvg_genes)} genes)...")
            
            adata_hvg = adata[:, hvg_genes].copy()
            
            try:
                model = Model.load(model_path, adata=adata_hvg)
                print(f"    ✓ Success!")
                print(f"\n    Model trained on {len(hvg_genes)} HVGs")
                print(f"    Genes: {hvg_genes[:10]}")
                return hvg_genes
                
            except Exception as e:
                print(f"    ❌ Failed: {e}")
        
        print(f"\n  ❌ Cannot determine training genes")
        return None
        
    except Exception as e:
        print(f"  ❌ Failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# MAIN DIAGNOSTIC
# ============================================================================

def main():
    """Run all diagnostic methods"""
    print("="*70)
    print("MODEL GENE EXTRACTION DIAGNOSTIC")
    print("="*70)
    
    # Load h5ad first
    print(f"\nLoading h5ad: {H5AD_PATH}")
    adata = sc.read_h5ad(H5AD_PATH)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    
    if 'highly_variable' in adata.var.columns:
        n_hvg = adata.var['highly_variable'].sum()
        print(f"  HVG marked: {n_hvg:,}")
    else:
        print(f"  No HVG info in var")
    
    # Test each model
    for model_name, model_path in MODELS_TO_TEST.items():
        print(f"\n\n{'#'*70}")
        print(f"TESTING: {model_name}")
        print(f"Path: {model_path}")
        print(f"{'#'*70}")
        
        if not model_path.exists():
            print(f"❌ Model not found: {model_path}")
            continue
        
        model_type = 'SCANVI' if 'scanvi' in str(model_path).lower() else 'SCVI'
        
        # Try all methods
        results = {}
        
        results['method1'] = method1_get_var_names(model_path, model_type)
        results['method2'] = method2_registry_search(model_path, model_type)
        results['method3'] = method3_adata_manager(model_path, adata, model_type)
        results['method4'] = method4_direct_file(model_path)
        results['method5'] = method5_compare_h5ad(model_path, adata, model_type)
        
        # Summary
        print(f"\n{'='*70}")
        print(f"SUMMARY for {model_name}")
        print(f"{'='*70}")
        
        successful_methods = {k: v for k, v in results.items() if v is not None}
        
        if successful_methods:
            print(f"\n✓ {len(successful_methods)} method(s) succeeded:")
            for method, genes in successful_methods.items():
                print(f"\n  {method}:")
                print(f"    Genes: {len(genes)}")
                print(f"    First 10: {genes[:10]}")
                print(f"    Last 10: {genes[-10:]}")
                
                # Check overlap with adata
                overlap = len(set(genes) & set(adata.var_names))
                print(f"    Overlap with h5ad: {overlap}/{len(genes)} ({overlap/len(genes)*100:.1f}%)")
        else:
            print(f"\n❌ All methods failed")
            print(f"\nPOSSIBLE SOLUTIONS:")
            print(f"  1. Model may be corrupted")
            print(f"  2. Model saved with different scvi-tools version")
            print(f"  3. Need to use scvi_model instead of scanvi_model")
            print(f"  4. Train new model with explicit gene saving")
    
    print(f"\n{'='*70}")
    print(f"DIAGNOSTIC COMPLETE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()