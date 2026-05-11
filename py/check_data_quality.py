#!/usr/bin/env python3
"""
Quick Data Quality Check for cNMF Pipeline
===========================================
Checks for problematic gene/cell names that could cause TSV parsing errors.

Usage:
    python check_data_quality.py /path/to/your_data.h5ad
"""

import sys
import scanpy as sc
import numpy as np

def check_for_special_chars(names, name_type="names"):
    """Check for tabs, newlines, and other problematic characters."""
    names_str = [str(n) for n in names]
    
    issues = {
        'tabs': [i for i, n in enumerate(names_str) if '\t' in n],
        'newlines': [i for i, n in enumerate(names_str) if '\n' in n or '\r' in n],
        'unprintable': [i for i, n in enumerate(names_str) if any(not c.isprintable() and c not in [' ', '\t', '\n', '\r'] for c in n)]
    }
    
    total_issues = sum(len(v) for v in issues.values())
    
    print(f"\n{name_type.upper()} Quality Check:")
    print(f"  Total {name_type}: {len(names):,}")
    print(f"  With tabs (\\t): {len(issues['tabs']):,}")
    print(f"  With newlines (\\n): {len(issues['newlines']):,}")
    print(f"  With unprintable chars: {len(issues['unprintable']):,}")
    print(f"  Total issues: {total_issues:,}")
    
    if total_issues > 0:
        print(f"\n  ⚠️  FOUND {total_issues} PROBLEMATIC {name_type.upper()}")
        print(f"  Examples (first 5):")
        shown = 0
        for issue_type, indices in issues.items():
            for idx in indices[:5-shown]:
                name = names_str[idx]
                repr_name = repr(name)  # Shows special chars
                print(f"    - [{idx}] {repr_name}")
                shown += 1
                if shown >= 5:
                    break
            if shown >= 5:
                break
        return True
    else:
        print(f"  ✓ All {name_type} are clean!")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python check_data_quality.py /path/to/data.h5ad")
        sys.exit(1)
    
    h5ad_path = sys.argv[1]
    
    print("=" * 70)
    print("cNMF Data Quality Check")
    print("=" * 70)
    print(f"File: {h5ad_path}")
    
    try:
        print("\nLoading data...")
        adata = sc.read_h5ad(h5ad_path)
        print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        print(f"  Layers: {list(adata.layers.keys())}")
        print(f"  Has .raw: {'✓' if adata.raw is not None else '✗'}")
        
        # Check gene names
        has_gene_issues = check_for_special_chars(adata.var_names, "gene names")
        
        # Check cell names
        has_cell_issues = check_for_special_chars(adata.obs_names, "cell names/barcodes")
        
        # Summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        
        if has_gene_issues or has_cell_issues:
            print("⚠️  ISSUES FOUND!")
            print("\nRecommendation:")
            print("  Use v1.3-HOTFIX2 pipeline which automatically sanitizes names.")
            print("  The pipeline will replace tabs/newlines with underscores.")
            print("\nOr manually fix with:")
            print("  adata.var_names = [str(g).replace('\\t','_').replace('\\n','_') for g in adata.var_names]")
            print("  adata.obs_names = [str(c).replace('\\t','_').replace('\\n','_') for c in adata.obs_names]")
        else:
            print("✅ No issues found!")
            print("\nYour data is clean and ready for cNMF analysis.")
        
        print("=" * 70)
        
    except Exception as e:
        print(f"\n❌ Error loading file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
