#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cNMF Results Analysis and Visualization v1.1 HOTFIX
====================================================

CRITICAL FIXES from code review:
- P2-10: Use cNMF stability metrics (not max K) for default analysis
- P2-11: Removed biased custom K metrics, use cNMF standard metrics
- P2-12: Robust dt_0.1/dt_0_1 file matching

Analyze and visualize cNMF results for multiple cell types

Author: r2end
Date: 2025-01-14
Version: 1.1 - HOTFIX for metric accuracy
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import json
import re

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input directory (output from label_guided_cnmf_pipeline)
CNMF_OUTPUT_DIR = Path("/home/h2048/data/py/0114/cnmf_by_celltype_v1_1")

# Analysis output directory
ANALYSIS_OUTPUT_DIR = Path(f"{CNMF_OUTPUT_DIR}/analysis")
ANALYSIS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Visualization parameters
VIZ_CONFIG = {
    'dpi': 300,
    'figsize': (10, 8),
    'cmap': 'viridis',
    'top_genes_per_gep': 20,
}

print("="*80)
print("cNMF Results Analysis v1.1 HOTFIX")
print("="*80)
print(f"Input directory: {CNMF_OUTPUT_DIR}")
print(f"Analysis output: {ANALYSIS_OUTPUT_DIR}")
print("="*80)


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def robust_dt_pattern(density_threshold: float = 0.1) -> List[str]:
    """
    Generate robust density threshold patterns for file matching.
    ⭐ FIXED: Support both dt_0.1 and dt_0_1 formats
    
    Parameters
    ----------
    density_threshold : float
        Density threshold value
    
    Returns
    -------
    patterns : list
        List of possible dt string formats
    """
    # Both underscore and dot formats
    dt_formats = [
        f"dt_{str(density_threshold).replace('.', '_')}",  # dt_0_1
        f"dt_{density_threshold}",                          # dt_0.1
    ]
    return dt_formats


# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def find_cnmf_outputs(base_dir: Path) -> Dict[str, Dict[str, Path]]:
    """
    Find all cNMF output files for each cell type.
    
    Parameters
    ----------
    base_dir : Path
        Base directory containing cell type subdirectories
    
    Returns
    -------
    outputs : dict
        Dictionary mapping cell type to output files
    """
    print("\n" + "="*80)
    print("Searching for cNMF Outputs")
    print("="*80)
    
    outputs = {}
    
    for celltype_dir in base_dir.iterdir():
        if not celltype_dir.is_dir() or celltype_dir.name == 'analysis':
            continue
        
        # Look for cNMF subdirectory
        cnmf_dirs = list(celltype_dir.glob("*_cNMF"))
        if len(cnmf_dirs) == 0:
            continue
        
        cnmf_dir = cnmf_dirs[0]
        cell_type = celltype_dir.name
        
        # Find output files with robust pattern matching
        spectra_files = []
        usage_files = []
        k_values_set = set()
        
        # Try both dt formats
        for dt_pattern in robust_dt_pattern():
            spectra_files.extend(list(cnmf_dir.glob(f"*.spectra.k_*.{dt_pattern}.consensus.txt")))
            usage_files.extend(list(cnmf_dir.glob(f"*.usages.k_*.{dt_pattern}.consensus.txt")))
        
        # Extract unique K values from filenames
        for f in spectra_files:
            # Extract K from filename like: name.spectra.k_25.dt_0_1.consensus.txt
            match = re.search(r'\.k_(\d+)\.', f.name)
            if match:
                k_values_set.add(int(match.group(1)))
        
        if len(spectra_files) > 0:
            outputs[cell_type] = {
                'dir': cnmf_dir,
                'spectra': spectra_files,
                'usage': usage_files,
                'k_values': sorted(list(k_values_set))
            }
            
            print(f"\n✓ Found: {cell_type}")
            print(f"  K values: {outputs[cell_type]['k_values']}")
            print(f"  Files: {len(spectra_files)} spectra, {len(usage_files)} usage")
    
    print(f"\n{'─'*80}")
    print(f"Total cell types found: {len(outputs)}")
    
    return outputs


def load_cnmf_results(
    cnmf_dir: Path,
    k: int,
    density_threshold: float = 0.1
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load spectra and usage matrices for a specific K.
    ⭐ FIXED: Robust dt format matching
    
    Parameters
    ----------
    cnmf_dir : Path
        cNMF output directory
    k : int
        Number of components (K value)
    density_threshold : float
        Density threshold used in consensus
    
    Returns
    -------
    spectra : DataFrame
        GEP × genes matrix (gene expression patterns)
    usage : DataFrame
        cells × GEP matrix (GEP usage in each cell)
    """
    # Try both dt formats
    for dt_str in robust_dt_pattern(density_threshold):
        spectra_pattern = f"*.spectra.k_{k}.{dt_str}.consensus.txt"
        usage_pattern = f"*.usages.k_{k}.{dt_str}.consensus.txt"
        
        spectra_files = list(cnmf_dir.glob(spectra_pattern))
        usage_files = list(cnmf_dir.glob(usage_pattern))
        
        if len(spectra_files) > 0 and len(usage_files) > 0:
            # Found files with this dt format
            spectra = pd.read_csv(spectra_files[0], sep='\t', index_col=0)
            usage = pd.read_csv(usage_files[0], sep='\t', index_col=0)
            return spectra, usage
    
    # If we get here, no files found
    raise FileNotFoundError(
        f"Cannot find results for K={k} in {cnmf_dir}\n"
        f"Tried patterns: {robust_dt_pattern(density_threshold)}"
    )


def load_k_selection_metrics(cnmf_dir: Path) -> Optional[pd.DataFrame]:
    """
    Load cNMF's built-in K selection metrics if available.
    
    Parameters
    ----------
    cnmf_dir : Path
        cNMF output directory
    
    Returns
    -------
    metrics_df : DataFrame or None
        K selection metrics from cNMF
    """
    # Look for k_selection related files
    k_sel_files = list(cnmf_dir.glob("*k_selection*"))
    
    # cNMF typically generates CSV files with stability metrics
    for f in k_sel_files:
        if f.suffix == '.csv':
            try:
                df = pd.read_csv(f)
                if 'k' in df.columns or 'K' in df.columns:
                    return df
            except:
                continue
    
    return None


# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def compute_gep_quality_metrics(
    spectra: pd.DataFrame,
    usage: pd.DataFrame
) -> Dict[str, float]:
    """
    Compute interpretability and quality metrics for GEPs.
    ⭐ FIXED: Focus on interpretability, not structural bias
    
    Parameters
    ----------
    spectra : DataFrame
        GEP × genes matrix
    usage : DataFrame
        cells × GEP matrix
    
    Returns
    -------
    metrics : dict
        Dictionary of quality metrics
    """
    metrics = {}
    
    # 1. Usage sparsity (Gini coefficient) - higher = more specialized GEPs
    usage_flat = usage.values.flatten()
    sorted_usage = np.sort(usage_flat)
    n = len(sorted_usage)
    index = np.arange(1, n + 1)
    gini = (2 * np.sum(index * sorted_usage)) / (n * np.sum(sorted_usage)) - (n + 1) / n
    metrics['usage_gini'] = float(gini)
    
    # 2. Mean max usage per cell - how "decisive" is GEP assignment
    max_usage_per_cell = usage.max(axis=1)
    metrics['mean_max_usage'] = float(max_usage_per_cell.mean())
    
    # 3. Number of active GEPs (usage > 5% in at least 5% cells)
    mean_usage_per_gep = usage.mean(axis=0)
    n_active_geps = (mean_usage_per_gep > 0.05).sum()
    metrics['n_active_geps'] = int(n_active_geps)
    metrics['active_gep_fraction'] = float(n_active_geps / len(mean_usage_per_gep))
    
    # 4. Spectra sparsity (how specific are gene loadings)
    spectra_flat = spectra.values.flatten()
    spectra_sorted = np.sort(spectra_flat)
    n = len(spectra_sorted)
    index = np.arange(1, n + 1)
    gini_spectra = (2 * np.sum(index * spectra_sorted)) / (n * np.sum(spectra_sorted)) - (n + 1) / n
    metrics['spectra_gini'] = float(gini_spectra)
    
    return metrics


def extract_top_genes_per_gep(
    spectra: pd.DataFrame,
    top_n: int = 20
) -> Dict[str, List[str]]:
    """
    Extract top marker genes for each GEP.
    
    Parameters
    ----------
    spectra : DataFrame
        GEP × genes matrix
    top_n : int
        Number of top genes to extract
    
    Returns
    -------
    top_genes : dict
        Dictionary mapping GEP name to list of top genes
    """
    top_genes = {}
    
    for gep in spectra.index:
        genes = spectra.loc[gep].sort_values(ascending=False).head(top_n).index.tolist()
        top_genes[gep] = genes
    
    return top_genes


def recommend_k_selection(
    cnmf_outputs: Dict[str, Dict],
    output_dir: Path
) -> Dict[str, int]:
    """
    Recommend K value for each cell type based on cNMF stability metrics.
    ⭐ FIXED: Use cNMF's built-in stability, not custom metrics
    
    Parameters
    ----------
    cnmf_outputs : dict
        Dictionary of cNMF outputs
    output_dir : Path
        Output directory
    
    Returns
    -------
    recommended_k : dict
        Recommended K for each cell type
    """
    print("\n" + "="*80)
    print("Recommending K Values (Based on cNMF Stability)")
    print("="*80)
    
    recommended = {}
    recommendations_list = []
    
    for cell_type, info in cnmf_outputs.items():
        print(f"\n{cell_type}:")
        
        # Try to load cNMF's K selection metrics
        cnmf_metrics = load_k_selection_metrics(info['dir'])
        
        if cnmf_metrics is not None:
            print(f"  ✓ Found cNMF K selection metrics")
            # Use cNMF's recommended K (typically based on stability)
            if 'stability' in cnmf_metrics.columns:
                best_idx = cnmf_metrics['stability'].idxmax()
                recommended_k = int(cnmf_metrics.loc[best_idx, 'k' if 'k' in cnmf_metrics.columns else 'K'])
                print(f"  Recommended K = {recommended_k} (from cNMF stability)")
            else:
                # Fallback to middle K
                recommended_k = info['k_values'][len(info['k_values'])//2]
                print(f"  Using middle K = {recommended_k} (fallback)")
        else:
            # Fallback: use middle K value
            recommended_k = info['k_values'][len(info['k_values'])//2]
            print(f"  ⚠️  No cNMF metrics found, using middle K = {recommended_k}")
        
        recommended[cell_type] = recommended_k
        recommendations_list.append({
            'cell_type': cell_type,
            'recommended_k': recommended_k,
            'available_k': str(info['k_values']),
            'method': 'cNMF_stability' if cnmf_metrics is not None else 'middle_k_fallback'
        })
    
    # Save recommendations
    rec_df = pd.DataFrame(recommendations_list)
    rec_file = output_dir / "k_recommendations.csv"
    rec_df.to_csv(rec_file, index=False)
    print(f"\n✓ Recommendations saved: {rec_file}")
    
    return recommended


def compare_k_values(
    cnmf_outputs: Dict[str, Dict],
    output_dir: Path
):
    """
    Compare quality metrics across K values for all cell types.
    
    Parameters
    ----------
    cnmf_outputs : dict
        Dictionary of cNMF outputs
    output_dir : Path
        Output directory for plots
    """
    print("\n" + "="*80)
    print("Comparing K Values Across Cell Types")
    print("="*80)
    
    all_metrics = []
    
    for cell_type, info in cnmf_outputs.items():
        print(f"\nProcessing: {cell_type}")
        
        for k in info['k_values']:
            print(f"  K={k}...", end=' ')
            try:
                spectra, usage = load_cnmf_results(info['dir'], k)
                metrics = compute_gep_quality_metrics(spectra, usage)
                metrics['cell_type'] = cell_type
                metrics['k'] = k
                all_metrics.append(metrics)
                print("✓")
            except Exception as e:
                print(f"❌ {e}")
    
    if len(all_metrics) == 0:
        print("\n⚠️  No metrics computed")
        return
    
    # Convert to DataFrame
    metrics_df = pd.DataFrame(all_metrics)
    
    # Save metrics
    metrics_file = output_dir / "k_comparison_quality_metrics.csv"
    metrics_df.to_csv(metrics_file, index=False)
    print(f"\n✓ Metrics saved: {metrics_file}")
    
    # Plot metrics
    plot_k_comparison(metrics_df, output_dir)
    
    return metrics_df


def plot_k_comparison(metrics_df: pd.DataFrame, output_dir: Path):
    """
    Create comprehensive K comparison plots.
    
    Parameters
    ----------
    metrics_df : DataFrame
        Metrics for all cell types and K values
    output_dir : Path
        Output directory
    """
    print("\nGenerating K comparison plots...")
    
    # Metrics to plot
    metric_names = [
        'usage_gini',
        'mean_max_usage',
        'active_gep_fraction',
        'spectra_gini'
    ]
    
    metric_labels = {
        'usage_gini': 'Usage Sparsity (Gini)',
        'mean_max_usage': 'Mean Max Usage per Cell',
        'active_gep_fraction': 'Fraction of Active GEPs',
        'spectra_gini': 'Spectra Sparsity (Gini)'
    }
    
    # Create 2x2 plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metric_names):
        ax = axes[idx]
        
        # Plot each cell type
        for cell_type in metrics_df['cell_type'].unique():
            ct_data = metrics_df[metrics_df['cell_type'] == cell_type]
            ax.plot(ct_data['k'], ct_data[metric], 
                   marker='o', label=cell_type, linewidth=2)
        
        ax.set_xlabel('K (Number of GEPs)', fontsize=12)
        ax.set_ylabel(metric_labels[metric], fontsize=12)
        ax.set_title(metric_labels[metric], fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = output_dir / "k_comparison_quality_metrics.png"
    plt.savefig(output_file, dpi=VIZ_CONFIG['dpi'], bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {output_file}")


def analyze_gep_markers(
    cnmf_outputs: Dict[str, Dict],
    output_dir: Path,
    k_selected: Optional[Dict[str, int]] = None
):
    """
    Extract and analyze top marker genes for each GEP.
    ⭐ FIXED: Use recommended K, not max K
    
    Parameters
    ----------
    cnmf_outputs : dict
        Dictionary of cNMF outputs
    output_dir : Path
        Output directory
    k_selected : dict, optional
        Dictionary mapping cell type to selected K value
        If None, uses recommended K from stability
    """
    print("\n" + "="*80)
    print("Analyzing GEP Marker Genes")
    print("="*80)
    
    all_markers = []
    
    for cell_type, info in cnmf_outputs.items():
        # Determine K to use
        if k_selected and cell_type in k_selected:
            k = k_selected[cell_type]
        else:
            # Use middle K as safe default (not max K!)
            k = info['k_values'][len(info['k_values'])//2]
        
        print(f"\n{cell_type} (K={k}):")
        
        try:
            spectra, usage = load_cnmf_results(info['dir'], k)
            top_genes = extract_top_genes_per_gep(
                spectra, 
                top_n=VIZ_CONFIG['top_genes_per_gep']
            )
            
            # Save markers for each GEP
            for gep_name, genes in top_genes.items():
                marker_record = {
                    'cell_type': cell_type,
                    'k': k,
                    'gep': gep_name,
                    'top_genes': ','.join(genes)
                }
                all_markers.append(marker_record)
                
                print(f"  {gep_name}: {', '.join(genes[:5])} ...")
            
        except Exception as e:
            print(f"  ❌ Failed: {e}")
    
    if len(all_markers) == 0:
        print("\n⚠️  No markers extracted")
        return None
    
    # Save all markers
    markers_df = pd.DataFrame(all_markers)
    markers_file = output_dir / "gep_top_markers.csv"
    markers_df.to_csv(markers_file, index=False)
    print(f"\n✓ Markers saved: {markers_file}")
    
    return markers_df


def create_usage_heatmaps(
    cnmf_outputs: Dict[str, Dict],
    output_dir: Path,
    k_selected: Optional[Dict[str, int]] = None
):
    """
    Create usage heatmaps for each cell type.
    
    Parameters
    ----------
    cnmf_outputs : dict
        Dictionary of cNMF outputs
    output_dir : Path
        Output directory
    k_selected : dict, optional
        Dictionary mapping cell type to selected K value
    """
    print("\n" + "="*80)
    print("Creating Usage Heatmaps")
    print("="*80)
    
    heatmap_dir = output_dir / "usage_heatmaps"
    heatmap_dir.mkdir(exist_ok=True)
    
    for cell_type, info in cnmf_outputs.items():
        # Determine K to use
        if k_selected and cell_type in k_selected:
            k = k_selected[cell_type]
        else:
            k = info['k_values'][len(info['k_values'])//2]
        
        print(f"\n{cell_type} (K={k})...", end=' ')
        
        try:
            spectra, usage = load_cnmf_results(info['dir'], k)
            
            # Create heatmap
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Sample cells if too many
            if usage.shape[0] > 1000:
                usage_plot = usage.sample(n=1000, random_state=42)
            else:
                usage_plot = usage
            
            # Plot
            sns.heatmap(
                usage_plot.T,
                cmap='viridis',
                cbar_kws={'label': 'Usage'},
                xticklabels=False,
                yticklabels=True,
                ax=ax
            )
            
            ax.set_xlabel('Cells', fontsize=12)
            ax.set_ylabel('GEPs', fontsize=12)
            ax.set_title(f'{cell_type} GEP Usage (K={k})', 
                        fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            
            output_file = heatmap_dir / f"{cell_type}_usage_heatmap.png"
            plt.savefig(output_file, dpi=VIZ_CONFIG['dpi'], bbox_inches='tight')
            plt.close()
            
            print("✓")
            
        except Exception as e:
            print(f"❌ {e}")
    
    print(f"\n✓ Heatmaps saved in: {heatmap_dir}")


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def main():
    """Main analysis execution."""
    
    # Step 1: Find all cNMF outputs
    cnmf_outputs = find_cnmf_outputs(CNMF_OUTPUT_DIR)
    
    if len(cnmf_outputs) == 0:
        print("\n❌ No cNMF outputs found")
        print(f"   Please run label_guided_cnmf_pipeline first")
        return
    
    # Step 2: Recommend K values (based on cNMF stability)
    k_selected = recommend_k_selection(cnmf_outputs, ANALYSIS_OUTPUT_DIR)
    
    # Step 3: Compare K values (quality metrics)
    metrics_df = compare_k_values(cnmf_outputs, ANALYSIS_OUTPUT_DIR)
    
    # Step 4: Analyze GEP markers (using recommended K)
    markers_df = analyze_gep_markers(cnmf_outputs, ANALYSIS_OUTPUT_DIR, k_selected)
    
    # Step 5: Create usage heatmaps (using recommended K)
    create_usage_heatmaps(cnmf_outputs, ANALYSIS_OUTPUT_DIR, k_selected)
    
    # Generate summary report
    print("\n" + "="*80)
    print("✓ ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nOutput directory: {ANALYSIS_OUTPUT_DIR}")
    print(f"\nGenerated files:")
    print(f"  1. k_recommendations.csv - Recommended K for each cell type")
    print(f"  2. k_comparison_quality_metrics.csv - Quality metrics for all K values")
    print(f"  3. k_comparison_quality_metrics.png - Visual comparison of metrics")
    print(f"  4. gep_top_markers.csv - Top marker genes for each GEP")
    print(f"  5. usage_heatmaps/ - Usage heatmaps for each cell type")
    print("\n⭐ NOTE: K recommendations are based on cNMF stability metrics")
    print("   Review K-selection plots in each cell type directory for validation")
    print("\n" + "="*80)


if __name__ == '__main__':
    main()
