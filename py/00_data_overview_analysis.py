#!/usr/bin/env python3
"""
00_data_overview_analysis.py

跨组织单细胞数据概览与质量检查
- 统计每个组织的细胞数和细胞类型分布
- 生成数据质量报告
- 为后续分析提供基础信息

作者：临床-生信团队
日期：2025-11-05
版本：v1.0
"""

import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
warnings.filterwarnings('ignore')

# ==================== 配置部分 ====================

# ========== 输入输出配置 ==========
input_h5ad_path = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated.h5ad"
output_dir = "/home/h2048/data/py/1029/cross_tissue_analysis/00_data_overview"

# ========== 关键列名配置 ==========
tissue_key = "tissue_sampling_method"  # 组织部位列名
cell_type_key = "cell_type"             # 细胞类型列名
batch_key = "dataset"                   # 批次列名

# ========== 重点分析的细胞类型 ==========
target_cell_types = ["Epithelial", "T", "Myeloid"]

# ========== 可视化配置 ==========
dpi = 300
figure_format = "pdf"

verbose = True

# ==================== 主要功能 ====================

def log_msg(msg):
    """打印日志"""
    if verbose:
        print(msg)


def log_step(step_num, step_name):
    """打印步骤标题"""
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)


def load_and_check_data(h5ad_path):
    """读取数据并进行基本检查"""
    import scanpy as sc
    
    log_step(1, "Loading and Checking Data")
    
    log_msg(f"\n读取文件: {h5ad_path}")
    if not Path(h5ad_path).exists():
        raise FileNotFoundError(f"文件未找到: {h5ad_path}")
    
    adata = sc.read_h5ad(h5ad_path)
    
    log_msg(f"✓ 数据加载成功")
    log_msg(f"   细胞总数: {adata.n_obs:,}")
    log_msg(f"   基因总数: {adata.n_vars:,}")
    
    # 检查必需的列
    required_cols = [tissue_key, cell_type_key, batch_key]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        raise ValueError(f"缺少必需的列: {', '.join(missing_cols)}")
    
    log_msg(f"\n✓ 必需列检查通过")
    log_msg(f"   组织列: {tissue_key}")
    log_msg(f"   细胞类型列: {cell_type_key}")
    log_msg(f"   批次列: {batch_key}")
    
    # 显示所有可用列
    log_msg(f"\n所有可用的metadata列:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        log_msg(f"   - {col}: {n_unique} unique values")
    
    return adata


def analyze_tissue_distribution(adata, output_dir):
    """分析组织分布"""
    log_step(2, "Tissue Distribution Analysis")
    
    # 统计每个组织的细胞数
    tissue_counts = adata.obs[tissue_key].value_counts().sort_index()
    
    log_msg(f"\n组织分布:")
    for tissue, count in tissue_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {tissue}: {count:,} cells ({pct:.1f}%)")
    
    log_msg(f"\n总计: {len(tissue_counts)} 个组织部位")
    
    # 保存统计表
    tissue_df = pd.DataFrame({
        'Tissue': tissue_counts.index,
        'Cell_Count': tissue_counts.values,
        'Percentage': tissue_counts.values / adata.n_obs * 100
    })
    
    tissue_csv = output_dir / "tissue_distribution.csv"
    tissue_df.to_csv(tissue_csv, index=False)
    log_msg(f"\n✓ 组织分布已保存: {tissue_csv}")
    
    return tissue_df


def analyze_celltype_distribution(adata, output_dir):
    """分析细胞类型分布"""
    log_step(3, "Cell Type Distribution Analysis")
    
    # 统计每个细胞类型的细胞数
    celltype_counts = adata.obs[cell_type_key].value_counts().sort_index()
    
    log_msg(f"\n细胞类型分布:")
    for celltype, count in celltype_counts.items():
        pct = count / adata.n_obs * 100
        marker = "★" if celltype in target_cell_types else " "
        log_msg(f" {marker} {celltype}: {count:,} cells ({pct:.1f}%)")
    
    log_msg(f"\n总计: {len(celltype_counts)} 种细胞类型")
    log_msg(f"★ = 重点分析的细胞类型")
    
    # 保存统计表
    celltype_df = pd.DataFrame({
        'Cell_Type': celltype_counts.index,
        'Cell_Count': celltype_counts.values,
        'Percentage': celltype_counts.values / adata.n_obs * 100,
        'Is_Target': [ct in target_cell_types for ct in celltype_counts.index]
    })
    
    celltype_csv = output_dir / "celltype_distribution.csv"
    celltype_df.to_csv(celltype_csv, index=False)
    log_msg(f"\n✓ 细胞类型分布已保存: {celltype_csv}")
    
    return celltype_df


def analyze_tissue_celltype_matrix(adata, output_dir):
    """分析组织×细胞类型交叉分布"""
    log_step(4, "Tissue × Cell Type Matrix")
    
    # 创建交叉表
    cross_tab = pd.crosstab(
        adata.obs[tissue_key],
        adata.obs[cell_type_key],
        margins=True,
        margins_name='Total'
    )
    
    log_msg(f"\n组织×细胞类型矩阵:")
    log_msg(f"   行数（组织）: {len(cross_tab) - 1}")
    log_msg(f"   列数（细胞类型）: {len(cross_tab.columns) - 1}")
    
    # 保存完整矩阵
    matrix_csv = output_dir / "tissue_celltype_matrix.csv"
    cross_tab.to_csv(matrix_csv)
    log_msg(f"\n✓ 交叉矩阵已保存: {matrix_csv}")
    
    # 显示重点细胞类型的分布
    log_msg(f"\n重点细胞类型在各组织的分布:")
    for celltype in target_cell_types:
        if celltype in cross_tab.columns:
            log_msg(f"\n  {celltype}:")
            for tissue in cross_tab.index[:-1]:  # 排除Total行
                count = cross_tab.loc[tissue, celltype]
                total = cross_tab.loc[tissue, 'Total']
                pct = count / total * 100 if total > 0 else 0
                log_msg(f"    {tissue}: {count:,} ({pct:.1f}%)")
    
    # 创建百分比矩阵（每个组织内的细胞类型比例）
    cross_tab_pct = cross_tab.iloc[:-1, :-1].div(cross_tab['Total'].iloc[:-1], axis=0) * 100
    
    pct_csv = output_dir / "tissue_celltype_percentage.csv"
    cross_tab_pct.to_csv(pct_csv)
    log_msg(f"\n✓ 百分比矩阵已保存: {pct_csv}")
    
    return cross_tab, cross_tab_pct


def analyze_batch_distribution(adata, output_dir):
    """分析批次分布"""
    log_step(5, "Batch Distribution Analysis")
    
    # 批次总体分布
    batch_counts = adata.obs[batch_key].value_counts().sort_index()
    
    log_msg(f"\n批次分布:")
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")
    
    # 批次×组织交叉表
    batch_tissue_tab = pd.crosstab(
        adata.obs[batch_key],
        adata.obs[tissue_key]
    )
    
    log_msg(f"\n批次×组织分布:")
    log_msg(batch_tissue_tab.to_string())
    
    # 保存
    batch_csv = output_dir / "batch_distribution.csv"
    batch_tissue_tab.to_csv(batch_csv)
    log_msg(f"\n✓ 批次分布已保存: {batch_csv}")
    
    return batch_tissue_tab


def check_data_quality(adata, output_dir):
    """检查数据质量指标"""
    log_step(6, "Data Quality Check")
    
    quality_metrics = {}
    
    # 1. 检查是否有缺失值
    log_msg("\n1. 检查缺失值:")
    for col in [tissue_key, cell_type_key, batch_key]:
        n_missing = adata.obs[col].isna().sum()
        pct_missing = n_missing / adata.n_obs * 100
        log_msg(f"   {col}: {n_missing} ({pct_missing:.2f}%)")
        quality_metrics[f'{col}_missing'] = n_missing
    
    # 2. 检查每个组织-细胞类型组合的细胞数
    log_msg("\n2. 检查样本量充足性:")
    log_msg("   组织-细胞类型组合的最小细胞数要求: 100")
    
    insufficient_groups = []
    
    for tissue in adata.obs[tissue_key].unique():
        for celltype in target_cell_types:
            mask = (adata.obs[tissue_key] == tissue) & (adata.obs[cell_type_key] == celltype)
            n_cells = mask.sum()
            
            if n_cells < 100:
                insufficient_groups.append({
                    'Tissue': tissue,
                    'Cell_Type': celltype,
                    'Cell_Count': n_cells,
                    'Status': 'Insufficient' if n_cells < 100 else 'OK'
                })
                log_msg(f"   ⚠️  {tissue} - {celltype}: {n_cells} cells (< 100)")
            else:
                log_msg(f"   ✓  {tissue} - {celltype}: {n_cells} cells")
    
    if insufficient_groups:
        log_msg(f"\n   ⚠️  发现 {len(insufficient_groups)} 个样本量不足的组合")
        log_msg(f"   建议：这些组合的分析结果需谨慎解读")
        
        insuf_df = pd.DataFrame(insufficient_groups)
        insuf_csv = output_dir / "insufficient_sample_sizes.csv"
        insuf_df.to_csv(insuf_csv, index=False)
        log_msg(f"   详细信息已保存: {insuf_csv}")
    else:
        log_msg(f"\n   ✓ 所有组合的样本量充足")
    
    # 3. 检查embeddings
    log_msg("\n3. 检查降维结果:")
    for emb in ['X_pca', 'X_umap']:
        if emb in adata.obsm:
            log_msg(f"   ✓ {emb}: {adata.obsm[emb].shape}")
        else:
            log_msg(f"   ✗ {emb}: 未找到")
    
    # 4. 检查聚类结果
    log_msg("\n4. 检查聚类结果:")
    leiden_cols = [col for col in adata.obs.columns if 'leiden' in col.lower()]
    if leiden_cols:
        for col in leiden_cols:
            n_clusters = adata.obs[col].nunique()
            log_msg(f"   ✓ {col}: {n_clusters} clusters")
    else:
        log_msg(f"   ✗ 未找到leiden聚类结果")
    
    # 保存质量报告
    quality_report = []
    quality_report.append("="*70)
    quality_report.append("Data Quality Report")
    quality_report.append("="*70)
    quality_report.append(f"\nTotal cells: {adata.n_obs:,}")
    quality_report.append(f"Total genes: {adata.n_vars:,}")
    quality_report.append(f"\nMissing values:")
    for col in [tissue_key, cell_type_key, batch_key]:
        n_missing = quality_metrics.get(f'{col}_missing', 0)
        quality_report.append(f"  {col}: {n_missing}")
    quality_report.append(f"\nEmbeddings available:")
    for emb in ['X_pca', 'X_umap']:
        status = "Yes" if emb in adata.obsm else "No"
        quality_report.append(f"  {emb}: {status}")
    quality_report.append(f"\nClustering results:")
    if leiden_cols:
        for col in leiden_cols:
            n_clusters = adata.obs[col].nunique()
            quality_report.append(f"  {col}: {n_clusters} clusters")
    else:
        quality_report.append(f"  None found")
    
    if insufficient_groups:
        quality_report.append(f"\n⚠️  Warning: {len(insufficient_groups)} tissue-celltype combinations with <100 cells")
    else:
        quality_report.append(f"\n✓ All tissue-celltype combinations have sufficient cells")
    
    quality_report.append("\n" + "="*70)
    
    report_text = '\n'.join(quality_report)
    
    report_path = output_dir / "data_quality_report.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    log_msg(f"\n✓ 质量报告已保存: {report_path}")
    
    return quality_metrics


def generate_summary_statistics(adata, output_dir):
    """生成汇总统计"""
    log_step(7, "Summary Statistics")
    
    summary = {
        'Total_Cells': adata.n_obs,
        'Total_Genes': adata.n_vars,
        'N_Tissues': adata.obs[tissue_key].nunique(),
        'N_Cell_Types': adata.obs[cell_type_key].nunique(),
        'N_Batches': adata.obs[batch_key].nunique(),
        'Tissues': ', '.join(sorted(adata.obs[tissue_key].unique())),
        'Target_Cell_Types': ', '.join(target_cell_types),
    }
    
    # 添加重点细胞类型的细胞数
    for celltype in target_cell_types:
        n_cells = (adata.obs[cell_type_key] == celltype).sum()
        summary[f'{celltype}_Cells'] = n_cells
    
    # 转换为DataFrame
    summary_df = pd.DataFrame([summary]).T
    summary_df.columns = ['Value']
    
    log_msg("\n汇总统计:")
    log_msg(summary_df.to_string())
    
    # 保存
    summary_csv = output_dir / "summary_statistics.csv"
    summary_df.to_csv(summary_csv)
    log_msg(f"\n✓ 汇总统计已保存: {summary_csv}")
    
    return summary_df


def generate_overview_report(adata, output_dir):
    """生成数据概览总报告"""
    log_step(8, "Generating Overview Report")
    
    from datetime import datetime
    
    report = []
    report.append("="*70)
    report.append("Cross-Tissue Single Cell Data Overview")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    # 基本信息
    report.append("[ Basic Information ]")
    report.append(f"  Total cells: {adata.n_obs:,}")
    report.append(f"  Total genes: {adata.n_vars:,}")
    report.append("")
    
    # 组织分布
    report.append("[ Tissue Distribution ]")
    tissue_counts = adata.obs[tissue_key].value_counts().sort_index()
    for tissue, count in tissue_counts.items():
        pct = count / adata.n_obs * 100
        report.append(f"  {tissue}: {count:,} cells ({pct:.1f}%)")
    report.append("")
    
    # 细胞类型分布
    report.append("[ Cell Type Distribution ]")
    celltype_counts = adata.obs[cell_type_key].value_counts().sort_index()
    for celltype, count in celltype_counts.items():
        pct = count / adata.n_obs * 100
        marker = "★" if celltype in target_cell_types else " "
        report.append(f" {marker} {celltype}: {count:,} cells ({pct:.1f}%)")
    report.append("")
    report.append("  ★ = Target cell types for analysis")
    report.append("")
    
    # 重点细胞类型统计
    report.append("[ Target Cell Types by Tissue ]")
    for celltype in target_cell_types:
        report.append(f"\n  {celltype}:")
        for tissue in sorted(adata.obs[tissue_key].unique()):
            mask = (adata.obs[tissue_key] == tissue) & (adata.obs[cell_type_key] == celltype)
            n_cells = mask.sum()
            total = (adata.obs[tissue_key] == tissue).sum()
            pct = n_cells / total * 100 if total > 0 else 0
            report.append(f"    {tissue}: {n_cells:,} cells ({pct:.1f}% of tissue)")
    report.append("")
    
    # 批次信息
    report.append("[ Batch Information ]")
    batch_counts = adata.obs[batch_key].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
    report.append("")
    
    # 输出文件
    report.append("[ Output Files ]")
    report.append(f"  - Tissue distribution: {output_dir / 'tissue_distribution.csv'}")
    report.append(f"  - Cell type distribution: {output_dir / 'celltype_distribution.csv'}")
    report.append(f"  - Tissue×CellType matrix: {output_dir / 'tissue_celltype_matrix.csv'}")
    report.append(f"  - Batch distribution: {output_dir / 'batch_distribution.csv'}")
    report.append(f"  - Summary statistics: {output_dir / 'summary_statistics.csv'}")
    report.append(f"  - Quality report: {output_dir / 'data_quality_report.txt'}")
    report.append("")
    
    report.append("="*70)
    report.append("Data overview completed successfully")
    report.append("="*70)
    
    report_text = '\n'.join(report)
    
    # 保存
    report_path = output_dir / "data_overview_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    log_msg(f"\n总报告已保存: {report_path}")
    log_msg("\n" + report_text)


# ==================== 主程序 ====================

def main():
    """主函数"""
    
    print("\n" + "="*70)
    print("跨组织单细胞数据概览分析")
    print("="*70)
    
    start_time = time.time()
    
    # 环境检查
    log_msg("\n检查环境...")
    try:
        import scanpy as sc
        log_msg(f"   scanpy: {sc.__version__}")
    except ImportError as e:
        print(f"\n错误: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 创建输出目录
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    log_msg(f"\n输出目录: {output_dir_path}")
    
    # 主要分析流程
    try:
        # Step 1: 加载数据
        adata = load_and_check_data(input_h5ad_path)
        
        # Step 2: 组织分布
        tissue_df = analyze_tissue_distribution(adata, output_dir_path)
        
        # Step 3: 细胞类型分布
        celltype_df = analyze_celltype_distribution(adata, output_dir_path)
        
        # Step 4: 交叉分布
        cross_tab, cross_tab_pct = analyze_tissue_celltype_matrix(adata, output_dir_path)
        
        # Step 5: 批次分布
        batch_tab = analyze_batch_distribution(adata, output_dir_path)
        
        # Step 6: 质量检查
        quality_metrics = check_data_quality(adata, output_dir_path)
        
        # Step 7: 汇总统计
        summary_df = generate_summary_statistics(adata, output_dir_path)
        
        # Step 8: 总报告
        generate_overview_report(adata, output_dir_path)
        
    except Exception as e:
        print(f"\n❌ 分析过程中出错: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # 完成
    total_time = time.time() - start_time
    
    print("\n" + "="*70)
    print("✅ 数据概览分析完成")
    print("="*70)
    print(f"\n输出目录: {output_dir_path}")
    print(f"运行时间: {total_time:.1f} 秒")
    print(f"\n请查看以下文件:")
    print(f"  - 总报告: data_overview_summary.txt")
    print(f"  - 质量报告: data_quality_report.txt")
    print(f"  - 统计表格: *.csv")
    print()
    
    return adata


if __name__ == "__main__":
    try:
        adata = main()
    except KeyboardInterrupt:
        print("\n\n⚠️  用户中断", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n❌ 错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)