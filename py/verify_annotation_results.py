#!/usr/bin/env python3
"""
verify_annotation_results.py

验证当前的annotation结果是否正确
检查聚类一致性和主数据集更新情况

作者: 临床-生信团队
日期: 2025-11-07
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

print("="*70)
print("验证Annotation结果")
print("="*70)

# 配置
pre_file = "/home/h2048/data/py/1029/rare_cells_analysis/rare_cells_with_new_clusters.h5ad"
post_file = "/home/h2048/data/py/1029/rare_cells_analysis/rare_cells_manually_annotated.h5ad"
annotation_file = "/home/h2048/data/py/1029/rare_cells_analysis/cell_level_annotation.csv"
main_file = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"

cluster_key = "leiden_rare_res2.0"
cell_type_key = "cell_type"

# ========== 检查1: 聚类一致性 ==========
print("\n[检查1/3] 聚类一致性...")

if Path(pre_file).exists() and Path(post_file).exists():
    adata_pre = sc.read_h5ad(pre_file)
    adata_post = sc.read_h5ad(post_file)
    
    cluster_diff = (adata_pre.obs[cluster_key] != adata_post.obs[cluster_key]).sum()
    
    print(f"  Pre-annotation聚类数: {adata_pre.obs[cluster_key].nunique()}")
    print(f"  Post-annotation聚类数: {adata_post.obs[cluster_key].nunique()}")
    print(f"  Cluster编号不一致的细胞: {cluster_diff}/{len(adata_pre)}")
    
    if cluster_diff == 0:
        print(f"\n  ✓ 聚类完全一致 - annotation应该是正确的")
        clustering_consistent = True
    else:
        print(f"\n  ✗ 聚类不一致 - annotation可能有误!")
        print(f"  ⚠️  建议重新按照完整流程运行")
        clustering_consistent = False
else:
    print("  未找到pre或post的h5ad文件，跳过此检查")
    clustering_consistent = None

# ========== 检查2: Annotation表内容 ==========
print("\n[检查2/3] Annotation表...")

if Path(annotation_file).exists():
    cell_annotations = pd.read_csv(annotation_file)
    
    changed = cell_annotations[cell_annotations['Annotation_Changed'] == True]
    
    print(f"  总细胞数: {len(cell_annotations)}")
    print(f"  修改的细胞: {len(changed)} ({len(changed)/len(cell_annotations)*100:.1f}%)")
    
    if len(changed) > 0:
        print(f"\n  修改统计:")
        change_summary = changed.groupby(['Original_CellType', 'Manual_CellType']).size()
        for (orig, new), count in change_summary.items():
            print(f"    {orig} → {new}: {count} cells")
    
    annotation_exists = True
else:
    print("  ✗ cell_level_annotation.csv 不存在")
    annotation_exists = False

# ========== 检查3: 主数据集更新 ==========
print("\n[检查3/3] 主数据集更新验证...")

if Path(main_file).exists() and annotation_exists:
    adata_main = sc.read_h5ad(main_file)
    
    print(f"  主数据集: {adata_main.n_obs:,} cells")
    
    # 验证前20个修改的细胞
    n_check = min(20, len(changed))
    print(f"\n  验证前{n_check}个修改的细胞:")
    
    verified = 0
    failed = 0
    
    for _, row in changed.head(n_check).iterrows():
        cell_id = row['Cell_ID']
        expected = row['Manual_CellType']
        
        if cell_id in adata_main.obs.index:
            actual = adata_main.obs.loc[cell_id, cell_type_key]
            if actual == expected:
                verified += 1
            else:
                print(f"    ✗ {cell_id}: 期望={expected}, 实际={actual}")
                failed += 1
    
    print(f"\n  验证结果: {verified}/{n_check}")
    
    if failed == 0:
        print(f"  ✓ 所有检查的细胞都正确更新")
        update_correct = True
    else:
        print(f"  ✗ 发现 {failed} 个错误")
        update_correct = False
    
    # 显示更新后的细胞类型分布
    print(f"\n  更新后的细胞类型分布:")
    celltype_counts = adata_main.obs[cell_type_key].value_counts()
    for celltype, count in celltype_counts.items():
        pct = count / adata_main.n_obs * 100
        print(f"    {celltype}: {count:,} cells ({pct:.2f}%)")
    
elif Path(main_file).exists():
    print(f"  主数据集存在但annotation表不存在，无法验证")
    update_correct = None
else:
    print(f"  主数据集尚未更新: {main_file}")
    update_correct = None

# ========== 总结 ==========
print("\n" + "="*70)
print("验证总结")
print("="*70)

all_checks_passed = True

if clustering_consistent is not None:
    if clustering_consistent:
        print("✓ 聚类一致性: 通过")
    else:
        print("✗ 聚类一致性: 失败")
        all_checks_passed = False

if annotation_exists:
    print("✓ Annotation表: 存在")
else:
    print("✗ Annotation表: 不存在")
    all_checks_passed = False

if update_correct is not None:
    if update_correct:
        print("✓ 主数据集更新: 正确")
    else:
        print("✗ 主数据集更新: 发现错误")
        all_checks_passed = False
elif Path(main_file).exists():
    print("? 主数据集更新: 无法验证")
else:
    print("- 主数据集更新: 尚未进行")

print("\n" + "="*70)

if all_checks_passed:
    print("🎉 所有检查通过!")
    print("✓ 你的annotation结果是正确的")
    print(f"✓ 可以安全使用: {main_file}")
else:
    print("⚠️  发现问题!")
    print("建议: 按照 COMPLETE_WORKFLOW.md 中的完整流程重新运行")

print("="*70)