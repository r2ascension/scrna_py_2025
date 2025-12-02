import scanpy as sc
import pandas as pd
import numpy as np

print("="*70)
print("验证当前annotation结果")
print("="*70)

# 1. 读取两次生成的数据
print("\n[1/3] 读取数据...")

# pre_annotation生成的数据（如果存在）
pre_file = "/home/h2048/data/py/1029/rare_cells_analysis/rare_cells_with_new_clusters.h5ad"
# post_annotation生成的数据
post_file = "/home/h2048/data/py/1029/rare_cells_analysis/rare_cells_manually_annotated.h5ad"

from pathlib import Path

if Path(pre_file).exists() and Path(post_file).exists():
    adata_pre = sc.read_h5ad(pre_file)
    adata_post = sc.read_h5ad(post_file)
    
    # 检查聚类是否一致
    cluster_key = "leiden_rare_res2.0"
    
    # 比较每个细胞的cluster编号
    cluster_diff = (adata_pre.obs[cluster_key] != adata_post.obs[cluster_key]).sum()
    
    print(f"✓ 两次数据加载成功")
    print(f"\n聚类一致性检查:")
    print(f"  Pre-annotation聚类: {adata_pre.obs[cluster_key].nunique()} clusters")
    print(f"  Post-annotation聚类: {adata_post.obs[cluster_key].nunique()} clusters")
    print(f"  Cluster编号不一致的细胞: {cluster_diff}/{len(adata_pre)}")
    
    if cluster_diff == 0:
        print(f"\n✓ 聚类结果完全一致 - cell_level_annotation.csv 是正确的！")
    else:
        print(f"\n⚠️ 聚类结果不一致 - 需要重新从头开始！")
        print(f"   建议：删除所有中间文件，按照下面的完整流程重新运行")
else:
    print("未找到pre_annotation数据文件，跳过一致性检查")
    print("假设cell_level_annotation.csv是基于最后一次聚类结果的")

# 2. 验证cell_level_annotation.csv的内容
print("\n[2/3] 验证annotation表...")

annotation_file = "/home/h2048/data/py/1029/rare_cells_analysis/cell_level_annotation.csv"
cell_annotations = pd.read_csv(annotation_file)

changed = cell_annotations[cell_annotations['Annotation_Changed'] == True]
print(f"  修改的细胞: {len(changed)}/{len(cell_annotations)}")
print(f"\n  修改统计:")
change_summary = changed.groupby(['Original_CellType', 'Manual_CellType']).size()
for (orig, new), count in change_summary.items():
    print(f"    {orig} → {new}: {count} cells")

# 3. 如果已经更新了主数据集，验证结果
print("\n[3/3] 验证主数据集更新...")

main_file = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"

if Path(main_file).exists():
    adata_main = sc.read_h5ad(main_file)
    
    # 检查前10个应该更新的细胞
    print(f"\n  验证前10个修改的细胞:")
    verified = 0
    failed = 0
    
    for _, row in changed.head(10).iterrows():
        cell_id = row['Cell_ID']
        expected = row['Manual_CellType']
        
        if cell_id in adata_main.obs.index:
            actual = adata_main.obs.loc[cell_id, 'cell_type']
            if actual == expected:
                print(f"    ✓ {cell_id}: {expected}")
                verified += 1
            else:
                print(f"    ✗ {cell_id}: 期望={expected}, 实际={actual}")
                failed += 1
    
    if failed == 0:
        print(f"\n✓ 验证通过！所有更新都正确")
        print(f"✓ 可以安全使用 {main_file}")
    else:
        print(f"\n✗ 验证失败！发现 {failed} 个错误")
        print(f"✗ 建议：按照下面的完整流程重新运行")
else:
    print(f"  主数据集尚未更新")

print("\n" + "="*70)