#!/usr/bin/env python3
"""
apply_annotation_safe.py

安全应用annotation - 基于固定的聚类结果
保证不会重新聚类，避免cluster编号不一致问题

使用场景:
1. 已运行 rare_cells_subset_analysis.py (pre_annotation模式)
2. 已创建 Annotation.csv (基于DotPlot手动注释)
3. 想要安全地应用annotation到主数据集

作者: 临床-生信团队
日期: 2025-11-07
"""

import scanpy as sc
import pandas as pd
from pathlib import Path
from datetime import datetime

print("="*70)
print("安全应用Annotation - 基于固定聚类结果")
print("="*70)
print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ========== 配置 ==========
rare_cells_h5ad = "/home/h2048/data/py/1029/rare_cells_analysis/rare_cells_with_new_clusters.h5ad"
annotation_csv = "/home/h2048/data/py/1029/rare_cells_analysis/Annotation.csv"
main_data_h5ad = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated.h5ad"
output_h5ad = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
output_annotation_csv = "/home/h2048/data/py/1029/rare_cells_analysis/cell_level_annotation_safe.csv"

cluster_key = "leiden_rare_res2.0"
cell_type_key = "cell_type"

# ========== Step 1: 读取稀有细胞（带聚类）==========
print("[1/5] 读取pre_annotation生成的聚类结果...")

if not Path(rare_cells_h5ad).exists():
    print(f"\n❌ 错误: {rare_cells_h5ad} 不存在")
    print("   请先运行:")
    print("   1. 设置 RUN_MODE='pre_annotation'")
    print("   2. python rare_cells_subset_analysis.py")
    exit(1)

adata_rare = sc.read_h5ad(rare_cells_h5ad)

print(f"✓ 加载成功")
print(f"  稀有细胞数: {adata_rare.n_obs}")
print(f"  使用聚类: {cluster_key}")
print(f"  Cluster数: {adata_rare.obs[cluster_key].nunique()}")

# 显示原始细胞类型分布
print(f"\n  原始细胞类型分布:")
orig_counts = adata_rare.obs[cell_type_key].value_counts()
for ct, count in orig_counts.items():
    pct = count / adata_rare.n_obs * 100
    print(f"    {ct}: {count} cells ({pct:.1f}%)")

# ========== Step 2: 读取Annotation ==========
print(f"\n[2/5] 读取Annotation.csv...")

if not Path(annotation_csv).exists():
    print(f"\n❌ 错误: {annotation_csv} 不存在")
    print("   请先:")
    print("   1. 查看 figures/rare_cells_markers_dotplot.pdf")
    print("   2. 编辑 Annotation_template.csv")
    print("   3. 保存为 Annotation.csv (无表头)")
    exit(1)

annotations = pd.read_csv(annotation_csv, header=None, names=['Cluster', 'CellType'])

print(f"✓ 加载成功")
print(f"  注释条目: {len(annotations)}")

# 验证注释完整性
annotated_clusters = set(annotations['Cluster'].astype(str))
actual_clusters = set(adata_rare.obs[cluster_key].astype(str).unique())

missing = actual_clusters - annotated_clusters
extra = annotated_clusters - actual_clusters

if missing:
    print(f"\n⚠️  警告: 以下cluster缺少annotation: {', '.join(missing)}")
    print("   这些细胞将保持原始注释")

if extra:
    print(f"\n⚠️  警告: Annotation包含不存在的cluster: {', '.join(extra)}")
    print("   这些注释将被忽略")

# 显示映射
print(f"\n  Annotation映射:")
for _, row in annotations.iterrows():
    # 统计该cluster的细胞数
    n_cells = (adata_rare.obs[cluster_key].astype(str) == str(row['Cluster'])).sum()
    print(f"    Cluster {row['Cluster']} ({n_cells} cells) -> {row['CellType']}")

# ========== Step 3: 应用到稀有细胞 ==========
print(f"\n[3/5] 应用annotation...")

annotation_dict = dict(zip(annotations['Cluster'].astype(str), annotations['CellType']))
adata_rare.obs['manual_celltype'] = adata_rare.obs[cluster_key].astype(str).map(annotation_dict)

# 处理未匹配的
unmatched = adata_rare.obs['manual_celltype'].isna()
if unmatched.sum() > 0:
    print(f"  {unmatched.sum()} 个细胞未匹配，保持原始注释")
    adata_rare.obs.loc[unmatched, 'manual_celltype'] = adata_rare.obs.loc[unmatched, cell_type_key]

# 统计变化
changed_mask = adata_rare.obs[cell_type_key] != adata_rare.obs['manual_celltype']
n_changed = changed_mask.sum()

print(f"\n✓ 应用完成")
print(f"  需要修改: {n_changed}/{adata_rare.n_obs} ({n_changed/adata_rare.n_obs*100:.1f}%)")

# 显示手动注释后的分布
print(f"\n  手动注释后的细胞类型分布:")
manual_counts = adata_rare.obs['manual_celltype'].value_counts()
for ct, count in manual_counts.items():
    pct = count / adata_rare.n_obs * 100
    orig_count = orig_counts.get(ct, 0)
    if count != orig_count:
        diff = count - orig_count
        marker = f"({diff:+d})" if diff != 0 else ""
        print(f"    {ct}: {count} cells ({pct:.1f}%) {marker}")
    else:
        print(f"    {ct}: {count} cells ({pct:.1f}%)")

if n_changed > 0:
    print(f"\n  修改详情:")
    for orig in adata_rare.obs[cell_type_key].unique():
        for new in adata_rare.obs['manual_celltype'].unique():
            count = ((adata_rare.obs[cell_type_key] == orig) & 
                     (adata_rare.obs['manual_celltype'] == new) & 
                     changed_mask).sum()
            if count > 0:
                print(f"    {orig} → {new}: {count} cells")

# 生成细胞级注释表
cell_annotation_list = []
for cell_id in adata_rare.obs.index:
    cell_obs = adata_rare.obs.loc[cell_id]
    cell_annotation_list.append({
        'Cell_ID': cell_id,
        'Original_CellType': cell_obs[cell_type_key],
        'Manual_CellType': cell_obs['manual_celltype'],
        'Cluster': cell_obs[cluster_key],
        'Annotation_Changed': cell_obs[cell_type_key] != cell_obs['manual_celltype']
    })

cell_annotation_df = pd.DataFrame(cell_annotation_list)
cell_annotation_df.to_csv(output_annotation_csv, index=False)
print(f"\n✓ 细胞级注释表已保存: {output_annotation_csv}")

# ========== Step 4: 读取主数据集 ==========
print(f"\n[4/5] 读取主数据集...")

if not Path(main_data_h5ad).exists():
    print(f"\n❌ 错误: {main_data_h5ad} 不存在")
    exit(1)

adata_main = sc.read_h5ad(main_data_h5ad)

print(f"✓ 加载成功")
print(f"  主数据细胞数: {adata_main.n_obs:,}")
print(f"  主数据基因数: {adata_main.n_vars:,}")

# 显示主数据原始分布
print(f"\n  主数据原始细胞类型分布:")
main_orig_counts = adata_main.obs[cell_type_key].value_counts()
for ct, count in main_orig_counts.items():
    pct = count / adata_main.n_obs * 100
    print(f"    {ct}: {count:,} cells ({pct:.2f}%)")

# ========== Step 5: 更新主数据集 ==========
print(f"\n[5/5] 更新主数据集...")

n_updated = 0
n_not_found = 0

for cell_id in adata_rare.obs.index:
    old_type = adata_rare.obs.loc[cell_id, cell_type_key]
    new_type = adata_rare.obs.loc[cell_id, 'manual_celltype']
    
    if old_type != new_type:
        if cell_id in adata_main.obs.index:
            adata_main.obs.loc[cell_id, cell_type_key] = new_type
            n_updated += 1
        else:
            n_not_found += 1

print(f"\n✓ 更新完成")
print(f"  成功更新: {n_updated} cells")
if n_not_found > 0:
    print(f"  未找到: {n_not_found} cells (可能已被过滤)")

# 显示主数据更新后的分布
print(f"\n  主数据更新后的细胞类型分布:")
main_new_counts = adata_main.obs[cell_type_key].value_counts()
for ct, count in main_new_counts.items():
    pct = count / adata_main.n_obs * 100
    orig_count = main_orig_counts.get(ct, 0)
    if count != orig_count:
        diff = count - orig_count
        marker = f"({diff:+,d})"
        print(f"    {ct}: {count:,} cells ({pct:.2f}%) {marker}")
    else:
        print(f"    {ct}: {count:,} cells ({pct:.2f}%)")

# 保存
print(f"\n保存更新后的主数据集...")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata_main.write_h5ad(output_h5ad, compression='gzip')

file_size = Path(output_h5ad).stat().st_size / (1024**3)
print(f"✓ 保存成功: {output_h5ad} ({file_size:.2f} GB)")

# ========== 总结 ==========
print("\n" + "="*70)
print("✓ 完成!")
print("="*70)

print(f"\n关键输出文件:")
print(f"  1. 细胞级注释表: {output_annotation_csv}")
print(f"  2. 更新的主数据集: {output_h5ad}")

print(f"\n修改统计:")
print(f"  - 稀有细胞总数: {adata_rare.n_obs}")
print(f"  - 注释修改: {n_changed} cells")
print(f"  - 主数据集更新: {n_updated} cells")

print(f"\n下一步:")
print(f"  运行 verify_annotation_results.py 验证结果")

print("="*70)