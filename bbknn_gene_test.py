#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
原始Counts层面的基因表达检测（统一修订版）
- 优先使用 layers['counts'] 作为原始 counts，并进行最大值体检提示
- 支持 var_names 去重
- 支持大小写/别名列的稳健基因匹配（如 var['gene_symbols'] / 'symbol' / 'SYMBOL' 等）
- 打印 Dotplot 缺失基因在“大群/上皮”两层面的表达统计
- 输出 CSV 至指定目录
"""

import os
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse

# ==================== 配置 ====================

MAIN_H5AD_PATH = "/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1029/bbknn_celltype_analysis/Epithelial/output_optimized"
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "gene_expression_analysis.csv")

# 判定“原始性”的快速阈值：最大值 < 20 则提示可能不是原始counts
RAW_MAX_HINT_THRESHOLD = 20

# 目标亚群列/取值（若找不到该列，则退化为全体细胞）
CELLTYPE_COL_CANDIDATES = ["cell_type", "celltype", "CellType", "celltype_major"]
EPITHELIAL_LABEL_CANDIDATES = ["Epithelial", "epithelial", "上皮", "EPI", "epi"]

# Dotplot缺失 marker
missing_genes = ['ASCL2', 'ASRGL1', 'KRT19', 'KRT8', 'MYH11',
                 'NOTCH3', 'POU2F3', 'SERPINB3', 'SFTPC', 'SOX9', 'SPDEF']

# 扩展上皮 marker
extended_markers = {
    '基底细胞': ['KRT5', 'KRT14', 'TP63', 'KRT15'],
    '管腔上皮': ['KRT8', 'KRT18', 'KRT19', 'EPCAM'],
    '分泌细胞': ['SCGB1A1', 'MUC5AC', 'MUC5B', 'SPDEF'],
    '纤毛细胞': ['FOXJ1', 'RSPH1', 'PIFO'],
    '杯状细胞': ['MUC5AC', 'MUC5B', 'TFF3'],
    '肺泡细胞': ['SFTPC', 'SFTPB', 'AGER', 'PDPN'],
    '簇细胞': ['POU2F3', 'ASRGL1'],
    '神经内分泌': ['ASCL2', 'CHGA', 'SYP'],
    '祖细胞': ['SOX9', 'SOX2', 'TP63'],
    '其他重要': ['NOTCH3', 'SERPINB3', 'MYH11']
}

# ==================== 工具函数 ====================

def header(title: str):
    print("=" * 80)
    print(title)
    print("=" * 80)

def list_layers(adata):
    print("\nLayers:")
    if len(adata.layers.keys()) == 0:
        print("  (none)")
    for layer_name in adata.layers.keys():
        print(f"  - {layer_name}")

def max_of_matrix(mat):
    try:
        if issparse(mat):
            mv = mat.max()
            return float(mv) if np.isscalar(mv) else float(mv.max())
        else:
            return float(np.max(mat))
    except Exception:
        try:
            return float(mat.max())
        except Exception:
            return None

def choose_counts_layer(adata):
    """优先使用 layers['counts']；同时打印最大值体检提示。"""
    if 'counts' in adata.layers:
        print("\n✓ 使用 layers['counts'] 作为原始counts")
        vmax = max_of_matrix(adata.layers['counts'])
        if vmax is not None and vmax < RAW_MAX_HINT_THRESHOLD:
            print(f"⚠️ layers['counts'] 最大值≈{vmax:.2f}，疑似已归一化/对数化，可能不是原始整数counts，请确认。")
        return 'counts'
    else:
        print("\n⚠️ 未找到 layers['counts']，使用 adata.X")
        vmax = max_of_matrix(adata.X)
        if vmax is not None and vmax < RAW_MAX_HINT_THRESHOLD:
            print(f"⚠️ adata.X 最大值≈{vmax:.2f}，疑似已预处理（归一化/对数化），可能不是原始counts。")
        return None

def find_epi_mask(adata):
    """寻找上皮亚群掩码；找不到时返回全体 True 掩码并提示。"""
    for col in CELLTYPE_COL_CANDIDATES:
        if col in adata.obs.columns:
            col_values = adata.obs[col].astype(str).fillna("")
            for lab in EPITHELIAL_LABEL_CANDIDATES:
                mask = (col_values == lab)
                if mask.any():
                    print(f"\n提取上皮细胞... (obs['{col}'] == '{lab}')")
                    return mask
            # 若列存在但没有匹配标签，尝试包含匹配（如 “Epithelial - basal”）
            contains_mask = col_values.str.contains("|".join(EPITHELIAL_LABEL_CANDIDATES), case=False, regex=True)
            if contains_mask.any():
                print(f"\n提取上皮细胞... (obs['{col}'] 包含 'Epithelial/上皮')")
                return contains_mask
    print("\n⚠️ 未找到明确的上皮细胞标注列/标签，将使用全体细胞作为“上皮”以继续检测。")
    return pd.Series([True] * adata.n_obs, index=adata.obs_names)

def make_gene_lookup(adata):
    """
    构建稳健基因查询：
    1) 直接 var_names 精确匹配
    2) var 中可选别名列（gene_symbols/symbol/SYMBOL/Gene/gene_name/gene）
    3) 忽略大小写匹配
    返回 lookup(gene)->(var_name, idx) 或 (None, None)
    """
    names = adata.var_names.to_list()
    name_to_idx = {g: i for i, g in enumerate(names)}
    lower_to_name = {g.lower(): g for g in names}

    alt_cols = [c for c in ['gene_symbols', 'symbol', 'SYMBOL', 'Gene', 'gene_name', 'gene']
                if c in adata.var.columns]
    alt_map = {}
    for c in alt_cols:
        ser = adata.var[c].astype(str)
        for i, sym in enumerate(ser):
            if sym and sym != 'nan':
                # 多重大小写键 → 首次出现的主名
                alt_map.setdefault(sym, names[i])
                alt_map.setdefault(sym.upper(), names[i])
                alt_map.setdefault(sym.lower(), names[i])

    def lookup(query: str):
        if query in name_to_idx:
            return query, name_to_idx[query]
        if query in alt_map:
            nm = alt_map[query]
            return nm, name_to_idx[nm]
        ql = query.lower()
        if ql in lower_to_name:
            nm = lower_to_name[ql]
            return nm, name_to_idx[nm]
        if query.upper() in alt_map:
            nm = alt_map[query.upper()]
            return nm, name_to_idx[nm]
        if ql in alt_map:
            nm = alt_map[ql]
            return nm, name_to_idx[nm]
        return None, None

    return lookup

def extract_gene_vector(adata, idx, use_layer=None):
    """按列索引抽取表达向量，返回 dense 1D numpy array。"""
    if use_layer and use_layer in adata.layers:
        col = adata.layers[use_layer][:, idx]
    else:
        col = adata.X[:, idx]
    if issparse(col):
        return col.toarray().ravel()
    # 可能是 np.matrix
    arr = np.asarray(col).ravel()
    return arr

def analyze_gene_expression(adata, gene, lookup, use_layer=None, dataset_name="Dataset"):
    """分析单基因表达统计（>0 计为表达）。"""
    var_name, idx = lookup(gene)
    if var_name is None:
        return {
            'Gene': gene,
            'VarName': None,
            'Dataset': dataset_name,
            'Status': '不存在',
            'N_Cells_Expr': 0,
            'Pct_Cells': 0.0,
            'Mean_Counts': 0.0,
            'Median_Counts': 0.0,
            'Max_Counts': 0.0,
            'Total_Counts': 0.0
        }

    vec = extract_gene_vector(adata, idx, use_layer)
    n = vec.shape[0]
    n_cells_expr = int((vec > 0).sum())
    pct_cells = (n_cells_expr / n * 100.0) if n > 0 else 0.0

    pos = vec[vec > 0]
    mean_counts = float(pos.mean()) if pos.size > 0 else 0.0
    median_counts = float(np.median(pos)) if pos.size > 0 else 0.0
    max_counts = float(vec.max()) if vec.size > 0 else 0.0
    total_counts = float(vec.sum()) if vec.size > 0 else 0.0

    return {
        'Gene': gene,
        'VarName': var_name,
        'Dataset': dataset_name,
        'Status': '存在',
        'N_Cells_Expr': n_cells_expr,
        'Pct_Cells': pct_cells,
        'Mean_Counts': mean_counts,
        'Median_Counts': median_counts,
        'Max_Counts': max_counts,
        'Total_Counts': total_counts
    }

# ==================== 主流程 ====================

def main():
    header("原始Counts层面的基因表达检测（统一修订版）")

    # 合并所有基因
    all_markers = set(missing_genes)
    for genes in extended_markers.values():
        all_markers.update(genes)
    all_markers = sorted(all_markers)

    print(f"\n待检测基因总数: {len(all_markers)}")
    print(f"  - Dotplot缺失: {len(missing_genes)}")
    print(f"  - 扩展marker: {len(all_markers) - len(missing_genes)}")

    # 读取数据
    print("\n" + "=" * 80)
    print("读取大群数据...")
    print(f"路径: {MAIN_H5AD_PATH}")
    if not os.path.exists(MAIN_H5AD_PATH):
        raise FileNotFoundError(f"File not found: {MAIN_H5AD_PATH}")

    adata_main = ad.read_h5ad(MAIN_H5AD_PATH)
    # var_names 去重以避免极端重复导致的定位问题
    adata_main.var_names_make_unique()

    print(f"维度: {adata_main.n_obs:,} cells × {adata_main.n_vars:,} genes")
    list_layers(adata_main)

    # 选择原始 counts 层
    use_layer = choose_counts_layer(adata_main)

    # 提取上皮细胞
    mask_epi = find_epi_mask(adata_main)
    adata_epi = adata_main[mask_epi].copy()
    print(f"上皮细胞: {adata_epi.n_obs:,} ({adata_epi.n_obs / adata_main.n_obs * 100:.1f}%)")

    # 基因匹配器（主数据集/上皮子集各建一个，防止列映射差异）
    lookup_main = make_gene_lookup(adata_main)
    lookup_epi = make_gene_lookup(adata_epi)

    # 分析
    header("分析基因表达情况 (使用原始counts)")
    results = []
    for gene in all_markers:
        results.append(analyze_gene_expression(adata_main, gene, lookup_main, use_layer, "大群"))
        results.append(analyze_gene_expression(adata_epi, gene, lookup_epi, use_layer, "上皮"))

    df = pd.DataFrame(results)

    # 重点：Dotplot缺失的基因
    print("\n" + "=" * 80)
    print("重点关注: Dotplot缺失的基因")
    print("=" * 80)

    for gene in missing_genes:
        print(f"\n基因: {gene}")
        print("-" * 60)
        gene_data = df[df['Gene'] == gene]
        for _, row in gene_data.iterrows():
            dataset = row['Dataset']
            status = row['Status']
            if status == '不存在':
                print(f"  {dataset:8s}: ❌ 不存在")
            else:
                n_cells = row['N_Cells_Expr']
                pct = row['Pct_Cells']
                mean = row['Mean_Counts']
                total = row['Total_Counts']
                will_filter = ""
                if dataset == "上皮":
                    if n_cells < 3:
                        will_filter = "  ⚠️ 会被min_cells=3过滤!"
                    elif n_cells < 10:
                        will_filter = f"  ⚠️ 仅{n_cells}个细胞,接近过滤边缘"
                    else:
                        will_filter = "  ✓ 安全"
                print(f"  {dataset:8s}: {n_cells:>6,} 细胞 ({pct:>5.2f}%) | "
                      f"平均counts: {mean:>6.1f} | 总counts: {total:>10.0f}{will_filter}")

    # 按上皮亚型分类统计
    print("\n" + "=" * 80)
    print("按上皮亚型分类的marker表达情况")
    print("=" * 80)

    for celltype, genes in extended_markers.items():
        print(f"\n{celltype}:")
        print("-" * 60)
        for gene in genes:
            row = df[(df['Gene'] == gene) & (df['Dataset'] == '上皮')].head(1)
            if row.empty:
                continue
            row = row.iloc[0]
            if row['Status'] == '不存在':
                print(f"  {gene:12s}: ❌ 不存在")
            else:
                n_cells = row['N_Cells_Expr']
                pct = row['Pct_Cells']
                if n_cells == 0:
                    status_mark = "❌ 无表达"
                elif n_cells < 3:
                    status_mark = f"⚠️ 仅{n_cells}细胞 (会被过滤)"
                elif n_cells < 10:
                    status_mark = f"⚠️ {n_cells}细胞 (低)"
                elif pct < 1:
                    status_mark = f"△ {n_cells}细胞 ({pct:.2f}%)"
                else:
                    status_mark = f"✓ {n_cells}细胞 ({pct:.1f}%)"
                print(f"  {gene:12s}: {status_mark}")

    # 过滤风险评估（上皮）
    print("\n" + "=" * 80)
    print("过滤风险评估 (上皮细胞, min_cells=3)")
    print("=" * 80)

    epi_data = df[df['Dataset'] == '上皮'].copy()
    safe = epi_data[epi_data['N_Cells_Expr'] >= 10]
    warning = epi_data[(epi_data['N_Cells_Expr'] >= 3) & (epi_data['N_Cells_Expr'] < 10)]
    danger = epi_data[(epi_data['N_Cells_Expr'] > 0) & (epi_data['N_Cells_Expr'] < 3)]
    will_lost = epi_data[epi_data['Status'] == '不存在']
    no_expr = epi_data[(epi_data['Status'] == '存在') & (epi_data['N_Cells_Expr'] == 0)]

    print(f"\n✓ 安全 (≥10细胞): {len(safe)} 个基因")
    if len(safe) > 0:
        safe_genes = safe['Gene'].tolist()
        print(f"  {', '.join(safe_genes[:10])}" + (f"\n  ... 及其他 {len(safe_genes)-10} 个" if len(safe_genes) > 10 else ""))

    print(f"\n△ 警告 (3-9细胞): {len(warning)} 个基因")
    if len(warning) > 0:
        for _, row in warning.iterrows():
            marker = "★" if row['Gene'] in missing_genes else " "
            print(f"  {marker} {row['Gene']:12s}: {row['N_Cells_Expr']} 细胞")

    print(f"\n⚠️ 危险 (1-2细胞, 会被过滤): {len(danger)} 个基因")
    if len(danger) > 0:
        for _, row in danger.iterrows():
            marker = "★" if row['Gene'] in missing_genes else " "
            print(f"  {marker} {row['Gene']:12s}: {row['N_Cells_Expr']} 细胞 ← 会被min_cells=3过滤")

    print(f"\n❌ 原始数据中不存在: {len(will_lost)} 个基因")
    if len(will_lost) > 0:
        lost_genes = will_lost['Gene'].tolist()
        print(f"  {', '.join(lost_genes)}")

    print(f"\n○ 基因存在但无表达: {len(no_expr)} 个基因")
    if len(no_expr) > 0:
        for _, row in no_expr.head(5).iterrows():
            print(f"  {row['Gene']:12s}: 存在但无细胞表达")

    # 保存结果
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print("\n" + "=" * 80)
    print(f"✓ 详细结果已保存: {OUTPUT_CSV}")
    print("=" * 80)

    # 建议与摘要
    print("\n" + "=" * 80)
    print("💡 解决建议")
    print("=" * 80)

    n_will_filter = len(danger)
    n_missing_in_data = len(will_lost)
    n_missing_dotplot = len([g for g in missing_genes if g in danger['Gene'].tolist() or g in will_lost['Gene'].tolist()])

    print(f"""
检测摘要:
  - 检测基因总数: {len(all_markers)}
  - Dotplot缺失基因: {len(missing_genes)}
  - 会被min_cells=3过滤: {n_will_filter}
  - 数据中完全不存在: {n_missing_in_data}
  - Dotplot缺失且需处理: {n_missing_dotplot}

根据检测结果:
""".rstrip())

    if n_will_filter > 0:
        print("1. 【推荐】修改过滤参数：")
        print("   min_cells_per_gene = 0  # 或 1")
        print(f"   这将保留 {n_will_filter} 个低表达但重要的marker基因\n")

    if n_missing_in_data > 0:
        print("2. 【注意】以下基因在原始数据中不存在，无法恢复：")
        print(f"   {', '.join(will_lost['Gene'].tolist())}\n")

    # 定义 available_genes 以便快速作图
    available_genes = sorted([g for g in all_markers if lookup_epi(g)[0] is not None])

    print("3. 【快速方案】使用上皮子集绘制 dotplot：")
    print("   import scanpy as sc")
    print("   # 你已有 adata_main / adata_epi；若在独立会话，请先重新读取/构造 adata_epi")
    print("   sc.pl.dotplot(adata_epi, var_names=available_genes, groupby='leiden_bbknn')")

    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()
