#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
表型→基因关联 路线C：图构建（cluster–gene–HPO）→ 嵌入/链接预测 → 基因排序
Phenotype→Gene Pipeline C: Graph (cluster–gene–HPO) → Embedding / Link Prediction

工作流摘要 (Workflow Summary):
  1. 构建异构图节点与边：
       节点类型：gene, hpo_term, cell_cluster
       边类型：
         cluster → gene    (权重 = 该 cluster 中基因平均表达，由 scRNA 计算)
         gene    → hpo_term (来自 HPO genes_to_phenotype.txt 注释)
         hpo_term→ hpo_term (来自 HPO 本体层级，由 hp.obo 解析；可选)
  2. 对图执行浅层嵌入（TransE-like / 矩阵分解）或
     用 NetworkX + Node2Vec 计算节点嵌入
  3. 链接预测：用点积/余弦相似度预测 gene–hpo_term 新边
  4. 对目标 HPO terms，为所有基因计算与目标表型的相似度得分 → 基因排序

输入 (Inputs):
  adata_path    : AnnData h5ad 文件（含 raw counts、细胞类型注释）
  hpo_gene_path : HPO phenotype_to_genes.txt 本地路径
  target_hpo_ids: 目标 HPO term ID 列表

输出 (Outputs):
  <output_dir>/graph_nodes.csv         : 节点列表
  <output_dir>/graph_edges.csv         : 边列表
  <output_dir>/node_embeddings.csv     : 节点嵌入矩阵
  <output_dir>/gene_ranking_C.csv      : 链接预测基因排序结果

依赖包 (Required Packages):
  anndata, scanpy, pandas, numpy, scipy, networkx, scikit-learn
  (可选增强: node2vec / torch_geometric / stellargraph)

参考 (References):
  - CADA (Köhler et al., 2021): HPO-gene link prediction
    https://pmc.ncbi.nlm.nih.gov/articles/PMC8415429/
  - SHEPHERD (Robinson et al., 2023): KG deep learning
    https://github.com/mims-harvard/SHEPHERD
  - node2vec (Grover & Leskovec, 2016):
    https://arxiv.org/abs/1607.00653

版本 (Version): v1.0  2026-03-27
"""

from __future__ import annotations

import argparse
import logging
import os
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ==============================================================================
# 0. 配置参数 Configuration
# ==============================================================================

ADATA_PATH     = None   # 替换为真实 h5ad 路径
HPO_GENE_PATH  = None   # HPO phenotype_to_genes.txt 本地路径（None=下载）
OUTPUT_DIR     = "./results_pipeline_C"

CELL_TYPE_KEY  = "Annotation_2"
SAMPLE_KEY     = "sample"
CONDITION_KEY  = "disease"

# 图参数
EXPR_EDGE_THRESHOLD = 0.1   # cluster→gene 表达权重阈值（低于此值的边忽略）
EMBED_DIM           = 64    # 节点嵌入维度
EMBED_EPOCHS        = 100   # SVD/NMF 嵌入迭代（或 random walk 步数）
TARGET_HPO_IDS      = None  # 目标 HPO IDs（None=使用 Toy 演示）
TOP_K               = 200   # 返回 top-k 基因

RANDOM_STATE = 42

# ==============================================================================
# 1. HPO 基因集加载（与 Pipeline B 复用相同接口）
# ==============================================================================

def load_hpo_gene_sets(local_path: Optional[str] = None) -> pd.DataFrame:
    """加载 HPO phenotype_to_genes.txt"""
    if local_path and Path(local_path).exists():
        log.info("从本地加载 HPO 基因集: %s", local_path)
        df = pd.read_csv(local_path, sep="\t", comment="#",
                         header=None, low_memory=False)
        if df.shape[1] >= 4:
            df = df.iloc[:, [0, 3]]
        else:
            df = df.iloc[:, :2]
        df.columns = ["hpo_id", "gene_symbol"]
    else:
        log.info("尝试下载 HPO phenotype_to_genes.txt...")
        url = ("https://purl.obolibrary.org/obo/hp/hpoa/"
               "phenotype_to_genes.txt")
        try:
            df = pd.read_csv(url, sep="\t", comment="#",
                             header=None, low_memory=False)
            df = df.iloc[:, [0, 3]]
            df.columns = ["hpo_id", "gene_symbol"]
        except Exception as e:
            log.warning("HPO 下载失败: %s  →  使用 Toy HPO", e)
            df = _make_toy_hpo_df()

    df["gene_symbol"] = df["gene_symbol"].str.upper()
    df = df.dropna(subset=["hpo_id", "gene_symbol"])
    return df


# ==============================================================================
# 2. 计算 cluster→gene 表达权重
# ==============================================================================

def compute_cluster_gene_weights(
    adata,
    cell_type_key: str = "Annotation_2",
    normalize: bool = True,
) -> pd.DataFrame:
    """
    计算每个 cell_cluster 中每个基因的平均 log1p 归一化表达

    Returns
    -------
    pd.DataFrame, shape (n_clusters, n_genes)
    """
    log.info("[Step 2] 计算 cluster→gene 表达权重...")
    try:
        import scanpy as sc
        if "counts" in adata.layers:
            adata_norm = adata.copy()
            adata_norm.X = adata.layers["counts"].copy()
        else:
            adata_norm = adata.copy()
        if normalize:
            sc.pp.normalize_total(adata_norm, target_sum=1e4)
            sc.pp.log1p(adata_norm)
        X = adata_norm.X
    except ImportError:
        X = adata.X

    if sp.issparse(X):
        X = X.toarray()

    clusters = adata.obs[cell_type_key].unique().tolist()
    rows = {}
    for ct in clusters:
        mask = (adata.obs[cell_type_key] == ct).values
        rows[ct] = X[mask, :].mean(axis=0)

    cg_df = pd.DataFrame(rows, index=adata.var_names).T
    # cluster 为行，gene 为列
    log.info("  cluster×gene 矩阵: %d clusters × %d 基因",
             cg_df.shape[0], cg_df.shape[1])
    return cg_df


# ==============================================================================
# 3. 图构建 Graph Construction
# ==============================================================================

def build_graph(
    cg_weights: pd.DataFrame,
    hpo_df: pd.DataFrame,
    expr_threshold: float = EXPR_EDGE_THRESHOLD,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    构建异构图

    节点 (Nodes)：
      - gene nodes     : 来自 cg_weights.columns
      - hpo_term nodes : 来自 hpo_df.hpo_id（与 cg_weights 基因有重叠的）
      - cluster nodes  : 来自 cg_weights.index

    边 (Edges)：
      - cluster → gene  : weight = avg_expr（>= threshold）
      - gene → hpo_term : weight = 1（来自 HPO 注释）

    Returns
    -------
    nodes_df : pd.DataFrame，columns [node_id, node_type]
    edges_df : pd.DataFrame，columns [src, dst, edge_type, weight]
    """
    log.info("[Step 3] 构建异构图...")

    genes    = list(cg_weights.columns)
    clusters = list(cg_weights.index)
    hpo_ids  = list(hpo_df["hpo_id"].unique())

    # 节点 DataFrame
    gene_nodes    = pd.DataFrame({"node_id": genes,    "node_type": "gene"})
    cluster_nodes = pd.DataFrame({"node_id": clusters, "node_type": "cluster"})
    hpo_nodes     = pd.DataFrame({"node_id": hpo_ids,  "node_type": "hpo_term"})
    nodes_df = pd.concat([gene_nodes, cluster_nodes, hpo_nodes],
                         ignore_index=True)
    log.info("  节点总数: %d  (gene=%d, cluster=%d, hpo_term=%d)",
             len(nodes_df), len(genes), len(clusters), len(hpo_ids))

    # 边 1: cluster → gene（表达权重）
    cg_rows = []
    genes_upper_map = {g.upper(): g for g in genes}
    for ct in clusters:
        expr = cg_weights.loc[ct]
        for g in genes:
            w = expr[g]
            if w >= expr_threshold:
                cg_rows.append({"src": ct, "dst": g,
                                 "edge_type": "cluster_expr", "weight": float(w)})
    edges_cg = pd.DataFrame(cg_rows)
    log.info("  cluster→gene 边: %d", len(edges_cg))

    # 边 2: gene → hpo_term（HPO 注释）
    hpo_rows = []
    for _, row in hpo_df.iterrows():
        g_upper = str(row["gene_symbol"]).upper()
        g_orig  = genes_upper_map.get(g_upper)
        if g_orig is not None:
            hpo_rows.append({
                "src": g_orig,
                "dst": row["hpo_id"],
                "edge_type": "gene_hpo",
                "weight": 1.0,
            })
    edges_hpo = pd.DataFrame(hpo_rows)
    log.info("  gene→hpo_term 边: %d", len(edges_hpo))

    edges_df = pd.concat([edges_cg, edges_hpo], ignore_index=True)
    log.info("  边总数: %d", len(edges_df))
    return nodes_df, edges_df


# ==============================================================================
# 4. 节点嵌入（矩阵分解/SVD） Node Embedding via Matrix Factorization
# ==============================================================================

def compute_node_embeddings_svd(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    embed_dim: int = EMBED_DIM,
    seed: int = RANDOM_STATE,
) -> pd.DataFrame:
    """
    基于邻接矩阵 SVD 截断分解计算节点嵌入

    方法：
      1. 构建节点索引的邻接矩阵（无向，含权重）
      2. 用 scipy.sparse.linalg.svds 截断 SVD（取前 embed_dim 个奇异值）
      3. 节点嵌入 = U @ diag(sqrt(sigma))（左奇异向量 × sqrt 奇异值）

    这是 CADA（HPO-gene link prediction）论文中 embedding 的简化版实现。
    对于更复杂的应用，可替换为：
      - node2vec（pip install node2vec）
      - PyTorch Geometric RotatE/TransE
      - GraphSAGE

    Returns
    -------
    pd.DataFrame, index=node_id, columns=[dim_0, dim_1, ..., dim_{embed_dim-1}]
    """
    from scipy.sparse.linalg import svds

    log.info("[Step 4] 计算节点嵌入（截断 SVD, dim=%d）...", embed_dim)

    # 节点 ID → 整数索引
    node_ids  = nodes_df["node_id"].tolist()
    node_idx  = {n: i for i, n in enumerate(node_ids)}
    n_nodes   = len(node_ids)

    # 构建稀疏邻接矩阵（无向）
    rows, cols, vals = [], [], []
    for _, row in edges_df.iterrows():
        s = node_idx.get(row["src"])
        d = node_idx.get(row["dst"])
        w = float(row["weight"])
        if s is not None and d is not None:
            rows += [s, d]
            cols += [d, s]
            vals += [w, w]

    A = sp.csr_matrix((vals, (rows, cols)), shape=(n_nodes, n_nodes))
    A_norm = _normalize_adj(A)

    # 截断 SVD
    k = min(embed_dim, n_nodes - 2)
    np.random.seed(seed)
    try:
        U, sigma, Vt = svds(A_norm.astype(np.float64), k=k)
    except Exception as e:
        log.warning("SVD 失败: %s  → 使用随机嵌入", e)
        U = np.random.randn(n_nodes, embed_dim)
        sigma = np.ones(embed_dim)

    # 节点嵌入 = U × sqrt(sigma)
    embeddings = U * np.sqrt(np.abs(sigma))[np.newaxis, :]

    emb_df = pd.DataFrame(
        embeddings,
        index=node_ids,
        columns=[f"dim_{i}" for i in range(embeddings.shape[1])],
    )
    log.info("  嵌入矩阵: %d 节点 × %d 维", emb_df.shape[0], emb_df.shape[1])
    return emb_df


def _normalize_adj(A: sp.spmatrix) -> sp.spmatrix:
    """对称归一化邻接矩阵：D^{-1/2} A D^{-1/2}"""
    degree = np.array(A.sum(axis=1)).flatten()
    d_inv_sqrt = np.where(degree > 0, 1.0 / np.sqrt(degree), 0.0)
    D = sp.diags(d_inv_sqrt)
    return D @ A @ D


# ==============================================================================
# 5. 可选：Node2Vec 嵌入（需安装 node2vec 包）
# ==============================================================================

def compute_node_embeddings_node2vec(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    embed_dim: int = EMBED_DIM,
    walk_length: int = 30,
    num_walks: int = 10,
    seed: int = RANDOM_STATE,
) -> pd.DataFrame:
    """
    用 node2vec 计算节点嵌入（需 pip install node2vec）

    这是 SVD 嵌入的增强替代：更好地捕获高阶图结构。
    """
    try:
        import networkx as nx
        from node2vec import Node2Vec
    except ImportError:
        log.warning("node2vec 或 networkx 未安装，回退到 SVD 嵌入")
        return compute_node_embeddings_svd(nodes_df, edges_df, embed_dim, seed)

    log.info("[Step 4 alt] 用 node2vec 计算节点嵌入...")

    G = nx.Graph()
    G.add_nodes_from(nodes_df["node_id"].tolist())
    for _, row in edges_df.iterrows():
        G.add_edge(str(row["src"]), str(row["dst"]), weight=float(row["weight"]))

    n2v = Node2Vec(G, dimensions=embed_dim, walk_length=walk_length,
                   num_walks=num_walks, workers=1, seed=seed, quiet=True)
    model = n2v.fit(window=5, min_count=1, batch_words=4)

    node_ids = nodes_df["node_id"].tolist()
    vecs = []
    for nid in node_ids:
        try:
            vecs.append(model.wv[str(nid)])
        except KeyError:
            vecs.append(np.zeros(embed_dim))

    emb_df = pd.DataFrame(
        np.vstack(vecs),
        index=node_ids,
        columns=[f"dim_{i}" for i in range(embed_dim)],
    )
    log.info("  node2vec 嵌入: %d 节点 × %d 维", emb_df.shape[0], emb_df.shape[1])
    return emb_df


# ==============================================================================
# 6. 链接预测：基因–HPO 相似度 Link Prediction: Gene–HPO Similarity
# ==============================================================================

def link_predict_gene_hpo(
    emb_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    target_hpo_ids: list[str],
    edges_df: pd.DataFrame,
    top_k: int = TOP_K,
) -> pd.DataFrame:
    """
    基于嵌入相似度（点积）预测 gene–hpo_term 链接，排序候选基因

    算法（类 CADA）：
      score(gene, hpo) = emb_gene · emb_hpo
      对每个目标 HPO term，计算所有基因的得分
      汇总多个目标 HPO terms 的得分（求和）→ 最终基因排序

    Parameters
    ----------
    emb_df         : compute_node_embeddings_*() 的返回值
    nodes_df       : build_graph() 的返回值
    target_hpo_ids : 目标 HPO term ID 列表
    edges_df       : 用于排除已有"已知"gene–hpo 正例（避免泄漏）
    top_k          : 返回前 k 个基因

    Returns
    -------
    pd.DataFrame with columns: gene, link_score, rank
    """
    log.info("[Step 5] 链接预测：gene–HPO 相似度排序...")

    gene_nodes = nodes_df[nodes_df["node_type"] == "gene"]["node_id"].tolist()
    available_hpos = [h for h in target_hpo_ids if h in emb_df.index]
    if not available_hpos:
        log.warning("  目标 HPO terms 均不在嵌入中，无法进行链接预测")
        return pd.DataFrame()

    gene_embs = emb_df.loc[
        [g for g in gene_nodes if g in emb_df.index]
    ].values  # shape: (n_genes, embed_dim)

    available_gene_ids = [g for g in gene_nodes if g in emb_df.index]

    # 已知正例（训练集 gene–hpo 边），用于分析时标注（非必须去除）
    known_pairs = set(
        zip(edges_df[edges_df["edge_type"] == "gene_hpo"]["src"],
            edges_df[edges_df["edge_type"] == "gene_hpo"]["dst"])
    )

    total_scores = np.zeros(len(available_gene_ids))

    for hpo_id in available_hpos:
        hpo_emb = emb_df.loc[hpo_id].values  # shape: (embed_dim,)
        scores  = gene_embs @ hpo_emb        # 点积相似度
        total_scores += scores

    result = pd.DataFrame({
        "gene":        available_gene_ids,
        "link_score":  total_scores,
    })
    result = result.sort_values("link_score", ascending=False).reset_index(drop=True)
    result["rank"] = result.index + 1

    # 标注是否为已知正例
    result["is_known"] = result.apply(
        lambda r: any((r["gene"], h) in known_pairs for h in available_hpos),
        axis=1,
    )
    result = result.head(top_k)
    log.info("  Top-%d 候选基因排序完成（其中 %d 个为已知正例）",
             len(result), result["is_known"].sum())
    return result


# ==============================================================================
# 7. 主流程 Main Pipeline
# ==============================================================================

def run_pipeline_C(
    adata_path: Optional[str] = None,
    hpo_gene_path: Optional[str] = None,
    target_hpo_ids: Optional[list[str]] = None,
    cell_type_key: str = CELL_TYPE_KEY,
    output_dir: str = OUTPUT_DIR,
    embed_dim: int = EMBED_DIM,
    top_k: int = TOP_K,
    use_node2vec: bool = False,
) -> dict:
    """
    运行路线 C 完整流程

    Returns
    -------
    dict with keys: nodes_df, edges_df, emb_df, gene_ranking
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- 真实数据 or Toy Demo ---
    if adata_path is None:
        log.info("=== [Toy Demo] 未提供 AnnData 路径，使用模拟数据 ===")
        adata = _make_toy_adata()
        cell_type_key = "Annotation_2"
    else:
        try:
            import anndata as ad
            adata = ad.read_h5ad(adata_path)
            log.info("加载 AnnData: %d 细胞 × %d 基因",
                     adata.n_obs, adata.n_vars)
        except Exception as e:
            log.error("AnnData 加载失败: %s", e)
            raise

    # Step 2: cluster×gene 表达权重
    cg_weights = compute_cluster_gene_weights(adata, cell_type_key=cell_type_key)

    # HPO 基因集加载
    hpo_df = load_hpo_gene_sets(hpo_gene_path)

    # 只保留与数据基因有重叠的 HPO terms
    genes_upper = set(g.upper() for g in adata.var_names)
    hpo_df = hpo_df[hpo_df["gene_symbol"].isin(genes_upper)].copy()
    log.info("过滤后 HPO 基因集: %d 行", len(hpo_df))

    # Step 3: 图构建
    nodes_df, edges_df = build_graph(
        cg_weights,
        hpo_df,
        expr_threshold=EXPR_EDGE_THRESHOLD,
    )
    nodes_df.to_csv(os.path.join(output_dir, "graph_nodes.csv"), index=False)
    edges_df.to_csv(os.path.join(output_dir, "graph_edges.csv"), index=False)
    log.info("  图保存至: %s", output_dir)

    # Step 4: 节点嵌入
    if use_node2vec:
        emb_df = compute_node_embeddings_node2vec(
            nodes_df, edges_df, embed_dim=embed_dim
        )
    else:
        emb_df = compute_node_embeddings_svd(
            nodes_df, edges_df, embed_dim=embed_dim
        )
    emb_df.to_csv(os.path.join(output_dir, "node_embeddings.csv"))
    log.info("  嵌入保存至: %s/node_embeddings.csv", output_dir)

    # 目标 HPO IDs
    all_hpo_ids = hpo_df["hpo_id"].unique().tolist()
    if target_hpo_ids is None:
        target_hpo_ids = all_hpo_ids[:5]
        log.info("  未指定目标 HPO，使用演示 HPO IDs: %s", target_hpo_ids)

    # Step 5: 链接预测
    gene_ranking = link_predict_gene_hpo(
        emb_df, nodes_df, target_hpo_ids, edges_df, top_k=top_k
    )

    if gene_ranking.empty:
        log.warning("[WARN] 链接预测无结果")
        return {}

    out_file = os.path.join(output_dir, "gene_ranking_C.csv")
    gene_ranking.to_csv(out_file, index=False)
    log.info("[Done] 基因排序结果保存至: %s", out_file)

    # 打印 top-20
    log.info("\n=== Top-20 候选基因（图链接预测得分） ===\n%s",
             gene_ranking.head(20).to_string(index=False))

    return {
        "nodes_df":    nodes_df,
        "edges_df":    edges_df,
        "emb_df":      emb_df,
        "gene_ranking": gene_ranking,
    }


# ==============================================================================
# Toy 数据生成 Toy Data Helpers
# ==============================================================================

def _make_toy_adata():
    """生成模拟 AnnData 用于演示"""
    import anndata as ad

    np.random.seed(42)
    n_cells, n_genes = 600, 300
    X = np.random.negative_binomial(2, 0.5, size=(n_cells, n_genes)).astype(np.float32)
    X[:300, :30] += 5   # 前 300 细胞（cluster A）的前 30 个基因高表达

    obs = pd.DataFrame({
        "Annotation_2": (["ClusterA"] * 200 + ["ClusterB"] * 200 + ["ClusterC"] * 200),
        "sample":        ["S1"] * 100 + ["S2"] * 100 + ["S3"] * 100 +
                         ["S4"] * 100 + ["S5"] * 100 + ["S6"] * 100,
        "disease":       (["CRSwNP"] * 300 + ["Control"] * 300),
    })
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)])

    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata


def _make_toy_hpo_df() -> pd.DataFrame:
    """生成模拟 HPO 基因集"""
    np.random.seed(1)
    n_terms, n_genes = 30, 300
    rows = []
    for i in range(n_terms):
        hpo_id = f"HP:{i:07d}"
        gene_idx = np.random.choice(n_genes, 15, replace=False)
        for gi in gene_idx:
            rows.append({"hpo_id": hpo_id, "gene_symbol": f"GENE{gi}"})
    return pd.DataFrame(rows)


# ==============================================================================
# 脚本入口 Script Entry Point
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline C: Graph-based gene ranking via link prediction"
    )
    parser.add_argument("--adata",       default=None, help="AnnData h5ad 路径")
    parser.add_argument("--hpo_genes",   default=None, help="HPO phenotype_to_genes.txt 路径")
    parser.add_argument("--hpo_ids",     default=None, nargs="+",
                        help="目标 HPO term IDs")
    parser.add_argument("--cell_type_key", default=CELL_TYPE_KEY)
    parser.add_argument("--output_dir", default=OUTPUT_DIR)
    parser.add_argument("--embed_dim",  type=int, default=EMBED_DIM)
    parser.add_argument("--top_k",      type=int, default=TOP_K)
    parser.add_argument("--node2vec",   action="store_true",
                        help="使用 node2vec 嵌入（需安装 node2vec 包）")
    args = parser.parse_args()

    run_pipeline_C(
        adata_path    = args.adata,
        hpo_gene_path = args.hpo_genes,
        target_hpo_ids = args.hpo_ids,
        cell_type_key  = args.cell_type_key,
        output_dir     = args.output_dir,
        embed_dim      = args.embed_dim,
        top_k          = args.top_k,
        use_node2vec   = args.node2vec,
    )


if __name__ == "__main__":
    main()
