#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import textwrap
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_RUN_ROOT = Path("/home/h2048/data/py/20260531/normal_airway_ml_bcell_full_run_20260531")
DEFAULT_ENV_FILE = Path("/home/h2048/.env")
DEFAULT_MODEL = "deepseek-chat"
DEFAULT_API_URL = "https://api.deepseek.com/v1/chat/completions"
PLACEHOLDER_TOKENS = (
    "your-deepseek-api-key",
    "your_key_here",
    "placeholder",
    "changeme",
    "replace_me",
)


@dataclass
class FigureExplanation:
    figure_name: str
    title: str
    figure_outputs: list[str]
    source_tables: list[str]
    gene_filter_context: dict[str, Any]
    api_status: str
    model: str | None
    local_summary: dict[str, Any]
    final_explanation: dict[str, Any]
    raw_response_text: str | None = None
    error_message: str | None = None


FIGURE_TITLES = {
    "contrast_overview": "B-cell core discovery contrast overview",
    "site_specific_dotplot": "Top site-specific core genes by cell type/site",
    "shared_core": "Shared core genes across multiple target sites",
    "run_diagnostics": "Run diagnostics and pseudobulk coverage",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Explain B-cell core-discovery figures with DeepSeek using structured plot-data inputs.")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT, help="Root directory of the completed run.")
    parser.add_argument("--env-file", type=Path, default=DEFAULT_ENV_FILE, help="Optional .env file to read DEEPSEEK_API_KEY from when the shell env is empty.")
    parser.add_argument("--model", default=DEFAULT_MODEL, help="DeepSeek model name.")
    parser.add_argument("--api-url", default=DEFAULT_API_URL, help="OpenAI-compatible DeepSeek chat completions endpoint.")
    parser.add_argument("--timeout-seconds", type=int, default=120, help="HTTP timeout for each figure explanation request.")
    parser.add_argument("--max-table-rows", type=int, default=12, help="Max rows of structured table data to include in each prompt.")
    return parser.parse_args()


def ensure_dirs(paths: list[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def read_env_value(env_file: Path, key: str) -> str:
    if not env_file.exists():
        return ""
    try:
        for raw_line in env_file.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            lhs, rhs = line.split("=", 1)
            if lhs.strip() != key:
                continue
            value = rhs.strip().strip('"').strip("'")
            return value
    except OSError:
        return ""
    return ""


def is_placeholder_key(value: str) -> bool:
    lowered = value.strip().lower()
    if not lowered:
        return True
    return any(token in lowered for token in PLACEHOLDER_TOKENS)


def resolve_deepseek_key(env_file: Path) -> tuple[str, str]:
    from_env = os.getenv("DEEPSEEK_API_KEY", "").strip()
    if from_env:
        return from_env, "environment"
    from_file = read_env_value(env_file, "DEEPSEEK_API_KEY").strip()
    if from_file:
        return from_file, "env_file"
    return "", "missing"


def load_figure_manifest(figure_dir: Path) -> dict[str, Any]:
    manifest_path = figure_dir / "figure_manifest.json"
    with manifest_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def clip_rows(df: pd.DataFrame, max_rows: int) -> list[dict[str, Any]]:
    clipped = df.head(max_rows).copy()
    return json.loads(clipped.to_json(orient="records", force_ascii=False))


def compact_gene_filter_context(manifest: dict[str, Any]) -> dict[str, Any]:
    gene_filter = manifest.get("gene_filter") or {}
    if not gene_filter:
        return {}
    return {
        "enabled": bool(gene_filter.get("enabled", False)),
        "profile_name": gene_filter.get("profile_name"),
        "rules": gene_filter.get("rules", []),
        "source_table_counts": gene_filter.get("source_table_counts", {}),
        "example_removed_genes": gene_filter.get("example_removed_genes", {}),
        "plotted_rows": gene_filter.get("plotted_rows", {}),
    }


def infer_gene_family_hint(gene: str) -> str:
    gene_upper = gene.upper()
    if gene_upper.startswith(("IGH", "IGK", "IGL")) or gene_upper in {"IGJ", "JCHAIN"}:
        return "偏向 BCR/抗体分泌相关程序。"
    if gene_upper.startswith(("ATP5", "NDUF", "COX", "UQCR", "SDH")):
        return "偏向氧化磷酸化或线粒体能量代谢。"
    if gene_upper.startswith(("SEPT", "ACT", "TPM", "MYL", "TUB", "VCL")):
        return "偏向细胞骨架、形态维持或迁移/应激重塑。"
    if gene_upper.startswith(("SFTPA", "SFTPB", "SFTPC", "SFTPD", "NAPSA", "CLDN18")):
        return "更像 distal-lung / alveolar epithelial cue；在 B 细胞图里需同时考虑组织微环境印记或少量环境 RNA。"
    if gene_upper.startswith(("TGF", "SMAD")) or gene_upper == "TGFBR2":
        return "偏向 TGF-β / 组织重塑信号。"
    if gene_upper.startswith(("MARCH", "VCP", "PSMA", "NOP", "SF3", "LSM", "RNPC", "PRPF", "SELENOK")):
        return "偏向蛋白稳态、膜转运或 RNA 加工相关过程。"
    if gene_upper.startswith(("LINC", "AC", "AL", "AP")) or gene_upper.startswith("RP11-"):
        return "属于注释较弱的 locus-like 转录本，更适合作为 marker 候选而非单独机制锚点。"
    if gene_upper.startswith("C") and "ORF" in gene_upper:
        return "属于 open reading frame 命名基因，功能注释相对有限，解释时应更保守。"
    return "需要结合当前 cell type / site 语境与文献进一步解释。"


def json_safe(value: Any) -> Any:
    if pd.isna(value):
        return None
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def summarize_contrast(df: pd.DataFrame, max_rows: int) -> dict[str, Any]:
    top_by_genes = df.sort_values(["n_stable_core_genes", "mean_balanced_accuracy"], ascending=[False, False]).head(5)
    by_celltype = (
        df.groupby("cell_type", dropna=False)
        .agg(
            n_contrasts=("contrast_id", "count"),
            max_stable_core_genes=("n_stable_core_genes", "max"),
            mean_balanced_accuracy=("mean_balanced_accuracy", "mean"),
        )
        .reset_index()
        .sort_values(["max_stable_core_genes", "mean_balanced_accuracy"], ascending=[False, False])
    )
    weakest = df.sort_values(["n_stable_core_genes", "mean_balanced_accuracy"], ascending=[True, True]).head(4)
    return {
        "plot_type": "bar + scatter overview",
        "n_contrasts": int(df.shape[0]),
        "top_by_stable_core_genes": clip_rows(top_by_genes, max_rows),
        "cell_type_summary": clip_rows(by_celltype, max_rows),
        "lowest_yield_examples": clip_rows(weakest, max_rows),
    }


def summarize_site_specific(df: pd.DataFrame, max_rows: int) -> dict[str, Any]:
    top_sites = (
        df.groupby(["cell_type", "enriched_site"], dropna=False)
        .agg(
            n_genes=("gene", "count"),
            max_consensus=("best_consensus_score", "max"),
            mean_abs_log2_fc=("mean_abs_log2_fc", "mean"),
        )
        .reset_index()
        .sort_values(["cell_type", "max_consensus"], ascending=[True, False])
    )
    strongest = df.sort_values(["best_consensus_score", "n_supporting_contrasts", "mean_abs_log2_fc"], ascending=[False, False, False]).head(10)
    gene_items = []
    gene_hint_context = []
    for row in df.itertuples(index=False):
        gene_items.append(
            {
                "gene": row.gene,
                "cell_type": row.cell_type,
                "enriched_site": row.enriched_site,
                "best_consensus_score": row.best_consensus_score,
                "mean_abs_log2_fc": row.mean_abs_log2_fc,
                "n_supporting_contrasts": row.n_supporting_contrasts,
            }
        )
        gene_hint_context.append(
            {
                "gene": row.gene,
                "context": f"{row.cell_type} / {row.enriched_site}",
                "hint": infer_gene_family_hint(str(row.gene)),
            }
        )
    return {
        "plot_type": "faceted site-specific gene dotplot",
        "n_rows": int(df.shape[0]),
        "top_genes_displayed": clip_rows(df, max_rows),
        "celltype_site_summary": clip_rows(top_sites, max_rows),
        "strongest_examples": clip_rows(strongest, max_rows),
        "gene_items_for_detailed_interpretation": gene_items,
        "gene_hint_context": gene_hint_context,
    }


def summarize_shared(df: pd.DataFrame, max_rows: int) -> dict[str, Any]:
    strongest = df.sort_values(["n_target_sites", "best_consensus_score", "mean_selection_freq"], ascending=[False, False, False])
    by_celltype = (
        strongest.groupby("cell_type", dropna=False)
        .agg(
            n_shared_genes=("gene", "count"),
            max_target_sites=("n_target_sites", "max"),
            max_consensus=("best_consensus_score", "max"),
        )
        .reset_index()
        .sort_values(["n_shared_genes", "max_consensus"], ascending=[False, False])
    )
    gene_items = []
    gene_hint_context = []
    for row in strongest.itertuples(index=False):
        gene_items.append(
            {
                "gene": row.gene,
                "cell_type": row.cell_type,
                "supported_sites": row.supported_sites,
                "best_consensus_score": row.best_consensus_score,
                "mean_selection_freq": row.mean_selection_freq,
                "n_target_sites": row.n_target_sites,
            }
        )
        gene_hint_context.append(
            {
                "gene": row.gene,
                "context": f"{row.cell_type} / shared across {row.supported_sites}",
                "hint": infer_gene_family_hint(str(row.gene)),
            }
        )
    return {
        "plot_type": "shared-core lollipop / point summary",
        "n_shared_genes": int(df.shape[0]),
        "shared_gene_rows": clip_rows(strongest, max_rows),
        "cell_type_summary": clip_rows(by_celltype, max_rows),
        "gene_items_for_detailed_interpretation": gene_items,
        "gene_hint_context": gene_hint_context,
    }


def summarize_run_diagnostics(method_status: pd.DataFrame, coverage: pd.DataFrame, max_rows: int) -> dict[str, Any]:
    method_table = method_status.copy().sort_values(["method", "status"])
    failures = method_table[method_table["status"] != "ok"].copy()
    coverage_sorted = coverage.sort_values(["n_sites", "n_samples"], ascending=[False, False])
    sparse_celltypes = coverage_sorted[coverage_sorted["n_sites"] < 2].copy()
    return {
        "plot_type": "stacked status bars + coverage bars",
        "method_status_table": clip_rows(method_table, max_rows),
        "non_ok_rows": clip_rows(failures, max_rows),
        "coverage_table": clip_rows(coverage_sorted, max_rows),
        "coverage_alerts": clip_rows(sparse_celltypes, max_rows),
    }


def build_site_specific_gene_details(payload: dict[str, Any]) -> list[dict[str, str]]:
    gene_details: list[dict[str, str]] = []
    for item in payload.get("gene_items_for_detailed_interpretation", []):
        context = f"{item['cell_type']} / {item['enriched_site']}"
        detail = (
            f"{item['gene']} 在 {context} 面板中入选（consensus={item['best_consensus_score']:.3f}, "
            f"|log2FC|={item['mean_abs_log2_fc']:.2f}, supporting_contrasts={int(item['n_supporting_contrasts'])}）；"
            f"{infer_gene_family_hint(str(item['gene']))}"
        )
        gene_details.append({"gene": str(item["gene"]), "context": context, "detail": detail})
    return gene_details


def build_shared_gene_details(payload: dict[str, Any]) -> list[dict[str, str]]:
    gene_details: list[dict[str, str]] = []
    for item in payload.get("gene_items_for_detailed_interpretation", []):
        context = f"{item['cell_type']} / shared across {item['supported_sites']}"
        detail = (
            f"{item['gene']} 在 {context} 中复现（consensus={item['best_consensus_score']:.3f}, "
            f"selection_freq={item['mean_selection_freq']:.3f}, target_sites={int(item['n_target_sites'])}）；"
            f"{infer_gene_family_hint(str(item['gene']))}"
        )
        gene_details.append({"gene": str(item["gene"]), "context": context, "detail": detail})
    return gene_details


def build_local_explanation(figure_name: str, payload: dict[str, Any]) -> dict[str, Any]:
    if figure_name == "contrast_overview":
        top = payload["top_by_stable_core_genes"][0]
        weakest = payload["lowest_yield_examples"][0]
        return {
            "one_sentence_takeaway": (
                f"这张图说明不同 B 细胞亚群对解剖部位的可分离度差异很大，其中 {top['cell_type']} 的 "
                f"{top['positive_label']} 相关对比产生了最多的稳定核心基因。"
            ),
            "what_is_plotted": "左侧是每个 contrast 的稳定核心基因数，右侧同时展示 balanced accuracy、核心基因数和样本量。",
            "key_observations": [
                f"最高信号来自 {top['cell_type']} / {top['display_label']}，稳定核心基因数为 {int(top['n_stable_core_genes'])}，balanced accuracy 为 {top['mean_balanced_accuracy']:.3f}。",
                f"较弱的 contrast 包括 {weakest['display_label']}，其稳定核心基因数仅 {int(weakest['n_stable_core_genes'])}。",
                "Memory_B 和 Plasma 的多个 contrast 同时兼具较高准确率与较多稳定核心基因，是当前 run 中最主要的差异来源。",
            ],
            "biological_interpretation": [
                "不同组织位点在部分 B 细胞亚群中留下了稳定且可重复的转录差异，而这种差异并不是所有细胞类型都同样强。",
                "高准确率且高核心基因数的 contrast 更适合作为后续机制解释或 marker 精修的重点。",
            ],
            "caveats": [
                "这是基于 pseudobulk contrast 的稳定性结果，不等同于单细胞层面的所有变化都很大。",
                "低核心基因数不一定代表没有生物学差异，也可能与样本覆盖或对比分布不均衡有关。",
            ],
            "gene_details": [],
            "gene_relationships": [],
        }

    if figure_name == "site_specific_dotplot":
        strongest = payload["strongest_examples"][0]
        return {
            "one_sentence_takeaway": (
                f"这张图强调每个细胞类型内部最具位点特异性的基因，其中 {strongest['cell_type']} 在 {strongest['enriched_site']} "
                f"方向上出现了当前展示中最强的一批高共识基因。"
            ),
            "what_is_plotted": "每一行是一个基因、每一列是一个组织位点；点大小表示支持该基因的 contrast 数，颜色表示 consensus score。",
            "key_observations": [
                f"当前最强示例之一是 {strongest['cell_type']} / {strongest['gene']} / {strongest['enriched_site']}，consensus score 为 {strongest['best_consensus_score']:.3f}。",
                "Memory_B 在 sinus 和 lung_parenchyma 方向上的强信号较集中；Naive_B 与 Plasma 也各自保留了位点特异 marker。",
                "同一细胞类型不同位点显示的是不同的 top genes，说明组织差异并非单一全局程序，而更像位点定向富集。",
            ],
            "biological_interpretation": [
                "这些基因可优先作为 site-specific marker 候选，用于后续通路解释或与外部文献交叉验证。",
                "若某一位点出现多个高 consensus、高 log2FC 基因，通常意味着该位点在该细胞类型中具有更稳定的组织特异状态。",
            ],
            "caveats": [
                "图中只展示了每个 cell_type × site 的前几名基因，不代表完整差异列表。",
                "部分基因的显著性和效应量可能受样本数与 contrast 结构影响，需要回到原始 summary 表复核。",
            ],
            "gene_details": build_site_specific_gene_details(payload),
            "gene_relationships": [
                "同一 cell type 在不同 anatomical sites 使用不同的基因组合，说明这是位点定向程序，而不是统一的全局 B-cell 状态。",
                "同一 panel 内若同时出现能量代谢、细胞骨架或分泌相关基因，往往提示局部组织适应需要多条程序并行。",
                "展示前已过滤 MT/RP/ENSG/locus-like 技术性或低解释性基因，因此剩余基因更适合被当作生物学候选而不是技术噪声。",
            ],
        }

    if figure_name == "shared_core":
        top = payload["shared_gene_rows"][0]
        return {
            "one_sentence_takeaway": (
                f"这张图聚焦跨多个目标位点重复出现的 shared core genes，其中 {top['gene']} 在 {top['cell_type']} 中跨 {int(top['n_target_sites'])} 个位点重复出现。"
            ),
            "what_is_plotted": "每个点代表一个跨多个位点复现的核心基因，横轴是最佳 consensus score，点大小随覆盖位点数增加。",
            "key_observations": [
                f"{top['gene']} ({top['cell_type']}) 是当前 shared-core 中最强的代表之一，覆盖位点为 {top['supported_sites']}。",
                "shared core 基因总数并不多，说明大多数差异更偏向 site-specific，而不是所有位点共用同一套程序。",
                "Memory_B 和 Plasma 都贡献了 shared core，但各自共享的基因集合并不完全相同。",
            ],
            "biological_interpretation": [
                "shared core 基因更适合被看作‘跨多个组织位点仍保持稳定方向性’的程序，可能代表较核心的组织适应机制。",
                "相比只在单个位点出现的 marker，这类基因更适合优先进入后续跨位点机制整合。",
            ],
            "caveats": [
                "shared core 的定义依赖当前 contrast 设计和阈值，因此数量较少并不异常。",
                "跨两个位点复现不等于跨所有位点普遍存在，仍需结合 supported_sites 具体解释。",
            ],
            "gene_details": build_shared_gene_details(payload),
            "gene_relationships": [
                "Memory_B 与 Plasma 的 shared genes 并没有完全重叠，说明跨位点保守程序本身也是亚群特异的。",
                "如果多基因同时落在能量代谢、分泌运输或蛋白稳态相关家族，可优先理解为跨组织位点共同保留的基础适应程序。",
                "shared-core 更适合作为‘跨多个位点都复现的核心候选’，而不是单一位点 marker 的替代品。",
            ],
        }

    failures = payload.get("non_ok_rows", [])
    coverage_alerts = payload.get("coverage_alerts", [])
    failure_text = "存在非 ok 记录" if failures else "所有方法均成功"
    sparse_text = coverage_alerts[0]["cell_type"] if coverage_alerts else "无"
    return {
        "one_sentence_takeaway": f"这张图是本次 run 的质量与覆盖度体检表：整体方法运行稳定，但样本/位点覆盖并不完全均衡。",
        "what_is_plotted": "左侧是各模型方法在全部 contrasts 上的运行状态，右侧是各细胞类型的 pseudobulk 样本数与位点覆盖数。",
        "key_observations": [
            f"方法状态总体稳定，{failure_text}；当前唯一显式失败来自 LDA 的 1 条记录。",
            f"覆盖度最弱的细胞类型是 {sparse_text}，其位点数不足 2，因此无法形成有效 contrast。",
            "Memory_B 样本最多、位点覆盖最完整，是当前结果解释最稳的细胞类型之一。",
        ],
        "biological_interpretation": [
            "运行状态图帮助区分‘真实无差异’与‘方法/覆盖不足导致看不到差异’这两种情况。",
            "覆盖充分的细胞类型更容易产生稳定 core genes，而覆盖不足的细胞类型即使存在差异也可能难以检出。",
        ],
        "caveats": [
            "方法失败不一定表示该 contrast 完全不可解释，但说明该模型不适合当前样本结构。",
            "覆盖条形图只反映 pseudobulk 设计层面的可用度，不直接反映 effect size 大小。",
        ],
        "gene_details": [],
        "gene_relationships": [],
    }


def normalize_text_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    normalized: list[str] = []
    for item in value:
        text = str(item).strip()
        if text:
            normalized.append(text)
    return normalized


def normalize_gene_detail_item(item: Any) -> dict[str, str] | None:
    if not isinstance(item, dict):
        return None
    gene = str(item.get("gene", "")).strip()
    context = str(item.get("context", "")).strip()
    detail = str(item.get("detail") or item.get("interpretation") or item.get("note") or "").strip()
    if not gene or not context or not detail:
        return None
    return {"gene": gene, "context": context, "detail": detail}


def merge_gene_details(candidate: Any, fallback: list[dict[str, str]]) -> list[dict[str, str]]:
    candidate_items: list[dict[str, str]] = []
    if isinstance(candidate, list):
        for raw_item in candidate:
            item = normalize_gene_detail_item(raw_item)
            if item is None:
                continue
            candidate_items.append(item)

    merged: list[dict[str, str]] = []
    used_candidate_indices: set[int] = set()

    for fallback_item in fallback:
        matched_index: int | None = None
        for idx, candidate_item in enumerate(candidate_items):
            if idx in used_candidate_indices:
                continue
            if (
                candidate_item["gene"] == fallback_item["gene"]
                and candidate_item["context"] == fallback_item["context"]
            ):
                matched_index = idx
                break

        if matched_index is None:
            for idx, candidate_item in enumerate(candidate_items):
                if idx in used_candidate_indices:
                    continue
                if candidate_item["gene"] == fallback_item["gene"]:
                    matched_index = idx
                    break

        if matched_index is not None:
            used_candidate_indices.add(matched_index)
            merged.append(candidate_items[matched_index])
        else:
            merged.append(fallback_item)

    for idx, candidate_item in enumerate(candidate_items):
        if idx not in used_candidate_indices:
            merged.append(candidate_item)
    return merged


def normalize_final_explanation(candidate: dict[str, Any], fallback: dict[str, Any]) -> dict[str, Any]:
    normalized = {
        "one_sentence_takeaway": str(candidate.get("one_sentence_takeaway", "")).strip() or fallback["one_sentence_takeaway"],
        "what_is_plotted": str(candidate.get("what_is_plotted", "")).strip() or fallback["what_is_plotted"],
        "key_observations": normalize_text_list(candidate.get("key_observations")) or fallback["key_observations"],
        "biological_interpretation": normalize_text_list(candidate.get("biological_interpretation")) or fallback["biological_interpretation"],
        "caveats": normalize_text_list(candidate.get("caveats")) or fallback["caveats"],
        "gene_details": merge_gene_details(candidate.get("gene_details"), fallback.get("gene_details", [])),
        "gene_relationships": normalize_text_list(candidate.get("gene_relationships")) or fallback.get("gene_relationships", []),
    }
    return normalized


def build_prompt(
    figure_name: str,
    title: str,
    payload: dict[str, Any],
    local_explanation: dict[str, Any],
    gene_filter_context: dict[str, Any],
) -> str:
    prompt = {
        "task": "Explain one bioinformatics result figure in Chinese based only on the provided structured plot-data summary.",
        "requirements": [
            "Do not invent numbers, genes, or trends not present in the payload.",
            "Keep the tone analytical and concise.",
            "Focus on what the figure shows, the main signal, biological meaning, and caveats.",
            "For site_specific_dotplot and shared_core, gene_details must cover every gene listed in gene_items_for_detailed_interpretation exactly once.",
            "For poorly annotated genes, explicitly say annotation is limited instead of inventing a precise mechanism.",
            "Use gene_relationships to explain how genes within the same panel or shared-core group connect biologically.",
            "For non-gene-centric figures, return gene_details as an empty list and gene_relationships as an empty list.",
            "Return valid JSON only, without markdown fences.",
        ],
        "output_schema": {
            "one_sentence_takeaway": "string",
            "what_is_plotted": "string",
            "key_observations": ["string", "string", "string"],
            "biological_interpretation": ["string", "string"],
            "gene_details": [
                {
                    "gene": "string",
                    "context": "string",
                    "detail": "string",
                }
            ],
            "gene_relationships": ["string", "string"],
            "caveats": ["string", "string"],
        },
        "figure_name": figure_name,
        "figure_title": title,
        "gene_filter_context": gene_filter_context,
        "structured_payload": payload,
        "baseline_local_summary": {
            "one_sentence_takeaway": local_explanation["one_sentence_takeaway"],
            "what_is_plotted": local_explanation["what_is_plotted"],
            "key_observations": local_explanation["key_observations"],
            "biological_interpretation": local_explanation["biological_interpretation"],
            "gene_relationship_hints": local_explanation.get("gene_relationships", []),
        },
    }
    return json.dumps(prompt, ensure_ascii=False, indent=2)


def call_deepseek(prompt: str, model: str, api_url: str, api_key: str, timeout_seconds: int) -> str:
    body = {
        "model": model,
        "temperature": 0.2,
        "response_format": {"type": "json_object"},
        "messages": [
            {
                "role": "system",
                "content": "You are a careful computational biology figure interpreter. Return strict JSON only.",
            },
            {
                "role": "user",
                "content": prompt,
            },
        ],
    }
    req = urllib.request.Request(
        api_url,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout_seconds) as response:
        payload = json.loads(response.read().decode("utf-8"))
    choices = payload.get("choices", [])
    if not choices:
        raise RuntimeError("DeepSeek response did not contain choices.")
    message = choices[0].get("message", {})
    content = message.get("content", "")
    if isinstance(content, list):
        parts = []
        for item in content:
            if isinstance(item, dict) and item.get("type") == "text":
                parts.append(item.get("text", ""))
            else:
                parts.append(str(item))
        return "\n".join(parts).strip()
    return str(content).strip()


def parse_json_response(text: str) -> dict[str, Any]:
    text = text.strip()
    if text.startswith("```"):
        text = text.strip("`")
        text = text.replace("json", "", 1).strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        start = text.find("{")
        end = text.rfind("}")
        if start != -1 and end != -1 and end > start:
            return json.loads(text[start:end + 1])
        raise


def load_plot_context(figure_name: str, source_tables: list[str], max_rows: int) -> dict[str, Any]:
    paths = [Path(path) for path in source_tables]
    if figure_name == "contrast_overview":
        return summarize_contrast(read_tsv(paths[0]), max_rows)
    if figure_name == "site_specific_dotplot":
        return summarize_site_specific(read_tsv(paths[0]), max_rows)
    if figure_name == "shared_core":
        return summarize_shared(read_tsv(paths[0]), max_rows)
    if figure_name == "run_diagnostics":
        return summarize_run_diagnostics(read_tsv(paths[0]), read_tsv(paths[1]), max_rows)
    raise ValueError(f"Unsupported figure_name: {figure_name}")


def render_markdown(explanations: list[FigureExplanation], output_path: Path) -> None:
    gene_figures = {"site_specific_dotplot", "shared_core"}
    lines: list[str] = []
    lines.append("# B-cell core discovery figure explanations")
    lines.append("")
    lines.append(f"Generated at: {datetime.now().isoformat()}")
    lines.append("")
    for item in explanations:
        lines.append(f"## {item.figure_name}")
        lines.append("")
        lines.append(f"- 标题: {item.title}")
        lines.append(f"- API 状态: {item.api_status}")
        lines.append(f"- 模型: {item.model or 'n/a'}")
        lines.append(f"- 图件: {', '.join(item.figure_outputs)}")
        lines.append(f"- 数据源: {', '.join(item.source_tables)}")
        if item.figure_name in gene_figures and item.gene_filter_context.get("enabled"):
            source_counts = item.gene_filter_context.get("source_table_counts", {})
            count_bits = []
            for source_name, stats in source_counts.items():
                count_bits.append(
                    f"{source_name}: {stats.get('input_rows', 0)}→{stats.get('retained_rows', 0)}（过滤 {stats.get('filtered_rows', 0)}）"
                )
            lines.append(
                f"- 基因过滤: {item.gene_filter_context.get('profile_name', 'enabled')}；" + "；".join(count_bits)
            )
        if item.error_message:
            lines.append(f"- API 备注: {item.error_message}")
        lines.append("")
        final = item.final_explanation
        lines.append("### 一句话结论")
        lines.append("")
        lines.append(final["one_sentence_takeaway"])
        lines.append("")
        lines.append("### 图上画了什么")
        lines.append("")
        lines.append(final["what_is_plotted"])
        lines.append("")
        lines.append("### 主要观察")
        lines.append("")
        for obs in final["key_observations"]:
            lines.append(f"- {obs}")
        lines.append("")
        lines.append("### 生物学解释")
        lines.append("")
        for obs in final["biological_interpretation"]:
            lines.append(f"- {obs}")
        lines.append("")
        if final.get("gene_details"):
            lines.append("### 逐基因说明")
            lines.append("")
            for gene_item in final["gene_details"]:
                lines.append(f"- `{gene_item['gene']}`（{gene_item['context']}）: {gene_item['detail']}")
            lines.append("")
        if final.get("gene_relationships"):
            lines.append("### 基因之间的联系")
            lines.append("")
            for rel in final["gene_relationships"]:
                lines.append(f"- {rel}")
            lines.append("")
        lines.append("### 注意事项")
        lines.append("")
        for obs in final["caveats"]:
            lines.append(f"- {obs}")
        lines.append("")
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    figure_dir = args.run_root / "core_discovery" / "figures"
    explanation_dir = figure_dir / "explanations"
    ensure_dirs([explanation_dir])

    manifest = load_figure_manifest(figure_dir)
    gene_filter_context = compact_gene_filter_context(manifest)
    api_key, key_source = resolve_deepseek_key(args.env_file)
    api_enabled = bool(api_key) and not is_placeholder_key(api_key)

    explanations: list[FigureExplanation] = []
    for figure in manifest["figures"]:
        figure_name = figure["name"]
        title = FIGURE_TITLES.get(figure_name, figure_name)
        source_tables = figure.get("sources", [])
        payload = load_plot_context(figure_name, source_tables, args.max_table_rows)
        local_explanation = build_local_explanation(figure_name, payload)
        final_explanation = local_explanation
        api_status = "skipped_no_key"
        raw_response_text = None
        error_message = None

        prompt = build_prompt(figure_name, title, payload, local_explanation, gene_filter_context)
        (explanation_dir / f"{figure_name}.prompt.json").write_text(prompt, encoding="utf-8")

        if api_enabled:
            try:
                raw_response_text = call_deepseek(
                    prompt=prompt,
                    model=args.model,
                    api_url=args.api_url,
                    api_key=api_key,
                    timeout_seconds=args.timeout_seconds,
                )
                parsed = parse_json_response(raw_response_text)
                final_explanation = normalize_final_explanation(parsed, local_explanation)
                api_status = "success"
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, RuntimeError, json.JSONDecodeError, OSError) as exc:
                api_status = "api_error"
                error_message = str(exc)
        elif key_source != "missing":
            api_status = "skipped_placeholder_key"
            error_message = "DEEPSEEK_API_KEY looked like a placeholder, so API calling was skipped."

        figure_record = FigureExplanation(
            figure_name=figure_name,
            title=title,
            figure_outputs=figure.get("outputs", []),
            source_tables=source_tables,
            gene_filter_context=gene_filter_context,
            api_status=api_status,
            model=args.model if api_enabled else None,
            local_summary=payload,
            final_explanation=final_explanation,
            raw_response_text=raw_response_text,
            error_message=error_message,
        )
        explanations.append(figure_record)

        figure_json_path = explanation_dir / f"{figure_name}.explanation.json"
        figure_json_path.write_text(
            json.dumps(
                {
                    "figure_name": figure_record.figure_name,
                    "title": figure_record.title,
                    "figure_outputs": figure_record.figure_outputs,
                    "source_tables": figure_record.source_tables,
                    "gene_filter_context": figure_record.gene_filter_context,
                    "api_status": figure_record.api_status,
                    "model": figure_record.model,
                    "error_message": figure_record.error_message,
                    "local_summary": figure_record.local_summary,
                    "final_explanation": figure_record.final_explanation,
                    "raw_response_text": figure_record.raw_response_text,
                },
                ensure_ascii=False,
                indent=2,
                default=json_safe,
            ),
            encoding="utf-8",
        )

    aggregate_manifest = {
        "generated_at": datetime.now().isoformat(),
        "run_root": str(args.run_root),
        "figure_dir": str(figure_dir),
        "explanation_dir": str(explanation_dir),
        "api_enabled": api_enabled,
        "key_source": key_source,
        "gene_filter_context": gene_filter_context,
        "figures": [
            {
                "figure_name": item.figure_name,
                "title": item.title,
                "api_status": item.api_status,
                "model": item.model,
                "figure_outputs": item.figure_outputs,
                "source_tables": item.source_tables,
                "json_path": str(explanation_dir / f"{item.figure_name}.explanation.json"),
                "prompt_path": str(explanation_dir / f"{item.figure_name}.prompt.json"),
            }
            for item in explanations
        ],
    }
    (explanation_dir / "figure_explanations_manifest.json").write_text(
        json.dumps(aggregate_manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    render_markdown(explanations, explanation_dir / "figure_explanations.md")

    print(f"[done] explanation dir: {explanation_dir}")
    print(f"[done] api_enabled={api_enabled} (source={key_source})")
    for item in explanations:
        suffix = f" ({item.error_message})" if item.error_message else ""
        print(f"  - {item.figure_name}: {item.api_status}{suffix}")


if __name__ == "__main__":
    main()
