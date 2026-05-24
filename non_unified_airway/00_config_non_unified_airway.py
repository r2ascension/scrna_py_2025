#!/usr/bin/env python3
"""Discover primary inputs for the non-unified airway workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PY_ROOT = SCRIPT_DIR.parent
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from non_unified_airway.common import (
    DEFAULT_METHODS_DOC,
    DEFAULT_REPO_ROOT,
    discover_input_inventory,
    write_json,
    write_tsv,
)

DEFAULT_OUT_DIR = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/config"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--methods-doc", type=Path, default=DEFAULT_METHODS_DOC)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def run_config_discovery(repo_root: Path, methods_doc: Path, out_dir: Path) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    inventory = discover_input_inventory(repo_root=repo_root, methods_doc_path=methods_doc)
    records_df = pd.DataFrame(inventory["records"])
    write_tsv(out_dir / "discovered_inputs.tsv", records_df)

    lineage_df = records_df.loc[records_df["role"] == "lineage_output_dir"].copy()
    write_tsv(out_dir / "discovered_lineage_output_dirs.tsv", lineage_df)

    reference_df = records_df.loc[records_df["role"] == "allcells_reference"].copy()
    write_tsv(out_dir / "discovered_allcells_reference_candidates.tsv", reference_df)

    query_df = records_df.loc[records_df["role"] == "mapped_query"].copy()
    write_tsv(out_dir / "discovered_mapped_query_candidates.tsv", query_df)

    write_json(out_dir / "config_manifest.json", inventory)
    return inventory


def main() -> int:
    args = parse_args()
    inventory = run_config_discovery(
        repo_root=args.repo_root,
        methods_doc=args.methods_doc,
        out_dir=args.out_dir,
    )
    print(f"Wrote config manifest with {len(inventory['records'])} discovered path records to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
