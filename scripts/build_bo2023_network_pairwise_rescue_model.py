#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_bo2023_network_correlation_validation import (  # noqa: E402
    build_group_reference,
    select_group_discriminative_genes,
)
from scripts.run_bo2023_network_pairwise_correlation_validation import (  # noqa: E402
    build_pair_models,
    derive_training_confusion_pairs,
    pair_key,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)
from scripts.run_bo2023_loso_validation import read_vsd_matrix  # noqa: E402


DEFAULT_OUT = ROOT / "data" / "models" / "bo2023_saleem_network_pairwise_rescue_model.json"
DEFAULT_VALIDATION = (
    ROOT
    / "results"
    / "bo2023_network_pairwise_correlation_full_loso_819_rerun_20260526"
    / "validation_summary.json"
)


def full_group_reference(values: np.ndarray, labels: np.ndarray, groups: list[str]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    columns: list[np.ndarray] = []
    training: dict[str, np.ndarray] = {}
    for group in groups:
        idx = np.flatnonzero(labels == group)
        if not len(idx):
            raise ValueError(f"group {group} has no samples")
        training[group] = idx
        columns.append(values[:, idx].mean(axis=1, dtype=np.float64))
    return np.column_stack(columns), training


def main() -> int:
    parser = argparse.ArgumentParser(description="Build full-training Bo2023 Network pairwise rescue model.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--endpoint", default="SaleemNetworks")
    parser.add_argument("--global-top-n-genes", type=int, default=200)
    parser.add_argument("--gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=100)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=3)
    parser.add_argument("--pair-min-margin", type=float, default=0.002)
    parser.add_argument("--validation-summary", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet)
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.endpoint].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    samples = matrix.columns.astype(str).tolist()
    labels = ann.set_index("sample_id").reindex(samples)["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    values = matrix.to_numpy(dtype=np.float32)
    gene_symbols = matrix.index.astype(str).to_numpy()

    reference, training = full_group_reference(values, labels, groups)
    gene_pool, gene_audit = select_group_discriminative_genes(values, groups, training, args.gene_pool_size)
    global_genes = gene_pool[: args.global_top_n_genes]
    pairs, pair_errors = derive_training_confusion_pairs(
        values,
        labels,
        groups,
        training,
        reference,
        global_genes,
        args.max_pairs_per_truth,
        args.min_pair_errors,
    )
    pair_models, pair_audit = build_pair_models(
        values,
        training,
        reference,
        pairs,
        groups,
        gene_pool,
        args.pair_top_n_genes,
    )

    group_pos = {group: i for i, group in enumerate(groups)}
    pair_entries: list[dict[str, Any]] = []
    selected_counts = pair_errors[pair_errors["selected"]].copy()
    for left, right in sorted(pair_models):
        genes = pair_models[(left, right)]
        ref = reference[genes, :][:, [group_pos[left], group_pos[right]]]
        left_errors = int(
            selected_counts.loc[
                (selected_counts["truth_network"] == left)
                & (selected_counts["confused_as_network"] == right),
                "error_count",
            ].sum()
        )
        right_errors = int(
            selected_counts.loc[
                (selected_counts["truth_network"] == right)
                & (selected_counts["confused_as_network"] == left),
                "error_count",
            ].sum()
        )
        pair_entries.append(
            {
                "left_network": left,
                "right_network": right,
                "key": "||".join(pair_key(left, right)),
                "genes": gene_symbols[genes].astype(str).tolist(),
                "reference": {
                    left: ref[:, 0].astype(float).tolist(),
                    right: ref[:, 1].astype(float).tolist(),
                },
                "training_error_count": int(left_errors + right_errors),
            }
        )

    validation = {}
    if args.validation_summary.exists():
        validation = json.loads(args.validation_summary.read_text(encoding="utf-8"))

    model = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "endpoint": args.endpoint,
        "method": "top3_pairwise_rescue_over_discriminative_network_correlation",
        "parameters": {
            "global_top_n_genes": int(args.global_top_n_genes),
            "gene_pool_size": int(args.gene_pool_size),
            "pair_top_n_genes": int(args.pair_top_n_genes),
            "max_pairs_per_truth": int(args.max_pairs_per_truth),
            "min_pair_errors": int(args.min_pair_errors),
            "pair_min_margin": float(args.pair_min_margin),
            "candidate_k": 3,
        },
        "networks": groups,
        "n_training_samples": int(len(samples)),
        "n_pairs": int(len(pair_entries)),
        "pairs": pair_entries,
        "validation": validation,
        "notes": (
            "Pairwise rescue is applied only within the baseline Top3 Network beam. "
            "It may switch Top1 but does not change the Top3 candidate set."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(model, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"out": str(args.out), "n_pairs": len(pair_entries)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
