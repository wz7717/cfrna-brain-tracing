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

from benchmark.metrics import compute_multiclass_auc  # noqa: E402
from scripts.run_bo2023_lomo_network_normalization_screen import (  # noqa: E402
    apply_pairwise_top3,
    make_row,
    normalize_fold,
)
from scripts.run_bo2023_loso_validation import correlation_scores, read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import (  # noqa: E402
    build_group_reference,
    select_group_discriminative_genes,
)
from scripts.run_bo2023_network_pairwise_correlation_validation import (  # noqa: E402
    build_pair_models,
    derive_training_confusion_pairs,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_network_normalization_full_loso_819_20260604"


def summarize(detail: pd.DataFrame, probabilities: pd.DataFrame, groups: list[str]) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probabilities[groups])),
        "median_true_rank": float(detail["true_rank"].median()),
        "n_switches": int(detail["switched"].sum()),
    }


def paired_changes(base: pd.DataFrame, tested: pd.DataFrame) -> dict[str, int]:
    return {
        "top1_gains": int(((base["hit1"] == 0) & (tested["hit1"] == 1)).sum()),
        "top1_losses": int(((base["hit1"] == 1) & (tested["hit1"] == 0)).sum()),
        "top3_gains": int(((base["hit3"] == 0) & (tested["hit3"] == 1)).sum()),
        "top3_losses": int(((base["hit3"] == 1) & (tested["hit3"] == 0)).sum()),
    }


def paired_pvalue(gains: int, losses: int) -> float:
    import math

    n = gains + losses
    if n == 0:
        return 1.0
    tail = min(gains, losses)
    probability = sum(math.comb(n, k) for k in range(tail + 1)) / (2**n)
    return float(min(1.0, 2.0 * probability))


def main() -> int:
    parser = argparse.ArgumentParser(description="Full LOSO Network normalization screen for Bo2023.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--variants", nargs="+", default=["raw", "sample_zscore", "sample_rank_center"])
    parser.add_argument("--global-top-n-genes", type=int, default=200)
    parser.add_argument("--gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=100)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=3)
    parser.add_argument("--pair-min-margin", type=float, default=0.0)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    matrix = map_matrix_to_symbols(read_vsd_matrix(args.matrix), args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet, usecols=["No.", args.network_col])
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.network_col].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    samples = matrix.columns.astype(str).tolist()
    if args.max_samples is not None:
        samples = samples[: max(1, int(args.max_samples))]
    full_samples = matrix.columns.astype(str).tolist()
    sample_ann = ann.set_index("sample_id").reindex(full_samples)
    labels = sample_ann["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    raw_values = matrix.to_numpy(dtype=np.float32)

    detail_frames: dict[str, pd.DataFrame] = {}
    probability_frames: dict[str, pd.DataFrame] = {}
    for variant in args.variants:
        rows: list[dict[str, Any]] = []
        probs: list[dict[str, Any]] = []
        for fold_no, sample_id in enumerate(samples, start=1):
            heldout_idx = full_samples.index(sample_id)
            train_indices = np.asarray([j for j in range(len(full_samples)) if j != heldout_idx], dtype=int)
            values = normalize_fold(raw_values, train_indices, variant)
            reference, training = build_group_reference(values, labels, groups, heldout_idx)
            gene_pool, _ = select_group_discriminative_genes(values, groups, training, args.gene_pool_size)
            global_genes = gene_pool[: args.global_top_n_genes]
            pairs, _ = derive_training_confusion_pairs(
                values,
                labels,
                groups,
                training,
                reference,
                global_genes,
                args.max_pairs_per_truth,
                args.min_pair_errors,
            )
            pair_models, _ = build_pair_models(values, training, reference, pairs, groups, gene_pool, args.pair_top_n_genes)
            truth = str(labels[heldout_idx])
            scores = correlation_scores(reference, values[:, heldout_idx], global_genes)
            order, switched, switch_pair, switch_margin = apply_pairwise_top3(
                values,
                heldout_idx,
                reference,
                groups,
                scores,
                pair_models,
                args.pair_min_margin,
            )
            row, prob = make_row(
                variant,
                sample_id,
                "LOSO",
                truth,
                groups,
                order,
                scores,
                switched,
                switch_pair,
                switch_margin,
            )
            rows.append(row)
            probs.append(prob)
        detail_frames[variant] = pd.DataFrame(rows)
        probability_frames[variant] = pd.DataFrame(probs)

    metrics = {variant: summarize(detail_frames[variant], probability_frames[variant], groups) for variant in args.variants}
    base_variant = "raw" if "raw" in detail_frames else args.variants[0]
    changes = {variant: paired_changes(detail_frames[base_variant], detail_frames[variant]) for variant in args.variants if variant != base_variant}
    pvalues = {
        variant: {
            "top1": paired_pvalue(change["top1_gains"], change["top1_losses"]),
            "top3": paired_pvalue(change["top3_gains"], change["top3_losses"]),
        }
        for variant, change in changes.items()
    }
    best_variant = max(args.variants, key=lambda variant: (metrics[variant]["top1_accuracy"], metrics[variant]["top3_accuracy"]))
    base = metrics[base_variant]
    best = metrics[best_variant]
    improved = best["top1_accuracy"] > base["top1_accuracy"] and best["top3_accuracy"] >= base["top3_accuracy"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "full strict LOSO Network normalization screen with fold-local genes and fold-local pairwise Top3 rescue",
        "n_test_samples": int(len(samples)),
        "n_classes": int(len(groups)),
        "variants": args.variants,
        "baseline_variant": base_variant,
        "routes": metrics,
        "paired_changes_vs_baseline": changes,
        "paired_pvalues_vs_baseline": pvalues,
        "best_variant": best_variant,
        "meets_internal_adoption_rule": bool(improved),
        "decision": (
            f"{best_variant} improves Network Top1 without lowering Top3; treat as a candidate normalization route."
            if improved
            else "No normalization variant improves Network Top1 under the no-Top3-loss constraint; keep raw VSD route."
        ),
    }
    pd.DataFrame([{"variant": variant, **metric} for variant, metric in metrics.items()]).to_csv(
        args.outdir / "network_normalization_loso_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for variant, detail in detail_frames.items():
        detail.to_csv(args.outdir / f"{variant}_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
