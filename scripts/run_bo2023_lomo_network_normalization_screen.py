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
from scripts.run_bo2023_loso_validation import correlation_scores, read_vsd_matrix, softmax  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_network_pairwise_correlation_validation import (  # noqa: E402
    pair_key,
    derive_training_confusion_pairs,
    build_pair_models,
)
from scripts.run_bo2023_v2_loso_validation import DEFAULT_GENE_MAP, DEFAULT_MATRIX, DEFAULT_SAMPLE_INFO, map_matrix_to_symbols  # noqa: E402


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_lomo_network_normalization_screen_20260601"
VARIANTS = ["raw", "train_gene_zscore", "train_gene_robust_zscore", "sample_zscore", "sample_rank_center"]


def normalize_fold(raw_values: np.ndarray, train_indices: np.ndarray, variant: str) -> np.ndarray:
    values = raw_values.astype(np.float32, copy=True)
    if variant == "raw":
        return values
    if variant == "train_gene_zscore":
        mean = values[:, train_indices].mean(axis=1, keepdims=True, dtype=np.float64)
        std = values[:, train_indices].std(axis=1, keepdims=True, dtype=np.float64)
        return ((values - mean) / (std + 1e-6)).astype(np.float32)
    if variant == "train_gene_robust_zscore":
        train = values[:, train_indices]
        median = np.median(train, axis=1, keepdims=True)
        q75 = np.percentile(train, 75, axis=1, keepdims=True)
        q25 = np.percentile(train, 25, axis=1, keepdims=True)
        iqr = q75 - q25
        return ((values - median) / (iqr + 1e-6)).astype(np.float32)
    if variant == "sample_zscore":
        mean = values.mean(axis=0, keepdims=True, dtype=np.float64)
        std = values.std(axis=0, keepdims=True, dtype=np.float64)
        return ((values - mean) / (std + 1e-6)).astype(np.float32)
    if variant == "sample_rank_center":
        ranked = np.empty_like(values, dtype=np.float32)
        denom = max(1, values.shape[0] - 1)
        for j in range(values.shape[1]):
            order = np.argsort(values[:, j], kind="mergesort")
            ranks = np.empty(values.shape[0], dtype=np.float32)
            ranks[order] = np.arange(values.shape[0], dtype=np.float32)
            ranked[:, j] = (ranks / denom) - 0.5
        return ranked
    raise ValueError(f"Unknown normalization variant: {variant}")


def build_reference(values: np.ndarray, labels: np.ndarray, groups: list[str], train_indices: np.ndarray) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reference = np.zeros((values.shape[0], len(groups)), dtype=np.float32)
    training: dict[str, np.ndarray] = {}
    for j, group in enumerate(groups):
        idx = train_indices[labels[train_indices] == group]
        if idx.size == 0:
            raise ValueError(f"No training samples for label {group!r}")
        training[group] = idx
        reference[:, j] = values[:, idx].mean(axis=1, dtype=np.float64).astype(np.float32)
    return reference, training


def make_row(
    variant: str,
    sample_id: str,
    monkey_id: str,
    truth: str,
    groups: list[str],
    order: list[int],
    scores: np.ndarray,
    switched: bool,
    switch_pair: str,
    switch_margin: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    ranked = [groups[i] for i in order]
    true_rank = ranked.index(truth) + 1
    probs = softmax(scores)
    row = {
        "variant": variant,
        "sample_id": sample_id,
        "monkey_id": monkey_id,
        "label": truth,
        "pred_top1": ranked[0],
        "pred_top2": ranked[1],
        "pred_top3": ranked[2],
        "true_rank": int(true_rank),
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "switched": int(switched),
        "switch_pair": switch_pair,
        "switch_margin": float(switch_margin),
    }
    prob = {"variant": variant, "sample_id": sample_id, "monkey_id": monkey_id, "label": truth}
    prob.update({group: float(probs[i]) for i, group in enumerate(groups)})
    return row, prob


def apply_pairwise_top3(
    values: np.ndarray,
    sample_idx: int,
    reference: np.ndarray,
    groups: list[str],
    scores: np.ndarray,
    pair_models: dict[tuple[str, str], np.ndarray],
    min_margin: float,
) -> tuple[list[int], bool, str, float]:
    order = np.argsort(scores)[::-1].tolist()
    anchor = order[0]
    best_position = -1
    best_margin = float(min_margin)
    best_pair = ""
    for position in range(1, min(3, len(order))):
        challenger = order[position]
        key = pair_key(groups[anchor], groups[challenger])
        genes = pair_models.get(key)
        if genes is None:
            continue
        pair_scores = correlation_scores(reference, values[:, sample_idx], genes)
        margin = float(pair_scores[challenger] - pair_scores[anchor])
        if margin > best_margin:
            best_position = position
            best_margin = margin
            best_pair = f"{key[0]} <> {key[1]}"
    if best_position >= 1:
        order[0], order[best_position] = order[best_position], order[0]
        return order, True, best_pair, best_margin
    return order, False, "", 0.0


def summarize(detail: pd.DataFrame, probabilities: pd.DataFrame, groups: list[str]) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_by_monkey_top1": float(detail.groupby("monkey_id")["hit1"].mean().mean()),
        "macro_by_monkey_top3": float(detail.groupby("monkey_id")["hit3"].mean().mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probabilities[groups])),
        "median_true_rank": float(detail["true_rank"].median()),
        "n_switches": int(detail["switched"].sum()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="LOMO Network normalization screen for Bo2023.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--monkey-col", default="MonkeyID")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--variants", nargs="+", default=VARIANTS)
    parser.add_argument("--global-top-n-genes", type=int, default=200)
    parser.add_argument("--gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=100)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=3)
    parser.add_argument("--pair-min-margin", type=float, default=0.0)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = map_matrix_to_symbols(read_vsd_matrix(args.matrix), args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet, usecols=["No.", args.monkey_col, args.network_col])
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["monkey_id"] = ann[args.monkey_col].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.network_col].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    samples = matrix.columns.astype(str).tolist()
    sample_ann = ann.set_index("sample_id").reindex(samples)
    labels = sample_ann["endpoint_label"].to_numpy(dtype=str)
    monkeys = sample_ann["monkey_id"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    raw_values = matrix.to_numpy(dtype=np.float32)
    monkey_ids = sorted(set(monkeys))

    all_details: list[pd.DataFrame] = []
    all_probs: list[pd.DataFrame] = []
    metrics: dict[str, dict[str, Any]] = {}
    for variant in args.variants:
        rows: list[dict[str, Any]] = []
        probs: list[dict[str, Any]] = []
        for monkey_id in monkey_ids:
            test_indices = np.where(monkeys == monkey_id)[0]
            train_indices = np.where(monkeys != monkey_id)[0]
            values = normalize_fold(raw_values, train_indices, variant)
            reference, training = build_reference(values, labels, groups, train_indices)
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
            for sample_idx in test_indices:
                sample_id = samples[sample_idx]
                truth = str(labels[sample_idx])
                scores = correlation_scores(reference, values[:, sample_idx], global_genes)
                order, switched, switch_pair, switch_margin = apply_pairwise_top3(
                    values, sample_idx, reference, groups, scores, pair_models, args.pair_min_margin
                )
                row, prob = make_row(
                    variant,
                    sample_id,
                    str(monkey_id),
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
        detail = pd.DataFrame(rows)
        probability = pd.DataFrame(probs)
        metrics[variant] = summarize(detail, probability, groups)
        all_details.append(detail)
        all_probs.append(probability)

    detail_all = pd.concat(all_details, ignore_index=True)
    prob_all = pd.concat(all_probs, ignore_index=True)
    metric_frame = pd.DataFrame([{"variant": variant, **metric} for variant, metric in metrics.items()])
    metric_frame = metric_frame.sort_values(["top3_accuracy", "top1_accuracy"], ascending=[False, False])
    metric_frame.to_csv(args.outdir / "network_normalization_variant_metrics.csv", index=False, encoding="utf-8-sig")
    detail_all.to_csv(args.outdir / "network_normalization_variant_detail.csv", index=False, encoding="utf-8-sig")
    prob_all.to_csv(args.outdir / "network_normalization_variant_probabilities.csv", index=False, encoding="utf-8-sig")
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "leave-one-monkey-out Network-only normalization screen with fold-local pairwise Top3 rescue",
        "variants": args.variants,
        "best_by_top3": str(metric_frame.iloc[0]["variant"]),
        "routes": metrics,
    }
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
