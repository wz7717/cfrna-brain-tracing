#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
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
from scripts.run_bo2023_loso_validation import read_vsd_matrix, softmax  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import (  # noqa: E402
    build_group_reference,
    select_group_discriminative_genes,
)
from scripts.run_bo2023_network_pairwise_correlation_validation import (  # noqa: E402
    BASELINE_ROUTE,
    PAIR_TOP3_ROUTE,
    build_pair_models,
    derive_training_confusion_pairs,
    evaluate_pairwise_rescue,
    make_detail,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_network_regularized_centroid_screen_loso_819_20260604"
WEIGHTED_ROUTE = "network_fisher_weighted_correlation_top200"
SHRINK_PREFIX = "network_shrinkage_centroid"
NSC_PREFIX = "network_nearest_shrunken_centroid"


def weighted_correlation_scores(reference: np.ndarray, sample: np.ndarray, genes: np.ndarray, weights: np.ndarray) -> np.ndarray:
    ref = reference[genes, :].astype(float, copy=False)
    vec = sample[genes].astype(float, copy=False)
    w = np.asarray(weights, dtype=float)
    w = np.clip(w, 0.0, None)
    if float(w.sum()) <= 1e-12:
        w = np.ones_like(w, dtype=float)
    w = w / (float(w.sum()) + 1e-12)
    ref_mean = (ref * w[:, None]).sum(axis=0, keepdims=True)
    vec_mean = float((vec * w).sum())
    ref0 = ref - ref_mean
    vec0 = vec - vec_mean
    numerator = (w[:, None] * ref0 * vec0[:, None]).sum(axis=0)
    denominator = np.sqrt((w[:, None] * np.square(ref0)).sum(axis=0) * float((w * np.square(vec0)).sum()) + 1e-12)
    return np.nan_to_num(numerator / denominator, nan=0.0, posinf=0.0, neginf=0.0)


def unweighted_correlation_scores(reference: np.ndarray, sample: np.ndarray, genes: np.ndarray) -> np.ndarray:
    return weighted_correlation_scores(reference, sample, genes, np.ones(len(genes), dtype=float))


def build_shrinkage_reference(reference: np.ndarray, training: dict[str, np.ndarray], values: np.ndarray, groups: list[str], alpha: float) -> np.ndarray:
    all_train = np.concatenate([training[group] for group in groups])
    global_mean = values[:, all_train].mean(axis=1, dtype=np.float64)
    return ((1.0 - float(alpha)) * reference + float(alpha) * global_mean[:, None]).astype(np.float32)


def build_nsc_reference(
    reference: np.ndarray,
    training: dict[str, np.ndarray],
    values: np.ndarray,
    groups: list[str],
    threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    all_train = np.concatenate([training[group] for group in groups])
    global_mean = values[:, all_train].mean(axis=1, dtype=np.float64)
    within_sum = np.zeros(values.shape[0], dtype=float)
    within_denom = 0
    for group in groups:
        x = values[:, training[group]].astype(float, copy=False)
        mean = x.mean(axis=1)
        within_sum += np.square(x - mean[:, None]).sum(axis=1)
        within_denom += max(0, x.shape[1] - 1)
    pooled_sd = np.sqrt(within_sum / max(1, within_denom) + 1e-8)
    standardized_delta = (reference - global_mean[:, None]) / pooled_sd[:, None]
    shrunk_delta = np.sign(standardized_delta) * np.maximum(np.abs(standardized_delta) - float(threshold), 0.0)
    shrunk_reference = global_mean[:, None] + shrunk_delta * pooled_sd[:, None]
    active = np.any(np.abs(shrunk_delta) > 0.0, axis=1)
    return shrunk_reference.astype(np.float32), active


def detail_from_scores(route: str, sample_id: str, truth: str, groups: list[str], scores: np.ndarray) -> tuple[dict[str, Any], dict[str, Any]]:
    order = np.argsort(scores)[::-1].tolist()
    return make_detail(route, sample_id, truth, groups, order, scores)


def summarize(detail: pd.DataFrame, probabilities: pd.DataFrame, groups: list[str]) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probabilities[groups])),
        "median_true_rank": float(detail["true_rank"].median()),
        "n_switches": int(detail.get("switched", pd.Series(0, index=detail.index)).sum()),
    }


def paired_change(base: pd.DataFrame, tested: pd.DataFrame) -> dict[str, int]:
    return {
        "top1_gains": int(((base["hit1"] == 0) & (tested["hit1"] == 1)).sum()),
        "top1_losses": int(((base["hit1"] == 1) & (tested["hit1"] == 0)).sum()),
        "top3_gains": int(((base["hit3"] == 0) & (tested["hit3"] == 1)).sum()),
        "top3_losses": int(((base["hit3"] == 1) & (tested["hit3"] == 0)).sum()),
    }


def paired_pvalue(gains: int, losses: int) -> float:
    n = gains + losses
    if n == 0:
        return 1.0
    tail = min(gains, losses)
    probability = sum(math.comb(n, k) for k in range(tail + 1)) / (2**n)
    return float(min(1.0, 2.0 * probability))


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict LOSO screen of regularized Network centroid scorers.")
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
    parser.add_argument("--pair-min-margin", type=float, default=0.0)
    parser.add_argument("--shrinkage-alphas", default="0.05,0.10,0.20")
    parser.add_argument("--nsc-thresholds", default="0.25,0.50,0.75,1.00")
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    shrinkage_alphas = [float(x) for x in str(args.shrinkage_alphas).split(",") if str(x).strip()]
    nsc_thresholds = [float(x) for x in str(args.nsc_thresholds).split(",") if str(x).strip()]

    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet)
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.endpoint].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    samples = matrix.columns.astype(str).tolist()
    if args.max_samples is not None:
        samples = samples[: max(1, int(args.max_samples))]
    all_samples = matrix.columns.astype(str).tolist()
    labels = ann.set_index("sample_id").reindex(all_samples)["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    values = matrix.to_numpy(dtype=np.float32)

    routes = [BASELINE_ROUTE, PAIR_TOP3_ROUTE, WEIGHTED_ROUTE]
    routes.extend([f"{SHRINK_PREFIX}_alpha_{str(alpha).replace('.', 'p')}" for alpha in shrinkage_alphas])
    routes.extend([f"{NSC_PREFIX}_threshold_{str(thr).replace('.', 'p')}" for thr in nsc_thresholds])
    route_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in routes}
    probability_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in routes}
    fold_rows: list[dict[str, Any]] = []

    for fold_no, sample_id in enumerate(samples, start=1):
        heldout_idx = all_samples.index(sample_id)
        truth = str(labels[heldout_idx])
        reference, training = build_group_reference(values, labels, groups, heldout_idx)
        gene_pool, gene_audit = select_group_discriminative_genes(values, groups, training, args.gene_pool_size)
        global_genes = gene_pool[: args.global_top_n_genes]
        fisher_weights = gene_audit.head(len(global_genes))["fisher_score"].to_numpy(dtype=float)
        fisher_weights = np.sqrt(np.maximum(fisher_weights, 0.0))
        sample = values[:, heldout_idx]

        base_scores = unweighted_correlation_scores(reference, sample, global_genes)
        detail, prob = detail_from_scores(BASELINE_ROUTE, sample_id, truth, groups, base_scores)
        route_rows[BASELINE_ROUTE].append(detail)
        probability_rows[BASELINE_ROUTE].append(prob)

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
        detail, prob = evaluate_pairwise_rescue(
            PAIR_TOP3_ROUTE,
            sample_id,
            truth,
            sample,
            reference,
            groups,
            base_scores,
            pair_models,
            3,
            args.pair_min_margin,
        )
        route_rows[PAIR_TOP3_ROUTE].append(detail)
        probability_rows[PAIR_TOP3_ROUTE].append(prob)

        weighted_scores = weighted_correlation_scores(reference, sample, global_genes, fisher_weights)
        detail, prob = detail_from_scores(WEIGHTED_ROUTE, sample_id, truth, groups, weighted_scores)
        route_rows[WEIGHTED_ROUTE].append(detail)
        probability_rows[WEIGHTED_ROUTE].append(prob)

        for alpha in shrinkage_alphas:
            route = f"{SHRINK_PREFIX}_alpha_{str(alpha).replace('.', 'p')}"
            shrunk_reference = build_shrinkage_reference(reference, training, values, groups, alpha)
            scores = unweighted_correlation_scores(shrunk_reference, sample, global_genes)
            detail, prob = detail_from_scores(route, sample_id, truth, groups, scores)
            route_rows[route].append(detail)
            probability_rows[route].append(prob)

        for threshold in nsc_thresholds:
            route = f"{NSC_PREFIX}_threshold_{str(threshold).replace('.', 'p')}"
            nsc_reference, active = build_nsc_reference(reference, training, values, groups, threshold)
            active_genes = global_genes[active[global_genes]]
            genes = active_genes if len(active_genes) >= 20 else global_genes
            scores = unweighted_correlation_scores(nsc_reference, sample, genes)
            detail, prob = detail_from_scores(route, sample_id, truth, groups, scores)
            route_rows[route].append(detail)
            probability_rows[route].append(prob)

        fold_rows.append({"fold": fold_no, "sample_id": sample_id, "truth": truth, "n_pairs": int(len(pairs))})

    details = {route: pd.DataFrame(rows) for route, rows in route_rows.items()}
    probabilities = {route: pd.DataFrame(rows) for route, rows in probability_rows.items()}
    metrics = {route: summarize(details[route], probabilities[route], groups) for route in routes}
    changes = {route: paired_change(details[PAIR_TOP3_ROUTE], details[route]) for route in routes if route != PAIR_TOP3_ROUTE}
    pvalues = {
        route: {
            "top1": paired_pvalue(change["top1_gains"], change["top1_losses"]),
            "top3": paired_pvalue(change["top3_gains"], change["top3_losses"]),
        }
        for route, change in changes.items()
    }
    best_route = max(routes, key=lambda route: (metrics[route]["top1_accuracy"], metrics[route]["top3_accuracy"]))
    base = metrics[PAIR_TOP3_ROUTE]
    best = metrics[best_route]
    improved = best_route != PAIR_TOP3_ROUTE and best["top1_accuracy"] > base["top1_accuracy"] and best["top3_accuracy"] >= base["top3_accuracy"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "endpoint": args.endpoint,
        "validation_design": "strict outer LOSO regularized centroid screen; fold-local gene selection and fold-local confusion pairs",
        "n_test_samples": int(len(samples)),
        "n_classes": int(len(groups)),
        "global_top_n_genes": int(args.global_top_n_genes),
        "gene_pool_size": int(args.gene_pool_size),
        "shrinkage_alphas": shrinkage_alphas,
        "nsc_thresholds": nsc_thresholds,
        "baseline_route": PAIR_TOP3_ROUTE,
        "routes": metrics,
        "paired_changes_vs_pairwise_top3": changes,
        "paired_pvalues_vs_pairwise_top3": pvalues,
        "best_route": best_route,
        "meets_internal_adoption_rule": bool(improved),
        "decision": (
            f"{best_route} improves Network Top1 without lowering Top3 versus pairwise Top3; keep as a candidate for independent confirmation."
            if improved
            else "No regularized centroid route improves Network Top1 under the no-Top3-loss constraint; keep pairwise Top3 as formal route."
        ),
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, frame in details.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_regularized_centroid_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
