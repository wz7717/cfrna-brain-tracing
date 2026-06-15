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
from scripts.run_bo2023_lomo_network_normalization_screen import normalize_fold  # noqa: E402
from scripts.run_bo2023_loso_validation import correlation_scores, read_vsd_matrix, softmax  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import (  # noqa: E402
    build_group_reference,
    select_group_discriminative_genes,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_network_multiview_fusion_loso_819_20260604"
BASELINE_ROUTE = "raw_top200_pearson"


def zscore_scores(scores: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=float)
    std = float(scores.std())
    if std <= 1e-12:
        return np.zeros_like(scores, dtype=float)
    return (scores - float(scores.mean())) / std


def rank_fusion_score(score_list: list[np.ndarray], method: str) -> np.ndarray:
    n_classes = len(score_list[0])
    fused = np.zeros(n_classes, dtype=float)
    for scores in score_list:
        order = np.argsort(scores)[::-1]
        ranks = np.empty(n_classes, dtype=float)
        ranks[order] = np.arange(1, n_classes + 1, dtype=float)
        if method == "rrf":
            fused += 1.0 / (60.0 + ranks)
        elif method == "borda":
            fused += n_classes - ranks
        else:
            raise ValueError(f"unknown rank fusion method: {method}")
    return fused


def detail_from_scores(route: str, sample_id: str, truth: str, groups: list[str], scores: np.ndarray) -> tuple[dict[str, Any], dict[str, Any]]:
    order = np.argsort(scores)[::-1]
    ranked = [groups[int(j)] for j in order]
    true_rank = ranked.index(truth) + 1
    probs = softmax(scores)
    row = {
        "route": route,
        "sample_id": sample_id,
        "label": truth,
        "pred_top1": ranked[0],
        "pred_top2": ranked[1],
        "pred_top3": ranked[2],
        "true_rank": int(true_rank),
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
    }
    probability: dict[str, Any] = {"sample_id": sample_id, "label": truth}
    probability.update({group: float(probs[j]) for j, group in enumerate(groups)})
    return row, probability


def summarize(detail: pd.DataFrame, probabilities: pd.DataFrame, groups: list[str]) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probabilities[groups])),
        "median_true_rank": float(detail["true_rank"].median()),
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
    parser = argparse.ArgumentParser(description="Strict LOSO multiview fusion screen for Bo2023 Network tracing.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--endpoint", default="SaleemNetworks")
    parser.add_argument("--views", nargs="+", default=["raw", "sample_zscore", "sample_rank_center"])
    parser.add_argument("--gene-counts", default="50,100,200,500,1000")
    parser.add_argument("--fusion-routes", nargs="+", default=[
        "raw_top100_raw_top200_zmean",
        "raw_top50_raw_top100_raw_top200_zmean",
        "raw_top200_rank_top200_zmean",
        "raw_top200_zscore_top200_rank_top200_zmean",
        "raw_top200_zscore_top200_rank_top200_rrf",
        "raw_top200_zscore_top200_rank_top200_borda",
    ])
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    gene_counts = sorted({int(x) for x in str(args.gene_counts).split(",") if str(x).strip()})
    max_gene_count = max(gene_counts)
    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet)
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.endpoint].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    all_samples = matrix.columns.astype(str).tolist()
    samples = all_samples if args.max_samples is None else all_samples[: max(1, int(args.max_samples))]
    labels = ann.set_index("sample_id").reindex(all_samples)["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    raw_values = matrix.to_numpy(dtype=np.float32)

    base_routes = [f"{view}_top{n}_pearson" for view in args.views for n in gene_counts]
    routes = base_routes + list(args.fusion_routes)
    route_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in routes}
    probability_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in routes}
    fold_rows: list[dict[str, Any]] = []

    for fold_no, sample_id in enumerate(samples, start=1):
        heldout_idx = all_samples.index(sample_id)
        train_indices = np.asarray([j for j in range(len(all_samples)) if j != heldout_idx], dtype=int)
        truth = str(labels[heldout_idx])
        fold_scores: dict[str, np.ndarray] = {}
        for view in args.views:
            values = normalize_fold(raw_values, train_indices, view)
            reference, training = build_group_reference(values, labels, groups, heldout_idx)
            gene_order, _ = select_group_discriminative_genes(values, groups, training, max_gene_count)
            sample = values[:, heldout_idx]
            for n in gene_counts:
                route = f"{view}_top{n}_pearson"
                genes = gene_order[: min(n, len(gene_order))]
                scores = correlation_scores(reference, sample, genes)
                fold_scores[route] = scores
                row, prob = detail_from_scores(route, sample_id, truth, groups, scores)
                route_rows[route].append(row)
                probability_rows[route].append(prob)

        fusion_defs: dict[str, list[str]] = {
            "raw_top100_raw_top200_zmean": ["raw_top100_pearson", "raw_top200_pearson"],
            "raw_top50_raw_top100_raw_top200_zmean": ["raw_top50_pearson", "raw_top100_pearson", "raw_top200_pearson"],
            "raw_top200_rank_top200_zmean": ["raw_top200_pearson", "sample_rank_center_top200_pearson"],
            "raw_top200_zscore_top200_rank_top200_zmean": [
                "raw_top200_pearson",
                "sample_zscore_top200_pearson",
                "sample_rank_center_top200_pearson",
            ],
            "raw_top200_zscore_top200_rank_top200_rrf": [
                "raw_top200_pearson",
                "sample_zscore_top200_pearson",
                "sample_rank_center_top200_pearson",
            ],
            "raw_top200_zscore_top200_rank_top200_borda": [
                "raw_top200_pearson",
                "sample_zscore_top200_pearson",
                "sample_rank_center_top200_pearson",
            ],
        }
        for route in args.fusion_routes:
            components = [fold_scores[name] for name in fusion_defs[route] if name in fold_scores]
            if not components:
                continue
            if route.endswith("_rrf"):
                scores = rank_fusion_score(components, method="rrf")
            elif route.endswith("_borda"):
                scores = rank_fusion_score(components, method="borda")
            else:
                scores = np.mean([zscore_scores(component) for component in components], axis=0)
            row, prob = detail_from_scores(route, sample_id, truth, groups, scores)
            route_rows[route].append(row)
            probability_rows[route].append(prob)

        fold_rows.append({"fold": fold_no, "sample_id": sample_id, "truth": truth})

    details = {route: pd.DataFrame(rows) for route, rows in route_rows.items() if rows}
    probabilities = {route: pd.DataFrame(rows) for route, rows in probability_rows.items() if rows}
    metrics = {route: summarize(details[route], probabilities[route], groups) for route in details}
    base_route = BASELINE_ROUTE if BASELINE_ROUTE in details else next(iter(details))
    changes = {route: paired_change(details[base_route], details[route]) for route in details if route != base_route}
    pvalues = {
        route: {
            "top1": paired_pvalue(change["top1_gains"], change["top1_losses"]),
            "top3": paired_pvalue(change["top3_gains"], change["top3_losses"]),
        }
        for route, change in changes.items()
    }
    best_top1 = max(metrics, key=lambda route: (metrics[route]["top1_accuracy"], metrics[route]["top3_accuracy"]))
    best_top3 = max(metrics, key=lambda route: (metrics[route]["top3_accuracy"], metrics[route]["top1_accuracy"]))
    base = metrics[base_route]
    best = metrics[best_top1]
    improved = best["top1_accuracy"] > base["top1_accuracy"] and best["top3_accuracy"] >= base["top3_accuracy"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "endpoint": args.endpoint,
        "validation_design": "strict outer LOSO multiview Network fusion screen; all references and genes are fold-local",
        "n_test_samples": int(len(samples)),
        "n_classes": int(len(groups)),
        "views": args.views,
        "gene_counts": gene_counts,
        "baseline_route": base_route,
        "routes": metrics,
        "paired_changes_vs_baseline": changes,
        "paired_pvalues_vs_baseline": pvalues,
        "best_top1_route": best_top1,
        "best_top3_route": best_top3,
        "meets_internal_adoption_rule": bool(improved),
        "decision": (
            f"{best_top1} improves Top1 without lowering Top3 versus {base_route}; keep as candidate."
            if improved
            else f"No multiview route improves Top1 under the no-Top3-loss constraint versus {base_route}."
        ),
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, frame in details.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_multiview_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
