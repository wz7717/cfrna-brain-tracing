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

from scripts.run_bo2023_exact_region_gene_count_variants import (  # noqa: E402
    paired_changes,
    paired_pvalue,
    select_fold_gene_order,
    zscore_candidate_scores,
)
from scripts.run_bo2023_hierarchical_region_correlation_validation import rank_candidates  # noqa: E402
from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    correlation_scores,
    read_annotations,
    read_vsd_matrix,
)
from scripts.run_bo2023_resolution_tier_validation import (  # noqa: E402
    DEFAULT_NETWORK_DETAIL,
    build_resolution_groups,
    candidate_training_indices,
    distinct_ranked_groups,
    region_network_assignment,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_group_aggregate_scoring_loso_814_20260604"
BASELINE_ROUTE = "first_region_group_from_zfusion"
MAX_ROUTE = "group_max_region_score"
MAX_MEAN_ROUTE = "group_max_plus_0p10_mean_score"


def group_order_first(ranked_regions: list[str], annotations: dict[str, dict[str, Any]]) -> list[str]:
    return distinct_ranked_groups(ranked_regions, annotations)


def group_order_aggregate(
    candidates: list[str],
    scores: np.ndarray,
    regions: list[str],
    annotations: dict[str, dict[str, Any]],
    mean_weight: float,
) -> list[str]:
    region_pos = {region: j for j, region in enumerate(regions)}
    grouped: dict[str, list[float]] = {}
    for region in candidates:
        group = str(annotations[region]["resolution_group"])
        grouped.setdefault(group, []).append(float(scores[region_pos[region]]))
    rows = []
    for group, values in grouped.items():
        score = max(values) + float(mean_weight) * (sum(values) / len(values))
        rows.append((group, score))
    return [group for group, _ in sorted(rows, key=lambda item: item[1], reverse=True)]


def make_detail(route: str, sample_id: str, truth_group: str, ranked_groups: list[str]) -> dict[str, Any]:
    true_rank = ranked_groups.index(truth_group) + 1 if truth_group in ranked_groups else len(ranked_groups) + 1
    padded = ranked_groups[:3] + [""] * max(0, 3 - len(ranked_groups))
    return {
        "route": route,
        "sample_id": sample_id,
        "true_resolution_group": truth_group,
        "pred_group_top1": padded[0],
        "pred_group_top2": padded[1],
        "pred_group_top3": padded[2],
        "group_true_rank": int(true_rank),
        "group_hit1": int(true_rank == 1),
        "group_hit3": int(true_rank <= 3),
    }


def summarize(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "group_top1_hits": int(detail["group_hit1"].sum()),
        "group_top1_accuracy": float(detail["group_hit1"].mean()),
        "group_top3_hits": int(detail["group_hit3"].sum()),
        "group_top3_accuracy": float(detail["group_hit3"].mean()),
        "median_group_true_rank": float(detail["group_true_rank"].median()),
    }


def group_paired_changes(base: pd.DataFrame, tested: pd.DataFrame) -> dict[str, int]:
    return {
        "top1_gains": int(((base["group_hit1"] == 0) & (tested["group_hit1"] == 1)).sum()),
        "top1_losses": int(((base["group_hit1"] == 1) & (tested["group_hit1"] == 0)).sum()),
        "top3_gains": int(((base["group_hit3"] == 0) & (tested["group_hit3"] == 1)).sum()),
        "top3_losses": int(((base["group_hit3"] == 1) & (tested["group_hit3"] == 0)).sum()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict LOSO comparison of first-occurrence versus aggregate scoring for resolution groups."
    )
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--top50-weight", type=float, default=0.25)
    parser.add_argument("--group-mean-weight", type=float, default=0.10)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.15)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=8)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    matrix = map_matrix_to_symbols(read_vsd_matrix(args.matrix), args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    network_ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet, usecols=["No.", args.network_col])
    network_ann["sample_id"] = network_ann["No."].astype(str).str.strip()
    network_ann["endpoint_label"] = network_ann[args.network_col].fillna("NA").astype(str).str.strip()
    ann = ann.merge(network_ann[["sample_id", "endpoint_label"]], on="sample_id", how="left")
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()

    values = matrix.to_numpy(dtype=np.float32)
    sample_ids = matrix.columns.astype(str).tolist()
    sample_pos = {sample_id: j for j, sample_id in enumerate(sample_ids)}
    reference_all, regions, _, region_indices = build_region_reference(values, sample_ids, ann)
    region_pos = {region: j for j, region in enumerate(regions)}
    network_detail = pd.read_csv(args.network_detail).set_index("sample_id")

    region_counts = ann.groupby("region_id")["sample_id"].size()
    singletons = set(region_counts[region_counts < 2].index)
    selected = ann[~ann["region_id"].isin(singletons)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)
    if args.max_samples is not None:
        selected = selected.head(max(1, int(args.max_samples))).copy()

    rows: dict[str, list[dict[str, Any]]] = {BASELINE_ROUTE: [], MAX_ROUTE: [], MAX_MEAN_ROUTE: []}
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        network_top = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        candidates = sorted(training_ann.loc[training_ann["endpoint_label"].isin(network_top), "region_id"].unique())
        candidate_indices = np.asarray([region_pos[region] for region in candidates], dtype=int)
        training = candidate_training_indices(candidates, region_indices, heldout_idx)
        local_rows = select_fold_gene_order(values, candidates, region_indices, heldout_idx, args.local_top_n_genes)
        annotations, _ = build_resolution_groups(
            values,
            candidates,
            training,
            region_network_assignment(training_ann, candidates),
            local_rows,
            args.min_resolution_samples,
            args.min_merge_samples,
            args.min_pair_errors,
            args.min_confusion_rate,
            args.similarity_threshold,
            args.merge_similarity_threshold,
            args.max_group_size,
        )
        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)

        gene_order = select_fold_gene_order(values, candidates, region_indices, heldout_idx, 100)
        scores50 = correlation_scores(reference, sample, gene_order[:50])
        scores100 = correlation_scores(reference, sample, gene_order[:100])
        fused = (
            float(args.top50_weight) * zscore_candidate_scores(scores50, candidate_indices)
            + (1.0 - float(args.top50_weight)) * zscore_candidate_scores(scores100, candidate_indices)
        )
        ranked_regions = rank_candidates(fused, regions, candidate_indices)
        truth_group = str(annotations[truth_region]["resolution_group"]) if truth_region in annotations else "outside_network_beam"
        orders = {
            BASELINE_ROUTE: group_order_first(ranked_regions, annotations),
            MAX_ROUTE: group_order_aggregate(candidates, fused, regions, annotations, mean_weight=0.0),
            MAX_MEAN_ROUTE: group_order_aggregate(candidates, fused, regions, annotations, mean_weight=args.group_mean_weight),
        }
        for route, ranked_groups in orders.items():
            rows[route].append(make_detail(route, sample_id, truth_group, ranked_groups))
        fold_rows.append(
            {
                "fold": fold,
                "sample_id": sample_id,
                "truth_region": truth_region,
                "truth_group": truth_group,
                "network_top3_hit": int(str(row.endpoint_label) in network_top),
                "n_candidate_regions": int(len(candidates)),
                "n_candidate_groups": int(len(set(annotations[r]["resolution_group"] for r in candidates))),
            }
        )

    details = {route: pd.DataFrame(data) for route, data in rows.items()}
    metrics = {route: summarize(detail) for route, detail in details.items()}
    changes = {route: group_paired_changes(details[BASELINE_ROUTE], details[route]) for route in [MAX_ROUTE, MAX_MEAN_ROUTE]}
    pvalues = {
        route: {
            "top1": paired_pvalue(change["top1_gains"], change["top1_losses"]),
            "top3": paired_pvalue(change["top3_gains"], change["top3_losses"]),
        }
        for route, change in changes.items()
    }
    best_route = max([MAX_ROUTE, MAX_MEAN_ROUTE], key=lambda route: (
        metrics[route]["group_top1_accuracy"],
        metrics[route]["group_top3_accuracy"],
    ))
    base = metrics[BASELINE_ROUTE]
    best = metrics[best_route]
    improved = best["group_top1_accuracy"] > base["group_top1_accuracy"] and best["group_top3_accuracy"] >= base["group_top3_accuracy"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; pairwise-corrected Top3 Network beam; Top50/Top100 z-fusion; fold-local adaptive max8 resolution groups; group-ranking aggregation only",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singletons).sum()),
        "top50_weight": float(args.top50_weight),
        "group_mean_weight": float(args.group_mean_weight),
        "local_top_n_genes": int(args.local_top_n_genes),
        "min_confusion_rate": float(args.min_confusion_rate),
        "max_group_size": int(args.max_group_size),
        "routes": metrics,
        "paired_changes_vs_first_occurrence": changes,
        "paired_pvalues_vs_first_occurrence": pvalues,
        "best_route": best_route,
        "meets_internal_adoption_rule": bool(improved),
        "decision": (
            f"{best_route} improves Group Top1 without lowering Group Top3; keep as a deployable secondary group-ranking candidate."
            if improved
            else "Aggregate group scoring does not improve Group Top1 under the no-Top3-loss constraint; keep first-occurrence grouping."
        ),
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_group_aggregate_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
