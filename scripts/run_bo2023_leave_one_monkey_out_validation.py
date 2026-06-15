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
from scripts.run_bo2023_hierarchical_region_correlation_validation import rank_candidates  # noqa: E402
from scripts.run_bo2023_loso_validation import correlation_scores, read_annotations, read_vsd_matrix, softmax  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_network_pairwise_correlation_validation import (  # noqa: E402
    BASELINE_ROUTE as NETWORK_BASELINE_ROUTE,
    PAIR_TOP3_ROUTE as NETWORK_PAIR_TOP3_ROUTE,
    build_pair_models,
    derive_training_confusion_pairs,
    evaluate_pairwise_rescue,
    make_detail as make_network_detail,
)
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups, score_route  # noqa: E402
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_20260601"
EXACT_ROUTE = "top3_beam_local_top50_top100_zfusion_w0p25"
GROUP_ROUTE = "top3_network_beam_local_region_candidates"


def normalize_fold(raw_values: np.ndarray, train_indices: np.ndarray, variant: str) -> np.ndarray:
    values = raw_values.astype(np.float32, copy=True)
    if variant == "raw":
        return values
    if variant == "sample_zscore":
        mean = values.mean(axis=0, keepdims=True, dtype=np.float64)
        std = values.std(axis=0, keepdims=True, dtype=np.float64)
        return ((values - mean) / (std + 1e-6)).astype(np.float32)
    raise ValueError(f"Unsupported full-route normalization variant: {variant}")


def build_label_reference(
    values: np.ndarray,
    labels: np.ndarray,
    groups: list[str],
    train_indices: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    training: dict[str, np.ndarray] = {}
    reference = np.zeros((values.shape[0], len(groups)), dtype=np.float32)
    for j, group in enumerate(groups):
        idx = train_indices[labels[train_indices] == group]
        if idx.size == 0:
            raise ValueError(f"No training samples for label {group!r}")
        training[group] = idx
        reference[:, j] = values[:, idx].mean(axis=1, dtype=np.float64).astype(np.float32)
    return reference, training


def build_region_training(
    ann: pd.DataFrame,
    sample_pos: dict[str, int],
    train_sample_ids: set[str],
) -> dict[str, np.ndarray]:
    training: dict[str, np.ndarray] = {}
    for region, rows in ann[ann["sample_id"].isin(train_sample_ids)].groupby("region_id"):
        idx = np.asarray([sample_pos[str(x)] for x in rows["sample_id"].astype(str)], dtype=int)
        if idx.size:
            training[str(region)] = idx
    return training


def build_region_reference(
    values: np.ndarray,
    regions: list[str],
    training: dict[str, np.ndarray],
) -> np.ndarray:
    return np.column_stack(
        [values[:, training[region]].mean(axis=1, dtype=np.float64) for region in regions]
    ).astype(np.float32)


def zscore(scores: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=np.float64)
    std = float(scores.std())
    if std <= 1e-12:
        return np.zeros_like(scores)
    return (scores - float(scores.mean())) / std


def make_exact_detail(
    sample_id: str,
    monkey_id: str,
    truth_region: str,
    truth_network: str,
    network_top: list[str],
    ranked_regions: list[str],
    n_total_regions: int,
) -> dict[str, Any]:
    true_rank = ranked_regions.index(truth_region) + 1 if truth_region in ranked_regions else n_total_regions + 1
    padded = ranked_regions[:3] + [""] * max(0, 3 - len(ranked_regions))
    return {
        "route": EXACT_ROUTE,
        "sample_id": sample_id,
        "monkey_id": monkey_id,
        "label": truth_region,
        "true_network": truth_network,
        "network_beam": " | ".join(network_top),
        "network_top1_hit": int(network_top[0] == truth_network),
        "network_top3_hit": int(truth_network in network_top),
        "pred_top1": padded[0],
        "pred_top2": padded[1],
        "pred_top3": padded[2],
        "true_rank": int(true_rank),
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "n_candidate_regions": int(len(ranked_regions)),
    }


def summarize_network(detail: pd.DataFrame, probabilities: pd.DataFrame, groups: list[str]) -> dict[str, Any]:
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
        "n_switches": int(detail.get("switched", pd.Series(dtype=int)).sum()),
    }


def summarize_exact(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_by_monkey_top1": float(detail.groupby("monkey_id")["hit1"].mean().mean()),
        "macro_by_monkey_top3": float(detail.groupby("monkey_id")["hit3"].mean().mean()),
        "median_true_rank": float(detail["true_rank"].median()),
        "conditional_top1_given_network_top1": float(detail.loc[detail["network_top1_hit"] == 1, "hit1"].mean()),
        "conditional_top3_given_network_top1": float(detail.loc[detail["network_top1_hit"] == 1, "hit3"].mean()),
    }


def summarize_group(detail: pd.DataFrame) -> dict[str, Any]:
    low = detail[detail["pred_top1_resolution_tier"] == "low_resolution"]
    high = detail[detail["pred_top1_resolution_tier"] == "high_resolution"]
    return {
        "n": int(len(detail)),
        "exact_top1_hits": int(detail["hit1"].sum()),
        "exact_top1_accuracy": float(detail["hit1"].mean()),
        "exact_top3_hits": int(detail["hit3"].sum()),
        "exact_top3_accuracy": float(detail["hit3"].mean()),
        "group_top1_hits": int(detail["group_hit1"].sum()),
        "group_top1_accuracy": float(detail["group_hit1"].mean()),
        "group_top3_hits": int(detail["group_hit3"].sum()),
        "group_top3_accuracy": float(detail["group_hit3"].mean()),
        "macro_by_monkey_group_top1": float(detail.groupby("monkey_id")["group_hit1"].mean().mean()),
        "macro_by_monkey_group_top3": float(detail.groupby("monkey_id")["group_hit3"].mean().mean()),
        "median_exact_true_rank": float(detail["true_rank"].median()),
        "median_group_true_rank": float(detail["group_true_rank"].median()),
        "low_resolution_predictions": int(len(low)),
        "low_resolution_fraction": float(len(low) / len(detail)),
        "low_resolution_exact_top1": float(low["hit1"].mean()) if len(low) else float("nan"),
        "high_resolution_predictions": int(len(high)),
        "high_resolution_exact_top1": float(high["hit1"].mean()) if len(high) else float("nan"),
    }


def per_monkey_metrics(detail: pd.DataFrame, metric_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for monkey_id, frame in detail.groupby("monkey_id", sort=True):
        row: dict[str, Any] = {"monkey_id": monkey_id, "n": int(len(frame))}
        for col in metric_cols:
            row[f"{col}_mean"] = float(frame[col].mean())
            row[f"{col}_hits"] = int(frame[col].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Bo2023 formal-route leave-one-monkey-out validation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--monkey-col", default="MonkeyID")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--global-top-n-genes", type=int, default=200)
    parser.add_argument("--gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=100)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=3)
    parser.add_argument("--pair-min-margin", type=float, default=0.002)
    parser.add_argument("--normalization", choices=["raw", "sample_zscore"], default="raw")
    parser.add_argument("--exact-fusion-weight", type=float, default=0.25)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--group-min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.15)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=8)
    parser.add_argument("--max-monkeys", type=int, default=None, help="Optional deterministic prefix for smoke tests.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = map_matrix_to_symbols(read_vsd_matrix(args.matrix), args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    full_ann = pd.read_excel(
        args.sample_info,
        sheet_name=args.sample_sheet,
        usecols=["No.", args.monkey_col, args.network_col],
    )
    full_ann["sample_id"] = full_ann["No."].astype(str).str.strip()
    full_ann["monkey_id"] = full_ann[args.monkey_col].astype(str).str.strip()
    full_ann["endpoint_label"] = full_ann[args.network_col].fillna("NA").astype(str).str.strip()
    ann = ann.merge(full_ann[["sample_id", "monkey_id", "endpoint_label"]], on="sample_id", how="left")
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    ann = ann.dropna(subset=["monkey_id", "endpoint_label"])

    raw_values = matrix.to_numpy(dtype=np.float32)
    samples = matrix.columns.astype(str).tolist()
    sample_pos = {sample_id: idx for idx, sample_id in enumerate(samples)}
    sample_ann = ann.set_index("sample_id").reindex(samples)
    labels = sample_ann["endpoint_label"].to_numpy(dtype=str)
    monkeys = sample_ann["monkey_id"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    all_regions = sorted(ann["region_id"].unique().tolist())
    monkey_ids = sorted(ann["monkey_id"].unique().tolist())
    if args.max_monkeys is not None:
        monkey_ids = monkey_ids[: max(1, int(args.max_monkeys))]

    network_rows: dict[str, list[dict[str, Any]]] = {
        NETWORK_BASELINE_ROUTE: [],
        NETWORK_PAIR_TOP3_ROUTE: [],
    }
    network_prob_rows: dict[str, list[dict[str, Any]]] = {
        NETWORK_BASELINE_ROUTE: [],
        NETWORK_PAIR_TOP3_ROUTE: [],
    }
    exact_rows: list[dict[str, Any]] = []
    group_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []
    pair_audit_rows: list[pd.DataFrame] = []
    group_audit_rows: list[pd.DataFrame] = []

    for fold_no, monkey_id in enumerate(monkey_ids, start=1):
        test_indices = np.where(monkeys == monkey_id)[0]
        train_indices = np.where(monkeys != monkey_id)[0]
        values = normalize_fold(raw_values, train_indices, args.normalization)
        train_sample_ids = {samples[i] for i in train_indices}
        test_sample_ids = [samples[i] for i in test_indices]

        network_reference, network_training = build_label_reference(values, labels, groups, train_indices)
        gene_pool, _ = select_group_discriminative_genes(values, groups, network_training, args.gene_pool_size)
        global_genes = gene_pool[: args.global_top_n_genes]
        pairs, pair_errors = derive_training_confusion_pairs(
            values,
            labels,
            groups,
            network_training,
            network_reference,
            global_genes,
            args.max_pairs_per_truth,
            args.min_pair_errors,
        )
        pair_errors.insert(0, "fold", fold_no)
        pair_errors.insert(1, "heldout_monkey_id", monkey_id)
        pair_audit_rows.append(pair_errors)
        pair_models, pair_audit = build_pair_models(
            values,
            network_training,
            network_reference,
            pairs,
            groups,
            gene_pool,
            args.pair_top_n_genes,
        )
        pair_audit.insert(0, "fold", fold_no)
        pair_audit.insert(1, "heldout_monkey_id", monkey_id)
        pair_audit_rows.append(pair_audit)

        region_training = build_region_training(ann, sample_pos, train_sample_ids)
        train_ann = ann[ann["sample_id"].isin(train_sample_ids)].copy()
        fold_info: dict[str, Any] = {
            "fold": fold_no,
            "heldout_monkey_id": monkey_id,
            "n_test_samples": int(len(test_indices)),
            "n_train_samples": int(len(train_indices)),
            "n_train_regions": int(len(region_training)),
            "n_pairwise_network_pairs": int(len(pair_models)),
        }
        fold_exact_evaluable = 0
        fold_group_evaluable = 0

        for sample_id, heldout_idx in zip(test_sample_ids, test_indices):
            truth_network = str(labels[heldout_idx])
            truth_region = str(sample_ann.loc[sample_id, "region_id"])
            scores = correlation_scores(network_reference, values[:, heldout_idx], global_genes)
            order = np.argsort(scores)[::-1].tolist()
            base_detail, base_prob = make_network_detail(
                NETWORK_BASELINE_ROUTE, sample_id, truth_network, groups, order, scores
            )
            pair_detail, pair_prob = evaluate_pairwise_rescue(
                NETWORK_PAIR_TOP3_ROUTE,
                sample_id,
                truth_network,
                values[:, heldout_idx],
                network_reference,
                groups,
                scores,
                pair_models,
                3,
                args.pair_min_margin,
            )
            for detail, prob, route in [
                (base_detail, base_prob, NETWORK_BASELINE_ROUTE),
                (pair_detail, pair_prob, NETWORK_PAIR_TOP3_ROUTE),
            ]:
                detail["monkey_id"] = monkey_id
                prob["monkey_id"] = monkey_id
                network_rows[route].append(detail)
                network_prob_rows[route].append(prob)

            network_top = [pair_detail[f"pred_top{i}"] for i in [1, 2, 3]]
            if truth_region not in region_training:
                continue
            candidates = sorted(
                region
                for region in train_ann.loc[train_ann["endpoint_label"].isin(network_top), "region_id"].unique().tolist()
                if region in region_training
            )
            if not candidates:
                continue

            candidate_training = {region: region_training[region] for region in candidates}
            candidate_reference = build_region_reference(values, candidates, candidate_training)
            if len(candidates) >= 2:
                gene_order, _ = select_group_discriminative_genes(
                    values, candidates, candidate_training, max(100, args.local_top_n_genes)
                )
            else:
                gene_order = np.arange(values.shape[0], dtype=int)
            scores50 = correlation_scores(candidate_reference, values[:, heldout_idx], gene_order[: min(50, len(gene_order))])
            scores100 = correlation_scores(candidate_reference, values[:, heldout_idx], gene_order[: min(100, len(gene_order))])
            fused = args.exact_fusion_weight * zscore(scores50) + (1.0 - args.exact_fusion_weight) * zscore(scores100)
            ranked_exact = [candidates[i] for i in np.argsort(fused)[::-1].tolist()]
            exact_rows.append(
                make_exact_detail(
                    sample_id,
                    monkey_id,
                    truth_region,
                    truth_network,
                    network_top,
                    ranked_exact,
                    len(all_regions),
                )
            )
            fold_exact_evaluable += 1

            local_rows = gene_order[: min(args.local_top_n_genes, len(gene_order))]
            assignment: dict[str, str | None] = {}
            for region in candidates:
                nets = sorted(train_ann.loc[train_ann["region_id"] == region, "endpoint_label"].astype(str).unique().tolist())
                assignment[region] = nets[0] if len(nets) == 1 else None
            annotations, audit = build_resolution_groups(
                values,
                candidates,
                candidate_training,
                assignment,
                local_rows,
                args.min_resolution_samples,
                args.min_merge_samples,
                args.group_min_pair_errors,
                args.min_confusion_rate,
                args.similarity_threshold,
                args.merge_similarity_threshold,
                args.max_group_size,
            )
            scores_group = correlation_scores(candidate_reference, values[:, heldout_idx], local_rows)
            ranked_group = [candidates[i] for i in np.argsort(scores_group)[::-1].tolist()]
            group_detail = score_route(
                GROUP_ROUTE,
                sample_id,
                truth_region,
                truth_network,
                network_top,
                ranked_group,
                annotations,
                len(all_regions),
            )
            group_detail["monkey_id"] = monkey_id
            group_rows.append(group_detail)
            fold_group_evaluable += 1
            if not audit.empty:
                audit.insert(0, "fold", fold_no)
                audit.insert(1, "heldout_monkey_id", monkey_id)
                audit.insert(2, "sample_id", sample_id)
                group_audit_rows.append(audit)

        fold_info["n_exact_evaluable_samples"] = int(fold_exact_evaluable)
        fold_info["n_group_evaluable_samples"] = int(fold_group_evaluable)
        fold_rows.append(fold_info)

    network_details = {route: pd.DataFrame(rows) for route, rows in network_rows.items()}
    network_probs = {route: pd.DataFrame(rows) for route, rows in network_prob_rows.items()}
    exact_detail = pd.DataFrame(exact_rows)
    group_detail = pd.DataFrame(group_rows)

    network_metrics = {
        route: summarize_network(network_details[route], network_probs[route], groups)
        for route in network_details
    }
    exact_metrics = summarize_exact(exact_detail)
    group_metrics = summarize_group(group_detail)
    summary: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "leave-one-monkey-out; all Network, Exact Region, and resolution-group components rebuilt fold-locally",
        "normalization": args.normalization,
        "heldout_monkey_ids": monkey_ids,
        "n_monkeys": int(len(monkey_ids)),
        "n_network_test_samples": int(len(network_rows[NETWORK_PAIR_TOP3_ROUTE])),
        "n_exact_evaluable_samples": int(len(exact_detail)),
        "n_group_evaluable_samples": int(len(group_detail)),
        "n_exact_or_group_excluded_samples": int(len(network_details[NETWORK_PAIR_TOP3_ROUTE]) - len(exact_detail)),
        "network_parameters": {
            "global_top_n_genes": int(args.global_top_n_genes),
            "gene_pool_size": int(args.gene_pool_size),
            "pair_top_n_genes": int(args.pair_top_n_genes),
            "max_pairs_per_truth": int(args.max_pairs_per_truth),
            "min_pair_errors": int(args.min_pair_errors),
            "pair_min_margin": float(args.pair_min_margin),
        },
        "exact_region_parameters": {
            "route": EXACT_ROUTE,
            "top50_top100_zfusion_weight": float(args.exact_fusion_weight),
        },
        "resolution_group_parameters": {
            "local_top_n_genes": int(args.local_top_n_genes),
            "min_resolution_samples": int(args.min_resolution_samples),
            "min_merge_samples": int(args.min_merge_samples),
            "min_pair_errors": int(args.group_min_pair_errors),
            "min_confusion_rate": float(args.min_confusion_rate),
            "similarity_threshold": float(args.similarity_threshold),
            "merge_similarity_threshold": float(args.merge_similarity_threshold),
            "max_group_size": int(args.max_group_size),
        },
        "routes": {
            "network": network_metrics,
            "exact_region": {EXACT_ROUTE: exact_metrics},
            "resolution_group": {GROUP_ROUTE: group_metrics},
        },
    }

    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_monkey_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([{"route": route, **metric} for route, metric in network_metrics.items()]).to_csv(
        args.outdir / "network_route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame([{"route": EXACT_ROUTE, **exact_metrics}]).to_csv(
        args.outdir / "exact_region_route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame([{"route": GROUP_ROUTE, **group_metrics}]).to_csv(
        args.outdir / "resolution_group_route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, frame in network_details.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    for route, frame in network_probs.items():
        frame.to_csv(args.outdir / f"{route}_probabilities.csv", index=False, encoding="utf-8-sig")
    exact_detail.to_csv(args.outdir / f"{EXACT_ROUTE}_detail.csv", index=False, encoding="utf-8-sig")
    group_detail.to_csv(args.outdir / f"{GROUP_ROUTE}_detail.csv", index=False, encoding="utf-8-sig")
    per_monkey_metrics(network_details[NETWORK_PAIR_TOP3_ROUTE], ["hit1", "hit3"]).to_csv(
        args.outdir / "network_pairwise_top3_per_monkey_metrics.csv", index=False, encoding="utf-8-sig"
    )
    per_monkey_metrics(exact_detail, ["hit1", "hit3"]).to_csv(
        args.outdir / "exact_region_per_monkey_metrics.csv", index=False, encoding="utf-8-sig"
    )
    per_monkey_metrics(group_detail, ["group_hit1", "group_hit3", "hit1", "hit3"]).to_csv(
        args.outdir / "resolution_group_per_monkey_metrics.csv", index=False, encoding="utf-8-sig"
    )
    if pair_audit_rows:
        pd.concat(pair_audit_rows, ignore_index=True).to_csv(
            args.outdir / "fold_network_pairwise_audit.csv", index=False, encoding="utf-8-sig"
        )
    if group_audit_rows:
        pd.concat(group_audit_rows, ignore_index=True).to_csv(
            args.outdir / "fold_resolution_group_pair_merge_audit.csv", index=False, encoding="utf-8-sig"
        )
    (args.outdir / "validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
