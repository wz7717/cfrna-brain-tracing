#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_bo2023_hierarchical_region_correlation_validation import rank_candidates, score_ranked_route  # noqa: E402
from scripts.run_bo2023_loso_validation import build_region_reference, correlation_scores, read_annotations, read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import build_group_reference, select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import DEFAULT_NETWORK_DETAIL  # noqa: E402
from scripts.run_bo2023_v2_loso_validation import DEFAULT_GENE_MAP, DEFAULT_MATRIX, DEFAULT_SAMPLE_INFO, map_matrix_to_symbols  # noqa: E402


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_exact_region_top5_conservative_rerank_loso_814_20260529"
BASELINE_ROUTE = "official_top3_beam_top50_top100_zfusion_w0p25"
RERANK_ROUTE = "top5_beam_conservative_rank4_5_into_top3"
TOP50_WEIGHT = 0.25


def paired_changes(base: pd.DataFrame, tested: pd.DataFrame) -> dict[str, int]:
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


def summarize(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "median_true_rank": float(detail["true_rank"].median()),
        "n_reranked": int(detail.get("reranked", pd.Series(0, index=detail.index)).sum()),
    }


def zscore_subset(scores: np.ndarray, indices: np.ndarray) -> np.ndarray:
    out = np.zeros_like(scores, dtype=np.float64)
    vals = scores[indices].astype(float, copy=False)
    std = float(vals.std())
    out[indices] = 0.0 if std <= 1e-12 else (vals - float(vals.mean())) / std
    return out


def fused_region_scores(
    reference: np.ndarray,
    sample: np.ndarray,
    candidate_indices: np.ndarray,
    gene_order: np.ndarray,
) -> np.ndarray:
    rows50 = gene_order[: min(50, len(gene_order))]
    rows100 = gene_order[: min(100, len(gene_order))]
    scores50 = correlation_scores(reference, sample, rows50)
    scores100 = correlation_scores(reference, sample, rows100)
    return TOP50_WEIGHT * zscore_subset(scores50, candidate_indices) + (1.0 - TOP50_WEIGHT) * zscore_subset(scores100, candidate_indices)


def select_fold_gene_order(
    values: np.ndarray,
    candidate_regions: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    top_n: int = 100,
) -> np.ndarray:
    training = {}
    for region in candidate_regions:
        indices = region_indices[region]
        indices = indices[indices != heldout_idx]
        if len(indices):
            training[region] = indices
    if len(training) < 2:
        return np.arange(values.shape[0], dtype=int)
    rows, _ = select_group_discriminative_genes(values, sorted(training), training, top_n)
    return rows.astype(int)


def network_rank_top5(
    values: np.ndarray,
    labels: np.ndarray,
    networks: list[str],
    heldout_idx: int,
    sample: np.ndarray,
    corrected_top3: list[str],
    top_n_genes: int,
) -> list[str]:
    reference, training = build_group_reference(values, labels, networks, heldout_idx)
    rows, _ = select_group_discriminative_genes(values, networks, training, top_n_genes)
    scores = correlation_scores(reference, sample, rows)
    natural = [networks[int(j)] for j in np.argsort(scores)[::-1]]
    out = list(dict.fromkeys(corrected_top3 + [network for network in natural if network not in corrected_top3]))
    return out[:5]


def training_confusion_pairs(
    values: np.ndarray,
    reference: np.ndarray,
    regions: list[str],
    candidates: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    gene_order: np.ndarray,
    max_pairs_per_truth: int,
    min_pair_errors: int,
) -> set[tuple[str, str]]:
    candidate_indices = np.asarray([regions.index(region) for region in candidates], dtype=int)
    counts: Counter[tuple[str, str]] = Counter()
    for truth in candidates:
        truth_indices = region_indices[truth]
        truth_indices = truth_indices[truth_indices != heldout_idx]
        if len(truth_indices) <= 1:
            continue
        truth_pos = regions.index(truth)
        for sample_idx in truth_indices:
            local_reference = reference.copy()
            remaining = truth_indices[truth_indices != sample_idx]
            local_reference[:, truth_pos] = values[:, remaining].mean(axis=1, dtype=np.float64)
            scores = fused_region_scores(local_reference, values[:, sample_idx], candidate_indices, gene_order)
            ranked = rank_candidates(scores, regions, candidate_indices)
            if ranked and ranked[0] != truth:
                counts[(truth, ranked[0])] += 1
    selected: set[tuple[str, str]] = set()
    for truth in sorted(candidates):
        retained = 0
        rows = sorted(
            [(pred, count) for (label, pred), count in counts.items() if label == truth],
            key=lambda item: (-item[1], item[0]),
        )
        for pred, count in rows:
            if count >= min_pair_errors and retained < max_pairs_per_truth:
                selected.add(tuple(sorted((truth, pred))))
                retained += 1
    return selected


def build_pair_models(
    values: np.ndarray,
    pairs: set[tuple[str, str]],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    gene_pool: np.ndarray,
    top_n_genes: int,
) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    models: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    for left, right in sorted(pairs):
        left_idx = region_indices[left][region_indices[left] != heldout_idx]
        right_idx = region_indices[right][region_indices[right] != heldout_idx]
        if len(left_idx) < 2 or len(right_idx) < 2:
            continue
        left_values = values[gene_pool[:, None], left_idx]
        right_values = values[gene_pool[:, None], right_idx]
        delta = left_values.mean(axis=1) - right_values.mean(axis=1)
        pooled = (left_values.var(axis=1, ddof=1) + right_values.var(axis=1, ddof=1)) / 2.0
        effect = np.abs(delta) / np.sqrt(pooled + 1e-8)
        local = np.argsort(effect)[::-1][: min(top_n_genes, len(gene_pool))]
        genes = gene_pool[local].astype(int)
        pair_ref = np.column_stack(
            [
                values[genes[:, None], left_idx].mean(axis=1, dtype=np.float64),
                values[genes[:, None], right_idx].mean(axis=1, dtype=np.float64),
            ]
        )
        models[(left, right)] = (genes, pair_ref)
    return models


def pair_margin_for_challenger(
    sample: np.ndarray,
    challenger: str,
    incumbent: str,
    pair_models: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]],
) -> float | None:
    key = tuple(sorted((challenger, incumbent)))
    model = pair_models.get(key)
    if model is None:
        return None
    genes, pair_ref = model
    pair_scores = correlation_scores(pair_ref, sample[genes])
    if key[0] == challenger:
        return float(pair_scores[0] - pair_scores[1])
    return float(pair_scores[1] - pair_scores[0])


def conservative_rerank(
    baseline_ranked: list[str],
    top5_ranked: list[str],
    sample: np.ndarray,
    pair_models: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]],
    margin_threshold: float,
) -> tuple[list[str], bool, str, float]:
    ranked = baseline_ranked.copy()
    if len(ranked) < 3:
        return ranked, False, "", 0.0
    best: tuple[int, str, str, float] | None = None
    for challenger in top5_ranked[3:5]:
        if challenger in ranked[:3]:
            continue
        incumbent = ranked[2]
        margin = pair_margin_for_challenger(sample, challenger, incumbent, pair_models)
        if margin is None or margin < margin_threshold:
            continue
        candidate = (3, challenger, incumbent, margin)
        if best is None or margin > best[3]:
            best = candidate
    if best is None:
        return ranked, False, "", 0.0
    _, challenger, incumbent, margin = best
    ranked[2] = challenger
    tail = [region for region in top5_ranked if region not in ranked[:3]]
    ranked = ranked[:3] + tail + [region for region in baseline_ranked[3:] if region not in ranked[:3] and region not in tail]
    return ranked, True, f"{challenger} > {incumbent}", float(margin)


def score_route(route: str, sample_id: str, truth_region: str, truth_network: str, network_top: list[str], ranked: list[str], n_total_regions: int, reranked: bool = False, pair: str = "", margin: float = 0.0) -> dict[str, Any]:
    detail = score_ranked_route(
        route,
        sample_id,
        truth_region,
        truth_network,
        network_top[0],
        ranked,
        n_total_regions,
        int(network_top[0] == truth_network),
        int(truth_network in network_top[:3]),
    )
    detail["network_top5_hit"] = int(truth_network in network_top[:5])
    detail["reranked"] = int(reranked)
    detail["rerank_pair"] = pair
    detail["rerank_margin"] = float(margin)
    detail["network_beam"] = " | ".join(network_top)
    return detail


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [BASELINE_ROUTE, RERANK_ROUTE]
    labels = ["Current official\nTop3 beam fusion", "Top5 conservative\nrank4/5->Top3"]
    top1 = [metrics[route]["top1_accuracy"] for route in routes]
    top3 = [metrics[route]["top3_accuracy"] for route in routes]
    x = np.arange(len(routes))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.6, 5.4), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Exact Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Exact Top3", color="#009E73")
    ax.set_ylabel("Accuracy")
    ax.set_ylim(0, min(1.0, max(top3) + 0.14))
    ax.set_xticks(x, labels)
    ax.set_title("Bo2023 exact Region: conservative Top5 rerank, strict LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "exact_region_top5_conservative_rerank.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "exact_region_top5_conservative_rerank.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict LOSO conservative Top5-beam rerank for exact Bo2023 Region Top3.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--network-top-n-genes", type=int, default=200)
    parser.add_argument("--pair-gene-pool", type=int, default=300)
    parser.add_argument("--pair-top-n-genes", type=int, default=50)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=2)
    parser.add_argument("--margin-grid", type=float, nargs="+", default=[0.02, 0.05, 0.08, 0.12])
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
    labels = ann.set_index("sample_id").reindex(sample_ids)["endpoint_label"].to_numpy(dtype=str)
    networks = sorted(set(labels))
    network_detail = pd.read_csv(args.network_detail).set_index("sample_id")

    region_counts = ann.groupby("region_id")["sample_id"].size()
    singletons = set(region_counts[region_counts < 2].index)
    selected = ann[~ann["region_id"].isin(singletons)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)
    if args.max_samples is not None:
        selected = selected.head(max(1, int(args.max_samples))).copy()

    route_rows = {BASELINE_ROUTE: [], RERANK_ROUTE: []}
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        corrected_top3 = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        network_top5 = network_rank_top5(values, labels, networks, heldout_idx, sample, corrected_top3, args.network_top_n_genes)

        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)

        top3_candidates = sorted(training_ann.loc[training_ann["endpoint_label"].isin(corrected_top3), "region_id"].unique().tolist())
        top5_candidates = sorted(training_ann.loc[training_ann["endpoint_label"].isin(network_top5), "region_id"].unique().tolist())
        top3_indices = np.asarray([region_pos[region] for region in top3_candidates], dtype=int)
        top5_indices = np.asarray([region_pos[region] for region in top5_candidates], dtype=int)
        top3_genes = select_fold_gene_order(values, top3_candidates, region_indices, heldout_idx, 100)
        top5_genes = select_fold_gene_order(values, top5_candidates, region_indices, heldout_idx, 100)

        base_scores = fused_region_scores(reference, sample, top3_indices, top3_genes)
        baseline_ranked = rank_candidates(base_scores, regions, top3_indices)
        route_rows[BASELINE_ROUTE].append(
            score_route(BASELINE_ROUTE, sample_id, truth_region, truth_network, corrected_top3, baseline_ranked, len(regions))
        )

        top5_scores = fused_region_scores(reference, sample, top5_indices, top5_genes)
        top5_ranked = rank_candidates(top5_scores, regions, top5_indices)
        pairs = training_confusion_pairs(
            values, reference, regions, top5_candidates, region_indices, heldout_idx, top5_genes,
            args.max_pairs_per_truth, args.min_pair_errors,
        )
        pair_pool = select_fold_gene_order(values, top5_candidates, region_indices, heldout_idx, args.pair_gene_pool)
        pair_models = build_pair_models(values, pairs, region_indices, heldout_idx, pair_pool, args.pair_top_n_genes)

        # Fold-internal conservative threshold selection: choose the largest threshold that still has at least one
        # selected pair model; this minimizes trigger rate and keeps Top1 untouched by construction.
        selected_margin = max(float(x) for x in args.margin_grid)
        if not pair_models:
            selected_margin = float("inf")
        reranked, switched, pair, margin = conservative_rerank(
            baseline_ranked, top5_ranked, sample, pair_models, selected_margin
        )
        route_rows[RERANK_ROUTE].append(
            score_route(
                RERANK_ROUTE, sample_id, truth_region, truth_network, network_top5, reranked, len(regions),
                reranked=switched, pair=pair, margin=margin,
            )
        )
        fold_rows.append(
            {
                "fold": fold,
                "sample_id": sample_id,
                "truth_region": truth_region,
                "truth_network": truth_network,
                "network_top3_hit": int(truth_network in corrected_top3),
                "network_top5_hit": int(truth_network in network_top5),
                "n_top3_candidate_regions": int(len(top3_candidates)),
                "n_top5_candidate_regions": int(len(top5_candidates)),
                "n_pair_models": int(len(pair_models)),
                "selected_margin": float(selected_margin) if np.isfinite(selected_margin) else "inf",
                "reranked": int(switched),
            }
        )

    details = {route: pd.DataFrame(rows) for route, rows in route_rows.items()}
    metrics = {route: summarize(frame) for route, frame in details.items()}
    changes = {RERANK_ROUTE: paired_changes(details[BASELINE_ROUTE], details[RERANK_ROUTE])}
    pvalues = {
        RERANK_ROUTE: {
            "top1": paired_pvalue(changes[RERANK_ROUTE]["top1_gains"], changes[RERANK_ROUTE]["top1_losses"]),
            "top3": paired_pvalue(changes[RERANK_ROUTE]["top3_gains"], changes[RERANK_ROUTE]["top3_losses"]),
        }
    }
    base = metrics[BASELINE_ROUTE]
    tested = metrics[RERANK_ROUTE]
    meets_rule = (
        tested["top1_hits"] >= base["top1_hits"]
        and tested["top3_hits"] > base["top3_hits"]
        and changes[RERANK_ROUTE]["top1_losses"] == 0
    )
    decision = (
        "Candidate passes the conservative adoption screen: exact Top3 improved, exact Top1 was not reduced, and Top1 was never changed."
        if meets_rule
        else "Candidate does not pass the adoption screen; keep the current official Top3-beam fusion route."
    )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; current official Top3 beam baseline; Top5 beam only for conservative rank4/5-to-Top3 rescue",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singletons).sum()),
        "top50_weight": TOP50_WEIGHT,
        "margin_grid": [float(x) for x in args.margin_grid],
        "routes": metrics,
        "paired_changes": changes,
        "paired_pvalues": pvalues,
        "meets_adoption_rule": bool(meets_rule),
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_rerank_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
