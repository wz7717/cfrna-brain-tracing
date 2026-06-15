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

from scripts.run_bo2023_hierarchical_region_correlation_validation import (  # noqa: E402
    rank_candidates,
    score_ranked_route,
)
from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    correlation_scores,
    read_annotations,
    read_vsd_matrix,
)
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import DEFAULT_NETWORK_DETAIL  # noqa: E402
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_exact_region_gene_count_variants_loso_814_20260529"
BASELINE_ROUTE = "top3_beam_local_top200_baseline"
FUSION_WEIGHTS = [0.25, 0.50, 0.75]


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
    correct_network = detail[detail["network_top1_hit"] == 1]
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "median_true_rank": float(detail["true_rank"].median()),
        "conditional_top1_given_network_top1": float(correct_network["hit1"].mean()) if len(correct_network) else float("nan"),
        "conditional_top3_given_network_top1": float(correct_network["hit3"].mean()) if len(correct_network) else float("nan"),
    }


def zscore_candidate_scores(scores: np.ndarray, candidate_indices: np.ndarray) -> np.ndarray:
    out = np.zeros_like(scores, dtype=np.float64)
    candidate_scores = scores[candidate_indices].astype(np.float64, copy=False)
    std = float(candidate_scores.std())
    if std <= 1e-12:
        out[candidate_indices] = 0.0
    else:
        out[candidate_indices] = (candidate_scores - float(candidate_scores.mean())) / std
    return out


def select_fold_gene_order(
    values: np.ndarray,
    candidate_regions: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    max_genes: int,
) -> np.ndarray:
    training: dict[str, np.ndarray] = {}
    for region in candidate_regions:
        indices = region_indices[region]
        indices = indices[indices != heldout_idx]
        if len(indices):
            training[region] = indices
    if len(training) < 2:
        return np.arange(values.shape[0], dtype=int)
    rows, _ = select_group_discriminative_genes(values, sorted(training), training, max_genes)
    return rows


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]], routes: list[str]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [route.replace("top3_beam_local_top", "Top").replace("_genes", "") for route in routes]
    top1 = [metrics[route]["top1_accuracy"] for route in routes]
    top3 = [metrics[route]["top3_accuracy"] for route in routes]
    x = np.arange(len(routes))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11.2, 5.6), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Exact Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Exact Top3", color="#009E73")
    ax.set_ylabel("Accuracy")
    ax.set_ylim(0, min(1.0, max(top3) + 0.14))
    ax.set_xticks(x, labels)
    ax.set_title("Bo2023 exact Region: local gene-count variants, strict LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "exact_region_gene_count_variants.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "exact_region_gene_count_variants.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any], routes: list[str]) -> None:
    lines = [
        "# Exact Region gene-count variants strict LOSO",
        "",
        "## Design",
        "",
        "- Strict outer LOSO unchanged.",
        "- Network gate unchanged: pairwise-corrected SaleemNetworks Top3 beam.",
        "- Only the fold-local Region discriminative gene window is varied.",
        "",
        "## Results",
        "",
        "| Route | Exact Top1 | Exact Top3 | Median true rank |",
        "| --- | ---: | ---: | ---: |",
    ]
    for route in routes:
        m = summary["routes"][route]
        lines.append(
            f"| {route} | {m['top1_hits']}/{summary['n_test_samples']} ({m['top1_accuracy']:.1%}) "
            f"| {m['top3_hits']}/{summary['n_test_samples']} ({m['top3_accuracy']:.1%}) "
            f"| {m['median_true_rank']:.1f} |"
        )
    lines.extend(["", "## Decision", "", summary["decision"], ""])
    (outdir / "exact_region_gene_count_variants_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict LOSO exact Region comparison over local gene-count variants.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--gene-counts", type=int, nargs="+", default=[50, 100, 200, 300, 500, 1000])
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_counts = sorted({int(x) for x in args.gene_counts if int(x) > 0})
    max_genes = max(gene_counts)
    routes = [f"top3_beam_local_top{n}_genes" for n in gene_counts]
    fusion_routes = [
        f"top3_beam_local_top50_top100_zfusion_w{str(weight).replace('.', 'p')}"
        for weight in FUSION_WEIGHTS
        if 50 in gene_counts and 100 in gene_counts
    ]
    all_routes = routes + fusion_routes

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

    route_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in all_routes}
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        network_top = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        candidates = sorted(
            training_ann.loc[training_ann["endpoint_label"].isin(network_top), "region_id"].unique().tolist()
        )
        candidate_indices = np.asarray([region_pos[region] for region in candidates], dtype=int)

        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)
        gene_order = select_fold_gene_order(values, candidates, region_indices, heldout_idx, max_genes)

        fold_info = {
            "fold": fold,
            "sample_id": sample_id,
            "truth_region": truth_region,
            "truth_network": truth_network,
            "network_top3_hit": int(truth_network in network_top),
            "n_candidate_regions": int(len(candidates)),
        }
        score_by_n: dict[int, np.ndarray] = {}
        for n, route in zip(gene_counts, routes):
            rows = gene_order[: min(n, len(gene_order))]
            scores = correlation_scores(reference, sample, rows)
            score_by_n[n] = scores
            ranked = rank_candidates(scores, regions, candidate_indices)
            detail = score_ranked_route(
                route,
                sample_id,
                truth_region,
                truth_network,
                network_top[0],
                ranked,
                len(regions),
                int(network_top[0] == truth_network),
                int(truth_network in network_top),
            )
            route_rows[route].append(detail)
            fold_info[f"{route}_true_rank"] = int(detail["true_rank"])
        if 50 in score_by_n and 100 in score_by_n:
            z50 = zscore_candidate_scores(score_by_n[50], candidate_indices)
            z100 = zscore_candidate_scores(score_by_n[100], candidate_indices)
            for weight in FUSION_WEIGHTS:
                route = f"top3_beam_local_top50_top100_zfusion_w{str(weight).replace('.', 'p')}"
                fused = weight * z50 + (1.0 - weight) * z100
                ranked = rank_candidates(fused, regions, candidate_indices)
                detail = score_ranked_route(
                    route,
                    sample_id,
                    truth_region,
                    truth_network,
                    network_top[0],
                    ranked,
                    len(regions),
                    int(network_top[0] == truth_network),
                    int(truth_network in network_top),
                )
                route_rows[route].append(detail)
                fold_info[f"{route}_true_rank"] = int(detail["true_rank"])
        fold_rows.append(fold_info)

    details = {route: pd.DataFrame(rows) for route, rows in route_rows.items()}
    metrics = {route: summarize(frame) for route, frame in details.items()}
    base_route = "top3_beam_local_top200_genes" if "top3_beam_local_top200_genes" in details else routes[0]
    changes = {route: paired_changes(details[base_route], details[route]) for route in all_routes if route != base_route}
    pvalues = {
        route: {
            "top1": paired_pvalue(change["top1_gains"], change["top1_losses"]),
            "top3": paired_pvalue(change["top3_gains"], change["top3_losses"]),
        }
        for route, change in changes.items()
    }
    best_top1 = max(all_routes, key=lambda route: (metrics[route]["top1_accuracy"], metrics[route]["top3_accuracy"]))
    best_top3 = max(all_routes, key=lambda route: (metrics[route]["top3_accuracy"], metrics[route]["top1_accuracy"]))
    base = metrics[base_route]
    best = metrics[best_top1]
    decision = (
        f"{best_top1} improves exact Top1 over Top200 without lowering Top3; use as next confirmation candidate."
        if best["top1_accuracy"] > base["top1_accuracy"] and best["top3_accuracy"] >= base["top3_accuracy"]
        else "No gene-count variant clearly improves exact Top1 without a Top3 tradeoff; keep Top200 baseline for now."
    )
    summary: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; same Network Top3 beam; fold-local Region gene-count variants",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singletons).sum()),
        "gene_counts": gene_counts,
        "fusion_weights": FUSION_WEIGHTS,
        "baseline_route": base_route,
        "best_top1_route": best_top1,
        "best_top3_route": best_top3,
        "routes": metrics,
        "paired_changes_vs_top200": changes,
        "paired_pvalues_vs_top200": pvalues,
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_gene_count_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics, all_routes)
    write_report(args.outdir, summary, all_routes)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
