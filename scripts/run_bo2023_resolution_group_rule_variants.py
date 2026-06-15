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

from scripts.run_bo2023_hierarchical_region_correlation_validation import build_local_discriminative_rows, rank_candidates  # noqa: E402
from scripts.run_bo2023_loso_validation import build_region_reference, correlation_scores, read_annotations, read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import (  # noqa: E402
    DEFAULT_NETWORK_DETAIL,
    UnionFind,
    build_resolution_groups,
    candidate_training_indices,
    centroid_reference,
    region_network_assignment,
    score_route,
)
from scripts.run_bo2023_v2_loso_validation import DEFAULT_GENE_MAP, DEFAULT_MATRIX, DEFAULT_SAMPLE_INFO, map_matrix_to_symbols  # noqa: E402


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_resolution_group_rule_variants_loso_814_20260529"
BASELINE_ROUTE = "same_network_max4_baseline"
SAME_NETWORK_MAX6_ROUTE = "same_network_max6"
CROSS_NETWORK_MAX6_ROUTE = "top3_beam_cross_network_max6"


def build_cross_network_resolution_groups(
    values: np.ndarray,
    candidates: list[str],
    training: dict[str, np.ndarray],
    region_network: dict[str, str | None],
    rows: np.ndarray,
    min_resolution_samples: int,
    min_merge_samples: int,
    min_pair_errors: int,
    min_confusion_rate: float,
    similarity_threshold: float,
    merge_similarity_threshold: float,
    max_group_size: int,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    """Build fold-local groups allowing small cross-Network merges inside current Top3 beam."""
    reference = centroid_reference(values, candidates, training, rows)
    positions = {region: j for j, region in enumerate(candidates)}
    annotations: dict[str, dict[str, Any]] = {}
    for region in candidates:
        annotations[region] = {
            "region_id": region,
            "network_id": region_network.get(region),
            "resolution_group": region,
            "group_members": [region],
            "resolution_tier": "high_resolution",
            "resolution_reasons": [],
            "training_samples": int(len(training[region])),
            "nearest_centroid_corr": None,
        }
        if len(training[region]) < min_resolution_samples:
            annotations[region]["resolution_reasons"].append(f"training_n<{min_resolution_samples}")
        if region_network.get(region) is None:
            annotations[region]["resolution_reasons"].append("cross_network_region")

    # Nearest similarity is computed over the full Top3 beam candidate space.
    for j, region in enumerate(candidates):
        competitors = [i for i in range(len(candidates)) if i != j]
        nearest = float(np.max(correlation_scores(reference[:, competitors], reference[:, j]))) if competitors else float("nan")
        annotations[region]["nearest_centroid_corr"] = nearest
        if np.isfinite(nearest) and nearest >= similarity_threshold:
            annotations[region]["resolution_reasons"].append(f"nearest_centroid_corr>={similarity_threshold:.2f}")

    confusions: Counter[tuple[str, str]] = Counter()
    for truth in candidates:
        truth_indices = training[truth]
        if len(truth_indices) <= 1:
            continue
        truth_pos = positions[truth]
        for sample_idx in truth_indices:
            local_ref = reference.copy()
            local_ref[:, truth_pos] = values[rows[:, None], truth_indices[truth_indices != sample_idx]].mean(
                axis=1, dtype=np.float64
            )
            scores = correlation_scores(local_ref, values[rows, sample_idx])
            predicted = candidates[int(np.argmax(scores))]
            if predicted != truth:
                confusions[(truth, predicted)] += 1

    pair_rows: list[dict[str, Any]] = []
    uf = UnionFind(candidates)
    for left_i, left in enumerate(candidates):
        for right in candidates[left_i + 1 :]:
            errors = int(confusions[(left, right)] + confusions[(right, left)])
            pair_n = int(len(training[left]) + len(training[right]))
            rate = float(errors / pair_n) if pair_n else 0.0
            corr = float(correlation_scores(reference[:, [positions[left], positions[right]]], reference[:, positions[left]])[1])
            left_network = region_network.get(left)
            right_network = region_network.get(right)
            cross_network = bool(left_network != right_network)
            eligible = (
                len(training[left]) >= min_merge_samples
                and len(training[right]) >= min_merge_samples
                and errors >= min_pair_errors
                and rate >= min_confusion_rate
                and corr >= merge_similarity_threshold
            )
            merged = False
            if eligible:
                merged = uf.union(left, right, max_group_size)
            pair_rows.append(
                {
                    "left_region": left,
                    "right_region": right,
                    "left_network": left_network,
                    "right_network": right_network,
                    "cross_network": cross_network,
                    "pair_errors": errors,
                    "confusion_rate": rate,
                    "centroid_correlation": corr,
                    "eligible_for_merge": bool(eligible),
                    "merged": bool(merged and eligible),
                }
            )

    groups: dict[str, list[str]] = {}
    for region in candidates:
        groups.setdefault(uf.find(region), []).append(region)
    for members in groups.values():
        members = sorted(members)
        if len(members) <= 1:
            continue
        networks = sorted({str(region_network.get(region)) for region in members})
        prefix = "XNET" if len(networks) > 1 else networks[0]
        group_id = f"{prefix}::{' + '.join(members)}"
        reason = "fold_local_cross_network_merged_confusion_group" if len(networks) > 1 else "fold_local_merged_confusion_group"
        for region in members:
            annotations[region]["resolution_group"] = group_id
            annotations[region]["group_members"] = members
            annotations[region]["resolution_reasons"].append(reason)

    for annotation in annotations.values():
        reasons = list(dict.fromkeys(annotation["resolution_reasons"]))
        annotation["resolution_reasons"] = reasons
        annotation["resolution_tier"] = "low_resolution" if reasons else "high_resolution"
    return annotations, pd.DataFrame(pair_rows)


def summarize(detail: pd.DataFrame) -> dict[str, Any]:
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
        "median_group_true_rank": float(detail["group_true_rank"].median()),
        "mean_candidate_groups": float(detail["n_candidate_groups"].mean()),
    }


def paired_changes(base: pd.DataFrame, tested: pd.DataFrame) -> dict[str, int]:
    return {
        "top1_gains": int(((base["group_hit1"] == 0) & (tested["group_hit1"] == 1)).sum()),
        "top1_losses": int(((base["group_hit1"] == 1) & (tested["group_hit1"] == 0)).sum()),
        "top3_gains": int(((base["group_hit3"] == 0) & (tested["group_hit3"] == 1)).sum()),
        "top3_losses": int(((base["group_hit3"] == 1) & (tested["group_hit3"] == 0)).sum()),
    }


def paired_pvalue(gains: int, losses: int) -> float:
    n = gains + losses
    if n == 0:
        return 1.0
    tail = min(gains, losses)
    probability = sum(math.comb(n, k) for k in range(tail + 1)) / (2**n)
    return float(min(1.0, 2.0 * probability))


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [BASELINE_ROUTE, SAME_NETWORK_MAX6_ROUTE, CROSS_NETWORK_MAX6_ROUTE]
    labels = ["Baseline\nsame max4", "Same Network\nmax6", "Top3 beam\ncross max6"]
    top1 = [metrics[route]["group_top1_accuracy"] for route in routes]
    top3 = [metrics[route]["group_top3_accuracy"] for route in routes]
    x = np.arange(len(routes))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9.8, 5.6), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Top3", color="#009E73")
    ax.set_ylabel("Resolvable group accuracy")
    ax.set_ylim(0, min(1.0, max(top3) + 0.14))
    ax.set_xticks(x, labels)
    ax.set_title("Bo2023 resolution group rule variants: strict LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "resolution_group_rule_variants_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "resolution_group_rule_variants_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    base = summary["routes"][BASELINE_ROUTE]
    same = summary["routes"][SAME_NETWORK_MAX6_ROUTE]
    cross = summary["routes"][CROSS_NETWORK_MAX6_ROUTE]
    change_same = summary["paired_changes"][SAME_NETWORK_MAX6_ROUTE]
    change_cross = summary["paired_changes"][CROSS_NETWORK_MAX6_ROUTE]
    text = f"""# Resolvable Group 合并规则变体严格 LOSO 验证

## 设计

- 严格 LOSO 框架不变：每个测试样本完全留出，分组规则、局部基因和 reference 均只由训练折生成。
- 基线：同 Network 内合并，最大 group size `4`，对应当前 `{base['group_top1_accuracy']:.1%}` / `{base['group_top3_accuracy']:.1%}`。
- 变体 A：同 Network 内合并，但最大 group size 放宽到 `6`。
- 变体 B：在当前 Top3 Network beam 候选空间内允许小范围跨 Network 合并，最大 group size `6`；仍要求训练折反复混淆、`centroid_corr >= {summary['merge_similarity_threshold']:.2f}`、`confusion_rate >= {summary['min_confusion_rate']:.2f}`。

## 结果

| 路线 | Group Top1 | Group Top3 | 平均候选 group 数 |
| --- | ---: | ---: | ---: |
| 基线 same-network max4 | {base['group_top1_hits']}/{summary['n_test_samples']} ({base['group_top1_accuracy']:.1%}) | {base['group_top3_hits']}/{summary['n_test_samples']} ({base['group_top3_accuracy']:.1%}) | {base['mean_candidate_groups']:.1f} |
| Same-network max6 | {same['group_top1_hits']}/{summary['n_test_samples']} ({same['group_top1_accuracy']:.1%}) | {same['group_top3_hits']}/{summary['n_test_samples']} ({same['group_top3_accuracy']:.1%}) | {same['mean_candidate_groups']:.1f} |
| Top3-beam cross-network max6 | {cross['group_top1_hits']}/{summary['n_test_samples']} ({cross['group_top1_accuracy']:.1%}) | {cross['group_top3_hits']}/{summary['n_test_samples']} ({cross['group_top3_accuracy']:.1%}) | {cross['mean_candidate_groups']:.1f} |

## 配对变化

| 路线 | Top1 新增 / 丢失 | Top1 p | Top3 新增 / 丢失 | Top3 p |
| --- | ---: | ---: | ---: | ---: |
| Same-network max6 | {change_same['top1_gains']} / {change_same['top1_losses']} | {summary['paired_pvalues'][SAME_NETWORK_MAX6_ROUTE]['top1']:.3f} | {change_same['top3_gains']} / {change_same['top3_losses']} | {summary['paired_pvalues'][SAME_NETWORK_MAX6_ROUTE]['top3']:.3f} |
| Top3-beam cross-network max6 | {change_cross['top1_gains']} / {change_cross['top1_losses']} | {summary['paired_pvalues'][CROSS_NETWORK_MAX6_ROUTE]['top1']:.3f} | {change_cross['top3_gains']} / {change_cross['top3_losses']} | {summary['paired_pvalues'][CROSS_NETWORK_MAX6_ROUTE]['top3']:.3f} |

## 判定

{summary['decision']}
"""
    (outdir / "resolution_group_rule_variants_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict LOSO comparison of resolution group merge-rule variants.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.20)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
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

    rows: dict[str, list[dict[str, Any]]] = {BASELINE_ROUTE: [], SAME_NETWORK_MAX6_ROUTE: [], CROSS_NETWORK_MAX6_ROUTE: []}
    fold_rows: list[dict[str, Any]] = []
    audit_rows: list[pd.DataFrame] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        network_top = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        candidates = sorted(training_ann.loc[training_ann["endpoint_label"].isin(network_top), "region_id"].unique().tolist())
        local_rows = build_local_discriminative_rows(values, candidates, region_indices, heldout_idx, args.local_top_n_genes)
        if local_rows is None:
            local_rows = np.arange(values.shape[0], dtype=int)
        training = candidate_training_indices(candidates, region_indices, heldout_idx)
        assignment = region_network_assignment(training_ann, candidates)
        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)
        scores = correlation_scores(reference, sample, local_rows)
        indices = np.asarray([region_pos[region] for region in candidates], dtype=int)
        ranked = rank_candidates(scores, regions, indices)

        variants = [
            (
                BASELINE_ROUTE,
                lambda: build_resolution_groups(
                    values, candidates, training, assignment, local_rows,
                    args.min_resolution_samples, args.min_merge_samples, args.min_pair_errors,
                    args.min_confusion_rate, args.similarity_threshold, args.merge_similarity_threshold, 4
                ),
            ),
            (
                SAME_NETWORK_MAX6_ROUTE,
                lambda: build_resolution_groups(
                    values, candidates, training, assignment, local_rows,
                    args.min_resolution_samples, args.min_merge_samples, args.min_pair_errors,
                    args.min_confusion_rate, args.similarity_threshold, args.merge_similarity_threshold, 6
                ),
            ),
            (
                CROSS_NETWORK_MAX6_ROUTE,
                lambda: build_cross_network_resolution_groups(
                    values, candidates, training, assignment, local_rows,
                    args.min_resolution_samples, args.min_merge_samples, args.min_pair_errors,
                    args.min_confusion_rate, args.similarity_threshold, args.merge_similarity_threshold, 6
                ),
            ),
        ]
        fold_info: dict[str, Any] = {"fold": fold, "sample_id": sample_id, "truth_region": truth_region}
        for route, builder in variants:
            annotations, audit = builder()
            rows[route].append(score_route(route, sample_id, truth_region, truth_network, network_top, ranked, annotations, len(regions)))
            if not audit.empty:
                audit.insert(0, "fold", fold)
                audit.insert(1, "sample_id", sample_id)
                audit.insert(2, "route", route)
                audit_rows.append(audit)
            groups = {item["resolution_group"] for item in annotations.values()}
            merged = {
                group for group in groups
                if any(len(value["group_members"]) > 1 and value["resolution_group"] == group for value in annotations.values())
            }
            cross_merged = {
                group for group in groups
                if any(
                    len(value["group_members"]) > 1
                    and value["resolution_group"] == group
                    and len({assignment.get(member) for member in value["group_members"]}) > 1
                    for value in annotations.values()
                )
            }
            fold_info[f"{route}_n_groups"] = int(len(groups))
            fold_info[f"{route}_n_merged_groups"] = int(len(merged))
            fold_info[f"{route}_n_cross_merged_groups"] = int(len(cross_merged))
        fold_rows.append(fold_info)

    details = {route: pd.DataFrame(data) for route, data in rows.items()}
    metrics = {route: summarize(detail) for route, detail in details.items()}
    changes = {route: paired_changes(details[BASELINE_ROUTE], details[route]) for route in [SAME_NETWORK_MAX6_ROUTE, CROSS_NETWORK_MAX6_ROUTE]}
    pvalues = {
        route: {metric: paired_pvalue(change[f"{metric}_gains"], change[f"{metric}_losses"]) for metric in ["top1", "top3"]}
        for route, change in changes.items()
    }
    base = metrics[BASELINE_ROUTE]
    best_route = max([SAME_NETWORK_MAX6_ROUTE, CROSS_NETWORK_MAX6_ROUTE], key=lambda route: (
        metrics[route]["group_top1_accuracy"], metrics[route]["group_top3_accuracy"]
    ))
    best = metrics[best_route]
    improved = best["group_top1_accuracy"] > base["group_top1_accuracy"] and best["group_top3_accuracy"] >= base["group_top3_accuracy"]
    if improved:
        decision = (
            f"`{best_route}` 在严格 LOSO 下提高 Group Top1 且未降低 Top3；"
            "可作为下一轮固定规则确认候选，但正式接入前仍需独立确认。"
        )
    else:
        decision = "两种小范围放宽合并规则均未在 Group Top1 提升且 Top3 不下降的约束下超过当前基线。"
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; merge-rule variants only; same region ranking reused across variants",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singletons).sum()),
        "local_top_n_genes": int(args.local_top_n_genes),
        "min_resolution_samples": int(args.min_resolution_samples),
        "min_merge_samples": int(args.min_merge_samples),
        "min_pair_errors": int(args.min_pair_errors),
        "min_confusion_rate": float(args.min_confusion_rate),
        "similarity_threshold": float(args.similarity_threshold),
        "merge_similarity_threshold": float(args.merge_similarity_threshold),
        "routes": metrics,
        "paired_changes": changes,
        "paired_pvalues": pvalues,
        "best_route": best_route,
        "meets_internal_adoption_rule": bool(improved),
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_rule_variant_summary.csv", index=False, encoding="utf-8-sig")
    if audit_rows:
        pd.concat(audit_rows, ignore_index=True).to_csv(args.outdir / "fold_rule_variant_merge_audit.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
