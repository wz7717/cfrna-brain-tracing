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

from scripts.run_bo2023_group_specific_correlation_validation import resolution_group_training  # noqa: E402
from scripts.run_bo2023_hierarchical_region_correlation_validation import build_local_discriminative_rows, rank_candidates  # noqa: E402
from scripts.run_bo2023_loso_validation import build_region_reference, correlation_scores, read_annotations, read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import (  # noqa: E402
    BEAM_ROUTE,
    DEFAULT_NETWORK_DETAIL,
    build_resolution_groups,
    candidate_training_indices,
    distinct_ranked_groups,
    region_network_assignment,
    score_route,
)
from scripts.run_bo2023_v2_loso_validation import DEFAULT_GENE_MAP, DEFAULT_MATRIX, DEFAULT_SAMPLE_INFO, map_matrix_to_symbols  # noqa: E402


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_group_pairwise_rescue_nested_loso_814_20260527"
PAIR_TOP3_ROUTE = "top3_network_beam_group_pairwise_rescue_top3"
PAIR_TOP5_ROUTE = "top3_network_beam_group_pairwise_rescue_top5"


def pair_key(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((left, right)))


def summarize(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "group_top1_hits": int(detail["group_hit1"].sum()),
        "group_top1_accuracy": float(detail["group_hit1"].mean()),
        "group_top3_hits": int(detail["group_hit3"].sum()),
        "group_top3_accuracy": float(detail["group_hit3"].mean()),
        "median_group_true_rank": float(detail["group_true_rank"].median()),
        "n_switches": int(detail.get("switched", pd.Series(0, index=detail.index)).sum()),
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


def derive_training_group_confusion_pairs(
    values: np.ndarray,
    candidates: list[str],
    training: dict[str, np.ndarray],
    annotations: dict[str, dict[str, Any]],
    local_rows: np.ndarray,
    max_pairs_per_truth: int,
    min_pair_errors: int,
    min_group_samples: int,
) -> tuple[set[tuple[str, str]], pd.DataFrame]:
    """Discover repeated group-level mistakes using outer-fold training samples only."""
    group_training = resolution_group_training(annotations, training)
    counts: Counter[tuple[str, str]] = Counter()
    reference = np.column_stack(
        [values[local_rows[:, None], training[region]].mean(axis=1, dtype=np.float64) for region in candidates]
    )
    for truth_region in candidates:
        truth_group = str(annotations[truth_region]["resolution_group"])
        if len(group_training[truth_group]) < min_group_samples or len(training[truth_region]) <= 1:
            continue
        truth_position = candidates.index(truth_region)
        for sample_idx in training[truth_region]:
            local_reference = reference.copy()
            remaining = training[truth_region][training[truth_region] != sample_idx]
            local_reference[:, truth_position] = values[local_rows[:, None], remaining].mean(axis=1, dtype=np.float64)
            scores = correlation_scores(local_reference, values[local_rows, sample_idx])
            ranked_regions = [candidates[int(j)] for j in np.argsort(scores)[::-1]]
            predicted_group = distinct_ranked_groups(ranked_regions, annotations)[0]
            if predicted_group != truth_group:
                counts[(truth_group, predicted_group)] += 1

    selected: set[tuple[str, str]] = set()
    rows: list[dict[str, Any]] = []
    for truth_group in sorted(group_training):
        candidates_for_truth = sorted(
            [(predicted, count) for (truth, predicted), count in counts.items() if truth == truth_group],
            key=lambda item: (-item[1], item[0]),
        )
        retained = 0
        for predicted, count in candidates_for_truth:
            sufficient = (
                count >= min_pair_errors
                and retained < max_pairs_per_truth
                and len(group_training.get(truth_group, [])) >= min_group_samples
                and len(group_training.get(predicted, [])) >= min_group_samples
            )
            rows.append(
                {
                    "truth_group": truth_group,
                    "confused_as_group": predicted,
                    "error_count": int(count),
                    "truth_training_samples": int(len(group_training.get(truth_group, []))),
                    "predicted_training_samples": int(len(group_training.get(predicted, []))),
                    "selected": bool(sufficient),
                }
            )
            if sufficient:
                selected.add(pair_key(truth_group, predicted))
                retained += 1
    return selected, pd.DataFrame(rows)


def build_pair_models(
    values: np.ndarray,
    group_training: dict[str, np.ndarray],
    pairs: set[tuple[str, str]],
    gene_pool: np.ndarray,
    top_n_genes: int,
) -> tuple[dict[tuple[str, str], tuple[np.ndarray, np.ndarray]], pd.DataFrame]:
    models: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    audit: list[dict[str, Any]] = []
    for left, right in sorted(pairs):
        left_values = values[gene_pool[:, None], group_training[left]].astype(float, copy=False)
        right_values = values[gene_pool[:, None], group_training[right]].astype(float, copy=False)
        delta = left_values.mean(axis=1) - right_values.mean(axis=1)
        pooled = (left_values.var(axis=1, ddof=1) + right_values.var(axis=1, ddof=1)) / 2.0
        effect = np.abs(delta) / np.sqrt(pooled + 1e-8)
        selected_rows = np.argsort(effect)[::-1][: min(top_n_genes, len(gene_pool))]
        genes = gene_pool[selected_rows].astype(int)
        reference = np.column_stack(
            [
                values[genes[:, None], group_training[left]].mean(axis=1, dtype=np.float64),
                values[genes[:, None], group_training[right]].mean(axis=1, dtype=np.float64),
            ]
        )
        models[(left, right)] = (genes, reference)
        audit.append(
            {
                "left_group": left,
                "right_group": right,
                "left_training_samples": int(len(group_training[left])),
                "right_training_samples": int(len(group_training[right])),
                "n_genes": int(len(genes)),
                "mean_selected_effect": float(effect[selected_rows].mean()),
                "centroid_pair_correlation": float(correlation_scores(reference, reference[:, 0])[1]),
            }
        )
    return models, pd.DataFrame(audit)


def rescue_detail(
    route: str,
    base_detail: dict[str, Any],
    baseline_groups: list[str],
    pair_models: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]],
    sample: np.ndarray,
    candidate_k: int,
    min_margin: float,
) -> dict[str, Any]:
    ranked = baseline_groups.copy()
    anchor = ranked[0]
    selected_position = -1
    selected_margin = float(min_margin)
    selected_pair = ""
    for position in range(1, min(candidate_k, len(ranked))):
        challenger = ranked[position]
        key = pair_key(anchor, challenger)
        model = pair_models.get(key)
        if model is None:
            continue
        genes, reference = model
        pair_scores = correlation_scores(reference, sample[genes])
        if key[0] == anchor:
            margin = float(pair_scores[1] - pair_scores[0])
        else:
            margin = float(pair_scores[0] - pair_scores[1])
        if margin > selected_margin:
            selected_position = position
            selected_margin = margin
            selected_pair = f"{key[0]} <> {key[1]}"
    switched = selected_position >= 1
    if switched:
        ranked[0], ranked[selected_position] = ranked[selected_position], ranked[0]
    truth_group = str(base_detail["true_resolution_group"])
    true_rank = ranked.index(truth_group) + 1 if truth_group in ranked else int(base_detail["n_candidate_groups"]) + 1
    padded = ranked[:3] + [""] * max(0, 3 - len(ranked))
    detail = dict(base_detail)
    detail.update(
        {
            "route": route,
            "pred_group_top1": padded[0],
            "pred_group_top2": padded[1],
            "pred_group_top3": padded[2],
            "group_true_rank": int(true_rank),
            "group_hit1": int(true_rank == 1),
            "group_hit3": int(true_rank <= 3),
            "switched": int(switched),
            "switch_pair": selected_pair,
            "switch_margin": float(selected_margin if switched else 0.0),
        }
    )
    return detail


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [BEAM_ROUTE, PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE]
    labels = ["Baseline", "Pairwise Top3", "Pairwise Top5"]
    top1 = [metrics[route]["group_top1_accuracy"] for route in routes]
    top3 = [metrics[route]["group_top3_accuracy"] for route in routes]
    x = np.arange(len(routes))
    width = 0.38
    fig, ax = plt.subplots(figsize=(9.4, 5.6), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Top3", color="#009E73")
    ax.set_ylabel("Resolvable group accuracy")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, min(1.0, max(top3) + 0.14))
    ax.set_title("Bo2023 pair-specific group rescue: strict outer LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "group_pairwise_rescue_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "group_pairwise_rescue_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    base = summary["routes"][BEAM_ROUTE]
    top3 = summary["routes"][PAIR_TOP3_ROUTE]
    top5 = summary["routes"][PAIR_TOP5_ROUTE]
    change3 = summary["paired_changes"][PAIR_TOP3_ROUTE]
    change5 = summary["paired_changes"][PAIR_TOP5_ROUTE]
    text = f"""# 高频混淆 Resolvable Group Pair-specific Correlation 严格 LOSO 验证

## 设计

- 基线固定为当前最优 `Top3 Network beam -> Region-first resolvable group`，Group Top1 `{base['group_top1_accuracy']:.1%}`、Top3 `{base['group_top3_accuracy']:.1%}`。
- 每个外层 LOSO fold 内，先使用训练样本生成 resolution group 与局部 Region correlation 基线；随后仅在训练样本内部统计真实 group 被误判为其他 group 的重复错误对。
- 保留每个真实 group 最多 `{summary['max_pairs_per_truth']}` 个、训练错误不少于 `{summary['min_pair_errors']}` 次、每侧训练样本不少于 `{summary['min_group_samples']}` 个的高频混淆对。
- Pair-specific 模型仅针对被选错误组对，从训练折候选基因池中选择 Top `{summary['pair_top_n_genes']}` 个效应量基因。
- `Pairwise Top3` 仅允许当前 Top3 group 内挑战 Top1；`Pairwise Top5` 允许 Top5 内挑战，从而可能改变 Top3 集合。

## 结果

| 路径 | Group Top1 | Group Top3 | 切换次数 |
| --- | ---: | ---: | ---: |
| 当前基线 | {base['group_top1_hits']}/{summary['n_test_samples']} ({base['group_top1_accuracy']:.1%}) | {base['group_top3_hits']}/{summary['n_test_samples']} ({base['group_top3_accuracy']:.1%}) | 0 |
| Pairwise Top3 | {top3['group_top1_hits']}/{summary['n_test_samples']} ({top3['group_top1_accuracy']:.1%}) | {top3['group_top3_hits']}/{summary['n_test_samples']} ({top3['group_top3_accuracy']:.1%}) | {top3['n_switches']} |
| Pairwise Top5 | {top5['group_top1_hits']}/{summary['n_test_samples']} ({top5['group_top1_accuracy']:.1%}) | {top5['group_top3_hits']}/{summary['n_test_samples']} ({top5['group_top3_accuracy']:.1%}) | {top5['n_switches']} |

## 配对变化

| 路径 | Top1 新增 / 丢失 | Top1 p | Top3 新增 / 丢失 | Top3 p |
| --- | ---: | ---: | ---: | ---: |
| Pairwise Top3 | {change3['top1_gains']} / {change3['top1_losses']} | {summary['paired_pvalues'][PAIR_TOP3_ROUTE]['top1']:.3f} | {change3['top3_gains']} / {change3['top3_losses']} | {summary['paired_pvalues'][PAIR_TOP3_ROUTE]['top3']:.3f} |
| Pairwise Top5 | {change5['top1_gains']} / {change5['top1_losses']} | {summary['paired_pvalues'][PAIR_TOP5_ROUTE]['top1']:.3f} | {change5['top3_gains']} / {change5['top3_losses']} | {summary['paired_pvalues'][PAIR_TOP5_ROUTE]['top3']:.3f} |

## 判定

{summary['decision']}
"""
    (outdir / "group_pairwise_rescue_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Fold-local pair-specific rescue of repeatedly confused resolution groups.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--pair-gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=100)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=3)
    parser.add_argument("--min-group-samples", type=int, default=3)
    parser.add_argument("--pair-min-margin", type=float, default=0.0)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--resolution-min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.20)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=4)
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
    counts = ann.groupby("region_id")["sample_id"].size()
    singleton_regions = set(counts[counts < 2].index)
    selected = ann[~ann["region_id"].isin(singleton_regions)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)
    if args.max_samples is not None:
        selected = selected.head(max(1, int(args.max_samples))).copy()

    rows: dict[str, list[dict[str, Any]]] = {BEAM_ROUTE: [], PAIR_TOP3_ROUTE: [], PAIR_TOP5_ROUTE: []}
    confusion_audit: list[pd.DataFrame] = []
    model_audit: list[pd.DataFrame] = []
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
        local_rows = build_local_discriminative_rows(
            values, candidates, region_indices, heldout_idx, args.local_top_n_genes
        )
        if local_rows is None:
            local_rows = np.arange(values.shape[0], dtype=int)
        pair_pool = build_local_discriminative_rows(
            values, candidates, region_indices, heldout_idx, args.pair_gene_pool_size
        )
        if pair_pool is None:
            pair_pool = np.arange(values.shape[0], dtype=int)
        training = candidate_training_indices(candidates, region_indices, heldout_idx)
        annotations, _ = build_resolution_groups(
            values,
            candidates,
            training,
            region_network_assignment(training_ann, candidates),
            local_rows,
            args.min_resolution_samples,
            args.min_merge_samples,
            args.resolution_min_pair_errors,
            args.min_confusion_rate,
            args.similarity_threshold,
            args.merge_similarity_threshold,
            args.max_group_size,
        )
        selected_pairs, errors = derive_training_group_confusion_pairs(
            values,
            candidates,
            training,
            annotations,
            local_rows,
            args.max_pairs_per_truth,
            args.min_pair_errors,
            args.min_group_samples,
        )
        if not errors.empty:
            errors.insert(0, "fold", fold)
            errors.insert(1, "sample_id", sample_id)
            confusion_audit.append(errors)
        group_training = resolution_group_training(annotations, training)
        pair_models, audit = build_pair_models(
            values, group_training, selected_pairs, pair_pool, args.pair_top_n_genes
        )
        if not audit.empty:
            audit.insert(0, "fold", fold)
            audit.insert(1, "sample_id", sample_id)
            model_audit.append(audit)

        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)
        baseline_scores = correlation_scores(reference, sample, local_rows)
        candidate_indices = np.asarray([region_pos[region] for region in candidates], dtype=int)
        ranked_regions = rank_candidates(baseline_scores, regions, candidate_indices)
        base = score_route(
            BEAM_ROUTE, sample_id, truth_region, truth_network, network_top, ranked_regions, annotations, len(regions)
        )
        base["switched"] = 0
        base["switch_pair"] = ""
        base["switch_margin"] = 0.0
        rows[BEAM_ROUTE].append(base)
        baseline_groups = distinct_ranked_groups(ranked_regions, annotations)
        rows[PAIR_TOP3_ROUTE].append(
            rescue_detail(PAIR_TOP3_ROUTE, base, baseline_groups, pair_models, sample, 3, args.pair_min_margin)
        )
        rows[PAIR_TOP5_ROUTE].append(
            rescue_detail(PAIR_TOP5_ROUTE, base, baseline_groups, pair_models, sample, 5, args.pair_min_margin)
        )
        fold_rows.append(
            {
                "fold": fold,
                "sample_id": sample_id,
                "truth_region": truth_region,
                "n_candidate_regions": int(len(candidates)),
                "n_candidate_groups": int(len(group_training)),
                "n_selected_pairs": int(len(selected_pairs)),
                "n_pair_models": int(len(pair_models)),
            }
        )

    details = {route: pd.DataFrame(data) for route, data in rows.items()}
    metrics = {route: summarize(detail) for route, detail in details.items()}
    changes = {route: paired_changes(details[BEAM_ROUTE], details[route]) for route in [PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE]}
    pvalues = {
        route: {
            metric: paired_pvalue(changes[route][f"{metric}_gains"], changes[route][f"{metric}_losses"])
            for metric in ["top1", "top3"]
        }
        for route in changes
    }
    base = metrics[BEAM_ROUTE]
    best_route = max([PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE], key=lambda route: (
        metrics[route]["group_top1_accuracy"], metrics[route]["group_top3_accuracy"]
    ))
    best = metrics[best_route]
    improved = (
        best["group_top1_accuracy"] > base["group_top1_accuracy"]
        and best["group_top3_accuracy"] >= base["group_top3_accuracy"]
    )
    if improved:
        decision = (
            f"`{best_route}` 在当前内部严格 LOSO 中提高 Group Top1 且未降低 Group Top3；"
            "可进入固定规则确认轮，但在独立确认前不替换正式路径。"
        )
    else:
        decision = (
            "高频混淆组对的 pair-specific correlation 未在 Group Top1 提升且 Group Top3 不下降约束下超过当前基线，"
            "不接入正式路径。"
        )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO with fold-local resolution groups, training-confusion pair discovery and pair-specific group correlation",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singleton_regions).sum()),
        "local_top_n_genes": int(args.local_top_n_genes),
        "pair_gene_pool_size": int(args.pair_gene_pool_size),
        "pair_top_n_genes": int(args.pair_top_n_genes),
        "max_pairs_per_truth": int(args.max_pairs_per_truth),
        "min_pair_errors": int(args.min_pair_errors),
        "min_group_samples": int(args.min_group_samples),
        "pair_min_margin": float(args.pair_min_margin),
        "routes": metrics,
        "paired_changes": changes,
        "paired_pvalues": pvalues,
        "best_route": best_route,
        "meets_internal_adoption_rule": bool(improved),
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_detail.csv", index=False, encoding="utf-8-sig")
    if confusion_audit:
        pd.concat(confusion_audit, ignore_index=True).to_csv(
            args.outdir / "fold_training_group_confusion_pairs.csv", index=False, encoding="utf-8-sig"
        )
    if model_audit:
        pd.concat(model_audit, ignore_index=True).to_csv(
            args.outdir / "fold_group_pair_model_audit.csv", index=False, encoding="utf-8-sig"
        )
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
