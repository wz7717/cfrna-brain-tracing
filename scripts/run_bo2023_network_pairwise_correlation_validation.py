#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
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

from benchmark.metrics import compute_multiclass_auc  # noqa: E402
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


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_network_pairwise_correlation_full_loso_819"
BASELINE_ROUTE = "network_discriminative_correlation_top200"
PAIR_TOP3_ROUTE = "network_pairwise_correlation_rescue_top3"
PAIR_TOP5_ROUTE = "network_pairwise_correlation_rescue_top5"


def pair_key(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def derive_training_confusion_pairs(
    values: np.ndarray,
    labels: np.ndarray,
    groups: list[str],
    training: dict[str, np.ndarray],
    reference: np.ndarray,
    global_genes: np.ndarray,
    max_pairs_per_truth: int,
    min_pair_errors: int,
) -> tuple[set[tuple[str, str]], pd.DataFrame]:
    """Find repeatedly confused Network pairs using only the outer-fold training samples."""
    group_pos = {group: i for i, group in enumerate(groups)}
    selected_reference = reference[global_genes, :].astype(float, copy=False)
    counts: Counter[tuple[str, str]] = Counter()
    for truth in groups:
        truth_pos = group_pos[truth]
        truth_training = training[truth]
        if len(truth_training) <= 1:
            continue
        for sample_idx in truth_training:
            local_reference = selected_reference.copy()
            local_reference[:, truth_pos] = (
                selected_reference[:, truth_pos] * len(truth_training) - values[global_genes, sample_idx]
            ) / (len(truth_training) - 1)
            scores = correlation_scores(local_reference, values[global_genes, sample_idx])
            scores[truth_pos] = -np.inf
            predicted = groups[int(np.argmax(scores))]
            counts[(truth, predicted)] += 1

    rows: list[dict[str, object]] = []
    chosen: set[tuple[str, str]] = set()
    for truth in groups:
        truth_rows = sorted(
            [(predicted, count) for (actual, predicted), count in counts.items() if actual == truth],
            key=lambda item: (-item[1], item[0]),
        )
        retained = 0
        for predicted, count in truth_rows:
            selected = count >= min_pair_errors and retained < max_pairs_per_truth
            rows.append(
                {
                    "truth_network": truth,
                    "confused_as_network": predicted,
                    "error_count": int(count),
                    "selected": bool(selected),
                }
            )
            if selected:
                chosen.add(pair_key(truth, predicted))
                retained += 1
    return chosen, pd.DataFrame(rows)


def build_pair_models(
    values: np.ndarray,
    training: dict[str, np.ndarray],
    reference: np.ndarray,
    pairs: set[tuple[str, str]],
    groups: list[str],
    gene_pool: np.ndarray,
    top_n_genes: int,
) -> tuple[dict[tuple[str, str], np.ndarray], pd.DataFrame]:
    """Select fold-local high-effect genes for each selected confusion pair."""
    group_pos = {group: i for i, group in enumerate(groups)}
    models: dict[tuple[str, str], np.ndarray] = {}
    audit_rows: list[dict[str, object]] = []
    for left, right in sorted(pairs):
        left_values = values[gene_pool[:, None], training[left]]
        right_values = values[gene_pool[:, None], training[right]]
        mean_difference = left_values.mean(axis=1) - right_values.mean(axis=1)
        pooled_variance = (left_values.var(axis=1, ddof=1) + right_values.var(axis=1, ddof=1)) / 2.0
        effect = np.abs(mean_difference) / np.sqrt(pooled_variance + 1e-8)
        local_rows = np.argsort(effect)[::-1][: min(int(top_n_genes), len(gene_pool))]
        selected = gene_pool[local_rows].astype(int)
        models[(left, right)] = selected
        pair_reference = reference[selected, :][:, [group_pos[left], group_pos[right]]]
        audit_rows.append(
            {
                "left_network": left,
                "right_network": right,
                "n_genes": int(len(selected)),
                "mean_selected_effect": float(effect[local_rows].mean()),
                "centroid_pair_correlation": float(correlation_scores(pair_reference, pair_reference[:, 0])[1]),
            }
        )
    return models, pd.DataFrame(audit_rows)


def baseline_scores(reference: np.ndarray, sample: np.ndarray, genes: np.ndarray) -> np.ndarray:
    return correlation_scores(reference, sample, genes)


def make_detail(
    route: str,
    sample_id: str,
    truth: str,
    groups: list[str],
    ranked_positions: list[int],
    scores: np.ndarray,
    switched: bool = False,
    switch_pair: str = "",
    switch_margin: float = 0.0,
) -> tuple[dict[str, Any], dict[str, Any]]:
    ranked = [groups[j] for j in ranked_positions]
    true_rank = ranked.index(truth) + 1
    adjusted_scores = scores.copy()
    natural_order = np.argsort(scores)[::-1].tolist()
    if ranked_positions != natural_order:
        ordered_values = np.sort(scores)[::-1]
        adjusted_scores = np.empty_like(scores)
        for rank, position in enumerate(ranked_positions):
            adjusted_scores[position] = ordered_values[rank]
    probs = softmax(adjusted_scores)
    detail = {
        "route": route,
        "sample_id": sample_id,
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
    probability: dict[str, Any] = {"sample_id": sample_id, "label": truth}
    probability.update({group: float(probs[j]) for j, group in enumerate(groups)})
    return detail, probability


def evaluate_pairwise_rescue(
    route: str,
    sample_id: str,
    truth: str,
    sample: np.ndarray,
    reference: np.ndarray,
    groups: list[str],
    scores: np.ndarray,
    pair_models: dict[tuple[str, str], np.ndarray],
    candidate_k: int,
    min_margin: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    order = np.argsort(scores)[::-1].tolist()
    anchor = order[0]
    group_pos = {group: i for i, group in enumerate(groups)}
    best_position = -1
    best_margin = float(min_margin)
    best_pair = ""
    for position in range(1, min(candidate_k, len(order))):
        challenger = order[position]
        key = pair_key(groups[anchor], groups[challenger])
        genes = pair_models.get(key)
        if genes is None:
            continue
        pair_scores = correlation_scores(reference, sample, genes)
        margin = float(pair_scores[challenger] - pair_scores[anchor])
        if margin > best_margin:
            best_margin = margin
            best_position = position
            best_pair = f"{key[0]} <> {key[1]}"
    switched = best_position >= 1
    if switched:
        order[0], order[best_position] = order[best_position], order[0]
    return make_detail(
        route,
        sample_id,
        truth,
        groups,
        order,
        scores,
        switched=switched,
        switch_pair=best_pair,
        switch_margin=best_margin if switched else 0.0,
    )


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


def paired_change(base: pd.DataFrame, tested: pd.DataFrame) -> dict[str, int]:
    return {
        "top1_gains": int(((base["hit1"] == 0) & (tested["hit1"] == 1)).sum()),
        "top1_losses": int(((base["hit1"] == 1) & (tested["hit1"] == 0)).sum()),
        "top3_gains": int(((base["hit3"] == 0) & (tested["hit3"] == 1)).sum()),
        "top3_losses": int(((base["hit3"] == 1) & (tested["hit3"] == 0)).sum()),
    }


def export_plot(outdir: Path, metrics: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["Top200 baseline", "Pairwise Top3", "Pairwise Top5"]
    routes = [BASELINE_ROUTE, PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE]
    values = metrics.set_index("route").loc[routes, ["top1_accuracy", "top3_accuracy", "macro_auc"]].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), constrained_layout=True)
    for j, title in enumerate(["Top1 accuracy", "Top3 accuracy", "Macro AUC"]):
        bars = axes[j].bar(labels, values[:, j], color=["#0072B2", "#D55E00", "#009E73"])
        axes[j].set_title(title, fontweight="bold")
        axes[j].grid(axis="y", alpha=0.25)
        axes[j].set_axisbelow(True)
        axes[j].tick_params(axis="x", rotation=17)
        axes[j].set_ylim(0, min(1.0, float(values[:, j].max()) + 0.14))
        for bar, value in zip(bars, values[:, j]):
            label = f"{value:.1%}" if j < 2 else f"{value:.3f}"
            axes[j].text(bar.get_x() + bar.get_width() / 2, value + 0.02, label, ha="center", fontweight="bold")
    fig.suptitle("Bo2023 Network pairwise discriminative correlation: full LOSO", fontweight="bold")
    fig.savefig(outdir / "network_pairwise_correlation_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "network_pairwise_correlation_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    base = summary["routes"][BASELINE_ROUTE]
    top3 = summary["routes"][PAIR_TOP3_ROUTE]
    top5 = summary["routes"][PAIR_TOP5_ROUTE]
    change3 = summary["paired_changes"][PAIR_TOP3_ROUTE]
    change5 = summary["paired_changes"][PAIR_TOP5_ROUTE]
    text = f"""# Network 高频混淆对专用判别基因 Correlation 优化验证

## 设计

- 终点：`SaleemNetworks`（{summary['n_classes']} 类）；数据范围：Bo2023 全部 `{summary['n_test_samples']}` 个样本。
- 验证：严格外层 LOSO；每一折的测试样本不参与 reference、全局判别基因、训练混淆对或 pair 专用基因选择。
- 当前正式基线：训练折内 Top `{summary['global_top_n_genes']}` 判别基因 Pearson correlation。
- 优化：仅在训练折内部按基线的 leave-one-out 误判频率选取每个真实 Network 最多 `{summary['max_pairs_per_truth']}` 个高频混淆对，最低误判数 `{summary['min_pair_errors']}`；每对从 Top `{summary['gene_pool_size']}` 全局判别基因池中选择 Top `{summary['pair_top_n_genes']}` 专用基因。
- `Pairwise Top3`：仅允许基线 Top3 内的高频混淆候选挑战 Top1；Top3 候选集合不变。
- `Pairwise Top5`：允许基线 Top5 内候选挑战 Top1，可能改变 Top3 候选集合。

## 结果

| 路径 | Top1 | Top3 | Macro AUC | 发生校正 |
| --- | ---: | ---: | ---: | ---: |
| 正式 Top200 correlation 基线 | {base['top1_hits']}/{summary['n_test_samples']} ({base['top1_accuracy']:.1%}) | {base['top3_hits']}/{summary['n_test_samples']} ({base['top3_accuracy']:.1%}) | {base['macro_auc']:.3f} | - |
| 高频混淆对 Pairwise Top3 | {top3['top1_hits']}/{summary['n_test_samples']} ({top3['top1_accuracy']:.1%}) | {top3['top3_hits']}/{summary['n_test_samples']} ({top3['top3_accuracy']:.1%}) | {top3['macro_auc']:.3f} | {top3['n_switches']} |
| 高频混淆对 Pairwise Top5 | {top5['top1_hits']}/{summary['n_test_samples']} ({top5['top1_accuracy']:.1%}) | {top5['top3_hits']}/{summary['n_test_samples']} ({top5['top3_accuracy']:.1%}) | {top5['macro_auc']:.3f} | {top5['n_switches']} |

## 配对变化

| 路径 | Top1 新增命中 | Top1 丢失命中 | Top3 新增命中 | Top3 丢失命中 |
| --- | ---: | ---: | ---: | ---: |
| Pairwise Top3 | {change3['top1_gains']} | {change3['top1_losses']} | {change3['top3_gains']} | {change3['top3_losses']} |
| Pairwise Top5 | {change5['top1_gains']} | {change5['top1_losses']} | {change5['top3_gains']} | {change5['top3_losses']} |

## 判定

{summary['decision']}

## 限制

本次方法设计来自此前全量错误结构观察，因此虽然每一折计算没有使用当前测试样本，本轮结果仍属于回顾性严格 LOSO 比较；若准备替换正式路径，需要在新增独立样本或外部图谱上再确认一次。
"""
    (outdir / "network_pairwise_correlation_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Full LOSO validation of pair-specific Network correlation rescue.")
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
    parser.add_argument("--pair-min-margin", type=float, default=0.002)
    parser.add_argument("--n-samples", type=int, default=819)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet)
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.endpoint].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    samples = matrix.columns.astype(str).tolist()
    if args.n_samples != len(samples):
        raise ValueError("This validation route is defined for full LOSO; --n-samples must match all matrix samples.")
    labels = ann.set_index("sample_id").reindex(samples)["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    values = matrix.to_numpy(dtype=np.float32)
    gene_symbols = matrix.index.astype(str).to_numpy()

    route_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in [BASELINE_ROUTE, PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE]}
    probability_rows: dict[str, list[dict[str, Any]]] = {route: [] for route in route_rows}
    fold_pair_rows: list[pd.DataFrame] = []
    pair_model_rows: list[pd.DataFrame] = []
    for fold_no, sample_id in enumerate(samples, start=1):
        heldout_idx = fold_no - 1
        truth = labels[heldout_idx]
        reference, training = build_group_reference(values, labels, groups, heldout_idx)
        gene_pool, _ = select_group_discriminative_genes(values, groups, training, args.gene_pool_size)
        global_genes = gene_pool[: args.global_top_n_genes]
        pairs, pair_errors = derive_training_confusion_pairs(
            values,
            labels,
            groups,
            training,
            reference,
            global_genes,
            args.max_pairs_per_truth,
            args.min_pair_errors,
        )
        pair_errors.insert(0, "fold", fold_no)
        pair_errors.insert(1, "sample_id", sample_id)
        fold_pair_rows.append(pair_errors)
        pair_models, pair_audit = build_pair_models(
            values,
            training,
            reference,
            pairs,
            groups,
            gene_pool,
            args.pair_top_n_genes,
        )
        pair_audit.insert(0, "fold", fold_no)
        pair_audit.insert(1, "sample_id", sample_id)
        pair_model_rows.append(pair_audit)
        scores = baseline_scores(reference, values[:, heldout_idx], global_genes)
        base_order = np.argsort(scores)[::-1].tolist()
        base_detail, base_probability = make_detail(BASELINE_ROUTE, sample_id, truth, groups, base_order, scores)
        route_rows[BASELINE_ROUTE].append(base_detail)
        probability_rows[BASELINE_ROUTE].append(base_probability)
        for route, candidate_k in [(PAIR_TOP3_ROUTE, 3), (PAIR_TOP5_ROUTE, 5)]:
            detail, probability = evaluate_pairwise_rescue(
                route,
                sample_id,
                truth,
                values[:, heldout_idx],
                reference,
                groups,
                scores,
                pair_models,
                candidate_k,
                args.pair_min_margin,
            )
            route_rows[route].append(detail)
            probability_rows[route].append(probability)

    detail_frames = {route: pd.DataFrame(rows) for route, rows in route_rows.items()}
    probability_frames = {route: pd.DataFrame(rows) for route, rows in probability_rows.items()}
    metrics = {route: summarize(detail_frames[route], probability_frames[route], groups) for route in route_rows}
    paired = {
        route: paired_change(detail_frames[BASELINE_ROUTE], detail_frames[route])
        for route in [PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE]
    }
    base = metrics[BASELINE_ROUTE]
    best_route = max([PAIR_TOP3_ROUTE, PAIR_TOP5_ROUTE], key=lambda route: (metrics[route]["top1_accuracy"], metrics[route]["top3_accuracy"]))
    best = metrics[best_route]
    improved = best["top1_accuracy"] > base["top1_accuracy"] and best["top3_accuracy"] >= base["top3_accuracy"]
    decision = (
        f"`{best_route}` 在本轮回顾性全量严格 LOSO 中同时提高 Top1 且未降低 Top3，可进入独立队列确认；独立确认前不替换正式 Top200 路径。"
        if improved
        else "高频混淆对专用 correlation 未能在 Top1 提升且 Top3 不下降的约束下超过正式 Top200 路径，因此不接入正式输出。"
    )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "endpoint": args.endpoint,
        "validation_design": "retrospective full strict LOSO with fold-local high-confusion pair discovery and pair-specific correlation",
        "n_test_samples": int(len(samples)),
        "n_classes": int(len(groups)),
        "global_top_n_genes": int(args.global_top_n_genes),
        "gene_pool_size": int(args.gene_pool_size),
        "pair_top_n_genes": int(args.pair_top_n_genes),
        "max_pairs_per_truth": int(args.max_pairs_per_truth),
        "min_pair_errors": int(args.min_pair_errors),
        "pair_min_margin": float(args.pair_min_margin),
        "routes": metrics,
        "paired_changes": paired,
        "best_pairwise_route": best_route,
        "meets_adoption_rule": bool(improved),
        "decision": decision,
    }
    metrics_frame = pd.DataFrame([{"route": route, **values} for route, values in metrics.items()])
    metrics_frame.to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    for route, frame in detail_frames.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.concat(fold_pair_rows, ignore_index=True).to_csv(args.outdir / "fold_training_confusion_pairs.csv", index=False, encoding="utf-8-sig")
    pd.concat(pair_model_rows, ignore_index=True).to_csv(args.outdir / "fold_pair_model_audit.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics_frame)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
