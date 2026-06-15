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

from scripts.run_bo2023_hierarchical_region_correlation_validation import build_local_discriminative_rows, rank_candidates  # noqa: E402
from scripts.run_bo2023_loso_validation import build_region_reference, correlation_scores, read_annotations, read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_region_soft_fusion_pairwise_nested_validation import (  # noqa: E402
    fused_scores,
    network_probabilities,
    region_network_prior,
)
from scripts.run_bo2023_resolution_tier_validation import (  # noqa: E402
    BEAM_ROUTE,
    DEFAULT_NETWORK_DETAIL,
    build_resolution_groups,
    candidate_training_indices,
    region_network_assignment,
    score_route,
)
from scripts.run_bo2023_v2_loso_validation import DEFAULT_GENE_MAP, DEFAULT_MATRIX, DEFAULT_SAMPLE_INFO, map_matrix_to_symbols  # noqa: E402


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_group_network_probability_fusion_nested_loso_814_20260527"


def route_name(alpha: float) -> str:
    return f"group_network_probability_fusion_alpha_{alpha:.2f}".replace(".", "p")


def summarize(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "group_top1_hits": int(detail["group_hit1"].sum()),
        "group_top1_accuracy": float(detail["group_hit1"].mean()),
        "group_top3_hits": int(detail["group_hit3"].sum()),
        "group_top3_accuracy": float(detail["group_hit3"].mean()),
        "median_group_true_rank": float(detail["group_true_rank"].median()),
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


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]], alphas: list[float]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [BEAM_ROUTE] + [route_name(alpha) for alpha in alphas]
    labels = ["Baseline"] + [f"alpha={alpha:.2f}" for alpha in alphas]
    top1 = [metrics[route]["group_top1_accuracy"] for route in routes]
    top3 = [metrics[route]["group_top3_accuracy"] for route in routes]
    x = np.arange(len(routes))
    width = 0.38
    fig, ax = plt.subplots(figsize=(12.0, 5.8), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Top3", color="#009E73")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Resolvable group accuracy")
    ax.set_ylim(0, min(1.0, max(top3) + 0.13))
    ax.set_title("Bo2023 Network-probability fusion on resolvable groups: nested LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f"{bar.get_height():.1%}", ha="center", fontsize=9)
    fig.savefig(outdir / "group_network_probability_fusion_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "group_network_probability_fusion_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any], alphas: list[float]) -> None:
    base = summary["routes"][BEAM_ROUTE]
    rows = []
    for alpha in alphas:
        route = route_name(alpha)
        metric = summary["routes"][route]
        change = summary["paired_changes"][route]
        pvalue = summary["paired_pvalues"][route]
        rows.append(
            f"| `{alpha:.2f}` | {metric['group_top1_hits']}/{summary['n_test_samples']} ({metric['group_top1_accuracy']:.1%}) "
            f"| {metric['group_top3_hits']}/{summary['n_test_samples']} ({metric['group_top3_accuracy']:.1%}) "
            f"| {change['top1_gains']} / {change['top1_losses']} (`p={pvalue['top1']:.3f}`) "
            f"| {change['top3_gains']} / {change['top3_losses']} (`p={pvalue['top3']:.3f}`) |"
        )
    text = f"""# Network 概率融合 Resolvable Group 内部开发验证

## 设计

- 仍使用每个外层 LOSO 训练折内生成的 `resolvable Region group`，不改变分组阈值。
- 基线：Region-first group 排名，Group Top1 `{base['group_top1_accuracy']:.1%}`、Top3 `{base['group_top3_accuracy']:.1%}`。
- 融合：对候选 Region correlation 分数加入 `alpha * log(P(Network))` 后，再映射至首次出现的 group。
- `P(Network)` 由折内 Network correlation 生成，并按现有 pairwise-corrected Top3 排名对齐；测试样本标签未参与概率计算。
- 本轮比较多个固定 `alpha`，因此属于本数据集内部开发选择，不能单独作为正式替换证据。

## 结果

| alpha | Group Top1 | Group Top3 | Top1 新增 / 丢失 | Top3 新增 / 丢失 |
| ---: | ---: | ---: | ---: | ---: |
| Baseline | {base['group_top1_hits']}/{summary['n_test_samples']} ({base['group_top1_accuracy']:.1%}) | {base['group_top3_hits']}/{summary['n_test_samples']} ({base['group_top3_accuracy']:.1%}) | - | - |
{chr(10).join(rows)}

## 判定

{summary['decision']}
"""
    (outdir / "group_network_probability_fusion_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Network probability fusion on fold-local resolvable Region groups.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--network-top-n-genes", type=int, default=200)
    parser.add_argument("--alphas", default="0.02,0.05,0.10,0.20,0.40")
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.20)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=4)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    alphas = [float(value) for value in str(args.alphas).split(",") if str(value).strip()]
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
    labels = ann.set_index("sample_id").reindex(sample_ids)["endpoint_label"].to_numpy(dtype=str)
    networks = sorted(set(labels))
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

    rows: dict[str, list[dict[str, Any]]] = {BEAM_ROUTE: []}
    for alpha in alphas:
        rows[route_name(alpha)] = []
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        corrected_top3 = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        _, network_probability = network_probabilities(
            values, labels, networks, heldout_idx, heldout_idx,
            corrected_top3=corrected_top3, top_n_genes=args.network_top_n_genes
        )
        candidates = sorted(
            training_ann.loc[training_ann["endpoint_label"].isin(corrected_top3), "region_id"].unique().tolist()
        )
        local_rows = build_local_discriminative_rows(
            values, candidates, region_indices, heldout_idx, args.local_top_n_genes
        )
        if local_rows is None:
            local_rows = np.arange(values.shape[0], dtype=int)
        training = candidate_training_indices(candidates, region_indices, heldout_idx)
        annotations, _ = build_resolution_groups(
            values, candidates, training, region_network_assignment(training_ann, candidates), local_rows,
            args.min_resolution_samples, args.min_merge_samples, args.min_pair_errors, args.min_confusion_rate,
            args.similarity_threshold, args.merge_similarity_threshold, args.max_group_size
        )
        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)
        scores = correlation_scores(reference, sample, local_rows)
        indices = np.asarray([region_pos[region] for region in candidates], dtype=int)
        baseline_ranked = rank_candidates(scores, regions, indices)
        rows[BEAM_ROUTE].append(
            score_route(
                BEAM_ROUTE, sample_id, truth_region, truth_network, corrected_top3,
                baseline_ranked, annotations, len(regions)
            )
        )
        priors = region_network_prior(candidates, training_ann, network_probability)
        for alpha in alphas:
            route = route_name(alpha)
            adjusted = fused_scores(scores, regions, candidates, priors, alpha)
            ranked = rank_candidates(adjusted, regions, indices)
            rows[route].append(
                score_route(
                    route, sample_id, truth_region, truth_network, corrected_top3,
                    ranked, annotations, len(regions)
                )
            )
        fold_rows.append({
            "fold": fold,
            "sample_id": sample_id,
            "truth_region": truth_region,
            "network_top1": corrected_top3[0],
            "network_top1_probability": float(network_probability[corrected_top3[0]]),
            "truth_network_probability": float(network_probability[truth_network]),
            "n_candidate_regions": int(len(candidates)),
        })

    details = {route: pd.DataFrame(data) for route, data in rows.items()}
    metrics = {route: summarize(detail) for route, detail in details.items()}
    changes = {route: paired_changes(details[BEAM_ROUTE], details[route]) for route in rows if route != BEAM_ROUTE}
    pvalues = {
        route: {metric: paired_pvalue(change[f"{metric}_gains"], change[f"{metric}_losses"]) for metric in ["top1", "top3"]}
        for route, change in changes.items()
    }
    best_route = max([route_name(alpha) for alpha in alphas], key=lambda route: (
        metrics[route]["group_top1_accuracy"], metrics[route]["group_top3_accuracy"]
    ))
    base = metrics[BEAM_ROUTE]
    best = metrics[best_route]
    improved = (
        best["group_top1_accuracy"] > base["group_top1_accuracy"]
        and best["group_top3_accuracy"] >= base["group_top3_accuracy"]
    )
    best_alpha = alphas[[route_name(alpha) for alpha in alphas].index(best_route)]
    if improved:
        decision = (
            f"固定候选权重中，`alpha={best_alpha:.2f}` 同时提高 Group Top1 且未降低 Top3；"
            "该权重可进入下一步严格 OOF stacked ranker 对照，但由于本轮在同一数据集选择 alpha，不直接替换正式路径。"
        )
    else:
        decision = (
            "固定 Network 概率融合权重未在 Group Top1 提升且 Top3 不下降约束下超过 Region-first group 基线；"
            "不进入正式路径。"
        )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; fold-local resolution groups and fold-local network probabilities; retrospective alpha grid comparison",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singleton_regions).sum()),
        "local_top_n_genes": int(args.local_top_n_genes),
        "network_top_n_genes": int(args.network_top_n_genes),
        "alphas": alphas,
        "routes": metrics,
        "paired_changes": changes,
        "paired_pvalues": pvalues,
        "best_alpha": float(best_alpha),
        "proceed_to_stacked_ranker": bool(improved),
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **value} for route, value in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_probability_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics, alphas)
    write_report(args.outdir, summary, alphas)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
