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
from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    correlation_scores,
    read_annotations,
    read_vsd_matrix,
    softmax,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    choose_validation_samples,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_discriminative_correlation_unseen_confirmation"


def select_fold_discriminative_genes(
    values: np.ndarray,
    regions: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    top_n: int = 1000,
    min_region_train_samples: int = 3,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Select high between/within-region variance genes from training samples only."""
    means: list[np.ndarray] = []
    weights: list[int] = []
    within_sum = np.zeros(values.shape[0], dtype=float)
    within_denom = 0
    for region in regions:
        idx = region_indices[region]
        train_idx = idx[idx != heldout_idx]
        if len(train_idx) < min_region_train_samples:
            continue
        x = values[:, train_idx].astype(float, copy=False)
        mean = x.mean(axis=1)
        means.append(mean)
        weights.append(int(len(train_idx)))
        within_sum += np.square(x - mean[:, None]).sum(axis=1)
        within_denom += len(train_idx) - 1
    if not means or within_denom <= 0:
        raise ValueError("insufficient training regions to select discriminative genes")
    centroid_matrix = np.column_stack(means)
    region_weights = np.asarray(weights, dtype=float)
    overall = np.average(centroid_matrix, axis=1, weights=region_weights)
    between = np.average(np.square(centroid_matrix - overall[:, None]), axis=1, weights=region_weights)
    within = within_sum / float(within_denom)
    fisher_score = between / (within + 1e-8)
    selected = np.argsort(fisher_score)[::-1][: min(int(top_n), len(fisher_score))]
    audit = pd.DataFrame(
        {
            "gene_index": selected.astype(int),
            "fisher_score": fisher_score[selected],
            "between_variance": between[selected],
            "within_variance": within[selected],
        }
    )
    return selected.astype(int), audit


def score_route(
    sample_id: str,
    truth: str,
    sample: np.ndarray,
    reference: np.ndarray,
    regions: list[str],
    rows: np.ndarray | None,
    route: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    scores = correlation_scores(reference, sample, rows)
    order = np.argsort(scores)[::-1]
    ranked = [regions[int(j)] for j in order]
    true_rank = ranked.index(truth) + 1
    probs = softmax(scores)
    detail = {
        "route": route,
        "sample_id": sample_id,
        "label": truth,
        "pred_top1": ranked[0],
        "pred_top2": ranked[1],
        "pred_top3": ranked[2],
        "true_rank": true_rank,
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "top1_score": float(scores[int(order[0])]),
        "true_region_score": float(scores[regions.index(truth)]),
        "decision_margin": float(scores[int(order[0])] - scores[int(order[1])]),
    }
    probability: dict[str, Any] = {"sample_id": sample_id, "label": truth}
    probability.update({region: float(probs[j]) for j, region in enumerate(regions)})
    return detail, probability


def summarize(detail: pd.DataFrame, probability: pd.DataFrame, regions: list[str]) -> dict[str, float | int]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probability[regions])),
        "median_true_rank": float(detail["true_rank"].median()),
        "mean_true_rank": float(detail["true_rank"].mean()),
    }


def export_plot(outdir: Path, metrics: pd.DataFrame, target_top1: float, target_top3: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base = metrics.loc[metrics["route"] == "correlation_all_genes"].iloc[0]
    enhanced = metrics.loc[metrics["route"] == "discriminative_correlation_top1000"].iloc[0]
    labels = ["Current correlation", "Discriminative correlation"]
    values = np.asarray(
        [
            [base.top1_accuracy, base.top3_accuracy, base.macro_auc],
            [enhanced.top1_accuracy, enhanced.top3_accuracy, enhanced.macro_auc],
        ]
    )
    targets = [target_top1, target_top3, None]
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.7), constrained_layout=True)
    for j, title in enumerate(["Top1 accuracy", "Top3 accuracy", "Macro AUC"]):
        bars = axes[j].bar(labels, values[:, j], color=["#0072B2", "#009E73"], width=0.58)
        axes[j].set_title(title, fontweight="bold")
        axes[j].grid(axis="y", alpha=0.25)
        axes[j].set_axisbelow(True)
        axes[j].tick_params(axis="x", rotation=15)
        if targets[j] is not None:
            axes[j].axhline(targets[j], color="#CC3311", linestyle="--", label="Requested target")
            axes[j].legend(fontsize=8)
        axes[j].set_ylim(0, max(0.85 if j < 2 else 0.95, float(values[:, j].max()) + 0.12))
        for bar, value in zip(bars, values[:, j]):
            text = f"{value:.1%}" if j < 2 else f"{value:.3f}"
            axes[j].text(bar.get_x() + bar.get_width() / 2, value + 0.025, text, ha="center", fontweight="bold")
    fig.suptitle("Bo2023 unseen confirmation: discriminative correlation path", fontsize=14, fontweight="bold")
    fig.savefig(outdir / "discriminative_correlation_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "discriminative_correlation_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    base = summary["routes"]["correlation_all_genes"]
    new = summary["routes"]["discriminative_correlation_top1000"]
    text = f"""# Correlation 主路径判别基因优化验证

## 设计

- 测试队列：新的 `{summary['n_test_samples']}` 个未见样本，排除此前已检查的 `{summary['n_prior_samples_excluded']}` 个样本。
- 当前对照：全部 gene 的 VSD centroid Pearson correlation。
- 优化路径：每个 LOSO fold 仅用训练样本计算 Region 间方差 / Region 内方差比，选择 Top `{summary['top_n_genes']}` 个判别基因，再计算 Pearson correlation。
- 低样本 Region：少于 `{summary['min_region_train_samples']}` 个训练样本的 Region 不参与判别基因筛选，但仍保留为可预测候选。
- 本轮未使用 marker 或测试集调参。

## 结果

| 方法 | Top1 | Top3 | Macro AUC | 真实 Region 中位排名 |
| --- | ---: | ---: | ---: | ---: |
| 当前 correlation | {base['top1_hits']}/30 ({base['top1_accuracy']:.1%}) | {base['top3_hits']}/30 ({base['top3_accuracy']:.1%}) | {base['macro_auc']:.3f} | {base['median_true_rank']:.1f} |
| 判别基因 correlation | {new['top1_hits']}/30 ({new['top1_accuracy']:.1%}) | {new['top3_hits']}/30 ({new['top3_accuracy']:.1%}) | {new['macro_auc']:.3f} | {new['median_true_rank']:.1f} |
| 目标 | 50.0% | 70.0% | - | - |

## 判定

{summary['decision']}

## 下一步

{summary['next_step']}
"""
    (outdir / "discriminative_correlation_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict unseen validation of discriminative-gene Pearson correlation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260530)
    parser.add_argument("--top-n-genes", type=int, default=1000)
    parser.add_argument("--min-region-train-samples", type=int, default=3)
    parser.add_argument("--target-top1", type=float, default=0.50)
    parser.add_argument("--target-top3", type=float, default=0.70)
    parser.add_argument("--exclude-samples", type=Path, action="append", default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "fold_selected_genes").mkdir(parents=True, exist_ok=True)
    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    excluded_ids: set[str] = set()
    for path in args.exclude_samples:
        if path.exists():
            excluded_ids.update(pd.read_csv(path)["sample_id"].astype(str))
    selected, singleton_samples = choose_validation_samples(
        ann, matrix.columns.tolist(), args.n_samples, args.seed, excluded_ids
    )
    values = matrix.to_numpy(dtype=np.float32)
    genes = matrix.index.astype(str).to_numpy()
    full_reference, regions, region_counts, region_indices = build_region_reference(
        values, matrix.columns.tolist(), ann
    )
    sample_pos = {sample_id: j for j, sample_id in enumerate(matrix.columns)}
    region_pos = {region: j for j, region in enumerate(regions)}
    detail: dict[str, list[dict[str, Any]]] = {
        "correlation_all_genes": [],
        "discriminative_correlation_top1000": [],
    }
    probabilities: dict[str, list[dict[str, Any]]] = {key: [] for key in detail}
    fold_manifest: list[dict[str, Any]] = []

    for fold_no, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth = str(row.region_id)
        heldout_idx = sample_pos[sample_id]
        truth_j = region_pos[truth]
        reference = full_reference.copy()
        truth_train = region_indices[truth][region_indices[truth] != heldout_idx]
        reference[:, truth_j] = values[:, truth_train].mean(axis=1, dtype=np.float64).astype(np.float32)
        selected_rows, audit = select_fold_discriminative_genes(
            values,
            regions,
            region_indices,
            heldout_idx,
            top_n=args.top_n_genes,
            min_region_train_samples=args.min_region_train_samples,
        )
        audit["gene_symbol"] = genes[audit["gene_index"].to_numpy(dtype=int)]
        gene_file = args.outdir / "fold_selected_genes" / f"fold_{fold_no:02d}_{sample_id}_top_genes.csv"
        audit.to_csv(gene_file, index=False, encoding="utf-8-sig")
        for route, rows in [
            ("correlation_all_genes", None),
            ("discriminative_correlation_top1000", selected_rows),
        ]:
            output, probability = score_route(
                sample_id, truth, values[:, heldout_idx], reference, regions, rows, route
            )
            detail[route].append(output)
            probabilities[route].append(probability)
        fold_manifest.append(
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "truth_region": truth,
                "n_selected_genes": int(len(selected_rows)),
                "selected_gene_file": str(gene_file.relative_to(args.outdir)),
            }
        )
    detail_frames = {key: pd.DataFrame(value) for key, value in detail.items()}
    probability_frames = {key: pd.DataFrame(value) for key, value in probabilities.items()}
    route_results = {key: summarize(detail_frames[key], probability_frames[key], regions) for key in detail}
    base = detail_frames["correlation_all_genes"]
    enhanced = detail_frames["discriminative_correlation_top1000"]
    paired = {
        "top1_gains": int(((base["hit1"] == 0) & (enhanced["hit1"] == 1)).sum()),
        "top1_losses": int(((base["hit1"] == 1) & (enhanced["hit1"] == 0)).sum()),
        "top3_gains": int(((base["hit3"] == 0) & (enhanced["hit3"] == 1)).sum()),
        "top3_losses": int(((base["hit3"] == 1) & (enhanced["hit3"] == 0)).sum()),
    }
    achieved = (
        route_results["discriminative_correlation_top1000"]["top1_accuracy"] >= args.target_top1
        and route_results["discriminative_correlation_top1000"]["top3_accuracy"] >= args.target_top3
    )
    if achieved:
        decision = "判别基因 correlation 在本轮未见样本上达到指定 Top1/Top3 门槛，可进入更大规模严格复核；在复核完成前不直接覆盖生产默认。"
        next_step = "运行扩大的未见样本确认，并将折内选择出的稳定判别基因面板固化为 Bo2023 correlation reference 配置。"
    else:
        decision = "判别基因 correlation 未达到 Top1 50% / Top3 70% 门槛，因此不能宣称该目标已实现，也不替换正式默认路径。"
        next_step = "当前 110 个细粒度 Region 与每区样本量不足限制明显。应评估预注册的分层终点：先预测可重复的大区/网络，再在高可分辨 Region 子集内细分；若必须保持 110 类精确终点，需要新增训练样本或外部图谱，而非继续在同一批样本上调参。"
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict unseen LOSO with fold-local Fisher gene selection",
        "n_test_samples": int(len(selected)),
        "n_prior_samples_excluded": int(len(excluded_ids)),
        "seed": int(args.seed),
        "n_regions": int(len(regions)),
        "n_gene_symbols": int(matrix.shape[0]),
        "n_singleton_samples_excluded": int(len(singleton_samples)),
        "top_n_genes": int(args.top_n_genes),
        "min_region_train_samples": int(args.min_region_train_samples),
        "requested_target": {"top1": float(args.target_top1), "top3": float(args.target_top3)},
        "routes": route_results,
        "paired_changes": paired,
        "target_achieved": bool(achieved),
        "decision": decision,
        "next_step": next_step,
    }
    route_metrics = pd.DataFrame([{"route": key, **value} for key, value in route_results.items()])
    route_metrics.to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(args.outdir / "selected_test_samples.csv", index=False, encoding="utf-8-sig")
    for key, frame in detail_frames.items():
        frame.to_csv(args.outdir / f"{key}_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "fold_manifest.json").write_text(json.dumps(fold_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, route_metrics, args.target_top1, args.target_top3)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
