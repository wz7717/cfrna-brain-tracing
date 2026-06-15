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
    select_fold_signature,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_correlation_marker_routes_unseen_confirmation"


def build_stable_marker_signature(
    values: np.ndarray,
    genes: list[str],
    regions: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    topk_per_region: int = 30,
    min_region_train_samples: int = 3,
    min_consistency: float = 0.75,
    min_effect: float = 0.5,
) -> pd.DataFrame:
    """Select markers using training samples only and sample-level consistency."""
    train_mask = np.ones(values.shape[1], dtype=bool)
    train_mask[heldout_idx] = False
    train_values = values[:, train_mask].astype(float, copy=False)
    global_sd = train_values.std(axis=1) + 1e-6
    total_sum = train_values.sum(axis=1)
    training_total = int(train_mask.sum())
    rows: list[dict[str, Any]] = []
    for region in regions:
        idx = region_indices[region]
        in_idx = idx[idx != heldout_idx]
        if len(in_idx) < min_region_train_samples:
            continue
        outside_n = training_total - len(in_idx)
        if outside_n <= 0:
            continue
        inside = values[:, in_idx].astype(float, copy=False)
        inside_mean = inside.mean(axis=1)
        outside_mean = (total_sum - inside.sum(axis=1)) / outside_n
        effect = (inside_mean - outside_mean) / global_sd
        consistency = (inside > outside_mean[:, None]).mean(axis=1)
        quality = effect * consistency
        eligible = np.flatnonzero((effect >= min_effect) & (consistency >= min_consistency))
        if not len(eligible):
            continue
        chosen = eligible[np.argsort(quality[eligible])[::-1][:topk_per_region]]
        for i in chosen:
            rows.append(
                {
                    "region_id": region,
                    "gene_symbol": str(genes[int(i)]),
                    "gene_index": int(i),
                    "effect_size": float(effect[i]),
                    "consistency": float(consistency[i]),
                    "marker_quality": float(quality[i]),
                    "outside_mean": float(outside_mean[i]),
                    "scale": float(global_sd[i]),
                    "n_region_train_samples": int(len(in_idx)),
                }
            )
    return pd.DataFrame(rows)


def score_stable_marker_support(
    sample: np.ndarray,
    regions: list[str],
    signature_df: pd.DataFrame,
    min_markers: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """Score whether a sample shows the fold-local marker program of each region."""
    scores = np.full(len(regions), np.nan, dtype=float)
    counts = np.zeros(len(regions), dtype=int)
    if signature_df.empty:
        return scores, counts
    region_pos = {region: j for j, region in enumerate(regions)}
    for region, group in signature_df.groupby("region_id", sort=False):
        j = region_pos.get(str(region))
        if j is None or len(group) < min_markers:
            continue
        idx = group["gene_index"].to_numpy(dtype=int)
        z = (sample[idx].astype(float) - group["outside_mean"].to_numpy(dtype=float)) / group[
            "scale"
        ].to_numpy(dtype=float)
        scores[j] = float(np.clip(z, -4.0, 4.0).mean())
        counts[j] = int(len(group))
    return scores, counts


def rerank_correlation_topk(
    correlation: np.ndarray,
    marker_support: np.ndarray,
    candidate_k: int = 10,
    marker_weight: float = 0.20,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply a conservative marker adjustment only among correlation TopK candidates."""
    corr_order = np.argsort(correlation)[::-1]
    candidates = corr_order[:candidate_k]
    valid = np.isfinite(marker_support[candidates])
    marker_z = np.zeros(len(candidates), dtype=float)
    if valid.sum() >= 2:
        raw = marker_support[candidates][valid]
        std = float(raw.std())
        if std > 1e-8:
            marker_z[valid] = (raw - raw.mean()) / std
    corr_top = correlation[candidates]
    corr_std = max(float(corr_top.std()), 1e-8)
    adjusted_top = corr_top + marker_weight * corr_std * marker_z
    local_order = np.argsort(adjusted_top)[::-1]
    reranked_order = np.concatenate([candidates[local_order], corr_order[candidate_k:]])
    adjusted = correlation.copy()
    adjusted[candidates] = adjusted_top
    return reranked_order, adjusted, marker_z


def annotate_correlation_top1(
    corr_order: np.ndarray,
    marker_support: np.ndarray,
    marker_counts: np.ndarray,
    candidate_k: int = 10,
) -> tuple[str, int | None, float | None, int]:
    candidates = corr_order[:candidate_k]
    top1 = int(candidates[0])
    valid = np.isfinite(marker_support[candidates])
    n_valid = int(valid.sum())
    if not np.isfinite(marker_support[top1]):
        return "insufficient_marker_evidence", None, None, n_valid
    valid_candidates = candidates[valid]
    order = valid_candidates[np.argsort(marker_support[valid_candidates])[::-1]]
    rank = int(np.where(order == top1)[0][0]) + 1
    if rank <= min(3, n_valid):
        status = "supported"
    elif n_valid >= 8 and rank >= n_valid - 1:
        status = "contradicted"
    else:
        status = "inconclusive"
    return status, rank, float(marker_support[top1]), int(marker_counts[top1])


def result_row(
    method: str,
    sample_id: str,
    truth: str,
    regions: list[str],
    order: np.ndarray,
    scores: np.ndarray,
) -> tuple[dict[str, Any], dict[str, float | str]]:
    predictions = [regions[int(j)] for j in order]
    true_rank = predictions.index(truth) + 1
    probs = softmax(scores)
    row = {
        "method": method,
        "sample_id": sample_id,
        "label": truth,
        "pred_top1": predictions[0],
        "pred_top2": predictions[1],
        "pred_top3": predictions[2],
        "true_rank": true_rank,
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "top1_score": float(scores[int(order[0])]),
    }
    probability: dict[str, float | str] = {"sample_id": sample_id, "label": truth}
    probability.update({region: float(probs[j]) for j, region in enumerate(regions)})
    return row, probability


def summarize_route(detail: pd.DataFrame, probabilities: pd.DataFrame, regions: list[str]) -> dict[str, float | int]:
    auc = compute_multiclass_auc(
        detail["label"].astype(str).tolist(),
        probabilities[regions],
    )
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(auc),
        "median_true_rank": float(detail["true_rank"].median()),
        "mean_true_rank": float(detail["true_rank"].mean()),
    }


def export_comparison_plot(outdir: Path, route_metrics: pd.DataFrame, annotation_counts: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["Correlation", "Top10 rerank"]
    baseline = route_metrics.loc[route_metrics["route"] == "correlation_primary"].iloc[0]
    rerank = route_metrics.loc[route_metrics["route"] == "correlation_marker_rerank_top10"].iloc[0]
    values = np.asarray(
        [
            [baseline["top1_accuracy"], baseline["top3_accuracy"], baseline["macro_auc"]],
            [rerank["top1_accuracy"], rerank["top3_accuracy"], rerank["macro_auc"]],
        ],
        dtype=float,
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7), constrained_layout=True)
    colors = ["#0072B2", "#D55E00"]
    titles = ["Top1 accuracy", "Top3 accuracy", "Macro AUC"]
    for j, title in enumerate(titles):
        bars = axes[j].bar(labels, values[:, j], color=colors, width=0.58)
        axes[j].set_title(title, fontweight="bold")
        axes[j].grid(axis="y", alpha=0.25)
        axes[j].set_axisbelow(True)
        axes[j].tick_params(axis="x", rotation=15)
        limit = max(0.5 if j < 2 else 0.9, float(values[:, j].max()) * 1.2)
        axes[j].set_ylim(0, limit)
        delta = values[1, j] - values[0, j]
        delta_text = f"Delta: {delta * 100:+.1f} pp" if j < 2 else f"Delta: {delta:+.3f}"
        axes[j].text(0.5, 0.97, delta_text, transform=axes[j].transAxes, ha="center", va="top")
        for bar, value in zip(bars, values[:, j]):
            text = f"{value:.1%}" if j < 2 else f"{value:.3f}"
            axes[j].text(
                bar.get_x() + bar.get_width() / 2,
                value + limit * 0.025,
                text,
                ha="center",
                va="bottom",
                fontweight="bold",
            )
    annotation_text = ", ".join(
        f"{r.annotation_status}: {int(r.n)} (Top1 {r.top1_accuracy:.1%})"
        for r in annotation_counts.itertuples(index=False)
    )
    fig.suptitle("Bo2023 unseen cohort: conservative marker routes (n=30)", fontweight="bold", fontsize=14)
    fig.text(0.5, -0.01, f"Annotation keeps correlation Top1 unchanged. {annotation_text}", ha="center", fontsize=9)
    fig.savefig(outdir / "correlation_marker_routes_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "correlation_marker_routes_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    status_order = ["supported", "inconclusive", "contradicted"]
    plot_counts = annotation_counts.set_index("annotation_status").reindex(status_order).fillna(0).reset_index()
    fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    bars = ax.bar(
        plot_counts["annotation_status"],
        plot_counts["top1_accuracy"],
        color=["#009E73", "#E69F00", "#CC3311"],
        width=0.6,
    )
    ax.axhline(
        float(route_metrics.loc[route_metrics["route"] == "correlation_primary", "top1_accuracy"].iloc[0]),
        color="#0072B2",
        linestyle="--",
        label="Overall correlation Top1",
    )
    ax.set_ylim(0, 0.55)
    ax.set_ylabel("Correlation Top1 accuracy")
    ax.set_title("Marker annotation as confidence stratification", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    for bar, row in zip(bars, plot_counts.itertuples(index=False)):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            float(row.top1_accuracy) + 0.025,
            f"{float(row.top1_accuracy):.1%}\n(n={int(row.n)})",
            ha="center",
            va="bottom",
        )
    ax.legend()
    fig.savefig(outdir / "marker_annotation_confidence_stratification.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "marker_annotation_confidence_stratification.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(
    outdir: Path,
    metrics: pd.DataFrame,
    annotation_counts: pd.DataFrame,
    paired: dict[str, Any],
    summary: dict[str, Any],
) -> None:
    baseline = metrics.loc[metrics["route"] == "correlation_primary"].iloc[0]
    rerank = metrics.loc[metrics["route"] == "correlation_marker_rerank_top10"].iloc[0]
    annotation_lines = "\n".join(
        f"| `{r.annotation_status}` | {int(r.n)} | {int(r.top1_hits)}/{int(r.n)} ({r.top1_accuracy:.1%}) |"
        for r in annotation_counts.itertuples(index=False)
    )
    text = f"""# Correlation 结合稳定 marker 的两条保守路线测试

## 设计

- 数据：Bo2023 VSD 图谱；新的 30 个未见样本，排除此前使用过的 `{summary['n_prior_samples_excluded']}` 个测试样本。
- 对照：正式 V2 当前 `correlation` 主路径，在每折训练样本生成的 signature 基因集合上计算相关性。
- 稳定 marker：只在训练折内构建；要求 Region 至少 `{summary['min_region_train_samples']}` 个训练样本、方向一致性不低于 `{summary['min_consistency']:.0%}`、标准化效应量不低于 `{summary['min_effect']:.2f}`，每 Region 最多 `{summary['stable_topk_per_region']}` 个。
- 路线 A：`Correlation + marker rerank Top10`，只在 correlation Top10 内重排；marker 调整权重固定为 `{summary['rerank_marker_weight']:.2f}`。
- 路线 B：`Correlation + marker confidence annotation`，不修改 correlation 排名；依据 Top10 候选内 marker 支持排名标记 `supported`、`contradicted` 或 `inconclusive`。

## 路线 A：Top10 重排

| 方法 | Top1 | Top3 | Macro AUC | 真实 Region 中位排名 |
| --- | ---: | ---: | ---: | ---: |
| Correlation 主路径 | {int(baseline.top1_hits)}/30 ({baseline.top1_accuracy:.1%}) | {int(baseline.top3_hits)}/30 ({baseline.top3_accuracy:.1%}) | {baseline.macro_auc:.3f} | {baseline.median_true_rank:.1f} |
| Correlation + marker rerank Top10 | {int(rerank.top1_hits)}/30 ({rerank.top1_accuracy:.1%}) | {int(rerank.top3_hits)}/30 ({rerank.top3_accuracy:.1%}) | {rerank.macro_auc:.3f} | {rerank.median_true_rank:.1f} |
| 差值 | {rerank.top1_accuracy - baseline.top1_accuracy:+.1%} | {rerank.top3_accuracy - baseline.top3_accuracy:+.1%} | {rerank.macro_auc - baseline.macro_auc:+.3f} | {rerank.median_true_rank - baseline.median_true_rank:+.1f} |

Top1 配对变化：新增命中 `{paired['rerank_top1_gains']}` 个，丢失命中 `{paired['rerank_top1_losses']}` 个。Top3 配对变化：新增命中 `{paired['rerank_top3_gains']}` 个，丢失命中 `{paired['rerank_top3_losses']}` 个。

## 路线 B：Marker 支持注释

该路线保持 correlation 的 Top1/Top3 和 AUC 不变，仅评估 marker 注释是否能够区分可信预测。

| Marker 注释状态 | 样本数 | Correlation Top1 命中 |
| --- | ---: | ---: |
{annotation_lines}

## 判定

`{summary['decision']}`

`supported` 覆盖了大多数样本；即使该组准确率高于总体，当前分层仍过宽，不能据此直接拒绝或采纳某个预测。

## 下一步

`{summary['next_step']}`
"""
    (outdir / "correlation_marker_routes_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate conservative correlation plus stable-marker routes on Bo2023.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260528)
    parser.add_argument("--reference-topk-per-region", type=int, default=80)
    parser.add_argument("--stable-topk-per-region", type=int, default=30)
    parser.add_argument("--min-region-train-samples", type=int, default=3)
    parser.add_argument("--min-consistency", type=float, default=0.75)
    parser.add_argument("--min-effect", type=float, default=0.5)
    parser.add_argument("--min-markers-for-support", type=int, default=8)
    parser.add_argument("--rerank-candidate-k", type=int, default=10)
    parser.add_argument("--rerank-marker-weight", type=float, default=0.20)
    parser.add_argument("--exclude-samples", type=Path, action="append", default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "fold_stable_markers").mkdir(parents=True, exist_ok=True)
    raw_matrix = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw_matrix, args.gene_map)
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
    genes = matrix.index.astype(str).tolist()
    full_reference, regions, region_counts, region_indices = build_region_reference(
        values, matrix.columns.tolist(), ann
    )
    sample_pos = {sample_id: j for j, sample_id in enumerate(matrix.columns)}
    region_pos = {region: j for j, region in enumerate(regions)}
    detail_by_route: dict[str, list[dict[str, Any]]] = {
        "correlation_primary": [],
        "correlation_marker_rerank_top10": [],
    }
    probability_by_route: dict[str, list[dict[str, float | str]]] = {
        "correlation_primary": [],
        "correlation_marker_rerank_top10": [],
    }
    annotation_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    marker_manifest: list[dict[str, Any]] = []

    for fold_no, heldout in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(heldout.sample_id)
        truth = str(heldout.region_id)
        test_idx = sample_pos[sample_id]
        truth_j = region_pos[truth]
        reference = full_reference.copy()
        truth_train_idx = region_indices[truth][region_indices[truth] != test_idx]
        reference[:, truth_j] = values[:, truth_train_idx].mean(axis=1, dtype=np.float64).astype(np.float32)
        signature_mask, _ = select_fold_signature(
            genes, regions, reference, args.reference_topk_per_region
        )
        correlation = correlation_scores(reference, values[:, test_idx], signature_mask)
        corr_order = np.argsort(correlation)[::-1]
        stable_df = build_stable_marker_signature(
            values,
            genes,
            regions,
            region_indices,
            heldout_idx=test_idx,
            topk_per_region=args.stable_topk_per_region,
            min_region_train_samples=args.min_region_train_samples,
            min_consistency=args.min_consistency,
            min_effect=args.min_effect,
        )
        marker_scores, marker_counts = score_stable_marker_support(
            values[:, test_idx], regions, stable_df, min_markers=args.min_markers_for_support
        )
        rerank_order, rerank_scores, marker_z = rerank_correlation_topk(
            correlation,
            marker_scores,
            candidate_k=args.rerank_candidate_k,
            marker_weight=args.rerank_marker_weight,
        )
        baseline_row, baseline_prob = result_row(
            "correlation_primary", sample_id, truth, regions, corr_order, correlation
        )
        rerank_row, rerank_prob = result_row(
            "correlation_marker_rerank_top10", sample_id, truth, regions, rerank_order, rerank_scores
        )
        detail_by_route["correlation_primary"].append(baseline_row)
        detail_by_route["correlation_marker_rerank_top10"].append(rerank_row)
        probability_by_route["correlation_primary"].append(baseline_prob)
        probability_by_route["correlation_marker_rerank_top10"].append(rerank_prob)
        annotation_status, marker_rank, top1_marker_score, top1_marker_count = annotate_correlation_top1(
            corr_order, marker_scores, marker_counts, candidate_k=args.rerank_candidate_k
        )
        annotation_rows.append(
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "truth_region": truth,
                "correlation_top1": regions[int(corr_order[0])],
                "correlation_hit1": baseline_row["hit1"],
                "correlation_hit3": baseline_row["hit3"],
                "annotation_status": annotation_status,
                "top1_marker_rank_within_candidates": marker_rank,
                "top1_marker_support_score": top1_marker_score,
                "top1_marker_count": top1_marker_count,
            }
        )
        candidates = corr_order[: args.rerank_candidate_k]
        for local_rank, j in enumerate(candidates, start=1):
            candidate_rows.append(
                {
                    "fold": fold_no,
                    "sample_id": sample_id,
                    "truth_region": truth,
                    "candidate_region": regions[int(j)],
                    "correlation_rank": local_rank,
                    "correlation_score": float(correlation[int(j)]),
                    "marker_support_score": float(marker_scores[int(j)]) if np.isfinite(marker_scores[int(j)]) else np.nan,
                    "marker_count": int(marker_counts[int(j)]),
                    "marker_z_for_rerank": float(marker_z[local_rank - 1]),
                    "is_truth": int(regions[int(j)] == truth),
                }
            )
        marker_file = args.outdir / "fold_stable_markers" / f"fold_{fold_no:02d}_{sample_id}_stable_markers.csv"
        stable_df.to_csv(marker_file, index=False, encoding="utf-8-sig")
        marker_manifest.append(
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "truth_region": truth,
                "n_marker_pairs": int(len(stable_df)),
                "n_regions_with_marker_support": int(np.isfinite(marker_scores).sum()),
                "file": str(marker_file.relative_to(args.outdir)),
            }
        )

    detail_frames = {key: pd.DataFrame(rows) for key, rows in detail_by_route.items()}
    probability_frames = {key: pd.DataFrame(rows) for key, rows in probability_by_route.items()}
    summaries = {
        key: summarize_route(detail_frames[key], probability_frames[key], regions)
        for key in detail_frames
    }
    route_metrics = pd.DataFrame(
        [{"route": key, **value} for key, value in summaries.items()]
    )
    annotations = pd.DataFrame(annotation_rows)
    annotation_counts = (
        annotations.groupby("annotation_status", dropna=False)["correlation_hit1"]
        .agg(n="size", top1_hits="sum", top1_accuracy="mean")
        .reset_index()
    )
    baseline = detail_frames["correlation_primary"]
    rerank = detail_frames["correlation_marker_rerank_top10"]
    paired = {
        "rerank_top1_gains": int(((baseline["hit1"] == 0) & (rerank["hit1"] == 1)).sum()),
        "rerank_top1_losses": int(((baseline["hit1"] == 1) & (rerank["hit1"] == 0)).sum()),
        "rerank_top3_gains": int(((baseline["hit3"] == 0) & (rerank["hit3"] == 1)).sum()),
        "rerank_top3_losses": int(((baseline["hit3"] == 1) & (rerank["hit3"] == 0)).sum()),
    }
    baseline_metrics = summaries["correlation_primary"]
    rerank_metrics = summaries["correlation_marker_rerank_top10"]
    improves_primary = (
        rerank_metrics["top1_accuracy"] > baseline_metrics["top1_accuracy"]
        and rerank_metrics["top3_accuracy"] >= baseline_metrics["top3_accuracy"]
    )
    supported_row = annotation_counts[annotation_counts["annotation_status"] == "supported"]
    supported_enriches = (
        len(supported_row) > 0
        and float(supported_row.iloc[0]["top1_accuracy"]) > float(baseline_metrics["top1_accuracy"])
    )
    mean_marker_pairs = float(np.mean([row["n_marker_pairs"] for row in marker_manifest]))
    mean_supported_regions = float(
        np.mean([row["n_regions_with_marker_support"] for row in marker_manifest])
    )
    if improves_primary:
        decision = "Top10 重排在本轮未见队列上改善主指标，可进入另一批未见样本复核，但在复核通过前仍不替换 correlation 默认路径。"
    else:
        decision = "Top10 重排未满足超过 correlation 的条件，不切换主路径。"
    if supported_enriches:
        next_step = "保留 marker 作为置信度注释候选，继续在另一批未见样本验证 supported 分组是否持续富集正确预测；Top10 重排仅在重新校准 marker 权重后再测。"
    else:
        next_step = "两条 marker 路径均不升级；先重构稳定 marker 的邻近脑区区分能力和权重校准，再开启新的未见队列确认。"
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict LOSO on previously unseen test cohort; fold-local stable markers",
        "seed": int(args.seed),
        "n_test_samples": int(args.n_samples),
        "n_prior_samples_excluded": int(len(excluded_ids)),
        "n_regions": int(len(regions)),
        "n_gene_symbols": int(matrix.shape[0]),
        "n_singleton_samples_excluded": int(len(singleton_samples)),
        "stable_topk_per_region": int(args.stable_topk_per_region),
        "min_region_train_samples": int(args.min_region_train_samples),
        "min_consistency": float(args.min_consistency),
        "min_effect": float(args.min_effect),
        "min_markers_for_support": int(args.min_markers_for_support),
        "rerank_candidate_k": int(args.rerank_candidate_k),
        "rerank_marker_weight": float(args.rerank_marker_weight),
        "mean_fold_marker_pairs": mean_marker_pairs,
        "mean_regions_with_marker_support": mean_supported_regions,
        "routes": summaries,
        "paired_changes": paired,
        "decision": decision,
        "next_step": next_step,
    }
    route_metrics.to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    annotations.to_csv(args.outdir / "marker_confidence_annotations.csv", index=False, encoding="utf-8-sig")
    annotation_counts.to_csv(args.outdir / "annotation_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(candidate_rows).to_csv(args.outdir / "top10_candidate_marker_scores.csv", index=False, encoding="utf-8-sig")
    for key, frame in detail_frames.items():
        frame.to_csv(args.outdir / f"{key}_detail.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(args.outdir / "selected_test_samples.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "marker_fold_manifest.json").write_text(
        json.dumps(marker_manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.outdir / "validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    export_comparison_plot(args.outdir, route_metrics, annotation_counts)
    write_report(args.outdir, route_metrics, annotation_counts, paired, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
