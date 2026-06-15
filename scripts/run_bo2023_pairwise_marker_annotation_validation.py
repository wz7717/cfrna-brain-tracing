#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
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


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_pairwise_marker_annotation_unseen_confirmation"


def _stable_hash(value: str, seed: int) -> int:
    payload = f"{seed}:{value}".encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16)


def split_inner_training_indices(
    region_indices: dict[str, np.ndarray],
    sample_ids: list[str],
    heldout_idx: int,
    seed: int,
    calibration_fraction: float = 0.20,
    min_build_samples: int = 3,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    build_indices: dict[str, np.ndarray] = {}
    calibration: list[int] = []
    for region, original in region_indices.items():
        available = np.asarray([i for i in original if int(i) != heldout_idx], dtype=int)
        if len(available) >= min_build_samples + 1:
            ordered = sorted(available.tolist(), key=lambda i: _stable_hash(sample_ids[int(i)], seed))
            n_calibration = max(1, int(round(len(available) * calibration_fraction)))
            n_calibration = min(n_calibration, len(available) - min_build_samples)
            calibration.extend(ordered[:n_calibration])
            build_indices[region] = np.asarray(ordered[n_calibration:], dtype=int)
        else:
            build_indices[region] = available
    return build_indices, np.asarray(calibration, dtype=int)


def build_reference_from_indices(
    values: np.ndarray,
    regions: list[str],
    indices: dict[str, np.ndarray],
) -> np.ndarray:
    reference = np.zeros((values.shape[0], len(regions)), dtype=np.float32)
    for j, region in enumerate(regions):
        idx = indices[region]
        if not len(idx):
            raise ValueError(f"region {region} has no inner-build samples")
        reference[:, j] = values[:, idx].mean(axis=1, dtype=np.float64).astype(np.float32)
    return reference


def derive_confusion_neighbors(
    reference: np.ndarray,
    regions: list[str],
    build_indices: dict[str, np.ndarray],
    neighbors_per_region: int = 5,
    min_samples_per_region: int = 3,
) -> dict[str, list[str]]:
    similarities = np.corrcoef(reference.T)
    similarities = np.nan_to_num(similarities, nan=-np.inf)
    out: dict[str, list[str]] = {}
    for j, region in enumerate(regions):
        if len(build_indices[region]) < min_samples_per_region:
            out[region] = []
            continue
        ranked = np.argsort(similarities[j])[::-1]
        candidates = [
            regions[int(k)]
            for k in ranked
            if int(k) != j and len(build_indices[regions[int(k)]]) >= min_samples_per_region
        ]
        out[region] = candidates[:neighbors_per_region]
    return out


def build_pairwise_marker_signature(
    values: np.ndarray,
    genes: list[str],
    build_indices: dict[str, np.ndarray],
    neighbor_map: dict[str, list[str]],
    max_markers_per_pair: int = 6,
    max_markers_per_region: int = 18,
    min_effect: float = 0.80,
    min_consistency: float = 0.75,
) -> pd.DataFrame:
    """Build sparse directional markers for correlation-neighbor region pairs."""
    rows: list[dict[str, Any]] = []
    for region, competitors in neighbor_map.items():
        used_genes: set[int] = set()
        inside = values[:, build_indices[region]].astype(float, copy=False)
        if inside.shape[1] < 3:
            continue
        for competitor in competitors:
            outside = values[:, build_indices[competitor]].astype(float, copy=False)
            if outside.shape[1] < 3 or len(used_genes) >= max_markers_per_region:
                continue
            inside_mean = inside.mean(axis=1)
            outside_mean = outside.mean(axis=1)
            pooled_sd = np.sqrt((inside.var(axis=1) + outside.var(axis=1)) / 2.0) + 1e-6
            effect = (inside_mean - outside_mean) / pooled_sd
            midpoint = (inside_mean + outside_mean) / 2.0
            consistency_inside = (inside > midpoint[:, None]).mean(axis=1)
            consistency_outside = (outside < midpoint[:, None]).mean(axis=1)
            consistency = np.minimum(consistency_inside, consistency_outside)
            quality = effect * consistency
            eligible = np.flatnonzero((effect >= min_effect) & (consistency >= min_consistency))
            picked = 0
            for i in eligible[np.argsort(quality[eligible])[::-1]]:
                gene_idx = int(i)
                if gene_idx in used_genes:
                    continue
                rows.append(
                    {
                        "region_id": region,
                        "competitor_region": competitor,
                        "gene_symbol": str(genes[gene_idx]),
                        "gene_index": gene_idx,
                        "effect_size": float(effect[gene_idx]),
                        "consistency": float(consistency[gene_idx]),
                        "marker_quality": float(quality[gene_idx]),
                        "midpoint": float(midpoint[gene_idx]),
                        "scale": float(pooled_sd[gene_idx]),
                    }
                )
                used_genes.add(gene_idx)
                picked += 1
                if picked >= max_markers_per_pair or len(used_genes) >= max_markers_per_region:
                    break
    return pd.DataFrame(rows)


def pairwise_annotation_score(
    sample: np.ndarray,
    predicted_region: str,
    candidate_regions: list[str],
    marker_df: pd.DataFrame,
    min_markers_per_pair: int = 3,
) -> tuple[float, int, int]:
    if marker_df.empty:
        return float("nan"), 0, 0
    pred_markers = marker_df[marker_df["region_id"] == predicted_region]
    pair_scores: list[float] = []
    n_markers = 0
    for competitor in candidate_regions:
        if competitor == predicted_region:
            continue
        pair = pred_markers[pred_markers["competitor_region"] == competitor]
        if len(pair) < min_markers_per_pair:
            continue
        idx = pair["gene_index"].to_numpy(dtype=int)
        z = (sample[idx].astype(float) - pair["midpoint"].to_numpy(dtype=float)) / pair[
            "scale"
        ].to_numpy(dtype=float)
        pair_scores.append(float(np.clip(z, -4.0, 4.0).mean()))
        n_markers += int(len(pair))
    if not pair_scores:
        return float("nan"), 0, 0
    # The weakest local competitor evidence controls support for a Top1 call.
    return float(min(pair_scores)), len(pair_scores), n_markers


def calibrate_annotation_thresholds(
    calibration: pd.DataFrame,
    min_group_size: int = 8,
) -> dict[str, float | int]:
    valid = calibration[np.isfinite(calibration["annotation_score"])].copy()
    if len(valid) < min_group_size * 2:
        return {
            "support_threshold": float("inf"),
            "contradiction_threshold": float("-inf"),
            "n_valid_calibration": int(len(valid)),
            "calibration_accuracy": float("nan"),
        }
    overall = float(valid["hit1"].mean())
    min_size = min(min_group_size, max(3, len(valid) // 4))
    support_candidates: list[tuple[float, float, int]] = []
    for threshold in np.unique(valid["annotation_score"]):
        group = valid[valid["annotation_score"] >= threshold]
        if len(group) >= min_size:
            support_candidates.append((float(group["hit1"].mean()), float(threshold), int(len(group))))
    if support_candidates:
        support_candidates.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
        support_precision, support_threshold, _ = support_candidates[0]
        if support_precision <= overall:
            support_threshold = float("inf")
    else:
        support_threshold = float("inf")
    contradict_candidates: list[tuple[float, float, int]] = []
    for threshold in np.unique(valid["annotation_score"]):
        group = valid[valid["annotation_score"] <= threshold]
        if len(group) >= min_size:
            contradict_candidates.append((float((1 - group["hit1"]).mean()), float(threshold), int(len(group))))
    if contradict_candidates:
        contradict_candidates.sort(key=lambda item: (item[0], -item[1], item[2]), reverse=True)
        error_precision, contradiction_threshold, _ = contradict_candidates[0]
        if error_precision <= 1.0 - overall:
            contradiction_threshold = float("-inf")
    else:
        contradiction_threshold = float("-inf")
    return {
        "support_threshold": float(support_threshold),
        "contradiction_threshold": float(contradiction_threshold),
        "n_valid_calibration": int(len(valid)),
        "calibration_accuracy": overall,
    }


def apply_annotation(score: float, thresholds: dict[str, float | int]) -> str:
    if not np.isfinite(score):
        return "insufficient_marker_evidence"
    if score >= float(thresholds["support_threshold"]):
        return "supported"
    if score <= float(thresholds["contradiction_threshold"]):
        return "contradicted"
    return "inconclusive"


def evaluate_correlation_prediction(
    sample_id: str,
    truth: str,
    sample: np.ndarray,
    reference: np.ndarray,
    regions: list[str],
    signature_mask: np.ndarray,
) -> tuple[dict[str, Any], dict[str, Any], np.ndarray]:
    scores = correlation_scores(reference, sample, signature_mask)
    order = np.argsort(scores)[::-1]
    ranked_regions = [regions[int(j)] for j in order]
    true_rank = ranked_regions.index(truth) + 1
    probs = softmax(scores)
    detail = {
        "sample_id": sample_id,
        "label": truth,
        "pred_top1": ranked_regions[0],
        "pred_top2": ranked_regions[1],
        "pred_top3": ranked_regions[2],
        "true_rank": true_rank,
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
    }
    probability: dict[str, Any] = {"sample_id": sample_id, "label": truth}
    probability.update({region: float(probs[j]) for j, region in enumerate(regions)})
    return detail, probability, order


def export_annotation_plot(outdir: Path, overall_accuracy: float, annotation_summary: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    desired = ["supported", "inconclusive", "contradicted", "insufficient_marker_evidence"]
    table = annotation_summary.set_index("annotation_status").reindex(desired).fillna(0).reset_index()
    labels = ["supported", "inconclusive", "contradicted", "insufficient"]
    fig, ax = plt.subplots(figsize=(8.7, 5.0), constrained_layout=True)
    bars = ax.bar(labels, table["top1_accuracy"], color=["#009E73", "#E69F00", "#CC3311", "#999999"])
    ax.axhline(overall_accuracy, linestyle="--", color="#0072B2", linewidth=2, label="Overall correlation Top1")
    ax.set_ylim(0, max(0.55, float(table["top1_accuracy"].max()) + 0.15))
    ax.set_ylabel("Correlation Top1 accuracy")
    ax.set_title("Pairwise marker annotation on unseen samples", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    for bar, row in zip(bars, table.itertuples(index=False)):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            float(row.top1_accuracy) + 0.025,
            f"{float(row.top1_accuracy):.1%}\n(n={int(row.n)})",
            ha="center",
            va="bottom",
        )
    ax.legend()
    fig.savefig(outdir / "pairwise_marker_annotation_stratification.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "pairwise_marker_annotation_stratification.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any], groups: pd.DataFrame) -> None:
    group_lines = "\n".join(
        f"| `{r.annotation_status}` | {int(r.n)} | {int(r.top1_hits)}/{int(r.n)} ({r.top1_accuracy:.1%}) | {int(r.top3_hits)}/{int(r.n)} ({r.top3_accuracy:.1%}) |"
        for r in groups.itertuples(index=False)
    )
    text = f"""# 成对混淆 marker 的 annotation 优化验证

## 严格验证设计

- 测试队列：新的 `{summary['n_test_samples']}` 个 Bo2023 未见样本，排除了此前使用的 `{summary['n_prior_samples_excluded']}` 个测试样本。
- 外层预测：保留 `correlation` 为唯一排名路径，annotation 不改变 Top1/Top3。
- 折内建模：每个外层 fold 的训练样本再划分为 marker 构建集和 calibration 集；外层测试样本不参与 marker、混淆对或阈值设定。
- 成对 marker：仅针对构建集 centroid correlation 的近邻混淆脑区，按方向一致性与 pairwise 效应量筛选，每个 Region 最多 `{summary['max_markers_per_region']}` 个 marker。
- 阈值校准：每个 fold 仅使用内部 calibration 样本选择 `supported` 与 `contradicted` 阈值。

## Correlation 基线

| 指标 | 结果 |
| --- | ---: |
| Top1 | {summary['baseline']['top1_hits']}/{summary['n_test_samples']} ({summary['baseline']['top1_accuracy']:.1%}) |
| Top3 | {summary['baseline']['top3_hits']}/{summary['n_test_samples']} ({summary['baseline']['top3_accuracy']:.1%}) |
| Macro AUC | {summary['baseline']['macro_auc']:.3f} |

## Annotation 分层结果

| 状态 | 样本数 | Top1 命中 | Top3 命中 |
| --- | ---: | ---: | ---: |
{group_lines}

## 覆盖与判定

- 平均每折 marker 对数：`{summary['marker_coverage']['mean_marker_pairs']:.1f}`；平均可评分混淆方向数：`{summary['marker_coverage']['mean_pair_directions']:.1f}`。
- 产生有效 annotation score 的测试样本：`{summary['annotation']['n_with_marker_score']}/{summary['n_test_samples']}`。
- Supported 组相对总体 Top1 富集：`{summary['annotation']['supported_top1_delta_pp']:+.1f}` 个百分点。

{summary['decision']}

## 下一步

{summary['next_step']}
"""
    (outdir / "pairwise_marker_annotation_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict validation of pairwise marker confidence annotation for Bo2023.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260529)
    parser.add_argument("--reference-topk-per-region", type=int, default=80)
    parser.add_argument("--calibration-fraction", type=float, default=0.20)
    parser.add_argument("--neighbors-per-region", type=int, default=5)
    parser.add_argument("--max-markers-per-pair", type=int, default=6)
    parser.add_argument("--max-markers-per-region", type=int, default=18)
    parser.add_argument("--min-effect", type=float, default=0.80)
    parser.add_argument("--min-consistency", type=float, default=0.75)
    parser.add_argument("--candidate-k", type=int, default=10)
    parser.add_argument("--exclude-samples", type=Path, action="append", default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "fold_pairwise_markers").mkdir(parents=True, exist_ok=True)
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
    sample_ids = matrix.columns.tolist()
    sample_pos = {sample_id: j for j, sample_id in enumerate(sample_ids)}
    region_pos = {region: j for j, region in enumerate(regions)}
    detail_rows: list[dict[str, Any]] = []
    probability_rows: list[dict[str, Any]] = []
    manifest: list[dict[str, Any]] = []

    for fold_no, heldout in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(heldout.sample_id)
        truth = str(heldout.region_id)
        test_idx = sample_pos[sample_id]
        truth_j = region_pos[truth]
        reference = full_reference.copy()
        train_idx = region_indices[truth][region_indices[truth] != test_idx]
        reference[:, truth_j] = values[:, train_idx].mean(axis=1, dtype=np.float64).astype(np.float32)
        outer_mask, _ = select_fold_signature(genes, regions, reference, args.reference_topk_per_region)
        test_detail, probability, test_order = evaluate_correlation_prediction(
            sample_id, truth, values[:, test_idx], reference, regions, outer_mask
        )
        build_indices, calibration_indices = split_inner_training_indices(
            region_indices,
            sample_ids,
            test_idx,
            seed=args.seed + fold_no,
            calibration_fraction=args.calibration_fraction,
        )
        build_reference = build_reference_from_indices(values, regions, build_indices)
        inner_mask, _ = select_fold_signature(
            genes, regions, build_reference, args.reference_topk_per_region
        )
        neighbors = derive_confusion_neighbors(
            build_reference,
            regions,
            build_indices,
            neighbors_per_region=args.neighbors_per_region,
        )
        marker_df = build_pairwise_marker_signature(
            values,
            genes,
            build_indices,
            neighbors,
            max_markers_per_pair=args.max_markers_per_pair,
            max_markers_per_region=args.max_markers_per_region,
            min_effect=args.min_effect,
            min_consistency=args.min_consistency,
        )
        calibration_rows: list[dict[str, Any]] = []
        for idx in calibration_indices:
            cal_truth = str(ann.loc[ann["sample_id"] == sample_ids[int(idx)], "region_id"].iloc[0])
            cal_detail, _, cal_order = evaluate_correlation_prediction(
                sample_ids[int(idx)],
                cal_truth,
                values[:, int(idx)],
                build_reference,
                regions,
                inner_mask,
            )
            candidates = [regions[int(j)] for j in cal_order[1 : args.candidate_k]]
            score, n_pairs, n_markers = pairwise_annotation_score(
                values[:, int(idx)], cal_detail["pred_top1"], candidates, marker_df
            )
            calibration_rows.append(
                {
                    "annotation_score": score,
                    "hit1": cal_detail["hit1"],
                    "n_pairs": n_pairs,
                    "n_markers": n_markers,
                }
            )
        calibration_df = pd.DataFrame(calibration_rows)
        thresholds = calibrate_annotation_thresholds(calibration_df)
        test_candidates = [regions[int(j)] for j in test_order[1 : args.candidate_k]]
        test_score, n_pairs, n_markers = pairwise_annotation_score(
            values[:, test_idx], test_detail["pred_top1"], test_candidates, marker_df
        )
        status = apply_annotation(test_score, thresholds)
        test_detail.update(
            {
                "annotation_status": status,
                "annotation_score": test_score,
                "annotation_pairs": n_pairs,
                "annotation_markers": n_markers,
                "support_threshold": thresholds["support_threshold"],
                "contradiction_threshold": thresholds["contradiction_threshold"],
                "n_valid_calibration": thresholds["n_valid_calibration"],
                "calibration_accuracy": thresholds["calibration_accuracy"],
            }
        )
        detail_rows.append(test_detail)
        probability_rows.append(probability)
        marker_file = args.outdir / "fold_pairwise_markers" / f"fold_{fold_no:02d}_{sample_id}_pairwise_markers.csv"
        marker_df.to_csv(marker_file, index=False, encoding="utf-8-sig")
        manifest.append(
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "truth_region": truth,
                "n_inner_calibration_samples": int(len(calibration_indices)),
                "n_valid_calibration_scores": int(thresholds["n_valid_calibration"]),
                "n_marker_pairs": int(len(marker_df)),
                "n_pair_directions": int(marker_df[["region_id", "competitor_region"]].drop_duplicates().shape[0])
                if len(marker_df)
                else 0,
                "annotation_status": status,
                "marker_file": str(marker_file.relative_to(args.outdir)),
            }
        )

    detail = pd.DataFrame(detail_rows)
    probabilities = pd.DataFrame(probability_rows)
    baseline = {
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probabilities[regions])),
    }
    annotation_summary = (
        detail.groupby("annotation_status", dropna=False)
        .agg(n=("hit1", "size"), top1_hits=("hit1", "sum"), top1_accuracy=("hit1", "mean"), top3_hits=("hit3", "sum"), top3_accuracy=("hit3", "mean"))
        .reset_index()
    )
    supported = annotation_summary[annotation_summary["annotation_status"] == "supported"]
    supported_n = int(supported["n"].iloc[0]) if len(supported) else 0
    supported_acc = float(supported["top1_accuracy"].iloc[0]) if len(supported) else float("nan")
    supported_delta_pp = (
        (supported_acc - baseline["top1_accuracy"]) * 100 if np.isfinite(supported_acc) else float("nan")
    )
    if supported_n >= 5 and supported_delta_pp >= 10.0:
        decision = "成对 marker annotation 在本轮形成了更窄且明显富集正确预测的 supported 组，可作为实验性置信度标记进入下一批复核；仍不改变 correlation 排名。"
        next_step = "用另一批未见样本重复验证同一预设规则；若 supported 富集再次成立，再考虑在界面输出 marker-supported 标记和人工复核提示。"
    else:
        decision = "成对 marker annotation 未形成足够稳定的高可信 supported 分层，暂不接入正式输出。"
        next_step = "继续针对具体混淆对压缩 marker 集并检查 pair 覆盖，随后重新进行未见队列验证。"
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "outer unseen LOSO with inner build/calibration split for annotation only",
        "seed": int(args.seed),
        "n_test_samples": int(len(detail)),
        "n_prior_samples_excluded": int(len(excluded_ids)),
        "n_regions": int(len(regions)),
        "n_gene_symbols": int(matrix.shape[0]),
        "n_singleton_samples_excluded": int(len(singleton_samples)),
        "neighbors_per_region": int(args.neighbors_per_region),
        "max_markers_per_pair": int(args.max_markers_per_pair),
        "max_markers_per_region": int(args.max_markers_per_region),
        "min_effect": float(args.min_effect),
        "min_consistency": float(args.min_consistency),
        "baseline": baseline,
        "marker_coverage": {
            "mean_marker_pairs": float(np.mean([x["n_marker_pairs"] for x in manifest])),
            "mean_pair_directions": float(np.mean([x["n_pair_directions"] for x in manifest])),
        },
        "annotation": {
            "n_with_marker_score": int(np.isfinite(detail["annotation_score"]).sum()),
            "supported_n": supported_n,
            "supported_top1_accuracy": supported_acc,
            "supported_top1_delta_pp": supported_delta_pp,
        },
        "decision": decision,
        "next_step": next_step,
    }
    detail.to_csv(args.outdir / "pairwise_annotation_detail.csv", index=False, encoding="utf-8-sig")
    annotation_summary.to_csv(args.outdir / "annotation_summary.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(args.outdir / "selected_test_samples.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "fold_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_annotation_plot(args.outdir, baseline["top1_accuracy"], annotation_summary)
    write_report(args.outdir, summary, annotation_summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
