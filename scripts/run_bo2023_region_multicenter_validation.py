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
from sklearn.cluster import KMeans


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_bo2023_hierarchical_region_correlation_validation import (  # noqa: E402
    build_local_discriminative_rows,
    paired_binomial_pvalue,
    paired_changes,
    rank_candidates,
    score_ranked_route,
    summarize,
)
from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    correlation_scores,
    read_annotations,
    read_vsd_matrix,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_NETWORK_DETAIL = (
    ROOT
    / "results"
    / "bo2023_network_pairwise_correlation_full_loso_819_rerun_20260526"
    / "network_pairwise_correlation_rescue_top3_detail.csv"
)
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_region_multicenter_full_loso_814_20260527"
BEAM_SINGLE = "top3_beam_local_single_centroid"
BEAM_MULTI = "top3_beam_local_multicenter"
ORACLE_SINGLE = "oracle_network_local_single_centroid"
ORACLE_MULTI = "oracle_network_local_multicenter"


def candidate_training_indices(
    candidate_regions: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
) -> dict[str, np.ndarray]:
    return {
        region: region_indices[region][region_indices[region] != heldout_idx]
        for region in candidate_regions
    }


def selected_reference(
    values: np.ndarray,
    candidate_regions: list[str],
    training: dict[str, np.ndarray],
    rows: np.ndarray,
) -> np.ndarray:
    return np.column_stack(
        [values[rows[:, None], training[region]].mean(axis=1, dtype=np.float64) for region in candidate_regions]
    )


def resolution_labels(
    values: np.ndarray,
    candidate_regions: list[str],
    training: dict[str, np.ndarray],
    rows: np.ndarray,
    min_center_samples: int,
    similarity_threshold: float,
) -> dict[str, tuple[str, str, float]]:
    reference = selected_reference(values, candidate_regions, training, rows)
    labels: dict[str, tuple[str, str, float]] = {}
    for j, region in enumerate(candidate_regions):
        count = len(training[region])
        competitors = [i for i in range(len(candidate_regions)) if i != j]
        nearest = float(np.max(correlation_scores(reference[:, competitors], reference[:, j]))) if competitors else float("nan")
        reasons: list[str] = []
        if count < min_center_samples:
            reasons.append(f"training_n<{min_center_samples}")
        if np.isfinite(nearest) and nearest >= similarity_threshold:
            reasons.append(f"nearest_centroid_corr>={similarity_threshold:.2f}")
        status = "low_resolution" if reasons else "high_resolution"
        labels[region] = (status, ";".join(reasons), nearest)
    return labels


def multicenter_scores(
    values: np.ndarray,
    sample: np.ndarray,
    candidate_regions: list[str],
    regions: list[str],
    training: dict[str, np.ndarray],
    rows: np.ndarray,
    min_center_samples: int,
    min_cluster_size: int,
) -> tuple[np.ndarray, dict[str, int]]:
    scores = np.full(len(regions), -np.inf, dtype=float)
    region_pos = {region: j for j, region in enumerate(regions)}
    n_centers: dict[str, int] = {}
    sample_local = sample[rows]
    for region in candidate_regions:
        x = values[rows[:, None], training[region]].astype(float, copy=False)
        centers = x.mean(axis=1, keepdims=True)
        if x.shape[1] >= min_center_samples:
            model = KMeans(n_clusters=2, random_state=0, n_init=10)
            assignments = model.fit_predict(x.T)
            counts = np.bincount(assignments, minlength=2)
            if int(counts.min()) >= min_cluster_size:
                centers = model.cluster_centers_.T
        n_centers[region] = int(centers.shape[1])
        scores[region_pos[region]] = float(np.max(correlation_scores(centers, sample_local)))
    return scores, n_centers


def annotated_detail(
    route: str,
    sample_id: str,
    truth_region: str,
    truth_network: str,
    network_top: list[str],
    scores: np.ndarray,
    regions: list[str],
    candidate_regions: list[str],
    resolution: dict[str, tuple[str, str, float]],
    n_centers: dict[str, int],
) -> dict[str, Any]:
    detail = score_ranked_route(
        route,
        sample_id,
        truth_region,
        truth_network,
        " | ".join(network_top),
        rank_candidates(scores, regions, np.asarray([regions.index(region) for region in candidate_regions], dtype=int)),
        len(regions),
        int(network_top[0] == truth_network),
        int(truth_network in network_top),
    )
    predicted = str(detail["pred_top1"])
    pred_status, pred_reason, pred_similarity = resolution[predicted]
    if truth_region in resolution:
        truth_status, truth_reason, truth_similarity = resolution[truth_region]
    else:
        truth_status, truth_reason, truth_similarity = "outside_candidate_set", "outside_network_beam", float("nan")
    detail.update(
        {
            "pred_top1_resolution": pred_status,
            "pred_top1_resolution_reason": pred_reason,
            "pred_top1_nearest_centroid_corr": pred_similarity,
            "pred_top1_n_centers": int(n_centers.get(predicted, 1)),
            "truth_resolution": truth_status,
            "truth_resolution_reason": truth_reason,
            "truth_nearest_centroid_corr": truth_similarity,
            "truth_n_centers": int(n_centers.get(truth_region, 0)),
        }
    )
    return detail


def export_plot(outdir: Path, metrics: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [BEAM_SINGLE, BEAM_MULTI, ORACLE_SINGLE, ORACLE_MULTI]
    labels = ["Top3 beam\nsingle", "Top3 beam\nmulticenter", "Oracle\nsingle", "Oracle\nmulticenter"]
    values = metrics.set_index("route").loc[routes, ["top1_accuracy", "top3_accuracy"]].to_numpy(dtype=float)
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(10.6, 5.5), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, values[:, 0], width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, values[:, 1], width, label="Top3", color="#009E73")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Accuracy")
    ax.set_ylim(0, float(values.max()) + 0.13)
    ax.set_title("Bo2023 Region multicenter reference: strict LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "region_multicenter_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "region_multicenter_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    routes = summary["routes"]
    beam_single = routes[BEAM_SINGLE]
    beam_multi = routes[BEAM_MULTI]
    oracle_single = routes[ORACLE_SINGLE]
    oracle_multi = routes[ORACLE_MULTI]
    beam_change = summary["paired_changes"]["beam_multicenter_vs_single"]
    oracle_change = summary["paired_changes"]["oracle_multicenter_vs_single"]
    text = f"""# Region 折内多中心 Reference 严格 LOSO 验证

## 设计

- 数据：Bo2023 VSD；精确 Region `110` 类，严格可评估样本 `{summary['n_test_samples']}` 个；`{summary['n_singleton_samples_excluded']}` 个 singleton Region 不参与准确率评估。
- 候选入口：复用 pairwise-corrected `SaleemNetworks` Top3 beam。
- 局部基因：每个外层 LOSO fold 内，仅用训练样本与当前候选 Region 构建 Top `{summary['local_top_n_genes']}` 局部判别基因。
- 多中心规则：训练折内 Region 样本数不少于 `{summary['min_center_samples']}` 时尝试 `2` 个 KMeans 表达中心；若任一簇少于 `{summary['min_cluster_size']}` 个训练样本，则回退为单中心。
- 低分辨率标记：仅依据训练折判断；训练样本少于 `{summary['min_center_samples']}`，或在当前候选空间中与最近 Region centroid 的相关系数不低于 `{summary['similarity_threshold']:.2f}`，标记为 `low_resolution`。
- Oracle 对照：使用真实 Network 限定 Region 候选，用于观察二级 Region reference 本身的区分潜力，不属于可部署预测。

## 准确率

| 路径 | Top1 | Top3 | 真实 Region 中位排名 |
| --- | ---: | ---: | ---: |
| Top3 beam + 单中心 reference | {beam_single['top1_hits']}/{summary['n_test_samples']} ({beam_single['top1_accuracy']:.1%}) | {beam_single['top3_hits']}/{summary['n_test_samples']} ({beam_single['top3_accuracy']:.1%}) | {beam_single['median_true_rank']:.1f} |
| Top3 beam + 多中心 reference | {beam_multi['top1_hits']}/{summary['n_test_samples']} ({beam_multi['top1_accuracy']:.1%}) | {beam_multi['top3_hits']}/{summary['n_test_samples']} ({beam_multi['top3_accuracy']:.1%}) | {beam_multi['median_true_rank']:.1f} |
| Oracle Network + 单中心 reference | {oracle_single['top1_hits']}/{summary['n_test_samples']} ({oracle_single['top1_accuracy']:.1%}) | {oracle_single['top3_hits']}/{summary['n_test_samples']} ({oracle_single['top3_accuracy']:.1%}) | {oracle_single['median_true_rank']:.1f} |
| Oracle Network + 多中心 reference | {oracle_multi['top1_hits']}/{summary['n_test_samples']} ({oracle_multi['top1_accuracy']:.1%}) | {oracle_multi['top3_hits']}/{summary['n_test_samples']} ({oracle_multi['top3_accuracy']:.1%}) | {oracle_multi['median_true_rank']:.1f} |

## 多中心配对变化

| 比较 | Top1 新增 | Top1 丢失 | Top1 p | Top3 新增 | Top3 丢失 | Top3 p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Top3 beam 多中心 vs 单中心 | {beam_change['top1_gains']} | {beam_change['top1_losses']} | {summary['paired_pvalues']['beam_multicenter_vs_single']['top1']:.3f} | {beam_change['top3_gains']} | {beam_change['top3_losses']} | {summary['paired_pvalues']['beam_multicenter_vs_single']['top3']:.3f} |
| Oracle 多中心 vs 单中心 | {oracle_change['top1_gains']} | {oracle_change['top1_losses']} | {summary['paired_pvalues']['oracle_multicenter_vs_single']['top1']:.3f} | {oracle_change['top3_gains']} | {oracle_change['top3_losses']} | {summary['paired_pvalues']['oracle_multicenter_vs_single']['top3']:.3f} |

## 低分辨率候选

Top3 beam 多中心路线输出的 Top1 候选中，`{summary['resolution']['beam_multi_low_resolution_predictions']}/{summary['n_test_samples']}` 个被训练侧规则标记为 `low_resolution`；标记为 `high_resolution` 的预测 Top1 命中率为 `{summary['resolution']['beam_multi_high_resolution_top1']:.1%}`，`low_resolution` 组为 `{summary['resolution']['beam_multi_low_resolution_top1']:.1%}`。该分组用于提示结论粒度和人工复核需求，不改变总体准确率口径。

## 判定

{summary['decision']}
"""
    (outdir / "region_multicenter_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict LOSO validation of fold-local multicenter Region references.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--min-center-samples", type=int, default=8)
    parser.add_argument("--min-cluster-size", type=int, default=2)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = map_matrix_to_symbols(read_vsd_matrix(args.matrix), args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    net = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet, usecols=["No.", args.network_col])
    net["sample_id"] = net["No."].astype(str).str.strip()
    net["endpoint_label"] = net[args.network_col].fillna("NA").astype(str).str.strip()
    ann = ann.merge(net[["sample_id", "endpoint_label"]], on="sample_id", how="left")
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    values = matrix.to_numpy(dtype=np.float32)
    sample_ids = matrix.columns.astype(str).tolist()
    sample_pos = {sample_id: j for j, sample_id in enumerate(sample_ids)}
    reference_all, regions, _, region_indices = build_region_reference(values, sample_ids, ann)
    region_pos = {region: j for j, region in enumerate(regions)}
    network_detail = pd.read_csv(args.network_detail).set_index("sample_id")

    region_counts = ann.groupby("region_id")["sample_id"].size()
    singleton_regions = set(region_counts[region_counts < 2].index)
    selected = ann[~ann["region_id"].isin(singleton_regions)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)

    routes: dict[str, list[dict[str, Any]]] = {route: [] for route in [BEAM_SINGLE, BEAM_MULTI, ORACLE_SINGLE, ORACLE_MULTI]}
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id]
        network_top = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        beam_regions = sorted(training_ann.loc[training_ann["endpoint_label"].isin(network_top), "region_id"].unique().tolist())
        oracle_regions = sorted(training_ann.loc[training_ann["endpoint_label"] == truth_network, "region_id"].unique().tolist())

        reference = reference_all.copy()
        reference[:, region_pos[truth_region]] = values[
            :, region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        ].mean(axis=1, dtype=np.float64).astype(np.float32)

        row_data: dict[str, Any] = {"fold": fold, "sample_id": sample_id, "truth_region": truth_region}
        for scope, candidates, single_route, multi_route, network_labels in [
            ("beam", beam_regions, BEAM_SINGLE, BEAM_MULTI, network_top),
            ("oracle", oracle_regions, ORACLE_SINGLE, ORACLE_MULTI, [truth_network]),
        ]:
            local_rows = build_local_discriminative_rows(values, candidates, region_indices, heldout_idx, args.local_top_n_genes)
            if local_rows is None:
                local_rows = np.arange(values.shape[0], dtype=int)
            training = candidate_training_indices(candidates, region_indices, heldout_idx)
            resolution = resolution_labels(
                values,
                candidates,
                training,
                local_rows,
                args.min_center_samples,
                args.similarity_threshold,
            )
            single_scores = correlation_scores(reference, sample, local_rows)
            single_centers = {region: 1 for region in candidates}
            routes[single_route].append(
                annotated_detail(
                    single_route, sample_id, truth_region, truth_network, network_labels,
                    single_scores, regions, candidates, resolution, single_centers
                )
            )
            multi_scores, n_centers = multicenter_scores(
                values,
                sample,
                candidates,
                regions,
                training,
                local_rows,
                args.min_center_samples,
                args.min_cluster_size,
            )
            routes[multi_route].append(
                annotated_detail(
                    multi_route, sample_id, truth_region, truth_network, network_labels,
                    multi_scores, regions, candidates, resolution, n_centers
                )
            )
            row_data[f"{scope}_n_candidates"] = int(len(candidates))
            row_data[f"{scope}_n_multicenter_regions"] = int(sum(value == 2 for value in n_centers.values()))
            row_data[f"{scope}_truth_resolution"] = resolution.get(truth_region, ("outside_candidate_set", "", float("nan")))[0]
        fold_rows.append(row_data)

    details = {route: pd.DataFrame(rows) for route, rows in routes.items()}
    metrics = {route: summarize(frame) for route, frame in details.items()}
    change_sets = {
        "beam_multicenter_vs_single": paired_changes(details[BEAM_SINGLE], details[BEAM_MULTI]),
        "oracle_multicenter_vs_single": paired_changes(details[ORACLE_SINGLE], details[ORACLE_MULTI]),
    }
    pvalues = {
        name: {
            metric: paired_binomial_pvalue(changes[f"{metric}_gains"], changes[f"{metric}_losses"])
            for metric in ["top1", "top3"]
        }
        for name, changes in change_sets.items()
    }
    beam_multi = details[BEAM_MULTI]
    high = beam_multi[beam_multi["pred_top1_resolution"] == "high_resolution"]
    low = beam_multi[beam_multi["pred_top1_resolution"] == "low_resolution"]
    improved_beam = metrics[BEAM_MULTI]["top1_accuracy"] > metrics[BEAM_SINGLE]["top1_accuracy"] and metrics[BEAM_MULTI]["top3_accuracy"] >= metrics[BEAM_SINGLE]["top3_accuracy"]
    improved_oracle = metrics[ORACLE_MULTI]["top1_accuracy"] > metrics[ORACLE_SINGLE]["top1_accuracy"] and metrics[ORACLE_MULTI]["top3_accuracy"] >= metrics[ORACLE_SINGLE]["top3_accuracy"]
    if improved_beam and improved_oracle:
        decision = "折内多中心 reference 在 Top3 beam 与 Oracle 条件下均提高 Region Top1 且未降低 Top3，可进入独立样本确认；正式接入前仍保留当前单中心路线为默认。"
    elif improved_oracle:
        decision = "多中心 reference 改善了 Oracle 条件下的 Region 区分能力，但未稳定转化为 Top3 beam 可部署增益；暂不替换当前单中心 Region 路线。"
    else:
        decision = "本轮固定规则多中心 reference 未在 Oracle 或 Top3 beam 条件下形成同时改善 Top1/Top3 的证据，不接入正式候选路径。"
    summary: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; fold-local local genes, two-center Region reference and resolution flags",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singleton_regions).sum()),
        "local_top_n_genes": int(args.local_top_n_genes),
        "min_center_samples": int(args.min_center_samples),
        "min_cluster_size": int(args.min_cluster_size),
        "similarity_threshold": float(args.similarity_threshold),
        "routes": metrics,
        "paired_changes": change_sets,
        "paired_pvalues": pvalues,
        "resolution": {
            "beam_multi_high_resolution_predictions": int(len(high)),
            "beam_multi_low_resolution_predictions": int(len(low)),
            "beam_multi_high_resolution_top1": float(high["hit1"].mean()) if len(high) else float("nan"),
            "beam_multi_low_resolution_top1": float(low["hit1"].mean()) if len(low) else float("nan"),
        },
        "decision": decision,
    }
    metrics_frame = pd.DataFrame([{"route": route, **values} for route, values in metrics.items()])
    metrics_frame.to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    for route, frame in details.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_multicenter_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics_frame)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
