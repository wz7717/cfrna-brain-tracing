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

from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    correlation_scores,
    read_annotations,
    read_vsd_matrix,
)
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_NETWORK_DETAIL = (
    ROOT
    / "results"
    / "bo2023_network_discriminative_correlation_full_loso_819"
    / "network_discriminative_correlation_top200_detail.csv"
)
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_hierarchical_region_correlation_full_loso_814"


def rank_candidates(
    scores: np.ndarray,
    regions: list[str],
    candidate_indices: np.ndarray | None = None,
) -> list[str]:
    indices = np.arange(len(regions), dtype=int) if candidate_indices is None else np.asarray(candidate_indices, dtype=int)
    order = indices[np.argsort(scores[indices])[::-1]]
    return [regions[int(j)] for j in order]


def score_ranked_route(
    route: str,
    sample_id: str,
    truth_region: str,
    truth_network: str,
    predicted_network: str,
    ranked: list[str],
    n_total_regions: int,
    network_top1_hit: int,
    network_top3_hit: int,
) -> dict[str, Any]:
    true_rank = ranked.index(truth_region) + 1 if truth_region in ranked else n_total_regions + 1
    padded = ranked[:3] + [""] * max(0, 3 - len(ranked))
    return {
        "route": route,
        "sample_id": sample_id,
        "label": truth_region,
        "true_network": truth_network,
        "predicted_network": predicted_network,
        "network_top1_hit": int(network_top1_hit),
        "network_top3_hit": int(network_top3_hit),
        "pred_top1": padded[0],
        "pred_top2": padded[1],
        "pred_top3": padded[2],
        "true_rank": int(true_rank),
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "n_candidate_regions": int(len(ranked)),
    }


def summarize(detail: pd.DataFrame) -> dict[str, float | int]:
    correct_network = detail[detail["network_top1_hit"] == 1]
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "median_true_rank": float(detail["true_rank"].median()),
        "n_when_network_top1_correct": int(len(correct_network)),
        "conditional_top1_given_network_top1": float(correct_network["hit1"].mean()) if len(correct_network) else float("nan"),
        "conditional_top3_given_network_top1": float(correct_network["hit3"].mean()) if len(correct_network) else float("nan"),
    }


def build_local_discriminative_rows(
    values: np.ndarray,
    candidate_regions: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int,
    top_n_genes: int,
) -> np.ndarray | None:
    training: dict[str, np.ndarray] = {}
    for region in candidate_regions:
        indices = region_indices[region]
        indices = indices[indices != heldout_idx]
        if len(indices):
            training[region] = indices
    if len(training) < 2:
        return None
    rows, _ = select_group_discriminative_genes(values, sorted(training), training, top_n_genes)
    return rows


def export_plot(outdir: Path, metrics: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [
        "flat_region_correlation_all_genes",
        "hierarchical_top1_network_local_discriminative_top200",
        "hierarchical_top3_network_beam_local_discriminative_top200",
        "oracle_true_network_local_discriminative_top200",
    ]
    labels = ["Flat Region", "Top1 Network -> Local", "Top3 beam -> Local", "Oracle -> Local"]
    plot_df = metrics.set_index("route").loc[routes]
    values = plot_df[["top1_accuracy", "top3_accuracy"]].to_numpy(dtype=float)
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11.8, 5.4), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, values[:, 0], width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, values[:, 1], width, label="Top3", color="#009E73")
    ax.set_ylim(0, min(1.0, float(values.max()) + 0.15))
    ax.set_xticks(x, labels)
    ax.set_ylabel("Accuracy")
    ax.set_title("Bo2023 exact Region: hierarchical correlation strict LOSO (n=814)", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, height + 0.012, f"{height:.1%}", ha="center", fontsize=9)
    fig.savefig(outdir / "hierarchical_region_correlation_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "hierarchical_region_correlation_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    routes = summary["routes"]
    flat = routes["flat_region_correlation_all_genes"]
    hierarchical = routes["hierarchical_top1_network_region_correlation"]
    local = routes["hierarchical_top1_network_local_discriminative_top200"]
    beam = routes["hierarchical_top3_network_beam_region_correlation"]
    beam_local = routes["hierarchical_top3_network_beam_local_discriminative_top200"]
    oracle = routes["oracle_true_network_region_correlation"]
    oracle_local = routes["oracle_true_network_local_discriminative_top200"]
    network = summary["network_gate"]
    text = f"""# 分层 Region Correlation 严格 LOSO 验证

## 设计

- 数据集：Bo2023 VSD，精确 Region 共 `{summary['n_regions']}` 类。
- 可评估样本：`{summary['n_test_samples']}` 个；排除 `{summary['n_singleton_samples_excluded']}` 个 singleton Region 样本，因为留出后没有同 Region 训练参考。
- 一级预测：复用已经完成的全量严格 LOSO `{summary['network_route_description']}` 预测；每个测试样本的一级预测未使用自身表达构建参考。
- 二级预测：仅在训练折中出现于一级 Top1 Network 的 Region 中重新计算 Pearson correlation；少数跨 Network 的同名 Region 不被强行归入唯一上层类别。
- 局部判别变体：仅使用当前一级候选 Network 内、由训练折构建的 Top `{summary['local_top_n_genes']}` 判别基因重新计算 Pearson correlation。
- Top3 Network beam：将一级预测前三个 Network 在训练折中出现过的 Region 合并为候选集，再以该候选集内部折内 Top `{summary['local_top_n_genes']}` 判别基因进行 correlation 排名。
- Oracle 诊断：假设一级 Network 已知为真，再进行 Region correlation，用于估计二级区分本身的上限，不属于可部署路线。

## 结果

| 路径 | Top1 | Top3 | 真实 Region 中位排名 |
| --- | ---: | ---: | ---: |
| 全部 Region 直接 correlation | {flat['top1_hits']}/{summary['n_test_samples']} ({flat['top1_accuracy']:.1%}) | {flat['top3_hits']}/{summary['n_test_samples']} ({flat['top3_accuracy']:.1%}) | {flat['median_true_rank']:.1f} |
| Top1 Network 内 Region correlation | {hierarchical['top1_hits']}/{summary['n_test_samples']} ({hierarchical['top1_accuracy']:.1%}) | {hierarchical['top3_hits']}/{summary['n_test_samples']} ({hierarchical['top3_accuracy']:.1%}) | {hierarchical['median_true_rank']:.1f} |
| Top1 Network 内局部判别基因 correlation | {local['top1_hits']}/{summary['n_test_samples']} ({local['top1_accuracy']:.1%}) | {local['top3_hits']}/{summary['n_test_samples']} ({local['top3_accuracy']:.1%}) | {local['median_true_rank']:.1f} |
| Top3 Network beam 内 Region correlation | {beam['top1_hits']}/{summary['n_test_samples']} ({beam['top1_accuracy']:.1%}) | {beam['top3_hits']}/{summary['n_test_samples']} ({beam['top3_accuracy']:.1%}) | {beam['median_true_rank']:.1f} |
| Top3 Network beam 内局部判别基因 correlation | {beam_local['top1_hits']}/{summary['n_test_samples']} ({beam_local['top1_accuracy']:.1%}) | {beam_local['top3_hits']}/{summary['n_test_samples']} ({beam_local['top3_accuracy']:.1%}) | {beam_local['median_true_rank']:.1f} |
| Oracle 真实 Network 内 Region correlation | {oracle['top1_hits']}/{summary['n_test_samples']} ({oracle['top1_accuracy']:.1%}) | {oracle['top3_hits']}/{summary['n_test_samples']} ({oracle['top3_accuracy']:.1%}) | {oracle['median_true_rank']:.1f} |
| Oracle 真实 Network 内局部判别基因 correlation | {oracle_local['top1_hits']}/{summary['n_test_samples']} ({oracle_local['top1_accuracy']:.1%}) | {oracle_local['top3_hits']}/{summary['n_test_samples']} ({oracle_local['top3_accuracy']:.1%}) | {oracle_local['median_true_rank']:.1f} |

一级 Network gate 在这 `{summary['n_test_samples']}` 个可评估样本上的准确率为 Top1 `{network['top1_hits']}/{summary['n_test_samples']} ({network['top1_accuracy']:.1%})`、Top3 `{network['top3_hits']}/{summary['n_test_samples']} ({network['top3_accuracy']:.1%})`。

## 配对变化

- Top1 Network 内 Region correlation 相比全局 Region correlation：新增命中 `{summary['paired_changes']['hierarchical_vs_flat']['top1_gains']}`，丢失命中 `{summary['paired_changes']['hierarchical_vs_flat']['top1_losses']}`；Top3 新增 `{summary['paired_changes']['hierarchical_vs_flat']['top3_gains']}`，丢失 `{summary['paired_changes']['hierarchical_vs_flat']['top3_losses']}`。
- Top1 Network 内局部判别基因 correlation 相比全局 Region correlation：Top1 新增 `{summary['paired_changes']['local_vs_flat']['top1_gains']}`，丢失 `{summary['paired_changes']['local_vs_flat']['top1_losses']}`（配对二项检验 `p={summary['paired_binomial_pvalues']['local_vs_flat']['top1']:.3f}`）；Top3 新增 `{summary['paired_changes']['local_vs_flat']['top3_gains']}`，丢失 `{summary['paired_changes']['local_vs_flat']['top3_losses']}`（`p={summary['paired_binomial_pvalues']['local_vs_flat']['top3']:.3f}`）。
- Top3 Network beam 局部判别基因 correlation 相比 Top1 Network 局部判别路线：Top1 新增 `{summary['paired_changes']['beam_local_vs_top1_local']['top1_gains']}`，丢失 `{summary['paired_changes']['beam_local_vs_top1_local']['top1_losses']}`（`p={summary['paired_binomial_pvalues']['beam_local_vs_top1_local']['top1']:.3f}`）；Top3 新增 `{summary['paired_changes']['beam_local_vs_top1_local']['top3_gains']}`，丢失 `{summary['paired_changes']['beam_local_vs_top1_local']['top3_losses']}`（`p={summary['paired_binomial_pvalues']['beam_local_vs_top1_local']['top3']:.3f}`）。
- Top3 Network beam 局部判别基因 correlation 相比全局 Region correlation：Top1 新增 `{summary['paired_changes']['beam_local_vs_flat']['top1_gains']}`，丢失 `{summary['paired_changes']['beam_local_vs_flat']['top1_losses']}`（`p={summary['paired_binomial_pvalues']['beam_local_vs_flat']['top1']:.3f}`）；Top3 新增 `{summary['paired_changes']['beam_local_vs_flat']['top3_gains']}`，丢失 `{summary['paired_changes']['beam_local_vs_flat']['top3_losses']}`（`p={summary['paired_binomial_pvalues']['beam_local_vs_flat']['top3']:.3f}`）。
- 局部判别基因变体相比普通分层 correlation：Top1 新增 `{summary['paired_changes']['local_vs_hierarchical']['top1_gains']}`，丢失 `{summary['paired_changes']['local_vs_hierarchical']['top1_losses']}`；Top3 新增 `{summary['paired_changes']['local_vs_hierarchical']['top3_gains']}`，丢失 `{summary['paired_changes']['local_vs_hierarchical']['top3_losses']}`。

## 判定

{summary['decision']}

## 解释边界

分层路径的 Region 结论受到一级 Network gate 上限约束。Oracle 结果仅表示一级分类正确时二级 Region correlation 可达到的潜力，不能直接作为生产准确率。
"""
    (outdir / "hierarchical_region_correlation_report_cn.md").write_text(text, encoding="utf-8")


def paired_changes(base: pd.DataFrame, candidate: pd.DataFrame) -> dict[str, int]:
    return {
        "top1_gains": int(((base["hit1"] == 0) & (candidate["hit1"] == 1)).sum()),
        "top1_losses": int(((base["hit1"] == 1) & (candidate["hit1"] == 0)).sum()),
        "top3_gains": int(((base["hit3"] == 0) & (candidate["hit3"] == 1)).sum()),
        "top3_losses": int(((base["hit3"] == 1) & (candidate["hit3"] == 0)).sum()),
    }


def paired_binomial_pvalue(gains: int, losses: int) -> float:
    discordant = gains + losses
    if discordant == 0:
        return 1.0
    tail = min(gains, losses)
    probability = sum(math.comb(discordant, k) for k in range(tail + 1)) / (2**discordant)
    return float(min(1.0, 2.0 * probability))


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict LOSO validation of hierarchical Network-to-Region correlation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument(
        "--network-route-description",
        default="SaleemNetworks Top200 判别基因 Pearson correlation",
    )
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    network_ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet, usecols=["No.", args.network_col])
    network_ann["sample_id"] = network_ann["No."].astype(str).str.strip()
    network_ann["endpoint_label"] = network_ann[args.network_col].fillna("NA").astype(str).str.strip()
    ann = ann.merge(network_ann[["sample_id", "endpoint_label"]], on="sample_id", how="left")
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    if ann["endpoint_label"].isna().any():
        raise ValueError("Missing SaleemNetworks labels in aligned annotations.")

    values = matrix.to_numpy(dtype=np.float32)
    sample_ids = matrix.columns.astype(str).tolist()
    sample_pos = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    reference_all, regions, region_counts, region_indices = build_region_reference(values, sample_ids, ann)
    region_pos = {region: j for j, region in enumerate(regions)}
    network_detail = pd.read_csv(args.network_detail).set_index("sample_id")
    missing_network_predictions = sorted(set(sample_ids) - set(network_detail.index.astype(str)))
    if missing_network_predictions:
        raise ValueError(f"Missing strict-LOSO network predictions for {len(missing_network_predictions)} samples.")

    counts_by_region = ann.groupby("region_id")["sample_id"].size()
    singleton_regions = set(counts_by_region[counts_by_region < 2].index)
    selected = ann[~ann["region_id"].isin(singleton_regions)].copy()
    selected["sample_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sample_order").drop(columns="sample_order").reset_index(drop=True)

    routes: dict[str, list[dict[str, Any]]] = {
        "flat_region_correlation_all_genes": [],
        "hierarchical_top1_network_region_correlation": [],
        "hierarchical_top1_network_local_discriminative_top200": [],
        "hierarchical_top3_network_beam_region_correlation": [],
        "hierarchical_top3_network_beam_local_discriminative_top200": [],
        "oracle_true_network_region_correlation": [],
        "oracle_true_network_local_discriminative_top200": [],
    }
    fold_rows: list[dict[str, Any]] = []
    for fold_no, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        reference = reference_all.copy()
        truth_j = region_pos[truth_region]
        truth_training = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, truth_j] = values[:, truth_training].mean(axis=1, dtype=np.float64).astype(np.float32)
        scores = correlation_scores(reference, sample)

        network_row = network_detail.loc[sample_id]
        network_top = [str(network_row[f"pred_top{i}"]) for i in [1, 2, 3]]
        predicted_network = network_top[0]
        network_top1_hit = int(predicted_network == truth_network)
        network_top3_hit = int(truth_network in network_top)
        training_ann = ann[ann["sample_id"] != sample_id]

        flat_ranked = rank_candidates(scores, regions)
        routes["flat_region_correlation_all_genes"].append(
            score_ranked_route(
                "flat_region_correlation_all_genes",
                sample_id,
                truth_region,
                truth_network,
                predicted_network,
                flat_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )

        candidate_regions = sorted(
            training_ann.loc[training_ann["endpoint_label"] == predicted_network, "region_id"].unique().tolist()
        )
        candidate_indices = np.asarray([region_pos[region] for region in candidate_regions], dtype=int)
        hierarchical_ranked = rank_candidates(scores, regions, candidate_indices)
        routes["hierarchical_top1_network_region_correlation"].append(
            score_ranked_route(
                "hierarchical_top1_network_region_correlation",
                sample_id,
                truth_region,
                truth_network,
                predicted_network,
                hierarchical_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )

        local_rows = build_local_discriminative_rows(
            values, candidate_regions, region_indices, heldout_idx, args.local_top_n_genes
        )
        local_scores = scores if local_rows is None else correlation_scores(reference, sample, local_rows)
        local_ranked = rank_candidates(local_scores, regions, candidate_indices)
        routes["hierarchical_top1_network_local_discriminative_top200"].append(
            score_ranked_route(
                "hierarchical_top1_network_local_discriminative_top200",
                sample_id,
                truth_region,
                truth_network,
                predicted_network,
                local_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )

        beam_regions = sorted(
            training_ann.loc[training_ann["endpoint_label"].isin(network_top), "region_id"].unique().tolist()
        )
        beam_indices = np.asarray([region_pos[region] for region in beam_regions], dtype=int)
        beam_ranked = rank_candidates(scores, regions, beam_indices)
        routes["hierarchical_top3_network_beam_region_correlation"].append(
            score_ranked_route(
                "hierarchical_top3_network_beam_region_correlation",
                sample_id,
                truth_region,
                truth_network,
                " | ".join(network_top),
                beam_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )
        beam_local_rows = build_local_discriminative_rows(
            values, beam_regions, region_indices, heldout_idx, args.local_top_n_genes
        )
        beam_local_scores = scores if beam_local_rows is None else correlation_scores(reference, sample, beam_local_rows)
        beam_local_ranked = rank_candidates(beam_local_scores, regions, beam_indices)
        routes["hierarchical_top3_network_beam_local_discriminative_top200"].append(
            score_ranked_route(
                "hierarchical_top3_network_beam_local_discriminative_top200",
                sample_id,
                truth_region,
                truth_network,
                " | ".join(network_top),
                beam_local_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )

        oracle_regions = sorted(
            training_ann.loc[training_ann["endpoint_label"] == truth_network, "region_id"].unique().tolist()
        )
        oracle_indices = np.asarray([region_pos[region] for region in oracle_regions], dtype=int)
        oracle_ranked = rank_candidates(scores, regions, oracle_indices)
        routes["oracle_true_network_region_correlation"].append(
            score_ranked_route(
                "oracle_true_network_region_correlation",
                sample_id,
                truth_region,
                truth_network,
                truth_network,
                oracle_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )
        oracle_local_rows = build_local_discriminative_rows(
            values, oracle_regions, region_indices, heldout_idx, args.local_top_n_genes
        )
        oracle_local_scores = scores if oracle_local_rows is None else correlation_scores(reference, sample, oracle_local_rows)
        oracle_local_ranked = rank_candidates(oracle_local_scores, regions, oracle_indices)
        routes["oracle_true_network_local_discriminative_top200"].append(
            score_ranked_route(
                "oracle_true_network_local_discriminative_top200",
                sample_id,
                truth_region,
                truth_network,
                truth_network,
                oracle_local_ranked,
                len(regions),
                network_top1_hit,
                network_top3_hit,
            )
        )
        fold_rows.append(
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "truth_region": truth_region,
                "truth_network": truth_network,
                "predicted_network": predicted_network,
                "network_top1_hit": network_top1_hit,
                "network_top3_hit": network_top3_hit,
                "n_predicted_network_regions": int(len(candidate_regions)),
                "n_local_discriminative_genes": int(len(local_rows)) if local_rows is not None else 0,
                "beam_networks": " | ".join(network_top),
                "n_beam_regions": int(len(beam_regions)),
                "n_beam_local_discriminative_genes": int(len(beam_local_rows)) if beam_local_rows is not None else 0,
            }
        )

    details = {route: pd.DataFrame(rows) for route, rows in routes.items()}
    metrics = {route: summarize(frame) for route, frame in details.items()}
    flat = details["flat_region_correlation_all_genes"]
    hierarchical = details["hierarchical_top1_network_region_correlation"]
    local = details["hierarchical_top1_network_local_discriminative_top200"]
    beam_local = details["hierarchical_top3_network_beam_local_discriminative_top200"]
    network_top1 = int(pd.DataFrame(fold_rows)["network_top1_hit"].sum())
    network_top3 = int(pd.DataFrame(fold_rows)["network_top3_hit"].sum())
    change_sets = {
        "hierarchical_vs_flat": paired_changes(flat, hierarchical),
        "local_vs_flat": paired_changes(flat, local),
        "local_vs_hierarchical": paired_changes(hierarchical, local),
        "beam_local_vs_top1_local": paired_changes(local, beam_local),
        "beam_local_vs_flat": paired_changes(flat, beam_local),
    }
    paired_pvalues = {
        comparison: {
            metric: paired_binomial_pvalue(values[f"{metric}_gains"], values[f"{metric}_losses"])
            for metric in ["top1", "top3"]
        }
        for comparison, values in change_sets.items()
    }
    beam_top1_p = paired_pvalues["beam_local_vs_top1_local"]["top1"]
    beam_top3_p = paired_pvalues["beam_local_vs_top1_local"]["top3"]
    beam_beats_top1 = metrics["hierarchical_top3_network_beam_local_discriminative_top200"]["top1_accuracy"] > metrics[
        "hierarchical_top1_network_local_discriminative_top200"
    ]["top1_accuracy"] and metrics["hierarchical_top3_network_beam_local_discriminative_top200"]["top3_accuracy"] >= metrics[
        "hierarchical_top1_network_local_discriminative_top200"
    ]["top3_accuracy"]
    if beam_beats_top1 and beam_top1_p < 0.05 and beam_top3_p < 0.05:
        decision = (
            "Top3 Network beam 局部判别 correlation 相比 Top1 Network 局部判别路线的 Top1 与 Top3 "
            f"均获得配对改善（Top1 p={beam_top1_p:.3f}，Top3 p={beam_top3_p:.3f}）。"
            "该路线可作为下一版 Region 二级候选的拟接入算法，但属于回顾性严格 LOSO 优化，需独立队列确认后再切换默认路径。"
        )
    elif beam_beats_top1:
        decision = (
            "Top3 Network beam 局部判别 correlation 数值上超过 Top1 Network 局部判别路线，"
            f"但配对证据尚不足（Top1 p={beam_top1_p:.3f}，Top3 p={beam_top3_p:.3f}）。"
            "因此不直接替换正式 Region 二级候选路径。"
        )
    else:
        decision = (
            "Top3 Network beam 未在 Top1 与 Top3 同时超过 Top1 Network 局部判别路线，"
            f"不替换正式 Region 二级候选路径（Top1 p={beam_top1_p:.3f}，Top3 p={beam_top3_p:.3f}）。"
        )
    summary: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict full evaluable LOSO; reused strict-LOSO fold-local Network predictions",
        "network_prediction_source": str(args.network_detail),
        "network_route_description": str(args.network_route_description),
        "n_regions": int(len(regions)),
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singleton_regions).sum()),
        "singleton_regions_excluded": sorted(singleton_regions),
        "local_top_n_genes": int(args.local_top_n_genes),
        "network_gate": {
            "top1_hits": network_top1,
            "top1_accuracy": float(network_top1 / len(selected)),
            "top3_hits": network_top3,
            "top3_accuracy": float(network_top3 / len(selected)),
        },
        "routes": metrics,
        "paired_changes": change_sets,
        "paired_binomial_pvalues": paired_pvalues,
        "decision": decision,
    }
    metric_frame = pd.DataFrame([{"route": route, **values} for route, values in metrics.items()])
    metric_frame.to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_detail.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(args.outdir / "evaluated_samples.csv", index=False, encoding="utf-8-sig")
    for route, frame in details.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metric_frame)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
