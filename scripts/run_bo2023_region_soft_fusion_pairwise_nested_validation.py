#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
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
    softmax,
)
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


DEFAULT_NETWORK_DETAIL = (
    ROOT
    / "results"
    / "bo2023_network_pairwise_correlation_full_loso_819_rerun_20260526"
    / "network_pairwise_correlation_rescue_top3_detail.csv"
)
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_region_soft_fusion_pairwise_nested_loso_814_20260526"
BASELINE_ROUTE = "top3_beam_local_discriminative_top200_baseline"
SOFT_ROUTE = "top3_beam_network_soft_fusion"
PAIR_ROUTE = "top3_beam_network_soft_fusion_region_pairwise"


def stable_order(value: str) -> int:
    return int(hashlib.sha1(value.encode("utf-8")).hexdigest()[:12], 16)


def build_training_indices(
    region_indices: dict[str, np.ndarray],
    excluded: set[int],
) -> dict[str, np.ndarray]:
    return {
        region: indices[~np.isin(indices, np.asarray(sorted(excluded), dtype=int))]
        for region, indices in region_indices.items()
    }


def make_inner_split(
    ann: pd.DataFrame,
    sample_pos: dict[str, int],
    outer_idx: int,
) -> tuple[set[int], list[int]]:
    build: set[int] = set()
    calibration: list[int] = []
    for _, group in ann.groupby("region_id"):
        indices = [sample_pos[str(sample_id)] for sample_id in group["sample_id"] if sample_pos[str(sample_id)] != outer_idx]
        indices = sorted(indices, key=lambda idx: stable_order(str(idx)))
        n_calibration = min(max(1, int(round(len(indices) * 0.20))), max(0, len(indices) - 2))
        if n_calibration:
            calibration.extend(indices[:n_calibration])
            build.update(indices[n_calibration:])
        else:
            build.update(indices)
    return build, sorted(calibration)


def build_reference_from_indices(
    values: np.ndarray,
    regions: list[str],
    training: dict[str, np.ndarray],
) -> np.ndarray:
    columns: list[np.ndarray] = []
    for region in regions:
        idx = training[region]
        if not len(idx):
            columns.append(np.zeros(values.shape[0], dtype=np.float32))
        else:
            columns.append(values[:, idx].mean(axis=1, dtype=np.float64).astype(np.float32))
    return np.column_stack(columns)


def select_region_gene_pool(
    values: np.ndarray,
    regions: list[str],
    training: dict[str, np.ndarray],
    top_n: int,
) -> np.ndarray:
    eligible = {region: idx for region, idx in training.items() if len(idx) >= 2}
    rows, _ = select_group_discriminative_genes(values, sorted(eligible), eligible, top_n)
    return rows


def network_probabilities(
    values: np.ndarray,
    network_labels: np.ndarray,
    networks: list[str],
    heldout_idx: int,
    sample_idx: int,
    corrected_top3: list[str] | None = None,
    top_n_genes: int = 200,
) -> tuple[list[str], dict[str, float]]:
    reference, training = build_group_reference(values, network_labels, networks, heldout_idx)
    rows, _ = select_group_discriminative_genes(values, networks, training, top_n_genes)
    scores = correlation_scores(reference, values[:, sample_idx], rows)
    natural_order = [networks[int(j)] for j in np.argsort(scores)[::-1]]
    order = natural_order
    adjusted = scores.copy()
    if corrected_top3:
        order = corrected_top3 + [network for network in natural_order if network not in corrected_top3]
        sorted_scores = np.sort(scores)[::-1]
        pos = {network: j for j, network in enumerate(networks)}
        for rank, network in enumerate(order):
            adjusted[pos[network]] = sorted_scores[rank]
    probabilities = softmax(adjusted)
    return order, {network: float(probabilities[j]) for j, network in enumerate(networks)}


def build_network_model_from_indices(
    values: np.ndarray,
    network_labels: np.ndarray,
    networks: list[str],
    build_indices: set[int],
    top_n_genes: int = 200,
) -> tuple[np.ndarray, np.ndarray]:
    training = {
        network: np.asarray(
            [idx for idx in sorted(build_indices) if network_labels[idx] == network],
            dtype=int,
        )
        for network in networks
    }
    reference = np.column_stack(
        [values[:, training[network]].mean(axis=1, dtype=np.float64).astype(np.float32) for network in networks]
    )
    rows, _ = select_group_discriminative_genes(values, networks, training, top_n_genes)
    return reference, rows


def score_network_model(
    reference: np.ndarray,
    rows: np.ndarray,
    sample: np.ndarray,
    networks: list[str],
) -> tuple[list[str], dict[str, float]]:
    scores = correlation_scores(reference, sample, rows)
    order = [networks[int(j)] for j in np.argsort(scores)[::-1]]
    probabilities = softmax(scores)
    return order, {network: float(probabilities[j]) for j, network in enumerate(networks)}


def region_network_prior(
    candidate_regions: list[str],
    training_ann: pd.DataFrame,
    network_probability: dict[str, float],
) -> dict[str, float]:
    memberships = (
        training_ann[training_ann["region_id"].isin(candidate_regions)]
        .groupby("region_id")["endpoint_label"]
        .apply(lambda values: sorted(set(values.astype(str))))
    )
    return {
        region: max(network_probability.get(network, 1e-12) for network in memberships.get(region, []))
        for region in candidate_regions
    }


def fused_scores(
    region_scores: np.ndarray,
    regions: list[str],
    candidate_regions: list[str],
    priors: dict[str, float],
    alpha: float,
) -> np.ndarray:
    out = region_scores.copy()
    region_pos = {region: j for j, region in enumerate(regions)}
    for region in candidate_regions:
        out[region_pos[region]] += float(alpha) * math.log(max(priors.get(region, 1e-12), 1e-12))
    return out


def discover_region_pairs(
    calibration_predictions: pd.DataFrame,
    max_pairs_per_truth: int,
    min_pair_errors: int,
) -> set[tuple[str, str]]:
    errors = calibration_predictions[calibration_predictions["label"] != calibration_predictions["pred_top1"]]
    errors = errors[errors["pred_top1"].astype(str).ne("")]
    counts = Counter((str(row.label), str(row.pred_top1)) for row in errors.itertuples(index=False))
    selected: set[tuple[str, str]] = set()
    for truth in sorted(errors["label"].astype(str).unique().tolist()):
        candidates = sorted(
            [(pred, count) for (label, pred), count in counts.items() if label == truth],
            key=lambda item: (-item[1], item[0]),
        )
        retained = 0
        for predicted, count in candidates:
            if count >= min_pair_errors and retained < max_pairs_per_truth:
                selected.add(tuple(sorted((truth, predicted))))
                retained += 1
    return selected


def build_region_pair_models(
    values: np.ndarray,
    pairs: set[tuple[str, str]],
    training: dict[str, np.ndarray],
    gene_pool: np.ndarray,
    top_n_genes: int,
) -> tuple[dict[tuple[str, str], np.ndarray], pd.DataFrame]:
    models: dict[tuple[str, str], np.ndarray] = {}
    rows: list[dict[str, Any]] = []
    for left, right in sorted(pairs):
        if len(training.get(left, [])) < 2 or len(training.get(right, [])) < 2:
            continue
        left_values = values[gene_pool[:, None], training[left]]
        right_values = values[gene_pool[:, None], training[right]]
        delta = left_values.mean(axis=1) - right_values.mean(axis=1)
        pooled = (left_values.var(axis=1, ddof=1) + right_values.var(axis=1, ddof=1)) / 2.0
        effect = np.abs(delta) / np.sqrt(pooled + 1e-8)
        chosen = np.argsort(effect)[::-1][: min(top_n_genes, len(gene_pool))]
        models[(left, right)] = gene_pool[chosen].astype(int)
        rows.append(
            {
                "left_region": left,
                "right_region": right,
                "n_genes": int(len(chosen)),
                "mean_effect": float(effect[chosen].mean()),
            }
        )
    return models, pd.DataFrame(rows)


def apply_pairwise_adjustment(
    scores: np.ndarray,
    reference: np.ndarray,
    sample: np.ndarray,
    regions: list[str],
    candidate_regions: list[str],
    pair_models: dict[tuple[str, str], np.ndarray],
    beta: float,
) -> np.ndarray:
    out = scores.copy()
    region_pos = {region: j for j, region in enumerate(regions)}
    candidates = set(candidate_regions)
    for (left, right), genes in pair_models.items():
        if left not in candidates or right not in candidates:
            continue
        pair_scores = correlation_scores(reference, sample, genes)
        margin = float(pair_scores[region_pos[left]] - pair_scores[region_pos[right]])
        out[region_pos[left]] += float(beta) * margin / 2.0
        out[region_pos[right]] -= float(beta) * margin / 2.0
    return out


def evaluate_detail(
    route: str,
    sample_id: str,
    truth_region: str,
    truth_network: str,
    network_top: list[str],
    scores: np.ndarray,
    regions: list[str],
    candidate_regions: list[str],
) -> dict[str, Any]:
    candidate_indices = np.asarray([regions.index(region) for region in candidate_regions], dtype=int)
    ranked = rank_candidates(scores, regions, candidate_indices)
    return score_ranked_route(
        route,
        sample_id,
        truth_region,
        truth_network,
        " | ".join(network_top),
        ranked,
        len(regions),
        int(network_top[0] == truth_network),
        int(truth_network in network_top),
    )


def calibration_score(detail: pd.DataFrame) -> tuple[float, float]:
    return float(detail["hit1"].mean()), float(detail["hit3"].mean())


def export_plot(outdir: Path, metrics: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["Beam baseline", "Network soft fusion", "Soft fusion + Region pairwise"]
    values = metrics.set_index("route").loc[[BASELINE_ROUTE, SOFT_ROUTE, PAIR_ROUTE], ["top1_accuracy", "top3_accuracy"]].to_numpy()
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(10.4, 5.3), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, values[:, 0], width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, values[:, 1], width, label="Top3", color="#009E73")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, float(values.max()) + 0.12)
    ax.set_ylabel("Accuracy")
    ax.set_title("Bo2023 Region nested LOSO: soft fusion and pair-specific correlation", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "region_soft_fusion_pairwise_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "region_soft_fusion_pairwise_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    baseline = summary["routes"][BASELINE_ROUTE]
    soft = summary["routes"][SOFT_ROUTE]
    pair = summary["routes"][PAIR_ROUTE]
    change_soft = summary["paired_changes"][SOFT_ROUTE]
    change_pair = summary["paired_changes"][PAIR_ROUTE]
    text = f"""# Top3 Network 软融合与 Region Pairwise Correlation 嵌套 LOSO 验证

## 设计

- 数据：Bo2023 VSD，`{summary['n_test_samples']}` 个可评估 Region 样本，排除 `{summary['n_singleton_samples_excluded']}` 个 singleton Region 样本。
- 外层：逐一样本留出评估，外层测试样本不参与 reference、局部基因、混淆对或融合参数选择。
- 内层：每个外层训练折再划分构建集和校准集；在校准集上选择 `alpha`（Network 概率软融合权重）与 `beta`（Region pair-specific correlation 权重），并发现高频 Region 混淆对。
- 一级候选：使用严格 LOSO 产出的 pairwise-corrected SaleemNetworks Top3 排名。
- 当前基线：`Top3 Network beam + 局部 Top200 判别基因 Region correlation`，即上一轮的 `21.3% / 41.6%` 路线。

## 评分

`soft_score = local_region_correlation + alpha * log(network_probability)`

在 Region pairwise 路线中，再对内层训练发现的高频混淆 Region 对加入 pair-specific correlation margin，权重为内层选出的 `beta`。

## 结果

| 路径 | Top1 | Top3 | 真实 Region 中位排名 |
| --- | ---: | ---: | ---: |
| Top3 beam 局部 correlation 基线 | {baseline['top1_hits']}/{summary['n_test_samples']} ({baseline['top1_accuracy']:.1%}) | {baseline['top3_hits']}/{summary['n_test_samples']} ({baseline['top3_accuracy']:.1%}) | {baseline['median_true_rank']:.1f} |
| Network 概率软融合 | {soft['top1_hits']}/{summary['n_test_samples']} ({soft['top1_accuracy']:.1%}) | {soft['top3_hits']}/{summary['n_test_samples']} ({soft['top3_accuracy']:.1%}) | {soft['median_true_rank']:.1f} |
| 软融合 + Region pairwise correlation | {pair['top1_hits']}/{summary['n_test_samples']} ({pair['top1_accuracy']:.1%}) | {pair['top3_hits']}/{summary['n_test_samples']} ({pair['top3_accuracy']:.1%}) | {pair['median_true_rank']:.1f} |

## 相对当前基线的配对变化

| 路径 | Top1 新增 | Top1 丢失 | Top1 p | Top3 新增 | Top3 丢失 | Top3 p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Network 概率软融合 | {change_soft['top1_gains']} | {change_soft['top1_losses']} | {summary['paired_pvalues'][SOFT_ROUTE]['top1']:.3f} | {change_soft['top3_gains']} | {change_soft['top3_losses']} | {summary['paired_pvalues'][SOFT_ROUTE]['top3']:.3f} |
| 软融合 + Region pairwise | {change_pair['top1_gains']} | {change_pair['top1_losses']} | {summary['paired_pvalues'][PAIR_ROUTE]['top1']:.3f} | {change_pair['top3_gains']} | {change_pair['top3_losses']} | {summary['paired_pvalues'][PAIR_ROUTE]['top3']:.3f} |

## 判定

{summary['decision']}
"""
    (outdir / "region_soft_fusion_pairwise_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Nested LOSO validation of Region soft fusion and pair-specific correlation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--pair-gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=50)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=2)
    parser.add_argument("--alpha-grid", type=float, nargs="+", default=[0.0, 0.05, 0.10, 0.20, 0.35])
    parser.add_argument("--beta-grid", type=float, nargs="+", default=[0.0, 0.05, 0.10, 0.20])
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
    samples = matrix.columns.astype(str).tolist()
    sample_pos = {sample_id: j for j, sample_id in enumerate(samples)}
    _, regions, _, region_indices = build_region_reference(values, samples, ann)
    network_labels = ann.set_index("sample_id").reindex(samples)["endpoint_label"].to_numpy(dtype=str)
    networks = sorted(set(network_labels))
    network_detail = pd.read_csv(args.network_detail).set_index("sample_id")

    region_counts = ann.groupby("region_id")["sample_id"].size()
    singleton_regions = set(region_counts[region_counts < 2].index)
    selected = ann[~ann["region_id"].isin(singleton_regions)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)

    routes: dict[str, list[dict[str, Any]]] = {BASELINE_ROUTE: [], SOFT_ROUTE: [], PAIR_ROUTE: []}
    fold_rows: list[dict[str, Any]] = []
    pair_audits: list[pd.DataFrame] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        test_idx = sample_pos[sample_id]
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        corrected_top3 = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]

        build_indices, calibration_indices = make_inner_split(ann, sample_pos, test_idx)
        build_ann = ann[ann["sample_id"].map(sample_pos).isin(build_indices)].copy()
        build_training = build_training_indices(region_indices, set(range(len(samples))) - build_indices)
        build_regions = [region for region in regions if len(build_training[region])]
        build_reference = build_reference_from_indices(values, regions, build_training)
        build_gene_pool = select_region_gene_pool(values, build_regions, build_training, args.pair_gene_pool_size)
        build_local_rows = build_gene_pool[: args.local_top_n_genes]
        inner_network_reference, inner_network_rows = build_network_model_from_indices(
            values, network_labels, networks, build_indices
        )

        preliminary_rows: list[dict[str, Any]] = []
        calibration_data: list[dict[str, Any]] = []
        for cal_idx in calibration_indices:
            cal_sample_id = samples[cal_idx]
            cal_truth = str(ann.loc[ann["sample_id"] == cal_sample_id, "region_id"].iloc[0])
            cal_network = str(ann.loc[ann["sample_id"] == cal_sample_id, "endpoint_label"].iloc[0])
            net_order, net_prob = score_network_model(
                inner_network_reference, inner_network_rows, values[:, cal_idx], networks
            )
            top3 = net_order[:3]
            candidates = sorted(build_ann.loc[build_ann["endpoint_label"].isin(top3), "region_id"].unique().tolist())
            if cal_truth not in candidates or not candidates:
                preliminary_rows.append({"label": cal_truth, "pred_top1": ""})
                continue
            local_scores = correlation_scores(build_reference, values[:, cal_idx], build_local_rows)
            priors = region_network_prior(candidates, build_ann, net_prob)
            preliminary = evaluate_detail(
                "calibration", cal_sample_id, cal_truth, cal_network, top3, local_scores, regions, candidates
            )
            preliminary_rows.append(preliminary)
            calibration_data.append(
                {
                    "sample_id": cal_sample_id,
                    "sample_idx": cal_idx,
                    "truth_region": cal_truth,
                    "truth_network": cal_network,
                    "network_top": top3,
                    "candidates": candidates,
                    "local_scores": local_scores,
                    "priors": priors,
                }
            )
        pairs = discover_region_pairs(
            pd.DataFrame(preliminary_rows),
            args.max_pairs_per_truth,
            args.min_pair_errors,
        )
        pair_models_build, pair_audit = build_region_pair_models(
            values, pairs, build_training, build_gene_pool, args.pair_top_n_genes
        )
        if len(pair_audit):
            pair_audit.insert(0, "fold", fold)
            pair_audit.insert(1, "sample_id", sample_id)
            pair_audits.append(pair_audit)

        best_key = (float("-inf"), float("-inf"), float("-inf"), float("-inf"))
        selected_alpha = 0.0
        selected_beta = 0.0
        for alpha in args.alpha_grid:
            for beta in args.beta_grid:
                rows: list[dict[str, Any]] = []
                for item in calibration_data:
                    fused = fused_scores(item["local_scores"], regions, item["candidates"], item["priors"], alpha)
                    adjusted = apply_pairwise_adjustment(
                        fused,
                        build_reference,
                        values[:, int(item["sample_idx"])],
                        regions,
                        item["candidates"],
                        pair_models_build,
                        beta,
                    )
                    rows.append(
                        evaluate_detail(
                            "calibration",
                            str(item["sample_id"]),
                            str(item["truth_region"]),
                            str(item["truth_network"]),
                            list(item["network_top"]),
                            adjusted,
                            regions,
                            list(item["candidates"]),
                        )
                    )
                if not rows:
                    continue
                top1, top3 = calibration_score(pd.DataFrame(rows))
                key = (top1, top3, -float(alpha), -float(beta))
                if key > best_key:
                    best_key = key
                    selected_alpha = float(alpha)
                    selected_beta = float(beta)

        outer_training = build_training_indices(region_indices, {test_idx})
        outer_reference = build_reference_from_indices(values, regions, outer_training)
        outer_ann = ann[ann["sample_id"] != sample_id].copy()
        beam_regions = sorted(outer_ann.loc[outer_ann["endpoint_label"].isin(corrected_top3), "region_id"].unique().tolist())
        beam_rows = build_local_discriminative_rows(values, beam_regions, region_indices, test_idx, args.local_top_n_genes)
        baseline_scores = correlation_scores(outer_reference, values[:, test_idx], beam_rows)
        routes[BASELINE_ROUTE].append(
            evaluate_detail(BASELINE_ROUTE, sample_id, truth_region, truth_network, corrected_top3, baseline_scores, regions, beam_regions)
        )
        _, outer_prob = network_probabilities(
            values, network_labels, networks, test_idx, test_idx, corrected_top3=corrected_top3
        )
        outer_priors = region_network_prior(beam_regions, outer_ann, outer_prob)
        soft_scores = fused_scores(baseline_scores, regions, beam_regions, outer_priors, selected_alpha)
        routes[SOFT_ROUTE].append(
            evaluate_detail(SOFT_ROUTE, sample_id, truth_region, truth_network, corrected_top3, soft_scores, regions, beam_regions)
        )
        outer_pool = select_region_gene_pool(values, regions, outer_training, args.pair_gene_pool_size)
        pair_models_outer, _ = build_region_pair_models(
            values, pairs, outer_training, outer_pool, args.pair_top_n_genes
        )
        pair_scores = apply_pairwise_adjustment(
            soft_scores, outer_reference, values[:, test_idx], regions, beam_regions, pair_models_outer, selected_beta
        )
        routes[PAIR_ROUTE].append(
            evaluate_detail(PAIR_ROUTE, sample_id, truth_region, truth_network, corrected_top3, pair_scores, regions, beam_regions)
        )
        fold_rows.append(
            {
                "fold": fold,
                "sample_id": sample_id,
                "truth_region": truth_region,
                "truth_network": truth_network,
                "alpha": selected_alpha,
                "beta": selected_beta,
                "n_calibration": int(len(calibration_data)),
                "n_region_pairs": int(len(pair_models_outer)),
                "beam_regions": int(len(beam_regions)),
            }
        )

    details = {route: pd.DataFrame(rows) for route, rows in routes.items()}
    metrics = {route: summarize(frame) for route, frame in details.items()}
    baseline = details[BASELINE_ROUTE]
    changes = {route: paired_changes(baseline, details[route]) for route in [SOFT_ROUTE, PAIR_ROUTE]}
    pvalues = {
        route: {
            metric: paired_binomial_pvalue(change[f"{metric}_gains"], change[f"{metric}_losses"])
            for metric in ["top1", "top3"]
        }
        for route, change in changes.items()
    }
    pair = metrics[PAIR_ROUTE]
    base = metrics[BASELINE_ROUTE]
    decision = (
        "软融合加 Region pair-specific correlation 同时提高了当前 beam 基线的 Top1 与 Top3，"
        "可作为后续独立确认的候选路径。"
        if pair["top1_accuracy"] > base["top1_accuracy"] and pair["top3_accuracy"] > base["top3_accuracy"]
        else "本轮软融合加 Region pair-specific correlation 未同时超过当前 beam 基线，不接入正式候选路径。"
    )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "outer strict LOSO with fold-internal build/calibration for alpha, beta and Region confusion pairs",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singleton_regions).sum()),
        "alpha_grid": args.alpha_grid,
        "beta_grid": args.beta_grid,
        "local_top_n_genes": int(args.local_top_n_genes),
        "pair_gene_pool_size": int(args.pair_gene_pool_size),
        "pair_top_n_genes": int(args.pair_top_n_genes),
        "routes": metrics,
        "paired_changes": changes,
        "paired_pvalues": pvalues,
        "decision": decision,
    }
    metrics_frame = pd.DataFrame([{"route": route, **values} for route, values in metrics.items()])
    metrics_frame.to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    for route, frame in details.items():
        frame.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_tuning_detail.csv", index=False, encoding="utf-8-sig")
    if pair_audits:
        pd.concat(pair_audits, ignore_index=True).to_csv(args.outdir / "fold_region_pair_models.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics_frame)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
