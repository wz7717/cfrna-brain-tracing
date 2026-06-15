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

from scripts.run_bo2023_hierarchical_region_correlation_validation import (  # noqa: E402
    build_local_discriminative_rows,
    rank_candidates,
)
from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    correlation_scores,
    read_annotations,
    read_vsd_matrix,
)
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import (  # noqa: E402
    BEAM_ROUTE,
    DEFAULT_NETWORK_DETAIL,
    ORACLE_ROUTE,
    build_resolution_groups,
    candidate_training_indices,
    region_network_assignment,
    score_route,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_group_specific_correlation_nested_loso_814_20260527"
BEAM_GROUP_ROUTE = "top3_network_beam_group_specific_correlation"
ORACLE_GROUP_ROUTE = "oracle_network_group_specific_correlation"


def resolution_group_training(
    annotations: dict[str, dict[str, Any]],
    training: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    members: dict[str, list[str]] = {}
    for region, annotation in annotations.items():
        members.setdefault(str(annotation["resolution_group"]), []).append(region)
    return {
        group: np.unique(np.concatenate([training[region] for region in regions])).astype(int)
        for group, regions in members.items()
    }


def rank_groups_directly(
    values: np.ndarray,
    sample: np.ndarray,
    group_training: dict[str, np.ndarray],
    top_n_genes: int,
) -> tuple[list[str], np.ndarray]:
    groups = sorted(group_training)
    if len(groups) < 2:
        return groups, np.arange(values.shape[0], dtype=int)
    rows, _ = select_group_discriminative_genes(values, groups, group_training, top_n_genes)
    reference = np.column_stack(
        [values[rows[:, None], group_training[group]].mean(axis=1, dtype=np.float64) for group in groups]
    )
    scores = correlation_scores(reference, sample[rows])
    order = np.argsort(scores)[::-1]
    return [groups[int(j)] for j in order], rows


def score_direct_group_route(
    route: str,
    sample_id: str,
    truth_region: str,
    truth_network: str,
    network_values: list[str],
    ranked_groups: list[str],
    annotations: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    truth_group = (
        str(annotations[truth_region]["resolution_group"])
        if truth_region in annotations
        else "outside_network_beam"
    )
    true_rank = ranked_groups.index(truth_group) + 1 if truth_group in ranked_groups else len(ranked_groups) + 1
    padded = ranked_groups[:3] + [""] * max(0, 3 - len(ranked_groups))
    return {
        "route": route,
        "sample_id": sample_id,
        "label": truth_region,
        "true_network": truth_network,
        "network_beam": " | ".join(network_values),
        "network_top1_hit": int(network_values[0] == truth_network),
        "network_top3_hit": int(truth_network in network_values),
        "true_resolution_group": truth_group,
        "pred_group_top1": padded[0],
        "pred_group_top2": padded[1],
        "pred_group_top3": padded[2],
        "group_true_rank": int(true_rank),
        "group_hit1": int(true_rank == 1),
        "group_hit3": int(true_rank <= 3),
        "truth_resolution_tier": (
            annotations[truth_region]["resolution_tier"] if truth_region in annotations else "outside_network_beam"
        ),
        "n_candidate_groups": int(len(ranked_groups)),
    }


def summarize_groups(detail: pd.DataFrame) -> dict[str, Any]:
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
    discordant = gains + losses
    if discordant == 0:
        return 1.0
    tail = min(gains, losses)
    probability = sum(math.comb(discordant, k) for k in range(tail + 1)) / (2**discordant)
    return float(min(1.0, 2.0 * probability))


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    routes = [BEAM_ROUTE, BEAM_GROUP_ROUTE, ORACLE_ROUTE, ORACLE_GROUP_ROUTE]
    labels = ["Beam\nRegion-first group", "Beam\nGroup-specific", "Oracle\nRegion-first group", "Oracle\nGroup-specific"]
    top1 = [metrics[route]["group_top1_accuracy"] for route in routes]
    top3 = [metrics[route]["group_top3_accuracy"] for route in routes]
    x = np.arange(len(routes))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11.6, 5.8), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Top3", color="#009E73")
    ax.set_ylabel("Resolvable group accuracy")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, min(1.0, max(top3) + 0.14))
    ax.set_title("Bo2023 group-specific discriminative correlation: strict nested LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "group_specific_correlation_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "group_specific_correlation_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    metrics = summary["routes"]
    beam_base = metrics[BEAM_ROUTE]
    beam_test = metrics[BEAM_GROUP_ROUTE]
    oracle_base = metrics[ORACLE_ROUTE]
    oracle_test = metrics[ORACLE_GROUP_ROUTE]
    beam_change = summary["paired_changes"]["beam_group_specific_vs_region_first"]
    oracle_change = summary["paired_changes"]["oracle_group_specific_vs_region_first"]
    text = f"""# Group-specific 判别基因 Correlation 严格嵌套 LOSO 验证

## 设计

- 在每个外层 LOSO 训练折中，先用既有 Region 局部判别路线学习 `resolvable Region group`，测试样本不参与分组。
- 基线：按精确 Region correlation 排名，再将排名映射为首次出现的 group。
- 新路线：将训练折中的合并组直接作为标签，使用 `group-vs-group` Top `{summary['group_top_n_genes']}` 判别基因建立组 centroid 并对组直接排名。
- 一级候选仍固定复用 pairwise-corrected `SaleemNetworks` Top3 beam；Oracle 仅用于诊断。

## 结果

| 条件 | 评分路径 | Group Top1 | Group Top3 | 中位真实组排名 |
| --- | --- | ---: | ---: | ---: |
| Top3 beam | Region-first group baseline | {beam_base['group_top1_hits']}/{summary['n_test_samples']} ({beam_base['group_top1_accuracy']:.1%}) | {beam_base['group_top3_hits']}/{summary['n_test_samples']} ({beam_base['group_top3_accuracy']:.1%}) | {beam_base['median_group_true_rank']:.1f} |
| Top3 beam | Group-specific correlation | {beam_test['group_top1_hits']}/{summary['n_test_samples']} ({beam_test['group_top1_accuracy']:.1%}) | {beam_test['group_top3_hits']}/{summary['n_test_samples']} ({beam_test['group_top3_accuracy']:.1%}) | {beam_test['median_group_true_rank']:.1f} |
| Oracle Network | Region-first group baseline | {oracle_base['group_top1_hits']}/{summary['n_test_samples']} ({oracle_base['group_top1_accuracy']:.1%}) | {oracle_base['group_top3_hits']}/{summary['n_test_samples']} ({oracle_base['group_top3_accuracy']:.1%}) | {oracle_base['median_group_true_rank']:.1f} |
| Oracle Network | Group-specific correlation | {oracle_test['group_top1_hits']}/{summary['n_test_samples']} ({oracle_test['group_top1_accuracy']:.1%}) | {oracle_test['group_top3_hits']}/{summary['n_test_samples']} ({oracle_test['group_top3_accuracy']:.1%}) | {oracle_test['median_group_true_rank']:.1f} |

## 配对变化

| 比较 | Top1 新增/丢失 | Top1 p | Top3 新增/丢失 | Top3 p |
| --- | ---: | ---: | ---: | ---: |
| Top3 beam 新路线 vs 基线 | {beam_change['top1_gains']} / {beam_change['top1_losses']} | {summary['paired_pvalues']['beam_group_specific_vs_region_first']['top1']:.3f} | {beam_change['top3_gains']} / {beam_change['top3_losses']} | {summary['paired_pvalues']['beam_group_specific_vs_region_first']['top3']:.3f} |
| Oracle 新路线 vs 基线 | {oracle_change['top1_gains']} / {oracle_change['top1_losses']} | {summary['paired_pvalues']['oracle_group_specific_vs_region_first']['top1']:.3f} | {oracle_change['top3_gains']} / {oracle_change['top3_losses']} | {summary['paired_pvalues']['oracle_group_specific_vs_region_first']['top3']:.3f} |

## 判定

{summary['decision']}
"""
    (outdir / "group_specific_correlation_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict nested LOSO of group-specific Region correlation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--group-top-n-genes", type=int, default=200)
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
    singletons = set(counts[counts < 2].index)
    selected = ann[~ann["region_id"].isin(singletons)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)
    if args.max_samples is not None:
        selected = selected.head(max(1, int(args.max_samples))).copy()

    rows: dict[str, list[dict[str, Any]]] = {
        BEAM_ROUTE: [],
        BEAM_GROUP_ROUTE: [],
        ORACLE_ROUTE: [],
        ORACLE_GROUP_ROUTE: [],
    }
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        network_top = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)
        fold_info: dict[str, Any] = {"fold": fold, "sample_id": sample_id, "truth_region": truth_region}
        for scope, network_values, baseline_route, test_route in [
            ("beam", network_top, BEAM_ROUTE, BEAM_GROUP_ROUTE),
            ("oracle", [truth_network], ORACLE_ROUTE, ORACLE_GROUP_ROUTE),
        ]:
            candidates = sorted(
                training_ann.loc[training_ann["endpoint_label"].isin(network_values), "region_id"].unique().tolist()
            )
            local_rows = build_local_discriminative_rows(
                values, candidates, region_indices, heldout_idx, args.local_top_n_genes
            )
            if local_rows is None:
                local_rows = np.arange(values.shape[0], dtype=int)
            training = candidate_training_indices(candidates, region_indices, heldout_idx)
            annotations, _ = build_resolution_groups(
                values,
                candidates,
                training,
                region_network_assignment(training_ann, candidates),
                local_rows,
                args.min_resolution_samples,
                args.min_merge_samples,
                args.min_pair_errors,
                args.min_confusion_rate,
                args.similarity_threshold,
                args.merge_similarity_threshold,
                args.max_group_size,
            )
            baseline_scores = correlation_scores(reference, sample, local_rows)
            candidate_indices = np.asarray([region_pos[region] for region in candidates], dtype=int)
            ranked_regions = rank_candidates(baseline_scores, regions, candidate_indices)
            rows[baseline_route].append(
                score_route(
                    baseline_route, sample_id, truth_region, truth_network, network_values,
                    ranked_regions, annotations, len(regions)
                )
            )
            group_training = resolution_group_training(annotations, training)
            ranked_groups, group_rows = rank_groups_directly(
                values, sample, group_training, args.group_top_n_genes
            )
            rows[test_route].append(
                score_direct_group_route(
                    test_route, sample_id, truth_region, truth_network, network_values,
                    ranked_groups, annotations
                )
            )
            fold_info[f"{scope}_n_groups"] = int(len(group_training))
            fold_info[f"{scope}_n_group_genes"] = int(len(group_rows))
        fold_rows.append(fold_info)

    details = {route: pd.DataFrame(data) for route, data in rows.items()}
    metrics = {route: summarize_groups(detail) for route, detail in details.items()}
    changes = {
        "beam_group_specific_vs_region_first": paired_changes(details[BEAM_ROUTE], details[BEAM_GROUP_ROUTE]),
        "oracle_group_specific_vs_region_first": paired_changes(details[ORACLE_ROUTE], details[ORACLE_GROUP_ROUTE]),
    }
    pvalues = {
        comparison: {
            metric: paired_pvalue(change[f"{metric}_gains"], change[f"{metric}_losses"])
            for metric in ["top1", "top3"]
        }
        for comparison, change in changes.items()
    }
    beam_base = metrics[BEAM_ROUTE]
    beam_test = metrics[BEAM_GROUP_ROUTE]
    improve = (
        beam_test["group_top1_accuracy"] > beam_base["group_top1_accuracy"]
        and beam_test["group_top3_accuracy"] >= beam_base["group_top3_accuracy"]
    )
    if improve:
        decision = (
            "Group-specific 判别基因 correlation 在可部署 Top3 beam 下同时提高 Group Top1 且未降低 Top3；"
            "进入下一步 Network 概率加权融合验证，但在后续配对确认前不替换正式输出。"
        )
    else:
        decision = (
            "Group-specific 判别基因 correlation 未在可部署 Top3 beam 下同时改善 Group Top1/Top3；"
            "不进入正式路径，也不以此路线作为概率融合基础。"
        )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO; fold-local resolution groups and fold-local group-specific discriminative genes",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singletons).sum()),
        "local_top_n_genes": int(args.local_top_n_genes),
        "group_top_n_genes": int(args.group_top_n_genes),
        "routes": metrics,
        "paired_changes": changes,
        "paired_pvalues": pvalues,
        "proceed_to_probability_fusion": bool(improve),
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **value} for route, value in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
