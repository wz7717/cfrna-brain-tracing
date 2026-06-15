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
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_resolution_tier_nested_loso_814_20260527"
DEFAULT_MODEL_OUT = ROOT / "data" / "models" / "bo2023_region_resolution_groups.json"
BEAM_ROUTE = "top3_network_beam_local_region_candidates"
ORACLE_ROUTE = "oracle_network_local_region_candidates"


class UnionFind:
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}
        self.members = {item: {item} for item in items}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: str, right: str, max_size: int) -> bool:
        a = self.find(left)
        b = self.find(right)
        if a == b:
            return True
        if len(self.members[a]) + len(self.members[b]) > max_size:
            return False
        if len(self.members[a]) < len(self.members[b]):
            a, b = b, a
        self.parent[b] = a
        self.members[a].update(self.members[b])
        del self.members[b]
        return True


def candidate_training_indices(
    candidates: list[str],
    region_indices: dict[str, np.ndarray],
    heldout_idx: int | None,
) -> dict[str, np.ndarray]:
    training: dict[str, np.ndarray] = {}
    for region in candidates:
        indices = region_indices[region]
        training[region] = indices if heldout_idx is None else indices[indices != heldout_idx]
    return training


def centroid_reference(
    values: np.ndarray,
    candidates: list[str],
    training: dict[str, np.ndarray],
    rows: np.ndarray,
) -> np.ndarray:
    return np.column_stack(
        [values[rows[:, None], training[region]].mean(axis=1, dtype=np.float64) for region in candidates]
    )


def region_network_assignment(training_ann: pd.DataFrame, candidates: list[str]) -> dict[str, str | None]:
    assignment: dict[str, str | None] = {}
    for region in candidates:
        networks = sorted(
            training_ann.loc[training_ann["region_id"] == region, "endpoint_label"].dropna().astype(str).unique().tolist()
        )
        assignment[region] = networks[0] if len(networks) == 1 else None
    return assignment


def build_resolution_groups(
    values: np.ndarray,
    candidates: list[str],
    training: dict[str, np.ndarray],
    region_network: dict[str, str | None],
    rows: np.ndarray,
    min_resolution_samples: int,
    min_merge_samples: int,
    min_pair_errors: int,
    min_confusion_rate: float,
    similarity_threshold: float,
    merge_similarity_threshold: float,
    max_group_size: int,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    annotations: dict[str, dict[str, Any]] = {}
    audit_rows: list[dict[str, Any]] = []
    reference = centroid_reference(values, candidates, training, rows)
    positions = {region: j for j, region in enumerate(candidates)}

    for region in candidates:
        annotations[region] = {
            "region_id": region,
            "network_id": region_network.get(region),
            "resolution_group": region,
            "group_members": [region],
            "resolution_tier": "high_resolution",
            "resolution_reasons": [],
            "training_samples": int(len(training[region])),
            "nearest_centroid_corr": None,
        }
        if len(training[region]) < min_resolution_samples:
            annotations[region]["resolution_reasons"].append(f"training_n<{min_resolution_samples}")
        if region_network.get(region) is None:
            annotations[region]["resolution_reasons"].append("cross_network_region")

    network_values = sorted({network for network in region_network.values() if network is not None})
    for network in network_values:
        network_regions = sorted([region for region in candidates if region_network.get(region) == network])
        if len(network_regions) < 2:
            continue
        network_positions = [positions[region] for region in network_regions]
        network_ref = reference[:, network_positions]
        for j, region in enumerate(network_regions):
            competitors = [i for i in range(len(network_regions)) if i != j]
            nearest = float(np.max(correlation_scores(network_ref[:, competitors], network_ref[:, j])))
            annotations[region]["nearest_centroid_corr"] = nearest
            if nearest >= similarity_threshold:
                annotations[region]["resolution_reasons"].append(
                    f"nearest_centroid_corr>={similarity_threshold:.2f}"
                )

        confusions: Counter[tuple[str, str]] = Counter()
        for truth in network_regions:
            truth_indices = training[truth]
            if len(truth_indices) <= 1:
                continue
            truth_pos = network_regions.index(truth)
            for sample_idx in truth_indices:
                local_ref = network_ref.copy()
                local_ref[:, truth_pos] = values[rows[:, None], truth_indices[truth_indices != sample_idx]].mean(
                    axis=1, dtype=np.float64
                )
                scores = correlation_scores(local_ref, values[rows, sample_idx])
                predicted = network_regions[int(np.argmax(scores))]
                if predicted != truth:
                    confusions[(truth, predicted)] += 1

        pair_candidates: list[dict[str, Any]] = []
        for left_i, left in enumerate(network_regions):
            for right in network_regions[left_i + 1 :]:
                errors = int(confusions[(left, right)] + confusions[(right, left)])
                pair_n = int(len(training[left]) + len(training[right]))
                rate = float(errors / pair_n) if pair_n else 0.0
                corr = float(correlation_scores(
                    reference[:, [positions[left], positions[right]]], reference[:, positions[left]]
                )[1])
                eligible = (
                    len(training[left]) >= min_merge_samples
                    and len(training[right]) >= min_merge_samples
                    and errors >= min_pair_errors
                    and rate >= min_confusion_rate
                    and corr >= merge_similarity_threshold
                )
                pair_candidates.append(
                    {
                        "network_id": network,
                        "left_region": left,
                        "right_region": right,
                        "pair_errors": errors,
                        "confusion_rate": rate,
                        "centroid_correlation": corr,
                        "eligible_for_merge": bool(eligible),
                    }
                )
        uf = UnionFind(network_regions)
        for pair in sorted(
            pair_candidates,
            key=lambda x: (-int(x["eligible_for_merge"]), -x["pair_errors"], -x["centroid_correlation"]),
        ):
            merged = False
            if pair["eligible_for_merge"]:
                merged = uf.union(str(pair["left_region"]), str(pair["right_region"]), max_group_size)
            pair["merged"] = bool(merged and pair["eligible_for_merge"])
            audit_rows.append(pair)

        groups: dict[str, list[str]] = {}
        for region in network_regions:
            groups.setdefault(uf.find(region), []).append(region)
        for members in groups.values():
            members = sorted(members)
            if len(members) <= 1:
                continue
            group_id = f"{network}::{' + '.join(members)}"
            for region in members:
                annotations[region]["resolution_group"] = group_id
                annotations[region]["group_members"] = members
                annotations[region]["resolution_reasons"].append("fold_local_merged_confusion_group")

    for annotation in annotations.values():
        reasons = list(dict.fromkeys(annotation["resolution_reasons"]))
        annotation["resolution_reasons"] = reasons
        annotation["resolution_tier"] = "low_resolution" if reasons else "high_resolution"
    return annotations, pd.DataFrame(audit_rows)


def distinct_ranked_groups(ranked_regions: list[str], annotations: dict[str, dict[str, Any]]) -> list[str]:
    groups: list[str] = []
    for region in ranked_regions:
        group = str(annotations[region]["resolution_group"])
        if group not in groups:
            groups.append(group)
    return groups


def score_route(
    route: str,
    sample_id: str,
    truth_region: str,
    truth_network: str,
    network_top: list[str],
    ranked_regions: list[str],
    annotations: dict[str, dict[str, Any]],
    n_total_regions: int,
) -> dict[str, Any]:
    truth_in_candidates = truth_region in annotations
    true_rank = ranked_regions.index(truth_region) + 1 if truth_region in ranked_regions else n_total_regions + 1
    predicted = ranked_regions[0]
    pred_annotation = annotations[predicted]
    ranked_groups = distinct_ranked_groups(ranked_regions, annotations)
    true_group = str(annotations[truth_region]["resolution_group"]) if truth_in_candidates else "outside_network_beam"
    group_rank = ranked_groups.index(true_group) + 1 if true_group in ranked_groups else n_total_regions + 1
    padded_regions = ranked_regions[:3] + [""] * max(0, 3 - len(ranked_regions))
    padded_groups = ranked_groups[:3] + [""] * max(0, 3 - len(ranked_groups))
    return {
        "route": route,
        "sample_id": sample_id,
        "label": truth_region,
        "true_network": truth_network,
        "network_beam": " | ".join(network_top),
        "network_top1_hit": int(network_top[0] == truth_network),
        "network_top3_hit": int(truth_network in network_top),
        "pred_top1": padded_regions[0],
        "pred_top2": padded_regions[1],
        "pred_top3": padded_regions[2],
        "true_rank": int(true_rank),
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "true_resolution_group": true_group,
        "pred_group_top1": padded_groups[0],
        "pred_group_top2": padded_groups[1],
        "pred_group_top3": padded_groups[2],
        "group_true_rank": int(group_rank),
        "group_hit1": int(group_rank == 1),
        "group_hit3": int(group_rank <= 3),
        "pred_top1_resolution_tier": pred_annotation["resolution_tier"],
        "pred_top1_resolution_reasons": ";".join(pred_annotation["resolution_reasons"]),
        "pred_top1_group_members": " | ".join(pred_annotation["group_members"]),
        "truth_resolution_tier": annotations[truth_region]["resolution_tier"] if truth_in_candidates else "outside_network_beam",
        "truth_group_members": " | ".join(annotations[truth_region]["group_members"]) if truth_in_candidates else "",
        "n_candidate_regions": int(len(ranked_regions)),
        "n_candidate_groups": int(len(ranked_groups)),
    }


def summarize(detail: pd.DataFrame) -> dict[str, Any]:
    low = detail[detail["pred_top1_resolution_tier"] == "low_resolution"]
    high = detail[detail["pred_top1_resolution_tier"] == "high_resolution"]
    return {
        "n": int(len(detail)),
        "exact_top1_hits": int(detail["hit1"].sum()),
        "exact_top1_accuracy": float(detail["hit1"].mean()),
        "exact_top3_hits": int(detail["hit3"].sum()),
        "exact_top3_accuracy": float(detail["hit3"].mean()),
        "group_top1_hits": int(detail["group_hit1"].sum()),
        "group_top1_accuracy": float(detail["group_hit1"].mean()),
        "group_top3_hits": int(detail["group_hit3"].sum()),
        "group_top3_accuracy": float(detail["group_hit3"].mean()),
        "median_exact_true_rank": float(detail["true_rank"].median()),
        "median_group_true_rank": float(detail["group_true_rank"].median()),
        "low_resolution_predictions": int(len(low)),
        "low_resolution_fraction": float(len(low) / len(detail)),
        "low_resolution_exact_top1": float(low["hit1"].mean()) if len(low) else float("nan"),
        "high_resolution_predictions": int(len(high)),
        "high_resolution_exact_top1": float(high["hit1"].mean()) if len(high) else float("nan"),
    }


def export_plot(outdir: Path, metrics: dict[str, dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["Top3 beam\nExact Region", "Top3 beam\nResolvable Group", "Oracle\nExact Region", "Oracle\nResolvable Group"]
    top1 = [
        metrics[BEAM_ROUTE]["exact_top1_accuracy"],
        metrics[BEAM_ROUTE]["group_top1_accuracy"],
        metrics[ORACLE_ROUTE]["exact_top1_accuracy"],
        metrics[ORACLE_ROUTE]["group_top1_accuracy"],
    ]
    top3 = [
        metrics[BEAM_ROUTE]["exact_top3_accuracy"],
        metrics[BEAM_ROUTE]["group_top3_accuracy"],
        metrics[ORACLE_ROUTE]["exact_top3_accuracy"],
        metrics[ORACLE_ROUTE]["group_top3_accuracy"],
    ]
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11.3, 5.8), constrained_layout=True)
    bars1 = ax.bar(x - width / 2, top1, width, label="Top1", color="#0072B2")
    bars3 = ax.bar(x + width / 2, top3, width, label="Top3", color="#009E73")
    ax.set_ylabel("Accuracy")
    ax.set_ylim(0, min(1.0, max(top3) + 0.14))
    ax.set_xticks(x, labels)
    ax.set_title("Bo2023 hierarchical resolution endpoint: strict nested LOSO", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend()
    for bars in [bars1, bars3]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012, f"{bar.get_height():.1%}", ha="center")
    fig.savefig(outdir / "resolution_tier_hierarchical_endpoint_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "resolution_tier_hierarchical_endpoint_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    beam = summary["routes"][BEAM_ROUTE]
    oracle = summary["routes"][ORACLE_ROUTE]
    text = f"""# Resolvable Region Group 分层终点严格嵌套 LOSO 验证

## 设计

- 数据：Bo2023 VSD；精确 Region `110` 类；可评估样本 `{summary['n_test_samples']}` 个，排除 `{summary['n_singleton_samples_excluded']}` 个 singleton Region 样本。
- 一级入口：已完成严格 LOSO 的 pairwise-corrected `SaleemNetworks` Top3 beam。
- 二级精确候选：每个外层折内，使用候选 Region 的训练样本建立 Top `{summary['local_top_n_genes']}` 局部判别基因 correlation 排名。
- Resolution group：仅在同一 Network 内，使用外层训练折的内层 LOO 混淆和 centroid 相似性识别反复不可分 Region；满足 `pair_errors >= {summary['min_pair_errors']}`、`confusion_rate >= {summary['min_confusion_rate']:.2f}`、`centroid_corr >= {summary['merge_similarity_threshold']:.2f}` 的 Region 才允许合并，每组最多 `{summary['max_group_size']}` 个 Region。
- 低分辨率标记：训练样本少于 `{summary['min_resolution_samples']}`、邻近 centroid correlation 不低于 `{summary['similarity_threshold']:.2f}`、跨 Network 同名 Region 或进入合并组时标记为 `low_resolution`。该标记只提示人工复核，不改变精确 Region 排名。
- Oracle：用真实 Network 限定候选，观察在一级入口无误时 resolution group 的可达到表现；不是可部署准确率。

## 分层终点结果

| 候选范围 | Exact Region Top1 | Exact Region Top3 | Resolvable Group Top1 | Resolvable Group Top3 |
| --- | ---: | ---: | ---: | ---: |
| Top3 Network beam | {beam['exact_top1_hits']}/{summary['n_test_samples']} ({beam['exact_top1_accuracy']:.1%}) | {beam['exact_top3_hits']}/{summary['n_test_samples']} ({beam['exact_top3_accuracy']:.1%}) | {beam['group_top1_hits']}/{summary['n_test_samples']} ({beam['group_top1_accuracy']:.1%}) | {beam['group_top3_hits']}/{summary['n_test_samples']} ({beam['group_top3_accuracy']:.1%}) |
| Oracle true Network | {oracle['exact_top1_hits']}/{summary['n_test_samples']} ({oracle['exact_top1_accuracy']:.1%}) | {oracle['exact_top3_hits']}/{summary['n_test_samples']} ({oracle['exact_top3_accuracy']:.1%}) | {oracle['group_top1_hits']}/{summary['n_test_samples']} ({oracle['group_top1_accuracy']:.1%}) | {oracle['group_top3_hits']}/{summary['n_test_samples']} ({oracle['group_top3_accuracy']:.1%}) |

## 输出粒度

- Top3 beam 路线中，Top1 输出被标记为 `low_resolution` 的样本数为 `{beam['low_resolution_predictions']}/{summary['n_test_samples']}`（{beam['low_resolution_fraction']:.1%}）。
- `low_resolution` 输出的精确 Region Top1 命中率为 `{beam['low_resolution_exact_top1']:.1%}`；`high_resolution` 输出数为 `{beam['high_resolution_predictions']}`，精确 Top1 命中率为 `{beam['high_resolution_exact_top1']:.1%}`。
- 全训练集生产注释模型包含 `{summary['production_model']['n_low_resolution_regions']}/{summary['production_model']['n_region_network_entries']}` 个低分辨率 Region-Network 条目与 `{summary['production_model']['n_merged_groups']}` 个合并候选组。

## 判定

{summary['decision']}
"""
    (outdir / "resolution_tier_report_cn.md").write_text(text, encoding="utf-8")


def build_production_model(
    values: np.ndarray,
    ann: pd.DataFrame,
    region_indices: dict[str, np.ndarray],
    local_top_n_genes: int,
    args: argparse.Namespace,
) -> dict[str, Any]:
    entries: dict[str, dict[str, Any]] = {}
    group_records: list[dict[str, Any]] = []
    for network in sorted(ann["endpoint_label"].unique().tolist()):
        regions = sorted(ann.loc[ann["endpoint_label"] == network, "region_id"].unique().tolist())
        unique_regions = [
            region for region in regions
            if ann.loc[ann["region_id"] == region, "endpoint_label"].nunique() == 1
        ]
        cross_regions = sorted(set(regions) - set(unique_regions))
        if len(unique_regions) >= 2:
            training = {region: region_indices[region] for region in unique_regions}
            rows, _ = select_group_discriminative_genes(values, unique_regions, training, local_top_n_genes)
            annotations, audit = build_resolution_groups(
                values, unique_regions, training, {region: network for region in unique_regions}, rows,
                args.min_resolution_samples, args.min_merge_samples, args.min_pair_errors,
                args.min_confusion_rate, args.similarity_threshold, args.merge_similarity_threshold,
                args.max_group_size,
            )
            for region, value in annotations.items():
                entries[f"{network}||{region}"] = value
            if not audit.empty:
                grouped = {value["resolution_group"]: value["group_members"] for value in annotations.values()}
                for group_id, members in grouped.items():
                    if len(members) > 1:
                        group_records.append({"network_id": network, "group_id": group_id, "members": members})
        elif len(unique_regions) == 1:
            region = unique_regions[0]
            entries[f"{network}||{region}"] = {
                "region_id": region,
                "network_id": network,
                "resolution_group": region,
                "group_members": [region],
                "resolution_tier": "low_resolution",
                "resolution_reasons": ["single_region_network;exact_resolution_not_testable"],
                "training_samples": int(len(region_indices[region])),
                "nearest_centroid_corr": None,
            }
        for region in cross_regions:
            entries[f"{network}||{region}"] = {
                "region_id": region,
                "network_id": network,
                "resolution_group": region,
                "group_members": [region],
                "resolution_tier": "low_resolution",
                "resolution_reasons": ["cross_network_region"],
                "training_samples": int((ann["region_id"].eq(region) & ann["endpoint_label"].eq(network)).sum()),
                "nearest_centroid_corr": None,
            }
    unique_groups = {(row["network_id"], row["group_id"]) for row in group_records}
    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "training_design": "full Bo2023 training-set model for output annotation; validation metrics use fold-local models only",
        "parameters": {
            "local_top_n_genes": int(local_top_n_genes),
            "min_resolution_samples": int(args.min_resolution_samples),
            "min_merge_samples": int(args.min_merge_samples),
            "min_pair_errors": int(args.min_pair_errors),
            "min_confusion_rate": float(args.min_confusion_rate),
            "similarity_threshold": float(args.similarity_threshold),
            "merge_similarity_threshold": float(args.merge_similarity_threshold),
            "max_group_size": int(args.max_group_size),
        },
        "entries": entries,
        "groups": group_records,
        "summary": {
            "n_region_network_entries": int(len(entries)),
            "n_low_resolution_regions": int(sum(v["resolution_tier"] == "low_resolution" for v in entries.values())),
            "n_merged_groups": int(len(unique_groups)),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict nested LOSO validation of fold-local resolvable Region groups.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--network-detail", type=Path, default=DEFAULT_NETWORK_DETAIL)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.20)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=4)
    parser.add_argument("--max-samples", type=int, default=None, help="Optional deterministic prefix for smoke tests only.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--model-out", type=Path, default=DEFAULT_MODEL_OUT)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.model_out.parent.mkdir(parents=True, exist_ok=True)
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

    region_counts = ann.groupby("region_id")["sample_id"].size()
    singletons = set(region_counts[region_counts < 2].index)
    selected = ann[~ann["region_id"].isin(singletons)].copy()
    selected["sort_order"] = selected["sample_id"].map(sample_pos)
    selected = selected.sort_values("sort_order").drop(columns="sort_order").reset_index(drop=True)
    if args.max_samples is not None:
        selected = selected.head(max(1, int(args.max_samples))).copy()

    route_rows: dict[str, list[dict[str, Any]]] = {BEAM_ROUTE: [], ORACLE_ROUTE: []}
    audit_rows: list[pd.DataFrame] = []
    fold_rows: list[dict[str, Any]] = []
    for fold, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth_region = str(row.region_id)
        truth_network = str(row.endpoint_label)
        heldout_idx = sample_pos[sample_id]
        sample = values[:, heldout_idx]
        training_ann = ann[ann["sample_id"] != sample_id].copy()
        network_top = [str(network_detail.loc[sample_id, f"pred_top{i}"]) for i in [1, 2, 3]]
        fold_info: dict[str, Any] = {"fold": fold, "sample_id": sample_id, "truth_region": truth_region}

        reference = reference_all.copy()
        truth_train = region_indices[truth_region][region_indices[truth_region] != heldout_idx]
        reference[:, region_pos[truth_region]] = values[:, truth_train].mean(axis=1, dtype=np.float64)
        for scope, network_values, route in [
            ("beam", network_top, BEAM_ROUTE),
            ("oracle", [truth_network], ORACLE_ROUTE),
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
            assignment = region_network_assignment(training_ann, candidates)
            annotations, audit = build_resolution_groups(
                values, candidates, training, assignment, local_rows,
                args.min_resolution_samples, args.min_merge_samples, args.min_pair_errors,
                args.min_confusion_rate, args.similarity_threshold, args.merge_similarity_threshold,
                args.max_group_size,
            )
            scores = correlation_scores(reference, sample, local_rows)
            indices = np.asarray([region_pos[region] for region in candidates], dtype=int)
            ranked = rank_candidates(scores, regions, indices)
            detail = score_route(
                route, sample_id, truth_region, truth_network, network_values, ranked,
                annotations, len(regions)
            )
            route_rows[route].append(detail)
            if not audit.empty:
                audit.insert(0, "fold", fold)
                audit.insert(1, "sample_id", sample_id)
                audit.insert(2, "scope", scope)
                audit_rows.append(audit)
            groups = {item["resolution_group"] for item in annotations.values()}
            merged = {group for group in groups if any(len(v["group_members"]) > 1 and v["resolution_group"] == group for v in annotations.values())}
            fold_info[f"{scope}_n_candidate_regions"] = int(len(candidates))
            fold_info[f"{scope}_n_resolution_groups"] = int(len(groups))
            fold_info[f"{scope}_n_merged_groups"] = int(len(merged))
            fold_info[f"{scope}_truth_resolution_tier"] = (
                annotations[truth_region]["resolution_tier"] if truth_region in annotations else "outside_network_beam"
            )
        fold_rows.append(fold_info)

    details = {route: pd.DataFrame(rows) for route, rows in route_rows.items()}
    metrics = {route: summarize(frame) for route, frame in details.items()}
    production_model = build_production_model(values, ann, region_indices, args.local_top_n_genes, args)
    args.model_out.write_text(json.dumps(production_model, ensure_ascii=False, indent=2), encoding="utf-8")
    decision = (
        "正式输出保留 Network 主结论与精确 Region 候选排名；对训练侧不可分或证据不足的 "
        "Region 输出 low_resolution 警告及可解释候选组，转入人工复核，而不强制声明精确 Top1。"
    )
    summary: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "strict outer LOSO with fold-local inner-confusion resolution groups and fold-local local genes",
        "n_test_samples": int(len(selected)),
        "n_singleton_samples_excluded": int(ann["region_id"].isin(singletons).sum()),
        "local_top_n_genes": int(args.local_top_n_genes),
        "min_resolution_samples": int(args.min_resolution_samples),
        "min_merge_samples": int(args.min_merge_samples),
        "min_pair_errors": int(args.min_pair_errors),
        "min_confusion_rate": float(args.min_confusion_rate),
        "similarity_threshold": float(args.similarity_threshold),
        "merge_similarity_threshold": float(args.merge_similarity_threshold),
        "max_group_size": int(args.max_group_size),
        "routes": metrics,
        "production_model": production_model["summary"],
        "production_model_path": str(args.model_out),
        "decision": decision,
    }
    pd.DataFrame([{"route": route, **metric} for route, metric in metrics.items()]).to_csv(
        args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig"
    )
    for route, detail in details.items():
        detail.to_csv(args.outdir / f"{route}_detail.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_rows).to_csv(args.outdir / "fold_resolution_summary.csv", index=False, encoding="utf-8-sig")
    if audit_rows:
        pd.concat(audit_rows, ignore_index=True).to_csv(
            args.outdir / "fold_pair_merge_audit.csv", index=False, encoding="utf-8-sig"
        )
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, metrics)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    print(f"Production annotation model written to: {args.model_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
