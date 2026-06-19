#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import (  # noqa: E402
    apply_projector,
    compute_logcpm,
    load_projector_npz,
    map_index_to_symbols,
    read_bo2023_gene_matrix,
    read_gene_map,
    write_json,
)
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_SAMPLE_INFO, DEFAULT_VSD  # noqa: E402
from scripts.run_ahba_human_rnaseq_external_validation import (  # noqa: E402
    hit_any,
    mapping_for_ahba_label,
    split_pipe,
    summarize_bool,
)
from scripts.run_ahba_projected_vsd_external_validation import load_ahba_counts, read_bo_metadata  # noqa: E402
from scripts.run_bo2023_loso_validation import correlation_scores  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP  # noqa: E402
from scripts.run_bo2023_projected_vsd_formal_lomo import EXACT_ROUTE, GROUP_ROUTE  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups, distinct_ranked_groups  # noqa: E402
from scripts.run_bo2023_leave_one_monkey_out_validation import build_region_reference, zscore  # noqa: E402


DEFAULT_ZIP_DIR = ROOT / "data" / "ahba_human_rnaseq" / "raw_zips"
DEFAULT_PROJECTOR = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "bo2023_reference_projector_linear_full.npz"
)
DEFAULT_OUTDIR = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "ahba_external_formal_three_tier"
)


def build_training(metadata: pd.DataFrame, labels_col: str, sample_pos: dict[str, int]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for label, rows in metadata.groupby(labels_col):
        idx = np.asarray([sample_pos[str(sample)] for sample in rows["sample_id"].astype(str)], dtype=int)
        if len(idx):
            out[str(label)] = idx
    return out


def centroid_reference(values: np.ndarray, labels: list[str], training: dict[str, np.ndarray]) -> np.ndarray:
    return np.column_stack([values[:, training[label]].mean(axis=1, dtype=np.float64) for label in labels]).astype(
        np.float32
    )


def pad(values: list[str], k: int) -> list[str]:
    return values[:k] + [""] * max(0, k - len(values))


def candidate_regions_from_networks(metadata: pd.DataFrame, network_top: list[str], region_training: dict[str, np.ndarray]) -> list[str]:
    return sorted(
        region
        for region in metadata.loc[metadata["network_id"].isin(network_top), "region_id"].dropna().astype(str).unique()
        if region in region_training
    )


def group_hits(
    ranked_regions: list[str],
    annotations: dict[str, dict[str, Any]],
    allowed_regions: list[str],
) -> tuple[bool | None, bool | None, list[str], list[str]]:
    if not allowed_regions:
        return None, None, [], []
    allowed_groups = sorted(
        {
            str(annotations[region]["resolution_group"])
            for region in allowed_regions
            if region in annotations
        }
    )
    if not allowed_groups:
        return False, False, [], []
    ranked_groups = distinct_ranked_groups(ranked_regions, annotations)
    return (
        hit_any(ranked_groups[:1], allowed_groups),
        hit_any(ranked_groups[:3], allowed_groups),
        ranked_groups,
        allowed_groups,
    )


def summarize_detail(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for route, frame in detail.groupby("route", sort=True):
        supported = frame[frame["supported_for_accuracy"].astype(bool)]
        exact = supported[supported["allowed_bo2023_regions"].fillna("").astype(str).str.len() > 0]
        rows.append(
            {
                "route": route,
                "n_samples_total": int(len(frame)),
                "n_samples_supported_for_accuracy": int(len(supported)),
                "n_samples_exact_region_evaluable": int(len(exact)),
                "network_top1_accuracy_coarse": summarize_bool(supported["network_top1_hit"]),
                "network_top3_accuracy_coarse": summarize_bool(supported["network_top3_hit"]),
                "group_top1_accuracy_exact_mapped": summarize_bool(exact["group_top1_hit"]),
                "group_top3_accuracy_exact_mapped": summarize_bool(exact["group_top3_hit"]),
                "region_top1_accuracy_exact_mapped": summarize_bool(exact["region_top1_exact_hit"]),
                "region_top3_accuracy_exact_mapped": summarize_bool(exact["region_top3_exact_hit"]),
                "mean_candidate_regions": float(supported["n_candidate_regions"].mean()),
                "low_resolution_top1_fraction": float((supported["pred_top1_resolution_tier"] == "low_resolution").mean()),
            }
        )
    return pd.DataFrame(rows)


def summarize_special_labels(detail: pd.DataFrame) -> pd.DataFrame:
    labels = {"Str::Caudate", "Str::Putamen", "Ins::Insula"}
    frame = detail[detail["public_label"].isin(labels)].copy()
    rows = []
    for (route, label), group in frame.groupby(["route", "public_label"], sort=True):
        rows.append(
            {
                "route": route,
                "public_label": label,
                "n": int(len(group)),
                "network_top1_accuracy_coarse": summarize_bool(group["network_top1_hit"]),
                "network_top3_accuracy_coarse": summarize_bool(group["network_top3_hit"]),
                "group_top1_accuracy_exact_mapped": summarize_bool(group["group_top1_hit"]),
                "group_top3_accuracy_exact_mapped": summarize_bool(group["group_top3_hit"]),
                "region_top1_accuracy_exact_mapped": summarize_bool(group["region_top1_exact_hit"]),
                "region_top3_accuracy_exact_mapped": summarize_bool(group["region_top3_exact_hit"]),
                "top1_predictions": " | ".join(group["region_top1"].astype(str).value_counts().head(5).index.tolist()),
            }
        )
    return pd.DataFrame(rows)


def plot_formal_metrics(metrics: pd.DataFrame, out_path: Path) -> None:
    cols = [
        "network_top1_accuracy_coarse",
        "network_top3_accuracy_coarse",
        "group_top3_accuracy_exact_mapped",
        "region_top1_accuracy_exact_mapped",
        "region_top3_accuracy_exact_mapped",
    ]
    labels = ["Network Top1", "Network Top3", "Group Top3", "Exact Top1", "Exact Top3"]
    routes = metrics["route"].astype(str).tolist()
    values = metrics.set_index("route").loc[routes, cols].to_numpy(dtype=float)
    x = np.arange(len(cols))
    width = 0.36
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    colors = ["#3B6EA8", "#C76F2D"]
    for idx, route in enumerate(routes):
        offset = (idx - (len(routes) - 1) / 2) * width
        bars = ax.bar(x + offset, values[idx], width=width, label=route, color=colors[idx % len(colors)])
        ax.bar_label(bars, labels=[f"{v:.2f}" for v in values[idx]], padding=3, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Accuracy")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title("AHBA formal three-tier external validation")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="AHBA projected-VSD formal three-tier external validation.")
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIP_DIR)
    parser.add_argument("--bo-counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--bo-vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--projector", type=Path, default=DEFAULT_PROJECTOR)
    parser.add_argument("--global-top-n-genes", type=int, default=200)
    parser.add_argument("--network-gene-pool-size", type=int, default=1000)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--exact-fusion-weight", type=float, default=0.25)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--group-min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.15)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=8)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ahba_metadata, ahba_counts = load_ahba_counts(args.zip_dir)
    ahba_logcpm = compute_logcpm(ahba_counts)
    ahba_projected = apply_projector(load_projector_npz(args.projector), ahba_logcpm)

    gene_map = read_gene_map(args.gene_map)
    bo_counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.bo_counts, dtype="float32"), gene_map)
    bo_vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.bo_vsd, dtype="float32"), gene_map)
    common = sorted(set(bo_counts.index.astype(str)) & set(bo_vsd.index.astype(str)))
    bo_counts = bo_counts.loc[common]
    bo_vsd = bo_vsd.loc[common]
    bo_logcpm = compute_logcpm(bo_counts)

    bo_metadata = read_bo_metadata(args.sample_info, args.sample_sheet, args.region_col, args.network_col)
    samples = [sample for sample in bo_counts.columns.astype(str) if sample in set(bo_metadata["sample_id"])]
    bo_metadata = bo_metadata[bo_metadata["sample_id"].isin(samples)].copy()
    bo_counts = bo_counts.loc[:, samples]
    bo_vsd = bo_vsd.loc[:, samples]
    bo_logcpm = bo_logcpm.loc[:, samples]

    sample_pos = {sample: idx for idx, sample in enumerate(samples)}
    region_training = build_training(bo_metadata, "region_id", sample_pos)
    network_training = build_training(bo_metadata, "network_id", sample_pos)
    regions = sorted(region_training)
    networks = sorted(network_training)
    region_network = {
        region: (
            values[0]
            if len(values := sorted(bo_metadata.loc[bo_metadata["region_id"].eq(region), "network_id"].astype(str).unique())) == 1
            else None
        )
        for region in regions
    }

    route_spaces = {
        "logcpm_baseline": {
            "network_train": bo_logcpm.loc[common].to_numpy(dtype=np.float32),
            "network_query": ahba_logcpm.reindex(common).fillna(0.0).astype("float32"),
            "exact_train": bo_logcpm.loc[common].to_numpy(dtype=np.float32),
            "exact_query": ahba_logcpm.reindex(common).fillna(0.0).astype("float32"),
        },
        "projected_vsd": {
            "network_train": bo_vsd.loc[common].to_numpy(dtype=np.float32),
            "network_query": ahba_projected.reindex(common).fillna(0.0).astype("float32"),
            "exact_train": bo_vsd.loc[common].to_numpy(dtype=np.float32),
            "exact_query": ahba_projected.reindex(common).fillna(0.0).astype("float32"),
        },
        "hybrid_projected_network_logcpm_exact": {
            "network_train": bo_vsd.loc[common].to_numpy(dtype=np.float32),
            "network_query": ahba_projected.reindex(common).fillna(0.0).astype("float32"),
            "exact_train": bo_logcpm.loc[common].to_numpy(dtype=np.float32),
            "exact_query": ahba_logcpm.reindex(common).fillna(0.0).astype("float32"),
        },
    }

    detail_rows: list[dict[str, Any]] = []
    resolution_audits: list[pd.DataFrame] = []

    for route_name, route in route_spaces.items():
        network_train_values = route["network_train"]
        network_query_matrix = route["network_query"]
        exact_train_values = route["exact_train"]
        exact_query_matrix = route["exact_query"]
        network_reference = centroid_reference(network_train_values, networks, network_training)
        network_rows, _ = select_group_discriminative_genes(
            network_train_values,
            networks,
            network_training,
            args.network_gene_pool_size,
        )
        network_rows = network_rows[: args.global_top_n_genes]

        for meta in ahba_metadata.itertuples(index=False):
            sample_id = str(meta.ahba_sample_uid)
            if sample_id not in network_query_matrix.columns or sample_id not in exact_query_matrix.columns:
                continue
            label_map = mapping_for_ahba_label(str(meta.main_structure), str(meta.sub_structure))
            supported = bool(label_map["supported_for_accuracy"])
            allowed_networks = split_pipe(label_map["allowed_bo2023_networks"])
            allowed_regions = split_pipe(label_map["allowed_bo2023_regions"])
            network_sample = network_query_matrix[sample_id].to_numpy(dtype=np.float32)
            exact_sample = exact_query_matrix[sample_id].to_numpy(dtype=np.float32)

            network_scores = correlation_scores(network_reference, network_sample, network_rows)
            network_top = [networks[i] for i in np.argsort(network_scores)[::-1][:3].tolist()]
            candidates = candidate_regions_from_networks(bo_metadata, network_top, region_training)
            if len(candidates) < 1:
                continue

            candidate_training = {region: region_training[region] for region in candidates}
            if len(candidates) >= 2:
                local_rows, _ = select_group_discriminative_genes(
                    exact_train_values,
                    candidates,
                    candidate_training,
                    args.local_top_n_genes,
                )
            else:
                local_rows = np.arange(train_values.shape[0], dtype=int)
            local_rows = np.asarray(local_rows, dtype=int)
            annotations, audit = build_resolution_groups(
                exact_train_values,
                candidates,
                candidate_training,
                {region: region_network.get(region) for region in candidates},
                local_rows,
                args.min_resolution_samples,
                args.min_merge_samples,
                args.group_min_pair_errors,
                args.min_confusion_rate,
                args.similarity_threshold,
                args.merge_similarity_threshold,
                args.max_group_size,
            )
            if not audit.empty:
                audit = audit.copy()
                audit["route"] = route_name
                audit["sample_id"] = sample_id
                resolution_audits.append(audit)

            reference = build_region_reference(exact_train_values, candidates, candidate_training)
            rows50 = local_rows[: min(50, len(local_rows))]
            rows100 = local_rows[: min(100, len(local_rows))]
            scores50 = correlation_scores(reference, exact_sample, rows50)
            scores100 = correlation_scores(reference, exact_sample, rows100)
            fused = args.exact_fusion_weight * zscore(scores50) + (1.0 - args.exact_fusion_weight) * zscore(scores100)
            ranked_regions = [candidates[i] for i in np.argsort(fused)[::-1].tolist()]
            ranked_groups = distinct_ranked_groups(ranked_regions, annotations)
            group_hit1, group_hit3, _, allowed_groups = group_hits(ranked_regions, annotations, allowed_regions)
            padded_regions = pad(ranked_regions, 3)
            padded_groups = pad(ranked_groups, 3)
            pred_annotation = annotations[padded_regions[0]]
            detail_rows.append(
                {
                    "route": route_name,
                    "route_stage": f"{GROUP_ROUTE} -> {EXACT_ROUTE}",
                    "network_space": "projected_vsd" if route_name.startswith("hybrid") else route_name,
                    "exact_space": "logcpm_baseline" if route_name.startswith("hybrid") else route_name,
                    "sample_id": sample_id,
                    "ahba_donor": str(meta.ahba_donor),
                    "public_label": label_map["public_label"],
                    "public_main_structure": label_map["public_main_structure"],
                    "public_sub_structure": label_map["public_sub_structure"],
                    "supported_for_accuracy": supported,
                    "accuracy_level": label_map["accuracy_level"],
                    "allowed_bo2023_networks": label_map["allowed_bo2023_networks"],
                    "allowed_bo2023_regions": label_map["allowed_bo2023_regions"],
                    "allowed_resolution_groups": " | ".join(allowed_groups),
                    "network_top1": network_top[0],
                    "network_top2": network_top[1] if len(network_top) > 1 else "",
                    "network_top3": network_top[2] if len(network_top) > 2 else "",
                    "network_top1_hit": hit_any(network_top[:1], allowed_networks) if supported else None,
                    "network_top3_hit": hit_any(network_top, allowed_networks) if supported else None,
                    "region_top1": padded_regions[0],
                    "region_top2": padded_regions[1],
                    "region_top3": padded_regions[2],
                    "region_top1_exact_hit": hit_any(padded_regions[:1], allowed_regions) if supported and allowed_regions else None,
                    "region_top3_exact_hit": hit_any(padded_regions, allowed_regions) if supported and allowed_regions else None,
                    "group_top1": padded_groups[0],
                    "group_top2": padded_groups[1],
                    "group_top3": padded_groups[2],
                    "group_top1_hit": group_hit1 if supported and allowed_regions else None,
                    "group_top3_hit": group_hit3 if supported and allowed_regions else None,
                    "pred_top1_resolution_tier": pred_annotation["resolution_tier"],
                    "pred_top1_resolution_group": pred_annotation["resolution_group"],
                    "pred_top1_group_members": " | ".join(pred_annotation["group_members"]),
                    "n_candidate_regions": int(len(candidates)),
                    "n_local_genes": int(len(local_rows)),
                    "network_overlap_genes": int(len(network_rows)),
                    "exact_overlap_genes_top100": int(len(rows100)),
                }
            )

    detail = pd.DataFrame(detail_rows)
    metrics = summarize_detail(detail)
    special = summarize_special_labels(detail)
    detail.to_csv(args.outdir / "ahba_formal_three_tier_sample_detail.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(args.outdir / "ahba_formal_three_tier_metrics.csv", index=False, encoding="utf-8-sig")
    special.to_csv(args.outdir / "ahba_formal_three_tier_special_labels.csv", index=False, encoding="utf-8-sig")
    if resolution_audits:
        pd.concat(resolution_audits, ignore_index=True).to_csv(
            args.outdir / "ahba_formal_three_tier_resolution_audit.csv",
            index=False,
            encoding="utf-8-sig",
        )
    plot_formal_metrics(metrics, args.outdir / "ahba_formal_three_tier_accuracy.png")

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "AHBA human RNA-seq raw counts",
        "route": "Network Top3 beam -> resolution group -> local exact rerank",
        "metrics": metrics.to_dict(orient="records"),
        "special_labels": special.to_dict(orient="records"),
        "caution": "Cross-species external validation; exact metrics only for AHBA labels with stable Bo2023 region mappings.",
    }
    write_json(args.outdir / "ahba_formal_three_tier_summary.json", summary)

    lines = [
        "# AHBA formal three-tier external validation",
        "",
        "Route: Network Top3 beam -> resolution group -> local exact rerank.",
        "",
        "| route | network Top1 | network Top3 | group Top1 | group Top3 | exact Top1 | exact Top3 | exact n |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics.to_dict(orient="records"):
        lines.append(
            "| {route} | {network_top1_accuracy_coarse:.6f} | {network_top3_accuracy_coarse:.6f} | "
            "{group_top1_accuracy_exact_mapped:.6f} | {group_top3_accuracy_exact_mapped:.6f} | "
            "{region_top1_accuracy_exact_mapped:.6f} | {region_top3_accuracy_exact_mapped:.6f} | "
            "{n_samples_exact_region_evaluable} |".format(**row)
        )
    lines.extend(["", "Special labels: Caudate, Putamen, Insula are written to `ahba_formal_three_tier_special_labels.csv`."])
    (args.outdir / "ahba_formal_three_tier_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(args.outdir)
    print(metrics.to_string(index=False))
    print(special.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
