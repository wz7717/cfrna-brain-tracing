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
from scripts.run_ahba_human_rnaseq_external_validation import summarize_bool  # noqa: E402
from scripts.run_ahba_projected_vsd_external_validation import read_bo_metadata  # noqa: E402
from scripts.run_bo2023_loso_validation import correlation_scores  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups, distinct_ranked_groups  # noqa: E402
from scripts.run_bo2023_leave_one_monkey_out_validation import build_region_reference, zscore  # noqa: E402
from scripts.score_tcga_gbm_lgg_sample_tracing_with_mri_labels import NETWORK_TO_BROAD, sample_to_patient  # noqa: E402


DEFAULT_COUNTS_MATRIX = (
    ROOT
    / "data"
    / "tcga_brain_tumor_expression"
    / "tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv"
)
DEFAULT_MANIFEST = (
    ROOT
    / "data"
    / "tcga_brain_tumor_expression"
    / "tcga_gbm_lgg_gdc_star_counts_manifest.csv"
)
DEFAULT_LABELS = (
    ROOT
    / "results"
    / "brats_tcga_lgg_65_mri_truth_corrected_20260612"
    / "corrected_direct_overlap_mri_truth.csv"
)
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
    / "tcga_labeled_hybrid_formal_external"
)


def split_candidates(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return []
    out: list[str] = []
    for chunk in text.replace(";", "|").split("|"):
        item = chunk.strip()
        if not item:
            continue
        out.append(item.split(":", 1)[0].strip())
    return list(dict.fromkeys(out))


def norm(value: Any) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def hit(predicted: list[str], truth: list[str]) -> bool | None:
    truth_norm = {norm(x) for x in truth if norm(x)}
    if not truth_norm:
        return None
    return any(norm(x) in truth_norm for x in predicted)


def build_training(metadata: pd.DataFrame, label_col: str, sample_pos: dict[str, int]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for label, rows in metadata.groupby(label_col):
        idx = np.asarray([sample_pos[str(sample)] for sample in rows["sample_id"].astype(str)], dtype=int)
        if len(idx):
            out[str(label)] = idx
    return out


def centroid_reference(values: np.ndarray, labels: list[str], training: dict[str, np.ndarray]) -> np.ndarray:
    return np.column_stack([values[:, training[label]].mean(axis=1, dtype=np.float64) for label in labels]).astype(
        np.float32
    )


def load_tcga_counts(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = matrix.columns.astype(str).str.strip()
    matrix = matrix.groupby(matrix.index).mean()
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")


def candidate_regions(metadata: pd.DataFrame, network_top: list[str], region_training: dict[str, np.ndarray]) -> list[str]:
    return sorted(
        region
        for region in metadata.loc[metadata["network_id"].isin(network_top), "region_id"].dropna().astype(str).unique()
        if region in region_training
    )


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for route, frame in detail.groupby("route", sort=True):
        rows.append(
            {
                "route": route,
                "n_samples": int(len(frame)),
                "n_patients": int(frame["patient_barcode"].nunique()),
                "network_top1_accuracy": summarize_bool(frame["network_top1_hit"]),
                "network_top3_accuracy": summarize_bool(frame["network_top3_hit"]),
                "lobe_top1_accuracy": summarize_bool(frame["lobe_top1_hit"]),
                "lobe_top3_accuracy": summarize_bool(frame["lobe_top3_hit"]),
                "broad_top1_accuracy": summarize_bool(frame["broad_top1_hit"]),
                "broad_top3_accuracy": summarize_bool(frame["broad_top3_hit"]),
                "mean_candidate_regions": float(frame["n_candidate_regions"].mean()),
                "low_resolution_top1_fraction": float((frame["pred_top1_resolution_tier"] == "low_resolution").mean()),
            }
        )
    return pd.DataFrame(rows)


def plot_metrics(metrics: pd.DataFrame, out_path: Path) -> None:
    cols = ["network_top1_accuracy", "network_top3_accuracy", "lobe_top1_accuracy", "lobe_top3_accuracy"]
    labels = ["Network Top1", "Network Top3", "Lobe Top1", "Lobe Top3"]
    routes = metrics["route"].astype(str).tolist()
    values = metrics.set_index("route").loc[routes, cols].to_numpy(dtype=float)
    x = np.arange(len(cols))
    width = 0.8 / max(len(routes), 1)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    colors = ["#3B6EA8", "#C76F2D", "#2F7D59"]
    for idx, route in enumerate(routes):
        offset = (idx - (len(routes) - 1) / 2) * width
        bars = ax.bar(x + offset, values[idx], width=width, label=route, color=colors[idx % len(colors)])
        ax.bar_label(bars, labels=[f"{v:.2f}" for v in values[idx]], padding=3, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Accuracy")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title("TCGA/BraTS MRI-labeled hybrid external validation")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Hybrid projected-network/logCPM-exact external validation on labeled TCGA/BraTS.")
    parser.add_argument("--tcga-counts", type=Path, default=DEFAULT_COUNTS_MATRIX)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--bo-counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--bo-vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--projector", type=Path, default=DEFAULT_PROJECTOR)
    parser.add_argument("--region-meta", type=Path, default=ROOT / "bo2023_bulk_atlas_buildkit" / "01_region_level_meta.csv")
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

    labels = pd.read_csv(args.labels)
    labels["patient_barcode"] = labels["patient_barcode"].astype(str).str.upper().str.slice(0, 12)
    manifest = pd.read_csv(args.manifest)
    sample_project = (
        manifest[["sample_submitter_id", "project_id", "case_submitter_id"]]
        .drop_duplicates("sample_submitter_id")
        .set_index("sample_submitter_id")
        .to_dict("index")
    )

    tcga_counts = load_tcga_counts(args.tcga_counts)
    sample_ids = [
        sample
        for sample in tcga_counts.columns.astype(str)
        if sample_to_patient(sample) in set(labels["patient_barcode"])
    ]
    tcga_counts = tcga_counts.loc[:, sample_ids]
    tcga_logcpm = compute_logcpm(tcga_counts)
    tcga_projected = apply_projector(load_projector_npz(args.projector), tcga_logcpm)

    gene_map = read_gene_map(args.gene_map)
    bo_counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.bo_counts, dtype="float32"), gene_map)
    bo_vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.bo_vsd, dtype="float32"), gene_map)
    common = sorted(set(bo_counts.index.astype(str)) & set(bo_vsd.index.astype(str)))
    bo_counts = bo_counts.loc[common]
    bo_vsd = bo_vsd.loc[common]
    bo_logcpm = compute_logcpm(bo_counts)

    bo_metadata = read_bo_metadata(args.sample_info, args.sample_sheet, args.region_col, args.network_col)
    bo_samples = [sample for sample in bo_counts.columns.astype(str) if sample in set(bo_metadata["sample_id"])]
    bo_metadata = bo_metadata[bo_metadata["sample_id"].isin(bo_samples)].copy()
    bo_logcpm = bo_logcpm.loc[common, bo_samples]
    bo_vsd = bo_vsd.loc[common, bo_samples]

    sample_pos = {sample: idx for idx, sample in enumerate(bo_samples)}
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
    region_lobe = dict(
        zip(
            pd.read_csv(args.region_meta, usecols=["Region", "Lobe"])["Region"].astype(str),
            pd.read_csv(args.region_meta, usecols=["Region", "Lobe"])["Lobe"].astype(str).str.lower(),
        )
    )

    route_spaces = {
        "logcpm_baseline": {
            "network_train": bo_logcpm.to_numpy(dtype=np.float32),
            "network_query": tcga_logcpm.reindex(common).fillna(0.0).astype("float32"),
            "exact_train": bo_logcpm.to_numpy(dtype=np.float32),
            "exact_query": tcga_logcpm.reindex(common).fillna(0.0).astype("float32"),
        },
        "projected_vsd": {
            "network_train": bo_vsd.to_numpy(dtype=np.float32),
            "network_query": tcga_projected.reindex(common).fillna(0.0).astype("float32"),
            "exact_train": bo_vsd.to_numpy(dtype=np.float32),
            "exact_query": tcga_projected.reindex(common).fillna(0.0).astype("float32"),
        },
        "hybrid_projected_network_logcpm_exact": {
            "network_train": bo_vsd.to_numpy(dtype=np.float32),
            "network_query": tcga_projected.reindex(common).fillna(0.0).astype("float32"),
            "exact_train": bo_logcpm.to_numpy(dtype=np.float32),
            "exact_query": tcga_logcpm.reindex(common).fillna(0.0).astype("float32"),
        },
    }

    rows: list[dict[str, Any]] = []
    for route_name, route in route_spaces.items():
        network_ref = centroid_reference(route["network_train"], networks, network_training)
        network_rows, _ = select_group_discriminative_genes(
            route["network_train"],
            networks,
            network_training,
            args.network_gene_pool_size,
        )
        network_rows = network_rows[: args.global_top_n_genes]
        for sample_id in sample_ids:
            patient = sample_to_patient(sample_id)
            truth = labels[labels["patient_barcode"].eq(patient)].iloc[0].to_dict()
            network_truth = split_candidates(truth.get("corrected_network_candidates", ""))
            lobe_truth = split_candidates(truth.get("corrected_lobe_candidates", ""))
            broad_truth = split_candidates(truth.get("corrected_broad_candidates", ""))
            net_sample = route["network_query"][sample_id].to_numpy(dtype=np.float32)
            exact_sample = route["exact_query"][sample_id].to_numpy(dtype=np.float32)
            net_scores = correlation_scores(network_ref, net_sample, network_rows)
            network_top = [networks[i] for i in np.argsort(net_scores)[::-1][:3].tolist()]
            candidates = candidate_regions(bo_metadata, network_top, region_training)
            if len(candidates) < 1:
                continue
            candidate_training = {region: region_training[region] for region in candidates}
            if len(candidates) >= 2:
                local_rows, _ = select_group_discriminative_genes(
                    route["exact_train"],
                    candidates,
                    candidate_training,
                    args.local_top_n_genes,
                )
            else:
                local_rows = np.arange(route["exact_train"].shape[0], dtype=int)
            local_rows = np.asarray(local_rows, dtype=int)
            annotations, _ = build_resolution_groups(
                route["exact_train"],
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
            region_ref = build_region_reference(route["exact_train"], candidates, candidate_training)
            rows50 = local_rows[: min(50, len(local_rows))]
            rows100 = local_rows[: min(100, len(local_rows))]
            fused = args.exact_fusion_weight * zscore(correlation_scores(region_ref, exact_sample, rows50)) + (
                1.0 - args.exact_fusion_weight
            ) * zscore(correlation_scores(region_ref, exact_sample, rows100))
            ranked_regions = [candidates[i] for i in np.argsort(fused)[::-1].tolist()]
            ranked_groups = distinct_ranked_groups(ranked_regions, annotations)
            region_top3 = ranked_regions[:3]
            lobe_top3 = [region_lobe.get(region, "") for region in region_top3]
            broad_top3 = [NETWORK_TO_BROAD.get(network, network) for network in network_top]
            pred_top1 = annotations[region_top3[0]]
            rows.append(
                {
                    "route": route_name,
                    "sample_id": sample_id,
                    "patient_barcode": patient,
                    "project_id": sample_project.get(sample_id, {}).get("project_id", ""),
                    "segmentation_source": truth.get("segmentation_source", ""),
                    "truth_network_candidates": " | ".join(network_truth),
                    "truth_lobe_candidates": " | ".join(lobe_truth),
                    "truth_broad_candidates": " | ".join(broad_truth),
                    "network_top1": network_top[0],
                    "network_top3": " | ".join(network_top),
                    "region_top1": region_top3[0],
                    "region_top3": " | ".join(region_top3),
                    "region_group_top1": ranked_groups[0] if ranked_groups else "",
                    "region_group_top3": " | ".join(ranked_groups[:3]),
                    "pred_lobe_top1": lobe_top3[0] if lobe_top3 else "",
                    "pred_lobe_top3": " | ".join(lobe_top3),
                    "pred_broad_top1": broad_top3[0] if broad_top3 else "",
                    "pred_broad_top3": " | ".join(broad_top3),
                    "network_top1_hit": hit(network_top[:1], network_truth),
                    "network_top3_hit": hit(network_top[:3], network_truth),
                    "lobe_top1_hit": hit(lobe_top3[:1], lobe_truth),
                    "lobe_top3_hit": hit(lobe_top3[:3], lobe_truth),
                    "broad_top1_hit": hit(broad_top3[:1], broad_truth),
                    "broad_top3_hit": hit(broad_top3[:3], broad_truth),
                    "pred_top1_resolution_tier": pred_top1["resolution_tier"],
                    "pred_top1_resolution_group": pred_top1["resolution_group"],
                    "n_candidate_regions": int(len(candidates)),
                    "n_local_genes": int(len(local_rows)),
                }
            )

    detail = pd.DataFrame(rows)
    metrics = summarize(detail)
    detail.to_csv(args.outdir / "tcga_labeled_hybrid_formal_sample_detail.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(args.outdir / "tcga_labeled_hybrid_formal_metrics.csv", index=False, encoding="utf-8-sig")
    plot_metrics(metrics, args.outdir / "tcga_labeled_hybrid_formal_accuracy.png")
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "TCGA/BraTS glioma RNA-seq with corrected MRI-derived labels",
        "n_labeled_expression_samples": int(len(sample_ids)),
        "n_labeled_patients": int(len({sample_to_patient(x) for x in sample_ids})),
        "metrics": metrics.to_dict(orient="records"),
        "label_file": str(args.labels),
        "caution": "MRI truth exact labels are human atlas labels, not Bo2023 macaque region IDs; exact Bo2023 region accuracy is not reported.",
    }
    write_json(args.outdir / "tcga_labeled_hybrid_formal_summary.json", summary)
    lines = [
        "# TCGA/BraTS labeled hybrid formal external validation",
        "",
        "Route tested: projected VSD network Top3 beam + logCPM resolution/local exact rerank.",
        "",
        "| route | network Top1 | network Top3 | lobe Top1 | lobe Top3 | broad Top1 | broad Top3 | n |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics.to_dict(orient="records"):
        lines.append(
            "| {route} | {network_top1_accuracy:.6f} | {network_top3_accuracy:.6f} | "
            "{lobe_top1_accuracy:.6f} | {lobe_top3_accuracy:.6f} | "
            "{broad_top1_accuracy:.6f} | {broad_top3_accuracy:.6f} | {n_samples} |".format(**row)
        )
    lines.append("")
    lines.append("Exact Bo2023 region accuracy is not reported because MRI truth regions use human atlas labels.")
    (args.outdir / "tcga_labeled_hybrid_formal_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.outdir)
    print(metrics.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
