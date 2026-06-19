#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP  # noqa: E402


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
    / "ahba_external_projected_vsd"
)


def load_ahba_counts(zip_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata_rows: list[pd.DataFrame] = []
    matrix_parts: list[pd.DataFrame] = []
    for zip_path in sorted(zip_dir.glob("*.zip")):
        donor = zip_path.stem.replace("_rnaseq", "")
        with zipfile.ZipFile(zip_path) as zf:
            annot = pd.read_csv(zf.open("SampleAnnot.csv"))
            sample_ids = [f"{donor}|{sample}" for sample in annot["RNAseq_sample_name"].astype(str)]
            annot = annot.copy()
            annot["ahba_donor"] = donor
            annot["ahba_sample_uid"] = sample_ids
            counts = pd.read_csv(zf.open("RNAseqCounts.csv"), header=None)
        if counts.shape[1] - 1 != len(sample_ids):
            raise ValueError(f"{zip_path.name} counts/sample mismatch: {counts.shape[1] - 1} vs {len(sample_ids)}")
        counts = counts.rename(columns={0: "gene_symbol"})
        counts.columns = ["gene_symbol", *sample_ids]
        matrix_parts.append(counts.set_index("gene_symbol"))
        metadata_rows.append(annot)
    metadata = pd.concat(metadata_rows, ignore_index=True)
    matrix = pd.concat(matrix_parts, axis=1)
    matrix.index = matrix.index.astype(str).str.strip()
    matrix = matrix.groupby(matrix.index).mean()
    return metadata, matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float32")


def read_bo_metadata(path: Path, sheet: str, region_col: str, network_col: str) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", region_col, network_col])
    info = info.copy()
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region_id"] = info[region_col].astype(str).str.strip()
    info["network_id"] = info[network_col].astype(str).str.strip()
    return info.drop_duplicates("sample_id")


def build_centroids(matrix: pd.DataFrame, metadata: pd.DataFrame, label_col: str) -> pd.DataFrame:
    sample_ids = [sample for sample in matrix.columns.astype(str) if sample in set(metadata["sample_id"])]
    label_by_sample = metadata.set_index("sample_id").loc[sample_ids, label_col].astype(str)
    refs = []
    labels = []
    for label, label_samples in label_by_sample.groupby(label_by_sample).groups.items():
        cols = list(label_samples)
        if cols:
            refs.append(matrix.loc[:, cols].mean(axis=1))
            labels.append(str(label))
    return pd.concat(refs, axis=1, keys=labels).sort_index(axis=1)


def pearson_scores(reference: pd.DataFrame, sample: pd.Series, genes: list[str]) -> np.ndarray:
    ref = reference.loc[genes].to_numpy(dtype=np.float64)
    vec = sample.loc[genes].to_numpy(dtype=np.float64)
    ref0 = ref - ref.mean(axis=0, keepdims=True)
    vec0 = vec - vec.mean()
    denom = np.sqrt(np.square(ref0).sum(axis=0) * np.square(vec0).sum() + 1e-12)
    return np.nan_to_num((ref0 * vec0[:, None]).sum(axis=0) / denom, nan=0.0, posinf=0.0, neginf=0.0)


def rank_labels(reference: pd.DataFrame, sample: pd.Series, genes: list[str], k: int = 3) -> list[tuple[str, float]]:
    scores = pearson_scores(reference, sample, genes)
    order = np.argsort(-scores)[:k]
    labels = reference.columns.astype(str).to_numpy()
    return [(str(labels[i]), float(scores[i])) for i in order]


def evaluate_route(
    route_name: str,
    query_matrix: pd.DataFrame,
    network_reference: pd.DataFrame,
    region_reference: pd.DataFrame,
    ahba_metadata: pd.DataFrame,
    region_lobe: dict[str, str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    network_genes = sorted(set(query_matrix.index.astype(str)) & set(network_reference.index.astype(str)))
    region_genes = sorted(set(query_matrix.index.astype(str)) & set(region_reference.index.astype(str)))
    rows: list[dict[str, Any]] = []
    for meta in ahba_metadata.itertuples(index=False):
        sample_id = str(meta.ahba_sample_uid)
        if sample_id not in query_matrix.columns:
            continue
        label_map = mapping_for_ahba_label(str(meta.main_structure), str(meta.sub_structure))
        allowed_networks = split_pipe(label_map["allowed_bo2023_networks"])
        allowed_lobes = split_pipe(label_map["allowed_bo2023_lobes"])
        allowed_regions = split_pipe(label_map["allowed_bo2023_regions"])
        supported = bool(label_map["supported_for_accuracy"])
        sample = query_matrix[sample_id]
        network_top = rank_labels(network_reference, sample, network_genes, k=3)
        region_top = rank_labels(region_reference, sample, region_genes, k=3)
        network_labels = [label for label, _ in network_top]
        region_labels = [label for label, _ in region_top]
        region_lobes = [region_lobe.get(label, "") for label in region_labels]
        rows.append(
            {
                "route": route_name,
                "sample_id": sample_id,
                "ahba_donor": str(meta.ahba_donor),
                "public_label": label_map["public_label"],
                "public_main_structure": label_map["public_main_structure"],
                "public_sub_structure": label_map["public_sub_structure"],
                "supported_for_accuracy": supported,
                "accuracy_level": label_map["accuracy_level"],
                "allowed_bo2023_networks": label_map["allowed_bo2023_networks"],
                "allowed_bo2023_lobes": label_map["allowed_bo2023_lobes"],
                "allowed_bo2023_regions": label_map["allowed_bo2023_regions"],
                "network_top1": network_labels[0] if network_labels else "",
                "network_top2": network_labels[1] if len(network_labels) > 1 else "",
                "network_top3": network_labels[2] if len(network_labels) > 2 else "",
                "network_top1_score": network_top[0][1] if network_top else np.nan,
                "region_top1": region_labels[0] if region_labels else "",
                "region_top2": region_labels[1] if len(region_labels) > 1 else "",
                "region_top3": region_labels[2] if len(region_labels) > 2 else "",
                "region_top1_lobe": region_lobes[0] if region_lobes else "",
                "region_top1_score": region_top[0][1] if region_top else np.nan,
                "network_top1_hit": hit_any(network_labels[:1], allowed_networks) if supported else None,
                "network_top3_hit": hit_any(network_labels, allowed_networks) if supported else None,
                "region_lobe_top1_hit": hit_any(region_lobes[:1], allowed_lobes) if supported else None,
                "region_top1_exact_hit": hit_any(region_labels[:1], allowed_regions) if supported and allowed_regions else None,
                "region_top3_exact_hit": hit_any(region_labels, allowed_regions) if supported and allowed_regions else None,
                "network_overlap_genes": len(network_genes),
                "region_overlap_genes": len(region_genes),
            }
        )
    detail = pd.DataFrame(rows)
    supported_detail = detail[detail["supported_for_accuracy"].astype(bool)].copy()
    exact_detail = supported_detail[
        supported_detail["allowed_bo2023_regions"].fillna("").astype(str).str.len() > 0
    ].copy()
    metrics = {
        "route": route_name,
        "n_samples_total": int(len(detail)),
        "n_samples_supported_for_accuracy": int(len(supported_detail)),
        "n_samples_exact_region_evaluable": int(len(exact_detail)),
        "network_top1_accuracy_coarse": summarize_bool(supported_detail["network_top1_hit"]),
        "network_top3_accuracy_coarse": summarize_bool(supported_detail["network_top3_hit"]),
        "region_lobe_top1_accuracy_coarse": summarize_bool(supported_detail["region_lobe_top1_hit"]),
        "region_top1_accuracy_exact_mapped": summarize_bool(exact_detail["region_top1_exact_hit"]),
        "region_top3_accuracy_exact_mapped": summarize_bool(exact_detail["region_top3_exact_hit"]),
        "network_overlap_genes": int(len(network_genes)),
        "region_overlap_genes": int(len(region_genes)),
    }
    return detail, metrics


def summarize_by_label(detail_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (route, label), frame in detail_df.groupby(["route", "public_label"], sort=True):
        supported = frame[frame["supported_for_accuracy"].astype(bool)]
        exact = supported[supported["allowed_bo2023_regions"].fillna("").astype(str).str.len() > 0]
        rows.append(
            {
                "route": route,
                "public_label": label,
                "n": int(len(frame)),
                "n_supported": int(len(supported)),
                "n_exact_evaluable": int(len(exact)),
                "network_top1_accuracy_coarse": summarize_bool(supported["network_top1_hit"]),
                "network_top3_accuracy_coarse": summarize_bool(supported["network_top3_hit"]),
                "region_lobe_top1_accuracy_coarse": summarize_bool(supported["region_lobe_top1_hit"]),
                "region_top1_accuracy_exact_mapped": summarize_bool(exact["region_top1_exact_hit"]),
                "region_top3_accuracy_exact_mapped": summarize_bool(exact["region_top3_exact_hit"]),
            }
        )
    return pd.DataFrame(rows)


def plot_metric_comparison(metrics_df: pd.DataFrame, out_path: Path) -> None:
    metric_cols = [
        "network_top1_accuracy_coarse",
        "network_top3_accuracy_coarse",
        "region_lobe_top1_accuracy_coarse",
        "region_top1_accuracy_exact_mapped",
        "region_top3_accuracy_exact_mapped",
    ]
    labels = ["Network Top1", "Network Top3", "Lobe Top1", "Exact Top1", "Exact Top3"]
    routes = metrics_df["route"].astype(str).tolist()
    values = metrics_df.set_index("route").loc[routes, metric_cols].to_numpy(dtype=float)
    x = np.arange(len(metric_cols))
    width = 0.36
    fig, ax = plt.subplots(figsize=(10, 4.8))
    colors = ["#3B6EA8", "#C76F2D"]
    for idx, route in enumerate(routes):
        offset = (idx - (len(routes) - 1) / 2) * width
        bars = ax.bar(x + offset, values[idx], width=width, label=route, color=colors[idx % len(colors)])
        ax.bar_label(bars, labels=[f"{v:.2f}" for v in values[idx]], padding=3, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Accuracy")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title("AHBA external validation accuracy")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="AHBA raw-count external validation for Bo2023 projected VSD.")
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIP_DIR)
    parser.add_argument("--bo-counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--bo-vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--projector", type=Path, default=DEFAULT_PROJECTOR)
    parser.add_argument("--region-meta", type=Path, default=ROOT / "bo2023_bulk_atlas_buildkit" / "01_region_level_meta.csv")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    ahba_metadata, ahba_counts = load_ahba_counts(args.zip_dir)
    ahba_logcpm = compute_logcpm(ahba_counts)
    projected = apply_projector(load_projector_npz(args.projector), ahba_logcpm)

    gene_map = read_gene_map(args.gene_map)
    bo_counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.bo_counts, dtype="float32"), gene_map)
    bo_vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.bo_vsd, dtype="float32"), gene_map)
    common = sorted(set(bo_counts.index.astype(str)) & set(bo_vsd.index.astype(str)))
    bo_counts = bo_counts.loc[common]
    bo_vsd = bo_vsd.loc[common]
    bo_logcpm = compute_logcpm(bo_counts)

    bo_metadata = read_bo_metadata(args.sample_info, args.sample_sheet, args.region_col, args.network_col)
    region_meta = pd.read_csv(args.region_meta, usecols=["Region", "Lobe"])
    region_lobe = dict(zip(region_meta["Region"].astype(str), region_meta["Lobe"].astype(str)))

    route_refs = {
        "logcpm_baseline": {
            "query": ahba_logcpm,
            "network_reference": build_centroids(bo_logcpm, bo_metadata, "network_id"),
            "region_reference": build_centroids(bo_logcpm, bo_metadata, "region_id"),
        },
        "projected_vsd": {
            "query": projected,
            "network_reference": build_centroids(bo_vsd, bo_metadata, "network_id"),
            "region_reference": build_centroids(bo_vsd, bo_metadata, "region_id"),
        },
    }

    details: list[pd.DataFrame] = []
    metrics: list[dict[str, Any]] = []
    for route_name, route in route_refs.items():
        detail, route_metrics = evaluate_route(
            route_name,
            route["query"],
            route["network_reference"],
            route["region_reference"],
            ahba_metadata,
            region_lobe,
        )
        details.append(detail)
        metrics.append(route_metrics)

    detail_df = pd.concat(details, ignore_index=True)
    metrics_df = pd.DataFrame(metrics)
    label_df = summarize_by_label(detail_df)
    detail_df.to_csv(args.outdir / "ahba_projected_vsd_external_sample_detail.csv", index=False, encoding="utf-8-sig")
    metrics_df.to_csv(args.outdir / "ahba_projected_vsd_external_metrics.csv", index=False, encoding="utf-8-sig")
    label_df.to_csv(args.outdir / "ahba_projected_vsd_external_label_metrics.csv", index=False, encoding="utf-8-sig")
    plot_metric_comparison(metrics_df, args.outdir / "ahba_projected_vsd_external_accuracy.png")

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "AHBA human RNA-seq raw RNAseqCounts.csv",
        "n_ahba_samples": int(ahba_counts.shape[1]),
        "n_ahba_genes": int(ahba_counts.shape[0]),
        "routes": metrics,
        "caution": (
            "AHBA is human normal brain tissue; Bo2023 is macaque atlas. Treat accuracy as cross-species "
            "coarse anatomical transfer validation, with exact-region metrics only for stable mapped labels."
        ),
    }
    write_json(args.outdir / "ahba_projected_vsd_external_summary.json", summary)

    lines = [
        "# AHBA projected-VSD external validation",
        "",
        "Input: AHBA `RNAseqCounts.csv` raw counts from local zip files.",
        "",
        "| route | network Top1 | network Top3 | lobe Top1 | exact Top1 | exact Top3 | supported n | exact n |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics:
        lines.append(
            "| {route} | {network_top1_accuracy_coarse:.6f} | {network_top3_accuracy_coarse:.6f} | "
            "{region_lobe_top1_accuracy_coarse:.6f} | {region_top1_accuracy_exact_mapped:.6f} | "
            "{region_top3_accuracy_exact_mapped:.6f} | {n_samples_supported_for_accuracy} | "
            "{n_samples_exact_region_evaluable} |".format(**row)
        )
    lines.extend(
        [
            "",
            "Caution: this is cross-species external validation. Exact-region accuracy is only reported for AHBA labels with stable Bo2023 region mappings.",
        ]
    )
    (args.outdir / "ahba_projected_vsd_external_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.outdir)
    print(metrics_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
