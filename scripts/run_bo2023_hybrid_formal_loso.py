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

from core.reference_projection import (  # noqa: E402
    align_matrices,
    compute_logcpm,
    map_index_to_symbols,
    read_bo2023_gene_matrix,
    read_gene_map,
    write_json,
)
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_MODEL_GENES, DEFAULT_SAMPLE_INFO, DEFAULT_VSD, read_locked_model_genes  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP  # noqa: E402
from scripts.run_bo2023_projected_vsd_region_local_rerank import build_centroids, loo_project_rows  # noqa: E402
from scripts.run_bo2023_loso_validation import correlation_scores  # noqa: E402
from scripts.run_bo2023_leave_one_monkey_out_validation import build_region_reference, zscore  # noqa: E402
from scripts.run_bo2023_projected_vsd_formal_lomo import EXACT_ROUTE, GROUP_ROUTE  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups, score_route  # noqa: E402


DEFAULT_OUTDIR = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "formal_three_tier_loso_hybrid"
)


def read_metadata(path: Path, sheet: str, region_col: str, network_col: str) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", region_col, network_col])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region_id"] = info[region_col].astype(str).str.strip()
    info["network_id"] = info[network_col].astype(str).str.strip()
    info = info.drop_duplicates("sample_id").set_index("sample_id")
    return info[info["region_id"].ne("") & info["network_id"].ne("")].copy()


def network_row(sample_id: str, truth_network: str, network_top: list[str], n_genes: int) -> dict[str, Any]:
    true_rank = network_top.index(truth_network) + 1 if truth_network in network_top else len(network_top) + 1
    padded = network_top[:3] + [""] * max(0, 3 - len(network_top))
    return {
        "route_family": "hybrid_projected_network_logcpm_exact",
        "sample_id": sample_id,
        "label": truth_network,
        "pred_top1": padded[0],
        "pred_top2": padded[1],
        "pred_top3": padded[2],
        "hit1": int(padded[0] == truth_network),
        "hit3": int(truth_network in padded),
        "true_rank": int(true_rank),
        "n_network_genes": int(n_genes),
    }


def exact_row(sample_id: str, truth_region: str, truth_network: str, network_top: list[str], ranked: list[str], n_total_regions: int) -> dict[str, Any]:
    true_rank = ranked.index(truth_region) + 1 if truth_region in ranked else n_total_regions + 1
    padded = ranked[:3] + [""] * max(0, 3 - len(ranked))
    return {
        "route": EXACT_ROUTE,
        "route_family": "hybrid_projected_network_logcpm_exact",
        "sample_id": sample_id,
        "label": truth_region,
        "true_network": truth_network,
        "network_beam": " | ".join(network_top),
        "network_top1_hit": int(network_top[0] == truth_network),
        "network_top3_hit": int(truth_network in network_top),
        "pred_top1": padded[0],
        "pred_top2": padded[1],
        "pred_top3": padded[2],
        "true_rank": int(true_rank),
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "n_candidate_regions": int(len(ranked)),
    }


def summarize_network(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "route_family": "hybrid_projected_network_logcpm_exact",
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "median_true_rank": float(detail["true_rank"].median()),
    }


def summarize_exact(detail: pd.DataFrame) -> dict[str, Any]:
    return {
        "route_family": "hybrid_projected_network_logcpm_exact",
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "median_true_rank": float(detail["true_rank"].median()),
        "conditional_top1_given_network_top1": float(detail.loc[detail["network_top1_hit"] == 1, "hit1"].mean()),
        "conditional_top3_given_network_top1": float(detail.loc[detail["network_top1_hit"] == 1, "hit3"].mean()),
    }


def summarize_group(detail: pd.DataFrame) -> dict[str, Any]:
    low = detail[detail["pred_top1_resolution_tier"] == "low_resolution"]
    high = detail[detail["pred_top1_resolution_tier"] == "high_resolution"]
    return {
        "route_family": "hybrid_projected_network_logcpm_exact",
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Complete three-tier LOSO for hybrid projected-network/logCPM-exact route.")
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_MODEL_GENES)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--exact-fusion-weight", type=float, default=0.25)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--group-min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.15)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=8)
    parser.add_argument("--max-folds", type=int, default=0)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata(args.sample_info, args.sample_sheet, args.region_col, args.network_col)
    samples = [sample for sample in samples if sample in metadata.index]
    counts = counts.loc[genes, samples]
    vsd = vsd.loc[genes, samples]
    logcpm = compute_logcpm(counts)

    gene_to_idx = {gene: i for i, gene in enumerate(genes)}
    locked_genes = [gene for gene in read_locked_model_genes(args.locked_model_genes) if gene in gene_to_idx]
    locked_rows = np.asarray([gene_to_idx[gene] for gene in locked_genes], dtype=int)
    if len(locked_rows) < 20:
        locked_rows = np.arange(len(genes), dtype=int)

    region_labels = metadata.loc[samples, "region_id"].astype(str).to_numpy()
    network_labels = metadata.loc[samples, "network_id"].astype(str).to_numpy()
    networks = sorted(set(network_labels))
    all_regions = sorted(set(region_labels))
    vsd_values = vsd.to_numpy(dtype=np.float32)
    logcpm_values = logcpm.to_numpy(dtype=np.float32)

    network_rows: list[dict[str, Any]] = []
    exact_rows: list[dict[str, Any]] = []
    group_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []
    unsupported_rows: list[dict[str, Any]] = []
    max_folds = len(samples) if args.max_folds == 0 else min(args.max_folds, len(samples))

    for fold_no, sample_idx in enumerate(range(max_folds), start=1):
        train_idx = np.setdiff1d(np.arange(len(samples), dtype=int), np.asarray([sample_idx], dtype=int))
        sample_id = samples[sample_idx]
        truth_region = str(region_labels[sample_idx])
        truth_network = str(network_labels[sample_idx])
        train_regions = sorted(set(region_labels[train_idx]))
        region_training = {
            region: train_idx[region_labels[train_idx] == region]
            for region in train_regions
        }

        network_reference = build_centroids(vsd_values[locked_rows, :], network_labels, networks, train_idx)
        projected_locked = loo_project_rows(logcpm_values, vsd_values, locked_rows, sample_idx)
        network_scores = correlation_scores(network_reference, projected_locked)
        network_top = [networks[int(i)] for i in np.argsort(network_scores)[::-1][:3].tolist()]
        network_rows.append(network_row(sample_id, truth_network, network_top, len(locked_rows)))

        region_evaluable = truth_region in region_training
        fold_row = {
            "fold_no": fold_no,
            "sample_id": sample_id,
            "n_train_samples": int(len(train_idx)),
            "network_evaluable": True,
            "region_evaluable": bool(region_evaluable),
            "region_non_evaluable_reason": "" if region_evaluable else "truth_region_absent_from_training_fold",
            "n_candidate_regions": 0,
            "n_local_genes": 0,
        }
        if not region_evaluable:
            unsupported_rows.append(
                {
                    "fold_no": fold_no,
                    "sample_id": sample_id,
                    "truth_network": truth_network,
                    "truth_region": truth_region,
                    "reason": "truth_region_absent_from_training_fold",
                    "network_included_in_evaluation": True,
                    "resolution_group_included_in_evaluation": False,
                    "exact_region_included_in_evaluation": False,
                }
            )
            fold_rows.append(fold_row)
            continue

        candidates = sorted(
            region
            for region in train_regions
            if metadata.loc[samples[int(region_training[region][0])], "network_id"] in set(network_top)
        )
        if len(candidates) < 2:
            fold_row["region_evaluable"] = False
            fold_row["region_non_evaluable_reason"] = "fewer_than_two_candidate_regions_in_network_beam"
            unsupported_rows.append(
                {
                    "fold_no": fold_no,
                    "sample_id": sample_id,
                    "truth_network": truth_network,
                    "truth_region": truth_region,
                    "reason": "fewer_than_two_candidate_regions_in_network_beam",
                    "network_included_in_evaluation": True,
                    "resolution_group_included_in_evaluation": False,
                    "exact_region_included_in_evaluation": False,
                }
            )
            fold_rows.append(fold_row)
            continue
        candidate_training = {region: region_training[region] for region in candidates}
        gene_order, _ = select_group_discriminative_genes(
            logcpm_values,
            candidates,
            candidate_training,
            max(100, args.local_top_n_genes),
        )
        candidate_reference = build_region_reference(logcpm_values, candidates, candidate_training)
        sample_vec = logcpm_values[:, sample_idx]
        scores50 = correlation_scores(candidate_reference, sample_vec, gene_order[: min(50, len(gene_order))])
        scores100 = correlation_scores(candidate_reference, sample_vec, gene_order[: min(100, len(gene_order))])
        fused = args.exact_fusion_weight * zscore(scores50) + (1.0 - args.exact_fusion_weight) * zscore(scores100)
        ranked_exact = [candidates[i] for i in np.argsort(fused)[::-1].tolist()]
        exact_rows.append(exact_row(sample_id, truth_region, truth_network, network_top, ranked_exact, len(all_regions)))

        local_rows = gene_order[: min(args.local_top_n_genes, len(gene_order))]
        assignment: dict[str, str | None] = {}
        train_meta = metadata.loc[[samples[int(i)] for i in train_idx]]
        for region in candidates:
            nets = sorted(train_meta.loc[train_meta["region_id"] == region, "network_id"].astype(str).unique())
            assignment[region] = nets[0] if len(nets) == 1 else None
        annotations, _ = build_resolution_groups(
            logcpm_values,
            candidates,
            candidate_training,
            assignment,
            local_rows,
            args.min_resolution_samples,
            args.min_merge_samples,
            args.group_min_pair_errors,
            args.min_confusion_rate,
            args.similarity_threshold,
            args.merge_similarity_threshold,
            args.max_group_size,
        )
        group_scores = correlation_scores(candidate_reference, sample_vec, local_rows)
        ranked_group = [candidates[i] for i in np.argsort(group_scores)[::-1].tolist()]
        group_detail = score_route(
            GROUP_ROUTE,
            sample_id,
            truth_region,
            truth_network,
            network_top,
            ranked_group,
            annotations,
            len(all_regions),
        )
        group_detail["route_family"] = "hybrid_projected_network_logcpm_exact"
        group_rows.append(group_detail)
        fold_row["n_candidate_regions"] = int(len(candidates))
        fold_row["n_local_genes"] = int(len(local_rows))
        fold_rows.append(fold_row)

    network_df = pd.DataFrame(network_rows)
    exact_df = pd.DataFrame(exact_rows)
    group_df = pd.DataFrame(group_rows)
    network_summary = pd.DataFrame([summarize_network(network_df)])
    exact_summary = pd.DataFrame([summarize_exact(exact_df)])
    group_summary = pd.DataFrame([summarize_group(group_df)])

    network_df.to_csv(args.outdir / "hybrid_formal_loso_network_detail.csv", index=False)
    exact_df.to_csv(args.outdir / "hybrid_formal_loso_exact_region_detail.csv", index=False)
    group_df.to_csv(args.outdir / "hybrid_formal_loso_resolution_group_detail.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(args.outdir / "hybrid_formal_loso_fold_summary.csv", index=False)
    pd.DataFrame(unsupported_rows).to_csv(
        args.outdir / "hybrid_formal_loso_region_unsupported_samples.csv",
        index=False,
    )
    network_summary.to_csv(args.outdir / "hybrid_formal_loso_network_route_metrics.csv", index=False)
    exact_summary.to_csv(args.outdir / "hybrid_formal_loso_exact_region_route_metrics.csv", index=False)
    group_summary.to_csv(args.outdir / "hybrid_formal_loso_resolution_group_route_metrics.csv", index=False)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_design": "three-tier LOSO; Network is evaluated for every supported Network label, while resolution-group and exact-region metrics require the truth region to remain represented in the training fold",
        "evaluation_denominators": {
            "network_n": int(len(network_df)),
            "resolution_group_n": int(len(group_df)),
            "exact_region_n": int(len(exact_df)),
            "region_unsupported_n": int(len(unsupported_rows)),
            "region_unsupported_reason": "truth region absent from the training fold or fewer than two candidate regions in the Network Top3 beam",
        },
        "network_implementation": "fold-local logCPM-to-VSD projection with locked Network genes and correlation-ranked Network Top3; no pairwise rescue in this LOSO script",
        "region_implementation": "logCPM-compatible local resolution-group and exact-region reranking within the projected-VSD Network Top3 beam",
        "network": network_summary.to_dict(orient="records")[0],
        "resolution_group": group_summary.to_dict(orient="records")[0],
        "exact_region": exact_summary.to_dict(orient="records")[0],
    }
    write_json(args.outdir / "hybrid_formal_loso_summary.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
