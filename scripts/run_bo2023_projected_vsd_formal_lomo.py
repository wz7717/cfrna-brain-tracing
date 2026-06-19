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
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_SAMPLE_INFO, DEFAULT_VSD  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_network_pairwise_correlation_validation import (  # noqa: E402
    PAIR_TOP3_ROUTE,
    build_pair_models,
    derive_training_confusion_pairs,
    evaluate_pairwise_rescue,
)
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP  # noqa: E402
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups, score_route  # noqa: E402
from scripts.run_bo2023_hierarchical_region_correlation_validation import rank_candidates  # noqa: E402
from scripts.run_bo2023_loso_validation import correlation_scores  # noqa: E402
from scripts.run_bo2023_leave_one_monkey_out_validation import (  # noqa: E402
    build_label_reference,
    build_region_reference,
    build_region_training,
    make_exact_detail,
    summarize_exact,
    summarize_group,
    summarize_network,
    zscore,
)
from scripts.run_bo2023_projected_vsd_region_local_rerank import fit_project_rows  # noqa: E402


EXACT_ROUTE = "top3_beam_local_top50_top100_zfusion_w0p25"
GROUP_ROUTE = "top3_network_beam_local_region_candidates"


def read_metadata(path: Path, sheet: str, region_col: str, network_col: str, monkey_col: str) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", region_col, network_col, monkey_col])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region_id"] = info[region_col].astype(str).str.strip()
    info["endpoint_label"] = info[network_col].astype(str).str.strip()
    info["monkey_id"] = info[monkey_col].astype(str).str.strip()
    info = info.drop_duplicates("sample_id")
    return info[
        info["sample_id"].ne("")
        & info["region_id"].ne("")
        & info["endpoint_label"].ne("")
        & info["monkey_id"].ne("")
    ].copy()


def sample_probability(detail_prob: dict[str, Any], monkey_id: str, route_family: str) -> dict[str, Any]:
    out = dict(detail_prob)
    out["monkey_id"] = monkey_id
    out["route_family"] = route_family
    return out


def add_route_family(row: dict[str, Any], monkey_id: str, route_family: str, route: str | None = None) -> dict[str, Any]:
    out = dict(row)
    out["monkey_id"] = monkey_id
    out["route_family"] = route_family
    if route is not None:
        out["route"] = route
    return out


def per_monkey_metrics(detail: pd.DataFrame, metric_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (route_family, monkey_id), frame in detail.groupby(["route_family", "monkey_id"], sort=True):
        row: dict[str, Any] = {"route_family": route_family, "monkey_id": monkey_id, "n": int(len(frame))}
        for col in metric_cols:
            row[f"{col}_mean"] = float(frame[col].mean())
            row[f"{col}_hits"] = int(frame[col].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Projected-VSD formal three-tier LOMO validation.")
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--monkey-col", default="MonkeyID")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--global-top-n-genes", type=int, default=200)
    parser.add_argument("--gene-pool-size", type=int, default=1000)
    parser.add_argument("--pair-top-n-genes", type=int, default=100)
    parser.add_argument("--max-pairs-per-truth", type=int, default=2)
    parser.add_argument("--min-pair-errors", type=int, default=3)
    parser.add_argument("--pair-min-margin", type=float, default=0.002)
    parser.add_argument("--exact-fusion-weight", type=float, default=0.25)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--min-resolution-samples", type=int, default=8)
    parser.add_argument("--min-merge-samples", type=int, default=3)
    parser.add_argument("--group-min-pair-errors", type=int, default=2)
    parser.add_argument("--min-confusion-rate", type=float, default=0.15)
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    parser.add_argument("--merge-similarity-threshold", type=float, default=0.90)
    parser.add_argument("--max-group-size", type=int, default=8)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata(args.sample_info, args.sample_sheet, args.region_col, args.network_col, args.monkey_col)
    samples = [sample for sample in samples if sample in set(metadata["sample_id"])]
    counts = counts.loc[genes, samples]
    vsd = vsd.loc[genes, samples]
    logcpm = compute_logcpm(counts)

    sample_pos = {sample_id: idx for idx, sample_id in enumerate(samples)}
    sample_ann = metadata.set_index("sample_id").reindex(samples)
    network_labels = sample_ann["endpoint_label"].to_numpy(dtype=str)
    monkey_labels = sample_ann["monkey_id"].to_numpy(dtype=str)
    groups = sorted(set(network_labels))
    all_regions = sorted(sample_ann["region_id"].dropna().astype(str).unique().tolist())

    vsd_values = vsd.to_numpy(dtype=np.float32)
    logcpm_values = logcpm.to_numpy(dtype=np.float32)

    network_details: list[dict[str, Any]] = []
    network_probs: list[dict[str, Any]] = []
    exact_rows: list[dict[str, Any]] = []
    group_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []

    for fold_no, monkey_id in enumerate(sorted(set(monkey_labels)), start=1):
        test_indices = np.flatnonzero(monkey_labels == monkey_id)
        train_indices = np.flatnonzero(monkey_labels != monkey_id)
        train_sample_ids = {samples[i] for i in train_indices}
        test_sample_ids = [samples[i] for i in test_indices]

        train_ann = metadata[metadata["sample_id"].isin(train_sample_ids)].copy()
        region_training = build_region_training(metadata.rename(columns={"region_id": "region_id"}), sample_pos, train_sample_ids)

        projected_test = fit_project_rows(
            logcpm_values,
            vsd_values,
            np.arange(len(genes), dtype=int),
            train_indices,
            test_indices,
        )

        route_spaces = {
            "native_vsd": {
                "network_train_values": vsd_values,
                "network_test_values": vsd_values[:, test_indices],
                "exact_train_values": vsd_values,
                "exact_test_values": vsd_values[:, test_indices],
            },
            "projected_vsd": {
                "network_train_values": vsd_values,
                "network_test_values": projected_test,
                "exact_train_values": vsd_values,
                "exact_test_values": projected_test,
            },
            "logcpm_baseline": {
                "network_train_values": logcpm_values,
                "network_test_values": logcpm_values[:, test_indices],
                "exact_train_values": logcpm_values,
                "exact_test_values": logcpm_values[:, test_indices],
            },
            "hybrid_projected_network_logcpm_exact": {
                "network_train_values": vsd_values,
                "network_test_values": projected_test,
                "exact_train_values": logcpm_values,
                "exact_test_values": logcpm_values[:, test_indices],
            },
        }

        fold_info: dict[str, Any] = {
            "fold": fold_no,
            "heldout_monkey_id": monkey_id,
            "n_train_samples": int(len(train_indices)),
            "n_test_samples": int(len(test_indices)),
            "n_train_regions": int(len(region_training)),
        }

        for route_family, route_data in route_spaces.items():
            network_values = route_data["network_train_values"]
            network_test_values = route_data["network_test_values"]
            exact_values = route_data["exact_train_values"]
            exact_test_values = route_data["exact_test_values"]
            network_reference, network_training = build_label_reference(network_values, network_labels, groups, train_indices)
            gene_pool, _ = select_group_discriminative_genes(network_values, groups, network_training, args.gene_pool_size)
            global_genes = gene_pool[: args.global_top_n_genes]
            pairs, _ = derive_training_confusion_pairs(
                network_values,
                network_labels,
                groups,
                network_training,
                network_reference,
                global_genes,
                args.max_pairs_per_truth,
                args.min_pair_errors,
            )
            pair_models, _ = build_pair_models(
                network_values,
                network_training,
                network_reference,
                pairs,
                groups,
                gene_pool,
                args.pair_top_n_genes,
            )

            fold_exact = 0
            fold_group = 0
            for local_pos, sample_idx in enumerate(test_indices):
                sample_id = samples[int(sample_idx)]
                truth_network = str(network_labels[int(sample_idx)])
                truth_region = str(sample_ann.loc[sample_id, "region_id"])
                network_sample_vec = network_test_values[:, local_pos]
                exact_sample_vec = exact_test_values[:, local_pos]
                scores = correlation_scores(network_reference, network_sample_vec, global_genes)
                pair_detail, pair_prob = evaluate_pairwise_rescue(
                    PAIR_TOP3_ROUTE,
                    sample_id,
                    truth_network,
                    network_sample_vec,
                    network_reference,
                    groups,
                    scores,
                    pair_models,
                    3,
                    args.pair_min_margin,
                )
                network_details.append(add_route_family(pair_detail, monkey_id, route_family))
                network_probs.append(sample_probability(pair_prob, monkey_id, route_family))

                network_top = [pair_detail[f"pred_top{i}"] for i in [1, 2, 3]]
                if truth_region not in region_training:
                    continue
                candidates = sorted(
                    region
                    for region in train_ann.loc[
                        train_ann["endpoint_label"].isin(network_top), "region_id"
                    ].unique().tolist()
                    if region in region_training
                )
                if not candidates:
                    continue
                candidate_training = {region: region_training[region] for region in candidates}
                candidate_reference = build_region_reference(exact_values, candidates, candidate_training)
                if len(candidates) >= 2:
                    gene_order, _ = select_group_discriminative_genes(
                        exact_values,
                        candidates,
                        candidate_training,
                        max(100, args.local_top_n_genes),
                    )
                else:
                    gene_order = np.arange(exact_values.shape[0], dtype=int)
                scores50 = correlation_scores(
                    candidate_reference,
                    exact_sample_vec,
                    gene_order[: min(50, len(gene_order))],
                )
                scores100 = correlation_scores(
                    candidate_reference,
                    exact_sample_vec,
                    gene_order[: min(100, len(gene_order))],
                )
                fused = args.exact_fusion_weight * zscore(scores50) + (1.0 - args.exact_fusion_weight) * zscore(scores100)
                ranked_exact = [candidates[i] for i in np.argsort(fused)[::-1].tolist()]
                exact_rows.append(
                    add_route_family(
                        make_exact_detail(
                            sample_id,
                            monkey_id,
                            truth_region,
                            truth_network,
                            network_top,
                            ranked_exact,
                            len(all_regions),
                        ),
                        monkey_id,
                        route_family,
                        EXACT_ROUTE,
                    )
                )
                fold_exact += 1

                local_rows = gene_order[: min(args.local_top_n_genes, len(gene_order))]
                assignment: dict[str, str | None] = {}
                for region in candidates:
                    nets = sorted(train_ann.loc[train_ann["region_id"] == region, "endpoint_label"].astype(str).unique())
                    assignment[region] = nets[0] if len(nets) == 1 else None
                annotations, _ = build_resolution_groups(
                    exact_values,
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
                group_scores = correlation_scores(candidate_reference, exact_sample_vec, local_rows)
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
                group_rows.append(add_route_family(group_detail, monkey_id, route_family))
                fold_group += 1

            fold_info[f"{route_family}_n_exact_evaluable"] = int(fold_exact)
            fold_info[f"{route_family}_n_group_evaluable"] = int(fold_group)
        fold_rows.append(fold_info)

    network_df = pd.DataFrame(network_details)
    network_prob_df = pd.DataFrame(network_probs)
    exact_df = pd.DataFrame(exact_rows)
    group_df = pd.DataFrame(group_rows)

    network_summary_rows = []
    exact_summary_rows = []
    group_summary_rows = []
    for route_family in sorted(network_df["route_family"].unique()):
        nd = network_df[network_df["route_family"] == route_family].copy()
        npb = network_prob_df[network_prob_df["route_family"] == route_family].copy()
        ed = exact_df[exact_df["route_family"] == route_family].copy()
        gd = group_df[group_df["route_family"] == route_family].copy()
        network_summary_rows.append({"route_family": route_family, **summarize_network(nd, npb, groups)})
        exact_summary_rows.append({"route_family": route_family, **summarize_exact(ed)})
        group_summary_rows.append({"route_family": route_family, **summarize_group(gd)})

    network_summary = pd.DataFrame(network_summary_rows)
    exact_summary = pd.DataFrame(exact_summary_rows)
    group_summary = pd.DataFrame(group_summary_rows)

    pd.DataFrame(fold_rows).to_csv(args.outdir / "formal_lomo_fold_summary.csv", index=False)
    network_df.to_csv(args.outdir / "formal_lomo_network_detail.csv", index=False)
    exact_df.to_csv(args.outdir / "formal_lomo_exact_region_detail.csv", index=False)
    group_df.to_csv(args.outdir / "formal_lomo_resolution_group_detail.csv", index=False)
    network_summary.to_csv(args.outdir / "formal_lomo_network_route_metrics.csv", index=False)
    exact_summary.to_csv(args.outdir / "formal_lomo_exact_region_route_metrics.csv", index=False)
    group_summary.to_csv(args.outdir / "formal_lomo_resolution_group_route_metrics.csv", index=False)
    per_monkey_metrics(network_df, ["hit1", "hit3"]).to_csv(args.outdir / "formal_lomo_network_per_monkey_metrics.csv", index=False)
    per_monkey_metrics(exact_df, ["hit1", "hit3"]).to_csv(args.outdir / "formal_lomo_exact_region_per_monkey_metrics.csv", index=False)
    per_monkey_metrics(group_df, ["group_hit1", "group_hit3", "hit1", "hit3"]).to_csv(
        args.outdir / "formal_lomo_resolution_group_per_monkey_metrics.csv", index=False
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_design": "formal three-tier leave-one-monkey-out; route-specific network, resolution-group, and exact-region components rebuilt fold-locally",
        "routes": {
            "network": network_summary.to_dict(orient="records"),
            "resolution_group": group_summary.to_dict(orient="records"),
            "exact_region": exact_summary.to_dict(orient="records"),
        },
        "parameters": {
            "global_top_n_genes": args.global_top_n_genes,
            "gene_pool_size": args.gene_pool_size,
            "pair_top_n_genes": args.pair_top_n_genes,
            "pair_min_margin": args.pair_min_margin,
            "exact_route": EXACT_ROUTE,
            "resolution_group_route": GROUP_ROUTE,
        },
    }
    write_json(args.outdir / "formal_lomo_validation_summary.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
