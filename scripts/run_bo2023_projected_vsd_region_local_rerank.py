#!/usr/bin/env python
from __future__ import annotations

import argparse
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
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_SAMPLE_INFO, DEFAULT_VSD, read_locked_model_genes  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes  # noqa: E402
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP, DEFAULT_CLEANED_LOCKED_GENES  # noqa: E402
from scripts.run_bo2023_projected_vsd_loso import corr_scores  # noqa: E402


def read_metadata(path: Path, sheet: str, region_col: str, network_col: str, monkey_col: str) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", region_col, network_col, monkey_col])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region"] = info[region_col].astype(str).str.strip()
    info["network"] = info[network_col].astype(str).str.strip()
    info["monkey_id"] = info[monkey_col].astype(str).str.strip()
    info = info.drop_duplicates("sample_id").set_index("sample_id")
    return info[info["region"].ne("") & info["network"].ne("")]


def build_centroids(values: np.ndarray, labels: np.ndarray, groups: list[str], sample_indices: np.ndarray) -> np.ndarray:
    columns = []
    for group in groups:
        idx = sample_indices[labels[sample_indices] == group]
        columns.append(values[:, idx].mean(axis=1, dtype=np.float64))
    return np.column_stack(columns).astype(np.float32)


def rank_row(
    split: str,
    fold_no: int,
    heldout_id: str,
    sample_id: str,
    route: str,
    truth: str,
    network_top3: list[str],
    regions: list[str],
    scores: np.ndarray,
    n_genes: int,
) -> dict[str, Any]:
    order = np.argsort(scores)[::-1]
    ranked = [regions[int(i)] for i in order]
    true_rank = ranked.index(truth) + 1 if truth in ranked else len(ranked) + 1
    padded = ranked[:3] + [""] * max(0, 3 - len(ranked))
    return {
        "split": split,
        "fold_no": fold_no,
        "heldout_id": heldout_id,
        "sample_id": sample_id,
        "route": route,
        "label": truth,
        "network_beam": " | ".join(network_top3),
        "pred_top1": padded[0],
        "pred_top2": padded[1],
        "pred_top3": padded[2],
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "true_rank": int(true_rank),
        "score_top1": float(scores[int(order[0])]) if len(order) else float("nan"),
        "score_top2": float(scores[int(order[1])]) if len(order) > 1 else float("nan"),
        "decision_margin": float(scores[int(order[0])] - scores[int(order[1])]) if len(order) > 1 else float("nan"),
        "n_candidate_regions": int(len(regions)),
        "n_local_genes": int(n_genes),
        "abstained": 0,
        "abstain_reason": "",
    }


def abstain_row(
    split: str,
    fold_no: int,
    heldout_id: str,
    sample_id: str,
    route: str,
    truth: str,
    reason: str,
) -> dict[str, Any]:
    return {
        "split": split,
        "fold_no": fold_no,
        "heldout_id": heldout_id,
        "sample_id": sample_id,
        "route": route,
        "label": truth,
        "network_beam": "",
        "pred_top1": "",
        "pred_top2": "",
        "pred_top3": "",
        "hit1": 0,
        "hit3": 0,
        "true_rank": -1,
        "score_top1": float("nan"),
        "score_top2": float("nan"),
        "decision_margin": float("nan"),
        "n_candidate_regions": 0,
        "n_local_genes": 0,
        "abstained": 1,
        "abstain_reason": reason,
    }


def loo_project_rows(logcpm: np.ndarray, vsd: np.ndarray, rows: np.ndarray, sample_idx: int) -> np.ndarray:
    x = logcpm[rows, :]
    y = vsd[rows, :]
    mask = np.ones(x.shape[1], dtype=bool)
    mask[sample_idx] = False
    xt = x[:, mask].astype(np.float64)
    yt = y[:, mask].astype(np.float64)
    x_mean = xt.mean(axis=1)
    y_mean = yt.mean(axis=1)
    xc = xt - x_mean[:, None]
    yc = yt - y_mean[:, None]
    denom = np.square(xc).sum(axis=1)
    slope = np.divide((xc * yc).sum(axis=1), denom, out=np.zeros(len(rows), dtype=float), where=denom > 0)
    intercept = y_mean - slope * x_mean
    pred = slope * logcpm[rows, sample_idx].astype(np.float64) + intercept
    low = np.quantile(yt, 0.005, axis=1)
    high = np.quantile(yt, 0.995, axis=1)
    return np.clip(pred, low, high).astype(np.float32)


def fit_project_rows(logcpm: np.ndarray, vsd: np.ndarray, rows: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray) -> np.ndarray:
    x = logcpm[rows[:, None], train_idx[None, :]].astype(np.float64)
    y = vsd[rows[:, None], train_idx[None, :]].astype(np.float64)
    x_mean = x.mean(axis=1)
    y_mean = y.mean(axis=1)
    xc = x - x_mean[:, None]
    yc = y - y_mean[:, None]
    denom = np.square(xc).sum(axis=1)
    slope = np.divide((xc * yc).sum(axis=1), denom, out=np.zeros(len(rows), dtype=float), where=denom > 0)
    intercept = y_mean - slope * x_mean
    pred = slope[:, None] * logcpm[rows[:, None], test_idx[None, :]].astype(np.float64) + intercept[:, None]
    low = np.quantile(y, 0.005, axis=1)
    high = np.quantile(y, 0.995, axis=1)
    return np.clip(pred, low[:, None], high[:, None]).astype(np.float32)


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split, route), frame in detail.groupby(["split", "route"], sort=True):
        valid = frame[frame["abstained"] == 0]
        rows.append(
            {
                "split": split,
                "route": route,
                "n_total": int(len(frame)),
                "n_valid": int(len(valid)),
                "abstain_rate": float(frame["abstained"].mean()) if len(frame) else float("nan"),
                "exact_region_top1": float(valid["hit1"].mean()) if len(valid) else float("nan"),
                "exact_region_top3": float(valid["hit3"].mean()) if len(valid) else float("nan"),
                "median_true_rank": float(valid["true_rank"].median()) if len(valid) else float("nan"),
                "mean_candidate_regions": float(valid["n_candidate_regions"].mean()) if len(valid) else float("nan"),
                "median_local_genes": float(valid["n_local_genes"].median()) if len(valid) else float("nan"),
                "median_decision_margin": float(valid["decision_margin"].median()) if len(valid) else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Region-specific fold-local local reranking for projected VSD.")
    parser.add_argument("--splits", nargs="+", choices=["loso", "lomo"], default=["loso", "lomo"])
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--monkey-col", default="MonkeyID")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_CLEANED_LOCKED_GENES)
    parser.add_argument("--local-top-n-genes", type=int, default=200)
    parser.add_argument("--max-folds", type=int, default=0)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata(args.sample_info, args.sample_sheet, args.region_col, args.network_col, args.monkey_col)
    samples = [sample for sample in samples if sample in metadata.index]
    counts = counts.loc[genes, samples]
    vsd = vsd.loc[genes, samples]
    logcpm = compute_logcpm(counts)

    gene_to_idx = {gene: i for i, gene in enumerate(genes)}
    locked_genes = [gene for gene in read_locked_model_genes(args.locked_model_genes) if gene in gene_to_idx]
    locked_rows = np.asarray([gene_to_idx[gene] for gene in locked_genes], dtype=int)
    sample_to_idx = {sample: i for i, sample in enumerate(samples)}
    region_labels = metadata.loc[samples, "region"].astype(str).to_numpy()
    network_labels = metadata.loc[samples, "network"].astype(str).to_numpy()
    monkey_labels = metadata.loc[samples, "monkey_id"].astype(str).to_numpy()
    networks = sorted(set(network_labels))
    all_rows = np.arange(len(genes), dtype=int)
    rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []

    for split in args.splits:
        if split == "loso":
            folds = [(sample, np.asarray([sample_to_idx[sample]], dtype=int)) for sample in samples]
        else:
            folds = [
                (monkey, np.flatnonzero(monkey_labels == monkey))
                for monkey in sorted(set(monkey_labels))
            ]
        if args.max_folds:
            folds = folds[: args.max_folds]

        for fold_no, (heldout_id, test_idx) in enumerate(folds, start=1):
            train_idx = np.setdiff1d(np.arange(len(samples), dtype=int), test_idx, assume_unique=False)
            train_regions = sorted(set(region_labels[train_idx]))
            region_training: dict[str, np.ndarray] = {
                region: train_idx[region_labels[train_idx] == region]
                for region in train_regions
            }
            network_reference_vsd = build_centroids(vsd.to_numpy(dtype=np.float32)[locked_rows, :], network_labels, networks, train_idx)
            network_reference_logcpm = build_centroids(logcpm.to_numpy(dtype=np.float32)[locked_rows, :], network_labels, networks, train_idx)
            fold_rows.append(
                {
                    "split": split,
                    "fold_no": fold_no,
                    "heldout_id": heldout_id,
                    "n_train_samples": int(len(train_idx)),
                    "n_test_samples": int(len(test_idx)),
                    "n_train_regions": int(len(train_regions)),
                    "n_test_regions": int(len(set(region_labels[test_idx]))),
                    "n_train_networks": int(len(set(network_labels[train_idx]))),
                }
            )

            projected_locked = None
            if split == "lomo":
                projected_locked = fit_project_rows(
                    logcpm.to_numpy(dtype=np.float32),
                    vsd.to_numpy(dtype=np.float32),
                    locked_rows,
                    train_idx,
                    test_idx,
                )

            for local_test_pos, sample_idx in enumerate(test_idx):
                sample_id = samples[int(sample_idx)]
                truth_region = str(region_labels[int(sample_idx)])
                if truth_region not in region_training:
                    for route in [
                        "native_vsd_top3_network_local_region_genes",
                        "projected_vsd_top3_network_local_region_genes",
                        "hybrid_projected_network_logcpm_local_region_genes",
                        "logcpm_top3_network_local_region_genes",
                    ]:
                        rows.append(
                            abstain_row(
                                split,
                                fold_no,
                                str(heldout_id),
                                sample_id,
                                route,
                                truth_region,
                                "truth_region_absent_from_training_reference",
                            )
                        )
                    continue

                locked_vsd_sample = vsd.iloc[locked_rows, int(sample_idx)].to_numpy(dtype=np.float32)
                locked_logcpm_sample = logcpm.iloc[locked_rows, int(sample_idx)].to_numpy(dtype=np.float32)
                if split == "loso":
                    locked_projected_sample = loo_project_rows(
                        logcpm.to_numpy(dtype=np.float32),
                        vsd.to_numpy(dtype=np.float32),
                        locked_rows,
                        int(sample_idx),
                    )
                else:
                    locked_projected_sample = projected_locked[:, local_test_pos]

                route_specs = [
                    (
                        "native_vsd_top3_network_local_region_genes",
                        vsd.to_numpy(dtype=np.float32),
                        locked_vsd_sample,
                        network_reference_vsd,
                        "vsd",
                    ),
                    (
                        "projected_vsd_top3_network_local_region_genes",
                        vsd.to_numpy(dtype=np.float32),
                        locked_projected_sample,
                        network_reference_vsd,
                        "projected",
                    ),
                    (
                        "hybrid_projected_network_logcpm_local_region_genes",
                        logcpm.to_numpy(dtype=np.float32),
                        locked_projected_sample,
                        network_reference_vsd,
                        "logcpm",
                    ),
                    (
                        "logcpm_top3_network_local_region_genes",
                        logcpm.to_numpy(dtype=np.float32),
                        locked_logcpm_sample,
                        network_reference_logcpm,
                        "logcpm",
                    ),
                ]

                for route, selection_values, network_sample, network_reference, value_space in route_specs:
                    network_scores = corr_scores(network_reference, network_sample)
                    network_order = np.argsort(network_scores)[::-1]
                    top3_networks = [networks[int(i)] for i in network_order[:3]]
                    candidate_regions = [
                        region for region in train_regions
                        if metadata.loc[samples[int(region_training[region][0])], "network"] in set(top3_networks)
                    ]
                    if truth_region not in candidate_regions:
                        # Keep the candidate-beam restriction explicit. The true region can
                        # be missed by the network gate even if the local scorer is strong.
                        candidate_regions = candidate_regions
                    if len(candidate_regions) < 2:
                        rows.append(
                            abstain_row(
                                split,
                                fold_no,
                                str(heldout_id),
                                sample_id,
                                route,
                                truth_region,
                                "fewer_than_two_candidate_regions",
                            )
                        )
                        continue
                    local_rows, _ = select_group_discriminative_genes(
                        selection_values,
                        sorted(candidate_regions),
                        {region: region_training[region] for region in candidate_regions},
                        args.local_top_n_genes,
                    )
                    candidate_regions_sorted = sorted(candidate_regions)
                    ref = build_centroids(selection_values[local_rows, :], region_labels, candidate_regions_sorted, train_idx)
                    if value_space == "projected":
                        if split == "loso":
                            sample_vec = loo_project_rows(
                                logcpm.to_numpy(dtype=np.float32),
                                vsd.to_numpy(dtype=np.float32),
                                local_rows,
                                int(sample_idx),
                            )
                        else:
                            sample_vec = fit_project_rows(
                                logcpm.to_numpy(dtype=np.float32),
                                vsd.to_numpy(dtype=np.float32),
                                local_rows,
                                train_idx,
                                np.asarray([sample_idx], dtype=int),
                            )[:, 0]
                    else:
                        sample_vec = selection_values[local_rows, int(sample_idx)]
                    scores = corr_scores(ref, sample_vec)
                    rows.append(
                        rank_row(
                            split,
                            fold_no,
                            str(heldout_id),
                            sample_id,
                            route,
                            truth_region,
                            top3_networks,
                            candidate_regions_sorted,
                            scores,
                            len(local_rows),
                        )
                    )

    detail = pd.DataFrame(rows)
    split_tag = "_".join(args.splits)
    detail_path = args.outdir / f"bo2023_projected_vsd_region_local_rerank_{split_tag}_detail.csv"
    detail.to_csv(detail_path, index=False)
    folds_df = pd.DataFrame(fold_rows)
    folds_df.to_csv(args.outdir / f"bo2023_projected_vsd_region_local_rerank_{split_tag}_folds.csv", index=False)
    summary_df = summarize(detail)
    summary_df.to_csv(args.outdir / f"bo2023_projected_vsd_region_local_rerank_{split_tag}_route_summary.csv", index=False)
    write_json(
        args.outdir / f"bo2023_projected_vsd_region_local_rerank_{split_tag}_summary.json",
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "local_top_n_genes": int(args.local_top_n_genes),
            "network_gate": "route-specific Top3 SaleemNetworks centroid correlation using locked 200 genes",
            "region_rerank": "fold-local candidate-region discriminative genes selected from all common genes",
            "n_common_genes": int(len(genes)),
            "n_locked_network_genes": int(len(locked_rows)),
            "routes": summary_df.to_dict(orient="records"),
        },
    )
    print(f"Wrote local rerank outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
