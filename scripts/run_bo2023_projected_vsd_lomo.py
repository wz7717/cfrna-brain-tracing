#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import (  # noqa: E402
    align_matrices,
    apply_projector,
    compute_logcpm,
    fit_linear_projector,
    map_index_to_symbols,
    read_bo2023_gene_matrix,
    read_gene_map,
    write_json,
)
from scripts.build_bo2023_reference_projector import (  # noqa: E402
    DEFAULT_COUNTS,
    DEFAULT_GENE_MAP,
    DEFAULT_MODEL_GENES,
    DEFAULT_SAMPLE_INFO,
    DEFAULT_VSD,
    read_locked_model_genes,
)
from scripts.run_bo2023_projected_vsd_loso import corr_scores, make_rank_row  # noqa: E402


def read_metadata(path: Path, sheet: str, monkey_col: str, network_col: str) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", monkey_col, network_col])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["monkey_id"] = info[monkey_col].astype(str).str.strip()
    info["network"] = info[network_col].astype(str).str.strip()
    info = info.drop_duplicates("sample_id").set_index("sample_id")
    return info[info["monkey_id"].ne("") & info["network"].ne("")]


def build_centroids(values: pd.DataFrame, metadata: pd.DataFrame, train_samples: list[str], groups: list[str]) -> np.ndarray:
    cols = []
    for group in groups:
        group_samples = [sample for sample in train_samples if metadata.loc[sample, "network"] == group]
        cols.append(values[group_samples].mean(axis=1).to_numpy(dtype=np.float32))
    return np.column_stack(cols).astype(np.float32)


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict fold-local LOMO validation for projected Bo2023-like VSD.")
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--monkey-col", default="MonkeyID")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_MODEL_GENES)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata(args.sample_info, args.sample_sheet, args.monkey_col, args.network_col)
    samples = [sample for sample in samples if sample in metadata.index]
    locked_genes = read_locked_model_genes(args.locked_model_genes)
    selected_genes = [gene for gene in locked_genes if gene in set(genes)]
    if len(selected_genes) < 20:
        selected_genes = genes
    counts = counts.loc[selected_genes, samples]
    vsd = vsd.loc[selected_genes, samples]
    logcpm = compute_logcpm(counts)

    rows = []
    fold_rows = []
    monkey_ids = sorted(metadata.loc[samples, "monkey_id"].unique().tolist())
    for fold_no, monkey_id in enumerate(monkey_ids, start=1):
        test_samples = [sample for sample in samples if metadata.loc[sample, "monkey_id"] == monkey_id]
        train_samples = [sample for sample in samples if metadata.loc[sample, "monkey_id"] != monkey_id]
        groups = sorted(set(metadata.loc[train_samples, "network"].astype(str)))
        train_regions = int(metadata.loc[train_samples, "network"].nunique())
        test_regions = int(metadata.loc[test_samples, "network"].nunique())
        fold_rows.append(
            {
                "heldout_monkey_id": monkey_id,
                "n_train_samples": len(train_samples),
                "n_test_samples": len(test_samples),
                "n_train_regions": train_regions,
                "n_test_regions": test_regions,
            }
        )
        fit, _ = fit_linear_projector(logcpm[train_samples], vsd[train_samples])
        projected = apply_projector(fit, logcpm[test_samples])
        native_reference = build_centroids(vsd, metadata, train_samples, groups)
        logcpm_reference = build_centroids(logcpm, metadata, train_samples, groups)
        for sample in test_samples:
            truth = str(metadata.loc[sample, "network"])
            if truth not in groups:
                continue
            rows.append(make_rank_row(fold_no, sample, "native_vsd", truth, groups, corr_scores(native_reference, vsd[sample].to_numpy(dtype=np.float32)), len(selected_genes)))
            rows.append(make_rank_row(fold_no, sample, "projected_vsd", truth, groups, corr_scores(native_reference, projected[sample].to_numpy(dtype=np.float32)), len(selected_genes)))
            rows.append(make_rank_row(fold_no, sample, "logcpm_baseline", truth, groups, corr_scores(logcpm_reference, logcpm[sample].to_numpy(dtype=np.float32)), len(selected_genes)))

    detail = pd.DataFrame(rows)
    detail.to_csv(args.outdir / "bo2023_projected_vsd_lomo_detail.csv", index=False)
    fold_df = pd.DataFrame(fold_rows)
    fold_df.to_csv(args.outdir / "bo2023_projected_vsd_lomo_folds.csv", index=False)
    route_summary = (
        detail.groupby("route")
        .agg(
            n=("sample_id", "size"),
            network_top1=("hit1", "mean"),
            network_top3=("hit3", "mean"),
            median_true_rank=("true_rank", "median"),
            mean_decision_margin=("decision_margin", "mean"),
            median_decision_margin=("decision_margin", "median"),
            n_overlap_genes=("n_overlap_genes", "median"),
        )
        .reset_index()
    )
    route_summary.to_csv(args.outdir / "bo2023_projected_vsd_lomo_route_summary.csv", index=False)
    write_json(
        args.outdir / "bo2023_projected_vsd_lomo_summary.json",
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "strict leave-one-monkey-out SaleemNetworks validation",
            "n_monkey_folds": int(len(monkey_ids)),
            "n_evaluated_rows": int(len(detail)),
            "n_eval_genes": int(len(selected_genes)),
            "folds": fold_df.to_dict(orient="records"),
            "routes": route_summary.to_dict(orient="records"),
        },
    )
    print(f"Wrote LOMO outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
