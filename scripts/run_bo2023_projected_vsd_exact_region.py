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
    DEFAULT_SAMPLE_INFO,
    DEFAULT_VSD,
    read_locked_model_genes,
)
from scripts.run_bo2023_projected_vsd_loso import corr_scores, make_rank_row  # noqa: E402


DEFAULT_CLEANED_GENE_MAP = (
    ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv"
)
DEFAULT_CLEANED_LOCKED_GENES = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "bo2023_saleem_network_top200_cleaned_model_genes.csv"
)


def read_metadata(path: Path, sheet: str, region_col: str, monkey_col: str) -> pd.DataFrame:
    usecols = ["No.", region_col]
    if monkey_col:
        usecols.append(monkey_col)
    info = pd.read_excel(path, sheet_name=sheet, usecols=usecols)
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region"] = info[region_col].astype(str).str.strip()
    if monkey_col:
        info["monkey_id"] = info[monkey_col].astype(str).str.strip()
    else:
        info["monkey_id"] = ""
    info = info.drop_duplicates("sample_id").set_index("sample_id")
    return info[info["region"].ne("") & info["region"].ne("nan")]


def build_centroids(values: pd.DataFrame, metadata: pd.DataFrame, train_samples: list[str], regions: list[str]) -> np.ndarray:
    columns = []
    for region in regions:
        samples = [sample for sample in train_samples if metadata.loc[sample, "region"] == region]
        columns.append(values[samples].mean(axis=1).to_numpy(dtype=np.float32))
    return np.column_stack(columns).astype(np.float32)


def abstain_row(
    fold_no: int,
    sample: str,
    route: str,
    truth: str,
    reason: str,
    n_genes: int,
) -> dict[str, object]:
    return {
        "fold_no": fold_no,
        "sample_id": sample,
        "route": route,
        "label": truth,
        "pred_top1": "",
        "pred_top2": "",
        "pred_top3": "",
        "hit1": 0,
        "hit3": 0,
        "true_rank": -1,
        "score_top1": np.nan,
        "score_top2": np.nan,
        "decision_margin": np.nan,
        "n_overlap_genes": int(n_genes),
        "abstained": 1,
        "abstain_reason": reason,
    }


def valid_rank_row(*args: object, **kwargs: object) -> dict[str, object]:
    row = make_rank_row(*args, **kwargs)
    row["abstained"] = 0
    row["abstain_reason"] = ""
    return row


def route_summary(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for route, frame in detail.groupby("route", sort=True):
        valid = frame[frame["abstained"] == 0]
        rows.append(
            {
                "route": route,
                "n_total": int(len(frame)),
                "n_valid": int(len(valid)),
                "abstain_rate": float(frame["abstained"].mean()) if len(frame) else float("nan"),
                "exact_region_top1": float(valid["hit1"].mean()) if len(valid) else float("nan"),
                "exact_region_top3": float(valid["hit3"].mean()) if len(valid) else float("nan"),
                "median_true_rank": float(valid["true_rank"].median()) if len(valid) else float("nan"),
                "mean_decision_margin": float(valid["decision_margin"].mean()) if len(valid) else float("nan"),
                "median_decision_margin": float(valid["decision_margin"].median()) if len(valid) else float("nan"),
                "n_overlap_genes": float(valid["n_overlap_genes"].median()) if len(valid) else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Fold-local exact-region validation for projected Bo2023-like VSD.")
    parser.add_argument("--split", choices=["loso", "lomo"], required=True)
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--monkey-col", default="MonkeyID")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_CLEANED_LOCKED_GENES)
    parser.add_argument("--gene-panel", choices=["locked", "all"], default="locked")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata(args.sample_info, args.sample_sheet, args.region_col, args.monkey_col)
    samples = [sample for sample in samples if sample in metadata.index]

    if args.gene_panel == "locked":
        locked_genes = read_locked_model_genes(args.locked_model_genes)
        selected_genes = [gene for gene in locked_genes if gene in set(genes)]
        if len(selected_genes) < 20:
            raise ValueError("Fewer than 20 locked model genes overlap the Bo2023 common panel.")
    else:
        selected_genes = genes

    counts = counts.loc[selected_genes, samples]
    vsd = vsd.loc[selected_genes, samples]
    logcpm = compute_logcpm(counts)

    if args.split == "loso":
        folds = [(sample, [sample]) for sample in samples]
    else:
        folds = [
            (monkey_id, [sample for sample in samples if metadata.loc[sample, "monkey_id"] == monkey_id])
            for monkey_id in sorted(metadata.loc[samples, "monkey_id"].unique().tolist())
        ]

    rows: list[dict[str, object]] = []
    fold_rows: list[dict[str, object]] = []
    for fold_no, (fold_id, test_samples) in enumerate(folds, start=1):
        test_set = set(test_samples)
        train_samples = [sample for sample in samples if sample not in test_set]
        train_regions = sorted(set(metadata.loc[train_samples, "region"].astype(str)))
        test_regions = sorted(set(metadata.loc[test_samples, "region"].astype(str)))
        fold_rows.append(
            {
                "fold_no": fold_no,
                "heldout_id": fold_id,
                "n_train_samples": int(len(train_samples)),
                "n_test_samples": int(len(test_samples)),
                "n_train_regions": int(len(train_regions)),
                "n_test_regions": int(len(test_regions)),
                "n_test_regions_absent_from_train": int(len(set(test_regions) - set(train_regions))),
            }
        )
        fit, _ = fit_linear_projector(logcpm[train_samples], vsd[train_samples])
        projected = apply_projector(fit, logcpm[test_samples])
        native_reference = build_centroids(vsd, metadata, train_samples, train_regions)
        logcpm_reference = build_centroids(logcpm, metadata, train_samples, train_regions)

        for sample in test_samples:
            truth = str(metadata.loc[sample, "region"])
            if truth not in train_regions:
                for route in ["native_vsd", "projected_vsd", "logcpm_baseline"]:
                    rows.append(
                        abstain_row(
                            fold_no,
                            sample,
                            route,
                            truth,
                            "truth_region_absent_from_training_reference",
                            len(selected_genes),
                        )
                    )
                continue
            rows.append(
                valid_rank_row(
                    fold_no,
                    sample,
                    "native_vsd",
                    truth,
                    train_regions,
                    corr_scores(native_reference, vsd[sample].to_numpy(dtype=np.float32)),
                    len(selected_genes),
                )
            )
            rows.append(
                valid_rank_row(
                    fold_no,
                    sample,
                    "projected_vsd",
                    truth,
                    train_regions,
                    corr_scores(native_reference, projected[sample].to_numpy(dtype=np.float32)),
                    len(selected_genes),
                )
            )
            rows.append(
                valid_rank_row(
                    fold_no,
                    sample,
                    "logcpm_baseline",
                    truth,
                    train_regions,
                    corr_scores(logcpm_reference, logcpm[sample].to_numpy(dtype=np.float32)),
                    len(selected_genes),
                )
            )

    detail = pd.DataFrame(rows)
    prefix = f"bo2023_projected_vsd_exact_region_{args.split}"
    detail.to_csv(args.outdir / f"{prefix}_detail.csv", index=False)
    folds_df = pd.DataFrame(fold_rows)
    folds_df.to_csv(args.outdir / f"{prefix}_folds.csv", index=False)
    summary_df = route_summary(detail)
    summary_df.to_csv(args.outdir / f"{prefix}_route_summary.csv", index=False)
    write_json(
        args.outdir / f"{prefix}_summary.json",
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "split": args.split,
            "endpoint": "exact_region",
            "region_col": args.region_col,
            "gene_panel": args.gene_panel,
            "n_eval_genes": int(len(selected_genes)),
            "n_folds": int(len(folds)),
            "n_samples": int(len(samples)),
            "n_regions_total": int(metadata.loc[samples, "region"].nunique()),
            "routes": summary_df.to_dict(orient="records"),
            "folds": folds_df.to_dict(orient="records"),
        },
    )
    print(f"Wrote exact-region {args.split.upper()} outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
