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
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_GENE_MAP, DEFAULT_VSD  # noqa: E402
from scripts.build_bo2023_reference_projector import DEFAULT_MODEL_GENES, DEFAULT_SAMPLE_INFO, read_locked_model_genes  # noqa: E402


def corr_scores(reference: np.ndarray, sample: np.ndarray) -> np.ndarray:
    ref0 = reference - reference.mean(axis=0, keepdims=True)
    vec0 = sample - sample.mean()
    denom = np.sqrt(np.square(ref0).sum(axis=0) * np.square(vec0).sum() + 1e-12)
    return np.nan_to_num((ref0 * vec0[:, None]).sum(axis=0) / denom)


def read_labels(path: Path, sheet: str, network_col: str) -> pd.Series:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", network_col])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info[network_col] = info[network_col].astype(str).str.strip()
    return info.drop_duplicates("sample_id").set_index("sample_id")[network_col]


def build_centroids(values: pd.DataFrame, labels: pd.Series, train_samples: list[str], groups: list[str]) -> np.ndarray:
    cols = []
    for group in groups:
        group_samples = [sample for sample in train_samples if labels.loc[sample] == group]
        cols.append(values[group_samples].mean(axis=1).to_numpy(dtype=np.float32))
    return np.column_stack(cols).astype(np.float32)


def make_rank_row(
    fold_no: int,
    sample: str,
    route: str,
    truth: str,
    groups: list[str],
    scores: np.ndarray,
    n_genes: int,
) -> dict[str, object]:
    order = np.argsort(scores)[::-1]
    ranked = [groups[int(i)] for i in order]
    true_rank = ranked.index(truth) + 1 if truth in ranked else len(ranked) + 1
    padded = ranked[:3] + [""] * max(0, 3 - len(ranked))
    return {
        "fold_no": fold_no,
        "sample_id": sample,
        "route": route,
        "label": truth,
        "pred_top1": padded[0],
        "pred_top2": padded[1],
        "pred_top3": padded[2],
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "true_rank": int(true_rank),
        "score_top1": float(scores[int(order[0])]),
        "score_top2": float(scores[int(order[1])]) if len(order) > 1 else float("nan"),
        "decision_margin": float(scores[int(order[0])] - scores[int(order[1])]) if len(order) > 1 else float("nan"),
        "n_overlap_genes": int(n_genes),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Fold-local LOSO smoke/formal validation for projected Bo2023-like VSD.")
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--network-col", default="SaleemNetworks")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_MODEL_GENES)
    parser.add_argument("--max-folds", type=int, default=0, help="Use 0 for all folds.")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    labels = read_labels(args.sample_info, args.sample_sheet, args.network_col)
    samples = [sample for sample in samples if sample in labels.index and labels.loc[sample] not in {"", "nan", "None"}]
    locked_genes = read_locked_model_genes(args.locked_model_genes)
    selected_genes = [gene for gene in locked_genes if gene in set(genes)]
    if len(selected_genes) < 20:
        selected_genes = genes
    counts = counts.loc[selected_genes, samples]
    vsd = vsd.loc[selected_genes, samples]
    logcpm = compute_logcpm(counts)
    max_folds = len(samples) if args.max_folds == 0 else min(args.max_folds, len(samples))

    rows = []
    for fold_no, sample in enumerate(samples[:max_folds], start=1):
        train_samples = [s for s in samples if s != sample]
        truth = str(labels.loc[sample])
        groups = sorted(set(labels.loc[train_samples].astype(str)))
        if truth not in groups:
            continue
        fit, _ = fit_linear_projector(logcpm[train_samples], vsd[train_samples])
        projected = apply_projector(fit, logcpm[[sample]])
        native_reference = build_centroids(vsd, labels, train_samples, groups)
        logcpm_reference = build_centroids(logcpm, labels, train_samples, groups)
        rows.append(make_rank_row(fold_no, sample, "native_vsd", truth, groups, corr_scores(native_reference, vsd[sample].to_numpy(dtype=np.float32)), len(selected_genes)))
        rows.append(make_rank_row(fold_no, sample, "projected_vsd", truth, groups, corr_scores(native_reference, projected[sample].to_numpy(dtype=np.float32)), len(selected_genes)))
        rows.append(make_rank_row(fold_no, sample, "logcpm_baseline", truth, groups, corr_scores(logcpm_reference, logcpm[sample].to_numpy(dtype=np.float32)), len(selected_genes)))
    detail = pd.DataFrame(rows)
    detail.to_csv(args.outdir / "bo2023_projected_vsd_loso_detail.csv", index=False)
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
    route_summary.to_csv(args.outdir / "bo2023_projected_vsd_loso_route_summary.csv", index=False)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "fold-local SaleemNetworks LOSO using locked network-model genes when available",
        "n_folds": int(max_folds),
        "n_common_samples": int(len(samples)),
        "n_common_genes": int(len(genes)),
        "n_eval_genes": int(len(selected_genes)),
        "routes": route_summary.to_dict(orient="records"),
    }
    write_json(args.outdir / "bo2023_projected_vsd_loso_summary.json", summary)
    print(f"Wrote LOSO outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
