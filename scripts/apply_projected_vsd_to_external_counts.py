#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import compute_logcpm, load_projector_npz, apply_projector, write_json  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402


def read_external_counts(path: Path, gene_col: str | None) -> pd.DataFrame:
    sep = "," if path.suffix == ".csv" or path.name.endswith(".csv.gz") else "\t"
    df = pd.read_csv(path, sep=sep, compression="infer")
    if gene_col is None:
        gene_col = df.columns[0]
    if gene_col not in df.columns:
        raise ValueError(f"gene column {gene_col!r} not found in {path}")
    df[gene_col] = df[gene_col].astype(str).str.strip()
    value_cols = [col for col in df.columns if col != gene_col]
    matrix = df.set_index(gene_col)[value_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix = matrix.groupby(matrix.index.astype(str), sort=True).mean()
    matrix.index.name = "gene_symbol"
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description="Apply a Bo2023-trained projected-VSD projector to external raw counts.")
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--projector", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--gene-col", default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--run-network-tracing", action="store_true")
    parser.add_argument("--network-model", type=Path, default=None)
    parser.add_argument("--network-metadata", type=Path, default=None)
    parser.add_argument("--pairwise-rescue-model", type=Path, default=None)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    counts = read_external_counts(args.counts, args.gene_col)
    fit = load_projector_npz(args.projector)
    projected = apply_projector(fit, compute_logcpm(counts))
    projected_path = args.outdir / f"external_projected_vsd_{args.dataset}_matrix.tsv.gz"
    projected.to_csv(projected_path, sep="\t", compression="gzip")

    details = []
    if args.run_network_tracing:
        for sample_id in projected.columns.astype(str):
            expression = pd.DataFrame(
                {"gene_symbol": projected.index.astype(str), "tpm_value": projected[sample_id].to_numpy()}
            )
            trace_kwargs = {}
            if args.network_model is not None:
                trace_kwargs["model_path"] = args.network_model
            if args.network_metadata is not None:
                trace_kwargs["metadata_path"] = args.network_metadata
            if args.pairwise_rescue_model is not None:
                trace_kwargs["pairwise_rescue_path"] = args.pairwise_rescue_model
            traced = trace_network_expression(expression, **trace_kwargs)
            meta = traced.get("meta", {})
            for row in traced.get("results", [])[:3]:
                details.append({"sample_id": sample_id, **row, **{"n_overlap_genes": meta.get("n_overlap_genes")}})
        pd.DataFrame(details).to_csv(args.outdir / f"external_projected_vsd_{args.dataset}_detail.csv", index=False)

    overlap = int(pd.Index(fit.genes).isin(counts.index.astype(str)).sum())
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": args.dataset,
        "counts_path": str(args.counts),
        "projector_path": str(args.projector),
        "projected_matrix_path": str(projected_path),
        "n_input_genes": int(counts.shape[0]),
        "n_projector_genes": int(len(fit.genes)),
        "n_overlap_projector_genes": overlap,
        "overlap_fraction": float(overlap / max(len(fit.genes), 1)),
        "n_samples": int(counts.shape[1]),
        "interpretation_boundary": "cross-domain projected-space analysis; not native Bo2023 VSD",
    }
    write_json(args.outdir / f"external_projected_vsd_{args.dataset}_summary.json", summary)
    print(f"Wrote projected matrix to {projected_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
