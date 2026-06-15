#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402


DEFAULT_MATRIX = ROOT / "data" / "tcga_brain_tumor_expression" / "tcga_gbm_lgg_primary_tumor_project_mean_tpm.tsv"
DEFAULT_OUTDIR = ROOT / "results" / "tcga_gbm_lgg_project_centroid_tracing_20260602"
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"


def expression_frame(matrix: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": matrix.index.astype(str),
            "tpm_value": pd.to_numeric(matrix[sample_id], errors="coerce").fillna(0.0).to_numpy(),
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Score TCGA-GBM/LGG project-mean TPM profiles with formal tracing route.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--atlas-id", type=int, default=4)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    meta: dict[str, Any] = {}
    for sample_id in matrix.columns.astype(str).tolist():
        expr = expression_frame(matrix, sample_id)
        network_out = trace_network_expression(expr, min_overlap_fraction=0.20)
        region_out = trace_bo2023_secondary_regions(expr, network_out, str(args.db), int(args.atlas_id), topk=15)
        meta[sample_id] = {"network": network_out.get("meta", {}), "region": region_out.get("meta", {})}
        for row in network_out.get("results", [])[:10]:
            network_rows.append({"sample_id": sample_id, **row})
        for row in region_out.get("results", [])[:15]:
            region_rows.append({"sample_id": sample_id, **row})

    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    network_df.to_csv(args.outdir / "tcga_gbm_lgg_project_centroid_network_tracing.csv", index=False, encoding="utf-8-sig")
    region_df.to_csv(args.outdir / "tcga_gbm_lgg_project_centroid_region_tracing.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "tcga_gbm_lgg_project_centroid_tracing_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    summary_rows = []
    for sample_id in matrix.columns.astype(str).tolist():
        net_top = network_df[network_df["sample_id"].eq(sample_id)].sort_values("rank").head(3)
        reg_top = region_df[region_df["sample_id"].eq(sample_id)].sort_values("rank").head(5)
        summary_rows.append(
            {
                "sample_id": sample_id,
                "network_top3": " | ".join(net_top["network_id"].astype(str).tolist()),
                "network_top1": str(net_top.iloc[0]["network_id"]) if len(net_top) else "",
                "region_top5": " | ".join(reg_top["region_id"].astype(str).tolist()),
                "region_top1": str(reg_top.iloc[0]["region_id"]) if len(reg_top) else "",
                "network_overlap_genes": meta[sample_id]["network"].get("n_overlap_genes"),
                "network_overlap_fraction": meta[sample_id]["network"].get("overlap_fraction"),
                "region_overlap_genes": meta[sample_id]["region"].get("n_overlap_genes"),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.outdir / "tcga_gbm_lgg_project_centroid_tracing_summary.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
