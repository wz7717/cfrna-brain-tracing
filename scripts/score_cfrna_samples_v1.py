#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.reference_tracing.preprocessing import read_expression_file
from src.reference_tracing.scoring import run_reference_tracing


def main() -> int:
    parser = argparse.ArgumentParser(description="Score cfRNA samples against reference tracing v1.")
    parser.add_argument("--expr", required=True, help="Input TSV/CSV expression matrix with genes in rows.")
    parser.add_argument("--gene-col", required=True, help="Gene symbol column name.")
    parser.add_argument("--outdir", default="results/cfrna_tracing_v1/", help="Output directory.")
    parser.add_argument("--expression-type", default="raw counts", choices=["raw counts", "TPM", "CPM", "FPKM", "counts", "raw"], help="Input expression scale.")
    parser.add_argument("--db", default=None, help="Optional reference SQLite path.")
    args = parser.parse_args()
    expr = read_expression_file(args.expr)
    run_reference_tracing(expr, args.gene_col, args.expression_type, args.outdir, db_path=args.db, make_plots=True)
    print(f"Wrote cfRNA tracing outputs to {Path(args.outdir).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
