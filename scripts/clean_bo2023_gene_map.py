#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import clean_excel_date_gene_symbol  # noqa: E402


DEFAULT_IN = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.csv"
DEFAULT_OUT = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv"


def main() -> int:
    parser = argparse.ArgumentParser(description="Clean Excel-mangled date-like gene symbols in the Bo2023 gene map.")
    parser.add_argument("--input", type=Path, default=DEFAULT_IN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    frame = pd.read_csv(args.input, dtype=str)
    if "Gene.name" not in frame.columns:
        raise ValueError("input gene map must contain Gene.name")
    original = frame["Gene.name"].fillna("").astype(str).str.strip()
    cleaned = original.map(clean_excel_date_gene_symbol)
    changes = frame.loc[original.ne(cleaned), ["Gene.stable.ID", "Gene.name"]].copy()
    changes["Gene.name.cleaned"] = cleaned[original.ne(cleaned)].to_numpy()
    frame["Gene.name"] = cleaned
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)

    change_path = args.output.with_suffix(".changes.csv")
    changes.to_csv(change_path, index=False)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(args.input),
        "output": str(args.output),
        "changes": str(change_path),
        "n_rows": int(len(frame)),
        "n_cleaned_symbols": int(len(changes)),
        "cleaned_symbols": changes.to_dict(orient="records"),
    }
    args.output.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote cleaned gene map to {args.output}")
    print(f"n_cleaned_symbols={len(changes)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
