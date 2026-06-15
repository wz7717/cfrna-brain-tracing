#!/usr/bin/env python
"""Inspect downloaded reference atlas files without loading large matrices into memory."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference_atlases"
OUT = BASE / "00_manifest" / "reference_file_inspection.tsv"

TEXT_SUFFIXES = (".csv", ".tsv", ".txt", ".gct")


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def logical_name(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name


def count_delimited(path: Path, sep: str, skiprows: int = 0) -> tuple[int, int, list[str]]:
    row_count = 0
    header: list[str] = []
    with open_text(path) as handle:
        for _ in range(skiprows):
            next(handle, None)
        reader = csv.reader(handle, delimiter=sep)
        header = next(reader, [])
        for row_count, _ in enumerate(reader, start=1):
            pass
    return row_count, len(header), header[:10]


def inspect_table(path: Path) -> dict[str, str]:
    name = logical_name(path)
    skiprows = 2 if name.endswith(".gct") else 0
    sep = "," if name.endswith(".csv") else "\t"
    rows, cols, first_cols = count_delimited(path, sep, skiprows=skiprows)
    columns_lower = [c.lower() for c in first_cols]
    all_cols = []
    try:
        preview = pd.read_csv(path, sep=sep, skiprows=skiprows, nrows=0)
        all_cols = [str(c).lower() for c in preview.columns]
    except Exception:
        all_cols = columns_lower
    joined = " ".join(all_cols)
    return {
        "file_path": str(path),
        "file_type": "table",
        "file_size_bytes": str(path.stat().st_size),
        "row_count": str(rows),
        "column_count": str(cols),
        "first_10_columns": "|".join(first_cols),
        "whether_gene_column_exists": str(any(c in {"gene", "gene_id", "gene name", "gene_name", "name", "description"} for c in all_cols)),
        "whether_ensembl_id_exists": str("ensembl" in joined or any(c.startswith("ensg") for c in all_cols)),
        "whether_gene_symbol_exists": str(any(x in joined for x in ["gene symbol", "gene_symbol", "gene name", "gene_name", "symbol"])),
        "h5ad_can_open_backed": "",
        "n_obs": "",
        "n_vars": "",
        "obs_columns": "",
        "var_columns": "",
        "error": "",
    }


def inspect_h5ad(path: Path) -> dict[str, str]:
    row = {
        "file_path": str(path),
        "file_type": "h5ad",
        "file_size_bytes": str(path.stat().st_size),
        "row_count": "",
        "column_count": "",
        "first_10_columns": "",
        "whether_gene_column_exists": "",
        "whether_ensembl_id_exists": "",
        "whether_gene_symbol_exists": "",
        "h5ad_can_open_backed": "False",
        "n_obs": "",
        "n_vars": "",
        "obs_columns": "",
        "var_columns": "",
        "error": "",
    }
    try:
        import anndata as ad

        data = ad.read_h5ad(path, backed="r")
        row["h5ad_can_open_backed"] = "True"
        row["n_obs"] = str(data.n_obs)
        row["n_vars"] = str(data.n_vars)
        row["obs_columns"] = "|".join(map(str, list(data.obs.columns)[:50]))
        row["var_columns"] = "|".join(map(str, list(data.var.columns)[:50]))
        data.file.close()
    except Exception as exc:  # noqa: BLE001
        row["error"] = str(exc)
    return row


def inspect_file(path: Path) -> dict[str, str] | None:
    name = logical_name(path)
    if name.endswith(TEXT_SUFFIXES):
        try:
            return inspect_table(path)
        except Exception as exc:  # noqa: BLE001
            return {
                "file_path": str(path),
                "file_type": "table",
                "file_size_bytes": str(path.stat().st_size),
                "row_count": "",
                "column_count": "",
                "first_10_columns": "",
                "whether_gene_column_exists": "",
                "whether_ensembl_id_exists": "",
                "whether_gene_symbol_exists": "",
                "h5ad_can_open_backed": "",
                "n_obs": "",
                "n_vars": "",
                "obs_columns": "",
                "var_columns": "",
                "error": str(exc),
            }
    if name.endswith(".h5ad"):
        return inspect_h5ad(path)
    return None


def main() -> int:
    BASE.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in BASE.rglob("*"):
        if not path.is_file():
            continue
        if ".part" in path.suffixes:
            continue
        result = inspect_file(path)
        if result:
            rows.append(result)
    columns = [
        "file_path",
        "file_type",
        "file_size_bytes",
        "row_count",
        "column_count",
        "first_10_columns",
        "whether_gene_column_exists",
        "whether_ensembl_id_exists",
        "whether_gene_symbol_exists",
        "h5ad_can_open_backed",
        "n_obs",
        "n_vars",
        "obs_columns",
        "var_columns",
        "error",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT} ({len(rows)} inspected files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
