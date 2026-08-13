from __future__ import annotations

import io
from pathlib import Path
from typing import BinaryIO

import pandas as pd


GENE_COLUMN_CANDIDATES = ("gene_symbol", "gene", "symbol")
VALUE_COLUMN_GROUPS = (
    ("raw_counts", ("raw_counts", "raw_count", "counts", "count", "read_count", "readcount", "reads")),
    ("logcpm", ("logcpm", "log_cpm", "log2cpm", "log2_cpm")),
    ("logtpm_fallback", ("logtpm", "log_tpm", "log1p_tpm")),
    ("tpm_fallback", ("tpm_value", "tpm", "expression", "value")),
)


def normalize_query_expression(frame: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Normalize one uploaded expression table for the locked production route."""
    lower_map = {str(column).strip().lower(): column for column in frame.columns}
    gene_col = next((lower_map[name] for name in GENE_COLUMN_CANDIDATES if name in lower_map), None)
    value_col = None
    query_source = ""
    for source, candidates in VALUE_COLUMN_GROUPS:
        value_col = next((lower_map[name] for name in candidates if name in lower_map), None)
        if value_col is not None:
            query_source = source
            break
    if gene_col is None or value_col is None:
        raise ValueError(
            "Input must include gene_symbol/gene/symbol and one expression column: "
            "raw counts, logCPM, logTPM or TPM."
        )

    out = frame[[gene_col, value_col]].copy()
    out.columns = ["gene_symbol", "query_value"]
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip()
    out["query_value"] = pd.to_numeric(out["query_value"], errors="coerce")
    out = out.dropna(subset=["gene_symbol", "query_value"])
    out = out[out["gene_symbol"] != ""]
    if out.empty:
        raise ValueError("No valid expression rows were found.")
    out = out.groupby("gene_symbol", as_index=False)["query_value"].mean()

    if query_source == "raw_counts":
        out["read_count"] = out["query_value"].clip(lower=0)
        if float(out["read_count"].sum()) <= 0:
            raise ValueError("Raw counts must sum to a positive value.")
        return out[["gene_symbol", "read_count"]], query_source
    if query_source in {"logcpm", "logtpm_fallback"}:
        out["log_tpm"] = out["query_value"]
        return out[["gene_symbol", "log_tpm"]], query_source

    out["tpm_value"] = out["query_value"].clip(lower=0)
    return out[["gene_symbol", "tpm_value"]], query_source


def read_expression_file(
    source: str | Path | BinaryIO,
    *,
    filename: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """Read CSV/TSV/TXT/XLSX input and apply the shared query normalization."""
    if isinstance(source, (str, Path)):
        path = Path(source)
        name = path.name.lower()
        if name.endswith(".xlsx"):
            frame = pd.read_excel(path)
        else:
            raw = path.read_bytes()
            separator = "\t" if name.endswith((".tsv", ".txt")) or raw[:2048].count(b"\t") > raw[:2048].count(b",") else ","
            frame = pd.read_csv(io.BytesIO(raw), sep=separator)
    else:
        name = str(filename or getattr(source, "name", "input.csv")).lower()
        if name.endswith(".xlsx"):
            frame = pd.read_excel(source)
        else:
            raw = source.getvalue() if hasattr(source, "getvalue") else source.read()
            if isinstance(raw, str):
                raw = raw.encode("utf-8")
            separator = "\t" if name.endswith((".tsv", ".txt")) or raw[:2048].count(b"\t") > raw[:2048].count(b",") else ","
            frame = pd.read_csv(io.BytesIO(raw), sep=separator)
    return normalize_query_expression(frame)
