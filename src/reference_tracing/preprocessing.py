from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_expression_file(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def prepare_expression_matrix(
    expr: pd.DataFrame,
    gene_col: str,
    expression_type: str = "raw counts",
) -> tuple[pd.DataFrame, dict]:
    if gene_col not in expr.columns:
        raise ValueError(f"Gene column not found: {gene_col}")
    df = expr.copy()
    df[gene_col] = df[gene_col].astype(str).str.strip().str.upper()
    sample_cols = [c for c in df.columns if c != gene_col]
    if not sample_cols:
        raise ValueError("Expression matrix must contain at least one sample column.")
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=[gene_col])
    df = df[df[gene_col].str.len() > 0]
    df = df.groupby(gene_col, as_index=True)[sample_cols].median()
    raw_gene_count = df.shape[0]
    expression_type = expression_type.lower()
    if expression_type in {"raw counts", "counts", "raw"}:
        library = df.sum(axis=0).replace(0, np.nan)
        df = df.div(library, axis=1) * 1_000_000
        transform = "raw counts converted to CPM then log1p"
    else:
        transform = f"{expression_type} transformed with log1p"
    df = np.log1p(df.fillna(0))
    df.index.name = "gene_symbol"
    meta = {
        "sample_number": len(sample_cols),
        "gene_number": raw_gene_count,
        "detected_genes": int((df > 0).any(axis=1).sum()),
        "expression_transform": transform,
    }
    return df, meta


def zscore_rows(values: pd.DataFrame) -> pd.DataFrame:
    mean = values.mean(axis=1)
    std = values.std(axis=1).replace(0, np.nan)
    return values.sub(mean, axis=0).div(std, axis=0).fillna(0)


def percentile_series(values: pd.Series) -> pd.Series:
    if len(values) < 5:
        return pd.Series(np.nan, index=values.index)
    return values.rank(pct=True) * 100


def level_from_percentile(percentile: float | int | None, score: float | None = None) -> str:
    if percentile is None or pd.isna(percentile):
        if score is None or pd.isna(score):
            return "Low"
        if score >= 2.5:
            return "High"
        if score >= 1.0:
            return "Moderate"
        return "Low"
    if percentile >= 90:
        return "High"
    if percentile >= 75:
        return "Moderate"
    return "Low"
