from __future__ import annotations

import sqlite3
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from signature_builder import get_latest_sigset_id as _sig_latest
from .models import apply_value_transform

# ---------------------------------------------------------------------------
# 安全白名单：仅允许对已知表执行 PRAGMA / 格式化查询
# ---------------------------------------------------------------------------
_ALLOWED_TABLES = frozenset({
    'reference_expression',
    'region_gene_signature',
    'cfrna_expression',
    'cfrna_samples',
    'macaque_brain_atlas',
    'atlas_versions',
    'signature_sets',
    'signature_genes',
    'analysis_runs',
    'analysis_results',
    'source_tracing_results',
    'analysis_history',
    'disease_associations',
    'sample_qc',
})


def _validate_table_name(name: str) -> str:
    """校验表名在白名单内，防止 PRAGMA 注入，返回原始表名。"""
    if name not in _ALLOWED_TABLES:
        raise ValueError(f"Unknown or disallowed table name: {name!r}")
    return name


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    cur = conn.cursor()
    cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (name,),
    )
    return cur.fetchone() is not None


def _col_exists(conn: sqlite3.Connection, table: str, col: str) -> bool:
    _validate_table_name(table)
    cur = conn.cursor()
    cur.execute(f"PRAGMA table_info({table})")
    return col in {row[1] for row in cur.fetchall()}


def load_signature_genes(db_path: str, sigset_id: int) -> Optional[List[str]]:
    conn = sqlite3.connect(db_path)
    try:
        if not _table_exists(conn, 'signature_genes'):
            return None
        cur = conn.cursor()
        cur.execute(
            'SELECT DISTINCT gene_symbol FROM signature_genes WHERE sigset_id=?',
            (int(sigset_id),),
        )
        rows = cur.fetchall()
        return [r[0] for r in rows] if rows else None
    finally:
        conn.close()


def get_latest_sigset_id(db_path: str, atlas_id: int = 1) -> Optional[int]:
    return _sig_latest(db_path, atlas_id)


def load_marker_signature_genes(
    db_path: str, topk_per_region: int = 200
) -> Optional[List[str]]:
    conn = sqlite3.connect(db_path)
    try:
        if not _table_exists(conn, 'region_gene_signature'):
            return None
        has_atlas = _col_exists(conn, 'region_gene_signature', 'atlas_id')
        cols = (
            'region_id, gene_symbol, specificity_score, is_marker, atlas_id'
            if has_atlas
            else 'region_id, gene_symbol, specificity_score, is_marker'
        )
        df = pd.read_sql_query(
            f"SELECT {cols} FROM region_gene_signature WHERE is_marker = 1",
            conn,
        )
        if df.empty:
            return None
        if 'atlas_id' in df.columns:
            df = df[
                (df['atlas_id'].isna())
                | (df['atlas_id'].astype(float).astype(int) == 1)
            ]
        df = (
            df.sort_values(
                ['region_id', 'specificity_score'], ascending=[True, False]
            )
            .groupby('region_id')
            .head(int(topk_per_region))
        )
        return sorted(df['gene_symbol'].dropna().astype(str).unique().tolist())
    finally:
        conn.close()


def load_reference_matrix(
    db_path: str,
    cache: Dict,
    atlas_id: int,
    sigset_id: Optional[int],
    use_value: str,
    fallback_marker_topk: int = 200,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    cache_key = (
        int(atlas_id),
        int(sigset_id) if sigset_id is not None else None,
        use_value,
        int(fallback_marker_topk),
    )
    if cache_key in cache:
        return cache[cache_key]

    genes_filter = (
        load_signature_genes(db_path, sigset_id)
        if sigset_id is not None
        else None
    )
    if genes_filter is None:
        genes_filter = load_marker_signature_genes(
            db_path, topk_per_region=fallback_marker_topk
        )

    conn = sqlite3.connect(db_path)
    try:
        has_atlas = _col_exists(conn, 'reference_expression', 'atlas_id')
        select_sql = (
            'SELECT region_id, gene_symbol, avg_tpm, atlas_id FROM reference_expression'
            if has_atlas
            else 'SELECT region_id, gene_symbol, avg_tpm FROM reference_expression'
        )
        ref = pd.read_sql_query(select_sql, conn)
        if ref.empty:
            result = (np.array([], dtype=object), np.zeros((0, 0)), [])
            cache[cache_key] = result
            return result

        ref = ref.dropna(subset=['region_id', 'gene_symbol', 'avg_tpm'])
        ref['gene_symbol'] = ref['gene_symbol'].astype(str)
        if 'atlas_id' in ref.columns:
            ref = ref[
                (ref['atlas_id'].isna())
                | (ref['atlas_id'].astype(float).astype(int) == int(atlas_id))
            ]
        if genes_filter is not None:
            ref = ref[ref['gene_symbol'].isin(set(genes_filter))]

        mat = (
            ref.pivot_table(
                index='gene_symbol',
                columns='region_id',
                values='avg_tpm',
                aggfunc='mean',
            )
            .fillna(0.0)
        )
        mat = mat.loc[mat.sum(axis=1) > 0]
        regions = list(mat.columns)
        genes = mat.index.values.astype(object)
        A = apply_value_transform(mat.values.astype(float), use_value)
        result = (genes, A, regions)
        cache[cache_key] = result
        return result
    finally:
        conn.close()


def load_sample_vector(
    db_path: str,
    sample_id: str,
    genes: np.ndarray,
    use_value: str,
) -> np.ndarray:
    if genes.size == 0:
        return np.array([], dtype=float)

    conn = sqlite3.connect(db_path)
    try:
        col = 'tpm_value'
        if use_value == 'zscore' and _col_exists(conn, 'cfrna_expression', 'zscore_tpm'):
            col = 'zscore_tpm'
        elif use_value in ('log1p', 'zscore') and _col_exists(
            conn, 'cfrna_expression', 'log_tpm'
        ):
            col = 'log_tpm'

        df = pd.read_sql_query(
            f'SELECT gene_symbol, {col} AS value FROM cfrna_expression WHERE sample_id = ?',
            conn,
            params=[sample_id],
        )
        if df.empty:
            return np.array([], dtype=float)

        df = df.dropna(subset=['gene_symbol', 'value'])
        df['gene_symbol'] = df['gene_symbol'].astype(str)
        s = df.groupby('gene_symbol')['value'].mean()
        vals = (
            s.reindex(genes)
            .astype(float)
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
            .values
        )
        if col == 'tpm_value':
            vals = apply_value_transform(vals, use_value)
        elif col == 'log_tpm' and use_value == 'zscore':
            vals = (vals - vals.mean()) / (vals.std() + 1e-8)
        return vals
    finally:
        conn.close()
