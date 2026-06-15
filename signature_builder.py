from __future__ import annotations

import json
import sqlite3
from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd

DEFAULT_HOUSEKEEPING = {
    "ACTB", "GAPDH", "RPLP0", "RPL13A", "B2M",
    "HPRT1", "EEF1A1", "RPS18", "RPL19", "PPIA",
}
DEFAULT_BLOOD_BACKGROUND = {
    "HBB", "HBA1", "HBA2", "ALAS2", "GYPA",
    "PTPRC", "LYZ", "S100A8", "S100A9", "LCN2",
}


def ensure_signature_tables(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS signature_sets (
            sigset_id              INTEGER PRIMARY KEY AUTOINCREMENT,
            atlas_id               INTEGER,
            method                 TEXT NOT NULL,
            topk_per_region        INTEGER NOT NULL,
            remove_housekeeping    INTEGER DEFAULT 1,
            remove_blood_background INTEGER DEFAULT 1,
            created_at             TEXT DEFAULT (datetime('now')),
            params_json            TEXT
        )
    """)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS signature_genes (
            sigset_id   INTEGER NOT NULL,
            region_id   TEXT NOT NULL,
            gene_symbol TEXT NOT NULL,
            weight      REAL DEFAULT 1.0,
            PRIMARY KEY(sigset_id, region_id, gene_symbol)
        )
    """)
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_sig_genes "
        "ON signature_genes(sigset_id, region_id, gene_symbol)"
    )
    conn.commit()


def _has_col(conn: sqlite3.Connection, table: str, col: str) -> bool:
    cur = conn.cursor()
    cur.execute(f"PRAGMA table_info({table})")
    return col in {r[1] for r in cur.fetchall()}


def _compute_region_scores(
    genes: List[str],
    regions: List[str],
    X: np.ndarray,
    method: str,
) -> np.ndarray:
    """对每个 region 计算 gene 的特异性得分矩阵，返回 shape=(n_genes, n_regions)。"""
    n_regions = len(regions)
    score_mat = np.zeros((len(genes), n_regions), dtype=float)
    for j in range(n_regions):
        in_r = X[:, j]
        rest_mat = np.delete(X, j, axis=1)
        rest_mean = rest_mat.mean(axis=1)
        rest_std = rest_mat.std(axis=1) + 1e-8
        fold = in_r - rest_mean
        eff = (in_r - rest_mean) / rest_std
        if method in ('foldchange', 'fc'):
            score_mat[:, j] = fold
        elif method in ('effect_size', 'es'):
            score_mat[:, j] = eff
        else:  # hybrid_specificity (default)
            score_mat[:, j] = 0.65 * eff + 0.35 * fold
    return score_mat


def _select_signature_genes(
    genes: List[str],
    regions: List[str],
    score_mat: np.ndarray,
    X: np.ndarray,
    exclude: set,
    topk_per_region: int,
    min_log_expr_in_region: float,
    max_global_repeat: int,
) -> List[Tuple[str, str, float]]:
    """从得分矩阵中按策略挑选签名基因，返回 (region, gene, score) 三元组列表。"""
    sig_rows: List[Tuple[str, str, float]] = []
    global_counts: Dict[str, int] = {}

    for j, region in enumerate(regions):
        in_r = X[:, j]
        scores = score_mat[:, j]
        picked = 0
        for idx in np.argsort(scores)[::-1]:
            g = genes[int(idx)]
            if (
                g in exclude
                or in_r[idx] < float(min_log_expr_in_region)
                or global_counts.get(g, 0) >= int(max_global_repeat)
            ):
                continue
            sig_rows.append((region, g, float(scores[idx])))
            global_counts[g] = global_counts.get(g, 0) + 1
            picked += 1
            if picked >= int(topk_per_region):
                break
    return sig_rows


def _persist_signature_set(
    conn: sqlite3.Connection,
    atlas_id: int,
    method: str,
    topk_per_region: int,
    remove_housekeeping: bool,
    remove_blood_background: bool,
    sig_rows: List[Tuple[str, str, float]],
    params: Dict,
    min_log_expr_in_region: float,
    max_global_repeat: int,
) -> int:
    """将签名集写入数据库，返回新生成的 sigset_id。"""
    cur = conn.cursor()
    cur.execute(
        "INSERT INTO signature_sets"
        "(atlas_id, method, topk_per_region, remove_housekeeping,"
        " remove_blood_background, params_json)"
        " VALUES (?, ?, ?, ?, ?, ?)",
        (
            int(atlas_id),
            method,
            int(topk_per_region),
            1 if remove_housekeeping else 0,
            1 if remove_blood_background else 0,
            json.dumps(
                {
                    'min_log_expr_in_region': float(min_log_expr_in_region),
                    'max_global_repeat': int(max_global_repeat),
                    **params,
                },
                ensure_ascii=False,
            ),
        ),
    )
    sigset_id = int(cur.lastrowid)
    cur.executemany(
        "INSERT OR REPLACE INTO signature_genes"
        "(sigset_id, region_id, gene_symbol, weight)"
        " VALUES (?, ?, ?, ?)",
        [(sigset_id, r, g, w) for (r, g, w) in sig_rows],
    )
    conn.commit()
    return sigset_id


def build_signature_set(
    db_path: str,
    atlas_id: int = 1,
    method: str = 'hybrid_specificity',
    topk_per_region: int = 120,
    remove_housekeeping: bool = True,
    remove_blood_background: bool = True,
    extra_exclude_genes: Optional[Iterable[str]] = None,
    min_log_expr_in_region: float = 0.10,
    max_global_repeat: int = 3,
    params: Optional[Dict] = None,
) -> int:
    """构建并持久化一套区域特异性签名基因集。

    Args:
        db_path: SQLite 数据库路径。
        atlas_id: 图谱版本 ID。
        method: 得分策略，支持 'hybrid_specificity'（默认）、'foldchange'、'effect_size'。
        topk_per_region: 每个脑区最多保留的签名基因数。
        remove_housekeeping: 是否排除管家基因。
        remove_blood_background: 是否排除血液背景基因。
        extra_exclude_genes: 额外需要排除的基因集合。
        min_log_expr_in_region: 基因在该区域的最低 log 表达阈值。
        max_global_repeat: 同一基因最多出现在多少个区域的签名中。
        params: 额外参数，记录在 params_json 字段中。

    Returns:
        新创建的 sigset_id。
    """
    params = params or {}

    # 构建排除集
    exclude: set = set(extra_exclude_genes or [])
    if remove_housekeeping:
        exclude |= DEFAULT_HOUSEKEEPING
    if remove_blood_background:
        exclude |= DEFAULT_BLOOD_BACKGROUND

    conn = sqlite3.connect(db_path)
    try:
        ensure_signature_tables(conn)

        # 读取参考表达谱
        select_sql = (
            'SELECT region_id, gene_symbol, avg_tpm, atlas_id FROM reference_expression'
            if _has_col(conn, 'reference_expression', 'atlas_id')
            else 'SELECT region_id, gene_symbol, avg_tpm FROM reference_expression'
        )
        ref = pd.read_sql_query(select_sql, conn)
        if ref.empty:
            raise ValueError(
                'reference_expression is empty; cannot build signature set.'
            )
        if 'atlas_id' in ref.columns:
            ref = ref[
                (ref['atlas_id'].isna())
                | (ref['atlas_id'].astype(float).astype(int) == int(atlas_id))
            ]

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
        genes = mat.index.astype(str).tolist()
        regions = list(mat.columns)
        X = np.log1p(np.clip(mat.values.astype(float), 0, None))

        # 计算得分并挑选签名基因
        score_mat = _compute_region_scores(genes, regions, X, method)
        sig_rows = _select_signature_genes(
            genes, regions, score_mat, X,
            exclude, topk_per_region,
            min_log_expr_in_region, max_global_repeat,
        )

        # 持久化
        sigset_id = _persist_signature_set(
            conn, atlas_id, method, topk_per_region,
            remove_housekeeping, remove_blood_background,
            sig_rows, params, min_log_expr_in_region, max_global_repeat,
        )
        return sigset_id
    finally:
        conn.close()


def get_latest_sigset_id(db_path: str, atlas_id: int = 1) -> Optional[int]:
    """获取指定 atlas 最新的签名集 ID，不存在时返回 None。"""
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='signature_sets'"
        )
        if cur.fetchone() is None:
            return None
        cur.execute(
            "SELECT sigset_id FROM signature_sets "
            "WHERE atlas_id = ? OR atlas_id IS NULL "
            "ORDER BY sigset_id DESC LIMIT 1",
            (int(atlas_id),),
        )
        row = cur.fetchone()
        if row:
            return int(row[0])
        # fallback：不过滤 atlas_id
        cur.execute(
            "SELECT sigset_id FROM signature_sets ORDER BY sigset_id DESC LIMIT 1"
        )
        row = cur.fetchone()
        return int(row[0]) if row else None
    finally:
        conn.close()
