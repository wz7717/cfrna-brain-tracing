from __future__ import annotations

import sqlite3

import pandas as pd


def table_exists(db_path: str, table: str) -> bool:
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table,))
        return cur.fetchone() is not None
    finally:
        conn.close()


def get_atlas_options(db_path: str, species_mode: str | None = None, include_legacy: bool = True):
    if not table_exists(db_path, "atlas_versions"):
        return [(1, "default atlas (legacy)")] if include_legacy else []
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            "SELECT atlas_id, atlas_name, species, level, build_version FROM atlas_versions ORDER BY atlas_id DESC",
            conn,
        )
    finally:
        conn.close()
    if df.empty:
        return [(1, "default atlas (legacy)")] if include_legacy else []

    if species_mode:
        mode = str(species_mode).lower().strip()
        if mode == "human":
            tokens = ["human", "homo sapiens"]
        else:
            tokens = ["rhesus", "macaca", "macaque", "macaca mulatta", "macaca fascicularis"]
        species_series = df["species"].fillna("").astype(str).str.lower()
        df = df[species_series.apply(lambda x: any(t in x for t in tokens))].copy()
    if df.empty:
        return []

    out = []
    for _, r in df.iterrows():
        out.append(
            (
                int(r["atlas_id"]),
                f"{int(r['atlas_id'])} | {r['atlas_name']} | {r.get('species', '')} {r.get('level', '')} {r.get('build_version', '')}".strip(),
            )
        )
    return out


def get_sigset_options(db_path: str, atlas_id: int):
    if not table_exists(db_path, "signature_sets"):
        return [(None, "no signature set available (recommended to build one first)")]
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            "SELECT sigset_id, atlas_id, method, topk_per_region, created_at FROM signature_sets WHERE atlas_id = ? OR atlas_id IS NULL ORDER BY sigset_id DESC",
            conn,
            params=[int(atlas_id)],
        )
    finally:
        conn.close()
    if df.empty:
        return [(None, "no signature set available (recommended to build one first)")]
    out = []
    for _, r in df.iterrows():
        aid = r.get("atlas_id")
        aid_s = "legacy" if pd.isna(aid) else str(int(aid))
        out.append(
            (
                int(r["sigset_id"]),
                f"{int(r['sigset_id'])} | atlas:{aid_s} | {r['method']} | topK={int(r['topk_per_region'])} | {r['created_at']}",
            )
        )
    out.append((None, "do not use a signature set (not recommended; discrimination will be weaker)"))
    return out


def get_atlas_metadata(db_path: str, atlas_id: int) -> dict:
    if not table_exists(db_path, "atlas_versions"):
        return {}
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            """
            SELECT atlas_id, atlas_name, species, level, build_version,
                   gene_id_type, normalization, notes
            FROM atlas_versions
            WHERE atlas_id = ?
            """,
            conn,
            params=[int(atlas_id)],
        )
    finally:
        conn.close()
    if df.empty:
        return {}
    return df.iloc[0].to_dict()


def list_runs(db_path: str, limit: int = 200) -> pd.DataFrame:
    if not table_exists(db_path, "analysis_runs"):
        return pd.DataFrame()
    conn = sqlite3.connect(db_path)
    try:
        return pd.read_sql_query(
            "SELECT run_id, sample_id, atlas_id, sigset_id, method, created_at FROM analysis_runs ORDER BY created_at DESC LIMIT ?",
            conn,
            params=[int(limit)],
        )
    finally:
        conn.close()


def get_run_results(db_path: str, run_id: str) -> pd.DataFrame:
    if not table_exists(db_path, "analysis_results"):
        return pd.DataFrame()
    conn = sqlite3.connect(db_path)
    try:
        return pd.read_sql_query(
            "SELECT * FROM analysis_results WHERE run_id=? ORDER BY rank ASC",
            conn,
            params=[run_id],
        )
    finally:
        conn.close()


def get_system_metrics(db_path: str) -> dict:
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        n_samples = 0
        if table_exists(db_path, "cfrna_samples"):
            cur.execute("SELECT COUNT(*) FROM cfrna_samples")
            n_samples = int(cur.fetchone()[0])
        n_legacy = n_v2 = 0
        if table_exists(db_path, "source_tracing_results"):
            cur.execute("SELECT COUNT(*) FROM source_tracing_results")
            n_legacy = int(cur.fetchone()[0])
        if table_exists(db_path, "analysis_runs"):
            cur.execute("SELECT COUNT(*) FROM analysis_runs")
            n_v2 = int(cur.fetchone()[0])
        return {
            "n_samples": n_samples,
            "n_analyses": n_legacy + n_v2,
            "n_analyses_legacy": n_legacy,
            "n_analyses_v2": n_v2,
        }
    finally:
        conn.close()
