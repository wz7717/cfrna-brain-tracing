from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DB = ROOT / "data" / "processed_reference" / "v1" / "cfrna_brain_tracing_reference_v1.sqlite"


def connect(db_path: str | Path | None = None) -> sqlite3.Connection:
    return sqlite3.connect(str(db_path or DEFAULT_DB))


def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    row = conn.execute("select name from sqlite_master where type='table' and name=?", (table,)).fetchone()
    return row is not None


def read_table(table: str, db_path: str | Path | None = None) -> pd.DataFrame:
    with connect(db_path) as conn:
        if not table_exists(conn, table):
            return pd.DataFrame()
        return pd.read_sql_query(f"select * from {table}", conn)


def read_marker_candidates(db_path: str | Path | None = None) -> pd.DataFrame:
    markers = read_table("marker_candidates", db_path)
    if markers.empty:
        return markers
    markers["marker_type_normalized"] = markers["marker_type"].replace(
        {
            "brain_region": "brain_region_marker",
            "brain_celltype": "brain_celltype_marker",
            "tbi_injury": "tbi_injury_response",
        }
    )
    return markers
