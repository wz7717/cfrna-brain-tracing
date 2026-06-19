
from __future__ import annotations

import re
import sqlite3
from pathlib import Path
from typing import Dict

import pandas as pd


def sanitize_table_name(filename: str) -> str:
    stem = Path(filename).stem
    stem = re.sub(r'^\d+_', '', stem)
    stem = stem.lower()
    stem = re.sub(r'[^a-z0-9]+', '_', stem).strip('_')
    return f'bo2023_{stem}'


def sanitize_column_name(name: str, used: set[str] | None = None) -> str:
    used = used if used is not None else set()
    raw = str(name).strip().lower()
    raw = raw.replace('%', 'pct').replace('#', 'num')
    raw = re.sub(r'[^a-z0-9]+', '_', raw).strip('_') or 'col'
    cand = raw
    i = 2
    while cand in used:
        cand = f'{raw}_{i}'
        i += 1
    used.add(cand)
    return cand


def _normalize_dataframe(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    used: set[str] = set()
    mapping = {c: sanitize_column_name(c, used) for c in df.columns}
    out = df.rename(columns=mapping).copy()
    for col in out.columns:
        if out[col].dtype == 'object':
            out[col] = out[col].where(~out[col].isna(), None)
    return out, mapping


def ensure_metadata_tables(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS bo2023_buildkit_catalog (
            table_name TEXT PRIMARY KEY,
            source_file TEXT NOT NULL,
            source_rows INTEGER,
            source_cols INTEGER,
            imported_at TEXT DEFAULT (datetime('now')),
            notes TEXT
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS bo2023_buildkit_column_map (
            table_name TEXT NOT NULL,
            source_column TEXT NOT NULL,
            sqlite_column TEXT NOT NULL,
            PRIMARY KEY (table_name, source_column)
        )
        """
    )
    conn.commit()


def import_buildkit_dir(
    db_path: str | Path,
    buildkit_dir: str | Path,
    include_readme: bool = False,
) -> Dict[str, int]:
    db_path = str(db_path)
    buildkit_dir = Path(buildkit_dir)
    csv_files = sorted(buildkit_dir.glob('*.csv'))
    conn = sqlite3.connect(db_path)
    try:
        ensure_metadata_tables(conn)
        imported = 0
        for csv_path in csv_files:
            table = sanitize_table_name(csv_path.name)
            df = pd.read_csv(csv_path)
            clean, mapping = _normalize_dataframe(df)
            clean.to_sql(table, conn, if_exists='replace', index=False)
            conn.execute('DELETE FROM bo2023_buildkit_catalog WHERE table_name=?', (table,))
            conn.execute(
                'INSERT INTO bo2023_buildkit_catalog(table_name, source_file, source_rows, source_cols, notes) VALUES (?, ?, ?, ?, ?)',
                (table, csv_path.name, int(df.shape[0]), int(df.shape[1]), 'Imported from Bo2023 supplementary build kit'),
            )
            conn.execute('DELETE FROM bo2023_buildkit_column_map WHERE table_name=?', (table,))
            conn.executemany(
                'INSERT INTO bo2023_buildkit_column_map(table_name, source_column, sqlite_column) VALUES (?, ?, ?)',
                [(table, str(src), str(dst)) for src, dst in mapping.items()],
            )
            imported += 1

        _create_helper_views(conn)
        _create_helper_indexes(conn)
        conn.commit()
        return {'tables_imported': imported, 'views_created': 4}
    finally:
        conn.close()


def _safe_drop_view(conn: sqlite3.Connection, name: str) -> None:
    conn.execute(f'DROP VIEW IF EXISTS {name}')


def _create_helper_views(conn: sqlite3.Connection) -> None:
    _safe_drop_view(conn, 'bo2023_regions')
    conn.execute(
        """
        CREATE VIEW bo2023_regions AS
        SELECT
            region,
            roi173,
            neocortexregion,
            lobe,
            full_name,
            dictionary_lobe,
            regional_map
        FROM bo2023_region_annotation_joined
        """
    )
    _safe_drop_view(conn, 'bo2023_region_qc_overview')
    conn.execute(
        """
        CREATE VIEW bo2023_region_qc_overview AS
        SELECT
            region,
            neocortexregion,
            lobe,
            rin,
            total_reads,
            uniquely_mapped,
            uniquely_mapped_percent,
            pct_pf_reads_aligned,
            pct_mrna_bases,
            pct_usable_bases,
            pct_correct_strand_reads,
            median_5prime_to_3prime_bias,
            at_dropout,
            gc_dropout
        FROM bo2023_region_qc_summary
        """
    )
    _safe_drop_view(conn, 'bo2023_region_sample_coverage')
    conn.execute(
        """
        CREATE VIEW bo2023_region_sample_coverage AS
        SELECT
            region,
            roi173,
            neocortexregion,
            lobe,
            total_sample_count
        FROM bo2023_region_sample_presence_matrix
        """
    )
    _safe_drop_view(conn, 'bo2023_ct_gene_catalog')
    conn.execute(
        """
        CREATE VIEW bo2023_ct_gene_catalog AS
        SELECT * FROM bo2023_ct_related_genes_1005
        """
    )


def _create_helper_indexes(conn: sqlite3.Connection) -> None:
    statements = [
        'CREATE INDEX IF NOT EXISTS idx_bo2023_regions_region ON bo2023_region_annotation_joined(region)',
        'CREATE INDEX IF NOT EXISTS idx_bo2023_region_qc_region ON bo2023_region_qc_summary(region)',
        'CREATE INDEX IF NOT EXISTS idx_bo2023_region_presence_region ON bo2023_region_sample_presence_matrix(region)',
        'CREATE INDEX IF NOT EXISTS idx_bo2023_deg_region ON bo2023_region_deg_summary(region)',
        'CREATE INDEX IF NOT EXISTS idx_bo2023_manifest_output_file ON bo2023_manifest(output_file)',
    ]
    for sql in statements:
        try:
            conn.execute(sql)
        except sqlite3.OperationalError:
            pass


def list_imported_tables(db_path: str | Path) -> pd.DataFrame:
    conn = sqlite3.connect(str(db_path))
    try:
        return pd.read_sql_query('SELECT * FROM bo2023_buildkit_catalog ORDER BY table_name', conn)
    finally:
        conn.close()
