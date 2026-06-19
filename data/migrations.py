from __future__ import annotations

import sqlite3

# 允许执行 PRAGMA / ALTER TABLE 操作的表白名单
_ALLOWED_TABLES = frozenset({
    'cfrna_expression',
    'cfrna_samples',
    'macaque_brain_atlas',
    'reference_expression',
    'region_gene_signature',
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
    """校验表名在白名单内，防止 PRAGMA / ALTER TABLE 格式化注入。"""
    if name not in _ALLOWED_TABLES:
        raise ValueError(f"Unknown or disallowed table name: {name!r}")
    return name


def _table_exists(conn: sqlite3.Connection, table: str) -> bool:
    cur = conn.cursor()
    cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    )
    return cur.fetchone() is not None


def _col_exists(conn: sqlite3.Connection, table: str, col: str) -> bool:
    _validate_table_name(table)
    cur = conn.cursor()
    cur.execute(f"PRAGMA table_info({table})")
    return col in {row[1] for row in cur.fetchall()}


def _ensure_column(
    conn: sqlite3.Connection, table: str, col: str, decl: str
) -> None:
    if _table_exists(conn, table) and not _col_exists(conn, table, col):
        _validate_table_name(table)
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} {decl}")


def run_migrations(db_path: str) -> None:
    conn = sqlite3.connect(db_path)
    try:
        _ensure_column(conn, 'cfrna_expression', 'log_tpm', 'REAL')
        _ensure_column(conn, 'cfrna_expression', 'gene_id_type', 'TEXT')
        _ensure_column(conn, 'cfrna_expression', 'expression_unit', 'TEXT')

        for table in ('macaque_brain_atlas', 'reference_expression', 'region_gene_signature'):
            _ensure_column(conn, table, 'atlas_id', 'INTEGER DEFAULT 1')

        for table in ('reference_expression', 'region_gene_signature', 'macaque_brain_atlas'):
            if _table_exists(conn, table):
                conn.execute(
                    f"UPDATE {table} SET atlas_id = COALESCE(atlas_id, 1)"
                )

        for col, decl in [
            ('plasma_volume_ml', 'REAL'),
            ('sample_type', 'TEXT'),
            ('gene_id_type', 'TEXT'),
            ('brain_traceability', 'TEXT'),
            ('post_op_day', 'REAL'),
            ('surgery_region', 'TEXT'),
            ('surgery_side', 'TEXT'),
        ]:
            _ensure_column(conn, 'cfrna_samples', col, decl)

        if _table_exists(conn, 'reference_expression'):
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_ref_atlas_region_gene "
                "ON reference_expression(atlas_id, region_id, gene_symbol)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_ref_celltype_marker "
                "ON reference_expression(cell_type_marker, atlas_id)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_ref_atlas_celltype "
                "ON reference_expression(atlas_id, cell_type_marker)"
            )
        if _table_exists(conn, 'region_gene_signature'):
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_sig_atlas_region_gene "
                "ON region_gene_signature(atlas_id, region_id, gene_symbol)"
            )
        if _table_exists(conn, 'signature_genes'):
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_signature_genes_sigset_region_gene "
                "ON signature_genes(sigset_id, region_id, gene_symbol)"
            )
        if _table_exists(conn, 'cfrna_expression'):
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_cfrna_sample_gene "
                "ON cfrna_expression(sample_id, gene_symbol)"
            )
        conn.commit()
    finally:
        conn.close()
