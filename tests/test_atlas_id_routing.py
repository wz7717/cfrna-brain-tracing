from __future__ import annotations

import sqlite3

from core.reference_loader import load_marker_signature_genes
from database_init import CSFRNASourceDatabase


def test_marker_signature_genes_filter_by_requested_atlas(tmp_path):
    db_path = tmp_path / "markers.db"
    conn = sqlite3.connect(db_path)
    try:
        conn.execute(
            """
            CREATE TABLE region_gene_signature (
                region_id TEXT,
                gene_symbol TEXT,
                specificity_score REAL,
                is_marker INTEGER,
                atlas_id INTEGER
            )
            """
        )
        conn.executemany(
            "INSERT INTO region_gene_signature VALUES (?, ?, ?, ?, ?)",
            [
                ("R1", "ATLAS1_GENE", 10.0, 1, 1),
                ("R1", "ATLAS4_GENE", 10.0, 1, 4),
                ("R2", "ATLAS4_SECOND", 9.0, 1, 4),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    assert load_marker_signature_genes(str(db_path), atlas_id=1) == ["ATLAS1_GENE"]
    assert load_marker_signature_genes(str(db_path), atlas_id=4) == [
        "ATLAS4_GENE",
        "ATLAS4_SECOND",
    ]


def test_fresh_database_schema_includes_atlas_id_columns(tmp_path):
    db_path = tmp_path / "fresh.db"
    db = CSFRNASourceDatabase(str(db_path))
    conn = db.connect()
    try:
        db.create_database_schema()
        for table in [
            "macaque_brain_atlas",
            "reference_expression",
            "region_gene_signature",
        ]:
            columns = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
            assert "atlas_id" in columns
    finally:
        conn.close()
