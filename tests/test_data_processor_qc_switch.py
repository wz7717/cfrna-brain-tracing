from __future__ import annotations

import sqlite3

import pandas as pd

from data_processor import DataProcessor


def _create_database(path: str) -> None:
    with sqlite3.connect(path) as conn:
        conn.executescript(
            """
            CREATE TABLE cfrna_samples (
                sample_id TEXT PRIMARY KEY,
                subject_id TEXT,
                species TEXT,
                diagnosis TEXT,
                collection_date TEXT,
                qc_status TEXT,
                gene_id_type TEXT,
                brain_traceability TEXT
            );
            CREATE TABLE cfrna_expression (
                sample_id TEXT,
                gene_symbol TEXT NOT NULL,
                tpm_value REAL,
                read_count REAL,
                detected INTEGER,
                expression_unit TEXT,
                FOREIGN KEY (sample_id) REFERENCES cfrna_samples(sample_id)
            );
            CREATE TABLE sample_qc (sample_id TEXT);
            """
        )
        conn.execute(
            "INSERT INTO cfrna_samples VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            ("S1", "P1", "Macaca mulatta", "", "", None, None, None),
        )


def test_get_all_samples_is_independent_of_upload_qc_switch(tmp_path) -> None:
    db_path = tmp_path / "samples.db"
    _create_database(str(db_path))
    processor = DataProcessor(str(db_path))

    samples = processor.get_all_samples()

    assert samples["sample_id"].tolist() == ["S1"]


def test_data_processor_connections_enable_foreign_keys(tmp_path) -> None:
    db_path = tmp_path / "samples.db"
    _create_database(str(db_path))
    processor = DataProcessor(str(db_path))

    with processor._get_conn() as conn:
        assert conn.execute("PRAGMA foreign_keys").fetchone()[0] == 1


def test_save_expression_can_skip_legacy_qc_without_affecting_storage(tmp_path) -> None:
    db_path = tmp_path / "samples.db"
    _create_database(str(db_path))
    processor = DataProcessor(str(db_path))
    expression = pd.DataFrame(
        {
            "gene_symbol": ["A1BG", "A2M"],
            "tpm_value": [1.0, 2.0],
            "read_count": [10.0, 20.0],
            "detected": [1, 1],
            "expression_unit": ["logCPM_from_raw_counts"] * 2,
        }
    )

    processor.save_expression_data("S1", expression, run_qc=False)

    assert len(processor.get_sample_expression("S1")) == 2
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT count(*) FROM sample_qc").fetchone()[0] == 0
        assert conn.execute("SELECT qc_status FROM cfrna_samples WHERE sample_id='S1'").fetchone()[0] is None


def test_atomic_sample_save_commits_metadata_and_expression_together(tmp_path) -> None:
    db_path = tmp_path / "samples.db"
    _create_database(str(db_path))
    processor = DataProcessor(str(db_path))
    expression = pd.DataFrame(
        {
            "gene_symbol": ["A1BG", "A2M"],
            "tpm_value": [1.0, 2.0],
            "detected": [1, 1],
        }
    )

    processor.save_sample_with_expression(
        {"sample_id": "S2", "subject_id": "P2", "species": "Macaca mulatta", "qc_status": None},
        expression,
    )

    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT subject_id FROM cfrna_samples WHERE sample_id='S2'").fetchone()[0] == "P2"
        assert conn.execute("SELECT count(*) FROM cfrna_expression WHERE sample_id='S2'").fetchone()[0] == 2


def test_atomic_sample_save_rolls_back_metadata_when_expression_insert_fails(tmp_path) -> None:
    db_path = tmp_path / "samples.db"
    _create_database(str(db_path))
    processor = DataProcessor(str(db_path))
    invalid_expression = pd.DataFrame(
        {
            "gene_symbol": [None],
            "tpm_value": [1.0],
            "detected": [1],
        }
    )

    try:
        processor.save_sample_with_expression(
            {"sample_id": "S2", "subject_id": "P2", "species": "Macaca mulatta"},
            invalid_expression,
        )
    except sqlite3.IntegrityError:
        pass
    else:
        raise AssertionError("Expected the expression insert to fail")

    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT count(*) FROM cfrna_samples WHERE sample_id='S2'").fetchone()[0] == 0
        assert conn.execute("SELECT count(*) FROM cfrna_expression WHERE sample_id='S2'").fetchone()[0] == 0


def test_metadata_upsert_does_not_delete_existing_expression(tmp_path) -> None:
    db_path = tmp_path / "samples.db"
    _create_database(str(db_path))
    processor = DataProcessor(str(db_path))
    expression = pd.DataFrame(
        {"gene_symbol": ["A1BG"], "tpm_value": [1.0], "detected": [1]}
    )
    processor.save_expression_data("S1", expression, run_qc=False)

    processor.save_sample_metadata(
        {"sample_id": "S1", "subject_id": "UPDATED", "species": "Macaca mulatta"}
    )

    assert len(processor.get_sample_expression("S1")) == 1
    assert processor.get_all_samples().iloc[0]["subject_id"] == "UPDATED"
