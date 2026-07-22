from __future__ import annotations

import io
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from app.pages.tracing_page import (
    _locked_route_expression,
    _read_demo_expression,
    _run_locked_bo2023_route,
)
from app.pages.upload_page import _input_audit, _read_uploaded_file, _upload_mode_defaults
from core.network_tracing import load_network_model
from data_processor import DataProcessor


def test_input_audit_reports_raw_count_reads_and_model_overlap() -> None:
    model_genes = load_network_model()["genes"].astype(str)
    frame = pd.DataFrame(
        {
            "gene_symbol": [model_genes[0], model_genes[1], "NOT_IN_MODEL"],
            "tpm_value": [1.0, 0.0, 2.0],
            "read_count": [10, 20, 30],
            "expression_unit": ["logCPM_from_raw_counts"] * 3,
        }
    )

    audit = _input_audit(frame)

    assert audit["unit"] == "logCPM_from_raw_counts"
    assert audit["valid_genes"] == 3
    assert audit["nonzero_genes"] == 2
    assert audit["total_reads"] == 60
    assert audit["model_overlap"] == 2
    assert audit["model_genes"] == len(model_genes)


def test_input_audit_omits_total_reads_for_logcpm() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": ["GENE1", "GENE2"],
            "tpm_value": [3.0, 4.0],
            "expression_unit": ["logCPM", "logCPM"],
        }
    )

    audit = _input_audit(frame)

    assert audit["unit"] == "logCPM"
    assert audit["total_reads"] is None
    assert audit["valid_genes"] == 2
    assert audit["nonzero_genes"] == 2


def test_upload_defaults_do_not_fabricate_experimental_metadata() -> None:
    defaults = _upload_mode_defaults("rhesus")

    assert defaults["subject_id"] == ""
    assert defaults["age_years"] is None
    assert defaults["diagnosis"] == ""
    assert defaults["post_op_day"] is None
    assert defaults["sequencing_platform"] == ""
    assert defaults["total_reads"] is None
    assert defaults["mapping_rate"] is None


def test_upload_reader_accepts_uppercase_csv_extension() -> None:
    uploaded = io.BytesIO(b"gene_symbol,raw_counts\nA1BG,10\nA2M,20\n")
    uploaded.name = "SAMPLE.CSV"

    result = _read_uploaded_file(uploaded)

    assert result["gene_symbol"].tolist() == ["A1BG", "A2M"]


def test_standard_logcpm_column_is_consistent_across_upload_and_tracing_entries() -> None:
    payload = (
        "gene_symbol,logCPM\n"
        + "\n".join(f"GENE{i},{i / 10:.1f}" for i in range(1, 11))
        + "\n"
    ).encode("utf-8")
    upload_file = io.BytesIO(payload)
    upload_file.name = "sample.csv"
    tracing_file = io.BytesIO(payload)
    tracing_file.name = "sample.csv"

    uploaded = _read_uploaded_file(upload_file)
    processor = DataProcessor()
    is_valid, errors = processor.validate_expression_data(uploaded)
    processed = processor.preprocess_expression_data(uploaded, min_tpm=0.0)
    tracing_expression, tracing_source = _read_demo_expression(tracing_file)

    assert is_valid, errors
    assert processed["expression_unit"].unique().tolist() == ["logCPM"]
    assert tracing_source == "logcpm"
    assert processed.set_index("gene_symbol")["tpm_value"].to_dict() == (
        tracing_expression.set_index("gene_symbol")["log_tpm"].to_dict()
    )


def test_saved_logcpm_preserves_direct_three_tier_ranking(tmp_path) -> None:
    model_path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "models"
        / "bo2023_formal_region_logcpm_reference_matrix.npz"
    )
    with np.load(model_path, allow_pickle=False) as model:
        genes = model["genes"].astype(str)
        values = model["matrix"][:, 0].astype(float)

    payload = (
        "gene_symbol,logCPM\n"
        + "\n".join(f"{gene},{value:.9g}" for gene, value in zip(genes, values))
        + "\n"
    ).encode("utf-8")
    direct_file = io.BytesIO(payload)
    direct_file.name = "sample.csv"
    upload_file = io.BytesIO(payload)
    upload_file.name = "sample.csv"

    direct_expression, direct_source = _read_demo_expression(direct_file)
    uploaded = _read_uploaded_file(upload_file)
    processor_input = DataProcessor().preprocess_expression_data(uploaded, min_tpm=0.0)

    db_path = tmp_path / "samples.db"
    with sqlite3.connect(db_path) as conn:
        conn.executescript(
            """
            CREATE TABLE braintrace_samples (
                sample_id TEXT PRIMARY KEY,
                subject_id TEXT,
                species TEXT,
                qc_status TEXT,
                gene_id_type TEXT,
                brain_traceability TEXT
            );
            CREATE TABLE braintrace_expression (
                sample_id TEXT,
                gene_symbol TEXT NOT NULL,
                tpm_value REAL,
                read_count REAL,
                detected INTEGER,
                log_tpm REAL,
                zscore_tpm REAL,
                gene_id_type TEXT,
                expression_unit TEXT,
                FOREIGN KEY (sample_id) REFERENCES braintrace_samples(sample_id)
            );
            """
        )
    processor = DataProcessor(str(db_path))
    processor.save_sample_with_expression(
        {"sample_id": "LOGCPM1", "subject_id": "", "species": "Macaca mulatta"},
        processor_input,
    )
    stored_rows = processor.get_sample_expression("LOGCPM1")
    stored_expression, stored_source = _locked_route_expression(stored_rows)

    assert direct_source == "logcpm"
    assert stored_source == "stored_logCPM"
    assert len(direct_expression) == len(processor_input) == len(stored_rows) == len(genes)
    pd.testing.assert_frame_equal(direct_expression, stored_expression)

    direct_network, direct_regions = _run_locked_bo2023_route(direct_expression, atlas_id=None)
    stored_network, stored_regions = _run_locked_bo2023_route(stored_expression, atlas_id=None)

    assert stored_network["results"] == direct_network["results"]
    assert (
        stored_regions["meta"]["region_resolution_annotation"]["group_ranking"]
        == direct_regions["meta"]["region_resolution_annotation"]["group_ranking"]
    )
    assert stored_regions["results"] == direct_regions["results"]


def test_raw_counts_min_zero_keeps_zero_rows_and_recomputes_logcpm() -> None:
    raw = pd.DataFrame(
        {
            "gene_symbol": [f"GENE{i}" for i in range(10)],
            "raw_counts": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        }
    )

    processed = DataProcessor().preprocess_expression_data(raw, min_tpm=0.0)

    assert len(processed) == len(raw)
    assert processed["read_count"].sum() == raw["raw_counts"].sum()
    assert processed["expression_unit"].unique().tolist() == ["logCPM_from_raw_counts"]
