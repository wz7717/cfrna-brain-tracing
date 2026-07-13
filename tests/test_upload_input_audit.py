from __future__ import annotations

import io

import pandas as pd

from app.pages.upload_page import _input_audit, _read_uploaded_file, _upload_mode_defaults
from core.network_tracing import load_network_model


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
