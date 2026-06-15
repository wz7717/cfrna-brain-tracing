from __future__ import annotations

import pandas as pd

from src.reference_tracing.scoring import run_reference_tracing


def test_reference_tracing_toy_matrix(tmp_path):
    genes = [
        "RBFOX3",
        "SNAP25",
        "SYT1",
        "GFAP",
        "MBP",
        "HBA1",
        "HBA2",
        "HBB",
        "PTPRC",
        "LST1",
        "TYROBP",
        "NEUROD6",
        "OPALIN",
        "AQP4",
        "CLDN5",
    ]
    df = pd.DataFrame(
        {
            "gene_symbol": genes,
            "brain_high": [500, 450, 400, 120, 80, 2, 2, 2, 3, 2, 2, 300, 200, 100, 90],
            "rbc_high": [10, 8, 6, 4, 4, 900, 850, 950, 20, 10, 10, 4, 4, 4, 4],
            "immune_high": [10, 8, 6, 4, 4, 20, 10, 10, 900, 850, 800, 4, 4, 4, 4],
            "mixed": [100, 80, 70, 50, 30, 200, 180, 220, 180, 160, 150, 80, 70, 60, 50],
            "low": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        }
    )
    result = run_reference_tracing(df, "gene_symbol", "raw counts", tmp_path, make_plots=False)
    summary = result["sample_overall_tracing_summary"]
    assert not summary.empty
    assert {"brain_signal_score", "rbc_risk", "immune_risk", "overall_interpretation"}.issubset(summary.columns)
    assert (tmp_path / "sample_overall_tracing_summary.tsv").exists()
    assert (tmp_path / "cfrna_tracing_report.md").exists()
    assert "rbc_high" in set(summary["sample"])
