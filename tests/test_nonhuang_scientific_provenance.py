import hashlib
from pathlib import Path

from scripts import generate_nonhuang_scientific_provenance_artifacts as generator
from scripts.generate_nonhuang_scientific_provenance_artifacts import bh_adjust
from scripts.verify_nonhuang_scientific_provenance import run_checks


ROOT = Path(__file__).resolve().parents[1]


def test_nonhuang_scientific_provenance_arithmetic() -> None:
    payload = run_checks(ROOT)
    assert payload["status"] == "PASS"
    assert payload["n_checks"] >= 20


def test_four_test_bh_adjustment() -> None:
    adjusted = bh_adjust([0.03125, 0.375, 0.59375, 0.32421875])
    assert adjusted == [0.125, 0.5, 0.59375, 0.5]


def test_repository_text_source_hash_is_newline_canonical(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(generator, "REPO_ROOT", tmp_path)
    source = tmp_path / "source.json"
    source.write_bytes(b'{"value": 1}\n')
    lf_record = generator.source_record(source)
    source.write_bytes(b'{"value": 1}\r\n')
    crlf_record = generator.source_record(source)

    expected = hashlib.sha256(b'{"value": 1}\n').hexdigest()
    assert lf_record["sha256"] == crlf_record["sha256"] == expected
    assert lf_record["hash_mode"] == "utf8_lf_canonical_text"


def test_external_dataset_provenance_regressions() -> None:
    provenance = (ROOT / "DATA_PROVENANCE.md").read_text(encoding="utf-8")

    assert "10.7937/K9/TCIA.2017.GJQ7R0EF" in provenance
    assert "10.5281/zenodo.3718921" not in provenance
    assert "65 public training-set MRI cases" in provenance
    assert "CC BY 3.0" in provenance
    assert "59,453 raw gene rows" in provenance
    assert "72,108 expression rows" not in provenance

    detailed = (ROOT / "reproducibility" / "DATA_PROVENANCE.md").read_text(encoding="utf-8")
    assert "2 RNA-seq donors with a selected set of matched anatomical structures" in detailed
    assert "both with whole-brain coverage" not in detailed
    assert "6 total; 2 with whole-brain coverage" not in detailed
    assert "65-subject processed training set public (CC BY 3.0); separate 43-subject test arm controlled" in detailed
    assert "HRA007247" in detailed
    assert "underlying sequencing data controlled" in detailed
    assert "All 10 supplementary CSV files" not in detailed
    assert "scripts/generate_nonhuang_scientific_provenance_artifacts.py" in detailed

    trace_text = "\n".join(
        (ROOT / path).read_text(encoding="utf-8")
        for path in (
            "reproducibility/generate_all_csvs.py",
            "reproducibility/v4_p0_5_ahba_trace.csv",
            "reproducibility/v4_p0_5_ahba_trace_manuscript_aligned.csv",
        )
    )
    assert "whole-brain" not in trace_text
    assert "4 of 6 AHBA donors excluded" not in trace_text
    assert "selected matched" in trace_text

    ledger = (ROOT / "NONHUANG_SCIENTIFIC_CONFLICT_LEDGER.csv").read_text(encoding="utf-8")
    assert "2 RNA-seq donors with whole-brain coverage" not in ledger
    assert "2 RNA-seq donors with a selected set of matched anatomical structures" in ledger
