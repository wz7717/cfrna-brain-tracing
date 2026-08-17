from pathlib import Path

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
