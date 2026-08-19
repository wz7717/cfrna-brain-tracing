from __future__ import annotations

import json
from pathlib import Path

from scripts.verify_huang2025_run import EXPECTED, validate


def test_huang_verifier_rejects_changed_profile_count(tmp_path: Path) -> None:
    path = tmp_path / "summary.json"
    payload = {
        **EXPECTED,
        "patient_paired_analysis": "NOT_SUPPORTED",
        "synthetic_matched_admixture": "REMOVED_FROM_CANONICAL_ANALYSIS",
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    assert validate(path)["status"] == "PASS"
    payload["n_csf"] = 78
    path.write_text(json.dumps(payload), encoding="utf-8")
    assert validate(path)["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
