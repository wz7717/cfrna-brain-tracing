from __future__ import annotations

from pathlib import Path

from scripts.audit_scientific_freeze import check_hashes


def test_freeze_hash_check_rejects_changed_text_artifact(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.json"
    artifact.write_text('{"value": 1}\n', encoding="utf-8")
    manifest = {
        "exact_text_sha256": {"artifact.json": "0" * 64},
        "exact_binary_sha256": {},
        "semantic_provenance_baseline_sha256": {},
    }
    records, errors = check_hashes(tmp_path, manifest)
    assert records[0]["status"] == "FAIL"
    assert errors == ["frozen artifact hash changed: artifact.json"]
