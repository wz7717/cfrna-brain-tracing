from __future__ import annotations

import json

from scripts.sync_release_metadata import (
    FULL_REPRO_DOI,
    PROVENANCE_BODY_MARKERS,
    SOFTWARE_DOI,
    _check_release_gate,
    validate_manifest,
)


def test_release_manifest_rejects_substitute_reserved_doi() -> None:
    manifest = {
        "version": "0.1.17", "tag": "v0.1.17", "scientific_frozen": True,
        "software_version_doi": SOFTWARE_DOI, "full_repro_version_doi": "10.5281/zenodo.1",
        "creators": [{"name": "Author"}], "title": "BrainTrace", "description": "Stable description",
    }
    errors = validate_manifest(manifest)
    assert "full-repro DOI is not the reserved v0.1.17 draft DOI" in errors
    manifest["full_repro_version_doi"] = FULL_REPRO_DOI
    assert validate_manifest(manifest) == []


def test_pre_finalization_gate_defers_checksum_until_exact_archive_exists(tmp_path) -> None:
    gate = tmp_path / "release_gate.json"
    gate.write_text(
        json.dumps(
            {
                "status": "PASS",
                "scientific_drift": 0,
                "app_docker": "PASS",
                "repro_docker": "PASS",
                "github_actions": "GREEN",
                "checksum": "DEFERRED_UNTIL_FINAL_ARCHIVE",
            }
        ),
        encoding="utf-8",
    )
    assert _check_release_gate(gate) == []


def test_finalization_preserves_each_provenance_document_body_at_its_real_heading() -> None:
    assert PROVENANCE_BODY_MARKERS == {
        "DATA_PROVENANCE.md": "## Source Datasets",
        "reproducibility/DATA_PROVENANCE.md": "## 1. Primary Atlas Data",
    }


def test_finalization_can_reapply_generator_owned_release_state_rewrites() -> None:
    from scripts.sync_release_metadata import _replace_release_state_text

    assert _replace_release_state_text("before", "before", "after", "test") == "after"
    assert _replace_release_state_text("after", "before", "after", "test") == "after"
