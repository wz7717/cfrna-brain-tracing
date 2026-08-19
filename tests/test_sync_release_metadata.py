from __future__ import annotations

from scripts.sync_release_metadata import FULL_REPRO_DOI, SOFTWARE_DOI, validate_manifest


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
