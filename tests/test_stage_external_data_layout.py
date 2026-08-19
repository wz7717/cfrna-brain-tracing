from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.stage_external_data_layout import stage


def test_staging_creates_canonical_hardlink_without_legacy_fallback(tmp_path: Path) -> None:
    source = tmp_path / "source"
    legacy = source / "legacy" / "input.txt"
    legacy.parent.mkdir(parents=True)
    legacy.write_bytes(b"input")
    manifest = {
        "schema": "braintrace.external_input_manifest.v1",
        "canonical_external_root": "external_data",
        "profiles": {"full": {}},
        "inputs": [
            {
                "alias": "fixture",
                "storage": "external",
                "kind": "file",
                "canonical_path": "Fixture/input.txt",
                "legacy_path": "legacy/input.txt",
                "sha256": hashlib.sha256(b"input").hexdigest(),
                "required": True,
                "profiles": ["full"],
                "classification": "downloadable public",
                "role": "test",
                "accession_or_source": "test",
            }
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    destination = tmp_path / "external_data"
    payload = stage(source, destination, manifest_path)
    staged = destination / "Fixture" / "input.txt"
    assert payload["status"] == "PASS"
    assert staged.read_bytes() == b"input"
    assert staged.stat().st_ino == legacy.stat().st_ino
