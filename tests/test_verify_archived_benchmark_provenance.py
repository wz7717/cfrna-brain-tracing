from __future__ import annotations

import json
from pathlib import Path

from scripts.verify_archived_benchmark_provenance import EXPECTED_ORIGIN_SHA256, EXPECTED_STAGED_SHA256, validate


def test_archived_benchmark_provenance_rejects_changed_hash(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text("immutable benchmark\n", encoding="utf-8")
    provenance = tmp_path / "provenance.json"
    provenance.write_text(json.dumps({
        "origin_path": "external_source::manifest.json",
        "origin_sha256": EXPECTED_ORIGIN_SHA256,
        "staged_path": "reproducibility/formal_real_input_performance_manifest.json",
        "staged_sha256": EXPECTED_STAGED_SHA256,
    }), encoding="utf-8")
    payload = validate(manifest, provenance)
    assert payload["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
    assert any("manifest SHA-256" in error for error in payload["errors"])
