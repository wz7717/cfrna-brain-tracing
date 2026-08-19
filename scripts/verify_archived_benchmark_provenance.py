#!/usr/bin/env python
"""Verify that the immutable archived GSE189919 benchmark has not changed."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.provenance_hashes import sha256_utf8_lf_text

EXPECTED_ORIGIN_SHA256 = "9442D7901F47A9A0A826E4502A58B0C4D1ADE61FCA70912890D66FDECC8498A3"
EXPECTED_STAGED_SHA256 = "958D23E17D60FF3610780D41F4B9A83A55AEF10F851A7D61AE4A11E7826DDA5C"
EXPECTED_STAGED_PATH = "reproducibility/formal_real_input_performance_manifest.json"


def validate(manifest_path: Path, provenance_path: Path) -> dict[str, Any]:
    errors: list[str] = []
    observed_sha = sha256_utf8_lf_text(manifest_path) if manifest_path.is_file() else None
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        provenance = {}
        errors.append(f"cannot read archived benchmark provenance: {exc}")
    if observed_sha != EXPECTED_STAGED_SHA256:
        errors.append("archived benchmark manifest SHA-256 differs from the accepted immutable artifact")
    if provenance.get("staged_path") != EXPECTED_STAGED_PATH:
        errors.append("archived benchmark staged path changed")
    if provenance.get("origin_sha256") != EXPECTED_ORIGIN_SHA256:
        errors.append("archived benchmark origin SHA-256 changed")
    if provenance.get("staged_sha256") != EXPECTED_STAGED_SHA256:
        errors.append("archived benchmark staged SHA-256 changed")
    return {
        "schema": "braintrace.archived_benchmark_provenance_verification.v1",
        "status": "PASS" if not errors else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT",
        "archived_manifest_sha256": observed_sha,
        "accepted_staged_manifest_sha256": EXPECTED_STAGED_SHA256,
        "accepted_origin_manifest_sha256": EXPECTED_ORIGIN_SHA256,
        "provenance_origin_path": provenance.get("origin_path"),
        "provenance_staged_path": provenance.get("staged_path"),
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = validate(args.manifest, args.provenance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
