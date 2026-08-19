#!/usr/bin/env python
"""Verify frozen-scope Huang2025 rerun counts and provenance guardrails."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


EXPECTED = {"n_profiles": 159, "n_csf": 77, "n_plasma": 82, "n_traceable_outputs": 159}
EXPECTED_ROUTE = "projected_vsd_network_top3_logcpm_resolution_local_exact"


def validate(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    errors: list[str] = []
    for key, expected in EXPECTED.items():
        if payload.get(key) != expected:
            errors.append(f"{key} changed: {payload.get(key)!r}, expected {expected}")
    if payload.get("patient_paired_analysis") != "NOT_SUPPORTED":
        errors.append("patient pairing guardrail changed")
    if payload.get("synthetic_matched_admixture") != "REMOVED_FROM_CANONICAL_ANALYSIS":
        errors.append("synthetic matched-mixture guardrail changed")
    return {
        "schema": "braintrace.huang2025_run_verification.v1",
        "status": "PASS" if not errors else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT",
        "observed": {key: payload.get(key) for key in EXPECTED},
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = validate(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
