#!/usr/bin/env python
"""Compare a regenerated NPZ artifact to its frozen scientific counterpart."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.external_inputs import sha256_file


VOLATILE_METADATA_KEYS = frozenset({"created_at_utc", "generated_at_utc"})


def _metadata_equal(left: np.ndarray, right: np.ndarray) -> tuple[bool, list[str]]:
    """Compare JSON metadata while excluding only recorded generation timestamps."""

    try:
        expected = json.loads(str(left.item()))
        observed = json.loads(str(right.item()))
    except (AttributeError, TypeError, ValueError, json.JSONDecodeError):
        return np.array_equal(left, right), []
    if not isinstance(expected, dict) or not isinstance(observed, dict):
        return expected == observed, []
    ignored = sorted((set(expected) | set(observed)) & VOLATILE_METADATA_KEYS)
    for key in ignored:
        expected.pop(key, None)
        observed.pop(key, None)
    return expected == observed, ignored


def compare(canonical: Path, regenerated: Path, canonical_label: str, regenerated_label: str) -> dict[str, Any]:
    mismatches: list[dict[str, Any]] = []
    ignored_metadata_fields: list[str] = []
    with np.load(canonical, allow_pickle=False) as expected, np.load(regenerated, allow_pickle=False) as observed:
        expected_keys = sorted(expected.files)
        observed_keys = sorted(observed.files)
        if expected_keys != observed_keys:
            mismatches.append({"field": "keys", "expected": expected_keys, "observed": observed_keys})
        for key in sorted(set(expected_keys) & set(observed_keys)):
            left = expected[key]
            right = observed[key]
            if left.dtype != right.dtype or left.shape != right.shape:
                mismatches.append(
                    {
                        "field": key,
                        "reason": "dtype_or_shape",
                        "expected_dtype": str(left.dtype),
                        "observed_dtype": str(right.dtype),
                        "expected_shape": list(left.shape),
                        "observed_shape": list(right.shape),
                    }
                )
            elif key == "metadata":
                equal, ignored = _metadata_equal(left, right)
                ignored_metadata_fields.extend(ignored)
                if not equal:
                    mismatches.append({"field": key, "reason": "metadata_values"})
            elif not np.array_equal(left, right):
                mismatches.append({"field": key, "reason": "array_values"})
    return {
        "schema": "braintrace.regenerated_npz_comparison.v1",
        "status": "PASS" if not mismatches else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT",
        "canonical_artifact": canonical_label,
        "regenerated_artifact": regenerated_label,
        "canonical_sha256": sha256_file(canonical),
        "regenerated_sha256": sha256_file(regenerated),
        "comparison": "exact NPZ key/dtype/shape/value comparison; only created_at_utc/generated_at_utc metadata timestamps are ignored; serialization SHA may differ",
        "ignored_volatile_metadata_fields": sorted(set(ignored_metadata_fields)),
        "n_mismatch": len(mismatches),
        "mismatches": mismatches,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--regenerated", type=Path, required=True)
    parser.add_argument("--canonical-label", required=True)
    parser.add_argument("--regenerated-label", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = compare(args.canonical, args.regenerated, args.canonical_label, args.regenerated_label)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"], "n_mismatch": payload["n_mismatch"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
