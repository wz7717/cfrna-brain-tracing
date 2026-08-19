#!/usr/bin/env python
"""Regenerate deterministic current evidence and compare it to frozen artifacts.

This is intentionally a comparison gate, not a scientific remediation tool.
It writes only to the caller-provided audit location and never overwrites
canonical manuscript artifacts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.lomo_exact_f1 import (  # noqa: E402
    CANONICAL_FORMAL_PATH,
    cross3_summary_row,
    compute_lomo_exact_metrics,
    load_formal_predictions,
)
from core.resolution_group_baselines import (  # noqa: E402
    CANONICAL_PATHS,
    compute_baseline_record,
    load_formal_rows,
)


NUMERIC_TEXT_RE = re.compile(r"^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
NUMERIC_TOLERANCE = 5e-10  # Canonical class-level CSV serializes floats to 10 decimals.


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _compare(expected: Any, observed: Any, path: str, mismatches: list[dict[str, Any]]) -> None:
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(expected) != set(observed):
            mismatches.append({"path": path, "reason": "key_set", "expected": sorted(expected), "observed": sorted(observed)})
            return
        for key in expected:
            _compare(expected[key], observed[key], f"{path}.{key}", mismatches)
        return
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            mismatches.append({"path": path, "reason": "length", "expected": len(expected), "observed": len(observed)})
            return
        for index, (left, right) in enumerate(zip(expected, observed)):
            _compare(left, right, f"{path}[{index}]", mismatches)
        return
    if isinstance(expected, float) or isinstance(observed, float):
        try:
            if math.isclose(float(expected), float(observed), rel_tol=0, abs_tol=NUMERIC_TOLERANCE):
                return
        except (TypeError, ValueError):
            pass
    if isinstance(expected, str) and isinstance(observed, str) and NUMERIC_TEXT_RE.fullmatch(expected) and NUMERIC_TEXT_RE.fullmatch(observed):
        if math.isclose(float(expected), float(observed), rel_tol=0, abs_tol=NUMERIC_TOLERANCE):
            return
    if expected != observed:
        mismatches.append({"path": path, "reason": "value", "expected": expected, "observed": observed})


def _lomo_class_rows(metrics: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {
            "class": str(row["class"]),
            "support": str(row["support"]),
            "tp": str(row["tp"]),
            "fp": str(row["fp"]),
            "fn": str(row["fn"]),
            "predicted_positive": str(row["predicted_positive"]),
            "precision": f"{float(row['precision']):.10f}",
            "recall": f"{float(row['recall']):.10f}",
            "f1": f"{float(row['f1']):.10f}",
        }
        for row in metrics["classes"]
    ]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path) -> dict[str, Any]:
    output = output.resolve()
    generated = output.parent / ".regeneration_tmp"
    generated.mkdir(parents=True, exist_ok=True)
    mismatches: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []

    metrics = compute_lomo_exact_metrics(load_formal_predictions(CANONICAL_FORMAL_PATH))
    regenerated_lomo = _lomo_class_rows(metrics)
    generated_lomo_path = generated / "formal_lomo_exact_region_f1.csv"
    _write_csv(generated_lomo_path, regenerated_lomo)
    archived_lomo_path = ROOT / "reproducibility" / "formal_lomo_exact_region_f1.csv"
    archived_lomo = _read_csv(archived_lomo_path)
    before = len(mismatches)
    _compare(regenerated_lomo, archived_lomo, "formal_lomo_exact_region_f1", mismatches)
    comparisons.append(
        {
            "canonical_artifact": "reproducibility/formal_lomo_exact_region_f1.csv",
            "regenerated_artifact": "regenerated::formal_lomo_exact_region_f1.csv",
            "comparison_mode": "exact_schema_and_canonical_text_fields",
            "canonical_sha256": sha256_file(archived_lomo_path),
            "regenerated_sha256": sha256_file(generated_lomo_path),
            "n_mismatch": len(mismatches) - before,
        }
    )

    cross3_path = ROOT / "reproducibility" / "p1_cross1_5" / "cross3_f1_distribution_summary.csv"
    archived_cross3 = next(row for row in _read_csv(cross3_path) if row["endpoint"] == "LOMO_Exact")
    regenerated_cross3 = {key: str(value) for key, value in cross3_summary_row(metrics).items()}
    before = len(mismatches)
    _compare(regenerated_cross3, archived_cross3, "lomo_exact_cross3_summary", mismatches)
    comparisons.append(
        {
            "canonical_artifact": "reproducibility/p1_cross1_5/cross3_f1_distribution_summary.csv#LOMO_Exact",
            "regenerated_artifact": "regenerated::lomo_exact_cross3_summary",
            "comparison_mode": "exact_schema_with_numeric_tolerance_5e-10",
            "n_mismatch": len(mismatches) - before,
        }
    )

    baseline_path = ROOT / "reproducibility" / "formal_resolution_group_random_baselines.json"
    archived_baselines = json.loads(baseline_path.read_text(encoding="utf-8"))
    regenerated_baselines = []
    for endpoint in ("LOSO", "LOMO"):
        rows, _ = load_formal_rows(CANONICAL_PATHS[endpoint], endpoint)
        regenerated_baselines.append(compute_baseline_record(rows, endpoint))
    generated_baseline_path = generated / "formal_resolution_group_random_baselines.json"
    generated_baseline_path.write_text(json.dumps(regenerated_baselines, indent=2) + "\n", encoding="utf-8")
    archived_by_endpoint = {row["endpoint"]: row for row in archived_baselines["records"]}
    expected_baselines = [
        {key: archived_by_endpoint[row["endpoint"]][key] for key in row}
        for row in regenerated_baselines
    ]
    before = len(mismatches)
    _compare(regenerated_baselines, expected_baselines, "resolution_group_random_baselines", mismatches)
    comparisons.append(
        {
            "canonical_artifact": "reproducibility/formal_resolution_group_random_baselines.json#records",
            "regenerated_artifact": "regenerated::formal_resolution_group_random_baselines.json",
            "comparison_mode": "exact_schema_with_numeric_tolerance_5e-10",
            "canonical_sha256": sha256_file(baseline_path),
            "regenerated_sha256": sha256_file(generated_baseline_path),
            "n_mismatch": len(mismatches) - before,
        }
    )

    payload = {
        "schema": "braintrace.regeneration_comparison.v1",
        "status": "PASS" if not mismatches else "FAIL",
        "n_mismatch": len(mismatches),
        "comparisons": comparisons,
        "mismatches": mismatches,
        "policy": "Canonical staged inputs are regenerated into a temporary audit directory; no canonical scientific artifact is overwritten.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = run(args.output)
    print(json.dumps({"status": payload["status"], "n_mismatch": payload["n_mismatch"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
