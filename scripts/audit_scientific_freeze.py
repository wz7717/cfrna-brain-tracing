#!/usr/bin/env python
"""Audit frozen scientific artifacts against the accepted d8296f7 baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.lomo_exact_f1 import FORMAL_N, FORMAL_TOP1, FORMAL_TOP3, compute_lomo_exact_metrics, load_formal_predictions
from core.lomo_f1 import CANONICAL_FORMAL_PATH as LOMO_NETWORK_PATH, compute_lomo_network_metrics, load_formal_predictions as load_lomo_network_predictions
from core.model_lock import verify_locked_model_bundle
from core.provenance_hashes import sha256_utf8_lf_text
from core.resolution_group_baselines import EXPECTED as GROUP_EXPECTED, CANONICAL_PATHS, load_formal_rows
from scripts.verify_archived_benchmark_provenance import validate as validate_archived_benchmark
from scripts.verify_huang2025_run import validate as validate_huang
from scripts.verify_nonhuang_scientific_provenance import run_checks as verify_current_provenance


DEFAULT_MANIFEST = ROOT / "reproducibility" / "scientific_freeze_manifest.json"


def _raw_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def check_hashes(root: Path, manifest: dict[str, Any]) -> tuple[list[dict[str, str]], list[str]]:
    """Compare immutable artifact content and record semantic-only provenance files."""

    records: list[dict[str, str]] = []
    errors: list[str] = []
    for section, hasher in (
        ("exact_text_sha256", sha256_utf8_lf_text),
        ("exact_binary_sha256", _raw_sha256),
    ):
        entries = manifest.get(section, {})
        if not isinstance(entries, dict):
            errors.append(f"freeze manifest {section} is malformed")
            continue
        for relative, expected in sorted(entries.items()):
            path = root / str(relative)
            if not path.is_file():
                errors.append(f"missing frozen artifact: {relative}")
                continue
            observed = hasher(path)
            status = "PASS" if observed == expected else "FAIL"
            records.append({"path": str(relative), "mode": section, "status": status})
            if status != "PASS":
                errors.append(f"frozen artifact hash changed: {relative}")
    semantic = manifest.get("semantic_provenance_baseline_sha256", {})
    if not isinstance(semantic, dict):
        errors.append("freeze manifest semantic_provenance_baseline_sha256 is malformed")
    else:
        for relative in sorted(semantic):
            path = root / str(relative)
            status = "PASS" if path.is_file() else "FAIL"
            records.append({"path": str(relative), "mode": "semantic_provenance", "status": status})
            if status != "PASS":
                errors.append(f"missing semantic-provenance artifact: {relative}")
    return records, errors


def _require(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def check_endpoint_contracts(root: Path) -> tuple[list[dict[str, str]], list[str]]:
    """Recompute the named frozen outcomes, independent of artifact hashing."""

    checks: list[dict[str, str]] = []
    errors: list[str] = []
    try:
        verify_locked_model_bundle()
        checks.append({"name": "locked_model_bundle", "status": "PASS"})
    except Exception as exc:  # pragma: no cover - exercised by an altered lock artifact
        errors.append(f"locked model bundle: {exc}")

    exact = compute_lomo_exact_metrics(load_formal_predictions(root / "reproducibility" / "p2_publication_completeness" / "formal_lomo_exact_region_detail.csv"))["summary"]
    _require(errors, (exact["n_samples"], exact["top1_correct"], exact["top3_correct"]) == (FORMAL_N, FORMAL_TOP1, FORMAL_TOP3), "LOMO Exact 177/812 or 346/812 changed")
    exact_expected = {
        "micro_f1": 0.21798029556650247,
        "macro_f1": 0.20341533397979938,
        "sd_class_f1": 0.2159329084651328,
        "median_class_f1": 0.15384615384615385,
        "iqr_class_f1": 0.2857142857142857,
        "weighted_f1": 0.21352686247014377,
        "conditional_macro_f1_nonzero": 0.2858810099175559,
    }
    for name, expected in exact_expected.items():
        _require(errors, math.isclose(float(exact[name]), expected, abs_tol=1e-15), f"LOMO Exact {name} changed")
    _require(errors, exact["n_classes"] == 104 and exact["n_zero_f1_classes"] == 30, "LOMO Exact class universe or zero-F1 count changed")
    checks.append({"name": "lomo_exact", "status": "PASS" if not errors else "FAIL"})

    network = compute_lomo_network_metrics(load_lomo_network_predictions(LOMO_NETWORK_PATH))["summary"]
    _require(errors, (network["top1_correct"], network["top3_correct"]) == (455, 750), "LOMO Network 455/819 or 750/819 changed")
    checks.append({"name": "lomo_network", "status": "PASS" if not errors else "FAIL"})

    for endpoint, expected in GROUP_EXPECTED.items():
        rows, _ = load_formal_rows(CANONICAL_PATHS[endpoint], endpoint)
        observed = (sum(int(float(row["group_hit1"])) for row in rows), sum(int(float(row["group_hit3"])) for row in rows), len(rows))
        target = (expected["top1"], expected["top3"], expected["n"])
        _require(errors, observed == target, f"{endpoint} Group frozen numerator/denominator changed")
    checks.append({"name": "resolution_groups", "status": "PASS" if not errors else "FAIL"})

    with (root / "reproducibility" / "v4_p0_9_triple_ci.csv").open(newline="", encoding="utf-8-sig") as handle:
        endpoint_rows = {row["metric"]: row for row in csv.DictReader(handle)}
    frozen_counts = {
        "LOSO Network Top1": (483, 819), "LOSO Network Top3": (753, 819),
        "LOMO Network Top1": (455, 819), "LOMO Network Top3": (750, 819),
        "LOSO Exact Top1": (182, 814), "LOSO Exact Top3": (368, 814),
        "LOMO Exact Top1": (177, 812), "LOMO Exact Top3": (346, 812),
        "LOSO ResGroup Top1": (368, 814), "LOSO ResGroup Top3": (590, 814),
        "LOMO ResGroup Top1": (344, 812), "LOMO ResGroup Top3": (569, 812),
    }
    for endpoint, (correct, denominator) in frozen_counts.items():
        row = endpoint_rows.get(endpoint, {})
        _require(errors, (int(row.get("correct", -1)), int(row.get("n", -1))) == (correct, denominator), f"{endpoint} frozen numerator/denominator changed")
    checks.append({"name": "loso_lomo_endpoint_counts", "status": "PASS" if not errors else "FAIL"})

    provenance = verify_current_provenance(root)
    _require(errors, provenance.get("status") == "PASS", "current scientific provenance contract failed")
    checks.append({"name": "current_scientific_provenance", "status": provenance.get("status", "FAIL")})

    archived = validate_archived_benchmark(
        root / "reproducibility" / "formal_real_input_performance_manifest.json",
        root / "reproducibility" / "formal_real_input_performance_provenance.json",
    )
    _require(errors, archived.get("status") == "PASS", "archived GSE189919 benchmark provenance changed")
    checks.append({"name": "archived_gse189919_benchmark", "status": str(archived.get("status", "FAIL"))})

    huang = validate_huang(root / "reproducibility" / "huang_2025" / "huang_2025_canonical_summary.json")
    _require(errors, huang.get("status") == "PASS", "Huang2025 profile-count/provenance contract changed")
    checks.append({"name": "huang2025", "status": str(huang.get("status", "FAIL"))})
    return checks, errors


def audit(root: Path, manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    hash_records, errors = check_hashes(root, manifest)
    contract_checks, contract_errors = check_endpoint_contracts(root)
    errors.extend(contract_errors)
    return {
        "schema": "braintrace.scientific_freeze_audit.v1",
        "accepted_scientific_baseline": manifest.get("accepted_scientific_baseline"),
        "status": "PASS" if not errors else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT",
        "scientific_drift_count": len(errors),
        "immutable_artifact_checks": hash_records,
        "semantic_contract_checks": contract_checks,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit(args.repo_root.resolve(), args.manifest.resolve())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"], "scientific_drift_count": payload["scientific_drift_count"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
