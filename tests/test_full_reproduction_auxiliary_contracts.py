from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

from scripts.prepare_brats_tcga_lgg_bundle import run as stage_brats
from scripts.verify_gse189919_benchmark_run import validate as validate_benchmark
from scripts.verify_tcga_brats_truth_basis_run import EXPECTED, validate as validate_tcga


def test_brats_staging_generates_portable_relative_patient_audit(tmp_path: Path) -> None:
    archive = tmp_path / "brats.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        for index in range(65):
            patient = f"TCGA-XX-{index:04d}"
            root = f"Pre-operative_TCGA_LGG_NIfTI_and_Segmentations/{patient}"
            handle.writestr(f"{root}/{patient}_t1.nii.gz", b"t1")
            suffix = "GlistrBoost_ManuallyCorrected" if index == 0 else "GlistrBoost"
            handle.writestr(f"{root}/{patient}_{suffix}.nii.gz", b"seg")
    payload = stage_brats(archive, tmp_path / "stage")
    assert payload["status"] == "PASS"
    audit = (tmp_path / "stage" / "brats_tcga_lgg_training_65_patient_audit.csv").read_text(encoding="utf-8")
    assert "D:\\" not in audit
    assert "extracted/" in audit


def test_gse_benchmark_verifier_accepts_contract_not_hardware_timing(tmp_path: Path) -> None:
    manifest = {
        "status": "formal_real_input_performance_gate",
        "route": "projected_vsd_network_top3_logcpm_resolution_local_exact",
        "input": {"workload_samples": 51, "total_samples_in_header": 51},
        "preregistration": {"formal_warm_repeats": 3, "formal_execution_authorized": True},
        "cold": {"timing": {key: 1.0 for key in ("wall_total_seconds", "wall_seconds_per_sample", "sample_time_p50_seconds", "sample_time_p95_seconds")}},
        "warm": {
            "repeats": [
                {"timing": {"samples": 51, "wall_total_seconds": 1.0}}
                for _ in range(3)
            ],
            "aggregate": {
                "samples": 153,
                **{key: 1.0 for key in ("wall_total_seconds", "wall_seconds_per_sample", "sample_time_p50_seconds", "sample_time_p95_seconds")},
            },
        },
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    timing_path = tmp_path / "timing.csv"
    with timing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["phase"])
        writer.writeheader()
        writer.writerows([{"phase": "cold"} for _ in range(51)] + [{"phase": "warm"} for _ in range(153)])
    assert validate_benchmark(manifest_path, timing_path)["status"] == "PASS"


def test_tcga_truth_verifier_rejects_scientific_drift(tmp_path: Path) -> None:
    path = tmp_path / "truth.csv"
    fields = ["patient_barcode", "edema_voxels", "edema_network_dominant"] + [
        f"{truth}_{level}_top3_strict" for level, truth in EXPECTED
    ]
    rows: list[dict[str, object]] = []
    for index in range(65):
        row: dict[str, object] = {
            "patient_barcode": f"TCGA-XX-{index:04d}",
            "edema_voxels": 1 if index < 63 else (1 if index == 63 else 0),
            "edema_network_dominant": "network" if index < 63 else ("out_of_scope" if index == 63 else "network"),
        }
        for (level, truth), (correct, _n) in EXPECTED.items():
            eligible_index = index if truth != "edema" else index
            row[f"{truth}_{level}_top3_strict"] = eligible_index < correct
        rows.append(row)
    rows[63]["patient_barcode"] = "TCGA-HT-7680"
    rows[64]["patient_barcode"] = "TCGA-HT-7686"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    assert validate_tcga(path)["status"] == "PASS"
    rows[0]["center_network_top3_strict"] = False
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    assert validate_tcga(path)["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
