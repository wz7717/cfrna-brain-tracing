#!/usr/bin/env python
"""Stage the verified BraTS-TCGA-LGG source bundle for a full reproduction run.

The external source is mounted read-only.  This helper extracts it only into
the caller-owned audit directory, writes a portable relative-path patient
audit, and refuses archive members that could escape the requested directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import zipfile
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_extract(archive: zipfile.ZipFile, destination: Path) -> None:
    root = destination.resolve()
    for member in archive.infolist():
        candidate = (root / member.filename).resolve()
        if candidate != root and root not in candidate.parents:
            raise ValueError(f"unsafe archive member: {member.filename}")
    archive.extractall(destination)


def _single(paths: list[Path], description: str, patient: str) -> Path:
    if len(paths) != 1:
        raise ValueError(f"{patient}: expected one {description}, found {len(paths)}")
    return paths[0]


def build_patient_audit(extracted_root: Path, audit_path: Path) -> list[dict[str, str]]:
    patients = sorted(path for path in extracted_root.iterdir() if path.is_dir() and path.name.startswith("TCGA-"))
    if len(patients) != 65:
        raise ValueError(f"expected 65 patient directories, found {len(patients)}")
    rows: list[dict[str, str]] = []
    for patient_dir in patients:
        patient = patient_dir.name.upper()
        t1 = _single(sorted(patient_dir.glob("*_t1.nii.gz")), "T1 NIfTI", patient)
        manual = sorted(patient_dir.glob("*_GlistrBoost_ManuallyCorrected.nii.gz"))
        automatic = sorted(patient_dir.glob("*_GlistrBoost.nii.gz"))
        segmentation = _single(manual or automatic, "segmentation NIfTI", patient)
        rows.append(
            {
                "patient_barcode": patient,
                "t1_path": t1.relative_to(audit_path.parent).as_posix(),
                "preferred_segmentation": segmentation.relative_to(audit_path.parent).as_posix(),
                "has_manual_segmentation": "True" if manual else "False",
            }
        )
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_barcode", "t1_path", "preferred_segmentation", "has_manual_segmentation"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def run(archive_path: Path, outdir: Path) -> dict[str, Any]:
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite staging directory: {outdir}")
    if not zipfile.is_zipfile(archive_path):
        raise ValueError("BraTS source is not a valid ZIP archive")
    outdir.mkdir(parents=True)
    try:
        with zipfile.ZipFile(archive_path) as archive:
            if archive.testzip() is not None:
                raise ValueError("BraTS source ZIP failed integrity test")
            _safe_extract(archive, outdir / "extracted")
        extracted_root = outdir / "extracted" / "Pre-operative_TCGA_LGG_NIfTI_and_Segmentations"
        if not extracted_root.is_dir():
            raise ValueError("BraTS source archive has no expected NIfTI root")
        audit_path = outdir / "brats_tcga_lgg_training_65_patient_audit.csv"
        rows = build_patient_audit(extracted_root, audit_path)
        payload: dict[str, Any] = {
            "schema": "braintrace.brats_bundle_staging.v1",
            "status": "PASS",
            "archive": {
                "locator": "external_source::brats_training_bundle/Pre-operative_TCGA_LGG_NIfTI_and_Segmentations.zip",
                "sha256": sha256_file(archive_path),
                "bytes": archive_path.stat().st_size,
            },
            "extraction": {
                "destination": "generated::brats_staging/extracted",
                "patient_directories": len(rows),
                "manual_segmentation_cases": sum(row["has_manual_segmentation"] == "True" for row in rows),
                "automatic_segmentation_cases": sum(row["has_manual_segmentation"] == "False" for row in rows),
            },
            "patient_audit": {
                "locator": "generated::brats_staging/brats_tcga_lgg_training_65_patient_audit.csv",
                "sha256": sha256_file(audit_path),
            },
        }
        (outdir / "BRATS_BUNDLE_STAGING_AUDIT.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n"
        )
        return payload
    except Exception:
        shutil.rmtree(outdir, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Verified read-only BraTS nested ZIP.")
    parser.add_argument("--outdir", type=Path, required=True, help="New caller-owned staging directory.")
    args = parser.parse_args()
    payload = run(args.archive, args.outdir)
    print(json.dumps({"status": payload["status"], "patient_directories": payload["extraction"]["patient_directories"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
