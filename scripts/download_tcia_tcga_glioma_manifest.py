#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri"
NBIA_BASE = "https://services.cancerimagingarchive.net/nbia-api/services/v1"


def nbia_get_json(endpoint: str, params: dict[str, str]) -> list[dict[str, Any]]:
    params = {**params, "format": "json"}
    query = urllib.parse.urlencode(params)
    url = f"{NBIA_BASE}/{endpoint}?{query}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8")
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        preview = text[:300].replace("\n", " ")
        raise RuntimeError(f"NBIA did not return JSON for {url}; first bytes: {preview}") from exc


def normalize_patient_id(value: Any) -> str:
    text = str(value or "").strip().upper()
    if text.startswith("TCGA-") and len(text) >= 12:
        return text[:12]
    return text


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download TCIA/NBIA metadata manifests for TCGA glioma MRI collections without downloading images."
    )
    parser.add_argument(
        "--collections",
        nargs="+",
        default=["TCGA-GBM", "TCGA-LGG"],
        help="TCIA collection names to query.",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    all_patients: list[dict[str, Any]] = []
    all_series: list[dict[str, Any]] = []
    errors: dict[str, str] = {}
    for collection in args.collections:
        patients_path = args.outdir / f"{collection}_patients.csv"
        series_path = args.outdir / f"{collection}_series.csv"
        if args.skip_existing and patients_path.exists() and series_path.exists():
            patients = pd.read_csv(patients_path).to_dict("records")
            series = pd.read_csv(series_path).to_dict("records")
        else:
            try:
                patients = nbia_get_json("getPatient", {"Collection": collection})
                series = nbia_get_json("getSeries", {"Collection": collection})
            except (urllib.error.URLError, RuntimeError) as exc:
                errors[collection] = str(exc)
                continue
            pd.DataFrame(patients).to_csv(patients_path, index=False, encoding="utf-8-sig")
            pd.DataFrame(series).to_csv(series_path, index=False, encoding="utf-8-sig")

        for row in patients:
            row = dict(row)
            row["collection"] = collection
            row["patient_barcode"] = normalize_patient_id(row.get("PatientID"))
            all_patients.append(row)
        for row in series:
            row = dict(row)
            row["collection"] = collection
            row["patient_barcode"] = normalize_patient_id(row.get("PatientID"))
            all_series.append(row)

    patient_df = pd.DataFrame(all_patients).drop_duplicates()
    series_df = pd.DataFrame(all_series).drop_duplicates()
    patient_df.to_csv(args.outdir / "tcia_tcga_glioma_patients_manifest.csv", index=False, encoding="utf-8-sig")
    series_df.to_csv(args.outdir / "tcia_tcga_glioma_series_manifest.csv", index=False, encoding="utf-8-sig")
    summary = {
        "collections": args.collections,
        "n_patients": int(len(patient_df)),
        "n_series": int(len(series_df)),
        "collection_errors": errors,
        "patient_manifest": str(args.outdir / "tcia_tcga_glioma_patients_manifest.csv"),
        "series_manifest": str(args.outdir / "tcia_tcga_glioma_series_manifest.csv"),
        "note": "This script downloads metadata only. Use NBIA Data Retriever/TCIA portal or BraTS-TCGA NIfTI derivatives for image volumes.",
    }
    (args.outdir / "tcia_tcga_glioma_manifest_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
