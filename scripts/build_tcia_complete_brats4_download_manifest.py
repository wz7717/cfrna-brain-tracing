#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATCH = (
    ROOT
    / "results"
    / "tcga_rnaseq_tcia_mri_collection_match_20260605"
    / "tcga_rnaseq_patient_to_tcia_mri_match.csv"
)
DEFAULT_SERIES = (
    ROOT
    / "results"
    / "tcga_rnaseq_tcia_mri_collection_match_20260605"
    / "tcia_gbm_lgg_series_with_modality_flags.csv"
)
DEFAULT_OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri" / "download_manifests"


def modality_priority_column(modality: str) -> str:
    return {
        "FLAIR": "series_is_flair",
        "T1": "series_is_t1",
        "T1CE": "series_is_t1ce",
        "T2": "series_is_t2",
    }[modality]


def choose_series(sub: pd.DataFrame, modality: str) -> dict[str, Any] | None:
    flag = modality_priority_column(modality)
    candidates = sub[sub[flag].astype(bool)].copy()
    if candidates.empty:
        return None
    candidates["image_count_numeric"] = pd.to_numeric(candidates.get("Image Count", 0), errors="coerce").fillna(0)
    candidates["file_size_numeric"] = pd.to_numeric(candidates.get("File Size", 0), errors="coerce").fillna(0)
    # Prefer richer 3D/volumetric-looking series but keep the rule deterministic.
    candidates = candidates.sort_values(
        ["image_count_numeric", "file_size_numeric", "Series Date", "Series Instance UID"],
        ascending=[False, False, True, True],
    )
    row = candidates.iloc[0].to_dict()
    return {
        "patient_barcode": row.get("patient_barcode", ""),
        "mri_project_id": row.get("project_id_mri", ""),
        "modality": modality,
        "series_instance_uid": row.get("Series Instance UID", ""),
        "study_instance_uid": row.get("Study Instance UID", ""),
        "series_description": row.get("Series Description", ""),
        "protocol_name": row.get("Protocol Name", ""),
        "study_description": row.get("Study Description", ""),
        "series_date": row.get("Series Date", ""),
        "image_count": int(pd.to_numeric(row.get("Image Count", 0), errors="coerce") or 0),
        "file_size": int(pd.to_numeric(row.get("File Size", 0), errors="coerce") or 0),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a TCIA/NBIA series download manifest for complete four-modality patients.")
    parser.add_argument("--match", type=Path, default=DEFAULT_MATCH)
    parser.add_argument("--series", type=Path, default=DEFAULT_SERIES)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    match = pd.read_csv(args.match)
    series = pd.read_csv(args.series)
    complete = match[match["complete_brats4"].astype(bool)].copy()
    patients = sorted(complete["patient_barcode"].astype(str).unique())

    rows: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for patient in patients:
        sub = series[series["patient_barcode"].astype(str).eq(patient)].copy()
        for modality in ["T1", "T1CE", "T2", "FLAIR"]:
            selected = choose_series(sub, modality)
            if selected is None:
                missing.append({"patient_barcode": patient, "modality": modality})
            else:
                rows.append(selected)
    manifest = pd.DataFrame(rows)
    manifest.to_csv(args.outdir / "tcga_73_complete_brats4_selected_series_manifest.csv", index=False, encoding="utf-8-sig")
    manifest[["series_instance_uid"]].drop_duplicates().to_csv(
        args.outdir / "tcga_73_complete_brats4_series_uids.txt",
        index=False,
        header=False,
        encoding="utf-8",
    )
    summary = {
        "n_complete_brats4_patients": int(len(patients)),
        "n_selected_series": int(len(manifest)),
        "n_unique_series": int(manifest["series_instance_uid"].nunique()) if len(manifest) else 0,
        "expected_series": int(len(patients) * 4),
        "missing": missing,
        "modality_counts": manifest["modality"].value_counts().to_dict() if len(manifest) else {},
        "output_manifest": str(args.outdir / "tcga_73_complete_brats4_selected_series_manifest.csv"),
    }
    (args.outdir / "tcga_73_complete_brats4_selected_series_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
