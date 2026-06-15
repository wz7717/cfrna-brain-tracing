#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SELECTED = (
    ROOT
    / "data"
    / "tcia_tcga_glioma_mri"
    / "download_manifests"
    / "tcga_73_complete_brats4_selected_series_manifest.csv"
)
DEFAULT_OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri" / "tcia_data_retriever_manifests"


def main() -> int:
    parser = argparse.ArgumentParser(description="Export TCIA Data Retriever manifest spreadsheets for 73 complete four-modality patients.")
    parser.add_argument("--selected", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    selected = pd.read_csv(args.selected)
    # Data Retriever beta accepts CSV/TSV/XLSX. Keep standard DICOM/TCIA field names.
    series = selected.rename(
        columns={
            "patient_barcode": "Patient ID",
            "mri_project_id": "Collection",
            "study_instance_uid": "Study Instance UID",
            "series_instance_uid": "Series Instance UID",
            "series_description": "Series Description",
            "protocol_name": "Protocol Name",
            "image_count": "Image Count",
            "file_size": "File Size",
        }
    )[
        [
            "Collection",
            "Patient ID",
            "modality",
            "Study Instance UID",
            "Series Instance UID",
            "Series Description",
            "Protocol Name",
            "Image Count",
            "File Size",
        ]
    ].copy()
    series = series.rename(columns={"modality": "Selected Modality"})
    unique_series = series.drop_duplicates("Series Instance UID").copy()
    patient = (
        series.groupby(["Collection", "Patient ID"], as_index=False)
        .agg(
            n_selected_modality_slots=("Selected Modality", "count"),
            n_unique_series=("Series Instance UID", "nunique"),
            selected_modalities=("Selected Modality", lambda s: " | ".join(s.astype(str).tolist())),
            series_instance_uids=("Series Instance UID", lambda s: " | ".join(pd.unique(s.astype(str)).tolist())),
        )
        .sort_values(["Collection", "Patient ID"])
    )

    series_csv = args.outdir / "tcga_73_complete_brats4_tcia_data_retriever_series.csv"
    series_tsv = args.outdir / "tcga_73_complete_brats4_tcia_data_retriever_series.tsv"
    series_xlsx = args.outdir / "tcga_73_complete_brats4_tcia_data_retriever_series.xlsx"
    unique_csv = args.outdir / "tcga_73_complete_brats4_tcia_data_retriever_unique_series.csv"
    patient_csv = args.outdir / "tcga_73_complete_brats4_patient_download_audit.csv"
    series.to_csv(series_csv, index=False, encoding="utf-8-sig")
    series.to_csv(series_tsv, sep="\t", index=False, encoding="utf-8")
    series.to_excel(series_xlsx, index=False)
    unique_series.to_csv(unique_csv, index=False, encoding="utf-8-sig")
    patient.to_csv(patient_csv, index=False, encoding="utf-8-sig")

    # A simple one-UID-per-line helper is useful for manual portal cart searches or copy/paste.
    uid_txt = args.outdir / "tcga_73_complete_brats4_unique_series_instance_uids.txt"
    uid_txt.write_text("\n".join(unique_series["Series Instance UID"].astype(str).tolist()) + "\n", encoding="utf-8")

    summary: dict[str, Any] = {
        "n_patients": int(patient["Patient ID"].nunique()),
        "n_modality_slots": int(len(series)),
        "n_unique_series": int(unique_series["Series Instance UID"].nunique()),
        "by_collection": patient["Collection"].value_counts().to_dict(),
        "outputs": {
            "series_csv": str(series_csv),
            "series_tsv": str(series_tsv),
            "series_xlsx": str(series_xlsx),
            "unique_series_csv": str(unique_csv),
            "patient_audit_csv": str(patient_csv),
            "series_uid_txt": str(uid_txt),
        },
        "note": (
            "Use the CSV/XLSX with the new TCIA Data Retriever, or use the UID list to recreate a cart in the "
            "TCIA Radiology Portal after logging in for restricted TCGA-GBM/LGG access. A classic .tcia manifest "
            "for restricted collections normally must be exported by the portal session."
        ),
    }
    (args.outdir / "tcga_73_complete_brats4_tcia_data_retriever_manifest_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
