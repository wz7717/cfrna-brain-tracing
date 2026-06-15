#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "idc_index_py312"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import pandas as pd  # noqa: E402
from idc_index import IDCClient  # noqa: E402


DEFAULT_SELECTED = (
    ROOT
    / "data"
    / "tcia_tcga_glioma_mri"
    / "download_manifests"
    / "tcga_73_complete_brats4_selected_series_manifest.csv"
)
DEFAULT_OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri" / "download_manifests"


def sql_quote(values: list[str]) -> str:
    return ",".join("'" + value.replace("'", "''") + "'" for value in values)


def main() -> int:
    parser = argparse.ArgumentParser(description="Match selected TCIA series to IDC index and write download manifest.")
    parser.add_argument("--selected", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    selected = pd.read_csv(args.selected)
    uids = sorted(selected["series_instance_uid"].dropna().astype(str).unique())
    client = IDCClient()
    query = f"""
    SELECT
      collection_id,
      PatientID,
      StudyInstanceUID,
      SeriesInstanceUID,
      Modality,
      SeriesDescription,
      instanceCount,
      series_size_MB,
      aws_bucket,
      crdc_series_uuid,
      series_aws_url
    FROM index
    WHERE SeriesInstanceUID IN ({sql_quote(uids)})
    """
    idc = client.sql_query(query)
    idc = idc.rename(columns={"PatientID": "idc_patient_id"})
    merged = selected.merge(
        idc,
        left_on="series_instance_uid",
        right_on="SeriesInstanceUID",
        how="left",
    )
    merged["found_in_idc"] = merged["SeriesInstanceUID"].notna()
    unique_status = (
        merged[["series_instance_uid", "found_in_idc"]]
        .drop_duplicates()
        .groupby("found_in_idc")
        .size()
        .to_dict()
    )
    patient_ready = (
        merged.groupby("patient_barcode")
        .agg(
            n_selected_slots=("modality", "count"),
            n_found_slots=("found_in_idc", "sum"),
            n_unique_found_series=("SeriesInstanceUID", "nunique"),
            total_idc_size_mb=("series_size_MB", "sum"),
        )
        .reset_index()
    )
    patient_ready["all_selected_slots_found_in_idc"] = patient_ready["n_found_slots"].eq(patient_ready["n_selected_slots"])
    merged.to_csv(args.outdir / "tcga_73_complete_brats4_selected_series_idc_manifest.csv", index=False, encoding="utf-8-sig")
    patient_ready.to_csv(args.outdir / "tcga_73_complete_brats4_idc_patient_download_summary.csv", index=False, encoding="utf-8-sig")
    found_series = merged[merged["found_in_idc"]][["series_instance_uid"]].drop_duplicates()
    found_series.to_csv(args.outdir / "tcga_73_complete_brats4_idc_series_uids.txt", index=False, header=False, encoding="utf-8")
    summary: dict[str, Any] = {
        "idc_version": client.get_idc_version(),
        "n_selected_slots": int(len(selected)),
        "n_unique_selected_series": int(len(uids)),
        "n_unique_series_found_in_idc": int(idc["SeriesInstanceUID"].nunique()) if len(idc) else 0,
        "unique_series_found_status": {str(k): int(v) for k, v in unique_status.items()},
        "n_patients": int(patient_ready["patient_barcode"].nunique()),
        "n_patients_all_selected_slots_found": int(patient_ready["all_selected_slots_found_in_idc"].sum()),
        "total_idc_size_mb_selected_slots": float(merged["series_size_MB"].fillna(0).sum()),
        "total_idc_size_mb_unique_found_series": float(idc.drop_duplicates("SeriesInstanceUID")["series_size_MB"].sum())
        if len(idc)
        else 0.0,
        "outputs": {
            "idc_manifest": str(args.outdir / "tcga_73_complete_brats4_selected_series_idc_manifest.csv"),
            "patient_summary": str(args.outdir / "tcga_73_complete_brats4_idc_patient_download_summary.csv"),
            "series_uid_list": str(args.outdir / "tcga_73_complete_brats4_idc_series_uids.txt"),
        },
    }
    (args.outdir / "tcga_73_complete_brats4_idc_manifest_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
