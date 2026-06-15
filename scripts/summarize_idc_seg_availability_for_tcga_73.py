#!/usr/bin/env python
from __future__ import annotations

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


PATIENTS = (
    ROOT
    / "results"
    / "tcga_rnaseq_tcia_mri_collection_match_20260605"
    / "tcga_rnaseq_patients_recommended_for_download_and_segmentation.csv"
)
OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri" / "download_manifests"


def sql_list(values: list[str]) -> str:
    return ",".join("'" + v.replace("'", "''") + "'" for v in values)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    recommended = pd.read_csv(PATIENTS)
    complete = recommended[recommended["complete_brats4"].astype(bool)].copy()
    patients = sorted(complete["patient_barcode"].astype(str).unique())
    client = IDCClient()
    q = f"""
    SELECT collection_id, PatientID, StudyInstanceUID, SeriesInstanceUID,
           Modality, SeriesDescription, instanceCount, series_size_MB,
           series_aws_url, crdc_series_uuid
    FROM index
    WHERE PatientID IN ({sql_list(patients)})
      AND Modality IN ('SEG', 'ANN', 'OT', 'MR')
    ORDER BY PatientID, Modality, SeriesInstanceUID
    """
    rows = client.sql_query(q)
    rows.to_csv(OUTDIR / "tcga_73_complete_brats4_idc_seg_ann_ot_mr_availability.csv", index=False, encoding="utf-8-sig")
    patient = (
        rows.pivot_table(index="PatientID", columns="Modality", values="SeriesInstanceUID", aggfunc="count", fill_value=0)
        .reset_index()
        if len(rows)
        else pd.DataFrame(columns=["PatientID"])
    )
    all_patient_df = pd.DataFrame({"PatientID": patients})
    patient = all_patient_df.merge(patient, on="PatientID", how="left").fillna(0)
    for col in ["MR", "SEG", "ANN", "OT"]:
        if col not in patient.columns:
            patient[col] = 0
        patient[col] = patient[col].astype(int)
    patient["has_idc_mr"] = patient["MR"].gt(0)
    patient["has_idc_seg"] = patient["SEG"].gt(0)
    patient["has_idc_ann"] = patient["ANN"].gt(0)
    patient["has_idc_ot"] = patient["OT"].gt(0)
    patient.to_csv(OUTDIR / "tcga_73_complete_brats4_idc_patient_seg_availability.csv", index=False, encoding="utf-8-sig")
    summary: dict[str, Any] = {
        "idc_version": client.get_idc_version(),
        "n_complete_brats4_patients": int(len(patient)),
        "n_has_idc_mr": int(patient["has_idc_mr"].sum()),
        "n_has_idc_seg": int(patient["has_idc_seg"].sum()),
        "n_has_idc_ann": int(patient["has_idc_ann"].sum()),
        "n_has_idc_ot": int(patient["has_idc_ot"].sum()),
        "modality_series_counts": {
            col: int(rows[rows["Modality"].eq(col)]["SeriesInstanceUID"].nunique()) if len(rows) else 0
            for col in ["MR", "SEG", "ANN", "OT"]
        },
        "outputs": {
            "series_availability": str(OUTDIR / "tcga_73_complete_brats4_idc_seg_ann_ot_mr_availability.csv"),
            "patient_availability": str(OUTDIR / "tcga_73_complete_brats4_idc_patient_seg_availability.csv"),
        },
    }
    (OUTDIR / "tcga_73_complete_brats4_idc_seg_availability_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
