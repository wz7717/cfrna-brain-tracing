from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = (
    ROOT
    / "data"
    / "brats_tcga_lgg_training_65"
    / "extracted"
    / "Pre-operative_TCGA_LGG_NIfTI_and_Segmentations"
)
RNA_SUMMARY = (
    ROOT
    / "results"
    / "tcga_gbm_lgg_sample_mri_label_tracing_20260605"
    / "tcga_gbm_lgg_sample_mri_label_tracing_summary.csv"
)
OUTDIR = ROOT / "results" / "brats_tcga_lgg_training_65_audit_20260609"


def classify_file(name: str) -> str:
    lower = name.lower()
    if "glistrboost_manuallycorrected" in lower:
        return "segmentation_manual"
    if "glistrboost" in lower:
        return "segmentation_auto"
    if "_flair." in lower:
        return "flair"
    if "_t1gd." in lower:
        return "t1ce"
    if "_t1." in lower:
        return "t1"
    if "_t2." in lower:
        return "t2"
    return "other"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rna = pd.read_csv(RNA_SUMMARY, dtype=str)
    lgg_rna = rna.loc[rna["project_id"] == "TCGA-LGG"].copy()
    rna_patients = set(lgg_rna["patient_barcode"])

    rows: list[dict[str, object]] = []
    for patient_dir in sorted(DATA_ROOT.iterdir()):
        if not patient_dir.is_dir():
            continue
        files = sorted(patient_dir.glob("*.nii.gz"))
        by_type = {classify_file(path.name): path for path in files}
        patient = patient_dir.name
        rows.append(
            {
                "patient_barcode": patient,
                "n_nifti_files": len(files),
                "has_t1": "t1" in by_type,
                "has_t1ce": "t1ce" in by_type,
                "has_t2": "t2" in by_type,
                "has_flair": "flair" in by_type,
                "has_auto_segmentation": "segmentation_auto" in by_type,
                "has_manual_segmentation": "segmentation_manual" in by_type,
                "preferred_segmentation": (
                    str(by_type["segmentation_manual"])
                    if "segmentation_manual" in by_type
                    else str(by_type.get("segmentation_auto", ""))
                ),
                "t1_path": str(by_type.get("t1", "")),
                "t1ce_path": str(by_type.get("t1ce", "")),
                "t2_path": str(by_type.get("t2", "")),
                "flair_path": str(by_type.get("flair", "")),
                "matches_lgg_rnaseq": patient in rna_patients,
            }
        )

    audit = pd.DataFrame(rows)
    audit.to_csv(OUTDIR / "brats_tcga_lgg_training_65_patient_audit.csv", index=False)
    matched = audit.loc[audit["matches_lgg_rnaseq"]].merge(
        lgg_rna,
        on="patient_barcode",
        how="left",
        validate="one_to_many",
    )
    matched.to_csv(
        OUTDIR / "brats_tcga_lgg_training_65_rnaseq_matched_tracing.csv",
        index=False,
    )

    summary = {
        "n_brats_patients": int(len(audit)),
        "n_complete_four_modality": int(
            audit[["has_t1", "has_t1ce", "has_t2", "has_flair"]].all(axis=1).sum()
        ),
        "n_auto_segmentations": int(audit["has_auto_segmentation"].sum()),
        "n_manual_segmentations": int(audit["has_manual_segmentation"].sum()),
        "n_matching_lgg_rnaseq_patients": int(audit["matches_lgg_rnaseq"].sum()),
        "n_matching_lgg_rnaseq_samples": int(len(matched)),
        "data_root": str(DATA_ROOT),
        "outputs": {
            "patient_audit": str(
                OUTDIR / "brats_tcga_lgg_training_65_patient_audit.csv"
            ),
            "rnaseq_matched_tracing": str(
                OUTDIR / "brats_tcga_lgg_training_65_rnaseq_matched_tracing.csv"
            ),
        },
    }
    (OUTDIR / "brats_tcga_lgg_training_65_audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
