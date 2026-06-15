#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RNASEQ_MANIFEST = (
    ROOT / "data" / "tcga_brain_tumor_expression" / "tcga_gbm_lgg_gdc_star_counts_manifest.csv"
)
DEFAULT_GBM_DIGEST = ROOT / "data" / "tcia_tcga_glioma_mri" / "manifests" / "TCGA-GBM_nbia_digest.xlsx"
DEFAULT_LGG_DIGEST = ROOT / "data" / "tcia_tcga_glioma_mri" / "manifests" / "TCGA-LGG_nbia_digest.xlsx"
DEFAULT_OUTDIR = ROOT / "results" / "tcga_rnaseq_tcia_mri_collection_match_20260605"


BARCODE_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", re.IGNORECASE)


def patient_barcode(value: Any) -> str:
    text = str(value or "").strip().upper().replace("_", "-")
    match = BARCODE_RE.search(text)
    if match:
        return match.group(1).upper()
    parts = text.split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 and parts[0] == "TCGA" else text


def normalize_text(value: Any) -> str:
    return str(value or "").strip().upper()


def classify_series(row: pd.Series) -> dict[str, bool]:
    text = " ".join(
        normalize_text(row.get(col, ""))
        for col in ["Series Description", "Protocol Name", "Study Description"]
        if col in row.index
    )
    has_flair = "FLAIR" in text or "FLR" in text
    has_t2 = bool(re.search(r"(^|[^A-Z0-9])T2([^A-Z0-9]|$)", text)) or "FSE T2" in text
    has_t1 = bool(re.search(r"(^|[^A-Z0-9])T1([^A-Z0-9]|$)", text)) or "SPGR" in text or "MPRAGE" in text
    contrast_terms = ["POST", "GAD", "GD", "GADO", "CONTRAST", "T1C", "T1CE", "T1-GD", "T1GD"]
    precontrast_terms = ["PRE", "NONCONTRAST", "NON-CONTRAST", "WITHOUT"]
    has_contrast = any(term in text for term in contrast_terms)
    has_precontrast = any(term in text for term in precontrast_terms)
    has_t1ce = has_t1 and has_contrast and not ("PRE" in text and "POST" not in text)
    return {
        "series_is_flair": has_flair,
        "series_is_t1": has_t1,
        "series_is_t1ce": has_t1ce,
        "series_is_t2": has_t2,
        "series_mentions_contrast": has_contrast,
        "series_mentions_precontrast": has_precontrast,
    }


def load_rnaseq_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["patient_barcode"] = df["case_submitter_id"].map(patient_barcode)
    sample_patient = (
        df[["project_id", "case_submitter_id", "sample_submitter_id", "sample_type", "patient_barcode"]]
        .drop_duplicates()
        .copy()
    )
    return sample_patient


def load_tcia_digest(path: Path, project_id: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Metadata")
    df["project_id_mri"] = project_id
    df["patient_barcode"] = df["Patient ID"].map(patient_barcode)
    flags = df.apply(classify_series, axis=1, result_type="expand")
    return pd.concat([df, flags], axis=1)


def compact_unique(values: pd.Series, max_items: int = 30) -> str:
    items = [str(x).strip() for x in values.dropna().astype(str).tolist() if str(x).strip()]
    unique = sorted(set(items))
    if len(unique) > max_items:
        return " | ".join(unique[:max_items]) + f" | ...(+{len(unique) - max_items})"
    return " | ".join(unique)


def summarize_patient_mri(mri: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for patient, sub in mri.groupby("patient_barcode", sort=True):
        n_series = int(len(sub))
        n_studies = int(sub["Study Instance UID"].nunique()) if "Study Instance UID" in sub.columns else 0
        n_images = int(pd.to_numeric(sub.get("Image Count", 0), errors="coerce").fillna(0).sum())
        has_flair = bool(sub["series_is_flair"].any())
        has_t1 = bool(sub["series_is_t1"].any())
        has_t1ce = bool(sub["series_is_t1ce"].any())
        has_t2 = bool(sub["series_is_t2"].any())
        complete_brats4 = has_flair and has_t1 and has_t1ce and has_t2
        minimal_segmentation_ready = has_flair and (has_t1ce or has_t1) and has_t2
        if complete_brats4:
            priority = "A_complete_T1_T1ce_T2_FLAIR"
        elif minimal_segmentation_ready:
            priority = "B_minimal_FLAIR_T2_T1_or_T1ce"
        elif has_flair and (has_t1 or has_t1ce or has_t2):
            priority = "C_partial_multimodal"
        else:
            priority = "D_low_priority_review"
        rows.append(
            {
                "patient_barcode": patient,
                "mri_project_id": compact_unique(sub["project_id_mri"]),
                "n_mri_series": n_series,
                "n_mri_studies": n_studies,
                "n_mri_images": n_images,
                "has_flair": has_flair,
                "has_t1": has_t1,
                "has_t1ce": has_t1ce,
                "has_t2": has_t2,
                "complete_brats4": complete_brats4,
                "minimal_segmentation_ready": minimal_segmentation_ready,
                "download_priority": priority,
                "series_descriptions": compact_unique(sub.get("Series Description", pd.Series(dtype=str))),
                "protocol_names": compact_unique(sub.get("Protocol Name", pd.Series(dtype=str))),
                "study_descriptions": compact_unique(sub.get("Study Description", pd.Series(dtype=str))),
                "series_instance_uids": compact_unique(sub.get("Series Instance UID", pd.Series(dtype=str)), max_items=80),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Match downloaded TCGA RNA-seq samples to TCIA MRI collection metadata.")
    parser.add_argument("--rnaseq-manifest", type=Path, default=DEFAULT_RNASEQ_MANIFEST)
    parser.add_argument("--gbm-digest", type=Path, default=DEFAULT_GBM_DIGEST)
    parser.add_argument("--lgg-digest", type=Path, default=DEFAULT_LGG_DIGEST)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rnaseq = load_rnaseq_manifest(args.rnaseq_manifest)
    mri = pd.concat(
        [
            load_tcia_digest(args.gbm_digest, "TCGA-GBM"),
            load_tcia_digest(args.lgg_digest, "TCGA-LGG"),
        ],
        ignore_index=True,
    )
    mri_patient = summarize_patient_mri(mri)

    rnaseq_patient = (
        rnaseq.groupby(["patient_barcode", "project_id"], as_index=False)
        .agg(
            n_rnaseq_samples=("sample_submitter_id", "nunique"),
            rnaseq_sample_ids=("sample_submitter_id", lambda s: compact_unique(s, max_items=20)),
            rnaseq_sample_types=("sample_type", compact_unique),
        )
        .sort_values(["project_id", "patient_barcode"])
    )
    matched = rnaseq_patient.merge(mri_patient, on="patient_barcode", how="left")
    matched["has_mri_collection_match"] = matched["n_mri_series"].notna()
    bool_cols = ["has_flair", "has_t1", "has_t1ce", "has_t2", "complete_brats4", "minimal_segmentation_ready"]
    for col in bool_cols:
        matched[col] = matched[col].map(lambda value: bool(value) if pd.notna(value) else False)
    matched["n_mri_series"] = matched["n_mri_series"].fillna(0).astype(int)
    matched["n_mri_studies"] = matched["n_mri_studies"].fillna(0).astype(int)
    matched["n_mri_images"] = matched["n_mri_images"].fillna(0).astype(int)
    matched["download_priority"] = matched["download_priority"].fillna("E_no_mri_collection_match")

    mri.to_csv(args.outdir / "tcia_gbm_lgg_series_with_modality_flags.csv", index=False, encoding="utf-8-sig")
    mri_patient.to_csv(args.outdir / "tcia_gbm_lgg_patient_mri_summary.csv", index=False, encoding="utf-8-sig")
    matched.to_csv(args.outdir / "tcga_rnaseq_patient_to_tcia_mri_match.csv", index=False, encoding="utf-8-sig")
    matched[matched["has_mri_collection_match"]].to_csv(
        args.outdir / "tcga_rnaseq_patients_with_mri.csv", index=False, encoding="utf-8-sig"
    )
    matched[~matched["has_mri_collection_match"]].to_csv(
        args.outdir / "tcga_rnaseq_patients_without_mri.csv", index=False, encoding="utf-8-sig"
    )
    matched[matched["minimal_segmentation_ready"]].to_csv(
        args.outdir / "tcga_rnaseq_patients_recommended_for_download_and_segmentation.csv",
        index=False,
        encoding="utf-8-sig",
    )

    summary: dict[str, Any] = {
        "n_rnaseq_patients": int(len(rnaseq_patient)),
        "n_tcia_mri_patients": int(mri_patient["patient_barcode"].nunique()),
        "n_matched_patients": int(matched["has_mri_collection_match"].sum()),
        "matched_fraction": float(matched["has_mri_collection_match"].mean()) if len(matched) else None,
        "by_rnaseq_project": {},
        "download_priority_counts": matched["download_priority"].value_counts().to_dict(),
        "modality_counts_among_rnaseq_patients": {
            "has_flair": int(matched["has_flair"].sum()),
            "has_t1": int(matched["has_t1"].sum()),
            "has_t1ce": int(matched["has_t1ce"].sum()),
            "has_t2": int(matched["has_t2"].sum()),
            "complete_brats4": int(matched["complete_brats4"].sum()),
            "minimal_segmentation_ready": int(matched["minimal_segmentation_ready"].sum()),
        },
        "outputs": {
            "matched_table": str(args.outdir / "tcga_rnaseq_patient_to_tcia_mri_match.csv"),
            "recommended_for_download": str(
                args.outdir / "tcga_rnaseq_patients_recommended_for_download_and_segmentation.csv"
            ),
            "series_flags": str(args.outdir / "tcia_gbm_lgg_series_with_modality_flags.csv"),
        },
    }
    for project, sub in matched.groupby("project_id", sort=True):
        summary["by_rnaseq_project"][project] = {
            "n_rnaseq_patients": int(len(sub)),
            "n_matched_patients": int(sub["has_mri_collection_match"].sum()),
            "matched_fraction": float(sub["has_mri_collection_match"].mean()) if len(sub) else None,
            "complete_brats4": int(sub["complete_brats4"].sum()),
            "minimal_segmentation_ready": int(sub["minimal_segmentation_ready"].sum()),
        }
    (args.outdir / "tcga_rnaseq_tcia_mri_match_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
