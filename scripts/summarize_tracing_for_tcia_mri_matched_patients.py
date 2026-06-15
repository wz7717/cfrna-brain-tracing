#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TRACING = (
    ROOT
    / "results"
    / "tcga_gbm_lgg_sample_mri_label_tracing_20260605"
    / "tcga_gbm_lgg_sample_mri_label_tracing_summary.csv"
)
DEFAULT_MATCH = (
    ROOT
    / "results"
    / "tcga_rnaseq_tcia_mri_collection_match_20260605"
    / "tcga_rnaseq_patient_to_tcia_mri_match.csv"
)
DEFAULT_OUTDIR = ROOT / "results" / "tcga_tracing_tcia_mri_matched_patients_20260605"


def counts(series: pd.Series) -> dict[str, int]:
    return {str(k): int(v) for k, v in series.fillna("").astype(str).value_counts().items() if str(k)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize tracing outputs for RNA-seq patients matched to TCIA MRI metadata.")
    parser.add_argument("--tracing", type=Path, default=DEFAULT_TRACING)
    parser.add_argument("--match", type=Path, default=DEFAULT_MATCH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    tracing = pd.read_csv(args.tracing)
    match = pd.read_csv(args.match)
    matched_patients = match[match["has_mri_collection_match"].astype(bool)].copy()
    merged = tracing.merge(
        matched_patients[
            [
                "patient_barcode",
                "mri_project_id",
                "n_mri_series",
                "has_flair",
                "has_t1",
                "has_t1ce",
                "has_t2",
                "complete_brats4",
                "minimal_segmentation_ready",
                "download_priority",
                "series_instance_uids",
                "series_descriptions",
            ]
        ],
        on="patient_barcode",
        how="inner",
    )
    recommended = merged[merged["minimal_segmentation_ready"].astype(bool)].copy()
    complete = merged[merged["complete_brats4"].astype(bool)].copy()

    merged.to_csv(args.outdir / "tcga_tracing_for_156_mri_matched_patients.csv", index=False, encoding="utf-8-sig")
    recommended.to_csv(
        args.outdir / "tcga_tracing_for_105_segmentation_ready_patients.csv",
        index=False,
        encoding="utf-8-sig",
    )
    complete.to_csv(
        args.outdir / "tcga_tracing_for_73_complete_brats4_patients.csv",
        index=False,
        encoding="utf-8-sig",
    )

    summary: dict[str, Any] = {
        "interpretation": (
            "These are MRI-collection-matched patients, not tumor-location-labeled patients. "
            "Accuracy against lobe/network requires segmentation-derived or curated tumor location labels."
        ),
        "n_tracing_samples_total": int(len(tracing)),
        "n_mri_matched_tracing_samples": int(len(merged)),
        "n_mri_matched_patients": int(merged["patient_barcode"].nunique()),
        "n_segmentation_ready_samples": int(len(recommended)),
        "n_segmentation_ready_patients": int(recommended["patient_barcode"].nunique()),
        "n_complete_brats4_samples": int(len(complete)),
        "n_complete_brats4_patients": int(complete["patient_barcode"].nunique()),
        "mri_matched_project_counts": counts(merged["project_id"]),
        "mri_matched_network_top1_counts": counts(merged["network_top1"]),
        "mri_matched_broad_top1_counts": counts(merged["predicted_broad_top1"]),
        "mri_matched_region_top1_counts": counts(merged["region_top1"]),
        "segmentation_ready_network_top1_counts": counts(recommended["network_top1"]),
        "complete_brats4_network_top1_counts": counts(complete["network_top1"]),
        "outputs": {
            "mri_matched_tracing": str(args.outdir / "tcga_tracing_for_156_mri_matched_patients.csv"),
            "segmentation_ready_tracing": str(args.outdir / "tcga_tracing_for_105_segmentation_ready_patients.csv"),
            "complete_brats4_tracing": str(args.outdir / "tcga_tracing_for_73_complete_brats4_patients.csv"),
        },
    }
    (args.outdir / "tcga_tracing_tcia_mri_matched_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
