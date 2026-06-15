#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "mri_label_tools"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import nibabel as nib  # noqa: E402


AUDIT = ROOT / "results" / "brats_tcga_lgg_training_65_audit_20260609" / "brats_tcga_lgg_training_65_patient_audit.csv"
QC = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_qc_audit_20260612" / "patient_level_mri_atlas_truth_qc.csv"
ATLAS = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_evaluation_20260609" / "sri24_tzo116_resampled_to_brats.nii.gz"
OUTDIR = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_qc_audit_20260612" / "montages"
PATIENTS = ["TCGA-HT-7680", "TCGA-HT-7879", "TCGA-CS-6186", "TCGA-DU-7008", "TCGA-CS-6188", "TCGA-DU-8167"]


def normalize(image: np.ndarray) -> np.ndarray:
    nonzero = image[image > 0]
    if not len(nonzero):
        return image
    low, high = np.quantile(nonzero, [0.01, 0.99])
    return np.clip((image - low) / max(high - low, 1e-8), 0, 1)


def slice_axis(array: np.ndarray, axis: int, index: int) -> np.ndarray:
    return np.take(array, index, axis=axis).T


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(AUDIT).set_index("patient_barcode")
    qc = pd.read_csv(QC).set_index("patient_barcode")
    atlas = np.asanyarray(nib.load(str(ATLAS)).dataobj).astype(np.int16)
    for patient in PATIENTS:
        row = audit.loc[patient]
        q = qc.loc[patient]
        image = normalize(np.asanyarray(nib.load(row["flair_path"]).dataobj).astype(float))
        seg = np.asanyarray(nib.load(row["preferred_segmentation"]).dataobj) > 0
        coords = np.argwhere(seg)
        center = np.rint(coords.mean(axis=0)).astype(int)
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        for ax, axis, title in zip(axes, (0, 1, 2), ("Sagittal", "Coronal", "Axial")):
            base = slice_axis(image, axis, int(center[axis]))
            tumor = slice_axis(seg.astype(int), axis, int(center[axis]))
            labels = slice_axis(atlas, axis, int(center[axis]))
            ax.imshow(base, cmap="gray", origin="lower")
            if np.any(labels):
                ax.contour(labels > 0, levels=[0.5], colors="#2A78B8", linewidths=0.45, alpha=0.8)
            if np.any(tumor):
                ax.contour(tumor, levels=[0.5], colors="#E33D3D", linewidths=1.2)
            ax.set_title(title)
            ax.axis("off")
        fig.suptitle(
            f"{patient} | direct={q.direct_network} | filled={q.filled_network} | "
            f"coverage={q.atlas_direct_overlap_fraction:.1%}",
            fontsize=11,
        )
        fig.tight_layout()
        fig.savefig(OUTDIR / f"{patient}_flair_seg_atlas.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
    print(OUTDIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
