#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MRI_DIR = ROOT / "data" / "tcia_tcga_glioma_mri"
DEFAULT_OUT = DEFAULT_MRI_DIR / "tcia_tcga_glioma_mri_derived_labels.csv"


BARCODE_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", re.IGNORECASE)


def patient_barcode_from_path(path: Path) -> str:
    match = BARCODE_RE.search(str(path).replace("_", "-"))
    return match.group(1).upper() if match else path.stem.upper()


def load_nifti(path: Path) -> tuple[np.ndarray, np.ndarray]:
    try:
        import nibabel as nib
    except ImportError as exc:
        raise RuntimeError(
            "nibabel is required for NIfTI segmentation processing. Install nibabel or provide a prebuilt label CSV."
        ) from exc
    img = nib.load(str(path))
    return np.asarray(img.get_fdata()), np.asarray(img.affine, dtype=float)


def read_lut(path: Path | None) -> dict[int, dict[str, str]]:
    if path is None:
        return {}
    if path.suffix.lower() == ".json":
        raw = json.loads(path.read_text(encoding="utf-8"))
        rows = raw.values() if isinstance(raw, dict) else raw
        return {int(row["label"]): {k: str(v) for k, v in row.items()} for row in rows}
    df = pd.read_csv(path)
    if "label" not in df.columns:
        raise ValueError("atlas LUT must contain a 'label' column")
    return {
        int(row["label"]): {col: str(row[col]) for col in df.columns if pd.notna(row[col])}
        for _, row in df.iterrows()
    }


def top_label(mask: np.ndarray, atlas: np.ndarray, lut: dict[int, dict[str, str]]) -> dict[str, Any]:
    labels, counts = np.unique(atlas[mask].astype(int), return_counts=True)
    keep = labels > 0
    labels = labels[keep]
    counts = counts[keep]
    if len(labels) == 0:
        return {}
    order = np.argsort(counts)[::-1]
    label = int(labels[order[0]])
    info = dict(lut.get(label, {}))
    total = int(counts.sum())
    return {
        "atlas_label": label,
        "atlas_region": info.get("region", info.get("name", "")),
        "tumor_lobe": info.get("lobe", ""),
        "tumor_broad_anatomy": info.get("broad_anatomy", info.get("anatomy", "")),
        "tumor_network": info.get("network", ""),
        "atlas_overlap_voxels": int(counts[order[0]]),
        "atlas_overlap_fraction": float(counts[order[0]] / total) if total else float("nan"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build patient-level tumor/edema location labels from BraTS/TCIA segmentation NIfTI files."
    )
    parser.add_argument("--seg-dir", type=Path, default=DEFAULT_MRI_DIR / "segmentations")
    parser.add_argument("--seg-glob", default="**/*seg*.nii*")
    parser.add_argument("--atlas", type=Path, default=None, help="Atlas NIfTI in the same voxel space as segmentations.")
    parser.add_argument("--atlas-lut", type=Path, default=None, help="CSV/JSON with label, region, lobe, broad_anatomy, network.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    seg_paths = sorted(args.seg_dir.glob(args.seg_glob))
    if not seg_paths:
        raise FileNotFoundError(f"no segmentation files matched {args.seg_dir / args.seg_glob}")

    atlas = None
    lut: dict[int, dict[str, str]] = {}
    if args.atlas is not None:
        atlas, _ = load_nifti(args.atlas)
        lut = read_lut(args.atlas_lut)

    rows: list[dict[str, Any]] = []
    for seg_path in seg_paths:
        seg, affine = load_nifti(seg_path)
        tumor_mask = seg > 0
        edema_mask = np.isclose(seg, 2)
        if not np.any(tumor_mask):
            continue
        coords = np.argwhere(tumor_mask)
        center_voxel = coords.mean(axis=0)
        center_mm = np.append(center_voxel, 1.0) @ affine.T
        row: dict[str, Any] = {
            "patient_barcode": patient_barcode_from_path(seg_path),
            "segmentation_path": str(seg_path),
            "tumor_voxels": int(tumor_mask.sum()),
            "edema_voxels": int(edema_mask.sum()),
            "tumor_center_i": float(center_voxel[0]),
            "tumor_center_j": float(center_voxel[1]),
            "tumor_center_k": float(center_voxel[2]),
            "tumor_center_x": float(center_mm[0]),
            "tumor_center_y": float(center_mm[1]),
            "tumor_center_z": float(center_mm[2]),
            "tumor_hemisphere": "left" if float(center_mm[0]) < 0 else "right",
        }
        if atlas is not None:
            if atlas.shape != seg.shape:
                raise ValueError(f"atlas shape {atlas.shape} does not match {seg_path.name} shape {seg.shape}")
            row.update(top_label(tumor_mask, atlas, lut))
            edema_info = top_label(edema_mask, atlas, lut) if np.any(edema_mask) else {}
            for key, value in edema_info.items():
                row[f"edema_{key}"] = value
        rows.append(row)

    out = pd.DataFrame(rows).drop_duplicates("patient_barcode")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(json.dumps({"n_labeled_patients": int(len(out)), "output": str(args.out)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
