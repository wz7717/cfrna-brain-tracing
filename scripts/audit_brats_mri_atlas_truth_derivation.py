#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "mri_label_tools"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import nibabel as nib  # noqa: E402
from nibabel.orientations import aff2axcodes  # noqa: E402
from nibabel.processing import resample_from_to  # noqa: E402

from scripts.evaluate_brats_tcga_lgg_65_mri_truth import (  # noqa: E402
    anatomy_for_region,
    parse_tzo_lut,
)


AUDIT = ROOT / "results" / "brats_tcga_lgg_training_65_audit_20260609" / "brats_tcga_lgg_training_65_patient_audit.csv"
DERIVED = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_evaluation_20260609" / "brats_tcga_lgg_65_mri_truth_and_predictions.csv"
ATLAS = ROOT / "data" / "atlases" / "sri24" / "labels" / "sri24" / "tzo116plus.nii"
ATLAS_LUT = ROOT / "data" / "atlases" / "sri24" / "labels" / "sri24" / "SRI24-tzo116plus.txt"
SRI_T1 = ROOT / "data" / "atlases" / "sri24" / "anatomy" / "sri24" / "spgr.nii"
OUTDIR = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_qc_audit_20260612"


def resample_labels(reference: nib.spatialimages.SpatialImage) -> np.ndarray:
    image = nib.load(str(ATLAS))
    data = np.squeeze(np.asanyarray(image.dataobj)).astype(np.int16)
    image3d = nib.Nifti1Image(data, image.affine, image.header)
    result = resample_from_to(image3d, (reference.shape[:3], reference.affine), order=0)
    labels = np.rint(np.asanyarray(result.dataobj)).astype(np.int16)
    labels[(labels < 1) | (labels > 116)] = 0
    return labels


def top_dimension(labels: np.ndarray, info: dict[int, dict[str, str]], dimension: str) -> tuple[str, float]:
    labels = labels[labels > 0]
    if not len(labels):
        return "", float("nan")
    values: dict[str, int] = {}
    for label, count in zip(*np.unique(labels, return_counts=True)):
        name = info[int(label)][dimension]
        values[name] = values.get(name, 0) + int(count)
    ordered = sorted(values.items(), key=lambda item: (-item[1], item[0]))
    total = sum(values.values())
    return ordered[0][0], ordered[0][1] / total


def normalized_correlation(left: np.ndarray, right: np.ndarray) -> float:
    mask = (left > 0) & (right != 0)
    if mask.sum() < 100:
        return float("nan")
    x = left[mask].astype(float)
    y = right[mask].astype(float)
    x = np.clip(x, *np.quantile(x, [0.01, 0.99]))
    y = np.clip(y, *np.quantile(y, [0.01, 0.99]))
    return float(np.corrcoef(x, y)[0, 1])


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(AUDIT)
    prior = pd.read_csv(DERIVED).set_index("patient_barcode")
    lut = parse_tzo_lut(ATLAS_LUT)
    info = {}
    for label, region in lut.items():
        lobe, broad, network = anatomy_for_region(region)
        info[label] = {"region": region, "lobe": lobe, "broad": broad, "network": network}

    reference = nib.load(str(audit.iloc[0]["t1_path"]))
    atlas = resample_labels(reference)
    distance, indices = distance_transform_edt(atlas == 0, return_indices=True)
    filled = atlas[tuple(indices)]
    direct_flip = np.squeeze(np.asanyarray(nib.load(str(ATLAS)).dataobj)).astype(np.int16)[:, ::-1, :]
    direct_flip[(direct_flip < 1) | (direct_flip > 116)] = 0
    atlas_resample_equals_y_flip = bool(np.array_equal(atlas, direct_flip))

    sri = np.squeeze(np.asanyarray(nib.load(str(SRI_T1)).dataobj)).astype(float)
    rows = []
    for item in audit.itertuples(index=False):
        images = {
            "t1": nib.load(item.t1_path),
            "t1ce": nib.load(item.t1ce_path),
            "t2": nib.load(item.t2_path),
            "flair": nib.load(item.flair_path),
            "seg": nib.load(item.preferred_segmentation),
        }
        shapes_equal = len({tuple(image.shape[:3]) for image in images.values()}) == 1
        affines_equal = all(np.allclose(images["t1"].affine, image.affine, atol=1e-5) for image in images.values())
        seg = np.rint(np.asanyarray(images["seg"].dataobj)).astype(np.int16)
        tumor = seg > 0
        coords = np.argwhere(tumor)
        center_float = coords.mean(axis=0)
        center_round = np.rint(center_float).astype(int)
        center_inside_tumor = bool(tumor[tuple(center_round)])
        tumor_distances = distance[tumor]
        direct_labels = atlas[tumor]
        filled_labels = filled[tumor]
        direct_network, direct_network_fraction = top_dimension(direct_labels, info, "network")
        filled_network, filled_network_fraction = top_dimension(filled_labels, info, "network")
        direct_lobe, direct_lobe_fraction = top_dimension(direct_labels, info, "lobe")
        filled_lobe, filled_lobe_fraction = top_dimension(filled_labels, info, "lobe")
        direct_broad, direct_broad_fraction = top_dimension(direct_labels, info, "broad")
        filled_broad, filled_broad_fraction = top_dimension(filled_labels, info, "broad")
        prior_row = prior.loc[item.patient_barcode]
        t1 = np.asanyarray(images["t1"].dataobj).astype(float)
        correlations = {
            "identity": normalized_correlation(sri, t1),
            "flip_y": normalized_correlation(sri[:, ::-1, :], t1),
            "flip_x": normalized_correlation(sri[::-1, :, :], t1),
            "flip_xy": normalized_correlation(sri[::-1, ::-1, :], t1),
        }
        best_orientation = max(correlations, key=lambda key: correlations[key])
        rows.append(
            {
                "patient_barcode": item.patient_barcode,
                "segmentation_source": "manual_corrected" if bool(item.has_manual_segmentation) else "automatic",
                "all_modalities_same_shape": shapes_equal,
                "all_modalities_same_affine": affines_equal,
                "t1_shape": "x".join(map(str, images["t1"].shape[:3])),
                "t1_orientation": "".join(aff2axcodes(images["t1"].affine)),
                "seg_orientation": "".join(aff2axcodes(images["seg"].affine)),
                "seg_unique_labels": "|".join(map(str, np.unique(seg))),
                "tumor_voxels_recomputed": int(tumor.sum()),
                "tumor_voxels_matches_prior": int(tumor.sum()) == int(prior_row["tumor_voxels"]),
                "atlas_direct_overlap_fraction": float(np.mean(direct_labels > 0)),
                "nearest_fill_distance_median_vox": float(np.median(tumor_distances)),
                "nearest_fill_distance_p95_vox": float(np.quantile(tumor_distances, 0.95)),
                "nearest_fill_distance_max_vox": float(np.max(tumor_distances)),
                "fraction_fill_distance_gt_3vox": float(np.mean(tumor_distances > 3)),
                "fraction_fill_distance_gt_5vox": float(np.mean(tumor_distances > 5)),
                "center_i": int(center_round[0]),
                "center_j": int(center_round[1]),
                "center_k": int(center_round[2]),
                "center_inside_tumor": center_inside_tumor,
                "center_distance_to_tumor_vox": float(distance_transform_edt(~tumor)[tuple(center_round)]),
                "direct_network": direct_network,
                "direct_network_fraction_labeled_only": direct_network_fraction,
                "filled_network": filled_network,
                "filled_network_fraction_all_tumor": filled_network_fraction,
                "network_direct_vs_filled_match": direct_network == filled_network,
                "filled_network_matches_prior": filled_network == prior_row["whole_tumor_network_dominant"],
                "direct_lobe": direct_lobe,
                "direct_lobe_fraction_labeled_only": direct_lobe_fraction,
                "filled_lobe": filled_lobe,
                "filled_lobe_fraction_all_tumor": filled_lobe_fraction,
                "lobe_direct_vs_filled_match": direct_lobe == filled_lobe,
                "filled_lobe_matches_prior": filled_lobe == prior_row["whole_tumor_lobe_dominant"],
                "direct_broad": direct_broad,
                "direct_broad_fraction_labeled_only": direct_broad_fraction,
                "filled_broad": filled_broad,
                "filled_broad_fraction_all_tumor": filled_broad_fraction,
                "broad_direct_vs_filled_match": direct_broad == filled_broad,
                "filled_broad_matches_prior": filled_broad == prior_row["whole_tumor_broad_dominant"],
                "sri_t1_corr_identity": correlations["identity"],
                "sri_t1_corr_flip_y": correlations["flip_y"],
                "sri_t1_corr_flip_x": correlations["flip_x"],
                "sri_t1_corr_flip_xy": correlations["flip_xy"],
                "best_array_orientation_by_corr": best_orientation,
            }
        )

    frame = pd.DataFrame(rows)
    frame["qc_flag"] = (
        (~frame["all_modalities_same_shape"])
        | (~frame["all_modalities_same_affine"])
        | (~frame["tumor_voxels_matches_prior"])
        | (~frame["filled_network_matches_prior"])
        | (~frame["filled_lobe_matches_prior"])
        | (~frame["filled_broad_matches_prior"])
        | (frame["atlas_direct_overlap_fraction"] < 0.50)
        | (frame["fraction_fill_distance_gt_5vox"] > 0.10)
        | (~frame["center_inside_tumor"])
    )
    frame.to_csv(OUTDIR / "patient_level_mri_atlas_truth_qc.csv", index=False, encoding="utf-8-sig")
    disagreements = frame[
        (~frame["network_direct_vs_filled_match"])
        | (~frame["lobe_direct_vs_filled_match"])
        | (~frame["broad_direct_vs_filled_match"])
        | (~frame["center_inside_tumor"])
        | (frame["fraction_fill_distance_gt_5vox"] > 0.10)
    ].copy()
    disagreements.to_csv(OUTDIR / "truth_derivation_review_cases.csv", index=False, encoding="utf-8-sig")

    summary = {
        "n_patients": int(len(frame)),
        "n_manual_segmentations": int(frame["segmentation_source"].eq("manual_corrected").sum()),
        "n_same_shape_all_modalities": int(frame["all_modalities_same_shape"].sum()),
        "n_same_affine_all_modalities": int(frame["all_modalities_same_affine"].sum()),
        "segmentation_label_sets": frame["seg_unique_labels"].value_counts().to_dict(),
        "atlas_resample_equals_exact_y_axis_flip": atlas_resample_equals_y_flip,
        "t1_orientation": frame["t1_orientation"].value_counts().to_dict(),
        "median_direct_atlas_overlap_fraction": float(frame["atlas_direct_overlap_fraction"].median()),
        "minimum_direct_atlas_overlap_fraction": float(frame["atlas_direct_overlap_fraction"].min()),
        "n_direct_overlap_below_50pct": int((frame["atlas_direct_overlap_fraction"] < 0.50).sum()),
        "n_more_than_10pct_tumor_farther_than_5vox_from_label": int(
            (frame["fraction_fill_distance_gt_5vox"] > 0.10).sum()
        ),
        "n_center_outside_tumor": int((~frame["center_inside_tumor"]).sum()),
        "n_network_changed_by_nearest_fill": int((~frame["network_direct_vs_filled_match"]).sum()),
        "n_lobe_changed_by_nearest_fill": int((~frame["lobe_direct_vs_filled_match"]).sum()),
        "n_broad_changed_by_nearest_fill": int((~frame["broad_direct_vs_filled_match"]).sum()),
        "n_recomputed_truth_mismatch_prior": {
            "network": int((~frame["filled_network_matches_prior"]).sum()),
            "lobe": int((~frame["filled_lobe_matches_prior"]).sum()),
            "broad": int((~frame["filled_broad_matches_prior"]).sum()),
        },
        "best_sri_t1_array_orientation_by_correlation": frame[
            "best_array_orientation_by_corr"
        ].value_counts().to_dict(),
        "n_qc_flagged": int(frame["qc_flag"].sum()),
    }
    (OUTDIR / "mri_atlas_truth_qc_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
