#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "mri_label_tools"
for path in (ROOT, VENDOR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import nibabel as nib  # noqa: E402
from nibabel.processing import resample_from_to  # noqa: E402

from scripts.evaluate_brats_tcga_lgg_65_mri_truth import anatomy_for_region, parse_tzo_lut  # noqa: E402


AUDIT = ROOT / "results" / "brats_tcga_lgg_training_65_audit_20260609" / "brats_tcga_lgg_training_65_patient_audit.csv"
OLD_TRUTH = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_evaluation_20260609" / "brats_tcga_lgg_65_mri_truth_and_predictions.csv"
ATLAS = ROOT / "data" / "atlases" / "sri24" / "labels" / "sri24" / "tzo116plus.nii"
LUT = ROOT / "data" / "atlases" / "sri24" / "labels" / "sri24" / "SRI24-tzo116plus.txt"
OUTDIR = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_corrected_20260612"


def resample_atlas(reference: nib.spatialimages.SpatialImage) -> np.ndarray:
    image = nib.load(str(ATLAS))
    data = np.squeeze(np.asanyarray(image.dataobj)).astype(np.int16)
    result = resample_from_to(
        nib.Nifti1Image(data, image.affine, image.header),
        (reference.shape[:3], reference.affine),
        order=0,
    )
    labels = np.rint(np.asanyarray(result.dataobj)).astype(np.int16)
    labels[(labels < 1) | (labels > 116)] = 0
    return labels


def summarize(labels: np.ndarray, info: dict[int, dict[str, str]]) -> dict[str, Any]:
    labels = labels[labels > 0]
    dimensions: dict[str, dict[str, int]] = {
        "region": defaultdict(int),
        "lobe": defaultdict(int),
        "broad": defaultdict(int),
        "network": defaultdict(int),
    }
    for label, count in zip(*np.unique(labels, return_counts=True)):
        for dimension in dimensions:
            dimensions[dimension][info[int(label)][dimension]] += int(count)
    result: dict[str, Any] = {}
    for dimension, counts in dimensions.items():
        ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        total = sum(counts.values())
        fractions = [(name, count / total) for name, count in ordered]
        candidates = [name for name, fraction in fractions if fraction >= 0.20]
        if not candidates:
            candidates = [fractions[0][0]]
        result[f"{dimension}_dominant"] = fractions[0][0]
        result[f"{dimension}_dominant_fraction"] = fractions[0][1]
        result[f"{dimension}_candidates"] = " | ".join(candidates)
        result[f"{dimension}_distribution"] = " | ".join(
            f"{name}:{fraction:.4f}" for name, fraction in fractions
        )
    return result


def split(value: Any) -> set[str]:
    return {item.strip() for item in str(value or "").split("|") if item.strip()}


def metric(frame: pd.DataFrame, level: str, topk: int, tolerant: bool) -> dict[str, Any]:
    truth_column = f"corrected_{level}_{'candidates' if tolerant else 'dominant'}"
    pred_column = "network_prediction_top1" if topk == 1 else "network_prediction_top3"
    rows = frame.copy()
    if level == "network":
        rows = rows[~rows["corrected_network_dominant"].eq("out_of_scope")]
        predictions = rows[pred_column].map(split)
    else:
        network_to_level = {}
        for network, values in (
            ("Cingulate gyrus", {"cingulate"}),
            ("Frontal (agranular frontal motor areas)", {"frontal"}),
            ("Hippocampal formation", {"temporal"} if level == "lobe" else {"medial_temporal", "temporal"}),
            ("Lateral Prefrontal Cortex", {"frontal"}),
            ("Occipital/Temporal", {"occipital", "temporal"} if level == "lobe" else {"occipital_temporal", "occipital", "temporal"}),
            ("Operculum/Insula", {"frontal", "insula"} if level == "lobe" else {"insula_operculum", "insula"}),
            ("Orbitomedial Prefrontal Cortex (OMPFC)", {"frontal"}),
            ("Parietal, and Parieto-occipital region", {"parietal", "occipital"} if level == "lobe" else {"parietal", "parietal_occipital", "occipital"}),
            ("Subcortical", {"subcortical"}),
            ("Temporal", {"temporal"}),
        ):
            network_to_level[network] = values
        predictions = rows[pred_column].map(
            lambda value: set().union(*(network_to_level.get(item, set()) for item in split(value)))
        )
    truths = rows[truth_column].map(split)
    hits = [bool(prediction & truth) for prediction, truth in zip(predictions, truths)]
    return {
        "level": level,
        "top_k": topk,
        "rule": "tolerant" if tolerant else "strict",
        "n": len(hits),
        "correct": int(sum(hits)),
        "accuracy": float(np.mean(hits)) if hits else None,
    }


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(AUDIT)
    old = pd.read_csv(OLD_TRUTH).set_index("patient_barcode")
    lut = parse_tzo_lut(LUT)
    info = {}
    for label, region in lut.items():
        lobe, broad, network = anatomy_for_region(region)
        info[label] = {"region": region, "lobe": lobe, "broad": broad, "network": network}
    reference = nib.load(str(audit.iloc[0]["t1_path"]))
    atlas = resample_atlas(reference)

    rows = []
    for item in audit.itertuples(index=False):
        seg = np.rint(np.asanyarray(nib.load(item.preferred_segmentation).dataobj)).astype(np.int16)
        tumor = seg > 0
        direct = atlas[tumor]
        label_info = summarize(direct, info)
        old_row = old.loc[item.patient_barcode]
        row = {
            "patient_barcode": item.patient_barcode,
            "segmentation_source": "manual_corrected" if bool(item.has_manual_segmentation) else "automatic",
            "tumor_voxels": int(tumor.sum()),
            "direct_labeled_voxels": int((direct > 0).sum()),
            "direct_overlap_fraction": float(np.mean(direct > 0)),
            "network_prediction_top1": old_row["network_prediction_top1"],
            "network_prediction_top3": old_row["network_prediction_top3"],
            "old_network_dominant": old_row["whole_tumor_network_dominant"],
            "old_lobe_dominant": old_row["whole_tumor_lobe_dominant"],
            "old_broad_dominant": old_row["whole_tumor_broad_dominant"],
        }
        row.update({f"corrected_{key}": value for key, value in label_info.items()})
        row["network_changed_vs_old"] = row["corrected_network_dominant"] != row["old_network_dominant"]
        row["lobe_changed_vs_old"] = row["corrected_lobe_dominant"] != row["old_lobe_dominant"]
        row["broad_changed_vs_old"] = row["corrected_broad_dominant"] != row["old_broad_dominant"]
        row["network_evaluable"] = row["corrected_network_dominant"] != "out_of_scope"
        row["low_direct_coverage_flag"] = row["direct_overlap_fraction"] < 0.50
        rows.append(row)

    frame = pd.DataFrame(rows)
    frame.to_csv(OUTDIR / "corrected_direct_overlap_mri_truth.csv", index=False, encoding="utf-8-sig")
    metrics = pd.DataFrame(
        [
            metric(frame, level, topk, tolerant)
            for level in ("network", "lobe", "broad")
            for topk in (1, 3)
            for tolerant in (False, True)
        ]
    )
    metrics.to_csv(OUTDIR / "corrected_accuracy_metrics.csv", index=False, encoding="utf-8-sig")
    summary = {
        "method": "direct TZO116 overlap only; no unbounded nearest-label filling",
        "n_patients": int(len(frame)),
        "n_network_evaluable": int(frame["network_evaluable"].sum()),
        "n_out_of_scope_cerebellar": int((~frame["network_evaluable"]).sum()),
        "n_low_direct_coverage_below_50pct": int(frame["low_direct_coverage_flag"].sum()),
        "n_truth_changed_vs_old": {
            "network": int(frame["network_changed_vs_old"].sum()),
            "lobe": int(frame["lobe_changed_vs_old"].sum()),
            "broad": int(frame["broad_changed_vs_old"].sum()),
        },
        "metrics": metrics.to_dict("records"),
    }
    (OUTDIR / "corrected_truth_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
