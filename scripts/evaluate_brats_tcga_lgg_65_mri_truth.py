#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "mri_label_tools"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import nibabel as nib  # noqa: E402
from nibabel.processing import resample_from_to  # noqa: E402


AUDIT = ROOT / "results" / "brats_tcga_lgg_training_65_audit_20260609" / "brats_tcga_lgg_training_65_patient_audit.csv"
NETWORK_PRED = ROOT / "results" / "tcga_gbm_lgg_sample_mri_label_tracing_20260605" / "tcga_gbm_lgg_sample_network_tracing.csv"
ATLAS = ROOT / "data" / "atlases" / "sri24" / "labels" / "sri24" / "tzo116plus.nii"
ATLAS_LUT = ROOT / "data" / "atlases" / "sri24" / "labels" / "sri24" / "SRI24-tzo116plus.txt"
OUTDIR = ROOT / "results" / "brats_tcga_lgg_65_mri_truth_evaluation_20260609"


NETWORK_TO_LOBES = {
    "Cingulate gyrus": {"cingulate"},
    "Frontal (agranular frontal motor areas)": {"frontal"},
    "Hippocampal formation": {"temporal"},
    "Lateral Prefrontal Cortex": {"frontal"},
    "Occipital/Temporal": {"occipital", "temporal"},
    "Operculum/Insula": {"insula"},
    "Orbitomedial Prefrontal Cortex (OMPFC)": {"frontal"},
    "Parietal, and Parieto-occipital region": {"parietal", "occipital"},
    "Subcortical": {"subcortical"},
    "Temporal": {"temporal"},
}

NETWORK_TO_BROAD = {
    "Cingulate gyrus": {"cingulate"},
    "Frontal (agranular frontal motor areas)": {"frontal"},
    "Hippocampal formation": {"medial_temporal", "temporal"},
    "Lateral Prefrontal Cortex": {"frontal"},
    "Occipital/Temporal": {"occipital_temporal", "occipital", "temporal"},
    "Operculum/Insula": {"insula_operculum", "insula"},
    "Orbitomedial Prefrontal Cortex (OMPFC)": {"frontal"},
    "Parietal, and Parieto-occipital region": {"parietal_occipital", "parietal", "occipital"},
    "Subcortical": {"subcortical"},
    "Temporal": {"temporal"},
}


def parse_tzo_lut(path: Path) -> dict[int, str]:
    out: dict[int, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.strip().split()
        if len(parts) >= 2 and parts[0].isdigit():
            label = int(parts[0])
            if 1 <= label <= 116:
                out[label] = parts[1]
    return out


def anatomy_for_region(name: str) -> tuple[str, str, str]:
    base = name.removesuffix("_L").removesuffix("_R")
    if base.startswith("Precentral") or base.startswith("Supp_Motor"):
        return "frontal", "frontal", "Frontal (agranular frontal motor areas)"
    if base.startswith("Paracentral"):
        return "parietal", "parietal", "Parietal, and Parieto-occipital region"
    if base.startswith("Frontal_Inf_Oper") or base.startswith("Rolandic_Oper"):
        return "frontal", "insula_operculum", "Operculum/Insula"
    if base.startswith("Frontal"):
        if "_Orb" in base or "_Med" in base:
            return "frontal", "frontal", "Orbitomedial Prefrontal Cortex (OMPFC)"
        return "frontal", "frontal", "Lateral Prefrontal Cortex"
    if base.startswith("Olfactory") or base.startswith("Rectus"):
        return "frontal", "frontal", "Orbitomedial Prefrontal Cortex (OMPFC)"
    if base.startswith("Insula"):
        return "insula", "insula_operculum", "Operculum/Insula"
    if base.startswith("Cingulum"):
        return "cingulate", "cingulate", "Cingulate gyrus"
    if base.startswith(("Hippocampus", "ParaHippocampal")):
        return "temporal", "medial_temporal", "Hippocampal formation"
    if base.startswith("Amygdala"):
        return "subcortical", "subcortical", "Subcortical"
    if base.startswith(("Calcarine", "Cuneus", "Occipital", "Lingual")):
        return "occipital", "occipital", "Occipital/Temporal"
    if base.startswith("Fusiform"):
        return "temporal", "occipital_temporal", "Occipital/Temporal"
    if base.startswith(("Postcentral", "Parietal", "SupraMarginal", "Angular", "Precuneus")):
        return "parietal", "parietal", "Parietal, and Parieto-occipital region"
    if base.startswith(("Heschl", "Temporal")):
        return "temporal", "temporal", "Temporal"
    if base.startswith(("Caudate", "Putamen", "Pallidum", "Thalamus")):
        return "subcortical", "subcortical", "Subcortical"
    if base.startswith(("Cerebelum", "Vermis")):
        return "cerebellum", "cerebellum", "out_of_scope"
    raise ValueError(f"Unmapped TZO region: {name}")


def resample_atlas_to_brats(atlas_path: Path, reference_path: Path) -> np.ndarray:
    atlas_img = nib.load(str(atlas_path))
    atlas_data = np.squeeze(np.asanyarray(atlas_img.dataobj)).astype(np.int16)
    atlas_3d = nib.Nifti1Image(atlas_data, atlas_img.affine, atlas_img.header)
    reference = nib.load(str(reference_path))
    resampled = resample_from_to(atlas_3d, (reference.shape[:3], reference.affine), order=0)
    data = np.rint(np.asanyarray(resampled.dataobj)).astype(np.int16)
    data[(data < 1) | (data > 116)] = 0
    return data


def nearest_filled_labels(atlas: np.ndarray) -> np.ndarray:
    zero = atlas == 0
    indices = distance_transform_edt(zero, return_distances=False, return_indices=True)
    return atlas[tuple(indices)].astype(np.int16)


def summarize_mask(
    mask: np.ndarray,
    direct_atlas: np.ndarray,
    filled_atlas: np.ndarray,
    label_info: dict[int, dict[str, str]],
) -> dict[str, Any]:
    if not np.any(mask):
        return {}
    direct_fraction = float(np.mean(direct_atlas[mask] > 0))
    labels, counts = np.unique(filled_atlas[mask], return_counts=True)
    total = int(counts.sum())
    dimensions: dict[str, dict[str, int]] = {
        "region": defaultdict(int),
        "lobe": defaultdict(int),
        "broad": defaultdict(int),
        "network": defaultdict(int),
    }
    for label, count in zip(labels, counts):
        info = label_info[int(label)]
        dimensions["region"][info["region"]] += int(count)
        dimensions["lobe"][info["lobe"]] += int(count)
        dimensions["broad"][info["broad"]] += int(count)
        dimensions["network"][info["network"]] += int(count)

    out: dict[str, Any] = {"direct_atlas_overlap_fraction": direct_fraction}
    for dimension, values in dimensions.items():
        ordered = sorted(values.items(), key=lambda item: (-item[1], item[0]))
        fractions = [(name, count / total) for name, count in ordered]
        candidates = [name for name, fraction in fractions if fraction >= 0.20]
        if not candidates:
            candidates = [fractions[0][0]]
        out[f"{dimension}_dominant"] = fractions[0][0]
        out[f"{dimension}_dominant_fraction"] = float(fractions[0][1])
        out[f"{dimension}_candidates"] = " | ".join(candidates)
        out[f"{dimension}_distribution"] = " | ".join(f"{name}:{fraction:.4f}" for name, fraction in fractions)
    return out


def split_candidates(value: Any) -> set[str]:
    return {part.strip() for part in str(value or "").split("|") if part.strip()}


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(center - half), float(center + half)


def add_prediction_matches(row: dict[str, Any], predictions: list[str], prefix: str) -> None:
    if f"{prefix}_network_dominant" not in row:
        for level in ("network", "lobe", "broad"):
            for k in (1, 3):
                for rule in ("strict", "tolerant"):
                    row[f"{prefix}_{level}_top{k}_{rule}"] = None
        return
    truth_network_dominant = {str(row[f"{prefix}_network_dominant"])}
    truth_network_candidates = split_candidates(row[f"{prefix}_network_candidates"])
    truth_lobe_dominant = {str(row[f"{prefix}_lobe_dominant"])}
    truth_lobe_candidates = split_candidates(row[f"{prefix}_lobe_candidates"])
    truth_broad_dominant = {str(row[f"{prefix}_broad_dominant"])}
    truth_broad_candidates = split_candidates(row[f"{prefix}_broad_candidates"])
    for k in (1, 3):
        top = predictions[:k]
        pred_lobes = set().union(*(NETWORK_TO_LOBES.get(item, set()) for item in top))
        pred_broad = set().union(*(NETWORK_TO_BROAD.get(item, set()) for item in top))
        row[f"{prefix}_network_top{k}_strict"] = bool(set(top) & truth_network_dominant)
        row[f"{prefix}_network_top{k}_tolerant"] = bool(set(top) & truth_network_candidates)
        row[f"{prefix}_lobe_top{k}_strict"] = bool(pred_lobes & truth_lobe_dominant)
        row[f"{prefix}_lobe_top{k}_tolerant"] = bool(pred_lobes & truth_lobe_candidates)
        row[f"{prefix}_broad_top{k}_strict"] = bool(pred_broad & truth_broad_dominant)
        row[f"{prefix}_broad_top{k}_tolerant"] = bool(pred_broad & truth_broad_candidates)


def metric_table(df: pd.DataFrame, cohort: str, mask: pd.Series) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    sub = df.loc[mask]
    for truth_basis in ("whole_tumor", "core", "edema", "center"):
        if truth_basis == "center" and f"{truth_basis}_network_dominant" not in sub:
            continue
        for level in ("network", "lobe", "broad"):
            for k in (1, 3):
                for rule in ("strict", "tolerant"):
                    col = f"{truth_basis}_{level}_top{k}_{rule}"
                    if col not in sub:
                        continue
                    values = sub[col].dropna().astype(bool)
                    n = int(len(values))
                    correct = int(values.sum())
                    low, high = wilson_interval(correct, n)
                    rows.append(
                        {
                            "cohort": cohort,
                            "truth_basis": truth_basis,
                            "level": level,
                            "top_k": k,
                            "matching_rule": rule,
                            "n": n,
                            "correct": correct,
                            "accuracy": correct / n if n else np.nan,
                            "ci95_low": low,
                            "ci95_high": high,
                        }
                    )
    return rows


def center_summary(
    mask: np.ndarray,
    filled_atlas: np.ndarray,
    label_info: dict[int, dict[str, str]],
) -> dict[str, Any]:
    coords = np.argwhere(mask)
    center = np.rint(coords.mean(axis=0)).astype(int)
    label = int(filled_atlas[tuple(center)])
    info = label_info[label]
    return {
        "center_i": int(center[0]),
        "center_j": int(center[1]),
        "center_k": int(center[2]),
        "region_dominant": info["region"],
        "region_dominant_fraction": 1.0,
        "region_candidates": info["region"],
        "region_distribution": f'{info["region"]}:1.0000',
        "lobe_dominant": info["lobe"],
        "lobe_dominant_fraction": 1.0,
        "lobe_candidates": info["lobe"],
        "lobe_distribution": f'{info["lobe"]}:1.0000',
        "broad_dominant": info["broad"],
        "broad_dominant_fraction": 1.0,
        "broad_candidates": info["broad"],
        "broad_distribution": f'{info["broad"]}:1.0000',
        "network_dominant": info["network"],
        "network_dominant_fraction": 1.0,
        "network_candidates": info["network"],
        "network_distribution": f'{info["network"]}:1.0000',
    }


def make_accuracy_plot(metrics: pd.DataFrame, path: Path) -> None:
    primary = metrics[
        (metrics["cohort"] == "all_65")
        & (metrics["truth_basis"] == "whole_tumor")
        & (metrics["matching_rule"] == "tolerant")
    ].copy()
    primary["label"] = primary["level"].str.title() + " Top" + primary["top_k"].astype(str)
    order = ["Network Top1", "Network Top3", "Lobe Top1", "Lobe Top3", "Broad Top1", "Broad Top3"]
    primary["label"] = pd.Categorical(primary["label"], order, ordered=True)
    primary = primary.sort_values("label")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(primary))
    y = primary["accuracy"].to_numpy()
    err = np.vstack([y - primary["ci95_low"].to_numpy(), primary["ci95_high"].to_numpy() - y])
    bars = ax.bar(x, y, color=["#3B6EA8", "#3B6EA8", "#4F8F68", "#4F8F68", "#9B6A3C", "#9B6A3C"])
    ax.errorbar(x, y, yerr=err, fmt="none", ecolor="#222222", capsize=4, linewidth=1)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Accuracy")
    ax.set_xticks(x, primary["label"], rotation=20, ha="right")
    ax.grid(axis="y", alpha=0.2)
    for bar, value in zip(bars, y):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.1%}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_network_diagnostic_plot(truth: pd.DataFrame, path: Path) -> None:
    labels = list(NETWORK_TO_LOBES)
    true_counts = truth["whole_tumor_network_dominant"].value_counts().reindex(labels, fill_value=0)
    pred_counts = truth["network_prediction_top1"].value_counts().reindex(labels, fill_value=0)
    short = {
        "Cingulate gyrus": "Cingulate",
        "Frontal (agranular frontal motor areas)": "Motor frontal",
        "Hippocampal formation": "Hippocampal",
        "Lateral Prefrontal Cortex": "Lateral PFC",
        "Occipital/Temporal": "Occipital/Temporal",
        "Operculum/Insula": "Operculum/Insula",
        "Orbitomedial Prefrontal Cortex (OMPFC)": "OMPFC",
        "Parietal, and Parieto-occipital region": "Parietal/PO",
        "Subcortical": "Subcortical",
        "Temporal": "Temporal",
    }
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(y + 0.18, true_counts.to_numpy(), height=0.34, label="MRI dominant truth", color="#4F8F68")
    ax.barh(y - 0.18, pred_counts.to_numpy(), height=0.34, label="RNA-seq Top1", color="#B05A5A")
    ax.set_yticks(y, [short[item] for item in labels])
    ax.set_xlabel("Patients")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_confusion_matrix(truth: pd.DataFrame, outdir: Path) -> None:
    matrix = pd.crosstab(
        truth["whole_tumor_network_dominant"],
        truth["network_prediction_top1"],
        margins=True,
    )
    matrix.to_csv(outdir / "network_top1_confusion_matrix.csv", encoding="utf-8-sig")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate SRI24-derived MRI truth and evaluate 65 TCGA-LGG RNA-seq traces.")
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--network-predictions", type=Path, default=NETWORK_PRED)
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    parser.add_argument("--atlas-lut", type=Path, default=ATLAS_LUT)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    audit = pd.read_csv(args.audit)
    pred = pd.read_csv(args.network_predictions).sort_values(["patient_barcode", "rank"])
    pred["patient_barcode"] = pred["patient_barcode"].astype(str).str.upper()
    audit["patient_barcode"] = audit["patient_barcode"].astype(str).str.upper()

    lut = parse_tzo_lut(args.atlas_lut)
    label_info: dict[int, dict[str, str]] = {}
    for label, region in lut.items():
        lobe, broad, network = anatomy_for_region(region)
        label_info[label] = {"region": region, "lobe": lobe, "broad": broad, "network": network}
    missing = sorted(set(range(1, 117)) - set(label_info))
    if missing:
        raise ValueError(f"Missing TZO mappings: {missing}")

    atlas = resample_atlas_to_brats(args.atlas, Path(audit.iloc[0]["t1_path"]))
    filled = nearest_filled_labels(atlas)
    nib.save(
        nib.Nifti1Image(atlas, nib.load(str(audit.iloc[0]["t1_path"])).affine),
        str(args.outdir / "sri24_tzo116_resampled_to_brats.nii.gz"),
    )

    rows: list[dict[str, Any]] = []
    for _, audit_row in audit.iterrows():
        patient = audit_row["patient_barcode"]
        seg = np.rint(np.asanyarray(nib.load(str(audit_row["preferred_segmentation"])).dataobj)).astype(np.int16)
        masks = {
            "whole_tumor": seg > 0,
            "core": np.isin(seg, [1, 4]),
            "edema": seg == 2,
        }
        row: dict[str, Any] = {
            "patient_barcode": patient,
            "segmentation_source": "manual_corrected" if bool(audit_row["has_manual_segmentation"]) else "automatic",
            "segmentation_path": audit_row["preferred_segmentation"],
            "tumor_voxels": int(masks["whole_tumor"].sum()),
            "core_voxels": int(masks["core"].sum()),
            "edema_voxels": int(masks["edema"].sum()),
        }
        for prefix, mask in masks.items():
            info = summarize_mask(mask, atlas, filled, label_info)
            for key, value in info.items():
                row[f"{prefix}_{key}"] = value
        for key, value in center_summary(masks["whole_tumor"], filled, label_info).items():
            row[f"center_{key}"] = value

        p = pred[pred["patient_barcode"].eq(patient)].head(3)
        predictions = p["network_id"].astype(str).tolist()
        scores = p["score"].astype(float).tolist()
        row["sample_id"] = p.iloc[0]["sample_id"] if len(p) else ""
        row["network_prediction_top1"] = predictions[0] if predictions else ""
        row["network_prediction_top3"] = " | ".join(predictions)
        row["network_top1_score"] = scores[0] if scores else np.nan
        row["network_top1_top2_margin"] = scores[0] - scores[1] if len(scores) >= 2 else np.nan
        for prefix in ("whole_tumor", "core", "edema", "center"):
            add_prediction_matches(row, predictions, prefix)
        rows.append(row)

    truth = pd.DataFrame(rows)
    truth.to_csv(args.outdir / "brats_tcga_lgg_65_mri_truth_and_predictions.csv", index=False, encoding="utf-8-sig")

    metric_rows: list[dict[str, Any]] = []
    metric_rows.extend(metric_table(truth, "all_65", pd.Series(True, index=truth.index)))
    metric_rows.extend(metric_table(truth, "manual_only_62", truth["segmentation_source"].eq("manual_corrected")))
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(args.outdir / "brats_tcga_lgg_65_accuracy_metrics.csv", index=False, encoding="utf-8-sig")
    make_accuracy_plot(metrics, args.outdir / "brats_tcga_lgg_65_primary_accuracy.png")
    make_network_diagnostic_plot(truth, args.outdir / "brats_tcga_lgg_65_network_distribution.png")
    write_confusion_matrix(truth, args.outdir)

    primary = metrics[
        (metrics["cohort"] == "all_65")
        & (metrics["truth_basis"] == "whole_tumor")
        & (metrics["matching_rule"] == "tolerant")
    ]
    summary = {
        "n_patients": int(len(truth)),
        "n_manual_corrected": int(truth["segmentation_source"].eq("manual_corrected").sum()),
        "n_automatic": int(truth["segmentation_source"].eq("automatic").sum()),
        "atlas": "SRI24/TZO116+ v2.0",
        "truth_definition": "whole-tumor overlap; candidate labels require >=20% overlap; nearest-label fill for unlabeled white matter",
        "median_direct_atlas_overlap_fraction": float(truth["whole_tumor_direct_atlas_overlap_fraction"].median()),
        "primary_metrics": {
            f'{row.level}_top{int(row.top_k)}': {
                "n": int(row.n),
                "correct": int(row.correct),
                "accuracy": float(row.accuracy),
                "ci95": [float(row.ci95_low), float(row.ci95_high)],
            }
            for row in primary.itertuples()
        },
    }
    (args.outdir / "brats_tcga_lgg_65_evaluation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
