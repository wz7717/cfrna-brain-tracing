#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PREDICTIONS = (
    ROOT
    / "results"
    / "tcga_lgg_65_tumor_adapted_network_20260612"
    / "tcga_lgg_65_tumor_adapted_predictions.csv"
)
DEFAULT_TRUTH = (
    ROOT
    / "results"
    / "brats_tcga_lgg_65_mri_truth_corrected_20260612"
    / "corrected_direct_overlap_mri_truth.csv"
)
DEFAULT_OUTDIR = ROOT / "results" / "tcga_lgg_65_broad_anatomy_output_20260613"

# The primary mapping is used for a concise single-label report. Compatible
# labels retain the anatomical breadth of networks that span adjacent areas.
NETWORK_TO_PRIMARY_BROAD = {
    "Cingulate gyrus": "cingulate",
    "Frontal (agranular frontal motor areas)": "frontal",
    "Hippocampal formation": "medial_temporal",
    "Lateral Prefrontal Cortex": "frontal",
    "Occipital/Temporal": "occipital_temporal",
    "Operculum/Insula": "insula_operculum",
    "Orbitomedial Prefrontal Cortex (OMPFC)": "frontal",
    "Parietal, and Parieto-occipital region": "parietal_occipital",
    "Subcortical": "subcortical",
    "Temporal": "temporal",
}

NETWORK_TO_COMPATIBLE_BROAD = {
    "Cingulate gyrus": ("cingulate",),
    "Frontal (agranular frontal motor areas)": ("frontal",),
    "Hippocampal formation": ("medial_temporal", "temporal"),
    "Lateral Prefrontal Cortex": ("frontal",),
    "Occipital/Temporal": ("occipital_temporal", "occipital", "temporal"),
    "Operculum/Insula": ("insula_operculum", "insula"),
    "Orbitomedial Prefrontal Cortex (OMPFC)": ("frontal",),
    "Parietal, and Parieto-occipital region": (
        "parietal_occipital",
        "parietal",
        "occipital",
    ),
    "Subcortical": ("subcortical",),
    "Temporal": ("temporal",),
}


def split_labels(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split("|") if item.strip()]


def ordered_unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def primary_broad(value: Any) -> str:
    networks = split_labels(value)
    return NETWORK_TO_PRIMARY_BROAD.get(networks[0], "unmapped") if networks else ""


def compatible_broad(value: Any) -> str:
    labels: list[str] = []
    for network in split_labels(value):
        labels.extend(NETWORK_TO_COMPATIBLE_BROAD.get(network, ("unmapped",)))
    return " | ".join(ordered_unique(labels))


def add_broad_columns(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    routes: list[str] = []
    for column in frame.columns:
        if not column.endswith("_top1"):
            continue
        route = column.removesuffix("_top1")
        top3 = f"{route}_top3"
        if top3 not in frame.columns:
            continue
        routes.append(route)
        frame[f"{route}_broad_top1_primary"] = frame[column].map(primary_broad)
        frame[f"{route}_broad_top1_compatible"] = frame[column].map(compatible_broad)
        frame[f"{route}_broad_top3"] = frame[top3].map(compatible_broad)
    return frame, routes


def evaluate(frame: pd.DataFrame, routes: list[str]) -> pd.DataFrame:
    if "corrected_broad_dominant" not in frame or "corrected_broad_candidates" not in frame:
        return pd.DataFrame()
    rows = []
    dominant = frame["corrected_broad_dominant"].astype(str)
    candidates = frame["corrected_broad_candidates"].map(set_from_pipe)
    for route in routes:
        top1_primary = frame[f"{route}_broad_top1_primary"].astype(str)
        top1_compatible = frame[f"{route}_broad_top1_compatible"].map(set_from_pipe)
        top3 = frame[f"{route}_broad_top3"].map(set_from_pipe)
        metrics = {
            "top1_primary_strict": top1_primary.eq(dominant),
            "top1_compatible_strict": [
                truth in predicted for truth, predicted in zip(dominant, top1_compatible)
            ],
            "top1_compatible_tolerant": [
                bool(truth & predicted) for truth, predicted in zip(candidates, top1_compatible)
            ],
            "top3_strict": [truth in predicted for truth, predicted in zip(dominant, top3)],
            "top3_tolerant": [
                bool(truth & predicted) for truth, predicted in zip(candidates, top3)
            ],
        }
        for metric, hits in metrics.items():
            hits = pd.Series(hits, dtype=bool)
            rows.append(
                {
                    "route": route,
                    "metric": metric,
                    "n": len(hits),
                    "correct": int(hits.sum()),
                    "accuracy": float(hits.mean()),
                }
            )
    return pd.DataFrame(rows)


def set_from_pipe(value: Any) -> set[str]:
    return set(split_labels(value))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--truth", type=Path, default=DEFAULT_TRUTH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    output = pd.read_csv(args.predictions)
    output, routes = add_broad_columns(output)

    if args.truth.exists():
        truth = pd.read_csv(args.truth)
        truth_columns = [
            "patient_barcode",
            "corrected_broad_dominant",
            "corrected_broad_candidates",
            "corrected_broad_distribution",
            "direct_overlap_fraction",
            "low_direct_coverage_flag",
        ]
        output = output.drop(
            columns=[column for column in truth_columns[1:] if column in output],
            errors="ignore",
        ).merge(truth[truth_columns], on="patient_barcode", how="left", validate="one_to_one")

    args.outdir.mkdir(parents=True, exist_ok=True)
    output_path = args.outdir / "human_sample_broad_anatomy_predictions.csv"
    metrics_path = args.outdir / "broad_anatomy_accuracy_metrics.csv"
    output.to_csv(output_path, index=False, encoding="utf-8-sig")
    metrics = evaluate(output, routes)
    metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")
    print(f"patients={len(output)} routes={len(routes)}")
    print(output_path)
    print(metrics_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
