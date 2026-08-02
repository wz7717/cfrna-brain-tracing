#!/usr/bin/env python
"""Audit inflation of AHBA mapped-label Network Top3 under multi-label truth."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import pandas as pd


ROUTE = "hybrid_projected_network_logcpm_exact"
NETWORKS = 10
K = 3


def split_pipe(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [item.strip() for item in str(value).split("|") if item.strip()]


def hit(predicted: list[str], allowed: list[str]) -> int:
    return int(bool(set(predicted) & set(allowed)))


def precision_at_k(predicted: list[str], allowed: list[str], k: int = K) -> float:
    return len(set(predicted[:k]) & set(allowed)) / k


def recall_at_k(predicted: list[str], allowed: list[str], k: int = K) -> float:
    return len(set(predicted[:k]) & set(allowed)) / len(set(allowed))


def random_hit_probability(m: int, n: int = NETWORKS, k: int = K) -> float:
    if m <= 0:
        return 0.0
    if n - m < k:
        return 1.0
    return 1.0 - math.comb(n - m, k) / math.comb(n, k)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--detail", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.detail)
    frame = frame[
        frame["route"].eq(ROUTE) & frame["supported_for_accuracy"].astype(bool)
    ].copy()
    frame["allowed"] = frame["allowed_bo2023_networks"].map(split_pipe)
    frame["predicted"] = frame.apply(
        lambda row: [str(row.network_top1), str(row.network_top2), str(row.network_top3)],
        axis=1,
    )
    frame["n_allowed"] = frame["allowed"].map(len)
    frame["unique_label"] = frame["n_allowed"].eq(1)
    frame["locked_top1_hit"] = frame.apply(
        lambda row: hit(row.predicted[:1], row.allowed), axis=1
    )
    frame["locked_top3_any_hit"] = frame.apply(
        lambda row: hit(row.predicted[:3], row.allowed), axis=1
    )
    frame["locked_precision_at3"] = frame.apply(
        lambda row: precision_at_k(row.predicted, row.allowed), axis=1
    )
    frame["locked_recall_at3"] = frame.apply(
        lambda row: recall_at_k(row.predicted, row.allowed), axis=1
    )
    frame["random_top3_any_probability"] = frame["n_allowed"].map(random_hit_probability)
    frame["random_expected_precision_at3"] = frame["n_allowed"] / NETWORKS

    # Leave-one-out frequency baseline: for each sample, rank Networks by the
    # number of other samples for which each Network is an allowed label.
    universe = sorted({label for labels in frame["allowed"] for label in labels})
    frequency_predictions: list[list[str]] = []
    for idx, row in frame.iterrows():
        other = frame.drop(index=idx)
        counts = {
            label: sum(label in labels for labels in other["allowed"])
            for label in universe
        }
        frequency_predictions.append(
            sorted(universe, key=lambda label: (-counts[label], label))[:K]
        )
    frame["frequency_predicted"] = frequency_predictions
    frame["frequency_top1_hit"] = frame.apply(
        lambda row: hit(row.frequency_predicted[:1], row.allowed), axis=1
    )
    frame["frequency_top3_any_hit"] = frame.apply(
        lambda row: hit(row.frequency_predicted[:3], row.allowed), axis=1
    )
    frame["frequency_precision_at3"] = frame.apply(
        lambda row: precision_at_k(row.frequency_predicted, row.allowed), axis=1
    )

    export = frame.copy()
    for column in ("allowed", "predicted", "frequency_predicted"):
        export[column] = export[column].map(lambda values: " | ".join(values))
    export[
        [
            "sample_id", "ahba_donor", "public_label", "allowed",
            "n_allowed", "unique_label", "predicted", "locked_top1_hit",
            "locked_top3_any_hit", "locked_precision_at3", "locked_recall_at3",
            "random_top3_any_probability", "random_expected_precision_at3",
            "frequency_predicted", "frequency_top1_hit",
            "frequency_top3_any_hit", "frequency_precision_at3",
        ]
    ].to_csv(args.outdir / "ahba_network_multilabel_audit_per_sample.csv", index=False)

    rows = []
    for subset_name, subset in (
        ("all_supported", frame),
        ("unique_label", frame[frame["unique_label"]]),
        ("multi_label", frame[~frame["unique_label"]]),
    ):
        rows.append(
            {
                "subset": subset_name,
                "n": int(len(subset)),
                "mean_allowed_labels": float(subset["n_allowed"].mean()),
                "locked_top1_accuracy": float(subset["locked_top1_hit"].mean()),
                "locked_top3_any_hit": float(subset["locked_top3_any_hit"].mean()),
                "locked_precision_at3": float(subset["locked_precision_at3"].mean()),
                "locked_recall_at3": float(subset["locked_recall_at3"].mean()),
                "random_expected_top3_any_hit": float(
                    subset["random_top3_any_probability"].mean()
                ),
                "random_expected_precision_at3": float(
                    subset["random_expected_precision_at3"].mean()
                ),
                "frequency_top1_accuracy": float(subset["frequency_top1_hit"].mean()),
                "frequency_top3_any_hit": float(
                    subset["frequency_top3_any_hit"].mean()
                ),
                "frequency_precision_at3": float(
                    subset["frequency_precision_at3"].mean()
                ),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(args.outdir / "ahba_network_multilabel_audit_summary.csv", index=False)

    payload = {
        "route": ROUTE,
        "n_supported": int(len(frame)),
        "n_unique_label": int(frame["unique_label"].sum()),
        "n_multi_label": int((~frame["unique_label"]).sum()),
        "multi_label_fraction": float((~frame["unique_label"]).mean()),
        "metric_boundary": (
            "The historical Top3 value is an any-allowed-label hit rate under "
            "set-valued truth, not conventional single-label classification accuracy."
        ),
        "random_baseline": (
            "Analytic probability for three uniformly sampled distinct Networks "
            "to overlap a sample's allowed-label set."
        ),
        "frequency_baseline": (
            "Leave-one-sample-out ranking of Networks by allowed-label frequency; "
            "ties broken alphabetically."
        ),
        "input_sha256": sha256(args.detail),
        "summary": summary.to_dict(orient="records"),
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
