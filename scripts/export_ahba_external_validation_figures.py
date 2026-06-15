#!/usr/bin/env python
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "ahba_human_rnaseq_external_validation_20260603"


def main() -> int:
    metrics = json.loads((OUTDIR / "ahba_rnaseq_external_validation_metrics.json").read_text(encoding="utf-8"))
    by_label = pd.read_csv(OUTDIR / "ahba_rnaseq_external_validation_by_label.csv")
    sample = pd.read_csv(OUTDIR / "ahba_rnaseq_external_validation_sample_summary.csv")

    metric_rows = [
        ("Network Top1", metrics["network_top1_accuracy_coarse"]),
        ("Network Top3", metrics["network_top3_accuracy_coarse"]),
        ("Region lobe Top1", metrics["region_lobe_top1_accuracy_coarse"]),
        ("Exact-mapped Region Top1", metrics["region_top1_exact_accuracy_on_exact_mapped_labels"]),
        ("Exact-mapped Region Top3", metrics["region_top3_exact_accuracy_on_exact_mapped_labels"]),
    ]
    fig, ax = plt.subplots(figsize=(10, 5))
    labels = [x[0] for x in metric_rows]
    values = [x[1] for x in metric_rows]
    bars = ax.bar(labels, values, color=["#4C78A8", "#4C78A8", "#59A14F", "#F28E2B", "#F28E2B"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Accuracy")
    ax.set_title("AHBA human RNA-seq external validation")
    ax.tick_params(axis="x", rotation=25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.02, f"{value:.1%}", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(OUTDIR / "ahba_external_validation_accuracy_summary.png", dpi=180)
    plt.close(fig)

    supported = by_label[by_label["supported_for_accuracy"].astype(bool)].copy()
    supported = supported.sort_values("region_lobe_top1_accuracy", ascending=False)
    fig, ax = plt.subplots(figsize=(12, 7))
    plot_labels = supported["public_label"].astype(str)
    ax.barh(plot_labels, supported["region_lobe_top1_accuracy"], color="#59A14F")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Region lobe Top1 accuracy")
    ax.set_title("AHBA label-harmonized lobe accuracy by label")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUTDIR / "ahba_external_validation_by_label_lobe_accuracy.png", dpi=180)
    plt.close(fig)

    confusion = pd.crosstab(sample["public_major_anatomy"], sample["predicted_public_major_top1"], normalize="index")
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(confusion.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(confusion.shape[1]))
    ax.set_xticklabels(confusion.columns, rotation=45, ha="right")
    ax.set_yticks(range(confusion.shape[0]))
    ax.set_yticklabels(confusion.index)
    ax.set_title("AHBA true public anatomy vs predicted public anatomy")
    fig.colorbar(im, ax=ax, label="Row fraction")
    fig.tight_layout()
    fig.savefig(OUTDIR / "ahba_external_validation_public_anatomy_confusion.png", dpi=180)
    plt.close(fig)
    print(f"Figures written to {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
