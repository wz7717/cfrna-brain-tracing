#!/usr/bin/env python
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import export_bioinformatics_publication_figures as base


base.OUT = base.ROOT / "manuscript" / "figures_publication_20260613"
base.DATA_OUT = base.OUT / "source_data"


def figure1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.35))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    base.draw_box(ax, (0.03, 0.74), 0.25, 0.17, "Primate reference atlas",
                  "819 macaque RNA-seq samples\n110 regions; 10 Networks",
                  base.PALE_BLUE, base.BLUE)
    base.draw_box(ax, (0.375, 0.74), 0.25, 0.17, "Fold-local model",
                  "Top-200 discriminative genes\nPearson correlation",
                  "#EFF7F2", base.GREEN)
    base.draw_box(ax, (0.72, 0.74), 0.25, 0.17, "Boundary rescue",
                  "Top-3 pairwise correlation\nmargin = 0.002",
                  base.PALE_ORANGE, base.VERMILLION)
    base.arrow(ax, (0.28, 0.825), (0.375, 0.825))
    base.arrow(ax, (0.625, 0.825), (0.72, 0.825))

    endpoints = [
        (0.05, "Primary endpoint", "Network Top-1 and Top-3", "#DDEBF7", base.BLUE),
        (0.375, "Secondary endpoint", "Adaptive Region Group", "#E2F0D9", base.GREEN),
        (0.70, "Exploratory endpoint", "Exact Region in Top-3 beam", "#FCE4D6", base.ORANGE),
    ]
    for x, title, body, fill, edge in endpoints:
        base.draw_box(ax, (x, 0.45), 0.25, 0.16, title, body, fill, edge)
    base.arrow(ax, (0.845, 0.74), (0.175, 0.61))
    base.arrow(ax, (0.30, 0.53), (0.375, 0.53))
    base.arrow(ax, (0.625, 0.53), (0.70, 0.53))

    validations = [
        (0.015, "Internal", "LOSO and leave-one-monkey-out"),
        (0.265, "Normal human", "AHBA harmonized labels"),
        (0.515, "Paired glioma", "TCGA-LGG RNA-seq + BraTS MRI"),
        (0.765, "Liquid biopsy", "EV-RNA and CSF-RNA transfer stress tests"),
    ]
    for x, title, body in validations:
        base.draw_box(ax, (x, 0.12), 0.22, 0.17, title, body, "#F3F4F6",
                      body_width=22, body_size=6.1)
    ax.plot([0.175, 0.175], [0.45, 0.37], color=base.GRAY, lw=1.0)
    ax.plot([0.125, 0.875], [0.37, 0.37], color=base.GRAY, lw=1.0)
    for x in (0.125, 0.375, 0.625, 0.875):
        base.arrow(ax, (x, 0.37), (x, 0.29))

    ax.text(0.5, 0.965, "Hierarchical RNA source tracing and validation design",
            ha="center", va="top", fontsize=10, fontweight="bold")
    ax.text(0.5, 0.03,
            "Liquid-biopsy cohorts test transfer and confidence; only paired MRI supplies anatomical truth.",
            ha="center", va="bottom", fontsize=7, color=base.VERMILLION,
            fontweight="bold")
    base.save_figure(fig, "Figure1_study_design_20260613")


def figure4_ahba() -> None:
    coarse = pd.DataFrame(
        [
            ["Network Top-1", 76, 233, 0.3262],
            ["Network Top-3", 129, 233, 0.5536],
            ["Broad anatomy Top-1", 103, 233, 0.4421],
        ],
        columns=["metric", "correct", "n", "accuracy"],
    )
    exact = pd.DataFrame(
        [
            ["Exact Region Top-1", 9, 91, 0.0989],
            ["Exact Region Top-3", 27, 91, 0.2967],
            ["Low-resolution/review", 223, 242, 0.9215],
        ],
        columns=["metric", "correct", "n", "fraction"],
    )
    coarse.to_csv(base.DATA_OUT / "Figure4_AHBA_coarse_metrics.csv", index=False)
    exact.to_csv(base.DATA_OUT / "Figure4_AHBA_resolution_metrics.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25))
    colors = [base.BLUE, base.SKY, base.GREEN]
    axes[0].barh(np.arange(3), coarse["accuracy"], color=colors)
    axes[0].set_yticks(np.arange(3), coarse["metric"])
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 0.65)
    axes[0].set_xlabel("Accuracy")
    axes[0].set_title("A  Harmonized coarse-label validation", loc="left", fontweight="bold")
    for i, row in coarse.iterrows():
        axes[0].text(row["accuracy"] + 0.015, i, f'{row["accuracy"]:.1%}',
                     va="center", fontsize=7)

    axes[1].barh(np.arange(3), exact["fraction"],
                 color=[base.ORANGE, base.VERMILLION, base.GRAY])
    axes[1].set_yticks(np.arange(3), exact["metric"])
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 1.0)
    axes[1].set_xlabel("Fraction")
    axes[1].set_title("B  Resolution and review burden", loc="left", fontweight="bold")
    for i, row in exact.iterrows():
        if row["fraction"] > 0.85:
            axes[1].text(row["fraction"] - 0.02, i, f'{row["fraction"]:.1%}',
                         ha="right", va="center", fontsize=7, color="white",
                         fontweight="bold")
        else:
            axes[1].text(row["fraction"] + 0.02, i, f'{row["fraction"]:.1%}',
                         va="center", fontsize=7)
    fig.suptitle("AHBA cross-species validation", fontsize=10, fontweight="bold")
    fig.tight_layout()
    base.save_figure(fig, "Figure4_AHBA_external_validation_20260613")


def figure5_paired_glioma() -> None:
    domain = pd.DataFrame(
        [
            ["Cramer's V", 0.1580],
            ["Jensen-Shannon divergence", 0.026619],
        ],
        columns=["metric", "value"],
    )
    low_conf = pd.DataFrame(
        [["GBM", 0.384], ["LGG", 0.260]],
        columns=["group", "fraction"],
    )
    paired = pd.DataFrame(
        [
            ["Lobe Top-3", 0.8462, 0.8923],
            ["Broad anatomy Top-3", 0.7538, 0.8308],
            ["Network Top-3", 0.2188, 0.3594],
        ],
        columns=["endpoint", "strict", "tolerant"],
    )
    liquid = pd.DataFrame(
        [
            ["GSE228512\nserum EV", 0.659, 0.000, np.nan],
            ["GSE106804\ntumor EV", 0.231, 0.000, 0.211],
            ["GSE189919\nCSF", 0.650, 0.000, 0.875],
        ],
        columns=["cohort", "dominant_fraction_case_group", "ood_acceptance",
                 "cross_input_or_route_agreement"],
    )
    domain.to_csv(base.DATA_OUT / "Figure5_domain_shift_metrics.csv", index=False)
    low_conf.to_csv(base.DATA_OUT / "Figure5_low_confidence_fraction.csv", index=False)
    paired.to_csv(base.DATA_OUT / "Figure5_paired_baseline_metrics.csv", index=False)
    liquid.to_csv(base.DATA_OUT / "Figure5_liquid_biopsy_diagnostics.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.45))

    axes[0].bar(low_conf["group"], low_conf["fraction"],
                color=[base.VERMILLION, base.BLUE])
    axes[0].set_ylim(0, 0.5)
    axes[0].set_ylabel("Low-confidence fraction")
    axes[0].set_title("A  GBM-LGG domain shift", loc="left", fontweight="bold")
    for i, value in enumerate(low_conf["fraction"]):
        axes[0].text(i, value + 0.018, f"{value:.1%}", ha="center", fontsize=7)
    axes[0].text(0.5, 0.47, "Permutation P = 0.00060\nCramer's V = 0.158",
                 ha="center", va="top", fontsize=6.5)

    y = np.arange(len(paired))
    h = 0.34
    axes[1].barh(y - h / 2, paired["strict"], height=h,
                 color=base.BLUE, label="Strict")
    axes[1].barh(y + h / 2, paired["tolerant"], height=h,
                 color=base.ORANGE, label="Tolerant")
    axes[1].set_yticks(y, paired["endpoint"])
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Accuracy")
    axes[1].set_title("B  Paired MRI baseline", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, fontsize=6.5, loc="lower right")

    x = np.arange(len(liquid))
    axes[2].bar(x, liquid["dominant_fraction_case_group"],
                color=[base.BLUE, base.GREEN, base.ORANGE],
                label="Dominant Top-1 fraction")
    for i, value in enumerate(liquid["cross_input_or_route_agreement"]):
        if np.isfinite(value):
            axes[2].plot(i, value, marker="D", color=base.VERMILLION,
                         markersize=5)
    axes[2].plot([], [], marker="D", linestyle="none", color=base.VERMILLION,
                 label="Route/input agreement")
    axes[2].set_xticks(x, liquid["cohort"], rotation=20, ha="right")
    axes[2].set_ylim(0, 1.0)
    axes[2].set_ylabel("Fraction")
    axes[2].set_title("C  Liquid-biopsy transfer", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, fontsize=6.0, loc="upper center")
    axes[2].text(0.5, 0.72, "Separate adapted-model OOD accepted: 0/186",
                 transform=axes[2].transAxes, ha="center", va="center",
                 fontsize=6.3, fontweight="bold",
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85,
                       "pad": 1.5})

    fig.suptitle("Glioma validation and liquid-biopsy transfer limits",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    base.save_figure(fig, "Figure5_paired_glioma_validation_20260613")


def supplementary_adaptation() -> None:
    sensitivity = pd.DataFrame(
        [
            ["Baseline", 0.0469, 0.3594, 0.8308],
            ["Raw TPM", 0.0781, 0.4375, 0.7692],
            ["log1p", 0.1406, 0.4375, 0.6462],
            ["Harmonized", 0.2031, 0.4844, 0.6000],
            ["Harmonized + calibrated", 0.2031, 0.4688, 0.7077],
        ],
        columns=["route", "network_top1_tolerant", "network_top3_tolerant",
                 "broad_top3_tolerant"],
    )
    sensitivity.to_csv(base.DATA_OUT / "FigureS3_adaptation_sensitivity.csv",
                       index=False)
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    x = np.arange(len(sensitivity))
    ax.plot(x, sensitivity["network_top1_tolerant"], marker="o",
            color=base.VERMILLION, label="Network Top-1")
    ax.plot(x, sensitivity["network_top3_tolerant"], marker="s",
            color=base.BLUE, label="Network Top-3")
    ax.plot(x, sensitivity["broad_top3_tolerant"], marker="^",
            color=base.GREEN, label="Broad anatomy Top-3")
    ax.set_xticks(x, ["Baseline", "Raw TPM", "log1p", "Harmonized",
                      "Harmonized + calibrated"], rotation=20, ha="right")
    ax.set_ylim(0, 0.9)
    ax.set_ylabel("Tolerant accuracy")
    ax.set_title("Exploratory adaptation redistributed performance across endpoints",
                 fontweight="bold")
    ax.legend(frameon=False, ncol=3, loc="upper center")
    ax.text(0.01, 0.03,
            "Harmonized: transductive; calibrated: cohort-internal; "
            "quantile mapping requires tie-aware retraining.",
            transform=ax.transAxes, fontsize=7, color=base.VERMILLION)
    fig.tight_layout()
    base.save_figure(fig, "FigureS3_adaptation_sensitivity_20260613")


def write_captions() -> None:
    text = """# Figure legends and style

## Figure 1. Hierarchical RNA source-tracing framework and validation design
Network is the sole overall primary endpoint. Region Group and Exact Region are secondary and exploratory. Validation proceeds from internal macaque resampling to AHBA normal-human transfer, paired TCGA-LGG/BraTS evaluation, and disease-domain stress tests. cfRNA remains the intended downstream application.

## Figure 2. Internal performance across endpoint resolutions
LOSO and leave-one-monkey-out performance for Network, Region Group and Exact Region. Network provides the most stable level of anatomical inference.

## Figure 3. Boundary-rescue threshold audit and donor-level generalization
Retrospective threshold audit for pairwise rescue and donor-isolated Network performance.

## Figure 4. AHBA cross-species validation
(A) Harmonized Network and broad-anatomy accuracy in 233 supported samples. (B) Exact-region accuracy in 91 evaluable samples and the fraction requiring low-resolution interpretation or manual review.

## Figure 5. Glioma validation and liquid-biopsy transfer limits
(A) Low-confidence prediction fractions in TCGA GBM and LGG; the global distribution difference was assessed by permutation. (B) Locked strict and tolerant Top-3 accuracy in 65 paired TCGA-LGG/BraTS cases, with Network evaluated in 64 in-scope cases. (C) Dominant Top-1 fractions in the case groups of three liquid-biopsy cohorts. Diamonds denote adapted-route agreement in GSE106804 and TPM-versus-CPM locked-baseline agreement in GSE189919. In GSE189919, the case-control Top-1 distribution was not significantly different (100,000-label permutation P=0.178). A separately tested macaque-derived adapted-model OOD rule accepted none of 186 liquid-biopsy samples; this rule was not part of the locked baseline. These cohorts lacked imaging truth and were not used to calculate localization accuracy.

## Supplementary figures
Figure S1 reports leave-one-monkey-out performance by donor and endpoint resolution. Figure S2 audits pairwise-rescue switches over the tested margin thresholds. Figure S3 reports exploratory adaptation sensitivity; harmonized routes are transductive, the calibrated route is cohort-internal, and the current quantile mapping requires tie-aware retraining.

## Style specification
- Final width: 180 mm for multi-panel main figures.
- Typeface: Arial; minimum plotted text approximately 7 pt at final size.
- Palette: Okabe-Ito color-blind-safe colors.
- Outputs: 600 dpi PNG plus vector PDF and SVG.
- Figure source tables are stored in `source_data`.
"""
    (base.OUT / "figure_legends_and_style.md").write_text(text, encoding="utf-8")


def main() -> int:
    base.OUT.mkdir(parents=True, exist_ok=True)
    base.DATA_OUT.mkdir(parents=True, exist_ok=True)
    base.setup_style()
    figure1_workflow()
    base.figure2_primary_results()
    base.figure3_threshold_and_lomo()
    figure4_ahba()
    figure5_paired_glioma()
    base.supplementary_figures()
    supplementary_adaptation()
    write_captions()
    print(f"Publication figures written to {base.OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
