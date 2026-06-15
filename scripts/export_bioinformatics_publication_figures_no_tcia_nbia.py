#!/usr/bin/env python
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import export_bioinformatics_publication_figures as base


base.OUT = base.ROOT / "manuscript" / "figures_publication_no_TCIA_NBIA"
base.DATA_OUT = base.OUT / "source_data"


def figure1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    base.draw_box(ax, (0.03, 0.72), 0.25, 0.19, "Primate reference atlas", "819 macaque brain RNA-seq samples\n110 regions; 10 Networks", base.PALE_BLUE, base.BLUE)
    base.draw_box(ax, (0.375, 0.72), 0.25, 0.19, "Fold-local model", "Top200 discriminative genes\nPearson correlation", "#EFF7F2", base.GREEN)
    base.draw_box(ax, (0.72, 0.72), 0.25, 0.19, "Boundary rescue", "Top3 pairwise correlation\nmargin = 0.002", base.PALE_ORANGE, base.VERMILLION)
    base.arrow(ax, (0.28, 0.815), (0.375, 0.815))
    base.arrow(ax, (0.625, 0.815), (0.72, 0.815))

    base.draw_box(ax, (0.05, 0.40), 0.25, 0.17, "Primary endpoint", "Network Top1 and Top3", "#DDEBF7", base.BLUE)
    base.draw_box(ax, (0.375, 0.40), 0.25, 0.17, "Secondary endpoint", "Within-Network adaptive region group", "#E2F0D9", base.GREEN)
    base.draw_box(ax, (0.70, 0.40), 0.25, 0.17, "Exploratory endpoint", "Top3 Network beam exact region", "#FCE4D6", base.ORANGE)
    base.arrow(ax, (0.845, 0.72), (0.175, 0.57))
    base.arrow(ax, (0.30, 0.485), (0.375, 0.485))
    base.arrow(ax, (0.625, 0.485), (0.70, 0.485))

    validation_boxes = [
        (0.08, "Internal validation", "Strict LOSO and leave-one-monkey-out"),
        (0.39, "Cross-species", "AHBA normal-brain label harmonization"),
        (0.70, "Disease domain", "Ivy GAP and TCGA prediction distributions"),
    ]
    for x, title, body in validation_boxes:
        base.draw_box(ax, (x, 0.085), 0.22, 0.16, title, body, "#F3F4F6", body_width=22, body_size=6.3)
    ax.plot([0.175, 0.175], [0.40, 0.32], color=base.GRAY, lw=1.0)
    ax.plot([0.19, 0.81], [0.32, 0.32], color=base.GRAY, lw=1.0)
    for x in (0.19, 0.50, 0.81):
        base.arrow(ax, (x, 0.32), (x, 0.245))

    ax.text(0.5, 0.965, "Hierarchical brain-origin RNA tracing and validation design", ha="center", va="top", fontsize=10, fontweight="bold")
    ax.text(
        0.5,
        0.018,
        "cfRNA is a prospective application; the current evidence is derived primarily from tissue RNA-seq.",
        ha="center",
        va="bottom",
        fontsize=7,
        color=base.VERMILLION,
        fontweight="bold",
    )
    base.save_figure(fig, "Figure1_study_design_no_TCIA_NBIA")


def figure5_tumor_domain() -> None:
    ivy = pd.read_csv(base.IVY_DIST)
    ivy = ivy[ivy["endpoint"] == "network_top1"].copy()
    tcga = pd.read_csv(base.TCGA_SUMMARY)

    ivy["network_short"] = ivy["value"].map(base.short_network)
    ivy_pivot = ivy.pivot_table(index="structure_acronym", columns="network_short", values="fraction", fill_value=0)
    ivy_order = ["CT-reference-histology", "CTmvp-reference-histology", "CTpan-reference-histology", "IT-reference-histology", "LE-reference-histology"]
    ivy_pivot = ivy_pivot.reindex(ivy_order)
    ivy_pivot.to_csv(base.DATA_OUT / "Figure5_ivy_network_distribution.csv")

    tcga_counts = pd.crosstab(tcga["project_id"], tcga["network_top1"])
    tcga_frac = tcga_counts.div(tcga_counts.sum(axis=1), axis=0)
    tcga_frac.columns = [base.short_network(value) for value in tcga_frac.columns]
    tcga_frac.to_csv(base.DATA_OUT / "Figure5_tcga_network_distribution.csv")

    network_colors = {
        "Subcortical": base.BLUE,
        "Hippocampal": base.ORANGE,
        "OMPFC": base.GREEN,
        "Motor frontal": base.PURPLE,
        "Cingulate": base.VERMILLION,
        "Operculum/Insula": base.SKY,
        "Occipital/Temporal": base.YELLOW,
        "Parietal": "#999999",
    }

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 4.1), gridspec_kw={"width_ratios": [1.45, 0.8], "wspace": 0.35})
    left = np.zeros(len(ivy_pivot))
    for column in ivy_pivot.columns:
        values = ivy_pivot[column].to_numpy(float)
        ax_a.barh(np.arange(len(ivy_pivot)), values, left=left, color=network_colors.get(column, "#BDBDBD"), label=column, height=0.68)
        left += values
    ax_a.set_yticks(np.arange(len(ivy_pivot)), ["Cellular tumor", "Microvascular proliferation", "Pseudopalisading cells", "Infiltrating tumor", "Leading edge"])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 1)
    ax_a.set_xticks(np.linspace(0, 1, 6), [f"{value:.0%}" for value in np.linspace(0, 1, 6)])
    ax_a.set_xlabel("Fraction of samples")
    ax_a.set_title("Ivy GAP prediction distribution\n(not anatomical accuracy)")
    base.clean_axis(ax_a, grid=False)
    ax_a.grid(axis="x", color=base.LIGHT_GRAY, linewidth=0.6)
    base.panel_label(ax_a, "A")

    bottom = np.zeros(len(tcga_frac))
    for column in tcga_frac.columns:
        values = tcga_frac[column].to_numpy(float)
        ax_b.bar(np.arange(len(tcga_frac)), values, bottom=bottom, color=network_colors.get(column, "#BDBDBD"), width=0.62)
        bottom += values
    ax_b.set_xticks(np.arange(len(tcga_frac)), tcga_frac.index)
    ax_b.set_ylim(0, 1)
    ax_b.set_yticks(np.linspace(0, 1, 6), [f"{value:.0%}" for value in np.linspace(0, 1, 6)])
    ax_b.set_ylabel("Fraction of samples")
    ax_b.set_title("TCGA prediction distribution\n(no anatomical labels)", pad=8)
    base.clean_axis(ax_b)
    ax_b.text(-0.16, 1.08, "B", transform=ax_b.transAxes, fontsize=11, fontweight="bold", va="top")

    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncols=4, loc="lower center", bbox_to_anchor=(0.5, 0.015), columnspacing=1.2, handlelength=1.4)
    fig.suptitle("Glioma RNA-seq reveals disease-domain prediction shift", y=0.99, fontweight="bold")
    fig.subplots_adjust(bottom=0.24, top=0.80)
    base.save_figure(fig, "Figure5_glioma_domain_shift_no_TCIA_NBIA")


def write_captions() -> None:
    text = """# Publication figure plan and legends

## Main figures

### Figure 1. Hierarchical brain-origin RNA tracing and validation design
The macaque reference atlas is used to build fold-local discriminative-correlation models. Network is the prespecified primary endpoint; within-Network region group is secondary and exact region is exploratory. Internal validation uses strict leave-one-sample-out (LOSO) and leave-one-monkey-out (LOMO) designs. AHBA provides label-harmonized cross-species evidence, whereas Ivy GAP and TCGA characterize disease-domain behavior. cfRNA is a prospective application; the current evidence is derived primarily from tissue RNA-seq.

### Figure 2. Internal validation supports Network as the primary endpoint
(A) Network Top1 accuracy for the Top200 discriminative-correlation baseline and Top3-constrained pairwise rescue. The LOSO paired comparison used the same 819 held-out samples (two-sided paired P=0.0026). (B) Outcome of the 99 LOSO rescue switches. Forty-five corrected baseline errors, 20 introduced errors and 34 produced no net correctness change, yielding a net gain of 25 correct samples. (C) Top1 and Top3 accuracy across the endpoint hierarchy. Network is primary, Region Group secondary and Exact Region exploratory. Exact evaluable sample counts were 814 for LOSO and 812 for LOMO.

### Figure 3. Threshold selection and cross-individual robustness
(A) Post-hoc screen of the minimum pairwise-rescue margin using existing strict LOSO predictions without model refitting. The selected threshold of 0.002 maximized Top1 accuracy while preserving Top3. The red line denotes the number of switched samples. This panel documents retrospective threshold selection and is not an independent validation. (B) Network Top1 accuracy for each held-out monkey in LOMO validation. Point size represents held-out sample count and horizontal intervals are 95% Wilson confidence intervals. The dashed line denotes aggregate LOMO accuracy (53.2%).

### Figure 4. Cross-species transfer to normal human brain RNA-seq
(A) Label-harmonized AHBA accuracy. Network and broad-anatomy metrics were evaluated in 233 supported samples; Exact Region metrics were restricted to 91 exact-mapped samples. (B) Row-normalized broad-anatomy prediction pattern. AHBA labels and macaque-derived predictions are harmonized coarse categories and are not directly equivalent to the macaque internal-validation labels.

### Figure 5. Glioma RNA-seq reveals disease-domain prediction shift
(A) Ivy GAP Network Top1 prediction distributions across five glioblastoma microanatomical structures. Ivy labels describe tumor compartments rather than normal brain location; the panel is not an accuracy analysis. (B) TCGA-GBM/LGG Network Top1 prediction distributions. TCGA has no anatomical truth in this analysis, so the panel characterizes domain shift rather than localization accuracy.

## Supplementary figures

### Figure S1. LOMO performance across endpoint resolutions
Heatmap of Top1 and Top3 accuracy for Network, Region Group and Exact Region in each held-out monkey. Sample counts are shown below monkey identifiers.

### Figure S2. Pairwise rescue switch audit across margin thresholds
Top1 gains, Top1 losses and total switched samples across the retrospective margin screen. The vertical dashed line denotes the selected threshold of 0.002.

## Style specification

- Final width: 180 mm for multi-panel main figures.
- Typeface: Arial; minimum plotted text approximately 7 pt at final size.
- Color palette: Okabe-Ito color-blind-safe palette.
- Outputs: 600 dpi PNG plus vector PDF and SVG.
- All figure source tables are stored in `source_data`.
"""
    (base.OUT / "figure_legends_and_style.md").write_text(text, encoding="utf-8")


def main() -> int:
    base.OUT.mkdir(parents=True, exist_ok=True)
    base.DATA_OUT.mkdir(parents=True, exist_ok=True)
    base.setup_style()
    figure1_workflow()
    base.figure2_primary_results()
    base.figure3_threshold_and_lomo()
    base.figure4_ahba()
    figure5_tumor_domain()
    base.supplementary_figures()
    write_captions()
    print(f"Publication figures without TCIA/NBIA written to {base.OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
