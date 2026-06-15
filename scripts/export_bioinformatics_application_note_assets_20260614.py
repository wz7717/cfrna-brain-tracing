from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import export_bioinformatics_publication_figures as base


ROOT = Path(__file__).resolve().parents[1]
FIG_OUT = ROOT / "manuscript" / "figures_application_note_20260614"
TABLE_OUT = ROOT / "manuscript" / "tables_application_note_20260614"


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG_OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG_OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def export_main_figure() -> None:
    fig = plt.figure(figsize=(7.1, 6.8))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], wspace=0.34, hspace=0.42)

    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    base.draw_box(ax, (0.03, 0.72), 0.27, 0.17, "Input",
                  "Gene symbol + abundance\nOptional sample metadata",
                  base.PALE_BLUE, base.BLUE, body_width=20, body_size=5.6)
    base.draw_box(ax, (0.365, 0.72), 0.27, 0.17, "Correlation",
                  "200 genes\nPearson centroids",
                  "#EFF7F2", base.GREEN, body_width=18, body_size=5.6)
    base.draw_box(ax, (0.70, 0.72), 0.27, 0.17, "Rescue",
                  "Top-3 pairwise model\nmargin <= 0.002",
                  base.PALE_ORANGE, base.VERMILLION, body_width=20, body_size=5.6)
    base.arrow(ax, (0.30, 0.805), (0.365, 0.805))
    base.arrow(ax, (0.635, 0.805), (0.70, 0.805))
    base.draw_box(ax, (0.12, 0.38), 0.31, 0.18, "Interfaces",
                  "Streamlit\nCommand line",
                  "#F3F4F6", base.GRAY, body_width=16, body_size=5.8)
    base.draw_box(ax, (0.57, 0.38), 0.31, 0.18, "Outputs",
                  "Lobe | Broad\nNetwork | Region*",
                  "#F3F4F6", base.GRAY, body_width=18, body_size=5.8)
    base.arrow(ax, (0.835, 0.72), (0.725, 0.56))
    base.arrow(ax, (0.43, 0.47), (0.57, 0.47))
    ax.text(0.5, 0.17, "Confidence, margin, entropy, marker coverage,\n"
            "rescue status and out-of-scope warnings",
            ha="center", va="center", fontsize=6.5,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white",
                  "edgecolor": base.BLUE, "linewidth": 1})
    ax.text(0.0, 1.02, "A  Software workflow", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom")
    ax.text(0.5, 0.03, "*Exploratory; cerebellum is outside the current reference",
            ha="center", fontsize=5.8, color=base.VERMILLION)

    ax = fig.add_subplot(grid[0, 1])
    labels = ["LOSO", "LOMO"]
    top1 = [0.558, 0.532]
    top3 = [0.880, 0.867]
    x = np.arange(2)
    width = 0.34
    ax.bar(x - width / 2, top1, width, color=base.BLUE, label="Top1")
    ax.bar(x + width / 2, top3, width, color=base.ORANGE, label="Top3")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Network accuracy")
    ax.set_title("B  Internal validation", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=6.5)
    for xpos, values in ((x - width / 2, top1), (x + width / 2, top3)):
        for xi, value in zip(xpos, values):
            ax.text(xi, value + 0.025, f"{value:.1%}", ha="center", fontsize=6.3)

    ax = fig.add_subplot(grid[1, 0])
    metrics = ["AHBA\nNetwork Top1", "AHBA\nNetwork Top3",
               "Glioma broad\nTop3 strict", "Glioma broad\nTop3 tolerant"]
    values = [0.326, 0.554, 0.754, 0.831]
    colors = [base.SKY, base.BLUE, base.GREEN, base.ORANGE]
    bars = ax.barh(np.arange(4), values, color=colors)
    ax.set_yticks(np.arange(4), metrics)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Accuracy / candidate coverage")
    ax.set_title("C  External coarse-resolution evaluation", loc="left",
                 fontweight="bold")
    for bar, value in zip(bars, values):
        ax.text(value + 0.02, bar.get_y() + bar.get_height() / 2,
                f"{value:.1%}", va="center", fontsize=6.3)

    ax = fig.add_subplot(grid[1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    base.draw_box(ax, (0.04, 0.64), 0.40, 0.24, "With truth",
                  "Accuracy\nStrict + tolerant",
                  "#E2F0D9", base.GREEN, body_width=18, body_size=5.8)
    base.draw_box(ax, (0.56, 0.64), 0.40, 0.24, "Without truth",
                  "Transfer stress test\nNo accuracy claim",
                  base.PALE_ORANGE, base.VERMILLION, body_width=19, body_size=5.8)
    base.draw_box(ax, (0.16, 0.22), 0.68, 0.22, "Diagnostics",
                  "Top3 | margin | entropy | coverage\n"
                  "stability | scope warning",
                  "#F3F4F6", base.GRAY, body_width=32, body_size=5.8)
    base.arrow(ax, (0.24, 0.64), (0.38, 0.44))
    base.arrow(ax, (0.76, 0.64), (0.62, 0.44))
    ax.set_title("D  Interpretation policy", loc="left", fontweight="bold")

    fig.suptitle("cfRNA-BrainTrace: implementation and validation",
                 fontsize=11, fontweight="bold", y=0.99)
    save_figure(fig, "Figure1_cfRNA_BrainTrace_application_note")


def export_main_table() -> None:
    TABLE_OUT.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(
        [
            ["Input", "Gene-symbol expression table; optional metadata",
             "Marker overlap and non-zero coverage"],
            ["Scoring", "200-gene Pearson correlation plus pairwise rescue",
             "Ranked Network Top1/Top3 and rescue status"],
            ["Hierarchy", "Lobe, broad anatomy, Network and exact region",
             "Resolution-specific candidates and out-of-scope flags"],
            ["Diagnostics", "Correlation, display probability, margin and entropy",
             "Confidence and transfer-quality audit"],
            ["Interfaces", "Streamlit web app and command-line interface",
             "Interactive review and reproducible batch export"],
            ["Validation policy", "Accuracy only with independent anatomical truth",
             "Unlabelled biofluids reported as stress tests"],
        ],
        columns=["Component", "Implementation", "Principal output"],
    )
    table.to_csv(TABLE_OUT / "Table1_cfRNA_BrainTrace_features.csv",
                 index=False, encoding="utf-8-sig")
    lines = [
        "| " + " | ".join(table.columns) + " |",
        "| " + " | ".join(["---"] * len(table.columns)) + " |",
    ]
    lines.extend("| " + " | ".join(map(str, row)) + " |"
                 for row in table.itertuples(index=False, name=None))
    (TABLE_OUT / "Table1_cfRNA_BrainTrace_features.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    base.setup_style()
    export_main_figure()
    export_main_table()
    print(FIG_OUT)
    print(TABLE_OUT)


if __name__ == "__main__":
    main()
