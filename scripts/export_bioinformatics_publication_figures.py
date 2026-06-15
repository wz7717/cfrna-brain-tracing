#!/usr/bin/env python
from __future__ import annotations

import json
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "figures_publication"
DATA_OUT = OUT / "source_data"

LOSO_JSON = ROOT / "results" / "bo2023_network_pairwise_correlation_full_loso_819_margin0p002_20260604" / "validation_summary.json"
LOMO_JSON = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_margin0p002_20260604" / "validation_summary.json"
LOMO_MONKEY = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_margin0p002_20260604" / "network_pairwise_top3_per_monkey_metrics.csv"
MARGIN = ROOT / "results" / "bo2023_network_pairwise_margin_threshold_screen_20260604" / "pairwise_margin_threshold_metrics.csv"
AHBA_METRICS = ROOT / "results" / "ahba_human_rnaseq_external_validation_margin0p002_20260604" / "ahba_rnaseq_external_validation_metrics.json"
AHBA_SAMPLE = ROOT / "results" / "ahba_human_rnaseq_external_validation_margin0p002_20260604" / "ahba_rnaseq_external_validation_sample_summary.csv"
IVY_DIST = ROOT / "results" / "ivy_gap_anatomic_rnaseq_tracing_margin0p002_20260604" / "ivy_gap_anatomic_structure_prediction_distributions.csv"
TCGA_SUMMARY = ROOT / "results" / "tcga_gbm_lgg_sample_mri_label_tracing_20260605" / "tcga_gbm_lgg_sample_mri_label_tracing_summary.csv"
TCIA_MATCH = ROOT / "results" / "tcga_rnaseq_tcia_mri_collection_match_20260605" / "tcga_rnaseq_tcia_mri_match_summary.json"


# Color-blind-safe Okabe-Ito palette.
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
YELLOW = "#F0E442"
BLACK = "#222222"
GRAY = "#6B7280"
LIGHT_GRAY = "#E5E7EB"
PALE_BLUE = "#EAF4FA"
PALE_ORANGE = "#FCEFE8"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.titlesize": 10,
            "axes.linewidth": 0.7,
            "axes.edgecolor": BLACK,
            "axes.labelcolor": BLACK,
            "xtick.color": BLACK,
            "ytick.color": BLACK,
            "text.color": BLACK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def wilson(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return (math.nan, math.nan)
    p = successes / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / den
    return center - half, center + half


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def clean_axis(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.6)
        ax.set_axisbelow(True)


def accuracy_axis(ax: plt.Axes, ymax: float = 1.0) -> None:
    ax.set_ylim(0, ymax)
    ticks = np.arange(0, ymax + 0.001, 0.2)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{x:.0%}" for x in ticks])
    ax.set_ylabel("Accuracy")
    clean_axis(ax)


def annotate_bars(ax: plt.Axes, bars, values, dy: float = 0.018, fontsize: float = 7) -> None:
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            float(value) + dy,
            f"{float(value):.1%}",
            ha="center",
            va="bottom",
            fontsize=fontsize,
        )


def save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix, kwargs in (
        (".png", {"dpi": 600}),
        (".pdf", {}),
        (".svg", {}),
    ):
        fig.savefig(OUT / f"{stem}{suffix}", bbox_inches="tight", pad_inches=0.04, **kwargs)
    plt.close(fig)


def draw_box(ax, xy, width, height, title, body, facecolor, edgecolor=BLACK, body_width=30, body_size=6.8):
    box = patches.FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=0.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
    )
    ax.add_patch(box)
    x, y = xy
    ax.text(x + width / 2, y + height * 0.70, title, ha="center", va="center", fontsize=8.2, fontweight="bold")
    ax.text(
        x + width / 2,
        y + height * 0.36,
        "\n".join(textwrap.wrap(body, width=body_width)),
        ha="center",
        va="center",
        fontsize=body_size,
        color="#374151",
        linespacing=1.15,
    )


def arrow(ax, start, end, color=GRAY):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": color, "shrinkA": 2, "shrinkB": 2},
    )


def figure1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    draw_box(ax, (0.03, 0.72), 0.25, 0.19, "Primate reference atlas", "819 macaque brain RNA-seq samples\n110 regions; 10 Networks", PALE_BLUE, BLUE)
    draw_box(ax, (0.375, 0.72), 0.25, 0.19, "Fold-local model", "Top200 discriminative genes\nPearson correlation", "#EFF7F2", GREEN)
    draw_box(ax, (0.72, 0.72), 0.25, 0.19, "Boundary rescue", "Top3 pairwise correlation\nmargin = 0.002", PALE_ORANGE, VERMILLION)
    arrow(ax, (0.28, 0.815), (0.375, 0.815))
    arrow(ax, (0.625, 0.815), (0.72, 0.815))

    draw_box(ax, (0.05, 0.40), 0.25, 0.17, "Primary endpoint", "Network Top1 and Top3", "#DDEBF7", BLUE)
    draw_box(ax, (0.375, 0.40), 0.25, 0.17, "Secondary endpoint", "Within-Network adaptive region group", "#E2F0D9", GREEN)
    draw_box(ax, (0.70, 0.40), 0.25, 0.17, "Exploratory endpoint", "Top3 Network beam exact region", "#FCE4D6", ORANGE)
    arrow(ax, (0.845, 0.72), (0.175, 0.57))
    arrow(ax, (0.30, 0.485), (0.375, 0.485))
    arrow(ax, (0.625, 0.485), (0.70, 0.485))

    draw_box(ax, (0.015, 0.075), 0.22, 0.16, "Internal validation", "Strict LOSO and leave-one-monkey-out", "#F3F4F6", body_width=22, body_size=6.3)
    draw_box(ax, (0.265, 0.075), 0.22, 0.16, "Cross-species", "AHBA normal-brain label harmonization", "#F3F4F6", body_width=22, body_size=6.3)
    draw_box(ax, (0.515, 0.075), 0.22, 0.16, "Disease domain", "Ivy GAP and TCGA prediction distributions", "#F3F4F6", body_width=22, body_size=6.3)
    draw_box(ax, (0.765, 0.075), 0.22, 0.16, "Draft B validation", "TCGA RNA-seq plus TCIA MRI location truth", "#F3F4F6", body_width=22, body_size=6.3)
    ax.plot([0.175, 0.175], [0.40, 0.31], color=GRAY, lw=1.0)
    ax.plot([0.125, 0.875], [0.31, 0.31], color=GRAY, lw=1.0)
    for x in (0.125, 0.375, 0.625, 0.875):
        arrow(ax, (x, 0.31), (x, 0.235))

    ax.text(0.5, 0.965, "Hierarchical brain-origin RNA tracing and validation design", ha="center", va="top", fontsize=10, fontweight="bold")
    ax.text(
        0.5,
        0.015,
        "cfRNA is a prospective application; the current evidence is derived primarily from tissue RNA-seq.",
        ha="center",
        va="bottom",
        fontsize=7,
        color=VERMILLION,
        fontweight="bold",
    )
    save_figure(fig, "Figure1_study_design")


def build_endpoint_data() -> pd.DataFrame:
    loso = load_json(LOSO_JSON)
    lomo = load_json(LOMO_JSON)
    rows = [
        ["LOSO", "Network", 819, 457, 0.557997557997558, 721, 0.8803418803418803, "Primary"],
        ["LOSO", "Region Group", 814, round(814 * 0.43611793611793614), 0.43611793611793614, round(814 * 0.6928746928746928), 0.6928746928746928, "Secondary"],
        ["LOSO", "Exact Region", 814, round(814 * 0.22235872235872237), 0.22235872235872237, round(814 * 0.44471744471744473), 0.44471744471744473, "Exploratory"],
        ["LOMO", "Network", 819, lomo["routes"]["network"]["network_pairwise_correlation_rescue_top3"]["top1_hits"], lomo["routes"]["network"]["network_pairwise_correlation_rescue_top3"]["top1_accuracy"], lomo["routes"]["network"]["network_pairwise_correlation_rescue_top3"]["top3_hits"], lomo["routes"]["network"]["network_pairwise_correlation_rescue_top3"]["top3_accuracy"], "Primary"],
        ["LOMO", "Region Group", 812, lomo["routes"]["resolution_group"]["top3_network_beam_local_region_candidates"]["group_top1_hits"], lomo["routes"]["resolution_group"]["top3_network_beam_local_region_candidates"]["group_top1_accuracy"], lomo["routes"]["resolution_group"]["top3_network_beam_local_region_candidates"]["group_top3_hits"], lomo["routes"]["resolution_group"]["top3_network_beam_local_region_candidates"]["group_top3_accuracy"], "Secondary"],
        ["LOMO", "Exact Region", 812, lomo["routes"]["exact_region"]["top3_beam_local_top50_top100_zfusion_w0p25"]["top1_hits"], lomo["routes"]["exact_region"]["top3_beam_local_top50_top100_zfusion_w0p25"]["top1_accuracy"], lomo["routes"]["exact_region"]["top3_beam_local_top50_top100_zfusion_w0p25"]["top3_hits"], lomo["routes"]["exact_region"]["top3_beam_local_top50_top100_zfusion_w0p25"]["top3_accuracy"], "Exploratory"],
    ]
    return pd.DataFrame(rows, columns=["validation", "endpoint", "n", "top1_hits", "top1", "top3_hits", "top3", "role"])


def figure2_primary_results() -> None:
    loso = load_json(LOSO_JSON)
    lomo = load_json(LOMO_JSON)
    endpoint = build_endpoint_data()
    endpoint.to_csv(DATA_OUT / "Figure2_endpoint_metrics.csv", index=False)

    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95], hspace=0.42, wspace=0.34)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    # A: primary Network baseline versus rescue.
    labels = ["LOSO", "LOMO"]
    baseline = [
        loso["routes"]["network_discriminative_correlation_top200"]["top1_accuracy"],
        lomo["routes"]["network"]["network_discriminative_correlation_top200"]["top1_accuracy"],
    ]
    rescued = [
        loso["routes"]["network_pairwise_correlation_rescue_top3"]["top1_accuracy"],
        lomo["routes"]["network"]["network_pairwise_correlation_rescue_top3"]["top1_accuracy"],
    ]
    x = np.arange(2)
    width = 0.34
    b1 = ax_a.bar(x - width / 2, baseline, width, color=GRAY, label="Top200 baseline")
    b2 = ax_a.bar(x + width / 2, rescued, width, color=BLUE, label="Pairwise rescue")
    annotate_bars(ax_a, b1, baseline)
    annotate_bars(ax_a, b2, rescued)
    ax_a.set_xticks(x, labels)
    ax_a.set_title("Network Top1")
    accuracy_axis(ax_a, 0.8)
    ax_a.legend(frameon=False, loc="upper right")
    ax_a.text(0, 0.72, "P=0.0026", ha="center", fontsize=7)
    panel_label(ax_a, "A")

    # B: paired changes in LOSO.
    changes = [45, 20, 34]
    change_labels = ["Corrected\nerrors", "Introduced\nerrors", "No net\nbenefit"]
    colors = [GREEN, VERMILLION, LIGHT_GRAY]
    bars = ax_b.bar(change_labels, changes, color=colors)
    ax_b.set_ylabel("Number of switched samples")
    ax_b.set_title("LOSO rescue changes (99 switches)")
    clean_axis(ax_b)
    for bar, value in zip(bars, changes):
        ax_b.text(bar.get_x() + bar.get_width() / 2, value + 1, str(value), ha="center", fontsize=7)
    ax_b.text(0.98, 0.94, "Net +25 correct", transform=ax_b.transAxes, ha="right", va="top", color=GREEN, fontweight="bold")
    panel_label(ax_b, "B")

    # C: hierarchical endpoint summary.
    endpoint_order = ["Network", "Region Group", "Exact Region"]
    xpos = np.arange(3)
    width = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5]
    combos = [
        ("LOSO", "top1", "LOSO Top1", BLUE),
        ("LOSO", "top3", "LOSO Top3", SKY),
        ("LOMO", "top1", "LOMO Top1", ORANGE),
        ("LOMO", "top3", "LOMO Top3", YELLOW),
    ]
    for offset, (validation, col, legend, color) in zip(offsets, combos):
        sub = endpoint[endpoint["validation"] == validation].set_index("endpoint").loc[endpoint_order]
        vals = sub[col].to_numpy(float)
        bars = ax_c.bar(xpos + offset * width, vals, width, label=legend, color=color, edgecolor="white", linewidth=0.4)
        annotate_bars(ax_c, bars, vals, dy=0.012, fontsize=6)
    ax_c.set_xticks(xpos, ["Network\n(primary)", "Region Group\n(secondary)", "Exact Region\n(exploratory)"])
    ax_c.set_title("Hierarchical performance")
    accuracy_axis(ax_c, 1.0)
    ax_c.legend(frameon=False, ncols=4, loc="upper center", bbox_to_anchor=(0.5, -0.20))
    panel_label(ax_c, "C")

    fig.suptitle("Internal validation supports Network as the primary endpoint", y=0.995, fontweight="bold")
    save_figure(fig, "Figure2_internal_validation")


def figure3_threshold_and_lomo() -> None:
    margin = pd.read_csv(MARGIN)
    monkey = pd.read_csv(LOMO_MONKEY)
    margin.to_csv(DATA_OUT / "Figure3_margin_screen.csv", index=False)
    monkey.to_csv(DATA_OUT / "Figure3_lomo_per_monkey.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), gridspec_kw={"width_ratios": [1.05, 1.35]})
    ax_a, ax_b = axes

    # A: post-hoc threshold screen.
    x = margin["threshold"].to_numpy(float)
    top1 = margin["top1_accuracy"].to_numpy(float)
    switches = margin["n_switches"].to_numpy(float)
    ax_a.plot(x, top1, marker="o", color=BLUE, linewidth=1.4, markersize=3.5)
    ax_a.axhline(432 / 819, color=GRAY, linestyle="--", linewidth=0.9, label="Top200 baseline")
    best_idx = int(np.argmin(np.abs(x - 0.002)))
    ax_a.scatter([x[best_idx]], [top1[best_idx]], s=45, facecolor=ORANGE, edgecolor=BLACK, linewidth=0.6, zorder=5)
    ax_a.annotate(
        "Selected: 0.002\nTop1 55.8%",
        (x[best_idx], top1[best_idx]),
        xytext=(0.018, 0.575),
        arrowprops={"arrowstyle": "->", "lw": 0.8, "color": BLACK},
        fontsize=7,
        ha="left",
    )
    ax_a.set_xscale("symlog", linthresh=0.001)
    ax_a.set_xlabel("Minimum rescue margin")
    ax_a.set_ylabel("LOSO Network Top1")
    ax_a.set_ylim(0.51, 0.59)
    ax_a.set_yticks([0.52, 0.54, 0.56, 0.58], ["52%", "54%", "56%", "58%"])
    clean_axis(ax_a)
    ax_a.legend(frameon=False, loc="lower right")
    ax_a.set_title("Post-hoc threshold screen")
    panel_label(ax_a, "A")

    ax_a2 = ax_a.twinx()
    ax_a2.plot(x, switches, color=VERMILLION, linewidth=0.9, alpha=0.65)
    ax_a2.set_ylabel("Switched samples", color=VERMILLION)
    ax_a2.tick_params(axis="y", colors=VERMILLION)
    ax_a2.spines["top"].set_visible(False)

    # B: LOMO individual performance and sample size.
    monkey = monkey.sort_values("hit1_mean")
    y = np.arange(len(monkey))
    lo = []
    hi = []
    for hits, n in zip(monkey["hit1_hits"], monkey["n"]):
        a, b = wilson(int(hits), int(n))
        lo.append(a)
        hi.append(b)
    vals = monkey["hit1_mean"].to_numpy(float)
    xerr = np.vstack([vals - np.array(lo), np.array(hi) - vals])
    sizes = 25 + 100 * np.sqrt(monkey["n"].to_numpy(float) / monkey["n"].max())
    ax_b.errorbar(vals, y, xerr=xerr, fmt="none", ecolor="#9CA3AF", elinewidth=1.0, capsize=2)
    ax_b.scatter(vals, y, s=sizes, c=vals, cmap="viridis", vmin=0.35, vmax=0.75, edgecolor="white", linewidth=0.6, zorder=3)
    ax_b.axvline(436 / 819, color=BLACK, linestyle="--", linewidth=0.9, label="Aggregate 53.2%")
    ax_b.set_yticks(y, [f"{m}  (n={n})" for m, n in zip(monkey["monkey_id"], monkey["n"])])
    ax_b.set_xlim(0.25, 0.90)
    ax_b.set_xticks([0.3, 0.5, 0.7, 0.9], ["30%", "50%", "70%", "90%"])
    ax_b.set_xlabel("Network Top1 accuracy (95% Wilson CI)")
    ax_b.set_title("Leave-one-monkey-out heterogeneity")
    clean_axis(ax_b, grid=False)
    ax_b.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)
    ax_b.legend(frameon=False, loc="lower right")
    panel_label(ax_b, "B")

    fig.suptitle("Threshold selection and cross-individual robustness", y=1.02, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "Figure3_threshold_lomo")


def figure4_ahba() -> None:
    metrics = load_json(AHBA_METRICS)
    sample = pd.read_csv(AHBA_SAMPLE)
    supported = sample[sample["supported_for_accuracy"].astype(str).str.lower().eq("true")].copy()

    metric_df = pd.DataFrame(
        [
            ["Network Top1", 233, metrics["network_top1_accuracy_coarse"], BLUE],
            ["Network Top3", 233, metrics["network_top3_accuracy_coarse"], SKY],
            ["Broad anatomy Top1", 233, metrics["region_lobe_top1_accuracy_coarse"], GREEN],
            ["Exact Region Top1", 91, metrics["region_top1_exact_accuracy_on_exact_mapped_labels"], ORANGE],
            ["Exact Region Top3", 91, metrics["region_top3_exact_accuracy_on_exact_mapped_labels"], YELLOW],
        ],
        columns=["metric", "n", "accuracy", "color"],
    )
    metric_df.drop(columns="color").to_csv(DATA_OUT / "Figure4_ahba_metrics.csv", index=False)

    true_order = [
        "frontal lobe",
        "parietal lobe",
        "temporal lobe",
        "occipital lobe",
        "cingulate cortex",
        "insula",
        "parahippocampal cortex",
        "striatum",
        "globus pallidus",
    ]
    pred_order = ["frontal lobe", "cingulate cortex", "insula", "subcortical structures"]
    confusion = pd.crosstab(supported["public_major_anatomy"], supported["predicted_public_major_top1"], normalize="index")
    confusion = confusion.reindex(index=[x for x in true_order if x in confusion.index], columns=pred_order, fill_value=0)
    confusion.to_csv(DATA_OUT / "Figure4_ahba_broad_confusion.csv")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [0.9, 1.25]})
    bars = ax_a.barh(np.arange(len(metric_df)), metric_df["accuracy"], color=metric_df["color"])
    ax_a.set_yticks(np.arange(len(metric_df)), [f"{m}\n(n={n})" for m, n in zip(metric_df["metric"], metric_df["n"])])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 0.7)
    ax_a.set_xticks([0, 0.2, 0.4, 0.6], ["0%", "20%", "40%", "60%"])
    ax_a.set_xlabel("Label-harmonized accuracy")
    ax_a.set_title("Cross-species summary")
    clean_axis(ax_a, grid=False)
    ax_a.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)
    for bar, val in zip(bars, metric_df["accuracy"]):
        ax_a.text(val + 0.015, bar.get_y() + bar.get_height() / 2, f"{val:.1%}", va="center", fontsize=7)
    panel_label(ax_a, "A")

    data = confusion.to_numpy(float)
    im = ax_b.imshow(data, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax_b.set_xticks(np.arange(len(confusion.columns)), [x.replace(" structures", "\nstructures") for x in confusion.columns], rotation=30, ha="right")
    ax_b.set_yticks(np.arange(len(confusion.index)), confusion.index)
    ax_b.set_xlabel("Predicted broad anatomy")
    ax_b.set_ylabel("AHBA label")
    ax_b.set_title("Broad-anatomy prediction pattern")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i, j] >= 0.05:
                ax_b.text(j, i, f"{data[i, j]:.0%}", ha="center", va="center", fontsize=6.5, color="white" if data[i, j] > 0.55 else BLACK)
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.03)
    cbar.set_label("Row fraction")
    panel_label(ax_b, "B")

    fig.suptitle("Cross-species transfer to normal human brain RNA-seq", y=1.01, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "Figure4_ahba_external_validation")


def short_network(name: str) -> str:
    mapping = {
        "Subcortical": "Subcortical",
        "Hippocampal formation": "Hippocampal",
        "Orbitomedial Prefrontal Cortex (OMPFC)": "OMPFC",
        "Frontal (agranular frontal motor areas)": "Motor frontal",
        "Cingulate gyrus": "Cingulate",
        "Operculum/Insula": "Operculum/Insula",
        "Occipital/Temporal": "Occipital/Temporal",
        "Parietal, and Parieto-occipital region": "Parietal",
    }
    return mapping.get(name, name)


def figure5_tumor_and_mri() -> None:
    ivy = pd.read_csv(IVY_DIST)
    ivy = ivy[ivy["endpoint"] == "network_top1"].copy()
    tcga = pd.read_csv(TCGA_SUMMARY)
    match = load_json(TCIA_MATCH)

    # Ivy stacked distribution.
    ivy["network_short"] = ivy["value"].map(short_network)
    ivy_pivot = ivy.pivot_table(index="structure_acronym", columns="network_short", values="fraction", fill_value=0)
    ivy_order = ["CT-reference-histology", "CTmvp-reference-histology", "CTpan-reference-histology", "IT-reference-histology", "LE-reference-histology"]
    ivy_pivot = ivy_pivot.reindex(ivy_order)
    ivy_pivot.to_csv(DATA_OUT / "Figure5_ivy_network_distribution.csv")

    # TCGA stacked distribution.
    tcga_counts = pd.crosstab(tcga["project_id"], tcga["network_top1"])
    tcga_frac = tcga_counts.div(tcga_counts.sum(axis=1), axis=0)
    tcga_frac.columns = [short_network(x) for x in tcga_frac.columns]
    tcga_frac.to_csv(DATA_OUT / "Figure5_tcga_network_distribution.csv")

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.90], hspace=0.67, wspace=0.33)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    network_colors = {
        "Subcortical": BLUE,
        "Hippocampal": ORANGE,
        "OMPFC": GREEN,
        "Motor frontal": PURPLE,
        "Cingulate": VERMILLION,
        "Operculum/Insula": SKY,
        "Occipital/Temporal": YELLOW,
        "Parietal": "#999999",
    }

    left = np.zeros(len(ivy_pivot))
    for col in ivy_pivot.columns:
        vals = ivy_pivot[col].to_numpy(float)
        ax_a.barh(np.arange(len(ivy_pivot)), vals, left=left, color=network_colors.get(col, "#BDBDBD"), label=col, height=0.68)
        left += vals
    ax_a.set_yticks(np.arange(len(ivy_pivot)), ["Cellular tumor", "Microvascular proliferation", "Pseudopalisading cells", "Infiltrating tumor", "Leading edge"])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 1)
    ax_a.set_xticks(np.linspace(0, 1, 6), [f"{x:.0%}" for x in np.linspace(0, 1, 6)])
    ax_a.set_xlabel("Fraction of samples")
    ax_a.set_title("Ivy GAP: Network Top1 prediction distribution (not anatomical accuracy)")
    ax_a.legend(frameon=False, ncols=4, loc="upper center", bbox_to_anchor=(0.5, -0.19), columnspacing=1.2, handlelength=1.4)
    clean_axis(ax_a, grid=False)
    ax_a.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)
    panel_label(ax_a, "A")

    left = np.zeros(len(tcga_frac))
    for col in tcga_frac.columns:
        vals = tcga_frac[col].to_numpy(float)
        ax_b.bar(np.arange(len(tcga_frac)), vals, bottom=left, color=network_colors.get(col, "#BDBDBD"), width=0.62)
        left += vals
    ax_b.set_xticks(np.arange(len(tcga_frac)), tcga_frac.index)
    ax_b.set_ylim(0, 1)
    ax_b.set_yticks(np.linspace(0, 1, 6), [f"{x:.0%}" for x in np.linspace(0, 1, 6)])
    ax_b.set_ylabel("Fraction of samples")
    ax_b.set_title("TCGA Network Top1 distribution\n(no MRI truth yet)", pad=8)
    clean_axis(ax_b)
    panel_label(ax_b, "B")

    stages = ["RNA-seq\npatients", "MRI\nmatched", "Segmentation-\nready", "Complete\n4-modality"]
    values = [match["n_rnaseq_patients"], match["n_matched_patients"], match["modality_counts_among_rnaseq_patients"]["minimal_segmentation_ready"], match["modality_counts_among_rnaseq_patients"]["complete_brats4"]]
    colors = [GRAY, SKY, GREEN, BLUE]
    bars = ax_c.bar(stages, values, color=colors, width=0.68)
    ax_c.set_ylabel("Patients")
    ax_c.set_title("TCGA-TCIA validation cohort")
    clean_axis(ax_c)
    for bar, value in zip(bars, values):
        ax_c.text(bar.get_x() + bar.get_width() / 2, value + 14, f"{value}\n({value / values[0]:.1%})", ha="center", va="bottom", fontsize=7)
    ax_c.set_ylim(0, 930)
    panel_label(ax_c, "C")

    fig.suptitle("Glioma domain shift and construction of the MRI validation cohort", y=0.995, fontweight="bold")
    save_figure(fig, "Figure5_glioma_domain_mri_cohort")


def supplementary_figures() -> None:
    endpoint = build_endpoint_data()
    lomo = load_json(LOMO_JSON)
    monkey_exact = pd.read_csv(ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_margin0p002_20260604" / "exact_region_per_monkey_metrics.csv")
    monkey_group = pd.read_csv(ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_margin0p002_20260604" / "resolution_group_per_monkey_metrics.csv")
    monkey_network = pd.read_csv(LOMO_MONKEY)
    per = monkey_network[["monkey_id", "n", "hit1_mean", "hit3_mean"]].rename(columns={"hit1_mean": "Network Top1", "hit3_mean": "Network Top3"})
    per = per.merge(monkey_group[["monkey_id", "group_hit1_mean", "group_hit3_mean"]].rename(columns={"group_hit1_mean": "Region Group Top1", "group_hit3_mean": "Region Group Top3"}), on="monkey_id")
    per = per.merge(monkey_exact[["monkey_id", "hit1_mean", "hit3_mean"]].rename(columns={"hit1_mean": "Exact Region Top1", "hit3_mean": "Exact Region Top3"}), on="monkey_id")
    per.to_csv(DATA_OUT / "FigureS1_lomo_hierarchical_per_monkey.csv", index=False)

    data_cols = ["Network Top1", "Network Top3", "Region Group Top1", "Region Group Top3", "Exact Region Top1", "Exact Region Top3"]
    data = per.set_index("monkey_id")[data_cols].T.to_numpy(float)
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    im = ax.imshow(data, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(per)), [f"{m}\n(n={n})" for m, n in zip(per["monkey_id"], per["n"])])
    ax.set_yticks(np.arange(len(data_cols)), data_cols)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{data[i, j]:.0%}", ha="center", va="center", fontsize=6, color="white" if data[i, j] < 0.35 or data[i, j] > 0.75 else BLACK)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Accuracy")
    ax.set_title("LOMO performance across endpoint resolutions")
    fig.tight_layout()
    save_figure(fig, "FigureS1_lomo_hierarchical_heatmap")

    # Supplementary paired changes by threshold.
    margin = pd.read_csv(MARGIN)
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    threshold_x = np.arange(len(margin))
    ax.plot(threshold_x, margin["top1_gains_vs_top200"], marker="o", color=GREEN, label="Top1 gains")
    ax.plot(threshold_x, margin["top1_losses_vs_top200"], marker="o", color=VERMILLION, label="Top1 losses")
    ax.plot(threshold_x, margin["n_switches"], marker="o", color=GRAY, label="All switches")
    selected_idx = int(np.flatnonzero(np.isclose(margin["threshold"], 0.002))[0])
    ax.axvline(selected_idx, color=ORANGE, linestyle="--", linewidth=1.0)
    ax.set_xticks(threshold_x)
    ax.set_xticklabels([f"{value:g}" for value in margin["threshold"]], rotation=45, ha="right")
    ax.set_xlabel("Minimum rescue margin")
    ax.set_ylabel("Samples")
    ax.set_title("Pairwise rescue changes across the post-hoc threshold screen")
    clean_axis(ax)
    ax.legend(frameon=False, ncols=3)
    fig.tight_layout()
    save_figure(fig, "FigureS2_margin_switch_audit")


def write_captions() -> None:
    text = """# Publication figure plan and legends

## Main figures

### Figure 1. Hierarchical brain-origin RNA tracing and validation design
The macaque reference atlas is used to build fold-local discriminative-correlation models. Network is the prespecified primary endpoint; within-Network region group is secondary and exact region is exploratory. Internal validation uses strict leave-one-sample-out (LOSO) and leave-one-monkey-out (LOMO) designs. AHBA provides label-harmonized cross-species evidence, whereas Ivy GAP and TCGA characterize disease-domain behavior. TCGA-TCIA transcriptome-MRI validation is reserved for Draft B. cfRNA is a prospective application; the current evidence is derived primarily from tissue RNA-seq.

### Figure 2. Internal validation supports Network as the primary endpoint
(A) Network Top1 accuracy for the Top200 discriminative-correlation baseline and Top3-constrained pairwise rescue. The LOSO paired comparison used the same 819 held-out samples (two-sided paired P=0.0026). (B) Outcome of the 99 LOSO rescue switches. Forty-five corrected baseline errors, 20 introduced errors and 34 produced no net correctness change, yielding a net gain of 25 correct samples. (C) Top1 and Top3 accuracy across the endpoint hierarchy. Network is primary, Region Group secondary and Exact Region exploratory. Exact evaluable sample counts were 814 for LOSO and 812 for LOMO.

### Figure 3. Threshold selection and cross-individual robustness
(A) Post-hoc screen of the minimum pairwise-rescue margin using existing strict LOSO predictions without model refitting. The selected threshold of 0.002 maximized Top1 accuracy while preserving Top3. The red line denotes the number of switched samples. This panel documents retrospective threshold selection and is not an independent validation. (B) Network Top1 accuracy for each held-out monkey in LOMO validation. Point size represents held-out sample count and horizontal intervals are 95% Wilson confidence intervals. The dashed line denotes aggregate LOMO accuracy (53.2%).

### Figure 4. Cross-species transfer to normal human brain RNA-seq
(A) Label-harmonized AHBA accuracy. Network and broad-anatomy metrics were evaluated in 233 supported samples; Exact Region metrics were restricted to 91 exact-mapped samples. (B) Row-normalized broad-anatomy prediction pattern. AHBA labels and macaque-derived predictions are harmonized coarse categories and are not directly equivalent to the macaque internal-validation labels.

### Figure 5. Glioma domain shift and construction of the MRI validation cohort
(A) Ivy GAP Network Top1 prediction distributions across five glioblastoma microanatomical structures. Ivy labels describe tumor compartments rather than normal brain location; the panel is not an accuracy analysis. (B) TCGA-GBM/LGG Network Top1 distributions before MRI-label generation. (C) Attrition from 800 TCGA RNA-seq patients to the MRI-matched, segmentation-ready and complete four-modality validation cohorts. TCGA anatomical accuracy remains unreported until MRI-derived truth is available.

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
- Accuracy is shown as a proportion or percentage with a zero baseline unless a deliberately truncated diagnostic axis is clearly labeled.
- All figure source tables are stored in `source_data`.
"""
    (OUT / "figure_legends_and_style.md").write_text(text, encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    setup_style()
    figure1_workflow()
    figure2_primary_results()
    figure3_threshold_and_lomo()
    figure4_ahba()
    figure5_tumor_and_mri()
    supplementary_figures()
    write_captions()
    print(f"Publication figures written to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
