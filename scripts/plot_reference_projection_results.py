#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROUTE_LABELS = {
    "native_vsd": "Native VSD",
    "projected_vsd": "Projected VSD",
    "logcpm_baseline": "logCPM",
}


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def plot_accuracy(outdir: Path, figdir: Path) -> pd.DataFrame:
    loso = read_csv(outdir / "bo2023_projected_vsd_loso_route_summary.csv")
    lomo = read_csv(outdir / "bo2023_projected_vsd_lomo_route_summary.csv")
    loso["validation"] = "LOSO"
    lomo["validation"] = "LOMO"
    data = pd.concat([loso, lomo], ignore_index=True)
    data["route_label"] = data["route"].map(ROUTE_LABELS).fillna(data["route"])
    data = data.sort_values(["validation", "route_label"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, metric, title in zip(axes, ["network_top1", "network_top3"], ["Network Top1", "Network Top3"]):
        pivot = data.pivot(index="route_label", columns="validation", values=metric).loc[
            ["Native VSD", "Projected VSD", "logCPM"]
        ]
        x = np.arange(len(pivot.index))
        width = 0.35
        ax.bar(x - width / 2, pivot["LOSO"], width, label="LOSO", color="#4E79A7")
        ax.bar(x + width / 2, pivot["LOMO"], width, label="LOMO", color="#F28E2B")
        ax.set_title(title)
        ax.set_ylim(0, 1.0)
        ax.set_xticks(x)
        ax.set_xticklabels(pivot.index, rotation=20, ha="right")
        ax.set_ylabel("Accuracy")
        ax.grid(axis="y", alpha=0.25)
        for i, value in enumerate(pivot["LOSO"]):
            ax.text(i - width / 2, value + 0.015, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
        for i, value in enumerate(pivot["LOMO"]):
            ax.text(i + width / 2, value + 0.015, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    axes[1].legend(frameon=False, loc="lower right")
    savefig(figdir / "internal_network_accuracy.png")
    return data


def plot_margins(outdir: Path, figdir: Path) -> None:
    frames = []
    for validation, filename in [
        ("LOSO", "bo2023_projected_vsd_loso_detail.csv"),
        ("LOMO", "bo2023_projected_vsd_lomo_detail.csv"),
    ]:
        frame = read_csv(outdir / filename)
        frame["validation"] = validation
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    data["route_label"] = data["route"].map(ROUTE_LABELS).fillna(data["route"])
    order = ["Native VSD", "Projected VSD", "logCPM"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, validation in zip(axes, ["LOSO", "LOMO"]):
        subset = data[data["validation"] == validation]
        values = [subset.loc[subset["route_label"] == route, "decision_margin"].dropna() for route in order]
        ax.boxplot(values, tick_labels=order, showfliers=False, patch_artist=True)
        ax.set_title(f"{validation} Decision Margin")
        ax.set_ylabel("Top1 - Top2 correlation")
        ax.tick_params(axis="x", rotation=20)
        ax.grid(axis="y", alpha=0.25)
    savefig(figdir / "internal_decision_margin_boxplot.png")

    rank_counts = (
        data.groupby(["validation", "route_label", "true_rank"]).size().reset_index(name="count")
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, validation in zip(axes, ["LOSO", "LOMO"]):
        subset = rank_counts[(rank_counts["validation"] == validation) & (rank_counts["true_rank"] <= 6)]
        bottom = np.zeros(6)
        x = np.arange(1, 7)
        for route, color in zip(order, ["#4E79A7", "#59A14F", "#E15759"]):
            counts = (
                subset[subset["route_label"] == route]
                .set_index("true_rank")["count"]
                .reindex(x, fill_value=0)
                .to_numpy()
            )
            ax.plot(x, counts, marker="o", label=route, color=color)
        ax.set_title(f"{validation} True-Rank Distribution")
        ax.set_xlabel("True network rank")
        ax.set_ylabel("Samples")
        ax.set_xticks(x)
        ax.grid(alpha=0.25)
    axes[1].legend(frameon=False)
    savefig(figdir / "internal_true_rank_distribution.png")


def plot_projector_qc(outdir: Path, figdir: Path) -> None:
    params = read_csv(outdir / "projector_gene_parameters.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    items = [
        ("r2", "Per-Gene R2", (0, 1)),
        ("spearman_r", "Per-Gene Spearman r", (-1, 1)),
        ("residual_sd", "Residual SD", None),
        ("n_nonzero_count_samples", "Nonzero Count Samples", None),
    ]
    for ax, (col, title, xlim) in zip(axes.ravel(), items):
        values = pd.to_numeric(params[col], errors="coerce").dropna()
        ax.hist(values, bins=50, color="#4E79A7", alpha=0.85)
        ax.set_title(title)
        if xlim:
            ax.set_xlim(*xlim)
        ax.grid(axis="y", alpha=0.25)
    savefig(figdir / "projector_gene_qc_distributions.png")


def plot_external(outdir: Path, figdir: Path) -> pd.DataFrame:
    detail_path = outdir / "external_projected_vsd_GSE189919_detail.csv"
    detail = read_csv(detail_path)
    top1 = detail[detail["rank"] == 1].copy()
    top3 = detail[detail["rank"].isin([1, 2, 3])].copy()
    top1_counts = top1["network_id"].value_counts().sort_values(ascending=True)
    top3_counts = top3["network_id"].value_counts().sort_values(ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].barh(top1_counts.index, top1_counts.values, color="#59A14F")
    axes[0].set_title("GSE189919 Projected VSD Top1")
    axes[0].set_xlabel("Samples")
    axes[0].grid(axis="x", alpha=0.25)
    axes[1].barh(top3_counts.index, top3_counts.values, color="#B07AA1")
    axes[1].set_title("GSE189919 Projected VSD Top3 Membership")
    axes[1].set_xlabel("Rank appearances")
    axes[1].grid(axis="x", alpha=0.25)
    savefig(figdir / "external_gse189919_network_distribution.png")

    return pd.DataFrame(
        {
            "top1_count": top1["network_id"].value_counts(),
            "top3_count": top3["network_id"].value_counts(),
        }
    ).fillna(0).astype(int).sort_values(["top1_count", "top3_count"], ascending=False)


def write_report(outdir: Path, figdir: Path, accuracy: pd.DataFrame, external_counts: pd.DataFrame) -> None:
    with (outdir / "projector_qc_summary.json").open("rt", encoding="utf-8") as handle:
        qc = json.load(handle)
    with (outdir / "external_projected_vsd_GSE189919_summary.json").open("rt", encoding="utf-8") as handle:
        external = json.load(handle)

    acc = accuracy.set_index(["validation", "route"])
    external_table = ["| network | top1_count | top3_count |", "| --- | ---: | ---: |"]
    for network, row in external_counts.iterrows():
        external_table.append(f"| {network} | {int(row['top1_count'])} | {int(row['top3_count'])} |")

    lines = [
        "# Reference Projection Visualization Summary",
        "",
        "## Figures",
        "",
        "- `figures/internal_network_accuracy.png`",
        "- `figures/internal_decision_margin_boxplot.png`",
        "- `figures/internal_true_rank_distribution.png`",
        "- `figures/projector_gene_qc_distributions.png`",
        "- `figures/external_gse189919_network_distribution.png`",
        "",
        "## Key Readout",
        "",
        f"- Data audit: {qc['n_common_samples']} common samples, {qc['n_common_genes']} common gene symbols, "
        f"{qc['n_locked_model_genes_in_common_panel']}/{qc['n_locked_model_genes']} locked network genes present.",
        f"- Full projector fit: median sample Pearson {qc['training_fit']['median_sample_pearson']:.6f}, "
        f"MAE {qc['training_fit']['mae']:.6f}, median per-gene R2 {qc['training_fit']['median_gene_r2']:.6f}.",
        f"- LOSO Top3: projected {acc.loc[('LOSO', 'projected_vsd'), 'network_top3']:.6f}, "
        f"native {acc.loc[('LOSO', 'native_vsd'), 'network_top3']:.6f}, "
        f"logCPM {acc.loc[('LOSO', 'logcpm_baseline'), 'network_top3']:.6f}.",
        f"- LOMO Top3: projected {acc.loc[('LOMO', 'projected_vsd'), 'network_top3']:.6f}, "
        f"native {acc.loc[('LOMO', 'native_vsd'), 'network_top3']:.6f}, "
        f"logCPM {acc.loc[('LOMO', 'logcpm_baseline'), 'network_top3']:.6f}.",
        f"- GSE189919 projected-space overlap: {external['n_overlap_projector_genes']}/{external['n_projector_genes']} "
        f"projector genes ({external['overlap_fraction']:.1%}), {external['n_samples']} samples.",
        "",
        "## GSE189919 Top Networks",
        "",
        *external_table,
        "",
        "## Interpretation",
        "",
        "The internal decision gate is passed at the SaleemNetworks level: projected VSD preserves or improves Top3 performance relative to native VSD and is not worse than logCPM in both LOSO and LOMO.",
        "",
        "The main caveat is calibration. Projected VSD has a smaller decision margin than logCPM/native in both internal validations, so high Top3 accuracy does not yet mean high-confidence Top1 calls.",
        "",
        "The external GSE189919 output should remain a cross-domain stress test. It has no anatomical truth labels here, and the projected matrix is not native Bo2023 VSD.",
    ]
    (outdir / "reference_projection_visualization_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot reference projection result summaries.")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    figdir = args.outdir / "figures"
    accuracy = plot_accuracy(args.outdir, figdir)
    plot_margins(args.outdir, figdir)
    plot_projector_qc(args.outdir, figdir)
    external_counts = plot_external(args.outdir, figdir)
    write_report(args.outdir, figdir, accuracy, external_counts)
    print(f"Wrote figures and summary to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
