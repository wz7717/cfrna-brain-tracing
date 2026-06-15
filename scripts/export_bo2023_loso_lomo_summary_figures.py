#!/usr/bin/env python
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "bo2023_loso_lomo_summary_figures_20260601"

NETWORK_LOSO = ROOT / "results" / "bo2023_network_pairwise_correlation_full_loso_819_rerun_20260526" / "validation_summary.json"
EXACT_LOSO = ROOT / "data" / "models" / "bo2023_exact_region_validation_summary.json"
GROUP_LOSO = ROOT / "data" / "models" / "bo2023_region_resolution_adaptive_max8_validation_summary.json"
LOMO = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_20260601" / "validation_summary.json"
LOMO_NETWORK_PER_MONKEY = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_20260601" / "network_pairwise_top3_per_monkey_metrics.csv"
LOMO_EXACT_PER_MONKEY = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_20260601" / "exact_region_per_monkey_metrics.csv"
LOMO_GROUP_PER_MONKEY = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_20260601" / "resolution_group_per_monkey_metrics.csv"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def pct(value: float) -> str:
    return f"{value * 100:.1f}%"


def label_bars(ax: plt.Axes, bars) -> None:
    for bar in bars:
        value = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            pct(float(value)),
            ha="center",
            va="bottom",
            fontsize=9,
        )


def style_accuracy_axis(ax: plt.Axes) -> None:
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Accuracy")
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_yticklabels([pct(x) for x in np.linspace(0, 1, 6)])
    ax.grid(axis="y", alpha=0.24)
    ax.set_axisbelow(True)


def build_summary_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    network_loso = load_json(NETWORK_LOSO)
    exact_loso = load_json(EXACT_LOSO)
    group_loso = load_json(GROUP_LOSO)
    lomo = load_json(LOMO)

    loso_network_base = network_loso["routes"]["network_discriminative_correlation_top200"]
    loso_network_rescue = network_loso["routes"]["network_pairwise_correlation_rescue_top3"]
    loso_exact = exact_loso["routes"]["top3_beam_local_top50_top100_zfusion_w0p25"]
    loso_group = group_loso["routes"]["top3_network_beam_local_region_candidates"]

    lomo_network_base = lomo["routes"]["network"]["network_discriminative_correlation_top200"]
    lomo_network_rescue = lomo["routes"]["network"]["network_pairwise_correlation_rescue_top3"]
    lomo_exact = lomo["routes"]["exact_region"]["top3_beam_local_top50_top100_zfusion_w0p25"]
    lomo_group = lomo["routes"]["resolution_group"]["top3_network_beam_local_region_candidates"]

    endpoint_rows = [
        {
            "validation": "LOSO",
            "endpoint": "Network",
            "n": loso_network_rescue["n"],
            "top1": loso_network_rescue["top1_accuracy"],
            "top3": loso_network_rescue["top3_accuracy"],
            "notes": "Top3 pairwise rescue",
        },
        {
            "validation": "LOSO",
            "endpoint": "Exact Region",
            "n": loso_exact["n"],
            "top1": loso_exact["top1_accuracy"],
            "top3": loso_exact["top3_accuracy"],
            "notes": "Top3 Network beam + Top50/Top100 z-fusion w=0.25",
        },
        {
            "validation": "LOSO",
            "endpoint": "Region Group",
            "n": loso_group["n"],
            "top1": loso_group["group_top1_accuracy"],
            "top3": loso_group["group_top3_accuracy"],
            "notes": "same-Network adaptive max8",
        },
        {
            "validation": "LOMO",
            "endpoint": "Network",
            "n": lomo_network_rescue["n"],
            "top1": lomo_network_rescue["top1_accuracy"],
            "top3": lomo_network_rescue["top3_accuracy"],
            "notes": "Top3 pairwise rescue",
        },
        {
            "validation": "LOMO",
            "endpoint": "Exact Region",
            "n": lomo_exact["n"],
            "top1": lomo_exact["top1_accuracy"],
            "top3": lomo_exact["top3_accuracy"],
            "notes": "Top3 Network beam + Top50/Top100 z-fusion w=0.25",
        },
        {
            "validation": "LOMO",
            "endpoint": "Region Group",
            "n": lomo_group["n"],
            "top1": lomo_group["group_top1_accuracy"],
            "top3": lomo_group["group_top3_accuracy"],
            "notes": "same-Network adaptive max8",
        },
    ]
    network_rows = [
        {
            "validation": "LOSO",
            "route": "Top200 baseline",
            "n": loso_network_base["n"],
            "top1": loso_network_base["top1_accuracy"],
            "top3": loso_network_base["top3_accuracy"],
        },
        {
            "validation": "LOSO",
            "route": "Top3 pairwise rescue",
            "n": loso_network_rescue["n"],
            "top1": loso_network_rescue["top1_accuracy"],
            "top3": loso_network_rescue["top3_accuracy"],
        },
        {
            "validation": "LOMO",
            "route": "Top200 baseline",
            "n": lomo_network_base["n"],
            "top1": lomo_network_base["top1_accuracy"],
            "top3": lomo_network_base["top3_accuracy"],
        },
        {
            "validation": "LOMO",
            "route": "Top3 pairwise rescue",
            "n": lomo_network_rescue["n"],
            "top1": lomo_network_rescue["top1_accuracy"],
            "top3": lomo_network_rescue["top3_accuracy"],
        },
    ]
    return pd.DataFrame(endpoint_rows), pd.DataFrame(network_rows)


def plot_endpoint_summary(endpoint_df: pd.DataFrame) -> None:
    endpoints = ["Network", "Exact Region", "Region Group"]
    validations = ["LOSO", "LOMO"]
    colors = {"Top1": "#0072B2", "Top3": "#009E73"}
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), constrained_layout=True, sharey=True)
    x = np.arange(len(endpoints))
    width = 0.36
    for ax, validation in zip(axes, validations):
        sub = endpoint_df[endpoint_df["validation"] == validation].set_index("endpoint").loc[endpoints]
        bars1 = ax.bar(x - width / 2, sub["top1"], width, label="Top1", color=colors["Top1"])
        bars3 = ax.bar(x + width / 2, sub["top3"], width, label="Top3", color=colors["Top3"])
        ax.set_title(f"{validation}: formal hierarchical endpoints", fontweight="bold")
        ax.set_xticks(x, endpoints)
        style_accuracy_axis(ax)
        label_bars(ax, bars1)
        label_bars(ax, bars3)
        ax.legend(loc="upper right")
    fig.savefig(OUTDIR / "loso_lomo_endpoint_summary.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUTDIR / "loso_lomo_endpoint_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_network_rescue(network_df: pd.DataFrame) -> None:
    validations = ["LOSO", "LOMO"]
    routes = ["Top200 baseline", "Top3 pairwise rescue"]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2), constrained_layout=True, sharey=True)
    x = np.arange(len(routes))
    width = 0.36
    for ax, validation in zip(axes, validations):
        sub = network_df[network_df["validation"] == validation].set_index("route").loc[routes]
        bars1 = ax.bar(x - width / 2, sub["top1"], width, label="Top1", color="#0072B2")
        bars3 = ax.bar(x + width / 2, sub["top3"], width, label="Top3", color="#009E73")
        ax.set_title(f"{validation}: Network baseline vs rescue", fontweight="bold")
        ax.set_xticks(x, routes)
        style_accuracy_axis(ax)
        label_bars(ax, bars1)
        label_bars(ax, bars3)
        ax.legend(loc="upper right")
    fig.savefig(OUTDIR / "network_baseline_vs_pairwise_rescue.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUTDIR / "network_baseline_vs_pairwise_rescue.pdf", bbox_inches="tight")
    plt.close(fig)


def build_per_monkey_table() -> pd.DataFrame:
    network = pd.read_csv(LOMO_NETWORK_PER_MONKEY)[["monkey_id", "n", "hit1_mean", "hit3_mean"]].rename(
        columns={"hit1_mean": "network_top1", "hit3_mean": "network_top3"}
    )
    exact = pd.read_csv(LOMO_EXACT_PER_MONKEY)[["monkey_id", "hit1_mean", "hit3_mean"]].rename(
        columns={"hit1_mean": "exact_top1", "hit3_mean": "exact_top3"}
    )
    group = pd.read_csv(LOMO_GROUP_PER_MONKEY)[["monkey_id", "group_hit1_mean", "group_hit3_mean"]].rename(
        columns={"group_hit1_mean": "group_top1", "group_hit3_mean": "group_top3"}
    )
    return network.merge(exact, on="monkey_id").merge(group, on="monkey_id")


def plot_lomo_per_monkey(per_monkey: pd.DataFrame) -> None:
    monkey_ids = per_monkey["monkey_id"].astype(str).tolist()
    x = np.arange(len(monkey_ids))
    width = 0.26
    fig, ax = plt.subplots(figsize=(13.4, 5.7), constrained_layout=True)
    bars_n = ax.bar(x - width, per_monkey["network_top3"], width, label="Network Top3", color="#0072B2")
    bars_e = ax.bar(x, per_monkey["exact_top3"], width, label="Exact Region Top3", color="#D55E00")
    bars_g = ax.bar(x + width, per_monkey["group_top3"], width, label="Region Group Top3", color="#009E73")
    ax.set_title("LOMO per held-out monkey: Top3 accuracy", fontweight="bold")
    ax.set_xticks(x, monkey_ids)
    ax.set_xlabel("Held-out monkey")
    style_accuracy_axis(ax)
    ax.legend(loc="upper right", ncols=3)
    for bars in [bars_n, bars_e, bars_g]:
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.012,
                f"{bar.get_height() * 100:.0f}%",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )
    fig.savefig(OUTDIR / "lomo_per_monkey_top3.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUTDIR / "lomo_per_monkey_top3.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_lomo_heatmap(per_monkey: pd.DataFrame) -> None:
    metrics = ["network_top1", "network_top3", "exact_top1", "exact_top3", "group_top1", "group_top3"]
    labels = ["Network Top1", "Network Top3", "Exact Top1", "Exact Top3", "Group Top1", "Group Top3"]
    data = per_monkey.set_index("monkey_id")[metrics].T.to_numpy()
    fig, ax = plt.subplots(figsize=(11.8, 5.6), constrained_layout=True)
    im = ax.imshow(data, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_title("LOMO per held-out monkey: hierarchical endpoint heatmap", fontweight="bold")
    ax.set_xticks(np.arange(len(per_monkey)), per_monkey["monkey_id"].astype(str).tolist())
    ax.set_yticks(np.arange(len(labels)), labels)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, pct(float(data[i, j])), ha="center", va="center", color="black", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, shrink=0.84)
    cbar.set_label("Accuracy")
    fig.savefig(OUTDIR / "lomo_per_monkey_heatmap.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUTDIR / "lomo_per_monkey_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def write_markdown(endpoint_df: pd.DataFrame, network_df: pd.DataFrame, per_monkey: pd.DataFrame) -> None:
    lines = [
        "# Bo2023 LOSO vs LOMO validation summary",
        "",
        "## Formal endpoint summary",
        "",
        "| Validation | Endpoint | n | Top1 | Top3 | Notes |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ]
    for row in endpoint_df.itertuples(index=False):
        lines.append(f"| {row.validation} | {row.endpoint} | {row.n} | {pct(row.top1)} | {pct(row.top3)} | {row.notes} |")
    lines.extend(
        [
            "",
            "## Network route comparison",
            "",
            "| Validation | Route | n | Top1 | Top3 |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in network_df.itertuples(index=False):
        lines.append(f"| {row.validation} | {row.route} | {row.n} | {pct(row.top1)} | {pct(row.top3)} |")
    lines.extend(
        [
            "",
            "## LOMO per held-out monkey",
            "",
            "| Monkey | n | Network Top1 | Network Top3 | Exact Top1 | Exact Top3 | Group Top1 | Group Top3 |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in per_monkey.itertuples(index=False):
        lines.append(
            f"| {row.monkey_id} | {row.n} | {pct(row.network_top1)} | {pct(row.network_top3)} "
            f"| {pct(row.exact_top1)} | {pct(row.exact_top3)} | {pct(row.group_top1)} | {pct(row.group_top3)} |"
        )
    lines.extend(
        [
            "",
            "## Figures",
            "",
            "- `loso_lomo_endpoint_summary.png`",
            "- `network_baseline_vs_pairwise_rescue.png`",
            "- `lomo_per_monkey_top3.png`",
            "- `lomo_per_monkey_heatmap.png`",
            "",
        ]
    )
    (OUTDIR / "bo2023_loso_lomo_validation_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    endpoint_df, network_df = build_summary_tables()
    per_monkey = build_per_monkey_table()

    endpoint_df.to_csv(OUTDIR / "loso_lomo_endpoint_summary.csv", index=False, encoding="utf-8-sig")
    network_df.to_csv(OUTDIR / "network_baseline_vs_pairwise_rescue.csv", index=False, encoding="utf-8-sig")
    per_monkey.to_csv(OUTDIR / "lomo_per_monkey_summary.csv", index=False, encoding="utf-8-sig")

    plot_endpoint_summary(endpoint_df)
    plot_network_rescue(network_df)
    plot_lomo_per_monkey(per_monkey)
    plot_lomo_heatmap(per_monkey)
    write_markdown(endpoint_df, network_df, per_monkey)
    print(f"Outputs written to: {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
