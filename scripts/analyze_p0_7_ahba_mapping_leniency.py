#!/usr/bin/env python
"""P0-7 audit of AHBA mapped-label multiplicity and strict/lenient scoring.

Predictions are read from the locked, precomputed hybrid AHBA route.  This
script does not change the model or anatomical mapping; it only makes the
number and granularity of allowed labels explicit and rescoring transparent.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DETAIL = ROOT / "reports" / "p0_7_ahba_mapping_20260711" / "formal_three_tier" / "ahba_formal_three_tier_sample_detail.csv"
DEFAULT_GROUPS = ROOT / "data" / "models" / "bo2023_region_resolution_groups.json"
DEFAULT_OUTDIR = ROOT / "reports" / "p0_7_ahba_mapping_20260711" / "mapping_leniency_audit"
HYBRID = "hybrid_projected_network_logcpm_exact"


def split_pipe(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [item.strip() for item in str(value).split("|") if item.strip()]


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    z = 1.959963984540054
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return centre - half, centre + half


def summarize(frame: pd.DataFrame, endpoint: str, definition: str, hit_cols: tuple[str, str]) -> list[dict]:
    rows = []
    for name, col in zip(("Top1", "Top3"), hit_cols):
        values = pd.to_numeric(frame[col], errors="coerce").dropna().astype(int)
        hits, n = int(values.sum()), int(len(values))
        lo, hi = wilson(hits, n)
        rows.append({
            "endpoint": endpoint,
            "mapping_definition": definition,
            "metric": name,
            "n": n,
            "hits": hits,
            "accuracy": hits / n if n else float("nan"),
            "wilson95_low": lo,
            "wilson95_high": hi,
        })
    return rows


def plot_metrics(metrics: pd.DataFrame, outpath: Path) -> None:
    top3 = metrics[metrics["metric"].eq("Top3")].copy()
    definitions = ["lenient", "strict_unique_mapping"]
    endpoints = ["Network", "Resolution group", "Exact region"]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(endpoints))
    width = 0.34
    colors = ["#2d6ca2", "#c96f30"]
    for i, definition in enumerate(definitions):
        subset = top3[top3["mapping_definition"].eq(definition)].set_index("endpoint").reindex(endpoints)
        values = subset["accuracy"].to_numpy(dtype=float)
        bars = ax.bar(x + (i - 0.5) * width, values, width, label=definition.replace("_", " "), color=colors[i])
        ax.bar_label(bars, labels=[f"{100*v:.1f}%\nn={int(n)}" for v, n in zip(values, subset["n"])], padding=3, fontsize=8)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Top3 accuracy")
    ax.set_xticks(x, endpoints)
    ax.set_title("AHBA mapped-label leniency sensitivity (locked predictions)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detail", type=Path, default=DEFAULT_DETAIL)
    parser.add_argument("--groups", type=Path, default=DEFAULT_GROUPS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    detail = pd.read_csv(args.detail)
    detail = detail[detail["route"].eq(HYBRID)].copy()
    if detail.empty:
        raise ValueError(f"No {HYBRID} records in {args.detail}")
    groups = json.loads(args.groups.read_text(encoding="utf-8"))["entries"]
    region_to_group = {str(item["region_id"]): str(item["resolution_group"]) for item in groups.values()}

    detail["allowed_network_list"] = detail["allowed_bo2023_networks"].map(split_pipe)
    region_key_col = "allowed_bo2023_region_keys" if "allowed_bo2023_region_keys" in detail.columns else "allowed_bo2023_regions"
    detail["allowed_region_list"] = detail[region_key_col].map(split_pipe)
    detail["allowed_static_group_list"] = detail["allowed_region_list"].map(
        lambda items: sorted(
            {
                region_to_group[item.rsplit("::", 1)[-1]]
                for item in items
                if item.rsplit("::", 1)[-1] in region_to_group
            }
        )
    )
    detail["n_allowed_networks"] = detail["allowed_network_list"].map(len)
    detail["n_allowed_regions"] = detail["allowed_region_list"].map(len)
    detail["n_allowed_static_groups"] = detail["allowed_static_group_list"].map(len)
    detail["network_unique_mapping"] = detail["n_allowed_networks"].eq(1)
    detail["group_unique_mapping"] = detail["n_allowed_static_groups"].eq(1)
    detail["exact_unique_mapping"] = detail["n_allowed_regions"].eq(1)
    detail["supported_for_accuracy"] = detail["supported_for_accuracy"].astype(bool)
    detail["exact_mapped"] = detail["n_allowed_regions"].gt(0)

    audit_cols = [
        "sample_id", "ahba_donor", "public_label", "public_main_structure", "public_sub_structure", "supported_for_accuracy",
        "accuracy_level", "allowed_bo2023_networks", "allowed_bo2023_regions", region_key_col, "allowed_network_list", "allowed_region_list",
        "allowed_static_group_list", "n_allowed_networks", "n_allowed_regions", "n_allowed_static_groups", "network_unique_mapping",
        "group_unique_mapping", "exact_unique_mapping", "network_top1", "network_top2", "network_top3", "network_top1_hit",
        "network_top3_hit", "group_top1_hit", "group_top3_hit", "region_top1_exact_hit", "region_top3_exact_hit",
    ]
    detail[audit_cols].to_csv(args.outdir / "ahba_mapping_multiplicity_per_sample.csv", index=False)

    granular_rows = []
    for scope, frame in (("all", detail), ("supported", detail[detail["supported_for_accuracy"]]), ("exact_mapped", detail[detail["exact_mapped"]])):
        for metric in ("n_allowed_networks", "n_allowed_regions", "n_allowed_static_groups"):
            values = frame[metric].astype(float)
            granular_rows.append({
                "scope": scope, "quantity": metric, "n_samples": int(len(values)), "mean": float(values.mean()),
                "median": float(values.median()), "iqr_low": float(values.quantile(0.25)), "iqr_high": float(values.quantile(0.75)),
                "single_label_fraction": float((values == 1).mean()), "multi_label_fraction": float((values > 1).mean()),
            })
    granularity = pd.DataFrame(granular_rows)
    granularity.to_csv(args.outdir / "ahba_mapping_granularity_distribution.csv", index=False)

    metrics: list[dict] = []
    # Lenient = the original mapped-label scoring: a hit against any allowed mapping.
    metrics += summarize(detail[detail["supported_for_accuracy"]], "Network", "lenient", ("network_top1_hit", "network_top3_hit"))
    metrics += summarize(detail[detail["supported_for_accuracy"] & detail["network_unique_mapping"]], "Network", "strict_unique_mapping", ("network_top1_hit", "network_top3_hit"))
    metrics += summarize(detail[detail["exact_mapped"]], "Resolution group", "lenient", ("group_top1_hit", "group_top3_hit"))
    metrics += summarize(detail[detail["exact_mapped"] & detail["group_unique_mapping"]], "Resolution group", "strict_unique_mapping", ("group_top1_hit", "group_top3_hit"))
    metrics += summarize(detail[detail["exact_mapped"]], "Exact region", "lenient", ("region_top1_exact_hit", "region_top3_exact_hit"))
    metrics += summarize(detail[detail["exact_mapped"] & detail["exact_unique_mapping"]], "Exact region", "strict_unique_mapping", ("region_top1_exact_hit", "region_top3_exact_hit"))
    metric_df = pd.DataFrame(metrics)
    metric_df.to_csv(args.outdir / "ahba_strict_lenient_accuracy.csv", index=False)
    plot_metrics(metric_df, args.outdir / "ahba_strict_lenient_top3.png")

    methods = {
        "route": HYBRID,
        "lenient_definition": "Original mapped-label scoring: prediction is correct if any predicted label matches any allowed Bo2023 label.",
        "network_strict_definition": "Supported AHBA samples with exactly one allowed Bo2023 Network.",
        "group_strict_definition": "Exact-mapped AHBA samples whose allowed Bo2023 regions resolve to exactly one static Bo2023 resolution group in data/models/bo2023_region_resolution_groups.json.",
        "exact_strict_definition": "Exact-mapped AHBA samples with exactly one allowed Bo2023 region.",
        "important_caveat": "Strict and lenient results use the identical locked predictions. Strict subsets answer a narrower, lower-ambiguity label question; they are not a replacement for lenient mapped-label transfer metrics.",
    }
    (args.outdir / "P0_7_METHODS.json").write_text(json.dumps(methods, indent=2), encoding="utf-8")
    print(metric_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
