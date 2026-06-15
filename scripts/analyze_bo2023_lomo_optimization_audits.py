#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULT_DIR = ROOT / "results" / "bo2023_leave_one_monkey_out_formal_route_20260601"


def _rate_frame(df: pd.DataFrame, group_cols: list[str], metric_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for keys, frame in df.groupby(group_cols, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {col: key for col, key in zip(group_cols, keys)}
        row["n"] = int(len(frame))
        for col in metric_cols:
            row[f"{col}_hits"] = int(frame[col].sum())
            row[f"{col}_rate"] = float(frame[col].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def network_confusion(pair: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    top1 = (
        pair.groupby(["monkey_id", "label", "pred_top1"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["monkey_id", "n"], ascending=[True, False])
    )
    misses = pair[pair["hit3"] == 0].copy()
    top3_misses = (
        misses.groupby(["monkey_id", "label", "pred_top1", "pred_top2", "pred_top3"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["monkey_id", "n"], ascending=[True, False])
    )
    return top1, top3_misses


def pairwise_gains_losses(base: pd.DataFrame, pair: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = [
        "sample_id",
        "monkey_id",
        "label",
        "pred_top1",
        "pred_top2",
        "pred_top3",
        "true_rank",
        "hit1",
        "hit3",
    ]
    merged = base[cols].merge(
        pair[
            cols
            + [
                "switched",
                "switch_pair",
                "switch_margin",
            ]
        ],
        on=["sample_id", "monkey_id", "label"],
        suffixes=("_baseline", "_pairwise"),
    )
    merged["top1_change"] = merged["hit1_pairwise"] - merged["hit1_baseline"]
    merged["top3_change"] = merged["hit3_pairwise"] - merged["hit3_baseline"]
    merged["change_type"] = "unchanged"
    merged.loc[merged["top1_change"].gt(0), "change_type"] = "top1_gain"
    merged.loc[merged["top1_change"].lt(0), "change_type"] = "top1_loss"
    by_monkey = _rate_frame(merged, ["monkey_id", "change_type"], ["switched"])
    by_pair = (
        merged[merged["switched"].eq(1)]
        .groupby(["switch_pair", "change_type"], dropna=False)
        .agg(
            n=("sample_id", "size"),
            mean_margin=("switch_margin", "mean"),
            top1_delta=("top1_change", "sum"),
            top3_delta=("top3_change", "sum"),
        )
        .reset_index()
        .sort_values(["top1_delta", "n"], ascending=[False, False])
    )
    return merged, pd.concat([by_monkey.assign(summary_level="monkey"), by_pair.assign(summary_level="pair")], ignore_index=True)


def group_stability(group: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    truth_rows: list[dict[str, Any]] = []
    for region, frame in group.groupby("label", sort=True):
        groups = frame["true_resolution_group"].fillna("").astype(str)
        members = frame["truth_group_members"].fillna("").astype(str)
        top_groups = groups.value_counts()
        truth_rows.append(
            {
                "region_id": region,
                "n": int(len(frame)),
                "n_monkeys": int(frame["monkey_id"].nunique()),
                "n_distinct_true_groups": int(groups.nunique()),
                "dominant_true_group": str(top_groups.index[0]) if len(top_groups) else "",
                "dominant_true_group_n": int(top_groups.iloc[0]) if len(top_groups) else 0,
                "dominant_true_group_fraction": float(top_groups.iloc[0] / len(frame)) if len(top_groups) else 0.0,
                "n_distinct_truth_member_sets": int(members.nunique()),
                "group_top1_rate": float(frame["group_hit1"].mean()),
                "group_top3_rate": float(frame["group_hit3"].mean()),
                "exact_top1_rate": float(frame["hit1"].mean()),
                "exact_top3_rate": float(frame["hit3"].mean()),
            }
        )
    stability = pd.DataFrame(truth_rows).sort_values(
        ["dominant_true_group_fraction", "n"], ascending=[True, False]
    )
    failure = (
        group[group["group_hit3"].eq(0)]
        .groupby(["monkey_id", "label", "true_resolution_group", "pred_group_top1"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["monkey_id", "n"], ascending=[True, False])
    )
    bridge = (
        group[(group["hit3"].eq(0)) & (group["group_hit3"].eq(1))]
        .groupby(["monkey_id", "label", "true_resolution_group"], dropna=False)
        .size()
        .reset_index(name="n_exact_miss_group_hit")
        .sort_values(["monkey_id", "n_exact_miss_group_hit"], ascending=[True, False])
    )
    return stability, failure, bridge


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate optimization audits from Bo2023 LOMO formal-route results.")
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    args = parser.parse_args()

    base = pd.read_csv(args.result_dir / "network_discriminative_correlation_top200_detail.csv")
    pair = pd.read_csv(args.result_dir / "network_pairwise_correlation_rescue_top3_detail.csv")
    exact = pd.read_csv(args.result_dir / "top3_beam_local_top50_top100_zfusion_w0p25_detail.csv")
    group = pd.read_csv(args.result_dir / "top3_network_beam_local_region_candidates_detail.csv")

    network_top1_confusion, network_top3_misses = network_confusion(pair)
    pairwise_detail, pairwise_summary = pairwise_gains_losses(base, pair)
    group_stability_table, group_top3_failures, exact_miss_group_hit = group_stability(group)

    network_top1_confusion.to_csv(args.result_dir / "audit_network_top1_confusion_by_monkey.csv", index=False, encoding="utf-8-sig")
    network_top3_misses.to_csv(args.result_dir / "audit_network_top3_misses_by_monkey.csv", index=False, encoding="utf-8-sig")
    pairwise_detail.to_csv(args.result_dir / "audit_network_pairwise_gains_losses_detail.csv", index=False, encoding="utf-8-sig")
    pairwise_summary.to_csv(args.result_dir / "audit_network_pairwise_gains_losses_summary.csv", index=False, encoding="utf-8-sig")
    group_stability_table.to_csv(args.result_dir / "audit_resolution_group_stability_by_region.csv", index=False, encoding="utf-8-sig")
    group_top3_failures.to_csv(args.result_dir / "audit_resolution_group_top3_failures.csv", index=False, encoding="utf-8-sig")
    exact_miss_group_hit.to_csv(args.result_dir / "audit_exact_miss_group_hit_cases.csv", index=False, encoding="utf-8-sig")

    summary = {
        "network_top3_misses": int((pair["hit3"] == 0).sum()),
        "network_top1_pairwise_gains": int(((base["hit1"] == 0) & (pair["hit1"] == 1)).sum()),
        "network_top1_pairwise_losses": int(((base["hit1"] == 1) & (pair["hit1"] == 0)).sum()),
        "network_top3_pairwise_gains": int(((base["hit3"] == 0) & (pair["hit3"] == 1)).sum()),
        "network_top3_pairwise_losses": int(((base["hit3"] == 1) & (pair["hit3"] == 0)).sum()),
        "group_unstable_regions": int((group_stability_table["dominant_true_group_fraction"] < 1.0).sum()),
        "group_top3_failures": int((group["group_hit3"] == 0).sum()),
        "exact_miss_group_hit_samples": int(((group["hit3"] == 0) & (group["group_hit3"] == 1)).sum()),
    }
    (args.result_dir / "audit_optimization_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Audit outputs written to: {args.result_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
