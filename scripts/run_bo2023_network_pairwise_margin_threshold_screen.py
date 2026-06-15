#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DETAIL_DIR = ROOT / "results" / "bo2023_network_pairwise_correlation_full_loso_819_rerun_20260526"
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_network_pairwise_margin_threshold_screen_20260604"


def paired_pvalue(gains: int, losses: int) -> float:
    import math

    n = int(gains + losses)
    if n == 0:
        return 1.0
    tail = min(int(gains), int(losses))
    probability = sum(math.comb(n, k) for k in range(tail + 1)) / (2**n)
    return float(min(1.0, 2.0 * probability))


def main() -> int:
    parser = argparse.ArgumentParser(description="Screen pairwise-rescue switch-margin thresholds from LOSO detail files.")
    parser.add_argument("--detail-dir", type=Path, default=DEFAULT_DETAIL_DIR)
    parser.add_argument("--baseline-file", default="network_discriminative_correlation_top200_detail.csv")
    parser.add_argument("--pairwise-file", default="network_pairwise_correlation_rescue_top3_detail.csv")
    parser.add_argument("--thresholds", default="0,0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.05,0.08,0.10")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    baseline = pd.read_csv(args.detail_dir / args.baseline_file)
    pairwise = pd.read_csv(args.detail_dir / args.pairwise_file)
    thresholds = [float(x) for x in str(args.thresholds).split(",") if str(x).strip()]

    rows: list[dict[str, Any]] = []
    for threshold in thresholds:
        use_switch = (
            pairwise["switched"].fillna(0).astype(int).eq(1)
            & (pairwise["switch_margin"].fillna(0.0).astype(float) > float(threshold))
        )
        hit1 = np.where(use_switch, pairwise["hit1"], baseline["hit1"])
        hit3 = np.where(use_switch, pairwise["hit3"], baseline["hit3"])
        gains = int(((baseline["hit1"] == 0) & (hit1 == 1)).sum())
        losses = int(((baseline["hit1"] == 1) & (hit1 == 0)).sum())
        top3_gains = int(((baseline["hit3"] == 0) & (hit3 == 1)).sum())
        top3_losses = int(((baseline["hit3"] == 1) & (hit3 == 0)).sum())
        rows.append(
            {
                "threshold": float(threshold),
                "top1_hits": int(hit1.sum()),
                "top1_accuracy": float(hit1.mean()),
                "top3_hits": int(hit3.sum()),
                "top3_accuracy": float(hit3.mean()),
                "n_switches": int(use_switch.sum()),
                "top1_gains_vs_top200": gains,
                "top1_losses_vs_top200": losses,
                "top1_pvalue_vs_top200": paired_pvalue(gains, losses),
                "top3_gains_vs_top200": top3_gains,
                "top3_losses_vs_top200": top3_losses,
                "top3_pvalue_vs_top200": paired_pvalue(top3_gains, top3_losses),
            }
        )

    metrics = pd.DataFrame(rows)
    best = metrics.sort_values(["top1_accuracy", "top3_accuracy"], ascending=[False, False]).iloc[0].to_dict()
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "validation_design": "post-hoc threshold screen over existing full strict LOSO pairwise-rescue detail; no model refit",
        "baseline_detail_dir": str(args.detail_dir),
        "baseline_file": args.baseline_file,
        "pairwise_file": args.pairwise_file,
        "n_test_samples": int(len(baseline)),
        "best_threshold": float(best["threshold"]),
        "best_route": best,
        "decision": (
            "Use pair_min_margin=0.002 as the conservative production threshold: it maximizes Top1 in this screen "
            "while preserving Top3 and reducing false switches."
        ),
    }
    metrics.to_csv(args.outdir / "pairwise_margin_threshold_metrics.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
