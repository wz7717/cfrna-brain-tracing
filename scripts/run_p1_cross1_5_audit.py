#!/usr/bin/env python
"""Reproducible calculations supporting P1-CROSS1--5 manuscript repair."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def macro_metrics(path: Path) -> pd.DataFrame:
    payload = json.loads(path.read_text(encoding="utf-8"))
    frame = pd.DataFrame(payload["data"])
    frame = frame[frame["endpoint"].isin(
        ["LOSO_Network", "LOMO_Network", "LOSO_Exact", "LOMO_Exact"]
    )].copy()
    for column in ("n", "precision", "recall", "f1"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    rows = []
    for endpoint, part in frame.groupby("endpoint", sort=True):
        support = part["n"].sum()
        # In single-label multiclass evaluation micro-F1 equals accuracy.
        tp = (part["recall"] * part["n"]).sum()
        rows.append(
            {
                "endpoint": endpoint,
                "n_classes": int(len(part)),
                "n_samples": int(round(support)),
                "macro_f1": float(part["f1"].mean()),
                # Descriptive spread of the same evaluable class-level F1
                # values used for macro-F1.  Use sample SD (ddof=1) and make
                # the denominator explicit in the manuscript rather than
                # presenting a second, incompatible uncertainty quantity.
                "sd_class_f1": float(part["f1"].std(ddof=1)),
                "median_class_f1": float(part["f1"].median()),
                "weighted_f1": float(np.average(part["f1"], weights=part["n"])),
                "micro_f1": float(tp / support),
                "n_zero_f1_classes": int((part["f1"] == 0).sum()),
                "fraction_zero_f1_classes": float((part["f1"] == 0).mean()),
                "conditional_macro_f1_nonzero": float(
                    part.loc[part["f1"] > 0, "f1"].mean()
                ),
                "conditional_median_f1_nonzero": float(
                    part.loc[part["f1"] > 0, "f1"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--macro-class-data", type=Path, required=True)
    parser.add_argument("--lambda-friedman", type=Path, required=True)
    parser.add_argument("--lambda-sensitivity", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    macro = macro_metrics(args.macro_class_data)
    macro.to_csv(args.outdir / "cross3_f1_distribution_summary.csv", index=False)

    friedman = pd.read_csv(args.lambda_friedman)
    sensitivity = pd.read_csv(args.lambda_sensitivity)
    candidate = sensitivity[sensitivity["lambda"].isin([0.25, 0.50, 0.75])].copy()
    exact_top3_range_pp = (
        candidate["exact_top3"].max() - candidate["exact_top3"].min()
    ) * 100
    cross1 = {
        "candidate_lambdas": candidate[
            ["lambda", "exact_top1", "exact_top3"]
        ].to_dict(orient="records"),
        "exact_top3_range_percentage_points": float(exact_top3_range_pp),
        "friedman_top3": friedman.loc[
            friedman["lambda"].astype(str).eq("Friedman_test")
        ].iloc[0].to_dict(),
        "interpretation": (
            "p=0.764 indicates insufficient evidence to detect a lambda effect; "
            "it does not identify lambda=0.25 as statistically optimal. Lambda=0.25 "
            "is retained as a prespecified stability-specificity operating choice."
        ),
    }
    (args.outdir / "cross1_lambda_interpretation.json").write_text(
        json.dumps(cross1, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    logcpm = {
        "formula": "x_gs = ln(1 + 1e6 * c_gs / sum_h c_hs)",
        "normalization": "within-sample total-count library-size scaling only",
        "pseudocount": "1 CPM inside ln(1+CPM)",
        "edgeR_prior_count": (
            "No edgeR prior.count parameter is used; no TMM, quantile, batch, "
            "or target-cohort distribution normalization is fitted."
        ),
        "reranking": (
            "Within regions admitted by the fixed Network Top3 candidate filter, "
            "Pearson correlations are computed on fold-local discriminative gene "
            "orders. Exact score = 0.25*z(corr Top50) + 0.75*z(corr Top100); "
            "resolution-group ranking uses the fold-local panel correlation."
        ),
    }
    (args.outdir / "cross5_logcpm_reranking_definition.json").write_text(
        json.dumps(logcpm, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(macro.to_string(index=False))
    print(json.dumps(cross1, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
