#!/usr/bin/env python
"""Donor-clustered inference for the locked Bo2023 LOSO/LOMO route.

This analysis does not retrain, tune, or alter the tracing route.  It reads the
locked-route detail tables, joins the LOSO records to MonkeyID, and reports
sample-weighted, donor-macro, donor-bootstrap, and cluster-robust results.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core.external_inputs import resolve_alias  # noqa: E402
DEFAULT_INPUT = ROOT / "reports" / "p0_donor_cluster_20260711"
DEFAULT_OUTDIR = DEFAULT_INPUT / "donor_cluster_inference"
DEFAULT_SAMPLE_INFO = resolve_alias("bo2023_sample_metadata")
SEED = 20260711
N_BOOTSTRAP = 10_000


ENDPOINTS = (
    ("Network Top3", "hybrid_formal_loso_network_detail.csv", "formal_lomo_network_detail.csv", "hit3"),
    ("Resolution-group Top3", "hybrid_formal_loso_resolution_group_detail.csv", "formal_lomo_resolution_group_detail.csv", "group_hit3"),
    ("Exact-region Top3", "hybrid_formal_loso_exact_region_detail.csv", "formal_lomo_exact_region_detail.csv", "hit3"),
)


def monkey_map(sample_info: Path) -> pd.DataFrame:
    info = pd.read_excel(sample_info, sheet_name="mfas5_819samples_phenSet4", usecols=["No.", "MonkeyID"])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["monkey_id"] = info["MonkeyID"].astype(str).str.strip()
    return info[["sample_id", "monkey_id"]].drop_duplicates("sample_id")


def add_monkey_id(frame: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    if "monkey_id" in frame.columns:
        frame = frame.drop(columns=["monkey_id"])
    out = frame.merge(mapping, on="sample_id", how="left", validate="one_to_one")
    if out["monkey_id"].isna().any():
        missing = out.loc[out["monkey_id"].isna(), "sample_id"].head().tolist()
        raise ValueError(f"Missing MonkeyID for samples: {missing}")
    return out


def t_critical_975(df: int) -> float:
    try:
        from scipy.stats import t

        return float(t.ppf(0.975, df))
    except Exception:
        # df=8 in this analysis; conservative fallback if SciPy is unavailable.
        return 2.306


def two_sided_t_pvalue(t_value: float, df: int) -> float:
    try:
        from scipy.stats import t

        return float(2 * t.sf(abs(t_value), df))
    except Exception:
        return float("nan")


def cluster_robust_mean(values: np.ndarray, clusters: np.ndarray, bounds: tuple[float, float] | None = (0.0, 1.0)) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    clusters = np.asarray(clusters, dtype=str)
    n = len(values)
    labels = np.unique(clusters)
    g = len(labels)
    estimate = float(values.mean())
    scores = np.array([values[clusters == label].sum() - (clusters == label).sum() * estimate for label in labels])
    variance = (g / (g - 1)) * float(np.square(scores).sum()) / (n * n) if g > 1 else float("nan")
    se = math.sqrt(max(variance, 0.0))
    crit = t_critical_975(g - 1)
    low, high = estimate - crit * se, estimate + crit * se
    if bounds is not None:
        low, high = max(bounds[0], low), min(bounds[1], high)
    return {
        "cluster_n": int(g),
        "cluster_robust_se": se,
        "cluster_robust_ci_low": float(low),
        "cluster_robust_ci_high": float(high),
    }


def donor_bootstrap(values: np.ndarray, clusters: np.ndarray, rng: np.random.Generator, n_bootstrap: int) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    clusters = np.asarray(clusters, dtype=str)
    labels = np.unique(clusters)
    per_cluster = [values[clusters == label] for label in labels]
    per_cluster_means = np.array([x.mean() for x in per_cluster])
    weighted = np.empty(n_bootstrap, dtype=float)
    macro = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        draw = rng.integers(0, len(labels), size=len(labels))
        weighted[i] = np.concatenate([per_cluster[j] for j in draw]).mean()
        macro[i] = per_cluster_means[draw].mean()
    return {
        "donor_bootstrap_weighted_ci_low": float(np.quantile(weighted, 0.025)),
        "donor_bootstrap_weighted_ci_high": float(np.quantile(weighted, 0.975)),
        "donor_bootstrap_macro_ci_low": float(np.quantile(macro, 0.025)),
        "donor_bootstrap_macro_ci_high": float(np.quantile(macro, 0.975)),
    }


def exact_mcnemar(x: np.ndarray, y: np.ndarray) -> dict[str, float | int]:
    x = np.asarray(x, dtype=int)
    y = np.asarray(y, dtype=int)
    only_x = int(((x == 1) & (y == 0)).sum())
    only_y = int(((x == 0) & (y == 1)).sum())
    discordant = only_x + only_y
    if discordant == 0:
        p = 1.0
    else:
        from scipy.stats import binomtest

        p = float(binomtest(min(only_x, only_y), discordant, 0.5, alternative="two-sided").pvalue)
    return {"loso_only_hits": only_x, "lomo_only_hits": only_y, "discordant_pairs": discordant, "mcnemar_exact_p": p}


def bonferroni(p: np.ndarray) -> np.ndarray:
    return np.minimum(1.0, p * len(p))


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    adjusted = np.empty_like(p, dtype=float)
    running = 1.0
    m = len(p)
    for rank in range(m, 0, -1):
        idx = order[rank - 1]
        running = min(running, p[idx] * m / rank)
        adjusted[idx] = running
    return np.minimum(1.0, adjusted)


def endpoint_summary(name: str, loso: pd.DataFrame, lomo: pd.DataFrame, hit_col: str, rng: np.random.Generator, n_bootstrap: int) -> tuple[dict, pd.DataFrame]:
    merged = loso[["sample_id", "monkey_id", hit_col]].merge(
        lomo[["sample_id", "monkey_id", hit_col]], on="sample_id", suffixes=("_loso", "_lomo"), validate="one_to_one"
    )
    if not (merged["monkey_id_loso"] == merged["monkey_id_lomo"]).all():
        raise ValueError(f"Donor mismatch after LOSO/LOMO merge for {name}")
    merged = merged.rename(columns={"monkey_id_loso": "monkey_id"}).drop(columns=["monkey_id_lomo"])
    x = merged[f"{hit_col}_loso"].to_numpy(dtype=float)
    y = merged[f"{hit_col}_lomo"].to_numpy(dtype=float)
    clusters = merged["monkey_id"].to_numpy()
    diff = x - y
    donor = merged.groupby("monkey_id", sort=True).agg(
        n=("sample_id", "size"),
        loso_accuracy=(f"{hit_col}_loso", "mean"),
        lomo_accuracy=(f"{hit_col}_lomo", "mean"),
    ).reset_index()
    donor["loso_minus_lomo"] = donor["loso_accuracy"] - donor["lomo_accuracy"]
    donor.insert(0, "endpoint", name)

    summary = {
        "endpoint": name,
        "paired_n": int(len(merged)),
        "n_donors": int(merged["monkey_id"].nunique()),
        "loso_sample_weighted_accuracy": float(x.mean()),
        "lomo_sample_weighted_accuracy": float(y.mean()),
        "loso_donor_macro_accuracy": float(donor["loso_accuracy"].mean()),
        "lomo_donor_macro_accuracy": float(donor["lomo_accuracy"].mean()),
        "sample_weighted_difference_loso_minus_lomo": float(diff.mean()),
        "donor_macro_difference_loso_minus_lomo": float(donor["loso_minus_lomo"].mean()),
    }
    for prefix, values in (("loso", x), ("lomo", y), ("difference", diff)):
        bounds = None if prefix == "difference" else (0.0, 1.0)
        summary.update({f"{prefix}_{k}": v for k, v in cluster_robust_mean(values, clusters, bounds=bounds).items()})
        summary.update({f"{prefix}_{k}": v for k, v in donor_bootstrap(values, clusters, rng, n_bootstrap).items()})
    diff_cr = cluster_robust_mean(diff, clusters, bounds=None)
    if diff_cr["cluster_robust_se"]:
        t_value = float(diff.mean() / diff_cr["cluster_robust_se"])
        summary["difference_cluster_robust_t"] = t_value
        summary["difference_cluster_robust_p"] = two_sided_t_pvalue(t_value, int(diff_cr["cluster_n"]) - 1)
    else:
        summary["difference_cluster_robust_t"] = float("nan")
        summary["difference_cluster_robust_p"] = float("nan")
    summary.update(exact_mcnemar(x, y))
    return summary, donor


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    mapping = monkey_map(args.sample_info)
    rng = np.random.default_rng(args.seed)
    summaries, donor_tables = [], []
    for name, loso_file, lomo_file, hit_col in ENDPOINTS:
        loso = pd.read_csv(args.input_root / "formal_loso_locked_route" / loso_file)
        lomo = pd.read_csv(args.input_root / "formal_lomo_locked_route" / lomo_file)
        for frame_name, frame in (("LOSO", loso), ("LOMO", lomo)):
            if "route_family" in frame.columns:
                selected = frame["route_family"].eq("hybrid_projected_network_logcpm_exact")
                if not selected.any():
                    raise ValueError(f"Locked hybrid route is absent from {frame_name} {name} detail table")
                if frame_name == "LOSO":
                    loso = frame.loc[selected].copy()
                else:
                    lomo = frame.loc[selected].copy()
        loso = add_monkey_id(loso, mapping)
        lomo = add_monkey_id(lomo, mapping)
        summary, donor = endpoint_summary(name, loso, lomo, hit_col, rng, args.n_bootstrap)
        summaries.append(summary)
        donor_tables.append(donor)
    results = pd.DataFrame(summaries)
    for col in ("mcnemar_exact_p", "difference_cluster_robust_p"):
        p = results[col].to_numpy(dtype=float)
        results[f"{col}_bonferroni_3_endpoints"] = bonferroni(p)
        results[f"{col}_bh_fdr_3_endpoints"] = benjamini_hochberg(p)
    results.to_csv(args.outdir / "donor_clustered_loso_lomo_summary.csv", index=False)
    pd.concat(donor_tables, ignore_index=True).to_csv(args.outdir / "donor_level_endpoint_metrics.csv", index=False)
    methods = {
        "locked_route": "hybrid_projected_network_logcpm_exact",
        "resampling_unit": "monkey/donor; all samples from a selected donor retained together",
        "n_donors": 9,
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "cluster_robust_ci": "sandwich variance for a clustered mean with small-cluster correction G/(G-1); t critical value with G-1=8 degrees of freedom",
        "multiple_testing_family": "three prespecified hierarchical Top3 LOSO-vs-LOMO comparisons: Network, resolution group, exact region",
        "multiple_testing_methods": ["Bonferroni", "Benjamini-Hochberg FDR"],
        "mcnemar_role": "paired sample-level sensitivity analysis; it does not replace donor-clustered inference",
    }
    (args.outdir / "STATISTICAL_METHODS.json").write_text(json.dumps(methods, indent=2), encoding="utf-8")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
