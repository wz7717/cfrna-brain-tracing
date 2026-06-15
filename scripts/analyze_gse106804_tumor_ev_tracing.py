#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.network_tracing import trace_network_expression  # noqa: E402
from scripts.analyze_gse228512_ev_rna_tracing import (  # noqa: E402
    ADAPTED_MODEL,
    BASELINE_MODEL,
    NETWORK_TO_PRIMARY_BROAD,
    add_route,
    bh_adjust,
    model_scores,
    stratified_permutation_p,
)


DATA_DIR = ROOT / "data" / "external_validation" / "GSE106804"
OUTDIR = ROOT / "results" / "gse106804_tumor_ev_tracing_20260613"
GSE228_RESULTS = ROOT / "results" / "gse228512_ev_rna_tracing_20260613"


def load_data(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(
        data_dir / "GSE106804_Gene_counts.txt.gz", sep="\t", index_col="ORF"
    )
    counts.index = counts.index.astype(str).str.strip()
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = counts.groupby(level=0, sort=True).sum()

    raw = gzip.open(
        data_dir / "GSE106804_family.soft.gz",
        "rt",
        encoding="utf-8",
        errors="replace",
    ).read()
    rows = []
    for block in raw.split("^SAMPLE = ")[1:]:
        accession = re.search(r"!Sample_geo_accession = (.*)", block).group(1)
        title = re.search(r"!Sample_title = (.*)", block).group(1)
        characteristics = {}
        for item in re.findall(r"!Sample_characteristics_ch1 = (.*)", block):
            if ":" in item:
                key, value = item.split(":", 1)
                characteristics[key.strip().lower()] = value.strip()
        rows.append(
            {
                "sample_id": title,
                "geo_accession": accession,
                "group": (
                    "GBM"
                    if characteristics.get("group") == "Glioblastoma Multiforme patients"
                    else "Healthy"
                ),
                "biofluid": characteristics.get("biofluid", ""),
                "capture": characteristics.get("cell type", ""),
            }
        )
    metadata = pd.DataFrame(rows)
    if set(metadata["sample_id"]) != set(counts.columns):
        raise ValueError("Matrix and SOFT sample sets differ")
    metadata = (
        metadata.set_index("sample_id")
        .loc[counts.columns]
        .reset_index()
        .rename(columns={"index": "sample_id"})
    )
    library_sizes = counts.sum(axis=0)
    cpm = counts.divide(library_sizes.where(library_sizes > 0), axis=1) * 1_000_000.0
    return cpm, metadata


def production_predictions(cpm: pd.DataFrame, output: pd.DataFrame) -> pd.DataFrame:
    top1, top3, switched = [], [], []
    expression = np.log1p(cpm)
    for sample in expression.columns:
        frame = pd.DataFrame(
            {"gene_symbol": expression.index, "tpm_value": expression[sample].to_numpy()}
        )
        traced = trace_network_expression(frame, min_overlap_fraction=0.50)
        rows = traced["results"]
        top1.append(rows[0]["network_id"])
        top3.append(" | ".join(row["network_id"] for row in rows[:3]))
        switched.append(
            bool(traced["meta"].get("pairwise_rescue", {}).get("switched", False))
        )
    output["baseline_production_top1_network"] = top1
    output["baseline_production_top3_network"] = top3
    output["baseline_production_top1_broad"] = output[
        "baseline_production_top1_network"
    ].map(NETWORK_TO_PRIMARY_BROAD)
    output["baseline_production_top3_broad"] = output[
        "baseline_production_top3_network"
    ].map(
        lambda value: " | ".join(
            dict.fromkeys(
                NETWORK_TO_PRIMARY_BROAD[item.strip()]
                for item in str(value).split("|")
                if item.strip()
            )
        )
    )
    output["baseline_production_pairwise_switched"] = switched
    return output


def continuous_tests(
    output: pd.DataFrame, routes: list[str], iterations: int
) -> pd.DataFrame:
    rows = []
    labels = output["group"].to_numpy()
    strata = output["biofluid"].to_numpy()
    for route in routes:
        for metric in (
            "max_score",
            "raw_margin",
            "top1_probability",
            "probability_margin",
            "normalized_entropy",
        ):
            column = f"{route}_{metric}"
            gbm = output.loc[output.group.eq("GBM"), column].to_numpy(dtype=float)
            healthy = output.loc[output.group.eq("Healthy"), column].to_numpy(dtype=float)
            _, p = mannwhitneyu(gbm, healthy, alternative="two-sided")

            def statistic(values: np.ndarray, groups: np.ndarray) -> float:
                return float(values[groups == "GBM"].mean() - values[groups == "Healthy"].mean())

            permutation_p = stratified_permutation_p(
                output[column].to_numpy(dtype=float),
                labels.copy(),
                strata,
                statistic,
                iterations,
                seed=106804,
            )
            rows.append(
                {
                    "route": route,
                    "metric": metric,
                    "gbm_median": float(np.median(gbm)),
                    "healthy_median": float(np.median(healthy)),
                    "gbm_mean": float(gbm.mean()),
                    "healthy_mean": float(healthy.mean()),
                    "mean_difference_gbm_minus_healthy": float(gbm.mean() - healthy.mean()),
                    "mann_whitney_p": float(p),
                    "biofluid_stratified_permutation_p": float(permutation_p),
                }
            )
    frame = pd.DataFrame(rows)
    frame["mann_whitney_fdr"] = bh_adjust(frame["mann_whitney_p"])
    frame["biofluid_stratified_permutation_fdr"] = bh_adjust(
        frame["biofluid_stratified_permutation_p"]
    )
    return frame


def prediction_summaries(
    output: pd.DataFrame, routes: list[str], iterations: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    distributions, collapse, tests, top3 = [], [], [], []
    labels = output["group"].to_numpy()
    strata = output["biofluid"].to_numpy()
    for route in routes:
        column = f"{route}_top1_network"
        table = pd.crosstab(output["group"], output[column])
        if table.shape[1] > 1:
            chi2, p, _, _ = chi2_contingency(table)

            def statistic(_: np.ndarray, groups: np.ndarray) -> float:
                return float(chi2_contingency(pd.crosstab(groups, output[column]))[0])

            permutation_p = stratified_permutation_p(
                np.arange(len(output)),
                labels.copy(),
                strata,
                statistic,
                iterations,
                seed=804106,
            )
        else:
            chi2, p, permutation_p = 0.0, 1.0, 1.0
        tests.append(
            {
                "route": route,
                "chi_square": float(chi2),
                "chi_square_p": float(p),
                "biofluid_stratified_permutation_p": float(permutation_p),
                "active_networks": int(output[column].nunique()),
            }
        )
        for group, sub in output.groupby("group"):
            counts = sub[column].value_counts()
            proportions = counts / len(sub)
            hhi = float(np.square(proportions).sum())
            collapse.append(
                {
                    "route": route,
                    "group": group,
                    "n": int(len(sub)),
                    "active_networks": int(len(counts)),
                    "dominant_network": str(counts.index[0]),
                    "dominant_count": int(counts.iloc[0]),
                    "dominant_fraction": float(proportions.iloc[0]),
                    "hhi": hhi,
                    "effective_network_count": float(1.0 / hhi),
                    "collapse_flag_dominant_fraction_ge_0p70": bool(
                        proportions.iloc[0] >= 0.70
                    ),
                }
            )
            for network, count in counts.items():
                distributions.append(
                    {
                        "route": route,
                        "group": group,
                        "network": network,
                        "count": int(count),
                        "fraction": float(count / len(sub)),
                        "broad_anatomy": NETWORK_TO_PRIMARY_BROAD[network],
                    }
                )
            occurrences = Counter(
                item.strip()
                for value in sub[f"{route}_top3_network"]
                for item in str(value).split("|")
                if item.strip()
            )
            for network, count in occurrences.items():
                top3.append(
                    {
                        "route": route,
                        "group": group,
                        "network": network,
                        "samples_with_label_in_top3": int(count),
                        "fraction": float(count / len(sub)),
                        "broad_anatomy": NETWORK_TO_PRIMARY_BROAD[network],
                    }
                )
    tests = pd.DataFrame(tests)
    tests["chi_square_fdr"] = bh_adjust(tests["chi_square_p"])
    tests["biofluid_stratified_permutation_fdr"] = bh_adjust(
        tests["biofluid_stratified_permutation_p"]
    )
    return (
        pd.DataFrame(distributions),
        pd.DataFrame(top3),
        pd.DataFrame(collapse),
        tests,
    )


def compare_collapse(current: pd.DataFrame, previous_path: Path) -> pd.DataFrame:
    previous = pd.read_csv(previous_path)
    previous["dataset"] = "GSE228512_unselected_serum_EV"
    current = current.copy()
    current["dataset"] = "GSE106804_tumor_specific_EV"
    columns = [
        "dataset",
        "route",
        "group",
        "n",
        "active_networks",
        "dominant_network",
        "dominant_fraction",
        "hhi",
        "effective_network_count",
        "collapse_flag_dominant_fraction_ge_0p70",
    ]
    return pd.concat([previous[columns], current[columns]], ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument("--permutations", type=int, default=10000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cpm, metadata = load_data(args.data_dir)
    output = metadata.copy()

    baseline_scores, baseline_networks, baseline_genes, baseline_present = model_scores(
        cpm, BASELINE_MODEL, "logcpm"
    )
    output = add_route(output, "baseline_logcpm", baseline_scores, baseline_networks)
    output = production_predictions(cpm, output)

    adapted = np.load(ADAPTED_MODEL, allow_pickle=False)
    threshold = float(adapted["ood_max_correlation_threshold"][0])
    adapted_present = None
    for route, transform in (
        ("adapted_cpm", "cpm"),
        ("adapted_logcpm", "logcpm"),
        ("adapted_harmonized", "harmonized"),
    ):
        scores, networks, _, present = model_scores(cpm, ADAPTED_MODEL, transform)
        adapted_present = present
        output = add_route(output, route, scores, networks, threshold)

    score_routes = [
        "baseline_logcpm",
        "adapted_cpm",
        "adapted_logcpm",
        "adapted_harmonized",
    ]
    prediction_routes = ["baseline_production", *score_routes]
    continuous = continuous_tests(output, score_routes, args.permutations)
    distribution, top3, collapse, tests = prediction_summaries(
        output, prediction_routes, args.permutations
    )
    comparison = compare_collapse(
        collapse, GSE228_RESULTS / "prediction_collapse_metrics.csv"
    )

    output.to_csv(args.outdir / "gse106804_sample_predictions.csv", index=False, encoding="utf-8-sig")
    continuous.to_csv(args.outdir / "group_continuous_metric_tests.csv", index=False, encoding="utf-8-sig")
    distribution.to_csv(args.outdir / "network_broad_prediction_distribution.csv", index=False, encoding="utf-8-sig")
    top3.to_csv(args.outdir / "network_broad_top3_occurrence.csv", index=False, encoding="utf-8-sig")
    collapse.to_csv(args.outdir / "prediction_collapse_metrics.csv", index=False, encoding="utf-8-sig")
    tests.to_csv(args.outdir / "group_distribution_tests.csv", index=False, encoding="utf-8-sig")
    comparison.to_csv(
        args.outdir / "gse106804_vs_gse228512_collapse_comparison.csv",
        index=False,
        encoding="utf-8-sig",
    )

    summary = {
        "dataset": "GSE106804",
        "n_samples": int(len(output)),
        "group_counts": output.group.value_counts().to_dict(),
        "biofluid_counts": output.biofluid.value_counts().to_dict(),
        "biofluid_by_group": pd.crosstab(output.biofluid, output.group).to_dict(),
        "model_gene_coverage": {
            "baseline": {
                "present": int(baseline_present.sum()),
                "total": int(len(baseline_genes)),
                "fraction": float(baseline_present.mean()),
            },
            "adapted": {
                "present": int(adapted_present.sum()),
                "total": int(len(adapted_present)),
                "fraction": float(adapted_present.mean()),
            },
        },
        "ood_threshold": threshold,
        "ood_acceptance": {
            route: {
                group: {
                    "accepted": int(
                        output.loc[output.group.eq(group), f"{route}_ood_accepted"].sum()
                    ),
                    "n": int(output.group.eq(group).sum()),
                    "coverage": float(
                        output.loc[output.group.eq(group), f"{route}_ood_accepted"].mean()
                    ),
                }
                for group in ("GBM", "Healthy")
            }
            for route in ("adapted_cpm", "adapted_logcpm", "adapted_harmonized")
        },
        "baseline_production_pairwise_rescue": {
            "n_switched": int(output.baseline_production_pairwise_switched.sum()),
            "fraction_switched": float(output.baseline_production_pairwise_switched.mean()),
        },
        "interpretation": (
            "Tumor-specific EV enrichment is assessed by comparison with GSE228512. "
            "No localization accuracy is computed because patient-level tumor locations are unavailable."
        ),
    }
    (args.outdir / "analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
