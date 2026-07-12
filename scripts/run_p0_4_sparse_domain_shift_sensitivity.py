#!/usr/bin/env python
"""P0-4 locked-route sparse-expression sensitivity analysis.

The training reference, feature selection, hierarchy construction, and scoring
parameters are held fixed within each strict LOSO fold.  Only the held-out
query count vector is perturbed.  The perturbation combines abundance-weighted
gene retention (dropout) with binomial read-depth thinning, and is therefore a
simulation of sparse-expression stress rather than real cfRNA validation.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import align_matrices, compute_logcpm, map_index_to_symbols, read_bo2023_gene_matrix, read_gene_map
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_MODEL_GENES, DEFAULT_SAMPLE_INFO, DEFAULT_VSD, read_locked_model_genes
from scripts.run_bo2023_hybrid_formal_loso import (  # noqa: E402
    EXACT_ROUTE,
    GROUP_ROUTE,
    build_centroids,
    build_region_reference,
    exact_row,
    network_row,
    read_metadata,
    zscore,
)
from scripts.run_bo2023_loso_validation import correlation_scores
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups, score_route


DEFAULT_OUTDIR = ROOT / "reports" / "p0_4_sparse_domain_shift_20260711"
SCENARIOS = (
    ("baseline", 1.00, 1.00),
    ("mild", 0.50, 0.80),
    ("moderate", 0.20, 0.60),
    ("severe", 0.05, 0.40),
    ("extreme", 0.01, 0.20),
)


@dataclass(frozen=True)
class FoldCache:
    sample_id: str
    monkey_id: str
    truth_network: str
    truth_region: str
    network_reference: np.ndarray
    projection_slope: np.ndarray
    projection_intercept: np.ndarray
    projection_low: np.ndarray
    projection_high: np.ndarray
    candidates: list[str]
    candidate_reference: np.ndarray | None
    gene_order: np.ndarray | None
    local_rows: np.ndarray | None
    annotations: dict[str, dict[str, Any]] | None


def projection_parameters(logcpm: np.ndarray, vsd: np.ndarray, rows: np.ndarray, sample_idx: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = logcpm[rows, :]
    y = vsd[rows, :]
    mask = np.ones(x.shape[1], dtype=bool)
    mask[sample_idx] = False
    xt, yt = x[:, mask].astype(np.float64), y[:, mask].astype(np.float64)
    x_mean, y_mean = xt.mean(axis=1), yt.mean(axis=1)
    xc, yc = xt - x_mean[:, None], yt - y_mean[:, None]
    denom = np.square(xc).sum(axis=1)
    slope = np.divide((xc * yc).sum(axis=1), denom, out=np.zeros(len(rows), dtype=float), where=denom > 0)
    intercept = y_mean - slope * x_mean
    return slope, intercept, np.quantile(yt, 0.005, axis=1), np.quantile(yt, 0.995, axis=1)


def perturb_counts(counts: np.ndarray, depth_fraction: float, gene_retention: float, rng: np.random.Generator) -> np.ndarray:
    """Keep a fixed fraction of detected genes with abundance-weighted sampling, then thin reads."""
    original = np.asarray(np.rint(np.clip(counts, 0, None)), dtype=np.int64)
    positive = np.flatnonzero(original > 0)
    if not len(positive):
        return np.zeros_like(original)
    keep_n = max(1, int(round(len(positive) * gene_retention)))
    if keep_n == len(positive):
        kept = positive
    else:
        # Higher-abundance transcripts are more likely to be detected in sparse libraries.
        weights = np.sqrt(original[positive].astype(float))
        weights /= weights.sum()
        kept = rng.choice(positive, size=keep_n, replace=False, p=weights)
    result = np.zeros_like(original)
    if depth_fraction >= 1.0:
        result[kept] = original[kept]
    else:
        result[kept] = rng.binomial(original[kept], depth_fraction)
    return result


def normalized_entropy(scores: np.ndarray) -> float:
    shifted = scores - np.max(scores)
    p = np.exp(shifted)
    p /= p.sum()
    return float(-(p * np.log(np.maximum(p, 1e-12))).sum() / np.log(len(p)))


def donor_bootstrap(values: pd.DataFrame, metric: str, rng: np.random.Generator, n_bootstrap: int) -> tuple[float, float]:
    donor_values = [group[metric].to_numpy(dtype=float) for _, group in values.groupby("monkey_id", sort=True)]
    draws = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        idx = rng.integers(0, len(donor_values), size=len(donor_values))
        draws[i] = np.concatenate([donor_values[j] for j in idx]).mean()
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def scenario_summary(detail: pd.DataFrame, n_bootstrap: int, seed: int) -> pd.DataFrame:
    rows = []
    for scenario, frame in detail.groupby("scenario", sort=False):
        rng = np.random.default_rng(seed + list(detail["scenario"].drop_duplicates()).index(scenario))
        row: dict[str, Any] = {
            "scenario": scenario,
            "depth_fraction": float(frame["depth_fraction"].iloc[0]),
            "target_gene_retention": float(frame["target_gene_retention"].iloc[0]),
            "n_rows": int(len(frame)),
            "n_samples": int(frame["sample_id"].nunique()),
            "n_donors": int(frame["monkey_id"].nunique()),
        }
        for metric in (
            "network_hit1", "network_hit3", "group_hit3", "exact_hit3", "gene_coverage_fraction",
            "network_marker_coverage", "local_marker_coverage", "network_entropy", "network_margin",
            "low_confidence_warning",
        ):
            valid = frame[["monkey_id", metric]].dropna(subset=[metric]).copy()
            row[f"{metric}_n"] = int(len(valid))
            row[metric] = float(valid[metric].mean())
            low, high = donor_bootstrap(valid, metric, rng, n_bootstrap)
            row[f"{metric}_donor_bootstrap_ci_low"] = low
            row[f"{metric}_donor_bootstrap_ci_high"] = high
        rows.append(row)
    return pd.DataFrame(rows)


def plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    x = np.arange(len(summary))
    labels = summary["scenario"].tolist()
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for metric, label, ax in (
        ("network_hit3", "Network Top3 accuracy", axes[0, 0]),
        ("group_hit3", "Resolution-group Top3 accuracy", axes[0, 1]),
        ("exact_hit3", "Exact-region Top3 accuracy", axes[1, 0]),
        ("network_marker_coverage", "Detected Network-marker coverage", axes[1, 1]),
    ):
        y = summary[metric].to_numpy()
        lo = y - summary[f"{metric}_donor_bootstrap_ci_low"].to_numpy()
        hi = summary[f"{metric}_donor_bootstrap_ci_high"].to_numpy() - y
        ax.errorbar(x, y, yerr=np.vstack([lo, hi]), marker="o", capsize=4, color="#2d6ca2")
        ax.set_title(label)
        ax.set_ylim(0, 1.05)
        ax.set_xticks(x, labels, rotation=25, ha="right")
        ax.set_ylabel("Proportion")
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Locked-route sensitivity to simulated sparse expression")
    fig.savefig(outdir / "p0_4_sparse_sensitivity_curves.png", dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_CLEANED_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_MODEL_GENES)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--max-samples", type=int, default=0, help="Smoke-test limiter; 0 analyses all samples.")
    parser.add_argument("--seed", type=int, default=20260711)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--summarize-only", action="store_true", help="Regenerate summary/plot from an existing detail CSV.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.summarize_only:
        detail_path = args.outdir / "p0_4_sparse_sensitivity_sample_detail.csv"
        detail = pd.read_csv(detail_path)
        summary = scenario_summary(detail, args.n_bootstrap, args.seed)
        summary.to_csv(args.outdir / "p0_4_sparse_sensitivity_summary.csv", index=False)
        plot_summary(summary, args.outdir)
        print(summary.to_string(index=False))
        return 0

    gene_map = read_gene_map(args.gene_map)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.counts, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(args.vsd, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata(args.sample_info, args.sample_sheet, "Region", "SaleemNetworks")
    donor_info = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet, usecols=["No.", "MonkeyID"])
    donor_info["sample_id"] = donor_info["No."].astype(str).str.strip()
    donor_by_sample = donor_info.drop_duplicates("sample_id").set_index("sample_id")["MonkeyID"].astype(str).to_dict()
    samples = [s for s in samples if s in metadata.index and s in donor_by_sample]
    counts, vsd = counts.loc[genes, samples], vsd.loc[genes, samples]
    logcpm = compute_logcpm(counts)
    count_values, vsd_values, logcpm_values = counts.to_numpy(dtype=np.float32), vsd.to_numpy(dtype=np.float32), logcpm.to_numpy(dtype=np.float32)
    region_labels = metadata.loc[samples, "region_id"].astype(str).to_numpy()
    network_labels = metadata.loc[samples, "network_id"].astype(str).to_numpy()
    networks = sorted(set(network_labels))
    all_regions = sorted(set(region_labels))
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    locked_genes = [g for g in read_locked_model_genes(args.locked_model_genes) if g in gene_to_idx]
    locked_rows = np.asarray([gene_to_idx[g] for g in locked_genes], dtype=int)
    if len(locked_rows) < 20:
        raise RuntimeError("Too few locked Network genes after symbol mapping")

    records: list[dict[str, Any]] = []
    max_samples = len(samples) if args.max_samples == 0 else min(args.max_samples, len(samples))
    for sample_idx, sample_id in enumerate(samples[:max_samples]):
        train_idx = np.setdiff1d(np.arange(len(samples)), [sample_idx])
        truth_network, truth_region = str(network_labels[sample_idx]), str(region_labels[sample_idx])
        network_reference = build_centroids(vsd_values[locked_rows, :], network_labels, networks, train_idx)
        slope, intercept, clip_low, clip_high = projection_parameters(logcpm_values, vsd_values, locked_rows, sample_idx)
        train_regions = sorted(set(region_labels[train_idx]))
        region_training = {region: train_idx[region_labels[train_idx] == region] for region in train_regions}
        region_evaluable = truth_region in region_training

        # Cache all fold-local references before perturbing the held-out query.
        cached_by_beam: dict[tuple[str, ...], tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, dict[str, Any]]]] = {}
        original_detected = int((count_values[:, sample_idx] > 0).sum())
        for scenario_i, (scenario, depth_fraction, gene_retention) in enumerate(SCENARIOS):
            repeats = 1 if scenario == "baseline" else args.replicates
            for repeat in range(repeats):
                rng = np.random.default_rng(args.seed + sample_idx * 1000 + scenario_i * 100 + repeat)
                perturbed_counts = perturb_counts(count_values[:, sample_idx], depth_fraction, gene_retention, rng)
                total = int(perturbed_counts.sum())
                query_logcpm = np.log1p(perturbed_counts / total * 1_000_000.0) if total else np.zeros(len(genes), dtype=float)
                projected = np.clip(slope * query_logcpm[locked_rows] + intercept, clip_low, clip_high)
                network_scores = correlation_scores(network_reference, projected)
                order = np.argsort(network_scores)[::-1]
                network_top = [networks[int(i)] for i in order[:3]]
                network_margin = float(network_scores[order[0]] - network_scores[order[1]])
                row: dict[str, Any] = {
                    "sample_id": sample_id,
                    "monkey_id": str(donor_by_sample[sample_id]),
                    "scenario": scenario,
                    "depth_fraction": depth_fraction,
                    "target_gene_retention": gene_retention,
                    "replicate": repeat,
                    "truth_network": truth_network,
                    "truth_region": truth_region,
                    "gene_coverage_fraction": float((perturbed_counts > 0).sum() / max(original_detected, 1)),
                    "network_marker_coverage": float((perturbed_counts[locked_rows] > 0).mean()),
                    "network_entropy": normalized_entropy(network_scores),
                    "network_margin": network_margin,
                    "network_hit1": int(network_top[0] == truth_network),
                    "network_hit3": int(truth_network in network_top),
                    "group_hit3": np.nan,
                    "exact_hit3": np.nan,
                    "local_marker_coverage": np.nan,
                    "region_evaluable": bool(region_evaluable),
                }
                # The production route has no validated absolute correlation-margin abstention threshold.
                # Retain margin as a continuous diagnostic; the warning itself is coverage-only and prespecified.
                row["low_confidence_warning"] = int(row["network_marker_coverage"] < 0.50)
                if region_evaluable:
                    candidates = tuple(sorted(region for region in train_regions if metadata.loc[samples[int(region_training[region][0])], "network_id"] in set(network_top)))
                    if len(candidates) >= 2:
                        if candidates not in cached_by_beam:
                            candidate_training = {region: region_training[region] for region in candidates}
                            gene_order, _ = select_group_discriminative_genes(logcpm_values, list(candidates), candidate_training, 200)
                            candidate_reference = build_region_reference(logcpm_values, list(candidates), candidate_training)
                            local_rows = gene_order[: min(200, len(gene_order))]
                            train_meta = metadata.loc[[samples[int(i)] for i in train_idx]]
                            assignment = {region: (lambda nets: nets[0] if len(nets) == 1 else None)(sorted(train_meta.loc[train_meta["region_id"] == region, "network_id"].astype(str).unique())) for region in candidates}
                            annotations, _ = build_resolution_groups(logcpm_values, list(candidates), candidate_training, assignment, local_rows, 8, 3, 2, 0.15, 0.95, 0.90, 8)
                            cached_by_beam[candidates] = candidate_reference, gene_order, local_rows, annotations
                        candidate_reference, gene_order, local_rows, annotations = cached_by_beam[candidates]
                        scores50 = correlation_scores(candidate_reference, query_logcpm, gene_order[: min(50, len(gene_order))])
                        scores100 = correlation_scores(candidate_reference, query_logcpm, gene_order[: min(100, len(gene_order))])
                        fused = 0.25 * zscore(scores50) + 0.75 * zscore(scores100)
                        ranked_exact = [candidates[i] for i in np.argsort(fused)[::-1]]
                        exact = exact_row(sample_id, truth_region, truth_network, network_top, ranked_exact, len(all_regions))
                        group_scores = correlation_scores(candidate_reference, query_logcpm, local_rows)
                        ranked_group = [candidates[i] for i in np.argsort(group_scores)[::-1]]
                        group = score_route(GROUP_ROUTE, sample_id, truth_region, truth_network, network_top, ranked_group, annotations, len(all_regions))
                        row["group_hit3"] = int(group["group_hit3"])
                        row["exact_hit3"] = int(exact["hit3"])
                        row["local_marker_coverage"] = float((perturbed_counts[local_rows] > 0).mean()) if len(local_rows) else np.nan
                records.append(row)
        if (sample_idx + 1) % 100 == 0:
            print(f"processed {sample_idx + 1}/{len(samples)} samples", flush=True)

    detail = pd.DataFrame(records)
    summary = scenario_summary(detail, args.n_bootstrap, args.seed)
    detail.to_csv(args.outdir / "p0_4_sparse_sensitivity_sample_detail.csv", index=False)
    summary.to_csv(args.outdir / "p0_4_sparse_sensitivity_summary.csv", index=False)
    plot_summary(summary, args.outdir)
    methods = {
        "design": "strict LOSO held-out query perturbation; reference and all fold-local feature/group construction remain unperturbed and locked",
        "simulation": "abundance-weighted retention of detected genes followed by binomial read-depth thinning",
        "scenarios": [{"name": n, "depth_fraction": d, "target_gene_retention": r} for n, d, r in SCENARIOS],
        "replicates_per_nonbaseline_scenario": args.replicates,
        "seed": args.seed,
        "donor_inference": f"{args.n_bootstrap} donor bootstrap replicates; nine monkeys are resampled as clusters",
        "warning_rule": "predefined low-confidence warning if detected Network-marker coverage <0.50; score margin and entropy are reported as continuous diagnostics because this locked route has no externally calibrated absolute margin threshold",
        "interpretation_limit": "This is simulated sparse-expression sensitivity analysis, not clinical cfRNA localization validation.",
    }
    (args.outdir / "P0_4_METHODS.json").write_text(json.dumps(methods, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
