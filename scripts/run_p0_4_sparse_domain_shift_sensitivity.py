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
import hashlib
import json
import shutil
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


DEFAULT_OUTDIR = ROOT / "reports" / "p0_8_sparse_domain_shift_20260716"
BOOTSTRAP_SEED = 20260716
EXPECTED_INPUT_SHA256 = {
    "counts": "1FB3A512DA11AB0C327C07C114DA3B9C38CAB0A504682F2C7C036EEDB3C7561A",
    "vsd": "286AEAB66B21B7FA012FAC8CEAA24497894327E0736F9F6B200334C57089A1B3",
    "sample_info": "9A2FE2BEC1475F6AD613883D0FF5925B1E6BA36E800CAA922C35D4F8AE7D3645",
    "gene_map": "24E0545A478ED2643322F994627A0D1C8BFAC3061F6A755EAE49080B2B92A78A",
    "locked_model_genes": "FCAF3A1927AA0E7B55513C9E3486333D67742E3CEAF0A7FB650A47B1DDA15929",
}
EXPECTED_BASELINE = {
    "network_hit1": (483, 819),
    "network_hit3": (753, 819),
    "group_hit1": (368, 814),
    "group_hit3": (590, 814),
    "exact_hit1": (182, 814),
    "exact_hit3": (368, 814),
}
METRICS = (
    "network_hit1", "network_hit3", "group_hit1", "group_hit3", "exact_hit1", "exact_hit3",
    "gene_coverage_fraction", "network_marker_coverage", "local_marker_coverage", "network_entropy",
    "network_margin", "low_confidence_warning",
)
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def validate_input_contract(paths: dict[str, Path]) -> dict[str, dict[str, Any]]:
    audit: dict[str, dict[str, Any]] = {}
    for name, expected in EXPECTED_INPUT_SHA256.items():
        path = paths[name].resolve()
        observed = sha256_file(path)
        audit[name] = {"path": str(path), "sha256": observed, "expected_sha256": expected, "matched": observed == expected}
    failures = [name for name, row in audit.items() if not row["matched"]]
    if failures:
        raise RuntimeError(f"Frozen-input SHA-256 contract failed: {', '.join(failures)}")
    return audit


def read_detail_csv(path: Path) -> pd.DataFrame:
    """Read persisted detail without splitting numeric-looking donor IDs by chunk dtype."""
    return pd.read_csv(path, dtype={"monkey_id": "string"}, low_memory=False)


def repeat_summary(detail: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (scenario, repeat), frame in detail.groupby(["scenario", "replicate"], sort=False):
        for metric in METRICS:
            valid = frame[["monkey_id", metric]].dropna(subset=[metric])
            if valid.empty:
                continue
            donor_means = valid.groupby("monkey_id", sort=True)[metric].mean()
            common = {
                "scenario": scenario,
                "replicate": int(repeat),
                "depth_fraction": float(frame["depth_fraction"].iloc[0]),
                "target_gene_retention": float(frame["target_gene_retention"].iloc[0]),
                "metric": metric,
                "n_rows": int(len(valid)),
                "n_samples": int(frame.loc[valid.index, "sample_id"].nunique()),
                "n_donors": int(len(donor_means)),
            }
            rows.append({**common, "estimator": "sample_weighted", "value": float(valid[metric].mean())})
            rows.append({**common, "estimator": "donor_macro", "value": float(donor_means.mean())})
    return pd.DataFrame(rows)


def _cluster_bootstrap_draws(frame: pd.DataFrame, metric: str, n_bootstrap: int, seed: int) -> dict[str, np.ndarray]:
    """Resample donors inside each repeat, then average estimates across repeats."""
    repeat_draws: dict[str, list[np.ndarray]] = {"sample_weighted": [], "donor_macro": []}
    for repeat, repeat_frame in frame.groupby("replicate", sort=True):
        valid = repeat_frame[["monkey_id", metric]].dropna(subset=[metric])
        donor = valid.groupby("monkey_id", sort=True)[metric].agg(["sum", "count", "mean"])
        if donor.empty:
            continue
        rng = np.random.default_rng(np.random.SeedSequence([seed, int(repeat)]))
        indices = rng.integers(0, len(donor), size=(n_bootstrap, len(donor)))
        sums, counts, means = donor["sum"].to_numpy(), donor["count"].to_numpy(), donor["mean"].to_numpy()
        repeat_draws["sample_weighted"].append(sums[indices].sum(axis=1) / counts[indices].sum(axis=1))
        repeat_draws["donor_macro"].append(means[indices].mean(axis=1))
    if not repeat_draws["sample_weighted"]:
        return {}
    return {name: np.vstack(draws).mean(axis=0) for name, draws in repeat_draws.items()}


def across_repeat_summary(detail: pd.DataFrame, repeats: pd.DataFrame, n_bootstrap: int, seed: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    scenario_order = list(detail["scenario"].drop_duplicates())
    for scenario_i, scenario in enumerate(scenario_order):
        scenario_detail = detail[detail["scenario"] == scenario]
        for metric_i, metric in enumerate(METRICS):
            metric_repeats = repeats[(repeats["scenario"] == scenario) & (repeats["metric"] == metric)]
            if metric_repeats.empty:
                continue
            draws = _cluster_bootstrap_draws(
                scenario_detail, metric, n_bootstrap, seed + scenario_i * 1000 + metric_i * 10,
            )
            for estimator, estimator_frame in metric_repeats.groupby("estimator", sort=False):
                values = estimator_frame.sort_values("replicate")["value"].to_numpy(dtype=float)
                bootstrap = draws[estimator]
                rows.append({
                    "scenario": scenario,
                    "depth_fraction": float(scenario_detail["depth_fraction"].iloc[0]),
                    "target_gene_retention": float(scenario_detail["target_gene_retention"].iloc[0]),
                    "metric": metric,
                    "estimator": estimator,
                    "n_repeats": int(len(values)),
                    "n_donors": int(scenario_detail["monkey_id"].nunique()),
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=0)),
                    "min": float(values.min()),
                    "max": float(values.max()),
                    "mc_q025": float(np.quantile(values, 0.025)),
                    "mc_q975": float(np.quantile(values, 0.975)),
                    "donor_bootstrap_ci_low": float(np.quantile(bootstrap, 0.025)),
                    "donor_bootstrap_ci_high": float(np.quantile(bootstrap, 0.975)),
                })
    return pd.DataFrame(rows)


def validate_baseline(detail: pd.DataFrame, ontology: dict[str, int], full_run: bool) -> dict[str, Any]:
    baseline = detail[detail["scenario"] == "baseline"]
    observed = {
        metric: (int(baseline[metric].sum()), int(baseline[metric].notna().sum()))
        for metric in EXPECTED_BASELINE
    }
    contract = {
        "expected": {key: list(value) for key, value in EXPECTED_BASELINE.items()},
        "observed": {key: list(value) for key, value in observed.items()},
        "n_donors": int(baseline["monkey_id"].nunique()),
        "ontology": ontology,
        "enable_pairwise_rescue": False,
        "full_run": bool(full_run),
    }
    if not full_run:
        raise RuntimeError("Canonical baseline gate requires all 819 samples; --max-samples is not permitted for evidence output")
    failures = [key for key in EXPECTED_BASELINE if observed[key] != EXPECTED_BASELINE[key]]
    if contract["n_donors"] != 9:
        failures.append("n_donors")
    if ontology != {"regions": 110, "networks": 10, "total_labels": 120}:
        failures.append("ontology")
    if failures:
        raise RuntimeError(f"Canonical baseline gate failed: {', '.join(failures)}")
    contract["passed"] = True
    return contract


def plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    selected = summary[(summary["estimator"] == "donor_macro") & summary["metric"].isin(("network_hit3", "group_hit3", "exact_hit3", "network_marker_coverage"))]
    labels = list(selected["scenario"].drop_duplicates())
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for metric, label, ax in (
        ("network_hit3", "Network Top3 accuracy", axes[0, 0]),
        ("group_hit3", "Resolution-group Top3 accuracy", axes[0, 1]),
        ("exact_hit3", "Exact-region Top3 accuracy", axes[1, 0]),
        ("network_marker_coverage", "Detected Network-marker coverage", axes[1, 1]),
    ):
        frame = selected[selected["metric"] == metric].set_index("scenario").loc[labels]
        y = frame["mean"].to_numpy()
        lo = y - frame["donor_bootstrap_ci_low"].to_numpy()
        hi = frame["donor_bootstrap_ci_high"].to_numpy() - y
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
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--max-samples", type=int, default=0, help="Development-only limiter; cannot pass the evidence gate.")
    parser.add_argument("--seed", type=int, default=20260711)
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED)
    parser.add_argument("--n-bootstrap", type=int, default=50000)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--summarize-only", action="store_true", help="Regenerate summary/plot from an existing detail CSV.")
    args = parser.parse_args()
    if args.summarize_only:
        detail_path = args.outdir / "p0_4_sparse_sensitivity_sample_detail.csv"
        detail = read_detail_csv(detail_path)
        repeats = repeat_summary(detail)
        summary = across_repeat_summary(detail, repeats, args.n_bootstrap, args.bootstrap_seed)
        repeats.to_csv(args.outdir / "p0_8_sparse_sensitivity_per_repeat.csv", index=False)
        summary.to_csv(args.outdir / "p0_8_sparse_sensitivity_across_repeats.csv", index=False)
        plot_summary(summary, args.outdir)
        print(summary.to_string(index=False))
        return 0

    if args.outdir.exists():
        raise FileExistsError(f"Refusing to overwrite existing evidence directory: {args.outdir}")
    if args.replicates < 1:
        raise ValueError("--replicates must be at least 1")
    if args.n_bootstrap < 1:
        raise ValueError("--n-bootstrap must be at least 1")
    input_audit = validate_input_contract({
        "counts": args.counts,
        "vsd": args.vsd,
        "sample_info": args.sample_info,
        "gene_map": args.gene_map,
        "locked_model_genes": args.locked_model_genes,
    })
    args.outdir.mkdir(parents=True)

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
    ontology = {
        "regions": len(all_regions),
        "networks": len(networks),
        "total_labels": len(all_regions) + len(networks),
    }
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
                rng_seed = args.seed + sample_idx * 1000 + scenario_i * 100 + repeat
                rng = np.random.default_rng(rng_seed)
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
                    "rng_seed": rng_seed,
                    "truth_network": truth_network,
                    "truth_region": truth_region,
                    "gene_coverage_fraction": float((perturbed_counts > 0).sum() / max(original_detected, 1)),
                    "network_marker_coverage": float((perturbed_counts[locked_rows] > 0).mean()),
                    "network_entropy": normalized_entropy(network_scores),
                    "network_margin": network_margin,
                    "network_hit1": int(network_top[0] == truth_network),
                    "network_hit3": int(truth_network in network_top),
                    "group_hit3": np.nan,
                    "group_hit1": np.nan,
                    "exact_hit3": np.nan,
                    "exact_hit1": np.nan,
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
                        row["group_hit1"] = int(group["group_hit1"])
                        row["exact_hit3"] = int(exact["hit3"])
                        row["exact_hit1"] = int(exact["hit1"])
                        row["local_marker_coverage"] = float((perturbed_counts[local_rows] > 0).mean()) if len(local_rows) else np.nan
                records.append(row)
        if (sample_idx + 1) % 100 == 0:
            print(f"processed {sample_idx + 1}/{len(samples)} samples", flush=True)

    detail = pd.DataFrame(records)
    baseline_contract = validate_baseline(detail, ontology, max_samples == len(samples))
    repeats = repeat_summary(detail)
    summary = across_repeat_summary(detail, repeats, args.n_bootstrap, args.bootstrap_seed)
    detail.to_csv(args.outdir / "p0_4_sparse_sensitivity_sample_detail.csv", index=False)
    repeats.to_csv(args.outdir / "p0_8_sparse_sensitivity_per_repeat.csv", index=False)
    summary.to_csv(args.outdir / "p0_8_sparse_sensitivity_across_repeats.csv", index=False)
    plot_summary(summary, args.outdir)
    methods = {
        "design": "strict LOSO held-out query perturbation; reference and all fold-local feature/group construction remain unperturbed and locked",
        "simulation": "abundance-weighted retention of detected genes followed by binomial read-depth thinning",
        "scenarios": [{"name": n, "depth_fraction": d, "target_gene_retention": r} for n, d, r in SCENARIOS],
        "replicates_per_nonbaseline_scenario": args.replicates,
        "perturbation_seed": args.seed,
        "perturbation_seed_formula": "seed + sample_idx * 1000 + scenario_idx * 100 + repeat; saved as rng_seed on every detail row",
        "bootstrap_seed": args.bootstrap_seed,
        "donor_inference": f"{args.n_bootstrap} vectorized donor bootstrap draws; donors are resampled independently within repeat and estimates are then averaged across repeats",
        "estimators": ["sample_weighted", "donor_macro"],
        "monte_carlo_summary": ["mean", "population SD", "min", "max", "2.5th percentile", "97.5th percentile"],
        "frozen_route": {"canonical_regions": 110, "canonical_networks": 10, "total_labels": 120, "enable_pairwise_rescue": False},
        "warning_rule": "predefined low-confidence warning if detected Network-marker coverage <0.50; score margin and entropy are reported as continuous diagnostics because this locked route has no externally calibrated absolute margin threshold",
        "interpretation_limit": "This is simulated sparse-expression sensitivity analysis, not clinical cfRNA localization validation.",
    }
    (args.outdir / "P0_4_METHODS.json").write_text(json.dumps(methods, indent=2), encoding="utf-8")
    (args.outdir / "input_audit.json").write_text(json.dumps(input_audit, indent=2), encoding="utf-8")
    (args.outdir / "baseline_gate.json").write_text(json.dumps(baseline_contract, indent=2), encoding="utf-8")
    shutil.copy2(__file__, args.outdir / "executed_run_p0_4_sparse_domain_shift_sensitivity.py")
    output_files = sorted(path for path in args.outdir.iterdir() if path.is_file() and path.name != "SHA256SUMS.txt")
    (args.outdir / "SHA256SUMS.txt").write_text(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in output_files), encoding="utf-8",
    )
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
