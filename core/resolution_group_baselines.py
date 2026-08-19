"""Deterministic current-formal resolution-group random baselines."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np


RNG_SEED = 20260629
N_WEIGHTED_RANDOM_DRAWS = 10000
FORMAL_ROUTE_FAMILY = "hybrid_projected_network_logcpm_exact"
EXPECTED = {
    "LOSO": {"n": 814, "top1": 368, "top3": 590},
    "LOMO": {"n": 812, "top1": 344, "top3": 569},
}

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_PATHS = {
    "LOSO": ROOT
    / "reproducibility"
    / "p2_publication_completeness"
    / "formal_loso_resolution_group_detail.csv",
    "LOMO": ROOT
    / "reproducibility"
    / "p2_publication_completeness"
    / "formal_lomo_resolution_group_detail.csv",
}


def load_formal_rows(path: Path, endpoint: str) -> tuple[list[dict[str, str]], list[str]]:
    """Read only the current hybrid formal universe for one endpoint."""

    if endpoint not in EXPECTED:
        raise ValueError(f"Unsupported endpoint: {endpoint}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        fields = list(reader.fieldnames)
        rows = [
            row
            for row in reader
            if str(row.get("route_family", "")).strip() == FORMAL_ROUTE_FAMILY
        ]
    required = {"sample_id", "true_resolution_group", "group_hit1", "group_hit3", "n_candidate_groups"}
    if missing := required - set(fields):
        raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")
    expected = EXPECTED[endpoint]
    if len(rows) != expected["n"]:
        raise ValueError(f"{endpoint} formal group detail must contain {expected['n']} rows")
    if len({row["sample_id"] for row in rows}) != expected["n"]:
        raise ValueError(f"{endpoint} formal group detail has duplicate sample IDs")
    top1 = sum(int(float(row["group_hit1"])) for row in rows)
    top3 = sum(int(float(row["group_hit3"])) for row in rows)
    if (top1, top3) != (expected["top1"], expected["top3"]):
        raise ValueError(f"{endpoint} observed group counts do not match the locked endpoint")
    return rows, fields


def uniform_topk_baseline(
    rows: Iterable[dict[str, str]], top_k: int = 3
) -> tuple[float, str]:
    """Mean sample-specific min(TopK / n_candidate_groups, 1)."""

    selected = list(rows)
    values = [
        min(top_k / max(float(row["n_candidate_groups"]), 1.0), 1.0)
        for row in selected
    ]
    return float(sum(values) / len(values)), "sample_specific_uniform_top3_over_n_candidate_groups"


def weighted_topk_baseline(
    rows: Iterable[dict[str, str]], top_k: int = 3
) -> tuple[float, str]:
    """Global truth-prevalence sampling without replacement, as in the original generator."""

    selected = list(rows)
    labels = sorted({str(row["true_resolution_group"]) for row in selected})
    label_to_index = {label: index for index, label in enumerate(labels)}
    counts = np.array(
        [sum(row["true_resolution_group"] == label for row in selected) for label in labels],
        dtype=float,
    )
    weights = counts / counts.sum()
    truth_indices = np.array(
        [label_to_index[row["true_resolution_group"]] for row in selected], dtype=int
    )
    candidate_k = np.array(
        [
            min(top_k, max(int(float(row["n_candidate_groups"])), 1))
            for row in selected
        ],
        dtype=int,
    )
    rng = np.random.default_rng(RNG_SEED + top_k + len(selected))
    inclusion_by_k: dict[int, np.ndarray] = {}
    for k in sorted(set(int(value) for value in candidate_k)):
        draw_counts = np.zeros(len(labels), dtype=float)
        for _ in range(N_WEIGHTED_RANDOM_DRAWS):
            draw = rng.choice(len(labels), size=min(k, len(labels)), replace=False, p=weights)
            draw_counts[draw] += 1.0
        inclusion_by_k[k] = draw_counts / N_WEIGHTED_RANDOM_DRAWS
    probabilities = np.array(
        [inclusion_by_k[int(k)][int(index)] for k, index in zip(candidate_k, truth_indices)]
    )
    return (
        float(probabilities.mean()),
        "global_truth_prior_weighted_without_replacement; candidate_count_limits_k_only",
    )


def compute_baseline_record(rows: list[dict[str, str]], endpoint: str) -> dict[str, object]:
    uniform, uniform_formula = uniform_topk_baseline(rows)
    weighted, weighted_formula = weighted_topk_baseline(rows)
    n = len(rows)
    return {
        "endpoint": endpoint,
        "top_k": 3,
        "n_profiles": n,
        "observed_hits": sum(int(float(row["group_hit3"])) for row in rows),
        "observed_rate": sum(int(float(row["group_hit3"])) for row in rows) / n,
        "uniform_random_rate": uniform,
        "weighted_random_rate": weighted,
        "uniform_formula": uniform_formula,
        "weighted_formula": weighted_formula,
        "rng_seed": RNG_SEED,
        "n_weighted_random_draws": N_WEIGHTED_RANDOM_DRAWS,
    }
