#!/usr/bin/env python3
"""Reproducible calculations for BrainTrace round-4 P1-STAT1--4."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


OUT = Path("manuscript/calculations")
OUT.mkdir(parents=True, exist_ok=True)


def bh_adjust(labels, pvalues):
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    adjusted = [0.0] * len(pvalues)
    running = 1.0
    m = len(pvalues)
    for rank_index in range(m - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        running = min(running, pvalues[idx] * m / rank)
        adjusted[idx] = min(1.0, running)
    return [{"endpoint": x, "raw_p": p, "bh_q_m4": q}
            for x, p, q in zip(labels, pvalues, adjusted)]


labels = ["Network Top1", "Network Top3", "Resolution-group Top3", "Exact-region Top3"]
pvalues = [0.031250, 0.375000, 0.593750, 0.324219]
bh = bh_adjust(labels, pvalues)

# The MDE simulation used donor-level paired rates, so ICC was not a simulation
# input. This grid quantifies how a hypothetical sample-level analysis would
# lose information under within-donor correlation.
n_samples = 819
n_donors = 9
mean_cluster_size = n_samples / n_donors
pooled_network_top3 = (0.9194 + 0.9158) / 2
bernoulli_variance = pooled_network_top3 * (1 - pooled_network_top3)
icc_rows = []
for icc in [0.00, 0.01, 0.05, 0.10, 0.30]:
    design_effect = 1 + (mean_cluster_size - 1) * icc
    icc_rows.append({
        "assumed_icc": icc,
        "mean_cluster_size": mean_cluster_size,
        "design_effect": design_effect,
        "effective_sample_size": n_samples / design_effect,
        "total_bernoulli_variance": bernoulli_variance,
        "between_donor_component": icc * bernoulli_variance,
        "within_donor_sample_component": (1 - icc) * bernoulli_variance,
    })

# Normal-theory prediction interval for one additional perturbation repeat:
# mean +/- t_(.975,29) * s * sqrt(1 + 1/30). This is not a per-query clinical PI.
t_975_df29 = 2.045229642
multiplier = t_975_df29 * math.sqrt(1 + 1 / 30)
s15 = {
    "Mild": {"Network Top3": (91.11, .47), "Resolution-group Top3": (70.96, .83), "Exact-region Top3": (42.02, .94)},
    "Moderate": {"Network Top3": (83.84, .93), "Resolution-group Top3": (63.14, 1.24), "Exact-region Top3": (33.24, 1.39)},
    "Severe": {"Network Top3": (71.83, 1.16), "Resolution-group Top3": (51.96, 1.02), "Exact-region Top3": (24.50, 1.17)},
    "Extreme": {"Network Top3": (58.54, 1.35), "Resolution-group Top3": (39.84, 1.70), "Exact-region Top3": (16.88, 1.34)},
}
pi_rows = []
for scenario, endpoints in s15.items():
    for endpoint, (mean, sd) in endpoints.items():
        half = multiplier * sd
        pi_rows.append({
            "scenario": scenario,
            "endpoint": endpoint,
            "repeats": 30,
            "mean_percent": mean,
            "sd_percent": sd,
            "prediction_multiplier": multiplier,
            "future_repeat_pi95_low_percent": max(0, mean - half),
            "future_repeat_pi95_high_percent": min(100, mean + half),
        })

payload = {
    "stat1_bh_family": bh,
    "stat2_mde_design": {
        "simulation_unit": "nine paired donor-rate differences",
        "icc_input": None,
        "monte_carlo_draws": 10000,
        "alpha_two_sided": 0.05,
        "target_power": 0.80,
        "reported_mde_pp": "approximately 6",
        "icc_sensitivity": icc_rows,
    },
    "stat3_prediction_interval": {
        "target": "one additional perturbation repeat under the same scenario",
        "formula": "mean +/- t_(0.975,29) * SD * sqrt(1 + 1/30)",
        "not_target": "a clinical query or a donor-population prediction interval",
        "rows": pi_rows,
    },
    "stat4_interpretation": "n=9 supports direction-consistency language only; non-significance is not robustness, equivalence, or no difference.",
}
(OUT / "P1_STAT1-4_audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

with (OUT / "P1_STAT1-4_prediction_intervals.csv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=pi_rows[0].keys())
    w.writeheader()
    w.writerows(pi_rows)

print(json.dumps(payload, indent=2))
