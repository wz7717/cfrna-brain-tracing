"""
P2-7: Enhanced beta-binomial CR2 coverage simulation.
Adds p=0.56 grid cell, uses n_sim=5000.
"""
import numpy as np
from pathlib import Path
import json, csv

SEED = 20260805
N_SIM = 5000
N_BOOT = 999
DONOR_SIZES = [17, 219, 11, 17, 114, 79, 222, 98, 42]  # 9 donors, total 819
G = len(DONOR_SIZES)
TOTAL_N = sum(DONOR_SIZES)

rng = np.random.default_rng(SEED)
outdir = Path("D:/Download/文章改稿/github_main_sync/reproducibility/round5_analysis/p2_7_beta_binomial_enhanced")
outdir.mkdir(parents=True, exist_ok=True)

def cr2_interval(hits_per_donor, n_per_donor, df_satt=4.27):
    """Fast bias-reduced CR2 formula: rate +/- t_{df,0.975} * SE_CR2."""
    rates = hits_per_donor / n_per_donor
    n_total = sum(n_per_donor)
    theta_hat = sum(hits_per_donor) / n_total
    # CR2 variance: sum((n_i*(rate_i - theta_hat))^2) / (sum(n_i))^2 * G/(G-1)
    # But this is simplified. Use the actual CR2 formula:
    # SE_CR2 = sqrt( G/(G-1) * sum_i (w_i * (rate_i - theta_hat))^2 )
    # where w_i = n_i / sum(n_i)
    w = n_per_donor / n_total
    residuals = rates - theta_hat
    se_cr2 = np.sqrt(G / (G - 1) * np.sum((w * residuals) ** 2))
    t_val = 2.776  # t_{4.27, 0.975}
    lo = theta_hat - t_val * se_cr2
    hi = theta_hat + t_val * se_cr2
    return lo, hi

def percentile_bootstrap_ci(hits, n_vec, n_boot):
    g = len(n_vec)
    rates = np.full(n_boot, np.nan)
    for b in range(n_boot):
        idx = rng.integers(0, g, size=g)
        rates[b] = np.sum(hits[idx]) / np.sum(n_vec[idx])
    return np.percentile(rates, [2.5, 97.5])

def cr1_t8_interval(hits, n_vec):
    """CR1/G(G-1) with t_8 df."""
    rates = hits / n_vec
    n_total = np.sum(n_vec)
    theta = np.sum(hits) / n_total
    w = n_vec / n_total
    se_cr1 = np.sqrt(1 / (G - 1) * np.sum(w * (rates - theta) ** 2))
    t8 = 2.306
    return theta - t8 * se_cr1, theta + t8 * se_cr1

def beta_binomial_counts(n_vec, p, rho):
    shape = (1 - rho) / rho
    theta = rng.beta(p * shape, (1 - p) * shape, size=len(n_vec))
    return rng.binomial(n_vec, theta)

# Grid: p × ICC
p_grid = [0.45, 0.56, 0.72, 0.92]
icc_grid = [0.05, 0.10, 0.20, 0.30]
donor_sizes = np.array(DONOR_SIZES)

print(f"Running {len(p_grid) * len(icc_grid)} cells x {N_SIM} sims...")
results_rows = []

for cell_i, (p, rho) in enumerate([(p, r) for p in p_grid for r in icc_grid]):
    print(f"  Cell {cell_i+1}/{len(p_grid)*len(icc_grid)}: p={p}, ICC={rho}...")

    cov_pct = cov_cr1 = cov_cr2 = 0
    wid_pct = wid_cr1 = wid_cr2 = 0.0

    for s in range(N_SIM):
        hits = beta_binomial_counts(donor_sizes, p, rho)

        # Percentile bootstrap
        pct_lo, pct_hi = percentile_bootstrap_ci(hits, donor_sizes, N_BOOT)
        if pct_lo <= p <= pct_hi: cov_pct += 1
        wid_pct += pct_hi - pct_lo

        # CR1/t8
        cr1_lo, cr1_hi = cr1_t8_interval(hits, donor_sizes)
        if cr1_lo <= p <= cr1_hi: cov_cr1 += 1
        wid_cr1 += cr1_hi - cr1_lo

        # CR2/Satterthwaite
        cr2_lo, cr2_hi = cr2_interval(hits, donor_sizes)
        if cr2_lo <= p <= cr2_hi: cov_cr2 += 1
        wid_cr2 += cr2_hi - cr2_lo

    rate = cov_cr2 / N_SIM
    mcse = np.sqrt(rate * (1 - rate) / N_SIM)

    row = {
        "p": p, "ICC": rho,
        "n_sim": N_SIM, "n_donors": G, "total_n": TOTAL_N,
        "pct_cov": cov_pct / N_SIM, "pct_width": wid_pct / N_SIM,
        "cr1_cov": cov_cr1 / N_SIM, "cr1_width": wid_cr1 / N_SIM,
        "cr2_cov": cov_cr2 / N_SIM, "cr2_width": wid_cr2 / N_SIM,
        "cr2_mcse": mcse,
    }
    results_rows.append(row)

    print(f"    CR2 cov={row['cr2_cov']:.4f} (±{mcse:.6f}), Pct={row['pct_cov']:.4f}, CR1={row['cr1_cov']:.4f}")

# Save
keys = ["p","ICC","n_sim","n_donors","total_n","pct_cov","pct_width","cr1_cov","cr1_width","cr2_cov","cr2_width","cr2_mcse"]
with open(outdir / "beta_binomial_coverage_enhanced.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    w.writerows(results_rows)

# Summary
summary = {
    "description": "Enhanced beta-binomial CR2 coverage simulation (P2-7)",
    "grid": {"p": p_grid, "ICC": icc_grid},
    "n_sim": N_SIM, "n_boot": N_BOOT, "seed": SEED,
    "donor_sizes": DONOR_SIZES,
    "cr2_coverage_range": f"{min(r['cr2_cov'] for r in results_rows):.4f} - {max(r['cr2_cov'] for r in results_rows):.4f}",
    "worst_case": min(results_rows, key=lambda r: r['cr2_cov']),
}
with open(outdir / "summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nSaved to {outdir}")
print(f"CR2 range: {summary['cr2_coverage_range']}")
worst = summary['worst_case']
print(f"Worst case: p={worst['p']}, ICC={worst['ICC']}, CR2={worst['cr2_cov']:.4f} ± {worst['cr2_mcse']:.6f}")
