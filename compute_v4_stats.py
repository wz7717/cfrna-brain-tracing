#!/usr/bin/env python
"""Compute statistics needed for v4 P0 fixes.

P0-9:  Wilson / Clopper-Pearson / Agresti-Coull triple CI
P0-10: Lambda sensitivity summary (3 points available: 0.25, 0.50, 0.75)
P0-11: RF fair comparator summary
P0-12: TCGA/BraTS bootstrap CI
P0-13: macro-F1 and class imbalance analysis
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist

ROOT = Path(__file__).resolve().parent

# ─── P0-9: Triple CI (Wilson / Clopper-Pearson / Agresti-Coull) ───────────

def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return centre - half, centre + half

def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    if k == 0:
        lo = 0.0
    else:
        lo = float(beta_dist.ppf(alpha / 2, k, n - k + 1))
    if k == n:
        hi = 1.0
    else:
        hi = float(beta_dist.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi

def agresti_coull_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    n_tilde = n + z * z
    p_tilde = (k + z * z / 2) / n_tilde
    half = z * math.sqrt(p_tilde * (1 - p_tilde) / n_tilde)
    return max(0.0, p_tilde - half), min(1.0, p_tilde + half)

print("=" * 80)
print("P0-9: Wilson / Clopper-Pearson / Agresti-Coull Triple CI")
print("=" * 80)

donor_csv = ROOT / "reports/validation_recheck_20260713_canonical110/donor_cluster_inference/donor_clustered_loso_lomo_summary.csv"
donor_df = pd.read_csv(donor_csv)

# Extract k/n for each endpoint and method
ci_rows = []
for _, row in donor_df.iterrows():
    endpoint = row["endpoint"]
    n = int(row["paired_n"])
    for method, prefix in [("LOSO", "loso"), ("LOMO", "lomo")]:
        acc = float(row[f"{prefix}_sample_weighted_accuracy"])
        k = int(round(acc * n))
        w_lo, w_hi = wilson_ci(k, n)
        cp_lo, cp_hi = clopper_pearson_ci(k, n)
        ac_lo, ac_hi = agresti_coull_ci(k, n)
        ci_rows.append({
            "endpoint": endpoint,
            "method": method,
            "n": n,
            "hits": k,
            "accuracy": acc,
            "wilson_low": w_lo, "wilson_high": w_hi,
            "clopper_pearson_low": cp_lo, "clopper_pearson_high": cp_hi,
            "agresti_coull_low": ac_lo, "agresti_coull_high": ac_hi,
            "cp_minus_wilson_width": (cp_hi - cp_lo) - (w_hi - w_lo),
            "ac_minus_wilson_width": (ac_hi - ac_lo) - (w_hi - w_lo),
        })

ci_df = pd.DataFrame(ci_rows)
print(ci_df.to_string(index=False, float_format="%.6f"))
ci_df.to_csv(ROOT / "v4_p0_9_triple_ci.csv", index=False)
print(f"\nSaved to v4_p0_9_triple_ci.csv")

# ─── P0-13: macro-F1 and class imbalance ───────────────────────────────────

print("\n" + "=" * 80)
print("P0-13: macro-F1 and class imbalance analysis")
print("=" * 80)

from sklearn.metrics import f1_score

def compute_macro_f1(detail_path: str, truth_col: str, pred_col: str, hit_col: str, label: str) -> dict:
    df = pd.read_csv(detail_path)
    if truth_col not in df.columns or pred_col not in df.columns:
        return {"label": label, "error": f"missing columns: {truth_col} or {pred_col}"}
    truth = df[truth_col].astype(str)
    pred = df[pred_col].astype(str)
    classes = sorted(set(truth.unique()) | set(pred.unique()))
    macro_f1 = f1_score(truth, pred, average="macro", labels=classes, zero_division=0)
    weighted_f1 = f1_score(truth, pred, average="weighted", labels=classes, zero_division=0)
    # Per-class F1
    per_class = {}
    for cls in classes:
        tp = int(((truth == cls) & (pred == cls)).sum())
        fp = int(((truth != cls) & (pred == cls)).sum())
        fn = int(((truth == cls) & (pred != cls)).sum())
        n_cls = int((truth == cls).sum())
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        per_class[cls] = {"n": n_cls, "precision": precision, "recall": recall, "f1": f1}
    return {
        "label": label,
        "n_samples": len(df),
        "n_classes": len(classes),
        "macro_f1": macro_f1,
        "weighted_f1": weighted_f1,
        "per_class": per_class,
    }

# LOSO Network (columns: label, pred_top1, hit1, hit3)
loso_net = compute_macro_f1(
    str(ROOT / "reports/validation_recheck_20260713_canonical110/loso/hybrid_formal_loso_network_detail.csv"),
    "label", "pred_top1", "hit1", "LOSO Network Top1"
)
# LOSO Exact Region (columns: label, pred_top1, hit1, hit3)
loso_exact = compute_macro_f1(
    str(ROOT / "reports/validation_recheck_20260713_canonical110/loso/hybrid_formal_loso_exact_region_detail.csv"),
    "label", "pred_top1", "hit1", "LOSO Exact Region Top1"
)
# LOMO Network
lomo_net = compute_macro_f1(
    str(ROOT / "reproducibility/p2_publication_completeness/formal_lomo_network_detail.csv"),
    "truth", "pred_top1", "hit1", "LOMO Network Top1"
)
# LOMO Exact Region
lomo_exact = compute_macro_f1(
    str(ROOT / "reports/validation_recheck_20260713_canonical110/lomo/formal_lomo_exact_region_detail.csv"),
    "label", "pred_top1", "hit1", "LOMO Exact Region Top1"
)

for result in [loso_net, loso_exact, lomo_net, lomo_exact]:
    if "error" in result:
        print(f"\n{result['label']}: ERROR - {result['error']}")
        continue
    print(f"\n{result['label']}:")
    print(f"  macro-F1 = {result['macro_f1']:.4f}, weighted-F1 = {result['weighted_f1']:.4f}")
    print(f"  n_classes = {result['n_classes']}, n_samples = {result['n_samples']}")
    # Show smallest classes
    pc = result["per_class"]
    sorted_classes = sorted(pc.items(), key=lambda x: x[1]["n"])
    print("  Smallest 5 classes:")
    for cls, info in sorted_classes[:5]:
        print(f"    {cls}: n={info['n']}, F1={info['f1']:.4f}, P={info['precision']:.4f}, R={info['recall']:.4f}")
    print("  Largest 5 classes:")
    for cls, info in sorted_classes[-5:]:
        print(f"    {cls}: n={info['n']}, F1={info['f1']:.4f}, P={info['precision']:.4f}, R={info['recall']:.4f}")

# ─── P0-12: TCGA/BraTS bootstrap CI ────────────────────────────────────────

print("\n" + "=" * 80)
print("P0-12: TCGA/BraTS bootstrap CI")
print("=" * 80)

tcga_path = ROOT / "results/archived_tcga_brats/tcga_gbm_lgg_sample_mri_label_tracing_20260605/tcga_gbm_lgg_sample_mri_label_tracing_summary.csv"
tcga_df = pd.read_csv(tcga_path)
# Filter to only samples with MRI labels (has_mri_label=True)
labeled = tcga_df[tcga_df["has_mri_label"] == True].copy()
n_labeled = len(labeled)
print(f"Total TCGA samples: {len(tcga_df)}, with MRI labels: {n_labeled}")

# Compute accuracy and Wilson/CP/AC CIs for each metric
for metric, match_col in [
    ("Network Top1", "network_top1_match"),
    ("Network Top3", "network_top3_match"),
    ("Broad anatomy Top1", "broad_top1_match"),
    ("Broad anatomy Top3", "broad_top3_match"),
    ("Exact region Top1", "exact_region_top1_match"),
    ("Exact region Top5", "exact_region_top5_match"),
]:
    vals = pd.to_numeric(labeled[match_col], errors="coerce").dropna().astype(int)
    n = len(vals)
    k = int(vals.sum())
    if n == 0:
        print(f"  {metric}: no labeled samples")
        continue
    w_lo, w_hi = wilson_ci(k, n)
    cp_lo, cp_hi = clopper_pearson_ci(k, n)
    ac_lo, ac_hi = agresti_coull_ci(k, n)
    print(f"  {metric}: {k}/{n} = {k/n:.4f}, Wilson=[{w_lo:.4f}, {w_hi:.4f}], CP=[{cp_lo:.4f}, {cp_hi:.4f}], AC=[{ac_lo:.4f}, {ac_hi:.4f}]")

# ─── P0-10: Lambda sensitivity ─────────────────────────────────────────────

print("\n" + "=" * 80)
print("P0-10: Lambda sensitivity summary")
print("=" * 80)

lambda_results = []
for lam_dir, lam_label in [
    ("loso", "0.25 (default)"),
    ("loso_lambda_0p50", "0.50"),
    ("loso_lambda_0p75", "0.75"),
]:
    net_path = ROOT / f"reports/validation_recheck_20260713_canonical110/{lam_dir}/hybrid_formal_loso_network_detail.csv"
    exact_path = ROOT / f"reports/validation_recheck_20260713_canonical110/{lam_dir}/hybrid_formal_loso_exact_region_detail.csv"
    if not net_path.exists():
        print(f"  Lambda {lam_label}: network detail not found at {net_path}")
        continue
    net_df = pd.read_csv(net_path)
    n = len(net_df)
    net_top1 = float(net_df["hit1"].sum()) / n if "hit1" in net_df.columns else float("nan")
    net_top3 = float(net_df["hit3"].sum()) / n if "hit3" in net_df.columns else float("nan")

    exact_top1 = float("nan")
    exact_top3 = float("nan")
    exact_n = 0
    if exact_path.exists():
        exact_df = pd.read_csv(exact_path)
        exact_n = len(exact_df)
        if "hit1" in exact_df.columns:
            exact_top1 = float(exact_df["hit1"].sum()) / exact_n
        if "hit3" in exact_df.columns:
            exact_top3 = float(exact_df["hit3"].sum()) / exact_n

    lambda_results.append({
        "lambda": lam_label,
        "network_top1": net_top1,
        "network_top3": net_top3,
        "network_n": n,
        "exact_top1": exact_top1,
        "exact_top3": exact_top3,
        "exact_n": exact_n,
    })

lambda_df = pd.DataFrame(lambda_results)
print(lambda_df.to_string(index=False, float_format="%.4f"))
lambda_df.to_csv(ROOT / "v4_p0_10_lambda_sensitivity.csv", index=False)

# ─── P0-11: RF fair comparator summary ─────────────────────────────────────

print("\n" + "=" * 80)
print("P0-11: RF fair comparator summary")
print("=" * 80)

rf_pred_path = ROOT / "validation_runs/r08_rf_fair_comparator_20260717/full_outputs/full_predictions.csv"
rf_df = pd.read_csv(rf_pred_path)
print(f"RF predictions columns: {rf_df.columns.tolist()}")
print(f"RF predictions shape: {rf_df.shape}")

# Compute RF Top1/Top3 for Network and Exact
for endpoint_name in ["Network", "Exact"]:
    subset = rf_df[rf_df["endpoint"] == endpoint_name]
    if subset.empty:
        continue
    n = len(subset)
    top1_hits = int(subset["hit1"].sum())
    top3_hits = int(subset["hit3"].sum())
    w1_lo, w1_hi = wilson_ci(top1_hits, n)
    w3_lo, w3_hi = wilson_ci(top3_hits, n)
    cp1_lo, cp1_hi = clopper_pearson_ci(top1_hits, n)
    cp3_lo, cp3_hi = clopper_pearson_ci(top3_hits, n)
    print(f"\n  RF {endpoint_name}: n={n}, Top1={top1_hits}/{n}={top1_hits/n:.4f} (Wilson: [{w1_lo:.4f}, {w1_hi:.4f}]), Top3={top3_hits}/{n}={top3_hits/n:.4f} (Wilson: [{w3_lo:.4f}, {w3_hi:.4f}])")

    # Also compute macro-F1 for RF
    truth = subset["truth"].astype(str)
    pred = subset["pred_top1"].astype(str)
    classes = sorted(set(truth.unique()) | set(pred.unique()))
    macro_f1 = f1_score(truth, pred, average="macro", labels=classes, zero_division=0)
    weighted_f1 = f1_score(truth, pred, average="weighted", labels=classes, zero_division=0)
    print(f"    macro-F1={macro_f1:.4f}, weighted-F1={weighted_f1:.4f}, n_classes={len(classes)}")

print("\nRF contract key fields:")
contract_path = ROOT / "validation_runs/r08_rf_fair_comparator_20260717/full_outputs/full_contract.json"
with open(contract_path) as f:
    contract = json.load(f)
# Print results summary if available
for key in contract:
    if "result" in key.lower() or "summary" in key.lower() or "metric" in key.lower():
        print(f"  {key}: {contract[key]}")

print("\n✅ All computations complete.")
