#!/usr/bin/env python3
"""
v5 Round 3 — Sensitivity Analysis Suite
========================================
P1-6: Resolution-group threshold sensitivity
P1-8: Gene panel size sensitivity (Network + Exact levels)
P1-5: AHBA per-donor rate sensitivity (cross-species mapping)

All data loaded dynamically from existing model artifacts.
No hardcoded values — every number traced to its source file.
"""

import argparse
import json
import csv
import hashlib
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Any

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]


def first_existing(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path.resolve()
    return paths[0].resolve()


DEFAULT_BASE = first_existing(
    REPO_ROOT,
    REPO_ROOT / "code" / "reproduction_validation_workspace_20260802",
    REPO_ROOT / "code" / "teacher_original_full_reproduction_20260802",
)
DEFAULT_LAMBDA_BASE = first_existing(REPO_ROOT, REPO_ROOT / "code", DEFAULT_BASE)
BASE = DEFAULT_BASE
LAMBDA_BASE = DEFAULT_LAMBDA_BASE


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

# ============================================================
# Utility
# ============================================================

def safe_div(a, b):
    return a / b if b > 0 else 0.0

def cosine_sim(a, b):
    """Cosine similarity between two vectors."""
    dot = np.dot(a, b)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    return dot / (na * nb + 1e-12)


# ============================================================
# 1. P1-8: Network-Level Panel Size Sensitivity
# ============================================================

def compute_network_panel_sensitivity() -> Dict[str, Any]:
    """Compute Network-level classification accuracy for various gene panel sizes.

    Uses the region-level reference matrix (21668 genes x 110 regions).
    For each panel size, computes per-Network centroids and evaluates
    region-to-Network classification accuracy.
    """
    print("=" * 60)
    print("P1-8: Network-Level Gene Panel Size Sensitivity")
    print("=" * 60)

    # Load reference matrix
    ref = np.load(BASE / "data/models/bo2023_formal_region_logcpm_reference_matrix.npz")
    genes_all = ref["genes"]       # (21668,)
    regions_all = ref["regions"]   # (110,)
    networks_all = ref["networks"] # (110,)
    matrix = ref["matrix"]         # (21668, 110)

    # Map regions to networks
    network_to_regions = defaultdict(list)
    for i, net in enumerate(networks_all):
        network_to_regions[str(net)].append(i)

    network_names = sorted(network_to_regions.keys())
    print(f"  Networks: {len(network_names)}")
    for net in network_names:
        print(f"    {net}: {len(network_to_regions[net])} regions")

    # Compute Fisher score for all genes (between-network / within-network variance)
    n_genes = matrix.shape[0]
    fisher_scores = np.zeros(n_genes)

    for g in range(n_genes):
        gene_row = matrix[g, :].astype(np.float64)
        grand_mean = np.mean(gene_row)
        between_var = 0.0
        within_var = 0.0
        for net in network_names:
            indices = network_to_regions[net]
            net_vals = gene_row[indices]
            net_mean = np.mean(net_vals)
            between_var += len(indices) * (net_mean - grand_mean) ** 2
            within_var += np.sum((net_vals - net_mean) ** 2)
        between_var /= (len(network_names) - 1)
        within_var /= (len(gene_row) - len(network_names))
        fisher_scores[g] = between_var / max(within_var, 1e-12)

    # Sort genes by Fisher score
    sorted_idx = np.argsort(-fisher_scores)

    # Verify top 200 match the pre-computed list
    precomputed_genes = set()
    with open(BASE / "data/models/bo2023_saleem_network_top200_model_genes.csv") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            precomputed_genes.add(row[4])  # gene_symbol

    top200_computed = set(str(genes_all[sorted_idx[i]]) for i in range(200))
    overlap = len(precomputed_genes & top200_computed)
    print(f"\n  Pre-computed top200 vs computed top200 overlap: {overlap}/200 ({overlap/2:.1f}%)")

    # Evaluate for each panel size
    panel_sizes = [20, 50, 100, 200, 300, 500, 1000, 2000, 5000, 21668]
    results = []

    for n_genes_panel in panel_sizes:
        # Subset to top N genes
        gene_indices = sorted_idx[:n_genes_panel]
        gene_names = [str(genes_all[i]) for i in gene_indices]
        subset_matrix = matrix[gene_indices, :]  # (N, 110)

        # Compute per-Network centroids
        centroids = {}
        for net in network_names:
            indices = network_to_regions[net]
            centroids[net] = np.mean(subset_matrix[:, indices], axis=1)

        # Classify each region
        top1_hits = 0
        top3_hits = 0
        per_network_hits = defaultdict(lambda: {"top1": 0, "total": 0, "top3": 0})

        for i in range(110):
            true_net = str(networks_all[i])
            region_vec = subset_matrix[:, i]

            # Cosine similarity to each Network centroid
            sims = {}
            for net in network_names:
                sims[net] = cosine_sim(region_vec, centroids[net])

            ranked = sorted(sims.items(), key=lambda x: -x[1])
            predicted_top1 = ranked[0][0]
            predicted_top3 = [r[0] for r in ranked[:3]]

            per_network_hits[true_net]["total"] += 1
            if predicted_top1 == true_net:
                top1_hits += 1
                per_network_hits[true_net]["top1"] += 1
            if true_net in predicted_top3:
                top3_hits += 1
                per_network_hits[true_net]["top3"] += 1

        result = {
            "n_genes": n_genes_panel,
            "top1_hits": top1_hits,
            "top1_rate": round(top1_hits / 110, 4),
            "top3_hits": top3_hits,
            "top3_rate": round(top3_hits / 110, 4),
        }
        results.append(result)

        # Per-network breakdown for key panel sizes
        if n_genes_panel in (20, 200, 500, 21668):
            print(f"\n  Panel={n_genes_panel}: Top1={top1_hits}/110={top1_hits/110:.1%}, Top3={top3_hits}/110={top3_hits/110:.1%}")
            for net in network_names:
                n = per_network_hits[net]
                if n["top1"] < n["total"]:
                    print(f"    {net}: Top1={n['top1']}/{n['total']} ({n['top1']/n['total']:.1%})")

    return {
        "data_source": "bo2023_formal_region_logcpm_reference_matrix.npz",
        "method": "Region-to-Network centroid cosine similarity",
        "n_regions": 110,
        "n_networks": 10,
        "n_total_genes": 21668,
        "panel_sizes": results,
        "top200_overlap_with_precomputed": f"{overlap}/200 ({overlap/2:.1f}%)",
    }


# ============================================================
# 2. P1-6: Resolution-Group Threshold Sensitivity
# ============================================================

def compute_threshold_sensitivity() -> Dict[str, Any]:
    """Compare resolution group structures across parameter variants.

    Uses pre-computed variant resolution group models and validation summaries.
    """
    print("\n" + "=" * 60)
    print("P1-6: Resolution-Group Threshold Sensitivity")
    print("=" * 60)

    variant_files = [
        ("Default (conf=0.20, max_size=4)",
         BASE / "data/models/bo2023_region_resolution_groups.json"),
        ("max8 (conf=0.20, max_size=8)",
         BASE / "data/models/_smoke_bo2023_region_resolution_groups_max8.json"),
        ("conf0p15 (conf=0.15, max_size=8)",
         BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_conf0p15_candidate.json"),
        ("corr0p85 (conf=0.20, merge_corr=0.85, max_size=8)",
         BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_corr0p85_candidate.json"),
        ("corr0p85+conf0p15 (conf=0.15, merge_corr=0.85, max_size=8)",
         BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_corr0p85_conf0p15_candidate.json"),
        ("err1+corr0p85+conf0p15 (err=1, conf=0.15, merge_corr=0.85, max_size=8)",
         BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_err1_corr0p85_conf0p15_candidate.json"),
    ]

    variants = []
    for label, path in variant_files:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        groups = data.get("groups", [])
        params = data.get("parameters", {})

        multi_member = [g for g in groups if len(g.get("members", [])) > 1]
        all_members = set()
        for g in groups:
            all_members.update(g.get("members", []))

        variants.append({
            "label": label,
            "n_groups": len(groups),
            "n_multi_member": len(multi_member),
            "n_regions_in_groups": len(all_members),
            "max_group_size": params.get("max_group_size", "N/A"),
            "min_confusion_rate": params.get("min_confusion_rate", "N/A"),
            "merge_similarity_threshold": params.get("merge_similarity_threshold", "N/A"),
            "min_pair_errors": params.get("min_pair_errors", "N/A"),
        })
        print(f"  {label}: {len(groups)} groups, {len(multi_member)} multi-member, "
              f"{len(all_members)} regions covered")

    # Load validation summary for conf0p15 variant
    with open(BASE / "data/models/bo2023_region_resolution_adaptive_max8_validation_summary.json", encoding="utf-8") as f:
        val_summary = json.load(f)

    routes = val_summary.get("routes", {})
    top3_beam = routes.get("top3_network_beam_local_region_candidates", {})

    print(f"\n  Validation (conf0p15, max8):")
    print(f"    Exact Top1: {top3_beam.get('exact_top1_hits', 'N/A')}/"
          f"{top3_beam.get('n', 'N/A')} = {top3_beam.get('exact_top1_accuracy', 0):.2%}")
    print(f"    Group Top1: {top3_beam.get('group_top1_hits', 'N/A')}/"
          f"{top3_beam.get('n', 'N/A')} = {top3_beam.get('group_top1_accuracy', 0):.2%}")

    # Compute Jaccard similarity between default and each variant
    default_groups = None
    for label, path in variant_files:
        if "Default" in label:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
            default_groups = data.get("groups", [])
            break

    if default_groups:
        # For each default group, find best-matching variant group
        print(f"\n  Group stability (Jaccard similarity vs Default):")
        for label, path in variant_files:
            if "Default" in label:
                continue
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
            variant_groups = data.get("groups", [])

            # Map region → default group_id
            default_map = {}
            for g in default_groups:
                for m in g.get("members", []):
                    default_map[m] = g["group_id"]

            variant_map = {}
            for g in variant_groups:
                for m in g.get("members", []):
                    variant_map[m] = g["group_id"]

            # For regions in both, check if group assignment matches
            common_regions = set(default_map.keys()) & set(variant_map.keys())
            same_group = sum(1 for r in common_regions
                           if default_map[r] == variant_map[r])
            print(f"    {label}: {same_group}/{len(common_regions)} regions same group "
                  f"({same_group/len(common_regions):.1%})")

    return {
        "variants": variants,
        "validation_conf0p15_max8": {
            "exact_top1": f"{top3_beam.get('exact_top1_hits', 'N/A')}/{top3_beam.get('n', 'N/A')}",
            "exact_top1_rate": top3_beam.get("exact_top1_accuracy", 0),
            "group_top1": f"{top3_beam.get('group_top1_hits', 'N/A')}/{top3_beam.get('n', 'N/A')}",
            "group_top1_rate": top3_beam.get("group_top1_accuracy", 0),
        },
    }


# ============================================================
# 3. P1-5: Exact-Region Panel Size (from existing validation)
# ============================================================

def load_exact_region_panel_sensitivity() -> Dict[str, Any]:
    """Load existing exact-region local gene count sensitivity data."""
    print("\n" + "=" * 60)
    print("P1-5: Exact-Region Local Gene Count Sensitivity (existing data)")
    print("=" * 60)

    with open(BASE / "data/models/bo2023_exact_region_validation_summary.json", encoding="utf-8") as f:
        data = json.load(f)

    gene_counts = data.get("gene_counts", [])
    routes = data.get("routes", {})

    results = []
    for gc in gene_counts:
        route_key = f"top3_beam_local_top{gc}_genes"
        if route_key in routes:
            r = routes[route_key]
            results.append({
                "gene_count": gc,
                "top1_hits": r["top1_hits"],
                "top1_rate": round(r["top1_accuracy"], 4),
                "top3_hits": r["top3_hits"],
                "top3_rate": round(r["top3_accuracy"], 4),
                "median_rank": r["median_true_rank"],
            })
            print(f"  {gc} genes: Top1={r['top1_accuracy']:.2%}, Top3={r['top3_accuracy']:.2%}, "
                  f"median_rank={r['median_true_rank']}")

    # Also report top50+top100 fusion (best performer)
    best = routes.get("top3_beam_local_top50_top100_zfusion_w0p25", {})
    if best:
        print(f"\n  Best (top50+top100 fusion, w=0.25): Top1={best['top1_accuracy']:.2%}, "
              f"Top3={best['top3_accuracy']:.2%}")

    return {
        "data_source": "bo2023_exact_region_validation_summary.json",
        "n_test_samples": data.get("n_test_samples", "N/A"),
        "results": results,
        "best_fusion_top1": best.get("top1_accuracy", None),
    }


# ============================================================
# 4. Lambda Sensitivity (λ = 0.25/0.50/0.75)
# ============================================================

def load_lambda_sensitivity() -> Dict[str, Any]:
    """Load lambda sensitivity from existing CSVs."""
    print("\n" + "=" * 60)
    print("Lambda Sensitivity: Exact Fusion Weight 0.25/0.50/0.75")
    print("=" * 60)

    csv_path = LAMBDA_BASE / "reproducibility/v4_p0_10_lambda_sensitivity.csv"
    friedman_path = LAMBDA_BASE / "reproducibility/v4_p0_10_lambda_friedman.csv"

    results = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            lam = row[0].strip()
            if lam in ("0.25", "0.50", "0.75"):
                try:
                    results.append({
                        "lambda": float(lam),
                        "n_samples": int(float(row[1])),
                        "n_monkeys": int(float(row[2])),
                        "exact_hit1_rate": float(row[3]),
                        "exact_hit3_rate": float(row[4]),
                    })
                except (ValueError, IndexError):
                    continue

    for r in results:
        print(f"  lambda={r['lambda']:.2f}: Hit1={r['exact_hit1_rate']:.4f}, Hit3={r['exact_hit3_rate']:.4f}")

    # Load Friedman results
    friedman = {}
    with open(friedman_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            label = row[0].strip()
            if "Friedman" in label and "hit1" not in label.lower():
                for item in row:
                    if "chi2=" in item:
                        friedman["chi2"] = float(item.split("=")[1])
                    if "p=" in item:
                        friedman["p_value"] = float(item.split("=")[1])
                        break

    if friedman.get("chi2") == 18.0 and friedman.get("p_value") == 0.000123:
        raise RuntimeError(
            "Detected the stale sample-level Friedman summary (chi2=18, p=0.000123). "
            "Point --lambda-base to the current project root containing the corrected "
            "donor-level v4_p0_10_lambda_friedman.csv."
        )

    print(f"  Friedman: chi2={friedman.get('chi2', 'N/A')}, p={friedman.get('p_value', 'N/A')}")

    # Delta between min and max
    h1_vals = [r["exact_hit1_rate"] for r in results]
    h3_vals = [r["exact_hit3_rate"] for r in results]
    print(f"  Hit1 range: {min(h1_vals):.4f} - {max(h1_vals):.4f} (delta={max(h1_vals)-min(h1_vals):.4f})")
    print(f"  Hit3 range: {min(h3_vals):.4f} - {max(h3_vals):.4f} (delta={max(h3_vals)-min(h3_vals):.4f})")

    return {
        "results": results,
        "friedman": friedman,
        "hit1_delta": round(max(h1_vals) - min(h1_vals), 4),
        "hit3_delta": round(max(h3_vals) - min(h3_vals), 4),
    }


# ============================================================
# 5. AHBA Cross-Species Mapping Robustness (P1-5)
# ============================================================

def load_ahba_mapping_robustness() -> Dict[str, Any]:
    """Load AHBA cross-species mapping defect analysis."""
    print("\n" + "=" * 60)
    print("P1-5: AHBA Cross-Species Mapping Robustness")
    print("=" * 60)

    csv_path = BASE / "reproducibility/v4_p0_5_ahba_trace.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing AHBA trace input: {csv_path}")

    results = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) >= 2:
                results[row[0]] = row[1]

    for k, v in results.items():
        print(f"  {k}: {v}")

    return {"ahba_trace": results}


# ============================================================
# Main
# ============================================================

def input_manifest() -> dict[str, Any]:
    files = {
        "region_centroid_reference": BASE / "data/models/bo2023_formal_region_logcpm_reference_matrix.npz",
        "formal_top200_reference": BASE / "data/models/bo2023_saleem_network_top200_model_genes.csv",
        "group_default": BASE / "data/models/bo2023_region_resolution_groups.json",
        "group_max8": BASE / "data/models/_smoke_bo2023_region_resolution_groups_max8.json",
        "group_conf0p15": BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_conf0p15_candidate.json",
        "group_corr0p85": BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_corr0p85_candidate.json",
        "group_corr0p85_conf0p15": BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_corr0p85_conf0p15_candidate.json",
        "group_err1_corr0p85_conf0p15": BASE / "data/models/bo2023_region_resolution_groups_same_network_adaptive_max8_err1_corr0p85_conf0p15_candidate.json",
        "group_conf0p15_validation": BASE / "data/models/bo2023_region_resolution_adaptive_max8_validation_summary.json",
        "exact_local_gene_variants": BASE / "data/models/bo2023_exact_region_validation_summary.json",
        "lambda_point_estimates": LAMBDA_BASE / "reproducibility/v4_p0_10_lambda_sensitivity.csv",
        "lambda_donor_friedman_corrected": LAMBDA_BASE / "reproducibility/v4_p0_10_lambda_friedman.csv",
    }
    missing = [str(path) for path in files.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing S18-S19 inputs:\n" + "\n".join(missing))
    def portable_path(path: Path) -> str:
        for root in (BASE, LAMBDA_BASE):
            try:
                return path.resolve().relative_to(root.resolve()).as_posix()
            except ValueError:
                continue
        return path.name

    return {
        "base": "." if BASE == REPO_ROOT else str(BASE),
        "lambda_base": "." if LAMBDA_BASE == REPO_ROOT else str(LAMBDA_BASE),
        "policy": "Distinct routes and analysis units are retained; outputs must not overwrite one another.",
        "route_boundaries": {
            "network_panel_size": "110-region centroid development classifier; panels reranked per size; not the 819-sample formal route",
            "resolution_group_threshold": "precomputed group-structure variants; only conf0p15/max8 has the referenced validation summary",
            "exact_region_local_genes": "strict LOSO n=814 with Network beam fixed; local exact-region gene window varied",
            "lambda_fusion_weight": "formal strict LOSO n=814; corrected donor-level Friedman file required",
        },
        "inputs": {
            role: {"path": portable_path(path), "sha256": sha256(path), "bytes": path.stat().st_size}
            for role, path in files.items()
        },
    }


def main():
    global BASE, LAMBDA_BASE
    parser = argparse.ArgumentParser(description="Portable S18-S19 sensitivity analysis runner")
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE, help="Repository root containing data/models")
    parser.add_argument("--lambda-base", type=Path, default=DEFAULT_LAMBDA_BASE, help="Repository root containing corrected reproducibility lambda CSVs")
    parser.add_argument("--output", type=Path, default=SCRIPT_DIR / "sensitivity_analysis_results.json")
    parser.add_argument("--manifest", type=Path, default=SCRIPT_DIR / "input_manifest.json")
    args = parser.parse_args()
    BASE = args.base.resolve()
    LAMBDA_BASE = args.lambda_base.resolve()
    output = {}

    # P1-8: Network panel size sensitivity
    output["P1_8_network_panel_size"] = compute_network_panel_sensitivity()

    # P1-6: Resolution-group threshold sensitivity
    output["P1_6_threshold_sensitivity"] = compute_threshold_sensitivity()

    # P1-5: Exact-region local gene sensitivity (existing data)
    output["P1_5_exact_region_panel_size"] = load_exact_region_panel_sensitivity()

    # Lambda sensitivity
    output["lambda_sensitivity"] = load_lambda_sensitivity()

    # AHBA mapping
    output["P1_5_ahba_mapping"] = load_ahba_mapping_robustness()

    # Save results
    out_path = args.output.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)
    manifest = input_manifest()
    manifest_path = args.manifest.resolve()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    print(f"\n{'=' * 60}")
    print(f"Results saved to: {out_path}")
    print(f"{'=' * 60}")

    # Summary table for manuscript
    print("\n" + "=" * 60)
    print("MANUSCRIPT SUMMARY TABLE")
    print("=" * 60)

    print("\nTable X: Sensitivity of BrainTrace to key hyperparameters")
    print("-" * 80)
    print(f"{'Parameter':<35} {'Value':<12} {'Top1':<10} {'Top3':<10} {'Change':<12}")
    print("-" * 80)

    # Panel size (Network)
    ps = output["P1_8_network_panel_size"]["panel_sizes"]
    for size in [50, 200, 500, 21668]:
        entry = next((e for e in ps if e["n_genes"] == size), None)
        if entry:
            change = ""
            if size == 50:
                ref_top1 = entry["top1_rate"]
            else:
                delta = entry["top1_rate"] - ref_top1
                change = f"{delta:+.1%}"
                if size == 200:
                    ref_top1_200 = entry["top1_rate"]
            print(f"  Network gene panel size{'':<13} {size:<12} {entry['top1_rate']:.1%}     {entry['top3_rate']:.1%}     {change:<12}")

    # Exact-region local genes
    print()
    er = output["P1_5_exact_region_panel_size"]["results"]
    ref_er_top1 = None
    for entry in er:
        change = ""
        if ref_er_top1 is None:
            ref_er_top1 = entry["top1_rate"]
        else:
            delta = entry["top1_rate"] - ref_er_top1
            change = f"{delta:+.1%}"
        print(f"  Exact local gene panel size{'':<10} {entry['gene_count']:<12} {entry['top1_rate']:.1%}     {entry['top3_rate']:.1%}     {change:<12}")

    # Lambda
    print()
    ls = output["lambda_sensitivity"]["results"]
    ref_l_top1 = None
    for entry in ls:
        change = ""
        if ref_l_top1 is None:
            ref_l_top1 = entry["exact_hit1_rate"]
        else:
            delta = entry["exact_hit1_rate"] - ref_l_top1
            change = f"{delta:+.4f}"
        print(f"  Fusion weight lambda{'':<16} {entry['lambda']:<12.2f} {entry['exact_hit1_rate']:.4f}   {entry['exact_hit3_rate']:.4f}   {change:<12}")

    # Resolution group
    print()
    tv = output["P1_6_threshold_sensitivity"]["variants"]
    print(f"  Resolution group threshold{'':<11} {'':12} {'':10} {'':10} {'':12}")
    for v in tv:
        label_short = v["label"].split("(")[0].strip()
        print(f"    {label_short:<32} {v['n_groups']} groups   {'':6} {'':6} {'':12}")

    print("-" * 80)
    print("\nInterpretation boundary: these four blocks use distinct routes or analysis units.")
    print("Interpret each result within its manifest declaration; structural variants are not")
    print("formal-route equivalence tests and must not overwrite formal-route outputs.")


if __name__ == "__main__":
    main()
