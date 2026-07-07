#!/usr/bin/env python
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "reports" / "p0_hard_evidence_20260629"
MODEL_DIR = ROOT / "data" / "models"
RESULT_DIR = ROOT / "results" / "bo2023_correlation_marker_routes_unseen_confirmation"


def pct(x: float) -> str:
    return f"{x:.2%}"


def write_network_marker_summary() -> pd.DataFrame:
    genes = pd.read_csv(MODEL_DIR / "bo2023_saleem_network_top200_model_genes.csv")
    out = pd.DataFrame(
        [
            {
                "marker_set": "Network global discriminative genes",
                "n_genes": int(len(genes)),
                "selection_score": "Fisher-like between-Network / within-Network variance ratio",
                "fisher_score_min": float(genes["fisher_score"].min()),
                "fisher_score_median": float(genes["fisher_score"].median()),
                "fisher_score_max": float(genes["fisher_score"].max()),
                "between_variance_median": float(genes["between_variance"].median()),
                "within_variance_median": float(genes["within_variance"].median()),
                "source_file": "data/models/bo2023_saleem_network_top200_model_genes.csv",
            }
        ]
    )
    out.to_csv(OUTDIR / "marker_network_top200_summary.csv", index=False)
    genes.to_csv(OUTDIR / "marker_network_top200_genes.csv", index=False)
    return out


def write_pairwise_marker_summary() -> pd.DataFrame:
    data = json.loads((MODEL_DIR / "bo2023_saleem_network_pairwise_rescue_model.json").read_text(encoding="utf-8"))
    rows = []
    for pair in data.get("pairs", []):
        rows.append(
            {
                "pair_key": pair.get("key"),
                "left_network": pair.get("left_network"),
                "right_network": pair.get("right_network"),
                "n_pair_genes": len(pair.get("genes", [])),
                "first_10_genes": " | ".join(pair.get("genes", [])[:10]),
                "selection_scope": "high-confusion Network pair; Top1 rescue constrained to original Top3 beam",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "marker_pairwise_rescue_summary.csv", index=False)
    return out


def write_marker_methodology_audit() -> pd.DataFrame:
    pair_model = json.loads((MODEL_DIR / "bo2023_saleem_network_pairwise_rescue_model.json").read_text(encoding="utf-8"))
    marker_validation = json.loads((RESULT_DIR / "validation_summary.json").read_text(encoding="utf-8"))
    coverage = pd.read_csv(RESULT_DIR / "marker_coverage_summary.csv")
    cov = dict(zip(coverage["metric"], coverage["value"]))
    rows = [
        {
            "component": "Network global marker panel",
            "stage": "Network Top3 beam generation",
            "selection_rule": "Rank genes by Fisher-like between-Network / within-Network variance ratio on Bo2023 training data; retain Top 200.",
            "n_markers_or_features": 200,
            "fold_policy": "Validation scripts rebuild references fold-locally; production model stores locked Top200 marker order.",
            "cross_species_mapping": "Gene symbols are harmonized through cleaned macaque gene-symbol mapping before model fitting; AHBA external validation uses mapped human labels and shared gene symbols where present.",
            "cfRNA_degradation_tolerance": "Not empirically proven on anatomical-truth plasma cfRNA; handled as marker/gene overlap and coverage diagnostics, not as a claimed degradation-resistance validation.",
            "evidence_file": "marker_network_top200_genes.csv; marker_network_top200_summary.csv",
        },
        {
            "component": "Network pairwise rescue markers",
            "stage": "Top1 reordering within retained Network Top3 beam",
            "selection_rule": "Identify high-confusion Network pairs in training predictions; for each pair retain pair-specific discriminative genes.",
            "n_markers_or_features": int(pair_model["parameters"]["pair_top_n_genes"]),
            "fold_policy": "LOMO formal validation builds pair models fold-locally; production pairwise model records locked full-reference pair panels.",
            "cross_species_mapping": "Same cleaned gene-symbol space as Network model.",
            "cfRNA_degradation_tolerance": "No standalone cfRNA degradation claim; output remains constrained to Top3 beam and is accompanied by score-margin/coverage diagnostics.",
            "evidence_file": "marker_pairwise_rescue_summary.csv",
        },
        {
            "component": "Resolution-group local discriminative genes",
            "stage": "Resolution-group reranking inside Network Top3 beam",
            "selection_rule": "Within each fold/sample candidate set, rank candidate regions/groups by local discriminative score; use local Top 200 genes for group scoring.",
            "n_markers_or_features": 200,
            "fold_policy": "Fold-local/sample-local; held-out truth sample is excluded from training reference.",
            "cross_species_mapping": "Internal Bo2023 region hierarchy only; AHBA results reported only at mapped-label supported levels.",
            "cfRNA_degradation_tolerance": "No direct cfRNA validation; low marker coverage should trigger coarse/low-confidence interpretation.",
            "evidence_file": "formal_loso_rerun_20260629/*resolution_group_detail.csv; formal_three_tier_lomo_hybrid/*resolution_group_detail.csv",
        },
        {
            "component": "Exact-region local discriminative genes",
            "stage": "Exploratory exact-region reranking inside Network Top3 beam",
            "selection_rule": "Within candidate regions, compute local discriminative gene order; fuse Top50 and Top100 correlation z-scores with fixed exact-fusion weight 0.25.",
            "n_markers_or_features": "Top50/Top100 fusion from local gene order",
            "fold_policy": "Fold-local/sample-local; exact-region endpoint is exploratory.",
            "cross_species_mapping": "Internal Bo2023 exact regions; AHBA exact metrics only for exact-evaluable mapped subset.",
            "cfRNA_degradation_tolerance": "Not established; exact-region output should not be interpreted as cfRNA localization accuracy.",
            "evidence_file": "formal_loso_rerun_20260629/*exact_region_detail.csv; formal_three_tier_lomo_hybrid/*exact_region_detail.csv",
        },
        {
            "component": "Development marker annotation route",
            "stage": "Marker support/annotation only, not adopted as main reranking route",
            "selection_rule": f"Stable markers required min_consistency={marker_validation['min_consistency']}, min_effect={marker_validation['min_effect']}, min_markers_for_support={marker_validation['min_markers_for_support']}; Top10 marker rerank was tested but not adopted.",
            "n_markers_or_features": f"stable_topk_per_region={marker_validation['stable_topk_per_region']}; marker_pairs_per_fold_mean={cov.get('marker_pairs_per_fold_mean')}",
            "fold_policy": f"Strict unseen-sample check, n={marker_validation['n_test_samples']}, seed={marker_validation['seed']}; prior samples excluded={marker_validation['n_prior_samples_excluded']}.",
            "cross_species_mapping": "Bo2023 internal marker annotation route; not a cross-species marker validation.",
            "cfRNA_degradation_tolerance": "Useful as coverage/support diagnostic only; route decision says marker rerank did not replace correlation primary.",
            "evidence_file": "results/bo2023_correlation_marker_routes_unseen_confirmation/",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "marker_methodology_audit.csv", index=False)
    return out


def write_report(network: pd.DataFrame, pairwise: pd.DataFrame, audit: pd.DataFrame) -> None:
    net = network.iloc[0]
    pair_genes = pairwise["n_pair_genes"].describe()
    text = f"""# P0-5 Marker methodology evidence

This note completes the P0-5 reviewer request: marker selection criteria, marker counts, cross-species mapping policy and cfRNA degradation-tolerance boundaries.

## Marker panels and selection rules

1. **Network global marker panel.** The locked Network model uses the Top 200 genes ranked by a Fisher-like between-Network / within-Network variance ratio. The exported Top200 panel has median Fisher score {net['fisher_score_median']:.4f}, minimum {net['fisher_score_min']:.4f} and maximum {net['fisher_score_max']:.4f}. Full gene list: `marker_network_top200_genes.csv`.

2. **Pairwise Network rescue markers.** The pairwise rescue model contains {len(pairwise)} high-confusion Network pairs. Each pair stores pair-specific discriminative genes; median pair panel size is {pair_genes['50%']:.0f} genes. Rescue can reorder Top1 only within the original Network Top3 beam and cannot introduce a new Network outside the beam.

3. **Resolution-group local markers.** Resolution-group reranking uses fold-local/sample-local discriminative genes within the retained Network Top3 candidate set. The scoring route uses local Top 200 genes. Held-out samples are excluded from the fold training reference.

4. **Exact-region local markers.** Exact-region reranking is exploratory. It uses a local gene order inside the retained Network beam and fuses Top50 and Top100 correlation z-scores with fixed weight 0.25. This endpoint should not be described as deterministic localization.

5. **Development marker annotation route.** A separate stable-marker annotation/rerank experiment required region-marker consistency and effect-size thresholds, but it did not improve over the primary correlation route and was not adopted as the main route. It remains evidence for marker-support diagnostics.

## Cross-species gene mapping

The model uses cleaned macaque gene-symbol mapping before fitting. AHBA external validation is interpreted as mapped-label transfer and uses shared gene symbols where present. Cross-species label mapping is not treated as direct human-macaque exact anatomical equivalence.

## cfRNA degradation tolerance

No current file establishes true plasma cfRNA degradation tolerance against patient-level anatomical truth. The correct claim boundary is: marker overlap, coverage, score margin and entropy are diagnostics for reliability and domain transfer; they do not prove cfRNA degradation-resistant localization. Manuscript text should explicitly state this limitation.

## Evidence files

- `marker_methodology_audit.csv`
- `marker_network_top200_summary.csv`
- `marker_network_top200_genes.csv`
- `marker_pairwise_rescue_summary.csv`
- `results/bo2023_correlation_marker_routes_unseen_confirmation/marker_coverage_summary.csv`
- `results/bo2023_correlation_marker_routes_unseen_confirmation/validation_summary.json`

## Suggested manuscript wording

Marker selection was performed in the training reference only. For Network-level scoring, genes were ranked by a Fisher-like ratio of between-Network to within-Network variance, and the locked model retained the top 200 genes. Pairwise rescue used separate pair-specific panels for high-confusion Network pairs, but rescue was constrained to reorder Top1 within the original Top3 Network beam. Downstream resolution-group and exact-region reranking used fold-local discriminative gene orders within the retained Network beam; exact-region scoring fused Top50 and Top100 local-correlation z-scores and was treated as exploratory. Marker coverage and score diagnostics were reported to flag sparse or domain-shifted profiles. Because no anatomical-truth plasma cfRNA cohort is available, these marker panels should not be interpreted as empirically validated cfRNA degradation-tolerant markers.
"""
    (OUTDIR / "MARKER_METHODOLOGY_REPORT.md").write_text(text, encoding="utf-8")


def update_completion_checklist() -> None:
    path = OUTDIR / "P0_COMPLETION_CHECKLIST.md"
    if not path.exists():
        return
    text = path.read_text(encoding="utf-8")
    line = "- [x] Marker selection criteria, counts, cross-species mapping policy and cfRNA degradation boundary - marker_methodology_audit.csv; MARKER_METHODOLOGY_REPORT.md"
    if "Marker selection criteria" not in text:
        text = text.rstrip() + "\n" + line + "\n"
        path.write_text(text, encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    network = write_network_marker_summary()
    pairwise = write_pairwise_marker_summary()
    audit = write_marker_methodology_audit()
    write_report(network, pairwise, audit)
    update_completion_checklist()
    print(f"Wrote marker methodology evidence to {OUTDIR}")


if __name__ == "__main__":
    main()
