"""Derive the P1 E4 per-Network AHBA mapping-impact summary.

The primary AHBA Network endpoint is the set-valued, any-allowed-label Top1
hit.  Because AHBA truth can contain more than one allowed Network, the same
sample may contribute to more than one Network-specific denominator.  A
strict same-Network Top1 sensitivity is reported alongside it.  Neither
metric is treated as a donor-level inferential estimate (AHBA has two
donors); the percentage-point effect is defined against the pooled mapped
sample Top1 hit rate (165/223 = 73.99%).

The public repository stores the derived CSV.  The default input paths point
to the local validation workspace used for the revision package; callers can
provide equivalent files with ``--humanization`` and ``--sample-detail``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


NETWORKS = [
    "Cingulate gyrus",
    "Frontal (agranular frontal motor areas)",
    "Hippocampal formation",
    "Lateral Prefrontal Cortex",
    "Occipital/Temporal",
    "Operculum/Insula",
    "Orbitomedial Prefrontal Cortex (OMPFC)",
    "Parietal, and Parieto-occipital region",
    "Subcortical",
    "Temporal",
]

HIGH_DEFICIT = {
    "Cingulate gyrus",
    "Lateral Prefrontal Cortex",
    "Operculum/Insula",
    "Orbitomedial Prefrontal Cortex (OMPFC)",
    "Parietal, and Parieto-occipital region",
}

WORKSPACE = Path(__file__).resolve().parents[3]
DEFAULT_HUMANIZATION = (
    WORKSPACE
    / "code"
    / "reproduction_validation_workspace_20260802"
    / "validation_runs"
    / "r08_high_feasibility_20260717"
    / "outputs"
    / "humanization"
    / "network_humanization_enrichment.csv"
)
DEFAULT_SAMPLE_DETAIL = (
    WORKSPACE
    / "code"
    / "reproducibility"
    / "p0_bio3_ahba_formal"
    / "ahba_formal_three_tier_sample_detail.csv"
)
DEFAULT_OUTPUT = (
    WORKSPACE
    / "github_main_sync"
    / "reproducibility"
    / "p1_e4_ahba_network_mapping_impact.csv"
)
DEFAULT_MANIFEST = (
    WORKSPACE
    / "github_main_sync"
    / "reproducibility"
    / "p1_e4_ahba_network_mapping_impact_manifest.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_labels(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {item.strip() for item in str(value).split(" | ") if item.strip()}


def derive(humanization: Path, sample_detail: Path) -> tuple[pd.DataFrame, dict]:
    mapping = pd.read_csv(humanization)
    required_mapping = {
        "unit",
        "signature_rows",
        "unmapped_rows",
        "unmapped_row_fraction",
    }
    missing = required_mapping - set(mapping.columns)
    if missing:
        raise ValueError(f"humanization CSV is missing columns: {sorted(missing)}")
    mapping = mapping.set_index("unit")
    missing_networks = set(NETWORKS) - set(mapping.index)
    if missing_networks:
        raise ValueError(f"humanization CSV is missing Networks: {sorted(missing_networks)}")

    detail = pd.read_csv(sample_detail)
    required_detail = {
        "route",
        "supported_for_accuracy",
        "allowed_bo2023_networks",
        "network_top1",
        "network_top1_hit",
    }
    missing = required_detail - set(detail.columns)
    if missing:
        raise ValueError(f"AHBA sample-detail CSV is missing columns: {sorted(missing)}")
    eligible = detail[
        (detail["route"] == "hybrid_projected_network_logcpm_exact")
        & detail["supported_for_accuracy"].astype(bool)
    ].copy()
    if len(eligible) != 223:
        raise ValueError(f"expected 223 Network-evaluable AHBA samples, found {len(eligible)}")

    pooled_hits = int(eligible["network_top1_hit"].astype(bool).sum())
    pooled_n = int(len(eligible))
    pooled_accuracy = pooled_hits / pooled_n
    rows: list[dict[str, object]] = []
    for network in NETWORKS:
        truth_mask = eligible["allowed_bo2023_networks"].map(split_labels).map(
            lambda labels: network in labels
        )
        subset = eligible[truth_mask]
        n = int(len(subset))
        any_hits = int(subset["network_top1_hit"].astype(bool).sum())
        strict_hits = int((subset["network_top1"].astype(str) == network).sum())
        rows.append(
            {
                "network": network,
                "mapping_group": "deficit-heavy" if network in HIGH_DEFICIT else "lower-deficit",
                "signature_rows": int(mapping.loc[network, "signature_rows"]),
                "unmapped_rows": int(mapping.loc[network, "unmapped_rows"]),
                "unmapped_fraction": float(mapping.loc[network, "unmapped_row_fraction"]),
                "ahba_allowed_network_n": n,
                "ahba_any_allowed_top1_hits": any_hits,
                "ahba_any_allowed_top1_accuracy": any_hits / n if n else np.nan,
                "ahba_strict_network_top1_hits": strict_hits,
                "ahba_strict_network_top1_accuracy": strict_hits / n if n else np.nan,
                "effect_delta_vs_pooled_any_allowed_top1_pp": (
                    100 * (any_hits / n - pooled_accuracy) if n else np.nan
                ),
                "effect_delta_vs_pooled_strict_network_top1_pp": (
                    100 * (strict_hits / n - pooled_accuracy) if n else np.nan
                ),
            }
        )
    result = pd.DataFrame(rows)

    valid = result[result["ahba_allowed_network_n"] > 0]
    rho_any = float(
        np.corrcoef(
            valid["unmapped_fraction"].to_numpy(float),
            valid["ahba_any_allowed_top1_accuracy"].to_numpy(float),
        )[0, 1]
    )
    rho_strict = float(
        np.corrcoef(
            valid["unmapped_fraction"].to_numpy(float),
            valid["ahba_strict_network_top1_accuracy"].to_numpy(float),
        )[0, 1]
    )
    manifest = {
        "analysis": "P1 E4 per-Network AHBA mapping-impact audit",
        "primary_endpoint": "set-valued any-allowed-label Network Top1 hit",
        "sensitivity_endpoint": "strict same-Network Network Top1 hit",
        "truth_denominator": "AHBA samples whose allowed Network set contains the named Network; denominators overlap under multi-label truth",
        "pooled_any_allowed_top1": {"hits": pooled_hits, "n": pooled_n, "accuracy": pooled_accuracy},
        "n_networks_with_ahba_support": int((result["ahba_allowed_network_n"] > 0).sum()),
        "n_networks_without_ahba_support": int((result["ahba_allowed_network_n"] == 0).sum()),
        "pearson_rho_unmapped_fraction_vs_any_allowed_accuracy": rho_any,
        "pearson_rho_unmapped_fraction_vs_strict_accuracy": rho_strict,
        "multi_label_truth_note": "These are descriptive transfer metrics from two AHBA donors; no Network-level inferential P value or causal mapping-bias claim is made.",
        "source_sha256": {
            "network_humanization_enrichment.csv": sha256(humanization),
            "ahba_formal_three_tier_sample_detail.csv": sha256(sample_detail),
        },
    }
    return result, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--humanization", type=Path, default=DEFAULT_HUMANIZATION)
    parser.add_argument("--sample-detail", type=Path, default=DEFAULT_SAMPLE_DETAIL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    for path in (args.humanization, args.sample_detail):
        if not path.exists():
            raise FileNotFoundError(
                f"input not found: {path}; provide the validation input with the corresponding CLI option"
            )
    result, manifest = derive(args.humanization, args.sample_detail)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False, float_format="%.10g")
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(result.to_string(index=False))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
