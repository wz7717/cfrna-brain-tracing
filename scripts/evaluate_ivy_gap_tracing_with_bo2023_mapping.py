#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TRACING_DIR = ROOT / "results" / "ivy_gap_anatomic_rnaseq_tracing_20260603"
DEFAULT_MAPPING = ROOT / "docs" / "deliverables" / "bo2023_public_anatomical_label_mapping.csv"
DEFAULT_OUTDIR = ROOT / "results" / "ivy_gap_anatomic_mapping_interpretation_20260603"


def mode_fraction(values: pd.Series) -> tuple[str, int, float]:
    values = values.fillna("").astype(str)
    if values.empty:
        return "", 0, 0.0
    counts = values.value_counts()
    value = str(counts.index[0])
    count = int(counts.iloc[0])
    return value, count, float(count / len(values))


def main() -> int:
    parser = argparse.ArgumentParser(description="Interpret Ivy GAP tracing results with Bo2023 public anatomical mapping.")
    parser.add_argument("--tracing-dir", type=Path, default=DEFAULT_TRACING_DIR)
    parser.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    sample = pd.read_csv(args.tracing_dir / "ivy_gap_anatomic_tracing_sample_summary.csv")
    mapping = pd.read_csv(args.mapping)
    sample["region_top1"] = sample["region_top1"].astype(str)
    mapping["bo2023_region_id"] = mapping["bo2023_region_id"].astype(str)

    cols = [
        "bo2023_region_id",
        "bo2023_region_name",
        "bo2023_lobe",
        "bo2023_saleem_network",
        "resolution_tier",
        "resolution_group",
        "public_major_anatomy",
        "public_network_or_subdivision",
        "public_fine_anatomy_label",
        "mapping_level",
        "mapping_confidence",
        "ivy_gap_mapping_status",
        "do_not_use_as_exact_ground_truth",
        "caution_note",
    ]
    mapped = sample.merge(
        mapping[cols],
        left_on="region_top1",
        right_on="bo2023_region_id",
        how="left",
    )
    mapped["mapping_found_for_region_top1"] = mapped["bo2023_region_id"].notna()
    mapped.to_csv(args.outdir / "ivy_gap_sample_tracing_with_public_mapping.csv", index=False, encoding="utf-8-sig")

    rows = []
    for structure, sub in mapped.groupby("structure_acronym"):
        net, net_n, net_frac = mode_fraction(sub["network_top1"])
        region, region_n, region_frac = mode_fraction(sub["region_top1"])
        public_major, public_major_n, public_major_frac = mode_fraction(sub["public_major_anatomy"])
        public_fine, public_fine_n, public_fine_frac = mode_fraction(sub["public_fine_anatomy_label"])
        group, group_n, group_frac = mode_fraction(sub["resolution_group_top1"])
        low_frac = float((sub["top1_resolution_tier"].astype(str) == "low_resolution").mean())
        high_mapping_frac = float((sub["mapping_confidence"].astype(str) == "high").mean())
        rows.append(
            {
                "ivy_structure": structure,
                "n_samples": int(len(sub)),
                "dominant_network_top1": net,
                "dominant_network_top1_n": net_n,
                "dominant_network_top1_fraction": net_frac,
                "dominant_region_top1": region,
                "dominant_region_top1_n": region_n,
                "dominant_region_top1_fraction": region_frac,
                "dominant_public_major_anatomy": public_major,
                "dominant_public_major_anatomy_n": public_major_n,
                "dominant_public_major_anatomy_fraction": public_major_frac,
                "dominant_public_fine_anatomy": public_fine,
                "dominant_public_fine_anatomy_n": public_fine_n,
                "dominant_public_fine_anatomy_fraction": public_fine_frac,
                "dominant_resolution_group": group,
                "dominant_resolution_group_n": group_n,
                "dominant_resolution_group_fraction": group_frac,
                "top1_low_resolution_fraction": low_frac,
                "top1_high_confidence_mapping_fraction": high_mapping_frac,
                "accuracy_status": "not_computable_no_direct_ivy_to_normal_region_ground_truth",
            }
        )
    structure = pd.DataFrame(rows).sort_values("ivy_structure")
    structure.to_csv(args.outdir / "ivy_gap_structure_mapping_interpretation_summary.csv", index=False, encoding="utf-8-sig")

    all_metrics = {
        "n_samples": int(len(mapped)),
        "n_ivy_structures": int(mapped["structure_acronym"].nunique()),
        "region_top1_mapping_coverage": float(mapped["mapping_found_for_region_top1"].mean()),
        "low_resolution_top1_fraction": float((mapped["top1_resolution_tier"].astype(str) == "low_resolution").mean()),
        "manual_review_top1_fraction": float(mapped["top1_manual_review_recommended"].astype(bool).mean()),
        "high_confidence_public_mapping_fraction": float((mapped["mapping_confidence"].astype(str) == "high").mean()),
        "true_accuracy_status": "not_computable",
        "reason": (
            "Ivy GAP labels CT/CTmvp/CTpan/IT/LE are GBM microanatomic structures. "
            "The Bo2023 mapping table explicitly marks Ivy GAP as no direct normal brain region mapping, "
            "so exact or coarse brain-region accuracy cannot be computed without an additional externally justified label map."
        ),
        "recommended_report_term": "mapped prediction distribution / external tumor transcriptome stress test",
    }
    (args.outdir / "ivy_gap_mapping_interpretation_metrics.json").write_text(
        json.dumps(all_metrics, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    lines = [
        "# Ivy GAP tracing interpreted with Bo2023 public anatomical mapping",
        "",
        "## Accuracy status",
        "",
        "Formal accuracy is not computable for this dataset because Ivy GAP structure labels are tumor microanatomic structures, not normal Bo2023 brain region labels.",
        "",
        "## Computable proxy summaries",
        "",
        f"- Samples: {all_metrics['n_samples']}",
        f"- Region Top1 mapping coverage: {all_metrics['region_top1_mapping_coverage']:.1%}",
        f"- Top1 low-resolution/manual-review fraction: {all_metrics['low_resolution_top1_fraction']:.1%}",
        f"- Top1 high-confidence public anatomical mapping fraction: {all_metrics['high_confidence_public_mapping_fraction']:.1%}",
        "",
        "## By Ivy structure",
        "",
        structure.to_string(index=False),
        "",
        "Interpretation: use these values as distribution concentration by Ivy structure, not as accuracy.",
    ]
    (args.outdir / "ivy_gap_mapping_interpretation_report.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps(all_metrics, ensure_ascii=False, indent=2))
    print(structure.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
