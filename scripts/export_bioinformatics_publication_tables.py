from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "tables_publication"
FIG_DATA = ROOT / "manuscript" / "figures_publication" / "source_data"
AHBA_FORMAL = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "ahba_external_formal_three_tier"
    / "ahba_formal_three_tier_metrics.csv"
)


def wilson(hits: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = hits / n
    denominator = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denominator
    half_width = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return centre - half_width, centre + half_width


def estimate(hits: int, n: int) -> str:
    low, high = wilson(hits, n)
    return f"{100 * hits / n:.1f}% ({100 * low:.1f}-{100 * high:.1f})"


def hits_from_accuracy(metric: float, n: int) -> int:
    return int(round(float(metric) * int(n)))


def markdown_table(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    table1 = pd.DataFrame(
        [
            ["Bo2023 macaque atlas", "Macaque brain tissue RNA-seq", "819 samples; 9 monkeys", "110 regions; 10 Networks", "Model development and internal LOSO/LOMO validation"],
            ["AHBA RNA-seq", "Normal human brain tissue RNA-seq", "242 total; 233 supported", "Formal Network/resolution-group/exact-region endpoints; 91 exact-mapped", "Cross-species external validation"],
            ["Ivy GAP", "Human glioblastoma anatomic-structure RNA-seq", "122 samples", "5 tumor microanatomic structures", "Disease-domain prediction distribution; no location accuracy"],
            ["TCGA GBM/LGG", "Human glioma bulk tissue RNA-seq", "801 samples; 800 patients", "GBM and LGG projects", "Disease-domain prediction distribution"],
            ["TCIA-linked TCGA", "RNA-seq patients with MRI collections", "156 matched; 105 segmentation-ready; 73 complete BraTS4", "MRI location truth pending", "Draft B validation cohort assembly"],
        ],
        columns=["Dataset", "Material", "Evaluable size", "Label space", "Role"],
    )

    endpoint = pd.read_csv(FIG_DATA / "Figure2_endpoint_metrics.csv")
    role_by_endpoint = {
        "Network": "Primary",
        "Region Group": "Secondary",
        "Exact Region": "Exploratory",
    }
    table2_rows: list[list[str | int]] = []
    for row in endpoint.itertuples(index=False):
        n = int(row.n)
        top1_hits = int(getattr(row, "top1_hits", hits_from_accuracy(row.top1, n)))
        top3_hits = int(getattr(row, "top3_hits", hits_from_accuracy(row.top3, n)))
        table2_rows.append(
            [
                row.validation,
                row.endpoint,
                getattr(row, "role", role_by_endpoint.get(row.endpoint, "")),
                f"{top1_hits}/{n}",
                estimate(top1_hits, n),
                f"{top3_hits}/{n}",
                estimate(top3_hits, n),
            ]
        )
    table2 = pd.DataFrame(
        table2_rows,
        columns=["Validation", "Endpoint", "Analysis role", "Top1 hits/n", "Top1 accuracy (95% CI)", "Top3 hits/n", "Top3 accuracy (95% CI)"],
    )

    ahba_metrics = pd.read_csv(AHBA_FORMAL)
    ahba = ahba_metrics.loc[ahba_metrics["route"].eq("hybrid_projected_network_logcpm_exact")].iloc[0]
    mri = json.loads(
        (ROOT / "results" / "tcga_rnaseq_tcia_mri_collection_match_20260605" / "tcga_rnaseq_tcia_mri_match_summary.json").read_text(encoding="utf-8")
    )
    ahba_supported = int(ahba["n_samples_supported_for_accuracy"])
    ahba_exact_n = int(ahba["n_samples_exact_region_evaluable"])
    ahba_network_top1 = hits_from_accuracy(ahba["network_top1_accuracy_coarse"], ahba_supported)
    ahba_network_top3 = hits_from_accuracy(ahba["network_top3_accuracy_coarse"], ahba_supported)
    ahba_group_top1 = hits_from_accuracy(ahba["group_top1_accuracy_exact_mapped"], ahba_exact_n)
    ahba_group_top3 = hits_from_accuracy(ahba["group_top3_accuracy_exact_mapped"], ahba_exact_n)
    ahba_exact_top1 = hits_from_accuracy(ahba["region_top1_accuracy_exact_mapped"], ahba_exact_n)
    ahba_exact_top3 = hits_from_accuracy(ahba["region_top3_accuracy_exact_mapped"], ahba_exact_n)

    table3 = pd.DataFrame(
        [
            ["AHBA", "Network Top1", f"{ahba_network_top1}/{ahba_supported}", estimate(ahba_network_top1, ahba_supported), "Primary endpoint; harmonized cross-species labels"],
            ["AHBA", "Network Top3", f"{ahba_network_top3}/{ahba_supported}", estimate(ahba_network_top3, ahba_supported), "Primary endpoint; harmonized cross-species labels"],
            ["AHBA", "Resolution Group Top1", f"{ahba_group_top1}/{ahba_exact_n}", estimate(ahba_group_top1, ahba_exact_n), "Secondary endpoint; exact-mapped subset"],
            ["AHBA", "Resolution Group Top3", f"{ahba_group_top3}/{ahba_exact_n}", estimate(ahba_group_top3, ahba_exact_n), "Secondary endpoint; exact-mapped subset"],
            ["AHBA", "Exact Region Top1", f"{ahba_exact_top1}/{ahba_exact_n}", estimate(ahba_exact_top1, ahba_exact_n), "Exploratory; exact-mapped subset"],
            ["AHBA", "Exact Region Top3", f"{ahba_exact_top3}/{ahba_exact_n}", estimate(ahba_exact_top3, ahba_exact_n), "Exploratory; exact-mapped subset"],
            ["TCGA-TCIA", "MRI collection match", f"{mri['n_matched_patients']}/{mri['n_rnaseq_patients']}", f"{100 * mri['matched_fraction']:.1f}%", "Cohort coverage only; not tracing accuracy"],
            ["TCGA-TCIA", "Segmentation-ready", "105/800", "13.1%", "Draft B candidate cohort"],
            ["TCGA-TCIA", "Complete BraTS4", "73/800", "9.1%", "Draft B highest-priority cohort"],
        ],
        columns=["Dataset", "Metric", "Count", "Estimate (95% CI where applicable)", "Interpretation"],
    )

    table1.to_csv(OUT / "Table1_datasets_and_roles.csv", index=False, encoding="utf-8-sig")
    table2.to_csv(OUT / "Table2_internal_validation.csv", index=False, encoding="utf-8-sig")
    table3.to_csv(OUT / "Table3_external_validation_and_mri_linkage.csv", index=False, encoding="utf-8-sig")

    text = f"""# Publication-ready tables for Bioinformatics Draft A

## Table 1. Datasets, material and analytical roles

{markdown_table(table1)}

**Legend.** Network is the sole primary endpoint. Region Group is secondary and Exact Region is exploratory. cfRNA is a prospective application; all current model-building and validation data are tissue RNA-seq.

## Table 2. Hierarchical internal validation performance

{markdown_table(table2)}

**Legend.** Values are binomial proportions with two-sided Wilson 95% confidence intervals. LOSO denotes strict leave-one-sample-out validation; LOMO denotes leave-one-monkey-out validation. The main internal results use the P0 formal hard-evidence hybrid route.

## Table 3. Human external validation and MRI-linked cohort coverage

{markdown_table(table3)}

**Legend.** AHBA accuracy uses the P0 formal three-tier hybrid route on harmonized human-to-macaque labels. TCGA-TCIA rows describe cohort availability only. Tumor-location accuracy must not be reported until segmentation-derived or curated MRI location truth is available.

## Reporting rules

- Keep Network as the only primary endpoint in the title, abstract, main text and figures.
- Describe Region Group as secondary and Exact Region as exploratory.
- Do not interpret Ivy GAP or TCGA prediction distributions as location accuracy.
- State that cfRNA is prospective and that the present evidence is derived primarily from tissue RNA-seq.
"""
    (OUT / "publication_tables_and_legends.md").write_text(text, encoding="utf-8")
    print(f"Publication tables written to {OUT}")


if __name__ == "__main__":
    main()
