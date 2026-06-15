from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "tables_publication"
FIG_DATA = ROOT / "manuscript" / "figures_publication" / "source_data"


def wilson(hits: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = hits / n
    denominator = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denominator
    half_width = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return centre - half_width, centre + half_width


def estimate(hits: int, n: int) -> str:
    low, high = wilson(hits, n)
    return f"{100 * hits / n:.1f}% ({100 * low:.1f}-{100 * high:.1f})"


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
            ["AHBA RNA-seq", "Normal human brain tissue RNA-seq", "242 total; 233 supported", "Harmonized Network/broad anatomy; 91 exact-mapped", "Cross-species external validation"],
            ["Ivy GAP", "Human glioblastoma anatomic-structure RNA-seq", "122 samples", "5 tumor microanatomic structures", "Disease-domain prediction distribution; no location accuracy"],
            ["TCGA GBM/LGG", "Human glioma bulk tissue RNA-seq", "801 samples; 800 patients", "GBM and LGG projects", "Disease-domain prediction distribution"],
            ["TCIA-linked TCGA", "RNA-seq patients with MRI collections", "156 matched; 105 segmentation-ready; 73 complete BraTS4", "MRI location truth pending", "Draft B validation cohort assembly"],
        ],
        columns=["Dataset", "Material", "Evaluable size", "Label space", "Role"],
    )

    endpoint = pd.read_csv(FIG_DATA / "Figure2_endpoint_metrics.csv")
    table2_rows: list[list[str | int]] = []
    for row in endpoint.itertuples(index=False):
        table2_rows.append(
            [
                row.validation,
                row.endpoint,
                row.role,
                f"{row.top1_hits}/{row.n}",
                estimate(int(row.top1_hits), int(row.n)),
                f"{row.top3_hits}/{row.n}",
                estimate(int(row.top3_hits), int(row.n)),
            ]
        )
    table2 = pd.DataFrame(
        table2_rows,
        columns=["Validation", "Endpoint", "Analysis role", "Top1 hits/n", "Top1 accuracy (95% CI)", "Top3 hits/n", "Top3 accuracy (95% CI)"],
    )

    ahba = json.loads(
        (ROOT / "results" / "ahba_human_rnaseq_external_validation_margin0p002_20260604" / "ahba_rnaseq_external_validation_metrics.json").read_text(encoding="utf-8")
    )
    mri = json.loads(
        (ROOT / "results" / "tcga_rnaseq_tcia_mri_collection_match_20260605" / "tcga_rnaseq_tcia_mri_match_summary.json").read_text(encoding="utf-8")
    )
    table3 = pd.DataFrame(
        [
            ["AHBA", "Network Top1", "58/233", estimate(58, 233), "Primary endpoint; harmonized cross-species labels"],
            ["AHBA", "Network Top3", "129/233", estimate(129, 233), "Primary endpoint; harmonized cross-species labels"],
            ["AHBA", "Broad anatomy Top1", "103/233", estimate(103, 233), "Secondary coarse-anatomy evidence"],
            ["AHBA", "Exact Region Top1", "9/91", estimate(9, 91), "Exploratory; exact-mapped subset"],
            ["AHBA", "Exact Region Top3", "27/91", estimate(27, 91), "Exploratory; exact-mapped subset"],
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

**Legend.** Values are binomial proportions with two-sided Wilson 95% confidence intervals. LOSO denotes strict leave-one-sample-out validation; LOMO denotes leave-one-monkey-out validation. The pairwise-rescue margin of 0.002 was selected retrospectively in the LOSO threshold screen and requires independent or nested confirmation.

## Table 3. Human external validation and MRI-linked cohort coverage

{markdown_table(table3)}

**Legend.** AHBA accuracy uses harmonized human-to-macaque labels and is not strict macaque exact-region accuracy. TCGA-TCIA rows describe cohort availability only. Tumor-location accuracy must not be reported until segmentation-derived or curated MRI location truth is available.

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
