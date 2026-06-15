from __future__ import annotations

from pathlib import Path

import pandas as pd

from export_bioinformatics_publication_tables import estimate, markdown_table


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "tables_publication_no_TCIA_NBIA"
FIG_DATA = ROOT / "manuscript" / "figures_publication_no_TCIA_NBIA" / "source_data"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    table1 = pd.DataFrame(
        [
            ["Bo2023 macaque atlas", "Macaque brain tissue RNA-seq", "819 samples; 9 monkeys", "110 regions; 10 Networks", "Model development and internal LOSO/LOMO validation"],
            ["AHBA RNA-seq", "Normal human brain tissue RNA-seq", "242 total; 233 supported", "Harmonized Network/broad anatomy; 91 exact-mapped", "Cross-species external validation"],
            ["Ivy GAP", "Human glioblastoma anatomic-structure RNA-seq", "122 samples", "5 tumor microanatomic structures", "Disease-domain prediction distribution; no location accuracy"],
            ["TCGA GBM/LGG", "Human glioma bulk tissue RNA-seq", "801 samples; 800 patients", "GBM and LGG projects", "Disease-domain prediction distribution; no location accuracy"],
        ],
        columns=["Dataset", "Material", "Evaluable size", "Label space", "Role"],
    )

    endpoint = pd.read_csv(FIG_DATA / "Figure2_endpoint_metrics.csv")
    table2_rows = []
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

    table3 = pd.DataFrame(
        [
            ["Network Top1", "58/233", estimate(58, 233), "Primary endpoint; harmonized cross-species labels"],
            ["Network Top3", "129/233", estimate(129, 233), "Primary endpoint; harmonized cross-species labels"],
            ["Broad anatomy Top1", "103/233", estimate(103, 233), "Secondary coarse-anatomy evidence"],
            ["Exact Region Top1", "9/91", estimate(9, 91), "Exploratory; exact-mapped subset"],
            ["Exact Region Top3", "27/91", estimate(27, 91), "Exploratory; exact-mapped subset"],
        ],
        columns=["AHBA metric", "Hits/n", "Accuracy (95% CI)", "Interpretation"],
    )

    table1.to_csv(OUT / "Table1_datasets_and_roles.csv", index=False, encoding="utf-8-sig")
    table2.to_csv(OUT / "Table2_internal_validation.csv", index=False, encoding="utf-8-sig")
    table3.to_csv(OUT / "Table3_AHBA_external_validation.csv", index=False, encoding="utf-8-sig")

    text = f"""# Publication-ready tables for Bioinformatics Draft A

## Table 1. Datasets, material and analytical roles

{markdown_table(table1)}

**Legend.** Network is the sole primary endpoint. Region Group is secondary and Exact Region is exploratory. cfRNA is a prospective application; all current model-building and validation data are tissue RNA-seq.

## Table 2. Hierarchical internal validation performance

{markdown_table(table2)}

**Legend.** Values are binomial proportions with two-sided Wilson 95% confidence intervals. LOSO denotes strict leave-one-sample-out validation; LOMO denotes leave-one-monkey-out validation. The pairwise-rescue margin of 0.002 was selected retrospectively in the LOSO threshold screen and requires independent or nested confirmation.

## Table 3. AHBA human external-validation performance

{markdown_table(table3)}

**Legend.** AHBA accuracy uses harmonized human-to-macaque labels and is not strict macaque exact-region accuracy. Ivy GAP and TCGA are retained as unlabeled disease-domain analyses and are therefore summarized graphically rather than as accuracy estimates.

## Reporting rules

- Keep Network as the only primary endpoint in the title, abstract, main text and figures.
- Describe Region Group as secondary and Exact Region as exploratory.
- Do not interpret Ivy GAP or TCGA prediction distributions as location accuracy.
- State that cfRNA is prospective and that the present evidence is derived primarily from tissue RNA-seq.
- Keep this manuscript version focused on transcriptomic evidence and omit the excluded imaging-database workflow.
"""
    (OUT / "publication_tables_and_legends.md").write_text(text, encoding="utf-8")
    print(f"Publication tables without TCIA/NBIA written to {OUT}")


if __name__ == "__main__":
    main()
