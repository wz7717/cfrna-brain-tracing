from __future__ import annotations

from pathlib import Path

import pandas as pd

from export_bioinformatics_publication_tables import estimate, markdown_table


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "tables_publication_20260613"
FIG_DATA = ROOT / "manuscript" / "figures_publication_20260613" / "source_data"


def save_table(df: pd.DataFrame, stem: str) -> None:
    df.to_csv(OUT / f"{stem}.csv", index=False, encoding="utf-8-sig")
    (OUT / f"{stem}.md").write_text(markdown_table(df), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    table1 = pd.DataFrame(
        [
            ["Bo2023 macaque atlas", "Macaque brain tissue RNA-seq", "819 samples; 9 monkeys", "110 regions; 10 Networks", "Model development; LOSO and donor-isolated validation"],
            ["AHBA", "Normal human brain RNA-seq", "242 total; 233 supported", "Harmonized Network/broad anatomy; 91 exact-mapped", "Cross-species external validation"],
            ["TCGA GBM/LGG", "Human glioma bulk RNA-seq", "800 patients", "GBM and LGG projects", "Unpaired disease-domain shift analysis"],
            ["TCGA-LGG + BraTS", "Paired tumor RNA-seq and MRI", "65 patients; Network n=64", "Direct tumor-mask overlap with TZO116 labels", "Paired anatomical validation"],
            ["Ivy GAP", "Glioblastoma microanatomic RNA-seq", "122 samples", "Five tumor structures", "Exploratory disease-domain characterization"],
            ["GSE228512", "Serum extracellular-vesicle RNA", "85 GBM; 31 healthy", "No patient-level anatomical truth", "Liquid-biopsy transfer stress test"],
            ["GSE106804", "Tumor-associated EV capture RNA", "13 GBM; 6 healthy", "No patient-level anatomical truth", "Enrichment and transfer stress test"],
            ["GSE189919", "Cerebrospinal-fluid RNA", "40 medulloblastoma; 11 normal", "No patient-level anatomical truth", "Technical validation and algorithm audit"],
        ],
        columns=["Dataset", "Material", "Evaluable size", "Label space", "Role"],
    )

    endpoint = pd.read_csv(FIG_DATA / "Figure2_endpoint_metrics.csv")
    table2 = pd.DataFrame(
        [
            [
                row.validation, row.endpoint, row.role,
                f"{row.top1_hits}/{row.n}", estimate(int(row.top1_hits), int(row.n)),
                f"{row.top3_hits}/{row.n}", estimate(int(row.top3_hits), int(row.n)),
            ]
            for row in endpoint.itertuples(index=False)
        ],
        columns=["Validation", "Endpoint", "Analysis role", "Top-1 hits/n",
                 "Top-1 accuracy (95% CI)", "Top-3 hits/n",
                 "Top-3 accuracy (95% CI)"],
    )

    table3 = pd.DataFrame(
        [
            ["AHBA", "Network Top-1", "76/233", estimate(76, 233), "Latest 2026-06-13 summary; primary endpoint"],
            ["AHBA", "Network Top-3", "129/233", estimate(129, 233), "Harmonized cross-species labels"],
            ["AHBA", "Broad anatomy Top-1", "103/233", estimate(103, 233), "Supportive coarse-anatomy endpoint"],
            ["AHBA", "Exact Region Top-1", "9/91", estimate(9, 91), "Exploratory exact-mapped subset"],
            ["AHBA", "Exact Region Top-3", "27/91", estimate(27, 91), "Exploratory exact-mapped subset"],
            ["TCGA-LGG/BraTS", "Lobe Top-3 strict", "55/65", estimate(55, 65), "Supportive paired validation"],
            ["TCGA-LGG/BraTS", "Lobe Top-3 tolerant", "58/65", estimate(58, 65), "Supportive paired validation"],
            ["TCGA-LGG/BraTS", "Broad anatomy Top-3 strict", "49/65", estimate(49, 65), "Principal disease-scenario readout"],
            ["TCGA-LGG/BraTS", "Broad anatomy Top-3 tolerant", "54/65", estimate(54, 65), "Principal disease-scenario readout"],
            ["TCGA-LGG/BraTS", "Network Top-3 strict", "14/64", estimate(14, 64), "Secondary approximation"],
            ["TCGA-LGG/BraTS", "Network Top-3 tolerant", "23/64", estimate(23, 64), "Secondary approximation"],
        ],
        columns=["Cohort", "Metric", "Hits/n", "Accuracy (95% CI)", "Interpretation"],
    )

    table_s1 = pd.DataFrame(
        [
            ["Baseline", "Locked", "No", "No", "4.69%", "35.94%", "83.08%"],
            ["Adapted raw TPM", "Exploratory sensitivity", "No", "No", "7.81%", "43.75%", "76.92%"],
            ["Adapted log1p", "Exploratory sensitivity", "No", "No", "14.06%", "43.75%", "64.62%"],
            ["Adapted harmonized", "Transductive sensitivity; tie-aware retraining required", "Yes", "No", "20.31%", "48.44%", "60.00%"],
            ["Harmonized + calibrated", "Cohort-internal experiment; tie-aware retraining required", "Yes", "Yes", "20.31%", "46.88%", "70.77%"],
        ],
        columns=["Route", "Analysis role", "Uses target unlabeled distribution",
                 "Uses paired labels internally", "Network Top-1 tolerant",
                 "Network Top-3 tolerant", "Broad Top-3 tolerant"],
    )

    table_s2 = pd.DataFrame(
        [
            ["GSE228512", "Serum EV RNA", "85 GBM; 31 healthy", "65.9% GBM; 71.0% healthy", "0/116", "High entropy and near-uniform Top-1 probability; no localization accuracy"],
            ["GSE106804", "Tumor-associated EV RNA", "13 GBM; 6 healthy", "23.1% in GBM", "0/19", "Lower class dominance but lower margin and 21.1% adapted-route agreement"],
            ["GSE189919", "CSF RNA", "40 medulloblastoma; 11 normal", "65.0% MB; 45.5% normal", "0/51", "178/200 markers; 67/67 audit checks; Top1 distribution P=0.178; baseline TPM-CPM agreement 87.5%/90.9%"],
        ],
        columns=["Cohort", "Material", "Groups", "Baseline dominant Top-1 fraction",
                 "Separate adapted-model OOD accepted", "Interpretation"],
    )

    save_table(table1, "Table1_datasets_and_roles")
    save_table(table2, "Table2_internal_validation")
    save_table(table3, "Table3_external_validation")
    save_table(table_s1, "TableS1_domain_adaptation_sensitivity")
    save_table(table_s2, "TableS2_liquid_biopsy_stress_tests")
    print(f"Publication tables written to {OUT}")


if __name__ == "__main__":
    main()
