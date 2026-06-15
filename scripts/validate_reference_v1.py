#!/usr/bin/env python
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DB = ROOT / "data" / "processed_reference" / "v1" / "cfrna_brain_tracing_reference_v1.sqlite"
OUT = ROOT / "data" / "processed_reference" / "v1" / "validation"

CORE = {
    "neuron": ["RBFOX3", "SNAP25", "SYT1", "TUBB3", "MAP2"],
    "astrocyte": ["GFAP", "AQP4", "ALDH1L1", "SLC1A3"],
    "oligodendrocyte": ["MBP", "PLP1", "MOBP", "MAG", "MOG"],
    "OPC": ["PDGFRA", "VCAN", "CSPG4"],
    "microglia": ["P2RY12", "CX3CR1", "TMEM119", "AIF1"],
    "endothelial_vascular": ["CLDN5", "PECAM1", "VWF", "KDR"],
    "RBC_contamination": ["HBA1", "HBA2", "HBB", "ALAS2", "CA1", "CA2", "SLC4A1", "GYPA"],
    "immune_contamination": ["PTPRC", "LST1", "TYROBP", "CD74", "HLA-DRA", "S100A8", "S100A9"],
}

TYPE_MAP = {"brain_region": "brain_region_marker", "brain_celltype": "brain_celltype_marker", "tbi_injury": "tbi_injury_response"}
EXPECTED = ["brain_enriched", "brain_region_marker", "brain_celltype_marker", "tbi_injury_response", "contamination_marker"]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB)
    markers = pd.read_sql("select * from marker_candidates", con)
    contam = pd.read_sql("select * from contamination_markers", con)
    universe = pd.read_sql("select gene_symbol from gene_id_map", con)
    con.close()

    markers["marker_type_normalized"] = markers["marker_type"].replace(TYPE_MAP)
    cat_rows = []
    for cat in EXPECTED:
        if cat == "contamination_marker":
            count = len(contam)
        else:
            count = int((markers["marker_type_normalized"] == cat).sum())
        cat_rows.append({"marker_category": cat, "row_count": count, "exists": count > 0})
    cat_df = pd.DataFrame(cat_rows)
    cat_df.to_csv(OUT / "reference_v1_marker_category_summary.tsv", sep="\t", index=False)

    available = set(universe["gene_symbol"].astype(str).str.upper()) | set(markers["gene_symbol"].astype(str).str.upper()) | set(contam["gene_symbol"].astype(str).str.upper())
    rows = []
    for group, genes in CORE.items():
        for gene in genes:
            rows.append({"marker_group": group, "gene_symbol": gene, "present": gene in available})
    core_df = pd.DataFrame(rows)
    core_df.to_csv(OUT / "reference_v1_core_marker_check.tsv", sep="\t", index=False)

    missing_core = core_df[(~core_df["present"]) & (~core_df["marker_group"].str.contains("contamination"))]
    missing_contam = core_df[(~core_df["present"]) & (core_df["marker_group"].str.contains("contamination"))]
    brain_n = int(cat_df.loc[cat_df["marker_category"] == "brain_enriched", "row_count"].iloc[0])
    injury_n = int(cat_df.loc[cat_df["marker_category"] == "tbi_injury_response", "row_count"].iloc[0])
    lines = [
        "# Reference v1 Validation Report",
        "",
        f"- SQLite database: `{DB}`",
        f"- brain_enriched marker count: {brain_n}",
        f"- tbi_injury_response marker count: {injury_n}",
        "",
        "## Marker Categories",
        "",
        cat_df.to_csv(sep="\t", index=False),
        "",
        "## Core Marker Missingness",
        "",
        f"- Missing core brain markers: {', '.join(missing_core['gene_symbol']) if not missing_core.empty else 'None'}",
        f"- Missing contamination markers: {', '.join(missing_contam['gene_symbol']) if not missing_contam.empty else 'None'}",
        "",
        "## Interpretation Notes",
        "",
        "- brain_enriched marker count is considered potentially low if below 50; current count is " + str(brain_n) + ".",
        "- injury-state markers depend on Garza TBI count matrices; Garza metadata is limited and current TBI/control grouping was inferred from sample names.",
        "- HPA aggregate lacks detection fraction, so celltype markers are candidate markers rather than final gold-standard markers.",
        "- Outputs should be interpreted as evidence, signal, risk, and candidate source. They are not clinical decision or definitive origin calls.",
    ]
    (OUT / "reference_v1_validation_report.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote validation outputs to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
