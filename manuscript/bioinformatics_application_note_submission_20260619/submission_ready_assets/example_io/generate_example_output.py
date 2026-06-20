from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
GENES = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"
DB = ROOT / "cfrna_source_tracing.db"

sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import trace_bo2023_secondary_regions
from core.network_tracing import trace_network_expression


def main() -> None:
    genes = pd.read_csv(GENES)
    expression = genes[["gene_symbol"]].copy()
    expression["tpm_value"] = [round(1.0 + ((i * 37) % 200) / 10.0, 3) for i in range(len(expression))]
    expression.to_csv(OUT / "example_expression_input_full_200genes.csv", index=False)

    network = trace_network_expression(expression)
    (OUT / "example_network_output.json").write_text(
        json.dumps(network, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    pd.DataFrame(network.get("results", [])).to_csv(
        OUT / "example_network_ranked_candidates.csv",
        index=False,
    )

    region = trace_bo2023_secondary_regions(expression, network, db_path=str(DB), atlas_id=1)
    (OUT / "example_three_tier_output.json").write_text(
        json.dumps(region, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    exact = pd.DataFrame(region.get("results", []))
    exact.to_csv(OUT / "example_exact_region_ranked_candidates.csv", index=False)

    if not exact.empty:
        cols = [
            "resolution_group_rank",
            "resolution_group",
            "resolution_group_score",
            "network_id",
            "resolution_group_members",
            "group_plausibility_tier",
            "group_calibration_flags",
            "manual_review_recommended",
        ]
        groups = (
            exact.sort_values(["resolution_group_rank", "rank"])
            .drop_duplicates("resolution_group_rank")[cols]
            .sort_values("resolution_group_rank")
        )
    else:
        groups = pd.DataFrame()
    groups.to_csv(OUT / "example_resolution_group_ranked_candidates.csv", index=False)

    print("network_traceability", network["meta"].get("traceability"))
    print("network_overlap", network["meta"].get("n_overlap_genes"), network["meta"].get("overlap_fraction"))
    print("three_tier_traceability", region["meta"].get("traceability"))
    print("exact_candidates", len(region.get("results", [])))


if __name__ == "__main__":
    main()
