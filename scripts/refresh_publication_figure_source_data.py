from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATA = ROOT / "manuscript" / "figures_publication" / "source_data"
ARCHIVED_SOURCE_DATA = ROOT / "manuscript" / "figures_publication_20260613" / "source_data"
LOSO_LOMO = ROOT / "results" / "bo2023_loso_lomo_summary_figures_20260601"
P0_METRICS = ROOT / "reports" / "p0_hard_evidence_20260629" / "metric_summary_with_ci.csv"
AHBA = ROOT / "results" / "ahba_human_rnaseq_external_validation_20260603"
AHBA_FORMAL = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "ahba_external_formal_three_tier"
)
IVY = ROOT / "results" / "ivy_gap_anatomic_rnaseq_tracing_20260603"


def copy_loso_lomo_sources() -> None:
    copies = {
        "loso_lomo_endpoint_summary.csv": "Figure2_endpoint_metrics_route_development_legacy.csv",
        "lomo_per_monkey_summary.csv": "Figure3_lomo_per_monkey.csv",
    }
    for source_name, target_name in copies.items():
        source = LOSO_LOMO / source_name
        if source.exists():
            shutil.copyfile(source, SOURCE_DATA / target_name)

    margin_screen = ARCHIVED_SOURCE_DATA / "Figure3_margin_screen.csv"
    if margin_screen.exists():
        shutil.copyfile(margin_screen, SOURCE_DATA / "Figure3_margin_screen.csv")

    per_monkey = pd.read_csv(LOSO_LOMO / "lomo_per_monkey_summary.csv")
    figure_s1 = per_monkey.rename(
        columns={
            "network_top1": "Network Top1",
            "network_top3": "Network Top3",
            "group_top1": "Region Group Top1",
            "group_top3": "Region Group Top3",
            "exact_top1": "Exact Region Top1",
            "exact_top3": "Exact Region Top3",
        }
    )[
        [
            "monkey_id",
            "n",
            "Network Top1",
            "Network Top3",
            "Region Group Top1",
            "Region Group Top3",
            "Exact Region Top1",
            "Exact Region Top3",
        ]
    ]
    figure_s1.to_csv(SOURCE_DATA / "FigureS1_lomo_hierarchical_per_monkey.csv", index=False)


def refresh_figure2_formal_sources() -> None:
    metrics = pd.read_csv(P0_METRICS)
    endpoint_order = [
        ("formal_internal_network_loso", "LOSO", "Network", "Primary"),
        ("formal_internal_resolution_group_loso", "LOSO", "Region Group", "Secondary"),
        ("formal_internal_exact_region_loso", "LOSO", "Exact Region", "Exploratory"),
        ("formal_internal_network_lomo", "LOMO", "Network", "Primary"),
        ("formal_internal_resolution_group_lomo", "LOMO", "Region Group", "Secondary"),
        ("formal_internal_exact_region_lomo", "LOMO", "Exact Region", "Exploratory"),
    ]
    rows: list[dict[str, object]] = []
    for endpoint_key, validation, endpoint, role in endpoint_order:
        endpoint_metrics = metrics.loc[metrics["endpoint"].eq(endpoint_key)]
        top1 = endpoint_metrics.loc[endpoint_metrics["metric"].eq("Top1")].iloc[0]
        top3 = endpoint_metrics.loc[endpoint_metrics["metric"].eq("Top3")].iloc[0]
        rows.append(
            {
                "validation": validation,
                "endpoint": endpoint,
                "n": int(top1["n"]),
                "top1": float(top1["accuracy"]),
                "top3": float(top3["accuracy"]),
                "top1_hits": int(top1["hits"]),
                "top3_hits": int(top3["hits"]),
                "role": role,
                "notes": "P0 formal hard-evidence hybrid_projected_network_logcpm_exact",
            }
        )
    pd.DataFrame(rows).to_csv(SOURCE_DATA / "Figure2_endpoint_metrics.csv", index=False)


def refresh_ahba_sources() -> None:
    legacy = json.loads((AHBA / "ahba_rnaseq_external_validation_metrics.json").read_text(encoding="utf-8"))
    legacy_supported_n = int(legacy["n_samples_supported_for_accuracy"])
    legacy_exact_n = int(legacy["exact_region_evaluable_samples"])
    pd.DataFrame(
        [
            ["Network Top1", legacy_supported_n, legacy["network_top1_accuracy_coarse"]],
            ["Network Top3", legacy_supported_n, legacy["network_top3_accuracy_coarse"]],
            ["Broad anatomy Top1", legacy_supported_n, legacy["region_lobe_top1_accuracy_coarse"]],
            ["Exact Region Top1", legacy_exact_n, legacy["region_top1_exact_accuracy_on_exact_mapped_labels"]],
            ["Exact Region Top3", legacy_exact_n, legacy["region_top3_exact_accuracy_on_exact_mapped_labels"]],
        ],
        columns=["metric", "n", "accuracy"],
    ).to_csv(SOURCE_DATA / "Figure4_ahba_metrics_standard_legacy.csv", index=False)

    metrics = pd.read_csv(AHBA_FORMAL / "ahba_formal_three_tier_metrics.csv")
    hybrid = metrics.loc[metrics["route"].eq("hybrid_projected_network_logcpm_exact")].iloc[0]
    supported_n = int(hybrid["n_samples_supported_for_accuracy"])
    exact_n = int(hybrid["n_samples_exact_region_evaluable"])
    pd.DataFrame(
        [
            ["Network Top1", supported_n, hybrid["network_top1_accuracy_coarse"]],
            ["Network Top3", supported_n, hybrid["network_top3_accuracy_coarse"]],
            ["Resolution Group Top1", exact_n, hybrid["group_top1_accuracy_exact_mapped"]],
            ["Resolution Group Top3", exact_n, hybrid["group_top3_accuracy_exact_mapped"]],
            ["Exact Region Top1", exact_n, hybrid["region_top1_accuracy_exact_mapped"]],
            ["Exact Region Top3", exact_n, hybrid["region_top3_accuracy_exact_mapped"]],
        ],
        columns=["metric", "n", "accuracy"],
    ).to_csv(SOURCE_DATA / "Figure4_ahba_metrics.csv", index=False)

    sample = pd.read_csv(AHBA / "ahba_rnaseq_external_validation_sample_summary.csv")
    sample = sample[sample["public_major_anatomy"].notna() & sample["predicted_public_major_top1"].notna()]
    confusion = pd.crosstab(
        sample["public_major_anatomy"],
        sample["predicted_public_major_top1"],
        normalize="index",
    )
    confusion.index.name = "public_major_anatomy"
    confusion.reset_index().to_csv(SOURCE_DATA / "Figure4_ahba_broad_confusion_standard_legacy.csv", index=False)

    formal = pd.read_csv(AHBA_FORMAL / "ahba_formal_three_tier_sample_detail.csv")
    formal = formal.loc[
        formal["route"].eq("hybrid_projected_network_logcpm_exact")
        & formal["supported_for_accuracy"].astype(bool)
        & formal["allowed_bo2023_networks"].notna()
        & formal["network_top1"].notna()
    ].copy()
    formal["allowed_network_set"] = formal["allowed_bo2023_networks"].astype(str)
    formal_confusion = pd.crosstab(
        formal["allowed_network_set"],
        formal["network_top1"],
        normalize="index",
    )
    formal_confusion.index.name = "allowed_bo2023_networks"
    formal_confusion.reset_index().to_csv(SOURCE_DATA / "Figure4_ahba_broad_confusion.csv", index=False)


def refresh_ivy_sources() -> None:
    distribution = pd.read_csv(IVY / "ivy_gap_anatomic_structure_prediction_distributions.csv")
    network = distribution.loc[distribution["endpoint"].eq("network_top1")]
    pivot = network.pivot_table(
        index="structure_acronym",
        columns="value",
        values="fraction",
        aggfunc="sum",
        fill_value=0.0,
    )
    pivot.reset_index().to_csv(SOURCE_DATA / "Figure5_ivy_network_distribution.csv", index=False)


def main() -> None:
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    copy_loso_lomo_sources()
    refresh_figure2_formal_sources()
    refresh_ahba_sources()
    refresh_ivy_sources()
    print(f"Figure source data refreshed in {SOURCE_DATA}")


if __name__ == "__main__":
    main()
