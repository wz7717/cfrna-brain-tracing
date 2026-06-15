#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402
from core.region_resolution import annotate_region_candidates  # noqa: E402


DEFAULT_DATA_DIR = ROOT / "data" / "ivy_gap_anatomic_rnaseq"
DEFAULT_MATRIX = DEFAULT_DATA_DIR / "ivy_gap_anatomic_structure_tpm_gene_symbol_matrix.tsv"
DEFAULT_METADATA = DEFAULT_DATA_DIR / "ivy_gap_anatomic_structure_sample_metadata_122.csv"
DEFAULT_OUTDIR = ROOT / "results" / "ivy_gap_anatomic_rnaseq_tracing_20260603"
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"


def expression_frame(matrix: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": matrix.index.astype(str),
            "tpm_value": pd.to_numeric(matrix[sample_id], errors="coerce").fillna(0.0).to_numpy(),
        }
    )


def top_values(rows: list[dict[str, Any]], key: str, k: int) -> list[str]:
    return [str(row.get(key, "")) for row in sorted(rows, key=lambda r: int(r.get("rank", 999)))[:k]]


def add_distribution_rows(
    records: list[dict[str, Any]],
    df: pd.DataFrame,
    group_col: str,
    value_col: str,
    out_label: str,
) -> None:
    for group_value, sub in df.groupby(group_col, dropna=False):
        counts = sub[value_col].fillna("").astype(str).value_counts(dropna=False)
        denom = int(len(sub))
        for value, count in counts.items():
            records.append(
                {
                    group_col: group_value,
                    "endpoint": out_label,
                    "value": value,
                    "count": int(count),
                    "n_samples": denom,
                    "fraction": float(count / denom) if denom else 0.0,
                }
            )


def add_top3_membership_rows(
    records: list[dict[str, Any]],
    df: pd.DataFrame,
    group_col: str,
    list_col: str,
    out_label: str,
) -> None:
    for group_value, sub in df.groupby(group_col, dropna=False):
        counter: Counter[str] = Counter()
        for values in sub[list_col].fillna("").astype(str):
            for value in values.split(" | "):
                value = value.strip()
                if value:
                    counter[value] += 1
        denom = int(len(sub))
        for value, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
            records.append(
                {
                    group_col: group_value,
                    "endpoint": out_label,
                    "value": value,
                    "count": int(count),
                    "n_samples": denom,
                    "fraction": float(count / denom) if denom else 0.0,
                }
            )


def plot_distribution(df: pd.DataFrame, value_col: str, title: str, out_path: Path, top_n: int = 12) -> None:
    import matplotlib.pyplot as plt

    if df.empty:
        return
    grouped = (
        df.groupby(["structure_acronym", value_col])
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )
    keep_values = grouped.groupby(value_col)["count"].sum().sort_values(ascending=False).head(top_n).index
    grouped = grouped[grouped[value_col].isin(keep_values)]
    pivot = grouped.pivot_table(index="structure_acronym", columns=value_col, values="count", fill_value=0)
    pivot = pivot.div(pivot.sum(axis=1).replace(0, 1), axis=0)

    fig, ax = plt.subplots(figsize=(12, 6))
    pivot.plot(kind="bar", stacked=True, ax=ax, width=0.85)
    ax.set_title(title)
    ax.set_xlabel("Ivy GAP structure")
    ax.set_ylabel("Fraction of samples")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8)
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Trace Ivy GAP Anatomic Structures RNA-seq samples and summarize predictions by true Ivy structure."
    )
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--atlas-id", type=int, default=4)
    parser.add_argument("--topk-region", type=int, default=15)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    metadata = pd.read_csv(args.metadata)
    metadata["sample_id"] = metadata["sample_id"].astype(str)

    missing = sorted(set(metadata["sample_id"]) - set(matrix.columns.astype(str)))
    extra = sorted(set(matrix.columns.astype(str)) - set(metadata["sample_id"]))
    if missing or extra:
        raise ValueError(f"metadata/matrix sample mismatch: missing={missing[:5]}, extra={extra[:5]}")

    meta_by_sample = metadata.set_index("sample_id").to_dict(orient="index")
    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    sample_summary_rows: list[dict[str, Any]] = []
    trace_meta: dict[str, Any] = {}

    for i, sample_id in enumerate(metadata["sample_id"].tolist(), start=1):
        expr = expression_frame(matrix, sample_id)
        sample_meta = meta_by_sample[sample_id]
        network_out = trace_network_expression(expr, min_overlap_fraction=0.20)
        region_out = trace_bo2023_secondary_regions(
            expr,
            network_out,
            str(args.db),
            int(args.atlas_id),
            topk=int(args.topk_region),
        )
        region_out = annotate_region_candidates(region_out, network_out)
        trace_meta[sample_id] = {
            "network": network_out.get("meta", {}),
            "region": region_out.get("meta", {}),
        }

        common = {
            "sample_id": sample_id,
            "tumor_name": sample_meta.get("tumor_name", ""),
            "structure_acronym": sample_meta.get("structure_acronym", ""),
            "structure_name": sample_meta.get("structure_name", ""),
            "study_name": sample_meta.get("study_name", ""),
            "molecular_subtype": sample_meta.get("molecular_subtype", ""),
        }
        for row in network_out.get("results", [])[:10]:
            network_rows.append({**common, **row})
        for row in region_out.get("results", [])[: int(args.topk_region)]:
            region_rows.append({**common, **row})

        network_top3 = top_values(network_out.get("results", []), "network_id", 3)
        region_top3 = top_values(region_out.get("results", []), "region_id", 3)
        group_top3 = top_values(region_out.get("results", []), "resolution_group", 3)
        top_region = (region_out.get("results", []) or [{}])[0]
        sample_summary_rows.append(
            {
                **common,
                "network_top1": network_top3[0] if network_top3 else "",
                "network_top3": " | ".join(network_top3),
                "region_top1": region_top3[0] if region_top3 else "",
                "region_top3": " | ".join(region_top3),
                "resolution_group_top1": group_top3[0] if group_top3 else "",
                "resolution_group_top3": " | ".join(group_top3),
                "top1_resolution_tier": str(top_region.get("resolution_tier", "")),
                "top1_manual_review_recommended": bool(top_region.get("manual_review_recommended", False)),
                "top1_group_plausibility_tier": str(top_region.get("group_plausibility_tier", "")),
                "network_overlap_genes": network_out.get("meta", {}).get("n_overlap_genes"),
                "network_overlap_fraction": network_out.get("meta", {}).get("overlap_fraction"),
                "region_overlap_genes": region_out.get("meta", {}).get("n_overlap_genes"),
            }
        )
        print(f"[{i:03d}/{len(metadata)}] {sample_id} {common['structure_acronym']}")

    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    sample_summary = pd.DataFrame(sample_summary_rows)
    network_df.to_csv(args.outdir / "ivy_gap_anatomic_network_tracing_per_sample_top10.csv", index=False, encoding="utf-8-sig")
    region_df.to_csv(args.outdir / "ivy_gap_anatomic_region_tracing_per_sample_top15.csv", index=False, encoding="utf-8-sig")
    sample_summary.to_csv(args.outdir / "ivy_gap_anatomic_tracing_sample_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "ivy_gap_anatomic_tracing_meta.json").write_text(
        json.dumps(trace_meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    distribution_rows: list[dict[str, Any]] = []
    add_distribution_rows(distribution_rows, sample_summary, "structure_acronym", "network_top1", "network_top1")
    add_top3_membership_rows(distribution_rows, sample_summary, "structure_acronym", "network_top3", "network_top3_membership")
    add_distribution_rows(distribution_rows, sample_summary, "structure_acronym", "region_top1", "region_top1")
    add_top3_membership_rows(distribution_rows, sample_summary, "structure_acronym", "region_top3", "region_top3_membership")
    add_distribution_rows(
        distribution_rows,
        sample_summary,
        "structure_acronym",
        "resolution_group_top1",
        "resolution_group_top1",
    )
    add_top3_membership_rows(
        distribution_rows,
        sample_summary,
        "structure_acronym",
        "resolution_group_top3",
        "resolution_group_top3_membership",
    )
    add_distribution_rows(
        distribution_rows,
        sample_summary,
        "structure_acronym",
        "top1_resolution_tier",
        "top1_resolution_tier",
    )
    distribution = pd.DataFrame(distribution_rows)
    distribution.to_csv(args.outdir / "ivy_gap_anatomic_structure_prediction_distributions.csv", index=False, encoding="utf-8-sig")

    structure_summary = (
        sample_summary.groupby("structure_acronym")
        .agg(
            n_samples=("sample_id", "count"),
            n_tumors=("tumor_name", pd.Series.nunique),
            dominant_network_top1=("network_top1", lambda s: s.value_counts().index[0] if len(s) else ""),
            dominant_region_top1=("region_top1", lambda s: s.value_counts().index[0] if len(s) else ""),
            dominant_resolution_group_top1=(
                "resolution_group_top1",
                lambda s: s.value_counts().index[0] if len(s) else "",
            ),
            low_resolution_top1_fraction=(
                "top1_resolution_tier",
                lambda s: float((s == "low_resolution").mean()) if len(s) else 0.0,
            ),
            manual_review_top1_fraction=(
                "top1_manual_review_recommended",
                lambda s: float(pd.Series(s).astype(bool).mean()) if len(s) else 0.0,
            ),
        )
        .reset_index()
    )
    structure_summary.to_csv(args.outdir / "ivy_gap_anatomic_structure_summary.csv", index=False, encoding="utf-8-sig")

    plot_distribution(
        sample_summary,
        "network_top1",
        "Ivy GAP structures: Network Top1 distribution",
        args.outdir / "ivy_gap_structure_network_top1_distribution.png",
    )
    plot_distribution(
        sample_summary,
        "region_top1",
        "Ivy GAP structures: Exact Region Top1 distribution",
        args.outdir / "ivy_gap_structure_region_top1_distribution.png",
    )
    plot_distribution(
        sample_summary,
        "resolution_group_top1",
        "Ivy GAP structures: Resolution Group Top1 distribution",
        args.outdir / "ivy_gap_structure_resolution_group_top1_distribution.png",
    )

    run_summary = {
        "n_samples": int(len(sample_summary)),
        "n_structures": int(sample_summary["structure_acronym"].nunique()),
        "matrix": str(args.matrix),
        "metadata": str(args.metadata),
        "db": str(args.db),
        "atlas_id": int(args.atlas_id),
        "interpretation": (
            "Ivy GAP structure labels are GBM microanatomic structures, not Bo2023 normal brain regions. "
            "These outputs summarize prediction distributions by true Ivy structure and do not report normal-region accuracy."
        ),
        "outputs": {
            "sample_summary": str(args.outdir / "ivy_gap_anatomic_tracing_sample_summary.csv"),
            "network_per_sample_top10": str(args.outdir / "ivy_gap_anatomic_network_tracing_per_sample_top10.csv"),
            "region_per_sample_top15": str(args.outdir / "ivy_gap_anatomic_region_tracing_per_sample_top15.csv"),
            "structure_distribution": str(args.outdir / "ivy_gap_anatomic_structure_prediction_distributions.csv"),
            "structure_summary": str(args.outdir / "ivy_gap_anatomic_structure_summary.csv"),
        },
    }
    (args.outdir / "ivy_gap_anatomic_tracing_run_summary.json").write_text(
        json.dumps(run_summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(run_summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
