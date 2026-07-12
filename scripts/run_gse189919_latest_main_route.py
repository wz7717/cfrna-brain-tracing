#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import ROUTE_NAME, trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402


DATA_DIR = ROOT / "data" / "external_validation" / "GSE189919"
DEFAULT_OUTDIR = ROOT / "results" / "gse189919_latest_main_route_20260708"
DEFAULT_DB_PATH = ROOT / "cfrna_source_tracing.db"


def sample_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def write_json(path: Path, payload: dict[str, Any] | list[dict[str, Any]]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_default) + "\n",
        encoding="utf-8",
    )


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows._"
    text = frame.copy()
    text = text.astype(object).where(pd.notna(text), "")
    columns = [str(column) for column in text.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in text.itertuples(index=False):
        values = [str(value).replace("\n", " ") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def read_soft_metadata(path: Path) -> pd.DataFrame:
    raw = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    rows: list[dict[str, Any]] = []
    for block in raw.split("^SAMPLE = ")[1:]:
        accession = re.search(r"!Sample_geo_accession = (.*)", block)
        title = re.search(r"!Sample_title = (.*)", block)
        source = re.search(r"!Sample_source_name_ch1 = (.*)", block)
        characteristics: dict[str, str] = {}
        for item in re.findall(r"!Sample_characteristics_ch1 = (.*)", block):
            if ":" in item:
                key, value = item.split(":", 1)
                characteristics[key.strip().lower()] = value.strip()
        disease = characteristics.get("disease", "")
        rows.append(
            {
                "geo_accession": accession.group(1) if accession else "",
                "geo_title": title.group(1) if title else "",
                "sample_key": sample_key(title.group(1) if title else ""),
                "group": "MB" if disease == "Medulloblastoma" else "Normal",
                "molecular_subgroup": characteristics.get("sub group", "not available"),
                "biofluid": source.group(1) if source else "",
            }
        )
    return pd.DataFrame(rows)


def read_expression(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(data_dir / "GSE189919_count.csv.gz", index_col="Geneid")
    tpm = pd.read_csv(data_dir / "GSE189919_tpm_count.csv.gz", index_col="Geneid")
    for frame in (counts, tpm):
        frame.index = frame.index.astype(str).str.strip()
        frame.columns = frame.columns.astype(str)
    counts = counts.apply(pd.to_numeric, errors="raise")
    tpm = tpm.apply(pd.to_numeric, errors="raise")

    if not counts.index.equals(tpm.index):
        raise ValueError("Count and TPM gene orders differ")
    if not counts.columns.equals(tpm.columns):
        raise ValueError("Count and TPM sample orders differ")
    if counts.index.duplicated().any() or tpm.index.duplicated().any():
        raise ValueError("Duplicated gene symbols in official matrix")
    if counts.isna().any().any() or tpm.isna().any().any():
        raise ValueError("Missing expression value")
    if (counts < 0).any().any() or (tpm < 0).any().any():
        raise ValueError("Negative expression value")

    soft = read_soft_metadata(data_dir / "GSE189919_family.soft.gz")
    by_key = soft.set_index("sample_key")
    matrix_keys = pd.Index([sample_key(column) for column in counts.columns])
    missing = matrix_keys.difference(by_key.index)
    extra = by_key.index.difference(matrix_keys)
    if len(missing) or len(extra):
        raise ValueError(
            f"Matrix/SOFT mismatch: missing metadata={missing.tolist()}, "
            f"extra metadata={extra.tolist()}"
        )
    metadata = by_key.loc[matrix_keys].reset_index(drop=True)
    metadata.insert(0, "sample_id", counts.columns)
    return counts.astype(float), tpm.astype(float), metadata


def sample_frame_from_counts(counts: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": counts.index.astype(str),
            "read_count": counts[sample_id].to_numpy(dtype=float),
        }
    )


def top_networks(network_out: dict[str, Any], k: int = 3) -> str:
    return " | ".join(str(row.get("network_id", "")) for row in network_out.get("results", [])[:k])


def top_regions(region_out: dict[str, Any], key: str, k: int = 3) -> str:
    return " | ".join(str(row.get(key, "")) for row in region_out.get("results", [])[:k])


def run(args: argparse.Namespace) -> int:
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    counts, tpm, metadata = read_expression(args.data_dir)
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("Zero-sized count library")

    sample_rows: list[dict[str, Any]] = []
    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    group_rows: list[dict[str, Any]] = []
    full_records: list[dict[str, Any]] = []

    for i, meta_row in enumerate(metadata.itertuples(index=False), start=1):
        sample_id = str(meta_row.sample_id)
        expression = sample_frame_from_counts(counts, sample_id)
        network_out = trace_network_expression(
            expression,
            min_overlap_fraction=args.min_network_overlap,
            project_to_vsd=True,
        )
        region_out = trace_bo2023_secondary_regions(
            expression,
            network_out,
            str(args.db_path),
            int(args.atlas_id),
            topk=int(args.topk_regions),
        )

        network_meta = network_out.get("meta", {})
        region_meta = region_out.get("meta", {})
        projection_meta = network_meta.get("reference_projection", {})
        resolution_meta = region_meta.get("region_resolution_annotation", {})
        network_ranked = network_out.get("results", [])
        region_ranked = region_out.get("results", [])

        sample_rows.append(
            {
                "sample_id": sample_id,
                "geo_accession": meta_row.geo_accession,
                "geo_title": meta_row.geo_title,
                "group": meta_row.group,
                "molecular_subgroup": meta_row.molecular_subgroup,
                "biofluid": meta_row.biofluid,
                "route": ROUTE_NAME,
                "input_route": "official_raw_counts_to_logcpm",
                "network_traceability": network_meta.get("traceability"),
                "region_traceability": region_meta.get("traceability"),
                "library_size": float(library_sizes[sample_id]),
                "n_detected_count_genes": int((counts[sample_id] > 0).sum()),
                "n_detected_tpm_genes": int((tpm[sample_id] > 0).sum()),
                "n_network_model_genes": network_meta.get("n_model_genes"),
                "n_network_overlap_genes": network_meta.get("n_overlap_genes"),
                "network_overlap_fraction": network_meta.get("overlap_fraction"),
                "n_projector_genes": projection_meta.get("n_projector_genes"),
                "n_input_projector_overlap_genes": projection_meta.get("n_input_projector_overlap_genes"),
                "projector_output_scale": projection_meta.get("output_scale"),
                "network_top1": network_ranked[0].get("network_id") if network_ranked else "",
                "network_top2": network_ranked[1].get("network_id") if len(network_ranked) > 1 else "",
                "network_top3": network_ranked[2].get("network_id") if len(network_ranked) > 2 else "",
                "network_top3_beam": top_networks(network_out, 3),
                "network_pairwise_switched": bool(
                    network_meta.get("pairwise_rescue", {}).get("switched", False)
                ),
                "n_candidate_regions": region_meta.get("n_candidate_regions"),
                "n_region_overlap_genes": region_meta.get("n_overlap_genes"),
                "n_local_candidate_genes": region_meta.get("n_local_candidate_genes"),
                "region_top1": region_ranked[0].get("region_id") if region_ranked else "",
                "region_top2": region_ranked[1].get("region_id") if len(region_ranked) > 1 else "",
                "region_top3": region_ranked[2].get("region_id") if len(region_ranked) > 2 else "",
                "region_top3_list": top_regions(region_out, "region_id", 3),
                "resolution_group_top1": resolution_meta.get("top1_resolution_group"),
                "resolution_tier_top1": resolution_meta.get("top1_resolution_tier"),
                "group_plausibility_tier_top1": resolution_meta.get("top1_group_plausibility_tier"),
                "manual_review_recommended": resolution_meta.get("manual_review_recommended"),
                "region_error": region_meta.get("error", ""),
            }
        )

        for row in network_ranked:
            network_rows.append(
                {
                    "sample_id": sample_id,
                    "group": meta_row.group,
                    "rank": row.get("rank"),
                    "network_id": row.get("network_id"),
                    "score": row.get("score"),
                    "confidence": row.get("confidence"),
                }
            )
        for row in region_ranked:
            item = {"sample_id": sample_id, "group": meta_row.group, **dict(row)}
            region_rows.append(item)
        for row in resolution_meta.get("group_ranking", []) or []:
            item = {"sample_id": sample_id, "group": meta_row.group, **dict(row)}
            group_rows.append(item)

        full_records.append(
            {
                "sample_id": sample_id,
                "metadata": {
                    "geo_accession": meta_row.geo_accession,
                    "geo_title": meta_row.geo_title,
                    "group": meta_row.group,
                    "molecular_subgroup": meta_row.molecular_subgroup,
                    "biofluid": meta_row.biofluid,
                },
                "network_output": network_out,
                "region_output": region_out,
            }
        )
        print(f"[{i:02d}/{len(metadata):02d}] {sample_id}: {top_networks(network_out, 3)}")

    sample_df = pd.DataFrame(sample_rows)
    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    group_df = pd.DataFrame(group_rows)

    sample_df.to_csv(outdir / "gse189919_latest_main_route_sample_summary.csv", index=False)
    network_df.to_csv(outdir / "gse189919_latest_main_route_network_rankings.csv", index=False)
    region_df.to_csv(outdir / "gse189919_latest_main_route_region_rankings.csv", index=False)
    group_df.to_csv(outdir / "gse189919_latest_main_route_resolution_group_rankings.csv", index=False)
    metadata.to_csv(outdir / "gse189919_sample_metadata.csv", index=False)

    with (outdir / "gse189919_latest_main_route_full_outputs.jsonl").open("wt", encoding="utf-8") as handle:
        for record in full_records:
            handle.write(json.dumps(record, ensure_ascii=False, default=json_default) + "\n")

    group_distribution = (
        sample_df.groupby(["group", "network_top1"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["group", "n"], ascending=[True, False])
    )
    region_distribution = (
        sample_df.groupby(["group", "region_top1"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["group", "n"], ascending=[True, False])
    )
    group_distribution.to_csv(outdir / "gse189919_latest_main_route_network_top1_distribution.csv", index=False)
    region_distribution.to_csv(outdir / "gse189919_latest_main_route_region_top1_distribution.csv", index=False)

    script_snapshot = outdir / "run_gse189919_latest_main_route.py"
    script_snapshot.write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "route": ROUTE_NAME,
        "input_data_dir": str(args.data_dir.resolve()),
        "outdir": str(outdir),
        "db_path": str(args.db_path),
        "atlas_id": int(args.atlas_id),
        "samples": int(len(sample_df)),
        "groups": sample_df["group"].value_counts().to_dict(),
        "interpretation_boundary": (
            "GSE189919 CSF RNA-seq has no patient-level anatomical or MRI localization truth in this repository; "
            "results are transfer/projection stress-test outputs, not localization accuracy."
        ),
    }
    report_file_names = sorted(
        {
            *(path.name for path in outdir.iterdir() if path.is_file()),
            "GSE189919_LATEST_MAIN_ROUTE_REPORT.md",
            "manifest.json",
        }
    )

    report = [
        "# GSE189919 latest main-route recalculation",
        "",
        f"- Created UTC: {manifest['created_utc']}",
        f"- Route: `{ROUTE_NAME}`",
        "- Input used for main route: official raw counts, internally converted to logCPM.",
        "- Stage 1: reference-fitted projected-VSD SaleemNetworks Top3 beam.",
        "- Stages 2-3: logCPM-compatible resolution-group and exact-region reranking within the Network Top3 beam.",
        "- Interpretation boundary: GSE189919 has no patient-level anatomical/MRI truth here; these outputs are transfer/projection stress-test results, not localization-accuracy validation.",
        "",
        "## Cohort",
        "",
        markdown_table(sample_df["group"].value_counts().rename_axis("group").reset_index(name="n")),
        "",
        "## Network Top1 distribution",
        "",
        markdown_table(group_distribution),
        "",
        "## Region Top1 distribution",
        "",
        markdown_table(region_distribution),
        "",
        "## Output files",
        "",
        "\n".join(f"- `{name}`" for name in report_file_names),
        "",
    ]
    (outdir / "GSE189919_LATEST_MAIN_ROUTE_REPORT.md").write_text("\n".join(report), encoding="utf-8")
    manifest["files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    write_json(outdir / "manifest.json", manifest)

    print(f"Saved GSE189919 latest main-route outputs to: {outdir}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Recalculate GSE189919 with the current paper main route.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db-path", type=Path, default=DEFAULT_DB_PATH)
    parser.add_argument("--atlas-id", type=int, default=1)
    parser.add_argument("--topk-regions", type=int, default=30)
    parser.add_argument("--min-network-overlap", type=float, default=0.50)
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
