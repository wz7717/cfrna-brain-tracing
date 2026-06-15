#!/usr/bin/env python
"""Summarize reference atlas download and inspection quality."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference_atlases"
OUT = ROOT / "data" / "processed_reference" / "00_qc"


def classify_file(path: str, first_cols: str = "") -> str:
    text = f"{path} {first_cols}".lower()
    name = Path(path).name.lower()
    if any(x in name for x in ["manifest", "inspection", "failed_downloads"]):
        return "administrative"
    if any(x in text for x in ["microarrayexpression", "rnaseqtpm", "rnaseqcounts", "gene_median_tpm", "ntpm", "ncpm", "gene_count_matrix"]):
        return "expression_matrix"
    if any(x in text for x in ["sampleannot", "metadata", "annotation", "ontology", "probes", "genes.csv", "clusters", "cluster_types", "query.csv", "data_files"]):
        return "metadata_annotation"
    if name.startswith("transcript_"):
        return "transcript_level_large_table"
    return "unknown"


def classify_failed(row: pd.Series) -> str:
    url = str(row.get("source_url", "")).lower()
    name = str(row.get("file_name", "")).lower()
    reason = str(row.get("reason", "")).lower()
    if "transcript" in url or "h5ad" in url or "loom" in url or "expression_matrices" in url or "mapmycells" in url:
        return "needs_--large"
    if "braincellatlas" in url or "cellxgene" in url or "manual" in name:
        return "needs_manual_download"
    if "404" in reason or "not found" in reason:
        return "url_invalid_replace_entry"
    if "sra" in url or any(x in name for x in [".fastq", ".bam", ".cram", ".bigwig", ".bw"]):
        return "ignorable_raw"
    return "review"


def atlas_from_path(path: str) -> str:
    low = path.lower()
    if "gtex_v8" in low:
        return "GTEx_v8_bulk_RNAseq"
    if "allen_human_brain_atlas" in low:
        return "Allen_Human_Brain_Atlas_Hawrylycz2012"
    if "human_protein_atlas_brain" in low:
        return "Human_Protein_Atlas_Brain_Sjostedt2020"
    if "garza2023" in low:
        return "Garza2023_human_TBI_snRNAseq"
    if "allen_aging_dementia_tbi" in low:
        return "Allen_Aging_Dementia_TBI_GSE104687"
    if "siletti" in low:
        return "Siletti2023_adult_human_whole_brain"
    if "brain_cell_atlas" in low:
        return "Chen2024_Brain_Cell_Atlas"
    if "hodge2019" in low:
        return "Hodge2019_human_cortex_MTG"
    return "unknown"


def md_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_None._"
    view = df.head(max_rows).copy() if max_rows else df.copy()
    view = view.fillna("").astype(str)
    cols = list(view.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in view.iterrows():
        vals = [str(row[c]).replace("|", "\\|").replace("\n", " ") for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    if max_rows and len(df) > max_rows:
        lines.append(f"\n_Only first {max_rows} of {len(df)} rows shown._")
    return "\n".join(lines)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(BASE / "00_manifest" / "reference_atlas_manifest.tsv", sep="\t", dtype=str).fillna("")
    inspection = pd.read_csv(BASE / "00_manifest" / "reference_file_inspection.tsv", sep="\t", dtype=str).fillna("")
    failed_path = BASE / "99_logs" / "failed_downloads.tsv"
    failed = pd.read_csv(failed_path, sep="\t", dtype=str).fillna("") if failed_path.exists() else pd.DataFrame()

    files = []
    for path in BASE.rglob("*"):
        if path.is_file() and ".part" not in path.name:
            files.append({"file_path": str(path), "atlas_name": atlas_from_path(str(path)), "file_size_bytes": path.stat().st_size})
    local = pd.DataFrame(files)

    status_counts = manifest.pivot_table(index="atlas_name", columns="download_status", values="file_name", aggfunc="count", fill_value=0)
    size_counts = local.groupby("atlas_name", as_index=False).agg(local_file_count=("file_path", "count"), local_total_size_bytes=("file_size_bytes", "sum"))
    summary = size_counts.merge(status_counts.reset_index(), on="atlas_name", how="outer").fillna(0)
    summary["failed_url_count"] = 0
    if not failed.empty:
        fc = failed.groupby("atlas_name").size().rename("failed_url_count").reset_index()
        summary = summary.drop(columns=["failed_url_count"]).merge(fc, on="atlas_name", how="left").fillna({"failed_url_count": 0})

    inspection["atlas_name"] = inspection["file_path"].map(atlas_from_path)
    inspection["qc_file_class"] = [classify_file(p, c) for p, c in zip(inspection["file_path"], inspection.get("first_10_columns", ""))]
    inspection["row_count_num"] = pd.to_numeric(inspection["row_count"], errors="coerce")
    inspection["column_count_num"] = pd.to_numeric(inspection["column_count"], errors="coerce")
    inspection["qc_flag"] = ""
    inspection.loc[(inspection["qc_file_class"] == "expression_matrix") & (inspection["column_count_num"] < 2), "qc_flag"] = "abnormal_column_count"
    inspection.loc[(inspection["qc_file_class"] == "expression_matrix") & (inspection["row_count_num"] < 1000), "qc_flag"] = "abnormal_row_count"
    inspection.loc[(inspection["qc_file_class"] == "expression_matrix") & (inspection["whether_gene_symbol_exists"].astype(str) != "True"), "qc_flag"] += ";missing_gene_symbol"
    inspection.loc[(inspection["qc_file_class"] == "expression_matrix") & (inspection["whether_ensembl_id_exists"].astype(str) != "True"), "qc_flag"] += ";missing_ensembl_id"

    failed_summary = pd.DataFrame()
    if not failed.empty:
        failed["failure_category"] = failed.apply(classify_failed, axis=1)
        failed_summary = failed.groupby(["failure_category", "atlas_name"], as_index=False).size()

    summary.to_csv(OUT / "download_qc_summary.tsv", sep="\t", index=False)
    inspection.to_csv(OUT / "download_qc_file_details.tsv", sep="\t", index=False)
    if not failed.empty:
        failed.to_csv(OUT / "failed_downloads_classified.tsv", sep="\t", index=False)

    lines = [
        "# Download QC Summary",
        "",
        f"- Manifest rows: {len(manifest)}",
        f"- Inspected files: {len(inspection)}",
        f"- Local files under reference_atlases: {len(local)}",
        f"- Failed URL records: {len(failed)}",
        "",
        "## Atlas Summary",
        "",
        md_table(summary),
        "",
        "## Expression Matrices",
        "",
        md_table(inspection[inspection["qc_file_class"] == "expression_matrix"][["atlas_name", "file_path", "row_count", "column_count", "qc_flag"]], 120),
        "",
        "## Metadata / Annotation Files",
        "",
        md_table(inspection[inspection["qc_file_class"] == "metadata_annotation"][["atlas_name", "file_path", "row_count", "column_count"]], 80),
        "",
        "## Transcript-Level Large Tables Excluded From V1",
        "",
        md_table(inspection[inspection["file_path"].str.contains("transcript_", case=False, na=False)][["atlas_name", "file_path", "row_count", "column_count"]]),
        "",
        "## Files Missing Gene IDs",
        "",
        md_table(inspection[(inspection["qc_file_class"] == "expression_matrix") & (inspection["qc_flag"].str.contains("missing", na=False))][["atlas_name", "file_path", "qc_flag"]], 120),
        "",
        "## Failed Downloads",
        "",
        md_table(failed_summary) if not failed_summary.empty else "No failed downloads recorded.",
        "",
    ]
    (OUT / "download_qc_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT / 'download_qc_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
