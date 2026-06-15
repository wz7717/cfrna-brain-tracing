#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "data" / "tcga_brain_tumor_expression"
GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"


def gdc_post(endpoint: str, payload: dict[str, Any]) -> dict[str, Any]:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        endpoint,
        data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode("utf-8"))


def gdc_filters(projects: list[str], sample_types: list[str]) -> dict[str, Any]:
    content: list[dict[str, Any]] = [
        {
            "op": "in",
            "content": {"field": "cases.project.project_id", "value": projects},
        },
        {
            "op": "=",
            "content": {"field": "data_category", "value": "Transcriptome Profiling"},
        },
        {
            "op": "=",
            "content": {"field": "data_type", "value": "Gene Expression Quantification"},
        },
        {
            "op": "=",
            "content": {"field": "experimental_strategy", "value": "RNA-Seq"},
        },
        {
            "op": "=",
            "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"},
        },
        {
            "op": "=",
            "content": {"field": "access", "value": "open"},
        },
    ]
    if sample_types:
        content.append(
            {
                "op": "in",
                "content": {"field": "cases.samples.sample_type", "value": sample_types},
            }
        )
    return {"op": "and", "content": content}


def query_gdc_files(projects: list[str], sample_types: list[str]) -> list[dict[str, Any]]:
    fields = [
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "cases.case_id",
        "cases.submitter_id",
        "cases.project.project_id",
        "cases.samples.sample_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.diagnoses.primary_diagnosis",
        "cases.diagnoses.tumor_grade",
        "cases.diagnoses.classification_of_tumor",
    ]
    payload = {
        "filters": gdc_filters(projects, sample_types),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": 2000,
        "sort": "cases.project.project_id:asc,cases.submitter_id:asc",
    }
    data = gdc_post(GDC_FILES, payload)
    return data["data"]["hits"]


def flatten_file_record(record: dict[str, Any]) -> dict[str, Any]:
    case = record.get("cases", [{}])[0] or {}
    sample = (case.get("samples") or [{}])[0] or {}
    diagnosis = (case.get("diagnoses") or [{}])[0] or {}
    project = case.get("project") or {}
    return {
        "file_id": record.get("file_id", ""),
        "file_name": record.get("file_name", ""),
        "md5sum": record.get("md5sum", ""),
        "file_size": record.get("file_size", ""),
        "project_id": project.get("project_id", ""),
        "case_id": case.get("case_id", ""),
        "case_submitter_id": case.get("submitter_id", ""),
        "sample_id": sample.get("sample_id", ""),
        "sample_submitter_id": sample.get("submitter_id", ""),
        "sample_type": sample.get("sample_type", ""),
        "primary_diagnosis": diagnosis.get("primary_diagnosis", ""),
        "tumor_grade": diagnosis.get("tumor_grade", ""),
        "classification_of_tumor": diagnosis.get("classification_of_tumor", ""),
    }


def download_file(file_id: str, destination: Path, retries: int = 3) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = f"{GDC_DATA}/{file_id}"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=300) as resp, destination.open("wb") as out:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            return
        except (urllib.error.URLError, TimeoutError):
            if attempt == retries:
                raise
            time.sleep(2 * attempt)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_star_counts(path: Path, value_col: str) -> pd.Series:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        df = pd.read_csv(fh, sep="\t", comment="#")
    required = {"gene_id", "gene_name", value_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} missing columns: {sorted(missing)}")
    df = df[df["gene_id"].astype(str).str.startswith("ENSG")].copy()
    df["gene_symbol"] = df["gene_name"].astype(str).str.strip()
    df = df[df["gene_symbol"].ne("") & df["gene_symbol"].ne("nan")]
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce").fillna(0.0)
    return df.groupby("gene_symbol", sort=True)[value_col].mean()


def make_file_level_column(row: Any) -> str:
    sample = str(row.sample_submitter_id or row.case_submitter_id)
    file_short = str(row.file_id)[:8]
    return f"{sample}|{file_short}"


def write_file_level_matrix(
    metadata: pd.DataFrame,
    raw_dir: Path,
    out_path: Path,
    value_col: str,
) -> pd.DataFrame:
    series_by_sample: dict[str, pd.Series] = {}
    for row in metadata.itertuples(index=False):
        sample_name = make_file_level_column(row)
        path = raw_dir / str(row.file_name)
        series_by_sample[sample_name] = read_star_counts(path, value_col)
    matrix = pd.DataFrame(series_by_sample).fillna(0.0)
    matrix.index.name = "gene_symbol"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out_path, sep="\t", encoding="utf-8")
    return matrix


def write_sample_level_matrix(file_matrix: pd.DataFrame, metadata: pd.DataFrame, out_path: Path, reducer: str) -> pd.DataFrame:
    column_to_sample = {
        make_file_level_column(row): str(row.sample_submitter_id or row.case_submitter_id)
        for row in metadata.itertuples(index=False)
    }
    groups: dict[str, list[str]] = {}
    for column, sample in column_to_sample.items():
        groups.setdefault(sample, []).append(column)
    data: dict[str, pd.Series] = {}
    for sample, columns in sorted(groups.items()):
        sub = file_matrix[columns]
        if reducer == "mean":
            data[sample] = sub.mean(axis=1)
        elif reducer == "sum":
            data[sample] = sub.sum(axis=1)
        else:
            raise ValueError(f"unknown reducer: {reducer}")
    matrix = pd.DataFrame(data).fillna(0.0)
    matrix.index.name = "gene_symbol"
    matrix.to_csv(out_path, sep="\t", encoding="utf-8")
    return matrix


def write_project_centroids(sample_matrix: pd.DataFrame, metadata: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    sample_project = (
        metadata[["sample_submitter_id", "project_id"]]
        .drop_duplicates("sample_submitter_id")
        .set_index("sample_submitter_id")["project_id"]
        .to_dict()
    )
    data = {}
    for project in sorted(set(sample_project.values())):
        columns = [col for col in sample_matrix.columns if sample_project.get(col) == project]
        data[project] = sample_matrix[columns].mean(axis=1)
    centroids = pd.DataFrame(data).fillna(0.0)
    centroids.index.name = "gene_symbol"
    centroids.to_csv(out_path, sep="\t", encoding="utf-8")
    return centroids


def main() -> int:
    parser = argparse.ArgumentParser(description="Download TCGA-GBM/LGG GDC STAR-count expression matrices.")
    parser.add_argument("--projects", nargs="+", default=["TCGA-GBM", "TCGA-LGG"])
    parser.add_argument("--sample-types", nargs="+", default=["Primary Tumor"])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--max-files", type=int, default=None, help="Limit files for smoke testing.")
    parser.add_argument("--dry-run", action="store_true", help="Only query manifest and metadata.")
    parser.add_argument("--rebuild-only", action="store_true", help="Use existing manifest/raw files and rebuild matrices.")
    parser.add_argument("--skip-existing", action="store_true", default=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.outdir / "gdc_star_counts_raw"
    manifest_path = args.outdir / "tcga_gbm_lgg_gdc_star_counts_manifest.csv"
    if args.rebuild_only:
        metadata = pd.read_csv(manifest_path)
    else:
        records = query_gdc_files(args.projects, args.sample_types)
        metadata = pd.DataFrame([flatten_file_record(x) for x in records])
    metadata = metadata.sort_values(["project_id", "case_submitter_id", "sample_submitter_id", "file_id"]).reset_index(drop=True)
    if args.max_files is not None:
        metadata = metadata.head(max(1, int(args.max_files))).copy()
    metadata.to_csv(manifest_path, index=False, encoding="utf-8-sig")

    if args.dry_run:
        summary = {
            "projects": args.projects,
            "sample_types": args.sample_types,
            "n_files": int(len(metadata)),
            "projects_observed": metadata["project_id"].value_counts().to_dict() if len(metadata) else {},
            "sample_types_observed": metadata["sample_type"].value_counts().to_dict() if len(metadata) else {},
        }
        (args.outdir / "tcga_gbm_lgg_download_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0

    if not args.rebuild_only:
        downloaded = 0
        for row in metadata.itertuples(index=False):
            destination = raw_dir / str(row.file_name)
            if destination.exists() and args.skip_existing:
                expected = str(row.md5sum or "")
                if not expected or md5(destination) == expected:
                    continue
            download_file(str(row.file_id), destination)
            expected = str(row.md5sum or "")
            if expected and md5(destination) != expected:
                raise ValueError(f"MD5 mismatch for {destination.name}")
            downloaded += 1
            print(f"downloaded {downloaded}: {destination.name}", flush=True)

    missing_raw = [str(row.file_name) for row in metadata.itertuples(index=False) if not (raw_dir / str(row.file_name)).exists()]
    if missing_raw:
        raise FileNotFoundError(f"{len(missing_raw)} raw files are missing; first missing: {missing_raw[0]}")

    tpm_file = write_file_level_matrix(
        metadata,
        raw_dir,
        args.outdir / "tcga_gbm_lgg_primary_tumor_tpm_unstranded_file_level.tsv",
        "tpm_unstranded",
    )
    tpm = write_sample_level_matrix(
        tpm_file,
        metadata,
        args.outdir / "tcga_gbm_lgg_primary_tumor_tpm_unstranded_sample_mean.tsv",
        "mean",
    )
    counts_file = write_file_level_matrix(
        metadata,
        raw_dir,
        args.outdir / "tcga_gbm_lgg_primary_tumor_unstranded_counts_file_level.tsv",
        "unstranded",
    )
    counts = write_sample_level_matrix(
        counts_file,
        metadata,
        args.outdir / "tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv",
        "sum",
    )
    centroids = write_project_centroids(
        tpm,
        metadata,
        args.outdir / "tcga_gbm_lgg_primary_tumor_project_mean_tpm.tsv",
    )
    summary = {
        "projects": args.projects,
        "sample_types": args.sample_types,
        "n_files": int(len(metadata)),
        "n_gene_symbols_tpm": int(tpm.shape[0]),
        "n_gene_symbols_counts": int(counts.shape[0]),
        "n_files": int(len(metadata)),
        "n_file_level_columns": int(tpm_file.shape[1]),
        "n_samples": int(tpm.shape[1]),
        "n_project_centroids": int(centroids.shape[1]),
        "duplicate_file_count_over_samples": int(len(metadata) - tpm.shape[1]),
        "outputs": {
            "metadata": str(manifest_path),
            "tpm_file_level_matrix": str(args.outdir / "tcga_gbm_lgg_primary_tumor_tpm_unstranded_file_level.tsv"),
            "tpm_sample_mean_matrix": str(args.outdir / "tcga_gbm_lgg_primary_tumor_tpm_unstranded_sample_mean.tsv"),
            "count_file_level_matrix": str(args.outdir / "tcga_gbm_lgg_primary_tumor_unstranded_counts_file_level.tsv"),
            "count_sample_sum_matrix": str(args.outdir / "tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv"),
            "project_mean_tpm_matrix": str(args.outdir / "tcga_gbm_lgg_primary_tumor_project_mean_tpm.tsv"),
        },
    }
    (args.outdir / "tcga_gbm_lgg_download_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
