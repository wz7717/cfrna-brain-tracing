#!/usr/bin/env python
from __future__ import annotations

import argparse
import io
import json
import re
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "data" / "ivy_gap_anatomic_rnaseq"

SAMPLE_DETAILS_URL = "https://glioblastoma.alleninstitute.org/api/v2/gbm/rna_seq_samples_details.csv"
RAW_FPKM_MANIFEST_URL = "https://glioblastoma.alleninstitute.org/rnaseq/raw_fpkm.csv"
GENE_EXPRESSION_ZIP_URL = "https://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873915"

ANATOMIC_REFERENCE_STRUCTURES = {
    "CT-reference-histology",
    "CTmvp-reference-histology",
    "CTpan-reference-histology",
    "IT-reference-histology",
    "LE-reference-histology",
}


def download_url(url: str, destination: Path, retries: int = 3) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
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


def read_url_csv(url: str) -> pd.DataFrame:
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = resp.read()
    return pd.read_csv(io.BytesIO(data))


def safe_name(text: Any) -> str:
    value = str(text).strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def ensure_gene_annotation(outdir: Path, refresh: bool = False) -> pd.DataFrame:
    zip_path = outdir / "gene_expression_matrix_2014-11-25.zip"
    rows_path = outdir / "rows-genes.csv"
    if refresh or not rows_path.exists():
        if refresh or not zip_path.exists():
            download_url(GENE_EXPRESSION_ZIP_URL, zip_path)
        with zipfile.ZipFile(zip_path) as zf:
            with zf.open("rows-genes.csv") as src:
                rows_path.write_bytes(src.read())
    genes = pd.read_csv(rows_path)
    required = {"gene_entrez_id", "gene_symbol"}
    missing = required - set(genes.columns)
    if missing:
        raise ValueError(f"rows-genes.csv missing columns: {sorted(missing)}")
    genes = genes.copy()
    genes["gene_entrez_id"] = pd.to_numeric(genes["gene_entrez_id"], errors="coerce")
    genes = genes.dropna(subset=["gene_entrez_id", "gene_symbol"])
    genes["gene_entrez_id"] = genes["gene_entrez_id"].astype("int64").astype(str)
    genes["gene_symbol"] = genes["gene_symbol"].astype(str).str.strip()
    genes = genes[genes["gene_symbol"].ne("")]
    return genes.drop_duplicates("gene_entrez_id", keep="first")


def profile_path(raw_dir: Path, row: Any) -> Path:
    sample_id = safe_name(row.rna_well)
    structure = safe_name(row.structure_acronym)
    specimen = safe_name(row.specimen_name)
    return raw_dir / f"{sample_id}__{structure}__{specimen}.tsv"


def read_profile(path: Path, value_col: str, entrez_to_symbol: dict[str, str]) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    required = {"gene_id", value_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} missing columns: {sorted(missing)}")
    gene_id = pd.to_numeric(df["gene_id"], errors="coerce")
    df = df.loc[gene_id.notna()].copy()
    df["gene_entrez_id"] = gene_id.loc[gene_id.notna()].astype("int64").astype(str).to_numpy()
    df["gene_symbol"] = df["gene_entrez_id"].map(entrez_to_symbol)
    df = df.dropna(subset=["gene_symbol"])
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce").fillna(0.0)
    return df.groupby("gene_symbol", sort=True)[value_col].mean()


def write_matrix(
    manifest: pd.DataFrame,
    raw_dir: Path,
    out_path: Path,
    value_col: str,
    entrez_to_symbol: dict[str, str],
) -> pd.DataFrame:
    series_by_sample: dict[str, pd.Series] = {}
    for row in manifest.itertuples(index=False):
        sample_id = str(row.rna_well)
        series_by_sample[sample_id] = read_profile(profile_path(raw_dir, row), value_col, entrez_to_symbol)
    matrix = pd.DataFrame(series_by_sample).fillna(0.0)
    matrix.index.name = "gene_symbol"
    matrix.to_csv(out_path, sep="\t", encoding="utf-8")
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download Ivy GAP Anatomic Structures RNA-seq sample metadata and gene-level FPKM/TPM matrices."
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--refresh", action="store_true", help="Re-download metadata, manifests, and profiles.")
    parser.add_argument("--dry-run", action="store_true", help="Only write filtered metadata/manifest, do not download profiles.")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.outdir / "raw_gene_level_profiles"
    raw_dir.mkdir(parents=True, exist_ok=True)

    details = read_url_csv(SAMPLE_DETAILS_URL)
    raw_manifest = read_url_csv(RAW_FPKM_MANIFEST_URL)
    details.to_csv(args.outdir / "ivy_gap_rna_seq_samples_details_all.csv", index=False, encoding="utf-8-sig")
    raw_manifest.to_csv(args.outdir / "ivy_gap_raw_fpkm_manifest_all.csv", index=False, encoding="utf-8-sig")

    anatomic_details = details[details["structure_acronym"].isin(ANATOMIC_REFERENCE_STRUCTURES)].copy()
    anatomic_manifest = raw_manifest[raw_manifest["structure_acronym"].isin(ANATOMIC_REFERENCE_STRUCTURES)].copy()
    anatomic_details.to_csv(
        args.outdir / "ivy_gap_anatomic_structure_sample_metadata_122.csv",
        index=False,
        encoding="utf-8-sig",
    )
    anatomic_manifest.to_csv(
        args.outdir / "ivy_gap_anatomic_structure_raw_fpkm_manifest_122.csv",
        index=False,
        encoding="utf-8-sig",
    )

    summary: dict[str, Any] = {
        "source": "Ivy Glioblastoma Atlas Project",
        "sample_details_url": SAMPLE_DETAILS_URL,
        "raw_fpkm_manifest_url": RAW_FPKM_MANIFEST_URL,
        "gene_expression_zip_url": GENE_EXPRESSION_ZIP_URL,
        "all_metadata_samples": int(len(details)),
        "all_raw_profile_links": int(len(raw_manifest)),
        "anatomic_structure_samples": int(len(anatomic_details)),
        "anatomic_structure_profile_links": int(len(anatomic_manifest)),
        "selected_structures": sorted(ANATOMIC_REFERENCE_STRUCTURES),
        "downloaded_profiles": 0,
        "notes": [
            "BAM files are intentionally not downloaded.",
            "Raw profiles contain Entrez gene IDs; matrices are converted to gene symbols using rows-genes.csv.",
        ],
    }

    if args.dry_run:
        (args.outdir / "ivy_gap_anatomic_rnaseq_download_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0

    genes = ensure_gene_annotation(args.outdir, refresh=args.refresh)
    entrez_to_symbol = dict(zip(genes["gene_entrez_id"].astype(str), genes["gene_symbol"].astype(str)))

    downloaded = 0
    for row in anatomic_manifest.itertuples(index=False):
        path = profile_path(raw_dir, row)
        if args.refresh or not path.exists() or path.stat().st_size == 0:
            download_url(str(row.file_download_link), path)
        downloaded += 1

    fpkm = write_matrix(
        anatomic_manifest,
        raw_dir,
        args.outdir / "ivy_gap_anatomic_structure_fpkm_gene_symbol_matrix.tsv",
        "FPKM",
        entrez_to_symbol,
    )
    tpm = write_matrix(
        anatomic_manifest,
        raw_dir,
        args.outdir / "ivy_gap_anatomic_structure_tpm_gene_symbol_matrix.tsv",
        "TPM",
        entrez_to_symbol,
    )

    summary.update(
        {
            "downloaded_profiles": downloaded,
            "gene_annotation_rows": int(len(genes)),
            "fpkm_matrix_shape": [int(fpkm.shape[0]), int(fpkm.shape[1])],
            "tpm_matrix_shape": [int(tpm.shape[0]), int(tpm.shape[1])],
            "outputs": {
                "metadata": str(args.outdir / "ivy_gap_anatomic_structure_sample_metadata_122.csv"),
                "manifest": str(args.outdir / "ivy_gap_anatomic_structure_raw_fpkm_manifest_122.csv"),
                "raw_profiles_dir": str(raw_dir),
                "fpkm_matrix": str(args.outdir / "ivy_gap_anatomic_structure_fpkm_gene_symbol_matrix.tsv"),
                "tpm_matrix": str(args.outdir / "ivy_gap_anatomic_structure_tpm_gene_symbol_matrix.tsv"),
            },
        }
    )
    (args.outdir / "ivy_gap_anatomic_rnaseq_download_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
