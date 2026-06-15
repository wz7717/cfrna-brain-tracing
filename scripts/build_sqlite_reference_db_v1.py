#!/usr/bin/env python
"""Build a lightweight SQLite database for Streamlit reference queries."""

from __future__ import annotations

import re
import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V1 = ROOT / "data" / "processed_reference" / "v1"
REF = ROOT / "data" / "reference_atlases"
DB = V1 / "cfrna_brain_tracing_reference_v1.sqlite"


def clean_col(name: object) -> str:
    text = str(name).strip()
    text = re.sub(r"[^0-9A-Za-z_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    if not text:
        text = "col"
    if text[0].isdigit():
        text = "c_" + text
    return text[:120]


def make_unique(cols: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for col in cols:
        base = clean_col(col)
        idx = seen.get(base, 0)
        seen[base] = idx + 1
        out.append(base if idx == 0 else f"{base}_{idx}")
    return out


def read_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = make_unique(list(df.columns))
    return df


def write_table(conn: sqlite3.Connection, table: str, df: pd.DataFrame) -> None:
    df.to_sql(table, conn, if_exists="replace", index=False, chunksize=10_000)


def ensure_source_fields(df: pd.DataFrame, source_atlas: str, doi: str, use: str) -> pd.DataFrame:
    out = df.copy()
    if "source_atlas" not in out.columns:
        out["source_atlas"] = source_atlas
    if "doi" not in out.columns:
        out["doi"] = doi
    if "recommended_use" not in out.columns:
        out["recommended_use"] = use
    return out


def create_index(conn: sqlite3.Connection, table: str, column: str) -> None:
    cols = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
    if column in cols:
        conn.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_{column} ON {table} ({column})")


def main() -> int:
    V1.mkdir(parents=True, exist_ok=True)
    if DB.exists():
        DB.unlink()
    conn = sqlite3.connect(DB)
    manifest = read_tsv(REF / "00_manifest" / "reference_atlas_manifest.tsv")
    inspection = read_tsv(REF / "00_manifest" / "reference_file_inspection.tsv")
    write_table(conn, "atlas_manifest", manifest)
    write_table(conn, "file_inspection", inspection)

    table_specs = [
        ("gtex_median_tpm", V1 / "peripheral_background_gtex_median_tpm.tsv", "GTEx_v8_bulk_RNAseq", "10.1126/science.aaz1776", "peripheral_background"),
        ("brain_celltype_reference", V1 / "celltype_reference_hpa_single_nuclei_brain.tsv", "Human_Protein_Atlas_Brain_Sjostedt2020", "10.1126/science.aay5947", "celltype_reference"),
        ("injury_state_reference", V1 / "injury_state_reference_garza2023_tbi.tsv", "Garza2023_human_TBI_snRNAseq", "10.1016/j.celrep.2023.113395", "injury_state_reference"),
        ("gene_id_map", V1 / "reference_gene_universe.tsv", "processed_reference_v1", "", "gene_id_map"),
    ]
    for table, path, source, doi, use in table_specs:
        if path.exists():
            write_table(conn, table, ensure_source_fields(read_tsv(path), source, doi, use))

    region_frames = []
    region_hpa = V1 / "brain_region_reference_hpa.tsv"
    if region_hpa.exists():
        df = ensure_source_fields(read_tsv(region_hpa), "Human_Protein_Atlas_Brain_Sjostedt2020", "10.1126/science.aay5947", "brain_region_reference")
        df["reference_subtype"] = "hpa_brain_region"
        region_frames.append(df)
    region_allen = V1 / "brain_region_reference_allen_hba.tsv"
    if region_allen.exists():
        df = ensure_source_fields(read_tsv(region_allen), "Allen_Human_Brain_Atlas_Hawrylycz2012", "10.1038/nature11405", "brain_region_reference")
        df["reference_subtype"] = "allen_hba_rnaseq"
        region_frames.append(df)
    if region_frames:
        write_table(conn, "brain_region_reference", pd.concat(region_frames, ignore_index=True, sort=False))

    aging = V1 / "injury_state_reference_allen_aging_tbi.tsv"
    if aging.exists():
        df = read_tsv(aging)
        df.to_sql("injury_state_reference_metadata", conn, if_exists="replace", index=False, chunksize=10_000)

    marker_dir = V1 / "marker_candidates"
    marker_frames = []
    marker_map = {
        "brain_enriched_genes_vs_gtex.tsv": "brain_enriched",
        "brain_region_marker_candidates.tsv": "brain_region",
        "brain_celltype_marker_candidates.tsv": "brain_celltype",
        "tbi_injury_response_candidates.tsv": "tbi_injury",
    }
    for fname, marker_type in marker_map.items():
        path = marker_dir / fname
        if path.exists():
            df = read_tsv(path)
            df["marker_type"] = marker_type
            marker_frames.append(df)
    if marker_frames:
        write_table(conn, "marker_candidates", pd.concat(marker_frames, ignore_index=True, sort=False))
    contam = marker_dir / "blood_immune_rbc_contamination_genes.tsv"
    if contam.exists():
        write_table(conn, "contamination_markers", read_tsv(contam))

    for table in [
        "atlas_manifest",
        "file_inspection",
        "gtex_median_tpm",
        "brain_region_reference",
        "brain_celltype_reference",
        "injury_state_reference",
        "marker_candidates",
        "contamination_markers",
        "gene_id_map",
    ]:
        create_index(conn, table, "gene_symbol")
        create_index(conn, table, "source_atlas")
        create_index(conn, table, "region")
        create_index(conn, table, "celltype")
        create_index(conn, table, "atlas_name")
    conn.commit()
    conn.close()
    print(f"Wrote {DB}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
