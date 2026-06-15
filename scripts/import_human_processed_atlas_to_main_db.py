#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from signature_builder import build_signature_set


DB_PATH = ROOT / "cfrna_source_tracing.db"
PROCESSED = ROOT / "data" / "processed_reference" / "v1"


ATLAS_SPECS = [
    {
        "file": PROCESSED / "brain_region_reference_hpa.tsv",
        "atlas_name": "human_brain_region_hpa",
        "species": "human",
        "level": "L2",
        "build_version": "v1.0.0",
        "gene_id_type": "symbol",
        "normalization": "log1p_tpm",
        "notes": "Imported from processed human HPA brain-region reference.",
        "prefix": "HPA",
    },
    {
        "file": PROCESSED / "brain_region_reference_allen_hba.tsv",
        "atlas_name": "human_brain_region_allen_hba",
        "species": "human",
        "level": "L2",
        "build_version": "v1.0.0",
        "gene_id_type": "symbol",
        "normalization": "log1p_tpm",
        "notes": "Imported from processed Allen Human Brain Atlas RNA-seq reference.",
        "prefix": "AHBA",
    },
]


META_COLS = {"gene_symbol", "gene_symbol_raw", "ensembl_id", "source_atlas", "doi", "recommended_use"}


def clean_region_id(prefix: str, name: str) -> str:
    text = str(name).strip().upper()
    text = re.sub(r"[^0-9A-Z]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return f"{prefix}_{text}"[:80]


def infer_expression_class(value: float) -> str:
    if value >= 5:
        return "high"
    if value >= 2:
        return "medium"
    return "low"


def ensure_schema(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS atlas_versions (
            atlas_id INTEGER PRIMARY KEY,
            atlas_name TEXT,
            species TEXT,
            level TEXT,
            build_version TEXT,
            gene_id_type TEXT,
            normalization TEXT,
            created_at TEXT DEFAULT (datetime('now')),
            notes TEXT
        )
        """
    )
    conn.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_ref_atlas_region_gene
        ON reference_expression(atlas_id, region_id, gene_symbol)
        """
    )
    conn.commit()


def next_atlas_id(conn: sqlite3.Connection) -> int:
    cur = conn.cursor()
    cur.execute("SELECT COALESCE(MAX(atlas_id), 0) FROM atlas_versions")
    return int(cur.fetchone()[0]) + 1


def upsert_atlas_version(conn: sqlite3.Connection, spec: dict) -> int:
    cur = conn.cursor()
    cur.execute(
        """
        SELECT atlas_id FROM atlas_versions
        WHERE atlas_name=? AND species=? AND build_version=?
        """,
        (spec["atlas_name"], spec["species"], spec["build_version"]),
    )
    row = cur.fetchone()
    if row:
        atlas_id = int(row[0])
        cur.execute("DELETE FROM reference_expression WHERE atlas_id=?", (atlas_id,))
        conn.commit()
        return atlas_id

    atlas_id = next_atlas_id(conn)
    cur.execute(
        """
        INSERT INTO atlas_versions (
            atlas_id, atlas_name, species, level, build_version,
            gene_id_type, normalization, notes
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            atlas_id,
            spec["atlas_name"],
            spec["species"],
            spec["level"],
            spec["build_version"],
            spec["gene_id_type"],
            spec["normalization"],
            spec["notes"],
        ),
    )
    conn.commit()
    return atlas_id


def import_matrix(conn: sqlite3.Connection, spec: dict) -> tuple[int, int]:
    df = pd.read_csv(spec["file"], sep="\t")
    expr_cols = [c for c in df.columns if c not in META_COLS]
    atlas_id = upsert_atlas_version(conn, spec)
    rows = []
    for region_name in expr_cols:
        region_id = clean_region_id(spec["prefix"], region_name)
        vals = pd.to_numeric(df[region_name], errors="coerce").fillna(0.0)
        for gene_symbol, gene_name, ensembl_id, avg_tpm in zip(
            df["gene_symbol"].astype(str),
            df["gene_symbol_raw"].astype(str),
            df["ensembl_id"].fillna("").astype(str),
            vals,
        ):
            rows.append(
                (
                    gene_symbol,
                    gene_name,
                    ensembl_id or None,
                    None,
                    region_id,
                    str(region_name),
                    float(avg_tpm),
                    0.0,
                    float(avg_tpm),
                    1,
                    infer_expression_class(float(avg_tpm)),
                    None,
                    atlas_id,
                )
            )

    conn.executemany(
        """
        INSERT INTO reference_expression (
            gene_symbol, gene_name, ensembl_id, ncbi_id,
            region_id, region_name, avg_tpm, std_tpm, median_tpm,
            sample_count, expression_class, cell_type_marker, atlas_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )
    conn.commit()
    return atlas_id, len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Import processed human atlas matrices into the main cfRNA SQLite database.")
    parser.add_argument("--db", default=str(DB_PATH))
    parser.add_argument("--build-signatures", action="store_true", help="Build signature sets for imported atlases after import.")
    args = parser.parse_args()

    conn = sqlite3.connect(args.db)
    try:
        ensure_schema(conn)
        imported = []
        for spec in ATLAS_SPECS:
            if not spec["file"].exists():
                print(f"skip missing file: {spec['file']}")
                continue
            atlas_id, n_rows = import_matrix(conn, spec)
            imported.append((atlas_id, spec["atlas_name"], n_rows))
            print(f"imported atlas_id={atlas_id} name={spec['atlas_name']} rows={n_rows}")
    finally:
        conn.close()

    if args.build_signatures:
        for atlas_id, atlas_name, _ in imported:
            sigset_id = build_signature_set(
                str(args.db),
                atlas_id=atlas_id,
                method="hybrid_specificity",
                topk_per_region=80,
                remove_housekeeping=True,
                remove_blood_background=True,
            )
            print(f"built signature set sigset_id={sigset_id} for atlas_id={atlas_id} name={atlas_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
