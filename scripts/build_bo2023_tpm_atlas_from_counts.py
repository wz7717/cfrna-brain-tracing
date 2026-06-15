#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "bo2023 data"
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"
DEFAULT_COUNTS = DEFAULT_DATA_DIR / "mfas5_819samples_28415genes_featurecounts_counts.txt"
DEFAULT_SAMPLE_INFO = DEFAULT_DATA_DIR / "Information of sequenced samples_update_full878_filter819.xlsx"
DEFAULT_GENE_MAP = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.csv"


def clean_text(value: object) -> str:
    return str(value or "").strip()


def expression_class(value: float) -> str:
    if value <= 0:
        return "silent"
    if value < 1:
        return "low"
    if value < 10:
        return "medium"
    return "high"


def inspect_counts_header(path: Path) -> dict:
    columns = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    return {
        "n_columns": len(columns),
        "first_columns": columns[:12],
        "has_geneid": "Geneid" in columns,
        "has_length": "Length" in columns,
        "has_featurecounts_annotation": all(c in columns for c in ["Geneid", "Chr", "Start", "End", "Strand", "Length"]),
    }


def read_gene_lengths(path: Path | None, id_col: str, length_col: str) -> pd.Series | None:
    if path is None:
        return None
    lengths = pd.read_csv(path)
    missing = [c for c in [id_col, length_col] if c not in lengths.columns]
    if missing:
        raise ValueError(f"gene length file is missing required columns: {missing}")
    s = lengths[[id_col, length_col]].copy()
    s[id_col] = s[id_col].astype(str).str.strip()
    s[length_col] = pd.to_numeric(s[length_col], errors="coerce")
    s = s.dropna(subset=[id_col, length_col])
    s = s[s[length_col] > 0]
    return s.drop_duplicates(subset=[id_col]).set_index(id_col)[length_col]


def read_counts_and_lengths(
    counts_path: Path,
    gene_lengths_path: Path | None,
    gene_length_id_col: str,
    gene_length_col: str,
) -> tuple[pd.DataFrame, pd.Series, dict]:
    header = inspect_counts_header(counts_path)
    if header["has_featurecounts_annotation"]:
        raw = pd.read_csv(counts_path, sep="\t")
        raw["Geneid"] = raw["Geneid"].astype(str).str.strip()
        lengths = pd.to_numeric(raw["Length"], errors="coerce")
        counts = raw.drop(columns=["Chr", "Start", "End", "Strand", "Length"], errors="ignore")
        counts = counts.set_index("Geneid")
    else:
        counts = pd.read_csv(counts_path, sep="\t", index_col=0)
        counts.index = counts.index.astype(str).str.strip()
        lengths = read_gene_lengths(gene_lengths_path, gene_length_id_col, gene_length_col)
        if lengths is None:
            raise ValueError(
                "The counts matrix does not contain a Length column. "
                "Strict TPM cannot be computed from this file alone. "
                "Provide a gene length table with --gene-lengths, or rebuild the original featureCounts output "
                "with Geneid/Chr/Start/End/Strand/Length columns."
            )

    counts = counts.loc[~counts.index.duplicated(keep="first")]
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lengths = pd.Series(lengths, index=counts.index if len(lengths) == len(counts.index) else lengths.index)
    lengths.index = lengths.index.astype(str).str.strip()
    lengths = pd.to_numeric(lengths, errors="coerce")
    lengths = lengths.reindex(counts.index)
    valid = lengths.notna() & (lengths > 0)
    counts = counts.loc[valid].copy()
    lengths = lengths.loc[valid].copy()
    audit = {
        **header,
        "n_genes_with_counts": int(len(valid)),
        "n_genes_with_valid_length": int(valid.sum()),
        "n_genes_dropped_missing_length": int((~valid).sum()),
        "n_samples": int(counts.shape[1]),
    }
    return counts, lengths, audit


def compute_tpm(counts: pd.DataFrame, lengths_bp: pd.Series) -> pd.DataFrame:
    length_kb = lengths_bp.astype(float) / 1000.0
    rpk = counts.div(length_kb, axis=0)
    denom = rpk.sum(axis=0).replace(0, pd.NA)
    return rpk.div(denom, axis=1).fillna(0.0) * 1_000_000.0


def read_sample_annotation(path: Path, sheet: str, region_col: str) -> pd.DataFrame:
    ann = pd.read_excel(path, sheet_name=sheet)
    if "No." not in ann.columns:
        raise ValueError("sample annotation must contain column 'No.'")
    if region_col not in ann.columns:
        raise ValueError(f"sample annotation must contain region column {region_col!r}")
    ann = ann.copy()
    ann["No."] = ann["No."].astype(str).str.strip()
    ann[region_col] = ann[region_col].astype(str).str.strip()
    ann = ann.dropna(subset=["No.", region_col]).drop_duplicates(subset=["No."])
    return ann


def read_region_names(path: Path) -> dict[str, str]:
    try:
        abbr = pd.read_excel(path, sheet_name="abbreviations")
    except Exception:
        return {}
    if "Abbreviation" not in abbr.columns:
        return {}
    full_col = "Full name " if "Full name " in abbr.columns else "Full name"
    if full_col not in abbr.columns:
        return {}
    return {
        clean_text(row["Abbreviation"]): clean_text(row[full_col])
        for _, row in abbr.iterrows()
        if clean_text(row["Abbreviation"])
    }


def read_gene_mapping(path: Path, gene_ids: pd.Index) -> tuple[pd.DataFrame, dict]:
    mapping = pd.read_csv(path)
    required = ["Gene.stable.ID", "Gene.name"]
    missing = [c for c in required if c not in mapping.columns]
    if missing:
        raise ValueError(f"gene map is missing required columns: {missing}")
    mapping = mapping[required].copy()
    mapping["ensembl_id"] = mapping["Gene.stable.ID"].astype(str).str.strip()
    mapping["gene_symbol"] = mapping["Gene.name"].astype(str).str.strip()
    mapping = mapping[mapping["ensembl_id"].isin(set(gene_ids.astype(str)))]
    missing_ids = sorted(set(gene_ids.astype(str)) - set(mapping["ensembl_id"]))
    audit = {
        "matrix_gene_ids": int(len(gene_ids)),
        "mapped_gene_ids": int(mapping["ensembl_id"].nunique()),
        "missing_gene_ids": int(len(missing_ids)),
        "missing_examples": missing_ids[:20],
        "id_multi_symbol": int((mapping.groupby("ensembl_id")["gene_symbol"].nunique() > 1).sum()),
        "symbol_multi_id": int((mapping.groupby("gene_symbol")["ensembl_id"].nunique() > 1).sum()),
        "mapping_rule": (
            "Use Gene.stable.ID as ENSMFAG ID and Gene.name as gene_symbol. "
            "Duplicate symbols are aggregated by region while preserving one representative ENSMFAG ID."
        ),
    }
    mapping = mapping.drop_duplicates(subset=["ensembl_id", "gene_symbol"])
    return mapping[["ensembl_id", "gene_symbol"]], audit


def build_region_long_table(
    tpm: pd.DataFrame,
    ann: pd.DataFrame,
    mapping: pd.DataFrame,
    region_col: str,
    region_names: dict[str, str],
) -> pd.DataFrame:
    sample_to_region = ann.set_index("No.")[region_col].to_dict()
    matched_samples = [c for c in tpm.columns if c in sample_to_region]
    if not matched_samples:
        raise ValueError("No count matrix samples match sample annotation column 'No.'")
    tpm = tpm[matched_samples].copy()
    region_series = pd.Series({s: sample_to_region[s] for s in matched_samples})

    frames = []
    for region_id, sample_ids in region_series.groupby(region_series).groups.items():
        sub = tpm[list(sample_ids)]
        stat = pd.DataFrame(
            {
                "ensembl_id": sub.index,
                "region_id": region_id,
                "region_name": region_names.get(region_id, region_id),
                "avg_tpm": sub.mean(axis=1).values,
                "median_tpm": sub.median(axis=1).values,
                "std_tpm": sub.std(axis=1, ddof=0).values,
                "sample_count": int(len(sample_ids)),
            }
        )
        frames.append(stat)
    long_df = pd.concat(frames, ignore_index=True)
    long_df = long_df.merge(mapping, on="ensembl_id", how="inner")
    long_df["gene_name"] = long_df["gene_symbol"]

    # Multiple ENSMFAG IDs may map to one symbol; aggregate after TPM calculation.
    agg = (
        long_df.groupby(["gene_symbol", "region_id"], as_index=False)
        .agg(
            gene_name=("gene_name", "first"),
            ensembl_id=("ensembl_id", "first"),
            region_name=("region_name", "first"),
            avg_tpm=("avg_tpm", "mean"),
            median_tpm=("median_tpm", "mean"),
            std_tpm=("std_tpm", "mean"),
            sample_count=("sample_count", "max"),
        )
    )
    agg["expression_class"] = agg["avg_tpm"].map(expression_class)
    agg["cell_type_marker"] = ""
    agg["ncbi_id"] = ""
    return agg[
        [
            "gene_symbol",
            "gene_name",
            "ensembl_id",
            "ncbi_id",
            "region_id",
            "region_name",
            "avg_tpm",
            "std_tpm",
            "median_tpm",
            "sample_count",
            "expression_class",
            "cell_type_marker",
        ]
    ]


def next_atlas_id(conn: sqlite3.Connection) -> int:
    row = conn.execute("SELECT COALESCE(MAX(atlas_id), 0) + 1 FROM atlas_versions").fetchone()
    return int(row[0])


def import_to_sqlite(
    db_path: Path,
    region_df: pd.DataFrame,
    atlas_name: str,
    build_version: str,
    summary: dict,
    replace_existing: bool,
) -> int:
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        existing = cur.execute("SELECT atlas_id FROM atlas_versions WHERE atlas_name=?", (atlas_name,)).fetchone()
        if existing and not replace_existing:
            raise ValueError(f"atlas_name {atlas_name!r} already exists. Use --replace-existing to overwrite it.")
        if existing:
            atlas_id = int(existing[0])
            cur.execute("DELETE FROM reference_expression WHERE atlas_id=?", (atlas_id,))
            cur.execute("DELETE FROM macaque_brain_atlas WHERE atlas_id=?", (atlas_id,))
            cur.execute("DELETE FROM atlas_versions WHERE atlas_id=?", (atlas_id,))
        else:
            atlas_id = next_atlas_id(conn)

        cur.execute(
            """
            INSERT INTO atlas_versions
            (atlas_id, atlas_name, species, level, build_version, gene_id_type, normalization, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                atlas_id,
                atlas_name,
                "Macaca fascicularis",
                "region",
                build_version,
                "ensembl_id_with_gene_symbol",
                "TPM_from_featureCounts",
                json.dumps(summary, ensure_ascii=False),
            ),
        )

        region_meta = region_df[["region_id", "region_name"]].drop_duplicates()
        for _, row in region_meta.iterrows():
            cur.execute(
                """
                INSERT INTO macaque_brain_atlas
                (region_id, region_name, region_acronym, atlas_version, atlas_id)
                VALUES (?, ?, ?, ?, ?)
                """,
                (row["region_id"], row["region_name"], row["region_id"], build_version, atlas_id),
            )

        insert_df = region_df.copy()
        insert_df["atlas_id"] = atlas_id
        insert_df.to_sql("reference_expression", conn, if_exists="append", index=False)
        conn.commit()
    return atlas_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a Bo2023 gene x region TPM atlas from raw featureCounts counts. "
            "Strict TPM requires gene lengths."
        )
    )
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--gene-lengths", type=Path, default=None, help="Optional CSV with gene ID and length in bp.")
    parser.add_argument("--gene-length-id-col", default="Geneid")
    parser.add_argument("--gene-length-col", default="Length")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--atlas-name", default="Bo2023_WangLab_TPM_region")
    parser.add_argument("--build-version", default="819samples_featureCounts_TPM_region")
    parser.add_argument("--output", type=Path, default=ROOT / "bo2023_tpm_region_reference.csv")
    parser.add_argument("--summary-output", type=Path, default=ROOT / "bo2023_tpm_region_build_summary.json")
    parser.add_argument("--import-db", action="store_true")
    parser.add_argument("--replace-existing", action="store_true")
    parser.add_argument("--inspect-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    header = inspect_counts_header(args.counts)
    if args.inspect_only:
        print(json.dumps(header, ensure_ascii=False, indent=2))
        return

    counts, lengths, count_audit = read_counts_and_lengths(
        args.counts,
        args.gene_lengths,
        args.gene_length_id_col,
        args.gene_length_col,
    )
    ann = read_sample_annotation(args.sample_info, args.sample_sheet, args.region_col)
    mapping, mapping_audit = read_gene_mapping(args.gene_map, counts.index)
    region_names = read_region_names(args.sample_info)
    tpm = compute_tpm(counts, lengths)
    region_df = build_region_long_table(tpm, ann, mapping, args.region_col, region_names)
    summary = {
        "counts": str(args.counts),
        "sample_info": str(args.sample_info),
        "gene_map": str(args.gene_map),
        "gene_lengths": str(args.gene_lengths) if args.gene_lengths else "Length column from featureCounts table",
        "atlas_name": args.atlas_name,
        "build_version": args.build_version,
        "normalization": "TPM_from_featureCounts",
        "count_matrix": count_audit,
        "mapping": mapping_audit,
        "n_regions": int(region_df["region_id"].nunique()),
        "n_gene_symbols_after_duplicate_symbol_aggregation": int(region_df["gene_symbol"].nunique()),
        "n_expression_rows": int(len(region_df)),
        "source_matrix_note": "TPM computed from raw counts and gene lengths; values are valid TPM estimates.",
    }
    region_df.to_csv(args.output, index=False)
    args.summary_output.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

    if args.import_db:
        atlas_id = import_to_sqlite(args.db, region_df, args.atlas_name, args.build_version, summary, args.replace_existing)
        print(f"Imported {args.atlas_name} as atlas_id={atlas_id}")


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
