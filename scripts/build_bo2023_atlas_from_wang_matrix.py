#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from signature_builder import build_signature_set  # noqa: E402


DEFAULT_DATA_DIR = ROOT / "bo2023 data"
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"
DEFAULT_MATRIX = DEFAULT_DATA_DIR / "mfas5_819samples_23605genes_vsd4_rmbatch.xls"
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


def read_expression_matrix(path: Path) -> pd.DataFrame:
    # The Wang lab .xls file is a tab-separated text matrix, not a binary Excel file.
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix = matrix.loc[matrix.sum(axis=1) != 0]
    return matrix


def read_sample_annotation(path: Path, sheet: str) -> pd.DataFrame:
    ann = pd.read_excel(path, sheet_name=sheet)
    if "No." not in ann.columns:
        raise ValueError("sample annotation must contain column 'No.'")
    ann["No."] = ann["No."].astype(str).str.strip()
    ann = ann.dropna(subset=["No."]).drop_duplicates(subset=["No."])
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


def read_gene_mapping(path: Path, matrix_gene_ids: Iterable[str]) -> tuple[pd.DataFrame, dict]:
    mapping = pd.read_csv(path)
    required = {"Gene.stable.ID", "Gene.name"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"gene mapping file is missing columns: {sorted(missing)}")

    mapping = mapping[["Gene.stable.ID", "Gene.name", "Gene.type_ensembl", "Gene.type_Fig2B"]].copy()
    mapping["Gene.stable.ID"] = mapping["Gene.stable.ID"].astype(str).str.strip()
    mapping["Gene.name"] = mapping["Gene.name"].astype(str).str.strip()
    mapping = mapping.dropna(subset=["Gene.stable.ID"]).drop_duplicates(subset=["Gene.stable.ID"], keep="first")

    matrix_ids = set(str(x).strip() for x in matrix_gene_ids)
    covered = mapping[mapping["Gene.stable.ID"].isin(matrix_ids)].copy()
    missing_ids = sorted(matrix_ids - set(covered["Gene.stable.ID"]))

    id_multi_symbol = int(covered.groupby("Gene.stable.ID")["Gene.name"].nunique().gt(1).sum())
    symbol_multi_id = int(
        covered[covered["Gene.name"].ne("")]
        .groupby("Gene.name")["Gene.stable.ID"]
        .nunique()
        .gt(1)
        .sum()
    )
    stats = {
        "matrix_gene_ids": len(matrix_ids),
        "mapped_gene_ids": int(covered["Gene.stable.ID"].nunique()),
        "missing_gene_ids": len(missing_ids),
        "id_multi_symbol": id_multi_symbol,
        "symbol_multi_id": symbol_multi_id,
        "missing_examples": missing_ids[:20],
        "mapping_rule": (
            "Use Gene.stable.ID as the original ENSMFAG ID and Gene.name as gene_symbol. "
            "ENSMFAG -> gene_symbol is unique in this source; gene_symbol -> ENSMFAG is not guaranteed unique. "
            "Duplicate symbols are aggregated by region while preserving one representative ENSMFAG ID and exporting the full mapping audit table."
        ),
    }
    return covered, stats


def validate_sample_alignment(matrix: pd.DataFrame, ann: pd.DataFrame) -> dict:
    matrix_samples = [str(c).strip() for c in matrix.columns]
    ann_samples = set(ann["No."].astype(str))
    missing_in_ann = sorted(set(matrix_samples) - ann_samples)
    missing_in_matrix = sorted(ann_samples - set(matrix_samples))
    return {
        "matrix_samples": len(matrix_samples),
        "annotation_samples": len(ann_samples),
        "matched_samples": len(set(matrix_samples) & ann_samples),
        "missing_in_annotation": missing_in_ann[:20],
        "missing_in_matrix": missing_in_matrix[:20],
        "n_missing_in_annotation": len(missing_in_ann),
        "n_missing_in_matrix": len(missing_in_matrix),
    }


def build_region_annotation(ann: pd.DataFrame, region_col: str, sample_info: Path) -> pd.DataFrame:
    if region_col not in ann.columns:
        raise ValueError(f"sample annotation does not contain region column: {region_col}")
    region_names = read_region_names(sample_info)
    rows = []
    for region_id, g in ann.groupby(region_col, dropna=True):
        rid = clean_text(region_id)
        if not rid:
            continue
        lobe = clean_text(g["Lobe"].mode().iloc[0]) if "Lobe" in g.columns and not g["Lobe"].mode().empty else ""
        saleem = clean_text(g["SaleemNetworks"].mode().iloc[0]) if "SaleemNetworks" in g.columns and not g["SaleemNetworks"].mode().empty else ""
        neocortex = clean_text(g["NeocortexRegion"].mode().iloc[0]) if "NeocortexRegion" in g.columns and not g["NeocortexRegion"].mode().empty else ""
        roi173 = clean_text(g["roi173"].mode().iloc[0]) if "roi173" in g.columns and not g["roi173"].mode().empty else ""
        rows.append(
            {
                "region_id": rid,
                "region_name": region_names.get(rid, rid),
                "region_acronym": rid,
                "parent_region_id": lobe or saleem,
                "lobe": lobe,
                "saleem_network": saleem,
                "neocortex_flag": neocortex,
                "roi173": roi173,
                "sample_count": int(len(g)),
            }
        )
    return pd.DataFrame(rows).sort_values("region_id").reset_index(drop=True)


def aggregate_expression_by_region(
    matrix: pd.DataFrame,
    ann: pd.DataFrame,
    gene_map: pd.DataFrame,
    region_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ann_use = ann[["No.", region_col]].dropna().copy()
    ann_use["No."] = ann_use["No."].astype(str).str.strip()
    ann_use[region_col] = ann_use[region_col].astype(str).str.strip()
    ann_use = ann_use[ann_use["No."].isin(matrix.columns)]

    gene_meta = gene_map.set_index("Gene.stable.ID")
    meta = pd.DataFrame(index=matrix.index)
    meta["ensembl_id"] = meta.index.astype(str)
    meta["gene_symbol"] = meta["ensembl_id"].map(gene_meta["Gene.name"]).fillna(meta["ensembl_id"])
    meta["gene_symbol"] = meta["gene_symbol"].replace("", pd.NA).fillna(meta["ensembl_id"])
    meta["gene_name"] = meta["gene_symbol"]

    chunks = []
    for region_id, sample_rows in ann_use.groupby(region_col):
        samples = [s for s in sample_rows["No."].tolist() if s in matrix.columns]
        if not samples:
            continue
        sub = matrix[samples]
        frame = pd.DataFrame(
            {
                "ensembl_id": meta["ensembl_id"].values,
                "gene_symbol": meta["gene_symbol"].values,
                "gene_name": meta["gene_name"].values,
                "region_id": clean_text(region_id),
                "avg_tpm": sub.mean(axis=1).values,
                "median_tpm": sub.median(axis=1).values,
                "std_tpm": sub.std(axis=1).fillna(0.0).values,
                "sample_count": len(samples),
            }
        )
        chunks.append(frame)

    long_df = pd.concat(chunks, ignore_index=True)

    # Because several ENSMFAG IDs can share one symbol, aggregate to one row per
    # gene_symbol x region for compatibility with gene-symbol based cfRNA inputs.
    audit = (
        long_df.groupby("gene_symbol")["ensembl_id"]
        .agg(lambda x: ";".join(sorted(set(map(str, x)))))
        .reset_index(name="ensembl_id_list")
    )
    audit["n_ensembl_ids"] = audit["ensembl_id_list"].str.count(";") + 1

    agg = (
        long_df.groupby(["gene_symbol", "gene_name", "region_id"], as_index=False)
        .agg(
            avg_tpm=("avg_tpm", "mean"),
            median_tpm=("median_tpm", "median"),
            std_tpm=("std_tpm", "mean"),
            sample_count=("sample_count", "max"),
            ensembl_id=("ensembl_id", "first"),
        )
    )
    agg["expression_class"] = agg["avg_tpm"].map(expression_class)
    return agg, audit


def next_atlas_id(conn: sqlite3.Connection) -> int:
    cur = conn.cursor()
    cur.execute("SELECT COALESCE(MAX(atlas_id), 0) + 1 FROM atlas_versions")
    return int(cur.fetchone()[0])


def existing_atlas_id(conn: sqlite3.Connection, atlas_name: str, build_version: str) -> int | None:
    cur = conn.cursor()
    cur.execute(
        "SELECT atlas_id FROM atlas_versions WHERE atlas_name=? AND build_version=?",
        (atlas_name, build_version),
    )
    row = cur.fetchone()
    return int(row[0]) if row else None


def insert_atlas(
    db_path: Path,
    atlas_name: str,
    build_version: str,
    normalization: str,
    region_df: pd.DataFrame,
    expr_df: pd.DataFrame,
    notes: dict,
    replace: bool,
) -> int:
    conn = sqlite3.connect(str(db_path), timeout=30)
    try:
        cur = conn.cursor()
        atlas_id = existing_atlas_id(conn, atlas_name, build_version)
        if atlas_id is not None and not replace:
            raise ValueError(f"atlas already exists: atlas_id={atlas_id}. Use --replace to rebuild it.")
        if atlas_id is not None and replace:
            cur.execute("DELETE FROM reference_expression WHERE atlas_id=?", (atlas_id,))
            cur.execute("DELETE FROM macaque_brain_atlas WHERE atlas_id=?", (atlas_id,))
            cur.execute("DELETE FROM signature_genes WHERE sigset_id IN (SELECT sigset_id FROM signature_sets WHERE atlas_id=?)", (atlas_id,))
            cur.execute("DELETE FROM signature_sets WHERE atlas_id=?", (atlas_id,))
            cur.execute("DELETE FROM atlas_versions WHERE atlas_id=?", (atlas_id,))
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
                normalization,
                json.dumps(notes, ensure_ascii=False),
            ),
        )

        atlas_version_text = f"{atlas_name} | {build_version}"
        atlas_rows = []
        for r in region_df.itertuples(index=False):
            coords = json.dumps(
                {
                    "lobe": r.lobe,
                    "saleem_network": r.saleem_network,
                    "neocortex_flag": r.neocortex_flag,
                    "roi173": r.roi173,
                    "source_sample_count": int(r.sample_count),
                },
                ensure_ascii=False,
            )
            atlas_rows.append(
                (
                    str(r.region_id),
                    str(r.region_name),
                    str(r.region_acronym),
                    str(r.parent_region_id),
                    None,
                    None,
                    atlas_version_text,
                    coords,
                    atlas_id,
                )
            )
        cur.executemany(
            """
            INSERT INTO macaque_brain_atlas
            (region_id, region_name, region_acronym, parent_region_id, hemi, layer, atlas_version, coordinates, atlas_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            atlas_rows,
        )

        region_name_map = dict(zip(region_df["region_id"].astype(str), region_df["region_name"].astype(str)))
        rows = []
        for r in expr_df.itertuples(index=False):
            rows.append(
                (
                    str(r.gene_symbol),
                    str(r.gene_name),
                    str(r.ensembl_id),
                    None,
                    str(r.region_id),
                    region_name_map.get(str(r.region_id), str(r.region_id)),
                    float(r.avg_tpm),
                    float(r.std_tpm),
                    float(r.median_tpm),
                    int(r.sample_count),
                    str(r.expression_class),
                    None,
                    atlas_id,
                )
            )
        cur.executemany(
            """
            INSERT INTO reference_expression
            (gene_symbol, gene_name, ensembl_id, ncbi_id, region_id, region_name,
             avg_tpm, std_tpm, median_tpm, sample_count, expression_class,
             cell_type_marker, atlas_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            rows,
        )
        conn.commit()
        return atlas_id
    finally:
        conn.close()


def write_outputs(outdir: Path, region_df: pd.DataFrame, expr_df: pd.DataFrame, audit_df: pd.DataFrame, summary: dict) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    region_df.to_csv(outdir / "bo2023_region_annotation.csv", index=False, encoding="utf-8-sig")
    audit_df.to_csv(outdir / "bo2023_ensmfag_to_symbol_audit.csv", index=False, encoding="utf-8-sig")
    expr_df.to_csv(outdir / "bo2023_gene_by_region_long.csv", index=False, encoding="utf-8-sig")
    (outdir / "bo2023_build_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a Bo2023 Wang-lab macaque brain atlas from gene-by-sample matrices.")
    parser.add_argument("--db", default=str(DEFAULT_DB))
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX), help="VSD + batch-removed gene x sample matrix.")
    parser.add_argument("--counts", default=str(DEFAULT_COUNTS), help="Raw counts matrix, recorded in metadata only.")
    parser.add_argument("--sample-info", default=str(DEFAULT_SAMPLE_INFO))
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--gene-map", default=str(DEFAULT_GENE_MAP))
    parser.add_argument("--region-col", default="Region", help="Annotation column used to aggregate samples, e.g. Region, roi173, Lobe.")
    parser.add_argument("--atlas-name", default="Bo2023_WangLab_VSD_region")
    parser.add_argument("--build-version", default="819samples_vsd4_rmbatch_region")
    parser.add_argument("--normalization", default="VSD_batch_removed")
    parser.add_argument("--outdir", default=str(ROOT / "bo2023_atlas_build_output"))
    parser.add_argument("--dry-run", action="store_true", help="Build and export files without writing SQLite.")
    parser.add_argument("--replace", action="store_true", help="Replace an existing atlas with the same name and build_version.")
    parser.add_argument("--build-signature", action="store_true", help="Build a signature set after SQLite import.")
    parser.add_argument("--topk-per-region", type=int, default=80)
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    sample_info = Path(args.sample_info)
    gene_map_path = Path(args.gene_map)
    outdir = Path(args.outdir)

    matrix = read_expression_matrix(matrix_path)
    ann = read_sample_annotation(sample_info, args.sample_sheet)
    alignment = validate_sample_alignment(matrix, ann)
    if alignment["n_missing_in_annotation"] or alignment["n_missing_in_matrix"]:
        raise ValueError(f"sample IDs are not fully aligned: {alignment}")

    gene_map, mapping_stats = read_gene_mapping(gene_map_path, matrix.index)
    if mapping_stats["missing_gene_ids"]:
        raise ValueError(f"gene mapping is incomplete: {mapping_stats}")

    region_df = build_region_annotation(ann, args.region_col, sample_info)
    expr_df, audit_df = aggregate_expression_by_region(matrix, ann, gene_map, args.region_col)

    summary = {
        "matrix": str(matrix_path),
        "counts": str(Path(args.counts)),
        "sample_info": str(sample_info),
        "gene_map": str(gene_map_path),
        "region_col": args.region_col,
        "atlas_name": args.atlas_name,
        "build_version": args.build_version,
        "normalization": args.normalization,
        "alignment": alignment,
        "mapping": mapping_stats,
        "n_regions": int(region_df["region_id"].nunique()),
        "n_gene_symbols_after_duplicate_symbol_aggregation": int(expr_df["gene_symbol"].nunique()),
        "n_expression_rows": int(len(expr_df)),
        "n_symbols_with_multiple_ensembl_ids": int((audit_df["n_ensembl_ids"] > 1).sum()),
    }
    write_outputs(outdir, region_df, expr_df, audit_df, summary)

    print(json.dumps(summary, indent=2, ensure_ascii=False))

    if args.dry_run:
        print(f"dry-run complete; files written to {outdir}")
        return 0

    notes = {
        **summary,
        "source_matrix_note": "Input matrix is VSD + batch removed; values are stored in avg_tpm/median_tpm fields for schema compatibility but are not TPM.",
    }
    atlas_id = insert_atlas(
        Path(args.db),
        args.atlas_name,
        args.build_version,
        args.normalization,
        region_df,
        expr_df,
        notes,
        replace=args.replace,
    )
    print(f"imported atlas_id={atlas_id}")

    if args.build_signature:
        sigset_id = build_signature_set(
            str(args.db),
            atlas_id=atlas_id,
            method="hybrid_specificity",
            topk_per_region=args.topk_per_region,
            remove_housekeeping=True,
            remove_blood_background=True,
            params={"source": "Bo2023 Wang-lab VSD region atlas"},
        )
        print(f"built signature set sigset_id={sigset_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
