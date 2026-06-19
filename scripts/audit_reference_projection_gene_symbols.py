#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_reference_projection_20260616"
DEFAULT_GENE_MAP = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.csv"
DEFAULT_LOCKED_GENES = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"


DATE_ISO_RE = re.compile(r"^\d{4}-\d{2}-\d{2}(?: 00:00:00)?$")
DATE_SLASH_RE = re.compile(r"^\d{1,2}/\d{1,2}/\d{2,4}$")
ENSEMBL_MACAQUE_RE = re.compile(r"^ENSMFAG\d{11}$")
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9.-]{1,30}$")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


def classify_symbol(value: object) -> list[str]:
    text = "" if pd.isna(value) else str(value).strip()
    flags: list[str] = []
    if not text:
        return ["empty"]
    if text.lower() in {"nan", "none", "null"}:
        flags.append("null_like")
    if DATE_ISO_RE.match(text) or DATE_SLASH_RE.match(text):
        flags.append("date_like")
    if ENSMBL := ENSEMBL_MACAQUE_RE.match(text):
        flags.append("ensembl_fallback")
    if len(text) > 40:
        flags.append("very_long")
    if any(ch.isspace() for ch in text):
        flags.append("contains_whitespace")
    if not GENE_SYMBOL_RE.match(text) and not ENSEMBL_MACAQUE_RE.match(text):
        flags.append("nonstandard_symbol_pattern")
    return flags or ["ok"]


def excel_date_candidate(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    match = re.match(r"^\d{4}-(\d{2})-(\d{2})(?: 00:00:00)?$", text)
    if not match:
        return ""
    month = int(match.group(1))
    day = int(match.group(2))
    month_prefix = {
        3: "MARCH",
        9: "SEPT",
        12: "DEC",
    }.get(month)
    if not month_prefix:
        return ""
    return f"{month_prefix}{day}"


def load_locked_genes(path: Path) -> set[str]:
    if not path.exists():
        return set()
    frame = pd.read_csv(path)
    if "gene_symbol" not in frame.columns:
        return set()
    return set(frame["gene_symbol"].dropna().astype(str).str.strip())


def load_external_overlap(outdir: Path) -> set[str]:
    matrix_path = outdir / "external_projected_vsd_GSE189919_matrix.tsv.gz"
    if not matrix_path.exists():
        return set()
    # Only the first column is needed, but pandas has no cheap compressed index-only
    # reader. The matrix is small enough for this audit.
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0, compression="infer")
    return set(matrix.index.astype(str))


def load_gse189919_raw_symbols() -> set[str]:
    path = ROOT / "data" / "external_validation" / "GSE189919" / "GSE189919_count.csv.gz"
    if not path.exists():
        return set()
    frame = pd.read_csv(path, usecols=["Geneid"], compression="infer")
    return set(frame["Geneid"].dropna().astype(str).str.strip())


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit gene symbols used by the Bo2023 reference projector.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--locked-genes", type=Path, default=DEFAULT_LOCKED_GENES)
    args = parser.parse_args()

    params_path = args.outdir / "projector_gene_parameters.csv"
    count_map_path = args.outdir / "count_gene_symbol_mapping_audit.csv"
    vsd_map_path = args.outdir / "vsd_gene_symbol_mapping_audit.csv"
    common_panel_path = args.outdir / "common_gene_panel.csv"

    params = pd.read_csv(params_path, dtype={"gene_symbol": str})
    count_map = pd.read_csv(count_map_path, dtype=str)
    vsd_map = pd.read_csv(vsd_map_path, dtype=str)
    common = pd.read_csv(common_panel_path, dtype={"gene_symbol": str})
    gene_map = pd.read_csv(args.gene_map, dtype=str)
    locked = load_locked_genes(args.locked_genes)
    external_projected_genes = load_external_overlap(args.outdir)
    gse189919_raw_symbols = load_gse189919_raw_symbols()

    audit = params[["gene_symbol"]].copy()
    audit["flags"] = audit["gene_symbol"].map(lambda x: ";".join(classify_symbol(x)))
    audit["excel_date_candidate_symbol"] = audit["gene_symbol"].map(excel_date_candidate)
    audit["is_suspicious"] = audit["flags"].ne("ok")
    audit["is_date_like"] = audit["flags"].str.contains("date_like", regex=False)
    audit["is_ensembl_fallback"] = audit["flags"].str.contains("ensembl_fallback", regex=False)
    audit["in_locked_network_model"] = audit["gene_symbol"].isin(locked)
    audit["in_gse189919_projected_matrix"] = audit["gene_symbol"].isin(external_projected_genes)
    audit["candidate_symbol_in_gse189919_raw_counts"] = audit["excel_date_candidate_symbol"].isin(gse189919_raw_symbols)
    audit["in_common_panel"] = audit["gene_symbol"].isin(set(common["gene_symbol"].astype(str)))

    count_sources = (
        count_map.groupby("gene_symbol")["gene_id"]
        .agg(lambda s: ";".join(sorted(set(map(str, s)))))
        .rename("count_gene_ids")
    )
    vsd_sources = (
        vsd_map.groupby("gene_symbol")["gene_id"]
        .agg(lambda s: ";".join(sorted(set(map(str, s)))))
        .rename("vsd_gene_ids")
    )
    audit = audit.merge(count_sources, left_on="gene_symbol", right_index=True, how="left")
    audit = audit.merge(vsd_sources, left_on="gene_symbol", right_index=True, how="left")
    audit["n_count_gene_ids"] = audit["count_gene_ids"].fillna("").map(lambda x: 0 if not x else len(str(x).split(";")))
    audit["n_vsd_gene_ids"] = audit["vsd_gene_ids"].fillna("").map(lambda x: 0 if not x else len(str(x).split(";")))
    audit["multi_count_ids_per_symbol"] = audit["n_count_gene_ids"] > 1
    audit["multi_vsd_ids_per_symbol"] = audit["n_vsd_gene_ids"] > 1

    suspicious = audit[audit["is_suspicious"] | audit["multi_count_ids_per_symbol"] | audit["multi_vsd_ids_per_symbol"]].copy()
    suspicious.to_csv(args.outdir / "gene_symbol_suspicious_audit.csv", index=False)
    audit.to_csv(args.outdir / "gene_symbol_full_audit.csv", index=False)

    date_like = audit[audit["is_date_like"]].copy()
    ensembl_fallback = audit[audit["is_ensembl_fallback"]].copy()
    multi_count = audit[audit["multi_count_ids_per_symbol"]].copy()
    multi_vsd = audit[audit["multi_vsd_ids_per_symbol"]].copy()

    raw_gene_map_flags = gene_map.assign(
        Gene_name_flags=gene_map["Gene.name"].map(lambda x: ";".join(classify_symbol(x)))
    )
    raw_suspicious_gene_map = raw_gene_map_flags[raw_gene_map_flags["Gene_name_flags"].ne("ok")].copy()
    raw_suspicious_gene_map.to_csv(args.outdir / "source_gene_map_suspicious_symbols.csv", index=False)

    date_like_count = int(audit["is_date_like"].sum())
    date_like_locked_count = int(date_like["in_locked_network_model"].sum())
    if date_like_count == 0 and date_like_locked_count == 0:
        recommendation = (
            "Date-like Excel-mangled symbols are fixed in this projector and locked gene panel. Remaining suspicious "
            "symbols are ENSMFAG fallback IDs or multi-ID aggregations; keep them documented, but they are not the "
            "date-conversion bug. External biological interpretation should still report the fallback-ID fraction."
        )
    else:
        recommendation = (
            "Rebuild the projector and locked network model with a cleaned gene map before using all-gene projected "
            "values for external interpretation. Date-like symbols should be restored to likely Excel-mangled symbols "
            "such as MARCH*, SEPT*, or DEC* when confirmed, or otherwise replaced by stable ENSMFAG IDs. The current "
            "internal validation is mostly preserved conceptually, but date-like locked genes should be corrected "
            "before publishing or expanding external claims."
        )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_projector_genes": int(len(audit)),
        "n_ok_symbols": int((audit["flags"] == "ok").sum()),
        "n_suspicious_symbols": int(audit["is_suspicious"].sum()),
        "n_date_like_symbols": date_like_count,
        "n_ensembl_fallback_symbols": int(audit["is_ensembl_fallback"].sum()),
        "n_multi_count_ids_per_symbol": int(audit["multi_count_ids_per_symbol"].sum()),
        "n_multi_vsd_ids_per_symbol": int(audit["multi_vsd_ids_per_symbol"].sum()),
        "n_suspicious_in_locked_network_model": int(suspicious["in_locked_network_model"].sum()),
        "n_date_like_in_locked_network_model": date_like_locked_count,
        "n_suspicious_in_gse189919_projected_matrix": int(suspicious["in_gse189919_projected_matrix"].sum()),
        "n_raw_gene_map_suspicious_names": int(len(raw_suspicious_gene_map)),
        "n_date_like_with_candidate_in_gse189919_raw_counts": int(
            date_like["excel_date_candidate_symbol"].isin(gse189919_raw_symbols).sum()
        ),
        "date_like_examples": date_like["gene_symbol"].head(20).tolist(),
        "ensembl_fallback_examples": ensembl_fallback["gene_symbol"].head(20).tolist(),
        "recommendation": recommendation,
        "outputs": {
            "full_audit": str(args.outdir / "gene_symbol_full_audit.csv"),
            "suspicious_audit": str(args.outdir / "gene_symbol_suspicious_audit.csv"),
            "source_gene_map_suspicious": str(args.outdir / "source_gene_map_suspicious_symbols.csv"),
        },
    }
    write_json(args.outdir / "gene_symbol_audit_summary.json", summary)

    report_lines = [
        "# Gene Symbol Audit",
        "",
        "## Summary",
        "",
        f"- Projector genes: {summary['n_projector_genes']}",
        f"- OK symbols: {summary['n_ok_symbols']}",
        f"- Suspicious symbols: {summary['n_suspicious_symbols']}",
        f"- Date-like symbols: {summary['n_date_like_symbols']}",
        f"- ENSMFAG fallback symbols: {summary['n_ensembl_fallback_symbols']}",
        f"- Multi count IDs per symbol: {summary['n_multi_count_ids_per_symbol']}",
        f"- Multi VSD IDs per symbol: {summary['n_multi_vsd_ids_per_symbol']}",
        f"- Suspicious symbols in locked 200 network model: {summary['n_suspicious_in_locked_network_model']}",
        f"- Date-like symbols in locked 200 network model: {summary['n_date_like_in_locked_network_model']}",
        f"- Suspicious symbols in GSE189919 projected matrix: {summary['n_suspicious_in_gse189919_projected_matrix']}",
        f"- Date-like candidate symbols present in GSE189919 raw counts: {summary['n_date_like_with_candidate_in_gse189919_raw_counts']}",
        "",
        "## Date-Like Examples",
        "",
        *[f"- `{x}`" for x in summary["date_like_examples"]],
        "",
        "## Recommendation",
        "",
        summary["recommendation"],
    ]
    (args.outdir / "gene_symbol_audit_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Wrote gene-symbol audit to {args.outdir}")
    print(
        "date_like={n_date_like_symbols} suspicious_locked={n_suspicious_in_locked_network_model}".format(
            **summary
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
