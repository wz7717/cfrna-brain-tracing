#!/usr/bin/env python
"""Build first-version processed reference matrices for cfRNA brain tracing."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "reference_atlases"
OUT = ROOT / "data" / "processed_reference" / "v1"


def strip_ensembl_version(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    return re.sub(r"\.\d+$", "", text)


def clean_symbol(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return text.upper()


def numeric_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in {"gene_symbol", "gene_symbol_raw", "ensembl_id", "source_atlas", "doi", "recommended_use"}]


def collapse_gene_matrix(df: pd.DataFrame, method: str = "median") -> tuple[pd.DataFrame, int]:
    df = df[df["gene_symbol"].astype(str).str.len() > 0].copy()
    before = len(df)
    nums = numeric_cols(df)
    for col in nums:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    agg = {col: method for col in nums}
    agg["gene_symbol_raw"] = "first"
    agg["ensembl_id"] = "first"
    collapsed = df.groupby("gene_symbol", as_index=False).agg(agg)
    collapsed = collapsed[["gene_symbol", "gene_symbol_raw", "ensembl_id", *nums]]
    return collapsed, before - len(collapsed)


def write_matrix(df: pd.DataFrame, path: Path, source_atlas: str, doi: str, recommended_use: str) -> pd.DataFrame:
    df = df.copy()
    for col in ["source_atlas", "doi", "recommended_use"]:
        if col in df.columns:
            df = df.drop(columns=[col])
    df.insert(3, "source_atlas", source_atlas)
    df.insert(4, "doi", doi)
    df.insert(5, "recommended_use", recommended_use)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return df


def build_gtex(report: list[str]) -> pd.DataFrame:
    path = REF / "04_peripheral_background" / "gtex_v8" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
    df = pd.read_csv(path, sep="\t", skiprows=2)
    expr_cols = [c for c in df.columns if c not in {"Name", "Description"}]
    out = pd.DataFrame(
        {
            "gene_symbol": df["Description"].map(clean_symbol),
            "gene_symbol_raw": df["Description"],
            "ensembl_id": df["Name"].map(strip_ensembl_version),
        }
    )
    out[expr_cols] = np.log1p(df[expr_cols].apply(pd.to_numeric, errors="coerce"))
    out, dup = collapse_gene_matrix(out, "median")
    report.append(f"- GTEx median TPM: {out.shape[0]} genes x {len(expr_cols)} tissues; duplicate symbols collapsed={dup}.")
    return write_matrix(
        out,
        OUT / "peripheral_background_gtex_median_tpm.tsv",
        "GTEx_v8_bulk_RNAseq",
        "10.1126/science.aaz1776",
        "peripheral_background",
    )


def build_allen_hba(report: list[str]) -> pd.DataFrame:
    base = REF / "02_normal_human_brain_region" / "allen_human_brain_atlas"
    matrices = []
    for donor_dir in sorted(base.glob("H0351.*_rnaseq")):
        tpm = donor_dir / "RNAseqTPM.csv"
        meta = donor_dir / "SampleAnnot.csv"
        if not tpm.exists() or not meta.exists():
            continue
        sample_meta = pd.read_csv(meta)
        sample_names = sample_meta["RNAseq_sample_name"].astype(str).tolist()
        region_names = (
            sample_meta["main_structure"].astype(str).str.strip()
            + "|"
            + sample_meta["sub_structure"].astype(str).str.strip()
            + "|"
            + donor_dir.name.split("_")[0]
        ).tolist()
        df = pd.read_csv(tpm, header=None)
        if df.shape[1] - 1 != len(sample_names):
            report.append(f"- Allen HBA warning: {donor_dir.name} sample metadata length mismatch; using generic sample names.")
            sample_names = [f"sample_{i}" for i in range(df.shape[1] - 1)]
            region_names = sample_names
        sample_names = [f"{name}_{idx}" if sample_names.count(name) > 1 else name for idx, name in enumerate(sample_names)]
        sample_to_region = dict(zip(sample_names, region_names))
        df.columns = ["gene_symbol_raw", *sample_names]
        out = pd.DataFrame({"gene_symbol": df["gene_symbol_raw"].map(clean_symbol), "gene_symbol_raw": df["gene_symbol_raw"], "ensembl_id": ""})
        value_cols = sample_names
        out[value_cols] = np.log1p(df[value_cols].apply(pd.to_numeric, errors="coerce"))
        group_cols = {}
        for region in sorted(set(region_names)):
            cols = [c for c in value_cols if sample_to_region[c] == region]
            group_cols[region] = out[cols].mean(axis=1)
        grouped = pd.concat([out[["gene_symbol", "gene_symbol_raw", "ensembl_id"]], pd.DataFrame(group_cols)], axis=1)
        matrices.append(grouped)
    if not matrices:
        report.append("- Allen HBA RNA-seq: no usable RNAseqTPM.csv files found.")
        return pd.DataFrame()
    merged = matrices[0]
    for mat in matrices[1:]:
        merged = merged.merge(mat, on=["gene_symbol", "gene_symbol_raw", "ensembl_id"], how="outer")
    expr = numeric_cols(merged)
    region_base: dict[str, list[str]] = {}
    for col in expr:
        key = re.sub(r"\|H0351\.\d+$", "", col)
        region_base.setdefault(key, []).append(col)
    final = merged[["gene_symbol", "gene_symbol_raw", "ensembl_id"]].copy()
    for region, cols in region_base.items():
        final[region] = merged[cols].mean(axis=1)
    final, dup = collapse_gene_matrix(final, "median")
    report.append(f"- Allen HBA RNA-seq TPM: {final.shape[0]} genes x {len(numeric_cols(final))} region groups; duplicate symbols collapsed={dup}.")
    return write_matrix(
        final,
        OUT / "brain_region_reference_allen_hba.tsv",
        "Allen_Human_Brain_Atlas_Hawrylycz2012",
        "10.1038/nature11405",
        "brain_region_reference",
    )


def pivot_hpa_long(path: Path, group_col: str, value_col: str, report: list[str], label: str) -> pd.DataFrame:
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=500_000):
        needed = ["Gene", "Gene name", group_col, value_col]
        missing = [c for c in needed if c not in chunk.columns]
        if missing:
            raise ValueError(f"{path} missing columns: {missing}")
        sub = chunk[needed].copy()
        sub["gene_symbol"] = sub["Gene name"].map(clean_symbol)
        sub["gene_symbol_raw"] = sub["Gene name"]
        sub["ensembl_id"] = sub["Gene"].map(strip_ensembl_version)
        sub[value_col] = pd.to_numeric(sub[value_col], errors="coerce")
        chunks.append(sub.groupby(["gene_symbol", "gene_symbol_raw", "ensembl_id", group_col], as_index=False)[value_col].mean())
    long = pd.concat(chunks, ignore_index=True)
    long = long.groupby(["gene_symbol", "gene_symbol_raw", "ensembl_id", group_col], as_index=False)[value_col].mean()
    pivot = long.pivot_table(index=["gene_symbol", "gene_symbol_raw", "ensembl_id"], columns=group_col, values=value_col, aggfunc="mean").reset_index()
    pivot.columns = [str(c) for c in pivot.columns]
    for col in numeric_cols(pivot):
        pivot[col] = np.log1p(pd.to_numeric(pivot[col], errors="coerce"))
    pivot, dup = collapse_gene_matrix(pivot, "median")
    report.append(f"- {label}: {pivot.shape[0]} genes x {len(numeric_cols(pivot))} groups; duplicate symbols collapsed={dup}.")
    return pivot


def build_hpa_region(report: list[str]) -> pd.DataFrame:
    path = REF / "02_normal_human_brain_region" / "human_protein_atlas_brain" / "rna_brain_region_hpa.tsv" / "rna_brain_region_hpa.tsv"
    df = pivot_hpa_long(path, "Brain region", "nTPM", report, "HPA brain region nTPM")
    return write_matrix(df, OUT / "brain_region_reference_hpa.tsv", "Human_Protein_Atlas_Brain_Sjostedt2020", "10.1126/science.aay5947", "brain_region_reference")


def build_hpa_celltype(report: list[str]) -> pd.DataFrame:
    path = REF / "02_normal_human_brain_region" / "human_protein_atlas_brain" / "rna_single_nuclei_cluster_type.tsv" / "rna_single_nuclei_cluster_type.tsv"
    df = pivot_hpa_long(path, "Cluster type", "nCPM", report, "HPA single-nuclei brain cluster-type nCPM")
    return write_matrix(df, OUT / "celltype_reference_hpa_single_nuclei_brain.tsv", "Human_Protein_Atlas_Brain_Sjostedt2020", "10.1126/science.aay5947", "celltype_reference")


def build_garza(report: list[str]) -> pd.DataFrame:
    base = REF / "03_brain_injury_state" / "garza2023_tbi_gse209552"
    matrices = []
    for path in sorted(base.glob("GSE209552_*gene_count_matrix_2.csv")):
        label = "TBI" if "TBI_gene" in path.name else "hGPC"
        df = pd.read_csv(path, sep="\t")
        sample_cols = [c for c in df.columns if c not in {"Geneid", "Chr", "Start", "End", "Strand", "Length"}]
        counts = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        lib = counts.sum(axis=0).replace(0, np.nan)
        cpm = counts.div(lib, axis=1) * 1_000_000
        out = pd.DataFrame({"gene_symbol": df["Geneid"].map(clean_symbol), "gene_symbol_raw": df["Geneid"], "ensembl_id": ""})
        out[[f"{label}|{c}" for c in sample_cols]] = np.log1p(cpm)
        matrices.append(out)
    if not matrices:
        report.append("- Garza TBI: no count matrices found.")
        return pd.DataFrame()
    merged = matrices[0]
    for mat in matrices[1:]:
        merged = merged.merge(mat, on=["gene_symbol", "gene_symbol_raw", "ensembl_id"], how="outer")
    merged, dup = collapse_gene_matrix(merged, "median")
    report.append(f"- Garza 2023 TBI/hGPC raw counts converted to log1p CPM: {merged.shape[0]} genes x {len(numeric_cols(merged))} samples; duplicate symbols collapsed={dup}.")
    return write_matrix(merged, OUT / "injury_state_reference_garza2023_tbi.tsv", "Garza2023_human_TBI_snRNAseq", "10.1016/j.celrep.2023.113395", "injury_state_reference")


def build_allen_aging(report: list[str]) -> pd.DataFrame:
    base = REF / "03_brain_injury_state" / "allen_aging_dementia_tbi_gse104687"
    donors = pd.read_csv(base / "query.csv")
    files = pd.read_csv(base / "tbi_data_files.csv")
    merged = files.merge(donors, on="donor_id", how="left", suffixes=("", "_donor"))
    merged.insert(0, "source_atlas", "Allen_Aging_Dementia_TBI_GSE104687")
    merged.insert(1, "doi", "GSE104687")
    merged.insert(2, "recommended_use", "injury_state_reference")
    path = OUT / "injury_state_reference_allen_aging_tbi.tsv"
    merged.to_csv(path, sep="\t", index=False)
    report.append(f"- Allen Aging/Dementia/TBI: expression links and sample/pathology metadata retained for {len(merged)} rows; gene-level FPKM files were not downloaded in this no-redownload stage.")
    return merged


def build_gene_universe(mats: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for source, df in mats.items():
        if "gene_symbol" not in df.columns:
            continue
        sub = df[["gene_symbol", "gene_symbol_raw", "ensembl_id"]].drop_duplicates().copy()
        sub["source_reference"] = source
        rows.append(sub)
    uni = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["gene_symbol", "gene_symbol_raw", "ensembl_id", "source_reference"])
    flags = uni.groupby("gene_symbol")["source_reference"].apply(lambda x: ";".join(sorted(set(map(str, x))))).reset_index(name="sources")
    raw = uni.groupby("gene_symbol", as_index=False).agg(gene_symbol_raw=("gene_symbol_raw", "first"), ensembl_id=("ensembl_id", "first"))
    out = raw.merge(flags, on="gene_symbol", how="left")
    out.to_csv(OUT / "reference_gene_universe.tsv", sep="\t", index=False)
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report = ["# Reference Build Report", ""]
    mats: dict[str, pd.DataFrame] = {}
    mats["gtex"] = build_gtex(report)
    mats["allen_hba"] = build_allen_hba(report)
    mats["hpa_region"] = build_hpa_region(report)
    mats["garza2023_tbi"] = build_garza(report)
    build_allen_aging(report)
    mats["hpa_celltype"] = build_hpa_celltype(report)
    universe = build_gene_universe(mats)
    report.append(f"- Reference gene universe: {len(universe)} unique upper-case gene symbols.")
    report.extend(
        [
            "",
            "## Standardization",
            "",
            "- Expression matrices are gene x sample_or_group.",
            "- Gene symbols are upper-case in `gene_symbol`; original names are retained in `gene_symbol_raw`.",
            "- Ensembl versions are stripped where Ensembl IDs exist.",
            "- TPM/nTPM matrices are transformed with log1p.",
            "- Raw count matrices are converted to CPM and then transformed with log1p.",
            "- Duplicate gene symbols are collapsed by median.",
            "- FASTQ/BAM/bigWig files are not processed.",
            "",
            "## Limitations",
            "",
            "- Allen Aging/Dementia/TBI has downloaded metadata and FPKM links only, not the per-sample expression files.",
            "- Hodge 2019, Siletti ABC metadata, and Brain Cell Atlas remain manual or --large for first-version cell-level integration.",
        ]
    )
    (OUT / "reference_build_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote processed reference v1 to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
