#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path
import pandas as pd

def parse_gtf_gene_lengths(gtf_path: str) -> pd.DataFrame:
    gene_lengths = {}
    with open(gtf_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            start, end = int(parts[3]), int(parts[4])
            attrs = parts[8]
            gene_id = None
            gene_name = None
            for item in attrs.split(";"):
                item = item.strip()
                if item.startswith("gene_id "):
                    gene_id = item.split(" ", 1)[1].strip().strip('"')
                elif item.startswith("gene_name "):
                    gene_name = item.split(" ", 1)[1].strip().strip('"')
            if gene_id is None:
                continue
            gene_lengths.setdefault((gene_id, gene_name or gene_id), 0)
            gene_lengths[(gene_id, gene_name or gene_id)] += (end - start + 1)
    rows = [{"gene_id": k[0], "gene_name": k[1], "gene_length": v} for k, v in gene_lengths.items()]
    return pd.DataFrame(rows)

def load_featurecounts_table(path: Path):
    df = pd.read_csv(path, sep="\t", comment="#")
    count_col = df.columns[-1]
    meta = df.iloc[:, :6].copy()
    counts = df[["Geneid", count_col]].copy()
    counts.columns = ["gene_id", path.name.replace(".featureCounts.txt", "")]
    return meta, counts

def write_tsv_gz(df: pd.DataFrame, path: Path):
    with gzip.open(path, "wt", encoding="utf-8") as gz:
        df.to_csv(gz, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-dir", required=True)
    ap.add_argument("--sample-sheet", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    counts_dir = Path(args.counts_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted(counts_dir.glob("*.featureCounts.txt"))
    if not files:
        raise SystemExit("counts 目录里没找到 *.featureCounts.txt")

    metas = None
    merged = None
    for fp in files:
        meta, cnt = load_featurecounts_table(fp)
        if metas is None:
            metas = meta
            merged = cnt
        else:
            merged = merged.merge(cnt, on="gene_id", how="outer")

    merged = merged.fillna(0)
    gene_lengths = parse_gtf_gene_lengths(args.gtf)
    matrix = merged.merge(gene_lengths, on="gene_id", how="left")

    if "gene_name" not in matrix.columns:
        matrix["gene_name"] = matrix["gene_id"]
    if "gene_length" not in matrix.columns:
        raise SystemExit("没能从 GTF 解析出 gene_length。")

    sample_cols = [c for c in matrix.columns if c not in {"gene_id", "gene_name", "gene_length"}]
    counts_only = matrix[["gene_id", "gene_name"] + sample_cols].copy()
    write_tsv_gz(counts_only, outdir / "bo2023_gene_by_sample_counts.tsv.gz")

    # TPM
    tpm = matrix[["gene_id", "gene_name", "gene_length"]].copy()
    for col in sample_cols:
        rpk = matrix[col] / (matrix["gene_length"] / 1000.0)
        scale = rpk.sum() / 1e6
        tpm[col] = rpk / scale if scale > 0 else 0
    tpm_out = tpm.drop(columns=["gene_length"])
    write_tsv_gz(tpm_out, outdir / "bo2023_gene_by_sample_tpm.tsv.gz")

    # region mean TPM（需要样本注释）
    meta = pd.read_csv(args.sample_sheet, sep="\t")
    meta["run_id"] = meta["run_id"].astype(str)
    if "brain_region" in meta.columns and meta["brain_region"].fillna("").ne("").any():
        keep = meta[(meta.get("include", 1).astype(str) != "0") & meta["brain_region"].fillna("").ne("")]
        region_map = keep[["run_id", "brain_region"]].drop_duplicates()

        long = tpm_out.melt(id_vars=["gene_id", "gene_name"], var_name="run_id", value_name="TPM")
        long = long.merge(region_map, on="run_id", how="inner")
        region_mean = long.groupby(["gene_id", "gene_name", "brain_region"], as_index=False)["TPM"].mean()
        region_wide = region_mean.pivot(index=["gene_id", "gene_name"], columns="brain_region", values="TPM").reset_index()
        region_wide.columns.name = None
        write_tsv_gz(region_wide, outdir / "bo2023_gene_by_region_mean_tpm.tsv.gz")

    print("完成：")
    print(outdir / "bo2023_gene_by_sample_counts.tsv.gz")
    print(outdir / "bo2023_gene_by_sample_tpm.tsv.gz")
    if (outdir / "bo2023_gene_by_region_mean_tpm.tsv.gz").exists():
        print(outdir / "bo2023_gene_by_region_mean_tpm.tsv.gz")

if __name__ == "__main__":
    main()
