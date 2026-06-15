#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runinfo", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.runinfo)
    if "Run" not in df.columns:
        raise SystemExit("RunInfo.csv 里没找到 Run 列，请检查文件。")

    srr = df["Run"].dropna().astype(str).drop_duplicates().tolist()
    (outdir / "srr_list.txt").write_text("\n".join(srr) + "\n", encoding="utf-8")

    template = pd.DataFrame({
        "run_id": df["Run"].astype(str),
        "sample_name": df["SampleName"].astype(str) if "SampleName" in df.columns else df["Run"].astype(str),
        "library_name": df["LibraryName"].astype(str) if "LibraryName" in df.columns else "",
        "scientific_name": df["ScientificName"].astype(str) if "ScientificName" in df.columns else "",
        "spots": df["spots"] if "spots" in df.columns else "",
        "bases": df["bases"] if "bases" in df.columns else "",
        "platform": df["Platform"].astype(str) if "Platform" in df.columns else "",
        "monkey_id": "",
        "brain_region": "",
        "region_group": "",
        "include": 1,
        "notes": ""
    })
    template.to_csv(outdir / "sample_annotation_master.tsv", sep="\t", index=False)

    print(f"写出: {outdir / 'srr_list.txt'}")
    print(f"写出: {outdir / 'sample_annotation_master.tsv'}")

if __name__ == "__main__":
    main()
