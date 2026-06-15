#!/usr/bin/env python
from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_UNIQUE = (
    ROOT
    / "data"
    / "tcia_tcga_glioma_mri"
    / "tcia_data_retriever_manifests"
    / "tcga_73_complete_brats4_tcia_data_retriever_unique_series.csv"
)
DEFAULT_OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri" / "tcia_data_retriever_manifests"


def main() -> int:
    parser = argparse.ArgumentParser(description="Write a classic .tcia manifest for TCIA/NBIA Data Retriever.")
    parser.add_argument("--unique-series", type=Path, default=DEFAULT_UNIQUE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--include-annotation", default="false", choices=["true", "false"])
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.unique_series)
    uids = sorted(df["Series Instance UID"].dropna().astype(str).unique())
    manifest = args.outdir / "tcga_73_complete_brats4_unique_series.tcia"
    lines = [
        "downloadServerUrl=https://nbia.cancerimagingarchive.net/nbia-download/servlet/DownloadServlet",
        f"includeAnnotation={args.include_annotation}",
        "noOfrRetry=4",
        f"databasketId=manifest-tcga-73-complete-brats4-{int(time.time())}.tcia",
        "manifestVersion=3.0",
        "ListOfSeriesToDownload=",
        *uids,
        "",
    ]
    manifest.write_text("\n".join(lines), encoding="utf-8")
    template = args.outdir / "tcia_restricted_credentials_template.txt"
    if not template.exists():
        template.write_text("userName=YourUserName\npassWord=YourPassword\n", encoding="utf-8")
    print(f"wrote {manifest}")
    print(f"n_series={len(uids)}")
    print(f"credential_template={template}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
