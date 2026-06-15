from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "tcia_tcga_glioma_mri" / "download_manifests" / "tcga_73_complete_brats4_selected_series_manifest.csv"
OUTDIR = ROOT / "data" / "tcia_tcga_glioma_mri" / "tcia_data_retriever_manifests"
OUTPUT = OUTDIR / "tcga_73_complete_brats4_retriever_series_uid_only.csv"


def main() -> None:
    df = pd.read_csv(INPUT, dtype=str).fillna("")
    if "series_instance_uid" not in df.columns:
        raise ValueError(f"Missing series_instance_uid in {INPUT}")
    series = (
        df[["series_instance_uid"]]
        .drop_duplicates()
        .rename(columns={"series_instance_uid": "SeriesInstanceUID"})
        .sort_values("SeriesInstanceUID")
    )
    OUTDIR.mkdir(parents=True, exist_ok=True)
    series.to_csv(OUTPUT, index=False)
    print(f"wrote {OUTPUT}")
    print(f"unique_series={len(series)}")


if __name__ == "__main__":
    main()
