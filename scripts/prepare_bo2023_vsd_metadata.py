from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_metadata import normalize_bo2023_region_labels  # noqa: E402


SHEET = "mfas5_819samples_phenSet4"
REQUIRED_COLUMNS = [
    "No.",
    "MonkeyID",
    "Sample",
    "Side",
    "Batch",
    "batch2",
    "Region",
    "Lobe",
    "SaleemNetworks",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export the Bo2023 819-sample phenotype sheet for R."
    )
    parser.add_argument("workbook", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    metadata = pd.read_excel(args.workbook, sheet_name=SHEET)
    missing = [column for column in REQUIRED_COLUMNS if column not in metadata.columns]
    if missing:
        raise ValueError(f"Missing required metadata columns: {missing}")

    metadata = metadata.loc[:, REQUIRED_COLUMNS].copy()
    metadata = normalize_bo2023_region_labels(metadata, columns=["Region"])
    metadata["No."] = metadata["No."].astype(str)
    if metadata["No."].duplicated().any():
        duplicates = metadata.loc[metadata["No."].duplicated(), "No."].tolist()
        raise ValueError(f"Duplicate sample IDs in metadata: {duplicates[:10]}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    metadata.set_index("No.").to_csv(args.output_csv)
    print(f"Wrote {len(metadata)} samples to {args.output_csv}")


if __name__ == "__main__":
    main()
