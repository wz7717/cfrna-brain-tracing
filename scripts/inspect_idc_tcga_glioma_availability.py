#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "idc_index_py312"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

from idc_index import IDCClient  # noqa: E402


def main() -> int:
    client = IDCClient()
    patients = ["TCGA-02-0003", "TCGA-06-0130", "TCGA-CS-4941", "TCGA-DU-6400"]
    print("IDC version", client.get_idc_version())
    for patient in patients:
        df = client.sql_query(
            f"""
            SELECT collection_id, PatientID, SeriesInstanceUID, Modality, SeriesDescription
            FROM index
            WHERE PatientID = '{patient}'
            LIMIT 20
            """
        )
        print("\nPATIENT", patient, "rows", len(df))
        print(df.to_string(index=False)[:3000])
    collections = client.sql_query(
        """
        SELECT DISTINCT collection_id
        FROM index
        WHERE lower(collection_id) LIKE '%tcga%'
           OR lower(collection_id) LIKE '%gbm%'
           OR lower(collection_id) LIKE '%lgg%'
           OR lower(collection_id) LIKE '%brats%'
        ORDER BY collection_id
        """
    )
    print("\nCOLLECTIONS")
    print(collections.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
