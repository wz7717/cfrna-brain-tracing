from __future__ import annotations

from collections.abc import Mapping

import pandas as pd


# Bo et al. report 110 post-QC anatomical regions.  Two individual assay rows
# carry a parent Network that disagrees with every other assay for the same
# anatomical region.  Canonicalize the parent label at analysis time without
# modifying the source workbook.
CANONICAL_REGION_NETWORKS: Mapping[str, str] = {
    "10m": "Orbitomedial Prefrontal Cortex (OMPFC)",
    "V2": "Occipital/Temporal",
}


def normalize_bo2023_network_labels(
    metadata: pd.DataFrame,
    *,
    region_col: str = "region_id",
    network_col: str = "network_id",
) -> pd.DataFrame:
    """Return Bo2023 metadata with one canonical parent Network per region.

    The operation is deliberately performed on a copy so the source metadata
    remains an immutable audit input.
    """

    missing = {region_col, network_col} - set(metadata.columns)
    if missing:
        raise ValueError(f"Bo2023 metadata missing canonicalization columns: {sorted(missing)}")

    normalized = metadata.copy()
    regions = normalized[region_col].astype(str).str.strip()
    normalized[region_col] = regions
    normalized[network_col] = normalized[network_col].astype(str).str.strip()
    for region_id, canonical_network in CANONICAL_REGION_NETWORKS.items():
        normalized.loc[regions.eq(region_id), network_col] = canonical_network
    return normalized


def assert_unique_region_network_mapping(
    metadata: pd.DataFrame,
    *,
    region_col: str = "region_id",
    network_col: str = "network_id",
) -> None:
    """Fail when one anatomical region remains assigned to multiple Networks."""

    counts = metadata.groupby(region_col, dropna=False)[network_col].nunique(dropna=False)
    ambiguous = counts[counts > 1]
    if len(ambiguous):
        raise ValueError(
            "Bo2023 regions retain multiple parent Networks after canonicalization: "
            + ", ".join(map(str, ambiguous.index.tolist()))
        )
