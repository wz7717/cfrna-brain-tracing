from __future__ import annotations

from collections.abc import Mapping
import re

import pandas as pd


# Bo et al. report 110 post-QC anatomical regions.  Two individual assay rows
# carry a parent Network that disagrees with every other assay for the same
# anatomical region.  Canonicalize the parent label at analysis time without
# modifying the source workbook.
CANONICAL_REGION_NETWORKS: Mapping[str, str] = {
    "10m": "Orbitomedial Prefrontal Cortex (OMPFC)",
    "V2": "Occipital/Temporal",
}

# These are source-label corrections, not changes to the immutable Bo2023
# workbook.  Apply them only to fields that carry a canonical region identity;
# sample IDs, gene IDs, and provenance/audit fields must retain their source
# spelling.
CANONICAL_REGION_IDS: Mapping[str, str] = {
    "AI": "A1",
    "44563": "1-2",
}

CANONICAL_REGION_COLUMNS = frozenset(
    {
        "Region", "roi173", "regionLR", "regions", "region_ids", "region_labels",
        "candidate_regions", "region_id", "region_key",
        "bo2023_region_id", "display_region", "display_regions", "label", "class",
        "truth", "truth_region", "left_region", "right_region", "block_id", "group_id",
        "members", "group_members", "resolution_group", "resolution_group_members",
        "allowed_resolution_groups", "best_region_id", "region_top1", "region_top2",
        "region_top3", "region_top5", "region_top1_list", "region_top3_list",
        "exact_region_top3", "exact_top3", "region_group_top1", "region_group_top2",
        "region_group_top3", "region_group_top5",
        "group_top1", "group_top2", "group_top3", "pred_top1", "pred_top2", "pred_top3",
        "pred_top1_group_members", "pred_top1_resolution_group", "pred_group_top1",
        "pred_group_top2", "pred_group_top3", "true_resolution_group", "truth_group_members",
        "resolution_group_top1", "resolution_group_top2", "resolution_group_top3",
    }
)


def canonicalize_bo2023_region_text(value: object) -> object:
    """Canonicalize a region-bearing text value without touching other IDs."""

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    text = str(value)
    text = re.sub(r"(?<![A-Za-z0-9_])AI(?=_[LR]\b)", "A1", text)
    text = re.sub(r"(?<![A-Za-z0-9_])44563(?=_[LR]\b)", "1-2", text)
    text = re.sub(r"(?<![A-Za-z0-9_])AI(?![A-Za-z0-9_])", "A1", text)
    text = re.sub(r"(?<![A-Za-z0-9_])44563(?![A-Za-z0-9_])", "1-2", text)
    return text


def normalize_bo2023_region_labels(
    metadata: pd.DataFrame,
    *,
    columns: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Return metadata with source region IDs normalized at the canonical layer."""

    normalized = metadata.copy()
    selected = list(columns) if columns is not None else [
        column for column in normalized.columns if column in CANONICAL_REGION_COLUMNS
    ]
    for column in selected:
        if column in normalized.columns:
            normalized[column] = normalized[column].map(canonicalize_bo2023_region_text)
    return normalized


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
