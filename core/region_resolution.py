from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any
import math


DEFAULT_BO2023_REGION_RESOLUTION_MODEL = (
    Path(__file__).resolve().parents[1] / "data" / "models" / "bo2023_region_resolution_groups.json"
)


def load_region_resolution_model(
    model_path: Path = DEFAULT_BO2023_REGION_RESOLUTION_MODEL,
) -> dict[str, Any]:
    if not model_path.exists():
        return {}
    return json.loads(model_path.read_text(encoding="utf-8"))


def annotate_region_candidates(
    region_output: dict[str, Any],
    network_output: dict[str, Any],
    model_path: Path = DEFAULT_BO2023_REGION_RESOLUTION_MODEL,
) -> dict[str, Any]:
    """Add resolution warnings to Region candidates without modifying their order or scores."""
    annotated = deepcopy(region_output)
    rows = annotated.get("results", [])
    network_rows = network_output.get("results", [])
    model = load_region_resolution_model(model_path)
    if not rows or not network_rows or not model:
        return annotated

    primary_network = str(network_rows[0].get("network_id", ""))
    entries = model.get("entries", {})
    for row in rows:
        region_id = str(row.get("region_id", ""))
        entry = entries.get(f"{primary_network}||{region_id}")
        if entry is None:
            entry = {
                "resolution_tier": "low_resolution",
                "resolution_group": region_id,
                "group_members": [region_id],
                "resolution_reasons": ["outside_primary_network_resolution_map"],
            }
        row["resolution_tier"] = str(entry.get("resolution_tier", "low_resolution"))
        row["resolution_group"] = str(entry.get("resolution_group", region_id))
        row["resolution_group_members"] = " | ".join(map(str, entry.get("group_members", [region_id])))
        row["resolution_reasons"] = ";".join(map(str, entry.get("resolution_reasons", [])))
        row["group_plausibility_tier"] = str(entry.get("group_plausibility_tier", "not_calibrated"))
        row["group_calibration_flags"] = ";".join(map(str, entry.get("group_calibration_flags", [])))
        row["manual_review_recommended"] = row["resolution_tier"] == "low_resolution"

    group_scores: dict[str, dict[str, Any]] = {}
    for row in rows:
        group = str(row.get("resolution_group", row.get("region_id", "")))
        score = row.get("score")
        try:
            numeric_score = float(score)
        except (TypeError, ValueError):
            numeric_score = float("-inf")
        item = group_scores.setdefault(
            group,
            {
                "resolution_group": group,
                "best_region_id": str(row.get("region_id", "")),
                "group_members": str(row.get("resolution_group_members", "")),
                "best_region_score": numeric_score,
                "member_scores": [],
                "member_region_ids": [],
                "resolution_tier": str(row.get("resolution_tier", "low_resolution")),
                "manual_review_recommended": bool(row.get("manual_review_recommended", True)),
                "group_plausibility_tier": str(row.get("group_plausibility_tier", "not_calibrated")),
            },
        )
        item["member_scores"].append(numeric_score)
        item["member_region_ids"].append(str(row.get("region_id", "")))
        item["manual_review_recommended"] = bool(item["manual_review_recommended"]) or bool(
            row.get("manual_review_recommended", False)
        )
        if numeric_score > float(item["best_region_score"]):
            item["best_region_score"] = numeric_score
            item["best_region_id"] = str(row.get("region_id", ""))
            item["resolution_tier"] = str(row.get("resolution_tier", "low_resolution"))
            item["group_plausibility_tier"] = str(row.get("group_plausibility_tier", "not_calibrated"))

    group_rows = []
    for item in group_scores.values():
        scores = [float(x) for x in item.pop("member_scores") if math.isfinite(float(x))]
        if not scores:
            aggregate = float("-inf")
            mean_score = float("nan")
        else:
            # Conservative evidence pooling: reward the best matching region, with a small
            # mean-score contribution from additional returned members in the same group.
            mean_score = float(sum(scores) / len(scores))
            aggregate = float(max(scores) + 0.10 * mean_score)
        item["group_score"] = aggregate
        item["mean_member_score"] = mean_score
        item["n_returned_group_members"] = int(len(scores))
        group_rows.append(item)
    group_rows.sort(key=lambda item: float(item["group_score"]), reverse=True)
    for rank, item in enumerate(group_rows, start=1):
        item["rank"] = int(rank)

    top = rows[0]
    metadata = annotated.setdefault("meta", {})
    metadata["region_resolution_annotation"] = {
        "enabled": True,
        "model_path": str(model_path),
        "primary_network": primary_network,
        "top1_resolution_tier": top["resolution_tier"],
        "top1_resolution_group": top["resolution_group"],
        "top1_group_members": top["resolution_group_members"],
        "top1_resolution_reasons": top["resolution_reasons"],
        "top1_group_plausibility_tier": top["group_plausibility_tier"],
        "top1_group_calibration_flags": top["group_calibration_flags"],
        "manual_review_recommended": bool(top["manual_review_recommended"]),
        "interpretation": "Resolution status annotates Region granularity only and does not rerank exact Region candidates.",
        "group_ranking_method": "best_region_score_plus_0p10_mean_returned_member_score",
        "group_ranking": group_rows[:10],
    }
    return annotated
