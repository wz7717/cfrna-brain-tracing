#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"
DEFAULT_RESOLUTION_MODEL = ROOT / "data" / "models" / "bo2023_region_resolution_groups.json"
DEFAULT_OUT = ROOT / "docs" / "deliverables" / "bo2023_public_anatomical_label_mapping.csv"


LOBE_PUBLIC_LABELS = {
    "Frontal": ("frontal lobe", "frontal cortex"),
    "Temporal": ("temporal lobe", "temporal cortex"),
    "Parietal": ("parietal lobe", "parietal cortex"),
    "Occipital": ("occipital lobe", "occipital cortex"),
    "Cingulate": ("cingulate cortex", "cingulate gyrus"),
    "Insula": ("insula", "insular cortex"),
    "Subcortical": ("subcortical structures", "deep brain structures"),
}

NETWORK_PUBLIC_LABELS = {
    "Orbitomedial Prefrontal Cortex (OMPFC)": "orbitomedial prefrontal cortex",
    "Lateral Prefrontal Cortex": "lateral prefrontal cortex",
    "Frontal (agranular frontal motor areas)": "primary/premotor frontal cortex",
    "Cingulate gyrus": "cingulate cortex",
    "Parietal, and Parieto-occipital region": "parietal/parieto-occipital cortex",
    "Occipital/Temporal": "visual occipital-temporal cortex",
    "Temporal": "temporal cortex",
    "Operculum/Insula": "operculum/insula",
    "Subcortical": "subcortical structures",
    "Hippocampal formation": "hippocampal formation",
}

EXACT_PUBLIC_LABELS = {
    "HC": ("hippocampus", "hippocampal formation; hippocampus"),
    "amy": ("amygdala", "amygdala"),
    "thalamus": ("thalamus", "thalamus"),
    "cd": ("caudate", "caudate nucleus; caudate"),
    "pu": ("putamen", "putamen"),
    "GPeGPi": ("globus pallidus", "globus pallidus; external globus pallidus; internal globus pallidus"),
    "Pir": ("piriform cortex", "piriform cortex"),
    "cla": ("claustrum", "claustrum"),
    "V1": ("primary visual cortex", "primary visual cortex; V1; occipital cortex"),
    "V2": ("secondary visual cortex", "V2; visual area 2; occipital cortex"),
    "V3d": ("visual cortex area V3", "V3; dorsal V3; occipital cortex"),
    "V3v": ("visual cortex area V3", "V3; ventral V3; occipital cortex"),
    "V4": ("visual cortex area V4", "V4; visual area 4; occipital/temporal visual cortex"),
    "V4v": ("visual cortex area V4", "V4; ventral V4; occipital/temporal visual cortex"),
    "SII": ("secondary somatosensory cortex", "S2; secondary somatosensory cortex; parietal operculum"),
    "3a/b": ("primary somatosensory cortex", "primary somatosensory cortex; areas 3a/3b"),
    "F1": ("primary motor cortex", "M1; primary motor cortex"),
    "F2": ("dorsal premotor cortex", "PMd; dorsal premotor cortex"),
    "F3": ("supplementary motor area", "SMA; supplementary motor area"),
    "F4": ("ventral premotor cortex", "PMv; ventral premotor cortex"),
    "F5": ("ventral premotor cortex", "PMv; ventral premotor cortex"),
    "F6": ("pre-supplementary motor area", "pre-SMA; presupplementary motor area"),
    "F7": ("dorsal premotor cortex", "PMd; rostral dorsal premotor cortex"),
    "G": ("gustatory cortex", "gustatory cortex; insula"),
    "STG": ("superior temporal gyrus", "superior temporal gyrus; temporal cortex"),
    "TG": ("temporal pole", "temporal pole"),
    "TF": ("parahippocampal cortex", "parahippocampal cortex"),
    "TFO": ("parahippocampal cortex", "parahippocampal cortex"),
}


def load_atlas(db_path: Path, atlas_id: int) -> pd.DataFrame:
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            """
            SELECT region_id, region_name, coordinates
            FROM macaque_brain_atlas
            WHERE atlas_id = ?
            ORDER BY region_id
            """,
            conn,
            params=(int(atlas_id),),
        )
    finally:
        conn.close()
    if df.empty:
        raise ValueError(f"no macaque_brain_atlas rows for atlas_id={atlas_id}")
    coordinates = df["coordinates"].map(lambda value: json.loads(value) if isinstance(value, str) and value else {})
    df["bo2023_lobe"] = coordinates.map(lambda item: item.get("lobe", ""))
    df["bo2023_saleem_network"] = coordinates.map(lambda item: item.get("saleem_network", ""))
    df["bo2023_neocortex_flag"] = coordinates.map(lambda item: item.get("neocortex_flag", ""))
    df["bo2023_roi173"] = coordinates.map(lambda item: item.get("roi173", ""))
    df["bo2023_source_sample_count"] = coordinates.map(lambda item: item.get("source_sample_count", ""))
    return df.drop(columns=["coordinates"])


def load_resolution_entries(model_path: Path) -> dict[str, dict[str, Any]]:
    if not model_path.exists():
        return {}
    model = json.loads(model_path.read_text(encoding="utf-8"))
    return model.get("entries", {})


def is_prefrontal(region_id: str, region_name: str, network: str) -> bool:
    text = f"{region_id} {region_name} {network}".lower()
    return "prefrontal" in text or region_id in {
        "10m",
        "10o",
        "11l",
        "11m",
        "12l",
        "12m",
        "12o",
        "12r",
        "13a",
        "13b",
        "13l",
        "13m",
        "14c",
        "14r",
        "25",
        "32",
        "44",
        "45",
        "46d",
        "46v",
        "8A",
        "8Bd",
        "8Bm",
        "8Bs",
        "9d",
        "9m",
    }


def infer_mapping(row: pd.Series) -> dict[str, Any]:
    region_id = str(row["region_id"])
    region_name = str(row["region_name"])
    lobe = str(row["bo2023_lobe"])
    network = str(row["bo2023_saleem_network"])
    public_lobe, public_lobe_synonym = LOBE_PUBLIC_LABELS.get(lobe, (lobe.lower(), lobe.lower()))
    public_network = NETWORK_PUBLIC_LABELS.get(network, network.lower())

    public_fine = region_name
    synonyms = f"{region_id}; {region_name}; {public_network}; {public_lobe_synonym}"
    mapping_level = "lobe/network"
    confidence = "medium"
    rule = "mapped by Bo2023 lobe and Saleem network; exact cross-species area identity is not assumed"

    if region_id in EXACT_PUBLIC_LABELS:
        public_fine, synonyms = EXACT_PUBLIC_LABELS[region_id]
        mapping_level = "exact/common_anatomical_name"
        confidence = "high"
        rule = "mapped by widely used anatomical name shared across macaque and public brain datasets"
    elif is_prefrontal(region_id, region_name, network):
        public_fine = public_network
        synonyms = f"{region_id}; {region_name}; prefrontal cortex; {public_network}; frontal cortex"
        mapping_level = "network/subdivision"
        confidence = "medium"
        rule = "mapped to prefrontal subdivision by Saleem network and region name"
    elif lobe in {"Cingulate", "Insula"}:
        public_fine = public_network
        synonyms = f"{region_id}; {region_name}; {public_network}; {public_lobe_synonym}"
        mapping_level = "major_structure"
        confidence = "medium"
        rule = "mapped to major cortical structure by Bo2023 lobe/network"
    elif lobe == "Subcortical":
        public_fine = public_network
        synonyms = f"{region_id}; {region_name}; subcortical; {public_network}"
        mapping_level = "major_structure"
        confidence = "low" if region_id not in EXACT_PUBLIC_LABELS else confidence
        rule = "subcortical abbreviation lacks stable public-dataset synonym; use only as coarse structure unless manually curated"
    elif lobe in {"Temporal", "Parietal", "Occipital", "Frontal"}:
        public_fine = public_network
        synonyms = f"{region_id}; {region_name}; {public_network}; {public_lobe_synonym}"

    compatible_public_datasets = [
        "AHBA/Allen Human Brain Atlas anatomical labels",
        "GTEx tissue-level brain labels",
        "TCGA/CGGA tumor location labels when available",
        "Ivy GAP only as tumor microanatomic grouping, not normal-region ground truth",
    ]
    return {
        "public_major_anatomy": public_lobe,
        "public_network_or_subdivision": public_network,
        "public_fine_anatomy_label": public_fine,
        "public_label_synonyms": synonyms,
        "mapping_level": mapping_level,
        "mapping_confidence": confidence,
        "mapping_rule": rule,
        "compatible_public_dataset_labels": " | ".join(compatible_public_datasets),
        "ivy_gap_mapping_status": "no_direct_normal_brain_region_mapping",
        "do_not_use_as_exact_ground_truth": True,
        "caution_note": (
            "Use as anatomical label harmonization. Do not interpret public human/tumor labels as exact "
            "Bo2023 macaque region ground truth without dataset-specific manual review."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Export Bo2023-to-public anatomical label harmonization table.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--atlas-id", type=int, default=4)
    parser.add_argument("--resolution-model", type=Path, default=DEFAULT_RESOLUTION_MODEL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    atlas = load_atlas(args.db, args.atlas_id)
    resolution_entries = load_resolution_entries(args.resolution_model)
    rows: list[dict[str, Any]] = []
    for item in atlas.to_dict(orient="records"):
        row = pd.Series(item)
        mapping = infer_mapping(row)
        region_id = str(row["region_id"])
        network = str(row["bo2023_saleem_network"])
        resolution_entry = resolution_entries.get(f"{network}||{region_id}", {})
        rows.append(
            {
                "bo2023_region_id": region_id,
                "bo2023_region_name": row["region_name"],
                "bo2023_lobe": row["bo2023_lobe"],
                "bo2023_saleem_network": network,
                "bo2023_neocortex_flag": row["bo2023_neocortex_flag"],
                "bo2023_roi173": row["bo2023_roi173"],
                "bo2023_source_sample_count": row["bo2023_source_sample_count"],
                "resolution_tier": resolution_entry.get("resolution_tier", ""),
                "resolution_group": resolution_entry.get("resolution_group", region_id),
                "resolution_group_members": " | ".join(map(str, resolution_entry.get("group_members", [region_id]))),
                "group_plausibility_tier": resolution_entry.get("group_plausibility_tier", ""),
                **mapping,
            }
        )

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, encoding="utf-8-sig")
    summary = {
        "output": str(args.out),
        "n_rows": int(len(out)),
        "mapping_confidence_counts": out["mapping_confidence"].value_counts().to_dict(),
        "mapping_level_counts": out["mapping_level"].value_counts().to_dict(),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
