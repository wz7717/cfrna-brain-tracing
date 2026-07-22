from __future__ import annotations

import json
import sqlite3

import numpy as np
import pandas as pd
import pytest

from app.shared import DB_PATH
import core.bo2023_region_tracing as bo2023_region_tracing
from core.bo2023_region_tracing import (
    DEFAULT_BO2023_COUNTS,
    DEFAULT_BO2023_GENE_MAP,
    DEFAULT_BO2023_SAMPLE_INFO,
    ROUTE_NAME,
    trace_bo2023_secondary_regions,
)
from core.network_tracing import load_network_model


ROOT = DEFAULT_BO2023_COUNTS.parents[1]
STAGED_BO2023_INPUTS = (
    ROOT / "tests" / "controlled_data" / "bo2023" / "mfas5_819samples_28415genes_featurecounts_counts.txt",
    ROOT / "tests" / "controlled_data" / "bo2023" / "Information of sequenced samples_update_full878_filter819.xlsx",
    ROOT
    / "tests"
    / "controlled_data"
    / "bo2023"
    / "04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv",
)
DEFAULT_CONTROLLED_BO2023_INPUTS = (
    DEFAULT_BO2023_COUNTS,
    DEFAULT_BO2023_SAMPLE_INFO,
    DEFAULT_BO2023_GENE_MAP,
)


def _controlled_bo2023_inputs_or_skip() -> tuple:
    for paths in (STAGED_BO2023_INPUTS, DEFAULT_CONTROLLED_BO2023_INPUTS):
        if all(path.exists() for path in paths):
            return paths
    missing = [path for path in STAGED_BO2023_INPUTS + DEFAULT_CONTROLLED_BO2023_INPUTS if not path.exists()]
    if missing:
        pytest.skip(
            "Controlled Bo2023 raw expression inputs are not included in the public release: "
            + ", ".join(str(path) for path in missing)
        )
    raise AssertionError("unreachable")


def test_bo2023_secondary_region_tracing_uses_network_top3_beam(monkeypatch):
    counts_path, sample_info_path, gene_map_path = _controlled_bo2023_inputs_or_skip()

    def load_controlled_reference(db_path: str, atlas_id: int):
        return bo2023_region_tracing._load_raw_logcpm_reference_matrix(  # noqa: SLF001
            counts_path,
            sample_info_path,
            gene_map_path,
        )

    monkeypatch.setattr(bo2023_region_tracing, "_load_reference_matrix", load_controlled_reference)

    model = load_network_model()
    expression = pd.DataFrame(
        {
            "gene_symbol": model["genes"],
            "tpm_value": model["reference"][:, 0],
        }
    )
    network_output = {
        "results": [
            {"network_id": str(model["networks"][0]), "rank": 1},
            {"network_id": str(model["networks"][1]), "rank": 2},
            {"network_id": str(model["networks"][2]), "rank": 3},
        ]
    }

    out = trace_bo2023_secondary_regions(expression, network_output, DB_PATH, atlas_id=4, topk=8)

    assert out["meta"]["method"] == ROUTE_NAME
    assert out["meta"]["candidate_region_source"] == "SaleemNetworks Top3 beam"
    assert out["meta"]["n_scoring_genes_top50"] <= 50
    assert out["meta"]["n_scoring_genes_top100"] <= 100
    assert len(out["results"]) == 8
    assert {"region_id", "score", "top50_corr_component", "top100_corr_component"}.issubset(out["results"][0])


def test_bo2023_secondary_region_tracing_rejects_missing_network_beam():
    out = trace_bo2023_secondary_regions(
        pd.DataFrame({"gene_symbol": ["A"], "tpm_value": [1.0]}),
        {"results": []},
        DB_PATH,
        atlas_id=4,
    )

    assert out["results"] == []
    assert out["meta"]["traceability"] == "insufficient"


def test_production_reference_rejects_malformed_package_without_silent_fallback(monkeypatch):
    packaged = pd.DataFrame({"N1::R1": [1.0]}, index=["GENE1"])
    network_map = {"N1::R1": "N1"}

    monkeypatch.setattr(
        bo2023_region_tracing,
        "_load_packaged_region_reference_matrix",
        lambda *args, **kwargs: (packaged, network_map, "packaged_region_logcpm_reference"),
    )
    monkeypatch.setattr(
        bo2023_region_tracing,
        "_load_raw_logcpm_reference_matrix",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("production inference must not prefer workstation-only raw inputs")
        ),
    )

    matrix, mapping, source = bo2023_region_tracing._load_reference_matrix("unused.db", None)  # noqa: SLF001

    assert matrix.empty
    assert mapping == {}
    assert source.startswith("canonical_formal_region_assets_invalid:")


def test_committed_formal_assets_pass_all_canonical_invariants():
    matrix, mapping, _ = bo2023_region_tracing._load_packaged_region_reference_matrix()  # noqa: SLF001
    beams = bo2023_region_tracing._load_formal_beam_gene_panels()  # noqa: SLF001

    bo2023_region_tracing._validate_formal_beam_gene_panels(beams, matrix, mapping)  # noqa: SLF001

    assert matrix.shape[1] == 110
    assert len(set(mapping.values())) == 10
    assert len(beams) == 120
    assert bo2023_region_tracing.packaged_formal_region_assets_available()


def test_raw_metadata_loader_applies_canonical_parent_networks(monkeypatch, tmp_path):
    source = pd.DataFrame(
        {
            "No.": ["a", "b", "c", "d"],
            "Region": ["10m", "10m", "V2", "V2"],
            "SaleemNetworks": [
                "Lateral Prefrontal Cortex",
                "Orbitomedial Prefrontal Cortex (OMPFC)",
                "Parietal, and Parieto-occipital region",
                "Occipital/Temporal",
            ],
        }
    )
    monkeypatch.setattr(pd, "read_excel", lambda *args, **kwargs: source.copy())

    metadata = bo2023_region_tracing._read_bo2023_sample_metadata(tmp_path / "metadata.xlsx")  # noqa: SLF001

    assert metadata.loc[metadata["region_id"].eq("10m"), "network_id"].unique().tolist() == [
        "Orbitomedial Prefrontal Cortex (OMPFC)"
    ]
    assert metadata.loc[metadata["region_id"].eq("V2"), "network_id"].unique().tolist() == [
        "Occipital/Temporal"
    ]


def test_db_fallback_qualifies_and_canonicalizes_region_identity(tmp_path):
    db_path = tmp_path / "reference.db"
    with sqlite3.connect(db_path) as conn:
        conn.execute("CREATE TABLE reference_expression (atlas_id INTEGER, gene_symbol TEXT, region_id TEXT, avg_tpm REAL)")
        conn.execute("CREATE TABLE macaque_brain_atlas (atlas_id INTEGER, region_id TEXT, coordinates TEXT)")
        conn.execute("INSERT INTO reference_expression VALUES (1, 'GENE1', '10m', 1.0)")
        conn.executemany(
            "INSERT INTO macaque_brain_atlas VALUES (1, '10m', ?)",
            [
                (json.dumps({"saleem_network": "Lateral Prefrontal Cortex"}),),
                (json.dumps({"saleem_network": "Orbitomedial Prefrontal Cortex (OMPFC)"}),),
            ],
        )

    matrix, mapping, source = bo2023_region_tracing._load_db_reference_matrix(str(db_path), 1)  # noqa: SLF001

    key = "Orbitomedial Prefrontal Cortex (OMPFC)::10m"
    assert matrix.columns.tolist() == [key]
    assert mapping == {key: "Orbitomedial Prefrontal Cortex (OMPFC)"}
    assert source == "db_reference_expression_avg_tpm_fallback"


def test_canonical_reference_guard_rejects_112_region_semantics():
    networks = [f"N{i}" for i in range(10)]
    regions = [f"{networks[i % 10]}::R{i}" for i in range(110)]
    regions.extend(["N1::R0", "N2::R1"])
    matrix = pd.DataFrame([range(112)], index=["GENE1"], columns=regions)
    mapping = {region: region.split("::", 1)[0] for region in regions}

    with pytest.raises(ValueError, match="expected 110 regions"):
        bo2023_region_tracing._validate_canonical_region_reference(matrix, mapping)  # noqa: SLF001


def test_canonical_reference_guard_rejects_cross_network_duplicate_display_region():
    networks = [f"N{i}" for i in range(10)]
    regions = [f"{networks[i % 10]}::R{i}" for i in range(108)]
    regions.extend(["N0::DUP", "N1::DUP"])
    matrix = pd.DataFrame([range(110)], index=["GENE1"], columns=regions)
    mapping = {region: region.split("::", 1)[0] for region in regions}

    with pytest.raises(ValueError, match="more than one Network"):
        bo2023_region_tracing._validate_canonical_region_reference(matrix, mapping)  # noqa: SLF001


def test_formal_beam_guard_rejects_incomplete_beam_set():
    matrix, mapping, _ = bo2023_region_tracing._load_packaged_region_reference_matrix()  # noqa: SLF001
    beams = dict(bo2023_region_tracing._load_formal_beam_gene_panels())  # noqa: SLF001
    beams.pop(next(iter(beams)))

    with pytest.raises(ValueError, match="all 120 canonical Network beams"):
        bo2023_region_tracing._validate_formal_beam_gene_panels(beams, matrix, mapping)  # noqa: SLF001


def test_formal_group_and_exact_rankings_are_independent():
    rows = [
        {
            "region_id": "R1",
            "resolution_group": "G1",
            "resolution_group_members": "R1",
            "exact_local_score": 3.0,
            "group_local_score": 0.1,
            "resolution_tier": "high_resolution",
            "manual_review_recommended": False,
        },
        {
            "region_id": "R2",
            "resolution_group": "G2",
            "resolution_group_members": "R2",
            "exact_local_score": 1.0,
            "group_local_score": 0.9,
            "resolution_tier": "high_resolution",
            "manual_review_recommended": False,
        },
    ]

    exact = bo2023_region_tracing._rank_exact_regions(rows, topk=2)  # noqa: SLF001
    group_input = sorted(rows, key=lambda row: row["group_local_score"], reverse=True)
    groups = bo2023_region_tracing._rank_resolution_groups(group_input)  # noqa: SLF001

    assert [row["region_id"] for row in exact] == ["R1", "R2"]
    assert [row["resolution_group"] for row in groups] == ["G2", "G1"]
    assert groups[0]["group_score"] == 0.9


def test_formal_route_reports_and_executes_locked_panel_sizes(monkeypatch):
    genes = pd.Index([f"G{i}" for i in range(200)])
    matrix = pd.DataFrame(
        {
            "N1::R1": np.linspace(0.0, 1.0, 200),
            "N2::R2": np.linspace(1.0, 0.0, 200),
        },
        index=genes,
    )
    mapping = {"N1::R1": "N1", "N2::R2": "N2"}
    monkeypatch.setattr(
        bo2023_region_tracing,
        "_load_reference_matrix",
        lambda *args, **kwargs: (matrix, mapping, "controlled_reference"),
    )
    monkeypatch.setattr(
        bo2023_region_tracing,
        "_formal_beam_gene_order",
        lambda *args, **kwargs: (np.arange(200, dtype=int), "controlled_top200"),
    )
    monkeypatch.setattr(bo2023_region_tracing, "load_region_resolution_model", lambda: {})
    monkeypatch.setattr(bo2023_region_tracing, "_load_formal_beam_gene_panels", lambda: {})
    expression = pd.DataFrame({"gene_symbol": genes, "log_tpm": np.linspace(0.0, 1.0, 200)})
    network_output = {"results": [{"network_id": "N1"}, {"network_id": "N2"}, {"network_id": "N3"}]}

    out = trace_bo2023_secondary_regions(expression, network_output, "unused.db", None)

    assert out["meta"]["network_top_k"] == bo2023_region_tracing.NETWORK_TOP_K == 3
    assert out["meta"]["n_local_candidate_genes"] == bo2023_region_tracing.DEFAULT_LOCAL_TOP_N_GENES == 200
    assert out["meta"]["n_scoring_genes_top50"] == bo2023_region_tracing.EXACT_TOP50_GENE_COUNT == 50
    assert out["meta"]["n_scoring_genes_top100"] == bo2023_region_tracing.EXACT_TOP100_GENE_COUNT == 100
    assert out["meta"]["min_required_region_overlap_genes"] == bo2023_region_tracing.MIN_REGION_GENE_OVERLAP == 20


def test_formal_route_uses_locked_minimum_region_overlap(monkeypatch):
    genes = pd.Index([f"G{i}" for i in range(200)])
    matrix = pd.DataFrame(
        {"N1::R1": np.arange(200), "N2::R2": np.arange(200)[::-1]},
        index=genes,
    )
    mapping = {"N1::R1": "N1", "N2::R2": "N2"}
    monkeypatch.setattr(
        bo2023_region_tracing,
        "_load_reference_matrix",
        lambda *args, **kwargs: (matrix, mapping, "controlled_reference"),
    )
    expression = pd.DataFrame({"gene_symbol": genes[:19], "log_tpm": np.arange(19)})
    network_output = {"results": [{"network_id": "N1"}, {"network_id": "N2"}, {"network_id": "N3"}]}

    out = trace_bo2023_secondary_regions(expression, network_output, "unused.db", None)

    assert out["results"] == []
    assert out["meta"]["n_overlap_genes"] == 19
    assert out["meta"]["min_required_region_overlap_genes"] == 20
