from __future__ import annotations

import inspect
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import scripts.run_huang2025_external_candidate as huang


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reproducibility" / "huang_2025"


def test_parse_sample_id_exposes_no_inferred_patient_identity() -> None:
    csf = huang.parse_sample_id("GLI_CSF16")
    plasma = huang.parse_sample_id("GLI_plasma16")

    assert csf["disease_group"] == "glioma"
    assert csf["fluid"] == "CSF"
    assert plasma["fluid"] == "plasma"
    assert csf["tumor_status"] == plasma["tumor_status"] == "tumour"
    assert pd.isna(csf["patient_id"])
    assert pd.isna(plasma["patient_id"])
    assert csf["patient_id_status"] == plasma["patient_id_status"] == huang.PATIENT_ID_STATUS
    assert "patient_key" not in csf
    assert "patient_number" not in csf


def test_parse_sample_id_rejects_unknown_contract() -> None:
    with pytest.raises(ValueError, match="Unexpected Huang"):
        huang.parse_sample_id("GBM_CSF1")


def test_log2_rpm_conversion_is_exact_scale_change() -> None:
    source = np.array([0.0, 1.0, 2.0, 3.5])
    converted = huang.source_log2_rpm_to_log1p_rpm(source)
    np.testing.assert_allclose(converted, source * math.log(2.0), rtol=1e-7)


def test_log2_rpm_conversion_rejects_negative_values() -> None:
    with pytest.raises(ValueError, match="negative"):
        huang.source_log2_rpm_to_log1p_rpm(np.array([0.0, -0.1]))


def test_benjamini_hochberg_is_monotone_by_p_value() -> None:
    p_values = [0.04, 0.001, 0.02]
    adjusted = huang.benjamini_hochberg(p_values)
    ordered = [adjusted[index] for index in np.argsort(p_values)]
    assert ordered == sorted(ordered)
    assert all(p <= q <= 1.0 for p, q in zip(p_values, adjusted))


def test_tumour_control_test_is_profile_level_when_pairing_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(huang, "N_BOOTSTRAP", 20)
    rows = []
    for fluid in ("CSF", "plasma"):
        for status, offset in (("tumour", 1.0), ("control", 0.0)):
            for index in range(3):
                rows.append(
                    {
                        "fluid": fluid,
                        "tumor_status": status,
                        "atlas_fit_score": offset + index / 10,
                        "network_margin": offset + index / 20,
                        "network_entropy": offset + index / 30,
                    }
                )
    comparisons = huang.compare_tumour_control(pd.DataFrame(rows))
    assert len(comparisons) == 6
    assert comparisons["test"].eq("two-sided Mann-Whitney U (profile-level; pairing unavailable)").all()
    assert comparisons["n_tumour"].eq(3).all()
    assert comparisons["n_control"].eq(3).all()
    assert comparisons["bh_fdr"].between(0, 1).all()


def test_canonical_outputs_have_full_cohort_accounting() -> None:
    summary = json.loads((OUT / "huang_2025_canonical_summary.json").read_text(encoding="utf-8"))
    ledger = pd.read_csv(OUT / "huang_2025_sample_ledger.csv")
    comparisons = pd.read_csv(OUT / "huang_2025_tumour_control_comparisons.csv")

    assert summary["protocol_status"] == "huang_2025_provenance_remediated"
    assert (summary["n_profiles"], summary["n_csf"], summary["n_plasma"]) == (159, 77, 82)
    assert summary["n_traceable_outputs"] == 159
    assert summary["patient_paired_analysis"] == "NOT_SUPPORTED"
    assert summary["synthetic_matched_admixture"] == "REMOVED_FROM_CANONICAL_ANALYSIS"
    assert len(ledger) == 159
    assert ledger["fluid"].eq("CSF").sum() == 77
    assert ledger["fluid"].eq("plasma").sum() == 82
    assert ledger["patient_id"].isna().all()
    assert ledger["patient_id_status"].eq(huang.PATIENT_ID_STATUS).all()
    assert ledger["BrainTrace_output_available"].all()
    assert len(comparisons) == 6
    assert set(comparisons["n_tumour"].astype(int)) == {59, 64}
    assert set(comparisons["n_control"].astype(int)) == {18}
    assert comparisons["test"].eq("two-sided Mann-Whitney U (profile-level; pairing unavailable)").all()
    assert summary["minimum_bh_fdr"] == pytest.approx(comparisons["bh_fdr"].min())


def test_canonical_distribution_denominators_and_percentages_are_arithmetic() -> None:
    fluid = pd.read_csv(OUT / "huang_2025_fluid_summary.csv")
    for _, row in fluid.iterrows():
        assert row["OMPFC_top1_denominator"] == row["n_profiles"]
        assert row["OMPFC_top3_denominator"] == row["n_profiles"]
        assert row["OMPFC_top1_percent"] == pytest.approx(
            100 * row["OMPFC_top1_numerator"] / row["OMPFC_top1_denominator"]
        )
        assert row["OMPFC_top3_percent"] == pytest.approx(
            100 * row["OMPFC_top3_numerator"] / row["OMPFC_top3_denominator"]
        )


def test_runner_source_has_no_pseudopairing_or_name_substitution() -> None:
    source = inspect.getsource(huang)
    forbidden = (
        "patient_key",
        'replace("_CSF"',
        "paired_stability",
        "permutation_pvalue",
        "matched_csf_plasma",
    )
    assert all(token not in source for token in forbidden)
    assert "profile-level; pairing unavailable" in source


def test_bio1_bio4_audit_no_longer_generates_huang_pseudopairs() -> None:
    source = (ROOT / "scripts" / "run_p1_bio1_4_audit.py").read_text(encoding="utf-8")
    forbidden = (
        "analyze_huang_and_admixture",
        "patient_key",
        "matched_csf_plasma",
        "_CSF\"",
        "plasma_fraction",
    )
    assert all(token not in source for token in forbidden)


def test_canonical_output_directory_contains_no_retired_pseudopair_artifact() -> None:
    names = {path.name.lower() for path in OUT.iterdir() if path.is_file()}
    assert not any("paired_csf_plasma" in name or "admixture" in name for name in names)
    report = (OUT / "HUANG_2025_RESULTS.md").read_text(encoding="utf-8")
    assert "No patient-level CSF-plasma correspondence was assumed." in report
    assert "minimum p 0.304" not in report



def test_asset_record_uses_portable_paths(tmp_path: Path) -> None:
    internal = huang.asset_record(
        ROOT / "data/models/bo2023_saleem_network_top200_model.npz",
        "locked_network_model",
    )
    assert internal["path"] == "data/models/bo2023_saleem_network_top200_model.npz"
    assert internal["path_kind"] == "repository_relative"
    assert not Path(internal["path"]).is_absolute()

    external_path = tmp_path / "huang_external.csv"
    external_path.write_text("x\n", encoding="utf-8")

    external = huang.asset_record(external_path, "external_test")
    assert external["path"] == "huang_external.csv"
    assert external["path_kind"] == "external_basename"
    assert not Path(external["path"]).is_absolute()



def test_generated_huang_report_does_not_assert_patient_independence() -> None:
    report = (OUT / "HUANG_2025_RESULTS.md").read_text(encoding="utf-8")

    forbidden = (
        "independent fluid-specific cohorts",
        "## Independent tumour-control diagnostics",
        "independent-cohort distributions",
    )

    assert all(token not in report for token in forbidden)
    assert "patient-level dependence" in report
