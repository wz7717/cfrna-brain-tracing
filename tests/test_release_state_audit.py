from __future__ import annotations

from scripts.audit_release_state import _is_release_context, audit


def test_final_phase_has_no_stale_release_state_wording() -> None:
    payload = audit()
    assert payload["release_state"] == "final"
    assert payload["status"] == "PASS"
    assert payload["counts"]["STALE_CURRENT_RELEASE_STATE"] == 0
    assert payload["counts"]["STALE_CURRENT_REFERENCE"] == 0
    assert payload["counts"]["AMBIGUOUS"] == 0


def test_candidate_ranking_remains_scientific_context_with_release_metadata_on_same_line() -> None:
    assert not _is_release_context("candidate", "The current release is a candidate ranking tool.")
    assert _is_release_context("candidate", "The v0.1.18 patch candidate is unreleased.")
