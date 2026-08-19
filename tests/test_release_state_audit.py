from __future__ import annotations

from scripts.audit_release_state import audit


def test_engineering_phase_classifies_current_candidate_wording_without_stale_state() -> None:
    payload = audit()
    assert payload["release_state"] == "engineering_pre_finalization"
    assert payload["status"] == "PASS"
    assert payload["counts"]["STALE_CURRENT_RELEASE_STATE"] == 0
    assert payload["counts"]["STALE_CURRENT_REFERENCE"] == 0
    assert payload["counts"]["AMBIGUOUS"] == 0
