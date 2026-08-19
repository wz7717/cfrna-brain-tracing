from __future__ import annotations

from scripts.audit_public_provenance_paths import audit


def test_current_public_provenance_has_no_private_machine_paths() -> None:
    payload = audit()
    assert payload["status"] == "PASS"
    assert payload["CURRENT_PUBLIC_ABSOLUTE_LOCAL_PATH_MATCHES"] == 0
