from __future__ import annotations

from scripts.audit_public_provenance_paths import ROOT, audit, tracked_paths


def test_public_provenance_scope_includes_xml_artifacts() -> None:
    scanned = {path.relative_to(ROOT).as_posix() for path in tracked_paths()}
    assert "reproducibility/coverage.xml" in scanned


def test_current_public_provenance_has_no_private_machine_paths() -> None:
    payload = audit()
    assert payload["status"] == "PASS"
    assert payload["CURRENT_PUBLIC_ABSOLUTE_LOCAL_PATH_MATCHES"] == 0
