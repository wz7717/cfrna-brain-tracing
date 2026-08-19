from __future__ import annotations

import reproduce_all


def test_git_probe_is_safe_when_git_is_not_installed(monkeypatch) -> None:
    def missing_git(*_args, **_kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(reproduce_all.subprocess, "run", missing_git)
    assert reproduce_all._git_value("rev-parse", "HEAD") == "unavailable"


def test_environment_record_uses_container_build_provenance(monkeypatch) -> None:
    monkeypatch.setenv("BRAINTRACE_GIT_SHA", "0123456789abcdef")
    monkeypatch.setenv("BRAINTRACE_GIT_CLEAN", "true")
    monkeypatch.setattr(reproduce_all, "_git_value", lambda *_args: "unavailable")
    record = reproduce_all.environment_record()
    assert record["git_sha"] == "0123456789abcdef"
    assert record["git_clean"] is True
