from __future__ import annotations

from pathlib import Path

import scripts.audit_public_provenance_paths as public_paths
import scripts.validate_csv_assets as csv_assets


def _missing_git(*_args, **_kwargs):
    raise FileNotFoundError("git")


def test_public_path_audit_discovers_public_tree_without_git(monkeypatch, tmp_path: Path) -> None:
    (tmp_path / "README.md").write_text("portable source image\n", encoding="utf-8")
    monkeypatch.setattr(public_paths, "ROOT", tmp_path)
    monkeypatch.setattr(public_paths.subprocess, "run", _missing_git)
    assert list(public_paths.tracked_paths()) == [tmp_path / "README.md"]


def test_csv_validation_discovers_repository_csvs_without_git(monkeypatch, tmp_path: Path) -> None:
    tracked = tmp_path / "reproducibility" / "tracked.csv"
    tracked.parent.mkdir()
    tracked.write_text("column\nvalue\n", encoding="utf-8")
    ignored = tmp_path / "external_data" / "raw.csv"
    ignored.parent.mkdir()
    ignored.write_text("column\nvalue\n", encoding="utf-8")
    monkeypatch.setattr(csv_assets, "ROOT", tmp_path)
    monkeypatch.setattr(csv_assets.subprocess, "run", _missing_git)
    assert csv_assets.tracked_csv_paths() == [tracked]
