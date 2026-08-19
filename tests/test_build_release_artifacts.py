from __future__ import annotations

import zipfile
from pathlib import Path

from scripts.build_release_artifacts import build_zip


def test_deterministic_zip_uses_sorted_members_and_fixed_timestamp(tmp_path: Path) -> None:
    (tmp_path / "b.txt").write_text("b\n", encoding="utf-8")
    (tmp_path / "a.txt").write_text("a\n", encoding="utf-8")
    first = tmp_path / "first.zip"
    second = tmp_path / "second.zip"
    first_record = build_zip(tmp_path, ["b.txt", "a.txt"], first, "bundle")
    second_record = build_zip(tmp_path, ["a.txt", "b.txt"], second, "bundle")
    assert first_record["sha256"] == second_record["sha256"]
    with zipfile.ZipFile(first) as archive:
        assert archive.namelist() == ["bundle/a.txt", "bundle/b.txt"]
        assert archive.getinfo("bundle/a.txt").date_time == (1980, 1, 1, 0, 0, 0)
