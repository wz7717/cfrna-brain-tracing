from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

from scripts.build_release_artifacts import build_zip, member_payload


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


def test_git_blob_payload_ignores_checkout_line_endings(tmp_path: Path) -> None:
    def git(*args: str) -> None:
        subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "test@example.invalid")
    git("config", "user.name", "Release test")
    member = tmp_path / "member.txt"
    member.write_bytes(b"line\n")
    git("add", "member.txt")
    git("commit", "-qm", "fixture")
    member.write_bytes(b"line\r\n")

    archive_path = tmp_path / "canonical.zip"
    build_zip(
        tmp_path,
        ["member.txt"],
        archive_path,
        "bundle",
        payload_reader=lambda relative: member_payload(tmp_path, relative, set()),
    )
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.read("bundle/member.txt") == b"line\n"
