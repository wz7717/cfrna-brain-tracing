from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

import pytest

from scripts.build_release_artifacts import build_zip, enforce_public_release_content, member_payload


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
    member = tmp_path / "member.txt"
    member.write_bytes(b"line\r\n")

    def git_show(args: list[str], **_: object) -> subprocess.CompletedProcess[bytes]:
        assert args == ["git", "show", "HEAD:member.txt"]
        return subprocess.CompletedProcess(args, 0, stdout=b"line\n", stderr=b"")

    import scripts.build_release_artifacts as release_artifacts

    original_run = release_artifacts.subprocess.run
    release_artifacts.subprocess.run = git_show
    try:
        archive_path = tmp_path / "canonical.zip"
        build_zip(
            tmp_path,
            ["member.txt"],
            archive_path,
            "bundle",
            payload_reader=lambda relative: member_payload(tmp_path, relative, set()),
        )
    finally:
        release_artifacts.subprocess.run = original_run

    with zipfile.ZipFile(archive_path) as archive:
        assert archive.read("bundle/member.txt") == b"line\n"


def test_release_builder_rejects_submission_authoring_tool(tmp_path: Path) -> None:
    script = tmp_path / "scripts" / "update_main_manuscript.py"
    script.parent.mkdir()
    script.write_text("print('local submission utility')\n", encoding="utf-8")

    with pytest.raises(ValueError, match="public release content audit failed"):
        enforce_public_release_content(tmp_path)
