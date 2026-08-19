from __future__ import annotations

from pathlib import Path

from scripts.audit_lfs_materialization import is_lfs_pointer


def test_lfs_pointer_detection(tmp_path: Path) -> None:
    pointer = tmp_path / "pointer.bin"
    pointer.write_bytes(b"version https://git-lfs.github.com/spec/v1\noid sha256:abc\n")
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"\x89PNG\r\n\x1a\n")
    assert is_lfs_pointer(pointer) is True
    assert is_lfs_pointer(payload) is False
