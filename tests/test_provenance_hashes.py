from __future__ import annotations

from pathlib import Path

from core.provenance_hashes import sha256_utf8_lf_text


def test_repository_text_sha_is_line_ending_invariant(tmp_path: Path) -> None:
    path = tmp_path / "artifact.json"
    path.write_bytes(b'{"value": 1}\r\n')
    windows_digest = sha256_utf8_lf_text(path)
    path.write_bytes(b'{"value": 1}\n')
    assert sha256_utf8_lf_text(path) == windows_digest
