"""Portable SHA-256 helpers for repository-staged provenance artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_utf8_lf_text(path: Path) -> str:
    """Hash UTF-8 repository text after canonicalizing line endings.

    Raw external-file identity is intentionally handled elsewhere.  This helper
    is only for tracked/staged text, whose Git content is LF regardless of a
    contributor's checkout setting.
    """

    text = path.read_text(encoding="utf-8")
    canonical = text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")
    return hashlib.sha256(canonical).hexdigest().upper()
