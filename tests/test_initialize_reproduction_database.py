from __future__ import annotations

import sqlite3
from pathlib import Path

from scripts.initialize_reproduction_database import initialize


def test_reproduction_database_is_schema_only(tmp_path: Path) -> None:
    path = tmp_path / "reproduction.db"
    initialize(path)
    with sqlite3.connect(path) as conn:
        assert conn.execute("SELECT COUNT(*) FROM atlas_versions").fetchone()[0] == 0
        assert conn.execute("SELECT COUNT(*) FROM reference_expression").fetchone()[0] == 0
