#!/usr/bin/env python
"""Initialize an empty BrainTrace schema for a temporary reproduction database.

Unlike ``database_init.py``'s interactive/demo entry point, this command never
adds demonstration atlas rows or synthetic expression.  It exists so full
reproduction can build a clean Bo2023 database from verified raw inputs.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from database_init import CSFRNASourceDatabase  # noqa: E402


def initialize(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite reproduction database: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    database = CSFRNASourceDatabase(str(path))
    try:
        database.connect()
        database.create_database_schema()
    finally:
        database.close()
    with sqlite3.connect(path) as conn:
        tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        forbidden_population = sum(
            int(conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
            for table in ("macaque_brain_atlas", "reference_expression", "atlas_versions")
        )
    required = {"macaque_brain_atlas", "reference_expression", "atlas_versions", "signature_sets"}
    if not required.issubset(tables) or forbidden_population != 0:
        raise RuntimeError("temporary reproduction database schema is incomplete or unexpectedly populated")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    args = parser.parse_args()
    initialize(args.db)
    print("reproduction database schema initialized")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
