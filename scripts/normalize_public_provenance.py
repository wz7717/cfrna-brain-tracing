#!/usr/bin/env python
"""Normalize generated public provenance to portable symbolic locators.

This generator performs a metadata-only migration: source identities and all
scientific values are retained verbatim while machine-local path strings are
replaced by stable ``repository::`` or ``external_source::`` locators.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SIGNATURE_PROVENANCE = ROOT / "data" / "signatures" / "mb_cell_state_signatures.provenance.json"
PRIVATE_PATH_RE = re.compile(r"(?:[A-Za-z]:\\|C:\\Users\\|/Users/|/home/)", re.IGNORECASE)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def normalize_signature_provenance() -> None:
    payload = json.loads(SIGNATURE_PROVENANCE.read_text(encoding="utf-8"))
    if "signature_counts" not in payload or "selection" not in payload:
        raise ValueError("signature provenance does not have the expected generated schema")
    payload["source_file"] = "external_source::mb_literature_signature_sets/mb_literature_signature_gene_sets.csv"
    payload["output_file"] = "repository::data/signatures/mb_cell_state_signatures.csv"
    payload["provenance_locator_policy"] = "portable symbolic locators; source identity and generated values unchanged"
    payload["provenance_generator"] = "scripts/normalize_public_provenance.py"
    write_json(SIGNATURE_PROVENANCE, payload)


def verify() -> dict[str, Any]:
    paths = [SIGNATURE_PROVENANCE]
    matches: list[str] = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        if PRIVATE_PATH_RE.search(text):
            matches.append(path.relative_to(ROOT).as_posix())
    return {"status": "PASS" if not matches else "FAIL", "absolute_path_matches": matches}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if not args.verify_only:
        normalize_signature_provenance()
    payload = verify()
    print(json.dumps(payload))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
