#!/usr/bin/env python
"""Semantically classify release-state wording and stale current-version text."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "release" / "v0.1.17" / "release_manifest.json"
STATE_PATTERNS = re.compile(
    r"candidate|release candidate|unreleased|pending(?: release)?|future release|after release|once released|planned|patch candidate|not yet published|will be updated|release requires|待发布|尚未发布|发布后|候选版本",
    re.IGNORECASE,
)
CURRENT_OLD_VERSION = re.compile(r"(?:current(?: executable| immutable| software| release| metadata)?|当前[^\n]{0,24})(?:[^\n]{0,80})\bv0\.1\.(?:8|9|1[0-6])\b", re.IGNORECASE)
SCIENTIFIC_CANDIDATE_CONTEXT = re.compile(r"candidate(?:[-\s]+)(?:ranking|set)", re.IGNORECASE)


def _is_release_context(term: str, line: str) -> bool:
    if term.lower() == "candidate" and SCIENTIFIC_CANDIDATE_CONTEXT.search(line):
        return False
    return term.lower() != "candidate" or bool(re.search(r"release|version|unreleased|patch|provenance", line, re.IGNORECASE))


def _scan_paths() -> list[Path]:
    candidates = [
        ROOT / "README.md",
        ROOT / "DATA_PROVENANCE.md",
        ROOT / "reproducibility" / "DATA_PROVENANCE.md",
        ROOT / "RELEASE_CANDIDATE_v0.1.17_HUANG.md",
        ROOT / ".zenodo.json",
        ROOT / "CITATION.cff",
    ]
    candidates.extend(sorted((ROOT / "release").rglob("*")) if (ROOT / "release").exists() else [])
    return [path for path in candidates if path.is_file()]


def audit(manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    version = str(manifest["version"])
    phase = str(manifest["release_state"])
    final_phase = phase == "final"
    rows: list[dict[str, Any]] = []
    for path in _scan_paths():
        text = path.read_text(encoding="utf-8", errors="replace")
        relative = path.relative_to(ROOT).as_posix()
        for match in STATE_PATTERNS.finditer(text):
            line_number = text.count("\n", 0, match.start()) + 1
            line = text.splitlines()[line_number - 1]
            term = match.group(0)
            # "candidate" is also a central scientific term (candidate
            # ranking/candidate set).  It is release-state wording only when
            # the same line supplies release/version/provenance context.
            release_context = _is_release_context(term, line)
            current_version_context = f"v{version}" in line or relative.endswith("RELEASE_CANDIDATE_v0.1.17_HUANG.md")
            historical_context = bool(re.search(r"historical|previous|prior|immutable old", line, re.IGNORECASE))
            if not release_context:
                classification = "VALID_RELEASE_STATE"
            elif historical_context:
                classification = "VALID_HISTORICAL_REFERENCE"
            elif not final_phase:
                classification = "VALID_PRE_RELEASE_STATE"
            elif final_phase and current_version_context:
                classification = "STALE_CURRENT_RELEASE_STATE"
            else:
                classification = "AMBIGUOUS"
            rows.append({"path": relative, "line": line_number, "term": term, "classification": classification})
        for match in CURRENT_OLD_VERSION.finditer(text):
            line_number = text.count("\n", 0, match.start()) + 1
            line = text.splitlines()[line_number - 1]
            classification = "STALE_CURRENT_REFERENCE" if final_phase else "VALID_PRE_RELEASE_STATE"
            rows.append({"path": relative, "line": line_number, "term": line.strip(), "classification": classification})
    counts = {name: sum(row["classification"] == name for row in rows) for name in {
        "VALID_HISTORICAL_REFERENCE", "VALID_PRE_RELEASE_STATE", "VALID_RELEASE_STATE",
        "STALE_CURRENT_RELEASE_STATE", "STALE_CURRENT_REFERENCE", "AMBIGUOUS",
    }}
    status = "PASS" if not counts["STALE_CURRENT_RELEASE_STATE"] and not counts["STALE_CURRENT_REFERENCE"] and not counts["AMBIGUOUS"] else "FAIL"
    return {
        "schema": "braintrace.release_state_audit.v1",
        "version": version,
        "release_state": phase,
        "status": status,
        "counts": counts,
        "records": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-text", type=Path)
    args = parser.parse_args()
    payload = audit(args.manifest)
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    if args.output_text:
        lines = [f"status={payload['status']}", f"version={payload['version']}", f"release_state={payload['release_state']}"]
        lines.extend(f"{key}={value}" for key, value in sorted(payload["counts"].items()))
        args.output_text.parent.mkdir(parents=True, exist_ok=True)
        args.output_text.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload, ensure_ascii=False))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
