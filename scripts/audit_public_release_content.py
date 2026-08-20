#!/usr/bin/env python3
"""Fail closed when manuscript-production tooling enters the public tree."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
CODE_SUFFIXES = {".py", ".r", ".ps1", ".sh", ".js", ".ts", ".ipynb", ".qmd", ".rmd"}
FALLBACK_IGNORED_PARTS = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    ".venv",
    ".venv-repro",
    "__pycache__",
    "external_data",
    "manuscript_remediation",
    "submission_seal_audit",
}

PROHIBITED_ARTIFACT_NAME = re.compile(
    r"(?:main[_ -]?manuscript|supplement(?:ary)?[_ -]?(?:file|material)?|"
    r"cover[_ -]?letter).*(?:\.docx|\.docm|\.pdf)$",
    re.IGNORECASE,
)
PROHIBITED_SCRIPT_NAME = re.compile(
    r"(?:^|/)(?:(?:apply|update|modify|edit|write|build|create|generate)[^/]*"
    r"(?:manuscript|supplement|submission|docx|change[_-]?list|checklist)|"
    r"[^/]*(?:manuscript|supplement|submission)[^/]*"
    r"(?:update|edit|format|layout|docx|writer))\.(?:py|r|ps1|sh|js|ts)$",
    re.IGNORECASE,
)
PROHIBITED_REVISION_RENDER = re.compile(
    r"(?:^|/)figure\d+[a-z]?_revision/.*\.(?:png|jpe?g|svg|pdf)$",
    re.IGNORECASE,
)
DOCX_IMPORT = re.compile(r"(?:from\s+docx\s+import|import\s+docx)", re.IGNORECASE)
DOCUMENT_SAVE = re.compile(r"\.save\s*\(")
MAIN_OR_SUPPLEMENT_TARGET = re.compile(
    r"Main_Manuscript|Supplementary_File|MAIN_SOURCE|SUPP_SOURCE|"
    r"update_main|update_supplement|submission documents?|主稿|主文|补充材料",
    re.IGNORECASE,
)
TEXT_OUTPUT = re.compile(r"\.write_text\s*\(|open\s*\([^\n]{0,200}['\"]w", re.IGNORECASE)
IMAGE_OUTPUT = re.compile(r"\.(?:savefig|save)\s*\(", re.IGNORECASE)
MANUSCRIPT_MEDIA_SOURCE = re.compile(
    r"manuscript[/\\].*(?:_extracted_imgs|word[/\\]media)",
    re.IGNORECASE,
)

# Split the literals so this policy implementation cannot match its own source.
AUTHORING_AID_LITERALS = (
    "Suggested " + "manuscript wording",
    "Main " + "manuscript remediation",
    "Supplement " + "remediation",
    "all manuscript/" + "supplement wording",
    "修改" + "清单",
    "限长" + "压缩",
    "Word 修订" + "模式",
    "真实 Word " + "修订",
    "Publication-ready " + "tables for",
    "## Reporting " + "rules",
)


def _git_files(root: Path) -> list[str] | None:
    try:
        result = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
            cwd=root,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    if result.returncode != 0:
        return None
    return sorted(
        item.decode("utf-8")
        for item in result.stdout.split(b"\0")
        if item and (root / item.decode("utf-8")).is_file()
    )


def _fallback_files(root: Path) -> list[str]:
    rows: list[str] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if any(part in FALLBACK_IGNORED_PARTS or part.endswith(".egg-info") for part in relative.parts):
            continue
        rows.append(relative.as_posix())
    return sorted(rows)


def public_files(root: Path = ROOT) -> list[str]:
    """Return the intended public working tree, with a gitless-container fallback."""

    git_files = _git_files(root)
    return git_files if git_files is not None else _fallback_files(root)


def classify(relative: str, raw: bytes) -> list[dict[str, str]]:
    """Classify one public file under the manuscript-production exclusion policy."""

    normalized = relative.replace("\\", "/")
    findings: list[dict[str, str]] = []
    if PROHIBITED_ARTIFACT_NAME.search(Path(normalized).name):
        findings.append(
            {
                "path": normalized,
                "category": "PROHIBITED_MANUSCRIPT_ARTIFACT",
                "reason": "Main/Supplement/Cover Letter document is not public release material",
            }
        )
    if PROHIBITED_REVISION_RENDER.search(normalized):
        findings.append(
            {
                "path": normalized,
                "category": "PROHIBITED_MANUSCRIPT_REVISION_RENDER",
                "reason": "rendered figure-revision artifact belongs to submission production, not the public scientific tree",
            }
        )

    suffix = Path(normalized).suffix.lower()
    if suffix not in CODE_SUFFIXES:
        if suffix in {".json", ".md", ".txt", ".yml", ".yaml"}:
            text = raw.decode("utf-8", errors="replace")
            if MANUSCRIPT_MEDIA_SOURCE.search(text):
                findings.append(
                    {
                        "path": normalized,
                        "category": "PROHIBITED_MANUSCRIPT_SOURCE_METADATA",
                        "reason": "metadata points to a private word-processing media source",
                    }
                )
        return findings

    text = raw.decode("utf-8", errors="replace")
    if PROHIBITED_SCRIPT_NAME.search(normalized):
        findings.append(
            {
                "path": normalized,
                "category": "PROHIBITED_MANUSCRIPT_AUTHORING_SCRIPT",
                "reason": "script filename identifies manuscript writing, revision, layout or checklist tooling",
            }
        )
    if DOCX_IMPORT.search(text) and DOCUMENT_SAVE.search(text) and MAIN_OR_SUPPLEMENT_TARGET.search(text):
        findings.append(
            {
                "path": normalized,
                "category": "PROHIBITED_MAIN_SUPPLEMENT_DOCX_WRITER",
                "reason": "script imports python-docx, saves a document and targets Main/Supplement content",
            }
        )
    if TEXT_OUTPUT.search(text) and any(literal.lower() in text.lower() for literal in AUTHORING_AID_LITERALS):
        findings.append(
            {
                "path": normalized,
                "category": "PROHIBITED_MANUSCRIPT_WORDING_GENERATOR",
                "reason": "script writes manuscript-remediation or suggested-manuscript wording",
            }
        )
    if IMAGE_OUTPUT.search(text) and MANUSCRIPT_MEDIA_SOURCE.search(text):
        findings.append(
            {
                "path": normalized,
                "category": "PROHIBITED_MANUSCRIPT_FIGURE_LAYOUT_SCRIPT",
                "reason": "script composes a rendered figure from private word-processing media",
            }
        )

    unique: dict[tuple[str, str], dict[str, str]] = {}
    for finding in findings:
        unique[(finding["path"], finding["category"])] = finding
    return list(unique.values())


def audit(root: Path = ROOT, files: Iterable[str] | None = None) -> dict[str, object]:
    selected = list(files) if files is not None else public_files(root)
    findings: list[dict[str, str]] = []
    for relative in selected:
        path = root / relative
        if not path.is_file():
            continue
        findings.extend(classify(relative, path.read_bytes()))
    return {
        "schema": "braintrace.public_release_content_audit.v1",
        "status": "PASS" if not findings else "FAIL",
        "files_scanned": len(selected),
        "prohibited_count": len(findings),
        "findings": findings,
        "policy": (
            "Main/Supplement/Cover Letter artifacts and scripts that write, revise, format or plan "
            "their wording/layout are excluded from every public release surface."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit(args.root.resolve())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    print(json.dumps(payload, ensure_ascii=False))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
