#!/usr/bin/env python
"""Manifest-driven, one-shot BrainTrace release metadata finalization.

The script is intentionally unable to finalize an engineering pre-release.
``--finalize`` requires a separately recorded all-green release gate, then
updates every version-bearing release surface from one manifest in one run.
"""

from __future__ import annotations

import argparse
import json
import re
import tomllib
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "release" / "v0.1.17" / "release_manifest.json"
SOFTWARE_DOI = "10.5281/zenodo.22006038"
FULL_REPRO_DOI = "10.5281/zenodo.22005947"


def load_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_manifest(manifest: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    version = str(manifest.get("version", ""))
    if not re.fullmatch(r"0\.1\.\d+", version):
        errors.append("version must use the 0.1.x form")
    if manifest.get("tag") != f"v{version}":
        errors.append("tag must equal v<version>")
    if manifest.get("software_version_doi") != SOFTWARE_DOI:
        errors.append("software DOI is not the reserved v0.1.17 draft DOI")
    if manifest.get("full_repro_version_doi") != FULL_REPRO_DOI:
        errors.append("full-repro DOI is not the reserved v0.1.17 draft DOI")
    if manifest.get("scientific_frozen") is not True:
        errors.append("scientific_frozen must be true")
    if not isinstance(manifest.get("creators"), list) or not manifest["creators"]:
        errors.append("creators must be a non-empty authoritative list")
    if not str(manifest.get("title", "")).strip() or not str(manifest.get("description", "")).strip():
        errors.append("title and description are required")
    return errors


def _replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"could not uniquely update {label}")
    return updated


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _creators(manifest: dict[str, Any]) -> list[dict[str, str]]:
    return [{"name": str(record["name"])} for record in manifest["creators"]]


def _zenodo_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "title": manifest["title"],
        "description": manifest["description"],
        "creators": _creators(manifest),
        "keywords": list(manifest.get("keywords", [])),
        "license": manifest["license"],
        "version": manifest["version"],
        "related_identifiers": [
            {
                "identifier": manifest["repository"],
                "relation": "isSupplementTo",
                "scheme": "url",
            },
            {
                "identifier": f"https://doi.org/{manifest['software_version_doi']}",
                "relation": "isIdenticalTo",
                "scheme": "url",
            },
            {
                "identifier": f"https://doi.org/{manifest['full_repro_version_doi']}",
                "relation": "isSupplementTo",
                "scheme": "url",
            },
        ],
    }


def _release_block(manifest: dict[str, Any]) -> str:
    return (
        "<!-- BRAINTRACE:CURRENT_RELEASE:START -->\n"
        "## Current release\n\n"
        f"BrainTrace v{manifest['version']} is the current software release. "
        f"Software DOI: [https://doi.org/{manifest['software_version_doi']}](https://doi.org/{manifest['software_version_doi']}). "
        f"Full reproducibility DOI: [https://doi.org/{manifest['full_repro_version_doi']}](https://doi.org/{manifest['full_repro_version_doi']}).\n"
        "<!-- BRAINTRACE:CURRENT_RELEASE:END -->\n"
    )


def _data_provenance_header(manifest: dict[str, Any]) -> str:
    return (
        f"# Data Provenance - BrainTrace v{manifest['version']}\n\n"
        "This manifest documents every external dataset used to build, validate, and benchmark the "
        "BrainTrace hierarchical brain-origin candidate ranking tool. The current executable metadata and immutable "
        f"software release are v{manifest['version']}, archived at version DOI "
        f"[{manifest['software_version_doi']}](https://doi.org/{manifest['software_version_doi']}). "
        "Its full reproducibility archive, including the materialized required Git LFS payload, is archived separately at "
        f"[{manifest['full_repro_version_doi']}](https://doi.org/{manifest['full_repro_version_doi']}); "
        f"the persistent software concept DOI is [{manifest['software_concept_doi']}](https://doi.org/{manifest['software_concept_doi']}). "
        "Public source data are available under the cited accessions except where explicitly noted below. "
    )


def _release_notes(manifest: dict[str, Any]) -> str:
    return (
        f"# BrainTrace v{manifest['version']}\n\n"
        f"- Software DOI: [https://doi.org/{manifest['software_version_doi']}](https://doi.org/{manifest['software_version_doi']})\n"
        f"- Full reproducibility DOI: [https://doi.org/{manifest['full_repro_version_doi']}](https://doi.org/{manifest['full_repro_version_doi']})\n"
        f"- Model lock: `{manifest['model_lock_id']}`\n"
        "- Frozen scientific artifacts were compared against the accepted scientific baseline before release finalization.\n"
        "- External input verification is fail-closed and uses the canonical `external_data/` layout.\n"
    )


def _check_release_gate(gate_path: Path) -> list[str]:
    """Validate facts that exist before the one-shot metadata finalization.

    The exact release-archive checksum is intentionally excluded here: the
    archive can only be built from a clean tree *after* finalization.  It is
    verified as a separate, mandatory pre-merge release-artifact gate.
    """
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    expected = {
        "status": "PASS",
        "scientific_drift": 0,
        "app_docker": "PASS",
        "repro_docker": "PASS",
        "github_actions": "GREEN",
    }
    return [f"release gate {key!r} is not {value!r}" for key, value in expected.items() if gate.get(key) != value]


def finalize(root: Path, manifest_path: Path, gate_path: Path) -> list[Path]:
    manifest = load_manifest(manifest_path)
    errors = [*validate_manifest(manifest), *_check_release_gate(gate_path)]
    if errors:
        raise ValueError("; ".join(errors))
    if manifest.get("release_state") != "engineering_pre_finalization":
        raise ValueError("release manifest is not in engineering_pre_finalization state")

    version = str(manifest["version"])
    pyproject = root / "pyproject.toml"
    _write(pyproject, _replace_once(pyproject.read_text(encoding="utf-8"), r'^version = "[^"]+"$', f'version = "{version}"', "pyproject version"))
    cli = root / "core" / "cli.py"
    _write(cli, _replace_once(cli.read_text(encoding="utf-8"), r'^CURRENT_SOFTWARE_VERSION = "[^"]+"$', f'CURRENT_SOFTWARE_VERSION = "{version}"', "CLI version"))

    docker = root / "Dockerfile"
    docker_text = docker.read_text(encoding="utf-8")
    docker_text = re.sub(r"braintrace:v0\.1\.\d+", f"braintrace:v{version}", docker_text)
    docker_text = _replace_once(docker_text, r'LABEL org\.opencontainers\.image\.version="v[^"]+"', f'LABEL org.opencontainers.image.version="v{version}"', "application Docker version")
    _write(docker, docker_text)
    repro = root / "Dockerfile.repro"
    repro_text = repro.read_text(encoding="utf-8")
    repro_text = re.sub(r"braintrace-repro:v0\.1\.\d+", f"braintrace-repro:v{version}", repro_text)
    repro_text = _replace_once(repro_text, r"ARG BRAINTRACE_VERSION=[^\n]+", f"ARG BRAINTRACE_VERSION={version}", "repro Docker version")
    _write(repro, repro_text)

    readme = root / "README.md"
    readme_text = readme.read_text(encoding="utf-8")
    block = _release_block(manifest)
    if "<!-- BRAINTRACE:CURRENT_RELEASE:START -->" in readme_text:
        readme_text = re.sub(r"<!-- BRAINTRACE:CURRENT_RELEASE:START -->.*?<!-- BRAINTRACE:CURRENT_RELEASE:END -->\n", block, readme_text, flags=re.DOTALL)
    else:
        readme_text = _replace_once(readme_text, r"(?=## Interfaces)", block + "\n", "README current-release block")
    readme_text = re.sub(r"braintrace:v0\.1\.\d+", f"braintrace:v{version}", readme_text)
    _write(readme, readme_text)

    for relative in ("DATA_PROVENANCE.md", "reproducibility/DATA_PROVENANCE.md"):
        path = root / relative
        text = path.read_text(encoding="utf-8")
        marker = "## Source Datasets"
        index = text.find(marker)
        if index < 0:
            raise ValueError(f"could not find source-dataset section in {relative}")
        _write(path, _data_provenance_header(manifest) + "\n" + text[index:])

    _write(root / ".zenodo.json", json.dumps(_zenodo_payload(manifest), ensure_ascii=False, indent=2) + "\n")
    release_dir = root / "release" / f"v{version}"
    _write(release_dir / "RELEASE_NOTES.md", _release_notes(manifest))
    _write(release_dir / "RELEASE_CHECKLIST.md", "# Release checklist\n\n- [x] Version metadata synchronized from `release_manifest.json`.\n- [ ] Upload exact checked archives to the two reserved Zenodo drafts.\n- [ ] Verify remote file checksums and metadata before publication.\n")
    manifest["release_state"] = "final"
    _write(manifest_path, json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    return [pyproject, cli, docker, repro, readme, root / "DATA_PROVENANCE.md", root / "reproducibility" / "DATA_PROVENANCE.md", root / ".zenodo.json", release_dir / "RELEASE_NOTES.md", release_dir / "RELEASE_CHECKLIST.md", manifest_path]


def check(root: Path, manifest_path: Path) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    errors = validate_manifest(manifest)
    pyproject_authors = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8")).get("project", {}).get("authors", [])
    author_names = [str(author.get("name", "")) for author in pyproject_authors if isinstance(author, dict)]
    if author_names != [record["name"] for record in manifest.get("creators", []) if isinstance(record, dict)]:
        errors.append("release-manifest creator order does not match pyproject authors")
    return {
        "schema": "braintrace.release_metadata_sync_audit.v1",
        "version": manifest.get("version"),
        "release_state": manifest.get("release_state"),
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "finalization_ready": manifest.get("release_state") == "engineering_pre_finalization" and not errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--finalize", action="store_true")
    parser.add_argument("--release-gate-record", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check and args.finalize:
        parser.error("--check and --finalize are mutually exclusive")
    if args.finalize:
        if args.release_gate_record is None:
            parser.error("--finalize requires --release-gate-record")
        changed = finalize(ROOT, args.manifest.resolve(), args.release_gate_record.resolve())
        payload: dict[str, Any] = {"status": "PASS", "changed_files": [path.relative_to(ROOT).as_posix() for path in changed]}
    else:
        payload = check(ROOT, args.manifest.resolve())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload, ensure_ascii=False))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
