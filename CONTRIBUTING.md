# Contributing to BrainTrace

Thank you for your interest in contributing.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate  # or .venv\Scripts\Activate.ps1 on Windows
pip install -r requirements_reproducible.txt
pip install -e .
```

## Running tests

```bash
pip install -r requirements-repro.txt
python -m pytest tests/ -v
```

## Code style

- Python 3.11+; type hints are used throughout (`from __future__ import annotations`)
- Four-space indentation; 100-character line limit preferred
- Functions and classes carry NumPy-style docstrings
- Public API changes should update the CLI `--help` output and `README.md`

## Pull request checklist

- [ ] New or modified code is covered by tests in `tests/`
- [ ] `python -m pytest tests/` passes
- [ ] `python scripts/audit_public_release_content.py` passes
- [ ] `python scripts/generate_package_integrity.py --check` passes
- [ ] `README.md` is updated when the external interface changes
- [ ] Reproducibility scripts in `reproducibility/` still produce identical outputs (check `SHA256SUMS.txt`)
- [ ] No submission-document authoring, wording, layout or office-review artifact is included; see `PUBLIC_RELEASE_CONTENT_POLICY.md`

## Repository layout

| Path | Purpose |
|------|---------|
| `app/` | Streamlit web application |
| `core/` | Inference engine, model I/O, metadata |
| `data/models/` | Lightweight frozen model artifacts (`*.npz`, `*.csv`, `*.json`) |
| `scripts/` | Scientific build, validation, provenance and external-analysis scripts |
| `reproducibility/` | Manuscript reproducibility pipelines |
| `tests/` | Unit and integration tests |
| `docs/` | Supplementary documentation |

## Data policy

This repository includes lightweight model artifacts. External source datasets
(Bo2023, AHBA, TCGA, etc.) are **not** included; see `DATA_PROVENANCE.md` and
`SHA256SUMS.txt` for download instructions and integrity checksums.
